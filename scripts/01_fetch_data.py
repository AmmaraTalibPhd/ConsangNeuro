"""
fetch_data.py — Step 1 of the consanguinity / neurological disease pipeline.  (v5)

Builds the input table for the model:
  gene | disorder_group | ancestry | q (cumulative P/LP allele frequency) | ...

Sources (all public, no login needed):
  1. Genomics England PanelApp  -> green (diagnostic-grade) BIALLELIC genes per disorder group
  2. ClinVar variant_summary     -> Pathogenic / Likely pathogenic variants (GRCh38, >=1 review star)
  3. gnomAD v4 GraphQL API       -> ancestry-specific allele counts for those variants

v5 changes
  - Each gene's chromosome is recorded (column `chrom`). X-linked genes are kept in the data
    files but excluded from the recessive model downstream (PanelApp labels X-linked genes
    "biallelic" for affected females); they are used only for the X-linked burden estimate.

v3 changes
  - Super panels (e.g. 'Leukodystrophy, childhood onset', which bundles ~1200 genes) are skipped.
  - P/LP variants with allele frequency > 1% in any ancestry are excluded and listed in
    data/excluded_common_variants.csv (common variants mislabelled pathogenic inflated risk).

v2 changes
  - Genes are matched by stable IDs (HGNC ID for ClinVar, Ensembl ID for gnomAD), not by
    symbol, so renamed genes (AARS->AARS1, GBA->GBA1, MUT->MMUT ...) are no longer lost.
  - Epilepsy panels are matched with a broader keyword (PanelApp renamed the panel).
  - Every web request waits and retries if the internet connection drops, instead of crashing.
  - Cached empty "Gene not found" results from v1 are ignored and fetched again.

Run from the repository root (needs internet):
    python scripts/01_fetch_data.py
Re-run the same command at any time; finished genes are cached and skipped.
"""

import json
import os
import time
import urllib.request

import pandas as pd
import requests

import pathlib
import sys
DATA = pathlib.Path(__file__).resolve().parent.parent / "data"
# Safety: the files in data/ are the exact inputs used in the paper. A fresh download goes to
# data/fresh/ so the deposited inputs are never overwritten by accident. To replace them on purpose:
#     python scripts/01_fetch_data.py --replace-deposited
OUT = str(DATA if ("--replace-deposited" in sys.argv or "--repair-deposited" in sys.argv) else DATA / "fresh")
os.makedirs(OUT, exist_ok=True)
CACHE = os.path.join(str(DATA), "gnomad_cache")      # shared cache, safe to reuse
os.makedirs(CACHE, exist_ok=True)

# Disorder groups -> PanelApp panel name keywords (matched case-insensitively)
PANEL_GROUPS = {
    "ID_NDD": ["intellectual disability"],
    "Epilepsy": ["epilep"],
    "HSP": ["hereditary spastic paraplegia"],
    "Leukodystrophy": ["leukodystrophy", "white matter"],
}

ANCESTRIES = ["sas", "mid", "nfe", "afr", "eas", "amr"]  # gnomAD v4 genetic ancestry groups
PLP = {"Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"}
MAX_VARIANT_AF = 0.01   # drop P/LP variants with AF > 1% in any ancestry (too common for a
                        # severe recessive disease: likely misclassified or low-penetrance)
MAX_GENES_PER_PANEL = 500  # a disease-specific panel this large is almost certainly a super panel
MIN_STARS_STATUS = {
    "criteria provided, single submitter",
    "criteria provided, multiple submitters, no conflicts",
    "reviewed by expert panel",
    "practice guideline",
}


# ---------------------------------------------------------------- robust web requests
def web(method, url, max_wait_min=60, **kw):
    """requests.get/post that survives dropped connections and rate limits."""
    waited, delay = 0, 15
    while True:
        try:
            r = requests.request(method, url, timeout=120, **kw)
            if r.status_code == 429 or r.status_code >= 500:
                raise requests.exceptions.ConnectionError(f"HTTP {r.status_code}")
            return r
        except (requests.exceptions.ConnectionError, requests.exceptions.Timeout) as e:
            if waited >= max_wait_min * 60:
                raise
            print(f"\n  Connection problem ({str(e)[:80]}). Retrying in {delay}s...")
            time.sleep(delay)
            waited += delay
            delay = min(delay * 2, 300)


# ---------------------------------------------------------------- 1. PanelApp
def is_super_panel(p):
    # Only a true super panel. NOT "Component Of Super Panel", which is an ordinary panel
    # (e.g. Intellectual disability, Early onset or syndromic epilepsy).
    for t in (p.get("types") or []):
        name = (t.get("name") or "").strip().lower()
        slug = (t.get("slug") or "").strip().lower()
        if name == "super panel" or slug in ("superpanel", "super-panel"):
            return True
    return False


def chromosome(gene_data):
    try:
        for ver in gene_data["ensembl_genes"]["GRch38"].values():
            return str(ver["location"]).split(":")[0]
    except (KeyError, TypeError, AttributeError):
        return None


def ensembl_id(gene_data):
    try:
        for ver in gene_data["ensembl_genes"]["GRch38"].values():
            return ver["ensembl_id"]
    except (KeyError, TypeError, AttributeError):
        return None


def fetch_panelapp():
    base = "https://panelapp.genomicsengland.co.uk/api/v1/panels/"
    panels, url = [], base
    while url:
        r = web("GET", url).json()
        panels += r["results"]
        url = r.get("next")
    rows = []
    for group, keys in PANEL_GROUPS.items():
        hits = [p for p in panels if any(k in p["name"].lower() for k in keys)
                and not is_super_panel(p)]
        if not hits:
            print(f"  WARNING: no PanelApp panel matched for {group}")
        for p in hits:
            detail = web("GET", f"{base}{p['id']}/").json()
            n_green_bi = sum(1 for g in detail["genes"] if g.get("confidence_level") == "3"
                             and "BIALLELIC" in (g.get("mode_of_inheritance") or "").upper())
            if group != "ID_NDD" and n_green_bi > MAX_GENES_PER_PANEL:
                print(f"  SKIPPED (looks like a super panel, {n_green_bi} genes): {p['name']}")
                continue
            for g in detail["genes"]:
                moi = (g.get("mode_of_inheritance") or "").upper()
                if g.get("confidence_level") == "3" and "BIALLELIC" in moi and "MONOALLELIC" not in moi:
                    gd = g["gene_data"]
                    rows.append({"gene": gd["gene_symbol"], "hgnc_id": gd.get("hgnc_id"),
                                 "chrom": chromosome(gd),
                                 "ensembl_id": ensembl_id(gd), "disorder_group": group,
                                 "panel": p["name"], "panel_version": detail.get("version")})
            print(f"  PanelApp: {group:15s} <- {p['name']} (v{detail.get('version')})")
    df = pd.DataFrame(rows).drop_duplicates(["gene", "disorder_group"])
    missing = [g for g in PANEL_GROUPS if g not in set(df.disorder_group)]
    if missing:
        raise SystemExit(f"STOP: no genes found for {missing}. Send Claude the console output above.")
    df.to_csv(os.path.join(OUT, "genes_panelapp.csv"), index=False)
    print(f"PanelApp: {df.gene.nunique()} unique biallelic green genes "
          f"({df.ensembl_id.isna().sum()} rows without Ensembl ID)")
    print(df.groupby("disorder_group").gene.nunique().to_string())
    return df


# ---------------------------------------------------------------- 2. ClinVar
def _gzip_ok(path):
    import gzip
    try:
        with gzip.open(path, "rb") as f:
            while f.read(1 << 24):
                pass
        return True
    except (OSError, EOFError):
        return False


def download_resumable(url, path, tries=50):
    """Download a large file, resuming after dropped connections; verify the gzip at the end."""
    if os.path.exists(path) and _gzip_ok(path):
        return
    part = path + ".part"
    if os.path.exists(path):                       # truncated file from an earlier failed run
        os.replace(path, part)
    total = None
    for k in range(tries):
        done = os.path.getsize(part) if os.path.exists(part) else 0
        try:
            hdr = {"Range": f"bytes={done}-"} if done else {}
            with requests.get(url, headers=hdr, stream=True, timeout=120) as r:
                if r.status_code == 416:               # already complete
                    break
                r.raise_for_status()
                if done and r.status_code != 206:      # server ignored the range: start again
                    done = 0
                if r.status_code == 206:
                    total = int(r.headers["Content-Range"].split("/")[-1])
                else:
                    total = int(r.headers.get("Content-Length", 0)) or None
                with open(part, "ab" if done else "wb") as f:
                    for chunk in r.iter_content(1 << 20):
                        f.write(chunk); done += len(chunk)
                        if total:
                            print(f"  ClinVar download: {done / 1e6:,.0f} / {total / 1e6:,.0f} MB", end="\r", flush=True)
            if total is None or done >= total:
                break
        except requests.RequestException as e:
            print(f"\n  connection dropped ({e}); resuming in 30 s", flush=True)
            time.sleep(30)
    print()
    if not _gzip_ok(part):
        os.remove(part)
        raise SystemExit("ClinVar download is corrupt; please run the script again.")
    os.replace(part, path)

def fetch_clinvar(genes_df):
    path = os.path.join(str(DATA), "variant_summary.txt.gz")
    download_resumable("https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz", path)
    cols = ["GeneSymbol", "HGNC_ID", "ClinicalSignificance", "ReviewStatus", "Assembly", "Chromosome",
            "PositionVCF", "ReferenceAlleleVCF", "AlternateAlleleVCF", "Type", "VariationID"]
    cv = pd.read_csv(path, sep="\t", usecols=cols, low_memory=False, dtype=str)
    hgnc = set(genes_df.hgnc_id.dropna())
    symbols = set(genes_df.gene)
    in_set = cv.HGNC_ID.isin(hgnc) | cv.GeneSymbol.isin(symbols)
    cv = cv[(cv.Assembly == "GRCh38") & in_set
            & cv.ClinicalSignificance.isin(PLP) & cv.ReviewStatus.isin(MIN_STARS_STATUS)]
    cv = cv[cv.ReferenceAlleleVCF.ne("na") & cv.AlternateAlleleVCF.ne("na")]
    cv["variant_id"] = cv.Chromosome + "-" + cv.PositionVCF + "-" + cv.ReferenceAlleleVCF + "-" + cv.AlternateAlleleVCF
    cv = cv.drop_duplicates("variant_id")
    cv.to_csv(os.path.join(OUT, "clinvar_plp.csv"), index=False)
    print(f"ClinVar: {len(cv)} P/LP variants in {cv.HGNC_ID.nunique()} genes")
    return cv


# ---------------------------------------------------------------- 3. gnomAD
GQL = """
query($id: String!) {
  gene(gene_id: $id, reference_genome: GRCh38) {
    symbol
    variants(dataset: gnomad_r4) {
      variant_id
      exome  { ac an populations { id ac an } }
      genome { ac an populations { id ac an } }
    }
  }
}"""
GQL_SYMBOL = GQL.replace("gene_id: $id", "gene_symbol: $id")


def _has_variants(data):
    return bool(((data or {}).get("data") or {}).get("gene"))


FAILED = []   # genes gnomAD could not serve this run (not cached; fetched again next run)


def _permanent(msg):
    """Only 'gene not found'-type answers are real; overload / unknown errors are temporary."""
    return "not found" in (msg or "").lower()


def _usable_cache(path):
    if not os.path.exists(path):
        return None
    data = json.load(open(path))
    if _has_variants(data) and data["data"]["gene"].get("variants") is not None and \
            ("note" not in data or _permanent(data["note"])):
        return data
    return None                     # empty placeholder from a temporary error: fetch again


def _fetch(key, query, var, label):
    """One gnomAD gene query with cache; temporary errors are retried, never cached."""
    cache = os.path.join(CACHE, f"{key}.json")
    data = _usable_cache(cache)
    if data is not None:
        return data
    for attempt in range(8):
        data = web("POST", "https://gnomad.broadinstitute.org/api",
                   json={"query": query, "variables": {"id": var}}).json()
        if _has_variants(data):
            json.dump(data, open(cache, "w"))
            time.sleep(6)                                   # stay under the public API rate limit
            return data
        msg = (data.get("errors") or [{"message": "no data"}])[0]["message"]
        if _permanent(msg):
            data = {"data": {"gene": {"variants": []}}, "note": msg}
            json.dump(data, open(cache, "w"))
            return data
        wait = min(60 * 2 ** attempt, 900)
        print(f"\n  gnomAD busy for {label} ({msg}); waiting {wait // 60} min and trying again", flush=True)
        time.sleep(wait)
    return None                                             # still failing: not cached


def gnomad_gene(symbol, ens):
    old = _usable_cache(os.path.join(CACHE, f"{symbol}.json"))      # v1 cache keyed by symbol
    if old is not None and old["data"]["gene"].get("variants"):
        return old
    data = _fetch(ens, GQL, ens, symbol) if ens else None
    if ens and data is None:
        FAILED.append(symbol)
        return {"data": {"gene": {"variants": []}}}
    if data is None or "note" in data:
        # gnomAD (GENCODE v39) sometimes uses a different Ensembl ID than PanelApp, e.g. TBCE, SCO2:
        # fall back to the gene symbol before concluding that the gene is absent.
        data = _fetch(f"symbol_{symbol}", GQL_SYMBOL, symbol, symbol)
        if data is None:
            FAILED.append(symbol)
            return {"data": {"gene": {"variants": []}}}
        if "note" in data:
            print(f"\n  gnomAD: {symbol} not found by Ensembl ID or symbol: {data['note']}")
    return data


def ancestry_counts(v):
    """Combine exome + genome allele counts per ancestry group."""
    out = {}
    for src in ("exome", "genome"):
        d = v.get(src)
        if not d:
            continue
        for p in d["populations"]:
            if p["id"] in ANCESTRIES:  # skips sex-split ids like 'sas_XX'
                a = out.setdefault(p["id"], [0, 0])
                a[0] += p["ac"]; a[1] += p["an"]
    return out


def build_gene_freqs(genes_df, cv):
    rows, var_rows, excluded = [], [], []
    plp_all = set(cv.variant_id)
    genes = genes_df.drop_duplicates("gene")[["gene", "hgnc_id", "ensembl_id"]].sort_values("gene")
    n_plp = cv.groupby("HGNC_ID").size().to_dict()
    for i, g in enumerate(genes.itertuples(), 1):
        print(f"  gnomAD {i}/{len(genes)}: {g.gene:12s}", end="\r")
        ens = g.ensembl_id if isinstance(g.ensembl_id, str) else None
        variants = (gnomad_gene(g.gene, ens).get("data", {}).get("gene") or {}).get("variants") or []
        acc = {a: {"ac_sum": 0, "q": 0.0, "ans": [], "n": 0} for a in ANCESTRIES}
        for v in variants:
            if v["variant_id"] not in plp_all:
                continue
            counts = ancestry_counts(v)
            max_af = max((ac / an for ac, an in counts.values() if an > 0), default=0)
            if max_af > MAX_VARIANT_AF:
                excluded.append({"gene": g.gene, "variant_id": v["variant_id"], "max_ancestry_af": max_af})
                continue
            for a, (ac, an) in counts.items():
                var_rows.append({"gene": g.gene, "variant_id": v["variant_id"], "ancestry": a, "ac": ac, "an": an})
                if an > 0:
                    acc[a]["q"] += ac / an
                    acc[a]["ac_sum"] += ac
                    acc[a]["ans"].append(an)
                    acc[a]["n"] += ac > 0
        for a, d in acc.items():
            rows.append({"gene": g.gene, "ancestry": a, "q": min(d["q"], 0.5),
                         "ac_total": d["ac_sum"], "n_variants_observed": d["n"],
                         "an_median": int(pd.Series(d["ans"]).median()) if d["ans"] else 0,
                         "n_plp_clinvar": n_plp.get(g.hgnc_id, 0)})
    freqs = pd.DataFrame(rows).merge(genes_df[["gene", "disorder_group"]], on="gene")
    freqs.to_csv(os.path.join(OUT, "gene_freqs.csv"), index=False)
    pd.DataFrame(var_rows).to_csv(os.path.join(OUT, "variant_counts.csv"), index=False)
    pd.DataFrame(excluded).to_csv(os.path.join(OUT, "excluded_common_variants.csv"), index=False)
    print(f"\nExcluded {len(excluded)} P/LP variants with AF > {MAX_VARIANT_AF:.0%} "
          f"(see {OUT}/excluded_common_variants.csv)")
    empty = genes.gene[~genes.gene.isin(pd.DataFrame(var_rows).get("gene", pd.Series(dtype=str)))]
    print(f"\nWrote {OUT}/gene_freqs.csv  ({freqs.gene.nunique()} genes x {len(ANCESTRIES)} ancestries)")
    print(f"{len(empty)} genes had no ClinVar P/LP variant present in gnomAD (kept with q = 0).")


if __name__ == "__main__":
    if "--repair-deposited" in sys.argv:
        # Rebuild gnomAD frequencies for the paper's own gene list and ClinVar variants (no new PanelApp
        # or ClinVar download). Cached genes are reused; only genes missing from the cache are fetched.
        genes_df = pd.read_csv(os.path.join(OUT, "genes_panelapp.csv"))
        cv = pd.read_csv(os.path.join(OUT, "clinvar_plp.csv"), dtype=str)
        print(f"Repairing deposited data: {genes_df.gene.nunique()} genes, {len(cv)} ClinVar variants")
    else:
        genes_df = fetch_panelapp()
        cv = fetch_clinvar(genes_df)
    build_gene_freqs(genes_df, cv)
    if FAILED:
        raise SystemExit(f"\n{len(FAILED)} genes could not be fetched because gnomAD was overloaded: "
                         f"{', '.join(FAILED)}.\nRun the same command again later; finished genes are cached "
                         "and only these will be fetched. Do not use the output files until this message is gone.")
    print("Done. Next: python scripts/02_model.py")
