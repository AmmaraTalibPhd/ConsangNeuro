"""
06_plof_analysis.py — classification-independent test of the ancestry ranking (Table S14).

Question: is European predicted risk higher than South Asian risk because European pathogenic alleles
are genuinely more frequent, or because ClinVar has classified more of them?

Approach: repeat the first-cousin risk calculation using gnomAD high-confidence predicted loss-of-function
(pLoF) variants, chosen by LOFTEE annotation alone and NOT by ClinVar classification, and compare the
European/South Asian ratio with the ratio from ClinVar P/LP variants in the same genes.

Two steps
  1. FETCH (needs internet, ~1.5-2 h, resumable): queries gnomAD v4 per gene and writes
     data/plof_variant_counts.csv (variant x ancestry allele counts with LOFTEE annotation).
     Run:  python scripts/06_plof_analysis.py --fetch
  2. ANALYSE (seconds; uses the deposited CSV):
     Run:  python scripts/06_plof_analysis.py

pLoF definition: LOFTEE high-confidence ("HC"), no LOFTEE flags, consequence stop_gained /
frameshift_variant / splice_donor_variant / splice_acceptor_variant, PASS in at least one of exome or
genome, allele frequency <= 1% in every ancestry (same ceiling as the ClinVar analysis).

Gene sets
  Primary  : "LoF-mechanism genes" — genes in which at least one HC pLoF variant observed in gnomAD is
             classified P/LP in ClinVar (ClinVar is used only to choose GENES whose disease mechanism
             includes loss of function, never to choose VARIANTS).
  Secondary: all 948 autosomal genes.
"""
import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd
import requests

sys.path.insert(0, os.path.dirname(__file__))
import common as C  # noqa: E402

OUT = C.DATA / "plof_variant_counts.csv"
CACHE = C.DATA / "gnomad_lof_cache"
LOF_CSQ = {"stop_gained", "frameshift_variant", "splice_donor_variant", "splice_acceptor_variant"}

GQL = """
query($id: String!) {
  gene(gene_id: $id, reference_genome: GRCh38) {
    symbol
    variants(dataset: gnomad_r4) {
      variant_id
      consequence
      lof
      lof_filter
      lof_flags
      exome  { ac an filters populations { id ac an } }
      genome { ac an filters populations { id ac an } }
    }
  }
}"""
GQL_SYMBOL = GQL.replace("gene_id: $id", "gene_symbol: $id")


# ---------------------------------------------------------------- 1. fetch
def post(query, var, tries=30):
    for k in range(tries):
        try:
            r = requests.post("https://gnomad.broadinstitute.org/api",
                              json={"query": query, "variables": {"id": var}}, timeout=120)
            if r.status_code == 429:
                raise requests.RequestException("rate limited")
            r.raise_for_status()
            return r.json()
        except requests.RequestException as e:
            wait = min(60 * (k + 1), 600)
            print(f"\n  connection problem ({e}); retrying in {wait // 60} min", flush=True)
            time.sleep(wait)
    raise SystemExit("gnomAD could not be reached; re-run later (finished genes are cached).")


def fetch():
    CACHE.mkdir(exist_ok=True)
    genes = C.genes_table()
    genes = genes[~genes.is_x].drop_duplicates("gene")[["gene", "ensembl_id"]].sort_values("gene")
    rows, failed = [], []
    for i, g in enumerate(genes.itertuples(), 1):
        print(f"  gnomAD pLoF {i}/{len(genes)}: {g.gene:12s}", end="\r", flush=True)
        ens = g.ensembl_id if isinstance(g.ensembl_id, str) and g.ensembl_id.startswith("ENSG") else None
        f = CACHE / f"{ens or g.gene}.json"
        data = json.load(open(f)) if f.exists() else None
        if data is not None and "note" in data and "not found" not in data["note"].lower():
            data = None                                   # placeholder from a temporary error: fetch again
        if data is None:
            for attempt in range(8):
                data = post(GQL, ens) if ens else post(GQL_SYMBOL, g.gene)
                if (data.get("data") or {}).get("gene"):
                    json.dump(data, open(f, "w")); time.sleep(6)
                    break
                msg = (data.get("errors") or [{"message": "no data"}])[0]["message"]
                if "Cannot query field" in msg:
                    raise SystemExit(f"gnomAD API schema changed: {msg}. Please report this message.")
                if "not found" in msg.lower() and ens:
                    ens = None                       # gnomAD uses another Ensembl ID: retry by symbol
                    f = CACHE / f"symbol_{g.gene}.json"
                    if f.exists():
                        data = json.load(open(f)); break
                    data = None
                    continue
                if "not found" in msg.lower():
                    print(f"\n  gnomAD: {g.gene} not in gnomAD: {msg}")
                    data = {"data": {"gene": {"variants": []}}, "note": msg}
                    json.dump(data, open(f, "w"))
                    break
                wait = min(60 * 2 ** attempt, 900)
                print(f"\n  gnomAD busy for {g.gene} ({msg}); waiting {wait // 60} min and trying again", flush=True)
                time.sleep(wait)
                data = None
            if data is None:
                failed.append(g.gene)
                continue
        for v in ((data.get("data") or {}).get("gene") or {}).get("variants") or []:
            if v.get("consequence") not in LOF_CSQ:
                continue
            passed = any(v.get(src) and not v[src].get("filters") for src in ("exome", "genome"))
            counts = {}
            for src in ("exome", "genome"):
                d = v.get(src)
                if not d:
                    continue
                for p in d["populations"]:
                    if p["id"] in C.ANCESTRIES:
                        a = counts.setdefault(p["id"], [0, 0]); a[0] += p["ac"]; a[1] += p["an"]
            for anc, (ac, an) in counts.items():
                rows.append(dict(gene=g.gene, variant_id=v["variant_id"], consequence=v["consequence"],
                                 lof=v.get("lof"), lof_flags=v.get("lof_flags") or "", lof_filter=v.get("lof_filter") or "",
                                 pass_filter=passed, ancestry=anc, ac=ac, an=an))
    if failed:
        raise SystemExit(f"\n{len(failed)} genes could not be fetched because gnomAD was overloaded: {', '.join(failed)}.\n"
                         "Run the same command again later; finished genes are cached and only these will be fetched.")
    pd.DataFrame(rows).to_csv(OUT, index=False)
    print(f"\nWrote {OUT} ({len(rows)} rows)")


# ---------------------------------------------------------------- 2. analyse
def gene_q(v):
    w = v[v.an > 0]
    return (w.assign(af=w.ac / w.an).groupby(["gene", "ancestry"]).af.sum().clip(upper=0.5)
            .unstack().reindex(columns=C.ANCESTRIES).fillna(0))


def first_cousin(qmat, genes):
    return {a: 1e5 * C.group_risk(qmat.reindex(sorted(genes)).fillna(0)[a].values, 1 / 16) for a in C.ANCESTRIES}


def analyse():
    if not OUT.exists():
        raise SystemExit(f"{OUT} not found. Run first:  python scripts/06_plof_analysis.py --fetch")
    raw = pd.read_csv(OUT, keep_default_na=False, dtype={"lof_flags": str, "lof_filter": str})
    raw["pass_filter"] = raw.pass_filter.astype(str).str.lower() == "true"
    auto = C.genes_table(); auto = auto[~auto.is_x]
    grp = auto.groupby("gene").disorder_group.apply(set)
    v = raw[(raw.lof == "HC") & (raw.lof_flags == "") & raw.pass_filter & raw.consequence.isin(LOF_CSQ)].copy()
    v["af"] = v.ac / v.an.where(v.an > 0)
    v = v[v.variant_id.map(v.groupby("variant_id").af.max()) <= C.MAX_VARIANT_AF]
    plp = set(pd.read_csv(C.DATA / "clinvar_plp.csv", usecols=["variant_id"]).variant_id)
    v["clinvar_plp"] = v.variant_id.isin(plp)

    lof_genes = set(v[v.clinvar_plp & (v.ac > 0)].gene)
    all_genes = set(auto.gene)
    id_genes = {g for g in all_genes if "ID_NDD" in grp.get(g, set())}

    q_plof = gene_q(v)
    q_cv = C.gene_freqs().drop_duplicates(["gene", "ancestry"]).pivot(index="gene", columns="ancestry", values="q")
    q_cv_in_plof = gene_q(v[v.clinvar_plp])          # ClinVar-classified subset of the same pLoF variants

    rows = []
    for setname, genes in [("LoF-mechanism genes (primary)", lof_genes), ("All autosomal genes", all_genes)]:
        for grpname, gs in [("ID_NDD", genes & id_genes), ("All groups", genes)]:
            r_p, r_c, r_k = first_cousin(q_plof, gs), first_cousin(q_cv, gs), first_cousin(q_cv_in_plof, gs)
            for a in C.ANCESTRIES:
                rows.append(dict(gene_set=setname, disorder_group=grpname, n_genes=len(gs), ancestry=a,
                                 risk_hc_plof=r_p[a], risk_clinvar_plp=r_c[a], risk_clinvar_classified_plof=r_k[a]))
    S = pd.DataFrame(rows)

    # share of pLoF allele frequency that ClinVar has classified, per ancestry (primary gene set, ID/NDD)
    w = v[v.gene.isin(lof_genes & id_genes) & (v.an > 0)]
    share = w.groupby("ancestry").apply(lambda d: d.loc[d.clinvar_plp, "af"].sum() / d.af.sum() if d.af.sum() else np.nan)
    nvar = w[w.ac > 0].groupby("ancestry").apply(lambda d: d.clinvar_plp.mean())

    # European downsampling of pLoF counts
    rs = np.random.default_rng(11); ds = {}
    nfe = v[(v.ancestry == "nfe") & (v.an > 0)]
    G = sorted(lof_genes & id_genes)
    for tgt in ("sas", "mid"):
        N = int(v[(v.ancestry == tgt) & (v.an > 0)].groupby("variant_id").an.max().median())
        res = []
        for _ in range(300):
            n = np.minimum(nfe.an.values, N)
            ac = rs.hypergeometric(nfe.ac.values, (nfe.an - nfe.ac).values, n)
            q = pd.Series(ac / n).groupby(nfe.gene.values).sum().clip(upper=0.5).reindex(G).fillna(0).values
            res.append(1e5 * C.group_risk(q, 1 / 16))
        ds[tgt] = dict(target_AN=N, median=float(np.median(res)), lo=float(np.quantile(res, .025)), hi=float(np.quantile(res, .975)))

    # Sensitivity analyses: is the pLoF ratio driven by a few variants or genes?
    G_id = lof_genes & id_genes

    def ratio(vv, genes=G_id):
        r = first_cousin(gene_q(vv), genes)
        return dict(sas=r["sas"], mid=r["mid"], nfe=r["nfe"], afr=r["afr"], nfe_sas=r["nfe"] / r["sas"], nfe_mid=r["nfe"] / r["mid"])

    vmax = v.variant_id.map(v.groupby("variant_id").af.max())
    contrib = {a: pd.Series(C.risk(q_plof.reindex(sorted(G_id)).fillna(0)[a].values, 1 / 16), index=sorted(G_id)).sort_values(ascending=False)
               for a in ("sas", "nfe")}
    sens = [("Primary (HC pLoF, AF <= 1%)", ratio(v)),
            ("AF <= 0.5% in every ancestry", ratio(v[vmax <= 0.005])),
            ("AF <= 0.1% in every ancestry", ratio(v[vmax <= 0.001])),
            ("Stop-gained and frameshift only (no splice)", ratio(v[v.consequence.isin({"stop_gained", "frameshift_variant"})])),
            ("Without ClinVar-classified variants (unclassified pLoF only)", ratio(v[~v.clinvar_plp]))]
    for k in (10, 25):
        for a in ("sas", "nfe"):
            sens.append((f"Without top {k} genes by {C.ANC_LABEL[a]} pLoF risk", ratio(v, G_id - set(contrib[a].index[:k]))))
    SENS = pd.DataFrame([dict(analysis=n, **r) for n, r in sens])
    SENS.round(3).to_csv(C.SUPP / "TableS14c_plof_sensitivity.csv", index=False)
    both = q_plof.reindex(sorted(G_id)).fillna(0)
    both = both[(both.sas > 0) & (both.nfe > 0)]
    gene_cmp = dict(n_both=len(both), share_sas_higher=float((both.sas > both.nfe).mean()),
                    median_nfe_sas_ratio=float((both.nfe / both.sas).median()))

    S.round(2).to_csv(C.SUPP / "TableS14_plof_ancestry_comparison.csv", index=False)
    gl = q_plof.reindex(sorted(lof_genes)).add_prefix("q_plof_").join(q_cv.reindex(sorted(lof_genes)).add_prefix("q_clinvar_"))
    gl.round(6).to_csv(C.SUPP / "TableS14b_plof_gene_frequencies.csv")

    p = S[(S.gene_set.str.startswith("LoF")) & (S.disorder_group == "ID_NDD")].set_index("ancestry")
    K = {"n_hc_plof_variants": int(v.variant_id.nunique()), "n_lof_mechanism_genes": len(lof_genes),
         "n_lof_mechanism_ID_genes": len(lof_genes & id_genes),
         "ID_first_cousin_per100k": p[["risk_hc_plof", "risk_clinvar_plp", "risk_clinvar_classified_plof"]].round(1).to_dict(),
         "nfe_sas_ratio_plof": p.risk_hc_plof["nfe"] / p.risk_hc_plof["sas"],
         "nfe_sas_ratio_clinvar": p.risk_clinvar_plp["nfe"] / p.risk_clinvar_plp["sas"],
         "nfe_mid_ratio_plof": p.risk_hc_plof["nfe"] / p.risk_hc_plof["mid"],
         "nfe_mid_ratio_clinvar": p.risk_clinvar_plp["nfe"] / p.risk_clinvar_plp["mid"],
         "share_plof_af_classified": share.round(3).to_dict(),
         "share_observed_plof_variants_classified": nvar.round(3).to_dict(),
         "nfe_downsampled_plof_ID_per100k": ds,
         "sensitivity": SENS.set_index("analysis")[["nfe_sas", "nfe_mid"]].round(3).to_dict("index"),
         "gene_by_gene_plof": gene_cmp,
         "classified_share_of_plof_risk": {a: p.risk_clinvar_classified_plof[a] / p.risk_hc_plof[a] for a in C.ANCESTRIES}}
    json.dump(K, open(C.RESULTS / "plof_key_numbers.json", "w"), indent=1, default=float)

    # Figure S2
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.6))
    x = np.arange(len(C.FOCUS)); wd = 0.38
    axs[0].bar(x - wd / 2, [p.risk_clinvar_plp[a] for a in C.FOCUS], wd, color="#9ecae1", label="ClinVar P/LP variants")
    axs[0].bar(x + wd / 2, [p.risk_hc_plof[a] for a in C.FOCUS], wd, color="#3182bd", label="High-confidence pLoF (classification-independent)")
    axs[0].set_xticks(x, [C.ANC_LABEL[a] for a in C.FOCUS]); axs[0].set_ylabel("Predicted ID/NDD risk per 100,000\nfirst-cousin births")
    axs[0].legend(frameon=False, fontsize=8, loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=1)
    axs[0].tick_params(axis="x", labelsize=8); axs[1].tick_params(axis="x", labelsize=8)
    axs[1].bar(x, [100 * share.get(a, np.nan) for a in C.FOCUS], color=[C.ANC_COLOR[a] for a in C.FOCUS])
    axs[1].set_xticks(x, [C.ANC_LABEL[a] for a in C.FOCUS]); axs[1].set_ylabel("pLoF allele frequency classified\nP/LP in ClinVar (%)")
    for ax, l in zip(axs, "ab"):
        ax.set_title(l, loc="left", fontweight="bold"); ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(); fig.savefig(C.FIGURES / "FigureS2_plof.png", dpi=300); plt.close(fig)
    print(json.dumps(K, indent=1, default=float))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--fetch", action="store_true")
    if ap.parse_args().fetch:
        fetch()
    analyse()
