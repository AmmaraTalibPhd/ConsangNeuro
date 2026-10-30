"""
05_literature_variants.py — how many causal variants from published Pakistani families
would the model count?  (Table S10)

Input : data/external/riazuddin2017_variants.tsv (Riazuddin et al. 2017, Tables 1-3; GRCh37)
Needs : pip install pyliftover ; internet once, to download the GRCh37->GRCh38 chain file
Run   : python scripts/05_literature_variants.py
"""
import os
import re
import sys
import urllib.request

import pandas as pd
from pyliftover import LiftOver

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402

CHAIN_URL = ("https://raw.githubusercontent.com/broadinstitute/gatk/master/scripts/funcotator/"
             "data_sources/gnomAD/b37ToHg38.over.chain")
CHAIN = C.DATA / "b37ToHg38.over.chain"
if not CHAIN.exists():
    urllib.request.urlretrieve(CHAIN_URL, CHAIN)
lo = LiftOver(str(CHAIN))
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}

lit = pd.read_csv(C.EXT / "riazuddin2017_variants.tsv", sep="\t")
lit = lit[lit.chrom.astype(str) != "X"].copy()
lit["pos_b38"] = [(lambda x: x[0][1] + 1 if x else None)(lo.convert_coordinate(str(c), int(p) - 1))
                  for c, p in zip(lit.chrom, lit.pos_b37)]
genes = C.genes_table()
panel = set(genes.loc[~genes.is_x, "gene"])
cv = pd.read_csv(C.DATA / "clinvar_plp.csv", dtype=str)
cv = cv[cv.variant_id.notna()].copy()
cv["pos"] = cv.PositionVCF.astype(int)
vc = pd.read_csv(C.DATA / "variant_counts.csv")


def match(row):
    """Exact match: substitution at same position with consistent alleles (either strand),
    or deletion of the reported bases left-anchored one base upstream.
    Substitution with consistent alleles one base away = 'uncertain'."""
    cand = cv[(cv.Chromosome == str(row.chrom)) & ((cv.pos - row.pos_b38).abs() <= 3)]
    m_sub = re.match(r"c\.[\d+\-*]+([ACGT])>([ACGT])", row.cdna)
    m_del = re.match(r"c\.[\d_+\-]+del([ACGT]+)$", row.cdna)
    for x in cand.itertuples():
        _, p, ref, alt = x.variant_id.split("-")
        d = int(p) - row.pos_b38
        if m_sub and len(ref) == 1 and len(alt) == 1:
            a, b = m_sub.groups()
            if (ref, alt) in [(a, b), (COMP[a], COMP[b])]:
                if d == 0:
                    return "confident", x.variant_id
                if abs(d) == 1:
                    return "uncertain", x.variant_id
        if m_del and len(alt) == 1 and len(ref) == len(m_del.group(1)) + 1 and d == -1:
            deleted = ref[1:]
            rc = "".join(COMP[c] for c in reversed(m_del.group(1)))
            if deleted in (m_del.group(1), rc):
                return "confident", x.variant_id
    return "none", ""


res = lit.apply(lambda r: pd.Series(match(r), index=["clinvar_match", "variant_id_b38"]), axis=1)
lit = pd.concat([lit, res], axis=1)
lit["gene_in_panel"] = lit.gene.isin(panel)
for anc in ("sas", "nfe"):
    ac = vc[vc.ancestry == anc].groupby("variant_id").ac.sum()
    lit[f"{anc}_ac"] = lit.variant_id_b38.map(ac).fillna(0).astype(int)
lit.to_csv(C.SUPP / "TableS10_pakistani_family_variants.csv", index=False)

n = len(lit)
conf = (lit.clinvar_match == "confident").sum(); unc = (lit.clinvar_match == "uncertain").sum()
counted_lo = ((lit.clinvar_match == "confident") & (lit.sas_ac > 0)).sum()
counted_hi = ((lit.clinvar_match != "none") & (lit.sas_ac > 0)).sum()
print(f"{n} autosomal variants in {lit.family.nunique()} families")
print(f"gene in panel: {lit.gene_in_panel.sum()} ({100 * lit.gene_in_panel.mean():.0f}%)")
print(f"ClinVar P/LP: {conf}-{conf + unc} ({100 * conf / n:.0f}-{100 * (conf + unc) / n:.0f}%)")
print(f"counted by model: {counted_lo}-{counted_hi} ({100 * counted_lo / n:.0f}-{100 * counted_hi / n:.0f}%)")
step = pd.Series("counted", index=lit.index)
step[(lit.clinvar_match != "none") & (lit.sas_ac == 0)] = "classified but not seen in gnomAD South Asians"
step[lit.clinvar_match == "none"] = "on gene list but not classified in ClinVar"
step[~lit.gene_in_panel] = "gene not on diagnostic-grade list"
lit["failure_step"] = step
lit.to_csv(C.SUPP / "TableS10_pakistani_family_variants.csv", index=False)
print(step.value_counts().to_string())
