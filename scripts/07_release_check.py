"""
07_release_check.py — do the results depend on the ClinVar / PanelApp release? (Table S15)

Reruns scripts 02 and 03 on a newer download in data/fresh/ (made by `python scripts/01_fetch_data.py`)
and compares the headline numbers with those from the paper's inputs in data/.
Outputs go to release_check/ ; the comparison table to supplementary/TableS15_release_comparison.csv.
"""
import json
import os
import subprocess
import sys

import pandas as pd

import common as C

FRESH = C.ROOT / "data" / "fresh"
OUT = C.ROOT / "release_check"
if not (FRESH / "gene_freqs.csv").exists():
    raise SystemExit("data/fresh/gene_freqs.csv not found: run  python scripts/01_fetch_data.py  first.")
env = dict(os.environ, CONSANG_DATA=str(FRESH), CONSANG_OUT=str(OUT))
for s in ("02_model.py", "03_analyses.py"):
    subprocess.run([sys.executable, str(C.ROOT / "scripts" / s)], env=env, check=True, stdout=subprocess.DEVNULL)

A = json.load(open(C.ROOT / "results" / "key_numbers.json"))
B = json.load(open(OUT / "results" / "key_numbers.json"))
rows = [("Autosomal genes", A["n_genes"], B["n_genes"]),
        ("ClinVar P/LP variants (autosomal)", A["clinvar_plp_variants_autosomal"], B["clinvar_plp_variants_autosomal"])]
for a, b in zip(A["countries_all_groups"], B["countries_all_groups"]):
    rows += [(f"{a['country']}: risk per 100,000", a["per100k"], b["per100k"]),
             (f"{a['country']}: attributable fraction (%)", 100 * a["PAF"], 100 * b["PAF"])]
for anc in C.FOCUS:
    rows.append((f"ID/NDD first-cousin risk, {C.ANC_LABEL[anc]}, baseline", A["table2_ID_first_cousin"][anc]["Baseline"],
                 B["table2_ID_first_cousin"][anc]["Baseline"]))
for vs in A["nfe_to_sas_ratio_by_set"]:
    rows.append((f"European/South Asian ratio, {vs}", A["nfe_to_sas_ratio_by_set"][vs], B["nfe_to_sas_ratio_by_set"][vs]))
for vs, v in A["pakistan_births_2.5_50_97.5_by_set"].items():
    rows.append((f"Pakistan births per year, {vs}", v[1], B["pakistan_births_2.5_50_97.5_by_set"][vs][1]))
rows += [("Homozygous share of recessive diagnoses, Saudi cohort (%)", 100 * A["validation_homozygous"]["saudi_cohort_mix"],
          100 * B["validation_homozygous"]["saudi_cohort_mix"]),
         ("Newborn screening: pooled observed/expected", A["nbs_calibration"]["pooled_OE"], B["nbs_calibration"]["pooled_OE"]),
         ("Newborn screening: calibration slope", A["nbs_calibration"]["slope"], B["nbs_calibration"]["slope"]),
         ("South Asian lowest in every variant set", A["sas_lowest_in_every_set"], B["sas_lowest_in_every_set"]),
         ("Genes in every ancestry's top 25", ", ".join(A["top25_shared_all_four"]), ", ".join(B["top25_shared_all_four"]))]
df = pd.DataFrame(rows, columns=["quantity", "paper_release", "newer_release"])
num = pd.to_numeric(df.paper_release, errors="coerce"), pd.to_numeric(df.newer_release, errors="coerce")
df["change_pct"] = (100 * (num[1] / num[0] - 1)).round(2)
df.to_csv(C.SUPP / "TableS15_release_comparison.csv", index=False)
g = pd.read_csv(FRESH / "genes_panelapp.csv")
print("Newer release panels:", "; ".join(sorted({f"{p} v{v}" for p, v in zip(g.panel, g.panel_version)})))
print(df.to_string(index=False))
