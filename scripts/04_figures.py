"""
04_figures.py — publication figures 1, 2, 3, 6 and Figure S1.
(Figures 4 and 5 are written by 02_model.py.)  Run after 02 and 03:
    python scripts/04_figures.py
"""
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import FancyBboxPatch  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                     "axes.spines.top": False, "axes.spines.right": False})
K = json.load(open(C.RESULTS / "key_numbers.json"))
FOUR = [(a, C.ANC_LABEL[a]) for a in C.FOCUS]

# ---------------------------------------------------------------- Fig. 1 framework
fig, ax = plt.subplots(figsize=(10, 4.6)); ax.set_xlim(0, 100); ax.set_ylim(0, 46); ax.axis("off")


def box(x, y, w, h, title, body, fc):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.4,rounding_size=1.2", fc=fc, ec="#444", lw=0.8))
    ax.text(x + w / 2, y + h - 2.2, title, ha="center", va="top", fontsize=9.5, fontweight="bold")
    ax.text(x + w / 2, y + h - 6.2, body, ha="center", va="top", fontsize=8, linespacing=1.35)


def arrow(x1, y1, x2, y2):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1), arrowprops=dict(arrowstyle="-|>", color="#444", lw=1.1))


gp = K["genes_per_group"]
ax.text(1, 44.5, "a", fontsize=12, fontweight="bold"); ax.text(1, 20, "b", fontsize=12, fontweight="bold")
box(3, 26, 19, 15, "PanelApp", f"{K['n_genes']} autosomal recessive genes\nID/NDD {gp['ID_NDD']} · Epilepsy {gp['Epilepsy']}\n"
    f"Leukodystrophy {gp['Leukodystrophy']} · HSP {gp['HSP']}", "#e8f0fa")
box(27, 26, 19, 15, "ClinVar", f"{K['clinvar_plp_variants_autosomal']:,} P/LP variants\n≥1 review star\nAF ≤1% in every ancestry", "#e8f0fa")
box(51, 26, 19, 15, "gnomAD v4", "Allele counts in 6 ancestries\nSouth Asian · Middle Eastern\nEuropean · African · others", "#e8f0fa")
box(75, 26, 22, 15, "Gene frequency q", "q = Σ pathogenic allele freqs\nper gene and ancestry\n+ Beta uncertainty", "#fdf0e2")
for a, b in [(22.8, 27), (46.8, 51), (70.8, 75)]:
    arrow(a, 33.5, b, 33.5)
box(3, 2, 27, 15, "Static risk", "r(q,F) = q² + F q(1−q)\ngenotype-defined risk and\nattributable fraction by country", "#eef7ea")
box(36, 2, 27, 15, "Clinical comparison and data gap", "% homozygous (Saudi cohort)\nnewborn-screening incidence\nsample size, classification, pLoF", "#eef7ea")
box(69, 2, 28, 15, "Multigenerational model", "Selection s, purging, mutation u\nconsanguinity trends\ncarrier screening scenarios", "#eef7ea")
for x in (16.5, 49.5, 83):
    arrow(86, 25.4, x, 17.6)
fig.tight_layout(); fig.savefig(C.FIGURES / "Figure1_framework.png", dpi=300); plt.close(fig)

# ---------------------------------------------------------------- Fig. 2 validation
v = K["validation_homozygous"]
fig, axs = plt.subplots(1, 2, figsize=(10, 4.4), gridspec_kw=dict(width_ratios=[1, 1.4]))
vals = [100 * v["outbred"], 100 * v["saudi_cohort_mix"], 100 * v["first_cousin"], 100 * v["observed"]]
labs = ["Outbred\n(model)", "Saudi cohort\nmix (model)", "First cousins\n(model)", "Saudi cohort\n(observed)"]
bars = axs[0].bar(range(4), vals, color=["#bbbbbb", "#ff7f0e", "#ffbb78", "#444444"], width=0.6)
for b, x in zip(bars, vals):
    axs[0].text(b.get_x() + b.get_width() / 2, x + 1.5, f"{x:.1f}%" if x != 97 else "97%", ha="center", fontsize=8.5)
axs[0].set_xticks(range(4), labs, fontsize=8); axs[0].set_ylim(0, 110)
axs[0].axhline(97, ls="--", lw=0.8, color="#444444")
axs[0].set_ylabel("Recessive diagnoses that are homozygous (%)"); axs[0].set_title("a", loc="left", fontweight="bold")
t = pd.read_csv(C.SUPP / "TableS8_newborn_screening_validation.csv").dropna(subset=["pred_per100k"]).sort_values("pred_per100k")
ax = axs[1]
for i, x in enumerate(t.itertuples()):
    ax.plot([x.obs_lo, x.obs_hi], [i, i], color="#444", lw=1.2); ax.plot(x.obs_per100k, i, "o", color="#444", ms=5)
    ax.plot(x.pred_per100k, i, "D", color="#d62728" if x.disorder == "Phenylketonuria" else "#ff7f0e", ms=6)
names = [d.replace(" deficiency", "").replace("Glutaric aciduria type 1", "Glutaric aciduria 1")
         .replace("Maple syrup urine disease", "MSUD") for d in t.disorder]
ax.set_yticks(range(len(t)), names, fontsize=8)
ax.set_xscale("symlog", linthresh=1); ax.set_xlim(0, 100); ax.set_xticks([0, 1, 10, 100], ["0", "1", "10", "100"])
ax.set_xlabel("Incidence per 100,000 births")
ax.plot([], [], "o", color="#444", label="Observed, newborn screening (95% CI)"); ax.plot([], [], "D", color="#ff7f0e", label="Predicted by model")
ax.legend(frameon=False, fontsize=8, loc="lower center", bbox_to_anchor=(0.55, 1.04), ncol=2)
ax.set_title("b", loc="left", fontweight="bold", pad=18)
fig.tight_layout(); fig.savefig(C.FIGURES / "Figure2_validation.png", dpi=300); plt.close(fig)

# ---------------------------------------------------------------- Fig. 3 reference-data bias
S4 = pd.read_csv(C.SUPP / "TableS4_variant_set_sensitivity.csv")
S4 = S4[(S4.group == "ID_NDD") & np.isclose(S4.F, 0.0625)]
sets = ["Baseline", "AF<=0.5%", ">=2 stars", "Pathogenic only"]
set_labels = ["Baseline", "AF ≤0.5%", "≥2 review\nstars", "Pathogenic\nonly"]
gap = pd.read_csv(C.RESULTS / "database_gap.csv")
sb = pd.read_csv(C.RESULTS / "static_burden.csv"); sb = sb[np.isclose(sb.F, 0.0625)]
w = 0.2
fig, axs = plt.subplots(1, 3, figsize=(11, 3.8))
for i, (a, lab) in enumerate(FOUR):
    axs[0].bar(np.arange(4) + (i - 1.5) * w, [S4[(S4.variant_set == s) & (S4.ancestry == a)].per100k.iloc[0] for s in sets],
               w, label=lab, color=C.ANC_COLOR[a])
    axs[1].bar(np.arange(4) + (i - 1.5) * w, [100 * gap[(gap.ancestry == a) & (gap.disorder_group == g)].frac_genes_observed.iloc[0]
                                            for g in C.GROUPS], w, color=C.ANC_COLOR[a])
    axs[2].bar(np.arange(4) + (i - 1.5) * w, [100 * sb[(sb.ancestry == a) & (sb.disorder_group == g)]
                                            .eval("per100k_unobserved_genes / per100k").iloc[0] for g in C.GROUPS],
               w, color=C.ANC_COLOR[a])
axs[0].set_xticks(range(4), set_labels); axs[0].set_ylabel("Predicted ID/NDD risk per 100,000\nfirst-cousin births")
for ax in axs[1:]:
    ax.set_xticks(range(4), ["ID/NDD", "Epilepsy", "HSP", "Leuko."])
axs[1].set_ylabel("Genes with ≥1 P/LP allele in gnomAD (%)"); axs[1].set_ylim(0, 100)
axs[2].set_ylabel("Share of predicted risk from\nunobserved genes (%)")
for ax, l in zip(axs, "abc"):
    ax.set_title(l, loc="left", fontweight="bold")
h, l = axs[0].get_legend_handles_labels(); fig.legend(h, l, loc="lower center", ncol=4, frameon=False)
fig.tight_layout(rect=(0, 0.07, 1, 1)); fig.savefig(C.FIGURES / "Figure3_reference_gap.png", dpi=300); plt.close(fig)

# ---------------------------------------------------------------- Fig. 6 panel capture
gf = C.gene_freqs()
fig, ax = plt.subplots(figsize=(6.5, 4))
for a, lab in FOUR:
    x = np.sort(C.risk(gf[gf.ancestry == a].drop_duplicates("gene").q.values, 1 / 16))[::-1]
    ax.plot(np.arange(1, len(x) + 1), 100 * np.cumsum(x) / x.sum(), color=C.ANC_COLOR[a], label=lab)
for y in (50, 90):
    ax.axhline(y, ls="--", color="#999", lw=0.8); ax.text(1.05, y + 1.5, f"{y}%", color="#666", fontsize=8)
ax.set_xscale("log"); ax.set_ylim(0, 102)
ax.set_xlabel("Genes on screening panel (ranked by contribution to risk)"); ax.set_ylabel("Share of first-cousin risk captured (%)")
ax.legend(frameon=False, loc="upper left", fontsize=8, bbox_to_anchor=(0.02, 0.85))
fig.tight_layout(); fig.savefig(C.FIGURES / "Figure6_panel_capture.png", dpi=300); plt.close(fig)

# ---------------------------------------------------------------- Figure S1 population history


def sim(q0, F_now, F_cal, s=0.9, gens=8, decline=None, cap=1e-5):
    q = q0.copy(); p = 1 - q; aa = q ** 2 + F_cal * p * q
    qs = (p * q * (1 - F_cal) + (1 - s) * aa) / (1 - s * aa); u = np.clip((q - qs) / (1 - qs), 0, None); inc = []
    if cap:
        u = np.minimum(u, cap)
    for t_ in range(gens + 1):
        F = F_now if t_ == 0 else (0.0 if decline is None else F_now * np.exp(-decline * t_))
        p = 1 - q; aa = q ** 2 + F * p * q; inc.append(1e5 * (1 - np.prod(1 - aa)))
        qs = (p * q * (1 - F) + (1 - s) * aa) / (1 - s * aa); q = qs + u * (1 - qs)
    return np.array(inc)


H = [(0.0, "Above equilibrium (recent rise in consanguinity)", "#2ca02c"), (0.03, "At equilibrium (baseline)", "#1f77b4"),
     (0.06, "Below equilibrium (higher past consanguinity)", "#d62728")]
fig, axs = plt.subplots(2, 2, figsize=(10, 7.6)); yrs = np.arange(9) * 25
for j, a in enumerate(["sas", "mid"]):
    q0 = gf[gf.ancestry == a].drop_duplicates("gene").q.values
    for Fc, lab, col in H:
        axs[0, j].plot(yrs[1:], sim(q0, 0.03, Fc)[1:], "o-", color=col, ms=3, label=lab + ", u capped at 1e-5")
        axs[0, j].plot(yrs[1:], sim(q0, 0.03, Fc, cap=None)[1:], "o--", color=col, ms=3, mfc="white", label=lab + ", uncapped")
        axs[1, j].plot(yrs, sim(q0, 0.03, Fc, decline=0.35), "o-", color=col, ms=3)
        axs[1, j].plot(yrs, sim(q0, 0.03, Fc, decline=0.35, cap=None), "o--", color=col, ms=3, mfc="white")
    axs[0, j].set_title(f"{C.ANC_LABEL[a]}: consanguinity ceases", fontsize=9)
    axs[1, j].set_title(f"{C.ANC_LABEL[a]}: consanguinity halves every ~2 generations", fontsize=9)
    for i in range(2):
        axs[i, j].set_xlabel("Years from now"); axs[i, j].set_ylabel("Births with biallelic pathogenic\ngenotypes per 100,000")
for ax, l in zip(axs.flat, "abcd"):
    ax.text(-0.12, 1.05, l, transform=ax.transAxes, fontweight="bold", fontsize=11)
h, l = axs[0, 0].get_legend_handles_labels(); fig.legend(h, l, loc="lower center", ncol=2, frameon=False, fontsize=7.5)
fig.tight_layout(rect=(0, 0.11, 1, 1)); fig.savefig(C.FIGURES / "FigureS1_population_history.png", dpi=300); plt.close(fig)
print("Figures written to", C.FIGURES)
