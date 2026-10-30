"""
model.py — Step 2: static burden, database gap, screening panels and multigenerational dynamics.

Input : data/gene_freqs.csv   (from fetch_data.py)
Output: results/*.csv and figures/*.png  (one per planned manuscript figure)

    python scripts/02_model.py              # real data (run from the repository root)
    python scripts/02_model.py --synthetic  # pipeline test only; numbers are MEANINGLESS

X-linked genes are excluded here (see common.gene_freqs); the recessive model does not apply.

Model
-----
For a recessive gene with cumulative pathogenic allele frequency q and inbreeding coefficient F
(1/16 first cousins, 1/64 second cousins), the probability a birth is affected is

    r(q, F) = q^2 + F q (1 - q)

Risk is linear in F, so a population with a mix of marriage types needs only the mean
inbreeding coefficient F_bar. Over a disorder group (G genes):  R = 1 - prod_g (1 - r_g).
Attributable fraction (PAF) = (R(F_bar) - R(0)) / R(F_bar).

Dynamics (per gene, per generation), genotype frequencies under inbreeding:
    AA = p^2 + Fpq,  Aa = 2pq(1-F),  aa = q^2 + Fpq
aa genotypes are removed from reproduction with probability s (disease severity) and, once
screening starts, affected births are avoided with probability a = uptake x (gene on panel).
    q_sel = [pq(1-F) + (1-s_eff)(q^2+Fpq)] / [1 - s_eff(q^2+Fpq)],  s_eff = 1-(1-s)(1-a)
    q'    = q_sel + u (1 - q_sel)
u is calibrated so today's q is stationary under today's F (mutation-selection balance).
Purging falls out naturally: higher F exposes more aa, removing alleles faster.
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common  # noqa: E402

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

GROUPS = ["ID_NDD", "Epilepsy", "HSP", "Leukodystrophy"]
F_GRID = np.round(np.arange(0, 0.0626, 0.0025), 4)   # 0 .. 1/16
FOCUS_ANC = ["sas", "mid", "nfe", "afr"]
ANC_LABEL = {"sas": "South Asian", "mid": "Middle Eastern", "nfe": "European (non-Finnish)",
             "afr": "African/African American", "eas": "East Asian", "amr": "Admixed American"}
os.chdir(common.OUTROOT)
os.makedirs("results", exist_ok=True)
os.makedirs("figures", exist_ok=True)


# ------------------------------------------------------------------ data
def synthetic_data(seed=1):
    """Fake inputs with realistic shape, ONLY to test that the code runs."""
    rng = np.random.default_rng(seed)
    sizes = {"ID_NDD": 450, "Epilepsy": 180, "HSP": 60, "Leukodystrophy": 50}
    rows = []
    for grp, n in sizes.items():
        for i in range(n):
            base = rng.lognormal(np.log(3e-4), 1.4)
            for anc, (scale, an) in {"sas": (1.1, 60000), "mid": (1.2, 6000), "nfe": (1.0, 1_200_000),
                                     "afr": (0.9, 80000), "eas": (0.9, 50000), "amr": (1.0, 60000)}.items():
                q_true = min(base * rng.lognormal(np.log(scale), 0.8), 0.05)
                ac = rng.binomial(an, q_true)
                rows.append({"gene": f"{grp}_{i}", "disorder_group": grp, "ancestry": anc,
                             "q": ac / an, "ac_total": ac, "an_median": an})
    return pd.DataFrame(rows)


def load(synthetic):
    df = synthetic_data() if synthetic else common.gene_freqs(autosomal=True)
    missing = [g for g in GROUPS if g not in set(df.disorder_group)]
    if missing:
        print(f"WARNING: no genes for {missing} in gene_freqs.csv — re-run fetch_data.py")
    return df[df.disorder_group.isin(GROUPS)].copy()


# ------------------------------------------------------------------ core model
def risk(q, F):
    return q ** 2 + F * q * (1 - q)


def group_risk(q, F):
    return 1 - np.prod(1 - risk(np.asarray(q), F))


def posterior_q(q, an, rng, n, an_fill):
    """Beta(ac+0.5, an-ac+0.5) draws: frequency uncertainty, esp. for small samples (e.g. gnomAD 'mid').
    Genes with no P/LP variant in gnomAD (an == 0) are treated as 'observed zero times' in a sample
    of the ancestry's typical size (an_fill). This is how the reference-database gap enters the burden."""
    an = np.where(an > 0, an, an_fill).astype(float)
    k = q * an                      # effective allele count, so the posterior is centred on q
    return rng.beta(k[:, None] + 0.5, (an - k)[:, None] + 0.5, size=(len(q), n))


# ------------------------------------------------------------------ Fig 2: static burden
def static_burden(df, n_draws=500, seed=0):
    rng = np.random.default_rng(seed)
    out = []
    an_typ = df[df.an_median > 0].groupby("ancestry").an_median.median().to_dict()
    for (anc, grp), d in df.groupby(["ancestry", "disorder_group"]):
        draws = posterior_q(d.q.values, d.an_median.values, rng, n_draws, an_typ[anc])
        seen = d.q.values > 0
        for F in F_GRID:
            R_obs = group_risk(d.q.values, F)                       # observed alleles only
            R_uns = 1 - np.prod(1 - risk(draws[~seen], F), axis=0)  # genes with no observed allele
            R = 1 - np.prod(1 - risk(draws, F), axis=0)
            R0 = 1 - np.prod(1 - risk(draws, 0), axis=0)
            paf = (R - R0) / R
            out.append({"ancestry": anc, "disorder_group": grp, "F": F,
                        "per100k_observed_only": 1e5 * R_obs,
                        "per100k_unobserved_genes": 1e5 * np.median(R_uns),
                        "per100k": 1e5 * np.median(R), "per100k_lo": 1e5 * np.quantile(R, .025),
                        "per100k_hi": 1e5 * np.quantile(R, .975),
                        "PAF": np.median(paf), "PAF_lo": np.quantile(paf, .025), "PAF_hi": np.quantile(paf, .975)})
    res = pd.DataFrame(out)
    res.to_csv("results/static_burden.csv", index=False)

    fig, axes = plt.subplots(1, len(GROUPS), figsize=(16, 4), sharex=True)
    for ax, grp in zip(axes, GROUPS):
        for anc in FOCUS_ANC:
            s = res[(res.ancestry == anc) & (res.disorder_group == grp)]
            ax.plot(s.F, s.per100k, label=ANC_LABEL[anc])
            ax.fill_between(s.F, s.per100k_lo, s.per100k_hi, alpha=.15)
        for cF, cl in ((0.0079, "Turkey"), (0.0245, "Saudi Arabia"), (0.0271, "Qatar"), (0.0331, "Pakistan")):
            ax.axvline(cF, color="#999999", lw=0.6, ls=":")
            if grp == GROUPS[-1]:
                ax.text(cF + 0.0006, 0.97, cl, transform=ax.get_xaxis_transform(), fontsize=6.5, color="#555555",
                        rotation=90, ha="left", va="top")
        ax.set_title(grp); ax.set_xlabel("Mean inbreeding coefficient F")
    axes[0].set_ylabel("Births with biallelic pathogenic\ngenotypes per 100,000"); axes[0].legend(fontsize=8)
    fig.tight_layout(); fig.savefig("figures/Figure4_risk_by_inbreeding.png", dpi=200); plt.close(fig)
    return res


# ------------------------------------------------------------------ Fig 3: reference database gap
def database_gap(df):
    rows = []
    for (anc, grp), d in df.groupby(["ancestry", "disorder_group"]):
        rows.append({"ancestry": anc, "disorder_group": grp, "genes": len(d),
                     "genes_with_observed_PLP": int((d.ac_total > 0).sum()),
                     "frac_genes_observed": (d.ac_total > 0).mean(),
                     "median_AN": d.an_median.median(),
                     "risk_first_cousin_per100k": 1e5 * group_risk(d.q.values, 1 / 16)})
    gap = pd.DataFrame(rows)
    gap.to_csv("results/database_gap.csv", index=False)

    fig, ax = plt.subplots(figsize=(8, 4))
    piv = gap.pivot(index="disorder_group", columns="ancestry", values="frac_genes_observed").reindex(index=GROUPS, columns=FOCUS_ANC)
    piv.rename(columns=ANC_LABEL).plot.bar(ax=ax)
    ax.set_ylabel("Share of genes with ≥1 observed P/LP allele"); ax.set_xlabel("")
    ax.tick_params(axis="x", rotation=0); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig("figures/model_database_gap.png", dpi=200); plt.close(fig)
    return gap


# ------------------------------------------------------------------ Fig 6a: screening panel capture
def panel_capture(df, F=1 / 16):
    rows, curves = [], {}
    for anc in FOCUS_ANC:
        d = df[df.ancestry == anc].groupby("gene").q.first()
        r = np.sort(risk(d.values, F))[::-1]
        cum = np.cumsum(r) / r.sum() if r.sum() > 0 else np.zeros_like(r)
        curves[anc] = cum
        for target in (.5, .8, .9, .95):
            rows.append({"ancestry": anc, "F": F, "target_capture": target,
                         "genes_needed": int(np.searchsorted(cum, target) + 1), "genes_total": len(r)})
    res = pd.DataFrame(rows)
    res.to_csv("results/panel_capture.csv", index=False)

    fig, ax = plt.subplots(figsize=(6, 4))
    for anc, c in curves.items():
        ax.plot(np.arange(1, len(c) + 1), c, label=ANC_LABEL[anc])
    ax.set_xscale("log"); ax.set_xlabel("Genes on screening panel (ranked by risk)")
    ax.set_ylabel("Share of recessive neuro risk captured"); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig("figures/model_panel_capture.png", dpi=200); plt.close(fig)
    return res


# ------------------------------------------------------------------ Figs 5 & 6b: dynamics
U_CAP = 1e-5   # plausible upper limit for the per-gene rate of new pathogenic mutations per generation


def simulate(q0, F0, gens=8, s=0.9, F_decline=0.0, screen_start=None, uptake=0.0, panel=None, u_cap=U_CAP):
    """Returns incidence per 100k births per generation for one gene set (vector q0)."""
    q = np.asarray(q0, float).copy()
    p = 1 - q
    aa0 = q ** 2 + F0 * p * q
    q_sel0 = (p * q * (1 - F0) + (1 - s) * aa0) / (1 - s * aa0)
    u = np.clip((q - q_sel0) / (1 - q_sel0), 0, None)          # calibrated mutation input
    if u_cap is not None:
        u = np.minimum(u, u_cap)                                # genes above this are held by founder effects/drift, not mutation
    on_panel = np.ones_like(q, bool) if panel is None else panel
    inc = []
    for t in range(gens + 1):
        F = F0 * np.exp(-F_decline * t)
        p = 1 - q
        aa = q ** 2 + F * p * q
        a = uptake * on_panel if (screen_start is not None and t >= screen_start) else 0.0
        inc.append(1e5 * (1 - np.prod(1 - aa * (1 - a))))
        s_eff = 1 - (1 - s) * (1 - a)
        q_sel = (p * q * (1 - F) + (1 - s_eff) * aa) / (1 - s_eff * aa)
        q = q_sel + u * (1 - q_sel)
    return np.array(inc)


def dynamics(df, anc="sas", gens=8, s=0.9):
    d = df[df.ancestry == anc].groupby("gene").q.first()
    q0 = d.values
    order = np.argsort(risk(q0, 1 / 16))[::-1]
    panel_top = np.zeros(len(q0), bool)
    panel_top[order[:max(1, len(q0) // 10)]] = True               # top 10% of genes

    scen = {
        "Constant consanguinity (F=0.03)": dict(F0=0.03),
        "Consanguinity halves per ~2 generations": dict(F0=0.03, F_decline=0.35),
        "No consanguinity from now (F=0)": dict(F0=0.0),
        "Screening, full panel, 50% uptake": dict(F0=0.03, screen_start=1, uptake=.5),
        "Screening, top-10% panel, 50% uptake": dict(F0=0.03, screen_start=1, uptake=.5, panel=panel_top),
        "Decline + full-panel screening": dict(F0=0.03, F_decline=0.35, screen_start=1, uptake=.5),
    }
    # 'No consanguinity' needs calibration at today's F, then switch: emulate with instant decline
    scen["No consanguinity from now (F=0)"] = dict(F0=0.03, F_decline=50)
    rows = []
    for name, kw in scen.items():
        inc = simulate(q0, gens=gens, s=s, **kw)
        for t, v in enumerate(inc):
            rows.append({"ancestry": anc, "scenario": name, "generation": t, "years": 25 * t, "per100k": v})
    res = pd.DataFrame(rows)
    res.to_csv(f"results/dynamics_{anc}.csv", index=False)

    fig, ax = plt.subplots(figsize=(8, 4.5))
    for name, s_ in res.groupby("scenario", sort=False):
        ax.plot(s_.years, s_.per100k, marker="o", ms=3, label=name)
    ax.set_xlabel("Years from now (25-year generations)"); ax.set_ylabel("Births with biallelic pathogenic\ngenotypes per 100,000")
    ax.set_title(f"Recessive neurological disease, {ANC_LABEL[anc]} frequencies"); ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(f"figures/Figure5{'A' if anc == 'sas' else 'B'}_dynamics_{anc}.png", dpi=200); plt.close(fig)
    return res


def sensitivity_s(df, anc="sas"):
    """How the post-cessation rebound depends on disease severity s (selection against aa)."""
    q0 = df[df.ancestry == anc].groupby("gene").q.first().values
    rows = []
    for s in (0.5, 0.7, 0.9, 1.0):
        stop = simulate(q0, 0.03, gens=8, s=s, F_decline=50)
        decl = simulate(q0, 0.03, gens=8, s=s, F_decline=0.35)
        rows.append({"s": s, "gen0": stop[0], "stop_gen1": stop[1], "stop_gen8": stop[-1],
                     "stop_rebound_gen8_vs_gen1": stop[-1] / stop[1],
                     "decline_gen8_vs_gen0": decl[-1] / decl[0]})
    pd.DataFrame(rows).to_csv(f"results/sensitivity_selection_{anc}.csv", index=False)


# ------------------------------------------------------------------ main
if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--synthetic", action="store_true")
    args = ap.parse_args()
    df = load(args.synthetic)
    if args.synthetic:
        print("*** SYNTHETIC TEST DATA — results are not real ***")
    sb = static_burden(df)
    gap = database_gap(df)
    pc = panel_capture(df)
    for anc in ("sas", "mid"):
        dynamics(df, anc)
        sensitivity_s(df, anc)
    f1c = sb[np.isclose(sb.F, 0.0625)]
    print(f1c.pivot(index="disorder_group", columns="ancestry", values="per100k").round(1))
    print(pc[pc.target_capture == .9])
    print("Figures in figures/, tables in results/")
