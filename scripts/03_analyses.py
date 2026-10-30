"""
03_analyses.py — every analysis reported in the manuscript beyond 02_model.py.

Writes supplementary/S1..S9 tables and results/key_numbers.json (all numbers quoted in the text).
Run from the repository root after 02_model.py:
    python scripts/03_analyses.py
"""
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import binom, chi2, poisson

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402

rng = np.random.default_rng(1)
K = {}                                    # key numbers for the manuscript
gf = C.gene_freqs(autosomal=True)
genes = C.genes_table()
auto_genes = genes[~genes.is_x]
V = C.variants(autosomal=True)
an_typ = gf[gf.an_median > 0].groupby("ancestry").an_median.median().to_dict()
base = gf.drop_duplicates(["gene", "ancestry"])[["gene", "ancestry", "an_median"]]

K["n_genes"] = int(gf.gene.nunique())
K["genes_per_group"] = gf.drop_duplicates(["gene", "disorder_group"]).groupby("disorder_group").gene.nunique().to_dict()
cv = pd.read_csv(C.DATA / "clinvar_plp.csv", dtype=str)
K["clinvar_plp_variants_autosomal"] = int(cv[cv.HGNC_ID.isin(auto_genes.hgnc_id) | cv.GeneSymbol.isin(auto_genes.gene)].variant_id.nunique())
K["n_x_linked_excluded"] = len(C.x_linked_genes())


def gene_set(grp):
    return auto_genes.gene.unique() if grp == "All" else auto_genes[auto_genes.disorder_group == grp].gene.unique()


# ---------------------------------------------------------------- S4: four co-primary variant sets
rows, pak = [], {}
for name, sel in C.VARIANT_SETS.items():
    q = C.q_from_variants(V[sel(V)])
    for anc in C.ANCESTRIES:
        b = base[base.ancestry == anc]
        for grp in C.GROUPS + ["All"]:
            d = b[b.gene.isin(gene_set(grp))].merge(q[q.ancestry == anc][["gene", "q"]], on="gene", how="left").fillna({"q": 0})
            for F in (1 / 16, 0.0331, 0.03):
                R, R0 = C.group_risk(d.q, F), C.group_risk(d.q, 0)
                rows.append(dict(variant_set=name, ancestry=anc, group=grp, F=round(F, 4),
                                 per100k=round(1e5 * R, 1), PAF=round((R - R0) / R, 4)))
            if anc == "sas" and grp == "All":     # Pakistan: posterior draws at survey F
                an = np.where(d.an_median > 0, d.an_median, an_typ["sas"]).astype(float)
                k = d.q.values * an
                draws = rng.beta(k[:, None] + 0.5, (an - k)[:, None] + 0.5, size=(len(k), 1000))
                Rd = 1 - np.prod(1 - C.risk(draws, 0.494 / 16 + 0.142 / 64), axis=0)
                pak[name] = [int(round(C.PAKISTAN_BIRTHS * x)) for x in np.quantile(Rd, [0.025, 0.5, 0.975])]
S4 = pd.DataFrame(rows)
S4.to_csv(C.SUPP / "TableS4_variant_set_sensitivity.csv", index=False)
idf = S4[(S4.group == "ID_NDD") & (S4.F == 0.0625)]
K["table2_ID_first_cousin"] = idf.pivot(index="variant_set", columns="ancestry", values="per100k").to_dict()
K["nfe_to_sas_ratio_by_set"] = {s: round(d.set_index("ancestry").per100k["nfe"] / d.set_index("ancestry").per100k["sas"], 2)
                                for s, d in idf.groupby("variant_set")}
K["sas_lowest_in_every_set"] = bool(all(d.set_index("ancestry").per100k.idxmin() == "sas" for _, d in idf.groupby("variant_set")))
K["PAF_range_first_cousin_all_sets"] = [float(S4.PAF[S4.F == 0.0625].min()), float(S4.PAF[S4.F == 0.0625].max())]
K["pakistan_births_2.5_50_97.5_by_set"] = pak

# ---------------------------------------------------------------- validation 1: homozygosity (Monies 2017)
def homozygous_fraction(anc, F_cons, frac_cons):
    w = V[V.ancestry == anc].dropna(subset=["af"])
    g = w.groupby("gene").af.agg(q="sum", s2=lambda x: (x ** 2).sum())
    g["q"] = g.q.clip(upper=0.5)
    hom = lambda F: (g.s2 * (1 - F) + F * g.q).sum()
    aff = lambda F: (g.q ** 2 + F * g.q * (1 - g.q)).sum()
    return (frac_cons * hom(F_cons) + (1 - frac_cons) * hom(0)) / (frac_cons * aff(F_cons) + (1 - frac_cons) * aff(0))

F_cons = C.SAUDI_MEAN_F / C.SAUDI_CONSANG_RATE
K["validation_homozygous"] = {"outbred": homozygous_fraction("mid", 0, 0),
                              "saudi_cohort_mix": homozygous_fraction("mid", F_cons, C.MONIES_CONSANG_FRACTION),
                              "first_cousin": homozygous_fraction("mid", 1 / 16, 1),
                              "observed": C.MONIES_OBSERVED_HOMOZYGOUS, "F_consanguineous_couples": F_cons}

# ---------------------------------------------------------------- validation 2: newborn screening (S8)
nbs = pd.read_csv(C.EXT / "nbs_saudi_mohamed2020.csv")
mid = gf[gf.ancestry == "mid"].drop_duplicates("gene").set_index("gene")


def p_two_sided(o, e):
    if e == 0:
        return 1.0 if o == 0 else 0.0
    return min(1.0, 2 * min(poisson.cdf(o, e), poisson.sf(o - 1, e)))


out = []
for x in nbs.itertuples():
    have = [g for g in x.genes.split(";") if g in mid.index]
    p = 1 - np.prod([1 - C.risk(mid.loc[g, "q"], C.SAUDI_MEAN_F) for g in have]) if have else np.nan
    o = x.observed_cases
    lo = chi2.ppf(0.025, 2 * o) / 2 if o > 0 else 0.0
    hi = chi2.ppf(0.975, 2 * (o + 1)) / 2
    e = p * C.NBS_INFANTS if have else np.nan
    out.append(dict(disorder=x.disorder, genes_in_set=",".join(have) or "-", observed=o,
                    obs_per100k=1e5 * o / C.NBS_INFANTS, obs_lo=1e5 * lo / C.NBS_INFANTS, obs_hi=1e5 * hi / C.NBS_INFANTS,
                    pred_per100k=1e5 * p if have else np.nan, pred_cases=e,
                    exact_poisson_p=p_two_sided(o, e) if have else np.nan))
S8 = pd.DataFrame(out)
S8.round(4).to_csv(C.SUPP / "TableS8_newborn_screening_validation.csv", index=False)
m = S8.dropna(subset=["pred_cases"])
ex = m[m.disorder != "Phenylketonuria"]
K["nbs"] = {"n_disorders": len(m), "n_rejected_p05": int((m.exact_poisson_p < 0.05).sum()),
            "P_3plus_rejections_if_correct": float(1 - binom.cdf(2, len(m), 0.05)),
            "excl_PKU_observed": int(ex.observed.sum()), "excl_PKU_expected": float(ex.pred_cases.sum()),
            "excl_PKU_p": p_two_sided(int(ex.observed.sum()), float(ex.pred_cases.sum()))}

# ---------------------------------------------------------------- S1: countries
surv = pd.read_csv(C.EXT / "country_surveys.csv")
rows = []
for x in surv.itertuples():
    F = float(x.F_reported) if str(x.F_reported) not in ("", "nan") else x.first_cousin / 16 + float(x.other_consanguineous) / 64
    for grp in C.GROUPS + ["All"]:
        d = gf[gf.ancestry == x.ancestry] if grp == "All" else gf[(gf.ancestry == x.ancestry) & (gf.disorder_group == grp)]
        q = d.drop_duplicates("gene").q.values
        R, R0 = C.group_risk(q, F), C.group_risk(q, 0)
        rows.append(dict(country=x.country, ancestry=x.ancestry, mean_F=round(F, 4), group=grp,
                         per100k=round(1e5 * R, 1), PAF=round((R - R0) / R, 3), source=x.source))
S1 = pd.DataFrame(rows)
S1.to_csv(C.SUPP / "TableS1_country_burden.csv", index=False)
K["countries_all_groups"] = S1[S1.group == "All"][["country", "mean_F", "per100k", "PAF"]].to_dict("records")
qs = gf[gf.ancestry == "sas"].drop_duplicates("gene").q.values
K["pakistan_F_sensitivity_per100k"] = {lab: round(1e5 * C.group_risk(qs, F), 1) for lab, F in [
    ("other consanguineous as unrelated", 0.494 / 16), ("baseline", 0.494 / 16 + 0.142 / 64),
    ("other consanguineous as F=1/32", 0.494 / 16 + 0.142 / 32),
    ("10% of first cousins double first cousins", 0.494 * 0.9 / 16 + 0.494 * 0.1 / 8 + 0.142 / 64)]}

# ---------------------------------------------------------------- robustness of the attributable fraction (Pakistan)
F_pk = 0.494 / 16 + 0.142 / 64
K["PAF_if_frequencies_multiplied"] = {}
for k in (1, 2, 3, 5):
    R, R0 = C.group_risk(qs * k, F_pk), C.group_risk(qs * k, 0)
    K["PAF_if_frequencies_multiplied"][f"x{k}"] = {"per100k": round(1e5 * R, 1), "PAF": round((R - R0) / R, 4)}
rp, r0 = C.risk(qs, F_pk), qs ** 2
pafs = []
for _ in range(1000):
    pen = rng.uniform(0.2, 1, len(qs))
    R, R0 = 1 - np.prod(1 - pen * rp), 1 - np.prod(1 - pen * r0)
    pafs.append((R - R0) / R)
K["PAF_random_gene_penetrance_95pct"] = [float(x) for x in np.quantile(pafs, [0.025, 0.975])]

# ---------------------------------------------------------------- S5: dynamics under three histories x mutation models
def simulate(q0, F_now, F_cal, s, gens=8, decline=None, cap=1e-5, uptake=0.0, panel=None, cease=False):
    """Calibrate u at F_cal (optionally capped), then project from F_now. decline: rate of exp(-k t) decline."""
    q = q0.copy(); p = 1 - q
    aa = q ** 2 + F_cal * p * q
    qs = (p * q * (1 - F_cal) + (1 - s) * aa) / (1 - s * aa)
    u = np.clip((q - qs) / (1 - qs), 0, None)
    if cap:
        u = np.minimum(u, cap)
    on = np.ones_like(q, bool) if panel is None else panel
    inc = []
    for t in range(gens + 1):
        F = F_now if t == 0 else (0.0 if cease else (F_now if decline is None else F_now * np.exp(-decline * t)))
        a = uptake * on if t >= 1 else 0.0
        p = 1 - q; aa = q ** 2 + F * p * q
        inc.append(1e5 * (1 - np.prod(1 - aa * (1 - a))))
        se = 1 - (1 - s) * (1 - a)
        qs = (p * q * (1 - F) + (1 - se) * aa) / (1 - se * aa)
        q = qs + u * (1 - qs)
    return np.array(inc)


HISTORIES = [(0.0, "Above equilibrium (recent rise in consanguinity)"), (0.03, "At equilibrium (baseline)"),
             (0.06, "Below equilibrium (higher past consanguinity)")]
MUT_MODELS = [(1e-5, "Primary: u capped at 1e-5"), (None, "Sensitivity: uncapped"), (3e-6, "Additional sensitivity: u capped at 3e-6")]
rows = []
for anc in ["sas", "mid"]:
    q0 = gf[gf.ancestry == anc].drop_duplicates("gene").q.values
    top = np.zeros(len(q0), bool); top[np.argsort(C.risk(q0, 1 / 16))[::-1][:max(1, len(q0) // 10)]] = True
    for cap, caplab in MUT_MODELS:
        for Fc, lab in HISTORIES:
            for s in (0.5, 0.7, 0.9, 1.0):
                kw = dict(F_cal=Fc, s=s, cap=cap)
                const = simulate(q0, 0.03, **kw); dec = simulate(q0, 0.03, decline=0.35, **kw)
                ces = simulate(q0, 0.03, cease=True, **kw)
                full = simulate(q0, 0.03, uptake=.5, **kw); t10 = simulate(q0, 0.03, uptake=.5, panel=top, **kw)
                ds = simulate(q0, 0.03, decline=0.35, uptake=.5, **kw)
                rows.append(dict(ancestry=anc, mutation_model=caplab, history=lab, F_calibration=Fc, s=s,
                                 risk_now=round(const[0], 1), constant_F_gen8=round(const[-1], 1),
                                 declining_F_gen8=round(dec[-1], 1), cessation_gen1=round(ces[1], 1), cessation_gen8=round(ces[-1], 1),
                                 full_screening_gen8=round(full[-1], 1), top10_screening_gen8=round(t10[-1], 1),
                                 decline_plus_screening_gen8=round(ds[-1], 1),
                                 rebound_pct=round(100 * (ces[-1] / ces[1] - 1), 1),
                                 gradual_decline_reduction_pct=round(100 * (1 - dec[-1] / dec[0]), 1)))
S5 = pd.DataFrame(rows)
S5.to_csv(C.SUPP / "TableS5_dynamics_sensitivity.csv", index=False)
K["decline_reduction_range_pct"] = [float(S5.gradual_decline_reduction_pct.min()), float(S5.gradual_decline_reduction_pct.max())]
K["rebound_range_pct"] = [float(S5.rebound_pct.min()), float(S5.rebound_pct.max())]
cp = S5[S5.mutation_model.str.startswith("Primary")]
K["capped_decline_reduction_range_pct"] = [float(cp.gradual_decline_reduction_pct.min()), float(cp.gradual_decline_reduction_pct.max())]
K["capped_cessation_change_range_pct"] = [float(cp.rebound_pct.min()), float(cp.rebound_pct.max())]

# ---------------------------------------------------------------- S7: top genes; panel yield; screening effort
rows = []
for anc in C.FOCUS:
    d = gf[gf.ancestry == anc].drop_duplicates("gene").copy()
    d["risk"] = C.risk(d.q, 1 / 16); d["share"] = d.risk / d.risk.sum()
    for i, x in enumerate(d.sort_values("risk", ascending=False).head(25).itertuples(), 1):
        rows.append(dict(ancestry=C.ANC_LABEL[anc], rank=i, gene=x.gene, q=round(x.q, 5),
                         first_cousin_risk_per100k=round(1e5 * x.risk, 1), share_pct=round(100 * x.share, 2)))
S7 = pd.DataFrame(rows)
S7.to_csv(C.SUPP / "TableS7_top_genes.csv", index=False)
tops = {a: set(d.gene) for a, d in S7.groupby("ancestry")}
K["top25_shared_all_four"] = sorted(set.intersection(*tops.values()))
K["top25_overlap_sas_mid"] = len(tops["South Asian"] & tops["Middle Eastern"])

d = gf[gf.ancestry == "sas"].drop_duplicates("gene").copy()
Fpk = 0.494 / 16 + 0.142 / 64
d["r"] = C.risk(d.q, Fpk); d = d.sort_values("r", ascending=False)
K["pakistan_births_per_gene_by_rank_block"] = {f"{a + 1}-{b}": round(C.PAKISTAN_BIRTHS * d.r.iloc[a:b].sum() / (b - a))
                                               for a, b in [(0, 10), (10, 50), (50, 100), (100, 300), (300, len(d))]}
effort = {}
for lab, F in [("first_cousin_couples", 1 / 16), ("all_couples_pakistan", Fpk)]:
    x = np.sort(C.risk(d.q.values, F))[::-1]
    effort[lab] = {k: round(1 / (1 - np.prod(1 - x[:k]))) for k in (10, 50, 100, 300, len(x))}
K["couples_screened_per_affected_pregnancy"] = effort

# ---------------------------------------------------------------- S9: ascertainment versus founder genes
rows = []
for k in (0, 10, 25, 50):
    top = (gf[gf.ancestry == "nfe"].drop_duplicates("gene").assign(r=lambda x: C.risk(x.q, 1 / 16))
           .sort_values("r", ascending=False).gene.head(k))
    row = {"top_european_genes_removed": k}
    for anc in C.ANCESTRIES:
        dd = gf[(gf.ancestry == anc) & (gf.disorder_group == "ID_NDD") & (~gf.gene.isin(top))]
        row[f"{anc}_ID_NDD_per100k"] = round(1e5 * C.group_risk(dd.q.values, 1 / 16), 1)
    row["nfe_to_sas_ratio"] = round(row["nfe_ID_NDD_per100k"] / row["sas_ID_NDD_per100k"], 2)
    rows.append(row)
pd.DataFrame(rows).to_csv(C.SUPP / "TableS9_founder_gene_removal.csv", index=False)
wide = gf.drop_duplicates(["gene", "ancestry"]).pivot(index="gene", columns="ancestry", values="q")
wide.to_csv(C.SUPP / "TableS9b_gene_by_gene_frequencies.csv")
both = wide[(wide.nfe > 0) & (wide.sas > 0)]
K["gene_by_gene"] = {"n_both": len(both), "share_nfe_higher": float((both.nfe > both.sas).mean()),
                     "median_ratio": float((both.nfe / both.sas).median()),
                     "only_nfe": int(((wide.nfe > 0) & (wide.sas == 0)).sum()),
                     "only_sas": int(((wide.sas > 0) & (wide.nfe == 0)).sum())}

# ---------------------------------------------------------------- other supplementary tables
sb = pd.read_csv(C.RESULTS / "static_burden.csv")
sb.to_csv(C.SUPP / "TableS2_burden_by_ancestry_and_F.csv", index=False)
f1 = sb[np.isclose(sb.F, 0.0625)]
K["unobserved_share_mid"] = (f1[f1.ancestry == "mid"].set_index("disorder_group")
                             .eval("per100k_unobserved_genes / per100k").round(3).to_dict())
pd.read_csv(C.DATA / "excluded_common_variants.csv").to_csv(C.SUPP / "TableS3_excluded_common_variants.csv", index=False)
auto_genes.drop(columns="is_x").to_csv(C.SUPP / "TableS6_gene_list_panel_versions.csv", index=False)
gap = pd.read_csv(C.RESULTS / "database_gap.csv")
K["share_genes_observed_range"] = gap.groupby("ancestry").frac_genes_observed.agg(["min", "max"]).round(3).to_dict("index")

# ---------------------------------------------------------------- X-linked order-of-magnitude estimate
gx = C.gene_freqs(autosomal=False)
gx = gx[gx.gene.isin(C.x_linked_genes())].drop_duplicates(["gene", "ancestry"])
K["x_linked_affected_per100k_births"] = (gx.groupby("ancestry").q.sum() * 1e5 / 2).round(1).to_dict()

# ---------------------------------------------------------------- S11: sample-size equalisation (rarefaction)
ID = set(auto_genes[auto_genes.disorder_group == "ID_NDD"].gene)


def id_summary(df_q):
    q = df_q[df_q.gene.isin(ID)].groupby("gene").q.sum().clip(upper=0.5).reindex(sorted(ID)).fillna(0).values
    return 1e5 * C.group_risk(q, 1 / 16), float((q > 0).mean())


rows = []
for anc in C.FOCUS:
    w = V[(V.ancestry == anc) & (V.an > 0)].copy(); w["q"] = w.ac / w.an
    b, c = id_summary(w)
    rows.append(dict(dataset=f"{C.ANC_LABEL[anc]} (observed)", iteration=0, ID_NDD_per100k=b, share_genes_observed=c))
nfe = V[(V.ancestry == "nfe") & (V.an > 0)][["gene", "ac", "an"]]
rs = np.random.default_rng(7)
for tgt in ("sas", "mid"):
    N = int(V[(V.ancestry == tgt) & (V.an > 0)].groupby("variant_id").an.max().median())
    for it in range(300):
        n = np.minimum(nfe.an.values, N)
        ac = rs.hypergeometric(nfe.ac.values, (nfe.an - nfe.ac).values, n)
        b, c = id_summary(pd.DataFrame({"gene": nfe.gene.values, "q": ac / n}))
        rows.append(dict(dataset=f"European downsampled to {C.ANC_LABEL[tgt]} size", iteration=it + 1, ID_NDD_per100k=b, share_genes_observed=c))
S11 = pd.DataFrame(rows); S11.to_csv(C.SUPP / "TableS11_sample_size_equalisation.csv", index=False)
K["downsampling"] = S11.groupby("dataset").agg(burden_median=("ID_NDD_per100k", "median"), burden_lo=("ID_NDD_per100k", lambda x: x.quantile(.025)),
                                               burden_hi=("ID_NDD_per100k", lambda x: x.quantile(.975)), share_median=("share_genes_observed", "median"),
                                               share_lo=("share_genes_observed", lambda x: x.quantile(.025)),
                                               share_hi=("share_genes_observed", lambda x: x.quantile(.975))).round(3).to_dict("index")

# ---------------------------------------------------------------- S12: priors for genes with no observed allele
PRIORS = {"Zero contribution": None, "Beta(0.5,0.5) (baseline)": (0.5, 0.5), "Beta(1,1)": (1, 1), "Conservative Beta(0.1,1)": (0.1, 1)}
rows = []
rp2 = np.random.default_rng(3)
for name, pr in PRIORS.items():
    for anc in C.FOCUS:
        for grp in ["ID_NDD", "All"]:
            d = (gf[gf.ancestry == anc] if grp == "All" else gf[(gf.ancestry == anc) & (gf.disorder_group == grp)]).drop_duplicates("gene")
            an = np.where(d.an_median > 0, d.an_median, an_typ[anc]).astype(float); q = d.q.values; seen = q > 0
            if pr is None:
                draws = np.repeat(q[:, None], 300, axis=1)
            else:
                k = q * an
                draws = rp2.beta(np.where(seen, k + 0.5, pr[0])[:, None], np.where(seen, an - k + 0.5, an + pr[1])[:, None], size=(len(q), 300))
            R = 1 - np.prod(1 - C.risk(draws, 1 / 16), axis=0); R0 = 1 - np.prod(1 - C.risk(draws, 0), axis=0)
            rows.append(dict(prior=name, ancestry=anc, group=grp, per100k=float(np.median(R)) * 1e5, PAF=float(np.median((R - R0) / R))))
S12 = pd.DataFrame(rows); S12.round(4).to_csv(C.SUPP / "TableS12_unobserved_gene_priors.csv", index=False)
K["priors_ID_first_cousin"] = S12[S12.group == "ID_NDD"].pivot(index="prior", columns="ancestry", values="per100k").round(0).to_dict()
K["priors_PAF_range_all_groups"] = [float(S12[S12.group == "All"].PAF.min()), float(S12[S12.group == "All"].PAF.max())]

# ---------------------------------------------------------------- newborn screening calibration
import statsmodels.api as sm  # noqa: E402
m8 = S8.dropna(subset=["pred_cases"]).copy()
ci = [(o / e, (chi2.ppf(.025, 2 * o) / 2 if o > 0 else 0) / e, chi2.ppf(.975, 2 * (o + 1)) / 2 / e) if e > 0 else (np.nan,) * 3
      for o, e in zip(m8.observed, m8.pred_cases)]
m8[["obs_to_exp", "obs_to_exp_lo", "obs_to_exp_hi"]] = ci
m8.round(4).to_csv(C.SUPP / "TableS8_newborn_screening_validation.csv", index=False)
pos = m8[m8.pred_cases > 0]
glm = sm.GLM(pos.observed.values, sm.add_constant(np.log(pos.pred_cases.values)), family=sm.families.Poisson()).fit()
O, E = int(m8.observed.sum()), float(m8.pred_cases.sum())
K["nbs_calibration"] = {"pooled_obs": O, "pooled_exp": E, "pooled_OE": O / E,
                        "pooled_OE_CI": [chi2.ppf(.025, 2 * O) / 2 / E, chi2.ppf(.975, 2 * (O + 1)) / 2 / E],
                        "slope": float(glm.params[1]), "slope_CI": [float(x) for x in glm.conf_int()[1]],
                        "intercept": float(glm.params[0]), "intercept_CI": [float(x) for x in glm.conf_int()[0]],
                        "within_2fold": int(((pos.obs_to_exp >= .5) & (pos.obs_to_exp <= 2)).sum()), "n_with_expected": len(pos)}
nc = K["nbs_calibration"]
summ = pd.DataFrame([
    {"disorder": "Pooled (all 12 disorders)", "observed": O, "pred_cases": round(E, 2), "obs_to_exp": nc["pooled_OE"],
     "obs_to_exp_lo": nc["pooled_OE_CI"][0], "obs_to_exp_hi": nc["pooled_OE_CI"][1]},
    {"disorder": "Calibration model log E(Y) = alpha + beta log(Yhat): intercept alpha (95% CI)",
     "obs_to_exp": nc["intercept"], "obs_to_exp_lo": nc["intercept_CI"][0], "obs_to_exp_hi": nc["intercept_CI"][1]},
    {"disorder": "Calibration model: slope beta (95% CI); beta = 1 is correct relative calibration",
     "obs_to_exp": nc["slope"], "obs_to_exp_lo": nc["slope_CI"][0], "obs_to_exp_hi": nc["slope_CI"][1]}])
pd.concat([m8, summ], ignore_index=True).round(4).to_csv(C.SUPP / "TableS8_newborn_screening_validation.csv", index=False)

# ---------------------------------------------------------------- S13: plausibility of calibrated mutation input
rows = []
for anc in ["sas", "mid"]:
    q = gf[gf.ancestry == anc].drop_duplicates("gene").set_index("gene").q
    q = q[q > 0]; p = 1 - q; aa = q ** 2 + 0.03 * p * q
    qs = (p * q * 0.97 + 0.1 * aa) / (1 - 0.9 * aa); u = ((q - qs) / (1 - qs)).clip(lower=0)
    for g in q.index:
        rows.append(dict(ancestry=anc, gene=g, q=q[g], u_calibrated=u[g]))
S13 = pd.DataFrame(rows); S13.to_csv(C.SUPP / "TableS13_calibrated_mutation_input.csv", index=False)
K["u_calibrated"] = {}
for anc, d in S13.groupby("ancestry"):
    r_ = C.risk(d.q, 1 / 16)
    K["u_calibrated"][anc] = {"median": float(d.u_calibrated.median()), "p90": float(d.u_calibrated.quantile(.9)),
                              "share_genes_u_gt_1e-5": float((d.u_calibrated > 1e-5).mean()),
                              "share_risk_u_gt_1e-5": float(r_[d.u_calibrated > 1e-5].sum() / r_.sum())}


def simulate_capped(q0, cap, decline=None, s=0.9, gens=8):
    q = q0.copy(); p = 1 - q; aa = q ** 2 + 0.03 * p * q
    qs = (p * q * 0.97 + (1 - s) * aa) / (1 - s * aa); u = np.clip((q - qs) / (1 - qs), 0, None)
    if cap:
        u = np.minimum(u, cap)
    inc = []
    for t in range(gens + 1):
        F = 0.03 if t == 0 else (0.0 if decline is None else 0.03 * np.exp(-decline * t))
        p = 1 - q; aa = q ** 2 + F * p * q; inc.append(1e5 * (1 - np.prod(1 - aa)))
        qs = (p * q * (1 - F) + (1 - s) * aa) / (1 - s * aa); q = qs + u * (1 - qs)
    return np.array(inc)


K["capped_mutation_dynamics"] = {}
for anc in ["sas", "mid"]:
    q0 = gf[gf.ancestry == anc].drop_duplicates("gene").q.values
    for cap in (None, 1e-5, 3e-6):
        x, y = simulate_capped(q0, cap), simulate_capped(q0, cap, decline=0.35)
        K["capped_mutation_dynamics"][f"{anc}_cap_{cap}"] = {"cessation_change_pct": 100 * (x[-1] / x[1] - 1),
                                                             "gradual_decline_reduction_pct": 100 * (1 - y[-1] / y[0])}

json.dump(K, open(C.RESULTS / "key_numbers.json", "w"), indent=1, default=float)
print(json.dumps({k: K[k] for k in ["n_genes", "validation_homozygous", "nbs", "pakistan_births_2.5_50_97.5_by_set",
                                    "decline_reduction_range_pct", "rebound_range_pct", "nfe_to_sas_ratio_by_set"]},
                 indent=1, default=float))
