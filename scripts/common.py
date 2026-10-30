"""Shared paths, constants and helpers used by scripts 02-04."""
import os
import pathlib

import numpy as np
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parent.parent
# CONSANG_DATA / CONSANG_OUT let 07_release_check.py rerun the pipeline on a newer data release
DATA = pathlib.Path(os.environ.get("CONSANG_DATA", ROOT / "data"))
EXT = ROOT / "data" / "external"
OUTROOT = pathlib.Path(os.environ.get("CONSANG_OUT", ROOT))
RESULTS, FIGURES, SUPP = OUTROOT / "results", OUTROOT / "figures", OUTROOT / "supplementary"
for d in (RESULTS, FIGURES, SUPP):
    d.mkdir(parents=True, exist_ok=True)

GROUPS = ["ID_NDD", "Epilepsy", "HSP", "Leukodystrophy"]
ANCESTRIES = ["sas", "mid", "nfe", "afr", "eas", "amr"]
FOCUS = ["sas", "mid", "nfe", "afr"]
ANC_LABEL = {"sas": "South Asian", "mid": "Middle Eastern", "nfe": "European", "afr": "African",
             "eas": "East Asian", "amr": "Admixed American"}
ANC_COLOR = {"sas": "#1f77b4", "mid": "#ff7f0e", "nfe": "#2ca02c", "afr": "#d62728"}

MAX_VARIANT_AF = 0.01                  # P/LP variants above this AF in any ancestry are excluded
TWO_STAR = {"criteria provided, multiple submitters, no conflicts",
            "reviewed by expert panel", "practice guideline"}

# Validation cohort: Monies et al. 2017 (Saudi Arabia)
MONIES_CONSANG_FRACTION = 0.78
MONIES_OBSERVED_HOMOZYGOUS = 0.97
SAUDI_MEAN_F = 0.336 / 16 + 0.224 / 64      # El Mouzan 2008 survey
SAUDI_CONSANG_RATE = 0.56
NBS_INFANTS = 56632                          # Mohamed et al. 2020
PAKISTAN_BIRTHS = 6_909_545                  # UN WPP 2024, medium variant, 2025


def risk(q, F):
    """Probability that a birth is affected: q^2 + F q (1 - q)."""
    return q ** 2 + F * q * (1 - q)


def group_risk(q, F):
    return 1 - np.prod(1 - risk(np.asarray(q, float), F))


def genes_table():
    g = pd.read_csv(DATA / "genes_panelapp.csv")
    g["is_x"] = g["chrom"].astype(str).isin(["X", "Y"])
    return g


def x_linked_genes():
    g = genes_table()
    return set(g.loc[g.is_x, "gene"])


def gene_freqs(autosomal=True):
    """Gene x ancestry cumulative P/LP frequencies (one row per gene, ancestry and disorder group)."""
    df = pd.read_csv(DATA / "gene_freqs.csv")
    df = df[df.disorder_group.isin(GROUPS)]
    if autosomal:
        df = df[~df.gene.isin(x_linked_genes())]
    return df


def variants(autosomal=True):
    """Variant x ancestry allele counts joined to ClinVar annotations, common variants removed."""
    v = pd.read_csv(DATA / "variant_counts.csv")
    cv = pd.read_csv(DATA / "clinvar_plp.csv", dtype=str,
                     usecols=["variant_id", "ClinicalSignificance", "ReviewStatus"])
    v = v.merge(cv, on="variant_id", how="left")
    v["af"] = v.ac / v.an.where(v.an > 0)
    max_af = v.groupby("variant_id").af.max()
    v["max_af"] = v.variant_id.map(max_af)
    v = v[v.max_af <= MAX_VARIANT_AF]
    if autosomal:
        v = v[~v.gene.isin(x_linked_genes())]
    return v


VARIANT_SETS = {
    "Baseline": lambda v: v.max_af <= MAX_VARIANT_AF,
    "AF<=0.5%": lambda v: v.max_af <= 0.005,
    ">=2 stars": lambda v: v.ReviewStatus.isin(TWO_STAR),
    "Pathogenic only": lambda v: v.ClinicalSignificance == "Pathogenic",
}


def q_from_variants(v):
    """Cumulative frequency per gene and ancestry from a variant table."""
    return (v.dropna(subset=["af"]).groupby(["gene", "ancestry"]).af.sum()
             .clip(upper=0.5).rename("q").reset_index())
