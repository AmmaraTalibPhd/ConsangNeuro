#!/usr/bin/env bash
# Reproduce every result from the deposited data (script 05 downloads one small file).
# To re-download the inputs first, run: python scripts/01_fetch_data.py
set -e
python scripts/02_model.py
python scripts/03_analyses.py
python scripts/04_figures.py
python scripts/05_literature_variants.py   # downloads a small chain file once
# pLoF analysis (Table S14): needs data/plof_variant_counts.csv,
# created once with: python scripts/06_plof_analysis.py --fetch
if [ -f data/plof_variant_counts.csv ]; then python scripts/06_plof_analysis.py; fi
# Release check (Table S15): needs data/fresh/ from python scripts/01_fetch_data.py
if [ -f data/fresh/gene_freqs.csv ]; then python scripts/07_release_check.py; fi
