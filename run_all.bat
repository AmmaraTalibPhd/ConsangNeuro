@echo off
REM Reproduce every result from the deposited data (no internet needed).
REM To re-download the inputs first, run: python scripts\01_fetch_data.py
python scripts\02_model.py || exit /b 1
python scripts\03_analyses.py || exit /b 1
python scripts\04_figures.py || exit /b 1
python scripts\05_literature_variants.py || exit /b 1
REM pLoF analysis (Table S14): needs data\plof_variant_counts.csv, created once with: python scripts\06_plof_analysis.py --fetch
if exist data\plof_variant_counts.csv python scripts\06_plof_analysis.py || exit /b 1
REM Release check (Table S15): needs data\fresh\ from python scripts\01_fetch_data.py
if exist data\fresh\gene_freqs.csv python scripts\07_release_check.py || exit /b 1
