# Hand-entered inputs from published sources

- `country_surveys.csv`: national consanguinity surveys. Proportions are of all marriages. For Qatar the survey reported the mean inbreeding coefficient directly (includes double first cousins).
- `nbs_saudi_mohamed2020.csv`: Mohamed et al. (2020) Saudi Med J 41:703-708. Newborn screening, Prince Sultan Military Medical City, Riyadh, 2012-2017, 56,632 infants screened. Genes column maps each disorder to its causal genes (MUT is the older symbol of MMUT).
- Validation cohort (constants in `scripts/common.py`): Monies et al. (2017) Hum Genet 136:921-939: 1,019 tests; 78% of families with information consanguineous; 97% of autosomal recessive diagnoses homozygous.
- Births: UN World Population Prospects 2024, medium-variant projection for Pakistan, 2025: 6,909,545.
- `riazuddin2017_variants.tsv`: causal variants from Riazuddin et al. (2017) Mol Psychiatry 22:1604-1614, Tables 1-3 (GRCh37), as published in the PMC full text (PMC5658665).
