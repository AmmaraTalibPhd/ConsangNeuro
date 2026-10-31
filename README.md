# Consanguinity and recessive neurological disease
[![DOI](https://zenodo.org/badge/1440281643.svg)](https://doi.org/10.5281/zenodo.23822689)

Code and data for:

> Talib, A., Aşır, S., and Farman, M. **Genomic reference under-representation obscures recessive neurological disease burden in consanguineous populations.** Submitted to *The American Journal of Human Genetics* (2026). [DOI to add on publication]

Corresponding author: Ammara Talib (ammara.talib@neu.edu.tr), The Research Center for Science Technology and Engineering, Department of Biomedical Engineering, Near East University, Mersin 10, Türkiye.

The model estimates genotype-defined risk of autosomal recessive neurological disease (intellectual disability and neurodevelopmental disorders, epilepsy, hereditary spastic paraplegia and leukodystrophy) as a function of parental inbreeding, using 948 diagnostic-grade genes from Genomics England PanelApp, ClinVar pathogenic variants and ancestry-specific gnomAD v4.1 allele frequencies. It then asks whether current reference data are sufficient to quantify this burden: it separates the effects of sample size (downsampling) and variant classification (a comparison using high-confidence predicted loss-of-function variants selected independently of clinical classification), compares predictions with two clinical datasets, measures how many causal variants from Pakistani families the model would count, and projects risk over generations under changing consanguinity and carrier screening.

## Quick start: reproduce every result (about a minute)

```bash
pip install -r requirements.txt
bash run_all.sh            # macOS / Linux
run_all.bat                # Windows
```

This uses the deposited data in `data/` (script 05 downloads one small chain file from GitHub) and writes:

- `results/`: model outputs, plus `key_numbers.json` and `plof_key_numbers.json` with every number quoted in the manuscript
- `supplementary/`: Tables S1–S15 and Note S1
- `figures/`: Figures 1–6 and Figures S1–S2


## Files in the paper and the scripts that make them

Names follow the article: main figures are Figures 1–6 and main tables are Tables 1–2 (in the manuscript file); supplemental items are Tables S1–S15, Figures S1–S2 and Note S1.

| Item in the article | File | Made by |
|---|---|---|
| Figure 1 | `figures/Figure1_framework.png` | `04_figures.py` |
| Figure 2 | `figures/Figure2_validation.png` | `04_figures.py` |
| Figure 3 | `figures/Figure3_reference_gap.png` | `04_figures.py` |
| Figure 4 | `figures/Figure4_risk_by_inbreeding.png` | `02_model.py` |
| Figure 5A, 5B | `figures/Figure5A_dynamics_sas.png`, `figures/Figure5B_dynamics_mid.png` | `02_model.py` |
| Figure 6 | `figures/Figure6_panel_capture.png` | `04_figures.py` |
| Figure S1 | `figures/FigureS1_population_history.png` | `04_figures.py` |
| Figure S2 | `figures/FigureS2_plof.png` | `06_plof_analysis.py` |
| Table S1 | `supplementary/TableS1_country_burden.csv` | `03_analyses.py` |
| Table S2 | `supplementary/TableS2_burden_by_ancestry_and_F.csv` | `03_analyses.py` |
| Table S3 | `supplementary/TableS3_excluded_common_variants.csv` | `03_analyses.py` |
| Table S4 | `supplementary/TableS4_variant_set_sensitivity.csv` | `03_analyses.py` |
| Table S5 | `supplementary/TableS5_dynamics_sensitivity.csv` | `03_analyses.py` |
| Table S6 | `supplementary/TableS6_gene_list_panel_versions.csv` | `03_analyses.py` |
| Table S7 | `supplementary/TableS7_top_genes.csv` | `03_analyses.py` |
| Table S8 | `supplementary/TableS8_newborn_screening_validation.csv` | `03_analyses.py` |
| Table S9 | `supplementary/TableS9_founder_gene_removal.csv`, `TableS9b_gene_by_gene_frequencies.csv` | `03_analyses.py` |
| Table S10 | `supplementary/TableS10_pakistani_family_variants.csv` | `05_literature_variants.py` |
| Table S11 | `supplementary/TableS11_sample_size_equalisation.csv` | `03_analyses.py` |
| Table S12 | `supplementary/TableS12_unobserved_gene_priors.csv` | `03_analyses.py` |
| Table S13 | `supplementary/TableS13_calibrated_mutation_input.csv` | `03_analyses.py` |
| Table S14 | `supplementary/TableS14_plof_ancestry_comparison.csv`, `TableS14b_plof_gene_frequencies.csv`, `TableS14c_plof_sensitivity.csv` | `06_plof_analysis.py` |
| Table S15 | `supplementary/TableS15_release_comparison.csv` | `07_release_check.py` |
| Note S1 | `supplementary/NoteS1_X_linked.md` (exploratory X-linked estimate) | text; numbers from `03_analyses.py` |

## Repository layout

| Path | Contents |
|---|---|
| `scripts/01_fetch_data.py` | Downloads gene panels (PanelApp), pathogenic variants (ClinVar) and allele frequencies (gnomAD v4.1). Needs internet; 1–2 hours. **Not needed to reproduce the paper.** |
| `scripts/02_model.py` | Core model: risk by ancestry and inbreeding with uncertainty, reference-data gap, screening panels, multigenerational projections (mutation input capped at 10⁻⁵). Writes Figures 4 and 5. |
| `scripts/03_analyses.py` | Four co-primary variant sets, Pakistan estimate, comparison with clinical data, country estimates, population-history and mutation-cap sensitivity, gene-level tables, founder-gene removal, X-linked estimate, downsampling, priors for unobserved genes, newborn-screening calibration, calibrated mutation input. Writes Tables S1–S9 and S11–S13 and `key_numbers.json`. |
| `scripts/04_figures.py` | Figures 1, 2, 3, 6 and Figure S1. |
| `scripts/05_literature_variants.py` | Checks causal variants from published Pakistani families against the gene list, ClinVar and gnomAD (Table S10). Downloads a GRCh37→GRCh38 chain file once. |
| `scripts/06_plof_analysis.py` | Classification-independent comparison (Table S14, Figure S2) using gnomAD high-confidence predicted loss-of-function variants. `--fetch` downloads the LOFTEE-annotated variants once (internet, about 1.5–2 hours, resumable) into `data/plof_variant_counts.csv`; without `--fetch` it runs the analysis in seconds. |
| `scripts/07_release_check.py` | Reruns scripts 02–03 on the newer PanelApp/ClinVar download in `data/fresh/` and compares headline results with the primary analysis (Table S15). |
| `scripts/common.py` | Shared constants and helpers. |
| `data/` | The frozen inputs of the primary analysis (see below). |
| `data/fresh/` | The October 2026 PanelApp/ClinVar/gnomAD download used only for the release-robustness check (Table S15). |
| `data/external/` | Values hand-entered from published sources (consanguinity surveys, newborn screening, Pakistani family variants), with citations. |

## Data

| File | Description |
|---|---|
| `data/genes_panelapp.csv` | Green, biallelic genes per disorder group, with panel names and versions, HGNC and Ensembl IDs and chromosome. X-linked genes are kept here but excluded from the recessive model. |
| `data/clinvar_plp.csv` | ClinVar Pathogenic / Likely pathogenic variants (GRCh38, ≥1 review star) in these genes. |
| `data/variant_counts.csv` | gnomAD v4.1 allele counts (exome + genome) per variant and ancestry. |
| `data/gene_freqs.csv` | Cumulative pathogenic allele frequency per gene and ancestry. |
| `data/excluded_common_variants.csv` | P/LP variants excluded for allele frequency >1% in any ancestry. |
| `data/plof_variant_counts.csv` | gnomAD v4.1 predicted loss-of-function variants (consequence, LOFTEE annotation, filters) with allele counts per ancestry. |

**Primary analysis inputs (frozen):** PanelApp panels as listed in `data/genes_panelapp.csv` (Intellectual disability v11.31); ClinVar `variant_summary.txt.gz` downloaded September 2026; gnomAD v4.1 (GraphQL API, accessed September 2026). Genes are retrieved from gnomAD by Ensembl ID, or by gene symbol where gnomAD uses a different Ensembl ID (for example *TBCE*).

ClinVar and PanelApp are updated monthly, so a new download gives slightly different inputs; `07_release_check.py` shows that the conclusions do not change (Table S15). To reproduce the published numbers exactly, use the deposited `data/` files. The raw downloads (`variant_summary.txt.gz` and the per-gene gnomAD responses in `gnomad_cache/` and `gnomad_lof_cache/`) are not in this repository; `01_fetch_data.py` and `06_plof_analysis.py --fetch` regenerate them.

## Re-running from scratch (optional)

```bash
python scripts/01_fetch_data.py                       # writes to data/fresh/, leaves the paper's inputs untouched
python scripts/01_fetch_data.py --replace-deposited   # only if you deliberately want to replace data/
python scripts/06_plof_analysis.py --fetch            # re-download the loss-of-function variants
```

All downloads are resumable: re-run the same command if the connection drops or gnomAD is busy.

## Licences

Code: MIT (see `LICENSE`). Input data remain under the terms of their sources: Genomics England PanelApp, ClinVar (NCBI, public domain) and gnomAD (Broad Institute).

## Citation

Please cite the article above and the archived release: Talib, A., Aşır, S., and Farman, M. ConsangNeuro. Zenodo. https://doi.org/10.5281/zenodo.23822690

Code: https://github.com/AmmaraTalibPhd/ConsangNeuro
