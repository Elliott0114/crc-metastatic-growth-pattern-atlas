# Regenerative-like epithelial programme in colorectal liver metastases

This repository contains the core analysis code and panel-level source data for the manuscript:

> **Cross-platform mapping links a regenerative-like epithelial programme to replacement growth and lesion size in colorectal liver metastases**

The study integrates targeted in situ sequencing, single-cell RNA sequencing, bulk RNA sequencing and spatial transcriptomics to connect the NC3 neoplastic state with a source-defined regenerative epithelial cell (REC) programme and to test that programme across histopathological growth-pattern and lesion-size contexts.

## Repository scope

This is a reader-facing release of the analyses used in the current manuscript. It includes:

- patient-level analysis of neoplastic-state composition and tumour–liver proximity;
- cross-platform NC3-to-REC state mapping;
- source-only derivation of the fixed 50-gene REC programme;
- projection into GSE151165, the Latacz interface study, E-MTAB-12022, E-MTAB-12043 and GSE294385;
- cross-context gene-level and REC–hypoxia component analyses;
- the fixed REC-programme table, analysis specifications and dataset-role metadata;
- panel-level source data for all six main figures and five supplementary figures.

Only analyses supporting the current manuscript are retained; superseded exploratory and figure-development branches are excluded.

## Repository layout

| Path | Contents |
|---|---|
| `scripts/` | Core Python and R analysis scripts |
| `metadata/` | Source-specific analysis specifications, dataset roles and the GSE294385 sample manifest |
| `analysis_results/ogden_anchor_program/selected_anchor_program_top50.tsv` | Fixed 50-gene REC-programme definition used for projection |
| `source_data/figure_source_data/` | Panel-level values underlying the displayed figures |
| `source_data/supplementary_tables/` | Supplementary Tables 1 and 2 in TSV and XLSX formats |
| `environment.yml` | Conda environment used for the reported analyses |

## Reproducibility boundary

Raw public datasets are not redistributed. Several source studies provide large or study-specific objects, so this release separates data acquisition/reconstruction from the statistical analysis layer:

- the scripts in `scripts/` implement the reported patient-level analyses from the expected analysis-ready inputs;
- the exact expected input paths and source accessions are documented in [`docs/DATA_INPUTS.md`](docs/DATA_INPUTS.md);
- the fixed programme and panel-level source data are included directly so the analysed gene set and all displayed quantitative values can be inspected without downloading the raw datasets.

The repository does not claim a single turnkey raw-FASTQ-to-manuscript workflow across all platforms. Source-specific preprocessing must reproduce the analysis-ready files described in `docs/DATA_INPUTS.md` before the corresponding analysis script is run.

## Environment

Create the recorded environment from the repository root:

```bash
conda env create -f environment.yml
```

Run Python and R scripts with:

```bash
conda run -n crc-metastatic-growth python scripts/<script.py>
conda run -n crc-metastatic-growth Rscript scripts/<script.R>
```

All scripts resolve paths relative to the repository root. Random procedures use seeds recorded in the scripts and analysis specifications.

## Computational dependency order

The order below reflects file dependencies, not the inferential priority assigned to each dataset in the manuscript.

1. `scripts/analyze_hgp_interface_ecology.py`
2. `scripts/analyze_cross_platform_state_anchor.py`
3. `scripts/analyze_ogden_anchor_program.R`
4. `scripts/analyze_rec_program_hgp_scrna_projection.py`
5. `scripts/analyze_rec_program_hgp_spatial_projection.py`
6. `scripts/analyze_rec_program_cross_modal_concordance.py`
7. `scripts/analyze_gse294385_rec_program_extension.py`
8. `scripts/analyze_rec_program_two_axis_interpretation.py`
9. `scripts/analyze_rec_program_gse151165_bulk_projection.R`
10. `scripts/analyze_rec_program_gse151165_tissue_context.R`
11. `scripts/analyze_rec_program_latacz_interface_summary.R`
12. `scripts/analyze_rec_program_multisource_core.R`
13. `scripts/analyze_rec_program_reader_facing_specificity_R.R`
14. `scripts/prepare_rec_program_spatial_display.py`

The primary GSE151165 endpoint uses the source-derived 50-gene programme. Its position in this computational order only allows the same script to append secondary gene subsets created by later cross-context summaries; it does not alter the frozen primary score.

## Statistical unit

Patients are the inferential units throughout. Cells, spots, regions of interest and computational resamples are not treated as independent biological replicates. The source-specific specifications in `metadata/` define eligibility rules, endpoints and interpretation boundaries.

## Data access

The principal public accessions are GEO **GSE151165** and **GSE294385**, and ArrayExpress/BioStudies **E-MTAB-12022** and **E-MTAB-12043**. The targeted ISS, Ogden multiomic and Latacz interface resources are identified in `docs/DATA_INPUTS.md` and in the manuscript Data availability statement.

## Citation

Please cite this repository using [`CITATION.cff`](CITATION.cff). A versioned snapshot is available from the GitHub Releases page.

## Contact

Correspondence about the manuscript and code: Sheng Dai, Department of Colorectal Surgery, Sir Run Run Shaw Hospital, Zhejiang University School of Medicine.
