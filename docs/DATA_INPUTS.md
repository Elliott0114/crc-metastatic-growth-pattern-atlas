# Analysis-ready inputs

The repository contains derived panel-level source data but does not redistribute raw or large public-source files. Place analysis-ready inputs at the paths below, relative to the repository root.

| Resource | Expected path(s) | Public source or accession | Used by |
|---|---|---|---|
| Targeted CRLM HGP ISS panel 2 | `data_sources/EscrivaConde_BJC_2026_HGP_ISS/prepared_cells_panel_2.h5ad` | Source article and panel-2 TMAP: `https://crclm2.serve.scilifelab.se/iss_panel_2.tmap?dl=1` | ISS ecology and NC3 mapping |
| Ogden CRLM multiome | `data_sources/Ogden_2025_CRLM_multiome/CRCLM_multiome_GEX_decontaminated.h5ad`; `data_sources/Ogden_2025_CRLM_multiome/supplementary_tables/mmc3.xlsx`, `mmc6.xlsx`, `mmc8.xlsx` | Source publication and supplementary material | NC3-to-REC mapping and programme definition |
| GSE151165 bulk RNA-seq | `data_sources/GSE151165/GSE151165_RNA_seq.raw_read_count.xlsx` | GEO GSE151165 | Independent bulk HGP projection and tumour–liver context |
| Latacz tumour–liver interface summary | `data_sources/Latacz_2024_targeted_HGP_interface/SupplTablesLataczforupload.xlsx` | Source-article supplementary tables | Interface gene-level corroboration |
| E-MTAB-12022 single-cell RNA-seq | `analysis_results/e_mtab_12022_primary_reconstruction.h5ad` | ArrayExpress/BioStudies E-MTAB-12022 | Epithelial patient-pseudobulk projection |
| E-MTAB-12043 FFPE spatial RNA | `analysis_work/e_mtab_12043_spatial_full/<sample>/outs/`; `data_sources/E_MTAB_HGP_CRLM/E-MTAB-12043_L1_spatial_pilot_download_manifest.tsv` | ArrayExpress/BioStudies E-MTAB-12043 | Regional and boundary projection |
| GSE294385 spatial RNA | `data_sources/Liu_2026_GSE294385/extracted/<sample>/filtered_feature_bc_matrix/`; `data_sources/Liu_2026_GSE294385/visium_liver_meta_after_qc.tsv.gz` | GEO GSE294385 | Paired micro-/macrometastatic analysis |

## Notes on source reconstruction

- `prepared_cells_panel_2.h5ad` is an analysis-ready representation of the public targeted ISS panel and must retain the observation fields checked by the scripts.
- `e_mtab_12022_primary_reconstruction.h5ad` is the tumour-only reconstruction used in the manuscript. It must retain patient, HGP, broad cell-compartment and count-layer information required by `analyze_rec_program_hgp_scrna_projection.py`.
- E-MTAB-12043 inputs are six study-specific Space Ranger output directories. The proprietary Space Ranger software and the large raw sequencing/image files are not redistributed here.
- GSE294385 sample-to-patient and region availability are fixed in `metadata/gse294385_liver_sample_manifest.tsv`.

Each analysis script validates the columns, dimensions or coverage needed for its endpoint and writes derived results under `analysis_results/`.
