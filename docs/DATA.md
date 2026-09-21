# Input Data

`metadata/analysis_inputs.tsv` is the file-level input contract: location,
repository/data-root path, staged path, bytes, SHA-256, source and consuming
modules. Downloaded filenames alone are not version identifiers. `check --mode
full` requires the recorded bytes and hashes and stops on drift.
`upstream_sha256` preserves the source hash when a bundled metadata file needs
a path-only portability edit; `adaptation` records that edit. The checked
`sha256` is always the actual distributed input, not a substituted source hash.

## Quick Route

Everything required is bundled in `data/frozen/` and `assets/`. No network,
original project, old workbook, unpublished manuscript, or external matrix is
needed. `reference/` is the independent validation baseline, not an input to
statistical estimation or workbook content assembly.

## Full Route

Use a separate directory, referred to as `INPUTS` below. Preserve each external
row's `path` beneath it. Bundled definitions and metadata are staged from the
repository automatically. Do not copy an entire old project into `INPUTS`.

1. Obtain the sources below and extract the particular files listed in the
   manifest. Extraction is required for archive members; do not substitute an
   archive for the listed matrix.
2. Extract `emtab-matrix-inputs-2026-09-21.tar.gz` into `INPUTS`. This companion
   is a separately delivered release asset, deliberately outside Git.
3. Run `python reproduce.py check --mode full --data-dir INPUTS`.
4. Only after the check passes, run `python reproduce.py full --data-dir INPUTS
   --output-dir NEW_OUTPUT`.

Paths in the manifest are relative. Historical `analysis_results/` names in
the **input** contract identify explicitly frozen programme definitions or
the reconstructed E-MTAB matrix, not permission to use fitted results.
Large inputs may be file symlinks in a local data directory; they are still
checked individually. The runner never links entire source directories.

## Source Inventory

| Resource | Fixed acquisition and preparation |
| --- | --- |
| CRC Atlas | [CELLxGENE collection](https://cellxgene.cziscience.com/collections/3a844375-60ea-474b-adf1-98e76928baee), collection version `f7281b78-3e3a-4c89-b2e7-96155853871b`, dataset version `4a8b9568-965e-46b8-a427-baab6bf018e5`, schema 7.1.0. Download the [versioned H5AD](https://datasets.cellxgene.cziscience.com/4a8b9568-965e-46b8-a427-baab6bf018e5.h5ad) as `data_sources/CRC_Atlas_CZI_core/crc_atlas_core_czi.h5ad`. Analyses use `raw.X`; Atlas is an acquisition container, not an independent validation study. |
| Ogden | [Mendeley Data v1](https://data.mendeley.com/datasets/yd7bb3shn5/1), DOI `10.17632/yd7bb3shn5.1`: obtain `CRCLM_multiome_GEX_decontaminated.h5ad`. Obtain mmc3/mmc6/mmc8 workbooks from the [source article](https://doi.org/10.1016/j.xgen.2025.100881) / [PMC supplementary archive](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12278653/supplementaryFiles). |
| ISS panel 2 | Download the [author's panel-2 archive](https://crclm2.serve.scilifelab.se/iss_panel_2.tmap?dl=1), extract `files/prepared_cells_panel_2.h5ad`, and place it at the manifest path. The [source paper](https://doi.org/10.1038/s41416-026-03567-y) defines the published states and HGP labels; no reclustering is substituted. |
| E-MTAB-12022 / 12043 | [Single-cell accession](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-12022) and [spatial accession](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-12043). Use the companion reconstructed matrices described below, not an unrelated processed portal export. |
| GSE151165 | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151165), `GSE151165_RNA_seq.raw_read_count.xlsx`; retain all original sample columns and their order. |
| GSE294385 | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE294385). Extract the sample directories containing `filtered_feature_bc_matrix/{matrix.mtx.gz,features.tsv.gz,barcodes.tsv.gz}`. Obtain the regional annotation from the authors and prepare it as below; it is not included in this repository or the companion asset. The bundled sample manifest retains the source pairing. |
| Latacz interface | Supplementary coefficient workbook from [the article](https://doi.org/10.1007/s10585-024-10319-w). The source does not supply a patient count matrix here: this module reassembles published coefficients and does **not** refit the authors' patient model. |
| GSE159216 | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159216), deposited series matrix and GPL17586 annotation. This is microarray expression, not RNA-seq counts. |
| Nissen proteomics | `mmc4.xlsx` from the [source supplementary archive](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12335997/supplementaryFiles). Start from the authors' normalized protein measurements, not raw mass spectra. |
| MRTX1133 | [GSE307774](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE307774), sample H5 matrices and matching deposited metadata. |
| Plexin B2 | [GSE267981](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE267981), deposited single-cell matrices; the source has one library per acute condition, so these contrasts remain descriptive. |
| GATA6 / serial selection | [GSE290753](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE290753) and [GSE290752](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE290752), deposited counts and series metadata. Preserve the Ensembl GRCm38 release-101 GTF and recorded orthologue tables. |
| Programmes and annotation | MSigDB human `2025.1.Hs`, author-defined HRC/RSC/iCMS/C8 signatures, source regulons and fixed orthologue/symbol tables. Exact files and hashes are in the input manifest; attribution is in `THIRD_PARTY_NOTICES.md`. |

### GSE294385 Author Annotation

The source annotation is excluded because its redistribution terms have not
been confirmed. Obtain `Meta_data/visium_meta_after_qc.rds` directly from the
authors via their [repository](https://github.com/yliuup/CRC_micromets_ST/tree/b6c40903eeda67d812712f2bb77064041db38e6c/Meta_data),
pinned to commit `b6c40903eeda67d812712f2bb77064041db38e6c`. The RDS SHA-256 is
`7c6d6d3fe32daf0136176fec53aef180e147b4c4372bd3dd5847adb83181c609`.
Use the provided original-code exporter without changing the source labels:

```bash
Rscript --vanilla tools/prepare_gse294385_annotation.R /absolute/path/to/visium_meta_after_qc.rds /absolute/path/to/INPUTS
python reproduce.py check --mode full --data-dir /absolute/path/to/INPUTS
```

The exporter retains liver rows, their order, Layer1/2/3 labels and barcodes,
and the original patient mapping from the bundled sample manifest. It writes
`data_sources/Liu_2026_GSE294385/visium_liver_meta_after_qc.tsv.gz` beneath
`INPUTS`. Expected SHA-256:
`cfcab046fe13bbc26cfb2646e1c646a28cd19493219b3746f64cefe7b52aa1bc`.
This is not newly inferred pathology. Quick reproduction uses derived result
tables and does not need the source annotation. Full reproduction stops if it
is missing or differs from the recorded version.

Full reproduction starts from assay-appropriate deposited matrices. The clinical,
proteomic and Latacz exceptions above are explicit; "full" does not imply
reconstruction from sequencing reads, mass spectra or unavailable patient data.

## E-MTAB Reconstruction

These matrices were reconstructed previously in this project; the public code
starts from them and does not rerun or package the FASTQ-processing pipeline.

- E-MTAB-12022: tumour R1/R2 libraries for PT44, PT54, PT59, PT52, PT55 and PT61;
  GENCODE v32 transcriptome; kallisto 0.52.0; bustools 0.45.1; 10x v3 chemistry;
  archived `3M-february-2018` barcode allowlist. Primary retention requires at
  least 150 detected genes and at most 70% mitochondrial UMIs. The matrix has
  38,425 cells and 59,368 genes. Broad-compartment annotations use target-disjoint
  Trebo T0 reference panels. This is not a reconstruction of the source authors'
  manual cluster/doublet filtering or EmptyDrops barcode calls.
- E-MTAB-12043: Space Ranger 4.1.0 / Martian v4.0.15; GRCh38-2020-A reference;
  `Visium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A`; templated ligation,
  filtered probes, read-2 length 50, image-based tissue detection. The companion
  retains the six filtered feature matrices, spot coordinates, scale factors
  and sample/HGP metadata. Original processing summaries accompany the asset.
- The retained filename `E-MTAB-12043_L1_spatial_pilot_download_manifest.tsv`
  is used only for the source sample-to-patient/HGP mapping. Its legacy pilot
  acquisition columns do not define the full reconstructed matrix inputs.
- The assays share four patients. Eligibility is reapplied by the selected
  analysis modules; six available patients do not mean six eligible epithelial
  pseudobulks.

`metadata/companion_files.tsv` and `metadata/release_assets.tsv` record the
individual files and compressed archive. The companion is distributed through
the [manuscript release](https://github.com/Elliott0114/crc-metastatic-growth-pattern-atlas/releases/tag/manuscript-2026-09-21),
outside Git. Verify the archive against its accompanying `.sha256` file before
extracting it. It contains no GSE294385 source annotation.

## Networking And Resources

The tested platform is Linux. In addition to Conda packages, the retained R
readers use `gzip`, `zcat`, `awk` and `sha256sum`; the full preflight checks that
these commands are available. On another machine, provide the GNU gzip/coreutils
and an awk implementation before running. No installer is invoked by the runner.

All analysis runs are offline. Downloads must use direct connections with
`HTTP_PROXY`, `HTTPS_PROXY`, `ALL_PROXY` and lowercase equivalents unset.
No script installs dependencies or changes proxy settings. The Atlas alone is
30,875,155,333 bytes; allow additional space for GSE294385, working matrices and
outputs. Full runs are substantially slower and more memory-intensive than
quick redraws; the validation report records measured times.
