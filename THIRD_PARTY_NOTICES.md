# Third-Party Sources And Terms

The MIT license covers original code in this repository only. It does not
relicense third-party matrices, author annotations, gene sets, histology images,
software packages or the manuscript. No third-party software environment is
redistributed. Retain the source citations and original terms when reusing data.

| Material | Attribution and terms |
| --- | --- |
| Ogden processed matrix | Sam Ogden; Mendeley Data version 1, DOI `10.17632/yd7bb3shn5.1`; [CC BY 4.0 stated on the deposit](https://data.mendeley.com/datasets/yd7bb3shn5/1). The large matrix is downloaded separately. |
| MSigDB 2025.1.Hs Hallmark, GO BP and Reactome sets | Copyright Broad Institute, Massachusetts Institute of Technology and Regents of the University of California; [MSigDB CC BY 4.0 terms](https://www.gsea-msigdb.org/gsea/msigdb_license_terms.jsp). Fixed subsets/format changes retain source identity. This export does not add KEGG or BioCarta collections. |
| Reactome pathway definitions | [Reactome data are CC0](https://reactome.org/license). The MSigDB-distributed representation also retains the MSigDB attribution above. |
| HRC, RSC, iCMS, REC/C8 and regulon definitions | Author-defined signatures from the papers and source supplementary tables recorded in the definition registries. They remain attributed source material, not MIT code. Clinical/functional source workbooks are not republished as raw files here. |
| CRC Atlas, ISS, GEO, ArrayExpress, Ensembl and UCSC inputs | Obtain from the original deposits listed in `docs/DATA.md` and the file manifest. Public availability is not treated as an MIT grant. Consult the relevant deposit and original source terms; accession records and local hashes identify the exact data used. |
| E-MTAB reconstructed matrices and H&E images | Derived from the Fleischer et al. public E-MTAB-12022/12043 deposits. The [source article's rights section](https://doi.org/10.1186/s12943-023-01713-1) states CC BY 4.0 for article material and CC0 for article data, subject to individual credit-line exceptions. The companion and two retained source H&E views preserve provenance; check the deposited material's applicable terms rather than treating the code license as permission. |
| GSE294385 regional annotations | Not redistributed in this repository or its companion asset. Full reproduction requires the source `Meta_data/visium_meta_after_qc.rds` from the [author repository](https://github.com/yliuup/CRC_micromets_ST), commit `b6c40903eeda67d812712f2bb77064041db38e6c`; see `docs/DATA.md` for preparation and checksums. No explicit redistribution license was found there. The original-code exporter does not relicense the source annotation. Derived analysis results and the sample-level pairing manifest remain available for reproduction. |
| Figure 1 schematic artwork | AI-generated/edited with OpenAI imagegen in the original project; source notes are under `assets/` and `docs/IMAGE_SOURCES.md`. This artwork is illustrative, not a data-bearing image. No image model/version was exposed by the original interface. |
| Analysis dependencies | Installed separately under their respective licenses, including Python/R, NumPy, pandas, SciPy, AnnData, Matplotlib, edgeR, limma and metafor. `environment.yml` and `metadata/environment_versions.json` record the tested versions. |

The full repository contains derived scientific tables needed for redraws and
small frozen definitions/annotations. It excludes the source articles, private
submission materials, correspondence, credentials and large original matrices.
No statement here expands permission beyond a source's own terms. Before
reuse, readers should check the applicable terms for the companion matrices
and retained source H&E views with their
institutional requirements. This record is provenance, not a legal determination.
