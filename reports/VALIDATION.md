# Validation Report

Baseline: Biology Direct manuscript draft, 21 September 2026.
Tested with the unchanged `crc-metastatic-growth` environment on Linux;
exact Python/R package versions are in `metadata/environment_versions.json`.

| Route | Actual completion | Result |
| --- | --- | --- |
| Quick | 11 figures, 52 panels / 58 panel source tables, 8 workbooks / 127 sheets | 204 artifact checks and 144 manuscript checks passed |
| Full | 61 declared modules from matrix/metadata inputs, then all figures and workbooks | 205 artifact/completion checks and 144 manuscript checks passed |

Both CLI tests used relocated repositories and output paths containing spaces.
Quick received no external inputs. Full staged only the checksummed file allowlist,
including the independently extracted E-MTAB companion. File-access traces found
no undeclared original-project access. No frozen-result fallback was used.
The final verifier additionally rechecked the clean full output from a relocated
copy after integer-count and figure-index checks were tightened. Scientific and
plotting scripts were hash-identical to those used by the complete cold run.

Full clean CLI elapsed time: 114.9 minutes.
See `module_validation.tsv` for measured module times. Data acquisition time is
not included. The cold run used 180 files, including 37 bundled definitions
and annotations. For public distribution, the byte-identical GSE294385 source
annotation is external: the published contract has 36 bundled and 144 external
inputs. Scientific scripts, input hashes and statistical baselines are unchanged.

## Comparisons

PNG dimensions and RGBA pixels agree exactly with the manuscript baseline.
Panel data and workbook cells check identifiers, membership order, missingness
and values; floats use `atol=1e-8`, `rtol=1e-6`. Manuscript numbers additionally
check their reported rounding. The source indices are rebuilt with current file
hashes and checked for complete membership, within-figure panel order and source
order within worksheets. Figure groups in the index follow the current redraw
order instead of the manuscript's historical editing order. PDF/ZIP timestamps and
compressed-file byte differences are not treated as statistical differences.

The development replay also compared 417 intermediate
tables with the retained project results. 18 whole-file
differences are explained in `module_difference_explanations.json`; these are
unused historical side outputs, omitted unused columns, portable image paths or
runtime hashes. They do not relax the current manuscript artifact checks.
Existing internal reconstruction assertions remain in the scientific scripts.

## Historical REC44 Definition

Supplementary Table 4 / `REC44_cross_context` retains the pre-correction spatial
column whose epithelial markers included KRT20. The isolated historical branch
recomputed all 44 retained values from the matrices. Corrected main HGP analyses
remain target-disjoint and unchanged. This inherited table must not be presented
as using the corrected marker definition; see `docs/WORKFLOW.md`.

## Failure And Scope Checks

Six unit tests and four real CLI failure tests passed: missing bundled input,
checksum mismatch, missing full inputs without fallback, and preservation of an
existing output directory. The companion's individual file hashes were verified
after extraction. R libraries resolved from the target Conda environment; no
packages were installed or changed and `renv::restore()` was not used.

Full starts from assay-appropriate matrices. It does not rerun FASTQ processing,
raw mass-spectrometry processing, or unavailable Latacz patient-level models.
Latacz coefficients and fixed signature/annotation definitions retain their
original-source role. This is manuscript reproduction, not independent validation
of study design, pathological labels or clinical claims. Third-party material
retains its original terms. Publication checks are recorded separately in
`publication_validation.json`; they do not represent a second cold full run.
