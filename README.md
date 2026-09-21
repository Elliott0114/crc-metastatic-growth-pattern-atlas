# CRC Regenerative and Junction Programmes

Analysis code for **Claudin and polarity-complex transcription differ in their
associations with regenerative expression in colorectal liver metastases**,
Biology Direct manuscript draft, 21 September 2026.

This export targets five main figures, six supplementary figures (52 panels),
and eight supplementary workbooks (127 worksheets). It does not typeset the
manuscript or reproduce unused exploratory analyses. Original scientific
scripts retain their statistical methods, random seeds, patient eligibility,
gene definitions, and multiple-testing families.

## Run

Use the versions in `environment.yml`. This repository uses Conda, not renv.
Create a separate environment on a new machine:

```bash
conda env create -f environment.yml
conda activate crc-metastatic-growth
python reproduce.py check --mode quick
python reproduce.py quick --output-dir /absolute/path/to/quick-run
python reproduce.py verify --output-dir /absolute/path/to/quick-run
```

Quick reproduction uses bundled frozen numerical results and image assets;
it does **not** refit models. Every invocation uses a new output directory.

For count-matrix analysis, obtain the exact inputs described in
[DATA.md](docs/DATA.md), including the separately distributed E-MTAB companion:

```bash
python reproduce.py check --mode full --data-dir /absolute/path/to/inputs
python reproduce.py full --data-dir /absolute/path/to/inputs --output-dir /absolute/path/to/full-run
python reproduce.py verify --output-dir /absolute/path/to/full-run
```

The full route must not substitute frozen results for missing analysis inputs.
It starts at deposited or reconstructed matrices and metadata, not FASTQ.
Only paths declared in the input manifest are staged. Historical date-stamped
paths inside a run are stable internal identifiers, not dependencies on old
manuscript directories outside this repository.

## Contents

- `analysis/`: selected R/Python calculations in dependency order.
- `plotting/`: final 21 September figure implementation.
- `data/frozen/`: quick-route numerical inputs, annotations and image sources.
- `metadata/`: file checksums, run order, panel and worksheet recipes.
- `reference/`: independent manuscript baselines used only for validation.
- `docs/`: data acquisition, provenance and operating instructions.
- `reports/`: actual validation coverage and outstanding limitations.

Generated results include `figures/`, `tables/`, `source_data/`, `logs/`,
`run.json`, and `verification.json`. Numeric verification uses absolute
tolerance `1e-8` and relative tolerance `1e-6`; discrete identifiers, membership,
order and missingness must agree. Existing stricter scientific checks are retained.
PNG equality is environment-dependent and is reported separately from numbers.

Repository: [Elliott0114/crc-metastatic-growth-pattern-atlas](https://github.com/Elliott0114/crc-metastatic-growth-pattern-atlas).
See [Chinese instructions](docs/使用说明.md) and [release instructions](docs/GITHUB.md).
The GSE294385 source regional annotation is not redistributed; obtain it from
the authors for the full route. Quick reproduction does not require that file.

The [execution order](docs/RUN_ORDER.md), [panel/worksheet source map](metadata/reproduction_map.tsv),
and [measured validation report](reports/VALIDATION.md) document the selected workflow.
`CHECKSUMS.tsv` records every distributed file except the checksum list itself.

The inherited REC44 supplementary table uses a pre-correction spatial mask.
Its isolated compatibility branch is documented in [WORKFLOW.md](docs/WORKFLOW.md);
it does not replace the corrected main analysis.

## Citation And Terms

Author information follows the current draft; see [CITATION.cff](CITATION.cff).
The manuscript is not represented here as an accepted or published article.
Original code is MIT-licensed. Third-party data and annotations are **not**
relicensed under MIT; see [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
