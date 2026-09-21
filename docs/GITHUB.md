# Repository And Releases

The public repository is
[Elliott0114/crc-metastatic-growth-pattern-atlas](https://github.com/Elliott0114/crc-metastatic-growth-pattern-atlas).
The 21 September manuscript snapshot uses tag `manuscript-2026-09-21`.
Earlier commits and the `v1.0.0` release remain in the repository history.

```bash
git clone https://github.com/Elliott0114/crc-metastatic-growth-pattern-atlas.git
cd crc-metastatic-growth-pattern-atlas
git checkout manuscript-2026-09-21
```

Download `emtab-matrix-inputs-2026-09-21.tar.gz` and its `.sha256` file from
the [release page](https://github.com/Elliott0114/crc-metastatic-growth-pattern-atlas/releases/tag/manuscript-2026-09-21).
Check the archive before extracting it into the explicitly selected data root:

```bash
sha256sum -c emtab-matrix-inputs-2026-09-21.tar.gz.sha256
mkdir -p /absolute/path/to/inputs
tar -xzf emtab-matrix-inputs-2026-09-21.tar.gz -C /absolute/path/to/inputs
```

The GSE294385 source regional annotation is not uploaded. Obtain it from the
authors and follow [DATA.md](DATA.md). It is not needed for quick reproduction.

## Future Updates

Work in a separate checkout of this repository, never in the parent research
project. Review the staged files and third-party terms before committing.
Regenerate `CHECKSUMS.tsv` when distributed files change; do not alter input
hashes simply to accept unexpected data drift. Use normal commits and new
release tags, preserving published history. `.gitattributes` preserves recorded
file bytes across Git checkouts.

```bash
git add .
git diff --cached --stat
git diff --cached --check
git commit -m "Update manuscript reproduction release"
git push origin main
```

Keep raw matrices, FASTQs, caches, outputs, credentials, manuscripts and
submission correspondence outside Git. Companion matrices are Release assets,
not Git objects. MIT covers original code only; retain all source attribution
and the separate terms in `THIRD_PARTY_NOTICES.md`.
