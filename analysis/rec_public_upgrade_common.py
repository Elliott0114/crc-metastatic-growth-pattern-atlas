"""Shared I/O for the public REC extension; no outcome-dependent signature fitting."""
from __future__ import annotations

import csv
import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
VERSION = "rec_public_upgrade_2026-09-07"
META = ROOT / "metadata" / VERSION
OUT = ROOT / "analysis_results" / VERSION
MS = ROOT / "manuscript" / VERSION
DATA = ROOT / "data_sources" / VERSION
OLD = ROOT / "analysis_results/deep_biology_upgrade_2026-08-31/phase2_mechanistic_specificity"


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def soft_header(lines):
    header = {}
    for line in lines:
        if line.startswith("!series_matrix_table_begin"):
            break
        if line.startswith("!"):
            row = next(csv.reader([line], delimiter="\t"))
            header.setdefault(row[0], []).append(row[1:])
    return header
