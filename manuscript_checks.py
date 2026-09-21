"""Numeric statements in the manuscript, independent of document formatting."""
import numpy as np
import pandas as pd
import json
from pathlib import Path


def check_claims(output):
    def read(name):
        return pd.read_csv(output / "source_data/panels" / name, sep="\t")

    a, b = read("Figure_2_A.tsv"), read("Figure_2_B.tsv")
    rows = []

    def check(label, actual, expected):
        rows.append(dict(claim=label, actual=actual, expected=expected, passed=bool(actual == expected)))

    check("Discovery patients (Fig. 2A)", int(a.patient.nunique()), 14)
    check("Independent validation patients (Fig. 2B)", int(b.patient.nunique()), 28)
    check("Independent validation studies", int(b.study.nunique()), 5)
    check("Independent malignant cells", int(b.n_cells.sum()), 59095)
    check("Positive patient component differences", int((b.delta_z > 0).sum()), 23)
    check("Study patient counts, sorted", b.groupby("study").size().sort_values().tolist(), [3, 4, 4, 5, 12])
    d = read("Figure_3_D.tsv")
    check("Predefined junction members", int(d.gene.nunique()), 30)
    check("Study-by-member estimates", int(d.estimate.notna().sum()), 125)
    check("Unavailable study-by-member entries", int(d.estimate.isna().sum()), 25)
    check("Fisher-z display inversion", bool(np.allclose(d.rho, np.tanh(d.estimate), equal_nan=True, atol=1e-8, rtol=1e-6)), True)
    h = read("Figure_5_C.tsv")
    check("HRC coefficient pairs", int(h.gene.nunique()), 88)
    for key, expected in [("unadjusted", .559432621), ("source_adjusted", .137163329)]:
        actual = float(h[key].median())
        rows.append(dict(claim=f"Median HRC {key} coefficient", actual=actual, expected=expected, passed=bool(np.isclose(actual, expected, atol=1e-8, rtol=1e-6))))
    regional = read("Supplementary_Figure_6_A.tsv")
    counts = regional.groupby(["dataset", "normalization", "model"]).count_eligible.sum()
    check("Regional specifications", len(counts), 8)
    check("Eligible junction genes in every specification", list(map(int, counts)), [24] * 8)
    spatial = read("Supplementary_Figure_4_region_counts.tsv")
    check("HGP spatial patients", int(spatial.patient.nunique()), 6)
    check("Mutually exclusive displayed spatial spots", int(spatial.spots.sum()), 21629)
    registry = Path(__file__).resolve().parent / 'metadata/manuscript_numeric_claims.json'
    if registry.exists():
        for rule in json.loads(registry.read_text(encoding='utf-8')):
            kind, source = rule['source'].split(':', 1)
            path = (output / 'source_data/panels' if kind == 'panel' else output / 'workspace') / source
            frame = pd.read_csv(path, sep='\t')
            for column, value in rule['filters'].items():
                frame = frame[frame[column] == value]
            values = frame[rule['column']]
            if rule['method'] == 'one':
                if len(values) != 1:
                    raise ValueError(f"Non-unique numeric claim source: {rule['claim']}")
                actual = float(values.iloc[0])
            else:
                actual = float(getattr(values, rule['method'])())
            if rule.get('transform') == 'tanh':
                actual = float(np.tanh(actual))
            actual *= rule['scale']
            expected = float(rule['reported'])
            rows.append(dict(claim=rule['claim'], section=rule['section'], source=rule['source'],
                             actual=actual, expected=expected, rounding_atol=rule['rounding_atol'],
                             passed=abs(actual-expected) <= rule['rounding_atol']))
    return rows
