"""Reader summaries: fixed wording, values read from full-precision source tables."""
from functools import lru_cache
import json
import math

import pandas as pd


@lru_cache(maxsize=128)
def frame(path):
    return pd.read_csv(path, sep="\t")


def value(binding, workspace, identity):
    data = frame(workspace / binding["source"])
    if binding.get("operation") == "fraction_positive":
        for column, target in binding["group"].items():
            data = data[data[column] == target]
        if data.empty:
            raise ValueError(f"Empty source selection: {binding}")
        return float((data[binding["column"]] > 0).mean())
    row = data.iloc[binding["row"]]
    for column, expected in identity.items():
        if column in row and str(row[column]) != expected:
            raise ValueError(f"Summary source identity changed: {binding['source']}, {column}")
    return float(row[binding["column"]]) * binding.get("scale", 1)


def number(value):
    return format(value, ".4g") if math.isfinite(value) else "\u2014"


def rebuild_summary(table, workspace, root):
    if table == 1:
        recipe = json.loads((root / "metadata/programme_summary.json").read_text(encoding="utf-8"))
        coverage = frame(workspace / recipe["sources"][0])
        definitions = frame(workspace / recipe["sources"][1])
        lock = frame(workspace / recipe["sources"][2])
        contexts = frame(workspace / recipe["sources"][3])
        programs = ["FROZEN_REC_TOP50", "CANELLAS_CORE_HRC", "REACTOME_TIGHT_JUNCTION_INTERACTIONS", "HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT", "WHITE_RSC", "WHITE_CBC", "ICMS2_TEMPLATE_UP", "ICMS3_TEMPLATE_UP"]
        rows = [recipe["headers"]]
        for labels, program in zip(recipe["labels"][:9], programs, strict=True):
            ogden = coverage[(coverage.dataset == "Ogden_patient_state") & (coverage.set_id == program)].iloc[0]
            bulk = coverage[(coverage.dataset == "GSE151165") & (coverage.set_id == program)].iloc[0]
            measured = int(bulk.n_measured)
            if program == "REACTOME_TIGHT_JUNCTION_INTERACTIONS":
                measured = int(contexts[(contexts.dataset == "GSE151165") & (contexts.component == "Junction")].measured_variable.sum())
            rows.append(labels + [str(int(bulk.n_defined)), str(int(ogden.n_measured)), str(measured)])
        for labels, component in zip(recipe["labels"][9:13], ["Junction", "Claudin", "Polarity", "HRC"], strict=True):
            n_defined = definitions[definitions.component == component].gene.nunique()
            n_ogden = lock[(lock.component == component) & (lock.role == "target")].measured.sum()
            n_bulk = contexts[(contexts.dataset == "GSE151165") & (contexts.component == component)].measured_variable.sum()
            rows.append(labels + [str(n_defined), str(int(n_ogden)), str(int(n_bulk))])
        spatial_hrc = set(definitions.loc[definitions.component == "HRC", "gene"]) - {"CEACAM5"}
        rows.append(recipe["labels"][13] + [str(len(spatial_hrc)), "\u2014", "\u2014"])
        return rows
    book = json.loads((root / "metadata/summary_recipes.json").read_text(encoding="utf-8"))[str(table)]
    rows = [book["headers"]]
    for recipe in book["rows"]:
        cells = list(recipe["labels"])
        for template, bindings in zip(recipe["templates"], recipe["bindings"], strict=True):
            values = [number(value(binding, workspace, recipe["identity"])) for binding in bindings]
            cells.append(template.format(*values))
        rows.append(cells)
    return rows
