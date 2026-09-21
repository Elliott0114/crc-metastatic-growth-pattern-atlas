"""Assemble supplementary workbooks from declared tabular sources."""
from pathlib import Path
import json
import math

import pandas as pd


def source_frame(recipe, workspace):
    path = workspace / recipe["source"]
    frame = pd.read_csv(path, sep="\t")
    operation = recipe.get("operation", "select")
    if operation == "ma_bins":
        frame = pd.read_csv(path, sep="\t")
        curves = []
        for (dataset, normalization), group in frame[frame.count_eligible].groupby(["dataset", "normalization"]):
            group = group.copy()
            group["bin"] = pd.qcut(group.A, 20, duplicates="drop")
            curve = group.groupby("bin", observed=True)[["A", "M"]].median().reset_index()
            curve["bin"] = curve["bin"].astype(str)
            curve["dataset"] = dataset
            curve["normalization"] = normalization
            curves.append(curve)
        frame = pd.concat(curves, ignore_index=True)
    elif operation == "bulk_display":
        frame = pd.read_csv(path, sep="\t")
        definitions = pd.read_csv(workspace / "analysis_results/rec_junction_context_2026-09-16/definitions.tsv", sep="\t")
        genes = set(definitions.loc[definitions.component.isin(["Claudin", "Polarity", "F11R", "Junction"]), "gene"])
        column = frame.columns[0]
        numeric = frame[frame[column].isin(genes)].set_index(column).select_dtypes(include="number")
        sd = numeric.std(axis=1, ddof=1)
        numeric = numeric.loc[sd > 0]
        frame = numeric.sub(numeric.mean(axis=1), axis=0).div(sd[sd > 0], axis=0).reset_index()
    return frame


def typed_rows(frame, recipe):
    columns = recipe["columns"]
    frame = frame.iloc[recipe["rows"]][columns]
    rows = [columns]
    for record in frame.itertuples(index=False, name=None):
        converted = []
        for j, (value, kind) in enumerate(zip(record, recipe["types"], strict=True)):
            if pd.isna(value) or value == "":
                converted.append(recipe["missing_values"][j])
            elif kind == "number":
                value = float(value)
                converted.append(int(value) if value.is_integer() else value)
            elif kind == "bool":
                converted.append(str(value).lower() == "true")
            else:
                converted.append(str(value))
        rows.append(converted)
    return rows


def rebuilt_sheet(recipe, workspace, root):
    if recipe["kind"] == "summary":
        from summary_sources import rebuild_summary
        return rebuild_summary(recipe["table"], workspace, root)
    if recipe["kind"] == "provenance":
        if recipe["sheet"] == "Figure_panel_index":
            frame = pd.read_csv(workspace / "manuscript/rec_text_polish_2026-09-21/panel_manifest.tsv", sep="\t", keep_default_na=False)
            return [list(frame.columns)] + frame.values.tolist()
        from workflow_io import sha256
        all_recipes = json.loads((root / "metadata/workbook_recipes.json").read_text(encoding="utf-8"))
        rows = [["sheet", "source", "sha256"]]
        for other in all_recipes:
            if other["table"] != recipe["table"] or other["kind"] == "provenance":
                continue
            sources = other.get("sources", [other["source"]] if other.get("source") else [])
            for source in sources:
                path = workspace / source
                rows.append([other["sheet"], source, sha256(path)])
        return rows
    return typed_rows(source_frame(recipe, workspace), recipe)
