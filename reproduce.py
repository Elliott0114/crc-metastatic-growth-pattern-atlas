#!/usr/bin/env python3
"""Reproduce the 21 September 2026 manuscript in an isolated workspace."""
from pathlib import Path
import argparse
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
import time

from workflow_io import read_tsv, read_workbook, sha256, write_json, write_workbook

ROOT = Path(__file__).resolve().parent
PAPER = Path("manuscript/rec_text_polish_2026-09-21")


def clean_env():
    env = {k: v for k, v in os.environ.items() if k.lower() not in ("http_proxy", "https_proxy", "all_proxy")}
    env.update(OPENBLAS_NUM_THREADS="2", OMP_NUM_THREADS="2", MKL_NUM_THREADS="2", PYTHONDONTWRITEBYTECODE="1", MPLBACKEND="Agg")
    return env


def validate_files(rows, base, hint):
    errors = []
    for row in rows:
        relative = Path(row["path"])
        if relative.is_absolute() or ".." in relative.parts:
            errors.append(f"Unsafe manifest path: {row['path']}")
            continue
        path = base / relative
        if not path.is_file():
            errors.append(f"Missing input: {row['path']}; source: {row.get('source', hint)}")
        elif path.stat().st_size != int(row["bytes"]) or sha256(path) != row["sha256"]:
            errors.append(f"Checksum mismatch: {row['path']}; restore the recorded version ({row.get('source', hint)})")
    return errors


def check(mode, data_dir=None):
    errors = []
    from environment_check import inspect_environment
    versions, environment_errors = inspect_environment(ROOT, mode, clean_env())
    errors.extend(environment_errors)
    rows = read_tsv(ROOT / "metadata/frozen_files.tsv")
    errors.extend(validate_files(rows, ROOT, "the same repository release; see README.md"))
    integrity = ROOT / "CHECKSUMS.tsv"
    if integrity.exists():
        errors.extend(validate_files(read_tsv(integrity), ROOT, "the same repository release"))
    if mode == "full":
        if data_dir is None:
            errors.append("Full analysis requires --data-dir; see docs/DATA.md")
        elif not (ROOT / "metadata/analysis_inputs.tsv").is_file():
            errors.append("Full-analysis input manifest is missing")
        else:
            for row in read_tsv(ROOT / "metadata/analysis_inputs.tsv"):
                errors.extend(validate_files([row], ROOT if row["location"] == "bundled" else data_dir, "docs/DATA.md"))
    if errors:
        raise RuntimeError("\n".join(errors))
    return {"mode": mode, "bundled_files_checked": len(rows), "environment": versions, "passed": True}


def command_run(command, workspace, log, records):
    started = time.monotonic()
    print("Running", " ".join(command[1:]), flush=True)
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("w", encoding="utf-8") as stream:
        result = subprocess.run(command, cwd=workspace, env=clean_env(), stdout=stream, stderr=subprocess.STDOUT)
    records.append({"command": command[1:], "returncode": result.returncode, "seconds": round(time.monotonic() - started, 2), "log": log.name})
    if result.returncode:
        raise RuntimeError(f"Command failed: {command[1:]}; see {log}")


def stage_plotting(workspace, mode="quick"):
    target = workspace / PAPER
    shutil.copytree(ROOT / "plotting", target / "revision", dirs_exist_ok=True)
    shutil.copytree(ROOT / "assets", target / "assets", dirs_exist_ok=True)
    shutil.copy2(ROOT / "assets/Study_design_sources.tsv", target / "Study_design_sources.tsv")
    shutil.copy2(ROOT / "reference/baseline_panel_hashes.json", target / "baseline_panel_hashes.json")
    if mode == "full":
        shutil.copytree(ROOT / "reference/panels", target / "reference_panels")
    write_json(target / "validation_mode.json", dict(mode=mode))
    write_json(target / "panel_manifest.json", [])
    return target


def finalize_panel_index(target):
    import pandas as pd
    records = json.loads((target / 'panel_manifest.json').read_text(encoding='utf-8'))
    for row in records:
        row['panel_data_sha256'] = sha256(target / row['panel_source_data'])
        row['output_sha256'] = sha256(target / row['output_path'])
    write_json(target / 'panel_manifest.json', records)
    pd.DataFrame(records).to_csv(target / 'panel_manifest.tsv', sep='\t', index=False, encoding='utf-8')


def build_workbooks(output, workspace, mode):
    recipes = json.loads((ROOT / "metadata/workbook_recipes.json").read_text(encoding="utf-8"))
    books = {}
    for recipe in recipes:
        from workbook_sources import rebuilt_sheet
        rows = rebuilt_sheet(recipe, workspace, ROOT)
        books.setdefault(recipe["table"], []).append((recipe["sheet"], rows))
    for number, sheets in books.items():
        write_workbook(output / "tables" / f"Supplementary_Table_{number}.xlsx", sheets)


def verify(output):
    import numpy as np
    import pandas as pd
    state = json.loads((output / "run.json").read_text(encoding="utf-8"))
    checks = []
    from PIL import Image
    pixels = {r["figure"]: r for r in json.loads((ROOT / "reference/figure_pixels.json").read_text(encoding="utf-8"))}
    for row in json.loads((ROOT / "reference/figure_hashes.json").read_text(encoding="utf-8")):
        path = output / "figures" / f"{row['figure']}.png"
        pixel_equal = False
        if path.is_file():
            with Image.open(path) as im:
                rgba = im.convert("RGBA")
                ref = pixels[row["figure"]]
                pixel_equal = rgba.size == (ref["width"], ref["height"]) and hashlib.sha256(rgba.tobytes()).hexdigest() == ref["rgba_sha256"]
        checks.append(dict(item=str(path.relative_to(output)), passed=pixel_equal, png_bytes_equal=path.is_file() and sha256(path) == row["sha256"], method="RGBA pixels and dimensions; PNG bytes recorded separately"))
    for ref in sorted((ROOT / "reference/panels").glob("*.tsv")):
        path = output / "source_data/panels" / ref.name
        detail = ""
        try:
            a, b = pd.read_csv(ref, sep="\t"), pd.read_csv(path, sep="\t")
            pd.testing.assert_frame_equal(a, b, check_dtype=False, check_exact=False, atol=1e-8, rtol=1e-6)
            for column in a.columns:
                if pd.api.types.is_integer_dtype(a[column].dtype):
                    pd.testing.assert_series_equal(a[column], b[column], check_dtype=False, check_exact=True)
            good = True
        except (AssertionError, FileNotFoundError) as exc:
            good = False
            detail = str(exc)[:1800]
        checks.append(dict(item="panels/" + ref.name, passed=good, method="columns, ordering, missingness and values", detail=detail))
    recipes = json.loads((ROOT / "metadata/workbook_recipes.json").read_text(encoding="utf-8"))
    books = {i: dict(read_workbook(output / "tables" / f"Supplementary_Table_{i}.xlsx")) for i in range(1, 9)}
    for i, book in books.items():
        expected = [r["sheet"] for r in recipes if r["table"] == i]
        checks.append(dict(item=f"Table {i}/worksheet order", passed=list(book) == expected, method="exact worksheet names and order"))
    for row in recipes:
        ref = json.loads((ROOT / row["reference"]).read_text(encoding="utf-8"))
        actual = books[row["table"]].get(row["sheet"])
        same = cells_equal(actual, ref)
        if row["kind"] == "provenance":
            same = verify_provenance(row, actual, output / "workspace")
        checks.append(dict(item=f"Table {row['table']}/{row['sheet']}", passed=same, method="typed cell content" if row["kind"] != "provenance" else "source index rebuilt for this run"))
    from manuscript_checks import check_claims
    claims = check_claims(output)
    completed = state.get("analysis_steps_completed", 0)
    expected_steps = len(json.loads((ROOT / "metadata/analysis_steps.json").read_text(encoding="utf-8")))
    if state["mode"] == "full":
        checks.append(dict(item="analysis completion", passed=completed == expected_steps, method="all declared modules completed"))
    report = dict(mode=state["mode"], analysis_steps_completed=completed, models_refitted=state["mode"] == "full" and completed == expected_steps, figures=11, workbooks=8, worksheets=127,
                  passed=all(x["passed"] for x in checks + claims), checks=checks, manuscript_claims=claims,
                  tolerance={"atol": 1e-8, "rtol": 1e-6})
    write_json(output / "verification.json", report)
    if not report["passed"]:
        bad = [x.get("item", x.get("claim")) for x in checks + claims if not x["passed"]]
        raise RuntimeError(f"Verification failed: {bad}; see {output / 'verification.json'}")
    print(f"Verified {len(checks)} figure, panel and workbook checks", flush=True)
    state['status'] = 'verified'
    write_json(output / 'run.json', state)
    return report


def cells_equal(actual, expected):
    if actual is None or len(actual) != len(expected):
        return False
    for a, b in zip(actual, expected, strict=True):
        if len(a) != len(b):
            return False
        for x, y in zip(a, b, strict=True):
            if type(x) is bool or type(y) is bool:
                if type(x) is not type(y) or x != y:
                    return False
            elif type(x) is int and type(y) is int:
                if x != y:
                    return False
            elif isinstance(x, (int, float)) and isinstance(y, (int, float)):
                if not math.isclose(x, y, abs_tol=1e-8, rel_tol=1e-6):
                    return False
            elif x != y:
                return False
    return True


def verify_provenance(recipe, rows, workspace):
    if not rows:
        return False
    if recipe["sheet"] == "Added_source_index":
        if rows[0] != ["sheet", "source", "sha256"] or len(rows) < 2:
            return False
        recipes = json.loads((ROOT / "metadata/workbook_recipes.json").read_text(encoding="utf-8"))
        expected = []
        for other in recipes:
            if other['table'] != recipe['table'] or other['kind'] == 'provenance':
                continue
            for source in other.get('sources', [other['source']] if other.get('source') else []):
                expected.append([other['sheet'], source])
        if [row[:2] for row in rows[1:]] != expected:
            return False
        for _, source, digest in rows[1:]:
            path = workspace / source
            if not path.is_file() or sha256(path) != digest:
                return False
        return True
    if recipe["sheet"] == "Figure_panel_index":
        if len(rows) != 53 or "panel_id" not in rows[0]:
            return False
        records = [dict(zip(rows[0], row, strict=True)) for row in rows[1:]]
        expected = read_tsv(ROOT / "metadata/panel_manifest.tsv")
        if {(r['figure_id'], r['panel_id']) for r in records} != {(r['figure_id'], r['panel_id']) for r in expected}:
            return False
        for figure_id in {r['figure_id'] for r in expected}:
            if ([r['panel_id'] for r in records if r['figure_id'] == figure_id]
                    != [r['panel_id'] for r in expected if r['figure_id'] == figure_id]):
                return False
        for record in records:
            hashes = json.loads(record['source_hashes'])
            if set(hashes) != set(json.loads(record['source_result_paths'])):
                return False
            for source, digest in hashes.items():
                if not (workspace / source).is_file() or sha256(workspace / source) != digest:
                    return False
            panel = workspace / PAPER / record['panel_source_data']
            if not panel.is_file() or sha256(panel) != record['panel_data_sha256']:
                return False
            figure = workspace / PAPER / record['output_path']
            if not figure.is_file() or sha256(figure) != record['output_sha256']:
                return False
        return True
    return False


def reproduce(args):
    preflight = check(args.command, args.data_dir)
    output = args.output_dir
    if output.exists() and any(output.iterdir()):
        raise RuntimeError(f"Use a new --output-dir; an existing run is preserved: {output}")
    workspace = output / "workspace"
    workspace.mkdir(parents=True, exist_ok=True)
    records = []
    state = dict(mode=args.command, status="running", analysis_snapshot="crc-junction-q2-analysis-2026-09-19", figure_snapshot="crc-junction-text-polish-2026-09-21", analysis_steps_completed=0, jobs=records, preflight=preflight)
    write_json(output / "run.json", state)
    try:
        if args.command == "quick":
            shutil.copytree(ROOT / "data/frozen", workspace, dirs_exist_ok=True)
        else:
            from run_analysis import run_full
            state["analysis_steps_completed"] = run_full(ROOT, workspace, args.data_dir, output, records, command_run)
        target = stage_plotting(workspace, args.command)
        command_run([sys.executable, str(PAPER / "revision/redraw_all.py")], workspace, output / "logs/figures.log", records)
        finalize_panel_index(target)
        for folder, extension in (("figures", "pdf"), ("editable_figures", "svg"), ("review", "png")):
            (output / "figures").mkdir(exist_ok=True)
            for source in (target / folder).glob(f"*.{extension}"):
                shutil.copy2(source, output / "figures" / source.name)
        shutil.copytree(target / "source_data", output / "source_data", dirs_exist_ok=True)
        build_workbooks(output, workspace, args.command)
        state["status"] = "completed_unverified"
    except Exception as exc:
        state.update(status="failed", error=str(exc))
        raise
    finally:
        progress = output / "analysis_progress.json"
        if progress.exists():
            state["analysis_steps_completed"] = len(json.loads(progress.read_text(encoding="utf-8"))["completed"])
        write_json(output / "run.json", state)
    verify(output)
    state["status"] = "verified"
    write_json(output / "run.json", state)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("check", "quick", "full", "verify"))
    parser.add_argument("--mode", choices=("quick", "full"), default="quick", help="Mode to preflight with check")
    parser.add_argument("--data-dir", type=lambda s: Path(s).expanduser().resolve())
    parser.add_argument("--output-dir", type=lambda s: Path(s).expanduser().resolve(), default=ROOT / "outputs")
    args = parser.parse_args()
    try:
        if args.command == "check":
            print(json.dumps(check(args.mode, args.data_dir), indent=2))
        elif args.command == "verify":
            verify(args.output_dir)
        else:
            reproduce(args)
    except (RuntimeError, FileNotFoundError) as exc:
        parser.exit(1, str(exc) + "\n")


if __name__ == "__main__":
    main()
