"""Execute the count-matrix workflow with an explicit, file-level input contract."""
import json
import shutil
import sys
from pathlib import Path

from workflow_io import read_tsv, write_json


def run_full(root, workspace, data_dir, output, records, command_run):
    if any(workspace.iterdir()):
        raise RuntimeError("Full analysis requires an empty workspace")
    shutil.copytree(root / "analysis", workspace / "scripts")
    for row in read_tsv(root / "metadata/analysis_inputs.tsv"):
        relative = Path(row["workspace_path"])
        if relative.is_absolute() or ".." in relative.parts:
            raise RuntimeError(f"Unsafe staging path: {relative}")
        source = (root if row["location"] == "bundled" else data_dir) / row["path"]
        target = workspace / row["workspace_path"]
        target.parent.mkdir(parents=True, exist_ok=True)
        if row["location"] == "bundled":
            shutil.copy2(source, target)
        else:
            target.symlink_to(source.resolve())
    steps = json.loads((root / "metadata/analysis_steps.json").read_text(encoding="utf-8"))
    completed = []
    try:
        for step in steps:
            interpreter = [sys.executable, "-B"] if step["script"].endswith(".py") else ["Rscript", "--vanilla"]
            command_run(interpreter + ["scripts/" + step["script"]] + step.get("arguments", []), workspace, output / "logs" / (step["id"] + ".log"), records)
            completed.append(step["id"])
            write_json(output / "analysis_progress.json", dict(completed=completed, total=len(steps), current_status="running"))
        from materialize_sources import materialize
        materialize(root, workspace)
    except Exception:
        write_json(output / "analysis_progress.json", dict(completed=completed, total=len(steps), current_status="failed", frozen_results_used=False))
        raise
    write_json(output / "analysis_progress.json", dict(completed=completed, total=len(steps), current_status="completed", frozen_results_used=False))
    return len(completed)
