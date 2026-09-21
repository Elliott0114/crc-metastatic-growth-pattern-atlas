"""Record and check the packages actually used by the selected workflow."""
from importlib.metadata import PackageNotFoundError, version
import json
import platform
import shutil
import subprocess


def inspect_environment(root, mode, env):
    expected = json.loads((root / "metadata/environment_versions.json").read_text(encoding="utf-8"))
    errors, actual = [], {"python": platform.python_version()}
    wanted = expected["quick_python"] if mode == "quick" else expected["full_python"]
    if actual["python"] != expected["python"]:
        errors.append(f"Python version mismatch: {actual['python']}; expected {expected['python']}")
    for package, pinned in wanted.items():
        try:
            actual[package] = version(package)
        except PackageNotFoundError:
            errors.append(f"Missing Python package: {package}=={pinned}; see environment.yml")
            continue
        if actual[package] != pinned:
            errors.append(f"Version mismatch: {package}=={actual[package]}; expected {pinned}")
    if mode == "full":
        actual['system_commands'] = {}
        for command in ('gzip', 'zcat', 'awk', 'sha256sum'):
            available = shutil.which(command) is not None
            actual['system_commands'][command] = available
            if not available:
                errors.append(f"Missing system command: {command}; see docs/DATA.md")
        if shutil.which("Rscript") is None:
            errors.append("Rscript is missing; see environment.yml")
        else:
            packages = list(expected["r_packages"])
            code = 'cat("R\\t",as.character(getRversion()),"\\n",sep=""); for(p in c(' + ','.join(json.dumps(p) for p in packages) + ')) cat(p,"\\t",if(requireNamespace(p,quietly=TRUE)) as.character(packageVersion(p)) else "MISSING","\\n",sep="")'
            result = subprocess.run(["Rscript", "--vanilla", "-e", code], text=True, encoding="utf-8", capture_output=True, env=env)
            if result.returncode:
                errors.append("R environment check failed: " + result.stderr[-1500:])
            else:
                actual["R"] = dict(line.split("\t", 1) for line in result.stdout.splitlines() if "\t" in line)
                for package, pinned in {"R": expected["r_version"], **expected["r_packages"]}.items():
                    if actual["R"].get(package) != pinned:
                        errors.append(f"R version mismatch/missing: {package}={actual['R'].get(package)}; expected {pinned}")
    return actual, errors
