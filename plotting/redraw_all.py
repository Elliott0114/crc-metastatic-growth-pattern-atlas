"""Current figure entry point; use frozen inputs and the final Figure 1 framework."""
from pathlib import Path
import os
import subprocess
import sys

REVISION = Path(__file__).resolve().parent
env = {k:v for k,v in os.environ.items() if k.lower() not in ['http_proxy','https_proxy','all_proxy']}
for script, args in [('build_main_figures.py', ['2','3','5','4']),
                     ('build_supplementary_figures.py', []),
                     ('build_figure1_framework.py', [])]:
    subprocess.run([sys.executable, str(REVISION/script), *args], env=env, check=True)
