"""run_all.py - run every analysis script in order.   Usage: python run_all.py [--datasets mivia lv bi]
Figures and tables go to analysis/outputs/<validation|performance|ablation|bidomain>/ ."""
import os
import subprocess
import sys
from pathlib import Path

here = Path(__file__).resolve().parent
ds = sys.argv[1:] if sys.argv[1:2] == ["--datasets"] else []
env = dict(os.environ, MPLBACKEND="Agg")
for script, takes_ds in [("01_validation.py", True), ("02_performance.py", True),
                         ("03_ablation.py", False), ("05_reference_bugs.py", False)]: #, ("04_bidomain.py", False)]:
    print(f"\n{'=' * 78}\n{script}\n{'=' * 78}", flush=True)
    r = subprocess.run([sys.executable, "-u", str(here / script)] + (ds if takes_ds else []), env=env, cwd=here)
    if r.returncode:
        print(f"{script} failed (exit {r.returncode})"); sys.exit(r.returncode)
print("\nAll done. Outputs are in", here / "outputs")
