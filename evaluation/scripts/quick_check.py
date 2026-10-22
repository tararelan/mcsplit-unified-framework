# check_fidelity.py - unified vs reference node counts, before vs after the fix (one dataset)
# Usage: python check_fidelity.py [algo] [dataset]     e.g. python check_fidelity.py mcsplit bi
import sys
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

ALGO = sys.argv[1] if len(sys.argv) > 1 else "mcsplit"
DS = sys.argv[2] if len(sys.argv) > 2 else "bi"
RUNS = {"before fix (old)": "hpc/publication/old", "after fix (results)": "hpc/publication/results"}
TIMEOUT_S = 999

# find the repo root (folder containing hpc/publication)
root = next(p for p in [Path.cwd(), *Path.cwd().parents] if (p / "hpc/publication").exists())

fig, axes = plt.subplots(1, len(RUNS), figsize=(5 * len(RUNS), 4.6), sharex=True, sharey=True)

for ax, (run_name, folder) in zip(axes, RUNS.items()):
    path = root / folder / f"{ALGO}_{DS}_with_ref_all.csv"
    if not path.exists():
        ax.set_title(f"{run_name}\n(missing: {path.name})", fontsize=9)
        print(f"{run_name:20} MISSING {path}")
        continue

    d = pd.read_csv(path)
    for c in ["unified_nodes", "ref_nodes", "unified_time", "ref_time"]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    # completed pairs on both sides, sizes agree
    d = d[(d["unified_time"] <= TIMEOUT_S) & (d["ref_time"] <= TIMEOUT_S) & (d["match"] == "PASS")]
    d = d[(d["unified_nodes"] > 0) & (d["ref_nodes"] > 0)]

    exact = d["unified_nodes"] == d["ref_nodes"]
    ratio = d["ref_nodes"] / d["unified_nodes"]
    print(f"{run_name:20} pairs={len(d):6}  exact={exact.mean():6.1%}  "
          f"ref>uni={(ratio > 1).mean():5.1%}  ref<uni={(ratio < 1).mean():5.1%}  "
          f"median ratio={ratio.median():.3f}")

    ax.scatter(d.loc[exact, "unified_nodes"], d.loc[exact, "ref_nodes"],
               s=7, alpha=0.4, color="#0072B2", edgecolors="none", label=f"exact ({exact.sum()})")
    ax.scatter(d.loc[~exact, "unified_nodes"], d.loc[~exact, "ref_nodes"],
               s=7, alpha=0.4, color="#D55E00", edgecolors="none", label=f"differs ({(~exact).sum()})")
    lo = min(d["unified_nodes"].min(), d["ref_nodes"].min())
    hi = max(d["unified_nodes"].max(), d["ref_nodes"].max())
    ax.plot([lo, hi], [lo, hi], "k--", lw=0.8)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("unified nodes")
    ax.set_title(f"{run_name}  (exact {exact.mean():.1%})", fontsize=10)
    ax.legend(fontsize=8, loc="upper left")

axes[0].set_ylabel("reference nodes")
fig.suptitle(f"{ALGO} on {DS}: unified vs reference search nodes (completed, size-matching pairs)")
fig.tight_layout()
out = root / f"fidelity_{ALGO}_{DS}.png"
fig.savefig(out, dpi=150)
print("saved", out)
plt.show()