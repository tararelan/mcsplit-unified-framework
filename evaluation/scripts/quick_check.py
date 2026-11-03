# quick_check.py - unified vs reference for ALL algorithms (one dataset), two graphs:
#   1) search nodes   2) runtime
# Loops over the ALGOS list below. DSB and SymSplit result CSVs written by the older run_instance.sh have
# rows split in two by a stray "BIDOMAINS:" line; load() stitches those back together in memory and
# recomputes `match` (the CSV files themselves are never modified). Clean CSVs pass through unchanged.
# Usage: python quick_check.py [dataset] [results_folder]
#   e.g. python quick_check.py bi
#        python quick_check.py bi hpc/publication/old      (the before-fix runs)
import io
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

DS = sys.argv[1] if len(sys.argv) > 1 else "bi"
FOLDER = sys.argv[2] if len(sys.argv) > 2 else "hpc/publication/results"
# ALGOS = ["mcsplit", "down", "rl", "rl_pair", "rl_rate", "rl_down",
#          "dal", "ll", "dsb", "rrsplit", "symsplit"]
ALGOS = ["mcsplit", "down", "rl", "ll", "dsb"]
TIMEOUT_S = 999
# algorithms with no reference binary of their own are compared against another algorithm's reference
# (reference columns are taken from that algorithm's CSV, matched on instance pair)
REF_FROM = {"rl_pair": "rl", "rl_rate": "rl", "rl_down": "rl"}
MIN_TIME = 1e-3   # reference clocks are ms-resolution; floor so the log axes work

# find the repo root (folder containing hpc/publication)
root = next(p for p in [Path.cwd(), *Path.cwd().parents] if (p / "hpc/publication").exists())

NCOLS = 4
NROWS = -(-len(ALGOS) // NCOLS)


def new_grid():
    fig, axes = plt.subplots(NROWS, NCOLS, figsize=(4.2 * NCOLS, 3.8 * NROWS))
    return fig, axes.ravel()


def read_repaired(path):
    """Read a result CSV, joining rows split by a stray 'BIDOMAINS:' line (DSB / SymSplit)."""
    fixed, joined = [], 0
    for line in path.read_text().splitlines():
        if line.startswith("BIDOMAINS:") and fixed:
            # previous line ended with the first part of ref_size; this line carries the rest of that
            # field ("BIDOMAINS:<mean>") then ref_nodes, ref_time, bd mean, bd max, match
            fixed[-1] += "," + ",".join(line.split(",")[1:])
            joined += 1
        else:
            fixed.append(line)
    d = pd.read_csv(io.StringIO("\n".join(fixed)))
    if joined:
        # stored `match` compared a two-line string, so recompute it from the sizes
        us = pd.to_numeric(d["unified_size"], errors="coerce")
        rs = pd.to_numeric(d["ref_size"], errors="coerce")
        d["match"] = "MISSING_DATA"
        d.loc[us.notna() & rs.notna() & (us == rs), "match"] = "PASS"
        d.loc[us.notna() & rs.notna() & (us != rs), "match"] = "FAIL"
        print(f"{'':10} (stitched {joined} split rows)")
    return d


def load(algo):
    """Completed, size-matching pairs for one algorithm, or None if there is nothing to compare."""
    path = root / FOLDER / f"{algo}_{DS}_with_ref_all.csv"
    if not path.exists():
        print(f"{algo:10} MISSING {path.name}")
        return None
    d = read_repaired(path)
    ref_name = REF_FROM.get(algo)
    if ref_name:
        rpath = root / FOLDER / f"{ref_name}_{DS}_with_ref_all.csv"
        if not rpath.exists():
            print(f"{algo:10} MISSING {rpath.name} (needed for its reference)")
            return None
        r = read_repaired(rpath)[["instance_a", "instance_b", "ref_size", "ref_nodes", "ref_time"]]
        d = d.drop(columns=["ref_size", "ref_nodes", "ref_time"]).merge(r, on=["instance_a", "instance_b"], how="inner")
        us = pd.to_numeric(d["unified_size"], errors="coerce")
        rs = pd.to_numeric(d["ref_size"], errors="coerce")
        d["match"] = "MISSING_DATA"
        d.loc[us.notna() & rs.notna() & (us == rs), "match"] = "PASS"
        d.loc[us.notna() & rs.notna() & (us != rs), "match"] = "FAIL"
        print(f"{algo:10} compared against the {ref_name} reference ({len(d)} pairs joined)")
    for c in ["unified_nodes", "ref_nodes", "unified_time", "ref_time"]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d = d[(d["unified_time"] <= TIMEOUT_S) & (d["ref_time"] <= TIMEOUT_S) & (d["match"] == "PASS")]
    d = d[(d["unified_nodes"] > 0) & (d["ref_nodes"] > 0)]
    if d.empty:
        print(f"{algo:10} no reference data (no reference binary configured, or nothing completed)")
        return None
    return d


def empty_panel(ax, algo, why):
    ax.set_title(algo, fontsize=10)
    ax.text(0.5, 0.5, why, ha="center", va="center", fontsize=9, color="grey", transform=ax.transAxes)
    ax.set_xticks([]); ax.set_yticks([])


fig_n, axes_n = new_grid()   # search nodes
fig_t, axes_t = new_grid()   # runtime

for i, algo in enumerate(ALGOS):
    d = load(algo)
    if d is None:
        empty_panel(axes_n[i], algo, "no reference data")
        empty_panel(axes_t[i], algo, "no reference data")
        continue

    # ---- graph 1: search nodes ----
    ax = axes_n[i]
    exact = d["unified_nodes"] == d["ref_nodes"]
    ratio = d["ref_nodes"] / d["unified_nodes"]
    print(f"{algo:10} pairs={len(d):6}  nodes: exact={exact.mean():6.1%}  "
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
    ax.set_xlabel("unified nodes"); ax.set_ylabel(f"{REF_FROM.get(algo, '')} reference nodes".strip())
    ax.set_title(f"{algo}  (exact {exact.mean():.1%}, n={len(d)})", fontsize=10)
    ax.legend(fontsize=8, loc="upper left")

    # ---- graph 2: runtime ----
    ax = axes_t[i]
    ut = d["unified_time"].clip(lower=MIN_TIME)
    rt = d["ref_time"].clip(lower=MIN_TIME)
    tratio = rt / ut
    print(f"{'':10} {'':12} time:  median ref/unified={tratio.median():.2f}  "
          f"unified faster={(tratio > 1).mean():5.1%}")
    ax.scatter(ut, rt, s=7, alpha=0.4, color="#009E73", edgecolors="none")
    lo, hi = min(ut.min(), rt.min()), max(ut.max(), rt.max())
    ax.plot([lo, hi], [lo, hi], "k--", lw=0.8)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("unified time (s)"); ax.set_ylabel(f"{REF_FROM.get(algo, '')} reference time (s)".strip())
    ax.set_title(f"{algo}  (median ref/unified {tratio.median():.2f}x, n={len(d)})", fontsize=10)

# hide any unused panels
for k in range(len(ALGOS), len(axes_n)):
    axes_n[k].axis("off"); axes_t[k].axis("off")

fig_n.suptitle(f"{DS}: unified vs reference search nodes (completed, size-matching pairs)")
fig_t.suptitle(f"{DS}: unified vs reference runtime (completed, size-matching pairs)")
fig_n.tight_layout(); fig_t.tight_layout()
out_n = root / f"fidelity_nodes_{DS}.png"
out_t = root / f"fidelity_runtime_{DS}.png"
# fig_n.savefig(out_n, dpi=150)
# fig_t.savefig(out_t, dpi=150)
print("would save", out_n, "and", out_t)
plt.show()
