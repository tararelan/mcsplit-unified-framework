"""
07_pairwise.py - pairwise node-count and runtime scatter grids      (supplementary information)

For every dataset, on the shared medium set, one grid per metric (nodes, time). Each panel compares two
algorithms instance by instance: x = the later algorithm in the list, y = the earlier one, one dot per
pair of graphs, log-log, with the dashed line y = x. A dot below the line means the x-algorithm needed
MORE than the y-algorithm on that instance; a dot above it means LESS. Runs that hit the timeout sit on the
cap for time (and at whatever node count they reached for nodes).

Only the upper triangle is drawn (A vs B and B vs A carry the same information).

Usage:  python 07_pairwise.py [--datasets mivia lv bi] [--algos mcsplit rl ll dsb ...] [--include-unsound]
        (default: every algorithm with results; pass --algos to draw a smaller grid)
Output: analysis/outputs/performance/pairwise_{nodes,time}_<ds>.png
        --source reference draws the same grids from the reference binaries (pairwise_*_<ds>_ref.png)
"""
import argparse

import numpy as np

import config as C

FOLDER = "performance"
METRICS = [("nodes", "Nodes", "nodes"), ("time", "Time (s)", "time")]
DPI = 150          # these grids are large (n-1 x n-1 panels); 150 dpi keeps the files a sensible size


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datasets", nargs="+", default=list(C.DATASETS))
    ap.add_argument("--algos", nargs="+", default=None)
    ap.add_argument("--include-unsound", action="store_true")
    ap.add_argument("--source", choices=["unified", "reference"], default="unified",
                    help="unified = this project's solvers (default); reference = the original reference binaries")
    args = ap.parse_args()
    plt = C.style()

    for ds in args.datasets:
        ref = args.source == "reference"
        algos = C.available(ds, include_unsound=args.include_unsound, validation=ref)
        if args.algos:
            algos = [k for k in algos if k in args.algos]
        if len(algos) < 2:
            print(f"{ds}: fewer than two algorithms, skipped"); continue
        cls = C.classify(ds)
        med = cls[cls["cls"] == "medium"][["a", "b"]]
        data = {k: med.merge(C.load(k, ds), on=["a", "b"]).set_index(["a", "b"]) for k in algos}
        if ref:                                   # only algorithms that have a reference run on this dataset
            algos = [k for k in algos if data[k]["ref_done"].any()]
            if len(algos) < 2:
                print(f"{ds}: fewer than two algorithms with a reference, skipped"); continue
        idx = data[algos[0]].index
        for k in algos[1:]:
            idx = idx.intersection(data[k].index)
        n = len(algos)
        print(f"\n{C.DATASETS[ds]}: {len(idx)} medium pairs, {n} algorithms")

        for col, axlabel, fname in METRICS:
            src_col = ("ref_" + col) if ref else col
            vals = {k: np.clip(data[k].loc[idx, src_col].astype(float).values, 1e-3, None) for k in algos}
            lo = min(np.nanmin(v) for v in vals.values()) * 0.7
            hi = max(np.nanmax(v) for v in vals.values()) * 1.4
            panel = 2.6
            fig, axes = plt.subplots(n - 1, n - 1, figsize=(panel * (n - 1) + 0.6, panel * (n - 1) + 0.6), squeeze=False)
            for r in range(n - 1):
                for c in range(n - 1):
                    ax = axes[r][c]
                    i, j = r, c + 1                      # y = algos[i], x = algos[j], needs i < j
                    if i >= j:
                        ax.axis("off"); continue
                    ax.scatter(vals[algos[j]], vals[algos[i]], s=6, alpha=0.25, linewidths=0, color="#1f77b4")
                    ax.plot([lo, hi], [lo, hi], ls="--", lw=0.8, color="#888888")
                    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
                    ax.set_box_aspect(1)
                    ax.set_title(f"{C.label(algos[j])} vs {C.label(algos[i])}", fontsize=10)
                    ax.tick_params(labelsize=8)
                    if j == n - 1 or r == c:             # bottom panel of each column = the diagonal
                        ax.set_xlabel(f"{axlabel} ({C.label(algos[j])})", fontsize=9)
                    if c == r:                           # first panel of each row
                        ax.set_ylabel(f"{axlabel} ({C.label(algos[i])})", fontsize=9)
            fig.tight_layout()
            p = C.out_dir(FOLDER) / f"pairwise_{fname}_{ds}{'_ref' if ref else ''}.png"
            fig.savefig(p, dpi=DPI, bbox_inches="tight")
            plt.close(fig)
            print("  saved", p)


if __name__ == "__main__":
    main()
