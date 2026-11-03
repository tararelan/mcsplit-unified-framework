"""
04_bidomain.py - how many bidomains does the search create per node, and on which datasets?

Bidomains per search node (mean and maximum over the whole search) are recorded by the current pipeline
(hpc/publication, 21-column files). A "bidomain" is a group of still-unmatched vertex pairs that are
compatible with each other (Figure 13 of the dissertation); more bidomains per node means more tie-breaking
decisions for an implementation to resolve.

Outputs:
  bidomain_summary.csv         per dataset (baseline algorithm): pairs, mean/median/p90/max of the per-pair mean and max
  bidomain_mean_hist.png       histogram of each pair's MEAN bidomains per node, one panel per dataset
  bidomain_max_hist.png        same for the MAX (log x axis, and a linear version clipped at the 98th percentile)
  bidomain_vs_size.png         mean bidomains vs solution size (log-log, per-bin medians) - is it just instance size?
  bidomain_by_algorithm.csv/.png   median mean-bidomains per node for every algorithm and dataset

Usage:  python 04_bidomain.py [--baseline mcsplit] [--solved-only]
Output: analysis/outputs/bidomain/
"""
import argparse

import numpy as np
import pandas as pd

import config as C

FOLDER = "bidomain"


def pairs_with_stats(k, ds, solved_only):
    d = C.load(k, ds, need_bidomain=True)
    if d is None:
        return None
    d = d[(d["bd_mean"] > 0) & (d["nodes"] > 0)]
    return d[d["done"]] if solved_only else d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--baseline", default="mcsplit", help="algorithm used for the dataset comparison")
    ap.add_argument("--solved-only", action="store_true", help="only pairs the algorithm solved (default: all pairs)")
    args = ap.parse_args()
    plt = C.style()
    base = args.baseline
    ds_colors = {"mivia": "#2a78d6", "lv": "#1baf7a", "bi": "#eb6834"}

    frames = {ds: pairs_with_stats(base, ds, args.solved_only) for ds in C.DATASETS}
    frames = {ds: d for ds, d in frames.items() if d is not None and len(d)}
    if not frames:
        raise SystemExit(f"No bidomain statistics found for {C.label(base)}.")
    print(f"Baseline algorithm: {C.label(base)}  ({'solved pairs only' if args.solved_only else 'all pairs'})")
    for ds, d in frames.items():
        print(f"  {C.DATASETS[ds]:6} {len(d):6} pairs  <- {d.attrs['source']}")

    # ---- summary table
    rows = []
    for ds, d in frames.items():
        rows.append(dict(dataset=C.DATASETS[ds], pairs=len(d),
                         mean_of_means=round(d["bd_mean"].mean(), 2), median_of_means=round(d["bd_mean"].median(), 2),
                         p90_of_means=round(d["bd_mean"].quantile(0.9), 2),
                         mean_of_max=round(d["bd_max"].mean(), 1), median_of_max=d["bd_max"].median(),
                         p90_of_max=d["bd_max"].quantile(0.9), max_of_max=d["bd_max"].max()))
    t = pd.DataFrame(rows)
    t.to_csv(C.out_dir(FOLDER) / "bidomain_summary.csv", index=False)
    print("\nBidomains per search node:\n" + t.to_string(index=False))

    # ---- histograms
    def hist_by_dataset(col, title, fname, log_x=False, clip_pct=None):
        fig, axes = plt.subplots(1, len(frames), figsize=(4.6 * len(frames), 3.8), sharey=False)
        axes = np.atleast_1d(axes)
        for ax, (ds, d) in zip(axes, frames.items()):
            v = d[col].values.astype(float)
            dropped = 0
            if clip_pct:
                cut = np.percentile(v, clip_pct)
                dropped = int((v > cut).sum()); v = v[v <= cut]
            bins = np.logspace(np.log10(max(v.min(), 1e-2)), np.log10(v.max() + 1e-9), 40) if log_x \
                else np.linspace(v.min(), v.max(), 40)
            ax.hist(v, bins=bins, color=ds_colors.get(ds, "grey"), alpha=0.9, edgecolor=C.SURFACE, linewidth=0.5)
            ax.axvline(np.median(v), color=C.INK2, ls="--", lw=1)
            if log_x:
                ax.set_xscale("log")
            ax.set_xlabel(title); ax.set_ylabel("pairs")
            ax.set_title(f"{C.DATASETS[ds]}  (median {np.median(v):.2f}, n={len(v)}" + (f", {dropped} above p{clip_pct} hidden)" if dropped else ")"), fontsize=9.5)
        fig.suptitle(f"{C.label(base)}: {title} per pair")
        fig.tight_layout(); C.save(fig, FOLDER, fname); plt.close(fig)

    hist_by_dataset("bd_mean", "mean bidomains per node", "bidomain_mean_hist.png", clip_pct=98)
    hist_by_dataset("bd_max", "max bidomains at any node", "bidomain_max_hist.png", log_x=True)

    # ---- is it just instance size?
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for ds, d in frames.items():
        s = d[(d["size"] > 0)]
        ax.scatter(s["size"], s["bd_mean"], s=7, alpha=0.2, color=ds_colors[ds], edgecolors="none")
        edges = np.logspace(0, np.log10(s["size"].max() + 1), 15)
        mid = s.groupby(pd.cut(s["size"], edges), observed=True).agg(
            x=("size", "median"), y=("bd_mean", "median"), n=("bd_mean", "size"))
        mid = mid[mid["n"] >= 5]
        ax.plot(mid["x"], mid["y"], color=ds_colors[ds], lw=2.2, marker="o", ms=4, label=C.DATASETS[ds])
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("solution size (vertices)"); ax.set_ylabel("mean bidomains per node")
    ax.set_title(f"{C.label(base)}: bidomains vs solution size (lines = per-size-bin medians)")
    ax.legend()
    fig.tight_layout(); C.save(fig, FOLDER, "bidomain_vs_size.png"); plt.close(fig)

    # ---- every algorithm
    rows = []
    for ds in C.DATASETS:
        for k in C.ALL_ALGOS:
            d = pairs_with_stats(k, ds, args.solved_only)
            if d is not None and len(d) and not (k, ds) in C.UNSOUND:
                rows.append(dict(dataset=C.DATASETS[ds], key=k, algorithm=C.label(k), pairs=len(d),
                                 median_mean_bidomains=round(d["bd_mean"].median(), 3),
                                 median_max_bidomains=d["bd_max"].median()))
    if rows:
        a = pd.DataFrame(rows)
        a.to_csv(C.out_dir(FOLDER) / "bidomain_by_algorithm.csv", index=False)
        print("\nMedian mean-bidomains per node, by algorithm:")
        print(a.pivot(index="algorithm", columns="dataset", values="median_mean_bidomains").to_string())
        fig, ax = plt.subplots(figsize=(8.2, 0.42 * a["algorithm"].nunique() + 1.8))
        order = [C.label(k) for k in C.ALL_ALGOS if C.label(k) in set(a["algorithm"])]
        for i, name in enumerate(order):
            for j, ds in enumerate(C.DATASETS):
                r = a[(a["algorithm"] == name) & (a["dataset"] == C.DATASETS[ds])]
                if len(r):
                    ax.plot(r["median_mean_bidomains"].iloc[0], i + (j - 1) * 0.2, "o", ms=6, color=ds_colors[ds],
                            label=C.DATASETS[ds] if i == 0 else None)
        ax.set_yticks(range(len(order))); ax.set_yticklabels(order); ax.invert_yaxis()
        ax.grid(axis="y", visible=False)
        ax.set_xlabel("median of per-pair mean bidomains per node"); ax.legend(title="dataset")
        ax.set_title("Bidomains per node by algorithm")
        fig.tight_layout(); C.save(fig, FOLDER, "bidomain_by_algorithm.png"); plt.close(fig)


if __name__ == "__main__":
    main()
