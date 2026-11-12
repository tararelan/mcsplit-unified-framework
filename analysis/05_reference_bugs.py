"""
05_reference_bugs.py - evidence for the reference-implementation defects      (dissertation Section 4.1)

  DAL       solution-inflating reward bug: the reference reports a larger common subgraph than exists.
  RRSplit   on the directed BI instances the reference writes no result at all (it fails) on almost every pair.
  SymSplit  on BI the reference prints a size that is wrong in BOTH directions (too large and too small).

The three defects look different, so every algorithm x dataset is sorted into the same four outcomes. For every
pair that plain McSplit (exact, unified) solves, McSplit's size is the true optimum. The reference of the algorithm
under test then either
  - gives no output   (size and node count blank: it failed, or timed out; the McSplit reference does this on
                       well under 1% of the same pairs, so a high share points to a failure),
  - agrees            (reference size = optimum),
  - reports MORE      (impossible for an exact solver: inflation, the DAL signature),
  - reports LESS      (a finished exact solver cannot undershoot either).
The unified re-implementation of the same algorithm is checked against the same optimum.

MIVIA and LV are undirected, so RRSplit/SymSplit are expected to behave there: they are the control that shows
the defect is the direction handling.

Output (analysis/outputs/reference_bugs/):
  reference_bugs.csv      one row per algorithm x dataset: pairs, then no output / = / > / < the optimum (count and %),
                          the largest excess and shortfall, and the unified implementation's agreement with the optimum
  reference_bugs.tex      the same table for the paper
  reference_outcomes_<k>.png  one image per algorithm (dal, rrsplit, symsplit): stacked bars of the four outcomes,
                          one bar per dataset, exact counts beside each bar
  reference_bugs.png      reference size against true size (pairs where the reference produced output); points
                          above the diagonal are inflated, below are undershoots

Usage:  python 05_reference_bugs.py [--truth mcsplit]
"""
import argparse

import numpy as np
import pandas as pd

import config as C

FOLDER = "reference_bugs"
CASES = ["dal", "rrsplit", "symsplit"]        # algorithms whose reference is under test
DATASETS = ["mivia", "lv", "bi"]
OUTCOMES = ["equal", "larger", "smaller", "no output"]
OUTCOME_COLORS = {"equal": "#a9c7b8", "larger": "#e07a5f", "smaller": "#5b7fa6", "no output": "#dcdcd8"}


def compare(k, ds, truth_algo):
    """Per-pair table for one algorithm on one dataset: truth, reference size, unified size. None if no data."""
    d, t = C.load(k, ds), C.load(truth_algo, ds)
    if d is None or t is None:
        return None
    t = t[t["done"]][["a", "b", "size"]].rename(columns={"size": "truth"})
    m = d.merge(t, on=["a", "b"], how="inner").copy()
    m["excess"] = m["ref_size"] - m["truth"]          # NaN where the reference gave no output
    m["outcome"] = np.select([m["ref_size"].isna(), m["excess"] > 0, m["excess"] < 0],
                             ["no output", "larger", "smaller"], default="equal")
    return m


def summarise(k, ds, m):
    n = len(m)
    cnt = m["outcome"].value_counts()
    c = {o: int(cnt.get(o, 0)) for o in OUTCOMES}
    pct = lambda o: round(100 * c[o] / n, 1) if n else np.nan
    u = m[m["done"]]
    return dict(algorithm=C.label(k), dataset=C.DATASETS[ds], pairs=n,
                ref_no_output=c["no output"], ref_equal=c["equal"], ref_larger=c["larger"], ref_smaller=c["smaller"],
                pct_no_output=pct("no output"), pct_equal=pct("equal"), pct_larger=pct("larger"),
                pct_smaller=pct("smaller"),
                max_excess=int(m["excess"].max()) if c["larger"] else 0,
                max_shortfall=int(-m["excess"].min()) if c["smaller"] else 0,
                unified_pairs=len(u), unified_equal=int((u["size"] == u["truth"]).sum()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth", default="mcsplit", help="exact algorithm whose size is taken as the true optimum")
    args = ap.parse_args()
    plt = C.style()

    tables, rows = {}, []
    for k in CASES:
        for ds in DATASETS:
            m = compare(k, ds, args.truth)
            if m is None or m.empty:
                print(f"  [{ds}] {k}: no reference results, skipped"); continue
            tables[(k, ds)] = m
            rows.append(summarise(k, ds, m))
    if not rows:
        print("nothing to report"); return
    t = pd.DataFrame(rows)
    out = C.out_dir(FOLDER)
    t.to_csv(out / "reference_bugs.csv", index=False)
    print("\nreference size vs the true optimum (reference finished and McSplit finished):")
    print(t.to_string(index=False))

    tex = t[["algorithm", "dataset", "pairs", "pct_no_output", "pct_equal", "pct_larger", "pct_smaller",
             "max_excess", "max_shortfall", "unified_equal", "unified_pairs"]].rename(columns={
        "pct_no_output": "no output \\%", "pct_equal": "$=$ opt \\%", "pct_larger": "$>$ opt \\%",
        "pct_smaller": "$<$ opt \\%", "max_excess": "max excess", "max_shortfall": "max shortfall",
        "unified_equal": "unif. $=$ opt", "unified_pairs": "unif. pairs"})
    (out / "reference_bugs.tex").write_text(tex.to_latex(index=False, escape=False, float_format="%.1f"))
    print("  saved", (out / "reference_bugs.tex").relative_to(C.REPO))

    ks_ = [k for k in CASES if any((k, ds) in tables for ds in DATASETS)]
    # ---- figure 1: the four outcomes as stacked bars, one image per algorithm, one bar per dataset
    for k in ks_:
        tk = t[t["algorithm"] == C.label(k)].reset_index(drop=True)
        fig, ax = plt.subplots(figsize=(11, 0.8 * len(tk) + 2.2))
        left = np.zeros(len(tk))
        names = {"equal": "matches optimum", "larger": "too large", "smaller": "too small", "no output": "no output"}
        for o in OUTCOMES:
            v = tk["pct_" + o.replace(" ", "_")].fillna(0).values
            ax.barh(range(len(tk)), v, left=left, height=0.6, color=OUTCOME_COLORS[o], edgecolor="none",
                    zorder=3, label=names[o])
            for i, (x0, w) in enumerate(zip(left, v)):
                if w >= 6:
                    ax.text(x0 + w / 2, i, f"{w:.0f}%", ha="center", va="center", fontsize=C.FONT, zorder=4,
                            color="white" if o in ("larger", "smaller") else C.INK)
            left += v
        ax.set_yticks(range(len(tk)))
        ax.set_yticklabels([f"{d}  ({n:,})" for d, n in zip(tk["dataset"], tk["pairs"])])
        ax.invert_yaxis(); ax.set_xlim(0, 100); ax.set_ylim(len(tk) - 0.5, -0.5)
        ax.set_xticks(range(0, 101, 25)); ax.set_xticklabels([f"{x}%" for x in range(0, 101, 25)])
        ax.set_xlabel("% of pairs")
        ax.grid(False)
        for sp in ("left", "right", "top", "bottom"):
            ax.spines[sp].set_visible(False)
        ax.tick_params(axis="both", length=0)
        ax.set_title(f"{C.label(k)} reference vs true optimum", loc="left", pad=40)
        h, l = ax.get_legend_handles_labels()
        ax.legend(h, l, loc="lower left", bbox_to_anchor=(-0.02, 1.0), ncol=4, frameon=False,
                  handlelength=1.0, handleheight=1.0, columnspacing=1.4)
        fig.tight_layout(); C.save(fig, FOLDER, f"reference_outcomes_{k}.png"); plt.close(fig)

    ks = ks_
    # ---- figure 2: reference size vs true size, a row per algorithm, a column per dataset
    ds_ = [ds for ds in DATASETS if any((k, ds) in tables for k in ks)]
    fig, axes = plt.subplots(len(ks), len(ds_), figsize=(6.2 * len(ds_) + 1.0, 5.8 * len(ks) + 0.6), squeeze=False)
    rng = np.random.default_rng(0)
    for r, k in enumerate(ks):
        for c, ds in enumerate(ds_):
            ax = axes[r][c]
            m = tables.get((k, ds))
            if m is None:
                ax.axis("off"); continue
            m = m[m["outcome"] != "no output"]            # only pairs where the reference printed a size
            bad = m["excess"] != 0
            jx, jy = rng.uniform(-0.18, 0.18, (2, len(m)))   # sizes are integers; jitter so the points do not stack
            ax.scatter((m["truth"] + jx)[~bad], (m["ref_size"] + jy)[~bad], s=40, color=C.INK2, alpha=0.35,
                       edgecolors="none", label="reference = optimum")
            ax.scatter((m["truth"] + jx)[bad], (m["ref_size"] + jy)[bad], s=170 if C.marker_of(k) == "*" else 80,
                       color=C.ALGOS[k].color, marker=C.marker_of(k), alpha=0.8, edgecolors=C.SURFACE, linewidths=0.8,
                       label="reference differs from optimum")
            lim = [0, (max(m["truth"].max(), m["ref_size"].max()) if len(m) else 1) + 1]
            ax.plot(lim, lim, color=C.INK2, ls="--", lw=1.8, zorder=0)
            ax.set_xlim(lim); ax.set_ylim(lim)
            ax.set_box_aspect(1)
            full = tables[(k, ds)]
            ax.set_title(f"{C.label(k)} on {C.DATASETS[ds]}\n{int(bad.sum())} of {len(m)} differ; "
                         f"{int((full['outcome'] == 'no output').sum())} of {len(full)} no output")
            if r == len(ks) - 1:
                ax.set_xlabel("true optimum (McSplit)")
            if c == 0:
                ax.set_ylabel("reference size")
    h, l = axes[0][0].get_legend_handles_labels()
    if h:
        fig.legend(h, l, loc="lower center", ncol=2, bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout()
    C.save(fig, FOLDER, "reference_bugs.png")


if __name__ == "__main__":
    main()
