"""
03_ablation.py - what does each algorithm component contribute?          (dissertation Section 4.3)

For each parent algorithm the full algorithm is compared with its variants on the shared medium set of each
dataset (medium = defined in config.py), by instances solved, plus the counters that explain the result.

  McSplit   plain McSplit vs McSplit-down
  RL        McSplit+RL vs its Pair / Rate / Down variants
  LL        LL vs LSM-only, LUM-only
  DAL       DAL vs its RL-only / DAL-only policy
  DSB       DSB vs never / always using the bound
  RRSplit   RRSplit vs no maximality / no bound / no vertex equivalence
  SymSplit  SymSplit vs variable-only / value-only symmetry

McSplit and RL variants are complete algorithms and come from the normal results. LL / DSB / RRSplit /
SymSplit variants come from hpc/publication ({variant}_{ds}_all.csv); a dataset with no ablation file for
a variant simply has no point for it. RRSplit and SymSplit are skipped on BI (unsound reference).

Output:
  ablation_solved.csv       one row per dataset x algorithm: solved %, change vs the parent, mean bound/symmetry
                            pruned branches (millions, over every instance solved), median time per node
  ablation_dotplot_<ds>.png one file per dataset: solved % on the medium set grouped by parent;
                            filled = with the component, open = without; dashed line = McSplit   (Fig 25)
  ablation_family_<parent>.png  one image per family (mcsplit, rl, ll, dal, dsb, rrsplit, symsplit), one panel per
                            dataset side by side, a row per variant; x = contribution of the component in percentage points (--figure family or both)

Usage:  python 03_ablation.py [--datasets mivia lv bi] [--figure dataset|family|both]
Output: analysis/outputs/ablation/
"""
import argparse

import pandas as pd
from matplotlib.ticker import FuncFormatter

import config as C

FOLDER = "ablation"
FONT = C.FONT        # text size and DPI come from config.py (FONT, SAVE_DPI), shared by every script
ADDITIVE = ("mcsplit", "rl")      # parent stripped, variants added on top (filled dots); the others ablate the full algorithm
PARENTS = [p for p in ["mcsplit", "rl", "ll", "dal", "dsb", "rrsplit", "symsplit"] if p in C.ALGOS]   # "dal" only while its ALGOS line is not commented out


def load_variant(key, ds):
    """Results of one ablation variant on one dataset, or None."""
    if key in C.ALGOS:                                   # a full algorithm (rl_pair, down, ...)
        return C.load(key, ds)
    for tmpl in C.ABLATION_SOURCES:
        if tmpl.count("{ds}") == 0 and ds != "mivia":    # dataset-less names are the MIVIA run
            continue
        path = C.REPO / tmpl.format(a=key, ds=ds)
        if path.exists():
            return C._standardise(C._read_csv(path), key)
    return None


def ablation_table(ds):
    cls = C.classify(ds)
    medium = cls[cls["cls"] == "medium"][["a", "b"]]
    nm = len(medium)

    def row(key, lab, parent, kind, whole):
        d = whole.merge(medium, on=["a", "b"], how="inner")
        solved_all = whole[whole["done"]]                # every instance solved, not only medium
        u = d[d["done"] & (d["nodes"] > 0)]
        return dict(dataset=ds, parent=C.label(parent), kind=kind, key=key, algorithm=lab, medium_pairs=nm,
                    solved=int(d["done"].sum()), solved_pct=round(100 * d["done"].sum() / nm, 2),
                    bound_pruned_mean_M=round(solved_all["bound_pruned"].mean() / 1e6, 2),
                    sym_pruned_mean_M=round(solved_all["sym_pruned"].mean() / 1e6, 2),
                    median_us_per_node=round((u["time"] / u["nodes"] * 1e6).median(), 3))

    rows = []
    for parent in PARENTS:
        # UNSOUND concerns the REFERENCE binaries; the ablation uses only the unified results, so no skip here
        whole = C.load(parent, ds)
        if whole is None:
            print(f"  [{ds}] skipping {parent}: no results"); continue
        rows.append(row(parent, C.label(parent) + (" (plain)" if parent == "mcsplit" else " (full)"),
                        parent, "full", whole))
        for key, lab, p in C.VARIANTS:
            if p != parent:
                continue
            v = load_variant(key, ds)
            if v is None:
                print(f"  [{ds}] no results for variant {key}"); continue
            rows.append(row(key, lab, parent, "variant", v))
    t = pd.DataFrame(rows)
    if len(t):
        full = t[t["kind"] == "full"].set_index("parent")["solved_pct"]
        t["change_vs_parent_pp"] = (t["solved_pct"] - t["parent"].map(full)).round(2)
    return t, nm


def plot_by_dataset(t, sizes, plt):
    # ---- Fig 25: one file per dataset. In each group the filled dots are the algorithm WITH the component
    # (top rows), the open dot is the version without it. For McSplit and RL the parent is the stripped
    # version and the variants are what is added on top of it; for LL/DSB/RRSplit/SymSplit the parent is the
    # full algorithm and the variants are ablations of it.
    for ds in sizes:
        td = t[t["dataset"] == ds]
        groups = []                                       # list of rows (key, label, filled), top to bottom
        for p in PARENTS:
            sub = td[td["parent"] == C.label(p)]
            if sub.empty:
                continue
            par = sub[sub["kind"] == "full"].iloc[0]
            var = sub[sub["kind"] == "variant"]
            nice = lambda r: C.label(r["key"]) if r["key"] in C.ALGOS else r["algorithm"]
            if p in ADDITIVE:
                rows_ = [(r, nice(r), True) for _, r in var.iterrows()] + [(par, C.label(p), False)]
            else:
                rows_ = [(par, C.label(p), True)] + [(r, nice(r), False) for _, r in var.iterrows()]
            groups.append(rows_)
        n_rows = sum(len(g) for g in groups)
        fig, ax = plt.subplots(figsize=(13, 0.5 * n_rows + 0.6 * len(groups) + 2.4))
        y, ys, labels = 0.0, [], []
        mc = td[td["key"] == "mcsplit"]["solved_pct"]
        for gi, g in enumerate(groups):
            if gi:                                         # faint separator between groups
                ax.axhline(y - 0.9, color=C.GRID, lw=1.0, zorder=0)
            for r, lab, filled in g:
                color = C.color_of(r["key"])
                ax.plot(r["solved_pct"], y, marker=C.marker_of(r["key"]), ms=11, mfc=color if filled else C.SURFACE,
                        mec=color, mew=2.4, ls="none")
                ax.annotate(f"{r['solved_pct']:.1f}", (r["solved_pct"], y), xytext=(0, 9), textcoords="offset points",
                            ha="center", va="bottom", fontsize=FONT, color=C.INK2)   # a fixed gap above the marker
                ys.append(y); labels.append(lab); y += 1
            y += 0.8
        if len(mc):
            ax.axvline(mc.iloc[0], color=C.INK2, ls="--", lw=1.8, label=f"McSplit ({mc.iloc[0]:.1f}%)")
            ax.legend(loc="lower right")
        ax.set_yticks(ys); ax.set_yticklabels(labels); ax.invert_yaxis()
        ax.grid(axis="y", visible=False); ax.margins(x=0.08)
        ax.set_xlabel(f"instances solved on the medium set (%), {C.DATASETS[ds]}, n={sizes[ds]}")
        ax.set_title(f"{C.DATASETS[ds]} ablation: solved % by component\n"
                     "(filled = with the component, open = without it; marker shape = algorithm family)", fontsize=FONT)
        fig.tight_layout(); C.save(fig, FOLDER, f"ablation_dotplot_{ds}.png"); plt.close(fig)


def contribution(parent, r):
    """Contribution of the component in percentage points: positive = the component helps.
    Additive parents (McSplit, RL) are stripped and the variant is added on top; otherwise the variant removes it."""
    return r["change_vs_parent_pp"] if parent in ADDITIVE else -r["change_vs_parent_pp"]


def plot_by_family(t, plt):
    """One image per parent family (ablation_family_<parent>.png), with one panel per dataset side by side
    (MIVIA | LV | BI). Rows = the family's variants, shared across the panels. x = contribution of the component in
    percentage points (positive = it helps). Colour = variant, marker shape = family, as in the per-dataset
    figure. Filled = component added to a stripped parent (McSplit, RL); open = removed from the full algorithm.
    Every image uses the same x-range, so families can be compared. A variant with no result on a dataset has no dot."""
    datasets = [d for d in ("mivia", "lv", "bi") if d in set(t["dataset"])]
    panels = {}
    for p in PARENTS:
        full = t[(t["kind"] == "full") & (t["key"] == p)]
        var = t[(t["parent"] == C.label(p)) & (t["kind"] == "variant")]
        order = list(dict.fromkeys(var["key"]))
        if not full.empty and order:
            panels[p] = (full, var, order)
    if not panels:
        return
    xs = [contribution(p, r) for p, (_, v, _) in panels.items() for _, r in v.iterrows()] + [0.0]
    pad = 0.1 * (max(xs) - min(xs)) + 1
    xlim = (min(xs) - pad, max(xs) + pad)

    for p, (full, var, order) in panels.items():
        fig, axes = plt.subplots(1, len(datasets), sharey=True, squeeze=False,
                                 figsize=(5.6 * len(datasets) + 3.6, 1.0 * len(order) + 4.2))
        labs = [var[var["key"] == k]["algorithm"].iloc[0] for k in order]
        for ax, ds in zip(axes[0], datasets):
            ax.axvline(0, color=C.INK2, ls="--", lw=1.8, zorder=1)
            for yi, k in enumerate(order):
                r = var[(var["key"] == k) & (var["dataset"] == ds)]
                if r.empty:
                    continue
                x = contribution(p, r.iloc[0])
                col, mk = C.color_of(k), C.marker_of(k)
                ax.plot(x, yi, marker=mk, ms=17 if mk == "*" else 13, mfc=col if p in ADDITIVE else C.SURFACE,
                        mec=col, mew=2.4, ls="none", zorder=3)
                ax.annotate(f"{x:+.1f}%", (x, yi), xytext=(0, 13 if mk == "*" else 10), textcoords="offset points",
                            ha="center", va="bottom", fontsize=FONT, color=C.INK2, zorder=5,
                            bbox=dict(boxstyle="round,pad=0.12", fc=C.SURFACE, ec="none"))   # fixed gap above the marker
            ax.set_ylim(len(order) - 0.5, -0.5)
            ax.set_yticks(range(len(order))); ax.set_yticklabels(labs)
            ax.set_xlim(*xlim); ax.grid(axis="y", visible=False)
            f = full[full["dataset"] == ds]
            ax.set_title(C.DATASETS[ds] + (f"\n{C.label(p)}: {f['solved_pct'].iloc[0]:.1f}% solved" if len(f) else ""),
                         fontsize=FONT)
            ax.set_xlabel("change in % of instances solved\n(percentage points, vs " + ("the stripped parent" if p in ADDITIVE else "the full algorithm") + ")")
            ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:+.0f}%" if v else "0"))
        fig.suptitle(f"{C.label(p)} ablation: contribution of each component\n"
                     f"({'added to the stripped parent' if p in ADDITIVE else 'removed from the full algorithm'}; "
                     "right of the dashed line = the component helps)", x=0.01, ha="left", fontsize=FONT)
        fig.tight_layout(); C.save(fig, FOLDER, f"ablation_family_{p}.png"); plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datasets", nargs="+", default=list(C.DATASETS), choices=list(C.DATASETS))
    ap.add_argument("--figure", choices=["dataset", "family", "both"], default="both",
                    help="dataset = one figure per dataset (Fig 25); family = one image per family with all "
                         "datasets on each row; both (default)")
    args = ap.parse_args()
    plt = C.style()

    tables, sizes = [], {}
    for ds in args.datasets:
        if C.load("mcsplit", ds) is None:
            print(f"{ds}: no McSplit results, skipping"); continue
        t, nm = ablation_table(ds)
        if t.empty:
            continue
        sizes[ds] = nm
        tables.append(t)
        print(f"\n{C.DATASETS[ds]}: medium set {nm} pairs")
        print(t[["parent", "algorithm", "solved", "solved_pct", "change_vs_parent_pp", "bound_pruned_mean_M",
                 "sym_pruned_mean_M", "median_us_per_node"]].to_string(index=False))
    if not tables:
        return
    t = pd.concat(tables, ignore_index=True)
    t.to_csv(C.out_dir(FOLDER) / "ablation_solved.csv", index=False)

    if args.figure in ("dataset", "both"):
        plot_by_dataset(t, sizes, plt)
    if args.figure in ("family", "both"):
        plot_by_family(t, plt)


if __name__ == "__main__":
    main()
