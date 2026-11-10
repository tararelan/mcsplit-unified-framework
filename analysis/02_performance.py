"""
02_performance.py - how do the algorithms compare on each dataset?        (dissertation Sections 4.2 - 4.5)

Per dataset (MIVIA, LV, BI):
  1. difficulty classes: easy / medium / hard                     (Table 7)   -> classification.csv
  2. instances solved on the shared medium set, re-implementation
     and reference, with ranks                                    (Table 8)   -> solved_medium.csv
  3. cactus plots on the medium set: time and nodes vs solved     (Fig 23-24) -> cactus_time_<ds>.png, cactus_nodes_<ds>.png (+ _log versions)
  4. time per node on the medium set (re-implementation, reference)(Table 9)  -> time_per_node.csv
  5. pairwise "found a larger incumbent" matrix on hard instances (Table 10)  -> hard_incumbent_<ds>.csv/.png

Difficulty is defined by the CORE algorithms in config.py (every algorithm evaluated), so the classes do not
change when a new algorithm is added; every algorithm with results is then evaluated on those classes.
A pair missing from an algorithm's result file counts as unsolved by that algorithm.

Usage:  python 02_performance.py [--datasets mivia lv bi] [--include-unsound]
Output: analysis/outputs/performance/
"""
import argparse

import numpy as np
import pandas as pd

import config as C

FOLDER = "performance"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datasets", nargs="+", default=list(C.DATASETS))
    ap.add_argument("--include-unsound", action="store_true")
    args = ap.parse_args()
    plt = C.style()

    cls_rows, solved_rows, tpn_rows = [], [], []
    for ds in args.datasets:
        name = C.DATASETS[ds]
        algos = C.available(ds, include_unsound=args.include_unsound)
        if not algos:
            print(f"{name}: no result files found"); continue
        print(f"\n===== {name} =====")
        C.print_sources(ds, algos)
        cls = C.classify(ds)
        core = cls.attrs.get("core", [])
        n = len(cls)
        counts = cls["cls"].value_counts()
        print(f"1. difficulty classes (core algorithms: {', '.join(C.label(k) for k in core)})")
        for c in ("easy", "medium", "hard"):
            cls_rows.append(dict(dataset=name, difficulty=c, count=int(counts.get(c, 0)),
                                 pct=round(100 * counts.get(c, 0) / n, 2)))
            print(f"   {c:7} {counts.get(c, 0):6}  ({100 * counts.get(c, 0) / n:5.2f}%)")

        medium = cls[cls["cls"] == "medium"][["a", "b"]]
        hard = cls[cls["cls"] == "hard"][["a", "b"]]
        nm = len(medium)
        med = {k: C.load(k, ds).merge(medium, on=["a", "b"], how="inner") for k in algos}

        # ---- 2. solved on the medium set
        rows = []
        for k in algos:
            d = med[k]
            r = dict(dataset=name, key=k, algorithm=C.label(k), medium_pairs=nm,
                     solved=int(d["done"].sum()), solved_pct=round(100 * d["done"].sum() / nm, 2))
            has_ref = k not in C.REF_BORROW and (k, ds) not in C.UNSOUND and d["ref_done"].any()
            r["ref_solved"] = int(d["ref_done"].sum()) if has_ref else np.nan
            r["ref_solved_pct"] = round(100 * d["ref_done"].sum() / nm, 2) if has_ref else np.nan
            rows.append(r)
        t = pd.DataFrame(rows)
        t["rank"] = t["solved"].rank(ascending=False, method="min").astype(int)
        t["ref_rank"] = t["ref_solved"].rank(ascending=False, method="min")
        solved_rows.append(t)
        print(f"\n2. instances solved on the shared medium set (n={nm})")
        print(t[["algorithm", "solved", "solved_pct", "rank", "ref_solved", "ref_solved_pct", "ref_rank"]]
              .sort_values("rank").to_string(index=False))

        # ---- 3. cactus plots
        for metric, col, xlabel, fname in [("time", "time", "wall-clock time (s)", "cactus_time"),
                                           ("nodes", "nodes", "search nodes", "cactus_nodes")]:
            for log in (False, True):             # linear and log-x version of every cactus plot
                fig, ax = plt.subplots(figsize=(10.5, 5.6))
                for i, k in enumerate(algos):
                    d = med[k]
                    v = np.sort(d.loc[d["done"], col].dropna().values)
                    if len(v) == 0:
                        continue
                    a = C.ALGOS[k]
                    ax.plot(np.clip(v, 1e-3 if col == "time" else 1, None), np.arange(1, len(v) + 1),
                            color=a.color, ls=a.ls, marker=a.marker, lw=1.0, alpha=0.9, ms=3.2, mew=0,
                            label=f"{a.label} ({len(v)})",
                            markevery=(int(len(v) * i / (len(algos) * 8)), max(1, len(v) // 8)))
                if log:
                    ax.set_xscale("log")
                ax.set_ylim(bottom=0)
                ax.set_xlabel(xlabel + (" (log scale)" if log else "")); ax.set_ylabel("instances solved (medium set)")
                ax.set_title(f"{name}: {xlabel} vs instances solved, medium instances (n={nm})")
                ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), title="algorithm (solved)", title_fontsize=8.5)
                fig.tight_layout(); C.save(fig, FOLDER, f"{fname}{'_log' if log else ''}_{ds}.png"); plt.close(fig)

        # ---- 4. time per node on the medium set
        print(f"\n4. time per node on the medium set (microseconds per node)")
        rows = []
        for k in algos:
            d = med[k]
            u = d[d["done"] & (d["nodes"] > 0)]
            upn = u["time"] / u["nodes"] * 1e6
            r = dict(dataset=name, key=k, algorithm=C.label(k), n=len(u),
                     median_us=round(upn.median(), 3), mean_us=round(upn.mean(), 3))
            if k not in C.REF_BORROW and (k, ds) not in C.UNSOUND and d["ref_done"].any():
                q = d[d["ref_done"] & (d["ref_nodes"] > 0) & (d["ref_time"] > 0)]
                rpn = q["ref_time"] / q["ref_nodes"] * 1e6
                r.update(ref_n=len(q), ref_median_us=round(rpn.median(), 3), ref_mean_us=round(rpn.mean(), 3))
            rows.append(r)
        t = pd.DataFrame(rows)
        tpn_rows.append(t)
        print(t.drop(columns=["dataset", "key"]).to_string(index=False))

        # ---- 5. pairwise larger incumbent on hard instances
        if len(hard):
            hd = {k: C.load(k, ds).merge(hard, on=["a", "b"], how="inner").set_index(["a", "b"])["size"] for k in algos}
            m = pd.DataFrame(index=[C.label(k) for k in algos], columns=[C.label(k) for k in algos], dtype=float)
            for r_ in algos:
                for c_ in algos:
                    if r_ == c_:
                        continue
                    j = pd.concat([hd[r_], hd[c_]], axis=1, join="inner").dropna()
                    m.loc[C.label(r_), C.label(c_)] = 100 * (j.iloc[:, 0] > j.iloc[:, 1]).mean() if len(j) else np.nan
            m.round(1).to_csv(C.out_dir(FOLDER) / f"hard_incumbent_{ds}.csv")
            print(f"\n5. % of hard instances (n={len(hard)}) where the ROW algorithm found a larger incumbent than the COLUMN algorithm")
            print(m.round(0).to_string(float_format=lambda x: f"{x:.0f}"))
            fig, ax = plt.subplots(figsize=(0.75 * len(algos) + 2.5, 0.65 * len(algos) + 1.8))
            from matplotlib.colors import LinearSegmentedColormap
            cmap = LinearSegmentedColormap.from_list("blue", ["#fcfcfb", "#9ec5f4", "#256abf", "#0d366b"])
            ax.imshow(m.values.astype(float), cmap=cmap, vmin=0, vmax=30)   # capped so the top-down columns (see README) do not wash out the rest
            ax.set_xticks(range(len(algos))); ax.set_xticklabels(m.columns, rotation=45, ha="right")
            ax.set_yticks(range(len(algos))); ax.set_yticklabels(m.index)
            ax.grid(False)
            for i in range(len(algos)):
                for j in range(len(algos)):
                    v = m.values[i, j]
                    if not np.isnan(v):
                        ax.text(j, i, f"{v:.0f}", ha="center", va="center", fontsize=8,
                                color="white" if v > 17 else C.INK)
            ax.set_title(f"{name}: % of hard instances where row found a larger incumbent than column (n={len(hard)})", fontsize=9.5)
            fig.tight_layout(); C.save(fig, FOLDER, f"hard_incumbent_{ds}.png"); plt.close(fig)

    out = C.out_dir(FOLDER)
    pd.DataFrame(cls_rows).to_csv(out / "classification.csv", index=False)
    if solved_rows:
        pd.concat(solved_rows).to_csv(out / "solved_medium.csv", index=False)
        pd.concat(tpn_rows).to_csv(out / "time_per_node.csv", index=False)
    print("\nTables saved to", out.relative_to(C.REPO))


if __name__ == "__main__":
    main()
