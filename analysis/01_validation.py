"""
01_validation.py - does each re-implementation reproduce its reference?   (dissertation Section 4.1)

What it does, per dataset:
  A. PASS / FAIL / timed-out table                       (Table 6)   -> validation_classification.csv
  B. search-node fidelity: R^2, exact match, within 2x   (Table 3)   -> validation_fidelity.csv
  C. runtime fidelity: median unified/reference ratio    (Table 5)   -> validation_fidelity.csv
  D. figures: nodes and runtime, unified vs reference    (Fig 17-22) -> nodes_<ds>.png, runtime_<ds>.png

rl_pair / rl_rate / rl_down have no reference binary of their own, so they are not part of the validation
(they are judged in 02_performance and 03_ablation). McSplit-down is validated: its reference is the McSplit
binary run with --big-first.

Usage:  python 01_validation.py [--datasets mivia lv bi] [--include-unsound]
Output: analysis/outputs/validation/
"""
import argparse

import numpy as np
import pandas as pd

import config as C

FOLDER = "validation"


def complete_pairs(d):
    """Pairs where both programs finished, sizes agree, and both counted nodes."""
    return d[d["done"] & d["ref_done"] & (d["match"] == "PASS") & (d["nodes"] > 0) & (d["ref_nodes"] > 0)]


def r_squared(x, y):
    if len(x) < 3 or np.log10(x).std() == 0 or np.log10(y).std() == 0:
        return np.nan
    return float(np.corrcoef(np.log10(x), np.log10(y))[0, 1] ** 2)


def ref_note(k):
    return f" (vs {C.label(C.REF_BORROW[k])} ref.)" if k in C.REF_BORROW else ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datasets", nargs="+", default=list(C.DATASETS))
    ap.add_argument("--include-unsound", action="store_true",
                    help="also show RRSplit/SymSplit on BI (their references are unsound on directed graphs)")
    args = ap.parse_args()
    plt = C.style()

    class_rows, fid_rows = [], []
    for ds in args.datasets:
        name = C.DATASETS[ds]
        # only algorithms with a reference binary of their own: RL's variants (Pair/Rate/Down) have none,
        # so a comparison with RL's reference is not a fidelity check and is left out
        algos = [k for k in C.available(ds, include_unsound=args.include_unsound) if k not in C.REF_BORROW]
        frames = {k: C.load(k, ds) for k in algos}
        frames = {k: d for k, d in frames.items() if d is not None}
        if not frames:
            print(f"{name}: no result files found"); continue
        print(f"\n===== {name} =====")
        C.print_sources(ds, algos)

        # ---- A + B + C: tables
        for k, d in frames.items():
            vc = d["match"].value_counts()
            class_rows.append(dict(dataset=name, algorithm=C.label(k), key=k, pairs=len(d),
                                   PASS=int(vc.get("PASS", 0)), FAIL=int(vc.get("FAIL", 0)),
                                   timed_out_or_missing=int(vc.get("MISSING", 0)),
                                   pass_pct=round(100 * vc.get("PASS", 0) / len(d), 2),
                                   reference="own" if k not in C.REF_BORROW else f"{C.REF_BORROW[k]} (borrowed)"))
            p = complete_pairs(d)
            if p.empty:
                continue
            ratio_n = np.maximum(p["nodes"], p["ref_nodes"]) / np.minimum(p["nodes"], p["ref_nodes"])
            t = p[p["ref_time"] > 0]
            tr = t["time"] / t["ref_time"]
            tl = t[t["ref_time"] > 1.0]
            fid_rows.append(dict(
                dataset=name, algorithm=C.label(k), key=k, pairs=len(p),
                r2_nodes=round(r_squared(p["nodes"], p["ref_nodes"]), 3),
                exact_nodes_pct=round(100 * (p["nodes"] == p["ref_nodes"]).mean(), 1),
                within_2x_pct=round(100 * (ratio_n <= 2).mean(), 1),
                median_time_ratio=round(float(tr.median()), 3),
                median_time_ratio_ref_over_1s=round(float((tl["time"] / tl["ref_time"]).median()), 3) if len(tl) else np.nan,
                reference="own" if k not in C.REF_BORROW else f"{C.REF_BORROW[k]} (borrowed: comparison, not fidelity)"))

        cls = pd.DataFrame([r for r in class_rows if r["dataset"] == name])
        fid = pd.DataFrame([r for r in fid_rows if r["dataset"] == name])
        print("\nA. validation classification (PASS = sizes identical; timed out/missing = a reference or run timed out)")
        print(cls[["algorithm", "PASS", "FAIL", "timed_out_or_missing", "pass_pct"]].to_string(index=False))
        print("\nB/C. fidelity on pairs both programs completed (ratio = unified time / reference time)")
        print(fid[["algorithm", "pairs", "r2_nodes", "exact_nodes_pct", "within_2x_pct",
                   "median_time_ratio", "median_time_ratio_ref_over_1s"]].to_string(index=False))

        # ---- D: figures
        ks = [k for k in frames if not complete_pairs(frames[k]).empty]
        if ks:
            fig, axes = C.grid(len(ks))
            for ax, k in zip(axes, ks):
                p = complete_pairs(frames[k])
                exact = (p["nodes"] == p["ref_nodes"])
                ax.scatter(p["nodes"], p["ref_nodes"], s=7, alpha=0.35, color=C.ALGOS[k].color, edgecolors="none")
                lo, hi = min(p["nodes"].min(), p["ref_nodes"].min()), max(p["nodes"].max(), p["ref_nodes"].max())
                ax.plot([lo, hi], [lo, hi], color=C.INK2, ls="--", lw=0.9)
                ax.set_xscale("log"); ax.set_yscale("log")
                ax.set_xlabel("unified nodes"); ax.set_ylabel("reference nodes")
                ax.set_title(f"{C.label(k)}{ref_note(k)}\nR²={r_squared(p['nodes'], p['ref_nodes']):.3f}, "
                             f"exact {exact.mean():.1%}, n={len(p)}", fontsize=9.5)
            fig.suptitle(f"{name}: search-tree fidelity (dashed line = perfect agreement)")
            C.outer_labels(fig, axes)
            fig.tight_layout(); C.save(fig, FOLDER, f"nodes_{ds}.png"); plt.close(fig)

            fig, axes = C.grid(len(ks))
            for ax, k in zip(axes, ks):
                p = complete_pairs(frames[k])
                ut, rt = p["time"].clip(lower=1e-3), p["ref_time"].clip(lower=1e-3)
                ax.scatter(ut, rt, s=7, alpha=0.35, color=C.ALGOS[k].color, edgecolors="none")
                lo, hi = min(ut.min(), rt.min()), max(ut.max(), rt.max())
                ax.plot([lo, hi], [lo, hi], color=C.INK2, ls="--", lw=0.9)
                ax.set_xscale("log"); ax.set_yscale("log")
                ax.set_xlabel("unified time (s)"); ax.set_ylabel("reference time (s)")
                ax.set_title(f"{C.label(k)}{ref_note(k)}\nmedian unified/ref = "
                             f"{(p['time'] / p['ref_time'].where(p['ref_time'] > 0)).median():.2f}, n={len(p)}", fontsize=9.5)
            fig.suptitle(f"{name}: runtime fidelity (dashed line = equal time; reference clocks are ms-resolution)")
            C.outer_labels(fig, axes)
            fig.tight_layout(); C.save(fig, FOLDER, f"runtime_{ds}.png"); plt.close(fig)

    out = C.out_dir(FOLDER)
    pd.DataFrame(class_rows).to_csv(out / "validation_classification.csv", index=False)
    pd.DataFrame(fid_rows).to_csv(out / "validation_fidelity.csv", index=False)
    print("\nTables saved to", out.relative_to(C.REPO))


if __name__ == "__main__":
    main()
