"""
06_contested.py - why does an algorithm lose instances that McSplit solves?     (dissertation Section 4.2)

Contested instances: medium instances McSplit solves within the timeout but algorithm X does not.
For each one we ask: had X searched McSplit's tree at X's own per-node cost, would it have finished
within the timeout?
  no  -> "overhead only"  (X's tree is not the problem; its per-node work is)
  yes -> "worse search"   (X must have explored a larger tree than McSplit's)

X's per-node cost (seconds/node) is taken two ways:
  per-instance : X's own time/nodes on that instance. A timed-out run still reports the nodes it
                 reached, so this is the rate X was actually running at.            <- the one to report
  global       : median time/nodes over the medium instances X solved (the Table 9 figure); a
                 cruder check that assumes every instance costs the median.

Not meaningful for the top-down algorithms (down, rl_down): their node count on a failed run is not
comparable with McSplit's bottom-up tree, so read their rows with care.

Usage:  python 06_contested.py [--datasets mivia lv bi] [--include-unsound]
Output: analysis/outputs/performance/contested.csv
"""
import argparse

import pandas as pd

import config as C

FOLDER = "performance"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datasets", nargs="+", default=list(C.DATASETS))
    ap.add_argument("--include-unsound", action="store_true")
    args = ap.parse_args()

    rows = []
    for ds in args.datasets:
        algos = C.available(ds, include_unsound=args.include_unsound)
        if "mcsplit" not in algos:
            print(f"{ds}: no McSplit results, skipped"); continue
        cls = C.classify(ds)
        med = cls[cls["cls"] == "medium"][["a", "b"]]

        def get(k):
            return med.merge(C.load(k, ds), on=["a", "b"]).set_index(["a", "b"])

        M = get("mcsplit")
        mdone = M.done.astype(bool)
        print(f"\n{C.DATASETS[ds]}: {len(M)} medium pairs, McSplit solves {int(mdone.sum())}, timeout {C.TIMEOUT_S}s")
        print(f"  {'algorithm':16s} {'contested':>9s} | {'per-inst overhead':>20s} {'worse':>6s} | "
              f"{'global overhead':>18s} {'worse':>6s} | {'solves, McSplit not':>19s}")
        for k in algos:
            if k == "mcsplit":
                continue
            X = get(k).reindex(M.index)
            xdone = X.done.fillna(False).astype(bool)
            cont = mdone & ~xdone
            n = int(cont.sum())
            gain = int((xdone & ~mdone).sum())
            if n == 0:
                rows.append(dict(dataset=ds, algo=k, contested=0, overhead_inst=0, worse_inst=0,
                                 overhead_global=0, worse_global=0, solves_mcsplit_not=gain))
                print(f"  {C.label(k):16s} {0:9d}"); continue
            glob = (X.time[xdone] / X.nodes[xdone]).median()
            o_i = int((M.nodes[cont] * (X.time[cont] / X.nodes[cont]) > C.TIMEOUT_S).sum())
            o_g = int((M.nodes[cont] * glob > C.TIMEOUT_S).sum())
            rows.append(dict(dataset=ds, algo=k, contested=n, overhead_inst=o_i, worse_inst=n - o_i,
                             overhead_global=o_g, worse_global=n - o_g, solves_mcsplit_not=gain))
            print(f"  {C.label(k):16s} {n:9d} | {o_i:9d} ({100 * o_i / n:4.0f}%) {n - o_i:6d} | "
                  f"{o_g:9d} ({100 * o_g / n:4.0f}%) {n - o_g:6d} | {gain:19d}")

    if rows:
        out = C.out_dir(FOLDER) / "contested.csv"
        pd.DataFrame(rows).to_csv(out, index=False)
        print("\n  saved", out.relative_to(C.REPO))


if __name__ == "__main__":
    main()
