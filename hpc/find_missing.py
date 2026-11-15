"""
find_missing.py - list, per algorithm and dataset, the instance pairs that have no usable result.

A pair is MISSING when its row in hpc/publication/<algo>_<ds>[_with_ref]_all.csv is absent, or present but blank
(unified_nodes empty: the run was killed, e.g. out of memory, and wrote an empty row). Stray lines in a file
(e.g. "BIDOMAINS:..." rows) are ignored; only rows for real instance pairs count.

Usage (from the repo root):   python hpc/find_missing.py [--datasets mivia lv bi] [--algos mcsplit dsb ...]
Output: hpc/missing/<algo>_<dataset>.txt, one "pathA pathB" line per missing pair (same format as the
        hpc/instances_*.txt lists the submit scripts build), plus a summary table. Lists for pairs that are
        no longer missing are removed from hpc/missing/ so the folder always reflects the current files.
Then rerun with:              bash hpc/resubmit_missing.sh <algo> <dataset> [--with-ref] [--swap] [--mem=40G]
"""
import argparse
import re
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
PUB = REPO / "hpc" / "publication"
OUT = REPO / "hpc" / "missing"
ALGOS = ["mcsplit", "down", "rl", "rl_pair", "rl_rate", "rl_down", "ll", "dal", "dsb", "rrsplit", "symsplit"]


def expected_pairs(ds):
    """The pair list the submit scripts build for each dataset, as (pathA, pathB) strings relative to the repo."""
    if ds == "bi":
        files = sorted(p.as_posix() for p in (REPO / "instances/SIP/biochemicalReactions").glob("*.txt"))
        files = [str(Path(f).relative_to(REPO).as_posix()) for f in files]
        return [(files[i], files[j]) for i in range(len(files)) for j in range(i + 1, len(files))]
    if ds == "lv":
        files = [p for p in (REPO / "instances/SIP/LV").glob("g*") if not p.name.startswith("._")]
        files = sorted(files, key=lambda p: int(re.sub(r"\D", "", p.name) or 0))
        files = [p.relative_to(REPO).as_posix() for p in files]
        return [(files[i], files[j]) for i in range(len(files)) for j in range(i + 1, len(files))]
    if ds == "mivia":
        out = []
        for a in sorted((REPO / "instances/MIVIA").rglob("*.A0[0-9]")):
            b = Path(str(a).replace(".A", ".B", 1))
            if b.exists():
                out.append((a.relative_to(REPO).as_posix(), b.relative_to(REPO).as_posix()))
        return out
    raise SystemExit(f"unknown dataset {ds}")


def published_file(algo, ds):
    for name in (f"{algo}_{ds}_with_ref_all.csv", f"{algo}_{ds}_all.csv"):
        if (PUB / name).exists():
            return PUB / name
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datasets", nargs="+", default=["mivia", "lv", "bi"])
    ap.add_argument("--algos", nargs="+", default=ALGOS)
    args = ap.parse_args()
    OUT.mkdir(exist_ok=True)
    rows = []
    for ds in args.datasets:
        pairs = expected_pairs(ds)
        keyed = {(Path(a).name, Path(b).name): (a, b) for a, b in pairs}
        for algo in args.algos:
            f = published_file(algo, ds)
            target = OUT / f"{algo}_{ds}.txt"
            if f is None:
                rows.append((ds, algo, len(pairs), "no file", "-")); continue
            d = pd.read_csv(f, skipinitialspace=True, dtype=str)
            d.columns = d.columns.str.strip()
            d["k"] = list(zip(d["instance_a"].str.strip(), d["instance_b"].str.strip()))
            d = d[d["k"].isin(keyed)]
            ok = d[pd.to_numeric(d["unified_nodes"], errors="coerce").notna()]["k"]
            have, seen = set(ok), set(d["k"])
            absent = [k for k in keyed if k not in seen]
            blank = [k for k in keyed if k in seen and k not in have]
            missing = [keyed[k] for k in keyed if k not in have]
            if missing:
                target.write_text("".join(f"{a} {b}\n" for a, b in missing), encoding="utf-8", newline="\n")
            elif target.exists():
                target.unlink()
            rows.append((ds, algo, len(pairs), len(missing), f"{len(blank)} blank + {len(absent)} absent"))
    print(f"{'dataset':8s} {'algorithm':9s} {'pairs':>6s} {'missing':>8s}  detail")
    for ds, algo, n, m, det in rows:
        print(f"{ds:8s} {algo:9s} {n:6d} {str(m):>8s}  {det}")
    print(f"\nLists written to {OUT.relative_to(REPO)}/")


if __name__ == "__main__":
    main()
