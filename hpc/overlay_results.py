"""
overlay_results.py - fill the blank/absent pairs of the published result files with rerun rows from the cluster.

The cluster's results folders hold only the rerun pairs for most algorithms, so the merged cluster file is PARTIAL.
Copying it over hpc/publication/ would throw away the good rows. This script instead takes each rerun file and, for
each pair, only fills in a pair that is absent or blank (unified_nodes empty) in the published file. Pairs that already
have a valid row are left untouched. Blank rows in the rerun file are never copied. Rows are never removed, except a
blank row that is replaced by its valid rerun row (and any duplicate blank rows for that pair).

Usage (from the repo root):
    python hpc/overlay_results.py NAME [NAME ...] [--dry-run] [--incoming hpc/rerun_incoming]
where NAME is a file stem without "_all.csv", e.g.  down_lv_with_ref  mcsplit_lv_with_ref  rl_pair_lv  rl_rate_lv
For each NAME it reads  hpc/rerun_incoming/NAME_all.csv  and updates  hpc/publication/NAME_all.csv  (same name only),
after copying the original to hpc/publication/backup_before_overlay/NAME_all.<timestamp>.csv.
"""
import argparse
import csv
import shutil
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
PUB = REPO / "hpc" / "publication"
NODES = 5  # column index of unified_nodes


def parse(line):
    """(key, fields) for a real result row, else None (stray lines such as BIDOMAINS rows are kept verbatim)."""
    if "BIDOMAINS" in line:
        return None
    try:
        f = next(csv.reader([line.strip("\r\n")], skipinitialspace=True))
    except Exception:
        return None
    if len(f) < 6 or not f[0].strip() or not f[1].strip():
        return None
    return (f[0].strip(), f[1].strip()), f


def valid(f):
    try:
        float(f[NODES])
        return True
    except (ValueError, IndexError):
        return False


def read(path):
    with open(path, encoding="utf-8", newline="") as fh:
        lines = fh.read().splitlines(keepends=True)
    return lines[0], lines[1:]


def stats(rows):
    ks = {}
    for k, f in rows:
        ks[k] = ks.get(k, False) or valid(f)
    return len(ks), sum(ks.values()), sum(1 for v in ks.values() if not v)


def overlay(name, incoming_dir, dry):
    inc, tgt = incoming_dir / f"{name}_all.csv", PUB / f"{name}_all.csv"
    for p in (inc, tgt):
        if not p.exists():
            print(f"{name}: SKIPPED, missing {p}")
            return
    h_in, l_in = read(inc)
    h_tg, l_tg = read(tgt)
    if h_in.strip() != h_tg.strip():
        print(f"{name}: SKIPPED, header differs between {inc.name} and {tgt.name}")
        return

    new, skipped_blank = {}, 0
    for line in l_in:
        p = parse(line)
        if p is None:
            continue
        k, f = p
        if valid(f):
            new[k] = line if line.endswith("\n") else line + "\n"
        else:
            skipped_blank += 1

    parsed = [parse(l) for l in l_tg]
    rows_before = [p for p in parsed if p]
    valid_have = {k for k, f in rows_before if valid(f)}
    n0, v0, b0 = stats(rows_before)

    out, done, filled, kept = [], set(), 0, 0
    for line, p in zip(l_tg, parsed):
        if p is None:
            out.append(line)
            continue
        k, f = p
        if k in valid_have:
            out.append(line)
        elif k in new:
            if k not in done:
                out.append(new[k])
                done.add(k)
                filled += 1
        else:
            out.append(line)
    seen = {p[0] for p in parsed if p}
    added = [k for k in new if k not in seen]
    for k in added:
        out.append(new[k])
    kept = sum(1 for k in new if k in valid_have)

    rows_after = [p for p in (parse(l) for l in out) if p]
    n1, v1, b1 = stats(rows_after)
    assert n1 >= n0 and v1 >= v0, "refusing: pair or valid-row count would drop"
    print(f"{name}: pairs {n0} -> {n1}, valid {v0} -> {v1}, blank {b0} -> {b1} | "
          f"filled {filled}, added {len(added)}, rerun-valid-but-already-had-valid {kept}, rerun-blank-ignored {skipped_blank}")
    if dry:
        return
    bdir = PUB / "backup_before_overlay"
    bdir.mkdir(exist_ok=True)
    shutil.copy2(tgt, bdir / f"{name}_all.{time.strftime('%Y%m%d-%H%M%S')}.csv")
    with open(tgt, "w", encoding="utf-8", newline="") as fh:
        fh.write(h_tg)
        fh.writelines(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("names", nargs="+")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--incoming", default="hpc/rerun_incoming")
    a = ap.parse_args()
    inc = Path(a.incoming)
    inc = inc if inc.is_absolute() else REPO / inc
    for n in a.names:
        overlay(n.replace("_all.csv", ""), inc, a.dry_run)
    if a.dry_run:
        print("(dry run: nothing written)")


if __name__ == "__main__":
    main()
