"""Mutation tests for the direction-aware RRSplit / SymSplit proofs.

Each mutant makes RRSplit or SymSplit ignore edge direction in one place,
as the reference implementations do. VeriPB must reject every one on at
least one directed instance. The unmutated solver must pass them all,
otherwise a "caught" result means nothing.

    python3 test/mutants_directed.py -n 40
    python3 test/mutants_directed.py -n 40 --csv mutants_directed.csv   # also save the counts

Instances are directed random graphs plus directed "stars": hubs with leaf
groups whose arcs point randomly in or out, so twins exist and so do leaves
that look alike ignoring direction but are not twins.
"""
import argparse
import os
import random
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
import types
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "src")
sys.path.insert(0, HERE)
import check_proofs  # noqa: E402
import gen  # noqa: E402

check_proofs.FAMILIES["dstar"] = (["-i"], True, False)

RR = ["rrsplit", "rrsplit_nomax", "rrsplit_noveq", "rrsplit_nobound"]
SYM = ["symsplit", "symsplit_valonly", "symsplit_varonly"]

# name, file, old text, new text, how many times old occurs, algorithms
MUTANTS = [
    ("twin test compares adjacency, not direction", "graph.cpp",
     "g.adjmat[i][node_nneg] != g.adjmat[node_neg][node_nneg]",
     "(g.adjmat[i][node_nneg] != 0) != (g.adjmat[node_neg][node_nneg] != 0)", 1, RR + SYM),
    ("twin test drops the check on the edge between the twins", "graph.cpp",
     "if (equiv && g.adjmat[i][node_neg] != g.adjmat[node_neg][i])",
     "if (equiv && false)", 1, RR + SYM),
    ("RRSplit doesn't split bidomains by edge direction", "rrsplit.cpp",
     "if (rr_multiway && left_len && right_len)",
     "if (false && left_len && right_len)", 2, RR),
    ("SymSplit doesn't split bidomains by edge direction", "symsplit.cpp",
     "if (sym_multiway && left_len && right_len)",
     "if (false && left_len && right_len)", 2, SYM),
]


def directed_star(seed, out, hubs, per):
    r = random.Random(seed)
    edges, n = set(), hubs
    for i in range(1, hubs):
        a = r.randrange(i)
        edges.add((a, i) if r.random() < 0.5 else (i, a))
    for h in range(hubs):
        # a group of leaves with one shared orientation, so they look like
        # twins; sometimes one arc joins two of them, and swapping those two
        # would reverse it, so they are not twins. Plus a few loose leaves.
        out_arcs = r.random() < 0.5
        group = []
        for _ in range(r.randint(2, per)):
            edges.add((h, n) if out_arcs else (n, h))
            group.append(n)
            n += 1
        if len(group) >= 2 and r.random() < 0.6:
            a, b = r.sample(group, 2)
            edges.add((a, b))
        for _ in range(r.randint(0, 2)):
            edges.add((h, n) if r.random() < 0.5 else (n, h))
            n += 1
    gen.write_dimacs(out, n, sorted(edges))
    return n


def make_instances(count, workdir):
    out = []
    for i in range(count):
        r = random.Random(f"dstar-{i}")
        a, b = (os.path.join(workdir, f"dstar{i}-{s}.dimacs") for s in ("pattern", "target"))
        na = directed_star(2 * i + 100, a, r.randint(2, 3), r.randint(2, 4))
        nb = directed_star(2 * i + 101, b, r.randint(2, 4), r.randint(2, 5))
        if na > nb:
            a, b = b, a
        out.append(("dstar", f"dstar{i}", a, b, "-d"))
    out += [("directed",) + inst for inst in check_proofs.make_instances("directed", count, workdir)]
    return out


def build(workdir, name, file, old, new, times):
    src = os.path.join(workdir, name, "src")
    shutil.copytree(SRC, src)
    path = os.path.join(src, file)
    text = open(path, newline="").read()
    if text.count(old) != times:
        return None, f"stale: text found {text.count(old)} times in {file}, expected {times}"
    open(path, "w", newline="").write(text.replace(old, new))
    binary = os.path.join(workdir, name, "mcsp")
    flags = shlex.split(os.environ.get("CXXFLAGS", "-O2 -std=c++17"))
    sources = sorted(os.path.join(src, f) for f in os.listdir(src) if f.endswith(".cpp"))
    p = subprocess.run([os.environ.get("CXX", "c++")] + flags + ["-o", binary] + sources
                       + shlex.split(os.environ.get("LDFLAGS", "")), capture_output=True, text=True)
    if p.returncode != 0:
        return None, "build failed: " + p.stderr[-300:]
    return binary, None


def run_jobs(solver, args, workdir, algos, instances):
    opts = types.SimpleNamespace(solver=solver, veripb=args.veripb, cake=None, timeout=args.timeout, keep=None)
    jobs = [(a,) + inst for a in algos for inst in instances]
    with ThreadPoolExecutor(args.jobs) as pool:
        errs = list(pool.map(lambda j: check_proofs.check(opts, workdir, *j), jobs))
    return jobs, errs


def kind(err):
    """Why a job failed: what the proof checker or the solver said."""
    if err is None:
        return "verified"
    for prefix, name in (("VeriPB", "VeriPB rejected the proof"), ("solver timed out", "solver timed out"),
                         ("solver failed", "solver crashed"), ("CakePB", "CakePB rejected the proof")):
        if err.startswith(prefix):
            return name
    return "other"


def tally(jobs, errs):
    """(algorithm, family) -> [jobs, failed]. jobs are (algo, family, name, pattern, target, fmt)."""
    t = {}
    for j, e in zip(jobs, errs):
        c = t.setdefault((j[0], j[1]), [0, 0])
        c[0] += 1
        c[1] += e is not None
    return t


def show(t, algos, families, label, ok_word):
    """One line per algorithm x family. label says what the number counts (verified / caught)."""
    for a in algos:
        for f in families:
            if (a, f) not in t:
                continue
            total, bad = t[(a, f)]
            n = total - bad if ok_word == "verified" else bad
            print(f"  {a:18} {f:9} {label} {n}/{total}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--solver", default=os.path.join(HERE, "..", "bin", "mcsp"), help="unmutated solver")
    ap.add_argument("--veripb", default="veripb")
    ap.add_argument("-n", "--count", type=int, default=30, help="instances per family")
    ap.add_argument("-j", "--jobs", type=int, default=os.cpu_count())
    ap.add_argument("--timeout", type=int, default=30)
    ap.add_argument("--csv", help="also write every (stage, algorithm, family) count to this CSV file")
    ap.add_argument("--no-controls", action="store_true",
                    help="only run the algorithms each mutant is meant to break (faster); by default the other "
                         "algorithms are run too, to show the mutant does not leak into them")
    args = ap.parse_args()

    t0 = time.time()
    workdir = tempfile.mkdtemp(prefix="mutants_directed_")
    instances = make_instances(args.count, workdir)
    families = list(dict.fromkeys(i[0] for i in instances))
    ALL = RR + SYM
    print("instances: " + ", ".join(f"{f} {sum(1 for i in instances if i[0] == f)}" for f in families)
          + f"   algorithms: {len(ALL)} ({len(RR)} RRSplit, {len(SYM)} SymSplit)")
    rows = []                                              # for --csv

    jobs, errs = run_jobs(args.solver, args, workdir, ALL, instances)
    bad = [e for e in errs if e]
    print(f"\nbaseline (unmutated): {len(errs) - len(bad)}/{len(errs)} verified   [{time.time() - t0:.0f}s]")
    base = tally(jobs, errs)
    show(base, ALL, families, "verified", "verified")
    for (a, f), (n, b_) in base.items():
        rows.append(dict(stage="baseline", mutant="", algorithm=a, family=f, role="baseline", jobs=n, failed=b_))
    if bad:
        print("  first failure:", bad[0])
        print("baseline must pass before the mutants mean anything")
        sys.exit(2)

    survived = 0
    for i, (name, file, old, new, times, algos) in enumerate(MUTANTS):
        print(f"\nmutant {i + 1}: {name}   [{file}]")
        binary, err = build(workdir, f"m{i}", file, old, new, times)
        if err:
            print(f"?? {name}: {err}")
            survived += 1
            continue
        controls = [] if args.no_controls else [a for a in ALL if a not in algos]
        jobs, errs = run_jobs(binary, args, workdir, algos + controls, instances)
        hit = [(j, e) for j, e in zip(jobs, errs) if j[0] in algos]
        caught = sum(e is not None for _, e in hit)
        t = tally(jobs, errs)
        print(f" affected algorithms (the mutant changes their code; VeriPB should reject):")
        show(t, algos, families, "caught", "caught")
        kinds = {}
        for _, e in hit:
            if e is not None:
                kinds[kind(e)] = kinds.get(kind(e), 0) + 1
        if kinds:
            print("  how they were caught: " + ", ".join(f"{k} {v}" for k, v in sorted(kinds.items())))
        if controls:
            print(" other algorithms (code the mutant does not touch; should still verify):")
            show(t, controls, families, "verified", "verified")
            leaked = sum(t[k][1] for k in t if k[0] in controls)
            if leaked:
                print(f"  !! {leaked} control job(s) failed: the mutant leaks into algorithms it should not affect")
        for (a, f), (n, b_) in t.items():
            rows.append(dict(stage="mutant", mutant=name, algorithm=a, family=f,
                             role="affected" if a in algos else "control", jobs=n, failed=b_))
        print(f"{'ok' if caught else 'SURVIVED'} {name}: caught on {caught}/{len(hit)}")
        survived += caught == 0

    print(f"\n{len(MUTANTS) - survived}/{len(MUTANTS)} mutants caught, {time.time() - t0:.0f}s in total")
    if args.csv:
        pd_cols = ["stage", "mutant", "algorithm", "family", "role", "jobs", "failed"]
        with open(args.csv, "w") as f:
            f.write(",".join(pd_cols) + "\n")
            for r in rows:
                f.write(",".join('"%s"' % str(r[c]) if c == "mutant" else str(r[c]) for c in pd_cols) + "\n")
        print("wrote", args.csv)
    shutil.rmtree(workdir, ignore_errors=True)
    sys.exit(1 if survived else 0)


if __name__ == "__main__":
    main()
