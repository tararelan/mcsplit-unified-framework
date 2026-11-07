"""Mutation tests for the direction-aware RRSplit / SymSplit proofs.

Each mutant makes RRSplit or SymSplit ignore edge direction in one place,
as the reference implementations do. VeriPB must reject every one on at
least one directed instance. The unmutated solver must pass them all,
otherwise a "caught" result means nothing.

    python3 test/mutants_directed.py -n 40

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
        for _ in range(r.randint(2, per)):
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


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--solver", default=os.path.join(HERE, "..", "bin", "mcsp"), help="unmutated solver")
    ap.add_argument("--veripb", default="veripb")
    ap.add_argument("-n", "--count", type=int, default=30, help="instances per family")
    ap.add_argument("-j", "--jobs", type=int, default=os.cpu_count())
    ap.add_argument("--timeout", type=int, default=30)
    args = ap.parse_args()

    workdir = tempfile.mkdtemp(prefix="mutants_directed_")
    instances = make_instances(args.count, workdir)

    jobs, errs = run_jobs(args.solver, args, workdir, RR + SYM, instances)
    bad = [e for e in errs if e]
    print(f"baseline (unmutated): {len(errs) - len(bad)}/{len(errs)} verified")
    if bad:
        print("  first failure:", bad[0])
        print("baseline must pass before the mutants mean anything")
        sys.exit(2)

    survived = 0
    for i, (name, file, old, new, times, algos) in enumerate(MUTANTS):
        binary, err = build(workdir, f"m{i}", file, old, new, times)
        if err:
            print(f"?? {name}: {err}")
            survived += 1
            continue
        jobs, errs = run_jobs(binary, args, workdir, algos, instances)
        caught = sum(e is not None for e in errs)
        print(f"{'ok' if caught else 'SURVIVED'} {name}: caught on {caught}/{len(jobs)}")
        survived += caught == 0
    shutil.rmtree(workdir, ignore_errors=True)
    sys.exit(1 if survived else 0)


if __name__ == "__main__":
    main()
