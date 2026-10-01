"""Checks the solver's VeriPB proofs on generated and stored instances.

For each algorithm and instance, runs the solver with -p, checks the proof
with VeriPB, and checks that the certified optimum matches the solver's
answer. With --cake, proofs on undirected unlabelled instances are also
elaborated and checked by CakePB's verified MCIS checker (cake_pb_mcis),
which builds the model from the LAD files itself.

    python3 test/check_proofs.py                      # everything, default sizes
    python3 test/check_proofs.py -n 3 -A symsplit     # quick check of one algorithm
    python3 test/check_proofs.py --cake ~/cakepb-dev/graph/cake_pb_mcis
"""
import argparse
import os
import random
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import gen  # noqa: E402

GENERAL = ["mcsplit", "rl", "ll", "ll_lsm", "ll_lum", "dal", "dal_rl", "dal_dal",
           "dsb", "dsb_always", "dsb_never"]
# these ignore direction and labels, so only run on undirected unlabelled instances
UNDIRECTED_ONLY = ["rrsplit", "rrsplit_nomax", "rrsplit_noveq", "rrsplit_nobound",
                   "symsplit", "symsplit_valonly", "symsplit_varonly"]

# family -> (solver flags, works with UNDIRECTED_ONLY, LAD so CakePB can check it)
FAMILIES = {
    "random": ([], True, True),
    "directed": (["-i"], False, False),
    "labelled": (["-x"], False, False),
    "binary": ([], True, False),
    "tree": ([], True, True),
    "star": ([], True, True),
    "stored": ([], True, True),
}


def make_instances(family, count, workdir):
    """Yields (name, pattern file, target file, format flag) for a family."""
    if family == "stored":
        inst = os.path.join(HERE, "instances")
        for f in sorted(os.listdir(inst)):
            if f.endswith("-pattern.lad"):
                yield f[:-len("-pattern.lad")], os.path.join(inst, f), os.path.join(inst, f.replace("pattern", "target")), "-l"
        return
    for i in range(count):
        r = random.Random(f"{family}-{i}")
        a, b = (os.path.join(workdir, f"{family}{i}-{s}") for s in ("pattern", "target"))
        if family in ("random", "directed", "labelled"):
            n1 = r.randint(5, 12)
            n2 = n1 + r.randint(0, 4)
            opts = {"directed": family == "directed", "labels": 3 if family == "labelled" else 0}
            gen.random_graph(n1, r.randint(1, 8) / 10, 2 * i, a + ".dimacs", **opts)
            gen.random_graph(n2, r.randint(1, 8) / 10, 2 * i + 1, b + ".dimacs", **opts)
            fmt = "-d"
        elif family == "binary":
            n1 = r.randint(8, 17)
            gen.binary_graph(n1, r.randint(2, 7) / 10, 2 * i, a + ".bin")
            gen.binary_graph(n1 + r.randint(0, 4), r.randint(2, 7) / 10, 2 * i + 1, b + ".bin")
            yield f"{family}{i}", a + ".bin", b + ".bin", ""
            continue
        elif family == "tree":
            n1 = r.randint(8, 15)
            gen.tree_graph(n1, r.randint(0, 2), 2 * i + 5000, a + ".dimacs")
            gen.tree_graph(n1 + r.randint(0, 6), r.randint(0, 3), 2 * i + 5001, b + ".dimacs")
            fmt = "-d"
        else:
            gen.star_graph(r.randint(3, 5), r.randint(2, 5), r.randint(0, 2), 2 * i + 9000, a + ".dimacs")
            gen.star_graph(r.randint(3, 6), r.randint(2, 5), r.randint(0, 2), 2 * i + 9001, b + ".dimacs")
            fmt = "-d"
        # the solver swaps the graphs if the first is bigger, but CakePB doesn't
        sizes = [int(open(x + ".dimacs").readline().split()[2]) for x in (a, b)]
        if sizes[0] > sizes[1]:
            os.replace(a + ".dimacs", a + ".tmp")
            os.replace(b + ".dimacs", a + ".dimacs")
            os.replace(a + ".tmp", b + ".dimacs")
        if FAMILIES[family][2]:
            gen.dimacs_to_lad(a + ".dimacs", a + ".lad")
            gen.dimacs_to_lad(b + ".dimacs", b + ".lad")
            yield f"{family}{i}", a + ".lad", b + ".lad", "-l"
        else:
            yield f"{family}{i}", a + ".dimacs", b + ".dimacs", fmt


def run(cmd, timeout):
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        return p.returncode, p.stdout + p.stderr
    except subprocess.TimeoutExpired:
        return None, ""


def check(args, workdir, algo, family, name, pattern, target, fmt):
    """Returns None if everything verified, or a description of the failure."""
    flags = FAMILIES[family][0]
    prefix = os.path.join(workdir, f"{algo}-{name}")
    cmd = [args.solver, "-q", "-A", algo, "-p", prefix] + ([fmt] if fmt else []) + flags + [pattern, target]
    code, out = run(cmd, args.timeout)
    if code is None:
        return "solver timed out"
    if code != 0:
        return f"solver failed: {out.strip()[-300:]}"
    size = int(out.split()[0])
    # the proof's objective counts unmapped vertices of the smaller graph
    with open(prefix + ".opb") as f:
        n_g = next(line for line in f if line.startswith("min:")).count("_null")
    expected = f"s VERIFIED BOUNDS {n_g - size} <= obj <= {n_g - size}"
    elab = prefix + ".elab.pbp"
    use_cake = args.cake and FAMILIES[family][2]
    code, out = run([args.veripb] + (["--elaborate", elab] if use_cake else []) + [prefix + ".opb", prefix + ".pbp"],
                    args.timeout * 10)
    verdict = next((line for line in out.splitlines() if line.startswith("s ")), out.strip()[-300:])
    if verdict != expected:
        return f"VeriPB: {verdict!r}, expected {expected!r}"
    if use_cake:
        code, out = run([args.cake, pattern, target, elab], args.timeout * 10)
        if f"s VERIFIED MAX CIS SIZE |CIS| = {size}" not in out:
            return f"CakePB: {out.strip()[-300:]!r}, expected |CIS| = {size}"
    if not args.keep:
        for ext in (".opb", ".pbp", ".elab.pbp"):
            if os.path.exists(prefix + ext):
                os.remove(prefix + ext)
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--solver", default=os.path.join(HERE, "..", "bin", "mcsp"))
    ap.add_argument("--veripb", default="veripb")
    ap.add_argument("--cake", help="path to cake_pb_mcis, to also check with CakePB")
    ap.add_argument("-A", "--algorithm", action="append", help="algorithms to check (default: all)")
    ap.add_argument("-F", "--family", action="append", choices=list(FAMILIES), help="instance families (default: all)")
    ap.add_argument("-n", "--count", type=int, default=8, help="instances per generated family")
    ap.add_argument("-j", "--jobs", type=int, default=os.cpu_count())
    ap.add_argument("--timeout", type=int, default=60, help="solver timeout (checkers get 10x)")
    ap.add_argument("--keep", help="keep instances and proofs in this directory")
    args = ap.parse_args()

    algos = args.algorithm or GENERAL + UNDIRECTED_ONLY
    families = args.family or list(FAMILIES)
    workdir = args.keep or tempfile.mkdtemp(prefix="check_proofs_")
    os.makedirs(workdir, exist_ok=True)

    jobs = []
    for family in families:
        instances = list(make_instances(family, args.count, workdir))
        for algo in algos:
            if algo in UNDIRECTED_ONLY and not FAMILIES[family][1]:
                continue
            for inst in instances:
                jobs.append((algo, family) + inst)

    with ThreadPoolExecutor(args.jobs) as pool:
        results = list(pool.map(lambda j: (j, check(args, workdir, *j)), jobs))

    failures = 0
    summary = {}
    for (algo, family, name, pattern, target, fmt), err in results:
        ok, total = summary.get((algo, family), (0, 0))
        summary[(algo, family)] = (ok + (err is None), total + 1)
        if err:
            failures += 1
            flags = " ".join(FAMILIES[family][0] + ([fmt] if fmt else []))
            print(f"FAIL {algo} {name}: {err}\n     {args.solver} -q -A {algo} {flags} {pattern} {target}")
    for (algo, family), (ok, total) in summary.items():
        print(f"{'ok  ' if ok == total else 'FAIL'} {algo:18} {family:9} {ok}/{total}")
    print(f"{len(results) - failures}/{len(results)} verified" + (" (with CakePB)" if args.cake else ""))
    if not args.keep:
        shutil.rmtree(workdir, ignore_errors=True)
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
