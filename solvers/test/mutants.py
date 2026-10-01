"""Mutation tests for proof logging: each mutant makes the solver prune
unsoundly, or makes the proof logger justify something wrongly, and must
be caught (a proof that fails to verify, or a wrong answer) on at least one
instance. A mutant that survives means the proofs aren't checking what
they should be.

    python3 test/mutants.py
    CXXFLAGS="-O2 -std=c++17 -I/opt/homebrew/include" LDFLAGS="-L/opt/homebrew/lib -largp" python3 test/mutants.py

The mutations are exact text replacements, so after changing the code they
can go stale; a mutant whose text isn't found is reported as such.
"""
import argparse
import os
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

# name, file, old text, new text, algorithms to run it with
MUTANTS = [
    ("mcsplit prunes one too early", "mcsplit.cpp",
     "if (bound <= (int)incumbent.size() || bound < (int)goal) {",
     "if (bound <= (int)incumbent.size() + 1 || bound < (int)goal) {",
     ["mcsplit"]),
    ("RRSplit skips the value after the previous twin's", "rrsplit.cpp",
     "w_lower = a.w;", "w_lower = a.w + 1;", ["rrsplit", "rrsplit_nomax"]),
    ("SymSplit skips the value after the previous twin's", "symsplit.cpp",
     "w_lower = a.w;", "w_lower = a.w + 1;", ["symsplit", "symsplit_varonly"]),
    ("SymSplit value symmetry ignores twin classes", "symsplit.cpp",
     "if (h_eqn_classes[arr[start_idx + i]] == w_class && arr[start_idx + i] < w)",
     "if (arr[start_idx + i] < w)", ["symsplit", "symsplit_valonly"]),
    ("SymSplit maximality doesn't rule out bottom", "symsplit.cpp",
     "proof->maximal_match_with_twins(v, w, true);", "proof->maximal_match_with_twins(v, w, false);",
     ["symsplit", "symsplit_varonly"]),
    ("LUM doesn't rule out the leaf being unmapped", "proof.cpp",
     'std::vector<long long> ruled_out{rule_out_value(v, w, -1, "red")};',
     "std::vector<long long> ruled_out;", ["ll", "dal", "symsplit_valonly"]),
    ("lex order reversed", "proof.cpp",
     "lex_weight[g_order[i]] = dec_mul(lex_weight[g_order[i + 1]], B);",
     "lex_weight[g_order[g_n - 1 - i]] = dec_mul(lex_weight[g_order[g_n - 2 - i]], B);",
     ["rrsplit", "symsplit"]),
    ("value symmetry swaps with the wrong vertex", "proof.cpp",
     "b = h_order[b];", "b = h_order[(b + 1) % h_n];", ["symsplit", "symsplit_valonly"]),
    ("maximality rewrite unmaps the holder instead", "proof.cpp",
     '<< match_var(a, w) << " -> 0 " << value(a, z) << " -> 1 ;\\n";',
     '<< match_var(a, w) << " -> 0 " << value(a, -1) << " -> 1 ;\\n";',
     ["rrsplit", "symsplit", "ll"]),
    ("McSplit bound ignores a bidomain's excess", "proof.cpp",
     "if (i == bd_idx || bd.left_len <= bd.right_len)",
     "if (i == bd_idx || bd.left_len <= bd.right_len + 1)", ["rrsplit_nomax"]),
]


def build(workdir, name, file, old, new):
    src = os.path.join(workdir, name, "src")
    shutil.copytree(SRC, src)
    path = os.path.join(src, file)
    text = open(path).read()
    if text.count(old) != 1:
        return None, f"stale: text found {text.count(old)} times in {file}"
    open(path, "w").write(text.replace(old, new))
    binary = os.path.join(workdir, name, "mcsp")
    cxx = os.environ.get("CXX", "c++")
    flags = shlex.split(os.environ.get("CXXFLAGS", "-O2 -std=c++17"))
    ldflags = shlex.split(os.environ.get("LDFLAGS", ""))
    sources = sorted(os.path.join(src, f) for f in os.listdir(src) if f.endswith(".cpp"))
    p = subprocess.run([cxx] + flags + ["-o", binary] + sources + ldflags, capture_output=True, text=True)
    if p.returncode != 0:
        return None, "build failed: " + p.stderr[-300:]
    return binary, None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--veripb", default="veripb")
    ap.add_argument("-n", "--count", type=int, default=8, help="instances per generated family")
    ap.add_argument("-j", "--jobs", type=int, default=os.cpu_count())
    ap.add_argument("--timeout", type=int, default=30)
    args = ap.parse_args()

    workdir = tempfile.mkdtemp(prefix="mutants_")
    families = ["random", "tree", "star", "stored"]
    instances = [(f,) + inst for f in families for inst in check_proofs.make_instances(f, args.count, workdir)]
    survived = 0
    for i, (name, file, old, new, algos) in enumerate(MUTANTS):
        binary, err = build(workdir, f"m{i}", file, old, new)
        if err:
            print(f"?? {name}: {err}")
            survived += 1
            continue
        opts = types.SimpleNamespace(solver=binary, veripb=args.veripb, cake=None, timeout=args.timeout, keep=None)
        jobs = [(algo,) + inst for algo in algos for inst in instances]
        with ThreadPoolExecutor(args.jobs) as pool:
            caught = sum(r is not None for r in pool.map(lambda j: check_proofs.check(opts, workdir, *j), jobs))
        print(f"{'ok' if caught else 'SURVIVED'} {name}: caught on {caught}/{len(jobs)}")
        survived += caught == 0
    shutil.rmtree(workdir, ignore_errors=True)
    sys.exit(1 if survived else 0)


if __name__ == "__main__":
    main()
