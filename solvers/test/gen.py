"""Random graph generators for testing proof logging.

    gen.py random N P SEED OUT [--labels K] [--directed]   DIMACS, G(n, p)
    gen.py binary N P SEED OUT [--directed]                MIVIA binary, random 16-bit labels
    gen.py tree N EXTRA SEED OUT                           random tree plus EXTRA edges (many leaves)
    gen.py star HUBS LEAVES EXTRA SEED OUT                 hubs with many leaves (large twin classes)
    gen.py lad IN OUT                                      convert an undirected DIMACS graph to LAD
"""
import argparse
import random
import struct


def write_dimacs(out, n, edges, labels=None):
    with open(out, "w") as f:
        f.write(f"p edge {n} {len(edges)}\n")
        for i, j in edges:
            f.write(f"e {i + 1} {j + 1}\n")
        for i, l in enumerate(labels or []):
            f.write(f"n {i + 1} {l}\n")


def random_graph(n, p, seed, out, labels=0, directed=False):
    r = random.Random(seed)
    edges = [(i, j) for i in range(n) for j in range(n)
             if i != j and (directed or j > i) and r.random() < p]
    write_dimacs(out, n, edges, [r.randrange(labels) for _ in range(n)] if labels else None)


def binary_graph(n, p, seed, out, directed=False):
    r = random.Random(seed)
    word = lambda x: struct.pack("<H", x)
    adj = [[] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j and (directed or j > i) and r.random() < p:
                adj[i].append((j, r.randrange(65536)))
    with open(out, "wb") as f:
        f.write(word(n))
        for _ in range(n):
            f.write(word(r.randrange(65536)))
        for i in range(n):
            f.write(word(len(adj[i])))
            for j, l in adj[i]:
                f.write(word(j) + word(l))


def tree_graph(n, extra, seed, out):
    r = random.Random(seed)
    edges = set()
    for i in range(1, n):
        j = r.randrange(i) if r.random() < 0.6 else r.randrange(max(1, i // 3))
        edges.add((j, i))
    extra = min(extra, n * (n - 1) // 2 - len(edges))
    while extra > 0:
        i, j = r.randrange(n), r.randrange(n)
        if i != j and (min(i, j), max(i, j)) not in edges:
            edges.add((min(i, j), max(i, j)))
            extra -= 1
    write_dimacs(out, n, sorted(edges))


def star_graph(hubs, per, extra, seed, out):
    r = random.Random(seed)
    edges = set()
    n = hubs
    for i in range(1, hubs):
        edges.add((r.randrange(i), i))
    for h in range(hubs):
        for _ in range(r.randint(1, per)):
            edges.add((h, n))
            n += 1
    extra = min(extra, hubs * (hubs - 1) // 2 - (hubs - 1))
    while extra > 0:
        i, j = r.randrange(hubs), r.randrange(hubs)
        if i != j and (min(i, j), max(i, j)) not in edges:
            edges.add((min(i, j), max(i, j)))
            extra -= 1
    write_dimacs(out, n, sorted(edges))


def dimacs_to_lad(inp, out):
    n, adj = 0, []
    for line in open(inp):
        parts = line.split()
        if parts and parts[0] == "p":
            n = int(parts[2])
            adj = [set() for _ in range(n)]
        elif parts and parts[0] == "e":
            i, j = int(parts[1]) - 1, int(parts[2]) - 1
            if i != j:
                adj[i].add(j)
                adj[j].add(i)
        elif parts and parts[0] == "n":
            raise ValueError("LAD can't hold vertex labels")
    with open(out, "w") as f:
        f.write(f"{n}\n")
        for i in range(n):
            f.write(" ".join(map(str, [len(adj[i])] + sorted(adj[i]))) + "\n")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="kind", required=True)
    a = sub.add_parser("random")
    for name, t in (("n", int), ("p", float), ("seed", int), ("out", str)):
        a.add_argument(name, type=t)
    a.add_argument("--labels", type=int, default=0)
    a.add_argument("--directed", action="store_true")
    a = sub.add_parser("binary")
    for name, t in (("n", int), ("p", float), ("seed", int), ("out", str)):
        a.add_argument(name, type=t)
    a.add_argument("--directed", action="store_true")
    a = sub.add_parser("tree")
    for name, t in (("n", int), ("extra", int), ("seed", int), ("out", str)):
        a.add_argument(name, type=t)
    a = sub.add_parser("star")
    for name, t in (("hubs", int), ("leaves", int), ("extra", int), ("seed", int), ("out", str)):
        a.add_argument(name, type=t)
    a = sub.add_parser("lad")
    a.add_argument("inp")
    a.add_argument("out")
    args = ap.parse_args()
    if args.kind == "random":
        random_graph(args.n, args.p, args.seed, args.out, args.labels, args.directed)
    elif args.kind == "binary":
        binary_graph(args.n, args.p, args.seed, args.out, args.directed)
    elif args.kind == "tree":
        tree_graph(args.n, args.extra, args.seed, args.out)
    elif args.kind == "star":
        star_graph(args.hubs, args.leaves, args.extra, args.seed, args.out)
    else:
        dimacs_to_lad(args.inp, args.out)
