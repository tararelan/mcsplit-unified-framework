"""
config.py - shared settings, colours and data loaders for every script in analysis/.

Everything an analysis script needs to know lives here, so the scripts themselves stay short:
  * ALGOS      every algorithm: label, family, colour, line style, marker
  * DATASETS   mivia / lv / bi: where each dataset's result CSVs live
  * load()     one loader that reads any of the result-file layouts into ONE standard table
  * helpers    medium-instance classification, plotting style, output folders

Standard table columns (one row per instance pair):
  a, b, algo                         instance names, algorithm key
  size, edges, nodes, time           unified solver: solution size / edges / search nodes / seconds
  aborted                            1 = hit the timeout
  nodes_to_best, time_to_best, cut, bound_pruned, sym_pruned
  bd_mean, bd_max                    bidomain count per search node (mean, max), if the run recorded it
  ref_size, ref_nodes, ref_time      reference binary (blank if it timed out / has no reference)
  done, ref_done                     solved within the timeout (unified / reference)
  match                              PASS (both done, same size) | FAIL (both done, sizes differ) | MISSING
"""
import io
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------- paths
TIMEOUT_S = 1000      # per-instance timeout used for every run
EASY_S = 10           # "easy" = every algorithm solves it within this many seconds


def _find_repo():
    here = Path(__file__).resolve()
    for p in [here.parent, *here.parents]:
        if (p / "hpc").is_dir() and (p / "solvers").is_dir():
            return p
    raise SystemExit("config.py: could not find the repo root (a folder containing hpc/ and solvers/)")


REPO = _find_repo()
OUT = Path(__file__).resolve().parent / "outputs"


def out_dir(name):
    d = OUT / name
    d.mkdir(parents=True, exist_ok=True)
    return d


# --------------------------------------------------------------------------- algorithms
@dataclass(frozen=True)
class Algo:
    key: str        # the -A flag / file prefix
    label: str      # name used in figures and tables
    family: str     # baseline | rl | ll | dsb | rr | sym
    color: str      # hue = family
    ls: str         # line style  -> tells variants of one family apart
    marker: str     # marker      -> same


# Colour = family (hue), line style + marker = variant within the family. Six hues + a neutral grey
# keep every pair distinguishable (checked adjacent-pair for colour-vision deficiency, light surface).
_GREY, _BLUE, _ORANGE, _AQUA, _YELLOW, _MAGENTA, _VIOLET = (
    "#52514e", "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#6250d6")

ALGOS = {a.key: a for a in [
    Algo("mcsplit", "McSplit",            "baseline", "#3d3c39", "-",  "o"),
    Algo("down",    "McSplit-down",       "baseline", "#a67c4e", "--", "o"),
    Algo("rl",      "McSplit+RL",         "rl",       _BLUE,    "-",  "s"),
    Algo("rl_pair", "McSplit+RL Pair",    "rl",       "#0b2f7a", "--", "s"),
    Algo("rl_rate", "McSplit+RL Rate",    "rl",       "#5bbcf2", ":",  "s"),
    Algo("rl_down", "McSplit+RL-down",    "rl",       "#12a3c7", "-.", "s"),
    Algo("ll",      "McSplit+LL",         "ll",       _AQUA,    "-",  "^"),
    Algo("dal",     "McSplit+DAL",        "dal",      _ORANGE,  "-",  "*"),   # DAL: comment out to drop it everywhere
    Algo("dsb",     "McSplit+DSB",        "dsb",      _YELLOW,  "-",  "D"),
    Algo("rrsplit", "RRSplit",            "rr",       _VIOLET,  "-",  "p"),
    Algo("symsplit", "SymSplit",          "sym",      _MAGENTA, "-",  "h"),
]}
ALL_ALGOS = list(ALGOS)

# Marker shape = algorithm family, the same in every figure: circle McSplit / down, square the RL family,
# triangle LL, star DAL, diamond DSB, pentagon RRSplit, hexagon SymSplit. (Colour and line style tell variants apart.)
FAMILY_MARKER = {"baseline": "o", "rl": "s", "ll": "^", "dal": "*", "dsb": "D", "rr": "p", "sym": "h"}


def marker_of(key):
    """Marker of an algorithm, or of an ablation variant (the shape of its parent's family)."""
    if key in ALGOS:
        return ALGOS[key].marker
    parent = next(p for k, _, p in VARIANTS if k == key)
    return ALGOS[parent].marker

# Ablation variants without their own ALGOS entry: same hue family as the parent, clearly different shade.
VARIANT_COLORS = {
    "ll_lsm": "#0d6e4b",  "ll_lum": "#7fe0b8",
    "dal_rl": "#a8370b", "dal_dal": "#f6b08c",
    "dsb_never": "#a86a00", "dsb_always": "#f7d36b",
    "rrsplit_nomax": "#3a2c9e", "rrsplit_nobound": "#a79bf2", "rrsplit_noveq": "#8a3fc0",
    "symsplit_varonly": "#b8316c", "symsplit_valonly": "#f4a9c8",
}


def color_of(key):
    return ALGOS[key].color if key in ALGOS else VARIANT_COLORS[key]


def label(key):
    return ALGOS[key].label if key in ALGOS else key


# Algorithms with no reference binary of their own are compared against another algorithm's reference.
# (That is a comparison with RL, not a fidelity check - the validation outputs say so.)
REF_BORROW = {"rl_pair": "rl", "rl_rate": "rl", "rl_down": "rl"}

# The set of algorithms that defines easy / medium / hard: every algorithm evaluated in the paper
# (RRSplit and SymSplit are left out on BI because their reference is unsound there). Hard = none of them
# solves it, so adding a stronger algorithm shrinks the hard class. The RL variants (Pair/Rate/Down) have no
# reference of their own; their correctness is checked separately (proof logging).
_ALL = ["mcsplit", "down", "rl", "rl_pair", "rl_rate", "rl_down", "ll", "dsb"]
CORE = {
    "mivia": _ALL + ["rrsplit", "symsplit"],
    "lv":    _ALL + ["rrsplit", "symsplit"],
    "bi":    list(_ALL),
}

# Algorithms whose reference is known to be unsound on a dataset: shown nowhere unless include_unsound=True.
UNSOUND = {("rrsplit", "bi"), ("symsplit", "bi")}

# Ablation variants (Table 11 of the dissertation): key, label, parent algorithm
VARIANTS = [
    ("ll_lsm",           "LSM only",              "ll"),
    ("ll_lum",           "LUM only",              "ll"),
    ("dal_rl",           "RL policy only",        "dal"),   # DAL: comment out to drop it from the ablation
    ("dal_dal",          "DAL policy only",       "dal"),   # DAL
    ("dsb_never",        "No DSB",                "dsb"),
    ("dsb_always",       "Always DSB",            "dsb"),
    ("rrsplit_nomax",    "No maximality",         "rrsplit"),
    ("rrsplit_nobound",  "No bound",              "rrsplit"),
    ("rrsplit_noveq",    "No vertex equivalence", "rrsplit"),
    ("symsplit_varonly", "Variable symmetry only", "symsplit"),
    ("symsplit_valonly", "Value symmetry only",   "symsplit"),
    # these are full algorithms with their own results (not hpc/ablation): RL's variants and McSplit's down
    ("rl_pair",          "Pair reward",           "rl"),
    ("rl_rate",          "Rate reward",           "rl"),
    ("rl_down",          "Down",                  "rl"),
    ("ll",               "LL",                    "rl"),   # LL builds on RL (also the parent of its own LSM/LUM group)
    ("down",             "Down",                  "mcsplit"),
]

# --------------------------------------------------------------------------- datasets
DATASETS = {"mivia": "MIVIA", "lv": "LV", "bi": "BI"}

# Where each dataset's per-algorithm result CSVs are, in order of preference (first file that exists wins).
# kind tells the loader what layout to expect:
#   pub  hpc/publication  21-column "with_ref" files (current pipeline, has bidomain columns)
#   val  hpc/validation   original validation files (unified + reference in one row)
#   xds  hpc/cross-dataset-evaluation  unified columns; reference may be in a separate *_reference_all.csv
SOURCES = {
    "bi": [("hpc/publication/results", "{a}_bi_with_ref_all.csv", "pub"),
           ("hpc/publication/old", "{a}_bi_with_ref_all.csv", "pub"),
           ("hpc/cross-dataset-evaluation/results", "{a}_bi_all.csv", "xds")],
    "mivia": [("hpc/validation/results", "{a}_all.csv", "val"),
              ("hpc/publication/results", "{a}_mivia_with_ref_all.csv", "pub"),
              ("hpc/publication/old", "{a}_mivia_with_ref_all.csv", "pub")],
    "lv": [("hpc/publication/old", "{a}_lv_with_ref_all.csv", "pub"),
           ("hpc/cross-dataset-evaluation/results", "{a}_lv_all.csv", "xds"),
           ("hpc/publication/results", "{a}_lv_with_ref_all.csv", "pub")],
}
# Ablation runs (unified solver only). Tried in order; the old MIVIA run has no dataset in its file names.
ABLATION_SOURCES = ["hpc/ablation/results/{a}_{ds}_all.csv",
                    "hpc/ablation/results/{ds}/{a}_all.csv",
                    "hpc/ablation/results/{a}_all.csv"]       # MIVIA only (used for ds == "mivia")

# --------------------------------------------------------------------------- loading
_RENAME = {
    "instance_a": "a", "instance_b": "b",
    "unified_size": "size", "unified_edges": "edges", "unified_nodes": "nodes", "unified_time": "time",
    "unified_aborted": "aborted", "unified_nodes_to_best": "nodes_to_best",
    "unified_time_to_best": "time_to_best", "unified_cut_branches": "cut",
    "unified_bound_pruned": "bound_pruned", "unified_sym_pruned": "sym_pruned",
    "unified_bidomain_mean": "bd_mean", "unified_bidomain_max": "bd_max",
    "ref_bidomain_mean": "ref_bd_mean", "ref_bidomain_max": "ref_bd_max",
}
_NUMERIC = ["size", "edges", "nodes", "time", "aborted", "nodes_to_best", "time_to_best", "cut",
            "bound_pruned", "sym_pruned", "bd_mean", "bd_max", "ref_size", "ref_nodes", "ref_time",
            "ref_bd_mean", "ref_bd_max"]
_cache = {}


def _read_csv(path):
    """Read a result CSV. Older DSB/SymSplit files have rows split in two by a stray 'BIDOMAINS:' line
    (the reference printed an extra line); those are stitched back together here, in memory."""
    fixed = []
    for line in Path(path).read_text().splitlines():
        if line.startswith("BIDOMAINS:") and fixed:
            fixed[-1] += "," + ",".join(line.split(",")[1:])
        else:
            fixed.append(line)
    df = pd.read_csv(io.StringIO("\n".join(fixed)), skipinitialspace=True)
    df.columns = [c.strip() for c in df.columns]
    return df


def _add_reference_file(df, path, algo, ds):
    """lv files keep the reference in a separate headerless file: a, b, algo, size, nodes, time."""
    rpath = Path(path).with_name(f"{algo}_{ds}_reference_all.csv")
    if not rpath.exists():
        return df
    r = pd.read_csv(rpath, header=None, skipinitialspace=True,
                    names=["instance_a", "instance_b", "algo", "r_size", "r_nodes", "r_time"])
    df = df.merge(r[["instance_a", "instance_b", "r_size", "r_nodes", "r_time"]],
                  on=["instance_a", "instance_b"], how="left")
    for new, old in [("ref_size", "r_size"), ("ref_nodes", "r_nodes"), ("ref_time", "r_time")]:
        if new not in df.columns:
            df[new] = np.nan
        df[new] = pd.to_numeric(df[new], errors="coerce").fillna(pd.to_numeric(df[old], errors="coerce"))
    return df.drop(columns=["r_size", "r_nodes", "r_time"])


def _match(df):
    """Dissertation rule (Section 4.1): PASS = the reference finished and both sizes are identical;
    FAIL = both finished (within the timeout) but the sizes differ; MISSING = anything else
    (the reference timed out, or the re-implementation timed out with a different incumbent)."""
    m = pd.Series("MISSING", index=df.index)
    m[df["ref_done"] & (df["size"] == df["ref_size"])] = "PASS"
    m[df["ref_done"] & df["done"] & (df["size"] != df["ref_size"])] = "FAIL"
    return m


def _standardise(df, algo):
    df = df.rename(columns=_RENAME)
    for c in _NUMERIC:
        df[c] = pd.to_numeric(df[c], errors="coerce") if c in df.columns else np.nan
    df["algo"] = algo
    df["a"], df["b"] = df["a"].astype(str), df["b"].astype(str)
    # a reference that timed out leaves size/nodes blank (its time is written as 0)
    df["ref_done"] = df["ref_size"].notna() & df["ref_nodes"].notna() & (df["ref_nodes"] > 0)
    for c in ("ref_size", "ref_nodes", "ref_time"):
        df.loc[~df["ref_done"], c] = np.nan
    df["done"] = (df["aborted"] == 0) if df["aborted"].notna().any() else (df["time"] < TIMEOUT_S)
    df["match"] = _match(df)
    keep = ["a", "b", "algo"] + _NUMERIC + ["done", "ref_done", "match"]
    return df[keep].drop_duplicates(subset=["a", "b"], keep="first").reset_index(drop=True)


def load(algo, ds, need_bidomain=False):
    """Standard table for one algorithm on one dataset, or None if there is no result file.
    need_bidomain=True skips files that carry no bidomain statistics."""
    key = (algo, ds, need_bidomain)
    if key in _cache:
        return _cache[key]
    result = None
    for folder, tmpl, kind in SOURCES[ds]:
        path = REPO / folder / tmpl.format(a=algo)
        if not path.exists():
            continue
        df = _read_csv(path)
        if kind == "xds":
            df = _add_reference_file(df, path, algo, ds)
        df = _standardise(df, algo)
        if need_bidomain and df["bd_mean"].notna().sum() == 0:
            continue
        df.attrs["source"] = str(path.relative_to(REPO))
        result = df
        break
    _cache[key] = result
    return result


def load_validation(algo, ds):
    """Like load(), but algorithms in REF_BORROW take their reference columns from another algorithm's
    reference results (matched on the instance pair) and `match` is recomputed against that."""
    df = load(algo, ds)
    if df is None or algo not in REF_BORROW:
        return df
    src = load(REF_BORROW[algo], ds)
    if src is None:
        return None
    r = src[["a", "b", "ref_size", "ref_nodes", "ref_time", "ref_done"]].rename(columns={"ref_done": "rd"})
    df = df.drop(columns=["ref_size", "ref_nodes", "ref_time", "ref_done"]).merge(r, on=["a", "b"], how="left")
    df["ref_done"] = df["rd"].eq(True)
    df = df.drop(columns="rd")
    df["match"] = _match(df)
    df.attrs["source"] = f"{algo}: unified from its own file; reference from {REF_BORROW[algo]}"
    return df


def available(ds, algos=None, include_unsound=False):
    """Algorithms (in ALGOS order) that have results for a dataset."""
    out = []
    for k in (algos or ALL_ALGOS):
        if not include_unsound and (k, ds) in UNSOUND:
            continue
        if load(k, ds) is not None:
            out.append(k)
    return out


def print_sources(ds, algos):
    """Say which file each algorithm's numbers come from, so mixed sources are never hidden."""
    print("   data files:")
    for k in algos:
        d = load(k, ds)
        if d is not None:
            print(f"     {k:9} <- {d.attrs.get('source', '?')}")


# --------------------------------------------------------------------------- difficulty classes
def classify(ds, core=None):
    """Easy / medium / hard for every instance pair of a dataset, as in the dissertation:
    easy   = every core algorithm solves it within EASY_S seconds
    hard   = no core algorithm solves it within TIMEOUT_S
    medium = everything else.
    Returns a DataFrame (a, b, cls) over the pairs present for ALL core algorithms."""
    wanted = core or CORE[ds]
    core = [k for k in wanted if load(k, ds) is not None]
    missing = [k for k in wanted if k not in core]
    if missing:
        print(f"WARNING: classify({ds}): no results for {', '.join(missing)} - "
              f"classes use only {', '.join(core)}. Rerun once the files exist.")
    wide = None
    for k in core:
        d = load(k, ds)[["a", "b", "done", "time"]].copy()
        d[k + "_easy"] = d["done"] & (d["time"] <= EASY_S)
        d[k + "_done"] = d["done"]
        d = d[["a", "b", k + "_easy", k + "_done"]]
        wide = d if wide is None else wide.merge(d, on=["a", "b"], how="inner")
    easy = wide[[k + "_easy" for k in core]].all(axis=1)
    hard = ~wide[[k + "_done" for k in core]].any(axis=1)
    wide["cls"] = np.where(easy, "easy", np.where(hard, "hard", "medium"))
    wide.attrs["core"] = core
    return wide[["a", "b", "cls"]]


# --------------------------------------------------------------------------- plotting
SAVE_DPI = 450   # resolution of every saved figure
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e1"


def style():
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
        "axes.edgecolor": INK2, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
        "text.color": INK, "axes.titlecolor": INK, "axes.titlesize": 10.5, "axes.labelsize": 9.5,
        "xtick.labelsize": 8.5, "ytick.labelsize": 8.5, "legend.fontsize": 8.5, "legend.frameon": False,
        "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
        "grid.linewidth": 0.6, "axes.axisbelow": True, "lines.linewidth": 1.8, "figure.dpi": 110,
    })
    return plt


def grid(n, panel=3.7):
    """Square panels in a balanced grid: 1-3 panels in one row, 4 in 2x2, 5-6 in 2x3, 7-8 in 2x4, more in rows
    of 4. Returns (fig, axes); axes in reading order, spare cells are hidden."""
    import matplotlib.pyplot as plt
    ncols = n if n <= 3 else 2 if n == 4 else 3 if n <= 6 else 4
    nrows = -(-n // ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(panel * ncols + 0.4, panel * nrows + 0.9), squeeze=False)
    flat = axes.ravel()
    for ax in flat[n:]:
        ax.axis("off")
    for ax in flat[:n]:
        ax.set_box_aspect(1)
    fig._grid_shape = (nrows, ncols, n)
    return fig, flat[:n]


def outer_labels(fig, axes):
    """Keep axis labels only on the left column / bottom row of a grid() figure (call after labelling)."""
    nrows, ncols, n = fig._grid_shape
    for i, ax in enumerate(axes):
        r, c = divmod(i, ncols)
        if c != 0:
            ax.set_ylabel("")
        if r != nrows - 1 and i + ncols < n:
            ax.set_xlabel("")


def save(fig, folder, name):
    p = out_dir(folder) / name
    fig.savefig(p, dpi=SAVE_DPI, bbox_inches="tight")
    print("  saved", p.relative_to(REPO))
    return p
