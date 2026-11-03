# analysis/

One place for every figure and table in the dissertation (and the paper). It replaces the many scripts in
`evaluation/scripts/` and `hpc/validation/scripts/`, which are left untouched.

```
python run_all.py                       # everything, all datasets  -> outputs/
python 01_validation.py --datasets bi   # one script, one dataset
```

Needs `python3`, `pandas`, `numpy`, `matplotlib`. Run from anywhere: the repo root is found automatically.

## What is where

| File | Question it answers | Dissertation | Output (`outputs/…`) |
|---|---|---|---|
| `config.py` | Shared settings: algorithms, colours, data locations, loaders. **Edit this, not the scripts.** | – | – |
| `01_validation.py` | Does each re-implementation reproduce its reference (solution size, search nodes, runtime)? | 4.1, Tables 3, 5, 6, Fig 15–22 | `validation/` |
| `02_performance.py` | How do the algorithms compare on each dataset? Easy/medium/hard, solved on the medium set, cactus plots, time per node, hard-instance incumbents | 4.2–4.5, Tables 7–10, Fig 23–24 | `performance/` |
| `03_ablation.py` | What does each component (LSM/LUM, DSB, RRSplit reductions, SymSplit symmetries) contribute? | 4.3, Table 11, Fig 25 | `ablation/` |
| `04_bidomain.py` | How many bidomains per search node, per dataset and algorithm? | 4.1.3 | `bidomain/` |
| `run_all.py` | Runs 01–04 in order | – | – |

Each script prints its tables to the terminal and saves them as CSV next to the figures. Every script starts
with a docstring saying what it does and what it writes.

## Algorithms and colours

Colour = algorithm family; line style and marker = variant within the family. All defined in `ALGOS` in `config.py`.

| Family | Algorithms | Colour |
|---|---|---|
| McSplit | McSplit (solid), McSplit-down (dashed) | grey |
| RL | McSplit+RL (solid), RL Pair (dashed), RL Rate (dotted), RL-down (dash-dot) | blue |
| LL | McSplit+LL | green |
| DAL | McSplit+DAL (star marker) | orange |
| DSB | McSplit+DSB | yellow |
| RRSplit | RRSplit | violet |
| SymSplit | SymSplit | magenta |

To add an algorithm: add one line to `ALGOS` (and to `CORE` in `config.py` only if it should define easy/medium/hard).

## Where the data comes from

`config.py` → `SOURCES` lists, per dataset, the result folders in order of preference; the first file that exists
wins, and every script prints which file each algorithm's numbers came from. To use newer results for a dataset,
move its folder to the top of that dataset's list.

| Dataset | Preferred source | Falls back to |
|---|---|---|
| BI | `hpc/publication/results` (current pipeline) | `hpc/publication/old`, `hpc/cross-dataset-evaluation` |
| MIVIA | `hpc/validation/results` (dissertation runs) | `hpc/publication/*` |
| LV | `hpc/cross-dataset-evaluation/results` (dissertation runs) | `hpc/publication/*` |
| Ablation (MIVIA) | `hpc/ablation/results` | – |

The loader reads every file layout into one standard table (documented at the top of `config.py`), stitches
together DSB/SymSplit rows that an older `run_instance.sh` split in two (in memory; the files are not changed), and
recomputes PASS / FAIL / MISSING with the dissertation's rule.

## Rules baked in (all in `config.py`)

- **PASS** = the reference finished and both solution sizes are identical. **FAIL** = both finished within 1000 s but
  sizes differ. **MISSING** = anything else (a timeout on either side).
- **Easy** = every core algorithm solves it within 10 s; **hard** = no core algorithm solves it; **medium** = the rest.
  Core algorithms = McSplit, RL, LL, DSB, RRSplit, SymSplit (BI: without RRSplit and SymSplit).
- RRSplit and SymSplit on BI are left out by default because their references are unsound on directed graphs
  (`--include-unsound` shows them).
- `rl_pair`, `rl_rate`, `rl_down` have no reference binary, so they are compared with the **RL** reference. That is
  a comparison with RL, not a fidelity check, and the outputs say so.

## Things to know when reading the outputs

- McSplit-down and RL-down search top-down and only hold an incumbent when a search finishes, so on timed-out (hard)
  instances their reported solution size is far below the others'. The hard-instance incumbent matrix
  (`02_performance.py`) reflects that; it is a property of the strategy, not a bug.
- A pair missing from an algorithm's result file counts as unsolved by that algorithm.
- The bidomain counters are recorded in every run, so reported times include that (small) overhead.
