# The University of St Andrews

## Understanding and Evaluating Variants of the McSplit Algorithm for the Maximum Common Subgraph Problem

### Tara Relan (250016623)
### Supervisor: Dr. Ruth Hoffmann
### August 11th 2026

## Contents

- [Abstract](#abstract)
- [Acknowledgements](#acknowledgements)
- [Reference Implementations](#reference-implementations)
- [User Manual](#user-manual)
  - [Requirements](#requirements)
  - [Building](#building)
  - [Running a Single Instance](#running-a-single-instance)
    - [Options](#options)
    - [Available Algorithms](#available-algorithms)
    - [Examples](#examples)
    - [Output Format](#output-format)
  - [Running Experiments on HPC](#running-experiments-on-hpc-hypatia)
    - [1. Sync files to Hypatia](#1-sync-files-to-hypatia)
    - [2. Build on Hypatia](#2-build-on-hypatia)
    - [3. Submit scripts](#3-submit-scripts)
    - [4. Monitoring progress](#4-monitoring-progress)
    - [5. Collect and merge results](#5-collect-and-merge-results)
    - [6. Copy results back to local machine](#6-copy-results-back-to-local-machine)
  - [Adding a New Algorithm](#adding-a-new-algorithm)

<!-- ## Abstract

The Maximum Common Subgraph problem is a fundamental, NP-hard graph-matching problem with applications in cheminformatics, malware detection, and pattern recognition. McSplit and its many variants --- incorporating reinforcement-learning-guided branching, tighter bounding, and structural symmetry reduction --- represent the current state of the art, but each has been evaluated independently, under different hardware, timeouts, and benchmark suites, making it unclear how much of any claimed improvement is attributable to algorithm design rather than experimental setup.

This dissertation presents open-source C++ re-implementations of seven McSplit-family algorithms within a shared search skeleton, isolating each algorithm's distinguishing components --- branching heuristic, bound computation, and symmetry or redundancy reduction --- for direct, controlled comparison. Each re-implementation is validated against its original reference for solution correctness, search-tree fidelity, and runtime fidelity across synthetic, directed biochemical, and large heterogeneous benchmark datasets.

This validation uncovered two previously undocumented defects in published reference implementations --- a solution-inflating reward bug in McSplit+DAL, and a direction-awareness omission rendering RRSplit and SymSplit unsound on directed graphs --- both traced to specific source-level causes and corroborated across independent implementations. This shows that prior evaluations of these algorithms, including their own published comparisons, may not be fully reliable, and demonstrates the necessity of unified, comparable re-implementation. Beyond validation, cross-dataset evaluation shows no single algorithm dominates: structural pruning outperforms learned branching on synthetic and heterogeneous graphs, while reinforcement-learning-guided branching becomes competitive on directed biochemical networks, where repeated reaction motifs provide an exploitable reward signal. An accompanying ablation study isolates each component's individual contribution, offering practical, structure-aware guidance for solver selection. -->

## Acknowledgements

This repository contains reimplementations of the following algorithms:

- **McSplit** - C. McCreesh and P. Prosser and J. Trimble. [A Partitioning Algorithm for Maximum Common Subgraph Problems](https://doi.org/10.24963/ijcai.2017/99). IJCAI 2017.
- **McSplit+RL** - Y. Liu and C.-M. Li and H. Jiang and K. He. [A Learning Based Branch and Bound for Maximum Common Subgraph Related Problems](https://doi.org/10.1609/aaai.v34i03.5619). AAAI 2020.
- **McSplit+LL** - J. Zhou and K. He and J. Zheng and C.-M. Li and Y. Liu. [A Strengthened Branch and Bound Algorithm for the Maximum Common (Connected) Subgraph Problem](https://doi.org/10.24963/ijcai.2022/265). IJCAI 2022.
- **McSplit+DAL** - Y. Liu and J. Zhao and C.-M. Li and H. Jiang and K. He. [Hybrid Learning with New Value Function for the Maximum Common Induced Subgraph Problem](https://doi.org/10.1609/aaai.v37i4.25519). AAAI 2023.
- **McSplit+DSB** - B. W. Kothalawala and H. Koehler and Q. Wang. [Learning to Bound for Maximum Common Subgraph Algorithms](https://doi.org/10.4230/LIPIcs.CP.2025.22). CP 2025.
- **RRSplit** - K. Yu and K. Wang and C. Long and L.V.S. Lakshmanan and R. Cheng. [Fast Maximum Common Subgraph Search: A Redundancy-Reduced Backtracking Approach](https://doi.org/10.1145/3725404). ACM. 2024.
- **SymSplit** - B. W. Kothalawala and H. Koehler and M. Farhan. [Accelerating Maximum Common Subgraph Computation by Exploiting Symmetries](https://doi.org/10.1145/3802005). ACM 2026.

All implementations were re-written in C++ based on the published papers and validated against the authors' reference binaries. Reference binaries are included in `reference/` for validation purposes only and remain the property of their respective authors.

## Reference Implementations

Original author code used for validation (Section 4.1 of the dissertation) was sourced as follows:

| Algorithm | Source |
|-----------|--------|
| McSplit | [ijcai2017-partitioning-common-subgraph](https://github.com/jamestrimble/ijcai2017-partitioning-common-subgraph) |
| McSplit+RL | [McSplit-RL](https://github.com/JHL-HUST/McSplit-RL) |
| McSplit+LL | [McSplit-LL](https://github.com/JHL-HUST/McSplit-LL) |
| McSplit+DAL | [SIGMOD25-MCSS](https://github.com/KaiqiangYu/SIGMOD25-MCSS) (`./McSplitDAL` subfolder) |
| McSplit+DSB | [mcsplit-dsb](https://github.com/BuddhiWathsala/mcsplit-dsb) |
| RRSplit | [SIGMOD25-MCSS](https://github.com/KaiqiangYu/SIGMOD25-MCSS) (`./RRSplit` subfolder) |
| SymSplit | [symsplit](https://github.com/mcsolver/symsplit) |

Your local `reference/` folder should be laid out as:

```
reference/
├── mcsplit/
│   └── mcsp      # compiled reference binary
├── mcsplit-rl/
│   └── mcsp+RL
├── mcsplit-ll/
│   └── mcsp+ll
├── mcsplit-dal/
│   └── mcspDAL
├── mcsplit-dsb/
│   └── bin/
│       └── run.o
├── rrsplit/
│   └── mcsp
└── symsplit/
    └── bin/
        └── run.o
```

These paths match what `hpc/run_instance.sh` calls in `--with-ref` mode - if a reference binary isn't at the expected path, that algorithm's ref-mode comparison will silently return empty reference fields rather than erroring, so it's worth checking each path exists after building.

## User Manual

This manual describes how to build and run the unified MCS solver framework. The framework implements seven algorithms (McSplit, McSplit+RL, McSplit+LL, McSplit+DAL, McSplit+DSB, RRSplit, SymSplit) and their ablation variants within a single binary.

---

### Requirements

- g++ 11.4.1 or later with C++17 support
- The [MIVIA/BI/LV instance set](https://github.com/tararelan/mcsplit-unified-framework/releases/tag/instances-v1), unzipped into `instances/`
- GNU Make
- POSIX-compatible shell (bash)
- Boost library (for HPC submission scripts only)

---

### Building

From the project root:

```bash
make -C solvers
```

This produces `solvers/bin/mcsp`. To rebuild from scratch:

```bash
rm -f solvers/bin/mcsp
make -C solvers
```

---

### Running a Single Instance

```bash
./solvers/bin/mcsp -A <algorithm> -t <timeout> [-l] -q <graph1> <graph2>
```

#### Options

| Flag | Long form | Description |
|------|-----------|-------------|
| `-q` | `--quiet` | Suppress per-instance output; print only the solution size |
| `-d` | `--dimacs` | Read graphs in DIMACS format (default: MIVIA binary) |
| `-l` | `--lad` | Read graphs in LAD format (default: MIVIA binary) |
| `-i` | `--directed` | Treat graphs as directed |
| `-a` | `--labelled` | Use both vertex and edge labels |
| `-x` | `--vertex-labelled` | Use vertex labels only (no edge labels) |
| `-t <seconds>` | `--timeout` | Abort search after this many seconds (0 = no limit) |
| `-A <algo>` | `--algorithm` | Algorithm to run (see algorithm table below) |

#### Available Algorithms

| Flag | Algorithm |
|------|-----------|
| `mcsplit` | McSplit (baseline) |
| `rl` | McSplit+RL |
| `ll` | McSplit+LL (full: LSM + LUM) |
| `ll_lsm` | McSplit+LL with LSM only |
| `ll_lum` | McSplit+LL with LUM only |
| `dal` | McSplit+DAL (hybrid alternating policy) |
| `dal_rl` | McSplit+DAL with RL policy locked |
| `dal_dal` | McSplit+DAL with DAL policy locked |
| `dsb` | McSplit+DSB (RL-gated bound) |
| `dsb_always` | McSplit+DSB with bound always on |
| `dsb_never` | McSplit+DSB with bound always off |
| `rrsplit` | RRSplit (all three reductions) |
| `rrsplit_noveq` | RRSplit without vertex-equivalence reduction |
| `rrsplit_nomax` | RRSplit without maximality reduction |
| `rrsplit_nobound` | RRSplit without tighter bound |
| `symsplit` | SymSplit (variable + value symmetry) |
| `symsplit_varonly` | SymSplit with variable symmetry only |
| `symsplit_valonly` | SymSplit with value symmetry only |

#### Examples

**MIVIA binary format:**
```bash
./solvers/bin/mcsp -A mcsplit -t 1000 -q \
    instances/MIVIA/mcs10/bvg/b03/mcs10_b03_s40.A00 \
    instances/MIVIA/mcs10/bvg/b03/mcs10_b03_s40.B00
```

**SIP LAD format (BI - biochemical reactions):**
```bash
./solvers/bin/mcsp -A rrsplit -t 1000 -l -q \
    instances/SIP/biochemicalReactions/001.txt \
    instances/SIP/biochemicalReactions/002.txt
```

**SIP LAD format (LV):**
```bash
./solvers/bin/mcsp -A dal_dal -t 0 -l -q \
    instances/SIP/LV/g2 \
    instances/SIP/LV/g5
```

#### Output Format

One line, space-separated:
```
<unified_size> <unified_edges> <unified_nodes> <unified_time> <unified_aborted> <unified_nodes_to_best> <unified_time_to_best> <unified_cut_branches> <unified_bound_pruned> <unified_sym_pruned>
```

| Field | Description |
|-------|-------------|
| `unified_size` | Solution size (MCS vertex count) |
| `unified_edges` | Edges in the common subgraph |
| `unified_nodes` | Search tree nodes explored |
| `unified_time` | Wall-clock time in seconds |
| `unified_aborted` | 1 if timed out before proving optimality, 0 if completed |
| `unified_nodes_to_best` | Nodes explored before best incumbent found |
| `unified_time_to_best` | Time at which best incumbent was found (seconds) |
| `unified_cut_branches` | Total branches pruned |
| `unified_bound_pruned` | Branches pruned by the algorithm's bound |
| `unified_sym_pruned` | Branches pruned by symmetry or redundancy reduction |

When run with `--with-ref` on HPC (see below), four extra fields are appended: `ref_size`, `ref_nodes`, `ref_time`, `match` (`PASS`/`FAIL`/`MISSING_DATA`).

---

### Running Experiments on HPC (Hypatia)

#### 1. Sync files to Hypatia

```bash
rsync -avu --progress \
  "[source-directory]/" \
  [account]@hypatia.st-andrews.ac.uk:work/[target-directory]/
```

#### 2. Connect to Hypatia

```bash
ssh [account]@hypatia.st-andrews.ac.uk
```

#### 3. Build on Hypatia

```bash
cd ~/work/[target-directory]/solvers && rm -f bin/mcsp && make
```

#### 4. Submit scripts

Submit one algorithm at a time so as not to hit the QOS limits. Each script accepts an optional `--with-ref` flag to also run the corresponding reference binary and record a PASS/FAIL/MISSING_DATA comparison alongside the unified solver's results.

```bash
cd ~/work/[target-directory]

# plain mode
bash hpc/submit_mivia.sh mcsplit
bash hpc/submit_bi.sh dsb
bash hpc/submit_lv.sh symsplit

# with reference comparison
bash hpc/submit_mivia.sh mcsplit --with-ref
bash hpc/submit_bi.sh dsb --with-ref
bash hpc/submit_lv.sh symsplit --with-ref
```

Wait for each array to finish (`squeue -u [account]`) before submitting the next, to stay within QOS limits.

#### 5. Monitoring progress

```bash
squeue -u [account]                                    # job status
ls hpc/results/<algo>_<dataset>/ | wc -l               # instances completed (plain mode)
ls hpc/results/<algo>_<dataset>_with_ref/ | wc -l      # instances completed (ref mode)
head -5 -q hpc/results/<algo>_<dataset>/*.csv          # spot-check output
```

e.g. `ls hpc/results/mcsplit_mivia/ | wc -l` or `ls hpc/results/dsb_bi_with_ref/ | wc -l`.

#### 6. Collect and merge results

Per-instance CSVs (one file per instance pair, no header) are merged into a single file with the correct header using `hpc/merge_results.sh`:

```bash
bash hpc/merge_results.sh <algo> <dataset> [--with-ref]
```

e.g.:
```bash
bash hpc/merge_results.sh dsb bi                   # → hpc/results/dsb_bi_all.csv
bash hpc/merge_results.sh symsplit lv --with-ref   # → hpc/results/symsplit_lv_with_ref_all.csv
```

Plain-mode header (13 columns):
```
instance_a,instance_b,algo,unified_size,unified_edges,unified_nodes,unified_time,unified_aborted,unified_nodes_to_best,unified_time_to_best,unified_cut_branches,unified_bound_pruned,unified_sym_pruned
```

Ref-mode header (17 columns):
```
instance_a,instance_b,algo,unified_size,unified_edges,unified_nodes,unified_time,unified_aborted,unified_nodes_to_best,unified_time_to_best,unified_cut_branches,unified_bound_pruned,unified_sym_pruned,ref_size,ref_nodes,ref_time,match
```

#### 7. Copy results back to local machine

```bash
scp [account]@hypatia.st-andrews.ac.uk:work/[target-directory]/hpc/results/*_all.csv \
    "[source-directory]/"
```

---

### Adding a New Algorithm

1. Implement the algorithm in a new `.cpp` file in `solvers/src/`, following the pattern of existing implementations (entry point returning `std::vector<VtxPair>`).
2. Register the dispatch in `solvers/src/main.cpp` - add an `else if (algo == "my_algo")` branch.
3. Add the file to the compiler command in `solvers/Makefile`.
4. If a reference binary exists, add it under `reference/my_algo/` (see [Reference Implementations](#reference-implementations) for the expected layout) and add a case for it in `hpc/run_instance.sh` under the `WITH_REF=1` block.
5. Rebuild and validate:

```bash
make -C solvers
bash hpc/submit_mivia.sh my_algo --with-ref
```