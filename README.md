# The University of St Andrews

## Understanding and Evaluating Variants of the McSplit Algorithm for the Maximum Common Subgraph Problem

### Tara Relan (250016623)
### Supervisor: Dr. Ruth Hoffmann
### August 11th 2026

## Contents

- [Abstract](#abstract)
- [Acknowledgements](#acknowledgements)
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
    - [4. Collect results](#4-collect-results)
	- [5. Copy results back to local machine](#5-copy-results-back-to-local-machine)
  - [Adding a New Algorithm](#adding-a-new-algorithm)

## Abstract

## Acknowledgements

This repository contains reimplementations of the following algorithms:

- **McSplit** — McCreesh, Prosser & Trimble. [A Partitioning Algorithm for Maximum Common Subgraph Problems](https://doi.org/10.24963/ijcai.2017/99). IJCAI 2017.
- **McSplit+RL** — Liu, Li, Jiang & He. A Learning Based Branch and Bound for Maximum Common Subgraph Related Problems. AAAI 2020.
- **McSplit+LL** — Zhou, He, Zheng, Li & Liu. A Strengthened Branch and Bound Algorithm for the Maximum Common (Connected) Subgraph Problem. IJCAI 2020.
- **McSplit+DAL** — Liu, Zhao, Li, Jiang & He. [Hybrid Learning with New Value Function for the Maximum Common Induced Subgraph Problem](https://doi.org/10.1609/aaai.v37i4.25519). AAAI 2023.
- **McSplit+DSB** — Kothalawala, Koehler & Wang. [Learning to Bound for Maximum Common Subgraph Algorithms](https://doi.org/10.4230/LIPIcs.CP.2025.22). CP 2025.
- **RRSplit** — Yu, Wang, Long, Lakshmanan & Cheng. Fast Maximum Common Subgraph Search: A Redundancy-Reduced Backtracking Approach. 2024.
- **SymSplit** — Kothalawala, Koehler & Farhan. [Accelerating Maximum Common Subgraph Computation by Exploiting Symmetries](https://doi.org/10.1145/3802005). SIGMOD 2026.

All implementations were written from scratch in C++ based on the published papers and validated against the authors' reference binaries. Reference binaries are included in `reference/` for validation purposes only and remain the property of their respective authors.

## User Manual

This manual describes how to build and run the unified MCS solver framework. The framework implements seven algorithms (McSplit, McSplit+RL, McSplit+LL, McSplit+DAL, McSplit+DSB, RRSplit, SymSplit) and their ablation variants within a single binary.

---

### Requirements

- g++ 11.4.0 or later with C++17 support
- instances.zip unzipped
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
| `unified_cut_branches` | Branches pruned by the bound |
| `unified_bound_pruned` | Branches pruned by algorithm-specific bound |
| `unified_sym_pruned` | Branches pruned by symmetry or redundancy reduction |

---

### Running Experiments on HPC (Hypatia)

#### 1. Sync files to Hypatia

```bash
rsync -avu --progress \
  --exclude='instances/SIP/images-CVIU11' \
  --exclude='instances/SIP/images-PR15' \
  --exclude='instances/SIP/meshes-CVIU11' \
  --exclude='instances/SIP/phase' \
  --exclude='instances/SIP/scalefree' \
  --exclude='instances/SIP/si' \
  "[source-directory]/" \
  [account]@hypatia.st-andrews.ac.uk:work/[target-directory]/
```

#### 2. Build on Hypatia

```bash
cd ~/work/[target-directory]/solvers && rm -f bin/mcsp && make
```

#### 3. Submit scripts

You should submit one algorithm at a time so as not to hit the QOS limits.

```bash
cd ~/work/[target-directory]
bash hpc/submit_mivia.sh mcsplit
# wait for completion (squeue -u [account]), then submit next
bash hpc/submit_bi.sh dsb
bash hpc/submit_lv.sh symsplit
```

To include a reference binary comparison:
```bash
bash hpc/submit_mivia.sh mcsplit --with-ref
bash hpc/submit_bi.sh dsb --with-ref
bash hpc/submit_lv.sh symsplit --with-ref
```

rm -rf hpc/results/symsplit_lv_reference/

bash -n hpc/job_array.slurm 2>&1 | head  # won't fully validate since placeholders aren't substituted, but catches gross syntax issues
bash hpc/submit_bi.sh dsb              # plain mode
bash hpc/submit_bi.sh dsb --with-ref   # ref mode

To see progress
squeue -u tr77
ls hpc/results/mcsplit_mivia/ | wc -l
head -5 -q hpc/results/mcsplit_mivia/*.csv
ls hpc/results/dsb_bi_with_ref/ | wc -l
head -5 -q hpc/results/dsb_bi_with_ref/*.csv
ls hpc/results/symsplit_lv_reference/ | wc -l
head -5 -q hpc/results/symsplit_lv_reference/*.csv
head -5 -q hpc/results/dsb_bi/*.csv
ls hpc/results/mcsplit_mivia/ | wc -l

The output would be one line, space-separated:
```
<unified_size> <unified_edges> <unified_nodes> <unified_time> <unified_aborted> <unified_nodes_to_best> <unified_time_to_best> <unified_cut_branches> <unified_bound_pruned> <unified_sym_pruned> <ref_size> <ref_nodes> <ref_time> <match>
```

| Field | Description |
|-------|-------------|
| `ref_size` | Reference solution size (MCS vertex count) |
| `ref_nodes` | Reference nodes explored in search tree |
| `ref_time` | Reference wall-clock time |

#### 4. Collect results

```bash
for ALGO in mcsplit rl ll ll_lsm ll_lum dal dal_rl dal_dal \
            dsb dsb_always dsb_never \
            rrsplit rrsplit_noveq rrsplit_nomax rrsplit_nobound \
            symsplit symsplit_valonly symsplit_varonly; do
    # MIVIA
    cat hpc/results/${ALGO}_mivia/*.csv > hpc/results/${ALGO}_mivia_all.csv
    # BI
    cat hpc/results/${ALGO}_bi/*.csv > hpc/results/${ALGO}_bi_all.csv
    # LV
    cat hpc/results/${ALGO}_lv/*.csv > hpc/results/${ALGO}_lv_all.csv
done
```

bash hpc/merge_results.sh dsb bi              # plain mode → dsb_bi_all.csv
bash hpc/merge_results.sh symsplit lv --with-ref   # ref mode   → dsb_bi_with_ref_all.csv

cat hpc/results/mcsplit_mivia/*.csv > hpc/results/mcsplit_mivia_all.csv

Note that merging the CSVs do not have any headings - you will have to add them yourself.

instance_a,instance_b,algo,unified_size,unified_edges,unified_nodes,unified_time,unified_aborted,unified_nodes_to_best,unified_time_to_best,unified_cut_branches,unified_bound_pruned,unified_sym_pruned

instance_a,instance_b,algo,unified_size,unified_edges,unified_nodes,unified_time,unified_aborted,nified_nodes_to_best,unified_time_to_best,unified_cut_branches,unified_bound_pruned,unified_sym_pruned,ref_size,ref_nodes,ref_time,match

#### 5. Copy results back to local machine

```bash
scp [account]@hypatia.st-andrews.ac.uk:work/[target-directory]/hpc/results/*_all.csv \
    "[source-directory]/"
```

---

### Adding a New Algorithm

1. Implement the algorithm in a new `.cpp` file in `solvers/src/`, following the pattern of existing implementations (entry point returning `std::vector<VtxPair>`).
2. Register the dispatch in `solvers/src/main.cpp` - add an `else if (algo == "my_algo")` branch.
3. Add the file to the compiler command in `solvers/Makefile`.
4. If a reference binary exists, add a case for it in `hpc/run_instance.sh` under the `WITH_REF=1` block.
5. Rebuild and validate:

```bash
make -C solvers
bash hpc/submit_mivia.sh my_algo --with-ref
```