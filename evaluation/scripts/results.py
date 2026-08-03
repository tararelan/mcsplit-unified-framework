"""
4.2 lv Benchmark Results
4.2.1 Instance Difficulty Classification - easy/medium/hard counts
4.2.2 Instances Solved on Shared Moderate Set - solve table + cactus plot (time)
4.2.3 Per-Node Overhead Analysis - nodes vs time scatter
4.2.4 Hard Instance Analysis - incumbent sizes, flat result motivation
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

pd.set_option('display.width', 200)
pd.set_option('display.float_format', '{:.3f}'.format)

ALGOS = ['mcsplit', 'rl', 'll', 'dsb'] #, 'rrsplit', 'symsplit']
DATA_DIR = Path('hpc/cross-dataset-evaluation/results')

LABELS = {
    'mcsplit': 'McSplit', 'rl': 'McSplit+RL', 'll': 'McSplit+LL',
    'dal': 'McSplit+DAL', 'dsb': 'McSplit+DSB',
    'rrsplit': 'RRSplit', 'symsplit': 'SymSplit',
}
COLORS = {
    'mcsplit': '#555555', 'rl': '#e07b39', 'll': '#e0a800',
    'dal': '#c94040', 'dsb': '#7a5c9e',
    'rrsplit': '#2e86ab', 'symsplit': '#3a9e5f',
}

COLS = [
    'instance_a', 'instance_b', 'algo',
    'unified_size', 'unified_edges', 'unified_nodes', 'unified_time',
    'unified_aborted', 'unified_nodes_to_best',
    'unified_time_to_best', 'unified_cut_branches', 'unified_bound_pruned',
    'unified_sym_pruned',
    'ref_size', 'ref_nodes', 'ref_time', 'match'
]

TIMEOUT = 1000.0
EASY_THRESHOLD = 10.0
inst_key = ['instance_a', 'instance_b']

# ── Load ──────────────────────────────────────────────────────────────────────
dfs = []
for algo in ALGOS:
    path = DATA_DIR / f'{algo}_bi_all.csv'
    if not path.exists():
        print(f'WARNING: {path} not found, skipping {algo}')
        continue
    df = pd.read_csv(path, skipinitialspace=True)
    df = df.reindex(columns=COLS)
    dfs.append(df)
    print(len(df))

all_data = pd.concat(dfs, ignore_index=True)
for col in ['unified_size', 'unified_aborted', 'unified_time', 'unified_nodes',
            'ref_size', 'ref_time', 'ref_nodes']:
    all_data[col] = pd.to_numeric(all_data[col], errors='coerce')

# ── 4.2.1 Difficulty Classification ──────────────────────────────────────────
print('=== 4.2.1 Instance Difficulty Classification ===\n')
print('Easy:   all algorithms solve within 10s')
print('Hard:   no algorithm solves within 1000s')
print('Medium: everything else (primary comparison set)\n')

all_data['solve_time'] = all_data.apply(
    lambda r: r['unified_time'] if r['unified_aborted'] == 0 else np.nan, axis=1
)

time_pivot = all_data.pivot_table(
    index=inst_key, columns='algo', values='solve_time', aggfunc='first'
)
abort_pivot = all_data.pivot_table(
    index=inst_key, columns='algo', values='unified_aborted', aggfunc='first'
)
abort_any = (abort_pivot.fillna(1).astype(int) == 1).astype(int)
all_instances = all_data[inst_key].drop_duplicates().set_index(inst_key).index

def classify(idx):
    if idx not in abort_any.index:
        return 'hard'
    row_abort = abort_any.loc[idx]
    row_time  = time_pivot.loc[idx] if idx in time_pivot.index else pd.Series(dtype=float)
    n_solved = (row_abort == 0).sum()
    n_data   = row_abort.notna().sum()
    n_fast   = (row_time <= EASY_THRESHOLD).sum()
    if n_solved == 0:
        return 'hard'
    if n_fast == n_data:
        return 'easy'
    return 'medium'

difficulty = pd.Series(
    [classify(idx) for idx in all_instances],
    index=all_instances, name='difficulty'
)

print('Difficulty distribution:')
print(difficulty.value_counts())
# difficulty.value_counts().to_csv("evaluation/results/lv_difficulty_distribution.csv")

medium_instances = set(difficulty[difficulty == 'medium'].index)
hard_instances   = set(difficulty[difficulty == 'hard'].index)

med_mask  = all_data.set_index(inst_key).index.isin(medium_instances)
hard_mask = all_data.set_index(inst_key).index.isin(hard_instances)
medium_data = all_data[med_mask].copy()
hard_data   = all_data[hard_mask].copy()

# ── 4.2.2 Instances Solved on Shared Moderate Set ────────────────────────────
print(f'\n=== 4.2.2 Instances Solved on Shared Moderate Set (n={len(medium_instances)}) ===\n')

rows = []
total = len(medium_instances)
for algo in ALGOS:
    sub = medium_data[medium_data['algo'] == algo]
    if sub.empty:
        continue
    solved = int((sub['unified_aborted'] == 0).sum())
    rows.append({
        'Algorithm': LABELS[algo],
        'Solved':    solved,
        'Solved %':  round(100 * solved / total, 2),
    })

solved_df = pd.DataFrame(rows).set_index('Algorithm')
print(solved_df)
# solved_df.to_csv("evaluation/results/lv_solved_difficulty.csv")

# Cactus plot: time vs instances solved
fig, ax = plt.subplots(figsize=(9, 5))
for algo in ALGOS:
    sub = medium_data[
        (medium_data['algo'] == algo) &
        (medium_data['unified_aborted'] == 0)
    ]['unified_time'].sort_values().values
    if len(sub) == 0:
        continue
    ax.plot(sub, np.arange(1, len(sub) + 1),
            label=LABELS[algo], color=COLORS[algo], linewidth=1.8)

ax.set_xlabel('Wall-clock time (s)', fontsize=11)
ax.set_ylabel('Instances solved', fontsize=11)
# ax.set_xscale('log')
ax.set_title('Time versus Instances Solved', fontsize=11)
ax.legend(fontsize=9, loc='upper left')
ax.set_ylim(bottom=0)
ax.grid(True, alpha=0.3)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
# plt.savefig('evaluation/results/lv_cactus_time_with.png', bbox_inches='tight', dpi=500)
plt.show()
print('Saved lv_cactus_time_with.png')

# Cactus plot: nodes vs instances solved
fig, ax = plt.subplots(figsize=(9, 5))
for algo in ALGOS:
    sub = medium_data[
        (medium_data['algo'] == algo) &
        (medium_data['unified_aborted'] == 0)
    ]['unified_nodes'].sort_values().values
    if len(sub) == 0:
        continue
    ax.plot(sub, np.arange(1, len(sub) + 1),
            label=LABELS[algo], color=COLORS[algo], linewidth=1.8)

ax.set_xlabel('Nodes', fontsize=11)
ax.set_ylabel('Instances solved', fontsize=11)
# ax.set_xscale('log')
ax.set_title('Nodes versus Instances Solved', fontsize=11)
ax.legend(fontsize=9, loc='upper left')
ax.set_ylim(bottom=0)
ax.grid(True, alpha=0.3)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
# plt.savefig('evaluation/results/lv_cactus_nodes.png', bbox_inches='tight', dpi=500)
plt.show()
print('Saved lv_cactus_nodes.png')

# Cactus plot: time vs instances solved
fig, ax = plt.subplots(figsize=(9, 5))
for algo in ALGOS:
    sub = medium_data[
        (medium_data['algo'] == algo) &
        (medium_data['unified_aborted'] == 0)
    ]['unified_time'].sort_values().values
    if len(sub) == 0:
        continue
    ax.plot(sub, np.arange(1, len(sub) + 1),
            label=LABELS[algo], color=COLORS[algo], linewidth=1.8)

ax.set_xlabel('Wall-clock time (s)', fontsize=11)
ax.set_ylabel('Instances solved', fontsize=11)
ax.set_xscale('log')
# ax.set_yscale('log')
ax.set_title('Time versus Instances Solved', fontsize=11)
ax.legend(fontsize=9, loc='upper left')
ax.set_ylim(bottom=0)
ax.grid(True, alpha=0.3)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
# plt.savefig('evaluation/results/lv_cactus_time_log.png', bbox_inches='tight', dpi=500)
plt.show()
print('Saved lv_cactus_time_log.png')

# Cactus plot: nodes vs instances solved
fig, ax = plt.subplots(figsize=(9, 5))
for algo in ALGOS:
    sub = medium_data[
        (medium_data['algo'] == algo) &
        (medium_data['unified_aborted'] == 0)
    ]['unified_nodes'].sort_values().values
    if len(sub) == 0:
        continue
    ax.plot(sub, np.arange(1, len(sub) + 1),
            label=LABELS[algo], color=COLORS[algo], linewidth=1.8)

ax.set_xlabel('Nodes', fontsize=11)
ax.set_ylabel('Instances solved', fontsize=11)
ax.set_xscale('log')
# ax.set_yscale('log')
ax.set_title('Nodes versus Instances Solved', fontsize=11)
ax.legend(fontsize=9, loc='upper left')
ax.set_ylim(bottom=0)
ax.grid(True, alpha=0.3)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
# plt.savefig('evaluation/results/lv_cactus_nodes_log.png', bbox_inches='tight', dpi=500)
plt.show()
print('Saved lv_cactus_nodes_log.png')

# ── 4.2.3 Per-Node Overhead Analysis ─────────────────────────────────────────
print('\n=== 4.2.3 Per-Node Overhead Analysis ===\n')
print('Time per node (microseconds) - proxy for per-node cost.')
print('Higher = more work done per branching decision.\n')

completed = medium_data[medium_data['unified_aborted'] == 0].copy()
completed['time_per_node_us'] = (
    completed['unified_time'] / completed['unified_nodes'] * 1e6
)

overhead = completed.groupby('algo')['time_per_node_us'].agg(['median', 'mean'])
overhead.index = [LABELS[a] for a in overhead.index]
overhead.columns = ['Median µs/node', 'Mean µs/node']
print(overhead.round(3))
# overhead.to_csv("evaluation/results/overhead.csv")

completed = medium_data[medium_data['unified_aborted'] == 0].copy()
completed['time_per_node_us'] = completed['unified_time'] / completed['unified_nodes'] * 1e6

overhead = completed.groupby('algo')['time_per_node_us'].median().reindex(ALGOS)

fig, ax = plt.subplots(figsize=(8, 4))
ax.bar([LABELS[a] for a in ALGOS], overhead.values, color=[COLORS[a] for a in ALGOS])
ax.set_ylabel('Median time per node (µs)', fontsize=10)
ax.set_title('Per-Node Overhead (medium instances, lv)', fontsize=10)
ax.tick_params(axis='x', rotation=20)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
# plt.savefig('evaluation/results/lv_nodes_vs_time.png', bbox_inches='tight', dpi=500)
plt.show()
print('Saved lv_nodes_vs_time.png')

# ── 4.2.4 Hard Instance Analysis ─────────────────────────────────────────────
print(f'\n=== 4.2.4 Hard Instance Analysis (n={len(hard_instances)} hard instances) ===\n')
print('For hard instances, unified_size = best incumbent found before timeout.\n')

rows = []
for algo in ALGOS:
    sub = hard_data[hard_data['algo'] == algo]
    if sub.empty:
        continue
    sizes = sub['unified_size'].dropna()
    rows.append({
        'Algorithm':       LABELS[algo],
        'Hard instances':  len(sub),
        'Mean incumbent':  round(sizes.mean(), 2),
        'Median incumbent': round(sizes.median(), 2),
        'Max incumbent':   int(sizes.max()) if len(sizes) > 0 else np.nan,
    })

hard_df = pd.DataFrame(rows).set_index('Algorithm')
print(hard_df)
# hard_df.to_csv("evaluation/results/lv_hard_instances.csv")

# Check whether any algorithm consistently finds larger incumbents
# Pivot incumbent sizes and compute pairwise dominance
pivot = hard_data.pivot_table(
    index=inst_key, columns='algo', values='unified_size', aggfunc='first'
)

print('\nPairwise: fraction of hard instances where algo A finds strictly larger incumbent than algo B:')
print('(Only instances where both completed with a valid incumbent)\n')

header = [''] + [LABELS[a] for a in ALGOS]
rows_dom = []
for a in ALGOS:
    row = [LABELS[a]]
    for b in ALGOS:
        if a == b:
            row.append('-')
            continue
        if a not in pivot.columns or b not in pivot.columns:
            row.append('N/A')
            continue
        mask = pivot[a].notna() & pivot[b].notna()
        if mask.sum() == 0:
            row.append('N/A')
            continue
        frac = (pivot.loc[mask, a] > pivot.loc[mask, b]).mean()
        row.append(f'{frac:.2f}')
    rows_dom.append(row)

dom_df = pd.DataFrame(rows_dom, columns=header).set_index('')
print(dom_df)
# dom_df.to_csv("evaluation/results/lv_hard_dominates.csv")
print('\nNo algorithm consistently dominates on lv hard instances.')
print('This motivates cross-dataset evaluation on BI and LV (Section 4.4),')
print('where structural properties allow algorithms to meaningfully separate.')

all_data['unified_nodes_to_best'] = pd.to_numeric(all_data['unified_nodes_to_best'], errors='coerce')
all_data['unified_time_to_best'] = pd.to_numeric(all_data['unified_time_to_best'], errors='coerce')
all_data['unified_nodes'] = pd.to_numeric(all_data['unified_nodes'], errors='coerce')

# for dataset_name, data, instances in [
#     ('lv', all_data, medium_data),
# ]:
#     sub = data[data['algo'].isin(['rrsplit', 'symsplit'])].copy()
    
#     # Pivot nodes_to_best and time_to_best
#     nodes_pivot = sub.pivot_table(
#         index=inst_key, columns='algo', 
#         values='unified_nodes_to_best', aggfunc='first'
#     )
#     time_pivot = sub.pivot_table(
#         index=inst_key, columns='algo',
#         values='unified_time_to_best', aggfunc='first'
#     )
    
#     # Only instances where both have data
#     mask = nodes_pivot['rrsplit'].notna() & nodes_pivot['symsplit'].notna()
#     nodes_pivot = nodes_pivot[mask]
#     time_pivot = time_pivot[mask]
    
#     print(f'\n=== {dataset_name} ===')
#     print(f'Instances compared: {mask.sum()}')
    
#     print('\nNodes to best incumbent:')
#     print(f'  SymSplit median:  {nodes_pivot["symsplit"].median():.0f}')
#     print(f'  RRSplit median:   {nodes_pivot["rrsplit"].median():.0f}')
#     print(f'  SymSplit finds incumbent in fewer nodes: '
#           f'{(nodes_pivot["symsplit"] < nodes_pivot["rrsplit"]).sum()}/{mask.sum()} instances')
    
#     print('\nTime to best incumbent:')
#     print(f'  SymSplit median:  {time_pivot["symsplit"].median():.3f}s')
#     print(f'  RRSplit median:   {time_pivot["rrsplit"].median():.3f}s')
#     print(f'  SymSplit finds incumbent faster: '
#           f'{(time_pivot["symsplit"] < time_pivot["rrsplit"]).sum()}/{mask.sum()} instances')
    
#     # Total nodes explored comparison on instances both solved
#     total_nodes = sub.pivot_table(
#         index=inst_key, columns='algo',
#         values='unified_nodes', aggfunc='first'
#     )
    
#     # Both solved (not aborted)
#     abort_pivot = sub.pivot_table(
#         index=inst_key, columns='algo',
#         values='unified_aborted', aggfunc='first'
#     )
#     both_solved = (abort_pivot['rrsplit'] == 0) & (abort_pivot['symsplit'] == 0)
#     total_nodes_solved = total_nodes[both_solved]
    
#     print(f'\nTotal nodes explored (instances both solved, n={both_solved.sum()}):')
#     print(f'  SymSplit median:  {total_nodes_solved["symsplit"].median():.0f}')
#     print(f'  RRSplit median:   {total_nodes_solved["rrsplit"].median():.0f}')
#     print(f'  SymSplit explores fewer total nodes: '
#           f'{(total_nodes_solved["symsplit"] < total_nodes_solved["rrsplit"]).sum()}/{both_solved.sum()} instances')
    
#     # Also check time per node to confirm the per-node cost explanation
#     time_solved = sub[sub['unified_aborted'] == 0].pivot_table(
#         index=inst_key, columns='algo',
#         values='unified_time', aggfunc='first'
#     )
#     nodes_solved = total_nodes_solved
#     common = time_solved.index.intersection(nodes_solved.index)
    
#     sym_tpn = (time_solved.loc[common, 'symsplit'] / nodes_solved.loc[common, 'symsplit'] * 1e6).median()
#     rr_tpn  = (time_solved.loc[common, 'rrsplit']  / nodes_solved.loc[common, 'rrsplit']  * 1e6).median()
#     print(f'\nMedian time per node on solved instances (µs):')
#     print(f'  SymSplit: {sym_tpn:.3f}')
#     print(f'  RRSplit:  {rr_tpn:.3f}')