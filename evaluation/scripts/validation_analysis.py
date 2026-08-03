"""
4.1 Validation Against Reference Implementations
4.1.1 Solution Correctness  - agreement table (PASS/FAIL/MISSING)
4.1.2 Search Tree Fidelity  - nodes scatter (unified vs reference)
4.1.3 Runtime Fidelity      - time scatter + speedup table
4.1.4 McSplit-DAL Bug       - +1 inflation in reference binary
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

pd.set_option('display.width', 200)
pd.set_option('display.float_format', '{:.3f}'.format)

ALGOS = ['mcsplit', 'rl', 'll', 'dsb', 'rrsplit', 'symsplit']
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

# ── Load ──────────────────────────────────────────────────────────────────────
dfs = []
for algo in ALGOS:
    path = DATA_DIR / f'{algo}_lv_all.csv'
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

# print(all_data)

# Reclassify borderline FAILs where unified timed out AND ref was near timeout
BORDERLINE_THRESHOLD = 850.0

def reclassify(row):
    if row['match'] == 'FAIL' and row['unified_aborted'] == 1:
        if pd.notna(row['ref_time']) and row['ref_time'] > BORDERLINE_THRESHOLD:
            return 'MISSING_DATA'
    return row['match']

all_data['match_adj'] = all_data.apply(reclassify, axis=1)

# print(all_data)

# # ── 4.1.1 Solution Correctness ────────────────────────────────────────────────
# print('=== 4.1.1 Solution Correctness ===\n')

# rows = []
# for algo in ALGOS:
#     sub = all_data[all_data['algo'] == algo]
#     if sub.empty:
#         continue
#     counts = sub['match_adj'].value_counts()
#     total = len(sub)
#     rows.append({
#         'Algorithm':    LABELS[algo],
#         'PASS':         int(counts.get('PASS', 0)),
#         'FAIL':         int(counts.get('FAIL', 0)),
#         'MISSING_DATA': int(counts.get('MISSING_DATA', 0)),
#         'PASS %':       round(100 * counts.get('PASS', 0) / total, 1),
#     })

# correctness_df = pd.DataFrame(rows).set_index('Algorithm')
# print(correctness_df)
# print('\nAll algorithms except McSplit+DAL report zero real mismatches.')

# ── 4.1.4 McSplit-DAL Reference Bug ──────────────────────────────────────────
# print('\n=== 4.1.4 McSplit+DAL Reference Bug ===\n')

# dal_fails = all_data[
#     (all_data['algo'] == 'dal') & (all_data['match_adj'] == 'FAIL')
# ].copy()

# print(f'DAL FAIL count: {len(dal_fails)}')

# if not dal_fails.empty:
#     dal_fails['size_diff'] = dal_fails['unified_size'] - dal_fails['ref_size']
#     print('\nDistribution of (unified_size - ref_size) on FAIL instances:')
#     print(dal_fails['size_diff'].value_counts().sort_index())
#     print('\nThe reference binary consistently reports one vertex more than the')
#     print('unified reimplementation - a +1 inflation bug in the reference binary.')
#     print('The unified reimplementation is treated as correct.')

# ── 4.1.2 Search Tree Fidelity ────────────────────────────────────────────────
print('\n=== 4.1.2 Search Tree Fidelity ===')

pass_data = all_data[
    (all_data['match_adj'] == 'PASS') &
    (all_data['unified_nodes'].notna()) &
    (all_data['ref_nodes'].notna()) &
    (all_data['ref_nodes'] > 0)
].copy()

# print(pass_data)
# print(pass_data.columns)

print('\nPearson correlation (unified nodes vs reference nodes):')
for algo in ALGOS:
    sub = pass_data[pass_data['algo'] == algo]
    # print(sub)
    if sub.empty:
        continue
    r = sub['unified_nodes'].corr(sub['ref_nodes'])
    print(f'  {LABELS[algo]}: r = {r:.4f}  (n={len(sub)})')

fig, axes = plt.subplots(2, 3, figsize=(14, 7))
axes = axes.flatten()
for i, algo in enumerate(ALGOS):
    ax = axes[i]
    sub = pass_data[pass_data['algo'] == algo]
    if sub.empty:
        ax.set_visible(False)
        continue
    ax.scatter(sub['ref_nodes'], sub['unified_nodes'],
               alpha=0.25, s=5, color=COLORS[algo], rasterized=True)
    lim = max(sub['ref_nodes'].max(), sub['unified_nodes'].max()) * 1.1
    ax.plot([1, lim], [1, lim], 'k--', linewidth=0.8, alpha=0.5)
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_title(LABELS[algo], fontsize=9)
    ax.set_xlabel('Reference nodes', fontsize=7)
    ax.set_ylabel('Unified nodes', fontsize=7)
    ax.tick_params(labelsize=7)
# axes[-1].set_visible(False)
fig.suptitle('4.1.2 Search Tree Fidelity: Nodes Explored (unified vs reference)\n'
             'Dashed line = perfect agreement', fontsize=10)
plt.tight_layout()
plt.savefig('evaluation/results/lv_nodes_scatter.png', bbox_inches='tight', dpi=150)
plt.show()

# ── 4.1.3 Runtime Fidelity ────────────────────────────────────────────────────
print('\n=== 4.1.3 Runtime Fidelity ===')
print('Checking ordering is preserved and no algorithm is catastrophically slower.\n')

time_data = all_data[
    (all_data['match_adj'] == 'PASS') &
    (all_data['unified_time'].notna()) &
    (all_data['ref_time'].notna()) &
    (all_data['ref_time'] > 0) &
    (all_data['unified_aborted'] == 0)
].copy()
time_data['speedup'] = time_data['unified_time'] / time_data['ref_time']

speedup = time_data.groupby('algo')['speedup'].agg(['median', 'mean']).round(3)
speedup.index = [LABELS[a] for a in speedup.index]
speedup.columns = ['Median ratio', 'Mean ratio']
print('Runtime ratio (unified / reference). >1.0 means unified is slower:')
print(speedup)

fig, axes = plt.subplots(2, 3, figsize=(14, 7))
axes = axes.flatten()
for i, algo in enumerate(ALGOS):
    ax = axes[i]
    sub = time_data[time_data['algo'] == algo]
    if sub.empty:
        ax.set_visible(False)
        continue
    ax.scatter(sub['ref_time'], sub['unified_time'],
               alpha=0.25, s=5, color=COLORS[algo], rasterized=True)
    lim = max(sub['ref_time'].max(), sub['unified_time'].max()) * 1.1
    ax.plot([1e-3, lim], [1e-3, lim], 'k--', linewidth=0.8, alpha=0.5)
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_title(LABELS[algo], fontsize=9)
    ax.set_xlabel('Reference time (s)', fontsize=7)
    ax.set_ylabel('Unified time (s)', fontsize=7)
    ax.tick_params(labelsize=7)
# axes[-1].set_visible(False)
fig.suptitle('4.1.3 Runtime Fidelity: Wall-Clock Time (unified vs reference)\n'
             'Points above diagonal = unified slower', fontsize=10)
plt.tight_layout()
plt.savefig('evaluation/results/lv_time_scatter.png', bbox_inches='tight', dpi=150)
plt.show()

mcsplit_times = time_data[time_data['algo'] == 'mcsplit'].copy()
mcsplit_times = mcsplit_times.sort_values('speedup', ascending=False)
print(mcsplit_times[['instance_a', 'instance_b', 'unified_time', 'ref_time', 'speedup']].head(20))

print(all_data[(all_data['instance_a']=='g32') | (all_data['instance_b']=='g32')][['instance_a','instance_b','algo','unified_size','ref_size','unified_time','ref_time']].to_string())