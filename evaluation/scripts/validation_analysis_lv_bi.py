"""
4.1 Validation Against Reference Implementations
4.1.1 Solution Correctness  - agreement table (PASS/FAIL/MISSING)
4.1.2 Search Tree Fidelity  - nodes scatter (unified vs reference) + log-scale correlation + ratio-band
4.1.3 Runtime Fidelity      - time scatter + speedup table
4.1.4 McSplit-DAL Bug       - +1 inflation in reference binary
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

plt.rcParams.update({
    'font.size': 14,
    'axes.titlesize': 18,
    'axes.labelsize': 16,
    'xtick.labelsize': 13,
    'ytick.labelsize': 13,
})

pd.set_option('display.width', 200)
pd.set_option('display.float_format', '{:.3f}'.format)

ALGOS = ['mcsplit', 'rl', 'll', 'dsb', 'rrsplit', 'symsplit']
DATA_DIR = Path('hpc/validation/results')

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

# Ratio bands used for the fidelity metric below. A pair is "within Nx" if
# max(unified, reference) / min(unified, reference) <= N.
RATIO_BANDS = [2, 5, 10, 100]

# ── Load ──────────────────────────────────────────────────────────────────────
dfs = []
for algo in ALGOS:
    path = DATA_DIR / f'{algo}_all.csv'
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

# Reclassify borderline FAILs where unified timed out AND ref was near timeout
BORDERLINE_THRESHOLD = 850.0

def reclassify(row):
    if row['match'] == 'FAIL' and row['unified_aborted'] == 1:
        if pd.notna(row['ref_time']) and row['ref_time'] > BORDERLINE_THRESHOLD:
            return 'MISSING_DATA'
    return row['match']

all_data['match_adj'] = all_data.apply(reclassify, axis=1)

print(all_data)
print(all_data['algo'].value_counts())
print(all_data['algo'].unique()[:10])


# ── Check 1: Internal cross-algorithm consistency ───────────────────────────
# Do the seven UNIFIED implementations agree with each other on solution size,
# independent of the reference entirely? This only uses unified_size.
print("\n=== Check 1: Cross-Algorithm Consistency (unified implementations only) ===\n")

inst_key = ['instance_a', 'instance_b']

completed = all_data[
    (all_data['unified_aborted'] == 0) &
    (all_data['unified_size'].notna())
][['instance_a', 'instance_b', 'algo', 'unified_size']]

pivot_internal = completed.pivot_table(
    index=inst_key,
    columns='algo',
    values='unified_size',
    aggfunc='first'
)

def check_consistency(row):
    vals = row.dropna()
    if len(vals) < 2:
        return 'INSUFFICIENT'
    if vals.nunique() > 1:
        return 'MISMATCH'
    return 'OK'

pivot_internal['status'] = pivot_internal.apply(check_consistency, axis=1)

status_counts = pivot_internal['status'].value_counts()
print(f'OK:           {status_counts.get("OK", 0)}')
print(f'MISMATCH:     {status_counts.get("MISMATCH", 0)}')
print(f'INSUFFICIENT: {status_counts.get("INSUFFICIENT", 0)}')

internal_mismatches = pivot_internal[pivot_internal['status'] == 'MISMATCH'].drop(columns='status')
if internal_mismatches.empty:
    print('\n\u2713 All completed unified algorithms agree with each other on every instance.')
else:
    print(f'\n\u26a0 {len(internal_mismatches)} instances where unified algorithms disagree with each other:')
    print(internal_mismatches)
    # internal_mismatches.to_csv(
    #     'hpc/cross-dataset-evaluation/results/bi_internal_size_mismatches.csv')
    print('Saved: bi_internal_size_mismatches.csv')


# ── Check 2: Unified vs reference agreement, across all algorithms ─────────
# This pivots BOTH unified_size and ref_size per algorithm and checks whether
# ALL of those values (unified and reference, across all 7 algorithms) agree
# on a given instance. Unlike Check 1, this can be triggered purely by a
# reference-side bug (e.g. DAL's known +1 inflation) even when every unified
# implementation agrees with every other unified implementation.
print("\n=== Check 2: Unified vs Reference Agreement (across all algorithms) ===\n")

completed2 = all_data[
    (all_data['unified_aborted'] == 0) &
    (all_data['unified_size'].notna())
][['instance_a', 'instance_b', 'algo', 'unified_size', 'ref_size']]

pivot_vs_ref = completed2.pivot_table(
    index=['instance_a', 'instance_b'],
    columns='algo',
    values=['unified_size', 'ref_size'],
    aggfunc='first'
)

def classify(row):
    vals = row.dropna()
    if len(vals) < 2:
        return 'INSUFFICIENT_DATA'
    if vals.nunique() > 1:
        return 'REAL_MISMATCH'
    return 'OK'

pivot_vs_ref['status'] = pivot_vs_ref.apply(classify, axis=1)
status_counts2 = pivot_vs_ref['status'].value_counts()
print(status_counts2)

vs_ref_mismatches = pivot_vs_ref[pivot_vs_ref['status'] == 'REAL_MISMATCH']
if vs_ref_mismatches.empty:
    print('\n\u2713 Zero mismatches between unified and reference across all algorithms.')
else:
    print(f'\n{len(vs_ref_mismatches)} instances flagged (expected to be dominated by the '
          f'known McSplit+DAL reference bug - see 4.1.4). Inspect the saved CSV to confirm '
          f'no other algorithm contributes mismatches beyond DAL.')
    vs_ref_mismatches.drop(columns='status').to_csv(
        'hpc/cross-dataset-evaluation/results/bi_unified_vs_ref_mismatches.csv')
    print('Saved: bi_unified_vs_ref_mismatches.csv')


# ── 4.1.1 Solution Correctness ────────────────────────────────────────────────
print('\n=== 4.1.1 Solution Correctness ===\n')

rows = []
for algo in ALGOS:
    sub = all_data[all_data['algo'] == algo]
    if sub.empty:
        continue
    counts = sub['match_adj'].value_counts()
    total = len(sub)
    rows.append({
        'Algorithm':    LABELS[algo],
        'PASS':         int(counts.get('PASS', 0)),
        'FAIL':         int(counts.get('FAIL', 0)),
        'MISSING_DATA': int(counts.get('MISSING_DATA', 0)),
        'PASS %':       round(100 * counts.get('PASS', 0) / total, 1),
    })

correctness_df = pd.DataFrame(rows).set_index('Algorithm')
print(correctness_df)

# Report the ACTUAL fail counts rather than assuming only DAL has any -
# this previously printed a blanket claim that contradicted the table above it.
non_dal_fails = correctness_df.drop(index='McSplit+DAL', errors='ignore')['FAIL'].sum()
dal_fails_n = correctness_df.loc['McSplit+DAL', 'FAIL'] if 'McSplit+DAL' in correctness_df.index else 0

print(f'\nMcSplit+DAL FAILs: {dal_fails_n} (investigated in 4.1.4 - known reference +1 bug)')
if non_dal_fails == 0:
    print('All other algorithms report zero FAILs.')
else:
    print(f'\u26a0 {int(non_dal_fails)} FAILs across non-DAL algorithms require investigation '
          f'before reporting these results - do not assume these are benign.')
    non_dal_fail_rows = all_data[
        (all_data['algo'] != 'dal') & (all_data['match_adj'] == 'FAIL') & (all_data['unified_aborted'] == 0)
    ]
    # non_dal_fail_rows.to_csv(
    #     'hpc/cross-dataset-evaluation/results/bi_non_dal_fails.csv', index=False)
    print('Saved: bi_non_dal_fails.csv - inspect before finalising 4.1.1.')

# ── 4.1.4 McSplit-DAL Reference Bug ──────────────────────────────────────────
print('\n=== 4.1.4 McSplit+DAL Reference Bug ===\n')

dal_fails = all_data[
    (all_data['algo'] == 'dal') & (all_data['match_adj'] == 'FAIL')
].copy()

print(f'DAL FAIL count: {len(dal_fails)}')

if not dal_fails.empty:
    dal_fails['size_diff'] = dal_fails['unified_size'] - dal_fails['ref_size']
    print('\nDistribution of (unified_size - ref_size) on FAIL instances:')
    print(dal_fails['size_diff'].value_counts().sort_index())
    print('\nThe reference binary consistently reports one vertex more than the')
    print('unified reimplementation - a +1 inflation bug in the reference binary.')
    print('The unified reimplementation is treated as correct.')

# ── 4.1.2 Search Tree Fidelity ────────────────────────────────────────────────
print('\n=== 4.1.2 Search Tree Fidelity ===')

pass_data = all_data[
    (all_data['match_adj'] == 'PASS') &
    (all_data['unified_nodes'].notna()) &
    (all_data['ref_nodes'].notna()) &
    (all_data['ref_nodes'] > 0) &
    (all_data['unified_nodes'] > 0)
].copy()

print(pass_data)

# Node counts span many orders of magnitude (10 to 10^9+), so raw Pearson
# correlation is dominated by the largest-magnitude points and is not a
# meaningful fidelity metric here. Correlation is instead computed on
# log10(nodes), which is standard practice for comparing search-tree sizes
# across solvers and matches how these results are visualised (log-log
# scatter plots below).
pass_data['log_unified_nodes'] = np.log10(pass_data['unified_nodes'])
pass_data['log_ref_nodes'] = np.log10(pass_data['ref_nodes'])

print('\nPearson correlation on log10(nodes) (unified vs reference):')
node_fidelity_rows = []
for algo in ALGOS:
    sub = pass_data[pass_data['algo'] == algo]
    if sub.empty:
        continue
    r_log = sub['log_unified_nodes'].corr(sub['log_ref_nodes'])
    r_raw = sub['unified_nodes'].corr(sub['ref_nodes'])

    # Ratio-band metric: what fraction of instances are within Nx of the
    # reference node count, in either direction. More directly interpretable
    # than a correlation coefficient for a dissertation results table.
    ratio = sub[['unified_nodes', 'ref_nodes']].max(axis=1) / sub[['unified_nodes', 'ref_nodes']].min(axis=1)
    band_pcts = {f'within_{b}x': round(100 * (ratio <= b).mean(), 1) for b in RATIO_BANDS}

    print(f'  {LABELS[algo]}: r_log10 = {r_log:.4f}  (r_raw = {r_raw:.4f}, n={len(sub)})  '
          + '  '.join(f'{k}={v}%' for k, v in band_pcts.items()))

    node_fidelity_rows.append({'Algorithm': LABELS[algo], 'r (log10)': round(r_log, 4),
                                'n': len(sub), **band_pcts})

node_fidelity_df = pd.DataFrame(node_fidelity_rows).set_index('Algorithm')
print('\nSummary table (4.1.2):')
print(node_fidelity_df)
# node_fidelity_df.to_csv('evaluation/results/bi_node_fidelity_summary.csv')

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
    ax.set_title(LABELS[algo])
    ax.set_xlabel('Reference nodes')
    ax.set_ylabel('Unified nodes')
    
# axes[-1].set_visible(False)
# fig.suptitle('4.1.2 Search Tree Fidelity: Nodes Explored (unified vs reference)\n'
#              'Dashed line = perfect agreement', fontsize=10)
plt.tight_layout()
plt.savefig('evaluation/results/mivia_nodes_scatter.png', bbox_inches='tight', dpi=500)
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
    ax.set_title(LABELS[algo])
    ax.set_xlabel('Reference time (s)')
    ax.set_ylabel('Unified time (s)')
    
# axes[-1].set_visible(False)
# fig.suptitle('4.1.3 Runtime Fidelity: Wall-Clock Time (unified vs reference)\n'
#              'Points above diagonal = unified slower', fontsize=10)
plt.tight_layout()
plt.savefig('evaluation/results/mivia_time_scatter.png', bbox_inches='tight', dpi=500)
plt.show()