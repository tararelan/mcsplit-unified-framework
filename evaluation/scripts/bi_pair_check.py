"""
Paired undirected-vs-directed fidelity comparison.

Rather than a synthetic/instrumented test, this uses the same instance pairs,
run under both undirected and directed search, compared directly. This avoids
any confound from the two datasets covering different instance subsets
(directed runs may have more timeouts, fewer completed pairs, etc.) and gives
a real-data answer to "does fidelity actually drop on directed instances
specifically" rather than a synthetic proxy for the mechanism.

Also computes a weak indirect proxy for branching intensity using columns
you already collect (cut_branches / nodes), since this doesn't require any
new instrumentation - it's already in every row of your existing CSVs.

Run this after validation_analysis_fixed.py has loaded both datasets, or
adjust DIR_UNDIRECTED / DIR_DIRECTED below to point at your two results
folders directly.
"""

import pandas as pd
import numpy as np
from pathlib import Path

ALGOS = ['mcsplit', 'rl', 'll', 'dsb']
LABELS = {
    'mcsplit': 'McSplit', 'rl': 'McSplit+RL', 'll': 'McSplit+LL',
    'dal': 'McSplit+DAL', 'dsb': 'McSplit+DSB',
    'rrsplit': 'RRSplit', 'symsplit': 'SymSplit',
}

DIR_UNDIRECTED = Path('hpc/cross-dataset-evaluation/results/old')
DIR_DIRECTED = Path('hpc/cross-dataset-evaluation/results')

COLS = [
    'instance_a', 'instance_b', 'algo',
    'unified_size', 'unified_edges', 'unified_nodes', 'unified_time',
    'unified_aborted', 'unified_nodes_to_best',
    'unified_time_to_best', 'unified_cut_branches', 'unified_bound_pruned',
    'unified_sym_pruned',
    'ref_size', 'ref_nodes', 'ref_time', 'match'
]

def load(data_dir):
    dfs = []
    for algo in ALGOS:
        path = data_dir / f'{algo}_bi_all.csv'
        if not path.exists():
            print(f'WARNING: {path} not found, skipping {algo}')
            continue
        df = pd.read_csv(path, skipinitialspace=True)
        df = df.reindex(columns=COLS)
        dfs.append(df)
    return pd.concat(dfs, ignore_index=True)

undirected = load(DIR_UNDIRECTED)
directed = load(DIR_DIRECTED)

for df, col_set in [(undirected, ['unified_size', 'unified_aborted', 'unified_time',
                                   'unified_nodes', 'unified_cut_branches',
                                   'ref_size', 'ref_time', 'ref_nodes']),
                     (directed, ['unified_size', 'unified_aborted', 'unified_time',
                                  'unified_nodes', 'unified_cut_branches',
                                  'ref_size', 'ref_time', 'ref_nodes'])]:
    for c in col_set:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors='coerce')

# Only compare pairs where BOTH settings completed successfully on BOTH sides
def clean(df):
    return df[
        (df['unified_aborted'] == 0) &
        (df['unified_nodes'].notna()) & (df['unified_nodes'] > 0) &
        (df['ref_nodes'].notna()) & (df['ref_nodes'] > 0) &
        (df['match'] == 'PASS')
    ][['instance_a', 'instance_b', 'algo', 'unified_nodes', 'ref_nodes', 'unified_cut_branches']]

u_clean = clean(undirected).rename(columns={
    'unified_nodes': 'u_nodes_undirected', 'ref_nodes': 'ref_nodes_undirected',
    'unified_cut_branches': 'u_cut_undirected'})
d_clean = clean(directed).rename(columns={
    'unified_nodes': 'u_nodes_directed', 'ref_nodes': 'ref_nodes_directed',
    'unified_cut_branches': 'u_cut_directed'})

paired = pd.merge(u_clean, d_clean, on=['instance_a', 'instance_b', 'algo'], how='inner')

print(f'Instance pairs completed under BOTH undirected and directed search: {len(paired)}\n')

print('=== Paired fidelity comparison (same instances, both settings) ===\n')
rows = []
for algo in ALGOS:
    sub = paired[paired['algo'] == algo]
    if len(sub) < 5:
        print(f'{LABELS.get(algo, algo)}: insufficient paired data (n={len(sub)}), skipping')
        continue

    r_und = np.log10(sub['u_nodes_undirected']).corr(np.log10(sub['ref_nodes_undirected']))
    r_dir = np.log10(sub['u_nodes_directed']).corr(np.log10(sub['ref_nodes_directed']))

    ratio_und = sub[['u_nodes_undirected', 'ref_nodes_undirected']].max(axis=1) / \
                sub[['u_nodes_undirected', 'ref_nodes_undirected']].min(axis=1)
    ratio_dir = sub[['u_nodes_directed', 'ref_nodes_directed']].max(axis=1) / \
                sub[['u_nodes_directed', 'ref_nodes_directed']].min(axis=1)

    within2x_und = round(100 * (ratio_und <= 2).mean(), 1)
    within2x_dir = round(100 * (ratio_dir <= 2).mean(), 1)

    # Weak indirect proxy for branching intensity: cut_branches per node.
    # Higher = more pruning decisions relative to nodes explored. This is NOT
    # a direct bidomain count, just something already in your existing data
    # that might move in the direction the branching-factor mechanism predicts.
    cut_ratio_und = (sub['u_cut_undirected'] / sub['u_nodes_undirected']).median()
    cut_ratio_dir = (sub['u_cut_directed'] / sub['u_nodes_directed']).median()

    print(f'{LABELS.get(algo, algo)} (n={len(sub)} paired instances):')
    print(f'  r_log10:        undirected={r_und:.4f}  directed={r_dir:.4f}  (delta={r_dir - r_und:+.4f})')
    print(f'  within_2x:      undirected={within2x_und}%  directed={within2x_dir}%  (delta={within2x_dir - within2x_und:+.1f}pp)')
    print(f'  median cut/node: undirected={cut_ratio_und:.3f}  directed={cut_ratio_dir:.3f}  (delta={cut_ratio_dir - cut_ratio_und:+.3f})')
    print()

    rows.append({
        'Algorithm': LABELS.get(algo, algo), 'n_paired': len(sub),
        'r_log10_undirected': round(r_und, 4), 'r_log10_directed': round(r_dir, 4),
        'delta_r': round(r_dir - r_und, 4),
        'within_2x_undirected': within2x_und, 'within_2x_directed': within2x_dir,
        'delta_within_2x_pp': round(within2x_dir - within2x_und, 1),
        'cut_per_node_undirected': round(cut_ratio_und, 3), 'cut_per_node_directed': round(cut_ratio_dir, 3),
    })

summary = pd.DataFrame(rows).set_index('Algorithm')
print('=== Summary table ===')
print(summary)
summary.to_csv('evaluation/results/bi_paired_directed_vs_undirected.csv')
print('\nSaved: evaluation/results/bi_paired_directed_vs_undirected.csv')