import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

plt.rcParams.update({
    'font.size': 14,
    'axes.titlesize': 18,
    'figure.titlesize': 24,
    'axes.labelsize': 16,
    'xtick.labelsize': 13,
    'ytick.labelsize': 13,
})

pd.set_option('display.width', 200)

ALGOS = ['mcsplit', 'rl', 'll', 'dsb', 'rrsplit', 'symsplit']
DATA_DIR = Path('hpc/cross-dataset-evaluation/results')

LABELS = {
    'mcsplit': 'McSplit', 'rl': 'McSplit+RL', 'll': 'McSplit+LL',
    'dsb': 'McSplit+DSB', 'rrsplit': 'RRSplit', 'symsplit': 'SymSplit',
}

COLS = [
    'instance_a', 'instance_b', 'algo',
    'unified_size', 'unified_edges', 'unified_nodes', 'unified_time',
    'unified_aborted', 'unified_nodes_to_best', 'unified_time_to_best',
    'unified_cut_branches', 'unified_bound_pruned', 'unified_sym_pruned',
    'ref_size', 'ref_nodes', 'ref_time', 'match'
]

INF_VAL = 1000

dfs = []
for algo in ALGOS:
    path = DATA_DIR / f'{algo}_lv_all.csv'
    if not path.exists():
        print(f'WARNING: {path} not found, skipping {algo}')
        continue
    df = pd.read_csv(path, skipinitialspace=True)
    df = df.reindex(columns=COLS)
    dfs.append(df)
    print(f'{algo}: {len(df)} rows')

all_data = pd.concat(dfs, ignore_index=True)
for col in ['unified_size', 'unified_aborted', 'unified_time', 'unified_nodes',
            'ref_size', 'ref_time', 'ref_nodes']:
    all_data[col] = pd.to_numeric(all_data[col], errors='coerce')

def get_metric(algo, metric, data):
    sub = data[data['algo'] == algo].set_index(['instance_a', 'instance_b'])
    return sub[metric]

def plot_grid(metric, ylabel, out_path):
    n = len(ALGOS)
    fig, axes = plt.subplots(n, n, figsize=(24, 24), constrained_layout=True)
    for i, algo_row in enumerate(ALGOS):
        for j, algo_col in enumerate(ALGOS):
            ax = axes[i, j]
            if j <= i:
                ax.set_visible(False)
                continue
            m_row = get_metric(algo_row, metric, all_data)
            m_col = get_metric(algo_col, metric, all_data)
            shared = m_row.index.intersection(m_col.index)
            if len(shared) == 0:
                ax.set_visible(False)
                continue
            x = m_col.loc[shared].fillna(INF_VAL).clip(upper=INF_VAL)
            y = m_row.loc[shared].fillna(INF_VAL).clip(upper=INF_VAL)
            ax.scatter(x, y, alpha=0.2, s=6, rasterized=True)
            ax.plot([1, INF_VAL], [1, INF_VAL], 'k--', linewidth=0.8, alpha=0.5)
            ax.set_xscale('log'); ax.set_yscale('log')
            ax.set_title(f'{LABELS[algo_col]} vs {LABELS[algo_row]}')
            ax.set_xlabel(f'{ylabel} ({LABELS[algo_col]})')
            ax.set_ylabel(f'{ylabel} ({LABELS[algo_row]})')
            ax.tick_params(labelsize=7)
    fig.suptitle(f'Pairwise {ylabel} Comparison (LV)')
    # plt.tight_layout()
    plt.savefig(out_path, dpi=500)
    # plt.show()
    print(f'Saved to {out_path}')

plot_grid('ref_nodes', 'Nodes', 'evaluation/results/lv_pairwise_nodes_grid_ref.png')
plot_grid('ref_time', 'Time', 'evaluation/results/lv_pairwise_time_grid_ref.png')