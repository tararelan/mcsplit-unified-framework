import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

BASE_ALGOS = ['mcsplit', 'rl', 'll', 'dal', 'dsb', 'rrsplit', 'symsplit']

ALL_ALGOS = [
    'mcsplit',
    'rl',
    'll', 'll_lum', 'll_lsm',
    'dal', 'dal_rl', 'dal_dal',
    'dsb', 'dsb_never', 'dsb_always',
    'rrsplit', 'rrsplit_nobound', 'rrsplit_nomax', 'rrsplit_noveq',
    'symsplit', 'symsplit_varonly', 'symsplit_valonly',
]

GROUPS = {
    'McSplit': ['mcsplit'],
    'McSplit+RL':         ['rl'],
    'McSplit+LL':         ['ll', 'll_lsm', 'll_lum'],
    'McSplit+DAL':        ['dal', 'dal_rl', 'dal_dal'],
    'McSplit+DSB':        ['dsb', 'dsb_never', 'dsb_always'],
    'RRSplit':            ['rrsplit', 'rrsplit_nobound', 'rrsplit_nomax', 'rrsplit_noveq'],
    'SymSplit':           ['symsplit', 'symsplit_varonly', 'symsplit_valonly'],
}

LABELS = {
    'mcsplit': 'McSplit', 'rl': 'McSplit+RL',
    'll': 'McSplit+LL', 'll_lsm': 'LSM only', 'll_lum': 'LUM only',
    'dal': 'McSplit+DAL', 'dal_rl': 'RL policy only', 'dal_dal': 'DAL policy only',
    'dsb': 'McSplit+DSB', 'dsb_never': 'Never applied', 'dsb_always': 'Always applied',
    'rrsplit': 'RRSplit', 'rrsplit_nobound': 'No bound',
    'rrsplit_nomax': 'No maximality', 'rrsplit_noveq': 'No vertex-equivalence',
    'symsplit': 'SymSplit', 'symsplit_varonly': 'Variable symmetry only', 'symsplit_valonly': 'Value symmetry only',
}

COLOURS = {
    'McSplit': '#555555', 'McSplit+RL': "#ff8e43",
    'McSplit+LL': "#e0a800ce", 'McSplit+DAL': '#c94040',
    'McSplit+DSB': '#7a5c9e', 'RRSplit': '#2e86ab', 'SymSplit': '#3a9e5f',
}

TIMEOUT = 1000.0
EASY_THRESHOLD = 10.0
BORDERLINE_THRESHOLD = 850.0

inst_key = ['instance_a', 'instance_b']

def load_files(data_directory, cols, algos=ALL_ALGOS):
    dfs = []
    for algo in algos:
        path = data_directory / f'{algo}_all.csv'
        if not path.exists():
            print(f'WARNING: {path} not found, skipping {algo}')
            continue
        df = pd.read_csv(path, header=None, names=cols)

        dfs.append(df)

    all_data = pd.concat(dfs, ignore_index=True)
    
    return all_data

def reclassify_instances(row):
    if row['match'] == 'FAIL' and row['unified_aborted'] == 1:
        if pd.notna(row['ref_time']) and row['ref_time'] > BORDERLINE_THRESHOLD:
            return 'MISSING_DATA'
    return row['match']

def cactus_plot():
    pass

def scatter_plot(data, x, y, x_label, y_label, filename, title, algos=BASE_ALGOS):
    fig, axes = plt.subplots(2, 4, figsize=(14, 7))
    axes = axes.flatten()
    for i, algo in enumerate(algos):
        ax = axes[i]
        sub = data[data['algo'] == algo]
        algo_name = LABELS[algo]
        if sub.empty:
            ax.set_visible(False)
            continue
        ax.scatter(x, y,
                alpha=0.25, s=5, color=COLOURS[algo_name], rasterized=True)
        lim = max(x.max(), y.max()) * 1.1
        ax.plot([1, lim], [1, lim], 'k--', linewidth=0.8, alpha=0.5)
        ax.set_xscale('log'); ax.set_yscale('log')
        ax.set_title(LABELS[algo], fontsize=9)
        ax.set_xlabel(x_label, fontsize=7)
        ax.set_ylabel(y_label, fontsize=7)
        ax.tick_params(labelsize=7)
    axes[-1].set_visible(False)
    fig.suptitle(title, fontsize=10)
    plt.tight_layout()
    plt.savefig(f'evaluation/results/{filename}.png', bbox_inches='tight', dpi=150)
    plt.show()

def classify_instances(idx, data):
    data['solve_time'] = data.apply(
        lambda r: r['unified_time'] if r['unified_aborted'] == 0 else np.nan, axis=1
    )

    time_pivot = data.pivot_table(
        index=inst_key, columns='algo', values='solve_time', aggfunc='first'
    )
    abort_pivot = data.pivot_table(
        index=inst_key, columns='algo', values='unified_aborted', aggfunc='first'
    )
    abort_any = (abort_pivot.fillna(1).astype(int) == 1).astype(int)

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

def medium_instances_summary():
    pass

def easy_instances_summary():
    pass

def hard_instances_summary():
    pass

def solved_instances():
    pass

def symmetry_pruned_summary():
    pass

def overhead_analysis():
    pass