"""
4.4 Cross-Dataset Evaluation (BI and LV)
Produces:
  - Difficulty distributions
  - Solved counts table
  - Cactus plots (time and nodes)
  - Hard instance incumbent analysis (delta distribution)
  - Symmetry-pruned branches comparison across datasets
"""

# ── 1. Setup ──────────────────────────────────────────────────────────────────
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

pd.set_option('display.width', 200)
pd.set_option('display.float_format', '{:.3f}'.format)

ALGOS = ['mcsplit', 'rl', 'll', 'dal', 'dsb', 'rrsplit', 'symsplit']
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
    'unified_aborted', 'unified_root_ub', 'unified_nodes_to_best',
    'unified_time_to_best', 'unified_cut_branches', 'unified_bound_pruned',
    'unified_sym_pruned', 'unified_conflicts',
]

TIMEOUT = 1000.0
EASY_THRESHOLD = 10.0
inst_key = ['instance_a', 'instance_b']

# ── 2. Load helper ────────────────────────────────────────────────────────────
def load_dataset(tag):
    dfs = []
    for algo in ALGOS:
        path = DATA_DIR / f'{algo}_{tag}_all.csv'
        if not path.exists():
            print(f'WARNING: {path} not found, skipping {algo}')
            continue
        df = pd.read_csv(path, header=None, names=COLS)
        dfs.append(df)
    if not dfs:
        return pd.DataFrame()
    data = pd.concat(dfs, ignore_index=True)
    for col in ['unified_size', 'unified_aborted', 'unified_time',
                'unified_nodes', 'unified_bound_pruned', 'unified_sym_pruned']:
        data[col] = pd.to_numeric(data[col], errors='coerce')
    return data

# ── 3. Difficulty classification helper ───────────────────────────────────────
def classify_difficulty(data):
    data = data.copy()
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
    all_instances = data[inst_key].drop_duplicates().set_index(inst_key).index

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

    return pd.Series(
        [classify(idx) for idx in all_instances],
        index=all_instances, name='difficulty'
    )

# ── 4. Analysis helper functions ──────────────────────────────────────────────
def solved_table(medium_data, medium_instances):
    total = len(medium_instances)
    rows = []
    for algo in ALGOS:
        sub = medium_data[medium_data['algo'] == algo]
        if sub.empty:
            rows.append({'Algorithm': LABELS[algo], 'Solved': np.nan, 'Solved %': np.nan})
            continue
        solved = (sub['unified_aborted'] == 0).sum()
        rows.append({'Algorithm': LABELS[algo], 'Solved': int(solved),
                     'Solved %': round(100 * solved / total, 1)})
    return pd.DataFrame(rows).set_index('Algorithm')

def cactus_time(medium_data, title, fname):
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
    ax.set_title(title, fontsize=11)
    ax.legend(fontsize=9, loc='upper left')
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(fname, bbox_inches='tight')
    plt.show()
    print(f'Saved {fname}')

def hard_instance_analysis(data, difficulty, dataset_name):
    """For hard instances: what incumbent size did each algo find before timeout?"""
    hard_instances = set(difficulty[difficulty == 'hard'].index)
    hard_mask = data.set_index(inst_key).index.isin(hard_instances)
    hard_data = data[hard_mask].copy()

    print(f'\n--- Hard Instance Analysis: {dataset_name} ---')
    print(f'Total hard instances: {len(hard_instances)}\n')

    # Incumbent sizes found before timeout
    rows = []
    for algo in ALGOS:
        sub = hard_data[hard_data['algo'] == algo]
        if sub.empty:
            continue
        # All hard instances are aborted; unified_size = best incumbent found
        sizes = sub['unified_size'].dropna()
        rows.append({
            'Algorithm': LABELS[algo],
            'Mean incumbent': round(sizes.mean(), 2),
            'Median incumbent': round(sizes.median(), 2),
        })
    incumbent_df = pd.DataFrame(rows).set_index('Algorithm')
    print('Incumbent sizes found on hard instances (before timeout):')
    print(incumbent_df)

    # Delta analysis: for each instance where algos disagree,
    # how much larger is the best algo's incumbent vs others?
    # Pivot incumbent sizes
    pivot = hard_data.pivot_table(
        index=inst_key, columns='algo', values='unified_size', aggfunc='first'
    )

    # Pairwise delta between symsplit (best expected) and each other algo
    if 'symsplit' in pivot.columns:
        print('\nDelta: SymSplit incumbent - other algorithm (on hard instances where they differ):')
        for algo in ALGOS:
            if algo == 'symsplit' or algo not in pivot.columns:
                continue
            mask = pivot['symsplit'].notna() & pivot[algo].notna()
            delta = pivot.loc[mask, 'symsplit'] - pivot.loc[mask, algo]
            nonzero = delta[delta != 0]
            if len(nonzero) == 0:
                print(f'  vs {LABELS[algo]}: no differences')
                continue
            print(f'  vs {LABELS[algo]}: mean={nonzero.mean():.2f}, '
                  f'median={nonzero.median():.1f}, '
                  f'SymSplit better in {(nonzero > 0).sum()}/{len(nonzero)} cases')

    return incumbent_df

def sym_pruned_summary(medium_data, dataset_name):
    """Symmetry-pruned branches for RRSplit/SymSplit — shows if mechanisms fire."""
    sym_algos = ['rrsplit', 'symsplit']
    sub = medium_data[
        medium_data['algo'].isin(sym_algos) & (medium_data['unified_aborted'] == 0)
    ]
    result = sub.groupby('algo')['unified_sym_pruned'].agg(['mean', 'median', 'sum'])
    result.index = [LABELS.get(a, a) for a in result.index]
    result.columns = ['Mean sym-pruned', 'Median sym-pruned', 'Total sym-pruned']
    print(f'\nSymmetry-pruned branches ({dataset_name}, medium instances):')
    print(result.round(1))
    return result

# ─────────────────────────────────────────────────────────────────────────────
# ── 5. BI Dataset ─────────────────────────────────────────────────────────────
# ─────────────────────────────────────────────────────────────────────────────
print('=' * 60)
print('4.4.1 Results on BI (Biochemical Reactions)')
print('=' * 60)

bi_data = load_dataset('bi')

if not bi_data.empty:
    bi_difficulty = classify_difficulty(bi_data)
    print(f'\nBI difficulty distribution:')
    print(bi_difficulty.value_counts())

    bi_medium = set(bi_difficulty[bi_difficulty == 'medium'].index)
    bi_med_mask = bi_data.set_index(inst_key).index.isin(bi_medium)
    bi_medium_data = bi_data[bi_med_mask].copy()

    print(f'\nBI solved counts (medium instances, n={len(bi_medium)}):')
    print(solved_table(bi_medium_data, bi_medium))

    cactus_time(bi_medium_data,
                'Cactus Plot: BI Dataset (Biochemical Reactions)',
                'evaluation/results/bi_cactus_time.png')

    bi_incumbent = hard_instance_analysis(bi_data, bi_difficulty, 'BI')
    bi_sym = sym_pruned_summary(bi_medium_data, 'BI')
else:
    print('No BI data loaded — check file paths.')

# ─────────────────────────────────────────────────────────────────────────────
# ── 6. LV Dataset ─────────────────────────────────────────────────────────────
# ─────────────────────────────────────────────────────────────────────────────
print('\n' + '=' * 60)
print('4.4.2 Results on LV (Biological Networks)')
print('=' * 60)

lv_data = load_dataset('lv')

if not lv_data.empty:
    lv_difficulty = classify_difficulty(lv_data)
    print(f'\nLV difficulty distribution:')
    print(lv_difficulty.value_counts())

    lv_medium = set(lv_difficulty[lv_difficulty == 'medium'].index)
    lv_med_mask = lv_data.set_index(inst_key).index.isin(lv_medium)
    lv_medium_data = lv_data[lv_med_mask].copy()

    print(f'\nLV solved counts (medium instances, n={len(lv_medium)}):')
    print(solved_table(lv_medium_data, lv_medium))

    cactus_time(lv_medium_data,
                'Cactus Plot: LV Dataset (Biological Networks)',
                'evaluation/results/lv_cactus_time.png')

    lv_incumbent = hard_instance_analysis(lv_data, lv_difficulty, 'LV')
    lv_sym = sym_pruned_summary(lv_medium_data, 'LV')
else:
    print('No LV data loaded — check file paths.')

# ─────────────────────────────────────────────────────────────────────────────
# ── 7. Comparison with MIVIA: symmetry-pruned branches across datasets ─────────
# ─────────────────────────────────────────────────────────────────────────────
print('\n' + '=' * 60)
print('4.4.4 Comparison with MIVIA: Symmetry Mechanism Activity')
print('=' * 60)

# Load MIVIA medium data for the sym comparison
# (reuse validation results for base algos)
VALIDATION_DIR = Path('../../hpc/validation/results')
VAL_COLS = COLS + ['ref_size', 'ref_nodes', 'ref_time', 'match']

mivia_dfs = []
for algo in ['rrsplit', 'symsplit']:
    path = VALIDATION_DIR / f'{algo}_all.csv'
    if not path.exists():
        continue
    df = pd.read_csv(path, header=None, names=VAL_COLS)
    df['algo'] = algo
    mivia_dfs.append(df[COLS])

if mivia_dfs:
    mivia_data = pd.concat(mivia_dfs, ignore_index=True)
    for col in ['unified_sym_pruned', 'unified_aborted']:
        mivia_data[col] = pd.to_numeric(mivia_data[col], errors='coerce')

    # Use all completed instances for comparison
    mivia_completed = mivia_data[mivia_data['unified_aborted'] == 0]

    datasets_sym = {}
    for algo in ['rrsplit', 'symsplit']:
        sub_mivia = mivia_completed[mivia_completed['algo'] == algo]['unified_sym_pruned']
        datasets_sym[(algo, 'MIVIA')] = sub_mivia.mean()
        if not bi_data.empty:
            sub_bi = bi_data[
                (bi_data['algo'] == algo) & (bi_data['unified_aborted'] == 0)
            ]['unified_sym_pruned']
            datasets_sym[(algo, 'BI')] = sub_bi.mean()
        if not lv_data.empty:
            sub_lv = lv_data[
                (lv_data['algo'] == algo) & (lv_data['unified_aborted'] == 0)
            ]['unified_sym_pruned']
            datasets_sym[(algo, 'LV')] = sub_lv.mean()

    comparison_df = pd.Series(datasets_sym).unstack(level=1).round(1)
    comparison_df.index = [LABELS.get(a, a) for a in comparison_df.index]
    print('\nMean symmetry-pruned branches per completed instance:')
    print(comparison_df)
    print('\nHigher values = symmetry mechanisms fire more = algorithm is doing more useful work.')

    # Bar chart comparison
    fig, ax = plt.subplots(figsize=(7, 4))
    datasets = [d for d in ['MIVIA', 'BI', 'LV'] if d in comparison_df.columns]
    x = np.arange(len(datasets))
    width = 0.35

    for i, algo in enumerate(['rrsplit', 'symsplit']):
        label = LABELS.get(algo, algo)
        if label in comparison_df.index:
            vals = [comparison_df.loc[label, d] for d in datasets]
            ax.bar(x + i * width, vals, width, label=label,
                   color=COLORS[algo], alpha=0.85)

    ax.set_xticks(x + width / 2)
    ax.set_xticklabels(datasets, fontsize=10)
    ax.set_ylabel('Mean symmetry-pruned branches', fontsize=10)
    ax.set_title('Symmetry Mechanism Activity Across Datasets', fontsize=10)
    ax.legend(fontsize=9)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig('evaluation/results/sym_pruned_comparison.png', bbox_inches='tight')
    plt.show()
    print('Saved sym_pruned_comparison.png')
