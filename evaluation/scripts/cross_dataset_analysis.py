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
    'unified_aborted', 'unified_nodes_to_best',
    'unified_time_to_best', 'unified_cut_branches', 'unified_bound_pruned',
    'unified_sym_pruned'
]

TIMEOUT = 1000.0
EASY_THRESHOLD = 10.0
inst_key = ['instance_a', 'instance_b']

# ── 2. Load helper ────────────────────────────────────────────────────────────
def load_dataset(tag):
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
    # ax.set_xscale('log')
    ax.set_title(title, fontsize=11)
    ax.legend(fontsize=9, loc='upper left')
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(fname, bbox_inches='tight')
    plt.show()
    print(f'Saved {fname}')


def cactus_nodes(medium_data, title, fname):
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
    ax.set_ylabel('Instances Solved', fontsize=11)
    # ax.set_xscale('log')
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
    hard_instances = set(difficulty[difficulty == 'hard'].index)
    hard_mask = data.set_index(inst_key).index.isin(hard_instances)
    hard_data = data[hard_mask].copy()

    print(f'\n--- Hard Instance Analysis: {dataset_name} ---')
    print(f'Total hard instances: {len(hard_instances)}\n')

    # Per-algorithm hard instance counts and incumbent sizes
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
            'Max incumbent':   int(sizes.max()) if len(sizes) > 0 else np.nan,
        })

    hard_df = pd.DataFrame(rows).set_index('Algorithm')
    print(hard_df)

    # Pairwise dominance matrix
    pivot = hard_data.pivot_table(
        index=inst_key, columns='algo', values='unified_size', aggfunc='first'
    )

    print('\nPairwise: fraction of hard instances where algo A finds strictly larger incumbent than algo B:')
    print('(Only instances where both have a valid incumbent)\n')

    header = [''] + [LABELS[a] for a in ALGOS if a in pivot.columns]
    rows_dom = []
    for a in ALGOS:
        if a not in pivot.columns:
            continue
        row = [LABELS[a]]
        for b in ALGOS:
            if b not in pivot.columns:
                continue
            if a == b:
                row.append('-')
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

    # SymSplit vs RRSplit on non-equal instances (replicates SymSplit paper Table 2)
    if 'symsplit' in pivot.columns and 'rrsplit' in pivot.columns:
        disagree = pivot['symsplit'] != pivot['rrsplit']
        both_valid = pivot['symsplit'].notna() & pivot['rrsplit'].notna()
        mask = disagree & both_valid
        if mask.sum() > 0:
            sym_wins = (pivot.loc[mask, 'symsplit'] > pivot.loc[mask, 'rrsplit']).mean()
            print(f'\nSymSplit > RRSplit among disagreeing instances: {sym_wins:.2%} '
                  f'(n={mask.sum()}) — cf. SymSplit paper Table 2 (>84%)')

def sym_pruned_summary(medium_data, dataset_name):
    """Symmetry-pruned branches for RRSplit/SymSplit - shows if mechanisms fire."""
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
                'BI Dataset (Biochemical Reactions) Time versus Instances',
                'evaluation/results/bi_cactus_time.png')
    
    cactus_nodes(bi_medium_data,
                'BI Dataset (Biochemical Reactions) Nodes versus Instances',
                'evaluation/results/bi_cactus_nodes.png')
    bi_incumbent = hard_instance_analysis(bi_data, bi_difficulty, 'BI')
    bi_sym = sym_pruned_summary(bi_medium_data, 'BI')
else:
    print('No BI data loaded - check file paths.')

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
                'LV Dataset Time versus Instances',
                'evaluation/results/lv_cactus_time.png')
    
    cactus_nodes(lv_medium_data,
                'LV Dataset Nodes versus Instances',
                'evaluation/results/lv_cactus_nodes.png')

    lv_incumbent = hard_instance_analysis(lv_data, lv_difficulty, 'LV')
    lv_sym = sym_pruned_summary(lv_medium_data, 'LV')
else:
    print('No LV data loaded - check file paths.')

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


# Compare nodes_to_best and time_to_best between SymSplit and RRSplit
# on instances where both completed (medium) or both timed out (hard)

for dataset_name, data, instances in [
    ('BI', bi_data, bi_medium),
    ('LV', lv_data, lv_medium),
]:
    sub = data[data['algo'].isin(['rrsplit', 'symsplit'])].copy()
    
    # Pivot nodes_to_best and time_to_best
    nodes_pivot = sub.pivot_table(
        index=inst_key, columns='algo', 
        values='unified_nodes_to_best', aggfunc='first'
    )
    time_pivot = sub.pivot_table(
        index=inst_key, columns='algo',
        values='unified_time_to_best', aggfunc='first'
    )
    
    # Only instances where both have data
    mask = nodes_pivot['rrsplit'].notna() & nodes_pivot['symsplit'].notna()
    nodes_pivot = nodes_pivot[mask]
    time_pivot = time_pivot[mask]
    
    print(f'\n=== {dataset_name} ===')
    print(f'Instances compared: {mask.sum()}')
    
    print('\nNodes to best incumbent:')
    print(f'  SymSplit median:  {nodes_pivot["symsplit"].median():.0f}')
    print(f'  RRSplit median:   {nodes_pivot["rrsplit"].median():.0f}')
    print(f'  SymSplit finds incumbent in fewer nodes: '
          f'{(nodes_pivot["symsplit"] < nodes_pivot["rrsplit"]).sum()}/{mask.sum()} instances')
    
    print('\nTime to best incumbent:')
    print(f'  SymSplit median:  {time_pivot["symsplit"].median():.3f}s')
    print(f'  RRSplit median:   {time_pivot["rrsplit"].median():.3f}s')
    print(f'  SymSplit finds incumbent faster: '
          f'{(time_pivot["symsplit"] < time_pivot["rrsplit"]).sum()}/{mask.sum()} instances')
    
    # Total nodes explored comparison on instances both solved
    total_nodes = sub.pivot_table(
        index=inst_key, columns='algo',
        values='unified_nodes', aggfunc='first'
    )
    
    # Both solved (not aborted)
    abort_pivot = sub.pivot_table(
        index=inst_key, columns='algo',
        values='unified_aborted', aggfunc='first'
    )
    both_solved = (abort_pivot['rrsplit'] == 0) & (abort_pivot['symsplit'] == 0)
    total_nodes_solved = total_nodes[both_solved]
    
    print(f'\nTotal nodes explored (instances both solved, n={both_solved.sum()}):')
    print(f'  SymSplit median:  {total_nodes_solved["symsplit"].median():.0f}')
    print(f'  RRSplit median:   {total_nodes_solved["rrsplit"].median():.0f}')
    print(f'  SymSplit explores fewer total nodes: '
          f'{(total_nodes_solved["symsplit"] < total_nodes_solved["rrsplit"]).sum()}/{both_solved.sum()} instances')
    
    # Also check time per node to confirm the per-node cost explanation
    time_solved = sub[sub['unified_aborted'] == 0].pivot_table(
        index=inst_key, columns='algo',
        values='unified_time', aggfunc='first'
    )
    nodes_solved = total_nodes_solved
    common = time_solved.index.intersection(nodes_solved.index)
    
    sym_tpn = (time_solved.loc[common, 'symsplit'] / nodes_solved.loc[common, 'symsplit'] * 1e6).median()
    rr_tpn  = (time_solved.loc[common, 'rrsplit']  / nodes_solved.loc[common, 'rrsplit']  * 1e6).median()
    print(f'\nMedian time per node on solved instances (µs):')
    print(f'  SymSplit: {sym_tpn:.3f}')
    print(f'  RRSplit:  {rr_tpn:.3f}')


for dataset_name, data, medium_instances in [
    ('BI', bi_data, bi_medium),
    ('LV', lv_data, lv_medium),
]:
    med_mask = data.set_index(inst_key).index.isin(medium_instances)
    medium_data_ds = data[med_mask].copy()
    completed = medium_data_ds[medium_data_ds['unified_aborted'] == 0].copy()
    completed['time_per_node_us'] = (
        completed['unified_time'] / completed['unified_nodes'] * 1e6
    )
    overhead = completed.groupby('algo')['time_per_node_us'].agg(['median', 'mean'])
    overhead.index = [LABELS.get(a, a) for a in overhead.index]
    overhead.columns = ['Median µs/node', 'Mean µs/node']
    print(f'\n=== Per-Node Overhead: {dataset_name} (medium instances) ===')
    print(overhead.round(3))