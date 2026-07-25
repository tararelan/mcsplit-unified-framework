"""
4.3 Ablation Study
Produces:
  - Dot plot: all ablation variants vs McSplit baseline
  - Symmetry-pruned branches on MIVIA (mechanisms rarely fire)
  - Bound-pruned branches comparison across families
"""

# ── 1. Setup ──────────────────────────────────────────────────────────────────
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from pathlib import Path
from utils import BASE_ALGOS, ALL_ALGOS, GROUPS, LABELS, COLOURS as GROUP_COLORS, TIMEOUT, EASY_THRESHOLD

pd.set_option('display.width', 200)
pd.set_option('display.float_format', '{:.3f}'.format)

VALIDATION_DIR = Path('hpc/validation/results')
ABLATION_DIR   = Path('hpc/ablation/results')

dfs = []

COLS = [
    'instance_a', 'instance_b', 'algo',
    'unified_size', 'unified_edges', 'unified_nodes', 'unified_time',
    'unified_aborted', 'unified_nodes_to_best',
    'unified_time_to_best', 'unified_cut_branches', 'unified_bound_pruned',
    'unified_sym_pruned',
]

# Base algorithms from validation results
for algo in BASE_ALGOS:
    path = VALIDATION_DIR / f'{algo}_all.csv'
    if not path.exists():
        print(f'WARNING: {path} not found, skipping {algo}')
        continue
    df = pd.read_csv(path, header=None, names=COLS + ['ref_size', 'ref_nodes', 'ref_time', 'match'])
    df['algo'] = algo
    dfs.append(df[COLS])

# Ablation variants
ablation_variants = [a for a in ALL_ALGOS if a not in BASE_ALGOS]
for algo in ablation_variants:
    path = ABLATION_DIR / f'{algo}_all.csv'
    if not path.exists():
        print(f'WARNING: {path} not found, skipping {algo}')
        continue
    df = pd.read_csv(path, header=None, skipinitialspace=True)
    df = df.iloc[:, :15]  # take first 15 columns
    df.columns = COLS
    df = df[df['algo'] == algo].copy()  # drop header rows
    dfs.append(df)

all_data = pd.concat(dfs, ignore_index=True)

print(all_data['algo'].unique())

for col in ['unified_size', 'unified_aborted', 'unified_time', 'unified_nodes',
            'unified_bound_pruned', 'unified_sym_pruned', 'unified_cut_branches']:
    all_data[col] = pd.to_numeric(all_data[col], errors='coerce')

inst_key = ['instance_a', 'instance_b']

# ── 3. Difficulty Classification (using base algorithms only) ─────────────────
base_data = all_data[all_data['algo'].isin(BASE_ALGOS)].copy()
base_data['solve_time'] = base_data.apply(
    lambda r: r['unified_time'] if r['unified_aborted'] == 0 else np.nan, axis=1
)

time_pivot = base_data.pivot_table(
    index=inst_key, columns='algo', values='solve_time', aggfunc='first'
)
abort_pivot = base_data.pivot_table(
    index=inst_key, columns='algo', values='unified_aborted', aggfunc='first'
)
abort_any = (abort_pivot.fillna(1).astype(int) == 1).astype(int)
all_instances = base_data[inst_key].drop_duplicates().set_index(inst_key).index

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
print('Difficulty distribution (base algorithms):')
print(difficulty.value_counts())


medium_instances = set(difficulty[difficulty == 'medium'].index)


print('Medium instances count:', len(medium_instances))
print('Example medium instance:', list(medium_instances)[:2])

lum = all_data[all_data['algo'] == 'll_lum']
print('ll_lum total rows:', len(lum))
print('ll_lum example instance:', list(zip(lum['instance_a'].head(2), lum['instance_b'].head(2))))

# Check how many ll_lum instances are in medium_instances
lum_idx = set(zip(lum['instance_a'], lum['instance_b']))
print('ll_lum instances in medium set:', len(lum_idx & medium_instances))


med_mask = all_data.set_index(inst_key).index.isin(medium_instances)
medium_data = all_data[med_mask].copy()

# ── 4. Solved count helper ────────────────────────────────────────────────────
def medium_summary(data, algos):
    rows = []
    total = len(medium_instances)
    for algo in algos:
        sub = data[data['algo'] == algo]
        if sub.empty:
            rows.append({'algo': algo, 'solved': np.nan, 'solved_pct': np.nan})
            continue
        solved = (sub['unified_aborted'] == 0).sum()
        rows.append({'algo': algo, 'solved': solved, 'solved_pct': 100 * solved / total})
    return pd.DataFrame(rows).set_index('algo')

# ── 5. Full ablation summary table ───────────────────────────────────────────
print('\n=== Full Ablation Summary (medium instances) ===\n')
full_summary = medium_summary(medium_data, ALL_ALGOS)
print(full_summary.round(2))

BASELINE_PCT = full_summary.loc['mcsplit', 'solved_pct']

# ── 6. Dot plot ───────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(9, 7))

y = 0
yticks, yticklabels = [], []
group_midpoints = {}

for group, algos in reversed(list(GROUPS.items())):
    group_start = y
    for algo in algos:
        pct = full_summary.loc[algo, 'solved_pct'] if algo in full_summary.index else np.nan
        if np.isnan(pct):
            yticks.append(y); yticklabels.append(LABELS[algo]); y += 1; continue
        color = GROUP_COLORS[group]
        ax.plot(pct, y, 'o', color=color, markersize=7, zorder=3)
        ax.plot([BASELINE_PCT, pct], [y, y], color=color, alpha=0.3, linewidth=1, zorder=2)
        yticks.append(y); yticklabels.append(LABELS[algo]); y += 1
    group_midpoints[group] = (group_start + y - 1) / 2
    y += 0.8

ax.axvline(BASELINE_PCT, color='#555555', linestyle='--', linewidth=1.2,
           label=f'McSplit baseline ({BASELINE_PCT:.1f}%)', zorder=1)
ax.set_yticks(yticks); ax.set_yticklabels(yticklabels, fontsize=9)
ax.set_xlabel('Solved (% of medium instances)', fontsize=10)
ax.set_title('Ablation Study - Medium Instances (MIVIA)', fontsize=11, fontweight='bold')
ax.legend(fontsize=8, loc='upper left')
ax.set_xlim(74, 101)
ax.grid(axis='x', linestyle=':', alpha=0.5)
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

ax2 = ax.twinx()
ax2.set_ylim(ax.get_ylim())
ax2.set_yticks(list(group_midpoints.values()))
ax2.set_yticklabels(list(group_midpoints.keys()), fontsize=8.5)
ax2.tick_params(length=0)
ax2.spines['top'].set_visible(False)
ax2.spines['left'].set_visible(False)
ax2.spines['bottom'].set_visible(False)
for label, group in zip(ax2.get_yticklabels(), group_midpoints.keys()):
    label.set_color(GROUP_COLORS[group])

plt.tight_layout()
plt.savefig('evaluation/results/ablation_dotplot.png', bbox_inches='tight')
plt.show()
print('Saved ablation_dotplot.png')

# ── 7. Symmetry-pruned branches on MIVIA ─────────────────────────────────────
print('\n=== Symmetry-Pruned Branches (medium instances, MIVIA) ===\n')
print('Shows how often symmetry mechanisms fire - expected to be low on synthetic graphs.\n')

sym_algos = ['rrsplit', 'rrsplit_noveq', 'symsplit', 'symsplit_varonly', 'symsplit_valonly']
sym_data = medium_data[medium_data['algo'].isin(sym_algos) & (medium_data['unified_aborted'] == 0)]

sym_summary = sym_data.groupby('algo')['unified_sym_pruned'].agg(['mean', 'median', 'sum'])
sym_summary.index = [LABELS.get(a, a) for a in sym_summary.index]
sym_summary.columns = ['Mean pruned', 'Median pruned', 'Total pruned']
print(sym_summary.round(1))

# ── 8. Bound-pruned branches across families ──────────────────────────────────
print('\n=== Bound-Pruned Branches by Algorithm Family (medium instances, MIVIA) ===\n')

bound_algos = ['mcsplit', 'rl', 'll', 'dal', 'dsb', 'rrsplit', 'symsplit']
bound_data = medium_data[
    medium_data['algo'].isin(bound_algos) & (medium_data['unified_aborted'] == 0)
]

bound_summary = bound_data.groupby('algo')['unified_bound_pruned'].agg(['mean', 'median'])
bound_summary.index = [LABELS.get(a, a) for a in bound_summary.index]
bound_summary.columns = ['Mean bound-pruned', 'Median bound-pruned']
print(bound_summary.round(0))

fig, ax = plt.subplots(figsize=(8, 4))
medians = bound_data.groupby('algo')['unified_bound_pruned'].median().reindex(bound_algos)
colors  = ['#555555', '#e07b39', '#e0a800', '#c94040', '#7a5c9e', '#2e86ab', '#3a9e5f']
bars = ax.bar([LABELS[a] for a in bound_algos], medians.values, color=colors)
ax.set_ylabel('Median bound-pruned branches', fontsize=10)
ax.set_title('Bound-Pruned Branches by Algorithm (medium instances, MIVIA)', fontsize=10)
ax.tick_params(axis='x', rotation=20)
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig('evaluation/results/ablation_bound_pruned.png', bbox_inches='tight')
plt.show()
print('Saved ablation_bound_pruned.png')