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

ALGOS = ['mcsplit', 'rl', 'll', 'dal', 'dsb', 'rrsplit', 'symsplit']
DATA_DIR = Path('hpc/cross-dataset-evaluation/results')

COLS = [
    'instance_a', 'instance_b', 'algo',
    'unified_size', 'unified_edges', 'unified_nodes', 'unified_time',
    'unified_aborted', 'unified_nodes_to_best',
    'unified_time_to_best', 'unified_cut_branches', 'unified_bound_pruned',
    'unified_sym_pruned',
    'ref_size', 'ref_nodes', 'ref_time', 'match'
]

# ── Load data for all algorithms ────────────────────────────────────────
dfs = []
for algo in ALGOS:
    path = DATA_DIR / f'{algo}_lv_all.csv'
    if not path.exists():
        print(f'WARNING: {path} not found, skipping {algo}')
        continue
    df = pd.read_csv(path, skipinitialspace=True)
    df = df.reindex(columns=COLS)
    dfs.append(df)

all_data = pd.concat(dfs, ignore_index=True)
for col in ['unified_size', 'unified_aborted', 'unified_time', 'unified_nodes',
            'ref_size', 'ref_time', 'ref_nodes']:
    all_data[col] = pd.to_numeric(all_data[col], errors='coerce')

inst_key = ['instance_a', 'instance_b']

# ── Build per-instance consensus from every completed size (unified + reference,
# across all seven algorithms), then flag instances where DAL's reference is
# exactly one larger than that consensus ─────────────────────────────────────
completed = all_data[
    (all_data['unified_aborted'] == 0) & all_data['unified_size'].notna()
][inst_key + ['algo', 'unified_size']].rename(columns={'unified_size': 'size'})
completed['source'] = completed['algo'] + '_unified'

ref_completed = all_data[
    all_data['ref_size'].notna()
][inst_key + ['algo', 'ref_size']].rename(columns={'ref_size': 'size'})
ref_completed['source'] = ref_completed['algo'] + '_ref'

everything = pd.concat([completed, ref_completed], ignore_index=True)

def analyse_instance(group):
    dal_ref_row = group[group['source'] == 'dal_ref']
    if dal_ref_row.empty:
        return None
    dal_ref_size = dal_ref_row['size'].iloc[0]

    others = group[group['source'] != 'dal_ref']
    if len(others) < 2:  # need enough corroborating sources to trust a consensus
        return None

    consensus_sizes = others['size'].mode()
    if len(consensus_sizes) != 1:
        return None  # no clean majority agreement, skip
    consensus = consensus_sizes.iloc[0]

    agreement = (others['size'] == consensus).mean()
    if agreement < 0.75:  # require most other sources to agree before trusting consensus
        return None

    return pd.Series({
        'consensus_size': consensus,
        'dal_ref_size': dal_ref_size,
        'diff': dal_ref_size - consensus,
        'n_corroborating': len(others),
    })

results = everything.groupby(inst_key).apply(analyse_instance).dropna(how='all').reset_index()

# Instances where DAL's reference is exactly +1 over consensus
plus_one = results[results['diff'] == 1].copy()
print(f'Instances where DAL reference = consensus + 1: {len(plus_one)}')
print(f'Instances where DAL reference matches consensus: {(results["diff"] == 0).sum()}')
print(f'Other discrepancies (diff != 0, 1): {((results["diff"] != 0) & (results["diff"] != 1)).sum()}')

# ── Pull node/time data for the flagged instances, DAL only, for the scatter ──
dal_data = all_data[all_data['algo'] == 'dal'].copy()
flagged = dal_data.merge(plus_one[inst_key], on=inst_key, how='inner')

plot_data = flagged[
    flagged['unified_nodes'].notna() & flagged['ref_nodes'].notna() &
    (flagged['unified_nodes'] > 0) & (flagged['ref_nodes'] > 0) &
    flagged['unified_time'].notna() & flagged['ref_time'].notna() &
    (flagged['unified_time'] > 0) & (flagged['ref_time'] > 0) &
    (flagged['unified_aborted'] == 0)
].copy()

print(f'Plottable +1 instances: {len(plot_data)}')

# ── Plot ───────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))

ax = axes[0]
ax.scatter(plot_data['ref_nodes'], plot_data['unified_nodes'],
           alpha=0.4, s=60, color='#c94040', edgecolors='none')
lim = max(plot_data['ref_nodes'].max(), plot_data['unified_nodes'].max()) * 1.1
lo = min(plot_data['ref_nodes'].min(), plot_data['unified_nodes'].min()) * 0.9
ax.plot([lo, lim], [lo, lim], 'k--', linewidth=0.8, alpha=0.5)
# ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlabel('Reference nodes')
ax.set_ylabel('Unified nodes')
ax.set_title('Nodes Explored')

ax = axes[1]
ax.scatter(plot_data['ref_time'], plot_data['unified_time'],
           alpha=0.4, s=60, color='#c94040', edgecolors='none')
lim = max(plot_data['ref_time'].max(), plot_data['unified_time'].max()) * 1.1
lo = min(plot_data['ref_time'].min(), plot_data['unified_time'].min()) * 0.9
ax.plot([lo, lim], [lo, lim], 'k--', linewidth=0.8, alpha=0.5)
# ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlabel('Reference time (s)')
ax.set_ylabel('Unified time (s)')
ax.set_title('Wall-Clock Time')

# fig.suptitle(f'Instances where DAL Reference Size = Consensus + 1 (n={len(plot_data)})\n'
#              'Dashed line = perfect agreement', fontsize=11)
plt.tight_layout()
# plt.savefig('evaluation/results/bi_dal_bug_analysis.png', bbox_inches='tight', dpi=500)
plt.show()