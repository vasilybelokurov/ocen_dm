#!/usr/bin/env python3
"""Plot completed rung-0/rung-1 comparisons from preserved summaries."""
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

root = Path(__file__).resolve().parents[1]
output = root/'plots'
output.mkdir(exist_ok=True)
data_output = root/'results/plot_data'
data_output.mkdir(parents=True, exist_ok=True)
labels = [['rung0_K1', 'rung0_K2_nfw', 'rung0_K2_cored'],
          ['rung1_K1_slice', 'rung1_K2_nfw', 'rung1_K2_cored_slice']]
records, hashes = [], {}
for group in labels:
    row = []
    for label in group:
        path = root/'results/fits'/label/'summary.json'
        row.append(json.loads(path.read_text()))
        hashes[str(path.relative_to(root))] = hashlib.sha256(path.read_bytes()).hexdigest()
    records.append(row)
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
names = ['HST radial', 'HST tangential', 'MUSE LOS', 'Gaia radial', 'Gaia tangential']
for rung, color in [(0, '#858b93'), (1, '#276cad')]:
    x = np.arange(3)+(rung-.5)*.34
    bars = axes[0].bar(x, [d['chi2_ml_total'] for d in records[rung]], .34, color=color, label=f'Rung {rung}')
    axes[0].bar_label(bars, fmt='%.1f', fontsize=9, padding=3)
    axes[1].bar(np.arange(5)+(rung-.5)*.34,
                [v['chi2'] for v in records[rung][-1]['chi2_ml'].values()], .34, color=color)
axes[0].set(xticks=np.arange(3), xticklabels=['K1', 'NFW', 'Core'], ylabel=r'Total $\chi^2$ (89 measurements)',
            title='All three families improve with constant anisotropy', ylim=(0, 850))
axes[1].set(xticks=np.arange(5), xticklabels=names, ylabel=r'$\chi^2$ contribution',
            title='Cored model: where the improvement occurs')
axes[1].tick_params(axis='x', labelrotation=30)
axes[0].legend()
for ax in axes:
    ax.grid(axis='y', alpha=.2); ax.set_axisbelow(True)
fig.savefig(output/'rung0_rung1_chi2.png', dpi=180)
(data_output/'rung0_rung1_chi2_provenance.json').write_text(json.dumps(hashes, indent=2)+'\n')
