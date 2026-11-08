#!/usr/bin/env python3
"""Compare UFF phonon bands with gradually added forcefield components."""
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys, os

files = {
    'Bonds only':              'diamond_phonon_bands_THz_PBC_uff_noangles_nodihedrals_noinversions.npz',
    'Bonds+Angles':            'diamond_phonon_bands_THz_PBC_uff_nodihedrals_noinversions.npz',
    'Bonds+Angles+Dihedrals':  'diamond_phonon_bands_THz_PBC_uff.npz',
}
colors  = ['r',  'g',  'b' ]
styles  = ['-',  '--', ':' ]
widths  = [1.0,  1.4,  1.8 ]

fig, (ax, ax2) = plt.subplots(1, 2, figsize=(16, 6))

kdist = None
gamma_freqs = {}
for (label, fname), col, ls, lw in zip(files.items(), colors, styles, widths):
    if not os.path.exists(fname):
        print(f"WARNING: {fname} not found, skipping")
        continue
    d = np.load(fname)
    if kdist is None: kdist = d['kdist']
    freqs = d['freqs']
    gamma_freqs[label] = np.sort(freqs[0])
    for i in range(freqs.shape[1]):
        ax.plot(kdist, freqs[:, i], color=col, ls=ls, lw=lw, alpha=0.7,
                label=label if i == 0 else '')
    # also plot on ax2 with y-range clipped to [-5, 80]
    for i in range(freqs.shape[1]):
        ax2.plot(kdist, freqs[:, i], color=col, ls=ls, lw=lw, alpha=0.7,
                 label=label if i == 0 else '')

ax.axhline(0, color='k', lw=0.5, ls='--')
ax.set_ylabel('Frequency (THz)')
ax.set_xlabel('k-path')
ax.set_title('UFF diamond phonons: all components')
ax.set_xlim(kdist[0], kdist[-1])
ax.legend(loc='upper right', fontsize=8)

ax2.axhline(0, color='k', lw=0.5, ls='--')
ax2.set_ylim(-5, 75)
ax2.set_ylabel('Frequency (THz)')
ax2.set_xlabel('k-path')
ax2.set_title('UFF diamond phonons: clipped [-5,75 THz]')
ax2.set_xlim(kdist[0], kdist[-1])
ax2.legend(loc='upper right', fontsize=8)

plt.tight_layout()
outpng = 'diamond_uff_comparison_THz.png'
plt.savefig(outpng, dpi=150)
print(f"UFF comparison plot saved to: {outpng}")

print("\nΓ-point optical frequencies (THz):")
for label, gf in gamma_freqs.items():
    print(f"  {label:35s}: {gf}")
