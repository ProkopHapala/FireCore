#!/usr/bin/env python3
"""Standalone phonon band plotter (no ASan issues)."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys

if len(sys.argv) < 2:
    print("Usage: python3 plot_phonon_bands.py <npz_file> [--unit THz|cm-1]")
    sys.exit(1)

npz_path = sys.argv[1]
unit_label = 'THz'
if '--unit' in sys.argv:
    idx = sys.argv.index('--unit')
    if idx + 1 < len(sys.argv):
        unit_label = sys.argv[idx + 1]

data = np.load(npz_path)
kdist = data['kdist']
freqs = data['freqs']
super_n = int(data['super_n']) if 'super_n' in data.files else 3
pbc = bool(data['pbc']) if 'pbc' in data.files else False
asr = bool(data['asr']) if 'asr' in data.files else False

mode_tag = 'PBC' if pbc else ('asr' if asr else 'cluster')

print(f"Loaded: {npz_path}")
print(f"super_n={super_n}, pbc={pbc}, asr={asr}")
print(f"freqs shape: {freqs.shape}")
print(f"freqs max: {np.max(freqs):.2f} {unit_label}")
print(f"freqs min: {np.min(freqs):.2f} {unit_label}")
print(f"Gamma point (sorted): {np.sort(freqs[0])}")

fig, ax = plt.subplots(figsize=(8, 6))
for imode in range(freqs.shape[1]):
    ax.plot(kdist, freqs[:, imode], 'b-', lw=0.8)

ax.set_ylabel(f'Frequency ({unit_label})')
ax.set_xlabel('k-path')
ax.set_title(f'Diamond phonons ({unit_label}, {super_n}x{super_n}x{super_n} {mode_tag})')
ax.set_xlim(kdist[0], kdist[-1])
ax.set_ylim(0, np.max(freqs) * 1.05)

plt.tight_layout()
outpng = npz_path.replace('.npz', '.png')
plt.savefig(outpng, dpi=150)
print(f"Saved: {outpng}")
