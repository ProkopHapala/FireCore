#!/usr/bin/env python3
"""
Compare diamond phonon bands with and without angle terms.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys

if len(sys.argv) < 3:
    print("Usage: python3 plot_angle_comparison.py <full_ff.npz> <bonds_only.npz>")
    sys.exit(1)

full_file = sys.argv[1]
bond_file = sys.argv[2]

data_full = np.load(full_file)
data_bond = np.load(bond_file)

kdist_full = data_full['kdist']
freqs_full = data_full['freqs']
kdist_bond = data_bond['kdist']
freqs_bond = data_bond['freqs']

unit = 'THz' if 'THz' in full_file else 'cm-1'

fig, ax = plt.subplots(figsize=(10, 6))

# Plot full forcefield (dotted, thicker)
for imode in range(freqs_full.shape[1]):
    ax.plot(kdist_full, freqs_full[:, imode], 'b--', lw=1.5, alpha=0.7, label='Full FF' if imode == 0 else '')

# Plot bond-only (solid, thinner)
for imode in range(freqs_bond.shape[1]):
    ax.plot(kdist_bond, freqs_bond[:, imode], 'r-', lw=0.8, alpha=0.7, label='Bond-only' if imode == 0 else '')

ax.set_ylabel(f'Frequency ({unit})')
ax.set_xlabel('k-path')
ax.set_title(f'Diamond phonons: Full FF vs Bond-only ({unit})')
ax.set_xlim(kdist_full[0], kdist_full[-1])
ax.legend()
plt.tight_layout()

outpng = f'diamond_angle_comparison_{unit}.png'
plt.savefig(outpng, dpi=150)
print(f"Comparison plot saved to: {outpng}")

# Print Γ-point comparison
print(f"\nΓ-point comparison ({unit}):")
print(f"  Full FF:    {np.sort(freqs_full[0])}")
print(f"  Bond-only:  {np.sort(freqs_bond[0])}")
print(f"  Difference: {np.sort(freqs_full[0]) - np.sort(freqs_bond[0])}")
