#!/usr/bin/env python3
"""
Compare MMFF phonon bands with different forcefield components.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys

if len(sys.argv) < 4:
    print("Usage: python3 plot_mmff_comparison.py <bondsonly.npz> <bondangle.npz> <full.npz>")
    sys.exit(1)

bond_file = sys.argv[1]
bondangle_file = sys.argv[2]
full_file = sys.argv[3]

data_bond = np.load(bond_file)
data_bondangle = np.load(bondangle_file)
data_full = np.load(full_file)

kdist = data_bond['kdist']
freqs_bond = data_bond['freqs']
freqs_bondangle = data_bondangle['freqs']
freqs_full = data_full['freqs']

unit = 'THz' if 'THz' in bond_file else 'cm-1'

fig, ax = plt.subplots(figsize=(12, 6))

# Plot bond-only (red, solid, thin)
for imode in range(freqs_bond.shape[1]):
    ax.plot(kdist, freqs_bond[:, imode], 'r-', lw=0.8, alpha=0.7, label='Bond-only' if imode == 0 else '')

# Plot bond+angle (green, dashed, medium)
for imode in range(freqs_bondangle.shape[1]):
    ax.plot(kdist, freqs_bondangle[:, imode], 'g--', lw=1.2, alpha=0.7, label='Bond+Angle' if imode == 0 else '')

# Plot full (blue, dotted, thick)
for imode in range(freqs_full.shape[1]):
    ax.plot(kdist, freqs_full[:, imode], 'b:', lw=1.8, alpha=0.7, label='Full (Bond+Angle+PiSigma+PiPiI)' if imode == 0 else '')

ax.set_ylabel(f'Frequency ({unit})')
ax.set_xlabel('k-path')
ax.set_title(f'Diamond phonons: MMFF forcefield components ({unit})')
ax.set_xlim(kdist[0], kdist[-1])
ax.legend(loc='upper right')
plt.tight_layout()

outpng = f'diamond_mmff_comparison_{unit}.png'
plt.savefig(outpng, dpi=150)
print(f"MMFF comparison plot saved to: {outpng}")

# Print Γ-point comparison
print(f"\nΓ-point comparison ({unit}):")
print(f"  Bond-only:         {np.sort(freqs_bond[0])}")
print(f"  Bond+Angle:         {np.sort(freqs_bondangle[0])}")
print(f"  Full:               {np.sort(freqs_full[0])}")
print(f"\n  Bond+Angle - Bond-only:   {np.sort(freqs_bondangle[0]) - np.sort(freqs_bond[0])}")
print(f"  Full - Bond+Angle:       {np.sort(freqs_full[0]) - np.sort(freqs_bondangle[0])}")
