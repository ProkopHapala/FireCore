#!/usr/bin/env python3
"""
Compare MMFF phonon bands with parameter scaling (default, bond scaled, angle scaled).
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys

if len(sys.argv) < 4:
    print("Usage: python3 plot_mmff_scaling.py <default.npz> <bond_scaled.npz> <angle_scaled.npz>")
    sys.exit(1)

default_file = sys.argv[1]
bond_file = sys.argv[2]
angle_file = sys.argv[3]

data_default = np.load(default_file)
data_bond = np.load(bond_file)
data_angle = np.load(angle_file)

kdist = data_default['kdist']
freqs_default = data_default['freqs']
freqs_bond = data_bond['freqs']
freqs_angle = data_angle['freqs']

unit = 'THz' if 'THz' in default_file else 'cm-1'

fig, ax = plt.subplots(figsize=(12, 6))

# Plot default (black, solid, medium)
for imode in range(freqs_default.shape[1]):
    ax.plot(kdist, freqs_default[:, imode], 'k-', lw=1.0, alpha=0.6, label='Default' if imode == 0 else '')

# Plot bond scaled (red, dashed, medium)
for imode in range(freqs_bond.shape[1]):
    ax.plot(kdist, freqs_bond[:, imode], 'r--', lw=1.0, alpha=0.6, label='Bond scaled (1.25x)' if imode == 0 else '')

# Plot angle scaled (blue, dotted, medium)
for imode in range(freqs_angle.shape[1]):
    ax.plot(kdist, freqs_angle[:, imode], 'b:', lw=1.0, alpha=0.6, label='Angle scaled (1.25x)' if imode == 0 else '')

ax.set_ylabel(f'Frequency ({unit})')
ax.set_xlabel('k-path')
ax.set_title(f'Diamond phonons: MMFF parameter scaling ({unit})')
ax.set_xlim(kdist[0], kdist[-1])
ax.legend(loc='upper right')
plt.tight_layout()

outpng = f'diamond_mmff_scaling_{unit}.png'
plt.savefig(outpng, dpi=150)
print(f"MMFF scaling plot saved to: {outpng}")

# Print Γ-point comparison
print(f"\nΓ-point comparison ({unit}):")
print(f"  Default:           {np.sort(freqs_default[0])}")
print(f"  Bond scaled (1.25x): {np.sort(freqs_bond[0])}")
print(f"  Angle scaled (1.25x): {np.sort(freqs_angle[0])}")
print(f"\n  Bond - Default:     {np.sort(freqs_bond[0]) - np.sort(freqs_default[0])}")
print(f"  Angle - Default:    {np.sort(freqs_angle[0]) - np.sort(freqs_default[0])}")
