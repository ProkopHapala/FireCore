#!/usr/bin/python3 -u

import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description='Parse arguments')
parser.add_argument('--filename', type=str,   required=True, help="data file to plot")
parser.add_argument('--Emin',     type=float, required=True, help="min color limit")
parser.add_argument('--Emax',     type=float, required=True, help="max color limit")
parser.add_argument('--title',    type=str,   required=True, help="title for plot")
parser.add_argument('--savefile', type=str,   required=True, help="file where to save the plot")
parser.add_argument('--dmax',     type=float, required=True, help="max distance")
args = parser.parse_args()

with open(args.filename) as f:
    lines = f.readlines()

# Split data into blocks by blank lines
blocks = []
current_block = []
for line in lines:
    line = line.strip()
    if not line:
        if current_block:
            blocks.append(np.array(current_block, dtype=float))
            current_block = []
        continue
    if line.startswith("#"):
        continue
    current_block.append(line.split())
if current_block:
    blocks.append(np.array(current_block, dtype=float))

# Infer dimensions
n_theta = len(blocks)          # number of angle values
n_r = len(blocks[0])           # number of distance values
ns = (n_theta, n_r)

# Build Theta, R, Z arrays
Theta = np.zeros(ns)
R = np.zeros(ns)
Z = np.zeros(ns)
for i, block in enumerate(blocks):
    Theta[i, :] = block[:, 0]
    R[i, :]     = block[:, 1]
    Z[i, :]     = block[:, 2]

# Convert theta to radians and shift
Theta_rad = np.radians(Theta)

# Clamp values
Z = np.clip( Z, a_min=args.Emin, a_max=args.Emax )

# Create figure and axis
fig, ax = plt.subplots(subplot_kw=dict(projection='polar'), figsize=(6, 6))

# Contour plot
cax = ax.contourf(Theta_rad, R, Z, 100, cmap='bwr', vmin=args.Emin, vmax=args.Emax)

# Polar limits
ax.set_theta_zero_location("W")   # zero points to the left
ax.set_theta_direction(1)         # angles increase counter-clockwise
ax.set_thetamin(-90)
ax.set_thetamax(90)
ax.set_rmin(0.0)
ax.set_rmax(args.dmax)
ax.set_rlabel_position(180)   # place radial labels on the left axis
ax.text(np.pi, ax.get_rmax()*0.2, "distance (Ang)", rotation=90, va="center", ha="center")

# Colorbar
cb = fig.colorbar(cax, ax=ax, shrink=0.7, aspect=20, pad=-0.07)
cb.set_label("Energy (kcal/mol)")

# Title
ax.set_title(args.title)

# Save figure
fig.savefig(args.savefile, bbox_inches='tight', dpi=300)
