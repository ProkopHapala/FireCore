import numpy as np
import matplotlib.pyplot as plt

filename = 'MYFILENAME'  # file to plot
MIN, MAX = MYMIN, MYMAX  # color limits
title    = 'MYTITLE'
savefile = 'MYPNG'

def read_block_data(fname):
    """
    Reads data with blank-line-separated blocks and returns
    Theta, R, Z arrays, along with inferred shape ns.
    """
    with open(fname) as f:
        lines = f.readlines()

    # Split data into blocks by blank lines
    blocks = []
    current_block = []
    for line in lines:
        if line.strip() == "":
            if current_block:
                blocks.append(np.array(current_block, dtype=float))
                current_block = []
        else:
            current_block.append(line.split())
    if current_block:
        blocks.append(np.array(current_block, dtype=float))

    # Infer ns
    n_r = len(blocks)       # number of radial points
    n_theta = len(blocks[0])  # number of angular points
    ns = (n_theta, n_r)

    # Build Theta, R, Z arrays
    Z = np.zeros(ns)
    Theta = np.zeros(ns)
    R = np.zeros(ns)

    for j, block in enumerate(blocks):
        Theta[:, j] = block[:, 0]
        R[:, j] = block[:, 1]
        Z[:, j] = block[:, 2]

    return Theta, R, Z, ns

def plot_polar(fname, title=None, clim=None, cmap='bwr', savefile=None):
    Theta, R, Z, ns = read_block_data(fname)

    # Convert theta to radians and shift
    Theta_rad = np.radians(Theta)

    # Clamp values
    if clim is not None:
        Z = np.clip(Z, clim[0], clim[1])

    # Create figure and axis
    fig, ax = plt.subplots(subplot_kw=dict(projection='polar'), figsize=(6, 6))

    # Contour plot
    cax = ax.contourf(Theta_rad, R, Z, 100, cmap=cmap, vmin=clim[0], vmax=clim[1])

    # Polar limits
    ax.set_theta_zero_location("W")   # zero points to the left
    ax.set_theta_direction(1)         # angles increase counter-clockwise
    ax.set_thetamin(-90)
    ax.set_thetamax(90)
    ax.set_rmin(0.0)
    ax.set_rmax(np.max(R))
    ax.set_rlabel_position(180)   # place radial labels on the left axis
    ax.text(np.pi, ax.get_rmax()*0.2, "distance (Ang)", rotation=90, va="center", ha="center")
    
    # Colorbar
    cb = fig.colorbar(cax, ax=ax, shrink=0.7, aspect=20, pad=-0.07)
    cb.set_label("Energy (kcal/mol)")

    # Title
    if title is not None:
        ax.set_title(title)

    # Save figure if requested
    if savefile is not None:
        fig.savefig(savefile, bbox_inches='tight', dpi=300)

plot_polar(filename, title=title, clim=(MIN, MAX), savefile=savefile)
