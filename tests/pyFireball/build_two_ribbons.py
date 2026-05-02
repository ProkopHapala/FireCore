#!/usr/bin/env python3
"""Build two-ribbon unit cell for H-transfer NEB.

Creates a unit cell with two ribbons (N-passivated and NH-passivated)
stacked along y-direction with H-bond distance ~3.5 Å.
"""
import sys
import os
import numpy as np

sys.path.append("../../")
from doc.Topics.Kekule_Topology.GrapheneRibbonBuilder import GrapheneRibbonBuilder
from pyBall import plotUtils

ELEM_MAP = {'H': 1, 'C': 6, 'N': 7, 'O': 8}
ELEM_MAP_INV = {v: k for k, v in ELEM_MAP.items()}

def build_ribbon(passivation, width_chains, length_cells, Lx, a_CC=1.42):
    """Build a single ribbon geometry."""
    xa_nom = a_CC * np.cos(np.pi / 6)
    b = GrapheneRibbonBuilder(a_CC=a_CC)
    scale_x = Lx / (2.0 * xa_nom)
    pos2d, elems, bonds = b.build_zigzag_ribbon(
        width_chains=width_chains, length_cells=length_cells,
        passivation=passivation, scale_x=scale_x)
    atypes = np.array([ELEM_MAP[e] for e in elems], dtype=np.int32)
    return np.array(pos2d), atypes, elems

def build_two_ribbon_cell(width_chains=4, length_cells=1, Lx=2.4, a_CC=1.42, L_Hb=2.0):
    """Build two-ribbon unit cell with N and NH passivation for hydrogen bonding.
    
    Args:
        width_chains: ribbon width in chains
        length_cells: length in unit cells
        Lx: lattice constant along x
        a_CC: C-C bond length
        L_Hb: hydrogen bond length between ribbons (Å)
    
    Returns:
        apos: (natoms, 3) atomic positions
        atypes: (natoms,) atomic types
        elems: (natoms,) element names
        lvs: (3,3) lattice vectors
    """
    # Build N-passivated ribbon (bottom)
    pos2d_N, atypes_N, elems_N = build_ribbon('N', width_chains, length_cells, Lx, a_CC)
    
    # Build NH-passivated ribbon (top)
    pos2d_NH, atypes_NH, elems_NH = build_ribbon('NH', width_chains, length_cells, Lx, a_CC)
    
    # Convert to 3D
    apos_N = np.zeros((len(atypes_N), 3))
    apos_N[:, 0:2] = pos2d_N
    
    apos_NH = np.zeros((len(atypes_NH), 3))
    apos_NH[:, 0:2] = pos2d_NH
    
    # Center each ribbon in y
    apos_N[:, 1] -= np.mean(apos_N[:, 1])
    apos_NH[:, 1] -= np.mean(apos_NH[:, 1])
    
    # Calculate ribbon spans in y
    y_min_N = np.min(apos_N[:, 1])
    y_max_N = np.max(apos_N[:, 1])
    y_span_N = y_max_N - y_min_N
    
    y_min_NH = np.min(apos_NH[:, 1])
    y_max_NH = np.max(apos_NH[:, 1])
    y_span_NH = y_max_NH - y_min_NH
    
    # Shift NH ribbon so: y_min_NH = y_max_N + L_Hb
    shift_y = (y_max_N + L_Hb) - y_min_NH
    apos_NH[:, 1] += shift_y
    
    # Recalculate extents after shift
    y_min_NH = np.min(apos_NH[:, 1])
    y_max_NH = np.max(apos_NH[:, 1])
    
    print(f"Ribbon N: y_span={y_span_N:.2f} Å, extent [{y_min_N:.2f}, {y_max_N:.2f}]")
    print(f"Ribbon NH: y_span={y_span_NH:.2f} Å, extent [{y_min_NH:.2f}, {y_max_NH:.2f}]")
    print(f"Shifted NH ribbon by {shift_y:.3f} Å")
    print(f"H-bond gap: {y_min_NH - y_max_N:.2f} Å (target L_Hb={L_Hb} Å)")
    
    # Calculate cell size: Ly = y_span_N + y_span_NH + 2 * L_Hb
    Ly = y_span_N + y_span_NH + 2 * L_Hb
    
    # Combine ribbons
    apos = np.vstack([apos_N, apos_NH])
    atypes = np.concatenate([atypes_N, atypes_NH])
    elems = list(elems_N) + list(elems_NH)
    
    # Add z-separation (vacuum)
    z_center = 10.0
    apos[:, 2] = z_center
    
    # Lz is vacuum in z
    Lz = 20.0
    
    # Center in cell (y-direction centered on the combined structure)
    apos[:, 1] -= np.mean(apos[:, 1])
    apos[:, 2] -= np.mean(apos[:, 2])
    apos[:, 2] += 0.5 * Lz
    
    # Lattice vectors: periodic along x, vacuum along y, z
    lvs = np.array([[Lx, 0.0, 0.0], [0.0, Ly, 0.0], [0.0, 0.0, Lz]])
    
    return apos, atypes, elems, lvs

def plot_geometry(apos, atypes, elems, lvs, title="Two-ribbon unit cell"):
    """Plot the geometry using plotGeometry from plotUtils."""
    plotUtils.plotGeometry(
        apos=apos, 
        atypes=atypes, 
        lvs=lvs, 
        bond_dist=1.8, 
        bBondLabels=True, 
        axes=(0, 1), 
        title=title, 
        fname='two_ribbons_geometry.png',
        figsize=(10, 6)
    )
    print(f"Saved plot: two_ribbons_geometry.png")

def save_xyz(apos, elems, fname):
    """Save geometry to XYZ file."""
    with open(fname, 'w') as f:
        f.write(f"{len(elems)}\n")
        f.write("Two-ribbon unit cell\n")
        for elem, pos in zip(elems, apos):
            f.write(f"{elem} {pos[0]:.6f} {pos[1]:.6f} {pos[2]:.6f}\n")
    print(f"Saved XYZ: {fname}")

def main():
    # Parameters
    width_chains = 4
    length_cells = 1
    Lx = 2.4
    a_CC = 1.42
    L_Hb = 2.0  # Hydrogen bond length between ribbons
    
    # Build geometry
    apos, atypes, elems, lvs = build_two_ribbon_cell(
        width_chains=width_chains, length_cells=length_cells, 
        Lx=Lx, a_CC=a_CC, L_Hb=L_Hb)
    
    print(f"Built two-ribbon cell:")
    print(f"  Natoms: {len(elems)}")
    print(f"  Elements: {dict(zip(*np.unique(elems, return_counts=True)))}")
    print(f"  Lattice vectors: {lvs[0,0]:.2f} x {lvs[1,1]:.2f} x {lvs[2,2]:.2f} Å")
    
    # Plot
    plot_geometry(apos, atypes, elems, lvs, title=f"Two-ribbon cell (Lx={Lx} Å)")
    
    # Save XYZ
    save_xyz(apos, elems, 'two_ribbons.xyz')
    
    # Save GenFormat for DFTB+
    enameset = sorted(set(elems))
    ename_to_idx = {e: i+1 for i, e in enumerate(enameset)}
    with open('two_ribbons.gen', 'w') as f:
        f.write(f"  {len(elems)}  S\n")
        f.write('  ' + ' '.join(enameset) + '\n')
        for i, (elem, pos) in enumerate(zip(elems, apos)):
            idx = ename_to_idx[elem]
            f.write(f"  {i+1} {idx}   {pos[0]:.10f}   {pos[1]:.10f}   {pos[2]:.10f}\n")
        f.write('  0.000000000  0.000000000  0.000000000\n')
        for row in lvs:
            f.write(f'  {row[0]:.10f}  {row[1]:.10f}  {row[2]:.10f}\n')
    print("Saved GenFormat: two_ribbons.gen")

if __name__ == '__main__':
    main()
