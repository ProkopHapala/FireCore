"""
XPBD_2D_utils.py - Utility functions for 2D XPBD simulator

Contains:
- LiveViz2D: Real-time matplotlib visualization
- Molecule setup functions (chain, triangle, from_xyz)
- Momentum computation utilities
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

# Add repo root to import pyBall.AtomicSystem (AtomicSystem uses relative imports)
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

from XPBD_2D import (build_neighs_bk_from_bonds_2d, make_bk_slots_2d, make_stiffness_from_bonds_2d)
from pyBall.AtomicSystem import AtomicSystem


class LiveViz2D:
    """Lightweight live 2D updater that keeps view persistent (like LivePortViz but 2D)."""
    def __init__(self, elems=None, view_scale=None):
        self.plt = plt
        self.elems = elems
        self.view_scale = view_scale  # Fixed viewport scale if provided
        plt.ion()
        self.fig, self.ax = plt.subplots(figsize=(10, 8))
        self.ax.set_aspect('equal')
        self.ax.grid(True, alpha=0.3)
        self.scatter_atoms = self.ax.scatter([], [], s=500, facecolors=(1.0, 0.0, 0.0, 0.25), edgecolors=(1.0, 0.0, 0.0, 0.8), linewidths=1.0, zorder=5)
        self.scatter_ports = self.ax.scatter([], [], s=120, c='k', marker='+', linewidths=1.2, zorder=6)
        self.lc_bonds = LineCollection([], colors=[(0.0, 0.6, 0.0, 0.6)], linewidths=1.5, zorder=2)
        self.ax.add_collection(self.lc_bonds)
        self.lc_atom2port = LineCollection([], colors=[(0.1, 0.3, 1.0, 0.7)], linewidths=0.8, zorder=3)
        self.ax.add_collection(self.lc_atom2port)
        self.lc_port2neigh = LineCollection([], colors=[(1.0, 0.0, 1.0, 0.7)], linewidths=0.8, zorder=3)
        self.ax.add_collection(self.lc_port2neigh)
        self.info_text = self.ax.text(0.02, 0.98, '', transform=self.ax.transAxes,
                                      verticalalignment='top', fontsize=9,
                                      bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        # Note: blitting currently disabled in test script; keep artists non-animated so they draw.
        self.scatter_atoms.set_animated(False)
        self.scatter_ports.set_animated(False)
        self.lc_bonds.set_animated(False)
        self.lc_atom2port.set_animated(False)
        self.lc_port2neigh.set_animated(False)
        self.info_text.set_animated(False)
        self.ax.set_xlabel('X')
        self.ax.set_ylabel('Y')
        self.ax.set_title('XPBD_2D')
        self._last_pos = None
        self._inited_view = False
        self.fig.show()

    def artists(self):
        return [self.scatter_atoms, self.scatter_ports, self.lc_bonds, self.lc_atom2port, self.lc_port2neigh, self.info_text]

    def init_view(self, pos, *, margin=2.0):
        if self._inited_view:
            return
        pos = np.asarray(pos)
        if self.view_scale is not None:
            # Use fixed scale centered on the initial positions
            cx = (pos[:, 0].min() + pos[:, 0].max()) * 0.5
            cy = (pos[:, 1].min() + pos[:, 1].max()) * 0.5
            half_scale = self.view_scale * 0.5
            xmin, xmax = cx - half_scale, cx + half_scale
            ymin, ymax = cy - half_scale, cy + half_scale
        else:
            xmin, xmax = pos[:, 0].min() - margin, pos[:, 0].max() + margin
            ymin, ymax = pos[:, 1].min() - margin, pos[:, 1].max() + margin
        self.ax.set_xlim(xmin, xmax)
        self.ax.set_ylim(ymin, ymax)
        self._inited_view = True

    def update(self, pos, neighs=None, nnode=None, title="", info="", port_local=None, port_n=None, rot=None):
        self.ax.set_title(title)
        self._last_pos = np.asarray(pos)
        self.scatter_atoms.set_offsets(pos)
        self.init_view(pos)

        if neighs is not None and nnode is not None:
            segs = []
            for i in range(int(nnode)):
                for k in range(4):
                    j = int(neighs[i, k])
                    if j >= 0:
                        segs.append([pos[i], pos[j]])
            self.lc_bonds.set_segments(segs)
        else:
            self.lc_bonds.set_segments([])

        if (port_local is not None) and (port_n is not None) and (nnode is not None):
            pps = []
            seg_atom2port = []
            seg_port2neigh = []
            if rot is None:
                for i in range(int(nnode)):
                    for k in range(4):
                        j = int(neighs[i, k]) if neighs is not None else -1
                        if j < 0:
                            continue
                        tip = pos[i] + port_local[i, k]
                        pps.append(tip)
                        seg_atom2port.append([pos[i], tip])
                        seg_port2neigh.append([tip, pos[j]])
            else:
                for i in range(int(nnode)):
                    zi = rot[i]
                    for k in range(4):
                        j = int(neighs[i, k]) if neighs is not None else -1
                        if j < 0:
                            continue
                        p = port_local[i, k]
                        pr = np.array([zi[0] * p[0] - zi[1] * p[1], zi[1] * p[0] + zi[0] * p[1]], dtype=np.float32)
                        tip = pos[i] + pr
                        pps.append(tip)
                        seg_atom2port.append([pos[i], tip])
                        seg_port2neigh.append([tip, pos[j]])
            self.scatter_ports.set_offsets(np.array(pps, dtype=np.float32) if len(pps) else np.empty((0, 2), dtype=np.float32))
            self.lc_atom2port.set_segments(seg_atom2port)
            self.lc_port2neigh.set_segments(seg_port2neigh)
        else:
            self.scatter_ports.set_offsets(np.empty((0, 2), dtype=np.float32))
            self.lc_atom2port.set_segments([])
            self.lc_port2neigh.set_segments([])
        self.info_text.set_text(info)

        # NOTE: For performance (blitting), we do NOT call plt.pause() here.
        # The caller (FuncAnimation) drives the GUI loop.

        return self.artists()


def setup_from_mol(sim, mol, *, k_bond=200.0, perturbation=0.0, perturb_rot=0.0, bAllNodes=False, seed=0):
    """Setup simulator from an already loaded AtomicSystem."""
    n_atoms = len(mol.apos)
    
    if mol.bonds is None or len(mol.bonds) == 0:
        mol.findBonds()
    
    bonds = [(int(b[0]), int(b[1])) for b in mol.bonds]
    neighs, bks = build_neighs_bk_from_bonds_2d(n_atoms, bonds)
    
    # Determine nodes (atoms with >1 bond) vs capping atoms
    if bAllNodes:
        nnode = n_atoms
    else:
        nnode = sum(1 for i in range(n_atoms) if np.sum(neighs[i] >= 0) > 1)
    
    # Reorder: nodes first, then capping
    node_mask = np.zeros(n_atoms, dtype=bool)
    for i in range(n_atoms):
        if bAllNodes or np.sum(neighs[i] >= 0) > 1:
            node_mask[i] = True
    
    old_to_new = {}
    new_idx = 0
    for i in range(n_atoms):
        if node_mask[i]:
            old_to_new[i] = new_idx
            new_idx += 1
    for i in range(n_atoms):
        if not node_mask[i]:
            old_to_new[i] = new_idx
            new_idx += 1
    
    new_pos = np.zeros((n_atoms, 2), dtype=np.float32)
    for old_i, new_i in old_to_new.items():
        new_pos[new_i] = [mol.apos[old_i, 0], mol.apos[old_i, 1]]  # Take x,y only
    
    new_bonds = []
    for (i, j) in bonds:
        ni = int(old_to_new[int(i)])
        nj = int(old_to_new[int(j)])
        if ni == nj:
            continue
        if ni < nj:
            new_bonds.append((ni, nj))
        else:
            new_bonds.append((nj, ni))
    new_bonds = sorted(set(new_bonds))

    new_neighs, new_bks = build_neighs_bk_from_bonds_2d(n_atoms, new_bonds)
    
    bkSlots = make_bk_slots_2d(new_neighs, nnode=nnode, natoms=n_atoms)
    stiffness = make_stiffness_from_bonds_2d(n_atoms, new_neighs, k_bond=k_bond)
    
    sim.upload_topology(new_neighs, bkSlots, stiffness, nnode=nnode)
    
    # FIX: Only use first nnode entries for node arrays
    neighs_nodes = new_neighs[:nnode]  # Only nodes have ports
    stiffness_nodes = stiffness[:nnode]  # Only nodes have stiffness
    
    port_local = np.zeros((nnode, 4, 2), dtype=np.float32)
    port_n = np.zeros((nnode,), dtype=np.uint8)
    
    for i in range(nnode):
        n_ports = 0
        for k in range(4):
            j = new_neighs[i, k]
            if j < 0:
                continue
            dx = new_pos[j, 0] - new_pos[i, 0]
            dy = new_pos[j, 1] - new_pos[i, 1]
            port_local[i, k] = [dx, dy]
            n_ports += 1
        port_n[i] = n_ports
    
    sim.upload_ports(port_local, port_n, nnode=nnode)
    
    if perturbation > 0:
        np.random.seed(int(seed))
        new_pos[:nnode] += np.random.randn(nnode, 2) * perturbation
    rot = np.zeros((n_atoms, 2), dtype=np.float32)
    rot[:, 0] = 1.0
    if perturb_rot > 0:
        np.random.seed(int(seed) + 10007)
        ang = (np.random.randn(nnode).astype(np.float32)) * float(perturb_rot)
        rot[:nnode] = rot_from_angles(ang)
    sim.upload_state(new_pos, rot=rot)
    
    return {
        'neighs': neighs_nodes,  # Only nodes, not all atoms
        'bks': new_bks[:nnode],  # Only nodes
        'bkSlots': bkSlots,
        'stiffness': stiffness_nodes,  # Only nodes
        'port_local': port_local,
        'port_n': port_n,
        'bond_length': 1.0,
        'nnode': nnode,
        'elems': [mol.enames[old_to_new[i]] for i in range(n_atoms)] if mol.enames else None
    }


def setup_from_xyz(sim, xyz_path, k_bond=200.0, perturbation=0.0, perturb_rot=0.0, bAllNodes=False, seed=0):
    mol = AtomicSystem(fname=xyz_path)
    return setup_from_mol(sim, mol, k_bond=k_bond, perturbation=perturbation, perturb_rot=perturb_rot, bAllNodes=bAllNodes, seed=seed)


def rot_from_angles(ang):
    ang = np.asarray(ang, dtype=np.float32)
    return np.stack([np.cos(ang), np.sin(ang)], axis=-1).astype(np.float32)


def compute_port_error(pos, rot, neighs, bks, port_local, nnode, port_n=None):
    nnode = int(nnode)
    max_err = 0.0
    sum2 = 0.0
    nerr = 0
    for i in range(nnode):
        zi = rot[i]
        for k in range(4):
            j = int(neighs[i, k])
            if j < 0:
                continue
            pi = port_local[i, k]
            pri = np.array([zi[0] * pi[0] - zi[1] * pi[1], zi[1] * pi[0] + zi[0] * pi[1]], dtype=np.float32)
            di = (pos[i] + pri) - pos[j]
            err = float(np.sqrt(di[0] * di[0] + di[1] * di[1]))
            if err > max_err:
                max_err = err
            sum2 += err * err
            nerr += 1
    rms = float(np.sqrt(sum2 / max(1, nerr)))
    return max_err, rms


def attach_picker_2d(viz, sim, *, pick_radius=0.5, verbose=0):
    """Attach matplotlib callbacks to pick and drag atoms (2D).

    Dragging directly overwrites atom position on GPU; no silent error handling.
    """
    if viz is None:
        raise ValueError('attach_picker_2d: viz is None')
    if sim is None:
        raise ValueError('attach_picker_2d: sim is None')

    pick = {"idx": None, "active": False, "mouse": np.array([0.0, 0.0], dtype=np.float32)}

    def on_press(event):
        if event.inaxes != viz.ax or event.xdata is None or event.ydata is None:
            return
        if viz._last_pos is None:
            return
        mouse_xy = np.array([event.xdata, event.ydata], dtype=np.float32)
        d2 = np.sum((viz._last_pos[:, :2] - mouse_xy) ** 2, axis=1)
        i_min = int(np.argmin(d2))
        if d2[i_min] <= float(pick_radius) ** 2:
            pick["idx"] = i_min
            pick["active"] = True
            pick["mouse"] = mouse_xy
            if int(verbose) > 0:
                print(f"[DEBUG] pick press idx={i_min} d={np.sqrt(float(d2[i_min])):.4f} mouse={mouse_xy}")
        else:
            if int(verbose) > 0:
                print(f"[DEBUG] pick miss closest idx={i_min} d={np.sqrt(float(d2[i_min])):.4f} mouse={mouse_xy}")

    def on_release(event):
        if pick["active"] and int(verbose) > 0:
            print(f"[DEBUG] pick release idx={pick['idx']}")
        pick["active"] = False
        pick["idx"] = None

    def on_motion(event):
        if not pick["active"]:
            return
        if event.xdata is None or event.ydata is None:
            return
        pick["mouse"] = np.array([event.xdata, event.ydata], dtype=np.float32)
        ia = int(pick["idx"])
        sim.set_atom_pos(ia, pick["mouse"])
        if viz._last_pos is not None:
            viz._last_pos[ia, 0] = pick["mouse"][0]
            viz._last_pos[ia, 1] = pick["mouse"][1]

    viz.fig.canvas.mpl_connect('button_press_event', on_press)
    viz.fig.canvas.mpl_connect('button_release_event', on_release)
    viz.fig.canvas.mpl_connect('motion_notify_event', on_motion)
    return pick


def create_chain_topology(n_atoms, bond_length=1.0):
    """Create a linear chain of atoms with bonds between neighbors."""
    bonds = [(i, i+1) for i in range(n_atoms - 1)]
    return bonds


def setup_chain_molecule(sim, n_atoms, nnode, bond_length=1.0, perturbation=0.0):
    """Setup initial state for a chain molecule."""
    bonds = create_chain_topology(nnode, bond_length)
    neighs, bks = build_neighs_bk_from_bonds_2d(n_atoms, bonds)
    bkSlots = make_bk_slots_2d(neighs, nnode=nnode, natoms=n_atoms)
    stiffness = make_stiffness_from_bonds_2d(n_atoms, neighs, k_bond=200.0)
    
    sim.upload_topology(neighs, bkSlots, stiffness)
    
    port_local = np.zeros((nnode, 4, 2), dtype=np.float32)
    port_n = np.zeros((nnode,), dtype=np.uint8)
    
    for i in range(nnode):
        n_ports = 0
        for k in range(4):
            j = neighs[i, k]
            if j < 0:
                continue
            if j > i:
                port_local[i, n_ports] = [bond_length * 0.5, 0.0]
            else:
                port_local[i, n_ports] = [-bond_length * 0.5, 0.0]
            n_ports += 1
        port_n[i] = n_ports
    
    sim.upload_ports(port_local, port_n, nnode=nnode)
    
    pos = np.zeros((n_atoms, 2), dtype=np.float32)
    for i in range(n_atoms):
        pos[i] = [i * bond_length, 0.0]
    
    if perturbation > 0:
        pos[:nnode] += np.random.randn(nnode, 2) * perturbation
    
    sim.upload_state(pos)
    
    return {
        'neighs': neighs,
        'bkSlots': bkSlots,
        'stiffness': stiffness,
        'port_local': port_local,
        'port_n': port_n,
        'bond_length': bond_length
    }


def setup_triangle_molecule(sim, bond_length=1.0, perturbation=0.0):
    """Setup a triangle of 3 rigid nodes."""
    n_atoms = 6
    nnode = 3
    bonds = [(0, 1), (1, 2), (2, 0)]
    
    neighs, bks = build_neighs_bk_from_bonds_2d(n_atoms, bonds)
    bkSlots = make_bk_slots_2d(neighs, nnode=nnode, natoms=n_atoms)
    stiffness = make_stiffness_from_bonds_2d(n_atoms, neighs, k_bond=300.0)
    
    sim.upload_topology(neighs, bkSlots, stiffness)
    
    port_local = np.zeros((nnode, 4, 2), dtype=np.float32)
    port_n = np.full((nnode,), 2, dtype=np.uint8)
    
    angles = [np.pi/2, -np.pi/6, -5*np.pi/6]
    for i in range(nnode):
        a1 = angles[i]
        a2 = angles[(i+1) % 3]
        port_local[i, 0] = [bond_length * 0.5 * np.cos(a1), bond_length * 0.5 * np.sin(a1)]
        port_local[i, 1] = [bond_length * 0.5 * np.cos(a2), bond_length * 0.5 * np.sin(a2)]
    
    sim.upload_ports(port_local, port_n, nnode=nnode)
    
    pos = np.zeros((n_atoms, 2), dtype=np.float32)
    for i in range(nnode):
        pos[i] = [bond_length * np.cos(angles[i]), bond_length * np.sin(angles[i])]
    
    if perturbation > 0:
        pos[:nnode] += np.random.randn(nnode, 2) * perturbation
    
    sim.upload_state(pos)
    
    return {
        'neighs': neighs,
        'bkSlots': bkSlots,
        'stiffness': stiffness,
        'port_local': port_local,
        'port_n': port_n,
        'bond_length': bond_length,
        'nnode': nnode
    }


def compute_momentum_2d(pos, vel, omega=None):
    """Compute total linear and angular momentum."""
    n = pos.shape[0]
    mass = 1.0
    
    P = np.sum(vel, axis=0) * mass
    
    L = 0.0
    for i in range(n):
        r = pos[i]
        p = vel[i] * mass
        L += r[0] * p[1] - r[1] * p[0]
    
    if omega is not None:
        I = 0.4 * mass
        L += np.sum(omega) * I
    
    return P, L
