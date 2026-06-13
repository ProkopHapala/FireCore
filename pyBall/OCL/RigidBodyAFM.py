import os
import numpy as np
from .RigidBodyDynamics import RigidBodyDynamics, _reqs_to_plq, _load_npy_legacy

def sample_gridff_single_atom(scan_positions, gridff_path, sub_xyz, atom_req=(1.487, 0.0006808, 0.0, 0.0), atom_mass=1.008, alpha_morse=1.5, debug=False, grid_p0=None, grid_step=None):
    """
    Samples the GridFF at the given scan_positions using a single test atom.
    scan_positions: (N, 3) array of coordinates.
    atom_req: tuple of (R, E, Q, H) for the test atom. Default is H atom.
    grid_p0: Grid origin (x0, y0, z0). If None, defaults to (0,0,0).
    grid_step: Grid spacing (dx, dy, dz). If None, calculated from lattice vectors.
    """
    scan_positions = np.asarray(scan_positions, dtype=np.float32)
    n_bodies = len(scan_positions)
    
    try:
        grid = np.load(gridff_path)
    except Exception:
        grid = _load_npy_legacy(gridff_path)
        
    with open(sub_xyz, 'r') as f:
        lines = f.readlines()
        comment = lines[1].strip()
        lvec = None
        for prefix in ["lvec:", "lvs"]:
            if prefix in comment:
                idx = comment.find(prefix) + len(prefix)
                parts = comment[idx:].split()
                try:
                    vals = [float(v) for v in parts if v.strip()]
                    if len(vals) >= 9:
                        lvec = np.array(vals[:9]).reshape(3,3).astype(np.float32)
                        break
                except ValueError:
                    pass
        if lvec is None:
            raise ValueError(f"Substrate lattice vectors missing in {sub_xyz}")
            
    # Use provided grid_p0 and grid_step, or calculate from lattice
    if grid_p0 is None:
        grid_p0 = (0.0, 0.0, 0.0)
    if grid_step is None:
        ax = float(np.linalg.norm(lvec[0]))
        ay = float(np.linalg.norm(lvec[1]))
        az = float(np.linalg.norm(lvec[2]))
        grid_step = (ax / grid.shape[0], ay / grid.shape[1], az / grid.shape[2])

    rbd = RigidBodyDynamics(debug=debug)
    rbd.realloc(n_bodies=n_bodies, num_atoms=1)
    rbd.enames = ['TestAtom']
    rbd.atom_types_assigned = ['TestAtom']
    
    reqs = np.array([atom_req], dtype=np.float32)
    rbd.atom_REQ = reqs
    rbd.atom_masses = np.array([atom_mass], dtype=np.float32)
    rbd.mass_physical = float(atom_mass)
    rbd.mass_trans = float(atom_mass)
    rbd.mass_rot = float(atom_mass)
    
    atom_plq_single = _reqs_to_plq(reqs, alpha=alpha_morse)
    atom_plq = np.repeat(atom_plq_single[None, :, :], n_bodies, axis=0).reshape(n_bodies, 4)
    rbd.atom_PLQ = atom_plq.copy()
    
    pos4 = np.zeros((n_bodies, 4), dtype=np.float32)
    pos4[:, :3] = scan_positions
    pos4[:, 3] = float(atom_mass)
    
    quat4 = np.zeros((n_bodies, 4), dtype=np.float32)
    quat4[:, 3] = 1.0
    zero4 = np.zeros((n_bodies, 4), dtype=np.float32)
    
    Iinv_relax = np.eye(3, dtype=np.float32)
    atom_body = np.zeros((n_bodies, 1, 3), dtype=np.float32)
    
    rbd.upload_state(pos4, quat4, zero4, zero4, rbd.mass_trans, 1.0 / rbd.mass_trans, np.repeat(Iinv_relax[None, :, :], n_bodies, axis=0), atom_body, atom_PLQ=atom_plq)
    rbd.init_gridff(grid, grid_p0=grid_p0, grid_step=grid_step)
    
    # Run 1 step with dt=0.0 to evaluate forces without moving
    rbd.run_gridff(num_steps=1, dt=0.0)
    
    outputs = rbd.download_selected(('atom_force',))
    forces = outputs['atom_force'][:, 0, :3]
    energies = outputs['atom_force'][:, 0, 3]
    
    return forces, energies

def plot_gridff_diagnostics(gridff_path, sub_xyz, ax, ay, az, iz_slices=[48], iy_slice=None, save_path='grid_diagnostics.png'):
    """
    Diagnostic tool to plot GridFF channels and overlay substrate atoms.
    """
    import matplotlib.pyplot as plt
    try:
        grid = np.load(gridff_path)
    except Exception:
        grid = _load_npy_legacy(gridff_path)
    
    nx, ny, nz, nch = grid.shape
    dx, dy, dz = ax/nx, ay/ny, az/nz
    
    def read_xyz(fname):
        with open(fname, 'r') as f:
            n = int(f.readline())
            f.readline() # comment
            pos = []
            enames = []
            for _ in range(n):
                parts = f.readline().split()
                enames.append(parts[0])
                pos.append([float(x) for x in parts[1:4]])
        return np.array(pos), enames

    sub_apos, sub_enames = read_xyz(sub_xyz)
    colors = ['purple' if e in ['Ca', 'Na'] else 'green' for e in sub_enames]
    
    n_slices = len(iz_slices)
    fig, axs = plt.subplots(nch, n_slices + 1, figsize=(5*(n_slices+1), 4*nch))
    if nch == 1: axs = axs[None, :]

    names = ['Pauli (P)', 'London (L)', 'Electrostatic (Q)', 'Hydrogen (H)'][:nch]

    for i in range(nch):
        # XY Slices
        for j, iz in enumerate(iz_slices):
            z_val = iz * dz
            im = axs[i, j].imshow(grid[:, :, iz, i].T, extent=[0, ax, 0, ay], origin='lower', cmap='bwr', aspect='equal')
            # Overlay atoms near this z
            mask = np.abs(sub_apos[:, 2] - z_val) < 1.0
            axs[i, j].scatter(sub_apos[mask, 0], sub_apos[mask, 1], c=np.array(colors)[mask], s=20, alpha=0.5)
            axs[i, j].set_title(f"{names[i]} XY at z={z_val:.2f} A")
            plt.colorbar(im, ax=axs[i, j])

        # XZ Slice
        if iy_slice is None: iy_slice = ny // 2
        y_val = iy_slice * dy
        im_xz = axs[i, -1].imshow(grid[:, iy_slice, :, i].T, extent=[0, ax, 0, az], origin='lower', cmap='bwr', aspect='equal')
        mask_y = np.abs(sub_apos[:, 1] - y_val) < 2.0
        axs[i, -1].scatter(sub_apos[mask_y, 0], sub_apos[mask_y, 2], c=np.array(colors)[mask_y], s=20, alpha=0.5)
        axs[i, -1].set_title(f"{names[i]} XZ at y={y_val:.2f} A")
        axs[i, -1].set_ylim(0, 15)
        plt.colorbar(im_xz, ax=axs[i, -1])

    plt.tight_layout()
    plt.savefig(save_path)
    print(f"Saved GridFF diagnostics to {save_path}")

class RigidBodyAFM:
    """
    High-level class for simulating AFM scanning using GPU rigid body dynamics.
    The molecule is attached via a harmonic spring to a moving "tip" (anchor point).
    """
    def __init__(self, mol_path, gridff_path, sub_xyz, type_map=None, debug=False,
                 anchor_idx=None, anchor_k=0.0, mass_trans=1.0, mass_rot=None):
        self.mol_path = mol_path
        self.gridff_path = gridff_path
        self.sub_xyz = sub_xyz
        self.type_map = type_map or {}
        self.debug = debug
        self.rbd = None
        self.anchor_idx = anchor_idx
        self.anchor_k = anchor_k
        self.mass_trans = mass_trans
        self.mass_rot = mass_rot

    def prepare(self, n_bodies=1, initial_positions=None, initial_quats=None):
        self.rbd = RigidBodyDynamics.from_xyz_and_grid(
            self.mol_path, self.gridff_path, self.sub_xyz,
            n_bodies=n_bodies,
            body_positions=initial_positions,
            quats=initial_quats,
            type_map=self.type_map,
            debug=self.debug,
            mass_trans=self.mass_trans,
            mass_rot=self.mass_rot,
        )
        self.n_bodies = n_bodies
        # Initial anchors setup
        if self.anchor_idx is not None and self.anchor_k > 0.0:
            outputs = self.rbd.download_selected(('atom_positions',))
            world_atoms = outputs['atom_positions']
            anchors = np.zeros((self.n_bodies, self.rbd.num_atoms, 4), dtype=np.float32)
            anchors[:, :, 3] = -1.0 # default no anchor
            # Set anchor for specific atom
            anchors[:, self.anchor_idx, :3] = world_atoms[:, self.anchor_idx, :3]
            anchors[:, self.anchor_idx, 3] = self.anchor_k
            self.rbd.update_anchors(anchors.reshape(self.rbd.total_atoms, 4))
    
    def set_anchor_positions(self, tip_positions):
        """
        tip_positions: (n_bodies, 3) or (3,) array of new tip locations.
        """
        tip_positions = np.asarray(tip_positions, dtype=np.float32)
        if tip_positions.ndim == 1:
            tip_positions = np.repeat(tip_positions[None, :], self.n_bodies, axis=0)
        
        # We need the current anchors to update just the coordinates of the anchor atom
        # Actually, it's simpler: just create the array from scratch, assuming we only anchor anchor_idx
        anchors = np.zeros((self.n_bodies, self.rbd.num_atoms, 4), dtype=np.float32)
        anchors[:, :, 3] = -1.0
        anchors[:, self.anchor_idx, :3] = tip_positions
        anchors[:, self.anchor_idx, 3] = self.anchor_k
        self.rbd.update_anchors(anchors.reshape(self.rbd.total_atoms, 4))

    def relax_to_constraint(self, nsteps=1000, dt=0.05, fconv=1e-3, tconv=1e-3, chunk=100,
                              lin_damp=0.92, ang_damp=0.88, force_scale=1.0, torque_scale=1.0):
        """
        Relax the system for the current tip positions.
        """
        converged = np.zeros(self.n_bodies, dtype=bool)
        
        for i in range(0, nsteps, chunk):
            nrun = min(chunk, nsteps - i)
            self.rbd.run_gridff(nrun, dt, lin_damp=lin_damp, ang_damp=ang_damp,
                                force_scale=force_scale, torque_scale=torque_scale)
            outputs = self.rbd.download_selected(('body_force', 'body_torque'))
            bf = outputs['body_force'][:, :3]
            bt = outputs['body_torque'][:, :3]
            f_norm = np.linalg.norm(bf, axis=1)
            t_norm = np.linalg.norm(bt, axis=1)
            converged = (f_norm < fconv) & (t_norm < tconv)
            if np.all(converged):
                break
        
        # Download full state
        return self.rbd.download_outputs(), converged
