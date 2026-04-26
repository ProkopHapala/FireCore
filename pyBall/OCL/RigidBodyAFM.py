import os
import numpy as np
from .RigidBodyDynamics import RigidBodyDynamics

class RigidBodyAFM:
    """
    High-level class for simulating AFM scanning using GPU rigid body dynamics.
    The molecule is attached via a harmonic spring to a moving "tip" (anchor point).
    """
    def __init__(self, mol_path, gridff_path, sub_xyz, type_map=None, debug=False,
                 anchor_idx=None, anchor_k=0.0):
        self.mol_path = mol_path
        self.gridff_path = gridff_path
        self.sub_xyz = sub_xyz
        self.type_map = type_map or {}
        self.debug = debug
        self.rbd = None
        self.anchor_idx = anchor_idx
        self.anchor_k = anchor_k

    def prepare(self, n_bodies=1, initial_positions=None, initial_quats=None):
        self.rbd = RigidBodyDynamics.from_xyz_and_grid(
            self.mol_path, self.gridff_path, self.sub_xyz,
            n_bodies=n_bodies,
            body_positions=initial_positions,
            quats=initial_quats,
            type_map=self.type_map,
            debug=self.debug
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

    def relax_to_constraint(self, nsteps=1000, dt=0.05, fconv=1e-3, tconv=1e-3, chunk=100):
        """
        Relax the system for the current tip positions.
        """
        converged = np.zeros(self.n_bodies, dtype=bool)
        
        for i in range(0, nsteps, chunk):
            nrun = min(chunk, nsteps - i)
            self.rbd.run_gridff(nrun, dt)
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
