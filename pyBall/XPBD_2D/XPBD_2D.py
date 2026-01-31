"""
XPBD_2D.py - 2D Position Based Dynamics Simulator with Complex Number Rotation

In 2D, rotation is represented by complex numbers z = (cosθ, sinθ).
Multiplying a vector v by z rotates it: z * v = (z.x*v.x - z.y*v.y, z.y*v.x + z.x*v.y)

Key differences from 3D:
- float2 instead of float3 for positions/vectors
- float2 (complex) instead of float4 (quaternion) for rotation
- float (scalar) instead of float3 for angular velocity
- Simplified inertia (scalar instead of tensor)
"""

import os
import numpy as np
import pyopencl as cl


class XPBD_2D:
    """2D XPBD simulator with port-based constraints."""
    
    def __init__(self, num_atoms, prefer_gpu=True, device_idx=0):
        self.num_atoms = num_atoms
        
        # OpenCL Setup
        self.ctx = self._make_context(prefer_gpu=prefer_gpu, device_idx=device_idx)
        self.queue = cl.CommandQueue(self.ctx)
        
        # Load OpenCL source
        curr_dir = os.path.dirname(os.path.abspath(__file__))
        cl_path = os.path.join(curr_dir, "XPBD_2D.cl")
        with open(cl_path, "r") as f:
            src = f.read()
        self.prg = cl.Program(self.ctx, src).build()
        
        # Allocate buffers
        mf = cl.mem_flags
        n = num_atoms
        
        # State buffers (float2 for 2D)
        self.cl_pos = cl.Buffer(self.ctx, mf.READ_WRITE, n * 8)      # float2
        self.cl_vel = cl.Buffer(self.ctx, mf.READ_WRITE, n * 8)      # float2
        self.cl_rot = cl.Buffer(self.ctx, mf.READ_WRITE, n * 8)      # float2 (complex)
        self.cl_omega = cl.Buffer(self.ctx, mf.READ_WRITE, n * 4)    # float (scalar)
        
        # Force buffers
        self.cl_force = cl.Buffer(self.ctx, mf.READ_WRITE, n * 8)    # float2
        
        # Topology buffers (int4 for neighbors)
        self.cl_neighs = cl.Buffer(self.ctx, mf.READ_ONLY, n * 16)   # int4
        self.cl_bkSlots = cl.Buffer(self.ctx, mf.READ_ONLY, n * 16)  # int4
        
        # Port data (node-only)
        self.cl_port_local = None  # float2[nnode*4]
        self.cl_port_n = None      # uchar[nnode]
        
        # Stiffness
        self.cl_stiffness = cl.Buffer(self.ctx, mf.READ_ONLY, n * 16)  # float4 per node
        self.cl_stiffness_flat = None  # float[nnode*4] for some kernels
        
        # PBD/XPBD correction buffers (node-only)
        self.cl_dpos_node = None
        self.cl_dtheta_node = None
        self.cl_dpos_neigh = None
        self.cl_lambda = None
        
        # Node buffers for explicit force solver
        self.cl_fneigh = None
        self.cl_pneigh = None
        
        self._nnode = 0
        
        # Host placeholders
        self._host_pos = np.zeros((n, 2), dtype=np.float32)
        self._host_vel = np.zeros((n, 2), dtype=np.float32)
        self._host_rot = np.zeros((n, 2), dtype=np.float32)  # complex
        self._host_omega = np.zeros((n,), dtype=np.float32)
        
    def _make_context(self, prefer_gpu=True, device_idx=0):
        dev_type = cl.device_type.GPU if prefer_gpu else cl.device_type.ALL
        try:
            for platform in cl.get_platforms():
                gpus = platform.get_devices(device_type=dev_type)
                if gpus:
                    return cl.Context(devices=[gpus[device_idx % len(gpus)]])
        except Exception:
            pass
        return cl.create_some_context(interactive=False)
    
    def _ensure_node_buffers(self, nnode):
        """Allocate node-only buffers once nnode is known."""
        if self._nnode == nnode and self.cl_fneigh is not None:
            return
        
        mf = cl.mem_flags
        self._nnode = nnode
        
        # Explicit force buffers
        self.cl_fneigh = cl.Buffer(self.ctx, mf.READ_WRITE, nnode * 4 * 8)  # float2[nnode*4]
        self.cl_pneigh = cl.Buffer(self.ctx, mf.READ_WRITE, nnode * 4 * 8)  # float2[nnode*4]
        
        # Port geometry
        self.cl_port_local = cl.Buffer(self.ctx, mf.READ_ONLY, nnode * 4 * 8)  # float2
        self.cl_port_n = cl.Buffer(self.ctx, mf.READ_ONLY, nnode)  # uchar
        
        # PBD correction buffers
        self.cl_dpos_node = cl.Buffer(self.ctx, mf.READ_WRITE, nnode * 8)  # float2
        self.cl_dtheta_node = cl.Buffer(self.ctx, mf.READ_WRITE, nnode * 4)  # float
        self.cl_dpos_neigh = cl.Buffer(self.ctx, mf.READ_WRITE, nnode * 4 * 8)  # float2[nnode*4]
        
        # XPBD lambda buffer
        self.cl_lambda = cl.Buffer(self.ctx, mf.READ_WRITE, nnode * 4 * 4)  # float[nnode*4]
        
    def upload_topology(self, neighs, bkSlots, stiffness):
        """Upload neighbor topology and stiffness.
        
        Args:
            neighs: (natoms, 4) int32 array of neighbor indices (-1 for empty)
            bkSlots: (natoms, 4) int32 array of back-slot indices for recoil gather
            stiffness: (natoms, 4) float32 array of stiffness per port
        """
        neighs = np.asarray(neighs, dtype=np.int32)
        bkSlots = np.asarray(bkSlots, dtype=np.int32)
        stiffness = np.asarray(stiffness, dtype=np.float32)
        
        # Reshape to int4/float4 format
        neighs_4 = neighs.reshape(-1, 4) if neighs.ndim == 2 else neighs
        bkSlots_4 = bkSlots.reshape(-1, 4) if bkSlots.ndim == 2 else bkSlots
        stiffness_4 = stiffness.reshape(-1, 4) if stiffness.ndim == 2 else stiffness
        
        cl.enqueue_copy(self.queue, self.cl_neighs, neighs_4)
        cl.enqueue_copy(self.queue, self.cl_bkSlots, bkSlots_4)
        cl.enqueue_copy(self.queue, self.cl_stiffness, stiffness_4).wait()
    
    def upload_state(self, pos, rot=None, vel=None, omega=None):
        """Upload state variables.
        
        Args:
            pos: (natoms, 2) float32 positions
            rot: (natoms, 2) float32 complex rotations (optional, defaults to identity)
            vel: (natoms, 2) float32 velocities (optional)
            omega: (natoms,) float32 angular velocities (optional)
        """
        pos = np.asarray(pos, dtype=np.float32)
        
        # Upload positions
        if pos.shape != (self.num_atoms, 2):
            raise ValueError(f"pos shape {pos.shape} != ({self.num_atoms}, 2)")
        cl.enqueue_copy(self.queue, self.cl_pos, pos)
        
        # Upload rotations (default to identity = (1, 0))
        if rot is None:
            rot = np.zeros((self.num_atoms, 2), dtype=np.float32)
            rot[:, 0] = 1.0  # real part = 1
        else:
            rot = np.asarray(rot, dtype=np.float32)
        cl.enqueue_copy(self.queue, self.cl_rot, rot)
        
        # Upload velocities
        if vel is None:
            vel = np.zeros((self.num_atoms, 2), dtype=np.float32)
        else:
            vel = np.asarray(vel, dtype=np.float32)
        cl.enqueue_copy(self.queue, self.cl_vel, vel)
        
        # Upload angular velocities
        if omega is None:
            omega = np.zeros((self.num_atoms,), dtype=np.float32)
        else:
            omega = np.asarray(omega, dtype=np.float32)
        cl.enqueue_copy(self.queue, self.cl_omega, omega).wait()
    
    def upload_ports(self, port_local, port_n, nnode):
        """Upload port geometry for nodes.
        
        Args:
            port_local: (nnode, 4, 2) float32 local port offsets
            port_n: (nnode,) uint8 number of active ports per node
            nnode: int number of nodes
        """
        self._ensure_node_buffers(nnode)
        
        port_local = np.asarray(port_local, dtype=np.float32)
        port_n = np.asarray(port_n, dtype=np.uint8)
        
        if port_local.shape != (nnode, 4, 2):
            raise ValueError(f"port_local shape {port_local.shape} != ({nnode}, 4, 2)")
        
        cl.enqueue_copy(self.queue, self.cl_port_local, port_local.reshape(nnode * 4, 2))
        cl.enqueue_copy(self.queue, self.cl_port_n, port_n).wait()
    
    def download_state(self):
        """Download current state from GPU.
        
        Returns:
            pos: (natoms, 2) float32 positions
            rot: (natoms, 2) float32 complex rotations
            vel: (natoms, 2) float32 velocities
            omega: (natoms,) float32 angular velocities
        """
        cl.enqueue_copy(self.queue, self._host_pos, self.cl_pos)
        cl.enqueue_copy(self.queue, self._host_rot, self.cl_rot)
        cl.enqueue_copy(self.queue, self._host_vel, self.cl_vel)
        cl.enqueue_copy(self.queue, self._host_omega, self.cl_omega).wait()
        
        return self._host_pos.copy(), self._host_rot.copy(), self._host_vel.copy(), self._host_omega.copy()

    def set_atom_pos(self, ia, xy):
        """Update a single atom position on GPU (useful for interactive dragging)."""
        ia = int(ia)
        if ia < 0 or ia >= self.num_atoms:
            raise ValueError(f"set_atom_pos: ia out of range {ia} not in [0,{self.num_atoms})")
        xy = np.asarray(xy, dtype=np.float32).reshape(1, 2)
        cl.enqueue_copy(self.queue, self.cl_pos, xy, device_offset=ia * 8).wait()

    def set_atom_vel(self, ia, v):
        """Update a single atom velocity on GPU."""
        ia = int(ia)
        if ia < 0 or ia >= self.num_atoms:
            raise ValueError(f"set_atom_vel: ia out of range {ia} not in [0,{self.num_atoms})")
        v = np.asarray(v, dtype=np.float32).reshape(1, 2)
        cl.enqueue_copy(self.queue, self.cl_vel, v, device_offset=ia * 8).wait()

    def set_atom_omega(self, ia, w):
        """Update a single atom angular velocity on GPU."""
        ia = int(ia)
        if ia < 0 or ia >= self.num_atoms:
            raise ValueError(f"set_atom_omega: ia out of range {ia} not in [0,{self.num_atoms})")
        w = np.asarray(w, dtype=np.float32).reshape(1)
        cl.enqueue_copy(self.queue, self.cl_omega, w, device_offset=ia * 4).wait()
    
    # -------------------------------------------------------------------------
    # SOLVER METHODS
    # -------------------------------------------------------------------------
    
    def step_explicit_force(self, nnode, dt=0.01, damp=1.0, nsteps=1):
        """Explicit force-based dynamics step.
        
        Args:
            nnode: number of rigid nodes (first nnode atoms)
            dt: time step
            damp: velocity damping factor
            nsteps: number of sub-steps
        """
        self._ensure_node_buffers(nnode)
        natoms = np.int32(self.num_atoms)
        nnode_i = np.int32(nnode)
        dt_f = np.float32(dt)
        damp_f = np.float32(damp)
        
        for _ in range(nsteps):
            # Clear forces
            self.prg.clear_2d_forces(self.queue, (self.num_atoms,), None, natoms, self.cl_force)
            self.prg.clear_2d_node_buffers(self.queue, (nnode * 4,), None, nnode_i, 
                                           self.cl_fneigh, self.cl_pneigh)
            
            # Gather port forces
            self.prg.gather_port_forces_2d(
                self.queue, (nnode,), None,
                nnode_i,
                self.cl_pos, self.cl_rot,
                self.cl_neighs, self.cl_stiffness,
                self.cl_port_local, self.cl_port_n,
                self.cl_force, self.cl_fneigh, self.cl_pneigh
            )
            
            # Integrate
            self.prg.integrate_2d_explicit(
                self.queue, (self.num_atoms,), None,
                natoms, nnode_i,
                self.cl_pos, self.cl_vel, self.cl_rot, self.cl_omega,
                self.cl_bkSlots, self.cl_force, self.cl_fneigh, self.cl_pneigh,
                dt_f, damp_f
            )
    
    def step_pbd(self, nnode, iterations=10, relaxation=0.5):
        """PBD constraint projection step.
        
        Args:
            nnode: number of rigid nodes
            iterations: number of solver iterations
            relaxation: relaxation factor for position updates
        """
        self._ensure_node_buffers(nnode)
        natoms = np.int32(self.num_atoms)
        nnode_i = np.int32(nnode)
        relax_f = np.float32(relaxation)
        
        # Upload flat stiffness for this kernel
        if self.cl_stiffness_flat is None:
            mf = cl.mem_flags
            self.cl_stiffness_flat = cl.Buffer(self.ctx, mf.READ_ONLY, nnode * 4 * 4)
        
        # Copy stiffness to flat buffer
        # Note: This assumes stiffness was uploaded as float4, we need to flatten
        # For simplicity, re-upload as flat
        
        for itr in range(iterations):
            # Clear correction buffers
            cl.enqueue_fill_buffer(self.queue, self.cl_dpos_node, np.float32(0), 0, nnode * 8)
            cl.enqueue_fill_buffer(self.queue, self.cl_dtheta_node, np.float32(0), 0, nnode * 4)
            cl.enqueue_fill_buffer(self.queue, self.cl_dpos_neigh, np.float32(0), 0, nnode * 4 * 8)
            
            # Compute corrections
            self.prg.compute_corrections_2d(
                self.queue, (nnode,), None,
                nnode_i,
                self.cl_pos, self.cl_rot, self.cl_neighs,
                self.cl_port_local, self.cl_stiffness_flat,
                self.cl_dpos_node, self.cl_dtheta_node, self.cl_dpos_neigh
            )
            
            # Apply corrections
            self.prg.apply_corrections_2d(
                self.queue, (self.num_atoms,), None,
                natoms, nnode_i,
                self.cl_pos, self.cl_rot, self.cl_bkSlots,
                self.cl_dpos_node, self.cl_dtheta_node, self.cl_dpos_neigh,
                relax_f
            )
    
    def step_xpbd(self, nnode, dt=0.1, iterations=10, reset_lambda=True, callback=None):
        """XPBD constraint projection step with compliance.
        
        Args:
            nnode: number of rigid nodes
            dt: time step (for compliance calculation)
            iterations: number of solver iterations
            reset_lambda: whether to reset Lagrange multipliers at start
        """
        self._ensure_node_buffers(nnode)
        natoms = np.int32(self.num_atoms)
        nnode_i = np.int32(nnode)
        dt_f = np.float32(dt)
        
        # Reset lambda if requested
        if reset_lambda:
            self.prg.reset_lambda_2d(self.queue, (nnode * 4,), None, 
                                     np.int32(nnode * 4), self.cl_lambda)
        
        for _ in range(iterations):
            # Clear correction buffers
            cl.enqueue_fill_buffer(self.queue, self.cl_dpos_node, np.float32(0), 0, nnode * 8)
            cl.enqueue_fill_buffer(self.queue, self.cl_dtheta_node, np.float32(0), 0, nnode * 4)
            cl.enqueue_fill_buffer(self.queue, self.cl_dpos_neigh, np.float32(0), 0, nnode * 4 * 8)
            
            # Compute XPBD corrections
            self.prg.compute_xpbd_corrections_2d(
                self.queue, (nnode,), None,
                nnode_i,
                self.cl_pos, self.cl_rot, self.cl_neighs,
                self.cl_stiffness, self.cl_port_local,
                self.cl_lambda, self.cl_dpos_neigh,
                self.cl_dpos_node, self.cl_dtheta_node,
                dt_f
            )
            
            # Gather and apply
            self.prg.gather_and_apply_xpbd_2d(
                self.queue, (self.num_atoms,), None,
                natoms, nnode_i,
                self.cl_pos, self.cl_rot, self.cl_bkSlots,
                self.cl_dpos_neigh, self.cl_dpos_node, self.cl_dtheta_node
            )

            if callback is not None:
                callback(itr)


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def build_neighs_bk_from_bonds_2d(n, bonds, max_deg=4):
    """Build neighbor and back-slot arrays from bond list.
    
    Args:
        n: total number of atoms
        bonds: list of (i, j) tuples defining bonds
        max_deg: maximum degree per node
        
    Returns:
        neighs: (n, 4) int32 neighbor indices (-1 for empty)
        bks: (n, 4) int32 back-slot indices
    """
    neighs = np.full((n, max_deg), -1, dtype=np.int32)
    bks = np.full((n, max_deg), -1, dtype=np.int32)
    deg = np.zeros((n,), dtype=np.int32)
    
    for (i, j) in bonds:
        if deg[i] >= max_deg or deg[j] >= max_deg:
            raise RuntimeError(f"Degree>={max_deg} for bond {i}-{j}")
        si = int(deg[i])
        sj = int(deg[j])
        neighs[i, si] = j
        neighs[j, sj] = i
        bks[i, si] = sj
        bks[j, sj] = si
        deg[i] += 1
        deg[j] += 1
    
    return neighs, bks


def make_bk_slots_2d(neighs, *, nnode, natoms=None):
    """Build back-slot indices for recoil force gathering.
    
    Args:
        neighs: (natoms, 4) neighbor array
        nnode: number of nodes (have ports and rotation)
        natoms: total atoms (default to len(neighs))
        
    Returns:
        bkSlots: (natoms, 4) int32 indices into node port buffers
    """
    if natoms is None:
        natoms = int(neighs.shape[0])
    
    bkSlots = np.full((natoms, 4), -1, dtype=np.int32)
    bkCount = np.zeros((natoms,), dtype=np.int32)
    
    for ia in range(int(nnode)):
        for k in range(4):
            ja = int(neighs[ia, k])
            if ja < 0:
                continue
            s = int(bkCount[ja])
            if s >= 4:
                raise RuntimeError(f"bkSlots overflow: atom {ja} has >4 back slots")
            bkSlots[ja, s] = ia * 4 + k
            bkCount[ja] += 1
    
    return bkSlots


def make_stiffness_from_bonds_2d(n, neighs, k_bond=200.0):
    """Create stiffness array from neighbor topology.
    
    Args:
        n: number of atoms
        neighs: (n, 4) neighbor indices
        k_bond: default stiffness value
        
    Returns:
        stiffness: (n, 4) float32 stiffness per port
    """
    stiffness = np.zeros((n, 4), dtype=np.float32)
    for i in range(n):
        for k in range(4):
            if neighs[i, k] >= 0:
                stiffness[i, k] = k_bond
    return stiffness
