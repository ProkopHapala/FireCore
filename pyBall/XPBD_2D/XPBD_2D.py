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
        
        # Mass buffer (float per atom)
        self.cl_mass = cl.Buffer(self.ctx, mf.READ_WRITE, n * 4)     # float[natoms]
        # Initialize with unit masses
        cl.enqueue_fill_buffer(self.queue, self.cl_mass, np.float32(1.0), 0, n * 4)
        
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

        # Heavy-ball buffers
        self.cl_hb_pos = None   # float2[natoms]
        self.cl_hb_rot = None   # float2[nnode]
        
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

        # Heavy-ball buffers
        self.cl_hb_pos = cl.Buffer(self.ctx, mf.READ_WRITE, self.num_atoms * 8)  # float2[natoms]
        self.cl_hb_rot = cl.Buffer(self.ctx, mf.READ_WRITE, nnode * 8)  # float2[nnode]
        
        # Momentum buffers (for C++ style heavy-ball mixing)
        self.cl_mom_pos = cl.Buffer(self.ctx, mf.READ_WRITE, self.num_atoms * 8)  # float2[natoms]
        self.cl_mom_rot = cl.Buffer(self.ctx, mf.READ_WRITE, nnode * 8)  # float2[nnode]
        
    def upload_topology(self, neighs, bkSlots, stiffness, nnode=None):
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

        if nnode is not None:
            nnode_i = int(nnode)
            if nnode_i < 0 or nnode_i > self.num_atoms:
                raise ValueError(f"upload_topology: nnode out of range {nnode_i} not in [0,{self.num_atoms}]")
            self._ensure_node_buffers(nnode_i)
            if self.cl_stiffness_flat is None or getattr(self, "_stiffness_flat_nnode", None) != nnode_i:
                mf = cl.mem_flags
                self.cl_stiffness_flat = cl.Buffer(self.ctx, mf.READ_ONLY, nnode_i * 4 * 4)
                self._stiffness_flat_nnode = nnode_i
            stiff_flat = np.asarray(stiffness_4[:nnode_i], dtype=np.float32).reshape(nnode_i * 4)
            cl.enqueue_copy(self.queue, self.cl_stiffness_flat, stiff_flat).wait()
    
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

    def download_lambda(self, nnode):
        """Download XPBD lambda buffer from GPU for diagnostics.
        
        Args:
            nnode: number of nodes (determines buffer size)
            
        Returns:
            lambda_arr: (nnode, 4) float32 Lagrange multipliers per port
        """
        if self.cl_lambda is None or nnode <= 0:
            return np.zeros((nnode, 4), dtype=np.float32)
        host_lambda = np.zeros((nnode * 4,), dtype=np.float32)
        cl.enqueue_copy(self.queue, host_lambda, self.cl_lambda).wait()
        return host_lambda.reshape(nnode, 4).copy()

    def _ensure_constraint_debug_buffers(self, nnode):
        """Allocate buffers for detailed constraint diagnostics per iteration."""
        size = nnode * 4 * 4  # nnode * MAX_DEGREE * sizeof(float) for scalar values
        size2 = nnode * 4 * 8  # nnode * MAX_DEGREE * sizeof(float2) for float2 values
        if getattr(self, 'cl_dbg_C', None) is None or self.cl_dbg_C.size < size:
            self.cl_dbg_C = cl.Buffer(self.ctx, cl.mem_flags.READ_WRITE, size=size * 5)  # C, lambda, dtheta, K, alpha
            self.cl_dbg_dpos_i = cl.Buffer(self.ctx, cl.mem_flags.READ_WRITE, size=size2)
            self.cl_dbg_dpos_j = cl.Buffer(self.ctx, cl.mem_flags.READ_WRITE, size=size2)
            self.cl_dbg_r_world = cl.Buffer(self.ctx, cl.mem_flags.READ_WRITE, size=size2)
            self.cl_dbg_n = cl.Buffer(self.ctx, cl.mem_flags.READ_WRITE, size=size2)

    def download_constraint_diagnostics(self, nnode):
        """Download detailed per-constraint diagnostic data.

        Returns dict with arrays shaped (nnode, 4) or (nnode, 4, 2):
            C: constraint violation
            lambda: accumulated Lagrange multiplier
            dtheta: angular correction
            K: stiffness
            alpha: compliance
            dpos_i: position correction on node i (float2)
            dpos_j: position correction on neighbor j (float2)
            r_world: rotated port vector (float2)
            n: constraint normal (float2)
        """
        n = nnode * 4
        total_size = n * 5  # 5 scalar fields
        host_data = np.zeros(total_size, dtype=np.float32)
        host_dpos_i = np.zeros(n * 2, dtype=np.float32)
        host_dpos_j = np.zeros(n * 2, dtype=np.float32)
        host_r_world = np.zeros(n * 2, dtype=np.float32)
        host_n = np.zeros(n * 2, dtype=np.float32)

        # Read all diagnostic buffers
        cl.enqueue_copy(self.queue, host_data, self.cl_dbg_C).wait()
        cl.enqueue_copy(self.queue, host_dpos_i, self.cl_dbg_dpos_i).wait()
        cl.enqueue_copy(self.queue, host_dpos_j, self.cl_dbg_dpos_j).wait()
        cl.enqueue_copy(self.queue, host_r_world, self.cl_dbg_r_world).wait()
        cl.enqueue_copy(self.queue, host_n, self.cl_dbg_n).wait()

        # Split scalar data into separate arrays
        host_C = host_data[0:n]
        host_lambda = host_data[n:2*n]
        host_dtheta = host_data[2*n:3*n]
        host_K = host_data[3*n:4*n]
        host_alpha = host_data[4*n:5*n]

        def reshape_scalar(arr):
            return arr.reshape(nnode, 4).copy()

        def reshape_float2(arr):
            return arr.reshape(nnode, 4, 2).copy()

        return {
            'C': reshape_scalar(host_C),
            'lambda': reshape_scalar(host_lambda),
            'dtheta': reshape_scalar(host_dtheta),
            'K': reshape_scalar(host_K),
            'alpha': reshape_scalar(host_alpha),
            'dpos_i': reshape_float2(host_dpos_i),
            'dpos_j': reshape_float2(host_dpos_j),
            'r_world': reshape_float2(host_r_world),
            'n': reshape_float2(host_n),
        }

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

    def set_atom_mass(self, ia, mass):
        """Update a single atom mass on GPU."""
        ia = int(ia)
        if ia < 0 or ia >= self.num_atoms:
            raise ValueError(f"set_atom_mass: ia out of range {ia} not in [0,{self.num_atoms})")
        m = np.asarray(mass, dtype=np.float32).reshape(1)
        cl.enqueue_copy(self.queue, self.cl_mass, m, device_offset=ia * 4).wait()

    def get_atom_mass(self, ia):
        """Get a single atom mass from GPU."""
        ia = int(ia)
        if ia < 0 or ia >= self.num_atoms:
            raise ValueError(f"get_atom_mass: ia out of range {ia} not in [0,{self.num_atoms})")
        m = np.zeros(1, dtype=np.float32)
        cl.enqueue_copy(self.queue, m, self.cl_mass, device_offset=ia * 4).wait()
        return float(m[0])
    
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
    
    def step_pbd(self, nnode, dt=0.1, iterations=10, relax=0.7, bmix_pos=None, bmix_rot=None, bmix=0.0, reset_hb=True, callback=None):
        """PBD constraint projection step with mass-based diagonal term and proper momentum buffers.
        
        Uses Projective Dynamics formula: xi_cor = (xi_pred * a + sum_j Kij * xj) / (a + sum_j Kij)
        where a = M/dt^2. This allows pinning atoms by setting their mass very high.
        
        Heavy-ball momentum mixing uses explicit momentum buffers following C++ reference style:
        - p_next = p_corr + d_prev * bmix
        - d_next = p_next - p_prev
        
        Args:
            nnode: number of rigid nodes
            dt: time step (for mass diagonal term M/dt^2)
            iterations: number of constraint solver iterations
            relax: relaxation factor for constraint projection
            bmix_pos: heavy-ball mixing factor for position (0 disables)
            bmix_rot: heavy-ball mixing factor for rotation (0 disables)
            reset_hb: whether to reset heavy-ball and momentum state at start
            callback: optional callback(itr) called after each inner iteration
        """
        self._ensure_node_buffers(nnode)
        natoms = np.int32(self.num_atoms)
        nnode_i = np.int32(nnode)
        relax_f = np.float32(relax)
        # Backward compatibility: if caller used legacy bmix, apply to both pos/rot unless explicit
        if bmix_pos is None:
            bmix_pos = bmix
        if bmix_rot is None:
            bmix_rot = bmix
        bmix_pos_f = np.float32(bmix_pos)
        bmix_rot_f = np.float32(bmix_rot)
        dt_f = np.float32(dt)

        if reset_hb:
            self.prg.init_hb_pos_2d(self.queue, (self.num_atoms,), None, natoms, self.cl_pos, self.cl_hb_pos)
            self.prg.init_hb_rot_2d(self.queue, (nnode,), None, nnode_i, self.cl_rot, self.cl_hb_rot)
            # Initialize momentum buffers to zero
            cl.enqueue_fill_buffer(self.queue, self.cl_mom_pos, np.float32(0.0), 0, self.num_atoms * 8)
            cl.enqueue_fill_buffer(self.queue, self.cl_mom_rot, np.float32(0.0), 0, nnode * 8)

        for itr in range(int(iterations)):
            cl.enqueue_fill_buffer(self.queue, self.cl_dpos_node, np.float32(0), 0, nnode * 8)
            cl.enqueue_fill_buffer(self.queue, self.cl_dtheta_node, np.float32(0), 0, nnode * 4)
            cl.enqueue_fill_buffer(self.queue, self.cl_dpos_neigh, np.float32(0), 0, nnode * 4 * 8)

            self.prg.compute_corrections_2d(
                self.queue, (nnode,), None,
                nnode_i,
                self.cl_pos, self.cl_rot, self.cl_neighs,
                self.cl_port_local, self.cl_stiffness_flat,
                self.cl_mass,
                self.cl_dpos_node, self.cl_dtheta_node, self.cl_dpos_neigh,
                dt_f
            )

            self.prg.apply_corrections_2d(
                self.queue, (self.num_atoms,), None,
                natoms, nnode_i,
                self.cl_pos, self.cl_rot, self.cl_bkSlots,
                self.cl_dpos_node, self.cl_dtheta_node, self.cl_dpos_neigh,
                relax_f, self.cl_hb_pos, self.cl_hb_rot,
                self.cl_mom_pos, self.cl_mom_rot,
                bmix_pos_f, bmix_rot_f
            )

            if callback is not None:
                callback(itr)
    
    def step_pbd_md(self, nnode, dt=0.01, iterations=10, relax=0.7, damp=0.98, callback=None):
        """PBD molecular dynamics step with real velocity integration.
        
        Unlike pbd_relax which uses heavy-ball mixing for convergence acceleration,
        this method computes real physical velocities: v = (x_new - x_old)/dt
        
        Args:
            nnode: number of rigid nodes
            dt: time step for velocity integration
            iterations: number of constraint solver iterations
            relax: relaxation factor for constraint projection
            damp: velocity damping factor
            callback: optional callback(itr) called after each inner iteration
        """
        self._ensure_node_buffers(nnode)
        natoms = np.int32(self.num_atoms)
        nnode_i = np.int32(nnode)
        relax_f = np.float32(relax)
        dt_f = np.float32(dt)
        
        # Save old positions for velocity computation
        cl.enqueue_copy(self.queue, self.cl_hb_pos, self.cl_pos)
        # Only copy first nnode rotations (hb_rot is sized for nnode, rot is sized for natoms)
        cl.enqueue_copy(self.queue, self.cl_hb_rot, self.cl_rot, byte_count=nnode*8)
        
        # Constraint solver iterations (Jacobi-style)
        for itr in range(int(iterations)):
            cl.enqueue_fill_buffer(self.queue, self.cl_dpos_node, np.float32(0), 0, nnode * 8)
            cl.enqueue_fill_buffer(self.queue, self.cl_dtheta_node, np.float32(0), 0, nnode * 4)
            cl.enqueue_fill_buffer(self.queue, self.cl_dpos_neigh, np.float32(0), 0, nnode * 4 * 8)
            
            self.prg.compute_corrections_2d(
                self.queue, (nnode,), None,
                nnode_i,
                self.cl_pos, self.cl_rot, self.cl_neighs,
                self.cl_port_local, self.cl_stiffness_flat,
                self.cl_mass,  # Pass mass buffer
                self.cl_dpos_node, self.cl_dtheta_node, self.cl_dpos_neigh,
                dt_f  # Pass dt for M/dt^2 calculation
            )
            
            # Apply corrections without heavy-ball mixing (that's only for relaxation)
            self.prg.apply_corrections_2d(
                self.queue, (self.num_atoms,), None,
                natoms, nnode_i,
                self.cl_pos, self.cl_rot, self.cl_bkSlots,
                self.cl_dpos_node, self.cl_dtheta_node, self.cl_dpos_neigh,
                relax_f, self.cl_hb_pos, self.cl_hb_rot,
                self.cl_mom_pos, self.cl_mom_rot,
                np.float32(0.0), np.float32(0.0)  # bmix=0 for MD
            )
            
            if callback is not None:
                callback(itr)
        
        # Compute real velocities: v = (x_new - x_old)/dt
        self.prg.compute_velocities_from_positions(
            self.queue, (self.num_atoms,), None,
            natoms, self.cl_pos, self.cl_hb_pos, self.cl_vel, dt_f, np.float32(damp)
        )
        
        # For nodes, also compute angular velocities from rotation changes
        self.prg.compute_angular_velocities_from_rotations(
            self.queue, (nnode,), None,
            nnode_i, self.cl_rot, self.cl_hb_rot, self.cl_omega, dt_f, np.float32(damp)
        )

    def step_xpbd(self, nnode, dt=0.1, iterations=10, reset_lambda=True, bmix_pos=None, bmix_rot=None, bmix=0.0, reset_hb=True, callback=None):
        """XPBD constraint projection step with compliance.
        
        Args:
            nnode: number of rigid nodes
            dt: time step (for compliance calculation)
            iterations: number of solver iterations
            reset_lambda: whether to reset Lagrange multipliers at start
            bmix_pos: heavy-ball mixing factor for position (0 disables)
            bmix_rot: heavy-ball mixing factor for rotation (0 disables)
            reset_hb: whether to reset heavy-ball and momentum state at start
            callback: optional callback(itr) called after each inner iteration
        """
        self._ensure_node_buffers(nnode)
        natoms = np.int32(self.num_atoms)
        nnode_i = np.int32(nnode)
        dt_f = np.float32(dt)
        bmix_pos_f = np.float32(bmix_pos)
        bmix_rot_f = np.float32(bmix_rot)

        # Reset lambda if requested
        if reset_lambda:
            self.prg.reset_lambda_2d(self.queue, (nnode * 4,), None, 
                                     np.int32(nnode * 4), self.cl_lambda)

        if reset_hb:
            self.prg.init_hb_pos_2d(self.queue, (self.num_atoms,), None, natoms, self.cl_pos, self.cl_hb_pos)
            self.prg.init_hb_rot_2d(self.queue, (nnode,), None, nnode_i, self.cl_rot, self.cl_hb_rot)
            # Initialize momentum buffers to zero
            cl.enqueue_fill_buffer(self.queue, self.cl_mom_pos, np.float32(0.0), 0, self.num_atoms * 8)
            cl.enqueue_fill_buffer(self.queue, self.cl_mom_rot, np.float32(0.0), 0, nnode * 8)

        for itr in range(int(iterations)):
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
            
            # Gather and apply with momentum buffers
            self.prg.gather_and_apply_xpbd_2d(
                self.queue, (self.num_atoms,), None,
                natoms, nnode_i,
                self.cl_pos, self.cl_rot, self.cl_bkSlots,
                self.cl_dpos_neigh, self.cl_dpos_node, self.cl_dtheta_node,
                self.cl_hb_pos, self.cl_hb_rot,
                self.cl_mom_pos, self.cl_mom_rot,
                bmix_pos_f, bmix_rot_f
            )

            if callback is not None:
                callback(itr)

    def step_xpbd_debug(self, nnode, dt=0.1, iterations=10, reset_lambda=True, bmix_pos=None, bmix_rot=None, bmix=0.0, reset_hb=True, max_debug_steps=5):
        """XPBD constraint projection step with DEBUG output to buffers for diagnostics.

        Args:
            nnode: number of rigid nodes
            dt: time step (for compliance calculation)
            iterations: number of solver iterations
            reset_lambda: whether to reset Lagrange multipliers at start
            bmix: heavy-ball mixing factor (0 disables)
            reset_hb: whether to reset heavy-ball state at start
            max_debug_steps: number of steps to capture detailed diagnostics
        """
        self._ensure_node_buffers(nnode)
        self._ensure_constraint_debug_buffers(nnode)
        natoms = np.int32(self.num_atoms)
        nnode_i = np.int32(nnode)
        dt_f = np.float32(dt)
        if bmix_pos is None:
            bmix_pos = bmix
        if bmix_rot is None:
            bmix_rot = bmix
        bmix_pos_f = np.float32(bmix_pos)
        bmix_rot_f = np.float32(bmix_rot)

        if reset_lambda:
            self.prg.reset_lambda_2d(self.queue, (nnode * 4,), None,
                                     np.int32(nnode * 4), self.cl_lambda)

        if reset_hb:
            self.prg.init_hb_pos_2d(self.queue, (self.num_atoms,), None, natoms, self.cl_pos, self.cl_hb_pos)
            self.prg.init_hb_rot_2d(self.queue, (nnode,), None, nnode_i, self.cl_rot, self.cl_hb_rot)
            cl.enqueue_fill_buffer(self.queue, self.cl_mom_pos, np.float32(0.0), 0, self.num_atoms * 8)
            cl.enqueue_fill_buffer(self.queue, self.cl_mom_rot, np.float32(0.0), 0, nnode * 8)

        for itr in range(int(iterations)):
            # Clear correction buffers
            cl.enqueue_fill_buffer(self.queue, self.cl_dpos_node, np.float32(0), 0, nnode * 8)
            cl.enqueue_fill_buffer(self.queue, self.cl_dtheta_node, np.float32(0), 0, nnode * 4)
            cl.enqueue_fill_buffer(self.queue, self.cl_dpos_neigh, np.float32(0), 0, nnode * 4 * 8)

            # Compute XPBD corrections with DEBUG - write to buffers instead of printf
            self.prg.compute_xpbd_corrections_2d_debug(
                self.queue, (nnode,), None,
                nnode_i,
                self.cl_pos, self.cl_rot, self.cl_neighs,
                self.cl_stiffness, self.cl_port_local,
                self.cl_lambda, self.cl_dpos_neigh,
                self.cl_dpos_node, self.cl_dtheta_node,
                dt_f,
                np.int32(itr),
                np.int32(max_debug_steps),
                self.cl_dbg_C,
                self.cl_dbg_dpos_i,
                self.cl_dbg_dpos_j,
                self.cl_dbg_r_world,
                self.cl_dbg_n
            )

            # Gather and apply
            self.prg.gather_and_apply_xpbd_2d(
                self.queue, (self.num_atoms,), None,
                natoms, nnode_i,
                self.cl_pos, self.cl_rot, self.cl_bkSlots,
                self.cl_dpos_neigh, self.cl_dpos_node, self.cl_dtheta_node,
                self.cl_hb_pos, self.cl_hb_rot,
                self.cl_mom_pos, self.cl_mom_rot,
                bmix_pos_f, bmix_rot_f
            )


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
