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

VERBOSE = 0

def set_verbose(level:int):
    global VERBOSE
    VERBOSE = int(level)


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

        # Kernel cache dict must exist before caching kernels
        self._krn = {}
        # Cache all kernels to avoid RepeatedKernelRetrieval warnings
        self._init_cached_kernels()
        # Convenience cached kernel accessors (avoid repeated dict lookups)
        self._k_init_hb_pos_2d = self._k('init_hb_pos_2d')
        self._k_init_hb_rot_2d = self._k('init_hb_rot_2d')
        self._k_compute_corrections_2d = self._k('compute_corrections_2d')
        self._k_apply_corrections_2d = self._k('apply_corrections_2d')
        self._k_compute_collision_cluster = self._k('compute_collision_cluster')
        self._k_compute_cluster_fused_2d = self._k('compute_cluster_fused_2d') if 'compute_cluster_fused_2d' in self._krn else None
        self._k_update_bboxes_2d = self._k('update_bboxes_2d')
        self._k_build_local_topology_2d = self._k('build_local_topology_2d')
        self._k_reset_lambda_2d = self._k('reset_lambda_2d')
        self._k_compute_xpbd_corrections_2d = self._k('compute_xpbd_corrections_2d') if 'compute_xpbd_corrections_2d' in self._krn else None
        self._k_gather_and_apply_xpbd_2d = self._k('gather_and_apply_xpbd_2d') if 'gather_and_apply_xpbd_2d' in self._krn else None
        self._k_compute_velocities = self._k('compute_velocities_from_positions')
        self._k_compute_omega = self._k('compute_angular_velocities_from_rotations')
        self._k_reset_lambda_debug = self._krn.get('compute_xpbd_corrections_2d_debug', None)
        
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
        
        # Mass buffer (float2 per atom: M for translation, I for rotation)
        self.cl_mass = cl.Buffer(self.ctx, mf.READ_WRITE, n * 8)     # float2[natoms] -> (M, I)
        # Initialize with unit masses and default inertia
        mass_init = np.ones((n, 2), dtype=np.float32)
        cl.enqueue_copy(self.queue, self.cl_mass, mass_init)

        # Collision radius buffer (float per atom)
        self.cl_radius = cl.Buffer(self.ctx, mf.READ_WRITE, n * 4)   # float[natoms]
        rad_init = np.ones((n,), dtype=np.float32) * np.float32(0.5)
        cl.enqueue_copy(self.queue, self.cl_radius, rad_init)
        
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

        # Cluster collision solver buffers (lazy init)
        self._cluster_inited = False
        self.group_size = 64
        self.cluster_stride = self.group_size
        self.max_ghosts = 128
        self.num_groups = (int(self.num_atoms) + self.group_size - 1) // self.group_size
        self.cl_bboxes_min = None
        self.cl_bboxes_max = None
        self.cl_ghost_indices = None
        self.cl_ghost_counts = None
        self.cl_cluster_counts = None
        self.cl_neighs_local = None
        self.cl_prev_pos_cluster = None
        self.cl_mom_pos_cluster = None
        
        self._nnode = 0
        
        # Host placeholders
        self._host_pos = np.zeros((n, 2), dtype=np.float32)
        self._host_vel = np.zeros((n, 2), dtype=np.float32)
        self._host_rot = np.zeros((n, 2), dtype=np.float32)  # complex
        self._host_omega = np.zeros((n,), dtype=np.float32)

        self._host_radius = np.zeros((n,), dtype=np.float32)
        
    def _init_cached_kernels(self):
        """Cache all OpenCL kernels to avoid retrieval overhead."""
        kernel_names = [
            # Basic init/clear kernels
            'clear_2d_forces', 'clear_2d_node_buffers', 'init_hb_pos_2d', 'init_hb_rot_2d',
            'init_mom_pos_2d', 'init_mom_rot_2d', 'reset_lambda_2d',
            # Force-based kernels
            'gather_port_forces_2d', 'integrate_2d_explicit',
            # PBD/XPBD kernels
            'compute_corrections_2d', 'apply_corrections_2d',
            'compute_xpbd_corrections_2d', 'gather_and_apply_xpbd_2d',
            'compute_xpbd_corrections_2d_debug',
            # Cluster collision kernels
            'update_bboxes_2d', 'build_local_topology_2d', 'compute_collision_cluster', 'compute_cluster_fused_2d',
            # MD kernels
            'compute_velocities_from_positions', 'compute_angular_velocities_from_rotations',
        ]
        for name in kernel_names:
            try:
                self._krn[name] = cl.Kernel(self.prg, name)
            except cl.RuntimeError:
                # Kernel may not exist in the program (e.g., debug kernels)
                pass
        
    def _k(self, name):
        """Get cached kernel by name."""
        if name not in self._krn:
            # Lazy cache if not pre-initialized
            self._krn[name] = cl.Kernel(self.prg, name)
        return self._krn[name]
        
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
        """Allocate node-only buffers once nnode is known.
        dpos_node is sized for max(nnode, natoms) so collision corrections can write per-atom.
        """
        if self._nnode == nnode and self.cl_fneigh is not None:
            return
        
        mf = cl.mem_flags
        self._nnode = nnode
        size_nodes = max(nnode, self.num_atoms)
        
        # Explicit force buffers
        self.cl_fneigh = cl.Buffer(self.ctx, mf.READ_WRITE, nnode * 4 * 8)  # float2[nnode*4]
        self.cl_pneigh = cl.Buffer(self.ctx, mf.READ_WRITE, nnode * 4 * 8)  # float2[nnode*4]
        
        # Port geometry
        self.cl_port_local = cl.Buffer(self.ctx, mf.READ_ONLY, nnode * 4 * 8)  # float2
        self.cl_port_n = cl.Buffer(self.ctx, mf.READ_ONLY, nnode)  # uchar
        
        # PBD correction buffers
        self.cl_dpos_node = cl.Buffer(self.ctx, mf.READ_WRITE, size_nodes * 8)  # float2[max(nnode,natoms)]
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

    def upload_radius(self, radius):
        radius = np.asarray(radius, dtype=np.float32)
        if radius.shape != (self.num_atoms,):
            raise ValueError(f"upload_radius: radius shape {radius.shape} != ({self.num_atoms},)")
        cl.enqueue_copy(self.queue, self.cl_radius, radius).wait()

    def download_radius(self):
        cl.enqueue_copy(self.queue, self._host_radius, self.cl_radius).wait()
        return self._host_radius.copy()

    def download_dpos_node(self):
        """Download dpos_node buffer (collision/position corrections) from GPU.
        
        Returns:
            dpos: (max(nnode, natoms), 2) float32 position corrections
        """
        size_nodes = max(self._nnode, self.num_atoms)
        out = np.zeros((size_nodes, 2), dtype=np.float32)
        cl.enqueue_copy(self.queue, out, self.cl_dpos_node).wait()
        return out.copy()

    def _ensure_cluster_buffers(self):
        if self._cluster_inited and self.cl_bboxes_min is not None:
            return
        mf = cl.mem_flags
        ng = int(self.num_groups)
        n = int(self.num_atoms)
        self.cl_bboxes_min       = cl.Buffer(self.ctx, mf.READ_WRITE, ng * 16)  # float4[ng]
        self.cl_bboxes_max       = cl.Buffer(self.ctx, mf.READ_WRITE, ng * 16)  # float4[ng]
        self.cl_ghost_indices    = cl.Buffer(self.ctx, mf.READ_WRITE, ng * int(self.max_ghosts) * 4)  # int
        self.cl_ghost_counts     = cl.Buffer(self.ctx, mf.READ_WRITE, ng * 4)  # int
        self.cl_cluster_counts   = cl.Buffer(self.ctx, mf.READ_ONLY,  ng * 4)  # int
        self.cl_neighs_local     = cl.Buffer(self.ctx, mf.READ_WRITE, n * 16)  # int4[natoms]
        self.cl_prev_pos_cluster = cl.Buffer(self.ctx, mf.READ_WRITE, n * 8)  # float2[natoms]
        self.cl_mom_pos_cluster  = cl.Buffer(self.ctx, mf.READ_WRITE, n * 8)   # float2[natoms]
        # Default cluster_counts: sequential packing
        counts = np.empty((ng,), dtype=np.int32)
        for g in range(ng):
            start = g * int(self.group_size)
            remaining = n - start
            if remaining < 0: remaining = 0
            if remaining > self.group_size: remaining = int(self.group_size)
            counts[g] = remaining
        cl.enqueue_copy(self.queue, self.cl_cluster_counts, counts)
        cl.enqueue_fill_buffer(self.queue, self.cl_ghost_counts, np.int32(0), 0, ng * 4)
        self._cluster_inited = True

    def set_cluster_config(self, *, group_size=None, num_groups=None, cluster_counts=None, cluster_stride=None):
        """Configure cluster tiling used by cluster-local kernels.

        Args:
            group_size: local workgroup size (<= GROUP_SIZE in kernel source)
            num_groups: number of groups (defaults to ceil(natoms/group_size))
            cluster_counts: int32[ng] number of real atoms in each group (packed at start)
        """
        if group_size is not None:
            self.group_size = int(group_size)
        if cluster_stride is not None:
            self.cluster_stride = int(cluster_stride)
        else:
            self.cluster_stride = int(self.group_size)
        if num_groups is not None:
            self.num_groups = int(num_groups)
        else:
            self.num_groups = (int(self.num_atoms) + int(self.group_size) - 1) // int(self.group_size)

        # Force reinit of cluster buffers on next use
        self._cluster_inited = False
        self.cl_bboxes_min = None
        self.cl_bboxes_max = None
        self.cl_ghost_indices = None
        self.cl_ghost_counts = None
        self.cl_cluster_counts = None
        self.cl_neighs_local = None
        self.cl_prev_pos_cluster = None
        self.cl_mom_pos_cluster = None

        self._ensure_cluster_buffers()

        ng = int(self.num_groups)
        if cluster_counts is not None:
            cc = np.asarray(cluster_counts, dtype=np.int32).reshape(ng)
            cl.enqueue_copy(self.queue, self.cl_cluster_counts, cc)
        cl.enqueue_fill_buffer(self.queue, self.cl_ghost_counts, np.int32(0), 0, ng * 4).wait()

    def debug_clusters(self, *, ng_print=8, atoms_per_group=16, ghosts_per_group=16):
        """Dump cluster_counts, ghost_counts, ghost_indices, and neighs_local for inspection."""
        self._ensure_cluster_buffers()
        ng = int(self.num_groups)
        n = int(self.num_atoms)
        cc = np.zeros((ng,), dtype=np.int32)
        gc = np.zeros((ng,), dtype=np.int32)
        gi = np.zeros((ng, int(self.max_ghosts)), dtype=np.int32)
        nl = np.zeros((n, 4), dtype=np.int32)
        cl.enqueue_copy(self.queue, cc, self.cl_cluster_counts).wait()
        cl.enqueue_copy(self.queue, gc, self.cl_ghost_counts).wait()
        cl.enqueue_copy(self.queue, gi, self.cl_ghost_indices).wait()
        cl.enqueue_copy(self.queue, nl, self.cl_neighs_local).wait()
        print(f"[DEBUG] cluster_counts (ng={ng}):")
        for g in range(min(int(ng_print), ng)):
            start = g * int(self.group_size)
            gslice = gi[g, :int(min(ghosts_per_group, gc[g]))]
            print(f"  grp[{g:3d}] c_count={cc[g]:3d} g_count={gc[g]:3d} ghosts={gslice.tolist()}")
            for i in range(min(int(atoms_per_group), int(cc[g]))):
                ia = start + i
                if ia >= n:
                    break
                nl_i = nl[ia]
                print(f"    atom[{ia:4d}] nl=({nl_i[0]:4d},{nl_i[1]:4d},{nl_i[2]:4d},{nl_i[3]:4d})")

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
            self.cl_dbg_C       = cl.Buffer(self.ctx, cl.mem_flags.READ_WRITE, size=size * 5)  # C, lambda, dtheta, K, alpha
            self.cl_dbg_dpos_i  = cl.Buffer(self.ctx, cl.mem_flags.READ_WRITE, size=size2)
            self.cl_dbg_dpos_j  = cl.Buffer(self.ctx, cl.mem_flags.READ_WRITE, size=size2)
            self.cl_dbg_r_world = cl.Buffer(self.ctx, cl.mem_flags.READ_WRITE, size=size2)
            self.cl_dbg_n       = cl.Buffer(self.ctx, cl.mem_flags.READ_WRITE, size=size2)

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
        host_data    = np.zeros(total_size, dtype=np.float32)
        host_dpos_i  = np.zeros(n * 2, dtype=np.float32)
        host_dpos_j  = np.zeros(n * 2, dtype=np.float32)
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

    def set_atom_mass(self, ia, mass, inertia=None):
        """Update a single atom mass on GPU.
        
        Args:
            ia: atom index
            mass: translational mass (M), or tuple (M, I) if inertia not provided
            inertia: rotational inertia (I), optional
        """
        ia = int(ia)
        if ia < 0 or ia >= self.num_atoms:
            raise ValueError(f"set_atom_mass: ia out of range {ia} not in [0,{self.num_atoms})")
        
        if hasattr(mass, '__len__'):
            # mass is a tuple/array (M, I)
            m = np.asarray([mass[0], mass[1]], dtype=np.float32)
        elif inertia is not None:
            m = np.asarray([mass, inertia], dtype=np.float32)
        else:
            # Keep existing inertia, update only mass
            current = np.zeros(2, dtype=np.float32)
            cl.enqueue_copy(self.queue, current, self.cl_mass, device_offset=ia * 8).wait()
            m = np.asarray([mass, current[1]], dtype=np.float32)
        
        cl.enqueue_copy(self.queue, self.cl_mass, m, device_offset=ia * 8).wait()

    def get_atom_mass(self, ia):
        """Get a single atom mass from GPU.
        
        Returns:
            (M, I) tuple of (translational_mass, rotational_inertia)
        """
        ia = int(ia)
        if ia < 0 or ia >= self.num_atoms:
            raise ValueError(f"get_atom_mass: ia out of range {ia} not in [0,{self.num_atoms})")
        m = np.zeros(2, dtype=np.float32)
        cl.enqueue_copy(self.queue, m, self.cl_mass, device_offset=ia * 8).wait()
        return (float(m[0]), float(m[1]))
    
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
            self._k('clear_2d_forces')(self.queue, (self.num_atoms,), None, natoms, self.cl_force)
            self._k('clear_2d_node_buffers')(self.queue, (nnode * 4,), None, nnode_i,   self.cl_fneigh, self.cl_pneigh)
            
            # Gather port forces
            self._k('gather_port_forces_2d')(
                self.queue, (nnode,), None,
                nnode_i,
                self.cl_pos, self.cl_rot,
                self.cl_neighs, self.cl_stiffness,
                self.cl_port_local, self.cl_port_n,
                self.cl_force, self.cl_fneigh, self.cl_pneigh
            )
            
            # Integrate
            self._k('integrate_2d_explicit')(
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
            self._k('init_hb_pos_2d')(self.queue, (self.num_atoms,), None, natoms, self.cl_pos, self.cl_hb_pos)
            self._k('init_hb_rot_2d')(self.queue, (nnode,), None, nnode_i, self.cl_rot, self.cl_hb_rot)
            # Initialize momentum buffers to zero
            cl.enqueue_fill_buffer(self.queue, self.cl_mom_pos, np.float32(0.0), 0, self.num_atoms * 8)
            cl.enqueue_fill_buffer(self.queue, self.cl_mom_rot, np.float32(0.0), 0, nnode * 8)

        for itr in range(int(iterations)):
            cl.enqueue_fill_buffer(self.queue, self.cl_dpos_node, np.float32(0), 0, nnode * 8)
            cl.enqueue_fill_buffer(self.queue, self.cl_dtheta_node, np.float32(0), 0, nnode * 4)
            cl.enqueue_fill_buffer(self.queue, self.cl_dpos_neigh, np.float32(0), 0, nnode * 4 * 8)

            self._k_compute_corrections_2d(
                self.queue, (nnode,), None,
                nnode_i,
                self.cl_pos, self.cl_rot, self.cl_neighs,
                self.cl_port_local, self.cl_stiffness_flat,
                self.cl_mass,
                self.cl_dpos_node, self.cl_dtheta_node, self.cl_dpos_neigh,
                dt_f
            )

            self._k_apply_corrections_2d(
                self.queue, (natoms,), None,
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
            
            self._k_compute_corrections_2d(
                self.queue, (nnode,), None,
                nnode_i,
                self.cl_pos, self.cl_rot, self.cl_neighs,
                self.cl_port_local, self.cl_stiffness_flat,
                self.cl_mass,  # Pass mass buffer
                self.cl_dpos_node, self.cl_dtheta_node, self.cl_dpos_neigh,
                dt_f  # Pass dt for M/dt^2 calculation
            )
            
            # Apply corrections without heavy-ball mixing (that's only for relaxation)
            self._k_apply_corrections_2d(
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
        self._k_compute_velocities(
            self.queue, (self.num_atoms,), None,
            natoms, self.cl_pos, self.cl_hb_pos, self.cl_vel, dt_f, np.float32(damp)
        )
        
        # For nodes, also compute angular velocities from rotation changes
        self._k_compute_omega(
            self.queue, (nnode,), None,
            nnode_i, self.cl_rot, self.cl_hb_rot, self.cl_omega, dt_f, np.float32(damp)
        )

    def step_pbd_cluster_relax(self, *, dt=0.1, outer_iters=1, inner_iters=10, k_coll=1.0, bmix=0.8, reset_hb=True, bbox_margin=0.5, callback=None, verbose=0):
        """Cluster-local Jacobi relaxation with collisions.

        Notes:
        - This currently solves *collisions only* (no port constraints) in a tiled local-memory manner.
        - Bonded neighbors are excluded from collision resolution using the remapped `neighs_local`.
        - Uses heavy-ball mixing on positions via `prev_pos_cluster` and `mom_pos_cluster`.
        """
        self._ensure_cluster_buffers()
        natoms = int(self.num_atoms)
        ng = int(self.num_groups)
        dt_f = np.float32(dt)
        natoms_i = np.int32(self.num_atoms)
        relax_coll_f = np.float32(1.0)
        bmix_f = np.float32(bmix)

        # Reset heavy-ball state
        if reset_hb:
            cl.enqueue_copy(self.queue, self.cl_prev_pos_cluster, self.cl_pos)
            cl.enqueue_fill_buffer(self.queue, self.cl_mom_pos_cluster, np.float32(0.0), 0, natoms * 8)

        # Heuristic margin for ghost discovery: use (2*Rmax)^2
        # We compute Rmax on host (cheap enough) because radius rarely changes.
        rad = self.download_radius()
        rmax = float(np.max(rad)) if rad.size else 0.0
        margin_sq = np.float32((2.0 * rmax + float(bbox_margin)) ** 2)
        bbox_margin_f = np.float32(bbox_margin)

        # Work sizes
        global_size = (ng * int(self.group_size),)
        local_size = (int(self.group_size),)

        for it in range(int(outer_iters)):
            # 1) AABB per cluster
            self._k_update_bboxes_2d(
                self.queue, global_size, local_size,
                self.cl_pos, self.cl_radius,
                self.cl_bboxes_min, self.cl_bboxes_max,
                cl.LocalMemory(int(self.group_size) * 16),
                cl.LocalMemory(int(self.group_size) * 16),
                np.int32(natoms)
            )

            # 2) Ghost discovery + neighbor remap (local indices)
            self._k_build_local_topology_2d(
                self.queue, global_size, local_size,
                self.cl_pos,
                self.cl_bboxes_min, self.cl_bboxes_max,
                self.cl_neighs,
                self.cl_ghost_indices, self.cl_ghost_counts, self.cl_cluster_counts,
                self.cl_neighs_local,
                np.int32(natoms), np.int32(ng),
                margin_sq, bbox_margin_f,
                np.int32(self.cluster_stride)
            )

            # DEBUG: Print neighs_local mapping if verbose > 1
            if int(verbose) > 1:
                host_neighs_local = np.empty((self.num_atoms, 4), dtype=np.int32)
                cl.enqueue_copy(self.queue, host_neighs_local, self.cl_neighs_local).wait()
                print(f"[DEBUG] neighs_local (local indices including ghosts):")
                for ii in range(min(16, self.num_atoms)):
                    nl = host_neighs_local[ii]
                    print(f"  atom[{ii:3d}] -> ({nl[0]:4d}, {nl[1]:4d}, {nl[2]:4d}, {nl[3]:4d})")

            # 3) Cluster-local Jacobi collision solve - inner loop in Python for callback support
            for inner in range(int(inner_iters)):
                # Clear dpos_node over all atoms (collision writes natoms-wide)
                cl.enqueue_fill_buffer(self.queue, self.cl_dpos_node, np.float32(0), 0, self.num_atoms * 8)

                self._k_compute_collision_cluster(
                    self.queue, global_size, local_size,
                    self.cl_pos, self.cl_radius, self.cl_mass,
                    self.cl_neighs_local,
                    self.cl_ghost_indices, self.cl_ghost_counts, self.cl_cluster_counts,
                    self.cl_dpos_node,
                    np.int32(natoms),
                    dt_f, np.float32(k_coll), np.int32(self.cluster_stride)
                )

                # Apply collisions: treat all atoms as nodes so dpos_node is applied
                self._k_apply_corrections_2d(
                    self.queue, (self.num_atoms,), None,
                    natoms_i, natoms_i,  # nnode=natoms so dpos_node applies without double-weighting capping
                    self.cl_pos, self.cl_rot, self.cl_bkSlots,
                    self.cl_dpos_node, self.cl_dtheta_node, self.cl_dpos_neigh,
                    relax_coll_f, self.cl_prev_pos_cluster, self.cl_hb_rot,  # hb_rot unused here
                    self.cl_mom_pos_cluster, self.cl_mom_rot,               # mom_rot unused here
                    bmix_f, np.float32(0.0)
                )
                if callback is not None:
                    callback(it * int(inner_iters) + inner)

    def step_pbd_cluster_fused(self, nnode, *, dt=0.1, outer_iters=1, inner_iters=10, k_coll=1.0, relax=0.7, bbox_margin=0.5, callback=None, verbose=0):
        """Cluster-local fused solver: collisions + port corrections in one kernel."""
        self._ensure_node_buffers(nnode)
        self._ensure_cluster_buffers()
        if self._k_compute_cluster_fused_2d is None:
            raise RuntimeError('step_pbd_cluster_fused: compute_cluster_fused_2d kernel not available')

        natoms = int(self.num_atoms)
        ng = int(self.num_groups)
        nnode_i = np.int32(nnode)
        natoms_i = np.int32(natoms)
        dt_f = np.float32(dt)
        bbox_margin_f = np.float32(bbox_margin)

        rad = self.download_radius()
        rmax = float(np.max(rad)) if rad.size else 0.0
        margin_sq = np.float32((2.0 * rmax + float(bbox_margin)) ** 2)

        global_size = (ng * int(self.group_size),)
        local_size = (int(self.group_size),)

        for it in range(int(outer_iters)):
            self._k_update_bboxes_2d(
                self.queue, global_size, local_size,
                self.cl_pos, self.cl_radius,
                self.cl_bboxes_min, self.cl_bboxes_max,
                cl.LocalMemory(int(self.group_size) * 16),
                cl.LocalMemory(int(self.group_size) * 16),
                np.int32(natoms)
            )

            self._k_build_local_topology_2d(
                self.queue, global_size, local_size,
                self.cl_pos,
                self.cl_bboxes_min, self.cl_bboxes_max,
                self.cl_neighs,
                self.cl_ghost_indices, self.cl_ghost_counts, self.cl_cluster_counts,
                self.cl_neighs_local,
                np.int32(natoms), np.int32(ng),
                margin_sq, bbox_margin_f
            )

            for inner in range(int(inner_iters)):
                cl.enqueue_fill_buffer(self.queue, self.cl_dpos_node, np.float32(0), 0, self.num_atoms * 8)
                cl.enqueue_fill_buffer(self.queue, self.cl_dtheta_node, np.float32(0), 0, nnode * 4)
                cl.enqueue_fill_buffer(self.queue, self.cl_dpos_neigh, np.float32(0), 0, nnode * 4 * 8)

                self._k_compute_cluster_fused_2d(
                    self.queue, global_size, local_size,
                    nnode_i,
                    self.cl_pos, self.cl_rot,
                    self.cl_radius, self.cl_mass,
                    self.cl_neighs_local,
                    self.cl_ghost_indices, self.cl_ghost_counts, self.cl_cluster_counts,
                    self.cl_port_local, self.cl_stiffness_flat,
                    self.cl_dpos_node, self.cl_dtheta_node, self.cl_dpos_neigh,
                    natoms_i,
                    dt_f, np.float32(k_coll), np.int32(self.cluster_stride)
                )

                self._k_apply_corrections_2d(
                    self.queue, (self.num_atoms,), None,
                    natoms_i, nnode_i,
                    self.cl_pos, self.cl_rot, self.cl_bkSlots,
                    self.cl_dpos_node, self.cl_dtheta_node, self.cl_dpos_neigh,
                    np.float32(relax), self.cl_hb_pos, self.cl_hb_rot,
                    self.cl_mom_pos, self.cl_mom_rot,
                    np.float32(0.0), np.float32(0.0)
                )

                if callback is not None:
                    callback(it * int(inner_iters) + inner)

    def step_pbd_cluster_relax_ports(self, nnode, *, dt=0.1, outer_iters=1, inner_iters=10, port_iters=1, k_coll=1.0, relax_ports=0.7, bmix_coll=0.8, bmix_pos_ports=None, bmix_rot_ports=None, reset_hb=True, bbox_margin=0.5, callback=None, verbose=0):
        """Split-path relaxation: cluster-local collisions, then global port constraints.

        Notes:
        - Collisions are solved in tiled local memory via the cluster kernels.
        - Ports are solved using existing global Jacobi PD kernels (`compute_corrections_2d` + `apply_corrections_2d`).
        - Heavy-ball state is kept independent:
          - collisions: `prev_pos_cluster` + `mom_pos_cluster`
          - ports: `hb_pos/hb_rot` + `mom_pos/mom_rot`
        """
        self._ensure_node_buffers(nnode)
        self._ensure_cluster_buffers()

        if int(port_iters) != 1:
            raise ValueError(f"step_pbd_cluster_relax_ports: port_iters must be 1 (got {port_iters}); port iterations are controlled by inner_iters in the outer Jacobi loop")

        natoms = int(self.num_atoms)
        ng = int(self.num_groups)
        nnode_i = np.int32(nnode)
        natoms_i = np.int32(self.num_atoms)
        dt_f = np.float32(dt)
        relax_f = np.float32(relax_ports)

        if bmix_pos_ports is None:
            bmix_pos_ports = 0.0
        if bmix_rot_ports is None:
            bmix_rot_ports = 0.0
        bmix_pos_ports_f = np.float32(bmix_pos_ports)
        bmix_rot_ports_f = np.float32(bmix_rot_ports)

        # Reset heavy-ball state for both solvers
        if reset_hb:
            # Collisions HB state
            cl.enqueue_copy(self.queue, self.cl_prev_pos_cluster, self.cl_pos)
            cl.enqueue_fill_buffer(self.queue, self.cl_mom_pos_cluster, np.float32(0.0), 0, natoms * 8)
            # Ports HB state
            self._k_init_hb_pos_2d(self.queue, (self.num_atoms,), None, natoms_i, self.cl_pos, self.cl_hb_pos)
            self._k_init_hb_rot_2d(self.queue, (nnode,), None, nnode_i, self.cl_rot, self.cl_hb_rot)
            cl.enqueue_fill_buffer(self.queue, self.cl_mom_pos, np.float32(0.0), 0, self.num_atoms * 8)
            cl.enqueue_fill_buffer(self.queue, self.cl_mom_rot, np.float32(0.0), 0, nnode * 8)

        # Heuristic margin for ghost discovery: use (2*Rmax)^2
        rad = self.download_radius()
        rmax = float(np.max(rad)) if rad.size else 0.0
        margin_sq = np.float32((2.0 * rmax + float(bbox_margin)) ** 2)
        bbox_margin_f = np.float32(bbox_margin)

        # Work sizes
        global_size = (ng * int(self.group_size),)
        local_size = (int(self.group_size),)

        for it in range(int(outer_iters)):
            # 1) AABB per cluster
            self._k_update_bboxes_2d(
                self.queue, global_size, local_size,
                self.cl_pos, self.cl_radius,
                self.cl_bboxes_min, self.cl_bboxes_max,
                cl.LocalMemory(int(self.group_size) * 16),
                cl.LocalMemory(int(self.group_size) * 16),
                np.int32(natoms)
            )

            # 2) Ghost discovery + neighbor remap (local indices)
            self._k_build_local_topology_2d(
                self.queue, global_size, local_size,
                self.cl_pos,
                self.cl_bboxes_min, self.cl_bboxes_max,
                self.cl_neighs,
                self.cl_ghost_indices, self.cl_ghost_counts, self.cl_cluster_counts,
                self.cl_neighs_local,
                np.int32(natoms), np.int32(ng),
                margin_sq, bbox_margin_f
            )

            # 3) Cluster-local Jacobi collision solve - inner loop in Python for callback support
            for inner in range(int(inner_iters)):
                # Clear correction buffers (collision writes natoms; port uses nnode portions)
                cl.enqueue_fill_buffer(self.queue, self.cl_dpos_node, np.float32(0), 0, self.num_atoms * 8)
                cl.enqueue_fill_buffer(self.queue, self.cl_dtheta_node, np.float32(0), 0, nnode * 4)
                cl.enqueue_fill_buffer(self.queue, self.cl_dpos_neigh, np.float32(0), 0, nnode * 4 * 8)

                self._k_compute_collision_cluster(
                    self.queue, global_size, local_size,
                    self.cl_pos, self.cl_radius, self.cl_mass,
                    self.cl_neighs_local,
                    self.cl_ghost_indices, self.cl_ghost_counts, self.cl_cluster_counts,
                    self.cl_dpos_node,
                    np.int32(natoms),
                    dt_f, np.float32(k_coll)
                )

                self._k_compute_corrections_2d(
                    self.queue, (nnode,), None,
                    nnode_i,
                    self.cl_pos, self.cl_rot, self.cl_neighs,
                    self.cl_port_local, self.cl_stiffness_flat,
                    self.cl_mass,
                    self.cl_dpos_node, self.cl_dtheta_node, self.cl_dpos_neigh,
                    dt_f
                )

                self._k_apply_corrections_2d(
                    self.queue, (self.num_atoms,), None,
                    natoms_i, nnode_i,
                    self.cl_pos, self.cl_rot, self.cl_bkSlots,
                    self.cl_dpos_node, self.cl_dtheta_node, self.cl_dpos_neigh,
                    relax_f, self.cl_hb_pos, self.cl_hb_rot,
                    self.cl_mom_pos, self.cl_mom_rot,
                    bmix_pos_ports_f, bmix_rot_ports_f
                )

                if callback is not None:
                    callback(it * int(inner_iters) + inner)

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
            self._k_reset_lambda_2d(self.queue, (nnode * 4,), None, 
                                    np.int32(nnode * 4), self.cl_lambda)

        if reset_hb:
            self._k_init_hb_pos_2d(self.queue, (self.num_atoms,), None, natoms, self.cl_pos, self.cl_hb_pos)
            self._k_init_hb_rot_2d(self.queue, (nnode,), None, nnode_i, self.cl_rot, self.cl_hb_rot)
            # Initialize momentum buffers to zero
            cl.enqueue_fill_buffer(self.queue, self.cl_mom_pos, np.float32(0.0), 0, self.num_atoms * 8)
            cl.enqueue_fill_buffer(self.queue, self.cl_mom_rot, np.float32(0.0), 0, nnode * 8)

        for itr in range(int(iterations)):
            # Clear correction buffers
            cl.enqueue_fill_buffer(self.queue, self.cl_dpos_node, np.float32(0), 0, nnode * 8)
            cl.enqueue_fill_buffer(self.queue, self.cl_dtheta_node, np.float32(0), 0, nnode * 4)
            cl.enqueue_fill_buffer(self.queue, self.cl_dpos_neigh, np.float32(0), 0, nnode * 4 * 8)
            
            # Compute XPBD corrections
            self._k_compute_xpbd_corrections_2d(
                self.queue, (nnode,), None,
                nnode_i,
                self.cl_pos, self.cl_rot, self.cl_neighs,
                self.cl_stiffness, self.cl_port_local,
                self.cl_mass,
                self.cl_lambda, self.cl_dpos_neigh,
                self.cl_dpos_node, self.cl_dtheta_node,
                dt_f
            )
            
            # Gather and apply with momentum buffers
            self._k_gather_and_apply_xpbd_2d(
                self.queue, (natoms,), None,
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
            self._k_reset_lambda_2d(self.queue, (nnode * 4,), None,
                                     np.int32(nnode * 4), self.cl_lambda)

        if reset_hb:
            self._k_init_hb_pos_2d(self.queue, (self.num_atoms,), None, natoms, self.cl_pos, self.cl_hb_pos)
            self._k_init_hb_rot_2d(self.queue, (nnode,), None, nnode_i, self.cl_rot, self.cl_hb_rot)
            cl.enqueue_fill_buffer(self.queue, self.cl_mom_pos, np.float32(0.0), 0, self.num_atoms * 8)
            cl.enqueue_fill_buffer(self.queue, self.cl_mom_rot, np.float32(0.0), 0, nnode * 8)

        for itr in range(int(iterations)):
            # Clear correction buffers
            cl.enqueue_fill_buffer(self.queue, self.cl_dpos_node, np.float32(0), 0, nnode * 8)
            cl.enqueue_fill_buffer(self.queue, self.cl_dtheta_node, np.float32(0), 0, nnode * 4)
            cl.enqueue_fill_buffer(self.queue, self.cl_dpos_neigh, np.float32(0), 0, nnode * 4 * 8)

            # Compute XPBD corrections with DEBUG - write to buffers instead of printf
            self._k_compute_xpbd_corrections_2d_debug(
                self.queue, (nnode,), None,
                nnode_i,
                self.cl_pos, self.cl_rot, self.cl_neighs,
                self.cl_stiffness_flat, self.cl_port_local, self.cl_mass,
                self.cl_lambda,
                self.cl_dpos_neigh, self.cl_dpos_node, self.cl_dtheta_node,
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
