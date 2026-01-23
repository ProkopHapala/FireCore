import os
import numpy as np
import pyopencl as cl


class XPDB_new:
    """Simplified Position Based Dynamics with Tiled Jacobi Solver"""
    
    def __init__(self, num_atoms, group_size=64, prefer_gpu=True, device_idx=0):
        self.num_atoms = num_atoms
        self.group_size = group_size
        self.num_groups = (num_atoms + group_size - 1) // group_size
        self.max_ghosts = 128
        self.n_max_bonded = 16  # Must match CL define

        # OpenCL Setup
        self.ctx = self._make_context(prefer_gpu=prefer_gpu, device_idx=device_idx)
        self.queue = cl.CommandQueue(self.ctx)

        # Load Source
        curr_dir = os.path.dirname(os.path.abspath(__file__))
        with open(os.path.join(curr_dir, "XPDB_new.cl"), "r") as f:
            src = f.read()
        self.prg = cl.Program(self.ctx, src).build()

        # Buffers
        mf = cl.mem_flags

        # Atom Data
        self.cl_pos = cl.Buffer(self.ctx, mf.READ_WRITE, num_atoms * 16)  # float4
        self.cl_pred_pos = cl.Buffer(self.ctx, mf.READ_WRITE, num_atoms * 16)  # float4
        self.cl_vel = cl.Buffer(self.ctx, mf.READ_WRITE, num_atoms * 16)  # float4 (kept for parity/debug)
        self.cl_params = cl.Buffer(self.ctx, mf.READ_ONLY, num_atoms * 16)  # .x=rad, .w=mass

        # Bounding Boxes
        self.cl_bboxes_min = cl.Buffer(self.ctx, mf.READ_WRITE, self.num_groups * 16)
        self.cl_bboxes_max = cl.Buffer(self.ctx, mf.READ_WRITE, self.num_groups * 16)

        # Ghost Data
        self.cl_ghost_indices = cl.Buffer(self.ctx, mf.READ_WRITE, self.num_groups * self.max_ghosts * 4)
        self.cl_ghost_counts = cl.Buffer(self.ctx, mf.READ_WRITE, self.num_groups * 4)

        # Fixed-size Bond Buffers (keep global topology stable, rebuild local each frame)
        self.cl_bond_indices_global = cl.Buffer(self.ctx, mf.READ_WRITE, num_atoms * self.n_max_bonded * 4)  # int
        self.cl_bond_indices_local  = cl.Buffer(self.ctx, mf.READ_WRITE, num_atoms * self.n_max_bonded * 4)  # int
        self.cl_bond_lengths = cl.Buffer(self.ctx, mf.READ_WRITE, num_atoms * self.n_max_bonded * 4)  # float
        self.cl_bond_stiffness = cl.Buffer(self.ctx, mf.READ_WRITE, num_atoms * self.n_max_bonded * 4)  # float

        # Debug Buffers
        self.cl_debug_force_bond = cl.Buffer(self.ctx, mf.WRITE_ONLY, num_atoms * 16)
        self.cl_debug_force_coll = cl.Buffer(self.ctx, mf.WRITE_ONLY, num_atoms * 16)

        # Host placeholders
        self._host_pos = np.zeros((self.num_atoms, 4), dtype=np.float32)

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

    def upload_data(self, pos, vel, radius, mass):
        """Upload atom positions, velocities, radii, and masses"""
        # Pack params: x=radius, w=mass
        params = np.zeros((self.num_atoms, 4), dtype=np.float32)
        params[:, 0] = radius
        params[:, 3] = mass

        # Pack positions float4
        pos4 = np.zeros((self.num_atoms, 4), dtype=np.float32)
        pos4[:, :3] = pos

        # Pack velocity float4
        vel4 = np.zeros((self.num_atoms, 4), dtype=np.float32)
        vel4[:, :3] = vel

        cl.enqueue_copy(self.queue, self.cl_pos, pos4)
        cl.enqueue_copy(self.queue, self.cl_pred_pos, pos4)
        cl.enqueue_copy(self.queue, self.cl_vel, vel4)
        cl.enqueue_copy(self.queue, self.cl_params, params)

    def upload_bonds_fixed(self, bonds_adj):
        """
        Upload bond topology using fixed-size buffers.
        bonds_adj: list of lists, each element is [(neighbor_idx, length, stiffness), ...]
        """
        bond_indices = np.full((self.num_atoms, self.n_max_bonded), -1, dtype=np.int32)
        bond_lengths = np.zeros((self.num_atoms, self.n_max_bonded), dtype=np.float32)
        bond_stiffness = np.zeros((self.num_atoms, self.n_max_bonded), dtype=np.float32)

        for i, b_list in enumerate(bonds_adj):
            for k, b in enumerate(b_list[:self.n_max_bonded]):
                bond_indices[i, k] = b[0]
                bond_lengths[i, k] = b[1]
                bond_stiffness[i, k] = b[2]

        cl.enqueue_copy(self.queue, self.cl_bond_indices_global, bond_indices)
        cl.enqueue_copy(self.queue, self.cl_bond_lengths, bond_lengths)
        cl.enqueue_copy(self.queue, self.cl_bond_stiffness, bond_stiffness)

    def upload_bonds_fixed_from_arrays(self, bond_indices, bond_lengths, bond_stiffness):
        """
        Upload prebuilt bond index/length/stiffness arrays directly to global buffers.
        bond_indices: (num_atoms, n_max_bonded) int array, -1 for unused
        bond_lengths: (num_atoms, n_max_bonded) float array
        bond_stiffness: (num_atoms, n_max_bonded) float array
        """
        cl.enqueue_copy(self.queue, self.cl_bond_indices_global, bond_indices)
        cl.enqueue_copy(self.queue, self.cl_bond_lengths, bond_lengths)
        cl.enqueue_copy(self.queue, self.cl_bond_stiffness, bond_stiffness)

    def update_bboxes(self, margin_sq=4.0):
        """Update bounding boxes for all clusters"""
        loc_mem = cl.LocalMemory(self.group_size * 16)  # float4

        self.prg.update_bboxes(
            self.queue, (self.num_groups * self.group_size,), (self.group_size,),
            self.cl_pos, self.cl_params,
            self.cl_bboxes_min, self.cl_bboxes_max,
            loc_mem, loc_mem,
            np.int32(self.num_atoms)
        )

    def build_local_topology(self, Rmax, coll_scale=2.0, bbox_scale=2.0):
        """Build ghost lists and re-index bonds for tiled solver"""
        coll_margin = coll_scale * Rmax
        bbox_margin = bbox_scale * Rmax
        
        # Ensure bboxes are up to date
        self.update_bboxes(margin_sq=coll_margin**2)
        
        # Global size must be total threads: num_groups * group_size
        self.prg.build_local_topology(
            self.queue, (self.num_groups * self.group_size,), (self.group_size,),
            self.cl_pos, self.cl_bboxes_min, self.cl_bboxes_max, self.cl_bond_indices_global,
            self.cl_ghost_indices, self.cl_ghost_counts, self.cl_bond_indices_local,
            np.int32(self.num_atoms), np.int32(self.num_groups),
            np.float32(coll_margin**2), np.float32(bbox_margin)
        )

    def solve_cluster_jacobi(self, dt=0.01, iterations=10, k_coll=1000.0, omega=0.8, momentum_beta=0.0):
        """Run tiled Jacobi solver"""
        self.prg.solve_cluster_jacobi(
            self.queue, (self.num_groups * self.group_size,), (self.group_size,),
            self.cl_pos, self.cl_pred_pos, self.cl_params,
            self.cl_bond_indices_local, self.cl_bond_lengths, self.cl_bond_stiffness,
            self.cl_ghost_indices, self.cl_ghost_counts,
            np.int32(self.num_atoms), np.int32(iterations),
            np.float32(dt), np.float32(k_coll), np.float32(omega), np.float32(momentum_beta),
            self.cl_debug_force_bond, self.cl_debug_force_coll
        )

    def get_positions(self, out=None):
        """Download positions from GPU"""
        if out is None:
            out = self._host_pos
        cl.enqueue_copy(self.queue, out, self.cl_pos).wait()
        return out[:, :3]

    def get_debug_forces(self):
        """Download debug force buffers"""
        f_bond = np.empty((self.num_atoms, 4), dtype=np.float32)
        f_coll = np.empty((self.num_atoms, 4), dtype=np.float32)
        cl.enqueue_copy(self.queue, f_bond, self.cl_debug_force_bond)
        cl.enqueue_copy(self.queue, f_coll, self.cl_debug_force_coll).wait()
        return f_bond[:, :3], f_coll[:, :3]
