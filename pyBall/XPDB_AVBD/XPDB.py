import os

import numpy as np
import pyopencl as cl


class XPDB:
    def __init__(self, num_atoms, group_size=64, prefer_gpu=True, device_idx=0):
        self.num_atoms = num_atoms
        self.group_size = group_size
        self.num_groups = (num_atoms + group_size - 1) // group_size
        self.max_neighbors = 64 # Must match CL define

        # OpenCL Setup (prefer explicit GPU selection to avoid interactive prompts)
        self.ctx = self._make_context(prefer_gpu=prefer_gpu, device_idx=device_idx)
        self.queue = cl.CommandQueue(self.ctx)

        # Load Source
        curr_dir = os.path.dirname(os.path.abspath(__file__))
        with open(os.path.join(curr_dir, "XPDB.cl"), "r") as f:
            src = f.read()
        self.prg = cl.Program(self.ctx, src).build()

        # Buffers
        mf = cl.mem_flags

        # Atom Data
        self.cl_pos = cl.Buffer(self.ctx, mf.READ_WRITE, num_atoms * 16) # float4
        self.cl_vel = cl.Buffer(self.ctx, mf.READ_WRITE, num_atoms * 16)
        self.cl_params = cl.Buffer(self.ctx, mf.READ_ONLY, num_atoms * 16) # .x=rad, .w=mass

        # Solver Temp Buffers (Double Buffer)
        self.cl_pred_pos = cl.Buffer(self.ctx, mf.READ_WRITE, num_atoms * 16)
        self.cl_iter_pos_A = cl.Buffer(self.ctx, mf.READ_WRITE, num_atoms * 16)
        self.cl_iter_pos_B = cl.Buffer(self.ctx, mf.READ_WRITE, num_atoms * 16)

        # Bounding Boxes
        self.cl_bboxes_min = cl.Buffer(self.ctx, mf.READ_WRITE, self.num_groups * 16)
        self.cl_bboxes_max = cl.Buffer(self.ctx, mf.READ_WRITE, self.num_groups * 16)

        # Collision List
        self.cl_coll_neighbors = cl.Buffer(self.ctx, mf.READ_WRITE, num_atoms * self.max_neighbors * 4) # int
        self.cl_coll_count = cl.Buffer(self.ctx, mf.READ_WRITE, num_atoms * 4) # int

        # Host placeholders for Bond topology (need to upload once)
        self.num_bonds = 0

        # Host scratch for readback to avoid per-frame allocation
        self._host_pos = np.zeros((self.num_atoms, 4), dtype=np.float32)
        self._host_diag = {
            "bond": np.zeros((self.num_atoms, 2), dtype=np.float32),
            "coll": np.zeros((self.num_atoms, 2), dtype=np.float32),
        }

        # Diagnostics buffers (lazy init)
        self.cl_diag_bond = None
        self.cl_diag_coll = None
        # Bond buffers (may stay empty if no bonds)
        self.cl_bond_start = None
        self.cl_bond_count = None
        self.cl_bond_neighbors = None
        self.cl_bond_lengths = None
        self.cl_bond_stiffness = None

    def _ensure_diag_buffers(self):
        if self.cl_diag_bond is None:
            mf = cl.mem_flags
            self.cl_diag_bond = cl.Buffer(self.ctx, mf.READ_WRITE, self.num_atoms * 8)  # float2
            self.cl_diag_coll = cl.Buffer(self.ctx, mf.READ_WRITE, self.num_atoms * 8)  # float2

    def _make_context(self, prefer_gpu=True, device_idx=0):
        dev_type = cl.device_type.GPU if prefer_gpu else cl.device_type.ALL
        try:
            for platform in cl.get_platforms():
                gpus = platform.get_devices(device_type=dev_type)
                if gpus:
                    return cl.Context(devices=[gpus[device_idx % len(gpus)]])
        except Exception:
            pass
        # Fallback to PyOpenCL default chooser (may prompt if DISPLAY present)
        return cl.create_some_context(interactive=False)

    def upload_data(self, pos, vel, radius, mass):
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
        cl.enqueue_copy(self.queue, self.cl_vel, vel4)
        cl.enqueue_copy(self.queue, self.cl_params, params)

    def upload_bonds(self, bonds_adj):
        """
        bonds_adj: list of lists or similar structure.
        We need flattened CSR-like format: start, count, neighbors, length, stiffness
        """
        start = []
        count = []
        neighbors = []
        lengths = []
        stiffness = []

        current_offset = 0
        for i in range(self.num_atoms):
            start.append(current_offset)
            b_list = bonds_adj[i] # Expecting [(j, length, k), ...]
            count.append(len(b_list))
            for b in b_list:
                neighbors.append(b[0])
                lengths.append(b[1])
                stiffness.append(b[2])
            current_offset += len(b_list)

        self.cl_bond_start = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=np.array(start, dtype=np.int32))
        self.cl_bond_count = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=np.array(count, dtype=np.int32))
        self.cl_bond_neighbors = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=np.array(neighbors, dtype=np.int32))
        self.cl_bond_lengths = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=np.array(lengths, dtype=np.float32))
        self.cl_bond_stiffness = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=np.array(stiffness, dtype=np.float32))

    def set_empty_bonds(self):
        """Prepare zero-length bond buffers to allow stepping without bonds."""
        zero_i = np.zeros(self.num_atoms, dtype=np.int32)
        zero_f = np.zeros(1, dtype=np.float32)
        mf = cl.mem_flags
        self.cl_bond_start = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=zero_i)
        self.cl_bond_count = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=zero_i)
        self.cl_bond_neighbors = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=np.zeros(1, dtype=np.int32))
        self.cl_bond_lengths = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=zero_f)
        self.cl_bond_stiffness = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=zero_f)

    def update_verlet(self, margin_sq=4.0):
        # 1. Update BBoxes
        # Local memory size for reduction
        loc_mem = cl.LocalMemory(self.group_size * 16) # float4

        self.prg.update_bboxes(
            self.queue, (self.num_groups * self.group_size,), (self.group_size,),
            self.cl_pos, self.cl_params,
            self.cl_bboxes_min, self.cl_bboxes_max,
            loc_mem, loc_mem, # Two local buffers passed
            np.int32(self.num_atoms)
        )

        # 2. Update List
        self.prg.update_verlet_list(
            self.queue, (self.num_atoms,), None,
            self.cl_pos, self.cl_params,
            self.cl_bboxes_min, self.cl_bboxes_max,
            self.cl_coll_neighbors, self.cl_coll_count,
            np.int32(self.num_atoms), np.int32(self.num_groups),
            np.float32(margin_sq), np.float32(2.0) # bbox margin
        )

    def step(
        self,
        dt=0.01,
        iterations=5,
        k_coll=1000.0,
        omega=0.8,
        pd_scale=1.0,
        gravity=(0.0, -9.81, 0.0),
        box_min=None,
        box_max=None,
        box_k=2000.0,
        bottom_y=None,
        bottom_c=0.0,
        bottom_x0=0.0,
        bottom_a=0.0,
        clamp_z=False,
    ):
        g = np.float32(gravity)
        bmin = np.float32(box_min if box_min is not None else (-1e9, -1e9, -1e9))
        bmax = np.float32(box_max if box_max is not None else (1e9, 1e9, 1e9))
        by = np.float32(bottom_y if bottom_y is not None else bmin[1])

        # Ensure bond buffers exist (may be empty)
        if self.cl_bond_start is None:
            self.set_empty_bonds()

        # 1. Predict
        self.prg.predict(
            self.queue, (self.num_atoms,), None,
            self.cl_pos, self.cl_vel, self.cl_params, self.cl_pred_pos,
            np.int32(self.num_atoms), np.float32(dt),
            g[0], g[1], g[2],
            bmin[0], bmin[1], bmin[2],
            bmax[0], bmax[1], bmax[2],
            np.float32(box_k), by, np.float32(bottom_c),
            np.float32(bottom_x0), np.float32(bottom_a),
        )

        # Init Iteration: Copy curr_pos (or pred_pos) to A
        cl.enqueue_copy(self.queue, self.cl_iter_pos_A, self.cl_pos) # Start from current

        # 2. Solver Loop
        for i in range(iterations):
            # Input: A, Output: B
            self.prg.solve_simple_jacobi(
                self.queue, (self.num_atoms,), None,
                # Inputs
                self.cl_iter_pos_A, self.cl_pred_pos, self.cl_params,
                # Output
                self.cl_iter_pos_B,
                self.cl_bond_start, self.cl_bond_count, self.cl_bond_neighbors, self.cl_bond_lengths, self.cl_bond_stiffness,
                self.cl_coll_neighbors, self.cl_coll_count,

                np.int32(self.num_atoms), np.float32(dt),
                np.float32(k_coll), np.float32(omega), np.float32(pd_scale)
            )
            # Swap A and B
            self.cl_iter_pos_A, self.cl_iter_pos_B = self.cl_iter_pos_B, self.cl_iter_pos_A

        # 3. Integrate (using A as final result because of swap)
        self.prg.integrate(
            self.queue, (self.num_atoms,), None,
            self.cl_pos, self.cl_iter_pos_A, self.cl_vel, self.cl_pos,
            np.int32(self.num_atoms), np.float32(dt)
        )
        if clamp_z:
            self.prg.clamp_z0(self.queue, (self.num_atoms,), None, self.cl_pos, self.cl_vel, np.int32(self.num_atoms))

    def get_positions(self, out=None):
        if out is None:
            out = self._host_pos
        cl.enqueue_copy(self.queue, out, self.cl_pos).wait()
        return out[:, :3]

    def diagnostics(self):
        """Compute per-atom min/max bond violation and collision penetration, return summary."""
        if self.cl_bond_start is None:
            self.set_empty_bonds()
        self._ensure_diag_buffers()
        self.prg.diag_constraints(
            self.queue, (self.num_atoms,), None,
            self.cl_pos, self.cl_params,
            self.cl_bond_start, self.cl_bond_count, self.cl_bond_neighbors, self.cl_bond_lengths,
            self.cl_coll_neighbors, self.cl_coll_count,
            self.cl_diag_bond, self.cl_diag_coll,
            np.int32(self.num_atoms)
        )
        cl.enqueue_copy(self.queue, self._host_diag["bond"], self.cl_diag_bond)
        cl.enqueue_copy(self.queue, self._host_diag["coll"], self.cl_diag_coll).wait()
        bond = self._host_diag["bond"]
        coll = self._host_diag["coll"]
        return {
            "bond_min": float(np.min(bond[:, 0])) if bond.size else 0.0,
            "bond_max": float(np.max(bond[:, 1])) if bond.size else 0.0,
            "coll_min": float(np.min(coll[:, 0])) if coll.size else 0.0,
            "coll_max": float(np.max(coll[:, 1])) if coll.size else 0.0,
        }