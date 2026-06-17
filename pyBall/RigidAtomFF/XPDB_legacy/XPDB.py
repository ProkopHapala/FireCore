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
        self.cl_iter_pos_prev = cl.Buffer(self.ctx, mf.READ_WRITE, num_atoms * 16)
        self.cl_iter_pos_prev2 = cl.Buffer(self.ctx, mf.READ_WRITE, num_atoms * 16)

        # Bounding Boxes
        self.cl_bboxes_min = cl.Buffer(self.ctx, mf.READ_WRITE, self.num_groups * 16)
        self.cl_bboxes_max = cl.Buffer(self.ctx, mf.READ_WRITE, self.num_groups * 16)

        # Collision List
        self.cl_coll_neighbors = cl.Buffer(self.ctx, mf.READ_WRITE, num_atoms * self.max_neighbors * 4) # int
        self.cl_coll_count = cl.Buffer(self.ctx, mf.READ_WRITE, num_atoms * 4) # int

        # Host placeholders for Bond topology (need to upload once)
        self.num_bonds = 0
        self.num_bond_entries = 0

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

        # Bond residual scratch (lazy init)
        self.cl_bond_abs_residual = None
        self._host_bond_abs_residual = None

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

        self.num_bond_entries = int(current_offset)

        self.cl_bond_start = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=np.array(start, dtype=np.int32))
        self.cl_bond_count = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=np.array(count, dtype=np.int32))
        self.cl_bond_neighbors = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=np.array(neighbors, dtype=np.int32))
        self.cl_bond_lengths = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=np.array(lengths, dtype=np.float32))
        self.cl_bond_stiffness = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=np.array(stiffness, dtype=np.float32))

        # Invalidate bond residual scratch (size depends on bond entries)
        self.cl_bond_abs_residual = None
        self._host_bond_abs_residual = None

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
        self.num_bond_entries = 0
        self.cl_bond_abs_residual = None
        self._host_bond_abs_residual = None

    def _ensure_bond_residual_scratch(self):
        if self.num_bond_entries <= 0:
            self.cl_bond_abs_residual = None
            self._host_bond_abs_residual = None
            return
        if (self.cl_bond_abs_residual is None) or (self._host_bond_abs_residual is None) or (len(self._host_bond_abs_residual) != self.num_bond_entries):
            mf = cl.mem_flags
            self.cl_bond_abs_residual = cl.Buffer(self.ctx, mf.READ_WRITE, self.num_bond_entries * 4)
            self._host_bond_abs_residual = np.empty(self.num_bond_entries, dtype=np.float32)

    def bond_residual_norms_from_pos(self, cl_pos_buffer, eps=1e-20):
        """Compute bond residual norms using OpenCL kernel bond_residuals_csr.

        Returns dict {linf, l2, n} where:
        - linf = max |dist - L0|
        - l2   = sqrt(mean(|dist - L0|^2))  (RMS)
        - n    = number of bond entries evaluated
        """
        if self.cl_bond_start is None:
            self.set_empty_bonds()
        if self.num_bond_entries <= 0:
            return {"linf": 0.0, "l2": 0.0, "n": 0}

        self._ensure_bond_residual_scratch()
        self.prg.bond_residuals_csr(
            self.queue, (self.num_atoms,), None,
            cl_pos_buffer,
            self.cl_bond_start, self.cl_bond_count, self.cl_bond_neighbors, self.cl_bond_lengths,
            self.cl_bond_abs_residual,
            np.int32(self.num_atoms)
        )
        cl.enqueue_copy(self.queue, self._host_bond_abs_residual, self.cl_bond_abs_residual).wait()

        r = self._host_bond_abs_residual
        linf = float(np.max(r))
        l2 = float(np.sqrt(np.mean(r * r) + float(eps)))
        return {"linf": linf, "l2": l2, "n": int(r.size)}

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

    def build_neighbor_list(self, margin=1.0):
        """Update collision neighbor list on GPU."""
        margin_sq = np.float32(margin * margin)
        self.prg.update_verlet_list(
            self.queue, (self.num_atoms,), None,
            self.cl_pos, self.cl_params,
            self.cl_bboxes_min, self.cl_bboxes_max,
            self.cl_coll_neighbors, self.cl_coll_count,
            np.int32(self.num_atoms), np.int32(self.num_groups),
            margin_sq, np.float32(2.0) # bbox margin
        )

    def step(
        self,
        dt=0.01,
        iterations=5,
        k_coll=1000.0,
        omega=0.8,
        momentum_beta=0.0,
        cheby_enable=False,
        cheby_rho=0.99,
        cheby_delay=5,
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
        cl.enqueue_copy(self.queue, self.cl_iter_pos_prev, self.cl_iter_pos_A)
        cl.enqueue_copy(self.queue, self.cl_iter_pos_prev2, self.cl_iter_pos_A)

        omega_k = 1.0

        # 2. Solver Loop
        for i in range(iterations):
            # Input: A, Output: B
            self.prg.solve_simple_jacobi(
                self.queue, (self.num_atoms,), None,
                # Inputs
                self.cl_iter_pos_A, self.cl_iter_pos_prev, self.cl_pred_pos, self.cl_params,
                # Output
                self.cl_iter_pos_B,
                self.cl_bond_start, self.cl_bond_count, self.cl_bond_neighbors, self.cl_bond_lengths, self.cl_bond_stiffness,
                self.cl_coll_neighbors, self.cl_coll_count,

                np.int32(self.num_atoms), np.float32(dt),
                np.float32(k_coll), np.float32(omega), np.float32(momentum_beta)
            )
            if cheby_enable:
                # Update Chebyshev omega
                if i < cheby_delay:
                    omega_k_next = 1.0
                elif i == cheby_delay:
                    omega_k_next = 2.0 / (2.0 - cheby_rho * cheby_rho)
                else:
                    omega_k_next = 4.0 / (4.0 - cheby_rho * cheby_rho * omega_k)
                self.prg.cheby_mix(
                    self.queue, (self.num_atoms,), None,
                    self.cl_iter_pos_prev2, self.cl_iter_pos_B, self.cl_iter_pos_B,
                    np.float32(omega_k_next)
                )
                omega_k = omega_k_next
            # Rotate buffers: prev2 <- prev, prev <- curr, curr <- B
            self.cl_iter_pos_prev2, self.cl_iter_pos_prev, self.cl_iter_pos_A, self.cl_iter_pos_B = (
                self.cl_iter_pos_prev,
                self.cl_iter_pos_A,
                self.cl_iter_pos_B,
                self.cl_iter_pos_prev2,
            )

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

    def diagnostics_from_pos(self, cl_pos_buffer):
        """Compute per-atom min/max bond violation and collision penetration from a specific CL buffer."""
        if self.cl_bond_start is None:
            self.set_empty_bonds()
        
        # Temporary buffers for diagnostics output
        out_bond = np.empty((self.num_atoms, 2), dtype=np.float32)
        out_coll = np.empty((self.num_atoms, 2), dtype=np.float32)
        cl_out_bond = cl.Buffer(self.ctx, cl.mem_flags.WRITE_ONLY, self.num_atoms * 8)
        cl_out_coll = cl.Buffer(self.ctx, cl.mem_flags.WRITE_ONLY, self.num_atoms * 8)

        self.prg.diag_constraints(
            self.queue, (self.num_atoms,), None,
            cl_pos_buffer, self.cl_params,
            self.cl_bond_start, self.cl_bond_count, self.cl_bond_neighbors, self.cl_bond_lengths,
            self.cl_coll_neighbors, self.cl_coll_count,
            cl_out_bond, cl_out_coll,
            np.int32(self.num_atoms)
        )

        cl.enqueue_copy(self.queue, out_bond, cl_out_bond)
        cl.enqueue_copy(self.queue, out_coll, cl_out_coll).wait()

        return {
            "bond_min": float(np.min(out_bond[:, 0])) if out_bond.size else 0.0,
            "bond_max": float(np.max(out_bond[:, 1])) if out_bond.size else 0.0,
            "coll_min": float(np.min(out_coll[:, 0])) if out_coll.size else 0.0,
            "coll_max": float(np.max(out_coll[:, 1])) if out_coll.size else 0.0,
        }

    def diagnostics(self):
        """Compute per-atom min/max bond violation and collision penetration, return summary."""
        return self.diagnostics_from_pos(self.cl_pos)

    def ensure_neighbors(self, margin=1.0):
        """Thin wrapper to rebuild neighbor list with a given margin."""
        return self.build_neighbor_list(margin=margin)

    def reset_iteration_buffers(self, src_buffer=None, reset_pred=True):
        """
        Seed iteration buffers (A, prev, prev2) from src_buffer (default cl_pos).
        Optionally seed cl_pred_pos as well.
        """
        src = self.cl_pos if src_buffer is None else src_buffer
        cl.enqueue_copy(self.queue, self.cl_iter_pos_A, src).wait()
        cl.enqueue_copy(self.queue, self.cl_iter_pos_prev, self.cl_iter_pos_A).wait()
        cl.enqueue_copy(self.queue, self.cl_iter_pos_prev2, self.cl_iter_pos_A).wait()
        if reset_pred:
            cl.enqueue_copy(self.queue, self.cl_pred_pos, src).wait()

    def jacobi_iteration(
        self,
        dt,
        k_coll,
        omega,
        momentum_beta,
        cheby_state=None,
        pred_buffer=None,
        collect_diag=False,
        out_pos_host=None,
    ):
        """
        Single Jacobi iteration over current iteration buffers.
        - pred_buffer: CL buffer to use as predicted positions (defaults to cl_pred_pos).
        - cheby_state: dict with keys {rho, delay, omega_k, iter}; updated in-place.
        - collect_diag: if True, returns diagnostics for the current buffer A after rotation.
        - out_pos_host: optional preallocated host array to receive positions (float4).
        Returns (diag or None, out_pos_host or None).
        """
        pred = self.cl_pred_pos if pred_buffer is None else pred_buffer
        self.prg.solve_simple_jacobi(
            self.queue, (self.num_atoms,), None,
            self.cl_iter_pos_A, self.cl_iter_pos_prev, pred, self.cl_params,
            self.cl_iter_pos_B,
            self.cl_bond_start, self.cl_bond_count, self.cl_bond_neighbors, self.cl_bond_lengths, self.cl_bond_stiffness,
            self.cl_coll_neighbors, self.cl_coll_count,
            np.int32(self.num_atoms), np.float32(dt),
            np.float32(k_coll), np.float32(omega), np.float32(momentum_beta)
        )
        if cheby_state is not None:
            i = cheby_state.get("iter", 0)
            omega_k = cheby_state.get("omega_k", 1.0)
            rho = cheby_state.get("rho", 0.99)
            delay = cheby_state.get("delay", 0)
            if i < delay:
                omega_k_next = 1.0
            elif i == delay:
                omega_k_next = 2.0 / (2.0 - rho * rho)
            else:
                omega_k_next = 4.0 / (4.0 - rho * rho * omega_k)
            self.prg.cheby_mix(
                self.queue, (self.num_atoms,), None,
                self.cl_iter_pos_prev2, self.cl_iter_pos_B, self.cl_iter_pos_B,
                np.float32(omega_k_next)
            )
            cheby_state["omega_k"] = omega_k_next
            cheby_state["iter"] = i + 1
        # Rotate buffers
        self.cl_iter_pos_prev2, self.cl_iter_pos_prev, self.cl_iter_pos_A, self.cl_iter_pos_B = (
            self.cl_iter_pos_prev,
            self.cl_iter_pos_A,
            self.cl_iter_pos_B,
            self.cl_iter_pos_prev2,
        )

        diag = None
        if collect_diag:
            diag = self.diagnostics_from_pos(self.cl_iter_pos_A)
        if out_pos_host is not None:
            cl.enqueue_copy(self.queue, out_pos_host, self.cl_iter_pos_A).wait()
        return diag, out_pos_host

    def run_jacobi_iterations(
        self,
        iterations,
        dt=0.01,
        k_coll=1000.0,
        omega=0.8,
        momentum_beta=0.0,
        cheby_enable=False,
        cheby_rho=0.99,
        cheby_delay=5,
        collect_diagnostics=True,
    ):
        """
        Run Jacobi iterations without integration. Assumes cl_pos/vel/params and neighbor lists are ready.
        Returns (coll_series, bond_series) if collect_diagnostics else None.
        """
        # Initialize buffers to current position
        cl.enqueue_copy(self.queue, self.cl_pred_pos, self.cl_pos).wait()
        cl.enqueue_copy(self.queue, self.cl_iter_pos_A, self.cl_pos).wait()
        cl.enqueue_copy(self.queue, self.cl_iter_pos_prev, self.cl_iter_pos_A).wait()
        cl.enqueue_copy(self.queue, self.cl_iter_pos_prev2, self.cl_iter_pos_A).wait()

        coll_series = []
        bond_series = []
        omega_k = 1.0

        for i in range(iterations):
            self.prg.solve_simple_jacobi(
                self.queue, (self.num_atoms,), None,
                self.cl_iter_pos_A, self.cl_iter_pos_prev, self.cl_pred_pos, self.cl_params,
                self.cl_iter_pos_B,
                self.cl_bond_start, self.cl_bond_count, self.cl_bond_neighbors, self.cl_bond_lengths, self.cl_bond_stiffness,
                self.cl_coll_neighbors, self.cl_coll_count,
                np.int32(self.num_atoms), np.float32(dt),
                np.float32(k_coll), np.float32(omega), np.float32(momentum_beta)
            )
            if cheby_enable:
                if i < cheby_delay:
                    omega_k_next = 1.0
                elif i == cheby_delay:
                    omega_k_next = 2.0 / (2.0 - cheby_rho * cheby_rho)
                else:
                    omega_k_next = 4.0 / (4.0 - cheby_rho * cheby_rho * omega_k)
                self.prg.cheby_mix(
                    self.queue, (self.num_atoms,), None,
                    self.cl_iter_pos_prev2, self.cl_iter_pos_B, self.cl_iter_pos_B,
                    np.float32(omega_k_next)
                )
                omega_k = omega_k_next
            # Rotate buffers
            self.cl_iter_pos_prev2, self.cl_iter_pos_prev, self.cl_iter_pos_A, self.cl_iter_pos_B = (
                self.cl_iter_pos_prev,
                self.cl_iter_pos_A,
                self.cl_iter_pos_B,
                self.cl_iter_pos_prev2,
            )
            if collect_diagnostics:
                diag = self.diagnostics_from_pos(self.cl_iter_pos_A)
                coll_violation = max(-diag["coll_min"], 0.0)
                bond_violation = max(abs(diag["bond_min"]), abs(diag["bond_max"]))
                coll_series.append(coll_violation)
                bond_series.append(bond_violation)

        if collect_diagnostics:
            return np.array(coll_series), np.array(bond_series)
        return None

    ## --------------------------------------------------------
    ## -------- Gauss-Seidel solver (with coloring) 
    ## --------------------------------------------------------

    def init_GS(self, num_atoms, group_size=64):
        self.num_atoms = num_atoms
        self.group_size = group_size
        self.ctx = cl.create_some_context()
        self.queue = cl.CommandQueue(self.ctx)
        
        # Load Kernel
        curr_dir = os.path.dirname(os.path.abspath(__file__))
        with open("GS_Solver.cl", "r") as f:
            src = f.read()
        self.prg = cl.Program(self.ctx, src).build()
        
        # Allocate basic buffers
        mf = cl.mem_flags
        self.cl_pos = cl.Buffer(self.ctx, mf.READ_WRITE, num_atoms * 16)
        self.cl_pred = cl.Buffer(self.ctx, mf.READ_WRITE, num_atoms * 16)
        self.cl_params = cl.Buffer(self.ctx, mf.READ_ONLY, num_atoms * 16)
        self.cl_colors = cl.Buffer(self.ctx, mf.READ_ONLY, num_atoms * 4)

    def upload_topology_GS(self, bonds, colors):
        # bonds is list of lists: [[(neighbor_idx, len, k), ...], ...]
        flat_neigh = []
        flat_len = []
        flat_stiff = []
        start = []
        count = []
        
        offset = 0
        for i, b_list in enumerate(bonds):
            start.append(offset)
            count.append(len(b_list))
            for b in b_list:
                flat_neigh.append(b[0])
                flat_len.append(b[1])
                flat_stiff.append(b[2])
            offset += len(b_list)

        mf = cl.mem_flags
        self.cl_bond_start = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=np.array(start, dtype=np.int32))
        self.cl_bond_count = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=np.array(count, dtype=np.int32))
        self.cl_bond_neigh = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=np.array(flat_neigh, dtype=np.int32))
        self.cl_bond_len = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=np.array(flat_len, dtype=np.float32))
        self.cl_bond_k = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=np.array(flat_stiff, dtype=np.float32))
        
        # Upload Colors
        cl.enqueue_copy(self.queue, self.cl_colors, np.array(colors, dtype=np.int32))

    def solve_GS(self, dt=0.016, inner_iters=10):
        # Assuming cl_pos and cl_pred are already set up (omitted prediction step for brevity)
        
        self.prg.solve_gauss_seidel_local(
            self.queue, (self.num_atoms,), (self.group_size,),
            self.cl_pos, self.cl_pred, self.cl_params,
            self.cl_bond_start, self.cl_bond_count, self.cl_bond_neigh, self.cl_bond_len, self.cl_bond_k,
            self.cl_colors, np.int32(2), # 2 Colors (Red/Black)
            np.int32(inner_iters), np.float32(dt), np.float32(1.0) # Omega
        )

    ## --------------------------------------------------------
    ## -------- Tiled Topology & Solver (BBox based)
    ## --------------------------------------------------------

    def init_tiled(self, max_ghosts=128):
        self.max_ghosts = max_ghosts
        mf = cl.mem_flags
        # Bond indices in int4 format for the tiled kernel
        self.cl_global_bond_indices = cl.Buffer(self.ctx, mf.READ_ONLY, self.num_atoms * 16)
        self.cl_local_bond_indices  = cl.Buffer(self.ctx, mf.READ_WRITE, self.num_atoms * 16)
        self.cl_global_bond_L4      = cl.Buffer(self.ctx, mf.READ_ONLY, self.num_atoms * 16)
        self.cl_global_bond_K4      = cl.Buffer(self.ctx, mf.READ_ONLY, self.num_atoms * 16)
        
        # Ghost management
        max_ghosts_total = self.num_groups * max_ghosts
        self.cl_ghost_indices = cl.Buffer(self.ctx, mf.READ_WRITE, max_ghosts_total * 4)
        self.cl_ghost_counts  = cl.Buffer(self.ctx, mf.READ_WRITE, self.num_groups * 4)

        # Debug buffers
        self.cl_debug_force_bond = cl.Buffer(self.ctx, mf.WRITE_ONLY, self.num_atoms * 16)
        self.cl_debug_force_coll = cl.Buffer(self.ctx, mf.WRITE_ONLY, self.num_atoms * 16)
        # For neighbor inspection: 64 slots per atom
        self.cl_debug_neigh_indices = cl.Buffer(self.ctx, mf.WRITE_ONLY, self.num_atoms * 64 * 4)
        self.cl_debug_neigh_pos     = cl.Buffer(self.ctx, mf.WRITE_ONLY, self.num_atoms * 64 * 16)

    def upload_global_bonds(self, bonds):
        """
        Convert list-of-lists bonds to int4 global indices (max 4 per atom)
        and also upload float4 lengths and stiffness.
        """
        b_indices = np.full((self.num_atoms, 4), -1, dtype=np.int32)
        b_L4 = np.zeros((self.num_atoms, 4), dtype=np.float32)
        b_K4 = np.zeros((self.num_atoms, 4), dtype=np.float32)
        for i, b_list in enumerate(bonds):
            for k, b in enumerate(b_list[:4]):
                b_indices[i, k] = b[0]
                b_L4[i, k] = b[1]
                b_K4[i, k] = b[2]
        cl.enqueue_copy(self.queue, self.cl_global_bond_indices, b_indices).wait()
        cl.enqueue_copy(self.queue, self.cl_global_bond_L4, b_L4).wait()
        cl.enqueue_copy(self.queue, self.cl_global_bond_K4, b_K4).wait()

    def build_tiled_topology(self, Rmax, coll_scale=2.0, bbox_scale=2.0):
        """
        Launch build_local_topology kernel using Rmax instead of explicit margins.
        - Collision search margin: coll_scale * Rmax (squared passed to kernels)
        - BBox padding: bbox_scale * Rmax
        """
        coll_margin = coll_scale * Rmax
        bbox_margin = bbox_scale * Rmax
        # Ensure bboxes are up to date
        self.update_verlet(margin_sq=coll_margin**2) # Using existing update_verlet which computes bboxes
        # Global size must be total threads: num_groups * group_size
        self.prg.build_local_topology(
            self.queue, (self.num_groups * self.group_size,), (self.group_size,),
            self.cl_pos, self.cl_bboxes_min, self.cl_bboxes_max, self.cl_global_bond_indices,
            self.cl_ghost_indices, self.cl_ghost_counts, self.cl_local_bond_indices,
            np.int32(self.num_atoms), np.int32(self.num_groups),
            np.float32(coll_margin**2), np.float32(bbox_margin)
        )

    def solve_tiled_jacobi(self, dt=0.01, iterations=10, k_coll=1000.0, omega=0.8, momentum_beta=0.0):
        """
        Tiled Jacobi solver using solve_cluster_jacobi kernel.
        """
        # Global size must be multiple of group_size
        self.prg.solve_cluster_jacobi(
            self.queue, (self.num_groups * self.group_size,), (self.group_size,),
            self.cl_pos, self.cl_pred_pos, self.cl_params,
            self.cl_local_bond_indices, self.cl_global_bond_L4, self.cl_global_bond_K4,
            self.cl_ghost_indices, self.cl_ghost_counts,
            np.int32(self.num_atoms), np.int32(iterations),
            np.float32(dt), np.float32(k_coll), np.float32(omega), np.float32(momentum_beta),
            self.cl_debug_force_bond, self.cl_debug_force_coll,
            self.cl_debug_neigh_indices, self.cl_debug_neigh_pos
        )

    ## --------------------------------------------------------
    ## -------- Gauss-Seidel Block (keep all in local memory)
    ## --------------------------------------------------------

    def init_GSblock(self, num_atoms, wg_size=64, num_colors=8):
        self.wg_size = wg_size
        self.num_colors = num_colors
        self.block_cap = wg_size * num_colors
        
        # Padding: Num atoms must be multiple of BLOCK_CAPACITY for simplicity
        # or we handle bounds carefully.
        self.padded_size = ((num_atoms + self.block_cap - 1) // self.block_cap) * self.block_cap
        
        self.ctx = cl.create_some_context()
        self.queue = cl.CommandQueue(self.ctx)
        
        with open("BlockGS.cl", "r") as f:
            src = f.read()
        self.prg = cl.Program(self.ctx, src).build()
        
        mf = cl.mem_flags
        self.cl_pos = cl.Buffer(self.ctx, mf.READ_WRITE, self.padded_size * 16)
        self.cl_pred = cl.Buffer(self.ctx, mf.READ_WRITE, self.padded_size * 16)
        self.cl_params = cl.Buffer(self.ctx, mf.READ_ONLY, self.padded_size * 16)

    def upload_data_GSblock(self, pos, mass):
        # Pad data
        pad_len = self.padded_size - len(pos)
        pos_padded = np.pad(pos, ((0, pad_len), (0,0)), 'constant')
        mass_padded = np.pad(mass, (0, pad_len), 'constant')
        
        pos4 = np.zeros((self.padded_size, 4), dtype=np.float32)
        pos4[:, :3] = pos_padded
        
        params = np.zeros((self.padded_size, 4), dtype=np.float32)
        params[:, 3] = mass_padded
        
        cl.enqueue_copy(self.queue, self.cl_pos, pos4)
        cl.enqueue_copy(self.queue, self.cl_pred, pos4) # pred = pos initially
        cl.enqueue_copy(self.queue, self.cl_params, params)

    def upload_topology_GSblock(self, bonds_adj):
        # Flatten logic similar to before
        start = []
        count = []
        neigh = []
        L = []
        K = []
        
        offset = 0
        # We must iterate up to PADDED size
        for i in range(self.padded_size):
            start.append(offset)
            if i < len(bonds_adj):
                b_list = bonds_adj[i]
                count.append(len(b_list))
                for b in b_list:
                    neigh.append(b[0])
                    L.append(b[1])
                    K.append(b[2])
                offset += len(b_list)
            else:
                count.append(0)

        self.cl_b_start = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=np.array(start, dtype=np.int32))
        self.cl_b_count = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=np.array(count, dtype=np.int32))
        self.cl_b_neigh = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=np.array(neigh, dtype=np.int32))
        self.cl_b_len = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=np.array(L, dtype=np.float32))
        self.cl_b_k = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=np.array(K, dtype=np.float32))

    def solve_GS_block(self, dt=0.016, inner_iters=20):
        # Global work size = Padded Size / Num Colors
        # Because we launch 1 thread per (K atoms)
        # Wait, actually:
        # In Kernel: Global ID corresponds to lid + grp*BLOCK_CAP
        # The kernel logic uses `get_group_id` * BLOCK_CAPACITY.
        # We need to launch (PaddedSize / BLOCK_CAPACITY) * WG_SIZE threads.
        # i.e., Total atoms / Num Colors.
        
        gws = int(self.padded_size / self.num_colors)
        lws = int(self.wg_size)
        
        self.prg.solve_block_gs(
            self.queue, (gws,), (lws,),
            self.cl_pos, self.cl_pred, self.cl_params,
            self.cl_b_start, self.cl_b_count, self.cl_b_neigh, self.cl_b_len, self.cl_b_k,
            np.int32(self.padded_size), np.int32(inner_iters), np.float32(dt), np.float32(1.0)
        )