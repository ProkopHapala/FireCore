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
        self.cl_pred_pos = cl.Buffer(self.ctx, mf.READ_WRITE, num_atoms * 16)  # compatibility buffer
        self.cl_vel = cl.Buffer(self.ctx, mf.READ_WRITE, num_atoms * 16)  # float4 (kept for parity/debug)
        self.cl_params = cl.Buffer(self.ctx, mf.READ_ONLY, num_atoms * 16)  # .x=rad, .w=mass

        # Persistent prev positions + momentum (for heavy-ball mixing)
        self.cl_prev_pos = cl.Buffer(self.ctx, mf.READ_WRITE, num_atoms * 16)  # float4
        self.cl_momentum = cl.Buffer(self.ctx, mf.READ_WRITE, num_atoms * 16)  # float4

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

        # Bond residual scratch (lazy init)
        self.cl_bond_abs_residual = None
        self._host_bond_abs_residual = None

        # Debug Buffers
        self.cl_debug_force_bond = cl.Buffer(self.ctx, mf.WRITE_ONLY, num_atoms * 16)
        self.cl_debug_force_coll = cl.Buffer(self.ctx, mf.WRITE_ONLY, num_atoms * 16)

        # Host placeholders
        self._host_pos = np.zeros((self.num_atoms, 4), dtype=np.float32)

        # Rigid-atom (quat) solver buffers are allocated lazily
        self._rigid_inited = False
        self._rigid_host_pos = None
        self._rigid_host_quat = None
        self._rigid_host_vel = None
        self._rigid_host_omega = None

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

    def _init_rigid_buffers(self):
        if self._rigid_inited:
            return
        mf = cl.mem_flags
        n = int(self.num_atoms)

        self.cl_rpos_A = cl.Buffer(self.ctx, mf.READ_WRITE, n * 16)
        self.cl_rpos_B = cl.Buffer(self.ctx, mf.READ_WRITE, n * 16)
        self.cl_rvel_A = cl.Buffer(self.ctx, mf.READ_WRITE, n * 16)
        self.cl_rvel_B = cl.Buffer(self.ctx, mf.READ_WRITE, n * 16)
        self.cl_rquat_A = cl.Buffer(self.ctx, mf.READ_WRITE, n * 16)
        self.cl_rquat_B = cl.Buffer(self.ctx, mf.READ_WRITE, n * 16)
        self.cl_romega_A = cl.Buffer(self.ctx, mf.READ_WRITE, n * 16)
        self.cl_romega_B = cl.Buffer(self.ctx, mf.READ_WRITE, n * 16)

        self.cl_rneighs = cl.Buffer(self.ctx, mf.READ_ONLY, n * 16)     # int4
        self.cl_rbkneighs = cl.Buffer(self.ctx, mf.READ_ONLY, n * 16)   # int4
        self.cl_rbLs = cl.Buffer(self.ctx, mf.READ_ONLY, n * 16)        # float4
        self.cl_rbKs = cl.Buffer(self.ctx, mf.READ_ONLY, n * 16)        # float4

        # Port-based rigid solver buffers
        self.cl_ratom_types = cl.Buffer(self.ctx, mf.READ_ONLY, n * 1)   # uchar per atom
        self.cl_rglobal_ports = cl.Buffer(self.ctx, mf.READ_WRITE, n * 4 * 16)  # float4[ natoms*4 ]
        self.cl_rdelta_p = cl.Buffer(self.ctx, mf.READ_WRITE, n * 16)    # float4
        self.cl_rdelta_rot = cl.Buffer(self.ctx, mf.READ_WRITE, n * 16)  # float4
        self.cl_rpos_tmp = cl.Buffer(self.ctx, mf.READ_WRITE, n * 16)    # float4
        self.cl_rquat_tmp = cl.Buffer(self.ctx, mf.READ_WRITE, n * 16)   # float4

        # Explicit force-based rigid dynamics buffers
        self.cl_rforce  = cl.Buffer(self.ctx, mf.READ_WRITE, n * 16)     # float4
        self.cl_rfneigh = cl.Buffer(self.ctx, mf.READ_WRITE, n * 4 * 16) # float4[ natoms*4 ] recoil forces
        self.cl_rpneigh = cl.Buffer(self.ctx, mf.READ_WRITE, n * 4 * 16) # float4[ natoms*4 ] port lever arms r_arm

        # Explicit force-based port geometry (per atom, local frame, length-scaled)
        self.cl_rport_local = cl.Buffer(self.ctx, mf.READ_ONLY, n * 4 * 16) # float4[ natoms*4 ]
        self.cl_rport_n     = cl.Buffer(self.ctx, mf.READ_ONLY, n * 1)      # uchar[ natoms ]

        # Port template tables (type->dirs, type->nports)
        # N_PORT_TYPES = 4, dirs per type = 4 float4 (xyz dir, w unused)
        self.cl_port_dirs = cl.Buffer(self.ctx, mf.READ_ONLY, 4 * 4 * 16)
        self.cl_port_ns   = cl.Buffer(self.ctx, mf.READ_ONLY, 4 * 1)

        self._rigid_host_pos = np.empty((n, 4), dtype=np.float32)
        self._rigid_host_quat = np.empty((n, 4), dtype=np.float32)
        self._rigid_host_vel = np.empty((n, 4), dtype=np.float32)
        self._rigid_host_omega = np.empty((n, 4), dtype=np.float32)
        self._rigid_host_types = np.zeros((n,), dtype=np.uint8)

        self._upload_default_port_tables()

        self._rigid_inited = True

    def _upload_default_port_tables(self):
        # type 0: none
        # type 1: sp1 (2 ports)
        # type 2: sp2 (3 ports)
        # type 3: sp3 (4 ports)
        port_ns = np.array([0, 2, 3, 4], dtype=np.uint8)
        port_dirs = np.zeros((4, 4, 4), dtype=np.float32)
        # sp1
        port_dirs[1, 0, :3] = (1.0, 0.0, 0.0)
        port_dirs[1, 1, :3] = (-1.0, 0.0, 0.0)
        # sp2
        port_dirs[2, 0, :3] = (1.0, 0.0, 0.0)
        port_dirs[2, 1, :3] = (-0.5, 0.8660254, 0.0)
        port_dirs[2, 2, :3] = (-0.5, -0.8660254, 0.0)
        port_dirs[2, 3, :3] = (0.0, 0.0, 1.0)  # optional pi
        # sp3
        s = 0.57735026919
        port_dirs[3, 0, :3] = (s, s, s)
        port_dirs[3, 1, :3] = (s, -s, -s)
        port_dirs[3, 2, :3] = (-s, s, -s)
        port_dirs[3, 3, :3] = (-s, -s, s)
        cl.enqueue_copy(self.queue, self.cl_port_dirs, port_dirs.reshape(-1, 4))
        cl.enqueue_copy(self.queue, self.cl_port_ns, port_ns).wait()

    def upload_port_tables(self, port_dirs, port_ns):
        self._init_rigid_buffers()
        pd = np.array(port_dirs, dtype=np.float32, copy=False)
        pn = np.array(port_ns, dtype=np.uint8, copy=False)
        if pd.shape != (4, 4, 4):
            raise ValueError(f'upload_port_tables: port_dirs.shape={pd.shape} expected (4,4,4)')
        if pn.shape != (4,):
            raise ValueError(f'upload_port_tables: port_ns.shape={pn.shape} expected (4,)')
        cl.enqueue_copy(self.queue, self.cl_port_dirs, pd.reshape(-1, 4))
        cl.enqueue_copy(self.queue, self.cl_port_ns, pn).wait()

    def upload_rigid_state(self, pos, mass, *, quat=None, vel=None, omega=None):
        self._init_rigid_buffers()
        n = int(self.num_atoms)
        if pos.shape[0] != n:
            raise ValueError(f'upload_rigid_state: pos.shape={pos.shape} expected ({n},3)')
        if mass.shape[0] != n:
            raise ValueError(f'upload_rigid_state: mass.shape={mass.shape} expected ({n},)')

        pos4 = np.zeros((n, 4), dtype=np.float32)
        pos4[:, :3] = pos[:, :3]
        pos4[:, 3] = 1.0 / mass

        if vel is None:
            vel4 = np.zeros((n, 4), dtype=np.float32)
        else:
            vel4 = np.zeros((n, 4), dtype=np.float32)
            vel4[:, :3] = vel[:, :3]

        if omega is None:
            om4 = np.zeros((n, 4), dtype=np.float32)
        else:
            om4 = np.zeros((n, 4), dtype=np.float32)
            om4[:, :3] = omega[:, :3]

        if quat is None:
            q4 = np.zeros((n, 4), dtype=np.float32)
            q4[:, 3] = 1.0
        else:
            q4 = np.array(quat, dtype=np.float32, copy=False)
            if q4.shape != (n, 4):
                raise ValueError(f'upload_rigid_state: quat.shape={q4.shape} expected ({n},4)')

        cl.enqueue_copy(self.queue, self.cl_rpos_A, pos4)
        cl.enqueue_copy(self.queue, self.cl_rpos_B, pos4)
        cl.enqueue_copy(self.queue, self.cl_rvel_A, vel4)
        cl.enqueue_copy(self.queue, self.cl_rvel_B, vel4)
        cl.enqueue_copy(self.queue, self.cl_rquat_A, q4)
        cl.enqueue_copy(self.queue, self.cl_rquat_B, q4)
        cl.enqueue_copy(self.queue, self.cl_romega_A, om4)
        cl.enqueue_copy(self.queue, self.cl_romega_B, om4).wait()

    def upload_rigid_topology(self, neighs, bkNeighs, bLs, bKs):
        self._init_rigid_buffers()
        n = int(self.num_atoms)
        neighs = np.array(neighs, dtype=np.int32, copy=False)
        bkNeighs = np.array(bkNeighs, dtype=np.int32, copy=False)
        bLs = np.array(bLs, dtype=np.float32, copy=False)
        bKs = np.array(bKs, dtype=np.float32, copy=False)
        if neighs.shape != (n, 4):
            raise ValueError(f'upload_rigid_topology: neighs.shape={neighs.shape} expected ({n},4)')
        if bkNeighs.shape != (n, 4):
            raise ValueError(f'upload_rigid_topology: bkNeighs.shape={bkNeighs.shape} expected ({n},4)')
        if bLs.shape != (n, 4):
            raise ValueError(f'upload_rigid_topology: bLs.shape={bLs.shape} expected ({n},4)')
        if bKs.shape != (n, 4):
            raise ValueError(f'upload_rigid_topology: bKs.shape={bKs.shape} expected ({n},4)')
        cl.enqueue_copy(self.queue, self.cl_rneighs, neighs)
        cl.enqueue_copy(self.queue, self.cl_rbkneighs, bkNeighs)
        cl.enqueue_copy(self.queue, self.cl_rbLs, bLs)
        cl.enqueue_copy(self.queue, self.cl_rbKs, bKs).wait()

    def upload_rigid_atom_types(self, atom_types):
        self._init_rigid_buffers()
        n = int(self.num_atoms)
        at = np.array(atom_types, dtype=np.uint8, copy=False)
        if at.shape != (n,):
            raise ValueError(f'upload_rigid_atom_types: atom_types.shape={at.shape} expected ({n},)')
        cl.enqueue_copy(self.queue, self.cl_ratom_types, at).wait()

    def upload_rigid_ports_local(self, port_local, port_n):
        self._init_rigid_buffers()
        n = int(self.num_atoms)
        pl = np.array(port_local, dtype=np.float32, copy=False)
        pn = np.array(port_n, dtype=np.uint8, copy=False)
        if pl.shape != (n, 4, 4):
            raise ValueError(f'upload_rigid_ports_local: port_local.shape={pl.shape} expected ({n},4,4)')
        if pn.shape != (n,):
            raise ValueError(f'upload_rigid_ports_local: port_n.shape={pn.shape} expected ({n},)')
        cl.enqueue_copy(self.queue, self.cl_rport_local, pl.reshape(n*4, 4))
        cl.enqueue_copy(self.queue, self.cl_rport_n, pn).wait()

    def download_buffer(self, clbuf, shape, *, dtype=np.float32):
        out = np.empty(shape, dtype=dtype)
        cl.enqueue_copy(self.queue, out, clbuf).wait()
        return out

    def rigid_projective_step(self, *, nnode, dt=0.1, iterations=10):
        """Projective dynamics (position-based) step: project_ports + jacobi_solve_rigid ping-pong."""
        self._init_rigid_buffers()
        nnode = int(nnode)
        if nnode < 0 or nnode > self.num_atoms:
            raise ValueError(f'rigid_projective_step: nnode={nnode} out of range [0,{self.num_atoms}]')
        natoms = np.int32(self.num_atoms)
        nnode_i = np.int32(nnode)
        dt_f = np.float32(float(dt))

        pos_in = self.cl_rpos_A
        quat_in = self.cl_rquat_A
        pos_out = self.cl_rpos_tmp
        quat_out = self.cl_rquat_tmp

        for _ in range(int(iterations)):
            if nnode > 0:
                self.prg.project_ports(self.queue, (nnode,), None, nnode_i, pos_in, quat_in, self.cl_rglobal_ports, self.cl_rbLs, self.cl_ratom_types, self.cl_port_dirs, self.cl_port_ns)
            self.prg.jacobi_solve_rigid(
                self.queue, (self.num_atoms,), None,
                natoms, nnode_i,
                pos_in, quat_in,
                self.cl_rglobal_ports,
                pos_out, quat_out,
                self.cl_rneighs, self.cl_rbkneighs,
                self.cl_rbKs, self.cl_rbLs, self.cl_ratom_types, self.cl_port_dirs, self.cl_port_ns,
                dt_f
            )
            pos_in, pos_out = pos_out, pos_in
            quat_in, quat_out = quat_out, quat_in

        self.cl_rpos_A = pos_in
        self.cl_rquat_A = quat_in

    def rigid_force_step(self, *, nnode, dt=0.1, iterations=10, relax_p=1.0, relax_q=1.0):
        """Force/impulse backup step: integrate_and_project + solve_ports_xpbd + apply_deltas_rigid."""
        self._init_rigid_buffers()
        nnode = int(nnode)
        if nnode < 0 or nnode > self.num_atoms:
            raise ValueError(f'rigid_force_step: nnode={nnode} out of range [0,{self.num_atoms}]')
        natoms = np.int32(self.num_atoms)
        nnode_i = np.int32(nnode)
        dt_f = np.float32(float(dt))
        rp = np.float32(float(relax_p))
        rq = np.float32(float(relax_q))

        for _ in range(int(iterations)):
            # relaxation-style force step (no inertial integration inside constraint iterations)
            if nnode > 0:
                self.prg.project_ports(self.queue, (nnode,), None, nnode_i, self.cl_rpos_A, self.cl_rquat_A, self.cl_rglobal_ports, self.cl_rbLs, self.cl_ratom_types, self.cl_port_dirs, self.cl_port_ns)
            self.prg.solve_ports_xpbd(
                self.queue, (self.num_atoms,), None,
                natoms, nnode_i,
                self.cl_rpos_A, self.cl_rquat_A, self.cl_rglobal_ports,
                self.cl_rneighs, self.cl_rbkneighs,
                self.cl_rbLs, self.cl_rbKs,
                self.cl_rdelta_p, self.cl_rdelta_rot,
                self.cl_ratom_types,
                self.cl_port_dirs, self.cl_port_ns,
                dt_f
            )
            self.prg.apply_deltas_rigid(
                self.queue, (self.num_atoms,), None,
                natoms, nnode_i,
                self.cl_rpos_A, self.cl_rquat_A,
                self.cl_rdelta_p, self.cl_rdelta_rot,
                rp, rq
            )

    def rigid_force_explicit_step(self, *, nnode, dt=0.01, nsteps=1, damp=1.0):
        """Explicit force-based rigid dynamics using port springs + recoil gather, with omega+quat integration."""
        self._init_rigid_buffers()
        nnode = int(nnode)
        if nnode < 0 or nnode > self.num_atoms:
            raise ValueError(f'rigid_force_explicit_step: nnode={nnode} out of range [0,{self.num_atoms}]')
        natoms = np.int32(self.num_atoms)
        nnode_i = np.int32(nnode)
        dt_f = np.float32(float(dt))
        damp_f = np.float32(float(damp))
        for _ in range(int(nsteps)):
            self.prg.clear_rigid_forces(self.queue, (self.num_atoms,), None, natoms, self.cl_rforce, self.cl_rfneigh, self.cl_rpneigh)
            self.prg.gather_port_forces(
                self.queue, (self.num_atoms,), None,
                natoms, nnode_i,
                self.cl_rpos_A, self.cl_rquat_A,
                self.cl_rneighs, self.cl_rbkneighs,
                self.cl_rbKs,
                self.cl_rport_local,
                self.cl_rport_n,
                self.cl_rforce, self.cl_rfneigh, self.cl_rpneigh
            )
            self.prg.integrate_rigid_explicit(
                self.queue, (self.num_atoms,), None,
                natoms, nnode_i,
                self.cl_rpos_A, self.cl_rvel_A, self.cl_rquat_A, self.cl_romega_A,
                self.cl_rneighs, self.cl_rbkneighs,
                self.cl_rforce, self.cl_rfneigh, self.cl_rpneigh,
                dt_f, damp_f
            )

    def rigid_jacobi_bk(self, *, nnode, dt=0.1, iterations=10, k_rot=50.0):
        self._init_rigid_buffers()
        nnode = int(nnode)
        if nnode < 0 or nnode > self.num_atoms:
            raise ValueError(f'rigid_jacobi_bk: nnode={nnode} out of range [0,{self.num_atoms}]')
        dt = float(dt)
        k_rot = float(k_rot)
        natoms = np.int32(self.num_atoms)
        nnode_i = np.int32(nnode)
        dt_f = np.float32(dt)
        k_rot_f = np.float32(k_rot)

        pos_in = self.cl_rpos_A
        pos_out = self.cl_rpos_B
        vel_in = self.cl_rvel_A
        vel_out = self.cl_rvel_B
        quat_in = self.cl_rquat_A
        quat_out = self.cl_rquat_B
        omega_in = self.cl_romega_A
        omega_out = self.cl_romega_B

        for _ in range(int(iterations)):
            self.prg.solve_rigid_bk_symmetry(
                self.queue, (self.num_atoms,), None,
                natoms, nnode_i,
                pos_in, quat_in, vel_in, omega_in,
                pos_out, quat_out, vel_out, omega_out,
                self.cl_rneighs, self.cl_rbkneighs, self.cl_rbLs, self.cl_rbKs,
                dt_f, k_rot_f
            )
            pos_in, pos_out = pos_out, pos_in
            vel_in, vel_out = vel_out, vel_in
            quat_in, quat_out = quat_out, quat_in
            omega_in, omega_out = omega_out, omega_in

        self.cl_rpos_A = pos_in
        self.cl_rpos_B = pos_out
        self.cl_rvel_A = vel_in
        self.cl_rvel_B = vel_out
        self.cl_rquat_A = quat_in
        self.cl_rquat_B = quat_out
        self.cl_romega_A = omega_in
        self.cl_romega_B = omega_out

    def download_rigid_state(self):
        self._init_rigid_buffers()
        cl.enqueue_copy(self.queue, self._rigid_host_pos, self.cl_rpos_A)
        cl.enqueue_copy(self.queue, self._rigid_host_vel, self.cl_rvel_A)
        cl.enqueue_copy(self.queue, self._rigid_host_quat, self.cl_rquat_A)
        cl.enqueue_copy(self.queue, self._rigid_host_omega, self.cl_romega_A).wait()
        return self._rigid_host_pos, self._rigid_host_quat, self._rigid_host_vel, self._rigid_host_omega

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
        cl.enqueue_copy(self.queue, self.cl_prev_pos, pos4)
        cl.enqueue_copy(self.queue, self.cl_vel, vel4)
        cl.enqueue_copy(self.queue, self.cl_params, params)
        cl.enqueue_fill_buffer(self.queue, self.cl_momentum, np.float32(0.0), 0, self.num_atoms * 16)

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

        # Invalidate bond residual scratch
        self.cl_bond_abs_residual = None
        self._host_bond_abs_residual = None

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

        # Invalidate bond residual scratch
        self.cl_bond_abs_residual = None
        self._host_bond_abs_residual = None

    def reset_prev_pos(self, src_buffer=None):
        """Seed cl_prev_pos from src_buffer (default cl_pos) and clear momentum."""
        src = self.cl_pos if src_buffer is None else src_buffer
        cl.enqueue_copy(self.queue, self.cl_prev_pos, src).wait()
        cl.enqueue_fill_buffer(self.queue, self.cl_momentum, np.float32(0.0), 0, self.num_atoms * 16).wait()

    def _seed_prev_from_curr(self):
        cl.enqueue_copy(self.queue, self.cl_prev_pos, self.cl_pos)

    def _ensure_bond_residual_scratch(self):
        n = int(self.num_atoms * self.n_max_bonded)
        if (self.cl_bond_abs_residual is None) or (self._host_bond_abs_residual is None) or (self._host_bond_abs_residual.size != n):
            mf = cl.mem_flags
            self.cl_bond_abs_residual = cl.Buffer(self.ctx, mf.READ_WRITE, n * 4)
            self._host_bond_abs_residual = np.empty(n, dtype=np.float32)

    def bond_residual_norms(self, eps=1e-20):
        """Compute bond residual norms using fixed-slot global topology.

        Returns dict {linf, l2, n} where:
        - linf = max |dist - L0| over all slots
        - l2   = sqrt(mean(|dist - L0|^2)) over all slots
        - n    = num_atoms * n_max_bonded
        Unused slots are written as 0 by the kernel.
        """
        self._ensure_bond_residual_scratch()
        self.prg.bond_residuals_fixed_global(
            self.queue, (self.num_atoms,), None,
            self.cl_pos,
            self.cl_bond_indices_global,
            self.cl_bond_lengths,
            self.cl_bond_abs_residual,
            np.int32(self.num_atoms)
        )
        cl.enqueue_copy(self.queue, self._host_bond_abs_residual, self.cl_bond_abs_residual).wait()
        r = self._host_bond_abs_residual
        linf = float(np.max(r))
        l2 = float(np.sqrt(np.mean(r * r) + float(eps)))
        return {"linf": linf, "l2": l2, "n": int(r.size)}

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
        """Default/global tiled Jacobi solver (one CL iteration per outer loop)."""
        for _ in range(int(iterations)):
            self._seed_prev_from_curr()
            self.prg.solve_cluster_jacobi_step(
                self.queue, (self.num_groups * self.group_size,), (self.group_size,),
                self.cl_pos, self.cl_prev_pos, self.cl_momentum, self.cl_params,
                self.cl_bond_indices_local, self.cl_bond_lengths, self.cl_bond_stiffness,
                self.cl_ghost_indices, self.cl_ghost_counts,
                np.int32(self.num_atoms),
                np.float32(dt), np.float32(k_coll), np.float32(omega)
            )

    def solve_cluster_jacobi_local(self, dt=0.01, inner_iters=10, k_coll=1000.0, omega=0.8, momentum_beta=0.0):
        """Local multi-iteration tiled Jacobi solver (specialized).

        WARNING: This runs `inner_iters` iterations inside one kernel launch and therefore
        does NOT synchronize with other work-groups between iterations. Use only for
        isolated systems where cross-cluster collisions are disabled/irrelevant.
        Momentum is synchronized with global `cl_prev_pos` at kernel entry/exit.
        """
        self._seed_prev_from_curr()
        self.prg.solve_cluster_jacobi(
            self.queue, (self.num_groups * self.group_size,), (self.group_size,),
            self.cl_pos, self.cl_prev_pos, self.cl_momentum, self.cl_params,
            self.cl_bond_indices_local, self.cl_bond_lengths, self.cl_bond_stiffness,
            self.cl_ghost_indices, self.cl_ghost_counts,
            np.int32(self.num_atoms), np.int32(inner_iters),
            np.float32(dt), np.float32(k_coll), np.float32(omega)
        )

    def solve_cluster_jacobi_local_nocoll(self, dt=0.01, inner_iters=10, omega=0.8, momentum_beta=0.0):
        """Local multi-iteration Jacobi WITHOUT collisions (specialized).

        WARNING: Intended only for isolated small systems (e.g. single molecule) where
        collisions are disabled and where the system fits within a work-group.
        """
        self._seed_prev_from_curr()
        self.prg.solve_cluster_jacobi_nocoll(
            self.queue, (self.num_groups * self.group_size,), (self.group_size,),
            self.cl_pos, self.cl_prev_pos, self.cl_momentum, self.cl_params,
            self.cl_bond_indices_local, self.cl_bond_lengths, self.cl_bond_stiffness,
            np.int32(self.num_atoms), np.int32(inner_iters),
            np.float32(dt), np.float32(omega)
        )

    def solve_cluster_jacobi_step(self, dt=0.01, k_coll=1000.0, omega=0.8, momentum_beta=0.0):
        """One Jacobi iteration using persistent prev positions (for benchmarking/trajectory recording)."""
        self._seed_prev_from_curr()
        self.prg.solve_cluster_jacobi_step(
            self.queue, (self.num_groups * self.group_size,), (self.group_size,),
            self.cl_pos, self.cl_prev_pos, self.cl_momentum, self.cl_params,
            self.cl_bond_indices_local, self.cl_bond_lengths, self.cl_bond_stiffness,
            self.cl_ghost_indices, self.cl_ghost_counts,
            np.int32(self.num_atoms),
            np.float32(dt), np.float32(k_coll), np.float32(omega)
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
