import os
import numpy as np
import pyopencl as cl


def build_neighs_bk_from_bonds(n, bonds, max_deg=4):
    neighs = np.full((n, max_deg), -1, dtype=np.int32)
    bks = np.full((n, max_deg), -1, dtype=np.int32)
    deg = np.zeros((n,), dtype=np.int32)
    for (i, j) in bonds:
        if deg[i] >= max_deg or deg[j] >= max_deg:
            raise RuntimeError(f"Degree>={max_deg} for bond {i}-{j}")
        si = int(deg[i]); sj = int(deg[j])
        neighs[i, si] = j
        neighs[j, sj] = i
        bks[i, si] = sj
        bks[j, sj] = si
        deg[i] += 1
        deg[j] += 1
    return neighs, bks


def make_bk_slots_clustered(neighs, *, group_size, nnode_per_group, natoms=None):
    if natoms is None:
        natoms = int(neighs.shape[0])
    group_size = int(group_size)
    nnode_per_group = int(nnode_per_group)
    if (natoms % group_size) != 0:
        raise ValueError(f'make_bk_slots_clustered: natoms={natoms} not multiple of group_size={group_size}')
    ng = natoms // group_size
    if nnode_per_group < 0 or nnode_per_group > group_size:
        raise ValueError(f'make_bk_slots_clustered: nnode_per_group={nnode_per_group} out of range [0,{group_size}]')

    neighs = np.array(neighs, dtype=np.int32, copy=False)
    if neighs.shape != (natoms, 4):
        raise ValueError(f'make_bk_slots_clustered: neighs.shape={neighs.shape} expected ({natoms},4)')

    bkSlots = np.full((natoms, 4), -1, dtype=np.int32)
    bkCount = np.zeros((natoms,), dtype=np.int32)
    for ig in range(ng):
        abase = ig * group_size
        inode_base = ig * nnode_per_group
        for il in range(nnode_per_group):
            ia = abase + il
            inode = inode_base + il
            for k in range(4):
                ja = int(neighs[ia, k])
                if ja < 0:
                    continue
                s = int(bkCount[ja])
                if s >= 4:
                    raise RuntimeError(f"bkSlots overflow: atom {ja} has >4 back slots (from group {ig} node {il})")
                bkSlots[ja, s] = inode * 4 + k
                bkCount[ja] += 1
    return bkSlots


def make_exclusions_1st_2nd(neighs, *, max2=4):
    """Build (excl1,excl2) int4 arrays for each atom.

    excl1: up to 4 first neighbors (from neighs)
    excl2: up to 4 second neighbors (neighbors of neighbors, excluding self and first neighbors)

    This is a debug/smoke-test helper. It is not a full chemical exclusion builder.
    """
    neighs = np.asarray(neighs, dtype=np.int32)
    n = int(neighs.shape[0])
    if neighs.shape != (n, 4):
        raise ValueError(f"make_exclusions_1st_2nd: neighs.shape={neighs.shape} expected ({n},4)")

    excl1 = np.full((n, 4), -1, dtype=np.int32)
    excl2 = np.full((n, 4), -1, dtype=np.int32)

    for i in range(n):
        n1 = [int(x) for x in neighs[i] if int(x) >= 0]
        n1 = n1[:4]
        for k, j in enumerate(n1):
            excl1[i, k] = j

        s2 = []
        s1set = set(n1)
        for j in n1:
            for k in range(4):
                t = int(neighs[j, k])
                if t < 0 or t == i:
                    continue
                if t in s1set:
                    continue
                if t in s2:
                    continue
                s2.append(t)
                if len(s2) >= max2:
                    break
            if len(s2) >= max2:
                break
        for k, t in enumerate(s2[:4]):
            excl2[i, k] = int(t)

    return excl1, excl2


def make_revSlot_clustered(neighs, *, group_size, nnode_per_group, natoms=None):
    if natoms is None:
        natoms = int(neighs.shape[0])
    group_size = int(group_size)
    nnode_per_group = int(nnode_per_group)
    if (natoms % group_size) != 0:
        raise ValueError(f"make_revSlot_clustered: natoms={natoms} not multiple of group_size={group_size}")
    neighs = np.asarray(neighs, dtype=np.int32)
    if neighs.shape != (natoms, 4):
        raise ValueError(f"make_revSlot_clustered: neighs.shape={neighs.shape} expected ({natoms},4)")

    ng = natoms // group_size
    nnode_tot = ng * nnode_per_group
    revSlot = np.full((nnode_tot * 4,), -1, dtype=np.int32)

    for ig in range(ng):
        abase = ig * group_size
        inode_base = ig * nnode_per_group
        for il in range(nnode_per_group):
            ia = abase + il
            inode = inode_base + il
            i4 = inode * 4
            for k in range(4):
                ja = int(neighs[ia, k])
                if ja < 0:
                    continue
                jl = ja % group_size
                if jl >= nnode_per_group:
                    continue
                kk = -1
                for t in range(4):
                    if int(neighs[ja, t]) == ia:
                        kk = t
                        break
                if kk < 0:
                    raise RuntimeError(f"make_revSlot_clustered: missing reciprocal neigh: {ia} -> {ja} (slot {k})")
                jnode = (ja // group_size) * nnode_per_group + jl
                revSlot[i4 + k] = jnode * 4 + kk

    return revSlot


_Unbound = object()


class RRsp3:
    def __init__(self, num_atoms, group_size=64, prefer_gpu=True, device_idx=0, max_ghosts=128, build_options=None):
        self.num_atoms = int(num_atoms)
        self.group_size = int(group_size)
        if (self.num_atoms % self.group_size) != 0:
            raise ValueError(f"RRsp3: num_atoms must be multiple of group_size; got num_atoms={self.num_atoms} group_size={self.group_size}")

        self.num_groups = self.num_atoms // self.group_size
        self.max_ghosts = int(max_ghosts)

        self.ctx = self._make_context(prefer_gpu=prefer_gpu, device_idx=device_idx)
        self.queue = cl.CommandQueue(self.ctx)

        curr_dir = os.path.dirname(os.path.abspath(__file__))
        with open(os.path.join(curr_dir, "RRsp3.cl"), "r") as f:
            src = f.read()
        if build_options is None:
            build_options = []
        self.build_options = list(build_options)
        self.build_options.append(f"-DGROUP_SIZE={self.group_size}")
        self.build_options.append(f"-DMAX_GHOSTS={self.max_ghosts}")
        self.prg = cl.Program(self.ctx, src).build(options=self.build_options)
        self._kernels = {}

        mf = cl.mem_flags
        n = self.num_atoms
        ng = self.num_groups

        self.cl_pos = cl.Buffer(self.ctx, mf.READ_WRITE, n * 16)   # float4
        self.cl_quat = cl.Buffer(self.ctx, mf.READ_WRITE, n * 16)  # float4
        self.cl_radius = cl.Buffer(self.ctx, mf.READ_ONLY, n * 4)  # float

        self.cl_neighs = cl.Buffer(self.ctx, mf.READ_ONLY, n * 16)   # int4 global
        self.cl_excl1 = cl.Buffer(self.ctx, mf.READ_ONLY, n * 16)    # int4 global
        self.cl_excl2 = cl.Buffer(self.ctx, mf.READ_ONLY, n * 16)    # int4 global

        self.cl_bboxes_min = cl.Buffer(self.ctx, mf.READ_WRITE, ng * 16)
        self.cl_bboxes_max = cl.Buffer(self.ctx, mf.READ_WRITE, ng * 16)

        self.cl_ghost_indices = cl.Buffer(self.ctx, mf.READ_WRITE, ng * self.max_ghosts * 4)
        self.cl_ghost_counts = cl.Buffer(self.ctx, mf.READ_WRITE, ng * 4)

        self.cl_neighs_local = cl.Buffer(self.ctx, mf.READ_WRITE, n * 16)
        self.cl_excl1_local = cl.Buffer(self.ctx, mf.READ_WRITE, n * 16)
        self.cl_excl2_local = cl.Buffer(self.ctx, mf.READ_WRITE, n * 16)

        self.cl_dpos_coll = cl.Buffer(self.ctx, mf.READ_WRITE, n * 16)

        # Solver momentum buffers (heavy-ball)
        self.cl_dpos_mom = cl.Buffer(self.ctx, mf.READ_WRITE, n * 16)
        self.cl_dquat_mom = cl.Buffer(self.ctx, mf.READ_WRITE, n * 16)
        self.reset_momentum()

        # Dynamics buffers (linear + angular velocities, and previous state)
        self.cl_vel = cl.Buffer(self.ctx, mf.READ_WRITE, n * 16)       # float4 (vx,vy,vz,0)
        self.cl_omega = cl.Buffer(self.ctx, mf.READ_WRITE, n * 16)     # float4 (wx,wy,wz,0)  (only nodes are used)
        self.cl_pos_prev = cl.Buffer(self.ctx, mf.READ_WRITE, n * 16)  # float4
        self.cl_quat_prev = cl.Buffer(self.ctx, mf.READ_WRITE, n * 16) # float4
        self.reset_dynamics()

        self.cl_fixmask = cl.Buffer(self.ctx, mf.READ_ONLY, n * 4)   # int
        cl.enqueue_fill_buffer(self.queue, self.cl_fixmask, np.int32(0), 0, n * 4).wait()

        self._nnode_per_group = 0
        self._nnode_tot = 0
        self.cl_dpos_node = None
        self.cl_drot_node = None
        self.cl_dpos_neigh = None
        self.cl_tips = None
        self.cl_revSlot = None
        self.cl_port_local = None
        self.cl_Kflat = None
        self.cl_bkSlots = None
        self._tips_valid = False
        self._revSlot_valid = False

    def reset_dynamics(self):
        zero = np.float32(0.0)
        cl.enqueue_fill_buffer(self.queue, self.cl_vel, zero, 0, self.num_atoms * 16)
        cl.enqueue_fill_buffer(self.queue, self.cl_omega, zero, 0, self.num_atoms * 16).wait()

    def reset_momentum(self):
        """Clear solver momentum buffers (dpos_mom, dquat_mom). Call this before starting a new time step or relaxation."""
        zero = np.float32(0.0)
        cl.enqueue_fill_buffer(self.queue, self.cl_dpos_mom, zero, 0, self.num_atoms * 16)
        cl.enqueue_fill_buffer(self.queue, self.cl_dquat_mom, zero, 0, self.num_atoms * 16).wait()

    def _make_context(self, prefer_gpu=True, device_idx=0):
        dev_type = cl.device_type.GPU if prefer_gpu else cl.device_type.ALL
        try:
            for platform in cl.get_platforms():
                devs = platform.get_devices(device_type=dev_type)
                if devs:
                    return cl.Context(devices=[devs[device_idx % len(devs)]])
        except Exception:
            pass
        return cl.create_some_context(interactive=False)

    def _ensure_node_buffers(self, nnode_per_group):
        nnode_per_group = int(nnode_per_group)
        if nnode_per_group < 0 or nnode_per_group > self.group_size:
            raise ValueError(f"_ensure_node_buffers: nnode_per_group={nnode_per_group} out of range")
        nnode_tot = self.num_groups * nnode_per_group
        if (self._nnode_tot == nnode_tot) and (self.cl_dpos_node is not None):
            return

        mf = cl.mem_flags
        self._nnode_per_group = nnode_per_group
        self._nnode_tot = nnode_tot
        self.cl_dpos_node = cl.Buffer(self.ctx, mf.READ_WRITE, nnode_tot * 16)
        self.cl_drot_node = cl.Buffer(self.ctx, mf.READ_WRITE, nnode_tot * 16)
        self.cl_dpos_neigh = cl.Buffer(self.ctx, mf.READ_WRITE, nnode_tot * 4 * 16)
        self.cl_tips = cl.Buffer(self.ctx, mf.READ_WRITE, nnode_tot * 4 * 16)
        self.cl_revSlot = cl.Buffer(self.ctx, mf.READ_ONLY, nnode_tot * 4 * 4)  # int per port
        self.cl_port_local = cl.Buffer(self.ctx, mf.READ_ONLY, nnode_tot * 4 * 16)
        self.cl_Kflat = cl.Buffer(self.ctx, mf.READ_ONLY, nnode_tot * 4 * 4)
        self.cl_bkSlots = cl.Buffer(self.ctx, mf.READ_ONLY, self.num_atoms * 16)
        self._tips_valid = False
        self._revSlot_valid = False

    def upload_state(self, pos3, inv_mass, *, quat=None, nan_padding=True):
        pos3 = np.asarray(pos3, dtype=np.float32)
        inv_mass = np.asarray(inv_mass, dtype=np.float32)
        if pos3.shape != (self.num_atoms, 3):
            raise ValueError(f"upload_state: pos3.shape={pos3.shape} expected ({self.num_atoms},3)")
        if inv_mass.shape != (self.num_atoms,):
            raise ValueError(f"upload_state: inv_mass.shape={inv_mass.shape} expected ({self.num_atoms},)")
        pos4 = np.zeros((self.num_atoms, 4), dtype=np.float32)
        pos4[:, :3] = pos3
        pos4[:, 3] = inv_mass
        if nan_padding:
            # Make padding clearly invalid to catch accidental use (inv_mass==0)
            m0 = (inv_mass <= 1e-12)
            pos4[m0, 0] = np.nan
            pos4[m0, 1] = np.nan
            pos4[m0, 2] = np.nan

        if quat is None:
            quat4 = np.zeros((self.num_atoms, 4), dtype=np.float32)
            quat4[:, 3] = 1.0
        else:
            quat4 = np.asarray(quat, dtype=np.float32)
            if quat4.shape != (self.num_atoms, 4):
                raise ValueError(f"upload_state: quat.shape={quat4.shape} expected ({self.num_atoms},4)")

        cl.enqueue_copy(self.queue, self.cl_pos, pos4)
        cl.enqueue_copy(self.queue, self.cl_quat, quat4).wait()

    def upload_fixmask(self, fixmask):
        fixmask = np.asarray(fixmask, dtype=np.int32)
        if fixmask.shape != (self.num_atoms,):
            raise ValueError(f"upload_fixmask: fixmask.shape={fixmask.shape} expected ({self.num_atoms},)")
        cl.enqueue_copy(self.queue, self.cl_fixmask, fixmask).wait()

    def upload_radius(self, radius):
        r = np.asarray(radius, dtype=np.float32)
        if r.shape != (self.num_atoms,):
            raise ValueError(f"upload_radius: radius.shape={r.shape} expected ({self.num_atoms},)")
        cl.enqueue_copy(self.queue, self.cl_radius, r).wait()

    def download_radius(self):
        return self.download(self.cl_radius, (self.num_atoms,), np.float32)

    def download_bboxes(self):
        bmin = self.download(self.cl_bboxes_min, (self.num_groups, 4), np.float32)
        bmax = self.download(self.cl_bboxes_max, (self.num_groups, 4), np.float32)
        return bmin, bmax

    def download_ghosts(self):
        gi = self.download(self.cl_ghost_indices, (self.num_groups, self.max_ghosts), np.int32)
        gc = self.download(self.cl_ghost_counts, (self.num_groups,), np.int32)
        return gi, gc

    def download_neighs(self):
        return self.download(self.cl_neighs, (self.num_atoms, 4), np.int32)

    def download_neighs_local(self):
        return self.download(self.cl_neighs_local, (self.num_atoms, 4), np.int32)

    def download_excl_local(self):
        e1 = self.download(self.cl_excl1_local, (self.num_atoms, 4), np.int32)
        e2 = self.download(self.cl_excl2_local, (self.num_atoms, 4), np.int32)
        return e1, e2

    def download_bkSlots(self):
        if self.cl_bkSlots is None:
            raise RuntimeError("download_bkSlots: bkSlots buffer not allocated; call upload_cluster_ports first")
        return self.download(self.cl_bkSlots, (self.num_atoms, 4), np.int32)

    def download_port_local(self, *, nnode_per_group):
        nnode_per_group = int(nnode_per_group)
        self._ensure_node_buffers(nnode_per_group)
        pl = self.download(self.cl_port_local, (self._nnode_tot * 4, 4), np.float32)
        return pl.reshape(self._nnode_tot, 4, 4)

    def _get_kernel(self, name, default=_Unbound):
        if name not in self._kernels:
            try:
                self._kernels[name] = getattr(self.prg, name)
            except AttributeError:
                if default is not _Unbound:
                    return default
                raise
        return self._kernels[name]

    def download_dpos_coll(self):
        return self.download(self.cl_dpos_coll, (self.num_atoms, 4), np.float32)

    def download_dpos_node(self, *, nnode_per_group):
        nnode_per_group = int(nnode_per_group)
        self._ensure_node_buffers(nnode_per_group)
        return self.download(self.cl_dpos_node, (self._nnode_tot, 4), np.float32)

    def download_drot_node(self, *, nnode_per_group):
        nnode_per_group = int(nnode_per_group)
        self._ensure_node_buffers(nnode_per_group)
        return self.download(self.cl_drot_node, (self._nnode_tot, 4), np.float32)

    def download_dpos_neigh(self, *, nnode_per_group):
        nnode_per_group = int(nnode_per_group)
        self._ensure_node_buffers(nnode_per_group)
        return self.download(self.cl_dpos_neigh, (self._nnode_tot * 4, 4), np.float32).reshape(self._nnode_tot, 4, 4)

    def upload_neighs_and_exclusions(self, neighs, excl1, excl2):
        neighs = np.asarray(neighs, dtype=np.int32)
        excl1 = np.asarray(excl1, dtype=np.int32)
        excl2 = np.asarray(excl2, dtype=np.int32)
        if neighs.shape != (self.num_atoms, 4):
            raise ValueError(f"upload_neighs: neighs.shape={neighs.shape}")
        if excl1.shape != (self.num_atoms, 4) or excl2.shape != (self.num_atoms, 4):
            raise ValueError(f"upload_excl: excl1.shape={excl1.shape} excl2.shape={excl2.shape}")
        cl.enqueue_copy(self.queue, self.cl_neighs, neighs)
        cl.enqueue_copy(self.queue, self.cl_excl1, excl1)
        cl.enqueue_copy(self.queue, self.cl_excl2, excl2).wait()

    def upload_cluster_ports(self, port_local_atoms, K_atoms, *, nnode_per_group):
        """Upload node-only port geometry and stiffness.

        port_local_atoms: (natoms,4,4) float32, only first nnode_per_group per group are used.
        K_atoms: (natoms,4) float32, only nodes used.
        """
        nnode_per_group = int(nnode_per_group)
        self._ensure_node_buffers(nnode_per_group)

        plA = np.asarray(port_local_atoms, dtype=np.float32)
        kkA = np.asarray(K_atoms, dtype=np.float32)
        if plA.shape != (self.num_atoms, 4, 4):
            raise ValueError(f"upload_cluster_ports: port_local_atoms.shape={plA.shape} expected ({self.num_atoms},4,4)")
        if kkA.shape != (self.num_atoms, 4):
            raise ValueError(f"upload_cluster_ports: K_atoms.shape={kkA.shape} expected ({self.num_atoms},4)")

        nnode_tot = self._nnode_tot
        pl = np.zeros((nnode_tot, 4, 4), dtype=np.float32)
        kk = np.zeros((nnode_tot, 4), dtype=np.float32)

        for ig in range(self.num_groups):
            abase = ig * self.group_size
            inode_base = ig * nnode_per_group
            for il in range(nnode_per_group):
                ia = abase + il
                inode = inode_base + il
                pl[inode, :, :] = plA[ia, :, :]
                kk[inode, :] = kkA[ia, :]

        cl.enqueue_copy(self.queue, self.cl_port_local, pl.reshape(nnode_tot * 4, 4))
        cl.enqueue_copy(self.queue, self.cl_Kflat, kk.reshape(nnode_tot * 4)).wait()

    def upload_bkSlots(self, bkSlots):
        if self.cl_bkSlots is None:
            raise RuntimeError("upload_bkSlots: call upload_cluster_ports first (allocates bkSlots buffer)")
        bk = np.asarray(bkSlots, dtype=np.int32)
        if bk.shape != (self.num_atoms, 4):
            raise ValueError(f"upload_bkSlots: bkSlots.shape={bk.shape} expected ({self.num_atoms},4)")
        cl.enqueue_copy(self.queue, self.cl_bkSlots, bk).wait()

    def upload_revSlot(self, revSlot, *, nnode_per_group):
        nnode_per_group = int(nnode_per_group)
        self._ensure_node_buffers(nnode_per_group)
        if self.cl_revSlot is None:
            raise RuntimeError("upload_revSlot: revSlot buffer not allocated; call upload_cluster_ports first")
        rs = np.asarray(revSlot, dtype=np.int32)
        if rs.shape != (self._nnode_tot * 4,):
            raise ValueError(f"upload_revSlot: revSlot.shape={rs.shape} expected ({self._nnode_tot * 4},)")
        cl.enqueue_copy(self.queue, self.cl_revSlot, rs).wait()
        self._revSlot_valid = True

    def run_bboxes_and_topology(self, *, bbox_margin=0.5):
        """Run only broad-phase kernels to populate bboxes and ghost lists (no integration).

        This is used by the visual debugger to show AABBs + halo mapping immediately after upload.
        """
        natoms = np.int32(self.num_atoms)
        ng = np.int32(self.num_groups)
        bbox_margin_f = np.float32(float(bbox_margin))

        rad = self.download(self.cl_radius, (self.num_atoms,), np.float32)
        rmax = float(np.max(rad)) if rad.size else 0.0
        margin_sq = np.float32((2.0 * rmax + float(bbox_margin)) ** 2)

        global_size = (self.num_groups * self.group_size,)
        local_size = (self.group_size,)

        self._get_kernel('update_bboxes_rigid')(
            self.queue, global_size, local_size,
            self.cl_pos, self.cl_radius,
            self.cl_bboxes_min, self.cl_bboxes_max,
            cl.LocalMemory(self.group_size * 16),
            cl.LocalMemory(self.group_size * 16),
            natoms
        )

        self._get_kernel('build_local_topology_rigid')(
            self.queue, global_size, local_size,
            self.cl_pos,
            self.cl_bboxes_min, self.cl_bboxes_max,
            self.cl_neighs,
            self.cl_excl1, self.cl_excl2,
            self.cl_ghost_indices, self.cl_ghost_counts,
            self.cl_neighs_local,
            self.cl_excl1_local, self.cl_excl2_local,
            natoms, ng,
            margin_sq, bbox_margin_f
        ).wait()

    def step_cluster(self, *, nnode_per_group, dt=0.1, k_coll=50.0, relaxation=0.5, bbox_margin=0.5, momentum_beta=0.0, port_kernel='current', rot_mass_scale=1.0, n_rot_substeps=5, rot_eps=0.0, theta_max=0.0):
        nnode_per_group = int(nnode_per_group)
        self._ensure_node_buffers(nnode_per_group)

        natoms = np.int32(self.num_atoms)
        ng = np.int32(self.num_groups)
        dt_f = np.float32(float(dt))
        relax_f = np.float32(float(relaxation))
        kcoll_f = np.float32(float(k_coll))
        bbox_margin_f = np.float32(float(bbox_margin))
        beta_f = np.float32(float(momentum_beta))
        rms_f = np.float32(float(rot_mass_scale))
        rot_eps_f = np.float32(float(rot_eps))
        theta_max_f = np.float32(float(theta_max))

        # heuristic ghost margin
        rad = self.download(self.cl_radius, (self.num_atoms,), np.float32)
        rmax = float(np.max(rad)) if rad.size else 0.0
        margin_sq = np.float32((2.0 * rmax + float(bbox_margin)) ** 2)

        global_size = (self.num_groups * self.group_size,)
        local_size = (self.group_size,)

        self._get_kernel('update_bboxes_rigid')(
            self.queue, global_size, local_size,
            self.cl_pos, self.cl_radius,
            self.cl_bboxes_min, self.cl_bboxes_max,
            cl.LocalMemory(self.group_size * 16),
            cl.LocalMemory(self.group_size * 16),
            natoms
        )

        self._get_kernel('build_local_topology_rigid')(
            self.queue, global_size, local_size,
            self.cl_pos,
            self.cl_bboxes_min, self.cl_bboxes_max,
            self.cl_neighs,
            self.cl_excl1, self.cl_excl2,
            self.cl_ghost_indices, self.cl_ghost_counts,
            self.cl_neighs_local,
            self.cl_excl1_local, self.cl_excl2_local,
            natoms, ng,
            margin_sq, bbox_margin_f
        )

        cl.enqueue_fill_buffer(self.queue, self.cl_dpos_coll, np.float32(0.0), 0, self.num_atoms * 16)
        cl.enqueue_fill_buffer(self.queue, self.cl_dpos_node, np.float32(0.0), 0, self._nnode_tot * 16)
        cl.enqueue_fill_buffer(self.queue, self.cl_drot_node, np.float32(0.0), 0, self._nnode_tot * 16)
        cl.enqueue_fill_buffer(self.queue, self.cl_dpos_neigh, np.float32(0.0), 0, self._nnode_tot * 4 * 16)

        self._get_kernel('compute_collision_cluster_rigid')(
            self.queue, global_size, local_size,
            self.cl_pos, self.cl_radius,
            self.cl_excl1_local, self.cl_excl2_local,
            self.cl_ghost_indices, self.cl_ghost_counts,
            self.cl_dpos_coll,
            natoms,
            kcoll_f
        )

        if port_kernel in (None, 'current', 'rigid'):
            kfun = self._get_kernel('compute_ports_cluster_rigid', default=None)
            massless_rot = 0
        elif port_kernel in ('orid', 'orig', 'original'):
            # Backward-compat: kernel was renamed from _orid -> _orig
            kfun = self._get_kernel('compute_ports_cluster_rigid_orig', default=None)
            if kfun is None:
                kfun = self._get_kernel('compute_ports_cluster_rigid_orid', default=None)
            massless_rot = 0
        elif port_kernel in ('substep', 'substep_optimized'):
            kfun = self._get_kernel('compute_ports_cluster_rigid_substep_optimized', default=None)
            massless_rot = 1
        elif port_kernel in ('shapematch', 'shape_match'):
            kfun = self._get_kernel('compute_ports_cluster_rigid_shapematch', default=None)
            massless_rot = 1
        elif port_kernel in ('eigen', 'q_eigen'):
            kfun = self._get_kernel('compute_ports_cluster_rigid_eigen', default=None)
            massless_rot = 1
        else:
            raise ValueError(f"RRsp3.step_cluster: unknown port_kernel={port_kernel!r}; use 'current','orig','substep_optimized','shapematch','eigen'")
        if kfun is None:
            raise RuntimeError(f"RRsp3.step_cluster: requested port_kernel={port_kernel!r} not present in built program")

        if port_kernel in ('eigen', 'q_eigen'):
            # Two-pass massless eigen with symmetric tips:
            # If tips not precomputed (first step or buffer recreated), compute from current state
            if not self._tips_valid:
                self._get_kernel('compute_tips')(
                    self.queue, (self.num_atoms,), None,
                    natoms, np.int32(nnode_per_group),
                    self.cl_pos, self.cl_quat,
                    self.cl_port_local,
                    self.cl_tips
                ).wait()
                self._tips_valid = True
            # Pass 1: compute optimal rotations (geometric solve)
            self._get_kernel('compute_optimal_rotation_eigen')(
                self.queue, global_size, local_size,
                self.cl_pos, self.cl_quat, self.cl_radius,
                self.cl_neighs_local,
                self.cl_ghost_indices, self.cl_ghost_counts,
                self.cl_port_local,
                self.cl_Kflat,
                self.cl_drot_node,
                self.cl_dquat_mom,  # temp storage for quat_opt (unused now, kept for compat)
                natoms, np.int32(nnode_per_group),
                dt_f
            )
            # Pass 2: symmetric recoil using precomputed tips (reads tips, writes dpos)
            self._get_kernel('compute_ports_cluster_rigid_eigen_tips')(
                self.queue, global_size, local_size,
                self.cl_pos,
                self.cl_neighs_local,
                self.cl_ghost_indices, self.cl_ghost_counts,
                self.cl_Kflat,
                self.cl_tips,
                self.cl_dpos_node,
                self.cl_dpos_neigh,
                natoms, np.int32(nnode_per_group),
                dt_f
            )
        else:
            args = [
                self.queue, global_size, local_size,
                self.cl_pos, self.cl_quat, self.cl_radius,
                self.cl_neighs_local,
                self.cl_ghost_indices, self.cl_ghost_counts,
                self.cl_port_local,
                self.cl_Kflat,
                self.cl_dpos_node,
                self.cl_drot_node,
                self.cl_dpos_neigh,
                natoms, np.int32(nnode_per_group),
                dt_f, np.int32(0)
            ]
            if port_kernel in (None, 'current', 'rigid'):
                args.append(rms_f)
            elif port_kernel in ('substep', 'substep_optimized'):
                args.append(np.int32(int(n_rot_substeps)))
                args.append(rot_eps_f)
                args.append(theta_max_f)
            # orig kernel now takes two extra trailing args
            if port_kernel in ('orid', 'orig', 'original'):
                args.append(None)      # quat_opt = null
                args.append(np.int32(0))  # skip_rotation = 0
            kfun(*args)

        # apply_corrections now takes optional port_local and tips for massless_rot tip caching
        apply_args = [
            self.queue, (self.num_atoms,), None,
            natoms, np.int32(nnode_per_group),
            self.cl_pos, self.cl_quat,
            self.cl_fixmask,
            self.cl_bkSlots,
            self.cl_dpos_node, self.cl_drot_node, self.cl_dpos_neigh,
            self.cl_dpos_coll,
            self.cl_dpos_mom, self.cl_dquat_mom,
            relax_f, beta_f, np.int32(massless_rot)
        ]
        if massless_rot:
            apply_args.append(self.cl_port_local)
            apply_args.append(self.cl_tips)
        else:
            apply_args.append(None)
            apply_args.append(None)
        self._get_kernel('apply_corrections_rigid_ports')(*apply_args).wait()

    def step_dynamics(self, *, nnode_per_group, dt=0.1, k_coll=50.0, relaxation=1.0, bbox_margin=0.5, port_kernel='current', rot_mass_scale=1.0, n_rot_substeps=5, rot_eps=0.0, theta_max=0.0, damp=1.0):
        """Leapfrog/PBD-style dynamics with rotational DOFs.

        Scheme:
        1) Predict pos/quat using stored vel/omega (semi-implicit / leapfrog style)
        2) Project constraints using existing RRsp3 kernels (collision + ports + apply_corrections)
        3) Update vel/omega from (new - old)/dt

        This keeps NaN-padding diagnostics intact: padding atoms (invM==0) remain NaN.
        Pinning is via fixmask, not invm.
        """
        nnode_per_group = int(nnode_per_group)
        self._ensure_node_buffers(nnode_per_group)

        natoms = np.int32(self.num_atoms)
        ng = np.int32(self.num_groups)
        dt_f = np.float32(float(dt))
        relax_f = np.float32(float(relaxation))
        kcoll_f = np.float32(float(k_coll))
        bbox_margin_f = np.float32(float(bbox_margin))
        rms_f = np.float32(float(rot_mass_scale))
        rot_eps_f = np.float32(float(rot_eps))
        theta_max_f = np.float32(float(theta_max))
        damp_f = np.float32(float(damp))

        # Clear solver momentum buffers to ensure dynamics does not depend on previous relaxation state
        self.reset_momentum()

        # Predict (store previous state inside kernel)
        self._tips_valid = False
        self._get_kernel('predict_dynamics')(
            self.queue, (self.num_atoms,), None,
            natoms, np.int32(nnode_per_group),
            self.cl_pos, self.cl_quat,
            self.cl_fixmask,
            self.cl_vel, self.cl_omega,
            self.cl_pos_prev, self.cl_quat_prev,
            dt_f
        ).wait()

        # heuristic ghost margin
        rad = self.download(self.cl_radius, (self.num_atoms,), np.float32)
        rmax = float(np.max(rad)) if rad.size else 0.0
        margin_sq = np.float32((2.0 * rmax + float(bbox_margin)) ** 2)

        global_size = (self.num_groups * self.group_size,)
        local_size = (self.group_size,)

        self._get_kernel('update_bboxes_rigid')(
            self.queue, global_size, local_size,
            self.cl_pos, self.cl_radius,
            self.cl_bboxes_min, self.cl_bboxes_max,
            cl.LocalMemory(self.group_size * 16),
            cl.LocalMemory(self.group_size * 16),
            natoms
        )

        self._get_kernel('build_local_topology_rigid')(
            self.queue, global_size, local_size,
            self.cl_pos,
            self.cl_bboxes_min, self.cl_bboxes_max,
            self.cl_neighs,
            self.cl_excl1, self.cl_excl2,
            self.cl_ghost_indices, self.cl_ghost_counts,
            self.cl_neighs_local,
            self.cl_excl1_local, self.cl_excl2_local,
            natoms, ng,
            margin_sq, bbox_margin_f
        )

        cl.enqueue_fill_buffer(self.queue, self.cl_dpos_coll, np.float32(0.0), 0, self.num_atoms * 16)
        cl.enqueue_fill_buffer(self.queue, self.cl_dpos_node, np.float32(0.0), 0, self._nnode_tot * 16)
        cl.enqueue_fill_buffer(self.queue, self.cl_drot_node, np.float32(0.0), 0, self._nnode_tot * 16)
        cl.enqueue_fill_buffer(self.queue, self.cl_dpos_neigh, np.float32(0.0), 0, self._nnode_tot * 4 * 16)

        self._get_kernel('compute_collision_cluster_rigid')(
            self.queue, global_size, local_size,
            self.cl_pos, self.cl_radius,
            self.cl_excl1_local, self.cl_excl2_local,
            self.cl_ghost_indices, self.cl_ghost_counts,
            self.cl_dpos_coll,
            natoms,
            kcoll_f
        )

        # ports (same selection logic as step_cluster)
        if port_kernel in (None, 'current', 'rigid'):
            kfun = self._get_kernel('compute_ports_cluster_rigid', default=None)
            massless_rot = 0
        elif port_kernel in ('orid', 'orig', 'original'):
            kfun = self._get_kernel('compute_ports_cluster_rigid_orig', default=None)
            if kfun is None:
                kfun = self._get_kernel('compute_ports_cluster_rigid_orid', default=None)
            massless_rot = 0
        elif port_kernel in ('substep', 'substep_optimized'):
            kfun = self._get_kernel('compute_ports_cluster_rigid_substep_optimized', default=None)
            massless_rot = 1
        elif port_kernel in ('shapematch', 'shape_match'):
            kfun = self._get_kernel('compute_ports_cluster_rigid_shapematch', default=None)
            massless_rot = 1
        elif port_kernel in ('eigen', 'q_eigen'):
            kfun = self._get_kernel('compute_ports_cluster_rigid_eigen', default=None)
            massless_rot = 1
        else:
            raise ValueError(f"RRsp3.step_dynamics: unknown port_kernel={port_kernel!r}; use 'current','orig','substep_optimized','shapematch','eigen'")
        if kfun is None:
            raise RuntimeError(f"RRsp3.step_dynamics: requested port_kernel={port_kernel!r} not present in built program")

        if port_kernel in ('eigen', 'q_eigen'):
            if not self._tips_valid:
                self._get_kernel('compute_tips')(
                    self.queue, (self.num_atoms,), None,
                    natoms, np.int32(nnode_per_group),
                    self.cl_pos, self.cl_quat,
                    self.cl_port_local,
                    self.cl_tips
                ).wait()
                self._tips_valid = True
            self._get_kernel('compute_optimal_rotation_eigen')(
                self.queue, global_size, local_size,
                self.cl_pos, self.cl_quat, self.cl_radius,
                self.cl_neighs_local,
                self.cl_ghost_indices, self.cl_ghost_counts,
                self.cl_port_local,
                self.cl_Kflat,
                self.cl_drot_node,
                self.cl_dquat_mom,
                natoms, np.int32(nnode_per_group),
                dt_f
            )
            self._get_kernel('compute_ports_cluster_rigid_eigen_tips')(
                self.queue, global_size, local_size,
                self.cl_pos,
                self.cl_neighs_local,
                self.cl_ghost_indices, self.cl_ghost_counts,
                self.cl_Kflat,
                self.cl_tips,
                self.cl_dpos_node,
                self.cl_dpos_neigh,
                natoms, np.int32(nnode_per_group),
                dt_f
            )
        else:
            args = [
                self.queue, global_size, local_size,
                self.cl_pos, self.cl_quat, self.cl_radius,
                self.cl_neighs_local,
                self.cl_ghost_indices, self.cl_ghost_counts,
                self.cl_port_local,
                self.cl_Kflat,
                self.cl_dpos_node,
                self.cl_drot_node,
                self.cl_dpos_neigh,
                natoms, np.int32(nnode_per_group),
                dt_f, np.int32(0)
            ]
            if port_kernel in (None, 'current', 'rigid'):
                args.append(rms_f)
            elif port_kernel in ('substep', 'substep_optimized'):
                args.append(np.int32(int(n_rot_substeps)))
                args.append(rot_eps_f)
                args.append(theta_max_f)
            if port_kernel in ('orid', 'orig', 'original'):
                args.append(None)
                args.append(np.int32(0))
            kfun(*args)

        apply_args = [
            self.queue, (self.num_atoms,), None,
            natoms, np.int32(nnode_per_group),
            self.cl_pos, self.cl_quat,
            self.cl_fixmask,
            self.cl_bkSlots,
            self.cl_dpos_node, self.cl_drot_node, self.cl_dpos_neigh,
            self.cl_dpos_coll,
            self.cl_dpos_mom, self.cl_dquat_mom,
            relax_f, np.float32(0.0), np.int32(massless_rot)
        ]
        if massless_rot:
            apply_args.append(self.cl_port_local)
            apply_args.append(self.cl_tips)
        else:
            apply_args.append(None)
            apply_args.append(None)
        self._get_kernel('apply_corrections_rigid_ports')(*apply_args).wait()

        # Update velocities from pos/quat deltas
        self._get_kernel('update_velocities_dynamics')(
            self.queue, (self.num_atoms,), None,
            natoms, np.int32(nnode_per_group),
            self.cl_pos, self.cl_quat,
            self.cl_fixmask,
            self.cl_vel, self.cl_omega,
            self.cl_pos_prev, self.cl_quat_prev,
            dt_f, damp_f
        ).wait()

    def compute_cluster_deltas(self, *, nnode_per_group, dt=0.1, k_coll=50.0, bbox_margin=0.5, skip_ports=False):
        nnode_per_group = int(nnode_per_group)
        if not skip_ports:
            self._ensure_node_buffers(nnode_per_group)

        natoms = np.int32(self.num_atoms)
        ng = np.int32(self.num_groups)
        dt_f = np.float32(float(dt))
        kcoll_f = np.float32(float(k_coll))
        bbox_margin_f = np.float32(float(bbox_margin))

        rad = self.download(self.cl_radius, (self.num_atoms,), np.float32)
        rmax = float(np.max(rad)) if rad.size else 0.0
        margin_sq = np.float32((2.0 * rmax + float(bbox_margin)) ** 2)

        global_size = (self.num_groups * self.group_size,)
        local_size = (self.group_size,)

        self._get_kernel('update_bboxes_rigid')(
            self.queue, global_size, local_size,
            self.cl_pos, self.cl_radius,
            self.cl_bboxes_min, self.cl_bboxes_max,
            cl.LocalMemory(self.group_size * 16),
            cl.LocalMemory(self.group_size * 16),
            natoms
        )

        self._get_kernel('build_local_topology_rigid')(
            self.queue, global_size, local_size,
            self.cl_pos,
            self.cl_bboxes_min, self.cl_bboxes_max,
            self.cl_neighs,
            self.cl_excl1, self.cl_excl2,
            self.cl_ghost_indices, self.cl_ghost_counts,
            self.cl_neighs_local,
            self.cl_excl1_local, self.cl_excl2_local,
            natoms, ng,
            margin_sq, bbox_margin_f
        )

        cl.enqueue_fill_buffer(self.queue, self.cl_dpos_coll, np.float32(0.0), 0, self.num_atoms * 16)
        if not skip_ports:
            cl.enqueue_fill_buffer(self.queue, self.cl_dpos_node, np.float32(0.0), 0, self._nnode_tot * 16)
            cl.enqueue_fill_buffer(self.queue, self.cl_drot_node, np.float32(0.0), 0, self._nnode_tot * 16)
            cl.enqueue_fill_buffer(self.queue, self.cl_dpos_neigh, np.float32(0.0), 0, self._nnode_tot * 4 * 16)

        self._get_kernel('compute_collision_cluster_rigid')(
            self.queue, global_size, local_size,
            self.cl_pos, self.cl_radius,
            self.cl_excl1_local, self.cl_excl2_local,
            self.cl_ghost_indices, self.cl_ghost_counts,
            self.cl_dpos_coll,
            natoms,
            kcoll_f
        )

        if not skip_ports:
            self._get_kernel('compute_ports_cluster_rigid')(
                self.queue, global_size, local_size,
                self.cl_pos, self.cl_quat, self.cl_radius,
                self.cl_neighs_local,
                self.cl_ghost_indices, self.cl_ghost_counts,
                self.cl_port_local,
                self.cl_Kflat,
                self.cl_dpos_node,
                self.cl_drot_node,
                self.cl_dpos_neigh,
                natoms, np.int32(nnode_per_group),
                dt_f, np.int32(0), np.float32(1.0)
            ).wait()

        dpos_coll = self.download(self.cl_dpos_coll, (self.num_atoms, 4), np.float32)
        if skip_ports:
            return dpos_coll, None, None, None
        dpos_node = self.download(self.cl_dpos_node, (self._nnode_tot, 4), np.float32)
        drot_node = self.download(self.cl_drot_node, (self._nnode_tot, 4), np.float32)
        dpos_neigh = self.download(self.cl_dpos_neigh, (self._nnode_tot * 4, 4), np.float32)
        return dpos_coll, dpos_node, drot_node, dpos_neigh

    def apply_cluster_corrections(self, *, nnode_per_group, relaxation=0.5, momentum_beta=0.0):
        nnode_per_group = int(nnode_per_group)
        self._ensure_node_buffers(nnode_per_group)
        natoms = np.int32(self.num_atoms)
        relax_f = np.float32(float(relaxation))
        beta_f = np.float32(float(momentum_beta))
        self._get_kernel('apply_corrections_rigid_ports')(
            self.queue, (self.num_atoms,), None,
            natoms, np.int32(nnode_per_group),
            self.cl_pos, self.cl_quat,
            self.cl_fixmask,
            self.cl_bkSlots,
            self.cl_dpos_node, self.cl_drot_node, self.cl_dpos_neigh,
            self.cl_dpos_coll,
            self.cl_dpos_mom, self.cl_dquat_mom,
            relax_f, beta_f, np.int32(0),
            None, None
        ).wait()

    def download(self, buf, shape, dtype):
        out = np.empty(shape, dtype=dtype)
        cl.enqueue_copy(self.queue, out, buf).wait()
        return out

    def download_pos_quat(self):
        pos4 = self.download(self.cl_pos, (self.num_atoms, 4), np.float32)
        quat4 = self.download(self.cl_quat, (self.num_atoms, 4), np.float32)
        return pos4, quat4


def linear_momentum_from_positions(pos_old, pos_new, m, dt):
    pos_old = np.asarray(pos_old, dtype=np.float32)
    pos_new = np.asarray(pos_new, dtype=np.float32)
    m = np.asarray(m, dtype=np.float32)
    if pos_old.shape != pos_new.shape:
        raise ValueError(f"linear_momentum_from_positions: pos_old.shape={pos_old.shape} pos_new.shape={pos_new.shape}")
    if pos_old.shape[1] != 3:
        raise ValueError(f"linear_momentum_from_positions: expected pos[:,3], got {pos_old.shape}")
    if m.shape != (pos_old.shape[0],):
        raise ValueError(f"linear_momentum_from_positions: m.shape={m.shape} expected ({pos_old.shape[0]},)")
    dt = float(dt)
    if dt <= 0:
        raise ValueError(f"linear_momentum_from_positions: dt={dt} must be >0")
    v = (pos_new - pos_old) / dt
    return np.sum(v * m[:, None], axis=0)


def angular_momentum_from_positions(pos_old, pos_new, m, dt, origin=None):
    pos_old = np.asarray(pos_old, dtype=np.float32)
    pos_new = np.asarray(pos_new, dtype=np.float32)
    m = np.asarray(m, dtype=np.float32)
    dt = float(dt)
    if origin is None:
        origin = np.zeros((3,), dtype=np.float32)
    origin = np.asarray(origin, dtype=np.float32)
    if origin.shape != (3,):
        raise ValueError(f"angular_momentum_from_positions: origin.shape={origin.shape} expected (3,)")
    v = (pos_new - pos_old) / dt
    r = pos_new - origin[None, :]
    p = v * m[:, None]
    return np.sum(np.cross(r, p), axis=0)
