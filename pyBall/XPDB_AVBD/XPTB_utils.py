import os

import numpy as np


def invert_permutation(perm):
    perm = np.asarray(perm, dtype=np.int32)
    if perm.ndim != 1:
        raise ValueError(f"invert_permutation: perm.ndim={perm.ndim} expected 1")
    inv = np.empty_like(perm)
    inv[perm] = np.arange(perm.size, dtype=np.int32)
    return inv


def apply_permutation_to_bonds(bonds, perm_inv):
    """Map bonds (old indices) to new indices using perm_inv (new_idx_of_old)."""
    perm_inv = np.asarray(perm_inv, dtype=np.int32)
    out = []
    for (i, j) in bonds:
        i = int(i); j = int(j)
        ni = int(perm_inv[i])
        nj = int(perm_inv[j])
        if ni == nj:
            continue
        if ni < nj:
            out.append((ni, nj))
        else:
            out.append((nj, ni))
    return sorted(set(out))


def pack_molecules_contiguous(
    molecules,
    *,
    group_size=64,
    nodes_first=True,
    pad_to_group=True,
):
    """Pack multiple molecules into contiguous groups of size group_size.

    This is a host-side helper to prepare data for OpenCL kernels which assume:
    - atoms are in contiguous clusters (workgroups)
    - within each cluster: node atoms first, then caps, then optional padding

    Parameters
    ----------
    molecules : list of dict
        Each dict must provide:
        - 'elems': list[str] length n
        - 'pos': (n,3) float
        - 'bonds': list[(i,j)] within molecule indexing [0..n)
        - 'nnode': int (number of node atoms within molecule)
    group_size : int
        Workgroup size used by kernels.
    nodes_first : bool
        If True, keep node atoms first within each molecule.
    pad_to_group : bool
        If True, add dummy atoms (- padding) so each molecule occupies exactly group_size.

    Returns
    -------
    packed : dict
        - 'elems': list[str] length natoms_total
        - 'pos': (natoms_total,3) float32
        - 'group_id': (natoms_total,) int32
        - 'is_padding': (natoms_total,) bool
        - 'mol_of_atom': (natoms_total,) int32
        - 'perm_old_to_new': list[np.ndarray] per molecule
        - 'perm_new_to_old': list[np.ndarray] per molecule
        - 'bonds': list[(i,j)] in packed global indexing
        - 'nnode_group': (ngroups,) int32
        - 'natoms_total': int
        - 'group_size': int
    """
    if int(group_size) <= 0:
        raise ValueError(f"pack_molecules_contiguous: invalid group_size={group_size}")

    elems_all = []
    pos_all = []
    group_id_all = []
    is_pad_all = []
    mol_of_atom_all = []
    bonds_all = []
    nnode_group = []
    perm_old_to_new_list = []
    perm_new_to_old_list = []

    base = 0
    for imol, mol in enumerate(molecules):
        if not isinstance(mol, dict):
            raise ValueError(f"pack_molecules_contiguous: molecule[{imol}] is not dict")
        elems = mol.get('elems', None)
        pos = mol.get('pos', None)
        bonds = mol.get('bonds', None)
        nnode = int(mol.get('nnode', -1))
        if elems is None or pos is None or bonds is None:
            raise ValueError(f"pack_molecules_contiguous: molecule[{imol}] must contain elems,pos,bonds")
        pos = np.asarray(pos, dtype=np.float32)
        n = int(pos.shape[0])
        if pos.shape != (n, 3):
            raise ValueError(f"pack_molecules_contiguous: molecule[{imol}] pos.shape={pos.shape} expected (n,3)")
        if len(elems) != n:
            raise ValueError(f"pack_molecules_contiguous: molecule[{imol}] len(elems)={len(elems)} != n={n}")
        if nnode < 0 or nnode > n:
            raise ValueError(f"pack_molecules_contiguous: molecule[{imol}] nnode={nnode} out of range [0,{n}]")

        if pad_to_group and (n > int(group_size)):
            raise ValueError(f"pack_molecules_contiguous: molecule[{imol}] n={n} > group_size={group_size}; increase group_size or split molecule")

        # local reorder within molecule: nodes first, then caps
        # Rules:
        # - if mol provides explicit 'perm', we trust it (after validation)
        # - otherwise, if nodes_first, we infer node_mask from bond degree (>1), and fail if count != nnode
        # - otherwise identity
        perm_override = mol.get('perm', None)
        if perm_override is not None:
            perm = np.asarray(perm_override, dtype=np.int32)
            if perm.shape != (n,):
                raise ValueError(f"pack_molecules_contiguous: molecule[{imol}] perm.shape={perm.shape} expected ({n},)")
            if np.any(perm < 0) or np.any(perm >= n):
                raise ValueError(f"pack_molecules_contiguous: molecule[{imol}] perm has out-of-range entries")
            if len(set(int(x) for x in perm.tolist())) != n:
                raise ValueError(f"pack_molecules_contiguous: molecule[{imol}] perm is not a permutation")
        elif nodes_first:
            deg = np.zeros((n,), dtype=np.int32)
            for (i, j) in bonds:
                i = int(i); j = int(j)
                if i < 0 or i >= n or j < 0 or j >= n:
                    raise ValueError(f"pack_molecules_contiguous: molecule[{imol}] bond {(i,j)} out of range [0,{n})")
                if i == j:
                    continue
                deg[i] += 1
                deg[j] += 1
            node_mask = deg > 1
            nnode_inferred = int(np.sum(node_mask))
            if nnode_inferred != nnode:
                raise ValueError(f"pack_molecules_contiguous: molecule[{imol}] inferred nnode={nnode_inferred} from degree>1 but provided nnode={nnode}; provide mol['perm'] if you want a custom node/cap split")
            perm = np.concatenate([np.where(node_mask)[0], np.where(~node_mask)[0]]).astype(np.int32, copy=False)
        else:
            perm = np.arange(n, dtype=np.int32)

        perm_inv = invert_permutation(perm)
        perm_old_to_new_list.append(perm_inv)
        perm_new_to_old_list.append(perm)

        pos_loc = pos[perm]
        elems_loc = [elems[int(i)] for i in perm.tolist()]
        bonds_loc = apply_permutation_to_bonds(bonds, perm_inv)

        # append real atoms
        ng = int(base // group_size)
        elems_all += elems_loc
        pos_all.append(pos_loc)
        group_id_all += [ng] * n
        is_pad_all += [False] * n
        mol_of_atom_all += [imol] * n
        for (i, j) in bonds_loc:
            bonds_all.append((base + int(i), base + int(j)))

        # padding
        n_pad = 0
        if pad_to_group:
            n_pad = int(group_size) - n
        if n_pad > 0:
            elems_all += ['X'] * n_pad
            pos_all.append(np.zeros((n_pad, 3), dtype=np.float32))
            group_id_all += [ng] * n_pad
            is_pad_all += [True] * n_pad
            mol_of_atom_all += [imol] * n_pad

        nnode_group.append(nnode)
        base += n + n_pad

    pos_all = np.vstack(pos_all) if len(pos_all) else np.zeros((0, 3), dtype=np.float32)
    group_id_all = np.asarray(group_id_all, dtype=np.int32)
    is_pad_all = np.asarray(is_pad_all, dtype=bool)
    mol_of_atom_all = np.asarray(mol_of_atom_all, dtype=np.int32)
    nnode_group = np.asarray(nnode_group, dtype=np.int32)

    if pos_all.shape[0] != len(elems_all):
        raise RuntimeError(f"pack_molecules_contiguous: internal size mismatch pos={pos_all.shape[0]} elems={len(elems_all)}")
    if (pos_all.shape[0] % int(group_size)) != 0:
        raise RuntimeError(f"pack_molecules_contiguous: natoms_total={pos_all.shape[0]} not divisible by group_size={group_size}")

    return {
        'elems': elems_all,
        'pos': pos_all.astype(np.float32, copy=False),
        'group_id': group_id_all,
        'is_padding': is_pad_all,
        'mol_of_atom': mol_of_atom_all,
        'perm_old_to_new': perm_old_to_new_list,
        'perm_new_to_old': perm_new_to_old_list,
        'bonds': bonds_all,
        'nnode_group': nnode_group,
        'natoms_total': int(pos_all.shape[0]),
        'group_size': int(group_size),
    }


def as_unit(v):
    v = np.array(v, dtype=np.float32)
    n = float(np.linalg.norm(v))
    if not np.isfinite(n) or n <= 0.0:
        raise ValueError(f"direction must be non-zero finite vector, got {v}")
    return v / n


def deform_shift_atom(pos, atom_idx=0, direction=(1.0, 0.0, 0.0), shift=2.0):
    d = as_unit(direction) * float(shift)
    out = np.array(pos, dtype=np.float32, copy=True)
    if atom_idx < 0 or atom_idx >= len(out):
        raise IndexError(f"atom_idx out of range: {atom_idx} vs n={len(out)}")
    out[atom_idx] += d
    return out


def deform_scale_along_direction(pos, direction=(1.0, 0.0, 0.0), scale=1.2, origin=None, atom_indices=None):
    d = as_unit(direction)
    out = np.array(pos, dtype=np.float32, copy=True)
    if origin is None:
        origin = out.mean(axis=0)
    origin = np.array(origin, dtype=np.float32)
    if atom_indices is None:
        sel = np.arange(len(out), dtype=np.int32)
    else:
        sel = np.array(atom_indices, dtype=np.int32)
        if np.any(sel < 0) or np.any(sel >= len(out)):
            raise IndexError(f"atom_indices out of range: min={sel.min()} max={sel.max()} n={len(out)}")

    x = out[sel] - origin[None, :]
    s = np.dot(x, d)
    out[sel] = origin[None, :] + x + (float(scale) - 1.0) * s[:, None] * d[None, :]
    return out


def make_h2o_geometry(add_angle=True, k_oh=500.0, k_hh=200.0):
    """Return (pos, bonds_adj).

    bonds_adj is list-of-lists: bonds_adj[i] = [(j, L0, K), ...]

    Atom order: 0=O, 1=H1, 2=H2.
    """
    r = 0.9572
    ang = np.deg2rad(104.52)
    h1 = np.array([r, 0.0, 0.0], dtype=np.float32)
    h2 = np.array([r * np.cos(ang), r * np.sin(ang), 0.0], dtype=np.float32)
    o = np.array([0.0, 0.0, 0.0], dtype=np.float32)
    pos = np.stack([o, h1, h2], axis=0).astype(np.float32)

    def dist(i, j):
        return float(np.linalg.norm(pos[i] - pos[j]))

    bonds = [[] for _ in range(3)]
    L01 = dist(0, 1)
    L02 = dist(0, 2)
    bonds[0].append((1, L01, float(k_oh))); bonds[1].append((0, L01, float(k_oh)))
    bonds[0].append((2, L02, float(k_oh))); bonds[2].append((0, L02, float(k_oh)))

    if add_angle:
        L12 = dist(1, 2)
        bonds[1].append((2, L12, float(k_hh))); bonds[2].append((1, L12, float(k_hh)))

    return pos, bonds


def bonds_to_max_L0(bonds_adj, default=1.0):
    Lmax = 0.0
    for bl in bonds_adj:
        for (j, L0, K) in bl:
            Lmax = max(Lmax, float(L0))
    return Lmax if Lmax > 0.0 else float(default)


def ensure_outdir(path):
    os.makedirs(path, exist_ok=True)
    return path


def print_run_header(tag, params):
    keys = sorted(params.keys())
    s = " ".join([f"{k}={params[k]}" for k in keys])
    print(f"[RUN] {tag} {s}")


def plot_residual_series(res, title="", noshow=False, outpath=None):
    import matplotlib.pyplot as plt

    res = np.asarray(res, dtype=np.float32)
    plt.figure(figsize=(6, 4))
    plt.semilogy(res[:, 0], np.maximum(res[:, 1], 1e-12), label="L_inf")
    if res.shape[1] > 2:
        plt.semilogy(res[:, 0], np.maximum(res[:, 2], 1e-12), label="L2")
    plt.xlabel("Jacobi iteration")
    plt.ylabel("bond residual")
    plt.title(title)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    if outpath is not None:
        plt.savefig(outpath, dpi=150)
    if not noshow:
        plt.show()
    plt.close()


def load_xyz(fname):
    with open(fname, 'r') as f:
        lines_raw = [l.rstrip('\n') for l in f.readlines()]
    if len(lines_raw) < 2:
        raise ValueError(f"load_xyz: file too short '{fname}'")
    n = int(lines_raw[0].strip().split()[0])
    if n <= 0:
        raise ValueError(f"load_xyz: invalid nAtoms={n} in '{fname}'")
    comment = lines_raw[1].strip()
    i0 = 2
    if comment.startswith('-lvs'):
        i0 = 2
    while i0 < len(lines_raw) and (lines_raw[i0].strip() == ''):
        i0 += 1
    if i0 + n > len(lines_raw):
        raise ValueError(f"load_xyz: not enough atom lines in '{fname}' (need {n}, have {len(lines_raw)-i0})")
    elems = []
    xyz = np.zeros((n, 3), dtype=np.float32)
    q = np.zeros((n,), dtype=np.float32)
    for i in range(n):
        ws = lines_raw[i0 + i].split()
        elems.append(ws[0])
        xyz[i, 0] = float(ws[1]); xyz[i, 1] = float(ws[2]); xyz[i, 2] = float(ws[3])
        if len(ws) > 4:
            q[i] = float(ws[4])
    return elems, xyz, q


def masses_from_elems(elems):
    m = np.zeros((len(elems),), dtype=np.float32)
    for i, e in enumerate(elems):
        if e == 'H':
            m[i] = 1.0
        elif e == 'C':
            m[i] = 12.0
        elif e == 'N':
            m[i] = 14.0
        elif e == 'O':
            m[i] = 16.0
        else:
            raise ValueError(f"Unknown element '{e}'")
    return m


def quat_mul(a, b):
    out = np.empty_like(a)
    ax, ay, az, aw = a[:, 0], a[:, 1], a[:, 2], a[:, 3]
    bx, by, bz, bw = b[:, 0], b[:, 1], b[:, 2], b[:, 3]
    out[:, 0] = aw*bx + ax*bw + ay*bz - az*by
    out[:, 1] = aw*by - ax*bz + ay*bw + az*bx
    out[:, 2] = aw*bz + ax*by - ay*bx + az*bw
    out[:, 3] = aw*bw - ax*bx - ay*by - az*bz
    return out


def normalize_quat(q):
    norms = np.linalg.norm(q, axis=1)
    norms[norms == 0.0] = 1.0
    return q / norms[:, None]


def perturb_state(pos, quat, pos_scale, rot_scale, rng):
    pos_out = pos.copy()
    quat_out = quat.copy()
    if pos_scale > 0:
        pos_out += rng.normal(scale=pos_scale, size=pos_out.shape)
    if rot_scale > 0:
        axes = rng.normal(size=(len(quat_out), 3))
        axis_norm = np.linalg.norm(axes, axis=1) + 1e-12
        axes /= axis_norm[:, None]
        angles = rng.normal(scale=rot_scale, size=(len(quat_out),))
        half = 0.5 * angles
        sin_half = np.sin(half)
        dq = np.zeros_like(quat_out)
        dq[:, :3] = axes * sin_half[:, None]
        dq[:, 3] = np.cos(half)
        quat_out = normalize_quat(quat_mul(dq, quat_out))
    return pos_out, quat_out


def write_xyz_with_ports(fname, elems, pos, pneigh, port_n):
    n = pos.shape[0]
    nports = int(np.sum(port_n))
    with open(fname, 'a') as f:
        f.write(f"{n + nports}\n")
        f.write("\n")
        for i in range(n):
            f.write(f"{elems[i]} {pos[i,0]:.6f} {pos[i,1]:.6f} {pos[i,2]:.6f}\n")
        ip = 0
        for i in range(n):
            pi = pos[i]
            for k in range(int(port_n[i])):
                tip = pi + pneigh[i, k, :3]
                f.write(f"X {tip[0]:.6f} {tip[1]:.6f} {tip[2]:.6f}\n")
                ip += 1


def write_pdb_trajectory(filename, frames, symbols, bonds):
    with open(filename, 'w') as f:
        for model_idx, coords in enumerate(frames):
            f.write(f"MODEL     {model_idx + 1:4}\n")
            for i, (sym, pos) in enumerate(zip(symbols, coords)):
                f.write(f"HETATM{i+1:5} {sym:4} UNK     1    {pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}  1.00  0.00          {sym:>2}\n")
            f.write("ENDMDL\n")
        for b1, b2 in bonds:
            f.write(f"CONECT{b1+1:5}{b2+1:5}\n")
        f.write("END\n")


class LivePortViz:
    """Lightweight live 3D updater that keeps camera/view persistent."""
    def __init__(self, elems):
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import proj3d
        self.plt = plt
        self.proj3d = proj3d
        self.elems = elems
        plt.ion()
        self.fig = plt.figure(figsize=(6, 6))
        self.ax = self.fig.add_subplot(111, projection='3d')
        self.sc = self.ax.scatter([], [], [], s=60, c='k')
        self.lines = []
        self.quiver = None
        self.labels = []
        for _ in elems:
            self.labels.append(self.ax.text(0, 0, 0, '', zorder=10))
        self.ax.set_xlabel('x'); self.ax.set_ylabel('y'); self.ax.set_zlabel('z')
        self.fig.canvas.draw()
        self.fig.show()
        self._last_pos = None
        self._last_proj = None

    def ensure_lines(self, total_ports):
        while len(self.lines) < total_ports:
            ln, = self.ax.plot([0, 0], [0, 0], [0, 0], '-', c='C0', lw=1)
            self.lines.append(ln)

    def update(self, pos, pneigh, port_n, force=None, title=""):
        self.ax.set_title(title)
        self._last_pos = np.asarray(pos, dtype=np.float32)
        self.sc._offsets3d = (pos[:, 0], pos[:, 1], pos[:, 2])
        for i, lab in enumerate(self.labels):
            lab.set_position((pos[i, 0], pos[i, 1]))
            lab.set_3d_properties(pos[i, 2], zdir='z')
            lab.set_text(self.elems[i])
        total_ports = int(np.sum(port_n))
        self.ensure_lines(total_ports)
        ip = 0
        for i in range(pos.shape[0]):
            pi = pos[i]
            for k in range(int(port_n[i])):
                tip = pi + pneigh[i, k, :3]
                self.lines[ip].set_data_3d([pi[0], tip[0]], [pi[1], tip[1]], [pi[2], tip[2]])
                ip += 1
        for j in range(ip, len(self.lines)):
            self.lines[j].set_data_3d([0, 0], [0, 0], [0, 0])
        if self.quiver is not None:
            self.quiver.remove()
            self.quiver = None
        if force is not None:
            f = force[:, :3]
            self.quiver = self.ax.quiver(pos[:, 0], pos[:, 1], pos[:, 2], f[:, 0], f[:, 1], f[:, 2], length=0.1, normalize=True, color='r')
        self.plt.pause(0.001)


def attach_picker_3d(viz, sim, *, pick_radius_px=20, verbose=0):
    if viz is None:
        raise ValueError('attach_picker_3d: viz is None')
    if sim is None:
        raise ValueError('attach_picker_3d: sim is None')
    if not hasattr(viz, '_last_pos'):
        raise ValueError('attach_picker_3d: viz has no _last_pos')

    pick = {"idx": None, "active": False, "mouse": np.array([0.0, 0.0], dtype=np.float32), "mouse3": np.array([0.0, 0.0, 0.0], dtype=np.float32), "_zproj": None}

    def _project_all(pos):
        ax = viz.ax
        M = ax.get_proj()
        x2, y2, z2 = viz.proj3d.proj_transform(pos[:, 0], pos[:, 1], pos[:, 2], M)
        pts2 = np.vstack([x2, y2]).T
        pts_px = ax.transData.transform(pts2)
        return pts_px, z2

    def on_press(event):
        if event.inaxes != viz.ax:
            return
        if viz._last_pos is None:
            return
        pos = np.asarray(viz._last_pos, dtype=np.float32)
        if pos.size == 0:
            return
        pts_px, z2 = _project_all(pos)
        mouse_px = np.array([event.x, event.y], dtype=np.float32)
        d2 = np.sum((pts_px - mouse_px[None, :]) ** 2, axis=1)
        i_min = int(np.argmin(d2))
        if d2[i_min] <= float(pick_radius_px) ** 2:
            pick["idx"] = i_min
            pick["active"] = True
            pick["mouse"] = mouse_px
            pick["_zproj"] = float(z2[i_min])
            pick["mouse3"] = pos[i_min].copy()
            if int(verbose) > 0:
                print(f"[DEBUG] pick3d press idx={i_min} dpx={np.sqrt(float(d2[i_min])):.2f}")

    def on_release(event):
        if pick["active"] and int(verbose) > 0:
            print(f"[DEBUG] pick3d release idx={pick['idx']}")
        pick["active"] = False
        pick["idx"] = None
        pick["_zproj"] = None

    def on_motion(event):
        if not pick["active"]:
            return
        if event.inaxes != viz.ax:
            return
        if pick.get("_zproj", None) is None:
            return
        # Convert mouse pixels -> projection coords -> data coords at fixed projected depth
        mouse_px = np.array([event.x, event.y], dtype=np.float32)
        pick["mouse"] = mouse_px
        ax = viz.ax
        xy_proj = ax.transData.inverted().transform(mouse_px)
        try:
            x, y, z = viz.proj3d.inv_transform(float(xy_proj[0]), float(xy_proj[1]), float(pick["_zproj"]), ax.get_proj())
        except Exception as e:
            raise RuntimeError(f"attach_picker_3d: proj3d.inv_transform failed: {e}")
        pick["mouse3"] = np.array([x, y, z], dtype=np.float32)
        ia = int(pick["idx"])
        sim.set_atom_pos(ia, pick["mouse3"])

    viz.fig.canvas.mpl_connect('button_press_event', on_press)
    viz.fig.canvas.mpl_connect('button_release_event', on_release)
    viz.fig.canvas.mpl_connect('motion_notify_event', on_motion)
    return pick


def plot_state_with_ports(elems, pos, pneigh, port_n, force=None, *, title=""):
    if not hasattr(plot_state_with_ports, '_viz'):
        plot_state_with_ports._viz = LivePortViz(elems)
    viz = plot_state_with_ports._viz
    viz.update(pos, pneigh, port_n, force=force, title=title)
