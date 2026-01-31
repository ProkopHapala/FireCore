import os

import numpy as np


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
        self.plt = plt
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

    def ensure_lines(self, total_ports):
        while len(self.lines) < total_ports:
            ln, = self.ax.plot([0, 0], [0, 0], [0, 0], '-', c='C0', lw=1)
            self.lines.append(ln)

    def update(self, pos, pneigh, port_n, force=None, title=""):
        self.ax.set_title(title)
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


def plot_state_with_ports(elems, pos, pneigh, port_n, force=None, *, title=""):
    if not hasattr(plot_state_with_ports, '_viz'):
        plot_state_with_ports._viz = LivePortViz(elems)
    viz = plot_state_with_ports._viz
    viz.update(pos, pneigh, port_n, force=force, title=title)
