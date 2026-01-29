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
