#!/usr/bin/env python3
import os, sys, argparse
import numpy as np

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_ROOT     = os.path.realpath(os.path.join(_THIS_DIR, '..', '..'))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from pyBall.OCL.AFM import AFMulator, compute_df


def _parse_int3(s):
    p = [int(x) for x in s.split(',')]
    if len(p) != 3:
        raise argparse.ArgumentTypeError(f"expected 3 ints as a,b,c got '{s}'")
    return tuple(p)


def _parse_float3(s):
    p = [float(x) for x in s.split(',')]
    if len(p) != 3:
        raise argparse.ArgumentTypeError(f"expected 3 floats as a,b,c got '{s}'")
    return tuple(p)


def _ensure_finite(name, arr, abs_max=1e6):
    if not np.isfinite(arr).all():
        bad = np.where(~np.isfinite(arr))
        raise RuntimeError(f"{name}: non-finite values present; first bad index={tuple(int(x[0]) for x in bad)}")
    m = float(np.max(np.abs(arr))) if arr.size else 0.0
    if m > abs_max:
        raise RuntimeError(f"{name}: abs_max={m:.3e} exceeds threshold {abs_max:.3e}")


def _bytes_to_gb(nbytes):
    return nbytes / (1024.0**3)


def _estimate_grid_bytes(nx, ny, nz, components=4, dtype_bytes=4):
    return int(nx) * int(ny) * int(nz) * int(components) * int(dtype_bytes)


def _estimate_scan_bytes(nx_s, ny_s, nz, components=4, dtype_bytes=4):
    return int(nx_s) * int(ny_s) * int(nz) * int(components) * int(dtype_bytes)


def getDfWeight(n, dz=0.1):
    x = np.linspace(-1, 1, n+1)
    y = np.sqrt(np.maximum(0.0, 1 - x*x))
    dy = (y[1:] - y[:-1]) / (dz*n)
    fpi = (n-2)**2
    prefactor = -1.0 * (1.0 + fpi*(2/np.pi)) / (fpi + 1.0)
    return dy*prefactor, (x[1:]+x[:-1])*0.5


def Fz2df(F, dz=0.1, k0=1800.0, f0=30000.0, n=4, units=16.0217656):
    """Giessibl df via convolution (ppafm-style discretization).

    Assumes:
    - F is in eV/Å
    - dz is in Å
    Returns df in Hz.
    """
    dz = float(dz)
    if dz <= 0:
        raise ValueError(f"Fz2df: dz must be >0, got {dz}")
    W, xs = getDfWeight(int(n), dz=dz)
    if F.shape[2] < len(W):
        raise ValueError(f"Fz2df: need nz>=n, got nz={F.shape[2]} n={len(W)}")
    dFconv = np.apply_along_axis(lambda m: np.convolve(m, W, mode='valid'), axis=2, arr=F)
    return (dFconv*float(units)*float(f0)/float(k0)).astype(np.float32), xs


def _valid_heights(heights, n):
    """Heights array matching convolution(mode='valid') output length."""
    n = int(n)
    nz = int(len(heights))
    nzv = nz - (n - 1)
    if nzv <= 0:
        raise ValueError(f"_valid_heights: nz={nz} too small for n={n}")
    i0 = n//2
    i1 = i0 + nzv
    return heights[i0:i1]


def _scan_extent_xy(scan_p0, scan_da, scan_db, nx, ny):
    dx = float(scan_da[0]); dy = float(scan_db[1])
    x0 = float(scan_p0[0]) - 0.5*dx
    x1 = float(scan_p0[0]) + dx*(nx - 0.5)
    y0 = float(scan_p0[1]) - 0.5*dy
    y1 = float(scan_p0[1]) + dy*(ny - 0.5)
    return (x0, x1), (y0, y1)


def _pick_cell_origin_xy(p_xy, a1, a2):
    M = np.array([[a1[0], a2[0]], [a1[1], a2[1]]], dtype=float)
    det = float(np.linalg.det(M))
    if abs(det) < 1e-8:
        return np.zeros(2)
    uv = np.linalg.solve(M, np.array([p_xy[0], p_xy[1]], dtype=float))
    iu = np.floor(uv[0]); iv = np.floor(uv[1])
    return a1*iu + a2*iv


def plot_slices_with_overlay(Fz, heights, label, fname, x_ext, y_ext, atoms_xy, lvec_xy, save_dir, cmap_mode='bwr', sym=True, atoms_top_xy=None, plot_iz=None, lvec_tiling=None, flip_x=False, flip_y=False):
    """Plot each z-slice as a separate PNG with lattice overlay and atom dots."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    os.makedirs(save_dir, exist_ok=True)
    nx, ny, nz = Fz.shape

    use_xy_grid = isinstance(x_ext, np.ndarray) and isinstance(y_ext, np.ndarray)

    ks = list(range(len(heights)))
    if plot_iz is not None:
        kiz = int(plot_iz)
        if kiz < 0: kiz = len(heights) + kiz
        if (kiz < 0) or (kiz >= len(heights)):
            raise ValueError(f"plot_iz={plot_iz} out of range for nz={len(heights)}")
        ks = [kiz]

    for k in ks:
        h = heights[k]
        fig, ax = plt.subplots(1, 1, figsize=(8, 6))
        data = Fz[:, :, k]
        if flip_x: data = data[::-1, :]
        if flip_y: data = data[:, ::-1]
        if sym:
            vabs = max(float(np.percentile(np.abs(data), 99)), 1e-6)
            norm = plt.Normalize(vmin=-vabs, vmax=vabs)
        else:
            norm = None
        if use_xy_grid:
            # X,Y are scan-point centers with shape (nx,ny); use shading='nearest' to accept same-shaped C
            im = ax.pcolormesh(x_ext, y_ext, data, cmap=cmap_mode, shading='nearest', norm=norm)
        else:
            im = ax.imshow(data, origin='lower', cmap=cmap_mode, aspect='equal',
                           extent=[x_ext[0], x_ext[1], y_ext[0], y_ext[1]], norm=norm)
        ax.set_xlabel('x (Å)')
        ax.set_ylabel('y (Å)')
        ax.set_title(f'{label} h={h:.2f}Å')
        plt.colorbar(im, ax=ax, shrink=0.5, label='Fz (eV/Å)')

        if lvec_xy is not None:
            o, a1, a2 = lvec_xy
            if lvec_tiling is None:
                corners = np.array([o, o+a1, o+a1+a2, o+a2, o])
                ax.plot(corners[:, 0], corners[:, 1], color='lime', lw=0.5, label='unit cell')
            else:
                ix0, ix1, iy0, iy1 = lvec_tiling
                for iyc in range(int(iy0), int(iy1)+1):
                    for ixc in range(int(ix0), int(ix1)+1):
                        oo = o + a1*ixc + a2*iyc
                        corners = np.array([oo, oo+a1, oo+a1+a2, oo+a2, oo])
                        ax.plot(corners[:, 0], corners[:, 1], color='lime', lw=0.35)

        if atoms_xy is not None and len(atoms_xy) > 0:
            ax.scatter(atoms_xy[:, 0], atoms_xy[:, 1], c='lime', s=0.5, alpha=1.0, linewidths=0, zorder=5)

        if atoms_top_xy is not None and len(atoms_top_xy) > 0:
            ax.scatter(atoms_top_xy[:, 0], atoms_top_xy[:, 1], facecolors='none', edgecolors='yellow', s=28.0, alpha=0.95, linewidths=0.9, zorder=6)

        ax.set_aspect('equal')
        if not use_xy_grid:
            ax.set_xlim(x_ext[0], x_ext[1])
            ax.set_ylim(y_ext[0], y_ext[1])
        plt.tight_layout()
        out_path = os.path.join(save_dir, f'{fname}_z{k:03d}.png')
        plt.savefig(out_path, dpi=120, bbox_inches='tight')
        plt.close()
        print(f"  Saved {out_path}")


def plot_Fz_curve(Fz, heights, label, fname, center_idx, save_dir):
    """Plot F(z) curve at a specific (x,y) pixel."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    ix, iy = center_idx
    fz_profile = Fz[ix, iy, :]

    fig, ax = plt.subplots(1, 1, figsize=(6, 4))
    ax.plot(heights, fz_profile, 'b-', lw=1.5, marker='o', markersize=3)
    ax.set_xlabel('Height above top atom (Å)')
    ax.set_ylabel('Fz (eV/Å)')
    ax.set_title(f'{label} at pixel ({ix},{iy})')
    ax.grid(True, alpha=0.3)
    ax.axhline(0, color='k', lw=0.5)
    plt.tight_layout()
    out_path = os.path.join(save_dir, fname)
    plt.savefig(out_path, dpi=120, bbox_inches='tight')
    plt.close()
    print(f"  Saved {out_path}")


'''

### RUN LIKE THIS:

python3 afm_morse_pbc.py   --spacing 0.3,0.3,0.3   --scan_step 0.3,0.3,0.0   --z_start 6.0 --z_end 2.0 --dz 0.1   --df_n 4 --plot_sym --plot_cmap bwr --plot_exp

'''

def main():
    ap = argparse.ArgumentParser(description="AFM imaging (Morse + point charges) using relax.cl: build primitive-cell 3D FF image + PP relaxation (FIRE). Uses halo atom replication for FF projection and periodic image sampling for scan." )
    #ap.add_argument('--mol',     type=str,  default=os.path.join(_ROOT, 'tests', 'tMMFF', 'trj_ditetraceno_conf1_hamaker_final.mol2'))
    ap.add_argument('--mol',     type=str,  default=os.path.join(_ROOT, 'tests', 'tMMFF', 'trj_ditetraceno_conf2_hamaker_final.mol2'))
    ap.add_argument('--params',  type=str,  default=os.path.join(_ROOT, 'cpp', 'common_resources', 'ElementTypes.dat'))
    ap.add_argument('--outdir',  type=str,  default=os.path.join(_THIS_DIR, 'out_afm_morse_pbc'))
    ap.add_argument('--pbc',     type=_parse_int3, default=(1, 1, 0),       help='PBC replication radius as rx,ry,rz. Example 1,1,0 -> (3,3,1) copies (central).')
    ap.add_argument('--spacing', type=_parse_float3, default=(0.2, 0.2, 0.2), help='Voxel spacing (dx,dy,dz) in Angstrom. If set, overrides --grid_n.')
    ap.add_argument('--grid_n',  type=_parse_int3, default=(0, 0, 0),       help='Force-field grid resolution nx,ny,nz (used if --spacing not given).')
    ap.add_argument('--nxy',     type=_parse_int3, default=(120, 120, 0),   help='Scan resolution nx,ny,unused.')
    ap.add_argument('--scan_step',type=_parse_float3, default=(0.0, 0.0, 0.0), help='Scan step (dx,dy,dz) in Angstrom. If 0,0,0 uses same as spacing.')
    ap.add_argument('--margin',  type=float, default=2.0,  help='XY padding for auto grid box [A].')
    ap.add_argument('--margin_z',type=float, default=2.0, help='Bottom z padding for primitive-cell grid [A].')
    ap.add_argument('--z_top',   type=float, default=14.0, help='Extra z above molecule for grid box [A].')
    ap.add_argument('--z_start', type=float, default=6.0,  help='Probe start height above top atom [A].')
    ap.add_argument('--z_end',   type=float, default=2.0,  help='Probe end height above top atom [A].')
    ap.add_argument('--dz',      type=float, default=0.10, help='Probe z step [A].')
    ap.add_argument('--df_n',    type=int,   default=4,    help='Giessibl convolution n (A = n*dz).')
    ap.add_argument('--k0',      type=float, default=1800.0, help='Cantilever stiffness [N/m] for df conversion.')
    ap.add_argument('--f0',      type=float, default=30000.0, help='Cantilever frequency [Hz] for df conversion.')
    ap.add_argument('--df_units', type=float, default=0.0, help='Deprecated (kept for CLI compat). Not used when force is in eV/Å.')
    ap.add_argument('--plot_cmap', type=str, default='bwr', choices=['bwr', 'gray', 'magma'], help='Colormap for per-slice plots.')
    ap.add_argument('--plot_sym', action='store_true', help='Use symmetric +/- normalization (useful for debugging).')
    ap.add_argument('--plot_exp', action='store_true', help='Also export grayscale (experimental-style) slices for df.')
    ap.add_argument('--exp_only', action='store_true', help='If set, output only experimental (gray) slices and skip standard colormap plots.')
    ap.add_argument('--tag',        type=str,   default='',   help='Extra tag appended to output filenames.')
    ap.add_argument('--scan_frac',  type=_parse_float3, default=(0.0, 1.0, 0.0),   help='Scan in x,y uses [f0,f1] fraction of bbox span.')
    ap.add_argument('--scan_cells', type=_parse_float3, default=(3.0, 3.0, 0.0), help='Scan size in unit cells along a,b. Example 3,3,0 scans 3x3 cells (periodicity check).')
    ap.add_argument('--scan_cart',  type=int,   default=1, help='1=use rectangular scan grid in lab XY (imshow), 0=use lattice-aligned scan (skew pcolormesh).')
    ap.add_argument('--plot_iz',    type=int,   default=-1, help='Which z-slice index to plot (-1=all, -2=second last, etc.). Useful for fast visual debugging.')
    ap.add_argument('--top_z_tol',  type=float, default=0.5, help='Top atom selection tolerance [A] for overlay circles (within tol of max z).')
    ap.add_argument('--phase_align',type=int,   default=0, help='If 1, auto-shift scan_p0 in XY so strongest feature in selected slice lands on a lattice corner (integer u,v). Reruns scan and plots with tag suffix _ph.')
    ap.add_argument('--flip_x',     type=int,   default=0, help='If 1, flip AFM images along x before plotting (visual debug).')
    ap.add_argument('--flip_y',     type=int,   default=1, help='If 1, flip AFM images along y before plotting (visual debug).')
    ap.add_argument('--shift_xy',   type=float, nargs=3, default=(0.0, 0.0, 0.0), metavar=('DX','DY','DZ'), help='Manual XY shift [Ang] applied to plotting frame (z ignored). Pass three floats separated by spaces.')
    ap.add_argument('--raw', action='store_true', help='Also compute raw (no PP relax) scan.')
    ap.add_argument('--noplot', action='store_true', help='Do not generate png plots.')

    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    rx, ry, rz = args.pbc
    nx_s, ny_s, _ = args.nxy
    spacing = args.spacing
    use_spacing = spacing[0] > 0 or spacing[1] > 0 or spacing[2] > 0

    dz = float(args.dz)
    if dz <= 0:
        raise ValueError('--dz must be >0')
    dtip = -dz
    nz = int(round((args.z_start - args.z_end) / dz)) + 1
    if nz < 2:
        raise ValueError(f"Invalid z range: z_start={args.z_start} z_end={args.z_end} dz={dz} -> nz={nz}")

    tag = 'morse'
    tag += f"_pbc{rx}{ry}{rz}"
    if args.tag:
        tag += f"_{args.tag}"

    print(f"=== AFM Morse PBC ===")
    print(f"mol={args.mol}")
    print(f"params={args.params}")
    print(f"outdir={args.outdir}")
    print(f"pbc_radius={args.pbc}")
    if use_spacing:
        print(f"spacing={spacing}")
    else:
        print(f"grid_n={args.grid_n}")
    print(f"scan_nxy=({nx_s},{ny_s})  scan_step={args.scan_step}")
    print(f"z: start={args.z_start} end={args.z_end} dz={dz} -> nz={nz}")

    afm = AFMulator(use_morse=True)
    mol_prim = afm.load_molecule(args.mol)
    mol = mol_prim

    if (rx, ry, rz) != (0, 0, 0):
        if getattr(mol_prim, 'lvec', None) is None:
            raise RuntimeError('Requested --pbc != 0,0,0 but molecule has no lattice vectors (mol.lvec is None)')
        mol = mol_prim.clonePBC_images_central(nPBC=(rx, ry, rz))
        afm.mol = mol
        print(f"PBC halo atoms: natoms {mol.natoms}  lvec(primitive)=\n{mol.lvec}")

    if not os.path.exists(args.params):
        raise FileNotFoundError(args.params)
    afm.assign_params(params_path=args.params)

    # NOTE: atoms_arr will be shifted by setup_grid/setup_grid_lvec(); for scan box use primitive cell vectors
    if getattr(mol_prim, 'lvec', None) is None:
        raise RuntimeError('This script expects mol.lvec for primitive-cell grid + periodic scan')
    a = mol_prim.lvec[0, :3].astype(float)
    b = mol_prim.lvec[1, :3].astype(float)
    a1 = mol_prim.lvec[0, :2].astype(float)
    a2 = mol_prim.lvec[1, :2].astype(float)
    la = float(np.linalg.norm(a))
    lb = float(np.linalg.norm(b))
    if la <= 0 or lb <= 0:
        raise RuntimeError(f"Invalid lattice vectors norms |a|={la} |b|={lb}")

    if use_spacing:
        sx, sy, sz = spacing
        if sx <= 0: sx = 0.1
        if sy <= 0: sy = 0.1
        if sz <= 0: sz = 0.1
        nx_g = max(8, int(np.ceil(la / sx)))
        ny_g = max(8, int(np.ceil(lb / sy)))
        apos_z = mol.apos[:, 2]
        span_z = float(apos_z.max() - apos_z.min())
        nz_g = max(8, int(np.ceil((span_z + float(args.margin_z) + float(args.z_top)) / sz)))
        print(f"Computed primitive grid_n from spacing: ({nx_g},{ny_g},{nz_g}) |a|={la:.3f} |b|={lb:.3f}")
    else:
        nx_g, ny_g, nz_g = args.grid_n
        if nx_g <= 0 or ny_g <= 0 or nz_g <= 0:
            raise ValueError('For primitive-cell grid you must specify either --spacing or valid --grid_n')

    vram_bytes = _estimate_grid_bytes(nx_g, ny_g, nz_g, components=4, dtype_bytes=4)
    print(f"[mem] FF grid {nx_g}x{ny_g}x{nz_g} float4 -> {vram_bytes/1e9:.3f} GB")

    afm.setup_grid_lvec(n=(nx_g, ny_g, nz_g), margin_z=float(args.margin_z), z_top=float(args.z_top))
    afm.make_forcefield()

    sca, scb, _ = args.scan_cells
    sca = float(sca); scb = float(scb)
    if sca <= 0 or scb <= 0:
        raise ValueError(f"scan_cells must be >0 in a,b, got {args.scan_cells}")

    scan_step = args.scan_step
    if args.scan_cart:
        # rectangular scan in lab XY
        Lx = float(a[0]) * sca
        Ly = float(b[1]) * scb
        if (Lx <= 0) or (Ly <= 0):
            raise RuntimeError(f"Invalid cart scan size from lattice projections Lx={Lx} Ly={Ly} (a={a} b={b})")
        if scan_step[0] > 0 and scan_step[1] > 0:
            dx = float(scan_step[0]); dy = float(scan_step[1])
            nx_s = int(np.ceil(Lx / dx)) + 1
            ny_s = int(np.ceil(Ly / dy)) + 1
        dx = Lx / max(nx_s-1, 1)
        dy = Ly / max(ny_s-1, 1)
        scan_da = np.array([dx, 0.0, 0.0], dtype=np.float32)
        scan_db = np.array([0.0, dy, 0.0], dtype=np.float32)
        print(f"Rectangular scan grid: ({nx_s},{ny_s}) dx={dx:.3f} dy={dy:.3f}")
    else:
        # lattice-aligned scan (skew)
        if scan_step[0] > 0 and scan_step[1] > 0:
            da_len = float(scan_step[0])
            db_len = float(scan_step[1])
            nx_s = int(np.ceil((sca*la) / da_len)) + 1
            ny_s = int(np.ceil((scb*lb) / db_len)) + 1
            scan_da = (a * (da_len / la)).astype(np.float32)
            scan_db = (b * (db_len / lb)).astype(np.float32)
            print(f"Computed scan nxy from scan_step along a,b: ({nx_s},{ny_s}) da_len={da_len:.3f} db_len={db_len:.3f}")
        else:
            scan_da = (a * (sca / max(nx_s-1, 1))).astype(np.float32)
            scan_db = (b * (scb / max(ny_s-1, 1))).astype(np.float32)

    mol_z = afm.mol_z
    z0_tip = mol_z + float(args.z_start) + abs(float(afm.dpos0[2]))

    shift_xy = np.array(args.shift_xy[:2], dtype=np.float32)  # reuse later (plotting frame only)
    scan_p0 = np.array([0.0, 0.0, z0_tip], dtype=np.float32)

    probe_heights = np.array([float(args.z_start) + k*dtip for k in range(nz)], dtype=np.float32)

    print(f"scan_p0={scan_p0}  da={scan_da}  db={scan_db}")
    print(f"mol_z={mol_z:.3f}  z0_tip={z0_tip:.3f}")
    print(f"probe_heights: {probe_heights[0]:.2f} .. {probe_heights[-1]:.2f} A (above top atom)")

    np.save(os.path.join(args.outdir, f'probe_heights_{tag}.npy'), probe_heights)

    scan_bytes = _estimate_scan_bytes(nx_s, ny_s, nz, components=4, dtype_bytes=4)
    print(f"[mem] Scan buffers {nx_s}x{ny_s}x{nz} float4 -> {scan_bytes/1e9:.3f} GB (per FEs array)")

    if args.raw:
        print("--- RAW scan (no PP relax) ---")
        FEs_raw, _ = afm.get_raw_FE(nxy=(nx_s, ny_s), nz=nz, dtip=dtip, scan_p0=scan_p0, scan_da=scan_da, scan_db=scan_db)
        Fz_raw = FEs_raw[:, :, :, 2]
        _ensure_finite('Fz_raw', Fz_raw)
        print(f"Fz_raw stats: min={Fz_raw.min():.4g} max={Fz_raw.max():.4g} mean={Fz_raw.mean():.4g}")
        np.save(os.path.join(args.outdir, f'FEs_raw_{tag}.npy'), FEs_raw)
        np.save(os.path.join(args.outdir, f'Fz_raw_{tag}.npy'),  Fz_raw)

    def _do_relaxed(scan_p0_):
        FEs_relax_, pts_ = afm.run_scan(nxy=(nx_s, ny_s), nz=nz, dtip=dtip,
                                      scan_p0=scan_p0_, scan_da=scan_da, scan_db=scan_db)
        Fz_relax_ = FEs_relax_[:, :, :, 2]
        _ensure_finite('Fz_relax', Fz_relax_, abs_max=1e6)
        return FEs_relax_, Fz_relax_, pts_

    print("--- RELAXED scan (PP FIRE) ---")
    FEs_relax, Fz_relax, pts = _do_relaxed(scan_p0)
    print(f"Fz_relax stats: min={Fz_relax.min():.6g} max={Fz_relax.max():.6g} mean={Fz_relax.mean():.6g}")

    if int(args.phase_align) != 0:
        plot_iz0 = int(args.plot_iz)
        if plot_iz0 == -1: plot_iz0 = nz - 2
        if plot_iz0 < 0: plot_iz0 = nz + plot_iz0
        if (plot_iz0 < 0) or (plot_iz0 >= nz):
            raise ValueError(f"phase_align: plot_iz={args.plot_iz} out of range for nz={nz}")

        sl = Fz_relax[:, :, plot_iz0]
        # ignore borders to avoid wrap artifacts
        bx = max(3, nx_s//50)
        by = max(3, ny_s//50)
        sub = sl[bx:nx_s-bx, by:ny_s-by]
        mi = np.unravel_index(np.argmax(sub), sub.shape)
        ixm, iym = int(mi[0] + bx), int(mi[1] + by)
        pxy = pts[ixm, iym, :2].astype(float)

        # map XY -> fractional u,v in lattice basis
        M = np.array([[a1[0], a2[0]], [a1[1], a2[1]]], dtype=float)
        det = float(np.linalg.det(M))
        if abs(det) < 1e-12:
            raise RuntimeError(f"phase_align: singular lattice det={det}")
        uv = np.linalg.solve(M, pxy)
        fu = float(uv[0] - np.floor(uv[0]))
        fv = float(uv[1] - np.floor(uv[1]))
        shift = a1*fu + a2*fv
        print(f"phase_align: slice k={plot_iz0} peak@({ixm},{iym}) pos={pxy} uv={uv} frac=({fu:.6f},{fv:.6f}) shift_xy={shift}")

        scan_p0_ph = scan_p0.copy()
        scan_p0_ph[0] -= float(shift[0])
        scan_p0_ph[1] -= float(shift[1])

        tag = tag + "_ph"
        print("--- RELAXED scan (phase-aligned rerun) ---")
        FEs_relax, Fz_relax, pts = _do_relaxed(scan_p0_ph)
        print(f"Fz_relax(ph) stats: min={Fz_relax.min():.6g} max={Fz_relax.max():.6g} mean={Fz_relax.mean():.6g}")
        scan_p0 = scan_p0_ph

    # Giessibl convolution df (smooth)
    df_relax, xs = Fz2df(Fz_relax, dz=abs(dz), k0=args.k0, f0=args.f0, n=args.df_n)
    probe_heights_df = _valid_heights(probe_heights, args.df_n)
    _ensure_finite('df_relax', df_relax, abs_max=1e8)
    print(f"df_relax stats: min={df_relax.min():.4g} max={df_relax.max():.4g} mean={df_relax.mean():.4g}")

    np.save(os.path.join(args.outdir, f'FEs_relax_{tag}.npy'), FEs_relax)
    np.save(os.path.join(args.outdir, f'Fz_relax_{tag}.npy'),  Fz_relax)
    np.save(os.path.join(args.outdir, f'df_relax_{tag}.npy'),  df_relax)
    np.save(os.path.join(args.outdir, f'scan_pts_{tag}.npy'),   pts)

    if not args.noplot:
        ix = np.arange(nx_s, dtype=np.float32)[:, None]
        iy = np.arange(ny_s, dtype=np.float32)[None, :]
        X = scan_p0[0] + ix*scan_da[0] + iy*scan_db[0]
        Y = scan_p0[1] + ix*scan_da[1] + iy*scan_db[1]

        apos_shifted = mol_prim.apos[:, :3] + afm.mol_shift[None, :]
        atoms0_xy = apos_shifted[:, :2]
        zmax = float(apos_shifted[:, 2].max())
        top_sel = np.where(apos_shifted[:, 2] > (zmax - float(args.top_z_tol)))[0]
        atoms0_top_xy = atoms0_xy[top_sel, :] if len(top_sel) else None
        a1 = mol_prim.lvec[0, :2].astype(float)
        a2 = mol_prim.lvec[1, :2].astype(float)
        o  = np.zeros(2)

        if args.scan_cart:
            x_ext = [float(X.min()) + shift_xy[0], float(X.max()) + shift_xy[0]]
            y_ext = [float(Y.min()) + shift_xy[1], float(Y.max()) + shift_xy[1]]
            X_plot, Y_plot = x_ext, y_ext
        else:
            X_plot, Y_plot = X + shift_xy[0], Y + shift_xy[1]

        o_shift = o + shift_xy
        lvec_xy = (o_shift, a1, a2)

        # Replicate atoms across the plotted region (so you can check alignment everywhere, not only in one cell)
        if args.scan_cart:
            M = np.array([[a1[0], a2[0]], [a1[1], a2[1]]], dtype=float)
            det = float(np.linalg.det(M))
            if abs(det) < 1e-12:
                raise RuntimeError(f"Singular lattice in XY det={det}")
            invM = np.linalg.inv(M)
            corners = np.array([
                [x_ext[0], y_ext[0]],
                [x_ext[1], y_ext[0]],
                [x_ext[0], y_ext[1]],
                [x_ext[1], y_ext[1]],
            ], dtype=float)
            uv = (invM @ corners.T).T
            umin, umax = float(np.min(uv[:, 0])), float(np.max(uv[:, 0]))
            vmin, vmax = float(np.min(uv[:, 1])), float(np.max(uv[:, 1]))
            iu0, iu1 = int(np.floor(umin)) - 1, int(np.ceil(umax)) + 1
            iv0, iv1 = int(np.floor(vmin)) - 1, int(np.ceil(vmax)) + 1
            atoms_xy = []
            atoms_top_xy = []
            for iv in range(iv0, iv1+1):
                for iu in range(iu0, iu1+1):
                    sh = a1*iu + a2*iv
                    atoms_xy.append(atoms0_xy + sh[None, :])
                    if atoms0_top_xy is not None:
                        atoms_top_xy.append(atoms0_top_xy + sh[None, :])
            atoms_xy = np.concatenate(atoms_xy, axis=0) if len(atoms_xy) else atoms0_xy
            atoms_top_xy = np.concatenate(atoms_top_xy, axis=0) if len(atoms_top_xy) else atoms0_top_xy
            lvec_tiling = (iu0, iu1, iv0, iv1)
        else:
            atoms_xy = atoms0_xy + shift_xy[None, :]
            atoms_top_xy = atoms0_top_xy + shift_xy[None, :] if atoms0_top_xy is not None else None
            lvec_tiling = None

        plot_iz = None if int(args.plot_iz) == -1 else int(args.plot_iz)
        if plot_iz is not None:
            # df_relax has shorter z dimension due to convolution(mode='valid')
            df_shift = int(args.df_n)//2
            plot_iz_df = int(plot_iz) - df_shift
            if plot_iz_df < 0: plot_iz_df = 0
            if plot_iz_df >= len(probe_heights_df): plot_iz_df = len(probe_heights_df)-1
        else:
            plot_iz_df = None

        slices_dir = os.path.join(args.outdir, f'slices_{tag}')
        os.makedirs(slices_dir, exist_ok=True)
        print(f"Plotting z-slices to {slices_dir} ...")
        if not args.exp_only:
            plot_slices_with_overlay(-Fz_relax, probe_heights, f'Fz relax [{tag}]', f'Fz_relax_{tag}', X_plot, Y_plot, atoms_xy, lvec_xy, slices_dir, cmap_mode=args.plot_cmap, sym=args.plot_sym or (args.plot_cmap == 'bwr'), atoms_top_xy=atoms_top_xy, plot_iz=plot_iz, lvec_tiling=lvec_tiling, flip_x=bool(args.flip_x), flip_y=bool(args.flip_y))
            plot_slices_with_overlay(-df_relax, probe_heights_df, f'df relax [{tag}]', f'df_relax_{tag}', X_plot, Y_plot, atoms_xy, lvec_xy, slices_dir, cmap_mode=args.plot_cmap, sym=args.plot_sym or (args.plot_cmap == 'bwr'), atoms_top_xy=atoms_top_xy, plot_iz=plot_iz_df, lvec_tiling=lvec_tiling, flip_x=bool(args.flip_x), flip_y=bool(args.flip_y))
        if args.plot_exp:
            plot_slices_with_overlay(-df_relax, probe_heights_df, f'df relax exp [{tag}]', f'df_relax_exp_{tag}', X_plot, Y_plot, atoms_xy, lvec_xy, slices_dir, cmap_mode='gray', sym=False, atoms_top_xy=atoms_top_xy, plot_iz=plot_iz_df, lvec_tiling=lvec_tiling, flip_x=bool(args.flip_x), flip_y=bool(args.flip_y))

        center_idx = (nx_s // 2, ny_s // 2)
        plot_Fz_curve(Fz_relax, probe_heights, 'Fz relax', f'Fz_curve_{tag}.png', center_idx, args.outdir)
        plot_Fz_curve(df_relax, probe_heights_df, 'df relax', f'df_curve_{tag}.png', center_idx, args.outdir)

        if args.raw:
            raw_path = os.path.join(args.outdir, f'Fz_raw_{tag}.npy')
            if not os.path.exists(raw_path):
                print(f"Raw file not found for current tag (likely due to phase_align): {raw_path}")
                print("Recomputing RAW scan for current scan_p0/tag ...")
                FEs_raw, _ = afm.get_raw_FE(nxy=(nx_s, ny_s), nz=nz, dtip=dtip, scan_p0=scan_p0, scan_da=scan_da, scan_db=scan_db)
                Fz_raw = FEs_raw[:, :, :, 2]
                _ensure_finite('Fz_raw', Fz_raw)
                np.save(os.path.join(args.outdir, f'FEs_raw_{tag}.npy'), FEs_raw)
                np.save(raw_path, Fz_raw)
            else:
                Fz_raw = np.load(raw_path)
            df_raw, xs = Fz2df(Fz_raw, dz=abs(dz), k0=args.k0, f0=args.f0, n=args.df_n)
            _ensure_finite('df_raw', df_raw, abs_max=1e8)
            np.save(os.path.join(args.outdir, f'df_raw_{tag}.npy'), df_raw)
            if not args.exp_only:
                plot_slices_with_overlay(-Fz_raw, probe_heights, f'Fz raw [{tag}]', f'Fz_raw_{tag}', X_plot, Y_plot, atoms_xy, lvec_xy, slices_dir, cmap_mode=args.plot_cmap, sym=args.plot_sym or (args.plot_cmap == 'bwr'), atoms_top_xy=atoms_top_xy, plot_iz=plot_iz, lvec_tiling=lvec_tiling, flip_x=bool(args.flip_x), flip_y=bool(args.flip_y))
                plot_slices_with_overlay(-df_raw, probe_heights_df, f'df raw [{tag}]', f'df_raw_{tag}', X_plot, Y_plot, atoms_xy, lvec_xy, slices_dir, cmap_mode=args.plot_cmap, sym=args.plot_sym or (args.plot_cmap == 'bwr'), atoms_top_xy=atoms_top_xy, plot_iz=plot_iz_df, lvec_tiling=lvec_tiling, flip_x=bool(args.flip_x), flip_y=bool(args.flip_y))
            if args.plot_exp:
                plot_slices_with_overlay(-df_raw, probe_heights_df, f'df raw exp [{tag}]', f'df_raw_exp_{tag}', X_plot, Y_plot, atoms_xy, lvec_xy, slices_dir, cmap_mode='gray', sym=False, atoms_top_xy=atoms_top_xy, plot_iz=plot_iz_df, lvec_tiling=lvec_tiling, flip_x=bool(args.flip_x), flip_y=bool(args.flip_y))
            plot_Fz_curve(Fz_raw, probe_heights, 'Fz raw', f'Fz_curve_raw_{tag}.png', center_idx, args.outdir)
            plot_Fz_curve(df_raw, probe_heights_df, 'df raw', f'df_curve_raw_{tag}.png', center_idx, args.outdir)

    print("=== DONE ===")


if __name__ == '__main__':
    main()
