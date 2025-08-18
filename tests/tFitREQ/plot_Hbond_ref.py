import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys
import argparse
import fnmatch


'''

Directory structure like this:

/home/prokop/Desktop/CARBSIS/PEOPLE/Paolo/HbondFit_small_mols_2025_08_15/

.
├── confs
│   ├── b3lyp
│   └── wb97m
├── data
│   ├── dftbp
│   ├── gfn0-xtb
│   ├── gfn1-xtb
│   ├── gfn2-xtb
│   ├── gfn-ff
│   ├── pm6
│   ├── pm6-d3
│   ├── pm6-d3h4
│   ├── pm6-dh2
│   ├── pm6-dhp
│   ├── pm6-ml
│   └── pm7
├── maps
└── ref
    ├── b3lyp
    └── wb97m


'''

# METHODS: reorder/comment out to select what to plot
METHODS = [
    'DFT-b3lyp', 'DFT-wb97m',
    'dftbp',
    'gfn0-xtb', 'gfn1-xtb', 'gfn2-xtb', 'gfn-ff',
    'pm6', 'pm6-d3', 'pm6-d3h4', 'pm6-dh2', 'pm6-dhp', 'pm6-ml', 'pm7',
]

# map folder names to display names
LABELS = {'b3lyp': 'DFT-b3lyp', 'wb97m': 'DFT-wb97m'}
def label(name):
    return LABELS.get(name, name)

def pick_file(files, pattern=None, index=None):
    names = [Path(p).name for p in files]
    if pattern:
        if any(c in pattern for c in '*?[]'):
            for p, n in zip(files, names):
                if fnmatch.fnmatch(n, pattern):
                    return p
        else:
            for p, n in zip(files, names):
                if pattern in n:
                    return p
    if index is not None and 0 <= index < len(files):
        return files[index]
    return files[0] if files else None


# -----------------------------
# Metrics and helper utilities
# -----------------------------
def _ref_shift(grid):
    if grid.size == 0:
        return grid, np.nan
    try:
        last = grid[-1, :] if grid.ndim == 2 and grid.shape[0] > 0 else np.array([np.nan])
        ref = np.nanmin(last[np.isfinite(last)]) if np.any(np.isfinite(last)) else 0.0
    except Exception:
        ref = 0.0
    gs = grid - ref
    mloc = float(np.nanmin(gs)) if np.any(np.isfinite(gs)) else np.nan
    if np.isfinite(mloc) and mloc > 0: mloc = 0.0
    return gs, mloc

def load_shifted(items):
    out = []  # (name, a, d, gs, mloc)
    mins = []
    for name, filepath in items:
        a, d, g = load_energy_map_gnuplot(filepath)
        gs, mloc = _ref_shift(g)
        out.append((name, a, d, gs, mloc))
        mins.append(mloc)
    return out, mins

def softclamp_rmse(model, ref, y1=0.0, y2=1.0):
    if model.shape != ref.shape:
        return np.nan
    diff = model - ref
    mask = np.isfinite(diff)
    if not np.any(mask):
        return np.nan
    e = diff.copy()
    m = mask & (e > y1)
    if np.any(m):
        z = (e[m] - y1) / (y2 - y1) if (y2 != y1) else np.zeros_like(e[m])
        e[m] = y1 + (y2 - y1) * (1.0 - 1.0 / (1.0 + z))
    return float(np.sqrt(np.nanmean((e[mask])**2)))

def rmse(a, b):
    if a.shape != b.shape:
        return np.nan
    m = np.isfinite(a) & np.isfinite(b)
    if not np.any(m): return np.nan
    return float(np.sqrt(np.mean((a[m]-b[m])**2)))

def extract_min_curves(a, d, g, rmax=None):
    if a.size == 0 or d.size == 0 or g.size == 0:
        return np.full_like(a, np.nan, dtype=float), np.full_like(a, np.nan, dtype=float)
    nA = a.shape[0]
    rmin = np.full(nA, np.nan)
    emin = np.full(nA, np.nan)
    for j in range(nA):
        col = g[:, j]
        if np.any(np.isfinite(col)):
            idx = int(np.nanargmin(col))
            emin[j] = col[idx]
            rmin[j] = d[idx]
            if (rmax is not None) and np.isfinite(rmin[j]) and (rmin[j] > rmax):
                rmin[j] = np.nan; emin[j] = np.nan
    return rmin, emin

def plot_min_lines(loaded, titles=None, panel=(3,5), show=True, savepath=None, suptitle=None, rmax=None):
    # loaded: list of (name, a, d, g, mloc)
    fig, (axR, axE) = plt.subplots(1, 2, figsize=(panel[0]*2, panel[1]), constrained_layout=True)
    curves = {}
    for i, (name, a, d, g, mloc) in enumerate(loaded):
        rmin, emin = extract_min_curves(a, d, g, rmax=rmax)
        curves[name] = (a, rmin, emin)
        lab = titles[i] if titles else name
        axR.plot(a, rmin, label=lab)
        axE.plot(a, emin, label=lab)
    axR.set_title('r_min(angle)')
    axE.set_title('E_min(angle)')
    for ax in (axR, axE):
        ax.set_xlabel('Angle (deg)')
        ax.grid(True, linestyle='--', alpha=0.4)
    axR.set_ylabel('Distance (Å)')
    axE.set_ylabel('Energy')
    axR.legend(ncol=1, fontsize=8)
    axE.legend(ncol=1, fontsize=8)
    if suptitle: fig.suptitle(suptitle)
    if savepath:
        Path(savepath).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(savepath, dpi=150, bbox_inches='tight')
    if show: plt.show()
    else: plt.close(fig)
    return curves

def align_and_rmse(a_ref, y_ref, a, y):
    a_ref = np.asarray(a_ref); y_ref = np.asarray(y_ref)
    a = np.asarray(a); y = np.asarray(y)
    mref = np.isfinite(y_ref)
    m = np.isfinite(y)
    if not np.any(mref) or not np.any(m):
        return np.nan
    if a.shape == a_ref.shape and np.allclose(a, a_ref):
        return rmse(y_ref, y)
    try:
        yi = np.interp(a_ref[mref], a[m], y[m])
    except Exception:
        return np.nan
    return float(np.sqrt(np.mean((yi - y_ref[mref])**2)))


def load_energy_map_gnuplot(filename):
    """
    Load hydrogen bond energy data from file.
    Supports gnuplot block format or flat triplets.
    Returns: angles, distances, 2D grid with shape (len(distances), len(angles)).
    Robust to comments, extra tokens, and NaNs.
    """
    rows = []
    with open(filename, "r") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s[0] in "#;%@":
                continue
            parts = s.split()
            if len(parts) < 3:
                continue
            try:
                a = float(parts[0]); d = float(parts[1]); e = float(parts[2])
            except ValueError:
                continue
            rows.append((a, d, e))

    if not rows:
        return np.array([]), np.array([]), np.full((0, 0), np.nan)

    angles = np.unique([r[0] for r in rows])
    distances = np.unique([r[1] for r in rows])
    grid = np.full((len(distances), len(angles)), np.nan)

    ai = {a:i for i,a in enumerate(angles)}
    di = {d:i for i,d in enumerate(distances)}
    for a, d, e in rows:
        grid[di[d], ai[a]] = e

    return angles, distances, grid

def plot_energy_map(angles, distances, energy_grid, title="Energy Map"):
    """
    Plot hydrogen bond energy as a colormap.
    NaNs will be left blank.
    """
    fig, ax = plt.subplots(figsize=(5,8))
    mesh = ax.pcolormesh(angles, distances, energy_grid, 
                         cmap="seismic", shading="auto")
    ax.set_xlabel("Angle (deg)")
    ax.set_ylabel("Distance (Å)")
    ax.set_title(title)
    plt.colorbar(mesh, ax=ax, label="Energy")
    plt.show()

def plot_all_methods(base_dirs, titles=None, panel=None, cmap="seismic", limit=None, share_scale=True, non_uniform=False, show=True, savepath=None, suptitle=None):
    """
    Plot multiple methods as subplots.
    base_dirs: dict {name: path_to_datfile}
    titles: list of titles in order (if None, use keys)
    """
    # keep insertion order (matches ORDER or --methods)
    items = list(base_dirs.items())
    if limit is not None:
        items = items[:limit]
    n = len(items)
    nrows, ncols = 1, max(n, 1)
    if panel is None:
        panel = (5, 3)
    figsize = (panel[0]*ncols, panel[1]*nrows)

    # Preload, shift by asymptotic reference (last row min), and compute mins
    loaded, mins = load_shifted(items)
    vmin_global = vmax_global = None
    if share_scale and len(mins) > 0:
        mfin = [m for m in mins if np.isfinite(m)]
        if mfin:
            vmin_global = float(np.min(mfin))
            if vmin_global > 0:
                vmin_global = 0.0
            vmax_global = -vmin_global

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False, constrained_layout=True)
    mesh = None
    for idx, (name, a, d, g, mloc) in enumerate(loaded):
        row, col = divmod(idx, ncols)
        ax = axes[row, col]

        ok = (a.size > 1) and (d.size > 1) and (g.size > 0)
        if ok:
            grid = np.ma.masked_invalid(g)
            if grid.mask.all() if hasattr(grid, 'mask') else False:
                ok = False
        if ok:
            try:
                # decide color scale
                if share_scale and (vmin_global is not None):
                    vmin_i, vmax_i = vmin_global, vmax_global
                else:
                    if np.isfinite(mloc):
                        vmin_i = mloc if mloc <= 0 else 0.0
                        vmax_i = -vmin_i
                    else:
                        vmin_i = vmax_i = None
                if non_uniform:
                    mesh = ax.pcolormesh(a, d, grid, cmap=cmap, shading="auto", vmin=vmin_i, vmax=vmax_i)
                else:
                    mesh = ax.imshow(grid, origin='lower', aspect='auto', cmap=cmap, vmin=vmin_i, vmax=vmax_i)
                    # label ticks by physical values while keeping uniform spacing
                    xt = np.linspace(0, grid.shape[1]-1, min(6, grid.shape[1])).astype(int)
                    yt = np.linspace(0, grid.shape[0]-1, min(6, grid.shape[0])).astype(int)
                    ax.set_xticks(xt); ax.set_yticks(yt)
                    ax.set_xticklabels([f"{a[i]:.0f}" for i in xt])
                    ax.set_yticklabels([f"{d[i]:.2f}" for i in yt])
            except Exception:
                ok = False
        if not ok:
            ax.text(0.5, 0.5, "No data", ha='center', va='center', transform=ax.transAxes)
            ax.set_xticks([]); ax.set_yticks([])

        ax.set_title(titles[idx] if titles else name, fontsize=10)
        ax.set_xlabel("Angle (deg)")
        ax.set_ylabel("Distance (Å)")

    # Remove unused panels
    for j in range(n, nrows*ncols):
        row, col = divmod(j, ncols)
        axes[row, col].axis("off")

    if suptitle:
        fig.suptitle(suptitle)
    if mesh is not None:
        fig.colorbar(mesh, ax=axes, orientation="vertical", fraction=0.02, pad=0.02, label="Energy")
    if savepath:
        Path(savepath).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(savepath, dpi=150, bbox_inches='tight')
    if show:
        plt.show()
    else:
        plt.close(fig)

    return loaded

if __name__ == "__main__":
    DEFAULT_ROOT = Path("/home/prokop/Desktop/CARBSIS/PEOPLE/Paolo/HbondFit_small_mols_2025_08_15")

    # Examples:
    #   python tests/tFitREQ/plot_Hbond_ref.py
    #   python tests/tFitREQ/plot_Hbond_ref.py --list
    #   python tests/tFitREQ/plot_Hbond_ref.py --methods DFT-b3lyp dftbp gfn2-xtb --file 'C4H3*D1-y.dat'
    #   python tests/tFitREQ/plot_Hbond_ref.py --non-uniform
    #   python tests/tFitREQ/plot_Hbond_ref.py --multi 100
    #   python tests/tFitREQ/plot_Hbond_ref.py --root /path/to/dataset --outdir plots2
    #   python tests/tFitREQ/plot_Hbond_ref.py --methods pm6 pm7 --file '...-y.dat' --no-share-scale
    #   python tests/tFitREQ/plot_Hbond_ref.py --panel 4 3
    #   python tests/tFitREQ/plot_Hbond_ref.py --common-scale  # force same vmin/vmax across methods
    #   python tests/tFitREQ/plot_Hbond_ref.py --line-plot --methods DFT-wb97m gfn2-xtb pm7 --file 'H2O-D1_CH2O-A1-z' --print-metrics
    #   python tests/tFitREQ/plot_Hbond_ref.py --metric-plot softclamp --multi 30 --methods DFT-wb97m gfn2-xtb pm7 --print-metrics

    p = argparse.ArgumentParser(description='Plot 2D H-bond maps across methods')
    p.add_argument('--root', type=Path, default=DEFAULT_ROOT, help='Dataset root containing data/ and ref/')
    p.add_argument('--limit', type=int, default=None, help='Max number of methods to plot (default: all)')
    p.add_argument('--cmap', type=str, default='seismic', help='Colormap')
    p.add_argument('--non-uniform', action='store_true', help='Use non-uniform axes (pcolormesh)')
    p.add_argument('--no-share-scale', action='store_true', help='Disable shared color scale (deprecated, use --common-scale)')
    p.add_argument('--common-scale', action='store_true', help='Use common symmetric color scale across methods')
    p.add_argument('--methods', nargs='*', default=None, help='Subset of method names to include')
    p.add_argument('--list', action='store_true', help='List available methods and exit')
    p.add_argument('--file', type=str, default="H2O-D1_CH2O-A1-z", help='Filename substring or glob (e.g. C4H3*D1-y.dat) to select per method')
    p.add_argument('--file-index', type=int, default=None, help='Index of .dat file (fallback if --file not matched)')
    p.add_argument('--multi', type=int, default=0, help='Plot first N files across methods; save to --outdir and do not show')
    p.add_argument('--outdir', type=str, default='plots', help='Output directory (under --root) for saved figures')
    p.add_argument('--panel', type=float, nargs=2, metavar=('W','H'), default=(3,5), help='Panel size in inches (W H), default 3 5')
    # metrics + reference for comparisons
    p.add_argument('--line-plot', action='store_true', help='Plot r_min(angle) and E_min(angle) lines instead of 2D map')
    p.add_argument('--ref-method', type=str, default='DFT-wb97m', help='Reference method name for metrics (display name)')
    p.add_argument('--clamp-start', type=float, default=0.0, help='Soft-clamp start (y1) for 2D error metric')
    p.add_argument('--clamp-limit', type=float, default=1.0, help='Soft-clamp limit (y2) for 2D error metric')
    p.add_argument('--print-metrics', action='store_true', help='Print metrics vs reference (2D soft-clamp RMSE or 1D RMSEs)')
    p.add_argument('--min-rmax', type=float, default=None, help='Max r allowed for minima; if r_min>this, set to NaN (filter repulsive angles)')
    p.add_argument('--metric-plot', '--metric_plot', type=str, default=None, help='Compute and plot metric across first --multi systems (e.g., softclamp/softclam)')
    p.add_argument('--metric-show', action='store_true', help='Show metric figures interactively (in addition to saving)')
    args = p.parse_args()

    base = args.root / "data"
    ref  = args.root / "ref"

    # discover methods -> lists of files
    discovered = []
    files_map = {}
    for folder in [base, ref]:
        if not folder.exists():
            continue
        for sub in sorted(folder.iterdir()):
            if not sub.is_dir():
                continue
            disp = label(sub.name)
            discovered.append(disp)
            datfiles = sorted(sub.rglob("*.dat"))
            if datfiles:
                files_map[disp] = [str(p) for p in datfiles]

    if args.list:
        print("Methods:")
        for m in sorted(set(discovered), key=lambda x: (METHODS.index(x) if x in METHODS else 1e9, x)):
            print(m)
        sys.exit(0)

    # choose methods and order
    if args.methods:
        use = [m for m in args.methods if m in files_map]
    else:
        use = [m for m in METHODS if m in files_map]
        if not use:
            use = sorted(files_map.keys())

    if not use:
        print(f"No .dat files found under '{base}' or '{ref}'. Use --root to set dataset path or --list to inspect methods.")
        sys.exit(1)

    lim = None if (args.limit is None or args.limit <= 0) else args.limit
    pan = tuple(args.panel) if args.panel else None
    common = args.common_scale or (not args.no_share_scale)

    # MULTI: plot first N files and save figures
    if args.multi:
        # Metric plot mode over multiple systems
        metric = (args.metric_plot or '').lower()
        if metric in ('softclamp','softclam'):
            outdir = args.root / args.outdir
            counts = [len(files_map[m]) for m in use if m in files_map]
            if not counts:
                print("No files to plot in selected methods."); sys.exit(1)
            max_i = min(counts) if args.multi <= 0 else min(min(counts), args.multi)
            stems = []
            methods_eval = [m for m in use if m != args.ref_method]
            metrics = {m: [] for m in methods_eval}
            for i in range(max_i):
                filemap = { m: files_map[m][i] for m in use if len(files_map[m]) > i }
                first = next(iter(filemap.values()), None)
                stems.append(Path(first).stem if first else f"idx_{i:03d}")
                loaded, _ = load_shifted(list(filemap.items()))
                name2grid = {name: g for (name, a, d, g, mloc) in loaded}
                if args.ref_method not in name2grid:
                    for m in methods_eval: metrics[m].append(np.nan)
                    continue
                gref = name2grid[args.ref_method]
                for m in methods_eval:
                    g = name2grid.get(m, None)
                    val = softclamp_rmse(g, gref, args.clamp_start, args.clamp_limit) if g is not None else np.nan
                    metrics[m].append(val)
            # plot lines per method over systems
            fig, ax = plt.subplots(1,1, figsize=(max(8, 0.35*len(stems)), 4), constrained_layout=True)
            xs = np.arange(len(stems))
            for m in methods_eval:
                ax.plot(xs, metrics[m], label=m)
            ax.set_xticks(xs); ax.set_xticklabels(stems, rotation=90)
            ax.set_ylabel('SoftClamp RMSE')
            ax.grid(True, linestyle='--', alpha=0.4)
            ax.legend(ncol=2, fontsize=8)
            outdir.mkdir(parents=True, exist_ok=True)
            p_lines = outdir / f"metrics_softclamp_lines.png"
            fig.savefig(p_lines, dpi=150, bbox_inches='tight')
            print(f"Saved: {p_lines}")
            if args.metric_show:
                plt.show()
            else:
                plt.close(fig)
            # aggregate plot (mean over systems)
            agg = {m: float(np.nanmean(metrics[m])) if len(metrics[m])>0 else np.nan for m in methods_eval}
            fig, ax = plt.subplots(1,1, figsize=(max(6, 0.25*len(methods_eval)), 3.5), constrained_layout=True)
            ms = list(methods_eval); vals = [agg[m] for m in ms]
            ax.bar(ms, vals)
            ax.set_ylabel('Mean SoftClamp RMSE'); ax.set_title('Aggregate over systems')
            ax.grid(True, axis='y', linestyle='--', alpha=0.4)
            p_agg = outdir / f"metrics_softclamp_agg.png"
            fig.savefig(p_agg, dpi=150, bbox_inches='tight')
            print(f"Saved: {p_agg}")
            if args.metric_show:
                plt.show()
            else:
                plt.close(fig)
            if args.print_metrics:
                print("SoftClamp RMSE per system:")
                for m in methods_eval:
                    print(m, '->', ', '.join(f"{v:.4f}" if np.isfinite(v) else 'nan' for v in metrics[m]))
                print("Aggregate mean:")
                for m in methods_eval:
                    print(f"{m}: {agg[m]:.4f}")
            sys.exit(0)

        outdir = args.root / args.outdir
        # common length assumption; use min length across methods
        counts = [len(files_map[m]) for m in use if m in files_map]
        if not counts:
            print("No files to plot in selected methods.")
            sys.exit(1)
        max_i = min(counts) if args.multi <= 0 else min(min(counts), args.multi)
        for i in range(max_i):
            filemap = { m: files_map[m][i] for m in use if len(files_map[m]) > i }
            first = next(iter(filemap.values())) if filemap else None
            stem = Path(first).stem if first else f"idx_{i:03d}"
            fname = Path(first).name if first else f"idx_{i:03d}.dat"
            savepath = outdir / f"{i:03d}_{stem}{'_lines' if args.line_plot else ''}.png"
            print(f"Plotting: {fname}")
            if args.line_plot:
                items = list(filemap.items())
                loaded, _ = load_shifted(items)
                curves = plot_min_lines(loaded, titles=use, panel=pan, show=False, savepath=str(savepath), suptitle=fname, rmax=args.min_rmax)
                if args.print_metrics:
                    if args.ref_method in curves:
                        aR, rR, eR = curves[args.ref_method]
                        for m in curves:
                            if m == args.ref_method: continue
                            a, r, e = curves[m]
                            r_rmse = align_and_rmse(aR, rR, a, r)
                            e_rmse = align_and_rmse(aR, eR, a, e)
                            print(f"RMSE[{m} vs {args.ref_method}] r_min={r_rmse:.4f} E_min={e_rmse:.4f}")
                    else:
                        print(f"Reference '{args.ref_method}' not in curves; skip metrics.")
            else:
                loaded = plot_all_methods(filemap, titles=use, panel=pan, cmap=args.cmap,
                                          limit=lim, share_scale=common,
                                          non_uniform=args.non_uniform, show=False, savepath=str(savepath), suptitle=fname)
                if args.print_metrics and loaded:
                    name2grid = {name: g for (name, a, d, g, mloc) in loaded}
                    if args.ref_method in name2grid:
                        gref = name2grid[args.ref_method]
                        for name, g in name2grid.items():
                            if name == args.ref_method: continue
                            sc = softclamp_rmse(g, gref, args.clamp_start, args.clamp_limit)
                            print(f"SoftClamp-RMSE[{name} vs {args.ref_method}] = {sc:.4f}")
                    else:
                        print(f"Reference '{args.ref_method}' not found among plotted methods; skip metrics.")
        sys.exit(0)

    # SINGLE: pick file by name/pattern (or index fallback) and show
    filemap = {}
    for m in use:
        files = files_map.get(m, [])
        chosen = pick_file(files, pattern=args.file, index=args.file_index)
        if chosen is not None:
            filemap[m] = chosen
    # suptitle: show picked filename (use first available)
    fname = Path(next(iter(filemap.values()))).name if filemap else None
    if fname:
        print(f"Plotting: {fname}")
    if args.line_plot:
        items = list(filemap.items())
        loaded, _ = load_shifted(items)
        curves = plot_min_lines(loaded, titles=use, panel=pan, show=True, savepath=None, suptitle=fname, rmax=args.min_rmax)
        if args.print_metrics:
            if args.ref_method in curves:
                aR, rR, eR = curves[args.ref_method]
                for m in curves:
                    if m == args.ref_method: continue
                    a, r, e = curves[m]
                    r_rmse = align_and_rmse(aR, rR, a, r)
                    e_rmse = align_and_rmse(aR, eR, a, e)
                    print(f"RMSE[{m} vs {args.ref_method}] r_min={r_rmse:.4f} E_min={e_rmse:.4f}")
            else:
                print(f"Reference '{args.ref_method}' not in curves; skip metrics.")
    else:
        loaded = plot_all_methods(filemap, titles=use, panel=pan, cmap=args.cmap,
                                  limit=lim, share_scale=common,
                                  non_uniform=args.non_uniform, show=True, savepath=None, suptitle=fname)
        if args.print_metrics and loaded:
            name2grid = {name: g for (name, a, d, g, mloc) in loaded}
            if args.ref_method in name2grid:
                gref = name2grid[args.ref_method]
                for name, g in name2grid.items():
                    if name == args.ref_method: continue
                    sc = softclamp_rmse(g, gref, args.clamp_start, args.clamp_limit)
                    print(f"SoftClamp-RMSE[{name} vs {args.ref_method}] = {sc:.4f}")
            else:
                print(f"Reference '{args.ref_method}' not found among plotted methods; skip metrics.")


# if __name__ == "__main__":
#     file_path = "/mnt/data/C4H3NO2-A1_C4H3NO2-D1-y.dat"
#     angles, distances, energy_grid = load_energy_map(file_path)
#     plot_energy_map(angles, distances, energy_grid, title="B3LYP")