
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import argparse

# Add pyBall to path
sys.path.append(os.path.abspath('../../'))
import pyBall.FireCore as fc
from pyBall.FireballOCL.OCL_Hamiltonian import OCL_Hamiltonian as ocl

def _get_bcna_grid(parser, nz1, nz2, nz3, isorp=0):
    """Return (hx, hy, nx, ny) for the bcna tables."""
    _, d3c = parser.load_species_data([nz1, nz2, nz3])
    key = ('bcna', 1, isorp, nz1, nz2, nz3)
    rec = d3c.get(key)
    if rec is None:
        raise RuntimeError(f"Missing bcna record for ({nz1},{nz2},{nz3}) isorp={isorp}")
    hx = rec['xmax'] / (rec['numx'] - 1)
    hy = rec['ymax'] / (rec['numy'] - 1)
    nx = rec['numx']
    ny = rec['numy']
    return hx, hy, nx, ny


def _dump_bcna_table_stats(parser, nz1, nz2, nz3, isorp=0, itheta=1, i_nz=0, max_print=8):
    _, d3c = parser.load_species_data([nz1, nz2, nz3])
    key = ('bcna', itheta, isorp, nz1, nz2, nz3)
    rec = d3c.get(key)
    if rec is None:
        print(f"[PY-TABLE] Missing bcna record {key}")
        return
    data = rec['data']
    if data.ndim != 3:
        print(f"[PY-TABLE] Unexpected data.ndim={data.ndim} for {key}")
        return
    if i_nz < 0 or i_nz >= data.shape[2]:
        print(f"[PY-TABLE] i_nz out of range: i_nz={i_nz} shape={data.shape}")
        return
    sl = data[:, :, i_nz]
    print(f"[PY-TABLE] key={key} shape={data.shape} slice(i_nz={i_nz}) vmin/vmax={np.nanmin(sl): .9e}/{np.nanmax(sl): .9e}")
    vals = []
    for iy in range(min(max_print, sl.shape[0])):
        vals.append(float(sl[iy, 0]))
    print(f"[PY-TABLE] first-column y=0..{min(max_print, sl.shape[0]) - 1} at x=0: {vals}")



def plot_scan_3c(root, nz1, nz2, nz3, interaction, isorp, rs, applyRotation=True, mode="full", verbosity=0, run_ocl_raw=False):
    print(f"Scanning 3-center {root} (interaction={interaction}, isorp={isorp})...")
    
    fs = []
    ps = []
    xs = []
    ys = []
    dRj_list = []
    dRk_list = []
    
    in1, in2, in3 = 1, 1, 1
    hx = hy = nx = ny = None
    
    for r in rs:
        # Arrange in a triangle for 3-center scan
        # Basis 1 at (0,0,0), Basis 2 at (r, 0, 0), NA at (r/2, r/2, 0)
        dRj = np.array([r, 0, 0], dtype=np.float64)
        dRk = np.array([r*0.5, r*0.5, 0], dtype=np.float64)
        dRj_list.append(dRj)
        dRk_list.append(dRk)
        # Fortran geometry: y=|r2-r1|, x=|rna-(r1+r2)/2|
        r21 = dRj
        y = np.linalg.norm(r21)
        rnabc = dRk - 0.5*(np.array([0,0,0])+dRj)
        x = np.linalg.norm(rnabc)
        xs.append(x); ys.append(y)
        
        # Fortran / reference: raw vs full
        if mode != "raw":
            mat_f = fc.scanHamPiece3c(interaction, isorp, in1, in2, in3, dRj, dRk, applyRotation=applyRotation)
            fs.append(mat_f[0, 0])
            
            mat_p = ham.scanHamPiece3c(root, nz1, nz2, nz3, dRj, dRk, applyRotation=applyRotation)
            if mat_p is not None:
                ps.append(mat_p[0]) # first element of hlist
            else:
                ps.append(np.nan)
            
    xs = np.array(xs, dtype=np.float64)
    ys = np.array(ys, dtype=np.float64)
    dRjs = np.ascontiguousarray(dRj_list, dtype=np.float64)
    dRks = np.ascontiguousarray(dRk_list, dtype=np.float64)

    if mode == "raw":
        npoints = dRjs.shape[0]

        raw_fortran = fc.scanHamPiece3c_raw_batch(interaction, isorp, in1, in2, in3, dRjs, dRks)
        fs = raw_fortran[:, :, 0]   # [npoints, 5]

        if not run_ocl_raw:
            raise RuntimeError("raw mode now compares Fortran vs OpenCL; pass --debug-ocl-raw 1")

        raw_ocl = ham.scanHamPiece3c_raw_batch(root, nz1, nz2, nz3, dRjs.astype(np.float32), dRks.astype(np.float32))
        if raw_ocl is None:
            raise RuntimeError("OpenCL raw batch returned None")
        ps = raw_ocl[:, :, 0]       # [npoints, 5]

        if int(args.debug_geom) > 0:
            _dump_bcna_table_stats(ham.parser, nz1, nz2, nz3, isorp=isorp, itheta=1, i_nz=0, max_print=8)

        hx, hy, nx, ny = _get_bcna_grid(ham.parser, nz1, nz2, nz3, isorp)

        nprint = int(args.debug_geom)
        if nprint <= 0:
            nprint = npoints
        nprint = min(nprint, npoints)

        print(f"[RAW-CMP] npoints={npoints} hx={hx:.6f} hy={hy:.6f} nx={nx} ny={ny}")
        print("[RAW-CMP] columns: ip x y F0 O0 d0")
        for ip in range(nprint):
            f0 = float(fs[ip, 0])
            o0 = float(ps[ip, 0])
            d0 = abs(f0 - o0)
            print(f"[RAW-CMP] ip={ip:03d} x={xs[ip]:.6f} y={ys[ip]:.6f}  F0={f0: .9e}  O0={o0: .9e}  d0={d0: .3e}")

        diff = np.abs(fs - ps)
        vmax_diff = np.nanmax(diff)
        worst_idx = np.unravel_index(np.nanargmax(diff), diff.shape)
        worst_r = rs[worst_idx[0]]
        worst_comp = worst_idx[1]
        worst_f = fs[worst_idx]
        worst_p = ps[worst_idx]
        if np.isnan(diff).any():
            nan_locs = np.argwhere(np.isnan(diff))
            print(f"[3c_raw:{root}] WARNING NaN in diff at indices (r_idx, comp): {nan_locs.tolist()}")
            print("Fortran raw bcna_0k:\n", fs)
            print("OpenCL raw bcna_0k:\n", ps)
        print(f"[3c_raw:{root}] SUMMARY max|F-P|={vmax_diff: .6e} at r={worst_r:.4f} L{worst_comp} "
              f"(F={worst_f: .6e}, O={worst_p: .6e})")
        print(f"[3c_raw:{root}] vmin/vmax Fortran: {np.nanmin(fs): .6e}/{np.nanmax(fs): .6e}  "
              f"OpenCL: {np.nanmin(ps): .6e}/{np.nanmax(ps): .6e}  max|diff|={vmax_diff: .6e}")
        for k in range(fs.shape[1]):
            vmin_f = np.nanmin(fs[:, k])
            vmax_f = np.nanmax(fs[:, k])
            vmin_p = np.nanmin(ps[:, k])
            vmax_p = np.nanmax(ps[:, k])
            vmax_d = np.nanmax(diff[:, k])
            print(f"  [L{k}] Fortran vmin/vmax {vmin_f: .6e}/{vmax_f: .6e}  "
                  f"OpenCL vmin/vmax {vmin_p: .6e}/{vmax_p: .6e}  max|diff|={vmax_d: .6e}")

        plt.figure(figsize=(8, 5))
        colors = ['r','g','b','m','c']
        for k in range(fs.shape[1]):
            plt.plot(rs, fs[:,k], ':', lw=1.5, color=colors[k], label=f'Fortran L{k}')
            plt.plot(rs, ps[:,k], '-', lw=0.5, color=colors[k], label=f'OpenCL  L{k}')
        plt.title(f"3-Center RAW bcna: {root} ({nz1}-{nz2}-{nz3})")
        plt.xlabel("Distance Parameter [A]")
        plt.ylabel("bcna_0k value")
        plt.legend()
        plt.grid(True)
        fname = f"scan_3c_{root}_raw.png"
        plt.savefig(fname)
        print(f"Saved plot {fname}")
    else:
        fs = np.array(fs, dtype=np.float64)
        ps = np.array(ps, dtype=np.float64)
        diff = np.abs(fs - ps)
        vmax_diff = np.nanmax(diff)
        if np.isnan(diff).any():
            nan_locs = np.argwhere(np.isnan(diff))
            print(f"[3c:{root}] WARNING NaN in diff at indices (r_idx): {nan_locs.tolist()}")
            print("Fortran curves (raw):\n", fs)
            print("PyOCL curves (raw):\n", ps)
        print(f"[3c:{root}] vmin/vmax Fortran: {np.nanmin(fs): .6e}/{np.nanmax(fs): .6e}  "
              f"PyOCL: {np.nanmin(ps): .6e}/{np.nanmax(ps): .6e}  max|diff|={vmax_diff: .6e}")
        # Debug geometry sampling
        print(f"[3c:{root}] sample x,y (first 5):")
        for i in range(min(5, len(xs))):
            r_val = float(rs[i])
            x_val = float(xs[i])
            y_val = float(ys[i])
            f_val = float(fs[i])
            p_val = float(ps[i])
            print(f"  r={r_val:.4f}  x={x_val:.6f}  y={y_val:.6f}  F={f_val: .6e}  P={p_val: .6e}")
        diff_threshold = 1e-4
        if vmax_diff > diff_threshold:
            print(f"[WARN] 3c {root} diff exceeds threshold {diff_threshold: .2e}; dumping curves")
            print("Fortran curves:\n", fs)
            print("PyOCL curves:\n", ps)

        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(rs, fs, 'g:', lw=1.0, label='Fortran')
        ax.plot(rs, ps, 'k-', lw=0.5, label='PyOpenCL')
        ax.set_title(f"3-Center Scan: {root} ({nz1}-{nz2}-{nz3})")
        ax.set_xlabel("Distance Parameter [A]")
        ax.set_ylabel("Value")
        ax.legend()
        ax.grid(True)
        fname = f"scan_3c_{root}.png"
        plt.savefig(fname)
        print(f"Saved plot {fname}")

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--fdata", default="./Fdata")
    ap.add_argument("--nz1", type=int, default=1)
    ap.add_argument("--nz2", type=int, default=1)
    ap.add_argument("--mode", choices=["full","raw"], default="full", help="full=Legendre+recover+rotate via Fortran/OCL; raw=bare bcna_0k first component")
    ap.add_argument("--rmin", type=float, default=0.001)
    ap.add_argument("--rmax", type=float, default=4.0)
    ap.add_argument("--rnum", type=int, default=100)
    ap.add_argument("--debug-geom", type=int, default=1, help="print first few (r,x,y) with Fortran/Fdata raw values")
    ap.add_argument("--verbosity",  type=int, default=4, help="Fireball verbosity for debug output")
    ap.add_argument("--debug-ocl-raw", type=int, default=1, help="Invoke PyOpenCL raw kernel to print its debug info")
    args = ap.parse_args()

    fdata_dir = args.fdata
    print(f"Using Fdata from: {fdata_dir}")
    nz1, nz2 = args.nz1, args.nz2
    
    # 1. Initialize Fortran
    print("Initializing Fortran...")
    fc.initialize(atomType=np.array([nz1, nz2], dtype=np.int32),  atomPos=np.array([[0,0,0],[0,0,1.0]], dtype=np.float64))
    fc.setVerbosity(args.verbosity, 0)
    
    # 2. Initialize PyOpenCL
    print("Initializing PyOpenCL...")
    ham = ocl(fdata_dir)
    ham.prepare_splines([nz1, nz2])
    ham.prepare_data_3c([nz1, nz2])
    
    rs_3c = np.linspace(args.rmin, args.rmax, args.rnum)
    plot_scan_3c('bcna', nz1, nz2, nz1, 1, 0, rs_3c, applyRotation=False, mode=args.mode, verbosity=args.verbosity, run_ocl_raw=args.debug_ocl_raw)
    
    print("Verification scans complete.")
