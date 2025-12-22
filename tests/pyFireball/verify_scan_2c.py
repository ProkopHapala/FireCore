
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

# Add pyBall to path
sys.path.append(os.path.abspath('../../'))
import pyBall.FireCore as fc
from pyBall.FireballOCL.OCL_Hamiltonian import OCL_Hamiltonian as ocl

def _extract_components(mat):
    """Extract requested orbital pairs from a 4x4 block (s, py, pz, px)."""
    # Indices in Ortega convention: 0:s, 1:py, 2:pz, 3:px
    pairs = [
        ("s-s",   (0,0)),
        ("s-px",  (0,3)),
        ("s-py",  (0,1)),
        ("s-pz",  (0,2)),
        ("px-px", (3,3)),
        ("py-py", (1,1)),
        ("pz-pz", (2,2)),
        ("pz-py", (2,1)),
        ("px-pz", (3,2)),
        ("py-pz", (1,2)),
    ]
    out = []
    for _, (i,j) in pairs:
        if i < mat.shape[0] and j < mat.shape[1]:
            out.append(mat[i,j])
        else:
            out.append(np.nan)
    return np.array(out), pairs

# Report only failing components (> threshold) with per-component extrema
def _report_diffs(tag, xs, fs, ps, pairs, threshold=1e-4):
    diff = np.abs(fs - ps)
    vmax_diff = np.nanmax(diff)
    # Always print per-component extrema and max|diff|
    for i, (label, _) in enumerate(pairs):
        vmin_f, vmax_f = np.nanmin(fs[:, i]), np.nanmax(fs[:, i])
        vmin_p, vmax_p = np.nanmin(ps[:, i]), np.nanmax(ps[:, i])
        vmax_d = np.nanmax(diff[:, i])
        print(f"  [{label}] Fortran vmin/vmax {vmin_f: .6e}/{vmax_f: .6e}  PyOCL vmin/vmax {vmin_p: .6e}/{vmax_p: .6e}  max|diff|={vmax_d: .6e}")
    failing = []
    for i, (label, _) in enumerate(pairs):
        vmin_f, vmax_f = np.nanmin(fs[:, i]), np.nanmax(fs[:, i])
        vmin_p, vmax_p = np.nanmin(ps[:, i]), np.nanmax(ps[:, i])
        vmax_d = np.nanmax(diff[:, i])
        if vmax_d > threshold:
            failing.append((label, vmax_d, (vmin_f, vmax_f), (vmin_p, vmax_p), i))
    print(f"[2c:{tag}] max|diff|={vmax_diff: .6e}; failing components (>{threshold}):")
    if not failing:
        print("  none")
        return
    for label, vmax_d, (vmin_f, vmax_f), (vmin_p, vmax_p), _ in failing:
        print(f"  [{label}] max|diff|={vmax_d: .6e}  Fortran vmin/vmax {vmin_f: .6e}/{vmax_f: .6e}  PyOCL vmin/vmax {vmin_p: .6e}/{vmax_p: .6e}")
    # Dump only failing curves
    print("[WARN] dumping curves for failing components:")
    for label, _, _, _, idx in failing:
        print(f"  curves {label} (x, Fortran, PyOCL):")
        for x, vf, vp in zip(xs, fs[:, idx], ps[:, idx]):
            print(f"    x={x: .6f}  F={vf: .6e}  P={vp: .6e}")

def plot_scan_2c(root, nz1, nz2, interaction, isub, rs, applyRotation=True):
    print(f"Scanning 2-center {root} (interaction={interaction}, isub={isub})...")
    
    fs = []
    ps = []
    
    # Fortran expects species indices (1..nspecies). 
    # Since we initialize with one species, it's index 1.
    in1, in2 = 1, 1
    # Fortran conventions (see fortran/ASSEMBLERS/assemble_2c.f90 and INTERACTIONS/doscentros.f90):
    # - overlap/kinetic: in3=in2
    # - vna_atom (interaction=4): in3=in1
    # - vna_ontopl (2) and vna_ontopr (3): in3=in2
    in3 = in2
    if interaction == 4:
        in3 = in1
    
    for r in rs:
        dR = np.array([0, 0, r], dtype=np.float64)
        # Fortran
        mat_f = fc.scanHamPiece2c(interaction, isub, in1, in2, in3, dR, applyRotation=applyRotation)
        
        # PyOpenCL
        mat_p = ham.scanHamPiece2c(root, nz1, nz2, dR, applyRotation=applyRotation)
        
        if mat_p is not None:
            v_f, pairs = _extract_components(mat_f)
            v_p, _     = _extract_components(mat_p)
            fs.append(v_f)
            ps.append(v_p)
        else:
            # preserve length of vector for downstream stats
            if mat_f is not None:
                v_f, pairs = _extract_components(mat_f)
            else:
                v_f = np.full(10, np.nan)
                pairs = None
            fs.append(v_f)
            ps.append(np.full_like(v_f, np.nan))
            
    fs = np.array(fs)
    ps = np.array(ps)

    _report_diffs(root, rs, fs, ps, pairs)
    
    plt.figure(figsize=(8, 5))
    plt.plot(rs, fs[:, 0], 'r:', lw=1.5, label='Fortran (s-s)')
    plt.plot(rs, ps[:, 0], 'r-', lw=0.5, label='PyOpenCL (s-s)')
    # Indices must match _extract_components() order:
    # 0:s-s, 1:s-px, 2:s-py, 3:s-pz, 4:px-px, 5:py-py, 6:pz-pz, ...
    plt.plot(rs, fs[:, 1], 'g:', lw=1.5, label='Fortran (s-px)')
    plt.plot(rs, ps[:, 1], 'g-', lw=0.5, label='PyOpenCL (s-px)')
    plt.plot(rs, fs[:, 2], 'm:', lw=1.5, label='Fortran (s-py)')
    plt.plot(rs, ps[:, 2], 'm-', lw=0.5, label='PyOpenCL (s-py)')
    plt.plot(rs, fs[:, 3], 'b:', lw=1.5, label='Fortran (s-pz)')
    plt.plot(rs, ps[:, 3], 'b-', lw=0.5, label='PyOpenCL (s-pz)')
    plt.plot(rs, fs[:, 4], 'c:', lw=1.5, label='Fortran (px-px)')
    plt.plot(rs, ps[:, 4], 'c-', lw=0.5, label='PyOpenCL (px-px)')
    plt.plot(rs, fs[:, 5], 'y:', lw=1.5, label='Fortran (py-py)')
    plt.plot(rs, ps[:, 5], 'y-', lw=0.5, label='PyOpenCL (py-py)')
    plt.plot(rs, fs[:, 6], 'k:', lw=1.5, label='Fortran (pz-pz)')
    plt.plot(rs, ps[:, 6], 'k-', lw=0.5, label='PyOpenCL (pz-pz)')
    # Cross terms (may be zero by symmetry for this geometry but included for completeness)
    plt.plot(rs, fs[:, 7], color='0.5', linestyle=':', lw=1.0, label='Fortran (pz-py)')
    plt.plot(rs, ps[:, 7], color='0.5', linestyle='-', lw=0.5, label='PyOpenCL (pz-py)')
    plt.plot(rs, fs[:, 8], color='0.3', linestyle=':', lw=1.0, label='Fortran (px-pz)')
    plt.plot(rs, ps[:, 8], color='0.3', linestyle='-', lw=0.5, label='PyOpenCL (px-pz)')
    plt.plot(rs, fs[:, 9], color='0.7', linestyle=':', lw=1.0, label='Fortran (py-pz)')
    plt.plot(rs, ps[:, 9], color='0.7', linestyle='-', lw=0.5, label='PyOpenCL (py-pz)')
    plt.title(f"2-Center Scan: {root} ({nz1}-{nz2})")
    plt.xlabel("Distance [A]")
    plt.ylabel("Value")
    plt.legend()
    plt.grid(True)
    plt.savefig(f"scan_2c_{root}.png")
    print(f"Saved plot scan_2c_{root}.png")

def plot_scan_2c_angular(root, nz1, nz2, interaction, isub, r, thetas):
    print(f"Scanning 2-center angular {root} (r={r})...")
    
    fs = []
    ps = []
    
    in1, in2 = 1, 1
    in3 = in2
    if interaction == 4:
        in3 = in1
    
    pairs = [
        ("s-px",   (0,3)),
        ("s-py",   (0,1)),
        ("s-pz",   (0,2)),
        ("px-px",  (3,3)),
        ("py-py",  (1,1)),
        ("pz-pz",  (2,2)),
        ("pz-py",  (2,1)),
        ("px-pz",  (3,2)),
        ("py-pz",  (1,2)),
    ]
    for th in thetas:
        dR = np.array([0, r*np.sin(th), r*np.cos(th)], dtype=np.float64)
        # Fortran
        mat_f = fc.scanHamPiece2c(interaction, isub, in1, in2, in3, dR, applyRotation=True)
        # PyOpenCL
        mat_p = ham.scanHamPiece2c(root, nz1, nz2, dR, applyRotation=True)
        
        if mat_p is not None:
            v_f, _ = _extract_components(mat_f)
            v_p, _ = _extract_components(mat_p)
            fs.append(v_f); ps.append(v_p)
        else:
            v_f, _ = _extract_components(mat_f) if mat_f is not None else (np.full(10, np.nan), None)
            fs.append(v_f); ps.append(np.full_like(v_f, np.nan))
    fs = np.array(fs); ps = np.array(ps)
    _report_diffs(f"ang:{root}", thetas, fs, ps, pairs)

    plt.figure(figsize=(10, 6))
    colors = ['r','g','b','k','m','c','0.5','0.3','0.7']
    labels = ["s-px","s-py","s-pz","px-px","py-py","pz-pz","pz-py","px-pz","py-pz"]
    for i, lab in enumerate(labels):
        plt.plot(thetas, fs[:, i], linestyle=':', linewidth=1.5, color=colors[i], label=f'Fortran ({lab})')
        plt.plot(thetas, ps[:, i], linestyle='-', linewidth=0.5, color=colors[i], label=f'PyOpenCL ({lab})')
    plt.title(f"2-Center Angular Scan: {root} (r={r})")
    plt.xlabel("Theta [rad]")
    plt.ylabel("Value")
    plt.legend()
    plt.grid(True)
    plt.savefig(f"scan_2c_angular_{root}.png")
    print(f"Saved plot scan_2c_angular_{root}.png")

if __name__ == "__main__":
    fdata_dir = "./Fdata"
    print(f"Using Fdata from: {fdata_dir}")
    nz1, nz2 = 6, 6  # Use Carbon for p-orbitals
    
    # 1. Initialize Fortran
    print("Initializing Fortran...")
    poss = [[0,0,0],[0,0,1.0]]
    #poss = [[0,0,0],[0,0.8,0.8]]
    fc.initialize(atomType=np.array([nz1, nz2], dtype=np.int32), atomPos=np.array(poss, dtype=np.float64))
    
    # 2. Initialize PyOpenCL
    print("Initializing PyOpenCL...")
    ham = ocl(fdata_dir)
    ham.prepare_splines([nz1, nz2])
    ham.prepare_data_3c([nz1, nz2])
    
    rs = np.linspace(0.1, 5.0, 50)
    
    # 2-Center Scans
    # interaction: 1=Overlap, 13=Kinetic, 4=Vna (atom)
    plot_scan_2c('overlap', nz1, nz2, 1, 0, rs, applyRotation=False)
    plot_scan_2c('kinetic', nz1, nz2, 13, 0, rs, applyRotation=False)
    plot_scan_2c('vna',     nz1, nz2, 4, 0, rs, applyRotation=False)
    plot_scan_2c('vnl',     nz1, nz2, 5, 0, rs, applyRotation=False)

    # Vna decomposed tables (Fortran interactions 2 and 3)
    plot_scan_2c('vna_ontopl_00', nz1, nz2, 2, 0, rs, applyRotation=False)
    plot_scan_2c('vna_ontopr_00', nz1, nz2, 3, 0, rs, applyRotation=False)
    
    # Angular Scan
    thetas = np.linspace(0, 2*np.pi, 50)
    plot_scan_2c_angular('overlap', nz1, nz2, 1, 0, 1.5, thetas)
    plot_scan_2c_angular('kinetic', nz1, nz2, 13, 0, 1.5, thetas)
    plot_scan_2c_angular('vna',     nz1, nz2, 4, 0, 1.5, thetas)
    plot_scan_2c_angular('vnl',     nz1, nz2, 5, 0, 1.5, thetas)
    plot_scan_2c_angular('vna_ontopl_00', nz1, nz2, 2, 0, 1.5, thetas)
    plot_scan_2c_angular('vna_ontopr_00', nz1, nz2, 3, 0, 1.5, thetas)
    print("2-center verification complete.")
