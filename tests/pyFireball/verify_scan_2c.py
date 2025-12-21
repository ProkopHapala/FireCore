
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

def plot_scan_2c(root, nz1, nz2, interaction, isub, rs, applyRotation=True):
    print(f"Scanning 2-center {root} (interaction={interaction}, isub={isub})...")
    
    fs = []
    ps = []
    
    # Fortran expects species indices (1..nspecies). 
    # Since we initialize with one species, it's index 1.
    in1, in2 = 1, 1
    
    for r in rs:
        dR = np.array([0, 0, r], dtype=np.float64)
        # Fortran
        mat_f = fc.scanHamPiece2c(interaction, isub, in1, in2, in2, dR, applyRotation=applyRotation)
        
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

    # DEBUG: print extrema and max diff to catch sign issues
    diff = np.abs(fs - ps)
    vmax_diff = np.max(diff)
    if np.isnan(diff).any():
        nan_locs = np.argwhere(np.isnan(diff))
        print(f"[2c:{root}] WARNING NaN in diff at indices (r_idx, component): {nan_locs.tolist()}")
        print("Fortran curves (raw):\n", fs)
        print("PyOCL curves (raw):\n", ps)
    print(f"[2c:{root}] vmin/vmax Fortran: {fs.min(): .6e}/{fs.max(): .6e}  PyOCL: {ps.min(): .6e}/{ps.max(): .6e}  max|diff|={vmax_diff: .6e}")
    # Component-wise extrema for clarity
    for i, (label, _) in enumerate(pairs):
        vmin_f, vmax_f = fs[:, i].min(), fs[:, i].max()
        vmin_p, vmax_p = ps[:, i].min(), ps[:, i].max()
        vmax_d = diff[:, i].max()
        print(f"  [{label}] Fortran vmin/vmax {vmin_f: .6e}/{vmax_f: .6e}  PyOCL vmin/vmax {vmin_p: .6e}/{vmax_p: .6e}  max|diff|={vmax_d: .6e}")
    diff_threshold = 1e-4
    if vmax_diff > diff_threshold:
        print(f"[WARN] 2c {root} diff exceeds threshold {diff_threshold: .2e}; dumping curves for debugging")
        print("Fortran curves:\n", fs)
        print("PyOCL curves:\n", ps)
    
    plt.figure(figsize=(8, 5))
    plt.plot(rs, fs[:, 0], 'r:', lw=1.5, label='Fortran (s-s)')
    plt.plot(rs, ps[:, 0], 'r-', lw=0.5, label='PyOpenCL (s-s)')
    plt.plot(rs, fs[:, 3], 'g:', lw=1.5, label='Fortran (s-px)')
    plt.plot(rs, ps[:, 3], 'g-', lw=0.5, label='PyOpenCL (s-px)')
    plt.plot(rs, fs[:, 2], 'm:', lw=1.5, label='Fortran (s-py)')
    plt.plot(rs, ps[:, 2], 'm-', lw=0.5, label='PyOpenCL (s-py)')
    plt.plot(rs, fs[:, 1], 'b:', lw=1.5, label='Fortran (s-pz)')
    plt.plot(rs, ps[:, 1], 'b-', lw=0.5, label='PyOpenCL (s-pz)')
    plt.plot(rs, fs[:, 4], 'c:', lw=1.5, label='Fortran (px-px)')
    plt.plot(rs, ps[:, 4], 'c-', lw=0.5, label='PyOpenCL (px-px)')
    plt.plot(rs, fs[:, 5], 'y:', lw=1.5, label='Fortran (py-py)')
    plt.plot(rs, ps[:, 5], 'y-', lw=0.5, label='PyOpenCL (py-py)')
    plt.plot(rs, fs[:, 6], 'k:', lw=1.5, label='Fortran (pz-pz)')
    plt.plot(rs, ps[:, 6], 'k-', lw=0.5, label='PyOpenCL (pz-pz)')
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
    
    for th in thetas:
        dR = np.array([0, r*np.sin(th), r*np.cos(th)], dtype=np.float64)
        # Fortran
        mat_f = fc.scanHamPiece2c(interaction, isub, in1, in2, in2, dR, applyRotation=True)
        # PyOpenCL
        mat_p = ham.scanHamPiece2c(root, nz1, nz2, dR, applyRotation=True)
        
        if mat_p is not None:
            v_f = []
            v_p = []
            # s-px, s-py, s-pz, px-px, py-py, pz-pz
            idxs = [(0,3), (0,1), (0,2), (3,3), (1,1), (2,2)]
            for (i,j) in idxs:
                v_f.append(mat_f[i,j] if i < mat_f.shape[0] and j < mat_f.shape[1] else np.nan)
                v_p.append(mat_p[i,j] if i < mat_p.shape[0] and j < mat_p.shape[1] else np.nan)
            fs.append(v_f)
            ps.append(v_p)
        else:
            fs.append([np.nan]*6)
            ps.append([np.nan]*6)
             
    fs = np.array(fs)
    ps = np.array(ps)
    
    # Debug prints for angular components
    diff = np.abs(fs - ps)
    if np.isnan(diff).any():
        nan_locs = np.argwhere(np.isnan(diff))
        print(f"[2c_ang:{root}] WARNING NaN in diff at indices (theta_idx, component): {nan_locs.tolist()}")
        print("Fortran angular curves (raw):\n", fs)
        print("PyOCL angular curves (raw):\n", ps)
    print(f"[2c_ang:{root}] vmin/vmax Fortran: {fs.min(): .6e}/{fs.max(): .6e}  PyOCL: {ps.min(): .6e}/{ps.max(): .6e}  max|diff|={np.nanmax(diff): .6e}")
    labels = ["s-px","s-py","s-pz","px-px","py-py","pz-pz"]
    for i, lab in enumerate(labels):
        vmin_f, vmax_f = np.nanmin(fs[:, i]), np.nanmax(fs[:, i])
        vmin_p, vmax_p = np.nanmin(ps[:, i]), np.nanmax(ps[:, i])
        vmax_d = np.nanmax(diff[:, i])
        print(f"  [{lab}] Fortran vmin/vmax {vmin_f: .6e}/{vmax_f: .6e}  PyOCL vmin/vmax {vmin_p: .6e}/{vmax_p: .6e}  max|diff|={vmax_d: .6e}")
    if np.nanmax(diff) > 1e-4:
        print(f"[WARN] 2c angular {root} diff exceeds threshold; dumping curves")
        print("Fortran angular curves:\n", fs)
        print("PyOCL angular curves:\n", ps)

    plt.figure(figsize=(10, 6))
    colors = ['r','g','b','k','m','c']
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
    fc.initialize(atomType=np.array([nz1, nz2], dtype=np.int32), 
                  atomPos=np.array([[0,0,0],[0,0,1.0]], dtype=np.float64))
    
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
    
    # Angular Scan
    thetas = np.linspace(0, 2*np.pi, 50)
    plot_scan_2c_angular('overlap', nz1, nz2, 1, 0, 1.5, thetas)
    print("2-center verification complete.")
