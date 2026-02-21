import sys
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

sys.path.append("/home/prokophapala/git/FireCore-fitREQH/")
from pyBall import FitREQ_PN as fit

def load_atom_types(filepath):
    type_names = {}
    idx = 0
    try:
        with open(filepath, "r") as f:
            for line in f:
                line = line.strip()
                if line.startswith("#") or not line: continue
                parts = line.split()
                if len(parts) > 1:
                    type_names[idx] = parts[0]
                    idx += 1
    except Exception as e:
        print(f"Warning: could not load atom types from {filepath}: {e}")
    return type_names

def get_atom_labels(types, host, type_names, mode="chem"):
    # modes: "numeric" (e.g. 0T38), "chem" (e.g. O_3), "chem_no_unders" (e.g. O3)
    labels = []
    for i, (t, h) in enumerate(zip(types, host)):
        if mode == "numeric":
            lbl = f"{i}T{t}"
            if h >= 0:
                lbl += f"(E{h})"
            labels.append(lbl)
        else:
            chem_name = type_names.get(t, f"T{t}")
            if mode == "chem_no_unders":
                chem_name = chem_name.replace("_", "")
            labels.append(chem_name)
    return labels

def plot_2d_maps(angles, distances, E_dft, E_model, E_diff, min_a, min_d, Emin, save_path="plot_2d_maps.png"):
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    # Normalization based on Emin
    vmin, vmax = Emin, -Emin
    diff_scale = 0.25
    vmin_diff, vmax_diff = diff_scale * Emin, -diff_scale * Emin
    
    im0 = axes[0].imshow(E_dft, origin='lower', aspect='auto', cmap='coolwarm', vmin=vmin, vmax=vmax)
    axes[0].set_title('DFT Reference')
    
    im1 = axes[1].imshow(E_model, origin='lower', aspect='auto', cmap='coolwarm', vmin=vmin, vmax=vmax)
    axes[1].set_title('Model')
    
    im2 = axes[2].imshow(E_diff, origin='lower', aspect='auto', cmap='bwr', vmin=vmin_diff, vmax=vmax_diff)
    axes[2].set_title(f'Difference (Model - DFT) [scale={diff_scale}]')
    
    for ax in axes:
        # X is Angle
        xticks = np.linspace(0, len(angles)-1, 5, dtype=int)
        ax.set_xticks(xticks)
        ax.set_xticklabels([f"{angles[i]:.1f}" for i in xticks])
        
        # Y is Distance
        yticks = np.linspace(0, len(distances)-1, 5, dtype=int)
        ax.set_yticks(yticks)
        ax.set_yticklabels([f"{distances[i]:.2f}" for i in yticks])
        
        ax.set_xlabel("Angle [deg]")
        ax.set_ylabel("Distance [A]")
        
        # Plot crosshair (X=Angle_idx, Y=Dist_idx)
        ax.plot([min_a], [min_d], 'k+', markersize=15, markeredgewidth=2)
        
    plt.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.04)
    plt.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.04)
    plt.colorbar(im2, ax=axes[2], fraction=0.046, pad=0.04)
    
    plt.tight_layout()
    plt.savefig(save_path)
    plt.close()

def plot_1d_slices(angles, distances, E_dft, E_base, E_hcorr, E_sr, min_a, min_d, Emin, save_path="plot_1d_slices.png"):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    ang_val = angles[min_a]
    dist_val = distances[min_d]
    
    y_lims = [1.2 * Emin, -1.2 * Emin]
    dist_xlim = [min(distances), 6.0]
    
    # Distance Scan (fixed angle -> slice over dist)
    b_d = E_base[:, min_a]
    h_d = E_hcorr[:, min_a]
    s_d = E_sr[:, min_a]
    dft_d = E_dft[:, min_a]
    
    ax1.plot(distances, dft_d, 'k-', linewidth=2, label='DFT')
    ax1.plot(distances, b_d, color='gray', linestyle=':', linewidth=1, label='Baseline')
    ax1.fill_between(distances, b_d, b_d + h_d, color='cyan', alpha=0.5, label='+ Hcorr')
    ax1.fill_between(distances, b_d + h_d, b_d + h_d + s_d, color='orange', alpha=0.5, label='+ SR (E-pair)')
    
    ax1.axvline(dist_val, color='gray', linestyle=':')
    ax1.set_ylim(y_lims)
    ax1.set_xlim(dist_xlim)
    ax1.axhline(0, color='k', linewidth=0.5)
    ax1.set_title(f"Distance Scan @ Angle = {ang_val:.1f}")
    ax1.set_xlabel("Distance [A]")
    ax1.set_ylabel("Energy [kcal/mol]")
    ax1.legend()
    
    # Angular Scan (fixed distance -> slice over angle)
    b_a = E_base[min_d, :]
    h_a = E_hcorr[min_d, :]
    s_a = E_sr[min_d, :]
    dft_a = E_dft[min_d, :]
    
    ax2.plot(angles, dft_a, 'k-', linewidth=2, label='DFT')
    ax2.plot(angles, b_a, color='gray', linestyle=':', linewidth=1, label='Baseline')
    ax2.fill_between(angles, b_a, b_a + h_a, color='cyan', alpha=0.5)
    ax2.fill_between(angles, b_a + h_a, b_a + h_a + s_a, color='orange', alpha=0.5)
    
    ax2.axvline(ang_val, color='gray', linestyle=':')
    ax2.set_ylim(y_lims)
    ax2.axhline(0, color='k', linewidth=0.5)
    ax2.set_title(f"Angular Scan @ Distance = {dist_val:.2f}")
    ax2.set_xlabel("Angle [deg]")
    
    plt.tight_layout()
    plt.savefig(save_path)
    plt.close()

def plot_3d_geom(pos, n0, host, labels, save_path="plot_3d_geom.png"):
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    # Fragments
    f1_idx = list(range(n0))
    f2_idx = list(range(n0, len(pos)))
    
    f1 = pos[:n0]
    f2 = pos[n0:]
    
    ax.scatter(f1[:,0], f1[:,1], f1[:,2], c='blue', s=100, label='Fragment 1')
    ax.scatter(f2[:,0], f2[:,1], f2[:,2], c='red', s=100, label='Fragment 2')
    
    for i, lbl in enumerate(labels):
        ax.text(pos[i,0], pos[i,1], pos[i,2], lbl, color='black', fontsize=9, zorder=10)
    
    # Plot distances between dummy/atoms across fragments
    for i in f1_idx:
        for j in f2_idx:
            d = np.linalg.norm(pos[i] - pos[j])
            if d < 3.5: # only show relatively close ones to avoid clutter
                ax.plot([pos[i,0], pos[j,0]], [pos[i,1], pos[j,1]], [pos[i,2], pos[j,2]], 'k:', alpha=0.5)
                mid = (pos[i] + pos[j]) / 2.0
                ax.text(mid[0], mid[1], mid[2], f"{d:.2f}", color='green', fontsize=8)
                
    ax.set_title("Minimum Configuration Geometry")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.legend()
    plt.savefig(save_path)
    plt.close()

def plot_interaction_matrix(matrix, labels, title, save_path):
    fig, ax = plt.subplots(figsize=(10, 8))
    # We want a symmetric colormap around zero
    vmax = np.max(np.abs(matrix))
    if vmax < 1e-6: vmax = 1.0 # avoid zero division
    im = ax.imshow(matrix, cmap='bwr', vmin=-vmax, vmax=vmax)
    
    ax.set_xticks(np.arange(len(labels)))
    ax.set_yticks(np.arange(len(labels)))
    ax.set_xticklabels(labels, rotation=90)
    ax.set_yticklabels(labels)
    ax.set_title(title)
    
    plt.colorbar(im, ax=ax)
    plt.tight_layout()
    plt.savefig(save_path)
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Static visualization for FitREQ_PN surfaces")
    parser.add_argument("--input", type=str, default="wb97m-split/H2O-A1_H2O-D1-z.xyz", help="Input XYZ path (relative to script dir or absolute)")
    parser.add_argument("--dof_selection", type=str, default="data/dofSelection_run.dat", help="DOF selection file (relative to script dir or absolute)")
    parser.add_argument("--atom_types", type=str, default="data/AtomTypes.dat", help="AtomTypes.dat path")
    parser.add_argument("--element_types", type=str, default="data/ElementTypes.dat", help="ElementTypes.dat path")
    parser.add_argument("--label_mode", type=str, default="chem_no_unders", choices=["chem_no_unders", "chem", "numeric"], help="Label style")
    parser.add_argument("--lepairs", type=float, default=0.8, help="Epair distance from host")
    parser.add_argument("--kmorse", type=float, default=1.8, help="Morse curvature")
    args = parser.parse_args()

    here = os.path.dirname(__file__)
    # Resolve paths relative to script dir if not absolute
    def resolve(p):
        return p if os.path.isabs(p) else os.path.join(here, p)

    fname = resolve(args.input)
    dof_selection = resolve(args.dof_selection)
    atom_types_file = resolve(args.atom_types)
    element_types_file = resolve(args.element_types)
    label_mode = args.label_mode

    type_names = load_atom_types(atom_types_file)
    
    fit.setVerbosity(verbosity=1)
    fit.setModel(ivdW=4, iCoul=1, iHbond=3, Epairs=1, iEpairs=2, kMorse=args.kmorse, Lepairs=args.lepairs, bPN=True)
    fit.loadTypes(fEtypes=element_types_file,  fAtypes=atom_types_file)
    fit.loadDOFSelection(fname=dof_selection)
    
    n_total = fit.loadXYZ(fname, bAddEpairs=True, bOutXYZ=False, bEvalOnlyCorrections=False, bAppend=False)
    fit.getBuffs()
    
    Erefs, x0s = fit.read_xyz_data(fname)
    Eerr, Es_model, _ = fit.getEs(bDOFtoTypes=True, bEs=True)
    _, Es_Coul, Es_vdW, Es_Hcorr, Es_Epairs = fit.getEs_components()
    
    Es_base = Es_Coul + Es_vdW
    
    Gref, seq, axis, distances, angles = fit.parse_xyz_mapping(fname)
    
    shape = Gref.shape  # this is (len(distances), len(angles)) if correctly parsed
    Gmodel = np.full(shape, np.nan)
    Gbase = np.full(shape, np.nan)
    Ghcorr = np.full(shape, np.nan)
    Gsr = np.full(shape, np.nan)
    
    nmap = min(len(Es_model), len(seq))
    for k in range(nmap):
        idist, iang = seq[k]
        Gmodel[idist, iang] = Es_model[k]
        Gbase[idist, iang] = Es_base[k]
        Ghcorr[idist, iang] = Es_Hcorr[k]
        Gsr[idist, iang] = Es_Epairs[k]
        
    Gdiff = Gmodel - Gref
    
    # Find minimum in DFT (Gref)
    valid = ~np.isnan(Gref)
    min_idx = np.unravel_index(np.nanargmin(Gref), Gref.shape)
    min_d, min_a = min_idx
    Emin = Gref[min_d, min_a]
    
    print("="*60)
    print(f"Global min at angle_idx={min_a} (angle={angles[min_a]:.1f}), dist_idx={min_d} (dist={distances[min_d]:.2f})")
    print(f"Min Energy: {Emin}")
    
    # Find corresponding sample index
    isamp = -1
    for k in range(nmap):
        if seq[k] == (min_d, min_a):
            isamp = k
            break
            
    if isamp >= 0:
        pos, types, Qs, host, n0 = fit.getSampleGeom(isamp)
        labels = get_atom_labels(types, host, type_names, mode=label_mode)
        
        print("="*60)
        print(f"Fragment Separation: n0 = {n0}, total particles = {len(pos)}")
        print(f"Fragment 1 size: {n0}, Fragment 2 size: {len(pos) - n0}")
        print("Particle Labels & Epair Hosts:")
        for i, lbl in enumerate(labels):
            print(f"  {i}: {lbl} (host={host[i]})")
            
        plot_3d_geom(pos, n0, host, labels, save_path="plot_3d_geom.png")
        
        # Evaluate Sample Pairs
        pair_out = fit.evalSamplePairs(isamp, len(pos))
        
        components = [("vdW_Morse", 0), ("Coulomb", 1), ("Hcorr", 2), ("SR_Epairs", 3)]
        
        print("="*60)
        for name, comp_idx in components:
            mat = pair_out[:, :, comp_idx]
            plot_interaction_matrix(mat, labels, f"{name} Interaction Matrix", save_path=f"plot_mat_{name}.png")
            print(f"\n--- {name} Interaction Matrix ---")
            
            # Print header
            header = "     " + "".join([f"{l:>8}" for l in labels])
            print(header)
            for i in range(len(labels)):
                row_str = f"{labels[i]:>4} " + "".join([f"{mat[i,j]:8.4f}" for j in range(len(labels))])
                print(row_str)
                
    else:
        print("Could not find sample for minimum!")
        
    # Gref is [dist, angle]. We want X=angle, Y=dist.
    plot_2d_maps(angles, distances, Gref, Gmodel, Gdiff, min_a, min_d, Emin, save_path="plot_2d_maps.png")
    plot_1d_slices(angles, distances, Gref, Gbase, Ghcorr, Gsr, min_a, min_d, Emin, save_path="plot_1d_slices.png")

if __name__ == "__main__":
    main()
