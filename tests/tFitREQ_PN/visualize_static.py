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
    # Shift reference to its minimum (match plot_compare_combined)
    Eref_raw = np.array(E_dft, copy=True)
    ref_min = np.nanmin(Eref_raw)
    Eref = Eref_raw - ref_min
    Emod = np.array(E_model, copy=True) if E_model is not None else None
    if Emod is not None:
        Emod = Emod - ref_min
    Ediff = np.array(E_diff, copy=True)

    # Shared vmin/vmax from reference
    vmin = float(np.nanmin(Eref)) if np.any(np.isfinite(Eref)) else None
    if vmin is not None and vmin > 0: vmin = 0.0
    vmax = -vmin if vmin is not None else None

    # Difference clamp (symmetric, small span)
    E_MARGIN = 1.0  # kcal/ev units consistent with inputs (kcal here)

    # Zero line from raw DFT: angle-dependent E=0 crossing mapped to pixel indices
    zero_r_pixels = []
    r_vec = np.asarray(distances, dtype=float)
    for i_ang in range(len(angles)):
        row = Eref_raw[:, i_ang]
        r0 = np.nan
        for j in range(len(r_vec)-1):
            if not (np.isfinite(row[j]) and np.isfinite(row[j+1])):
                continue
            if row[j] > 0 and row[j+1] <= 0:
                r0 = r_vec[j] - row[j]*(r_vec[j+1]-r_vec[j])/(row[j+1]-row[j])
                break
        if np.isnan(r0):
            zero_r_pixels.append(np.nan)
        else:
            zero_r_pixels.append(int(np.argmin(np.abs(r_vec - r0))))

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    im0 = axes[0].imshow(Eref, origin='lower', aspect='auto', cmap='bwr', vmin=vmin, vmax=vmax)
    axes[0].set_title('DFT Reference')

    im1 = axes[1].imshow(Emod, origin='lower', aspect='auto', cmap='bwr', vmin=vmin, vmax=vmax)
    axes[1].set_title('Model')
    axes[1].plot(range(len(angles)), zero_r_pixels, linestyle=':', color='k', alpha=0.7, label='E_DFT=0')
    axes[1].legend(fontsize=8)

    im2 = axes[2].imshow(Ediff, origin='lower', aspect='auto', cmap='bwr', vmin=-E_MARGIN, vmax=E_MARGIN)
    axes[2].set_title('Difference (Model - DFT)')
    axes[2].plot(range(len(angles)), zero_r_pixels, linestyle=':', color='k', alpha=0.7)

    for ax in axes:
        xticks = np.linspace(0, len(angles)-1, 5, dtype=int)
        ax.set_xticks(xticks)
        ax.set_xticklabels([f"{angles[i]:.1f}" for i in xticks])

        yticks = np.linspace(0, len(distances)-1, 5, dtype=int)
        ax.set_yticks(yticks)
        ax.set_yticklabels([f"{distances[i]:.2f}" for i in yticks])

        ax.set_xlabel("Angle [deg]")
        ax.set_ylabel("Distance [A]")

        # Crosshair at global minimum
        ax.plot([min_a], [min_d], 'k+', markersize=15, markeredgewidth=2)

    plt.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.04)
    plt.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.04)
    plt.colorbar(im2, ax=axes[2], fraction=0.046, pad=0.04)

    plt.tight_layout()
    plt.savefig(save_path)
    plt.close()

def plot_min_lines_vs_angle(Epanel, distances):
    """For each angle column, find Emin and corresponding Rmin (distance)."""
    n_angles = Epanel.shape[1]
    Emin = np.full(n_angles, np.nan)
    Rmin = np.full(n_angles, np.nan)
    for i in range(n_angles):
        col = Epanel[:, i]
        if np.all(np.isnan(col)):
            continue
        j = np.nanargmin(col)
        Emin[i] = col[j]
        Rmin[i] = distances[j]
    return Rmin, Emin


def model_min_with_components(Gmodel, Gbase, Ghcorr, Gsr, distances):
    """Per angle: take argmin of total model, then sample base/H/SR at that same distance.
    Returns (Rmin_model, E_model, E_base, E_base_H, E_base_H_SR).
    """
    n_angles = Gmodel.shape[1]
    rmin = np.full(n_angles, np.nan)
    e_tot = np.full(n_angles, np.nan)
    e_base = np.full(n_angles, np.nan)
    e_base_h = np.full(n_angles, np.nan)
    e_full = np.full(n_angles, np.nan)
    for i in range(n_angles):
        col = Gmodel[:, i]
        if np.all(np.isnan(col)):
            continue
        j = np.nanargmin(col)
        rmin[i] = distances[j]
        e_tot[i] = Gmodel[j, i]
        e_base[i] = Gbase[j, i]
        e_base_h[i] = Gbase[j, i] + Ghcorr[j, i]
        e_full[i] = Gbase[j, i] + Ghcorr[j, i] + Gsr[j, i]
    return rmin, e_tot, e_base, e_base_h, e_full


def plot_1d_slices(angles, distances, E_dft, E_base, E_hcorr, E_sr, min_a, min_d, Emin, mode="rigid", save_path="plot_1d_slices.png"):
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
    
    if mode == "rigid":
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
    else:
        # Emin(angle) using model total argmin; sample base/H/SR at the same distance
        r_ref, e_ref = plot_min_lines_vs_angle(E_dft, distances)
        r_model, e_model, e_base, e_base_h, e_full = model_min_with_components(E_base + E_hcorr + E_sr, E_base, E_hcorr, E_sr, distances)

        ax2.plot(angles, e_ref, 'k-', linewidth=2, label='DFT (min per angle)')
        ax2.plot(angles, e_full, color='orange', linewidth=1.5, label='Model Emin')
        ax2.plot(angles, e_base, color='gray', linestyle=':', linewidth=1,  label='Base @ model Rmin')
        ax2.fill_between(angles, e_base, e_base_h, color='cyan',   alpha=0.5, label='+ Hcorr @ model Rmin')
        ax2.fill_between(angles, e_base_h, e_full, color='orange', alpha=0.3, label='+ SR @ model Rmin')
        ax2.axhline(0, color='k', linewidth=0.5)
        ax2.set_title("Emin vs Angle (model argmin per angle)")
        ax2.set_xlabel("Angle [deg]")
        ax2.set_ylabel("Energy [kcal/mol]")

        # Rmin overlay
        ax2b = ax2.twinx()
        ax2b.plot(angles, r_ref, 'k--', linewidth=1, label='Rmin (DFT)')
        ax2b.plot(angles, r_model, color='orange', linestyle='--', linewidth=1, label='Rmin (model)')
        ax2b.set_ylabel("Rmin [A]")
        ax2b.set_ylim(1.5, 3.0)
        # Combine legends
        lines, labels   = ax2.get_legend_handles_labels()
        lines2, labels2 = ax2b.get_legend_handles_labels()
        ax2b.legend(lines + lines2, labels + labels2, loc='upper right')

    
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
    parser.add_argument("--input",         type=str,   default="wb97m-split/H2O-A1_H2O-D1-y.xyz", help="Input XYZ path (relative to script dir or absolute)")
    #parser.add_argument("--dof_selection",type=str,   default="data/dofSelection_run.dat", help="DOF selection file (relative to script dir or absolute)")
    parser.add_argument("--dof_selection", type=str,   default="data/dofSelection_H2O.dat", help="DOF selection file (relative to script dir or absolute)")
    parser.add_argument("--atom_types",    type=str,   default="data/AtomTypes.dat", help="AtomTypes.dat path")
    parser.add_argument("--element_types", type=str,   default="data/ElementTypes.dat", help="ElementTypes.dat path")
    parser.add_argument("--label_mode",    type=str,   default="chem_no_unders", choices=["chem_no_unders", "chem", "numeric"], help="Label style")
    parser.add_argument("--lepairs",       type=float, default=0.8, help="Epair distance from host")
    parser.add_argument("--kmorse",        type=float, default=1.8, help="Morse curvature")
    parser.add_argument("--slice_mode",    type=str,   default="emin", choices=["rigid", "emin"], help="Angular slice mode: rigid uses fixed distance, emin plots per-angle minima Rmin/Emin")
    parser.add_argument("--init_relax",    type=int,   default=0, help="Do a single PN step with dt=0 to apply DOF xstart before evaluation")
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

    if args.init_relax:
        # One step, dt=0 ensures only DOFsToTypes/regularization applied (no movement)
        fit.run_PN(ialg=0, iparallel=0, nstep=1, Fmax=0.0, dt=0.0001, max_step=0.0, damping=0.0)

    Erefs, x0s = fit.read_xyz_data(fname)
    Eerr, Es_model, _ = fit.getEs(bDOFtoTypes=True, bEs=True)
    _, Es_Coul, Es_vdW, Es_Hcorr, Es_Epairs = fit.getEs_components()
    
    Es_base = Es_Coul + Es_vdW
    
    Gref, seq, axis, distances, angles = fit.parse_xyz_mapping(fname)
    
    shape  = Gref.shape  # this is (len(distances), len(angles)) if correctly parsed
    Gmodel = np.full(shape, np.nan)
    Gbase  = np.full(shape, np.nan)
    Ghcorr = np.full(shape, np.nan)
    Gsr   = np.full(shape, np.nan)
    
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
    plot_1d_slices(angles, distances, Gref, Gbase, Ghcorr, Gsr, min_a, min_d, Emin, mode=args.slice_mode, save_path="plot_1d_slices.png")

if __name__ == "__main__":
    main()
