import numpy as np
import os
import sys
import itertools
import matplotlib.pyplot as plt

# We will dynamically bind 'fit' based on initialization
fit = None

def init(use_pn=False):
    global fit
    if use_pn:
        import pyBall.FitREQ_PN as fit_module
    else:
        import pyBall.FitREQ as fit_module
    fit = fit_module
    return fit

def parse_dof_file(fname):
    dofs = []
    hyper = {}
    with open(fname, 'r') as f:
        for line in f:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith('#HYPER'):
                parts = stripped.split()[1:]
                for p in parts:
                    if '=' in p:
                        k, v = p.split('=', 1)
                        hyper[k] = float(v)
                continue
            if stripped.startswith('#'):
                continue
            parts = stripped.split()
            if len(parts) >= 10:
                dofs.append({
                    'typename': parts[0],
                    'comp': int(parts[1]),
                    'Min': float(parts[2]),
                    'Max': float(parts[3]),
                    'xlo': float(parts[4]),
                    'xhi': float(parts[5]),
                    'Klo': float(parts[6]),
                    'Khi': float(parts[7]),
                    'K0': float(parts[8]),
                    'xstart': float(parts[9]),
                    'invMass': float(parts[10]) if len(parts)>10 else 1.0
                })
    return dofs, hyper

def write_dof_file(fname, dofs, hyper=None):
    with open(fname, 'w') as f:
        if hyper:
            kv = " ".join(f"{k}={v}" for k, v in hyper.items())
            f.write(f"#HYPER {kv}\n")
        f.write("# typename comp Min Max xlo xhi Klo Khi K0 xstart invMass\n")
        for d in dofs:
            f.write(f"{d['typename']} {d['comp']} {d['Min']} {d['Max']} {d['xlo']} {d['xhi']} {d['Klo']} {d['Khi']} {d['K0']} {d['xstart']} {d['invMass']}\n")

def get_dof_values(dofs):
    return np.array([d['xstart'] for d in dofs])

def set_dof_values(dofs, values):
    for i, d in enumerate(dofs):
        d['xstart'] = values[i]

TYPES_LOADED = False

def setup_fitreq(dof_file, xyz_path, kMorse=1.7, Lepairs=0.9, hScale=0.0, softClamp=True, sc_start=4.0, sc_max=6.0, imodel=7, Regularize=1):
    global TYPES_LOADED
    
    if hasattr(fit, 'setup'):
        fit.setVerbosity(0)
        fit.setup(imodel=imodel, EvalJ=1, WriteJ=1, Regularize=Regularize, SaveJustElementXYZ=-1)
        if softClamp:
            fit.setGlobalParams(kMorse=kMorse, Lepairs=Lepairs, hScale=hScale, softClamp_start=sc_start, softClamp_max=sc_max)
        else:
            fit.setGlobalParams(kMorse=kMorse, Lepairs=Lepairs, hScale=hScale)
    else:
        fit.setVerbosity(verbosity=0)
        fit.setModel(ivdW=4, iCoul=1, iHbond=3, Epairs=1, iEpairs=2, kMorse=kMorse, Lepairs=Lepairs, bPN=True)
        
    if not TYPES_LOADED:
        fit.loadTypes()
        TYPES_LOADED = True
    fit.loadDOFSelection(fname=dof_file)
    fit.loadXYZ(xyz_path, bAddEpairs=True)
    fit.getBuffs()
    return fit.loadDOFnames(dof_file)

def evaluate_current(bFs=False):
    if hasattr(fit, 'getEs') and 'bOmp' in fit.getEs.__code__.co_varnames:
        Eerr, Es, Fs = fit.getEs(bOmp=False, bDOFtoTypes=True, bEs=True, bFs=bFs)
    else:
        Eerr, Es, Fs = fit.getEs(bDOFtoTypes=True, bEs=True, bFs=bFs)
    return Eerr, Es, Fs

def relax_tethered(dofs, kMorse=1.7, Lepairs=0.9, nstep=100, Fmax=1e-8, dt=0.05, damping=0.0):
    tmp_dof = "_tmp_dofs_tether.dat"
    write_dof_file(tmp_dof, dofs)
    fit.loadDOFSelection(tmp_dof)
    fit.getBuffs()
    
    if hasattr(fit, 'run'):
        fit.run(iparallel=0, ialg=1, nstep=nstep, Fmax=Fmax, dt=dt, damping=damping, max_step=0.1, bClamp=True)
    else:
        fit.run_PN(ialg=1, iparallel=0, nstep=nstep, Fmax=Fmax, dt=dt, max_step=0.1, damping=damping)
        
    fit.getBuffs()
    
    vals = np.array(fit.DOFs)
    set_dof_values(dofs, vals)
    return vals

def compute_fd_hessian(dofs, h=1e-4):
    n = len(dofs)
    H = np.zeros((n, n))
    base_vals = get_dof_values(dofs)
    
    for i in range(n):
        v_plus = base_vals.copy()
        v_plus[i] += h
        for j in range(n): fit.DOFs[j] = v_plus[j]
        _, _, F_plus = evaluate_current(bFs=True)
        
        v_minus = base_vals.copy()
        v_minus[i] -= h
        for j in range(n): fit.DOFs[j] = v_minus[j]
        _, _, F_minus = evaluate_current(bFs=True)
        
        H[i, :] = (F_plus - F_minus) / (2 * h)
        
    for j in range(n): fit.DOFs[j] = base_vals[j]
    
    H = (H + H.T) / 2.0
    eigvals, eigvecs = np.linalg.eigh(H)
    return H, eigvals, eigvecs

def run_grid_search(xyz_path, base_dof_file, out_dir, kMorse_vals, Lepairs_vals, weight_alpha_vals, hScale_vals):
    os.makedirs(out_dir, exist_ok=True)
    Gref, seq, axis, distances, angles = fit.parse_xyz_mapping(xyz_path)
    Erefs, x0s = fit.read_xyz_data(xyz_path)
    base_dofs, _ = parse_dof_file(base_dof_file)
    dof_names = setup_fitreq(base_dof_file, xyz_path, imodel=7, Regularize=1)
    
    if hasattr(fit, 'setVerbosity') and 'idebug' in fit.setVerbosity.__code__.co_varnames:
        fit.setVerbosity(1, 0, 0, 0, 0, 0, 0)
    else:
        fit.setVerbosity(1)
        
    results = []
    run_idx = 0
    for kM, Lep, w_alpha, hScale in itertools.product(kMorse_vals, Lepairs_vals, weight_alpha_vals, hScale_vals):
        print(f"Run {run_idx}: kMorse={kM}, Lepairs={Lep}, weight_alpha={w_alpha}, hScale={hScale}")
        dofs = [dict(d) for d in base_dofs]
        write_dof_file("tmp_grid_dofs.dat", dofs)
        fit.loadDOFSelection("tmp_grid_dofs.dat")
        
        if hasattr(fit, 'setGlobalParams'):
            fit.setGlobalParams(kMorse=kM, Lepairs=Lep, softClamp_start=4.0, softClamp_max=6.0, hScale=hScale)
        else:
            fit.setModel(ivdW=4, iCoul=1, iHbond=3, Epairs=1, iEpairs=2, kMorse=kM, Lepairs=Lep, bPN=True)
            
        weight_func = lambda E: fit.exp_weight_func(E, a=1.0, alpha=w_alpha)
        if hasattr(fit, 'split_and_weight_curves'):
            weights0, lens = fit.split_and_weight_curves(Erefs, x0s, n_before_min=100, weight_func=weight_func, EminMin=-0.02)
            fit.setWeights(weights0)
        
        if hasattr(fit, 'setVerbosity') and 'idebug' in fit.setVerbosity.__code__.co_varnames:
            fit.setVerbosity(0, 0, 0, 0, 0, 0, 0)
        else:
            fit.setVerbosity(0)
            
        if hasattr(fit, 'run'):
            fit.run(iparallel=0, ialg=1, nstep=150, Fmax=1e-8, dt=0.05, damping=0.0, max_step=0.1, bClamp=True)
        else:
            fit.run_PN(ialg=1, iparallel=0, nstep=150, Fmax=1e-8, dt=0.05, max_step=0.1, damping=0.0)
            
        fit.getBuffs()
        final_dofs = np.array(fit.DOFs)
        Eerr, Es, Fs = evaluate_current(bFs=False)
        
        if hasattr(fit, 'compute_model_grid'):
            run_params = dict(nstep=150, Fmax=1e-8, dt=0.05, max_step=0.1, damping=0.0)
            Gmodel = fit.compute_model_grid(xyz_path, seq, Gref.shape, do_fit=False, bAddEpairs=True, run_params=run_params, bOutXYZ=False)
            
            prefix = f"{out_dir}/run_{run_idx}_kM{kM}_L{Lep}_wa{w_alpha}_hS{hScale}"
            title = f"kM={kM} L={Lep} wa={w_alpha} hS={hScale} Err={Eerr:.4e}"
            params_lines = [f"kMorse={kM:.3f}", f"Lepairs={Lep:.3f}", f"hScale={hScale:.3f}", f"w_alpha={w_alpha:.3f}"]
            params_lines += [f"{name}={val:.4f}" for name, val in zip(dof_names, final_dofs)]
            params_text = "\n".join(params_lines)
            
            fit.plot_compare_combined(Gref, Gmodel, angles, distances, title, save_path=f"{prefix}.png", kcal=True, params_text=params_text)
        else:
            # We don't have plot_compare_combined in PN right now
            prefix = f"{out_dir}/run_{run_idx}_kM{kM}_L{Lep}_wa{w_alpha}_hS{hScale}"
        
        for i, d in enumerate(dofs):
            d['xstart'] = final_dofs[i]
        
        hyper_meta = {'kMorse': kM, 'Lepairs': Lep, 'weight_alpha': w_alpha, 'hScale': hScale}
        write_dof_file(f"{prefix}.dat", dofs, hyper=hyper_meta)
        
        results.append({'kMorse': kM, 'Lepairs': Lep, 'weight_alpha': w_alpha, 'error': Eerr, 'dofs': final_dofs})
        run_idx += 1
        
    np.savez(f"{out_dir}/ensemble_data.npz", 
             kMorse=np.array([r['kMorse'] for r in results]),
             Lepairs=np.array([r['Lepairs'] for r in results]),
             weight_alpha=np.array([r['weight_alpha'] for r in results]),
             error=np.array([r['error'] for r in results]),
             dofs=np.array([r['dofs'] for r in results]),
             dof_names=dof_names)
             
    print("Grid search complete.")

def plot_string_results(alphas, energies_init, energies_relax, chain_relaxed, gradients, param_names, out_prefix):
    fig, axes = plt.subplots(3, 1, figsize=(10, 12), sharex=True)
    
    # --- Plot 1: The Energy Barrier (Loss Landscape) ---
    ax = axes[0]
    ax.plot(alphas, energies_init, 'k--', label='Linear Interpolation (Rigid)')
    if energies_relax is not None and len(energies_relax) > 0:
        ax.plot(alphas, energies_relax, 'b-o', label='Restrained Relaxation')
    ax.set_ylabel("Fitting Error (Objective)")
    ax.set_title("Energy Barrier Between Solutions")
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # --- Plot 2: Parameter Evolution ---
    ax = axes[1]
    # Scale parameters that are too large or too small for [-1, 1] plot
    dofs_plot = np.copy(chain_relaxed)
    scales = np.ones(len(param_names))
    for j in range(len(param_names)):
        max_abs = np.max(np.abs(dofs_plot[:, j]))
        if max_abs > 1.0:
            scales[j] = 0.1
        elif max_abs > 0.0 and max_abs < 0.1:
            scales[j] = 10.0
        dofs_plot[:, j] *= scales[j]
        
    for j in range(len(param_names)):
        label = param_names[j]
        if scales[j] != 1.0:
            label += f" (*{scales[j]:.1f})"
        ax.plot(alphas, dofs_plot[:, j], '.-', label=label)
        
    ax.set_ylabel("Parameter Value")
    ax.set_title("Parameter Evolution along Relaxed Path")
    ax.legend(bbox_to_anchor=(1.04, 1), loc="upper left", fontsize='small', ncol=2)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(-1.1, 1.1)
    
    # --- Plot 3: Variational Forces (Gradients) ---
    ax = axes[2]
    for i, name in enumerate(param_names):
        ax.plot(alphas, gradients[:, i], '.-', label=f"dE/d({name})")
    ax.axhline(0, color='k', linestyle=':')
    ax.set_ylabel("Gradient (Force)")
    ax.set_xlabel("Interpolation Progress (α)")
    ax.set_title("Gradients along the path")
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f"{out_prefix}_string.png", dpi=150)
    plt.close()

def plot_hessian_analysis(H, eigvals, eigvecs, param_names, alpha, out_prefix):
    n = len(param_names)
    
    # 1. Softest Mode Composition
    soft_mode = eigvecs[:, 0]
    plt.figure(figsize=(8, 4))
    plt.bar(param_names, soft_mode)
    plt.title(f"Lowest Frequency Mode (alpha={alpha:.2f})\nEigval = {eigvals[0]:.2e}")
    plt.xticks(rotation=45, ha='right')
    plt.ylabel("Component")
    plt.tight_layout()
    plt.savefig(f"{out_prefix}_hessian_softmode_alpha_{alpha:.2f}.png")
    plt.close()
    
    # 2. Parameter Correlation Matrix
    # pseudo-inverse of Hessian ~ covariance matrix of parameters around minimum
    H_inv = np.linalg.pinv(H, rcond=1e-5)
    
    # Convert covariance to correlation matrix: C_ij / sqrt(C_ii * C_jj)
    # add small epsilon to avoid div by zero
    std_dev = np.sqrt(np.diag(H_inv)) + 1e-12
    corr = H_inv / np.outer(std_dev, std_dev)
    
    plt.figure(figsize=(8, 6))
    plt.imshow(corr, cmap='bwr', vmin=-1, vmax=1)
    plt.colorbar(label='Correlation Coefficient')
    plt.xticks(range(n), param_names, rotation=45, ha='right')
    plt.yticks(range(n), param_names)
    plt.title(f"Parameter Correlation Heatmap (alpha={alpha:.2f})")
    plt.tight_layout()
    plt.savefig(f"{out_prefix}_hessian_corr_alpha_{alpha:.2f}.png")
    plt.close()

def run_string_scan(xyz_path, dofA, dofB, n_points, kmorse, lepairs, springK, nstep, out_prefix, do_relax, do_hessian):
    import copy
    
    dofsA, hyperA = parse_dof_file(dofA)
    dofsB, hyperB = parse_dof_file(dofB)
    valsA = get_dof_values(dofsA)
    valsB = get_dof_values(dofsB)
    
    def hval(h, key, default):
        return h.get(key, default)
        
    kA, kB = hval(hyperA, 'kMorse', kmorse), hval(hyperB, 'kMorse', kmorse)
    lA, lB = hval(hyperA, 'Lepairs', lepairs), hval(hyperB, 'Lepairs', lepairs)
    wA, wB = hval(hyperA, 'weight_alpha', 1.0), hval(hyperB, 'weight_alpha', 1.0)
    hSA, hSB = hval(hyperA, 'hScale', 0.0), hval(hyperB, 'hScale', 0.0)
    
    dof_names = setup_fitreq(dofA, xyz_path, kMorse=kA, Lepairs=lA, hScale=hSA)
    Erefs, x0s = fit.read_xyz_data(xyz_path)
    
    alphas = np.linspace(0, 1, n_points)
    chain_initial = [(1 - a) * valsA + a * valsB for a in alphas]
    
    chain_relaxed = []
    energies_init = []
    energies_relax = []
    gradients = []
    
    kM_path, Lp_path, wa_path, hS_path = [], [], [], []
    
    dof_tmp = copy.deepcopy(dofsA)

    print(f"Running String Scan ({n_points} points). Relax: {do_relax}")
    for i, p_target in enumerate(chain_initial):
        a = alphas[i]
        print(f"Step {i+1}/{n_points} (alpha={a:.2f})")
        
        kM = (1-a)*kA + a*kB
        Lp = (1-a)*lA + a*lB
        wa = (1-a)*wA + a*wB
        hS = (1-a)*hSA + a*hSB
        
        kM_path.append(kM); Lp_path.append(Lp); wa_path.append(wa); hS_path.append(hS)
        
        if hasattr(fit, 'setGlobalParams'):
            fit.setGlobalParams(kMorse=kM, Lepairs=Lp, hScale=hS, softClamp_start=4.0, softClamp_max=6.0)
        else:
            fit.setModel(ivdW=4, iCoul=1, iHbond=3, Epairs=1, iEpairs=2, kMorse=kM, Lepairs=Lp, bPN=True)
            
        weight_func = lambda E: fit.exp_weight_func(E, a=1.0, alpha=wa)
        if hasattr(fit, 'split_and_weight_curves'):
            weights0, lens = fit.split_and_weight_curves(Erefs, x0s, n_before_min=100, weight_func=weight_func, EminMin=-0.02)
            fit.setWeights(weights0)
        
        set_dof_values(dof_tmp, p_target)
        for j in range(len(p_target)): fit.DOFs[j] = p_target[j]
        E_init, _, grads_init = evaluate_current(bFs=True)
        energies_init.append(E_init)
        
        if do_relax:
            for d in dof_tmp: d['K0'] = springK
            p_final = relax_tethered(dof_tmp, kMorse=kM, Lepairs=Lp, nstep=nstep)
            E_final, _, grads_final = evaluate_current(bFs=True)
            
            chain_relaxed.append(p_final)
            energies_relax.append(E_final)
            gradients.append(grads_final)
        else:
            chain_relaxed.append(p_target)
            gradients.append(grads_init)
            
    chain_relaxed = np.array(chain_relaxed)
    energies_init = np.array(energies_init)
    energies_relax = np.array(energies_relax) if do_relax else None
    gradients = np.array(gradients)
    
    np.savez(f"{out_prefix}_trajectory.npz", 
             alphas=alphas, 
             chain_initial=np.array(chain_initial),
             chain_relaxed=chain_relaxed, 
             energies_init=energies_init, 
             energies_relax=energies_relax if do_relax else np.array([]),
             param_names=dof_names,
             kMorse_path=np.array(kM_path),
             Lepairs_path=np.array(Lp_path),
             weight_alpha_path=np.array(wa_path),
             hScale_path=np.array(hS_path))
             
    plot_string_results(alphas, energies_init, energies_relax, chain_relaxed, gradients, dof_names, out_prefix)
    
    if do_hessian:
        print("\nComputing Hessians...")
        for h_alpha in [0.0, 0.5, 1.0]:
            idx = np.argmin(np.abs(alphas - h_alpha))
            a = alphas[idx]
            print(f"Hessian at alpha ~ {a:.2f}")
            
            if hasattr(fit, 'setGlobalParams'):
                fit.setGlobalParams(kMorse=kM_path[idx], Lepairs=Lp_path[idx], hScale=hS_path[idx], softClamp_start=4.0, softClamp_max=6.0)
            else:
                fit.setModel(ivdW=4, iCoul=1, iHbond=3, Epairs=1, iEpairs=2, kMorse=kM_path[idx], Lepairs=Lp_path[idx], bPN=True)
                
            weight_func = lambda E: fit.exp_weight_func(E, a=1.0, alpha=wa_path[idx])
            if hasattr(fit, 'split_and_weight_curves'):
                weights0, lens = fit.split_and_weight_curves(Erefs, x0s, n_before_min=100, weight_func=weight_func, EminMin=-0.02)
                fit.setWeights(weights0)
            
            set_dof_values(dof_tmp, chain_relaxed[idx])
            H, eigvals, eigvecs = compute_fd_hessian(dof_tmp)
            
            print(f"  Eigenvalues: {eigvals}")
            
            plt.figure(figsize=(6,4))
            plt.plot(eigvals, 'o-')
            plt.axhline(0, color='k', linestyle='--')
            plt.title(f"Hessian Eigenvalues (alpha={alphas[idx]:.2f})")
            plt.ylabel("Eigenvalue")
            plt.xlabel("Mode Index")
            plt.grid(True, alpha=0.3)
            plt.savefig(f"{out_prefix}_hessian_eigs_alpha_{alphas[idx]:.2f}.png")
            plt.close()
            
            plot_hessian_analysis(H, eigvals, eigvecs, dof_names, a, out_prefix)

def run_dr(traj_file, out_prefix):
    try:
        import pacmap
        HAS_PACMAP = True
    except ImportError:
        HAS_PACMAP = False
        
    try:
        import umap
        HAS_UMAP = True
    except ImportError:
        HAS_UMAP = False
        
    from sklearn.decomposition import PCA

    if not os.path.exists(traj_file):
        raise FileNotFoundError(f"File {traj_file} not found.")
        
    data = np.load(traj_file)

    if 'chain_initial' in data and 'chain_relaxed' in data:
        is_path = True
        chain_initial = data['chain_initial']
        chain_relaxed = data['chain_relaxed']
        energies_init = data['energies_init']
        energies_relax = data['energies_relax']
        X = np.vstack([chain_initial, chain_relaxed])
        E = np.concatenate([energies_init, energies_relax]) if len(energies_relax) > 0 else energies_init
        n_init = len(chain_initial)
    elif 'dofs' in data and 'error' in data:
        is_path = False
        X = np.asarray(data['dofs'])
        E = np.asarray(data['error'])
        chain_initial = X
        chain_relaxed = np.zeros((0, X.shape[1]))
        n_init = len(chain_initial)
    else:
        raise ValueError(f"Unsupported NPZ format: missing chain_* or dofs/error in {traj_file}")
    
    if len(E) != len(X):
        if len(E) < len(X):
            pad_val = E[-1] if len(E) > 0 else 0.0
            E = np.concatenate([E, np.full(len(X) - len(E), pad_val)])
        else:
            E = E[:len(X)]
            
    E_log = np.log(E - E.min() + 1e-6)
    
    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X)
    
    plt.figure(figsize=(8,6))
    sc = plt.scatter(X_pca[:, 0], X_pca[:, 1], c=E_log, cmap='viridis', s=50)
    if is_path and len(X) >= 2:
        plt.plot(X_pca[:n_init, 0], X_pca[:n_init, 1], 'k--', alpha=0.5, label='Initial Path')
        if len(chain_relaxed) > 0:
            plt.plot(X_pca[n_init:, 0], X_pca[n_init:, 1], 'b-', alpha=0.5, label='Relaxed Path')
    plt.colorbar(sc, label="Log(Error)")
    plt.title(f"PCA Projection (Explained Var: {pca.explained_variance_ratio_.sum():.2f})")
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    if is_path:
        plt.legend()
    plt.savefig(f"{out_prefix}_pca.png", dpi=150)
    plt.close()
    
    if HAS_UMAP and len(X) >= 5:
        n_neighbors = min(15, len(X) - 1)
        reducer = umap.UMAP(n_components=2, n_neighbors=n_neighbors, min_dist=0.1)
        X_umap = reducer.fit_transform(X)
        
        plt.figure(figsize=(8,6))
        sc = plt.scatter(X_umap[:, 0], X_umap[:, 1], c=E_log, cmap='viridis', s=50)
        if is_path and len(X) >= 2:
            plt.plot(X_umap[:n_init, 0], X_umap[:n_init, 1], 'k--', alpha=0.5, label='Initial Path')
            if len(chain_relaxed) > 0:
                plt.plot(X_umap[n_init:, 0], X_umap[n_init:, 1], 'b-', alpha=0.5, label='Relaxed Path')
        plt.colorbar(sc, label="Log(Error)")
        plt.title(f"UMAP Projection")
        plt.xlabel("UMAP 1")
        plt.ylabel("UMAP 2")
        if is_path:
            plt.legend()
        plt.savefig(f"{out_prefix}_umap.png", dpi=150)
        plt.close()
        
    if HAS_PACMAP and len(X) >= 5:
        n_neighbors = min(10, len(X) - 1)
        reducer = pacmap.PaCMAP(n_components=2, n_neighbors=n_neighbors, MN_ratio=0.5, FP_ratio=2.0)
        X_pacmap = reducer.fit_transform(X, init="pca")
        
        plt.figure(figsize=(8,6))
        sc = plt.scatter(X_pacmap[:, 0], X_pacmap[:, 1], c=E_log, cmap='viridis', s=50)
        if is_path and len(X) >= 2:
            plt.plot(X_pacmap[:n_init, 0], X_pacmap[:n_init, 1], 'k--', alpha=0.5, label='Initial Path')
            if len(chain_relaxed) > 0:
                plt.plot(X_pacmap[n_init:, 0], X_pacmap[n_init:, 1], 'b-', alpha=0.5, label='Relaxed Path')
        plt.colorbar(sc, label="Log(Error)")
        plt.title(f"PacMAP Projection")
        plt.xlabel("PacMAP 1")
        plt.ylabel("PacMAP 2")
        if is_path:
            plt.legend()
        plt.savefig(f"{out_prefix}_pacmap.png", dpi=150)
        plt.close()

    print(f"Saved DR plots to {out_prefix}_*.png")

def run_compare_endpoints(name, dof_file, xyz_path, Gref, seq, distances, angles, kmorse, lepairs, nstep, out_prefix, do_relax=False):
    print(f"\n--- Processing Endpoint {name} (Relax={do_relax}) ---")
    setup_fitreq(dof_file, xyz_path, kMorse=kmorse, Lepairs=lepairs,  softClamp=True, sc_start=4.0, sc_max=6.0, imodel=7, Regularize=1)
    dofs, hyper = parse_dof_file(dof_file)
    if do_relax:
        relax_tethered(dofs, kMorse=kmorse, Lepairs=lepairs, nstep=nstep)
        
    if hasattr(fit, 'compute_model_grid'):
        run_params = dict(nstep=nstep, Fmax=1e-8, dt=0.05, max_step=0.1, damping=0.0)
        Gmodel = fit.compute_model_grid(xyz_path, seq, Gref.shape, do_fit=False, bAddEpairs=True, run_params=run_params, bOutXYZ=False)
        prefix = f"{out_prefix}_{name}_{'relax' if do_relax else 'rigid'}"
        title  = f"Endpoint {name} ({'Relaxed' if do_relax else 'Rigid'})"
        fit.plot_compare(Gref, Gmodel, angles, distances, title, save_prefix=prefix + "_2d", line=True, kcal=True, save_fmt="png")
    else:
        print(f"Warning: compute_model_grid and plot_compare are not available in PN backend. Skipping grid plot for endpoint {name}.")


def run_cluster_ensemble(data_path, out_prefix):
    try:
        import umap
        HAS_UMAP = True
    except ImportError:
        HAS_UMAP = False

    if not os.path.exists(data_path):
        raise FileNotFoundError(f"Data file {data_path} not found.")

    data = np.load(data_path)
    dofs = data['dofs']
    kMorse = data['kMorse']
    Lepairs = data['Lepairs']
    weight_alpha = data['weight_alpha']
    error = data['error']
    dof_names = data['dof_names']

    log_err = np.log10(error)

    if not HAS_UMAP:
        print("UMAP not installed. Falling back to PCA.")
        from sklearn.decomposition import PCA
        reducer = PCA(n_components=2)
        proj = reducer.fit_transform(dofs)
    else:
        reducer = umap.UMAP(n_components=2, n_neighbors=5, min_dist=0.1, random_state=42)
        proj = reducer.fit_transform(dofs)

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()

    props = [
        (log_err, "Log10(Error)", "viridis"),
        (kMorse, "kMorse", "plasma"),
        (Lepairs, "Lepairs", "coolwarm"),
        (weight_alpha, "Weight Alpha", "magma")
    ]

    for ax, (cdata, title, cmap) in zip(axes, props):
        sc = ax.scatter(proj[:, 0], proj[:, 1], c=cdata, cmap=cmap, s=80, alpha=0.8, edgecolor='k')
        ax.set_title(f"Colored by {title}")
        fig.colorbar(sc, ax=ax, label=title)
        
        if title == "Log10(Error)":
            min_idx = np.argmin(cdata)
            max_idx = np.argmax(cdata)
            ax.annotate(f"Min Err ({min_idx})", (proj[min_idx, 0], proj[min_idx, 1]), xytext=(5,5), textcoords='offset points')
            ax.annotate(f"Max Err ({max_idx})", (proj[max_idx, 0], proj[max_idx, 1]), xytext=(5,5), textcoords='offset points')
            
    plt.suptitle("Ensemble Projection (UMAP/PCA) of Optimized Parameters", fontsize=16)
    plt.tight_layout()
    plt.savefig(f"{out_prefix}_all.png", dpi=150)
    print(f"Saved {out_prefix}_all.png")

def getLJ( r, R0=3.0, E0=1.0 ):
    u = R0/r
    u6 = u**6
    return E0*u6*(u6-2.0)

def getMorse( r, R0=3.0, E0=1.0, k=1.6 ):
    e = np.exp( -k*(r-R0) )
    return E0*e*(e-2.0)

def getSR2( r,  E0=1.0, w=1.0 ):
    u = r/w
    return E0*np.exp(-u*u)

def makeCircle( R=3.0 ):
    angs = np.linspace(0.0,2.0*np.pi,100)
    ps = np.zeros( (len(angs),2) )
    ps[:,0] = R*np.cos(angs)
    ps[:,1] = R*np.sin(angs)
    return ps

def _plot_directional( Es, Emax, xs, ys, Rcirc=None, apos=None, ps_circ=None ):
    plt.imshow( Es, origin='lower', vmin=-Emax, vmax=Emax, cmap='bwr', extent=[xs[0],xs[-1],ys[0],ys[-1]] )
    plt.colorbar()
    if apos is not None:
        plt.plot( apos[:,0], apos[:,1], '+k'  )
    if Rcirc is not None:
        ps_circ = makeCircle(Rcirc)    
    if ps_circ is not None:
        plt.plot( ps_circ[:,0], ps_circ[:,1], '--', lw=0.5, c='k' )
    plt.xlim( xs[0], xs[-1] )
    plt.ylim( ys[0], ys[-1] )
    plt.xlabel("x [A]")
    plt.ylabel("y [A]")
    
def run_directional_barrier(L, R, E, k, bMorse, out_file):
    xs = np.linspace( 0.0,8.0,300)
    ys = np.linspace(-3.0,3.0,300)
    Xs,Ys = np.meshgrid(xs,ys)
    Rs = np.sqrt(Xs**2 + Ys**2)

    E0  = E
    R0  = R
    Le  = L

    if bMorse:
        EHb = 1.0
        E_uncorr = getMorse(Rs, R0, E0, k)
    else:
        EHb = 10.0
        E_uncorr = getLJ(Rs, R0, E0)
        
    Rse    = np.sqrt((Xs-Le)**2 + Ys**2)
    E_corr = E_uncorr + getSR2(Rse, -EHb*2.0, w=0.9 )

    ps_circ = makeCircle(R0)
    apos = np.array( [[0.0,0.0],[Le,0.0]] )
    
    plt.figure(figsize=(15,5))
    plt.subplot(1,3,1); _plot_directional( E_uncorr , E0, xs, ys, apos=apos, ps_circ=ps_circ )
    plt.subplot(1,3,2); _plot_directional( E_corr   , E0, xs, ys, apos=apos, ps_circ=ps_circ )
    plt.subplot(1,3,3); _plot_directional( E_corr   , EHb*0.2, xs, ys, apos=apos, ps_circ=ps_circ )
    plt.grid()
    plt.savefig( out_file, bbox_inches='tight' )
    plt.close()

