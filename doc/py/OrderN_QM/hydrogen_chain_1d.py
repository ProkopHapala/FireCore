import argparse
import numpy as np
import matplotlib.pyplot as plt

# unlimited line lenght in numpy when printing
np.set_printoptions(linewidth=np.inf)

from OrderN import (
    find_neighbor_pairs,
    build_matrix_from_pairs,
    plot_matrix,
    canonical_reference,
    plot_spectrum_and_map,
)
import OMM
from FOE import foe_stochastic_density
from GF import (
    greens_function_probing,
    greens_function_random,
    get_contour_poles,
    pole_plot_data,
)


def plot_probe_response(probe_data, probe_idx=0, pole_indices=None, title_prefix="Probe"):
    """Single plot: probe comb, selected pole responses, and integrated rho_probe."""
    if probe_data is None:
        return
    eta = probe_data["eta"]
    rho_probe = probe_data["rho_probe"]
    responses = probe_data["responses"]
    weights = probe_data["weights"]
    if pole_indices is None or len(pole_indices) == 0:
        pole_indices = list(range(min(3, len(responses))))
    pole_indices = [i for i in pole_indices if 0 <= i < len(responses)]
    fig, ax = plt.subplots(figsize=(10, 4))
    x = np.arange(len(eta))
    ax.plot(x, eta.real * 1.0, label=f"{title_prefix} comb idx={probe_idx}", color="k", ls="--", lw=1.5)
    ax.plot(x, rho_probe * -1.0, label="rho_probe (sum all poles)", color="k", lw=3, alpha=0.7)
    npoles = len(pole_indices)
    for j, idx in enumerate(pole_indices):
        resp = responses[idx] * npoles * -2.0
        color = f"C{j%10}"
        ax.plot(x, resp.imag, color=color, lw=1.8, label=f"Im resp pole {idx} w={weights[idx]:.2e}")
    ax.set_xlabel("site")
    ax.set_ylabel("amplitude")
    ax.set_title(f"{title_prefix} responses (imag parts of selected poles)")
    ax.legend(fontsize=8, ncol=2)
    plt.tight_layout()


def plot_probe_set(probes, title="Probes"):
    if not probes:
        return
    fig, ax = plt.subplots(figsize=(10, 3))
    for i, p in enumerate(probes):
        ax.plot(np.arange(len(p)), p, label=f"probe {i}")
    ax.set_title(title)
    ax.set_xlabel("site")
    ax.set_ylabel("amp")
    ax.legend(ncol=3, fontsize=8)
    plt.tight_layout()


def build_dimer_chain(n_cells=10, cell_a=2.5, dimer_distance=0.74):
    """Return positions (N,3) and cell indices for a 1D chain of dimers."""
    atoms = []
    cell_ids = []
    for ic in range(n_cells):
        x0 = ic * cell_a
        atoms.append([x0, 0.0, 0.0])
        atoms.append([x0 + dimer_distance, 0.0, 0.0])
        cell_ids.extend([ic, ic])
    pos = np.array(atoms)
    cell_ids = np.array(cell_ids, dtype=int)
    print(f"#DEBUG build_dimer_chain n_atoms={len(atoms)} cell_a={cell_a} d={dimer_distance}")
    return pos, cell_ids


def build_hamiltonian(pos, cell_ids, onsite=0.0, decay=1.0, t0=-1.0, wrap_ends=False, conf_amp=0.0):
    # include both intra-dimer and inter-cell neighbors
    rmax = 1.05 * max(np.diff(np.sort(np.unique(pos[:,0]))).max(), 1e-6)
    pairs = find_neighbor_pairs(pos, cell_ids=cell_ids, tol=1.05, max_distance=rmax)
    if wrap_ends and len(pos) > 2:
        inter_r = pos[2,0] - pos[1,0] if pos.shape[0] >= 3 else rmax
        pairs.append((len(pos)-1, 0, inter_r))
    H = build_matrix_from_pairs(len(pos), pairs, diag_value=onsite, prefactor=t0, decay=decay)
    if conf_amp != 0.0:
        xs = pos[:,0]
        x0, x1 = xs.min(), xs.max()
        L = max(x1 - x0, 1e-9)
        xc = 0.5*(x0 + x1)
        shape = 0.5 * (1.0 - np.cos(np.pi * (xs - xc) / (0.5*L)))  # 0 at center, 1 at ends
        np.fill_diagonal(H, onsite + conf_amp * shape)
    return H


def build_overlap(pos, cell_ids, diag=1.0, decay=1.0, s0=1.0, wrap_ends=False):
    rmax = 1.05 * max(np.diff(np.sort(np.unique(pos[:,0]))).max(), 1e-6)
    pairs = find_neighbor_pairs(pos, cell_ids=cell_ids, tol=1.05, max_distance=rmax)
    if wrap_ends and len(pos) > 2:
        inter_r = pos[2,0] - pos[1,0] if pos.shape[0] >= 3 else rmax
        pairs.append((len(pos)-1, 0, inter_r))
    S = build_matrix_from_pairs(len(pos), pairs, diag_value=diag, prefactor=s0, decay=decay)
    return S


def plot_hamiltonian(H, title="1D dimer chain Hamiltonian"):
    plot_matrix(H, title=title, symmetric=True)


def plot_overlap(S, title="1D overlap matrix"):
    plot_matrix(S, title=title, symmetric=True)



def plot_geometry(pos, cell_a, title="1D geometry"):
    fig, ax = plt.subplots(figsize=(7, 2.5))
    ax.scatter(pos[:, 0], np.zeros_like(pos[:, 0]), c="k", s=30, label="H")
    # Draw bonds within cell and to next cell (visual only)
    for i in range(0, len(pos), 2):
        ax.plot(pos[i : i + 2, 0], [0, 0], "-", c="tab:blue", lw=2, label="dimer" if i == 0 else None)
    # guide lines for cells
    xmax = pos[:, 0].max()
    for x in np.arange(0, xmax + 1e-6, cell_a):
        ax.axvline(x, color="gray", ls="--", lw=0.5)
    ax.set_xlabel("x [Å]")
    ax.set_yticks([])
    ax.set_title(title)
    ax.legend(loc="upper right")
    plt.tight_layout()

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="1D hydrogen dimer chain tight-binding reference")
    ap.add_argument("--n_cells",  type=int, default=16)
    ap.add_argument("--cell_a",   type=float, default=2.0)
    ap.add_argument("--dimer_distance", type=float, default=0.74)
    # Defaults aligned to Harrison-like scale in eV/Å
    ap.add_argument("--onsite",          type=float, default=-13.6)   # H 1s energy in eV
    ap.add_argument("--decay",           type=float, default=0.60)   # decay length in Å
    ap.add_argument("--t0",              type=float, default=-24.5)  # ~-18 eV at 0.74 Å with decay 0.60
    ap.add_argument("--s0",              type=float, default=1.0)    # overlap scale
    ap.add_argument("--TS_decay_factor", type=float, default=1.2)
    ap.add_argument("--conf_amp",        type=float, default=4.0, help="cosine confinement amplitude added to onsite")
    ap.add_argument("--wrap_ends",       type=int, default=1, help="periodic end-to-end coupling")
    ap.add_argument("--plot_geom",       type=int, default=1, help="show geometry scatter")
    ap.add_argument("--plot_mats",       type=int, default=1, help="show H and S matrices")
    ap.add_argument("--plot_map",        type=int, default=1, help="show MO coefficient map")
    ap.add_argument("--plot_density",    type=int, default=1, help="show Mulliken density")
    ap.add_argument("--no_show",         type=int, default=0, help="skip plt.show (for batch)")
    ap.add_argument("--verbosity",       type=int, default=3, help="printing verbosity; states printed if >2")
    
    ap.add_argument("--orderN_method",   type=str, default="OMM", help="choose one of: FOE | GFcomb | GFrand | OMM")
    ap.add_argument("--foe_solver",      type=str, default="cholesky", help="FOE solver: jacobi|solve|cholesky")
    ap.add_argument("--foe_npoly",       type=int, default=200, help="FOE Chebyshev order")
    ap.add_argument("--foe_beta",        type=float, default=10.0, help="FOE beta (smaller = smoother)")
    ap.add_argument("--foe_no_ortho",    type=int, default=0, help="disable orthogonalization in FOE")
    ap.add_argument("--foe_diag_every",  type=int, default=1, help="if >0 print RMSE/maxdiff vs reference every k Chebyshev step")
    ap.add_argument("--foe_mu",          type=float, default=None, help="override Fermi level for FOE (defaults to HOMO/LUMO mid)")
    ap.add_argument("--foe_nrand",       type=int, default=100, help="number of random probes for FOE trace")
    ap.add_argument("--foe_beta_scaled", type=int, default=0, help="set to 1 if beta already in scaled units (skip automatic scaling)")
    
    ap.add_argument("--gf_n_poles",      type=int, default=64, help="number of contour poles")
    ap.add_argument("--gf_probe_dist",   type=int, default=6, help="probing distance (coloring step)")
    ap.add_argument("--gf_margin_e",     type=float, default=1.0, help="contour bottom offset below min eigen (energy units)")
    ap.add_argument("--gf_probe_idx",    type=int, default=10, help="which comb-probe/color to visualize (0-based)")
    ap.add_argument("--gf_plot_poles",   type=str, default="0,7,15", help="comma-separated pole indices to overlay in probe plot")
    ap.add_argument("--gfr_nrand",       type=int, default=64, help="random probes for GF")
    ap.add_argument("--gfr_probe_idx",   type=int, default=0, help="which probe to visualize (0-based)")
    ap.add_argument("--omm_iter",        type=int, default=8, help="OMM Jacobi iterations")
    ap.add_argument("--omm_damp",        type=float, default=0.2, help="OMM Jacobi damping")
    ap.add_argument("--omm_inertia",     type=float, default=1.0, help="OMM inertial diagonal (stabilizer)")
    ap.add_argument("--omm_max_step",    type=float, default=0.1, help="OMM max coefficient step (abs), set 0 to disable")
    ap.add_argument("--omm_support",     type=int, default=2, help="OMM support in +/- cells")
    ap.add_argument("--omm_plot",        type=int, default=1, help="plot OMM orbitals after orthogonalization")
    ap.add_argument("--omm_init",        type=str, default="const", help="OMM init mode: const | rand")
    ap.add_argument("--omm_rand_scale",  type=float, default=0.1, help="OMM random init scale")
    ap.add_argument("--omm_rand_seed",   type=int, default=None, help="OMM random init seed")
    ap.add_argument("--omm_plot_sel",    type=str, default="0,1", help="comma-separated orbital indices to plot; None=all")
    args = ap.parse_args()
    gf_poles_to_plot = [int(s) for s in args.gf_plot_poles.split(",") if s.strip()!=""]
    method = args.orderN_method.strip().lower()
    if method not in ("foe", "gfcomb", "gfrand", "omm"):
        raise ValueError("orderN_method must be one of: FOE | GFcomb | GFrand | OMM")

    pos, cell_ids = build_dimer_chain(n_cells=args.n_cells, cell_a=args.cell_a, dimer_distance=args.dimer_distance)
    H = build_hamiltonian(pos, cell_ids, onsite=args.onsite, decay=args.decay, t0=args.t0, wrap_ends=bool(args.wrap_ends), conf_amp=args.conf_amp)
    S = build_overlap(pos, cell_ids, diag=1.0, decay=args.decay*args.TS_decay_factor, s0=args.s0, wrap_ends=bool(args.wrap_ends))
    print(f"#DEBUG |H| max={np.abs(H).max():.3f} min={np.abs(H).min():.3f}")
    print(f"#DEBUG |S| max={np.abs(S).max():.3f} min={np.abs(S).min():.3f}")

    if method == "omm":
        C_list, errors, total_abs, max_abs, masks = OMM.run_omm(
            S,
            cell_ids,
            support_cells=args.omm_support,
            wrap_ends=bool(args.wrap_ends),
            n_iter=args.omm_iter,
            damping=args.omm_damp,
            inertia=args.omm_inertia,
            max_step=None if args.omm_max_step == 0 else args.omm_max_step,
            init_mode=args.omm_init,
            rand_scale=args.omm_rand_scale,
            rand_seed=args.omm_rand_seed,
            verbosity=args.verbosity,
        )
        print(f"#INFO OMM total_abs_err={total_abs:.6e} max_abs_err={max_abs:.6e}")
        if args.verbosity >= 3:
            Rmat = OMM._errors_to_matrix(errors, len(masks))
            print("#DEBUG OMM R-matrix\n", Rmat)
        if args.omm_plot:
            C_norm = OMM.normalize_orbitals(C_list, masks, S)
            dense = OMM.coeffs_to_dense(C_norm, masks, n_basis=S.shape[0])
            fig, ax = plt.subplots(figsize=(8, 4))
            xs = np.arange(S.shape[0])
            sel = None
            if args.omm_plot_sel:
                sel = [int(s) for s in args.omm_plot_sel.split(",") if s.strip()!=""]
                sel = [i for i in sel if 0 <= i < dense.shape[0]]
            indices = sel if sel is not None else range(dense.shape[0])
            for i in indices:
                ax.plot(xs, dense[i], label=f"orb {i}")
            ax.set_xlabel("site index")
            ax.set_ylabel("coeff")
            ax.set_title("OMM orbitals after orthogonalization")
            ax.legend(ncol=2, fontsize=8)
            plt.tight_layout()
            if args.no_show:
                plt.close(fig)
            else:
                plt.show()
        exit(0)

    if args.plot_geom:
        plot_geometry(pos, cell_a=args.cell_a, title="1D hydrogen dimer chain geometry")

    if args.plot_mats:
        fig, axes = plt.subplots(1, 2, figsize=(10, 5))
        vmax_h = np.abs(H).max()
        im0 = axes[0].imshow(H, origin="lower", cmap="bwr", vmin=-vmax_h, vmax=vmax_h)
        axes[0].set_title("Hamiltonian")
        axes[0].set_xlabel("j"); axes[0].set_ylabel("i")
        vmax_s = np.abs(S).max()
        im1 = axes[1].imshow(S, origin="lower", cmap="bwr", vmin=-vmax_s, vmax=vmax_s)
        axes[1].set_title("Overlap")
        axes[1].set_xlabel("j"); axes[1].set_ylabel("i")
        fig.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.04)
        fig.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.04)
        plt.tight_layout()

    # canonical O(N^3) reference
    e, C, rho = canonical_reference(H, S, n_occ=len(pos)//2)
    homo = e[len(pos)//2 - 1]
    lumo = e[len(pos)//2]
    gap = lumo - homo
    e_fermi = 0.5 * (homo + lumo)
    valence_width = homo - e[0]
    print(f"#DEBUG eigen min/max {e.min():.3f}/{e.max():.3f} gap={gap:.4f} homo={homo:.4f} lumo={lumo:.4f} valence_width={valence_width:.4f}")
    print(f"#INFO  BandGap={gap:.4f}  HOMO={homo:.4f}  LUMO={lumo:.4f}  ValenceWidth={valence_width:.4f}")
    print(f"#INFO  Fermi(mid HOMO/LUMO) = {e_fermi:.4f}")
    # report representative hoppings/overlaps
    if len(pos) >= 4:
        t_intra = H[0, 1]
        s_intra = S[0, 1]
        t_inter = H[1, 2]
        s_inter = S[1, 2]
        print(f"#INFO  T on/off ( {t_intra:.4f} , {t_inter:.4f} )")
        print(f"#INFO  S on/off ( {s_intra:.4f} , {s_inter:.4f} )")
    # dump eigenstates (controlled by verbosity)
    if args.verbosity > 2:
        for idx, ei in enumerate(e):
            coeffs = " ".join(f"{c:.3f}" for c in C[:, idx])
            print(f"#STATE {idx:02d} E={ei:.6f}  coeffs: {coeffs}")

    if args.plot_map:
        info = f"gap={gap:.3f}\nHOMO={homo:.3f}\nLUMO={lumo:.3f}\nvalence width={valence_width:.3f}"
        plot_spectrum_and_map(e, C, n_occ=len(pos)//2, positions=pos[:,0], plot_spectrum=False, info_text=info)

    if args.plot_density:
        fig, ax = plt.subplots(figsize=(7, 2.5))
        xs = pos[:,0]
        order = np.argsort(xs)
        ax.plot(xs[order], rho[order], ":k", lw=2, label="density (Mulliken)")
        ax.set_xlabel("x [Å]")
        ax.set_ylabel("ρ")
        ax.set_title("Total density along chain")
        ax.legend()
        plt.tight_layout()

    # FOE (order-N) stochastic density estimate for comparison
    if method == "foe":
        rho_foe, foe_info = foe_stochastic_density( H, S, n_poly=args.foe_npoly, n_random=args.foe_nrand, beta=args.foe_beta,
            mu=args.foe_mu if args.foe_mu is not None else e_fermi,
            jacobi_steps=12,
            solver=args.foe_solver,
            orthogonalize=not bool(args.foe_no_ortho),
            ref_rho=rho,
            diag_every=args.foe_diag_every,
            diag_prefix="#INFO FOE_iter",
            beta_is_scaled=bool(args.foe_beta_scaled),
        )
        print(f"#DEBUG FOE rho min/max {rho_foe.min():.4f}/{rho_foe.max():.4f} info={foe_info}")
        # per-site report
        xs = pos[:,0]
        order = np.argsort(xs)
        print("#INFO  Site  x  rho_ref  rho_FOE  diff")
        for idx in order:
            diff = rho_foe[idx] - rho[idx]
            print(f"#INFO  {idx:3d}  {xs[idx]:7.3f}  {rho[idx]:9.5f}  {rho_foe[idx]:9.5f}  {diff:9.5f}")
        diff_abs = np.abs(rho_foe - rho)
        print(f"#INFO  Summary rho_ref min/max=({rho.min():.5f},{rho.max():.5f}) "
              f"rho_FOE min/max=({rho_foe.min():.5f},{rho_foe.max():.5f}) "
              f"max|diff|={diff_abs.max():.5f}")
        if args.plot_density:
            ax.plot(xs[order], rho_foe[order], "-r", lw=0.5, label="density (FOE)")
            ax.legend()

    # Deterministic Green-function probing (no stochastic)
    if method == "gfcomb":
        z_nodes, w_nodes = get_contour_poles(emin=e.min(), e_fermi=e_fermi, n_poles=args.gf_n_poles, margin_e=args.gf_margin_e)
        rho_gf, probe_data = greens_function_probing( H, S, mu=e_fermi, emin=e.min(), n_poles=args.gf_n_poles, probing_distance=args.gf_probe_dist, solver="dense", return_probe_idx=args.gf_probe_idx, store_probes=True )
        diff_gf = rho_gf - rho
        print("#INFO  Method    Site   Charge    Diff (GF probing)")
        for i in order:
            print(f"#INFO  GF-Probe {i:02d}  {rho_gf[i]:9.5f}  {diff_gf[i]:+9.5f}")
        print(f"#INFO  GF Summary rho_ref min/max=({rho.min():.5f},{rho.max():.5f}) rho_GF min/max=({rho_gf.min():.5f},{rho_gf.max():.5f}) max|diff|={np.abs(diff_gf).max():.5f}")
        if args.plot_density:
            ax.plot(xs[order], rho_gf[order], "-g", lw=0.5, label="density (GF comb)")
            ax.legend()
        # Pole plot vs spectrum
        pole_data = pole_plot_data(e, z_nodes)
        fig_p, ax_p = plt.subplots(figsize=(6,4))
        ax_p.scatter(pole_data["eigs"], np.zeros_like(pole_data["eigs"]), c="k", s=20, label="eigenvalues")
        ax_p.scatter(pole_data["pole_re"], pole_data["pole_im"], c="tab:orange", s=30, label="poles (upper arc)")
        ax_p.axvline(e_fermi, color="tab:green", ls="--", label="Fermi")
        ax_p.set_xlabel("Energy")
        ax_p.set_ylabel("Im(z)")
        ax_p.set_title("Contour poles vs spectrum")
        ax_p.legend()
        if not args.no_show:
            plot_probe_response(probe_data, probe_idx=args.gf_probe_idx, pole_indices=gf_poles_to_plot, title_prefix="GF-probe (comb)")
            plot_probe_set(probe_data.get("probes"), title="GF comb probes")

    # Random-probe Green-function with per-probe visualization
    if method == "gfrand":
        rho_gfr, probe_data = greens_function_random( H, S, mu=e_fermi, emin=e.min(), n_poles=args.gf_n_poles, n_random=args.gfr_nrand, solver="dense", return_probe_idx=args.gfr_probe_idx, store_probes=True )
        diff_gfr = rho_gfr - rho
        print("#INFO  Method    Site   Charge    Diff (GF random)")
        for i in order:
            print(f"#INFO  GF-Rand {i:02d}  {rho_gfr[i]:9.5f}  {diff_gfr[i]:+9.5f}")
        print(f"#INFO  GF-Rand Summary rho_ref min/max=({rho.min():.5f},{rho.max():.5f}) rho_GF min/max=({rho_gfr.min():.5f},{rho_gfr.max():.5f}) max|diff|={np.abs(diff_gfr).max():.5f}")
        if args.plot_density:
            ax.plot(xs[order], rho_gfr[order], "-m", lw=0.5, label="density (GF random)")
            ax.legend()

        if probe_data and not args.no_show:
            plot_probe_response(probe_data, probe_idx=args.gfr_probe_idx, title_prefix="GF-probe (random)")
            plot_probe_set(probe_data.get("probes"), title="GF random probes")

    if not args.no_show:
        plt.show()
