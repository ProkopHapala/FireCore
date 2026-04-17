"""
STM analysis CLI: Compare transport methods and generate STM images.

PURPOSE:
  Main user-facing script for STM transport analysis using pre-computed DOS files.
  Supports multiple transport methods (PDOS approximation, NEGF Caroli, NEGF response),
  scanning protocols (height scan, lateral scan, 2D xy scan), and MO overlap diagnostics.

PHYSICAL BACKGROUND:

1) Transport Methods:

   a) PDOS approximation (Phase 1):
      Uses Tersoff-Hamann-like approximation with Wolfsberg-Helmholtz hopping:
        V_μν ≈ κ * (H_μμ + H_νν)/2 * S_μν
        I ≈ Σ_μν PDOS_μ^tip(E) |V_μν|² PDOS_ν^sample(E)
      This is a fast perturbative approach valid for weak coupling.

   b) NEGF Caroli (combined-system):
        G^r(E) = [(E+iη)S - H - Σ_L - Σ_R]^{-1}
        T(E) = Tr[Γ_L G^r Γ_R G^a]
      Full quantum transport with wide-band limit self-energies on lead orbitals.

   c) NEGF response metric (deterministic probe):
        A = (E+iη)S - H - Σ_L - Σ_R
        x = solve(A, b)  # b = source vector on lead L
        resp = x† Γ_measure x
      Single-solve probe of coupling strength (not a transmission coefficient).

2) Inter-System Coupling:
   Vacuum exponential radial decay + Slater-Koster angular dependence:
        H_μν(r) = A_μν * exp(-β(r - r₀)) * SK_angular(l,m,n)
   Parameters:
     - β ~1.0 Å⁻¹ (vacuum decay constant)
     - r₀ ~3.0 Å (reference distance)
     - A_μν: orbital-pair strengths (A_ss, A_sp, A_ppsig, A_pppi)

3) MO Overlap Diagnostic:
   For "featureless metallic tip" approximation:
        v_smp(s) = Σ_{i∈tip} H_{TS}(i,:)
        t_n(s) = v_smp(s) · C_{:,n}
        I(s) = |t_n(s)|²
   Computes signed amplitude and intensity for a selected MO (e.g., HOMO).
   The signed amplitude reveals nodal planes (sign flips), while intensity matches STM measurement.

4) Scanning Protocols:
   - Height scan: I(z) at fixed lateral position (tests tunneling decay)
   - Lateral scan: I(s) at constant z (diagonal, x, or y direction)
   - 2D xy scan: I(x,y) at constant z (STM image)

USAGE EXAMPLES:

  Basic NEGF analysis with response metric:
    cd tests/pyFireball
    python stm_analysis.py --method negf --negf_mode response --E -0.8

  1D lateral scan at constant z:
    python stm_analysis.py --method negf --scan_xy --zxy 5.0 --xy_dir diag

  MO overlap diagnostic (requires DOS with eigenvectors):
    python stm_analysis.py --method negf --scan_xy --scan_mo --mo 42

  Full 2D STM image:
    python stm_analysis.py --method negf --z2d 5.0 --xy_range 7.0

PREREQUISITES:
  Pre-computed DOS files (run once per molecule):
    python stm_compute_dos.py --xyz H2.xyz --out export/stm_phase1/dos_H.npz
    python stm_compute_dos.py --xyz PTCDA.xyz --out export/stm_phase1/dos_PTCDA.npz

  Note: For MO overlap diagnostic, DOS files must contain 'eps' and 'C' arrays
  (regenerate with updated stm_compute_dos.py after recent changes).

KEY CLI OPTIONS:
  --method {pdos,negf}       Transport method (pdos=approx, negf=full)
  --solver {direct,iter}     NEGF solver (direct=exact, iter=GMRES)
  --negf_mode {caroli,response}  Transmission estimator
  --E VALUE                 Single energy (eV); default uses Fermi level
  --source {L,R}            Response excitation lead
  --measure {L,R}           Response measurement lead
  --beta, --r0              Vacuum decay parameters
  --A_ss, --A_sp, --A_ppsig, --A_pppi  Coupling strengths
  --scan_xy                 Enable 1D lateral scan
  --scan_mo                 Enable MO overlap diagnostic
  --z2d, --zxy             Heights for 2D and lateral scans
  --xy_range, --xy_step     Lateral scan parameters

OUTPUTS:
  All outputs saved to export/stm_phase1/ (or custom --export_dir):
    - dos_plot.png: PDOS curves for tip and sample
    - scan1d.png: Height scan I(z)
    - scan1d_xy.png: Lateral scan I(s)
    - scan1d_xy_mo.png: MO overlap signed amplitude and intensity
    - stm2d.png: 2D xy scan I(x,y)
  Corresponding .npz files with numerical data.
"""
import argparse, os, sys, time
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

_THIS_DIR  = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.normpath(os.path.join(_THIS_DIR, "..", ".."))
if _REPO_ROOT not in sys.path: sys.path.insert(0, _REPO_ROOT)

from pyBall.FireballOCL.STM import (
    compute_hopping, stm_current,
    project_pdos_to_grid,
    build_inter_system_blocks_fdata, build_inter_system_blocks_exp_sk, assemble_combined_HS,
    select_atoms_by_xy_radius, select_lead_orbitals_by_atoms,
    negf_current, negf_current_iterative, negf_response_at_energy,
    coupling_vec_tip_to_sample, mo_overlap_amplitude
)
from pyBall.FireballOCL.STM_utils import (
    set_export_dir, save_plot, plot_dos, plot_atoms,
    DEFAULT_EXPORT_DIR
)

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--tip",       default="H2",  choices=["H2","CO"])
    p.add_argument("--method",    default="pdos", choices=["pdos","negf"])  # pdos: Phase1 approx; negf: combined-system Caroli
    p.add_argument("--bias",      type=float, default=0.1)
    p.add_argument("--z2d",       type=float, default=5.0)
    p.add_argument("--xy_range",  type=float, default=7.0)
    p.add_argument("--xy_step",   type=float, default=0.5)
    p.add_argument("--scan_xy",   action="store_true")      # extra 1D lateral scan at constant z
    p.add_argument("--zxy",       type=float, default=5.0)   # z for lateral scan
    p.add_argument("--xy_dir",    default="diag", choices=["diag","x","y"])  # scan direction
    p.add_argument("--scan_mo",   action="store_true")      # MO overlap diagnostic along lateral scan
    p.add_argument("--mo",        type=int, default=None)    # explicit MO index (0-based in eps/C)
    p.add_argument("--mo_win",    type=float, default=0.05)  # energy window (eV) around E to select closest MO
    p.add_argument("--reduce_tip", default="sum", choices=["sum","norm"])  # reduce H_TS over tip orbitals
    p.add_argument("--no2d",      action="store_true")
    p.add_argument("--noshow",    action="store_true")
    p.add_argument("--export_dir", default=None)
    p.add_argument("--gammaL",    type=float, default=2.0)   # WBL coupling strength (eV) for tip lead
    p.add_argument("--gammaR",    type=float, default=2.0)   # WBL coupling strength (eV) for sample lead
    p.add_argument("--eta",       type=float, default=1e-4)  # small imaginary for numerical stability
    p.add_argument("--solver",    default="direct", choices=["direct","iter"])  # NEGF solver: direct inverse vs iterative GMRES
    p.add_argument("--negf_mode", default="response", choices=["caroli","response"])  # transmission estimator
    p.add_argument("--E",         type=float, default=None)  # single energy (eV); default uses E_F
    p.add_argument("--source",    default="L", choices=["L","R"])                # response excitation lead
    p.add_argument("--measure",   default="R", choices=["L","R","mol"])         # response measurement projector
    p.add_argument("--nrand",     type=int, default=8)       # Hutchinson samples for iterative trace
    p.add_argument("--tol",       type=float, default=1e-6)  # GMRES tolerance
    p.add_argument("--maxiter",   type=int, default=400)     # GMRES maxiter
    p.add_argument("--rcut",      type=float, default=6.0)   # inter-system coupling cutoff (Å)
    p.add_argument("--beta",      type=float, default=1.0)   # exponential decay for tunneling coupling (Å^-1)
    p.add_argument("--r0",        type=float, default=3.0)   # reference distance for exponential (Å)
    p.add_argument("--A_ss",      type=float, default=-1.0)  # SK channel amplitudes (arb/eV)
    p.add_argument("--A_sp",      type=float, default=-1.0)
    p.add_argument("--A_ppsig",   type=float, default=-1.0)
    p.add_argument("--A_pppi",    type=float, default=+1.0)
    p.add_argument("--overlap_scale", type=float, default=0.0)  # S_TS = overlap_scale * H_TS ; 0 -> S_TS=0
    p.add_argument("--zmax",      type=float, default=10.0)  # max z for 1D scan
    p.add_argument("--E0",        type=float, default=None)  # integration window start (eV)
    p.add_argument("--E1",        type=float, default=None)  # integration window end (eV)
    args = p.parse_args()

    # Set output directory
    export_dir = args.export_dir if args.export_dir else DEFAULT_EXPORT_DIR
    set_export_dir(export_dir)
    os.makedirs(export_dir, exist_ok=True)
    results = {}

    # ═══ LOAD DOS ════════════════════════════════════════════════════════════
    tip_file = os.path.join(export_dir, f"dos_{args.tip}.npz")
    smp_file = os.path.join(export_dir, "dos_PTCDA.npz")
    for f in [tip_file, smp_file]:
        if not os.path.exists(f):
            raise FileNotFoundError(f"Missing DOS file: {f}\nRun stm_compute_dos.py first.")

    tip = np.load(tip_file, allow_pickle=True)
    smp = np.load(smp_file, allow_pickle=True)
    print(f"Loaded tip DOS: {tip_file}  shape A={tip['A'].shape}")
    print(f"Loaded smp DOS: {smp_file}  shape A={smp['A'].shape}")

    E_grid = smp['E_grid']; E_F = float(smp['E_fermi']); dE = float(E_grid[1]-E_grid[0])
    A_smp = smp['A']
    A_tip_raw = tip['A']; E_tip = tip['E_grid']; E_F_tip = float(tip['E_fermi'])
    A_tip = np.array([np.interp(E_grid, E_tip-E_F_tip+E_F, A_tip_raw[:,m], left=0, right=0)
                      for m in range(A_tip_raw.shape[1])]).T

    print(f"\nE_F = {E_F:.4f} eV  (sample HOMO)")
    print(f"E_grid: [{E_grid[0]:.3f}, {E_grid[-1]:.3f}]  nE={len(E_grid)}")
    E_single = float(args.E) if args.E is not None else float(E_F)

    # ═══ T3: Sum-rule + broadening ═══════════════════════════════════════════
    print("\n" + "="*65)
    print("T3: Sum-rule and broadening check")
    print("="*65)
    sr_smp = (A_smp*dE).sum(0); sr_tip = (A_tip*dE).sum(0)
    print(f"  Sample sum-rule: mean={sr_smp.mean():.3f} min={sr_smp.min():.3f} max={sr_smp.max():.3f}")
    print(f"  Tip    sum-rule: mean={sr_tip.mean():.3f} min={sr_tip.min():.3f} max={sr_tip.max():.3f}")
    results['T3_sumrule_smp'] = 'PASSED' if sr_smp.mean()>0.05 else 'FAILED'
    results['T3_sumrule_tip'] = 'PASSED' if sr_tip.mean()>0.01 else 'FAILED'

    contact_orbs = tip['contact_orbs']
    non_contact  = np.setdiff1d(np.arange(A_tip.shape[1]), contact_orbs)
    if len(contact_orbs)>0 and len(non_contact)>0:
        w_c  = (A_tip[:,contact_orbs]*dE).sum()  / len(contact_orbs)
        w_nc = (A_tip[:,non_contact]*dE).sum()   / len(non_contact)
        print(f"  Tip broadening: contact={w_c:.4f}  non-contact={w_nc:.4f}")
        results['T3_broadening'] = 'PASSED' if w_c >= w_nc*0.5 else 'FAILED'
    else: results['T3_broadening'] = 'N/A'

    # DOS plots
    fig, axes = plt.subplots(1,2,figsize=(14,5))
    fig.suptitle(f"DOS — {args.tip} tip | PTCDA sample (Γ_tip={float(tip['gamma']):.2f}eV, Γ_smp={float(smp['gamma']):.2f}eV)")
    plot_dos(axes[0], E_grid, A_tip, tip['eigen'], f"Tip ({args.tip}) PDOS", E_F)
    plot_dos(axes[1], E_grid, A_smp, smp['eigen'], "Sample (PTCDA) PDOS",   E_F)
    plt.tight_layout(); save_plot(fig, "dos_plot", export_dir, dpi=130)

    # ═══ Transport method selection ══════════════════════════════════════════
    use_negf = (args.method == 'negf')

    # ═══ T4: Hopping / Coupling at reference position ═══════════════════════
    print("\n" + "="*65)
    print("T4: Coupling at z=5Å above PTCDA centre")
    print("="*65)
    tip_types  = tip['atomTypes']; tip_pos  = tip['atomPos'];   H_tip_diag = np.diag(tip['H'])
    smp_types  = smp['atomTypes']; smp_pos  = smp['atomPos'];   H_smp_diag = np.diag(smp['H'])
    o2a_tip    = tip['orb2atom'].astype(int); o2a_smp = smp['orb2atom'].astype(int)

    if not use_negf:
        Vsq_ref, R_ref = compute_hopping(H_tip_diag, o2a_tip, tip_types, tip_pos,
                                         H_smp_diag, o2a_smp, smp_types, smp_pos,
                                         np.array([0,0,5.0]))
        print(f"  |V|² shape={Vsq_ref.shape}  max={Vsq_ref.max():.3e}  mean={Vsq_ref[Vsq_ref>0].mean():.3e}")
        results['T4_hopping'] = 'PASSED' if Vsq_ref.max()>0 else 'FAILED'

        fig2,ax2 = plt.subplots(figsize=(8,4))
        r_f = R_ref.ravel(); v_f = Vsq_ref.ravel()
        ax2.semilogy(r_f[np.argsort(r_f)], v_f[np.argsort(r_f)]+1e-30, '.', ms=1.5, alpha=0.4)
        ax2.set_xlabel("r (Å)"); ax2.set_ylabel("|V_μν|² (eV²)"); ax2.set_title("Hopping vs distance (z=5Å)")
        plt.tight_layout(); save_plot(fig2, "hopping_vs_dist", export_dir, dpi=120)
    else:
        # Build inter-system blocks from Fdata SK tables (sp-only 4x4, sliced by norb_per)
        fdata_dir = os.path.join(_THIS_DIR, 'Fdata')
        tip_pos_ref = tip_pos + np.array([0.0, 0.0, 5.0])
        pairs, Hb, Sb = build_inter_system_blocks_fdata(tip_types, tip_pos_ref, tip['norb_per'], smp_types, smp_pos, smp['norb_per'],
                                                        fdata_dir=fdata_dir, rcut=args.rcut, root_H='kinetic', root_S='overlap',
                                                        exp_beta=args.beta, exp_r0=3.0)
        print(f"  Inter-system pairs={len(pairs)}  (rcut={args.rcut}Å)")
        # Coupling norm diagnostic
        if len(Hb) > 0:
            mx = max(float(np.max(np.abs(b))) for b in Hb)
            results['T4_hopping'] = 'PASSED' if mx>0 else 'FAILED'
            print(f"  max|H_TS block| = {mx:.3e}")
        else:
            results['T4_hopping'] = 'FAILED'
            print("  ERROR: no coupling pairs found (rcut too small or missing Fdata tables?)")

    # ═══ T5: 1D height scan ══════════════════════════════════════════════════
    print("\n" + "="*65)
    print("T5: 1D height scan I(z)")
    print("="*65)
    z_arr = np.linspace(3.0, float(args.zmax), 36); I_arr = np.zeros(len(z_arr)); T_arr = np.zeros((len(z_arr), len(E_grid)))
    t0=time.time()
    if not use_negf:
        for iz,z in enumerate(z_arr):
            Vsq,_ = compute_hopping(H_tip_diag,o2a_tip,tip_types,tip_pos,
                                    H_smp_diag,o2a_smp,smp_types,smp_pos, np.array([0,0,z]))
            I,T   = stm_current(A_tip, A_smp, Vsq, E_grid, E_F, args.bias)
            I_arr[iz]=I; T_arr[iz]=T
            print(f"  z={z:.2f}Å  I={I:.3e}  T_max={T.max():.3e}")
    else:
        # Lead orbital selection: tip = all orbitals; sample = outer ring by xy radius
        leadL = np.arange(tip['H'].shape[0], dtype=int)
        smp_atoms_outer = select_atoms_by_xy_radius(smp_pos, rmin=0.8*np.max(np.linalg.norm(smp_pos[:,:2]-smp_pos[:,:2].mean(0),axis=1)))
        leadR_s = select_lead_orbitals_by_atoms(smp['orb2atom'], smp_atoms_outer)
        leadR = tip['H'].shape[0] + leadR_s
        E0 = args.E0 if args.E0 is not None else E_F
        E1 = args.E1 if args.E1 is not None else (E_F + args.bias)
        for iz,z in enumerate(z_arr):
            tip_pos_z = tip_pos + np.array([0.0, 0.0, z])
            pairs, Hb, Sb = build_inter_system_blocks_exp_sk(tip_types, tip_pos_z, tip['norb_per'], smp_types, smp_pos, smp['norb_per'],
                                                             rcut=args.rcut, beta=args.beta, r0=args.r0,
                                                             A_ss=args.A_ss, A_sp=args.A_sp, A_pp_sig=args.A_ppsig, A_pp_pi=args.A_pppi,
                                                             overlap_scale=args.overlap_scale)
            Hc, Sc = assemble_combined_HS(tip['H'], tip['S'], tip['orb2atom'], tip['norb_per'],
                                          smp['H'], smp['S'], smp['orb2atom'], smp['norb_per'],
                                          tip_pos_z, smp_pos, pairs, Hb, Sb, tip_offset=np.array([0,0,z]))
            if args.negf_mode == 'response':
                resp = negf_response_at_energy(E_single, Hc, Sc, leadL, leadR, gammaL=args.gammaL, gammaR=args.gammaR, eta=args.eta, source=args.source, measure=args.measure)
                I_arr[iz] = resp
                T_arr[iz, :] = 0.0
                print(f"  z={z:.2f}Å  resp(E={E_single:.3f})={resp:.3e}")
            else:
                # Caroli: single-energy only unless you explicitly want full E_grid (not supported here)
                if args.solver == 'iter':
                    resp = negf_current_iterative(np.array([E_single]), Hc, Sc, leadL, leadR, gammaL=args.gammaL, gammaR=args.gammaR, eta=args.eta,
                                                  E0=E_single, E1=E_single, nrand=args.nrand, tol=args.tol, maxiter=args.maxiter, seed=0, reuse=True)[1][0]
                else:
                    resp = negf_current(np.array([E_single]), Hc, Sc, leadL, leadR, gammaL=args.gammaL, gammaR=args.gammaR, eta=args.eta,
                                        E0=E_single, E1=E_single)[1][0]
                I_arr[iz] = resp
                T_arr[iz, :] = 0.0
                print(f"  z={z:.2f}Å  T(E={E_single:.3f})={resp:.3e}")
    print(f"  1D scan done in {time.time()-t0:.1f}s")

    logI = np.log(I_arr+1e-40); coef = np.polyfit(z_arr, logI, 1); kappa = -coef[0]
    print(f"  Decay constant κ = {kappa:.3f} Å⁻¹  (typical STM: ~1–2 Å⁻¹)")
    results['T5_decay']  = 'PASSED' if 0.1<kappa<10 else 'FAILED'
    results['T5_1dscan'] = 'PASSED'

    np.savez(os.path.join(export_dir,"scan1d.npz"), z=z_arr, I=I_arr, T=T_arr, E_grid=E_grid, kappa=kappa, method=args.method, solver=args.solver, negf_mode=args.negf_mode, E=E_single)
    fig3,axes3 = plt.subplots(1,2,figsize=(13,4))
    axes3[0].semilogy(z_arr, I_arr+1e-40, 'o-', lw=2)
    axes3[0].semilogy(z_arr, np.exp(coef[1]+coef[0]*z_arr), 'r--', label=f'κ={kappa:.2f}Å⁻¹')
    axes3[0].legend(); axes3[0].set_xlabel("z (Å)"); axes3[0].set_ylabel("I (arb)"); axes3[0].set_title(f"I(z) [{args.tip} tip]")
    im3 = axes3[1].imshow(T_arr, origin='lower', aspect='auto',
                          extent=[E_grid[0],E_grid[-1],z_arr[0],z_arr[-1]], cmap='hot')
    axes3[1].axvline(E_F, color='w', lw=0.8, ls='--'); axes3[1].set_xlabel("E (eV)"); axes3[1].set_ylabel("z (Å)")
    axes3[1].set_title("T(E,z) map"); plt.colorbar(im3,ax=axes3[1])
    plt.tight_layout(); save_plot(fig3, "scan1d", export_dir, dpi=130)

    # ═══ Optional: 1D lateral scan at constant z ═════════════════════════════
    if args.scan_xy:
        print("\n" + "="*65)
        print(f"T5b: 1D lateral scan at z={args.zxy:.2f}Å  dir={args.xy_dir}")
        print("="*65)
        t0=time.time()
        s_arr = np.arange(-args.xy_range, args.xy_range+args.xy_step*0.5, args.xy_step)
        Ixy = np.zeros(len(s_arr), dtype=float)
        t_mo = np.zeros(len(s_arr), dtype=float)
        t2_mo = np.zeros(len(s_arr), dtype=float)

        # Pick MO index on sample (PTCDA)
        mo_idx = None
        if args.scan_mo:
            if ('C' not in smp) or ('eps' not in smp):
                raise KeyError("Sample dos_PTCDA.npz missing 'C'/'eps'. Re-run stm_compute_dos.py to regenerate with eigenvectors.")
            eps = smp['eps']
            if args.mo is not None:
                mo_idx = int(args.mo)
            else:
                # choose closest MO to E_single
                mo_idx = int(np.argmin(np.abs(eps - float(E_single))))
                if abs(float(eps[mo_idx]) - float(E_single)) > float(args.mo_win):
                    print(f"[WARN] closest MO eps[{mo_idx}]={eps[mo_idx]:.4f} far from E={E_single:.4f} by {abs(eps[mo_idx]-E_single):.3f}eV")
            print(f"  MO selected: index={mo_idx}  eps={float(eps[mo_idx]):.4f} eV")
        for i,s in enumerate(s_arr):
            if args.xy_dir == 'diag':
                x = float(s); y = float(s)
            elif args.xy_dir == 'x':
                x = float(s); y = 0.0
            else:
                x = 0.0; y = float(s)
            if not use_negf:
                Vsq,_ = compute_hopping(H_tip_diag,o2a_tip,tip_types,tip_pos,
                                        H_smp_diag,o2a_smp,smp_types,smp_pos, np.array([x,y,float(args.zxy)]))
                Ixy[i],_ = stm_current(A_tip, A_smp, Vsq, E_grid, E_F, args.bias)
            else:
                tip_pos_xy = tip_pos + np.array([x, y, float(args.zxy)])
                pairs, Hb, Sb = build_inter_system_blocks_exp_sk(tip_types, tip_pos_xy, tip['norb_per'], smp_types, smp_pos, smp['norb_per'],
                                                                 rcut=args.rcut, beta=args.beta, r0=args.r0,
                                                                 A_ss=args.A_ss, A_sp=args.A_sp, A_pp_sig=args.A_ppsig, A_pp_pi=args.A_pppi,
                                                                 overlap_scale=args.overlap_scale)

                if args.scan_mo:
                    v_s = coupling_vec_tip_to_sample(tip['H'].shape[0], smp['H'].shape[0], tip['norb_per'], smp['norb_per'], pairs, Hb, reduce_tip=args.reduce_tip)
                    tt, tt2 = mo_overlap_amplitude(v_s, smp['C'], mo_idx)
                    t_mo[i] = tt
                    t2_mo[i] = tt2

                Hc, Sc = assemble_combined_HS(tip['H'], tip['S'], tip['orb2atom'], tip['norb_per'],
                                              smp['H'], smp['S'], smp['orb2atom'], smp['norb_per'],
                                              tip_pos_xy, smp_pos, pairs, Hb, Sb, tip_offset=np.array([x,y,float(args.zxy)]))
                if args.negf_mode == 'response':
                    Ixy[i] = negf_response_at_energy(E_single, Hc, Sc, leadL, leadR, gammaL=args.gammaL, gammaR=args.gammaR, eta=args.eta, source=args.source, measure=args.measure)
                else:
                    Ixy[i] = negf_current(np.array([E_single]), Hc, Sc, leadL, leadR, gammaL=args.gammaL, gammaR=args.gammaR, eta=args.eta,
                                          E0=E_single, E1=E_single)[1][0]
            if (i % max(1, len(s_arr)//10)) == 0:
                print(f"  s={s:7.3f}Å  I={Ixy[i]:.3e}")
        print(f"  lateral scan done in {time.time()-t0:.2f}s")
        np.savez(os.path.join(export_dir,"scan1d_xy.npz"), s=s_arr, I=Ixy, z=float(args.zxy), dir=args.xy_dir, method=args.method, solver=args.solver, negf_mode=args.negf_mode, E=E_single)
        fig5,ax5 = plt.subplots(1,1,figsize=(7,4))
        ax5.plot(s_arr, Ixy, 'o-', lw=1.5)
        ax5.set_xlabel("s (Å)")
        ax5.set_ylabel("I (arb)")
        ax5.set_title(f"1D lateral scan z={args.zxy:.2f}Å dir={args.xy_dir} E={E_single:.3f}eV")
        plt.tight_layout(); save_plot(fig5, "scan1d_xy", export_dir, dpi=130)
        print(f"  → {os.path.join(export_dir, 'scan1d_xy.png')}")

        if args.scan_mo:
            np.savez(os.path.join(export_dir,"scan1d_xy_mo.npz"), s=s_arr, t=t_mo, t2=t2_mo, I=Ixy, mo=mo_idx, eps=float(smp['eps'][mo_idx]), z=float(args.zxy), dir=args.xy_dir, reduce_tip=args.reduce_tip)
            fig6,axs6 = plt.subplots(2,1,figsize=(7,6), sharex=True)
            axs6[0].plot(s_arr, t_mo, 'o-', lw=1.2)
            axs6[0].axhline(0.0, color='k', lw=0.8)
            axs6[0].set_ylabel("t(s) signed")
            axs6[0].set_title(f"MO overlap amplitude (signed)  mo={mo_idx}  eps={float(smp['eps'][mo_idx]):.3f}eV  reduce_tip={args.reduce_tip}")
            axs6[1].plot(s_arr, t2_mo, 'o-', lw=1.2)
            axs6[1].set_ylabel("|t|^2")
            axs6[1].set_xlabel("s (Å)")
            plt.tight_layout(); save_plot(fig6, "scan1d_xy_mo", export_dir, dpi=130)
            print(f"  → {os.path.join(export_dir, 'scan1d_xy_mo.png')}")

    # ═══ T6: 2D xy scan ══════════════════════════════════════════════════════
    if not args.no2d:
        print("\n" + "="*65)
        print(f"T6: 2D xy scan I(x,y) at z={args.z2d}Å")
        print("="*65)
        x_arr = np.arange(-args.xy_range, args.xy_range+args.xy_step*0.5, args.xy_step)
        y_arr = np.arange(-args.xy_range, args.xy_range+args.xy_step*0.5, args.xy_step)
        I2d = np.zeros((len(x_arr), len(y_arr))); t0=time.time()
        if not use_negf:
            for ix,x in enumerate(x_arr):
                for iy,y in enumerate(y_arr):
                    Vsq,_ = compute_hopping(H_tip_diag,o2a_tip,tip_types,tip_pos,
                                            H_smp_diag,o2a_smp,smp_types,smp_pos, np.array([x,y,args.z2d]))
                    I2d[ix,iy],_ = stm_current(A_tip, A_smp, Vsq, E_grid, E_F, args.bias)
        else:
            leadL = np.arange(tip['H'].shape[0], dtype=int)
            smp_atoms_outer = select_atoms_by_xy_radius(smp_pos, rmin=0.8*np.max(np.linalg.norm(smp_pos[:,:2]-smp_pos[:,:2].mean(0),axis=1)))
            leadR_s = select_lead_orbitals_by_atoms(smp['orb2atom'], smp_atoms_outer)
            leadR = tip['H'].shape[0] + leadR_s
            for ix,x in enumerate(x_arr):
                for iy,y in enumerate(y_arr):
                    tip_pos_xy = tip_pos + np.array([float(x), float(y), float(args.z2d)])
                    pairs, Hb, Sb = build_inter_system_blocks_exp_sk(tip_types, tip_pos_xy, tip['norb_per'], smp_types, smp_pos, smp['norb_per'],
                                                                     rcut=args.rcut, beta=args.beta, r0=args.r0,
                                                                     A_ss=args.A_ss, A_sp=args.A_sp, A_pp_sig=args.A_ppsig, A_pp_pi=args.A_pppi,
                                                                     overlap_scale=args.overlap_scale)
                    Hc, Sc = assemble_combined_HS(tip['H'], tip['S'], tip['orb2atom'], tip['norb_per'],
                                                  smp['H'], smp['S'], smp['orb2atom'], smp['norb_per'],
                                                  tip_pos_xy, smp_pos, pairs, Hb, Sb, tip_offset=np.array([x,y,args.z2d]))
                    if args.negf_mode == 'response':
                        I2d[ix,iy] = negf_response_at_energy(E_single, Hc, Sc, leadL, leadR, gammaL=args.gammaL, gammaR=args.gammaR, eta=args.eta, source=args.source, measure=args.measure)
                    else:
                        if args.solver == 'iter':
                            I2d[ix,iy] = negf_current_iterative(np.array([E_single]), Hc, Sc, leadL, leadR, gammaL=args.gammaL, gammaR=args.gammaR, eta=args.eta,
                                                               E0=E_single, E1=E_single, nrand=args.nrand, tol=args.tol, maxiter=args.maxiter, seed=0, reuse=True)[1][0]
                        else:
                            I2d[ix,iy] = negf_current(np.array([E_single]), Hc, Sc, leadL, leadR, gammaL=args.gammaL, gammaR=args.gammaR, eta=args.eta,
                                                     E0=E_single, E1=E_single)[1][0]
        print(f"  2D scan done {time.time()-t0:.1f}s  I=[{I2d.min():.2e},{I2d.max():.2e}]")
        results['T6_2dscan'] = 'PASSED'
        np.savez(os.path.join(export_dir,"scan2d.npz"), x=x_arr, y=y_arr, I=I2d, z=args.z2d, method=args.method, solver=args.solver, negf_mode=args.negf_mode, E=E_single)

        fig4,axes4 = plt.subplots(1,2,figsize=(14,6))
        fig4.suptitle(f"2D STM Image  z={args.z2d}Å  tip={args.tip}  {args.method}/{args.negf_mode}  E={E_single:.3f}eV")
        ext=[x_arr[0],x_arr[-1],y_arr[0],y_arr[-1]]
        im4a = axes4[0].imshow(I2d.T, origin='lower', extent=ext, cmap='hot', aspect='equal')
        plt.colorbar(im4a,ax=axes4[0]); axes4[0].set_title("I(x,y)")
        plot_atoms(axes4[0], smp_pos, smp_types, color='green', ms=3)
        axes4[0].set_xlabel("x(Å)"); axes4[0].set_ylabel("y(Å)")
        I2d_log = np.log10(I2d+I2d[I2d>0].min()*0.01)
        im4b = axes4[1].imshow(I2d_log.T, origin='lower', extent=ext, cmap='RdYlBu_r', aspect='equal')
        plt.colorbar(im4b,ax=axes4[1]); axes4[1].set_title("log₁₀ I(x,y)")
        axes4[1].set_xlabel("x(Å)"); axes4[1].set_ylabel("y(Å)")
        plt.tight_layout(); save_plot(fig4, f"scan2d_{args.negf_mode}", export_dir, dpi=140)
        print(f"  → {os.path.join(export_dir, f'scan2d_{args.negf_mode}.png')}")
    else:
        results['T6_2dscan'] = 'SKIPPED'; print("\nT6: Skipped (--no2d)")

    # ═══ T7: PDOS real-space projection ══════════════════════════════════════
    print("\n" + "="*65)
    print("T7: PDOS real-space projection at E_F (GridProjector)")
    print("="*65)
    try:
        FDATA_BASIS = os.path.join(_THIS_DIR, "Fdata", "basis")
        if not os.path.isdir(FDATA_BASIS):
            raise FileNotFoundError(f"Basis dir not found: {FDATA_BASIS}")

        dens3d, gspec = project_pdos_to_grid(A_smp, E_grid, E_F, smp_types, smp_pos, smp['orb2atom'].astype(int),
                                              FDATA_BASIS, dE_window=0.3, step=0.2, margin=3.5)
        if dens3d is not None and dens3d.max()>0:
            print(f"  Grid shape={dens3d.shape}  range=[{dens3d.min():.3e},{dens3d.max():.3e}]")
            np.savez(os.path.join(export_dir,"pdos_grid.npz"),dens=dens3d,origin=gspec['origin'],step=0.2)

            proj_xy = dens3d.max(axis=2); proj_xz = dens3d[:,dens3d.shape[1]//2,:]
            fig5,axes5 = plt.subplots(1,2,figsize=(14,6))
            fig5.suptitle(f"PDOS@E_F±0.3eV — real-space (PTCDA)")
            step=0.2; pmin=gspec['origin']; ngrid=gspec['ngrid']
            ext_xy=[pmin[0],pmin[0]+ngrid[0]*step, pmin[1],pmin[1]+ngrid[1]*step]
            im5a = axes5[0].imshow(proj_xy.T, origin='lower', extent=ext_xy, cmap='hot')
            plot_atoms(axes5[0], smp_pos, smp_types, color='green', ms=2)
            axes5[0].set_title("Max-proj along z"); axes5[0].set_xlabel("x(Å)"); axes5[0].set_ylabel("y(Å)")
            plt.colorbar(im5a,ax=axes5[0])
            ext_xz=[pmin[0],pmin[0]+ngrid[0]*step, pmin[2],pmin[2]+ngrid[2]*step]
            im5b = axes5[1].imshow(proj_xz.T, origin='lower', extent=ext_xz, cmap='hot')
            axes5[1].set_title("x-z slice (y=centre)"); axes5[1].set_xlabel("x(Å)"); axes5[1].set_ylabel("z(Å)")
            plt.colorbar(im5b,ax=axes5[1])
            plt.tight_layout(); save_plot(fig5, "pdos_projection", export_dir, dpi=150)
            results['T7_pdos'] = 'PASSED'
        else:
            print("  WARNING: GridProjector returned empty grid"); results['T7_pdos'] = 'FAILED'
    except Exception as e:
        import traceback; traceback.print_exc()
        print(f"  T7 ERROR: {e}"); results['T7_pdos'] = f'ERROR'

    # ═══ SUMMARY ═════════════════════════════════════════════════════════════
    print("\n" + "="*65)
    print("TEST SUMMARY")
    print("="*65)
    npass=nfail=0
    for k,v in results.items():
        icon = "✓" if v=='PASSED' else ("✗" if v=='FAILED' else "~")
        print(f"  {icon}  {k:30s}  {v}")
        if v=='PASSED': npass+=1
        elif v=='FAILED': nfail+=1
    print(f"\n  {npass} passed  |  {nfail} failed  |  {len(results)-npass-nfail} other")
    print(f"  Output files: {export_dir}/")

    if not args.noshow:
        matplotlib.use('TkAgg'); plt.show()

if __name__=="__main__": main()
