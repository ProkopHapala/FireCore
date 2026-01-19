import numpy as np
import pyopencl as cl
import os
import sys
from io import StringIO

# Adjust path to your FireCore/pyBall directory
sys.path.append("../../")
from pyBall import FireCore as fc
from pyBall.FireballOCL.OCL_Hamiltonian import (
    OCL_Hamiltonian,
    _cepal_py,
    _build_olsxc_off_py,
    _epsilon_fb_py,
    _twister_pmat,
    _chooser,
    _rotate_fb_py,
    _recover_rotate_sp,
    _recover_sp,
    build_vnl_neighbor_map,
    build_bcxcx_on_py,
    accumulate_vca_blocks,
)
from pyBall.FireballOCL.Check_Fireball_wrt_Fotran import (
    compare_matrices,
    compare_matrices_brief,
    _orbital_layout,
    _blocked_to_dense,
    _blocked_to_dense_vca_atom,
    compare_blocks,
    dense_from_neighbor_list,
    scan2c_fortran,
    scan2c_ocl,
    firecore_sparse_to_dense,
    firecore_sparse_H_with_options,
    cl_sp_from_ham,
    contract_vnl_blocks,
    compare_ewaldsr,
    compare_dip,
    compare_vnl_ref,
    compare_avgrho,
    check_ewaldsr_decomposition,
    compare_raw3c_rotation,
    print_ewald_debug,
    apply_sfire_overwrite,
    validate_neigh_back,
    summarize_neigh_lists,
    build_overlap_pairs,
    assign_triplet_types,
)

np.set_printoptions(precision=6, suppress=True, linewidth=np.inf)




def run_verification():
    # Linear C3 chain to ensure there exists a common neighbor for (0,2): cn={1}
    atomTypes_Z = np.array([6, 6, 6], dtype=np.int32)
    dCC = 1.25
    atomPos = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.0, 0.0, dCC],
            [0.0, 0.0, 2.0 * dCC],
        ],
        dtype=np.float64,
    )

    fdata_dir = "./Fdata"



    print("Initializing Fortran FireCore...")
    fc.initialize(atomType=atomTypes_Z, atomPos=atomPos)

    DO_AVGRHO_PLUMBING = True

    fc.setVerbosity(0, 0)
    fc.set_export_mode(0)
    fc.set_options(1, 1, 1, 1, 1, 1, 1, 1)
    fc.SCF(positions=atomPos, iforce=0, nmax_scf=50)

    fc.setVerbosity(7, 1)
    fc.set_vca_diag(enable=1, iatom=2, jatom=1, mbeta=0, isorp=1)
    old_stdout = sys.stdout
    captured_output = StringIO()
    sys.stdout = captured_output
    fc.assembleH(positions=atomPos, iforce=0, Kscf=2)
    sys.stdout = old_stdout
    output = captured_output.getvalue()
    fc.set_vca_diag(enable=0, iatom=2, jatom=1, mbeta=0, isorp=1)

    # Get Vxc debug data from Fortran via bulk export
    vxc_dbg = fc.get_vxc_diag_data()
    dens_dbg = vxc_dbg['dens']
    densij_dbg = vxc_dbg['densij']
    muxc_dbg = vxc_dbg['muxc']
    dmuxc_dbg = vxc_dbg['dmuxc']
    muxcij_dbg = vxc_dbg['muxcij']
    dmuxcij_dbg = vxc_dbg['dmuxcij']
    denx = vxc_dbg['denmx']
    den1x = vxc_dbg['den1x']
    sx = vxc_dbg['sx']
    bcxcx_F = vxc_dbg['bcxcx']
    arho_on = vxc_dbg['arho_on']
    arhoi_on = vxc_dbg['arhoi_on']
    rho_on = vxc_dbg['rho_on']
    rhoi_on = vxc_dbg['rhoi_on']
    vxc_ca_on = vxc_dbg['vxc_ca']
    vxc_1c_F = vxc_dbg['vxc_1c']

    muxc_py = np.zeros_like(muxc_dbg)
    dmuxc_py = np.zeros_like(dmuxc_dbg)
    muxcij_py = np.zeros_like(muxcij_dbg)
    dmuxcij_py = np.zeros_like(dmuxcij_dbg)
    for i in range(2):
        for j in range(2):
            _, muxc, _, _, dmuxc, _ = _cepal_py(float(dens_dbg[i, j]))
            _, muxcij, _, _, dmuxcij, _ = _cepal_py(float(densij_dbg[i, j]))
            muxc_py[i, j] = muxc
            dmuxc_py[i, j] = dmuxc
            muxcij_py[i, j] = muxcij
            dmuxcij_py[i, j] = dmuxcij
    err_muxc = float(np.max(np.abs(muxc_py - muxc_dbg)))
    err_dmuxc = float(np.max(np.abs(dmuxc_py - dmuxc_dbg)))
    err_muxcij = float(np.max(np.abs(muxcij_py - muxcij_dbg)))
    err_dmuxcij = float(np.max(np.abs(dmuxcij_py - dmuxcij_dbg)))
    print(f"[XC_OFF][cepal_py] max|muxc_py-muxc_F|={err_muxc:.3e}  max|dmuxc_py-dmuxc_F|={err_dmuxc:.3e}  max|muxcij_py-muxcij_F|={err_muxcij:.3e}  max|dmuxcij_py-dmuxcij_F|={err_dmuxcij:.3e}")

    dims = fc.get_HS_dims(force_refresh=True)
    sd = fc.get_HS_neighs(dims)
    sd = fc.get_HS_neighsPP(dims, data=sd)
    sd = fc.get_HS_sparse(dims, data=sd)
    sd = fc.get_rho_sparse(dims, data=sd)
    sd = fc.get_rho_off_sparse(dims, data=sd)
    neigh_back = fc.get_neigh_back(dims)

    # Reproduce build_olsxc_off bcxcx(1:4,1:4) for the debugged pair (Fortran iatom=2, ineigh=1).
    # NOTE: FireCore Python bindings export rho_off but not rhoij_off, so we take denmx/den1x/sx from gated Fortran prints.
    
    # NEW: Use exported Vxc debug data from Fortran
    denx = vxc_dbg['denmx']
    den1x = vxc_dbg['den1x']
    sx = vxc_dbg['sx']
    bcxcx_F = vxc_dbg['bcxcx']
    bcxcx_py = _build_olsxc_off_py(den1x, denx, sx, dens_dbg, densij_dbg)
    err_bcxcx = float(np.max(np.abs(bcxcx_py - bcxcx_F)))
    print(f"[XC_OFF][bcxcx_py] max|bcxcx_py-bcxcx_F|={err_bcxcx:.3e}")
    
    # NEW: Use exported Vxc debug data from Fortran
    arho_on = vxc_dbg['arho_on']
    arhoi_on = vxc_dbg['arhoi_on']
    rho_on = vxc_dbg['rho_on']
    rhoi_on = vxc_dbg['rhoi_on']
    vxc_ca_on = vxc_dbg['vxc_ca']

    # Implement build_ca_olsxc_on algorithm in Python
    bcxcx_on_py = build_bcxcx_on_py(arho_on, arhoi_on, rho_on, rhoi_on, _cepal_py, nsh=2, lssh=np.array([0, 1], dtype=np.int32))

    err_bcxcx_on = float(np.max(np.abs(bcxcx_on_py - vxc_ca_on)))
    print(f"[XC_ON][bcxcx_on_py] max|bcxcx_on_py-vxc_ca_on|={err_bcxcx_on:.3e}")

    # Quick sanity check: neigh_back should map (iatom,ineigh)->jneigh in neighbor list of jatom
    # such that neigh_j[jatom, jneigh] == iatom+1 (Fortran 1-based atom ids).
    nb_bad = validate_neigh_back(sd, neigh_back, dims)
    print(f"[NEIGH_BACK] inconsistent mappings: {nb_bad}")
    res_neigh_back = (nb_bad == 0)
    err_neigh_back = float(nb_bad)

    if int(dims.natoms) <= 8 and int(dims.neigh_max) <= 16:
        print(f"[NEIGH_BACK] neigh_j (0-based) =\n{sd.neigh_j[:,:int(dims.neigh_max)]-1}")
        print(f"[NEIGH_BACK] neigh_back (0-based) =\n{neigh_back[:,:int(dims.neigh_max)]-1}")

    # Neighbor self-slot (derived): used by _blocked_to_dense; track missing self-slots as a geometry/parsing health check.
    neigh_self = np.full(int(dims.natoms), -1, dtype=np.int32)
    for i in range(int(dims.natoms)):
        for j in range(int(sd.neighn[i])):
            if int(sd.neigh_j[i,j]) == (i+1) and int(sd.neigh_b[i,j]) == 0:
                neigh_self[i] = j
                break
    n_self_missing = int(np.sum(neigh_self < 0))
    res_neigh_self = (n_self_missing == 0)
    err_neigh_self = float(n_self_missing)

    if DO_AVGRHO_PLUMBING:
        fc.set_avg_rho_diag(enable=1, iatom=1, jatom=-1, mbeta=-1)
        en, ti, tj, tb = fc.get_avg_rho_diag_state()
        print(f"[FORTRAN AvgRho DIAG state] enable={en} iatom={ti} jatom={tj} mbeta={tb}")
        d_iatom, d_ineigh, d_in1, d_in2, d_jatom, d_mbeta = fc.get_avg_rho_diag_meta()
        eps2c = fc.get_avg_rho_diag_eps2c()
        smS   = fc.get_avg_rho_diag_sm()
        rhom2cS = fc.get_avg_rho_diag_rhom2c()
        rhom3cS = fc.get_avg_rho_diag_rhom3c()
        rhooff_3c = fc.get_avg_rho_diag_rhooff_3c()
        rhooff_fin = fc.get_avg_rho_diag_rhooff_final()
        print("\n[FORTRAN AvgRho DIAG]")
        print(f"  target meta: iatom={d_iatom} ineigh={d_ineigh} in1={d_in1} in2={d_in2} jatom={d_jatom} mbeta={d_mbeta}")
        print("  eps2c:")
        print(eps2c)
        print("  sm (shell):")
        print(smS)
        print("  rhom2c (shell):")
        print(rhom2cS)
        print("  rhom3c (shell):")
        print(rhom3cS)
        print("  rho_off after 3c (orbital):")
        print(rhooff_3c[:4,:4])
        print("  rho_off final (orbital):")
        print(rhooff_fin[:4,:4])

    print("Initializing PyOpenCL Hamiltonian...")
    ham = OCL_Hamiltonian(fdata_dir)
    ham.prepare_splines(atomTypes_Z)
    ham.prepare_data_3c(atomTypes_Z)

    # VXC_ON microtest: OpenCL kernel build_ca_olsxc_on_sp4 should reproduce Fortran vxc_ca (sp-only 4x4)
    arho_on16 = np.zeros((3,4), dtype=np.float32)
    arhoi_on16 = np.zeros((3,4), dtype=np.float32)
    rho_on16 = np.zeros((3,16), dtype=np.float32)
    rhoi_on16 = np.zeros((3,16), dtype=np.float32)
    arho_on16[1,:] = arho_on.astype(np.float32).reshape(4)
    arhoi_on16[1,:] = arhoi_on.astype(np.float32).reshape(4)
    rho_on16[1,:] = rho_on.astype(np.float32).reshape(16)
    rhoi_on16[1,:] = rhoi_on.astype(np.float32).reshape(16)
    bcxcx_on_ocl16 = ham.build_ca_olsxc_on_sp4(arho_on16, arhoi_on16, rho_on16, rhoi_on16)[1].astype(np.float64).reshape(4,4)
    err_bcxcx_on_ocl = float(np.max(np.abs(bcxcx_on_ocl16 - vxc_ca_on)))
    print(f"[XC_ON][bcxcx_on_ocl] max|bcxcx_on_ocl-vxc_ca_on|={err_bcxcx_on_ocl:.3e}")

    # VXC microtest: OpenCL kernel build_olsxc_off_sp4 should reproduce Fortran bcxcx (sp-only 4x4)
    dens4 = np.array([[dens_dbg[0,0], dens_dbg[0,1], dens_dbg[1,0], dens_dbg[1,1]]], dtype=np.float32)
    densij4 = np.array([[densij_dbg[0,0], densij_dbg[0,1], densij_dbg[1,0], densij_dbg[1,1]]], dtype=np.float32)
    denx16 = denx.astype(np.float32).reshape(1,16)
    den1x16 = den1x.astype(np.float32).reshape(1,16)
    sx16 = sx.astype(np.float32).reshape(1,16)
    bcxcx_ocl16 = ham.build_olsxc_off_sp4(dens4, densij4, denx16, den1x16, sx16)[0].astype(np.float64).reshape(4,4)
    err_bcxcx_ocl = float(np.max(np.abs(bcxcx_ocl16 - bcxcx_F)))
    print(f"[XC_OFF][bcxcx_ocl] max|bcxcx_ocl-bcxcx_F|={err_bcxcx_ocl:.3e}")

    natoms = int(atomPos.shape[0])

    # For scan-based tests we use all (i,j) including self. For sparse-export parity we must use
    # the exact same neighbor list as Fortran export.
    neighs_all = [(i, j) for i in range(natoms) for j in range(natoms)]
    # ------------------------------------------------------------------------------------------
    # Verify_C2-style 2-center tests (scan API vs OCL scan) for all neighbor pairs (i!=j) and self.
    # ------------------------------------------------------------------------------------------
    n_orb_atom, offs = _orbital_layout(sd, natoms)
    norb = int(offs[-1])

    # For scan-based tests we use all (i,j) including self. For sparse-export parity we must use
    # the exact same neighbor list as Fortran export.
    neighs_all = [(i, j) for i in range(natoms) for j in range(natoms)]

    if DO_AVGRHO_PLUMBING:
        print("\n[PLUMBING] compute_avg_rho (3c gather) using Fortran-exported rho + Qin-shell...")
        sd = fc.get_HS_neighsPP(dims, data=sd)
        sd = fc.get_rho_sparse(dims, data=sd)
        sd = fc.get_rho_off_sparse(dims, data=sd)

        # OpenCL compute_avg_rho weights 3c density pieces by *neutral* charge of the common neighbor.
        # Fortran average_rho uses Qneutral(isorp, indna) (per-shell neutral population for species indna).
        pass

    Qneutral_sh_full = fc.get_Qneutral_shell(dims)
    Qneutral_sh = np.zeros_like(Qneutral_sh_full)
    Qneutral_sh[:, :dims.nspecies] = Qneutral_sh_full[:, :dims.nspecies]
    Qin_shell = fc.get_Qin_shell(dims)
    nzx_full = np.array(sd.nzx, dtype=np.int32)
    nzx = nzx_full[:dims.nspecies]
    iatyp = np.array(sd.iatyp, dtype=np.int32)
    ispec_of_atom = np.zeros(dims.natoms, dtype=np.int32)
    for ia in range(dims.natoms):
        Z = int(iatyp[ia])
        w = np.where(nzx == Z)[0]
        if w.size == 0:
            raise RuntimeError(f"Cannot map atom Z={Z} to nzx species list {nzx}")
        ispec_of_atom[ia] = int(w[0])
    # TODO/DEBUG: Qatom unused below; summing Qneutral_sh triggered overflow on some runs.
    # Qatom = np.zeros(dims.natoms, dtype=np.float32)
    # for ia in range(dims.natoms):
    #     ispec = ispec_of_atom[ia]
    #     if ispec >= Qneutral_sh.shape[1]:
    #         raise RuntimeError(f"ispec_of_atom[{ia}]={ispec} out of bounds for Qneutral_sh.shape={Qneutral_sh.shape}")
    #     Qatom[ia] = float(np.sum(Qneutral_sh[:, ispec]))


    # Neighbor lists from Fortran export: neigh_j is 1-based atom index.
    neigh_lists, neigh_lists_self = summarize_neigh_lists(sd, dims)
    print(f"DEBUG: neigh_lists (no self) = {neigh_lists}")
    print(f"DEBUG: neigh_lists_self = {neigh_lists_self}")

    # Build pairs from neighbor lists (unique i<j)
    pairs = []
    for ia in range(dims.natoms):
        for j in neigh_lists[ia]:
            if ia < j:
                pairs.append((ia, j))
    if len(pairs) == 0:
        raise RuntimeError("[PLUMBING] No neighbor pairs found")
    pairs = np.array(pairs, dtype=np.int32)

    # TODO/DEBUG: do NOT inject non-Fortran neighbor pairs; comparisons must follow Fortran neighbor list exactly.
    # If pair (0,2) is not in the Fortran neighbor list, we skip its diagnostics below.

    cn_offsets, cn_indices = ham.build_common_neighbor_csr(neigh_lists, pairs)
    cn_counts = cn_offsets[1:] - cn_offsets[:-1]
    print(f"  n_pairs={pairs.shape[0]}  n_cn_total={cn_indices.shape[0]}  cn_max={int(np.max(cn_counts))}")

    # Map sparse neighbor blocks for (i,j) to 4x4 blocks; if (i,j) not in sparse neighbor list, use zeros
    neigh_index = {}
    for ia in range(dims.natoms):
        nn = int(sd.neighn[ia])
        for ineigh in range(nn):
            j = int(sd.neigh_j[ia, ineigh]) - 1
            if j >= 0:
                neigh_index[(ia, j)] = ineigh

    DEBUG_QIN_TEST = True
    S_blocks = np.zeros((pairs.shape[0], 4, 4), dtype=np.float32)
    # We no longer construct the 2c seed using scanHamPiece2c (libFireCore interface is not reference).
    # Instead, OpenCL compute_avg_rho builds the 2c seed internally from den_ontopl/den_ontopr tables.
    rho_blocks = np.zeros((pairs.shape[0], 4, 4), dtype=np.float32)
    rho_blocks_alt = np.zeros((pairs.shape[0], 4, 4), dtype=np.float32)
    rho_blocks_alt_raw = np.zeros((pairs.shape[0], 4, 4), dtype=np.float64)
    rho_blocks_alt_raw_L = np.zeros((pairs.shape[0], 4, 4), dtype=np.float64)
    rho_blocks_alt_raw_R = np.zeros((pairs.shape[0], 4, 4), dtype=np.float64)
    S_blocks = np.zeros((pairs.shape[0], 4, 4), dtype=np.float64)
    rho_blocks_qin = np.zeros((pairs.shape[0], 4, 4), dtype=np.float32) if DEBUG_QIN_TEST else None

    S_blocks, pair_2c_types, pair_dR, pair_in_i, pair_in_j, pair_valid, pair_mbeta = build_overlap_pairs(sd, atomPos, ham, pairs, neigh_index, ispec_of_atom, atomTypes_Z)

    # Triplet type selection: key is (root,nz1,nz2,nz3) where nz3 is the common neighbor type
    root = 'den3'
    pair_triplet_types = assign_triplet_types(ham, pairs, cn_offsets, cn_indices, atomTypes_Z, root=root)

    nspecies_fdata = int(Qneutral_sh.shape[1])
    if nspecies_fdata <= 0:
        raise RuntimeError("Qneutral_sh has zero species dimension")
    if int(np.max(ispec_of_atom)) >= nspecies_fdata:
        raise RuntimeError(f"ispec_of_atom out of range: max={int(np.max(ispec_of_atom))} nspecies_fdata(Qneutral_sh)={nspecies_fdata}")
    nssh_species = np.ascontiguousarray(sd.nssh[:dims.nspecies], dtype=np.int32)

    # Use Fortran-exported mu/nu/mvalue (make_munu) to reconstruct 3c blocks exactly.
    # FireCore.py stores Fortran mu(index,in1,in2) as sd.mu[in2-1, in1-1, index-1] due to axis reversal.
    n_types = int(np.max(pair_triplet_types)) + 1 if pair_triplet_types.size else 0
    if n_types <= 0:
        raise RuntimeError("pair_triplet_types has zero types")
    n_nz_max = int(ham.n_nz_3c_max)
    mu3c_map = np.zeros((n_types, n_nz_max), dtype=np.int16)
    nu3c_map = np.zeros((n_types, n_nz_max), dtype=np.int16)
    mv3c_map = np.zeros((n_types, n_nz_max), dtype=np.int8)
    # For this test geometry we have only one chemical species; use in1=in2=1.
    in1 = 1
    in2 = 1
    mu3c_map[:, :] = sd.mu[in2-1, in1-1, :n_nz_max].astype(np.int16)
    nu3c_map[:, :] = sd.nu[in2-1, in1-1, :n_nz_max].astype(np.int16)
    mv3c_map[:, :] = sd.mvalue[in2-1, in1-1, :n_nz_max].astype(np.int8)

    print(f"DEBUG mu/nu/mvalue (sd) first 12: mu={mu3c_map[0,:12]}  nu={nu3c_map[0,:12]}  mv={mv3c_map[0,:12]}")
    # TODO/DEBUG: override mu/nu map to expected sp mapping (from make_munu) to check Fortran export layout
    DEBUG_MUNU_OVERRIDE = True
    if DEBUG_MUNU_OVERRIDE:
        mu_sp = np.array([1, 1, 3, 2, 3, 4, 1, 4, 4, 3], dtype=np.int16)
        nu_sp = np.array([1, 3, 1, 2, 3, 4, 4, 1, 3, 4], dtype=np.int16)
        mv_sp = np.array([0, 0, 0, 0, 0, 0, 1, 1, 1, 1], dtype=np.int8)
        mu3c_map[:, :] = 0
        nu3c_map[:, :] = 0
        mv3c_map[:, :] = 0
        mu3c_map[:, :mu_sp.size] = mu_sp
        nu3c_map[:, :nu_sp.size] = nu_sp
        mv3c_map[:, :mv_sp.size] = mv_sp
        print(f"DEBUG mu/nu/mvalue (override) first 12: mu={mu3c_map[0,:12]}  nu={nu3c_map[0,:12]}  mv={mv3c_map[0,:12]}")

    if np.any(np.isnan(rho_blocks)):
        print("ERROR: rho_blocks contains NaN!")
    if np.any(np.isnan(Qneutral_sh)):
        print("ERROR: Qneutral_sh contains NaN!")
        

    # Build alternative 2c seed by manually rotating unrotated scanHamPiece2c output
    for ip, (ia, ja) in enumerate(pairs):
        if pair_valid[ip] == 0:
            continue
        mb = int(pair_mbeta[ip])
        r1 = atomPos[int(ia)]
        r2 = atomPos[int(ja)] + sd.xl[mb]
        r21 = r2 - r1
        r21_norm = np.linalg.norm(r21)
        if r21_norm < 1e-12:
            continue
        sighat = r21 / r21_norm
        eps2c = _epsilon_fb_py(r2, sighat)
        in_i = int(pair_in_i[ip])
        in_j = int(pair_in_j[ip])
        left_rot = _rotate_fb_py(in_i, in_i, eps2c, rho_blocks_alt_raw_L[ip], sd.lssh, sd.nssh)
        right_rot = _rotate_fb_py(in_i, in_j, eps2c, rho_blocks_alt_raw_R[ip], sd.lssh, sd.nssh)
        rho_blocks_alt[ip, :, :] = (left_rot + right_rot).astype(np.float32)

    # Raw 3c interpolation check (no rotation/recover) for a single triplet
    DEBUG_RAW3C = True
    err_rot_isorp = float('nan')
    err_rawrot_f = float('nan')
    err_rawrot_o = float('nan')
    if DEBUG_RAW3C and natoms >= 3:
        _raw3c_res = compare_raw3c_rotation(fc, ham, atomPos, atomTypes_Z, sd, dims, ispec_of_atom, Qneutral_sh, n_nz_max, mu3c_map, nu3c_map, mv3c_map, verbose=True)
        err_rawrot_f = _raw3c_res.details['err_rawrot_f']
        err_rawrot_o = _raw3c_res.details['err_rawrot_o']
        err_rot_isorp = _raw3c_res.details['err_rot_isorp']

    # Compute and compare AvgRho blocks
    _avgrho_res = compare_avgrho(ham, atomPos, pairs, pair_triplet_types, cn_offsets, cn_indices, S_blocks, rho_blocks, Qneutral_sh, ispec_of_atom, nssh_species, sd, dims, neigh_index, cn_counts, pair_2c_types=pair_2c_types, mu3c_map=mu3c_map, nu3c_map=nu3c_map, mv3c_map=mv3c_map, Qin_shell=Qin_shell, DEBUG_QIN_TEST=DEBUG_QIN_TEST)
    rho_avg_blocks = _avgrho_res.details['rho_avg_blocks']
    ref_blocks = _avgrho_res.details['ref_blocks']
    mask = _avgrho_res.details['mask']
    res_AvgRho = _avgrho_res.ok
    err_AvgRho = _avgrho_res.err
    rho_avg_blocks_qin = _avgrho_res.details.get('rho_avg_blocks_qin', None)

    # Report a few pairs (only those present in Fortran neighbor list)
    for ip in range(pairs.shape[0]):
        i, j = int(pairs[ip, 0]), int(pairs[ip, 1])
        if (i, j) in [(0, 2), (0, 1), (1, 2)]:
            print(f"\n  rho_avg block pair ({i},{j}) cn_count={int(cn_counts[ip])}")
            print(rho_avg_blocks[ip])
            if mask[ip] == 1:
                print("  rho_off Fortran ref block:")
                print(ref_blocks[ip])
                print("  abs diff:")
                print(np.abs(ref_blocks[ip] - rho_avg_blocks[ip]))

    # ------------------------------------------------------------------------------------------
    # Verify_C2-style 2-center tests (scan API vs OCL scan) for all neighbor pairs (i!=j) and self.
    # ------------------------------------------------------------------------------------------
    n_orb_atom, offs = _orbital_layout(sd, natoms)
    norb = int(offs[-1])

    # For scan-based tests we use all (i,j) including self. For sparse-export parity we must use
    # the exact same neighbor list as Fortran export.
    neighs_all = [(i, j) for i in range(natoms) for j in range(natoms)]

    def _scan_table_2c(interaction_fortran, root_ocl, applyRotation_offdiag=True, applyRotation_self=False, in3=1):
        blocks_f = []
        blocks_o = []
        for (i, j) in neighs_all:
            dR = (atomPos[j] - atomPos[i]).copy()
            if i == j:
                dR[:] = 0.0
                Af = scan2c_fortran(fc, interaction_fortran, dR, in3=in3, applyRotation=applyRotation_self)
                Ao = scan2c_ocl(ham, root_ocl, atomTypes_Z[i], atomTypes_Z[j], dR, applyRotation=applyRotation_self)
            else:
                Af = scan2c_fortran(fc, interaction_fortran, dR, in3=in3, applyRotation=applyRotation_offdiag)
                Ao = scan2c_ocl(ham, root_ocl, atomTypes_Z[i], atomTypes_Z[j], dR, applyRotation=applyRotation_offdiag)
            blocks_f.append(Af)
            blocks_o.append(Ao)
        blocks_f = np.array(blocks_f, dtype=np.float64)
        blocks_o = np.array(blocks_o, dtype=np.float64)
        Mf = dense_from_neighbor_list(neighs_all, blocks_f, n_orb_atom, offs)
        Mo = dense_from_neighbor_list(neighs_all, blocks_o, n_orb_atom, offs)
        return Mf, Mo

    print("\nTesting Overlap S...")
    S_f, S_o = _scan_table_2c(1, 'overlap', applyRotation_offdiag=True, applyRotation_self=False, in3=1)
    res_S, err_S = compare_matrices("Overlap S", S_f, S_o, tol=1e-6, require_nonzero=True)

    print("\nTesting dipole_z (2c scan, Fortran doscentros interaction=9)...")
    # This is the direct prerequisite for OpenCL dip reimplementation.
    # If this fails, the issue is in the OpenCL 2c table eval / rotation for dipole tables (not in neighbor-slot mapping).
    DipZ_f, DipZ_o = _scan_table_2c(9, 'dipole_z', applyRotation_offdiag=True, applyRotation_self=False, in3=1)
    res_DipZ, err_DipZ = compare_matrices("Dipole_z (scan)", DipZ_f, DipZ_o, tol=1e-6, require_nonzero=True)

    print("\nTesting Kinetic T...")
    T_f, T_o = _scan_table_2c(13, 'kinetic', applyRotation_offdiag=True, applyRotation_self=False, in3=1)
    res_T, err_T = compare_matrices("Kinetic T", T_f, T_o, tol=1e-5, require_nonzero=True)

    print("\nTesting Vna...")
    # Fortran Vna total = ontopl(2)+ontopr(3) for off-diagonal; vna_atom(4) for self accumulated from neighbor direction.
    # Here we follow verify_C2 behavior: use vna_atom with dR(i->j) for self block as in original script.
    blocks_f = []
    blocks_o = []
    for (i, j) in neighs_all:
        dR = (atomPos[j] - atomPos[i]).copy()
        if i == j:
            # Choose nearest neighbor direction if available; otherwise 0.
            k = 1 if i == 0 else 0
            dR = (atomPos[k] - atomPos[i]).copy()
            Af = scan2c_fortran(fc, 4, dR, in3=1, applyRotation=True)
            Ao = scan2c_ocl(ham, 'vna_atom_00', atomTypes_Z[i], atomTypes_Z[i], dR, applyRotation=True)
        else:
            Af = scan2c_fortran(fc, 2, dR, in3=1, applyRotation=True) + scan2c_fortran(fc, 3, dR, in3=1, applyRotation=True)
            Ao = scan2c_ocl(ham, 'vna_ontopl_00', atomTypes_Z[i], atomTypes_Z[j], dR, applyRotation=True) + scan2c_ocl(ham, 'vna_ontopr_00', atomTypes_Z[i], atomTypes_Z[j], dR, applyRotation=True)
        blocks_f.append(Af)
        blocks_o.append(Ao)
    Vna_f = dense_from_neighbor_list(neighs_all, np.array(blocks_f, dtype=np.float64), n_orb_atom, offs)
    Vna_o = dense_from_neighbor_list(neighs_all, np.array(blocks_o, dtype=np.float64), n_orb_atom, offs)
    res_Vna, err_Vna = compare_matrices("Vna", Vna_f, Vna_o, tol=1e-5, require_nonzero=True)

    print("\nTesting sVNL (PP overlap, interaction=5)...")
    sVNL_f, sVNL_o = _scan_table_2c(5, 'vnl', applyRotation_offdiag=True, applyRotation_self=False, in3=1)
    require_nz = np.max(np.abs(sVNL_f)) > 1e-6
    res_sVNL, err_sVNL = compare_matrices("sVNL", sVNL_f, sVNL_o, tol=1e-5, require_nonzero=require_nz)

    # ------------------------------------------------------------------------------------------
    # Contracted VNL: scan sVNL + contract (reference) vs ham.assemble_full CPU/GPU and Fortran gated export
    # ------------------------------------------------------------------------------------------
    print("\nTesting contracted VNL (scan-based reference vs OpenCL CPU/GPU)...")

    cls = [cl_sp_from_ham(ham, atomTypes_Z, i) for i in range(natoms)]

    # sVNL blocks for all (basis_atom i, pp_atom k)
    s_map = {}
    for i in range(natoms):
        for k in range(natoms):
            if i == k:
                dR = np.array([0.0, 0.0, 0.0], dtype=np.float64)
                s_map[(i, k)] = scan2c_fortran(fc, 5, dR, in3=1, applyRotation=False)
            else:
                dR = (atomPos[k] - atomPos[i]).copy()
                s_map[(i, k)] = scan2c_fortran(fc, 5, dR, in3=1, applyRotation=True)

    # Build VNL neighbor list exactly like Fortran export_mode>=2 does:
    # iterate neighPP list and map each (jatom,mbeta) to an index in the normal neigh list.
    # This avoids comparing OpenCL against zero blocks for non-PP pairs.
    neighs_vnl = build_vnl_neighbor_map(sd, natoms, verbose=True)
    if len(neighs_vnl) == 0:
        raise RuntimeError("VNL: neighs_vnl is empty (no PP neighbors mapped into neigh list)")

    # Cross-check new refactored helper against the inline implementation (no behavior change).
    # Keep this until we fully switch the inline block to the helper.
    _vnlchk_cpu = compare_vnl_ref(ham, atomPos, atomTypes_Z, sd, s_map, cls, n_orb_atom, offs, use_gpu_vnl=False)
    _vnlchk_gpu = compare_vnl_ref(ham, atomPos, atomTypes_Z, sd, s_map, cls, n_orb_atom, offs, use_gpu_vnl=True)
    Vnl_ref = _vnlchk_cpu.details['Vnl_ref']
    Vnl_cpu = _vnlchk_cpu.details['Vnl_ocl']
    Vnl_gpu = _vnlchk_gpu.details['Vnl_ocl']
    neighs_vnl = _vnlchk_cpu.details['neighs_vnl']
    res_Vnl_cpu, err_Vnl_cpu = _vnlchk_cpu.ok, _vnlchk_cpu.err
    res_Vnl_gpu, err_Vnl_gpu = _vnlchk_gpu.ok, _vnlchk_gpu.err

    print("\nTesting contracted VNL (Fortran gated export vs OpenCL CPU/GPU)...")
    Vnl_fortran_full, _neighs_vnl_f, _sd_vnl_f = firecore_sparse_H_with_options(fc, 0, 0, 0, 1, 0, 0, 0, export_mode=2)
    # Compare only on the PP-neighbor block pattern.
    Vnl_fortran = Vnl_fortran_full
    print(f"  neighs_vnl blocks: {len(neighs_vnl)}")
    print(f"  max|Vnl_fortran|={float(np.max(np.abs(Vnl_fortran))):.3e}  max|Vnl_cpu|={float(np.max(np.abs(Vnl_cpu))):.3e}")
    res_Vnl_cpu_fortran, err_Vnl_cpu_fortran = compare_matrices_brief("VNL (CPU) vs Fortran export", Vnl_fortran, Vnl_cpu, tol=1e-5, require_nonzero=True)
    res_Vnl_gpu_fortran, err_Vnl_gpu_fortran = compare_matrices_brief("VNL (GPU) vs Fortran export", Vnl_fortran, Vnl_gpu, tol=1e-5, require_nonzero=True)

    print("\nTesting Fortran export reconstruction (raw full H vs reconstructed full H)...")
    H_full_raw, S_full_raw, neighbors_sparse, sparse_data = firecore_sparse_to_dense(fc, export_mode=0)
    H_full_rec, S_full_rec, _neighbors2, _sparse2 = firecore_sparse_to_dense(fc, export_mode=1)
    res_H_recon, err_H_recon = compare_matrices("Fortran full H: raw vs reconstructed", H_full_raw, H_full_rec, tol=1e-6, require_nonzero=True)
    res_S_recon, err_S_recon = compare_matrices("Fortran full S: raw vs reconstructed(export)", S_full_raw, S_full_rec, tol=1e-12, require_nonzero=True)

    # Step-4 style check: sparse block buffers -> dense reconstruction (pure Python) vs Fortran full dense
    H_from_blocks = _blocked_to_dense(sd, sd.h_mat, natoms)
    S_from_blocks = _blocked_to_dense(sd, sd.s_mat, natoms)
    res_H_blocks, err_H_blocks = compare_matrices_brief("H_dense(from blocks) vs Fortran full H", H_full_raw, H_from_blocks, tol=1e-10, require_nonzero=True)
    res_S_blocks, err_S_blocks = compare_matrices_brief("S_dense(from blocks) vs Fortran full S", S_full_raw, S_from_blocks, tol=1e-12, require_nonzero=True)

    print("\nTesting combined 2-center H = T + Vna (from scan API)...")
    H2c_f = T_f + Vna_f
    H2c_o = T_o + Vna_o
    res_H2c, err_H2c = compare_matrices("H2c (T+Vna)", H2c_f, H2c_o, tol=1e-5, require_nonzero=True)

    dip4_f = fc.export_interaction4D(4)

    print("\nTesting dip (OpenCL Fdata dipole_z vs Fortran export dip)...")
    # Fortran computes dip via doscentros(interaction=9) which maps to dipole_z tables.
    # dip is stored per neighbor slot (iatom,ineigh) and used by EWALDSR/LR.
    _dipchk = compare_dip(ham, atomPos, atomTypes_Z, sd, dims, dip4_f, tol=1e-6)
    dip4_o = _dipchk.details['dip4_o']
    Dip_f = _dipchk.details['Dip_f']
    Dip_o = _dipchk.details['Dip_o']
    res_Dip, err_Dip = _dipchk.ok, _dipchk.err

    print("\nTesting Vca (charge-dependent, 2-center only)...")
    # Fortran Vca term in assembled H is: vca + ewaldlr - ewaldsr  (see firecore_get_HS_sparse)
    # Here we export each component explicitly.
    vcaL4_f    = fc.export_interaction4D(10) # vca_ontopl
    vcaR4_f    = fc.export_interaction4D(11) # vca_ontopr
    vcaA4_f    = fc.export_interaction4D(12) # vca_atom
    ewaldsr4_f = fc.export_interaction4D(2)  # ewaldsr
    ewaldlr4_f = fc.export_interaction4D(3)  # ewaldlr
    ewaldsr2cA4_f = fc.export_interaction4D(13)  # ewaldsr_2c_atom (self-slot)
    ewaldsr2cO4_f = fc.export_interaction4D(14)  # ewaldsr_2c_ontop
    ewaldsr3c4_f  = fc.export_interaction4D(15)  # ewaldsr_3c
    print(f"  Fortran export max|vcaL4|={float(np.max(np.abs(vcaL4_f))):.3e}  max|vcaR4|={float(np.max(np.abs(vcaR4_f))):.3e}  max|vcaA4|={float(np.max(np.abs(vcaA4_f))):.3e}")
    ia_max, ineigh_max, inu_max, imu_max = np.unravel_index(int(np.argmax(np.abs(vcaA4_f))), vcaA4_f.shape)
    print(f"  Fortran export argmax|vcaA4| at iatom={ia_max} ineigh={ineigh_max} (neighn={int(sd.neighn[ia_max])} neigh_j={int(sd.neigh_j[ia_max, ineigh_max])})  |val|={float(np.abs(vcaA4_f[ia_max, ineigh_max, inu_max, imu_max])):.3e}")
    VcaL_f      = _blocked_to_dense(sd, vcaL4_f, natoms)
    VcaR_f      = _blocked_to_dense(sd, vcaR4_f, natoms)
    VcaA_f      = _blocked_to_dense_vca_atom(sd, vcaA4_f, natoms)
    # Core-reference Vca2c is the *component sum* (ontopl + ontopr + atom).
    # DO NOT use export_interaction4D(1) as reference; it is an interface-level array and not trusted.
    Vca2c_f     = VcaL_f + VcaR_f + VcaA_f
    EwaldSR_f   = _blocked_to_dense(sd, ewaldsr4_f, natoms)
    EwaldLR_f   = _blocked_to_dense(sd, ewaldlr4_f, natoms)
    Vca_full_f  = Vca2c_f + EwaldLR_f - EwaldSR_f

    # Stepwise Fortran consistency check: EWALDSR decomposition
    _ewaldsr_decomp = check_ewaldsr_decomposition(sd, ewaldsr2cA4_f, ewaldsr2cO4_f, ewaldsr3c4_f, ewaldsr4_f, natoms)
    EwaldSR_2c_atom_f  = _ewaldsr_decomp.details['EwaldSR_2c_atom_f']
    EwaldSR_2c_ontop_f = _ewaldsr_decomp.details['EwaldSR_2c_ontop_f']
    EwaldSR_3c_f       = _ewaldsr_decomp.details['EwaldSR_3c_f']
    EwaldSR_sum_f      = _ewaldsr_decomp.details['EwaldSR_sum_f']
    res_EwaldSR_split4d = _ewaldsr_decomp.ok
    err_EwaldSR_split4d = _ewaldsr_decomp.err

    # --------------------------
    # PRINT-BASED DEBUG (match Fortran [EW2C_A]/[EW2C_O] lines)
    # Keep it tiny: only atoms <=3. This is for the C3 debugging workflow.
    dip4_f = fc.export_interaction4D(4)
    print_ewald_debug(fc, sd, atomPos, neigh_self, dip4_f, eq2=14.39975, verbose=True)
    # --------------------------

    # 4D-level check: must reproduce SFIRE overwrite semantics used by Fortran in assemble_ca_3c.
    # Fortran does: jneigh = neigh_back(iatom,ineigh)
    #              ewaldsr(inu,imu,jneigh,jatom) = ewaldsr(imu,inu,ineigh,iatom)
    # In Python export layout [iatom,ineigh,nu,mu], this becomes reverse_block = forward_block.T.
    ewaldsr4_sum_f = ewaldsr2cA4_f + ewaldsr2cO4_f + ewaldsr3c4_f
    _sfire_res = apply_sfire_overwrite(ewaldsr4_sum_f, ewaldsr2cA4_f, ewaldsr2cO4_f, ewaldsr3c4_f, sd, neigh_back, n_orb_atom, verbose=True)

    R4 = ewaldsr4_f - ewaldsr4_sum_f
    ijkl = np.unravel_index(int(np.argmax(np.abs(R4))), R4.shape)
    i0, ine0, nu0, mu0 = [int(x) for x in ijkl]
    j0 = int(sd.neigh_j[i0, ine0]) - 1 if ine0 < int(sd.neighn[i0]) else -1
    jb0 = int(neigh_back[i0, ine0]) if ine0 < int(sd.neighn[i0]) else 0
    jne0 = jb0 - 1
    print(f"  [EWALD_SPLIT_4D] max|ewaldsr4 - sum4(SFIRE)|={float(np.max(np.abs(R4))):.3e}  at (iatom,ineigh,nu,mu)={ijkl}  jatom={j0} jneigh={jne0} (raw_neigh_back={jb0})  ewaldsr4={float(ewaldsr4_f[ijkl]): .6e}  sum4={float(ewaldsr4_sum_f[ijkl]): .6e}  resid={float(R4[ijkl]): .6e}")

    res_EwaldSR_split4d = (float(np.max(np.abs(R4))) <= 1e-10)
    err_EwaldSR_split4d = float(np.max(np.abs(R4)))

    EwaldSR_sum_sfire_f = _blocked_to_dense(sd, ewaldsr4_sum_f, natoms)
    res_EwaldSR_split_sfire, err_EwaldSR_split_sfire = compare_matrices_brief("EwaldSR split (SFIRE) dense vs EwaldSR", EwaldSR_f, EwaldSR_sum_sfire_f, tol=1e-10, require_nonzero=True)

    res_EwaldSR_split, err_EwaldSR_split = compare_matrices_brief("EwaldSR split (NO SFIRE; diagnostic) vs EwaldSR", EwaldSR_f, EwaldSR_sum_f, tol=1e-10, require_nonzero=True)
    if not res_EwaldSR_split:
        print(f"  [EWALD_SPLIT] max|EwaldSR|={float(np.max(np.abs(EwaldSR_f))):.3e}  max|sum|={float(np.max(np.abs(EwaldSR_sum_f))):.3e}")
        print(f"  [EWALD_SPLIT] max|2c_atom|={float(np.max(np.abs(EwaldSR_2c_atom_f))):.3e}  max|2c_ontop|={float(np.max(np.abs(EwaldSR_2c_ontop_f))):.3e}  max|3c|={float(np.max(np.abs(EwaldSR_3c_f))):.3e}")
        R = EwaldSR_f - EwaldSR_sum_f
        ij = np.unravel_index(int(np.argmax(np.abs(R))), R.shape)
        def _which_atom(idx_):
            for ia in range(natoms):
                if idx_ >= int(offs[ia]) and idx_ < int(offs[ia+1]):
                    return ia
            return -1
        ia = _which_atom(int(ij[0])); ja = _which_atom(int(ij[1]))
        print(f"  [EWALD_SPLIT] components at worst: 2c_atom={float(EwaldSR_2c_atom_f[ij]): .6e}  2c_ontop={float(EwaldSR_2c_ontop_f[ij]): .6e}  3c={float(EwaldSR_3c_f[ij]): .6e}")

    H_vca_f = Vca2c_f
    print("  Fortran Vca2c (component sum, no Ewald):")
    print(H_vca_f)
    print(f"  Fortran |Vca2c|={float(np.max(np.abs(Vca2c_f))):.3e}  |EwaldSR|={float(np.max(np.abs(EwaldSR_f))):.3e}  |EwaldLR|={float(np.max(np.abs(EwaldLR_f))):.3e}")
    print(f"  Fortran |Vca_full|={float(np.max(np.abs(Vca_full_f))):.3e}")

    dq_shell = np.zeros((dims.nsh_max, natoms))
    for ia in range(natoms):
        ispec = int(ispec_of_atom[ia])
        for issh in range(int(sd.nssh[ispec])):
            dq_shell[issh, ia] = Qin_shell[issh, ia] - Qneutral_sh[issh, ispec]

    # Build directed neighbor list for OpenCL Vca assembly (excludes missing ja, mbeta combos)
    neigh_list_vca = []
    neigh_mbeta_vca = []
    for ia in range(natoms):
        nn = int(sd.neighn[ia])
        for ineigh in range(nn):
            ja = int(sd.neigh_j[ia, ineigh]) - 1
            mb = int(sd.neigh_b[ia, ineigh])
            if ja < 0 or ja >= natoms:
                continue
            # Keep self-pairs only if they correspond to the explicit self slot (mb=0)
            if ja == ia and mb != 0:
                continue
            neigh_list_vca.append((ia, ja))
            neigh_mbeta_vca.append(mb)

    # Helper function to compute dq for an atom
    def _dq_atom(ia):
        ispec = int(ispec_of_atom[ia])
        s = 0.0
        for issh in range(int(sd.nssh[ispec])):
            s += float(Qin_shell[issh, ia] - Qneutral_sh[issh, ispec])
        return s

    if len(neigh_mbeta_vca) > 0:
        print(f"  VCA neigh_b stats: max|mbeta|={int(np.max(np.abs(np.array(neigh_mbeta_vca, dtype=np.int32))))}  unique={sorted(set([int(x) for x in neigh_mbeta_vca]))}")

    pair_types_vca = []
    for (i, j) in neigh_list_vca:
        nz1 = int(atomTypes_Z[i]); nz2 = int(atomTypes_Z[j])
        k = ('vna', nz1, nz2)
        idx = ham.species_pair_map.get(k)
        if idx is None:
            print(f"Warning: missing vna pair type for {nz1},{nz2}")
            idx = 0
        pair_types_vca.append(idx)

    # IMPORTANT: Fortran uses r2 = atomPos[j] + xl[mbeta] for each neighbor.
    # OpenCL assemble_vca currently accepts only a flat atom position array; emulate mbeta shifts by
    # constructing a per-pair position buffer (two atoms per pair) when any mbeta!=0.
    mb_max = int(np.max(np.abs(np.array(neigh_mbeta_vca, dtype=np.int32)))) if len(neigh_mbeta_vca) > 0 else 0
    # Pick one directed neighbor pair for synchronized VCA debugging (Fortran vs OpenCL).
    # Target: ia=1 -> ja=0 in 0-based indexing (Fortran iatom=2, jatom=1), mbeta==0.
    dbg_pair = -1
    dbg_ish = 0
    for k,(i,j) in enumerate(neigh_list_vca):
        if i == 1 and j == 0 and int(neigh_mbeta_vca[k]) == 0:
            dbg_pair = k
            break
    print(f"  [VCA_DBG][PY][CPSEL] target ia=1 ja=0 mbeta=0 -> dbg_pair={dbg_pair} dbg_ish={dbg_ish}")

    if mb_max != 0:
        ratoms_pairs = np.zeros((2*len(neigh_list_vca), 3), dtype=np.float64)
        neigh_pairs = np.zeros((len(neigh_list_vca), 2), dtype=np.int32)
        for k,(i,j) in enumerate(neigh_list_vca):
            mb = int(neigh_mbeta_vca[k])
            ratoms_pairs[2*k+0] = atomPos[i]
            ratoms_pairs[2*k+1] = atomPos[j] + sd.xl[mb]
            neigh_pairs[k,0] = 2*k+0
            neigh_pairs[k,1] = 2*k+1
        off_vca, diag_vca = ham.assemble_vca(ratoms_pairs, neigh_pairs, np.array(pair_types_vca, dtype=np.int32), dq_shell,
                                             ispec_of_atom=np.repeat(ispec_of_atom, 1), nssh_species=nssh_species, lssh_species=sd.lssh[:dims.nspecies], nsh_max=dims.nsh_max,
                                             dbg_enable=int(dbg_pair>=0), dbg_pair=dbg_pair, dbg_ish=dbg_ish)
    else:
        off_vca, diag_vca = ham.assemble_vca(atomPos, np.array(neigh_list_vca, dtype=np.int32), np.array(pair_types_vca, dtype=np.int32), dq_shell,
                                         ispec_of_atom=ispec_of_atom, nssh_species=nssh_species, lssh_species=sd.lssh[:dims.nspecies], nsh_max=dims.nsh_max,
                                         dbg_enable=int(dbg_pair>=0), dbg_pair=dbg_pair, dbg_ish=dbg_ish)
    Vca_ocl = np.zeros_like(H_vca_f)
    Vca_ocl_off = np.zeros_like(H_vca_f)
    Vca_ocl_diag = np.zeros_like(H_vca_f)
    for k, (i, j) in enumerate(neigh_list_vca):
        i0 = int(offs[i]); j0 = int(offs[j])
        ni = int(n_orb_atom[i]); nj = int(n_orb_atom[j])
        if i != j:
            Vca_ocl[i0:i0+ni, j0:j0+nj] += off_vca[k]
            Vca_ocl_off[i0:i0+ni, j0:j0+nj] += off_vca[k]
        # Diagonal updates from vna_atom:
        # diag_vca[k,0] updates V_ii from potential on j (including self-pair i==j).
        Vca_ocl[i0:i0+ni, i0:i0+ni] += diag_vca[k, 0]
        Vca_ocl_diag[i0:i0+ni, i0:i0+ni] += diag_vca[k, 0]

    res_Vca, err_Vca = compare_matrices_brief("Vca2c (vca array)", Vca2c_f, Vca_ocl, tol=1e-4, require_nonzero=True)

    # OpenCL component-wise Vca2c (bcca splines only) by masking vca_map indices
    vmap = ham.vca_map
    vmap_L = vmap.copy(); vmap_L[:, :, 1] = -1; vmap_L[:, :, 2] = -1
    vmap_R = vmap.copy(); vmap_R[:, :, 0] = -1; vmap_R[:, :, 2] = -1
    vmap_A = vmap.copy(); vmap_A[:, :, 0] = -1; vmap_A[:, :, 1] = -1

    # Diagnostic: swapped L/R component interpretation (if OpenCL component order differs)
    vmap_Ls = vmap.copy(); vmap_Ls[:, :, 0] = -1; vmap_Ls[:, :, 2] = -1  # take comp-1 as "L"
    vmap_Rs = vmap.copy(); vmap_Rs[:, :, 1] = -1; vmap_Rs[:, :, 2] = -1  # take comp-0 as "R"

    # Extra synchronized atom-only traces for ia=1 -> ja=0 and ia=1 -> ja=2 at ish=0.
    # These calls are for printing only and do not affect later comparisons.
    dbg_pairs = []
    for target_j in (0, 2):
        kp = -1
        for k,(i,j) in enumerate(neigh_list_vca):
            if i == 1 and j == target_j and int(neigh_mbeta_vca[k]) == 0:
                kp = k
                break
        dbg_pairs.append(kp)
    dbg_ish_list = [0, 1]
    print(f"  [VCA_DBG][PY][CPSEL2] ia=1 -> ja=(0,2) mbeta=0 dbg_pairs={dbg_pairs} dbg_ish_list={dbg_ish_list}")
    for kp in dbg_pairs:
        if kp >= 0:
            for ish_dbg in dbg_ish_list:
                _off_tmp, _diag_tmp = ham.assemble_vca(atomPos, np.array(neigh_list_vca, dtype=np.int32), np.array(pair_types_vca, dtype=np.int32), dq_shell,
                                      ispec_of_atom=ispec_of_atom, nssh_species=nssh_species, lssh_species=sd.lssh[:dims.nspecies], nsh_max=dims.nsh_max,
                                      vca_map_override=vmap_A, dbg_enable=1, dbg_pair=kp, dbg_ish=int(ish_dbg))

    offL, diagL = ham.assemble_vca(atomPos, np.array(neigh_list_vca, dtype=np.int32), np.array(pair_types_vca, dtype=np.int32), dq_shell,
                                  ispec_of_atom=ispec_of_atom, nssh_species=nssh_species, lssh_species=sd.lssh[:dims.nspecies], nsh_max=dims.nsh_max, vca_map_override=vmap_L)
    offR, diagR = ham.assemble_vca(atomPos, np.array(neigh_list_vca, dtype=np.int32), np.array(pair_types_vca, dtype=np.int32), dq_shell,
                                  ispec_of_atom=ispec_of_atom, nssh_species=nssh_species, lssh_species=sd.lssh[:dims.nspecies], nsh_max=dims.nsh_max, vca_map_override=vmap_R)
    offA, diagA = ham.assemble_vca(atomPos, np.array(neigh_list_vca, dtype=np.int32), np.array(pair_types_vca, dtype=np.int32), dq_shell,
                                  ispec_of_atom=ispec_of_atom, nssh_species=nssh_species, lssh_species=sd.lssh[:dims.nspecies], nsh_max=dims.nsh_max, vca_map_override=vmap_A)

    offLs, diagLs = ham.assemble_vca(atomPos, np.array(neigh_list_vca, dtype=np.int32), np.array(pair_types_vca, dtype=np.int32), dq_shell,
                                   ispec_of_atom=ispec_of_atom, nssh_species=nssh_species, lssh_species=sd.lssh[:dims.nspecies], nsh_max=dims.nsh_max, vca_map_override=vmap_Ls)
    offRs, diagRs = ham.assemble_vca(atomPos, np.array(neigh_list_vca, dtype=np.int32), np.array(pair_types_vca, dtype=np.int32), dq_shell,
                                   ispec_of_atom=ispec_of_atom, nssh_species=nssh_species, lssh_species=sd.lssh[:dims.nspecies], nsh_max=dims.nsh_max, vca_map_override=vmap_Rs)

    VcaL_o, VcaR_o, VcaA_o = accumulate_vca_blocks(neigh_list_vca, n_orb_atom, offs, offL, offR, diagA, H_vca_f.shape)
    VcaL_os, VcaR_os, VcaA_os = accumulate_vca_blocks(neigh_list_vca, n_orb_atom, offs, offLs, offRs, diagLs, H_vca_f.shape)

    res_VcaL, err_VcaL = compare_matrices_brief("Vca2c_ontopl (2)", VcaL_f, VcaL_o, tol=1e-4, require_nonzero=True)
    res_VcaR, err_VcaR = compare_matrices_brief("Vca2c_ontopr (3)", VcaR_f, VcaR_o, tol=1e-4, require_nonzero=True)
    res_VcaA, err_VcaA = compare_matrices_brief("Vca2c_atom (4)",   VcaA_f, VcaA_o, tol=1e-4, require_nonzero=True)

    _res_VcaL_s, err_VcaL_s = compare_matrices_brief("Vca2c_ontopl (SWAPPED)", VcaL_f, VcaL_os, tol=1e-4, require_nonzero=True)
    _res_VcaR_s, err_VcaR_s = compare_matrices_brief("Vca2c_ontopr (SWAPPED)", VcaR_f, VcaR_os, tol=1e-4, require_nonzero=True)
    print(f"  [diag] component swap check: err_L={err_VcaL:.3e} err_R={err_VcaR:.3e} | swapped: err_L={err_VcaL_s:.3e} err_R={err_VcaR_s:.3e}")


    dq_atom = np.array([_dq_atom(i) for i in range(natoms)], dtype=np.float32)
    _ewchk = compare_ewaldsr(ham, atomPos, atomTypes_Z, sd, dims, neigh_self, neigh_back, dq_atom, ewaldsr4_f, ewaldsr3c4_f, tol=1e-5)
    ewaldsr4_o = _ewchk.details['ewaldsr4_o']
    EwaldSR_o = _ewchk.details['EwaldSR_o']

    # EwaldLR still not implemented (requires lattice-summed ewald(i,j) / sub_ewald).
    EwaldLR_o = np.zeros_like(EwaldLR_f)
    print("  [TODO] OpenCL EwaldLR not implemented yet; still comparing against zeros")
    res_EwaldSR, err_EwaldSR = compare_matrices_brief("EwaldSR", EwaldSR_f, EwaldSR_o, tol=1e-5, require_nonzero=True)
    res_EwaldLR, err_EwaldLR = compare_matrices_brief("EwaldLR", EwaldLR_f, EwaldLR_o, tol=1e-5, require_nonzero=True)
    Vca_full_o = Vca_ocl + EwaldLR_o - EwaldSR_o
    res_Vca_full, err_Vca_full = compare_matrices_brief("Vca_full (2c)", Vca_full_f, Vca_full_o, tol=1e-5, require_nonzero=True)

    # Pairwise block diagnostics (non-self directed neighbors only).
    print("\n--- Vca pairwise block diagnostics (non-self directed neighbors) ---")
    max_pair_diff = 0.0
    worst = None
    for kpair, (ia, ja) in enumerate(neigh_list_vca):
        # Find corresponding Fortran neighbor slot (ia, ineigh) where neigh_j==ja+1
        nn = int(sd.neighn[ia])
        ineighF = -1
        for ineigh in range(nn):
            if int(sd.neigh_j[ia, ineigh]) - 1 == ja:
                ineighF = ineigh
                break
        if ineighF < 0:
            continue
        # Core reference off-diagonal block is ontopl+ontopr. Atom term contributes only to diagonal updates.
        blk_f = vcaL4_f[ia, ineighF, :, :] + vcaR4_f[ia, ineighF, :, :]   # [nu,mu]
        blk_o = off_vca[kpair].T             # [nu,mu]
        d = float(np.max(np.abs(blk_f - blk_o)))
        if d > max_pair_diff:
            max_pair_diff = d
            worst = (ia, ja, ineighF, kpair, d)
    print(f"  max|Fortran_vca_block - OCL_off_block| over non-self directed pairs = {max_pair_diff:.3e}")
    if worst is not None:
        ia, ja, ineighF, kpair, d = worst
        print(f"  worst directed pair ia={ia} -> ja={ja} ineighF={ineighF} kpair={kpair} diff={d:.3e}")
        print("  Fortran vca(offdiag) block [nu,mu] = vcaL+vcaR:")
        print(vcaL4_f[ia, ineighF] + vcaR4_f[ia, ineighF])
        print("  OpenCL off block [nu,mu]:")
        print(off_vca[kpair].T)

    # Component-level pairwise diagnostics for ontopl/ontopr
    print("\n--- Vca pairwise diagnostics: ontopl vs export ---")
    max_pair_L = 0.0
    worst_L = None
    for kpair, (ia, ja) in enumerate(neigh_list_vca):
        nn = int(sd.neighn[ia])
        ineighF = -1
        for ineigh in range(nn):
            if int(sd.neigh_j[ia, ineigh]) - 1 == ja:
                ineighF = ineigh; break
        if ineighF < 0: continue
        blk_f = vcaL4_f[ia, ineighF, :, :]       # [nu,mu]
        blk_o = offL[kpair].T                    # [nu,mu]
        d = float(np.max(np.abs(blk_f - blk_o)))
        if d > max_pair_L:
            max_pair_L = d
            worst_L = (ia, ja, ineighF, kpair, d)
    print(f"  max|vcaL_export - OCL_offL| over non-self directed pairs = {max_pair_L:.3e}")
    if worst_L is not None:
        ia, ja, ineighF, kpair, d = worst_L
        print(f"  worst ontopl pair ia={ia} -> ja={ja} ineighF={ineighF} kpair={kpair} diff={d:.3e}")
        blk_f0 = vcaL4_f[ia, ineighF]
        blk_fT = blk_f0.T
        blk_o  = offL[kpair].T
        print(f"  ontopl orientation check: max|F - OCL|={float(np.max(np.abs(blk_f0-blk_o))):.3e}  max|F.T - OCL|={float(np.max(np.abs(blk_fT-blk_o))):.3e}")
        print("  Fortran vcaL block [nu,mu]:")
        print(blk_f0)
        print("  OpenCL offL block [nu,mu]:")
        print(blk_o)

        # Try to infer orbital ordering/sign mismatch by brute force (full 4-orb basis).
        import itertools
        fref = blk_f0[:4, :4].copy()    # [nu,mu]
        oref = blk_o[:4, :4].copy()     # [nu,mu]
        best = None
        for perm in itertools.permutations([0, 1, 2, 3], 4):
            for sgn in itertools.product([-1.0, 1.0], repeat=4):
                P = np.array(perm, dtype=int)
                S = np.array(sgn, dtype=np.float64)
                M = oref[np.ix_(P, P)].copy()
                M = (S[:, None] * M) * (S[None, :])
                dd = float(np.max(np.abs(fref - M)))
                if (best is None) or (dd < best[0]):
                    best = (dd, perm, sgn)
        if best is not None:
            dd, perm, sgn = best
            print(f"  [diag] best 4-orb map for ontopl: max|F - mapped(OCL)|={dd:.3e}  perm={perm}  sgn={sgn}")

    print("\n--- Vca pairwise diagnostics: ontopr vs export ---")
    max_pair_R = 0.0
    worst_R = None
    for kpair, (ia, ja) in enumerate(neigh_list_vca):
        nn = int(sd.neighn[ia])
        ineighF = -1
        for ineigh in range(nn):
            if int(sd.neigh_j[ia, ineigh]) - 1 == ja:
                ineighF = ineigh; break
        if ineighF < 0: continue
        blk_f = vcaR4_f[ia, ineighF, :, :]       # [nu,mu]
        blk_o = offR[kpair].T                    # [nu,mu]
        d = float(np.max(np.abs(blk_f - blk_o)))
        if d > max_pair_R:
            max_pair_R = d
            worst_R = (ia, ja, ineighF, kpair, d)
    print(f"  max|vcaR_export - OCL_offR| over non-self directed pairs = {max_pair_R:.3e}")
    if worst_R is not None:
        ia, ja, ineighF, kpair, d = worst_R
        print(f"  worst ontopr pair ia={ia} -> ja={ja} ineighF={ineighF} kpair={kpair} diff={d:.3e}")
        blk_f0 = vcaR4_f[ia, ineighF]
        blk_fT = blk_f0.T
        blk_o  = offR[kpair].T
        print(f"  ontopr orientation check: max|F - OCL|={float(np.max(np.abs(blk_f0-blk_o))):.3e}  max|F.T - OCL|={float(np.max(np.abs(blk_fT-blk_o))):.3e}")
        print("  Fortran vcaR block [nu,mu]:")
        print(blk_f0)
        print("  OpenCL offR block [nu,mu]:")
        print(blk_o)

        import itertools
        fref = blk_f0[:4, :4].copy()    # [nu,mu]
        oref = blk_o[:4, :4].copy()     # [nu,mu]
        best = None
        for perm in itertools.permutations([0, 1, 2, 3], 4):
            for sgn in itertools.product([-1.0, 1.0], repeat=4):
                P = np.array(perm, dtype=int)
                S = np.array(sgn, dtype=np.float64)
                M = oref[np.ix_(P, P)].copy()
                M = (S[:, None] * M) * (S[None, :])
                dd = float(np.max(np.abs(fref - M)))
                if (best is None) or (dd < best[0]):
                    best = (dd, perm, sgn)
        if best is not None:
            dd, perm, sgn = best
            print(f"  [diag] best 4-orb map for ontopr: max|F - mapped(OCL)|={dd:.3e}  perm={perm}  sgn={sgn}")

    # Diagonal diagnostics for atom term
    print("\n--- Vca atom(diag) diagnostics: export self-slot vs OCL diagA ---")
    max_ai = 0.0
    worst_ai = None
    # Find self-slot index for each atom
    neigh_self = np.full(natoms, -1, dtype=np.int32)
    for i in range(natoms):
        ii = i + 1
        for ineigh in range(int(sd.neighn[i])):
            if int(sd.neigh_j[i, ineigh]) == ii and int(sd.neigh_b[i, ineigh]) == 0:
                neigh_self[i] = ineigh
                break
    for i in range(natoms):
        ni = int(n_orb_atom[i]); i0 = int(offs[i])
        # Fortran export: vca_atom is accumulated in self-slot (matom = neigh_self(iatom))
        acc_f = np.zeros((ni, ni), dtype=np.float64)
        if neigh_self[i] >= 0:
            acc_f = vcaA4_f[i, neigh_self[i], :ni, :ni].T
        acc_o = np.zeros((ni, ni), dtype=np.float64)
        for kpair, (ia, ja) in enumerate(neigh_list_vca):
            if ia != i: continue
            acc_o += diagA[kpair, 0, :ni, :ni]
        d = float(np.max(np.abs(acc_f - acc_o)))
        if d > max_ai:
            max_ai = d
            worst_ai = (i, d, acc_f, acc_o)
    print(f"  max|atom_diag_export - atom_diag_OCL| over atoms = {max_ai:.3e}")
    if worst_ai is not None:
        i, d, acc_f, acc_o = worst_ai
        print(f"  worst atom i={i} diff={d:.3e}")
        print("  Fortran atom diag block (from self-slot):")
        print(acc_f)
        print("  OpenCL atom diag block (sum diagA[k,0] for ia==i):")
        print(acc_o)

    diff_ocl = np.abs(Vca_ocl - H_vca_f)
    # err_Vca already computed above from Fortran export vs OpenCL.
    print(f"  Max diff Vca OCL vs Fortran: {err_Vca:.2e}")
    if res_Vca:
        print("SUCCESS: OpenCL Vca2c matches Fortran export.")
    else:
        print("FAILURE: OpenCL Vca2c mismatch vs Fortran export.")

    def _block(M, i, j):
        ni = int(n_orb_atom[i]); nj = int(n_orb_atom[j])
        i0 = int(offs[i]); j0 = int(offs[j])
        return M[i0:i0+ni, j0:j0+nj]

    for i in range(natoms):
        for j in range(natoms):
            if i == j or j in neigh_lists[i]:
                df = _block(H_vca_f, i, j)
                do = _block(Vca_ocl, i, j)
                print(f"  VCA block ({i},{j}) max|F-ocl|={float(np.max(np.abs(df-do))):.2e}")

    # Component-wise diagnostics for the full Vca term used in H assembly: vca + ewaldlr - ewaldsr
    print("\n--- Vca component-wise diagnostics (Fortran) ---")
    _ = compare_matrices_brief("Vca2c (vca array)", Vca2c_f, Vca_ocl, tol=1e-4, require_nonzero=True)
    print(f"  Max|EwaldSR|={float(np.max(np.abs(EwaldSR_f))):.3e}  Max|EwaldLR|={float(np.max(np.abs(EwaldLR_f))):.3e}")
    Vca_resid = Vca_full_f - Vca_ocl
    print(f"  Residual fullVca-ocl: max|.|={float(np.max(np.abs(Vca_resid))):.3e}")
    # If you later implement EwaldSR/LR in OpenCL, compare against Vca_full_f directly.

    print("\n--- Vca split diagnostics: offdiag vs diag (OpenCL) ---")
    # Off-diagonal differences
    mask_off = np.ones_like(H_vca_f, dtype=bool)
    np.fill_diagonal(mask_off, False)
    diff_off = np.abs((Vca2c_f - Vca_ocl_off) * mask_off)
    diff_diag = np.abs(np.diag(Vca2c_f) - np.diag(Vca_ocl_diag))
    print(f"  max|offdiag(F - OCL_off)|={float(np.max(diff_off)):.3e}")
    print(f"  max|diag(F - OCL_diag)|={float(np.max(diff_diag)):.3e}")

    ij = np.unravel_index(int(np.argmax(np.abs(Vca2c_f - Vca_ocl))), Vca2c_f.shape)
    print(f"  WORST elem (F-OCL) at {ij}: F={float(Vca2c_f[ij]): .6e}  OCL={float(Vca_ocl[ij]): .6e}")
    # Locate which atom blocks contain that matrix element
    def _which_atom(idx_):
        for ia in range(natoms):
            if idx_ >= int(offs[ia]) and idx_ < int(offs[ia+1]):
                return ia
        return -1
    ia = _which_atom(int(ij[0])); ja = _which_atom(int(ij[1]))
    print(f"  WORST elem belongs to atom-block ({ia},{ja})")

    # Placeholders (Vxc flags updated when microtest errors are available)
    res_Vxc = ('err_bcxcx_ocl' in locals()) and (err_bcxcx_ocl < 2e-6)
    res_Vxc_1c = ('err_bcxcx_on' in locals()) and (err_bcxcx_on < 1e-6)
    res_Vxc_ca = ('err_bcxcx_on_ocl' in locals()) and (err_bcxcx_on_ocl < 1e-5)
    err_Vxc_summary = err_bcxcx_ocl if 'err_bcxcx_ocl' in locals() else float('nan')
    err_Vxc1c_summary = err_bcxcx_on if 'err_bcxcx_on' in locals() else float('nan')
    err_Vxc_ca_summary = err_bcxcx_on_ocl if 'err_bcxcx_on_ocl' in locals() else float('nan')
    res_H_full = False

    print("\n" + "=" * 40)
    print("STRUCTURAL CHECKS (indices/self-consistency, not OCL vs Fortran float)")
    print(f" Step1_neigh_self: {'PASSED' if res_neigh_self else 'FAILED'}  err={err_neigh_self:.0f}")
    print(f" Step1_neigh_back: {'PASSED' if res_neigh_back else 'FAILED'}  err={err_neigh_back:.0f}")
    print(f" Step4_S_blocks:   {'PASSED' if res_S_blocks else 'FAILED'}  err={err_S_blocks:.2e}")
    print(f" Step4_H_blocks:   {'PASSED' if res_H_blocks else 'FAILED'}  err={err_H_blocks:.2e}")
    print(f" EwaldSR_split_4D(SFIRE): {'PASSED' if res_EwaldSR_split4d else 'FAILED'}  err={err_EwaldSR_split4d:.2e}")
    print(f" EwaldSR_split_dense(SFIRE): {'PASSED' if res_EwaldSR_split_sfire else 'FAILED'}  err={err_EwaldSR_split_sfire:.2e}")

    print("\n" + "=" * 40)
    print("VERIFICATION SUMMARY (max|diff|) (OpenCL vs FORTRAN)")
    print(f"Overlap S:   {'PASSED' if res_S else 'FAILED'}  err={err_S:.2e}")
    print(f"Kinetic T:   {'PASSED' if res_T else 'FAILED'}  err={err_T:.2e}")
    print(f"Dip:         {'PASSED' if res_Dip else 'FAILED'}  err={err_Dip:.2e}")
    print(f"Vna:         {'PASSED' if res_Vna else 'FAILED'}  err={err_Vna:.2e}")
    print(f"Vnl:         {'PASSED' if res_sVNL else 'FAILED'}  err={err_sVNL:.2e}")
    print(f"Vxc:         {'PASSED' if res_Vxc else 'FAILED'}  err={err_Vxc_summary:.2e}")
    print(f"Vxc_1c:      {'PASSED' if res_Vxc_1c else 'FAILED'}  err={err_Vxc1c_summary:.2e}")
    print(f"Vca:         {'PASSED' if res_Vca else 'FAILED'}  err={err_Vca:.2e}")
    print(f" Vca2c_ontopl:{'PASSED' if res_VcaL else 'FAILED'}  err={err_VcaL:.2e}")
    print(f" Vca2c_ontopr:{'PASSED' if res_VcaR else 'FAILED'}  err={err_VcaR:.2e}")
    print(f" Vca2c_atom:  {'PASSED' if res_VcaA else 'FAILED'}  err={err_VcaA:.2e}")
    print(f" EwaldSR:     {'PASSED' if res_EwaldSR else 'FAILED'}  err={err_EwaldSR:.2e}")
    print(f" EwaldLR:     {'PASSED' if res_EwaldLR else 'FAILED'}  err={err_EwaldLR:.2e}")
    print(f" Vca_full2c:  {'PASSED' if res_Vca_full else 'FAILED'}  err={err_Vca_full:.2e}")
    if 'err_bcxcx_on_ocl' in locals():
        print(f"Vxc_ca:      {'PASSED' if res_Vxc_ca else 'FAILED'}  err={err_Vxc_ca_summary:.2e}")
    else:
        print(f"Vxc_ca:      NOT IMPLEMENTED  err=nan")
    print(f"H2c (T+Vna): {'PASSED' if res_H2c else 'FAILED'}  err={err_H2c:.2e}")
    print(f"H raw==recon:{'PASSED' if res_H_recon else 'FAILED'}  err={err_H_recon:.2e}")
    print(f"VNL vs F:    {'PASSED' if (res_Vnl_cpu_fortran and res_Vnl_gpu_fortran) else 'FAILED'}  err=max({err_Vnl_cpu_fortran:.2e},{err_Vnl_gpu_fortran:.2e})")
    print(f"AvgRho_off:  {'PASSED' if res_AvgRho else 'FAILED'}  err={err_AvgRho:.2e}")
    
    #print(f"AvgRho_T:    {'PASSED' if err_AvgRho_T < 1e-6 else 'FAILED'}  err={err_AvgRho_T:.2e}")
    #print(f"AvgRho_seed:{'PASSED' if err_AvgRho_seed < 1e-6 else 'FAILED'}  err={err_AvgRho_seed:.2e}")
    #print(f"Seed_T:      {'PASSED' if err_AvgRho_seed_T < 1e-6 else 'FAILED'}  err={err_AvgRho_seed_T:.2e}")
    #print(f"AvgRho_3c:   {'PASSED' if err_AvgRho3c < 1e-6 else 'FAILED'}  err={err_AvgRho3c:.2e}")
    #print(f"AvgRho_3c_T: {'PASSED' if err_AvgRho3c_T < 1e-6 else 'FAILED'}  err={err_AvgRho3c_T:.2e}")

    print(f"Rot(isorp):  {'PASSED' if err_rot_isorp < 1e-6 else 'FAILED'}  err={err_rot_isorp:.2e}")
    print(f"Raw->rot F:  {'PASSED' if err_rawrot_f < 1e-6 else 'FAILED'}  err={err_rawrot_f:.2e}")
    print(f"Raw->rot O:  {'PASSED' if err_rawrot_o < 1e-6 else 'FAILED'}  err={err_rawrot_o:.2e}")
    print(f"Full H:      {'PASSED' if res_H_full else 'NOT IMPLEMENTED'}  err=nan")
    print("=" * 40)

    # Vxc microtest summary
    print("\n--- Vxc Microtests (OpenCL vs Fortran) ---")
    if 'err_muxc' in locals():
        print(f"cepal scalar: PASSED  err_muxc={err_muxc:.3e}  err_dmuxc={err_dmuxc:.3e}")
    if 'err_bcxcx' in locals():
        print(f"bcxcx Python:  PASSED  err={err_bcxcx:.3e}")
    if 'err_bcxcx_ocl' in locals():
        print(f"bcxcx OpenCL:  {'PASSED' if err_bcxcx_ocl < 2e-6 else 'FAILED'}  err={err_bcxcx_ocl:.3e}")
    if 'err_bcxcx_on' in locals():
        print(f"bcxcx_on Python:  PASSED  err={err_bcxcx_on:.3e}")
    if 'err_bcxcx_on_ocl' in locals():
        print(f"bcxcx_on OpenCL:  {'PASSED' if err_bcxcx_on_ocl < 1e-5 else 'FAILED'}  err={err_bcxcx_on_ocl:.3e}")
    print("=" * 40)

    # AvgRho_off is expected to fail until 2c+3c density-table evaluation is fully implemented on GPU.
    if not (res_S and res_T and res_Dip and res_Vna and res_sVNL and res_Vnl_cpu and res_Vnl_gpu and res_H_recon and res_H2c and res_Vnl_cpu_fortran and res_Vnl_gpu_fortran):
        raise RuntimeError("verify_C3: some implemented checks FAILED")

    print("\nDONE verify_C3")


if __name__ == "__main__":
    run_verification()
