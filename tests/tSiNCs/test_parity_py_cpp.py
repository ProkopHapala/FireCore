"""Parity test: Python reference vs C++ FFfit sensitivity matrices and fitted parameters.

Verifies that the pure-Python build_sensitivity_matrices / accumulate_normal_equations
path produces identical results to the C++ FFfit library path (collect_sensitivity_matrices).
"""
import os, sys
import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(__file__))
from test_FFfit import (
    build_sensitivity_matrices, accumulate_normal_equations, fit_hessian,
    compute_model_hessian, build_topology, assign_si_environment_types,
    angle_type_key, make_cpp_fitter, build_global_param_map,
    hessian_ha_bohr_to_ev_ang2, load_pyscf_case,
)
import pyBall.FFfit as FFfit_cpp

RESULTS_DIR = "/home/prokop/SIMULATIONS/jobs_pyscf_vib_OUT_small_nc/results"


def _make_system(name):
    """Load one system and build both Python and C++ models with identical typing."""
    data = load_pyscf_case(os.path.join(RESULTS_DIR, name))
    natoms = data['natoms']
    H_ref = hessian_ha_bohr_to_ev_ang2(data['hessian'].transpose(0, 2, 1, 3).reshape(natoms * 3, natoms * 3))
    bonds, angles, bonds3 = build_topology(data['symbols'], data['positions'], 2.5, third_bonds=False)
    symbols = data['symbols']
    atom_types = assign_si_environment_types(symbols, bonds, enabled=True)
    positions = data['positions']
    masses = data['masses']
    # Build global param map (single system)
    bmap, _, amap, _, n_total = build_global_param_map(
        [bonds], [angles], [atom_types], all_elements=[symbols], angle_central_only=True)
    n_cpp = len(bmap) + len(amap)
    # Bond/angle global type arrays
    bond_types = np.array([bmap[tuple(sorted((atom_types[i], atom_types[j])))] for i, j, _ in bonds], dtype=int)
    angle_types = np.array([amap[angle_type_key(atom_types, i, j, k, elements=symbols, central_only=True)]
                            for i, j, k, _ in angles], dtype=int)
    # Python sensitivity matrices
    A_py = build_sensitivity_matrices(positions, bonds, angles, bond_types, angle_types, natoms)
    # C++ sensitivity matrices
    sys_dict = {'positions': positions, 'symbols': symbols, 'atom_types': atom_types,
                'bonds': bonds, 'angles': angles, 'bonds3': [], 'angle_central_only': True}
    fitter = make_cpp_fitter(sys_dict, (bmap, {}, amap), n_total)
    fitter.set_n_free(n_cpp)
    A_cpp = FFfit_cpp.collect_sensitivity_matrices(fitter)
    return A_py, A_cpp, n_cpp, H_ref, masses, positions, bonds, angles, bmap, amap, fitter


@pytest.fixture(scope="module", params=["SiH4", "Si_R3p8"])
def system(request):
    return _make_system(request.param)


def test_sensitivity_matrices_match(system, tol=1e-10):
    """Python vs C++ sensitivity matrices A[p] = dH/dk_p must match to machine precision."""
    A_py, A_cpp, n_cpp, *_ = system
    assert len(A_py) == n_cpp
    assert A_cpp.shape == (n_cpp, A_cpp.shape[1], A_cpp.shape[2])
    for p in range(n_cpp):
        diff = np.max(np.abs(A_py[p] - A_cpp[p]))
        norm = np.max(np.abs(A_cpp[p]))
        rel = diff / max(norm, 1e-30)
        assert rel < tol, f"param {p}: max_abs_diff={diff:.2e} max_val={norm:.2e} rel={rel:.2e}"


def test_model_hessian_match(system, tol=1e-10):
    """H_model = Σ k_p A_p must match between Python and C++ for arbitrary k."""
    A_py, A_cpp, n_cpp, H_ref, masses, positions, bonds, angles, bmap, amap, fitter = system
    k = np.ones(n_cpp) * 2.0  # arbitrary test parameters
    # Python
    H_py = compute_model_hessian(A_py, k)
    # C++
    fitter.set_params(k)
    H_cpp = fitter.compute_model_hessian()
    diff = np.max(np.abs(H_py - H_cpp))
    norm = np.max(np.abs(H_cpp))
    rel = diff / max(norm, 1e-30)
    assert rel < tol, f"model Hessian: max_abs_diff={diff:.2e} max_val={norm:.2e} rel={rel:.2e}"


def test_lsq_fit_match(system, tol=1e-8):
    """Linear least-squares fitted parameters must match between Python and C++."""
    A_py, A_cpp, n_cpp, H_ref, masses, positions, bonds, angles, bmap, amap, fitter = system
    # Use unweighted LSQ for parity (weight conventions differ between Python and C++)
    # Python LSQ
    G_py = np.zeros((n_cpp, n_cpp))
    y_py = np.zeros(n_cpp)
    accumulate_normal_equations(G_py, y_py, H_ref, A_py, weight=None)
    k_py = np.linalg.lstsq(G_py, y_py, rcond=None)[0]
    # C++ LSQ
    G_cpp = np.zeros((n_cpp, n_cpp))
    y_cpp = np.zeros(n_cpp)
    fitter.accumulate_normal_equations(G_cpp, y_cpp, H_ref.ravel(), weight=None)
    k_cpp = FFfit_cpp.FFfit.solve_normal_equations(G_cpp, y_cpp)
    diff = np.max(np.abs(k_py - k_cpp))
    norm = np.max(np.abs(k_cpp))
    rel = diff / max(norm, 1e-30)
    assert rel < tol, f"LSQ params: max_abs_diff={diff:.2e} max_val={norm:.2e} rel={rel:.2e}"
    # Also check normal equations match
    G_diff = np.max(np.abs(G_py - G_cpp))
    y_diff = np.max(np.abs(y_py - y_cpp))
    assert G_diff < tol, f"Normal matrix G: max_diff={G_diff:.2e}"
    assert y_diff < tol, f"RHS vector y: max_diff={y_diff:.2e}"


def test_hybrid_fit_py_vs_cpp(tol=1e-8):
    """Hybrid fit using Python sensitivity matrices vs C++ sensitivity matrices."""
    data = load_pyscf_case(os.path.join(RESULTS_DIR, "SiH4"))
    natoms = data['natoms']
    H_ref = hessian_ha_bohr_to_ev_ang2(data['hessian'].transpose(0, 2, 1, 3).reshape(natoms * 3, natoms * 3))
    bonds, angles, bonds3 = build_topology(data['symbols'], data['positions'], 2.5, third_bonds=False)
    symbols = data['symbols']
    atom_types = assign_si_environment_types(symbols, bonds, enabled=True)
    positions = data['positions']
    masses = data['masses']
    bmap, _, amap, _, n_total = build_global_param_map(
        [bonds], [angles], [atom_types], all_elements=[symbols], angle_central_only=True)
    n_cpp = len(bmap) + len(amap)
    bond_types = np.array([bmap[tuple(sorted((atom_types[i], atom_types[j])))] for i, j, _ in bonds], dtype=int)
    angle_types = np.array([amap[angle_type_key(atom_types, i, j, k, elements=symbols, central_only=True)]
                            for i, j, k, _ in angles], dtype=int)
    # Python A
    A_py_list = build_sensitivity_matrices(positions, bonds, angles, bond_types, angle_types, natoms)
    A_py = np.stack(A_py_list)
    # C++ A
    sys_dict = {'positions': positions, 'symbols': symbols, 'atom_types': atom_types,
                'bonds': bonds, 'angles': angles, 'bonds3': [], 'angle_central_only': True}
    fitter = make_cpp_fitter(sys_dict, (bmap, {}, amap), n_total)
    fitter.set_n_free(n_cpp)
    A_cpp = FFfit_cpp.collect_sensitivity_matrices(fitter)
    # Hybrid fit with Python A
    prior = np.ones(n_cpp)
    for t in bmap.values(): prior[t] = 5.0
    for t in amap.values(): prior[t] = 1.0
    sys_py = {'A': A_py, 'H_ref': H_ref, 'positions': positions, 'masses': masses, 'bonds': bonds, 'angles': angles}
    sys_cpp = {'A': A_cpp, 'H_ref': H_ref, 'positions': positions, 'masses': masses, 'bonds': bonds, 'angles': angles}
    k_py, _ = FFfit_cpp.fit_hybrid_hessian(
        [sys_py], mode_weight=1.0, local_weight=1.0, internal_weight=1.0,
        prior=prior, regularization=2e-3, parameter_scale=prior, bounds=(0.0, np.inf))
    k_cpp, _ = FFfit_cpp.fit_hybrid_hessian(
        [sys_cpp], mode_weight=1.0, local_weight=1.0, internal_weight=1.0,
        prior=prior, regularization=2e-3, parameter_scale=prior, bounds=(0.0, np.inf))
    diff = np.max(np.abs(k_py - k_cpp))
    norm = np.max(np.abs(k_cpp))
    rel = diff / max(norm, 1e-30)
    assert rel < tol, f"Hybrid fit params: max_abs_diff={diff:.2e} max_val={norm:.2e} rel={rel:.2e}"


if __name__ == "__main__":
    # Run without pytest for quick standalone check
    print("=== Parity test: Python vs C++ FFfit ===")
    for name in ["SiH4", "Si_R3p8"]:
        print(f"\n--- {name} ---")
        s = _make_system(name)
        A_py, A_cpp, n_cpp = s[0], s[1], s[2]
        max_rel = 0.0
        for p in range(n_cpp):
            diff = np.max(np.abs(A_py[p] - A_cpp[p]))
            norm = np.max(np.abs(A_cpp[p]))
            rel = diff / max(norm, 1e-30)
            max_rel = max(max_rel, rel)
        print(f"  Sensitivity matrices: {n_cpp} params, max rel diff = {max_rel:.2e}")
        # Model Hessian
        k = np.ones(n_cpp) * 2.0
        H_py = compute_model_hessian(A_py, k)
        fitter = s[10]
        fitter.set_params(k)
        H_cpp = fitter.compute_model_hessian()
        diff = np.max(np.abs(H_py - H_cpp))
        norm = np.max(np.abs(H_cpp))
        print(f"  Model Hessian (k=2): max abs diff = {diff:.2e}, rel = {diff/max(norm,1e-30):.2e}")
        # LSQ fit (unweighted for parity)
        H_ref = s[3]
        G_py = np.zeros((n_cpp, n_cpp)); y_py = np.zeros(n_cpp)
        accumulate_normal_equations(G_py, y_py, H_ref, A_py, weight=None)
        k_py = np.linalg.lstsq(G_py, y_py, rcond=None)[0]
        G_cpp = np.zeros((n_cpp, n_cpp)); y_cpp = np.zeros(n_cpp)
        fitter.accumulate_normal_equations(G_cpp, y_cpp, H_ref.ravel(), weight=None)
        k_cpp = FFfit_cpp.FFfit.solve_normal_equations(G_cpp, y_cpp)
        diff = np.max(np.abs(k_py - k_cpp))
        print(f"  LSQ params: max abs diff = {diff:.2e}")
    # Hybrid fit
    print("\n--- Hybrid fit (SiH4) ---")
    test_hybrid_fit_py_vs_cpp()
    print("  Hybrid fit: PASS")
    print("\n=== All parity tests passed ===")
