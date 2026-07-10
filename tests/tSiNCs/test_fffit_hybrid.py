import numpy as np
import pytest

from pyBall import FFfit
from pyBall.FFfit_utils import assign_si_environment_types, dihedral_angle, dihedral_energy_gradient, subtype_shrinkage_rows, family_mean_prior_rows, build_cross_param_maps, compute_cross_sensitivity
from pyBall.FFfit_plots import cluster_1d


def make_triatomic_problem():
    positions = np.array([[0.0, 0.0, 0.0], [1.25, 0.0, 0.0], [1.65, 1.05, 0.2]])
    masses = np.array([28.085, 28.085, 1.008])
    bonds = [(0, 1, 1.25), (1, 2, np.linalg.norm(positions[2] - positions[1]))]
    angles = [(0, 1, 2, 1.9)]
    B, _ = FFfit.build_wilson_matrix(positions, bonds, angles)
    A = np.stack((B[:2].T @ B[:2], B[2:].T @ B[2:]))
    return positions, masses, bonds, angles, B, A


def test_rigid_and_internal_coordinate_spaces_are_physical():
    positions, masses, bonds, angles, B, _ = make_triatomic_problem()
    rigid, vibration = FFfit.rigid_and_vibrational_bases(positions, masses)
    U, s = FFfit.internal_coordinate_basis(B, masses)
    F, info = FFfit.internal_hessian_projection(B.T @ B, B, masses)
    cartesian_rigid = rigid / np.sqrt(np.repeat(masses, 3))[:, None]
    assert rigid.shape == (9, 6)
    assert vibration.shape == (9, 3)
    assert U.shape == (9, 3)
    assert np.all(s > 0.0)
    assert np.linalg.norm(B @ cartesian_rigid) < 1e-12
    assert np.linalg.norm(U.T @ rigid) < 1e-12
    assert F.shape == (3, 3)
    assert info['rank'] == 3
    assert np.linalg.norm(F - F.T) < 1e-12


def test_hybrid_fit_exactly_recovers_representable_hessian():
    positions, masses, bonds, angles, _, A = make_triatomic_problem()
    k_ref = np.array([4.2, 0.85])
    H_ref = np.einsum('p,pij->ij', k_ref, A)
    system = {'A': A, 'H_ref': H_ref, 'positions': positions, 'masses': masses, 'bonds': bonds, 'angles': angles}
    k, diag = FFfit.fit_hybrid_hessian([system], mode_weight=1.0, local_weight=1.0, internal_weight=1.0,
                                       regularization=0.0, bounds=(0.0, np.inf))
    assert np.allclose(k, k_ref, rtol=1e-10, atol=1e-10)
    assert diag['relative_residual'] < 1e-12
    assert diag['rank'] == 2
    assert diag['systems'][0]['n_modes'] == 3
    assert diag['systems'][0]['internal_rank'] == 3


def test_wilson_projection_uses_dimensionless_bond_scaling_and_recovers_spring_constant():
    positions, masses, bonds, angles, B, _ = make_triatomic_problem()
    labels = [('bond', 0, 1), ('bond', 1, 2), ('angle', 0, 1, 2)]
    scale = FFfit.dimensionless_wilson_scale(positions, labels)
    k_bond = 7.3
    H = k_bond * np.outer(B[0], B[0])
    F, _ = FFfit.internal_hessian_projection(H, B[:1], masses, coordinate_scale=scale[:1])
    r0 = np.linalg.norm(positions[1] - positions[0])
    assert F[0, 0] / r0**2 == pytest.approx(k_bond, rel=1e-12, abs=1e-12)


def test_internal_hybrid_component_exactly_recovers_representable_curvature():
    positions, masses, bonds, angles, _, A = make_triatomic_problem()
    H_ref = 3.1*A[0] + 0.6*A[1]
    X, y, info = FFfit.assemble_hybrid_hessian_system(A, H_ref, positions, masses, bonds, angles,
                                                       mode_weight=0.0, local_weight=0.0, internal_weight=1.0)
    assert info['component_rows']['internal'] == info['internal_rank'] * (info['internal_rank'] + 1) // 2
    assert np.linalg.norm(X @ np.array([3.1, 0.6]) - y) < 1e-12


def test_coupled_regularization_shrinks_subtype_to_parent():
    X = np.array([[1.0, 0.0]])
    y = np.array([8.0])
    R = np.array([[20.0, -20.0]])
    k, _ = FFfit.solve_regularized_lsq(X, y, coupling_matrix=R, coupling_target=np.zeros(1), bounds=(-np.inf, np.inf))
    assert k[0] == pytest.approx(8.0, abs=1e-12)
    assert k[1] == pytest.approx(8.0, abs=1e-12)


def test_subtype_shrinkage_rows_leave_the_family_mean_free():
    names = ['bond:Si-Si', 'bond:Si-SiH', 'bond:Si-SiH2']
    R, target, groups = subtype_shrinkage_rows(names, strength=3.0, parameter_scale=np.full(3, 5.0))
    assert R.shape == (3, 3)
    assert np.allclose(R @ np.full(3, 8.0), target)
    assert np.linalg.norm(R @ np.array([8.0, 9.0, 7.0])) > 0.0
    assert groups['bond:Si-Si'] == [0, 1, 2]


def test_family_mean_prior_targets_only_the_subtype_mean():
    names = ['bond:Si-Si', 'bond:Si-SiH', 'bond:Si-SiH2']
    R, target = family_mean_prior_rows(names, {'bond:Si-Si': 10.0}, strength=2.0)
    assert np.allclose(R @ np.array([8.0, 10.0, 12.0]), target)
    assert np.allclose(R @ np.array([10.0, 10.0, 10.0]), target)


def test_cross_valence_sensitivities_are_symmetric_and_rigid_invariant():
    positions, masses, bonds, angles, B, _ = make_triatomic_problem()
    ss = {('H', 'Si', 'Si'): 0}
    sb = {('H', 'Si', 'Si'): 1}
    A = compute_cross_sensitivity(positions, bonds, angles, ['Si', 'Si', 'H'], ss, sb)
    scale = FFfit.dimensionless_wilson_scale(positions, [('bond', 0, 1), ('bond', 1, 2), ('angle', 0, 1, 2)])
    expected_ss = np.outer(scale[0]*B[0], scale[1]*B[1]) + np.outer(scale[1]*B[1], scale[0]*B[0])
    expected_sb = np.outer((scale[0]*B[0] + scale[1]*B[1])/np.sqrt(2.0), B[2]) + np.outer(B[2], (scale[0]*B[0] + scale[1]*B[1])/np.sqrt(2.0))
    rigid, _ = FFfit.rigid_and_vibrational_bases(positions, masses)
    cartesian_rigid = rigid / np.sqrt(np.repeat(masses, 3))[:, None]
    for p, expected in ((0, expected_ss), (1, expected_sb)):
        assert np.allclose(A[p], expected)
        assert np.allclose(A[p], A[p].T)
        assert np.linalg.norm(A[p] @ cartesian_rigid) < 1e-12


def test_cross_parameter_maps_are_contiguous_when_both_terms_are_enabled():
    angles = [[(0, 1, 2, 1.9), (2, 1, 3, 1.9)]]
    symbols = [['H', 'Si', 'Si', 'H']]
    ss, sb, npar = build_cross_param_maps(angles, symbols, offset=5, stretch_stretch=True, stretch_bend=True)
    assert sorted(list(ss.values()) + list(sb.values())) == list(range(5, npar))


def test_cluster_indices_refer_to_the_original_unsorted_values():
    clusters = cluster_1d(np.array([10.0, 1.0, 1.1, 10.1]), gap_threshold=2.0)
    assert {tuple(sorted(indices.tolist())) for _, _, indices in clusters} == {(0, 3), (1, 2)}


def test_cpp_normal_equation_solver_fails_loudly_when_singular():
    with pytest.raises(np.linalg.LinAlgError, match="singular"):
        FFfit.FFfit.solve_normal_equations(np.array([[1.0, 1.0], [1.0, 1.0]]), np.array([1.0, 1.0]))


def test_regularized_bounded_fit_stays_nonnegative_for_incomplete_model():
    positions, masses, bonds, angles, B, A = make_triatomic_problem()
    missing = np.outer(B[0] + 0.4*B[2], B[1] - 0.3*B[2])
    missing = 0.5 * (missing + missing.T)
    H_ref = 3.0*A[0] + 0.7*A[1] + 0.25*missing
    system = {'A': A, 'H_ref': H_ref, 'positions': positions, 'masses': masses, 'bonds': bonds, 'angles': angles}
    k, diag = FFfit.fit_hybrid_hessian([system], mode_weight=1.0, local_weight=1.0, internal_weight=1.0,
                                       prior=np.array([4.0, 1.0]), parameter_scale=np.array([4.0, 1.0]),
                                       regularization=1e-2, bounds=(0.0, np.inf))
    assert np.all(np.isfinite(k))
    assert np.all(k >= 0.0)
    assert np.isfinite(diag['condition'])
    assert 0.0 < diag['relative_residual'] < 0.2


def test_uff_dihedral_gradient_matches_finite_difference():
    positions = np.array([[0.2, -0.4, 0.7], [0.0, 0.0, 0.0], [1.1, 0.1, -0.2], [1.8, 0.9, 0.6]])
    E, grad = dihedral_energy_gradient(positions, d=1, n=3)
    grad_fd = np.zeros(12)
    h = 1e-6
    for i in range(12):
        pos_p = positions.copy(); pos_p.flat[i] += h
        pos_m = positions.copy(); pos_m.flat[i] -= h
        grad_fd[i] = (dihedral_energy_gradient(pos_p, d=1, n=3)[0] - dihedral_energy_gradient(pos_m, d=1, n=3)[0]) / (2.0*h)
    assert np.isfinite(E)
    assert -np.pi <= dihedral_angle(positions) <= np.pi
    assert np.allclose(grad.ravel(), grad_fd, rtol=2e-6, atol=2e-7)


def test_si_environment_typing_uses_bonded_hydrogen_count():
    symbols = ['Si', 'Si', 'Si', 'Si'] + ['H']*7
    bonds = [(0, 4, 1.0), (0, 5, 1.0), (0, 6, 1.0), (0, 7, 1.0),
             (1, 8, 1.0), (1, 9, 1.0), (2, 10, 1.0), (0, 1, 2.0), (1, 2, 2.0), (2, 3, 2.0)]
    labels = assign_si_environment_types(symbols, bonds, enabled=True)
    assert labels[:4] == ['SiH3', 'SiH2', 'SiH', 'Si']
    assert labels[4:] == ['H']*7
