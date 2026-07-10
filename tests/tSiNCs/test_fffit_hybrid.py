import numpy as np

from pyBall import FFfit
from pyBall.FFfit_utils import assign_si_environment_types, dihedral_angle, dihedral_energy_gradient


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
