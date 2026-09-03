#!/usr/bin/env python3
"""Fail-loud: a harmonic spectrum is the Hessian at that energy's own minimum. No sqrt(max(λ,0))."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_ROOT))

from pyBall.FFfit_utils import (
    assert_harmonic_spectrum_at_minimum, cartesian_modes_from_hessian, signed_frequencies_cm1,
)

L2 = _ROOT / "tests/tSiNCs/OUT_dftb_vs_mmff"
WRONG = L2 / "WRONG_at_DFTB_geometry"
CASES = ("octahedron_C", "truncated_octahedron_C", "rhombic_dodecahedron_C")


def _ok_spectrum(n=30, n_rigid=6):
    w = np.linspace(80.0, 2200.0, n)
    w[:n_rigid] = np.array([0.0, 0.01, -0.02, 0.03, -0.04, 0.05])[:n_rigid]
    return w


def test_assert_passes_six_rigid_and_positive_vibs():
    st = assert_harmonic_spectrum_at_minimum(_ok_spectrum())
    if st["n_imag"] != 0 or st["n_rigid"] != 6:
        raise AssertionError(st)


def test_assert_fails_large_imaginary():
    w = _ok_spectrum()
    w[10] = -80.0
    with pytest.raises(RuntimeError, match="imaginary modes"):
        assert_harmonic_spectrum_at_minimum(w)


def test_assert_fails_missing_rotations():
    """WRONG default MMFF at DFTB q: 0 large imag, only 3 rigid, first vib ~18 cm⁻¹."""
    w = _ok_spectrum()
    w[:6] = [0.0, 0.0, 0.0, 18.0, 22.0, 30.0]
    with pytest.raises(RuntimeError, match="rigid-body"):
        assert_harmonic_spectrum_at_minimum(w)


def test_complex_imaginary_as_negative():
    w = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 200.0, 0.0 + 80.0j])
    with pytest.raises(RuntimeError, match="imaginary modes"):
        assert_harmonic_spectrum_at_minimum(w)


def test_indefinite_hessian_cartesian_modes_fails():
    n = 4
    masses = np.ones(n)
    H = np.eye(3 * n)
    H[0, 0] = -2.0
    with pytest.raises(RuntimeError, match="imaginary modes"):
        cartesian_modes_from_hessian(H, masses)


def test_signed_frequencies_keep_negative_curvature():
    n = 3
    masses = np.ones(n)
    H = np.diag(np.concatenate([[-1.0], np.ones(8)]))
    om = signed_frequencies_cm1(H, masses)
    if not np.any(om < -10.0):
        raise AssertionError(f"signed frequencies hid the negative eigenvalue: {om}")


@pytest.mark.parametrize("name", CASES)
def test_l2_ownmin_mmff_is_a_spectrum(name):
    p = L2 / name / "mmff_ownmin_protocol.json"
    if not p.is_file():
        pytest.skip(f"no own-min protocol for {name}")
    for kind in ("default", "scaled"):
        f = np.load(L2 / name / f"mmff_{kind}_frequencies_cm1.npy")
        assert_harmonic_spectrum_at_minimum(f, ctx=f"{name} MMFF {kind} own-min: ")


@pytest.mark.parametrize("name", CASES)
def test_wrong_at_dftb_geometry_is_not_a_spectrum(name):
    d = WRONG / name
    if not (d / "mmff_scaled_frequencies_cm1.npy").is_file():
        pytest.skip(f"no quarantine files for {name}")
    f_s = np.load(d / "mmff_scaled_frequencies_cm1.npy")
    with pytest.raises(RuntimeError, match="imaginary modes|rigid-body"):
        assert_harmonic_spectrum_at_minimum(f_s, ctx=f"{name} WRONG scaled: ")
    f_d = np.load(d / "mmff_default_frequencies_cm1.npy")
    with pytest.raises(RuntimeError, match="rigid-body"):
        assert_harmonic_spectrum_at_minimum(f_d, ctx=f"{name} WRONG default: ")
