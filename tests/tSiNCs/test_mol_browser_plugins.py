#!/usr/bin/env python3
"""Tests for VispyMolBrowser plugin system and vibration spectrum helpers."""
import os
import sys

import numpy as np
import pytest

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, _ROOT)

PASSIVATION_DIR = os.path.join(os.path.dirname(__file__), 'fixtures', 'si_1nm_passivation', '02_sphere_13A_H_standard')


def test_is_viewer_listable_basename():
    from pyBall.io.crystal_npz import is_viewer_listable_basename
    assert is_viewer_listable_basename('01_crystal.npz')
    assert is_viewer_listable_basename('02_relaxed.npz')
    assert not is_viewer_listable_basename('04_hessian.npz')
    assert not is_viewer_listable_basename('05_spectrum.npz')
    assert not is_viewer_listable_basename('foo_hessian.npz')


@pytest.mark.skipif(not os.path.isfile(os.path.join(PASSIVATION_DIR, '04_hessian.npz')), reason='passivation fixture missing')
def test_solve_normal_modes_from_hessian_npz():
    from pyBall.nanocrystal_pipeline import solve_normal_modes_from_hessian_npz
    hess = os.path.join(PASSIVATION_DIR, '04_hessian.npz')
    nm = solve_normal_modes_from_hessian_npz(hess)
    assert nm['natoms'] > 0
    assert nm['modes_cart'].shape == (3 * nm['natoms'], 3 * nm['natoms'])
    assert nm['vib_indices'].size > 0
    assert np.all(np.isfinite(nm['omegas_cm_vib']))


@pytest.mark.skipif(not os.path.isfile(os.path.join(PASSIVATION_DIR, '05_spectrum.npz')), reason='passivation fixture missing')
def test_load_spectrum_npz():
    from pyBall.nanocrystal_pipeline import load_spectrum_npz
    spec = load_spectrum_npz(os.path.join(PASSIVATION_DIR, '05_spectrum.npz'))
    assert spec['omega_centers'].size == spec['hist'].size
    assert spec['omegas_modes'].size == 3 * int(np.load(os.path.join(PASSIVATION_DIR, '04_hessian.npz'))['natoms'])


def test_default_plugin_registry():
    from pyBall.GUI.mol_browser_plugins import default_plugin_registry
    reg = default_plugin_registry()
    assert len(reg.plugins) >= 1
    assert reg.plugins[0].plugin_id == 'vibration_spectrum'


@pytest.mark.skipif(not os.path.isdir(PASSIVATION_DIR), reason='passivation fixture missing')
def test_vibration_plugin_relevant():
    from pyBall.GUI.mol_browser_plugins import MolBrowserContext, VibrationSpectrumPlugin
    plugin = VibrationSpectrumPlugin()
    ctx = MolBrowserContext(directory=PASSIVATION_DIR, selected_path=os.path.join(PASSIVATION_DIR, '02_relaxed.npz'), molecule=None, viewer=None)
    assert plugin.is_relevant(ctx)


@pytest.mark.skipif(not os.path.isdir(PASSIVATION_DIR), reason='passivation fixture missing')
def test_scan_hides_analysis_npz():
    from pyBall.GUI.VispyMolBrowser import scan_molecule_entries
    entries = scan_molecule_entries(PASSIVATION_DIR)
    basenames = [os.path.basename(p) for p in entries]
    assert '04_hessian.npz' not in basenames
    assert '05_spectrum.npz' not in basenames
    assert '02_relaxed.npz' in basenames
