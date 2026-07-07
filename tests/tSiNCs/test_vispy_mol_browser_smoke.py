#!/usr/bin/env python3
import os
import sys

import pytest

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, _ROOT)

FIXTURE_DIR = os.path.join(os.path.dirname(__file__), 'fixtures', 'npz_viewer')


@pytest.fixture(scope='module', autouse=True)
def ensure_fixtures():
    if not os.path.isfile(os.path.join(FIXTURE_DIR, '01_init.npz')):
        import subprocess
        mjs = os.path.join(FIXTURE_DIR, 'bootstrap_fixtures.mjs')
        if os.path.isfile(mjs):
            subprocess.check_call(['node', mjs], cwd=FIXTURE_DIR)
        else:
            subprocess.check_call([sys.executable, os.path.join(FIXTURE_DIR, 'bootstrap_fixtures.py')])


def test_import_vispy_mol_browser():
    from pyBall.GUI.VispyMolBrowser import VispyMolBrowser, load_molecule_entry  # noqa: F401


def test_load_molecule_entry_offscreen():
    os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
    from pyBall.GUI.VispyMolBrowser import load_molecule_entry, scan_molecule_entries
    entries = scan_molecule_entries(FIXTURE_DIR)
    assert len(entries) >= 2
    mol = load_molecule_entry(os.path.join(FIXTURE_DIR, '01_init.npz'))
    assert mol.natoms > 0
    assert mol.bonds_ij is not None
