#!/usr/bin/env python3
import os
import sys
import time

import numpy as np
import pytest

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, _ROOT)

from pyBall.io.crystal_npz import (
    bonds_ij_to_segments, load_crystal_npz, load_npy_crystal, load_topology_npz,
)

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


def test_load_crystal_npz_shapes():
    path = os.path.join(FIXTURE_DIR, '01_init.npz')
    d = load_crystal_npz(path)
    assert d['pos'].shape[1] == 3
    assert d['Z'].shape[0] == d['pos'].shape[0]
    assert d['bonds_ij'] is not None
    assert d['bonds_ij'].shape[1] == 2


def test_load_topology_npz_keys():
    path = os.path.join(FIXTURE_DIR, '03_topology.npz')
    t = load_topology_npz(path)
    for k in ('group_bbox_min', 'group_bbox_max', 'icolGroup'):
        assert k in t
    assert t['group_bbox_min'].shape[0] == t['group_bbox_max'].shape[0]


def test_load_npy_crystal():
    d = load_npy_crystal(os.path.join(FIXTURE_DIR, 'npy_fast'))
    assert d['pos'].shape[0] == d['Z'].shape[0]


def test_crystal_npz_load_timing():
    path = os.path.join(FIXTURE_DIR, 'bench_500.npz')
    t0 = time.perf_counter()
    d = load_crystal_npz(path, mmap=True)
    pos = np.asarray(d['pos'])
    Z = np.asarray(d['Z'])
    bonds = np.asarray(d['bonds_ij'])
    segs = bonds_ij_to_segments(pos, bonds)
    elapsed_ms = (time.perf_counter() - t0) * 1000.0
    assert pos.shape[0] == 500
    assert segs.shape[0] == bonds.shape[0] * 2
    assert elapsed_ms < 200.0, f'load+bond build took {elapsed_ms:.1f} ms'


def test_validate_topology_allows_relaxed_pos_drift():
    from pyBall.io.crystal_npz import validate_topology_crystal_parity
    crystal = {'pos': np.ones((2, 3)), 'Z': np.array([6, 6], dtype=np.int32), 'bonds_ij': np.array([[0, 1]], dtype=np.int32)}
    topology = {'pos': np.zeros((2, 3)), 'Z': crystal['Z'], 'neigh_idx': np.full((2, 4), -1, dtype=np.int32), 'bond_k': np.zeros((2, 4)), 'bond_l0': np.zeros((2, 4)), 'n_bond': np.int32(1)}
    validate_topology_crystal_parity(crystal, topology)


def test_build_hessian_from_linear_topology_2atom():
    from pyBall.FTIR import build_hessian_from_linear_topology
    pos = np.array([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]])
    neigh = np.array([[-1, 1, -1, -1], [-1, -1, -1, -1]], dtype=np.int32)
    k = np.zeros((2, 4)); l0 = np.zeros((2, 4))
    k[0, 1] = 5.0; l0[0, 1] = 1.4
    sc = np.zeros((2, 4), dtype=np.int32); sc[0, 1] = 1
    H = build_hessian_from_linear_topology(pos, neigh, k, l0, stick_class=sc)
    assert H.shape == (6, 6)
    assert np.max(np.abs(H - H.T)) < 1e-10
    w = np.linalg.eigvalsh(H)
    assert w.min() > -1e-10


def test_missing_pos_raises():
    import tempfile
    with tempfile.NamedTemporaryFile(suffix='.npz', delete=False) as f:
        tmp = f.name
    try:
        np.savez_compressed(tmp, Z=np.array([6], dtype=np.int32))
        with pytest.raises(KeyError, match='pos'):
            load_crystal_npz(tmp, mmap=False)
    finally:
        os.unlink(tmp)
