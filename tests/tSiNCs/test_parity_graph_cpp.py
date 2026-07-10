"""Parity test: Python reference vs C++ for graph algorithms, Wilson matrix, and dihedral Hessian.

Verifies that the C++ implementations of:
  - bond_graph_distances (BFS all-pairs shortest paths)
  - local_hessian_mask (graph-distance mask)
  - find_3rd_neighbor_bonds (1-4 pairs)
  - enumerate_dihedrals (torsion enumeration)
  - build_wilson_matrix (Wilson B matrix)
  - dihedral_hessian (12x12 FD Hessian)
produce identical results to the Python reference implementations.
"""
import os, sys
import numpy as np
from collections import deque

sys.path.insert(0, os.path.dirname(__file__))
from test_FFfit import (
    shortest_path_distances, build_3rd_neighbor_bonds, build_dihedrals,
    build_topology, dihedral_hessian, dihedral_energy_gradient,
    load_pyscf_case, hessian_ha_bohr_to_ev_ang2,
)
import pyBall.FFfit as FFfit_cpp
from pyBall.FFfit import build_wilson_matrix as py_build_wilson_matrix

RESULTS_DIR = "/home/prokop/SIMULATIONS/jobs_pyscf_vib_OUT_small_nc/results"

def _data_available(name):
    return os.path.exists(os.path.join(RESULTS_DIR, name, 'masses.npy'))

def _load_system(name="Si2H6"):
    """Load a small test system."""
    data = load_pyscf_case(os.path.join(RESULTS_DIR, name), geometry_only=True)
    bonds, angles, bonds3 = build_topology(data['symbols'], data['positions'], 2.5, third_bonds=False)
    return data, bonds, angles


def _make_random_topology(natoms=10, nbonds=15, seed=42):
    """Generate a random connected bond topology for testing."""
    rng = np.random.default_rng(seed)
    bond_pairs = []
    # Ensure connectivity: chain first
    for i in range(natoms - 1):
        bond_pairs.append((i, i + 1))
    # Add random extra bonds
    while len(bond_pairs) < nbonds:
        i, j = rng.integers(0, natoms, size=2)
        if i != j:
            pair = (min(i, j), max(i, j))
            if pair not in bond_pairs:
                bond_pairs.append(pair)
    positions = rng.uniform(-5, 5, size=(natoms, 3))
    return bond_pairs, positions, natoms


class TestBondGraphDistances:
    def test_small_random(self):
        bp, pos, natoms = _make_random_topology()
        bp_arr = np.array(bp, dtype=np.int32).ravel()
        # Python reference
        dist_py = shortest_path_distances(bp, natoms)
        # C++
        dist_cpp = FFfit_cpp.FFfit.bond_graph_distances(bp_arr, natoms)
        np.testing.assert_array_equal(dist_py, dist_cpp)

    def test_disconnected(self):
        natoms = 6
        bp = [(0,1), (1,2), (3,4)]  # 5 is isolated, {0,1,2} and {3,4} components
        bp_arr = np.array(bp, dtype=np.int32).ravel()
        dist_py = shortest_path_distances(bp, natoms)
        dist_cpp = FFfit_cpp.FFfit.bond_graph_distances(bp_arr, natoms)
        np.testing.assert_array_equal(dist_py, dist_cpp)

    def test_real_system(self):
        if not _data_available("Si2H6"): print("Skip: no data"); return
        data, bonds, angles = _load_system("Si2H6")
        natoms = data['natoms']
        bp = [(i, j) for i, j, _ in bonds]
        bp_arr = np.array(bp, dtype=np.int32).ravel()
        dist_py = shortest_path_distances(bp, natoms)
        dist_cpp = FFfit_cpp.FFfit.bond_graph_distances(bp_arr, natoms)
        np.testing.assert_array_equal(dist_py, dist_cpp)


class TestLocalHessianMask:
    def test_small_random(self):
        bp, pos, natoms = _make_random_topology()
        # Python reference
        adj = [[] for _ in range(natoms)]
        for i, j in bp:
            adj[i].append(j); adj[j].append(i)
        max_d = 2
        block = np.zeros((natoms, natoms), dtype=bool)
        for src in range(natoms):
            dist = np.full(natoms, -1, dtype=int)
            dist[src] = 0
            queue = deque([src])
            while queue:
                u = queue.popleft()
                if dist[u] >= max_d: continue
                for v in adj[u]:
                    if dist[v] < 0:
                        dist[v] = dist[u] + 1
                        queue.append(v)
            block[src, dist >= 0] = dist[dist >= 0] <= max_d
        mask_py = np.repeat(np.repeat(block, 3, axis=0), 3, axis=1)
        # C++
        mask_cpp = FFfit_cpp.FFfit.local_hessian_mask(bp, natoms, max_graph_distance=max_d)
        np.testing.assert_array_equal(mask_py, mask_cpp)

    def test_real_system(self):
        if not _data_available("Si2H6"): print("Skip: no data"); return
        data, bonds, angles = _load_system("Si2H6")
        natoms = data['natoms']
        bp = [(i, j) for i, j, _ in bonds]
        # Python
        adj = [[] for _ in range(natoms)]
        for i, j in bp:
            adj[i].append(j); adj[j].append(i)
        max_d = 2
        block = np.zeros((natoms, natoms), dtype=bool)
        for src in range(natoms):
            dist = np.full(natoms, -1, dtype=int)
            dist[src] = 0
            queue = deque([src])
            while queue:
                u = queue.popleft()
                if dist[u] >= max_d: continue
                for v in adj[u]:
                    if dist[v] < 0:
                        dist[v] = dist[u] + 1
                        queue.append(v)
            block[src, dist >= 0] = dist[dist >= 0] <= max_d
        mask_py = np.repeat(np.repeat(block, 3, axis=0), 3, axis=1)
        # C++
        mask_cpp = FFfit_cpp.FFfit.local_hessian_mask(bp, natoms, max_graph_distance=max_d)
        np.testing.assert_array_equal(mask_py, mask_cpp)


class TestFind3rdNeighborBonds:
    def test_small_random(self):
        bp, pos, natoms = _make_random_topology()
        symbols = ['Si'] * natoms
        # Python reference
        bonds3_py = build_3rd_neighbor_bonds(symbols, pos, bp, max_dist=None)
        pairs3_py = set((i, j) for i, j, _ in bonds3_py)
        # C++
        bp_arr = np.array(bp, dtype=np.int32).ravel()
        pairs3_cpp = FFfit_cpp.FFfit.find_3rd_neighbor_bonds(bp_arr, natoms, pos, max_dist=0.0)
        pairs3_cpp_set = set(pairs3_cpp)
        assert pairs3_py == pairs3_cpp_set, f"Mismatch: py={pairs3_py} cpp={pairs3_cpp_set}"

    def test_real_system(self):
        if not _data_available("Si2H6"): print("Skip: no data"); return
        data, bonds, angles = _load_system("Si2H6")
        natoms = data['natoms']
        bp = [(i, j) for i, j, _ in bonds]
        bonds3_py = build_3rd_neighbor_bonds(data['symbols'], data['positions'], bp, max_dist=None)
        pairs3_py = set((i, j) for i, j, _ in bonds3_py)
        bp_arr = np.array(bp, dtype=np.int32).ravel()
        pairs3_cpp = FFfit_cpp.FFfit.find_3rd_neighbor_bonds(bp_arr, natoms, data['positions'], max_dist=0.0)
        pairs3_cpp_set = set(pairs3_cpp)
        assert pairs3_py == pairs3_cpp_set, f"Mismatch: py={pairs3_py} cpp={pairs3_cpp_set}"


class TestEnumerateDihedrals:
    def test_small_random(self):
        bp, pos, natoms = _make_random_topology()
        # Python reference
        diheds_py = build_dihedrals(['C'] * natoms, pos, [(i, j, 1.0) for i, j in bp], dihedral=True)
        tuples_py = set((i, j, k, l) for i, j, k, l, _, _ in diheds_py)
        # C++
        bp_arr = np.array(bp, dtype=np.int32).ravel()
        diheds_cpp = FFfit_cpp.FFfit.enumerate_dihedrals(bp_arr, natoms)
        tuples_cpp = set(diheds_cpp)
        assert tuples_py == tuples_cpp, f"Mismatch: py={tuples_py} cpp={tuples_cpp}"

    def test_real_system(self):
        if not _data_available("Si2H6"): print("Skip: no data"); return
        data, bonds, angles = _load_system("Si2H6")
        natoms = data['natoms']
        bp = [(i, j) for i, j, _ in bonds]
        diheds_py = build_dihedrals(data['symbols'], data['positions'], bonds, dihedral=True)
        tuples_py = set((i, j, k, l) for i, j, k, l, _, _ in diheds_py)
        bp_arr = np.array(bp, dtype=np.int32).ravel()
        diheds_cpp = FFfit_cpp.FFfit.enumerate_dihedrals(bp_arr, natoms)
        tuples_cpp = set(diheds_cpp)
        assert tuples_py == tuples_cpp, f"Mismatch: py={tuples_py} cpp={tuples_cpp}"

    def test_reversed_bond_endpoints_are_equivalent(self):
        bp, pos, natoms = _make_random_topology()
        reversed_bp = [(j, i) for i, j in bp]
        py = build_dihedrals(['C'] * natoms, pos, [(i, j, 1.0) for i, j in reversed_bp], dihedral=True)
        cpp = FFfit_cpp.FFfit.enumerate_dihedrals(np.array(reversed_bp, dtype=np.int32).ravel(), natoms)
        assert set((i, j, k, l) for i, j, k, l, _, _ in py) == set(cpp)


class TestBuildWilsonMatrix:
    def test_real_system(self):
        if not _data_available("Si2H6"): print("Skip: no data"); return
        data, bonds, angles = _load_system("Si2H6")
        natoms = data['natoms']
        positions = data['positions']
        # Python reference
        B_py, labels_py = py_build_wilson_matrix(positions, bonds, angles)
        # C++ — need to set up FFfit instance
        fitter = FFfit_cpp.FFfit()
        fitter.set_geometry(positions)
        for i, j, l0 in bonds:
            fitter.add_bond(i, j, l0)
        for i, j, k, theta0 in angles:
            fitter.add_angle(i, j, k, theta0)
        B_cpp = fitter.build_wilson_matrix_cpp()
        np.testing.assert_allclose(B_py, B_cpp, atol=1e-12)


class TestDihedralHessian:
    def test_random_geometry(self):
        rng = np.random.default_rng(123)
        pos4 = rng.uniform(-3, 3, size=(4, 3))
        # Python reference
        H_py = dihedral_hessian(pos4, d=1, n=3, h=1e-5)
        # C++ — set up FFfit instance with 4 atoms
        fitter = FFfit_cpp.FFfit()
        fitter.set_geometry(pos4)
        H_cpp = fitter.dihedral_hessian_cpp(0, 1, 2, 3)
        np.testing.assert_allclose(H_py, H_cpp, atol=1e-8)

    def test_real_system(self):
        """Test dihedral Hessian on a real Si-H system with a known dihedral."""
        if not _data_available("Si5H12"): print("Skip: no data"); return
        data, bonds, angles = _load_system("Si5H12")
        natoms = data['natoms']
        bp = [(i, j) for i, j, _ in bonds]
        diheds = build_dihedrals(data['symbols'], data['positions'], bonds, dihedral=True)
        if not diheds:
            print("No dihedrals in Si5H12, skipping")
            return
        i, j, k, l, d, n = diheds[0]
        pos4 = data['positions'][[i, j, k, l]]
        # Python
        H_py = dihedral_hessian(pos4, d=d, n=n, h=1e-5)
        # C++
        fitter = FFfit_cpp.FFfit()
        fitter.set_geometry(data['positions'])
        H_cpp = fitter.dihedral_hessian_cpp(i, j, k, l)
        np.testing.assert_allclose(H_py, H_cpp, atol=1e-8)


class TestDihedralBatch:
    def test_batch_vs_single(self):
        """Verify batch dihedral_dHdk_batch matches sum of per-dihedral dHdk calls."""
        rng = np.random.default_rng(42)
        natoms = 10
        positions = rng.uniform(-5, 5, size=(natoms, 3))
        # Build a simple bond graph: chain 0-1-2-3-4-5-6-7-8-9
        bonds = [(i, i+1, 1.5) for i in range(natoms - 1)]
        bp = [(i, i+1) for i in range(natoms - 1)]
        diheds = build_dihedrals(['Si']*natoms, positions, bonds, dihedral=True)
        if not diheds:
            print("No dihedrals generated, skipping")
            return
        fitter = FFfit_cpp.FFfit()
        fitter.set_geometry(positions)
        # Per-dihedral accumulation
        n3 = natoms * 3
        A_single = np.zeros((n3, n3), dtype=np.float64)
        for (i, j, k, l, d, n) in diheds:
            A_single += fitter.dihedral_dHdk_cpp(i, j, k, l)
        # Batch
        dh_flat = np.array([(i, j, k, l) for (i, j, k, l, d, n) in diheds], dtype=np.int32)
        A_batch = fitter.dihedral_dHdk_batch_cpp(dh_flat)
        np.testing.assert_allclose(A_single, A_batch, atol=1e-12)
        print(f"  Batch vs single: {len(diheds)} dihedrals, max diff={np.max(np.abs(A_single - A_batch)):.2e}")

    def test_batch_typed_vs_python(self):
        """Verify batch_typed matches Python compute_dihedral_sensitivity."""
        rng = np.random.default_rng(99)
        natoms = 8
        positions = rng.uniform(-5, 5, size=(natoms, 3))
        symbols = ['Si'] * natoms
        bonds = [(i, i+1, 1.5) for i in range(natoms - 1)]
        diheds = build_dihedrals(symbols, positions, bonds, dihedral=True)
        if not diheds:
            print("No dihedrals, skipping")
            return
        # Simple type map: all dihedrals share one type (Si-Si-Si-Si)
        from pyBall.FFfit_utils import dihedral_type_key
        type_map = {}
        type_idx_list = []
        for (i, j, k, l, d, n) in diheds:
            key = dihedral_type_key(symbols, i, j, k, l)
            if key not in type_map:
                type_map[key] = len(type_map)
            type_idx_list.append(type_map[key])
        n_types = len(type_map)
        # Python reference
        from pyBall.FFfit_utils import compute_dihedral_sensitivity
        A_py = compute_dihedral_sensitivity(positions, symbols, diheds, type_map, h=1e-5)
        # C++ batch
        fitter = FFfit_cpp.FFfit()
        fitter.set_geometry(positions)
        dh_flat = np.array([(i, j, k, l) for (i, j, k, l, d, n) in diheds], dtype=np.int32)
        ti = np.array(type_idx_list, dtype=np.int32)
        A_cpp = fitter.dihedral_dHdk_batch_typed_cpp(dh_flat, ti, n_types)
        # Compare
        for p in range(n_types):
            np.testing.assert_allclose(A_py[p], A_cpp[p], atol=1e-6)
        print(f"  Batch typed vs Python: {len(diheds)} dihedrals, {n_types} types, max diff={max(np.max(np.abs(A_py[p] - A_cpp[p])) for p in range(n_types)):.2e}")


if __name__ == "__main__":
    # Run all tests
    import traceback
    test_classes = [TestBondGraphDistances, TestLocalHessianMask,
                    TestFind3rdNeighborBonds, TestEnumerateDihedrals,
                    TestBuildWilsonMatrix, TestDihedralHessian,
                    TestDihedralBatch]
    n_pass = 0; n_fail = 0
    for cls in test_classes:
        for method_name in dir(cls):
            if method_name.startswith("test_"):
                try:
                    getattr(cls(), method_name)()
                    print(f"  PASS: {cls.__name__}.{method_name}")
                    n_pass += 1
                except Exception as e:
                    print(f"  FAIL: {cls.__name__}.{method_name}: {e}")
                    traceback.print_exc()
                    n_fail += 1
    print(f"\n=== {n_pass} passed, {n_fail} failed ===")
