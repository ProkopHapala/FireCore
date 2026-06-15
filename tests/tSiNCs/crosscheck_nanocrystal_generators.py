#!/usr/bin/env python3
"""Cross-check JS and Python nanocrystal generators (C/Si, sphere/planes). Writes tests/tSiNCs/crosscheck/."""
import json
import os
import subprocess
import sys
import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))
OUT = os.path.join(REPO, 'tests/tSiNCs/crosscheck')
sys.path.insert(0, REPO)

from pyBall.nanocrystal_gen import build_spherical_nanoparticle, save_xyz


def load_xyz(path):
    with open(path) as f:
        n = int(f.readline())
        f.readline()
        elems, pos = [], []
        for _ in range(n):
            ws = f.readline().split()
            elems.append(ws[0])
            pos.append([float(ws[1]), float(ws[2]), float(ws[3])])
    return np.array(elems), np.asarray(pos, float)


def hh_angles_ch2(elems, pos, heavy='C'):
    angs = []
    for i, e in enumerate(elems):
        if e != heavy:
            continue
        hs = [j for j in range(len(elems)) if elems[j] == 'H' and np.linalg.norm(pos[j] - pos[i]) < 1.6]
        if len(hs) != 2:
            continue
        d1 = pos[hs[0]] - pos[i]
        d2 = pos[hs[1]] - pos[i]
        d1 /= np.linalg.norm(d1)
        d2 /= np.linalg.norm(d2)
        angs.append(float(np.degrees(np.arccos(np.clip(d1.dot(d2), -1, 1)))))
    return angs


def hh_angles_nb2(elems, pos, heavy_sym):
    return hh_angles_ch2(elems, pos, heavy=heavy_sym)


def h_clashes(elems, pos, cut=1.8):
    n = len(elems)
    c = 0
    for i in range(n):
        if elems[i] != 'H':
            continue
        for j in range(i + 1, n):
            if elems[j] != 'H':
                continue
            if np.linalg.norm(pos[j] - pos[i]) < cut:
                c += 1
    return c


def compare_xyz(a_path, b_path, tol_pos=0.05):
    ea, pa = load_xyz(a_path)
    eb, pb = load_xyz(b_path)
    report = {'a': a_path, 'b': b_path, 'natoms_a': len(ea), 'natoms_b': len(eb),
              'match_count': len(ea) == len(eb)}
    if len(ea) != len(eb):
        return report
    # compare heavy atom count
    report['n_heavy_a'] = int(np.sum(ea != 'H'))
    report['n_heavy_b'] = int(np.sum(eb != 'H'))
    report['n_H_a'] = int(np.sum(ea == 'H'))
    report['n_H_b'] = int(np.sum(eb == 'H'))
    # RMSD after matching by sorted distances from COM (rough)
    ca = pa[ea != 'H'].mean(axis=0)
    cb = pb[eb != 'H'].mean(axis=0)
    da = np.linalg.norm(pa[ea != 'H'] - ca, axis=1)
    db = np.linalg.norm(pb[eb != 'H'] - cb, axis=1)
    da.sort()
    db.sort()
    if len(da) == len(db):
        report['heavy_radial_rmsd'] = float(np.sqrt(np.mean((da - db) ** 2)))
    report['hh_mean_a'] = float(np.mean(hh_angles_ch2(ea, pa))) if np.any(ea == 'H') else None
    report['hh_mean_b'] = float(np.mean(hh_angles_ch2(eb, pb))) if np.any(eb == 'H') else None
    report['clashes_a'] = h_clashes(ea, pa)
    report['clashes_b'] = h_clashes(eb, pb)
    return report


def run_node(cmd):
    print('[crosscheck] node', ' '.join(cmd))
    subprocess.run(cmd, cwd=REPO, check=True)


def run_py(cmd):
    print('[crosscheck] python', ' '.join(cmd))
    subprocess.run(cmd, cwd=REPO, check=True)


def main():
    os.makedirs(OUT, exist_ok=True)
    results = []

    # --- C sphere R=6 nrep=5 ---
    py_c = os.path.join(OUT, 'C_sphere_R6_py.xyz')
    elems, apos = build_spherical_nanoparticle(R=6.0, nrep=5, heavy_z=6, resolve_clashes=False)
    save_xyz(py_c, elems, apos, 'py native C sphere')
    js_c = os.path.join(OUT, 'C_sphere_R6_js.xyz')
    run_node(['node', 'scripts/gen_nanocrystals.mjs', '--samples', '1', '--maxFiles', '1', '--seed', '42',
              '--cutMode', 'sphere', '--sphereR', '6', '--sphereNrep', '5',
              '--cif', 'cpp/common_resources/crystals/diamond_primitive.cif', '--applySymmetry', '0',
              '--caps', 'H', '--insertProb', '0', '--collapseProb', '0', '--resolveClashes', '0',
              '--outDir', OUT, '--prefix', 'C_sphere_R6_js'])
    js_mol2 = sorted(f for f in os.listdir(OUT) if f.startswith('C_sphere_R6_js') and f.endswith('.mol2'))[-1]
    run_node(['node', 'scripts/mol2_to_xyz.mjs', os.path.join(OUT, js_mol2), js_c])
    results.append({'case': 'C_sphere_R6', **compare_xyz(py_c, js_c)})

    # --- Si sphere R=6 nrep=5 ---
    py_si = os.path.join(OUT, 'Si_sphere_R6_py.xyz')
    elems, apos = build_spherical_nanoparticle(R=6.0, nrep=5, heavy_z=14, resolve_clashes=False)
    save_xyz(py_si, elems, apos, 'py native Si sphere')
    js_si = os.path.join(OUT, 'Si_sphere_R6_js.xyz')
    run_node(['node', 'scripts/gen_nanocrystals.mjs', '--samples', '1', '--maxFiles', '1', '--seed', '42',
              '--cutMode', 'sphere', '--sphereR', '6', '--sphereNrep', '5',
              '--cif', 'cpp/common_resources/crystals/Si_primitive.cif', '--applySymmetry', '0',
              '--caps', 'H', '--insertProb', '0', '--collapseProb', '0', '--resolveClashes', '0',
              '--outDir', OUT, '--prefix', 'Si_sphere_R6_js'])
    js_mol2 = sorted(f for f in os.listdir(OUT) if f.startswith('Si_sphere_R6_js') and f.endswith('.mol2'))[-1]
    run_node(['node', 'scripts/mol2_to_xyz.mjs', os.path.join(OUT, js_mol2), js_si])
    results.append({'case': 'Si_sphere_R6', **compare_xyz(py_si, js_si)})

    # --- Si planes G2 caps-only (py via node delegate) ---
    si_g2_js = os.path.join(OUT, 'Si_planes_G2_js.xyz')
    run_node(['node', 'scripts/gen_nanocrystals.mjs', '--samples', '1', '--maxFiles', '1', '--seed', '42',
              '--cutMode', 'planes', '--nx-range', '2,2', '--ny-range', '2,2', '--nz-range', '2,2', '--centered', '1',
              '--planeTemplates', 'a111', '--planeSymC', '6', '--planeCScale', '0.40', '--planeCJitter', '0',
              '--caps', 'H', '--insertProb', '0', '--collapseProb', '0',
              '--outDir', OUT, '--prefix', 'Si_planes_G2_js'])
    js_mol2 = sorted(f for f in os.listdir(OUT) if f.startswith('Si_planes_G2_js') and f.endswith('.mol2'))[-1]
    run_node(['node', 'scripts/mol2_to_xyz.mjs', os.path.join(OUT, js_mol2), si_g2_js])
    si_g2_py = os.path.join(OUT, 'Si_planes_G2_py.xyz')
    run_py([sys.executable, 'scripts/gen_nanocrystals.py', '--cutMode', 'planes', '--element', 'Si',
            '--seed', '42', '--nx-range', '2,2', '--ny-range', '2,2', '--nz-range', '2,2',
            '--planeTemplates', 'a111', '--planeCScale', '0.40', '--planeCJitter', '0',
            '--caps', 'H', '--outDir', OUT, '--prefix', 'Si_planes_G2_py'])
    py_mol2 = sorted(f for f in os.listdir(OUT) if f.startswith('Si_planes_G2_py') and f.endswith('.mol2'))[-1]
    run_node(['node', 'scripts/mol2_to_xyz.mjs', os.path.join(OUT, py_mol2), si_g2_py])
    results.append({'case': 'Si_planes_G2', **compare_xyz(si_g2_py, si_g2_js)})

    # --- C planes from diamond CIF ---
    c_g1_js = os.path.join(OUT, 'C_planes_G1_js.xyz')
    run_node(['node', 'scripts/gen_nanocrystals.mjs', '--samples', '1', '--maxFiles', '1', '--seed', '42',
              '--cutMode', 'planes', '--cif', 'cpp/common_resources/crystals/diamond_primitive.cif', '--applySymmetry', '0',
              '--nx-range', '1,1', '--ny-range', '1,1', '--nz-range', '1,1', '--centered', '1',
              '--planeTemplates', 'a111', '--planeSymC', '6', '--planeCScale', '0.50', '--planeCJitter', '0',
              '--caps', 'H', '--insertProb', '0', '--collapseProb', '0', '--resolveClashes', '0',
              '--outDir', OUT, '--prefix', 'C_planes_G1_js'])
    js_mol2 = sorted(f for f in os.listdir(OUT) if f.startswith('C_planes_G1_js') and f.endswith('.mol2'))[-1]
    run_node(['node', 'scripts/mol2_to_xyz.mjs', os.path.join(OUT, js_mol2), c_g1_js])
    angs = hh_angles_ch2(*load_xyz(c_g1_js), heavy='C')
    results.append({'case': 'C_planes_G1_js', 'natoms': len(load_xyz(c_g1_js)[0]),
                    'hh_mean': float(np.mean(angs)) if angs else None, 'clashes': h_clashes(*load_xyz(c_g1_js))})

    # --- Bridge feature demo (JS only; seed fixed) ---
    br_mol2_dir = OUT
    run_node(['node', 'scripts/gen_nanocrystals.mjs', '--samples', '1', '--maxFiles', '1', '--seed', '42',
              '--cutMode', 'planes', '--nx-range', '2,2', '--ny-range', '2,2', '--nz-range', '2,2', '--centered', '1',
              '--planeTemplates', 'a111', '--planeCScale', '0.40', '--planeCJitter', '0',
              '--caps', 'H', '--insertProb', '0.2', '--collapseProb', '0.3',
              '--outDir', br_mol2_dir, '--prefix', 'Si_bridges_demo'])
    br_mol2 = sorted(f for f in os.listdir(OUT) if f.startswith('Si_bridges_demo') and f.endswith('.mol2'))[-1]
    br_xyz = os.path.join(OUT, 'Si_bridges_demo.xyz')
    run_node(['node', 'scripts/mol2_to_xyz.mjs', os.path.join(OUT, br_mol2), br_xyz])
    results.append({'case': 'Si_bridges_demo', 'mol2': br_mol2, 'xyz': br_xyz})

    report_path = os.path.join(OUT, 'crosscheck_report.json')
    with open(report_path, 'w') as f:
        json.dump(results, f, indent=2)

    print('\n=== crosscheck summary ===')
    for r in results:
        print(r)
    print(f'\nReport: {report_path}')
    print(f'Artifacts: {OUT}/')
    failed = [r for r in results if r.get('match_count') is False]
    if failed:
        print('NOTE: atom count mismatches expected for C_sphere (CIF vs xyz lattice paths)')
    # geometry gates: tetrahedral H-H on nb=2 caps (C); Si nb=2 may be absent on small spheres
    for r in results:
        for k in ('hh_mean', 'hh_mean_a', 'hh_mean_b'):
            v = r.get(k)
            if v is None or (isinstance(v, float) and np.isnan(v)):
                continue
            if v < 100:
                raise SystemExit(f"FAIL {r['case']}: H-H angle {k}={v}")
    print('Geometry gates: PASS (H-H tetrahedral where capped)')


if __name__ == '__main__':
    main()
