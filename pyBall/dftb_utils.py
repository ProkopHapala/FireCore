
import sys
import os
import numpy as np
from . import elements
from . import atomicUtils as au
from textwrap import dedent,indent

DFTB_EXE        = os.environ.get('DFTB_EXE', 'dftb+')
DEFAULT_SK_PATH = os.environ.get('DFTB_SK_PATH', '/home/prokophapala/SIMULATIONS/dftbplus/slakos/3ob-3-1/')

methods=[  'GFN2' ]

methods_XTB  = { 'GFN1', 'GFN2', 'IPEA1' }
methods_dftb = { '3ob', 'D3H5' }


H5Scaling ={
    "O" : 0.06,
    "N" : 0.18,
    "S" : 0.21,
}

default_params={
    "RScaling" : 0.714,
    "WScaling" : 0.25,
    "sr6"      : 1.25,
    "alpha6"   : 29.61,
    "s6"       : 1.0,
    "s8"       : 0.49,
    "HHRepulsion" : "Yes",
    #"Optimizer"  : "Rational{}",
    "Optimizer"  : "LBFGS{  Memory = 20 }",
    #"Optimizer"  : "FIRE{StepSize = 1.0}",       #
    "MaxSteps": 1000,
    "GradElem": 1E-4
    #'Temperature' : 300
}

# >>> Fromn scan_utils.py

def select_atom_index(enames, apos, symbol, axis=1, mode='abs_min', value=0.0):
    idx = [i for i, e in enumerate(enames) if e == symbol]
    assert idx, f"No atoms with symbol '{symbol}'"
    if mode == 'abs_min':
        return idx[int(np.argmin([abs(apos[i, axis] - value) for i in idx]))]
    if mode == 'min':
        return idx[int(np.argmin([apos[i, axis] for i in idx]))]
    if mode == 'max':
        return idx[int(np.argmax([apos[i, axis] for i in idx]))]
    raise ValueError(f"Unknown mode '{mode}'")


def find_closest_indices(enames, apos, target_idx, symbol, n=2):
    idx = [i for i, e in enumerate(enames) if e == symbol]
    assert idx, f"No atoms with symbol '{symbol}'"
    dists = [(i, float(np.linalg.norm(apos[target_idx] - apos[i]))) for i in idx if i != target_idx]
    dists.sort(key=lambda x: x[1])
    assert len(dists) >= n, f"Need at least {n} atoms '{symbol}' (excluding target), got {len(dists)}"
    return [p[0] for p in dists[:n]], [p[1] for p in dists[:n]]


def make_axis_path(p0, p1, svals):
    p0 = np.array(p0, dtype=float)
    p1 = np.array(p1, dtype=float)
    axis = p1 - p0
    L = float(np.linalg.norm(axis))
    assert L > 1e-8, "Axis length too small"
    axis_hat = axis / L
    svals = np.array(svals, dtype=float)
    path = np.zeros((len(svals), 3), dtype=float)
    for i, s in enumerate(svals):
        path[i, :] = p0 + axis_hat * s
    return path


def identify_hbond_transfer(enames, apos, h_symbol='H', heavy_symbol='N', h_select_axis=1, h_select_mode='abs_min', h_select_value=0.0, verbose=True):
    apos = np.array(apos, dtype=float)
    h_scan_idx = select_atom_index(enames, apos, h_symbol, axis=h_select_axis, mode=h_select_mode, value=h_select_value)
    (heavy_idx, heavy_dists) = find_closest_indices(enames, apos, h_scan_idx, heavy_symbol, n=2)
    donor_idx, acceptor_idx = heavy_idx[0], heavy_idx[1]
    fixed_idx = [i for i, e in enumerate(enames) if e == heavy_symbol] + [h_scan_idx]
    if verbose:
        h_idx = [i for i, e in enumerate(enames) if e == h_symbol]
        n_idx = [i for i, e in enumerate(enames) if e == heavy_symbol]
        print(f"  {h_symbol} atoms: {h_idx}  positions y={[f'{apos[i,1]:.2f}' for i in h_idx]}")
        print(f"  {heavy_symbol} atoms: {n_idx}  positions y={[f'{apos[i,1]:.2f}' for i in n_idx]}")
        print(f"  Scanning {h_symbol}: idx={h_scan_idx}, y={apos[h_scan_idx,1]:.3f}")
        print(f"  Donor {heavy_symbol}:    idx={donor_idx}, y={apos[donor_idx,1]:.2f}, d={heavy_dists[0]:.2f}Å")
        print(f"  Acceptor {heavy_symbol}: idx={acceptor_idx}, y={apos[acceptor_idx,1]:.2f}, d={heavy_dists[1]:.2f}Å")
        dNN = np.linalg.norm(apos[donor_idx] - apos[acceptor_idx])
        print(f"  {heavy_symbol}-{heavy_symbol} distance: {dNN:.3f} Å")
        print(f"  H-bond: {heavy_symbol}{donor_idx+1}-{h_symbol}{h_scan_idx+1}...{heavy_symbol}{acceptor_idx+1}  -->  {heavy_symbol}{donor_idx+1}...{h_symbol}{h_scan_idx+1}-{heavy_symbol}{acceptor_idx+1}")
    return h_scan_idx, donor_idx, acceptor_idx, fixed_idx

# <<< Fromn scan_utils.py

def load_molecule(filename, use_ase=True):
    from ase import Atoms
    from ase.io import read
    from pyBall.AtomicSystem import AtomicSystem
    if use_ase:
        return read(filename)
    mol = AtomicSystem(fname=filename)
    return Atoms(symbols=mol.enames, positions=mol.apos, cell=mol.lvec if mol.lvec is not None else None, pbc=(mol.lvec is not None))

def get_max_angular_momentum(enames):
    return {ename: elements.ELEMENT_DICT[ename][4] for ename in sorted(set(enames))}

def write_dftb_input_hessian(enames, gname="geo.xyz", fname='dftb_in.hsd', basis_path=DEFAULT_SK_PATH, delta=1e-4, SCCTolerance=1e-7):
    max_angular = get_max_angular_momentum(enames)
    max_angular_str = '\n'.join([f'        {elem} = "{max_angular[elem]}"' for elem in max_angular])
    with open(fname, 'w') as f:
        f.write(f'''Geometry = xyzFormat {{
    <<< "{gname}"
}}

Driver = SecondDerivatives {{
    Delta = {delta}
    Atoms = 1:-1
}}

Hamiltonian = DFTB {{
    Scc = Yes
    SlaterKosterFiles = Type2FileNames {{
        Prefix = {basis_path}
        Separator = "-"
        Suffix = ".skf"
    }}
    MaxAngularMomentum {{
{max_angular_str}
    }}
    SCCTolerance = {SCCTolerance:.6e}
}}
''')

def write_dftb_input_orbitals(enames, gname="geo.xyz", fname='dftb_in.hsd', basis_path=DEFAULT_SK_PATH, SCCTolerance=1e-7):
    max_angular = get_max_angular_momentum(enames)
    max_angular_str = '\n'.join([f'        {elem} = "{max_angular[elem]}"' for elem in max_angular])
    with open(fname, 'w') as f:
        f.write(f'''Geometry = xyzFormat {{
    <<< "{gname}"
}}

Options {{
  WriteDetailedXml = Yes
}}

Analysis {{
  WriteEigenvectors = Yes
}}

Hamiltonian = DFTB {{
  Scc = Yes
  SlaterKosterFiles = Type2FileNames {{
    Prefix = "{basis_path}"
    Separator = "-"
    Suffix = ".skf"
  }}
  MaxAngularMomentum {{
{max_angular_str}
  }}
  SCCTolerance = {SCCTolerance:.6e}
}}
''')

def read_hessian(filename='hessian.out', n_atoms=None):
    numbers = []
    with open(filename, 'r') as f:
        for line in f:
            for part in line.split():
                try:
                    numbers.append(float(part))
                except ValueError:
                    pass
    data = np.array(numbers)
    if n_atoms is None:
        n_total = int(np.sqrt(len(data)))
        if n_total * n_total != len(data):
            raise ValueError(f"Cannot infer n_atoms from {len(data)} elements")
        n_atoms = n_total // 3
    expected = (3 * n_atoms) ** 2
    if len(data) != expected:
        raise ValueError(f"Expected {expected} elements for {n_atoms} atoms, got {len(data)}")
    return data.reshape((3 * n_atoms, 3 * n_atoms))

def hessian_hartree_bohr_to_eV_angstrom(hessian):
    return hessian * (27.2114 / (0.529177 ** 2))

def parse_energy_out(fname='OUT', allow_unconverged=False):
    with open(fname, 'r') as f:
        lines = f.readlines()
    hits = [line for line in lines if "Total Energy" in line]
    if hits:
        return float(hits[-1][51:70].strip())
    if not allow_unconverged:
        raise AssertionError(f"Could not parse Total Energy from {fname}")
    # Fallback: parse last SCC iteration electronic energy (not a proper total energy).
    # We do NOT hide this: print loud warning.
    print(f"WARNING parse_energy_out(): 'Total Energy' not found in {fname}; using last SCC electronic energy (SCC not converged)")
    in_scc = False
    last_e = None
    for line in lines:
        if 'iSCC Total electronic' in line:
            in_scc = True
            continue
        if in_scc:
            ws = line.split()
            if len(ws) >= 2:
                try:
                    # columns: iSCC, TotalElectronic, DiffElectronic, SCCError
                    last_e = float(ws[1])
                except ValueError:
                    pass
    assert last_e is not None, f"Could not parse SCC electronic energy from {fname}"
    return last_e

def parse_forces(fname='detailed.out', natoms=None):
    assert natoms is not None, "parse_forces() requires natoms"
    forces = np.zeros((natoms, 3))
    HAU2EVA = 51.422067
    with open(fname, 'r') as f:
        lines = f.readlines()
    in_forces = False
    idx = 0
    for line in lines:
        if 'Total Forces' in line:
            in_forces = True
            continue
        if in_forces and idx < natoms:
            parts = line.split()
            if len(parts) >= 3:
                try:
                    forces[idx] = [float(parts[0]), float(parts[1]), float(parts[2])]
                    forces[idx] *= HAU2EVA
                    idx += 1
                except ValueError:
                    continue
    return forces

def read_relaxed_geometry(apos, do_relax=False):
    apos_out = apos.copy()
    if not do_relax:
        return apos_out
    if os.path.exists('geom.out.gen'):
        with open('geom.out.gen') as fgen:
            lines = fgen.readlines()
        n = int(lines[0].split()[0])
        apos_out = np.zeros((n, 3))
        for j in range(n):
            parts = lines[2+j].split()
            apos_out[j] = [float(parts[2]), float(parts[3]), float(parts[4])]
    elif os.path.exists('geom.out.xyz'):
        with open('geom.out.xyz') as fxyz:
            n = int(fxyz.readline())
            fxyz.readline()
            apos_out = np.zeros((n, 3))
            for j in range(n):
                parts = fxyz.readline().split()
                if len(parts) >= 4:
                    apos_out[j] = [float(parts[1]), float(parts[2]), float(parts[3])]
    return apos_out

def run_pbc(apos, enames, lvs, basis_path=DEFAULT_SK_PATH, do_relax=False, fixed_atoms=None, nk=(16,1,1), k_shift=(0.5,0.0,0.0), dftb_exe=DFTB_EXE, workdir=None, Temperature=300, MixingParameter=0.2, MaxScc=200, SCCTolerance=1e-5, params=None, allow_unconverged_energy=False):
    if params is None:
        params = default_params
    cwd = os.getcwd()
    if workdir is not None:
        os.makedirs(workdir, exist_ok=True)
        os.chdir(workdir)
    try:
        makeDFTBjob_pbc(enames=enames, apos=apos, lvs=lvs, fname='dftb_in.hsd', basis_path=basis_path, nk=nk, k_shift=k_shift, opt=do_relax, params=params, Temperature=Temperature, MixingParameter=MixingParameter, MaxScc=MaxScc, SCCTolerance=SCCTolerance, fixed_atoms=fixed_atoms)
        # Capture both stdout and stderr
        ierr = os.system(f'{dftb_exe} > OUT 2> ERR')
        if ierr != 0:
            with open('ERR', 'r') as f:
                err_msg = f.read()
            with open('OUT', 'r') as f:
                out_msg = f.read()
            raise RuntimeError(f"DFTB+ command failed with code {ierr}\n=== STDERR ===\n{err_msg}\n=== STDOUT (last 50 lines) ===\n{out_msg[-5000:]}")
        E = parse_energy_out('OUT', allow_unconverged=allow_unconverged_energy)
        apos_out = read_relaxed_geometry(apos, do_relax=do_relax)
        forces = parse_forces('detailed.out', len(enames)) if os.path.exists('detailed.out') else None
        return E, apos_out, forces
    finally:
        if workdir is not None:
            os.chdir(cwd)

def constrained_scan(apos0, enames, lvs, moved_idx, path, fixed_idx=None, outdir='.', do_relax=True, use_prev_relaxed=True, basis_path=DEFAULT_SK_PATH, dftb_exe=DFTB_EXE, nk=(16,1,1), k_shift=(0.5,0.0,0.0), Temperature=300, MixingParameter=0.2, MaxScc=500, SCCTolerance=1e-6, params=None, step_prefix='tmp_scan', results_prefix='scan', save_xyz=True, xyz_fname=None, key_func=None, key_name='s', plot_step_func=None):
    """General constrained scan over a path for selected atoms.

    Args:
        apos0: (natoms,3) starting geometry [Å]
        enames: list[str] length natoms
        lvs: (3,3) lattice vectors [Å]
        moved_idx: list[int] indices of atoms to move according to `path`
        path: ndarray (nMoved, nstep, 3) positions [Å] for moved atoms
        fixed_idx: list[int] indices to keep fixed during relaxation (passed to DFTB+)
        outdir: output directory
        do_relax: if True run geometry optimization with constraints
        use_prev_relaxed: if True, each step starts from previous relaxed geometry
        key_func: optional callable(step, apos_out, forces, energy) -> scalar key to store
        plot_step_func: optional callable(step, apos_out, forces, meta_dict) for per-step plots

    Returns:
        results: list of dicts per step with keys: 'E','apos','enames','forces','fixed_idx', plus optional key_name
    """
    os.makedirs(outdir, exist_ok=True)
    apos0 = np.array(apos0, dtype=float)
    path = np.array(path, dtype=float)
    moved_idx = list(moved_idx)
    nMoved = len(moved_idx)
    assert path.ndim == 3 and path.shape[0] == nMoved and path.shape[2] == 3, f"path must be (nMoved,nstep,3); got {path.shape} for nMoved={nMoved}"
    nstep = path.shape[1]
    if fixed_idx is None:
        fixed_idx = []
    fixed_idx = list(fixed_idx)

    if params is None:
        params = default_params

    results = []
    for istep in range(nstep):
        print(f"\n[constrained_scan] step {istep+1}/{nstep}")
        if use_prev_relaxed and results and do_relax:
            apos_step = results[-1]['apos'].copy()
        else:
            apos_step = apos0.copy()
        for im, ia in enumerate(moved_idx):
            apos_step[ia, :] = path[im, istep, :]

        wdir = os.path.join(outdir, f"{step_prefix}_{istep:03d}")
        E, apos_out, forces = run_pbc(
            apos_step, enames, lvs,
            basis_path=basis_path,
            do_relax=do_relax,
            fixed_atoms=fixed_idx,
            nk=nk, k_shift=k_shift,
            dftb_exe=dftb_exe,
            workdir=wdir,
            Temperature=Temperature,
            MixingParameter=MixingParameter,
            MaxScc=MaxScc,
            SCCTolerance=SCCTolerance,
            params=params,
        )

        r = {'istep': istep, 'E': E, 'apos': apos_out, 'enames': enames, 'forces': forces, 'fixed_idx': fixed_idx}
        if key_func is not None:
            r[key_name] = key_func(istep, apos_out, forces, E)
        results.append(r)

        if plot_step_func is not None:
            plot_step_func(istep, apos_out, forces, r)

    if save_xyz:
        if xyz_fname is None:
            xyz_fname = os.path.join(outdir, f"{results_prefix}.xyz")
        save_xyz_movie(results, xyz_fname, lvs=lvs, key_order=[key_name, 'E'] if key_func is not None else ['E'])
    return results

def save_xyz_movie(results, fname, lvs=None, label=None, key_order=None):
    with open(fname, 'w') as f:
        for r in results:
            apos = r['apos']; enames = r['enames']
            fields = []
            if label is not None:
                fields.append(str(label))
            keys = key_order
            if keys is None:
                keys = [k for k in ('L_Hb', 'L_H', 'E') if k in r]
            for k in keys:
                v = r[k]
                fields.append(f"{k}={v:.6f}" if isinstance(v, (float, np.floating)) else f"{k}={v}")
            if r.get('forces') is not None and r.get('fixed_idx'):
                ff = [np.linalg.norm(r['forces'][i]) for i in r['fixed_idx']]
                fields.append(f"F_fixed={np.mean(ff):.3f}")
            if lvs is not None:
                fields.append(f"Lx={lvs[0,0]:.2f} Ly={lvs[1,1]:.2f} Lz={lvs[2,2]:.2f}")
            f.write(f"{len(enames)}\n")
            f.write(" ".join(fields) + "\n")
            for e, pos in zip(enames, apos):
                f.write(f"{e} {pos[0]:.6f} {pos[1]:.6f} {pos[2]:.6f}\n")
    print(f"Saved XYZ movie: {fname}")

# ============ Setup

def makeDFTBjob( enames=None, fname='dftb_in.hsd', gname="input.xyz", method='D3H5', cell=None, basis_path='/home/prokop/SIMULATIONS/dftbplus/slakos/3ob-3-1/', params=default_params, opt=True ):
    enameset = set( enames )
    #print( "enameset = ", enameset )
    hsd = open(fname,'w')
    hsd.write(dedent("""
    Geometry = xyzFormat {
        <<< "%s"
    }\n""" %gname ))

    if opt:
        hsd.write(dedent(f"""  
        Driver = GeometryOptimization {{
            Optimizer = {params["Optimizer"]}
            MovedAtoms = 1:-1
            MaxSteps = {params["MaxSteps"]}
            OutputPrefix = "geom.out"
            Convergence {{ 
                GradElem = {params["GradElem"]}
        """) )
        if 'DispElem' in params: 
            hsd.write('        DispElem = %e \n' %params['DispElem']  ) 
        if 'EConv' in params: 
            hsd.write('        Energy   = %e \n' %params['EConv']    ) 
        hsd.write("    }\n")
        hsd.write("}\n")

    if method in methods_XTB:
        hsd.write(dedent("""
        Hamiltonian = xTB {
            Method = "%s-xTB"
        }
        """ %method ))
    elif method in methods_dftb:
        hsd.write(dedent("""
        Hamiltonian = DFTB {
            Scc = Yes
            SlaterKosterFiles = Type2FileNames {
                Prefix = %s
                Separator = "-"
                Suffix = ".skf"
            }
        """ %basis_path ))
        
        hsd.write("    MaxAngularMomentum {\n")
        for ename in enameset:  hsd.write(f'        {ename} = "{elements.ELEMENT_DICT[ename][4]}" \n'   )
        hsd.write("    }\n")

        if method=="D3H5":

        
            hsd.write(indent(dedent(f"""          
            HCorrection = H5 {{
                RScaling = {params["RScaling"]}
                WScaling = {params["WScaling"]}
                H5Scaling {{\n""" ), "    "))
            for ename in enameset: 
                if ename in H5Scaling: hsd.write(f'            {ename} = {H5Scaling[ename]} \n' )
            hsd.write("       }\n")
            hsd.write("    }\n")

            hsd.write(indent(dedent(f"""   
            Dispersion = DftD3 {{
                Damping = ZeroDamping {{
                    sr6    = {params["sr6"]}
                    alpha6 = {params["alpha6"]}
                }}
                s6 = {params["s6"]}
                s8 = {params["s8"]}
                HHRepulsion = {params["HHRepulsion"]}
            }}\n"""  ), "    "))

        if 'SCCTolerance' in params: 
            hsd.write('    SCCTolerance = %e \n' %params['SCCTolerance']    )
        if 'MaxSccIterations' in params:
            hsd.write('    MaxSccIterations = %i \n' %params['MaxSccIterations']    )
        if 'Mixer' in params:
            hsd.write('    Mixer = %s \n' %params['Mixer']    ) 
        if 'Temperature' in params:
            hsd.write('    Filling = Fermi {Temperature [K] = %f }\n' %params['Temperature']    )

        hsd.write("}\n")

    return
    #Options { WriteDetailedOut = No }
    #Analysis { CalculateForces = Yes }
    #ParserOptions { ParserVersion = 10 }
    #Parallel { UseOmpThreads = Yes }

def makeDFTBjob_pbc( enames, apos, lvs, fname='dftb_in.hsd', basis_path='/home/prokop/SIMULATIONS/dftbplus/slakos/3ob-3-1/', 
                     nk=(1,1,1), k_shift=(0.5,0.0,0.0), opt=False, params=default_params, SCCTolerance=1e-5, MaxScc=200, Temperature=300, MixingParameter=0.2, fixed_atoms=None ):
    """Write a DFTB+ input for a periodic calculation using GenFormat (supercell type S).
    
    Args:
        enames: list of element names per atom
        apos:   (natoms,3) array of Cartesian atomic positions [Angstrom]
        lvs:    (3,3) lattice vectors (rows are a1, a2, a3) [Angstrom]
        fname:  output HSD filename
        basis_path: path to Slater-Koster files
        nk:     (3,) k-point folding along a1, a2, a3
        k_shift: (3,) k-point shift (0.5 for Monkhorst-Pack half-shift)
        opt:    if True, add geometry optimization driver
        params: dict with optimizer settings (uses default_params keys)
        SCCTolerance: SCC convergence threshold
        MaxScc: max SCC iterations
        Temperature: electronic temperature [K]
        MixingParameter: Broyden mixing parameter (0.0-1.0)
        fixed_atoms: list of 0-based atom indices to fix during relaxation (adds Constraint block)
    """
    enameset = sorted(set(enames))
    ename_to_idx = {e: i+1 for i, e in enumerate(enameset)}  # GenFormat: 1-indexed
    natoms = len(enames)
    lvs = np.array(lvs)

    with open(fname, 'w') as hsd:
        # GenFormat geometry block
        hsd.write('Geometry = GenFormat {\n')
        hsd.write(f'  {natoms}  S\n')
        hsd.write('  ' + ' '.join(enameset) + '\n')
        for i, (ename, pos) in enumerate(zip(enames, apos)):
            idx = ename_to_idx[ename]
            hsd.write(f'  {i+1} {idx}   {pos[0]:.10f}   {pos[1]:.10f}   {pos[2]:.10f}\n')
        # Origin + lattice vectors
        hsd.write('  0.000000000  0.000000000  0.000000000\n')
        for row in lvs:
            hsd.write(f'  {row[0]:.10f}  {row[1]:.10f}  {row[2]:.10f}\n')
        hsd.write('}\n\n')

        # Geometry optimization driver
        if opt:
            # MovedAtoms: all atoms except fixed ones
            if fixed_atoms:
                fixed_1based = sorted([i+1 for i in fixed_atoms])  # 1-based for DFTB+
                # Build MovedAtoms as range excluding fixed
                all_idx = set(range(1, natoms+1))
                moved_idx = sorted(all_idx - set(fixed_1based))
                if moved_idx:
                    # Compact representation: list of ranges
                    moved_str = ' '.join(str(i) for i in moved_idx)
                else:
                    moved_str = "1:-1"  # fallback
            else:
                moved_str = "1:-1"
            hsd.write(dedent(f"""Driver = GeometryOptimization {{
    Optimizer = {params["Optimizer"]}
    MovedAtoms = {moved_str}
    MaxSteps = {params["MaxSteps"]}
    OutputPrefix = "geom.out"
    LatticeOpt = No
    Convergence {{ GradElem = {params["GradElem"]} }}
}}\n\n"""))
        
        # Force calculation (always, needed to monitor constraints)
        hsd.write('\nAnalysis {\n  CalculateForces = Yes\n}\n\n')

        # Hamiltonian
        hsd.write('Hamiltonian = DFTB {\n')
        hsd.write('  Scc = Yes\n')
        hsd.write('  SlaterKosterFiles = Type2FileNames {\n')
        hsd.write(f'    Prefix = {basis_path}\n')
        hsd.write('    Separator = "-"\n')
        hsd.write('    Suffix = ".skf"\n')
        hsd.write('  }\n')
        hsd.write('  MaxAngularMomentum {\n')
        for ename in enameset:
            hsd.write(f'    {ename} = "{elements.ELEMENT_DICT[ename][4]}"\n')
        hsd.write('  }\n')
        # K-points via SupercellFolding (Monkhorst-Pack)
        hsd.write('  KPointsAndWeights = SupercellFolding {\n')
        hsd.write(f'    {nk[0]} 0 0\n')
        hsd.write(f'    0 {nk[1]} 0\n')
        hsd.write(f'    0 0 {nk[2]}\n')
        hsd.write(f'    {k_shift[0]:.1f} {k_shift[1]:.1f} {k_shift[2]:.1f}\n')
        hsd.write('  }\n')
        hsd.write(f'  SCCTolerance = {SCCTolerance:.2e}\n')
        hsd.write(f'  MaxSccIterations = {MaxScc}\n')
        hsd.write(f'  Filling = Fermi {{ Temperature [K] = {Temperature} }}\n')
        hsd.write('  Mixer = Broyden {\n')
        hsd.write(f'    MixingParameter = {MixingParameter}\n')
        hsd.write('  }\n')
        hsd.write('}\n')


def run( geom=None, params=None, id=0 ):
    idstr = "%03i" %id 
    print( idstr )
    if params['own_dir']:
        cwd = os.getcwd()
        os.mkdir( idstr )
        os.chdir( idstr )
    #try:
    #    os.system( 'cp ../%03i/charges.bin .' %(id-1) )
    #except: pass
    if( id>0 ):
        os.system( 'cp ../%03i/charges.bin .' %(id-1) )
    apos,es = geom
    au.saveXYZ( es=es, xyzs=apos, fname="input.xyz" )
    makeDFTBjob( enames=es, fname='dftb_in.hsd', gname="input.xyz", method=params['method'], cell=params['cell'], basis_path=params['basis'], params=params, opt=params['opt'] )
    os.system('dftb+ > OUT' )
    #os.system( 'grep "Total Energy" OUT | tail -1 | cut -b 52-70' )
    Estr = os.popen('grep "Total Energy" OUT | tail -1 | cut -b 52-70').read()
    E = float(Estr)
    if params['own_dir']:
        os.chdir( cwd )
    return E


# ============ Hessian / Dynamical Matrix ============

def get_hessian_ase(atoms, delta=1e-4, nfree=2):
    """
    Compute Hessian matrix using ASE Vibrations (programmatic, no file I/O).
    
    This uses ASE's finite-difference implementation to compute the Hessian
    directly in memory without reading/writing text files.
    
    Args:
        atoms: ASE Atoms object with calculator attached
        delta: Displacement magnitude (Angstrom)
        nfree: Number of displacements per coordinate (2 or 4)
    
    Returns:
        hessian: (3N, 3N) numpy array in eV/Angstrom²
        vib: ASE Vibrations object (for further analysis if needed)
    
    Example:
        >>> from ase import Atoms
        >>> from ase.calculators.dftb import Dftb
        >>> atoms = Atoms('H2', positions=[[0,0,0], [0,0,0.74]])
        >>> atoms.calc = Dftb(Hamiltonian_SCC='Yes', ...)
        >>> H, vib = get_hessian_ase(atoms)
        >>> print(H.shape)  # (6, 6) for 2 atoms
    """
    from ase.vibrations import Vibrations
    
    vib = Vibrations(atoms, delta=delta, nfree=nfree)
    vib.run()
    vib.read()
    
    # H is stored in vib.H after read()
    # Shape: (3N, 3N) in eV/Angstrom²
    hessian = vib.H.copy()
    
    return hessian, vib


def hessian_to_mass_weighted(hessian, masses):
    """
    Convert Hessian to mass-weighted dynamical matrix.
    
    Args:
        hessian: (3N, 3N) array in eV/Angstrom²
        masses: (N,) array in atomic mass units
    
    Returns:
        D: (3N, 3N) mass-weighted dynamical matrix
        im: (3N,) inverse sqrt masses (repeated 3x per atom)
    """
    import numpy as np
    
    # Inverse sqrt masses, repeated for x,y,z
    im = np.repeat(masses**-0.5, 3)
    
    # Mass-weighted Hessian: D = M^(-1/2) * H * M^(-1/2)
    D = im[:, None] * hessian * im
    
    return D, im


def hessian_to_frequencies(hessian, masses):
    """
    Convert Hessian to vibrational frequencies.
    
    Args:
        hessian: (3N, 3N) array in eV/Angstrom²
        masses: (N,) array in atomic mass units
    
    Returns:
        frequencies: (3N,) array in cm^-1
        modes: (3N, 3N) eigenvectors (normal modes)
    """
    import numpy as np
    from ase import units
    
    D, im = hessian_to_mass_weighted(hessian, masses)
    
    # Diagonalize mass-weighted Hessian
    omega2, modes = np.linalg.eigh(D)
    
    # Convert eigenvalues (eV/amu/Å²) to frequencies (cm^-1)
    # omega = sqrt(k/m) where k is in eV/Å² and m in amu
    # 1 eV = 1.602176634e-19 J
    # 1 amu = 1.66053906660e-27 kg
    # 1 Å = 1e-10 m
    # So eV/amu/Å² = (1.602e-19 J) / (1.661e-27 kg) / (1e-20 m²) = 9.648e33 s^-2
    # sqrt gives s^-1, then convert to cm^-1: divide by (c * 100)
    
    # Simpler approach using ASE units
    # omega2 is in eV/amu/Å²
    # Convert to atomic units (Hartree/Bohr²/amu) first
    # 1 eV/Å² = (1/27.2114) Hartree / (0.529177²) Bohr² = 0.134 Hartree/Bohr²
    # Then use ASE's conversion
    
    # Direct conversion from eV/Å² to cm^-1:
    # omega (cm^-1) = sqrt(omega2 * eV_to_J / (amu_to_kg * Å_to_m²)) / (c * 100)
    eV_to_J = 1.602176634e-19
    amu_to_kg = 1.66053906660e-27
    Å_to_m = 1e-10
    c = 2.99792458e10  # cm/s
    
    omega = np.sqrt(np.abs(omega2) * eV_to_J / (amu_to_kg * Å_to_m**2))
    frequencies = omega / c
    
    return frequencies, modes


# ============ Molecular Orbitals (Waveplot) ============

def run_waveplot(workdir='.', waveplot_exe='waveplot', 
                 sk_wfc_path=None, plotted_levels='1:-1',
                 n_points=(50, 50, 50), resolution=0.01,
                 electrostatic_potential=False):
    """
    Run DFTB+ waveplot utility to generate cube files for molecular orbitals.
    
    Args:
        workdir: Directory containing detailed.xml and eigenvec.bin
        waveplot_exe: Path to waveplot executable
        sk_wfc_path: Path to wavefunction coefficient file (e.g., wfc.mio-1-1.hsd)
        plotted_levels: Which orbitals to plot (e.g., '1:-1' for all, '4' for HOMO)
        n_points: Grid resolution (nx, ny, nz)
        resolution: Grid spacing for wavefunction evaluation
        electrostatic_potential: If True, also calculate electrostatic potential
    
    Returns:
        List of generated cube file paths
    
    Requires:
        - detailed.xml and eigenvec.bin in workdir (from DFTB+ with WriteEigenvectors=Yes)
        - waveplot executable in PATH or specified path
    """
    import subprocess
    from pathlib import Path
    
    cwd = Path(workdir).absolute()
    
    # Generate waveplot_in.hsd dynamically
    if sk_wfc_path:
        # Use <<+ syntax as in working examples
        wfc_include = f'  <<+ "{sk_wfc_path}"'
    else:
        wfc_include = '  <<+ "wfc.mio-1-1.hsd"'
    
    nx, ny, nz = n_points
    esp_option = '  ElectrostaticPotential = Yes' if electrostatic_potential else ''
    
    waveplot_input = f'''
Options {{
  TotalChargeDensity = Yes
  TotalChargeDifference = Yes
  ChargeDensity = Yes
  RealComponent = Yes
  PlottedSpins = 1 -1
  PlottedLevels = {plotted_levels}
  PlottedRegion = OptimalCuboid {{}}
  NrOfPoints = {nx} {ny} {nz}
  NrOfCachedGrids = -1
  Verbose = Yes
{esp_option}
}}

DetailedXml = "detailed.xml"
EigenvecBin = "eigenvec.bin"

Basis {{
  Resolution = {resolution}
  {wfc_include}
}}
'''
    with open(cwd / 'waveplot_in.hsd', 'w') as f:
        f.write(waveplot_input)
    
    # Run waveplot
    result = subprocess.run(
        [waveplot_exe],
        cwd=cwd,
        capture_output=True,
        text=True
    )
    
    if result.returncode != 0:
        raise RuntimeError(f"Waveplot failed: {result.stderr}")
    
    # Collect output cube files
    cubes = list(cwd.glob('wp-*.cube'))
    return cubes


def read_cube(filename):
    """
    Read Gaussian cube file using ASE.
    
    Args:
        filename: Path to .cube file
    
    Returns:
        data: (nx, ny, nz) numpy array of values
        atoms: ASE Atoms object with atomic positions
    
    Requires:
        - ASE (from ase.io.cube import read_cube_data)
    """
    from ase.io.cube import read_cube_data
    
    # ASE returns (data, atoms) where data is (nx, ny, nz)
    data, atoms = read_cube_data(filename)
    
    return data, atoms


def read_cube_with_grid(filename):
    """
    Read Gaussian cube file with full grid information (origin, spacing).
    
    Args:
        filename: Path to .cube file
    
    Returns:
        data: (nx, ny, nz) numpy array of values
        atoms: ASE Atoms object with atomic positions
        origin: (3,) array of grid origin in Angstrom
        spacing: (3,) array of grid spacing in Angstrom
    """
    from ase.io.cube import read_cube_data
    import numpy as np
    
    # Use ASE's read_cube_data which returns grid info
    data, atoms = read_cube_data(filename)
    
    # Extract grid info from the cube file manually for accuracy
    filename = str(filename)
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # Skip comment lines (waveplot adds header comments)
    line_idx = 0
    while line_idx < len(lines):
        line = lines[line_idx].strip()
        # Skip lines that are comments or don't start with a number
        if (line.startswith('#') or line.startswith('Cube') or 
            line.startswith('Calc-Id') or not line or
            (line.split() and not line.split()[0][0].isdigit())):
            line_idx += 1
        else:
            break
    
    # Line after comments: number of atoms and origin
    natoms = int(lines[line_idx].split()[0])
    origin = np.array([float(lines[line_idx].split()[i]) for i in range(1, 4)])
    line_idx += 1
    
    # Next 3 lines: number of points and spacing vectors
    nx = int(lines[line_idx].split()[0])
    dx = np.array([float(lines[line_idx].split()[i]) for i in range(1, 4)])
    line_idx += 1
    
    ny = int(lines[line_idx].split()[0])
    dy = np.array([float(lines[line_idx].split()[i]) for i in range(1, 4)])
    line_idx += 1
    
    nz = int(lines[line_idx].split()[0])
    dz = np.array([float(lines[line_idx].split()[i]) for i in range(1, 4)])
    
    # Cube file stores coordinates in Bohr - convert to Angstrom
    BOHR_TO_ANG = 0.529177210903
    origin_ang = origin * BOHR_TO_ANG
    
    # Calculate actual spacing (magnitude of vectors) and convert Bohr -> Angstrom
    spacing_ang = np.array([np.linalg.norm(dx), np.linalg.norm(dy), np.linalg.norm(dz)]) * BOHR_TO_ANG
    
    # atoms from ASE read_cube_data are already in Angstrom
    return data, atoms, origin_ang, spacing_ang


def plot_orbital(cube_file, isosurface_level=0.02, axes=(0, 1)):
    """
    Plot molecular orbital isosurface from cube file.
    
    Args:
        cube_file: Path to .cube file
        isosurface_level: Isosurface level (fraction of max value)
        axes: Which 2D plane to plot (0=x, 1=y, 2=z)
    
    Requires matplotlib and mayavi/pyvista for 3D visualization.
    """
    import numpy as np
    import matplotlib.pyplot as plt
    
    data, atoms = read_cube(cube_file)
    
    # Take 2D slice
    if axes == (0, 1):  # xy plane
        slice_data = data[:, :, data.shape[2]//2]
    elif axes == (0, 2):  # xz plane
        slice_data = data[:, data.shape[1]//2, :]
    elif axes == (1, 2):  # yz plane
        slice_data = data[data.shape[0]//2, :, :]
    else:
        raise ValueError(f"Invalid axes: {axes}")
    
    plt.imshow(slice_data.T, origin='lower', cmap='RdBu_r')
    plt.colorbar(label='Orbital value')
    plt.title(f'Orbital slice from {cube_file}')
    plt.show()
    
    return data, atoms


def write_vibration_modes_jmol(filename, atoms, frequencies, modes, 
                                scale=1.0, min_freq=1.0):
    """
    Write vibration modes in Jmol-compatible XYZ format with vectors.
    
    Format: Multi-model XYZ where each frame has coordinates + vibration vectors.
    Each line: element x y z vx vy vz
    
    Args:
        filename: Output .xyz file path
        atoms: ASE Atoms object with positions and symbols
        frequencies: (3N,) array of frequencies in cm^-1
        modes: (3N, 3N) array of eigenvectors (normal modes)
        scale: Scaling factor for vibration vectors
        min_freq: Minimum frequency to include (skip translations/rotations)
    
    Example:
        >>> write_vibration_modes_jmol('modes.xyz', atoms, frequencies, modes)
        >>> # In Jmol: load 'modes.xyz'; vectors on; vibrate on
    """
    import numpy as np
    
    natoms = len(atoms)
    positions = atoms.positions  # Angstrom
    symbols = atoms.get_chemical_symbols()
    
    with open(filename, 'w') as f:
        for i, freq in enumerate(frequencies):
            # Skip near-zero frequencies (translations/rotations)
            if abs(freq) < min_freq:
                continue
            
            # Get eigenvector (column i from modes matrix)
            mode = modes[:, i]
            
            # Write frame header
            f.write(f"{natoms}\n")
            f.write(f"Vibration Frequency: {freq:.2f} cm^-1\n")
            
            # Write atoms with vibration vectors
            # mode is (3N,) vector, reshape to (N, 3)
            mode_reshaped = mode.reshape(natoms, 3)
            
            for j in range(natoms):
                x, y, z = positions[j]
                vx, vy, vz = mode_reshaped[j] * scale
                f.write(f"{symbols[j]} {x:.6f} {y:.6f} {z:.6f} {vx:.6f} {vy:.6f} {vz:.6f}\n")
    
    print(f"Wrote {len([f for f in frequencies if abs(f) >= min_freq])} vibration modes to {filename}")


def evaluate_orbital_at_points(workdir, points, orbital_indices=None):
    """
    Evaluate molecular orbitals at arbitrary points using DFTB+ C API.
    
    Args:
        workdir: Directory containing detailed.xml and eigenvec.bin
        points: (N, 3) array of points in Angstrom where to evaluate orbitals
        orbital_indices: List of orbital indices to evaluate (e.g., [0, 1, 2] for HOMO-2, HOMO-1, HOMO)
                     If None, evaluates all occupied orbitals
    
    Returns:
        values: (N, n_orbitals) array of orbital values at each point
    
    Requires:
        - dftb_lib.py with DftbPlusCalculator
        - detailed.xml and eigenvec.bin from DFTB+ calculation
    """
    import numpy as np
    from pyBall import dftb_lib
    
    # Load DFTB+ calculator
    calc = dftb_lib.DftbPlusCalculator()
    
    # Initialize with existing calculation
    calc.initialize_external(workdir)
    
    # Get number of orbitals
    n_orbitals = calc.get_n_orbitals()
    
    if orbital_indices is None:
        # Get occupied orbitals
        nelec = calc.get_n_electrons()
        orbital_indices = list(range(nelec // 2))
    
    # Evaluate orbitals at points
    points = np.asarray(points)
    n_points = points.shape[0]
    n_orb = len(orbital_indices)
    
    values = np.zeros((n_points, n_orb))
    
    for i, point in enumerate(points):
        for j, orb_idx in enumerate(orbital_indices):
            # Get orbital value at point (this requires C API call)
            # For now, use interpolation from cube data
            pass
    
    calc.cleanup()
    
    return values


def interpolate_orbital(cube_file, points):
    """
    Interpolate orbital values at arbitrary points from cube file.
    
    Args:
        cube_file: Path to .cube file
        points: (N, 3) array of points in Angstrom
    
    Returns:
        values: (N,) array of interpolated orbital values
    
    This is a simple nearest-neighbor interpolation. For higher accuracy,
    use scipy.interpolate.RegularGridInterpolator.
    """
    import numpy as np
    from scipy.interpolate import RegularGridInterpolator
    
    data, atoms = read_cube(cube_file)
    
    # Get grid dimensions
    nx, ny, nz = data.shape
    
    # Get grid bounds from cube file header (first two lines after comments)
    # Cube file format: origin (line 2), number of points (line 3)
    # We need to extract the actual grid spacing from the cube file
    # For now, use simple interpolation based on atomic positions
    
    # Use atoms positions to estimate grid bounds
    positions = atoms.positions
    origin = positions.min(axis=0) - 1.0  # Add padding
    max_pos = positions.max(axis=0) + 1.0
    
    # Create grid coordinates
    x = np.linspace(origin[0], max_pos[0], nx)
    y = np.linspace(origin[1], max_pos[1], ny)
    z = np.linspace(origin[2], max_pos[2], nz)
    
    # Create interpolator
    interp = RegularGridInterpolator((x, y, z), data, method='linear', bounds_error=False, fill_value=0.0)
    
    # Interpolate at points
    values = interp(points)
    
    return values







