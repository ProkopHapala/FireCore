
import sys
import os
import numpy as np
from . import elements
from . import atomicUtils as au
from textwrap import dedent,indent

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







