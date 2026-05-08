#!/usr/bin/env python3
"""
DFTBplusParser - Parse DFTB+ output files for waveplot OpenCL reimplementation.

Handles:
- detailed.xml: Geometry, occupations, k-points
- eigenvec.bin: Binary eigenvector data (Fortran unformatted)
- Basis HSD: Slater-type orbital parameters
"""

import numpy as np
import struct
import xml.etree.ElementTree as ET
from pathlib import Path


class DFTBplusParser:
    """Parser for DFTB+ waveplot input files."""
    
    def __init__(self, work_dir='.', verbosity=0):
        self.work_dir = Path(work_dir)
        self.verbosity = verbosity
        self.geometry = None
        self.eigenvectors = None
        self.basis = None
        self.kpoints = None
        self.occupations = None
        self.t_real = True
        self.identity = None
        
    def parse_detailed_xml(self, xml_path=None):
        """Parse detailed.xml for geometry, occupations, k-points."""
        if xml_path is None:
            xml_path = self.work_dir / 'detailed.xml'
        
        tree = ET.parse(xml_path)
        root = tree.getroot()
        
        # Get identity
        self.identity = int(root.find('Identity').text)
        
        # Parse geometry
        geo = root.find('Geometry')
        atoms = []
        species_names = []
        species_dict = {}
        
        for atom in geo.findall('Atom'):
            coords = [float(x) for x in atom.get('coords').split()]
            species = atom.get('species')
            if species not in species_dict:
                species_dict[species] = len(species_names)
                species_names.append(species)
            atoms.append({
                'coords': np.array(coords),
                'species': species_dict[species],
                'species_name': species
            })
        
        natoms = len(atoms)
        coords = np.array([a['coords'] for a in atoms])
        species = np.array([a['species'] for a in atoms], dtype=np.int32)
        
        # Check for periodic
        t_periodic = False
        lat_vecs = None
        periodic = root.find('Periodic')
        if periodic is not None:
            t_periodic = True
            lat_vecs = np.zeros((3, 3))
            for i, vec in enumerate(periodic.findall('LatVec')):
                lat_vecs[i] = [float(x) for x in vec.text.split()]
        
        self.geometry = {
            'natoms': natoms,
            'coords': coords,
            'species': species,
            'species_names': species_names,
            't_periodic': t_periodic,
            'lat_vecs': lat_vecs
        }
        
        # Parse k-points
        kp_node = root.find('KPoints')
        if kp_node is not None:
            nkp = int(kp_node.get('number'))
            kpoints = np.zeros((3, nkp))
            for i, kp in enumerate(kp_node.findall('KPoint')):
                kpoints[:, i] = [float(x) for x in kp.text.split()]
            self.kpoints = kpoints
        else:
            self.kpoints = np.zeros((3, 1))
            
        # Parse occupations
        occ_node = root.find('Occupations')
        if occ_node is not None:
            nspin = int(occ_node.get('nSpin'))
            nstate = int(occ_node.get('nStates'))
            nkpt = int(occ_node.get('nKPoints'))
            self.occupations = np.zeros((nstate, nkpt, nspin))
            for spin in range(nspin):
                for ik in range(nkpt):
                    occ_str = occ_node.get(f'spin{spin+1}_k{ik+1}')
                    if occ_str:
                        self.occupations[:, ik, spin] = [float(x) for x in occ_str.split()]
        
        # Real or complex
        real_node = root.find('Real')
        if real_node is not None:
            self.t_real = real_node.text.lower() == 'true'
            
        if self.verbosity > 0:
            print(f"[DFTBplusParser] Parsed detailed.xml:")
            print(f"  Identity: {self.identity}")
            print(f"  Atoms: {natoms}, Species: {len(species_names)}")
            print(f"  K-points: {self.kpoints.shape[1]}")
            print(f"  Real: {self.t_real}")
            
        return self.geometry
    
    def parse_eigenvec_bin(self, bin_path=None, nstates=None, nkpoints=None, nspin=None, norb=None):
        """Parse binary eigenvector file (Fortran unformatted)."""
        if bin_path is None:
            bin_path = self.work_dir / 'eigenvec.bin'
        
        with open(bin_path, 'rb') as f:
            # Read identity (Fortran record marker + int + record marker)
            rec_start = struct.unpack('i', f.read(4))[0]
            identity = struct.unpack('i', f.read(4))[0]
            rec_end = struct.unpack('i', f.read(4))[0]
            
            if self.verbosity > 0:
                print(f"[DFTBplusParser] eigenvec.bin identity: {identity}")
            
            # Determine dimensions from detailed.xml if available
            if nstates is None and self.occupations is not None:
                nstates = self.occupations.shape[0]
                nkpoints = self.occupations.shape[1]
                nspin = self.occupations.shape[2]
            
            if norb is None and self.basis is not None:
                norb = sum(self.basis['norb_per_atom'])
            
            if any(x is None for x in [nstates, nkpoints, nspin, norb]):
                raise ValueError("Must provide nstates, nkpoints, nspin, norb or parse detailed.xml and basis first")
            
            # Read eigenvectors
            if self.t_real:
                eigenvecs = np.zeros((norb, nstates, nkpoints, nspin), dtype=np.float64)
            else:
                eigenvecs = np.zeros((norb, nstates, nkpoints, nspin), dtype=np.complex128)
            
            for ispin in range(nspin):
                for ik in range(nkpoints):
                    for istate in range(nstates):
                        # Read record markers and data
                        rec_start = struct.unpack('i', f.read(4))[0]
                        if self.t_real:
                            data = np.fromfile(f, dtype=np.float64, count=norb)
                        else:
                            data = np.fromfile(f, dtype=np.float64, count=2*norb)
                            data = data[0::2] + 1j * data[1::2]
                        rec_end = struct.unpack('i', f.read(4))[0]
                        eigenvecs[:, istate, ik, ispin] = data
            
        self.eigenvectors = eigenvecs
        
        if self.verbosity > 0:
            print(f"[DFTBplusParser] Parsed eigenvec.bin: shape {eigenvecs.shape}")
            
        return eigenvecs
    
    def parse_basis_hsd(self, hsd_path=None):
        """Parse basis definition from HSD format."""
        if hsd_path is None:
            # Try to find in work_dir
            candidates = list(self.work_dir.glob('*.hsd')) + list(self.work_dir.glob('wfc.*.hsd'))
            if candidates:
                hsd_path = candidates[0]
        
        if hsd_path is None:
            raise FileNotFoundError("Could not find basis HSD file")
            
        hsd_path = Path(hsd_path)
        
        # Simple HSD parser (not full implementation, just basis extraction)
        species_basis = []
        current_species = None
        current_orbitals = []
        
        with open(hsd_path, 'r') as f:
            content = f.read()
            
        # Parse basis block
        import re
        
        # Find Basis block
        basis_match = re.search(r'Basis\s*{([^}]*)}', content, re.DOTALL)
        if not basis_match:
            raise ValueError("No Basis block found in HSD")
        
        basis_content = basis_match.group(1)
        
        # Parse species entries
        species_pattern = r'(\w+)\s*{([^}]*)}'
        for match in re.finditer(species_pattern, basis_content):
            species_name = match.group(1)
            species_block = match.group(2)
            
            # Parse atomic number
            an_match = re.search(r'AtomicNumber\s*=\s*(\d+)', species_block)
            atomic_number = int(an_match.group(1)) if an_match else 0
            
            # Parse resolution
            res_match = re.search(r'Resolution\s*=\s*(\S+)', species_block)
            resolution = float(res_match.group(1)) if res_match else 0.02
            
            # Parse orbitals
            orbitals = []
            orb_pattern = r'Orbital\s*{([^}]*)}'
            for orb_match in re.finditer(orb_pattern, species_block):
                orb_block = orb_match.group(1)
                
                # Parse angular momentum
                l_match = re.search(r'AngularMomentum\s*=\s*(\d+)', orb_block)
                l = int(l_match.group(1)) if l_match else 0
                
                # Parse occupation
                occ_match = re.search(r'Occupation\s*=\s*(\S+)', orb_block)
                occ = float(occ_match.group(1)) if occ_match else 0.0
                
                # Parse cutoff
                cut_match = re.search(r'Cutoff\s*=\s*(\S+)', orb_block)
                cutoff = float(cut_match.group(1)) if cut_match else 10.0
                
                # Parse exponents
                exp_match = re.search(r'Exponents\s*=\s*{([^}]*)}', orb_block)
                if exp_match:
                    exponents = [float(x) for x in exp_match.group(1).split()]
                else:
                    exponents = []
                
                # Parse coefficients
                coef_match = re.search(r'Coefficients\s*=\s*{([^}]*)}', orb_block)
                if coef_match:
                    all_coefs = [float(x) for x in coef_match.group(1).split()]
                    nexp = len(exponents)
                    npow = len(all_coefs) // nexp if nexp > 0 else 0
                    coefficients = np.array(all_coefs).reshape(nexp, npow).T  # (nPow, nAlpha)
                else:
                    coefficients = np.array([])
                
                orbitals.append({
                    'l': l,
                    'occupation': occ,
                    'cutoff': cutoff,
                    'exponents': np.array(exponents),
                    'coefficients': coefficients,  # (nPow, nAlpha)
                    'nAlpha': len(exponents),
                    'nPow': coefficients.shape[0] if coefficients.size > 0 else 0
                })
            
            species_basis.append({
                'name': species_name,
                'atomic_number': atomic_number,
                'resolution': resolution,
                'orbitals': orbitals
            })
        
        self.basis = {
            'species': species_basis,
            'nspecies': len(species_basis)
        }
        
        # Compute orbital layout per atom
        if self.geometry is not None:
            norb_per_atom = []
            for ia in range(self.geometry['natoms']):
                ispec = self.geometry['species'][ia]
                orb_count = sum(2*orb['l'] + 1 for orb in species_basis[ispec]['orbitals'])
                norb_per_atom.append(orb_count)
            
            orb_offsets = np.zeros(self.geometry['natoms'] + 1, dtype=np.int32)
            orb_offsets[1:] = np.cumsum(norb_per_atom)
            
            self.basis['norb_per_atom'] = np.array(norb_per_atom, dtype=np.int32)
            self.basis['orb_offsets'] = orb_offsets
            self.basis['total_norb'] = orb_offsets[-1]
        
        if self.verbosity > 0:
            print(f"[DFTBplusParser] Parsed basis: {len(species_basis)} species")
            for sp in species_basis:
                print(f"  {sp['name']} (Z={sp['atomic_number']}): {len(sp['orbitals'])} orbitals")
        
        return self.basis
    
    def load_all(self, work_dir=None):
        """Load all DFTB+ output files."""
        if work_dir is not None:
            self.work_dir = Path(work_dir)
        
        self.parse_detailed_xml()
        self.parse_basis_hsd()
        self.parse_eigenvec_bin()
        
        return {
            'geometry': self.geometry,
            'basis': self.basis,
            'eigenvectors': self.eigenvectors,
            'kpoints': self.kpoints,
            'occupations': self.occupations,
            't_real': self.t_real
        }


BOHR2ANG = 0.5291772109  # Bohr -> Angstrom


def parse_skf_xml(skf_path):
    """
    Parse SK file XML documentation to extract multi-zeta basis parameters.

    SK files contain embedded XML documentation with the actual basis parameters
    used in DFTB+ calculations (multi-zeta STO basis sets).

    Args:
        skf_path: Path to SK file (e.g., H-O.skf)

    Returns:
        dict with basis parameters for each atom:
        {
            'atom1': {
                'shells': list of shell types (e.g., ['1s']),
                'exponents': list of exponents (multi-zeta),
                'power': power parameter,
                'wavefunction': wavefunction parameter
            },
            'atom2': { ... }
        }
    """
    import xml.etree.ElementTree as ET
    import re

    with open(skf_path, 'r') as f:
        content = f.read()

    # Find XML documentation block
    match = re.search(r'<Documentation>(.*?)</Documentation>', content, re.DOTALL)
    if not match:
        raise ValueError("No XML documentation found in SK file")

    xml_content = match.group(1)
    root = ET.fromstring("<Documentation>" + xml_content + "</Documentation>")

    basis_info = {}

    # Extract basis information for each atom
    for basis in root.findall('.//Basis'):
        atom_num = basis.get('atom')
        shells = basis.find('Shells').text.strip().split()
        exponents = [float(x) for x in basis.find('Exponents').text.split()]
        power = int(basis.find('Power').text)
        wavefunction = [float(x) for x in basis.find('Wavefunction').text.split()]

        basis_info[f'atom{atom_num}'] = {
            'shells': shells,
            'exponents': exponents,
            'power': power,
            'wavefunction': wavefunction
        }

    return basis_info


def parse_wfc_hsd(wfc_path):
    """
    Parse wfc.*.hsd file to extract multi-zeta STO basis parameters.

    The wfc.*.hsd files are in atomic units (Bohr). This function keeps
    them in Bohr units for STO evaluation, and the evaluate_sto_1d function
    should be called with r in Bohr.

    Args:
        wfc_path: Path to wfc.*.hsd file

    Returns:
        dict with structure:
        {
            'H': {
                'AtomicNumber': 1,
                'orbitals': [
                    {
                        'AngularMomentum': 0,
                        'Occupation': 1.0,
                        'Cutoff': 5.0,  # in Bohr
                        'Exponents': [0.5, 1.0, 2.0],  # in Bohr^-1
                        'Coefficients': [[...], [...], [...]]  # nPow x nAlpha
                    }
                ]
            },
            ...
        }
    """
    import re

    with open(wfc_path, 'r') as f:
        content = f.read()

    # Parse using brace counting to handle nested structures
    basis_data = {}
    i = 0
    n = len(content)

    while i < n:
        # Skip whitespace
        while i < n and content[i].isspace():
            i += 1

        if i >= n:
            break

        # Find species name (before '=')
        if content[i] == '}':
            i += 1
            continue

        # Extract species name (until we hit '=', '{', or whitespace)
        species_start = i
        while i < n and content[i] not in ('=', '{') and not content[i].isspace():
            i += 1
        species_name = content[species_start:i].strip()

        # Skip whitespace
        while i < n and content[i].isspace():
            i += 1

        # Skip optional '=' (some wfc files use 'Species = {', others use 'Species {')
        if i < n and content[i] == '=':
            i += 1
            # Skip whitespace after '='
            while i < n and content[i].isspace():
                i += 1

        # Expect '{'
        if i >= n or content[i] != '{':
            print(f"Warning: Expected '{{' after species {species_name}")
            continue
        i += 1

        # Extract species block (count braces)
        brace_count = 1
        block_start = i
        while i < n and brace_count > 0:
            if content[i] == '{':
                brace_count += 1
            elif content[i] == '}':
                brace_count -= 1
            i += 1

        species_block = content[block_start:i-1]

        # Parse species block
        atomic_num = None
        orbitals = []

        # Extract atomic number
        atomic_match = re.search(r'AtomicNumber\s*=\s*(\d+)', species_block)
        if atomic_match:
            atomic_num = int(atomic_match.group(1))

        # Extract orbital blocks using brace counting
        orbital_idx = 0
        while orbital_idx < len(species_block):
            # Find next Orbital
            orbital_match = species_block.find('Orbital', orbital_idx)
            if orbital_match == -1:
                break

            # Find '{' after Orbital (handle both 'Orbital = {' and 'Orbital {' formats)
            brace_start = species_block.find('{', orbital_match)
            if brace_start == -1:
                break

            # Extract orbital block
            brace_count = 1
            block_start = brace_start + 1
            j = block_start
            while j < len(species_block) and brace_count > 0:
                if species_block[j] == '{':
                    brace_count += 1
                elif species_block[j] == '}':
                    brace_count -= 1
                j += 1

            orbital_block = species_block[block_start:j-1]
            orbital_idx = j

            # Parse orbital block
            l_match = re.search(r'AngularMomentum\s*=\s*(\d+)', orbital_block)
            l = int(l_match.group(1)) if l_match else 0

            occ_match = re.search(r'Occupation\s*=\s*([\d.]+)', orbital_block)
            occ = float(occ_match.group(1)) if occ_match else 0.0

            cutoff_match = re.search(r'Cutoff\s*=\s*([\d.]+)', orbital_block)
            cutoff = float(cutoff_match.group(1)) if cutoff_match else 5.0

            # Extract exponents (handle both 'Exponents = {' and 'Exponents {' formats)
            exp_match = re.search(r'Exponents\s*=?\s*\{([^}]+)\}', orbital_block, re.DOTALL)
            if exp_match:
                exps = [float(x) for x in exp_match.group(1).split()]
            else:
                exps = []

            # Extract coefficients (handle both 'Coefficients = {' and 'Coefficients {' formats)
            coeff_match = re.search(r'Coefficients\s*=?\s*\{([^}]+)\}', orbital_block, re.DOTALL)
            if coeff_match:
                coeff_lines = [line.strip() for line in coeff_match.group(1).split('\n') if line.strip()]
                coeffs = []
                for line in coeff_lines:
                    coeffs.extend([float(x) for x in line.split()])
                # Fortran reads: reshape(coeffs, [nPow, nAlpha]) with column-major order
                # i.e. fills aa(1,1),aa(2,1),...,aa(nPow,1),aa(1,2),... 
                n_alpha = len(exps)
                n_pow = len(coeffs) // n_alpha
                coeffs = np.array(coeffs).reshape(n_pow, n_alpha, order='F')  # (nPow, nAlpha)
            else:
                coeffs = np.array([[1.0]])

            orbitals.append({
                'AngularMomentum': l,
                'Occupation': occ,
                'Cutoff': cutoff,
                'Exponents': exps,
                'Coefficients': coeffs
            })

        basis_data[species_name] = {
            'AtomicNumber': atomic_num,
            'orbitals': orbitals
        }

    return basis_data


def evaluate_sto_1d(r, l, exps, coeffs):
    """
    Evaluate multi-zeta STO function at distance r using numpy vectorization.

    STO(r) = sum_{i=1}^{nAlpha} exp(-exps[i] * r) * sum_{j=1}^{nPow} coeffs[j,i] * r^{l + j - 1}

    This matches the waveplot implementation in slater.F90.
    All units are in Bohr (atomic units).

    Args:
        r: distance array in Bohr (can be scalar or array)
        l: angular momentum
        exps: array of exponents in Bohr^-1 (nAlpha,)
        coeffs: coefficient matrix (nPow, nAlpha)

    Returns:
        STO values at r (same shape as r)
    """
    r = np.asarray(r, dtype=float)
    n_alpha = len(exps)
    n_pow = coeffs.shape[0]  # coeffs shape: (nPow, nAlpha)

    # Build r^l, r^(l+1), ..., r^(l+nPow-1) = pows array
    # Matching Fortran: pows(1)=r^ll, pows(ii+1)=pows(ii)*rr
    # For l=0 at r=0: rTmp=1 (Fortran avoids 0^0)
    r0 = np.where((l == 0) & (r < 1e-12), 1.0, r ** l)
    r_powers = r0[:, None] * r[:, None] ** np.arange(n_pow)[None, :]  # (N, nPow)
    # Fix: first column should be r^l, not r^(l+1)
    # Actually: pows[0]=r^l, pows[1]=r^l*r=r^(l+1), etc.
    # r0 = r^l already, so: r_powers[:,j] = r0 * r^j = r^(l+j)
    # But that gives r^l * r^0=r^l for j=0: correct!

    # Compute exponentials: exp(-exps[i] * r)
    exp_terms = np.exp(-np.outer(r, exps))  # (N, nAlpha)

    # sum_j coeffs[j,i] * r_powers[:,j] for each i -> (N, nAlpha)
    radial_sum = r_powers @ coeffs   # (N, nPow) @ (nPow, nAlpha) = (N, nAlpha)

    # sum over i: sto[n] = sum_i exp_terms[n,i] * radial_sum[n,i]
    sto = (exp_terms * radial_sum).sum(axis=1)  # (N,)

    return sto


def evaluate_sto_2d(x_grid, y_grid, l, exps, coeffs, origin=(0, 0)):
    """
    Evaluate multi-zeta STO function on 2D grid (s-orbital radial part).

    For s-orbitals (l=0), this is just the radial function evaluated at distance from origin.

    Args:
        x_grid, y_grid: 2D meshgrid arrays
        l: angular momentum
        exps: array of exponents (nAlpha,)
        coeffs: coefficient matrix (nPow, nAlpha)
        origin: (x0, y0) center of the orbital

    Returns:
        STO values on 2D grid
    """
    r = np.sqrt((x_grid - origin[0])**2 + (y_grid - origin[1])**2)
    return evaluate_sto_1d(r, l, exps, coeffs)


def parse_basis_hsd_ang(hsd_path):
    """
    Parse waveplot_in.hsd Basis block and return species_list ready for load_basis_sto().

    All parameters in the HSD are in Bohr (DFTB+ internal units). This function
    converts them to Angstrom for use with OCL kernels that work in Angstrom:
      alpha_Ang  = alpha_Bohr  / BOHR2ANG
      cutoff_Ang = cutoff_Bohr * BOHR2ANG
      res_Ang    = res_Bohr    * BOHR2ANG
      coeff_Ang  = coeff_Bohr  / BOHR2ANG^l   (r^l carries units)

    Handles <<+ "filename.hsd" include directive to use official wfc files.

    Returns list of species dicts compatible with GridProjector.load_basis_sto().
    """
    import re
    B = BOHR2ANG
    with open(hsd_path, 'r') as f:
        content = f.read()

    # Check for <<+ include directive (use official wfc file)
    include_match = re.search(r'<<\+\s*"([^"]+)"', content)
    if include_match:
        wfc_path = Path(hsd_path).parent / include_match.group(1)
        if wfc_path.exists():
            print(f"[parse_basis_hsd_ang] Using included wfc file: {wfc_path}")
            basis_data = parse_wfc_hsd(str(wfc_path))
            # Convert from Bohr to Angstrom for OCL kernels
            species_list = []
            for sp_name, sp_data in basis_data.items():
                orbitals = []
                for orb in sp_data['orbitals']:
                    l = orb['AngularMomentum']
                    exps_b = np.asarray(orb['Exponents'], dtype=np.float64)  # in Bohr^-1
                    coeffs_b = np.asarray(orb['Coefficients'], dtype=np.float64)  # (nPow, nAlpha) in Bohr
                    cutoff_b = orb['Cutoff']  # in Bohr
                    nPow = coeffs_b.shape[0]
                    # Each power term j contributes r^(l+j), so scale by B^(l+j)
                    scale_factors = np.array([B ** (l + j) for j in range(nPow)])
                    coeffs_scaled = coeffs_b / scale_factors[:, None]  # (nPow, 1) broadcasting
                    orbitals.append({
                        'l': l,
                        'cutoff': cutoff_b * B,
                        'exponents': exps_b / B,
                        'coefficients': coeffs_scaled,
                    })
                species_list.append({
                    'name': sp_name,
                    'atomic_number': sp_data['AtomicNumber'],
                    'orbitals': orbitals,
                    'resolution': 0.04 * B,  # Default resolution, convert to Angstrom
                })
            # Parse top-level Resolution from original HSD
            res_match = re.search(r'Basis\s*\{[^}]*Resolution\s*=\s*(\S+)', content, re.DOTALL)
            res_bohr = float(res_match.group(1)) if res_match else 0.04
            for sp in species_list:
                sp['resolution'] = res_bohr * B
            return species_list
        else:
            print(f"[WARNING] wfc file not found: {wfc_path}, falling back to inline HSD parsing")
    else:
        print(f"[parse_basis_hsd_ang] No <<+ include directive found, using inline HSD parsing")

    # Parse top-level Resolution
    res_match = re.search(r'Basis\s*\{[^}]*Resolution\s*=\s*(\S+)', content, re.DOTALL)
    res_bohr = float(res_match.group(1)) if res_match else 0.04

    # Find Basis block (handles nested braces via simple depth tracking)
    start = content.find('Basis')
    depth = 0; basis_content = ''
    for i, ch in enumerate(content[start:], start):
        if ch == '{': depth += 1
        elif ch == '}':
            depth -= 1
            if depth == 0: basis_content = content[start:i+1]; break

    species_list = []
    # Match species blocks: 'Name = {' or 'Name{' — skip HSD keywords
    _SKIP = {'Basis', 'Orbital', 'Options', 'PlottedRegion', 'Box', 'Origin',
             'Resolution', 'AtomicNumber', 'AngularMomentum', 'Occupation',
             'Cutoff', 'Exponents', 'Coefficients', 'NrOfPoints', 'RealComponent',
             'Verbose', 'PlottedLevels', 'PlottedKPoints', 'PlottedSpins', 'GroundState'}
    sp_pat = re.compile(r'\b([A-Z][a-zA-Z0-9]*)\s*=\s*\{')  # species: capital letter, then = {
    i = basis_content.find('{') + 1  # skip opening Basis {
    while i < len(basis_content) - 1:
        m = sp_pat.search(basis_content, i)
        if not m: break
        sp_name = m.group(1)
        if sp_name in _SKIP: i = m.end(); continue
        # Extract species block (the { starting at m.end()-1)
        b_start = m.end() - 1
        depth = 0; b_end = b_start
        for j, ch in enumerate(basis_content[b_start:], b_start):
            if ch == '{': depth += 1
            elif ch == '}':
                depth -= 1
                if depth == 0: b_end = j; break
        sp_block = basis_content[b_start:b_end+1]
        i = b_end + 1

        an_m = re.search(r'AtomicNumber\s*=\s*(\d+)', sp_block)
        atomic_number = int(an_m.group(1)) if an_m else 0

        orbitals = []
        orb_pat = re.compile(r'Orbital\s*=\s*\{')
        oi = 0
        while True:
            om = orb_pat.search(sp_block, oi)
            if not om: break
            od = 0; oe = om.start()
            for j, ch in enumerate(sp_block[om.start():], om.start()):
                if ch == '{': od += 1
                elif ch == '}':
                    od -= 1
                    if od == 0: oe = j; break
            orb_block = sp_block[om.start():oe+1]
            oi = oe + 1

            l_m   = re.search(r'AngularMomentum\s*=\s*(\d+)', orb_block)
            cut_m = re.search(r'Cutoff\s*=\s*(\S+)', orb_block)
            exp_m = re.search(r'Exponents\s*=\s*\{([^}]*)\}', orb_block)
            cof_m = re.search(r'Coefficients\s*=\s*\{([^}]*)\}', orb_block)

            l      = int(l_m.group(1)) if l_m else 0
            cutoff_b = float(cut_m.group(1)) if cut_m else 10.0
            exps_b = np.array([float(x) for x in exp_m.group(1).split()]) if exp_m else np.array([1.0])
            if cof_m:
                all_c = np.array([float(x) for x in cof_m.group(1).split()])
                nexp  = len(exps_b)
                npow  = len(all_c) // nexp
                coef_b = all_c.reshape(nexp, npow).T  # (nPow, nAlpha)
            else:
                coef_b = np.ones((1, len(exps_b)))

            # Convert Bohr -> Angstrom with proper power term scaling
            nPow = coef_b.shape[0]
            scale_factors = np.array([B ** (l + j) for j in range(nPow)])
            coef_scaled = coef_b / scale_factors[:, None]
            orbitals.append({
                'l':            l,
                'cutoff':       cutoff_b * B,
                'exponents':    exps_b   / B,
                'coefficients': coef_scaled,
            })

        species_list.append({
            'name':          sp_name,
            'atomic_number': atomic_number,
            'resolution':    res_bohr * B,
            'orbitals':      orbitals,
        })

    return species_list


def compute_sto_radial(r, aa, alpha, ll):
    """
    Compute Slater-type orbital radial part analytically.
    
    Formula: R_l(r) = sum_{i=1}^{nAlpha} [ sum_{j=1}^{nPow} aa(j,i) * r^{l+j-1} ] * exp(-alpha_i * r)
    
    Args:
        r: radial distance (scalar or array)
        aa: (nPow, nAlpha) coefficients
        alpha: (nAlpha,) exponents
        ll: angular momentum
    
    Returns:
        STO radial value
    """
    aa = np.asarray(aa)
    alpha = np.asarray(alpha)
    r = np.asarray(r)
    
    nAlpha = len(alpha)
    nPow = aa.shape[0] if aa.ndim > 1 else 1
    
    result = np.zeros_like(r, dtype=np.float64)
    
    for i_alpha in range(nAlpha):
        exp_term = np.exp(-alpha[i_alpha] * r)
        
        for i_pow in range(nPow):
            coef = aa[i_pow, i_alpha] if aa.ndim > 1 else aa[i_alpha]
            power = ll + i_pow
            
            # Handle r=0 case for l=0, p=1 (r^0 = 1)
            if power == 0:
                r_pow = np.ones_like(r)
            else:
                r_pow = r ** power
            
            result += coef * r_pow * exp_term
    
    return result


def precompute_sto_grid(aa, alpha, ll, cutoff, resolution, n_nodes=None, dr=None):
    """
    Precompute STO values on uniform grid with spline second derivatives.
    
    Args:
        aa: (nPow, nAlpha) coefficients
        alpha: (nAlpha,) exponents
        ll: angular momentum
        cutoff: cutoff radius
        resolution: grid spacing (used if n_nodes/dr not provided)
        n_nodes: optional number of grid points (overrides cutoff/resolution)
        dr: optional grid spacing (overrides resolution)
    
    Returns:
        grid_values: (nNodes,) STO values
        grid_d2: (nNodes,) second derivatives for cubic spline
        dr: grid spacing
    """
    if n_nodes is None or dr is None:
        n_nodes = int(np.ceil(cutoff / resolution)) + 2
        dr = resolution
    
    r = np.arange(n_nodes) * dr
    
    # Evaluate STO at each grid point
    grid_values = compute_sto_radial(r, aa, alpha, ll)
    grid_values = grid_values.astype(np.float32)
    
    # Compute second derivatives (natural cubic spline)
    grid_d2 = _spline_d2_uniform(grid_values, dr)
    
    return grid_values, grid_d2, dr


def _spline_d2_uniform(y, h):
    """Natural cubic spline second derivatives for uniform grid."""
    n = len(y)
    if n < 3:
        return np.zeros(n, dtype=np.float32)
    
    a = np.ones(n-3, dtype=np.float64)
    b = np.full(n-2, 4.0, dtype=np.float64)
    c = np.ones(n-3, dtype=np.float64)
    rhs = np.zeros(n-2, dtype=np.float64)
    rhs[:] = 6.0 * (y[2:] - 2.0*y[1:-1] + y[:-2]) / (h*h)
    
    # Thomas algorithm
    for i in range(1, n-2):
        w = a[i-1] / b[i-1]
        b[i] -= w * c[i-1]
        rhs[i] -= w * rhs[i-1]
    
    d2_inner = np.zeros(n-2, dtype=np.float64)
    d2_inner[-1] = rhs[-1] / b[-1]
    for i in range(n-4, -1, -1):
        d2_inner[i] = (rhs[i] - c[i] * d2_inner[i+1]) / b[i]
    
    d2 = np.zeros(n, dtype=np.float32)
    d2[1:-1] = d2_inner.astype(np.float32)
    
    return d2


# ================================================================
# Standalone parsers for DFTB+ output files (compare_waveplot_lib.py)
# These handle the actual DFTB+ XML format (lowercase tags)
# ================================================================

def parse_detailed_xml_custom(xml_path):
    """
    Parse DFTB+ detailed.xml (actual format with lowercase tags).
    
    Returns dict with:
        species_names: list of unique species in order
        species_per_atom: 0-based index into species_names (natoms,)
        coords_bohr: (natoms, 3) array in Bohr
        nstates, norb, nkpoints, nspin
        occupations: (nstates, nkpoints, nspin)
    """
    import xml.etree.ElementTree as ET
    tree = ET.parse(xml_path)
    root = tree.getroot()
    
    # Parse geometry
    geo = root.find('geometry')
    typenames = geo.find('typenames').text.strip().split()
    # Remove quotes if present
    species_names = [s.strip('"\'') for s in typenames]
    
    # Parse typesandcoordinates: format: type_index x y z (Bohr)
    typesandcoords = geo.find('typesandcoordinates').text.strip().split()
    nlines = len(typesandcoords) // 4
    coords_bohr = []
    species_per_atom = []
    for i in range(nlines):
        idx = int(typesandcoords[i*4]) - 1  # 1-based -> 0-based
        x = float(typesandcoords[i*4 + 1])
        y = float(typesandcoords[i*4 + 2])
        z = float(typesandcoords[i*4 + 3])
        coords_bohr.append([x, y, z])
        species_per_atom.append(idx)
    
    coords_bohr = np.array(coords_bohr)
    species_per_atom = np.array(species_per_atom, dtype=np.int32)
    
    # Parse dimensions
    nstates = int(root.find('nrofstates').text)
    norb = int(root.find('nroforbitals').text)
    nkpoints = int(root.find('nrofkpoints').text)
    nspin = int(root.find('nrofspins').text)
    
    # Parse occupations
    occ = root.find('occupations')
    occ_data = occ.find('spin1').find('k1').text.strip().split()
    occupations_1d = np.array([float(x) for x in occ_data])
    # Reshape to (nstates, nkpoints, nspin) for compatibility
    occupations = occupations_1d.reshape(nstates, nkpoints, nspin)
    
    natoms = coords_bohr.shape[0]
    return {
        'species_names': species_names,
        'species_per_atom': species_per_atom,
        'coords_bohr': coords_bohr,
        'natoms': natoms,
        'nstates': nstates,
        'norb': norb,
        'nkpoints': nkpoints,
        'nspin': nspin,
        'occupations': occupations,
    }


def parse_eigenvec_bin_custom(bin_path, nstates, norb, nkpoints=1, nspin=1):
    """
    Parse DFTB+ eigenvec.bin (flat binary format, no Fortran record markers).
    
    Format: [4-byte identity int] [nstates * norb * 8-byte float64]
    
    Returns: evecs[nstates, norb] (float64)
    """
    import struct
    with open(bin_path, 'rb') as f:
        raw = f.read()
    
    # Read identity (first 4 bytes)
    identity = struct.unpack_from('i', raw, 0)[0]
    
    # Read eigenvectors (rest of file)
    evecs = np.frombuffer(raw[4:], dtype=np.float64).reshape(nstates, norb).copy()
    
    return evecs


def read_cube(path):
    """
    Read Gaussian cube file.
    
    Returns: (grid_data, origin, step, nPoints)
        grid_data: (nx, ny, nz) array
        origin: (3,) array in Bohr
        step: (3,) array in Bohr
        nPoints: (3,) tuple
    """
    import numpy as np
    with open(path, 'r') as f:
        lines = f.readlines()
    
    # Line 1: comment
    # Line 2: comment
    # Line 3: natoms origin_x origin_y origin_z
    natoms = int(lines[2].split()[0])
    origin = np.array([float(x) for x in lines[2].split()[1:4]])
    
    # Lines 4-6: nPoints and step vectors
    nPoints = []
    step = np.zeros(3)
    for i in range(3):
        parts = lines[3+i].split()
        nPoints.append(int(parts[0]))
        step[i] = float(parts[1+i])  # diagonal element only
    
    nPoints = tuple(nPoints)
    
    # Lines 7-7+natoms: atomic numbers and coordinates
    # Skip for now
    
    # Remaining lines: grid data (starts after atom section)
    # Atom section is lines 7-9 (0-indexed: 6-8), so data starts from index 6+natoms
    data_lines = lines[6+natoms:]
    data = []
    for line in data_lines:
        data.extend([float(x) for x in line.split()])
    
    grid_data = np.array(data).reshape(nPoints)
    
    # Parse atom section
    atoms = []
    for i in range(natoms):
        line = lines[6+i].split()
        z = int(line[0])
        x = float(line[1])
        y = float(line[2])
        z_coord = float(line[3])
        atoms.append((z, x, y, z_coord))
    
    return grid_data, origin, step, nPoints, atoms


def build_wp_basis(species_list_ang, species_names_ordered):
    """
    Convert Å-normalized basis (from parse_basis_hsd_ang) back to Bohr for libwaveplot.
    
    Args:
        species_list_ang: list from parse_basis_hsd_ang (Å units)
        species_names_ordered: list of species names in atom order
    
    Returns:
        basis: list of dicts with angMoms, cutoffs (Bohr), occupations, stos
        resoln_b: resolution in Bohr
    """
    B = BOHR2ANG
    sp_by_name = {sp['name']: sp for sp in species_list_ang}
    basis = []
    for sp_name in species_names_ordered:
        sp = sp_by_name[sp_name]
        angMoms   = [orb['l'] for orb in sp['orbitals']]
        cutoffs_b = [orb['cutoff'] / B for orb in sp['orbitals']]
        occs      = [1.0] * len(sp['orbitals'])
        stos = []
        for orb in sp['orbitals']:
            l = orb['l']
            alpha_b = list(np.array(orb['exponents']) * B)       # Å^-1 -> Bohr^-1
            coeffs = np.array(orb['coefficients'])
            nPow = coeffs.shape[0]
            # Reverse the power term scaling: multiply by B^(l+j) for each row
            scale_factors = np.array([B ** (l + j) for j in range(nPow)])
            coef_b = coeffs * scale_factors[:, None]
            aa = coef_b.tolist() if hasattr(coef_b, 'tolist') else list(coef_b)
            stos.append({'alpha': alpha_b, 'aa': aa})
        basis.append({'angMoms': angMoms, 'cutoffs': cutoffs_b,
                      'occupations': occs, 'stos': stos})
    resoln_b = species_list_ang[0]['resolution'] / B
    return basis, resoln_b


def evec_to_kernel_coeffs(evec_row, natoms, species_per_atom, species_names, species_list_ang):
    """
    Convert eigenvector row [norb] -> (natoms,4) kernel coeffs [px,py,pz,s].
    
    Args:
        evec_row: [norb] array
        natoms: number of atoms
        species_per_atom: [natoms] 0-based species indices
        species_names: list of species names
        species_list_ang: list from parse_basis_hsd_ang
    
    Returns:
        coeffs: (natoms, 4) float32 array [px, py, pz, s]
    """
    sp_by_name = {sp['name']: sp for sp in species_list_ang}
    c = np.zeros((natoms, 4), dtype=np.float32)
    offset = 0
    for ia in range(natoms):
        si = species_per_atom[ia]
        sp_name = species_names[si]
        sp = sp_by_name[sp_name]
        for orb in sp['orbitals']:
            l = orb['l']
            nm = 2 * l + 1
            chunk = evec_row[offset:offset+nm]
            if l == 0:
                c[ia, 3] = chunk[0]  # s -> slot 3
            elif l == 1:
                c[ia, 1] = chunk[0]  # py -> slot 1
                c[ia, 2] = chunk[1]  # pz -> slot 2
                c[ia, 0] = chunk[2]  # px -> slot 0
            offset += nm
    return c


def mask_sto_coefficients(species_list_ang, species_name, orbital_idx, active_pow, active_alpha):
    """
    Debug utility: Mask all coefficients to zero except one specific (nPow, nAlpha) term.
    
    For testing individual STO components against SLATER.F90 reference.
    Creates a deep copy of species_list_ang with masked coefficients.
    
    Args:
        species_list_ang: Original species list from parse_basis_hsd_ang
        species_name: Species to mask (e.g., 'H', 'O')
        orbital_idx: Which orbital to mask (0 for first s-orbital, etc.)
        active_pow: Which power term to keep non-zero (0, 1, 2, ...)
        active_alpha: Which exponent to keep non-zero (0, 1, 2, ...)
    
    Returns:
        masked_species_list: Deep copy with masked coefficients
        original_coeff: The original coefficient value at (active_pow, active_alpha)
    """
    import copy
    
    masked_list = copy.deepcopy(species_list_ang)
    original_coeff = None
    
    for sp in masked_list:
        if sp['name'] == species_name:
            if orbital_idx < len(sp['orbitals']):
                orb = sp['orbitals'][orbital_idx]
                coeffs = orb['coefficients'].copy()
                
                # Store original value before masking
                if active_pow < coeffs.shape[0] and active_alpha < coeffs.shape[1]:
                    original_coeff = coeffs[active_pow, active_alpha]
                
                # Mask all to zero
                coeffs[:, :] = 0.0
                
                # Set only active term to 1.0 (unity for testing)
                if active_pow < coeffs.shape[0] and active_alpha < coeffs.shape[1]:
                    coeffs[active_pow, active_alpha] = 1.0
                
                orb['coefficients'] = coeffs
                print(f"[mask_sto_coefficients] {species_name} orbital {orbital_idx}: "
                      f"masked all coeffs to 0, set ({active_pow},{active_alpha})=1.0 "
                      f"(was {original_coeff:.6f})")
    
    return masked_list, original_coeff


def test_single_sto_component(species_list_ang, species_name, orbital_idx, pow_idx, alpha_idx, 
                               r_max_ang=5.0, n_points=100):
    """
    Evaluate a single STO component (one term: r^(l+pow) * exp(-alpha*r)) 
    and compare with analytical formula.
    
    For debugging: each coefficient corresponds to one term in the STO expansion.
    This isolates individual terms for perfect parity testing.
    
    Args:
        species_list_ang: Species list (will be masked)
        species_name: Species to test
        orbital_idx: Which orbital (0 for s, etc.)
        pow_idx: Power index (j in r^(l+j))
        alpha_idx: Exponent index
        r_max_ang: Max radius in Angstrom
        n_points: Number of radial points
    
    Returns:
        r_ang, values: Radial grid and evaluated STO values
    """
    from .DFTBplusParser import evaluate_sto_1d
    
    # Mask to only test this component
    masked_list, orig_coeff = mask_sto_coefficients(
        species_list_ang, species_name, orbital_idx, pow_idx, alpha_idx
    )
    
    # Get the orbital parameters
    sp = None
    for s in masked_list:
        if s['name'] == species_name:
            sp = s
            break
    
    if sp is None:
        raise ValueError(f"Species {species_name} not found")
    
    orb = sp['orbitals'][orbital_idx]
    l = orb['l']
    exps = orb['exponents']
    coeffs = orb['coefficients']
    
    # Create radial grid
    r_ang = np.linspace(0, r_max_ang, n_points)
    r_bohr = r_ang / 0.5291772109
    
    # Evaluate (masked coefficients will give single term)
    values = evaluate_sto_1d(r_bohr, l, exps, coeffs)
    
    # Expected analytical form for this single component:
    # R(r) = r^(l + pow_idx) * exp(-exps[alpha_idx] * r)
    # (coefficient is 1.0 due to masking)
    alpha_bohr = exps[alpha_idx] * 0.5291772109  # convert back to Bohr^-1
    expected = r_bohr**(l + pow_idx) * np.exp(-alpha_bohr * r_bohr)
    
    # Handle r=0 case
    if l == 0 and pow_idx == 0:
        expected[0] = 1.0  # r^0 = 1 at r=0
    
    print(f"[test_single_sto_component] {species_name} l={l} pow={pow_idx} alpha_idx={alpha_idx}")
    print(f"  alpha = {alpha_bohr:.4f} Bohr^-1")
    print(f"  max|evaluated| = {np.abs(values).max():.6f}")
    print(f"  max|expected| = {np.abs(expected).max():.6f}")
    print(f"  RMS diff = {np.sqrt(np.mean((values - expected)**2)):.2e}")
    
    return r_ang, values, expected
