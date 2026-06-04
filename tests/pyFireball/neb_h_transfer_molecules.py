#!/usr/bin/env python3
"""NEB calculation for H-transfer between hydro and dehydro pyrazine molecules.

Correct implementation:
1. Molecules placed using bounding boxes with 2A gap (H-bond distance)
2. Hydrogens on nitrogen atoms in hydro molecule move toward nitrogens in dehydro
3. DFTB+ proper cube file generation using waveplot with manually generated detailed.xml

Periodic boundary conditions along y with gamma-point only.
"""
import sys
import os
import numpy as np
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import dftb_utils as dftbu
from pyBall import atomicUtils as au

# Parameters
Ly = 25.0  # Lattice constant along y (periodic direction)
Lx = 20.0  # Cell size x (vacuum)
Lz = 20.0  # Cell size z (vacuum)
sk_set = '3ob-3-1'
L_Hb = 2.0  # Hydrogen bond distance (A)
n_images = 7  # Number of NEB images (excluding endpoints)
spring_k = 5.0  # Spring constant (eV/A^2)

ELEM_MAP = {'H': 1, 'C': 6, 'N': 7, 'O': 8}
ELEM_Z = {'H': 1, 'C': 6, 'N': 7, 'O': 8}

def read_xyz(fname):
    """Read XYZ file. Returns (apos, atypes, comment)."""
    with open(fname, 'r') as f:
        lines = f.readlines()
    natoms = int(lines[0].strip())
    comment = lines[1].strip()
    apos = np.zeros((natoms, 3))
    atypes = []
    for i in range(natoms):
        parts = lines[2 + i].split()
        atypes.append(parts[0])
        apos[i, 0] = float(parts[1])
        apos[i, 1] = float(parts[2])
        apos[i, 2] = float(parts[3])
    return apos, atypes, comment

def rotate_around_y(apos, angle_deg):
    """Rotate positions around y axis by angle_deg (positive = x->z)."""
    theta = np.radians(angle_deg)
    c, s = np.cos(theta), np.sin(theta)
    R = np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])
    return apos @ R.T

def save_xyz(apos, atypes, fname, comment=""):
    """Save geometry to XYZ file."""
    with open(fname, 'w') as f:
        f.write(f"{len(atypes)}\n")
        f.write(f"{comment}\n")
        for elem, pos in zip(atypes, apos):
            f.write(f"{elem} {pos[0]:.6f} {pos[1]:.6f} {pos[2]:.6f}\n")
    print(f"Saved: {fname}")

def save_xyz_movie(images, atypes, fname, comments=None):
    """Save multiple geometries as XYZ movie (multi-frame XYZ file)."""
    if comments is None:
        comments = [f"Frame {i}" for i in range(len(images))]
    with open(fname, 'w') as f:
        for apos, comment in zip(images, comments):
            f.write(f"{len(atypes)}\n")
            f.write(f"{comment}\n")
            for elem, pos in zip(atypes, apos):
                f.write(f"{elem} {pos[0]:.6f} {pos[1]:.6f} {pos[2]:.6f}\n")
    print(f"Saved XYZ movie: {fname} ({len(images)} frames)")

def build_molecular_cell(hydro_file, dehydro_file, Ly, Lx, Lz, L_Hb=2.0, rotate_y_deg=0.0, use_NN_distance=False):
    """Build periodic cell with hydro and dehydro molecules.

    Places molecules along y-axis. Two modes:
    1. Bounding box mode (default): uses L_Hb gap between bounding boxes
    2. N-N distance mode (use_NN_distance=True): places molecules so donor N and acceptor N are ~3.03 A apart
       (gives H...N ~2 A when H is 1.03 A from donor N)

    If rotate_y_deg is non-zero, only the dehydro molecule is rotated around y axis
    (to prevent CH groups from hindering H-bonds).
    """
    # Read molecules
    apos_hydro, elems_hydro, _ = read_xyz(hydro_file)
    apos_dehydro, elems_dehydro, _ = read_xyz(dehydro_file)

    # Apply rotation to dehydro molecule only
    if rotate_y_deg != 0.0:
        apos_dehydro = rotate_around_y(apos_dehydro, rotate_y_deg)

    # Center each molecule in y
    apos_hydro = apos_hydro.copy()
    apos_dehydro = apos_dehydro.copy()
    apos_hydro[:, 1] -= np.mean(apos_hydro[:, 1])
    apos_dehydro[:, 1] -= np.mean(apos_dehydro[:, 1])

    if use_NN_distance:
        # 2-Hbond-aware placement with PBC in y:
        # We want TWO N..N distances to be D≈3.03 A, with one pair across the boundary.
        # Enforce sign pattern on wrapped dy: one +D, one -D.
        N_h = [i for i,e in enumerate(elems_hydro) if e == 'N']
        N_d = [i for i,e in enumerate(elems_dehydro) if e == 'N']
        if len(N_h) != 2 or len(N_d) != 2:
            raise ValueError(f"use_NN_distance requires exactly 2 N in each molecule; got hydro={len(N_h)} dehydro={len(N_d)}")
        # sort by y (lo, hi)
        N_h = sorted(N_h, key=lambda i: apos_hydro[i,1])
        N_d = sorted(N_d, key=lambda i: apos_dehydro[i,1])
        y_h = np.array([apos_hydro[i,1] for i in N_h])
        y_d = np.array([apos_dehydro[i,1] for i in N_d])
        D = 3.03

        def wrap(dy):
            return dy - np.round(dy / Ly) * Ly

        best = (1e9, 0.0, None)
        # two assignments: direct (lo->lo, hi->hi) and cross (lo->hi, hi->lo)
        for shift_y in np.linspace(0.0, Ly, 4001):
            yy = y_d + shift_y
            # direct
            dy_lo = wrap(yy[0]-y_h[0])
            dy_hi = wrap(yy[1]-y_h[1])
            # cross
            dy_lo_c = wrap(yy[1]-y_h[0])
            dy_hi_c = wrap(yy[0]-y_h[1])

            # cost with desired sign pattern (+D for upper, -D for lower)
            c_dir = (dy_lo + D)**2 + (dy_hi - D)**2
            c_crs = (dy_lo_c + D)**2 + (dy_hi_c - D)**2

            if c_dir <= c_crs:
                cost, assign = c_dir, [(N_h[0], N_d[0], dy_lo), (N_h[1], N_d[1], dy_hi)]
            else:
                cost, assign = c_crs, [(N_h[0], N_d[1], dy_lo_c), (N_h[1], N_d[0], dy_hi_c)]

            if cost < best[0]:
                best = (cost, shift_y, assign)

        _, shift_y, assign = best
        apos_dehydro[:,1] += shift_y
        print(f"NN-PBC placement (fixed Ly={Ly:.3f}): shift_y={shift_y:.3f} A")
        for ih, idd, dy in assign:
            print(f"  N_h[{ih}] -> N_d[{idd}]  dy_min_image={dy:.3f}  target={'-' if dy<0 else '+'}{D:.3f}")
    else:
        # Bounding box placement
        y_min_hydro = np.min(apos_hydro[:, 1])
        y_max_hydro = np.max(apos_hydro[:, 1])
        y_min_dehydro = np.min(apos_dehydro[:, 1])
        y_max_dehydro = np.max(apos_dehydro[:, 1])
        shift_y = (y_max_hydro + L_Hb) - y_min_dehydro
        apos_dehydro[:, 1] += shift_y
        Ly = (y_max_hydro - y_min_hydro) + (y_max_dehydro - y_min_dehydro) + 2*L_Hb
        print(f"Bounding box placement: shift_y={shift_y:.3f} A, gap={L_Hb} A")

    # Calculate cell size (for NN placement Ly already set)
    if not use_NN_distance:
        pass

    # Combine molecules
    apos = np.vstack([apos_hydro, apos_dehydro])
    elems = elems_hydro + elems_dehydro
    mol_id = np.array([0]*len(elems_hydro) + [1]*len(elems_dehydro), dtype=np.int32)

    # Center in cell
    apos[:, 1] -= np.mean(apos[:, 1])
    apos[:, 1] += 0.5 * Ly
    apos[:, 0] += Lx / 2.0
    apos[:, 2] += Lz / 2.0

    lvs = np.array([[Lx, 0.0, 0.0], [0.0, Ly, 0.0], [0.0, 0.0, Lz]])

    print(f"Built molecular cell: Natoms={len(elems)}, Ly={Ly:.3f} A")
    return apos, elems, lvs, mol_id

def find_hydrogens_on_nitrogen(apos, elems):
    """Find H atoms bonded to N atoms and the N atoms they bond to.

    Returns:
        h_indices: list of H atom indices bonded to N
        n_donor_indices: list of N atom indices that the H atoms bond to
        n_acceptor_indices: list of N atom indices in the other molecule
    """
    # Find N atoms
    n_indices = [i for i, e in enumerate(elems) if e == 'N']
    h_indices = [i for i, e in enumerate(elems) if e == 'H']

    print(f"Found {len(h_indices)} H atoms at indices: {h_indices}")
    print(f"Found {len(n_indices)} N atoms at indices: {n_indices}")

    # For each H, find closest N. If distance < 1.5 A, consider it bonded
    h_on_n = []
    n_donors = []
    for h_idx in h_indices:
        h_pos = apos[h_idx]
        min_dist = float('inf')
        closest_n = None
        for n_idx in n_indices:
            n_pos = apos[n_idx]
            dist = np.linalg.norm(h_pos - n_pos)
            if dist < min_dist:
                min_dist = dist
                closest_n = n_idx
        if min_dist < 1.5:  # N-H bond distance threshold
            h_on_n.append(h_idx)
            n_donors.append(closest_n)
            print(f"  H[{h_idx}] bonded to N[{closest_n}] at distance {min_dist:.3f} A")

    # Separate into two groups by y position (hydro = lower y, dehydro = higher y)
    if len(n_indices) >= 2:
        n_y = [apos[i, 1] for i in n_indices]
        # Sort by y position
        sorted_n = sorted(zip(n_indices, n_y), key=lambda x: x[1])
        mid = len(sorted_n) // 2
        n_lower = [idx for idx, y in sorted_n[:mid]]  # Hydro molecule (lower y)
        n_upper = [idx for idx, y in sorted_n[mid:]]    # Dehydro molecule (upper y)
    else:
        n_lower = n_indices[:len(n_indices)//2]
        n_upper = n_indices[len(n_indices)//2:]

    print(f"  N in hydro (lower y): {n_lower}")
    print(f"  N in dehydro (upper y): {n_upper}")

    return h_on_n, n_donors, n_lower, n_upper

def _pbc_vec_y(d, Ly_cell):
    d = d.copy(); d[1] -= np.round(d[1] / Ly_cell) * Ly_cell; return d

def find_hydrogen_bonds(apos, elems, lvs, mol_id, Xelems=('N','O','F'), r_XH=1.3, r_HY=(1.3, 2.6), r_XY=6.0):
    """Detect hydrogen bonds X–H...Y with PBC along y.

    Donor X is chosen ONLY among electronegative atoms in the same molecule as H using direct (non-PBC) distance.
    Acceptor Y is chosen among electronegative atoms in the other molecule using minimum-image along y.

    Returns list of dicts: {X,H,Y,imgY, dXH,dHY,dXY, score}.
    """
    Ly_cell = lvs[1, 1]
    Xidx = [i for i,e in enumerate(elems) if e in Xelems]
    Hidx = [i for i,e in enumerate(elems) if e == 'H']
    hbonds = []
    for h in Hidx:
        # donor must be in same molecule as H; do NOT use PBC for covalent bond
        bestX, dXH = None, 1e9
        for x in Xidx:
            if mol_id[x] != mol_id[h]: continue
            r = np.linalg.norm(apos[h] - apos[x])
            if r < dXH:
                dXH, bestX = r, x
        if bestX is None or dXH > r_XH: continue

        # acceptor: other molecule, minimum-image along y
        bestY, bestImgY, bestdXY, bestdHY = None, 0, 1e9, 1e9
        for y in Xidx:
            if mol_id[y] == mol_id[bestX]: continue
            dy = apos[y,1] - apos[bestX,1]
            imgY = int(-np.round(dy / Ly_cell))
            y_pos = apos[y].copy(); y_pos[1] += imgY * Ly_cell
            x_pos = apos[bestX].copy()
            # unwrap donor close to acceptor image
            x_un = x_pos.copy();
            dx = x_un[1] - y_pos[1]
            x_un[1] -= np.round(dx / Ly_cell) * Ly_cell
            dXY = np.linalg.norm(y_pos - x_un)
            if dXY > r_XY: continue

            h_pos = apos[h].copy();
            dhx = h_pos[1] - x_un[1]
            h_pos[1] -= np.round(dhx / Ly_cell) * Ly_cell
            dHY = np.linalg.norm(y_pos - h_pos)

            # prefer acceptors that give H...Y in window, but don't require it (phenazine initial may be longer)
            pen = 0.0
            if dHY < r_HY[0]: pen += (r_HY[0] - dHY)
            if dHY > r_HY[1]: pen += 0.2*(dHY - r_HY[1])
            score = dXY + 0.1*dHY + pen
            if score < bestdXY + 0.1*bestdHY:
                bestY, bestImgY, bestdXY, bestdHY = y, imgY, dXY, dHY

        if bestY is None: continue
        hbonds.append({'X':bestX,'H':h,'Y':bestY,'imgY':bestImgY,'dXH':dXH,'dHY':bestdHY,'dXY':bestdXY,'score':bestdXY + 0.1*bestdHY})

    hbonds.sort(key=lambda b: b['score'])
    return hbonds

def create_initial_final_states(apos, elems, lvs, mol_id):
    """Create initial and final states for H-transfer based on detected H-bonds."""
    Ly_cell = lvs[1, 1]
    hb = find_hydrogen_bonds(apos, elems, lvs, mol_id)
    if len(hb) == 0:
        raise ValueError('No hydrogen bonds (X–H...Y) detected')

    # group by donor X, keep best acceptor per donor-H
    usedH=set(); usedX=set(); pairs=[]
    for b in hb:
        if b['H'] in usedH: continue
        if b['X'] in usedX: continue
        usedH.add(b['H']); usedX.add(b['X']); pairs.append(b)
        if len(pairs) == 2: break
    if len(pairs) == 0:
        raise ValueError('Failed to select H-bond pairs')

    NH_BOND = 1.03
    apos_initial = apos.copy(); save_xyz(apos_initial, elems, 'initial_mol.xyz', 'Initial: H on donor X')
    apos_final   = apos.copy()

    h_to_move=[]; n_targets=[]; n_donors=[]; n_target_images=[]
    print(f"\nH-transfer pairs (Ly={Ly_cell:.3f} A):")
    for b in pairs:
        x,h,y,imgY = b['X'], b['H'], b['Y'], b['imgY']
        x_pos = apos[x].copy()
        y_pos = apos[y].copy(); y_pos[1] += imgY * Ly_cell
        # unwrap donor close to acceptor? keep as is; construct vector from acceptor->donor in same image as acceptor
        x_un = x_pos.copy();
        dx = x_un[1] - y_pos[1]
        x_un[1] -= np.round(dx / Ly_cell) * Ly_cell
        vec = x_un - y_pos; vec /= np.linalg.norm(vec)
        h_final = y_pos + vec * NH_BOND
        h_final[1] = h_final[1] % Ly_cell
        apos_final[h] = h_final
        print(f"  X[{x}]--H[{h}]...Y[{y}] imgY={imgY:+d}  dXH={b['dXH']:.3f} dHY={b['dHY']:.3f} dXY={b['dXY']:.3f}  -> H_y {apos[h,1]:.3f}->{h_final[1]:.3f}")
        h_to_move.append(h); n_targets.append(y); n_donors.append(x); n_target_images.append(imgY)

    save_xyz(apos_final, elems, 'final_mol.xyz', 'Final: H on acceptor Y')
    return apos_initial, apos_final, h_to_move, n_targets, n_donors, n_target_images

def generate_orbital_cubes_pyscf(apos, elems, work_dir, grid_resolution=0.2, margin=4.0,
                                   orbital_indices=None, basis='sto-3g'):
    """Generate cube files for molecular orbitals using pySCF.

    This is a fallback when DFTB+ waveplot/libwaveplot have library issues.
    Computes orbitals for the isolated molecule using pySCF.

    Args:
        apos: (natoms, 3) positions in Angstrom
        elems: list of element names
        work_dir: working directory for output
        grid_resolution: grid spacing in Angstrom
        margin: margin around molecule in Angstrom
        orbital_indices: list of orbital indices to plot (None = HOMO-1, HOMO, LUMO, LUMO+1)
        basis: pySCF basis set (default 'sto-3g' for speed, use 'ccpvdz' for accuracy)

    Returns:
        list of generated cube file paths
    """
    import pyscf
    from pyscf import gto, scf, tools

    print(f"  Generating orbital cubes using pySCF ({basis})...")

    # Build pySCF molecule
    atom_str = ""
    for elem, pos in zip(elems, apos):
        atom_str += f"{elem} {pos[0]:.10f} {pos[1]:.10f} {pos[2]:.10f};"

    mol = gto.M(atom=atom_str[:-1], basis=basis, verbose=0)

    # Run SCF
    mf = scf.RHF(mol)
    mf.kernel()

    # Determine which orbitals to plot
    nmo = mf.mo_coeff.shape[1]
    nocc = int(mol.nelectron // 2)
    if orbital_indices is None:
        # Plot HOMO-1, HOMO, LUMO, LUMO+1
        orbital_indices = [max(0, nocc - 2), max(0, nocc - 1),
                           min(nmo - 1, nocc), min(nmo - 1, nocc + 1)]

    # Generate grid
    min_pos = np.min(apos, axis=0) - margin
    max_pos = np.max(apos, axis=0) + margin
    ngrid = np.ceil((max_pos - min_pos) / grid_resolution).astype(np.int32)
    ngrid = np.maximum(ngrid, 2)

    # Generate cube files for each orbital
    cube_files = []
    os.makedirs(work_dir, exist_ok=True)

    for imo in orbital_indices:
        cube_name = f"orbital_{imo:03d}.cube"
        cube_path = os.path.join(work_dir, cube_name)

        # pySCF's cubegen
        tools.cubegen.orbital(mol, cube_path, mf.mo_coeff[:, imo],
                              nx=ngrid[0], ny=ngrid[1], nz=ngrid[2])

        e_mo = mf.mo_energy[imo]
        cube_files.append(cube_path)
        print(f"    Generated {cube_name} (E={e_mo:.6f} Hartree = {e_mo*27.2114:.3f} eV)")

    return cube_files


def _write_cube_file(fname, apos, elems, grid_spec, data):
    """Write Gaussian cube file (helper for any grid-based method)."""
    natoms = len(elems)
    origin = grid_spec['origin']
    dA = grid_spec['dA']
    dB = grid_spec['dB']
    dC = grid_spec['dC']
    nx, ny, nz = grid_spec['ngrid']

    BOHR = 1.8897259886
    origin_bohr = origin * BOHR
    dA_bohr = dA * BOHR
    dB_bohr = dB * BOHR
    dC_bohr = dC * BOHR

    with open(fname, 'w') as f:
        f.write("Orbital from DFTB+\n")
        f.write("Generated by neb_h_transfer_molecules.py\n")
        f.write(f"{natoms} {origin_bohr[0]:.6f} {origin_bohr[1]:.6f} {origin_bohr[2]:.6f}\n")
        f.write(f"{nx} {dA_bohr[0]:.6f} {dA_bohr[1]:.6f} {dA_bohr[2]:.6f}\n")
        f.write(f"{ny} {dB_bohr[0]:.6f} {dB_bohr[1]:.6f} {dB_bohr[2]:.6f}\n")
        f.write(f"{nz} {dC_bohr[0]:.6f} {dC_bohr[1]:.6f} {dC_bohr[2]:.6f}\n")

        elem_z = {'H': 1, 'C': 6, 'N': 7, 'O': 8}
        for i, (elem, pos) in enumerate(zip(elems, apos)):
            z = elem_z.get(elem, 0)
            pos_bohr = pos * BOHR
            f.write(f"{z} 0.0 {pos_bohr[0]:.6f} {pos_bohr[1]:.6f} {pos_bohr[2]:.6f}\n")

        count = 0
        for ix in range(nx):
            for iy in range(ny):
                for iz in range(nz):
                    f.write(f" {data[ix, iy, iz]:.6e}")
                    count += 1
                    if count % 6 == 0:
                        f.write("\n")
        if count % 6 != 0:
            f.write("\n")

def _generate_detailed_xml(workdir, elems, apos):
    """Generate minimal detailed.xml for waveplot from geometry."""
    import xml.etree.ElementTree as ET

    # Count electrons to estimate nstates
    elem_z = {'H': 1, 'C': 6, 'N': 7, 'O': 8}
    nelectrons = sum(elem_z.get(e, 0) for e in elems)
    nstates = nelectrons // 2  # Closed shell

    # Count orbitals
    norb = 0
    for e in elems:
        if e == 'H':
            norb += 1
        elif e in ['C', 'N']:
            norb += 4

    # Get unique species
    species_names = []
    species_per_atom = []
    for e in elems:
        if e not in species_names:
            species_names.append(e)
        species_per_atom.append(species_names.index(e))

    root = ET.Element('detailed')
    geometry = ET.SubElement(root, 'geometry')
    typenames = ET.SubElement(geometry, 'typenames')
    typenames.text = ' '.join(f'"{s}"' for s in species_names)

    typesandcoords = ET.SubElement(geometry, 'typesandcoordinates')
    coords_text = []
    for sp_idx, pos in zip(species_per_atom, apos):
        pos_bohr = pos * 1.889726133873
        coords_text.append(f"{sp_idx + 1} {pos_bohr[0]:.10f} {pos_bohr[1]:.10f} {pos_bohr[2]:.10f}")
    typesandcoords.text = '\n' + '\n'.join(coords_text) + '\n'

    ET.SubElement(root, 'nrofstates').text = str(nstates)
    ET.SubElement(root, 'nroforbitals').text = str(norb)
    ET.SubElement(root, 'nrofkpoints').text = '1'
    ET.SubElement(root, 'nrofspins').text = '1'

    occupations = ET.SubElement(root, 'occupations')
    spin1 = ET.SubElement(occupations, 'spin1')
    k1 = ET.SubElement(spin1, 'k1')
    occ_values = []
    for i in range(nstates):
        occ_values.append("2.0000000000" if i < nstates // 2 else "0.0000000000")
    k1.text = ' '.join(occ_values)

    ET.ElementTree(root).write(os.path.join(workdir, 'detailed.xml'), encoding='utf-8', xml_declaration=True)


def run_dftb_calculation(apos, elems, lvs, temp_dir, do_forces=True, do_orbitals=False, opt=False, fixed_atoms=None,
                         Temperature=800, MixingParameter=0.05, MaxScc=2000, SCCTolerance=1e-6, reuse_charges=False):
    """Run DFTB+ calculation using dftb_utils with gamma point only.

    Notes on SCC stability:
    - Temperature controls Fermi broadening (helps near-degeneracies)
    - MixingParameter too large can cause SCC divergence; default here is conservative
    - reuse_charges=False avoids reusing a bad initial guess between scan points
    """
    cwd = os.getcwd()
    os.makedirs(temp_dir, exist_ok=True)
    os.chdir(temp_dir)

    try:
        # Write periodic DFTB+ input with gamma point only
        dftbu.makeDFTBjob_pbc(
            enames=elems, apos=apos, lvs=lvs, fname='dftb_in.hsd',
            sk_set=sk_set,
            nk=(1, 1, 1), k_shift=(0.0, 0.0, 0.0),
            opt=opt, params=dftbu.default_params,
            Temperature=Temperature, MixingParameter=MixingParameter,
            MaxScc=MaxScc, SCCTolerance=SCCTolerance,
            fixed_atoms=fixed_atoms
        )

        # Remove ParserOptions block that causes version incompatibility
        with open('dftb_in.hsd', 'r') as f:
            content = f.read()
        content = content.replace('ParserOptions {\n  ParserVersion = 15\n}\n', '')
        with open('dftb_in.hsd', 'w') as f:
            f.write(content)

        # Add analysis options
        with open('dftb_in.hsd', 'a') as f:
            f.write("\nAnalysis {\n")
            f.write("  CalculateForces = Yes\n")
            if do_orbitals:
                f.write("  WriteEigenvectors = Yes\n")
            f.write("}\n")
            f.write("\nOptions {\n")
            f.write("  WriteDetailedOut = Yes\n")
            f.write("}\n")

        # Run DFTB+
        if not reuse_charges:
            if os.path.exists('charges.bin'): os.remove('charges.bin')
        dftb_path = '/home/prokophapala/miniconda3/bin/dftb+'
        ierr = os.system(f'{dftb_path} > OUT 2> ERR')
        if ierr != 0:
            raise RuntimeError(f"DFTB+ failed with exit code {ierr}")
        
        # Generate cube files using correct waveplot binary if eigenvec.bin exists
        if do_orbitals and os.path.exists('eigenvec.bin'):
            try:
                # Use correct waveplot binary (dftbplus version, not dftb-asi)
                wp_path = '/home/prokophapala/opt/dftbplus/bin/waveplot'
                if os.path.exists(wp_path):
                    # Generate minimal detailed.xml for waveplot
                    _generate_detailed_xml('.', elems, apos)
                    # Run waveplot with correct binary
                    import subprocess
                    result = subprocess.run([wp_path], cwd='.', capture_output=True, text=True)
                    if result.returncode == 0:
                        cube_files = [f for f in os.listdir('.') if f.endswith('.cube')]
                        print(f"  Generated {len(cube_files)} cube files with waveplot")
                    else:
                        print(f"  waveplot stderr: {result.stderr[:200]}")
                else:
                    print(f"  waveplot not found at {wp_path}")
            except Exception as e:
                print(f"  waveplot failed: {e}")

        # Parse energy from detailed.out — last column is already eV
        E = None
        if os.path.exists('detailed.out'):
            for line in open('detailed.out'):
                if line.startswith('Total energy:'):
                    E = float(line.split()[-2])   # eV column (second to last token)
                    break
        if E is None:
            raise ValueError("Could not parse energy from DFTB+ output")

        # Parse forces if requested
        forces = None
        if do_forces and os.path.exists('detailed.out'):
            forces = parse_forces_from_detailed_out('detailed.out', len(elems))

        # Parse orbital energies from band.out
        orbitals = None
        if do_orbitals and os.path.exists('band.out'):
            orbitals = parse_orbitals_from_band_out('band.out')

        # Generate cube files if requested (for isolated molecules only, using pySCF)
        cube_files = []
        if do_orbitals:
            try:
                # Use pySCF for cube generation (avoids DFTB+ library issues)
                cube_files = generate_orbital_cubes_pyscf(
                    apos, elems, '.', grid_resolution=0.2, margin=4.0
                )
            except Exception as e:
                print(f"  Cube file generation failed: {e}")
                import traceback
                traceback.print_exc()

        return E, forces, orbitals, cube_files

    except Exception as e:
        print(f"  ERROR: {e}")
        return None, None, None, None
    finally:
        os.chdir(cwd)

def parse_forces_from_detailed_out(fname, natoms):
    """Parse forces from DFTB+ detailed.out file."""
    forces = np.zeros((natoms, 3))
    with open(fname, 'r') as f:
        lines = f.readlines()
    in_forces = False
    atom_idx = 0
    for line in lines:
        if 'Total Forces' in line:
            in_forces = True
            continue
        if in_forces and atom_idx < natoms:
            parts = line.split()
            if len(parts) >= 3:
                try:
                    forces[atom_idx] = [float(parts[0]), float(parts[1]), float(parts[2])]
                    atom_idx += 1
                except ValueError:
                    continue
    return forces

def parse_orbitals_from_band_out(fname):
    """Parse orbital energies from DFTB+ band.out file."""
    energies = []
    with open(fname, 'r') as f:
        lines = f.readlines()
    for line in lines:
        parts = line.split()
        if len(parts) >= 2:
            try:
                energies.append(float(parts[1]))
            except ValueError:
                continue
    return np.array(energies)

def linear_interpolate(apos1, apos2, n_points, lvs=None):
    """Create linear interpolation between two geometries.

    If lvs is provided, uses PBC-aware minimum-image displacement along y (lvs[1,1])
    so moving atoms always take the short path, not the long way round the cell.
    """
    images = []
    if lvs is not None:
        Ly_cell = lvs[1, 1]
        # displacement wrapped to [-Ly/2, +Ly/2] — minimum image
        d = apos2 - apos1
        d[:, 1] -= np.round(d[:, 1] / Ly_cell) * Ly_cell
    else:
        d = apos2 - apos1
    for i in range(n_points):
        t = i / (n_points - 1) if n_points > 1 else 0
        apos = apos1 + t * d
        if lvs is not None:
            apos[:, 1] = apos[:, 1] % Ly_cell   # wrap y back into [0, Ly)
        images.append(apos.copy())
    return images

def compute_neb_forces(images, energies, forces_list):
    """Compute NEB forces for all images."""
    n_images = len(images)
    neb_forces = []
    for i in range(n_images):
        f_real = forces_list[i].copy()
        if i == 0 or i == n_images - 1:
            neb_forces.append(np.zeros_like(f_real))
            continue
        tau = images[i+1] - images[i-1]
        tau_norm = np.linalg.norm(tau)
        if tau_norm > 1e-10:
            tau = tau / tau_norm
        f_spring = spring_k * (np.linalg.norm(images[i+1] - images[i]) -
                               np.linalg.norm(images[i] - images[i-1])) * tau
        f_perp = f_real - np.dot(f_real.flatten(), tau.flatten()) * tau
        f_neb = f_perp + f_spring
        neb_forces.append(f_neb)
    return neb_forces

def run_neb_calculation(hydro_file, dehydro_file):
    """Run NEB calculation for hydrogen transfer between pyrazine molecules."""
    print("Building molecular cell...")
    apos, elems, lvs, mol_id = build_molecular_cell(hydro_file, dehydro_file, Ly, Lx, Lz, L_Hb)

    print("\nCreating initial and final states...")
    apos_initial, apos_final, h_to_move, n_targets, n_donors, n_target_images = create_initial_final_states(apos, elems, lvs, mol_id)

    print("\nSetting up NEB images...")
    n_total = n_images + 2
    images = linear_interpolate(apos_initial, apos_final, n_total, lvs=lvs)

    for i, img in enumerate(images):
        save_xyz(img, elems, f'neb_image_{i:02d}_init_mol.xyz', f"NEB image {i} initial")

    print(f"\nRunning NEB optimization with {n_total} images...")
    print("(Note: Using gamma-point only, periodic along y)")

    energies = []
    forces_list = []
    orbitals_list = []
    cube_files_list = []

    for i, img in enumerate(images):
        print(f"\nCalculating image {i}/{n_total-1}...")
        temp_dir = f'temp_neb_mol_{i}'
        E, forces, orbitals, cube_files = run_dftb_calculation(
            img, elems, lvs, temp_dir, do_forces=True, do_orbitals=True)
        if E is not None:
            energies.append(E)
            forces_list.append(forces)
            orbitals_list.append(orbitals)
            cube_files_list.append(cube_files)
            print(f"  E = {E:.6f} eV")
            if orbitals is not None and len(orbitals) > 1:
                print(f"  HOMO-LUMO gap: {(orbitals[0] - orbitals[-1]):.6f} eV")
            if cube_files:
                print(f"  Generated {len(cube_files)} cube files")
        else:
            print("  FAILED")

    if len(energies) == n_total:
        print("\nNEB calculation completed successfully!")
        print(f"Energy barrier: {max(energies) - min(energies):.6f} eV")

        plt.figure(figsize=(10, 6))
        plt.plot(range(n_total), energies, 'bo-', linewidth=2, markersize=8)
        plt.xlabel('NEB Image', fontsize=14)
        plt.ylabel('Energy (eV)', fontsize=14)
        plt.title('Hydrogen Transfer Energy Profile (Pyrazine)', fontsize=16)
        plt.grid(True, alpha=0.3)
        plt.savefig('neb_energy_profile_mol.png', dpi=300, bbox_inches='tight')
        print("Saved energy profile: neb_energy_profile_mol.png")

        comments = [f"Image {i}: E = {energies[i]:.6f} eV" for i in range(n_total)]
        save_xyz_movie(images, elems, 'neb_h_transfer_movie.xyz', comments)

        np.savez('neb_results_mol.npz',
                 energies=energies,
                 images=images,
                 elems=elems,
                 lvs=lvs)
        print("Saved results: neb_results_mol.npz")
    else:
        print(f"\nNEB calculation incomplete: {len(energies)}/{n_total} images succeeded")

def _parse_fixed_atoms(s, natoms=None):
    if s is None or s.strip() == "": return []
    idx = [int(x) for x in s.replace(',', ' ').split()]
    if natoms is not None:
        for i in idx:
            if i < 0 or i >= natoms: raise ValueError(f"fixed atom index out of range: {i}")
    return sorted(set(idx))

def _read_xyz_last_frame(fname):
    with open(fname, 'r') as f:
        frames = []
        while True:
            line = f.readline()
            if not line: break
            n = int(line.strip())
            f.readline()
            apos = np.zeros((n, 3)); elems = []
            for i in range(n):
                p = f.readline().split(); elems.append(p[0]); apos[i] = (float(p[1]), float(p[2]), float(p[3]))
            frames.append((apos, elems))
    if len(frames) == 0: raise ValueError(f"no frames in {fname}")
    return frames[-1]

def run_rigid_scan(hydro_file, dehydro_file, n_points=11, Temperature=800, MixingParameter=0.05, MaxScc=2000, SCCTolerance=1e-6, rotate_y_deg=0.0, use_NN_distance=False):
    """Run rigid scan: single-point energies along interpolated path."""
    print("Building molecular cell...")
    apos, elems, lvs, mol_id = build_molecular_cell(hydro_file, dehydro_file, Ly, Lx, Lz, L_Hb, rotate_y_deg, use_NN_distance)

    print("\nCreating initial and final states...")
    apos_initial, apos_final, h_to_move, n_targets, n_donors, n_target_images = create_initial_final_states(apos, elems, lvs, mol_id)

    print(f"\nSetting up rigid scan with {n_points} points...")
    images = linear_interpolate(apos_initial, apos_final, n_points, lvs=lvs)

    comments = [f"Image {i}/{n_points-1}" for i in range(n_points)]
    save_xyz_movie(images, elems, 'rigid_scan_movie.xyz', comments)

    print(f"\nRunning single-point DFTB+ calculations...")
    energies = []
    for i, img in enumerate(images):
        print(f"  Point {i+1}/{n_points}...", end=' ', flush=True)
        temp_dir = f'temp_scan_{i:03d}'
        E, _, _, _ = run_dftb_calculation(img, elems, lvs, temp_dir, do_forces=False, do_orbitals=False,
                                          Temperature=Temperature, MixingParameter=MixingParameter, MaxScc=MaxScc, SCCTolerance=SCCTolerance,
                                          reuse_charges=False)
        if E is not None:
            energies.append(E)
            print(f"E = {E:.6f} eV")
        else:
            print("FAILED")

    if len(energies) == n_points:
        print(f"\nRigid scan completed successfully!")
        print(f"Energy barrier: {max(energies) - min(energies):.6f} eV")
        plt.figure(figsize=(10, 6))
        reaction_coord = np.linspace(0, 1, n_points)
        plt.plot(reaction_coord, energies, 'bo-', linewidth=2, markersize=8)
        plt.xlabel('Reaction Coordinate', fontsize=14)
        plt.ylabel('Energy (eV)', fontsize=14)
        plt.title('Hydrogen Transfer Rigid Scan (Pyrazine)', fontsize=16)
        plt.grid(True, alpha=0.3)
        plt.savefig('rigid_scan_energy.png', dpi=300, bbox_inches='tight')
        print("Saved energy profile: rigid_scan_energy.png")

        np.savez('rigid_scan_results.npz', energies=np.array(energies), reaction_coord=reaction_coord, images=images, elems=elems, lvs=lvs,
                 fixed_donors=np.array(n_donors), moving_H=np.array(h_to_move), acceptors=np.array(n_targets))
        print("Saved results: rigid_scan_results.npz")
    else:
        print(f"\nRigid scan incomplete: {len(energies)}/{n_points} points succeeded")

def run_relax_scan(hydro_file, dehydro_file, n_points=11, Temperature=800, MixingParameter=0.05, MaxScc=2000, SCCTolerance=1e-6, fixed_extra=None, rotate_y_deg=0.0, use_NN_distance=False):
    """Relax scan: optimize each image with fixed atoms (default donor X + moving H)."""
    print("Building molecular cell...")
    apos, elems, lvs, mol_id = build_molecular_cell(hydro_file, dehydro_file, Ly, Lx, Lz, L_Hb, rotate_y_deg, use_NN_distance)

    print("\nCreating initial and final states...")
    apos_initial, apos_final, h_to_move, n_targets, n_donors, n_target_images = create_initial_final_states(apos, elems, lvs, mol_id)

    print(f"\nSetting up relax scan with {n_points} points...")
    images0 = linear_interpolate(apos_initial, apos_final, n_points, lvs=lvs)

    fixed0 = sorted(set((fixed_extra or []) + list(h_to_move) + list(n_donors)))
    print(f"Fixed atoms (0-based): {fixed0}")

    energies = []
    images_relaxed = []
    for i, img in enumerate(images0):
        print(f"\n  Relax point {i+1}/{n_points} ...")
        temp_dir = f'temp_relax_{i:03d}'
        E, _, _, _ = run_dftb_calculation(img, elems, lvs, temp_dir, do_forces=False, do_orbitals=False,
                                          opt=True, fixed_atoms=fixed0,
                                          Temperature=Temperature, MixingParameter=MixingParameter, MaxScc=MaxScc, SCCTolerance=SCCTolerance,
                                          reuse_charges=False)
        if E is None:
            print("  FAILED")
            continue
        energies.append(E)
        xyz = os.path.join(temp_dir, 'geom.out.xyz')
        if os.path.exists(xyz):
            apos_last, _ = _read_xyz_last_frame(xyz)
            images_relaxed.append(apos_last)
        else:
            images_relaxed.append(img.copy())
        print(f"  E = {E:.6f} eV")

    if len(energies) == n_points:
        comments = [f"Relax {i}: E={energies[i]:.6f} eV" for i in range(n_points)]
        save_xyz_movie(images_relaxed, elems, 'relax_scan_movie.xyz', comments)
        plt.figure(figsize=(10, 6))
        rc = np.linspace(0, 1, n_points)
        plt.plot(rc, energies, 'bo-', linewidth=2, markersize=8)
        plt.xlabel('Reaction Coordinate', fontsize=14)
        plt.ylabel('Energy (eV)', fontsize=14)
        plt.title('Hydrogen Transfer Relax Scan (fixed donor N + H)', fontsize=16)
        plt.grid(True, alpha=0.3)
        plt.savefig('relax_scan_energy.png', dpi=300, bbox_inches='tight')
        np.savez('relax_scan_results.npz', energies=np.array(energies), reaction_coord=rc, images=images_relaxed, elems=elems, lvs=lvs, fixed=np.array(fixed0))
        print("\nRelax scan completed: relax_scan_energy.png, relax_scan_movie.xyz, relax_scan_results.npz")
    else:
        print(f"\nRelax scan incomplete: {len(energies)}/{n_points} points succeeded")

def run_dry_run(hydro_file, dehydro_file, rotate_y_deg=0.0, use_NN_distance=False, n_points=11):
    """Geometry-only check: build cell, detect H-bonds, save initial/final states and interpolated movie without DFTB."""
    print("Building molecular cell...")
    apos, elems, lvs, mol_id = build_molecular_cell(hydro_file, dehydro_file, Ly, Lx, Lz, L_Hb, rotate_y_deg, use_NN_distance)

    print("\nCreating initial and final states...")
    apos_initial, apos_final, h_to_move, n_targets, n_donors, n_target_images = create_initial_final_states(apos, elems, lvs, mol_id)

    save_xyz(apos_initial, elems, 'dry_run_initial.xyz', 'Dry run: initial state')
    save_xyz(apos_final, elems, 'dry_run_final.xyz', 'Dry run: final state')

    # Interpolate and save movie (same as rigid_scan)
    Ly_cell = lvs[1, 1]
    images = linear_interpolate(apos_initial, apos_final, n_points, lvs=lvs)
    comments = [f"Frame {i}/{n_points-1}" for i in range(n_points)]
    save_xyz_movie(images, elems, 'dry_run_movie.xyz', comments)
    print(f"Saved interpolated movie: dry_run_movie.xyz ({n_points} frames)")

    print(f"\nGeometry check (Ly={Ly_cell:.3f} A):")
    print(f"  Natoms: {len(elems)}")
    print(f"  Fixed donor X: {n_donors}")
    print(f"  Moving H: {h_to_move}")
    print(f"  Acceptors Y: {n_targets}")

    # Check min distances
    def min_dist_pbc(a, b):
        d = a - b; d[1] -= np.round(d[1] / Ly_cell) * Ly_cell; return np.linalg.norm(d)

    min_HH = 1e9; min_HX = 1e9; min_XX = 1e9
    for i in range(len(elems)):
        for j in range(i+1, len(elems)):
            d = min_dist_pbc(apos[i], apos[j])
            if elems[i] == 'H' and elems[j] == 'H':
                min_HH = min(min_HH, d)
            elif elems[i] == 'H' or elems[j] == 'H':
                min_HX = min(min_HX, d)
            else:
                min_XX = min(min_XX, d)

    print(f"  Min H-H distance: {min_HH:.3f} A")
    print(f"  Min H-X distance: {min_HX:.3f} A")
    print(f"  Min X-X distance: {min_XX:.3f} A")
    print(f"\nSaved: dry_run_initial.xyz, dry_run_final.xyz, dry_run_movie.xyz")

def main():
    import argparse
    parser = argparse.ArgumentParser(description='NEB for H-transfer between pyrazine molecules')
    parser.add_argument('--hydro', type=str, required=True, help='Path to hydro-pyrazine XYZ')
    parser.add_argument('--dehydro', type=str, required=True, help='Path to dehydro-pyrazine XYZ')
    parser.add_argument('--mode', type=str, default='neb', choices=['single', 'neb', 'rigid_scan', 'relax_scan', 'dry_run'])
    parser.add_argument('--n_points', type=int, default=11, help='Number of points (for rigid_scan/relax_scan)')
    parser.add_argument('--Temperature', type=float, default=800, help='Fermi broadening temperature [K]')
    parser.add_argument('--MixingParameter', type=float, default=0.05, help='Broyden mixing parameter')
    parser.add_argument('--MaxScc', type=int, default=2000, help='Maximum SCC iterations')
    parser.add_argument('--SCCTolerance', type=float, default=1e-6, help='SCC tolerance')
    parser.add_argument('--fixed', type=str, default='', help='Extra fixed atom indices (0-based), e.g. "0,1,2"')
    parser.add_argument('--rotate_y', type=float, default=0.0, help='Rotate dehydro molecule around y axis by degrees (to avoid peripheral H clashes)')
    parser.add_argument('--use_NN_distance', action='store_true', help='Place molecules based on N-N distance (~3.03 A for H...N ~2 A) instead of bounding boxes')
    parser.add_argument('--Ly', type=float, default=25.0, help='Lattice constant along y')
    parser.add_argument('--Lx', type=float, default=20.0, help='Cell size x')
    parser.add_argument('--Lz', type=float, default=20.0, help='Cell size z')
    parser.add_argument('--L_Hb', type=float, default=2.0, help='H-bond distance (A)')
    args = parser.parse_args()

    global Ly, Lx, Lz, L_Hb
    Ly = args.Ly
    Lx = args.Lx
    Lz = args.Lz
    L_Hb = args.L_Hb

    if args.mode == 'single':
        print("Building molecular cell...")
        apos, elems, lvs, mol_id = build_molecular_cell(args.hydro, args.dehydro, Ly, Lx, Lz, L_Hb)
        print("\nRunning single-point DFTB+ calculation with molecular orbitals...")
        temp_dir = 'single_point_calc'
        E, forces, orbitals, cube_files = run_dftb_calculation(
            apos, elems, lvs, temp_dir, do_forces=True, do_orbitals=True,
            Temperature=args.Temperature, MixingParameter=args.MixingParameter, MaxScc=args.MaxScc, SCCTolerance=args.SCCTolerance,
            reuse_charges=False)
        if E is not None:
            print(f"Total Energy: {E:.6f} eV")
            if orbitals is not None:
                print(f"Found {len(orbitals)} orbital energies")
                print(f"HOMO: {orbitals[-1]:.6f} eV, LUMO: {orbitals[0]:.6f} eV")
            if cube_files:
                print(f"Generated cube files: {[os.path.basename(c) for c in cube_files]}")
    elif args.mode == 'rigid_scan':
        run_rigid_scan(args.hydro, args.dehydro, n_points=args.n_points,
                       Temperature=args.Temperature, MixingParameter=args.MixingParameter, MaxScc=args.MaxScc, SCCTolerance=args.SCCTolerance,
                       rotate_y_deg=args.rotate_y, use_NN_distance=args.use_NN_distance)
    elif args.mode == 'relax_scan':
        run_relax_scan(args.hydro, args.dehydro, n_points=args.n_points,
                       Temperature=args.Temperature, MixingParameter=args.MixingParameter, MaxScc=args.MaxScc, SCCTolerance=args.SCCTolerance,
                       fixed_extra=_parse_fixed_atoms(args.fixed), rotate_y_deg=args.rotate_y, use_NN_distance=args.use_NN_distance)
    elif args.mode == 'dry_run':
        run_dry_run(args.hydro, args.dehydro, rotate_y_deg=args.rotate_y, use_NN_distance=args.use_NN_distance, n_points=args.n_points)
    else:
        run_neb_calculation(args.hydro, args.dehydro)

if __name__ == '__main__':
    main()
