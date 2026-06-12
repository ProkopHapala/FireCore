#!/usr/bin/env python3
"""
Phonon band structure for diamond primitive cell using supercell force constants.

Workflow (phonopy-compatible):
1. Build NxNxN supercell (odd N => central primitive cell at R=0)
2. Displace only central-cell atoms; read forces on all supercell atoms (getPhononPhiBlocks)
3. Extract real-space Phi(0,R) blocks
4. Bloch sum D(k) and diagonalize

Use --pbc for periodic supercell (recommended: correct boundary coordination).
Cluster mode (no PBC) gives correct Gamma optical modes but acoustic branches can go
imaginary at finite k due to under-coordinated surface atoms in 3x3x3.
"""
import os, sys, tempfile, argparse
import numpy as np

sys.path.append("../../")
from pyBall import MMFF

MMFF.setVerbosity(0,0)

parser = argparse.ArgumentParser(description='Diamond phonon band structure')
parser.add_argument('--unit', choices=['cm-1', 'THz'], default='cm-1', help='Frequency unit (default: cm-1)')
parser.add_argument('--super-n', type=int, default=3, help='NxNxN supercell (odd, default: 3)')
parser.add_argument('--pbc', action='store_true', help='PBC on supercell (removes negatives but optical modes ~20x too soft — do not use)')
parser.add_argument('--asr', action='store_true', help='Phonopy FC symmetrization (fixes cluster acoustic negatives, needs phonopy)')
parser.add_argument('--plot', action='store_true', help='Generate PNG plot (may crash under ASan due to matplotlib/ft2font)')
parser.add_argument('--plot-only', type=str, default=None, help='Plot from cached .npz (skip MMFF FD); path to diamond_phonon_bands_*.npz')
parser.add_argument('--dx', type=float, default=1e-4, help='Finite-difference displacement (Bohr)')
parser.add_argument('--bondsOnly', action='store_true', help='Disable angles, PiSigma, and PiPiI terms (bond-only)')
parser.add_argument('--uff', action='store_true', help='Use UFF forcefield instead of MMFF')
parser.add_argument('--no-bonds', action='store_true', help='Disable bond terms')
parser.add_argument('--no-angles', action='store_true', help='Disable angle terms (MMFF) or angle terms (UFF)')
parser.add_argument('--no-pisigma', action='store_true', help='Disable PiSigma terms (MMFF only)')
parser.add_argument('--no-pipii', action='store_true', help='Disable PiPiI terms (MMFF only)')
parser.add_argument('--no-dihedrals', action='store_true', help='Disable dihedral terms (UFF only)')
parser.add_argument('--no-inversions', action='store_true', help='Disable inversion terms (UFF only)')
parser.add_argument('--fix-r0', action='store_true', help='Override UFF C_3-C_3 r0 to match actual diamond lattice constant')
parser.add_argument('--scale-bond', type=float, default=None, help='Scale MMFF bond stiffness by factor (e.g., 1.25 for +25%)')
parser.add_argument('--scale-angle', type=float, default=None, help='Scale MMFF angle stiffness by factor (e.g., 1.25 for +25%)')
args = parser.parse_args()

SUPER_N = args.super_n
dx = args.dx
mass_C = 12.0107
unit_label = args.unit

if args.plot_only is not None:
    data = np.load(args.plot_only)
    kdist = data['kdist']
    freqs = data['freqs']
    super_n = int(data['super_n']) if 'super_n' in data.files else SUPER_N
    pbc = bool(data['pbc']) if 'pbc' in data.files else False
    asr = bool(data['asr']) if 'asr' in data.files else False
    mode_tag = 'PBC' if pbc else ('asr' if asr else 'cluster')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8, 6))
    for imode in range(freqs.shape[1]):
        ax.plot(kdist, freqs[:, imode], 'b-', lw=0.8)
    ax.set_ylabel(f'Frequency ({unit_label})')
    ax.set_xlabel('k-path')
    ax.set_title(f'Diamond phonons ({unit_label}, {super_n}x{super_n}x{super_n} {mode_tag})')
    ax.set_xlim(kdist[0], kdist[-1])
    plt.tight_layout()
    outpng = f'diamond_phonon_bands_{unit_label}_{mode_tag}_plotonly.png'
    plt.savefig(outpng, dpi=150)
    print(f"Band structure saved to: {outpng}")
    raise SystemExit(0)
assert SUPER_N % 2 == 1, "SUPER_N must be odd for clean central cell"
Nc = SUPER_N // 2
RMAX = Nc  # extract Phi for |R| <= Nc (full supercell range)


def read_primitive_xyz(path):
    with open(path) as f:
        natoms = int(f.readline().strip())
        parts = f.readline().strip().split()
        lvec = np.array(list(map(float, parts[1:]))).reshape(3, 3)
        symbols, pos = [], []
        for i in range(natoms):
            tok = f.readline().strip().split()
            symbols.append(tok[0])
            pos.append([float(tok[1]), float(tok[2]), float(tok[3])])
    return np.array(pos), lvec, symbols


def extract_phi_blocks(phi, sc_cell, sc_ia, inds_disp, n_prim, rmax=None):
    """Build Phi(R) from rectangular phi matrix (n_sc*3 x n_disp*3)."""
    if rmax is None:
        rmax = max(max(abs(c[0]), abs(c[1]), abs(c[2])) for c in sc_cell)
    Phi_blocks = {}
    for p_idx, ip in enumerate(inds_disp):
        ia_i = sc_ia[ip]
        n_sc = len(sc_cell)
        for j in range(n_sc):
            R = sc_cell[j]
            if max(abs(R[0]), abs(R[1]), abs(R[2])) > rmax:
                continue
            Rkey = tuple(R)
            ia_j = sc_ia[j]
            if Rkey not in Phi_blocks:
                Phi_blocks[Rkey] = np.zeros((n_prim, n_prim, 3, 3))
            Phi_blocks[Rkey][ia_i, ia_j] = phi[j*3:(j+1)*3, p_idx*3:(p_idx+1)*3]
    return Phi_blocks


def reciprocal_lattice(cell):
    vol = np.dot(cell[0], np.cross(cell[1], cell[2]))
    recip = np.zeros((3, 3))
    recip[0] = 2 * np.pi * np.cross(cell[1], cell[2]) / vol
    recip[1] = 2 * np.pi * np.cross(cell[2], cell[0]) / vol
    recip[2] = 2 * np.pi * np.cross(cell[0], cell[1]) / vol
    return recip


def apply_phonopy_asr(phi, sc_pos, sc_cell, sc_ia, central_atoms, n_prim, prim_lvec, prim_pos, level=10):
    """Symmetrize supercell FC via phonopy (removes cluster-boundary acoustic artifacts)."""
    from phonopy.harmonic.force_constants import symmetrize_force_constants
    n_sc = len(sc_cell)
    sc_ref = []
    for iz in range(SUPER_N):
        for iy in range(SUPER_N):
            for ix in range(SUPER_N):
                R = ix * prim_lvec[0] + iy * prim_lvec[1] + iz * prim_lvec[2]
                for ia, p in enumerate(prim_pos):
                    sc_ref.append(p + R)
    sc_ref = np.array(sc_ref)
    perm = [int(np.argmin(np.linalg.norm(sc_ref - p, axis=1))) for p in sc_pos]
    fc = np.zeros((n_sc, n_sc, 3, 3))
    for p_idx, i0 in enumerate(central_atoms):
        for j in range(n_sc):
            fc[perm[i0], perm[j]] = phi[j*3:(j+1)*3, p_idx*3:(p_idx+1)*3]
    symmetrize_force_constants(fc, level=level)
    Phi_blocks = {}
    for p_idx, i0 in enumerate(central_atoms):
        for j in range(n_sc):
            Rkey = tuple(sc_cell[j])
            if max(abs(Rkey[0]), abs(Rkey[1]), abs(Rkey[2])) > RMAX:
                continue
            Phi_blocks.setdefault(Rkey, np.zeros((n_prim, n_prim, 3, 3)))
            Phi_blocks[Rkey][sc_ia[i0], sc_ia[j]] = fc[perm[i0], perm[j]]
    return Phi_blocks


def solve_bands(Phi_blocks, prim_lvec, masses, kpts_frac, unit_label, convention='signed'):
    n_prim = len(masses)
    dim = 3 * n_prim
    recip = reciprocal_lattice(prim_lvec)
    nk = len(kpts_frac)
    freqs = np.zeros((nk, dim))
    for ik, qf in enumerate(kpts_frac):
        k = qf[0] * recip[0] + qf[1] * recip[1] + qf[2] * recip[2]
        Dk = np.zeros((dim, dim), dtype=complex)
        for Rkey, Phi_R in Phi_blocks.items():
            R_cart = Rkey[0] * prim_lvec[0] + Rkey[1] * prim_lvec[1] + Rkey[2] * prim_lvec[2]
            phase = np.exp(1j * np.dot(k, R_cart))
            for i in range(n_prim):
                for j in range(n_prim):
                    block = Phi_R[i, j] * phase
                    Dk[i*3:(i+1)*3, j*3:(j+1)*3] += block / np.sqrt(masses[i] * masses[j])
        Dk = 0.5 * (Dk + Dk.conj().T)
        eigvals = np.linalg.eigvalsh(Dk)  # [eV/Å^2/amu]
        fTHz = thz_from_eig_massweighted(eigvals)
        if unit_label == 'cm-1':
            freqs[ik, :] = fTHz * 33.356
        else:
            freqs[ik, :] = fTHz
    return freqs

def thz_from_eig_massweighted(eig_eVA2_per_amu):
    eV=1.602176634e-19
    Ang=1e-10
    amu=1.66053906660e-27
    conv = (eV/(Ang*Ang))/amu
    w = np.sign(eig_eVA2_per_amu)*np.sqrt(np.abs(eig_eVA2_per_amu)*conv)
    f = w/(2*np.pi)*1e-12
    return f

def gamma_reports(Phi_blocks, masses_amu):
    n_prim=len(masses_amu); dim=3*n_prim
    K = np.zeros((dim,dim),dtype=float)
    for _,Phi_R in Phi_blocks.items():
        for i in range(n_prim):
            for j in range(n_prim):
                K[i*3:(i+1)*3, j*3:(j+1)*3] += Phi_R[i,j]
    K = 0.5*(K+K.T)
    lam = np.sort(np.linalg.eigvalsh(K))
    w = np.repeat(1.0/np.sqrt(np.array(masses_amu,dtype=float)),3)
    D = (K*w[None,:])*w[:,None]
    mu = np.sort(np.linalg.eigvalsh(0.5*(D+D.T)))
    fTHz = thz_from_eig_massweighted(mu)
    print("\nΓ-point reporting (always):")
    print("stiffness eig λ  [eV/Å^2]     min..max:", lam[0], lam[-1])
    print("mass-weighted eig μ [eV/Å^2/amu] min..max:", mu[0], mu[-1])
    print("freq             [THz]         sorted:", fTHz)


xyz_path = os.path.join(os.path.dirname(__file__), "..", "..", "cpp/common_resources/crystals/diamond_primitive.xyz")
prim_pos, prim_lvec, prim_sym = read_primitive_xyz(xyz_path)
n_prim = len(prim_pos)
print(f"Primitive cell: {n_prim} atoms")

sc_pos, sc_cell, sc_ia = [], [], []
for iz in range(SUPER_N):
    for iy in range(SUPER_N):
        for ix in range(SUPER_N):
            R = ix * prim_lvec[0] + iy * prim_lvec[1] + iz * prim_lvec[2]
            for ia, p in enumerate(prim_pos):
                sc_pos.append(p + R)
                sc_cell.append((ix - Nc, iy - Nc, iz - Nc))
                sc_ia.append(ia)
sc_pos = np.array(sc_pos)
n_sc = len(sc_pos)
sc_lvec = SUPER_N * prim_lvec
print(f"Supercell: {SUPER_N}x{SUPER_N}x{SUPER_N} = {n_sc} atoms, PBC={args.pbc}")

with tempfile.NamedTemporaryFile(mode='w', suffix='.xyz', delete=False) as f:
    tmp_xyz = f.name
    f.write(f"{n_sc}\n")
    f.write(f"lvs  {sc_lvec[0,0]:.6f} {sc_lvec[0,1]:.6f} {sc_lvec[0,2]:.6f}  {sc_lvec[1,0]:.6f} {sc_lvec[1,1]:.6f} {sc_lvec[1,2]:.6f}  {sc_lvec[2,0]:.6f} {sc_lvec[2,1]:.6f} {sc_lvec[2,2]:.6f}\n")
    for sym, p in zip([prim_sym[sc_ia[i]] for i in range(n_sc)], sc_pos):
        f.write(f"{sym}  {p[0]:.6f}  {p[1]:.6f}  {p[2]:.6f}\n")

nPBC = (1, 1, 1) if args.pbc else (0, 0, 0)
print(f"\nInitializing MMFF (nPBC={nPBC})...")
MMFF.init(xyz_name=tmp_xyz, nPBC=nPBC, bEpairs=False, bMMFF=not args.uff, bUFF=args.uff)

if args.uff:
    # UFF mode
    do_bond = -1 if args.no_bonds else 1
    do_angle = -1 if args.no_angles else 1
    do_dihedral = -1 if args.no_dihedrals else 1
    do_inversion = -1 if args.no_inversions else 1
    MMFF.setSwitchesUFF(DoBond=do_bond, DoAngle=do_angle, DoDihedral=do_dihedral, DoInversion=do_inversion, DoAssemble=1, SubtractBondNonBond=-1, ClampNonBonded=-1)
    print(f"UFF mode: bonds={do_bond}, angles={do_angle}, dihedrals={do_dihedral}, inversions={do_inversion}")
    if args.fix_r0:
        r0_actual = np.linalg.norm(prim_pos[1] - prim_pos[0])
        print(f"Overriding C_3-C_3 r0: UFF default ~1.514 => actual diamond bond {r0_actual:.6f} Ang")
        MMFF.setBondParamsByType('C_3', 'C_3', r0=r0_actual)
else:
    # MMFF mode
    if args.bondsOnly:
        MMFF.setSwitches(PBC=-1, NonBonded=-1, Angles=-1, PiSigma=-1, PiPiI=-1)
        print("Bond-only mode: Angles, PiSigma, PiPiI OFF")
    else:
        do_angles = -1 if args.no_angles else 1
        do_pisigma = -1 if args.no_pisigma else 1
        do_pipii = -1 if args.no_pipii else 1
        if args.pbc:
            MMFF.setSwitches(NonBonded=-1, Angles=do_angles, PiSigma=do_pisigma, PiPiI=do_pipii)
            print(f"Explicit supercell phonon mode: PBC kept ON, nonbonded OFF, angles={do_angles}, PiSigma={do_pisigma}, PiPiI={do_pipii}")
        else:
            MMFF.setSwitches(PBC=-1, NonBonded=-1, Angles=do_angles, PiSigma=do_pisigma, PiPiI=do_pipii)
            print(f"Cluster mode: PBC OFF, nonbonded OFF, angles={do_angles}, PiSigma={do_pisigma}, PiPiI={do_pipii}")
    
    # Apply MMFF parameter scaling if requested (using generalized API)
    if args.scale_bond is not None:
        # Scale all bonds by reading current values and multiplying
        MMFF.getBuffs()
        current_k = MMFF.bKs[MMFF.bKs != 0].mean()  # average non-zero bond stiffness
        MMFF.setBondParamsByType('C_3', 'C_3', k=current_k * args.scale_bond, forcefield='MMFF')
    if args.scale_angle is not None:
        # Scale all angles by reading current values and multiplying
        MMFF.getBuffs()
        current_k = MMFF.apars[:, 1].mean()  # average Kss
        MMFF.setAngleParamsByType('C_3', 'C_3', 'C_3', k=current_k * args.scale_angle, forcefield='MMFF')

inds_total = np.arange(n_sc, dtype=np.int32)
central_atoms = [i for i, c in enumerate(sc_cell) if c == (0, 0, 0)]
inds_disp = np.array(central_atoms, dtype=np.int32)
print(f"Central cell atoms: {list(central_atoms)}")
print(f"Computing Phi blocks (displace {len(inds_disp)} atoms, read forces on {n_sc})...")
phi = MMFF.getPhononPhiBlocks(inds_total, inds_disp, dx=dx)
if np.isnan(phi).any() or np.isinf(phi).any():
    raise ValueError("NaN/Inf in Phi matrix")
print(f"Phi matrix shape: {phi.shape}, norm: {np.linalg.norm(phi):.6e}")

Phi_blocks = extract_phi_blocks(phi, sc_cell, sc_ia, central_atoms, n_prim, rmax=RMAX)
print(f"Extracted {len(Phi_blocks)} Phi(0,R) blocks (|R|_inf <= {RMAX})")
if args.asr:
    try:
        Phi_blocks = apply_phonopy_asr(phi, sc_pos, sc_cell, sc_ia, central_atoms, n_prim, prim_lvec, prim_pos)
        print(f"Applied phonopy FC symmetrization (level=10)")
    except ImportError:
        raise ImportError("--asr requires phonopy (pip install phonopy)")

# Translational ASR diagnostic: sum_{j,R} Phi[i,j,R] should be ~0
trans = np.zeros((n_prim, 3, 3))
for i in range(n_prim):
    for Rkey, Pr in Phi_blocks.items():
        for j in range(n_prim):
            trans[i] += Pr[i, j]
print(f"Translational ASR violation norm: {np.linalg.norm(trans):.3e}")

os.unlink(tmp_xyz)

recip = reciprocal_lattice(prim_lvec)
print(f"Reciprocal lattice vectors:\n{recip}")

pts_frac = {
    'Γ': np.array([0.0, 0.0, 0.0]),
    'X': np.array([0.5, 0.5, 0.0]),
    'W': np.array([0.5, 0.25, 0.75]),
    'K': np.array([0.375, 0.375, 0.75]),
    'L': np.array([0.5, 0.5, 0.5]),
    'U': np.array([0.25, 0.25, 0.75]),
}
path_segments = [('Γ', 'X'), ('X', 'W'), ('W', 'K'), ('K', 'Γ'), ('Γ', 'L'), ('L', 'U'), ('U', 'W')]
npts_per_seg = 40
kpts, k_labels, k_frac = [], [], []
for seg in path_segments:
    p1, p2 = seg
    k1 = pts_frac[p1][0] * recip[0] + pts_frac[p1][1] * recip[1] + pts_frac[p1][2] * recip[2]
    k2 = pts_frac[p2][0] * recip[0] + pts_frac[p2][1] * recip[1] + pts_frac[p2][2] * recip[2]
    for i in range(npts_per_seg):
        t = i / (npts_per_seg - 1) if npts_per_seg > 1 else 0
        kpts.append(k1 + t * (k2 - k1))
        k_frac.append(pts_frac[p1] + t * (pts_frac[p2] - pts_frac[p1]))
    k_labels.append((len(kpts) - npts_per_seg, p1))
k_labels.append((len(kpts) - 1, p2))
kpts = np.array(kpts)
k_frac = np.array(k_frac)
masses = np.full(n_prim, mass_C)

print(f"\nComputing phonon frequencies at {len(kpts)} k-points...")
freqs = solve_bands(Phi_blocks, prim_lvec, masses, k_frac, unit_label, convention='signed')
gamma_reports(Phi_blocks, masses)
neg_count = int(np.sum(freqs < -0.01))
fs = np.sort(freqs, axis=1)
print(f"Negative frequency count (signed, < -0.01 {unit_label}): {neg_count}")
print(f"Lowest 3 bands min: {fs[:, :3].min(axis=0)}")
print(f"Γ-point (sorted): {np.sort(freqs[0])}")

if neg_count > 0 and not args.pbc and not args.asr:
    print("\nWARNING: cluster 3x3x3 (no PBC) gives imaginary acoustic branches at finite k (surface under-coordination).")
    print("         Re-run with --asr (phonopy FC symmetrization) to fix.")

kdist = np.zeros(len(kpts))
for i in range(1, len(kpts)):
    kdist[i] = kdist[i-1] + np.linalg.norm(kpts[i] - kpts[i-1])

mode_tag = 'PBC' if args.pbc else ('asr' if args.asr else 'cluster')
if args.uff:
    mode_tag += '_uff'
if args.bondsOnly:
    mode_tag += '_bondsonly'
if args.no_bonds:
    mode_tag += '_nobonds'
if args.no_angles:
    mode_tag += '_noangles'
if args.no_pisigma:
    mode_tag += '_nopisigma'
if args.no_pipii:
    mode_tag += '_nopipii'
if args.no_dihedrals:
    mode_tag += '_nodihedrals'
if args.no_inversions:
    mode_tag += '_noinversions'
if args.plot:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8, 6))
    for imode in range(3 * n_prim):
        ax.plot(kdist, freqs[:, imode], 'b-', lw=0.8)
    for pos, lab in k_labels:
        ax.axvline(x=kdist[pos], color='k', lw=0.5)
    ax.set_xticks([kdist[pos] for pos, _ in k_labels])
    ax.set_xticklabels([lab for _, lab in k_labels])
    ax.set_ylabel(f'Frequency ({unit_label})')
    ax.set_xlabel('k-path')
    ax.set_title(f'Diamond phonons ({unit_label}, {SUPER_N}x{SUPER_N}x{SUPER_N} {mode_tag})')
    ax.set_xlim(kdist[0], kdist[-1])
    plt.tight_layout()
    outpng = f'diamond_phonon_bands_{args.unit}_{mode_tag}.png'
    plt.savefig(outpng, dpi=150)
    print(f"\nBand structure saved to: {outpng}")

np.savez(f'diamond_phonon_bands_{args.unit}_{mode_tag}.npz',
         kdist=kdist, freqs=freqs, kpts=kpts, k_frac=k_frac,
         prim_lvec=prim_lvec, recip=recip, super_n=SUPER_N, pbc=args.pbc, asr=args.asr)
print(f"Data saved to: diamond_phonon_bands_{args.unit}_{mode_tag}.npz")
