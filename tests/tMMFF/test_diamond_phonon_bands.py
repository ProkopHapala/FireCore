#!/usr/bin/env python3
"""
Phonon band structure for diamond primitive cell using supercell force constants.

Method:
1. Build a NxNxN supercell from the primitive cell
2. Compute the full supercell Hessian (no PBC to avoid wrap-around)
3. Extract real-space force-constant blocks K(0,L) for the central cell
4. For each k-point along a path, build D(k) via Bloch phase sum
5. Diagonalize D(k) to get phonon frequencies omega(k)
"""
import os, sys, tempfile, argparse
import numpy as np

sys.path.append("../../")
from pyBall import MMFF

# ========================================================================
# CLI
# ========================================================================
parser = argparse.ArgumentParser(description='Diamond phonon band structure')
parser.add_argument('--unit', choices=['cm-1', 'THz'], default='cm-1', help='Frequency unit (default: cm-1)')
args = parser.parse_args()

# ========================================================================
# CONFIG
# ========================================================================
xyz_path = os.path.join(os.path.dirname(__file__), "..", "..", "cpp/common_resources/crystals/diamond_primitive.xyz")
SUPER_N = 3          # NxNxN supercell (odd required for central cell)
dx = 1e-4            # finite-difference displacement (Bohr)
mass_C = 12.0107     # amu
# Conversion: sqrt(Hartree/(amu*Bohr^2)) -> cm^-1
freq_conv_cm = 16.25
# 1 THz = 33.356 cm^-1
freq_conv = freq_conv_cm / 33.356 if args.unit == 'THz' else freq_conv_cm
unit_label = args.unit

# ========================================================================
# 1. Read primitive cell
# ========================================================================
def read_primitive_xyz(path):
    with open(path) as f:
        natoms = int(f.readline().strip())
        lvs_line = f.readline().strip()
        parts = lvs_line.split()
        # Format: "lvs  a1x a1y a1z  a2x a2y a2z  a3x a3y a3z"
        vals = list(map(float, parts[1:]))
        lvec = np.array(vals).reshape(3, 3)
        symbols = []
        pos = []
        for i in range(natoms):
            line = f.readline().strip()
            tok = line.split()
            symbols.append(tok[0])
            pos.append([float(tok[1]), float(tok[2]), float(tok[3])])
    return np.array(pos), lvec, symbols

prim_pos, prim_lvec, prim_sym = read_primitive_xyz(xyz_path)
n_prim = len(prim_pos)
print(f"Primitive cell: {n_prim} atoms")
print(f"Lattice vectors (Bohr):\n{prim_lvec}")

# ========================================================================
# 2. Build supercell (odd N => central cell at index (Nc,Nc,Nc))
# ========================================================================
assert SUPER_N % 2 == 1, "SUPER_N must be odd for clean central cell"
Nc = SUPER_N // 2  # center index

sc_pos = []
sc_cell = []       # record (ix,iy,iz) for each atom
sc_ia = []         # record primitive atom index for each atom

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
print(f"Supercell size: {SUPER_N}x{SUPER_N}x{SUPER_N} = {n_sc} atoms")

# Write temporary XYZ file for supercell
# Supercell lattice vectors (no PBC)
sc_lvec = SUPER_N * prim_lvec
with tempfile.NamedTemporaryFile(mode='w', suffix='.xyz', delete=False) as f:
    tmp_xyz = f.name
    f.write(f"{n_sc}\n")
    f.write(f"lvs  {sc_lvec[0,0]:.6f} {sc_lvec[0,1]:.6f} {sc_lvec[0,2]:.6f}  {sc_lvec[1,0]:.6f} {sc_lvec[1,1]:.6f} {sc_lvec[1,2]:.6f}  {sc_lvec[2,0]:.6f} {sc_lvec[2,1]:.6f} {sc_lvec[2,2]:.6f}\n")
    for sym, p in zip([prim_sym[sc_ia[i]] for i in range(n_sc)], sc_pos):
        f.write(f"{sym}  {p[0]:.6f}  {p[1]:.6f}  {p[2]:.6f}\n")

# ========================================================================
# 3. Compute supercell Hessian (NO PBC to avoid wrap-around!)
# ========================================================================
print("\nInitializing MMFF for supercell (no PBC)...")
MMFF.init(
    xyz_name=tmp_xyz,
    nPBC=(0, 0, 0),   # cluster/supercell without PBC
    bEpairs=False,
    bMMFF=True,
)
print("Supercell Hessian calculation...")
inds = np.arange(n_sc, dtype=np.int32)
H_sc = MMFF.getHessian3Nx3N(inds, dx=dx)
H_sc = 0.5 * (H_sc + H_sc.T)

if np.isnan(H_sc).any() or np.isinf(H_sc).any():
    raise ValueError("NaN/Inf in supercell Hessian")

print(f"Supercell Hessian shape: {H_sc.shape}")
print(f"Supercell Hessian norm: {np.linalg.norm(H_sc):.6e}")

# ========================================================================
# 4. Extract K(0, L) blocks from supercell Hessian
# ========================================================================
# We only extract blocks where the displaced atom is in the CENTRAL cell.
# Find atoms in central cell (cell index = (0,0,0))
central_atoms = [i for i, c in enumerate(sc_cell) if c == (0, 0, 0)]
print(f"Central cell atoms: {central_atoms}")

# Build map: (cell_ix, cell_iy, cell_iz, prim_ia) -> supercell atom index
atom_map = {}
for i, (c, ia) in enumerate(zip(sc_cell, sc_ia)):
    atom_map[(c[0], c[1], c[2], ia)] = i

# Extract Phi[i,j,R] where i is in central cell, j is in cell R
# Phi has shape (n_prim, n_prim, 3, 3) for each R vector
# We store as dict: R_tuple -> (n_prim, n_prim, 3, 3) array
Phi_blocks = {}

for i0 in central_atoms:
    ia_i = sc_ia[i0]        # primitive atom index of i
    for j in range(n_sc):
        R = sc_cell[j]      # cell of atom j relative to central
        ia_j = sc_ia[j]     # primitive atom index of j
        Rkey = R
        if Rkey not in Phi_blocks:
            Phi_blocks[Rkey] = np.zeros((n_prim, n_prim, 3, 3))
        # Extract 3x3 block from Hessian
        block = H_sc[i0*3:(i0+1)*3, j*3:(j+1)*3]
        Phi_blocks[Rkey][ia_i, ia_j] = block

print(f"Extracted {len(Phi_blocks)} unique R vectors")

# Clean up temp file
os.unlink(tmp_xyz)

# ========================================================================
# 5. Define k-path for diamond (FCC primitive cell)
# ========================================================================
# Reciprocal lattice vectors for primitive cell
# b_i satisfy a_j · b_i = 2π δ_ji  =>  b_i = 2π (a_j × a_k) / (a1 · (a2 × a3))
vol = np.dot(prim_lvec[0], np.cross(prim_lvec[1], prim_lvec[2]))
recip = np.zeros((3, 3))
recip[0] = 2 * np.pi * np.cross(prim_lvec[1], prim_lvec[2]) / vol
recip[1] = 2 * np.pi * np.cross(prim_lvec[2], prim_lvec[0]) / vol
recip[2] = 2 * np.pi * np.cross(prim_lvec[0], prim_lvec[1]) / vol
print(f"Reciprocal lattice vectors:\n{recip}")

# High-symmetry points in Cartesian (units of 2π/a for conventional cell, a=3.567A)
# For primitive FCC, use standard points in reciprocal space
# Here we express k-points directly in Cartesian (inverse Bohr)
a0 = np.linalg.norm(prim_lvec[0] + prim_lvec[1] - prim_lvec[2])  # conventional cubic edge
# Actually easier: define points in units of recip vectors, then convert

def k_from_frac(kf):
    """Convert fractional reciprocal coords to Cartesian."""
    return kf[0] * recip[0] + kf[1] * recip[1] + kf[2] * recip[2]

# Standard FCC BZ high-symmetry points in fractional reciprocal coordinates
# (i.e. coefficients of b1, b2, b3)
pts_frac = {
    'Γ': np.array([0.0, 0.0, 0.0]),
    'X': np.array([0.5, 0.5, 0.0]),
    'W': np.array([0.5, 0.25, 0.75]),
    'K': np.array([0.375, 0.375, 0.75]),
    'L': np.array([0.5, 0.5, 0.5]),
    'U': np.array([0.25, 0.25, 0.75]),
}

# k-path segments (standard diamond dispersion path)
path_segments = [
    ('Γ', 'X'),
    ('X', 'W'),
    ('W', 'K'),
    ('K', 'Γ'),
    ('Γ', 'L'),
    ('L', 'U'),
    ('U', 'W'),
]

npts_per_seg = 40

kpts = []
k_labels = []
k_frac = []

for seg in path_segments:
    p1, p2 = seg
    k1 = k_from_frac(pts_frac[p1])
    k2 = k_from_frac(pts_frac[p2])
    for i in range(npts_per_seg):
        t = i / (npts_per_seg - 1) if npts_per_seg > 1 else 0
        kpts.append(k1 + t * (k2 - k1))
        k_frac.append(pts_frac[p1] + t * (pts_frac[p2] - pts_frac[p1]))
    k_labels.append((len(kpts) - npts_per_seg, p1))
k_labels.append((len(kpts) - 1, p2))

kpts = np.array(kpts)
nk = len(kpts)

# ========================================================================
# 6. Build D(k) and diagonalize for each k-point
# ========================================================================
dim = 3 * n_prim
freqs = np.zeros((nk, dim))

masses = np.full(n_prim, mass_C)

print(f"\nComputing phonon frequencies at {nk} k-points...")
for ik, k in enumerate(kpts):
    Dk = np.zeros((dim, dim), dtype=complex)
    for Rkey, Phi_R in Phi_blocks.items():
        R_cart = Rkey[0] * prim_lvec[0] + Rkey[1] * prim_lvec[1] + Rkey[2] * prim_lvec[2]
        phase = np.exp(1j * np.dot(k, R_cart))
        for i in range(n_prim):
            for j in range(n_prim):
                block = Phi_R[i, j] * phase
                Dk[i*3:(i+1)*3, j*3:(j+1)*3] += block / np.sqrt(masses[i] * masses[j])

    # Enforce Hermitian (numerical noise may break it slightly)
    Dk = 0.5 * (Dk + Dk.conj().T)
    eigvals = np.linalg.eigvalsh(Dk)
    freqs[ik, :] = np.sign(eigvals) * np.sqrt(np.abs(eigvals)) * freq_conv

print("Done.")

# ========================================================================
# 7. Plot phonon band structure
# ========================================================================
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

kdist = np.zeros(nk)
for i in range(1, nk):
    dk = np.linalg.norm(kpts[i] - kpts[i-1])
    kdist[i] = kdist[i-1] + dk

fig, ax = plt.subplots(figsize=(8, 6))
for imode in range(dim):
    ax.plot(kdist, freqs[:, imode], 'b-', lw=0.8)

# Vertical lines at segment boundaries
for pos, lab in k_labels:
    ax.axvline(x=kdist[pos], color='k', lw=0.5)
ax.set_xticks([kdist[pos] for pos, _ in k_labels])
ax.set_xticklabels([lab for _, lab in k_labels])

ax.set_ylabel(f'Frequency ({unit_label})')
ax.set_xlabel('k-path')
ax.set_title(f'Diamond phonon dispersion ({unit_label}, supercell {SUPER_N}x{SUPER_N}x{SUPER_N})')
ax.set_xlim(kdist[0], kdist[-1])
ax.set_ylim(bottom=0)

plt.tight_layout()
outpng = 'diamond_phonon_bands.png'
plt.savefig(outpng, dpi=150)
print(f"\nBand structure saved to: {outpng}")
outpng_unit = f'diamond_phonon_bands_{args.unit}.png'
os.rename(outpng, outpng_unit)
outpng = outpng_unit
print(f"Renamed to: {outpng}")

# Also save data
np.savez(f'diamond_phonon_bands_{args.unit}.npz',
         kdist=kdist, freqs=freqs, kpts=kpts,
         prim_lvec=prim_lvec, recip=recip)
print(f"Data saved to: diamond_phonon_bands_{args.unit}.npz")

# Print Gamma point frequencies (first point in our path)
sorted_gamma = np.sort(np.abs(freqs[0]))
print(f"\nΓ-point frequencies ({unit_label}):")
for f in sorted_gamma:
    print(f"  {f:10.3f}")
