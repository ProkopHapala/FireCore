"""
Ewald2D  —  2D Fourier representation of the electrostatic potential
         above / inside a periodic ionic slab.

The potential is expressed as a Cartesian product of in-plane plane waves
cos(G·ρ) and out-of-plane real exponentials exp(−|G|z):

    φ(ρ,z) = φ₀(z) + Σ_{G≠0} (2π/(A|G|)) Σ_i q_i e^{iG·(ρ−ρ_i)} e^{−|G||z−z_i|}

See:  doc/Topics/OnSurfaceAssembly/Ewald_2D.md  for the full derivation.

Units:  Å for lengths, elementary charge for q  →  potential in e/Å.
        Multiply by 14.4 to get eV / e.
"""

import numpy as np

# ============================================================
#  Reciprocal lattice
# ============================================================

def make_reciprocal_2d(a, b):
    """Compute area and 2D reciprocal vectors from real-space lattice vectors.

    Parameters
    ----------
    a, b : array-like (2,)
        In-plane lattice vectors.

    Returns
    -------
    area : float
        |a × b|
    b1, b2 : ndarray (2,)
        Reciprocal vectors satisfying  a_i · b_j = 2π δ_ij.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    area = abs(a[0]*b[1] - a[1]*b[0])
    assert area > 1e-14, f"Degenerate lattice: area={area}"
    b1 = 2*np.pi/area * np.array([ b[1], -b[0]])
    b2 = 2*np.pi/area * np.array([-a[1],  a[0]])
    return area, b1, b2


def generate_G_vectors(b1, b2, n_harm):
    """Generate G = h*b1 + k*b2 for |h|,|k| ≤ n_harm, excluding G=0.

    Returns
    -------
    hs, ks : ndarray (N,)   integer indices
    Gx, Gy : ndarray (N,)   Cartesian components
    Gn     : ndarray (N,)   |G|
    """
    hs, ks = [], []
    for h in range(-n_harm, n_harm+1):
        for k in range(-n_harm, n_harm+1):
            if h == 0 and k == 0: continue
            hs.append(h); ks.append(k)
    hs = np.array(hs); ks = np.array(ks)
    Gx = hs*b1[0] + ks*b2[0]
    Gy = hs*b1[1] + ks*b2[1]
    Gn = np.hypot(Gx, Gy)
    return hs, ks, Gx, Gy, Gn

# ============================================================
#  Structure factors / coefficients
# ============================================================

def compute_C_G(Gx, Gy, Gn, rx, ry, rz, q, area):
    """Vacuum-side coefficients (centrosymmetric / cosine basis).

    C_G = (4π/(A|G|)) Σ_i q_i cos(G·ρ_i) exp(|G| z_i)

    Returns ndarray (N_G,).
    """
    Gdotr = Gx[:,None]*rx[None,:] + Gy[:,None]*ry[None,:]   # (N_G, N_ions)
    prefactor = 4*np.pi / (area * Gn)                        # (N_G,)
    C = prefactor * np.sum(q[None,:] * np.cos(Gdotr) * np.exp(Gn[:,None]*rz[None,:]), axis=1)
    return C


def compute_w_per_ion(Gx, Gy, Gn, rx, ry, q, area):
    """Per-ion complex weights for the full (interior) potential.

    w[g,i] = (2π/(A|G|)) q_i exp(−iG·ρ_i)

    Returns complex ndarray (N_G, N_ions).
    """
    Gdotr = Gx[:,None]*rx[None,:] + Gy[:,None]*ry[None,:]
    prefactor = 2*np.pi / (area * Gn)  # (N_G,)
    w = prefactor[:,None] * q[None,:] * np.exp(-1j * Gdotr)
    return w

# ============================================================
#  Potential evaluation
# ============================================================

def eval_potential_vacuum(C, Gx, Gy, Gn, X, Y, z):
    """Potential on an XY grid at a single height z ≥ z_max.

    φ = φ₀ + Σ_G C_G cos(G·ρ) exp(−|G|z).   (φ₀ = 0 for neutral cell)

    Parameters
    ----------
    C : (N_G,)   vacuum-side coefficients
    X, Y : 2-D arrays (grid)
    z : float    height
    """
    phase = Gx[:,None,None]*X[None,:,:] + Gy[:,None,None]*Y[None,:,:]
    decay = np.exp(-Gn * z)  # (N_G,)
    phi = np.sum(C[:,None,None] * np.cos(phase) * decay[:,None,None], axis=0)
    return phi


def eval_potential_full_2d(w, Gx, Gy, Gn, rz, q, area, X, Y, Z):
    """Full potential on a 2-D grid (e.g. XZ plane).

    φ = φ₀(z) + Σ_G Σ_i w[g,i] e^{iG·ρ} e^{−|G||z−z_i|}

    X, Y, Z : 2-D arrays of identical shape.
    """
    N_ions = len(q)
    phi = np.zeros_like(X, dtype=float)
    # G=0 term
    for i in range(N_ions):
        phi += -(2*np.pi/area) * q[i] * np.abs(Z - rz[i])
    # G≠0 terms
    N_G = len(Gx)
    for g in range(N_G):
        phase = np.exp(1j * (Gx[g]*X + Gy[g]*Y))
        for i in range(N_ions):
            decay = np.exp(-Gn[g] * np.abs(Z - rz[i]))
            phi += np.real(w[g,i] * phase * decay)
    return phi


def eval_potential_full_1d(w, Gx, Gy, Gn, rz, q, area, x0, y0, z_arr):
    """Full potential along a vertical line at (x0, y0)."""
    N_ions = len(q)
    phi = np.zeros(len(z_arr), dtype=float)
    for i in range(N_ions):
        phi += -(2*np.pi/area) * q[i] * np.abs(z_arr - rz[i])
    N_G = len(Gx)
    for g in range(N_G):
        phase = np.exp(1j * (Gx[g]*x0 + Gy[g]*y0))
        for i in range(N_ions):
            decay = np.exp(-Gn[g] * np.abs(z_arr - rz[i]))
            phi += np.real(w[g,i] * phase * decay)
    return phi

# ============================================================
#  Brute-force reference (direct Coulomb sum)
# ============================================================

def eval_potential_brute(rx, ry, rz, q, a_vec, b_vec, x0, y0, z_arr, N_rep=20):
    """Direct Coulomb sum over periodic images, circular-shell ordered.

    φ(r) = Σ_{n,m} Σ_i q_i / |r − R_{nm} − r_i|

    Converges for charge-neutral unit cells when summed over expanding circles.
    """
    a_vec = np.asarray(a_vec, dtype=float)
    b_vec = np.asarray(b_vec, dtype=float)
    nm_pairs = []
    for n in range(-N_rep, N_rep+1):
        for m in range(-N_rep, N_rep+1):
            Rx = n*a_vec[0] + m*b_vec[0]
            Ry = n*a_vec[1] + m*b_vec[1]
            nm_pairs.append((Rx*Rx + Ry*Ry, n, m))
    nm_pairs.sort()
    phi = np.zeros(len(z_arr), dtype=float)
    for (_, n, m) in nm_pairs:
        Rx = n*a_vec[0] + m*b_vec[0]
        Ry = n*a_vec[1] + m*b_vec[1]
        for i in range(len(q)):
            dx = x0 - (Rx + rx[i])
            dy = y0 - (Ry + ry[i])
            dz = z_arr - rz[i]
            phi += q[i] / np.sqrt(dx*dx + dy*dy + dz*dz)
    return phi

# ============================================================
#  Charge-density Fourier reconstruction
# ============================================================

def reconstruct_charge_xy(Gx, Gy, rx, ry, q, area, X_grid, Y_grid):
    """Reconstruct in-plane charge density from truncated Fourier series.

    ρ_approx(ρ) = Σ_G S_G cos(G·ρ)
    where S_G = (1/A) Σ_i q_i cos(G·ρ_i).
    """
    rho = np.zeros_like(X_grid)
    for g in range(len(Gx)):
        S_G = (1.0/area) * np.sum(q * np.cos(Gx[g]*rx + Gy[g]*ry))
        rho += S_G * np.cos(Gx[g]*X_grid + Gy[g]*Y_grid)
    return rho

# ============================================================
#  High-level convenience: setup from ion arrays
# ============================================================

class Ewald2D:
    """Container that precomputes reciprocal lattice and coefficients.

    Parameters
    ----------
    a_vec, b_vec : (2,)   in-plane lattice vectors
    rx, ry, rz   : (N,)   ion coordinates
    q             : (N,)   ion charges
    n_harm        : int    G-vector half-width
    """

    def __init__(self, a_vec, b_vec, rx, ry, rz, q, n_harm=3):
        self.a_vec = np.asarray(a_vec, dtype=float)
        self.b_vec = np.asarray(b_vec, dtype=float)
        self.rx = np.asarray(rx, dtype=float)
        self.ry = np.asarray(ry, dtype=float)
        self.rz = np.asarray(rz, dtype=float)
        self.q  = np.asarray(q,  dtype=float)
        self.n_harm = n_harm
        Q_tot = np.sum(self.q)
        if abs(Q_tot) > 1e-6:
            print(f"WARNING Ewald2D: cell not neutral, Q_tot={Q_tot:.6f}")
        self.area, self.b1, self.b2 = make_reciprocal_2d(self.a_vec, self.b_vec)
        self.hs, self.ks, self.Gx, self.Gy, self.Gn = generate_G_vectors(self.b1, self.b2, n_harm)
        self.C_G = compute_C_G(self.Gx, self.Gy, self.Gn, self.rx, self.ry, self.rz, self.q, self.area)
        self.w   = compute_w_per_ion(self.Gx, self.Gy, self.Gn, self.rx, self.ry, self.q, self.area)

    @classmethod
    def from_AtomicSystem(cls, sys, n_harm=3):
        """Construct from a pyBall.AtomicSystem (or compatible) object.

        Expects sys.apos (N,3), sys.qs (N,), sys.lvec (3,3).
        The first two rows of lvec are the 2D lattice vectors (xy components).
        """
        assert sys.lvec is not None, "AtomicSystem must have lattice vectors (lvec)"
        assert sys.qs   is not None, "AtomicSystem must have charges (qs)"
        a_vec = sys.lvec[0, :2]
        b_vec = sys.lvec[1, :2]
        rx = sys.apos[:, 0]
        ry = sys.apos[:, 1]
        rz = sys.apos[:, 2]
        q  = sys.qs
        return cls(a_vec, b_vec, rx, ry, rz, q, n_harm=n_harm)

    @property
    def N_G(self):
        return len(self.Gx)

    @property
    def z_max(self):
        return np.max(self.rz)

    @property
    def z_min(self):
        return np.min(self.rz)

    @property
    def Lx(self):
        return np.linalg.norm(self.a_vec)

    @property
    def Ly(self):
        return np.linalg.norm(self.b_vec)

    def phi_vacuum_xy(self, X, Y, z):
        """Potential on XY grid at height z (vacuum side)."""
        return eval_potential_vacuum(self.C_G, self.Gx, self.Gy, self.Gn, X, Y, z)

    def phi_full_2d(self, X, Y, Z):
        """Full potential on 2D grid (e.g. XZ slice)."""
        return eval_potential_full_2d(self.w, self.Gx, self.Gy, self.Gn, self.rz, self.q, self.area, X, Y, Z)

    def phi_full_1d(self, x0, y0, z_arr):
        """Full potential along vertical line at (x0, y0)."""
        return eval_potential_full_1d(self.w, self.Gx, self.Gy, self.Gn, self.rz, self.q, self.area, x0, y0, z_arr)

    def phi_brute_1d(self, x0, y0, z_arr, N_rep=20):
        """Brute-force Coulomb sum along vertical line."""
        return eval_potential_brute(self.rx, self.ry, self.rz, self.q, self.a_vec, self.b_vec, x0, y0, z_arr, N_rep=N_rep)

    def charge_density_xy(self, X, Y):
        """Reconstructed in-plane charge density."""
        return reconstruct_charge_xy(self.Gx, self.Gy, self.rx, self.ry, self.q, self.area, X, Y)

    def make_xy_grid(self, Ng=200):
        """Create XY meshgrid covering one unit cell."""
        xv = np.linspace(0, self.Lx, Ng)
        yv = np.linspace(0, self.Ly, Ng)
        return np.meshgrid(xv, yv)

    def make_xz_grid(self, Ng=200, z_lo=None, z_hi=None, y_fixed=0.0):
        """Create XZ meshgrid for a vertical slice at fixed y."""
        if z_lo is None: z_lo = self.z_min - 1.0
        if z_hi is None: z_hi = self.z_max + 4.0
        xv = np.linspace(0, self.Lx, Ng)
        zv = np.linspace(z_lo, z_hi, Ng)
        X, Z = np.meshgrid(xv, zv)
        Y = np.full_like(X, y_fixed)
        return X, Y, Z

    def print_info(self):
        print(f"Ewald2D: N_ions={len(self.q)}  Q_tot={np.sum(self.q):.6f}")
        print(f"  a = ({self.a_vec[0]:.4f}, {self.a_vec[1]:.4f})  |a| = {self.Lx:.4f}")
        print(f"  b = ({self.b_vec[0]:.4f}, {self.b_vec[1]:.4f})  |b| = {self.Ly:.4f}")
        print(f"  Area = {self.area:.4f} Å²")
        print(f"  b1 = ({self.b1[0]:.4f}, {self.b1[1]:.4f})  |b1| = {np.linalg.norm(self.b1):.4f}")
        print(f"  b2 = ({self.b2[0]:.4f}, {self.b2[1]:.4f})  |b2| = {np.linalg.norm(self.b2):.4f}")
        print(f"  n_harm = {self.n_harm}  N_G = {self.N_G}")
        print(f"  z_range = [{self.z_min:.4f}, {self.z_max:.4f}]")
