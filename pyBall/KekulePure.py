#!/usr/bin/env python3
"""
Pure-Python NumPy Kekule bond-order optimizer.

Optimizes pi-bond orders for sp2 atoms using a flat-array, vectorized
relaxation.  The sigma skeleton is taken from AtomicSystem.bonds; only
the extra pi contribution is optimized.

Algorithm:
    1. Atom valence: parabolic penalty around atom-specific n_pi
       (sum of pi orders at each atom).  sp2 atoms have n_pi=1, sp3 have 0.
    2. Aromatic bond energy: parabola centered at pi=0.5.
    3. Localization snap: three piecewise parabolas around 0, 0.5, and 1.
    4. Gradient-descent relaxation.
"""

import numpy as np


def _to_int_array(a):
    a = np.asarray(a)
    if a.dtype != np.int32:
        a = a.astype(np.int32)
    return a


class KekulePure:
    """
    Kekule pi-bond order optimizer.

    Parameters
    ----------
    system : AtomicSystem
        Must have `bonds` set.  Each bond is assumed to carry one sigma bond.
    n_pi : array-like of int, optional
        Target number of pi electrons per atom (e.g. 1 for sp2, 0 for sp3).
        If None, it is inferred from `pi_atoms`.
    pi_atoms : array-like of bool, optional
        Mask of atoms that participate in the pi system.  Only used when
        `n_pi` is not given.  If both are None, all atoms are treated as sp2
        (n_pi = 1).
    bonds : array-like of (int, int), optional
        Subset of bonds to optimize.  If None, all `system.bonds` are used.
    Kval : float
        Stiffness of the atom valence (pi-electron count) penalty.
    Kloc : float
        Stiffness of the three-piece localization snap potential (minima at
        pi=0, 0.5, 1).  Usually turned on after an aromatic pre-relaxation.
    Karo : float
        Aromatic stabilization energy.  A positive value makes pi=0.5 (total
        bond order 1.5) the minimum of the bond energy.
    Kbound : float
        Stiffness of the hard [0,1] bounds.
    """

    def __init__(self, system, n_pi=None, pi_atoms=None, bonds=None,
                 Kval=1.0, Kloc=0.0, Karo=0.3, Kbound=1.0, bRandomStart=False,
                 allow_aromatic=False):
        self.system = system
        self.natom = system.natoms
        self.allow_aromatic = allow_aromatic

        # bonds to optimize: (nbond, 2) int array
        if bonds is None:
            bonds = system.bonds
        self.bonds = _to_int_array(bonds)
        self.nbond = len(self.bonds)

        # target number of pi electrons per atom
        if n_pi is None:
            if pi_atoms is None:
                pi_atoms = np.ones(self.natom, dtype=bool)
            pi_atoms = np.asarray(pi_atoms, dtype=bool)
            n_pi = np.where(pi_atoms, 1.0, 0.0)
        self.n_pi = np.asarray(n_pi, dtype=float)

        self.Kval = Kval
        self.Kloc = Kloc
        self.Karo = Karo
        self.Kbound = Kbound

        # current pi bond orders
        if bRandomStart:
            self.bo = np.random.rand(self.nbond)
        else:
            self.bo = np.full(self.nbond, 0.5)
            # Random start is strongly recommended for discrete Kekule patterns
            # because the uniform 0.5 start rounds to 1 and all bonds get
            # pushed toward double bonds.
            self.bo = np.random.rand(self.nbond)

        # split atom indices for projection
        self.i0 = self.bonds[:, 0]
        self.i1 = self.bonds[:, 1]

        # incidence (atom x bond) for quadratic solves
        self._A = np.zeros((self.natom, self.nbond), dtype=float)
        ib = np.arange(self.nbond, dtype=np.int32)
        np.add.at(self._A, (self.i0, ib), 1.0)
        np.add.at(self._A, (self.i1, ib), 1.0)
        self._AT = self._A.T
        self._ATA = self._AT @ self._A

    def project_valence(self):
        """Return atom valences = sum of connected pi bond orders."""
        val = np.zeros(self.natom, dtype=float)
        np.add.at(val, self.i0, self.bo)
        np.add.at(val, self.i1, self.bo)
        return val

    def _localization_target(self, bo):
        """Return target bond order for the three-piece snap parabolas.

        Switching thresholds are 0.25 and 0.75.  If allow_aromatic is False,
        the middle basin is suppressed and the target is the nearest integer.
        """
        target = np.empty_like(bo)
        mask_0 = bo < 0.25
        mask_1 = bo > 0.75
        mask_a = ~(mask_0 | mask_1)
        target[mask_0] = 0.0
        target[mask_1] = 1.0
        if self.allow_aromatic:
            target[mask_a] = 0.5
        else:
            target[mask_a] = np.where(bo[mask_a] > 0.5, 1.0, 0.0)
        return target

    def set_random_bonds(self):
        self.bo = np.random.rand(self.nbond)

    def eval(self):
        """Evaluate forces and energy.  Returns total energy."""
        bo = self.bo

        # atom valence: E = Kval * (n_pi - sum_j pi_ij)^2
        val = self.project_valence()
        d_val = val - self.n_pi
        f_atom = -2.0 * self.Kval * d_val
        E_val = self.Kval * np.sum(d_val * d_val)

        # aromatic parabola centered at 0.5: E = Karo * (pi - 0.5)^2
        d_aro = bo - 0.5
        f_bond = -2.0 * self.Karo * d_aro
        E_aro = self.Karo * np.sum(d_aro * d_aro)

        # localization snap parabolas around 0 / 0.5 / 1
        target = self._localization_target(bo)
        d_loc = bo - target
        f_bond += -2.0 * self.Kloc * d_loc
        E_loc = self.Kloc * np.sum(d_loc * d_loc)

        # hard bounds [0,1]
        d_min = np.minimum(bo - 0.0, 0.0)
        d_max = np.maximum(bo - 1.0, 0.0)
        f_bond -= self.Kbound * (d_min + d_max)
        E_bound = 0.5 * self.Kbound * np.sum(d_min * d_min + d_max * d_max)

        # project atom forces back onto the bonds
        f_bond += f_atom[self.i0] + f_atom[self.i1]

        self._f_bond = f_bond
        self._f_atom = f_atom
        return E_val + E_aro + E_loc + E_bound

    def step(self, dt):
        """Gradient-descent step.  Returns squared force norm."""
        self.eval()
        self.bo += dt * self._f_bond
        self.bo = np.clip(self.bo, 0.0, 1.0)
        return float(np.sum(self._f_bond * self._f_bond))

    def relax(self, dt=0.1, nmax=2000, tol=1e-6, verbose=False):
        """Relax bond orders.  Returns F2 convergence metric."""
        F2 = 1.0
        for it in range(nmax):
            F2 = self.step(dt)
            if verbose and (it % 100 == 0 or it == nmax - 1):
                E = self.eval()
                print(f"iter {it:4d}  E={E:.6f}  F2={F2:.6e}")
            if F2 < tol:
                break
        return F2

    def relax_multistart(self, ntrials=20, dt=0.1, nmax=2000, tol=1e-6,
                         Kloc_final=2.0, aromatic_penalty=0.5):
        """
        Run multiple random starts and return the best discrete Kekule pattern.

        The schedule first enforces atom valence (Kval=10, Kloc=0), then
        localizes bonds to 0/1 (Kloc=Kloc_final, Kval=1).  A small penalty for
        aromatic bonds prevents spurious all-aromatic solutions for non-aromatic
        systems.
        """
        Kval_orig = self.Kval
        Kloc_orig = self.Kloc
        best = None
        best_score = 1e300
        for t in range(ntrials):
            self.set_random_bonds()
            # stage 1: satisfy valence
            self.Kval = 10.0
            self.Kloc = 0.0
            self.relax(dt=dt, nmax=nmax, tol=tol)
            # stage 2: localize
            self.Kval = Kval_orig
            self.Kloc = Kloc_final
            self.relax(dt=dt * 0.2, nmax=nmax * 2, tol=tol)
            E = self.eval()
            n_aromatic = np.sum(self.classify() == 1)
            score = E + aromatic_penalty * n_aromatic
            if score < best_score:
                best_score = score
                best = self.bo.copy()
        self.bo = best
        self.Kval = Kval_orig
        self.Kloc = Kloc_orig
        return best_score

    def pi_bond_orders(self):
        """Optimized pi bond orders, shape (nbond,)."""
        return self.bo.copy()

    def total_bond_orders(self):
        """Total bond orders = sigma (1) + pi."""
        return 1.0 + self.bo

    def snap(self, tol=0.15):
        """Round pi bond orders to the nearest discrete value (0, 0.5, 1).

        If allow_aromatic is False, only 0 or 1 are produced.
        """
        if self.allow_aromatic:
            bo = self.bo.copy()
            bo[bo < 0.25 - tol] = 0.0
            bo[bo > 0.75 + tol] = 1.0
            mask = (bo >= 0.25 - tol) & (bo <= 0.75 + tol)
            bo[mask] = 0.5
            self.bo = bo
        else:
            self.bo = np.round(self.bo)
        return self.bo

    def classify(self, tol=0.05):
        """Classify each pi bond as integer codes: 0 single, 1 aromatic, 2 double."""
        bo = self.bo
        out = np.empty(self.nbond, dtype=np.int8)
        if self.allow_aromatic:
            out[:] = 1
            out[bo < 0.25 - tol] = 0
            out[bo > 0.75 + tol] = 2
        else:
            out[:] = 0
            out[bo > 0.5] = 2
        return out

    def make_bond_style(self, tol=0.05):
        """
        Return line widths and colors for plotting bonds.

        Returns
        -------
        lws : (nbond,) array
        colors : (nbond,) array of color strings
        """
        cls = self.classify(tol=tol)
        lws = np.ones(self.nbond, dtype=float)
        colors = np.empty(self.nbond, dtype=object)
        colors[:] = 'k'
        lws[cls == 1] = 2.0
        colors[cls == 1] = 'green'
        lws[cls == 2] = 3.0
        colors[cls == 2] = 'k'
        return lws, colors

    def _solve_constrained_kkt(self, free, target=None, Karo=None, Kloc=None):
        """Solve the KKT system for free bonds with A x = n_pi enforced exactly.

        Returns updated full bond-order vector and the atom Lagrange multipliers.
        """
        if Karo is None: Karo = self.Karo
        if Kloc is None: Kloc = self.Kloc
        if target is None:
            target = np.full(self.nbond, 0.5)
        target = np.asarray(target, dtype=float)
        # fixed bonds are already stored in self.bo
        fixed = ~free
        # residual right-hand side for the atom constraints after fixed bonds
        rhs_atoms = self.n_pi.copy()
        if np.any(fixed):
            rhs_atoms -= self._A[:, fixed] @ self.bo[fixed]
        nfree = int(np.sum(free))
        if nfree == 0:
            return self.bo.copy(), np.zeros(self.natom)
        A_free = self._A[:, free]
        AT_free = A_free.T
        H = 2.0 * (Karo + Kloc) * np.eye(nfree)
        # [ H  A^T ] [ x ] = [ 2*(Karo*0.5 + Kloc*target) ]
        # [ A  0   ] [ l ]   [ rhs_atoms                ]
        top = np.block([[H, AT_free], [A_free, np.zeros((self.natom, self.natom))]])
        rhs_x = 2.0 * (Karo * 0.5 + Kloc * target[free])
        rhs = np.concatenate([rhs_x, rhs_atoms])
        try:
            sol = np.linalg.solve(top, rhs)
        except np.linalg.LinAlgError:
            sol, *_ = np.linalg.lstsq(top, rhs, rcond=None)
        x_new = self.bo.copy()
        x_new[free] = sol[:nfree]
        lam = sol[nfree:]
        return x_new, lam

    def solve_constrained(self, target=None, Karo=None, Kloc=None, max_iter=20,
                          tol=1e-9, clip=True):
        """Solve constrained QP: min energy s.t. A*bo = n_pi and 0 <= bo <= 1.

        Uses an active-set method: solve the unconstrained KKT, clip violating
        bonds to [0,1], fix them, and re-solve the reduced KKT until the active
        set stops changing.
        """
        free = np.ones(self.nbond, dtype=bool)
        for _ in range(max_iter):
            x_new, _ = self._solve_constrained_kkt(free, target=target, Karo=Karo, Kloc=Kloc)
            if clip:
                viol_lo = x_new < 0.0
                viol_hi = x_new > 1.0
                viol = viol_lo | viol_hi
                if np.any(viol):
                    x_new[viol_lo] = 0.0
                    x_new[viol_hi] = 1.0
                    free[viol] = False
                    self.bo = x_new
                    continue
            self.bo = x_new
            break
        # enforce constraint check loudly
        err = self._A @ self.bo - self.n_pi
        max_err = float(np.max(np.abs(err))) if err.size else 0.0
        if max_err > tol:
            raise RuntimeError(f"KekulePure.solve_constrained(): atom-sum constraint not satisfied, max|A@bo-n_pi|={max_err:.3e}")
        return self.bo

    def solve_snap(self, niter=25, Karo=None, Kloc=None):
        """Alternating constrained solve: update snap targets, then solve QP."""
        if Karo is None: Karo = self.Karo
        if Kloc is None: Kloc = self.Kloc
        for _ in range(niter):
            target = self._localization_target(self.bo)
            self.solve_constrained(target=target, Karo=Karo, Kloc=Kloc)
        return self.bo

    def solve_quadratic(self, Kval=None, Karo=None, Kloc=None, target=None, clip=True):
        """Backward-compatible wrapper: use the constrained KKT solver."""
        return self.solve_constrained(target=target, Karo=Karo, Kloc=Kloc)

    def kkt_matrix(self, target=None, Karo=None, Kloc=None):
        """Build the full KKT matrix for plotting/diagnostics.

        Blocks:
            [ 2(Karo+Kloc)I   A^T ]  <- top-left diagonal (bonds), top-right coupling
            [ A               0  ]  <- bottom-left coupling, bottom-right zeros
        """
        if Karo is None: Karo = self.Karo
        if Kloc is None: Kloc = self.Kloc
        H = 2.0 * (Karo + Kloc) * np.eye(self.nbond)
        return np.block([[H, self._AT], [self._A, np.zeros((self.natom, self.natom))]])

    def plot_kkt_matrix(self, target=None, Karo=None, Kloc=None, ax=None,
                        mode='signed', cmap='seismic', fname=None, show=None):
        """Visualize the KKT matrix with labelled blocks.

        Parameters
        ----------
        mode : 'signed' or 'logabs'
            'signed' shows matrix values with a diverging colormap and a
            symmetric color scale (vmin=-vmax).
            'logabs' shows log10(|M_ij|) with a sequential colormap; the colormap
            is reversed so the largest magnitudes are dark.
        cmap : str or Colormap
            'seismic' for signed mode, 'magma' or 'inferno' for logabs mode.
        """
        import matplotlib.pyplot as plt
        M = self.kkt_matrix(target=target, Karo=Karo, Kloc=Kloc)
        if show is None:
            show = (ax is None)
        if (ax is None):
            fig, ax = plt.subplots(figsize=(7, 7))
        if mode == 'signed':
            vmax = np.max(np.abs(M))
            im = ax.imshow(M, cmap=cmap, vmin=-vmax, vmax=vmax, origin='lower')
            plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label='value')
        elif mode == 'logabs':
            Mlog = np.log10(np.abs(M) + 1e-30)
            im = ax.imshow(Mlog, cmap=cmap, origin='lower')
            im.set_cmap(plt.cm.get_cmap(cmap).reversed())
            plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04,
                         label=r'$\log_{10}(|M_{ij}|)$')
        else:
            raise ValueError("mode must be 'signed' or 'logabs'")
        ax.axvline(self.nbond - 0.5, color='red', lw=1)
        ax.axhline(self.nbond - 0.5, color='red', lw=1)
        ax.text(self.nbond * 0.35, self.nbond * 0.35, r'$2(K_{aro}+K_{loc})I$',
                color='red', ha='center', va='center')
        ax.text(self.nbond + self.natom * 0.5, self.nbond * 0.35, r'$A^T$',
                color='red', ha='center', va='center')
        ax.text(self.nbond * 0.35, self.nbond + self.natom * 0.5, r'$A$',
                color='red', ha='center', va='center')
        ax.text(self.nbond + self.natom * 0.5, self.nbond + self.natom * 0.5,
                '0', color='red', ha='center', va='center')
        ax.set_xlabel('variables (bonds | lambdas)')
        ax.set_ylabel('variables (bonds | lambdas)')
        ax.set_title(f'KKT matrix ({mode})')
        if fname is not None:
            plt.tight_layout()
            plt.savefig(fname)
        if show:
            plt.tight_layout()
            plt.show()
        return ax


def make_pi_mask(system, elements=None):
    """
    Return boolean mask of atoms considered sp2/pi atoms.

    Parameters
    ----------
    elements : set of str, optional
        Element symbols treated as pi atoms.  Default is uppercase C/N/O.
    """
    if elements is None:
        elements = {'C', 'N', 'O'}
    return np.array([e in elements for e in system.enames], dtype=bool)


def make_n_pi(system, sp2=None, sp3=None):
    """
    Return n_pi target per atom from element case.

    Uppercase element symbols (C, N, O, ...) are treated as sp2 with one
    pi electron.  Lowercase symbols (c, n, o, ...) are treated as sp3 with
    zero pi electrons.  If the system has an `_enames_original` attribute
    (e.g. from the ASCII parser) it is used, otherwise the current `enames`
    are used.  If the system has no case information, all atoms are assumed
    to be sp2.
    """
    if sp2 is None:
        sp2 = {'C', 'N', 'O'}
    if sp3 is None:
        sp3 = {'c', 'n', 'o'}
    names = np.asarray(getattr(system, '_enames_original', system.enames), dtype=str)
    namesU = np.char.upper(names)
    namesL = np.char.lower(names)
    is_lower = (names == namesL)
    is_sp2 = (~is_lower) & np.isin(namesU, list(sp2))
    is_sp3 = is_lower & np.isin(namesL, list(sp3))
    n_pi = np.where(is_sp2, 1.0, 0.0)
    n_pi = np.where(is_sp3, 0.0, n_pi)
    # Default rule for ambiguous cases (e.g. elements outside {C,N,O}):
    n_pi = np.where((~is_sp2) & (~is_sp3) & np.isin(namesU, list(sp2)) & (~is_lower), 1.0, n_pi)
    return n_pi


def optimize_pi_bonds(system, n_pi=None, pi_atoms=None, bonds=None, dt=0.1,
                      nmax=2000, tol=1e-6, Kval=1.0, Kloc=0.0, Karo=0.3,
                      Kbound=1.0):
    """
    Convenience function: create and relax a KekulePure instance.

    Returns the optimizer instance.
    """
    k = KekulePure(system, n_pi=n_pi, pi_atoms=pi_atoms, bonds=bonds,
                   Kval=Kval, Kloc=Kloc, Karo=Karo, Kbound=Kbound)
    k.relax(dt=dt, nmax=nmax, tol=tol)
    return k
