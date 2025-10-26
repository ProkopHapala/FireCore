import numpy as np

from .MMFF import MMFF, initAtomProperties


def _law_of_cosines(r_ab, r_bc, cos_theta):
    term = r_ab * r_ab + r_bc * r_bc - 2.0 * r_ab * r_bc * cos_theta
    if term < 0.0:
        term = 0.0
    return np.sqrt(term, dtype=np.float64)


class MMFFL(MMFF):
    """Linearized MM force-field builder.

    Converts angular/planar restraints into equivalent harmonic bonds to enable
    linear solvers. Reuses MMFF topology detection, then augments the system
    with extra bonds and optional pi dummy atoms.
    """

    def __init__(self, *, L_pi=1.0, two_pi_dummies=False, Kang=0.0, Kpi_host=0.0, Kpi_orth=0.0, verbosity=1, lone_pairs_pi=False, align_pi_vectors=False):
        super().__init__(bTorsion=False, verbosity=verbosity)
        self.L_pi = float(L_pi)
        self.two_pi_dummies = bool(two_pi_dummies)
        self.K_ang = float(Kang)
        self.K_pi_host = float(Kpi_host)
        self.K_pi_orth = float(Kpi_orth)
        self.lone_pairs_pi = bool(lone_pairs_pi)
        self.align_pi_vectors = bool(align_pi_vectors)
        self.linear_bonds = []      # [(i, j, l0, k, tag)]
        self.pi_dummies = []        # [{'index': idx, 'host': ia, 'sign': +/-1.0}]
        self._dummy_start = 0
        self._next_dummy = 0

    def reset_linearization_state(self):
        self.linear_bonds = []
        self.pi_dummies = []
        self._dummy_start = self.natoms
        self._next_dummy = self.natoms

    def build_linearized(self, mol, atom_types=None, *, bUFF=False):
        if atom_types is None:
            atom_types = self.atom_types
        if mol.ngs is None:
            mol.neighs()
        npi_list, nep_list, is_node = initAtomProperties(mol, atom_types, self.capping_atoms)
        if isinstance(npi_list, np.ndarray): npi_list = npi_list.tolist()
        if isinstance(nep_list, np.ndarray): nep_list = nep_list.tolist()
        if isinstance(is_node, np.ndarray): is_node = is_node.tolist()
        mol.npi_list = npi_list
        mol.nep_list = nep_list
        mol.isNode = is_node

        super().toMMFFsp3_loc(
            mol=mol,
            atom_types=atom_types,
            bRealloc=True,
            bEPairs=False,
            bUFF=bUFF,
            lone_pairs_pi=self.lone_pairs_pi,
            align_pi_vectors=self.align_pi_vectors,
        )
        self.reset_linearization_state()
        self._build_angle_bonds(mol, atom_types)
        host_npi = getattr(mol, "npi_list", npi_list)
        if isinstance(host_npi, np.ndarray):
            host_npi = host_npi.tolist()
        host_npi = [int(x) for x in host_npi]
        self._build_pi_dummies(atom_types, host_npi)
        self._sync_pi_block()
        return self

    def _sync_pi_block(self):
        if self.pipos is None:
            return
        block = self.pipos[:self.nnode]
        self.apos[self.natoms:self.natoms + self.nnode, :3] = block
        self.apos[self.natoms:self.natoms + self.nnode, 3] = 0.0

    def _build_angle_bonds(self, mol, atom_types):
        angs, triplets = mol.findAngles(ngs=mol.ngs)
        for (ja, ib, jc) in triplets:
            if ja == jc:
                continue
            r_ab = self._equilibrium_bond_length(ib, ja)
            r_bc = self._equilibrium_bond_length(ib, jc)
            if r_ab <= 0.0 or r_bc <= 0.0:
                continue
            atype = atom_types.get(mol.enames[ib])
            if atype is None:
                continue
            theta = np.deg2rad(atype.Ass)
            l0 = _law_of_cosines(r_ab, r_bc, np.cos(theta))
            k = self.K_ang if self.K_ang != 0.0 else float(atype.Kss)
            self._register_linear_bond(ja, jc, l0, k, ("angle", (ja, ib, jc)))

    def _build_pi_dummies(self, atom_types, npi_list):
        n_per_pi = 2 if self.two_pi_dummies else 1
        valid_hosts = []
        for ia in range(self.nnode):
            npi = int(npi_list[ia])
            if npi <= 0:
                continue
            norm = np.linalg.norm(self.pipos[ia])
            if not np.isfinite(norm) or norm < 1e-6:
                continue
            valid_hosts.append(ia)
        total = len(valid_hosts) * n_per_pi
        self._extend_storage(total)
        for ia in valid_hosts:
            pi_dir = self.pipos[ia]
            pi_dir = pi_dir / np.linalg.norm(pi_dir)
            host_pos = self.apos[ia, :3]
            for idx_side, sign in enumerate(([+1.0, -1.0] if self.two_pi_dummies else [+1.0])):
                idx = self._next_dummy
                self._next_dummy += 1
                pos = host_pos + pi_dir * (self.L_pi * sign)
                self.apos[idx, :3] = pos.astype(np.float32)
                self.apos[idx, 3] = 0.0
                self.fapos[idx, :] = 0.0
                self.atypes[idx] = -1
                self.REQs[idx, :] = 0.0
                self.neighs[idx, :] = -1
                self.neighs[idx, 0] = ia
                self.neighCell[idx, :] = 0
                self.pi_dummies.append({'index': idx, 'host': ia, 'sign': sign})
                self._register_linear_bond(ia, idx, self.L_pi, self.K_pi_host, ("pi-host", (ia, idx)))
                for nb in self.neighs[ia]:
                    if nb < 0 or nb >= self.natoms:
                        continue
                    r_ab = self._equilibrium_bond_length(ia, nb)
                    if r_ab <= 0.0:
                        continue
                    l0 = np.sqrt(self.L_pi * self.L_pi + r_ab * r_ab)
                    self._register_linear_bond(idx, nb, l0, self.K_pi_orth, ("pi-orth", (idx, nb)))
        if self._next_dummy != self._dummy_start + total:
            raise ValueError("Pi dummy bookkeeping mismatch")

    def _extend_storage(self, n_extra):
        self._dummy_start = self.natoms
        self._next_dummy = self.natoms
        if n_extra <= 0:
            return
        old_atoms = self.natoms
        old_vecs = self.nvecs
        new_atoms = old_atoms + n_extra
        new_vecs = new_atoms + self.nnode
        old_pi_block = self.apos[old_atoms:old_atoms + self.nnode, :3].copy()
        old_pi_force = self.fapos[old_atoms:old_atoms + self.nnode].copy()
        self.apos = self._grow_matrix(self.apos, (new_vecs, 4), fill=0.0)
        self.fapos = self._grow_matrix(self.fapos, (new_vecs, 4), fill=0.0)
        self.atypes = self._grow_vector(self.atypes, new_atoms, fill=-1)
        self.neighs = self._grow_matrix(self.neighs, (new_atoms, 4), fill=-1, dtype=np.int32)
        self.neighCell = self._grow_matrix(self.neighCell, (new_atoms, 4), fill=0, dtype=np.int32)
        self.REQs = self._grow_matrix(self.REQs, (new_atoms, 4), fill=0.0)
        self.pipos = self._grow_matrix(self.pipos, (new_atoms, 3), fill=0.0)
        if hasattr(self, 'back_neighs') and self.back_neighs is not None:
            self.back_neighs = self._grow_matrix(self.back_neighs, (new_vecs, 4), fill=-1, dtype=np.int32)
        self.apos[new_atoms:new_atoms + self.nnode, :4] = 0.0
        self.apos[new_atoms:new_atoms + self.nnode, :3] = old_pi_block
        self.fapos[new_atoms:new_atoms + self.nnode, :] = old_pi_force
        self.ncap += n_extra
        self.natoms = new_atoms
        self.nvecs = new_vecs
        self.nDOFs = self.nvecs * 3

    @staticmethod
    def _grow_matrix(src, shape, *, fill=0.0, dtype=None):
        out = np.full(shape, fill, dtype=dtype or src.dtype)
        sh0 = min(src.shape[0], shape[0])
        sh1 = min(src.shape[1], shape[1])
        out[:sh0, :sh1] = src[:sh0, :sh1]
        return out

    @staticmethod
    def _grow_vector(src, size, *, fill=0, dtype=None):
        out = np.full(size, fill, dtype=dtype or src.dtype)
        out[:min(src.shape[0], size)] = src[:min(src.shape[0], size)]
        return out

    def _equilibrium_bond_length(self, i, j):
        if i < self.nnode:
            neighs = self.neighs[i]
            for idx, nb in enumerate(neighs):
                if nb == j and self.bLs is not None:
                    return float(self.bLs[i, idx])
        if j < self.nnode:
            neighs = self.neighs[j]
            for idx, nb in enumerate(neighs):
                if nb == i and self.bLs is not None:
                    return float(self.bLs[j, idx])
        pa = self.apos[i, :3]
        pb = self.apos[j, :3]
        return float(np.linalg.norm(pa - pb))

    def _register_linear_bond(self, i, j, l0, k, tag):
        self.linear_bonds.append((int(i), int(j), float(l0), float(k), tag))
