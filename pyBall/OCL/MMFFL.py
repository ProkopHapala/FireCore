import numpy as np
import os

from .MMFF import MMFF, initAtomProperties
from . import MMparams


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

    def __init__(self, *, L_pi=1.0, two_pi_dummies=False, Kang=0.0, Kpi_host=0.0, Kpi_orth=0.0, Kpi_align=0.0, L_epair=0.5, Kep_host=0.0, Kep_orth=0.0, Kep_pair=0.0, verbosity=1, lone_pairs_pi=False, align_pi_vectors=False, reorder_nodes_first=True):
        super().__init__(bTorsion=False, verbosity=verbosity, reorder_nodes_first=reorder_nodes_first)
        self.L_pi = float(L_pi)
        self.two_pi_dummies = bool(two_pi_dummies)
        self.K_ang = float(Kang)
        self.K_pi_host = float(Kpi_host)
        self.K_pi_orth = float(Kpi_orth)
        self.K_pi_align = float(Kpi_align)
        self.L_epair = float(L_epair)
        self.K_ep_host = float(Kep_host)
        self.K_ep_orth = float(Kep_orth)
        self.K_ep_pair = float(Kep_pair)
        self.lone_pairs_pi = bool(lone_pairs_pi)
        self.align_pi_vectors = bool(align_pi_vectors)
        self.linear_bonds = []      # [(i, j, l0, k, tag)]
        self.pi_dummies = []        # [{'index': idx, 'host': ia, 'sign': +/-1.0}]
        self.ep_dummies = []        # [{'index': idx, 'host': ia, 'slot': k}]
        self._dummy_start = 0
        self._next_dummy = 0

    def reset_linearization_state(self):
        self.linear_bonds = []
        self.pi_dummies = []
        self.ep_dummies = []
        self._dummy_start = self.natoms
        self._next_dummy = self.natoms

    @staticmethod
    def _default_params_path():
        base_path = os.path.dirname(os.path.abspath(__file__))
        return os.path.join(base_path, "../../cpp/common_resources/")

    @staticmethod
    def _load_params(path=None):
        if path is None:
            path = MMFFL._default_params_path()
        return MMparams.MMFFparams(path)

    @staticmethod
    def _resolve_type_name_table(ename, params: MMparams.MMFFparams):
        if ename in params.atom_types_map:
            return ename
        trial = ename + "_"
        if trial in params.atom_types_map:
            return trial
        # fallback: first type whose element matches
        for name, at in params.atom_types_map.items():
            if at.element_name == ename:
                return name
        raise ValueError(f"Cannot resolve atom type for ename='{ename}' using AtomTypes.dat")

    @staticmethod
    def assign_type_names(mol, *, type_source="table", params_path=None, bSimple=True, b141=True, bConj=False, bCumulene=False):
        """Return list[str] of type names (same length as mol atoms).

        type_source:
            - 'table' : use AtomTypes.dat mapping directly (ename->type key)
            - 'uff'   : run UFF topology inference to assign types
        """
        if mol.enames is None:
            raise ValueError("mol.enames is None; cannot assign types")
        params = MMFFL._load_params(params_path)
        if type_source == "table":
            return [MMFFL._resolve_type_name_table(str(e), params) for e in mol.enames]
        if type_source != "uff":
            raise ValueError(f"Unknown type_source='{type_source}'")
        from .UFFbuilder import UFF_Builder
        uff = UFF_Builder(mol, bSimple=bSimple, b141=b141, bConj=bConj, bCumulene=bCumulene)
        uff.assign_uff_types()
        out = []
        for ia in range(len(mol.enames)):
            it = int(uff.builder_atoms[ia].type)
            out.append(uff.params.atypes[it].name)
        return out

    @staticmethod
    def _count_nbond_for_atom(ng):
        if ng is None:
            return 0
        if isinstance(ng, dict):
            return len(ng)
        return sum(1 for x in ng if int(x) >= 0)

    @staticmethod
    def report_types_and_valence(mol, type_names, *, params_path=None, do_assert=True, atoms_assert=("C", "N", "O")):
        params = MMFFL._load_params(params_path)
        if mol.ngs is None:
            mol.neighs()
        for ia, tname in enumerate(type_names):
            at = params.atom_types_map.get(tname)
            if at is None:
                raise ValueError(f"Type '{tname}' missing in AtomTypes.dat")
            nbond = MMFFL._count_nbond_for_atom(mol.ngs[ia])
            nepair = int(getattr(at, "nepair", 0))
            npi = int(getattr(at, "valence", 0)) - nbond
            if npi < 0:
                npi = 0
            print(f"TYPE ia={ia:4d} ename={mol.enames[ia]} type={tname:8s} nbond={nbond} npi={npi} nepair={nepair}")
            elem = str(getattr(at, "element_name", ""))
            if do_assert and (elem in atoms_assert):
                if (nbond + npi + nepair) != 4:
                    raise ValueError(f"Octet-rule check failed ia={ia} elem={elem} nbond+npi+nepair={nbond+npi+nepair} != 4")

    def build_topology(self, mol, *, type_source="table", params_path=None, add_angle=True, add_pi=False, two_pi_dummies=False, add_pi_align=False, add_epair=False, add_epair_pairs=False, bUFF=False, report=True):
        """Build linearized topology for XPDB.

        Returns dict with:
            - apos        : (n_all,3)
            - type_names  : list[str] (n_all)
            - is_dummy    : (n_all,) bool
            - bonds_primary : list[(i,j)] among real atoms
            - bonds_angle   : list[(i,j)]
            - bonds_pi      : list[(i,j)]
            - bonds_epair   : list[(i,j)]
            - bonds_linear  : list[(i,j,l0,k,tag)] (from self.linear_bonds)
        """
        if mol.ngs is None:
            mol.neighs()
        type_names_real = self.assign_type_names(mol, type_source=type_source, params_path=params_path)
        if report:
            self.report_types_and_valence(mol, type_names_real, params_path=params_path)

        params = self._load_params(params_path)
        atom_types = params.atom_types_map

        enames_save = list(mol.enames)
        mol.enames = list(type_names_real)
        try:
            self.two_pi_dummies = bool(two_pi_dummies)
            self.build_linearized(mol, atom_types=atom_types, bUFF=bUFF, include_angle=bool(add_angle), include_pi=bool(add_pi))
            if add_pi and add_pi_align and mol.bonds is not None and len(mol.bonds) > 0:
                self._build_pi_align_bonds(mol.bonds)
            if add_epair:
                host_npi = getattr(mol, "npi_list", [0] * len(type_names_real))
                host_nep = getattr(mol, "nep_list", [0] * len(type_names_real))
                host_npi = [int(x) for x in (host_npi.tolist() if isinstance(host_npi, np.ndarray) else host_npi)]
                host_nep = [int(x) for x in (host_nep.tolist() if isinstance(host_nep, np.ndarray) else host_nep)]
                self._build_epair_dummies(atom_types, host_npi, host_nep)
                if add_epair_pairs:
                    self._build_epair_pair_bonds()
        finally:
            mol.enames = enames_save

        n_real = len(type_names_real)
        n_all = int(self.natoms)
        is_dummy = np.zeros(n_all, dtype=bool)
        if n_all > n_real:
            is_dummy[n_real:n_all] = True

        type_names_all = list(type_names_real)
        if n_all > n_real:
            type_names_all += ["Pi"] * (n_all - n_real)
            # epair dummies override names
            for d in self.ep_dummies:
                idx = int(d["index"])
                if idx < len(type_names_all):
                    type_names_all[idx] = "E"

        bonds_primary = []
        if mol.bonds is not None:
            for a, b in mol.bonds:
                ia = int(a); ib = int(b)
                bonds_primary.append((min(ia, ib), max(ia, ib)))

        bonds_angle = []
        bonds_pi = []
        bonds_pi_align = []
        bonds_ep = []
        bonds_ep_pair = []
        for (i, j, l0, k, tag) in self.linear_bonds:
            kind = tag[0] if isinstance(tag, tuple) else str(tag)
            ij = (min(int(i), int(j)), max(int(i), int(j)))
            if kind == "angle":
                bonds_angle.append(ij)
            elif kind.startswith("pi-align"):
                bonds_pi_align.append(ij)
            elif kind.startswith("pi"):
                bonds_pi.append(ij)
            elif kind.startswith("ep-pair"):
                bonds_ep_pair.append(ij)
            elif kind.startswith("ep"):
                bonds_ep.append(ij)

        return {
            "apos": self.apos[:n_all, :3].copy(),
            "type_names": type_names_all,
            "is_dummy": is_dummy,
            "bonds_primary": sorted(set(bonds_primary)),
            "bonds_angle": sorted(set(bonds_angle)),
            "bonds_pi": sorted(set(bonds_pi)),
            "bonds_pi_align": sorted(set(bonds_pi_align)),
            "bonds_epair": sorted(set(bonds_ep)),
            "bonds_epair_pair": sorted(set(bonds_ep_pair)),
            "bonds_linear": list(self.linear_bonds),
            "n_real": n_real,
            "n_all": n_all,
        }

    def build_linearized(self, mol, atom_types=None, *, bUFF=False, include_angle=True, include_pi=True):
        if atom_types is None:
            atom_types = self.atom_types
        if not hasattr(mol, "natoms"):
            mol.natoms = int(len(mol.apos))
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
        if include_angle:
            self._build_angle_bonds(mol, atom_types)
        if include_pi:
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

    def _build_epair_dummies(self, atom_types, npi_list, nep_list):
        if self.nnode <= 0:
            return
        # count total electron-pair dummies
        hosts = []
        for ia in range(self.nnode):
            nep = int(nep_list[ia]) if ia < len(nep_list) else 0
            if nep <= 0:
                continue
            hosts.append((ia, nep))
        total = sum(nep for _, nep in hosts)
        self._extend_storage(total)
        for ia, nep in hosts:
            nb_ids = [int(nb) for nb in self.neighs[ia] if int(nb) >= 0 and int(nb) < self.natoms]
            nb = len(nb_ids)
            npi = int(npi_list[ia]) if ia < len(npi_list) else 0
            dirs = self._epair_dirs(ia, nb_ids, npi, nep)
            host_pos = self.apos[ia, :3]
            for k, d in enumerate(dirs):
                idx = self._next_dummy
                self._next_dummy += 1
                pos = host_pos + d * self.L_epair
                self.apos[idx, :3] = pos.astype(np.float32)
                self.apos[idx, 3] = 0.0
                self.fapos[idx, :] = 0.0
                self.atypes[idx] = -1
                self.REQs[idx, :] = 0.0
                self.neighs[idx, :] = -1
                self.neighs[idx, 0] = ia
                self.neighCell[idx, :] = 0
                self.ep_dummies.append({'index': idx, 'host': ia, 'slot': k})
                self._register_linear_bond(ia, idx, self.L_epair, self.K_ep_host, ("ep-host", (ia, idx)))
                for nbj in nb_ids:
                    r_ab = self._equilibrium_bond_length(ia, nbj)
                    if r_ab <= 0.0:
                        continue
                    l0 = np.sqrt(self.L_epair * self.L_epair + r_ab * r_ab)
                    self._register_linear_bond(idx, nbj, l0, self.K_ep_orth, ("ep-orth", (idx, nbj)))
        if self._next_dummy != self._dummy_start + total:
            raise ValueError("Epair dummy bookkeeping mismatch")

    def _build_epair_pair_bonds(self):
        if not self.ep_dummies or self.K_ep_pair == 0.0:
            return
        by_host = {}
        for d in self.ep_dummies:
            by_host.setdefault(int(d["host"]), []).append(int(d["index"]))
        for host, idxs in by_host.items():
            if len(idxs) < 2:
                continue
            for a in range(len(idxs)):
                for b in range(a + 1, len(idxs)):
                    ia = idxs[a]; ib = idxs[b]
                    pos_a = self.apos[ia, :3]
                    pos_b = self.apos[ib, :3]
                    l0 = float(np.linalg.norm(pos_a - pos_b))
                    self._register_linear_bond(ia, ib, l0, self.K_ep_pair, ("ep-pair", (host, ia, ib)))

    def _build_pi_align_bonds(self, primary_bonds):
        if not self.pi_dummies or self.K_pi_align == 0.0:
            return
        host_map = {}
        for d in self.pi_dummies:
            host_map.setdefault(int(d["host"]), []).append(int(d["index"]))
        if primary_bonds is None:
            return
        if hasattr(primary_bonds, "__len__") and len(primary_bonds) == 0:
            return
        for a, b in primary_bonds:
            ia = int(a); ib = int(b)
            list_a = host_map.get(ia)
            list_b = host_map.get(ib)
            if not list_a or not list_b:
                continue
            for idx_a in list_a:
                for idx_b in list_b:
                    pos_a = self.apos[idx_a, :3]
                    pos_b = self.apos[idx_b, :3]
                    l0 = float(np.linalg.norm(pos_a - pos_b))
                    self._register_linear_bond(idx_a, idx_b, l0, self.K_pi_align, ("pi-align", (ia, ib, idx_a, idx_b)))

    def _epair_dirs(self, ia, nb_ids, npi, nep):
        # Mirror AtomicSystem.make_epair_geom shapes (but return directions, do not mutate mol)
        pos = self.apos[ia, :3]
        vs = []
        for j in nb_ids[:3]:
            v = self.apos[j, :3] - pos
            n = np.linalg.norm(v)
            if n > 1e-8:
                vs.append((v / n).astype(np.float32))
        nb = len(vs)
        out = []
        if nep <= 0:
            return out
        if npi == 0:
            if nb == 3:
                base = np.cross(vs[1] - vs[0], vs[2] - vs[0])
                n = np.linalg.norm(base)
                if n > 1e-8:
                    base = base / n
                else:
                    base = np.array([0.0, 0.0, 1.0], dtype=np.float32)
                if np.dot(base, (vs[0] + vs[1] + vs[2])) > 0:
                    base = -base
                out.append(base.astype(np.float32))
            elif nb == 2:
                m_c = vs[0] + vs[1]
                nmc = np.linalg.norm(m_c)
                if nmc > 1e-8:
                    m_c = m_c / nmc
                else:
                    m_c = np.array([0.0, 0.0, 1.0], dtype=np.float32)
                m_b = np.cross(vs[0], vs[1])
                nmb = np.linalg.norm(m_b)
                if nmb > 1e-8:
                    m_b = m_b / nmb
                else:
                    m_b = np.array([0.0, 1.0, 0.0], dtype=np.float32)
                cc = 0.57735026919
                cb = 0.81649658092
                ep1 = (-cc * m_c + cb * m_b).astype(np.float32)
                ep2 = (-cc * m_c - cb * m_b).astype(np.float32)
                ep1 /= max(np.linalg.norm(ep1), 1e-8)
                ep2 /= max(np.linalg.norm(ep2), 1e-8)
                out += [ep1, ep2]
        elif npi == 1:
            if nb == 2:
                m_c = vs[0] + vs[1]
                nmc = np.linalg.norm(m_c)
                if nmc > 1e-8:
                    m_c = m_c / nmc
                else:
                    m_c = np.array([0.0, 0.0, 1.0], dtype=np.float32)
                out.append((-m_c).astype(np.float32))
            elif nb == 1:
                v1 = vs[0]
                m_b = self.pipos[ia] if (self.pipos is not None and ia < self.pipos.shape[0]) else np.array([0.0, 0.0, 1.0], dtype=np.float32)
                nmb = np.linalg.norm(m_b)
                if nmb > 1e-8:
                    m_b = m_b / nmb
                else:
                    m_b = np.array([0.0, 0.0, 1.0], dtype=np.float32)
                m_c = np.cross(v1, m_b)
                nmc = np.linalg.norm(m_c)
                if nmc > 1e-8:
                    m_c = m_c / nmc
                else:
                    # pick arbitrary perpendicular
                    m_c = np.cross(v1, np.array([0.0, 1.0, 0.0], dtype=np.float32))
                    m_c /= max(np.linalg.norm(m_c), 1e-8)
                out.append((-0.5 * v1 + 0.86602540378 * m_c).astype(np.float32))
                out.append((-0.5 * v1 - 0.86602540378 * m_c).astype(np.float32))
        # clamp to requested count
        if len(out) > nep:
            out = out[:nep]
        # normalize
        normed = []
        for d in out:
            n = np.linalg.norm(d)
            if n > 1e-8:
                normed.append((d / n).astype(np.float32))
        return normed

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
        self.apars = self._grow_matrix(self.apars, (new_atoms, 4), fill=0.0)
        self.bLs = self._grow_matrix(self.bLs, (new_atoms, 4), fill=0.0)
        self.bKs = self._grow_matrix(self.bKs, (new_atoms, 4), fill=0.0)
        self.Ksp = self._grow_matrix(self.Ksp, (new_atoms, 4), fill=0.0)
        self.Kpp = self._grow_matrix(self.Kpp, (new_atoms, 4), fill=0.0)
        self.angles = self._grow_tensor(self.angles, (new_atoms, 6, 3), fill=0.0)
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
    def _grow_tensor(src, shape, *, fill=0.0, dtype=None):
        out = np.full(shape, fill, dtype=dtype or src.dtype)
        slices = tuple(slice(0, min(src.shape[i], shape[i])) for i in range(len(shape)))
        out[slices] = src[slices]
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
