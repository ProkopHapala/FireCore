import os
import numpy as np
from  .. import MMFF as mmff
from . import MMparams as mmparams
from .. import atomicUtils as au

class Atom:
    def __init__(self, type_idx, iconf):
        self.type = type_idx
        self.iconf = iconf
        self.REQ = np.zeros(4)

class Conformer:
    def __init__(self, nbond):
        self.nbond = nbond

class Bond:
    def __init__(self, atoms, order=1):
        self.atoms = atoms
        self.order = order
        self.l0 = 0.0
        self.k = 0.0

class Angle:
    def __init__(self):
        self.atoms = (0, 0, 0)
        self.bonds = (0, 0)
        self.k = 0.0
        self.C0 = 0.0
        self.C1 = 0.0
        self.C2 = 0.0
        self.C3 = 0.0

class Dihedral:
    def __init__(self):
        self.atoms = (0, 0, 0, 0)
        self.bonds = (0, 0, 0)
        self.k = 0.0
        self.d = 0.0
        self.n = 0.0

class Inversion:
    def __init__(self):
        self.atoms = (0, 0, 0, 0)
        self.bonds = (0, 0, 0)
        self.k = 0.0
        self.C0 = 0.0
        self.C1 = 0.0
        self.C2 = 0.0

class UFF_Builder:

    def __init__(self, mol, bSimple=False, b141=False, bConj=False, bCumulene=False):
        self.mol = mol
        base_path = os.path.dirname(os.path.abspath(__file__))
        data_path = os.path.join(base_path, "../../cpp/common_resources/")
        self.params = mmparams.MMFFparams(data_path)

        # Monkey-patch AtomicSystem instance to be compatible with this builder
        self.builder_atoms = []
        for i, ename in enumerate(self.mol.enames):
            type_idx = self.params.getAtomType(ename, bErr=False)
            if type_idx == -1:
                # If the direct lookup fails, try to find a default type (e.g., "H_" for "H")
                type_idx = self.params.getAtomType(ename + "_", bErr=False)
            if type_idx == -1:
                 # If that also fails, try the element name itself as a fallback
                for j, atype in enumerate(self.params.atypes):
                    if atype.element_name == ename:
                        type_idx = j
                        break
            if type_idx == -1:
                raise Exception(f"Cannot find a default atom type for element {ename}")
            self.builder_atoms.append(Atom(type_idx, i))
            print(f"Atom {i} ({ename}) initial type: {self.params.atypes[type_idx].name}")

        def make_neighs_bonds():
            if self.mol.bonds is None:
                self.mol.findBonds()
            ngs = au.neigh_bonds(len(self.mol.apos), self.mol.bonds)
            max_neighs = 4
            self.neighs = np.full((len(self.mol.apos), max_neighs), -1, dtype=int)
            self.neighBs = np.full((len(self.mol.apos), max_neighs), -1, dtype=int)
            for i, ng in enumerate(ngs):
                for j, (neighbor, bond_idx) in enumerate(ng.items()):
                    self.neighs[i, j] = neighbor
                    self.neighBs[i,j] = bond_idx
            return self.neighs, self.neighBs
        self.mol.make_neighs_bonds = make_neighs_bonds

        def get_bond_by_atoms(i, j):
            if self.bonds is None:
                return None
            for ib, b in enumerate(self.bonds):
                if (b.atoms[0] == i and b.atoms[1] == j) or (b.atoms[0] == j and b.atoms[1] == i):
                    return ib
            return None
        self.mol.get_bond_by_atoms = get_bond_by_atoms

        if self.mol.bonds is None:
            self.mol.findBonds()
        bonds_ = self.mol.bonds
        self.bonds = []
        for b in bonds_:
            self.bonds.append(Bond(b))
        ngs = au.neigh_bonds(len(self.mol.apos), bonds_)
        self.mol.confs = [Conformer(len(n)) for n in ngs]

        self.bSimple = bSimple
        self.b141 = b141
        self.bConj = bConj
        self.bCumulene = bCumulene

    def build(self):
        self.assign_uff_types()
        self.assign_uff_params()
        return self.get_arrays()

    def assign_uff_types(self):
        # This is a Python implementation of the C++ assignUFFtypes logic
        # It will modify self.mol.atom_types, self.mol.bond_orders etc.

        # init working arrays
        tol = 0.05
        set_atom = np.zeros(len(self.builder_atoms), dtype=bool)
        set_bond = np.zeros(len(self.bonds), dtype=bool)
        BOs = np.full(len(self.bonds), -1.0, dtype=np.float64)
        BOs_int = np.full(len(self.bonds), -1, dtype=np.int32)

        # make neighbor list
        neighs, _ = self.mol.make_neighs_bonds()

        # assign bond orders and atom types for trivial cases
        self.assign_uff_types_trivial(neighs, BOs, BOs_int, set_atom, set_bond)
        print("set_atom after trivial:", np.where(~set_atom)[0])

        # assign bond orders and atom types for nitro groups
        self.assign_uff_types_nitro(neighs, BOs, BOs_int, set_atom, set_bond)
        print("set_atom after nitro:", np.where(~set_atom)[0])

        # find a (hopefully) valid limit resonance structure
        self.assign_uff_types_treewalk(neighs, BOs_int)
        print("set_atom after treewalk:", np.where(~set_atom)[0])

        if self.bSimple:
            # assign resonant atoms according to the "simple" rule: atoms with two sp2 neighbors (or heteroatoms) are resonant
            self.assign_uff_types_simplerule(tol, neighs, BOs, set_atom, set_bond)
            print("set_atom after simplerule:", np.where(~set_atom)[0])
        else:
            # assign resonant atoms according to the Huckel rule for aromaticity
            self.assign_uff_types_findrings(tol, neighs, BOs, BOs_int, set_atom, set_bond)
            print("set_atom after findrings:", np.where(~set_atom)[0])

        # assign the rest of atom types
        self.assign_uff_types_assignrest(neighs, BOs, BOs_int, set_atom, set_bond)
        print("set_atom after assignrest:", np.where(~set_atom)[0])

        # try to assign the rest of double bonds
        self.assign_uff_types_fixsaturation(neighs, BOs, BOs_int, set_bond, tol)

        # exception to avoid cumulenes
        if self.bCumulene:
            self.assign_uff_types_cumulene(neighs, BOs, BOs_int, set_bond)

        # manually change sp3 nitrogen and oxygen to "resonant" when they are bonded to an sp2 atom (conjugation)
        if self.bConj:
            self.assign_uff_types_conjugation(neighs, BOs)

        # some check before assigning parameters
        self.assign_uff_types_checks(neighs, BOs, BOs_int, set_atom, set_bond, tol)

        # exception for amide groups
        self.assign_uff_types_amide(neighs, BOs)

        # store bond orders
        self.mol.bond_orders = BOs

    def assign_uff_params(self):
        # This is a Python implementation of the C++ assignUFFparams logic
        # It will calculate and store bond, angle, dihedral, and inversion parameters
        neighs, _ = self.mol.make_neighs_bonds()
        self.assign_uff_params_vdws()
        self.assign_uff_params_bonds()
        self.assign_uff_params_angles(neighs)
        self.assign_uff_params_dihedrals(neighs)
        self.assign_uff_params_inversions(neighs)

    def get_arrays(self):
        # This will be the equivalent of the C++ toUFF method
        # It will return a dictionary of numpy arrays for the OpenCL runner
        uff_data = {
            "atype": np.array([a.type for a in self.builder_atoms], dtype=np.int32),
            "REQs": np.array([a.REQ for a in self.builder_atoms], dtype=np.float32),
            "bonAtoms": np.array([b.atoms for b in self.bonds], dtype=np.int32),
            "bonParams": np.array([[b.k, b.l0] for b in self.bonds], dtype=np.float32),
            "angAtoms": np.array([list(a.atoms) + [-1] for a in self.mol.angles], dtype=np.int32),
            "angParams": np.array([[a.k, a.C0, a.C1, a.C2, a.C3] for a in self.mol.angles], dtype=np.float32),
            "dihAtoms": np.array([d.atoms for d in self.mol.dihedrals], dtype=np.int32),
            "dihParams": np.array([[d.k, d.d, d.n] for d in self.mol.dihedrals], dtype=np.float32),
            "invAtoms": np.array([i.atoms for i in self.mol.inversions], dtype=np.int32),
            "invParams": np.array([[i.k, i.C0, i.C1, i.C2] for i in self.mol.inversions], dtype=np.float32),
            "neighs": self.neighs,
            "neighBs": self.neighBs,
        }
        return uff_data

    # --- Helper methods for assign_uff_types ---

    def assign_uff_types_trivial(self, neighs, BOs, BOs_int, set_atom, set_bond):
        """
        Assigns UFF atom types and bond orders for simple, unambiguous cases.

        This function handles atoms and bonds that can be typed without ambiguity,
        such as hydrogens, sp3 carbons, and nitrile groups. It serves as the
        first pass in the typing process, resolving the easy cases before more
        complex algorithms are used for resonant or ambiguous structures.
        """
        for ia, ai in enumerate(self.builder_atoms):
            ci = self.mol.confs[ai.iconf]
            if self.params.atypes[ai.type].name.startswith('H'):
                ai.type = self.params.getAtomType("H_")
                set_atom[ia] = True
                ja = neighs[ia, 0]
                ib = self.mol.get_bond_by_atoms(ia, ja)
                if ib is not None:
                    BOs[ib] = 1.0
                    BOs_int[ib] = 1
                    set_bond[ib] = True
            else:
                if self.params.atypes[ai.type].name.startswith('C') and ci.nbond == 4:
                    ai.type = self.params.getAtomType("C_3")
                    set_atom[ia] = True
                    for i in range(4):
                        ja = neighs[ia, i]
                        if ja < 0: continue
                        ib = self.mol.get_bond_by_atoms(ia, ja)
                        if ib is not None:
                            BOs[ib] = 1.0
                            BOs_int[ib] = 1
                            set_bond[ib] = True
                elif self.params.atypes[ai.type].name.startswith('N') and ci.nbond == 3:
                    for i in range(3):
                        ja = neighs[ia, i]
                        if ja < 0: continue
                        ib = self.mol.get_bond_by_atoms(ia, ja)
                        if ib is not None:
                            BOs_int[ib] = 1
                elif self.params.atypes[ai.type].name.startswith('N') and ci.nbond == 1:
                    ai.type = self.params.getAtomType("N_1")
                    set_atom[ia] = True
                    ja = neighs[ia, 0]
                    ib = self.mol.get_bond_by_atoms(ia, ja)
                    if ib is not None:
                        BOs[ib] = 3.0
                        BOs_int[ib] = 3
                        set_bond[ib] = True

                    aj = self.builder_atoms[ja]
                    if self.params.atypes[aj.type].name.startswith('C'):
                        aj.type = self.params.getAtomType("C_1")
                        set_atom[ja] = True
                        cj = self.mol.confs[aj.iconf]
                        for j in range(cj.nbond):
                            ka = neighs[ja, j]
                            if ka < 0: continue
                            if ka == ia: continue
                            jb = self.mol.get_bond_by_atoms(ja, ka)
                            if jb is not None:
                                BOs[jb] = 1.0
                                BOs_int[jb] = 1
                                set_bond[jb] = True
                elif self.params.atypes[ai.type].name.startswith('O') and ci.nbond == 2:
                    for i in range(2):
                        ja = neighs[ia, i]
                        if ja < 0: continue
                        ib = self.mol.get_bond_by_atoms(ia, ja)
                        if ib is not None:
                            BOs_int[ib] = 1
                elif self.params.atypes[ai.type].name.startswith('O') and ci.nbond == 1:
                    ja = neighs[ia, 0]
                    ib = self.mol.get_bond_by_atoms(ia, ja)
                    if ib is not None:
                        BOs[ib] = 2.0
                        BOs_int[ib] = 2

    def assign_uff_types_nitro(self, neighs, BOs, BOs_int, set_atom, set_bond):
        """
        Identifies and assigns UFF types for nitro groups.

        This function specifically searches for nitrogen atoms bonded to two
        oxygen atoms (a nitro group). It assigns resonant types (N_R, O_R)
        to these atoms and sets the N-O bonds to a resonant bond order of 1.5.
        """
        for ia, ai in enumerate(self.builder_atoms):
            if self.params.atypes[ai.type].name.startswith('H'):
                continue
            ci = self.mol.confs[ai.iconf]
            if self.params.atypes[ai.type].name.startswith('N') and ci.nbond == 3:
                n_oxy = 0
                for i in range(3):
                    ja = neighs[ia, i]
                    if ja < 0: continue
                    aj = self.builder_atoms[ja]
                    if self.params.atypes[aj.type].name.startswith('O'):
                        n_oxy += 1

                if n_oxy == 2:
                    ai.type = self.params.getAtomType("N_R")
                    set_atom[ia] = True
                    for i in range(3):
                        ja = neighs[ia, i]
                        if ja < 0: continue
                        aj = self.builder_atoms[ja]
                        ib = self.mol.get_bond_by_atoms(ia, ja)
                        if ib is None: continue

                        if self.params.atypes[aj.type].name.startswith('O'):
                            aj.type = self.params.getAtomType("O_R")
                            set_atom[ja] = True
                            BOs[ib] = 1.5
                            set_bond[ib] = True
                        else:
                            BOs[ib] = 1.0
                            set_bond[ib] = True

    def _assign_uff_types_checkall(self, neighs, BOs_int):
        while True:
            changed = False
            for ia, atom in enumerate(self.builder_atoms):
                if self.params.atypes[atom.type].name.startswith('H'):
                    continue

                conf = self.mol.confs[atom.iconf]

                found_unset = False
                for i in range(conf.nbond):
                    ja = neighs[ia, i]
                    if ja < 0: continue
                    ib = self.mol.get_bond_by_atoms(ia, ja)
                    if ib is not None and BOs_int[ib] < 0:
                        found_unset = True
                        break
                if not found_unset:
                    continue

                nset = 0
                val = 0
                for i in range(conf.nbond):
                    ja = neighs[ia, i]
                    if ja < 0: continue
                    ib = self.mol.get_bond_by_atoms(ia, ja)
                    if ib is not None and BOs_int[ib] > 0:
                        nset += 1
                        val += BOs_int[ib]

                atom_valence = self.params.atypes[atom.type].valence
                if (nset == conf.nbond - 1) or (conf.nbond - nset == atom_valence - val):
                    changed = True
                    if nset == conf.nbond - 1:
                        for i in range(conf.nbond):
                            ja = neighs[ia, i]
                            if ja < 0: continue
                            ib = self.mol.get_bond_by_atoms(ia, ja)
                            if ib is not None and BOs_int[ib] < 0:
                                BOs_int[ib] = atom_valence - val
                                break
                    elif (conf.nbond - nset) == (atom_valence - val):
                        for i in range(conf.nbond):
                            ja = neighs[ia, i]
                            if ja < 0: continue
                            ib = self.mol.get_bond_by_atoms(ia, ja)
                            if ib is not None and BOs_int[ib] < 0:
                                BOs_int[ib] = 1
            if not changed:
                break

        return not np.any(BOs_int < 0)

    def _assign_uff_types_checkatom(self, ia, neighs, BOs_int):
        atom = self.builder_atoms[ia]
        conf = self.mol.confs[atom.iconf]
        for i in range(conf.nbond):
            ja = neighs[ia, i]
            if ja < 0: continue
            ib = self.mol.get_bond_by_atoms(ia, ja)
            if ib is not None and BOs_int[ib] < 0:
                return False
        return True

    def assign_uff_types_treewalk(self, neighs, BOs_int):
        """
        Resolves ambiguous resonance structures using a tree-like search.

        When the initial typing rules fail to assign all bond orders, this function
        attempts to find a valid Lewis structure by exploring different
        possibilities for double bonds. It iteratively tries assigning a double
        bond to an unassigned bond and then checks if this choice leads to a
        fully determined, valid structure using the `_assign_uff_types_checkall`
        helper. This is a brute-force approach to find a valid resonance form.
        """
        if self._assign_uff_types_checkall(neighs, BOs_int):
            return

        bos_int_tmp0 = BOs_int.copy()
        natoms = len(self.builder_atoms)

        for ia1 in range(natoms):
            a1 = self.builder_atoms[ia1]
            if self.params.atypes[a1.type].name.startswith('H'): continue
            if self._assign_uff_types_checkatom(ia1, neighs, bos_int_tmp0): continue

            c1 = self.mol.confs[a1.iconf]
            for i1 in range(c1.nbond):
                ja1 = neighs[ia1, i1]
                if ja1 < 0: continue
                ib1 = self.mol.get_bond_by_atoms(ia1, ja1)
                if ib1 is None or bos_int_tmp0[ib1] > 0: continue

                bos_int_tmp1 = bos_int_tmp0.copy()
                bos_int_tmp1[ib1] = 2
                if self._assign_uff_types_checkall(neighs, bos_int_tmp1):
                    np.copyto(BOs_int, bos_int_tmp1)
                    return

                for ia2 in range(natoms):
                    a2 = self.builder_atoms[ia2]
                    if self.params.atypes[a2.type].name.startswith('H'): continue
                    if self._assign_uff_types_checkatom(ia2, neighs, bos_int_tmp1): continue

                    c2 = self.mol.confs[a2.iconf]
                    for i2 in range(c2.nbond):
                        ja2 = neighs[ia2, i2]
                        if ja2 < 0: continue
                        ib2 = self.mol.get_bond_by_atoms(ia2, ja2)
                        if ib2 is None or bos_int_tmp1[ib2] > 0: continue

                        bos_int_tmp2 = bos_int_tmp1.copy()
                        bos_int_tmp2[ib2] = 2
                        if self._assign_uff_types_checkall(neighs, bos_int_tmp2):
                            np.copyto(BOs_int, bos_int_tmp2)
                            return

                        for ia3 in range(natoms):
                            a3 = self.builder_atoms[ia3]
                            if self.params.atypes[a3.type].name.startswith('H'): continue
                            if self._assign_uff_types_checkatom(ia3, neighs, bos_int_tmp2): continue

                            c3 = self.mol.confs[a3.iconf]
                            for i3 in range(c3.nbond):
                                ja3 = neighs[ia3, i3]
                                if ja3 < 0: continue
                                ib3 = self.mol.get_bond_by_atoms(ia3, ja3)
                                if ib3 is None or bos_int_tmp2[ib3] > 0: continue

                                bos_int_tmp3 = bos_int_tmp2.copy()
                                bos_int_tmp3[ib3] = 2
                                if self._assign_uff_types_checkall(neighs, bos_int_tmp3):
                                    np.copyto(BOs_int, bos_int_tmp3)
                                    return

        raise Exception("ERROR: UFF treewalk failed to find a valid resonance structure.")

    def assign_uff_types_simplerule(self, tol, neighs, BOs, set_atom, set_bond):
        """
        Assigns resonant types based on a simple neighbor-based heuristic.

        This function identifies atoms that are likely part of a resonant system
        by checking their local environment. The rule is that any sp2-like atom
        (C with 3 neighbors, N with >1, O with 2) that is bonded to at least
        two other sp2-like atoms is considered resonant. This is a fast,
        heuristic-based alternative to explicit ring finding for identifying
        conjugated systems.
        """
        for ia, atom in enumerate(self.builder_atoms):
            if self.params.atypes[atom.type].name.startswith('H'):
                continue

            conf = self.mol.confs[atom.iconf]
            atom_type_name = self.params.atypes[atom.type].name

            is_candidate = ( (atom_type_name.startswith('C') and conf.nbond == 3) or
                             (atom_type_name.startswith('N') and conf.nbond > 1) or
                             (atom_type_name.startswith('O') and conf.nbond == 2) )

            if is_candidate:
                n_neighbors = 0
                for i in range(conf.nbond):
                    ja = neighs[ia, i]
                    if ja < 0: continue
                    neigh_atom = self.builder_atoms[ja]
                    if self.params.atypes[neigh_atom.type].name.startswith('H'):
                        continue

                    neigh_conf = self.mol.confs[neigh_atom.iconf]
                    neigh_type_name = self.params.atypes[neigh_atom.type].name

                    is_neigh_candidate = ( (neigh_type_name.startswith('C') and neigh_conf.nbond == 3) or
                                           (neigh_type_name.startswith('N') and neigh_conf.nbond > 1) or
                                           (neigh_type_name.startswith('O') and neigh_conf.nbond == 2) )
                    if is_neigh_candidate:
                        n_neighbors += 1

                if n_neighbors > 1:
                    if set_atom[ia]:
                        if not self.params.atypes[atom.type].name.endswith('R'):
                            print(f"WARNING SIMPLERULE: atom {ia+1} would be set to resonant but it has already a type of {self.params.atypes[atom.type].name}")
                        continue

                    if atom_type_name.startswith('C'):
                        atom.type = self.params.getAtomType("C_R")
                    elif atom_type_name.startswith('N'):
                        atom.type = self.params.getAtomType("N_R")
                    elif atom_type_name.startswith('O'):
                        atom.type = self.params.getAtomType("O_R")
                    else:
                        raise ValueError(f"ERROR SIMPLERULE: atom {ia} type {atom.type} not recognized")

                    set_atom[ia] = True

                    for i in range(conf.nbond):
                        ja = neighs[ia, i]
                        if ja < 0: continue
                        neigh_atom = self.builder_atoms[ja]
                        if self.params.atypes[neigh_atom.type].name.startswith('H'):
                            continue

                        neigh_conf = self.mol.confs[neigh_atom.iconf]
                        neigh_type_name = self.params.atypes[neigh_atom.type].name

                        is_neigh_candidate = ( (neigh_type_name.startswith('C') and neigh_conf.nbond == 3) or
                                               (neigh_type_name.startswith('N') and neigh_conf.nbond > 1) or
                                               (neigh_type_name.startswith('O') and neigh_conf.nbond == 2) )

                        if is_neigh_candidate:
                            ib = self.mol.get_bond_by_atoms(ia, ja)
                            if ib is None: continue
                            if set_bond[ib]:
                                if abs(BOs[ib] - 1.5) > tol:
                                    print(f"WARNING SIMPLERULE: bond {ib+1} between atoms {ia+1} and {ja+1} would be set to 1.5 but it has already a bond order of {BOs[ib]}")
                                continue
                            BOs[ib] = 1.5
                            set_bond[ib] = True

    def _assign_uff_types_checkaroma(self, n, nb, BOs, BOs_int, tol):
        for i in range(1, n + 1):
            atom = self.builder_atoms[nb[i]]
            conf = self.mol.confs[atom.iconf]
            if self.params.atypes[atom.type].name.startswith('C') and conf.nbond != 3:
                return False

        npi = 0
        for i in range(1, n + 1):
            atom = self.builder_atoms[nb[i]]
            conf = self.mol.confs[atom.iconf]
            atom_type_name = self.params.atypes[atom.type].name
            if (atom_type_name.startswith('N') and conf.nbond == 3) or \
               (atom_type_name.startswith('O') and conf.nbond == 2):
                npi += 2

        npi_save = npi

        for i in range(1, n + 1):
            ib1 = self.mol.get_bond_by_atoms(nb[i], nb[i % n + 1])
            ib2 = self.mol.get_bond_by_atoms(nb[i-1 if i > 1 else n], nb[i])
            if (ib1 is not None and BOs_int[ib1] == 2) or \
               (ib2 is not None and BOs_int[ib2] == 2):
                npi += 1

        if npi == 6:
            return True

        npi = npi_save

        for i in range(1, n + 1):
            ib1 = self.mol.get_bond_by_atoms(nb[i], nb[i % n + 1])
            ib2 = self.mol.get_bond_by_atoms(nb[i-1 if i > 1 else n], nb[i])
            if (ib1 is not None and (BOs_int[ib1] == 2 or abs(BOs[ib1] - 1.5) < tol)) or \
               (ib2 is not None and (BOs_int[ib2] == 2 or abs(BOs[ib2] - 1.5) < tol)):
                npi += 1

        return npi == 6

    def _assign_uff_types_setaroma(self, n, nb, neighs, BOs, set_atom, set_bond, tol):
        for i in range(1, n + 1):
            ia = nb[i]
            atom = self.builder_atoms[ia]
            if set_atom[ia]:
                print(f"WARNING SETAROMA: atom {self.params.atypes[atom.type].name[0]} {ia+1} would be set to resonant but it has already a type of {atom.type}")
                continue

            atom_type_name = self.params.atypes[atom.type].name
            if atom_type_name.startswith('C'):
                atom.type = self.params.getAtomType("C_R")
            elif atom_type_name.startswith('N'):
                atom.type = self.params.getAtomType("N_R")
            elif atom_type_name.startswith('O'):
                atom.type = self.params.getAtomType("O_R")
            else:
                raise ValueError(f"ERROR SETAROMA: atom {ia+1} type {atom.type} not recognized")
            set_atom[ia] = True

        for i in range(1, n + 1):
            ia = nb[i]
            ja = nb[i % n + 1]
            ib = self.mol.get_bond_by_atoms(ia, ja)
            if ib is None: continue
            if set_bond[ib] and abs(BOs[ib] - 1.5) > tol:
                a1 = self.builder_atoms[ia]
                a2 = self.builder_atoms[ja]
                print(f"WARNING SETAROMA: bond {ib+1} between atoms {self.params.atypes[a1.type].name[0]} {ia+1} and {self.params.atypes[a2.type].name[0]} {ja+1} would be set to 1.5 but it has already a bond order of {BOs[ib]}")
                continue
            BOs[ib] = 1.5
            set_bond[ib] = True

        for i in range(1, n + 1):
            ia = nb[i]
            atom = self.builder_atoms[ia]
            if self.params.atypes[atom.type].name.startswith('H'): continue
            conf = self.mol.confs[atom.iconf]
            for i_neigh in range(conf.nbond):
                ja = neighs[ia, i_neigh]
                if ja < 0: continue

                is_in_ring = False
                for j in range(1, n + 1):
                    if ja == nb[j]:
                        is_in_ring = True
                        break
                if is_in_ring: continue

                neigh_atom = self.builder_atoms[ja]
                if self.params.atypes[neigh_atom.type].name.startswith('H'): continue

                neigh_conf = self.mol.confs[neigh_atom.iconf]
                ib = self.mol.get_bond_by_atoms(ia, ja)
                if ib is None: continue

                if self.params.atypes[neigh_atom.type].name.startswith('O') and neigh_conf.nbond == 1:
                    if set_atom[ja] and not self.params.atypes[neigh_atom.type].name.endswith('R'):
                        print(f"WARNING SETAROMA: carbonyl atom {self.params.atypes[neigh_atom.type].name[0]} {ja+1} would be set to resonant but it has already a type of {neigh_atom.type}")
                        continue
                    if set_bond[ib] and abs(BOs[ib] - 1.5) > tol:
                        print(f"WARNING SETAROMA: bond {ib+1} between atoms {self.params.atypes[atom.type].name[0]} {ia+1} and {self.params.atypes[neigh_atom.type].name[0]} {ja+1} would be set to 1.5 but it has already a bond order of {BOs[ib]}")
                        continue

                    neigh_atom.type = self.params.getAtomType("O_R")
                    set_atom[ja] = True
                    BOs[ib] = 1.5
                    set_bond[ib] = True

    def assign_uff_types_findrings(self, tol, neighs, BOs, BOs_int, set_atom, set_bond):
        """
        Finds 5- and 6-membered rings and assigns aromatic types based on Huckel's rule.

        This function performs a brute-force search for small rings in the molecule.
        For each ring found, it uses the `_assign_uff_types_checkaroma` helper to
        determine if the ring is aromatic (i.e., has 6 pi electrons). If a ring
        is determined to be aromatic, the `_assign_uff_types_setaroma` helper is
        used to assign resonant types and bond orders to the atoms and bonds
        in the ring.
        """
        natoms = len(self.builder_atoms)
        explored = np.zeros(natoms, dtype=bool)

        while True:
            changed = False
            for ia in range(natoms):
                a1 = self.builder_atoms[ia]
                if self.params.atypes[a1.type].name.startswith('H'): continue
                c1 = self.mol.confs[a1.iconf]

                for i1 in range(c1.nbond):
                    nb2 = neighs[ia, i1]
                    if nb2 < 0: continue
                    a2 = self.builder_atoms[nb2]
                    if self.params.atypes[a2.type].name.startswith('H'): continue
                    c2 = self.mol.confs[a2.iconf]

                    for i2 in range(c2.nbond):
                        nb3 = neighs[nb2, i2]
                        if nb3 < 0 or nb3 == ia: continue
                        a3 = self.builder_atoms[nb3]
                        if self.params.atypes[a3.type].name.startswith('H'): continue
                        c3 = self.mol.confs[a3.iconf]

                        for i3 in range(c3.nbond):
                            nb4 = neighs[nb3, i3]
                            if nb4 < 0 or nb4 == nb2: continue
                            a4 = self.builder_atoms[nb4]
                            if self.params.atypes[a4.type].name.startswith('H'): continue
                            c4 = self.mol.confs[a4.iconf]

                            for i4 in range(c4.nbond):
                                nb5 = neighs[nb4, i4]
                                if nb5 < 0 or nb5 == nb3: continue
                                a5 = self.builder_atoms[nb5]
                                if self.params.atypes[a5.type].name.startswith('H'): continue
                                c5 = self.mol.confs[a5.iconf]

                                for i5 in range(c5.nbond):
                                    nb6 = neighs[nb5, i5]
                                    if nb6 < 0 or nb6 == nb4: continue

                                    if nb6 == ia:
                                        ring = [0, ia, nb2, nb3, nb4, nb5]
                                        if any(explored[i] for i in ring[1:]): continue

                                        if self._assign_uff_types_checkaroma(5, ring, BOs, BOs_int, tol):
                                            changed = True
                                            for i in ring[1:]: explored[i] = True
                                            self._assign_uff_types_setaroma(5, ring, neighs, BOs, set_atom, set_bond, tol)

                                    a6 = self.builder_atoms[nb6]
                                    if self.params.atypes[a6.type].name.startswith('H'): continue
                                    c6 = self.mol.confs[a6.iconf]
                                    for i6 in range(c6.nbond):
                                        nb7 = neighs[nb6, i6]
                                        if nb7 < 0 or nb7 == nb5: continue

                                        if nb7 == ia:
                                            ring = [0, ia, nb2, nb3, nb4, nb5, nb6]
                                            if any(explored[i] for i in ring[1:]): continue

                                            if self._assign_uff_types_checkaroma(6, ring, BOs, BOs_int, tol):
                                                changed = True
                                                for i in ring[1:]: explored[i] = True
                                                self._assign_uff_types_setaroma(6, ring, neighs, BOs, set_atom, set_bond, tol)
            if not changed:
                break

    def assign_uff_types_assignrest(self, neighs, BOs, BOs_int, set_atom, set_bond):
        """
        Assigns UFF types to any remaining un-typed atoms.

        After the more specific typing rules (trivial, nitro, resonance) have run,
        this function handles any atoms that are still without a UFF type. It
        assigns types based on the element and the number of bonds (e.g., a
        3-bonded carbon becomes C_2, a 2-bonded oxygen becomes O_3). This acts as a
        catch-all for common, non-resonant structures.
        """
        for ia, atom in enumerate(self.builder_atoms):
            if set_atom[ia]:
                continue
            if self.params.atypes[atom.type].name.startswith('H'):
                continue

            conf = self.mol.confs[atom.iconf]
            atom_type_name = self.params.atypes[atom.type].name

            if atom_type_name.startswith('C'):
                if conf.nbond == 2:
                    atom.type = self.params.getAtomType("C_1")
                    set_atom[ia] = True
                elif conf.nbond == 3:
                    atom.type = self.params.getAtomType("C_2")
                    set_atom[ia] = True
            elif atom_type_name.startswith('N'):
                if conf.nbond == 2:
                    atom.type = self.params.getAtomType("N_2")
                    set_atom[ia] = True
                elif conf.nbond == 3:
                    atom.type = self.params.getAtomType("N_3")
                    set_atom[ia] = True
                    for i in range(3):
                        ja = neighs[ia, i]
                        if ja < 0: continue
                        ib = self.mol.get_bond_by_atoms(ia, ja)
                        if ib is not None:
                            BOs[ib] = 1.0
                            BOs_int[ib] = 1
                            set_bond[ib] = True
            elif atom_type_name.startswith('O'):
                if conf.nbond == 1:
                    atom.type = self.params.getAtomType("O_2")
                    set_atom[ia] = True
                    ja = neighs[ia, 0]
                    if ja >= 0:
                        ib = self.mol.get_bond_by_atoms(ia, ja)
                        if ib is not None:
                            BOs[ib] = 2.0
                            BOs_int[ib] = 2
                            set_bond[ib] = True

                        neigh_atom = self.builder_atoms[ja]
                        if not self.params.atypes[neigh_atom.type].name.startswith('H'):
                            neigh_conf = self.mol.confs[neigh_atom.iconf]
                            if self.params.atypes[neigh_atom.type].name.startswith('C') and not set_atom[ja] and neigh_conf.nbond == 3:
                                neigh_atom.type = self.params.getAtomType("C_2")
                                set_atom[ja] = True

                elif conf.nbond == 2:
                    atom.type = self.params.getAtomType("O_3")
                    set_atom[ia] = True
                    for i in range(2):
                        ja = neighs[ia, i]
                        if ja < 0: continue
                        ib = self.mol.get_bond_by_atoms(ia, ja)
                        if ib is not None:
                            BOs[ib] = 1.0
                            BOs_int[ib] = 1
                            set_bond[ib] = True

    def assign_uff_types_fixsaturation(self, neighs, BOs, BOs_int, set_bond, tol):
        """
        Iteratively corrects bond orders to satisfy atom valences.

        This function runs in a loop to resolve any remaining un-set bond orders.
        It identifies atoms where some, but not all, bond orders are known. If an
        atom has only one remaining unknown bond, the order of that bond is set
        to whatever value is needed to satisfy the atom's total valence. If
        multiple bonds are unknown, it checks if setting them all to 1.0 would
        satisfy the valence. This process repeats until all bond orders are set.
        """
        while True:
            changed = False
            for ia, atom in enumerate(self.builder_atoms):
                if self.params.atypes[atom.type].name.startswith('H'):
                    continue

                conf = self.mol.confs[atom.iconf]

                all_bonds_set = True
                for i in range(conf.nbond):
                    ja = neighs[ia, i]
                    if ja < 0: continue
                    ib = self.mol.get_bond_by_atoms(ia, ja)
                    if ib is not None and BOs[ib] < 0.0:
                        all_bonds_set = False
                        break
                if all_bonds_set:
                    continue

                nset = 0
                val = 0.0
                for i in range(conf.nbond):
                    ja = neighs[ia, i]
                    if ja < 0: continue
                    ib = self.mol.get_bond_by_atoms(ia, ja)
                    if ib is not None and BOs[ib] > 0.0:
                        nset += 1
                        val += BOs_int[ib]

                remaining_valence = self.params.atypes[atom.type].valence - val

                if nset == conf.nbond - 1:
                    if abs(remaining_valence - round(remaining_valence)) > tol or abs(remaining_valence) < tol:
                        raise ValueError(f"ERROR FIXSATURATION: atom {self.params.atypes[atom.type].name[0]} {ia+1} would have a valence of {remaining_valence}")

                    changed = True
                    for i in range(conf.nbond):
                        ja = neighs[ia, i]
                        if ja < 0: continue
                        ib = self.mol.get_bond_by_atoms(ia, ja)
                        if ib is not None and BOs[ib] < 0.0:
                            BOs[ib] = remaining_valence
                            BOs_int[ib] = int(round(remaining_valence))
                            set_bond[ib] = True
                            break
                elif int(round(remaining_valence)) == conf.nbond - nset:
                    changed = True
                    for i in range(conf.nbond):
                        ja = neighs[ia, i]
                        if ja < 0: continue
                        ib = self.mol.get_bond_by_atoms(ia, ja)
                        if ib is not None and BOs[ib] < 0.0:
                            BOs[ib] = 1.0
                            BOs_int[ib] = 1
                            set_bond[ib] = True

            if not changed:
                break

    def assign_uff_types_cumulene(self, neighs, BOs, BOs_int, set_bond):
        """
        Applies an ad-hoc rule to handle specific cumulene-like structures.

        This function identifies a very specific bonding pattern of four sp and sp2
        carbons (C(sp2)-C(sp)=C(sp)-C(sp2)) that can be misinterpreted by the general
        typing rules. It corrects the bond orders in this specific case to avoid
        energetically unfavorable cumulene structures in the final force field.
        """
        for ia1, a1 in enumerate(self.builder_atoms):
            if not self.params.atypes[a1.type].name.startswith('C'): continue
            c1 = self.mol.confs[a1.iconf]
            if c1.nbond != 2: continue

            for i1 in range(2):
                ia2 = neighs[ia1, i1]
                if ia2 < 0: continue
                a2 = self.builder_atoms[ia2]
                if not self.params.atypes[a2.type].name.startswith('C'): continue
                c2 = self.mol.confs[a2.iconf]
                if c2.nbond != 2: continue

                ib1 = self.mol.get_bond_by_atoms(ia1, ia2)
                if ib1 is None: continue

                ia4 = neighs[ia1, 1-i1]
                if ia4 < 0: continue
                a4 = self.builder_atoms[ia4]
                if not self.params.atypes[a4.type].name.startswith('C'): continue
                c4 = self.mol.confs[a4.iconf]
                if c4.nbond != 3: continue

                ib4 = self.mol.get_bond_by_atoms(ia1, ia4)
                if ib4 is None: continue

                for i2 in range(2):
                    ia3 = neighs[ia2, i2]
                    if ia3 < 0 or ia3 == ia1: continue
                    a3 = self.builder_atoms[ia3]
                    if not self.params.atypes[a3.type].name.startswith('C'): continue
                    c3 = self.mol.confs[a3.iconf]
                    if c3.nbond != 3: continue

                    for i3 in range(3):
                        if neighs[ia3, i3] == ia4:
                            ib2 = self.mol.get_bond_by_atoms(ia2, ia3)
                            ib3 = self.mol.get_bond_by_atoms(ia3, ia4)
                            if ib2 is None or ib3 is None: continue

                            BOs[ib1] = 3.0; BOs_int[ib1] = 3; set_bond[ib1] = True
                            BOs[ib2] = 1.0; BOs_int[ib2] = 1; set_bond[ib2] = True
                            BOs[ib3] = 2.0; BOs_int[ib3] = 2; set_bond[ib3] = True
                            BOs[ib4] = 1.0; BOs_int[ib4] = 1; set_bond[ib4] = True
                            break

    def assign_uff_types_conjugation(self, neighs, BOs):
        """
        Handles conjugation between sp3 heteroatoms and resonant systems.

        This function identifies sp3-hybridized nitrogen or oxygen atoms that are
        directly bonded to an atom already identified as part of a resonant
        system. In such cases, the lone pair on the heteroatom can participate
        in conjugation. The function updates the heteroatom's type to resonant
        (N_R or O_R) and sets the connecting bond order to 1.5.
        """
        for ia, atom in enumerate(self.builder_atoms):
            if self.params.atypes[atom.type].name.startswith('H'):
                continue

            conf = self.mol.confs[atom.iconf]
            atom_type_name = self.params.atypes[atom.type].name

            if (atom_type_name.startswith('N') and conf.nbond == 3) or \
               (atom_type_name.startswith('O') and conf.nbond == 2):
                for i in range(conf.nbond):
                    ja = neighs[ia, i]
                    if ja < 0: continue

                    neigh_atom = self.builder_atoms[ja]
                    if self.params.atypes[neigh_atom.type].name.endswith('R'):
                        if atom_type_name.startswith('N'):
                            atom.type = self.params.getAtomType("N_R")
                        elif atom_type_name.startswith('O'):
                            atom.type = self.params.getAtomType("O_R")

                        ib = self.mol.get_bond_by_atoms(ia, ja)
                        if ib is not None:
                            BOs[ib] = 1.5

    def assign_uff_types_checks(self, neighs, BOs, BOs_int, set_atom, set_bond, tol):
        """
        Performs sanity checks on the final assigned types and bond orders.

        After all typing rules have been applied, this function verifies that the
        resulting structure is self-consistent. It checks that all atoms and bonds
        have been assigned a type, that resonant atoms have resonant bonds, that
        sp2 atoms have exactly one double bond, and that resonant bonds have
        valid integer bond orders in their limiting resonance structure.
        """
        if not np.all(set_atom):
            unassigned = np.where(~set_atom)[0]
            raise ValueError(f"ERROR CHECKS: atoms {unassigned} are not set")
        if not np.all(set_bond):
            unassigned = np.where(~set_bond)[0]
            raise ValueError(f"ERROR CHECKS: bonds {unassigned} are not set")
        if np.any(BOs < 0.0):
            unassigned = np.where(BOs < 0.0)[0]
            raise ValueError(f"ERROR CHECKS: order for bonds {unassigned} is not set")
        if np.any(BOs_int < 0):
            unassigned = np.where(BOs_int < 0)[0]
            raise ValueError(f"ERROR CHECKS: integer order for bonds {unassigned} is not set")

        for ia, atom in enumerate(self.builder_atoms):
            atom_type_name = self.params.atypes[atom.type].name
            if atom_type_name.endswith('R'):
                conf = self.mol.confs[atom.iconf]
                found = False
                for i in range(conf.nbond):
                    ja = neighs[ia, i]
                    if ja < 0: continue
                    ib = self.mol.get_bond_by_atoms(ia, ja)
                    if ib is not None and abs(BOs[ib] - 1.5) < tol:
                        found = True
                        break
                if not found:
                    raise ValueError(f"ERROR CHECKS: atom {atom_type_name[0]} {ia+1} is resonant but it has no 1.5 bonds")
            elif atom_type_name.endswith('2'):
                conf = self.mol.confs[atom.iconf]
                found_double = False
                for i in range(conf.nbond):
                    ja = neighs[ia, i]
                    if ja < 0: continue
                    ib = self.mol.get_bond_by_atoms(ia, ja)
                    if ib is not None and abs(BOs[ib] - 2.0) < tol:
                        if found_double:
                            raise ValueError(f"ERROR CHECKS: atom {atom_type_name[0]} {ia+1} is sp2 but it has more than one double bond")
                        found_double = True

                if not found_double:
                    found_resonant = False
                    for i in range(conf.nbond):
                        ja = neighs[ia, i]
                        if ja < 0: continue
                        ib = self.mol.get_bond_by_atoms(ia, ja)
                        if ib is not None and abs(BOs[ib] - 1.5) < tol:
                            found_resonant = True
                            break
                    if found_resonant:
                        print(f"WARNING CHECKS: atom {atom_type_name[0]} {ia+1} is sp2 and it has no double bonds, only 1.5 bonds")
                    else:
                        raise ValueError(f"ERROR CHECKS: atom {atom_type_name[0]} {ia+1} is sp2 but it has no double bonds")

        for ib, bo in enumerate(BOs):
            if abs(bo - 1.5) < tol:
                if BOs_int[ib] not in [1, 2]:
                    raise ValueError(f"ERROR CHECKS: bond {ib+1} is 1.5 but in the limit resonance structure is neither single nor double")

    def assign_uff_types_amide(self, neighs, BOs):
        """
        Identifies and assigns special types for amide groups.

        This function specifically looks for the amide functional group (a carbonyl
        carbon bonded to a nitrogen). Due to the special resonance of the amide
        bond, it assigns resonant types to the C, O, and N atoms and sets the
        C-O bond order to 1.5 and the C-N bond order to 1.41, which is a UFF
        convention for amides.
        """
        for ia, atom in enumerate(self.builder_atoms):
            if self.params.atypes[atom.type].name.startswith('H'):
                continue

            conf = self.mol.confs[atom.iconf]
            atom_type_name = self.params.atypes[atom.type].name

            if atom_type_name.startswith('C') and conf.nbond == 3:
                found_o = False
                for i in range(conf.nbond):
                    ja = neighs[ia, i]
                    if ja < 0: continue
                    neigh_atom = self.builder_atoms[ja]
                    if self.params.atypes[neigh_atom.type].name.startswith('H'): continue
                    neigh_conf = self.mol.confs[neigh_atom.iconf]
                    if self.params.atypes[neigh_atom.type].name.startswith('O') and neigh_conf.nbond == 1:
                        found_o = True
                        break
                if not found_o:
                    continue

                found_n = False
                for i in range(conf.nbond):
                    ja = neighs[ia, i]
                    if ja < 0: continue
                    neigh_atom = self.builder_atoms[ja]
                    if self.params.atypes[neigh_atom.type].name.startswith('H'): continue
                    neigh_conf = self.mol.confs[neigh_atom.iconf]
                    if self.params.atypes[neigh_atom.type].name.startswith('N') and neigh_conf.nbond == 3:
                        found_n = True
                        break
                if not found_n:
                    continue

                atom.type = self.params.getAtomType("C_R")
                for i in range(conf.nbond):
                    ja = neighs[ia, i]
                    if ja < 0: continue
                    neigh_atom = self.builder_atoms[ja]
                    if self.params.atypes[neigh_atom.type].name.startswith('H'): continue
                    neigh_conf = self.mol.confs[neigh_atom.iconf]

                    neigh_atom_type_name = self.params.atypes[neigh_atom.type].name
                    if neigh_atom_type_name.startswith('O') and neigh_conf.nbond == 1:
                        neigh_atom.type = self.params.getAtomType("O_R")
                        ib = self.mol.get_bond_by_atoms(ia, ja)
                        if ib is not None:
                            BOs[ib] = 1.5
                    elif neigh_atom_type_name.startswith('N') and neigh_conf.nbond == 3:
                        neigh_atom.type = self.params.getAtomType("N_R")
                        ib = self.mol.get_bond_by_atoms(ia, ja)
                        if ib is not None:
                            BOs[ib] = 1.41

    # --- Helper methods for assign_uff_params ---

    def assign_uff_params_vdws(self):
        for atom in self.builder_atoms:
            self.params.assignRE(atom.type, atom.REQ, True)

    def assign_uff_params_bonds(self):
        for i, bond in enumerate(self.bonds):
            bond.l0 = self.assign_uff_params_calcrij(i)
            ei = self.params.elementOfAtomType(self.builder_atoms[bond.atoms[0]].type)
            ej = self.params.elementOfAtomType(self.builder_atoms[bond.atoms[1]].type)
            bond.k = 0.5 * 28.7989689090648 * ei.Quff * ej.Quff / (bond.l0**3)

    def assign_uff_params_angles(self, neighs):
        self.mol.angles = []
        for j in range(len(self.builder_atoms)):
            aj = self.builder_atoms[j]
            tj = self.params.atypes[aj.type]
            ct = np.cos(np.deg2rad(tj.Ass))
            st2 = np.sin(np.deg2rad(tj.Ass))**2
            if self.params.atypes[aj.type].name[0] == 'H':
                continue
            cj = self.mol.confs[aj.iconf]
            for in1 in range(cj.nbond - 1):
                i = neighs[j, in1]
                ei = self.params.elementOfAtomType(self.builder_atoms[i].type)
                for in2 in range(in1 + 1, cj.nbond):
                    k = neighs[j, in2]
                    ek = self.params.elementOfAtomType(self.builder_atoms[k].type)
                    a = Angle()
                    a.atoms = (i, j, k)
                    a.bonds = (self.mol.get_bond_by_atoms(i, j), self.mol.get_bond_by_atoms(j, k))
                    rij = self.assign_uff_params_calcrij(a.bonds[0])
                    rjk = self.assign_uff_params_calcrij(a.bonds[1])
                    rik = np.sqrt(rij**2 + rjk**2 - 2.0 * rij * rjk * ct)
                    kappa = 28.7989689090648 * ei.Quff * ek.Quff / (rik**5) * (3.0 * rij * rjk * st2 - rik**2 * ct)
                    if self.params.atypes[aj.type].name[2] in ['1', '2', 'R']:
                        if self.params.atypes[aj.type].name[2] == '1':
                            a.k = kappa
                            a.C0 = 1.0
                            a.C1 = 1.0
                            a.C2 = 0.0
                            a.C3 = 0.0
                        elif self.params.atypes[aj.type].name[2] in ['2', 'R']:
                            a.k = kappa / 9.0
                            a.C0 = 1.0
                            a.C1 = 0.0
                            a.C2 = 0.0
                            a.C3 = -1.0
                    elif self.params.atypes[aj.type].name[2] == '3':
                        a.k = kappa
                        a.C2 = 1.0 / (4.0 * st2)
                        a.C1 = -4.0 * a.C2 * ct
                        a.C0 = a.C2 * (2.0 * ct**2 + 1.0)
                        a.C3 = 0.0
                    self.mol.angles.append(a)

    def assign_uff_params_dihedrals(self, neighs):
        self.mol.dihedrals = []
        for i1 in range(len(self.builder_atoms)):
            a1 = self.builder_atoms[i1]
            if self.params.atypes[a1.type].name[0] == 'H':
                n1 = 1
            else:
                c1 = self.mol.confs[a1.iconf]
                n1 = c1.nbond
            for in1 in range(n1):
                i2 = neighs[i1, in1]
                a2 = self.builder_atoms[i2]
                if self.params.atypes[a2.type].name[0] == 'H' or self.params.atypes[a2.type].name[2] == '1':
                    continue
                c2 = self.mol.confs[a2.iconf]
                for in2 in range(c2.nbond):
                    i3 = neighs[i2, in2]
                    if i3 != i1:
                        a3 = self.builder_atoms[i3]
                        if self.params.atypes[a3.type].name[0] == 'H' or self.params.atypes[a3.type].name[2] == '1':
                            continue
                        c3 = self.mol.confs[a3.iconf]
                        for in3 in range(c3.nbond):
                            i4 = neighs[i3, in3]
                            if i4 != i2 and i4 > i1:
                                a4 = self.builder_atoms[i4]
                                d = Dihedral()
                                d.atoms = (i1, i2, i3, i4)
                                d.bonds = (self.mol.get_bond_by_atoms(i1, i2), self.mol.get_bond_by_atoms(i2, i3), self.mol.get_bond_by_atoms(i3, i4))
                                e2 = self.params.elementOfAtomType(a2.type)
                                e3 = self.params.elementOfAtomType(a3.type)
                                if self.params.atypes[a2.type].name[2] == '3' and self.params.atypes[a3.type].name[2] == '3':
                                    d.k = np.sqrt(e2.Vuff * e3.Vuff)
                                    d.d = 1
                                    d.n = 3
                                    if self.params.atypes[a2.type].name[0] in ['O', 'S'] and self.params.atypes[a3.type].name[0] in ['O', 'S']:
                                        d.k = 4.1840 / 60.2214076 / 1.602176634
                                        d.n = 2
                                        d.k *= 2.0 if self.params.atypes[a2.type].name[0] == 'O' else 6.8
                                        d.k *= 2.0 if self.params.atypes[a3.type].name[0] == 'O' else 6.8
                                        d.k = np.sqrt(d.k)
                                elif (self.params.atypes[a2.type].name[2] == '3' and self.params.atypes[a3.type].name[2] in ['2', 'R']) or (self.params.atypes[a2.type].name[2] in ['2', 'R'] and self.params.atypes[a3.type].name[2] == '3'):
                                    d.k = 4.1840 / 60.2214076 / 1.602176634
                                    d.d = -1
                                    d.n = 6
                                    if (self.params.atypes[a2.type].name[2] == '3' and self.params.atypes[a2.type].name[0] in ['O', 'S']) or (self.params.atypes[a3.type].name[2] == '3' and self.params.atypes[a3.type].name[0] in ['O', 'S']):
                                        ib = self.mol.get_bond_by_atoms(i2, i3)
                                        b = self.bonds[ib]
                                        d.k = 5.0 * np.sqrt(e2.Uuff * e3.Uuff) * (1.0 + 4.18 * np.log(b.order))
                                        d.d = 1
                                        d.n = 2
                                    if self.params.atypes[a2.type].name[2] == '3':
                                        found = any(self.params.atypes[self.builder_atoms[neighs[i3, i]].type].name[2] in ['2', 'R'] for i in range(c3.nbond))
                                        if found:
                                            d.k = 2.0 * 4.1840 / 60.2214076 / 1.602176634
                                            d.d = 1
                                            d.n = 3
                                    elif self.params.atypes[a3.type].name[2] == '3':
                                        found = any(self.params.atypes[self.builder_atoms[neighs[i2, i]].type].name[2] in ['2', 'R'] for i in range(c2.nbond))
                                        if found:
                                            d.k = 2.0 * 4.1840 / 60.2214076 / 1.602176634
                                            d.d = 1
                                            d.n = 3
                                elif self.params.atypes[a2.type].name[2] in ['2', 'R'] and self.params.atypes[a3.type].name[2] in ['2', 'R']:
                                    ib = self.mol.get_bond_by_atoms(i2, i3)
                                    b = self.bonds[ib]
                                    d.k = 5.0 * np.sqrt(e2.Uuff * e3.Uuff) * (1.0 + 4.18 * np.log(b.order))
                                    d.d = -1
                                    d.n = 2
                                else:
                                    raise ValueError(f"Dihedral case not found for atoms {self.params.atypes[a1.type].name} {self.params.atypes[a2.type].name} {self.params.atypes[a3.type].name} {self.params.atypes[a4.type].name}")
                                d.k = 0.5 * d.k / ((c2.nbond - 1) * (c3.nbond - 1))
                                self.mol.dihedrals.append(d)

    def assign_uff_params_inversions(self, neighs):
        self.mol.inversions = []
        for i1 in range(len(self.builder_atoms)):
            a1 = self.builder_atoms[i1]
            if self.params.atypes[a1.type].name[0] == 'H':
                continue
            c1 = self.mol.confs[a1.iconf]
            if c1.nbond != 3:
                continue
            i2 = neighs[i1, 0]
            i3 = neighs[i1, 1]
            i4 = neighs[i1, 2]
            self.assign_uff_params_assigninversion(i1, i2, i3, i4)
            self.assign_uff_params_assigninversion(i1, i4, i2, i3)
            self.assign_uff_params_assigninversion(i1, i3, i4, i2)

    def assign_uff_params_assigninversion(self, i1, i2, i3, i4):
        a1 = self.builder_atoms[i1]
        a2 = self.builder_atoms[i2]
        a3 = self.builder_atoms[i3]
        a4 = self.builder_atoms[i4]

        i = Inversion()
        i.atoms = (i1, i2, i3, i4)
        i.bonds = (self.mol.get_bond_by_atoms(i1, i2), self.mol.get_bond_by_atoms(i1, i3), self.mol.get_bond_by_atoms(i1, i4))
        if self.params.atypes[a1.type].name[0] == 'C' and self.params.atypes[a1.type].name[2] in ['2', 'R']:
            if any(self.params.atypes[self.builder_atoms[at].type].name[0] == 'O' and self.params.atypes[self.builder_atoms[at].type].name[2] == '2' for at in [i2, i3, i4]):
                i.k = 50.0 * 4.1840 / 60.2214076 / 1.602176634
            else:
                i.k = 6.0 * 4.1840 / 60.2214076 / 1.602176634
            i.C0 = 1.0
            i.C1 = -1.0
            i.C2 = 0.0
        elif self.params.atypes[a1.type].name[0] == 'N' and self.params.atypes[a1.type].name[2] in ['2', 'R']:
            i.k = 6.0 * 4.1840 / 60.2214076 / 1.602176634
            i.C0 = 1.0
            i.C1 = -1.0
            i.C2 = 0.0
        elif self.params.atypes[a1.type].name[0] == 'N' and self.params.atypes[a1.type].name[2] == '3':
            i.k = 0.0
            i.C0 = 0.0
            i.C1 = 0.0
            i.C2 = 0.0
        elif self.params.atypes[a1.type].name[0] in ['P', '3']:
            omega0 = np.deg2rad(84.4339)
            i.C0 = 4.0 * np.cos(omega0)**2 - np.cos(2.0 * omega0)
            i.C1 = -4.0 * np.cos(omega0)
            i.C2 = 1.0
            i.k = 22.0 * 4.1840 / 60.2214076 / 1.602176634 / (i.C0 + i.C1 + i.C2)
        else:
            raise ValueError(f"Inversion case not found for atoms {self.params.atypes[a1.type].name} {self.params.atypes[a2.type].name} {self.params.atypes[a3.type].name} {self.params.atypes[a4.type].name}")
        i.k /= 3.0
        self.mol.inversions.append(i)

    def assign_uff_params_calcrij(self, ib):
        bond = self.bonds[ib]
        ti = self.params.atypes[self.builder_atoms[bond.atoms[0]].type]
        tj = self.params.atypes[self.builder_atoms[bond.atoms[1]].type]
        ei = self.params.elementOfAtomType(self.builder_atoms[bond.atoms[0]].type)
        ej = self.params.elementOfAtomType(self.builder_atoms[bond.atoms[1]].type)
        rbo = -0.1332 * (ti.Ruff + tj.Ruff) * np.log(bond.order)
        ren = ti.Ruff * tj.Ruff * (np.sqrt(-ei.Eaff) - np.sqrt(-ej.Eaff))**2 / (-ei.Eaff * ti.Ruff - ej.Eaff * tj.Ruff)
        rij = ti.Ruff + tj.Ruff + rbo - ren
        return rij
