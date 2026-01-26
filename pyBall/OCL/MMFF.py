import numpy as np
import math
import os
from . import MMparams

# Constants (Define these based on your simulation requirements)
Lepair = 1.5    # Example bond length for electron pairs
Kepair = 100.0  # Example bond stiffness for electron pairs
deg2rad = np.pi / 180.0

verbosity = 0

PI_HINTS = {
    # Carbon families
    "C_2": 1, "C_R": 1, "C_ene": 1, "C_ald": 1, "C_COO": 1,
    "C_1": 2, "C_yne": 2,
    # Oxygen families
    "O_2": 1, "O_R": 1, "O_ald": 1, "O_pCOO": 1,
    # Nitrogen families
    "N_2": 1, "N_R": 1, "N_1": 2,
}

def initAtomProperties(mol, atom_types, capping_atoms={'H'}, bPrint=True):
    """Initialize atom properties (npi, nep, isNode) based on atom types.
    
    Args:
        mol: AtomicSystem object
        atom_types: Dictionary of AtomType objects
        capping_atoms: Set of element symbols for capping atoms (default: {'H'})
    """
    natoms = mol.natoms
    npi_list = np.zeros(natoms, dtype=np.int32)
    nep_list = np.zeros(natoms, dtype=np.int32)
    isNode   = np.ones(natoms,  dtype=np.int32)  # Default to node
    for i in range(natoms):
        atom_name = mol.enames[i]
        if atom_name in atom_types:
            at = atom_types[atom_name]
            npi_list[i] = at.npi
            nep_list[i] = at.nepair
            if bPrint:
                if npi_list[i] > 0:
                    print(f"DEBUG initAtomProperties: atom {i} {atom_name} npi={npi_list[i]} nep={nep_list[i]}")
                else:
                    hinted = PI_HINTS.get(atom_name, 0)
                    if hinted > 0:
                        npi_list[i] = hinted
                        print(f"DEBUG initAtomProperties: atom {i} {atom_name} using hint npi={hinted}")
        # Mark as capping if element is in capping set
        if atom_name in capping_atoms:
            isNode[i] = 0
        elif atom_name in atom_types:
            at = atom_types[atom_name]
            if at.element_name in capping_atoms:
                isNode[i] = 0
            elif bPrint:
                print(f"DEBUG initAtomProperties: atom {i} {atom_name} element={at.element_name} marked node")
        elif bPrint:
            print(f"DEBUG initAtomProperties: atom {i} {atom_name} not in atom_types -> treated as node")
    return npi_list, nep_list, isNode

class MMFF:
    """
    Represents the MMFF (Merck Molecular Force Field) parameters and configurations.
    """
    def __init__(self, bTorsion=False, verbosity=1, *, reorder_nodes_first: bool = True):
        """
        Initializes the MMFF instance.

        Parameters:
        - bTorsion (bool): Flag to indicate if torsions are considered.
        - verbosity (int): Level of verbosity for logging.
        """
        self.bTorsion = bTorsion
        self.natoms = 0
        self.nnode  = 0
        self.ncap   = 0
        self.nvecs  = 0
        self.ntors  = 0
        self.nDOFs  = 0

        # Load element_types and atom_types if needed for UFF calculations
        if not hasattr(self, 'element_types') or not hasattr(self, 'atom_types'):
            base_path = os.path.dirname(os.path.abspath(__file__))
            data_path = os.path.join(base_path, "../../cpp/common_resources/")
            self.element_types = MMparams.read_element_types(os.path.join(data_path, 'ElementTypes.dat'))
            self.atom_types = MMparams.read_atom_types(os.path.join(data_path, 'AtomTypes.dat'), self.element_types)

        # Initialize arrays
        self.apos      = None       # [nvecs,  4 ] Positions
        self.fapos     = None       # [nvecs,  4 ] Forces
        self.atypes    = None       # [natoms    ]  Atom types 
        self.neighs    = None       # [natoms, 4 ]  Neighbor indices 
        self.neighCell = None       # [natoms, 4 ] Neighbor cell indices
        self.REQs      = None       # [natoms, 4 ] Non-covalent parameters
        self.apars     = None       # [nnode,  4 ] Angle parameters
        self.bLs       = None       # [nnode,  4 ] Bond lengths
        self.bKs       = None       # [nnode,  4 ] Bond stiffness
        self.Ksp       = None       # [nnode,  4 ] Pi-sigma stiffness
        self.Kpp       = None       # [nnode,  4 ] Pi-pi stiffness
        self.angles    = None       # [nnode, 6, 3 ] Angles between bonds
        # self.tors2atom = None     # [ntors]        Torsion atom indices
        # self.torsParams = None    # [ntors, 4    ] Torsion parameters
        # self.constr = None        # [natoms, 4   ] Constraints
        # self.constrK = None       # [natoms, 4   ] Constraint stiffness
        self.invLvec    = None      # [nSystems, 9 ] Inverse lattice vectors
        self.pbc_shifts = None      # [nSystems, 3 ] PBC shifts
        self.nPBC       = None
        self.npbc       = None
        self.capping_atoms = {'H'}
        self.reorder_nodes_first = bool(reorder_nodes_first)
        self.dt         = 0.01
        self.damp       = 0.1
        self.Flimit     = 10.0
        # Pi-orbitals and electron pairs
        self.pipos = None         # [natoms, 3] Pi-orbital orientations

        # Optional debug hook for topology generation parity checks.
        # If set to a callable(str), it will receive line-by-line debug records.
        self.debug_topo = None

    def _dbg_topo(self, line: str):
        cb = getattr(self, 'debug_topo', None)
        if cb is None:
            return
        cb(str(line))

    def realloc(self, nnode, ncap, ntors=0, nPBC=(0,0,0) ):
        """
        Reallocates memory for MMFF parameters based on the system size.

        Parameters:
        - nnode (int): Number of nodes (configurations).
        - ncap (int): Number of capping atoms.
        - ntors (int): Number of torsions.
        """
        self.npbc   = (nPBC[0]*2+1)*(nPBC[1]*2+1)*(nPBC[2]*2+1)
        self.nPBC   = nPBC
        self.nnode  = nnode
        self.ncap   = ncap
        self.natoms = nnode + ncap
        self.nvecs  = self.natoms + nnode  # Including pi-orbitals
        self.ntors  = ntors
        self.nDOFs  = self.nvecs * 3
        if verbosity > 0:
            print(f"MMFF::realloc() natoms {self.natoms} nnode {self.nnode} ncap {self.ncap} nvecs {self.nvecs} ntors {self.ntors} nDOFs {self.nDOFs} nPBC {self.nPBC} npbc {self.npbc}")

        # Initialize arrays with default values
        self.apos       = np.full((self.nvecs, 4), 0.0, dtype=np.float32)
        self.fapos      = np.full((self.nvecs, 4), 0.0, dtype=np.float32)
        self.atypes     = np.full( self.natoms,    -1,  dtype=np.int32)
        self.neighs     = np.full((self.natoms,4), -1,  dtype=np.int32)
        self.neighCell  = np.full((self.natoms,4),  0,  dtype=np.int32)
        self.REQs       = np.full((self.natoms,4), 0.0, dtype=np.float32)
        self.apars      = np.full((self.nnode, 4), 0.0, dtype=np.float32)
        self.bLs        = np.full((self.nnode, 4), 0.0, dtype=np.float32)
        self.bKs        = np.full((self.nnode, 4), 0.0, dtype=np.float32)
        self.Ksp        = np.full((self.nnode, 4), 0.0, dtype=np.float32)
        self.Kpp        = np.full((self.nnode, 4), 0.0, dtype=np.float32)
        self.angles     = np.full((self.nnode, 6, 3), 0.0, dtype=np.float32)
        # self.tors2atom  = np.full( self.ntors,      0.0,dtype=np.int32)
        # self.torsParams = np.full((self.ntors, 4),  0.0, dtype=np.float32)
        # self.constr     = np.full((self.natoms, 4), 0.0, dtype=np.float32)
        # self.constrK    = np.full((self.natoms, 4), 0.0, dtype=np.float32)
        self.invLvec    = np.full((3,3), 0.0,  dtype=np.float32)  # Assuming single system; modify as needed
        self.pbc_shifts = np.full((self.npbc, 3), 0.0,  dtype=np.float32)
        self.pipos = np.zeros((self.natoms, 3), dtype=np.float32)

    def make_back_neighs(self, b_cap_neighs=True):
        """Populate back-neighbor indices similar to `MMFFsp3_loc::makeBackNeighs()`."""
        self.back_neighs = np.full((self.nvecs, 4), -1, dtype=np.int32)
        # Node atoms contribute entries for their sigma bonds and associated pi recoils.
        pi_offset = self.nnode * 4
        for ia in range(self.nnode):
            for ib_idx in range(4):
                ja = int(self.neighs[ia, ib_idx])
                if ja < 0 or ja >= self.natoms: continue
                packed = ia * 4 + ib_idx
                empty = np.where(self.back_neighs[ja] < 0)[0]
                if empty.size == 0: raise ValueError(f"Atom {ja} has >4 back-neighbors when adding {ia}")
                self.back_neighs[ja, empty[0]] = packed
                if ja < self.nnode:
                    pi_idx = self.natoms + ja
                    empty_pi = np.where(self.back_neighs[pi_idx] < 0)[0]
                    if empty_pi.size == 0: raise ValueError(f"Pi orbital {ja} has >4 back-neighbors when adding {ia}")
                    self.back_neighs[pi_idx, empty_pi[0]] = packed + pi_offset
        if b_cap_neighs:
            for ia in range(self.nnode, self.natoms):
                if self.back_neighs[ia, 0] < 0: raise ValueError(f"Capping atom {ia} missing back-neighbor")
                # mirror CPU behavior: set first neighbor to bonded node index
                self.neighs[ia, :] = -1
                self.neighs[ia, 0] = self.back_neighs[ia, 0] // 4

    def back_neighs_as_atoms(self):
        if not hasattr(self, 'back_neighs'): return None
        arr = np.full_like(self.back_neighs, -1)
        mask = self.back_neighs >= 0
        arr[mask] = self.back_neighs[mask] // 4
        return arr

    def _ensure_node_first(self, mol, atom_types):
        natoms = mol.apos.shape[0]
        perm_existing = getattr(mol, "perm_nodes_first", None)
        if perm_existing is not None and len(perm_existing) == natoms:
            npi_list = getattr(mol, "npi_list", None)
            nep_list = getattr(mol, "nep_list", None)
            is_node  = getattr(mol, "isNode", None)
            if npi_list is None or len(npi_list) != natoms or nep_list is None or len(nep_list) != natoms or is_node is None or len(is_node) != natoms:
                npi_list, nep_list, is_node = initAtomProperties(mol, atom_types, self.capping_atoms, bPrint=False)
                mol.npi_list = list(npi_list)
                mol.nep_list = list(nep_list)
                mol.isNode   = list(is_node)
            if mol.ngs is None: mol.neighs()
            return list(mol.npi_list), list(mol.nep_list), list(mol.isNode)

        npi_raw, nep_raw, is_node_raw = initAtomProperties(mol, atom_types, self.capping_atoms, bPrint=False)
        node_indices = [i for i, flag in enumerate(is_node_raw) if flag > 0]
        cap_indices  = [i for i, flag in enumerate(is_node_raw) if flag <= 0]
        perm         = node_indices + cap_indices

        if perm == list(range(natoms)):
            mol.perm_nodes_first = perm
            mol.perm_inverse = list(range(natoms))
            mol.npi_list     = list(int(x) for x in npi_raw)
            mol.nep_list     = list(int(x) for x in nep_raw)
            mol.isNode       = list(int(x) for x in is_node_raw)
            if mol.ngs is None:  mol.neighs()
            return list(mol.npi_list), list(mol.nep_list), list(mol.isNode)

        perm_arr = np.array(perm, dtype=np.int32)
        inv_perm = np.empty_like(perm_arr)
        inv_perm[perm_arr] = np.arange(natoms, dtype=np.int32)

        mol.apos = mol.apos[perm_arr].copy()
        if mol.atypes is not None: mol.atypes = mol.atypes[perm_arr].copy()
        if mol.qs     is not None: mol.qs     = mol.qs[perm_arr].copy()
        if mol.Rs     is not None: mol.Rs     = mol.Rs[perm_arr].copy()
        if mol.enames is not None: mol.enames = [mol.enames[i] for i in perm_arr]
        if mol.bonds  is not None:
            remapped_bonds = []
            for a, b in mol.bonds: remapped_bonds.append((int(inv_perm[int(a)]), int(inv_perm[int(b)])))
            mol.bonds = remapped_bonds
        npi_list     = np.asarray(npi_raw, dtype=np.int32)[perm_arr].tolist()
        nep_list     = np.asarray(nep_raw, dtype=np.int32)[perm_arr].tolist()
        is_node      = np.asarray(is_node_raw, dtype=np.int32)[perm_arr].tolist()
        mol.npi_list = npi_list
        mol.nep_list = nep_list
        mol.isNode   = is_node
        mol.perm_nodes_first = perm_arr.tolist()
        mol.perm_inverse = inv_perm.tolist()
        mol.ngs = None
        mol.neighs()
        return npi_list, nep_list, is_node

    def toMMFFsp3_loc(self, mol, atom_types, bRealloc=True, bEPairs=False, bUFF=False, *, lone_pairs_pi=False, align_pi_vectors=False):
        """
        Converts an AtomicSystem to the MMFFsp3_loc representation.

        Parameters:
        - mol (AtomicSystem): The atomic system containing all data.
        - atom_types (dict): Dictionary mapping atom names to AtomType instances.
        - bRealloc (bool): Flag to indicate if reallocation is needed.
        - bEPairs (bool): Flag to include electron pairs.
        - bUFF (bool): Flag to use UFF parameters.
        """
        ang0s = [109.5 * np.pi / 180.0, 120.0 * np.pi / 180.0, 180.0 * np.pi / 180.0]  # in radians

        natom  = len(mol.apos)

        if self.reorder_nodes_first:
            npi_list, nep_list, isNode = self._ensure_node_first(mol, atom_types)
        else:
            npi_list, nep_list, isNode = initAtomProperties(mol, atom_types, self.capping_atoms, bPrint=False)
            if isinstance(npi_list, np.ndarray): npi_list = npi_list.tolist()
            if isinstance(nep_list, np.ndarray): nep_list = nep_list.tolist()
            if isinstance(isNode, np.ndarray): isNode = isNode.tolist()
            mol.npi_list = list(npi_list)
            mol.nep_list = list(nep_list)
            mol.isNode_original = list(int(x) for x in isNode)
            isNode = [1] * natom
            mol.isNode = list(isNode)
        REQs     = getattr(mol, 'REQs', MMparams.generate_REQs_from_atom_types(mol, atom_types))
        #self.REQs[:,:] = REQs[:,:] # self.REQs not yet allocated

        nCmax  = len(mol.ngs)  # Assuming 'ngs' is a list of neighbor lists per atom
        ngs    = mol.ngs
        nngs   = np.zeros(len(ngs), dtype=np.int32)
        
        if isinstance(ngs[0], dict):
            if verbosity > 0:
                print( "!!!!!!!!!!! CONVERT FROM DICT TO ARRAY !!!!!!!!!!!!" )
            ngs_ = np.empty((nCmax,4), dtype=np.int32)
            ngbs = np.empty((nCmax,4), dtype=np.int32)
            ngs_[:,:]=-1
            for ia, ng in enumerate(ngs):
                i = 0
                nngs[ia] = len(ng)
                for k, v in ng.items():
                    ngbs[ia,i] = v # neigh bond
                    ngs_[ia,i] = k # neigh atom
                    i += 1
            ngs = ngs_
        else:
            for ia, row in enumerate(ngs):
                if isinstance(row, dict):
                    nngs[ia] = len(row)
                else:
                    nngs[ia] = sum(1 for nb in row if nb >= 0)

        for ia in range(natom):
            at_name = mol.enames[ia]
            at_type = atom_types.get(at_name)
            if at_type is None:
                continue
            nbond = int(nngs[ia])
            conf_ne = int(getattr(at_type, 'nepair', 0))
            val_total = int(getattr(at_type, 'valence', 4)) + conf_ne
            npi_calc = val_total - nbond - conf_ne
            if npi_calc < 0:
                npi_calc = 0
            if npi_calc > 2:
                npi_calc = 2
            npi_list[ia] = int(npi_calc)
            nep_list[ia] = conf_ne

        if lone_pairs_pi:
            augmented = np.array(npi_list, dtype=np.int32)
            for ia in range(natom):
                if augmented[ia] > 0:
                    continue
                if int(nep_list[ia]) <= 0:
                    continue
                neigh_container = ngs[ia]
                if isinstance(neigh_container, dict):
                    neigh_ids = list(neigh_container.keys())
                else:
                    neigh_ids = [int(nb) for nb in neigh_container if int(nb) >= 0]
                for ja in neigh_ids:
                    if ja < 0 or ja >= natom:
                        continue
                    if augmented[ja] > 0:
                        augmented[ia] = 1
                        break
            npi_list = augmented.tolist()

        
        npi_total = sum(npi_list) 
        ne_total  = sum(nep_list) 
        if not bEPairs:
            ne_total = 0
        
        # count isNode > 0
        if self.reorder_nodes_first:
            nnode = isNode.count(1)
            ncap  = natom - nnode
        else:
            nnode = natom
            ncap  = 0
        nb    = len(mol.bonds)
        
        if verbosity > 0:
            print(f"MMFF::toMMFFsp3_loc() natom {natom} nnode {nnode} ncap {ncap} nb {nb} npi {npi_total} ne {ne_total}")

        ntors = 0
        if bRealloc:
            self.realloc(nnode=nnode, ncap=ncap + ne_total, ntors=ntors)
            self.REQs[:,:] = REQs[:,:]

        # Assign atom types and positions
        etyp = atom_types.get("E", None)
        if etyp is None:
            raise ValueError("atom_types does not contain key 'E'.")


        # Initialize neighbors
        self.neighs[:] = -1  # Set all neighbors to -1
        self.bLs[:] = 0.0
        self.bKs[:] = 0.0
        self.Ksp[:] = 0.0
        self.Kpp[:] = 0.0

        # Assign positions and types
        for ia in range(natom):
            A_pos        = mol.apos[ia]
            A_type_index = mol.atypes[ia]
            A_ename      = mol.enames[ia]
            atom_type    = atom_types.get(A_ename, None)
            if atom_type is None:
                raise ValueError(f"Atom type '{A_ename}' not found in atom_types.")

            self.apos  [ia,:3] = A_pos.astype(np.float32)
            self.atypes[ia]    = A_type_index

            ngi   = ngs[ia]
            nbond = nngs[ia]

            if isNode[ia] <= 0:
                self.neighs[ia,:] =  ngi
            else:
                # Assign parameters
                conf_index = 0  # Default configuration index, no longer using iconf
                # Try to get neighbors for this atom, if available
                conf_npi = max(int(getattr(atom_type, 'npi', 0)), int(npi_list[ia]))
                conf_ne  = int(nep_list[ia])
                if conf_npi <= 0:
                    val_total = int(getattr(atom_type, 'valence', 4)) + conf_ne
                    npi_calc = val_total - nbond - conf_ne
                    if npi_calc < 0:
                        npi_calc = 0
                    if npi_calc > 2:
                        npi_calc = 2
                    conf_npi = npi_calc

                if conf_npi > 2:
                    print(f"ERROR in MM::Builder::toMMFFsp3_loc(): atom[{ia}].conf.npi({conf_npi}) > 2 => exit()")
                    self.printAtomConf(ia, mol)
                    exit(0)

                # Angle parameters
                ang0 = atom_type.Ass * deg2rad
                ang0 *= 0.5
                self.apars[ia, 0] = np.cos(ang0)    # ssC0
                self.apars[ia, 1] = np.sin(ang0)    # ssC1
                self.apars[ia, 2] = atom_type.Kss * 4.0  # ssK
                self.apars[ia, 3] = np.sin(atom_type.Ass * deg2rad)  # piC0


                if verbosity > 0: print( "nbond", nbond )
                # Setup neighbors - ngs is already populated above
                hs = np.zeros((4, 3), dtype=np.float32)
                print(f"DEBUG MMFF: ia={ia} {A_ename} nbond={nbond} conf_npi={conf_npi} conf_ne={conf_ne}")
                for k in range(nbond):
                    ja = ngi[k]  # ja is the atom index of the neighbor
                    if verbosity > 2: print( "ia,ja", ia,ja )
                    if ja < 0 or ja >= natom:
                        continue
                        
                    # Find the bond between atoms ia and ja
                    bond_index = ngbs[ia][k]

                    bond = mol.bonds[bond_index]
                    Aj = mol.apos[ja]
                    jtyp_ename = mol.enames[ja]
                    jtyp = atom_types.get(jtyp_ename, None)
                    if jtyp is None:
                        raise ValueError(f"Atom type '{jtyp_ename}' not found in atom_types.")

                    hs_k = Aj - A_pos
                    norm = np.linalg.norm(hs_k)
                    if norm != 0:
                        hs_k /= norm
                    else:
                        hs_k = np.array([0.0, 0.0, 0.0], dtype=np.float32)
                    hs[k] = hs_k.astype(np.float32)

                    # Assign neighbor
                    self.neighs[ia, k] = ja

                    if bUFF:
                        # Assign bond parameters using UFF
                        rij, kij = self.assignBondParamsUFF(ibond_index, ia, ja, mol, atom_types)
                        self.bLs[ia, k] = rij
                        self.bKs[ia, k] = kij
                    else:
                        rij, kij = self.assignBondParamsSimple(ia, ja, mol, atom_types)
                        self.bLs[ia, k] = rij
                        self.bKs[ia, k] = kij

                    # Assign Ksp
                    if conf_npi > 0 or conf_ne > 0:
                        self.Ksp[ia, k] = atom_type.Ksp
                    else:
                        self.Ksp[ia, k] = 0.0

                    # Assign Kpp
                    #nej  = nep_list[ja]
                    #npij = npi_list[ja]
                    if atom_type.Kpp > 0 and jtyp.Kpp > 0:
                        self.Kpp[ia, k] = np.sqrt(atom_type.Kpp * jtyp.Kpp)
                    else:
                        self.Kpp[ia, k] = 0.0

                # Setup configuration geometry based on bonds
                if nbond <= 1:
                    self.pipos[ia] = np.array([0.0, 0.0, 0.0], dtype=np.float32)
                
                hs_filled = self.makeConfGeom(nb=nbond, npi=conf_npi)
                if nbond < 4:
                    for idx in range(nbond, 4):
                        hs_filled[idx] = hs[idx]
                if conf_npi > 0 and hasattr(mol, 'get_atomi_pi_direction'):
                    pi_vec = mol.get_atomi_pi_direction(ia)
                else:
                    pi_vec = hs_filled[3]
                norm_pi = np.linalg.norm(pi_vec)
                if norm_pi >= 1e-6:
                    pi_vec = (pi_vec / norm_pi).astype(np.float32)
                else:
                    pi_vec = np.array([0.0, 0.0, 0.0], dtype=np.float32)
                self.pipos[ia] = pi_vec
                print(f"DEBUG MMFF: init pi_vec ia={ia} -> {self.pipos[ia]} (raw={hs_filled[3]})")

                if conf_npi > 0:
                    neigh_ids = [int(x) for x in self.neighs[ia] if int(x) >= 0]
                    hs_raw = hs_filled[3]
                    self._dbg_topo(
                        f"PI_INIT ia={ia:4d} nb={int(nbond)} npi={int(conf_npi)} neigh={neigh_ids} "
                        f"raw=({hs_raw[0]: .6f},{hs_raw[1]: .6f},{hs_raw[2]: .6f}) "
                        f"pi=({self.pipos[ia][0]: .6f},{self.pipos[ia][1]: .6f},{self.pipos[ia][2]: .6f})"
                    )

                if bEPairs:
                    # Generate electron pairs
                    for k in range(nbond, nbond + conf_ne):
                        if k >= 4:
                            print(f"ERROR: Atom {ia} has more than 4 neighbors; additional electron pairs are ignored.")
                            exit(0)
                        ie = nnode + ncap + k - nbond
                        if ie >= self.nvecs:
                            print(f"ERROR: Electron pair index {ie} exceeds nvecs {self.nvecs}.")
                            exit(0)
                        self.neighs[ia, k] = ie
                        self.apos[ie] = self.apos[ia] + hs_filled[k - nbond] * Lepair
                        self.atypes[ie] = atom_types["E"].name  # Assuming 'E' is the type index for electron pairs
                        self.bKs[ia, k] = Kepair
                        self.bLs[ia, k] = Lepair
                        if conf_npi > 0:
                            self.Ksp[ia, k] = atom_type.Ksp
                        else:
                            self.Ksp[ia, k] = 0.0

        mol.npi_list = [int(x) for x in npi_list]
        mol.nep_list = [int(x) for x in nep_list]

        self._propagate_pi_dirs(mol, npi_list)
        if align_pi_vectors:
            self._align_pi_signs(mol, npi_list)

        # Normalize pi-orbital vectors and provide a safe fallback if zero
        if self.pipos is not None:
            for i in range(self.nnode):
                v = self.pipos[i]
                n = np.linalg.norm(v)
                if not np.isfinite(n) or n < 1e-8:
                    print(f"WARNING MMFF: pi_dir unresolved for node {i}, leaving zero vector")
                    self.pipos[i] = np.array([0.0, 0.0, 0.0], dtype=np.float32)
                else:
                    self.pipos[i] = (v / n).astype(np.float32)
                if npi_list[i] > 0:
                    print(f"DEBUG MMFF: final pi_vec ia={i} -> {self.pipos[i]} (norm={n})")
                    self._dbg_topo(
                        f"PI_FINAL ia={i:4d} npi={int(npi_list[i])} "
                        f"pi=({self.pipos[i][0]: .6f},{self.pipos[i][1]: .6f},{self.pipos[i][2]: .6f}) norm={float(n): .6f}"
                    )

        # Write pi-orbital unit vectors into the tail of apos (indices natoms : natoms+nnode)
        # This matches kernels expecting apos[iav + nAtoms] to hold pi orientation for each node atom
        if self.pipos is not None and self.pipos.shape[0] >= self.nnode:
            self.apos[self.natoms:self.natoms + self.nnode, :3] = self.pipos[:self.nnode, :]
            self.apos[self.natoms:self.natoms + self.nnode,  3] = 0.0

        # Prepare simple back-neighbor placeholder (atom indices). GPU expects indices into fneigh; left for future mapping.
        self.make_back_neighs()

    def _propagate_pi_dirs(self, mol, npi_list, *, min_norm=0.7, max_iter=4):
        if self.pipos is None:
            return
        ngs = getattr(mol, 'ngs', None)
        if ngs is None:
            return
        nnode = min(self.nnode, len(npi_list))
        for _ in range(max_iter):
            updated = False
            for ia in range(nnode):
                if npi_list[ia] <= 0:
                    continue
                v_host = self.pipos[ia]
                host_norm = np.linalg.norm(v_host)
                if host_norm >= min_norm:
                    continue
                acc = np.zeros(3, dtype=np.float32)
                neigh_container = ngs[ia] if ia < len(ngs) else []
                if isinstance(neigh_container, dict):
                    neigh_ids = list(neigh_container.keys())
                else:
                    neigh_ids = [int(nb) for nb in neigh_container if int(nb) >= 0]
                for ja in neigh_ids:
                    if ja < 0 or ja >= nnode:
                        continue
                    vj = self.pipos[ja]
                    norm_j = np.linalg.norm(vj)
                    if norm_j < 1e-6:
                        continue
                    vj_unit = (vj / norm_j).astype(np.float32)
                    if host_norm > 1e-6 and np.dot(v_host, vj_unit) < 0.0:
                        vj_unit = -vj_unit
                    elif np.linalg.norm(acc) > 1e-6 and np.dot(acc, vj_unit) < 0.0:
                        vj_unit = -vj_unit
                    acc += vj_unit
                acc_norm = np.linalg.norm(acc)
                if acc_norm >= min_norm:
                    self.pipos[ia] = (acc / acc_norm).astype(np.float32)
                    print(f"DEBUG MMFF: propagate pi_vec ia={ia} -> {self.pipos[ia]} using neighbors {neigh_ids}")
                    self._dbg_topo(
                        f"PI_PROP ia={ia:4d} neigh={neigh_ids} "
                        f"acc=({acc[0]: .6f},{acc[1]: .6f},{acc[2]: .6f}) acc_norm={float(acc_norm): .6f} "
                        f"pi=({self.pipos[ia][0]: .6f},{self.pipos[ia][1]: .6f},{self.pipos[ia][2]: .6f})"
                    )
                    updated = True
            if not updated:
                break


    def _align_pi_signs(self, mol, npi_list, *, max_iter=6):
        if self.pipos is None:
            return
        ngs = getattr(mol, 'ngs', None)
        if ngs is None:
            return
        nnode = min(self.nnode, len(npi_list))
        for _ in range(max_iter):
            flipped = False
            for ia in range(nnode):
                if npi_list[ia] <= 0:
                    continue
                v_host = self.pipos[ia]
                host_norm = np.linalg.norm(v_host)
                if host_norm < 1e-6:
                    continue
                neigh_container = ngs[ia] if ia < len(ngs) else []
                if isinstance(neigh_container, dict):
                    neigh_ids = list(neigh_container.keys())
                else:
                    neigh_ids = [int(nb) for nb in neigh_container if int(nb) >= 0]
                baseline = None
                for ja in neigh_ids:
                    if ja < 0 or ja >= nnode or npi_list[ja] <= 0:
                        continue
                    vj = self.pipos[ja]
                    norm_j = np.linalg.norm(vj)
                    if norm_j < 1e-6:
                        continue
                    baseline = (vj / norm_j).astype(np.float32)
                    if np.dot(v_host, baseline) < 0.0:
                        self.pipos[ia] = -v_host
                        flipped = True
                        self._dbg_topo(
                            f"PI_FLIP ia={ia:4d} ja={int(ja):4d} "
                            f"dot={float(np.dot(v_host, baseline)): .6f} "
                            f"pi=({self.pipos[ia][0]: .6f},{self.pipos[ia][1]: .6f},{self.pipos[ia][2]: .6f})"
                        )
                    break
            if not flipped:
                break


    def assignBondParamsSimple(self, ia, ja, mol, atom_types ):
        ti = atom_types[mol.enames[ia]]
        tj = atom_types[mol.enames[ja]]
        rij = ti.Ruff + tj.Ruff
        kij = np.sqrt( tj.Kss * ti.Kss )
        return rij, kij
        

    def assignBondParamsUFF(self, ib, ai, aj, mol, atom_types):
        """
        Calculate bond parameters for a bond in the molecule using UFF.

        Args:
        - ib (int): Bond index.
        - ai (int): Index of the first atom in the bond.
        - aj (int): Index of the second atom in the bond.
        - mol (AtomicSystem): Molecular system.
        - atom_types (dict): Dictionary of atom types.

        Returns:
        - (float, float): Tuple of bond length (rij) and bond stiffness (kij).
        """
        bond = mol.bonds[ib]
        
        # Load element types if not already loaded
        if not hasattr(self, 'element_types'):
            base_path = os.path.dirname(os.path.abspath(__file__))
            data_path = os.path.join(base_path, "../../cpp/common_resources/")
            self.element_types = MMparams.read_element_types(os.path.join(data_path, 'ElementTypes.dat'))
        
        # Use getattr to get isNode or create it if not exist
        capping_atoms = ['H', 'F', 'Cl', 'Br', 'I']
        isNode = getattr(mol, 'isNode', 
                         np.array([0 if atom_types[e].element_name in capping_atoms else 1 
                                  for e in mol.enames]))
        
        # Use getattr to get npi_list or create it if not exist
        npi_list = getattr(mol, 'npi_list', 
                          np.array([atom_types[e].npi for e in mol.enames]))
        
        # Get the number of pi electrons for each atom in the bond
        if isNode[ai]:
            npi = npi_list[ai]
        else:
            npi = 0
        if isNode[aj]:
            npj = npi_list[aj]
        else:
            npj = 0

        # Get atom types for both atoms
        ti = atom_types[mol.enames[ai]]
        tj = atom_types[mol.enames[aj]]
        
        # Get element properties
        element_i = self.element_types[ti.element_name]
        element_j = self.element_types[tj.element_name]
        
        # Get electronegativities from ElementType
        Ei = element_i.Eaff  # Electronegativity
        Ej = element_j.Eaff  # Electronegativity
        
        # Use Qbase from AtomType for charges
        Qi = ti.Qbase
        Qj = tj.Qbase

        # Bond order calculation
        BO = 1 + min(npi, npj)
        
        # UFF natural bond radius from AtomType
        ri = ti.Ruff
        rj = tj.Ruff

        # Compute rBO (bond order correction) and rEN (electronegativity correction)
        if BO > 0:
            rBO = -0.1332 * (ri + rj) * np.log(BO)
        else:
            rBO = 0.0
            
        denominator = Ei * ri + Ej * rj
        if denominator != 0.0 and Ei >= 0 and Ej >= 0:
            rEN = ri * rj * ((np.sqrt(Ei) - np.sqrt(Ej))**2) / denominator
        else:
            # Default to zero if we can't calculate properly
            rEN = 0.0

        # Calculate the natural bond length
        rij = ri + rj + rBO - rEN

        # Calculate force constant in kcal/(mol*Angstrom^2)
        Zi = Qi  # Formal charge
        Zj = Qj  # Formal charge
        kij = 664.12 * ((Zi * Zj) / (rij**3)) * ri * rj

        # Debug information (commented out)
        # print(f"bondUFF[{ti.name},{tj.name},{BO}] kij={kij} rij={rij}({ri},{rj}|{rBO},{rEN}) k={kij}({Qi},{Qj}) E=({Ei},{Ej}) {ti.name} {tj.name}")

        return (rij, kij)

    def makeConfGeom(self, nb, npi):
        """
        Sets up the configuration geometry based on bond vectors.

        Parameters:
        - nb (int): Number of bonds.
        - npi (int): Number of pi-orbitals.
        - hs (np.ndarray): Preallocated array with shape (4, 3) to store bond direction vectors.

        Returns:
        - np.ndarray: The filled hs array with shape (4, 3).
        """
        hs = np.zeros((4, 3), dtype=np.float32)
        if nb == 3:
            # Defined by 3 sigma bonds
            cross_prod = np.cross(hs[1] - hs[0], hs[2] - hs[0])
            norm = np.linalg.norm(cross_prod)
            if norm != 0:
                cross_prod /= norm
            else:
                cross_prod = np.array([0.0, 0.0, 0.0], dtype=np.float32)

            if npi == 0:
                # sp3 no-pi
                dot_product = np.dot(cross_prod, hs[0] + hs[1] + hs[2])
                if dot_product > 0:
                    cross_prod *= -1
                hs[3] = cross_prod
            else:
                hs[3] = cross_prod

        elif nb == 2:
            # Defined by 2 sigma bonds
            cross_prod = np.cross(hs[0], hs[1])
            norm = np.linalg.norm(cross_prod)
            if norm != 0:
                cross_prod /= norm
            else:
                cross_prod = np.array([0.0, 0.0, 0.0], dtype=np.float32)

            if npi == 0:
                # -CH2- like sp3 no-pi => 109.5 degrees
                cb = 0.81649658092  # sqrt(2/3)
                cc = 0.57735026919  # sqrt(1/3)
                hs[2] = cross_prod * cc + hs[0] * cb
                hs[3] = cross_prod * cc - hs[0] * cb
            elif npi == 1:
                # =CH- like sp1-pi
                hs[2] = cross_prod
                hs[3] = hs[0]
            else:
                # #C- like sp2-pi
                hs[2] = -cross_prod
                hs[3] = hs[0]

        elif nb == 1:
            # Defined by 1 sigma bond
            norm_c = np.linalg.norm(hs[0])
            if norm_c != 0:
                c = hs[0] / norm_c
            else:
                c = np.array([0.0, 0.0, 0.0], dtype=np.float32)

            # Find an orthogonal vector to c
            if abs(c[0]) < abs(c[1]):
                if abs(c[0]) < abs(c[2]):
                    ortho1 = np.array([1.0, 0.0, 0.0], dtype=np.float32)
                else:
                    ortho1 = np.array([0.0, 0.0, 1.0], dtype=np.float32)
            else:
                if abs(c[1]) < abs(c[2]):
                    ortho1 = np.array([0.0, 1.0, 0.0], dtype=np.float32)
                else:
                    ortho1 = np.array([0.0, 0.0, 1.0], dtype=np.float32)

            a = np.cross(c, ortho1)
            norm_a = np.linalg.norm(a)
            if norm_a != 0:
                a /= norm_a
            else:
                a = np.array([0.0, 0.0, 0.0], dtype=np.float32)

            b = np.cross(c, a)
            norm_b = np.linalg.norm(b)
            if norm_b != 0:
                b /= norm_b
            else:
                b = np.array([0.0, 0.0, 0.0], dtype=np.float32)

            if npi == 0:
                # -CH3 like sp3 no-pi => 109.5 degrees
                ca = 0.81649658092  # sqrt(2/3)
                cb = 0.47140452079  # sqrt(2/9)
                cc = -0.33333333333  # 1/3
                hs[1] = c * cc + b * (cb * 2)
                hs[2] = c * cc - b * cb + a * ca
                hs[3] = c * cc - b * cb - a * ca
            elif npi == 1:
                # =CH2 like sp2 1-pi => 60 degrees
                ca = 0.86602540378  # cos(30 degrees)
                cc = -0.5           # sqrt(1/8)
                hs[1] = c * cc + a * ca
                hs[2] = c * cc - a * ca
                hs[3] = b
            else:
                # #CH like sp2-pi
                hs[1] = c * -1
                hs[2] = b
                hs[3] = a

        elif nb == 0:
            # Defined by 0 sigma bonds
            # Assume hs[0] is already defined
            norm_c = np.linalg.norm(hs[0])
            if norm_c != 0:
                c = hs[0] / norm_c
            else:
                c = np.array([0.0, 0.0, 0.0], dtype=np.float32)

            # Find an orthogonal vector to c
            if abs(c[0]) < abs(c[1]):
                if abs(c[0]) < abs(c[2]):
                    ortho1 = np.array([1.0, 0.0, 0.0], dtype=np.float32)
                else:
                    ortho1 = np.array([0.0, 0.0, 1.0], dtype=np.float32)
            else:
                if abs(c[1]) < abs(c[2]):
                    ortho1 = np.array([0.0, 1.0, 0.0], dtype=np.float32)
                else:
                    ortho1 = np.array([0.0, 0.0, 1.0], dtype=np.float32)

            a = np.cross(c, ortho1)
            norm_a = np.linalg.norm(a)
            if norm_a != 0:
                a /= norm_a
            else:
                a = np.array([0.0, 0.0, 0.0], dtype=np.float32)

            b = np.cross(c, a)
            norm_b = np.linalg.norm(b)
            if norm_b != 0:
                b /= norm_b
            else:
                b = np.array([0.0, 0.0, 0.0], dtype=np.float32)

            if npi == 0:
                # CH4 like sp3 no-pi => 109.5 degrees
                ca = 0.81649658092  # sqrt(2/3)
                cb = 0.47140452079  # sqrt(2/9)
                cc = -0.33333333333  # 1/3
                hs[1] = c * cc + b * (cb * 2)
                hs[2] = c * cc - b * cb + a * ca
                hs[3] = c * cc - b * cb - a * ca
                hs[0] = c  # Assuming hs[0] is the original bond direction
        return hs

        # After parameters and orientations are prepared, populate the pi-orbital vectors in the apos buffer


    def printAtomConf(self, ia, mol):
        """
        Prints the configuration of a specific atom for debugging.

        Parameters:
        - ia (int): Index of the atom.
        - mol (AtomicSystem): The atomic system containing all data.
        """
        print(f"printAtomConf({ia}) neighs: {self.neighs[ia]}", end="")
        if ia < self.nnode:
            print(f" apars: {self.apars[ia]} bLs: {self.bLs[ia]} bKs: {self.bKs[ia]} Ksp: {self.Ksp[ia]} Kpp: {self.Kpp[ia]} pipos: {self.pipos[ia]}")
        else:
            print()

    def printArrays(self):
        print("MMFF::printArrays(): natoms: {self.natoms}  nvecs: {self.nvecs}  nnode: {self.nnode}  ncap: {self.ncap}  ntors: {self.ntors}")
        print("apos",   self.apos )
        print("REQs",   self.REQs )
        print("neighs", self.neighs )
        print("bLs",    self.bLs )
        print("bKs",    self.bKs )
        print("apars",  self.apars )
        print("Ksp",    self.Ksp )
        print("Kpp",    self.Kpp )

# IF MAIN
if __name__ == "__main__":

    from ..atomicUtils import AtomicSystem

    # Define AtomType instances with npi and ne
    atom_types = {
        "C": AtomType(name="C", Kss=300.0, Asp=109.5, Ksp=100.0, Kpp=150.0, npi=0, ne=0),
        "H": AtomType(name="H", Kss=200.0, Asp=109.5, Ksp=50.0,  Kpp=75.0,  npi=0, ne=0),
        "O": AtomType(name="O", Kss=350.0, Asp=109.5, Ksp=120.0, Kpp=180.0, npi=0, ne=0),
        "E": AtomType(name="E", Kss=0.0,   Asp=0.0,   Ksp=0.0,   Kpp=0.0,   npi=0, ne=1)  # Electron pairs
    }

    mol = AtomicSystem(
        apos=np.array([[0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, 1.0]], dtype=np.float32),
        atypes=np.array([0, 1, 1, 1], dtype=np.int32),  # Example type indices
        enames=["C", "H", "H", "H"],
        lvec=np.identity(3, dtype=np.float32),
        qs=np.array([0.0, 0.0, 0.0, 0.0], dtype=np.float32),
        Rs=np.array([[1.5, 0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0, 0.0]], dtype=np.float32),
        bonds=[
            Bond(i=0, j=1, l0=1.0, k=300.0),
            Bond(i=0, j=2, l0=1.0, k=300.0),
            Bond(i=0, j=3, l0=1.0, k=300.0)
        ],
        ngs=[[1, 2, 3, -1]],  # Neighbor lists for each atom
        dihedrals=[],          # No dihedrals in this simple example
        iconf=[0, -1, -1, -1], # Atom 0 is in configuration 0, others are not
        npi_list=[0, 0, 0, 0],
        nep_list=[0, 0, 0, 0]
    )

    # Initialize MMFF instance
    mmff = MMFF(bTorsion=True, verbosity=1)

    # Convert AtomicSystem to MMFFsp3_loc representation
    mmff.toMMFFsp3_loc(
        mol=mol,
        atom_types=atom_types,
        bRealloc=True,
        bEPairs=True,
        bUFF=False
    )

    # Optionally, print atom configurations for debugging
    mmff.printAtomConf(0, mol)  # Replace 0 with desired atom index
