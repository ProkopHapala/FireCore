import numpy as np
import math

# Constants (Define these based on your simulation requirements)
Lepair = 1.5    # Example bond length for electron pairs
Kepair = 100.0  # Example bond stiffness for electron pairs
deg2rad = np.pi / 180.0

class AtomType:
    """
    Represents the parameters for each atom type.
    """
    def __init__(self, name, Ruff, Quff, Eaff, Kss, Ass, Ksp, Kpp, npi=0, ne=0):
        """
        Initializes the AtomType.

        Parameters:
        - name (str): Name of the atom type (e.g., 'C', 'H').
        - Kss (float): Stiffness for sigma-sigma interactions.
        - Ass (float): Angle parameter in degrees.
        - Ksp (float): Stiffness for pi-sigma interactions.
        - Kpp (float): Stiffness for pi-pi interactions.
        - npi (int): Number of pi-orbitals.
        - ne (int): Number of electron pairs.
        """
        self.name = name
        self.Kss = Kss
        self.Ass = Ass
        self.Ksp = Ksp
        self.Kpp = Kpp
        self.npi = npi
        self.ne = ne
        self.Ruff = Ruff
        self.Quff = Quff
        self.Eaff = Eaff

class Bond:
    """
    Represents a bond between two atoms.
    """
    def __init__(self, i, j, l0, k):
        """
        Initializes the Bond.

        Parameters:
        - i (int): Index of the first atom.
        - j (int): Index of the second atom.
        - l0 (float): Equilibrium bond length.
        - k (float): Bond stiffness.
        """
        self.i = i
        self.j = j
        self.l0 = l0
        self.k = k

    def get_neighbor_atom(self, atom_index):
        """
        Given one atom index, returns the index of the other atom in the bond.

        Parameters:
        - atom_index (int): Index of one atom in the bond.

        Returns:
        - int: Index of the other atom.
        """
        return self.j if atom_index == self.i else self.i

class Dihedral:
    """
    Represents a dihedral (torsion) between four atoms.
    """
    def __init__(self, atoms, a0, k, n):
        """
        Initializes the Dihedral.

        Parameters:
        - atoms (tuple of int): Tuple of four atom indices defining the dihedral.
        - a0 (float): Equilibrium dihedral angle in radians.
        - k (float): Dihedral stiffness.
        - n (int): Periodicity of the dihedral.
        """
        self.atoms = atoms  # Tuple of four atom indices
        self.a0 = a0        # Equilibrium angle in radians
        self.k = k          # Stiffness
        self.n = n          # Periodicity

class MMFF:
    """
    Represents the MMFF (Merck Molecular Force Field) parameters and configurations.
    """
    def __init__(self, bTorsion=False, verbosity=1):
        """
        Initializes the MMFF instance.

        Parameters:
        - bTorsion (bool): Flag to indicate if torsions are considered.
        - verbosity (int): Level of verbosity for logging.
        """
        self.bTorsion = bTorsion
        self.verbosity = verbosity
        self.natoms = 0
        self.nnode  = 0
        self.ncap   = 0
        self.nvecs  = 0
        self.ntors  = 0
        self.nDOFs  = 0

        # Initialize arrays
        self.apos = None          # [nvecs,  3 ] Positions
        self.fapos = None         # [nvecs,  3 ] Forces
        self.atypes = None        # [natoms    ]  Atom types 
        self.neighs = None        # [natoms, 4 ]  Neighbor indices 
        self.neighCell = None     # [natoms, 4 ] Neighbor cell indices
        self.REQs = None          # [natoms, 4 ] Non-covalent parameters
        self.apars = None         # [nnode,  4 ] Angle parameters
        self.bLs = None           # [nnode,  4 ] Bond lengths
        self.bKs = None           # [nnode,  4 ] Bond stiffness
        self.Ksp = None           # [nnode,  4 ] Pi-sigma stiffness
        self.Kpp = None           # [nnode,  4 ] Pi-pi stiffness
        self.angles = None        # [nnode, 6, 3 ] Angles between bonds
        # self.tors2atom = None     # [ntors]        Torsion atom indices
        # self.torsParams = None    # [ntors, 4    ] Torsion parameters
        # self.constr = None        # [natoms, 4   ] Constraints
        # self.constrK = None       # [natoms, 4   ] Constraint stiffness
        self.invLvec = None       # [nSystems, 9 ] Inverse lattice vectors
        self.pbc_shifts = None    # [nSystems, 3 ] PBC shifts
        self.nPBC = None
        self.npbc = None

        # Pi-orbitals and electron pairs
        self.pipos = None         # [natoms, 3] Pi-orbital orientations

    def realloc(self, nnode, ncap, ntors=0, nPBC=(0,0,0) ):
        """
        Reallocates memory for MMFF parameters based on the system size.

        Parameters:
        - nnode (int): Number of nodes (configurations).
        - ncap (int): Number of capping atoms.
        - ntors (int): Number of torsions.
        """
        self.npbc = (nPBC[0]*2+1)*(nPBC[1]*2+1)*(nPBC[2]*2+1)
        self.nPBC = nPBC
        self.nnode  = nnode
        self.ncap   = ncap
        self.natoms = nnode + ncap
        self.nvecs  = self.natoms + nnode  # Including pi-orbitals
        self.ntors  = ntors
        self.nDOFs  = self.nvecs * 3

        # Initialize arrays with default values
        self.apos       = np.full((self.nvecs, 3), 0.0, dtype=np.float32)
        self.fapos      = np.full((self.nvecs, 3), 0.0, dtype=np.float32)
        self.atypes     = np.full( self.natoms,    -1,  dtype=np.int32)
        self.neighs     = np.full((self.natoms,4), -1,  dtype=np.int32)
        self.neighCell  = np.full((self.natoms,4), -1,  dtype=np.int32)
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

    def countPiE(self, atomic_system):
        """
        Counts the total number of pi-orbitals and electron pairs.

        Parameters:
        - atomic_system (AtomicSystem): The atomic system.

        Returns:
        - (int, int): Tuple of total pi-orbitals and electron pairs.
        """
        npi = sum(atomic_system.npi_list) if hasattr(atomic_system, 'npi_list') else 0
        ne = sum(atomic_system.nep_list) if hasattr(atomic_system, 'nep_list') else 0
        return npi, ne

    def make_back_neighs(self):
        """
        Creates back neighbors based on neighbor lists.
        """
        self.back_neighs = np.full_like(self.neighs, -1)
        for ia in range(self.natoms):
            for k in range(4):
                ib = self.neighs[ia, k]
                if ib >= 0 and ib < self.natoms:
                    for kk in range(4):
                        if self.neighs[ib, kk] == ia:
                            self.back_neighs[ia, k] = ib
                            break

    def toMMFFsp3_loc(self, atomic_system, AtomTypeDict, bRealloc=True, bEPairs=True, bUFF=False):
        """
        Converts an AtomicSystem to the MMFFsp3_loc representation.

        Parameters:
        - atomic_system (AtomicSystem): The atomic system containing all data.
        - AtomTypeDict (dict): Dictionary mapping atom names to AtomType instances.
        - bRealloc (bool): Flag to indicate if reallocation is needed.
        - bEPairs (bool): Flag to include electron pairs.
        - bUFF (bool): Flag to use UFF parameters.
        """
        ang0s = [109.5 * np.pi / 180.0, 120.0 * np.pi / 180.0, 180.0 * np.pi / 180.0]  # in radians

        natom = len(atomic_system.apos)
        nCmax = len(atomic_system.ngs)  # Assuming 'ngs' is a list of neighbor lists per atom

        ngs  = atomic_system.ngs
        nngs = np.zeros(len(ngs), dtype=np.int32)
        isNode = atomic_system.isNode
        
        if isinstance(ngs[0], dict):
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

        
        npi_total, ne_total = self.countPiE(atomic_system)
        if not bEPairs:
            ne_total = 0
        
        # count isNode > 0
        nnode = len(atomic_system.isNode[atomic_system.isNode > 0])
        ncap = natom - nnode
        nb = len(atomic_system.bonds)
        
        if self.verbosity > 0:
            print(f"MM::Builder::toMMFFsp3_loc() nnode {nnode} ncap {ncap} npi {npi_total} ne {ne_total}")

        ntors = 0
        if bRealloc:
            self.realloc(nnode=nnode, ncap=ncap + ne_total, ntors=ntors)

        # Assign atom types and positions
        etyp = AtomTypeDict.get("E", None)
        if etyp is None:
            raise ValueError("AtomTypeDict does not contain key 'E'.")

        # Initialize neighbors
        self.neighs[:] = -1  # Set all neighbors to -1
        self.bLs[:] = 0.0
        self.bKs[:] = 0.0
        self.Ksp[:] = 0.0
        self.Kpp[:] = 0.0

        

        # Assign positions and types
        for ia in range(natom):
            A_pos = atomic_system.apos[ia]
            A_type_index = atomic_system.atypes[ia]
            A_ename = atomic_system.enames[ia]
            atom_type = AtomTypeDict.get(A_ename, None)
            if atom_type is None:
                raise ValueError(f"Atom type '{A_ename}' not found in AtomTypeDict.")

            self.apos[ia] = A_pos.astype(np.float32)
            self.atypes[ia] = A_type_index

            ngi   = ngs[ia]
            nbond = nngs[ia]

            if isNode[ia] <= 0:
                self.neighs[ia,:] =  ngi
            else:
                # Assign parameters
                conf_index = 0  # Default configuration index, no longer using iconf
                # Try to get neighbors for this atom, if available
                conf_npi = atom_type.npi
                conf_ne  = atom_type.ne

                if conf_npi > 2:
                    print(f"ERROR in MM::Builder::toMMFFsp3_loc(): atom[{ia}].conf.npi({conf_npi}) > 2 => exit()")
                    self.printAtomConf(ia, atomic_system)
                    exit(0)

                # Angle parameters
                ang0 = atom_type.Ass * deg2rad
                ang0 *= 0.5
                self.apars[ia, 0] = np.cos(ang0)    # ssC0
                self.apars[ia, 1] = np.sin(ang0)    # ssC1
                self.apars[ia, 2] = atom_type.Kss * 4.0  # ssK
                self.apars[ia, 3] = np.sin(atom_type.Ass * deg2rad)  # piC0


                if self.verbosity > 0: print( "nbond", nbond )
                # Setup neighbors - ngs is already populated above
                hs = np.zeros((4, 3), dtype=np.float32)
                for k in range(nbond):
                    ja = ngi[k]  # ja is the atom index of the neighbor
                    print( "ia,ja", ia,ja )
                    if ja < 0 or ja >= natom:
                        continue
                        
                    # Find the bond between atoms ia and ja
                    bond_index = ngbs[ia][k]

                        
                    bond = atomic_system.bonds[bond_index]
                    Aj = atomic_system.apos[ja]
                    jtyp_ename = atomic_system.enames[ja]
                    jtyp = AtomTypeDict.get(jtyp_ename, None)
                    if jtyp is None:
                        raise ValueError(f"Atom type '{jtyp_ename}' not found in AtomTypeDict.")

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
                        rij, kij = self.assignBondParamsUFF(ibond_index, ia, ja, atomic_system, AtomTypeDict)
                        self.bLs[ia, k] = rij
                        self.bKs[ia, k] = kij
                    else:
                        rij, kij = self.assignBondParamsSimple(ia, ja, atomic_system, AtomTypeDict)
                        self.bLs[ia, k] = rij
                        self.bKs[ia, k] = kij

                    # Assign Ksp
                    if conf_npi > 0 or conf_ne > 0:
                        self.Ksp[ia, k] = atom_type.Ksp
                    else:
                        self.Ksp[ia, k] = 0.0

                    # Assign Kpp
                    nej  = atomic_system.nep_list[ja]
                    npij = atomic_system.npi_list[ja]
                    if atom_type.Kpp > 0 and jtyp.Kpp > 0:
                        self.Kpp[ia, k] = np.sqrt(atom_type.Kpp * jtyp.Kpp)
                    else:
                        self.Kpp[ia, k] = 0.0

                # Setup configuration geometry based on bonds
                
                hs_filled = self.makeConfGeom(nb=nbond, npi=conf_npi)
                # hs_filled is a (4,3) array; assign pi-orbital orientations
                if nbond < 4:
                    # Fill remaining hs with existing directions
                    for idx in range(nbond, 4):
                        hs_filled[idx] = hs[idx]  # Existing bond directions or zeros
                self.pipos[ia] = hs_filled[3]  # Assign the fourth direction as pi orientation

                if bEPairs:
                    # Generate electron pairs
                    for k in range(nbond, nbond + conf_ne):
                        if k >= 4:
                            print(f"WARNING: Atom {ia} has more than 4 neighbors; additional electron pairs are ignored.")
                            break
                        ie = nnode + ncap + k - nbond
                        if ie >= self.nvecs:
                            print(f"ERROR: Electron pair index {ie} exceeds nvecs {self.nvecs}.")
                            continue
                        self.neighs[ia, k] = ie
                        self.apos[ie] = self.apos[ia] + hs_filled[k - nbond] * Lepair
                        self.atypes[ie] = AtomTypeDict["E"].name  # Assuming 'E' is the type index for electron pairs
                        self.bKs[ia, k] = Kepair
                        self.bLs[ia, k] = Lepair
                        if conf_npi > 0:
                            self.Ksp[ia, k] = atom_type.Ksp
                        else:
                            self.Ksp[ia, k] = 0.0
   
    def assignBondParamsSimple(self, ia, ja, atomic_system, AtomTypeDict ):
        ti = AtomTypeDict[atomic_system.enames[ia]]
        tj = AtomTypeDict[atomic_system.enames[ja]]
        rij = ti.Ruff + tj.Ruff
        kij = np.sqrt( tj.Kss * ti.Kss )
        return rij, kij
        

    def assignBondParamsUFF(self, ib, ai, aj, atomic_system, AtomTypeDict):
        """
        Assigns bond parameters using UFF (Universal Force Field).

        Parameters:
        - ib (int): Index of the bond in the atomic system's bond list.
        - atomic_system (AtomicSystem): The atomic system containing all data.
        - AtomTypeDict (dict): Dictionary mapping atom names to AtomType instances.

        Returns:
        - (float, float): Tuple of bond length (rij) and bond stiffness (kij).
        """
        bond = atomic_system.bonds[ib]
        isNode = atomic_system.isNode
        # ai = bond.i
        # aj = bond.j
        if isNode[ai]:
            npi = atomic_system.npi_list[ai]
        else:
            npi = 0
        if isNode[aj]:
            npj = atomic_system.npi_list[aj]
        else:
            npj = 0

        ti = AtomTypeDict[atomic_system.enames[ai]]
        tj = AtomTypeDict[atomic_system.enames[aj]]
        # Assuming ElementType is not needed; use Eaff directly
        Ei = ti.Eaff   # Replace with actual Eaff if available
        Ej = tj.Eaff   # Replace with actual Eaff if available
        Qi = ti.Quff   # Replace with actual Quff if available
        Qj = tj.Quff   # Replace with actual Quff if available

        BO = 1 + min(npi, npj)
        ri = ti.Ruff  # Replace with actual Ruff if available
        rj = tj.Ruff  # Replace with actual Ruff if available

        # Compute rBO and rEN
        if BO > 0:
            rBO = -0.1332 * (ri + rj) * np.log(BO)
        else:
            rBO = 0.0
        denominator = Ei * ri + Ej * rj
        if denominator != 0:
            rEN = ri * rj * (np.sqrt(Ei) - np.sqrt(Ej))**2 / denominator
        else:
            rEN = 0.0
        rij = ri + rj + rBO - rEN

        # Compute kij
        kij = 28.79898 * Qi * Qj / (rij**3)  # Adjust constants as needed
        print( "kij", kij, Qi, Qj, rij )

        if self.verbosity > 1:
            print(f"bondUFF[{ti.name},{tj.name},{BO}] r={rij}({ri},{rj}|{rBO},{rEN}) k={kij}({Qi},{Qj}) E=({Ei},{Ej}) {ti.name} {tj.name} {getattr(ti, 'element', 'N/A')} {getattr(tj, 'element', 'N/A')}")

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

    def printAtomConf(self, ia, atomic_system):
        """
        Prints the configuration of a specific atom for debugging.

        Parameters:
        - ia (int): Index of the atom.
        - atomic_system (AtomicSystem): The atomic system containing all data.
        """
        print(f"Atom {ia} configuration:")
        print(f"Neighbors: {self.neighs[ia]}")
        if ia < self.nnode:
            print(f"apars: {self.apars[ia]}")
            print(f"bLs: {self.bLs[ia]}")
            print(f"bKs: {self.bKs[ia]}")
            print(f"Ksp: {self.Ksp[ia]}")
            print(f"Kpp: {self.Kpp[ia]}")
            print(f"pipos: {self.pipos[ia]}")

# IF MAIN
if __name__ == "__main__":

    from ..atomicUtils import AtomicSystem

    # Define AtomType instances with npi and ne
    AtomTypeDict = {
        "C": AtomType(name="C", Kss=300.0, Asp=109.5, Ksp=100.0, Kpp=150.0, npi=0, ne=0),
        "H": AtomType(name="H", Kss=200.0, Asp=109.5, Ksp=50.0,  Kpp=75.0,  npi=0, ne=0),
        "O": AtomType(name="O", Kss=350.0, Asp=109.5, Ksp=120.0, Kpp=180.0, npi=0, ne=0),
        "E": AtomType(name="E", Kss=0.0,   Asp=0.0,   Ksp=0.0,   Kpp=0.0,   npi=0, ne=1)  # Electron pairs
    }

    atomic_system = AtomicSystem(
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
        atomic_system=atomic_system,
        AtomTypeDict=AtomTypeDict,
        bRealloc=True,
        bEPairs=True,
        bUFF=False
    )

    # Optionally, print atom configurations for debugging
    mmff.printAtomConf(0, atomic_system)  # Replace 0 with desired atom index
