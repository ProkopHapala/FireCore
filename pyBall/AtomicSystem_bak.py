#!/usr/bin/python

from random import random
import numpy as np
from . import elements
#import elements
#import numpy as np
import copy
from . import atomicUtils as au

class AtomicSystem( ):

    def __init__(self,fname=None, apos=None, atypes=None, enames=None, lvec=None, qs=None, Rs=None, bonds=None, ngs=None, bReadN=True ) -> None:
        self.apos    = apos
        self.atypes  = atypes
        self.enames  = enames
        self.qs      = qs
        self.Rs      = Rs
        self.bonds   = bonds
        self.ngs     = ngs 
        self.lvec    = lvec
        self.aux_labels = None
        if fname is not None:
            ext = fname.split('.')[-1]
            #print( f"AtomicSystem.__init__({fname}) ext=", ext  )
            if( 'mol' == ext ):
                self.apos,self.atypes,self.enames,self.qs,self.bonds = loadMol(fname=fname, bReadN=bReadN )
            elif ( 'mol2' == ext ):
                self.apos,self.atypes,self.enames,self.qs,self.bonds, self.lvec = loadMol2(fname=fname, bReadN=bReadN )
            elif ( 'xyz' == ext ):
                self.apos,self.atypes,self.enames,self.qs, comment = load_xyz(fname=fname, bReadN=bReadN )
                if comment is not None:
                    if comment[:3] == 'lvs':      
                        self.lvec = string_to_matrix( comment, nx=3,ny=3, bExactSize=False )
                        #print( f"AtomicSystem.__init__({fname}) lvec=\n", self.lvec   )
                #print( f"AtomicSystem.__init__({fname}) comment=", comment  )
            else:
                self.apos,self.atypes,self.enames,self.qs = loadAtomsNP(fname=fname , bReadN=bReadN )

    def saveXYZ(self, fname, mode="w", blvec=True, comment="", ignore_es=None, bQs=True, other_lines=None ):
        if blvec and (self.lvec is not None):
            #print( self.lvec )
            comment= ( "lvs %6.3f %6.3f %6.3f   %6.3f %6.3f %6.3f   %6.3f %6.3f %6.3f" %(self.lvec[0,0],self.lvec[0,1],self.lvec[0,2],  self.lvec[1,0],self.lvec[1,1],self.lvec[1,2],  self.lvec[2,0],self.lvec[2,1],self.lvec[2,2]   ) ) + comment
        qs = self.qs
        if(not bQs): qs=None
        saveXYZ( self.enames, self.apos, fname, qs=qs, Rs=self.Rs, mode=mode, comment=comment, ignore_es=ignore_es, other_lines=other_lines )

    def save_mol(self, fname, title="Avogadro"):
        """
        Save the current AtomicSystem in MDL MOL V2000 format (i.e. a ".mol" file).

        The MOL file format has the following structure:
        1. A title line (up to 80 characters) – here prefixed by two spaces.
        2. A blank line.
        3. A counts line, for example:
                "  3  2  0  0  0  0  0  0  0 0999 V2000"
            where the first number is the number of atoms (right justified in 3 columns),
            the second number is the number of bonds (3 columns), followed by seven fields (each "  0"),
            then a field " 0999" and the literal " V2000".
        4. An ATOM block: one line per atom in fixed‐width format.
            In MOL V2000 the typical atom line (columns) is as follows:
            - Columns 1–3: Atom number (3-digit integer, right justified)
            - Columns 4–12: x coordinate (10.4f)
            - Columns 13–22: y coordinate (10.4f)
            - Columns 23–32: z coordinate (10.4f)
            - Columns 34–36: Atom symbol (3-character string, left justified)
            - Then 12 fields of 3 characters each (usually zeros)
        5. A BOND block: one line per bond.
            Each bond line contains:
            - Columns 1–3: Bond number (3-digit integer, right justified)
            - Columns 4–6: First atom number (3-digit integer)
            - Columns 7–9: Second atom number (3-digit integer)
            - Columns 10–12: Bond type (3-digit integer)
            - Columns 13–15: 0 (3-digit integer)
            - Columns 16–18: 0 (3-digit integer)
            - Columns 19–21: 0 (3-digit integer)
            - Columns 22–24: 0 (3-digit integer)
        6. A termination line: "M  END"

        Parameters:
        fname : str
                The output filename.
        title : str, optional
                The title for the molecule (default "Avogadro").

        Returns:
        None.
        """
        with open(fname, "w") as fout:
            # --- Title line (with two leading spaces) ---
            fout.write("  " + title + "\n")
            # --- Blank line ---
            fout.write("\n")
            n_atoms = len(self.apos)
            n_bonds = len(self.bonds) if self.bonds is not None else 0
            # --- Counts line ---
            # The counts line: atom count (3d), bond count (3d),
            # then 7 fields of "  0", then " 0999 V2000"
            counts_line = f"  {n_atoms:>3d}{n_bonds:>3d}  0  0  0  0  0  0  0 0999 V2000"
            fout.write(counts_line + "\n")
            
            # --- Atom block ---
            fout.write("\n")
            for i in range(n_atoms):
                atom_id = i + 1  # MOL format uses 1-based indexing
                x, y, z = self.apos[i]
                # Use the element name as the atom symbol.
                symbol = self.enames[i]
                # Build the atom line using fixed-width formatting.
                # Here, we format:
                #   Atom number: 3d right justified
                #   x, y, z: each 10.4f (total width 10, with 4 decimal places)
                #   Atom symbol: left aligned in 3 characters
                #   Then 12 fields of 3 characters each set to 0.
                atom_line = f"{atom_id:>3d} {x:10.4f}{y:10.4f}{z:10.4f} {symbol:<3s}" + "  0"*12
                fout.write(atom_line + "\n")
            
            # --- Bond block ---
            fout.write("\n")
            for i, bond in enumerate(self.bonds):
                bond_id = i + 1
                # Assume bond is a tuple (i, j) with 0-based indices; convert to 1-based.
                a1 = bond[0] + 1
                a2 = bond[1] + 1
                # Bond type is set to 1; then 4 fields of 0.
                bond_line = f"{bond_id:>3d}{a1:>4d}{a2:>4d}{1:>4d}" + "  0"*4
                fout.write(bond_line + "\n")
            
            # --- Termination line ---
            fout.write("M  END\n")


    def save_mol2(self, fname, comment=""):
        """
        Save the current AtomicSystem in MOL2 format.

        The MOL2 file will have the following sections:
        - @<TRIPOS>MOLECULE: a header with molecule name and counts.
        - @<TRIPOS>ATOM: one line per atom including atom id, element name,
                            coordinates, atom type, substructure id, residue name,
                            and charge.
        - @<TRIPOS>BOND: one line per bond including bond id, indices of the two
                            atoms (1–based indexing) and bond type (default "1").
        
        Parameters:
        fname   : str
                    The output filename.
        comment : str, optional
                    A comment string to include in the MOL2 file header.
                    
        Returns:
        None.
        """
        with open(fname, "w") as fout:
            # Write the MOLECULE section.
            fout.write("@<TRIPOS>MOLECULE\n")
            # Use a default molecule name or comment.
            molecule_name = "Molecule"
            fout.write(molecule_name + "\n")
            n_atoms = len(self.apos)
            n_bonds = len(self.bonds) if self.bonds is not None else 0
            # MOL2 counts: atoms, bonds, (and 0 0 0 for other fields)
            fout.write(f"{n_atoms:>3d}{n_bonds:>3d} 0 0 0\n")
            # Add required SMALL and GASTEIGER lines
            fout.write("SMALL\nGASTEIGER\n\n")
            
            # Write the ATOM section.
            fout.write("@<TRIPOS>ATOM\n")
            for i in range(n_atoms):
                atom_id   = i + 1  # MOL2 uses 1-based indexing.
                atom_type = self.enames[i]
                ename     = atom_type.replace('.','_')
                ename     = ename.split('_')[0]
                x, y, z   = self.apos[i]
                substructure = 1  # Default substructure id.
                residue = "UNL1"   # Standard residue name for unknown ligand.
                # Use self.qs if available and has the correct length, otherwise 0.0.
                charge = self.qs[i] if (self.qs is not None and len(self.qs) == n_atoms) else 0.0
                # Format: atom_id, ename, x, y, z, atom_type, substructure, residue, charge.
                fout.write("{:>7d} {:<8s} {:>9.4f} {:>9.4f} {:>9.4f} {:<5s} {:>3d}  {:<7s} {:>10.4f}\n".format(atom_id, ename, x, y, z, atom_type, substructure, residue, charge))
            
            # Write the BOND section.
            fout.write("@<TRIPOS>BOND\n")
            if self.bonds is not None:
                for i, bond in enumerate(self.bonds):
                    bond_id = i + 1
                    # bond is assumed to be a tuple (i, j) with 0-based indices.
                    # Convert to 1-based indices.
                    a1 = bond[0] + 1
                    a2 = bond[1] + 1
                    bond_type = 1
                    fout.write("{:>6d} {:>5d} {:>5d} {:>4d}\n".format(bond_id, a1, a2, bond_type ))

    def toLines(self):
        #lines = []
        #for i,pos in enumerate(self.apos):
        #    lines.append(  "%s %3.5f %3.5f %3.5f\n" %(self.enames[i], pos[0],pos[1],pos[2]) )
        return geomLines( self.apos, self.enames )

    def toXYZ(self, fout, comment="#comment", ignore_es=None, other_lines=None, bHeader=False ):
        writeToXYZ( fout, self.enames, self.apos, qs=self.qs, Rs=self.Rs, bHeader=bHeader, comment=comment, ignore_es=ignore_es, other_lines=other_lines )

    def print(self):
        print( len(self.atypes), len(self.enames), len(self.apos) )
        for i in range(len(self.apos)):
            print( "[%i] %i=%s p(%10.5f,%10.5f,%10.5f)" %( i, self.atypes[i],self.enames[i], self.apos[i,0], self.apos[i,1], self.apos[i,2] ), end =" " )
            if(self.aux_labels is not None): print(self.aux_labels[i], end =" ")
            print("")

    def getValenceElectrons( self ):
        return  np.array( [ elements.ELEMENT_DICT[e][9] for e in self.enames ] )

    def subtractValenceE(self, f0=-1.0, f=+1.0 ):
        self.qs[:] = self.qs[:]*f0 + self.getValenceElectrons()*f       

    def printBonds(self):
        print("AtomicSystem.printBonds():")
        if self.bonds is None:
            print("No bonds defined")
            return
        for i, (a, b) in enumerate(self.bonds):
            print(f"[{i}] ({a},{b}) ({self.enames[a]},{self.enames[b]})")

    def printNeighs(self):
        print("AtomicSystem.printNeighs():")
        if self.neighs is None:
            print("No neighs defined")
            return
        for i, ngi in enumerate(self.ngs):
            print(f"ngs[{i}]: ", end="")
            for j,ia in enumerate(ngi):
                print(ia, end=" ")
            print("")

    def findBonds(self, Rcut=3.0, RvdwCut=1.5, RvdWs=None, byRvdW=True ):
        if self.atypes is None:
            self.atypes = [ elements.ELEMENT_DICT[e][0] for e in self.enames ]
        self.bonds, rs = findBondsNP( self.apos, self.atypes, Rcut=Rcut, RvdwCut=RvdwCut, RvdWs=RvdWs, byRvdW=byRvdW )
        return self.bonds, rs

    def findHBonds(self, Rb=1.5, Rh=2.5, angMax=60.0, typs1={"H"}, typs2=neg_types_set, bPrint=False, bHbase=False ):
        return findHBondsNP( self.apos, atypes=self.enames, Rb=Rb, Rh=Rh, angMax=angMax, typs1=typs1, typs2=typs2, bPrint=bPrint,  bHbase=bHbase )

    def findBondsOfAtom(self, ia, bAtom=False ):
        if bAtom: 
            return [ b[1] for b in self.bonds if(b[0]==ia) ] + [ b[0] for b in self.bonds if(b[1]==ia) ] 
        else:
            return [i for i,b in enumerate(self.bonds) if (b[0]==ia) or (b[1]==ia) ]

    def neighs( self, bBond=True ):
        if(self.bonds is None):
            self.findBonds()
        self.ngs = neigh_bonds( len(self.apos), self.bonds )
        return self.ngs

    def find_groups(self):
        if self.ngs is None: self.neighs()
        ngs = self.ngs
        #print( ngs )
        groups = { }
        for inod in range(len(self.apos)):
            if len(ngs[inod]) > 1: groups[inod] = [inod]
        for inod,g in groups.items():
            inod = g[0] 
            g += [ ia for ia in ngs[inod].keys() if ia not in groups ] 
        return groups

    def select_by_ename( self, elist ):
        return [ i for i,e in enumerate(self.enames) if e in elist ]

    def getNeighsOfType( self, selection, typ='N'):
        if self.ngs is None: self.neighs()
        return findNeighsOfType( selection, self.enames, self.ngs, typ=typ ) 

    def select_by_neighType( self, neighs, typ='N', neighTyps={'H':(1,2)} ):
        return findTypeNeigh_( self.enames, neighs, typ=typ, neighTyps=neighTyps )

    def findAngles(self, select=None, ngs=None, ):
        if ngs is None:
            ngs = self.neighs()
        return findAngles( self.apos, select=select, neighs=ngs )

    def findDihedral( self, select=None, ngs=None, neighTyp={'H'} ):
        if ngs is None:
            ngs = self.neighs()
        return findDihedral( self.apos, self.enames, ngs, select=select, neighTyp=neighTyp ) 

    def findCOG(self, apos, byBox=False ):
        return findCOG( apos, byBox=byBox )
    
    def projectAlongBondDir( self, i0, i1 ):
        return projectAlongBondDir( self.apos, i0, i1 )

    def store_bond_lengths(self):
        bond_lengths = {}
        bonds = self.findBonds()  # Get all bonds in the system
        for bond in bonds:
            ia,ja = bond
            if ia>ja: ia,ja = ja,ia
            length = np.linalg.norm(self.apos[ia]-self.apos[ja])
            bond_lengths[(ia,ja)] = length
        self.bond_legths = bond_lengths
        return bond_lengths

    def restore_bond_length(self, ij, L=None ):
        ia,ja= ij
        d = self.apos[ja] - self.apos[ia]
        Lnow = np.linalg.norm(d)
        if L is None:
            if ia>ja: i,j = ja,ia
            else:     i,j = ia,ja
            L = self.bond_lengths[(i,j)]
        f = L / Lnow
        self.apos[ia] = self.apos[ja] + d * f


    def clonePBC(self,nPBC=(1,1,1) ):
        nx,ny,nz= nPBC
        nxyz=nx*ny*nz
        na = len(self.apos)
        apos   = np.zeros((na*nxyz,3))
        #print( "clonePBC ", na, len(self.atypes) )
        if self.atypes is not None: 
            atypes = np.zeros(na*nxyz,np.int32)
        else:
            atypes = None

        if self.enames is not None: 
            enames = []
        else:
            enames = None

        if self.qs is not None: 
            qs = np.zeros(na*nxyz) 
        else:
            qs = None

        #print( nxyz, na, apos.shape, atypes.shape )
        if( nxyz > 1 ):
            lvec   = np.array([ self.lvec[0,:]*nx,self.lvec[1,:]*ny,self.lvec[2,:]*nz ]) 
            i0=0
            for iz in range(nz):
                for iy in range(ny):
                    for ix in range(nx):
                        shift = self.lvec[0,:]*ix  + self.lvec[1,:]*iy + self.lvec[2,:]*iz
                        apos  [i0:i0+na,:] = self.apos[:,:] + shift[None,:]
                        if atypes is not None: atypes[i0:i0+na  ] = self.atypes
                        if qs     is not None: qs    [i0:i0+na  ] = self.qs    
                        if enames is not None: enames[i0:i0+na  ] = self.enames
                        #if enames is not None: enames += self.enames
                        i0+=na
        else:
            lvec=self.lvec
            apos  [:,:] = self.apos[:,:]
            if atypes is not None: atypes[:] = self.atypes[:]
            if qs     is not None: qs    [:] = self.qs    [:]  
            if enames is not None: enames[:] = self.enames[:]

        return AtomicSystem(apos=apos, atypes=atypes, enames=enames, lvec=lvec, qs=qs ) 

    def symmetrized(self, d=0.1 ):
        # def atoms_symmetrized( atypes, apos, lvec, qs=None, REQs=None, d=0.1):
        atypes, apos, qs, REQs, ws = atoms_symmetrized( self.atypes, self.apos, self.lvec, qs=self.qs, d=d );
        enames = iz2enames( atypes )
        return AtomicSystem( apos=apos, atypes=atypes, enames=enames, lvec=self.lvec.copy(), qs=qs ), ws 

    def selectSubset(self, inds ):
        if self.atypes is not None: 
                atypes = self.atypes[inds]
        else:
            atypes = None

        if self.enames is not None: 
            enames = [ self.enames[i] for i in inds ]
        else:
            enames = None

        if self.qs is not None: 
            qs = self.qs[inds]
        else:
            qs = None

        lvec=self.lvec
        apos  = self.apos[inds,:]

        return AtomicSystem(apos=apos, atypes=atypes, enames=enames, lvec=lvec, qs=qs ) 

    def selectBondedCluster( self, s ):
        na = len(self.apos)
        if self.bonds is None: self.findBonds()
        s     = selectBondedCluster( s, self.bonds )
        ins  = [ i for i in range(na) if (i in s) ]
        outs = [ i for i in range(na) if (i not in s) ] 
        return ins,outs

    def makeRotMat( self, ip1, ip2, _0=1 ):
        fw  = self.apos[ip1[1]-_0]-self.apos[ip1[0]-_0]
        up  = self.apos[ip2[1]-_0]-self.apos[ip2[0]-_0]
        return makeRotMat( fw, up )

    def orient_mat(self, rot, p0=None, bCopy=False ):
        apos=self.apos  
        if(bCopy): apos=apos.copy()
        if p0  is not None: apos[:,:]-=p0[None,:]
        if rot is not None: mulpos( apos, rot )
        return apos

    def orient_vs(self, fw, up, p0=None, trans=None, bCopy=False ):
        if fw is None:
            rot = None
        else:
            rot = makeRotMat( fw, up )
            if trans is not None: rot=rot[trans,:]
        return self.orient_mat( rot, p0, bCopy )

    def orient( self, i0, b1, b2, _0=1, trans=None, bCopy=False ):
        #print( "orient i0 ", i0, " ip1 ", ip1, " ip2 ",ip2 )
        # p0  = self.apos[i0-_0]
        # fw  = self.apos[ip1[1]-_0]-self.apos[ip1[0]-_0]
        # up  = self.apos[ip2[1]-_0]-self.apos[ip2[0]-_0]
        p0, fw, up = makeVectros( self.apos, i0, b1, b2, _0=_0 )
        return self.orient_vs( fw, up, p0, trans=trans, bCopy=bCopy )
    
    def orientPCA(self, perm=None):
        orientPCA(self.apos, perm=perm )

    def shift(self, vec, sel=None ):
        if sel is None: 
            self.apos[:,0] += vec[0]
            self.apos[:,1] += vec[1]
            self.apos[:,2] += vec[2]
        else:
            self.apos[sel,0] += vec[0]
            self.apos[sel,1] += vec[1]
            self.apos[sel,2] += vec[2]

    def rotate_ax(self, ang, ax=(0,1), p0=None ):
        rot = makeRotMatAng( ang, ax=ax )
        if p0  is not None: self.apos[:,:]-=p0[None,:]
        mulpos( self.apos, rot )
        if p0  is not None: self.apos[:,:]+=p0[None,:]

    def delete_atoms(self, lst ):
        st = set(lst)
        if( self.apos   is not None ): self.apos   =  np.delete( self.apos,   lst, axis=0 )
        if( self.atypes is not None ): self.atypes =  np.delete( self.atypes, lst )
        if( self.qs     is not None ): self.qs     =  np.delete( self.qs,     lst )
        if( self.Rs     is not None ): self.Rs     =  np.delete( self.Rs,     lst )
        if( self.enames is not None ): self.enames =  np.delete( self.enames, lst )
        if( self.aux_labels is not None ): self.aux_labels = [ v for i,v in enumerate(self.aux_labels) if i not in st ] 


    def preinitialize_atomic_properties(self):
        """
        Preinitialize per-atom arrays for an AtomicSystem.
        
        This function assumes that the system’s atypes (or enames) have been set.
        It uses the global 'elements.ELEMENTS' (a list of lists) to set default values:
        - qs: set to the element’s default valence electron count (column index 9)
        - Rs: set to the element’s van der Waals radius (column index 7)
        - aux_labels: set to a default label (simply the atom’s index as a string)
        
        Parameters:
        atomicSystem (AtomicSystem): an instance of AtomicSystem.
        
        Raises:
        ValueError: if atomicSystem.atypes is not defined.
        """
        natoms = len(self.apos)
        if self.atypes is None:   raise ValueError("The system does not have atypes defined. Please initialize the system’s atypes (or enames) first.")
        if self.qs is None:   # Assume atypes is an array of atomic numbers (e.g. 6 for carbon, etc.)
            qs = []
            for z in self.atypes:   # our ELEMENTS list is zero-based: for atomic number z, use ELEMENTS[z-1]
                qs.append(elements.ELEMENTS[z-1][9])
            self.qs = np.array(qs)
        if self.Rs is None:  # For each atom, use the vdW radius (column index 7)
            Rs = []
            for z in self.atypes: Rs.append(elements.ELEMENTS[z-1][7])
            self.Rs = np.array(Rs)
        # Initialize aux_labels if not defined.
        if self.aux_labels is None: self.aux_labels = [str(i) for i in range(natoms)]
        self.neighs()
        # (If you have other arrays you want to preinitialize, do it here.)
        #print(f"Pre-initialized atomic properties for {natoms} atoms.")

        
    def check_atomic_properties(atomicSystem):
        """
        Check that the per-atom arrays (qs, Rs, aux_labels) are defined and
        have the correct length. If not, raise an error telling the user
        to run preinitialize_atomic_properties().
        """
        natoms = len(atomicSystem.apos)
        if (atomicSystem.qs is None or len(atomicSystem.qs) != natoms or
            atomicSystem.Rs is None or len(atomicSystem.Rs) != natoms or
            atomicSystem.aux_labels is None or len(atomicSystem.aux_labels) != natoms):
            raise ValueError("Not all per-atom arrays are initialized correctly. Please call preinitialize_atomic_properties() on your system.")
                            
                            
    # Example: modify append_atoms() to check rather than auto-initialize
    def append_atoms(self, B, pre="A"):
        # Ensure both systems have been pre-initialized:
        self.check_atomic_properties()
        B.check_atomic_properties()
        
        # Number of atoms in self and in B
        nA = len(self.apos)
        nB = len(B.apos)
        
        self.apos   = np.append(self.apos,   B.apos, axis=0)
        self.atypes = np.append(self.atypes, B.atypes)
        self.qs     = np.append(self.qs,     B.qs)
        self.Rs     = np.append(self.Rs,     B.Rs)
        self.enames = np.append(self.enames, B.enames)
        
        # Extend the aux_labels list:
        self.aux_labels.extend(B.aux_labels)


    def remap( self, lst ):
        dct = {   key:value for (value,key) in enumerate(self.aux_labels) }
        return [ dct.get(key,-1) for key in lst ]


    def attach_group( self, G,  i0, i1, iup,   bond,  up=(0.,0.,1.),  _0=1, pre="A"  ): 
        """
        Attach an end–group (G) to the backbone (self) at a specified bond.
        
        The attachment is done in two steps:
        1. **Internal Orientation of the Group:**  
            The group is reoriented in its own frame by calling:
                G.orient(i0, (i0, i1), iup, _0=_0)
            - *i0*: the index (or indices) for the pivot atom in the group. This atom
                    is moved to the attachment position.
            - *(i0, i1)*: a tuple defining a bond in the group that determines the
                        forward (direction) vector. The forward vector is computed as
                        the difference between the positions of the atom at i1 and i0.
                        The atom corresponding to i1 is then deleted (replaced) in the group.
            - *iup*: a tuple (or list) of two indices that defines the up vector in the group.
                    The up vector is computed (typically as the difference between the
                    positions of the atoms provided) and is used to fix the rotation about
                    the forward axis.
        
        2. **Alignment to the Backbone:**  
            The backbone provides the attachment bond (given by `bond`) and a backbone
            up vector (given by `up`). A rotation matrix is computed with:
                rot = rotmat_from_points(self.apos, ifw=bond, up=up, _0=_0)
            This matrix aligns the backbone’s forward vector (computed from the bond) with
            the group’s forward vector. The group is then rotated by this matrix (via
                G.orient_mat(rot)
            ) and translated so that the pivot atom of the group coincides with the backbone’s
            attachment position.
        
        Parameters:
        G      : AtomicSystem
                The end–group to attach. It must have its atoms pre‐oriented as per the
                expected coordinate system.
        i0     : int or iterable
                The index (or indices) of the pivot atom in G (1-based indexing if _0=1).
        i1     : int
                The index (1-based) of the atom in G used to define the forward vector.
                This atom will be removed after orientation.
        iup    : tuple (i_up0, i_up1)
                A pair of indices (1-based) in G whose difference defines the up vector.
        bond   : tuple (i_backbone1, i_backbone2)
                A pair of atom indices (1-based) in the backbone that define the bond where
                the end–group is attached. The forward vector on the backbone is computed as
                the vector from i_backbone1 to i_backbone2.
        up     : 3-tuple or array, optional (default=(0.,0.,1.))
                The up vector for the backbone. This is used to fix the rotation about the
                forward axis.
        _0     : int, optional (default=1)
                An offset to account for whether the provided indices are 0-based or 1-based.
        pre    : str, optional (default="A")
                A prefix for labeling the atoms that come from the group.
        
        After executing, the group G is reoriented, rotated, and translated so that its
        pivot atom is placed at the backbone’s attachment site. The atom used for forward
        definition (i1) is deleted, and the group’s atoms are appended to the backbone.
        """
        up  = np.array( up )
        rot = rotmat_from_points( self.apos, ifw=bond, up=up, _0=1 );   
        rot = rot.transpose()
        p0  = self.apos[bond[0]-_0]
        
        if( G.aux_labels is None ): G.aux_labels = [ pre+str(i) for  i in range(len(G.apos)) ]

        G.orient( i0,(i0,i1),iup, _0=_0 )
        G.orient_mat( rot ); 
        G.apos[:,:]+=p0[None,:]
        G.delete_atoms( [i1-_0] )

        self.append_atoms( G )

    def reindex_bonds(self, old_to_new_map, to_remove=None ):
        #print ("self.reindex_bonds: old_to_new_map \n", old_to_new_map)
        #print ("self.reindex_bonds: to_remove      \n", to_remove)
        self.bonds = reindex_bonds( self.bonds, old_to_new_map, to_remove )
        self.ngs   = None

    def extract_marker_pairs(self, markerX, markerY, remove=True):
        """Legacy method that combines finding marker pairs without removal."""
        pairs = self.find_marker_pairs(markerX, markerY)
        return pairs

    def find_marker_pairs(self, markerX, markerY):
        """
        Find marker pairs in this system based on element types and bonding information.
        For each atom with element equal to markerX, look for a bonded neighbor with element equal to markerY.
        Returns a list of tuples (iX, iY) where iX is the index of a markerX atom and iY is the index of its bonded markerY neighbor.
        """
        if self.ngs is None:
            self.neighs(bBond=True)
        mks = []
        for i, ename in enumerate(self.enames):
            if ename == markerX:    # pivot atom marker-X
                ngi = self.ngs[i]    
                for j in ngi:
                    if self.enames[j] == markerY:
                        i2=j   # index of a markerY-typed bonded to the pivot atom
                    else:
                        i3=j   # atom to which pivot atom is bonded, but is not a marker (i.e. Anchor-atom)  
                mks.append( (i, i2, i3))                
        return mks

    def ensure_numpy_arrays(self):
        """Ensure position arrays are numpy arrays."""
        if not isinstance(self.apos, np.ndarray):
            self.apos = np.array(self.apos)
        
    def filter_system(self, mask, bInverted=False):
        """Create a new filtered system without marker atoms.
        Parameters:
            mask : set Set of atom indices to keep ( or remove if inverted)
            bInverted : bool Invert the mask
        Returns:
            tuple : (filtered arrays, old_to_new index_map)
        """
        n = len(self.apos)
        #if bInverted: mask = set(range(len(self.apos))).difference(mask)
        if bInverted: mask = [i for i in range(n) if i not in mask]
        #print ("filter_system: mask ", mask)
        #keep_mask = np.array([i not in to_remove for i in range(len(self.apos))])
        # filtered = {
        #     'apos':   self.apos  [mask],
        #     'atypes': self.atypes[mask],
        #     'enames': self.enames[mask],
        #     'qs':     self.qs    [mask] if self.qs is not None else None,
        #     'Rs':     self.Rs    [mask] if self.Rs is not None else None,
        #     'aux_labels': [label for i, label in enumerate(self.aux_labels) if i not in mask] if self.aux_labels is not None else None
        # }

        old_to_new = make_reindex( n, mask, bInverted=False)    

        if self.bonds is not  None:
            bonds = [ b for b in self.bonds if (b[0] in mask) and (b[1] in mask) ]
            #print( "reindex_bonds() bonds \n", bonds )
            #bonds = [ (old_to_new[b[0]], old_to_new[b[1]]) for b in bonds if b[0] in old_to_new and b[1] in old_to_new ]
            bonds = [ (old_to_new[b[0]], old_to_new[b[1]]) for b in bonds  ]
        else:
            bonds = None

        filtered = AtomicSystem(
            apos  =self.apos  [mask], 
            atypes=self.atypes[mask],
            enames=self.enames[mask],
            lvec  =self.lvec,
            qs    =self.qs[mask] if self.qs is not None else None,
            Rs    =self.Rs[mask] if self.Rs is not None else None,
            #ngs   =self.ngs[mask] if self.ngs is not None else None,
            bonds = bonds
        )
          
        return filtered, old_to_new
        
    #def merge_arrays(self, other, other_bonds, offset):
    def merge_arrays(self, other, offset=None):
        """Merge arrays from other system into self.
        
        Parameters:
            other_arrays : dict
                Filtered arrays from other system
            other_bonds : ndarray
                Reindexed bonds from other system
            offset : int
                Offset for bond indices
        """
        # Merge main arrays
        # self.apos   = np.concatenate([self.apos,   other_arrays['apos'  ]], axis=0)
        # self.atypes = np.concatenate([self.atypes, other_arrays['atypes']])
        # self.enames = np.concatenate([self.enames, other_arrays['enames']])
        # if self.qs         is not None and other_arrays['qs']         is not None: self.qs         = np.concatenate([self.qs, other_arrays['qs']])
        # if self.Rs         is not None and other_arrays['Rs']         is not None: self.Rs         = np.concatenate([self.Rs, other_arrays['Rs']])
        # if self.aux_labels is not None and other_arrays['aux_labels'] is not None: 
        #     self.aux_labels = np.concatenate([self.aux_labels, other_arrays['aux_labels']])
        # else:
        #     self.aux_labels = None
            
        # # Merge bonds
        # if other_bonds is not None:
        #     adjusted_bonds = other_bonds + offset
        #     if self.bonds is not None:
        #         self.bonds = np.array(self.bonds)
        #         self.bonds = np.concatenate([self.bonds, adjusted_bonds], axis=0)
        #     else:
        #         self.bonds = adjusted_bonds

        if offset is None: offset = len(self.apos)

        self.apos   = np.concatenate([self.apos,   other.apos], axis=0)
        self.atypes = np.concatenate([self.atypes, other.atypes])
        self.enames = np.concatenate([self.enames, other.enames])
        if self.qs         is not None and other.qs         is not None: self.qs         = np.concatenate([self.qs, other.qs])
        if self.Rs         is not None and other.Rs         is not None: self.Rs         = np.concatenate([self.Rs, other.Rs])
        if self.aux_labels is not None and other.aux_labels is not None: 
            self.aux_labels = np.concatenate([self.aux_labels, other.aux_labels])
        else:
            self.aux_labels = None
            
        # Merge bonds
        if other.bonds is not None:
            adjusted_bonds = np.array(other.bonds) + offset
            if self.bonds is not None:
                self.bonds = np.array(self.bonds)
                self.bonds = np.concatenate([self.bonds, adjusted_bonds], axis=0)
            else:
                self.bonds = adjusted_bonds

    def add_bond(self, b):
        #print ("add_bond", b)
        #print( "bonds.shape", self.bonds.shape )
        self.bonds = np.concatenate((self.bonds, np.array([b])), axis=0)

    def merge_geometries(self, other, group_mk, backbone_mk ):
        """Merge another AtomicSystem into this one using the provided group marker pair for alignment.
        
        This implementation appends the atoms and bonds from 'other' into self, adjusting indices appropriately.
        The process follows these steps:
        1. Find neighbors of marker atoms in both systems
        2. Remove marker atoms from the group system
        3. Merge the remaining atoms and bonds
        4. Create new bonds between fragments based on marker neighbors
        
        Parameters:
            other : AtomicSystem
                   The system to merge into this one
            group_marker_pair : tuple
                   (iX, iY) marker pair from the group system used for alignment
        """
        # Ensure numpy arrays
        self.ensure_numpy_arrays()
        other.ensure_numpy_arrays()
        
        removed = set(group_mk[:2])
        other_filtered, old_to_new = other.filter_system( removed, bInverted=True            )     # Filter group system without markers 
        #other_bonds                = other.reindex_removed_bonds(other.bonds, removed, old_to_new)   # Reindex group bonds
        #other_bonds = reindex_bonds( other.bonds, old_to_new, to_remove=removed )
        
        # Merge arrays with offset
        offset = len(self.apos)
        #self.merge_arrays(other_filtered, other_bonds, offset)
        self.merge_arrays(other_filtered, offset)

        # # Create new bonds between fragments
        # self.create_fragment_bonds(backbone_neighbors, group_neighbors, old_to_new, offset)

        #self.create_bond_reindexed( (group_mk[2],backbone_mk[2]), old_to_new, offset )
        #print( "old_to_new", old_to_new )
        i2 = old_to_new[group_mk[2]]

        #print( "-----BEFORE self.bonds", self.bonds )
        self.add_bond( (backbone_mk[2], i2 + offset) )
        #print( "-----AFTER self.bonds", self.bonds )

        # Clear neighbor list since it needs to be rebuilt
        self.ngs = None

    def compute_group_orientation(self, G, backbone_pair, group_pair, _0=1):
        """Compute the orientation transformation for attaching a group to the backbone.
        
        Parameters:
            G : AtomicSystem
                The group to orient
            backbone_pair : tuple
                (iX, iY) indices of marker atoms in backbone
            group_pair : tuple
                (iX, iY) indices of marker atoms in group
            _0 : int
                Offset for index conversion
                
        Returns:
            tuple: (R, X_b, A2)
                R - rotation matrix
                X_b - translation point (backbone marker position)
                A2 - attachment point in group (for translation)
        """
        iX_b, iY_b,_  = backbone_pair
        iX_g, iY_g,_ = group_pair
        
        # Compute frames for backbone and group
        X_b, A1, M_target = compute_attachment_frame_from_indices(self.apos, iX_b, iY_b, self, bFlipFw=False, _0=_0)
        X_g, A2, M_group = compute_attachment_frame_from_indices(G.apos, iX_g, iY_g, G, bFlipFw=True, _0=_0)
        
        # Compute rotation matrix R = M_target @ (M_group)ᵀ
        R = M_target @ M_group.T
        
        return R, X_b, A2

    def delete_atoms(self, to_remove):
            rem = sorted(to_remove, reverse=True)
            for idx in rem:
                self.apos   = np.delete(self.apos,   idx, axis=0)
                self.atypes = np.delete(self.atypes, idx)
                self.enames = np.delete(self.enames, idx)
                if self.qs is not None:
                    self.qs = np.delete(self.qs, idx)
                if self.Rs is not None:
                    self.Rs = np.delete(self.Rs, idx)
                if self.aux_labels is not None:
                    self.aux_labels = np.delete(self.aux_labels, idx)

    def attach_group_by_marker(self, G, markerX="Xe", markerY="He", _0=1):
        """Attach an end–group G to this backbone using marker atoms and connectivity.
        Steps:
          1. Find marker pairs in both the backbone and the group.
          2. Ensure the group has exactly one marker pair.
          3. Orient and transform the group.
          4. Merge the transformed group.
          5. Remove marker atoms from the backbone.
          6. Update the neighbor list.
        """
        # 1. Find marker pairs
        backbone_inds = self.find_marker_pairs(markerX, markerY)
        group_inds    = G.find_marker_pairs(markerX, markerY)
        #print( "backbone_inds ", backbone_inds )
        #print( "group_inds    ", group_inds )
        if not backbone_inds:    raise ValueError(f"No marker pair ({markerX}, {markerY}) found in backbone")
        if len(group_inds) != 1: raise ValueError(f"Group must have exactly one marker pair, found {len(group_inds)}")
            
        # 2. Get orientation transformation
        R, X_b, A2 = self.compute_group_orientation(G, backbone_inds[0], group_inds[0], _0)
        
        # 3. Make a deep copy of G and transform it
        G_copy = copy.deepcopy(G)
        G_copy.apos = (R @ (G_copy.apos - A2).T).T + X_b
        
        # 4. Merge the transformed geometries
        #self.merge_geometries(G_copy, group_inds[0], backbone_inds[0] )
        
        # 5. Remove marker atoms from the backbone

        gind = group_inds[0]
        bind = backbone_inds[0]

        for ii, bind in enumerate(backbone_inds):

            bind = self.find_marker_pairs(markerX, markerY)[0]

            # 2. Get orientation transformation
            R, X_b, A2 = self.compute_group_orientation(G, bind, gind, _0)
            
            # 3. Make a deep copy of G and transform it
            G_copy = copy.deepcopy(G)
            G_copy.apos = (R @ (G_copy.apos - A2).T).T + X_b

            self.merge_geometries(G_copy, gind, bind )
            to_remove = set( [bind[0], bind[1]] )

            old_to_new = {}
            new_idx = 0
            for old_idx in range(len(self.apos)):
                if old_idx not in to_remove:
                    old_to_new[old_idx] = new_idx
                    new_idx += 1
            self.reindex_bonds(old_to_new, to_remove)

            self.delete_atoms( to_remove )
            #rem = sorted(to_remove, reverse=True)
            # for idx in rem:
            #     self.apos   = np.delete(self.apos,   idx, axis=0)
            #     self.atypes = np.delete(self.atypes, idx)
            #     self.enames = np.delete(self.enames, idx)
            #     if self.qs is not None:
            #         self.qs = np.delete(self.qs, idx)
            #     if self.Rs is not None:
            #         self.Rs = np.delete(self.Rs, idx)
            #     if self.aux_labels is not None:
            #         self.aux_labels = np.delete(self.aux_labels, idx)

        # 6. Update neighbor list
        self.neighs(bBond=True)