
## Task

* I'm trying to move my C/C++ OpenCL interface to python pyOpenCL interface fro Molecular dynamics simulations on GPU.
* Our approach is a bit special because it simulates multiple replicas of the same system, at the same time in order to exploit the GPU parallelism efficintly (the systems are small like 100 atoms, an we have like 100-1000 replicas, using >10000 GPU cores).
* Therefore we shoudl do the following steps:
  1. load atomic system from .xyz file using AtomicSystem() class implemented in atomicUtils.py and create bodning topology for it 
  2. process AtomicSystem() object to produce MMFF parameters

Now I will give you few relevent pieces of code from which you can better understand what exactly should be implemented. 

* You should read it and realize what needs to be implemented in the new python class MolecularDynamics.py using py OpenCL.
* Perhaps as intermediante we should create new class MMFF which implements MMFFforcefield class.

## class AtomicSystem in atomicUtils.py

This class we will use to read a molecule from .xyz file and build bodning topology and neighbor list

```python
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
            if( 'xyz' == ext ):
                self.apos,self.atypes,self.enames,self.qs, comment = load_xyz(fname=fname, bReadN=bReadN )
                if comment is not None:
                    if comment[:3] == 'lvs':      
                        self.lvec = string_to_matrix( comment, nx=3,ny=3, bExactSize=False )
                        #print( f"AtomicSystem.__init__({fname}) lvec=\n", self.lvec   )
                #print( f"AtomicSystem.__init__({fname}) comment=", comment  )
            else:
                self.apos,self.atypes,self.enames,self.qs = loadAtomsNP(fname=fname , bReadN=bReadN )

    def findBonds(self, Rcut=3.0, RvdwCut=1.2, RvdWs=None, byRvdW=True ):
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
```

## MMFFsp3_loc 

This is C++ class which implements MMFFforcefield. The data arrays correspond to the data arrays passed to OpenCL kernels. However to OpenCL krenel we pack multiple replicas of the system into the same buffer.

```C++
class MMFFsp3_loc : public NBFF { public:
    static constexpr const int nneigh_max = 4; // maximum number of neighbors

    // === inherited from NBFF
    // int natoms=0;        // [natoms] // from Atoms
    //int   * atypes  =0; // [natom]  atom types
    //Vec3d * apos  =0;     // [natoms] // from Atoms
    //Vec3d * vapos = 0;    // [natom]  velocities of atoms
    //Vec3d * fapos  =0;    // [natoms] // from NBFF
    //Quat4i* neighs =0;    // [natoms] // from NBFF
    //Quat4i* neighCell=0;  // [natoms] // from NBFF
    //Quat4d* REQs =0;      // [nnode]  // from NBFF
    //bool    bPBC=false;       // from NBFF
    //Vec3i   nPBC;         // from NBFF 
    //Mat3d   lvec;         // from NBFF
    //double  Rdamp  = 1.0; // from NBFF

    // dimensions of the system
    int  nDOFs=0,nnode=0,ncap=0,nvecs=0,ntors=0;
    double Etot,Eb,Ea, Eps,EppT,EppI;  // total energy, bond energy, angle energy, pi-sigma energy, pi-pi torsion energy, pi-pi interaction energy

    double *  DOFs = 0;   // degrees of freedom
    double * fDOFs = 0;   // forces

    bool doBonds  =true; // compute bonds
    bool doNeighs =true; // compute neighbors
    bool doPiPiI  =true; // compute pi-pi interaction parallel
    bool doPiPiT  =true; // compute pi-pi torsion     perpendicular
    bool doPiSigma=true; // compute pi-sigma interaction
    bool doAngles =true; // compute angles
    //bool doEpi    =true; // compute pi-electron interaction

    // bool    bCollisionDamping        = false; // if true we use collision damping
    // bool    bCollisionDampingAng     = false;
    // bool    bCollisionDampingNonBond = false;  // if true we use collision damping for non-bonded interactions
    // double  damping_medium           = 1.0;   // cdamp       = 1 -(damping_medium     /ndampstep     )
    // double  collisionDamping         = 0.0;   // col_damp    =     collisionDamping   /(dt*ndampstep )
    // double  collisionDamping_ang     = 0.0;   // col_damp_NB =     collisionDamping_NB/(dt*ndampstep )
    // double  collisionDamping_NB      = 0.0;   // col_damp_NB =     collisionDamping_NB/(dt*ndampstep )
    // int     ndampstep                = 10;    // how many steps it takes to decay velocity to to 1/e of the initial value
    // double  col_damp_dRcut1          = -0.2;   // non-covalent collision damping interaction goes between 1.0 to 0.0 on interval  [ Rvdw , Rvdw+col_damp_dRcut ]
    // double  col_damp_dRcut2          =  0.3;   // non-covalent collision damping interaction goes between 1.0 to 0.0 on interval  [ Rvdw , Rvdw+col_damp_dRcut ]
    // double col_damp      = 0.0;  //  collisionDamping   /(dt*ndampstep );
    // double col_damp_NB   = 0.0;  //  collisionDamping_NB/(dt*ndampstep );
    // double col_dampAng   = 0.0;  //  collisionDamping   /(dt*ndampstep );

    CollisionDamping colDamp;

    Vec3d cvf=Vec3dZero; // <f|f>, <v|v>, <v|f> 

    bool bEachAngle = false; // if true we compute angle energy for each angle separately, otherwise we use common parameters for all angles
    bool bTorsion   = false; // if true we compute torsion energy
    
    //                           c0     Kss    Ksp    c0_e
    Quat4d default_NeighParams{ -1.0,   1.0,   1.0,   -1.0 }; // default parameters for neighbors, c0 is cosine of equilibrium angle, Kss is bond stiffness, Ksp is pi-sigma stiffness, c0_e is cos of equilibrium angle for pi-electron interaction

    // Dynamical Varaibles;
    //Vec3d *   apos=0;   // [natom]
    //Vec3d *  fapos=0;   // [natom]
    Vec3d *  pipos=0;   // [nnode]
    Vec3d * fpipos=0;   // [nnode]

    // Aux Dynamil
    Vec3d * fneigh  =0;  // [nnode*4]     temporary store of forces on atoms form neighbors (before assembling step)
    Vec3d * fneighpi=0;  // [nnode*4]     temporary store of forces on pi    form neighbors (before assembling step)

    //Quat4i*  neighs =0;   // [natoms] // from NBFF
    Quat4i*  bkneighs=0;   // [natoms]  inverse neighbors
    
    Quat4d*  apars __attribute__((aligned(64))) =0;  // [nnode] angle parameters
    Quat4d*  bLs   __attribute__((aligned(64))) =0;  // [nnode] bond lengths
    Quat4d*  bKs   __attribute__((aligned(64))) =0;  // [nnode] bond stiffness
    Quat4d*  Ksp   __attribute__((aligned(64))) =0;  // [nnode] stiffness of pi-alignment
    Quat4d*  Kpp   __attribute__((aligned(64))) =0;  // [nnode] stiffness of pi-planarization

    Vec3d*   angles=0; // [nnode*6]  angles between bonds

    Quat4i*  tors2atom  __attribute__((aligned(64))) =0; // [ntors]  torsion atoms
    Quat4d*  torsParams __attribute__((aligned(64))) =0; // [ntors]  torsion parameters

    Quat4d* constr  __attribute__((aligned(64))) = 0; // [natom]  constraints
    Quat4d* constrK __attribute__((aligned(64))) = 0; // [natom]  constraints
    //Vec3d * vapos  __attribute__((aligned(64))) = 0; // [natom]  velocities of atoms

    Mat3d   invLvec; // inverse lattice vectors

    bool    bAngleCosHalf         = true;   // if true we use evalAngleCosHalf() instead of evalAngleCos() to compute anglular energy
    // these are defined in ForceFiled.h
    //bool    bSubtractAngleNonBond = false;  // if true we subtract angle energy from non-bonded energy
    //bool    bSubtractBondNonBond  = false;  // if true we subtract bond energy from non-bonded energy

    //int itr_DBG=0;

// =========================== Functions

// reallcoate MMFFsp3_loc
void realloc( int nnode_, int ncap_, int ntors_=0 ){
    nnode=nnode_; ncap=ncap_; ntors=ntors_;
    natoms= nnode + ncap; 
    nvecs = natoms+nnode;  // each atom as also pi-orientiation (like up-vector)
    nDOFs = nvecs*3;
    //printf( "MMFFsp3::realloc() natom(%i,nnode=%i,ncap=%i), npi=%i, nbond=%i \n", natoms, nnode, ncap, npi, nbonds );
    int ipi0=natoms;
    
    _realloc0(  DOFs    , nDOFs , (double)NAN );
    _realloc0( fDOFs    , nDOFs , (double)NAN );
    apos   = (Vec3d*) DOFs ;
    fapos  = (Vec3d*)fDOFs;
    pipos  = apos  + ipi0;
    fpipos = fapos + ipi0;
    // ---- Aux
    _realloc0( fneigh  , nnode*4, Vec3dNAN );
    _realloc0( fneighpi, nnode*4, Vec3dNAN );
    // ----- Params [natom]
    _realloc0( atypes    , natoms, -1 );
    _realloc0( neighs    , natoms, Quat4iMinusOnes );
    _realloc0( neighCell , natoms, Quat4iMinusOnes );
    _realloc0( bkneighs  , natoms, Quat4iMinusOnes);
    _realloc0( apars     , nnode, Quat4dNAN );
    _realloc0( bLs       , nnode, Quat4dNAN );
    _realloc0( bKs       , nnode, Quat4dNAN );
    _realloc0( Ksp       , nnode, Quat4dNAN );
    _realloc0( Kpp       , nnode, Quat4dNAN );

    // Additional:
    // Angles
    _realloc0( angles, nnode*6, Vec3dNAN );   // 6=4*3/2
    // Torsions
    _realloc0( tors2atom,  ntors, Quat4iZero );
    _realloc0( torsParams, ntors, Quat4dNAN  ); 

    _realloc0( constr    , natoms, Quat4dOnes*-1. );
    _realloc0( constrK   , natoms, Quat4dOnes*-1. );
}
```

## MMFFBuilder::toMMFFsp3_loc() 

This is function we use in C++ to assing forcefield parameters to MMFFsp3_loc (like bond stiffnes,angle stiffnes, electron pairs, etc...)

```C++
void toMMFFsp3_loc( MMFFsp3_loc& ff, bool bRealloc=true, bool bEPairs=true, bool bUFF=false ){

        //double c0s[3]{-0.33333,-0.5,-1.0}; // cos(angle)   sp1 sp2 sp3
        double ang0s[3]{ 109.5 *M_PI/180.0, 120.0*M_PI/180.0, 180.0*M_PI/180.0 }; // cos(angle)   sp1 sp2 sp3

        int nAmax = atoms.size();
        int nCmax = confs.size();
        int npi,ne; ne=countPiE( npi, 0,nCmax );
        if(!bEPairs) ne=0;
        int nconf = nCmax;
        int ncap  = nAmax - nconf;
        int nb    = bonds.size();
        if(verbosity>0)printf(  "MM::Builder::toMMFFsp3_loc() nconf %i ncap %i npi %i ne %i \n", nconf, ncap, npi, ne  );
        int ntors=0; if(ff.bTorsion) ntors=dihedrals.size();
        if(bRealloc)ff.realloc( nconf, ncap+ne, ntors );
        Vec3d hs[4];
        int ie0=nconf+ncap;
        int iie = 0;

        if(ff.bTorsion){
            for(int i=0; i<ff.ntors; i++){
                const Dihedral& dih=dihedrals[i];
                ff.tors2atom [i]=dih.atoms;
                //ff.torsParams[i]=Quat4d{cos(dih.a0),sin(dih.a0),dih.k,dih.n};
                ff.torsParams[i]=Quat4d{cos(dih.a0),sin(dih.a0),dih.k, (double)dih.n}; 
            }
        }

        //params->printAtomTypeDict();
        int etyp=-1;  if(params) etyp=params->atomTypeDict["E"];
        for(int i=0; i<ff.natoms; i++){ ff.neighs[i]=Quat4i{-1,-1,-1,-1}; };
        for(int i=0; i<ff.nnode;  i++){ ff.bLs[i]=Quat4dZero, ff.bKs[i]=Quat4dZero, ff.Ksp[i]=Quat4dZero, ff.Kpp[i]=Quat4dZero; }; // back neighbors
        for(int ia=0; ia<nAmax; ia++ ){
            const Atom& A =  atoms[ia];
            ff.apos  [ia] = A.pos;
            ff.atypes[ia] = A.type;
            AtomType& atyp = params->atypes[A.type];

            if(A.iconf>=0){
                // Prepare params and orientation
                AtomConf& conf = confs[A.iconf];
                //printf( "Builder::toMMFFsp3_loc() [%i] ", ia ); conf.print(); printf("\n");
                int npi_neigh = countAtomPiNeighs(ia);
                //assignSp3Params( A.type, conf.nbond, conf.npi, conf.ne, npi_neigh, ff.NeighParams[ia] );

                if( conf.npi>2 ){ printf("ERROR in MM::Builder::toMMFFsp3_loc(): atom[%i].conf.npi(%i)>2 => exit() \n", ia, conf.npi); printAtomConf(ia); exit(0); }
                //double ang0 = ang0s[conf.npi];
                double ang0   = atyp.Ass*deg2rad;
                ang0 *= 0.5;
                ff.apars[ia].x = cos(ang0);    // ssC0    // cos(angle) for angles (sigma-siamg)
                ff.apars[ia].y = sin(ang0);
                //ff.apars[ia].z = 1.0;          // ssK     // stiffness  for angles
                //ff.apars[ia].w = 0;            // piC0    // angle0 for orthogonalization sigma-pi 

                ff.apars[ia].z = atyp.Kss*4.0;   // ssK     // stiffness  for angles    ... ToDo: check if K/4 or K*4
                ff.apars[ia].w = sin(atyp.Asp*deg2rad);  // piC0    // angle0 for orthogonalization sigma-pi 

                // setup ff neighbors                
                int*     ngs  = ff.neighs[ia].array;
                double*  bL   = ff.bLs[ia].array;
                double*  bK   = ff.bKs[ia].array;
                double*  Ksp  = ff.Ksp[ia].array;
                double*  Kpp  = ff.Kpp[ia].array;

                // --- Generate Bonds
                //printf( "ia[%i nbond=%i \n", ia, conf.nbond  );
                for(int k=0; k<conf.nbond; k++){
                    int ib = conf.neighs[k];
                    const Bond& B = bonds[ib];
                    
                    int ja = B.getNeighborAtom(ia);

                    const Atom& Aj =  atoms[ja];
                    AtomType& jtyp = params->atypes[Aj.type];
                    hs[k]  = atoms[ja].pos - A.pos;
                    hs[k].normalize();
                    ngs[k] = ja;
                    bL [k]=B.l0;
                    bK [k]=B.k;
                    if(bUFF){
                        Vec2d bLK = assignBondParamsUFF( ib );
                        bL [k]=bLK.x;
                        bK [k]=bLK.y;
                    }
                    //Ksp[k]=0;
                    if( (conf.npi>0)||(conf.ne>0) ){ Ksp[k]= atyp.Ksp;}else{ Ksp[k]=0; }
                    int nej  = getAtom_ne (ja);
                    int npij = getAtom_npi(ja);
                    Kpp[k]   = sqrt( atyp.Kpp * jtyp.Kpp );
                }

                makeConfGeom( conf.nbond, conf.npi, hs );

                if(bEPairs){ // --- Generate electron pairs
                    int ns = conf.nbond+conf.ne;
                    for(int k=conf.nbond; k<ns; k++){    
                        int ie=ie0+iie;
                        ngs [k] =ie;
                        //printf( "atom[%i|%i] ie %i \n", ia, k, ie );
                        ff.apos  [ie] = atoms[ia].pos + hs[k]*Lepair;
                        ff.atypes[ie] = etyp;
                        bK [k]=Kepair;
                        bL [k]=Lepair;
                        //Ksp[k]=0;                angs[iang].x = cos(ang0);
                        if( conf.npi>0 ){ Ksp[k]=atyp.Ksp; }else{ Ksp[k]=0; }  // only electron on atoms without pi-orbital are conjugted with pi-orbitas on neighboring atoms
                        iie++;
                    }
                }
                ff.pipos[ia] = hs[N_NEIGH_MAX-1]; // Pi orientation

            } // if(A.iconf>=0){
        }

        ff.bPBC = bPBC;
        ff.makeBackNeighs();
        //if( bPBC ){ ff.initPBC(); updatePBC( ff.pbcShifts ); }
        //if(verbosity>0)
        printf(  "MM::Builder::toMMFFsp3_loc() DONE \n"  );
    }
```

## MolWorld_sp3_multi (multi system buffers)

There are some relevant function used to broadcast data of one system to all the replicas. We have separate pack_system() and upload_system() in C++ for fine-control and performance. But in our pyOpenCL we target more simplicity, therefore we can combine it.

```C++

void realloc( int nSystems_ ){
    printf("MolWorld_sp3_multi::realloc() \n");
    nSystems=nSystems_;
    printf( "MolWorld_sp3_multi::realloc() Systems %i nAtoms %i nnode %i \n", nSystems, ffl.natoms,  ffl.nnode );
    ocl.initAtomsForces( nSystems, ffl.natoms,  ffl.nnode, npbc+1 );
    //printf( "MolWorld_sp3_multi::realloc() Systems %i nAtoms %i nnode %i nvecs %i \n", nSystems, ocl.nAtoms, ocl.nnode, ocl.nvecs );
    // --- dynamical
    _realloc ( atoms,     ocl.nvecs*nSystems  );
    _realloc0( aforces,   ocl.nvecs*nSystems  , Quat4fZero );
    _realloc0( avel,      ocl.nvecs*nSystems  , Quat4fZero );
    _realloc0( cvfs,      ocl.nvecs*nSystems  , Quat4fZero );
    _realloc0( constr,    ocl.nAtoms*nSystems , Quat4fOnes*-1. );
    _realloc0( constrK,   ocl.nAtoms*nSystems , Quat4fOnes*-1. );
    _realloc0( bboxes,   nSystems, cl_Mat3{cl_float4{-1e+8,-1e+8,-1e+8,-1e+8,},cl_float4{+1e+8,+1e+8,+1e+8,+1e+8,}, cl_float4{-1.,-1.,-1.,-1.} }   );
    // --- params
    _realloc( neighs,    ocl.nAtoms*nSystems );
    _realloc( neighCell, ocl.nAtoms*nSystems );
    _realloc( bkNeighs,    ocl.nvecs*nSystems );
    _realloc( bkNeighs_new,ocl.nvecs*nSystems );
    _realloc( REQs,      ocl.nAtoms*nSystems );
    _realloc( MMpars,    ocl.nnode*nSystems  );
    _realloc( BLs,       ocl.nnode*nSystems  );
    _realloc( BKs,       ocl.nnode*nSystems  );
    _realloc( Ksp,       ocl.nnode*nSystems  );
    _realloc( Kpp,       ocl.nnode*nSystems  );

    _realloc( lvecs,     nSystems  );
    _realloc( ilvecs,    nSystems  );
    _realloc( MDpars,    nSystems  );
    _realloc0( TDrive,   nSystems, Quat4f{0.0,-1.0,0.0,0.0} );
    _realloc( pbcshifts, ocl.npbc*nSystems );

    _realloc( fire,      nSystems  );

    // ToDo : it may be good to bind buffer directly in p_cpu buffer inside   OCLsystem::newBuffer()
}

void pack_system( int isys, MMFFsp3_loc& ff, bool bParams=0, bool bForces=false, bool bVel=false, bool blvec=true, float l_rnd=-1 ){

    int i0n   = isys * ocl.nnode;
    int i0a   = isys * ocl.nAtoms;
    int i0v   = isys * ocl.nvecs;
    int i0pbc = isys*ocl.npbc;

    if(blvec){
        if(npbc==0){ pbcshifts[isys].f=Vec3fZero; };
        pack( npbc,  pbc_shifts, pbcshifts+i0pbc );

        Mat3_to_cl( ff.   lvec,  lvecs[isys] );
        Mat3_to_cl( ff.invLvec, ilvecs[isys] );
    }
    
    pack( ff.nvecs, ff.apos, atoms+i0v );
    for(int i=0; i<ocl.nAtoms; i++){ 
        Quat4f a=atoms[i+i0v]; 
        //a.w=-1.0;
        a.w = ffl.constr[i].w; 
        constr [i+i0a] = a; 
        constrK[i+i0a] = (Quat4f)ffl.constrK[i]; 
        //printf(  "atom[%3i|sys=%i](%8.5f,%8.5f,%8.5f) constr(%8.5f,%8.5f,%8.5f|%8.5f)\n", i, isys, ff.apos[i].x,ff.apos[i].y,ff.apos[i].z,     a.x,a.y,a.z,a.w );
    }  // contrains
    //Mat3d lvec0=ff.lvec;  
    if(bForces){ pack( ff.nvecs, ff.fapos, aforces+i0v ); }
    //if(bVel   ){ pack( ff.nvecs, opt.vel,  avel   +i0v ); }
    if(bParams){
        copy    ( ff.natoms, ff.neighCell, neighCell+i0a );
        copy    ( ff.natoms, ff.neighs,    neighs   +i0a );
        //copy_add( ff.natoms, ff.neighs,    neighs   +i0a,           0      );
        copy    ( ff.natoms, ff.bkneighs,  bkNeighs_new +i0v           );
        copy    ( ff.nnode,  ff.bkneighs,  bkNeighs_new +i0v+ff.natoms );
        copy_add( ff.natoms, ff.bkneighs,   bkNeighs +i0v,           i0n*8  );
        copy_add( ff.nnode , ff.bkneighs,   bkNeighs +i0v+ff.natoms, i0n*8 + 4*ff.nnode );
        pack    ( ff.nnode , ff.apars,      MMpars   +i0n );
        pack    ( ff.nnode , ff.bLs,        BLs      +i0n );
        pack    ( ff.nnode , ff.bKs,        BKs      +i0n );
        pack    ( ff.nnode , ff.Ksp,        Ksp      +i0n );
        pack    ( ff.nnode , ff.Kpp,        Kpp      +i0n );
        pack    ( nbmol.natoms, nbmol.REQs, REQs     +i0a );
    }
}

void upload_sys( int isys, bool bParams=false, bool bForces=0, bool bVel=true, bool blvec=true ){
    //printf("MolWorld_sp3_multi::upload() \n");
    int i0v   = isys * ocl.nvecs;
    int i0a   = isys * ocl.nAtoms;
    int err=0;
    err|= ocl.upload( ocl.ibuff_atoms,  atoms,    ocl.nvecs,  i0v );
    err|= ocl.upload( ocl.ibuff_constr,  constr,  ocl.nAtoms, i0a );
    err|= ocl.upload( ocl.ibuff_constrK, constrK, ocl.nAtoms, i0a );
    err|= ocl.upload( ocl.ibuff_bboxes,  bboxes, 1, isys  );
    if(bForces){ err|= ocl.upload( ocl.ibuff_aforces, aforces, ocl.nvecs, i0v ); }
    if(bVel   ){ 
        err|= ocl.upload( ocl.ibuff_avel,    avel,    ocl.nvecs, i0v ); 
        err|= ocl.upload( ocl.ibuff_cvf,     cvfs,    ocl.nvecs, i0v ); 
    }
    if(blvec){
        int i0pbc = isys * ocl.npbc;
        err|= ocl.upload( ocl.ibuff_lvecs,     lvecs     , 1, isys  );
        err|= ocl.upload( ocl.ibuff_ilvecs,    ilvecs    , 1, isys  );
        err|= ocl.upload( ocl.ibuff_pbcshifts, pbcshifts , ocl.npbc, i0pbc );
    }
    if(bParams){
        int i0n   = isys * ocl.nnode;
        int i0bk  = isys * ocl.nbkng;
        err|= ocl.upload( ocl.ibuff_neighs,    neighs    , ocl.nAtoms, i0a  );
        err|= ocl.upload( ocl.ibuff_neighCell, neighCell , ocl.nAtoms, i0a  );
        err|= ocl.upload( ocl.ibuff_bkNeighs,  bkNeighs  , ocl.nbkng, i0bk  );
        err|= ocl.upload( ocl.ibuff_bkNeighs_new, bkNeighs_new, ocl.nbkng, i0bk  );
        //err|= ocl.upload( ocl.ibuff_neighForce, neighForce  );
        err|= ocl.upload( ocl.ibuff_REQs,   REQs,  ocl.nAtoms, i0a  );
        err|= ocl.upload( ocl.ibuff_MMpars, MMpars, ocl.nnode, i0n );
        err|= ocl.upload( ocl.ibuff_BLs, BLs, ocl.nnode, i0n );
        err|= ocl.upload( ocl.ibuff_BKs, BKs, ocl.nnode, i0n );
        err|= ocl.upload( ocl.ibuff_Ksp, Ksp, ocl.nnode, i0n );
        err|= ocl.upload( ocl.ibuff_Kpp, Kpp, ocl.nnode, i0n );
    }
    err |= ocl.finishRaw(); 
    OCL_checkError(err, "MolWorld_sp2_multi::upload().finish");
}
```

## MolWorld_sp3_multi::run_ocl_opt()

This is the main driver function which runs the optimization loop on GPU. Gradually executing each krenel. We run muptiple iterations in innter loop (nPerVFs-times) before we download the results and check the convergence. 
This is performance optimization. 

```C++
int run_ocl_opt( int niter, double Fconv=1e-6 ){ 
    double F2conv = Fconv*Fconv;
    picked2GPU( ipicked,  1.0 );
    int err=0;
    if( task_MMFF==0)setup_MMFFf4_ocl();

    //int nPerVFs = _min(1,niter);
    nPerVFs = _min(10,niter);
    //int nPerVFs = _min(50,niter);
    int nVFs    = niter/nPerVFs;
    long T0     = getCPUticks();
    int niterdone=0;
    double F2=0;

    bool dovdW=true;
    //bool dovdW=false;
    ocl.bSubtractVdW=dovdW;
    
    if(ocl.nGroupTot<=0){ bGroups = false; };
    bool bExplore = false;

    for(int i=0; i<nVFs; i++){
        if(bAnimManipulation){ animate(); }

        bExplore = updateMultiExploring( Fconv );

        bool bGroupDrive = bGroups && bExplore;
        for(int j=0; j<nPerVFs; j++){
            {    
                if( bGroupDrive )err |= task_GroupUpdate->enque_raw();
                if(dovdW)[[likely]]{
                    if(bSurfAtoms)[[likely]]{
                        if  (bGridFF)[[likely]]{ 
                            err |= task_NBFF_Grid ->enque_raw();   //OCL_checkError(err, "task_NBFF_Grid->enque_raw(); ");
                        }else { 
                            //printf( "task_NBFF(), task_SurfAtoms() \n" );
                            err |= task_NBFF     ->enque_raw();  //OCL_checkError(err, "MolWorld_sp3_multi::run_ocl_opt().task_NBFF()" ); 
                            err |= task_SurfAtoms->enque_raw();  //OCL_checkError(err, "MolWorld_sp3_multi::run_ocl_opt().task_SurfAtoms()" );
                        }
                    }else{ 
                        err |= task_NBFF      ->enque_raw();     //OCL_checkError(err, "task_NBFF->enque_raw();");
                    }
                }
                err |= task_MMFF->enque_raw();    //OCL_checkError(err, "task_MMFF->enque_raw()");

                if( bGroupDrive ) err |= task_GroupForce->enque_raw();
                err |= task_move->enque_raw();    //OCL_checkError(err, "task_move->enque_raw()");
            }
            niterdone++;
            nloop++;
        }
        
        F2 = evalVFs( Fconv );
        if( F2<F2conv  ){ 
            double t=(getCPUticks()-T0)*tick2second;
            //printf( "run_omp_ocl(nSys=%i|iPara=%i) CONVERGED in %i/%i nsteps |F|=%g time=%g[ms]\n", nSystems, iParalel, itr,niter_max, sqrt(F2max), T1*1000 );
            if(verbosity>0)
            printf( "run_ocl_opt(nSys=%i|iPara=%i,bSurfAtoms=%i,bGridFF=%i,bExplore=%i) CONVERGED in %i/%i steps, |F|(%g)<%g time %g [ms]( %g [us/step]) bGridFF=%i \n", nSystems, iParalel, bSurfAtoms, bGridFF, bExplore, niterdone,niter, sqrt(F2), Fconv, t*1000, t*1e+6/niterdone, bGridFF ); 
            return niterdone; 
        }
    }
    if(bGroups){
        err|= ocl.download( ocl.ibuff_gcenters, gcenters );
        err|= ocl.download( ocl.ibuff_gfws, gfws );
        err|= ocl.download( ocl.ibuff_gups, gups );
    }
    err|= ocl.download( ocl.ibuff_atoms,    atoms   );
    err|= ocl.download( ocl.ibuff_aforces,  aforces );
    err|= ocl.finishRaw();  
    OCL_checkError(err, "run_ocl_opt().finishRaw()");
    checkBordersOfBbox();
    double t=(getCPUticks()-T0)*tick2second;
    if(verbosity>0)printf( "run_ocl_opt(nSys=%i|iPara=%i,bSurfAtoms=%i,bGridFF=%i,bExplore=%i) NOT CONVERGED in %i steps, |F|(%g)>%g time %g [ms]( %g [us/step]) bGridFF=%i iSysFMax=%i dovdW=%i \n", nSystems, iParalel, bSurfAtoms, bGridFF, bExplore,  niter, sqrt(F2), Fconv, t*1000, t*1e+6/niterdone, bGridFF, iSysFMax, dovdW ); 
    return niterdone;
}
```

it uses OpenCL tansks created here:

```C++
void setup_MMFFf4_ocl(){
    printf("MolWorld_sp3_multi::setup_MMFFf4_ocl() \n");
    ocl.nDOFs.x=ff.natoms;
    ocl.nDOFs.y=ff.nnode;
    if(!task_cleanF)   task_cleanF = ocl.setup_cleanForceMMFFf4 ( ffl.natoms, ffl.nnode       );
    if(!task_move  )   task_move   = ocl.setup_updateAtomsMMFFf4( ffl.natoms, ffl.nnode       );
    if(!task_MMFF  )   task_MMFF   = ocl.setup_getMMFFf4        ( ffl.natoms, ffl.nnode, bPBC );
    Vec3i nPBC_=nPBC; if(!bPBC){ nPBC_=Vec3iZero; }; 
    if((!task_NBFF_Grid)&&bGridFF ){ task_NBFF_Grid = ocl.setup_getNonBond_GridFF( ffl.natoms, ffl.nnode, nPBC_ ); } 
    if(!task_NBFF                 ){ task_NBFF      = ocl.setup_getNonBond       ( ffl.natoms, ffl.nnode, nPBC_ ); }
    if(!task_cleanF)task_cleanF = ocl.setup_cleanForceMMFFf4 ( ffl.natoms, ffl.nnode       );
}
```

The implementation is here

```C++

    OCLtask* setup_getNonBond( int na, int nNode, Vec3i nPBC_, OCLtask* task=0){
        printf("setup_getNonBond(na=%i,nnode=%i) \n", na, nNode);
        if(task==0) task = getTask("getNonBond");
        //int nloc = 1;
        //int nloc = 2;
        //int nloc = 4;
        //int nloc = 8;
        //int nloc = 16;
        int nloc = 32;
        //int nloc = 64;
        task->local.x  = nloc;
        task->global.x = na + nloc-(na%nloc); // round up to multiple of nloc
        task->global.y = nSystems;

        useKernel( task->ikernel );
        nDOFs.x=na; 
        nDOFs.y=nNode;
        //nDOFs.x=bPBC; 
        v2i4( nPBC_, nPBC );
        // ------- Maybe We do-not need to do this every frame ?
        int err=0;
        err |= _useArg   ( nDOFs );               // 1 
        // Dynamical
        err |= useArgBuff( ibuff_atoms      ); // 2  
        err |= useArgBuff( ibuff_aforces    ); // 3
        // parameters
        err |= useArgBuff( ibuff_REQs      );  // 4
        err |= useArgBuff( ibuff_neighs    );  // 5
        err |= useArgBuff( ibuff_neighCell );  // 6
        err |= useArgBuff( ibuff_lvecs     );  // 7
        //err |= _useArg( cl_lvec          );  // 7
        err |= _useArg( nPBC               );  // 8
        err |= _useArg( GFFparams          );  // 9
        OCL_checkError(err, "setup_getNonBond");
        return task;
        // const int4 ns,                  // 1
        // // Dynamical
        // __global float4*  atoms,        // 2
        // __global float4*  forces,       // 3
        // // Parameters
        // __global float4*  REQKs,        // 4
        // __global int4*    neighs,       // 5
        // __global int4*    neighCell,    // 6
        // const int4 nPBC,                // 7
        // const cl_Mat3 lvec,             // 8
        // float R2damp                    // 9
    }

    OCLtask* setup_getNonBond_GridFF( int na, int nNode, Vec3i nPBC_, OCLtask* task=0){
        printf("setup_getNonBond_GridFF(na=%i,nnode=%i) itex_FE_Paul=%i itex_FE_Lond=%i itex_FE_Coul=%i\n", na, nNode,  itex_FE_Paul,itex_FE_Lond,itex_FE_Coul );
        if((itex_FE_Paul<0)||(itex_FE_Paul<=0)||(itex_FE_Paul<=0)){ printf( "ERROR in setup_getNonBond_GridFF() GridFF textures not initialized(itex_FE_Paul=%i itex_FE_Lond=%i itex_FE_Coul=%i) => Exit() \n", itex_FE_Paul,itex_FE_Lond,itex_FE_Coul ); exit(0); }
        if(task==0) task = getTask("getNonBond_GridFF");        
        //int nloc = 1;
        //int nloc = 4;
        //int nloc = 8;
        int nloc = 32;
        //int nloc = 64;
        task->local.x  = nloc;
        task->global.x = na + nloc-(na%nloc);
        task->global.y = nSystems;
        grid_shift0_p0 = grid_p0 + grid_shift0;

        useKernel( task->ikernel );
        nDOFs.x=na; 
        nDOFs.y=nNode; 

        v2i4( nPBC_, nPBC );
        // ------- Maybe We do-not need to do this every frame ?
        int err=0;
        err |= _useArg   ( nDOFs );           // 1
        // Dynamical
        err |= useArgBuff( ibuff_atoms     ); // 2
        err |= useArgBuff( ibuff_aforces   ); // 3
        // parameters
        err |= useArgBuff( ibuff_REQs      );  // 4
        err |= useArgBuff( ibuff_neighs    );  // 5
        err |= useArgBuff( ibuff_neighCell );  // 6
        err |= useArgBuff( ibuff_lvecs     );  // 7
        err |= _useArg( nPBC               );  // 8
        err |= _useArg( GFFparams          );  // 9
        err |= useArgBuff( itex_FE_Paul    );  // 10
        err |= useArgBuff( itex_FE_Lond    );  // 11
        err |= useArgBuff( itex_FE_Coul    );  // 12   
        //err |= _useArg( cl_diGrid          );  // 13
        err |= _useArg( cl_grid_ilvec      );  // 13
        err |= _useArg( grid_shift0_p0     );  // 14
        OCL_checkError(err, "setup_getNonBond_GridFF");
        return task;
        // const int4 ns,                  // 1
        // // Dynamical
        // __global float4*  atoms,        // 2
        // __global float4*  forces,       // 3
        // // Parameters
        // __global float4*  REQKs,        // 4
        // __global int4*    neighs,       // 5
        // __global int4*    neighCell,    // 6
        // __global cl_Mat3* lvecs,        // 7
        // const int4 nPBC,                // 8
        // const float4 GFFparams,         // 9
        // // GridFF
        // __read_only image3d_t  FE_Paul, // 10
        // __read_only image3d_t  FE_Lond, // 11
        // __read_only image3d_t  FE_Coul, // 12
        // const cl_Mat3  grid_invd,       // 13
        // const float4   grid_p0          // 14
    }

    OCLtask* setup_getMMFFf4( int na, int nNode, bool bPBC=false, OCLtask* task=0){
        printf("setup_getMMFFf4(na=%i,nnode=%i) \n", na, nNode);
        if(task==0) task = getTask("getMMFFf4");
        int nloc = 32;
        task->local.x  = nloc;
        task->global.x = nNode + nloc-(nNode%nloc); // round up to multiple of nloc
        //task->global.x = nNode;
        task->global.y = nSystems;
        useKernel( task->ikernel );
        nDOFs.x=na; 
        nDOFs.y=nNode; 
        nDOFs.w=bPBC; 
        // ------- Maybe We do-not need to do this every frame ?
        int err=0;
        err |= _useArg   ( nDOFs );            // 1
        // Dynamical
        err |= useArgBuff( ibuff_atoms  );     // 2
        err |= useArgBuff( ibuff_aforces);     // 3
        err |= useArgBuff( ibuff_neighForce ); // 4
        // parameters
        err |= useArgBuff( ibuff_neighs    );  // 5
        err |= useArgBuff( ibuff_neighCell );  // 5
        err |= useArgBuff( ibuff_REQs   );     // 6
        err |= useArgBuff( ibuff_MMpars );     // 7
        err |= useArgBuff( ibuff_BLs    );     // 8
        err |= useArgBuff( ibuff_BKs    );     // 9
        err |= useArgBuff( ibuff_Ksp    );     // 10
        err |= useArgBuff( ibuff_Kpp    );     // 11
        err |= useArgBuff( ibuff_lvecs  );     // 12
        err |= useArgBuff( ibuff_ilvecs );     // 13
        err |= useArgBuff( ibuff_pbcshifts );  // 13
        err |= _useArg   ( npbc         );  
        err |= _useArg   ( bSubtractVdW ); 
        

        //err |= _useArg( cl_lvec    );        // 12
        //err |= _useArg( cl_invLvec );        // 13
        OCL_checkError(err, "setup_getMMFFf4");
        return task;
        // const int4 nDOFs,              // 1   (nAtoms,nnode)
        // // Dynamical
        // __global float4*  apos,        // 2    [natoms]
        // __global float4*  fapos,       // 3    [natoms]     
        // __global float4*  fneigh,      // 4    [nnode*4]
        // // parameters
        // __global int4*    neighs,       // 5  [nnode]  neighboring atoms
        // __global float4*  REQKs,        // 6  [natoms] non-boding parametes {R0,E0,Q} 
        // __global float4*  apars,        // 7 [nnode]  per atom forcefield parametrs {c0ss,Kss,c0sp}
        // __global float4*  bLs,          // 8 [nnode]  bond lengths  for each neighbor
        // __global float4*  bKs,          // 9 [nnode]  bond stiffness for each neighbor
        // __global float4*  Ksp,          // 10 [nnode]  stiffness of pi-alignment for each neighbor
        // __global float4*  Kpp,          // 11 [nnode]  stiffness of pi-planarization for each neighbor
        // const cl_Mat3 lvec,             // 12
        // const cl_Mat3 invLvec           // 13
    }

    OCLtask* setup_updateAtomsMMFFf4( int na, int nNode,  OCLtask* task=0 ){
        printf( "setup_updateAtomsMMFFf4() \n" );
        if(task==0) task = getTask("updateAtomsMMFFf4");
        int nvec = na+nNode;
        //task->local .x = 1;
        int nloc=32;
        task->local.x   = nloc;
        task->global.x  = nvec + nloc-(nvec%nloc);
        //task->global.x = nvec;
        task->global.y = nSystems;
        printf( "setup_updateAtomsMMFFf4() global.x=%i local.x=%i glob/loc=%g \n",  task->local.x, task->global.x, (task->global.x)/((float)task->local.x)  );
        //task->roundSizes();
        //if(n >=0  ) 
        nDOFs.x=na; 
        nDOFs.y=nNode; 
        //nDOFs.w=nMaxSysNeighs;
        nDOFs.w=0;
        int err=0;
        useKernel( task->ikernel  );
        err |= _useArg( nDOFs     );           // 1
        err |= useArgBuff( ibuff_atoms      ); // 2
        err |= useArgBuff( ibuff_avel       ); // 3
        err |= useArgBuff( ibuff_aforces    ); // 4
        err |= useArgBuff( ibuff_cvf        ); // 5
        err |= useArgBuff( ibuff_neighForce ); // 6
        err |= useArgBuff( ibuff_bkNeighs   ); // 7
        err |= useArgBuff( ibuff_constr     ); // 8
        err |= useArgBuff( ibuff_constrK    ); // 9
        err |= useArgBuff( ibuff_MDpars     ); // 10
        err |= useArgBuff( ibuff_TDrive     ); // 11
        err |= useArgBuff( ibuff_bboxes     ); // 12
        err |= useArgBuff( ibuff_sysneighs  ); // 13
        err |= useArgBuff( ibuff_sysbonds   ); // 14
        OCL_checkError(err, "setup_updateAtomsMMFFf4");
        return task;
        // const int4        n,            // 1 // (natoms,nnode) dimensions of the system
        // __global float4*  apos,         // 2 // positions of atoms  (including node atoms [0:nnode] and capping atoms [nnode:natoms] and pi-orbitals [natoms:natoms+nnode] )
        // __global float4*  avel,         // 3 // velocities of atoms 
        // __global float4*  aforce,       // 4 // forces on atoms
        // __global float4*  cvf,          // 5 // damping coefficients for velocity and force
        // __global float4*  fneigh,       // 6 // recoil forces on neighbors (and pi-orbitals)
        // __global int4*    bkNeighs,     // 7 // back neighbors indices (for recoil forces)
        // __global float4*  constr,       // 8 // constraints (x,y,z,K) for each atom
        // __global float4*  constrKs,     // 9 // constraints stiffness (kx,ky,kz,?) for each atom
        // __global float4*  MDparams,     // 10 // MD parameters (dt,damp,Flimit)
        // __global float4*  TDrives,       // 11 // Thermal driving (T,gamma_damp,seed,?)
        // __global cl_Mat3* bboxes        // 12 // bounding box (xmin,ymin,zmin)(xmax,ymax,zmax)(kx,ky,kz)
    }
```

## C++ interface to allocation of OpenCL buffers

```C++
    int initAtomsForces( int nSystems_, int nAtoms_, int nnode_, int npbc_ ){
        nSystems=nSystems_;
        nnode  = nnode_;
        nAtoms = nAtoms_;
        npi    = nnode_;
        nvecs  = nAtoms+npi;   // number of vectors (atoms and pi-orbitals)
        nbkng  = nnode*4*2;    // number of back-neighbors (4 neighbors per node, position and pi-orbital)
        ncap   = nAtoms-nnode; // number of capping atoms
        if(npbc_==0){ npbc=1; }; 
        npbc   = npbc_;
        printf( "initAtomsForces() nSystems %i nvecs %i natoms %i nnode %i nbkng %i \n", nSystems, nvecs, nAtoms, nnode, nbkng );
        printf( "initAtomsForces() nS*nvecs %i nS*natoms %i nS*nnode %i nS*nbkng %i \n", nSystems*nvecs,  nSystems*nAtoms, nSystems*nnode, nSystems*nbkng );
        if( (nSystems<=0)||(nAtoms<=0) ){ printf("ERROR in OCL_MM::initAtomsForces() invalid size nSystems=%i nAtoms=%i => Exit() \n", nSystems, nAtoms); exit(0); }
        ibuff_atoms      = newBuffer( "atoms",      nSystems*nvecs , sizeof(float4), 0, CL_MEM_READ_WRITE ); // atoms positions and velocities (x,y,z,m)
        ibuff_aforces    = newBuffer( "aforces",    nSystems*nvecs , sizeof(float4), 0, CL_MEM_READ_WRITE ); // atoms forces
        ibuff_REQs       = newBuffer( "REQs",       nSystems*nAtoms, sizeof(float4), 0, CL_MEM_READ_ONLY  ); // atoms parameters {R0,E0,Q}
        ibuff_neighs     = newBuffer( "neighs",     nSystems*nAtoms, sizeof(int4  ), 0, CL_MEM_READ_ONLY  );
        ibuff_neighCell  = newBuffer( "neighCell" , nSystems*nAtoms, sizeof(int4  ), 0, CL_MEM_READ_ONLY  );

        ibuff_constr     = newBuffer( "constr",     nSystems*nAtoms , sizeof(float4), 0, CL_MEM_READ_WRITE );
        ibuff_constrK    = newBuffer( "constrK",    nSystems*nAtoms , sizeof(float4), 0, CL_MEM_READ_WRITE );
        //ibuff_constr0    = newBuffer( "constr0",   nSystems*nAtoms , sizeof(float4), 0, CL_MEM_READ_WRITE );
        //ibuff_constrK    = newBuffer( "constrK",   nSystems*nAtoms , sizeof(float4), 0, CL_MEM_READ_WRITE );
        ibuff_bboxes     = newBuffer( "bboxes",      nSystems,        sizeof(cl_Mat3), 0, CL_MEM_READ_ONLY  );

        ibuff_bkNeighs     = newBuffer( "bkNeighs", nSystems*nvecs,     sizeof(int4  ), 0, CL_MEM_READ_ONLY  );     // back neighbors
        ibuff_bkNeighs_new = newBuffer( "bkNeighs_new", nSystems*nvecs, sizeof(int4  ), 0, CL_MEM_READ_ONLY  );   
        ibuff_avel       = newBuffer( "avel",       nSystems*nvecs,     sizeof(float4), 0, CL_MEM_READ_WRITE );     // atoms velocities (x,y,z,m)
        ibuff_cvf        = newBuffer( "cvf",        nSystems*nvecs ,    sizeof(float4), 0, CL_MEM_READ_WRITE );
        ibuff_neighForce = newBuffer( "neighForce", nSystems*nbkng,     sizeof(float4), 0, CL_MEM_READ_WRITE );

        ibuff_MMpars     = newBuffer( "MMpars",     nSystems*nnode,  sizeof(int4),   0, CL_MEM_READ_ONLY  );
        ibuff_BLs        = newBuffer( "BLs",        nSystems*nnode,  sizeof(float4), 0, CL_MEM_READ_ONLY  );
        ibuff_BKs        = newBuffer( "BKs",        nSystems*nnode,  sizeof(float4), 0, CL_MEM_READ_ONLY  );
        ibuff_Ksp        = newBuffer( "Ksp",        nSystems*nnode,  sizeof(float4), 0, CL_MEM_READ_ONLY  );
        ibuff_Kpp        = newBuffer( "Kpp",        nSystems*nnode,  sizeof(float4), 0, CL_MEM_READ_ONLY  );

        ibuff_MDpars     = newBuffer( "MDpars",     nSystems,        sizeof(float4),  0, CL_MEM_READ_ONLY  );
        ibuff_TDrive     = newBuffer( "TDrive",     nSystems,        sizeof(float4),  0, CL_MEM_READ_ONLY  );
        ibuff_lvecs      = newBuffer( "lvecs",      nSystems,        sizeof(cl_Mat3), 0, CL_MEM_READ_ONLY  );
        ibuff_ilvecs     = newBuffer( "ilvecs",     nSystems,        sizeof(cl_Mat3), 0, CL_MEM_READ_ONLY  );

        ibuff_pbcshifts  = newBuffer( "pbcshifts",  nSystems*npbc,   sizeof(float4), 0, CL_MEM_READ_ONLY  );

        ibuff_sysneighs  = newBuffer( "sysneighs",  nSystems*nMaxSysNeighs,  sizeof(int),    0, CL_MEM_READ_ONLY  );   // for each system contains array int[nMaxSysNeighs] of nearby other systems
        ibuff_sysbonds   = newBuffer( "sysbonds",   nSystems*nMaxSysNeighs,  sizeof(float4), 0, CL_MEM_READ_ONLY  );   // contains parameters of bonds (constrains) with neighbor systems   {Lmin,Lmax,Kpres,Ktens}

        // int nsamp_max = 1000; DEBUG
        // ibuff_samp_fs   = newBuffer( "samp_fs",   nsamp_max, sizeof(float4), 0, CL_MEM_READ_WRITE );   DEBUG
        // ibuff_samp_ps   = newBuffer( "samp_ps",   nsamp_max, sizeof(float4), 0, CL_MEM_READ_ONLY  );    DEBUG
        // ibuff_samp_REQ  = newBuffer( "samp_REQ",  nsamp_max, sizeof(float4), 0, CL_MEM_READ_ONLY  );    DEBUG

        return ibuff_atoms;
    }
```

## OpenCL kernel headers

```
// ======================================================================
//                          getMMFFf4()
// ======================================================================
// 1.  getMMFFf4() - computes bonding interactions between atoms and nodes and its neighbors (max. 4 neighbors allowed), the resulting forces on atoms are stored "fapos" array and recoil forces on neighbors are stored in "fneigh" array
//                   kernel run over all atoms and all systems in parallel to exploit GPU parallelism
//__attribute__((reqd_work_group_size(1,1,1)))
__kernel void getMMFFf4(
    const int4 nDOFs,               // 1   (nAtoms,nnode) dimensions of the system
    // Dynamical
    __global float4*  apos,         // 2  [natoms]     positions of atoms (including node atoms [0:nnode] and capping atoms [nnode:natoms] and pi-orbitals [natoms:natoms+nnode] )
    __global float4*  fapos,        // 3  [natoms]     forces on    atoms (just node atoms are evaluated)
    __global float4*  fneigh,       // 4  [nnode*4*2]  recoil forces on neighbors (and pi-orbitals)
    // parameters
    __global int4*    neighs,       // 5  [nnode]  neighboring atoms
    __global int4*    neighCell,    // 5  [nnode]  neighboring atom  cell index
    __global float4*  REQKs,        // 6  [natoms] non-boding parametes {R0,E0,Q} i.e. R0: van der Waals radii, E0: well depth and partial charge, Q: partial charge 
    __global float4*  apars,        // 7  [nnode]  per atom forcefield parametrs {c0ss,Kss,c0sp}, i.e. c0ss: cos(equlibrium angle/2) for sigma-sigma; Kss: stiffness of sigma-sigma angle; c0sp: is cos(equlibrium angle) for sigma-pi
    __global float4*  bLs,          // 8  [nnode]  bond length    between node and each neighbor
    __global float4*  bKs,          // 9  [nnode]  bond stiffness between node and each neighbor
    __global float4*  Ksp,          // 10 [nnode]  stiffness of pi-alignment for each neighbor     (only node atoms have pi-pi alignemnt interaction)
    __global float4*  Kpp,          // 11 [nnode]  stiffness of pi-planarization for each neighbor (only node atoms have pi-pi alignemnt interaction)
    __global cl_Mat3* lvecs,        // 12 lattice vectors         for each system
    __global cl_Mat3* ilvecs,       // 13 inverse lattice vectors for each system
    __global float4*  pbc_shifts,
    const int npbc,
    const int bSubtractVdW
){


// Assemble recoil forces from neighbors and  update atoms positions and velocities 
//__attribute__((reqd_work_group_size(1,1,1)))
__kernel void updateAtomsMMFFf4(
    const int4        n,            // 1 // (natoms,nnode) dimensions of the system
    __global float4*  apos,         // 2 // positions of atoms  (including node atoms [0:nnode] and capping atoms [nnode:natoms] and pi-orbitals [natoms:natoms+nnode] )
    __global float4*  avel,         // 3 // velocities of atoms 
    __global float4*  aforce,       // 4 // forces on atoms
    __global float4*  cvf,          // 5 // damping coefficients for velocity and force
    __global float4*  fneigh,       // 6 // recoil forces on neighbors (and pi-orbitals)
    __global int4*    bkNeighs,     // 7 // back neighbors indices (for recoil forces)
    __global float4*  constr,       // 8 // constraints (x,y,z,K) for each atom
    __global float4*  constrK,      // 9 // constraints stiffness (kx,ky,kz,?) for each atom
    __global float4*  MDparams,     // 10 // MD parameters (dt,damp,Flimit)
    __global float4*  TDrives,      // 11 // Thermal driving (T,gamma_damp,seed,?)
    __global cl_Mat3* bboxes,       // 12 // bounding box (xmin,ymin,zmin)(xmax,ymax,zmax)(kx,ky,kz)
    __global int*     sysneighs,    // 13 // // for each system contains array int[nMaxSysNeighs] of nearby other systems
    __global float4*  sysbonds      // 14 // // contains parameters of bonds (constrains) with neighbor systems   {Lmin,Lmax,Kpres,Ktens}
)

__kernel void printOnGPU(
    const int4        n,            // 1
    const int4        mask,         // 2
    __global float4*  apos,         // 3
    __global float4*  avel,         // 4
    __global float4*  aforce,       // 5
    __global float4*  fneigh,       // 6
    __global int4*    bkNeighs,     // 7
    __global float4*  constr        // 8
);

__kernel void cleanForceMMFFf4(
    const int4        n,           // 2
    __global float4*  aforce,      // 5
    __global float4*  fneigh       // 6
){

// ======================================================================
//                           getNonBond()
// ======================================================================
// Calculate non-bonded forces on atoms (icluding both node atoms and capping atoms), cosidering periodic boundary conditions
// It calculate Lenard-Jones, Coulomb and Hydrogen-bond forces between all atoms in the system
// it can be run in parallel for multiple systems, in order to efficiently use number of GPU cores (for small systems with <100 this is essential to get good performance)
// This is the most time consuming part of the forcefield evaluation, especially for large systems when nPBC>1
__attribute__((reqd_work_group_size(32,1,1)))
__kernel void getNonBond(
    const int4 ns,                  // 1 // (natoms,nnode) dimensions of the system
    // Dynamical
    __global float4*  atoms,        // 2 // positions of atoms  (including node atoms [0:nnode] and capping atoms [nnode:natoms] and pi-orbitals [natoms:natoms+nnode] )
    __global float4*  forces,       // 3 // forces on atoms
    // Parameters
    __global float4*  REQKs,        // 4 // non-bonded parameters (RvdW,EvdW,QvdW,Hbond)
    __global int4*    neighs,       // 5 // neighbors indices      ( to ignore interactions between bonded atoms )
    __global int4*    neighCell,    // 6 // neighbors cell indices ( to know which PBC image should be ignored  due to bond )
    __global cl_Mat3* lvecs,        // 7 // lattice vectors for each system
    const int4        nPBC,         // 8 // number of PBC images in each direction (x,y,z)
    const float4      GFFParams     // 9 // Grid-Force-Field parameters
);


// ======================================================================
//                           getNonBond_GridFF()
// ======================================================================
// Calculate non-bonded forces on atoms (icluding both node atoms and capping atoms), cosidering periodic boundary conditions
// It calculate Lenard-Jones, Coulomb and Hydrogen-bond forces between all atoms in the system, and interactions of these atoms rigid surface described by Grid-Force-Field (GFF) done by texture interpolation
// it can be run in parallel for multiple systems, in order to efficiently use number of GPU cores (for small systems with <100 this is essential to get good performance)
// This is the most time consuming part of the forcefield evaluation, especially for large systems when nPBC>1
// NOTE: this modified version of getNonBond() just by including Grid-Force-Field (GFF) forces
__attribute__((reqd_work_group_size(32,1,1)))
__kernel void getNonBond_GridFF(
    const int4 ns,                  // 1 // dimensions of the system (natoms,nnode,nvec)
    // Dynamical
    __global float4*  atoms,        // 2 // positions of atoms
    __global float4*  forces,       // 3 // forces on atoms
    // Parameters
    __global float4*  REQKs,        // 4 // parameters of Lenard-Jones potential, Coulomb and Hydrogen Bond (RvdW,EvdW,Q,H)
    __global int4*    neighs,       // 5 // indexes neighbors of atoms
    __global int4*    neighCell,    // 6 // indexes of cells of neighbor atoms
    __global cl_Mat3* lvecs,        // 7 // lattice vectors of the system
    const int4 nPBC,                // 8 // number of PBC images in each direction
    const float4  GFFParams,        // 9 // parameters of Grid-Force-Field (GFF) (RvdW,EvdW,Q,H)
    // GridFF
    __read_only image3d_t  FE_Paul, // 10 // Grid-Force-Field (GFF) for Pauli repulsion
    __read_only image3d_t  FE_Lond, // 11 // Grid-Force-Field (GFF) for London dispersion
    __read_only image3d_t  FE_Coul, // 12 // Grid-Force-Field (GFF) for Coulomb interaction
    const cl_Mat3  diGrid,          // 13 // inverse of grid spacing
    const float4   grid_p0          // 14 // origin of the grid
    //__global cl_Mat3* bboxes      // 15 // bounding box (xmin,ymin,zmin)(xmax,ymax,zmax)(kx,ky,kz)
);
```