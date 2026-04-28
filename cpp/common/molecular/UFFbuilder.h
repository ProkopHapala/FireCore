#ifndef UFFbuilder_h
#define UFFbuilder_h

//#include "MMFFBuilder.h"
#include "MMFFBuilderBase.h"
#include "UFF.h"

namespace MM {

class UFFBuilder : public BuilderBase {public:

    // ==== from BuilderBase ====
    //static int iDebug = 0;
    // std::vector<Atom>       atoms;
    // std::vector<Bond>       bonds;
    // std::vector<Angle>      angles;
    // std::vector<Dihedral>   dihedrals;
    // std::vector<Inversion>  inversions;
    // std::vector<AtomConf>   confs;
    // std::vector<Fragment>   frags;
    // //std::vector<int>  atom_neighs;
    // // --- other
    // bool  bPBC = false;
    // Mat3d lvec = Mat3dIdentity;
    // MMFFparams* params = 0;  // Each function which needs this can take it as parameter
    // std::unordered_set<int> capping_types;
    // // --------- Aux Params
    // // --- Types
    // int  itypHcap  =-1;
    // int  itypEpair =-1;
    // int  itypSigmaHole = -1;  // optional dedicated sigma-hole dummy type (e.g., E_h)
    // int  itype_min = 1; // type 0 is * (not valid type)
    // int  ignoreType=-1;
    // // --- Default params
    // Quat4d defaultREQ   {  1.5, 0.0, 0.0, 0.0    };
    // Bond   defaultBond  { -1, {-1,-1}, 1.5, 10.0 };
    // Angle  defaultAngle { -1, {-1,-1}, M_PI, 0.5 };
    // Bond   bondBrush = defaultBond;
    // // --- Caps
    // Atom   capAtom      = Atom{ (int)NeighType::H,     -1,-1, {0,0,0}, Atom::HcapREQ };
    // Atom   capAtomEpair = Atom{ (int)NeighType::epair, -1,-1, {0,0,0}, {0,0,0} };
    // Atom   capAtomPi    = Atom{ (int)NeighType::pi,    -1,-1, {0,0,0}, {0,0,0} };
    // Bond   capBond      = Bond{ -1,  {-1,-1},  1.07, 100/const_eVA2_Nm };
    // Vec3d  capUp        = Vec3d{0.0,0.0,1.0};
    // // --- bools
    // bool   bondByRvdW  = false;
    // bool   bDummyPi    = false;
    // bool   bDummyEpair = false;
    // bool   bAutoTypes  = true;
    // bool   bAddCaps    = true;
    // bool   bAddECaps   = false;

    // --- UFF-related functions from MMFFBuilder ---

    void toUFF( UFF& ff, bool bRealloc=true ){
        printf( "MM::Builder::toUFF() \n" );
    
        int natoms      = atoms.size();
        int nbonds      = bonds.size();
        int nangles     = angles.size();
        int ndihedrals  = dihedrals.size();
        int ninversions = inversions.size();
        if(verbosity>0) printf(  "MM::Builder::toUFF() natoms %i nbonds %i nangles %i ndihedrals %i ninversions %i\n", natoms, nbonds, nangles, ndihedrals, ninversions );
        if(bRealloc) ff.realloc( natoms, nbonds, nangles, ndihedrals, ninversions );
        
        for(int i=0; i<nbonds; i++ ){
            const Bond& B = bonds[i];
            ff.bonAtoms[i]  = Vec2i{B.atoms.x, B.atoms.y};
            ff.bonParams[i] = Vec2d{B.k, B.l0};
        }
        for(int i=0; i<nangles; i++ ){
            const Angle& A = angles[i];
            ff.angAtoms[i]  = Vec3i{A.atoms.x, A.atoms.y, A.atoms.z};
            ff.angParams[i] = double5{A.k, A.C0, A.C1, A.C2, A.C3};
        }
        for(int i=0; i<ndihedrals; i++ ){
            const Dihedral& D = dihedrals[i];
            ff.dihAtoms[i]  = Quat4i{D.atoms.x, D.atoms.y, D.atoms.z, D.atoms.w};
            ff.dihParams[i] = Vec3d{D.k, D.d, D.n};
        }
        for(int i=0; i<ninversions; i++ ){
            const Inversion& I = inversions[i];
            ff.invAtoms[i]  = Quat4i{I.atoms.x, I.atoms.y, I.atoms.z, I.atoms.w};
            ff.invParams[i] = Quat4d{I.k, I.C0, I.C1, I.C2};
        }
    
        for(int ia=0; ia<natoms; ia++ ){
            const Atom& A =  atoms[ia];
            ff.atypes[ia] = A.type;
            ff.apos [ia]  = A.pos; 
            ff.neighs[ia]=Quat4i{-1,-1,-1,-1};
            if(A.iconf>=0){
                AtomConf& conf = confs[A.iconf];
                int* ngs  = ff.neighs[ia].array;
                for(int k=0; k<conf.nbond; k++){
                    int ib = conf.neighs[k];
                    const Bond& B = bonds[ib];
                    int ja = B.getNeighborAtom(ia);
                    ngs[k] = ja;
                }
            }
        } 
    
        ff.makeNeighBs();         
        ff.bakeAngleNeighs();     
        ff.bakeDihedralNeighs();  
        ff.bakeInversionNeighs(); 
        ff.mapAtomInteractions(); 
        ff.printSizes();          
     
        ff.bPBC = bPBC;
        if(verbosity>0)printf("MM::Builder::toUFF DONE\n");
    
    }

    void assignUFFtypes_trivial( int* neighs, double* BOs, int* BOs_int, bool* set_atom, bool* set_bond ){

        for(int ia=0; ia<atoms.size(); ia++){
            Atom& Ai = atoms[ia];
            AtomConf& Ci = confs[Ai.iconf];
            if( params->atypes[Ai.type].name[0] == 'H' ){ // all hydrogens
                Ai.type = params->getAtomType("H_", true);
                set_atom[ia] = true;
                int ja = neighs[ia*4];
                int ib = getBondByAtoms( ia, ja );
                BOs[ib] = 1.0;
                BOs_int[ib] = 1;
                set_bond[ib] = true;
            } else {
                AtomConf& Ci = confs[Ai.iconf];
                if( params->atypes[Ai.type].name[0] == 'C' && Ci.nbond == 4 ){ // all sp3 carbons
                    Ai.type = params->getAtomType("C_3", true);
                    set_atom[ia] = true;
                    for(int in=ia*4; in<ia*4+4; in++){
                        int ja = neighs[in];
                        int ib = getBondByAtoms( ia, ja );
                        BOs[ib] = 1.0;
                        BOs_int[ib] = 1;
                        set_bond[ib] = true;
                    }    
                }else if( params->atypes[Ai.type].name[0] == 'N' && Ci.nbond == 3 ){ // all sp3 nitrogens
                    // not setting ufftype and final bond orders, as some may be resonant
                    for(int in=ia*4; in<ia*4+3; in++){
                        int ja = neighs[in];
                        int ib = getBondByAtoms( ia, ja );
                        BOs_int[ib] = 1;
                    }    
                }else if( params->atypes[Ai.type].name[0] == 'N' && Ci.nbond == 1 ){ // all sp1 nitrogens
                    Ai.type = params->getAtomType("N_1", true);
                    set_atom[ia] = true;
                    int ja = neighs[ia*4];
                    int ib = getBondByAtoms( ia, ja );
                    BOs[ib] = 3.0;
                    BOs_int[ib] = 3;
                    set_bond[ib] = true;
                    Atom& Aj = atoms[ja];
                    // nitrile group
                    if( params->atypes[Aj.type].name[0] == 'C' ){ 
                        Aj.type = params->getAtomType("C_1", true);
                        set_atom[ja] = true;
                        for(int jn=ja*4; jn<ja*4+2; jn++){
                            int ka = neighs[jn];
                            if( ka == ia ) continue;
                            int jb = getBondByAtoms( ja, ka );
                            BOs[jb] = 1.0;
                            BOs_int[jb] = 1;
                            set_bond[jb] = true;
                        }
                    }
                }else if( params->atypes[Ai.type].name[0] == 'O' && Ci.nbond == 2 ){ // all sp3 oxygens
                    // not setting ufftype and final bond orders, as some may be resonant
                    for(int in=ia*4; in<ia*4+2; in++){
                        int ja = neighs[in];
                        int ib = getBondByAtoms( ia, ja );
                        BOs_int[ib] = 1;
                    }    
                }else if( params->atypes[Ai.type].name[0] == 'O' && Ci.nbond == 1 ){ // all sp2 oxygens
                    // not setting ufftype and final bond orders, as some may be resonant
                    int ja = neighs[ia*4];
                    int ib = getBondByAtoms( ia, ja );
                    BOs_int[ib] = 2;
                }
            }
        }

    }

    void assignUFFtypes_nitro( int* neighs, double* BOs, int* BOs_int, bool* set_atom, bool* set_bond ){

        for(int ia=0; ia<atoms.size(); ia++){
            Atom& Ai = atoms[ia];
            if( params->atypes[Ai.type].name[0] == 'H' ) continue;
            AtomConf& Ci = confs[Ai.iconf];
            if( params->atypes[Ai.type].name[0] == 'N' && Ci.nbond == 3 ){ 
                int n = 0;
                for(int in=ia*4; in<ia*4+3; in++){
                    int ja = neighs[in];
                    Atom& Aj = atoms[ja];
                    if ( params->atypes[Aj.type].name[0] == 'O' ){n++;};
                }
                if( n == 2 ){
                    Ai.type = params->getAtomType("N_R", true);
                    set_atom[ia] = true;
                    for(int in=ia*4; in<ia*4+3; in++){
                        int ja = neighs[in];
                        Atom& Aj = atoms[ja];
                        if ( params->atypes[Aj.type].name[0] == 'O' ){
                            Aj.type = params->getAtomType("O_R", true);
                            set_atom[ja] = true;
                            int ib = getBondByAtoms( ia, ja );
                            BOs[ib] = 1.5;
                            set_bond[ib] = true;
                        }else{
                            int ib = getBondByAtoms( ia, ja );
                            BOs[ib] = 1.0;
                            set_bond[ib] = true;
                        }
                    }
                }
            }
        }

    }

    bool assignUFFtypes_checkall( int* neighs, int* BOs_int ){

        // check saturation
        for (; true; ) {
            bool changed = false;
            for(int ia=0; ia<atoms.size(); ia++){
                bool found = false;
                Atom& A = atoms[ia];
                if( params->atypes[A.type].name[0] == 'H' ) continue;
                AtomConf& C = confs[A.iconf];
                // skip already set atoms
                for(int in=ia*4; in<ia*4+C.nbond; in++){
                    int ja = neighs[in];
                    int ib = getBondByAtoms( ia, ja );
                    if( BOs_int[ib] < 0 ){ found = true; break; }
                }
                if( !found ) continue;
                // look for atoms that misses only one bond (or only single bonds)
                int nset = 0;
                int val = 0;
                for(int in=ia*4; in<ia*4+C.nbond; in++){
                    int ja = neighs[in];
                    int ib = getBondByAtoms( ia, ja );
                    if( BOs_int[ib] > 0 ){ nset++; val+=BOs_int[ib]; }
                }
                // if found, set it
                if( ( nset == C.nbond - 1 ) || ( C.nbond - nset == params->atypes[A.type].valence - val ) ){
                    changed = true;
                    if( nset == C.nbond-1 ){
                        for(int in=ia*4; in<ia*4+C.nbond; in++){
                            int ja = neighs[in];
                            int ib = getBondByAtoms( ia, ja );
                            if( BOs_int[ib] < 0 ){ BOs_int[ib] = params->atypes[A.type].valence - val; break; }
                        }
                    }else if( C.nbond - nset == params->atypes[A.type].valence - val ){
                        for(int in=ia*4; in<ia*4+C.nbond; in++){
                            int ja = neighs[in];
                            int ib = getBondByAtoms( ia, ja );
                            if( BOs_int[ib] < 0 ){ BOs_int[ib] = 1; }
                        }
                    }
                }
            }
            if (!changed) break;
        }

        // check if we are done
        bool done = true;
        // check if the bond orders are all set
        for (int ib = 0; ib < bonds.size(); ib++) {
            if ( BOs_int[ib] < 0 ) { done = false; break; }
        }
        return done;

    }

    bool assignUFFtypes_checkatom( int ia, int* neighs, int* BOs_int ){

        bool exclude = true;
        Atom& A = atoms[ia];
        AtomConf& C = confs[A.iconf];
        for(int in=ia*4; in<ia*4+C.nbond; in++){
            int ja = neighs[in];
            int ib = getBondByAtoms( ia, ja );
            if( BOs_int[ib] < 0 ){ exclude = false; break; }
        }
        return exclude;

    }

    void assignUFFtypes_treewalk( int* neighs, int* BOs_int ){

        // see if the problem is trivial
        bool done = assignUFFtypes_checkall ( neighs, &BOs_int[0] );
        if ( done ) return;

        // --- explore the graph
        // NB: at this point, only sp2 C, sp2 N and sp C (both triple bond and 2x double bond) are left to be assigned
        // TBD use arrays instead of vectors
        std::vector<int> BOs_int_tmp0(bonds.size());
        std::vector<int> BOs_int_tmp1(bonds.size());
        std::vector<int> BOs_int_tmp2(bonds.size());
        std::vector<int> BOs_int_tmp3(bonds.size());
        for (int ib = 0; ib < bonds.size(); ib++) { BOs_int_tmp0[ib] = BOs_int[ib]; }

        // pick an atom to attempt assigning a double bond
        for(int ia1=0; ia1<atoms.size(); ia1++){
            Atom& A1 = atoms[ia1];
            if( params->atypes[A1.type].name[0] == 'H' ) continue;
            if ( assignUFFtypes_checkatom( ia1, neighs, &BOs_int[0] ) ) continue;
            // pick a a neighbor
            AtomConf& C1 = confs[A1.iconf];
            for(int in1=ia1*4; in1<ia1*4+C1.nbond; in1++){
                int ib1 = getBondByAtoms( ia1, neighs[in1] );
                if( BOs_int_tmp0[ib1] > 0 ) continue;
                // try to set it up as a double bond
                for (int ib = 0; ib < bonds.size(); ib++) { BOs_int_tmp1[ib] = BOs_int_tmp0[ib]; }
                BOs_int_tmp1[ib1] = 2;
                // see if this solves the problem
                if ( assignUFFtypes_checkall ( neighs, &BOs_int_tmp1[0] ) ){ 
                    for (int ib = 0; ib < bonds.size(); ib++) { BOs_int[ib] = BOs_int_tmp1[ib]; } 
                    return; 
                }

                // otherwise, continue with the second nested attempt
                // pick an atom to attempt assigning a double bond
                for(int ia2=0; ia2<atoms.size(); ia2++){
                    Atom& A2 = atoms[ia2];
                    if( params->atypes[A2.type].name[0] == 'H' ) continue;
                    if ( assignUFFtypes_checkatom( ia2, neighs, &BOs_int[0] ) ) continue;
                    // pick a a neighbor
                    AtomConf& C2 = confs[A2.iconf];
                    for(int in2=ia2*4; in2<ia2*4+C2.nbond; in2++){
                        int ib2 = getBondByAtoms( ia2, neighs[in2] );
                        if( BOs_int_tmp1[ib2] > 0 ) continue;
                        // try to set it up as a double bond
                        for (int ib = 0; ib < bonds.size(); ib++) { BOs_int_tmp2[ib] = BOs_int_tmp1[ib]; }
                        BOs_int_tmp2[ib2] = 2;
                        // see if this solves the problem
                        if ( assignUFFtypes_checkall ( neighs, &BOs_int_tmp2[0] ) ){ 
                            for (int ib = 0; ib < bonds.size(); ib++) { BOs_int[ib] = BOs_int_tmp2[ib]; } 
                            return; 
                        }

                        // otherwise, continue with the third nested attempt
                        // pick an atom to attempt assigning a double bond
                        for(int ia3=0; ia3<atoms.size(); ia3++){
                            Atom& A3 = atoms[ia3];
                            if( params->atypes[A3.type].name[0] == 'H' ) continue;
                            if ( assignUFFtypes_checkatom( ia3, neighs, &BOs_int[0] ) ) continue;
                            // pick a a neighbor
                            AtomConf& C3 = confs[A3.iconf];
                            for(int in3=ia3*4; in3<ia3*4+C3.nbond; in3++){
                                int ib3 = getBondByAtoms( ia3, neighs[in3] );
                                if( BOs_int_tmp2[ib3] > 0 ) continue;
                                // try to set it up as a double bond
                                for (int ib = 0; ib < bonds.size(); ib++) { BOs_int_tmp3[ib] = BOs_int_tmp2[ib]; }
                                BOs_int_tmp3[ib3] = 2;
                                // see if this solves the problem
                                if ( assignUFFtypes_checkall ( neighs, &BOs_int_tmp3[0] ) ){ 
                                    for (int ib = 0; ib < bonds.size(); ib++) { BOs_int[ib] = BOs_int_tmp3[ib]; } 
                                    return; 
                                }
                            }
                        }
                    }
                }
            }
        }

        printf("ERROR TREEWALK: tree walk failed\n");
        printf("STOP\n");
        exit(0);

    }

    void assignUFFtypes_simplerule( double tol, int* neighs, double* BOs, bool* set_atom, bool* set_bond ){

        // simplest-est rule: atoms with two sp2 neighbors (or heteroatoms) are resonant
        for(int ia=0; ia<atoms.size(); ia++){
            Atom& A = atoms[ia];
            if( params->atypes[A.type].name[0] == 'H' ) continue;
            AtomConf& C = confs[A.iconf];
            if( ( params->atypes[A.type].name[0] == 'C' && C.nbond == 3 ) || // C_2
                ( params->atypes[A.type].name[0] == 'N' && C.nbond > 1 ) ||  // N_3 or N_2
                ( params->atypes[A.type].name[0] == 'O' && C.nbond == 2 ) ){ // O_3
                int n = 0;
                for(int in=ia*4; in<ia*4+C.nbond; in++){
                    int ja = neighs[in];
                    Atom& Aj = atoms[ja];
                    if( params->atypes[Aj.type].name[0] == 'H' ) continue;
                    AtomConf& Cj = confs[Aj.iconf];
                    if( ( params->atypes[Aj.type].name[0] == 'C' && Cj.nbond == 3 ) || // C_2
                        ( params->atypes[Aj.type].name[0] == 'N' && Cj.nbond > 1 ) ||  // N_3 or N_2
                        ( params->atypes[Aj.type].name[0] == 'O' && Cj.nbond == 2 ) ){ n++; }
                }
                if( n > 1 ){
                    // set atoms
                    if ( set_atom[ia] ) {
                        if ( params->atypes[A.type].name[2] != 'R' ) printf("WARNING SIMPLERULE: atom %i would be set to resonant but it has already a type of %s\n", ia+1, params->atypes[A.type].name);
                        continue;
                    }
                    if ( params->atypes[A.type].name[0] == 'C' ){ 
                        A.type = params->getAtomType("C_R", true); 
                    }else if ( params->atypes[A.type].name[0] == 'N' ){ 
                        A.type = params->getAtomType("N_R", true);
                    }else if ( params->atypes[A.type].name[0] == 'O' ){
                        A.type = params->getAtomType("O_R", true);
                    }else{
                        printf("ERROR SIMPLERULE: atom %i type %i not recognized\n", ia, A.type);
                        printf("STOP\n");
                        exit(0);
                    }
                    set_atom[ia] = true;
                    // set bonds
                    for(int in=ia*4; in<ia*4+C.nbond; in++){
                        int ja = neighs[in];
                        Atom& Aj = atoms[ja];
                        if( params->atypes[Aj.type].name[0] == 'H' ) continue;
                        AtomConf& Cj = confs[Aj.iconf];
                        if( ( params->atypes[Aj.type].name[0] == 'C' && Cj.nbond == 3 ) || // C_2
                            ( params->atypes[Aj.type].name[0] == 'N' && Cj.nbond > 1 ) ||  // N_3 or N_2
                            ( params->atypes[Aj.type].name[0] == 'O' && Cj.nbond == 2 ) ){ // O_3
                            int ib = getBondByAtoms( ia, neighs[in] );
                            if ( set_bond[ib] ){
                                if ( abs(BOs[ib]-1.5) > tol ) printf("WARNING SIMPLERULE: bond %i between atoms %i and %i would be set to 1.5 but it has already a bond order of %g\n", ib+1, ia+1, ja+1, BOs[ib]);
                                continue;
                            }
                            BOs[ib] = 1.5;
                            set_bond[ib] = true;
                        }
                    }
                }
            }
        }
    }

    void assignUFFtypes_findrings( double tol, int* neighs, double* BOs, int* BOs_int, bool* set_atom, bool* set_bond ){
        bool   explored[atoms.size()];
        for(int ia=0; ia<atoms.size(); ia++){explored[ia]=false;}
        std::vector<int>  nb(8);
        for (; true; ) {
            bool changed = false;
            for(int ia=0; ia<atoms.size(); ia++){
                nb[1] = ia; // 1st atom
                Atom& A1 = atoms[nb[1]];
                if ( params->atypes[A1.type].name[0] == 'H' ) continue;
                AtomConf& C1 = confs[A1.iconf];
                for(int in1=nb[1]*4; in1<nb[1]*4+C1.nbond; in1++){
                    nb[2] = neighs[in1]; // 2nd atom
                    Atom& A2 = atoms[nb[2]];
                    if ( params->atypes[A2.type].name[0] == 'H' ) continue;
                    AtomConf& C2 = confs[A2.iconf];
                    for(int in2=nb[2]*4; in2<nb[2]*4+C2.nbond; in2++){
                        nb[3] = neighs[in2]; // 3rd atom
                        if ( nb[3] == nb[1] ) continue;
                        Atom& A3 = atoms[nb[3]];
                        if ( params->atypes[A3.type].name[0] == 'H' ) continue;
                        AtomConf& C3 = confs[A3.iconf];
                        for(int in3=nb[3]*4; in3<nb[3]*4+C3.nbond; in3++){
                            nb[4] = neighs[in3]; // 4th atom
                            if ( nb[4] == nb[2] ) continue;
                            Atom& A4 = atoms[nb[4]];
                            if ( params->atypes[A4.type].name[0] == 'H' ) continue;
                            AtomConf& C4 = confs[A4.iconf];
                            for(int in4=nb[4]*4; in4<nb[4]*4+C4.nbond; in4++){
                                nb[5] = neighs[in4]; // 5th atom
                                if ( nb[5] == nb[3] ) continue;
                                Atom& A5 = atoms[nb[5]];
                                if ( params->atypes[A5.type].name[0] == 'H' ) continue;
                                AtomConf& C5 = confs[A5.iconf];
                                for(int in5=nb[5]*4; in5<nb[5]*4+C5.nbond; in5++){
                                    nb[6] = neighs[in5]; // 6th atom
                                    if ( nb[6] == nb[4] ) continue;

                                    // found 5-member ring
                                    if ( nb[6] == nb[1] ) {
                                        if ( explored[nb[1]] && explored[nb[2]] && explored[nb[3]] && explored[nb[4]] && explored[nb[5]] ) continue;
                                        nb[0] = nb[5];
                                        // count double bonds & heteroatoms
                                        if ( assignUFFtypes_checkaroma ( 5, &nb[0], &BOs[0], &BOs_int[0], tol ) ) {
                                            changed = true;
                                            for(int i=1; i<6; i++){explored[nb[i]] = true;}
                                            assignUFFtypes_setaroma ( 5, &nb[0], &neighs[0], &BOs[0], &set_atom[0], &set_bond[0], tol );
                                        }
                                    }

                                    Atom& A6 = atoms[nb[6]];
                                    if ( params->atypes[A6.type].name[0] == 'H' ) continue;
                                    AtomConf& C6 = confs[A6.iconf];
                                    for(int in6=nb[6]*4; in6<nb[6]*4+C6.nbond; in6++){
                                        nb[7] = neighs[in6]; // 7th atom
                                        if ( nb[7] == nb[5] ) continue;

                                        // found 6-member ring
                                        if ( nb[7] == nb[1] ) {
                                            if ( explored[nb[1]] && explored[nb[2]] && explored[nb[3]] && explored[nb[4]] && explored[nb[5]] && explored[nb[6]] ) continue;
                                            nb[0] = nb[6];
                                            // count double bonds & heteroatoms
                                            if ( assignUFFtypes_checkaroma ( 6, &nb[0], &BOs[0], &BOs_int[0], tol ) ) {
                                                changed = true;
                                                for(int i=1; i<7; i++){explored[nb[i]] = true;}
                                                assignUFFtypes_setaroma ( 6, &nb[0], &neighs[0], &BOs[0], &set_atom[0], &set_bond[0], tol );
                                            }
                                        }

                                    } // in6
                                } // in5
                            } // in4
                        } // in3
                    } // in2
                } // in1
            } // ia

            if (!changed) break;
        }

    }

    bool assignUFFtypes_checkaroma ( int n, int* nb, double* BOs, int* BOs_int, double tol ){

        // check that all carbons in the ring are sp2
        for(int i=1; i<n+1; i++){
            Atom& A = atoms[nb[i]];
            AtomConf& C = confs[A.iconf];
            if ( params->atypes[A.type].name[0] == 'C' && C.nbond != 3 ) return false;
        }
                
        int npi = 0;

        // count heteroatoms in the ring that can donate an electron pair        
        for(int i=1; i<n+1; i++){
            Atom& A = atoms[nb[i]];
            AtomConf& C = confs[A.iconf];
            if ( ( params->atypes[A.type].name[0] == 'N' && C.nbond == 3 ) || 
                 ( params->atypes[A.type].name[0] == 'O' && C.nbond == 2 ) ) npi += 2;
        }

        int npi_save = npi;

        // count pi electrons in the ring considering double bonds in the ring only
        for(int i=1; i<n+1; i++){
            int ib1 = getBondByAtoms( nb[i], nb[i+1] );
            int ib2 = getBondByAtoms( nb[i-1], nb[i] );
            if ( BOs_int[ib1] == 2 || BOs_int[ib2] == 2 ) npi += 1;
        }

        if ( npi == 6 ) return true;

        npi = npi_save;

        // count pi electrons in the ring considering also fused rings
        for(int i=1; i<n+1; i++){
            int ib1 = getBondByAtoms( nb[i], nb[i+1] );
            int ib2 = getBondByAtoms( nb[i-1], nb[i] );
            if ( BOs_int[ib1] == 2 || abs(BOs[ib1]-1.5) < tol || BOs_int[ib2] == 2 || abs(BOs[ib2]-1.5) < tol ) npi += 1;
        }

        if ( npi == 6 ) return true;

        return false;

    }

    void assignUFFtypes_setaroma ( int n, int* nb, int* neighs, double* BOs, bool* set_atom, bool* set_bond, double tol ){

        // set atoms
        for(int i=1; i<n+1; i++){
            Atom& A = atoms[nb[i]];
            if ( set_atom[nb[i]] ) {
                printf("WARNING SETAROMA: atom %c %i would be set to resonant but it has already a type of %i\n", params->atypes[A.type].name[0], nb[i]+1, A.type);
                continue;
            }
            if ( params->atypes[A.type].name[0] == 'C' ){ 
                A.type = params->getAtomType("C_R", true); 
            }else if ( params->atypes[A.type].name[0] == 'N' ){ 
                A.type = params->getAtomType("N_R", true);
            }else if ( params->atypes[A.type].name[0] == 'O' ){
                A.type = params->getAtomType("O_R", true);
            }else{
                printf("ERROR SIMPLERULE: atom %i type %i not recognized\n", nb[i]+1, A.type);
                printf("STOP\n");
                exit(0);
            }
            set_atom[nb[i]] = true;
        }
        // set ring bonds
        for(int i=1; i<n+1; i++){
            int ib = getBondByAtoms( nb[i], nb[i+1] );
            if ( set_bond[ib] && abs(BOs[ib]-1.5) > tol ) {
                Atom& A1 = atoms[nb[i]];
                Atom& A2 = atoms[nb[i+1]];
                printf("WARNING SETAROMA: bond %i between atoms %c %i and %c %i would be set to 1.5 but it has already a bond order of %g\n", ib+1, params->atypes[A1.type].name[0], nb[i]+1, params->atypes[A2.type].name[0], nb[i+1]+1, BOs[ib]);
                continue;
            }
            BOs[ib] = 1.5;
            set_bond[ib] = true;
        }

        // set other bonds
        for(int i=1; i<n+1; i++){
            Atom& A = atoms[nb[i]];
            if ( params->atypes[A.type].name[0] == 'H' ) continue;
            AtomConf& C = confs[A.iconf];
            for(int in=nb[i]*4; in<nb[i]*4+C.nbond; in++){
                bool ok = true;
                for(int j=1; j<n+1; j++){ if ( neighs[in] == nb[j] ) { ok = false; break; } }
                if ( !ok ) continue;
                Atom& Aj = atoms[neighs[in]];
                if ( params->atypes[Aj.type].name[0] == 'H' ) continue;
                AtomConf& Cj = confs[Aj.iconf];
                int ib = getBondByAtoms( nb[i], neighs[in] );
                if ( params->atypes[Aj.type].name[0] == 'O' && Cj.nbond == 1 ) { // delocalized carbonyl
                    if ( set_atom[neighs[in]] && params->atypes[Aj.type].name[2] != 'R' ) {
                        printf("WARNING SETAROMA: carbonyl atom %c %i would be set to resonant but it has already a type of %i\n", params->atypes[Aj.type].name[0], neighs[in]+1, Aj.type);
                        continue;
                    }
                    if ( set_bond[ib] && abs(BOs[ib]-1.5) > tol ) {
                        printf("WARNING SETAROMA: bond %i between atoms %c %i and %c %i would be set to 1.5 but it has already a bond order of %g\n", ib+1, params->atypes[A.type].name[0], nb[i]+1, params->atypes[Aj.type].name[0], neighs[in]+1, BOs[ib]);
                        continue;
                    }
                    Aj.type = params->getAtomType("O_R", true);
                    set_atom[neighs[in]] = true;
                    BOs[ib] = 1.5;
                    set_bond[ib] = true;
                } //else {
                  //  BOs[ib] = 1.0;
                  //  set_bond[ib] = true; }
            }
        }

    }

    void assignUFFtypes_assignrest( int* neighs, double* BOs, int* BOs_int, bool* set_atom, bool* set_bond ){

        for(int ia=0; ia<atoms.size(); ia++){
            if ( set_atom[ia] ) continue;
            Atom& A = atoms[ia];
            if ( params->atypes[A.type].name[0] == 'H' ) continue;
            AtomConf& C = confs[A.iconf];
            if ( params->atypes[A.type].name[0] == 'C' ) {
                if ( C.nbond == 2 ) { 
                    A.type = params->getAtomType("C_1", true); 
                    set_atom[ia] = true; 
                } else if ( C.nbond == 3 ) { 
                    A.type = params->getAtomType("C_2", true); 
                    set_atom[ia] = true; 
                }
            } else if ( params->atypes[A.type].name[0] == 'N' ) {
                if ( C.nbond == 2 ) { 
                    A.type = params->getAtomType("N_2", true); 
                    set_atom[ia] = true; 
                } else if ( C.nbond == 3 ) { 
                    A.type = params->getAtomType("N_3", true); 
                    set_atom[ia] = true;
                    for(int in=ia*4; in<ia*4+3; in++){
                        int ib = getBondByAtoms( ia, neighs[in] );
                        BOs[ib] = 1.0;
                        BOs_int[ib] = 1;
                        set_bond[ib] = true;
                    }
                }
            } else if ( params->atypes[A.type].name[0] == 'O' ) {
                if ( C.nbond == 1 ) { 
                    A.type = params->getAtomType("O_2", true); 
                    set_atom[ia] = true; 
                    int ib = getBondByAtoms( ia, neighs[ia*4] );
                    BOs[ib] = 2.0;
                    BOs_int[ib] = 2;
                    set_bond[ib] = true;
                    Atom& Aj = atoms[neighs[ia*4]];
                    if ( params->atypes[Aj.type].name[0] == 'H' ) continue;
                    AtomConf& Cj = confs[Aj.iconf];
                    if ( params->atypes[A.type].name[0] == 'C' && !set_atom[neighs[ia*4]] && Cj.nbond == 3 ) {
                        Aj.type = params->getAtomType("C_2", true);
                        set_atom[neighs[ia*4]] = true;
                    }
                } else if ( C.nbond == 2 ) { 
                    A.type = params->getAtomType("O_3", true); 
                    set_atom[ia] = true; 
                    for(int in=ia*4; in<ia*4+2; in++){
                        int ib = getBondByAtoms( ia, neighs[in] );
                        BOs[ib] = 1.0;
                        BOs_int[ib] = 1;
                        set_bond[ib] = true;
                    }
                }
            }
        }

    }

    void assignUFFtypes_fixsaturation( int* neighs, double* BOs, int* BOs_int, bool* set_bond, double tol ){

        for (; true; ) {
            bool changed = false;

            for(int ia=0; ia<atoms.size(); ia++){
                Atom& A = atoms[ia];
                if ( params->atypes[A.type].name[0] == 'H' ) continue;
                AtomConf& C = confs[A.iconf];
                // if the atom has all bonds set, skip
                bool ok = true;
                for(int in=ia*4; in<ia*4+C.nbond; in++){
                    int ib = getBondByAtoms( ia, neighs[in] );
                    if ( BOs[ib] < 0.0 ) { ok = false; break; }
                }
                if ( ok ) continue;
                // compute number of already set bonds and corresponding atom valence
                int nset = 0;
                double val = 0.0;
                for(int in=ia*4; in<ia*4+C.nbond; in++){
                    int ib = getBondByAtoms( ia, neighs[in] );
                    if ( BOs[ib] > 0.0 ) { 
                        nset++; 
                        val+=(double)BOs_int[ib]; 
                    }
                }
                val = (double)params->atypes[A.type].valence - val;
                // only one bond to be set
                if ( nset == C.nbond-1 ) {
                    // check for weird valence
                    if ( abs( val - (double)rint(val) ) > tol || abs(val) < tol ) {
                        printf("ERROR FIXSATURATION: atom %c %i would have a valence of %g\n", params->atypes[A.type].name[0], ia+1, val);
                        printf("STOP\n");
                        exit(0);
                    }
                    changed = true;
                    for(int in=ia*4; in<ia*4+C.nbond; in++){
                        int ib = getBondByAtoms( ia, neighs[in] );
                        if ( BOs[ib] < 0.0 ) {
                            BOs[ib] = val;
                            BOs_int[ib] = (int)rint(val);
                            set_bond[ib] = true;
                            break;
                        }
                    }
                // or multiple bonds to be set, but they are all single
                } else if ( (int)rint(val) == C.nbond-nset ) {
                        changed = true;
                        for(int in=ia*4; in<ia*4+C.nbond; in++){
                            int ib = getBondByAtoms( ia, neighs[in] );
                            if ( BOs[ib] < 0.0 ) {
                                BOs[ib] = 1.0;
                                BOs_int[ib] = 1;
                                set_bond[ib] = true;
                            }
                        }
                }
            }
                
            if( !changed ) break;
        }

    }

    void assignUFFtypes_cumulene( int* neighs, double* BOs, int* BOs_int, bool* set_bond ){

        // ad hoc exception for getting rid of cumulenes in the backbone, i.e.
        // ... C(ia3)                                           ... C(ia3)                          
        //       \\                              /                    \\                __            /     
        //        C(ia4) == C(ia1) == C(ia2) == C(ia3)	    ==>       C(ia4) -- C(ia1) == C(ia2) -- C(ia3) 
        //       /                              \\                     /                              \\           
        //                                        C(ia4) ...                                           C(ia4) ... 
        for(int ia1=0; ia1<atoms.size(); ia1++){
            Atom& A1 = atoms[ia1];
            if( params->atypes[A1.type].name[0] != 'C' ) continue; // looking for C
            AtomConf& C1 = confs[A1.iconf];
            if( C1.nbond != 2 ) continue; // looking for sp C
            for(int in1=ia1*4; in1<ia1*4+2; in1++){
                int ia2 = neighs[in1];
                Atom& A2 = atoms[ia2];
                if( params->atypes[A2.type].name[0] != 'C' ) continue; // looking for C
                AtomConf& C2 = confs[A2.iconf];
                if( C2.nbond != 2 ) continue; // looking for sp C
                int ib1 = getBondByAtoms( ia1, ia2 );
                int ia4;
                if(in1==ia1*4){ia4=neighs[ia1*4+1];}else{ia4=neighs[ia1*4];}
                Atom& A4 = atoms[ia4];
                if( params->atypes[A4.type].name[0] != 'C' ) continue; // looking for C
                AtomConf& C4 = confs[A4.iconf];
                if( C4.nbond != 3 ) continue; // looking for sp2 C
                int ib4 = getBondByAtoms( ia1, ia4 );
                for(int in2=ia2*4; in2<ia2*4+2; in2++){
                    int ia3 = neighs[in2];
                    if( ia3 == ia1 ) continue; // looking for the other neighbor of ia2
                    Atom& A3 = atoms[ia3];
                    if( params->atypes[A3.type].name[0] != 'C' ) continue; // looking for C
                    AtomConf& C3 = confs[A3.iconf];
                    if( C3.nbond != 3 ) continue; // looking for sp2 C
                    for(int in3=ia3*4; in3<ia3*4+3; in3++){
                        if( neighs[in3] != ia4 ) continue; // looking for the right ia3
                        int ib2 = getBondByAtoms( ia2, ia3 );
                        int ib3 = getBondByAtoms( ia3, ia4 );
                        BOs[ib1] = 3.0;
                        BOs[ib2] = 1.0;
                        BOs[ib3] = 2.0;
                        BOs[ib4] = 1.0;
                        BOs_int[ib1] = 3;
                        BOs_int[ib2] = 1;
                        BOs_int[ib3] = 2;
                        BOs_int[ib4] = 1;
                        set_bond[ib1] = true;
                        set_bond[ib2] = true;
                        set_bond[ib3] = true;
                        set_bond[ib4] = true;
                    }
                }
            }
        }
    }

    void assignUFFtypes_conjugation( int* neighs, double* BOs ){

        for(int ia=0; ia<atoms.size(); ia++){
            Atom& A = atoms[ia];
            if ( params->atypes[A.type].name[0] == 'H' ) continue;
            AtomConf& C = confs[A.iconf];
            if( (params->atypes[A.type].name[0]=='N'&&C.nbond==3) || (params->atypes[A.type].name[0]=='O'&&C.nbond==2) ){
                for(int in=ia*4; in<ia*4+C.nbond; in++){
                    int ja = neighs[in];
                    Atom& Aj = atoms[ja];
                    if ( params->atypes[Aj.type].name[2]=='R' ) {
                        if ( params->atypes[A.type].name[0] == 'N' ) { A.type = params->getAtomType("N_R", true); } 
                        else if ( params->atypes[A.type].name[0] == 'O' ) { A.type = params->getAtomType("O_R", true); }
                        int ib = getBondByAtoms( ia, ja );
                        BOs[ib] = 1.5;
                    }
                }
            }
        }

    }

    void assignUFFtypes_checks( int* neighs, double* BOs, int* BOs_int, bool* set_atom, bool* set_bond, double tol ){

        // sanity checks
        for(int ia=0; ia<atoms.size(); ia++){
            if(!set_atom[ia]){printf("ERROR CHECKS: atom %i is not set\n", ia+1);printf("STOP\n");exit(0);}
        }
        for(int ib=0; ib<bonds.size(); ib++){
            if(!set_bond[ib]){printf("ERROR CHECKS: bond %i is not set\n", ib+1);printf("STOP\n");exit(0);}
            if(BOs[ib]<0.0){printf("ERROR CHECKS: order for bond %i order is not set\n", ib+1);printf("STOP\n");exit(0);}
            if(BOs_int[ib]<0){printf("ERROR CHECKS: integer order for bond %i order is not set\n", ib+1);printf("STOP\n");exit(0);}
        }

        for(int ia=0; ia<atoms.size(); ia++){
            Atom& A = atoms[ia];
            // check that resonant atoms have at least one 1.5 bond
            if ( params->atypes[A.type].name[2]=='R' ) { 
                AtomConf& C = confs[A.iconf]; 
                bool found = false;
                for(int in=ia*4; in<ia*4+C.nbond; in++){
                    int ib = getBondByAtoms( ia, neighs[in] );
                    if ( abs(BOs[ib]-1.5) < tol ) { found = true; break; }
                }
                if ( !found ) {
                    printf("ERROR CHECKS: atom %c %i is resonant but it has no 1.5 bonds\n", params->atypes[A.type].name[0], ia+1);
                    for(int in=ia*4; in<ia*4+C.nbond; in++){
                        int ib = getBondByAtoms( ia, neighs[in] );
                        printf("  bond %i atoms %c %i %c %i bond order %g %i\n", ib+1, params->atypes[A.type].name[0], ia+1, params->atypes[atoms[neighs[in]].type].name[0], neighs[in]+1, BOs[ib], BOs_int[ib]);
                    }
                    printf("STOP\n");
                    exit(0);
                }
            // check that sp2 atoms have one localized double bond
            } else if ( params->atypes[A.type].name[2]=='2' ) { 
                AtomConf& C = confs[A.iconf]; 
                bool found = false;
                for(int in=ia*4; in<ia*4+C.nbond; in++){
                    int ib = getBondByAtoms( ia, neighs[in] );
                    if ( abs(BOs[ib]-2.0) < tol ) { 
                        if ( found ) {
                            printf("ERROR CHECKS: atom %c %i is sp2 but it has more than one double bond\n", params->atypes[A.type].name[0], ia+1);
                            for(int in=ia*4; in<ia*4+C.nbond; in++){
                                int ib = getBondByAtoms( ia, neighs[in] );
                                printf("  bond %i atoms %c %i %c %i bond order %g %i\n", ib+1, params->atypes[A.type].name[0], ia+1, params->atypes[atoms[neighs[in]].type].name[0], neighs[in]+1, BOs[ib], BOs_int[ib]);
                            }
                            printf("STOP\n");
                            exit(0);
                        }
                        found = true;
                    }
                }
                if ( !found ) {
                    found = false;
                    for(int in=ia*4; in<ia*4+C.nbond; in++){
                        int ib = getBondByAtoms( ia, neighs[in] );
                        if ( abs(BOs[ib]-1.5) < tol ) { found = true; break; }
                    }
                    if ( found ) {
                        printf("WARNING CHECKS: atom %c %i is sp2 and it has no double bonds, only 1.5 bonds\n", params->atypes[A.type].name[0], ia+1);
                    } else {
                        printf("ERROR CHECKS: atom %c %i is sp2 but it has no double bonds\n", params->atypes[A.type].name[0], ia+1);
                        for(int in=ia*4; in<ia*4+C.nbond; in++){
                            int ib = getBondByAtoms( ia, neighs[in] );
                            printf("  bond %i atoms %c %i %c %i bond order %g %i\n", ib+1, params->atypes[A.type].name[0], ia+1, params->atypes[atoms[neighs[in]].type].name[0], neighs[in]+1, BOs[ib], BOs_int[ib]);
                        }
                        printf("STOP\n");
                        exit(0);
                    }
                }
            }
        }

        // check that 1.5 bonds must be either 1 or 2 in the limit resonance structure
        for(int ib=0; ib<bonds.size(); ib++){
            if ( abs(BOs[ib]-1.5) < tol ) {
                if ( BOs_int[ib] != 1 && BOs_int[ib] != 2 ) {
                    printf("ERROR CHECKS: bond %i is 1.5 but in the limit resonance structure is neither single nor double\n", ib+1);
                    printf("STOP\n");
                    exit(0);
                }
            }
        }
    }

    void assignUFFtypes_amide( int* neighs, double* BOs ){

        for(int ia=0; ia<atoms.size(); ia++){
            Atom& A = atoms[ia];
            if ( params->atypes[A.type].name[0] == 'H' ) continue;
            AtomConf& C = confs[A.iconf];
            if( params->atypes[A.type].name[0] == 'C' && C.nbond == 3 ) {
                // look for a carbonyl oxygen
                bool found = false;
                for(int in=ia*4; in<ia*4+C.nbond; in++){
                    int ja = neighs[in];
                    Atom& Aj = atoms[ja];
                    if ( params->atypes[Aj.type].name[0] == 'H' ) continue;
                    AtomConf& Cj = confs[Aj.iconf];
                    if( params->atypes[Aj.type].name[0] == 'O' && Cj.nbond == 1 ) { found = true; break; }
                }
                if ( !found ) continue;
                // look for an amino nitrogen
                found = false;
                for(int in=ia*4; in<ia*4+C.nbond; in++){
                    int ja = neighs[in];
                    Atom& Aj = atoms[ja];
                    if ( params->atypes[Aj.type].name[0] == 'H' ) continue;
                    AtomConf& Cj = confs[Aj.iconf];
                    if( params->atypes[Aj.type].name[0] == 'N' && Cj.nbond == 3 ) { found = true; break; }
                }
                if ( !found ) continue;
                A.type = params->getAtomType("C_R", true);
                for(int in=ia*4; in<ia*4+C.nbond; in++){
                    int ja = neighs[in];
                    Atom& Aj = atoms[ja];
                    if ( params->atypes[Aj.type].name[0] == 'H' ) continue;
                    AtomConf& Cj = confs[Aj.iconf];
                    if( params->atypes[Aj.type].name[0] == 'O' && Cj.nbond == 1 ) {
                        Aj.type = params->getAtomType("O_R", true);
                        int ib = getBondByAtoms( ia, ja );
                        BOs[ib] = 1.5;
                    } else if( params->atypes[Aj.type].name[0] == 'N' && Cj.nbond == 3 ) {
                        Aj.type = params->getAtomType("N_R", true);
                        int ib = getBondByAtoms( ia, ja );
                        BOs[ib] = 1.41;
                    }
                }
            }
        }

    }

    // UFF atom-type assignement
    void assignUFFtypes( int* neighs=0, bool bCumulene=false, bool bDeallocNeighs=true, bool b141=true, bool bSimple=true, bool bConj=false ){ 

        // init working arrays
        double tol = 0.05;
        bool   set_atom[atoms.size()];
        bool   set_bond[bonds.size()];
        std::vector<double> BOs(bonds.size());
        std::vector<int>    BOs_int(bonds.size());
        for(int ia=0; ia<atoms.size(); ia++){set_atom[ia]=false;}
        for(int ib=0; ib<bonds.size(); ib++){BOs[ib]=-1.0; BOs_int[ib]=-1; set_bond[ib]=false;}

        bDeallocNeighs &= (neighs==0);
        makeNeighs ( neighs, 4 );

        // assign bond orders and atom types for trivial cases
        // e.g. sp3 carbons, hydrogens, etc...
        assignUFFtypes_trivial( neighs, &BOs[0], &BOs_int[0], &set_atom[0], &set_bond[0] );

        // assign bond orders and atom types for nitro groups
        assignUFFtypes_nitro( neighs, &BOs[0], &BOs_int[0], &set_atom[0], &set_bond[0] );

        // find a (hopefully) valid limit resonance structure
        assignUFFtypes_treewalk( neighs, &BOs_int[0] );

        if ( bSimple ) {
            // assign resonant atoms according to the "simple" rule: atoms with two sp2 neighbors (or heteroatoms) are resonant
            assignUFFtypes_simplerule( tol, neighs, &BOs[0], &set_atom[0], &set_bond[0] );
        } else {
            // assign resonant atoms according to the Huckel rule for aromaticity
            assignUFFtypes_findrings( tol, neighs, &BOs[0], &BOs_int[0], &set_atom[0], &set_bond[0] );
        }

        // assign the rest of atom types
        assignUFFtypes_assignrest( neighs, &BOs[0], &BOs_int[0], &set_atom[0], &set_bond[0] );

        // try to assign the rest of double bonds
        assignUFFtypes_fixsaturation( neighs, &BOs[0], &BOs_int[0], &set_bond[0], tol );

        // exception to avoid cumulenes
        if( bCumulene ) { assignUFFtypes_cumulene( neighs, &BOs[0], &BOs_int[0], &set_bond[0] ); }
  
        // manually change sp3 nitrogen and oxygen to "resonant" when they are bonded to an sp2 atom (conjugation)
        if( bConj ){ assignUFFtypes_conjugation( neighs, &BOs[0] ); }

        // some check before assigning parameters
        assignUFFtypes_checks( neighs, &BOs[0], &BOs_int[0], &set_atom[0], &set_bond[0], tol );

        // exception for amide groups
        assignUFFtypes_amide( neighs, &BOs[0] );

        // store bond orders
        for(int ib=0; ib<bonds.size(); ib++){
            Bond& B = bonds[ib];
            B.order = BOs[ib];
            //if(verbosity>0)printf( "bondOrder[%i](%i-%i)= %g\n", ib+1, B.atoms.i+1, B.atoms.j+1, B.order );            
        }

        // deallocate neighbor array
        if(bDeallocNeighs)delete [] neighs;

    }

    double assignUFFparams_calcrij ( int ib ){

        Bond& B = bonds[ib];   // reference
        //Bond* B_ = &bonds[ib]; // pointer

        AtomType& Ti = params->atypes[atoms[B.atoms.i].type];
        AtomType& Tj = params->atypes[atoms[B.atoms.j].type];
        const ElementType& Ei = *params->elementOfAtomType(atoms[B.atoms.i].type);
        const ElementType& Ej = *params->elementOfAtomType(atoms[B.atoms.j].type);
        double rBO = -0.1332 * ( Ti.Ruff + Tj.Ruff ) * log( B.order );
        double rEN = Ti.Ruff * Tj.Ruff * sq( sqrt(-Ei.Eaff) - sqrt(-Ej.Eaff) ) / ( -Ei.Eaff*Ti.Ruff -Ej.Eaff*Tj.Ruff );
        double rIJ = Ti.Ruff + Tj.Ruff + rBO - rEN;
        return rIJ;

    }

    void assignUFFparams_vdws( ){

        for( int ia=0; ia<atoms.size(); ia++){
            Atom& A = atoms[ia];
            params->assignRE( A.type, A.REQ, true );
        }

    }

    void assignUFFparams_bonds( ){

        for ( int ib=0; ib<bonds.size(); ib++ ) {
            Bond& B = bonds[ib];
            const ElementType& Ei = *params->elementOfAtomType(atoms[B.atoms.i].type);
            const ElementType& Ej = *params->elementOfAtomType(atoms[B.atoms.j].type);
            B.l0 = assignUFFparams_calcrij(ib);
            B.k = 0.5 * 28.7989689090648 * Ei.Quff * Ej.Quff / ( B.l0*sq(B.l0) );
        }

    }

    void assignUFFparams_angles( int* neighs ){

        angles.clear();
        for( int j=0; j<atoms.size(); j++){
            Atom& Aj = atoms[j];
            AtomType& Tj = params->atypes[Aj.type];
            double ct = cos(Tj.Ass*deg2rad);
            double st2 = sq(sin(Tj.Ass*deg2rad));
            if ( params->atypes[Aj.type].name[0] == 'H' ) continue;
            AtomConf& Cj = confs[Aj.iconf];
            for(int in1=j*4; in1<j*4+Cj.nbond-1; in1++){
                int i = neighs[in1];
                const ElementType& Ei = *params->elementOfAtomType(atoms[i].type);
                for(int in2=in1+1; in2<j*4+Cj.nbond; in2++){
                    int k = neighs[in2];
                    const ElementType& Ek = *params->elementOfAtomType(atoms[k].type);
                    Angle a;
                    a.atoms.x = i;
                    a.atoms.y = j;
                    a.atoms.z = k;
                    a.bonds.x = getBondByAtoms( i, j );
                    a.bonds.y = getBondByAtoms( j, k );
                    double rij = assignUFFparams_calcrij(a.bonds.x);
                    double rjk = assignUFFparams_calcrij(a.bonds.y);
                    double rik = sqrt( rij*rij + rjk*rjk - 2.0*rij*rjk*ct );
                    double kappa = 28.7989689090648 * Ei.Quff * Ek.Quff / (sq(rik)*sq(rik)*rik) * ( 3.0*rij*rjk*st2 - sq(rik)*ct );
                    if ( params->atypes[Aj.type].name[2]=='1' || params->atypes[Aj.type].name[2]=='2' || 
                    params->atypes[Aj.type].name[2]=='R' ) { // cosine/periodic
                        if ( params->atypes[Aj.type].name[2]=='1' ) {
                            //a.k = 0.5 * kappa;
                            a.k = kappa;
                            a.C0 = 1.0;
                            a.C1 = 1.0;
                            a.C2 = 0.0;
                            a.C3 = 0.0;
                        } else if ( params->atypes[Aj.type].name[2]=='2' || params->atypes[Aj.type].name[2]=='R' ) {
                            a.k = kappa / 9.0;
                            a.C0 = 1.0;
                            a.C1 = 0.0;
                            a.C2 = 0.0;
                            a.C3 = -1.0;
                        } 
                    } else if ( params->atypes[Aj.type].name[2]=='3' ) { // fourier
                        a.k = kappa;
                        a.C2 = 1.0 / ( 4.0 * st2 );
                        a.C1 = -4.0 * a.C2 * ct;
                        a.C0 = a.C2 * ( 2.0*sq(ct) + 1.0 );
                        a.C3 = 0.0;
                    }
                    angles.push_back(a);
                }
            }
        }

    }

    void assignUFFparams_dihedrals( int* neighs ){

        dihedrals.clear();
        for( int i1=0; i1<atoms.size(); i1++){
            Atom& A1 = atoms[i1];
            int n1;
            if ( params->atypes[A1.type].name[0] == 'H' ){ n1 = 1; } 
            else { AtomConf& C1 = confs[A1.iconf]; n1 = C1.nbond; }
            for(int in1=i1*4; in1<i1*4+n1; in1++){
                int i2 = neighs[in1];
                Atom& A2 = atoms[i2];
                if ( params->atypes[A2.type].name[0] == 'H' ) continue;
                if ( params->atypes[A2.type].name[2] == '1' ) continue; // Torsional potentials for central bonds involving sp-hybridized centers X-1 were assigned a value of zero
                AtomConf& C2 = confs[A2.iconf];
                for(int in2=i2*4; in2<i2*4+C2.nbond; in2++){
                    int i3 = neighs[in2];
                    if ( i3 != i1 ) {
                        Atom& A3 = atoms[i3];
                        if ( params->atypes[A3.type].name[0] == 'H' ) continue;
                        if ( params->atypes[A3.type].name[2] == '1' ) continue; // Torsional potentials for central bonds involving sp-hybridized centers X-1 were assigned a value of zero
                        AtomConf& C3 = confs[A3.iconf];
                        for(int in3=i3*4; in3<i3*4+C3.nbond; in3++){
                            int i4 = neighs[in3];
                            if ( i4 != i2 && i4 > i1 ) { // avoid 3-membered rings and double-counting
                                Atom& A4 = atoms[i4];
                                Dihedral d;
                                d.atoms.x = i1;
                                d.atoms.y = i2;
                                d.atoms.z = i3;
                                d.atoms.w = i4;
                                d.bonds.x = getBondByAtoms( i1, i2 );
                                d.bonds.y = getBondByAtoms( i2, i3 );
                                d.bonds.z = getBondByAtoms( i3, i4 );
                                const ElementType& E2 = *params->elementOfAtomType(atoms[i2].type);
                                const ElementType& E3 = *params->elementOfAtomType(atoms[i3].type);
                                // specific general case (a): * - sp3 - sp3 - *
                                if ( params->atypes[A2.type].name[2] == '3' && params->atypes[A3.type].name[2] == '3' ){
                                    d.k = sqrt( E2.Vuff * E3.Vuff );
                                    d.d = 1;
                                    d.n = 3;
                                    // special case of * - group 16 sp3 - group 16 sp3 - *
                                    if ( (params->atypes[A2.type].name[0]=='O'||params->atypes[A2.type].name[0]=='S') &&
                                        (params->atypes[A3.type].name[0]=='O'||params->atypes[A3.type].name[0]=='S') ){
                                            d.k = 4.1840/60.2214076/1.602176634; // 1 kcal/mol to eV
                                            d.n = 2;
                                            if ( params->atypes[A2.type].name[0] == 'O' ){ d.k *= 2.0; }
                                            else { d.k *= 6.8; }
                                            if ( params->atypes[A3.type].name[0] == 'O' ){ d.k *= 2.0; }
                                            else { d.k *= 6.8; }
                                            d.k = sqrt( d.k );
                                    }
                                } 
                                // specific general case (b): * - sp3 - sp2 - *
                                else if ( (params->atypes[A2.type].name[2]=='3'&&(params->atypes[A3.type].name[2]=='2'||params->atypes[A3.type].name[2]=='R')) ||
                                ((params->atypes[A2.type].name[2]=='2'||params->atypes[A2.type].name[2]=='R')&&params->atypes[A3.type].name[2]=='3') ){
                                    d.k = 4.1840/60.2214076/1.602176634; // 1 kcal/mol to eV
                                    d.d = -1;
                                    d.n = 6;
                                    // special case of * - group 16 sp3 - sp2 - *
                                    if ( (params->atypes[A2.type].name[2]=='3'&&(params->atypes[A2.type].name[0]=='O'||params->atypes[A2.type].name[0]=='S')) ||
                                        (params->atypes[A3.type].name[2]=='3'&&(params->atypes[A3.type].name[0]=='O'||params->atypes[A3.type].name[0]=='S')) ){
                                        int ib = getBondByAtoms( i2, i3 );
                                        Bond& B = bonds[ib];
                                        d.k = 5.0 * sqrt( E2.Uuff * E3.Uuff ) * ( 1.0 + 4.18 * log( B.order ) );
                                        d.d = 1;
                                        d.n = 2;
                                    }
                                    // special case of * - sp3 - sp2 bounded to another sp2 atom
                                    if ( params->atypes[A2.type].name[2]=='3' ){
                                        bool found = false;
                                        for(int in=i3*4; in<i3*4+C3.nbond; in++){
                                            Atom& A = atoms[neighs[in]];
                                            if ( params->atypes[A.type].name[2]=='2' || params->atypes[A.type].name[2]=='R' ) {
                                                found = true;
                                                break;
                                            }
                                        }
                                        if ( found ) {
                                            d.k = 2.0 * 4.1840/60.2214076/1.602176634; // 2 kcal/mol to eV
                                            d.d = 1;
                                            d.n = 3;
                                        }
                                    } else if ( params->atypes[A3.type].name[2]=='3' ){
                                        bool found = false;
                                        for(int in=i2*4; in<i2*4+C2.nbond; in++){
                                            Atom& A = atoms[neighs[in]];
                                            if ( params->atypes[A.type].name[2]=='2' || params->atypes[A.type].name[2]=='R' ) {
                                                found = true;
                                                break;
                                            }
                                        }
                                        if ( found ) {
                                            d.k = 2.0 * 4.1840/60.2214076/1.602176634; // 2 kcal/mol to eV
                                            d.d = 1;
                                            d.n = 3;
                                        }
                                    }
                                }
                                // specific general case (c): * - sp2 - sp2 - *
                                else if ( (params->atypes[A2.type].name[2]=='2'||params->atypes[A2.type].name[2]=='R') &&
                                        (params->atypes[A3.type].name[2]=='2'||params->atypes[A3.type].name[2]=='R') ){
                                    int ib = getBondByAtoms( i2, i3 );
                                    Bond& B = bonds[ib];
                                    d.k = 5.0 * sqrt( E2.Uuff * E3.Uuff ) * ( 1.0 + 4.18 * log( B.order ) );
                                    d.d = -1;
                                    d.n = 2;
                                } else {
                                    printf("ERROR: assignUFFparams_dihedrals: torsion case not found for atoms %s %s %s %s\n", params->atypes[A1.type].name, params->atypes[A2.type].name, params->atypes[A3.type].name, params->atypes[A4.type].name);
                                    printf("STOP\n");
                                    exit(0);
                                }
                                d.k = 0.5 * d.k / ( (double)(C2.nbond-1) * (double)(C3.nbond-1) );
                                dihedrals.push_back(d);
                            }
                        }
                    }
                }
            }
        }

    }       

    void assingUFFparams_assigninversion( int i1, int i2, int i3, int i4 ){

        Atom& A1 = atoms[i1];
        Atom& A2 = atoms[i2];
        Atom& A3 = atoms[i3];
        Atom& A4 = atoms[i4];

        Inversion i;
        i.atoms.x = i1;
        i.atoms.y = i2;
        i.atoms.z = i3;
        i.atoms.w = i4;
        i.bonds.x = getBondByAtoms( i1, i2 );
        i.bonds.y = getBondByAtoms( i1, i3 );
        i.bonds.z = getBondByAtoms( i1, i4 );
        // sp2 carbon
        if ( params->atypes[A1.type].name[0] == 'C' && (params->atypes[A1.type].name[2]=='2'||params->atypes[A1.type].name[2]=='R') ){
            // carbonyl
            if ( ( params->atypes[A2.type].name[0] == 'O' && params->atypes[A2.type].name[2] == '2' ) ||
                 ( params->atypes[A3.type].name[0] == 'O' && params->atypes[A3.type].name[2] == '2' ) ||
                 ( params->atypes[A4.type].name[0] == 'O' && params->atypes[A4.type].name[2] == '2' ) ){
                i.k = 50.0 * 4.1840/60.2214076/1.602176634; } // 50 kcal/mol to eV
            else { i.k = 6.0 * 4.1840/60.2214076/1.602176634; } // 6 kcal/mol to eV
            i.C0 = 1.0;
            i.C1 = -1.0;
            i.C2 = 0.0;
        }
        // sp2 nitrogen
        else if ( params->atypes[A1.type].name[0] == 'N' && (params->atypes[A1.type].name[2]=='2'||params->atypes[A1.type].name[2]=='R') ){
            i.k = 6.0 * 4.1840/60.2214076/1.602176634; // 6 kcal/mol to eV
            i.C0 = 1.0;
            i.C1 = -1.0;
            i.C2 = 0.0;
        }
        // sp3 nitrogen
        else if ( params->atypes[A1.type].name[0] == 'N' && params->atypes[A1.type].name[2]=='3' ){
            i.k = 0.0;
            i.C0 = 0.0;
            i.C1 = 0.0;
            i.C2 = 0.0;
        }
        // sp3 group 15
        else if ( params->atypes[A1.type].name[0] == 'P' || params->atypes[A1.type].name[0] == '3' ){
            double omega0 = 84.4339 * deg2rad;
            i.C0 = 4.0 * sq(cos(omega0)) - cos(2.0*omega0);
            i.C1 = -4.0 * cos(omega0);
            i.C2 = 1.0;
            i.k = 22.0 * 4.1840/60.2214076/1.602176634 / ( i.C0 + i.C1 + i.C2 ); // 22 kcal/mol to eV
        } else {
            printf("ERROR: assignUFFparams_inversions: improper case not found for atoms %s %s %s %s\n", params->atypes[A1.type].name, params->atypes[A2.type].name, params->atypes[A3.type].name, params->atypes[A4.type].name);
            printf("STOP\n");
            exit(0);
        }
        i.k = i.k / 3.0;
        inversions.push_back(i);

    }

    void assignUFFparams_inversions( int* neighs ){

        inversions.clear();
        for( int i1=0; i1<atoms.size(); i1++){
            Atom& A1 = atoms[i1];
            if ( params->atypes[A1.type].name[0] == 'H' ) continue;
            AtomConf& C1 = confs[A1.iconf]; 
            if ( C1.nbond != 3 ) continue;
            int i2 = neighs[i1*4];
            int i3 = neighs[i1*4+1];
            int i4 = neighs[i1*4+2];
            assingUFFparams_assigninversion( i1, i2, i3, i4 );
            assingUFFparams_assigninversion( i1, i4, i2, i3 );
            assingUFFparams_assigninversion( i1, i3, i4, i2 );
        }

    }

    void assignUFFparams( int* neighs=0, bool bDeallocNeighs=true ){

        bDeallocNeighs &= (neighs==0);
        makeNeighs ( neighs, 4 );

        // assign vdw parameters
        assignUFFparams_vdws();   

        // assign bond parameters
        assignUFFparams_bonds();

        // populate angle array & assign angle parameters
        assignUFFparams_angles( neighs );
   
        // assign torsion parameters
        assignUFFparams_dihedrals( neighs );
    
        // assign inversion parameters
        assignUFFparams_inversions( neighs );

        // deallocate neighbor array
        if(bDeallocNeighs)delete [] neighs;

        // write data
        //assignUFFparams_writedata();
          
    }
  
    void assignUFFparams_writedata( ){

        double tokcal = 60.2214076*1.602176634/4.1840;
        //double toeV = 4.1840/60.2214076/1.602176634;
        FILE* f=fopen( "mol.data.firecore", "w");
        fprintf( f, "\n" );
        fprintf( f, "\n" );
        fprintf( f, "%i atoms\n", atoms.size() );
        fprintf( f, "%i atom types\n", atoms.size() );
        fprintf( f, "%i bonds\n", bonds.size() );
        fprintf( f, "%i bond types\n", bonds.size() );
        fprintf( f, "%i angles\n", angles.size() );
        fprintf( f, "%i angle types\n", angles.size() );
        fprintf( f, "%i dihedrals\n", dihedrals.size() );
        fprintf( f, "%i dihedral types\n", dihedrals.size() );
        fprintf( f, "%i impropers\n", inversions.size() );
        fprintf( f, "%i improper types\n", inversions.size() );
        fprintf( f, "\n" );
        fprintf( f, "0.0 64.0 xlo xhi\n" );
        fprintf( f, "0.0 64.0 ylo yhi\n" );
        fprintf( f, "0.0 64.0 zlo zhi\n" );
        fprintf( f, "0.0 0.0 0.0 xy xz yz\n" );
        fprintf( f, "\n" );
        fprintf( f, "Masses\n" );
        fprintf( f, "\n" );
        for( int ia=0; ia<atoms.size(); ia++){
            Atom& A = atoms[ia];
            fprintf( f, "%i 1.0 # %s\n", ia+1, params->atypes[A.type].name );
        }
        fprintf( f, "\n" );
        fprintf( f, "Pair Coeffs\n" );
        fprintf( f, "\n" );
        for( int ia=0; ia<atoms.size(); ia++){
            Atom& A = atoms[ia];
            fprintf( f, "%i %f %f # %s %s\n", ia+1, sq(A.REQ.y)*tokcal, A.REQ.x*2.0*0.890898718140339, params->atypes[A.type].name, params->atypes[A.type].name );
        }
        if ( bonds.size() > 0 ){
            fprintf( f, "\n" );
            fprintf( f, "Bond Coeffs\n" );
            fprintf( f, "\n" );
            for( int ib=0; ib<bonds.size(); ib++){
                Bond& B = bonds[ib];
                fprintf( f, "%i %f %f # %s %s\n", ib+1, 0.5*B.k*tokcal, B.l0, params->atypes[atoms[B.atoms.i].type].name, params->atypes[atoms[B.atoms.j].type].name );
            }
        }
        if ( angles.size() > 0 ){
            fprintf( f, "\n" );
            fprintf( f, "Angle Coeffs\n" );
            fprintf( f, "\n" );
            bool bFourier = false;
            bool bCosine = false;
            for( int ia=0; ia<angles.size(); ia++){
                Angle& A = angles[ia];
                int i = A.atoms.j;
                if ( params->atypes[atoms[i].type].name[2] == '1' || params->atypes[atoms[i].type].name[2] == '2' || params->atypes[atoms[i].type].name[2] == 'R' ) { bCosine = true; }
                if ( params->atypes[atoms[i].type].name[2] == '3' ) { bFourier = true; }
            }
            for( int ia=0; ia<angles.size(); ia++){
                Angle& A = angles[ia];
                int i = A.atoms.j;
                if ( bFourier && bCosine ) { 
                    if ( params->atypes[atoms[i].type].name[2] == '1' ) { // cosine/periodic
                        fprintf( f, "%i cosine/periodic %f 1 1 # %s %s %s\n", ia+1, A.k*tokcal, params->atypes[atoms[A.atoms.x].type].name, params->atypes[atoms[A.atoms.y].type].name, params->atypes[atoms[A.atoms.z].type].name );
                    } else if ( params->atypes[atoms[i].type].name[2] == '2' || params->atypes[atoms[i].type].name[2] == 'R' ) { // cosine/periodic
                        fprintf( f, "%i cosine/periodic %f -1 3 # %s %s %s\n", ia+1, 9.0*A.k*tokcal, params->atypes[atoms[A.atoms.x].type].name, params->atypes[atoms[A.atoms.y].type].name, params->atypes[atoms[A.atoms.z].type].name );
                    } else if ( params->atypes[atoms[i].type].name[2] == '3' ) { // fourier
                        fprintf( f, "%i fourier %f %f %f %f # %s %s %s\n", ia+1, A.k*tokcal, A.C0, A.C1, A.C2, params->atypes[atoms[A.atoms.x].type].name, params->atypes[atoms[A.atoms.y].type].name, params->atypes[atoms[A.atoms.z].type].name );
                    }
                } else {
                    if ( params->atypes[atoms[i].type].name[2] == '1' ) { // cosine/periodic
                        fprintf( f, "%i %f 1 1 # %s %s %s\n", ia+1, A.k*tokcal, params->atypes[atoms[A.atoms.x].type].name, params->atypes[atoms[A.atoms.y].type].name, params->atypes[atoms[A.atoms.z].type].name );
                    } else if ( params->atypes[atoms[i].type].name[2] == '2' || params->atypes[atoms[i].type].name[2] == 'R' ) { // cosine/periodic
                        fprintf( f, "%i %f -1 3 # %s %s %s\n", ia+1, 9.0*A.k*tokcal, params->atypes[atoms[A.atoms.x].type].name, params->atypes[atoms[A.atoms.y].type].name, params->atypes[atoms[A.atoms.z].type].name );
                    } else if ( params->atypes[atoms[i].type].name[2] == '3' ) { // fourier
                        fprintf( f, "%i %f %f %f %f # %s %s %s\n", ia+1, A.k*tokcal, A.C0, A.C1, A.C2, params->atypes[atoms[A.atoms.x].type].name, params->atypes[atoms[A.atoms.y].type].name, params->atypes[atoms[A.atoms.z].type].name );
                    }
                }
            }
        }
        if ( dihedrals.size() > 0 ){
            fprintf( f, "\n" );
            fprintf( f, "Dihedral Coeffs\n" );
            fprintf( f, "\n" );
            for( int id=0; id<dihedrals.size(); id++){
                Dihedral& D = dihedrals[id];
                fprintf( f, "%i %f %i %i # %s %s %s %s\n", id+1, D.k*tokcal, D.d, D.n, params->atypes[atoms[D.atoms.x].type].name, params->atypes[atoms[D.atoms.y].type].name, params->atypes[atoms[D.atoms.z].type].name, params->atypes[atoms[D.atoms.w].type].name );
            }
        }
        if ( inversions.size() > 0 ){
            fprintf( f, "\n" );
            fprintf( f, "Improper Coeffs\n" );
            fprintf( f, "\n" );
            for( int ii=0; ii<inversions.size(); ii++){
                Inversion& I = inversions[ii];
                fprintf( f, "%i %f %f %f %f 0 # %s %s %s %s\n", ii+1, I.k*tokcal, I.C0, I.C1, I.C2, params->atypes[atoms[I.atoms.x].type].name, params->atypes[atoms[I.atoms.y].type].name, params->atypes[atoms[I.atoms.z].type].name, params->atypes[atoms[I.atoms.w].type].name );
            }
        }
        fprintf( f, "\n" );
        fprintf( f, "Atoms\n" );
        fprintf( f, "\n" );
        for( int ia=0; ia<atoms.size(); ia++){
            Atom& A = atoms[ia];
            fprintf( f, "%i 1 %i 0.0 %g %g %g\n", A.id+1, ia+1, A.pos.x, A.pos.y, A.pos.z );
        }
        if ( bonds.size() > 0 ){
            fprintf( f, "\n" );
            fprintf( f, "Bonds\n" );
            fprintf( f, "\n" );
            for( int ib=0; ib<bonds.size(); ib++){
                Bond& B = bonds[ib];
                fprintf( f, "%i %i %i %i\n", ib+1, ib+1, atoms[B.atoms.i].id+1, atoms[B.atoms.j].id+1 );
            }
        }
        if ( angles.size() > 0 ){
            fprintf( f, "\n" );
            fprintf( f, "Angles\n" );
            fprintf( f, "\n" );
            for( int ia=0; ia<angles.size(); ia++){
                Angle& A = angles[ia];
                fprintf( f, "%i %i %i %i %i\n", ia+1, ia+1, atoms[A.atoms.i].id+1, atoms[A.atoms.j].id+1, atoms[A.atoms.k].id+1 );
            }
        }
        if ( dihedrals.size() > 0 ){
            fprintf( f, "\n" );
            fprintf( f, "Dihedrals\n" );
            fprintf( f, "\n" );
            for( int id=0; id<dihedrals.size(); id++){
                Dihedral& D = dihedrals[id];
                fprintf( f, "%i %i %i %i %i %i\n", id+1, id+1, atoms[D.atoms.x].id+1, atoms[D.atoms.y].id+1, atoms[D.atoms.z].id+1, atoms[D.atoms.w].id+1 );
            }
        }
        if ( inversions.size() > 0 ){
            fprintf( f, "\n" );
            fprintf( f, "Impropers\n" );
            fprintf( f, "\n" );
            for( int ii=0; ii<inversions.size(); ii++){
                Inversion& I = inversions[ii];
                fprintf( f, "%i %i %i %i %i %i\n", ii+1, ii+1, atoms[I.atoms.x].id+1, atoms[I.atoms.y].id+1, atoms[I.atoms.z].id+1, atoms[I.atoms.w].id+1 );
            }
        }
        fclose(f);

    }

};

}

#endif
