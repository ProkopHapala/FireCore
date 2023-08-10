
import sys
import os
import numpy as np
#import psi4
from . import atomicUtils as au

resp = None
psi4 = None

#try:
#    import resp
#except Exception as e:
#    print("ERROR: Cannot open `resp` module for psi4")
#    print(e)


# https://psicode.org/psi4manual/master/api/psi4.core.Molecule.html
# https://github.com/psi4/psi4numpy/blob/master/Tutorials/01_Psi4NumPy-Basics/1b_molecule.ipynb    
#    * Example: Fitting Lennard-Jones Parameters from Potential Energy Scan

# see also : https://notebook.community/dsirianni/psi4_api-docs/01_interactive_Psi4_Basics


# ============ Setup


default_resp_options = {
'VDW_SCALE_FACTORS'  : [1.4, 1.6, 1.8, 2.0],
'VDW_POINT_DENSITY'  : 1.0,
'RESP_A'             : 0.0005,
'RESP_B'             : 0.1,
}

default_psi4_options = {
"geom_maxiter": 100,                # increase limit for geometry relaxation
"intrafrag_step_limit"    : 0.1,    # this helps with geometry relaxation convergence
"intrafrag_step_limit_min": 0.1,
"intrafrag_step_limit_max": 0.1,
"opt_coordinates" : "cartesian",
"step_type":  "nr"
}


def load_res():
    global resp
    if resp is None:
        import resp

# ============ Functions

def try_make_dirs( dname ):
    try:
        os.mkdir( dname )
    except:
        pass

def xyz2str(fname):
    ls = [ " ".join( l.split()[:4] ) for l in open(fname).readlines()[2:] ]
    #print(ls)
    return '\n'.join(ls)

def save_xyz_Q( fname, lines, Qs ):
    n = len(Qs)
    with open(fname,'w') as fout:
        fout.write("%i\n" %n)
        fout.write("#comment \n" )
        for i,Q in enumerate( Qs ):
            fout.write( "%s %10.5f \n" %(lines[i],  Q) )

#def preparemol(fname='input.xyz', geom_str='O 0. 0. 0.\n H 1. 0. 0.\n H 0. 1. 0.' ):
def preparemol(fname='input.xyz', geom_str=None ):
    load_res()
    if os.path.isfile(fname):
        geom_str  = xyz2str( fname )     #;print("geom>>%s<<" %geom )
    mol   = psi4.geometry( geom_str )
    mol.update_geometry()
    return mol 

def unpack_mol( mol ):
    na = mol.natom()
    apos = np.zeros((na,3))
    es   = [ str(mol.symbol(i)) for i in range(na) ]
    for i in range(na):
        p = mol.xyz(i)
        apos[i,0]=p[0]
        apos[i,1]=p[1]
        apos[i,2]=p[2]
    #print( "es  \n", es )
    #print( "pos \n", apos )
    return apos, es

def pack_mol( apos, es, ifrag_line=None ):
    print( "pack_mol es="+str(es)+" apos="+str(apos) )
    load_res()
    na = len(es)
    strs = [  "%s %g %g %g" %(es[i],apos[i,0],apos[i,1],apos[i,2]) for i in range(na) ]    #;print( strs )
    if ifrag_line is not None: strs.insert( ifrag_line, "--" )
    geom = "\n".join(strs)         ;print(geom)
    mol   = psi4.geometry( geom )
    return mol

def eval( geom, params=None, id=None ):
    load_res()
    pars = params.copy()
    method = pars['method'];   del pars['method'];
    basis  = pars['basis' ];   del pars['basis' ];
    bsse   = pars['bsse'  ];   del pars['bsse' ];  
     
    # ---- Fragments
    ifrag_line=None
    if 'ifrag_line' in params.keys():
        ifrag_line = params.get('ifrag_line'); del pars['ifrag_line']

    #print( method, basis  )
    method_basis=method+"/"+basis
    # ------ load geometry
    apos,es = geom
    mol = pack_mol( apos, es, ifrag_line=ifrag_line )
    mol.update_geometry()
    mol.symmetrize(1e-3)   # this heps prevent problems with symmetry : https://github.com/psi4/psi4webmo/issues/4
    psi4.set_options( pars )
        
    E = psi4.energy(method_basis, molecule=mol, bsse_type=bsse )

    return E

def relax( geom=None, params=None, fname=None ):
    load_res()
    if  geom is not None:
        apos,es = geom
        mol     = pack_mol( apos, es )
    elif fname is not None:
        smol = xyz2str( fname ) 
        mol  = psi4.geometry( smol )
    
    pars = params.copy()
    method = pars['method']
    basis  = pars['basis']
    del pars['method']
    del pars['basis']
    #print( method, basis  )
    method_basis=method+"/"+basis
    # ------ load geometry

    mol.update_geometry()
    mol.symmetrize(1e-3)   # this heps prevent problems with symmetry : https://github.com/psi4/psi4webmo/issues/4
    psi4.set_options( pars )
    psi4.optimize(method_basis, molecule=mol)
    #E = psi4.energy(method_basis, molecule=mol)
    return mol

def psi4resp( name, bRelax=True, indir="./input/", outdir="./output/", method='scf', basis='/STO-3G', resp_options=default_resp_options, psi4_options=default_psi4_options ):
    load_res()
    if resp is None:
        import resp
    method_basis = method+"/"+basis;    #print( method_basis )
    # ------ load geometry
    geom  = xyz2str( indir+name+".xyz")     #;print("geom>>%s<<" %geom )
    mol   = psi4.geometry( geom )
    mol.update_geometry()
    mol.symmetrize(1e-3)   # this heps prevent problems with symmetry : https://github.com/psi4/psi4webmo/issues/4

    psi4.set_options( psi4_options )
    
    # ------ relax geometry
    if bRelax:
        psi4.optimize(method_basis, molecule=mol)

    # ----- save output
    geom_lines = mol.save_string_xyz().split('\n')[1:]

    # ------------ Call for first stage fit
    resp_options['METHOD_ESP'] = method
    resp_options['BASIS_ESP' ] = basis
    Qs = resp.resp([mol], resp_options)   ;Q_esp=Qs[0]; Q_resp=Qs[1]
    print('ESP  Charges: ', Q_esp  )
    print('RESP Charges: ', Q_resp )
    save_xyz_Q( outdir+name+".xyz", geom_lines, Q_resp )

'''
def unpack_mol( mol, units=bohr2A ):
    apos= np.array([ a[1] for a in mol._atom ]) * units
    es  = np.array([ a[0] for a in mol._atom ])
    return apos, es

def pack_mol( apos, es ):
    return [ (es[i],apos[i]) for i in range(len(es))]

def printlist(lst):
    for item in lst: print( item )

def printObj(obj):
    printlist(dir(obj))

def saveAtoms(fname,atoms, unit=bohr2A ):
    apos = np.array([ a[1] for a in atoms ])*unit
    es   = [ a[0] for a in atoms ]
    au.saveXYZ( es, apos, fname )

def preparemol(fname='relaxed.xyz', conv_params=None, atoms='O 0 0 0; H 1 0 0; H 0 1 0' ):
    if conv_params is None: conv_params=default_conv_params
    if os.path.isfile(fname):
        print("found(%s) => no need for relaxation " %fname )
        mol = pyscf.M(atom=fname)
    else:
        h2o = pyscf.M(atom=atoms)
        h2o.verbose = verbosity
        calc = pyscf.scf.RHF(h2o)
        #calc = pyscf.dft.RKS(h2o); calc.xc = 'b3lyp'; #calc.init_guess='atom'
        mol = optimize(calc, maxsteps=1000, **conv_params)
        saveAtoms(fname,mol._atom)
    return mol 

def evalHf(inp):
    apos,es = inp            #;print( apos, es )
    m = pack_mol( apos, es)  #;print( m )
    mol=pyscf.M( atom=pack_mol( apos, es) )   #;print( " mol._atom \n", mol.atom )
    #mol.build()             #;print( mol._atom )
    mol.verbose = verbosity 
    out = pyscf.scf.UHF(mol).run()
    #if(verbosity>0): print( x, out.e_tot, out.energy_elec(), out.energy_nuc() )
    return out.e_tot*hartree2eV

def optHf(atoms, conv_params=None ):
    if conv_params is None: conv_params=default_conv_params
    print(atoms)
    job = pyscf.M(atom=atoms)
    job.SCF_max_cycle=100
    job.verbose = verbosity
    calc = pyscf.scf.RHF(job)
    mol = optimize(calc, maxsteps=1000, **conv_params)
    printlist(mol)
    return mol

'''




default_params_block='''
    scf_type df 
    opt_type MIN
    geom_maxiter 1000
    g_convergence qchem
    print_trajectory_xyz_file true
    opt_coordinates cartesian
    step_type nr
''' 

def extract_final_geom(fname):
    s = os.popen('grep -n "Final optimized geometry" %s | cut -b -10' %fname ).read()
    nstart = int( s.split(':')[0] )
    lines = []
    fin = open(fname,'r') 
    for i in range(nstart+5):
        next(fin)
    for l in fin:
        ws = l.split()
        if( len(ws)<4 ): break
        lines.append(l)
    fin.close()
    return lines

def write_geom( fname, lines, comment="#comment" ):
    fout = open(fname,'w')
    #print(len(lines))
    fout.write( "%i\n" %len(lines) )
    fout.write( comment )
    for l in lines:
        #print(l)
        fout.write(l)
    fout.close()

def extract_input_geom( fname ):
    s = os.popen('grep -n "molecule " %s | cut -b -10' %fname ).read()
    nstart = int( s.split(':')[0] )
    fin = open(fname,'r')
    lines = []
    for i in range(nstart+1):
        next(fin)
    for l in fin:
        ws = l.split()
        nw = len(ws)
        if(nw==4):
            lines.append( l )
        elif(nw==1):
            nhyphen = len( lines )
        else:    
            break
    fin.close()
    return lines, nhyphen 

def write_psi4_in( lines, nhyphen=None, mem='500MB', method='b3lyp', basis='CC-pVDZ', bsse="'cp'", params_block=None, fname='psi.in', q=0, multiplicity=1, opt=True ):
    fout = open(fname,'w')
    fout.write( 'memory '+mem+"\n" )
    fout.write( "molecule {\n" )
    fout.write( "%i %i\n" %(q,multiplicity) )
    il=0
    for l in lines:
        fout.write(l)
        if(il==nhyphen):
            fout.write('  --\n')
        il+=1
    fout.write( "units angstrom\n" )
    fout.write( "}\n" )
    fout.write( "set {\n" )
    fout.write( "    basis " + basis+"\n" )
    if params_block is not None:
        fout.write( params_block )
    fout.write( "}\n" )
    if (nhyphen is None):
        bsse=None
    if opt:
        if bsse is None:
            fout.write( "optimize( '%s')\n" %method )
        else:
            fout.write( "optimize( '%s', bsse_type=%s)\n" %(method,bsse) )
    else:
        if bsse is None:
             fout.write( "energy('%s')\n" %method )
        else:
             fout.write( "energy( '%s', bsse_type=%s )\n" %(method,bsse) )
    fout.close()
    