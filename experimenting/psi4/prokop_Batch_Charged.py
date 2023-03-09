import psi4
import resp
import os
import time

# ========= Setup

indir="./input/"
#outdir="./output/"

#bRelax=True
bRelax=False

#method='scf'     # Hartree-Fock
#method='pbe' 
#method='b3lyp'
#method='mp2'
#basis='sto-3g'
#basis='6-31+G'
#basis='cc-pvdz'

methods=[
#    'scf',
#    'pbe',
#    'b3lyp',
    'mp2'
]
basises=[
    # 'sto-3g',
    # '6-31+G',
    # '6-311+G*',
    # '6-311++G**',
    # '6-311++G(3df,3pd)',
     'cc-pvdz',
    # 'aug-cc-pvtz',
    # 'def2-QZVPPD',
]


resp_options = {
'VDW_SCALE_FACTORS'  : [1.4, 1.6, 1.8, 2.0],
'VDW_POINT_DENSITY'  : 1.0,
'RESP_A'             : 0.0005,
'RESP_B'             : 0.1,
}

psi4_options = {
"geom_maxiter": 100,                # increase limit for geometry relaxation
"intrafrag_step_limit"    : 0.1,    # this helps with geometry relaxation convergence
"intrafrag_step_limit_min": 0.1,
"intrafrag_step_limit_max": 0.1,
"opt_coordinates" : "cartesian",
"step_type":  "nr"
}




names=[
# ---- Minimal
#"HF",
#"HCCH",
#"HCN",
#"NH3",
#"H2O",
#"oxalate",
#"F-COOH",
#"NC-COOH",
"HCOOH",
#"OCH2",
#"formaldimine",
# ------- Other
#"acetic_acid",
#"acetimide",
#"acetamide",
#"amino_ethan",
#"hydroxo_ethylene",
#"amino_acetylene",
#"hydroxo_acetylene",
# ------- Heterocycles
#"furan",
#"pyrrole",
#"pyridine",
# ------- Quinone
#"NN-quinine",     # SCF NOT CONVERGE !!!
#"quinone",
#"hydroquinone",  
# ------- Aromatic
#"phenol",
#"aninline",
#"benzaldehyde",
#"benzoic_acid",
]


# ======== Functions

#def file2string(fname):
#    ls = [ l for l in open(fname) ]
#    return '\n'.join(ls[2:])

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
    #print("lines ", lines)
    #print("Qs ",Qs)
    with open(fname,'w') as fout:
        fout.write("%i\n" %n)
        fout.write("#comment \n" )
        for i,Q in enumerate( Qs ):
            fout.write( "%s %10.5f \n" %(lines[i],  Q) )

def process_molecule( name, bRelax=True, indir="./input/", outdir="./output/", method='scf', basis='/STO-3G' ):
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

# ======== Main

#psi4.core.be_quiet()

for method in methods:
    for basis in basises:
        outdir = method+"/"+basis+"/"
        try_make_dirs( method )
        try_make_dirs( outdir )
        for name in names:
            print( "# ======= Molecule: ", name, method, basis )
            psi4.core.set_output_file(outdir+name+'.log', False)
            t0 = time.time_ns()
            try:
                process_molecule( name, bRelax=bRelax, method=method, basis=basis, indir=indir, outdir=outdir )
            except Exception as e: 
                print(e)
            t = time.time_ns() - t0; print( "time: ", t*1e-9, "[s]" )