import psi4
import resp

# ========= Setup

fname='H2O.xyz'
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
#"HCOOH",
#"OCH2",
"formaldimine",
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

def file2string(fname):
    ls = [ l for l in open(fname) ]
    return '\n'.join(ls[2:])

def save_xyz_Q( fname, lines, Qs ):
    n = len(Qs)
    print("lines ", lines)
    print("Qs ",Qs)
    with open(fname,'w') as fout:
        fout.write("%i\n" %n)
        fout.write("#comment \n" )
        for i,Q in enumerate( Qs ):
            fout.write( "%s %10.5f \n" %(lines[i],  Q) )

def process_molecule( name, bRelax=True ):
    print( "# ======= Molecule: ", name )
    # ------ load geometry
    #geom = file2string("./molecules/"+name)
    geom  = file2string("./molecules/"+name+".xyz")
    mol   = psi4.geometry( geom )
    #mol  = psi4.Molecule.init_with_xyz(filen) 
    mol.update_geometry()
    
    psi4.set_options( psi4_options )

    if bRelax:
        #psi4.optimize('scf/cc-pvdz', molecule=mol)
        psi4.optimize('scf/STO-3G', molecule=mol)
        #psi4.optimize('scf/6-311pg_2df_2pd_', molecule=mol, symmetry='c1' )
        #psi4.optimize('scf/cc-pvdz', molecule=mol, engine='geometric' )
        #psi4.optimize('pbe/cc-pvdz', molecule=mol)
        #psi4.optimize('b3lyp/cc-pvdz', molecule=mol)
        #psi4.optimize('mp2/cc-pvdz', molecule=mol)
        #mol.save_xyz_file('relaxed/'+fname,1)

    # ----- save output
    geom_lines = mol.save_string_xyz().split('\n')[1:]

    # ------------ Call for first stage fit
    Qs = resp.resp([mol], resp_options)
    Q_esp  = Qs[0] 
    Q_resp = Qs[1]
    print('Electrostatic Potential Charges\n',            Q_esp)
    print('Restrained Electrostatic Potential Charges\n', Q_resp)

    save_xyz_Q( "relaxed/"+name+".xyz", geom_lines, Q_resp )

# ======== Main

#psi4.core.be_quiet()



for name in names:
    psi4.core.set_output_file(name+'.log', False)
    process_molecule( name, bRelax=True )


'''
# Change the value of the RESP parameter A
options['RESP_A'] = 0.001

# ----------- Add constraint for atoms fixed in second stage fit
constraint_charge = []
for i in range(4, 8):
    constraint_charge.append([charges1[1][i], [i+1]])
options['constraint_charge'] = constraint_charge
options['constraint_group'] = [[2, 3, 4]]
options['grid'] = ['1_%s_grid.dat' %mol.name()]
options['esp'] = ['1_%s_grid_esp.dat' %mol.name()]

# Call for second stage fit
charges2 = resp.resp([mol], options)

# Get RESP charges
print("\nStage Two:\n")
print('RESP Charges')
print(charges2[1])
'''
