import sys
import os
import time
import psi4
import resp

sys.path.append("../../")
from pyBall import psi4_utils as psi4u
from pyBall  import atomicUtils as au

# ========= Setup

indir="./input/"
#outdir="./output/"

bRelax=True
#bRelax=False

methods=[
#    'scf',
#    'pbe',
#    'b3lyp',
#    'mp2'
]
basises=[
#    'sto-3g',
#    '6-31+G',
#    '6-311+G*',
#    '6-311++G**',
#    '6-311++G(3df,3pd)',
#    'cc-pvdz',
#    'aug-cc-pvtz',
#    'def2-QZVPPD',
]

method_bas_pairs = [
('scf','sto-3g'), ('pbe','cc-pvdz')
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

# ======== Functions

def try_make_dirs( dname ):
    try:
        os.mkdir( dname )
    except:
        pass

# ======== Main

#psi4.core.be_quiet()

#names = [ f.split('.')[0] for f in os.listdir(indir) ]

#names = [ "backbone_pasivated-H", "backbone_pasivated-R" ]
names = [ "pyridine" ]


#names =["hexa_hb3_donor"]
print(names)

# method_bas_pair = []
# for method in methods:
#     for basis in basises:
#         method_bas_pairs.append( (method,basis) )

for method_bas in method_bas_pairs:
    method,basis = method_bas
    outdir = method+"/"+basis+"/"
    try_make_dirs( method )
    try_make_dirs( outdir )
    for name in names:
        print( "# ======= Molecule: ", name, method, basis )
        psi4.core.set_output_file(outdir+name+'.log', False)
        t0 = time.time_ns()
        try:
            psi4u.psi4resp( name, bRelax=bRelax, method=method, basis=basis, indir=indir, outdir=outdir, psi4_options=psi4_options, resp_options=resp_options )
        except Exception as e: 
            print(e)
        t = time.time_ns() - t0; print( "time: ", t*1e-9, "[s]" )