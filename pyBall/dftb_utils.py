
import sys
import os
import numpy as np
from . import elements
from . import atomicUtils as au
from textwrap import dedent,indent

methods=[  'GFN2' ]

methods_XTB  = { 'GFN1', 'GFN2', 'IPEA1' }
methods_dftb = { '3ob', 'D3H5' }


H5Scaling ={
    "O" : 0.06,
    "N" : 0.18,
    "S" : 0.21,
}

default_params={
    "RScaling" : 0.714,
    "WScaling" : 0.25,
    "sr6"      : 1.25,
    "alpha6"   : 29.61,
    "s6"       : 1.0,
    "s8"       : 0.49,
    "HHRepulsion" : "Yes",
    #"Optimizer"  : "Rational{}",
    "Optimizer"  : "LBFGS{  Memory = 20 }",
    #"Optimizer"  : "FIRE{StepSize = 1.0}",       #
    "MaxSteps": 1000,
    "GradAMax": 1E-4
    #'Temperature' : 300
}

# ============ Setup

def makeDFTBjob( enames=None, fname='dftb_in.hsd', gname="input.xyz", method='D3H5', cell=None, basis_path='/home/prokop/SIMULATIONS/dftbplus/slakos/3ob-3-1/', params=default_params, opt=True ):
    enameset = set( enames )
    #print( "enameset = ", enameset )
    hsd = open(fname,'w')
    hsd.write(dedent("""
    Geometry = xyzFormat {
        <<< "%s"
    }\n""" %gname ))

    if opt:
        hsd.write(dedent(f"""  
        Driver = GeometryOptimization {{
            Optimizer = {params["Optimizer"]}
            MovedAtoms = 1:-1
            MaxSteps = {params["MaxSteps"]}
            OutputPrefix = "geom.out"
            Convergence {{ GradAMax = {params["GradAMax"]} }}
        }}
        """))

    if method in methods_XTB:
        hsd.write(dedent("""
        Hamiltonian = xTB {
            Method = "%s-xTB"
        }
        """ %method ))
    elif method in methods_dftb:
        hsd.write(dedent("""
        Hamiltonian = DFTB {
            Scc = Yes
            SlaterKosterFiles = Type2FileNames {
                Prefix = %s
                Separator = "-"
                Suffix = ".skf"
            }
        """ %basis_path ))
    
        if method=="D3H5":
            hsd.write("    MaxAngularMomentum {\n")
            for ename in enameset:  hsd.write(f'        {ename} = "{elements.ELEMENT_DICT[ename][4]}" \n'   )
            hsd.write("    }")
        
            hsd.write(indent(dedent(f"""          
            HCorrection = H5 {{
                RScaling = {params["RScaling"]}
                WScaling = {params["WScaling"]}
                H5Scaling {{\n""" ), "    "))
            for ename in enameset: 
                if ename in H5Scaling: hsd.write(f'            {ename} = {H5Scaling[ename]} \n' )
            hsd.write("       }\n")
            hsd.write("    }\n")

            hsd.write(indent(dedent(f"""   
            Dispersion = DftD3 {{
                Damping = ZeroDamping {{
                    sr6    = {params["sr6"]}
                    alpha6 = {params["alpha6"]}
                }}
                s6 = {params["s6"]}
                s8 = {params["s8"]}
                HHRepulsion = {params["HHRepulsion"]}
            }}\n"""  ), "    "))

        if 'Temperature' in params:
            hsd.write('    Filling = Fermi {Temperature [K] = %f }\n' %params['Temperature']    )

        hsd.write("}\n")

    return
    #Options { WriteDetailedOut = No }
    #Analysis { CalculateForces = Yes }
    #ParserOptions { ParserVersion = 10 }
    #Parallel { UseOmpThreads = Yes }

def run( geom=None, params=None, id=0 ):
    idstr = "%03i" %id 
    print( idstr )
    if params['own_dir']:
        cwd = os.getcwd()
        os.mkdir( idstr )
        os.chdir( idstr )
    #try:
    #    os.system( 'cp ../%03i/charges.bin .' %(id-1) )
    #except: pass
    if( id>0 ):
        os.system( 'cp ../%03i/charges.bin .' %(id-1) )
    apos,es = geom
    au.saveXYZ( es=es, xyzs=apos, fname="input.xyz" )
    makeDFTBjob( enames=es, fname='dftb_in.hsd', gname="input.xyz", method=params['method'], cell=params['cell'], basis_path=params['basis'], params=params, opt=params['opt'] )
    os.system('dftb+ > OUT' )
    #os.system( 'grep "Total Energy" OUT | tail -1 | cut -b 52-70' )
    Estr = os.popen('grep "Total Energy" OUT | tail -1 | cut -b 52-70').read()
    E = float(Estr)
    if params['own_dir']:
        os.chdir( cwd )
    return E







