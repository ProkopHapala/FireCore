
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
    "GradElem": 1E-4
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
            Convergence {{ 
                GradElem = {params["GradElem"]}
        """) )
        if 'DispElem' in params: 
            hsd.write('        DispElem = %e \n' %params['DispElem']  ) 
        if 'EConv' in params: 
            hsd.write('        Energy   = %e \n' %params['EConv']    ) 
        hsd.write("    }\n")
        hsd.write("}\n")

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
        
        hsd.write("    MaxAngularMomentum {\n")
        for ename in enameset:  hsd.write(f'        {ename} = "{elements.ELEMENT_DICT[ename][4]}" \n'   )
        hsd.write("    }\n")

        if method=="D3H5":

        
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

        if 'SCCTolerance' in params: 
            hsd.write('    SCCTolerance = %e \n' %params['SCCTolerance']    )
        if 'MaxSccIterations' in params:
            hsd.write('    MaxSccIterations = %i \n' %params['MaxSccIterations']    )
        if 'Mixer' in params:
            hsd.write('    Mixer = %s \n' %params['Mixer']    ) 
        if 'Temperature' in params:
            hsd.write('    Filling = Fermi {Temperature [K] = %f }\n' %params['Temperature']    )

        hsd.write("}\n")

    return
    #Options { WriteDetailedOut = No }
    #Analysis { CalculateForces = Yes }
    #ParserOptions { ParserVersion = 10 }
    #Parallel { UseOmpThreads = Yes }

def makeDFTBjob_pbc( enames, apos, lvs, fname='dftb_in.hsd', basis_path='/home/prokop/SIMULATIONS/dftbplus/slakos/3ob-3-1/', 
                     nk=(1,1,1), k_shift=(0.5,0.0,0.0), opt=False, params=default_params, SCCTolerance=1e-5, MaxScc=200, Temperature=300, MixingParameter=0.2, fixed_atoms=None ):
    """Write a DFTB+ input for a periodic calculation using GenFormat (supercell type S).
    
    Args:
        enames: list of element names per atom
        apos:   (natoms,3) array of Cartesian atomic positions [Angstrom]
        lvs:    (3,3) lattice vectors (rows are a1, a2, a3) [Angstrom]
        fname:  output HSD filename
        basis_path: path to Slater-Koster files
        nk:     (3,) k-point folding along a1, a2, a3
        k_shift: (3,) k-point shift (0.5 for Monkhorst-Pack half-shift)
        opt:    if True, add geometry optimization driver
        params: dict with optimizer settings (uses default_params keys)
        SCCTolerance: SCC convergence threshold
        MaxScc: max SCC iterations
        Temperature: electronic temperature [K]
        MixingParameter: Broyden mixing parameter (0.0-1.0)
        fixed_atoms: list of 0-based atom indices to fix during relaxation (adds Constraint block)
    """
    enameset = sorted(set(enames))
    ename_to_idx = {e: i+1 for i, e in enumerate(enameset)}  # GenFormat: 1-indexed
    natoms = len(enames)
    lvs = np.array(lvs)

    with open(fname, 'w') as hsd:
        # GenFormat geometry block
        hsd.write('Geometry = GenFormat {\n')
        hsd.write(f'  {natoms}  S\n')
        hsd.write('  ' + ' '.join(enameset) + '\n')
        for i, (ename, pos) in enumerate(zip(enames, apos)):
            idx = ename_to_idx[ename]
            hsd.write(f'  {i+1} {idx}   {pos[0]:.10f}   {pos[1]:.10f}   {pos[2]:.10f}\n')
        # Origin + lattice vectors
        hsd.write('  0.000000000  0.000000000  0.000000000\n')
        for row in lvs:
            hsd.write(f'  {row[0]:.10f}  {row[1]:.10f}  {row[2]:.10f}\n')
        hsd.write('}\n\n')

        # Geometry optimization driver
        if opt:
            # MovedAtoms: all atoms except fixed ones
            if fixed_atoms:
                fixed_1based = sorted([i+1 for i in fixed_atoms])  # 1-based for DFTB+
                # Build MovedAtoms as range excluding fixed
                all_idx = set(range(1, natoms+1))
                moved_idx = sorted(all_idx - set(fixed_1based))
                if moved_idx:
                    # Compact representation: list of ranges
                    moved_str = ' '.join(str(i) for i in moved_idx)
                else:
                    moved_str = "1:-1"  # fallback
            else:
                moved_str = "1:-1"
            hsd.write(dedent(f"""Driver = GeometryOptimization {{
    Optimizer = {params["Optimizer"]}
    MovedAtoms = {moved_str}
    MaxSteps = {params["MaxSteps"]}
    OutputPrefix = "geom.out"
    LatticeOpt = No
    Convergence {{ GradElem = {params["GradElem"]} }}
}}\n\n"""))
        
        # Force calculation (always, needed to monitor constraints)
        hsd.write('\nAnalysis {\n  CalculateForces = Yes\n}\n\n')

        # Hamiltonian
        hsd.write('Hamiltonian = DFTB {\n')
        hsd.write('  Scc = Yes\n')
        hsd.write('  SlaterKosterFiles = Type2FileNames {\n')
        hsd.write(f'    Prefix = {basis_path}\n')
        hsd.write('    Separator = "-"\n')
        hsd.write('    Suffix = ".skf"\n')
        hsd.write('  }\n')
        hsd.write('  MaxAngularMomentum {\n')
        for ename in enameset:
            hsd.write(f'    {ename} = "{elements.ELEMENT_DICT[ename][4]}"\n')
        hsd.write('  }\n')
        # K-points via SupercellFolding (Monkhorst-Pack)
        hsd.write('  KPointsAndWeights = SupercellFolding {\n')
        hsd.write(f'    {nk[0]} 0 0\n')
        hsd.write(f'    0 {nk[1]} 0\n')
        hsd.write(f'    0 0 {nk[2]}\n')
        hsd.write(f'    {k_shift[0]:.1f} {k_shift[1]:.1f} {k_shift[2]:.1f}\n')
        hsd.write('  }\n')
        hsd.write(f'  SCCTolerance = {SCCTolerance:.2e}\n')
        hsd.write(f'  MaxSccIterations = {MaxScc}\n')
        hsd.write(f'  Filling = Fermi {{ Temperature [K] = {Temperature} }}\n')
        hsd.write('  Mixer = Broyden {\n')
        hsd.write(f'    MixingParameter = {MixingParameter}\n')
        hsd.write('  }\n')
        hsd.write('}\n')


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







