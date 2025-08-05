#!/usr/bin/env python3
"""
Test script for scanNonBond and scanNonBond2 kernels.

This script demonstrates the improved workflow with explicit file paths
and debugging capabilities for testing non-bonded force fields.
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'pyBall'))
from pyBall.OCL.MolecularDynamics import MolecularDynamics
from pyBall.plotUtils import plot1d

def exp_fe( x, b=1.6):
    y  = np.exp(-b*x)
    dy = -b*y
    return y, dy

def exp_pow_cubic( bMorse=1.6, n=5, x1=3.0, x2=5.0 ):
    # 1. Define the target function g(r) = exp(-b*(r-Ri)/n) and its derivative
    y1, dy1 = exp_fe(x1, bMorse/n)
    y2, dy2 = exp_fe(x2, bMorse/n)
    
    # 3. Set up and solve the 4x4 linear system for the coefficients of p(r)
    A = np.array([
        [   x1**3,   x1**2 ,x1 , 1 ], 
        [ 3*x1**2, 2*x1    ,1  , 0 ],
        [   x2**3,   x2**2 ,x2 , 1 ], 
        [ 3*x2**2, 2*x2    ,1  , 0 ]
    ])
    rhs            = np.array([y1, dy1, y2, dy2])
    coeffs         = np.linalg.solve(A, rhs)
    all_roots      = np.roots(coeffs)
    real_roots     = all_roots[np.isreal(all_roots)].real
    if len(real_roots) == 0: raise ValueError("Could not find a physical cutoff root for the given parameters.")
    Rc = np.max(real_roots)
    print(f"exp_pow_cubic() Rc: {Rc}", real_roots )
    print(f"exp_pow_cubic() coeffs: {coeffs}")
    return Rc, coeffs

def eval_exp_lin( x, b, n=5 ):
    p  = (1-x*b/n)
    dp = -b/n + (x*0.)
    #Rc = n/b
    mask = p<0
    p [mask]=0
    dp[mask]=0
    y  = p**n
    dy = n*p**(n-1)*dp 
    return y, dy

def eval_exp_pow_cubic( x, cs, n=5 ):
    dcs  = cs[:-1] * np.array([3, 2, 1])
    p    = np.polyval(cs , x)
    dp   = np.polyval(dcs, x)
    mask = p<0
    p [mask]=0
    dp[mask]=0
    y  = p**n
    dy = n*p**(n-1)*dp 
    return y, dy

def eval_Morse_pow_cubic(rs, R0, E0, cs, n=5):
    x = rs - R0
    p_base = np.polyval(cs, x)
    dcs = cs[:-1] * np.array([3, 2, 1])
    dp_base_dx = np.polyval(dcs, x)
    mask = p_base < 0
    p_base    [mask] = 0
    dp_base_dx[mask] = 0
    p = p_base**n
    dpdr = n * (p_base**(n-1)) * dp_base_dx
    E =        E0 * p*(p-2.0)
    F = -2.0 * E0 *   (p-1.0) * dpdr
    return E, F

def eval_Morse_exact(rs, R0, E0, b):
    x = rs - R0
    p = np.exp(-b * x)
    E = E0 * p * (p - 2.0)
    F =  2 * b *  E0*p*(p - 1.0)
    return E, F

def test_scan_kernels():
    """Test scanNonBond and scanNonBond2 kernels with improved workflow."""
    
    print("=== Testing scanNonBond and scanNonBond2 kernels ===")
    
    # Import after setting up paths

    
    # Initialize MolecularDynamics
    md = MolecularDynamics(nloc=32)
    
    # Define file paths
    base_path   = os.path.dirname(os.path.abspath(__file__))
    forces_path = os.path.join(base_path, '../../tests/tmp/data/cl/Forces.cl')
    kernel_path = os.path.join(base_path, '../../tests/tmp/data/cl/relax_multi_mini.cl')
    
    # Test parameters
    n  = 100  # number of test points
    na = 50  # number of atoms
    
    # Create test data
    np.random.seed(42)
    
    # Test positions (random points in 3D space)
    rs  = np.linspace(0.0, 10.0, n)
    pos = np.zeros((n, 4), dtype=np.float32)
    pos[:, 0] = rs
    
    # Atom positions
    apos = np.random.randn(na, 4).astype(np.float32)
    apos[:,0] = 0.0
    
    # Test atom parameters (REQH)
    REQH = np.array([3.0, 1.0, 0.0, 0.0], dtype=np.float32)  # R, E, Q, H
    #REQH = np.array([3.0, 1.0, 1.0, 0.0], dtype=np.float32)  # R, E, Q, H
    
    # Atom parameters
    REQs = np.random.rand(na, 4).astype(np.float32)
    REQs[:, 0] *= 2.0  # R values between 0-2
    REQs[:, 1] *= 0.5  # E values between 0-0.5
    REQs[:, 2] = 0.0   # Q values (neutral)
    REQs[:, 3] = 0.0   # H values
    
    # Force field parameters
    ffpar = np.array([0.1, 1.0, 0.0, 0.0], dtype=np.float32)  # R2damp, K, unused, unused
    
    # Parse Forces.cl to extract force functions
    print("\n=== Parsing Forces.cl ===")
    force_defs = md.parse_forces_cl(forces_path)
    print(f"Found {len(force_defs['functions'])} functions and {len(force_defs['macros'])} macros")
    
    #print(f"pos: \n", pos)
    #print(f"apos: \n", apos)
    
    # Test scanNonBond
    print("\n=== Testing scanNonBond ===")
    force = np.zeros((n, 4), dtype=np.float32)
    
    print("\n##################################")
    print("1. scanNonBond with getLJQH      ")
    print("##################################")
    
    R2damp = 1e-4
    bMorse = 1.6
    #Rc, morse_pcub = exp_pow_cubic(bMorse, 5, 3.0, 5.0)
    Rc, morse_pcub = exp_pow_cubic(bMorse, 5, 0.0, 2.0);   Rc+=3.0 


    Rcl5 = 5/bMorse

    # y_ref, dy_ref = exp_fe(rs, bMorse)
    # y_lin, dy_lin = eval_exp_lin      (rs, bMorse,     5 )
    # y_cub, dy_cub = eval_exp_pow_cubic(rs, morse_pcub, 5 )
    # fig, (ax1, ax2) = plot1d(rs, [y_ref, y_lin, y_cub], [-dy_ref, -dy_lin, -dy_cub], labels=["exp", "exp_lin", "exp_pow_cubic"],colors=['k','r','g'], figsize=(8,12))
    # #ax1.set_ylim(-1,1)
    # ax1.axvline( Rcl5, c='r', ls='--')
    # ax1.axvline( Rc,   c='g', ls='--')
    # #plt.show()
    
    morse_pcub_ = morse_pcub[::-1]

    #exit()


    # potentials={
    #     #  name    ffpar             code                             
    #     #"LJ"         :  ( [ R2damp ],          "fij = getLJQH(dp, REQH, ffpar.x);"             ),
    #     #"Morse"      :  ( [ R2damp, bMorse],   "fij = getMorseQH(dp, REQH, ffpar.y, ffpar.x);" ),
    #     "exp_r"      :  ( [ bMorse ],          "fij = exp_r     ( dp, ffpar.x );"               ),
    #     "exp_r_lin4" :  ( [ bMorse ],          "fij = exp_r_lin4( dp, ffpar.x );"               ),
    #     "exp_r_lin8" :  ( [ bMorse ],          "fij = exp_r_lin8( dp, ffpar.x );"               ),
    #     "exp_r_lin16":  ( [ bMorse ],          "fij = exp_r_lin16( dp, ffpar.x );"              ),
    #     "exp_r_cub4" :  ( [ *morse_pcub_, Rc], "fij = exp_r_cub4( dp, ffpar.lo, ffpar.hi.x );"  ),
    # }
  

    # inline float4 getMorse_lin5( float3 dp, float R0, float E0, float b ){
    potentials={
        #  name    ffpar             code                             
        "Morse"      :  ( [ bMorse ],          "fij = getMorse      ( dp, REQH.x, REQH.y, ffpar.x );"  ),
        "Morse_lin5" :  ( [ bMorse ],          "fij = getMorse_lin5 ( dp, REQH.x, REQH.y, ffpar.x );"  ),
        "Morse_lin9" :  ( [ bMorse ],          "fij = getMorse_lin9 ( dp, REQH.x, REQH.y, ffpar.x );"  ),
        "Morse_lin17":  ( [ bMorse ],          "fij = getMorse_lin17( dp, REQH.x, REQH.y, ffpar.x );"  ),
        "Morse_cub5" :  ( [ *morse_pcub_, Rc], "fij = getMorse_cub5 ( dp, REQH.x, REQH.y, ffpar.lo,ffpar.hi.x );"  ),
    }

    names    = []
    energies = []
    forces   = []
    for name,(ffpar, code) in potentials.items():
        print(f"{name}: {ffpar}")
        output_path   = os.path.join(base_path, f"tests/tmp/cl/tmp_{name}.cl")
        substitutions = { "macros": { "GET_FORCE_NONBOND": code } }
        md.preprocess_opencl_source(kernel_path, substitutions, output_path)
        md.load_program(kernel_path=output_path)
        fes = md.scanNonBond( pos=pos, force=force, REQH=REQH, ffpar=ffpar )

        names    .append(name)
        energies .append(fes[:,3])
        forces   .append(fes[:,0])
        
        print(f"   scanNonBond ({name}) completed cl_file={output_path}")
        #print(f"   Result forces: \n", forces)
        print(f"   result forces.shape {fes.shape} range: [{fes.min():.3f}, {fes.max():.3f}]")

    fig, (ax1, ax2) = plot1d( rs, energies, forces, labels=names, colors=['k', 'b','c','g',   'r'], figsize=(10,15) )
    ax1.set_ylim(-1,1)
    ax2.set_ylim(-1,1)
    ax1.axvline( Rcl5, c='r', ls='--')
    ax1.axvline( Rc,   c='g', ls='--')
    plt.show()

    # plt.
    # for name, forces in results.items():
    #     rs = pos[:,0]
    #     plt.plot(rs, forces[:,0],label=name+" F")
    #     plt.plot(rs, forces[:,1],label=name+" E")
    #     plt.ylim(-1,1)
    #     plt.grid(alpha=0.2)
    #     plt.legend()
    # plt.show()


    exit() # KEEP this until it works, then we can proceed
    
    print("\n##################################")
    print("2. scanNonBond with getMorseQH")
    print("##################################")
    
    output_path_morse = os.path.join(base_path, 'tests/tmp/data/cl/scanNonBond_getMorseQH.cl')
    
    substitutions_morse = {
        'macros': {
            'GET_FORCE_NONBOND': 'fij = getMorseQH(pos[i], REQH, ffpar);'
        }
    }
    
    md.load_preprocessed_program(kernel_path, substitutions_morse, output_path_morse)
    
    result_morse = md.scanNonBond(
        pos=pos,
        force=force,
        REQs=np.array([REQH] * n, dtype=np.float32),
        n=n,
        REQH=REQH,
        ffpar=ffpar,
        force_function='getMorseQH'
    )
    
    print(f"   scanNonBond (MorseQH) completed")
    print(f"   Result range: [{result_morse.min():.3f}, {result_morse.max():.3f}]")
    print(f"   Preprocessed file saved to: {output_path_morse}")
    


    print("\n##################################")
    print("3. scanNonBond2 with getLJQH")
    print("##################################")

    force2 = np.zeros((n, 4), dtype=np.float32)
    

    output_path2_ljqh = os.path.join(base_path, 'tests/tmp/data/cl/scanNonBond2_getLJQH.cl')
    
    # Variables match kernel context: p, lPos[j], lPar[j], REQH0, ffpar
    substitutions2_ljqh = {
        'macros': {
            'GET_FORCE_NONBOND': '''
                float3 dp = lPos[j].xyz - p;
                float4 REQH = lPar[j];
                REQH.x += REQH0.x;
                REQH.yzw *= REQH0.yzw;
                fij = getLJQH(dp, REQH, ffpar.x);
            '''
        }
    }
    
    md.load_preprocessed_program(kernel_path, substitutions2_ljqh, output_path2_ljqh)
    
    result2_ljqh = md.scanNonBond2(
        pos=pos,
        force=force2,
        apos=apos,
        REQs=REQs,
        n=n,
        na=na,
        REQH0=REQH,
        ffpar=ffpar,
        force_function='getLJQH'
    )
    
    print(f"   scanNonBond2 (LJQH) completed")
    print(f"   Result shape: {result2_ljqh.shape}")
    print(f"   Result range: [{result2_ljqh.min():.3f}, {result2_ljqh.max():.3f}]")
    print(f"   Preprocessed file saved to: {output_path2_ljqh}")
    

    print("\n##################################")
    print("4. scanNonBond2 with getMorseQH")
    print("##################################")
    
    output_path2_morse = os.path.join(base_path, 'tests/tmp/data/cl/scanNonBond2_getMorseQH.cl')
    
    # Variables match kernel context: p, lPos[j], lPar[j], REQH0, ffpar
    substitutions2_morse = {
        'macros': {
            'GET_FORCE_NONBOND': '''
                float3 dp = lPos[j].xyz - p;
                float4 REQH = lPar[j];
                REQH.x += REQH0.x;
                REQH.yzw *= REQH0.yzw;
                fij = getMorseQH(dp, REQH, ffpar.y, ffpar.x);
            '''
        }
    }
    
    md.load_preprocessed_program(kernel_path, substitutions2_morse, output_path2_morse)
    
    result2_morse = md.scanNonBond2(
        pos=pos,
        force=force2,
        apos=apos,
        REQs=REQs,
        n=n,
        na=na,
        REQH0=REQH,
        ffpar=ffpar,
        force_function='getMorseQH'
    )
    
    print(f"   scanNonBond2 (MorseQH) completed")
    print(f"   Result range: [{result2_morse.min():.3f}, {result2_morse.max():.3f}]")
    print(f"   Preprocessed file saved to: {output_path2_morse}")
    
    # Performance comparison
    print("\n=== Performance Test ===")
    import time
    
    # Time scanNonBond
    start = time.time()
    for _ in range(10):
        md.scanNonBond(
            pos=pos, force=force, REQs=np.array([REQH] * n, dtype=np.float32),
            n=n, REQH=REQH, ffpar=ffpar, force_function='getLJQH'
        )
    time1 = (time.time() - start) / 10
    
    # Time scanNonBond2
    start = time.time()
    for _ in range(10):
        md.scanNonBond2(
            pos=pos, force=force2, apos=apos, REQs=REQs,
            n=n, na=na, REQH0=REQH, ffpar=ffpar, force_function='getLJQH'
        )
    time2 = (time.time() - start) / 10
    
    print(f"scanNonBond average time: {time1*1000:.2f} ms")
    print(f"scanNonBond2 average time: {time2*1000:.2f} ms")
    print(f"Speedup: {time1/time2:.2f}x")
    
    print("\n=== All tests completed successfully! ===")
    print("\nPreprocessed files are available for debugging:")
    print(f"  - {output_path_ljqh}")
    print(f"  - {output_path_morse}")
    print(f"  - {output_path2_ljqh}")
    print(f"  - {output_path2_morse}")

if __name__ == "__main__":
    # run it like this:
    #   python -u -m pyBall.OCL.run_scanNonBond | tee OUT-scanNonBond

    test_scan_kernels()
