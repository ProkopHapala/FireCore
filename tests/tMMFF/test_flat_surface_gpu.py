
import sys
import os
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))
from pyBall import atomicUtils as au
from pyBall.AtomicSystem import AtomicSystem
from pyBall.OCL.MMFF import MMFF
from pyBall.OCL.MolecularDynamics import MolecularDynamics

def run_test_flat_surface():
    print("Testing Flat Surface Potential on GPU...")

    # 1. Create a simple system: Single atom (Oxygen)
    # We place it at various heights to scan the potential
    na = 1
    apos = np.zeros((na, 3), dtype=np.float64)
    apos[0] = [0.0, 0.0, 3.0] # Initial z
    
    atom_types = { "O": au.elements.ELEMENT_DICT["O"][0] } # Dummy
    # We don't really need full MMFF setup for a single atom if we just use run_getSurfFlat
    # But MolecularDynamics expects an MMFF object for realloc/pack_system usually.
    # Or we can use init_with_atoms.
    
    # Let's use init_with_atoms for simplicity of checking the kernel first
    md = MolecularDynamics()
    
    # Manually set params
    # R, E, Q
    REQs = np.zeros((na, 4), dtype=np.float32)
    REQs[0] = [1.5, 0.1, 0.0, 0.0] # R=1.5, E=0.1
    
    md.init_with_atoms(na=na, atoms=np.hstack([apos, np.zeros((na,1))]), REQs=REQs)
    md.setup_kernels() # Regenerate args including getSurfFlat
    
    # 2. Configure Surface Parameters
    # Plane at z=0, Normal=(0,0,1)
    surf_pos0 = np.array([0.0, 0.0, 0.0, 0.0], dtype=np.float32)
    surf_normal = np.array([0.0, 0.0, 1.0, 0.0], dtype=np.float32)
    
    # Surface REQ: R=1.5, E=0.1
    surf_REQ = np.array([1.5, 0.1, 0.0, 0.0], dtype=np.float32)
    
    # Hamaker LJ9-3 (mode=1), K=1.6 (not used for LJ93 but passed)
    # Param: x=K, y=mode
    surf_param_hamaker = np.array([1.6, 1.0, 0.0, 0.0], dtype=np.float32) 
    
    # Update kernel params
    md.kernel_params['surf_pos0'] = surf_pos0
    md.kernel_params['surf_normal'] = surf_normal
    md.kernel_params['surf_REQ'] = surf_REQ
    
    
    # 3. Scan Z
    zs = np.linspace(1.0, 6.0, 50)
    Es_hamaker = []
    Fs_hamaker = []
    
    print("Scanning Hamaker LJ 9-3...")
    md.kernel_params['surf_param'] = surf_param_hamaker
    
    for z in zs:
        # Update position
        atoms = np.zeros((na, 4), dtype=np.float32)
        atoms[0, :3] = [0.0, 0.0, z]
        md.toGPU('apos', atoms)
        
        # Clear forces (important!)
        # The kernel does fapos[iav] += ... so we must clear or assume 0
        # Actually init_with_atoms clears aforce. 
        # But we reuse buffer. We should clear it.
        md.toGPU('aforce', np.zeros_like(atoms))
        
        md.run_getSurfFlat()
        
        _, res_force = md.download_results()
        # res_force[0] is (fx, fy, fz, E)
        f = res_force[0]
        Fs_hamaker.append(f[2])
        Es_hamaker.append(f[3])
        
    # Analytical Hamaker Check
    # combineREQ: R = Ra+Rs, E = Ea*Es
    R_eff = REQs[0,0] + surf_REQ[0] # 1.5 + 1.5 = 3.0
    E_eff = REQs[0,1] * surf_REQ[1] # 0.1 * 0.1 = 0.01
    
    # U(z) = (E/2) * [ (R/z)^9 - 3*(R/z)^3 ]
    def hamaker_ana(z, R, E):
        ratio = R/z
        r3 = ratio**3
        r9 = r3**3
        return 0.5 * E * (r9 - 3*r3)
    
    Es_ana = [hamaker_ana(z, R_eff, E_eff) for z in zs]
    
    plt.figure(figsize=(10, 5))
    plt.subplot(1,2,1)
    plt.plot(zs, Es_hamaker, 'o', label='GPU Hamaker')
    plt.plot(zs, Es_ana, '-', label='Analytical')
    plt.xlabel('z (A)')
    plt.ylabel('Energy (eV)')
    plt.legend()
    plt.title('Hamaker Potential')
    
    # 4. Scan Morse
    print("Scanning Morse...")
    # Morse (mode=2), K=1.6
    surf_param_morse = np.array([1.6, 2.0, 0.0, 0.0], dtype=np.float32)
    md.kernel_params['surf_param'] = surf_param_morse
    
    Es_morse = []
    
    for z in zs:
        atoms = np.zeros((na, 4), dtype=np.float32)
        atoms[0, :3] = [0.0, 0.0, z]
        md.toGPU('apos', atoms)
        md.toGPU('aforce', np.zeros_like(atoms))
        
        md.run_getSurfFlat()
        
        _, res_force = md.download_results()
        f = res_force[0]
        Es_morse.append(f[3])
        
    # Analytical Morse
    # U(z) = E * [ exp(-K*(z-R)) - 1 ]^2 - E  ?
    # relax_multi_mini.cl implementation:
    # exp_term = exp( -K * (z - REQH.x) )
    # E = REQH.y * ( exp_term*exp_term - 2.0f*exp_term )
    # This is equivalent to E * ( (exp - 1)^2 - 1 )
    
    def morse_ana(z, R, E, K):
        ex = np.exp(-K * (z - R))
        return E * (ex*ex - 2*ex)
        
    Es_morse_ana = [morse_ana(z, R_eff, E_eff, 1.6) for z in zs]
    
    plt.subplot(1,2,2)
    plt.plot(zs, Es_morse, 'o', label='GPU Morse')
    plt.plot(zs, Es_morse_ana, '-', label='Analytical')
    plt.xlabel('z (A)')
    plt.ylabel('Energy (eV)')
    plt.legend()
    plt.title('Morse Potential')
    
    plt.savefig('test_flat_surface_gpu.png')
    print("Saved plot to test_flat_surface_gpu.png")
    
    # Check max error
    err_h = np.max(np.abs(np.array(Es_hamaker) - np.array(Es_ana)))
    err_m = np.max(np.abs(np.array(Es_morse) - np.array(Es_morse_ana)))
    print(f"Max Error Hamaker: {err_h:.2e}")
    print(f"Max Error Morse:   {err_m:.2e}")
    
    # float32 precision gives ~1e-5 for Hamaker (due to (R/z)^9 magnification near z≈R)
    if err_h < 5e-5 and err_m < 1e-5:
        print("PASS")
    else:
        print("FAIL (err_h=%e, err_m=%e)" % (err_h, err_m))

if __name__ == "__main__":
    run_test_flat_surface()
