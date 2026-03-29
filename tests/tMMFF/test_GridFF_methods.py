import os
import numpy as np
import matplotlib.pyplot as plt
import time
from pyBall.OCL.GridFF import GridFF_cl, GridShape
from pyBall import atomicUtils as au
from pyBall.AtomicSystem import AtomicSystem

def make_atoms_arrays( atoms=None, fname=None, bSymetrize=False, Element_Types_name="./data/ElementTypes.dat", bSqrtEvdw=True ): 
    if atoms is None:
        atoms = AtomicSystem( fname=fname )
    if bSymetrize:
        atoms, ws = atoms.symmetrized()
    REvdW = au.getVdWparams( atoms.atypes, fname=Element_Types_name )
    na = len(atoms.atypes)
    REQs=np.zeros( (na,4), dtype=np.float32 )
    xyzq=np.zeros( (na,4), dtype=np.float32 )
    xyzq[:,:3] = atoms.apos
    xyzq[:,3]  = atoms.qs
    REQs[:,0]  = REvdW[:,0]
    REQs[:,1]  = np.sqrt(REvdW[:,1]) if bSqrtEvdw else REvdW[:,1]
    REQs[:,2]  = atoms.qs
    REQs[:,3]  = 0.0
    return xyzq, REQs, atoms

def test_methods_convergence(fname="./data/xyz/NaCl_1x1_L1.xyz", Element_Types_name="../tMMFF/data/ElementTypes.dat"):
    print("Starting GridFF methods convergence test...")
    
    # 1. Setup atoms and grid
    xyzq, REQs, atoms = make_atoms_arrays(fname=fname, Element_Types_name=Element_Types_name)
    grid = GridShape(dg=(0.1, 0.1, 0.1), lvec=atoms.lvec)
    
    clgff = GridFF_cl()
    clgff.set_grid(grid)
    
    # 2. Generate Morse Potential (Reference for fitting)
    z0 = xyzq[:, 2].max()
    g0 = (-grid.Ls[0] * 0.5, -grid.Ls[1] * 0.5, z0)
    nPBC = (4, 4, 0) # Minimal PBC for fast test
    
    print("Generating Morse potential...")
    clgff.make_MorseFF(xyzq, REQs, nPBC=nPBC, lvec=atoms.lvec, g0=g0, GFFParams=(0.1, 1.5, 0.0, 0.0), bReturn=False)
    
    # 3. Test different methods
    methods = [
        {'name': 'MD',     'dt': 0.5,  'damp': 0.15},
        {'name': 'CG',     'dt': 0.5,  'damp': 0.0 }, # CG doesn't use dt/damp in the same way
        {'name': 'Jacobi', 'dt': 0.2,  'damp': 0.0 },
    ]
    
    results = {}
    
    for m in methods:
        print(f"\nTesting method: {m['name']}")
        # We use V_Paul_buff as it's already populated by make_MorseFF
        V_fit, trj = clgff.fit3D(clgff.V_Paul_buff, 
                                 nmaxiter=500, 
                                 dt=m['dt'], 
                                 damp=m['damp'], 
                                 method=m['name'], 
                                 bConvTrj=True, 
                                 bPrint=True)
        results[m['name']] = trj

    # 4. Plot results
    plt.figure(figsize=(10, 6))
    for name, trj in results.items():
        if trj is not None:
            plt.plot(trj[:, 0], trj[:, 1], label=f"{name} |F|")
            # plt.plot(trj[:, 0], trj[:, 2], '--', label=f"{name} |E|")

    plt.yscale('log')
    plt.xlabel('Iteration')
    plt.ylabel('Error')
    plt.title('Convergence Comparison: Morse Potential B-spline Fitting')
    plt.legend()
    plt.grid(True, which="both", ls="-", alpha=0.5)
    plt.savefig("fitting_convergence.png")
    print("\nResults saved to fitting_convergence.png")
    plt.show()

if __name__ == "__main__":
    test_methods_convergence()
