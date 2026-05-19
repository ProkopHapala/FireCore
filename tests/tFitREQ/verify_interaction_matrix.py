#!/usr/bin/env python3
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Add project root to path
sys.path.append(str(Path(__file__).parent.parent.parent))
from pyBall.OCL.NonBondFitting import setup_driver

def verify_consistency(drv, samples=None):
    if samples is None:
        samples = range(min(5, drv.n_samples))
    
    print(f"\nVerifying consistency for {len(samples)} samples...")
    
    # 1. Total energy from optimized kernel
    E_total_all = drv.evaluate_energies()
    
    results = []
    for iS in samples:
        E_total = E_total_all[iS]
        
        # 2. Interaction matrix
        M = drv.evaluate_interaction_matrix(iS)
        # Sum components
        E_sum = np.sum(M)
        
        diff = E_sum - E_total
        print(f"Sample {iS:3d}: E_total = {E_total:12.6f} eV | E_sum = {E_sum:12.6f} eV | Diff = {diff:12.6e} eV")
        results.append((E_total, E_sum, diff))
        
    return results

if __name__ == "__main__":
    ROOT = Path(__file__).parent.parent.parent
    
    # Setup driver
    model = 'ENERGY_MorseQ_PAIR'
    atypes = str(ROOT / 'cpp/common_resources/AtomTypes.dat')
    drv = setup_driver(model_name=model, atom_types_file=atypes, verbose=1)
    # Load some data
    data_path = ROOT / 'tests/tFitREQ/wb97m_input/'
    xyz_path = data_path / "HF-A1_HF-D1.xyz"
    
    if not xyz_path.exists():
        print(f"File not found: {xyz_path}")
        sys.exit(1)
        
    drv.load_data(str(xyz_path))
    drv.type_set = [('H' ,'H',0.85),('F' ,'H',-0.80)]
    drv.alphaMorse = 1.8
    drv.init_and_upload_energy_only()
    
    # Compile all kernels at once
    drv.compile_all_with_model(model, bPrint=True)
    
    # Verify
    results = verify_consistency(drv)
    
    # Check max diff
    max_diff = max(abs(r[2]) for r in results)
    if max_diff < 1e-5:
        print("\nSUCCESS: Interaction matrix sum matches total energy!")
    else:
        print(f"\nFAILURE: Large discrepancy found: {max_diff:e}")
        sys.exit(1)
