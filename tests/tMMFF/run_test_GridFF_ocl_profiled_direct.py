#!/usr/bin/env python3
import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time
import json

# Dictionary to store kernel execution times
kernel_times = {}

# Original imports
sys.path.append("../../")

from pyBall.tests import ocl_GridFF_new as gff
from pyBall.OCL.GridFF import GridFF_cl

# Save original methods
original_methods = {}

# Monkey patch the GridFF_cl class to add profiling
for name in dir(GridFF_cl):
    # Skip special methods and non-callable attributes
    if name.startswith('__') or not callable(getattr(GridFF_cl, name)):
        continue
    
    # Save the original method
    original_methods[name] = getattr(GridFF_cl, name)
    
    # Create a profiling wrapper
    def create_wrapper(method_name, original_method):
        def wrapper(self, *args, **kwargs):
            # Record start time
            start_time = time.time()
            
            # Call the original method
            result = original_method(self, *args, **kwargs)
            
            # Record end time
            end_time = time.time()
            
            # Calculate duration
            duration = (end_time - start_time) * 1000  # Convert to ms
            
            # Add to kernel times
            if method_name not in kernel_times:
                kernel_times[method_name] = []
            kernel_times[method_name].append(duration)
            
            return result
        return wrapper
    
    # Replace the method with the wrapper
    setattr(GridFF_cl, name, create_wrapper(name, original_methods[name]))

def test_gridFF_ocl():
    print("py======= test_gridFF_ocl() START")
    
    # ========= Setup
    
    bSymetrize = False
    print("bSymetrize ", bSymetrize)
    
    # ========= Atoms
    
    atoms = np.array([
        [2.0, 2.0, 0.0],
        [0.0, 0.0, 0.0]
    ])
    
    # ========= Atom Types
    
    atypes = np.array([17, 11], dtype=np.int32)
    print("Raw atom types:", atypes)
    
    # ========= Atom Charges
    
    charges = np.array([-0.7, 0.7])
    print("Raw atom charges:", charges)
    
    # Check total charge
    Qtot = np.sum(charges)
    Qabs = np.sum(np.abs(charges))
    print("Qtot ", Qtot, " Qabs ", Qabs)
    
    # ========= Van der Waals Parameters
    
    # Get parameters from ElementTypes.dat
    try_load_mmff()
    REVdW = mmff.getREVdW(atypes)
    print("Raw REvdW parameters from ElementTypes.dat:", REVdW)
    
    # ========= Combine Parameters
    
    # Combine R, E, Q parameters
    REQs = np.zeros((len(atoms), 4))
    REQs[:, 0] = REVdW[:, 0]  # R
    REQs[:, 1] = REVdW[:, 1]  # E
    REQs[:, 2] = charges      # Q
    
    # Print component-wise breakdown
    print("\nComponent-wise breakdown:")
    for i, atom in enumerate(atoms):
        print(f"Atom[{i}]:")
        print(f"  Type: {atypes[i]}")
        print(f"  Position: {atom}")
        print(f"  Charge: {charges[i]}")
        print(f"  REvdW raw: {REVdW[i]}")
        print(f"  Final REQ: {REQs[i]}")
    
    # ========= Grid Setup
    
    # Set up grid
    print("Atoms:", atoms)
    print("Debug :: REQs =", REQs)
    
    # Copy atoms to avoid modifying the original
    new_atoms = atoms.copy()
    print("New_Atoms:", new_atoms)
    
    # Set z0
    z0 = 0.0
    print("test_gridFF_ocl() z0= ", z0)
    
    # ========= Grid Parameters
    
    # Set grid parameters
    ng = [40, 40, 200]
    ng[0] = gff.next_nice(ng[0])
    ng[1] = gff.next_nice(ng[1])
    ng[2] = gff.next_nice(ng[2])
    
    # ========= Coulomb Potential
    
    print("!!!! Starting Coulomb potential calculation...")
    
    # Calculate Coulomb potential
    clgff = gff.clgff
    Vcoul = clgff.makeCoulombEwald_slab(new_atoms, REQs, ng, z0=z0, Lz_slab=20.0, bTranspose=True, bSaveQgrid=True)
    
    # ========= Morse Potential
    
    print("Starting Morse potential calculation...")
    
    # Calculate Morse potential
    lvec = np.array([[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 20.0]])
    nPBC = gff.autoPBC(lvec, Rcut=20.0, mask=(1, 1, 0))
    print("autoPBC(nPBC_mors): ", nPBC)
    
    # Set grid parameters for Morse potential
    ng = [40, 40, 200]
    ng[0] = gff.next_nice(ng[0])
    ng[1] = gff.next_nice(ng[1])
    ng[2] = gff.next_nice(ng[2])
    
    # Calculate Morse potential
    VPaul, VLond = clgff.make_MorseFF(new_atoms, REQs, ng, z0=z0, nPBC=nPBC)
    
    # ========= Save Results
    
    # Create output directory
    path = "./data/NaCl_1x1_L1"
    os.makedirs(path, exist_ok=True)
    print("test_gridFF_ocl() path = ", path)
    
    # Print min/max values
    print("Paul  min,max = ", VPaul.min(), VPaul.max())
    print("Lond  min,max = ", VLond.min(), VLond.max())
    print("Coul  min,max = ", Vcoul.min(), Vcoul.max())
    print("CoulB min,max = ", Vcoul.min(), Vcoul.max())
    
    # Save results
    print("test_gridFF_ocl() path = ", path)
    
    # Print final stats
    print("Final PLQ Coulomb layer stats:", Vcoul.min(), Vcoul.max(), np.sum(Vcoul))
    
    # Save to file
    full_name = path + "/Bspline_PLQd.npy"
    print("test_gridFF_ocl() - save Morse to: ", full_name)
    np.save(full_name, {"Paul": VPaul, "Lond": VLond, "Coul": Vcoul})
    
    print("py======= test_gridFF() DONE")

def try_load_mmff():
    global mmff
    if mmff is None:
        from pyBall import MMFF as mmff_
        mmff = mmff_

# Global variables
mmff = None

# Run the test
test_gridFF_ocl()

# Print profiling results
print("\n===== OpenCL Profiling Results =====")
print(f"Number of Methods: {len(kernel_times)}")

print("\nMethod Executions:")
for i, (name, times) in enumerate(sorted(kernel_times.items(), key=lambda x: sum(x[1]), reverse=True)):
    total_time = sum(times)
    avg_time = total_time / len(times)
    print(f"{i+1}. {name}")
    print(f"   Total Time: {total_time:.3f} ms")
    print(f"   Average Time: {avg_time:.3f} ms")
    print(f"   Executions: {len(times)}")

# Save profiling results to a file
output_dir = os.path.join(os.getcwd(), "firecore_profile_results")
os.makedirs(output_dir, exist_ok=True)

results = {
    "method_count": len(kernel_times),
    "methods": [{"name": name, "times": times, "total_time": sum(times), "avg_time": sum(times) / len(times), "executions": len(times)} for name, times in kernel_times.items()]
}

results_file = os.path.join(output_dir, "firecore_profiling_results.json")
with open(results_file, 'w') as f:
    json.dump(results, f, indent=2)

print(f"\nProfiling results saved to {results_file}")

# Plot method durations
# Extract method names and total times
names = [k for k in kernel_times.keys()]
total_times = [sum(kernel_times[k]) for k in names]

# Sort by total time
sorted_data = sorted(zip(names, total_times), key=lambda x: x[1], reverse=True)
names = [x[0] for x in sorted_data]
total_times = [x[1] for x in sorted_data]

# Create a bar chart
plt.figure(figsize=(12, 6))
bars = plt.bar(range(len(names)), total_times)

# Add labels and title
plt.xlabel('Method')
plt.ylabel('Total Duration (ms)')
plt.title('OpenCL Method Execution Times')
plt.xticks(range(len(names)), names, rotation=45, ha='right')

# Add duration values on top of bars
for bar in bars:
    height = bar.get_height()
    plt.text(bar.get_x() + bar.get_width()/2., height + 0.1,
            f'{height:.2f}', ha='center', va='bottom')

plt.tight_layout()

# Save to file
plot_file = os.path.join(output_dir, "firecore_method_durations.png")
plt.savefig(plot_file)
print(f"Method duration plot saved to {plot_file}")

# Restore original methods
for name, method in original_methods.items():
    setattr(GridFF_cl, name, method)
