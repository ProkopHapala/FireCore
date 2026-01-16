import numpy as np
import os
import sys

# Add the parent directory to the path to allow imports from pyBall
sys.path.append(os.path.join(os.path.dirname(__file__), '../..'))

from pyBall import eFF as eff
from pyBall.OCL import eFF_ocl as ocl

xyz_filename = os.path.normpath(os.path.join(os.path.dirname(__file__), "../../cpp/sketches_SDL/Molecular/data/H2O_fixcore.xyz"))
print("--- Comparing CPU and GPU eFF implementations for H2O_fixcore ---")

# ==========================
#      CPU Calculation
# ==========================
print("\n--- Running CPU Calculation ---")
eff.setVerbosity(4,1)
eff.processXYZ_e(xyz_filename, nstepMax=0)
eff.getBuffs()

# Store CPU forces
na = eff.na
ne = eff.ne
cpu_forces = np.zeros((na + ne, 4))
cpu_forces[:na, :3] = eff.aforce
cpu_forces[na:, :3] = eff.eforce
cpu_forces[na:, 3]  = eff.fsize

# ==========================
#      GPU Calculation
# ==========================
print("\n--- Running GPU Calculation ---")
gpu_forces = None
try:
    eff_gpu = ocl.EFF_OCL()
    eff_gpu.load_xyzs(xyz_filename)
    eff_gpu.realloc_buffers()
    eff_gpu.upload_data()
    eff_gpu.relax_systems(n_steps=1, dt=0.0, damping=0.0)
    gpu_forces = eff_gpu.get_forces()
except Exception as e:
    print(f"Could not run GPU calculation: {e}")

# ==========================
#      Comparison
# ==========================
print("\n--- Comparison ---")
if gpu_forces is not None:
    # The GPU kernel calculates forces for all particles in the workgroup (nloc)
    # We only care about the forces for the actual particles in the system.
    ntot = na + ne
    gpu_forces = gpu_forces[:ntot]

    if cpu_forces.shape != gpu_forces.shape:
        print(f"FAIL: Force array shapes do not match!")
        print(f"CPU shape: {cpu_forces.shape}, GPU shape: {gpu_forces.shape}")
    else:
        diff = cpu_forces - gpu_forces
        max_abs_diff = np.max(np.abs(diff))
        print(f"Max absolute difference: {max_abs_diff}")

        if max_abs_diff < 1e-4:
            print("PASS: Implementations match.")
        else:
            print("FAIL: Implementations do NOT match.")
            print("=========================================")
            print("CPU Forces:")
            print(cpu_forces)
            print("\nGPU Forces:")
            print(gpu_forces)
            print("\nDifference (CPU - GPU):")
            print(diff)
            print("=========================================")
            print("NOTE: The OpenCL kernel `localMD` likely has incorrect force accumulation logic.")
            print("      The test is correctly identifying a discrepancy. => modify `__kernel void localMD` in `/FireCore/cpp/common_resources/cl/eFF.cl`")
else:
    print("SKIPPED: GPU forces not available for comparison.")

