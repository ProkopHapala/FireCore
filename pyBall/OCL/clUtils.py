import sys
import os
import numpy as np
import ctypes
import pyopencl as cl
# import pyopencl.array as cl_array
# import pyopencl.cltypes as cltypes
# import matplotlib.pyplot as plt
import time

bytePerFloat = 4

FFT = None


import pyopencl as cl

def get_nvidia_device( what="nvidia"):
    platforms = cl.get_platforms()
    for platform in platforms:
        devices = platform.get_devices()
        for device in devices:
            if what in device.name.lower():
                # Create the OpenCL context and command queue
                ctx = cl.Context([device])
                queue = cl.CommandQueue(ctx)
                # Print information about the selected device
                print(f"Selected device: {device.name}")
                get_cl_info(device)
                return ctx, queue
    # If no NVIDIA device is found, return None
    print("pyOepnCL error: No {what} device found.")
    return None, None


def try_load_clFFT():
    global FFT
    if FFT is None:
        from gpyfft.fft import FFT as FFT_
        FFT = FFT_

def make_inds_pbc(n):
    return np.array([
    [ 0, 1,   2,   3   ],
    [ 0, 1,   2,   3-n ],
    [ 0, 1,   2-n, 3-n ],
    [ 0, 1-n, 2-n, 3-n ]], dtype=np.int32 );

def roundup_global_size(global_size, local_size):
    remainder = global_size % local_size
    if remainder == 0: return global_size
    return global_size + (local_size - remainder)

def roundup_global_size_3d( global_size,  local_size):
    return (
         roundup_global_size(global_size[0], local_size[0]),
         roundup_global_size(global_size[1], local_size[1]),
         roundup_global_size(global_size[2], local_size[2])
     )

def get_cl_info( device ):
    '''
    also use:
        nvidia-smi -q
    '''
    # Example for NVIDIA GPUs:
    # Pascal: 128 CUDA cores per SM
    # Turing: 64 CUDA cores per SM
    # Ampere: 128 CUDA cores per SM  10496 cores / 82 compute units

    print(f"Device Name: {device.name}")
    print(f"Max Compute Units: {device.max_compute_units}")
    print(f"Max Work Group Size: {device.max_work_group_size}")
    print(f"Global Memory Size: {device.global_mem_size / (1024*1024)} MB")
    print(f"Local Memory Size: {device.local_mem_size / 1024} KB")
    print(f"Max Clock Frequency: {device.max_clock_frequency} MHz")
              
    # Get local memory characteristics
    granularity      = cl.characterize.local_memory_access_granularity(device)
    bank_count       = cl.characterize.local_memory_bank_count(device)
    usable_local_mem = cl.characterize.usable_local_mem_size(device)
    # Print results
    print(f"Local Memory Access Granularity: {granularity} bytes")
    print(f"Number of Local Memory Banks: {bank_count}")
    print(f"Usable Local Memory Size: {usable_local_mem} bytes")

    # Retrieve various characteristics
    fast_math_options = cl.characterize.get_fast_inaccurate_build_options(device)
    simd_group_size   = cl.characterize.get_simd_group_size(device, 4)  # Assuming float4
    double_support    = cl.characterize.has_amd_double_support(device)
    #src_cache_support = cl.characterize.has_src_build_cache(device)

    # Print results
    print(f"Fast Math Options: {fast_math_options}")
    print(f"SIMD Group Size: {simd_group_size}")
    print(f"Double Precision Support: {double_support}")
    #print(f"Source Build Cache Support: {src_cache_support}")

def local_memory_per_workgroup( device, local_size=32, sp_per_cu=128 ):
    return device.local_mem_size /( sp_per_cu/local_size )

def local_memory_float_per_workgroup( device, local_size=32, sp_per_cu=128 ):
    return local_memory_per_workgroup( device, local_size=local_size, sp_per_cu=sp_per_cu )/4

def local_memory_float4_per_workgroup( device, local_size=32, sp_per_cu=128 ):
    return local_memory_per_workgroup( device, local_size=local_size, sp_per_cu=sp_per_cu )/(4*4)


################################################################################

def is_nice(n, allowed_factors={2, 3, 5}):
    """
    Checks if integer n factors completely into allowed_factors.
    For example, 60 (2²×3×5) is nice but 302 (2×151) is not.
    """
    temp = n
    for p in sorted(allowed_factors):
        while temp % p == 0:
            temp //= p
    return temp == 1

def next_nice(n, allowed_factors={2, 3, 5}):
    """
    Returns the smallest integer greater than or equal to n that is FFT-friendly.
    With desired voxel 0.1, for example, if n computes to 302 (2×151), then
    next_nice will bump it upward (to 320, since 320 = 2⁶×5 is acceptable).
    """
    original = n
    while not is_nice(n, allowed_factors):
        n += 1
    print(f"next_nice: for raw value {original} -> adjusted to {n}")
    return n

def adjust_dimensions(dim_vector, desired_voxel=0.1, allowed_factors={2, 3, 5}):
    """
    For each physical dimension L in dim_vector, compute:
         N_target = ceil(L / desired_voxel)
    then bump N_target to the nearest FFT-friendly integer.
    Returns a tuple (ns, new_voxels), where:
         ns         : list of grid sizes for each dimension,
         new_voxels : list of recalculated voxel sizes L/N.
    """
    ns = []
    new_voxels = []
    for L in dim_vector:
        N_target = int(np.ceil(L / desired_voxel))
        N = next_nice(N_target, allowed_factors)
        ns.append(N)
        new_voxels.append(L / N)
    return ns, new_voxels

class GridShape:
    def __init__(self, ns=None, dg=(0.0, 0.0, 0.0), lvec=None, Ls=None,
                 g0=(0.0, 0.0, 0.0), desired_voxel=0.1,
                 allowed_factors={2, 3, 5}):
        # Store lattice vectors so that downstream routines have access.
        self.lvec = lvec
        # If Ls isn’t provided, assume the lengths are given on the diagonal of lvec.
        if Ls is None:
            Ls = (lvec[0][0], lvec[1][1], lvec[2][2])
        # Either use provided ns or calculate them from Ls and desired_voxel.
        if ns is None:
            ns, dg_new = adjust_dimensions(Ls, desired_voxel, allowed_factors)
            dg = tuple(dg_new)
        self.ns = ns
        self.nxyz = np.prod(ns)
        self.Ls = Ls
        self.g0 = g0
        self.dg = dg
        self.dV = dg[0] * dg[1] * dg[2]
        self.V = self.dV * self.nxyz
    def __str__(self):
        return f"GridShape(ns={self.ns}, Ls={self.Ls}, g0={self.g0}, dg={self.dg})"

################################################################################################################



# # --- Add the helper functions here ---
# def is_nice(n, allowed_factors={2,3,5,7,11,13,17,19}):
#     """Checks if n factors completely into allowed_factors."""
#     temp = n
#     for p in sorted(allowed_factors):
#         while temp % p == 0:
#             temp //= p
#     return temp == 1

# def next_nice(n, allowed_factors={2,3,5,7,11,13,17,19}):
#     """Returns the smallest integer greater than or equal to n that is FFT-friendly."""
#     while not is_nice(n, allowed_factors):
#         n += 1
#     return n

# def adjust_dimensions(dim_vector, desired_voxel=0.17, allowed_factors={2,3,5,7,11,13,17,19}):
#     """
#     For each physical dimension L in dim_vector, computes the number of grid points as:
#          N_target = ceil(L/desired_voxel)
#     and then bumps N_target to the next FFT-friendly integer.
#     Returns the new dimensions as (ns, new_voxels) where:
#          ns   : a list of grid-point counts for each dimension,
#          new_voxels : a list of recalculated voxel sizes L/N.
#     """
#     ns = []
#     new_voxels = []
#     for L in dim_vector:
#         N_target = int(np.ceil(L / desired_voxel))
#         N = next_nice(N_target, allowed_factors)
#         ns.append(N)
#         new_voxels.append(L / N)
#     return ns, new_voxels


# class GridShape:
#     def __init__(self, ns=None, dg=(0.0, 0.0, 0.0), lvec=None, Ls=None, g0=(0.0,0.0,0.0), desired_voxel=0.17, allowed_factors={2,3,5,7,11,13,17,19}):
#         # Store lvec so that downstream routines (e.g. GridCL) have access.
#         self.lvec = lvec  
#         # If Ls isn’t provided, assume they come from the diagonal of lvec.
#         if Ls is None:
#             Ls = (lvec[0][0], lvec[1][1], lvec[2][2])
#         # Either use the provided ns or calculate them from Ls and desired_voxel.
#         if ns is None:
#             ns, dg_new = adjust_dimensions(Ls, desired_voxel, allowed_factors)
#             dg = tuple(dg_new)
#         self.ns = ns
#         self.nxyz = np.prod(ns)
#         self.Ls = Ls
#         self.g0 = g0
#         self.dg = dg
#         self.dV = dg[0] * dg[1] * dg[2]
#         self.V = self.dV * self.nxyz
#     def __str__(self):
#         return f"GridShape(ns={self.ns}, Ls={self.Ls}, g0={self.g0}, dg={self.dg})"
    
                                                           
# class GridShape:
#     def __init__(self, ns=None, dg=(0.0, 0.0, 0.0), lvec=None, Ls=None, g0=(0.0,0.0,0.0), desired_voxel=0.15):
#         self.lvec = lvec
#         # Use lvec to determine cell dimensions if Ls is not given
#         if Ls is None:
#             Ls = (lvec[0][0], lvec[1][1], lvec[2][2])
#         if ns is None:
#             ns = []
#             dg_new = []
#             for L in Ls:
#                 N, new_voxel = adjust_grid_size(L, desired_voxel)
#                 print("L =", L, " → N =", N, " new_voxel =", new_voxel)
#                 ns.append(N)
#                 dg_new.append(new_voxel)
#                 print(f"For L = {L:.4f}, chosen grid points: {N} with voxel size: {new_voxel:.4f}")
#             dg = dg_new
#         self.ns = ns
#         self.nxyz = np.prod(ns)
#         self.Ls = Ls
#         self.g0 = g0
#         self.dg = dg
#         self.dV = dg[0] * dg[1] * dg[2]
#         self.V = self.dV * self.nxyz

#     def __str__(self):
#         return f"GridShape(ns={self.ns}, Ls={self.Ls}, g0={self.g0}, dg={self.dg})"

# class GridShape:
#     def __init__(self, ns=None, dg=(0.0,0.0,0.0), lvec=None, Ls=None, g0=(0.0,0.0,0.0) ):
#         if lvec is None: lvec = np.array( [ [Ls[0],0.,0.], [0.,Ls[1],0.], [0.,0.,Ls[2]] ] )
#         if Ls   is None: Ls = (lvec[0][0],lvec[1][1],lvec[2][2])
#         if ns   is None: ns = (int((Ls[0]/dg[0])+0.5),int((Ls[1]/dg[1])+0.5),int((Ls[2]/dg[2])+0.5))
#         self.ns   = ns
#         self.nxyz = np.prod(ns)
#         self.Ls   = Ls
#         self.g0 = g0
#         self.dg   = dg
#         self.lvec = lvec
#         self.dV   = dg[0]*dg[1]*dg[2] 
#         self.V    = self.dV*self.nxyz
        

#     def __str__(self):
#         return f"GridShape(ns={self.ns}, Ls={self.Ls}, pos0={self.pos0}, dg={self.dg}, lvec={self.lvec})"

class GridCL:
    def __init__(self, gsh : GridShape ):
        self.nxyz = np.int32( gsh.ns[0]*gsh.ns[1]*gsh.ns[2] )
        self.ns = np.array( gsh.ns+[self.nxyz], dtype=np.int32   )
        self.g0 = np.array( gsh.g0+(0.,),   dtype=np.float32 )
        self.dg = np.array( gsh.dg+(0.,),   dtype=np.float32 )
        self.a  = np.array( ( *gsh.lvec[0], 0.), dtype=np.float32 )
        self.b  = np.array( ( *gsh.lvec[1], 0.), dtype=np.float32 )
        self.c  = np.array( ( *gsh.lvec[2], 0.), dtype=np.float32 )
        self.da = np.array( ( *(gsh.lvec[0]/gsh.ns[0]), 0.), dtype=np.float32 )
        self.db = np.array( ( *(gsh.lvec[1]/gsh.ns[1]), 0.), dtype=np.float32 )
        self.dc = np.array( ( *(gsh.lvec[2]/gsh.ns[2]), 0.), dtype=np.float32 )
        
    def __str__(self):
        return f"GridShape(ng={self.ng}, g0={self.g0}, dg={self.dg}\n   La={self.La},Lb={self.Lb},Lc={self.Lc}\n    da={self.da},db={self.db},dc={self.dc})"
