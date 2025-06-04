import os
import sys
import time
import numpy as np
import pyopencl as cl

# Dictionary to store kernel execution times
kernel_times = {}

def profile_kernel(prg, name, queue, global_size, local_size, *args, **kwargs):
    """Profile a kernel execution"""
    # Record start time
    start_time = time.time()
    
    # Call the kernel
    result = getattr(prg, name)(queue, global_size, local_size, *args, **kwargs)
    
    # Record end time
    end_time = time.time()
    
    # Calculate duration
    duration = (end_time - start_time) * 1000  # Convert to ms
    
    # Add to kernel times
    if name not in kernel_times:
        kernel_times[name] = []
    kernel_times[name].append(duration)
    
    return result

# Import the original GridFF module
from .GridFF import GridFF_cl as OriginalGridFF_cl
from .GridFF import GridShape

# Create a profiled version of the GridFF_cl class
class GridFF_cl(OriginalGridFF_cl):
    """Profiled version of the GridFF_cl class"""
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
    
    # Override methods that call OpenCL kernels
    
    def makeCoulombEwald_slab(self, *args, **kwargs):
        """Profiled version of makeCoulombEwald_slab"""
        return super().makeCoulombEwald_slab(*args, **kwargs)
    
    def make_MorseFF(self, *args, **kwargs):
        """Profiled version of make_MorseFF"""
        return super().make_MorseFF(*args, **kwargs)
    
    def fit3D(self, *args, **kwargs):
        """Profiled version of fit3D"""
        return super().fit3D(*args, **kwargs)

# Monkey patch the original GridFF_cl class to use our profiling function
original_enqueue_nd_range_kernel = cl.enqueue_nd_range_kernel

def profiling_enqueue_nd_range_kernel(queue, kernel, global_size, local_size=None, *args, **kwargs):
    """Wrapper for cl.enqueue_nd_range_kernel that profiles the kernel execution"""
    # Get the kernel name
    kernel_name = kernel.function_name if hasattr(kernel, 'function_name') else str(kernel)
    
    # Record start time
    start_time = time.time()
    
    # Call the original function
    event = original_enqueue_nd_range_kernel(queue, kernel, global_size, local_size, *args, **kwargs)
    
    # Record end time
    end_time = time.time()
    
    # Calculate duration
    duration = (end_time - start_time) * 1000  # Convert to ms
    
    # Add to kernel times
    if kernel_name not in kernel_times:
        kernel_times[kernel_name] = []
    kernel_times[kernel_name].append(duration)
    
    return event

# Replace the original function with our profiling function
cl.enqueue_nd_range_kernel = profiling_enqueue_nd_range_kernel
