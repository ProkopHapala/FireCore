#!/usr/bin/env python3

import inspect
import sys

try:
    import pyopencl as cl
    
    # Print PyOpenCL version
    print(f"PyOpenCL version: {cl.VERSION}")
    
    # Examine CommandQueue class
    print("\nCommandQueue methods:")
    for name, method in inspect.getmembers(cl.CommandQueue, predicate=inspect.isfunction):
        print(f"  {name}: {method}")
    
    # Examine low-level _cl module
    print("\n_cl module attributes:")
    for name in dir(cl._cl):
        if not name.startswith('__'):
            print(f"  {name}")
    
    # Look for kernel execution functions
    print("\nKernel execution functions:")
    for name in dir(cl):
        if 'kernel' in name.lower() or 'enqueue' in name.lower():
            print(f"  {name}")
    
    # Create a context and examine a real queue
    print("\nExamining a real CommandQueue instance:")
    platform = cl.get_platforms()[0]
    device = platform.get_devices()[0]
    ctx = cl.Context([device])
    queue = cl.CommandQueue(ctx, properties=cl.command_queue_properties.PROFILING_ENABLE)
    
    print(f"  Queue properties: {queue.properties if hasattr(queue, 'properties') else 'Not available'}")
    print(f"  Queue methods:")
    for name, method in inspect.getmembers(queue, predicate=inspect.ismethod):
        print(f"    {name}: {method}")
    
    # Check if enqueue_kernel exists
    if hasattr(queue, 'enqueue_nd_range_kernel'):
        print("\nFound enqueue_nd_range_kernel method on queue instance")
    else:
        print("\nNo enqueue_nd_range_kernel method found on queue instance")
    
    # Try to find the actual kernel execution function
    print("\nSearching for kernel execution function in cl module:")
    for name in dir(cl):
        if callable(getattr(cl, name)) and ('kernel' in name.lower() or 'enqueue' in name.lower()):
            func = getattr(cl, name)
            print(f"  {name}: {func}")
            print(f"    Signature: {inspect.signature(func) if hasattr(inspect, 'signature') else 'Not available'}")
    
    print("\nDone.")
    
except ImportError:
    print("PyOpenCL not available")
except Exception as e:
    print(f"Error: {e}")
    import traceback
    traceback.print_exc()
