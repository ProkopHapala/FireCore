#!/usr/bin/python

import numpy as np
import pyopencl as cl
import pyopencl.array as cla
from gpyfft.fft import FFT

# Create OpenCL context and command queue
ctx = cl.create_some_context(interactive=False)
queue = cl.CommandQueue(ctx)

# Define grid dimensions
nx, ny, nz = 16, 16, 16

# Create a uniform 3D array with random complex values
data_host = np.random.rand(nx, ny, nz).astype(np.complex64) + \
 1j * np.random.rand(nx, ny, nz).astype(np.complex64)

# Transfer data to GPU
data_gpu = cla.to_device(queue, data_host)

# Create FFT transform plan for the entire grid
transform = FFT(ctx, queue, data_gpu, axes=(0, 1, 2))

# Start computation and wait until it is finished
event, = transform.enqueue()
event.wait()

# Read back the transformed data from GPU to host
result_host = data_gpu.get()

# Output results (optional)
print("Transformed Data:")
print(result_host)