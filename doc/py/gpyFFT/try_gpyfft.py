#!/usr/bin/python

import numpy as np
import pyopencl as cl
import pyopencl.array as cla
from gpyfft.fft import FFT

# Create OpenCL context and command queue
ctx = cl.create_some_context(interactive=False)
queue = cl.CommandQueue(ctx)

# Prepare 3D data (4 batches of 16x16 complex numbers)
data_shape = (4, 16, 16)
data_host = np.random.rand(*data_shape).astype(np.complex64) + 1j * np.random.rand(*data_shape).astype(np.complex64)
data_gpu = cla.to_device(queue, data_host)

# Create FFT transform plan for batched inline 3D transform
transform = FFT(ctx, queue, data_gpu, axes=(1, 2))

# Start computation and wait until it is finished
event, = transform.enqueue()
event.wait()

# Read back the transformed data from GPU to host
result_host = data_gpu.get()

# Output results (optional)
print("Transformed Data:")
print(result_host)