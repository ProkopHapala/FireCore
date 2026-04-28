#!/usr/bin/python

import numpy as np
import pyopencl as cl
from pyfft.cl import Plan

# Create OpenCL context and command queue
ctx = cl.create_some_context(interactive=False)
queue = cl.CommandQueue(ctx)

# Prepare data
data = np.ones((16, 16), dtype=np.complex64)
gpu_data = cl.array.to_device(queue, data)

# Create FFT plan
plan = Plan(data.shape, queue=queue)

# Execute forward FFT
plan.execute(gpu_data.data)
result = gpu_data.get()

# Execute inverse FFT
plan.execute(gpu_data.data, inverse=True)
inverse_result = gpu_data.get()

# Check results
error = np.abs(np.sum(np.abs(data) - np.abs(inverse_result)) / data.size)
print("Error:", error < 1e-6)  # Should print True if correct
