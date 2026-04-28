import numpy as np
import pyopencl as cl
from gpyfft.fft import FFT

# Initialize OpenCL context and command queue
ctx = cl.create_some_context()
queue = cl.CommandQueue(ctx)

# Define grid dimensions
nx, ny, nz = 16, 16, 16
nxyz = nx * ny * nz

# Create a buffer on the GPU for input data (complex numbers)
input_buffer = cl.Buffer(ctx, cl.mem_flags.READ_WRITE,
                          size=nxyz * np.dtype(np.complex64).itemsize)

# Initialize data with zeros and set one voxel to 1.0 (e.g., at index 8)
input_data = np.zeros((nxyz,), dtype=np.complex64)
input_data[8] = 1.0 + 0.0j  # Set voxel value

# Copy initial data to GPU buffer
cl.enqueue_copy(queue, input_buffer, input_data)

# Create a dummy NumPy array for shape (needed for creating the FFT plan)
shape = (nx, ny, nz)
dummy_array = np.empty(shape, dtype=np.complex64)

# Create FFT transform plan for 3D data (in-place operation)
transform = FFT(ctx, queue, dummy_array, axes=(0, 1, 2))

# Enqueue the FFT operation using the input buffer
event, = transform.enqueue(input_buffer)
event.wait()

# Now input_buffer contains the transformed data in Fourier space
# To perform inverse FFT in-place:
inverse_transform = FFT(ctx, queue, dummy_array, axes=(0, 1, 2))
event_inverse, = inverse_transform.enqueue(input_buffer, inverse=True)
event_inverse.wait()

# Read back the result from input_buffer if needed
result_data = np.empty_like(input_data)
cl.enqueue_copy(queue, result_data, input_buffer)

print("Result after inverse FFT:")
print(result_data)