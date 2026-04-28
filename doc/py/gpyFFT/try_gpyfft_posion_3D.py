import numpy as np
import pyopencl as cl
import pyopencl.array as cla
from gpyfft.fft import FFT
import matplotlib.pyplot as plt

# Create OpenCL context and command queue
ctx = cl.create_some_context(interactive=False)
queue = cl.CommandQueue(ctx)

# Define grid dimensions
nx, ny, nz = 16, 16, 16

# Step 1: Initialize a 3D array of zeros
data_host = np.zeros((nx, ny, nz), dtype=np.complex64)

# Step 2: Set one voxel to 1.0 (e.g., at position (8, 8, 8))
data_host[8, 8, 8] = 1.0 + 0.0j

# Transfer data to GPU
data_gpu = cla.to_device(queue, data_host)

# Step 3: Create FFT transform plan for the entire grid
transform = FFT(ctx, queue, data_gpu, axes=(0, 1, 2))

# Start computation and wait until it is finished
event, = transform.enqueue()
event.wait()

# Read back the transformed data from GPU to host
fft_result_host = data_gpu.get()

# Step 4: Create k-space coordinates for multiplication
kx = np.fft.fftfreq(nx)[:, None, None] * nx
ky = np.fft.fftfreq(ny)[None, :, None] * ny
kz = np.fft.fftfreq(nz)[None, None, :] * nz

# Compute kx^2 + ky^2 + kz^2
k_squared = kx**2 + ky**2 + kz**2

# Avoid division by zero by setting the zero frequency component to a small number
k_squared[0, 0, 0] = np.finfo(np.float32).eps

# Multiply in Fourier space by 1/(kx^2 + ky^2 + kz^2)
fft_result_host /= k_squared

# Step 5: Create an inverse FFT transform plan (reuse transform)
inverse_transform = FFT(ctx, queue, data_gpu, axes=(0, 1, 2))

# Copy modified data back to GPU for inverse FFT
data_gpu.set(fft_result_host)

# Start inverse computation and wait until it is finished
#event, = inverse_transform.enqueue(inverse=True)
event, = inverse_transform.enqueue()
event.wait()

# Read back the transformed data from GPU to host (solution in real space)
solution_host = data_gpu.get()

# Step 6: Visualize a specific slice using imshow

# Choose a slice index (e.g., iz=8)
iz = 8

# Extract the slice from the solution (real part for visualization)
slice_data = np.abs(solution_host[:, :, iz])

# Plot using imshow
plt.figure(figsize=(10, 5))
plt.subplot(1, 2, 1); plt.imshow(data_host[:, :, iz].T.real,  origin='lower', cmap='viridis', aspect='auto'); plt.axis('equal'); plt.colorbar(); plt.title(f'Charge')
plt.subplot(1, 2, 2); plt.imshow(slice_data.T,                origin='lower', cmap='viridis', aspect='auto'); plt.axis('equal'); plt.colorbar(); plt.title(f'Potential')

# 
# plt.xlabel('X axis')
# plt.ylabel('Y axis')
plt.show()