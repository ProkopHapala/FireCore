import numpy as np
import matplotlib.pyplot as plt
import os

# Path to the data directory
data_dir = "/home/indranil/git/FireCore/tests/tMMFF/data"
nacl_dir = os.path.join(data_dir, "NaCl_1x1_L1")

# Load CPU potential files (before B-spline fitting)
cpu_paul_before = np.load(os.path.join(nacl_dir, "VPaul.npy"))
cpu_lond_before = np.load(os.path.join(nacl_dir, "VLond.npy"))
cpu_coul_before = np.load(os.path.join(nacl_dir, "VCoul.npy"))

# Load GPU potential files (before B-spline fitting)
gpu_paul_before = np.load(os.path.join(nacl_dir, "V_Paul_gpu_before.npy"))
gpu_lond_before = np.load(os.path.join(nacl_dir, "V_Lond_gpu_before.npy"))
gpu_coul_before = np.load(os.path.join(nacl_dir, "V_Coul_gpu_before.npy"))

# Load CPU potential files (after B-spline fitting)
cpu_paul_after = np.load(os.path.join(nacl_dir, "debug_BsplinePaul_pbc.npy"))
cpu_lond_after = np.load(os.path.join(nacl_dir, "debug_BsplineLond_pbc.npy"))
cpu_coul_after = np.load(os.path.join(nacl_dir, "debug_BsplineCoul_pbc.npy"))

# Load GPU potential files (after B-spline fitting)
gpu_paul_after = np.load(os.path.join(nacl_dir, "V_Paul_gpu_after.npy"))
gpu_lond_after = np.load(os.path.join(nacl_dir, "V_Lond_gpu_after.npy"))
gpu_coul_after = np.load(os.path.join(nacl_dir, "V_Coul_gpu_after.npy"))

# Function to ensure all arrays have the same shape
def ensure_same_shape(arrays):
    # Find minimum shape across all arrays
    min_shape = [float('inf'), float('inf'), float('inf')]
    for arr in arrays:
        for i in range(3):
            min_shape[i] = min(min_shape[i], arr.shape[i])
    
    # Slice all arrays to the minimum shape
    result = []
    for arr in arrays:
        result.append(arr[:min_shape[0], :min_shape[1], :min_shape[2]])
    
    return result

# Ensure all arrays have the same shape
(
    cpu_paul_before, cpu_lond_before, cpu_coul_before,
    gpu_paul_before, gpu_lond_before, gpu_coul_before,
    cpu_paul_after, cpu_lond_after, cpu_coul_after,
    gpu_paul_after, gpu_lond_after, gpu_coul_after
) = ensure_same_shape([
    cpu_paul_before, cpu_lond_before, cpu_coul_before,
    gpu_paul_before, gpu_lond_before, gpu_coul_before,
    cpu_paul_after, cpu_lond_after, cpu_coul_after,
    gpu_paul_after, gpu_lond_after, gpu_coul_after
])

# Print shapes and ranges for before fitting
print("\nDifference between CPU and GPU before fitting ")
print("CPU Paul shape:", cpu_paul_before.shape, "range:", cpu_paul_before.min(), cpu_paul_before.max())
print("GPU Paul shape:", gpu_paul_before.shape, "range:", gpu_paul_before.min(), gpu_paul_before.max())
print("CPU Lond shape:", cpu_lond_before.shape, "range:", cpu_lond_before.min(), cpu_lond_before.max())
print("GPU Lond shape:", gpu_lond_before.shape, "range:", gpu_lond_before.min(), gpu_lond_before.max())
print("CPU Coul shape:", cpu_coul_before.shape, "range:", cpu_coul_before.min(), cpu_coul_before.max())
print("GPU Coul shape:", gpu_coul_before.shape, "range:", gpu_coul_before.min(), gpu_coul_before.max())

# Print shapes and ranges for after fitting
print("\nDifference between CPU and GPU after fitting ")
print("CPU Paul shape:", cpu_paul_after.shape, "range:", cpu_paul_after.min(), cpu_paul_after.max())
print("GPU Paul shape:", gpu_paul_after.shape, "range:", gpu_paul_after.min(), gpu_paul_after.max())
print("CPU Lond shape:", cpu_lond_after.shape, "range:", cpu_lond_after.min(), cpu_lond_after.max())
print("GPU Lond shape:", gpu_lond_after.shape, "range:", gpu_lond_after.min(), gpu_lond_after.max())
print("CPU Coul shape:", cpu_coul_after.shape, "range:", cpu_coul_after.min(), cpu_coul_after.max())
print("GPU Coul shape:", gpu_coul_after.shape, "range:", gpu_coul_after.min(), gpu_coul_after.max())

# Calculate differences
diff_paul_before = cpu_paul_before - gpu_paul_before
diff_lond_before = cpu_lond_before - gpu_lond_before
diff_coul_before = cpu_coul_before - gpu_coul_before

diff_paul_after = cpu_paul_after - gpu_paul_after
diff_lond_after = cpu_lond_after - gpu_lond_after
diff_coul_after = cpu_coul_after - gpu_coul_after

# Calculate CPU before-after differences
diff_cpu_paul = cpu_paul_before - cpu_paul_after
diff_cpu_lond = cpu_lond_before - cpu_lond_after
diff_cpu_coul = cpu_coul_before - cpu_coul_after

# Calculate GPU before-after differences
diff_gpu_paul = gpu_paul_before - gpu_paul_after
diff_gpu_lond = gpu_lond_before - gpu_lond_after
diff_gpu_coul = gpu_coul_before - gpu_coul_after

# Plot 2D slices of the potentials (middle slice)
z_slice = cpu_paul_before.shape[0] // 2

# Plot before fitting
plt.figure(figsize=(15, 10))

# Paul potential
plt.subplot(3, 3, 1)
plt.imshow(cpu_paul_before[z_slice], origin='lower')
plt.colorbar()
plt.title('CPU Paul (Before Fitting)')

plt.subplot(3, 3, 2)
plt.imshow(gpu_paul_before[z_slice], origin='lower')
plt.colorbar()
plt.title('GPU Paul (Before Fitting)')

plt.subplot(3, 3, 3)
plt.imshow(diff_paul_before[z_slice], origin='lower', cmap='coolwarm')
plt.colorbar()
plt.title('Diff Paul (Before Fitting)')

# London potential
plt.subplot(3, 3, 4)
plt.imshow(cpu_lond_before[z_slice], origin='lower')
plt.colorbar()
plt.title('CPU London (Before Fitting)')

plt.subplot(3, 3, 5)
plt.imshow(gpu_lond_before[z_slice], origin='lower')
plt.colorbar()
plt.title('GPU London (Before Fitting)')

plt.subplot(3, 3, 6)
plt.imshow(diff_lond_before[z_slice], origin='lower', cmap='coolwarm')
plt.colorbar()
plt.title('Diff London (Before Fitting)')

# Coulomb potential
plt.subplot(3, 3, 7)
plt.imshow(cpu_coul_before[z_slice], origin='lower')
plt.colorbar()
plt.title('CPU Coulomb (Before Fitting)')

plt.subplot(3, 3, 8)
plt.imshow(gpu_coul_before[z_slice], origin='lower')
plt.colorbar()
plt.title('GPU Coulomb (Before Fitting)')

plt.subplot(3, 3, 9)
plt.imshow(diff_coul_before[z_slice], origin='lower', cmap='coolwarm')
plt.colorbar()
plt.title('Diff Coulomb (Before Fitting)')

plt.tight_layout()
plt.savefig('potential_comparison_before.png')
plt.show()

# Plot after fitting
plt.figure(figsize=(15, 10))

# Paul potential
plt.subplot(3, 3, 1)
plt.imshow(cpu_paul_after[z_slice], origin='lower')
plt.colorbar()
plt.title('CPU Paul (After Fitting)')

plt.subplot(3, 3, 2)
plt.imshow(gpu_paul_after[z_slice], origin='lower')
plt.colorbar()
plt.title('GPU Paul (After Fitting)')

plt.subplot(3, 3, 3)
plt.imshow(diff_paul_after[z_slice], origin='lower', cmap='coolwarm')
plt.colorbar()
plt.title('Diff Paul (After Fitting)')

# London potential
plt.subplot(3, 3, 4)
plt.imshow(cpu_lond_after[z_slice], origin='lower')
plt.colorbar()
plt.title('CPU London (After Fitting)')

plt.subplot(3, 3, 5)
plt.imshow(gpu_lond_after[z_slice], origin='lower')
plt.colorbar()
plt.title('GPU London (After Fitting)')

plt.subplot(3, 3, 6)
plt.imshow(diff_lond_after[z_slice], origin='lower', cmap='coolwarm')
plt.colorbar()
plt.title('Diff London (After Fitting)')

# Coulomb potential
plt.subplot(3, 3, 7)
plt.imshow(cpu_coul_after[z_slice], origin='lower')
plt.colorbar()
plt.title('CPU Coulomb (After Fitting)')

plt.subplot(3, 3, 8)
plt.imshow(gpu_coul_after[z_slice], origin='lower')
plt.colorbar()
plt.title('GPU Coulomb (After Fitting)')

plt.subplot(3, 3, 9)
plt.imshow(diff_coul_after[z_slice], origin='lower', cmap='coolwarm')
plt.colorbar()
plt.title('Diff Coulomb (After Fitting)')

plt.tight_layout()
plt.savefig('potential_comparison_after.png')
plt.show()

# Plot histograms of the differences before fitting
plt.figure(figsize=(15, 5))

plt.subplot(1, 3, 1)
plt.hist(diff_paul_before.flatten(), bins=100)
plt.title('Paul Difference Histogram (Before Fitting)')

plt.subplot(1, 3, 2)
plt.hist(diff_lond_before.flatten(), bins=100)
plt.title('London Difference Histogram (Before Fitting)')

plt.subplot(1, 3, 3)
plt.hist(diff_coul_before.flatten(), bins=100)
plt.title('Coulomb Difference Histogram (Before Fitting)')

plt.tight_layout()
plt.savefig('difference_histograms_before.png')
plt.show()

# Plot histograms of the differences after fitting
plt.figure(figsize=(15, 5))

plt.subplot(1, 3, 1)
plt.hist(diff_paul_after.flatten(), bins=100)
plt.title('Paul Difference Histogram (After Fitting)')

plt.subplot(1, 3, 2)
plt.hist(diff_lond_after.flatten(), bins=100)
plt.title('London Difference Histogram (After Fitting)')

plt.subplot(1, 3, 3)
plt.hist(diff_coul_after.flatten(), bins=100)
plt.title('Coulomb Difference Histogram (After Fitting)')

plt.tight_layout()
plt.savefig('difference_histograms_after.png')
plt.show()

# Advanced Analysis: Find the exact locations of maximum differences
print("\n=== DETAILED ANALYSIS OF MAXIMUM DIFFERENCES ===")

for name, cpu, gpu, diff in [
    ('Paul (before fitting)', cpu_paul_before, gpu_paul_before, diff_paul_before),
    ('London (before fitting)', cpu_lond_before, gpu_lond_before, diff_lond_before),
    ('Coulomb (before fitting)', cpu_coul_before, gpu_coul_before, diff_coul_before),
    ('Paul (after fitting)', cpu_paul_after, gpu_paul_after, diff_paul_after),
    ('London (after fitting)', cpu_lond_after, gpu_lond_after, diff_lond_after),
    ('Coulomb (after fitting)', cpu_coul_after, gpu_coul_after, diff_coul_after)
]:
    # Find max absolute difference and its location
    abs_diff = np.abs(diff)
    max_idx = np.unravel_index(np.argmax(abs_diff), abs_diff.shape)
    z, y, x = max_idx
    
    # Calculate relative error at this point
    eps = 1e-10  # To avoid division by zero
    rel_error = abs_diff[max_idx] / (np.abs(cpu[max_idx]) + eps)
    
    print(f"\n{name}:")
    print(f"  Max abs diff at (z={z}, y={y}, x={x}): {abs_diff[max_idx]:.6e}")
    print(f"  CPU value: {cpu[max_idx]:.6e}, GPU value: {gpu[max_idx]:.6e}")
    print(f"  Relative error: {rel_error:.6e}")
    
    # Calculate overall statistics
    mean_abs_diff = np.mean(abs_diff)
    median_abs_diff = np.median(abs_diff)
    std_abs_diff = np.std(abs_diff)
    
    print(f"  Mean abs diff: {mean_abs_diff:.6e}")
    print(f"  Median abs diff: {median_abs_diff:.6e}")
    print(f"  Std dev of abs diff: {std_abs_diff:.6e}")

# Plot 2D slices at the points of maximum difference
plt.figure(figsize=(18, 12))
for i, (name, cpu, gpu, diff) in enumerate([
    ('Paul (before fitting)', cpu_paul_before, gpu_paul_before, diff_paul_before),
    ('London (before fitting)', cpu_lond_before, gpu_lond_before, diff_lond_before),
    ('Coulomb (before fitting)', cpu_coul_before, gpu_coul_before, diff_coul_before),
    ('Paul (after fitting)', cpu_paul_after, gpu_paul_after, diff_paul_after),
    ('London (after fitting)', cpu_lond_after, gpu_lond_after, diff_lond_after),
    ('Coulomb (after fitting)', cpu_coul_after, gpu_coul_after, diff_coul_after)
]):
    # Find max difference location
    abs_diff = np.abs(diff)
    max_idx = np.unravel_index(np.argmax(abs_diff), abs_diff.shape)
    z, y, x = max_idx
    
    # Plot the slice at z where the max difference occurs
    plt.subplot(2, 3, i+1)
    plt.imshow(diff[z], origin='lower', cmap='coolwarm')
    plt.colorbar(label='CPU - GPU')
    plt.title(f"{name}\nSlice z={z}")
    plt.scatter(x, y, color='black', marker='x', s=50)  # Mark the max diff point
    
plt.tight_layout()
plt.savefig('max_difference_slices.png')
plt.show()

# Check for numerical precision issues
print("\n=== NUMERICAL PRECISION ANALYSIS ===")
for name, cpu, gpu in [
    ('Paul (before fitting)', cpu_paul_before, gpu_paul_before),
    ('London (before fitting)', cpu_lond_before, gpu_lond_before),
    ('Coulomb (before fitting)', cpu_coul_before, gpu_coul_before),
    ('Paul (after fitting)', cpu_paul_after, gpu_paul_after),
    ('London (after fitting)', cpu_lond_after, gpu_lond_after),
    ('Coulomb (after fitting)', cpu_coul_after, gpu_coul_after)
]:
    # Check for very small values that might cause precision issues
    small_vals_cpu = np.sum(np.abs(cpu) < 1e-10)
    small_vals_gpu = np.sum(np.abs(gpu) < 1e-10)
    
    # Check for very large values that might cause overflow
    large_vals_cpu = np.sum(np.abs(cpu) > 1e6)
    large_vals_gpu = np.sum(np.abs(gpu) > 1e6)
    
    print(f"\n{name} precision analysis:")
    print(f"  Small values (< 1e-10): CPU: {small_vals_cpu}, GPU: {small_vals_gpu}")
    print(f"  Large values (> 1e6): CPU: {large_vals_cpu}, GPU: {large_vals_gpu}")
    
    # Check if the differences follow a pattern related to magnitude
    diff = cpu - gpu
    abs_cpu = np.abs(cpu)
    
    # Group by magnitude bins and check average error
    bins = [0, 1e-6, 1e-3, 1e0, 1e3, np.inf]
    for i in range(len(bins)-1):
        mask = (abs_cpu >= bins[i]) & (abs_cpu < bins[i+1])
        if np.sum(mask) > 0:
            mean_diff = np.mean(np.abs(diff[mask]))
            print(f"  Mean abs diff for values {bins[i]:.1e}-{bins[i+1]:.1e}: {mean_diff:.6e} (count: {np.sum(mask)})")

# Check for spatial patterns in the differences
plt.figure(figsize=(15, 5))

# Create projections of the absolute differences along each axis
for i, (name, diff) in enumerate([
    ('Paul', diff_paul_after),
    ('London', diff_lond_after),
    ('Coulomb', diff_coul_after)
]):
    abs_diff = np.abs(diff)
    
    plt.subplot(1, 3, i+1)
    
    # Project along each axis
    z_proj = np.mean(abs_diff, axis=(1, 2))
    y_proj = np.mean(abs_diff, axis=(0, 2))
    x_proj = np.mean(abs_diff, axis=(0, 1))
    
    plt.plot(range(len(z_proj)), z_proj, label='Z projection')
    plt.plot(range(len(y_proj)), y_proj, label='Y projection')
    plt.plot(range(len(x_proj)), x_proj, label='X projection')
    
    plt.title(f"{name} Difference Projections")
    plt.xlabel("Index")
    plt.ylabel("Mean Absolute Difference")
    plt.legend()
    plt.grid(True)

plt.tight_layout()
plt.savefig('difference_projections.png')
plt.show()

# Compare the effect of fitting on CPU vs GPU
print("\n=== FITTING EFFECT ANALYSIS ===")
for name, before_cpu, after_cpu, before_gpu, after_gpu in [
    ('Paul', cpu_paul_before, cpu_paul_after, gpu_paul_before, gpu_paul_after),
    ('London', cpu_lond_before, cpu_lond_after, gpu_lond_before, gpu_lond_after),
    ('Coulomb', cpu_coul_before, cpu_coul_after, gpu_coul_before, gpu_coul_after)
]:
    # Calculate the change due to fitting
    cpu_change = after_cpu - before_cpu
    gpu_change = after_gpu - before_gpu
    
    # Calculate the difference in the change
    change_diff = cpu_change - gpu_change
    
    print(f"\n{name} fitting effect:")
    print(f"  CPU change range: {cpu_change.min():.6e} to {cpu_change.max():.6e}")
    print(f"  GPU change range: {gpu_change.min():.6e} to {gpu_change.max():.6e}")
    print(f"  Difference in change range: {change_diff.min():.6e} to {change_diff.max():.6e}")
    
    # Find where the fitting effect differs most
    abs_change_diff = np.abs(change_diff)
    max_idx = np.unravel_index(np.argmax(abs_change_diff), abs_change_diff.shape)
    z, y, x = max_idx
    
    print(f"  Max difference in fitting effect at (z={z}, y={y}, x={x}): {abs_change_diff[max_idx]:.6e}")
    print(f"  CPU fitting change: {cpu_change[max_idx]:.6e}, GPU fitting change: {gpu_change[max_idx]:.6e}")

# Plot histograms of the fitting effect differences
plt.figure(figsize=(15, 5))

for i, (name, before_cpu, after_cpu, before_gpu, after_gpu) in enumerate([
    ('Paul', cpu_paul_before, cpu_paul_after, gpu_paul_before, gpu_paul_after),
    ('London', cpu_lond_before, cpu_lond_after, gpu_lond_before, gpu_lond_after),
    ('Coulomb', cpu_coul_before, cpu_coul_after, gpu_coul_before, gpu_coul_after)
]):
    # Calculate the change due to fitting
    cpu_change = after_cpu - before_cpu
    gpu_change = after_gpu - before_gpu
    
    # Calculate the difference in the change
    change_diff = cpu_change - gpu_change
    
    plt.subplot(1, 3, i+1)
    plt.hist(change_diff.flatten(), bins=100)
    plt.title(f'{name} Fitting Effect Difference')

plt.tight_layout()
plt.savefig('fitting_effect_difference.png')
plt.show()


""""
difference between cpu before and after fitting

CPU Paul shape: (200, 40, 40) range: 1.2548612781138576e-24 36.97090288159046
CPU Paul shape: (200, 40, 40) range: -2.383758917275688e-16 71.93161877100597
CPU Lond shape: (200, 40, 40) range: -3.9108023748996485 -2.6377442055651597e-12
CPU Lond shape: (200, 40, 40) range: -6.301866207073345 -2.875900022601296e-12
CPU Coul shape: (200, 40, 40) range: -52.01220096563635 52.01220096563635
CPU Coul shape: (200, 40, 40) range: -71.19295590474475 71.19295590474474


Difference between CPU and GPU before fitting 

CPU Paul shape: (200, 40, 40) range: 1.2548612781138576e-24 36.97090288159046
GPU Paul shape: (200, 40, 40) range: 1.5480292e-24 46.34928
CPU Lond shape: (200, 40, 40) range: -3.9108023748996485 -2.6377442055651597e-12
GPU Lond shape: (200, 40, 40) range: -4.582958 -3.0270166e-12
CPU Coul shape: (200, 40, 40) range: -52.01220096563635 52.01220096563635
GPU Coul shape: (200, 40, 40) range: -55.474205 55.474205


Difference between CPU and GPU after fitting 
CPU Paul shape: (200, 40, 40) range: -2.383758917275688e-16 71.93161877100597
GPU Paul shape: (200, 40, 40) range: -2.3511054e-22 90.17905
CPU Lond shape: (200, 40, 40) range: -6.301866207073345 -2.875900022601296e-12
GPU Lond shape: (200, 40, 40) range: -7.3873377 -3.300213e-12
CPU Coul shape: (200, 40, 40) range: -71.19295590474475 71.19295590474474
GPU Coul shape: (200, 40, 40) range: -81.06157 81.061554
"""