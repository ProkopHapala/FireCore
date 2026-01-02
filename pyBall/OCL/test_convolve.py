import argparse
import numpy as np
import pyopencl as cl
import time
import matplotlib.pyplot as plt
from scipy.ndimage import convolve
from convolve import Conv3D


def select_preferred_device(preferred_vendor="NVIDIA"):
    """
    Prefer a device whose vendor string contains the preferred keyword. Return
    the platform and device (first match) or fall back to the first available.
    """
    platforms = cl.get_platforms()
    for platform in platforms:
        for device in platform.get_devices():
            vendor = device.vendor.upper()
            if preferred_vendor.upper() in vendor or preferred_vendor.upper() in device.name.upper():
                return platform, device
    # fallback to first device if preferred is not found
    for platform in platforms:
        devices = platform.get_devices()
        if devices:
            return platform, devices[0]
    raise RuntimeError("No OpenCL devices found.")

def get_test_data(shape):
    """Generate a smooth 3D function (Gaussian blob + sine wave)."""
    z, y, x = np.meshgrid(
        np.linspace(-3, 3, shape[0]),
        np.linspace(-3, 3, shape[1]),
        np.linspace(-3, 3, shape[2]),
        indexing='ij'
    )
    # Gaussian blob in center
    blob = np.exp(-(x**2 + y**2 + z**2))
    # Sine wave to test periodicity
    wave = 0.2 * np.sin(x*2.0) * np.cos(y*2.0)
    return (blob + wave).astype(np.float32)

def run_cpu_reference(data, pbc):
    """
    Run equivalent convolution on CPU using Scipy.
    Kernel: B-Spline separable [1/6, 2/3, 1/6]
    """
    k1d = np.array([1.0/6.0, 2.0/3.0, 1.0/6.0])
    
    # Scipy convolve doesn't support a sequence of modes.
    # We apply the 1D convolution sequentially for each axis.
    # NOTE: The axes in OpenCL kernel are (ns.x, ns.y, ns.z) mapped to (W, H, D).
    # In numpy, indexing is 'ij' (z, y, x) -> (D, H, W).
    # so:
    # OpenCL ns.x (W) -> numpy axis 2
    # OpenCL ns.y (H) -> numpy axis 1
    # OpenCL ns.z (D) -> numpy axis 0
    # Our pbc input is (Z, Y, X) as per CLI help, but we need to check 
    # if the GPU results actually use (X, Y, Z) order for bPBC.
    res = data.copy()
    # Scipy axis 0 = Z, 1 = Y, 2 = X
    # OpenCL bPBC.x -> X, bPBC.y -> Y, bPBC.z -> Z
    # Thus pbc[0]=Z, pbc[1]=Y, pbc[2]=X matches axes 0, 1, 2.
    for axis, p in enumerate(pbc):
        mode = 'wrap' if p else 'constant'
        # Reshape k1d to be a 3D kernel along the current axis
        k_shape = [1, 1, 1]
        k_shape[axis] = 3
        k_axis = k1d.reshape(k_shape)
        res = convolve(res, k_axis, mode=mode, cval=0.0)
    return res

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run 3D B-spline convolution on a preferred GPU.")
    parser.add_argument( "--grid-size", nargs=3, type=int,  default=[40,50,60], metavar=("Z", "Y", "X"), help="Grid dimensions (Z Y X).")
    parser.add_argument( "--ls-values", nargs="+", type=int,  default=[8], metavar="LS", help="Uniform local sizes to test along each axis.")
    parser.add_argument( "--pbc-axes", nargs=3, type=int, choices=(0,1), default=(1,1,1), help="Periodic boundary flags (Z Y X) as 1=periodic,0=zero.")
    parser.add_argument( "--random", type=int, default=1, help="Populate the grid with random values instead of the deterministic test pattern.")
    parser.add_argument( "--seed", type=int, default=42, help="Random seed when --random is enabled.")
    parser.add_argument( "--plot", type=int, default=1, help="Generate comparison plots for CPU vs GPU output.")
    args = parser.parse_args()

    # 1. Setup OpenCL
    platform, device = select_preferred_device()
    ctx      = cl.Context([device])
    queue    = cl.CommandQueue(ctx)
    
    print(f"Using Device: {device.name}")
    print("-" * 50)

    # 2. Test Parameters
    grid_shape = tuple(args.grid_size)
    pbc = tuple(bool(val) for val in args.pbc_axes)
    test_configs = [{"ls": (ls, ls, ls), "pbc": pbc} for ls in args.ls_values]

    if args.random:
        rng = np.random.default_rng(args.seed)
        data = rng.standard_normal(grid_shape).astype(np.float32)
    else:
        data = get_test_data(grid_shape)
    
    # 3. Iterate Configurations
    for cfg in test_configs:
        ls = cfg['ls']
        pbc = cfg['pbc']
        
        print(f"\nTesting Config: Tile={ls}, PBC={pbc}")
        
        # Initialize Module (Compiles Kernel with specific LS)
        bspline = Conv3D(ctx, queue, tile_size=ls)
        
        # Run GPU
        t0 = time.time()
        # OpenCL kernel expects bPBC as {nx, ny, nz, 0}
        # But wait, our pbc is (Z, Y, X).
        # Convolution3D_General( ns={nx,ny,nz,0}, ..., bPBC={px, py, pz, 0} )
        # So we MUST swap pbc order to (X, Y, Z) for the kernel call
        pbc_gpu = (pbc[2], pbc[1], pbc[0])
        res_gpu = bspline.run(data, pbc=pbc_gpu)
        queue.finish()
        gpu_time = (time.time() - t0) * 1000
        
        # Run CPU Reference
        res_cpu = run_cpu_reference(data, pbc)
        
        diff = np.abs(res_gpu - res_cpu)
        max_diff = np.max(diff)
        mean_diff = np.mean(diff)
        
        print(f"  GPU Time: {gpu_time:.4f} ms")
        print(f"  Max Error: {max_diff:.2e}")
        print(f"  Mean Error: {mean_diff:.2e}")
        
        if max_diff < 1e-5:
            print("  [PASS] Results match reference.")
        else:
            print("  [FAIL] Artifacts detected! Check corners:")
            print(f"    CPU corner: {res_cpu[0,0,0]:.6f}")
            print(f"    GPU corner: {res_gpu[0,0,0]:.6f}")

        if args.plot:
            midpoint = tuple(dim//2 for dim in grid_shape)
            slices = {
                "CPU": res_cpu[midpoint[0], :, :],
                "GPU": res_gpu[midpoint[0], :, :],
                "Diff": diff[midpoint[0], :, :]
            }
            fig, axes = plt.subplots(1, 3, figsize=(12, 4))
            for ax, (title, array) in zip(axes, slices.items()):
                im = ax.imshow(array, origin='lower', cmap='coolwarm')
                ax.set_title(f"{title} slice @ z={midpoint[0]}")
                plt.colorbar(im, ax=ax, shrink=0.8)
            plt.tight_layout()
            plt.show()