import pyopencl as cl
import numpy as np
import sys

print("--- PyOpenCL GPU Detection Test ---")

# 1. Find NVIDIA Platform and GPU Device
nvidia_platform = None
nvidia_device = None

try:
    print("Available OpenCL Platforms:")
    platforms = cl.get_platforms()
    if not platforms:
        print("!!! ERROR: No OpenCL platforms found!")
        sys.exit(1)

    for i, platform in enumerate(platforms):
        platform_name = platform.get_info(cl.platform_info.NAME)
        print(f"  [{i}] Platform Name: {platform_name}")
        if 'nvidia' in platform_name.lower():
            nvidia_platform = platform
            print(f"      (Selected as NVIDIA Platform)")

    if nvidia_platform is None:
        print("\n!!! ERROR: NVIDIA OpenCL platform not found!")
        sys.exit(1)

    print("\nAvailable GPU Devices on NVIDIA Platform:")
    devices = nvidia_platform.get_devices(device_type=cl.device_type.GPU)
    if not devices:
        print("!!! ERROR: No GPU devices found on NVIDIA platform!")
        sys.exit(1)

    for i, device in enumerate(devices):
        device_name = device.get_info(cl.device_info.NAME)
        print(f"  [{i}] Device Name: {device_name}")
        # Select the first available NVIDIA GPU
        if nvidia_device is None:
            nvidia_device = device
            print(f"      (Selected this device for context)")

    if nvidia_device is None:
        print("\n!!! ERROR: Could not select an NVIDIA GPU device!")
        sys.exit(1)

    print(f"\nUsing Device: {nvidia_device.get_info(cl.device_info.NAME)}")

except cl.Error as e:
    print(f"\n!!! OpenCL Error during platform/device selection: {e}")
    sys.exit(1)

# 2. Create Context and Command Queue
try:
    ctx = cl.Context(devices=[nvidia_device])
    # Add PROFILING_ENABLE if you want to try clGetEventProfilingInfo later
    queue = cl.CommandQueue(ctx, properties=cl.command_queue_properties.PROFILING_ENABLE)
    print("OpenCL Context and Command Queue created successfully.")
except cl.Error as e:
    print(f"!!! OpenCL Error creating context/queue: {e}")
    sys.exit(1)

# 3. Prepare Data
N = 1024 * 1024  # Size of vectors (e.g., 1 Million elements)
try:
    a_np = np.random.rand(N).astype(np.float32)
    b_np = np.random.rand(N).astype(np.float32)
    res_np = np.empty_like(a_np)
    print(f"Host data created (Size: {N} floats)")
except Exception as e:
     print(f"!!! Error creating numpy arrays: {e}")
     sys.exit(1)

# 4. Create Device Buffers and Copy Data
mf = cl.mem_flags
try:
    a_g = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=a_np)
    b_g = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=b_np)
    res_g = cl.Buffer(ctx, mf.WRITE_ONLY, res_np.nbytes)
    print("Device buffers created and data copied.")
    # Explicitly finish queue after initial copies (good practice)
    queue.finish()
except cl.Error as e:
    print(f"!!! OpenCL Error creating buffers or initial copy: {e}")
    sys.exit(1)

# 5. Define and Build Kernel
kernel_code = """
__kernel void simple_add(__global const float *a_g,
                         __global const float *b_g,
                         __global float *res_g)
{
  int gid = get_global_id(0);
  res_g[gid] = a_g[gid] + b_g[gid];
}
"""
try:
    prg = cl.Program(ctx, kernel_code).build()
    print("Kernel built successfully.")
except cl.Error as e:
    print(f"!!! OpenCL Error building kernel: {e}")
    sys.exit(1)

# 6. Execute Kernel
print("Executing kernel...")
try:
    # Get the kernel function
    simple_add_knl = prg.simple_add
    # Set global size; local size can often be None for OpenCL to choose
    global_size = (N,)
    local_size = None # Or set to a specific value like (256,) if needed
    # Launch the kernel
    kernel_event = simple_add_knl(queue, global_size, local_size, a_g, b_g, res_g)
    # Wait for the kernel event to complete (optional, but useful)
    kernel_event.wait()
    print("Kernel execution enqueued and waited for.")
except cl.Error as e:
    print(f"!!! OpenCL Error executing kernel: {e}")
    sys.exit(1)

# 7. Finish Queue (Ensure kernel execution completes and catch runtime errors)
print("Finishing command queue...")
try:
    queue.finish()
    print("Queue finished.")
except cl.Error as e:
    print(f"!!! OpenCL Error during queue.finish() (Kernel execution might have failed): {e}")
    sys.exit(1)

# 8. Read Result Back
print("Reading results back...")
try:
    read_event = cl.enqueue_copy(queue, res_np, res_g)
    # Wait for the read to complete before verification
    read_event.wait()
    print("Results read back.")
except cl.Error as e:
    print(f"!!! OpenCL Error reading results: {e}")
    sys.exit(1)

# 9. Verify Result (Optional)
print("Verifying results...")
try:
    expected_res = a_np + b_np
    if np.allclose(res_np, expected_res):
        print("--- Verification SUCCESSFUL! ---")
    else:
        diff = np.abs(res_np - expected_res)
        print(f"!!! Verification FAILED! Max difference: {np.max(diff)}")
except Exception as e:
    print(f"!!! Error during verification: {e}")

print("--- Test script finished. ---")