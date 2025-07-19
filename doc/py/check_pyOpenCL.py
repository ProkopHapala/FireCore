import pyopencl as cl
import numpy as np

# OpenCL kernel to add two arrays
kernel_code = """
__kernel void add_arrays(__global const float *a, 
                         __global const float *b, 
                         __global float *result) {
    int gid = get_global_id(0);
    result[gid] = a[gid] + b[gid];
}
"""

def main():
    # Get the list of OpenCL platforms
    platforms = cl.get_platforms()
    
    if not platforms:
        print("No OpenCL platforms found.")
        return

    for platform in platforms:
        print(f"Platform: {platform.name}")
        print(f"Version: {platform.version}")

        # Get the list of devices for each platform
        devices = platform.get_devices()
        for device in devices:
            print(f"  Device: {device.name}")
            print(f"  Type: {cl.device_type.to_string(device.type)}")
            print(f"  Max Compute Units: {device.max_compute_units}")
            print(f"  Global Memory Size: {device.global_mem_size / (1024 ** 2):.2f} MB")

            # Create a context and command queue
            context = cl.Context([device])
            queue = cl.CommandQueue(context)

            # Create input data
            n = 1024  # Size of the arrays
            a = np.random.rand(n).astype(np.float32)
            b = np.random.rand(n).astype(np.float32)
            result = np.empty_like(a)

            # Create OpenCL buffers
            mf = cl.mem_flags
            a_buf = cl.Buffer(context, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=a)
            b_buf = cl.Buffer(context, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=b)
            result_buf = cl.Buffer(context, mf.WRITE_ONLY, result.nbytes)

            # Build the kernel
            program = cl.Program(context, kernel_code).build()

            # Execute the kernel
            program.add_arrays(queue, a.shape, None, a_buf, b_buf, result_buf)

            # Read back the results
            cl.enqueue_copy(queue, result, result_buf)

            # Verify the results
            expected_result = a + b

            print("a:", a)
            print("b:", b)
            print("result          :", result )
            print("expected_result :", expected_result )
            if np.allclose(result, expected_result):
                print("Array addition successful!")
            else:
                print("Array addition failed!")

if __name__ == "__main__":
    main()