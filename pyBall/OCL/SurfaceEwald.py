"""
SurfaceEwald.py - OpenCL-accelerated 2D Ewald electrostatics

GPU implementation of 2D Ewald sum for surface electrostatics.
Reference: pyBall/Ewald2D.py (Python implementation)

Key features:
- Efficient complex multiplication for computing e^{iG·ρ}
- Precomputes z1_b1 = e^{i*b1·ρ}, z1_b2 = e^{i*b2·ρ}
- Uses powers: e^{ih*b1·ρ} = z1_b1^h
- Each thread evaluates one spatial point

Kernels:
- compute_ewald_coefficients: Compute C_G and w[g,i] coefficients
- eval_potential_vacuum: Evaluate potential for z > max(z_i)
- eval_potential_full: Evaluate potential for any z
- eval_potential_brute: Brute force Coulomb sum for validation
"""

import sys
import os
import numpy as np
import pyopencl as cl
import pyopencl.array as cl_array
import time

from . import clUtils as clu

COULOMB_CONST = 14.3996448915  # eV·Å/e²


class SurfaceEwaldCL:
    """
    OpenCL-accelerated 2D Ewald electrostatics evaluator.
    
    Usage:
        ew_cl = SurfaceEwaldCL()
        
        # Prepare data
        ion_data = np.column_stack([rx, ry, rz, q]).astype(np.float32)
        a_vec = np.array([ax, ay], dtype=np.float32)
        b_vec = np.array([bx, by], dtype=np.float32)
        
        # Compute coefficients
        ew_cl.prepare_system(ion_data, a_vec, b_vec, n_harm=3)
        
        # Evaluate potential
        X, Y = np.meshgrid(xv, yv)
        phi = ew_cl.eval_vacuum(X, Y, z_height)
    """
    
    def __init__(self, platform='nvidia'):
        """
        Initialize OpenCL context and compile kernels.
        
        Parameters:
            platform: 'nvidia', 'amd', or 'cpu' - which OpenCL platform to use
        """
        self.cl = cl
        
        # Get OpenCL context
        if platform == 'nvidia':
            self.ctx, self.queue = clu.get_nvidia_device(what="nvidia")
        elif platform == 'amd':
            self.ctx, self.queue = clu.get_amd_device()
        else:
            self.ctx, self.queue = clu.get_cpu_device()
        
        # Load and compile OpenCL program
        try:
            cl_path = os.path.normpath(os.path.join(
                os.path.dirname(os.path.abspath(__file__)),
                '../../cpp/common_resources/cl/Surface.cl'
            ))
            with open(cl_path, 'r') as f:
                cl_code = f.read()
            
            print(f"Compiling OpenCL program from {cl_path}...")
            self.prg = cl.Program(self.ctx, cl_code).build()
            print("OpenCL compilation successful")
            
        except Exception as e:
            print(f"Error compiling OpenCL program: {e}")
            import traceback
            traceback.print_exc()
            raise
        
        # System state
        self.ion_data = None          # float4 (rx, ry, rz, q)
        self.a_vec = None             # float2 (ax, ay)
        self.b_vec = None             # float2 (bx, by)
        self.b1 = None                # float2 reciprocal vector
        self.b2 = None                # float2 reciprocal vector
        self.area = None              # float
        self.n_harm = None            # int
        self.G_data = None            # float4 (h, k, Gn, 0)
        self.C_G = None               # float2 (real, imag)
        self.w = None                 # float2 per-ion weights
        
        # Buffers (allocated on demand)
        self.ion_buff = None
        self.G_buff = None
        self.b_vectors_buff = None
        self.C_G_buff = None
        self.w_buff = None
        
    def make_reciprocal_2d(self, a_vec, b_vec):
        """
        Compute reciprocal lattice vectors.
        
        Formula:
            area = |a × b|
            b1 = (2π/area) * (b_y, -b_x)
            b2 = (2π/area) * (-a_y, a_x)
        """
        a = np.asarray(a_vec, dtype=np.float64)
        b = np.asarray(b_vec, dtype=np.float64)
        area = abs(a[0]*b[1] - a[1]*b[0])
        
        b1 = (2*np.pi/area) * np.array([b[1], -b[0]])
        b2 = (2*np.pi/area) * np.array([-a[1], a[0]])
        
        return area, b1.astype(np.float32), b2.astype(np.float32)
    
    def generate_G_vectors(self, b1, b2, n_harm):
        """
        Generate G-vectors for |h|, |k| <= n_harm, excluding G=0.
        
        Returns:
            G_data: N_G x 4 array (h, k, Gn, 0)
        """
        G_list = []
        for h in range(-n_harm, n_harm+1):
            for k in range(-n_harm, n_harm+1):
                if h == 0 and k == 0:
                    continue
                Gx = h*b1[0] + k*b2[0]
                Gy = h*b1[1] + k*b2[1]
                Gn = np.sqrt(Gx**2 + Gy**2)
                G_list.append([h, k, Gn, 0])
        
        return np.array(G_list, dtype=np.float32)
    
    def prepare_system(self, ion_data, a_vec, b_vec, n_harm=3):
        """
        Prepare the Ewald system - compute all coefficients.
        
        Parameters:
            ion_data: (N_ions, 4) array - columns [rx, ry, rz, q]
            a_vec: (2,) array - lattice vector a
            b_vec: (2,) array - lattice vector b  
            n_harm: int - harmonic truncation (|h|,|k| <= n_harm)
        """
        print(f"Preparing Ewald2D system: N_ions={len(ion_data)}, n_harm={n_harm}")
        
        # Store system parameters
        self.ion_data = np.asarray(ion_data, dtype=np.float32)
        self.a_vec = np.asarray(a_vec, dtype=np.float32)
        self.b_vec = np.asarray(b_vec, dtype=np.float32)
        self.n_harm = n_harm
        
        # Compute reciprocal lattice
        self.area, self.b1, self.b2 = self.make_reciprocal_2d(a_vec, b_vec)
        print(f"  Area: {self.area:.4f} Å²")
        print(f"  Reciprocal b1: ({self.b1[0]:.4f}, {self.b1[1]:.4f})")
        print(f"  Reciprocal b2: ({self.b2[0]:.4f}, {self.b2[1]:.4f})")
        
        # Generate G-vectors
        self.G_data = self.generate_G_vectors(self.b1, self.b2, n_harm)
        N_G = len(self.G_data)
        print(f"  G-vectors: {N_G}")
        
        # Allocate buffers
        mf = cl.mem_flags
        self.ion_buff = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR,
                                  hostbuf=self.ion_data)
        self.G_buff = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR,
                                hostbuf=self.G_data)
        
        b_vectors = np.array([self.b1, self.b2], dtype=np.float32)
        self.b_vectors_buff = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR,
                                        hostbuf=b_vectors)
        
        # Output buffers for coefficients
        N_ions = len(self.ion_data)
        self.C_G = np.empty((N_G, 2), dtype=np.float32)
        self.w = np.empty((N_G, N_ions, 2), dtype=np.float32)
        
        self.C_G_buff = cl.Buffer(self.ctx, mf.WRITE_ONLY, self.C_G.nbytes)
        self.w_buff = cl.Buffer(self.ctx, mf.WRITE_ONLY, self.w.nbytes)
        
        # Launch coefficient computation kernel
        t0 = time.time()
        self.prg.compute_ewald_coefficients(
            self.queue, (N_G,), None,
            self.ion_buff, self.G_buff, self.b_vectors_buff,
            np.float32(self.area), np.int32(N_ions), np.int32(N_G),
            self.C_G_buff, self.w_buff
        )
        
        # Copy results back
        cl.enqueue_copy(self.queue, self.C_G, self.C_G_buff)
        cl.enqueue_copy(self.queue, self.w, self.w_buff)
        self.queue.finish()
        
        t1 = time.time()
        print(f"  Coefficient computation: {t1-t0:.3f} s")
        
        # Make C_G and w available as device buffers for evaluation
        self.C_G_buff = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR,
                                  hostbuf=self.C_G)
        self.w_buff = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR,
                                hostbuf=self.w.reshape(-1, 2))
    
    def eval_vacuum(self, X, Y, z, use_fast=False):
        """
        Evaluate potential on XY grid at fixed height z.
        Valid for z > max(z_i) (above all ions).
        
        Parameters:
            X, Y: 2D arrays - meshgrid coordinates
            z: float - height above surface
            use_fast: bool - use local memory optimization
            
        Returns:
            phi: 2D array - potential in eV/e
        """
        if self.C_G is None:
            raise ValueError("System not prepared. Call prepare_system() first.")
        
        # Flatten evaluation points
        x_flat = X.ravel().astype(np.float32)
        y_flat = Y.ravel().astype(np.float32)
        N_points = len(x_flat)
        
        eval_points = np.column_stack([x_flat, y_flat,
                                      np.full(N_points, z, dtype=np.float32),
                                      np.zeros(N_points, dtype=np.float32)])
        
        # Allocate buffers
        mf = cl.mem_flags
        eval_buff = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR,
                              hostbuf=eval_points)
        phi_out = np.empty(N_points, dtype=np.float32)
        phi_buff = cl.Buffer(self.ctx, mf.WRITE_ONLY, phi_out.nbytes)
        
        # Launch kernel
        N_G = len(self.G_data)
        t0 = time.time()
        
        if use_fast and self.n_harm <= 16:
            # Use local memory optimization for small n_harm
            local_size = min(64, N_points)
            global_size = ((N_points + local_size - 1) // local_size) * local_size
            self.prg.eval_potential_vacuum_fast(
                self.queue, (global_size,), (local_size,),
                eval_buff, self.C_G_buff, self.G_buff, self.b_vectors_buff,
                np.int32(N_points), np.int32(N_G), np.int32(self.n_harm),
                phi_buff
            )
        else:
            self.prg.eval_potential_vacuum(
                self.queue, (N_points,), None,
                eval_buff, self.C_G_buff, self.G_buff, self.b_vectors_buff,
                np.int32(N_points), np.int32(N_G), np.int32(self.n_harm),
                phi_buff
            )
        
        cl.enqueue_copy(self.queue, phi_out, phi_buff)
        self.queue.finish()
        t1 = time.time()
        
        print(f"Vacuum evaluation: {N_points} points in {t1-t0:.3f} s ({N_points/(t1-t0):.0f} pts/s)")
        
        return phi_out.reshape(X.shape)
    
    def eval_full(self, X, Y, Z):
        """
        Evaluate potential at arbitrary 3D positions.
        Valid for any z (inside or outside slab).
        
        Parameters:
            X, Y, Z: 2D arrays - coordinates (e.g., from meshgrid)
            
        Returns:
            phi: 2D array - potential in eV/e
        """
        if self.w is None:
            raise ValueError("System not prepared. Call prepare_system() first.")
        
        # Flatten evaluation points
        x_flat = X.ravel().astype(np.float32)
        y_flat = Y.ravel().astype(np.float32)
        z_flat = Z.ravel().astype(np.float32)
        N_points = len(x_flat)
        
        eval_points = np.column_stack([x_flat, y_flat, z_flat,
                                      np.zeros(N_points, dtype=np.float32)])
        
        # Allocate buffers
        mf = cl.mem_flags
        eval_buff = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR,
                              hostbuf=eval_points)
        phi_out = np.empty(N_points, dtype=np.float32)
        phi_buff = cl.Buffer(self.ctx, mf.WRITE_ONLY, phi_out.nbytes)
        
        # Launch kernel
        N_G = len(self.G_data)
        N_ions = len(self.ion_data)
        
        t0 = time.time()
        self.prg.eval_potential_full(
            self.queue, (N_points,), None,
            eval_buff, self.w_buff, self.ion_buff, self.G_buff,
            self.b_vectors_buff, np.float32(self.area),
            np.int32(N_points), np.int32(N_ions), np.int32(N_G),
            phi_buff
        )
        
        cl.enqueue_copy(self.queue, phi_out, phi_buff)
        self.queue.finish()
        t1 = time.time()
        
        print(f"Full evaluation: {N_points} points in {t1-t0:.3f} s ({N_points/(t1-t0):.0f} pts/s)")
        
        return phi_out.reshape(X.shape)
    
    def eval_brute(self, X, Y, Z, N_rep=20):
        """
        Brute force Coulomb sum over periodic images.
        Slow but exact (within N_rep shells).
        
        Parameters:
            X, Y, Z: 2D arrays - evaluation coordinates
            N_rep: int - number of PBC shells
            
        Returns:
            phi: 2D array - potential in eV/e
        """
        if self.ion_data is None:
            raise ValueError("System not prepared. Call prepare_system() first.")
        
        # Flatten evaluation points
        x_flat = X.ravel().astype(np.float32)
        y_flat = Y.ravel().astype(np.float32)
        z_flat = Z.ravel().astype(np.float32)
        N_points = len(x_flat)
        
        eval_points = np.column_stack([x_flat, y_flat, z_flat,
                                      np.zeros(N_points, dtype=np.float32)])
        
        # Allocate buffers
        mf = cl.mem_flags
        eval_buff = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR,
                              hostbuf=eval_points)
        
        a_vec_2d = self.a_vec.reshape(1, 2).astype(np.float32)
        b_vec_2d = self.b_vec.reshape(1, 2).astype(np.float32)
        a_buff = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR,
                           hostbuf=a_vec_2d)
        b_buff = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR,
                           hostbuf=b_vec_2d)
        
        phi_out = np.empty(N_points, dtype=np.float32)
        phi_buff = cl.Buffer(self.ctx, mf.WRITE_ONLY, phi_out.nbytes)
        
        # Launch kernel
        N_ions = len(self.ion_data)
        
        t0 = time.time()
        self.prg.eval_potential_brute(
            self.queue, (N_points,), None,
            eval_buff, self.ion_buff, a_buff, b_buff,
            np.int32(N_points), np.int32(N_ions), np.int32(N_rep),
            phi_buff
        )
        
        cl.enqueue_copy(self.queue, phi_out, phi_buff)
        self.queue.finish()
        t1 = time.time()
        
        print(f"Brute evaluation: {N_points} points in {t1-t0:.3f} s ({N_points/(t1-t0):.0f} pts/s)")
        
        return phi_out.reshape(X.shape)


def test_ewald_cl():
    """
    Test OpenCL Ewald implementation against Python reference.
    """
    print("="*60)
    print("Testing OpenCL Ewald2D")
    print("="*60)
    
    # Import Python reference
    import sys
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../..'))
    from pyBall.Ewald2D import Ewald2D
    
    # Create NaCl test system
    # NaCl lattice: a = 4.0 Å, ions at (0,0) and (2,2) in Å
    a = 4.0
    a_vec = np.array([a, 0.0])
    b_vec = np.array([0.0, a])
    
    # Ion positions (NaCl 1x1 unit cell)
    rx = np.array([0.0, 2.0])
    ry = np.array([0.0, 2.0])
    rz = np.array([0.0, 0.0])
    q = np.array([1.0, -1.0])
    
    ion_data = np.column_stack([rx, ry, rz, q])
    
    print(f"\nTest system:")
    print(f"  Lattice: a=({a_vec[0]:.1f}, {a_vec[1]:.1f}), b=({b_vec[0]:.1f}, {b_vec[1]:.1f})")
    print(f"  Ions: Na at (0,0,0) q=+1, Cl at (2,2,0) q=-1")
    
    # Initialize OpenCL Ewald
    print("\nInitializing OpenCL...")
    ew_cl = SurfaceEwaldCL()
    ew_cl.prepare_system(ion_data, a_vec, b_vec, n_harm=3)
    
    # Initialize Python reference
    print("\nInitializing Python reference...")
    ew_py = Ewald2D(a_vec, b_vec, rx, ry, rz, q, n_harm=3)
    
    # Test 1: Vacuum evaluation at z = 2.0
    print("\n" + "="*60)
    print("Test 1: Vacuum evaluation at z = 2.0 Å")
    print("="*60)
    
    xv = np.linspace(0, a, 20)
    yv = np.linspace(0, a, 20)
    X, Y = np.meshgrid(xv, yv)
    z_test = 2.0
    
    # OpenCL
    phi_cl = ew_cl.eval_vacuum(X, Y, z_test)
    
    # Python
    phi_py = ew_py.phi_vacuum_xy(X, Y, z_test)
    
    # Compare
    diff = phi_cl - phi_py
    rmse = np.sqrt(np.mean(diff**2))
    max_err = np.max(np.abs(diff))
    
    print(f"  RMSE: {rmse:.6e} eV")
    print(f"  Max error: {max_err:.6e} eV")
    
    if rmse < 1e-5:
        print("  ✓ PASS: Agreement within 1e-5 eV")
    else:
        print("  ✗ FAIL: Agreement worse than 1e-5 eV")
    
    # Test 2: Full evaluation at various z
    print("\n" + "="*60)
    print("Test 2: Full evaluation (1D line scan)")
    print("="*60)
    
    x0, y0 = 0.5, 0.5
    z_arr = np.linspace(-1.0, 5.0, 100)
    
    # OpenCL
    X_line = np.full((1, len(z_arr)), x0, dtype=np.float32)
    Y_line = np.full((1, len(z_arr)), y0, dtype=np.float32)
    Z_line = z_arr.reshape(1, -1).astype(np.float32)
    
    phi_cl_line = ew_cl.eval_full(X_line, Y_line, Z_line)[0, :]
    
    # Python
    phi_py_line = ew_py.phi_full_1d(x0, y0, z_arr)
    
    # Compare
    diff = phi_cl_line - phi_py_line
    rmse = np.sqrt(np.mean(diff**2))
    max_err = np.max(np.abs(diff))
    
    print(f"  RMSE: {rmse:.6e} eV")
    print(f"  Max error: {max_err:.6e} eV")
    
    if rmse < 1e-5:
        print("  ✓ PASS: Agreement within 1e-5 eV")
    else:
        print("  ✗ FAIL: Agreement worse than 1e-5 eV")
    
    # Test 3: Brute force (small test)
    print("\n" + "="*60)
    print("Test 3: Brute force validation")
    print("="*60)
    
    x_test = np.array([0.5], dtype=np.float32)
    y_test = np.array([0.5], dtype=np.float32)
    z_test_arr = np.array([2.0], dtype=np.float32)
    X_test = x_test.reshape(1, 1)
    Y_test = y_test.reshape(1, 1)
    Z_test = z_test_arr.reshape(1, 1)
    
    phi_brute_cl = ew_cl.eval_brute(X_test, Y_test, Z_test, N_rep=10)[0, 0]
    phi_brute_py = ew_py.phi_brute_1d(x_test[0], y_test[0], z_test_arr, N_rep=10)[0]
    
    diff_brute = abs(phi_brute_cl - phi_brute_py)
    print(f"  OpenCL brute: {phi_brute_cl:.6f} eV")
    print(f"  Python brute: {phi_brute_py:.6f} eV")
    print(f"  Difference: {diff_brute:.6e} eV")
    
    if diff_brute < 1e-4:
        print("  ✓ PASS: Brute force agreement within 1e-4 eV")
    else:
        print("  ✗ FAIL: Brute force agreement worse than 1e-4 eV")
    
    print("\n" + "="*60)
    print("All tests complete")
    print("="*60)
    
    return ew_cl, ew_py


if __name__ == "__main__":
    test_ewald_cl()
