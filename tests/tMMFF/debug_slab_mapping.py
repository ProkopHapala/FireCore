#!/usr/bin/env python3
import sys
import numpy as np
import pyopencl as cl
import os
sys.path.append('/home/indranil/git/FireCore')
from pyBall.OCL.GridFF import GridFF_cl
from pyBall.OCL.clUtils import GridShape

def main():
    """Test the slabPotential kernel with a known pattern to identify shifts"""
    print("=== Testing slabPotential kernel mapping ===")
    
    # Create a small test grid with a known pattern
    # Each position will have a unique value: z*10000 + y*100 + x
    nz, ny, nx = 10, 8, 6
    print(f"Creating test grid with shape: ({nz}, {ny}, {nx})")
    
    # Create test pattern in ZYX order
    host = np.zeros((nz, ny, nx), dtype=np.float32)
    for z in range(nz):
        for y in range(ny):
            for x in range(nx):
                host[z, y, x] = z*10000 + y*100 + x
    
    print("Input grid pattern: V[z,y,x] = z*10000 + y*100 + x")
    print(f"Example: V[1,2,3] = {host[1,2,3]}")
    
    # Initialize GridFF_cl with proper grid shape
    grid = GridFF_cl()
    
    # Create a GridShape object
    lvec = np.array([
        [nx * 0.1, 0, 0],
        [0, ny * 0.1, 0],
        [0, 0, nz * 0.1]
    ], dtype=np.float32)
    
    gsh = GridShape()
    gsh.init_from_lvec(lvec, (nx, ny, nz), origin=(-0.3, -0.4, -0.5))
    grid.set_grid(gsh)
    
    # Upload to GPU
    mf = cl.mem_flags
    buf_in = cl.Buffer(grid.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=host)
    
    # Run slabPotential
    print("\nRunning slabPotential kernel...")
    out = grid.slabPotential(buf_in, nz_slab=nz, dipol=0.0, bDownload=True)
    
    print(f"Output grid shape: {out.shape}")
    
    # Expected output: V_out[z,y,x] = V_in[nz-1-z, ny-1-y, nx-1-x]
    ref = np.zeros_like(host)
    for z in range(nz):
        for y in range(ny):
            for x in range(nx):
                ref[z, y, x] = host[nz-1-z, ny-1-y, nx-1-x]
    
    # Compare
    diff = out - ref
    max_diff = np.abs(diff).max()
    print(f"\nMax difference between expected and actual: {max_diff}")
    
    if max_diff > 1e-6:
        print("MISMATCH DETECTED!")
        # Find first few mismatches
        bad = np.where(np.abs(diff) > 1e-6)
        print("First few mismatches at (z,y,x):")
        for i in range(min(5, bad[0].size)):
            z, y, x = bad[0][i], bad[1][i], bad[2][i]
            print(f"  Position ({z},{y},{x}):")
            print(f"    Expected: {ref[z,y,x]}")
            print(f"    Actual:   {out[z,y,x]}")
            print(f"    Input at ({nz-1-z},{ny-1-y},{nx-1-x}): {host[nz-1-z, ny-1-y, nx-1-x]}")
            
            # Try to figure out where the value came from
            found = False
            for sz in range(nz):
                for sy in range(ny):
                    for sx in range(nx):
                        if abs(host[sz, sy, sx] - out[z, y, x]) < 1e-6:
                            print(f"    Output value matches input at ({sz},{sy},{sx})")
                            print(f"    Shift: ({sz-(nz-1-z)}, {sy-(ny-1-y)}, {sx-(nx-1-x)})")
                            found = True
                            break
                    if found: break
                if found: break
    else:
        print("Perfect match! No shift detected in slabPotential kernel.")

if __name__ == "__main__":
    main()
