from ase import Atoms
import numpy as np

a = 5.6413 / 2  # a = 2.82065 Å

nx,ny,nz = 13,12,3

positions = []
symbols   = []


Lx = a*nx
LxStep = Lx*0.5
ax = -1.0/nx


# Generate atomic positions
for iz in range(nz):
    for ix in range(nx):
        for iy in range(ny):
            # Calculate the position of the atom
            x = ix * a
            y = iy * a
            z = iz * a
            
            i = ix + iy + iz
            
            z_=z+ax*x
            x-=ax*z
            z=z_
            if x>LxStep:
               z+=a
               i+=1
            
            positions.append([x, y, z])
            # Alternate between Na and Cl based on the sum of indices
            if i % 2 == 0:
                symbols.append('Na')
            else:
                symbols.append('Cl')

# Create the ASE Atoms object
nacl = Atoms(symbols=symbols, positions=positions)

# Set the cell vectors
cell = [
    [nx * a, 0, 0],  # First cell vector (x-axis)
    [0, ny * a, 0],  # Second cell vector (y-axis)
    [0, 0, nz * a]   # Third cell vector (z-axis)
]
nacl.set_cell(cell)

# Center the structure and add vacuum along the z-axis
nacl.center(vacuum=5.0, axis=2)

# Save the structure to a file
nacl.write('NaCl_simple_cubic.xyz')
