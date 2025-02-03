
Q0 = 0.7

a = 5.6413 / 2  # a = 2.82065 Å

#nx,ny,nz = 13,12,3
nx,ny,nz = 15,8,3

#positions = []
#symbols   = []

Lx = a*nx
LxStep = Lx*0.5
ax = -1.0/nx

cell = [
    [nx * a, 0, 0],  # First cell vector (x-axis)
    [0, ny * a, 0],  # Second cell vector (y-axis)
    [0, 0, nz * a]   # Third cell vector (z-axis)
]




fout = open('NaCl_cubic.xyz','w')
fout.write(f'{nx*ny*nz}\n')
fout.write( f'lvs {cell[0][0]} {cell[0][1]} {cell[0][2]}   {cell[1][0]} {cell[1][1]} {cell[1][2]}   {cell[2][0]} {cell[2][1]} {cell[2][2]}\n' )
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
            
            #positions.append([x, y, z])
            # Alternate between Na and Cl based on the sum of indices
            if i % 2 == 0:
                S = 'Na'
                Q = +Q0
            else:
                S = 'Cl'
                Q = -Q0
            
            fout.write( f'{S} {x:10.5f} {y:10.5f} {z:10.5f} {Q} \n')    
            #symbols.append(S)

fout.close()

#from ase import Atoms
#import numpy as np
#nacl.set_cell(cell)
#nacl = Atoms(symbols=symbols, positions=positions)
#nacl.center(vacuum=5.0, axis=2)
#nacl.write('NaCl_simple_cubic.xyz')
