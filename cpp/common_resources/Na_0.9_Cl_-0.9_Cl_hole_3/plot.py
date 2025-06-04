import numpy as np
import matplotlib.pyplot as plt


# grid =np.load("Bspline_PLQd_ocl.npy")
# grid =np.load("old_Bspline_PLQd.npy")
grid =np.load("Bspline_PLQd.npy")
charge =np.load("Qgrid_ocl.npy")

print("Shape of grid:", grid.shape)  # expected (Nx, Ny, Nz, 3)

iz = 50
Ep_1=grid[:,:,iz,0] ; #print("Ep", Ep.shape)
El_1=grid[:,:,iz,1] ; #print("El", El.shape)
Eq_1=grid[:,:,iz,2] ; #print("Eq", Eq.shape)
# Eq_rot_1 = np.rot90(Eq_1,k=3)
Eq_rot_1 = Eq_1.transpose()

# Eq_rot_1 = Eq_1



ix = 0
Ep_2=grid[ix,:,:,0] ; 
El_2=grid[ix,:,:,1] ;
Eq_2=grid[ix,:,:,2] ; 
# Eq_rot_2 = np.rot90(Eq_2)
Eq_rot_2 = Eq_2

iy = 40
Ep_3=grid[:,iy,:,0] ; 
El_3=grid[:,iy,:,1] ;
Eq_3=grid[:,iy,:,2] ; 
# Eq_rot_3 = np.rot90(Eq_3)
Eq_rot_3 = Eq_3







'''
Transpose to the Morse and np.rot90 to the coulomb can make them align, for the case when the origin of grid is set to (0,0,0) all these operations are on top of the saved numpy array Bspline_PLQd.npy 

For the other case, when the origin of grid is set to (-Lx/2,-Ly/2,0) then same is working transpose to morse and np.rot90 to coulomb working fine .

However, abot the z axis the Morse and Coulomb are set differently; the Morse is fine if g0(,,z0) and Coulomb is fine if g0(,,0) they are not at same page regarding tihs, which may be due to the trancation in coulomb  

Due to this reason when the grid was computed for atoms having -ve coordinates shows problem in coulomb 


'''
x_points = np.array([160, 160, 160, 200])
y_points = np.array([200, 240, 280, 200])

x_points_shift = x_points -160
y_points_shift = y_points -160 


plt.figure(figsize=(15,5))
# plt.subplot(1,3,1); plt.imshow( Ep.transpose(), origin="lower" ); plt.colorbar(); plt.scatter(x_points, y_points, marker='o', color='red')
# plt.subplot(1,3,2); plt.imshow( El.transpose(), origin="lower" ); plt.colorbar(); plt.scatter(x_points, y_points, marker='o', color='red')
# plt.subplot(1,3,3); plt.imshow( Eq_rot, origin="lower" ); plt.colorbar(); plt.scatter(x_points, y_points, marker='o', color='red')

plt.subplot(1,3,1); plt.imshow( Ep_1.transpose(), origin="lower" ); plt.colorbar(); plt.scatter(x_points_shift, y_points_shift, marker='o', color='red')
plt.subplot(1,3,2); plt.imshow( El_1.transpose(), origin="lower" ); plt.colorbar(); plt.scatter(x_points_shift, y_points_shift, marker='o', color='red')
plt.subplot(1,3,3); plt.imshow( Eq_rot_1, origin="lower" ); plt.colorbar(); plt.scatter(x_points_shift, y_points_shift, marker='o', color='red')

plt.figure(figsize=(15,5))
# plt.subplot(1,3,1); plt.imshow( Ep_2.transpose(), origin="lower" ); plt.colorbar(); plt.scatter(x_points, y_points, marker='o', color='red')
# plt.subplot(1,3,2); plt.imshow( El_2.transpose(), origin="lower" ); plt.colorbar(); plt.scatter(x_points, y_points, marker='o', color='red')
# plt.subplot(1,3,3); plt.imshow( Eq_rot_2, origin="lower" ); plt.colorbar(); plt.scatter(x_points, y_points, marker='o', color='red')
plt.subplot(1,3,1); plt.imshow( Ep_2.transpose(), origin="lower" ); plt.colorbar(); plt.scatter(x_points_shift, y_points_shift, marker='o', color='red')
plt.subplot(1,3,2); plt.imshow( El_2.transpose(), origin="lower" ); plt.colorbar(); plt.scatter(x_points_shift, y_points_shift, marker='o', color='red')
plt.subplot(1,3,3); plt.imshow( Eq_rot_2, origin="lower" ); plt.colorbar(); plt.scatter(x_points_shift, y_points_shift, marker='o', color='red')

plt.figure(figsize=(15,5))
# plt.subplot(1,3,1); plt.imshow( Ep_3.transpose(), origin="lower" ); plt.colorbar(); plt.scatter(x_points, y_points, marker='o', color='red')
# plt.subplot(1,3,2); plt.imshow( El_3.transpose(), origin="lower" ); plt.colorbar(); plt.scatter(x_points, y_points, marker='o', color='red')
# plt.subplot(1,3,3); plt.imshow( Eq_rot_3, origin="lower" ); plt.colorbar(); plt.scatter(x_points, y_points, marker='o', color='red')
plt.subplot(1,3,1); plt.imshow( Ep_3.transpose(), origin="lower" ); plt.colorbar(); plt.scatter(x_points_shift, y_points_shift, marker='o', color='red')
plt.subplot(1,3,2); plt.imshow( El_3.transpose(), origin="lower" ); plt.colorbar(); plt.scatter(x_points_shift, y_points_shift, marker='o', color='red')
plt.subplot(1,3,3); plt.imshow( Eq_rot_3, origin="lower" ); plt.colorbar(); plt.scatter(x_points_shift, y_points_shift, marker='o', color='red')


plt.figure(figsize=(15,5))
print("Shape of charge:", charge.shape) 
iz1=2
iy1=40
ix1=0
plt.subplot(1,3,1); plt.imshow(charge[iz1,:,:], origin="lower"); plt.colorbar(); plt.title(f'XY plane (z={iz1})')
plt.subplot(1,3,2); plt.imshow(charge[:,iy1,:], origin="lower"); plt.colorbar(); plt.title(f'XZ plane (y={iy1})')
plt.subplot(1,3,3); plt.imshow(charge[:,:,ix1], origin="lower"); plt.colorbar(); plt.title(f'YZ plane (x={ix1})')


charge_sum_z = np.sum(charge, axis=0)  # Sum along z-axis (axis=0)

plt.figure(figsize=(10,8))
plt.imshow(charge_sum_z, origin="lower")
plt.colorbar(label='Summed Charge Density')
plt.title('Sum of all XY slices along z-axis')
plt.xlabel('X')
plt.ylabel('Y')
# plt.scatter(x_points, y_points, marker='o', color='red')

plt.show()
