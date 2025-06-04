import numpy as np
import matplotlib.pyplot as plt


# grid =np.load("Bspline_PLQd_ocl.npy")
# grid =np.load("old_Bspline_PLQd.npy")
grid =np.load("Bspline_PLQd.npy")

print("Shape of grid:", grid.shape)  # expected (Nx, Ny, Nz, 3)

iz=100
Ep=grid[:,:,iz,0] ; print("Ep", Ep.shape)
El=grid[:,:,iz,1] ; print("El", El.shape)
Eq=grid[:,:,iz,2] ; print("Eq", Eq.shape)

plt.figure(figsize=(15,5))
plt.subplot(1,3,1); plt.imshow( Ep, origin="lower"  ); plt.colorbar()
plt.subplot(1,3,2); plt.imshow( El, origin="lower" ); plt.colorbar()
plt.subplot(1,3,3); plt.imshow( Eq, origin="lower" ); plt.colorbar()

plt.show()