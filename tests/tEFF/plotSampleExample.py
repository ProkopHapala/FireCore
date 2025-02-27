import numpy as np
import matplotlib.pyplot as plt

# Define a vector function (example: E, fx, fy, fz)
def example_vector_function( points, w, freq=(1.0, 2.0, 3.0), amp=0.05 ):
    # deserialize input coordinates
    x = points[:,0]
    y = points[:,1]
    z = points[:,2]
    # prepare output array
    FEout = np.zeros( ( len(points), 4) )
    # our example function is  Gaussian  + Amp*sin(freq*p)
    gaussian = np.exp(-(x**2 + y**2 + z**2) / (w**2) )
    E  = gaussian          #+ amp*np.sin( freq[0]*x) + np.sin(  freq[1]*y) + np.sin(  freq[2]*z )
    fx = gaussian*x*2/w**2 #+ amp*np.cos( freq[0]*x ) *freq[0] * np.sin( freq[1]*y) + np.sin(  freq[2]*z )
    fy = gaussian*y*2/w**2 #+ amp*np.cos( freq[1]*y ) *freq[1] * np.sin( freq[0]*x) + np.sin(  freq[2]*z )
    fz = gaussian*z*2/w**2 #+ amp*np.cos( freq[2]*z ) *freq[2] * np.sin( freq[0]*x) + np.sin(  freq[1]*y )
    FEout[:,0] = fx
    FEout[:,1] = fy
    FEout[:,2] = fz
    FEout[:,3] = E
    return FEout

# ========== Evaluate the function on a grid ==========

nx, ny = 100, 100
x = np.linspace(-5, 5, nx)
y = np.linspace(-5, 5, ny)
X, Y = np.meshgrid(x, y)
z0= 1.0

# Serialize the 2D coordinates into an array of shape (nx, ny, 3)
points = np.zeros( (nx*ny, 3) )
points[:,0] = X.flatten()
points[:,1] = Y.flatten()
points[:,2] = z0

# Apply the vector function to each point in the grid
FEout = example_vector_function( points, 3.0 )

# Deserialize the results back into 2D grids
fx_grid = FEout[:,0].reshape(nx, ny)
fy_grid = FEout[:,1].reshape(nx, ny)
fz_grid = FEout[:,2].reshape(nx, ny)
E_grid  = FEout[:,3].reshape(nx, ny)

# ========== Plotting ==========

# Plot the components using imshow
fig, axes = plt.subplots(2, 2, figsize=(10, 10))

# Plot E
ax = axes[0, 0]
im = ax.imshow(E_grid, extent=(-5, 5, -5, 5), origin='lower', cmap='viridis')
ax.set_title('E')
fig.colorbar(im, ax=ax)

# Plot fx
ax = axes[0, 1]
im = ax.imshow(fx_grid, extent=(-5, 5, -5, 5), origin='lower', cmap='plasma')
ax.set_title('fx')
fig.colorbar(im, ax=ax)

# Plot fy
ax = axes[1, 0]
im = ax.imshow(fy_grid, extent=(-5, 5, -5, 5), origin='lower', cmap='plasma')
ax.set_title('fy')
fig.colorbar(im, ax=ax)

# Plot fz
ax = axes[1, 1]
im = ax.imshow(fz_grid, extent=(-5, 5, -5, 5), origin='lower', cmap='plasma')
ax.set_title('fz')
fig.colorbar(im, ax=ax)

plt.tight_layout()
plt.show()