import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import hsv_to_rgb

# ==========================================
# 1. DEFINE THE MOLECULE AND SDF
# ==========================================
# List of atoms: (x, y, radius)
atoms = [
    (0.0, 0.5, 0.4),   # Top atom
    (-0.5, -0.3, 0.3), # Bottom left
    (0.5, -0.3, 0.3)   # Bottom right
]

def smooth_min(dists, k=10.0):
    """
    LogSumExp smooth minimum. 
    Standard exact min() creates sharp Voronoi seams between atoms.
    Smooth min blends the fields naturally, giving continuous differentiable gradients.
    """
    # Exponentiate, sum, and log. 
    res = -np.log(np.sum(np.exp(-k * dists), axis=0)) / k
    return res

def compute_sdf(X, Y, atoms):
    """
    Computes the Signed Distance Field. 
    SDF(r) = min_i ( |r - r_i| - R_i )
    < 0 inside atoms, 0 at boundary, > 0 outside.
    """
    dists = []
    for (ax, ay, r) in atoms:
        dist = np.sqrt((X - ax)**2 + (Y - ay)**2) - r
        dists.append(dist)
    return smooth_min(np.array(dists), k=8.0)

# ==========================================
# 2. THE GEOMETRIC WARPING MATH
# ==========================================
# To warp a texture in a shader or plot, we map Physical Space (X,Y) -> Texture Space (U,V).
# If we want the physical space to look DENSE, the U,V coordinates must change RAPIDLY.
#
# Math: U(X) = X + grad(SDF) * f(SDF)
# We want the derivative of f(SDF) to be high near SDF=0, and 0 far away.
# f(SDF) = A * tanh(k * SDF) is perfect. 
# Its derivative is A*k at the surface (high density) and decays to 0 far away (uniform grid).

def warp_space(X, Y, sdf, grid_extent):
    # Compute Gradient of SDF using finite differences (vector calculus equivalent to Normal)
    dx = X[0, 1] - X[0, 0]
    dy = Y[1, 0] - Y[0, 0]
    dSDF_dY, dSDF_dX = np.gradient(sdf, dy, dx)
    
    # Normalize gradient (creates a unit vector field pointing away from molecule)
    grad_mag = np.sqrt(dSDF_dX**2 + dSDF_dY**2) + 1e-8
    nX = dSDF_dX / grad_mag
    nY = dSDF_dY / grad_mag

    # The mapping function parameters
    # 'density_boost' (A*k) controls how much denser the grid is at the surface.
    # 'falloff' (k) controls how quickly it returns to normal space.
    density_boost = 1.5
    falloff = 2.0
    A = density_boost / falloff

    # Calculate U, V
    U = X + nX * (A * np.tanh(falloff * sdf))
    V = Y + nY * (A * np.tanh(falloff * sdf))
    
    return U, V

# ==========================================
# 3. GENERATE THE COMPOSITE TEXTURE
# ==========================================
def generate_texture(U, V):
    """
    Creates a combined texture: Checkerboard + Gridlines + Gradient.
    """
    freq = 4.0 # Frequency of the grid
    
    # 1. Checkerboard
    check = (np.floor(U * freq) + np.floor(V * freq)) % 2
    check = 0.4 + 0.4 * check # Remap to 0.4 - 0.8 intensity
    
    # 2. Gridlines (using sine waves for smooth anti-aliased lines)
    thickness = 0.15
    grid_x = np.abs(np.sin(U * freq * np.pi))
    grid_y = np.abs(np.sin(V * freq * np.pi))
    lines = (grid_x < thickness) | (grid_y < thickness)
    
    # 3. Gradient (Color mapped by polar angle and distance in U,V space)
    angle = np.arctan2(V, U) / (2 * np.pi) + 0.5
    radius = np.sqrt(U**2 + V**2)
    
    # Construct HSV image: Hue = angle, Saturation = 1.0, Value = Checkerboard
    hsv = np.dstack((angle, np.ones_like(angle), check))
    rgb = hsv_to_rgb(hsv)
    
    # Overlay White Gridlines
    rgb[lines] = [1.0, 1.0, 1.0]
    
    return rgb

# ==========================================
# 4. EXECUTE AND PLOT
# ==========================================
# Create a high-res grid in Physical Space
extent = [-2, 2, -2, 2]
res = 800
x_lin = np.linspace(extent[0], extent[1], res)
y_lin = np.linspace(extent[2], extent[3], res)
X, Y = np.meshgrid(x_lin, y_lin)

# 1. Compute fields
sdf = compute_sdf(X, Y, atoms)
U, V = warp_space(X, Y, sdf, extent)

# 2. Generate Textures
# We generate the unwarped (raw X,Y) to compare with the warped (U,V)
img_unwarped = generate_texture(X, Y)
img_warped = generate_texture(U, V)

# 3. Plotting
fig, axes = plt.subplots(1, 2, figsize=(14, 7))

# Plot 1: Unwarped Space (What the grid natively looks like)
axes[0].imshow(img_unwarped, extent=extent, origin='lower')
axes[0].contour(X, Y, sdf, levels=[0], colors='white', linewidths=3) # The collision radius
axes[0].set_title("Uniform Space (Raw Grid)", fontsize=14)

# Plot 2: Warped Space (Grid warped around atoms)
axes[1].imshow(img_warped, extent=extent, origin='lower')
# Overlay SDF Contours to show how the grid aligns with the distance field
contours = axes[1].contour(X, Y, sdf, levels=np.linspace(-0.5, 2.5, 8), colors='black', alpha=0.5)
axes[1].contour(X, Y, sdf, levels=[0], colors='white', linewidths=3) # Collision radius
axes[1].set_title("Warped Space (Denser near atoms)", fontsize=14)

for ax in axes:
    ax.set_aspect('equal')
    ax.set_xlim([-1.5, 1.5])
    ax.set_ylim([-1.5, 1.5])

plt.tight_layout()
plt.show()