import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

# ==========================================
# 1. SETUP: Atoms on a Substrate (y=0)
# ==========================================
# (x, y, weight)
atoms_data = [
    (-0.7, 0.5, 0.25), 
    ( 0.5, 0.4, 0.15)
]

# Method of Images for a perfectly flat Substrate at y=0
sources = []
for (ax, ay, r) in atoms_data:
    sources.append((ax + 1j*ay, r))     # Real atom
    sources.append((ax - 1j*ay, r))     # Image atom

A = 1.5   # Strength of the vertical gradient (Height factor)
freq_u = 6.0 # Density of Horizontal "Height" shells
freq_v = 6.0 # Density of Vertical "Streamline" columns

# ==========================================
# 2. THE MATHEMATICAL FIELD (W = U + iV)
# ==========================================
def get_W(z):
    """
    W = -i*A*z + sum(R*log(z-z_i))
    U (Real) = Ay + Sum log|z-zi|  -> Controls HEIGHT/SHELLS
    V (Imag) = -Ax + Sum arg(z-zi) -> Controls COLUMN/RAYS
    """
    # Linear part rotated to be vertical
    w = -1j * A * z
    for pos, r in sources:
        w += r * np.log(z - pos)
    return w

def get_dW(z):
    """Complex Derivative for Newton Solver"""
    dw = -1j * A
    for pos, r in sources:
        dw += r / (z - pos)
    return dw

# ==========================================
# 3. ROBUST INVERSION (Finding the Red Dots)
# ==========================================
def find_grid_intersection(u_target, v_target):
    """
    Finds the Z point for a specific U, V.
    Uses a far-field guess to ensure we stay in the upper half-plane.
    """
    # Far field guess: Invert the linear part w = -iAz  => z = w / -iA
    z = (u_target + 1j * v_target) / (-1j * A)
    
    # Newton-Raphson refinement
    for _ in range(10):
        w_curr = get_W(z)
        dw = get_dW(z)
        step = (w_curr - (u_target + 1j * v_target)) / dw
        if np.abs(step) < 1e-9: break
        z -= step
    return z

# ==========================================
# 4. DATA GENERATION
# ==========================================
res = 600
x_range = np.linspace(-2, 2, res)
y_range = np.linspace(0.0, 2, res)
X, Y = np.meshgrid(x_range, y_range)
Z_pixels = X + 1j * Y

W_pix = get_W(Z_pixels)
U = np.real(W_pix)
V = np.imag(W_pix)

# Intersection points (Red Dots)
u_targets = np.arange(0.2, 3.5, 1.0/freq_u) # Horizontal tiers
v_targets = np.arange(-3.5, 3.5, 1.0/freq_v) # Vertical columns

sample_x, sample_y = [], []
for u in u_targets:
    for v in v_targets:
        z_root = find_grid_intersection(u, v)
        # Keep only points that are actually above the substrate and outside atoms
        if z_root.imag > 0.001:
            outside = True
            for ax, ay, r in atoms_data:
                if np.abs(z_root - (ax + 1j*ay)) < r:
                    outside = False; break
            if outside:
                sample_x.append(z_root.real)
                sample_y.append(z_root.imag)

# ==========================================
# 5. VISUALIZATION
# ==========================================
fig, ax = plt.subplots(figsize=(12, 8))

# A. Seamless Checkerboard
# We use sin() of the potential to create a smooth, branch-cut-free pattern
checker = (np.sin(U * freq_u * np.pi) > 0) ^ (np.sin(V * freq_v * np.pi) > 0)
ax.imshow(checker, extent=[-2, 2, 0, 2], origin='lower', cmap='binary', alpha=0.1)

# B. Grid Lines (The Manifolds)
# Cyan = Streamlines (Columns), Magenta = Isolines (Tiers)
ax.contour(X, Y, U, levels=u_targets, colors='magenta', linewidths=0.6, alpha=0.4)
ax.contour(X, Y, V, levels=v_targets, colors='cyan', linewidths=0.6, alpha=0.4)

# C. The Substrate (y=0)
ax.axhline(0, color='#1e293b', linewidth=5, zorder=5)
ax.fill_between([-2, 2], -0.1, 0, color='#f1f5f9', zorder=4)

# D. The Atoms
for ax_p, ay_p, r in atoms_data:
    circle = plt.Circle((ax_p, ay_p), r, color='#64748b', alpha=0.5, zorder=6)
    ax.add_patch(circle)
    ax.plot(ax_p, ay_p, 'k+', markersize=10, zorder=7)

# E. Intersection Samples (The Red Dots)
ax.scatter(sample_x, sample_y, s=15, color='#ef4444', edgecolors='white', linewidth=0.5, zorder=10)

ax.set_aspect('equal')
ax.set_xlim([-1.8, 1.8])
ax.set_ylim([-0.05, 1.8])
ax.set_axis_off()
ax.set_title("Corrected Substrate-Conformal Grid\nVertical Potential (Magenta) | Transverse Streamlines (Cyan)", fontsize=14)

plt.tight_layout()
plt.show()