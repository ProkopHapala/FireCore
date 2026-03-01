import numpy as np
import matplotlib.pyplot as plt

# ==========================================
# 1. Setup: atoms on surface
# ==========================================
atoms_data = [
    ( 0.0,  0.5, -0.05), 
    ( 1.0, 0.5, -0.05)
]

A = 1.2          # free-stream strength
m = 3            # decay power (1,2,3,...)
alpha = 0.05     # pole strength

# Method of images
sources = []
for (ax, ay, r) in atoms_data:
    sources.append(ax + 1j*ay)
    sources.append(ax + 1j*ay)

# ==========================================
# 2. Holomorphic Field
# ==========================================
def get_W(z):
    w = -1j*A*z
    for pos in sources:
        w += alpha / (z - pos)**m
    return w

def get_dW(z):
    dw = -1j*A
    for pos in sources:
        dw -= alpha*m / (z - pos)**(m+1)
    return dw

# ==========================================
# 3. Robust inversion
# ==========================================
def find_root(u_target, v_target):
    z = (u_target + 1j*v_target)/(-1j*A)
    for _ in range(20):
        f = get_W(z) - (u_target + 1j*v_target)
        df = get_dW(z)
        step = f/df
        z -= step
        if abs(step) < 1e-10:
            break
    return z

# ==========================================
# 4. Generate background field
# ==========================================
res = 600
x = np.linspace(-2,2,res)
y = np.linspace(0.0,2,res)
X,Y = np.meshgrid(x,y)
Z = X + 1j*Y

W = get_W(Z)
U = np.real(W)
V = np.imag(W)

# Targets
freq = 3.0
u_targets = np.arange(0.3,3.5,1/freq)
v_targets = np.arange(-3,3,1/freq)

# Sample points
sx, sy = [], []
for u in u_targets:
    for v in v_targets:
        z = find_root(u,v)
        if z.imag > 0:
            sx.append(z.real)
            sy.append(z.imag)

# ==========================================
# 5. Plot
# ==========================================
fig, ax = plt.subplots(figsize=(12,8))

ax.contour(X,Y,U,levels=u_targets,colors='magenta',linewidths=0.5)
ax.contour(X,Y,V,levels=v_targets,colors='cyan',linewidths=0.5)

ax.scatter(sx,sy,s=15,color='red')

# substrate
ax.axhline(0,color='black',linewidth=4)

# atoms
for ax_p, ay_p, r in atoms_data:
    circle = plt.Circle((ax_p,ay_p), r, color='gray', alpha=0.6)
    ax.add_patch(circle)

ax.set_aspect('equal')
ax.set_xlim([-1.8,1.8])
ax.set_ylim([0,1.8])
plt.show()