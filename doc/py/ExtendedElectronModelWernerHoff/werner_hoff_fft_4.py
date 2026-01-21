import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import art3d
from skimage import measure

def plot_intersecting_nodal_planes():
    # Grid Setup
    N = 60
    L = 4.0
    x = np.linspace(-L, L, N)
    X, Y, Z = np.meshgrid(x, x, x)

    # --- Scenario 1: Acetylene-like (Cylindrical Symmetry) ---
    # Phi_1 = p_x (Dumbbell along X) -> Node is Plane X=0
    # Phi_2 = p_y (Dumbbell along Y) -> Node is Plane Y=0
    
    # Let's add some "molecular" character (exponential decay)
    R = np.sqrt(X**2 + Y**2 + Z**2)
    radial = np.exp(-R)
    
    phi_1 = X * radial  # Node at X=0
    phi_2 = Y * radial  # Node at Y=0
    
    # Construct Complex Wavefunction (The Current-Carrying State)
    # Psi = px + i*py
    psi = phi_1 + 1j * phi_2
    rho = np.abs(psi)**2

    fig = plt.figure(figsize=(12, 6))
    
    # --- Visualization ---
    ax = fig.add_subplot(1, 2, 1, projection='3d')
    ax.set_title("Acetylene Model:\nIntersection of px and py Nodes")

    # 1. Plot Node of Phi_1 (Red Sheet) - The X=0 Plane
    # We find isosurface where phi_1 approx 0
    try:
        verts1, faces1, _, _ = measure.marching_cubes(phi_1, level=0.0)
        # Scale to match plot
        verts1 = verts1 * (2*L/N) - L
        mesh1 = art3d.Poly3DCollection(verts1[faces1], alpha=0.3, color='red', label='Node px')
        ax.add_collection3d(mesh1)
    except: pass

    # 2. Plot Node of Phi_2 (Blue Sheet) - The Y=0 Plane
    try:
        verts2, faces2, _, _ = measure.marching_cubes(phi_2, level=0.0)
        verts2 = verts2 * (2*L/N) - L
        mesh2 = art3d.Poly3DCollection(verts2[faces2], alpha=0.3, color='blue', label='Node py')
        ax.add_collection3d(mesh2)
    except: pass
    
    # 3. Plot The Vortex Core (Total Density ~ 0)
    # This should be the LINE where Red and Blue cross.
    try:
        verts_core, faces_core, _, _ = measure.marching_cubes(rho, level=0.001) # Small epsilon
        verts_core = verts_core * (2*L/N) - L
        mesh_core = art3d.Poly3DCollection(verts_core[faces_core], alpha=1.0, color='black')
        ax.add_collection3d(mesh_core)
    except: pass

    ax.set_xlim(-2, 2); ax.set_ylim(-2, 2); ax.set_zlim(-2, 2)
    ax.set_xlabel("X"); ax.set_ylabel("Y"); ax.set_zlabel("Z")

    # --- Scenario 2: Distorted/Twisted (Ethylene-like or General) ---
    # Phi_1 = p_x (Standard Node X=0)
    # Phi_2 = p_z (Node Z=0) -> Imagine mixing px and pz
    
    ax2 = fig.add_subplot(1, 2, 2, projection='3d')
    ax2.set_title("Mixed Symmetry:\nIntersection of px and pz")
    
    phi_A = X * radial
    phi_B = Z * radial # Orthogonal direction
    
    rho_mix = np.abs(phi_A + 1j * phi_B)**2

    # Plot Node A
    try:
        vertsA, facesA, _, _ = measure.marching_cubes(phi_A, level=0.0)
        vertsA = vertsA * (2*L/N) - L
        meshA = art3d.Poly3DCollection(vertsA[facesA], alpha=0.3, color='red')
        ax2.add_collection3d(meshA)
    except: pass

    # Plot Node B
    try:
        vertsB, facesB, _, _ = measure.marching_cubes(phi_B, level=0.0)
        vertsB = vertsB * (2*L/N) - L
        meshB = art3d.Poly3DCollection(vertsB[facesB], alpha=0.3, color='blue')
        ax2.add_collection3d(meshB)
    except: pass
    
    # Plot Intersection (Vortex)
    try:
        vertsC, facesC, _, _ = measure.marching_cubes(rho_mix, level=0.001)
        vertsC = vertsC * (2*L/N) - L
        meshC = art3d.Poly3DCollection(vertsC[facesC], alpha=1.0, color='black')
        ax2.add_collection3d(meshC)
    except: pass

    ax2.set_xlim(-2, 2); ax2.set_ylim(-2, 2); ax2.set_zlim(-2, 2)

    plt.show()

plot_intersecting_nodal_planes()