import numpy as np

def build_grid_2d(nx, ny, m=1.0, m_end=1000.0, l=1.0, k=1.0, k_diag=-1.0):
    """Build a 2D grid of points connected by springs"""
    np_total = (nx + 1) * (ny + 1)
    masses     = np.ones(np_total) * m
    masses[0]  = m_end
    masses[nx] = m_end
    
    # Create points
    points = np.zeros((np_total, 3))
    for iy in range(ny + 1):
        for ix in range(nx + 1):
            i = iy * (nx + 1) + ix
            points[i, 0] = ix * l
            points[i, 1] = -iy * l
    
    # Create bonds
    bonds = []
    ks = []
    fixed = [0, nx]  # Fixed points (0-based indexing)
    
    for iy in range(ny + 1):
        for ix in range(nx + 1):
            i = iy * (nx + 1) + ix
            # Horizontal bonds
            if ix < nx:
                bonds.append((i, i + 1))
                ks.append(k)
            # Vertical bonds
            if iy < ny:
                bonds.append((i, i + nx + 1))
                ks.append(k)
            # Diagonal bonds
            if k_diag > 0:
                if ix < nx and iy < ny:
                    bonds.append((i, i + nx + 2))
                    ks.append(k_diag)
                if ix > 0 and iy < ny:
                    bonds.append((i, i + nx))
                    ks.append(k_diag)
    
    return np.array(bonds), points, masses, np.array(ks), fixed

def make_pd_matrix(neighbs, bonds, masses, dt, ks):
    """Create the system matrix for projective dynamics"""
    np_total = len(masses)
    A = np.zeros((np_total, np_total))
    idt2 = 1.0 / (dt * dt)
    
    # Build the system matrix
    for i in range(np_total):
        A[i, i] = masses[i] * idt2
        for ib in neighbs[i]:
            k = ks[ib]
            i_, j_ = bonds[ib]
            A[i, i] += k
            if j_ > i:
                A[i, j_] = -k
                A[j_, i] = -k
            elif i_ > i:
                A[i, i_] = -k
                A[i_, i] = -k
    
    return A, masses * idt2

def build_neighbor_list(bonds, n_points):
    """Build list of neighboring bonds for each point"""
    neighbs = [[] for _ in range(n_points)]
    for i, (i_, j_) in enumerate(bonds):
        neighbs[i_].append(i)
        neighbs[j_].append(i)
    return neighbs

def solve_pd(points, velocity, bonds, masses, ks, dt=0.1, n_iter=100, gravity=np.array([0, -9.81, 0]), fixed_points=None, call_back=None):  
    """Solve the system using projective dynamics"""
    n_points = len(points)
    neighbs  = build_neighbor_list(bonds, n_points)
    A, Mt    = make_pd_matrix(neighbs, bonds, masses, dt, ks)
    
    pos = points.copy()

    pos_prev = points - velocity * dt
    
    t = 0
    # Main simulation loop
    for itr in range(n_iter):
        # Predict step
        #pos_pred = 2 * pos - pos_prev
        
        pos_pred = pos + velocity * dt

        # Add external forces (gravity)
        f_ext = np.zeros_like(points)
        f_ext[:, 1] = gravity[1]  # Apply gravity in y direction
        
        # Solve system (solve for each coordinate separately)
        pos_new = np.zeros_like(points)
        for i in range(3):  # For each coordinate (x, y, z)
            b             = Mt * pos_pred[:, i] +  f_ext[:, i]*dt*dt
            pos_new[:, i] = np.linalg.solve(A, b)
        
        # Update fixed points
        if fixed_points is not None:
            pos_new[fixed_points] = points[fixed_points]
        
        # Update positions
        pos_prev = pos.copy()
        pos      = pos_new.copy()

        # Update velocity
        velocity = (pos - pos_prev) / dt

        if call_back is not None:
            call_back(pos)

        print( f"iter:{itr} d={t} |v|={np.linalg.norm(velocity)}")

    return pos, velocity

# Example usage
if __name__ == "__main__":
    # Create a simple 5x5 grid
    nx, ny = 5, 5
    bonds, points, masses, ks, fixed = build_grid_2d(nx, ny)
    velocity = np.zeros_like(points)
    
    # Run simulation
    new_points, new_velocity = solve_pd(points, velocity, bonds, masses, ks, fixed_points=fixed)
