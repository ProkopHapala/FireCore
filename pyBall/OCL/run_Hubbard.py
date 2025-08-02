# run it like this:
#   python -u -m pyBall.OCL.HubbardSolver | tee OUT

import os
import numpy as np
import pyopencl as cl
import matplotlib.pyplot as plt

from .OpenCLBase import OpenCLBase
from .HubbardSolver import HubbardSolver, default_params, generate_xy_scan, generate_xV_scan, make_grid_sites, make_site_ring, find_closest_pTip, save_sites_to_txt, load_sites_from_txt, make_sparse_W, make_sparse_W_pbc, screened_coulomb, kBoltzmann, COULOMB_CONST

import matplotlib.pyplot as plt

# -----------------------------------------------------------------------------
#                     High-level convenience functions 
# -----------------------------------------------------------------------------

def solve_xy_scan(solver: HubbardSolver, posE: np.ndarray, extent, nxy=(10,10), params = default_params):
    """Run an XY scan and reshape outputs to (nx,ny)."""
    zTip             = params["zTip"]
    Vbias            = params["Vbias"]
    pTips            = generate_xy_scan(extent, nxy, zTip=zTip, Vbias=Vbias)
    Emin, iMin, Itot = solver.solve(posE, pTips, params=params)
    return Emin.reshape(nxy), iMin.reshape(nxy), Itot.reshape(nxy + (2,))

def solve_xV_scan(solver: HubbardSolver, posE: np.ndarray, p1, p2, Vrange=(0.0,1.0), nxV=(10,10), params=default_params ):
    """Run an X-V scan (position along a line vs bias).  Outputs shaped (nx,nV)."""
    pTips = generate_xV_scan(p1, p2, Vrange, nxV)
    pTips[:,2] = params["zTip"]
    Emin, iMin, Itot = solver.solve(posE, pTips, params=params)
    return Emin.reshape(nxV), iMin.reshape(nxV), Itot.reshape(nxV + (2,))

def solve_pSites(params, posE: np.ndarray, solver: HubbardSolver,  ):
    """
    Debugging function that runs the solver with tip positions matching the sites.
    Returns:
        Emin : (nSites,) ground-state energies
        iMin : (nSites,) occupation patterns
        Itot : (nSites,2) tunneling currents
    """
    # Create tips at site positions with given Vbias
    pTips = posE.copy()
    pTips[:,3] = params["Vbias"]  # Replace E0 with Vbias
    pTips[:,2] = params["zTip"]
    Emin, iMin, Itot = solver.solve(posE, pTips, params=params)
    for i in range(len(Emin)):
        print(f"Site {i}: Emin={Emin[i]:.6f}, iMin={iMin[i]:b}, Itot={Itot[i]} pTips={pTips[i]}")
    return Emin, iMin, Itot

def plot_sites(posE, ax=None, c="r", ms=2.0,angles=None):
    if ax is None: 
        ax = plt.gca()
    ax.plot(posE[:,0], posE[:,1], '.', color=c, markersize=ms)
    if angles is not None: # with  quiver
        ax.quiver(posE[:,0], posE[:,1], np.cos(angles), np.sin(angles), color=c)
    ax.set_aspect('equal')
    return ax

def plot2d( data, extent=None, title=None, ax=None, cmap=None, ps=None, c='g', ms=1):
    if ax is None:
        ax = plt.gca()
    im0 = ax.imshow(data, extent=extent, cmap=cmap, origin="lower");        
    ax.figure.colorbar(im0, ax=ax) 
    plot_sites(ps,ax, c=c, ms=ms); 
    ax.set_title(title);

def plot_site_maps_imshow(Esite_map, Tsite_map, occ_map, posE, tip_pos=None, title="", sz=3, axs=None, titles=["Esite", "Tsite", "Occupancy"]):
    """
    Plots Esite, Tsite, and occupancy for a single tip configuration on a grid using imshow.

    Args:
        Esite_map (np.ndarray): 2D array of on-site energies.
        Tsite_map (np.ndarray): 2D array of hopping amplitudes.
        occ_map (np.ndarray): 2D array of unpacked occupancies (0 or 1).
        posE (np.ndarray): The (N, 4) array of site positions to determine plot extent.
        tip_pos (np.ndarray, optional): The (x, y, z, V) position of the tip to mark on the plot.
        title (str, optional): A title for the figure.
    """
    # Determine the real-space extent of the site grid for imshow
    # Add half a grid cell to each side for better visualization of corner sites
    dx = (np.max(posE[:, 0]) - np.min(posE[:, 0])) / (Esite_map.shape[0] - 1)
    dy = (np.max(posE[:, 1]) - np.min(posE[:, 1])) / (Esite_map.shape[1] - 1)
    extent = [
        np.min(posE[:, 0]) - dx/2, np.max(posE[:, 0]) + dx/2,
        np.min(posE[:, 1]) - dy/2, np.max(posE[:, 1]) + dy/2
    ]

    if axs is None:
        fig, axs = plt.subplots(1, 3, figsize=((sz+.5)*3, sz))
    else:
        fig = axs[0].figure

    if title:
        fig.suptitle(title)
    # Plot Esite
    im0 = axs[0].imshow(Esite_map, cmap='magma', origin='lower', interpolation='nearest', extent=extent)
    if titles[0] is not None: axs[0].set_title(titles[0])
    fig.colorbar(im0, ax=axs[0])

    # Plot Tsite
    im1 = axs[1].imshow(Tsite_map, cmap='magma', origin='lower', interpolation='nearest', extent=extent)
    if titles[1] is not None: axs[1].set_title(titles[1])
    fig.colorbar(im1, ax=axs[1])

    # Plot Occupancy
    im2 = axs[2].imshow(occ_map, cmap='gray', origin='lower', interpolation='nearest', extent=extent)
    if titles[2] is not None: axs[2].set_title(titles[2])
    fig.colorbar(im2, ax=axs[2])

    # Mark tip position and site positions on all subplots
    for ax in axs:
        ax.plot(posE[:, 0], posE[:, 1], '.g', markersize=1, alpha=1.0)
        if tip_pos is not None: ax.plot(tip_pos[0], tip_pos[1], 'r+', markersize=8, markeredgewidth=1.5)
        #ax.set_xlabel("X (Angstrom)")
        #ax.set_ylabel("Y (Angstrom)")
        ax.set_aspect('equal')

    plt.tight_layout(rect=[0, 0, 1, 0.96])

def plot_sites_maps_imshow(Esite, Tsite, occ_bits, tip_indices, pTips, posE, nxy_sites, nSingle, sz=2 ):
    n_configs = len(tip_indices)
    fig, axs = plt.subplots(3, n_configs, figsize=(n_configs * (sz+.5), 3*sz), squeeze=False)
    for j, i_tip in enumerate(tip_indices):
        tip_pos = pTips[i_tip]
        Esite_map = Esite[i_tip, :].reshape(nxy_sites)
        Tsite_map = Tsite[i_tip, :].reshape(nxy_sites)
        #occ_bytes = occupation[i_tip]
        #occ_map = np.unpackbits(occ_bytes, bitorder='little')[:nSingle].reshape(nxy_sites)
        occ_map = occ_bits[i_tip][:nSingle].reshape(nxy_sites)
        titles = [f"Tip {i_tip}", None, None]
        plot_site_maps_imshow(Esite_map, Tsite_map, occ_map, posE, tip_pos, title=f"Tip {i_tip}", sz=sz, axs=axs[:, j], titles=titles)
    return axs
    
def plot_site_property(pos, prop_values, titles=None, nxy_sites=None, ms=500, sz=2, thisSites=None, cmap='viridis'):
    """
    Plots a floating-point property for each site at its specific location.
    This is a general version of plot_site_occupancy for continuous properties.
    For each configuration (e.g., tip position), it creates a subplot showing the
    property values for all sites using a color map.

    Args:
        pos (np.ndarray): The (N, 3) or (N, 4) array of site positions.
        prop_values (np.ndarray): A 2D array (n_configs, n_sites) of the property to plot.
        energy (np.ndarray, optional): 1D array of energies for each configuration to display in the title.
        nxy_sites (tuple): The (nx, ny) dimensions of the subplot grid.
        sz (int): Marker size for the scatter plot.
        thisSites (np.ndarray, optional): The (n_configs, 3) or (n_configs, 4) array of tip positions.
        title_prefix (str): A prefix for the subplot titles.
        cmap (str): Colormap for the scatter plot.
    """
    n_configs, n_sites = prop_values.shape if prop_values.ndim == 2 else (1, prop_values.shape[0])
    nx, ny = nxy_sites

    if n_configs != (nx * ny):
        print(f"Warning: Number of configurations ({n_configs}) does not match grid size ({nx*ny}). Truncating.")
        n_configs = min(n_configs, nx * ny)
    fig, axs = plt.subplots(ny, nx, figsize=(nx * sz, ny * sz))
    if n_configs == 1:
        axs = np.array([axs])
    axs = axs.flatten() # Make it easier to iterate

    vmin = np.min(prop_values)
    vmax = np.max(prop_values)

    for s in range(n_configs):
        ax = axs[s]
        if thisSites is not None:  ax.plot(thisSites[s, 0], thisSites[s, 1], 'r+', markersize=8, markeredgewidth=2, zorder=10)
        sc = ax.scatter(pos[:, 0], pos[:, 1], c=prop_values[s, :], s=ms, cmap=cmap, vmin=vmin, vmax=vmax, edgecolors='none')
        if titles is not None:  ax.set_title(titles[s])
        ax.set_aspect('equal')
        #ax.set_xlabel("X [A]")
        #ax.set_ylabel("Y [A]")
    fig.colorbar(sc, ax=axs.ravel().tolist(), orientation='vertical', fraction=0.02, pad=0.04)
    plt.tight_layout()

def plot_site_occupancy(pos,  occupation, energy, nxy_sites, nSingle, ms=500, sz=2, thisSites=None,):
    """
    Plots a binary occupancy pattern for each site at its specific location.
    This is a general version of plot_site_property for boolean properties.
    For each configuration (e.g., tip position), it creates a subplot showing the
    occupancy pattern for all sites using a color map.

    Args:
        pos (np.ndarray): The (N, 3) or (N, 4) array of site positions.
        occupation (np.ndarray): A 2D array (n_configs, n_sites) of the occupancy pattern to plot.
        energy (np.ndarray, optional): 1D array of energies for each configuration to display in the title.
        nxy_sites (tuple): The (nx, ny) dimensions of the subplot grid.
        sz (int): Marker size for the scatter plot.
        thisSites (np.ndarray, optional): The (n_configs, 3) or (n_configs, 4) array of tip positions.
        title_prefix (str): A prefix for the subplot titles.
        cmap (str): Colormap for the scatter plot.
    """
    print( "plot_site_occupancy() occupation.shape ", occupation.shape)

    if thisSites is None:
        thisSites = pos.copy()

    #bits_all = np.unpackbits(occupation.reshape(nTips, solver.occ_bytes), axis=1)[:, :nSingle].astype(bool)
    bits_all = np.unpackbits(occupation, axis=1)[:, :nSingle].astype(bool)

    fig, axs = plt.subplots(nxy_sites[1], nxy_sites[0], figsize=(nxy_sites[0]*sz, nxy_sites[1]*sz))
    #plt.subplots_adjust(wspace=0.3, hspace=0.3)
    for s in range(nSingle):
        ix = s // nxy_sites[0]
        jy = s % nxy_sites[0]
        ax = axs[ix, jy]
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_color('black')
            spine.set_linewidth(2)
        occ_config = bits_all[s].astype(bool)
        ax.plot(thisSites[s, 0], thisSites[s, 1], 'r+', markersize=8, markeredgewidth=2)
        ax.scatter(pos[:, 0][~occ_config], pos[:, 1][~occ_config], s=ms, facecolors='none', edgecolors='black', linewidths=2)
        ax.scatter(pos[:, 0][occ_config],  pos[:, 1][occ_config],  s=ms, color='black')
        ax.set_title(f"E={energy[s]:.3f}")
        ax.set_aspect('equal')
        #ax.axis('off')
    #plt.tight_layout()
    #plt.show()

def plot_occupancy_slice(occupation, nxy_scan, occ_bytes, slice_idx, axis='x'):
    """
    Plots a 1D slice of the occupancy array as a binary image.

    Args:
        occupation (np.ndarray): The 1D array of site occupancies (nTips * occ_bytes).
        nxy_scan (tuple): The (nx, ny) dimensions of the tip scan grid.
        occ_bytes (int): The number of bytes per tip position for occupancy.
        slice_idx (int): The index for the slice (e.g., ix if axis='x', iy if axis='y').
        axis (str): The axis along which to take the slice ('x' or 'y').
    """
    nx, ny = nxy_scan
    nTips  = nx * ny
    nSites = occ_bytes * 8
    
    occupancy_2d_bytes = occupation.reshape(ny, nx, occ_bytes)
    
    if axis == 'x':
        slice_data_bytes = occupancy_2d_bytes[:, slice_idx, :]
        title_suffix = f" (x_idx={slice_idx})"
    elif axis == 'y':
        slice_data_bytes = occupancy_2d_bytes[slice_idx, :, :]
        title_suffix = f" (y_idx={slice_idx})"
    else:
        print("Error: axis must be 'x' or 'y'")
        return

    binary_image_data = np.zeros((slice_data_bytes.shape[0], nSites), dtype=np.uint8)
    for i, byte_row in enumerate(slice_data_bytes):
        bits = np.unpackbits(byte_row, bitorder='little')
        binary_image_data[i, :] = bits

    plt.figure(figsize=(10, 5))
    plt.imshow(binary_image_data, cmap='gray', aspect='auto', interpolation='nearest')
    plt.title(f"Site Occupancy Slice{title_suffix}")
    plt.xlabel("Site Index (0 to nSites-1)")
    plt.ylabel(f"Tip Position Index along {'Y' if axis=='x' else 'X'} axis")
    plt.colorbar(label="Occupancy (0=empty, 1=occupied)")
    #plt.show()

def plot_occupancy_line(occupation, nSingle):
    nbyte = nSingle//8
    bits = np.unpackbits(occupation, axis=1, bitorder='little')
    print("plot_occupancy_line() bits.shape = ", bits.shape, nSingle, nbyte)
    bits = bits[:,:nSingle]
    plt.figure(figsize=(10, 5))
    plt.imshow(bits, cmap='gray', aspect='auto', interpolation='nearest')
    plt.title(f"Site Occupancy Slice")
    plt.xlabel("Site Index (0 to nSites-1)")
    plt.ylabel(f"Tip Position Index ")
    plt.colorbar(label="Occupancy (0=empty, 1=occupied)")
    #plt.show()

# New helper functions for demo_local_update
def solve_hopping_scan(solver: HubbardSolver, posE: np.ndarray, rots: np.ndarray, orbs: np.ndarray, extent, nxy=(10,10), params=default_params):
    """Run oriented hopping scan over XY grid and sum amplitudes over sites."""
    pTips = generate_xy_scan(extent, nxy, zTip=params["zTip"], Vbias=params["Vbias"])
    nSingle = posE.shape[0]
    Tout = solver.eval_oriented_hopping(posE, rots, orbs, pTips, bRealloc=True)
    return Tout.reshape(nxy+(nSingle,))

# -----------------------------------------------------------------------------
#                            Tests  and Demos 
# -----------------------------------------------------------------------------

def test_brute_force_solver():
    # trimer
    #posE = make_site_ring(3, 5.0, E0=params["E0"] )
    #save_sites_to_txt("trimer.txt", posE)

    # hexamer
    params={ 
        "tipDecay"  : 0.3,
        "tipRadius" : 3.0,
        "zTip"      : 5.0,
        "zMirror"   : -0.5,
        "Wcouple"   : 1.0,
        "E0"        : -0.10,
        #"Vbias"     : +0.3,
        #"Vbias"     : +0.2,
        "Vbias"     : +0.1,
        #"Vbias"     : -0.45,    
        #"Vbias"     : -0.55,
    }
    #outer,a1 = make_site_ring(3, 8.0, E0=params["E0"], ang0=np.pi*0.5 )
    #inner,a2 = make_site_ring(3, 5.2, E0=params["E0"], ang0=np.pi*1.5 )

    outer,a1 = make_site_ring(3, 0.86602540378*(2.0/3.0)*20.0, E0=params["E0"], ang0=np.pi*0.5 )
    inner,a2 = make_site_ring(3, 0.86602540378*(2.0/3.0)*10.0, E0=params["E0"], ang0=np.pi*1.5 )

    print( "|p1-p0| ",  np.sqrt( ((inner[0,:]-inner[1,:])**2).sum() ) )

    #outer = make_site_ring(3, 5.0, E0=params["E0"], ang0=0.0 )
    #inner = make_site_ring(3, 3.0, E0=params["E0"], ang0=np.pi )
    posE  = np.vstack([outer, inner])
    angles = np.hstack([a1,a2])
    print("posE.shape", posE.shape)
    #save_sites_to_txt("hexamer.txt", posE)
    #np.savetxt("hexamer.txt", posE)
    xya      = posE.copy()
    xya[:,2] = angles
    save_sites_to_txt("hexamer.txt", xya)
    #np.savetxt("hexamer.txt", xya)

    plot_sites(posE, angles=angles)
    #plt.show(); exit(0)
    
    # # ---- grid
    # params={ 
    #     "tipDecay"  : 0.3,
    #     "tipRadius" : 3.0,
    #     "zTip"      : 1.0,
    #     "zMirror"   : -0.5,
    #     "Wcouple"   : 1.0,
    #     "E0"        : -0.15,
    #     "Vbias"     : 0.2,
    # }
    # posE = make_grid_sites( nxy=(4,4), avec=(7.0,0.0), bvec=(0.0,5.0), E0=params["E0"] )
    # print("posE:\n", posE)
    # save_sites_to_txt("grid.txt", posE)

    solver = HubbardSolver()
    #solve_pSites(params, posE, solver)

    Rij, Wij = solver.eval_coupling_matrix(posE, params["zMirror"])
    #print("Rij:\n", Rij)
    #print("Wij:\n", Wij)

    extent=[-10.0,10.0,-10.0,10.0]
    # #posE = load_sites_from_txt("hexamer.txt")
    # Oriented hopping demo
    # Create random site rotations and orbitals for demonstration
    nSingle = posE.shape[0]
    # Simple rots: all zeros (no rotation)
    phis = np.zeros(nSingle)
    rots = np.vstack([np.cos(phis), np.sin(phis)]).T.astype(np.float32)
    # Simple orbitals: s-wave only with decay=1.0
    orbs = np.zeros((nSingle,4), dtype=np.float32)
    orbs[:,0] = 1.0  # C_s
    orbs[:,3] = 0.2  # decay
    Tsites = solve_hopping_scan(solver, posE, rots, orbs, extent=extent, nxy=(10,10), params=params)
    print("Tsites: shape: ", Tsites.shape, " min: ", Tsites.min(), " max: ", Tsites.max())
    # Plot hopping amplitude scan
    
    fig2, axs = plt.subplots( 1,nSingle, figsize=(nSingle*4,4) )
    for i in range(nSingle):
        plot2d(Tsites[i], extent=extent, title="Total Hopping Amplitude", ax=axs[i], cmap="inferno", ps=posE)
    # Also run XY scan for comparison plots
    
    plt.show(); exit()
    
    Emin, iMin, Itot = solve_xy_scan(solver, posE, extent=extent, nxy=(200,200), params=params)
    fig,axs = plt.subplots(2,2)
    # im0 = axs[0,0].imshow(Emin,        extent=extent);        fig.colorbar(im0, ax=axs[0,0]); plot_sites(posE,axs[0,0]); axs[0,0].set_title("Emin");
    # im1 = axs[0,1].imshow(iMin,        extent=extent);        fig.colorbar(im1, ax=axs[0,1]); plot_sites(posE,axs[0,1]); axs[0,1].set_title("iMin");
    # im2 = axs[1,0].imshow(Itot[:,:,0], extent=extent); fig.colorbar(im2, ax=axs[1,0]); plot_sites(posE,axs[1,0]); axs[1,0].set_title("I_occ");
    # im3 = axs[1,1].imshow(Itot[:,:,1], extent=extent); fig.colorbar(im3, ax=axs[1,1]); plot_sites(posE,axs[1,1]); axs[1,1].set_title("I_unocc");
    plot2d(Emin,        extent=extent, title="Emin",    ax=axs[0,0], cmap="viridis", ps=posE);
    plot2d(iMin,        extent=extent, title="iMin",    ax=axs[0,1], cmap="viridis", ps=posE);
    plot2d(Itot[:,:,0], extent=extent, title="I_occ",   ax=axs[1,0], cmap="viridis", ps=posE);
    plot2d(Itot[:,:,1], extent=extent, title="I_unocc", ax=axs[1,1], cmap="viridis", ps=posE);
    plt.tight_layout()
    plt.show()

def test_local_update_solver():
    # --- Example Usage for the new Local Update Solver ---
    
    # 1. Setup a larger system
    posE  = make_grid_sites(nxy=(8, 8), avec=(5.0, 0.0), bvec=(0.0, 5.0), E0=-0.15)
    nSite = posE.shape[0]
    print(f"Created a large system with {nSite} sites.")

    solver = HubbardSolver()

    # 2. Pre-calculate the dense interaction matrix
    _, Wij_dense = solver.eval_coupling_matrix(posE, zMirror=-0.5)

    # 3. Convert the dense matrix to the required sparse format
    nMaxNeigh = 16
    W_sparse = make_sparse_W(Wij_dense, nMaxNeigh=nMaxNeigh)
    print(f"Converted dense Wij to sparse format with nMaxNeigh={nMaxNeigh}.")

    # 4. Pre-calculate Esite and Thop for a few tip positions
    pTips = generate_xy_scan(extent=[-5, 5, -5, 5], nxy=(3, 1), zTip=3.0, Vbias=0.2)
    nTips = pTips.shape[0]

    # This part would typically involve a more complex physical model,
    # but for demonstration, we'll construct it manually.
    Esite = np.zeros((nTips, nSite), dtype=np.float32)
    Tsite = np.zeros((nTips, nSite), dtype=np.float32)
    for i in range(nTips):
        tip_pos   = pTips[i, :3]
        tip_vbias = pTips[i, 3]
        for j in range(nSite):
            dist = np.linalg.norm(tip_pos - posE[j, :3])
            # Simplified model for Esite and Thop
            Esite[i, j] = posE[j, 3] + tip_vbias * (3.0 / dist) # On-site energy
            Tsite[i, j] = np.exp(-0.3 * dist)                 # Hopping
    
    # 5. Run the local update solver
    print("Running the parallel local update solver...")
    # solver_params = {
    #     "kT": 0.001,
    #     "nIter": 20000,
    #     "solverMode": 2, # Annealing
    #     "seed": 42
    # }
    E_out, state_out, Itot_out = solver.solve_local_updates(Esite, Tsite, W_sparse )
    
    # 6. Print results
    for i in range(nTips):
        print(f"--- Tip {i} ---")
        print(f"  Final Energy: {E_out[i]:.6f}")
        print(f"  Final State Mask: {state_out[i]:0{nSite}b}") # Print as binary string
        print(f"  Final Currents: {Itot_out[i]}")

def make_site_lattice( lvecs, nxy=(6,6), p0=(0.0,0.0,0.0) ):
    nx, ny = nxy
    nSite = nx * ny
    pos = np.zeros((nSite, 3), dtype=np.float32)
    idx = 0
    for iy in range(ny):
        for ix in range(nx):
            pos[idx, 0] = ix * lvecs[0,0] + p0[0]
            pos[idx, 1] = iy * lvecs[1,1] + p0[1]
            pos[idx, 2] = p0[2]
            idx += 1
    return pos

def test_site_coupling( nxy=(6,6), cutoff=8.0, a=5.0, b=5.0, nMaxNeigh=8 ):
    # Lattice vectors define the size of the periodic cell
    lvs = np.array([ [a, 0.0, 0.0], [0.0, b, 0.0] ], dtype=np.float32)
    pos = make_site_lattice(lvs, nxy=nxy, p0=(0.0,0.0,0.0))
    lvs_pbc = lvs.copy(); lvs_pbc[1,:] *= nxy[1]; lvs_pbc[0,:] *= nxy[0] 
    W_val, W_idx, nNeigh = make_sparse_W_pbc( pos, lvs_pbc, cutoff, W_func=screened_coulomb, nMaxNeigh=nMaxNeigh )
    print("nNeigh:", nNeigh)
    for i in range(len(nNeigh)):
        print(f"Site {i}: nNeigh: {nNeigh[i]} ",end="")
        for j in range(nNeigh[i]):
            print(f"  {W_idx[i,j]}: {W_val[i,j]:.2f}", end=" ")
        print()

def demo_precalc_scan(solver: HubbardSolver=None, nxy_sites=(6, 6), nxy_scan=(100, 100), Vbias=0.1):
    """
    Demonstrates the usage of the precalc_esite_thop kernel by running a 2D scan.

    This function visualizes the fundamental interaction landscape between the tip and
    the individual sites before any many-body effects are considered. It plots maps of:
    - The minimum on-site energy induced by the tip.
    - The maximum on-site energy induced by the tip.
    - The minimum hopping amplitude (i.e., the site the tip is closest to).

    Args:
        solver (HubbardSolver): An instance of your solver class.
        nxy_sites (tuple): The (nx, ny) dimensions of the grid of sample sites.
        nxy_scan (tuple): The (nx, ny) resolution of the tip scan.
        Vbias (float): The bias voltage to apply to the tip.
    """
    if solver is None:
        solver = HubbardSolver()

    print("--- Running pre-calculation scan demo ---")

    # 1. Define the system of sites
    posE    = make_grid_sites(nxy=nxy_sites, avec=(5.0, 0.0), bvec=(0.0, 5.0), E0=-0.15)
    nSingle = posE.shape[0]
    print(f"Created a {nxy_sites[0]}x{nxy_sites[1]} grid of {nSingle} sites.")

    # 2. Define the 2D grid of tip positions for the scan
    extent = [-15.0, 15.0, -15.0, 15.0]
    pTips  = generate_xy_scan(extent, nxy=nxy_scan, zTip=3.0, Vbias=Vbias)
    nTips  = pTips.shape[0]
    print(f"Generated a {nxy_scan[0]}x{nxy_scan[1]} grid of {nTips} tip positions.")

    # 3. Define the physical model parameters (e.g., simple monopole)
    multipoleCoefs = np.zeros((nSingle, 1), dtype=np.float32)
    multipoleCoefs[:, 0] = 1.0  # All sites are simple monopoles
    
    params = {
        "Rtip": 3.0,
        "zMirror": -0.5,
        "Thop_decay": 0.5, # A higher decay for more localized hopping
        "Thop_amp": 1.0
    }

    # 4. Call the pre-calculation kernel
    print("Calling precalc_esite_thop kernel...")
    Esite,Tsite = solver.precalc_esite_thop(posE, pTips, multipoleCoefs=multipoleCoefs, params=params)
    # Output shape is (nTips, nSingle, 2)
    print(f"Esite min,max : {np.min(Esite), np.max(Esite)}")
    print(f"Tsite min,max : {np.min(Tsite), np.max(Tsite)}")
    
    # 5. Perform reduction in Python/NumPy to get the desired maps
    # For each tip position (axis 0), find the min/max over all sites (axis 1)
    min_E_map    = np.min(Esite, axis=1) # Min Esite
    max_E_map    = np.max(Esite, axis=1) # Max Esite
    max_Thop_map = np.max(Tsite, axis=1) # Max Thop    # Hopping is always positive, so min is fine. Use max for "strongest hopping".
    
    # 6. Reshape the 1D arrays back into 2D maps for plotting
    nx, ny = nxy_scan
    min_E_map    = min_E_map.reshape((nx, ny))
    max_E_map    = max_E_map.reshape((nx, ny))
    max_Thop_map = max_Thop_map.reshape((nx, ny))

    # 7. Plot the results using the provided helper functions
    fig, axs = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle(f"Tip-Site Interaction Maps (Vbias = {Vbias:.2f} V)")

    plot2d(min_E_map,    extent=extent, title="Minimum Site Energy (E_min)", ax=axs[0], cmap="viridis", ps=posE)
    plot2d(max_E_map,    extent=extent, title="Maximum Site Energy (E_max)", ax=axs[1], cmap="magma",   ps=posE)
    plot2d(max_Thop_map, extent=extent, title="Maximum Hopping (T_max)",     ax=axs[2], cmap="inferno", ps=posE)
           
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.show()

def demo_local_update(
        solver: HubbardSolver=None, 
        nxy_sites =(8,8), 
        avec      =(15.04, 0.0), 
        bvec      =(0.0, 11.57), 
        #avec      =(5.0, 0.0), 
        #bvec      =(0.0, 5.0), 
        nxy_scan  =(200, 200), 
        extent    =[-20.0, 20.0, -20.0, 20.0], 
        E0=    -0.15, 
        zTip=3.0, 
        Vbias=0.5, 
        cutoff=20.0, 
        W_amplitude=0.5, 
        T=1.0, 
        nIter=100000, 
        solverMode=2, 
        Efermi=0.0,
        #fig_dir=None,
        fig_dir="./figs",
    ):
    """
    Demonstrates the usage of the solve_local_updates kernel by running a 2D scan with Monte Carlo optimization.
    
    This function extends demo_precalc_scan by adding the many-body effects through Monte Carlo optimization.
    It calculates site energies and hoppings, then runs the local update solver for each tip position.
    The function plots maps of:
    - The minimum energy after Monte Carlo optimization
    - The maximum hopping amplitude after optimization
    - The total charge (number of occupied sites) for each tip position
    
    Args:
        solver (HubbardSolver): An instance of your solver class.
        nxy_sites (tuple): The (nx, ny) dimensions of the grid of sample sites.
        nxy_scan (tuple): The (nx, ny) resolution of the tip scan.
        Vbias (float): The bias voltage to apply to the tip.
        cutoff (float): Cutoff distance for site-site interactions.
        W_amplitude (float): Strength of the Coulomb interaction between sites.
        T (float): Temperature for the Monte Carlo simulation in Kelvin.
        nIter (int): Number of Monte Carlo iterations.
    """
    if solver is None:
        solver = HubbardSolver()

    print("--- Running local update Monte Carlo solver demo ---")

    # 1. Define the system of sites
    posE    = make_grid_sites(nxy=nxy_sites, avec=avec, bvec=bvec, E0=E0)
    nSingle = posE.shape[0]
    print(f"Created a {nxy_sites[0]}x{nxy_sites[1]} grid of {nSingle} sites.")

    # 2. Define the 2D grid of tip positions for the scan
    #extent = [-15.0, 15.0, -15.0, 15.0]
    pTips = generate_xy_scan(extent=extent, nxy=nxy_scan, zTip=zTip, Vbias=Vbias)
    nTips = pTips.shape[0]
    print(f"Generated a {nxy_scan[0]}x{nxy_scan[1]} grid of {nTips} tip positions.")

    # 3. Define the physical model parameters
    # 3a. Multipole coefficients (simple monopole)
    multipoleCoefs = np.zeros((nSingle, 1), dtype=np.float32)
    multipoleCoefs[:, 0] = 1.0 # All sites are simple monopoles
    
    # 3b. Pre-calculation kernel parameters
    params_precalc = {
        "Rtip": 3.0,
        "zMirror": -0.5,
        "Thop_decay": 0.1, # A higher decay for more localized hopping
        "Thop_amp":   1.0
    }
    
    # 3c. Generate site-site interactions (W matrix)
    # Periodic boundary conditions for the lattice
    lvs = np.array([[avec[0], avec[1], 0.0], [bvec[0], bvec[1], 0.0]], dtype=np.float32)
    lvs_pbc = lvs.copy()
    lvs_pbc[0,:] *= nxy_sites[0]
    lvs_pbc[1,:] *= nxy_sites[1]
    
    # Function to calculate Coulomb interaction with screening
    def screened_coulomb(r):
        return W_amplitude * np.exp(-r/cutoff) / (r + 0.1)
    
    # Create sparse W matrix with periodic boundary conditions
    nMaxNeigh = 8  # Maximum number of neighbors per site
    # Extract only the position components (x,y,z) from posE for the PBC calculation
    pos_xyz = posE[:, :3]  # Only take x, y, z components, not the energy E0
    W_val, W_idx, nNeigh = make_sparse_W_pbc(pos_xyz, lvs_pbc, cutoff, W_func=screened_coulomb, nMaxNeigh=nMaxNeigh)
    print(f"demo_local_update() Created sparse W matrix with cutoff {cutoff} and max {nMaxNeigh} neighbors per site.")


    plt.imshow(W_val)
    
    # 3d. Local update solver parameters
    params_solver = {
        "kT":    kBoltzmann * T,  # Temperature in energy units
        "nIter": nIter,      # Number of Monte Carlo iterations
        "solverMode": solverMode,     # Use Metropolis algorithm (0=deterministic, 2=simulated annealing)
        "seed":  np.random.randint(0, 2**31)  # Random seed
    }

    # Create figure directory if it doesn't exist
    os.makedirs(fig_dir, exist_ok=True)
    print(f"Saving figures to: {os.path.abspath(fig_dir)}")

    # 4. Call the pre-calculation kernel
    print("demo_local_update() Calling precalc_esite_thop kernel...")
    Esite, Tsite = solver.precalc_esite_thop(posE, pTips, multipoleCoefs=multipoleCoefs, params=params_precalc)

    #W_val *=0.05
    #W_val *=1.0
    #W_val *=0.5
    #W_val *=0.1
    #W_val *=0.05
    #W_val *=0.03
    #W_val *=0.02
    #W_val *=0.01
    #W_val *=0.0
    Esite += Efermi

    #print("demo_local_update() W_val: \n", W_val)
    #print("demo_local_update() W_idx: \n", W_idx)

    print(f"demo_local_update() Kernel finished. Sites: {nSingle}, Tips: {nTips}")
    print(f"demo_local_update() Esite shape: {Esite.shape}, min/max: {np.min(Esite):.3f}/{np.max(Esite):.3f}")
    print(f"demo_local_update() Tsite shape: {Tsite.shape}, min/max: {np.min(Tsite):.3f}/{np.max(Tsite):.3f}")

    # 5. Run the local update Monte Carlo solver
    print(f"demo_local_update() Running local update solver with T={T}K, nIter={nIter}...")
    
    # Allocate buffers for the solver
    solver.realloc_local_update_buffers(nSingle, nTips, nMaxNeigh)

    #bNoCoupling = True
    bNoCoupling = False

    #energy, current, occupation = solver.solve_local_updates( W_sparse=(W_val, W_idx, nNeigh), Esite=Esite, Tsite=Tsite, nTips=nTips, nSite=nSingle, nMaxNeigh=nMaxNeigh, params=params_solver, initMode=0, bNoCoupling=bNoCoupling )
    # energy, current, occupation = solver.solve_mc(W_sparse=(W_val, W_idx, nNeigh), Esite=Esite, Tsite=Tsite, nTips=nTips, nSite=nSingle, params=params_solver, bRealloc=True, initMode=3,  nxy_scan=nxy_scan,                  
    #     nLocalIter=200, prob_params=( 0.1, 0.0, 0.5, 0.0) )

    nGlobalSteps = 100

    # solve_mc_2phase(self, W_sparse, Esite, Tsite, nTips, nSite, nx, nGlobalSteps=100, nLocalIter=100, prob_params=(0.1, 0.6, 0.3), # (p_best, p_neighbor, p_random), bAlloc=True, bFinalize=True):
    # phase 1 - exploration
    solver                              .solve_mc_2phase(W_sparse=(W_val, W_idx, nNeigh), Esite=Esite, Tsite=Tsite, nTips=nTips, nSite=nSingle, nx=nxy_scan[0], nGlobalSteps=nGlobalSteps, nLocalIter=50,  prob_params=( 0.1, 0.0, 0.9, 0.0), bAlloc=True, bFinalize=False )
    # phase 2 - exploitation
    energy, current, occupation = solver.solve_mc_2phase(W_sparse=(W_val, W_idx, nNeigh), Esite=Esite, Tsite=Tsite, nTips=nTips, nSite=nSingle, nx=nxy_scan[0], nGlobalSteps=nGlobalSteps, nLocalIter=150, prob_params=( 0.1, 0.5, 0.5, 0.0), bAlloc=False, bFinalize=True )

    print( "demo_local_update() energy.shape: ", energy.shape )
    #print( "demo_local_update() current.shape: ", current.shape )
    #print( "demo_local_update() occupation.shape: ", occupation.shape, solver.occ_bytes )

    # --- Plotting ---

    # Reshape occupation array for easier indexing
    occupation_reshaped = occupation.reshape(nTips, solver.occ_bytes)

    #site_tip_indices = find_closest_pTip(posE, pTips, nxy_scan)

    #nSpec = 5
    #spec_tips = np.zeros((nSpec,4))
    #spec_tips[:,0] = np.linspace(-10.0, 10.0, nSpec, endpoint=True)

    bRow = True
    if bRow:
        middle_row = nxy_sites[1] // 2  # Integer division to get middle row index
        spec_sites = posE[middle_row*nxy_sites[0]:(middle_row+1)*nxy_sites[0], :]
    else:
        middle_col = nxy_sites[0] // 2  # Integer division to get middle column index
        spec_sites = posE[middle_col::nxy_sites[0], :]  # Step by row length to get column
    spec_inds = find_closest_pTip( spec_sites, pTips)

    # # 1. Plot occupancy, energy and tunneling patterns for tips closest to each site (original scatter plot)
    # print("Plotting occupancy patterns for tip-on-site configurations...")
    # plot_site_occupancy(posE, occupation_reshaped[site_tip_indices], energy[site_tip_indices], nxy_sites, nSingle)
    # plt.suptitle("Occupancy Patterns (Tip Closest to Each Site)")
    # print("Plotting on-site energy property for tip-on-site configurations...")
    # plot_site_property(posE, Esite[site_tip_indices, :], titles=[f"E={e:.3f}" for e in energy[site_tip_indices]], nxy_sites=nxy_sites, cmap='plasma', thisSites=pTips[site_tip_indices])
    # plt.suptitle("On-Site Energy (Esite)")
    # plot_site_property(posE, Tsite[site_tip_indices, :], titles=[f"T={t:.3f}" for t in energy[site_tip_indices]], nxy_sites=nxy_sites, cmap='plasma', thisSites=pTips[site_tip_indices])
    # plt.suptitle("Tunneling (Tsite)")

    # 3. Plot detailed maps for a single, central tip position using imshow
    # i_tip_center = nTips // 2
    # tip_pos_center = pTips[i_tip_center]
    # Esite_map_center = Esite[i_tip_center, :].reshape(nxy_sites)
    # Tsite_map_center = Tsite[i_tip_center, :].reshape(nxy_sites)
    # occ_bytes_center = occupation_reshaped[i_tip_center]
    # occ_map_center = np.unpackbits(occ_bytes_center)[:nSingle].reshape(nxy_sites)
    # plot_site_maps_imshow(  Esite_map_center, Tsite_map_center, occ_map_center, posE=posE, tip_pos=tip_pos_center, title=f"Site Maps for Tip at ({tip_pos_center[0]:.2f}, {tip_pos_center[1]:.2f})" )


    bits = np.unpackbits(occupation_reshaped, axis=1, bitorder='little')   # shape (nTips, 8*occ_bytes)

    #isite0 = 0
    #itip0  = site_tip_indices[isite0]
    # print_site_maps(Esite_map, Tsite_map, occ_map, total_energy, total_current, title=""):
    itip0 = spec_inds[0]
    #print_site_maps( Esite[itip0], Tsite[itip0], bits[itip0], total_energy=energy[itip0], total_current=current[itip0], title=f"Site Maps for iTip# {itip0} near site# {isite0} at pos({pTips[itip0][0]:.2f}, {pTips[itip0][1]:.2f})" )
    #print_site_maps( Esite[itip0], Tsite[itip0], bits[itip0], total_energy=energy[itip0], total_current=current[itip0], title=f"Site Maps for iTip# {itip0} at pos({ pTips[itip0][0]:.2f}, {pTips[itip0][1]:.2f})" )


    # 4. Plot detailed maps for a subset of tip-on-site configurations using imshow
    print("Plotting detailed maps for a subset of tip-on-site configurations...")

    axs = plot_sites_maps_imshow( Esite, Tsite, bits, tip_indices=spec_inds, pTips=pTips, posE=posE, nxy_sites=nxy_sites, nSingle=nSingle )
    axs[0,0].plot( spec_sites[:,0], spec_sites[:,1], '+g', ms=10 )
    plt.suptitle("Property maps for 5 closest tip-on-site configurations")
    plt.savefig(os.path.join(fig_dir, f"site_maps_imshow_V{Vbias:.2f}.png"))
    #plt.close()

    # Calculate total charge for each tip position by summing only the relevant bits
    total_charge = bits[:, :nSingle].sum(axis=1)

    # Reshape results into 2D maps
    nx, ny = nxy_scan
    energy_map       = energy       .reshape((nx, ny))
    current_map_occ  = current[:, 0].reshape((nx, ny))   # Occupied sites current
    current_map_unoc = current[:, 1].reshape((nx, ny))  # Unoccupied sites current
    charge_map       = total_charge .reshape((nx, ny))
    
    print(f"Energy range: {np.min(energy):.3f} to {np.max(energy):.3f}")
    print(f"Charge range: {np.min(total_charge)} to {np.max(total_charge)} sites")
    print(f"Current (occupied) range: {np.min(current[:, 0]):.3f} to {np.max(current[:, 0]):.3f}")
    
    # 5. Plot final summary maps (energy, charge, etc.)
    fig, axs = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle(f"Monte Carlo Optimization Results (T={T}K, nIter={nIter}, W={W_amplitude})", fontsize=16)
    
    cmap='plasma'
    plot2d(energy_map.T,       extent=extent, title="Energy (optimized)",              ax=axs[0, 0], cmap=cmap, ps=posE )
    plot2d(charge_map.T,       extent=extent, title="Total Charge (# occupied sites)", ax=axs[0, 1], cmap=cmap, ps=posE )
    plot2d(current_map_occ.T,  extent=extent, title="Current (occupied sites)",        ax=axs[1, 0], cmap=cmap, ps=posE )
    plot2d(current_map_unoc.T, extent=extent, title="Current (unoccupied sites)",      ax=axs[1, 1], cmap=cmap, ps=posE )
    
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(os.path.join(fig_dir, f"monte_carlo_results_V{Vbias:.2f}.png"))
    #plt.close()

    plt.show()


if __name__ == "__main__":

    # run it like this:
    #   python -u -m pyBall.OCL.run_Hubbard | tee OUT-Hubbard

    #test_brute_force_solver()
    #test_site_coupling()
    #demo_precalc_scan()
    demo_local_update()

    
