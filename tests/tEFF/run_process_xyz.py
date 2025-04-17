import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time

sys.path.append("../../")
from pyBall import eFF as eff

def extract_blocks(xyz_file):
    """Extract parameters from XYZ file comments (lines starting with #)
    Returns:
        dict: Dictionary of extracted parameters with NaNs for missing values
        (e.g. {'ang': [...], 'dist': [...], 'Etot': [...]})
    """
    all_keys = set()
    records = []
    # First pass: collect all keys and raw records
    with open(xyz_file) as f:
        for line in f:
            if line.startswith('#'):
                parts = line[1:].strip().split()
                record = {}
                for i in range(0, len(parts)-1, 2):
                    key = parts[i]
                    val = float(parts[i+1])
                    record[key] = val
                    all_keys.add(key)
                records.append(record)
    # Initialize params with NaN-filled arrays
    nrec= len(records)
    params = {key: np.full(nrec, np.nan) for key in all_keys}
    # Second pass: fill values
    for i, record in enumerate(records):
        for key, val in record.items():
            params[key][i] = val
    return params, nrec

def plot_energy_landscape( Xs, Ys, Es, Espan=None):
    """Plot energy landscape from XYZ file using imshow (simple and robust)."""
    #params, nrec = extract_blocks(xyz_file)
    xs, idx = np.unique ( Xs, return_inverse=True )
    ys, idy  = np.unique( Ys, return_inverse=True )
    energy_grid = np.full((len(xs), len(ys)), np.nan)
    # Fill the grid
    for i in range(len(Xs)):
        energy_grid[ idx[i], idy[i]] = Es[i]
    # Plot
    plt.figure(figsize=(10,8))
    if Espan is not None: 
        vmin = np.nanmin(energy_grid)
        vmax = vmin + Espan
    else:
        vmin = None
        vmax = None
    plt.imshow(energy_grid, extent=[ys.min(), ys.max(), xs.max(), xs.min()], aspect='auto', cmap='inferno', vmin=vmin, vmax=vmax)
    plt.colorbar(label='Total Energy ')
    plt.xlabel('Distance (Å)')
    plt.ylabel('Angle (rad)')
    #plt.title('Potential Energy Surface')
    
#plot_energy_landscape("export/scan_data/angdistscan_H2O.xyz")

params, nrec = extract_blocks("export/scan_data/angdistscan_H2O.xyz")
plot_energy_landscape( params['ang'], params['dist'], params['Etot'], Espan=5.0 )
plt.savefig("map2D_referece.png")

outEs = np.zeros((nrec,5))

with open("processXYZ.xyz", "w") as f: f.write("")
#eff.processXYZ( "export/scan_data/distscan_H2O.xyz", bOutXYZ=True, outEs );


eff.initOpt( dt=0.005, damping=0.005, f_limit=1000.0)
#eff.processXYZ( "export/scan_data/angdistscan_H2O.xyz", bOutXYZ=True, outEs=outEs, bCoreElectrons=False );
eff.processXYZ( "export/scan_data/angdistscan_H2O.xyz", bOutXYZ=True, outEs=outEs, bCoreElectrons=True, nstepMax=1000, dt=0.001, Fconv=1e-3, ialg=2 );
#eff.processXYZ( "export/scan_data/angdistscan_H2O.xyz", bOutXYZ=True, outEs=outEs, bCoreElectrons=True, nstepMax=0 );
#print(outEs)

plot_energy_landscape( params['ang'], params['dist'], outEs[:,0] )
plt.savefig("map2d_eFF.png")

plt.show()