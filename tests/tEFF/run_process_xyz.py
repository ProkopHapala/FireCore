from ctypes import wstring_at
import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import time

sys.path.append("../../")
from pyBall import eFF as eff

# can be used to create _ee file
elementPath = "export/scan_data/angdistscan_H2O.xyz"
elementPath_e = "export/scan_data/angdistscan_H2O_ee.xyz"
# elementPath = "export/scan_data/angdistscan_CH4.xyz"
# elementPath_e = "export/scan_data/angdistscan_CH4_ee.xyz"

def count_mask_lines( fgo_file, mask='#iconf' ):
    line = None
    nline = 0
    with open(fgo_file) as f:
        for l in f:
            if mask in l:
                line = l
                nline += 1
    return nline, line
                

def extract_nae( xyz_file ):
    with open(xyz_file) as f:
        ws = f.readline().strip().split()
        na= int(ws[0])
        ne= int(ws[1])
        return na, ne

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


if __name__ == "__main__":
    print("#=========== RUN /home/prokophapala/git/FireCore/tests/tEFF/run_process_xyz.py")
    theta0 = np.array([1.125, 0.9, -0.2])
    # theta0 = np.array([ 1.160836,  0.874741, -0.044889])
    # theta0 = np.array([ 1.153868 , 0.871132, -0.042265])
    eff.setVerbosity(0,0)
    print("verbos")
    atomParams = np.array([
    #  Q   sQ   sP   cP
    [ 0.,  1.0, 1.00, 0.0 ], # 0
    [ 1.,  0.0, 0.00, 0.0 ], # 1 H
    [ 0.,  1.0, 1.00, 1.0 ], # 2 He
    [ 1.,  0.0, 0.10, 1.0 ], # 3 Li
    [ 2.,  0.0, 0.10, 1.0 ], # 4 Be
    [ 3.,  0.0, 0.10, 1.0 ], # 5 B
    [ 4.,  0.0, 0.10, 1.0 ], # 6 C
    [ 5.,  0.0, 0.10, 1.0 ], # 7 N
    [ 6.,  0.0, 0.2, 1.0 ], # 8 O
    [ 7.,  0.0, 0.10, 1.0 ], # 9 F
    ], dtype=np.float64)
    eff.setAtomParams( atomParams )
    print("set atom par")
    params, nrec = extract_blocks(elementPath)
    plot_energy_landscape( params['ang'], params['dist'], params['Etot'], Espan=5.0 )
    plt.title("Before relaxetion")
    plt.savefig("map2D_referece.png")

    outEs = np.zeros((nrec,5))
    # apos = np.zeros((nrec,,3))
    # epos = np.zeros((nrec,4))

    with open("processXYZ.xyz", "w") as f: f.write("")
    #eff.processXYZ( "export/scan_data/distscan_H2O.xyz", bOutXYZ=True, outEs );

    eff.initOpt( dt=0.005, damping=0.005, f_limit=1000.0)

    ## ---- Previous scan

    bCoreElectrons = False
    eff.setSwitches( coreCoul=1 )
    #eff.setSwitches( coreCoul=0 )
    eff.preAllocateXYZ(elementPath, Rfac=-1.35, bCoreElectrons=bCoreElectrons )
    eff.getBuffs()
    eff.info()
    #eff.aPars[0,2]=1
    eff.esize[:]=0.7

    #eff.processXYZ( "export/scan_data/angdistscan_H2O.xyz", bOutXYZ=True, outEs=outEs, bCoreElectrons=False );
    #eff.processXYZ( "export/scan_data/angdistscan_H2O.xyz", bOutXYZ=True, outEs=outEs, bCoreElectrons=True, nstepMax=1000, dt=0.001, Fconv=1e-3, ialg=2 );
    #eff.processXYZ( "export/scan_data/angdistscan_H2O.xyz", bOutXYZ=True, outEs=outEs, bCoreElectrons=bCoreElectrons, bChangeCore=False, bChangeEsize=False, nstepMax=0 );
    eff.processXYZ( elementPath, outEs=outEs, bCoreElectrons=bCoreElectrons, bChangeCore=False, bChangeEsize=True, nstepMax=0, dt=0.005, Fconv=1e-3, ialg=2, xyz_out=elementPath_e  )
    #print(outEs)
    plot_energy_landscape( params['ang'], params['dist'], outEs[:,0] )
    plt.title("After relaxetion")
    plt.savefig("map2d_eFF.png")

    print("#=========== DONE /home/gabriel/git/FireCore/tests/tEFF/run_process_xyz.py")
    plt.show()


    '''
    ## ---- New scan over different voltages
    eff.preAllocateXYZ("export/scan_data/angdistscan_H2O.xyz", Rfac=-1.35, bCoreElectrons=True)
    eff.getBuffs()
    radii = np.linspace(0.07,0.10,0.15,0.20,0.25)  # example radii range
    energy_map = np.zeros((nrec, len(radii)))
    for j, rcore in enumerate(radii):
        eff.aPars[0,2] = rcore  # set core radius
        outEs = np.zeros((nrec,5))
        eff.processXYZ("export/scan_data/angdistscan_H2O.xyz", Rfac=-1.35, outEs=outEs, bCoreElectrons=True, bChangeCore=True, nstepMax=0)
        energy_map[:,j] = outEs[:,0]
    np.save("energy_map_radii.npy", energy_map)
    emap = energy_map[:,0]
    plot_energy_landscape(params['ang'], params['dist'], emap)
    plt.savefig("map2D_radius0.png")
    plt.show()
    '''

