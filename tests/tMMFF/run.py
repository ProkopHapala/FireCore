import sys
import numpy as np
import os
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import MMFF as mmff


########## 1D Scan #####################################################
def scanPlot(nscan=1000, span=(0.0,8.0), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0.0), label="E", saveFig=None, saveData=None):
    
    ts = np.linspace(span[0], span[1], nscan, endpoint=False)
  
    poss = np.zeros((nscan, 3))
    poss[:, 0] = p0[0] + ts * dir[0]
    poss[:, 1] = p0[1] + ts * dir[1]
    poss[:, 2] = p0[2] + ts * dir[2]

   
    Es, Fs, Ps = mmff.scan(poss, bF=True, bP=True)
    

    if saveData is not None:
        np.savetxt(saveData, np.column_stack((ts, Es)), header="ts\tEnergy", comments="# ")

    plt.title(label)
    plt.plot(ts, Es, '-', lw=0.5, label=label)
    plt.xlabel(f"Scaned along ({dir[0]}_{dir[1]}_{dir[2]}) direction ")
    plt.ylabel(f"Scaned Energy (eV)")
    
    # Optionally, save the figure to a file.
    if saveFig is not None:
        plt.savefig(saveFig)
    plt.show()



########## 2D Scan #####################################################
def scanPlot2D(nscan1=1000, nscan2=1000, span1=(0.0,4.0), span2=(0.0,4.0),
               p0=(0.0,0.0,0.0), dir1=(1.0,0.0,0.0), dir2=(0.0,1.0,0.0),
               label="E_2D", saveFig=None, saveData=None):
    # Create linspace arrays for both scan directions
    t1 = np.linspace(span1[0], span1[1], nscan1, endpoint=False)
    t2 = np.linspace(span2[0], span2[1], nscan2, endpoint=False)
    
    # Generate a meshgrid for these parameters.
    # 'indexing="ij"' ensures that T1 corresponds to the first dimension and T2 to the second.
    T1, T2 = np.meshgrid(t1, t2, indexing="ij")
    
    # Prepare an array for positions with shape (nscan1*nscan2, 3)
    poss = np.zeros((nscan1*nscan2, 3))
    
    # Each scanned position is the starting point plus contributions along two directions:
    # p = p0 + t1*dir1 + t2*dir2 for each point in the grid.
    poss[:, 0] = p0[0] + T1.ravel()*dir1[0] + T2.ravel()*dir2[0]
    poss[:, 1] = p0[1] + T1.ravel()*dir1[1] + T2.ravel()*dir2[1]
    poss[:, 2] = p0[2] + T1.ravel()*dir1[2] + T2.ravel()*dir2[2]
    
    # Call the scan function using the computed positions.
    Es, Fs, Ps = mmff.scan(poss, bF=True, bP=True)
    
    # Reshape the energies into a 2D grid matching the T1, T2 shape.
    Egrid = Es.reshape(nscan1, nscan2)
    
    # Create a contour plot for the 2D energy scan.
    plt.figure()
    cp = plt.contourf(T1, T2, Egrid, levels=20, cmap="viridis")
    plt.colorbar(cp)
    plt.title(label)
    plt.xlabel(f"Scan parameter along ({dir1[0]}_{dir1[1]}_{dir1[2]}) direction")
    plt.ylabel(f"Scan parameter along ({dir2[0]}_{dir2[1]}_{dir2[2]}) direction")
    # plt.show()

        # Optionally, save the figure.
    if saveFig is not None:
        plt.savefig(saveFig)
    
    # Optionally, save the 2D scan data to a text file that you can later plot with gnuplot.
    # The file will contain three columns: t1, t2, and Energy.
    if saveData is not None:
        data = np.column_stack((T1.ravel(), T2.ravel(), Es.ravel()))
        np.savetxt(saveData, data, header="t1\tt2\tEnergy", comments="# ")
    plt.show()

#======== Body ###########

mmff.setVerbosity( verbosity=1, idebug=1 )

#mmff.init( xyz_name="data/xyz/pyridine", surf_name="data/NaCl_1x1_L2" )    
#mmff.init( xyz_name="data/xyz/nHexadecan_dicarboxylic", bMMFF=True  )     
# mmff.init( xyz_name="data/xyz/O", surf_name="data/xyz/NaCl_1x1_L3" )  
#mmff.init( xyz_name="data/xyz/H2O", surf_name="data/xyz/NaCl_1x1_L3" )    
# mmff.init( xyz_name="data/xyz/PTCDA", surf_name="data/xyz/NaCl_8x8_L3" )    
# mmff.init( xyz_name="data/xyz/PTCDA", surf_name="data/xyz/NaCl_8x8_L3_copy" )

# mmff.init( xyz_name="data/xyz/C.iz0", surf_name="data/xyz/C.iz0" )
# mmff.init( xyz_name="data/xyz/C.iz0", surf_name="data/xyz/Na.iz0" )
# mmff.init( xyz_name="data/xyz/C.iz0", surf_name="data/xyz/Cl.iz0" )
# mmff.init( xyz_name="data/xyz/C.iz0", surf_name="data/xyz/NaCl.iz0" )
# mmff.init( xyz_name="data/xyz/C_Cl_iz0", surf_name="data/xyz/NaCl.iz0" )

# mmff.init( xyz_name="data/xyz/Na+Cl-", surf_name="data/xyz/NaCl_8x8_L3_Coulumb" )

# mmff.init( xyz_name="data/xyz/PTCDA", surf_name="data/xyz/NaCl_coulomb.iz0" )

mmff.init( xyz_name="data/xyz/PTCDA_paolo", surf_name="data/xyz/NaCl.ptcda" )


# mmff.init( xyz_name="data/xyz/molNaCl_Na.iz0", surf_name="data/xyz/NaCl_coulomb.iz0" )
# mmff.init( xyz_name="data/xyz/molNaCl_Cl.iz0", surf_name="data/xyz/NaCl_coulomb.iz0" )
# mmff.init( xyz_name="data/xyz/molNaCl.ix0.iy0", surf_name="data/xyz/NaCl_coulomb.iz0" )

# mmff.init( xyz_name="data/xyz/2_atom_NaCl_mol", surf_name="data/xyz/2_atom_NaCl" )




'''
Debug :: REQs = [[1.4915     0.03606831 0.         0.        ]] For Na single atom substrate  Na.iz0.xyz
Debug :: REQs = [[1.9735     0.09921518 0.         0.        ]] For Cl single atom substrate  Cl.iz0.xyz


MolWorld_sp3::scan_rigid()[ia=0] pos(  0.0000,  0.0000,  7.1569) REQ(  1.9255,      0.06747763,  0.0000,  0.0000) PLQd(     21.77017471,      1.21202304,      0.00000000,      0.00000000)
MolWorld_sp3::scan_rigid()[ia=0] pos(  0.0000,  0.0000,  7.1569) REQ(  1.9255,      0.06747763,  0.0000,  0.0000) PLQd(     21.77017471,      1.21202304,      0.00000000,      0.00000000)

'''

mmff.getBuffs()

#print( "ffflags ", mmff.ffflags )

mmff.setSwitches( NonBonded=-1, MMFF=-1, SurfAtoms=0, GridFF=1 )


mmff.PLQs[:,0] = 0.0  # delete Pauli
mmff.PLQs[:,1] = 0.0  # delete London
# mmff.PLQs[:,2 ] = 0.0 # delete Coulomb (charges)

# scanPlot( nscan=10, span=(0.0,4.0), dir=(1.0,0.0,0.0), p0=(0.0,0.0,0.0),  label="E_x", saveFig="E_x_scan.png", saveData="E_x_scan.dat")
# scanPlot( nscan=1000, span=(0.0,4.0), dir=(0.0,1.0,0.0), p0=(0.0,0.0,0.0),  label="E_y", saveFig="E_y_scan.png", saveData="E_y_scan.dat" )

# scanPlot( nscan=100, span=(3.3,8.3), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0.0), label="Atop_C", saveFig="E_z_scan_Atop_C.png", saveData="E_z_scan_Atop_C.dat" )
# scanPlot( nscan=100, span=(2.9,8.3), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0.0), label="Atop_C on Na", saveFig="E_z_scan_Atop_Na.png", saveData="E_z_scan_Atop_Na.dat" )
# scanPlot( nscan=100, span=(3.3,8.3), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0.0), label="Atop_C on Cl ", saveFig="E_z_scan_Atop_Cl.png", saveData="E_z_scan_Atop_Cl.dat" )
# scanPlot( nscan=100, span=(0.8,6.3), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0.0), label="Atop_Na", saveFig="E_z_scan_Atop_Na_SYM.png", saveData="E_z_scan_Atop_Na_SYM.dat" )
# scanPlot( nscan=100, span=(2.3,8.3), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0.0), label="Atop_Na", saveFig="E_z_scan_Atop_Na_Ewald.png", saveData="E_z_scan_Atop_Na_Ewald.dat" )
# scanPlot( nscan=100, span=(7.9,18.3), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0.0), label="Atop_Na", saveFig="E_z_scan_Atop_Na_Ewald.png", saveData="E_z_scan_Atop_Na_Ewald.dat" )
# scanPlot( nscan=100, span=(3.1,8), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0.0), label="Atop_Cl", saveFig="E_z_scan_Atop_Cl_No_SYM.png", saveData="E_z_scan_Atop_Cl_No_SYM.dat" )

# scanPlot( nscan=100, span=(3,10), dir=(0.0,0.0,1.0), p0=(0.0,0.0,5.656854), label="PTCDA z0=0", saveFig="E_z_scan_PTCDA_Ewald.png", saveData="E_z_scan_PTCDA_Ewald.dat" )

# scanPlot( nscan=100, span=(3,10), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0), label="PTCDA", saveFig="E_z_scan_PTCDA_Ewald_new_Direct.png", saveData="E_z_scan_PTCDA_Ewald_new_Direct.dat" )
scanPlot( nscan=100, span=(3,10), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0), label="PTCDA", saveFig="E_z_scan_PTCDA_Ewald_new_trial.png", saveData="E_z_scan_PTCDA_Ewald_new_trial.dat" )
# scanPlot( nscan=100, span=(3,10), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0), label="PTCDA", saveFig="E_z_scan_PTCDA_Morse_new.png", saveData="E_z_scan_PTCDA_Morse_new.dat" )

# scanPlot( nscan=101, span=(3,8.3), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0.0), label="Atop_Na", saveFig="E_z_scan_Atop_Na_Ewald.png", saveData="E_z_scan_Atop_Na_Ewald.dat" )
# scanPlot( nscan=101, span=(3,8.3), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0.0), label="Atop_Cl", saveFig="E_z_scan_Atop_Cl_Ewald.png", saveData="E_z_scan_Atop_Cl_Ewald.dat" )

# scanPlot( nscan=100, span=(3.1,8.1), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0.0), label="Atop_Na", saveFig="E_z_scan_Atop_Na_Morse.png", saveData="E_z_scan_Atop_Na_Morse.dat" )
# scanPlot( nscan=100, span=(2.5,8.5), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0.0), label="Atop_Cl", saveFig="E_z_scan_Atop_Cl_Morse.png", saveData="E_z_scan_Atop_Cl_Morse.dat" )

# scanPlot( nscan=100, span=(0,20), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0.0), label="Atop_Na", saveFig="E_z_scan_2atom_NaCl_Ewald.png", saveData="E_z_scan_2atom_NaCl_Ewald.dat" )

# # For add pos0+ adjust(1.5) -4.9,4.5, 
# 2.35,8.5  
# # For sub pos0 only 8,12.157
#add pos0 only -3.5,4.5
# 1.9255
##################
# scanPlot2D(nscan1=41, nscan2=41,
#            span1=(0.0, 4.1), span2=(0.0, 4.1),
#            dir1=(1.0,0.0,0.0), dir2=(0.0,1.0,0.0),
#            p0=(0.0,0.0,3.4), label="E_xy for z=3.4",
#            saveFig="E_xy_scan_mol_NaCl.png", saveData="E_xy_scan_mol_NaCl.dat")

# scanPlot2D(nscan1=41, nscan2=61,
#            span1=(0.0, 4.1), span2=(2.3, 8.3),
#            dir1=(1.0,0.0,0.0), dir2=(0.0,0.0,1.0),
#            p0=(0.0,0.0,0), label="E_xz for Y=3.4",
#            saveFig="E_xz_scan_C.png", saveData="E_xz_scan_C.dat")

# scanPlot2D(nscan1=41, nscan2=61,
#            span1=(0.0, 4.1), span2=(2.3, 8.3),
#            dir1=(0.0,1.0,0.0), dir2=(0.0,0.0,1.0),
#            p0=(0.0,0.0,0), label="E_xz for Y=3.4",
#            saveFig="E_yz_scan_C.png", saveData="E_yz_scan_C.dat")           

# plt.legend()
# plt.grid()
# plt.show()
# print(mmff.PLQs.shape)



# Initialize the system. Make sure you call mmff.init() and then mmff.getBuffs()
# mmff.setVerbosity(verbosity=1, idebug=1)
# mmff.init(xyz_name="data/xyz/PTCDA", surf_name="data/xyz/NaCl_8x8_L3_copy")
# mmff.getBuffs()

# (Optionally) set up any switches, etc.
# mmff.setSwitches(NonBonded=-1, MMFF=-1, SurfAtoms=0, GridFF=1)

# Now, call our new saveXSF_geometry function.
# We assume that: 
#   - The grid data is stored in an array registered in the C++ global dictionary with key "BsplinePaul_pbc"
#   - mmff.natoms holds the total number of atoms of the substrate plus PTCDA.
#   - mmff.AtomType is a pointer or array of atom types (or use mmff.params.atypes if that’s your interface)
#   - mmff.fapos holds the geometry (atomic positions)

# It might be that the name "BsplinePaul_pbc" is not registered; adjust the key if needed.
# grid_array = mmff.getArrayPointer("BsplinePaul_pbc")
# if grid_array is None:
#     raise RuntimeError('Grid array "BsplinePaul_pbc" not found. Ensure the grid is built.')

# # Save to XSF file.
# ret = mmff.saveXSF_geometry("rigid_scan.xsf", grid_array, pitch=1, offset=0, 
#                               natoms=mmff.natoms,
#                               atypes=mmff.params.atypes if hasattr(mmff, "params") else None,
#                               fapos=mmff.fapos, bPrimCoord=True)
# if ret != 0:
#     print("XSF file saved successfully!")
# else:
#     print("Error saving XSF file.")

# # Optionally, run the dynamics or further analysis
# mmff.run()



# Save the grid and geometry to an XSF file:
# mmff.saveXSF("rigid_scan.xsf", mmff.getArrayPointer("BsplinePaul_pbc"), 1, 0, mmff.natoms, mmff.AtomType, mmff.fapos)
# mmff.saveXSF("rigid_scan.xsf", mmff.getArrayPointer("BsplinePaul_pbc"))

# E = 0.0
# mmff.run(omp=True, nstepMax=20)
# print( "E=", E )
#for i in range(200):
#mmff.addSnapshot()
# with open("gopt_trajectory.xyz", "w") as file:pass
# mmff.printDatabase()

#mmff.init( xyz_name="dataxyz//pyridine", surf_name="dataxyz/NaCl_1x1_L2", bMMFF=False, gridStep=-1 )  # without gridFF
#mmff.getBuffs()
#mmff.eval()
#mmff.relax(1000)
#print( "FORCES:\n mmff.fapos:\n ", mmff.fapos )
#mmff.plot(bForce=True, Fscale=10.0 )
#plt.show()
# exit(0)


# '''
# mmff.init()
# #mmff.init_params( "data/AtomTypes.dat", "data/BondTypes.dat", "data/AngleTypes.dat" )
# #mmff.insertSMILES("CC");
# #mmff.insertSMILES("C=C");
# #mmff.insertSMILES("C#C");
# #mmff.insertSMILES("C#CCN=C", True );
# #mmff.insertSMILES("C1#CCN=C1", True );
# #mmff.insertSMILES("C=C1NC#CC1CO", True, True );

# #mmff.initWithSMILES( "C=C1NC#CC1CO" )
# mmff.initWithSMILES( "C=C" )
# mmff.getBuffs()
# mmff.relax(1000)
# mmff.plot()
# plt.show()

# exit()
# '''

# '''
# # ======== Oritent Molecule
# xyzs,Zs,enames,qs = au.loadAtomsNP( "data/xyz/Benzene_deriv.xyz" )
# au.orient( 2, (5,2), (1,3), xyzs, bFlipXZ=True )
# au.saveXYZ( enames, xyzs, "data/xyz/Benzene_deriv_.xyz", qs=qs, Rs=None )
# plt.plot( xyzs[:,0],xyzs[:,1], "o" )
# plt.axis('equal')
# plt.show()
# exit()
# '''

# '''
# # ============== C2H4,xyz
# #mmff.initWithMolFile( "C2H4.xyz", bNonBonded=False, bOptimizer=True)
# #mmff.printBuffNames()
# #mmff.getBuffs() #;print( mmff.ndims )
# #mmff.eval()
# #mmff.relax(1000, bWriteTrj=True )
# #Es=mmff.scanRotation( [1,4,5], 0, 0,1, np.pi*2, 100, bWriteTrj=True)   ;print("Es=", Es)
# #plt.plot(Es)
# #print( "Es(Etot,Eb,Ea,Eps,EppT,EppI):", mmff.Es )
# #nsel = mmff.splitAtBond(6-1)  ;print( "split to:\n", mmff.selection[:nsel],"\n", mmff.selection[nsel:] )
# '''



# # ============== Benzene_deriv.xyz
# #mmff.initWithMolFile( "data/xyz/Benzene_deriv.xyz", bNonBonded=False, bOptimizer=True)


# # mmff.initWithMolFile( "data/xyz/PTCDA.xyz", bNonBonded=False, bOptimizer=True)
# # mmff.getBuffs() #;print( mmff.ndims )

# #nsel = mmff.splitAtBond(5)  ;print( "split to:\n", mmff.selection[:nsel],"\n", mmff.selection[nsel:] )
# #nsel = mmff.splitAtBond(6)  ;print( "split to:\n", mmff.selection[:nsel],"\n", mmff.selection[nsel:] )
# #nsel = mmff.splitAtBond(10)  ;print( "split to:\n", mmff.selection[:nsel],"\n", mmff.selection[nsel:] )
# #nsel = mmff.splitAtBond(2)  ;print( "split to:\n", mmff.selection[:nsel],"\n", mmff.selection[nsel:] )
# #nsel = mmff.splitAtBond(4)  ;print( "split to:\n", mmff.selection[:nsel],"\n", mmff.selection[nsel:] )
# #print( "nsel ", nsel, len(mmff.selection)-nsel )
# #Es=mmff.scanRotation( 1, 1,11, np.pi*2, 100, bWriteTrj=True, _0=1) ;plt.plot(Es) ;print("Es=", Es)

# #Es = mmff.scanBondRotation( 6, np.pi*2, 100, bWriteTrj=True );  plt.plot(Es)
# #Es = mmff.scanBondRotation( 2, np.pi*2, 100, bWriteTrj=True );  plt.plot(Es); plt.grid()


# #mmff.eval()
# #mmff.relax(1000, Ftol=1e-4, bWriteTrj=True )
# #Es=mmff.scanRotation( 1, 1,11, np.pi*2, 100, sel=[11,13,14,20]+[29,30,31,32], bWriteTrj=True, _0=1) ;plt.plot(Es) ;print("Es=", Es)

# # plt.figure()
# # mmff.plot()
# # #mmff.plot_selection( mmff.selection[:nsel] )
# # #mmff.plot_selection( [1,2,3] )



# file_path = "/home/indranil/git/FireCore/cpp/common_resources/NaCl_8x8_L3_copy/Bspline_PLQd.npy"
# grid = np.load(file_path)

# print("Shape of grid:", grid.shape)  # expected (Nx, Ny, Nz, 3)

# channel = 0  # select the potential channel to view (e.g., Morse potential)

# # For this example we assume the grids share the same shape.
# Nx, Ny, Nz, _ = grid.shape

# # Use fixed indices (you can change these to whichever slice you want)
# mid_z = 10 + Nz//2  # for horizontal (xy) slice at constant z
# mid_x = 10 + Nx//2  # for vertical_x (yz) slice at constant x
# mid_y = 10 + Ny//2  # for vertical_y (xz) slice at constant y

# # Extract slices from the grid
# slice_xy = grid[:, :, mid_z, channel]
# slice_yz = grid[mid_x, :, :, channel]
# slice_xz = grid[:, mid_y, :, channel]


# # ---------------------------
# # Retrieve grid parameters and atom positions from the C++ extension
gff_shift0, gff_pos0, gff_cell, gff_dCell, gff_natoms, gff_natoms_ = mmff.get_gridFF_info()
# print("gff_pos0:", gff_pos0)
print("gff_shift0:", gff_shift0)
# print("gff_cell:\n", gff_cell)
# print("gff_dCell:\n", gff_dCell)

# # For scaling purposes (we assume dCell is diagonal and defines the grid spacing)
# dx = gff_dCell[0, 0]
# dy = gff_dCell[1, 1]
# dz = gff_dCell[2, 2]

# # Compute physical extents for each slice.
# # Here our coordinate system will be shifted so that the grid origin becomes zero.
# extent_xy = [0, Nx * dx, 0, Ny * dy]     # x spans along axis0 and y axis1
# extent_yz = [0, Ny * dy, 0, Nz * dz]       # for a slice taken at constant x; horizontal axis: y, vertical: z
# extent_xz = [0, Nx * dx, 0, Nz * dz]       # for a slice taken at constant y; horizontal: x, vertical: z
# # extent_xy = [ gff_pos0[0], gff_pos0[0] + Nx * dx, gff_pos0[1], gff_pos0[1] + Ny * dy ]
# # extent_yz = [ gff_pos0[1], gff_pos0[1] + Ny * dy, gff_pos0[2], gff_pos0[2] + Nz * dz ]
# # extent_xz = [ gff_pos0[0], gff_pos0[0] + Nx * dx, gff_pos0[2], gff_pos0[2] + Nz * dz ]

# # ---------------------------
# # Get atom positions (substate and molecule)
# substrate_apos, molecule_apos = mmff.get_atom_positions()

# # --- Scaling and Shifting Atom Positions ---
# # The grid (and substrate) has a global origin gff_pos0.
# # For plotting in a coordinate system where the lower left is 0, we subtract gff_pos0.
# # Additionally, the molecule atoms are offset by gff_shift0.
# substrate_apos_scaled = substrate_apos #+ gff_pos0
# molecule_apos_scaled  = molecule_apos  - gff_shift0  #+ gff_pos0

# # ---------------------------
# # Plotting: Create a figure with three subplots for grid slices.
# fig, axs = plt.subplots(1, 3, figsize=(18, 6))

# # Horizontal slice (xy plane) at constant z
# im0 = axs[0].imshow(slice_xy.T, extent=extent_xy, origin="lower", cmap="viridis", aspect="auto")
# axs[0].set_title(f"Horizontal (xy) slice at z = {mid_z*dz}")
# axs[0].set_xlabel("x (physical dimesion)")
# axs[0].set_ylabel("y (physical dimension)")
# fig.colorbar(im0, ax=axs[0], fraction=0.046, pad=0.04)
# # Overlay atom positions (x-y)
# axs[0].scatter(substrate_apos_scaled[:, 0], substrate_apos_scaled[:, 1],
#                s=10, color='red', label='Substrate')
# axs[0].scatter(molecule_apos_scaled[:, 0], molecule_apos_scaled[:, 1],
#                s=10, color='blue', label='Molecule')
# axs[0].legend()

# # Vertical slice (yz plane) taken at constant x (mid_x)
# im1 = axs[1].imshow(slice_yz.T, extent=extent_yz, origin="lower", cmap="viridis", aspect="auto")
# axs[1].set_title(f"Vertical (yz) slice at x = {mid_x*dx}")
# axs[1].set_xlabel("y (physical dimension)")
# axs[1].set_ylabel("z (physical dimension)")
# fig.colorbar(im1, ax=axs[1], fraction=0.046, pad=0.04)
# # Overlay atom positions in y-z plane
# axs[1].scatter(substrate_apos_scaled[:, 1], substrate_apos_scaled[:, 2],
#                s=10, color='red')
# axs[1].scatter(molecule_apos_scaled[:, 1], molecule_apos_scaled[:, 2],
#                s=10, color='blue')

# # Vertical slice (xz plane) taken at constant y (mid_y)
# im2 = axs[2].imshow(slice_xz.T, extent=extent_xz, origin="lower", cmap="viridis", aspect="auto")
# axs[2].set_title(f"Vertical (xz) slice at y = {mid_y*dy}")
# axs[2].set_xlabel("x (physical dimesion)")
# axs[2].set_ylabel("z (physical dimension)")
# fig.colorbar(im2, ax=axs[2], fraction=0.046, pad=0.04)
# # Overlay atom positions in x-z plane
# axs[2].scatter(substrate_apos_scaled[:, 0], substrate_apos_scaled[:, 2],
#                s=10, color='red')
# axs[2].scatter(molecule_apos_scaled[:, 0], molecule_apos_scaled[:, 2],
#                s=10, color='blue')

# plt.tight_layout()
# plt.show()

exit()