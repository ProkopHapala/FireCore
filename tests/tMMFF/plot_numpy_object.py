# import numpy as np
# import matplotlib.pyplot as plt

# # Change these file names to your actual grid file names
# # file1 = "/home/indranil/git/FireCore/cpp/common_resources/NaCl_8x8_L3_copy/no_Shift_Bspline_PLQd.npy"
# file2 = "/home/indranil/git/FireCore/cpp/common_resources/NaCl_8x8_L3_copy/Bspline_PLQd.npy"
# file3 = "/home/indranil/git/FireCore/cpp/common_resources/NaCl_8x8_L3_copy/no_Shift_Bspline_PLQd.npy"
   
# # Load the grid force field arrays
# # grid1 = np.load(file1)
# grid2 = np.load(file2)
# grid3 = np.load(file3)

# print("Shapes:",  grid2.shape, grid3.shape)

# # Let’s assume the arrays are shaped (Nx, Ny, Nz, 3): one channel per potential
# # You can now compare a given channel, for example channel 0 (e.g., Morse potential)

# slice_idx = 4 # For example, take the first z-slice if Nz corresponds to z
# channel = 0

# # diff12 = np.abs(grid1[:,:,slice_idx,channel] - grid2[:,:,slice_idx,channel])
# # diff13 = np.abs(grid1[:,:,slice_idx,channel] - grid3[:,:,slice_idx,channel])
# # diff23 = np.abs(grid2[:,:,slice_idx,channel] - grid3[:,:,slice_idx,channel])
# diff2 = (grid2[:,:,slice_idx,channel])
# diff3 = (grid3[:,:,slice_idx,channel])


# plt.figure(figsize=(15, 5))

# # plt.subplot(1,3,1)
# # # plt.imshow(diff12, origin="lower", cmap="viridis")
# # plt.colorbar()
# # plt.title("Abs diff: File1 - File2")

# plt.subplot(1,3,2)
# plt.imshow(diff2, origin="lower", cmap="viridis")
# plt.colorbar()
# plt.title("Abs diff: File1 - File3")

# plt.subplot(1,3,3)
# plt.imshow(diff3, origin="lower", cmap="viridis")
# plt.colorbar()
# plt.title("Abs diff: File2 - File3")

# plt.tight_layout()
# plt.show()

# # Optionally, compute some summary statistics
# # print("File1 vs File2:")
# # print("   Mean difference:", np.mean(diff12))
# # print("   Max difference:", np.max(diff12))
# # print("File1 vs File3:")
# # print("   Mean difference:", np.mean(diff13))
# # print("   Max difference:", np.max(diff13))
# # print("File2 vs File3:")
# print("   Mean difference:", np.mean(diff23))
# print("   Max difference:", np.max(diff23))


# ### To see the PLQD numpy object

# # import numpy as np
# # import matplotlib.pyplot as plt

# # # Load the Bspline_PLQd.npy file
# # data = np.load("Bspline_PLQd.npy")
# # print("Shape:", data.shape)

# # # Plotting based on data dimensions:
# # # If the data is 2D, we can directly use imshow.
# # # For larger or different dimensions, adjust the plotting accordingly.
# # if data.ndim == 2:
# #     plt.imshow(data, cmap='viridis', aspect='auto')
# #     plt.colorbar(label='Value')
# #     plt.title("Bspline_PLQd.npy 2D Plot")
# # elif data.ndim == 1:
# #     plt.plot(data)
# #     plt.title("Bspline_PLQd.npy 1D Plot")
# # else:
# #     # For 3D arrays, pick a slice along the first dimension
# #     plt.imshow(data[0], cmap='viridis', aspect='auto')
# #     plt.colorbar(label='Value')
# #     plt.title("Bspline_PLQd.npy Slice Plot (first slice)")

# # plt.xlabel("X-axis")
# # plt.ylabel("Y-axis")
# # plt.tight_layout()
# # plt.show()



# import numpy as np
# import matplotlib.pyplot as plt
# import ipywidgets as widgets
# from IPython.display import display

# # Load the grid force field arrays (replace with your actual file paths)
# file2 = "/home/indranil/git/FireCore/cpp/common_resources/NaCl_8x8_L3_copy/Bspline_PLQd.npy"
# file3 = "/home/indranil/git/FireCore/cpp/common_resources/NaCl_8x8_L3_copy/no_Shift_Bspline_PLQd.npy"
# grid2 = np.load(file2)
# grid3 = np.load(file3)

# # Select the channel you're interested in
# channel = 0  # Example: Morse potential

# # Create interactive plot
# def plot_3d_slice(slice_idx):
#     fig, axes = plt.subplots(1, 2, figsize=(12, 5))

#     im1 = axes[0].imshow(grid2[:, :, slice_idx, channel], origin="lower", cmap="viridis")
#     axes[0].set_title(f"grid2, Slice {slice_idx}")
#     fig.colorbar(im1, ax=axes[0])

#     im2 = axes[1].imshow(grid3[:, :, slice_idx, channel], origin="lower", cmap="viridis")
#     axes[1].set_title(f"grid3, Slice {slice_idx}")
#     fig.colorbar(im2, ax=axes[1])

#     plt.tight_layout()
#     plt.show()


# # Create a slider to control the slice index
# slice_slider = widgets.IntSlider(
#     min=0, max=grid2.shape[2] - 1, step=1, description="Slice:", continuous_update=False
# )

# # Display the interactive plot and the slider
# widgets.interactive(plot_3d_slice, slice_idx=slice_slider)


# # --- Difference Plot (Interactive)---
# def plot_diff(slice_idx):
#     diff = grid2[:, :, slice_idx, channel] - grid3[:, :, slice_idx, channel]
#     plt.figure(figsize=(6, 5))
#     plt.imshow(diff, origin="lower", cmap="viridis")
#     plt.colorbar()
#     plt.title(f"Difference (grid2 - grid3), Slice {slice_idx}")
#     plt.show()
#     print("   Mean difference:", np.mean(diff))
#     print("   Max difference:", np.max(diff))



# # Create a slider to control the slice index for the difference
# slice_slider_diff = widgets.IntSlider(
#     min=0, max=grid2.shape[2] - 1, step=1, description="Slice:", continuous_update=False
# )

# # Display the interactive difference plot and slider
# widgets.interactive(plot_diff, slice_idx=slice_slider_diff)

# import numpy as np
# from vispy import app, scene, color

# # Load the grid force field arrays (adjust the file path to your grid file)
# file2 = "/home/indranil/git/FireCore/cpp/common_resources/NaCl_8x8_L3_copy/Bspline_PLQd.npy"
# grid = np.load(file2)

# print("Shape of grid:", grid.shape)  # e.g. (320, 320, 200, 3)

# channel = 0  # Choose a potential channel

# Nx, Ny, Nz, _ = grid.shape

# # Create coordinate arrays
# x = np.arange(Nx)
# y = np.arange(Ny)
# z = np.arange(Nz)
# X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
# data = grid[:,:,:,channel]  # extracting the data

# # Flatten the arrays
# X_flat = X.flatten()
# Y_flat = Y.flatten()
# Z_flat = Z.flatten()
# data_flat = data.flatten()

# # Normalize data values for colors
# norm = (data_flat - data_flat.min()) / (data_flat.max() - data_flat.min())
# # Fixed: use .map() method on the colormap object
# colors = color.get_colormaps()['viridis'].map(norm)

# # Create a VisPy canvas
# canvas = scene.SceneCanvas(keys='interactive', show=True, bgcolor='white')
# view = canvas.central_widget.add_view()
# view.camera = scene.cameras.TurntableCamera(fov=45)

# # Add scatter plot (GPU accelerated via OpenGL)
# scatter = scene.visuals.Markers()
# scatter.set_data(np.column_stack((X_flat, Y_flat, Z_flat)), face_color=colors, size=3)
# view.add(scatter)

# if __name__ == '__main__':
#     app.run()


############################ To plot the grid object slices ###########################################################
# import numpy as np
# import matplotlib.pyplot as plt

# # Replace with your actual file path.
# file_path = "/home/indranil/git/FireCore/cpp/common_resources/NaCl_8x8_L3_copy/no_Shift_Bspline_PLQd.npy"
# grid = np.load(file_path)
# print("Shape of grid:", grid.shape)  # Expected (Nx, Ny, Nz, 3)

# channel = 0  # Select the potential channel to view (e.g., Morse potential)

# # Get grid dimensions
# Nx, Ny, Nz, _ = grid.shape

# # Calculate the middle index for each dimension
# # mid_z = Nz // 2   # for horizontal (xy) slice at constant z
# # mid_x = Nx // 2   # for vertical_x (yz) slice at constant x
# # mid_y = Ny // 2   # for vertical_y (xz) slice at constant y

# mid_z = 199   # for horizontal (xy) slice at constant z
# mid_x = 308   # for vertical_x (yz) slice at constant x
# mid_y = 310   # for vertical_y (xz) slice at constant y

# # Extract slices:
# # Horizontal cut (xy plane) at constant z = mid_z.
# slice_xy = grid[:, :, mid_z, channel]
# # Vertical cut (yz plane) at constant x = mid_x.
# slice_yz = grid[mid_x, :, :, channel]
# # Vertical cut (xz plane) at constant y = mid_y.
# slice_xz = grid[:, mid_y, :, channel]

# # Create subplots for the three cuts.
# fig, axs = plt.subplots(1, 3, figsize=(18, 6))

# # Horizontal slice: xy plane.
# im0 = axs[0].imshow(slice_xy, origin="lower", cmap="viridis", aspect="auto")
# axs[0].set_title(f"Horizontal Cut (xy) at z = {mid_z}")
# axs[0].set_xlabel("x")
# axs[0].set_ylabel("y")
# fig.colorbar(im0, ax=axs[0], fraction=0.046, pad=0.04)

# # Vertical slice: yz plane.
# im1 = axs[1].imshow(slice_yz, origin="lower", cmap="viridis", aspect="auto")
# axs[1].set_title(f"Vertical Cut (yz) at x = {mid_x}")
# axs[1].set_xlabel("y")
# axs[1].set_ylabel("z")
# fig.colorbar(im1, ax=axs[1], fraction=0.046, pad=0.04)

# # Vertical slice: xz plane.
# im2 = axs[2].imshow(slice_xz, origin="lower", cmap="viridis", aspect="auto")
# axs[2].set_title(f"Vertical Cut (xz) at y = {mid_y}")
# axs[2].set_xlabel("x")
# axs[2].set_ylabel("z")
# fig.colorbar(im2, ax=axs[2], fraction=0.046, pad=0.04)

# plt.tight_layout()
# plt.show()


############################ To plot the grid objects slices and compaire the difference ###########################################################
# import numpy as np
# import matplotlib.pyplot as plt

# # Replace with your actual file paths.
# file_path1 = "/home/indranil/git/FireCore/cpp/common_resources/NaCl_8x8_L3_copy/no_Shift_Bspline_PLQd.npy"
# file_path2 = "/home/indranil/git/FireCore/cpp/common_resources/NaCl_8x8_L3/cpu/Bspline_PLQd.npy"

# grid1 = np.load(file_path1)
# grid2 = np.load(file_path2)

# print("Shape of grid1:", grid1.shape)  # expected (Nx, Ny, Nz, 3)
# print("Shape of grid2:", grid2.shape)  # expected (Nx, Ny, Nz, 3)

# channel = 0  # select the potential channel to view (e.g., Morse potential)

# # For this example we assume the grids share the same shape.
# Nx, Ny, Nz, _ = grid1.shape

# # Use fixed indices (you can change these to whichever slice you want)
# mid_z = 199  # for horizontal (xy) slice at constant z
# mid_x = 308  # for vertical_x (yz) slice at constant x
# mid_y = 310  # for vertical_y (xz) slice at constant y

# # Extract slices from grid1
# slice_xy_1 = grid1[:, :, mid_z, channel]
# slice_yz_1 = grid1[mid_x, :, :, channel]
# slice_xz_1 = grid1[:, mid_y, :, channel]

# # Extract slices from grid2
# slice_xy_2 = grid2[:, :, mid_z, channel]
# slice_yz_2 = grid2[mid_x, :, :, channel]
# slice_xz_2 = grid2[:, mid_y, :, channel]

# # Compute differences of the slices (grid2 - grid1)
# slice_xy_diff = slice_xy_2 - slice_xy_1
# slice_yz_diff = slice_yz_2 - slice_yz_1
# slice_xz_diff = slice_xz_2 - slice_xz_1

# # Create a figure with 3 rows and 3 columns:
# # Row 0: grid1 slices; Row 1: grid2 slices; Row 2: difference slices.
# fig, axs = plt.subplots(3, 3, figsize=(18, 15))

# # ---------------------------
# # Row 0: grid1
# # Horizontal slice: xy plane.
# im0 = axs[0, 0].imshow(slice_xy_1, origin="lower", cmap="viridis", aspect="auto")
# axs[0, 0].set_title(f"Grid1: Horizontal (xy) at z = {mid_z}")
# axs[0, 0].set_xlabel("x")
# axs[0, 0].set_ylabel("y")
# fig.colorbar(im0, ax=axs[0, 0], fraction=0.046, pad=0.04)

# # Vertical slice: yz plane.
# im1 = axs[0, 1].imshow(slice_yz_1, origin="lower", cmap="viridis", aspect="auto")
# axs[0, 1].set_title(f"Grid1: Vertical (yz) at x = {mid_x}")
# axs[0, 1].set_xlabel("y")
# axs[0, 1].set_ylabel("z")
# fig.colorbar(im1, ax=axs[0, 1], fraction=0.046, pad=0.04)

# # Vertical slice: xz plane.
# im2 = axs[0, 2].imshow(slice_xz_1, origin="lower", cmap="viridis", aspect="auto")
# axs[0, 2].set_title(f"Grid1: Vertical (xz) at y = {mid_y}")
# axs[0, 2].set_xlabel("x")
# axs[0, 2].set_ylabel("z")
# fig.colorbar(im2, ax=axs[0, 2], fraction=0.046, pad=0.04)

# # ---------------------------
# # Row 1: grid2
# # Horizontal slice: xy plane.
# im3 = axs[1, 0].imshow(slice_xy_2, origin="lower", cmap="viridis", aspect="auto")
# axs[1, 0].set_title(f"Grid2: Horizontal (xy) at z = {mid_z}")
# axs[1, 0].set_xlabel("x")
# axs[1, 0].set_ylabel("y")
# fig.colorbar(im3, ax=axs[1, 0], fraction=0.046, pad=0.04)

# # Vertical slice: yz plane.
# im4 = axs[1, 1].imshow(slice_yz_2, origin="lower", cmap="viridis", aspect="auto")
# axs[1, 1].set_title(f"Grid2: Vertical (yz) at x = {mid_x}")
# axs[1, 1].set_xlabel("y")
# axs[1, 1].set_ylabel("z")
# fig.colorbar(im4, ax=axs[1, 1], fraction=0.046, pad=0.04)

# # Vertical slice: xz plane.
# im5 = axs[1, 2].imshow(slice_xz_2, origin="lower", cmap="viridis", aspect="auto")
# axs[1, 2].set_title(f"Grid2: Vertical (xz) at y = {mid_y}")
# axs[1, 2].set_xlabel("x")
# axs[1, 2].set_ylabel("z")
# fig.colorbar(im5, ax=axs[1, 2], fraction=0.046, pad=0.04)

# # ---------------------------
# # Row 2: differences: (grid2 - grid1)
# # Horizontal slice difference.
# im6 = axs[2, 0].imshow(slice_xy_diff, origin="lower", cmap="bwr", aspect="auto")
# axs[2, 0].set_title(f"Difference (xy) at z = {mid_z}")
# axs[2, 0].set_xlabel("x")
# axs[2, 0].set_ylabel("y")
# fig.colorbar(im6, ax=axs[2, 0], fraction=0.046, pad=0.04)

# # Vertical slice difference: yz plane.
# im7 = axs[2, 1].imshow(slice_yz_diff, origin="lower", cmap="bwr", aspect="auto")
# axs[2, 1].set_title(f"Difference (yz) at x = {mid_x}")
# axs[2, 1].set_xlabel("y")
# axs[2, 1].set_ylabel("z")
# fig.colorbar(im7, ax=axs[2, 1], fraction=0.046, pad=0.04)

# # Vertical slice difference: xz plane.
# im8 = axs[2, 2].imshow(slice_xz_diff, origin="lower", cmap="bwr", aspect="auto")
# axs[2, 2].set_title(f"Difference (xz) at y = {mid_y}")
# axs[2, 2].set_xlabel("x")
# axs[2, 2].set_ylabel("z")
# fig.colorbar(im8, ax=axs[2, 2], fraction=0.046, pad=0.04)

# plt.tight_layout()
# plt.show()

'''''''''
## Problem Solved need to perform the git commit accordigly:

# I was assigned to prepare a detailed documentation, which is not compleated yet.

# Problem Solved 0: (Runing all the test and get familiar with the test run)
Tested all the test run from ~/git/FireCore/tests/tMMFF folder  and produced output and commented errors and also refined some results in the test run. All the outputs are there in a 
md file called ~/git/FireCore/tests/tMMFF/tMMFF_Test_Runing_QA.md, Under few refinement results the correspoding chages are mentiond. the commit message should follow that.  


# Problem Solved 1: (For the grid dimension error fft pan error)
We had a problem with the grid dimension; If the dimetion is not a well factoriable number then it will show the fft pan error, whcih was due to the fixed grid step size and 
as a result it was dependent on the size of the substrate, which is solved by making the grid step size adjustable to make the overall grid size number a weel factoriable number (2,3,5).
To do that we have introduced a new function called ''next_nice'' defined in the clUtils.py file(pyBall/OCL). Also called in GridFF.py file (pyBall/OCL). 


# Problem Solved 2: (Shift0 issue with the positiong of the molecule atoms , substrate atoms, and the grid postion with respect to the global coordiates)

"Prokop Hapala
Fri, Feb 7, 5:00PM (5 days ago)
to me

""The meaning of grid.pos0 and gridFF.shift0 is like this:

1)  grid.pos0 is used to shift the origin (lower left corner) of the Grid in the world-space after it is created projected. Therefore it is used in all the interpolation functions. 
2)  gridFF.shift0 is used just to shift the atoms of the substrate before they are projected to the grid. This also means we need to use it when calculating direct interaction of molecule with substrate atoms (not using grid bGridFF=false)

However it seems we are not useing it consistently in all functions 
 - when we calculate it on GPU using pyOpenCL by external python script we need to make sure we are using the same shift0 as we have inside MolWorld in C++ 
( both MMFF_lib and MolGUI use the same MolWorld and GridFF object so the shift0 should be same there)
 - when we calculate GridFF on CPU inside GridFF.h it seems we are also not using shift0. Probabaly I did this to make it easier to compare pyOpenCL (GPU) 
 and C++ (CPU) creation of GridFF, but I forgot that it will make it inconsistent with direct evaluation (bGridFF=false)" " "

    2 (a): 
    Now only the inconsistency in GUI is solved, we can see the effect of shift0 as defined in the ~/git/FireCore/cpp/common/molecular/MolWorld_sp3.h:965 while runing the GUI after calculating
    GRID using the ocl code ~/git/FireCore/tests/tMMFF/run_test_GridFF_ocl_new.py, What I mean to say is that we can see how the molecule (PTCDA) is aligning itself with the grid and substrate.
    There is no inconsistancy between the Grid and the substarate atoms. However, while generating the grid using run_test_GridFF_ocl_new.py -> ocl_GridFF_new.py, 
    the shift0 is not affecting the grid generation, what I mean by this is that even though we add some shift +/- some value, it will also tacking into account while generating the grid
    by the function find_top and the variable z0 which will incorporate the effect of shift0 in grid generation. Thefore, the grid and the substarte atoms are in same positon irrespective of 
    the shift0, for the ocl code to generate Grid.
    So we can conclude that the shift0 is affecting the GUI after generating the Grid using GPUcode and we can not make it hardcoded in the code itself, need to be adjusted in the GUI. 

    
    Not Solved yet:
    2 (b):
    There is some issue with the symmetrization of the substarte atoms, in the GUI, which is temporaryly fixed by not using the symmetrization function.
    Also how to calculate results like scan and all PES from the grid, that I do not know yet.

    2 (c):
    While generating the grid using cpu code, there is some issue with the simetrization of substrate may be, because it does not look like as the grid generated by the gpu code.
    However, I do not know what is the status of the shift0 here in  the cpu code. which can be adjested in 
        "but when we create the grid in CPU ( initGrid -> tryLoad_new -> makeGridFF_Bspline_d -> tryLoadGridFF_potentials -> evalBsplineRef -> evalGridFFPoint_Mors )
        evalMorsePBC_PLQ_sym  -> evalMorsePBC_PLQ->  evalGridFFPoint"
    
    2 (d):
    While using the CLI mode for the runing Firecore after grid generation usimg ocl code, we can not see the effect of shift0, which is not the case in the GUI. 
    To address this we need to define a function which can show the grid sclice,the substrate atoms and the molecules atoms. 

    This is partially solved; need to veryfy the aligment with prokop, what was done till now is that modified three files one function is defined in MMFF_lib.cpp, called get_grid_info () which will pass the grid information to he the python wrapper MMFF.py.
    Similarly other two functions to retrive the information about the substrate and the molecule atoms respectively. Then in run.py we are retrivng all the informations and loding the grid from the numpy file and ploting the slice of the grid with the substrate and the molecule atoms.
    The shift0 is also taken into account while plotting the atoms, retrived from the library which has a problem of recompileing the library after changing the value at MolWorld_sp3.h:965. This can be solved if we propagate the value of shift0 from the scan I mean the run.py 
    to the MolWorld_sp3.h:965, which is not done yet.  

    2 (e): 
    The rigid scan results are not dependent on the shift0. 
    It got a long story : Now it is dependent on the shift0. Now how it is propagating is as follows:

    while scaning in the run.py -> MMFF.py scan which will use the library MMFF_lib.cpp where this scan will invoke the rigid scan through scan_rigid , which is defined in the MolWorld_sp3.h:2465, which further invoke eval_no_omp in the same file at 2024. which further calls 
    addAtom from the GridFF.h:450 for the case of if GridFF is true otherwise call evalMorsePBC_sym. 

    The adAtom then call getForce_Bspline @354  in GridFF, there the actual calculation of interaction is done, where the atomic position is taken and projected into the global coordinate, through add/sub pos0. Now this pos0 need to be added or substracted depends on the 
    sign of pos0 if the z coordinate is positive then add pos0 otherwise substract pos0. 

    Now if we wan to add the shift0 along with this then also the sign of shift0 matters, Therefore to make it more genric we should take the value of shift0 as +ve number then add it with pos0. 

    I think it doesnot matter what ever the sign of the poss0 is  if the molecule's z coordinates is set as zero then it will always be p.add (pos0+shift0), ad shift0 should always be +ve. i.e lets say if the molecules is at z =0 and pos0=-5 ad shift0=2 the the final position will be 
    p.add (pos0+shift0) = p.add (-5+2) = -3. if the pos0 is +ve 5 then also p.add (pos0+shift0) = p.add (5+2) = 7. which we wanted 

    
    Now we need to do couple of fixes, this value of shift0 will be taken from the run.py or where the scan is initiated not from the GriddFF. 
    Therefore we need to create a function in the python wrapper which should propagate from pytho to final cpp grid calculations.

     For calculation of interaction when Gridd is true: 
     run.py --> MMFF.py --> MMFF_lib.cpp --> scan_rigid --> eval_no_omp --> GridFF.h:450 addAtom --> getForce_Bspline -->  p.add (pos0+shift0).

     For Direct calculation of interaction:
     run.py --> MMFF.py --> MMFF_lib.cpp --> scan_rigid --> eval_no_omp --> GridFF.h:1249 evalMorsePBC_sym --> evalMorsePBC  with only taking apos_ (symmetrised)-->  
     
     (i) mostly the dp0 needs to me modified to make the correction of Shift0 and the previous fix of p.add 
     
     (ii) there is also one fixed required, about the symmetrization fuction. 
    
    Also this should reflect in GUI also. 


    There is no need to chnage Shift0 no need atall and what was implemneted previously is fine, I mean the p.sub it should be the case not the addition. due to the vector representation of the position vector.
    Now the positito of the molecule is clear it should be the same as the position of the substrate along z direction , I mean the ztop of the substarte should be the same as the z coordinate of the molecule, then start scanning from the ztop of the substrate. we can add the 
    shift while scaning through the p0 parameter and ow 2D scan is also possible.

    We have comapared the Morse and Coulomb interaction for the PTCDA with the NaCl 8x8_L3 substrate, the results are the same with LAMMPS.

    Now we need to go in details with the GPU accelaration. 


    













 '''''''''



######################## This secetio is to see the grid slice along with the atomic positions #################################################################
# import sys
# import numpy as np
# import matplotlib.pyplot as plt

# sys.path.append("../../")
# from pyBall import atomicUtils as au
# from pyBall import MMFF as mmff


# mmff.init( xyz_name="data/xyz/PTCDA", surf_name="data/xyz/NaCl_8x8_L3_copy" ) 

# file_path = "/home/indranil/git/FireCore/cpp/common_resources/NaCl_8x8_L3_copy/Bspline_PLQd.npy"
# grid = np.load(file_path)

# print("Shape of grid:", grid.shape)  # expected (Nx, Ny, Nz, 3)

# channel = 0  # select the potential channel to view (e.g., Morse potential)

# # For this example we assume the grids share the same shape.
# Nx, Ny, Nz, _ = grid.shape

# # Use fixed indices (you can change these to whichever slice you want)
# mid_z = 199  # for horizontal (xy) slice at constant z
# mid_x = 308  # for vertical_x (yz) slice at constant x
# mid_y = 310  # for vertical_y (xz) slice at constant y

# # Extract slices from the grid
# slice_xy = grid[:, :, mid_z, channel]
# slice_yz = grid[mid_x, :, :, channel]
# slice_xz = grid[:, mid_y, :, channel]

# # Create a figure with 1 row and 3 columns
# fig, axs = plt.subplots(1, 3, figsize=(18, 6))

# # ---------------------------
# # Horizontal slice: xy plane.
# im0 = axs[0].imshow(slice_xy, origin="lower", cmap="viridis", aspect="auto")
# axs[0].set_title(f"Horizontal (xy) at z = {mid_z}")
# axs[0].set_xlabel("x")
# axs[0].set_ylabel("y")
# fig.colorbar(im0, ax=axs[0], fraction=0.046, pad=0.04)

# # Vertical slice: yz plane.
# im1 = axs[1].imshow(slice_yz.T, origin="lower", cmap="viridis", aspect="auto")
# axs[1].set_title(f"Vertical (yz) at x = {mid_x}")
# axs[1].set_xlabel("y")
# axs[1].set_ylabel("z")
# fig.colorbar(im1, ax=axs[1], fraction=0.046, pad=0.04)

# # Vertical slice: xz plane.
# im2 = axs[2].imshow(slice_xz.T, origin="lower", cmap="viridis", aspect="auto")
# axs[2].set_title(f"Vertical (xz) at y = {mid_y}")
# axs[2].set_xlabel("x")
# axs[2].set_ylabel("z")
# fig.colorbar(im2, ax=axs[2], fraction=0.046, pad=0.04)

# # ---------------------------
# # Call get_gridFF_info() BEFORE get_atom_positions()
# gff_shift0, gff_pos0, gff_natoms, gff_natoms_ = mmff.get_gridFF_info()

# # Get atom positions
# substrate_apos, molecule_apos = mmff.get_atom_positions()

# # --- Scaling and Shifting Atom Positions ---
# # Assuming grid origin is gff_pos0 and the grid spans from 0 to Nx, Ny, Nz
# # We need to scale and shift the atom positions to match the grid's coordinate system
# # This assumes that the grid's cell vectors are aligned with the Cartesian axes.

# print("gff_pos0:", gff_pos0)
# # Scale and shift the atom positions
# substrate_apos_scaled = substrate_apos - gff_pos0
# molecule_apos_scaled  = molecule_apos  - gff_pos0 - gff_shift0

# # --- Plotting Atom Positions ---
# # Horizontal (xy)
# axs[0].scatter(substrate_apos_scaled[:, 0], substrate_apos_scaled[:, 1], s=10, color='red', label='Substrate')
# axs[0].scatter(molecule_apos_scaled[:, 0], molecule_apos_scaled[:, 1], s=10, color='blue', label='Molecule')

# # Vertical (yz)
# axs[1].scatter(substrate_apos_scaled[:, 1], substrate_apos_scaled[:, 2], s=10, color='red')
# axs[1].scatter(molecule_apos_scaled[:, 1], molecule_apos_scaled[:, 2], s=10, color='blue')

# # Vertical (xz)
# axs[2].scatter(substrate_apos_scaled[:, 0], substrate_apos_scaled[:, 2], s=10, color='red')
# axs[2].scatter(molecule_apos_scaled[:, 0], molecule_apos_scaled[:, 2], s=10, color='blue')

# axs[0].legend()  # Add legend to the first subplot

# plt.tight_layout()
# plt.show()




############################## It is to see the scaled grid slice along with the atomic positions #################################################################
import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import MMFF as mmff


mmff.init( xyz_name="data/xyz/PTCDA", surf_name="data/xyz/NaCl_8x8_L3_copy" ) 


file_path = "/home/indranil/git/FireCore/cpp/common_resources/NaCl_8x8_L3_copy/Bspline_PLQd.npy"
grid = np.load(file_path)

print("Shape of grid:", grid.shape)  # expected (Nx, Ny, Nz, 3)

channel = 0  # select the potential channel to view (e.g., Morse potential)

# For this example we assume the grids share the same shape.
Nx, Ny, Nz, _ = grid.shape

# Use fixed indices (you can change these to whichever slice you want)
mid_z = 10 + Nz//2  # for horizontal (xy) slice at constant z
mid_x = 10 + Nx//2  # for vertical_x (yz) slice at constant x
mid_y = 10 + Ny//2  # for vertical_y (xz) slice at constant y

# Extract slices from the grid
slice_xy = grid[:, :, mid_z, channel]
slice_yz = grid[mid_x, :, :, channel]
slice_xz = grid[:, mid_y, :, channel]


# ---------------------------
# Retrieve grid parameters and atom positions from the C++ extension
gff_shift0, gff_pos0, gff_cell, gff_dCell, gff_natoms, gff_natoms_ = mmff.get_gridFF_info()
print("gff_pos0:", gff_pos0)
print("gff_shift0:", gff_shift0)
print("gff_cell:\n", gff_cell)
print("gff_dCell:\n", gff_dCell)

# For scaling purposes (we assume dCell is diagonal and defines the grid spacing)
dx = gff_dCell[0, 0]
dy = gff_dCell[1, 1]
dz = gff_dCell[2, 2]

# Compute physical extents for each slice.
# Here our coordinate system will be shifted so that the grid origin becomes zero.
extent_xy = [0, Nx * dx, 0, Ny * dy]     # x spans along axis0 and y axis1
extent_yz = [0, Ny * dy, 0, Nz * dz]       # for a slice taken at constant x; horizontal axis: y, vertical: z
extent_xz = [0, Nx * dx, 0, Nz * dz]       # for a slice taken at constant y; horizontal: x, vertical: z
# extent_xy = [ gff_pos0[0], gff_pos0[0] + Nx * dx, gff_pos0[1], gff_pos0[1] + Ny * dy ]
# extent_yz = [ gff_pos0[1], gff_pos0[1] + Ny * dy, gff_pos0[2], gff_pos0[2] + Nz * dz ]
# extent_xz = [ gff_pos0[0], gff_pos0[0] + Nx * dx, gff_pos0[2], gff_pos0[2] + Nz * dz ]

# ---------------------------
# Get atom positions (substate and molecule)
substrate_apos, molecule_apos = mmff.get_atom_positions()

# --- Scaling and Shifting Atom Positions ---
# The grid (and substrate) has a global origin gff_pos0.
# For plotting in a coordinate system where the lower left is 0, we subtract gff_pos0.
# Additionally, the molecule atoms are offset by gff_shift0.
substrate_apos_scaled = substrate_apos #+ gff_pos0
molecule_apos_scaled  = molecule_apos  - gff_shift0  #+ gff_pos0

# ---------------------------
# Plotting: Create a figure with three subplots for grid slices.
fig, axs = plt.subplots(1, 3, figsize=(18, 6))

# Horizontal slice (xy plane) at constant z
im0 = axs[0].imshow(slice_xy.T, extent=extent_xy, origin="lower", cmap="viridis", aspect="auto")
axs[0].set_title(f"Horizontal (xy) slice at z = {mid_z*dz}")
axs[0].set_xlabel("x (physical dimesion)")
axs[0].set_ylabel("y (physical dimension)")
fig.colorbar(im0, ax=axs[0], fraction=0.046, pad=0.04)
# Overlay atom positions (x-y)
axs[0].scatter(substrate_apos_scaled[:, 0], substrate_apos_scaled[:, 1],
               s=10, color='red', label='Substrate')
axs[0].scatter(molecule_apos_scaled[:, 0], molecule_apos_scaled[:, 1],
               s=10, color='blue', label='Molecule')
axs[0].legend()

# Vertical slice (yz plane) taken at constant x (mid_x)
im1 = axs[1].imshow(slice_yz.T, extent=extent_yz, origin="lower", cmap="viridis", aspect="auto")
axs[1].set_title(f"Vertical (yz) slice at x = {mid_x*dx}")
axs[1].set_xlabel("y (physical dimension)")
axs[1].set_ylabel("z (physical dimension)")
fig.colorbar(im1, ax=axs[1], fraction=0.046, pad=0.04)
# Overlay atom positions in y-z plane
axs[1].scatter(substrate_apos_scaled[:, 1], substrate_apos_scaled[:, 2],
               s=10, color='red')
axs[1].scatter(molecule_apos_scaled[:, 1], molecule_apos_scaled[:, 2],
               s=10, color='blue')

# Vertical slice (xz plane) taken at constant y (mid_y)
im2 = axs[2].imshow(slice_xz.T, extent=extent_xz, origin="lower", cmap="viridis", aspect="auto")
axs[2].set_title(f"Vertical (xz) slice at y = {mid_y*dy}")
axs[2].set_xlabel("x (physical dimesion)")
axs[2].set_ylabel("z (physical dimension)")
fig.colorbar(im2, ax=axs[2], fraction=0.046, pad=0.04)
# Overlay atom positions in x-z plane
axs[2].scatter(substrate_apos_scaled[:, 0], substrate_apos_scaled[:, 2],
               s=10, color='red')
axs[2].scatter(molecule_apos_scaled[:, 0], molecule_apos_scaled[:, 2],
               s=10, color='blue')

plt.tight_layout()
plt.show()