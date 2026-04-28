# README: UFF Throughput Simulations

This document describes the procedure for running and analyzing UFF (Universal Force Field) throughput simulations. The goal is to efficiently explore the conformational space of a molecule on various surfaces and to analyze the results.

## Workflow

The workflow consists of three main steps:
1.  **Prepare Surface Data**: Generate B-spline data (GridFF) for the surfaces.
2.  **Run Simulations**: Use the `run_throughput_UFF.sh` script to launch a series of simulations.
3.  **Analyze Results**: Process the output data and generate plots.

---

### Step 1: Prepare Surface Data (GridFF Generation)

To efficiently calculate the molecule-surface interactions, it is necessary to pre-generate data for B-spline interpolation (GridFF).

**Procedure:**

1.  **Copy Surfaces**: Copy the `.xyz` surface files from the `data/xyz/surfaces_for_throughput` directory to a temporary folder outside this repository.

2.  **Switch Branch**: Switch to the branch containing the grid generation tools:
    ```bash
    git switch debug/grid_generation_and_test
    ```

3.  **Generate GridFF Data**: Run the script to generate the grid for all surfaces.
    * Copy the folder `data/xyz/surfaces_for_throughput` to `data/xyz/` in the `debug/grid_generation_and_test` branch.
    ```bash
    python3 ../tMMFF/run_test_GridFF_ocl_new.py
    ```
    This script should process all `.xyz` files in the given folder and generate the corresponding GridFF files.

4.  **Save and Return**: Copy the generated files to a safe location and switch back to the main branch:
    ```bash
    git switch prokop_and_master
    ```

5.  **Copy Data**: Copy the generated GridFF files into the `data/` directory of this project.

---

### Step 2: Run Simulations (`run_throughput_UFF.sh`)

The main script `run_throughput_UFF.sh` orchestrates the execution of individual simulations with various parameters.

**Setting Parameters in the Script:**

Before running, you can modify the following parameters directly in the `run_throughput_UFF.sh` script:

*   **Interactions**:
    *   `dovdW`: Controls Non-bonded interactions. This should be enabled all the time.
    *   `doSurfAtoms`: A positive value (e.g., `1`) enables interaction with the surface atoms (explicit calculation). A negative value disables the surface.
    *   `bGridFF`: A positive value (use `6` for B-spline) activates the surface interaction calculation using the pre-generated GridFF grid.

*   **Surface Selection**:
    *   The `Ns` array allows you to specify which surface sizes to use. The script iterates through the values in this array and looks for a surface named `NaCl_NxN_L3` for each size `N`. Use only one size eG `Ns=(3)` for debugging.
    ```bash
    Ns=(3 4 5 6 7 8 9 10 11 12 13 14 15 16) # Example: uses surfaces from 3x3 to 16x16
    ```
    You can also use a different surface by modifying the `surf_name` variable. Use `Cl_hole` or `Na_hole` instead of `L3` to use surfaces with a hole. eG `surf_name="common_resources/xyz/surfaces_for_throughput/NaCl_${N}x${N}_Cl_hole"`

*   **Simulation Parameters**:
    *   `replicas`: Number of parallel replicas (systems).
    *   `perframes`: Number of MD steps in one "frame".
    *   `perVF`: Number of force field evaluations within a single step.
    To generate data for the article, use the following parameters:
    ```bash
    replicas=(5000)
    perframes=(100 20)
    perVF=(100 20)
    ```
    Previously we used two lines (perframes=100, perVF=100 and perframes=20, perVF=20)

**Execution:**

The script is run directly from the command line:
```bash
./run_throughput_UFF.sh
```

The script automatically recompiles the necessary C++/OCL library and then runs `run_throughput_UFF.py` with various parameter combinations.

---

### Important Code Modifications for Data Generation

Depending on the data you want to generate, you may need to make small modifications to the C++ code in `cpp/common/molecular/MolWorld_sp3_multi.h` before running the simulations.

#### For the `eval_vs_time` plot:

To generate data that tracks the discovery of unique minima over time, you must enable saving each new structure to the database.

1.  **File**: `cpp/common/molecular/MolWorld_sp3_multi.h`
2.  **Function**: `double evalVFs( double Fconv=1e-6 )`
3.  **Change**: Locate the following `if` statement:

    ```cpp
    // save
    if(0*bSaveToDatabase){  // use it for evaluation_vs_time, do not use for nb_evale_vs_surf_size
        int sameMember = database->addIfNewDescriptor(&ffu);
        if(sameMember==-1){
            // ... code to save structure
        }
    }
    ```

4.  **Action**: Change `0*bSaveToDatabase` to `1*bSaveToDatabase` (or simply `bSaveToDatabase`) to enable this block. This will save every newly found converged structure, which is required for this plot.

5.  **Simulation Time**: The total duration of the simulation is controlled internally within the `MDLoop` function in `cpp/common/molecular/MolWorld_sp3_multi.h`. There is time limit of 59.5s.

#### For the `eval_vs_surf_size` plot:

For this plot, you are interested in the final statistics of a simulation of a fixed duration, not the history of found minima.

1.  **Database Saving**: Ensure the code block mentioned above is **disabled** (i.e., it remains `if(0*bSaveToDatabase)`). This prevents writing a large number of intermediate structures.

2.  **Simulation Time**: The total duration of the simulation is controlled internally within the `MDLoop` function in `cpp/common/molecular/MolWorld_sp3_multi.h`. There is time limit of 9.5s.

--- 

### Step 3: Analyze Results

The script generates several types of output files.

**Output Files:**

*   `traj_UFF_000.xyz`: Trajectory of the zeroth replica. Useful for visual inspection.
*   `results/all/`: A directory containing detailed data for each individual simulation. The name of each `minima*.dat` file is encoded according to the parameters used.
*   `results/all/results.dat`: A summary file where the last line from each `minima*.dat` file is stored. This provides overview of the results of all simulations.
*   `minima*.dat`: These files (renamed and stored in `results/all/`) contain the history of minima found during a single simulation.

**Important Columns in `minima*.dat` (and `results.dat`):**

*   **`converged structures`**: The number of structures that reached convergence (satisfied the force criterion). May be listed duplicately, calculated in different ways.
*   **`uniquely converged structures`**: The number of unique converged structures.
*   **`time`**: The total simulation time.
*   **`total evaluations`**: The total number of MD steps performed.
*   **`avg steps to converge`**: The average number of steps required to converge a single structure.

**Generating Plots:**

To visualize the results, two main plots need to be generated:

1.  **`eval_vs_time`**:
    *   **Purpose**: Shows how the number of unique minima found changes over time.
    *   **Data**: Data from the `minima*.dat` files (stored in `results/all/`) are used for this plot. You need to extract the time and the number of unique minima from the individual lines of these files.

2.  **`eval_vs_surf_size`**:
    *   **Purpose**: Shows how the total number of force evaluations depends on the surface size.
    *   **Data**: Data from the summary file `results/all/results.dat` are used for this plot. You need to plot the dependency of the total number of evaluations (one of the columns) on the surface size `N` (obtained from the filename in the first column).