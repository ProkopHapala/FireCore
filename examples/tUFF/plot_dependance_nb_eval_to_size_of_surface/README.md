# Performance vs. Surface Size Analysis

This directory contains the Python script `eval_vs_surf_size.py`, which is used to analyze and plot the performance of different force field methods against the size of the surface being simulated.

## How to Use

Follow these steps to generate the analysis plots and data files.

### 1. Prepare the Input Data

- The script requires an input file named `uff_data.dat` to be present in this directory.
- This file should be populated with data from the main results file, located at `tUFF/results/all/results.dat`.
- Copy the relevant result lines from `results.dat` into `uff_data.dat`.
- **Important**: For each method `code` you want to analyze, you should include at least two corresponding result lines in `uff_data.dat`. If you are updating the data, it is best to replace old lines with new ones for the same code to ensure the analysis is current.

### 2. Run the Analysis Script

- Once the `uff_data.dat` file is prepared, run the script from within this directory:
    ```bash
    python3 eval_vs_surf_size.py
    ```

## Generated Output

The script will generate the following files in the current directory:

- **Performance Plots (`images/performance_*.png`):**
    Several `.png` image files are generated, showing the relationship between the number of surface atoms and the performance (either in "Evaluations per second" or "Nanoseconds per evaluation"). The plots are generated with various linear and logarithmic axis scales.

- **Performance Data (`performance_vs_number_atoms.dat`):**
    A tab-separated data file containing a summary of the best and worse performance metrics for each algorithm (`No Surface`, `No GridFF`, `GridFF`) at different surface sizes.

### Note on Article Data

The data generated in `performance_vs_number_atoms.dat` is significant as it was the basis for the article. The results from this file were used by Indranil for further analysis with his script.
