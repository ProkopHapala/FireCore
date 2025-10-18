# Plotting Script: Evaluation Count vs. Time

This directory contains a Python script `evaluation_vs_time.py` designed to generate plots that compare the performance of different molecular mechanics force field methods. It visualizes the number of evaluations and discovered structures against computation time.

## Usage Workflow

To generate the plots, follow these steps:

1.  **Prepare Data Files:**
    *   First, you need to obtain the results data. These are  located in a `tUFF/results/all` directory from your run.
    *   From that directory, copy two specific `.dat` summary files into this one (`tests/tUFF/plot_dependance_nb_eval_to_time/`).
    *   The required files are:
        *   One file for the "all atoms" method, which has the code `11000` in its filename (e.g., `minima__11000_...dat`). This is labeled as 'No GridFF' in the plot legend.
        *   One file for the "GridFF" method, which has the code `11110` in its filename (e.g., `minima__11110_...dat`). This is labeled as 'GridFF' in the plot legend.

2.  **Run the Script:**
    *   Navigate to this directory in your terminal.
    *   Execute the Python script:
        ```bash
        python3 evaluation_vs_time.py
        ```

3.  **Check the Output:**
    *   These plots will be saved in the current directory in both `.png` and `.pdf` formats (e.g., `evaluations_linX_linY.png`).

## Plotted Data Details

The script generates two main types of plots from the data files, with the data sourced from specific columns in the `.dat` files:

1.  **Evaluations vs. Time:** This plot shows how the total number of force field evaluations grows over time.
    *   **X-axis:** `Time` (Time in seconds)
    *   **Y-axis:** `nbEvaluation` (Total number of evaluations)

2.  **Structures vs. Time:** These plots show the number of discovered molecular structures over time. The script generates multiple versions of this plot with different axis scales (linear-linear, linear-log, log-linear, and log-log).
    *   **X-axis:** `Time` (Time in seconds)
    *   **Y-axis (solid line):** `Totally (conv+nonconv)` (Total number of all structures found, both converged and non-converged)
    *   **Y-axis (dashed line):** `Unique conv` (Number of unique converged structures)

For the article we used `structures_linX_linY`
