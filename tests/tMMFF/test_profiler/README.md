# FireCore Plug-and-Play Profiler

This is a true plug-and-play profiler designed to work with any Python code, with special support for OpenCL applications. It follows the philosophy of standard profilers like NVIDIA's profiler, allowing you to profile existing code without any modifications.

## Features

- **Zero-modification profiling**: Works with your existing code without requiring any changes
- **Python function profiling**: Tracks execution time of Python functions using cProfile
- **OpenCL information**: Shows details about available OpenCL platforms and devices
- **Simple shell wrapper**: Easy-to-use shell script for profiling any Python script
- **Detailed reports**: Creates text files with comprehensive profiling data

## Usage

### Quick Start

The easiest way to profile any Python script is to use the `profile_any.sh` shell script:

```bash
# Profile any Python script
./profile_any.sh path/to/your_script.py [script_args...]

# Example: Profile the grid generation script
./profile_any.sh ../run_test_GridFF_ocl_new.py
```

This will run your script with profiling enabled and save the results to the `profile_results` directory.

### Manual Profiling

For more control over the profiling options, you can use the `run_profiler.py` script directly:

```bash
# Basic usage
python run_profiler.py -- python your_script.py [script_args...]

# Enable cProfile (more detailed profiling)
python run_profiler.py --cprofile -- python your_script.py [script_args...]

# Show OpenCL device information
python run_profiler.py --opencl-info -- python your_script.py [script_args...]

# Change the output directory
python run_profiler.py --output my_results -- python your_script.py [script_args...]

# Limit the number of functions shown in the output
python run_profiler.py --cprofile --limit 50 -- python your_script.py [script_args...]

# Change the sort order (cumulative, time, calls, pcalls, name)
python run_profiler.py --cprofile --sort time -- python your_script.py [script_args...]
```

## Output

The profiler generates output files in the `profile_results` directory (or a directory you specify with `--output`):

- `profile_YYYYMMDD_HHMMSS_profile.txt`: Standard cProfile output showing function call statistics

## How It Works

The profiler uses standard Python profiling tools to collect information without modifying your code:

1. **cProfile**: Standard Python profiling for detailed function-level profiling
2. **Time measurement**: Simple timing for overall execution time
3. **OpenCL information**: Displays information about available OpenCL platforms and devices
4. **Subprocess execution**: Runs your command as a subprocess to avoid interfering with its execution

This approach allows you to profile your existing code without modifying it, making it easy to identify bottlenecks and optimize performance.

## Examples

### Profiling Grid Generation

```bash
# Profile grid generation with default settings
./profile_any.sh ../run_test_GridFF_ocl_new.py

# Profile with cProfile and save results to a specific directory
python run_profiler.py --cprofile --output grid_profiles -- python ../run_test_GridFF_ocl_new.py
```

### Analyzing Results

The profiling results show the time spent in each function, allowing you to identify bottlenecks:

1. Look for functions that take a long time to execute
2. Identify functions that are called frequently
3. Focus optimization efforts on the most time-consuming parts of your code

For example, in the grid generation code, you might find that the `fit3D` method takes a significant amount of time, making it a good candidate for optimization.
