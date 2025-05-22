#!/bin/bash

# profile_kernels.sh - OpenCL kernel profiler wrapper
# Usage: ./profile_kernels.sh your_script.py [args...]

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROFILER="${SCRIPT_DIR}/opencl_kernel_profiler.py"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"

# Check if the profiler exists
if [ ! -f "$PROFILER" ]; then
    echo "Error: OpenCL kernel profiler not found at $PROFILER"
    exit 1
fi

# Make sure the profiler is executable
chmod +x "$PROFILER"

# Add FireCore project root to PYTHONPATH
export PYTHONPATH="$PROJECT_ROOT:$PYTHONPATH"

# Print environment information
echo "===== Environment Information ====="
echo "Project Root: $PROJECT_ROOT"
echo "PYTHONPATH: $PYTHONPATH"
echo "================================="

# Run the profiler with the given script
python "$PROFILER" "$@"
