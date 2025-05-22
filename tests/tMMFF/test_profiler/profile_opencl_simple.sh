#!/bin/bash

# profile_opencl_simple.sh - Simple OpenCL profiler wrapper
# Usage: ./profile_opencl_simple.sh your_script.py [args...]

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROFILER="${SCRIPT_DIR}/opencl_profiler.py"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"

# Check if the profiler exists
if [ ! -f "$PROFILER" ]; then
    echo "Error: OpenCL profiler not found at $PROFILER"
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

# Run the profiler with the given command
python "$PROFILER" python "$@"
