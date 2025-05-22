#!/bin/bash

# profile_opencl.sh - Specialized profiler for OpenCL applications
# Usage: ./profile_opencl.sh [options] script.py [script_args...]

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DETAILED_PROFILER="${SCRIPT_DIR}/detailed_profiler.py"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"

# Check if the detailed profiler exists
if [ ! -f "$DETAILED_PROFILER" ]; then
    echo "Error: Detailed profiler not found at $DETAILED_PROFILER"
    exit 1
fi

# Parse options
PROFILER_ARGS=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --no-opencl)
            PROFILER_ARGS="$PROFILER_ARGS --no-opencl"
            shift
            ;;
        --no-trace)
            PROFILER_ARGS="$PROFILER_ARGS --no-trace"
            shift
            ;;
        --output)
            PROFILER_ARGS="$PROFILER_ARGS --output $2"
            shift 2
            ;;
        --limit)
            PROFILER_ARGS="$PROFILER_ARGS --limit $2"
            shift 2
            ;;
        --sort)
            PROFILER_ARGS="$PROFILER_ARGS --sort $2"
            shift 2
            ;;
        --help|-h)
            echo "Usage: $0 [options] script.py [script_args...]"
            echo ""
            echo "Options:"
            echo "  --no-opencl       Disable OpenCL profiling"
            echo "  --no-trace        Disable function tracing"
            echo "  --output DIR      Output directory for profiling results"
            echo "  --limit N         Limit the number of functions shown in output"
            echo "  --sort KEY        Sort order (time, cumulative, calls, name)"
            echo "  --help, -h        Show this help message"
            exit 0
            ;;
        *)
            # First non-option argument is the script to profile
            break
            ;;
    esac
done

# Check if a script was provided
if [ $# -eq 0 ]; then
    echo "Error: No script specified"
    echo "Usage: $0 [options] script.py [script_args...]"
    exit 1
fi

# Set OpenCL environment variables for better profiling
export PYOPENCL_COMPILER_OUTPUT=1
export OPENCL_PROFILE=1

# For NVIDIA
export CUDA_PROFILE=1
export CUDA_PROFILE_CSV=1

# For AMD
export AMD_OCL_BUILD_OPTIONS_APPEND='-cl-opt-disable'

# For Intel
export INTEL_OPENCL_PROFILE=1

# Add FireCore project root to PYTHONPATH
export PYTHONPATH="$PROJECT_ROOT:$PYTHONPATH"

# Print environment information
echo "===== Environment Information ====="
echo "Project Root: $PROJECT_ROOT"
echo "PYTHONPATH: $PYTHONPATH"
echo "OpenCL Profiling: Enabled"
echo "================================="

# Run the detailed profiler
python "$DETAILED_PROFILER" $PROFILER_ARGS "$@"
