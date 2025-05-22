#!/bin/bash

# profile_detailed.sh - Wrapper script for the detailed profiler
# Usage: ./profile_detailed.sh [options] script.py [script_args...]

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DETAILED_PROFILER="${SCRIPT_DIR}/detailed_profiler.py"

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
        --help|-h)
            echo "Usage: $0 [options] script.py [script_args...]"
            echo ""
            echo "Options:"
            echo "  --no-opencl       Disable OpenCL profiling"
            echo "  --no-trace        Disable function tracing"
            echo "  --output DIR      Output directory for profiling results"
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

# Run the detailed profiler
python "$DETAILED_PROFILER" $PROFILER_ARGS "$@"
