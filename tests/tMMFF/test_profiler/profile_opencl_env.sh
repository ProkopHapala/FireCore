#!/bin/bash

# profile_opencl_env.sh - OpenCL Environment Profiler wrapper
# Usage: ./profile_opencl_env.sh [options] [command]

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROFILER="${SCRIPT_DIR}/opencl_env_profiler.py"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"

# Check if the profiler exists
if [ ! -f "$PROFILER" ]; then
    echo "Error: OpenCL Environment profiler not found at $PROFILER"
    exit 1
fi

# Make sure the profiler is executable
chmod +x "$PROFILER"

# Add FireCore project root to PYTHONPATH
export PYTHONPATH="$PROJECT_ROOT:$PYTHONPATH"

# Parse command-line arguments
JOB="PLQ"
DATASET="NaCl_1x1_L1"
OUTPUT_DIR="profile_results"
COMMAND_ARGS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --job)
            JOB="$2"
            shift 2
            ;;
        --dataset)
            DATASET="$2"
            shift 2
            ;;
        --output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --help|-h)
            echo "Usage: $0 [options] [command]" 
            echo ""
            echo "Options:"
            echo "  --job JOB           Job type (PLQ, Morse, etc.) [default: PLQ]"
            echo "  --dataset DATASET   Dataset name [default: NaCl_1x1_L1]"
            echo "  --output DIR        Output directory for profiling results [default: profile_results]"
            echo "  --help, -h          Show this help message"
            echo ""
            echo "If a command is provided, it will be profiled. Otherwise, GridFF will be profiled."
            exit 0
            ;;
        *)
            COMMAND_ARGS+=("$1")
            shift
            ;;
    esac
done

# Print environment information
echo "===== Environment Information ====="
echo "Project Root: $PROJECT_ROOT"
echo "PYTHONPATH: $PYTHONPATH"
echo "Job: $JOB"
echo "Dataset: $DATASET"
echo "Output Directory: $OUTPUT_DIR"
echo "================================="

# Run the profiler
if [ ${#COMMAND_ARGS[@]} -eq 0 ]; then
    # Run GridFF profiling
    python "$PROFILER" --job "$JOB" --dataset "$DATASET" --output "$OUTPUT_DIR"
else
    # Run command profiling
    python "$PROFILER" --output "$OUTPUT_DIR" "${COMMAND_ARGS[@]}"
fi
