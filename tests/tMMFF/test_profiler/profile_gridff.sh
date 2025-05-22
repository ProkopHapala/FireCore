#!/bin/bash

# profile_gridff.sh - GridFF OpenCL profiler wrapper
# Usage: ./profile_gridff.sh [options]

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROFILER="${SCRIPT_DIR}/gridff_profiler.py"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
KERNEL_PATH="${PROJECT_ROOT}/cpp/common_resources/cl/GridFF.cl"

# Check if the profiler exists
if [ ! -f "$PROFILER" ]; then
    echo "Error: GridFF profiler not found at $PROFILER"
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
            echo "Usage: $0 [options]"
            echo ""
            echo "Options:"
            echo "  --job JOB           Job type (PLQ, Morse, etc.) [default: PLQ]"
            echo "  --dataset DATASET   Dataset name [default: NaCl_1x1_L1]"
            echo "  --output DIR        Output directory for profiling results [default: profile_results]"
            echo "  --help, -h          Show this help message"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# Print environment information
echo "===== Environment Information ====="
echo "Project Root: $PROJECT_ROOT"
echo "PYTHONPATH: $PYTHONPATH"
echo "Kernel Path: $KERNEL_PATH"
echo "Job: $JOB"
echo "Dataset: $DATASET"
echo "Output Directory: $OUTPUT_DIR"
echo "================================="

# Run the profiler
python "$PROFILER" --job "$JOB" --dataset "$DATASET" --output "$OUTPUT_DIR" --kernel-path "$KERNEL_PATH"
