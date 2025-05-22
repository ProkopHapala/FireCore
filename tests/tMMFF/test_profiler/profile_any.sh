#!/bin/bash

# Simple shell script to profile any Python script without modification
# Usage: ./profile_any.sh [script_path] [script_args...]

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROFILER="$SCRIPT_DIR/run_profiler.py"
OUTPUT_DIR="$SCRIPT_DIR/profile_results"

# Check if script path is provided
if [ $# -lt 1 ]; then
    echo "Usage: $0 [script_path] [script_args...]"
    echo "Example: $0 ../run_test_GridFF_ocl_new.py"
    exit 1
fi

# Get script path and arguments
SCRIPT_PATH="$1"
shift
SCRIPT_ARGS="$@"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Run the profiler
python "$PROFILER" --cprofile --opencl-info --output "$OUTPUT_DIR" -- python "$SCRIPT_PATH" $SCRIPT_ARGS
