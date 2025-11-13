#!/bin/bash
# Systematic dt sweep for H2O convergence testing
# Tests different time steps to find optimal value

MOLECULE="H2O"
STEPS=10000

# Array of dt values to test
DT_VALUES=(0.1 0.05 0.02 0.01 0.005 0.002 0.001)

echo "=========================================="
echo "  Time Step (dt) Sweep Test"
echo "=========================================="
echo "Molecule: ${MOLECULE}"
echo "Steps per test: ${STEPS}"
echo "dt values: ${DT_VALUES[@]}"
echo "=========================================="
echo ""

for dt in "${DT_VALUES[@]}"; do
    echo ""
    echo "=========================================="
    echo "  Testing dt = ${dt}"
    echo "=========================================="
    
    LOGFILE="log_${MOLECULE}_dt${dt}.txt"
    
    # Run simulation
    python3 run_throughput_UFF.py \
        --xyz_name data/xyz/${MOLECULE}.xyz \
        --nSys 1 \
        --bUFF 1 \
        --bGridFF 1 \
        --gridnPBC "(1,1,0)" \
        --loops 1 \
        --perframe ${STEPS} \
        --perVF 100 \
        --Fconv 1e-6 \
        --dt ${dt} \
        2>&1 | tee ${LOGFILE}
    
    # Analyze results
    echo ""
    echo "Analyzing dt=${dt}..."
    python3 analyze_and_plot.py ${LOGFILE} "${MOLECULE}_dt${dt}"
    
    echo ""
    echo "✓ Completed dt=${dt}"
    echo ""
done

echo ""
echo "=========================================="
echo "  All tests complete!"
echo "=========================================="
echo ""
echo "Generating comparison summary..."
python3 compare_dt_results.py ${MOLECULE}
