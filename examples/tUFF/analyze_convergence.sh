#!/bin/bash
# Script to extract and plot convergence data from H2O UFF test

LOGFILE="log_h2o_100k_steps.txt"
DATAFILE="convergence_data.txt"

echo "Extracting convergence data from $LOGFILE..."

# Extract step number and force magnitude
grep "DEBUG: isys=0 nbEval=" $LOGFILE | \
    awk '{print $3, $4}' | \
    sed 's/nbEval=//' | \
    sed 's/|F|=//' > $DATAFILE

echo "Data extracted to $DATAFILE"
echo "Total steps: $(wc -l < $DATAFILE)"

# Calculate statistics
echo ""
echo "=== Convergence Statistics ==="
awk '{
    sum+=$2; 
    sumsq+=$2*$2;
    if($2<min || NR==1) min=$2; 
    if($2>max || NR==1) max=$2;
    if($2<1e-5) count_1e5++;
    if($2<1e-6) count_1e6++;
} 
END {
    avg=sum/NR; 
    stddev=sqrt(sumsq/NR - avg*avg);
    print "Total steps:", NR;
    print "Min |F|:", min, "eV/Å";
    print "Max |F|:", max, "eV/Å";
    print "Avg |F|:", avg, "eV/Å";
    print "Std Dev:", stddev, "eV/Å";
    print "Steps below 1e-5:", count_1e5, "("100*count_1e5/NR"%)";
    print "Steps below 1e-6:", count_1e6, "("100*count_1e6/NR"%)";
}' $DATAFILE

# Find best convergence
echo ""
echo "=== Best Convergence ==="
sort -k2 -g $DATAFILE | head -10 | awk '{printf "Step %d: |F| = %.3e eV/Å\n", $1, $2}'

# Check last 1000 steps
echo ""
echo "=== Last 1000 Steps Statistics ==="
tail -1000 $DATAFILE | awk '{
    sum+=$2; 
    if($2<min || NR==1) min=$2; 
    if($2>max || NR==1) max=$2;
} 
END {
    print "Min:", min, "eV/Å";
    print "Max:", max, "eV/Å"; 
    print "Avg:", sum/NR, "eV/Å";
}'

echo ""
echo "To plot the data, run: python3 plot_convergence.py"
