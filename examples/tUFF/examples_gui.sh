#!/bin/bash

# Examples of running the UFF scanning GUI with different configurations

echo "=== UFF Scanning GUI Examples ==="
echo ""
echo "Choose an example to run:"
echo ""
echo "1) Basic - NaCl 6x6 surface with xylitol (bug jump vhen enough large surface scannerd to atom without charge)"
echo "2) Small surface - NaCl 3x3"
echo "3) Large surface - NaCl 10x10"
echo "4) Cl hole surface - NaCl 6x6 with Cl vacancy"
echo "5) Na hole surface - NaCl 6x6 with Na vacancy"
echo "6) Custom - specify your own parameters"
echo ""
read -p "Enter choice [1-6]: " choice

case $choice in
    1)
        echo "Running basic example with NaCl 6x6..."
        ./run_scan_gui.sh NaCl_6x6_L3
        ;;
    2)
        echo "Running with small NaCl 3x3 surface..."
        ./run_scan_gui.sh NaCl_3x3_L3
        ;;
    3)
        echo "Running with large NaCl 10x10 surface..."
        ./run_scan_gui.sh NaCl_10x10_L3
        ;;
    4)
        echo "Running with Cl hole surface..."
        ./run_scan_gui.sh NaCl_6x6_Cl_hole
        ;;
    5)
        echo "Running with Na hole surface..."
        ./run_scan_gui.sh NaCl_6x6_Na_hole
        ;;
    6)
        echo ""
        read -p "Enter surface name (e.g., NaCl_6x6_L3): " surf_name
        read -p "Enter molecule file (default: xylitol_WO_gridFF.xyz): " mol_file
        read -p "Enter initial scan atom index (default: 5): " scan_atom
        
        mol_file=${mol_file:-data/xyz/xylitol_WO_gridFF.xyz}
        scan_atom=${scan_atom:-5}
        
        echo "Running with custom parameters..."
        python3 -u scan_files/scan_gui.py \
            --surf "xyz/surfaces_for_throughput/$surf_name" \
            --molecule "$mol_file" \
            --scan_atom "$scan_atom" \
            --preset grid-only \
            --grid-ff
        ;;
    *)
        echo "Invalid choice. Exiting."
        exit 1
        ;;
esac

