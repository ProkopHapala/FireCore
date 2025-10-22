#!/usr/bin/env python3
"""
Automatic convergence analysis and plotting for UFF double precision tests
Usage: python3 analyze_and_plot.py <logfile> <molecule_name>
"""
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import re

def extract_convergence_data(logfile):
    """Extract convergence data from log file"""
    steps = []
    forces = []
    
    with open(logfile, 'r') as f:
        for line in f:
            if 'DEBUG: isys=0 nbEval=' in line:
                # Extract step and force
                match = re.search(r'nbEval=(\d+).*\|F\|=([\d.e+-]+)', line)
                if match:
                    steps.append(int(match.group(1)))
                    forces.append(float(match.group(2)))
    
    return np.array(steps), np.array(forces)

def print_statistics(steps, forces, molecule):
    """Print convergence statistics"""
    print(f"\n{'='*60}")
    print(f"  Convergence Analysis: {molecule}")
    print(f"{'='*60}\n")
    
    print(f"Total steps: {len(steps)}")
    print(f"Min |F|: {forces.min():.3e} eV/Å (step {steps[forces.argmin()]:.0f})")
    print(f"Max |F|: {forces.max():.3e} eV/Å")
    print(f"Final |F|: {forces[-1]:.3e} eV/Å")
    print(f"Mean |F|: {forces.mean():.3e} eV/Å")
    print(f"Median |F|: {np.median(forces):.3e} eV/Å")
    
    # Convergence milestones
    print(f"\n{'='*60}")
    print("  Convergence Milestones")
    print(f"{'='*60}\n")
    milestones = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6]
    for threshold in milestones:
        idx = np.where(forces < threshold)[0]
        if len(idx) > 0:
            first_step = steps[idx[0]]
            count = len(idx)
            percent = 100 * count / len(steps)
            status = "✓" if threshold <= 1e-4 else "✗"
            print(f"{status} |F| < {threshold:.0e}: First at step {first_step:6.0f}, achieved {count:5d} times ({percent:5.1f}%)")
        else:
            print(f"✗ |F| < {threshold:.0e}: Never achieved")
    
    # Last 10% statistics
    last_10pct = int(len(forces) * 0.1)
    last_forces = forces[-last_10pct:]
    print(f"\n{'='*60}")
    print(f"  Last {last_10pct} Steps Statistics")
    print(f"{'='*60}\n")
    print(f"Min: {last_forces.min():.3e} eV/Å")
    print(f"Max: {last_forces.max():.3e} eV/Å")
    print(f"Mean: {last_forces.mean():.3e} eV/Å")
    print(f"Std Dev: {last_forces.std():.3e} eV/Å")
    if last_forces.min() > 0:
        print(f"Oscillation range: {last_forces.max()/last_forces.min():.1f}x")

def create_plots(steps, forces, molecule, output_file):
    """Create convergence plots"""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # 1. Full convergence plot (log scale)
    ax = axes[0, 0]
    ax.semilogy(steps, forces, 'b-', linewidth=0.5, alpha=0.7)
    ax.axhline(y=1e-4, color='r', linestyle='--', linewidth=2, label='1e-4 (float limit)')
    ax.axhline(y=1e-6, color='g', linestyle='--', linewidth=2, label='1e-6 (target)')
    ax.set_xlabel('MD Steps', fontsize=12)
    ax.set_ylabel('|F| (eV/Å)', fontsize=12)
    ax.set_title(f'{molecule}: Full Convergence History (Log Scale)', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=10)
    
    # 2. Initial convergence (first 20%)
    ax = axes[0, 1]
    cutoff = int(len(steps) * 0.2)
    ax.plot(steps[:cutoff], forces[:cutoff], 'b-', linewidth=0.8)
    ax.axhline(y=1e-4, color='r', linestyle='--', linewidth=2, label='1e-4')
    ax.set_xlabel('MD Steps', fontsize=12)
    ax.set_ylabel('|F| (eV/Å)', fontsize=12)
    ax.set_title(f'Initial Convergence (First {cutoff} steps)', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=10)
    
    # 3. Final convergence (last 30%)
    ax = axes[1, 0]
    cutoff = int(len(steps) * 0.3)
    ax.semilogy(steps[-cutoff:], forces[-cutoff:], 'b-', linewidth=0.8)
    ax.axhline(y=1e-5, color='orange', linestyle='--', linewidth=2, label='1e-5')
    ax.axhline(y=1e-6, color='g', linestyle='--', linewidth=2, label='1e-6 (target)')
    ax.set_xlabel('MD Steps', fontsize=12)
    ax.set_ylabel('|F| (eV/Å)', fontsize=12)
    ax.set_title(f'Final Convergence (Last {cutoff} steps)', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=10)
    
    # 4. Histogram of forces
    ax = axes[1, 1]
    ax.hist(np.log10(forces), bins=100, edgecolor='black', alpha=0.7, color='steelblue')
    ax.axvline(x=np.log10(1e-4), color='r', linestyle='--', linewidth=2, label='1e-4')
    ax.axvline(x=np.log10(1e-6), color='g', linestyle='--', linewidth=2, label='1e-6')
    ax.set_xlabel('log10(|F|) [eV/Å]', fontsize=12)
    ax.set_ylabel('Count', fontsize=12)
    ax.set_title('Distribution of Forces', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='y')
    ax.legend(fontsize=10)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"\n{'='*60}")
    print(f"  Plot saved as: {output_file}")
    print(f"{'='*60}\n")

def main():
    if len(sys.argv) < 3:
        print("Usage: python3 analyze_and_plot.py <logfile> <molecule_name>")
        sys.exit(1)
    
    logfile = sys.argv[1]
    molecule = sys.argv[2]
    
    # Extract data
    print(f"\nExtracting convergence data from {logfile}...")
    steps, forces = extract_convergence_data(logfile)
    
    if len(steps) == 0:
        print("ERROR: No convergence data found in log file!")
        sys.exit(1)
    
    # Save raw data
    datafile = f"convergence_{molecule}.txt"
    np.savetxt(datafile, np.column_stack([steps, forces]), 
               fmt=['%d', '%.6e'], 
               header='Step Force(eV/A)')
    print(f"Data saved to: {datafile}")
    
    # Print statistics
    print_statistics(steps, forces, molecule)
    
    # Create plots
    plotfile = f"convergence_{molecule}.png"
    create_plots(steps, forces, molecule, plotfile)
    
    print(f"\n{'='*60}")
    print(f"  Analysis Complete!")
    print(f"{'='*60}")
    print(f"  Data file: {datafile}")
    print(f"  Plot file: {plotfile}")
    print(f"{'='*60}\n")

if __name__ == '__main__':
    main()
