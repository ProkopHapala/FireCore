#!/usr/bin/env python3
"""
Compare convergence results across different dt values
Creates comprehensive comparison plots and summary table
"""
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import glob
import re

def load_convergence_data(filename):
    """Load convergence data from file"""
    try:
        data = np.loadtxt(filename)
        return data[:, 0], data[:, 1]  # steps, forces
    except:
        return None, None

def extract_dt_from_filename(filename):
    """Extract dt value from filename"""
    match = re.search(r'dt([\d.]+)', filename)
    if match:
        dt_str = match.group(1).rstrip('.')  # Remove trailing dot
        return float(dt_str)
    return None

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 compare_dt_results.py <molecule_name>")
        sys.exit(1)
    
    molecule = sys.argv[1]
    
    # Find all convergence data files for this molecule
    pattern = f"convergence_{molecule}_dt*.txt"
    files = sorted(glob.glob(pattern))
    
    if len(files) == 0:
        print(f"No convergence data files found matching: {pattern}")
        sys.exit(1)
    
    print(f"\nFound {len(files)} dt test results:")
    
    # Collect data
    results = []
    for f in files:
        dt = extract_dt_from_filename(f)
        steps, forces = load_convergence_data(f)
        if steps is not None and dt is not None:
            results.append({
                'dt': dt,
                'steps': steps,
                'forces': forces,
                'min_force': forces.min(),
                'final_force': forces[-1],
                'mean_force': forces.mean(),
                'median_force': np.median(forces),
                'below_1e4': np.sum(forces < 1e-4),
                'below_1e5': np.sum(forces < 1e-5),
                'below_1e6': np.sum(forces < 1e-6),
            })
            print(f"  dt={dt}: min |F|={forces.min():.3e}, final |F|={forces[-1]:.3e}")
    
    # Sort by dt
    results = sorted(results, key=lambda x: x['dt'], reverse=True)
    
    # Print summary table
    print(f"\n{'='*100}")
    print(f"  Time Step Comparison Summary - {molecule}")
    print(f"{'='*100}")
    print(f"{'dt':>8} | {'Min |F|':>12} | {'Final |F|':>12} | {'Mean |F|':>12} | {'<1e-4':>6} | {'<1e-5':>6} | {'<1e-6':>6} | Status")
    print(f"{'-'*100}")
    
    best_dt = None
    best_final = float('inf')
    
    for r in results:
        status = ""
        if r['final_force'] < 1e-4:
            status = "✓ Excellent"
            if r['final_force'] < best_final:
                best_final = r['final_force']
                best_dt = r['dt']
        elif r['final_force'] < 1e-3:
            status = "✓ Good"
        elif r['final_force'] < 1e-2:
            status = "○ Fair"
        else:
            status = "✗ Poor"
        
        pct_1e4 = 100 * r['below_1e4'] / len(r['forces'])
        pct_1e5 = 100 * r['below_1e5'] / len(r['forces'])
        pct_1e6 = 100 * r['below_1e6'] / len(r['forces'])
        
        print(f"{r['dt']:8.4f} | {r['min_force']:12.3e} | {r['final_force']:12.3e} | {r['mean_force']:12.3e} | "
              f"{pct_1e4:5.1f}% | {pct_1e5:5.1f}% | {pct_1e6:5.1f}% | {status}")
    
    print(f"{'='*100}")
    if best_dt:
        print(f"  ★ BEST: dt={best_dt} achieved final |F|={best_final:.3e} eV/Å")
    print(f"{'='*100}\n")
    
    # Create comparison plots
    fig = plt.figure(figsize=(16, 12))
    
    # 1. All convergence curves (log scale)
    ax1 = plt.subplot(2, 3, 1)
    colors = plt.cm.viridis(np.linspace(0, 1, len(results)))
    for i, r in enumerate(results):
        ax1.semilogy(r['steps'], r['forces'], '-', linewidth=1.5, 
                     label=f"dt={r['dt']}", color=colors[i], alpha=0.8)
    ax1.axhline(y=1e-4, color='r', linestyle='--', linewidth=2, label='1e-4 target')
    ax1.axhline(y=1e-6, color='g', linestyle='--', linewidth=2, label='1e-6 target')
    ax1.set_xlabel('MD Steps', fontsize=11)
    ax1.set_ylabel('|F| (eV/Å)', fontsize=11)
    ax1.set_title('All Convergence Curves (Log Scale)', fontsize=12, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=9, loc='best')
    
    # 2. Final force vs dt
    ax2 = plt.subplot(2, 3, 2)
    dts = [r['dt'] for r in results]
    finals = [r['final_force'] for r in results]
    ax2.semilogy(dts, finals, 'o-', linewidth=2, markersize=10, color='steelblue')
    ax2.axhline(y=1e-4, color='r', linestyle='--', linewidth=2, label='1e-4 target')
    ax2.set_xlabel('Time Step (dt)', fontsize=11)
    ax2.set_ylabel('Final |F| (eV/Å)', fontsize=11)
    ax2.set_title('Final Force vs Time Step', fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=9)
    ax2.invert_xaxis()  # Larger dt on left
    
    # 3. Min force vs dt
    ax3 = plt.subplot(2, 3, 3)
    mins = [r['min_force'] for r in results]
    ax3.semilogy(dts, mins, 's-', linewidth=2, markersize=10, color='darkgreen')
    ax3.axhline(y=1e-4, color='r', linestyle='--', linewidth=2, label='1e-4 target')
    ax3.axhline(y=1e-6, color='g', linestyle='--', linewidth=2, label='1e-6 target')
    ax3.set_xlabel('Time Step (dt)', fontsize=11)
    ax3.set_ylabel('Best |F| (eV/Å)', fontsize=11)
    ax3.set_title('Best Force vs Time Step', fontsize=12, fontweight='bold')
    ax3.grid(True, alpha=0.3)
    ax3.legend(fontsize=9)
    ax3.invert_xaxis()
    
    # 4. Convergence rate (first 20% of steps)
    ax4 = plt.subplot(2, 3, 4)
    for i, r in enumerate(results):
        cutoff = int(len(r['steps']) * 0.2)
        ax4.semilogy(r['steps'][:cutoff], r['forces'][:cutoff], '-', 
                     linewidth=1.5, label=f"dt={r['dt']}", color=colors[i], alpha=0.8)
    ax4.set_xlabel('MD Steps', fontsize=11)
    ax4.set_ylabel('|F| (eV/Å)', fontsize=11)
    ax4.set_title('Initial Convergence (First 20%)', fontsize=12, fontweight='bold')
    ax4.grid(True, alpha=0.3)
    ax4.legend(fontsize=9, loc='best')
    
    # 5. Percentage below thresholds
    ax5 = plt.subplot(2, 3, 5)
    pct_1e4 = [100 * r['below_1e4'] / len(r['forces']) for r in results]
    pct_1e5 = [100 * r['below_1e5'] / len(r['forces']) for r in results]
    pct_1e6 = [100 * r['below_1e6'] / len(r['forces']) for r in results]
    
    x = np.arange(len(dts))
    width = 0.25
    ax5.bar(x - width, pct_1e4, width, label='< 1e-4', color='orange', alpha=0.8)
    ax5.bar(x, pct_1e5, width, label='< 1e-5', color='green', alpha=0.8)
    ax5.bar(x + width, pct_1e6, width, label='< 1e-6', color='blue', alpha=0.8)
    ax5.set_xlabel('Time Step (dt)', fontsize=11)
    ax5.set_ylabel('% of Steps Below Threshold', fontsize=11)
    ax5.set_title('Convergence Quality', fontsize=12, fontweight='bold')
    ax5.set_xticks(x)
    ax5.set_xticklabels([f"{dt:.3f}" for dt in dts], rotation=45)
    ax5.legend(fontsize=9)
    ax5.grid(True, alpha=0.3, axis='y')
    
    # 6. Mean force vs dt
    ax6 = plt.subplot(2, 3, 6)
    means = [r['mean_force'] for r in results]
    ax6.semilogy(dts, means, 'd-', linewidth=2, markersize=10, color='purple')
    ax6.set_xlabel('Time Step (dt)', fontsize=11)
    ax6.set_ylabel('Mean |F| (eV/Å)', fontsize=11)
    ax6.set_title('Average Force vs Time Step', fontsize=12, fontweight='bold')
    ax6.grid(True, alpha=0.3)
    ax6.invert_xaxis()
    
    plt.tight_layout()
    
    output_file = f"dt_comparison_{molecule}.png"
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"\n{'='*100}")
    print(f"  Comparison plot saved as: {output_file}")
    print(f"{'='*100}\n")

if __name__ == '__main__':
    main()
