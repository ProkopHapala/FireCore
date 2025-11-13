#!/usr/bin/env python3
"""
Plot convergence data from H2O UFF double precision test
"""
import numpy as np
import matplotlib.pyplot as plt

# Load data
data = np.loadtxt('convergence_data.txt')
steps = data[:, 0]
forces = data[:, 1]

# Create figure with multiple subplots
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# 1. Full convergence plot (log scale)
ax = axes[0, 0]
ax.semilogy(steps, forces, 'b-', linewidth=0.5, alpha=0.7)
ax.axhline(y=1e-4, color='r', linestyle='--', label='1e-4 (float limit)')
ax.axhline(y=1e-6, color='g', linestyle='--', label='1e-6 (target)')
ax.set_xlabel('MD Steps')
ax.set_ylabel('|F| (eV/Å)')
ax.set_title('Full Convergence History (Log Scale)')
ax.grid(True, alpha=0.3)
ax.legend()

# 2. Linear scale (first 20k steps)
ax = axes[0, 1]
mask = steps <= 20000
ax.plot(steps[mask], forces[mask], 'b-', linewidth=0.8)
ax.axhline(y=1e-4, color='r', linestyle='--', label='1e-4')
ax.set_xlabel('MD Steps')
ax.set_ylabel('|F| (eV/Å)')
ax.set_title('Initial Convergence (Linear Scale, First 20k)')
ax.grid(True, alpha=0.3)
ax.legend()

# 3. Last 30k steps (zoomed)
ax = axes[1, 0]
mask = steps >= (steps[-1] - 30000)
ax.semilogy(steps[mask], forces[mask], 'b-', linewidth=0.8)
ax.axhline(y=1e-5, color='orange', linestyle='--', label='1e-5')
ax.axhline(y=1e-6, color='g', linestyle='--', label='1e-6 (target)')
ax.set_xlabel('MD Steps')
ax.set_ylabel('|F| (eV/Å)')
ax.set_title('Final Convergence (Last 30k Steps)')
ax.grid(True, alpha=0.3)
ax.legend()

# 4. Histogram of forces
ax = axes[1, 1]
ax.hist(np.log10(forces), bins=100, edgecolor='black', alpha=0.7)
ax.axvline(x=np.log10(1e-4), color='r', linestyle='--', label='1e-4')
ax.axvline(x=np.log10(1e-6), color='g', linestyle='--', label='1e-6')
ax.set_xlabel('log10(|F|) [eV/Å]')
ax.set_ylabel('Count')
ax.set_title('Distribution of Forces')
ax.grid(True, alpha=0.3)
ax.legend()

plt.tight_layout()
plt.savefig('h2o_convergence_analysis.png', dpi=150, bbox_inches='tight')
print("Plot saved as: h2o_convergence_analysis.png")

# Print summary statistics
print("\n=== Summary Statistics ===")
print(f"Total steps: {len(steps)}")
print(f"Min |F|: {forces.min():.3e} eV/Å (step {steps[forces.argmin()]:.0f})")
print(f"Max |F|: {forces.max():.3e} eV/Å")
print(f"Final |F|: {forces[-1]:.3e} eV/Å")
print(f"Mean |F|: {forces.mean():.3e} eV/Å")
print(f"Median |F|: {np.median(forces):.3e} eV/Å")

# Convergence milestones
print("\n=== Convergence Milestones ===")
milestones = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6]
for threshold in milestones:
    idx = np.where(forces < threshold)[0]
    if len(idx) > 0:
        first_step = steps[idx[0]]
        count = len(idx)
        percent = 100 * count / len(steps)
        print(f"|F| < {threshold:.0e}: First at step {first_step:.0f}, achieved {count} times ({percent:.1f}%)")
    else:
        print(f"|F| < {threshold:.0e}: Never achieved")

# Check oscillation in last 10k steps
last_10k = forces[-10000:]
print(f"\n=== Last 10k Steps ===")
print(f"Min: {last_10k.min():.3e} eV/Å")
print(f"Max: {last_10k.max():.3e} eV/Å")
print(f"Mean: {last_10k.mean():.3e} eV/Å")
print(f"Std Dev: {last_10k.std():.3e} eV/Å")
print(f"Oscillation range: {last_10k.max()/last_10k.min():.1f}x")

plt.show()
