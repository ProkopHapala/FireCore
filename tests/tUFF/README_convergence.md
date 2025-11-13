# UFF Double Precision Convergence Testing

## Quick Start

### 1. Edit `run_multi.sh` to configure your test:

```bash
# Configuration for convergence test
MOLECULE="H2O"           # Change to your molecule name (e.g., "xylitol", "CH4")
LOGFILE="log_${MOLECULE}_convergence.txt"

# Adjust simulation parameters:
python3 run_throughput_UFF.py \
  --xyz_name data/xyz/${MOLECULE}.xyz \
  --nSys 1 \
  --bUFF 1 \
  --bGridFF 1 \
  --gridnPBC "(1,1,0)" \
  --loops 1 \
  --perframe 10000 \    # Number of MD steps
  --perVF 100 \
  --Fconv 1e-6 \        # Convergence target
  --dt 0.001            # Time step (smaller = more stable, slower)
```

### 2. Run the test:

```bash
./run_multi.sh
```

**That's it!** The script will:
- Run the simulation
- Automatically extract convergence data
- Generate statistics
- Create plots
- Save everything with the molecule name

## Output Files

After running, you'll get:

1. **`log_<MOLECULE>_convergence.txt`** - Full simulation log
2. **`convergence_<MOLECULE>.txt`** - Extracted data (step, force)
3. **`convergence_<MOLECULE>.png`** - 4-panel convergence plot

## Time Step Guidelines

| Target Precision | Recommended dt | Expected Steps |
|------------------|----------------|----------------|
| 1e-4 eV/Å       | 0.01           | ~10,000        |
| 1e-5 eV/Å       | 0.005          | ~50,000        |
| 1e-6 eV/Å       | 0.001          | ~100,000+      |

## Two-Stage Approach (Recommended for 1e-6)

For best efficiency to reach 1e-6:

**Stage 1:** Rough convergence
```bash
MOLECULE="H2O"
dt=0.01, perframe=10000  # Reaches ~1e-4 quickly
```

**Stage 2:** Fine convergence (restart from Stage 1 output)
```bash
dt=0.001, perframe=50000  # Final push to 1e-6
```

## Example: Testing Different Molecules

### H2O (simple, fast):
```bash
MOLECULE="H2O"
perframe=10000, dt=0.001
```

### Xylitol (complex, slower):
```bash
MOLECULE="xylitol"
perframe=50000, dt=0.01  # Stage 1
perframe=100000, dt=0.001  # Stage 2
```

## Troubleshooting

### Simulation crashes after ~2000 steps
- **Cause:** Cleanup bug (double-free)
- **Fix:** Use `--loops 1` with large `--perframe` (already configured)

### Forces oscillate, won't converge
- **Cause:** Time step too large
- **Fix:** Reduce `dt` (e.g., 0.01 → 0.001)

### Convergence too slow
- **Cause:** Time step too small
- **Fix:** Increase `dt` for initial stages, then decrease for final convergence

## What the Plots Show

1. **Top-left:** Full convergence history (log scale)
2. **Top-right:** Initial convergence (first 20% of steps)
3. **Bottom-left:** Final convergence (last 30% of steps)
4. **Bottom-right:** Distribution histogram

## Success Criteria

✅ **Double precision working:** Achieved < 1e-4 eV/Å  
✅ **Good convergence:** Achieved < 1e-5 eV/Å  
✅ **Excellent convergence:** Achieved < 1e-6 eV/Å
