---
name: visual-debugging
description: Use when creating diagnostic plots, visualizations, or headless visual tests for debugging
trigger:
  glob:
    - "**/tests/**/*"
    - "**/*test*.py"
    - "**/*debug*.py"
    - "**/*benchmark*.py"
    - "**/*plot*.rs"
---

## Authority

**SSOT:** `doc/TEST_DESIGN.md`. This skill covers **L2 visual** and headless rendering.

---

## Three levels (reminder)

| Level | Output |
|-------|--------|
| L0 | `assert` / `TopologyDiff` |
| L1 | `.out` / `.log` — skill:`running-tests` |
| **L2** | `.png` — this skill |

---

## Rules

1. **Reuse helpers:**
   - `tests/helpers/topology_test.py` — `TopologySnapshot`, `TopologyDiff`, `render_before_after`
   - `tests/helpers/parity.py` — `overlay_plot()`, `assert_parity()`
   - `tests/helpers/geometry.py` — bond/angle checks
   - `spammm/GUI/VispyUtils.py` — interactive 3D only (not headless tests)

2. **Integrate L2 into pytest** when an assertable core exists:
   - `--visual` or `--develop` → `visual_output_dir` fixture → `debug/<script>/<test_func>.png`
   - L0 assertions always run; PNG only when fixture is not `None`

3. **Pure visual demos** (no assertable core): `testplot_*.py`, run via `python tests/...`. Not collected by pytest.

4. **Test backend logic, not GUI widgets.** Simulate via API:
   - `graph.pick_atom(pos)`, `backend.add_ring(q, r)`, etc.

5. **Output location:** `debug/<script_stem>/` only. Never `/tmp/` or repo root. Report exact paths.

6. **Plot style:**
   - Reference: `ls=':'`, `lw=1.5`
   - Model: `ls='-'`, `lw=0.5`
   - Residual twin axis: `(model - ref) * 100`
   - RMSE/MaxErr box: upper-left, monospace

7. **No `plt.show()`** in library code. Agg backend for headless.

8. **Foreground execution.** Never hide output (`| tail`, `| head`, `&`). Full stdout visible.

9. **CRITICAL — equal axes for spatial maps:** any `imshow` of XY/XZ/YZ density or potential
   **MUST** use `aspect='equal'` so **1 Å along x = 1 Å along y/z**. Never `aspect='auto'`
   (distorts spheres → ellipses). Tip/molecule diagnostics: plot in **molecular frame**
   (before pad+roll), mark where 1D cuts are taken on the 2D map.

---

## Topology editing pattern

- **L0:** `TopologySnapshot` / `TopologyDiff.assert_counts` — stable `_id`, not array indices
- **L1:** `make_review` + `graph.format_table()` — skill:`running-tests`
- **L2:** `render_before_after` with diff coloring (green=added, red=removed, blue=new bonds)
- Annotations: cursor crosshair, `id:element` labels, selection halos

```bash
pytest tests/topology/test_editing_ops.py -s                 # L0
pytest tests/topology/test_editing_ops.py --develop -s       # L0+L1+L2
```

---

## Develop mode

When testing a **new** feature, use `--develop`: plots on by default so the user can review immediately without a second run.

---

## Physics-First Plot Design (MANDATORY)

**Think like a physicist or engineer when designing test plots.** A plot is not
a checkbox — it is a diagnostic instrument. Its purpose is to reveal whether the
code produces physically correct behavior across a meaningful range of
parameters. A plot that shows a constant, all-zero, or tautological result is a
**failed plot** — it tells you nothing and wastes the reviewer's time.

### What a diagnostic plot MUST do

1. **Illustrate a real physical or engineering phenomenon.** The plot should
   show a recognizable relationship from physics or engineering: a flash curve
   (quality vs enthalpy), a conversion-vs-recycle-ratio curve, an isentropic
   temperature rise, a pressure-drop-vs-flow curve, a phase diagram, a
   convergence trajectory. If you cannot name the physical principle the plot
   demonstrates, the plot is wrong.

2. **Sweep a relevant range of the independent variable.** Do not plot a
   single operating point. Sweep the parameter that matters:
   - Flash: sweep enthalpy (or temperature) across the two-phase region
   - Compressor: sweep pressure ratio from 1 to ~5
   - Recycle: sweep recycle ratio R/F from 0 to ~3
   - Valve: sweep opening from 0.01 to 1.0
   - Reactor: sweep feed composition or temperature
   The range must span the interesting transition (e.g. subcooled → two-phase →
   superheated), not just one side of it.

3. **Show a non-trivial trend.** The dependent variable must change visibly
   across the sweep. If the plot is flat, either:
   - The parameter range is wrong (too narrow, or on the wrong side of a
     transition)
   - The code has a bug (silent fallback, wrong units, dead code path)
   - The chosen variable doesn't couple to the physics (e.g. sweeping
     temperature when enthalpy is the relevant input)
   In all three cases, **fix the plot or fix the code** — do not ship a flat
   plot.

4. **Compare against a known reference when possible.** Overlay the
   analytical formula, a textbook curve, or a reference implementation:
   - Isentropic T_out = T_in * r^((gamma-1)/gamma) for a compressor
   - Flash quality x = (h - h_l) / (h_g - h_l) in the two-phase region
   - Overall conversion = 1 - (1-X_sp)/(1+rf*X_sp) for a recycle reactor
   - Poiseuille flow F = pi*r^4*dP / (8*mu*L) for laminar pipe flow
   The comparison reveals whether the code matches theory — this is the whole
   point of the plot.

5. **Label axes with physical units.** "Specific enthalpy (kJ/mol)", not
   "h". "Compression ratio P_out/P_in", not "r". "Interface height (m)",
   not "h_int". A reviewer should understand the plot without reading the
   code.

### Anti-patterns (DO NOT ship these)

- **Constant plot:** All values identical across the sweep. This means the
  swept variable doesn't affect the output — either a code bug or a wrong
  choice of variable. Example: plotting "liquid fraction vs feed temperature"
  when the flash always returns all-gas because enthalpy is computed wrong.

- **All-zero plot:** Every value is 0.0. This usually means a silent fallback
  (e.g. `UnsupportedModel` caught and ignored, leaving moles at zero). The
  plot hides the bug instead of revealing it.

- **Single-point "plot":** One dot on an empty canvas. Not a sweep, not a
  trend, not diagnostic. Always sweep at least 10-20 points.

- **Tautological plot:** The output is mathematically forced to equal the
  input (e.g. plotting "flow vs flow" or "T vs T"). Tells you nothing about
  the code's correctness.

- **Convergence-only plot:** "Iterations vs parameter" with no physical
  content. Convergence count is useful metadata but not the primary
  diagnostic — the physical result (temperature, composition, flow) is.

- **Binary indicator plot:** A scatter of 0/1 dots showing "converged vs not"
  with no physical variable on the y-axis. This is a status report, not a
  diagnostic plot. If convergence fails, show *why* — plot the residual or
  the physical quantity that is oscillating.

### Design checklist (run this before generating any plot)

1. **What physical principle does this plot demonstrate?** Name it in one
   sentence. If you can't, redesign.
2. **What is the independent variable, and what range does it span?** Does
   the range cross a transition (phase change, flow reversal, convergence
   threshold)?
3. **What is the expected trend?** Sketch it mentally before running the
   code. If the code's output doesn't match the expectation, investigate.
4. **Is there an analytical or textbook reference to overlay?** If yes,
   include it. The gap between model and reference IS the diagnostic.
5. **Are the axes labeled with units?** If not, fix.
6. **Does the plot show a non-trivial trend?** If it's flat or constant,
   stop and debug — do not ship.
