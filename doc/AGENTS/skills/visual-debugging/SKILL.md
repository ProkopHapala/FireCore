---
name: visual-debugging
description: Use when creating diagnostic plots, visualizations, or headless visual tests for debugging
trigger:
  glob:
    - "**/tests/**/*"
    - "**/*test*.py"
    - "**/*debug*.py"
    - "**/*benchmark*.py"
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
