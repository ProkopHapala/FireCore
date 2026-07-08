---
name: running-tests
description: Use when running tests or writing new test scripts for SPAMMM
trigger:
  glob:
    - "**/tests/**/*"
    - "**/*test*.py"
    - "**/pytest.ini"
    - "**/conftest.py"
    - "**/doc/TEST_DESIGN.md"
---

## Authority

**SSOT:** `doc/TEST_DESIGN.md`. This skill is the operational summary.

---

## Three review levels

| Level | Name | Mechanism |
|-------|------|-----------|
| **L0** | Automatic | `assert`, `TopologyDiff`, `ref_data/` |
| **L1** | Agentic | `.out` + `.log` in `debug/<script>/` |
| **L2** | Visual | `.png` via `--visual` or `--develop` |

**If you can produce a good reference, you don't need L1** — assert against `ref_data/`. See skill:`reference-data`.

---

## Run modes

| Mode | Command | Use when |
|------|---------|----------|
| **Routine** | `pytest -m "not slow"` | Regression — did we break established code? |
| **Develop** | `pytest path --develop -s` | New feature/module — agent writes AND reviews artifacts |

```bash
pytest -m "not slow"                                    # routine (L0 only)
pytest tests/topology/test_editing_ops.py --develop -s  # L0+L1+L2
pytest tests/topology/test_editing_ops.py --review -s   # L0+L1
pytest tests/topology/test_editing_ops.py --visual -s   # L0+L2
pytest -m "gpu and not slow"
pytest --update-refs tests/test_folded_relax.py
```

### Develop-mode agent contract

1. Run **in foreground** — **never** `| grep`, `| tail`, `| head`, `&`, or background.
2. Use `-s` so stdout is visible (user feedback during long runs).
3. After run, read every `REVIEW: debug/...` path from stdout.
4. Read **`.out` first** (curated evaluation). If issue spotted, read **`.log`** (trace).
5. Report findings to user; include `.png` paths for human L2 review.

---

## File classes

| Pattern | pytest? | Purpose |
|---------|---------|---------|
| `test_*.py` | Yes | L0 (+ optional L1/L2 via flags) |
| `testplot_*.py` | No | Pure visual demos: `python tests/testplot_foo.py` |
| `run_*.py` | No | CLI utilities |
| `helpers/` | No | Imported utilities |

---

## Writing tests — what goes where

**Before coding:** What would I need to see to know this is wrong?

| Channel | Content |
|---------|---------|
| **stdout** | Progress, key scalars, `REVIEW: debug/...` pointers |
| **`.out`** | Intent, `AtomicGraph.format_table()`, metrics, agent checklist |
| **`.log`** | Step-by-step trace, internal variables, `debug_print` at level 3 |

Use `make_review('test_foo')` fixture + `ReviewSession` from `tests/helpers/review.py`.

**L0 is mandatory.** Even `assert np.isfinite(E).all()` or `diff.assert_counts(...)`. Permissive asserts are OK; absent asserts are not.

**Molecular dumps:** `graph.format_table(pos=False, neighbors=True)` — one atom per line. Use `.mol2` when geometry matters.

**Arrays:** `debug_summarize_array(x)` → shape, dtype, min, max. Check `np.isfinite`. Invariants: skill:`numerical-parity`, skill:`forcefield-validation`.

---

## Artifacts

```
debug/<script_stem>/
  <test_func>.out
  <test_func>.log
  <test_func>.png          # --visual or --develop
```

Gitignored except `debug/README.md`. Agents read by explicit path — gitignore does not block filesystem access. Never add `debug/` to `.cursorignore`.

---

## Fixtures (`conftest.py`)

- `make_review` — factory → `ReviewSession`
- `visual_output_dir` — `debug/<script>/` when `--visual` or `--develop`
- `review_dir`, `review_enabled`, `develop_mode`
- `xyz`, `substrate`, `dat`, `update_refs`

---

## Markers

`slow`, `gpu`, `visual`, `review` — see `pytest.ini`.

---

## Reference data

Git-tracked `tests/ref_data/*.ref.json` + `*.ref.xyz` for established regression. Regenerate: `pytest --update-refs`. Details: skill:`reference-data`.

---

## Helpers

- `helpers/review.py` — `ReviewSession`, `review_trace`
- `helpers/topology_test.py` — `TopologySnapshot`, `TopologyDiff`, `render_before_after`
- `helpers/parity.py` — `assert_parity`, `overlay_plot`
- `helpers/geometry.py`, `helpers/scan.py`, `helpers/folded_rigid.py`

---

## Execution time

Each test < 1s default; >1s → `@pytest.mark.slow`. Default: `pytest -m "not slow"`.

---

## Verbosity

`spammm.globals.debug_print(level, msg)` gated by `SPAMMM_VERBOSITY` env or `DEBUG_PRINT_LEVEL`. `--develop` bumps levels via `set_develop_mode()`.
