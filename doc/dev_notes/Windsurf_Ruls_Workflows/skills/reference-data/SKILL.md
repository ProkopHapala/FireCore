---
name: reference-data
description: Use when creating or updating pytest reference files in tests/ref_data/
trigger:
  glob:
    - "**/tests/ref_data/**/*"
    - "**/*ref.json"
    - "**/*ref.xyz"
    - "**/conftest.py"
    - "**/doc/TEST_DESIGN.md"
---

## When to use references vs agent review

| Situation | Approach |
|-----------|----------|
| Established behavior — "did we break it?" | `tests/ref_data/` + L0 assert |
| New feature — no golden answer yet | L1 agent reads `.out` (skill:`running-tests`) |
| Reference is good enough | Promote to `ref_data/`, replace L1 with L0 |

**Rule:** If you can produce a tight reference, you don't need L1 agent judgment — assert `rmse(result, ref) < tol`.

Most active development has **no reference yet**. Agent judgment + good `.out`/`.log` artifacts is the strength of agentic debugging. References accumulate over time as features stabilize.

---

## File layout

```
tests/ref_data/
  h2o_nacl.ref.json       # scalar properties, tolerances, test linkage
  h2o_nacl.ref.xyz        # geometry snapshot
```

- **Git-tracked** — permanent regression data, not debug artifacts
- JSON includes `test_func`, `test_module` for bidirectional linking
- Compare with **physical tolerances**, not exact float equality

---

## Regenerating references

After intentional physics changes:

```bash
pytest tests/test_folded_relax.py --update-refs
```

Always commit updated `.ref.json` and `.ref.xyz` together.

---

## What to snapshot

**Good for references:**
- Final geometry (`.ref.xyz`) when relaxation/MD is stable
- Scalar invariants: `z_rel`, force norms, torque, energy, bond lengths
- Parity metrics: RMSE vs analytical/brute reference

**Poor for references (use L1 instead):**
- Subjective chemistry (Kekule aesthetics)
- Plot appearance
- Cases still under active development

---

## Creating a new reference

1. Run test in develop mode; verify `.out` and plots look correct
2. Identify **minimal** assertable scalars (not everything in the dump)
3. Save via test's `save_reference()` or `--update-refs` path
4. Set tolerances from observed Fp32/Fp64 noise (skill:`numerical-parity`)
5. Document in test docstring what the reference represents

---

## Helpers

- `tests/helpers/folded_rigid.py` — `save_reference`, `compare_to_reference`
- `tests/helpers/parity.py` — `assert_parity`, `rmse`

See `tests/TEST_RESULTS.md` for logged results.
