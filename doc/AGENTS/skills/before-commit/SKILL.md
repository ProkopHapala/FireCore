---
name: before-commit
description: Classify git changes before commit (git add vs skip); code and docs only, no data/output/cache. Use before commit or with `git commit -a`.
---

# Before Commit

## Defaults

1. User commits with **`git commit -a`** → modified **tracked** files only.
2. **Untracked** paths need explicit **`git add`**.
3. **Include:** source, scripts, tests (code), docs.
4. **Exclude:** data, computed output, caches, scratch, plots.
5. **Do not commit** unless asked.

## Workflow

```bash
git status --porcelain -u
git diff --stat
```

Scope to current task when user says "today" / "this work"; else full tree. Classify paths → report briefly. Propose `.gitignore` lines only if asked; do not edit unless asked.

## ADD

| Kind | Examples |
|------|----------|
| Source | `*.py` `*.cpp` `*.h` `*.cl` `*.mjs` `*.js` in code dirs |
| Tests / runners | `test_*` `run_*` `generate_*` scripts |
| Generators | `bootstrap_*` `generate_*` (code that builds local data) |
| Docs | `*.md` in `doc/`, `README.md`, `AGENTS.md`, `CODEMAP.md` |
| Build | `CMakeLists.txt`, `Makefile` — if part of the change |

## SKIP

| Kind | Examples |
|------|----------|
| Data / output | `*.npz` `*.npy` `*.xyz` `*.cube`, run status/meta JSON, fixture binaries |
| Scratch | `work_*` `out_*` temp dirs, generated geometry in test trees |
| Plots | `*.png` `*.bmp`, `plots_*` |
| Caches | `.pytest_cache/` `__pycache__/` `.*_cache*/` |
| Local / IDE | `.cursor/` `.clangd` `compile_commands.json` |
| Backups | `*_bak*` `* copy.*` |

**Gray:** small reference inputs (CIF, params) — skip unless user wants them. Fixture `README.md` — add. Unrelated modified tracked files — flag ( `-a` picks them up).

## Fixtures

Commit test code + generator scripts + fixture README. Do not commit generated artifacts under `fixtures/` or `ref_data/` unless user overrides.

## Report (keep short)

```markdown
## `git commit -a` (task-related modified tracked)
- …

## `git add` (untracked code & docs)
- …

## Skip
| path | reason |

## Watch (`-a` extras)
- …

## Command
git add …
```

Group by subsystem. Concrete paths. One `git add` block for ADD only.
