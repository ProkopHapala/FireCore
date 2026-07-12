---
name: cpp-memory-ownership
description: >-
  Debug and implement C++ heap memory ownership in FireCore: Use when fixing memory crashes double-free, segfault
---

# C++ memory ownership (FireCore)

## Core rule

**One owner per heap array.** Other objects hold **borrowed aliases** (`T*` copies). Only the owner may `_dealloc`. Track with `bOwn*` flags set from `_bindOrRealloc` return value.

We use **explicit `dealloc()` / `clearFFs()`**, not strict RAII. Destructors on some classes (`Atoms`) call `dealloc()`, but MMFF teardown is orchestrated by `MolWorld::clearFFs()`.

## Macros (`cpp/common/macroUtils.h`)

- `_alloc` / `_realloc` / `_dealloc` — standard heap arrays
- `_bindOrRealloc(n, from, arr)` → `false` = borrow `from`, `true` = allocated `arr`
- `_dealloc` nulls **only the reference passed** — other aliases stay live

## Debug tools

| Tool | Use for |
|------|---------|
| **ASAN** (`Build-asan`, see `firecore-cpp-build` skill) | double-free, UAF, overflow — **first choice** |
| **DEBUG_ALLOCATOR** (`#define` in `macroUtils.h`) | alloc site ledger, invalid-free detection; call `debugAllocator_init()` |
| `tests/tMMFF/test_asan_minimal.py` | init/clear/Hessian smoke test |

Python `getBuffs()` = non-owning views. After `MMFF.clear()`, old NumPy arrays are UAF.

## When investigating a crash

1. Reproduce with ASAN build.
2. Find pointer in stack trace; grep `_dealloc` and bind sites.
3. Check `bOwn*` at bind (`initNBmol`, `bindOrRealloc`, `setOptimizer`).
4. Verify `clearFFs()` frees **owner before borrowers** (`ffl` before `nbmol`/`opt`).
5. Watch interior pointers: `apos`⊂`DOFs`, UFF `fbon`⊂`fint` — never free views.

## Adding a new array

1. Pick single owner.
2. Use `_bindOrRealloc` if own-or-borrow; set `bOwn*`.
3. Free only in owner's `dealloc()`.
4. Run ASAN + `init → clear → init`.

## Full reference

**[doc/dev_notes/MMFF/Memory_Ownership_and_Deallocation.md](../../../doc/dev_notes/MMFF/Memory_Ownership_and_Deallocation.md)**
