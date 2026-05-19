# AGENTS.md

## Core Philosophy: Scientific Software

We develop rigorous scientific software where numerical correctness, physical consistency, and debuggability are paramount.

Code design priorities:

1. **Debuggability** — code must be transparent and inspectable.
2. **Simplicity** — avoid over-engineering, unnecessary abstractions, and special cases.
3. **Performance** — physics simulations require efficient low-level execution and minimal overhead.

## Rule 1 — Fail Loudly

- Silent fallbacks and catch-all exception handling are forbidden.
- Unexpected states, invalid assumptions, and numerical failures must terminate with explicit errors.
- Never hide errors behind automatic recovery logic.

## Rule 2 — Surgical Changes

- Modify only what is necessary for the task.
- Preserve existing style, formatting, and architecture.
- Do not perform unrelated cleanup or refactoring.
- If requirements or behavior are ambiguous, stop and ask.
- After significant steps, summarize:
  - what was changed,
  - what was verified,
  - what remains unresolved.

## Rule 3 — Read Before Write

- Before writing new code, inspect existing modules, functions, and utilities for reuse.
- Prefer extending existing functionality over duplicating logic.
- Avoid introducing isolated or redundant implementations.
- If proper generalization would require risky architectural changes, ask first.

## Rule 4 — Validation & Debugging

- Before coding, define how correctness will be verified:
  - analytical solutions,
  - conservation laws,
  - symmetry checks,
  - parity with reference implementations,
  - known physical limits.
- Never claim code works unless relevant tests were actually run successfully.
- Add strategic numerical sanity checks (`NaN`, `Inf`, divergence, unreasonable magnitudes).
- Use gated `debug_print` statements with verbosity levels.
- Run validation/tests immediately after modifications and show full stdout.
- Generate diagnostic plots (`.png`) for visual inspection when useful and report their location.
- Store test outputs, plots, benchmarks, and debug artifacts in organized numbered directories (e.g. `tests/003_case_name/`).

## Rule 5 — Performance

- Keep Python orchestration minimal; avoid Python loops in performance-critical paths.
- Push heavy computation into optimized C/C++/OpenCL/CUDA kernels.
- Prefer flat contiguous arrays and data-oriented memory layouts.
- Consider cache locality, memory bandwidth, and synchronization costs explicitly.
- On GPU:
  - prefer gather over scatter,
  - minimize branching and atomics,
  - use shared/local memory efficiently,
  - avoid unnecessary host-device transfers.

## Rule 6 — Concise Code

- Prefer simple direct code over excessive abstractions.
- Inline trivial logic instead of creating tiny wrapper functions.
- Avoid passing excessive numbers of arguments; group related state logically.
- Prefer compact readable code with minimal unnecessary whitespace.
- Keep control flow obvious and easy to inspect during debugging.