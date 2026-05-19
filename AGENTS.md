
# AI_AGENTS.md

## Core Philosophy: Scientific Focus

We develop rigorous scientific software. The ability to instantly spot bugs or numerical/physical inconsistencies is paramount.

Code architecture and design decisions are governed by three pillars, in order of priority:

1. **Debuggability:** Code must be transparent and inspectable.
2. **Simplicity:** No over-engineering. Clear, clean, and direct logic. Elegant design that avoids branching and special conditions.
3. **Performance:** Physics simulations require streamlined execution, low-level efficiency, and minimal overhead (e.g. avoid loops in Python).

## Rule 1 — Fail Loudly (No Silent Fallbacks)

* **Never Make Silent Fallbacks:** If an expected code path or numerical assumption fails, you must raise an explicit exception, print a precise error message, and abort execution immediately.
* **Zero Tolerance for Hidden Errors:** Silent fallbacks or catch-all passes make debugging complex simulations impossible. If something unexpected happens, crash loudly.

## Rule 2 — Surgical edits & simplicity

* **Surgical Changes:** Touch only what is strictly necessary to solve the problem. Match existing codebase conventions and formatting perfectly. Never "clean up" or modify adjacent code, comments, or formatting outside your explicit scope.
* **Simplicity First:** Write the minimum code required to solve the task. No speculative features, no defensive abstractions for large single-use code.
* **No Guessing:** If you encounter ambiguity or uncertainty, stop and ask for clarification. Never make assumptions that could lead to incorrect behavior.
* **Strict Checkpointing:** Summarize what was done, what is verified, and what remains after every significant step.

## Rule 3 — Avoid Code Duplication (Read Before Write)

* **Inventory First:** You must always thoroughly review the provided reference source-code files to identify existing functions, modules, and data structures available for reuse, before writing anything from scratch.
* **Composability Over Bloat:** We build integrated systems, not isolated scripts. Compose existing, well-tested modules to keep the codebase lean.
* **Generalization Over Duplication:** If you find an existing function that *almost* fits your needs, try to generalize it rather than duplicating code. However, do not make major changes that risk breaking backwards compatibility. If this dilemma arises, **report it to me immediately** and ask for approval.

## Rule 4 — Test-Driven Development

* **Physical & Analytical Parity:** Every functionality must have a rigorous validation test. Before coding, define how to prove correctness (e.g., via parity checks against reference code, validation against known analytical solutions for simple cases, physical conservation laws, or symmetry).
* **Verbose Debug Prints:** Place debug prints of intermediate results throughout the code. Use a gated `debug_print` function and a verbosity level system from global modules (e.g., `global.py`, `global.h`) for this.
* **Are Numbers Reasonable?:** Place several strategic checks throughout the calculation to verify if values are within a reasonable range, and are not `NaN`, infinity, or unexpected zeros.
* **Test on Completion:** Run the relevant test or validation script immediately after finishing any code modification. Do not use "| tail/head" when running the tests, I want to see full stdout output. 
* **Visual Validation:** I test visually - prepare diagnostic plots via `matplotlib` plots saved as `.png` files for my revision, and report where can I find them.
* **Structured Test Artifacts:** Assume we will run dozens of tests before a task is complete. Group all debugging, benchmarking, and testing outputs into a dedicated, clearly named, and numbered folder. Do not clutter the root directories. Explicitly report the location of these results.

## Rule 5 - Concise style

* Do not create function stubs or abstractions that do too little (e.g., 1–3 lines of code). If it is simple, inline it.
* Do not pass too many arguments over function interfaces; group them into structs/dicts if needed, or utilize globals and class properties.
* Prefer compact code with long lines and minimal whitespace.
* Prefer one-line constructs. Assume infinite line length; do not break function calls onto multiple lines just because of length.