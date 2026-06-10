
## Core Philosophy: Scientific & Performance Focus

We develop rigorous scientific software where debuggability, numerical correctness, and physical consistency are paramount. Follow these principles:

1. **Debuggability:** Code must be transparent, inspectable, and prioritize trace-ability over user experience (never hide issues).
2. **Simplicity:** Clear, clean, direct logic. Elegant design that avoids branching, excessive special conditions, and defensive abstractions.
3. **Performance:** Streamlined execution with minimal overhead (avoid Python loops), data-oriented memory layouts, and cache awareness.

## Rule 1 — Fail Loudly

- **No silent fallbacks**, catch-all passes, or try-except blocks that mask bugs. Unexpected states must terminate with explicit errors and full stack traces. See `general-debug-guidelines.md`.
- **Root Cause Identification:** Find and fix the fundamental issue. Never apply "quick-fixes" that hide the root cause.

## Rule 2 — Surgical Edits & Simplicity

- **Minimum Intervention:** Write only the code necessary to solve the task. Touch only what is required; never perform unrelated formatting, cleanup, or aesthetic edits on adjacent code.
- **No Guessing:** If requirements, behavior, or architecture are ambiguous, stop and ask for clarification.
- **Report Problems Immediately:** Encountered ambiguities, unexpected errors, or risky decisions must be written in the chat at once. Never apply silent workarounds or unilateral fixes.
- **Strict Checkpointing:** After every significant step, summarize what changed, what was verified, and what remains unresolved.
- **Preservation & Backups:** Create a backup copy before major module changes. Comment out old, experimental, or deprecated code using `#` or `//` instead of deleting it to allow instant reversion. Mark unfinished code clearly with `# TODO` or `# DEBUG`.

## Rule 3 — Reusable Architecture

- **Inventory First:** Thoroughly review reference source-code files to identify existing functions, modules, and data structures before writing anything from scratch. Use `CODEMAP.md` and `/doc/topical_audit/topical_audit.md` for guidance.
- **Composability Over Bloat:** Build integrated systems, not isolated scripts. Refactor into reusable functions in shared modules.
- **Separation of Concerns:**
   - Separate compute algorithms from plotting/diagnostics (no plotting in core libraries).
   - Separate GUI, CLI test scripts, and backend modules. Test scripts are thin wrappers that call functions from shared backend modules.
   - Consolidate related test scripts into one with CLI routing parameters for different execution paths.
- **Generalization Over Duplication:** Try to generalize an existing function if it almost fits your needs. If generalization requires risky major changes that threaten backward compatibility, **stop and report it immediately for approval.**

## Rule 4 — Test-Driven Development & Validation

- **Numerical Range Sanity:** Strategically place checks throughout calculations to ensure values are within reasonable limits and are not `NaN`, infinity, or unexpected zeros.
- **Test on Completion:** Run validation tests immediately after any code modification. Never claim code works unless tests run successfully.
- **Physical & Analytical Parity:** Define how correctness will be verified *before* coding via parity checks against reference code, known analytical solutions, physical conservation laws, symmetry checks, or known physical limits. See `numerical-parity/SKILL.md`.
- **Foreground Execution:** Run tests synchronously with full output. Never use background commands, pipes (`| tail`, `| head`). Full stdout must be visible.
- **Visual Review & Diagnostics:** Use shared utilities for plotting, debugging, and diagnostics instead of ad-hoc code. See `visual-debugging/SKILL.md` for `plotUtils.py`, `VispyUtils.py`, `TestUtils.py`, and `testUtils.h`.
- **Invoke Relevant Skills:** When task matches skill description (numerical-parity, visual-debugging, gpu-debugging, forcefield-validation, port-to-opencl), invoke the skill tool to get detailed guidance.

## Rule 5 — Performance Optimization

- **Minimal Orchestration:** Keep Python orchestration minimal. Push heavy computations into optimized C/C++/OpenCL/CUDA/Compute Shader kernels.
- **Memory Optimization:** Prefer flat, contiguous arrays and data-oriented layouts. Be explicit about dtypes and shapes. Preallocate and reuse buffers; avoid repeated allocations in hot paths.
- **Low-Level & GPU Kernel Guidelines:** Design around memory latency, prefer gather over scatter, minimize branching/atomics/synchronization, maximize shared/local memory usage, avoid unnecessary host-device transfers. See `port-to-opencl/SKILL.md`.

## Rule 6 — Concise Style

- **No Micro-Abstractions:** Do not create 1 line function stubs or wrappers. If it is simple, inline it.
- **Clean Interfaces:** Avoid passing excessive numbers of arguments. Group related state into structs/dicts, or utilize globals and class properties. Use default named arguments to avoid long call strings.
- **Compact Layout:** Prefer compact code with long lines and minimal empty lines or whitespace. Avoid line wrapping that disrupts the readability of expressions; assume infinite line length.
- **Naming & Comments:** Use short, clear variable names for math/physics symbols (e.g., `E_tot`, `T_ij`, `m_i`). Avoid comments that state the obvious; use them for intent, rationale, or math/physics derivations. Place inline comments behind the code line.
- **Language-Specific Rules:**
- **C++:** Use `printf` for debugging instead of `std::cout`. Prefer plain C arrays (`double*`) in hot paths.
- **Doxygen:** Document using `///`; avoid `/* ... */` formatting.
- **Parity Work:** When mirroring features across languages (e.g., Python $\leftrightarrow$ JS), explicitly cite the reference file and function in the comments.

## Practical Navigation, Compilation, testing Protocols

- **Repository Navigation:** Review `CODEMAP.md` for structure and build instructions.
- **Test Location:** Place all test scripts within `/test`.
- **Automation Scripts:** Use provided `run.sh`/`make.sh` scripts in the test directory; never invoke `make` directly if helpers exist. Run tests from inside the test directory to ensure paths are set.
