---
trigger: always_on
---

- Focus on debuggability, not user experience => do not try to hide issues.  
- Add debug prints to track variables and flow (we remove them later).  
- Never silently handle errors => fail loudly.  
- Avoid try-except; crashes with stack trace are helpful.  
- Comment out old code instead of deleting it.  
- Make small, testable changes and run often.  
- Log function entries and key conditions when helpful.  
- Mark experimental or unfinished code clearly (e.g. # TODO, # DEBUG).

---

## Rule 1 — Fail Loudly

- Silent fallbacks and catch-all exception handling are forbidden.
- Unexpected states, invalid assumptions, and numerical failures must terminate with explicit errors.
- Never hide errors behind automatic recovery logic or hardcoded outputs.

## Rule 2 — Surgical Changes

- Modify only what is necessary for the task.
- Make small, testable changes and validate frequently.
- After significant steps, summarize: what changed, what was verified, what remains unresolved.
- If requirements or behavior are ambiguous, stop and ask.
- Preserve existing style, formatting, and architecture.
- Never perform unrelated formatting, cleanup, or aesthetic edits.

## Rule 3 — Read Before Write

- Before writing new code, inspect existing modules, functions, and utilities for reuse.
- Prefer extending existing functionality over duplicating logic.
- Avoid introducing isolated or redundant implementations.
- If proper generalization would require risky architectural changes, ask first.
- Write general reusable function in shared modules, which are called from test script, do not implement substantial logic directly in the test scripts

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
- Never hide command output (`tail`, `head`, background execution, silent redirects, etc.).
- Generate diagnostic plots (`.png`) for visual inspection when useful and report their location.
- Store test outputs, plots, benchmarks, and debug artifacts in organized numbered directories (e.g. `tests/003_case_name/`).

## Rule 5 — Performance

- Keep Python orchestration minimal; avoid Python loops in performance-critical paths.
- Push heavy computation into optimized C/C++/OpenCL/CUDA kernels.
- Prefer flat contiguous arrays and data-oriented memory layouts.
- Preallocate and reuse buffers in performance-critical code.
- Consider cache locality, memory bandwidth, and synchronization costs explicitly.

On GPU:
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

---

### Core Principles
- **Concise, modular code**: small, testable, step-by-step changes; divide-and-conquer
- **Reuse existing functions**: avoid duplication (check first, report if adaptation needed)
- **Pure, data-oriented functions**: explicit inputs/outputs; default named args to avoid long call strings
- **Fail loudly**: no silent handling; assertions for invariants; crashes with stack trace preferred to masking
- **Comment out deprecated/experimental code** instead of deleting; mark with TODO/DEBUG

### Debugging First
- **Debuggability > UX**: do not hide issues
- Initially add debug prints for key variables and flow; remove later after debugging
- Avoid broad try/except as they mask bugs; prefer loud crashes with stack trace
- Make small, testable changes and run after every change
- Log function entries and key conditions when helpful to track flow
- Mark unfinished/experimental code clearly (e.g., # TODO, # DEBUG)

### Performance Guidelines
- Prefer data-oriented code that is cache-friendly and avoids overheads
- Preallocate and reuse buffers; avoid repeated allocation in hot paths
- Be explicit about dtypes/shapes; prefer contiguous memory where possible

### Style Guidelines (Cross-Language)
- Prefer concise/compact code; avoid bloated structures and unnecessary empty lines
- Short variable names OK (math/physics symbols like E, T, m) when locally clear
- Prefer one-liner expressions, assume unlimited line-width
- Avoid line wrapping that hurts readability of expressions
- Inline comments behind the code line for rationale
- **Doxygen**: use `///`; avoid `/* ... */`
- **C++**: use `printf` for debugging over `std::cout`; prefer plain C arrays (double*) in hot paths

### Visualization
- Separate compute vs plotting; no plotting in core algorithms
- Plotting optional via flags (e.g., --noPlot, --saveFig)
- `plt.show()` only in CLI/main, never in library code
- Prefer shared plotting helpers (e.g., plot_utils.py) to avoid duplication


## Repository Navigation

- Read `CODEMAP.md` to understand repository structure.
- `CODEMAP.md` — repository structure and subsystem overview
- `tests/` — primary reference for implemented functionality and usage examples
- Look in `/doc/` for technical details

## Agent Operating Guidelines

These rules apply to any autonomous assistant working inside this repo:

- **No fallbacks or silent fixes**: if required data/geometry is missing, throw a descriptive error immediately. Never synthesize “best-effort” defaults—the failure must be visible.
- **Debug-first mindset**: prioritize traceable crashes, add targeted logging only while investigating.
- **Preserve context**: never remove comments during active debugging; if you must replace logic, comment out the old block instead of deleting so we can revert instantly if needed.
- **Respect existing structure**: extend modules in-place rather than rewriting or rearranging unless explicitly authorized.
- **Use official scripts**: when commands are required, rely on the provided `run.sh`/`make.sh` helpers; never invoke `make` directly.
- **Document parity work**: when mirroring Python ↔ JS features, cite the reference file/function in comments so future maintainers can diff implementations quickly.
- **Never use background commands**: always run commands synchronously to show full output; never use `Background=True`, tail, or similar output-hiding methods.


---

# Core Philosophy

We develop rigorous scientific software. Debuggability, numerical correctness, physical consistency  are paramount.

Code design priorities:

1. **Debuggability:** Code must be transparent, inspectable, and prioritize trace-ability over user experience (never hide issues).
2. **Simplicity:** Clear, clean, direct logic. Elegant design that avoids branching, excessive special conditions, and defensive abstractions.
3. **Performance:** streamlined execution with minimal overhead (e.g. avoid loops in Python), minimal branching, cache aware data-oriented design

## Rule 1 — Fail Loudly

- **No Silent fallbacks:** 
   - Catch-all passes, try-except blocks that mask bugs, and automatic recovery logic are strictly forbidden. 
   - Unexpected states, invalid assumptions, or numerical failures must terminate with explicit errors and full stack traces.
- **Find the root cause** and fix it, never apply "quick-fixes" that hide root causes.

## Rule 2 — Surgical edits & simplicity

- **Surgical Changes:** 
   - Write minimum code necessary to solve the problem. No speculative features, no defensive abstractions.
   - Touch only what is necessary. Never "clean up" or modify adjacent code, comments, or formatting outside your explicit scope.
- **No Guessing:** If you encounter ambiguity or uncertainty, stop and ask for clarification.
- **Strict Checkpointing:** Summarize what was done, what is verified, and what remains after every significant step.
- **backup** 
    - before major chnages of module, make backup copy.
    - commnet-out code with `#` and `//` rather than delete it when trying different version

## Rule 3 — Reusable architecture

- **Inventory First:** 
   - Always thoroughly review the provided reference source-code files to find existing functions, modules, and data structures that can be reused or extended, before writing anything from scratch.
- **Composability Over Bloat:** 
   - We build integrated systems, not isolated scripts. Compose existing, well-tested modules to keep the codebase lean.
   - Separate plotting / diagnostics and other utility modules so they can be reused across different parts of the codebase.
   - Separate GUI, CLI test scripts, and backend modules. All complex logic must be factored into reusable/general functions modules.
   - consolidate testing script, testing scripts should have CLI with params and routing different exection paths called from backend)
- **Generalization Over Duplication:** 
   - If you find an existing function that *almost* fits your needs, try to generalize it rather than duplicating code. 
   - However, do not make major changes that risk breaking backwards compatibility. If this dilemma arises, **report it to me immediately** and ask for approval.

## Rule 4 — Test-Driven Development

- **Physical & Analytical Parity:** Every functionality must have a rigorous validation test. Before coding, define how to prove correctness (e.g., via parity checks against reference code, validation against known analytical solutions for simple cases, physical conservation laws, or symmetry).
- **Verbose Debug Prints:** 
   - Place debug prints of intermediate results throughout the code. 
   - Log function entries and key conditions when helpful to track flow
    -Use a gated `debug_print` function and a verbosity level system from global modules (e.g., `global.py`, `global.h`) for this.
- **Are Numbers Reasonable?:** Strategically place checks throughout the calculation to verify if values are within a reasonable range, and are not `NaN`, infinity, or unexpected zeros.
- **Test on Completion:** 
    - Run tests immediately after code modification. 
    - Never claim code works before you run the tests successfully.
    - Run test on foreground (blocking) with full output, do not use `&` or `| tail`, I need to see full stdout. 
- **Visual review:** Tests in python should produce diagnostic plots using `matplotlib` and save them as `.png`. I will review then to confirm success. Report to me folder/files where to find those .png plots.
- **Structured Test Artifacts:** Assume we will run dozens of tests before a task is complete. Group all debugging, benchmarking, and testing outputs into a dedicated, clearly named, and numbered folder. Do not clutter the root directories. Explicitly report the location of these results.

## Rule 5 — Performance

- Keep Python orchestration minimal; python is slow => avoid loops in Python.
- Prefer flat contiguous arrays and data-oriented memory layouts.
- Preallocate and reuse buffers; avoid repeated allocation in hot paths
- For Low-level perforance critical code (C/C++/OpenCL/CUDA/Compute Shaders)
   - "Gather" operations are typically faster "scatter", prefer "gather" desing where possible.
   - consider memory latency and cache efficiency.
   - Avoid branching, atomics, and unnecessary synchronization in GPU kernels.
   - Use shared/local memory and workgroups efficiently on GPU.

## Rule 6 - Concise style

- Do not create function stubs or abstractions that do too little (e.g., 1–3 lines of code). If it is simple, inline it.
- Do not pass too many arguments over function interfaces; group them into structs/dicts if needed, or utilize globals and class properties.
- Prefer compact code with long lines and minimal whitespace.
- Avoid line wrapping that hurts readability of expressions
- Prefer one-line constructs. Assume infinite line length; do not break function calls onto multiple lines just because of length.
- Prefer short variable names (e.g. math/physics symbols like E_tot, T_ij, m_i)
- **comments**
  - Avoid comments stating the obvious, use it to supplement information (intend, purpose, physics/math derivation, '#TODO'/ '#FIXME')
  - Inline comments behind the code line for rationale
- **Doxygen**: use `///`; avoid `/* ... */`
- **C++**: use `printf` for debugging over `std::cout`; prefer plain C arrays (double*) in hot paths

## Testing compiled modules (C++/Fortran)

- properly settings paths is critical, do not try to guess, follow prescribed protocols and scripts
- test should be in `/test` directory
- for running test with depend on C++/Fotran library `run.sh` inside the test directory which automatically recompile the program, set up paths to inputs and libraries
- Edit `run.sh` scripts if needed
- see skills debug_Fortran and debug_C



