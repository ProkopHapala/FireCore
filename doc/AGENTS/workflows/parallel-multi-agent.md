---
description: Distribute a software project across multiple parallel LLM coding agents without integration failures
---

# Parallel Multi-Agent Development Workflow

How to split a software project across multiple LLM coding agents working in
parallel, avoiding the interface mismatches and integration failures that
plague naive parallelization.

Based on experience from the PaperDB project (6 parallel agents, 13 interface
mismatches, documented in `docs/tasks/paperdb/task7_integration_gaps.md` and
`docs/topical_audit/paper_db_notes.md` §20).

---

## When to use this workflow

- A project is too large for one agent to complete in a single session
- The project can be decomposed into modules with clear boundaries
- You have multiple agent sessions available (Cursor, Windsurf, Devin, etc.)
- You want to parallelize development to save wall-clock time

## When NOT to use

- Project is small enough for one agent (always prefer single-agent if feasible)
- Modules are tightly coupled with no clean boundaries
- The interface between modules is not yet understood (design first, then split)

---

## Step-by-step procedure

### Step 1: Design the architecture first

Before splitting work, produce a design document covering:

- **Module decomposition**: what modules exist, what each owns
- **Data flow**: how modules call each other
- **Shared data model**: what types/structures cross module boundaries
- **Storage layer**: database schema, file formats, APIs

Write this down in a single markdown file. All agents will reference it.

```
docs/design/<project>_notes.md
```

### Step 2: Identify the critical-path task

Every parallel project has one task that **all others depend on**. Usually
this is the foundation layer: database schema, core models, repository/API
interface.

**This task MUST complete before any other agent starts coding.**

Identify it explicitly:
- Which task defines the interfaces others call?
- Which task creates the shared types/models?
- Which task sets up the project skeleton (pyproject.toml, package structure)?

### Step 3: Write the interface contract file

This is the **most important** deliverable. Before any parallel work starts,
produce a single Python file that defines:

1. **All shared types** (Pydantic models, dataclasses, TypedDicts)
2. **All interface protocols** (`typing.Protocol` or ABC) with exact method
   signatures — parameter names, types, return types
3. **All shared constants** (enums, config keys, etc.)

```python
# <project>/interfaces.py — THE contract
from typing import Protocol, Optional
from dataclasses import dataclass

# Shared types — all agents import these, never redefine
@dataclass
class Paper:
    paper_key: str
    doi: Optional[str] = None
    title: Optional[str] = None
    # ...

# Interface protocol — exact method signatures
class RepositoryProtocol(Protocol):
    def upsert_paper(self, paper_key: str, doi: str = None,
                     title: str = None, ...) -> int: ...
    def get_paper(self, paper_id: int) -> Optional[Paper]: ...
    # ... every method, with exact parameter names and types
```

**Rules for the contract file:**
- Owned by the architect only (not by any parallel agent)
- All agents `from <project>.interfaces import ...`
- If an agent needs a method not in the contract, they **request it** — never invent
- The contract file uses **keyword arguments**, not model-object parameters
  (downstream callers universally prefer kwargs; this was a major source of
  disagreement in PaperDB)

### Step 4: Create the method-name vocabulary table

In the task README, maintain a table listing every cross-module method:

| Method | Owner task | Caller tasks | Signature |
|--------|-----------|-------------|-----------|
| `upsert_paper(paper_key, doi?, title?, ...)` | Task 1 | Task 2, 5 | `-> int` |
| `add_paper_file(paper_id, path, sha256?, ...)` | Task 1 | Task 2, 5 | `-> int` |
| `search(query, limit=20)` | Task 3 | Task 4 | `-> list[SearchResult]` |

Every agent references this table. If a method is missing, the agent raises
the gap immediately rather than inventing a name.

### Step 5: Define file ownership map

Each file is owned by **exactly one task**. No agent modifies another agent's
files. Document this explicitly:

| File | Owner |
|------|-------|
| `<project>/interfaces.py` | Architect |
| `<project>/__init__.py` | Task 1 |
| `<project>/db/repository.py` | Task 1 |
| `<project>/search/fts.py` | Task 3 |
| ... | ... |

If an agent needs a change in another agent's file, they document it as a
**dependency note** and wait — they do not modify the file.

### Step 6: Complete the critical-path task FIRST

```
Phase 1: Critical-path task (foundation) — COMPLETE
    ↓ interfaces.py exists, repository.py exists, models exist
Phase 2: All parallel tasks start, importing interfaces.py
    ↓ each task codes against the real contract
Phase 3: Integration owner wires facade, runs end-to-end tests
```

**Do NOT skip Phase 1.** The cost of waiting is hours. The cost of 13
interface mismatches across 6 modules is days of debugging.

### Step 7: Assign an explicit integration owner

One agent (or the architect) is responsible for:

- Wiring the main API facade to all submodule implementations
- Running end-to-end integration tests
- Resolving interface mismatches discovered during integration
- Verifying that `NotImplementedError` stubs are replaced with real delegation

**This is not a shared responsibility.** Without a single owner, facade
wiring falls through the cracks (as happened in PaperDB).

### Step 8: Smoke test gate for each task

Before a task is marked "done", it must pass:

```bash
# 1. Import check — does the module import without errors?
python -c "from <project>.ingest.scanner import scan_folder"

# 2. Contract compliance — does the module call the interface protocol correctly?
#    Use a mock that implements the protocol and asserts calls.
python -m pytest tests/<project>/test_integration/mock_test.py

# 3. No NotImplementedError — if the task owns facade methods, they must be wired.
python -c "from <project> import ProjectAPI; api = ProjectAPI(); api.search('test')"
```

### Step 9: End-to-end integration test (Phase 3)

After all tasks are complete, the integration owner writes and runs:

```bash
# Full workflow test: create → process → search → retrieve
python -m pytest tests/<project>/test_integration/test_e2e.py
```

This test exercises the full system end-to-end, catching any remaining
mismatches that unit tests with mocks missed.

---

## Anti-patterns to avoid

### 1. No interface contract — "just code against the spec"

**What happens**: Each agent imagines a different API. 13 mismatches at
integration time. Days of debugging.

**Fix**: Always produce `interfaces.py` before parallel work starts.

### 2. Model-object parameters instead of keyword arguments

**What happens**: Foundation task implements `repo.upsert_paper(paper: Paper)`.
Downstream tasks call `repo.upsert_paper(paper_key=..., doi=...)`. Crash.

**Fix**: Use keyword arguments in repository/service APIs. Reserve model
objects for return types and internal logic.

### 3. `hasattr` guards for missing methods

**What happens**: `if hasattr(repo, 'touch_file'): repo.touch_file(...)` —
silently skips functionality when method doesn't exist. Bugs hidden.

**Fix**: If the method is in the contract, call it directly. Let it crash.
If optional, define it in the protocol with a default no-op.

### 4. All agents start simultaneously

**What happens**: Downstream agents code against imagined APIs. Foundation
task's actual implementation doesn't match anyone's assumptions.

**Fix**: Stagger. Foundation task completes first. Then parallel tasks start
with the real contract in hand.

### 5. Facade stubs left as `NotImplementedError`

**What happens**: Foundation task leaves stubs "for other agents to wire".
No agent does the wiring. System doesn't work end-to-end.

**Fix**: Integration owner is explicitly responsible for wiring all stubs.
This is a deliverable, not an afterthought.

### 6. Method names invented in isolation

**What happens**: Task 6 calls `repo.add_method()`, but Task 1 implemented
`repo.upsert_method()`. Runtime `AttributeError`.

**Fix**: Method-name vocabulary table in the task README. Every cross-module
method is listed with exact name and signature before coding starts.

### 7. No end-to-end test

**What happens**: Each task passes its unit tests with mock repositories.
Integration fails because mocks don't match the real Repository.

**Fix**: Integration owner writes end-to-end tests as a separate deliverable.

---

## Task decomposition heuristics

### How to split work across agents

**Good split criteria:**
- **By layer**: foundation/db → identity/ingest → search → CLI/MCP → synthesis
- **By ownership boundary**: each agent owns a directory or set of files
- **By dependency depth**: tasks with no dependencies start first

**Bad split criteria:**
- **By feature** (e.g. "agent 1 does search, agent 2 does tagging") — these
  often share models and repository methods, leading to conflicts
- **By file count** — doesn't respect logical boundaries
- **Random split** — maximizes coupling

### How many agents?

- **2-3 agents**: manageable, minimal coordination overhead
- **4-6 agents**: works well with a contract file + integration owner
- **7+ agents**: coordination overhead exceeds parallelization benefit;
  consider splitting into phases instead

### Dependency graph

Draw the dependency graph before assigning tasks:

```
Task 1 (Foundation) ──┬──> Task 2 (Identity)
                      ├──> Task 3 (Search)
                      ├──> Task 4 (CLI/MCP)
                      ├──> Task 5 (Extraction)
                      └──> Task 6 (Synthesis)

Phase 1: Task 1
Phase 2: Tasks 2-6 (parallel)
Phase 3: Integration owner
```

Tasks with no dependencies between them can run in parallel. Tasks with
dependencies must be staggered.

---

## Checklist (copy this into your task README)

### Before parallel work starts
- [ ] Architecture design document written
- [ ] Interface contract file (`interfaces.py`) created with all types + protocols
- [ ] Method-name vocabulary table in task README
- [ ] File ownership map defined
- [ ] Critical-path task identified
- [ ] Integration owner assigned

### Critical-path task (Phase 1)
- [ ] Foundation task completed: repository, models, package skeleton
- [ ] `interfaces.py` importable by all downstream tasks
- [ ] Smoke test: `python -c "from <project> import ProjectAPI"` works

### Each parallel task (Phase 2)
- [ ] Imports from `interfaces.py`, not from imagined APIs
- [ ] Uses keyword-argument Repository methods (not model-object params)
- [ ] Does not modify files owned by other tasks
- [ ] Does not use `hasattr` guards for contract methods
- [ ] Passes import check: `python -c "from <project>.<module> import <func>"`
- [ ] Passes mock-repo test: calls match the protocol signatures

### Integration (Phase 3)
- [ ] Facade wired: no `NotImplementedError` stubs remain
- [ ] End-to-end test written and passing
- [ ] All interface mismatches resolved
- [ ] `python -c "from <project> import ProjectAPI; api = ProjectAPI(); api.status()"` works
