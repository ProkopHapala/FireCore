---
name: doc-task-summary
description: After implementing — update file headers, README, topical audits; note duplication
---

## After You Implement — Update Docs

Same session as the code, while context is fresh. Style rules: `doc-audit` skill.

### 1. File / module headers (changed files only)

At top of each new or materially changed source file:

- **Essence** — one sentence: what problem this file solves (not a file list).
- **Design** — only non-obvious choices (SSOT, data flow, why this approach).
- **Open issues / caveats** — bullets for landmines: partial/deferred state, ordering deps, duplicate parallel implementations, strict/fail-loud APIs, import/export limits, out-of-scope TODOs.
- **Functions** — one-line purpose comment only where the name is not enough (language style per `doc-audit`: Python docstrings, OpenCL `//`).

Do not delete existing comments. Do not document every parameter. Do not refactor while documenting.

### 2. README.md in the folder you touched

- If missing: create one (1–3 sentences + bullet list of key files, see `doc-audit` skill for format)
- If exists: add/update entries for new or changed files
- Keep it a quick index, not a manual — for user-facing features add a **Tutorial** snippet and **fit knobs** table when non-obvious (see `spammm/surfaces/README.md` contact-surface section)

### 3. Topical audit (`doc/topical_audit.md`)

- New take on an existing topic: add a row to the implementations table
- Topic section doesn't exist yet: add one (see `doc-audit` skill for format)
- Mark old implementations as `deprecated` if superseded
- Cross-cutting caveats from headers → audit **Open Issues** when relevant
- Link to `doc/Takeways.md` when debugging pitfalls are non-obvious (sign conventions, buffer layouts, coordinate frames)

### 4. Duplication check

- Search for similar logic across languages — if found, note in topical audit or file caveats
- Don't consolidate now unless trivial — just record it

### 5. Full topical doc (when feature is multi-file or user asks)

Use `doc-audit` OKF format — example: `doc/Topics/AFM/ContactSurface_Static.md`:

| Section | Content |
|---------|---------|
| YAML frontmatter | `type`, `title`, `tags`, `timestamp` |
| Summary | One paragraph + optional mermaid data-flow |
| Tutorial | Library/CLI snippet, env flags, review artifacts |
| API reference | Tables: AFMulator methods, classes, kernels |
| Parity status | Metrics + test file pointers |
| Background | Physics, basis math, implementation plan (existing sections) |
| Pitfalls | Link to `doc/Takeways.md`, not duplicate full debug narratives |

Also update `CODEMAP.md` (file one-liners, doc paths, test scripts).

## Don't Over-Do It

- **Small edit**: caveat in header if behavior changed; README bullet if user-facing
- **New feature**: headers + caveats on new files; README update; audit row if new topic
- **Full documentation**: only when explicitly asked — use `doc-audit` + topical doc + Takeways for pitfalls

## Related Skills

- `doc-audit` — format and inline-comment conventions per language
- `doc-read-navigate` — read existing headers (caveats) before extending code
- `doc/Takeways.md` — cross-cutting pitfalls to link from topical docs and audits
