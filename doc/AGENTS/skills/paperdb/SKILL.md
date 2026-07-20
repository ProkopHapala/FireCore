# PaperDB CLI — Skill File for Coding Agents

PaperDB is a scientific paper knowledge base. It indexes PDFs, extracts equations/methods/tags via LLM, and provides search + context-pack assembly for coding agents.

## Quick Start

```bash
# All commands go through python -m paperdb.cli
# Set data directory (default: ~/paperdb/)
export PAPERDB_DATA=/path/to/paperdb_data

# Check what's in the database
python -m paperdb.cli status

# Search for papers
python -m paperdb.cli search "molecular dynamics GPU"
```

## Global Options

All commands accept these global flags **before** the subcommand:

- `--json` — output as JSON (for programmatic consumption)
- `--data-dir PATH` — override data directory (env: `PAPERDB_DATA`)
- `--db-path PATH` — override DB path (env: `PAPERDB_DB`)
- `--llm-config KEY` — LLM config key from `config/LLMs.toml` (env: `PAPERDB_LLM`)

**Always use `--json` when calling from code** — human-readable output uses rich tables that are hard to parse.

## Command Reference

### Ingesting Papers

#### `scan` — Index PDFs from a folder
```bash
python -m paperdb.cli scan /path/to/pdfs --recursive
python -m paperdb.cli --json scan /path/to/pdfs
```
Indexes PDFs by SHA-256 hash, creates paper records, deduplicates. Does NOT extract content — use `ingest` for that.

#### `add` — Add a single paper from path/URL/DOI
```bash
python -m paperdb.cli add /path/to/paper.pdf
python -m paperdb.cli add "https://arxiv.org/pdf/2401.12345"
python -m paperdb.cli add "10.1103/PhysRevE.2024.12345"
```
Fetches metadata from CrossRef (DOI) or arXiv API, downloads PDF if URL/DOI given.

#### `ingest` — Full processing pipeline
```bash
# Ingest a single paper by key
python -m paperdb.cli ingest --paper Macklin_2016_XPBD

# Ingest all papers from a folder (scan + process)
python -m paperdb.cli ingest --folder /path/to/pdfs

# Ingest all unprocessed papers in DB
python -m paperdb.cli ingest --all
```
Runs: PDF→Markdown conversion (Docling), equation extraction, method extraction, summarization (LLM), tag extraction (LLM), search unit building. Skips already-processed papers unless `--force` (via Python API).

#### `sync` — Scan + ingest new/changed papers
```bash
python -m paperdb.cli sync
python -m paperdb.cli sync --folder /path/to/pdfs
```

### Searching & Retrieval

#### `search` — Full-text search with ranking
```bash
# Basic search
python -m paperdb.cli search "constraint-based dynamics"

# With tag filters (prefix with ! to exclude)
python -m paperdb.cli search "collision detection" --tag solver:xpbd --tag domain:game_physics

# Year range filter
python -m paperdb.cli search "force fields" --year 2015-2025

# Show scoring breakdown
python -m paperdb.cli search "Ewald summation" --explain

# JSON output (recommended for agents)
python -m paperdb.cli --json search "Ewald summation" --limit 10
```

**JSON output format** (list of dicts):
```json
[
  {
    "id": 1,
    "paper_key": "Macklin_2016_XPBD",
    "title": "XPBD: Position-Based Simulation of Compliant Constraints",
    "year": 2016,
    "score": 8,
    "breakdown": {"title": 5, "abstract": 2, "fts": 1},
    "matching_units": [],
    "paper": {"id": 1, "paper_key": "...", "title": "...", "doi": "...", ...}
  }
]
```

#### `context` — Assemble a context pack for an LLM
```bash
# Default token budget 24000
python -m paperdb.cli context "GPU-friendly neighbor search algorithms"

# Custom budget and output to file
python -m paperdb.cli context "SPH methods for fluid simulation" --budget 16000 --out context.md

# Include specific content types
python -m paperdb.cli context "rigid body solvers" --include equations,methods
```
Returns a markdown document with: selected paper summaries, relevant equations (with LaTeX), method cards, comparison matrix, and bibliography. This is the **primary output for feeding context to coding agents**.

### Inspecting Papers

#### `inspect` — Full metadata for a paper
```bash
python -m paperdb.cli inspect Macklin_2016_XPBD
python -m paperdb.cli --json inspect Macklin_2016_XPBD
```
Accepts paper key, DOI, or numeric ID.

#### `get` — Get paper content in formats
```bash
python -m paperdb.cli get Macklin_2016_XPBD              # all formats
python -m paperdb.cli get Macklin_2016_XPBD --markdown    # markdown only
python -m paperdb.cli get Macklin_2016_XPBD --json        # structured JSON
python -m paperdb.cli get Macklin_2016_XPBD --bib         # BibTeX
```

#### `equations` — List extracted equations
```bash
python -m paperdb.cli equations Macklin_2016_XPBD
python -m paperdb.cli --json equations Macklin_2016_XPBD
```
Returns equations with `latex_raw`, `latex_normalized`, `equation_number`, `section_path`, `page_number`.

#### `method` — Show method cards for a paper
```bash
python -m paperdb.cli method Macklin_2016_XPBD
python -m paperdb.cli method Macklin_2016_XPBD --name "XPBD solver"
python -m paperdb.cli --json method Macklin_2016_XPBD
```
Returns method cards with `name`, `method_type` (source_algorithm/reconstructed_method), `purpose`, `complexity`, `confidence`, `card_json`.

#### `methods` — Search for methods across papers
```bash
python -m paperdb.cli methods "position-based dynamics" --limit 20
```
Runs a search and returns matching papers.

#### `related` — Find papers sharing tags
```bash
python -m paperdb.cli related Macklin_2016_XPBD --limit 10
```

### Tags & Taxonomy

#### `tags` — List or merge tags
```bash
# List all tags
python -m paperdb.cli tags
python -m paperdb.cli tags --category solver

# Merge two tags (alias → canonical)
python -m paperdb.cli tags --merge xpbd position_based_dynamics
```

### Topical Reviews

#### `topic` — Generate a topical review document
```bash
python -m paperdb.cli topic "molecular force fields"
python -m paperdb.cli topic "GPU collision methods" --out review.md
```
Multi-step: interprets query → searches papers → retrieves method cards → builds comparison matrix → synthesizes review via LLM. Requires `PAPERDB_LLM` to be set.

#### `compare` — Compare methods along axes
```bash
python -m paperdb.cli compare "neighbor search methods" --axes complexity,synchronization,spatial_structure
```

### Export & Migration

#### `export` — Export library as BibTeX
```bash
python -m paperdb.cli export --bibtex --out library.bib
```

#### `migrate` — Import from legacy DB or Mendeley
```bash
python -m paperdb.cli migrate --from /path/to/legacy.db
python -m paperdb.cli migrate --from-mendeley /path/to/library.bib
```

### Status & Re-processing

#### `status` — Database statistics
```bash
python -m paperdb.cli status
```
Returns counts: papers, files, search_units, tags, equations, methods, summaries, runs, topics, context_packs, with_markdown, with_summary, with_doi.

#### `reindex` — Re-run processing operations
```bash
python -m paperdb.cli reindex --re-summarize --re-tag
python -m paperdb.cli reindex --re-extract-equations --llm-config deepseek_chat
```

## Typical Agent Workflows

### Workflow 1: "Find papers about X and get their equations"
```bash
# 1. Search
python -m paperdb.cli --json search "Ewald summation 2D periodicity" --limit 10

# 2. For each result, get equations
python -m paperdb.cli --json equations <paper_key>

# 3. Or get full context pack
python -m paperdb.cli context "Ewald summation 2D periodicity" --out context.md
```

### Workflow 2: "Add a paper and process it"
```bash
# 1. Add from DOI/arXiv/URL/path
python -m paperdb.cli add "10.1103/PhysRevE.2024.12345"

# 2. Ingest (convert + extract + summarize + tag)
python -m paperdb.cli ingest --all

# 3. Verify
python -m paperdb.cli status
python -m paperdb.cli inspect <paper_key>
```

### Workflow 3: "Build a comparison of methods for problem X"
```bash
# 1. Generate topical review
python -m paperdb.cli topic "short-range interaction search on GPU" --out review.md

# 2. Or compare along specific axes
python -m paperdb.cli compare "short-range interaction search" --axes complexity,synchronization,spatial_structure
```

### Workflow 4: "Get everything needed to implement method from paper X"
```bash
# 1. Get method cards
python -m paperdb.cli --json method <paper_key>

# 2. Get equations
python -m paperdb.cli --json equations <paper_key>

# 3. Get full markdown
python -m paperdb.cli get <paper_key> --markdown

# 4. Or assemble a targeted context pack
python -m paperdb.cli context "<paper_key> method equations" --include equations,methods --out context.md
```

## Key Concepts

- **Paper key**: Human-readable identifier like `Macklin_2016_XPBD`. Used in all CLI commands. Also accepts DOI or numeric ID.
- **Search units**: FTS-indexed chunks at section/paragraph/equation/method level — not whole-paper. Enables precise retrieval.
- **Processing runs**: Track which operations ran on which paper, with which backend/config. Enables skip-if-equivalent logic.
- **Method cards**: Two types — `source_algorithm` (extracted from paper text) and `reconstructed_method` (LLM-reconstructed with assumptions, steps, I/O).
- **Context pack**: Assembled markdown with selected papers, equations, methods, comparison matrix, bibliography. The primary artifact for feeding to coding agents.
- **Tags**: Canonical names with categories (domain, solver, method, physical_system, etc.). Aliases map variant names to canonical.

## Environment Variables

| Variable | Default | Description |
|---|---|---|
| `PAPERDB_DATA` | `~/paperdb/` | Data directory (DB, papers, logs) |
| `PAPERDB_DB` | `$PAPERDB_DATA/papers.db` | SQLite database path |
| `PAPERDB_LLM` | first in LLMs.toml | LLM config key for extraction/synthesis |

## Python API (for programmatic access)

All CLI commands delegate to `paperdb.PaperDB`:
```python
from paperdb import PaperDB
db = PaperDB()  # uses env vars or defaults

results = db.search("molecular dynamics", limit=10)
for r in results:
    print(r['paper_key'], r['score'], r['title'])
    eqs = db.get_equations(r['paper_key'])
    methods = db.get_methods(r['paper_key'])
```
