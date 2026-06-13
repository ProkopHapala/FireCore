#!/usr/bin/env python3
"""
Classify all .md files in the repo into categories:
  llm_chat     - LLM conversation transcripts (USER/Gemini/ChatGPT/Claude/etc. headings)
  dev_notes    - Debug logs, session notes, dated progress entries, parity discussions
  report       - Implementation reports, findings, status summaries
  documentation- Proper reference docs (tutorials, API docs, design docs)
"""

import os, re, sys
from pathlib import Path

REPO = "/home/prokop/git/FireCore"

# LLM model name patterns (as markdown headings level 1-3)
LLM_HEADING_RE = re.compile(
    r'^#{1,3}\s*('  # heading level 1-3
    r'USER'
    r'|Gemini'
    r'|ChatGPT[^\w]'
    r'|Claude[^\w]'
    r'|DeepSeek[^\w]'
    r'|Kimi[^\w]'
    r'|SWE-[^\w]'
    r'|GPT-[^\w]'
    r'|Devin[^\w]'
    r'|Copilot[^\w]'
    r'|Post-mortem'
    r'|Cody[^\w]'
    r'|Tabby[^\w]'
    r'|OpenAI[^\w]'
    r'|gpt-[^\w]'
    r'|o1[^\w]'
    r'|o3[^\w]'
    r'|llama[^\w]'
    r'|mistral[^\w]'
    r'|codestral[^\w]'
    r'|perplexity[^\w]'
    r'|cursor[^\w]'
    r'|sweep-[^\w]'
    r')',
    re.MULTILINE | re.IGNORECASE
)

# Share links for LLM chat transcripts
SHARE_LINK_RE = re.compile(
    r'https://(chatgpt\.com/share|gemini\.google\.com/share|claude\.ai/share|chat\.deepseek\.com/share|www\.kimi\.com/share)',
    re.IGNORECASE
)

# Date patterns for dev notes / debug logs
# YYYY-MM-DD, YYYY/MM/DD, or Month DD, YYYY
DATE_RE = re.compile(
    r'(?:^#{1,4}\s*(\d{4}[-/]\d{2}[-/]\d{2})'  # heading with date
    r'|\b(\d{4}[-/]\d{2}[-/]\d{2})\s+(session notes|parity recap|progress update|debugging|testing|check|meeting|review|log)'
    r'|\b(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)[a-z]*\s+\d{1,2},?\s+\d{4}\b'
    r'|\b\d{4}[-/]\d{2}[-/]\d{2}\b.*(?:notes|recap|update|status|log|debug|fix|parity|test))',
    re.MULTILINE | re.IGNORECASE
)

# Keywords/phrases typical of dev notes / debug logs
DEV_NOTE_KEYWORDS = [
    "session notes", "parity recap", "progress update", "debugging", "root cause",
    "problems encountered", "fixes applied", "solutions applied", "next steps",
    "current status", "parity status", "test outcome", "key findings",
    "bottlenecks", "inefficiencies", "crashes", "workaround", "regression",
    "memory leak", "buffer overrun", "race condition", "segfault", "nan",
    "mismatch", "discrepancy", "magnitude mismatch", "scaling issue",
    "double-free", "out-of-bounds", "undefined behavior", "hot path",
    "TODO", "FIXME", "HACK", "WIP", "DEBUG", "WARN",
]

# Keywords for reports (more polished than dev notes)
REPORT_KEYWORDS = [
    "implementation report", "design document", "system architecture",
    "overview and physics", "algorithmic strategy", "what we did",
    "what we achieved", "validation results", "key takeaways",
    "reproducible runs", "usage", "calibrated parameters", "output directory",
    "formal fitting", "consolidated verification", "output:",
]

LLM_KEYWORDS = [
    "# USER\n", "# Gemini\n", "# ChatGPT", "# Claude", "# DeepSeek",
    "# Kimi", "# SWE-", "# GPT-", "# Copilot", "# Devin",
    "# Post-mortem", "# gpt-", "# o1", "# o3",
]


def score_content(text: str, filename: str) -> dict:
    """Score a markdown file for each category. Higher = stronger evidence."""
    scores = {"llm_chat": 0, "dev_notes": 0, "report": 0, "documentation": 0}
    lower = text.lower()

    # --- LLM Chat scoring ---
    if SHARE_LINK_RE.search(text):
        scores["llm_chat"] += 10

    llm_headings = LLM_HEADING_RE.findall(text)
    if llm_headings:
        scores["llm_chat"] += len(llm_headings) * 5

    # Count conversation-style turns (# USER + # <model>)
    user_turns = len(re.findall(r'^#{1,3}\s*USER\s*$', text, re.MULTILINE | re.IGNORECASE))
    if user_turns > 0:
        scores["llm_chat"] += user_turns * 3

    # --- Dev Notes scoring ---
    date_matches = len(DATE_RE.findall(text))
    if date_matches > 0:
        scores["dev_notes"] += date_matches * 3

    for kw in DEV_NOTE_KEYWORDS:
        if kw.lower() in lower:
            scores["dev_notes"] += 1

    # Tables with "Status" or "Result" columns are common in dev notes
    if re.search(r'\|\s*(Status|Result|PASS|FAIL)\s*\|', text):
        scores["dev_notes"] += 3

    # --- Report scoring ---
    for kw in REPORT_KEYWORDS:
        if kw.lower() in lower:
            scores["report"] += 1

    # --- Documentation scoring ---
    # Clean, structured docs tend to have many ## sections and few LLM names
    section_count = len(re.findall(r'^#{2,3}\s', text, re.MULTILINE))
    if section_count > 3:
        scores["documentation"] += section_count

    # Penalize documentation score if LLM markers present
    if scores["llm_chat"] > 0:
        scores["documentation"] -= scores["llm_chat"] * 2

    # Filename heuristics
    fname_lower = filename.lower()
    if "_discussion" in fname_lower or "discussion_" in fname_lower:
        scores["llm_chat"] += 5
    if "_report" in fname_lower or "report_" in fname_lower:
        scores["report"] += 10
    if fname_lower.startswith("notes_") or "_notes" in fname_lower:
        scores["dev_notes"] += 5
    if "_plan" in fname_lower or "todo" in fname_lower:
        scores["dev_notes"] += 3

    # Paths strongly suggest category
    if "/devnotes/" in fname_lower or "/dev_notes/" in fname_lower:
        scores["dev_notes"] += 2
    if "/topics/" in fname_lower:
        scores["report"] += 1  # topics are often reports/discussions
    if "/markdown/" in fname_lower or "/doc/" in fname_lower:
        scores["documentation"] += 1

    return scores


def is_config(filepath: str) -> bool:
    """Agent/IDE config files (skills, workflows, AGENTS, etc.)"""
    f = filepath.lower()
    if "/.windsurf/" in f or "/.cursor/" in f or "/.devin/" in f or "/.aider/" in f:
        return True
    if "/skills/" in f and f.endswith("skill.md"):
        return True
    if "/workflows/" in f and (f.endswith(".md") or f.endswith("_workflow")):
        return True
    basename = os.path.basename(f)
    if basename in ("agents.md", "agents_alternative.md"):
        return True
    return False


WIP_MARKERS = [
    "TODO", "FIXME", "HACK", "WIP", "DRAFT",
    "TEMPORARY", "PLACEHOLDER", "PENDING",
    "UNFINISHED", "INCOMPLETE", "ROUGH", "SKETCH",
    "FUTURE WORK", "TO BE DONE",
]


def doc_quality(filepath: str, text: str) -> str:
    """For files classified as 'documentation', determine if polished reference or draft/WIP."""
    f = filepath.lower()
    fname = os.path.basename(f)

    # README / MIGRATION_GUIDE / TUTORIAL / MANIFEST files are almost always reference docs
    ref_names = ["readme.md", "migration_guide.md", "modular_system_summary.md", "index.md"]
    if fname in ref_names:
        return "reference_doc"

    # Explicit WIP filename signals
    wip_patterns = ["_plan", "_draft", "_research", "_ideas", "_wip", "_todo", "_experimental", "_sketch"]
    for pat in wip_patterns:
        if pat in fname:
            return "draft_notes"

    # Dev-notes paths are inherently WIP-oriented
    if "/dev_notes/" in f or "/devnotes/" in f:
        return "draft_notes"

    # Count strong WIP markers in text
    lower = text.lower()
    marker_count = sum(1 for m in WIP_MARKERS if m.lower() in lower)

    # Count sections
    section_count = len(re.findall(r'^#{2,3}\s', text, re.MULTILINE))

    # Many sections with few markers = polished reference doc
    if section_count >= 8 and marker_count <= 2:
        return "reference_doc"

    # Moderate sections, no markers = reference
    if section_count >= 4 and marker_count == 0:
        return "reference_doc"

    # Many markers regardless of structure = draft
    if marker_count >= 3:
        return "draft_notes"

    # Some markers but poor structure = draft
    if marker_count > 0 and section_count < 4:
        return "draft_notes"

    # Very short and few sections = draft
    if section_count < 2 and len(text) < 3000:
        return "draft_notes"

    # Default: moderate quality goes to reference
    return "reference_doc"


def classify_file(filepath: str, text: str) -> str:
    scores = score_content(text, filepath)
    # Determine winner
    best = max(scores, key=scores.get)
    best_score = scores[best]

    # Special case: if llm_chat score >= 10, it's almost certainly an LLM chat
    if scores["llm_chat"] >= 10:
        return "llm_chat"

    # If dev_notes and report are close, and file is in Topics, prefer report for polished ones
    if best == "dev_notes" and scores["report"] >= scores["dev_notes"] - 2:
        # Check if it looks like a polished report (has Usage, outputs, etc.)
        if "usage" in text.lower() and "output" in text.lower():
            return "report"

    # Default threshold: needs some positive score to not be "documentation"
    if best_score <= 0:
        return "documentation"

    # If report wins but score is low, and no strong dev markers, call it documentation
    if best == "report" and best_score < 3 and scores["dev_notes"] < 2:
        return "documentation"

    return best


def main():
    # Read the list of all .md files
    with open("/tmp/all_md_files.txt") as f:
        paths = [p.strip() for p in f if p.strip()]

    categories = {
        "llm_chat": [],
        "dev_notes": [],
        "report": [],
        "config": [],
        "reference_doc": [],
        "draft_notes": [],
    }

    for p in paths:
        try:
            with open(p, "r", encoding="utf-8", errors="ignore") as f:
                text = f.read(15_000)
        except Exception as e:
            print(f"WARN: cannot read {p}: {e}", file=sys.stderr)
            continue

        # 1. Config filter (highest priority)
        if is_config(p):
            cat = "config"
        else:
            cat = classify_file(p, text)
            # 2. Split broad documentation into polished reference vs draft/WIP
            if cat == "documentation":
                cat = doc_quality(p, text)

        categories[cat].append(p)

    # Sort within each category
    for cat in categories:
        categories[cat].sort()

    # Print summary to stdout
    print("# Markdown File Classification")
    print("")
    total = sum(len(v) for v in categories.values())
    print(f"Total files: {total}")
    print("")

    for cat, files in categories.items():
        print(f"## {cat} ({len(files)})")
        for f in files:
            rel = os.path.relpath(f, REPO)
            print(f"- {rel}")
        print("")

    # Write structured markdown output
    out_path = os.path.join(REPO, "doc", "MarkdownFileClassification.md")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("# Markdown File Classification\n\n")
        f.write("This document lists all `.md` files in the repository and classifies them by type.\n")
        f.write("Generated automatically by `classify_md.py`.\n\n")
        f.write("## Categories\n\n")
        f.write("- **llm_chat** — Conversation transcripts with LLMs (USER/Gemini/ChatGPT/Claude/DeepSeek/Kimi/GPT/SWE headings, share links)\n")
        f.write("- **dev_notes** — Debug logs, session notes, dated progress entries, parity discussions, work-in-progress notes\n")
        f.write("- **report** — Polished implementation reports, findings, status summaries with usage instructions\n")
        f.write("- **config** — Agent/IDE configuration files (skills, workflows, AGENTS, .windsurf, .cursor, .devin)\n")
        f.write("- **reference_doc** — Polished reference documentation (tutorials, API docs, design docs, READMEs, user guides)\n")
        f.write("- **draft_notes** — Temporary / exploratory / WIP notes, plans, research sketches, unfinished docs\n\n")

        for cat, files in categories.items():
            f.write(f"## {cat} ({len(files)})\n\n")
            for path in files:
                rel = os.path.relpath(path, REPO)
                f.write(f"- `{rel}`\n")
            f.write("\n")

    print(f"\nWrote classification to: {out_path}")


if __name__ == "__main__":
    main()
