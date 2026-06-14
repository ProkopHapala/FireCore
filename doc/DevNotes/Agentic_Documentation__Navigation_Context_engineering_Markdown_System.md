# USER

This is quation about context engineering for agentic coding (vibecoding)

I'm insterested to make good docummentation/navigation system for agetnts to my large repo furecore, especially in branch prokop where I do many experimentation and it is like scratchpad of playgroun 
https://github.com/ProkopHapala/FireCore/tree/prokop
and this 
https://github.com/ProkopHapala/SimpleSimulationEngine

I want some general tranferable system which I can use with any coding agent (Cursor, Devin, Antigravity, opencode, codex, claude code, pi etc)

The standartized system is only AGENTS.md ?

I found this system which seems quite simple and robust, and it is actually siilar to what I was already doing myself
https://github.com/agent0ai/dox/blob/main/AGENTS.md

but I'm not sure if it is best to add all documentation in AGNETS.md or if there can be different type of files?

are there other similar alternative system developed for this purpose?

What I do not like that they all have the same name so it is difficult to navigate from vscode ctrl+P but I can deal with that. 

Also I do not want to overhem agents with too much context

---

# Gemini

Managing the context window for coding agents during rapid, scratchpad experimentation (like your `prokop` branch in `FireCore` or your work in `SimpleSimulationEngine`) requires a careful balance. Overloading an agent causes context drift and degrades reasoning, while providing too little information leads to hallucinations and broken environments.

---

### 1. Is `AGENTS.md` the Only Standard?

While `AGENTS.md` is emerging as the cross-tool, open standard (backed by the Agentic AI Foundation to prevent tool sprawl), it is part of a broader ecosystem of context files.

Different ecosystems and platforms rely on specific configurations:

| Tool / Ecosystem | Primary Context Surface File(s) |
| --- | --- |
| **Cross-Tool Standard** | `AGENTS.md` (Supported natively or via fallback rules by Cursor, Claude Code, OpenCode, Aider, etc.) |
| **Claude Code (Terminal)** | `CLAUDE.md` (Uses `/init`. You can point it directly to your global file with an `@AGENTS.md` import line). |
| **Cursor / Windsurf** | `.cursorrules` / `.windsurfrules` (Supports YAML frontmatter for file-pattern targeting). |
| **GitHub / Copilot Workspace** | `.github/copilot-instructions.md` |

### 2. Solving the VS Code Navigation Dilemma (`Ctrl+P`)

It can be frustrating when multiple context files across different projects are all named `AGENTS.md`, making fuzzy-search navigation a chore. To fix this without breaking your AI toolchain, **leverage symlinks, localized mapping, or IDE setting exclusions.**

* **The Symlink Hack:** Name your actual file something distinct and easy to type, like `FIRECORE_CONTEXT.md` or `SIM_ENGINE_AGENT.md`. Then, create a symlink at the root level of your project:
```bash
ln -s FIRECORE_CONTEXT.md AGENTS.md

```

Coding agents will follow the `AGENTS.md` symlink and read the text perfectly, while you can hit `Ctrl+P` and type `FIRECORE` to jump straight to your documentation.
* **The VS Code Search Filter:** Alternatively, keep the name `AGENTS.md` but use directory-based search typing. Press `Ctrl+P` and type `fire/age` or `simp/age` to immediately isolate the correct file by its parent directory path.

### 3. Preventing Context Overwhelm: The "Split Layout" Pattern

Putting every piece of project documentation into a single root file will quickly degrade the performance of an agent. Instead, organize your context layers by combining a core `AGENTS.md` file with a modular, **Split Layout Architecture**.

```
Your Repository Root
├── AGENTS.md              <-- Minimalist, high-priority operational layer
├── .agent/                <-- Dedicated subdirectory for deep reference
│   ├── ARCHITECTURE.md    <-- Mental models, simulation step flow, math equations
│   ├── CONVENTIONS.md     <-- Low-level pointer/memory math or naming style rules
│   └── MEMORY.md          <-- Dynamic scratchpad of current branch goals

```

#### Layer 1: The Operational File (`AGENTS.md`)

Keep this file concise (ideally under 300 lines). It should act as an operational blueprint rather than an academic reference manual, containing text an agent cannot infer from reading the source code directly.

* **Project Core & Constraints:** 2–3 sentences defining what the tool does (e.g., *"This is a low-level physical engine utilizing customized relaxed scans to avoid brute-force surface potentials"*).
* **Toolchain Commands:** Exact execution flags for compiling, building, and running tests (e.g., `cargo test -- --nocapture` or custom Python/C++ build commands).
* **Hard Boundaries:** Explicit "Never/Always" instructions.
* *Example:* `NEVER modify the raw kernel compilation logic without authorization.`
* *Example:* `ALWAYS handle raw memory allocation boundaries explicitly; do not wrap methods in generic try/catch blocks if they use custom pointer results.`

#### Layer 2: The Reference Directory (`.agent/`)

Instead of pasting full code snippets or massive architecture notes into the main file, point the agent to specialized supplementary materials only when required:

```markdown
## Deep Reference Maps
- For GPU Workgroup optimization guidelines and memory layout, see `.agent/CONVENTIONS.md`.
- For the active state of current experiments and feature flags in the prokop branch, see `.agent/MEMORY.md`.

```

Modern agents are capable of identifying these file paths in your root brief and will fetch the secondary files autonomously only if a task demands it.

### 4. Tailoring Context for Rapid Experimentation (`prokop` branch)

Because your `prokop` branch acts as a fluid scratchpad and playground for physical simulations, your context strategy should account for changing codebases.

* **Maintain a `MEMORY.md` file:** When working in an active playground branch, keep a running file inside your `.agent/` folder that lists your current short-term goals. Updated manually or via a quick prompt at the end of a coding session, it ensures continuity across agent lifetimes:
```markdown
### Current Active Context
- We are currently testing WebGPU workgroup sync issues inside the simulation loop.
- Status: The baseline integration works, but the memory bound checks are throwing edge-case errors.
- Next Task: Optimize the trajectory interpolation kernel.

```

* **Isolate Code Experiments:** Avoid writing long prose describing how your modules link together. Agents navigate file structures efficiently on their own. Instead, use your context files to highlight non-standard architecture choices, keeping your playground organized and your agents aligned with your current goals.

---

# DeepSeek-V4-flash

I think you've touched on something really interesting here. The question isn't just about *one* file, but about a whole *system* for navigating your repo, especially when you're in an experimental branch like `prokop`, which is a dynamic, playful space.

The ecosystem is moving beyond a single `AGENTS.md` file to more sophisticated, multi-file systems. The good news is that many of these new approaches solve exactly the pain points you're describing.

### 📄 AGENTS.md: The Universal Standard (Your Starting Point)

Yes, `AGENTS.md` is the widely adopted, open standard for this purpose. It's essentially a "README for AI agents" written in plain Markdown, placed at your repository's root.

*   **What to Put In**: Focus on the essential, durable information an agent needs to work effectively:
    *   **Build & Test Commands**: The exact, copy-pasteable commands to build, test, and lint the project.
    *   **Project Architecture**: A quick map of the core modules and how they interact.
    *   **Non-Obvious Rules**: Any specific constraints, boundaries (files to never touch), or code style rules that differ from the default.
    *   **Pointers**: Instead of copying everything, you can link to more detailed documentation.

A 2025 Princeton study found that a well-written `AGENTS.md` can reduce agent runtime by **28.6%** and token usage by **16.6%**. Crucially, human-written files are much more effective than auto-generated ones.

**The Multi-File "Tree" Strategy:**
The standard is powerful because it's hierarchical. Your repo can have a root `AGENTS.md` and then other `AGENTS.md` files in sub-directories. The coding agent will read the closest one to the file it's working on, allowing for tailored instructions for specific parts of your project.
*   This directly addresses your desire to keep instructions near the experimental code in your `prokop` branch and your `SimpleSimulationEngine` repo, preventing the root file from getting bloated.

### 🤔 Other Formats: The "Pre-Standardization" World

Because `AGENTS.md` is a relatively new standard, you'll still see other tool-specific files in the wild. They all serve a similar purpose, but `AGENTS.md` is the one designed for cross-agent compatibility.

| Tool | Looks For | Why It Exists |
| :--- | :--- | :--- |
| **Claude Code** | `CLAUDE.md` | Anthropic's native format, though the `CLAUDE.md` spec is aligned with `AGENTS.md` to be compatible going forward. |
| **Cursor** | `.cursorrules` | Cursor's original method, now also supports `AGENTS.md`. |
| **GitHub Copilot** | `.github/copilot-instructions.md` | Native format for GitHub Copilot. |
| **Gemini CLI** | `GEMINI.md` | Native format for Gemini CLI. |

While you *could* use these, using `AGENTS.md` is future-proof. It's a standard stewarded by the Linux Foundation's Agentic AI Foundation and is supported across many tools including Codex CLI, Cursor, Aider, and OpenHands.

### 🤯 Beyond the Single File: Modern Multi-File Systems

This is where the most interesting solutions for your "large repo" and "navigation" needs lie. These are systems that build on top of the `AGENTS.md` idea.

*   **The DOX System (You Found It)**: Exactly as you described, it's a lightweight, self-documenting pattern where an agent builds a hierarchy of `AGENTS.md` files as it works. It starts with a root `AGENTS.md` that tells the agent to scan the repo and create child `AGENTS.md` files for different areas. Its simplicity is its strength.

*   **`agents-md` (Composable Fragments)**: A more automated approach. You write your instructions across many small Markdown files (e.g., `docs/overview.agents.md`, `api/endpoints.agents.md`). A CLI tool then *automatically composes* them into canonical `AGENTS.md` files in the correct directories. This elegantly solves your "Ctrl+P" navigation problem: all your source fragments have unique, meaningful names, and the final `AGENTS.md` files are just builds, not files you edit directly.

*   **`agentic-docs-templates` (Structured Framework)**: This provides a holistic blueprint that goes beyond a single command file. It suggests a structured documentation system including `AGENTS.md`, `ARCHITECTURE.md`, and folders for `design-docs/`, `exec-plans/`, etc.. It uses a `bootstrap.md` prompt to have an AI agent analyze your existing project and auto-generate the entire documentation structure for you. This is ideal for bringing a large, messy repo up to speed.

*   **`AgentsMesh` (Universal Canonical Source)**: This is for a different problem: maintaining consistency across *different AI tools*. It lets you define all your rules, commands, and agents once in a canonical `.agentsmesh/` directory, and then it *generates* the native config files (`AGENTS.md`, `.cursorrules`, etc.) for every tool you use. You don't need this unless you're actively switching between many different coding agents.

These modern systems are designed to solve the exact challenges you're facing: scaling documentation, making it navigable for humans, and keeping it from becoming a burden.

### 💡 A Practical Path for FireCore & SimpleSimulationEngine

Given your desire for a robust yet practical system that works across agents, here's a path you can take.

**1. Foundation: Adopt the Standard**
*   **Start Simple**: Add a concise `AGENTS.md` to the root of your `FireCore` and `SimpleSimulationEngine` repositories.
*   **Seed the Agent**: Use an AI agent with a prompt like this to auto-generate your initial files:
    > "Analyze this project's codebase, tech stack, and development scripts. Using the `agentic-docs-templates/bootstrap.md` framework, generate a `AGENTS.md` file for the root, an `ARCHITECTURE.md`, and a structured documentation tree in a `docs/` folder based on the code's modules."

**2. Scale: Implement a Multi-File Pattern**
Choose between the **DOX** or **agents-md** approach for managing your `prokop` branch and larger repositories.

*   **For your `prokop` branch**, adopt the **DOX pattern**:
    *   Place a DOX-like `AGENTS.md` in the `prokop/` directory itself. This local file can state: "This is an experimental scratchpad. Code here may be unstable, but should adhere to core project patterns when possible." This gives the agent crystal-clear, context-specific instructions.

*   **For larger refactoring**, start with the **agents-md** pattern:
    *   Create source fragments for different concerns (e.g., `docs/build.agents.md`, `docs/api.agents.md`).
    *   Run the `agents-md compose` command to automatically build your `AGENTS.md` files.
    *   Consider using the `pre-commit` hook to keep them automatically generated and up-to-date.

**3. Optimize for Your Brain: Human Workflow Integration**
Your "Ctrl+P" concern is valid. Here's how to solve it:

*   **Use the `agents-md` Pattern**: As mentioned, you edit source fragments with meaningful names like `database.agents.md` and `api.agents.md`, which are easily searchable. The final `AGENTS.md` is a build artifact, not a source file.
*   **Create an Index File**: Keep a root `PROJECT_MAP.md` (a human-only file) as a simple table of contents for your documentation structure.
*   **Project-Specific Aliases**: Use VS Code's workspace settings (`keybindings.json`) to create keyboard shortcuts that jump to your key agent files.
    ```json
    {
      "key": "ctrl+shift+a",
      "command": "workbench.action.quickOpen",
      "args": "AGENTS.md"
    }
    ```

### 🧠 Mastering the Context Window: The Discipline of Engineering

This is the most crucial part of "vibe coding." You are right to be concerned about overwhelming agents with too much context. The agent's context window is like **RAM**: it's fast and powerful but finite and volatile. It is not meant for long-term storage. Using it as a database leads to failures that look like model problems but are actually memory architecture problems.

Here are the key strategies, which align perfectly with the multi-file systems discussed:

*   **"RAM vs. Disk" Model**: Put only the *current, active task's state* in the context window (e.g., the specific `AGENTS.md` and key files you're working on). Everything else (the bulk of your repo docs, historical decisions) should be "on disk," meaning it's *available to be retrieved* via tool calls, not loaded by default.
*   **Write to be Retrieved**: `AGENTS.md` files should be concise and highly scannable. Use clear headings and bullet points to help agents find what they need quickly. A 200-400 line file is a good target.
*   **Use Tool Calls to Fetch**: Instead of loading `ARCHITECTURE.md` into the main context, the agent should be instructed to *use a tool to read the file* only when it becomes relevant to the current task. This keeps the working memory clean.
*   **Compress and Summarize**: For long conversation histories or large tool outputs, have the agent generate a concise summary before starting a new task, and load only that summary into the context for the next step.

### 🧑‍💻 Conclusion: A Summary of Key Principles & Actions

Here are the core principles and actionable steps to take immediately.

**Core Principles**

*   **Perspective**: View the context window as a **constrained resource (RAM)**, not a storage disk.
*   **Strategy**: Keep documentation **layered** (global vs. local) and **lightweight** (concise, scannable).
*   **Action**: Give the agent the ability to **fetch details on-demand** via tool calls, rather than pre-loading everything.
*   **Control**: Structure your documentation to be **human-and-AI friendly**, using an approach like DOX or agents-md.

**Actionable Next Steps**
1.  **Audit & Initialize**: Use an AI agent with the `agentic-docs-templates/bootstrap.md` prompt to perform an initial audit and generate a structured documentation scaffold for both `FireCore` and `SimpleSimulationEngine`.
2.  **Seed the `prokop` Branch**: Place a concise DOX-like `AGENTS.md` file inside the `prokop/` directory with specific, experimental guidelines.
3.  **Update Your Workflow**: Run `agents-md compose` and set up the `pre-commit` hook to keep your `AGENTS.md` files automatically generated from source fragments.
4.  **Establish a Human Workflow**: Create an index file or VS Code aliases to navigate your documentation without the "Ctrl+P" friction.

By adopting these standards and practices, you'll transform your chaotic "vibe coding" into a structured, reliable system where AI agents become powerful, context-aware collaborators rather than unpredictable guessers. Good luck with FireCore, and feel free to reach out if you have further questions as you implement these systems!

---

# Grok

**AGENTS.md is the de facto standardized, transferable system for agentic coding context.** It works across most modern coding agents (Cursor, Claude Code, GitHub Copilot agent mode, Aider, OpenAI Codex, and others). It's intentionally simple, version-controlled, and predictable—exactly what you want for a large/experimental repo like FireCore (prokop branch) or SimpleSimulationEngine.

### Why AGENTS.md (and not just one big file)?
- **It's agent-first**: README.md is for humans (overviews, installation). AGENTS.md is the "README for agents" with build/test commands, architecture notes, coding conventions, gotchas, workflows, etc. Agents load it automatically at the start of sessions.
- **Hierarchical & scalable for large/experimental repos**: Place one at the root for global rules. Add nested `AGENTS.md` files in subdirectories (e.g., `src/physics/AGENTS.md`, `experiments/prokop/AGENTS.md`). Agents read the nearest one (plus parents) for context-specific guidance. This is perfect for monorepos or scratchpad-style branches—keeps context focused and avoids bloating a single file. OpenAI's own repos use dozens (e.g., 88 in one case).
- **Transferable**: Works with most tools without proprietary files. For Claude Code specifically, symlink or reference it (e.g., `ln -s AGENTS.md CLAUDE.md` or `@AGENTS.md` in prompts). Cursor and others support it alongside their own rules.

**Do not dump everything into one AGENTS.md**—that's a common pitfall that overwhelms context windows, increases token costs, and can even reduce agent performance (some studies show overly verbose context files hurt success rates). Keep it concise (ideally <150-500 lines total per file; focus on what the agent *can't infer* from code).

### Recommended Structure for Your Setup
Use a lightweight hierarchy + modular references to avoid navigation issues and context overload:

1. **Root AGENTS.md** (project-wide):
   - High-level purpose/architecture overview.
   - Build/test/run commands.
   - Global conventions (style, git workflow, experimentation rules for prokop branch).
   - Navigation/index: "For physics simulation → read @src/physics/AGENTS.md", "For prokop experiments → read @branches/prokop/AGENTS.md or specific subdirs".
   - Boundaries (e.g., "Don't touch production code in experiments without review").

2. **Nested AGENTS.md files** (domain/experiment-specific):
   - In key folders like `src/`, `experiments/`, `prokop/`, etc.
   - Local contracts, purpose, workflows, verification steps.
   - The dox system you linked (from agent0ai) is a solid example of this hierarchical approach—agents must walk the chain and update docs after changes.

3. **Supplementary files** (don't overload AGENTS.md):
   - `architecture.md` or `docs/OVERVIEW.md` for deeper human+agent docs.
   - `.cursor/rules/` (Cursor-specific) or similar tool folders for granular rules.
   - Use `@filename.md` references in AGENTS.md so agents pull in only what's needed.
   - For very experimental scratchpads: A dedicated `prokop/AGENTS.md` with loose rules ("fast iteration, document breaking changes").

This keeps files short, discoverable, and non-repetitive. Agents handle the hierarchy automatically.

### Addressing Your Concerns
- **Same filename navigation (Ctrl+P in VSCode)**: Inevitable with the standard, but hierarchy helps—deeper files are in subdirs so paths differ (e.g., `prokop/AGENTS.md`). Use a `.agents/` or `docs/` folder for non-hierarchical extras and reference them. Many devs live with it for the standardization benefit.
- **Overwhelming context**: Prioritize minimal, operational instructions. Update iteratively (let agents help refine after tasks). Test what works—some research suggests slim, focused files outperform exhaustive ones. Re-read relevant chain before edits.
- **Prokop branch as scratchpad**: Treat it with a strong local AGENTS.md emphasizing experimentation rules, quick validation, and "update docs on durable changes."

### Alternatives & Complements
- **Tool-specific**: CLAUDE.md (Claude), .cursorrules / .cursor/rules (Cursor), copilot-instructions.md, etc. Use AGENTS.md as the core and symlink/reference for compatibility.
- **Dox-like systems**: The one you found is great for strict hierarchy/contracts. Similar ideas exist in skills/best-practices repos.
- **Other patterns**: `.ai/docs/`, architecture.md + rules folders, or even semantic search tools/MCP for dynamic context (more advanced). But stick to AGENTS.md for portability.
- **No single alternative dominates**—AGENTS.md won because it's simple and cross-tool.

**Start simple**: Create a root AGENTS.md based on the dox example or agents.md spec. Let an agent (e.g., in Cursor) generate/refine an initial version for FireCore/prokop. Iterate based on real sessions. This will transfer well to any agent.

If you share more details about your repo structure or specific pain points (e.g., physics sims, experiments), I can help draft a sample AGENTS.md.

---

# Mistral

Your goal—**building a transferable, agent-friendly documentation/navigation system** for large, experimental repos like [FireCore/prokop](https://github.com/ProkopHapala/FireCore/tree/prokop) and [SimpleSimulationEngine](https://github.com/ProkopHapala/SimpleSimulationEngine)—is both timely and practical. Let’s break this down into actionable insights, alternatives, and best practices, tailored to your workflow and the constraints you mentioned.

## **1. Standardized Systems for Agent Context**
### **AGENTS.md: The De Facto Standard**
- **What it is**: A community-driven convention (popularized by [agent0ai/dox](https://github.com/agent0ai/dox/blob/main/AGENTS.md)) for providing agent-specific instructions, repo structure, and workflows in a single file.
- **Pros**:
  - Simple, widely recognized, and easy to adopt.
  - Works with most agents (Cursor, Devin, Claude Code, etc.) as they often look for `AGENTS.md` by default.
  - Centralized: All context in one place.
- **Cons**:
  - Can become unwieldy for large repos (e.g., FireCore/prokop with its experimental branches).
  - Hard to navigate in VSCode (`Ctrl+P`) if all files are named the same.

### **Alternatives/Complements to AGENTS.md**
| System/File          | Purpose                                                                 | Pros                                                                 | Cons                                                                 |
|----------------------|-------------------------------------------------------------------------|----------------------------------------------------------------------|----------------------------------------------------------------------|
| **AGENTS.md**        | General agent instructions, repo overview, workflows.                 | Universal, simple.                                                  | Can get too long; hard to modularize.                                |
| **CONTEXT.md**       | Focused on *contextual* info (e.g., key files, dependencies, quirks). | Separates "how to use" (AGENTS.md) from "what to know" (CONTEXT.md). | Less standardized; agents may not look for it.                      |
| **ARCHITECTURE.md**  | High-level design, module interactions, diagrams.                    | Great for complex repos like FireCore.                              | Overkill for small projects.                                         |
| **.agentignore**     | Specifies files/dirs to *exclude* from agent context.                  | Prevents context overload (e.g., ignore `node_modules/`).          | Not all agents support it yet.                                       |
| **TASKS.md**         | Common tasks/workflows (e.g., "How to run tests in prokop branch").   | Modular; agents can reference specific tasks.                      | Requires agents to parse multiple files.                             |
| **EXPERIMENTS.md**   | For branches like `prokop`: tracks experiments, hypotheses, results. | Keeps scratchpad work organized.                                    | Niche; may need custom agent prompts to use it.                       |
| **`/.vibe/` or `/ai/`** | Dedicated directory for agent files (e.g., `ai/AGENTS.md`).          | Avoids `Ctrl+P` clutter; allows multiple files (e.g., `ai/TASKS.md`). | Less discoverable; agents need to be told to look here.              |

---
## **2. Addressing Your Concerns**
### **A. Navigation in VSCode (`Ctrl+P`)**
- **Problem**: All `AGENTS.md` files blend together.
- **Solutions**:
  1. **Prefix files with repo/branch name**:
     - `firecore-prokop-AGENTS.md`
     - `simplesimulation-AGENTS.md`
     - *Downside*: Less portable (agents expect `AGENTS.md`).
  2. **Use a dedicated directory**:
     - `/ai/AGENTS.md` (for FireCore)
     - `/ai/EXPERIMENTS.md` (for prokop branch)
     - *Upside*: Cleaner `Ctrl+P`; agents can be prompted to "check the `/ai/` directory."
  3. **Leverage VSCode workspaces**:
     - Group related repos in a workspace; use `Ctrl+P` with workspace filters.

### **B. Avoiding Context Overload**
- **Problem**: Agents get overwhelmed with too much info.
- **Solutions**:
  1. **Modularize documentation**:
     - Split into `AGENTS.md` (core instructions) + `CONTEXT.md` (repo-specific details) + `TASKS.md` (workflows).
     - Example for FireCore/prokop:
       ```
       /ai/
       ├── AGENTS.md          # General agent instructions
       ├── CONTEXT.md         # Key files, dependencies, quirks
       ├── EXPERIMENTS.md     # Prokop branch-specific experiments
       └── TASKS.md           # Common tasks (e.g., "How to run tests")
       ```
  2. **Use `.agentignore`**:
     - Exclude noisy files/dirs (e.g., `build/`, `*.log`).
     - Example:
       ```
       # .agentignore
       node_modules/
       *.pyc
       /build/
       ```
  3. **Prioritize with headers**:
     - In `AGENTS.md`, use clear sections like:
       ```markdown
       ## 🚀 Quick Start (For Agents)
       ## 📁 Repo Structure
       ## ⚙️ Workflows
       ## 🔬 Prokop Branch Experiments (See EXPERIMENTS.md)
       ```

## **3. Transferable System Proposal**
Here’s a **scalable, agent-friendly** structure you can reuse across repos:

### **Core Files**
| File/Directory       | Purpose                                                                 | Example Content                                                                 |
|----------------------|-------------------------------------------------------------------------|---------------------------------------------------------------------------------|
| `/ai/AGENTS.md`      | **Primary entry point** for agents.                                     | Repo overview, how to navigate, where to find more info.                       |
| `/ai/CONTEXT.md`     | **Static context**: dependencies, key files, quirks.                  | "FireCore uses Eigen3 for linear algebra. Main entry: `src/main.cpp`."         |
| `/ai/TASKS.md`       | **Dynamic workflows**: how to build, test, debug.                       | "To run tests in prokop branch: `./run_multi.sh 2>1 | tee OUT-run-multi`." |
| `/ai/EXPERIMENTS.md` | **Branch-specific notes** (e.g., prokop).                               | "Experiment: Testing memory sanitizers. See `libasan.so` in `run_multi.sh`."   |
| `/ai/.agentignore`   | Files/dirs to exclude from agent context.                               | `build/`, `*.log`, `node_modules/`                                             |

### **Example `/ai/AGENTS.md` for FireCore/prokop**
```markdown
# FireCore Agent Guide
**Branch**: prokop (experimental)

## 🚀 Quick Start
- **Primary entry**: `src/main.cpp`
- **Build**: `mkdir build && cd build && cmake .. && make`
- **Run tests**: `./run_multi.sh 2>1 | tee OUT-run-multi`

## 📁 Repo Structure
- `/src/`: Core C++ libraries (e.g., `libMMFFmulti_lib.so`)
- `/scripts/`: Python helpers (e.g., `test_UFF_multi.py`)
- `/common_resources/`: Shared symlinks

## ⚙️ Workflows
- **Multiprocessing**: Set `OMP_NUM_THREADS` to control CPU usage.
- **Debugging**: Use `LD_PRELOAD=libasan.so` for memory sanitization.

## 🔬 Prokop Branch
- **Purpose**: Scratchpad for experiments.
- **Details**: See [EXPERIMENTS.md](./EXPERIMENTS.md)
- **Note**: Always leave 1 CPU core free (`ncpu=$(($ncpu - 1))`).

## 📌 Rules for Agents
1. **Do not modify**: `common_resources/` (symlinks only).
2. **Check first**: [CONTEXT.md](./CONTEXT.md) for dependencies.
3. **Log errors**: Use `2>1 | tee OUT-*` to capture stderr.
```

## **4. Other Systems to Explore**
| System               | Description                                                                 | Link                                                                 |
|----------------------|-----------------------------------------------------------------------------|----------------------------------------------------------------------|
| **Cursor Rules**     | Cursor-specific instructions (similar to AGENTS.md).                     | [Cursor Docs](https://docs.cursor.com/)                              |
| **Devin Config**     | Devin uses `devin.yaml` for project-specific settings.                     | [Devin GitHub](https://github.com/devin-ai/devin)                   |
| **OpenCode Rules**   | OpenCode looks for `.opencode/` directory.                                | [OpenCode Docs](https://github.com/open-code-ai/open-code)          |
| **Claude Code Rules**| Claude Code supports `CLAUDE.md` or `.claude/` directory.                  | [Claude Code Docs](https://github.com/anthropics/claude-code)        |
| **Vibe Rules**       | Vibe (Mistral) supports `VIBE.md` or `/vibe/` directory.                   | (Internal)                                                           |
| **AI Metadata**      | JSON/YAML files for structured agent context (e.g., `ai-metadata.json`). | [Example](https://github.com/abahmed/ai-metadata)                   |

## **5. Recommendations for Your Repos**
### **For FireCore/prokop**
1. **Start with `/ai/AGENTS.md`** as the primary entry point.
2. **Add `/ai/EXPERIMENTS.md`** to track your scratchpad work in the prokop branch.
3. **Use `/ai/.agentignore`** to exclude build artifacts and logs.
4. **Link to external docs**: If you have diagrams (e.g., architecture), reference them in `AGENTS.md`.

### **For SimpleSimulationEngine**
1. **Reuse the same `/ai/` structure** for consistency.
2. **Keep `AGENTS.md` lightweight** and offload details to `CONTEXT.md`/`TASKS.md`.

### **Agent-Specific Tweaks**
- **Cursor**: Add a `cursor.rules.md` file (Cursor will auto-detect it).
- **Devin**: Add a `devin.yaml` for project-specific settings.
- **Claude Code**: Add a `CLAUDE.md` or use `/ai/AGENTS.md` (Claude checks both).

## **6. Testing Your System**
1. **Prompt an agent**:
   - *"Read `/ai/AGENTS.md` and summarize the FireCore/prokop branch."*
   - *"How do I run tests in SimpleSimulationEngine? Check `/ai/TASKS.md`."*
2. **Iterate**:
   - If agents miss info, add it to `AGENTS.md` or link to the relevant file.
   - If agents get overwhelmed, split content further or use `.agentignore`.

## **7. Example Workflow for You**
1. **Set up `/ai/` in FireCore/prokop**:
   ```bash
   mkdir -p FireCore/ai
   touch FireCore/ai/{AGENTS.md,CONTEXT.md,TASKS.md,EXPERIMENTS.md,.agentignore}
   ```
2. **Populate `AGENTS.md`** with the template above.
3. **Add prokop-specific notes** to `EXPERIMENTS.md`:
   ```markdown
   # Prokop Branch Experiments
   - **Memory Sanitization**: Tested with `libasan.so` (see `run_multi.sh`).
   - **Multiprocessing**: Always leave 1 core free (`ncpu=$(($ncpu - 1))`).
   - **Logs**: Use `2>1 | tee OUT-*` to capture errors.
   ```
4. **Test with an agent**:
   - Ask it to *"Explain the prokop branch experiments in FireCore."*

## **8. Open Questions for You**
- Do you want to **automate** any part of this (e.g., a script to generate `CONTEXT.md` from repo structure)?
- Should we **standardize a template** for your `/ai/` directories across all repos?
- Are there **specific agents** you want to optimize for (e.g., Devin’s `devin.yaml`)?

---

# Kimi 2.6

I'll help you research context engineering systems for agentic coding. Let me search for the latest approaches and alternatives to AGENTS.md.


Great question — you're touching on a genuinely open problem in context engineering. The short answer: **AGENTS.md is not the only standardized system, and for a large experimental repo like yours, a multi-file approach is likely better than a single monolithic file.**

Here's what the landscape looks like as of mid-2026, and how to think about structuring context for your `prokop` branch specifically.

## The Current File Format Landscape

There are **multiple competing conventions**, not one standard:

| File | Tool(s) | Purpose |
|------|---------|---------|
| `AGENTS.md` | Codex, generic standard | Tool-agnostic project context |
| `CLAUDE.md` | Claude Code | Anthropic-specific instructions |
| `.cursorrules` | Cursor | IDE-embedded rules |
| `.github/copilot-instructions.md` | GitHub Copilot | Copilot agent instructions |
| `.codex/instructions.md` | Codex (legacy) | OpenAI-specific instructions |

The MSR 2026 study found that only **5% of repositories** have adopted any of these formats, and the most common pairing was actually `AGENTS.md` + `CLAUDE.md` coexisting in the same repo . So the field is very much in flux.

## The Multi-File Approach (Recommended for Your Case)

For a large repo with an experimental "scratchpad" branch, **don't put everything in one file**. The 2026 consensus has moved away from kitchen-sink documentation. Instead, use a **hierarchical context system**:

### 1. Root-level: `AGENTS.md` (the "contract")
Keep this under **250 lines** . It should contain:
- **Project identity** (1 paragraph: what FireCore is, stack, audience)
- **Non-negotiables** (5–10 hard, falsifiable rules with enforcement method)
- **Communication style** (language, tone, verbosity)
- **Principles** (3–7 rules that apply when no other rule applies)

**Example for your case:**
```markdown
## Project
FireCore is a physics simulation engine (C++/Python) focused on 
molecular dynamics and force-field calculations. The `prokop` 
branch is an experimental playground — code there may be 
prototype-quality and is not subject to the same stability 
requirements as `main`.

## Non-negotiables
- All new physics kernels must have a corresponding Python test 
  in `tests/prokop/`
- Experimental code in `prokop/` must be marked with 
  `# EXPERIMENTAL: [description]` header
- Never modify `main/` from the `prokop` branch without explicit 
  human approval
```

### 2. Directory-level: `.agents.md` or `AGENTS.md` in subfolders
This is where you solve the **Ctrl+P navigation problem**. Many tools (Claude Code, Cursor) will **walk up the directory tree** and aggregate context files. So you can have:

```
FireCore/
├── AGENTS.md                 # Root: global rules
├── prokop/
│   ├── AGENTS.md             # Branch-specific: "this is a playground"
│   ├── skymol/
│   │   └── AGENTS.md         # Subproject: specific force-field conventions
│   └── new_grid/
│       └── AGENTS.md         # Subproject: grid algorithm patterns
```

**The key insight:** Agents read context files recursively. A file in `prokop/skymol/AGENTS.md` will be loaded *in addition to* the root `AGENTS.md`, not instead of it. This lets you keep the root file stable while iterating on branch-specific context.

### 3. Tool-specific overlays (optional)
If you use multiple agents, maintain **one canonical source** and generate tool-specific files:

```
scripts/
├── generate-context.sh       # Syncs rules.base.md → AGENTS.md, CLAUDE.md, .cursorrules
context/
├── rules.base.md             # Your single source of truth
├── AGENTS.md                 # Generated (tool-agnostic)
├── CLAUDE.md                 # Generated + Claude-specific terminal behavior
└── .cursorrules              # Generated + Cursor-specific diff formatting
```

This avoids the "contradiction death spiral" that kills multi-tool projects .

## Alternative Systems & Emerging Approaches

Beyond the markdown-file approach, there are several alternatives worth knowing:

### 1. **Skills / Lazy-Loaded Context** (Claude Code, Cursor)
Instead of loading everything at session start, you define **skills** that the agent loads on-demand. A skill is a folder with a `SKILL.md` and supporting files.

```
.skills/
├── physics-kernels/
│   ├── SKILL.md              # "How to write a new physics kernel"
│   ├── template.cpp          # Boilerplate
│   └── example_test.py       # Reference test
├── build-system/
│   └── SKILL.md              # "How to modify CMake for new modules"
└── git-workflow/
    └── SKILL.md              # "How we commit in prokop branch"
```

The agent reads the skill description and decides whether to load the full content. This is **the best solution for not overloading context** .

### 2. **Ralph Loop / Git-as-Memory**
For your experimental branch, consider Geoffrey Huntley's approach: spawn a **fresh agent session per task**, let it read the PRD/code, commit, and terminate. Git becomes the persistent memory layer, not the context file . This is ideal for "scratchpad" work where you don't want accumulated context to bias the next experiment.

### 3. **MCP Servers (Model Context Protocol)**
Instead of static files, expose dynamic context through MCP servers. For example, an MCP server could:
- Query your `prokop` branch's recent commits to summarize current experiments
- Read your `SimpleSimulationEngine` repo to cross-reference patterns
- Generate on-the-fly context based on the specific task

This is more infrastructure but avoids stale context entirely.

### 4. **Spec-Driven Development (spec-kitty, etc.)**
For experimental features, write a **spec file** (acceptance criteria, testable requirements) and let the agent implement against it. The spec becomes the context, not a general rules file .

## Specific Recommendations for Your Setup

Given your two repos (`FireCore/prokop` as playground, `SimpleSimulationEngine` as reference), here's a practical structure:

### FireCore repo:
```
AGENTS.md                      # Global: "this is a physics engine"
prokop/
├── AGENTS.md                  # "This branch is experimental. Rules are relaxed."
├── EXPERIMENTS.md             # Human-readable log of what's being tried
├── .skills/                   # Lazy-loaded context
│   ├── force-field-dev/
│   ├── grid-algorithms/
│   └── cpp-python-binding/
└── [subproject dirs with their own .agents.md]
```

### SimpleSimulationEngine repo:
```
AGENTS.md                      # "Reference engine. Stable patterns live here."
.skills/
├── architecture-patterns/     # "How we structure simulation loops"
├── math-utilities/            # "Vector/matrix conventions"
└── rendering-integration/     # "OpenGL/Vulkan patterns"
```

### Cross-repo reference:
In `FireCore/prokop/AGENTS.md`, add:
```markdown
## Reference Codebase
For stable patterns and architectural decisions, consult 
`SimpleSimulationEngine` (sibling repo). Prefer copying patterns 
from there over inventing new ones, unless the experiment 
explicitly requires divergence.
```

## The "Ctrl+P Problem" — Practical Workarounds

You mentioned all files having the same name makes navigation hard. A few tricks:

1. **Use `.agents.md` (hidden) in subdirectories** — VS Code's Ctrl+P filters hidden files by default, but you can type `.agents` to find them.
2. **Symlink with descriptive names**: `prokop-experiments.agents.md → prokop/AGENTS.md` (for human navigation only; agents still read the canonical file).
3. **Use VS Code workspace settings** to exclude `AGENTS.md` from quick-open if you have too many, and navigate via the Explorer instead.

## Bottom Line

For your use case, **don't use a single `AGENTS.md`**. Use a **hierarchical, lazy-loaded system**:

- **Root `AGENTS.md`**: 250-line "constitution"
- **Branch `prokop/AGENTS.md`**: Experimental ground rules
- **Subproject `.agents.md`**: Specific technical conventions
- **`.skills/` folders**: Deep context loaded on-demand
- **Tool-specific files**: Generated from a single `rules.base.md` if you use multiple agents

This gives you the transferability you want (any agent can read markdown), avoids context overload, and lets your `prokop` branch breathe without polluting the global context.

---

# Perplexity

For your use case, **AGENTS.md is the best default backbone, but it should not hold everything**. The strongest pattern is a short root `AGENTS.md` plus more specific docs that the agent only reads when needed, because keeping the main file concise is explicitly recommended and nested `AGENTS.md` files can override by directory scope. [agents](https://agents.md)

## What to put where

Use the root `AGENTS.md` for stable, high-value guidance: repo layout, build/test commands, coding conventions, do-not rules, and what “done” means. [developers.openai](https://developers.openai.com/codex/learn/best-practices)
Put deeper or more specialized knowledge in separate markdown files, like `docs/architecture.md`, `docs/debugging.md`, `docs/dev-notes.md`, or task-specific playbooks, and link to them from `AGENTS.md` only when relevant. [philschmid](https://www.philschmid.de/writing-good-agents)
This avoids flooding every agent session with your entire scratchpad while still preserving useful context for harder tasks. [developers.openai](https://developers.openai.com/codex/learn/best-practices)

## Best structure for your repo

For a large repo like FireCore and especially a noisy experimentation branch, I’d use a **progressive disclosure** setup:  
- `AGENTS.md` at repo root for universal rules.  
- Additional `AGENTS.md` files inside major subprojects or packages when their workflows differ, because the closest file wins. [agents](https://agents.md)
- Separate support docs for architecture, experiments, and “how to think about this subsystem” notes, so the agent can fetch them only when needed. [philschmid](https://www.philschmid.de/writing-good-agents)

A practical example: root `AGENTS.md` says where the core docs live, what commands validate changes, and which directories are experimental; then `scratchpad/README.md` or `experiments.md` can explain the branch-specific playground conventions without forcing that context into every task. [developers.openai](https://developers.openai.com/codex/learn/best-practices)

## Other file systems

Yes, there are alternatives, but they are mostly tool-specific rather than universal. Common ones include `.cursorrules` or `.cursor/rules/*.mdc` for Cursor, `CLAUDE.md` for Claude Code, `.github/copilot-instructions.md` plus path-specific instruction files for GitHub Copilot, and `llms.txt` for documentation discovery rather than agent behavior. [deployhq](https://www.deployhq.com/blog/ai-coding-config-files-guide)
If you want something transferable across many agents, `AGENTS.md` is the closest thing to a shared standard, while the others are useful only when you specifically want to tune one tool. [agents](https://agents.md)
For repo navigation by agents, `llms.txt` is a good complement, because it maps important docs and can point agents toward canonical pages while excluding low-value or experimental material. [mintlify](https://www.mintlify.com/docs/ai/llmstxt)

## Keeping context small

Your concern about overfeeding agents is exactly right, and the best mitigation is to keep the main instructions short and push detail into linked docs that are read only when relevant. [philschmid](https://www.philschmid.de/writing-good-agents)
A good rule is: if a rule applies to almost every task, keep it in `AGENTS.md`; if it applies only to debugging, architecture work, tests, or a specific subsystem, move it out. [developers.openai](https://developers.openai.com/codex/learn/best-practices)
That gives you a stable entry point without turning the agent context into a junk drawer. [developers.openai](https://developers.openai.com/codex/concepts/customization)

## Recommended setup

For your two repos, I’d start with this layout:  
- `AGENTS.md`: short root-level agent handbook.  
- `docs/llms.txt`: doc map for humans and agents.  
- `docs/architecture.md`: high-level system map.  
- `docs/experiments/`: scratchpad and branch-specific notes.  
- Optional nested `AGENTS.md` files inside major subtrees if workflows diverge. [mintlify](https://www.mintlify.com/docs/ai/llmstxt)

That gives you one transferable pattern that works across Cursor, Devin, Codex, Claude Code, opencode, and similar tools, while still keeping the active context tight. [agents](https://agents.md)
Would you like a concrete `AGENTS.md` template tailored to FireCore and a separate doc tree for the prokop scratchpad branch?