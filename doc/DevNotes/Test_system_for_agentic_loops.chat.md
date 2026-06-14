# USER

Can you find some examples of writing agentic loops? ( Like in openClaw, Hermes agent or other AI harness) I'm thinking what I can use it for. I'm computational physicist and chemist. I have massive codebases and projects on GitHub but my work is more exploratory and my ropos disorganized. I have also lot of disorganized data, notes, projects exploratory unfinished stuff on my computer. I would like if some AI help me organize it finish or refactor. I often reimplement one code to different platform (eg fortran to c++, c++ to opencl or even to JavaScript webgpu with compute shaders. What I really need in my repos is define clearly tests so I can identify when I break something but it will also tremendous help with reimplementing (just check if the rests give same results like parity cheks. Maybe sometime we do not want exactly same results as new method can give potentially better results then the old, but must be qualitatively similar and qualitatively meaningful (like conservation law, symmetry, common sense). For example if I have many forcefields for water NaCl interaction, oxygen must sit on Na+ ion, and hydrogen toward Cl- ion in distance close to sum of atomic radius but the exact distance and energy is unclear unless we have accurate DFT reference. We must be able to define good systematic multilevel tests, but these tests must be concise and fast we have huge repos on GitHub

---

# Gemini

Building custom agentic loops is an excellent way to handle this kind of "scientific code archaeology" and cross-platform porting. Off-the-shelf AI tools often fail here because they expect standard web-dev frameworks and 100% deterministic unit tests, whereas you need **multilevel validation** (from bit-wise parity to physical sanity checks like conservation laws and structural configurations).

To achieve this, you don't necessarily need massive frameworks like LangChain or AutoGPT. In fact, for scientific computing, a lightweight, highly deterministic loop—similar to the architecture of projects like **openClaw** or **Hermes**—works best.

Here is a breakdown of how you can structure these agentic loops, with concrete examples tailored to your workflow.

## 1. The Core Architecture of a Scientific Coding Agent

An agentic loop is essentially a `while` loop that allows an LLM to reason, take an action (run a script, compile code), observe the output, and correct its course until a specific condition is met.

For your needs, a robust pattern consists of four steps: **Plan $\rightarrow$ Act $\rightarrow$ Test/Observe $\rightarrow$ Reflect**.

```
[Plan] -> [Act (Write/Port Code)] -> [Test (Compile & Phys Check)] -> [Reflect/Fix]
  ^                                                                        |
  |_______________________ Fail (Loop Max N) ______________________________|

```

### Example 1: The "Cross-Platform Porting & Parity" Loop

This loop specializes in porting code (e.g., C++ to WebGPU) and running differential tests (parity checks) between the legacy and new implementations.

```python
import os
import subprocess
from openai import OpenAI

client = OpenAI()

def agentic_porting_loop(source_code, target_platform, physical_constraints, max_iterations=5):
    # Context given to the LLM to act as a rigorous scientific developer
    system_prompt = f"You are an expert in computational physics and {target_platform}. Port the code and satisfy these physical rules: {physical_constraints}"
    
    # 1. INITIAL PLAN & CODE GENERATION
    user_prompt = f"Port this source code to {target_platform}:\n\n{source_code}"
    messages = [
        {"role": "system", "content": system_prompt},
        {"role": "user", "content": user_prompt}
    ]
    
    for i in range(max_iterations):
        print(f"\n--- Iteration {i+1} ---")
        
        # ACT: Agent generates or fixes the code
        response = client.chat.completions.create(model="gpt-4o", messages=messages)
        agent_output = response.choices[0].message.content
        messages.append({"role": "assistant", "content": agent_output})
        
        # Extract code block (helper function assumed here)
        generated_code = extract_code(agent_output) 
        write_to_disk(generated_code, target_platform)
        
        # OBSERVE: Compile and Run physical validations
        compile_success, compile_log = compile_code(target_platform)
        if not compile_success:
            print("Status: Compilation Failed. Feeding logs back to agent.")
            messages.append({"role": "user", "content": f"Compilation failed with error:\n{compile_log}\nFix the syntax/types."})
            continue
            
        run_success, run_log, physical_metrics = run_validation_suite()
        if not run_success:
            print("Status: Execution Runtime Error.")
            messages.append({"role": "user", "content": f"Runtime error:\n{run_log}\nFix the bug."})
            continue
            
        # REFLECT: Evaluate physical sanity checks
        passed, physics_report = evaluate_physics(physical_metrics)
        if passed:
            print("Status: SUCCESS! Code compiles and physical constraints are met.")
            return generated_code
        else:
            print("Status: Physical constraints violated. Forcing reflection.")
            messages.append({"role": "user", "content": f"Physical validation failed:\n{physics_report}\nAdjust the algorithm."})
            
    print("Loop reached max iterations without perfect convergence.")
    return None

```

## 2. Implementing Your Multilevel "Physical" Testing Framework

To make the agentic loop above work, you need a testing harness that understands that **numerical drift is expected, but physical violations are failures**. Standard testing frameworks don't handle this well, so you can write a lightweight Python evaluator that calculates **Qualitative & Physical Sanity**.

Instead of checking if $E_{new} == E_{old}$, your automated test harness should evaluate a multilevel validation checklist:

### Step 1: Bitwise/Delta Check (For direct ports like Fortran $\rightarrow$ C++)

* Accepts a strict tolerance: $|A_{i} - B_{i}| < \epsilon$.

### Step 2: Conservation Laws & Symmetries (For algorithmic changes)

* Is total energy conserved ($\Delta E \approx 0$ in NVE)?
* Is the force matrix strictly symmetric ($F_{ij} = -F_{ji}$)?

### Step 3: Structural/Geometrical Common Sense (The "Water-NaCl" Test)

* Compute the Radial Distribution Function ($g(r)$) or minimum distances during a short test run.
* **Assertion:** $r(Na^+ - O) < r(Na^+ - H)$ and the first peak falls within a physical threshold (e.g., $2.2 - 2.6$ Å).

Here is an example of a test validator script the agentic loop can execute:

```python
import numpy as np

def run_validation_suite():
    # 1. Run the compiled artifact (e.g., WebGPU/OpenCL headless script)
    # 2. Parse the output trajectory/matrices
    coordinates = np.load("test_output_coords.npy") # Shape: (N_atoms, 3)
    forces = np.load("test_output_forces.npy")
    
    report = []
    passed = True
    
    # Test A: Symmetry Check
    max_symmetry_error = np.max(np.abs(forces + forces.T)) # If matrix is F_ij
    if max_symmetry_error > 1e-4:
        report.append(f"FAIL: Newton's Third Law violated. Max F_ij + F_ji = {max_symmetry_error}")
        passed = False
        
    # Test B: Water-NaCl Orientation Check (Geometric Common Sense)
    # Dummy logic looking at atom indices: Na=0, O=1, H1=2, H2=3
    r_Na_O = np.linalg.norm(coordinates[0] - coordinates[1])
    r_Na_H1 = np.linalg.norm(coordinates[0] - coordinates[2])
    
    if r_Na_O >= r_Na_H1:
        report.append(f"FAIL: Physical orientation inverted! Na+ is closer to H than O. r_Na_O={r_Na_O:.2f}, r_Na_H1={r_Na_H1:.2f}")
        passed = False
        
    return passed, "\n".join(report)

```

## 3. Organizing the "Disorganized Chaos" (The Repo-Sweeper Loop)

For your scattered notes, unfinished snippets, and messy GitHub repos, you can deploy a **"Map-Reduce" Agentic Loop** that doesn't edit code immediately, but indexes it using physical metadata.

### The Metadata Extractor Loop

You can point a local loop at your directory to generate a structured map of your chaos. Have the agent read a file, evaluate its purpose, and output a strict JSON schema.

```python
import glob

def chaos_sweeper(root_dir):
    all_files = glob.glob(f"{root_dir}/**/*.*", recursive=True)
    
    for file_path in all_files:
        if any(ext in file_path for ext in ['.cpp', '.f90', '.cl', '.js', '.cu']):
            with open(file_path, 'r') as f:
                content = f.read()[:2000] # Grab the header/first 2000 characters
                
            prompt = f"""
            Analyze this raw file path and snippet: '{file_path}'.
            Categorize it by:
            1. Physical System (e.g., Water, NaCl, Lennard-Jones, DFT, Tight-Binding)
            2. Language/Platform (e.g., Fortran, OpenCL, WebGPU)
            3. Completeness (Scale 1-5: 1 is a broken snippet, 5 is a fully working solver)
            4. Does it have tests? (Yes/No)
            
            Return strictly as valid JSON.
            """
            # Call LLM here, parse JSON, and append to a centralized database (or markdown index)

```

By sweeping your hard drive with this metadata extractor, you instantly get an inventory of your unfinished work. You can then sort by `Completeness == 4` and target them for automated testing.

## 4. Practical Roadmap for Your Setup

Instead of trying to configure a massive, fragile system, you can build your own minimal harness step by step:

1. **Define the Fast Test Targets:** Pick *one* specific calculation (e.g., a short 100-step MD run of 64 water molecules). Write a reference output file using your trusted, slowest, most accurate CPU code (even if it's in old Fortran).
2. **Write the Physical Oracle:** Create a standalone Python script (`validate.py`) that checks if a newly generated output directory matches the reference qualitatively (energy drift bounds, structural conservation).
3. **Deploy the File-to-File Loop:** Write a simple Python script using the native OpenAI, Anthropic, or local Ollama SDK that reads a file, sends it to the LLM with compilation error logs, and writes the output back to disk.

This gives you a lightweight, deterministic environment where the AI handles the repetitive syntax adjustments of compute shaders, while your physical assertions prevent it from hallucinating the science.

What specific code-porting project or disorganized repository on your disk would be the best candidate to start testing this loop structure on?

---

# DeepSeek-V4-flash

Agentic loops are a natural fit for the exploratory, multi-repo, and multi-language challenges you've described. They shift from manually prompting for code to building automated systems that plan, execute, and verify tasks like cross-language porting or refactoring.

I have found several frameworks with clear "agentic loop" patterns that align well with your computational science and GPU acceleration goals. You can think of these loops as automated research assistants that work through a structured cycle of tasks, constantly evaluating and adjusting their approach.

### 🔄 Foundational Agentic Loops: Core Frameworks

These frameworks provide the basic building blocks for an AI agent that can think, use tools, and act:

*   **OpenClaw**: Its loop is structured as `intake → context assembly → model inference → tool execution → streaming replies → persistence`, treating each session as a single, serialized run to keep state consistent.
*   **Attractor's Coding Agent Loop Spec**: This is a **language-agnostic specification** that defines a programmable loop where a host application can inspect, steer, and control the agent at every step.
*   **Hermes Agent**: Features a **"Closed Learning Loop"** that automatically logs tasks, generates reusable skills, and stores them in long-term memory to improve future performance.
*   **Looplet**: A framework-agnostic, iterator-first library that gives you granular control. Instead of `agent.run()`, you use a `for step in loop(...):` to own your control flow.

### 🏗️ Agentic Engineering Harnesses: Orchestrators and Systems

These tools build upon the core loops to create full-fledged systems for managing complex software tasks:

*   **Agent Loop**: Uses an **Orchestrator agent** to break down goals and delegate tasks to specialized Specialist Agents (like `Architect`, `Builder`, `Tester`).
*   **Loop-Harness**: An autonomous "loop engineering" harness. Each loop wakes on a schedule, works in isolated git worktrees, and uses a **second agent to verify the work** of the first before merging.
*   **Whilly Orchestrator**: Runs a loop of `pick a pending task → hand it to an LLM agent → verify result → commit → next`, similar to the "Ralph Wiggum loop" but with more structure.
*   **SwarmWeaver**: A Python harness for long-running autonomous coding sessions across unlimited sessions, with multi-agent swarms for parallel execution.
*   **Lumpcode CLI**: Automates batched coding work (e.g., "migrate every component") by running a coding agent over similar units in isolated git branches.

### 🧪 Agentic Testing & Science Agents: Ensuring Correctness

For verification, you can implement property-based testing or use dedicated "science agents":

*   **Property-Based Testing**: Define high-level properties (e.g., "oxygen must sit on Na+ ion") and have the AI generate random data to check if they hold.
*   **AstraAI**: A structure-aware AI coding agent designed for large scientific HPC codebases like AMReX, integrating compiler tooling directly into the generation loop for deterministic, structure-preserving program transformation.

### 🌐 Porting & Refactoring Agents: Cross-Language Translation

Automate the reimplementation of code across platforms with AI-assisted tools:

*   **AutoSoftware**: An autonomous engineering platform where AI agents scan codebases, plan features, and write production code.
*   **RepoAI**: A multi-agent system for repository-level refactoring that uses **cross-family validation** (independent LLMs) to judge code quality.

### 📊 Comparative Framework Table

Here is a quick comparison to help you choose the right tool:

| **Framework** | **Core Loop Pattern** | **Target Use Case** | **Key Differentiator** |
| :--- | :--- | :--- | :--- |
| **OpenClaw** / **Hermes** | Stateful, memory-augmented loop | General-purpose & personal AI agents | Has persistent memory and learns from past interactions |
| **Attractor Spec** | Programmable, provider-aligned loop | Building custom, low-level coding agents | Gives you full, granular control over the loop |
| **Agent Loop** | Orchestrator with specialist agents | AI-driven software development teams | Uses a team of specialized AIs to tackle complex tasks |
| **Loop-Harness** | 2-agent verification loop | Autonomous, scheduled code maintenance | Has a built-in second-agent verification gate for safety |
| **SwarmWeaver** | Multi-agent, parallel execution loops | Long-running, large-scale codebase work | Supports parallel agents for tackling big tasks |
| **Whilly** | Plan-Execute-Verify-Repair | End-to-end issue to pull-request automation | Structured pipeline is excellent for systematic checks |
| **Lumpcode** | Batch processing loop over units | Repetitive coding chores & large refactors | Perfect for sweeping changes across many files |

These agentic loops can automate your cross-platform coding, organize your messy repositories, and finally enforce the rigorous testing your work deserves. Given your background, which of these areas feels most urgent to tackle first? I can help you think through a practical first step.

---

# Mistral

Here are **practical agentic loop examples** and how you can apply them to organize, test, and reimplement your scientific codebases:

## 🔁 **Core Agentic Loop Patterns**

| Pattern | How It Works | Best For | Example Frameworks |
|---------|-------------|----------|-------------------|
| **ReAct** | **Reason → Act → Observe → Repeat** (thought-action-observation cycle) | Exploratory tasks, debugging, unfamiliar codebases | LangChain AgentExecutor, Claude Code |
| **Reflexion** | ReAct + **explicit self-reflection** after failures (stores "lessons learned" in memory) | Trial-and-error tasks, debugging, creative problem-solving | [Reflexion paper](https://arxiv.org/abs/2303.11366) |
| **Ralph Loop** | **Task-based iteration**: Pick task → build → test → commit → repeat | Code implementation, test coverage, code review | [ralph-loop](https://github.com/agenticloops-ai/ralph-loop) |
| **Plan-and-Execute** | Plan full sequence → Execute → Re-plan on errors | Known domains, structured workflows | AutoGen, CrewAI |
| **Agentic RAG** | Retrieval **inside** reasoning cycle (decides when/what to retrieve dynamically) | Knowledge-heavy tasks, large codebases | Custom implementations |

## 🎯 **For Your Scientific Computing Needs**

### **1. Repository Organization & Refactoring**
**Use Case**: Clean up disorganized repos, finish unfinished projects, migrate between languages (Fortran→C++→OpenCL→WebGPU)

**Agentic Loop Application**:
```
Loop:
  1. AGENT: Scan repo structure, identify patterns, dead code, unfinished features
  2. AGENT: Generate inventory (files, dependencies, TODO comments, branches)
  3. HUMAN: Review inventory, set priorities
  4. AGENT: For each priority:
     - Write migration plan (language A → B)
     - Implement in small chunks
     - Run parity tests (see below)
     - Commit with descriptive message
  5. Repeat until inventory empty
```

**Tools**:
- **Hermes Agent**: Persistent memory across sessions, learns your codebase conventions
- **ralph-loop**: Pre-built prompts for `implement`, `test`, `review` cycles
- **OpenClaw**: Can execute shell commands, edit files, run builds directly on your machine

**Example Workflow** (using Hermes):
```bash
# Set up persistent agent for your repo
hermes init my-science-repo
hermes skill add code-organizer
hermes task "Inventory all Fortran files, categorize by: production, experimental, deprecated"
hermes task "Create migration plan: Fortran → C++ for production modules"
```

### **2. Multi-Level Testing Framework**

For scientific code, you need **hierarchical validation**:

| Level | Test Type | Example (Water-NaCl Forcefield) | Speed | When to Run |
|-------|-----------|----------------------------------|-------|-------------|
| **L0** | Syntax/Compilation | Does code compile? | ms | Every commit |
| **L1** | Unit Tests | Energy calculation for single water molecule | ms | Every commit |
| **L2** | **Parity Tests** | Same input → same output (old vs new implementation) | ms-s | Every PR |
| **L3** | **Qualitative Tests** | Oxygen near Na⁺, hydrogen near Cl⁻ (conservation checks) | s | Nightly |
| **L4** | **Reference Tests** | Compare to DFT/ experimental data | min-hrs | Weekly |
| **L5** | Regression | Full simulation suite | hrs | Release |

**Agentic Testing Loop**:
```python
# Pseudocode for your test runner
def test_implementation(old_code, new_code, test_cases):
    for case in test_cases:
        # L2: Parity check
        old_result = old_code(case.input)
        new_result = new_code(case.input)
        
        # For numerical: allow tolerance
        if not np.allclose(old_result, new_result, rtol=1e-5):
            # L3: Qualitative check
            if not qualitative_check(new_result):
                return FAIL
            else:
                log_warning("Numerical difference but qualitatively valid")
        
        # L3: Conservation laws
        if not check_conservation_laws(new_result):
            return FAIL
    
    return PASS
```

**Agent Integration**:
- Agent **automatically generates** L2 parity tests when you start a reimplementation
- Agent **suggests** L3 qualitative checks based on domain knowledge
- Agent **flags** when L4 reference tests deviate beyond expected tolerance

### **3. Reimplementation Workflow**

**Problem**: Porting Fortran → WebGPU compute shaders

**Agentic Loop**:
```
1. AGENT: Extract core algorithm from Fortran (identify hot loops)
2. AGENT: Generate equivalent WebGPU kernel
3. AGENT: Create parity test harness:
   - Generate random inputs
   - Run both implementations
   - Compare outputs with tolerance
4. AGENT: For qualitative validation:
   - Check physical constraints (energy ≥ 0)
   - Check symmetry properties
   - Check conservation laws
5. HUMAN: Review edge cases, approve
6. AGENT: Optimize WebGPU version
7. Repeat for next module
```

**Example with OpenClaw**:
```bash
# OpenClaw can run this autonomously
openclaw task "Port force_calc.f90 to WebGPU:
1. Read force_calc.f90
2. Identify the main computation loop
3. Write equivalent WGSL compute shader
4. Create test in tests/parity/force_calc_test.py
5. Run test suite
6. If fails, debug and fix
7. Commit changes"
```

### **4. Concrete Frameworks to Try**

| Framework | Strengths | Your Use Case |
|-----------|-----------|---------------|
| **[OpenClaw](https://openclaw.ai/)** | Runs on your machine, direct file/process access | Full repo management, code execution |
| **[Hermes Agent](https://hermes-agent.nousresearch.com/)** | Learning loop, persistent memory | Long-term codebase understanding |
| **[ralph-loop](https://github.com/agenticloops-ai/ralph-loop)** | Simple, task-focused | Structured refactoring sessions |
| **[Agentic AI Engineering](https://github.com/agenticloops-ai/agentic-ai-engineering)** | Tutorials, first principles | Learn to build custom loops |

**Quick Start with ralph-loop**:
```bash
git clone https://github.com/agenticloops-ai/ralph-loop
cd ralph-loop
# Copy the structure to your project
cp -r .ralph ~/your-science-repo/
cd ~/your-science-repo
# Edit .ralph/prompts/implement.md to add:
# "For scientific code: always add parity tests when refactoring"
./ralph-loop.sh 10  # Run 10 autonomous iterations
```

### **5. Testing Patterns for Scientific Code**

**Parity Test Generator** (Agent writes these automatically):
```python
# Agent-generated test for your water-NaCl forcefield
def test_na_cl_interaction():
    # Generate test case
    na_pos = np.array([0, 0, 0])
    cl_pos = np.array([3, 0, 0])
    water = WaterMolecule([0, 2, 0], [1, 2.5, 0], [0, 2, 1])

    # Run both implementations
    old_energy = old_forcefield.energy(na_pos, cl_pos, water)
    new_energy = new_forcefield.energy(na_pos, cl_pos, water)

    # L2: Numerical parity
    assert np.isclose(old_energy, new_energy, rtol=1e-4)

    # L3: Qualitative checks
    forces = new_forcefield.forces(na_pos, cl_pos, water)
    assert forces.on_oxygen[0] > 0  # Oxygen attracted to Na+
    assert forces.on_hydrogen[0] < 0  # Hydrogen attracted to Cl-
    
    # L3: Distance check
    o_na_dist = np.linalg.norm(water.oxygen_pos - na_pos)
    assert 2.0 < o_na_dist < 3.0  # Reasonable distance
```

**Conservation Law Test**:
```python
def test_energy_conservation(simulation):
    # Run simulation for N steps
    energies = [simulation.total_energy() for _ in range(1000)]
    
    # Check conservation (allow small drift)
    max_drift = max(energies) - min(energies)
    assert max_drift / np.mean(energies) < 1e-6
```

### **6. What You Should Do Next**

1. **Start Small**: Pick **one disorganized repo** and run an inventory agent
2. **Set Up Testing**: Create a `tests/parity/` directory with agent-generated tests
3. **Choose a Framework**:
   - **OpenClaw** if you want direct machine access
   - **Hermes** if you want persistent learning
   - **ralph-loop** if you want simple, task-focused iterations
4. **Define Your Qualitative Rules**: Document domain-specific checks (conservation laws, symmetries, physical constraints)
5. **Automate the Loop**: Let the agent handle the repetitive parts (test generation, parity checking, commit messages)

**Example First Command**:
```bash
# With OpenClaw
openclaw task "Analyze ~/code/quantum-chemistry:
1. List all files by language
2. Identify all TODO/FIXME comments
3. Find all test files
4. Create a report in REPO_AUDIT.md
5. Suggest next steps for organization"
```

Would you like me to elaborate on any specific pattern or help you set up a pilot agentic loop for one of your repos?

---

# Kimi 2.6

I'll search for examples of agentic loops and how they can be applied to scientific computing workflows.
Let me search for more specific applications to scientific computing and cross-platform code migration
Based on my research, I'll give you a comprehensive overview of agentic loops and how they apply directly to your computational physics/chemistry workflow. The core pattern is simple but powerful, and there are several mature frameworks you can use or adapt.

## The Core Agentic Loop Pattern

All modern agentic systems (OpenClaw, Hermes, Claude Code, Aider) share the same fundamental **ReAct loop** (Reason + Act):

```python
while True:
    response = llm.call(context)

    if response.is_text():
        send_reply(response.text)
        break

    if response.is_tool_call():
        result = execute_tool(response.tool_name, response.tool_params)
        context.add_message("tool_result", result)
        # loop continues — model sees the result and decides next action
```

This is the entire engine. The model either replies (done) or requests a tool (file read, shell command, test run, etc.), sees the result, and loops. What separates the frameworks is what they wrap around this loop: budgets, guardrails, memory, skills, sub-agents, and tool registries. 

## Frameworks You Can Use Today

### 1. **OpenClaw** — Self-hosted, ReAct-based, heartbeat monitoring
- Uses the **Pi Agent framework** (TypeScript) with streaming
- **Skills**: Markdown files with YAML frontmatter that load on-demand
- **Heartbeat**: Periodic autonomous checks (e.g., every 30 min) via `HEARTBEAT.md` — perfect for long-running HPC jobs or test suites
- Memory: SQLite + Markdown files, embedding search via `sqlite-vec`
- Good for: Local control, privacy, custom tool integration 

### 2. **Hermes Agent** — Most sophisticated loop architecture
- **IterationBudget**: Thread-safe counter with mutex, separate budgets for parent/subagents (default 90 for parent, 50 per subagent)
- **Three-tier prompts**: Stable (cached), Context (project files), Volatile (memory, timestamps)
- **Self-improving loop**: After tasks, it extracts patterns, creates skills, and refines them
- **Tool registry + toolsets**: Global registry for all tools, but curatorial toolsets control blast radius
- **Subagents**: `delegate_task` spawns isolated children with separate budgets and tool grants
- Good for: Complex, multi-step scientific workflows that need to learn from repetition 

### 3. **Claude Code / Aider** — Immediate practical use
- **Claude Code**: Autonomous agent that searches codebase, plans multi-step tasks, runs tests, fixes failures in a loop. Sub-agents for parallel work. MCP protocol for external tools.
- **Aider**: Open-source, git-first, model-agnostic. Repo-map via tree-sitter. Every edit is a git commit. 

---

## How This Maps to Your Scientific Computing Needs

### **1. Codebase Organization & Refactoring**

For your disorganized repos and exploratory code, you can use the **agentic loop** as a "scientific code archaeologist":

```python
# Pseudocode for a scientific codebase organizer agent
tools = [
    "list_dir", "read_file", "grep_code", "run_test", 
    "git_log", "identify_dependencies", "extract_physics_invariants"
]

while not task_complete:
    response = llm.call(context + physics_prompt)
    
    if response.tool == "identify_dependencies":
        # Map which functions depend on which forcefields/integrators
        result = analyze_call_graph(response.params)
    elif response.tool == "extract_physics_invariants":
        # Find conservation laws, symmetries, boundary conditions
        result = parse_physics_semantics(response.params)
    elif response.tool == "run_test":
        # Quick smoke test to see if code still runs
        result = execute_quick_test(response.params)
```

**Key insight from Hermes**: Use **progressive disclosure** — keep a compact index of your physics modules in the stable prompt, load full `SKILL.md` files on demand when the agent needs to work on a specific forcefield or integrator. 

### **2. Cross-Platform Migration (Fortran → C++ → OpenCL → WebGPU)**

This is where the **subagent pattern** shines. From Hermes's architecture:

> "Subagents are processes, not threads — separate budgets, separate tool grants, separate SessionDB rows, summarized return channel." 

You can spawn parallel subagents for different migration paths:

```python
# Parent agent: "Migrate water forcefield to WebGPU"
parent_task = "Migrate TIP4P/2005 water model from Fortran to WebGPU compute shaders"

# Subagent 1: Analyze Fortran original
child_1 = delegate_task(
    task="Identify all array layouts, loop orders, and numerical precision in the Fortran source",
    toolset=["read_file", "grep", "parse_fortran"],
    budget=20
)

# Subagent 2: Draft C++ reference
child_2 = delegate_task(
    task="Create C++ reference implementation with identical numerical results",
    toolset=["read_file", "write_file", "compile", "run_test"],
    budget=30
)

# Subagent 3: Design OpenCL/WebGPU kernel
child_3 = delegate_task(
    task="Design compute shader structure for water force evaluation",
    toolset=["read_file", "write_file", "validate_wgsl"],
    budget=25
)

# Parent synthesizes results, checks parity
```

### **3. Systematic Multilevel Testing (Your Core Need)**

This is the most important application for you. Based on the **LARA-HPC** validation-driven approach from computational chemistry literature, you can build a **dry-run → validate → execute** loop: 

```python
class ScientificTestHarness:
    def __init__(self):
        self.test_levels = [
            "unit_test",           # Individual functions
            "parity_test",         # Bit-exact or tolerance match
            "invariant_test",      # Conservation laws, symmetries
            "qualitative_test",    # Physical common sense
            "convergence_test",    # Systematic mesh/k-point convergence
            "dft_reference"        # Expensive ground truth
        ]
    
    def agentic_loop(self, task):
        while True:
            response = llm.call(context + test_strategy_prompt)
            
            if response.is_text():
                return response.text
            
            if response.tool == "define_test":
                # Agent proposes a test based on physics invariants
                test = self.generate_test(response.params)
                result = self.dry_run(test)  # Cheap validation
                context.add_tool_result(result)
                
            elif response.tool == "run_parity_check":
                # Compare against reference implementation
                result = self.run_with_tolerance(
                    old_code=response.params["reference"],
                    new_code=response.params["target"],
                    tolerance=self.get_tolerance(response.params["level"])
                )
                context.add_tool_result(result)
                
            elif response.tool == "check_invariants":
                # Conservation laws, symmetry checks
                result = self.verify_physics_invariants(response.params)
                context.add_tool_result(result)
```

**Test hierarchy for your NaCl/water example**:

| Level | Test | Speed | Purpose |
|-------|------|-------|---------|
| 1 | **Unit**: Na+ coord. number | ~ms | Sanity check |
| 2 | **Parity**: Exact energy match vs old code | ~s | Regression detection |
| 3 | **Invariant**: Energy conservation (NVE) | ~s | Dynamics validity |
| 4 | **Qualitative**: O→Na+, H→Cl- orientation | ~s | Physical common sense |
| 5 | **Convergence**: Cutoff/radius convergence | ~min | Systematic error |
| 6 | **DFT Reference**: QM energy for dimer | ~hours | Ground truth |

The agent can **start with level 1-4** (fast, always run in CI), and only trigger level 5-6 when changes affect core algorithms.

### **4. Handling "Better but Different" Results**

Your key insight — that new methods may give *better* results, not identical ones — is handled by **qualitative test contracts**:

```python
def validate_migration(old_result, new_result, test_type):
    if test_type == "exact_parity":
        return abs(old_result - new_result) < 1e-10
    elif test_type == "conservation_law":
        # New method should conserve energy better, not worse
        return new_result.drift <= old_result.drift * 1.1
    elif test_type == "physical_common_sense":
        # O must be closer to Na+ than H is
        return new_result.d_O_Na < new_result.d_H_Na
    elif test_type == "convergence":
        # Must converge to same limit as old method
        return check_same_limit(old_result, new_result)
```

The agent loop can negotiate these: if parity fails but conservation improves, it flags for human review rather than auto-rejecting.

## Practical Implementation Path

Based on the **Hermes build-your-own checklist** (which is excellent for scientific use): 

**Phase 1 — The Loop (1-2 days)**
```python
class SciAgent:
    def run_conversation(self, messages, tools, budget=50):
        for i in range(budget):
            response = llm.call(messages, tools)
            if response.is_text():
                return response.text
            # Execute tool, append result
            result = self.execute_tool(response.tool_call)
            messages.append({"role": "tool", "content": result})
```

**Phase 2 — Scientific Tool Registry**
Register tools specific to your domain:
- `parse_physics_invariants`: Extract symmetries from code comments/math
- `dry_run_test`: Compile and run without full execution (like LARA-HPC's BigDFT dry-run)
- `parity_check`: Run both versions, compare with tolerance
- `convergence_sweep`: Systematically vary parameters, check monotonicity
- `git_snapshot`: Save state before destructive changes

**Phase 3 — Skills for Your Domain**
Create `SKILL.md` files for recurring tasks:
- `migrate_fortran_to_cpp`: "Preserve array layout, watch for 1-based vs 0-based indices, maintain double precision..."
- `validate_forcefield`: "Check O-Na+ distance < H-Na+ distance, energy minimum near sum of ionic radii..."
- `webgpu_compute_shader`: "Use workgroup size 64, shared memory for neighbor lists, avoid divergence..."

**Phase 4 — Memory & Learning**
- `MEMORY.md`: "User prefers OpenCL over CUDA for portability"
- `USER.md`: "Works on water/NaCl systems, cares about energy conservation, uses TIP4P/2005"
- Procedural memory: After successful migrations, the agent writes new skills

## Why This Matters for You Specifically

1. **Exploratory → Production**: Your disorganized repos are perfect for an agent with **heartbeat monitoring** — it can periodically scan for "orphaned" experiments, suggest consolidations, or flag when a branch hasn't been tested against main.

2. **Platform Migration**: The subagent pattern lets you parallelize the tedious work of reimplementing kernels while maintaining a single "integration agent" that checks parity.

3. **Test-Driven Refactoring**: Instead of writing tests after the fact, the agent can **infer tests from physics** — conservation laws, symmetries, known asymptotic behaviors — and make them fast enough to run on every commit.

4. **DFT Reference Integration**: The agent can manage the expensive ground-truth calculations as a separate tier, only running them when cheaper tests pass, and caching results.

The frameworks are mature enough to use today. I'd recommend starting with **Aider** (open-source, git-first, model-agnostic) for immediate codebase work, and graduating to **Hermes** or **OpenClaw** if you need the learning loop, subagents, and heartbeat monitoring for long-running HPC workflows.

---

# Claude Sonet 4.6

This is a rich and very relevant landscape for your situation. Let me synthesize what's out there and map it directly to your computational physics/chemistry use case.

## The Agentic Loop Landscape

### Core Pattern: The "Ralph Loop"

The key insight of the Ralph/RALF loop is storing state in the **external filesystem and git**, not in the context window. When an agent writes code, it's saved to files; when it commits, history accumulates. Context can reset at any time, but the filesystem and git are persistent. Each new agent starts "fresh" while fully inheriting the previous agent's results. This is exactly what you need for massive codebases.

The loop runs like this each iteration: implement → typecheck → run tests → lint → commit only after all checks pass. This ensures the agent doesn't submit destructive code.

### The Three Main Harnesses

**Claude Code** is your daily driver at your desk — a coding assistant that lives next to your code, in your terminal, designed for work where you want to be in the loop on every decision. **Hermes** and **OpenClaw** handle scheduled tasks on the go, handling background work while you're away from the keyboard.

Hermes Agent is an open-source AI agent runtime that gives a language model tools, memory, skills, files, terminal access, web access, scheduling, and communication channels. It supports provider flexibility across OpenRouter, Anthropic, OpenAI, local/custom endpoints.

## How This Maps to Your Specific Problems

### 1. Cross-platform porting (Fortran → C++ → OpenCL → WebGPU)

This is the **killer use case** for Ralph-style loops. The loop works well when: success is easy to validate (tests are green), failure is acceptable (complex cases get flagged and deferred), scope is tightly bounded per iteration, and lessons learned carry forward.

For your porting work, the loop would be: port a function → run parity tests → commit if pass → flag if qualitatively wrong → iterate.

### 2. Your multilevel test philosophy

This is exactly right and well-supported. Without tests, the agent moves forward without actually checking if things work. Your instinct about **qualitative correctness** (oxygen sits near Na⁺, H toward Cl⁻, distances near sum of atomic radii, energy conservation) is actually *better* for agentic loops than exact numerical parity — because it gives the agent a fast, robust oracle.

The levels you'd want:
- **Level 0** (seconds): Compile, dimension checks, conservation laws, symmetry checks
- **Level 1** (minutes): Qualitative geometry (your NaCl/water example), sign checks, monotonicity
- **Level 2** (hours): Quantitative regression vs DFT reference, RDF, PMF curves
- **Level 3** (manual/gated): Physics review — "does this make sense?"

### 3. Disorganized repos and data

The proposed human-led AI agent team has the **human as lead Architect** who designs features, breaks projects into manageable testable modules, and creates gold standard test sets. **Thinking agents** (expensive reasoning models) help understand the codebase and plan. **Coding agents** (fast cheaper models) perform actual coding and testing.

For your exploratory/unfinished repos, the pattern would be: Thinking Agent reads the repo → writes a `SUMMARY.md` + task list → Coding Agent works through it iteratively.

### 4. Agentic loop for test generation itself

A **Test Coverage Loop** works like this: point the agent at your coverage metrics, it finds uncovered lines, writes tests, and iterates until coverage meets your goal. You could seed it with your physical constraints (conservation laws, symmetry, known limits) and let it generate the test suite.

## Concrete Starting Point for You

Given you're a computational physicist with Fortran/C++/OpenCL/WebGPU codebases, I'd suggest:

1. **Claude Code + Ralph Loop** for active porting tasks — you stay in the loop, review physics
2. **Hermes** for background repo organization — runs overnight, generates `README.md`, tags unfinished work, finds duplicated implementations
3. A `PHYSICS_CONSTRAINTS.md` in each repo — your qualitative test oracle (conservation laws, symmetry, chemical common sense) — which becomes the agent's first-principles test spec

Want me to draft a `PHYSICS_CONSTRAINTS.md` template for your forcefield work, or sketch out what a Ralph-style loop for a Fortran→C++ port with parity checks would look like in practice?

---

# ChatGPT (5.5?)

I think you are thinking in exactly the direction where agentic loops are most valuable. Many demos like OpenClaw or Hermes agents are focused on "make a website" or "solve LeetCode", but for computational science the more interesting application is **long-running autonomous scientific software engineering**.

Your description also highlights something important: the AI should not just write code. It needs a **scientific verification harness**. In your domain, the tests are not just `assert x == 5`; they are physical invariants, reference comparisons, qualitative behavior, and multi-scale validation.

Let's first look at the architecture of agentic loops, then map them to your work.

---

# 1. Simple coding agent loop (OpenClaw style)

The simplest agent loop is:

```
while not done:

    understand task

    inspect codebase

    propose modification

    edit files

    run tests

    inspect failures

    reason about errors

    iterate
```

In pseudocode:

```python
while budget_left:

    state = read_repository()

    plan = LLM("""
    Current state:
    {state}

    Goal:
    {goal}

    What is the next smallest useful step?
    """)

    execute(plan.actions)

    result = run_tests()

    if result.success:
        commit()
        break

    feedback = analyze_logs(result)
```

This is already extremely useful for:

* refactoring C++,
* porting Fortran → C++,
* adding unit tests,
* fixing build systems,
* improving documentation.

# 2. More advanced "scientific coding agent"

For your work I would add a much richer loop:

```
              +----------------+
              |    Scientist   |
              | high-level goal|
              +--------+-------+
                       |
                       v
             +---------+---------+
             |  Planner agent    |
             +---------+---------+
                       |
             creates subtasks
                       |
      +----------------+----------------+
      |                                 |
      v                                 v

Code agent                    Knowledge agent

read code                     read papers
add functions                 inspect notes
refactor                      classify data

      \                                 /
       \                               /
        -------- Verification ----------
                    |
            run physical tests
                    |
          update confidence model
```

The critical thing is the **verification agent**.

# 3. Your most important missing piece: hierarchical tests

I actually think this is the biggest opportunity in your repositories.

You need tests at several levels:

## Level 0 — numerical sanity

Very cheap.

Example:

```text
Energy finite?
Forces finite?
No NaN?
No atoms exploded?
```

## Level 1 — symmetry tests

Examples:

```
Translation invariance:

E(R + Δ) = E(R)

Rotation invariance:

E(R) = E(U R)
```

For forces:

```
sum(F_i) = 0
```

for isolated systems.

## Level 2 — physical qualitative tests

This is exactly your NaCl-water example.

Not:

```
E = -0.3567 eV
```

but:

```
Water over Na+

Expected:
    O atom closest to Na
    distance 2-3 Å
    dipole points toward ion
    stable minimum exists
```

Another:

```
Water over Cl-

Expected:
    H atoms point toward Cl
    H-Cl distance ~2-3 Å
```

The test returns a score:

```
orientation score: 0.95
distance score: 0.83
binding energy sign: PASS
```

## Level 3 — reference comparisons

Against old implementation:

```
Fortran
    |
    | same geometry
    |
    v
C++ implementation

compare:
    energy difference < tolerance
    force cosine similarity > 0.99
```

This is extremely useful for ports.

## Level 4 — high-level scientific behavior

This is where many AI coding systems fail because they don't understand science.

Examples:

AFM simulation:

```
Does the image have expected symmetry?
Are bright spots above atoms?
Does contrast invert at expected height?
```

Molecular dynamics:

```
Energy conservation
Temperature distribution
Diffusion coefficient order of magnitude
```

Electronic structure:

```
Correct ordering of orbitals
Charge conservation
Symmetry of wavefunctions
```

# 4. The loop should optimize confidence, not "pass/fail"

I would actually avoid the standard software engineering mindset:

```
all tests green = good
```

For scientific software it is often wrong.

A better model is:

```
Scientific confidence score

Geometry tests:      0.92
Symmetry tests:      1.00
Reference parity:    0.97
Energy accuracy:     0.80
Performance:         0.95

Total confidence: 0.93
```

because, as you said, a new method can be **better** than the old one.

# 5. The agent should be allowed to build its own tests

This is the part where AI becomes really powerful.

Imagine giving it:

```
Goal:
Make the OpenCL version of MMFF equivalent
to the C++ version.
```

The agent does:

```
1. Read C++ implementation.

2. Identify all force terms:
    - Coulomb
    - Lennard-Jones
    - hydrogen bonds
    - pi interactions

3. Generate minimal molecular examples:
    - water dimer
    - NaCl cluster
    - benzene dimer

4. Create test inputs.

5. Run C++ reference.

6. Run OpenCL implementation.

7. Compare results.

8. Investigate discrepancies.
```

This is very close to what companies are trying to do with coding agents, but your domain has a huge advantage: physics gives you a lot of constraints.

# 6. Repository archaeology agent (I think this would help you enormously)

Your situation:

> massive codebases, disorganized repositories, old projects, notes, data

This is almost a perfect agent task.

Loop:

```
for each directory:

    classify content:
        code?
        paper?
        data?
        old experiment?

    extract metadata:
        language
        purpose
        dependencies
        status

    create summary.md
```

Then a second pass:

```
find duplicates:
    "this OpenCL kernel and this C++ file
     implement the same potential"

find dead code

find undocumented algorithms

build dependency graph
```

The final product would be something like:

```
FireCore/
    CODEMAP.md

Modules:
    MMFF
      status: active
      tests: 34
      coverage: 80%
      GPU parity: incomplete

    DFT
      status: mature
      Fortran origin
      C++ rewrite exists

Old experiments:
    water_surface_scan_2018/
       likely obsolete
```

# 7. I would actually design a "scientific CI"

A normal CI says:

```
Does it compile?
Do unit tests pass?
```

Your CI should say:

```
Scientific health report

Build:
    PASS

Performance:
    +12%

Physical invariants:
    PASS

Water/NaCl tests:
    PASS (0.93)

AFM image similarity:
    PASS (0.87)

Reference agreement:
    WARNING
    H-bond force differs 5%
```

# My impression

Based on the type of work you do (large old Fortran/C++/OpenCL codes, many exploratory branches, lots of data), I suspect the highest return is **not** an "AI programmer" that writes new algorithms. It is an **AI scientific software archaeologist + verification engineer**.

I would probably start with a repository-wide protocol such as:

```
AGENTS.md
CODEMAP.md
TESTS.md
PHYSICS_TESTS.md
EXPERIMENTS.md
```

where the agent continuously expands those files.

Given your previous work on AutoCrunchCoder and code dependency extraction, you are actually already halfway toward building exactly this kind of agent. The missing component is not code generation; it is the **closed loop of inspection → hypothesis → experiment → verification → documentation**.

---

# USER


look in @Test_system_for_agentic_loops.md 
and sumarize it

where the models agree, disagree, where some has inovative idea, what seems to you as a best path forward based on these information and your own judgement

---

# Kimi 2.6

# Summary of [Test_system_for_agentic_loops.md](cci:7://file:///home/prokop/git/FireCore/doc/DevNotes/Test_system_for_agentic_loops.md:0:0-0:0)

This document captures six AI model responses (Gemini, DeepSeek, Mistral, Kimi, Claude, ChatGPT) to your question about building agentic loops for computational physics/chemistry — specifically for cross-platform porting, multilevel scientific testing, and organizing disorganized repositories.

## Where the Models Agree

- **Core loop pattern**: All converge on the same fundamental cycle — *Plan/Reason → Act (code/edit) → Observe/Test → Reflect/Fix* — with iteration until success or budget exhaustion.
- **Scientific testing is hierarchical**: Universal consensus on multilevel validation — from fast syntactic/compilation checks (L0), through parity/regression (L1-L2), to physical invariants and qualitative "common sense" (L3-L4), up to expensive DFT/reference comparisons (L5). All agree exact numerical parity is often the wrong goal; qualitative physical correctness (energy conservation, symmetry, orientation) is more robust.
- **Cross-platform porting is the killer use case**: All identify Fortran → C++ → OpenCL/WebGPU as the ideal application for these loops.
- **Repository archaeology matters**: All recognize your disorganized repos as a solvable problem via automated inventory, metadata extraction, and duplicate/dead-code detection.
- **Tooling recommendations overlap**: OpenClaw, Hermes, Claude Code/Aider, and ralph-loop are mentioned consistently as viable starting points.

## Where the Models Disagree

- **Framework choice**:
  - **Gemini & ChatGPT**: Favor minimal custom Python harnesses over heavy frameworks; build your own lightweight loop with direct API calls.
  - **DeepSeek & Kimi**: Emphasize feature-rich orchestrators (Hermes for learning/memory, SwarmWeaver for parallelism, Loop-Harness for 2-agent verification).
  - **Claude**: Recommends a hybrid — Claude Code for daily active work, Hermes/OpenClaw for background/scheduled tasks.
  - **Mistral**: Advocates framework-specific adoption (ralph-loop for task focus, Hermes for persistent learning).
- **State persistence**:
  - **Claude** insists state lives in *filesystem + git*, not context window (Ralph Loop philosophy).
  - **Kimi & Hermes** emphasize persistent memory/skills databases (SQLite, Markdown memory, procedural learning).
- **Test philosophy nuance**:
  - **ChatGPT** is the most distinct: explicitly rejects binary pass/fail. Advocates for a *continuous confidence score* model (e.g., "Total confidence: 0.93") because new methods can legitimately outperform old ones.
  - **Gemini & Mistral** stick closer to traditional layered pass/fail with tolerance bands.

## Innovative Ideas

| Model | Innovation |
|-------|-----------|
| **Gemini** | "Physical Oracle" concept — a standalone Python `validate.py` that evaluates qualitative & physical sanity separately from numerical parity. Also the "Map-Reduce" chaos sweeper loop that indexes repos by physical system, language, and completeness before touching code. |
| **DeepSeek** | Comparative framework table and categorization into *Foundational Loops* / *Orchestrators* / *Testing Agents* / *Porting Agents*. Introduces "cross-family validation" (independent LLMs judging code quality). |
| **Mistral** | Formal 5-level test hierarchy (L0-L5) with explicit speed/purpose mapping. Also integrates ReAct, Reflexion, Plan-and-Execute, and Agentic RAG as distinct patterns with specific best-use cases. |
| **Kimi** | Hermes-style "progressive disclosure" — keep a compact physics index in the stable prompt, load full `SKILL.md` files on demand. Subagent pattern for parallel migration paths (Fortran analyzer, C++ drafter, WebGPU designer running concurrently). |
| **Claude** | Human-as-Architect / Thinking Agent / Coding Agent team structure. "Test Coverage Loop" — agent points itself at coverage metrics and writes tests until target is met. |
| **ChatGPT** | **Confidence-score model** instead of binary pass/fail. "Scientific CI" — a health report with graded scores per category (physical invariants, performance, reference agreement). Most explicit about the agent *generating its own tests* from physics invariants. |

## Best Path Forward (My Assessment)

Given your actual workflow in FireCore (Fortran DFT core, C++ apps, OpenCL/WebGPU ports, parity validation, disorganized exploratory code), here is what I would prioritize:

### 1. Start with the Minimal Loop, Not a Framework
Gemini and ChatGPT are correct here. Don't adopt Hermes or LangChain yet. Build a **~50-line Python loop** that:
- Reads a task from a markdown file
- Calls an LLM API with a system prompt containing your physics conventions
- Executes the returned shell/Python commands
- Captures stdout/stderr fully (per your rule: *never hide output*)
- Appends results back to the conversation
- Repeats until a test passes or budget expires

This is immediately useful and avoids framework lock-in.

### 2. Invest Heavily in the Test Harness First
All models agree this is the missing piece. I recommend ChatGPT's **confidence-score model** over binary pass/fail, because FireCore already has cases where the GPU implementation is *more* correct than the CPU reference (e.g., PME parity fixes). A test should report:

```
parity:      0.97 (tolerance 1e-5)
conservation: 1.00
symmetry:    0.99
physical:    0.93 (O→Na+ distance OK, but H-Cl slightly off)
```

This maps directly to your `AGENTS.md` rules about physical correctness and avoids false rejections of legitimate improvements.

### 3. Implement the "Repo Archaeologist" First
Before any porting, run the inventory loop Claude and ChatGPT describe. FireCore already has `CODEMAP.md` and `DOX_AGENTS.md` — extend this with an automated sweep that:
- Scans `tests/` and `cpp/` and `pyBall/` for untested modules
- Identifies duplicate implementations (e.g., same forcefield in C++ and OpenCL)
- Flags files with `TODO`/`FIXME` or debug exports
- Generates a prioritized task list

This gives you immediate value without risking code breakage.

### 4. Use Subagents for Porting (Kimi/Hermes pattern)
When you do port Fortran → OpenCL or C++ → WebGPU, use the subagent pattern:
- **Analyzer subagent**: reads Fortran, extracts loop order, array layout, precision
- **Reference subagent**: writes a minimal C++ test harness with the same inputs/outputs
- **Port subagent**: writes the GPU kernel
- **Verifier subagent**: runs parity and physical checks, reports confidence scores

Run these in sequence, not parallel, to keep deterministic — but delegate each phase as a bounded task with its own budget.

### 5. Adopt Claude Code / Aider for Daily Work, Custom Loop for Batch Work
Claude's hybrid recommendation is practical. Use **Aider** (git-first, model-agnostic, tree-sitter repo map) for active editing sessions where you review every change. Use your custom loop for overnight batch tasks (inventory sweeps, regression testing, documentation updates).

### 6. Keep State in Git + Markdown, Not Context Windows
Follow Claude's Ralph Loop principle. Every agent output should be a file edit or git commit. Conversation state should be reconstructible from `AGENTS.md`, `CODEMAP.md`, and commit history — not from an LLM's memory. This is the only way to survive context-window limits on a codebase the size of FireCore.

### Concrete First Step
1. Write `tests/agentic_loop.py` — a minimal loop that reads `task.md`, calls an LLM, runs commands, logs output.
2. Seed it with one bounded task: *"Inventory all `.cl` files in `cpp/common/OpenCL/` and `pyBall/FireballOCL/cl/` and report which kernels lack corresponding CPU reference tests."*
3. Run it, review the output, commit the result.

This gives you a working loop and actionable repo intelligence within an hour, with zero framework dependencies.


---

# USER

OK, now we must focus on most important part of these agentic loops which is define the goal and test if the goal was achived.
I can imagine several typical situations for which I would like to use such agentic loops
1) reimplementation and parity checking
2) performance optimization

I often reimplement same computiing method into different code (e.g. fotran to C++, C++ to OpenCL or rust, python to javascript etc), Or I have some 3rd-party opesource software like DFTB+, Grimme xTB, Fireball, SIESTA, and I try to exctract the minimal computational core and reimplement it in my software. Or I want to implement new method e.g. different formulation of forcefield, approximate quantum mechanical method, and I want to make sure it give qualitatively reaosnable results or how much it differen form some accurate reference. My way is typically imporving speed (often on GPU, using OpenCL kernels in single-precission) rather than imporving accuracy.

Distinguish human-vosial test and LLM-test base test and automatic test
1) what can be checked automatically (NaNs, RMSE parity) [cheap]
2) what need inteligent analysis - typically LLM can analyse text output of the code [medium cost]
3) human has better visual intiutun and can debug GUI and review .xyz files (trajectroy), .cube files in VESTA, and .png or .svg or .html p5.js plots [costly]

we should do as much as possible by atumatic check [cheap],  then by LLM [medium], and then by huan only rarely [costly]

we should assume thet the debug cycle:   edit_code -> compile -> test_run -> -review_result -> reflect -> plan next move
that this loop will be exectuted mutiple times unti the goal is achived, therefore economy of attention is important

I can imagine there can be multiple levels of testing
1) first it must compile errors and not crash (no segfaults, OpenCL device freeze etc)
2) then we should check if it gives reasonable numbers (no NaN, no infinity, also the number can have some reasonable natural range). User should either provide some reasonable range, or agent should infer by his physical intution what are the physical ranges involved

3) check mapping - very ofthn the problem is in topology, that the neihbor map indexing, is somehow inconsistent, we should either automatically check if with reference (that angles point to same atoms, that we do not have uninitialized arrays of neighborlist, and angles description, whe should have proper padding of empty slost e.g. by -1

5) rigorous fixed accuracy test

it is rather easy to compute energy and force for single point, but unles we copy exactly same method to differenc cde this typically is difficult to jude, the force or energy can be different between methods yet they mayconverge to same physically meaningfull geometry

If we reimplement the same method we can also accurately compare the rigi scans (wehere we move some atom, scan bond lengh, angle, dihedral, scan molecule rigidly over substrate etc) and plot the profile. LLM can compute RMSE, and see numbers in text but it should also produce some plots which can be easly and quickly reviwd by human, these debug artifaces must be well described and documented and LLM should tell use where exactly it is

4) we should relax the system and see if it produce reasonable geometry

then we use some simplest posible test examples and check if ti give some reaonsalbe qualitative results. Typically forst for QM is H2 molecule with 2 s orbitlas, then N2 molecule with sp basis, then maybe CO molecule (two different atoms, polarization), or H2O. in classical forcefielw swe typically start frm H2O, NH3 (tetrahedral), then formaldehyde and formic acid, oxalate, benzene or pyridine, which have pi-orbitals and sp2 hybridization (need either torsions and dihedrasl in UFF or pi orbitlans and pi-sigma orthogonalization pi-pi alignement in our. For non covalent forcefiled we can check carboxylic acid dimer. For interaction with surface we typically use fisrt NaCl as it is cubic, we start with scanning single atom/ion with selected LJ/mosre parameter and charge . We can have qualitative results like H2O absorb on Na+ site by oxygen, and then on Cl side it turn by hydrogen.

these qualitative tests are typically less automatic than the nex rigorous tests but are esential for undestating what happens.

we should ditinguis single point clauclation - where the forces are going, and optimization with forces

5) then comes performance testing - we can try different version of kernel, memory layour, strategy of reading data from global memory, workgrou psize, different strategy to write to global memory (coalesced) using the assembly strategy for recoild forces to avoid atomic writes or write collision from different threads, or to avoid need for costly GPU reduction

now plese elaborate in these

---


# Kimi  2.6 / devin

Here is a structured elaboration of the testing hierarchy, mapped to your specific FireCore workflows and the principle of **economy of attention**.

# Agentic Testing Hierarchy for FireCore

## Principle: Escalating Cost of Attention

```
Automatic (cheap, fast, run every iteration)
    ↓
LLM Analysis (medium, parse text/logs/plot descriptions)
    ↓
Human Review (expensive, visual, 3D, GUI — minimize)
```

Every artifact must be **self-describing**: filename, plot title, and caption should make the issue obvious without reading code.

## Level 0: Compile & Runtime Sanity (Automatic)

| Check | Method | Example in FireCore |
|-------|--------|---------------------|
| Compiles | `subprocess.run(["make","-f","makefile"])` or `run.sh` | `tests/tUFF/run.sh` |
| No segfault | Run with `timeout 30` | MD for 10 steps on H2O |
| No NaN/Inf | `np.isnan(arr).any()` / `np.isinf(arr).any()` | Forces, energies, densities |
| No GPU freeze | OpenCL event timeout / `CL_COMPLETE` check | `pyBall/FireballOCL` kernels |
| Output files exist | `os.path.exists()` | `.xyz`, `.npy`, `.png` |

**Fail behavior**: Hard stop. Log stderr verbatim into the agent's context.

## Level 1: Reasonable Number Ranges (Automatic + LLM Inference)

Two sub-cases: user-provided bounds, or LLM-inferred physical bounds.

### Automatic (user-provided)
```python
def check_physical_bounds(energy, forces, coords):
    assert energy > -1e6 and energy < 1e6          # eV scale
    assert np.all(np.abs(forces) < 1e4)             # eV/Å, non-exploded
    assert np.all(np.abs(coords) < 1e3)             # Å scale
```

### LLM-inferred (when user doesn't specify)
The agent reads the method name and infers:
- **Lennard-Jones ε**: typically 0.01–0.5 kcal/mol → expect well depth ~meV
- **Coulomb Na⁺–Cl⁻**: ~2.4 Å, energy ~−120 kcal/mol in vacuum
- **H₂ bond**: 0.74 Å, binding ~4.5 eV
- **Water dipole**: ~1.85 D → partial charges ~±0.4e

If values violate these, the agent flags immediately without human review.

## Level 2: Topology & Mapping Integrity (Automatic + Reference Comparison)

This is where most FireCore parity bugs live (per your memory: neighCell zero-shift, fneigh buffer underallocation, uninitialized arrays).

### Automatic checks against reference

| Check | How | FireCore Example |
|-------|-----|------------------|
| Neighbor count consistency | `len(neighs_new) == len(neighs_ref)` per atom | `MMFF.py` neighCell rebuild |
| Self-neighbor present | `any(nj==iatom+1 and nb==0 for ...)` | Sparse HS blocks |
| Angle atom indices valid | `all(0 <= idx < natoms for idx in angle_atoms)` | UFF angle/dihedral arrays |
| Padding correct | `assert padded_unused == -1 or 0` | `nang_max`, `ndih_max` padding |
| Array shapes match | `arr_new.shape == arr_ref.shape` | `fneigh`, `avel`, `cvf` buffers |
| Index mapping bijective | `set(idxs_new) == set(idxs_ref)` | After node reordering |

### Topology diff report
The agent prints a concise text diff:
```
Atom 7: ref_neighs=[3,5,8]  new_neighs=[3,5,8,8]  ← DUPLICATE
Angle 12: ref_atoms=(2,5,9)  new_atoms=(2,5,10)  ← INDEX SHIFT
```
This is pure text — LLM can parse and diagnose without human.

## Level 3: Single-Point Quantitative Parity (Automatic / LLM)

For **exact reimplementation** (Fortran→C++, C++→OpenCL).

### What to compare
| Quantity | Tolerance | FireCore Path |
|----------|-----------|---------------|
| Energy | `rtol=1e-5` (double), `1e-3` (single) | `firecore_get_energy()` |
| Forces per atom | `atol=1e-4 eV/Å` | `firecore_get_forces()` |
| Hamiltonian blocks | `atol=1e-6` per block | `get_HS_sparse()` |
| Density matrix | `atol=1e-5` | `get_rho()` |
| AO grid projection | `atol=1e-4` | `STM_utils.py` vs `Grid.cl` |

### Rigorous scan parity
For forces/energy that depend sensitively on geometry, a **rigid scan** is more discriminating than a single point:

```python
def rigid_scan_parity(ref_func, new_func, geometry, scan_dim, steps):
    # Scan bond length, angle, or molecule-surface distance
    # Return energies, forces along scan coordinate
```

**Example scans for FireCore:**
- **H–H distance**: 0.5–2.0 Å (Morss/DFT)
- **H₂O bend angle**: 90–120° (forcefield angular term)
- **Na⁺–H₂O distance**: 2.0–5.0 Å (nonbonded LJ+Q)
- **Benzene over NaCl(001)**: 2.0–6.0 Å z-scan (π–surface interaction)

**Artifact**: `scan_parity_<name>.png` with two curves (CPU dotted, GPU solid — per your existing convention in `tests/tEFF/run_scan_batch.sh`). Caption in plot: max |dE|, RMSE, max |dF|.

**Automatic**: RMSE, max abs diff.  
**LLM**: Reads the RMSE, decides if within tolerance.  
**Human**: Glances at plot shape (smooth? symmetric? minimum at expected position?).

## Level 4: Qualitative Physical Sanity (Automatic + LLM + Sparse Human)

These tests validate that the method produces *physically meaningful* behavior, even when exact parity is impossible (new approximation, different formulation).

### Tiered qualitative test suite for FireCore

#### QM / Tight-Binding
| System | Expected Behavior | Automatic Check | LLM Review |
|--------|-------------------|-----------------|------------|
| **H₂** (1s basis) | Bonding σ orbital, E_bind ~4.5 eV | `E_min < 0`, `dE/dr=0` at ~0.74 Å | Plot E(r) — minimum visible? |
| **N₂** (sp basis) | Triple bond, short r_eq ~1.1 Å | `r_eq < 1.2` Å | Verify orbital ordering |
| **CO** (polarization) | Dipole moment ~0.1 D (C→O) | `mu > 0` (sign) | Check charge asymmetry |
| **H₂O** (sp³, bent) | Angle ~104.5°, dipole ~1.85 D | `angle > 100°`, `angle < 110°` | Plot geometry |
| **NH₃** (tetrahedral) | Pyramidal, lone pair | `angle_HNH ~107°` | Check inversion barrier |

#### Classical Forcefields (UFF / MMFF)
| System | Expected Behavior | Automatic Check |
|--------|-------------------|-----------------|
| **H₂O** | Bent, O–H ~0.96 Å, H–O–H ~104.5° | Relax → geometry bounds |
| **NH₃** | Trigonal pyramidal, H–N–H ~107° | Relax → angle bounds |
| **Formaldehyde** (H₂CO) | Planar C=O, sp² carbon | Dihedral H–C–O–H ≈ 0° |
| **Formic acid** | Planar, strong H-bond dimer | Dimer E_bind > monomer |
| **Benzene / pyridine** | Planar ring, π-system | All ring dihedrals ~0° or 180° |
| **Carboxylic acid dimer** | Two H-bonds, cyclic dimer | O–H···O distance 1.6–1.8 Å |

#### Surface / Nonbonded
| System | Expected Behavior | Automatic Check | Human Review |
|--------|-------------------|-----------------|--------------|
| **NaCl(001)** | Cubic, alternating Na⁺/Cl⁻ | Coordination check | Visual in VESTA |
| **H₂O on Na⁺** | O closest to Na⁺, dipole pointing away | `r(Na–O) < r(Na–H)` | `.xyz` snapshot |
| **H₂O on Cl⁻** | H atoms point toward Cl⁻ | `r(Cl–H) < r(Cl–O)` | `.xyz` snapshot |
| **Rare gas on NaCl** | Physisorption, no chemisorption | E_bind < 0.5 eV | Plot E(z) |

### How the agent validates qualitatively

**Automatic**: Hard assertions on geometry ranges, energy sign, distance ordering.  
**LLM**: Reads the geometry output and confirms physical narrative:
```
"Water relaxes with O at (0.1, 0.1, 2.3) Å above Na⁺ site. 
Na–O distance = 2.31 Å, Na–H distances = 3.12 Å, 3.15 Å. 
Conclusion: O is closer to Na⁺ than H — physical orientation CORRECT."
```
**Human**: Only if LLM flags ambiguity (e.g., "O and H both 2.5 Å from Na — cannot determine orientation").

## Level 5: Relaxation & Dynamics (Automatic + LLM + Human for Trajectory)

Single-point checks are insufficient. The system must find stable minima and not explode.

### Tests
| Test | What | Automatic | LLM | Human |
|------|------|-----------|-----|-------|
| **Geometry relaxation** | Minimize 100 steps, check convergence | `max|F| < 0.01 eV/Å`, `E_final < E_initial` | Read convergence log | Review `.xyz` if NaN |
| **Conservation (NVE)** | 1000-step MD, no thermostat | `ΔE/E < 1e-4` | Check drift plot | Review if drift > 1e-3 |
| **Symmetry preservation** | Relax symmetric input | RMSD from input symmetry < 0.1 Å | — | VESTA if broken |
| **No cluster dissociation** | H₂O dimer → stays dimer | `max_pair_dist < 5 Å` | — | `.xyz` if split |

### Artifact: `relax_parity_<system>.png`
Plot E vs step for CPU (dotted) and GPU (solid). Same color, different line style — your existing convention.

## Level 6: Performance Optimization (Automatic + LLM Analysis)

When the goal is speed, not parity.

### What to measure
| Metric | How | Automatic |
|--------|-----|-----------|
| Kernel time | OpenCL event profiling | `clGetEventProfilingInfo` |
| Memory bandwidth | Bytes read/written per kernel | Compute from buffer sizes |
| Occupancy | Workgroup size, register usage | `clGetKernelWorkGroupInfo` |
| Host-device transfer time | `queue.finish()` timing | Separate upload/compute/download |
| End-to-end throughput | Steps/second for N atoms | `N_steps / total_time` |

### Optimization strategy loop
The agent iterates:
1. **Profile**: Identify bottleneck (memory-bound? compute-bound? atomic collision?)
2. **Hypothesize**: "Coalesced global reads" or "avoid atomics via recoil assembly"
3. **Implement**: Edit kernel memory layout or force accumulation strategy
4. **Validate**: Run same test suite (Levels 0–5) to ensure correctness preserved
5. **Measure**: Compare throughput, reject if slower or incorrect

### Performance artifact: `benchmark_<kernel>_<config>.txt`
```
System: H2O_64box
Kernel: evalMMFFf4_local
Workgroup: 64
Shared mem: 8192 B
Avoid atomics: True
Time/step: 0.023 ms
Throughput: 2.8M atom-steps/s
Speedup vs baseline: 1.45x
Parity max|dF|: 2.3e-5 eV/Å  ← from Level 3
```

## The Complete Debug Cycle

For one agentic iteration (e.g., porting Fortran LJ kernel to OpenCL):

```
1. EDIT: Agent writes .cl kernel
2. COMPILE: Run build script → log (Level 0)
   └─ FAIL → feed compiler error to LLM, goto 1
3. TEST RUN: Execute on H2O test (Level 1)
   └─ NaN/Inf → feed output to LLM, goto 1
4. TOPOLOGY CHECK: Compare neigh list vs CPU ref (Level 2)
   └─ Mismatch → feed diff to LLM, goto 1
5. PARITY SCAN: Rigid z-scan Na⁺–H₂O (Level 3)
   └─ RMSE > 1e-3 → plot diff, feed to LLM, goto 1
6. QUALITATIVE: Relax H₂O on NaCl(001) (Level 4)
   └─ O→Na⁺ check FAIL → feed geometry to LLM, goto 1
7. CONSERVATION: 1000-step NVE (Level 5)
   └─ ΔE > 1e-4 → feed energy drift plot to LLM, goto 1
8. PERFORMANCE: Benchmark vs CPU (Level 6, if goal)
   └─ Slower → profile, feed metrics to LLM, goto 1
9. COMMIT: All levels pass → git commit
```

Each level gates the next. No human is involved until Level 4 ambiguity or Level 5 trajectory review.

## Artifact Naming Convention

All artifacts must be discoverable without reading code:

```
test_H2O_sp/
  ├── scan_bond_HH.png              # Level 3: H–H rigid scan
  ├── scan_bend_HOH.png             # Level 3: H–O–H angle scan
  ├── relax_NaCl_001.xyz            # Level 4: relaxed geometry
  ├── relax_NaCl_001.png            # Level 4: before/after overlay
  ├── md_NVE_1000_E.png             # Level 5: energy vs time
  ├── md_NVE_1000_T.png             # Level 5: temperature vs time
  ├── benchmark_wg32.txt            # Level 6: profiling
  └── REPORT.md                     # LLM summary of all checks
```

`REPORT.md` is generated by the agent for human review:
```markdown
## Test: LJ_kernel_OpenCL_port
### Level 0: PASS (compile, no crash)
### Level 1: PASS (E=-0.234 eV, Fmax=1.2e-2, all finite)
### Level 2: PASS (neigh count match: 34/34 atoms)
### Level 3: RMSE scan = 2.1e-4 eV (PASS)
### Level 4: H2O on Na+ → O closest (PASS), r_NaO=2.31Å
### Level 5: NVE drift = 3e-5 (PASS)
### Level 6: 1.45x speedup vs CPU
**STATUS: READY FOR MERGE**
```

## What You Should Do Next

1. **Pick one reference implementation** (e.g., existing C++ MMFF) and one target (e.g., OpenCL `evalMMFFf4_local`).
2. **Write the 6-level test harness as a Python script** in `tests/tMMFF/agentic_validation.py` that executes Levels 0–3 automatically, generates plots, and produces `REPORT.md`.
3. **Run it once manually** to establish the baseline (CPU reference outputs).
4. **Then hand the loop to the agent**: "Port this function. The validation script is at `tests/tMMFF/agentic_validation.py`. Run it after each edit. Do not ask me until Level 4 fails."

This gives you an immediate, concrete loop that respects attention economy and maps directly to your existing FireCore testing conventions.

---

# Gemini

Building a resilient, cost-effective testing pyramid within an agentic loop requires balancing the rigorous determinism of scientific computing with the creative, iterative reasoning of LLMs.

Because you are focusing on high-speed GPU architectures (OpenCL/WebGPU) where single-precision floats, memory alignment, and device freezes are frequent hurdles, your agentic harness must prioritize **fast failures first** to protect your expensive loop iterations and human focus.

Here is an architectural elaboration on how to structure the testing hierarchy, map your qualitative test cases, and manage the economics of debug cycles.

## 1. The Economics of the Debug Cycle

To prevent the agent from burning API costs or wasting your time reviewing broken artifacts, we categorize the evaluation stack into three strict financial and temporal tiers:

| Tier | Type | Speed / Cost | Primary Responsibility |
| --- | --- | --- | --- |
| **Tier 1** | **Automatic Code Assertions** | < Seconds / Millicents | Syntax, device execution, NaNs/Infs, array boundaries, bitwise parity ($\Delta < \epsilon$). |
| **Tier 2** | **LLM-Driven Log Analysis** | ~ Seconds / Cents | Evaluating semantic text output, verifying physical boundaries, reviewing compiler warnings, assessing structural matrices ($F_{ij} = -F_{ji}$). |
| **Tier 3** | **Human-Visual Intuition** | ~ Minutes / High Cost | Spotting qualitative defects in geometries (`.xyz`), potential energy surface curves (`.png`/`.svg`), and visual web dashboards. |

### Rule of the Loop:

The agent must **never** escalate to a higher tier if a lower tier fails. If an OpenCL kernel causes a device crash (Tier 1), the loop immediately resets, adapts the source code, and re-runs. The human is only alerted when all Tier 1 and Tier 2 criteria pass.

## 2. Elaborated Multi-Level Testing Strategy

### Level 1: Execution & Runtime Integrity (Tier 1)

Before counting a single physical value, you must guarantee code stability. For OpenCL/WebGPU, this is notoriously fragile.

* **Compilation Verification:** Intercept compiler error logs. The loop reads the error line and injects it back to the prompt for fixing.
* **Crash & Freeze Guards:** Wrap execution in native OS timeouts. An infinite loop or GPU workgroup deadlock in OpenCL can freeze the display driver. The agentic harness must track subprocess execution time and treat a timeout as a catastrophic failure requiring kernel workgroup redistribution.
* **NaN & Infinity Sanitization:** A single uninitialized array element or a division by zero in a force kernel will propagate `NaN` values across your whole trajectory. The automated harness must explicitly scan the raw output text/binary arrays for `NaN` or `Inf`.

### Level 2: Physical Range Intuition & Topology Mapping (Tier 1 & Tier 2)

Scientific code bugs often result in numbers that look clean (no NaNs) but are physically absurd.

* **Inferring Ranges:** If you do not provide explicit physical boundaries, the agent must inspect your units (e.g., eV vs. Hartree, Ångströms vs. Bohr) and establish a bounding box. For example: *"This is a localized atomic orbital calculation in eV. If any diagonal element of the Hamiltonian matrix exceeds $\pm 100\text{ eV}$, fail."*
* **Topology Alignment Checks:** When mapping data structures (especially moving from object-oriented C++ to flat arrays in OpenCL/Compute Shaders), index shifting is a common source of bugs. The agent should write a fast validation routine to check the neighborhood topology:
* Ensure padding values (like empty slots in a neighbor list) are explicitly set to `-1` and not left as raw memory remnants or zeros (which would falsely point to atom 0).
* Verify that if atom $A$ counts atom $B$ as a neighbor, the reverse is true or appropriately handled by your recoil strategy.

### Level 3: Qualitative Rigidity Scans & Human Debug Artifacts (Tier 2 & Tier 3)

When exact single-point energy parity cannot be achieved (e.g., extracting a core from DFTB+ or xTB into your own framework where basis sets or grid implementations vary slightly), **differential scanning** is your strongest tool.

Rather than looking at one static number, the agent executes a **Rigid Scan** (e.g., stretching a bond, rotating a dihedral, or sliding a molecule over an NaCl substrate) and saves the energy profile.

```
       Energy (eV)
         ^
         |      _--_   <-- Legacy Method (Reference)
         |     /    \
         |    /      \
         |   /........\ <-- New GPU Method (Qualitative Match)
         +------------------------> Coordinate (r)

```

* **LLM Analysis (Tier 2):** The LLM calculates the Root Mean Square Error (RMSE) of the curve slopes. It checks if the location of the global minimum ($r_e$) aligns within a tight physical tolerance (e.g., $\pm 0.05\text{ Å}$), even if the absolute baseline energy is shifted.
* **Human Artifact Generation (Tier 3):** If the curve matches qualitatively but needs final human verification, the loop must automatically plot a `.png` or generate an interactive HTML page using `p5.js` or `Plotly`. The agent must output a clear message:
> 📊 **Debug Artifact Generated:** I have rendered the rigid scan profile for the $H_2O\text{--}NaCl$ interaction at `outputs/plots/water_nacl_scan.html`. The absolute energy difference is $0.12\text{ eV}$ due to single-precision rounding, but the adsorption minimum matches perfectly at $2.4\text{ Å}$.

### Level 4: Geometric Relaxation & System Minimum Verification (Tier 1 & Tier 2)

A forcefield or quantum method might look correct during a single-point evaluation but fall apart when atoms begin to move.

* **Single-Point Directionality:** The harness checks if the force vectors point in the correct direction (e.g., forces should oppose the displacement vector during a bond stretch).
* **System Relaxation:** The agent runs a geometry optimization loop (e.g., Conjugate Gradient or BFGS) on **minimalist reference systems**. It checks if the final relaxed coordinates pass clear geometric rules:

| System | Classification | Expected Qualitative Geometric Outcome |
| --- | --- | --- |
| **$H_2$** | QM (2s Basis) | Converges to a stable bond length near $0.74\text{ Å}$. |
| **$N_2$** | QM (sp Basis) | Correctly handles $\sigma$ and $\pi$ bonding; triple bond distance logic. |
| **$H_2O$** | Classical / QM | Converges to a non-linear, bent structure ($\approx 104.5^\circ$). |
| **$NH_3$** | Classical / QM | Maintains a stable pyramidal structure. |
| **Benzene** | Sp2 / $\pi$-systems | Remains strictly planar. Planarity checks calculated via out-of-plane dihedral angles ($\theta \approx 0^\circ$). |
| **$H_2O$ on $NaCl$** | Surface Adsorption | **Oxygen** self-orients directly toward the $Na^+$ ion; **Hydrogens** pivot toward the $Cl^-$ ions at a distance matching the sum of their physical atomic radii. |

## 3. Level 5: The Performance Optimization Loop

Once your code passes the multi-level physical validity pipeline, the agentic loop switches its goal from **Accuracy** to **Throughput**.

The loop architecture changes focus here: it keeps the functional validation suite fixed as a regression guard while iterating on code variations to maximize performance metrics (GFLOPS, bandwidth utilization, or execution time).

```
   [Generate Kernel Variant] 
              |
              v
     [Run Accuracy Suite]  --- (Failed) ---> [Discard & Log Assembly Constraint]
              |
           (Passed)
              v
   [Measure Execution Time] ---> [Update Best Performance Record]

```

### Strategies for the Agent to Optimize GPU Compute:

* **Memory Layout Transformation:** The agent can rewrite structures from Array-of-Structures (AoS) to Structure-of-Arrays (SoA) to ensure coalesced global memory access, which is critical for maximizing performance on single-precision GPU architectures.
* **Workgroup Size Sweeping:** The loop can systematically benchmark the identical kernel across varying local workgroup sizes (e.g., 32, 64, 128, 256) to identify the execution sweet spot for a given hardware profile.
* **Avoiding Atomic Write Collisions (Recoil Forces):** Instead of using expensive global atomic additions (`atom_add`) when calculating pairwise forces ($F_{ij}$ and $F_{ji}$ updates hitting the same memory locations), the agent can re-architect the loop execution. It can configure the threads to compute interactions using an asymmetric strategy or assign specific spatial tiles to unique workgroups, avoiding costly GPU reduction phases entirely.

## 4. Conceptual Implementation: The Master Orchestrator

This unified framework shows how you can tie these levels together into a single automated script:

```python
class AgenticValidationHarness:
    def __init__(self, target_script, reference_data):
        self.script = target_script
        self.ref = reference_data

    def execute_tier1_checks(self):
        # 1. Compile and catch system crashes/device lockups
        status, logs = self.compile_kernel()
        if not status: return False, "CRASH_COMPILATION_ERROR", logs
        
        # 2. Run simulation and extract output matrices
        run_status, arrays = self.run_simulation_binary()
        if not run_status: return False, "CRASH_RUNTIME_ERROR", arrays
        
        # 3. Check for NaNs and basic array topology bounds
        if np.isnan(arrays['forces']).any():
            return False, "VALIDATION_NAN_DETECTED", "Forces contain NaN values."
        if (arrays['neighbor_list'] == 0).any(): # Assuming 0-index risk instead of -1
            return False, "VALIDATION_TOPOLOGY_ERROR", "Neighbor list contains unpadded 0 indices."
            
        return True, "TIER1_PASSED", arrays

    def execute_tier2_checks(self, arrays):
        # 1. Physical Conservation Assertions
        symmetry_error = np.max(np.abs(arrays['forces'] + arrays['forces'].T))
        if symmetry_error > 1e-4:
            return False, f"PHYSICS_VIOLATION: Newton's 3rd Law failed. Max error: {symmetry_error}"
            
        # 2. Run rigid scan calculation and check qualitative curve matching
        rmse_slope, tracking_report = evaluate_rigid_scan_profile(arrays, self.ref)
        if rmse_slope > 0.1: # Qualitative mismatch threshold
            return False, f"PHYSICS_VIOLATION: Rigid scan curve shape mismatch. Report: {tracking_report}"
            
        return True, "SUCCESS"

    def generate_tier3_artifacts(self, arrays):
        # Generate clean visual feedback so the human can verify optimization at a glance
        generate_p5js_coordinate_view(arrays['coordinates'])
        generate_plotly_curve(arrays['energy_profile'])
        print("Human review files updated at ./debug_dashboard.html")

```

By framing your goals around this multi-tiered architecture, you give the AI agent a clear, logical sandbox to safely explore alternative code paths. The system will independently resolve structural alignment errors and compiler issues, leaving you to step in only for high-level visual and physical validation.

---

# DeepSeek-V4-flahs

Excellent. You've pinpointed the hardest and most critical part of any agentic loop: **goal definition and verifiable success criteria**. Without that, the agent is just a fancy text generator spinning in circles.

Let me elaborate on how to structure these tests for your two primary use cases, respecting the economy of attention (cheap automatic → medium LLM → costly human). I'll also show how to embed them into the debug cycle.

## 1. Core Principle: Tests as Stateful Checkpoints

The agent loop must maintain a **test suite that grows with the project**. Each test has:
- **Type**: `compile`, `crash`, `numeric_sanity`, `parity`, `geometry_relax`, `performance`
- **Cost**: cheap (`C`), medium (`M`), expensive (`E`)
- **Expected outcome**: e.g., `max_force < 1e-6`, `RMSE < 1e-5 eV`, `speedup > 2.0x`
- **Reference data / tolerances**

The agent runs the cheapest tests first, and only if they pass does it proceed to more expensive ones.

## 2. Concrete Test Hierarchy for Your Use Cases

### Use Case A: Reimplementation & Parity Checking

You have a **reference implementation** (old Fortran, DFTB+, etc.) and a **new implementation** (C++/OpenCL/Rust). Goals:
1. **Identical results** (bit-exact or floating-point tolerant) – for exact reimplementation.
2. **Physically equivalent** (same minima, different energy offsets) – for reformulated methods.

#### Test Levels (in order of increasing cost)

| Level | Test | Automatic? | What it verifies | Example |
|-------|------|------------|------------------|---------|
| 0 | Compilation & linking | C – auto | No syntax errors | `make` succeeds |
| 1 | No crash on trivial input | C – auto | No segfault, device freeze | Run single H2O → exit code 0 |
| 2 | No NaNs / infinities in outputs | C – auto | Basic numeric stability | Check energy, forces for `isfinite()` |
| 3 | Physical range sanity | C – auto (with user/LLM supplied bounds) | Values within plausible physics | Bond force in H2: `-100..100 eV/Å`; Density: `0..100 e/Å³` |
| 4 | Neighbor list / topology check | C – auto | No off-by-one, uninitialized slots | Compare to reference: same number of neighbors; no index = -1 used |
| 5 | Single-point parity (identical method) | C – auto | Numerically identical to reference | `max|E_new - E_ref| < 1e-8 eV` ; `max|F_new - F_ref| < 1e-6 eV/Å` |
| 6 | Rigid scan (bond/angle/dihedral) | M – LLM generates plots, human reviews quickly | Same shape, location of minima | Scan H–H in H2; plot both curves; RMSE; visual check of plots |
| 7 | Geometry relaxation (energy min) | M – LLM checks convergence & final structure | Reasonable final geometry | H2O → bond angle ~104.5°, bond length ~0.96Å |
| 8 | Qualitative behavior (surfaces, adsorption) | H – human reviews .xyz/.cube | Chemistry makes sense | Water oxygen binds to Na⁺, hydrogen to Cl⁻ |
| 9 | Performance scaling (multiple sizes) | C – auto | Speed vs reference | `time_new / time_ref < 0.5` for N=100 atoms |

**Important**: For level 5, if methods are **different** (e.g., DFT vs forcefield), you skip exact parity and jump to levels 6–8 with relaxed tolerances (e.g., RMSE < 0.1 eV).

### Use Case B: Performance Optimization (OpenCL/GPU kernels)

You keep the same algorithm, but change memory layout, workgroup size, atomics vs reduction, etc.

| Level | Test | Cost | What it verifies |
|-------|------|------|------------------|
| 0–4 | Same as above (compile, no crash, no NaNs) | C | Correctness unchanged |
| P1 | Single-run timing on small system | C | Measurable speed (not regression) |
| P2 | Multiple runs + statistical timing (median, stddev) | C | Stable performance |
| P3 | Scaling test – weak / strong scaling | M (LLM plots) | Efficiency with problem size |
| P4 | Profiling – occupancy, memory bandwidth, warp divergence | M (LLM reads profiling output) | Identify bottlenecks |
| P5 | Compare to baseline (same hardware) – speedup | C | Improvement target (e.g., 2x) |
| P6 | Energy / power consumption (if available) | C | Efficiency for laptops |

The agent should **automatically revert** if P0–4 fail, and only propose new kernel variants if P5 is not met.

## 3. How the Agent Distinguishes Test Types & Infers Ranges

### For `physical range` (level 3):
- **User provides**: a simple YAML file per project, e.g.:
  ```yaml
  physical_limits:
    energy: [-1e6, 1e6] eV  # per atom
    force: [-100, 100] eV/Å
    charge_density: [0, 10] e/Å³
    bond_length_HH: [0.5, 2.0] Å
  ```
- **If missing, LLM infers**: Based on system description ("H2 molecule", "water box", "NaCl crystal") the LLM queries its internal knowledge or a reference table (e.g., typical bond lengths, atomic radii). It then **asks the user for confirmation** once, storing the bounds for future runs.

### For `qualitative geometry` (level 8):
- **LLM-generated test**: After relaxation, the agent writes a short Python script that:
  1. Reads final `.xyz`.
  2. Computes distances / angles (e.g., O–Na distance, H–Cl distance, H–O–H angle).
  3. Checks against simple rules: `O_Na_dist < 3.0 Å` and `H_Cl_dist < 3.0 Å` for water on NaCl.
  4. If rules pass → auto OK. If borderline → flag for human.
- **Human review aid**: The agent also generates a **comparison image** showing reference vs new geometry (using e.g., `matplotlib` or `ase.visualize`) and saves it as `debug/geom_check.png` with a caption describing what to look for. This is **cheap for human** because they just glance.

## 4. The Debug Cycle in Practice: Agentic Loop Pseudocode

```python
state = Goal("reimplement DFTB+ SCF loop in OpenCL")
state.add_test_suite("parity", tests=[compile, no_crash, no_nan, range_sanity, neighbor_topology])
state.add_test_suite("exact_parity", tests=[single_point_rmse < 1e-6])
state.add_test_suite("geometry", tests=[relax_to_minima, check_H2O_angles])

while not state.goal_achieved:
    action = agent.plan(state)   # e.g., edit kernel.cl
    agent.execute(action)        # compiles
    test_results = []
    for test in state.get_next_cheapest_test():
        result = test.run()
        test_results.append(result)
        if not result.passed:
            agent.log_failure(test, result)
            break   # no point running more expensive tests
    if all passed for current suite:
        state.advance_to_next_suite()
    agent.reflect(test_results, state)   # decides next edit
```

The **key** is `get_next_cheapest_test()` – it returns the first test in the hierarchy that hasn't passed yet. This prevents wasting time on expensive geometry relaxations if the code doesn't even compile.

## 5. Tools to Embed Into the Agent

To make this work, your agent needs a few helper tools that you can implement once:

- **`run_safe(command, timeout=30)`** – catches segfaults, device resets.
- **`check_nan_inf(output_file)`** – scans logs for `nan`, `inf`, `1.#QNAN`.
- **`compare_vectors(ref_file, new_file, tol)`** – RMSE, max deviation.
- **`plot_scan(x_vals, ref_y, new_y, xlabel, ylabel, output_path)`** – generates PNG and HTML (for p5.js).
- **`check_topology(neighbor_list, max_atoms)`** – looks for indices out of range, negative padding values used in computation.
- **`relax_and_geometry_check(input_structure, force_field, expected_rules)`** – wraps an optimizer and then applies rule‑based checks.

The agent should be able to **invoke these as functions** (like tools) – not generate the code each time.

## 6. Handling the "Different Method, Same Physics" Case

For your reformulated forcefields or new QM methods, **absolute parity is impossible**. Instead, define:

- **Conservation laws**: `dE/dt = power` – integrate forces, check energy drift < threshold.
- **Symmetries**: In NaCl crystal, forces on symmetry-equivalent atoms should be equal (up to numeric noise).
- **Common sense**: H₂O relaxation should not invert structure (oxygen between hydrogens).
- **Reference benchmark**: Use a small set of molecules (H2, N2, CO, H2O, formic acid dimer) and compute RMSE against high‑level DFT. The agent can **stop** if RMSE > 0.1 eV/atom, but only **warn** if RMSE is larger than expected for that method class.

You can store a **method‑specific tolerance matrix**:
```json
{
  "method": "my_new_dft",
  "reference": "cc-pVTZ",
  "rmse_energy_eV_per_atom": 0.05,
  "rmse_force_eV_per_A": 0.2
}
```

The agent checks against this, not against zero.

## 7. What the Agent Should Report to You (Human)

After each loop iteration, the agent should produce a **compact summary**:

```
Iteration 12:
  [OK] Compiles
  [OK] No crash
  [OK] No NaN
  [OK] Neighbors valid
  [FAIL] Single-point parity: max|F| = 0.12 eV/Å (tolerance 0.001)
    → Suspect: atomic forces on O mismatch. See debug/forces_O.png
  Next plan: Check neighbor list indexing for oxygen. Suggested edit: kernel.cl line 142.

[Human action required?] Not yet – will auto-retry up to 3 times.
```

Only when the agent gets stuck (e.g., same failure after 3 attempts) it should **pause and ask for human review**, attaching the relevant plots and `.xyz` files.

## 8. Concrete Example: Water on NaCl (Qualitative Test)

You want to ensure water adsorbs with oxygen toward Na⁺. The agent can:

1. Run a short MD or relaxation starting from water 3 Å above the surface.
2. Write a Python script that computes the **O–Na distance** and **H–Cl distances**.
3. Check: `min(O_Na) < 3.0 Å` **and** `min(H_Cl) < 3.0 Å`.
4. If true → **auto-pass**.
5. If false, generate a `.png` showing final geometry with atoms colored by element, and **ask human** "Is this physically reasonable? (y/n)". The human answers in one word, and the agent stores that answer for future runs.

This keeps human cost minimal while still catching weird failures.

## Final Advice

Start by **manually writing a small suite of level 0–5 tests** for one of your existing code pairs (e.g., Fortran vs C++). Then let an agentic loop (like `Loop-Harness` or `Whilly`) run overnight trying to fix the parity. The agent will rapidly learn which tests are most informative and where your codebases diverge.

The real win is **not having to re‑write the same test logic** for each new port – you just tell the agent "use the existing `water_on_nacl` qualitative test", and it knows what to measure and how to present it.

Would you like me to sketch the actual **YAML test definition format** and a **Python class for the agentic test runner** that you could plug into one of the frameworks mentioned earlier?

---

# Mistral 

Here’s a **hierarchical testing framework** for your agentic loops, optimized for **economy of attention** (automatic → LLM → human). This is tailored for reimplementation and performance optimization in computational chemistry/physics.

## 🎯 **Testing Pyramid: 8 Levels of Validation**

*(Cheapest at bottom, most expensive at top. Agent should **automate L1-L5**, use LLM for **L6-L7**, and only escalate to human for **L8** or failures)*

| Level | Type | What to Check | Automation | LLM Role | Human Role | Cost |
|-------|------|---------------|------------|----------|------------|------|
| **L1** | **Stability** | Compiles? No segfaults? No OpenCL freezes? | ✅ Full | None | None | **$** |
| **L2** | **Sanity** | No NaN/Inf? Values in physical ranges? | ✅ Full | Infer ranges | Define ranges | **$** |
| **L3** | **Topology** | Neighbor maps consistent? Proper padding? | ✅ Full | Generate test cases | Define rules | **$** |
| **L4** | **Parity** | RMSE vs reference? Same method, diff impl. | ✅ Full | Compare outputs | Define tolerance | **$** |
| **L5** | **Rigid Scans** | Energy/force profiles vs reference | ✅ Full | Generate plots, compute RMSE | Review plots | **$$** |
| **L6** | **Relaxation** | Converges to reasonable geometry? | ⚠️ Partial | Analyze trajectory | Visualize in VESTA | **$$$** |
| **L7** | **Qualitative** | H₂O on Na⁺? H toward Cl⁻? | ❌ No | Suggest test cases, analyze | Visual review | **$$$$** |
| **L8** | **Performance** | Kernel speed, memory layout, workgroup size | ❌ No | Suggest optimizations | Profile, decide | **$$$$** |


## 🔬 **Detailed Breakdown by Level**

### **🟢 L1: Stability Tests** *(Fully Automatic)*
**Goal**: Code doesn’t crash.

```python
# test_stability.py
import subprocess
import signal
import sys

def test_compilation():
    """L1: Does it compile?"""
    result = subprocess.run(["make"], capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, f"Compilation failed:\n{result.stderr}"

def test_no_crash(executable, input_file):
    """L1: Does it run without segfault/OpenCL freeze?"""
    def timeout_handler(signum, frame):
        raise TimeoutError("Process hung (OpenCL freeze?)")

    signal.signal(signal.SIGALRM, timeout_handler)
    signal.alarm(30)  # 30 second timeout

    try:
        result = subprocess.run([executable, input_file],
                              capture_output=True, text=True, timeout=25)
        assert result.returncode == 0, f"Crash detected:\n{result.stderr}"
    finally:
        signal.alarm(0)
```

**Agent Action**: Run before every commit. If fails → **stop loop**, fix immediately.

### **🟢 L2: Sanity Checks** *(Fully Automatic + LLM Range Inference)*
**Goal**: Outputs are physically plausible.

```python
# test_sanity.py
import numpy as np

# === User-defined physical ranges ===
PHYSICAL_RANGES = {
    "energy": {"H2": (-100, 0), "H2O": (-1000, 0), "NaCl": (-500, 50)},  # Hartree
    "force": {"H2": (-10, 10), "H2O": (-50, 50)},  # Hartree/Bohr
    "distance": {"bond": (0.5, 3.0), "vdw": (2.0, 5.0)},  # Å
    "angle": (0, 180),  # degrees
}

# === LLM can infer these if not provided ===
def infer_physical_range(molecule, property):
    """LLM generates reasonable ranges (call once, cache result)"""
    # This would be a prompt to LLM:
    # "What is a reasonable range for the bond energy of H2 in Hartree?"
    # Returns: (-4.5, -4.0)  # Actual H2 dissociation energy ~4.5 eV
    return PHYSICAL_RANGES.get(property, {}).get(molecule, (-1e6, 1e6))

def test_no_nan(output):
    """L2: Check for NaN/Inf"""
    assert not np.any(np.isnan(output)), "NaN detected in output"
    assert not np.any(np.isinf(output)), "Inf detected in output"

def test_physical_range(output, molecule, property):
    """L2: Check if values are in physical range"""
    low, high = infer_physical_range(molecule, property)
    assert np.all(output >= low), f"{property} below physical minimum ({low})"
    assert np.all(output <= high), f"{property} above physical maximum ({high})"
```

**Agent Action**:
1. Run automatically after L1 passes
2. If LLM needs to infer range → **cache result** for future runs
3. If out of range → **flag for L6/L7 review**

### **🟢 L3: Topology/Mapping Tests** *(Fully Automatic)*
**Goal**: Neighbor lists, atom indices, and array padding are consistent.

```python
# test_topology.py
import numpy as np

def test_neighbor_map_consistency(old_map, new_map, cutoff=5.0):
    """L3: Compare neighbor maps between implementations"""
    # Check same number of neighbors
    assert old_map.shape == new_map.shape, "Neighbor map shape mismatch"

    # Check same atoms are neighbors (order may differ)
    for i in range(len(old_map)):
        old_neighbors = set(old_map[i][old_map[i] != -1])
        new_neighbors = set(new_map[i][new_map[i] != -1])
        assert old_neighbors == new_neighbors, f"Neighbor mismatch for atom {i}"

def test_padding(old_array, new_array, padding_value=-1):
    """L3: Check arrays are properly padded"""
    # Find padding positions (where value == padding_value)
    old_pad = (old_array == padding_value)
    new_pad = (new_array == padding_value)

    # Same positions should be padded
    assert np.array_equal(old_pad, new_pad), "Padding pattern mismatch"

    # No uninitialized values (non-padding, non-valid)
    valid_range = (0, 1000000)  # Example: atom indices
    assert np.all((new_array[~new_pad] >= valid_range[0]) &
                  (new_array[~new_pad] <= valid_range[1])), "Uninitialized values detected"
```

**Agent Action**:
- Generate test cases with **random but valid** geometries
- Compare neighbor maps, angle/dihedral definitions between implementations
- **Critical for**: Forcefield reimplementations, QM basis set mappings

### **🟢 L4: Parity Tests** *(Fully Automatic)*
**Goal**: Same method, different implementation → same results (within tolerance).

```python
# test_parity.py
import numpy as np

def test_energy_parity(ref_energy, new_energy, rtol=1e-6, atol=1e-8):
    """L4: Energy parity check"""
    assert np.isclose(ref_energy, new_energy, rtol=rtol, atol=atol), \
        f"Energy mismatch: {ref_energy} vs {new_energy} (diff={abs(ref_energy-new_energy)})"

def test_force_parity(ref_forces, new_forces, rtol=1e-5):
    """L4: Force parity check (component-wise)"""
    assert np.allclose(ref_forces, new_forces, rtol=rtol), \
        f"Force mismatch: max diff = {np.max(np.abs(ref_forces - new_forces))}"

def test_full_parity(ref_output, new_output, properties):
    """L4: Check all properties in output"""
    for prop in properties:
        if prop in ["energy", "total_energy"]:
            test_energy_parity(ref_output[prop], new_output[prop])
        elif prop in ["forces", "gradients"]:
            test_force_parity(ref_output[prop], new_output[prop])
        elif prop in ["stress", "virial"]:
            test_stress_parity(ref_output[prop], new_output[prop])
```

**Agent Action**:
- Run on **identical inputs** for both implementations
- Use **tighter tolerances** for double precision, **looser** for single (GPU)
- If fails → **check L3 (topology) first**, then L2 (sanity)

### **🟡 L5: Rigid Scan Tests** *(Automatic + LLM Plot Generation)*
**Goal**: Energy/force profiles match when scanning degrees of freedom.

```python
# test_rigid_scans.py
import numpy as np
import matplotlib.pyplot as plt

def rigid_scan_bond_length(molecule, bond_index, scan_range=(0.5, 3.0), steps=20):
    """L5: Scan bond length, compare energy profiles"""
    distances = np.linspace(*scan_range, steps)
    ref_energies = []
    new_energies = []

    for d in distances:
        # Modify molecule geometry
        mol = molecule.copy()
        mol.set_bond_length(bond_index, d)

        # Compute with both implementations
        ref_energies.append(ref_implementation.energy(mol))
        new_energies.append(new_implementation.energy(mol))

    # Compute RMSE
    rmse = np.sqrt(np.mean((np.array(ref_energies) - np.array(new_energies))**2))

    # Generate plot for human review
    plt.figure(figsize=(10, 6))
    plt.plot(distances, ref_energies, 'b-', label='Reference')
    plt.plot(distances, new_energies, 'r--', label='New')
    plt.xlabel('Bond Length (Å)')
    plt.ylabel('Energy (Hartree)')
    plt.title(f'Bond Scan RMSE: {rmse:.2e}')
    plt.legend()
    plt.savefig(f'scan_bond_{bond_index}.png')
    plt.close()

    return rmse

def test_scan_parity(rmse, max_rmse=1e-4):
    """L5: Check if scan profile is close enough"""
    assert rmse < max_rmse, f"Scan RMSE {rmse} exceeds threshold {max_rmse}"
```

**Agent Action**:
1. **Automatically generate** scans for:
   - Bond lengths (H₂, O-H in H₂O)
   - Bond angles (H-O-H in H₂O)
   - Dihedrals (torsion in butane)
   - Rigid molecule over surface (H₂O over NaCl)
2. **Save plots** with clear labels in `debug/plots/`
3. **LLM analyzes** the RMSE and suggests if difference is acceptable
4. **Human reviews** plots only if RMSE > threshold

**Example Test Cases**:
```python
# For forcefield validation
test_cases = [
    ("H2", "bond_length", (0.5, 1.5)),  # H-H bond
    ("H2O", "bond_angle", (90, 120)),     # H-O-H angle
    ("H2O", "scan_over_NaCl", (2.0, 5.0)), # Distance from surface
    ("benzene", "dihedral", (0, 180)),   # Ring puckering
]
```

### **🟡 L6: Relaxation Tests** *(Partial Automation + LLM Analysis)*
**Goal**: Does the new implementation relax to a **physically reasonable** geometry?

```python
# test_relaxation.py
import numpy as np
from ase.io import write  # For .xyz output

def test_relaxation(molecule, max_steps=100, tol=1e-4):
    """L6: Relax structure and check convergence"""
    # Run relaxation with new implementation
    traj = new_implementation.relax(molecule, max_steps=max_steps, tol=tol)

    # Save trajectory for human review
    write(f"debug/trajectories/{molecule.name}_relax.xyz", traj)

    # Check convergence
    energies = [frame.get_potential_energy() for frame in traj]
    forces = [np.linalg.norm(frame.get_forces()) for frame in traj]

    assert len(energies) < max_steps, "Did not converge in max steps"
    assert forces[-1] < tol, f"Final force norm {forces[-1]} > tolerance {tol}"

    # LLM can analyze trajectory
    analysis = analyze_trajectory(traj, molecule)
    return analysis

def analyze_trajectory(traj, molecule):
    """LLM analyzes if relaxation is physically reasonable"""
    # This would be a prompt to LLM:
    prompt = f"""
    Analyze this relaxation trajectory for {molecule.name}:
    - Initial geometry: {traj[0].get_positions()}
    - Final geometry: {traj[-1].get_positions()}
    - Energy change: {traj[-1].get_potential_energy() - traj[0].get_potential_energy()}
    - Final forces: {np.linalg.norm(traj[-1].get_forces())}

    Is this relaxation physically reasonable? Consider:
    1. Are bond lengths in expected ranges?
    2. Are angles in expected ranges?
    3. Did the molecule dissociate unexpectedly?
    4. For {molecule.name}, what is the expected equilibrium geometry?
    """

    # LLM response would be cached and used for pass/fail
    return "PASS"  # or "FAIL: Bond length too short"
```

**Agent Action**:
- Save `.xyz` files in `debug/trajectories/`
- LLM **generates a summary** of the relaxation
- Human **only reviews** if LLM flags as suspicious

### **🟠 L7: Qualitative Physical Tests** *(LLM + Human Visual Review)*
**Goal**: Does the implementation capture **essential physics**?

| Test Case | What to Check | LLM Role | Human Role |
|-----------|---------------|----------|------------|
| **H₂ (2s)** | Binding energy ~4.5 eV, bond length ~0.74 Å | Suggest test, analyze output | Visualize MO diagram |
| **N₂ (sp)** | Triple bond, bond length ~1.1 Å | Suggest test, analyze output | Review bond order |
| **CO** | Polar bond (C⁻-O⁺), dipole moment | Suggest test, analyze output | Visualize electron density |
| **H₂O** | Tetrahedral, O-H ~0.96 Å, H-O-H ~104.5° | Suggest test, analyze output | Review in VESTA |
| **NaCl** | Cubic, Na⁺-Cl⁻ distance ~2.8 Å | Suggest test, analyze output | Visualize crystal |
| **H₂O + NaCl** | O on Na⁺, H on Cl⁻ | **Generate test case**, analyze | Visualize in VESTA |
| **Carboxylic acid dimer** | H-bonding pattern | Suggest test, analyze output | Visualize dimer |
| **Benzene** | Planar, C-C ~1.4 Å, D6h symmetry | Suggest test, analyze output | Review structure |

```python
# test_qualitative.py
def test_h2o_nacl_interaction():
    """L7: Check if H2O adsorbs correctly on NaCl"""
    # Set up system: NaCl surface + H2O
    nacl = NaClSurface((3, 3, 1))
    h2o = H2OMolecule()
    
    # Place H2O at various positions
    test_positions = [
        nacl.get_na_position() + [0, 0, 2.0],  # Above Na+
        nacl.get_cl_position() + [0, 0, 2.0],  # Above Cl-
    ]
    
    results = []
    for pos in test_positions:
        h2o.set_position(pos)
        energy = new_implementation.energy(nacl + h2o)
        forces = new_implementation.forces(h2o)
        results.append({
            "position": pos,
            "energy": energy,
            "forces_on_O": forces[0],  # Oxygen
            "forces_on_H": forces[1:], # Hydrogens
        })
    
    # LLM analyzes results
    analysis = llm_analyze_h2o_nacl(results)
    
    # Save for human review
    with open("debug/qualitative/h2o_nacl_analysis.md", "w") as f:
        f.write(analysis)
    
    # Human reviews in VESTA
    write("debug/structures/h2o_above_na.xyz", nacl + h2o)
    
    return analysis

def llm_analyze_h2o_nacl(results):
    """LLM checks if physics is correct"""
    prompt = f"""
    Analyze H2O adsorption on NaCl:

    Above Na+ (position {results[0]['position']}):
    - Energy: {results[0]['energy']}
    - Force on O: {results[0]['forces_on_O']}
    - Forces on H: {results[0]['forces_on_H']}

    Above Cl- (position {results[1]['position']}):
    - Energy: {results[1]['energy']}
    - Force on O: {results[1]['forces_on_O']}
    - Forces on H: {results[1]['forces_on_H']}

    Expected behavior:
    1. Above Na+: Oxygen should be attracted (negative force z-component)
       Hydrogens should be repelled (positive force z-component)
    2. Above Cl-: Hydrogens should be attracted (negative force z-component)
       Oxygen should be repelled (positive force z-component)
    3. Energy should be lower above Na+ than above Cl-

    Is this physically correct? If not, what's wrong?
    """
    # LLM response
    return "PASS: Correct adsorption behavior"
```

**Agent Action**:
- **LLM generates** the test cases (e.g., "test H₂O on NaCl")
- **LLM analyzes** the numerical output
- **Human reviews** only the **visualizations** (VESTA, plots)
- **Save all artifacts** in `debug/qualitative/`

### **🔴 L8: Performance Optimization** *(LLM Suggestions + Human Decision)*
**Goal**: Optimize GPU kernels, memory layout, workgroup sizes.

```python
# test_performance.py
import time
import numpy as np

def benchmark_kernel(kernel_func, input_data, warmup=10, runs=100):
    """L8: Benchmark a kernel"""
    # Warmup
    for _ in range(warmup):
        kernel_func(input_data)

    # Measure
    times = []
    for _ in range(runs):
        start = time.perf_counter()
        kernel_func(input_data)
        times.append(time.perf_counter() - start)

    avg_time = np.mean(times)
    std_time = np.std(times)
    return {"avg": avg_time, "std": std_time, "min": min(times), "max": max(times)}

def test_kernel_optimizations():
    """L8: Compare different kernel versions"""
    # Test different workgroup sizes
    workgroup_sizes = [(8, 8, 8), (16, 16, 4), (32, 8, 2), (64, 4, 1)]

    results = {}
    for wg in workgroup_sizes:
        kernel = create_kernel(workgroup_size=wg)
        results[str(wg)] = benchmark_kernel(kernel, test_data)

    # LLM suggests best option
    best = llm_suggest_best_kernel(results)
    return best

def llm_suggest_best_kernel(results):
    """LLM analyzes performance data and suggests optimizations"""
    prompt = f"""
    Kernel performance results (seconds):
    {results}

    Suggest:
    1. Which workgroup size is fastest?
    2. Is the performance scaling as expected with workgroup size?
    3. Should we try other sizes? Which ones?
    4. Are there likely memory access patterns causing slowdowns?
    """
    # LLM response
    return "(16, 16, 4)"  # Example
```

**Agent Action**:
- **LLM suggests** optimizations (workgroup size, memory layout)
- **Human profiles** with `nvprof`, `nsight`, or AMD tools
- **Human decides** final implementation

## 🤖 **Agentic Loop Integration**

### **Reimplementation Workflow**
```python
# agent_reimplement.py
from agentic_loop import Agent

def reimplementation_loop(source_code, target_language, test_cases):
    agent = Agent(
        goal=f"Reimplement {source_code} in {target_language} with parity to reference",
        max_iterations=50,
    )

    while not agent.goal_achieved():
        # === L1: Stability ===
        if not agent.test_compilation():
            agent.fix_compilation_errors()
            continue

        # === L2: Sanity ===
        if not agent.test_sanity():
            agent.fix_nan_inf()
            continue

        # === L3: Topology ===
        if not agent.test_topology():
            agent.fix_neighbor_maps()
            continue

        # === L4: Parity ===
        if not agent.test_parity(test_cases):
            agent.debug_parity_issues()
            continue

        # === L5: Rigid Scans ===
        if not agent.test_rigid_scans():
            agent.analyze_scan_differences()
            continue

        # === L6: Relaxation ===
        if not agent.test_relaxation():
            agent.analyze_relaxation()
            continue

        # === L7: Qualitative ===
        if not agent.test_qualitative():
            # Generate plots for human
            agent.generate_qualitative_artifacts()
            agent.request_human_review("Check debug/qualitative/")
            continue

        # === L8: Performance ===
        agent.optimize_performance()

    return agent.get_final_code()
```

### **Performance Optimization Workflow**
```python
# agent_optimize.py
def optimization_loop(kernel, test_data):
    agent = Agent(
        goal=f"Optimize {kernel} for speed on GPU",
        max_iterations=20,
    )

    while not agent.goal_achieved():
        # === L1-L4: Must still pass ===
        if not agent.test_stability() or not agent.test_parity():
            agent.revert_changes()
            continue

        # === L8: Performance ===
        perf = agent.benchmark_kernel()
        if perf["avg"] < agent.best_performance:
            agent.best_performance = perf["avg"]
            agent.save_best_kernel()
        else:
            agent.revert_changes()

        # LLM suggests next optimization
        next_optimization = agent.llm_suggest_optimization(perf)
        agent.apply_optimization(next_optimization)

    return agent.get_best_kernel()
```

## 📁 **Directory Structure for Debug Artifacts**

```
your_repo/
├── src/
│   ├── fortran_original/
│   ├── cpp_reimplementation/
│   └── opencl_gpu/
│
├── tests/
│   ├── L1_stability/
│   ├── L2_sanity/
│   ├── L3_topology/
│   ├── L4_parity/
│   ├── L5_rigid_scans/
│   ├── L6_relaxation/
│   └── L7_qualitative/
│
└── debug/  # Generated by agent for human review
    ├── plots/          # L5: Energy profiles, force curves
    │   ├── scan_bond_0.png
    │   └── scan_angle_HOH.png
    │
    ├── trajectories/   # L6: Relaxation paths (.xyz)
    │   ├── H2O_relax.xyz
    │   └── benzene_relax.xyz
    │
    ├── structures/     # L7: For VESTA review
    │   ├── H2O_above_Na.xyz
    │   └── carboxylic_dimer.xyz
    │
    ├── qualitative/    # L7: LLM analysis
    │   ├── h2o_nacl_analysis.md
    │   └── benzene_symmetry_check.md
    │
    └── performance/    # L8: Benchmark data
        ├── kernel_benchmark.json
        └── nvprof_output.txt
```

## 🎯 **Economy of Attention: When to Escalate**

| Situation | Automatic | LLM | Human | Action |
|-----------|-----------|-----|-------|--------|
| Compilation error | ✅ | ❌ | ❌ | Fix immediately |
| NaN/Inf detected | ✅ | ❌ | ❌ | Check input, fix |
| Topology mismatch | ✅ | ❌ | ❌ | Debug neighbor maps |
| Parity RMSE > threshold | ❌ | ✅ | ❌ | LLM analyzes differences |
| Scan profile weird | ❌ | ✅ | ⚠️ | LLM flags, human reviews plot |
| Relaxation suspicious | ❌ | ✅ | ✅ | LLM + human review trajectory |
| Qualitative test fails | ❌ | ✅ | ✅ | LLM + human review visualization |
| Performance plateau | ❌ | ✅ | ✅ | LLM suggests, human profiles |

## 🚀 **Recommended Test Cases by Domain**

### **Quantum Mechanics (QM)**
| Molecule | Basis | Test | Expected |
|----------|-------|------|----------|
| H₂ | STO-3G | Binding energy | ~3.5 eV |
| H₂ | 6-31G* | Binding energy | ~4.5 eV |
| He | STO-3G | Energy | ~-2.86 Hartree |
| H₂O | STO-3G | Geometry | O-H ~0.96 Å, ∠HOH ~104.5° |
| N₂ | STO-3G | Bond length | ~1.1 Å |
| CO | STO-3G | Dipole moment | ~0.1 D (small, C⁻-O⁺) |

### **Classical Forcefields**
| System | Test | Expected |
|--------|------|----------|
| H₂O | Single point | O-H ~0.96 Å, ∠HOH ~104.5° |
| NH₃ | Single point | N-H ~1.01 Å, ∠HNH ~107° |
| Formaldehyde (H₂C=O) | Single point | C=O ~1.21 Å, C-H ~1.10 Å |
| Formic acid (HCOOH) | Single point | C=O ~1.20 Å, O-H ~0.96 Å |
| Benzene | Single point | Planar, C-C ~1.40 Å |
| Pyridine | Single point | Planar, C-N ~1.34 Å |
| Carboxylic acid dimer | Single point | H-bond length ~1.7-2.0 Å |
| NaCl (bulk) | Single point | Lattice constant ~5.64 Å |
| H₂O + NaCl | Adsorption | O on Na⁺, H on Cl⁻ |

### **Non-Covalent Interactions**
| System | Test | Expected |
|--------|------|----------|
| H₂O dimer | Binding energy | ~5 kcal/mol |
| Benzene dimer | Binding energy | ~2 kcal/mol (π-π stacking) |
| Na⁺ + Cl⁻ | Binding energy | ~-7.9 kcal/mol (gas phase) |
| H₂O + Na⁺ | Binding energy | ~-25 kcal/mol |
| H₂O + Cl⁻ | Binding energy | ~-15 kcal/mol |

## 💡 **Key Principles for Your Agentic Loop**

1. **Fail Fast**: Run L1-L4 **automatically on every iteration**. Only proceed if all pass.
2. **Cache LLM Calls**: Physical ranges, test case suggestions → **cache in `tests/config.yaml`**
3. **Generate Artifacts**: Plots, trajectories, structures → **save in `debug/` for human review**
4. **Escalate Smartly**: Only ask human when **LLM can’t decide** or **visual review is needed**
5. **Document Decisions**: Every human review → **add note in `debug/REVIEW.md`**
6. **Parallelize Tests**: Run L1-L4 in parallel (they’re independent)
7. **Prioritize Tests**: Run **fastest tests first** (L1-L2 before L5-L7)
## **Next Steps for You**

1. **Pick one reimplementation project** (e.g., Fortran → OpenCL for a specific forcefield)
2. **Set up the 8-level test hierarchy** for that project
3. **Start with L1-L4 automation** (these are easiest)
4. **Add L5-L7 gradually** as you identify edge cases
5. **Use an agent framework** (OpenClaw, Hermes, or custom) to run the loop

Would you like me to:
- **Generate a starter `test_*.py` template** for one of your specific codes?
- **Design the YAML config** for your test hierarchy?
- **Write the agent loop controller** in Python?
- **Create a specific test case** (e.g., H₂O + NaCl) with all 8 levels?

---

# Kimi 2.6 chat

I'll search for more specific examples of agentic testing frameworks and scientific code migration patterns, then elaborate on your excellent framework.
Your framework is excellent and maps directly to what cutting-edge research is doing. Let me elaborate each layer with concrete implementations, drawing from recent work on **agentic translation loops** (Ralph loops for scientific code migration) and **property-based testing** for physics invariants. 

## The Core Principle: Goal-Conditioned Agentic Loops

The breakthrough insight from the 2026 land-surface model translation paper is the **Ralph loop** — an orchestration pattern that intercepts the agent's exit signal and re-injects the driving prompt until a verifiable goal condition is met. This addresses "agentic laziness" (agents declaring success prematurely) and statelessness across sessions. 

```bash
# The Ralph loop for scientific code migration
while true; do
    cat PROMPT_translate.md | claude --dangerously-skip-permissions
    # Loop forces re-evaluation of all plan.md tasks before exit
    # plan.md and CHANGELOG.md provide external memory
done
```

For your workflow, this becomes:

```python
class ScientificAgentLoop:
    def __init__(self):
        self.plan = "plan.md"          # What needs to be done
        self.changelog = "CHANGELOG.md" # What was attempted, what failed
        self.memory = "MEMORY.md"       # Physics invariants learned
        self.budget = IterationBudget(max=90)  # From Hermes
    
    def run(self, goal: str):
        while not self.goal_achieved(goal) and self.budget.remaining > 0:
            # 1. ASSESS: Read current state
            context = self.read_plan() + self.read_changelog() + self.read_memory()
            
            # 2. ACT: Generate code/tests/analysis
            action = llm.generate(context, goal)
            
            # 3. EXECUTE: Run the action (compile, test, etc.)
            result = self.execute(action)
            
            # 4. VERIFY: Check goal condition
            status = self.verify(result, goal)
            
            # 5. UPDATE: Write to external memory
            self.update_changelog(action, result, status)
            
            # 6. DECIDE: Continue, repair, or escalate
            if status == "PASS":
                self.advance_plan()
            elif status == "FAIL":
                self.spawn_repair_agent(action, result)
            elif status == "NEEDS_HUMAN":
                self.flag_for_human_review()
```

## Level 1: Compile & Runtime Sanity (100% Automatic, ~ms)

**Goal**: Code compiles, runs without crashing, produces finite numbers in physically reasonable ranges.

```python
class Level1RuntimeTest:
    """Cheap, automatic, must pass on every iteration"""
    
    def test_compile(self, source_files, build_cmd):
        """Automatic: Compiler does the work"""
        result = subprocess.run(build_cmd, capture_output=True)
        assert result.returncode == 0, f"Compile failed: {result.stderr}"
        return True
    
    def test_no_nan_inf(self, executable, test_input):
        """Automatic: Runtime check"""
        output = run(executable, test_input)
        assert not np.isnan(output).any(), "NaN detected"
        assert not np.isinf(output).any(), "Inf detected"
        return True
    
    def test_reasonable_range(self, output, property_name):
        """
        Can be automatic if user provides ranges, 
        or LLM-inferred from physical intuition
        """
        # User-provided: "Lennard-Jones epsilon should be 0.1-10 kJ/mol"
        # LLM-inferred: "Water O-H bond length should be ~0.96 Å"
        
        ranges = self.get_ranges(property_name)  # From MEMORY.md
        assert ranges['min'] <= output <= ranges['max'], \
            f"{property_name} = {output} outside reasonable range [{ranges['min']}, {ranges['max']}]"
        return True
```

**Physical Range Inference** (LLM task, done once, cached):
```python
# Agent queries LLM: "What is reasonable range for TIP4P/2005 O-H bond length?"
# LLM responds: "0.95-0.96 Å based on known water models"
# Stored in MEMORY.md for future reuse
```

## Level 2: Topology & Mapping (95% Automatic, ~ms-s)

**Goal**: Atoms, bonds, angles, neighbor lists point to the right things. No uninitialized arrays. Proper padding.

This is where most reimplementation bugs hide. The agent must verify **structural isomorphism** between reference and target.

```python
class Level2TopologyTest:
    def test_neighbor_list_consistency(self, ref_nbl, new_nbl):
        """
        Automatic: Compare topology maps
        Critical for Fortran→C++ where indexing changes (1-based→0-based)
        """
        assert len(ref_nbl) == len(new_nbl), "Atom count mismatch"
        
        for i in range(len(ref_nbl)):
            ref_neighbors = set(ref_nbl[i])
            new_neighbors = set(new_nbl[i])
            assert ref_neighbors == new_neighbors, \
                f"Atom {i}: neighbor mismatch {ref_neighbors} vs {new_neighbors}"
        
        # Check padding consistency (e.g., -1 for empty slots)
        assert (new_nbl[new_nbl == -1]).all() == (ref_nbl[ref_nbl == -1]).all()
        return True
    
    def test_angle_indices(self, ref_angles, new_angles, atom_map):
        """
        Angles must point to same atoms after reindexing
        """
        for ref_angle, new_angle in zip(ref_angles, new_angles):
            # atom_map: old index → new index (for reordered systems)
            mapped = [atom_map[i] for i in ref_angle]
            assert set(mapped) == set(new_angle), \
                f"Angle atom mismatch: {ref_angle} → {mapped} vs {new_angle}"
        return True
    
    def test_no_uninitialized(self, code_ast):
        """
        Static analysis: Check all arrays are initialized before use
        Critical for OpenCL where uninitialized local memory causes silent errors
        """
        # Use tree-sitter or clang AST to find declaration-use chains
        uninitialized = static_analysis.find_uninitialized_arrays(code_ast)
        assert len(uninitialized) == 0, f"Uninitialized arrays: {uninitialized}"
        return True
```

**Key insight**: Topology errors are **silent killers** — the code runs but gives wrong physics. The agent must generate **mapping verification tests** automatically by comparing the graph structure of the reference and reimplemented systems.

## Level 3: Parity & RMSE (80% Automatic, ~s)

**Goal**: Rigid scans, single-point calculations match reference within tolerance. But "better" methods are allowed to differ systematically.

This is the core of your reimplementation workflow. The Ralph loop paper uses **golden I/O files** — recorded inputs/outputs from the reference code that serve as oracles. 

```python
class Level3ParityTest:
    def __init__(self, reference_executable, tolerance_strategy):
        self.ref = reference_executable
        self.tol = tolerance_strategy  # "exact", "relative", "systematic_shift"
    
    def test_rigid_scan_parity(self, scan_type, coordinates_range):
        """
        Generate scan (e.g., O-H bond length from 0.8 to 1.2 Å)
        Compare energy profiles
        """
        # Automatic: Generate scan points
        scan_points = np.linspace(coordinates_range['min'], 
                                  coordinates_range['max'], 50)
        
        ref_energies = []
        new_energies = []
        
        for r in scan_points:
            coords = self.generate_geometry(scan_type, r)
            ref_energies.append(self.run_reference(coords))
            new_energies.append(self.run_target(coords))
        
        # Automatic: Compute RMSE
        rmse = np.sqrt(np.mean((np.array(ref_energies) - np.array(new_energies))**2))
        
        # LLM Analysis: Determine if difference is acceptable
        if self.tol == "exact":
            assert rmse < 1e-6, f"RMSE {rmse} exceeds exact tolerance"
        elif self.tol == "relative":
            relative_rmse = rmse / np.mean(np.abs(ref_energies))
            assert relative_rmse < 1e-4, f"Relative RMSE {relative_rmse} exceeds 1e-4"
        elif self.tol == "systematic_shift":
            # "Better" method may have constant offset but same shape
            correlation = np.corrcoef(ref_energies, new_energies)[0,1]
            assert correlation > 0.999, f"Correlation {correlation} too low"
            # LLM checks if offset is physically reasonable (e.g., basis set difference)
        
        return {
            'rmse': rmse,
            'correlation': correlation if self.tol == 'systematic_shift' else None,
            'plot_path': self.generate_plot(scan_points, ref_energies, new_energies)
        }
    
    def generate_plot(self, x, ref, new):
        """LLM generates plot, human reviews quickly"""
        fig, ax = plt.subplots()
        ax.plot(x, ref, 'b-', label='Reference', linewidth=2)
        ax.plot(x, new, 'r--', label='Reimplementation', linewidth=2)
        ax.set_xlabel('Scan coordinate (Å)')
        ax.set_ylabel('Energy (kJ/mol)')
        ax.legend()
        path = f"/debug_artifacts/scan_{uuid4()}.png"
        fig.savefig(path)
        return path
```

**Debug Artifacts**: The agent must produce **well-documented plots** with clear filenames (`scan_OH_bond_LJ_vs_Morse_2024-06-14.png`) and a summary report pointing humans to exactly where to look.

## Level 4: Physical Invariants (60% Automatic, ~s-min)

**Goal**: Conservation laws, symmetries, detailed balance hold. These are **physics contracts** that any valid implementation must satisfy, regardless of method.

```python
class Level4InvariantTest:
    def test_energy_conservation_nve(self, trajectory, tolerance=1e-3):
        """
        Automatic: NVE ensemble must conserve total energy
        """
        E_total = trajectory.potential + trajectory.kinetic
        drift = np.max(E_total) - np.min(E_total)
        relative_drift = drift / np.mean(np.abs(E_total))
        
        assert relative_drift < tolerance, \
            f"Energy drift {relative_drift} exceeds {tolerance}"
        return True
    
    def test_symmetry(self, system, operation):
        """
        Automatic: Apply symmetry operation, check invariance
        """
        E_before = self.compute_energy(system)
        system_transformed = operation.apply(system)  # rotation, reflection
        E_after = self.compute_energy(system_transformed)
        
        assert abs(E_before - E_after) < 1e-10, \
            f"Symmetry broken: {E_before} vs {E_after}"
        return True
    
    def test_detailed_balance(self, transition_matrix):
        """
        Automatic: For Monte Carlo / kinetic methods
        """
        # π_i P(i→j) = π_j P(j→i)
        for i in range(n_states):
            for j in range(i+1, n_states):
                lhs = equilibrium_prob[i] * transition_matrix[i,j]
                rhs = equilibrium_prob[j] * transition_matrix[j,i]
                assert abs(lhs - rhs) < tolerance, \
                    f"Detailed balance violated for {i}→{j}"
        return True
    
    def test_virial_theorem(self, trajectory):
        """
        Automatic: 2<T> = -<<Σ r·F> for Coulomb systems
        """
        # ... implementation
        pass
```

**Property-Based Testing** insight: These invariants can be tested with **randomized inputs** (Hypothesis/QuickCheck style) rather than hand-crafted cases. The agent generates random configurations and checks that invariants hold across the space. 

## Level 5: Qualitative Geometry (40% Automatic, ~min)

**Goal**: Does H2O bind to Na+ with O-down? Does benzene stay planar? These require **physical intuition** — LLM territory.

```python
class Level5QualitativeTest:
    def test_water_on_nacl(self, optimized_geometry):
        """
        LLM Analysis: Check qualitative binding geometry
        """
        # Automatic: Extract distances and angles
        d_O_Na = distance(optimized_geometry.O, optimized_geometry.Na)
        d_H_Na = distance(optimized_geometry.H1, optimized_geometry.Na)
        d_H_Cl = distance(optimized_geometry.H1, optimized_geometry.Cl)
        
        # LLM evaluates: "Is O closer to Na+ than H is?"
        # This is a qualitative check, not exact distance
        qualitative_checks = {
            'O_closer_to_Na': d_O_Na < d_H_Na,
            'H_toward_Cl': d_H_Cl < d_O_Na,
            'reasonable_distance': 2.0 < d_O_Na < 3.5  # Å, from ionic radii
        }
        
        # LLM generates natural language report
        report = f"""
        Water on NaCl Surface - Qualitative Analysis:
        - O-Na distance: {d_O_Na:.2f} Å (expected: 2.0-3.5 Å)
        - H-Na distance: {d_H_Na:.2f} Å (should be > O-Na)
        - H-Cl distance: {d_H_Cl:.2f} Å (should be < O-Na for H-bonding to Cl-)
        
        Verdict: {'PASS' if all(qualitative_checks.values()) else 'FAIL'}
        """
        
        # Human reviews this report + geometry visualization
        return report, self.visualize_geometry(optimized_geometry)
```

**Test System Hierarchy** (from your notes):
| System | What It Tests | Method |
|--------|--------------|--------|
| H2 (2 s orbitals) | Minimal QM basis | Single-point |
| N2 (sp basis) | Bond formation | Optimization |
| CO (polarization) | Different atoms | Dipole moment |
| H2O | 3-body angles, lone pairs | Geometry |
| NH3 | Tetrahedral (sp3) | Inversion barrier |
| Formaldehyde | sp2, π-bond | Planarity |
| Benzene/Pyridine | π-system, aromaticity | Delocalization |
| Carboxylic dimer | H-bonding | Dimer energy |
| NaCl surface | Ionic, non-covalent | Adsorption geometry |

The agent should **progressively test** — start with H2, only advance to benzene if H2 passes.

## Level 6: DFT Reference (10% Automatic, ~hours)

**Goal**: Ultimate ground truth. Expensive, run only when cheaper tests pass.

```python
class Level6DFTReference:
    def __init__(self):
        self.dft_cache = {}  # Cache results to avoid recomputation
    
    def test_against_dft(self, system, method_under_test):
        """
        Run only if Levels 1-5 pass
        """
        cache_key = hash(system) + hash(method_under_test)
        if cache_key in self.dft_cache:
            return self.dft_cache[cache_key]
        
        # Expensive: Run DFT calculation
        dft_energy = self.run_dft(system)  # Hours on HPC
        
        # Compare with method under test
        test_energy = method_under_test.compute(system)
        
        # LLM evaluates: Is difference acceptable for this level of theory?
        error = abs(dft_energy - test_energy)
        
        # Store in cache
        self.dft_cache[cache_key] = error
        return error
```

**Caching is critical**: The agent must never re-run DFT for the same system. The `MEMORY.md` tracks what has been validated.

## The Multi-Agent Architecture

For your workflow, I recommend a **specialized multi-agent system** similar to the code review architecture, but adapted for scientific computing: 

```python
class ScientificAgentSwarm:
    def __init__(self):
        self.agents = {
            'translator': TranslationAgent(),      # Fortran → C++/OpenCL
            'test_generator': TestGeneratorAgent(), # Creates tests from physics
            'parity_checker': ParityAgent(),       # Compares outputs
            'invariant_verifier': InvariantAgent(), # Checks physics laws
            'performance_optimizer': PerformanceAgent(), # GPU kernel tuning
            'human_liaison': HumanReviewAgent(),   # Prepares visual artifacts
        }
        self.consensus = ConsensusAgent()
    
    def migrate_module(self, source_code, target_platform, goal):
        """
        Main workflow: Translate, test, repair, optimize
        """
        # Phase 1: Translation (subagent with budget)
        translation = self.agents['translator'].run(
            source=source_code, 
            target=target_platform,
            budget=30
        )
        
        # Phase 2: Generate tests (from physics invariants)
        tests = self.agents['test_generator'].generate(
            code=translation,
            physics_domain='forcefield',  # or 'QM', 'MD', etc.
            levels=[1, 2, 3, 4]  # Start cheap
        )
        
        # Phase 3: Run tests in order
        for level in [1, 2, 3, 4, 5]:
            results = self.run_level(level, tests, translation)
            if results['status'] == 'FAIL':
                # Spawn repair agent
                repair = self.agents['translator'].repair(
                    code=translation,
                    test_results=results,
                    budget=20
                )
                # Re-run from Level 1
                continue
        
        # Phase 4: Performance optimization (if tests pass)
        if target_platform in ['OpenCL', 'CUDA', 'WebGPU']:
            optimized = self.agents['performance_optimizer'].tune(
                code=translation,
                metrics=['memory_coalescing', 'workgroup_size', 'atomics_avoidance']
            )
            # Re-verify parity after optimization
        
        # Phase 5: Human review preparation
        artifacts = self.agents['human_liaison'].prepare(
            plots=self.collect_plots(),
            trajectories=self.collect_trajectories(),
            summary=self.generate_summary()
        )
        
        return optimized, artifacts
```

## Performance Optimization Agent (Your Specific Need)

For GPU kernel tuning, the agent needs **platform-specific skills**:

```python
class PerformanceAgent:
    def __init__(self):
        self.skills = {
            'opencl': 'SKILL_opencl_optimization.md',
            'cuda': 'SKILL_cuda_optimization.md',
            'webgpu': 'SKILL_webgpu_optimization.md',
        }
    
    def tune(self, kernel_code, platform):
        """
        Strategies from your notes:
        1. Memory layout (SOA vs AOS)
        2. Workgroup size (typically 64 for GPUs)
        3. Global memory coalescing
        4. Avoid atomics (use reduction strategies)
        5. Single precision vs double precision tradeoffs
        """
        strategies = self.load_skill(self.skills[platform])
        
        # Agentic loop: Try strategy → benchmark → compare → keep if better
        best_kernel = kernel_code
        best_time = float('inf')
        
        for strategy in strategies:
            candidate = self.apply_strategy(kernel_code, strategy)
            compile_time = self.benchmark_compile(candidate)
            if compile_time == 'FAIL':
                continue
            
            run_time = self.benchmark_run(candidate)
            parity = self.check_parity(candidate, best_kernel)  # Must still be correct!
            
            if run_time < best_time and parity['pass']:
                best_kernel = candidate
                best_time = run_time
                self.log(f"Strategy {strategy} improved: {best_time} ms")
        
        return best_kernel
```

## Economy of Attention: The Debug Loop

Your insight about the **edit→compile→test→review→reflect→plan** loop is the heart of the system. The key is **who does the work at each stage**:

| Stage | Who | Cost | Key Optimization |
|-------|-----|------|----------------|
| Edit Code | LLM (guided by human PRD) | Medium | Good prompts, clear goals |
| Compile | Compiler (automatic) | Cheap | Lint/typecheck in CI |
| Test Run | Test harness (automatic) | Cheap | Parallel, cached |
| Review Result | **Automatic** for text/numbers | Cheapest | NaN checks, RMSE |
| Review Result | **LLM** for pattern analysis | Medium | Parse logs, flag anomalies |
| Review Result | **Human** for plots/trajectories | Costly | Only when LLM flags |
| Reflect | LLM (diagnose failure) | Medium | Error log → root cause |
| Plan Next | LLM (strategy) | Medium | Based on failure mode |

**The agent must produce debug artifacts that are self-documenting**:

```
/debug_artifacts/
├── 2024-06-14_scan_OH_bond/
│   ├── scan_profile.png          # Human reviews this
│   ├── rmse_analysis.txt         # LLM reads this
│   ├── topology_check.json       # Automatic
│   └── conservation_drift.png    # Human reviews if drift > threshold
```
## Implementation Roadmap

**Week 1: The Loop Skeleton**
- Implement Ralph loop with `plan.md`, `CHANGELOG.md`, `MEMORY.md`
- Connect to Aider or Claude Code for file editing
- Add basic compile/run tools

**Week 2: Test Registry**
- Register your physics invariants as testable properties
- Implement Level 1-3 tests (compile, topology, parity)
- Cache reference outputs (golden I/O files)

**Week 3: LLM Analysis Layer**
- Add LLM-based qualitative checks (Level 4-5)
- Generate plots and reports automatically
- Build human review queue (only flagged items)

**Week 4: Multi-Agent & Performance**
- Spawn subagents for parallel platform migration
- Add performance optimization loop
- Integrate DFT reference pipeline (cached, expensive)

**Week 5: Self-Improvement**
- After successful migrations, extract patterns into `SKILL.md`
- Build procedural memory of what works for your specific codebase
- Tune iteration budgets based on task complexity

The framework you described is essentially a **physics-aware, cost-sensitive agentic testing system**. The research is converging on exactly this architecture — the key differentiator is your domain-specific knowledge of what constitutes "reasonable" in computational chemistry, which must be encoded in the agent's skills and memory.

![Agentic Testing Framework](sandbox:///mnt/agents/output/agentic_testing_framework.png)

---

# Claude Sonet 4.6

This is a really well-thought-out framework. Let me elaborate it into a structured system you can actually implement.Click any level for a deeper dive. Now let me elaborate the full framework in detail.

## The debug loop economy

Your loop is `edit → compile → test_run → review_result → reflect → plan`. Every iteration has a cost in time and attention. The framework above exists to **fail fast at the cheapest possible level** — you never want to wait 10 minutes for a rigid scan to tell you something a NaN check could have caught in 2 seconds.

The agent should always run levels sequentially and abort upward on first failure.

## Level 0 — Compile & crash

The agent checks: did it compile? Did it run without segfault, OpenCL device hang, or infinite loop? For OpenCL specifically: did `clCreateKernel`, `clSetKernelArg`, and `clEnqueueNDRangeKernel` all return `CL_SUCCESS`? Did the device come back after the kernel finished?

This is fully automatic, zero physics knowledge needed, and should block everything else.

## Level 1 — Numerical sanity

The agent reads the raw text output and checks a handful of invariants you provide in a `PHYSICS_CONSTRAINTS.md` or inline in the test spec:

**Automatic checks (no LLM needed):**
- No NaN or Inf anywhere in energy, force tensor, stress
- Total energy in a plausible range (e.g. for a water dimer: −0.5 to +0.5 eV)
- Forces pass a finite-difference consistency check: $F_i \approx -(E(+h) - E(-h)) / 2h$ for each atom/direction
- NVE trajectory: total energy drift < 0.01% over 100 steps
- Virial stress trace gives reasonable pressure

The crucial insight you identified: **the agent should infer physical ranges from the system type**. You give it `system: water_dimer, method: TIP4P` and it should know bond lengths are 0.95 Å, H-O-H angle ~104°, dimer binding ~−0.21 eV. This is exactly the kind of thing a system prompt in a `PHYSICS_CONSTRAINTS.md` can encode.

## Level 2 — Topology & mapping

This is underappreciated and is often where ports silently break. The agent diffs your neighbor list, angle list, and dihedral list against a reference (the original Fortran, or a Python reference build) and checks:

- All pairs within cutoff are present — no missing neighbors
- No pair appears twice
- Empty slots are padded with −1 (not 0, which is a valid atom index!)
- Angle entries `[i, j, k]` point to the same physical atoms in both implementations
- OpenCL workgroup padding doesn't leave uninitialized memory in an active slot

This is an automatic diff — no physics reasoning needed. It just compares index arrays.

## Level 3 — Qualitative physics

This is where you need the LLM, not just a comparator. The test cases are your "canonical molecules" — ordered from simplest to most demanding:

**For QM methods:** H₂ (minimal, s-only) → N₂ (triple bond, homo-nuclear) → CO (polar, LUMO on C) → H₂O (angle, dipole) → benzene (π delocalization)

**For classical forcefields:** H₂O (geometry) → NH₃ (inversion, pyramidal) → formaldehyde (sp² C) → benzene (aromatic planarity, UFF torsions or your π model) → carboxylic acid dimer (H-bond geometry)

**For surface interactions:** single ion/atom scan with LJ+Coulomb → H₂O on NaCl (your exact example: O toward Na⁺, H toward Cl⁻)

For each, the agent:
1. Runs the calculation
2. Reads text output
3. Checks the qualitative assertions: "Is the O–Na distance less than O–Cl? Is it between 2.0 and 3.5 Å? Is the binding energy negative?"
4. Writes a 3-line verdict: PASS/FAIL with reasoning

The LLM is doing physics reasoning here, not just pattern matching — that's what makes it medium-cost but high-value.

## Level 4 — Parity & rigid scans

For same-method ports: rigid bond/angle/dihedral scans produce 1D energy profiles. The agent computes RMSE between the two implementations. For a faithful port, RMSE should be near machine epsilon times the energy scale.

For new-method comparisons: you don't want RMSE — you want **qualitative shape agreement**. The agent checks: same minimum location (±0.1 Å)? Same barrier height order of magnitude? Correct sign of asymmetry (CO vs OC on surface)?

The agent should always produce a simple SVG or HTML plot of the scan — one curve per method, labelled. You spend 5 seconds looking at it. This is the "cheap human" step: you're not debugging, you're approving.

**For the agentic loop, the plot path is the key artifact.** The agent should always end its iteration by printing: `Scan plot written to: debug/scan_bondlength_OH_iter3.svg` — explicit, datestamped, never overwriting previous iterations.

## Level 5 — Performance

This only runs after Level 3 passes — you never optimize broken physics.

The agent sweeps:
- Workgroup sizes: 64, 128, 256, 512
- Memory layout: AoS (atom-of-structures) vs SoA (structure-of-arrays) — critical for coalesced reads in OpenCL
- Reduction strategy for recoil forces: atomic add (simple, slow at high occupancy) vs assembly/local memory reduction (fast, complex)
- Single vs double precision (for your GPU work this matters enormously — single is often 8–16× faster on consumer hardware)

The agent reads the profiler output (or `cl_event` timing), builds a table, and flags the winning configuration. Roofline comparison (achieved GFlop/s vs theoretical peak) tells you whether you're memory-bound or compute-bound — which determines where to focus next.

## The `PHYSICS_CONSTRAINTS.md` file

This is the single most important artifact you can write. It's the oracle that both the LLM and the automatic test runner share. Structure:

```markdown
## System: water_monomer
method_family: [TIP3P, TIP4P, SPC/E, our_pi_model]
bond_OH: [0.94, 1.02] Å
angle_HOH: [100, 110] deg  
dipole: [1.7, 1.9] D
energy_range: [-0.5, 0.1] eV  # relative to dissociated atoms

## System: NaCl_surface + H2O
O_near_Na: true          # O must be the closest water atom to Na+
H_toward_Cl: true        # H must point within 45° of Cl- 
adsorption_energy: [-0.8, -0.2] eV
O_Na_distance: [2.0, 3.2] Å
```

The agent reads this file at the start of every iteration. New methods don't need to match exact numbers — they need to satisfy the qualitative assertions.

## How the agentic loop uses all of this

The loop's state lives in `debug/progress.md` and git. Each iteration:

1. Agent reads `progress.md` + `git log --oneline -10`
2. Runs Level 0–2 (all automatic, seconds)
3. If pass: runs Level 3 (LLM analysis, up to a minute)
4. If pass: generates Level 4 scan plots, writes them to `debug/`
5. Updates `progress.md`: what passed, what failed, what changed
6. Commits with a descriptive message
7. If all pass: optionally runs Level 5 performance sweep
8. Loop exits when all levels pass, or flags for human review

The **economy of attention** you mentioned is encoded directly: you only look at plots when Level 3 passes. You only do a human debug session when the LLM says "I see the energy profile but the minimum is at 1.4 Å which seems too short for O–Na; please check atom index mapping in the angle list."

Want me to draft the actual `PHYSICS_CONSTRAINTS.md` template for your NaCl/water forcefield work, or sketch the `progress.md` structure and loop script?

---

# ChatGPT (5.5?)

I think you are converging on a very important idea: **scientific agentic loops are not primarily a coding problem; they are a problem of designing a hierarchy of verifiers.**

The agent should not ask a human after every iteration. The whole point is to have a very cheap inner loop that can run hundreds of times:

```
edit code
    |
compile
    |
run minimal tests
    |
automatic checks
    |
LLM interpretation
    |
reflect
    |
modify code
    |
repeat
```

Only after a significant milestone should the agent produce artifacts for human review.

A good way to think about it is **cost of attention**:

| Level | Who evaluates                | Cost                 | Frequency       |
| ----- | ---------------------------- | -------------------- | --------------- |
| L0    | Compiler/runtime             | microseconds–seconds | every iteration |
| L1    | Numerical scripts            | milliseconds–seconds | every iteration |
| L2    | LLM reading logs             | seconds              | occasionally    |
| L3    | LLM generating plots/reports | tens of seconds      | milestone       |
| L4    | Human visual inspection      | minutes–hours        | rarely          |

The goal is to push as much as possible down to L0–L2.

# I. Level 0 — Build and runtime survival

The first criterion is simply:

**"Does the program exist and execute?"**

Checks:

### Compilation

```
✓ no syntax errors
✓ no linker errors
✓ OpenCL kernels compile
✓ no missing symbols
```

### Runtime stability

```
✓ no segmentation fault
✓ no CUDA/OpenCL device error
✓ no illegal memory access
✓ no timeout
✓ no memory leak
```

### Basic numerical health

For every output array:

```
energy:
    finite?
    not NaN?
    not ±inf?

forces:
    finite?
    max(|F|) < F_max

coordinates:
    not exploded
```

A useful idea is to maintain expected physical ranges:

Example:

```
bond length:
    0.5 Å < r < 5 Å

energy:
    -1000 eV < E < 1000 eV

force:
    |F| < 1000 eV/Å
```

The agent can infer rough ranges, but the user should be able to override them.

# II. Level 1 — Topology and bookkeeping validation

In computational chemistry, a huge fraction of bugs are not numerical, but **mapping errors**.

Especially after:

```
Fortran -> C++
C++ -> OpenCL
CPU -> GPU
AoS -> SoA memory layout
```

Typical failures:

```
atom 15 becomes atom 16

neighbor list shifted by one

uninitialized padding

wrong periodic image

angles point to wrong atoms
```

These should be explicitly tested.

## Neighbor list checks

Example:

Reference:

```
atom 5:

neighbors:
[1, 3, 8, -1]
```

GPU:

```
atom 5:

neighbors:
[1, 8, 3, 52791]
```

Automatic checker should report:

```
FAIL

atom 5:
  expected:
    [1,3,8,-1]
  got:
    [1,8,3,52791]

problem:
  wrong order
  uninitialized padding
```

---

## Geometric consistency

For every bond:

```
distance(i,j) = distance(j,i)
```

For angles:

```
angle(i,j,k) uses existing bonds
```

For periodic systems:

```
minimum image convention consistent
```

These are extremely cheap and should run constantly.

---

# III. Level 2 — Single-point numerical parity

This is the most important test for direct reimplementation.

Example:

```
Fortran Fireball
       |
       | H2O geometry
       |
       v
OpenCL FireCore
```

Compare:

```
total energy

atomic forces

charges

orbital energies

Hamiltonian matrix
```

Metrics:

```
Energy RMSE

Force RMSE

Force direction cosine

max error
```

For forces, direction is often more important than magnitude:

```
cos(F_ref, F_test) > 0.999
```

---

But one should distinguish:

## Exact parity

For the same algorithm:

```
Fortran double precision

vs

OpenCL single precision
```

you may require:

```
Energy difference < 1e-5 eV

Force error < 1e-4 eV/Å
```

## Equivalent physics

Different formulation:

```
UFF

vs

new pi-orbital force field
```

then parity is impossible.

The question becomes:

```
Do they predict the same chemistry?
```

# IV. Level 3 — Energy landscape comparison

This is probably one of the most valuable tests in your field.

Do not test one geometry.

Test a controlled scan.

Examples:


## Bond scan

```
H-H distance:

0.5 Å
0.6 Å
...
3.0 Å
```

Compare:

```
E(r)
F(r)
```

Metrics:

```
RMSE(E)

minimum position

equilibrium distance

curvature near minimum
```

## Angle scan

Example:

```
H-O-H angle:

60° to 180°
```

Compare:

```
energy profile
```

## Dihedral scan

For molecules with sp² chemistry:

```
benzene
pyridine
formaldehyde
formic acid
oxalate
```

Check:

```
planarity

torsion barriers
```

## Surface scans

For NaCl:

Single atom:

```
Na+ over Na site

Na+ over Cl site

Na+ over bridge
```

Then:

```
E(x,y,z)
```

For water:

Expected qualitative behavior:

```
Water over Na+:

O points down

distance ~2-3 Å


Water over Cl−:

H atoms point down
```

The exact adsorption energy is not important at this stage.

# V. Level 4 — Geometry optimization tests

This is the first truly "chemistry-aware" level.

Single point tells:

```
Where are forces pointing?
```

Optimization tells:

```
Where does the system want to go?
```

Simple molecules:

## Quantum methods

Progression:

```
H2
 ↓
N2
 ↓
CO
 ↓
H2O
```

Each introduces more complexity:

H2:

```
two atoms
bond formation
```

N2:

```
multiple bond
sp basis
```

CO:

```
heteronuclear polarization
```

H2O:

```
angular geometry
lone pairs
```

## Classical force fields

Progression:

```
H2O

NH3

methane

formaldehyde

formic acid

benzene

pyridine

carboxylic acid dimer
```

Tests:

```
bond lengths

bond angles

planarity

hydrogen bonds

pi alignment
```

The result should not be a huge trajectory.

The agent should summarize:

```
Optimization summary:

initial energy:
 -12.5 eV

final energy:
 -15.2 eV

steps:
 56

max force:
 0.02 eV/Å

Geometry checks:

H-O-H angle:
 104.6°
 expected:
 100-110°
 PASS

O-H bond:
 0.97 Å
 expected:
 0.9-1.1 Å
 PASS
```

# VI. Level 5 — LLM-assisted interpretation

This is where an LLM is valuable.

It can read:

```
logs

tables

RMSE reports

optimization summaries
```

and reason:

> The energy is correct but forces differ strongly, suggesting a sign error or missing derivative term.

or:

> The minimum is shifted by 0.5 Å, suggesting incorrect equilibrium bond length or missing electrostatic contribution.

or:

> The OpenCL implementation matches small molecules but fails larger systems, suggesting race conditions or memory corruption.

The LLM should not inspect raw million-line logs.

The code should produce **compact diagnostic summaries**.

# VII. Human visual review

Only when the agent has already reduced the problem.

Good artifacts:

```
plots/
    H2_scan.png
    NaCl_xy_scan.png
    force_error.png

structures/
    H2O_relaxed.xyz
    benzene.xyz

fields/
    density.cube

reports/
    test_summary.md
```

The agent should say:

```
Review required.

Interesting artifacts:

plots/NaCl_scan.png

Reason:
The minimum is displaced by 0.8 Å from reference.

Hypothesis:
periodic image convention is wrong.
```

The human should never search through thousands of files.

# VIII. Performance optimization loops

This is a slightly different agent.

The goal is not "correctness" but:

```
maximize performance
subject to
accuracy >= threshold
```

So the loop is:

```
modify kernel
        |
compile
        |
run correctness benchmark
        |
if failed:
    reject
        |
else:
    run timing benchmark
        |
measure:
    GFLOPS
    memory bandwidth
    occupancy
    scaling
        |
keep best variant
```

For OpenCL this includes exploring:

```
memory layout:
    AoS vs SoA

precision:
    double vs float

work-group size:
    32, 64, 128, 256

local memory caching

vector types:
    float4, float8

loop unrolling

avoiding atomics

reordering computation

reducing global memory traffic
```

# IX. The key insight: "scientific test pyramid"

I would formulate the whole system as:

```
              Human scientific judgment
                       ▲
                       |
              Visual inspection
                       ▲
                       |
              LLM physical reasoning
                       ▲
                       |
              Energy landscape tests
                       ▲
                       |
              Single-point parity
                       ▲
                       |
              Topology checks
                       ▲
                       |
              Numerical sanity
                       ▲
                       |
              Compilation/runtime
```

The lower levels run **thousands of times**.

The upper levels run **only when necessary**.

I think this hierarchy is actually the missing piece that makes AI agents useful for computational science. Current software agents are very good at "run pytest until green", but in scientific computing **the notion of correctness is much richer**: symmetry, conservation laws, qualitative chemistry, energy landscapes, and physical plausibility. Designing these tests in a reusable way (e.g. `TESTS_PHYSICS.md` plus a library of small benchmark systems) could easily become one of the most valuable assets in your entire code ecosystem.

---

# Kimi 2.6 devin















# Summary of LLM Responses (lines 1335–3901)

The document contains responses from **Kimi/devin** (which is actually my prior response to you), **Gemini**, **DeepSeek-V4-flash**, **Mistral**, **Kimi 2.6 chat**, **Claude Sonet 4.6**, and **ChatGPT** on the testing hierarchy for scientific agentic loops.

---

## Where the Models Agree

**1. Three-tier attention economy**
Universal consensus: automatic scripts (cheap, fast) → LLM text/log analysis (medium cost) → human visual review (expensive, minimize). No model suggests human review for compilation errors or NaN checks.

**2. Hierarchical gating**
All models structure tests as a pyramid where cheap checks *gate* expensive ones. If compilation fails, you do not run scans. If topology is wrong, you do not check physics.

**3. Topology/mapping is critically underappreciated**
Every model identifies neighbor-list indexing, padding values, angle atom mappings, and 1-based vs 0-based shifts as the dominant source of silent bugs in cross-platform ports. All agree this should be tested early and automatically.

**4. Rigid scans > single points**
All agree that comparing energy/force profiles along bond lengths, angles, or surface distances is far more discriminating than a single geometry comparison. All suggest generating plots as artifacts.

**5. Performance optimization must preserve correctness**
All treat performance as a separate loop or final level with a hard requirement: the accuracy suite (levels 0–4) must still pass after every kernel optimization.

**6. Canonical test system progression**
All converge on the same sequence: H₂ (minimal) → N₂ (sp, multiple bonds) → CO (polar, heteronuclear) → H₂O (angles, lone pairs) → NH₃ (tetrahedral) → formaldehyde (sp²) → benzene/pyridine (π-systems) → carboxylic dimer (H-bonds) → NaCl surface (ionic/nonbonded).

---

## Where the Models Disagree

**Level count & ordering** (mostly cosmetic, but reveals priority):
- **Kimi/devin**: 6 levels (compile → sanity → topology → parity → qualitative → performance)
- **Mistral**: 8 levels (adds rigid scans and relaxation as separate levels)
- **DeepSeek**: 9 levels (0–8, with explicit "different method" skip logic)
- **ChatGPT**: 9 levels (0–8) framed as a scientific test pyramid
- **Claude**: ~5 core levels but explicitly separates LLM interpretation and human review as orthogonal overlays

**Topology placement**:
- **Kimi/devin, Mistral, Claude, ChatGPT**: Topology early (Level 1–2)
- **DeepSeek**: Topology at Level 4, after physical range sanity — slightly later

**Qualitative physics ordering**:
- **Kimi/devin**: Before relaxation (Level 4)
- **Mistral/DeepSeek**: After relaxation (L7/L8) — you need optimized geometry before judging qualitative binding
- **Claude**: Collapses qualitative and parity into adjacent levels

**Architecture**:
- **Kimi 2.6 chat** pushes a **multi-agent swarm** (translator, test_generator, parity_checker, performance_optimizer, human_liaison) — heavy orchestration.
- **Gemini, Claude, ChatGPT** favor a **single lightweight loop** with external state in markdown/git — simpler, more robust.
- **DeepSeek** sits in between: single loop but with a rich tool registry.

**Physical range inference**:
- **Mistral, Kimi 2.6 chat**: LLM infers ranges once, caches in `config.yaml` / `MEMORY.md`
- **Claude**: User writes `PHYSICS_CONSTRAINTS.md` as the shared oracle; agent reads it but does not generate it
- **DeepSeek**: YAML file per project, user-provided preferred, LLM-inferred with confirmation

---

## Innovative Ideas

| Model | Innovation |
|-------|-----------|
| **Gemini** | Explicit `AgenticValidationHarness` class with `execute_tier1/2/3` methods. Treats OpenCL device freezes as catastrophic failures requiring kernel redistribution. Suggests `generate_p5js_coordinate_view` for interactive human review. |
| **DeepSeek** | **Method-specific tolerance matrix** (JSON): `{"method": "my_new_dft", "rmse_energy_eV_per_atom": 0.05}`. This cleanly handles "different method, same physics." Also `get_next_cheapest_test()` — the loop queries which test to run next based on cost and state. |
| **Mistral** | Most concrete code for every level. Explicit cost symbols (`$` to `$$$$`). Parallelize L1–L4. `debug/REVIEW.md` documents every human decision for future agent learning. Comprehensive test catalog with expected values (e.g., H₂O dimer binding ~5 kcal/mol). |
| **Kimi 2.6 chat** | **Ralph loop** with `plan.md` + `CHANGELOG.md` + `MEMORY.md` as external state. **Golden I/O files** — recorded reference inputs/outputs that serve as oracles. **Property-based testing** (Hypothesis-style randomized invariants). 5-week implementation roadmap. |
| **Claude Sonet 4.6** | **`PHYSICS_CONSTRAINTS.md`** as the single most important artifact — a shared oracle for both automatic runner and LLM. Finite-difference force consistency check (`F_i ≈ -(E(+h)−E(−h))/2h`). Datestamped plot paths that never overwrite. `debug/progress.md` for loop state. |
| **ChatGPT** | Most crystalline formulation: **Scientific Test Pyramid** visualization. Distinguishes "Exact parity" vs "Equivalent physics" explicitly. **Force direction cosine** (`cos(F_ref, F_test) > 0.999`) as more important than force magnitude. LLM should read *compact diagnostic summaries*, never raw million-line logs. |

---

## Best Path Forward (My Assessment)

Given your actual FireCore conventions (`AGENTS.md`, `run.sh` scripts, parity tests in `tests/tEFF/`, CPU-dotted/GPU-solid plot convention) and your stated preferences (minimal dependencies, loud failures, reusable code, no ad-hoc scripts), here is what I recommend:

### 1. Adopt Claude's `PHYSICS_CONSTRAINTS.md` as the Central Oracle
This is the highest-leverage single file. It lives alongside your existing `AGENTS.md` hierarchy and encodes your domain knowledge in a format both automatic scripts and LLMs can read. Example:
```markdown
## System: water_monomer
method_family: [TIP4P, our_pi_model]
bond_OH: [0.94, 1.02] Å
angle_HOH: [100, 110] deg
energy_range: [-0.5, 0.1] eV

## System: NaCl_001 + H2O
O_near_Na: true
O_Na_distance: [2.0, 3.2] Å
```

### 2. Implement a 6-Level Python Harness in `tests/common/`
Do not create a new framework. Write a reusable `tests/common/agentic_harness.py` with classes matching your existing test style:
- `Level0RuntimeTest` — compile, timeout, segfault, NaN/Inf
- `Level1TopologyTest` — neighbor map diff, padding check, index bijection
- `Level2ParityTest` — single-point RMSE, force cosine
- `Level3ScanTest` — rigid bond/angle/surface scans, plot generation
- `Level4RelaxationTest` — geometry optimization, convergence, qualitative bounds
- `Level5PerformanceTest` — kernel timing, workgroup sweep (only if L0-L4 pass)

Each level returns `(passed: bool, report: str, artifacts: list[Path])`.

### 3. Use Golden I/O Files (Kimi 2.6 chat idea)
For every reference implementation (Fortran/C++), store serialized reference outputs in `tests/ref_data/<system>/`:
```
tests/ref_data/H2O_TIP4P/
  ├── energy.npy
  ├── forces.npy
  ├── scan_OH_bond.npy
  └── neigh_list.npy
```
The harness loads these automatically. You generate them once manually; the agent never regenerates them.

### 4. Artifact Naming Convention (consistent across all models)
All models converged on this. Use it:
```
debug/2024-06-14_H2O_OpenCL_port/
  ├── L0_compile.log
  ├── L1_topology_diff.txt
  ├── L2_parity_report.md
  ├── L3_scan_OH_bond.png       # CPU dotted, GPU solid
  ├── L4_relax_H2O_NaCl.xyz
  ├── L4_relax_geometry.png
  └── REPORT.md                 # LLM-generated summary
```

### 5. Performance Loop is Separate, with Regression Guard
Following DeepSeek and ChatGPT: the performance agent must run the full L0–L3 harness after *every* kernel variant. If correctness fails, reject immediately. Use OpenCL event profiling for timing, not wall-clock where possible.

### 6. Start Minimal — One Module, One Test
Do not build the multi-agent swarm. Pick one concrete task: e.g., validate the OpenCL `evalMMFFf4_local` kernel against the C++ reference. Write the harness for that one case. Once it works, generalize to `tests/common/`. This respects your "small, testable, step-by-step changes" rule.

### Concrete First File
Create `tests/common/scientific_harness.py` with:
```python
class ScientificTestHarness:
    def __init__(self, ref_dir, physics_constraints_path):
        self.ref = ReferenceData(ref_dir)
        self.constraints = PhysicsConstraints(physics_constraints_path)
    
    def run_level(self, level: int, target_func) -> TestResult:
        ...
```

Then a thin script in `tests/tMMFF/test_agentic_parity.py` that calls it. Run it manually once to generate golden I/O, then hand the loop to the agent.

This gives you a working, reusable system within one day that aligns with all seven models' recommendations while respecting FireCore's existing conventions and your personal rules about minimal dependencies and loud failures.
