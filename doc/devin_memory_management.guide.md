# Devin Desktop Memory Management

## Problem

The Codeium/Windsurf/Devin language server binary (`language_server`) has a known memory leak in large mixed-language monorepos (C++, Python, Rust). It can consume all available RAM (32 GB on this system) and crash the IDE.

## Solutions Implemented

### 1. IDE Settings (Devin Desktop)

**Location:** `~/.config/Devin/User/settings.json`

**Memory caps added:**
```json
"python.analysis.diagnosticMode": "openFilesOnly",
"python.analysis.indexing": false,
"C_Cpp.intelliSenseEngine": "Tag Parser",
"C_Cpp.workspaceParsingPriority": "low",
"C_Cpp.maxConcurrentThreads": 2,
"rust-analyzer.workspace.discoverPackages": false,
"rust-analyzer.cargo.sysroot": "discover",
"typescript.tsserver.maxTsServerMemory": 2048
```

**What this does:**
- Python: Only analyzes open files, no background indexing
- C++: Uses lightweight Tag Parser instead of full IntelliSense, limits to 2 threads
- Rust: Disables automatic workspace discovery
- TypeScript: Caps memory at 2 GB

### 2. `.codeiumignore`

**Location:** `/home/prokop/git/FireCore/.codeiumignore`

**Content:**
```
**/target/
**/.venv/
**/build/
**/dist/
**/__pycache__/
**/*.o
**/*.a
**/*.so
**/.git/
**/Build*/
**/Build-dbg/
**/Build-test/
```

**What this does:** Prevents the AI indexer from scanning build artifacts, object files, and cache directories.

### 3. Systemd Watchdog (Automatic Memory Killer)

**Location:** `~/.config/systemd/user/devin-memory-watchdog.service`

**Current threshold:** 10 GB (kills `language_server` if RSS > 10 GB)

**What it does:**
- Runs automatically on login (systemd user service)
- Every 15 seconds, checks if `language_server` exceeds threshold
- If yes, kills the process and logs to `/tmp/devin-watchdog.log`
- Devin Desktop automatically respawns the language server cleanly

**Status check:**
```bash
systemctl --user status devin-memory-watchdog.service
```

**View logs:**
```bash
cat /tmp/devin-watchdog.log
```

## Adjusting Thresholds

If the watchdog is too aggressive (killing too often) or too passive (memory still spikes):

### Change watchdog threshold

Edit `~/.config/systemd/user/devin-memory-watchdog.service`:
```ini
# Current: 10 GB (10485760 KB)
ExecStart=... [ "$rss" -gt 10485760 ] ...

# To change to 15 GB:
ExecStart=... [ "$rss" -gt 15728640 ] ...
```

Then reload:
```bash
systemctl --user daemon-reload
systemctl --user restart devin-memory-watchdog.service
```

### Change TypeScript server cap

Edit `~/.config/Devin/User/settings.json`:
```json
"typescript.tsserver.maxTsServerMemory": 3072  # 3 GB instead of 2 GB
```

Reload Devin Desktop (`Ctrl+Shift+P` → "Developer: Reload Window").

## Disabling the Watchdog

If you want to disable the automatic memory killer:

```bash
systemctl --user stop devin-memory-watchdog.service
systemctl --user disable devin-memory-watchdog.service
```

## Rebrand History

The IDE has rebranded multiple times, each creating a new config directory:
- Codeium → `~/.config/Codeium/`
- Windsurf → `~/.config/Windsurf/`
- Devin → `~/.config/Devin/` (current)

Settings in old directories are **not** loaded. Always edit `~/.config/Devin/User/settings.json` for current IDE.

## Quick Diagnostics

Check what's eating RAM:
```bash
ps aux --sort=-%mem | head -20
```

Check language_server memory specifically:
```bash
ps aux | grep language_server
```
