
# USER

Do you have any idea how to solve this problem with winsurf language server eating all my memory? I want some genral hard solution which basically kills it immadiately (not windsurf just ht elanguage server process) when it try to use more memory than certain level. It is better when it crashes than wen it freeze my computer so I cannot do anything and must turn it of by power button
---
**Summary**
Language server suddnely peak to ~20GB RAM usage freezing my desktop (can hardly move mouse) and then crashes
**What’s happening**
This happend repeately. Several times a day, whenever I open lager number of files lin large repo (pyhon, OpenCL, markdown). I know my repo is large. But can I limit the memory usage. I cannot allow single process to take 20GB of RAM alone (I have 32GB in total)
I tried this recomented by Windsurf suppoer robot, not sure what it should do but it did not helped
```
prokop@GTX3090:~$ echo "fs.inotify.max_user_watches=524288" | sudo tee -a /etc/sysctl.conf
fs.inotify.max_user_watches=524288
prokop@GTX3090:~$ echo "fs.inotify.max_user_instances=1024" | sudo tee -a /etc/sysctl.conf
fs.inotify.max_user_instances=1024
```
**Expected behavior**
The language server should either use less memory (keep belo some trashold e.g. 5GB) or crash, but not freeze my system. I do not mind if it crash, but I hate when it freze my system so I cannot eve move a mouse.
**Environment**
- OS: Linux Ubutu 2024 LTS
- Device: Desktop
- Windsurf version: 2.2.17
**Errors / Logs / Diagnosis**
this is what I see in Taks manager:
/usr/share/windsurf/resources/app/extensions/windsurf/bin/language_server_linux_x64 --api_server_url https://server.codeium.com --run_child --enable_lsp --extension_server_port 39689 --ide_name windsurf --random_port --inference_api_server_url https://inference.codeium.com --database_dir /home/prokop/.codeium/windsurf/database/9c0694567290725d9dcba14ade58e297 --enable_index_service --enable_local_search --search_max_workspace_file_count 1000 --indexed_files_retention_period_days 30 --workspace_id file_home_prokop_git_FireCore --sentry_telemetry --sentry_environment stable --codeium_dir .codeium/windsurf --extensions_dir /home/prokop/.windsurf/extensions --parent_pipe_path /tmp/server_8970f235f38775be --windsurf_version 2.2.17 --stdin_initial_metadata --detect_proxy=false

---

# DeepSeek-V4

You want a **hard solution** that kills the Windsurf language server process immediately if it exceeds a memory threshold (e.g., 5 GB), before it locks up your 32 GB system.  
The language server is a separate binary:  
`/usr/share/windsurf/resources/app/extensions/windsurf/bin/language_server_linux_x64`

The most reliable and immediate method is to use **Linux cgroups v2 memory limits** – the kernel will kill the process the moment it goes over the limit. No polling delay, no system freeze.  
Below I provide a step‑by‑step solution that:

1. Creates a cgroup with a memory limit.
2. Automatically moves any running language server process into that cgroup.
3. Starts a lightweight systemd service that keeps the cgroup enforced even when the language server restarts.

## Step 1 – Create a cgroup with a memory limit

Ubuntu 24.04 uses **cgroup v2** by default. We’ll create a child cgroup under the user slice.

```bash
# Create a cgroup named "limit_lsp"
sudo mkdir /sys/fs/cgroup/limit_lsp

# Set memory limit to 5 GB (use 5G, 5000M, or whatever you prefer)
echo "5G" | sudo tee /sys/fs/cgroup/limit_lsp/memory.max

# (Optional) also limit swap usage to prevent swapping out – same value is fine
echo "5G" | sudo tee /sys/fs/cgroup/limit_lsp/memory.swap.max
```

> **How it works:** When a process inside this cgroup tries to allocate more than 5 GB of memory, the kernel immediately kills it (SIGKILL). The rest of Windsurf continues running.

## Step 2 – Automatically put the language server into that cgroup

The language server is spawned by Windsurf. We need a daemon that repeatedly picks up its PID and writes it to `cgroup.procs`.  
Create a simple script `/usr/local/bin/limit_lsp.sh`:

```bash
#!/bin/bash
CGROUP="/sys/fs/cgroup/limit_lsp"

while true; do
    # Find any language server process (the one with --run_child flag)
    pid=$(pgrep -f "language_server_linux_x64.*--run_child" | head -1)
    if [[ -n "$pid" ]]; then
        # Check if the process is already in our cgroup
        if ! grep -q "^$pid$" "$CGROUP/cgroup.procs" 2>/dev/null; then
            echo "$pid" | sudo tee "$CGROUP/cgroup.procs" >/dev/null
            echo "Moved PID $pid into memory‑limited cgroup"
        fi
    fi
    sleep 2   # check every 2 seconds (adjust if you want even faster)
done
```

Make it executable:

```bash
sudo chmod +x /usr/local/bin/limit_lsp.sh
```

## Step 3 – Run the script as a systemd user service (auto‑start & resilient)

Create a systemd **user** service (so it runs under your login, not needing root for every echo).

First, make the script usable without password for the `tee` command. Add a sudoers entry (safe because only that script writes to that one cgroup file):

```bash
echo "your_username ALL=(ALL) NOPASSWD: /usr/bin/tee /sys/fs/cgroup/limit_lsp/cgroup.procs" | sudo tee /etc/sudoers.d/limit-lsp
```

Replace `your_username` with your actual login name.

Now create the service file `~/.config/systemd/user/limit-lsp.service`:

```bash
mkdir -p ~/.config/systemd/user
cat > ~/.config/systemd/user/limit-lsp.service <<EOF
[Unit]
Description=Limit Windsurf language server memory via cgroup
After=default.target

[Service]
ExecStart=/usr/local/bin/limit_lsp.sh
Restart=always
RestartSec=5
StandardOutput=journal

[Install]
WantedBy=default.target
EOF
```

Enable and start the service:

```bash
systemctl --user daemon-reload
systemctl --user enable limit-lsp.service
systemctl --user start limit-lsp.service
```

Check that it’s running:

```bash
systemctl --user status limit-lsp.service
```

## Step 4 – Verify it works

1. Open Windsurf and load your large repository.
2. Monitor the memory limit in another terminal:
   ```bash
   cat /sys/fs/cgroup/limit_lsp/memory.current
   ```
   Watch how big the language server grows.
3. When it hits 5 GB, the process will be **killed instantly** – you’ll see a message in the Windsurf output that the language server crashed, but your desktop remains responsive.  
   (The service will move the new language server process into the cgroup again after a restart.)

## Alternative simpler (but slower) solution – polling + kill

If you don’t want to use cgroups, you can use a **bash script** that runs every 5 seconds and kills the process if it exceeds 5 GB.  
This is less elegant (lag between detection and kill) but very easy:

```bash
#!/bin/bash
MEM_LIMIT_KB=$((5 * 1024 * 1024))   # 5 GB in KB

while true; do
    pid=$(pgrep -f "language_server_linux_x64.*--run_child")
    if [[ -n "$pid" ]]; then
        # Get RSS memory in KB (from /proc)
        mem=$(awk '/^VmRSS:/ {print $2}' /proc/$pid/status 2>/dev/null)
        if [[ -n "$mem" && $mem -gt $MEM_LIMIT_KB ]]; then
            kill -9 $pid
            echo "Killed language server (PID $pid) – used ${mem}KB > limit"
        fi
    fi
    sleep 5
done
```

Run it as a background script or a systemd service (similar to the cgroup method but without the sudoers trick).  
This will work but may let the process spike to 5.2 GB before the next check. The cgroup method is **precise and immediate**.

## Important notes

- The language server is **not** Windsurf itself – killing it only removes the advanced code intelligence (completions, indexing). Windsurf will keep running, often with a notification that the language server crashed. That’s exactly what you want.
- If you also want to **see** when it gets killed, check `journalctl --user -u limit-lsp.service` or look at kernel messages (`dmesg | tail -20`) for OOM or cgroup kill events.
- The cgroup limit applies only to the language server process. Windsurf’s main process and other extensions are unaffected.
- After a system reboot, the cgroup `/sys/fs/cgroup/limit_lsp` is **gone**. You need to re‑create it: add a simple startup service or just run the `mkdir` + `echo` commands once after boot. Alternatively, use a systemd **system** service (not user) to create the cgroup at boot. I can give you that snippet if you want persistence.

**Recommendation:** Use the cgroup method. It’s the “hard” immediate solution you asked for and will completely prevent your computer from freezing.

---

# Grok

**The best hard solution on modern Ubuntu (cgroup v2) is to enforce a memory.max limit on the language server process.** When it hits the limit, the kernel's OOM killer will immediately kill processes inside that cgroup (usually the leaky one), preventing system-wide freeze. This is exactly what you want: crash the LS instead of locking up your desktop.

### Quickest Manual Method (for testing)
1. Install tools if needed: `sudo apt install cgroup-tools`

2. Create a cgroup:
   ```bash
   sudo mkdir -p /sys/fs/cgroup/windsurf-ls
   ```

3. Set a hard memory limit (e.g., **5GB** — adjust as needed; use `5G` or bytes):
   ```bash
   echo "5G" | sudo tee /sys/fs/cgroup/windsurf-ls/memory.max
   # Optional: softer throttling point before hard kill
   echo "4G" | sudo tee /sys/fs/cgroup/windsurf-ls/memory.high
   ```

4. Find the language server PID (the long binary you posted):
   ```bash
   pgrep -f language_server_linux_x64
   ```

5. Move it into the cgroup:
   ```bash
   echo <PID> | sudo tee /sys/fs/cgroup/windsurf-ls/cgroup.procs
   ```

Repeat step 4-5 whenever it restarts (it will probably restart automatically).

To make it persistent and automatic, use a wrapper or systemd scope.

### Better: Launch via systemd-run (transient scope with limits)
This is cleaner for processes launched by an IDE:

```bash
systemd-run --scope --user -u windsurf-ls.scope \
  -p MemoryMax=5G \
  -p MemoryHigh=4G \
  -- /usr/share/windsurf/resources/app/extensions/windsurf/bin/language_server_linux_x64 [all the original args]
```

But since Windsurf launches it, you may need to create a wrapper script.

### Recommended Persistent Wrapper Script
Create `~/bin/windsurf-ls-wrapper` (make executable):

```bash
#!/bin/bash
# Wrapper for Windsurf language server with hard memory limit

CGROUP="/sys/fs/cgroup/windsurf-ls"

# Ensure cgroup exists and has limit
if [ ! -d "$CGROUP" ]; then
    sudo mkdir -p "$CGROUP"
    echo "5G" | sudo tee "$CGROUP/memory.max" > /dev/null
    echo "4G" | sudo tee "$CGROUP/memory.high" > /dev/null
fi

# Run the real binary and move it to cgroup
exec "$@" &
PID=$!
echo "$PID" | sudo tee "$CGROUP/cgroup.procs" > /dev/null
wait $PID
```

Then make it executable: `chmod +x ~/bin/windsurf-ls-wrapper`

You'd need to configure Windsurf (if possible) to use this wrapper instead of the direct binary, or use `cgexec`:

```bash
sudo cgcreate -g memory:windsurf-ls   # if using cgtools
cgexec -g memory:windsurf-ls /path/to/language_server_linux_x64 ...
```

### Other Quick Protections
- **System-wide soft guard**: Add to `/etc/security/limits.conf` (or a file in `/etc/security/limits.d/`):
  ```
  * hard as 6000000   # ~6GB virtual memory per process (rough)
  ```
  Then log out/in. This uses ulimit on address space.

- Monitor & auto-kill script (as backup):
  A simple loop with `pmap` or `ps` that kills if RSS > threshold (many examples online).

### Also Address the Root Cause
This is a known memory leak in Windsurf/Codeium language server on large repos.

- Create a `.codeiumignore` file in your repo root to exclude big/irrelevant dirs (node_modules, build, .git, large data files, etc.).
- Work in subfolders instead of the entire massive repo root.
- Check Windsurf settings for indexing limits or disable heavy features (local search, etc.).
- Update to latest version if available.

The cgroup `memory.max` approach gives you the immediate hard kill you asked for — the process dies cleanly (from the kernel's perspective) instead of OOMing the whole desktop. Start with 4-6G depending on your tolerance.

---

# Implementation Report (May 13, 2026)

## What Was Implemented

Following the DeepSeek-V4 cgroup solution, the persistent memory limit was set up on the system:

### Files Created
1. **`/usr/local/bin/limit_lsp.sh`** - Daemon script that auto-detects and moves the Windsurf language server into the memory-limited cgroup
2. **`~/.config/systemd/user/limit-lsp.service`** - Systemd user service that auto-starts the daemon script
3. **`/etc/sudoers.d/limit-lsp`** - Sudoers entry allowing passwordless writes to the cgroup.procs file
4. **`/usr/local/bin/create_limit_lsp_cgroup.sh`** - Boot-time script to create the cgroup directory and set memory limits
5. **`/etc/systemd/system/create-limit-lsp-cgroup.service`** - Systemd system service that runs the cgroup creation script at boot

### Cgroup Configuration
- **Cgroup path**: `/sys/fs/cgroup/limit_lsp`
- **Memory limit**: 16 GB (initially 5 GB, changed to 8 GB, then 14 GB, then 16 GB)
- **Swap limit**: 16 GB
- **Check interval**: 2 seconds (script polls every 2 seconds for new language server processes)

### Service Status
- **Service name**: `limit-lsp.service` (user service)
- **Status**: Active and running
- **Auto-start**: Enabled (starts automatically on user login)
- **Restart policy**: Always (restarts if script dies)

## Current Configuration

**Memory limit**: 16 GB (hard kill when exceeded)
**Username**: prokop (used in sudoers entry)
**Boot-time persistence**: Enabled via system service `create-limit-lsp-cgroup.service`

## How to Reconfigure

### Change Memory Limit
```bash
# Change to new value (e.g., 10G) - takes effect immediately
echo "10G" | sudo tee /sys/fs/cgroup/limit_lsp/memory.max
echo "10G" | sudo tee /sys/fs/cgroup/limit_lsp/memory.swap.max

# To make it persistent across reboots, also edit the boot-time script:
sudo nano /usr/local/bin/create_limit_lsp_cgroup.sh
# Change the "16G" values to your desired limit in both memory.max and memory.swap.max lines
```

### Change Check Interval
Edit `/usr/local/bin/limit_lsp.sh` and modify the `sleep 2` line to desired seconds.

### Monitor Current Memory Usage
```bash
cat /sys/fs/cgroup/limit_lsp/memory.current
```

### View Service Logs
```bash
journalctl --user -u limit-lsp.service -f
```

## How to Undo (Remove the Solution)

### Stop and Disable the Services
```bash
# Stop and disable user service
systemctl --user stop limit-lsp.service
systemctl --user disable limit-lsp.service

# Stop and disable boot-time system service
sudo systemctl stop create-limit-lsp-cgroup.service
sudo systemctl disable create-limit-lsp-cgroup.service
```

### Remove Files
```bash
# Remove the daemon script
sudo rm /usr/local/bin/limit_lsp.sh

# Remove the boot-time cgroup creation script
sudo rm /usr/local/bin/create_limit_lsp_cgroup.sh

# Remove the user service file
rm ~/.config/systemd/user/limit-lsp.service

# Remove the system service file
sudo rm /etc/systemd/system/create-limit-lsp-cgroup.service

# Remove the sudoers entry
sudo rm /etc/sudoers.d/limit-lsp

# Remove the cgroup (will be recreated on reboot anyway)
sudo rmdir /sys/fs/cgroup/limit_lsp
```

### Reload Systemd
```bash
# Reload user systemd
systemctl --user daemon-reload

# Reload system systemd
sudo systemctl daemon-reload
```

## Important Notes

### After System Reboot
The cgroup directory `/sys/fs/cgroup/limit_lsp` is removed on reboot. A system service `create-limit-lsp-cgroup.service` has been implemented to automatically recreate the cgroup at boot time with the configured memory limits.

**Boot-time persistence is now fully implemented.** The system service runs before the user service and ensures the cgroup exists before the daemon tries to move processes into it.

To verify the boot-time service is working:
```bash
sudo systemctl status create-limit-lsp-cgroup.service
```

### What Happens When Limit Is Exceeded
- The kernel immediately kills the language server process (SIGKILL)
- Windsurf continues running but shows a notification that the language server crashed
- The daemon script will detect the new language server process when it restarts and move it into the cgroup again
- Your desktop remains responsive (no freeze)

## Diagnostics and Verification

### Check Service Status
```bash
systemctl --user status limit-lsp.service
```

**Expected output** (example from successful implementation):
```
● limit-lsp.service - Limit Windsurf language server memory via cgroup
     Loaded: loaded (/home/prokop/.config/systemd/user/limit-lsp.service; enabled; preset: enabled)
     Active: active (running) since Wed 2026-05-13 10:41:24 CEST; 4min 58s ago
   Main PID: 9094 (limit_lsp.sh)
      Tasks: 2 (limit: 38330)
     Memory: 756.0K (peak: 2.7M)
        CPU: 3.666s
  CGroup: /user.slice/user-1000.slice/user@1000.service/app.slice/limit-lsp.service
          ├─ 9094 /bin/bash /usr/local/bin/limit_lsp.sh
          └─10987 sleep 2

May 13 10:41:24 GTX3090 systemd[2096]: Started limit-lsp.service - Limit Windsurf language server memory via cgroup.
May 13 10:41:24 GTX3090 sudo[9101]:   prokop : PWD=/home/prokop ; USER=root ; COMMAND=/usr/bin/tee /sys/fs/cgroup/limit_lsp/cgroup.procs
May 13 10:41:24 GTX3090 sudo[9101]: pam_unix(sudo:session): session opened for user root(uid=0) by prokop(uid=1000)
May 13 10:41:24 GTX3090 sudo[9101]: pam_unix(sudo:session): session closed for user root(uid=0)
May 13 10:41:24 GTX3090 limit_lsp.sh[9094]: Moved PID 6816 into memory-limited cgroup
```

**Key indicators of success:**
- `Active: active (running)` - service is running
- `enabled` - will auto-start on login
- `Memory: 756.0K` - script itself uses very little memory
- `Moved PID XXXX into memory-limited cgroup` - language server was successfully moved

### Check if Service is Enabled for Auto-Start
```bash
systemctl --user is-enabled limit-lsp.service
```
Should output: `enabled`

### Check Current Memory Usage in Cgroup
```bash
cat /sys/fs/cgroup/limit_lsp/memory.current
```
Output is in bytes. Example: `389251072` (~371 MB)

### Check Cgroup Limits
```bash
cat /sys/fs/cgroup/limit_lsp/memory.max
cat /sys/fs/cgroup/limit_lsp/memory.swap.max
```
Should output: `8589934592` (8 GB in bytes) or similar

### View Service Logs
```bash
# Recent logs
journalctl --user -u limit-lsp.service --since "5 minutes ago"

# Follow logs in real-time
journalctl --user -u limit-lsp.service -f
```

### Check if Language Server is in Cgroup
```bash
# Find language server PID
pgrep -f "language_server_linux_x64.*--run_child"

# Check if it's in the cgroup
cat /sys/fs/cgroup/limit_lsp/cgroup.procs
```

The PID from pgrep should appear in cgroup.procs if the service is working correctly.