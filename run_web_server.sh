#!/usr/bin/env bash

# Usage: ./run_web_server.sh [PORT] [DOCROOT]
# Defaults: PORT=8000, DOCROOT=.

set -u # Error on undefined variables

PORT="${1:-8000}"
DOCROOT="${2:-.}"

# --- STEP 1: KILL OLD PROCESSES ---
# We use lsof to find the PID. 
# We use 'kill -9' to FORCE kill (SIGKILL). This cannot be ignored by the process.

echo "Checking port ${PORT}..."
PIDS=$(lsof -tiTCP:"${PORT}" -sTCP:LISTEN 2>/dev/null)

if [ -n "$PIDS" ]; then
  echo "⚠️  Found process(es) ${PIDS} on port ${PORT}. Nuking them..."
  kill -9 ${PIDS} 2>/dev/null || true
  
  # Wait loop: Ensure it is actually dead before proceeding
  while lsof -tiTCP:"${PORT}" -sTCP:LISTEN >/dev/null; do
    echo "Waiting for port ${PORT} to free up..."
    sleep 0.1
  done
  echo "✅ Port ${PORT} is now free."
else
  echo "✅ Port ${PORT} is clear."
fi

# --- STEP 2: START NO-CACHE SERVER ---
echo "Starting NO-CACHE web server on port ${PORT} in ${DOCROOT}"
cd "${DOCROOT}"

# We feed the python script via standard input (<<EOF) so it's all in one file.
python3 - <<EOF
import http.server
import socketserver
import sys

PORT = ${PORT}

class NoCacheHTTPRequestHandler(http.server.SimpleHTTPRequestHandler):
    def end_headers(self):
        # FORCE THE BROWSER TO NEVER CACHE
        self.send_header("Cache-Control", "no-cache, no-store, must-revalidate")
        self.send_header("Pragma", "no-cache")
        self.send_header("Expires", "0")
        super().end_headers()

# Allow address reuse so you don't get "Address already in use" errors immediately after restart
socketserver.TCPServer.allow_reuse_address = True

with socketserver.TCPServer(("", PORT), NoCacheHTTPRequestHandler) as httpd:
    print(f"🐍 Serving HTTP on 0.0.0.0:{PORT} (Cache disabled)")
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\nStopping server...")
        sys.exit(0)
EOF