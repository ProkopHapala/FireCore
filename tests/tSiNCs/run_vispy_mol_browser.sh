#!/bin/bash
# Launch Vispy molecular browser on NPZ viewer fixtures.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
export PYTHONPATH="${ROOT}:${PYTHONPATH:-}"
FIX="${ROOT}/tests/tSiNCs/fixtures/npz_viewer"
if [[ ! -f "${FIX}/01_init.npz" ]]; then
  python3 "${FIX}/bootstrap_fixtures.py"
fi
DIR="${1:-${FIX}}"
exec python3 -m pyBall.GUI.VispyMolBrowser --dir "${DIR}"
