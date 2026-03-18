#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: ./zenodo_download.sh [--all_models]

Downloads data/models from Zenodo record 18312238.

Default behavior (no flags):
  - If ./dataset exists: warn and skip dataset download/extraction
  - If ./models exists:  warn and skip models download/extraction
  - Otherwise:
      * download dataset.zip -> extract to ./dataset
      * download two model files into ./models:
          - MAD-SURF.model
          - MAD-SURF_fewshot.model

With --all_models:
  - download models.zip -> extract to ./models (unless ./models already exists)
EOF
}

ALL_MODELS=0

for arg in "$@"; do
  case "$arg" in
    --all_models)
      ALL_MODELS=1
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Error: unknown argument '$arg'" >&2
      usage >&2
      exit 2
      ;;
  esac
done


download_file() {
  local url="$1"
  local output_file="$2"

  echo "Downloading $output_file"
  if command -v curl >/dev/null 2>&1; then
    curl -fL \
      --retry 10 \
      --retry-delay 10 \
      --retry-connrefused \
      --continue-at - \
      --progress-bar \
      "$url" -o "$output_file"
  else
    wget -c --tries=10 --waitretry=10 --progress=bar:force:noscroll "$url" -O "$output_file"
  fi
}

extract_zip() {
  local zip_file="$1"
  local output_dir="$2"

  if ! command -v unzip >/dev/null 2>&1; then
    echo "Error: unzip not found. Install it (e.g., sudo apt-get install unzip)." >&2
    exit 1
  fi

  echo "Extracting $(basename "$zip_file") into $output_dir"
  unzip -q "$zip_file" -d "$output_dir"
  echo "Done."
}

warn_skip_if_dir_exists() {
  local dir="$1"
  local label="$2"
  if [[ -d "$dir" ]]; then
    echo "Warning: '$dir' already exists; skipping ${label} download/extraction." >&2
    return 0
  fi
  return 1
}

# Use /records/ as primary; fall back to /record/
ZENODO_RECORD_ID="18312238"
ZENODO_BASE_URL_PRIMARY="https://zenodo.org/records/${ZENODO_RECORD_ID}/files"
ZENODO_BASE_URL_FALLBACK="https://zenodo.org/record/${ZENODO_RECORD_ID}/files"

fetch_from_zenodo() {
  local filename="$1"
  local outpath="$2"

  local url_primary="${ZENODO_BASE_URL_PRIMARY}/${filename}?download=1"
  local url_fallback="${ZENODO_BASE_URL_FALLBACK}/${filename}?download=1"

  if ! download_file "$url_primary" "$outpath"; then
    echo "Primary URL failed for ${filename}; trying fallback..." >&2
    download_file "$url_fallback" "$outpath"
  fi
}

REPO_DIR="$(pwd)"
DATASET_DIR="${REPO_DIR}/dataset"
MODELS_DIR="${REPO_DIR}/models"

# --- Dataset ---
if ! warn_skip_if_dir_exists "$DATASET_DIR" "dataset"; then
  tmp_zip="${REPO_DIR}/dataset.zip"
  fetch_from_zenodo "dataset.zip" "$tmp_zip"
  extract_zip "$tmp_zip" "$REPO_DIR"
  rm -f "$tmp_zip"
fi

# --- Models ---
if warn_skip_if_dir_exists "$MODELS_DIR" "models"; then
  echo "Done."
  exit 0
fi

mkdir -p "$MODELS_DIR"

if [[ "$ALL_MODELS" -eq 1 ]]; then
  tmp_zip="${REPO_DIR}/models.zip"
  fetch_from_zenodo "models.zip" "$tmp_zip"
  extract_zip "$tmp_zip" "$REPO_DIR"
  rm -f "$tmp_zip"
else
  # Download only the two core model files into ./models
  fetch_from_zenodo "MAD-SURF.model" "${MODELS_DIR}/MAD-SURF.model"
  fetch_from_zenodo "MAD-SURF_fewshot.model" "${MODELS_DIR}/MAD-SURF_fewshot.model"
fi

echo "Downloads finished successfully."

