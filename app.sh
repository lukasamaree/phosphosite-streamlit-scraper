#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

PYTHON="${PYTHON:-.venv312/bin/python}"
PORT="${PORT:-8523}"

exec "$PYTHON" -m streamlit run streamlit_phospho_scraper.py \
  --server.port "$PORT" \
  --server.address 127.0.0.1 \
  --server.headless true \
  --browser.gatherUsageStats false \
  --server.fileWatcherType none
