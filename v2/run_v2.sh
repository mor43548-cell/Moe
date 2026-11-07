#!/usr/bin/env bash
# Reproduce version 2 (lotto hypothesis). Design fixed in v2/DESIGN.md before testing.
set -euo pipefail
cd "$(dirname "$0")/.."
python3 v2/scripts/01_merge.py
python3 scripts/validate_nist.py
python3 v2/scripts/02_analysis.py
python3 v2/scripts/03_report.py
node v2/scripts/render_pdf.js || echo "PDF rendering skipped"
