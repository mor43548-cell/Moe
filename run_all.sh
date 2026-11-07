#!/usr/bin/env bash
# Reproduce the whole study from the 7 raw files. Requires: python3 (numpy, pandas, scipy, matplotlib, openpyxl), node + playwright (PDF only).
set -euo pipefail
cd "$(dirname "$0")"
python3 scripts/01_merge.py
python3 scripts/validate_nist.py          # optional: pass a directory containing the 'nistrng' package for a cross-check
python3 scripts/02_analysis.py
python3 scripts/03_report.py
node scripts/render_pdf.js || echo "PDF rendering skipped (node/playwright not available)"
