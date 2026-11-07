# Randomness evaluation of the Y/T/W/E/R/Q table

**Result:** the success criterion was a 100% pass rate (15/15 in NIST SP 800-22). The data **does not meet it**.
NIST: 2 PASS, 10 FAIL, 3 NOT APPLICABLE. Full report: `report/randomness_report.pdf` (Hebrew), with the same content in `.html` and `.xlsx`.

## Inputs
- 6 files were received: Y, T, W, E, R, Q (`data/raw/`). The request mentioned 7; only these 6 were analysed.
- Merged table (exported before any analysis): `data/merged/merged_table.csv` / `.xlsx` (328 rows × 6), with checks in `merge_checks.json`.
- Structure: every row holds 6 distinct numbers from 1..37 with Q<W<E<R<T<Y in every row, i.e. a sorted 6-of-37 set.

## Pipeline (run in order)
| script | what it does | output |
|---|---|---|
| `scripts/01_merge.py` | merges the files into one table, integrity checks, SHA-256 | `data/merged/` |
| `scripts/02_validate_sp800_22.py` | validates the SP 800-22 implementation against the spec's worked examples and the reference p-values for e (10^6 bits) | `output/machine_output/02_validation*` |
| `scripts/03_analysis.py` | χ²/KS, correlogram, runs tests, Shannon entropy, 15 NIST tests | `output/machine_output/03_*`, `output/figures/` |
| `scripts/04_report.py` (+ `scripts/tools/html_to_pdf.js`) | scorecard, checklist, HTML/PDF/XLSX report | `report/` |

```
pip install -r requirements.txt
python scripts/01_merge.py && python scripts/02_validate_sp800_22.py && python scripts/03_analysis.py && python scripts/04_report.py
node scripts/tools/html_to_pdf.js   # optional, needs playwright
```

## Notes for an external reviewer
- The official NIST STS 2.1.2 could not be downloaded (csrc.nist.gov was blocked in the build environment).
  `scripts/sp800_22.py` implements the 15 tests from the spec with STS constants. Validation: 59/61 reference values match.
  The 2 mismatches are the two short DFT examples printed in the spec. The full-length e DFT reference does match.
- To re-run with the official STS: `output/machine_output/bitstream.txt` is the exact 11,808-bit ASCII 0/1 input
  (each number as 6-bit binary of its value, in row-major order Y,T,W,E,R,Q).
- Under this encoding, even a perfectly uniform 1..37 source gives P(bit=1)=0.4189, not 0.5.
  The observed fraction is 0.4204, so the bit-level NIST failures are expected from the encoding alone.
