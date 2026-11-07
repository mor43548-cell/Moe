# Randomness assessment of 7 variables (Y, T, W, R, ABC, Q, E)

The package an external reviewer needs: input data, code, unedited machine output and the structured report.

## Result against the pre-defined criterion (pass 100% of the randomness tests): NOT MET

| Block | Result |
|---|---|
| Uniformity, chi-square (p > 0.01) | 1/7 variables pass (ABC only). Proportion 0.1429 against a target of ~0.99 |
| Uniformity, discrete KS (p > 0.01) | 1/7 (ABC only) |
| Runs test, original values (p > 0.05) | up/down: 6/7 pass (Q fails, p = 0.0086). Above/below median: 6/7 pass (R fails, p = 0.0490) |
| Autocorrelation, lags 1-20 | max \|r\| ≤ 0.161. 10/140 coefficients fall outside the 95% band (about 7 expected by chance). All bars are inside the band only for W and ABC |
| Shannon entropy ≥ 99.9% of maximum | 0/7 variables (best: ABC 0.9938) |
| NIST SP 800-22 (15 tests) | PASS 0/15, FAIL 9/15, N/A 6/15 (sequence shorter than the standard's minimum) |

The full checklist is in `report/Randomness_Report.pdf` (sections, scorecard, figures, verbal summary, appendices) and in
`report/Randomness_Report.xlsx` (every table plus the raw machine output).

## Layout
```
data/raw/            the 7 input files, unchanged
data/merged/         merged table (CSV + XLSX), step 1
scripts/01_merge.py  step 1: merge
scripts/nist_sp800_22.py  Python implementation of the 15 NIST SP 800-22 Rev.1a tests
scripts/validate_nist.py  checks the implementation against the worked examples in the NIST document
scripts/02_analysis.py    uniformity, independence, entropy, NIST
scripts/03_report.py + render_pdf.js   report (Excel, HTML, PDF)
results/             machine output (CSV/TXT/JSON). machine_output.txt is the full program log
results/bitstreams/  the exact bit sequences fed to NIST, so another tool (for example the official STS) can re-test them
figures/             frequency histograms, correlograms
run_all.sh           reproduces everything
```

## Input integrity (SHA-256)
```
14a52efa0a5313d59be44d5bca23cc7c19e06364b2ca8107ccc7d04e1f4e46ea  Y.csv
88725a1abd9ea971f70bbbf3e89ae0ecde3153c65f639fd0487da631443fdbf9  T.csv
2d141bafebd0bb30c74a13049ef1218819cf3bd4820e34cdabf83ece32e78ca1  W.csv
705506bfd70059e3bf88ddbf59efe4f8299361879c330af21a1dccbfa77cb976  R.csv
c8d7ca17846d48eaec1ee2c464d942678e8de10a1b65a3eeb48ec2760b3e526b  ABC.csv
91adb601c658f0e742d675b159dbe80c18c7fbca16b40bf789f03ff0c456f78b  Q.csv
7d7a1d823c1a5b009221d4bfc380f30b921e8e87455ea0fc7161b2866e027515  E.csv
```

## Disclosures for the reviewer
1. **NIST tool.** The official NIST STS (C) could not be downloaded in the analysis environment because the host was blocked.
   The suite was implemented in Python from SP 800-22 Rev.1a. Validation: 12 of 13 worked examples in the document are
   reproduced to 6 decimals. The DFT example is **not** reproduced, so the DFT result is flagged (`results/validation_nist.txt`).
   The bit streams are saved in `results/bitstreams/` so they can be re-run in the official STS.
2. **Tests not runnable.** Binary Matrix Rank, Overlapping Template, Maurer, Linear Complexity, Random Excursions and
   Random Excursions Variant need between 38,912 and 1,000,000 bits (and ≥ 500 cycles). The data give 10,824 bits and
   10 cycles. They are reported as N/A, not as PASS.
3. **Uniformity domain** = the integers from the observed minimum to the observed maximum of each variable. The true theoretical range is unknown.
4. **Binary encoding** = value − min in ⌈log₂k⌉ bits. When k is not a power of 2, this encoding lowers the share of 1-bits even for perfectly uniform data.
5. **Correction made during the work.** The first run used `scipy.stats.kstest`, which is wrong for tied discrete data
   because it inflates D. It was replaced by the discrete KS (Conover 1972). Only ABC changed: p went from ≈7.6e-10 to 0.719.
6. **Entropy at N=328.** The plug-in estimator is biased downward. Even exactly uniform data would give a ratio of about 0.987-0.995,
   so the 99.9% threshold is almost unreachable at this sample size. The threshold was not changed.
7. Passing a test means only that H0 is not rejected. It does not prove randomness.

## Reproduce
```
pip install -r requirements.txt   # tested: Python 3.13, numpy 2.5.3, scipy 1.18.1, pandas 3.0.5, matplotlib 3.11.2, openpyxl 3.1.5
./run_all.sh
```
