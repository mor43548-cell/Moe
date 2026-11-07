"""Step 1: merge the 7 raw files into one table (row-aligned by position) and export."""
import pandas as pd, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
COLS = ["Y", "T", "W", "R", "ABC", "Q", "E"]
frames = []
for c in COLS:
    s = pd.read_csv(ROOT / "data/raw" / f"{c}.csv")
    assert list(s.columns) == [c], (c, s.columns)
    frames.append(s[c])
df = pd.concat(frames, axis=1)
assert df.notna().all().all(), "missing values"
assert len({len(f) for f in frames}) == 1, "length mismatch"
df.insert(0, "row", range(1, len(df) + 1))
df.to_csv(ROOT / "data/merged/merged_7_files.csv", index=False)
df.to_excel(ROOT / "data/merged/merged_7_files.xlsx", index=False)
print(df.shape); print(df.describe().T.to_string())
