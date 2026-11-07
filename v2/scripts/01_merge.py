"""V2 step 1: one merged table of the 7 files, columns labelled under the lotto hypothesis
(6 main numbers 1..37 sorted ascending = N1..N6, strong number 1..7). Original names kept in the header."""
import pandas as pd
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
MAP = [("Q", "N1_smallest"), ("W", "N2"), ("E", "N3"), ("R", "N4"), ("T", "N5"), ("Y", "N6_largest"), ("ABC", "Strong")]
cols = {}
for orig, new in MAP:
    s = pd.read_csv(ROOT / "data/raw" / f"{orig}.csv")[orig]
    cols[f"{new} ({orig})"] = s
df = pd.DataFrame(cols); df.insert(0, "draw", range(1, len(df) + 1))
main = df.iloc[:, 1:7].to_numpy()
assert ((main[:, 1:] > main[:, :-1]).all()), "rows not strictly increasing"
assert main.min() >= 1 and main.max() <= 37 and df.iloc[:, 7].between(1, 7).all()
df.to_csv(ROOT / "v2/data/merged_lotto_table.csv", index=False)
df.to_excel(ROOT / "v2/data/merged_lotto_table.xlsx", index=False)
print(df.shape); print(df.head().to_string(index=False))
print("structure check: every row has 6 distinct main numbers in 1..37 (strictly increasing) and Strong in 1..7: OK")
