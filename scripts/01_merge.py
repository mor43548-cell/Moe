"""Step 1: merge the uploaded column files into ONE table and export it.

No value is altered. Column order = upload order (Y, T, W, E, R, Q).
Row index = original line order in the files (1..N).
"""
import hashlib, json, pathlib
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
OUT = ROOT / "data" / "merged"
COLS = ["Y", "T", "W", "E", "R", "Q"]

frames = []
for c in COLS:
    df = pd.read_csv(RAW / f"{c}.csv")
    assert list(df.columns) == [c], (c, df.columns)
    frames.append(df)
lens = {c: len(f) for c, f in zip(COLS, frames)}
assert len(set(lens.values())) == 1, lens
table = pd.concat(frames, axis=1)
table.index = range(1, len(table) + 1)
table.index.name = "row"

OUT.mkdir(parents=True, exist_ok=True)
table.to_csv(OUT / "merged_table.csv")
table.to_excel(OUT / "merged_table.xlsx")

v = table[COLS].to_numpy()
asc = table[["Q", "W", "E", "R", "T", "Y"]].to_numpy()
checks = {
    "files_received": len(COLS),
    "files_announced_by_user": 7,
    "rows": int(len(table)),
    "rows_per_file": lens,
    "missing_values": int(table.isna().sum().sum()),
    "all_integer": bool((table.dtypes == "int64").all()),
    "global_min": int(v.min()),
    "global_max": int(v.max()),
    "per_column_min": {c: int(table[c].min()) for c in COLS},
    "per_column_max": {c: int(table[c].max()) for c in COLS},
    "rows_with_6_distinct_values": int(sum(len(set(r)) == 6 for r in v)),
    "rows_strictly_ascending_Q<W<E<R<T<Y": int(sum(all(r[i] < r[i + 1] for i in range(5)) for r in asc)),
    "duplicate_rows": int(table.duplicated().sum()),
    "sha256_raw": {c: hashlib.sha256((RAW / f"{c}.csv").read_bytes()).hexdigest() for c in COLS},
}
(OUT / "merge_checks.json").write_text(json.dumps(checks, indent=2))
print(json.dumps(checks, indent=2))
print(table.head(10).to_string())
