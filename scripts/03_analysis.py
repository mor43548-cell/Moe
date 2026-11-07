"""Step 3: run exactly the requested measures/tests on the merged table.

Objects analysed (all derived from data/merged/merged_table.csv, no value altered):
  * each column Y,T,W,E,R,Q as a sequence of 328 values (draw order = file row order)
  * POOLED: all 6x328 = 1968 numbers (every number of every row)
  * STREAM: row-major sequence Y,T,W,E,R,Q,Y,T,... (1968 numbers) for correlogram/runs
  * BITS: STREAM with each number written as its own 6-bit binary (37 < 64): 11,808 bits,
          input to NIST SP 800-22 (which is defined on bits only).

Thresholds are the ones given by the requester; they are applied, never modified.
"""
import json, math, pathlib, sys
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import sp800_22 as NIST

ROOT = pathlib.Path(__file__).resolve().parents[1]
FIG = ROOT / "output" / "figures"; FIG.mkdir(parents=True, exist_ok=True)
MO = ROOT / "output" / "machine_output"; MO.mkdir(parents=True, exist_ok=True)
COLS = ["Y", "T", "W", "E", "R", "Q"]
K = 37                      # alphabet 1..37 (observed global min 1, max 37)
LAGS = 20
plt.rcParams.update({"figure.dpi": 150, "font.size": 9})

tbl = pd.read_csv(ROOT / "data" / "merged" / "merged_table.csv", index_col="row")
n_rows = len(tbl)
stream = tbl[COLS].to_numpy().ravel()          # row-major
pooled = stream
R = {"n_rows": n_rows, "columns": COLS, "alphabet": [1, K]}

# ---------------------------------------------------------------- 1. Uniformity
def chi2_uniform(x):
    obs = np.bincount(x, minlength=K + 1)[1:]
    exp = np.full(K, len(x) / K)
    chi2, p = stats.chisquare(obs, exp)
    return dict(n=int(len(x)), chi2=float(chi2), df=K - 1, p_value=float(p),
                expected_per_value=float(exp[0]), observed=obs.tolist(),
                min_obs=int(obs.min()), max_obs=int(obs.max()))

def ks_uniform(x):
    # KS against the discrete uniform CDF on 1..37.  For discrete data the statistic is
    # D = max_v |ECDF(v) - F(v)| over the support points v = 1..37 (scipy.stats.kstest assumes a
    # continuous F and, with ties, inflates D by the jump size, so it is not used here).
    # p-value: Kolmogorov distribution for sample size n (exact, scipy.stats.kstwo); for a
    # discrete F this p-value is conservative (Conover 1972).
    n = len(x)
    ecdf = np.cumsum(np.bincount(x, minlength=K + 1)[1:]) / n
    F = np.arange(1, K + 1) / K
    D = float(np.max(np.abs(ecdf - F)))
    return dict(D=D, p_value=float(stats.kstwo.sf(D, n)), n=n,
                method="discrete KS: D over support points; Kolmogorov p (conservative for discrete F)")

U = {"pooled": {"chi2": chi2_uniform(pooled), "ks": ks_uniform(pooled)}}
for c in COLS:
    U[c] = {"chi2": chi2_uniform(tbl[c].to_numpy()), "ks": ks_uniform(tbl[c].to_numpy())}
R["uniformity"] = U

# histogram: pooled
obs = np.array(U["pooled"]["chi2"]["observed"]); e = U["pooled"]["chi2"]["expected_per_value"]
fig, ax = plt.subplots(figsize=(9, 3.6))
ax.bar(range(1, K + 1), obs, color="#4C72B0", edgecolor="white")
ax.axhline(e, color="black", lw=1.2, label=f"expected under uniformity = {e:.2f}")
band = 1.96 * math.sqrt(n_rows * (6 / K) * (1 - 6 / K))  # count of a number ~ Binomial(328, 6/37)
ax.axhspan(e - band, e + band, color="grey", alpha=0.2, label="±1.96·sd, Binomial(328, 6/37)")
ax.set_xticks(range(1, K + 1)); ax.set_xlabel("number"); ax.set_ylabel("count")
ax.set_title(f"Frequency histogram, all {len(pooled)} numbers (6 x {n_rows})   "
             f"chi2={U['pooled']['chi2']['chi2']:.2f}, df=36, p={U['pooled']['chi2']['p_value']:.4f}")
ax.legend(loc="upper right", fontsize=8); fig.tight_layout(); fig.savefig(FIG / "fig1_histogram_pooled.png"); plt.close(fig)

# histogram: per column
fig, axes = plt.subplots(2, 3, figsize=(11, 5.5), sharex=True)
for ax, c in zip(axes.ravel(), COLS):
    o = np.array(U[c]["chi2"]["observed"]); ec = U[c]["chi2"]["expected_per_value"]
    ax.bar(range(1, K + 1), o, color="#4C72B0", edgecolor="white")
    ax.axhline(ec, color="black", lw=1)
    ax.set_title(f"{c}: chi2 p={U[c]['chi2']['p_value']:.2e}", fontsize=9)
    ax.set_xticks([1, 10, 20, 30, 37])
fig.suptitle("Frequency histogram per column (line = uniform expectation 328/37)")
fig.tight_layout(); fig.savefig(FIG / "fig2_histogram_columns.png"); plt.close(fig)

# ---------------------------------------------------------------- 2. Independence
def acf(x, nlags):
    x = np.asarray(x, float); x = x - x.mean(); d = np.sum(x * x)
    return np.array([1.0] + [np.sum(x[k:] * x[:-k]) / d for k in range(1, nlags + 1)])

def runs_updown(x):
    d = np.diff(np.asarray(x, float)); d = d[d != 0]          # ties removed (standard handling)
    s = np.sign(d); m = len(s) + 1                             # effective length
    runs = 1 + int(np.sum(s[1:] != s[:-1]))
    mu = (2 * m - 1) / 3; var = (16 * m - 29) / 90
    z = (runs - mu) / math.sqrt(var)
    return dict(runs=runs, expected=mu, z=z, p_value=float(2 * stats.norm.sf(abs(z))),
                ties_removed=int(len(x) - 1 - len(d)))

def runs_median(x):
    x = np.asarray(x, float); med = np.median(x); b = x[x != med] > med
    n1, n2 = int(b.sum()), int((~b).sum()); n = n1 + n2
    runs = 1 + int(np.sum(b[1:] != b[:-1]))
    mu = 2 * n1 * n2 / n + 1; var = 2 * n1 * n2 * (2 * n1 * n2 - n) / (n * n * (n - 1))
    z = (runs - mu) / math.sqrt(var)
    return dict(median=float(med), n_above=n1, n_below=n2, runs=runs, expected=mu, z=z,
                p_value=float(2 * stats.norm.sf(abs(z))), dropped_equal_to_median=int(len(x) - n))

def runs_identical(x):
    x = np.asarray(x); n = len(x)
    rep = int(np.sum(x[1:] == x[:-1]))
    p = np.bincount(x)[1:] / n; q = float(np.sum(p ** 2))      # P(two independent draws equal)
    expct = (n - 1) * q
    # two-sided exact binomial on the number of consecutive identical pairs
    pv = float(stats.binomtest(rep, n - 1, q).pvalue)
    return dict(consecutive_identical_pairs=rep, expected=expct, p_value=pv)

ind = {}
band_col = 1.96 / math.sqrt(n_rows)
for c in COLS:
    x = tbl[c].to_numpy()
    a = acf(x, LAGS)
    ind[c] = dict(acf=a.tolist(), band95=band_col,
                  lags_outside_band=[k for k in range(1, LAGS + 1) if abs(a[k]) > band_col],
                  max_abs_acf=float(np.max(np.abs(a[1:]))),
                  runs_updown=runs_updown(x), runs_median=runs_median(x), runs_identical=runs_identical(x))
a = acf(stream, 30); band_s = 1.96 / math.sqrt(len(stream))
ind["stream"] = dict(acf=a.tolist(), band95=band_s,
                     lags_outside_band=[k for k in range(1, 31) if abs(a[k]) > band_s],
                     max_abs_acf=float(np.max(np.abs(a[1:]))),
                     runs_updown=runs_updown(stream), runs_median=runs_median(stream))
R["independence"] = ind

fig, axes = plt.subplots(2, 3, figsize=(11, 5.5), sharey=True)
for ax, c in zip(axes.ravel(), COLS):
    a = np.array(ind[c]["acf"])
    ax.vlines(range(LAGS + 1), 0, a, color="#4C72B0", lw=2); ax.plot(range(LAGS + 1), a, "o", ms=3, color="#4C72B0")
    ax.axhspan(-band_col, band_col, color="grey", alpha=0.25); ax.axhline(0, color="black", lw=.6)
    ax.set_title(f"{c}: {len(ind[c]['lags_outside_band'])}/{LAGS} lags outside ±{band_col:.3f}", fontsize=9)
    ax.set_xlabel("lag (draws)"); ax.set_ylim(-0.3, 1.05)
fig.suptitle("Correlogram per column (95% band = ±1.96/√328)")
fig.tight_layout(); fig.savefig(FIG / "fig3_correlogram_columns.png"); plt.close(fig)

fig, ax = plt.subplots(figsize=(9, 3.4))
a = np.array(ind["stream"]["acf"])
ax.vlines(range(31), 0, a, color="#4C72B0", lw=2); ax.plot(range(31), a, "o", ms=3, color="#4C72B0")
ax.axhspan(-band_s, band_s, color="grey", alpha=0.25); ax.axhline(0, color="black", lw=.6)
ax.set_xlabel("lag (positions in row-major stream Y,T,W,E,R,Q,...)")
ax.set_title(f"Correlogram of the row-major stream (n={len(stream)}, band ±{band_s:.3f}): "
             f"{len(ind['stream']['lags_outside_band'])}/30 lags outside")
fig.tight_layout(); fig.savefig(FIG / "fig4_correlogram_stream.png"); plt.close(fig)

# ---------------------------------------------------------------- 3. Entropy
def shannon(x, k):
    p = np.bincount(x)[np.bincount(x) > 0] / len(x)
    h = float(-np.sum(p * np.log2(p)))
    return dict(H_bits=h, H_max=math.log2(k), ratio=h / math.log2(k))

bits = np.array([int(b) for v in stream for b in format(int(v), "06b")], dtype=np.int8)
ent = {"pooled": shannon(pooled, K)}
for c in COLS:
    ent[c] = shannon(tbl[c].to_numpy(), K)
ent["bits"] = shannon(bits, 2)
R["entropy"] = ent

# ---------------------------------------------------------------- 4. NIST SP 800-22 on BITS
popc = [bin(v).count("1") for v in range(1, K + 1)]
R["bit_encoding"] = dict(
    rule="each number written as 6-bit unsigned binary of its own value (1 -> 000001, 37 -> 100101), row-major order",
    n_bits=int(len(bits)), ones=int(bits.sum()), ones_fraction=float(bits.mean()),
    ones_fraction_if_numbers_perfectly_uniform_1_37=float(np.mean(popc) / 6),
    note="A perfectly uniform 1..37 source does NOT give P(bit=1)=1/2 under fixed 6-bit encoding; "
         "the expected one-fraction is the value above.")
nist = []
for f in NIST.ALL:
    r = f(bits)
    ps = r["p_values"]
    r["n_pvalues"] = len(ps)
    r["n_pvalues_ge_0.01"] = int(sum(p >= 0.01 for p in ps))
    r["n_pvalues_gt_0.05"] = int(sum(p > 0.05 for p in ps))
    r["min_p"] = float(min(ps)) if ps else None
    r["sequences_m"] = 1
    r["proportion"] = (1.0 if all(p >= 0.01 for p in ps) else 0.0) if ps else None
    r["status"] = ("NOT APPLICABLE" if not r["applicable"] else
                   "PASS" if all(p >= 0.01 for p in ps) else "FAIL")
    r["n_meets_recommended_min"] = (None if r["min_n_recommended"] is None else len(bits) >= r["min_n_recommended"])
    nist.append(r)
R["nist"] = nist
R["nist_proportion_interval_m1"] = [0.99 - 3 * math.sqrt(0.99 * 0.01 / 1), 0.99 + 3 * math.sqrt(0.99 * 0.01 / 1)]

(MO / "03_results.json").write_text(json.dumps(R, indent=2))
np.savetxt(MO / "bitstream.txt", bits[None, :], fmt="%d", delimiter="")

# ---------------------------------------------------------------- machine-readable stdout
L = []
P = L.append
P("=== 1. UNIFORMITY ===")
for k in ["pooled"] + COLS:
    u = U[k]
    P(f"{k:7s} chi2={u['chi2']['chi2']:.4f} df=36 p={u['chi2']['p_value']:.6g} | KS D={u['ks']['D']:.4f} p={u['ks']['p_value']:.6g}")
P("=== 2. INDEPENDENCE ===")
for k in COLS + ["stream"]:
    d = ind[k]
    s = f"{k:7s} max|acf|={d['max_abs_acf']:.4f} lags_outside={d['lags_outside_band']} | updown runs={d['runs_updown']['runs']} exp={d['runs_updown']['expected']:.2f} p={d['runs_updown']['p_value']:.6g} | median runs={d['runs_median']['runs']} exp={d['runs_median']['expected']:.2f} p={d['runs_median']['p_value']:.6g}"
    if "runs_identical" in d:
        s += f" | identical pairs={d['runs_identical']['consecutive_identical_pairs']} exp={d['runs_identical']['expected']:.2f} p={d['runs_identical']['p_value']:.6g}"
    P(s)
P("=== 3. ENTROPY ===")
for k, v in ent.items():
    P(f"{k:7s} H={v['H_bits']:.6f} Hmax={v['H_max']:.6f} ratio={v['ratio']:.6f}")
P("=== 4. NIST SP 800-22 (n=%d bits, m=1 sequence) ===" % len(bits))
P(json.dumps(R["bit_encoding"]))
for r in nist:
    P(f"{r['name']:40s} status={r['status']:15s} params={r['params']} p={['%.6f' % p for p in r['p_values']][:20]}{' ...' if len(r['p_values'])>20 else ''} "
      f"#p>=0.01={r['n_pvalues_ge_0.01']}/{r['n_pvalues']} note={r['note']}")
out = "\n".join(L)
(MO / "03_analysis_stdout.txt").write_text(out)
print(out)
