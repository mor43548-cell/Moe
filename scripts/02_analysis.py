"""
Step 2: the analyses requested, in the order requested.
  1. Uniformity     - chi-square frequency test + KS test per column, proportion, frequency histograms
  2. Independence   - autocorrelation (correlogram), runs test (up/down and above/below median)
  3. Entropy        - Shannon entropy relative to the maximum
  4. NIST SP 800-22 - 15 tests on the binary representation of the data

All fixed choices (domain, encoding, lags, alpha) are declared at the top and were fixed before
any result was computed. Every number written to results/ is the unmodified program output.
"""
import json, math, sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
import nist_sp800_22 as N

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"; FIG = ROOT / "figures"
RES.mkdir(exist_ok=True); FIG.mkdir(exist_ok=True)

# ---------------- fixed design choices ----------------
COLS = ["Y", "T", "W", "R", "ABC", "Q", "E"]
ALPHA_NIST = 0.01         # NIST significance level ("p > 0.01")
ALPHA_STRICT = 0.05       # "aspire to p > 0.05"; mandatory threshold for the runs test as specified
MAX_LAG = 20              # autocorrelation lags 1..20
CONF_Z = 1.959964         # 95% confidence band for the correlogram
TARGET_PROP = 0.99        # NIST target proportion
ENTROPY_TARGET = 0.999    # >= 99.9% of maximum entropy
# Domain of each variable for the uniformity hypothesis: integers from observed min to observed max.
# Binary encoding for NIST: value - min, written in ceil(log2(k)) bits, most significant bit first.
# Combined stream: row-major interleave of the 7 columns in the order above.

log = []
def out(s=""):
    print(s); log.append(str(s))

df = pd.read_csv(ROOT / "data/merged/merged_7_files.csv")
n_obs = len(df)
out(f"Merged table: {n_obs} rows x {len(COLS)} columns {COLS}")

def nist_prop_bounds(m, p=TARGET_PROP):
    half = 3 * math.sqrt(p * (1 - p) / m)
    return p - half, min(1.0, p + half)

# =========================================================
# 1. UNIFORMITY
# =========================================================
out("\n=== 1. UNIFORMITY ===")
uni_rows, counts_by_col = [], {}
for c in COLS:
    x = df[c].to_numpy()
    lo, hi = int(x.min()), int(x.max()); k = hi - lo + 1
    cats = np.arange(lo, hi + 1)
    obs = np.array([(x == v).sum() for v in cats])
    exp = np.full(k, n_obs / k)
    chi2, p_chi = stats.chisquare(obs, exp)
    # Discrete KS (Conover 1972): D = max over the support points of |F_n(v) - F_0(v)|.
    # scipy.stats.kstest is NOT used: its continuous formula counts the jump at tied values into D.
    # p-value from the continuous Kolmogorov distribution, which is conservative for a discrete null.
    ks_D = float(np.max(np.abs(np.cumsum(obs) / n_obs - np.arange(1, k + 1) / k)))
    ks = type("KS", (), dict(statistic=ks_D, pvalue=float(stats.kstwo.sf(ks_D, n_obs))))
    counts_by_col[c] = (cats, obs, n_obs / k)
    uni_rows.append(dict(column=c, min=lo, max=hi, k=k, N=n_obs, expected_per_value=n_obs / k,
                         min_observed=int(obs.min()), max_observed=int(obs.max()), empty_values=int((obs == 0).sum()),
                         chi2=chi2, df=k - 1, p_chi2=p_chi, KS_D=ks.statistic, p_KS=ks.pvalue,
                         chi2_p_gt_0_01=p_chi > ALPHA_NIST, chi2_p_gt_0_05=p_chi > ALPHA_STRICT,
                         KS_p_gt_0_01=ks.pvalue > ALPHA_NIST, KS_p_gt_0_05=ks.pvalue > ALPHA_STRICT))
uni = pd.DataFrame(uni_rows)
lo_b, hi_b = nist_prop_bounds(len(COLS))
prop_chi = float((uni.p_chi2 >= ALPHA_NIST).mean()); prop_ks = float((uni.p_KS >= ALPHA_NIST).mean())
out(uni.to_string(index=False))
out(f"Proportion of columns passing chi2 at alpha=0.01: {prop_chi:.4f}  (NIST acceptable range for m=7: [{lo_b:.4f}, {hi_b:.4f}])")
out(f"Proportion of columns passing KS   at alpha=0.01: {prop_ks:.4f}")
uni.to_csv(RES / "1_uniformity.csv", index=False)

fig, axes = plt.subplots(2, 4, figsize=(18, 8))
for ax, c in zip(axes.flat, COLS):
    cats, obs, e = counts_by_col[c]
    ax.bar(cats, obs, color="#4C72B0", edgecolor="white", width=0.85, label="observed")
    ax.axhline(e, color="#C44E52", lw=2, ls="--", label=f"uniform expectation = {e:.1f}")
    r = uni.set_index("column").loc[c]
    ax.set_title(f"{c}   chi2 p = {r.p_chi2:.3g}   KS p = {r.p_KS:.3g}", fontsize=11)
    ax.set_xlabel("value"); ax.set_ylabel("count"); ax.legend(fontsize=8)
axes.flat[-1].axis("off")
fig.suptitle("Frequency histograms (observed counts vs. flat uniform expectation), N = 328 per column", fontsize=14)
fig.tight_layout(); fig.savefig(FIG / "1_frequency_histograms.png", dpi=150); plt.close(fig)

# =========================================================
# 2. INDEPENDENCE
# =========================================================
out("\n=== 2. INDEPENDENCE ===")
def acf(x, L):
    x = np.asarray(x, float); x = x - x.mean(); d = np.sum(x * x)
    return np.array([1.0] + [np.sum(x[:-k] * x[k:]) / d for k in range(1, L + 1)])

def runs_up_down(x):
    d = np.diff(np.asarray(x, float)); s = np.sign(d); s = s[s != 0]
    n = len(s) + 1
    R = 1 + int(np.sum(s[1:] != s[:-1]))
    mu = (2 * n - 1) / 3.0; var = (16 * n - 29) / 90.0
    z = (R - mu) / math.sqrt(var)
    return R, mu, var, z, 2 * stats.norm.sf(abs(z)), int(np.sum(np.diff(x) == 0))

def runs_median(x):
    x = np.asarray(x, float); med = np.median(x)
    s = x[x != med] > med
    n1, n2 = int(s.sum()), int((~s).sum()); n = n1 + n2
    R = 1 + int(np.sum(s[1:] != s[:-1]))
    mu = 2 * n1 * n2 / n + 1; var = 2 * n1 * n2 * (2 * n1 * n2 - n) / (n * n * (n - 1))
    z = (R - mu) / math.sqrt(var)
    return R, mu, var, z, 2 * stats.norm.sf(abs(z)), med, n1, n2, int(np.sum(x == med))

band = CONF_Z / math.sqrt(n_obs)
acf_tab, ind_rows = {}, []
for c in COLS:
    x = df[c].to_numpy()
    r = acf(x, MAX_LAG); acf_tab[c] = r
    out_band = [k for k in range(1, MAX_LAG + 1) if abs(r[k]) > band]
    R1, mu1, v1, z1, p1, ties1 = runs_up_down(x)
    R2, mu2, v2, z2, p2, med, n1, n2, ties2 = runs_median(x)
    ind_rows.append(dict(column=c, acf_lag1=r[1], max_abs_acf_lags1_20=float(np.max(np.abs(r[1:]))),
                         band_95=band, lags_outside_band=len(out_band), which_lags_outside=" ".join(map(str, out_band)),
                         runs_updown_R=R1, runs_updown_expected=mu1, runs_updown_z=z1, runs_updown_p=p1,
                         runs_updown_zero_diffs_dropped=ties1,
                         runs_median_R=R2, runs_median_expected=mu2, runs_median_z=z2, runs_median_p=p2,
                         median=med, n_above=n1, n_below=n2, ties_at_median_dropped=ties2,
                         updown_p_gt_0_05=p1 > ALPHA_STRICT, median_p_gt_0_05=p2 > ALPHA_STRICT))
ind = pd.DataFrame(ind_rows)
acf_df = pd.DataFrame({c: acf_tab[c] for c in COLS}); acf_df.index.name = "lag"
out(ind.to_string(index=False))
out(f"Expected number of lags outside the 95% band among 20 under independence: 1.0 per column")
ind.to_csv(RES / "2_independence_runs.csv", index=False)
acf_df.to_csv(RES / "2_autocorrelation_lags_0_20.csv")
prop_ud = float((ind.runs_updown_p >= ALPHA_NIST).mean()); prop_md = float((ind.runs_median_p >= ALPHA_NIST).mean())

fig, axes = plt.subplots(2, 4, figsize=(18, 8))
for ax, c in zip(axes.flat, COLS):
    r = acf_tab[c]; lags = np.arange(MAX_LAG + 1)
    cols_ = ["#4C72B0" if (k == 0 or abs(r[k]) <= band) else "#C44E52" for k in lags]
    ax.bar(lags, r, color=cols_, width=0.6)
    ax.axhspan(-band, band, color="grey", alpha=0.2, label=f"95% band = ±{band:.3f}")
    ax.axhline(0, color="black", lw=0.8)
    ax.set_ylim(-0.4, 1.05); ax.set_title(f"{c}   lag-1 r = {r[1]:.3f}", fontsize=11)
    ax.set_xlabel("lag"); ax.set_ylabel("autocorrelation"); ax.legend(fontsize=8)
axes.flat[-1].axis("off")
fig.suptitle("Correlograms, lags 0-20 (red = outside the 95% confidence band)", fontsize=14)
fig.tight_layout(); fig.savefig(FIG / "2_correlograms.png", dpi=150); plt.close(fig)

# =========================================================
# 3. ENTROPY + binary encoding (encoding is also the NIST input)
# =========================================================
out("\n=== 3. ENTROPY ===")
streams, ent_rows = {}, []
for c in COLS:
    x = df[c].to_numpy()
    lo, hi = int(x.min()), int(x.max()); k = hi - lo + 1
    nb = int(math.ceil(math.log2(k)))
    bits = np.array([(int(v - lo) >> (nb - 1 - j)) & 1 for v in x for j in range(nb)], dtype=np.int64)
    streams[c] = bits
    p = np.bincount(x - lo, minlength=k) / n_obs; pnz = p[p > 0]
    H = float(-np.sum(pnz * np.log2(pnz))); Hmax = math.log2(k)
    p1 = bits.mean(); Hb = float(-(p1 * math.log2(p1) + (1 - p1) * math.log2(1 - p1)))
    # reference only (not a test): ones-fraction the encoding would give if values were exactly uniform on k
    ideal_p1 = sum(bin(v).count("1") for v in range(k)) / (k * nb)
    mm_ratio = 1 - (k - 1) / (2 * n_obs * math.log(2)) / Hmax  # Miller-Madow expected plug-in ratio, uniform data
    ent_rows.append(dict(column=c, k=k, H_bits=H, H_max=Hmax, H_ratio=H / Hmax, meets_99_9=H / Hmax >= ENTROPY_TARGET,
                         bits_per_value=nb, n_bits=len(bits), ones_fraction=p1, H_binary=Hb, H_binary_meets_99_9=Hb >= ENTROPY_TARGET,
                         ref_ones_fraction_if_exactly_uniform=ideal_p1, ref_expected_H_ratio_if_uniform_N328=mm_ratio))
comb = np.concatenate([np.array([(int(df[c].iloc[i] - df[c].min()) >> (len_ - 1 - j)) & 1 for c, len_ in
                                 [(cc, int(math.ceil(math.log2(df[cc].max() - df[cc].min() + 1)))) for cc in COLS]
                                 for j in range(len_)]) for i in range(n_obs)]).astype(np.int64)
streams["COMBINED"] = comb
p1 = comb.mean(); Hb = float(-(p1 * math.log2(p1) + (1 - p1) * math.log2(1 - p1)))
ent_rows.append(dict(column="COMBINED", n_bits=len(comb), ones_fraction=p1, H_binary=Hb, H_binary_meets_99_9=Hb >= ENTROPY_TARGET))
ent = pd.DataFrame(ent_rows)
out(ent.to_string(index=False))
ent.to_csv(RES / "3_entropy.csv", index=False)
for c, b in streams.items():
    (RES / "bitstreams").mkdir(exist_ok=True)
    (RES / "bitstreams" / f"{c}.txt").write_text("".join(map(str, b.tolist())) + "\n")

# =========================================================
# 4. NIST SP 800-22
# =========================================================
out("\n=== 4. NIST SP 800-22 (15 tests) ===")
rows = []
for s, bits in streams.items():
    for r in N.run_all(bits):
        if r["applicable"]:
            for i, p in enumerate(r["p_values"]):
                rows.append(dict(stream=s, n_bits=len(bits), test=r["test"], sub=i + 1, p_value=p,
                                 result="SUCCESS" if p >= ALPHA_NIST else "FAILURE",
                                 p_gt_0_05=p > ALPHA_STRICT, note=json.dumps(r["stats"], default=float)))
        else:
            rows.append(dict(stream=s, n_bits=len(bits), test=r["test"], sub=0, p_value=np.nan, result="N/A",
                             p_gt_0_05=None, note=r["reason"]))
nist = pd.DataFrame(rows)
nist.to_csv(RES / "4_nist_machine_output_all_streams.csv", index=False)
out(nist.drop(columns="note").to_string(index=False))

# all 148 aperiodic templates (m=9) for the non-overlapping template test (as the reference suite does)
tmpl = N.aperiodic_templates(9)
assert len(tmpl) == 148, len(tmpl)
nt_rows = []
for s, bits in streams.items():
    for t in tmpl:
        r = N.non_overlapping_template(bits, template=tuple(t))
        nt_rows.append(dict(stream=s, template="".join(map(str, t)), p_value=r["p_values"][0] if r["applicable"] else np.nan))
nt = pd.DataFrame(nt_rows)
nt.to_csv(RES / "4_nist_nonoverlapping_148_templates.csv", index=False)
nt_sum = nt.groupby("stream").p_value.agg(templates="count", passing=lambda v: int((v >= ALPHA_NIST).sum()), min_p="min")
out("\nNon-overlapping template, all 148 aperiodic m=9 templates:"); out(nt_sum.to_string())

order = [t.__name__ for t in N.ALL_TESTS]
names = ["Frequency", "BlockFrequency", "Runs", "LongestRun", "BinaryMatrixRank", "DFT", "NonOverlappingTemplate",
         "OverlappingTemplate", "MaurerUniversal", "LinearComplexity", "Serial", "ApproximateEntropy",
         "CumulativeSums", "RandomExcursions", "RandomExcursionsVariant"]
score = []
for t in names:
    sub = nist[nist.test == t]
    combined = sub[sub.stream == "COMBINED"]
    cols = sub[sub.stream != "COMBINED"]
    if (combined.result == "N/A").all():
        score.append(dict(test=t, combined_p_values="", combined_min_p=np.nan, proportion="", m_applicable=0,
                          prop_lower_bound=np.nan, status="N/A", reason=combined.note.iloc[0]))
        continue
    cp = combined.p_value.tolist()
    app = cols[cols.result != "N/A"]
    per_stream = app.groupby("stream").p_value.min()  # a stream passes a test only if all its sub-p-values >= 0.01
    m = len(per_stream); k_pass = int((per_stream >= ALPHA_NIST).sum())
    lb, _ = nist_prop_bounds(m) if m else (np.nan, np.nan)
    prop = k_pass / m if m else np.nan
    ok = (min(cp) >= ALPHA_NIST) and (m > 0 and prop >= lb)
    score.append(dict(test=t, combined_p_values=" / ".join((f"{p:.6f}" if p >= 1e-4 else f"{p:.3e}") for p in cp), combined_min_p=min(cp),
                      proportion=f"{k_pass}/{m} = {prop:.4f}" if m else "no column stream applicable",
                      m_applicable=m, prop_lower_bound=lb, status="PASS" if ok else "FAIL",
                      reason=("" if t != "DFT" else "implementation not validated against spec example (see validation)")
                      + ("" if m == len(COLS) else f"; {len(COLS) - m} column stream(s) below minimum length")))
sc = pd.DataFrame(score)
sc.to_csv(RES / "4_nist_scorecard.csv", index=False)
n_pass = int((sc.status == "PASS").sum()); n_fail = int((sc.status == "FAIL").sum()); n_na = int((sc.status == "N/A").sum())
out("\nSCORECARD"); out(sc.to_string(index=False))
out(f"\nNIST battery: PASS {n_pass}/15, FAIL {n_fail}/15, N/A (not runnable at this length) {n_na}/15")

summary = dict(N=n_obs, columns=COLS, proportion_bounds_m7=[lo_b, hi_b],
               uniformity=dict(prop_chi2=prop_chi, prop_ks=prop_ks),
               independence=dict(prop_runs_updown=prop_ud, prop_runs_median=prop_md),
               nist=dict(pass_=n_pass, fail=n_fail, na=n_na))
(RES / "summary.json").write_text(json.dumps(summary, indent=2))
(RES / "machine_output.txt").write_text("\n".join(log) + "\n")
