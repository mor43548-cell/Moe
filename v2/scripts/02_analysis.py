"""
V2 analysis: the requested test set, applied under the lotto hypothesis exactly as fixed in v2/DESIGN.md.
Every number written to v2/results is the unmodified program output.
"""
import json, math, sys
from math import comb
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import nist_sp800_22 as N

V2 = ROOT / "v2"; RES = V2 / "results"; FIG = V2 / "figures"
RES.mkdir(exist_ok=True); FIG.mkdir(exist_ok=True)

# ---------------- fixed design (v2/DESIGN.md) ----------------
NMAX, K, SMAX = 37, 6, 7
ADJ = 36 / 31            # chi-square scale factor for 6-without-replacement sampling
M_BLOCKS, BLOCK = 8, 41  # 328 = 8 x 41
MC_SIMS, SEED = 20000, 20261107
ALPHA, ALPHA_STRICT, MAX_LAG, Z95, ENT_T, TARGET = 0.01, 0.05, 20, 1.959964, 0.999, 0.99

log = []
def out(s=""):
    print(s); log.append(str(s))

df = pd.read_csv(V2 / "data/merged_lotto_table.csv")
main = df.iloc[:, 1:7].to_numpy()          # sorted ascending
strong = df.iloc[:, 7].to_numpy()
series_names = list(df.columns[1:8])
D = len(df)
out(f"Draws: {D}; main numbers pooled: {main.size}; strong numbers: {len(strong)}")

def bounds(m):
    h = 3 * math.sqrt(TARGET * (1 - TARGET) / m); return TARGET - h, min(1.0, TARGET + h)

def counts_main(m):
    return np.bincount(m.ravel() - 1, minlength=NMAX)

def chi2_stat(c, n_items):
    e = n_items / len(c); return float(np.sum((c - e) ** 2 / e))

def ks_D(c):
    k = len(c); return float(np.max(np.abs(np.cumsum(c) / c.sum() - np.arange(1, k + 1) / k)))

def main_chi_p(m):
    c = counts_main(m); x2 = chi2_stat(c, m.size)
    return x2, float(stats.chi2.sf(x2 * ADJ, NMAX - 1)), float(stats.chi2.sf(x2, NMAX - 1))

# ---------------- Monte Carlo under exact H0 ----------------
rng = np.random.default_rng(SEED)
mc_x2, mc_ks, mc_sx2, mc_sks = [], [], [], []
for _ in range(MC_SIMS // 500):
    keys = rng.random((500, D, NMAX))
    draws = np.argpartition(keys, K, axis=2)[:, :, :K] + 1
    for d in draws:
        c = counts_main(d); mc_x2.append(chi2_stat(c, d.size)); mc_ks.append(ks_D(c))
    s = rng.integers(1, SMAX + 1, size=(500, D))
    for v in s:
        c = np.bincount(v - 1, minlength=SMAX); mc_sx2.append(chi2_stat(c, D)); mc_sks.append(ks_D(c))
mc_x2, mc_ks, mc_sx2, mc_sks = map(np.array, (mc_x2, mc_ks, mc_sx2, mc_sks))
mcp = lambda arr, obs: float((1 + np.sum(arr >= obs - 1e-12)) / (1 + len(arr)))
out(f"Monte Carlo check of the chi-square scale: mean X² under H0 = {mc_x2.mean():.3f} (theory 36*31/36 = 31.000); "
    f"adjusted-asymptotic vs MC 95th percentile: {stats.chi2.ppf(0.95, 36) / ADJ:.3f} vs {np.percentile(mc_x2, 95):.3f}")

# =========================================================
# 1. UNIFORMITY
# =========================================================
out("\n=== 1. UNIFORMITY ===")
cm = counts_main(main); cs = np.bincount(strong - 1, minlength=SMAX)
x2m, p_m_adj, p_m_raw = main_chi_p(main)
ksm = ks_D(cm)
x2s = chi2_stat(cs, D); p_s = float(stats.chi2.sf(x2s, SMAX - 1)); kss = ks_D(cs)
uni = pd.DataFrame([
    dict(variable="Main numbers pooled (1..37)", n=int(main.size), k=NMAX, expected_per_value=main.size / NMAX,
         min_count=int(cm.min()), max_count=int(cm.max()), chi2=x2m, df=NMAX - 1,
         p_chi2_primary=p_m_adj, p_chi2_unadjusted=p_m_raw, p_chi2_MC=mcp(mc_x2, x2m), KS_D=ksm, p_KS_MC=mcp(mc_ks, ksm)),
    dict(variable="Strong number (1..7)", n=int(D), k=SMAX, expected_per_value=D / SMAX,
         min_count=int(cs.min()), max_count=int(cs.max()), chi2=x2s, df=SMAX - 1,
         p_chi2_primary=p_s, p_chi2_unadjusted=p_s, p_chi2_MC=mcp(mc_sx2, x2s), KS_D=kss, p_KS_MC=mcp(mc_sks, kss)),
])
out(uni.to_string(index=False))
uni.to_csv(RES / "1_uniformity.csv", index=False)
pd.DataFrame({"number": range(1, NMAX + 1), "count": cm, "expected": main.size / NMAX}).to_csv(RES / "1_counts_main.csv", index=False)
pd.DataFrame({"strong": range(1, SMAX + 1), "count": cs, "expected": D / SMAX}).to_csv(RES / "1_counts_strong.csv", index=False)

blk = []
for b in range(M_BLOCKS):
    sl = slice(b * BLOCK, (b + 1) * BLOCK)
    _, pm, _ = main_chi_p(main[sl])
    c = np.bincount(strong[sl] - 1, minlength=SMAX); ps = float(stats.chi2.sf(chi2_stat(c, BLOCK), SMAX - 1))
    blk.append(dict(block=b + 1, draws=f"{b * BLOCK + 1}-{(b + 1) * BLOCK}", p_main_chi2=pm, p_strong_chi2=ps))
blk = pd.DataFrame(blk)
lb = bounds(M_BLOCKS)[0]
prop_main = float((blk.p_main_chi2 >= ALPHA).mean()); prop_strong = float((blk.p_strong_chi2 >= ALPHA).mean())
out(blk.to_string(index=False))
out(f"Proportion (blocks with p >= 0.01): main {prop_main:.4f}, strong {prop_strong:.4f}; NIST interval for m=8: >= {lb:.4f}")
blk.to_csv(RES / "1_uniformity_blocks.csv", index=False)

fig, axes = plt.subplots(1, 2, figsize=(18, 5.5), gridspec_kw=dict(width_ratios=[3, 1]))
axes[0].bar(range(1, NMAX + 1), cm, color="#4C72B0", width=0.8, label="observed")
axes[0].axhline(main.size / NMAX, color="#C44E52", ls="--", lw=2, label=f"expected = {main.size / NMAX:.1f}")
sd = math.sqrt(D * (K / NMAX) * (1 - K / NMAX))
axes[0].axhspan(main.size / NMAX - 1.96 * sd, main.size / NMAX + 1.96 * sd, color="grey", alpha=0.15, label="±1.96 SD of a single count under H0")
axes[0].set_xticks(range(1, NMAX + 1)); axes[0].set_xlabel("main number"); axes[0].set_ylabel("count")
axes[0].set_title(f"Main numbers pooled (n = {main.size})   chi² p = {p_m_adj:.4f}   KS p = {uni.p_KS_MC[0]:.4f}"); axes[0].legend(fontsize=9)
axes[1].bar(range(1, SMAX + 1), cs, color="#55A868", width=0.8)
axes[1].axhline(D / SMAX, color="#C44E52", ls="--", lw=2)
axes[1].set_title(f"Strong number   chi² p = {p_s:.4f}"); axes[1].set_xlabel("strong number")
fig.suptitle("Frequency histograms under the lotto hypothesis", fontsize=14)
fig.tight_layout(); fig.savefig(FIG / "1_frequency_histograms.png", dpi=150); plt.close(fig)

# =========================================================
# 2. INDEPENDENCE (series over draws)
# =========================================================
out("\n=== 2. INDEPENDENCE ===")
def acf(x, L):
    x = np.asarray(x, float) - np.mean(x); d = np.sum(x * x)
    return np.array([1.0] + [np.sum(x[:-k] * x[k:]) / d for k in range(1, L + 1)])
def runs_up_down(x):
    s = np.sign(np.diff(np.asarray(x, float))); s = s[s != 0]; n = len(s) + 1
    R = 1 + int(np.sum(s[1:] != s[:-1])); mu = (2 * n - 1) / 3; var = (16 * n - 29) / 90
    z = (R - mu) / math.sqrt(var); return R, mu, z, 2 * stats.norm.sf(abs(z))
def runs_median(x):
    x = np.asarray(x, float); med = np.median(x); s = x[x != med] > med
    n1, n2 = int(s.sum()), int((~s).sum()); n = n1 + n2; R = 1 + int(np.sum(s[1:] != s[:-1]))
    mu = 2 * n1 * n2 / n + 1; var = 2 * n1 * n2 * (2 * n1 * n2 - n) / (n * n * (n - 1))
    z = (R - mu) / math.sqrt(var); return R, mu, z, 2 * stats.norm.sf(abs(z))
band = Z95 / math.sqrt(D)
ser = {nm: df[nm].to_numpy() for nm in series_names}
acfs, rows = {}, []
for nm, x in ser.items():
    r = acf(x, MAX_LAG); acfs[nm] = r
    ob = [k for k in range(1, MAX_LAG + 1) if abs(r[k]) > band]
    R1, m1, z1, p1 = runs_up_down(x); R2, m2, z2, p2 = runs_median(x)
    rows.append(dict(series=nm, acf_lag1=r[1], max_abs_acf=float(np.max(np.abs(r[1:]))), lags_outside_band=len(ob),
                     which_lags=" ".join(map(str, ob)), runs_updown_R=R1, runs_updown_expected=m1, runs_updown_p=p1,
                     runs_median_R=R2, runs_median_expected=m2, runs_median_p=p2))
ind = pd.DataFrame(rows)
out(ind.to_string(index=False)); out(f"95% band = ±{band:.4f}")
ind.to_csv(RES / "2_independence.csv", index=False)
pd.DataFrame(acfs).rename_axis("lag").to_csv(RES / "2_acf_lags_0_20.csv")
lb7 = bounds(7)[0]
prop_ud = float((ind.runs_updown_p >= ALPHA).mean()); prop_md = float((ind.runs_median_p >= ALPHA).mean())

fig, axes = plt.subplots(2, 4, figsize=(18, 8))
for ax, nm in zip(axes.flat, series_names):
    r = acfs[nm]; lags = np.arange(MAX_LAG + 1)
    ax.bar(lags, r, width=0.6, color=["#4C72B0" if k == 0 or abs(r[k]) <= band else "#C44E52" for k in lags])
    ax.axhspan(-band, band, color="grey", alpha=0.2, label=f"95% band ±{band:.3f}"); ax.axhline(0, color="k", lw=0.8)
    ax.set_ylim(-0.4, 1.05); ax.set_title(f"{nm}   lag-1 r = {r[1]:.3f}", fontsize=10); ax.set_xlabel("lag"); ax.legend(fontsize=8)
axes.flat[-1].axis("off")
fig.suptitle("Correlograms of the 7 series over draws (red = outside the 95% band)", fontsize=14)
fig.tight_layout(); fig.savefig(FIG / "2_correlograms.png", dpi=150); plt.close(fig)

# =========================================================
# 3. bits for NIST (rejection sampling) + ENTROPY
# =========================================================
out("\n=== 3. ENTROPY ===")
def colex_rank(c):
    return sum(comb(int(v) - 1, i + 1) for i, v in enumerate(sorted(c)))
assert colex_rank([1, 2, 3, 4, 5, 6]) == 0 and colex_rank(range(32, 38)) == comb(37, 6) - 1
ranks = np.array([colex_rank(r) for r in main])
keep_m = ranks < 2 ** 21; keep_s = (strong - 1) < 4
block_bits = [[] for _ in range(M_BLOCKS)]
for i in range(D):
    bb = []
    if keep_m[i]:
        bb += [(int(ranks[i]) >> (20 - j)) & 1 for j in range(21)]
    if keep_s[i]:
        bb += [((int(strong[i]) - 1) >> (1 - j)) & 1 for j in range(2)]
    block_bits[i // BLOCK] += bb
stream = np.array(sum(block_bits, []), dtype=np.int64)
(RES / "bitstream_v2.txt").write_text("".join(map(str, stream.tolist())) + "\n")
out(f"Bit stream: {len(stream)} bits; main draws used {int(keep_m.sum())}/{D} (expected {D * 2**21 / comb(37, 6):.1f}); "
    f"strong used {int(keep_s.sum())}/{D} (expected {D * 4 / 7:.1f})")

def H(c):
    p = c[c > 0] / c.sum(); return float(-np.sum(p * np.log2(p)))
p1 = stream.mean(); Hb = float(-(p1 * math.log2(p1) + (1 - p1) * math.log2(1 - p1)))
ent = pd.DataFrame([
    dict(variable="Main numbers pooled", H=H(cm), H_max=math.log2(NMAX), ratio=H(cm) / math.log2(NMAX),
         ref_expected_ratio_if_H0=1 - (NMAX - 1) / (2 * main.size * math.log(2)) / math.log2(NMAX)),
    dict(variable="Strong number", H=H(cs), H_max=math.log2(SMAX), ratio=H(cs) / math.log2(SMAX),
         ref_expected_ratio_if_H0=1 - (SMAX - 1) / (2 * D * math.log(2)) / math.log2(SMAX)),
    dict(variable=f"NIST bit stream (ones fraction {p1:.4f})", H=Hb, H_max=1.0, ratio=Hb,
         ref_expected_ratio_if_H0=1 - 1 / (2 * len(stream) * math.log(2))),
])
ent["meets_99_9"] = ent.ratio >= ENT_T
out(ent.to_string(index=False)); ent.to_csv(RES / "3_entropy.csv", index=False)

# =========================================================
# 4. NIST SP 800-22
# =========================================================
out("\n=== 4. NIST SP 800-22 ===")
rows = []
streams = {"WHOLE": stream} | {f"B{b + 1}": np.array(block_bits[b], dtype=np.int64) for b in range(M_BLOCKS)}
for s, bits in streams.items():
    for r in N.run_all(bits):
        if r["applicable"]:
            for i, p in enumerate(r["p_values"]):
                rows.append(dict(stream=s, n_bits=len(bits), test=r["test"], sub=i + 1, p_value=p,
                                 result="SUCCESS" if p >= ALPHA else "FAILURE", note=json.dumps(r["stats"], default=float)))
        else:
            rows.append(dict(stream=s, n_bits=len(bits), test=r["test"], sub=0, p_value=np.nan, result="N/A", note=r["reason"]))
nist = pd.DataFrame(rows); nist.to_csv(RES / "4_nist_machine_output.csv", index=False)
out(nist.drop(columns="note").to_string(index=False))
tmpl = N.aperiodic_templates(9)
nt = pd.DataFrame([dict(stream=s, template="".join(map(str, t)),
                        p_value=(lambda r: r["p_values"][0] if r["applicable"] else np.nan)(N.non_overlapping_template(b, template=tuple(t))))
                   for s, b in streams.items() for t in tmpl])
nt.to_csv(RES / "4_nist_nonoverlapping_148.csv", index=False)

names = ["Frequency", "BlockFrequency", "Runs", "LongestRun", "BinaryMatrixRank", "DFT", "NonOverlappingTemplate",
         "OverlappingTemplate", "MaurerUniversal", "LinearComplexity", "Serial", "ApproximateEntropy",
         "CumulativeSums", "RandomExcursions", "RandomExcursionsVariant"]
sc = []
for t in names:
    sub = nist[nist.test == t]; w = sub[sub.stream == "WHOLE"]; bl = sub[(sub.stream != "WHOLE") & (sub.result != "N/A")]
    if (w.result == "N/A").all():
        sc.append(dict(test=t, p_values="", proportion="", status="N/A", reason=w.note.iloc[0])); continue
    cp = w.p_value.tolist(); per = bl.groupby("stream").p_value.min(); m = len(per); kp = int((per >= ALPHA).sum())
    lbm = bounds(m)[0] if m else np.nan
    ok = min(cp) >= ALPHA and m > 0 and kp / m >= lbm
    note = "DFT implementation not validated against spec example" if t == "DFT" else ""
    if m < M_BLOCKS:
        extra = f"proportion over {m} blocks (others below minimum length)" if m else "no block long enough for a proportion"
        note += ("; " if note else "") + extra
    sc.append(dict(test=t, p_values=" / ".join(f"{p:.6f}" if p >= 1e-4 else f"{p:.3e}" for p in cp),
                   proportion=f"{kp}/{m} = {kp / m:.4f} (>= {lbm:.4f})" if m else "—",
                   status="PASS" if ok else "FAIL", reason=note))
sc = pd.DataFrame(sc); sc.to_csv(RES / "4_nist_scorecard.csv", index=False)
P = {s: int((sc.status == s).sum()) for s in ("PASS", "FAIL", "N/A")}
out("\nSCORECARD"); out(sc.to_string(index=False))
out(f"\nNIST battery: PASS {P['PASS']}/15, FAIL {P['FAIL']}/15, N/A {P['N/A']}/15")

json.dump(dict(prop_bounds_m8=bounds(8), prop_bounds_m7=bounds(7), prop_main=prop_main, prop_strong=prop_strong,
               prop_runs_ud=prop_ud, prop_runs_md=prop_md, nist=P, stream_bits=len(stream),
               used_main=int(keep_m.sum()), used_strong=int(keep_s.sum()),
               mc_mean_x2=float(mc_x2.mean())), open(RES / "summary.json", "w"), indent=2)
(RES / "machine_output.txt").write_text("\n".join(log) + "\n")
