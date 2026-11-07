"""
NIST SP 800-22 Rev.1a statistical test suite - Python implementation of the 15 tests.

Formulas follow NIST SP 800-22 Rev.1a (Rukhin et al., 2010), section 2 and 3, and the
reference C code (sts-2.1.2) where the document leaves details open.

Applicability policy (applied identically to every stream, decided BEFORE looking at results):
a test is only evaluated when the input-size requirement stated in SP 800-22 section 2.x.7
is met. Otherwise the test returns status "N/A" with the reason, and no p-value.
"""
from __future__ import annotations
import math
import numpy as np
from scipy.special import gammaincc, erfc
from scipy.stats import norm

ALPHA = 0.01


def igamc(a, x):
    return float(gammaincc(a, x))


def _res(name, p, applicable=True, reason="", **stats):
    if not applicable:
        return {"test": name, "p_values": [], "applicable": False, "reason": reason, "stats": stats}
    ps = p if isinstance(p, list) else [p]
    return {"test": name, "p_values": [float(v) for v in ps], "applicable": True, "reason": reason, "stats": stats}


def _bits(e):
    e = np.asarray(e, dtype=np.int64)
    assert set(np.unique(e)).issubset({0, 1})
    return e


# 2.1 Frequency (Monobit)
def frequency(e):
    e = _bits(e); n = len(e)
    if n < 100:
        return _res("Frequency", None, False, f"n={n} < 100")
    s = int(np.sum(2 * e - 1))
    s_obs = abs(s) / math.sqrt(n)
    return _res("Frequency", float(erfc(s_obs / math.sqrt(2))), n=n, S_n=s, s_obs=s_obs)


# 2.2 Frequency within a block
def block_frequency(e, M=None):
    e = _bits(e); n = len(e)
    if n < 100:
        return _res("BlockFrequency", None, False, f"n={n} < 100")
    if M is None:  # smallest M with M>=20, M>0.01n, N<100 (SP 800-22 2.2.7)
        M = 20
        while not (M > 0.01 * n and n // M < 100):
            M += 1
    N = n // M
    blocks = e[: N * M].reshape(N, M)
    pi = blocks.mean(axis=1)
    chi2 = 4.0 * M * float(np.sum((pi - 0.5) ** 2))
    return _res("BlockFrequency", igamc(N / 2.0, chi2 / 2.0), n=n, M=M, N=N, chi2=chi2)


# 2.3 Runs
def runs(e):
    e = _bits(e); n = len(e)
    if n < 100:
        return _res("Runs", None, False, f"n={n} < 100")
    pi = e.mean()
    tau = 2.0 / math.sqrt(n)
    if abs(pi - 0.5) >= tau:  # frequency prerequisite failed -> p = 0 (spec 2.3.4 step 2)
        return _res("Runs", 0.0, n=n, pi=pi, tau=tau, prerequisite="failed |pi-1/2|>=tau")
    V = 1 + int(np.sum(e[1:] != e[:-1]))
    p = float(erfc(abs(V - 2 * n * pi * (1 - pi)) / (2 * math.sqrt(2 * n) * pi * (1 - pi))))
    return _res("Runs", p, n=n, pi=pi, V_obs=V)


# 2.4 Longest run of ones in a block
def longest_run(e):
    e = _bits(e); n = len(e)
    if n < 128:
        return _res("LongestRun", None, False, f"n={n} < 128")
    if n < 6272:
        M, K, lo, pis = 8, 3, 1, [0.21484375, 0.3671875, 0.23046875, 0.1875]  # exact values, as in sts-2.1.2
    elif n < 750000:
        M, K, lo, pis = 128, 5, 4, [0.1174035788, 0.242955959, 0.249363483, 0.17517706, 0.102701071, 0.112398847]
    else:
        M, K, lo, pis = 10000, 6, 10, [0.0882, 0.2092, 0.2483, 0.1933, 0.1208, 0.0675, 0.0727]
    N = n // M
    nu = np.zeros(K + 1)
    for b in e[: N * M].reshape(N, M):
        best = cur = 0
        for x in b:
            cur = cur + 1 if x else 0
            best = max(best, cur)
        idx = min(max(best - lo, 0), K)
        nu[idx] += 1
    pis = np.array(pis)
    chi2 = float(np.sum((nu - N * pis) ** 2 / (N * pis)))
    return _res("LongestRun", igamc(K / 2.0, chi2 / 2.0), n=n, M=M, N=N, nu=nu.tolist(), chi2=chi2)


def _gf2_rank(m):
    m = m.copy() % 2
    rows, cols = m.shape
    r = 0
    for c in range(cols):
        piv = None
        for i in range(r, rows):
            if m[i, c]:
                piv = i; break
        if piv is None:
            continue
        m[[r, piv]] = m[[piv, r]]
        for i in range(rows):
            if i != r and m[i, c]:
                m[i] ^= m[r]
        r += 1
        if r == rows:
            break
    return r


# 2.5 Binary matrix rank
def matrix_rank(e, Mq=32):
    e = _bits(e); n = len(e)
    need = 38 * Mq * Mq
    if n < need:
        return _res("BinaryMatrixRank", None, False, f"n={n} < 38*M*Q={need}")
    N = n // (Mq * Mq)
    F = {Mq: 0, Mq - 1: 0, "other": 0}
    for k in range(N):
        r = _gf2_rank(e[k * Mq * Mq:(k + 1) * Mq * Mq].reshape(Mq, Mq).astype(np.uint8))
        F[r if r in (Mq, Mq - 1) else "other"] += 1

    def prob(r):
        prod = 1.0
        for i in range(r):
            prod *= ((1 - 2.0 ** (i - Mq)) ** 2) / (1 - 2.0 ** (i - r))
        return 2.0 ** (r * (2 * Mq - r) - Mq * Mq) * prod
    p32, p31 = prob(Mq), prob(Mq - 1); p30 = 1 - p32 - p31
    chi2 = ((F[Mq] - p32 * N) ** 2 / (p32 * N) + (F[Mq - 1] - p31 * N) ** 2 / (p31 * N)
            + (F["other"] - p30 * N) ** 2 / (p30 * N))
    return _res("BinaryMatrixRank", math.exp(-chi2 / 2.0), n=n, N=N, F=[F[Mq], F[Mq - 1], F["other"]], chi2=chi2)


# 2.6 Discrete Fourier Transform (spectral)
def dft(e, gate=True):
    e = _bits(e); n = len(e)
    if gate and n < 1000:
        return _res("DFT", None, False, f"n={n} < 1000")
    X = 2 * e - 1
    mod = np.abs(np.fft.fft(X))[: n // 2]
    T = math.sqrt(math.log(1 / 0.05) * n)
    N0 = 0.95 * n / 2.0
    N1 = int(np.sum(mod < T))
    d = (N1 - N0) / math.sqrt(n * 0.95 * 0.05 / 4.0)
    return _res("DFT", float(erfc(abs(d) / math.sqrt(2))), n=n, N0=N0, N1=N1, d=d)


def aperiodic_templates(m):
    out = []
    for v in range(2 ** m):
        b = [(v >> (m - 1 - i)) & 1 for i in range(m)]
        if all(b[k:] != b[:-k] for k in range(1, m)):
            out.append(b)
    return out


# 2.7 Non-overlapping template matching
def non_overlapping_template(e, template=(0, 0, 0, 0, 0, 0, 0, 0, 1), N=8, gate=True):
    e = _bits(e); n = len(e)
    B = list(template); m = len(B)
    M = n // N
    if gate and not (M > 0.01 * n and N <= 100 and n >= 100):
        return _res("NonOverlappingTemplate", None, False, f"n={n}: requires N<=100 and M>0.01n")
    W = []
    for j in range(N):
        blk = e[j * M:(j + 1) * M].tolist()
        i = c = 0
        while i <= M - m:
            if blk[i:i + m] == B:
                c += 1; i += m
            else:
                i += 1
        W.append(c)
    mu = (M - m + 1) / 2.0 ** m
    var = M * (1 / 2.0 ** m - (2 * m - 1) / 2.0 ** (2 * m))
    chi2 = float(sum((w - mu) ** 2 for w in W) / var)
    return _res("NonOverlappingTemplate", igamc(N / 2.0, chi2 / 2.0), n=n, template="".join(map(str, B)),
                N=N, M=M, W=W, mu=mu, chi2=chi2)


# 2.8 Overlapping template matching
def overlapping_template(e, m=9, M=1032, K=5):
    e = _bits(e); n = len(e)
    if n < 10 ** 6:
        return _res("OverlappingTemplate", None, False, f"n={n} < 1,000,000")
    N = n // M
    B = [1] * m
    nu = np.zeros(K + 1)
    for j in range(N):
        blk = e[j * M:(j + 1) * M].tolist()
        c = sum(1 for i in range(M - m + 1) if blk[i:i + m] == B)
        nu[min(c, K)] += 1
    pis = np.array([0.364091, 0.185659, 0.139381, 0.100571, 0.070432, 0.139865])
    chi2 = float(np.sum((nu - N * pis) ** 2 / (N * pis)))
    return _res("OverlappingTemplate", igamc(K / 2.0, chi2 / 2.0), n=n, N=N, nu=nu.tolist(), chi2=chi2)


# 2.9 Maurer's universal statistical test
def maurer_universal(e):
    e = _bits(e); n = len(e)
    if n < 387840:
        return _res("MaurerUniversal", None, False, f"n={n} < 387,840")
    table = [(387840, 6), (904960, 7), (2068480, 8), (4654080, 9), (10342400, 10), (22753280, 11),
             (49643520, 12), (107560960, 13), (231669760, 14), (496435200, 15), (1059061760, 16)]
    L = 6
    for lim, l in table:
        if n >= lim:
            L = l
    exp = {6: (5.2177052, 2.954), 7: (6.1962507, 3.125), 8: (7.1836656, 3.238), 9: (8.1764248, 3.311),
           10: (9.1723243, 3.356), 11: (10.170032, 3.384), 12: (11.168765, 3.401), 13: (12.168070, 3.410),
           14: (13.167693, 3.416), 15: (14.167488, 3.419), 16: (15.167379, 3.421)}
    Q = 10 * 2 ** L
    Kb = n // L - Q
    T = np.zeros(2 ** L, dtype=np.int64)
    blocks = e[: (Q + Kb) * L].reshape(Q + Kb, L).dot(1 << np.arange(L - 1, -1, -1))
    for i in range(Q):
        T[blocks[i]] = i + 1
    s = 0.0
    for i in range(Q, Q + Kb):
        s += math.log2(i + 1 - T[blocks[i]]); T[blocks[i]] = i + 1
    fn = s / Kb
    ev, var = exp[L]
    c = 0.7 - 0.8 / L + (4 + 32 / L) * (Kb ** (-3 / L)) / 15
    sigma = c * math.sqrt(var / Kb)
    return _res("MaurerUniversal", float(erfc(abs((fn - ev) / (math.sqrt(2) * sigma)))), n=n, L=L, Q=Q, K=Kb, fn=fn)


def _berlekamp_massey(s):
    n = len(s); b = [0] * n; c = [0] * n; b[0] = c[0] = 1
    L, m = 0, -1
    for N in range(n):
        d = s[N]
        for i in range(1, L + 1):
            d ^= c[i] & s[N - i]
        if d:
            t = c[:]
            for i in range(0, n - N + m):
                c[N - m + i] ^= b[i]
            if L <= N // 2:
                L, m, b = N + 1 - L, N, t
    return L


# 2.10 Linear complexity
def linear_complexity(e, M=500):
    e = _bits(e); n = len(e)
    if n < 10 ** 6:
        return _res("LinearComplexity", None, False, f"n={n} < 1,000,000")
    N = n // M; K = 6
    mu = M / 2.0 + (9 + (-1) ** (M + 1)) / 36.0 - (M / 3.0 + 2 / 9.0) / 2.0 ** M
    nu = np.zeros(7)
    for j in range(N):
        L = _berlekamp_massey(e[j * M:(j + 1) * M].tolist())
        T = (-1) ** M * (L - mu) + 2 / 9.0
        idx = 0 if T <= -2.5 else 6 if T > 2.5 else int(math.floor(T + 3.5)) if T > -2.5 else 0
        idx = min(max(idx, 0), 6)
        nu[idx] += 1
    pis = np.array([0.010417, 0.03125, 0.125, 0.5, 0.25, 0.0625, 0.020833])
    chi2 = float(np.sum((nu - N * pis) ** 2 / (N * pis)))
    return _res("LinearComplexity", igamc(K / 2.0, chi2 / 2.0), n=n, M=M, N=N, nu=nu.tolist(), chi2=chi2)


def _psi2(e, m):
    if m <= 0:
        return 0.0
    n = len(e)
    ext = np.concatenate([e, e[: m - 1]])
    idx = np.zeros(n, dtype=np.int64)
    for k in range(m):
        idx = (idx << 1) | ext[k:k + n]
    counts = np.bincount(idx, minlength=2 ** m)
    return (2.0 ** m / n) * float(np.sum(counts.astype(float) ** 2)) - n


# 2.11 Serial
def serial(e, m=None):
    e = _bits(e); n = len(e)
    lim = int(math.floor(math.log2(n))) - 2
    if m is None:
        m = lim - 1  # largest m satisfying m < floor(log2 n) - 2 (SP 800-22 2.11.7)
    if m < 3 or m >= lim:
        return _res("Serial", None, False, f"m={m} must satisfy 3 <= m < floor(log2 n)-2 = {lim}")
    p0, p1, p2 = _psi2(e, m), _psi2(e, m - 1), _psi2(e, m - 2)
    d1 = p0 - p1; d2 = p0 - 2 * p1 + p2
    return _res("Serial", [igamc(2 ** (m - 2), d1 / 2.0), igamc(2 ** (m - 3), d2 / 2.0)], n=n, m=m, del1=d1, del2=d2)


def _phi(e, m):
    if m == 0:
        return 0.0
    n = len(e)
    ext = np.concatenate([e, e[: m - 1]])
    idx = np.zeros(n, dtype=np.int64)
    for k in range(m):
        idx = (idx << 1) | ext[k:k + n]
    c = np.bincount(idx, minlength=2 ** m).astype(float) / n
    c = c[c > 0]
    return float(np.sum(c * np.log(c)))


# 2.12 Approximate entropy
def approximate_entropy(e, m=None):
    e = _bits(e); n = len(e)
    lim = int(math.floor(math.log2(n))) - 5
    if m is None:
        m = lim - 1  # largest m satisfying m < floor(log2 n) - 5 (SP 800-22 2.12.7)
    if m < 1 or m >= lim:
        return _res("ApproximateEntropy", None, False, f"m={m} must satisfy m < floor(log2 n)-5 = {lim}")
    apen = _phi(e, m) - _phi(e, m + 1)
    chi2 = 2.0 * n * (math.log(2) - apen)
    return _res("ApproximateEntropy", igamc(2 ** (m - 1), chi2 / 2.0), n=n, m=m, ApEn=apen, chi2=chi2)


# 2.13 Cumulative sums
def cumulative_sums(e):
    e = _bits(e); n = len(e)
    if n < 100:
        return _res("CumulativeSums", None, False, f"n={n} < 100")
    X = 2 * e - 1
    out = []
    for seq in (X, X[::-1]):
        z = int(np.max(np.abs(np.cumsum(seq))))
        s1 = 0.0
        k = int((-n / z + 1) / 4)
        while k <= (n / z - 1) / 4:
            s1 += norm.cdf((4 * k + 1) * z / math.sqrt(n)) - norm.cdf((4 * k - 1) * z / math.sqrt(n)); k += 1
        s2 = 0.0
        k = int((-n / z - 3) / 4)
        while k <= (n / z - 1) / 4:
            s2 += norm.cdf((4 * k + 3) * z / math.sqrt(n)) - norm.cdf((4 * k + 1) * z / math.sqrt(n)); k += 1
        out.append(float(1.0 - s1 + s2))
    return _res("CumulativeSums", out, n=n, labels=["forward", "backward"])


def _cycles(e):
    S = np.concatenate([[0], np.cumsum(2 * e - 1), [0]])
    zeros = np.where(S == 0)[0]
    return S, zeros, len(zeros) - 1


# 2.14 Random excursions
def random_excursions(e):
    e = _bits(e); n = len(e)
    S, zeros, J = _cycles(e)
    if n < 10 ** 6 or J < 500:
        return _res("RandomExcursions", None, False, f"n={n} (requires >= 1,000,000); J={J} cycles (requires >= 500)", J=J)
    states = [-4, -3, -2, -1, 1, 2, 3, 4]
    nu = {x: np.zeros(6) for x in states}
    for a, b in zip(zeros[:-1], zeros[1:]):
        cyc = S[a + 1:b]
        for x in states:
            nu[x][min(int(np.sum(cyc == x)), 5)] += 1

    def pi(k, x):
        x = abs(x)
        if k == 0:
            return 1 - 1 / (2 * x)
        if k < 5:
            return (1 / (4 * x * x)) * (1 - 1 / (2 * x)) ** (k - 1)
        return (1 / (2 * x)) * (1 - 1 / (2 * x)) ** 4
    ps = []
    for x in states:
        chi2 = sum((nu[x][k] - J * pi(k, x)) ** 2 / (J * pi(k, x)) for k in range(6))
        ps.append(igamc(2.5, chi2 / 2.0))
    return _res("RandomExcursions", ps, n=n, J=J, states=states)


# 2.15 Random excursions variant
def random_excursions_variant(e):
    e = _bits(e); n = len(e)
    S, zeros, J = _cycles(e)
    if n < 10 ** 6 or J < 500:
        return _res("RandomExcursionsVariant", None, False, f"n={n} (requires >= 1,000,000); J={J} cycles (requires >= 500)", J=J)
    states = [x for x in range(-9, 10) if x != 0]
    ps = []
    for x in states:
        xi = int(np.sum(S == x))
        ps.append(float(erfc(abs(xi - J) / math.sqrt(2 * J * (4 * abs(x) - 2)))))
    return _res("RandomExcursionsVariant", ps, n=n, J=J, states=states)


ALL_TESTS = [frequency, block_frequency, runs, longest_run, matrix_rank, dft, non_overlapping_template,
             overlapping_template, maurer_universal, linear_complexity, serial, approximate_entropy,
             cumulative_sums, random_excursions, random_excursions_variant]


def run_all(e):
    return [t(e) for t in ALL_TESTS]
