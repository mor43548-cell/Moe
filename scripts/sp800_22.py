"""NIST SP 800-22 Rev.1a statistical test suite -- the 15 tests, implemented
directly from the specification text (sections 2.1-2.15), with the STS 2.1.2
default parameters where the spec leaves a choice.

Every function returns a dict with:
  name, p_values (list), applicable (bool), note (str), params (dict)
A test is reported `applicable=False` ONLY where the specification itself
says the test is not to be applied (Runs pre-test is reported as p=0 per spec;
Random Excursions J<500; Universal n<387,840 because no (L,Q) is defined).
Inputs below NIST's *recommended* minimum n are still computed and the
recommendation is reported separately (min_n_recommended).

Validated by scripts/02_validate_sp800_22.py against the worked examples in
SP 800-22 and the Appendix B reference p-values for the expansion of e.
"""
import math
import numpy as np
from scipy.special import gammaincc, erfc
from scipy.stats import norm

ALPHA = 0.01


def igamc(a, x):
    return float(gammaincc(a, x))


def _res(name, p, params, applicable=True, note="", min_n=None, stats=None):
    p = [float(v) for v in (p if isinstance(p, (list, tuple, np.ndarray)) else [p])]
    return dict(name=name, p_values=p, applicable=applicable, note=note,
                params=params, min_n_recommended=min_n, stats=stats or {})


# 2.1 ---------------------------------------------------------------------
def frequency(e):
    n = len(e)
    s = np.sum(2 * e.astype(int) - 1)
    sobs = abs(s) / math.sqrt(n)
    return _res("Frequency (Monobit)", erfc(sobs / math.sqrt(2)), {}, min_n=100,
                stats={"S_n": int(s), "ones": int(e.sum()), "n": n})


# 2.2 ---------------------------------------------------------------------
def block_frequency(e, M=128):
    n = len(e)
    N = n // M
    blocks = e[: N * M].reshape(N, M)
    pi = blocks.mean(axis=1)
    chi2 = 4 * M * np.sum((pi - 0.5) ** 2)
    return _res("Block Frequency", igamc(N / 2, chi2 / 2), {"M": M, "N": N}, min_n=100,
                stats={"chi2": float(chi2)})


# 2.3 ---------------------------------------------------------------------
def runs(e):
    n = len(e)
    pi = e.mean()
    tau = 2 / math.sqrt(n)
    if abs(pi - 0.5) >= tau:
        return _res("Runs", 0.0, {}, note=f"frequency pre-test failed (|pi-1/2|={abs(pi-0.5):.5f} >= tau={tau:.5f}); per spec 2.3.4(2) P-value = 0.0000",
                    min_n=100, stats={"pi": float(pi)})
    v = 1 + int(np.sum(e[1:] != e[:-1]))
    p = erfc(abs(v - 2 * n * pi * (1 - pi)) / (2 * math.sqrt(2 * n) * pi * (1 - pi)))
    return _res("Runs", p, {}, min_n=100, stats={"V_obs": v, "pi": float(pi)})


# 2.4 ---------------------------------------------------------------------
def longest_run(e):
    n = len(e)
    if n < 128:
        return _res("Longest Run of Ones", [], {}, applicable=False, note="n<128", min_n=128)
    if n < 6272:
        M, K, edges = 8, 3, (1, 4)
        pi = [0.2148, 0.3672, 0.2305, 0.1875]
    elif n < 750000:
        M, K, edges = 128, 5, (4, 9)
        pi = [0.1174, 0.2430, 0.2493, 0.1752, 0.1027, 0.1124]
    else:
        M, K, edges = 10000, 6, (10, 16)
        pi = [0.0882, 0.2092, 0.2483, 0.1933, 0.1208, 0.0675, 0.0727]
    N = n // M
    v = np.zeros(K + 1, dtype=int)
    for b in e[: N * M].reshape(N, M):
        run = best = 0
        for bit in b:
            run = run + 1 if bit else 0
            best = max(best, run)
        idx = min(max(best, edges[0]), edges[1]) - edges[0]
        v[idx] += 1
    pi = np.array(pi)
    chi2 = np.sum((v - N * pi) ** 2 / (N * pi))
    return _res("Longest Run of Ones", igamc(K / 2, chi2 / 2), {"M": M, "K": K, "N": N}, min_n=128,
                stats={"chi2": float(chi2), "v": v.tolist()})


# 2.5 ---------------------------------------------------------------------
def _gf2_rank(mat):
    m = mat.copy().astype(np.uint8)
    rows, cols = m.shape
    r = 0
    for c in range(cols):
        piv = None
        for i in range(r, rows):
            if m[i, c]:
                piv = i
                break
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


def binary_matrix_rank(e, M=32, Q=32):
    n = len(e)
    N = n // (M * Q)
    if N == 0:
        return _res("Binary Matrix Rank", [], {"M": M, "Q": Q, "N": 0}, applicable=False,
                    note="fewer bits than one 32x32 matrix", min_n=38 * M * Q)
    ranks = [_gf2_rank(e[i * M * Q:(i + 1) * M * Q].reshape(M, Q)) for i in range(N)]
    FM = sum(r == M for r in ranks)
    FM1 = sum(r == M - 1 for r in ranks)
    rest = N - FM - FM1
    def _prob(r):  # exact rank probabilities, as computed by STS 2.1.2 (rank.c)
        pr = 2.0 ** (r * (Q + M - r) - M * Q)
        for i in range(r):
            pr *= (1 - 2.0 ** (i - Q)) * (1 - 2.0 ** (i - M)) / (1 - 2.0 ** (i - r))
        return pr
    p32, p31 = _prob(M), _prob(M - 1)
    p30 = 1 - p32 - p31
    chi2 = (FM - p32 * N) ** 2 / (p32 * N) + (FM1 - p31 * N) ** 2 / (p31 * N) + (rest - p30 * N) ** 2 / (p30 * N)
    note = "" if N >= 38 else f"N={N} matrices < 38 recommended by SP 800-22 2.5.7"
    return _res("Binary Matrix Rank", math.exp(-chi2 / 2), {"M": M, "Q": Q, "N": N}, note=note,
                min_n=38 * M * Q, stats={"chi2": chi2, "F_M": FM, "F_M-1": FM1, "rest": rest})


# 2.6 ---------------------------------------------------------------------
def dft(e):
    n = len(e)
    X = 2 * e.astype(float) - 1
    S = np.fft.fft(X)
    Mod = np.abs(S[: n // 2])
    T = math.sqrt(math.log(1 / 0.05) * n)
    N0 = 0.95 * n / 2
    N1 = int(np.sum(Mod < T))
    d = (N1 - N0) / math.sqrt(n * 0.95 * 0.05 / 4)
    return _res("Discrete Fourier Transform (Spectral)", erfc(abs(d) / math.sqrt(2)), {}, min_n=1000,
                stats={"N1": N1, "N0": N0, "d": d})


# 2.7 ---------------------------------------------------------------------
def aperiodic_templates(m):
    out = []
    for v in range(2 ** m):
        b = [(v >> (m - 1 - i)) & 1 for i in range(m)]
        if all(b[k:] != b[: m - k] for k in range(1, m)):
            out.append(b)
    return out


def non_overlapping_template(e, m=9, N=8, templates=None):
    n = len(e)
    M = n // N
    mu = (M - m + 1) / 2 ** m
    var = M * (1 / 2 ** m - (2 * m - 1) / 2 ** (2 * m))
    temps = templates if templates is not None else aperiodic_templates(m)
    ps, worst = [], None
    for B in temps:
        B = np.array(B, dtype=e.dtype)
        W = []
        for j in range(N):
            blk = e[j * M:(j + 1) * M]
            i = cnt = 0
            while i <= M - m:
                if np.array_equal(blk[i:i + m], B):
                    cnt += 1
                    i += m
                else:
                    i += 1
            W.append(cnt)
        chi2 = sum((w - mu) ** 2 for w in W) / var
        ps.append(igamc(N / 2, chi2 / 2))
    return _res("Non-overlapping Template Matching", ps, {"m": m, "N": N, "M": M, "templates": len(temps)},
                min_n=None, note="one P-value per aperiodic template")


# 2.8 ---------------------------------------------------------------------
def overlapping_template(e, m=9, M=1032, K=5, pi=None):
    n = len(e)
    N = n // M
    B = np.ones(m, dtype=e.dtype)
    v = np.zeros(K + 1, dtype=int)
    for j in range(N):
        blk = e[j * M:(j + 1) * M]
        cnt = sum(1 for i in range(M - m + 1) if np.array_equal(blk[i:i + m], B))
        v[min(cnt, K)] += 1
    # corrected class probabilities used by STS 2.1.2 (overlappingTemplateMatchings.c)
    pi = np.array(pi if pi is not None else [0.364091, 0.185659, 0.139381, 0.100571, 0.070432, 0.139865])
    chi2 = np.sum((v - N * pi) ** 2 / (N * pi))
    note = "" if n >= 10 ** 6 else f"n={n} < 10^6 recommended by SP 800-22 2.8.7 (N={N} blocks)"
    return _res("Overlapping Template Matching", igamc(K / 2, chi2 / 2), {"m": m, "M": M, "N": N, "K": K},
                note=note, min_n=10 ** 6, stats={"chi2": float(chi2), "v": v.tolist()})


# 2.9 ---------------------------------------------------------------------
_UNIV = [(387840, 6), (904960, 7), (2068480, 8), (4654080, 9), (10342400, 10),
         (22753280, 11), (49643520, 12), (107560960, 13), (231669760, 14),
         (496435200, 15), (1059061760, 16)]
_UNIV_EV = {6: (5.2177052, 2.954), 7: (6.1962507, 3.125), 8: (7.1836656, 3.238),
            9: (8.1764248, 3.311), 10: (9.1723243, 3.356), 11: (10.170032, 3.384),
            12: (11.168765, 3.401), 13: (12.168070, 3.410), 14: (13.167693, 3.416),
            15: (14.167488, 3.419), 16: (15.167379, 3.421)}


def universal(e, L=None, Q=None):
    n = len(e)
    if L is None:
        L = None
        for nmin, l in _UNIV:
            if n >= nmin:
                L = l
        if L is None:
            return _res("Maurer's Universal Statistical", [], {}, applicable=False,
                        note=f"n={n} < 387,840: SP 800-22 2.9.7 defines no (L,Q) below this length",
                        min_n=387840)
    Q = Q or 10 * 2 ** L
    K = n // L - Q
    T = np.zeros(2 ** L, dtype=np.int64)
    w = e[: (Q + K) * L].reshape(Q + K, L).astype(np.int64) @ (1 << np.arange(L - 1, -1, -1))
    for i in range(Q):
        T[w[i]] = i + 1
    s = 0.0
    for i in range(Q, Q + K):
        s += math.log2(i + 1 - T[w[i]])
        T[w[i]] = i + 1
    fn = s / K
    ev, var = _UNIV_EV[L]
    c = 0.7 - 0.8 / L + (4 + 32 / L) * K ** (-3 / L) / 15
    sigma = c * math.sqrt(var / K)
    return _res("Maurer's Universal Statistical", erfc(abs((fn - ev) / (math.sqrt(2) * sigma))),
                {"L": L, "Q": Q, "K": K}, min_n=387840, stats={"fn": fn})


# 2.10 --------------------------------------------------------------------
def _berlekamp_massey(s):
    n = len(s)
    b = np.zeros(n, dtype=np.uint8); c = np.zeros(n, dtype=np.uint8)
    b[0] = c[0] = 1
    L, m = 0, -1
    for N in range(n):
        d = s[N]
        for i in range(1, L + 1):
            d ^= c[i] & s[N - i]
        if d:
            t = c.copy()
            c[N - m:] ^= b[: n - (N - m)]
            if L <= N // 2:
                L, m, b = N + 1 - L, N, t
    return L


def linear_complexity(e, M=500):
    n = len(e)
    N = n // M
    mu = M / 2 + (9 + (-1) ** (M + 1)) / 36 - (M / 3 + 2 / 9) / 2 ** M
    v = np.zeros(7, dtype=int)
    for j in range(N):
        L = _berlekamp_massey(e[j * M:(j + 1) * M].astype(np.uint8))
        T = (-1) ** M * (L - mu) + 2 / 9
        if T <= -2.5: v[0] += 1
        elif T <= -1.5: v[1] += 1
        elif T <= -0.5: v[2] += 1
        elif T <= 0.5: v[3] += 1
        elif T <= 1.5: v[4] += 1
        elif T <= 2.5: v[5] += 1
        else: v[6] += 1
    pi = np.array([0.01047, 0.03125, 0.125, 0.5, 0.25, 0.0625, 0.020833])  # STS 2.1.2 constants
    chi2 = np.sum((v - N * pi) ** 2 / (N * pi))
    note = "" if N >= 200 else f"N={N} blocks < 200 recommended by SP 800-22 2.10.7"
    return _res("Linear Complexity", igamc(3, chi2 / 2), {"M": M, "N": N}, note=note,
                min_n=10 ** 6, stats={"chi2": float(chi2), "v": v.tolist()})


# 2.11 / 2.12 helpers -----------------------------------------------------
def _pattern_counts(e, m):
    if m <= 0:
        return np.zeros(1)
    n = len(e)
    ext = np.concatenate([e, e[: m - 1]]).astype(np.int64)
    idx = np.zeros(n, dtype=np.int64)
    for k in range(m):
        idx = (idx << 1) | ext[k:k + n]
    return np.bincount(idx, minlength=2 ** m)


def _psi2(e, m):
    if m <= 0:
        return 0.0
    n = len(e)
    c = _pattern_counts(e, m)
    return (2 ** m / n) * float(np.sum(c.astype(float) ** 2)) - n


def serial(e, m=None):
    n = len(e)
    if m is None:
        m = int(math.floor(math.log2(n))) - 3  # largest m with m < floor(log2 n) - 2
    p0, p1, p2 = _psi2(e, m), _psi2(e, m - 1), _psi2(e, m - 2)
    d1 = p0 - p1
    d2 = p0 - 2 * p1 + p2
    return _res("Serial", [igamc(2 ** (m - 2), d1 / 2), igamc(2 ** (m - 3), d2 / 2)], {"m": m},
                note="P-value1 (del psi^2), P-value2 (del^2 psi^2)", stats={"d1": d1, "d2": d2})


def approximate_entropy(e, m=None):
    n = len(e)
    if m is None:
        m = int(math.floor(math.log2(n))) - 6  # largest m with m < floor(log2 n) - 5
    def phi(mm):
        if mm == 0:
            return 0.0
        c = _pattern_counts(e, mm) / n
        c = c[c > 0]
        return float(np.sum(c * np.log(c)))
    ap = phi(m) - phi(m + 1)
    chi2 = 2 * n * (math.log(2) - ap)
    return _res("Approximate Entropy", igamc(2 ** (m - 1), chi2 / 2), {"m": m},
                stats={"ApEn": ap, "chi2": chi2})


# 2.13 --------------------------------------------------------------------
def _cusum_p(n, z):
    s1 = 0.0
    for k in range(int((-n / z + 1) / 4), int((n / z - 1) / 4) + 1):
        s1 += norm.cdf((4 * k + 1) * z / math.sqrt(n)) - norm.cdf((4 * k - 1) * z / math.sqrt(n))
    s2 = 0.0
    for k in range(int((-n / z - 3) / 4), int((n / z - 1) / 4) + 1):
        s2 += norm.cdf((4 * k + 3) * z / math.sqrt(n)) - norm.cdf((4 * k + 1) * z / math.sqrt(n))
    return 1 - s1 + s2


def cumulative_sums(e):
    n = len(e)
    X = 2 * e.astype(int) - 1
    zf = int(np.max(np.abs(np.cumsum(X))))
    zr = int(np.max(np.abs(np.cumsum(X[::-1]))))
    return _res("Cumulative Sums (Cusum)", [_cusum_p(n, zf), _cusum_p(n, zr)], {}, min_n=100,
                note="forward, reverse", stats={"z_forward": zf, "z_reverse": zr})


# 2.14 / 2.15 -------------------------------------------------------------
def _walk(e):
    S = np.cumsum(2 * e.astype(int) - 1)
    Sp = np.concatenate([[0], S, [0]])
    zeros = np.flatnonzero(Sp == 0)
    J = len(zeros) - 1
    return Sp, zeros, J


def random_excursions(e):
    n = len(e)
    Sp, zeros, J = _walk(e)
    xs = [-4, -3, -2, -1, 1, 2, 3, 4]
    req = max(0.005 * math.sqrt(n), 500)
    if J < req:
        return _res("Random Excursions", [], {"J": J}, applicable=False,
                    note=f"J={J} cycles < {req:.0f}: per SP 800-22 2.14.4(6)/2.14.7 the test is not applied",
                    min_n=10 ** 6, stats={"J": J})
    ps = []
    for x in xs:
        ax = abs(x)
        pi = [1 - 1 / (2 * ax)] + [(1 / (4 * ax * ax)) * (1 - 1 / (2 * ax)) ** (k - 1) for k in range(1, 5)] + \
             [(1 / (2 * ax)) * (1 - 1 / (2 * ax)) ** 4]
        v = np.zeros(6, dtype=int)
        for a, b in zip(zeros[:-1], zeros[1:]):
            cnt = int(np.sum(Sp[a:b + 1] == x))
            v[min(cnt, 5)] += 1
        pi = np.array(pi)
        chi2 = np.sum((v - J * pi) ** 2 / (J * pi))
        ps.append(igamc(5 / 2, chi2 / 2))
    return _res("Random Excursions", ps, {"J": J, "states": xs}, min_n=10 ** 6, stats={"J": J})


def random_excursions_variant(e):
    n = len(e)
    Sp, zeros, J = _walk(e)
    xs = [x for x in range(-9, 10) if x != 0]
    req = max(0.005 * math.sqrt(n), 500)
    if J < req:
        return _res("Random Excursions Variant", [], {"J": J}, applicable=False,
                    note=f"J={J} cycles < {req:.0f}: per SP 800-22 2.15.4(5)/2.15.7 the test is not applied",
                    min_n=10 ** 6, stats={"J": J})
    ps = []
    for x in xs:
        xi = int(np.sum(Sp == x))
        ps.append(erfc(abs(xi - J) / math.sqrt(2 * J * (4 * abs(x) - 2))))
    return _res("Random Excursions Variant", ps, {"J": J, "states": xs}, min_n=10 ** 6, stats={"J": J})


ALL = [frequency, block_frequency, runs, longest_run, binary_matrix_rank, dft,
       non_overlapping_template, overlapping_template, universal, linear_complexity,
       serial, approximate_entropy, cumulative_sums, random_excursions, random_excursions_variant]
