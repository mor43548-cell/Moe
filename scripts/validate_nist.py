"""
Validation of scripts/nist_sp800_22.py against the worked examples printed in
NIST SP 800-22 Rev.1a, section 2.x.8 ("Example"), which use the first 100 bits of
the binary expansion of pi (generated here from its hexadecimal expansion, not typed in).

Also an optional cross-check against the independent open-source package `nistrng`
(BSD-3) on the parameter-free tests, if it is importable (pass its directory as argv[1]).
Output: results/validation_nist.txt
"""
import sys, math
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import nist_sp800_22 as N

ROOT = Path(__file__).resolve().parents[1]
lines = []
def out(s=""):
    print(s); lines.append(s)

PI_HEX = "243F6A8885A308D313198A2E03707344A4093822299F31D008"
eps100 = np.array([int(c) for c in ("11" + "".join(f"{int(h, 16):04b}" for h in PI_HEX))[:100]])

def raw(fn, e, *a, **k):
    """Run a test bypassing the minimum-length gate (only for the 10/20-bit spec examples)."""
    return fn(e, *a, **k)

# --- direct formula checks for short examples (gate bypassed by calling formulas) ---
def freq_raw(e):
    n = len(e); s = abs(np.sum(2 * e - 1)) / math.sqrt(n); return float(N.erfc(s / math.sqrt(2)))
def runs_raw(e):
    n = len(e); pi = e.mean(); V = 1 + int(np.sum(e[1:] != e[:-1]))
    return float(N.erfc(abs(V - 2 * n * pi * (1 - pi)) / (2 * math.sqrt(2 * n) * pi * (1 - pi))))
LR128 = "11001100000101010110110001001100111000000000001001001101010100010001001111010110100000001101011111001100111001101101100010110010"
def b(s): return np.array([int(c) for c in s])

checks = [
    ("2.1.4 Frequency, e=1011010101", freq_raw(b("1011010101")), 0.527089),
    ("2.1.8 Frequency, 100 bits of pi", N.frequency(eps100)["p_values"][0], 0.109599),
    ("2.2.8 BlockFrequency M=10, 100 bits of pi", N.block_frequency(eps100, M=10)["p_values"][0], 0.706438),
    ("2.3.4 Runs, e=1001101011", runs_raw(b("1001101011")), 0.147232),
    ("2.3.8 Runs, 100 bits of pi", N.runs(eps100)["p_values"][0], 0.500798),
    ("2.4.8 LongestRun M=8, 128-bit example", N.longest_run(b(LR128))["p_values"][0], 0.180609),
    ("2.6.8 DFT, 100 bits of pi", N.dft(eps100, gate=False)["p_values"][0], 0.168669),
    ("2.7.4 NonOverlapping B=001 N=2, e=10100100101110010110",
     N.non_overlapping_template(b("10100100101110010110"), template=(0, 0, 1), N=2, gate=False)["p_values"][0], 0.344154),
    ("2.11.4 Serial m=3, e=0011011101 (P1)", None, 0.808792),
    ("2.11.4 Serial m=3, e=0011011101 (P2)", None, 0.670320),
    ("2.12.8 ApEn m=2, 100 bits of pi", None, 0.235301),
    ("2.13.8 CuSum forward, 100 bits of pi", N.cumulative_sums(eps100)["p_values"][0], 0.219194),
    ("2.13.8 CuSum backward, 100 bits of pi", N.cumulative_sums(eps100)["p_values"][1], 0.114866),
]
# Serial / ApEn short examples: call internal helpers directly (gate would reject n=10 / n=100 with m given)
e10 = b("0011011101")
p0, p1, p2 = N._psi2(e10, 3), N._psi2(e10, 2), N._psi2(e10, 1)
checks[8] = (checks[8][0], N.igamc(2, (p0 - p1) / 2), 0.808792)
checks[9] = (checks[9][0], N.igamc(1, (p0 - 2 * p1 + p2) / 2), 0.670320)
apen = N._phi(eps100, 2) - N._phi(eps100, 3)
checks[10] = (checks[10][0], N.igamc(2, 2 * 100 * (math.log(2) - apen) / 2), 0.235301)

out("VALIDATION AGAINST NIST SP 800-22 Rev.1a WORKED EXAMPLES")
out(f"{'example':48s} {'computed':>10s} {'NIST':>10s}  match(6dp)")
ok_all = True
for name, got, exp in checks:
    ok = abs(got - exp) < 5e-6
    ok_all &= ok
    out(f"{name:48s} {got:10.6f} {exp:10.6f}  {'YES' if ok else 'NO'}")
out(f"ALL SPEC EXAMPLES REPRODUCED: {ok_all}")
out("NOTE: a 'NO' row is reported as-is. For DFT the implementation follows the formulas of SP 800-22")
out("sections 2.6.4/3.6 literally (moduli of the first n/2 DFT coefficients, T=sqrt(ln(1/0.05)n), N0=0.95n/2);")
out("the counts N1 printed in the spec examples could not be reproduced from those formulas, so DFT results")
out("in this study carry the flag 'implementation not validated against spec example'.")

# --- optional cross-check with nistrng ---
if len(sys.argv) > 1:
    sys.path.insert(0, sys.argv[1])
    try:
        from nistrng.sp800_22r1a import (MonobitTest, RunsTest, LongestRunOnesInABlockTest,
                                         DiscreteFourierTransformTest, CumulativeSumsTest)
        rng = np.random.default_rng(12345)
        out(""); out("CROSS-CHECK vs independent package nistrng 1.2.3 (random reference stream, n=20000)")
        e = rng.integers(0, 2, 20000)
        pairs = [("Frequency", MonobitTest(), N.frequency), ("Runs", RunsTest(), N.runs),
                 ("LongestRun", LongestRunOnesInABlockTest(), N.longest_run),
                 ("DFT", DiscreteFourierTransformTest(), N.dft),
                 ("CumulativeSums", CumulativeSumsTest(), N.cumulative_sums)]
        for name, t, f in pairs:
            r = t.run(e.astype(np.int8))
            theirs = np.atleast_1d(r.score).astype(float).tolist()
            mine = f(e)["p_values"]
            out(f"{name:16s} ours={['%.6f' % x for x in mine]}  nistrng={['%.6f' % x for x in theirs]}")
    except Exception as ex:  # report, do not hide
        out(f"nistrng cross-check could not run: {ex!r}")

(ROOT / "results").mkdir(exist_ok=True)
(ROOT / "results/validation_nist.txt").write_text("\n".join(lines) + "\n")
