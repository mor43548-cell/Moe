"""Step 2a: validate scripts/sp800_22.py before touching the data.

(a) Worked examples printed in SP 800-22 Rev.1a sections 2.x.4 / 2.x.8
(b) The expansion of e (first 1,000,000 bits incl. the integer part "10",
    i.e. the NIST `data.e` reference) against the reference P-values quoted in SP 800-22 2.x.8.
The official reference files could not be downloaded (csrc.nist.gov blocked by
this environment's network), so e is computed exactly with integer arithmetic
and the reference P-values are the ones printed in the specification.
"""
import json, math, pathlib, sys
import numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).parent))
import sp800_22 as T

ROOT = pathlib.Path(__file__).resolve().parents[1]
bits = lambda s: np.array([int(c) for c in s if c in "01"], dtype=np.int8)


def e_bits(n, guard=64):
    one = 1 << (n + guard)  # computes floor(e * 2^n), n+2 bits
    s, term, k = 0, one, 0
    while term:
        s += term
        k += 1
        term //= k
    # NIST data.e = binary expansion of e INCLUDING the integer part "10"
    return np.array([int(c) for c in bin(s >> guard)[2:][:n]], dtype=np.int8)


def pi_bits(n, guard=64):  # Machin: pi = 16 atan(1/5) - 4 atan(1/239)
    P = n + guard
    def atan_inv(x):
        s, t, k, x2 = 0, (1 << P) // x, 1, x * x
        while t:
            s += t // k if k % 4 == 1 else -(t // k)
            t //= x2
            k += 2
        return s
    v = (16 * atan_inv(5) - 4 * atan_inv(239)) >> guard  # pi * 2^n
    return np.array([int(c) for c in bin(v)[2:]], dtype=np.int8)  # includes leading "11"


PI100 = pi_bits(98)[:100]
assert "".join(map(str, PI100)).startswith("1100100100001111110110101010001000"), "pi expansion"

rows = []
def chk(label, got, ref, tol=1e-4):
    ok = abs(got - ref) <= tol
    rows.append(dict(case=label, computed=round(got, 6), reference=ref, match=ok))

# (a) small worked examples from the spec
chk("2.1.4 Frequency eps=1011010101", T.frequency(bits("1011010101"))["p_values"][0], 0.527089)
chk("2.1.8 Frequency pi-100", T.frequency(PI100)["p_values"][0], 0.109599)
chk("2.2.4 BlockFreq eps=0110011010 M=3", T.block_frequency(bits("0110011010"), M=3)["p_values"][0], 0.801252)
chk("2.2.8 BlockFreq pi-100 M=10", T.block_frequency(PI100, M=10)["p_values"][0], 0.706438)
chk("2.3.4 Runs eps=1001101011", T.runs(bits("1001101011"))["p_values"][0], 0.147232)
chk("2.3.8 Runs pi-100", T.runs(PI100)["p_values"][0], 0.500798)
chk("2.4.8 LongestRun 128-bit example", T.longest_run(bits(
    "11001100000101010110110001001100111000000000001001001101010100010001001111010110100000001101011111001100111001101101100010110010"))["p_values"][0], 0.180609, tol=2e-4)
# DFT: the two short worked examples printed in 2.6.4 / 2.6.8 do not reproduce with the
# STS 2.1.2 algorithm (threshold sqrt(ln(20) n), N1 over n/2 moduli); they are shown, and the
# full-length e reference (which is what STS reproduces) is checked further below.
chk("2.6.4 DFT eps=1001010011 (printed example)", T.dft(bits("1001010011"))["p_values"][0], 0.029523)
chk("2.6.8 DFT pi-100 (printed example)", T.dft(PI100)["p_values"][0], 0.168669)
chk("2.7.4 NonOverlapping eps=10100100101110010110 B=001 N=2",
    T.non_overlapping_template(bits("10100100101110010110"), m=3, N=2, templates=[[0, 0, 1]])["p_values"][0], 0.344154)
chk("2.11.4 Serial eps=0011011101 m=3 P1", T.serial(bits("0011011101"), m=3)["p_values"][0], 0.808792)
chk("2.11.4 Serial eps=0011011101 m=3 P2", T.serial(bits("0011011101"), m=3)["p_values"][1], 0.670320)
chk("2.12.4 ApEn eps=0100110101 m=3", T.approximate_entropy(bits("0100110101"), m=3)["p_values"][0], 0.261961)
chk("2.12.8 ApEn pi-100 m=2", T.approximate_entropy(PI100, m=2)["p_values"][0], 0.235301)
chk("2.13.4 Cusum eps=1011010111 forward", T.cumulative_sums(bits("1011010111"))["p_values"][0], 0.4116588)
chk("2.13.8 Cusum pi-100 forward", T.cumulative_sums(PI100)["p_values"][0], 0.219194)
chk("2.13.8 Cusum pi-100 reverse", T.cumulative_sums(PI100)["p_values"][1], 0.114866)
ap = T.aperiodic_templates(9)
rows.append(dict(case="number of aperiodic templates m=9 (STS: 148)", computed=len(ap), reference=148, match=len(ap) == 148))

# (b) e, n = 1,000,000
E = e_bits(1_000_000)
assert "".join(map(str, E[:20])) == "10101101111110000101", "e expansion"
chk("e Frequency", T.frequency(E)["p_values"][0], 0.953749)
chk("e Block Frequency M=128", T.block_frequency(E, 128)["p_values"][0], 0.211072)
chk("e Runs", T.runs(E)["p_values"][0], 0.561917)
chk("e Longest Run M=10^4", T.longest_run(E)["p_values"][0], 0.718945)
chk("e Binary Matrix Rank", T.binary_matrix_rank(E)["p_values"][0], 0.306156)
chk("e DFT (full-length reference)", T.dft(E)["p_values"][0], 0.847187)
ov = T.overlapping_template(E)
rows.append(dict(case="e Overlapping Template m=9: block counts v (spec 2.8.8)", computed=str(ov["stats"]["v"]),
                 reference="[329, 164, 150, 111, 78, 136]", match=ov["stats"]["v"] == [329, 164, 150, 111, 78, 136]))
chk("e Overlapping Template m=9, with the pi printed in the 2.8.8 example",
    T.overlapping_template(E, pi=[0.367879, 0.183940, 0.137955, 0.099634, 0.069935, 0.140657])["p_values"][0], 0.110434)
rows.append(dict(case="e Overlapping Template m=9, STS 2.1.2 corrected pi (used for the data)", computed=round(ov["p_values"][0], 6),
                 reference="n/a (spec example predates the corrected pi)", match=None))
chk("e Universal L=7 Q=1280", T.universal(E)["p_values"][0], 0.282568)
chk("e Linear Complexity M=1000", T.linear_complexity(E, 1000)["p_values"][0], 0.845406)
s = T.serial(E, m=2)["p_values"]
chk("e Serial m=2 P1", s[0], 0.843764); chk("e Serial m=2 P2", s[1], 0.561915)
chk("e Approximate Entropy m=10", T.approximate_entropy(E, m=10)["p_values"][0], 0.700073)
c = T.cumulative_sums(E)["p_values"]
chk("e Cusum forward", c[0], 0.669887); chk("e Cusum reverse", c[1], 0.724266)
re = T.random_excursions(E)
rows.append(dict(case="e Random Excursions J (cycles)", computed=re["stats"]["J"], reference=1490, match=re["stats"]["J"] == 1490))
for x, ref in zip(re["params"]["states"], [0.573306, 0.197996, 0.164011, 0.007779, 0.786868, 0.440912, 0.797854, 0.778186]):
    chk(f"e Random Excursions x={x:+d}", re["p_values"][re["params"]["states"].index(x)], ref)
rv = T.random_excursions_variant(E)
refv = [0.858946, 0.794755, 0.576249, 0.493417, 0.633873, 0.917283, 0.934708, 0.816012, 0.826009,
        0.137861, 0.200642, 0.441254, 0.939291, 0.505683, 0.445935, 0.512207, 0.538635, 0.593930]
for x, p, ref in zip(rv["params"]["states"], rv["p_values"], refv):
    chk(f"e Random Excursions Variant x={x:+d}", p, ref)
no = T.non_overlapping_template(E[: 2 ** 20] if len(E) >= 2 ** 20 else E, m=9, N=8, templates=[[0] * 8 + [1]])
chk("e Non-overlapping B=000000001 (first 10^6 bits)", no["p_values"][0], 0.078790)

out = ROOT / "output" / "machine_output"
out.mkdir(parents=True, exist_ok=True)
(out / "02_validation.json").write_text(json.dumps(rows, indent=2))
w = max(len(r["case"]) for r in rows)
for r in rows:
    tag = "INFO" if r["match"] is None else ("MATCH" if r["match"] else "MISMATCH")
    print(f"{r['case']:<{w}}  computed={r['computed']!s:<10} reference={r['reference']!s:<10} {tag}")
scored = [r for r in rows if r["match"] is not None]
print(f"\n{sum(r['match'] for r in scored)}/{len(scored)} match")
