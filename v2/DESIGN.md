# V2 design: written and committed BEFORE any v2 test was run

## Hypothesis under test (H0)
Each row is one draw of a 6/37 lotto with a strong number. The 6 main numbers are a uniformly random 6-subset of {1..37},
stored sorted ascending (N1 = Q, N2 = W, N3 = E, N4 = R, N5 = T, N6 = Y). The strong number (ABC) is uniform on {1..7}
and independent of the main numbers. Draws are independent of each other.

## Test set and thresholds: unchanged from the request (same tests, same thresholds as v1)
Only the representation of the data changes, so that each test's null hypothesis matches H0:

| Requested test | v1 applied it to | v2 applies it to | Why |
|---|---|---|---|
| Frequency chi² / KS | each column on its own [min,max] | (a) all 1,968 main numbers pooled, on 1..37; (b) Strong, on 1..7 | under H0 the sorted columns are order statistics and are non-uniform by construction; each number 1..37 is drawn with the same probability 6/37 |
| Proportion | over 7 columns | over m = 8 consecutive blocks of 41 draws (328 = 8 × 41) | NIST proportion = share of sub-sequences that pass |
| Histogram | per column | pooled main numbers 1..37; Strong 1..7 | as above |
| Autocorrelation, Runs (up/down and above/below median) | each column as a series over draws | identical: the 7 series N1..N6, Strong over draws | under H0 each series is i.i.d. across draws, so these tests are valid without uniformity |
| Shannon entropy | per column | pooled main numbers (max log2 37); Strong (max log2 7); NIST bit stream (max 1) | |
| NIST SP 800-22 (15 tests) | value−min in ⌈log2 k⌉ bits (biased when k is not a power of 2) | unbiased bits: see below | so that NIST's H0 (i.i.d. fair bits) is exactly equivalent to the lotto H0 |

### Chi-square under lotto sampling
Main numbers are drawn without replacement within a draw, so the pooled counts are not multinomial. Their covariance is
(31/36) times the multinomial covariance. Hence X²·(36/31) ~ χ²(36) asymptotically. **Primary p-value = this adjusted
asymptotic p.** The unadjusted χ²(36) p-value and a Monte Carlo p-value under the exact H0 are also reported.
Strong number: ordinary multinomial χ²(6). KS (discrete, Conover): the p-value comes from Monte Carlo under the exact H0,
because no closed form exists for lotto sampling. Monte Carlo: 20,000 replicates, seed 20261107.

### Bits for NIST (unbiased under H0)
- Main draw → its combinatorial rank r = Σ C(c_i − 1, i), i = 1..6 (colex), which is uniform on 0..C(37,6)−1 = 0..2,324,783 under H0.
  If r < 2^21, emit the 21 bits of r (MSB first). Otherwise discard the draw for the bit stream (rejection sampling, ≈ 9.8% expected).
- Strong s → s − 1, uniform on 0..6. If < 4, emit 2 bits; otherwise discard (≈ 3/7 expected).
- Stream = draw by draw: main bits, then strong bits. Under H0 every emitted bit is an i.i.d. fair bit.
- p-value per test on the whole stream. Proportion over the 8 blocks of 41 draws (bits from each block).
- Same NIST implementation and the same applicability rules as v1 (scripts/nist_sp800_22.py, α = 0.01).

### Status rules (same as v1)
- chi² / KS / NIST: p > 0.01 required (p > 0.05 shown).
- Runs test: p > 0.05 required.
- Proportion: inside the NIST interval 0.99 ± 3·√(0.99·0.01/m).
- Entropy ≥ 99.9% of maximum.
- Correlogram: every lag 1..20 inside ±1.96/√N.
- Battery: 15/15 PASS.
