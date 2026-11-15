# -*- coding: utf-8 -*-
"""
Dynamics of the 328-row 6/37 table: transition rules, transition function f,
equations over the field F_37, cycle structure of f.
Exact arithmetic only. No simulation, no random numbers, no seeds, no resampling.
"""
import csv, math, hashlib, itertools
from fractions import Fraction
import numpy as np
from scipy.stats import binom, chi2, norm

OUT = []
def P(s=""):
    print(s); OUT.append(str(s))

C = math.comb
p_ = 37                      # the field characteristic
NMAIN, K = 37, 6
OMEGA = C(37, 6)             # 2,324,784
OMEGA_FULL = OMEGA * 7

# ---------------------------------------------------------------- S0: data
path = "/home/user/Moe/v2/data/merged_lotto_table.csv"
raw = open(path, "rb").read()
rows = list(csv.DictReader(open(path, encoding="utf-8")))
cols = ["N1_smallest (Q)", "N2 (W)", "N3 (E)", "N4 (R)", "N5 (T)", "N6_largest (Y)"]
S = np.array([[int(r[c]) for c in cols] for r in rows], dtype=np.int64)   # (328,6) sorted
AB = np.array([int(r["Strong (ABC)"]) for r in rows], dtype=np.int64)
N = len(S)
dates = [l.strip() for l in open("/home/user/Moe/data/raw/Date.csv", encoding="utf-8").read().splitlines()[1:1+N]]

P("="*100)
P("S0. DATA AND ORIENTATION")
P("="*100)
P(f"file                : {path}")
P(f"sha256              : {hashlib.sha256(raw).hexdigest()}")
P(f"rows N              : {N}   columns: 6 main (sorted) + 1 strong")
P(f"strictly increasing in every row (Q<W<E<R<T<Y) : {bool(np.all(np.diff(S,axis=1)>0))}")
P(f"main range          : [{S.min()},{S.max()}]   strong range: [{AB.min()},{AB.max()}]")
P(f"row 1 date          : {dates[0]}      row {N} date: {dates[-1]}")
P(f"time direction      : row index increases => date DECREASES (reverse chronological)")
P(f"  => the pair (row n -> row n+1) is a step BACKWARD in time. A forward-time")
P(f"     transition function is (row n+1 -> row n), i.e. the inverse relation.")
P("     Every test below is reported for the row order as stored; for a bijective f")
P("     the two directions are related by f -> f^-1 (same cycle type), for a")
P("     non-injective f they are different objects.")

# ------------------------------------------------- S1: state space, combinatorics
P(""); P("="*100)
P("S1. STATE SPACE (combinatorics / discrete mathematics)")
P("="*100)
states = [tuple(r) for r in S.tolist()]
states_full = [tuple(r)+(int(a),) for r,a in zip(S.tolist(), AB.tolist())]
P(f"|Omega| main = C(37,6)            = {OMEGA:,}")
P(f"|Omega| main x strong             = {OMEGA_FULL:,}")
P(f"observed states                   = {N}  ->  coverage {N/OMEGA:.3e} of Omega")
P(f"distinct main states              = {len(set(states))}  (duplicates: {N-len(set(states))})")
P(f"distinct full states              = {len(set(states_full))}  (duplicates: {N-len(set(states_full))})")
exp_coll = Fraction(C(N,2), OMEGA)
P(f"expected coincident pairs under H0 (exact) = C({N},2)/|Omega| = {exp_coll} = {float(exp_coll):.4f}")
P(f"  => all states distinct is the EXPECTED outcome, not an anomaly, and it is also")
P(f"     exactly what makes a deterministic rule unfalsifiable at the state level (S2).")

# ------------------------------------- S2: is a transition function well defined?
P(""); P("="*100)
P("S2. WELL-DEFINEDNESS OF A TRANSITION FUNCTION (definition check)")
P("="*100)
P("Definition used: f is a transition function for the data iff for all m,n:")
P("  pi(state_m) = pi(state_n)  =>  successor(m) = successor(n),   where pi is whatever")
P("  the rule reads from the state. A function cannot send one input to two outputs.")
P("")
T = N - 1   # 327 transitions
P(f"number of transitions T = {T}")
P("")
P("(a) full state as input")
succ = {}
conflict = 0
for n in range(T):
    k = states_full[n]
    if k in succ and succ[k] != states_full[n+1]: conflict += 1
    succ[k] = states_full[n+1]
P(f"    repeated inputs with different outputs: {conflict}")
P(f"    => the relation state_n -> state_(n+1) IS a function on the 328 observed states,")
P(f"       for the trivial reason that no input repeats. Any assignment of successors to")
P(f"       328 distinct points is a function, so EXISTENCE of f carries zero information.")
P(f"       Validity verdict: the claim 'a transition function f exists' is true and VACUOUS")
P(f"       (it is interpolation on 328 points, not a dynamical law).")
P("")
P("(b) rules driven by a projection pi of the state (autonomous rules)")
P("    For each pi: #distinct values, #values whose successor is not unique,")
P("    minimum number of equations any such rule must violate.")
hdr = f"    {'pi(state)':<34}{'distinct':>9}{'conflicting':>12}{'min violations':>15}  verdict"
P(hdr)
def test_pi(name, keys, outs):
    d = {}
    for k,o in zip(keys,outs): d.setdefault(k,[]).append(o)
    distinct = len(d)
    confl = sum(1 for v in d.values() if len(set(v))>1)
    # minimum violated equations = sum over groups (size - largest identical block)
    minv = sum(len(v)-max(v.count(u) for u in set(v)) for v in d.values())
    verdict = "FALSIFIED" if confl else "not falsified"
    P(f"    {name:<34}{distinct:>9}{confl:>12}{minv:>15}  {verdict}")
    return confl, minv
names = ["Q=N1","W=N2","E=N3","R=N4","T=N5","Y=N6"]
for j in range(6):
    test_pi(f"N{j+1} ({names[j]})", [int(S[n,j]) for n in range(T)], [int(S[n+1,j]) for n in range(T)])
test_pi("strong ABC", [int(AB[n]) for n in range(T)], [int(AB[n+1]) for n in range(T)])
rank_r = (S.sum(axis=1) - 21).astype(int)          # Gale/Young rank, see S6
test_pi("rank r = sum-21", [int(rank_r[n]) for n in range(T)], [int(rank_r[n+1]) for n in range(T)])
test_pi("sum mod 37", [int(S[n].sum()%37) for n in range(T)], [tuple(S[n+1]) for n in range(T)])
test_pi("(min,max)=(Q,Y)", [(int(S[n,0]),int(S[n,5])) for n in range(T)], [tuple(S[n+1]) for n in range(T)])
test_pi("#odd numbers in the draw", [int((S[n]%2).sum()) for n in range(T)], [tuple(S[n+1]) for n in range(T)])
test_pi("multiset of residues mod 7", [tuple(sorted(S[n]%7)) for n in range(T)], [tuple(S[n+1]) for n in range(T)])
P("")
P("    Reading: every rule of the form 'next value = g(one coordinate)' , 'next state =")
P("    g(rank)' , 'next state = g(sum mod 37)' etc. is FALSIFIED by the data, for every g,")
P("    of any degree and any complexity. This is a theorem about the data, not a test with")
P("    a p-value: one input with two outputs refutes the whole family.")

# ---------------------------------------- S3: equations over the field F_37
P(""); P("="*100)
P("S3. EQUATIONS OVER THE FIELD (algebra)")
P("="*100)
def is_prime(m): return m>1 and all(m%i for i in range(2,int(m**.5)+1))
P(f"37 prime: {is_prime(37)}  =>  Z/37 = F_37 is a field (order 37, characteristic 37).")
P("Encoding used: v in {1..37} -> v mod 37 in F_37 (37 -> 0). This is a bijection, so no")
P("information is lost and no value is privileged.")
P("Terminology check on the strong number: ABC in {1..7} lives in F_7, a different field.")
P(f"The joint state would live in Z/37 x Z/7 = Z/259, and 259 = 7*37 is not a prime power,")
P("so NO field of order 259 exists: Z/259 is a ring with zero divisors. 'The equations over")
P("the field' can therefore only mean F_37 for the 6 main numbers (or F_7 for ABC alone).")
P("")
X = (S[:-1] % 37).astype(np.int64)      # inputs  (327,6)
Y = (S[1:]  % 37).astype(np.int64)      # outputs (327,6)

def gf_rank(A, mod=37):
    A = (A % mod).astype(np.int64).copy(); r = 0; rowsn, colsn = A.shape
    for c in range(colsn):
        piv = None
        for i in range(r, rowsn):
            if A[i, c] % mod: piv = i; break
        if piv is None: continue
        A[[r, piv]] = A[[piv, r]]
        inv = pow(int(A[r, c]), mod-2, mod)
        A[r] = (A[r] * inv) % mod
        nz = np.nonzero(A[r+1:, c])[0]
        if len(nz): A[r+1+nz] = (A[r+1+nz] - np.outer(A[r+1+nz, c], A[r])) % mod
        r += 1
        if r == rowsn: break
    return r

def monomials(d, nv=6):
    out = [()]
    for deg in range(1, d+1):
        out += list(itertools.combinations_with_replacement(range(nv), deg))
    return out

P("Model family: f_j(x) = polynomial over F_37 in the 6 input coordinates, total degree <= d,")
P("one polynomial per output coordinate j. Equations: 327 per coordinate (1962 in total).")
P("A consistent solution exists iff rank(M) = rank([M|y]).")
P("")
P(f"    {'d':>2}{'unknowns/coord':>16}{'equations/coord':>17}{'rank M':>8}  consistency per output coordinate j=1..6")
cons_table = []
for d in [1,2,3,4]:
    mons = monomials(d)
    M = np.ones((T, len(mons)), dtype=np.int64)
    for c, mon in enumerate(mons):
        col = np.ones(T, dtype=np.int64)
        for v in mon: col = (col * X[:, v]) % 37
        M[:, c] = col
    rM = gf_rank(M)
    verd = []
    for j in range(6):
        Aug = np.concatenate([M, Y[:, [j]]], axis=1)
        verd.append("consistent" if gf_rank(Aug) == rM else "INCONSISTENT")
    cons_table.append((d, len(mons), rM, verd))
    P(f"    {d:>2}{len(mons):>16}{T:>17}{rM:>8}  " + " ".join(f"j{j+1}:{v[:4]}" for j,v in enumerate(verd)))
P("")
P("    d=1 is the affine model x -> Ax+b : 6*6+6 = 42 unknowns against 1962 equations.")
P("    All degrees 1..4 are overdetermined (unknowns < 327) and ALL are inconsistent:")
P("    no polynomial transition function of total degree <= 4 over F_37 reproduces the data.")
P("")
P("    Identifiability threshold (why d<=4 is the whole testable range):")
P(f"    {'d':>3}{'unknowns C(d+6,6)':>20}{'equations':>11}{'status':>18}{'parameter bits':>16}")
for d in range(1,8):
    u = C(d+6,6)
    st = "overdetermined" if u < T else ("exactly determined" if u==T else "UNDERDETERMINED")
    bits = 6*u*math.log2(37)
    P(f"    {d:>3}{u:>20}{T:>11}{st:>18}{bits:>16.0f}")
data_bits = T*math.log2(OMEGA)
P(f"    information in the 327 transitions = 327*log2(C(37,6)) = {data_bits:.0f} bits")
P(f"    => from d=5 the model has more free parameters than equations, so a solution exists")
P(f"       for counting reasons alone, and from d=4 the parameters already cost about as many")
P(f"       bits as the data contains ({6*C(10,6)*math.log2(37):.0f} vs {data_bits:.0f}): no compression, no law.")
P("")
P("    Vacuity theorem (Lagrange interpolation over a finite field):")
P("    every map F_37^6 -> F_37^6 is given by a unique reduced polynomial (all exponents <= 36).")
P(f"    reduced coefficient count per output coordinate = 37^6 = {37**6:,}")
P("    So 'there exists a polynomial f over F_37 with f(state_n)=state_(n+1)' is TRUE for ANY")
P("    table whose inputs are distinct. It is a property of finite fields, not a finding about")
P("    this data. Only a degree (or sparsity) bound makes the statement falsifiable - and under")
P("    every falsifiable bound (d<=4) the answer here is NO.")
P("")
P("    Exhaustive search over the diagonal subfamily  N_j(n+1) = a*N_j(n) + b  (mod 37),")
P("    all 37*37 = 1369 pairs (a,b) per coordinate, 327 equations each.")
P("    The null probability of a single match is NOT 1/37: N_j is an order statistic, so its")
P("    marginal pmf is p_j(v) = C(v-1,j-1)C(37-v,6-j)/C(37,6). For a rule g the exact match")
P("    probability under H0 (independent draws) is q = sum_v p_j(v) * p_j(g(v)).")
P("    Consecutive transitions share a draw, so the 327 indicators are dependent; the exact")
P("    independent test uses the 164 disjoint transitions (rows 1-2, 3-4, ...).")
P("")
pmf = {}
for j in range(1,7):
    pmf[j] = [0.0]*38
    for v in range(1,38):
        pmf[j][v] = C(v-1,j-1)*C(37-v,6-j)/OMEGA
def dec(u): return 37 if u % 37 == 0 else u % 37
disj = list(range(0, T, 2))
P(f"    {'coord':<6}{'family':<14}{'best (a,b)':>11}{'match/327':>11}{'q_exact':>9}{'E[match]':>9}"
  f"{'match/164 disj':>15}{'exact binom p':>15}{'verdict':>14}")
alpha_bonf = 0.05/(1332*6)
diag_rows = []
for j in range(6):
    for fam, avals in (("affine a!=0", range(1,37)), ("constant a=0", [0])):
        best = None
        for a in avals:
            for b in range(37):
                img = (a*X[:,j] + b) % 37
                m = int(np.count_nonzero(img == Y[:,j]))
                q = sum(pmf[j+1][v]*pmf[j+1][dec(a*(v%37)+b)] for v in range(1,38))
                md = int(np.count_nonzero(img[disj] == Y[disj,j]))
                pv = float(binom.sf(md-1, len(disj), q)) if q>0 else 1.0
                if best is None or pv < best[0]: best = (pv, (a,b), m, q, md)
        pv, ab, m, q, md = best
        sig = "significant" if (fam.startswith("affine") and pv < alpha_bonf) else "no"
        P(f"    N{j+1:<5}{fam:<14}{str(ab):>11}{m:>11}{q:>9.4f}{T*q:>9.1f}{md:>15}{pv:>15.4g}{sig:>14}")
        diag_rows.append((j+1, fam, ab, m, q, md, pv))
P(f"    Bonferroni threshold for the 1332*6 = 7992 non-degenerate affine rules at 0.05 : {alpha_bonf:.3g}")
P("    The a=0 row is the degenerate member 'next value = constant b'. It reads nothing from the")
P("    past, so it is not a transition rule; its best b is just the modal value of column j")
P("    (e.g. Y=37 in 50 rows, Q=1 in 67 rows) and that is a marginal frequency property,")
P("    already covered by the uniformity block of the existing report, not dynamics.")

# ------------------------------------------------- S4: cycle structure of f
P(""); P("="*100)
P("S4. CYCLE STRUCTURE OF f (dynamical systems)")
P("="*100)
P("Definitions: x is periodic iff f^k(x)=x for some k>=1; the cycle structure of a map on a")
P("finite set is the multiset of cycle lengths of its functional graph; the data give one")
P("ORBIT SEGMENT, not the map.")
P("")
P(f"observed orbit segment length            : {N} states, {T} edges")
P(f"repeated states (=> cycles observed)     : {N-len(set(states_full))}")
P(f"fixed points observed (f(x)=x)           : {sum(1 for n in range(T) if states_full[n]==states_full[n+1])}")
P(f"periodic points observed                 : 0")
P(f"functional graph of the observed partial map: one directed simple PATH on {N} vertices")
P(f"  -> connected components 1, cycles 0, cycle lengths {{}} (empty multiset),")
P(f"     transient length >= {N}, rho-shape (Floyd) not reached.")
P("")
nrem = OMEGA_FULL - T
log10_perm = float((math.lgamma(nrem+1))/math.log(10))
log10_fun  = float((OMEGA_FULL - T)*math.log10(OMEGA_FULL))
P("Identifiability of the cycle structure (exact counts of consistent extensions):")
P(f"  maps f: Omega->Omega agreeing with the observed path      = |Omega|^(|Omega|-{T}) , log10 = {log10_fun:.3e}")
P(f"  permutations of Omega containing the path (contract it)    = (|Omega|-{T})! , log10 = {log10_perm:.3e}")
P("  Every cycle type of Omega that contains at least one cycle of length >= 328 is realised by")
P("  at least one such extension, and so is every cycle type reachable by a non-injective f.")
P(f"  => the cycle structure of f is NOT identified by the data. Any statement of the form")
P(f"     'f has a cycle of length L' is unconstrained by these {N} rows.")
P("")
P("Even the best case (f linear and invertible over F_37) is out of reach:")
glo = 1
for i in range(6): glo *= (37**6 - 37**i)
P(f"  |GL_6(F_37)| = {glo:,}")
P(f"  maximal order of an element (Singer cycle) = 37^6 - 1 = {37**6-1:,}")
P(f"  factorisation of 37^6-1 : ", )
m = 37**6-1; fac = {}; dd = 2
while dd*dd <= m:
    while m % dd == 0: fac[dd] = fac.get(dd,0)+1; m//=dd
    dd += 1
if m>1: fac[m] = fac.get(m,0)+1
P("    " + " * ".join(f"{k}^{v}" if v>1 else str(k) for k,v in sorted(fac.items())))
P(f"  the orbit of a point under such an f has length dividing 37^6-1, i.e. up to ~2.6e9 steps,")
P(f"  so {N} consecutive states can never exhibit the period. Lower bound from the data on the")
P(f"  period of any deterministic f consistent with it: period > {N} (or f is not injective).")
P("")
P("Terminology correction for the findings already in the repository:")
P("  * 'periodicity of 6' (lag-6 autocorrelation 0.7285) is a period of the REPRESENTATION")
P("    (6 sorted numbers written in a row), not a period of a map: it is not an f^k(x)=x")
P("    statement. Quantified exactly in S5.")
P("  * 'J=8 cycles' in NIST Random Excursions are zero-crossing cycles of a random walk on the")
P("    bit stream, a different object again.")
P("  * 'pass ratio / yachas maavar' in the NIST tables is NIST's PROPORTION of passing")
P("    subsequences. It is not a transition rule and not a transition matrix.")

# ------------------- S5: the "periodicity of 6" quantified exactly under H0
P(""); P("="*100)
P("S5. THE LAG-6 'PERIODICITY' - EXACT PREDICTION UNDER H0 (artifact of the representation)")
P("="*100)
mu = [None]+[Fraction(38*i,7) for i in range(1,7)]
pmf = {j: [Fraction(0)]+[Fraction(C(v-1,j-1)*C(37-v,6-j), OMEGA) for v in range(1,38)] for j in range(1,7)}
for j in range(1,7):
    m = sum(v*pmf[j][v] for v in range(1,38))
    assert m == mu[j], (j, m, mu[j])
EN2 = {j: sum(v*v*pmf[j][v] for v in range(1,38)) for j in range(1,7)}
def joint_E(i, j):                      # E[N_i N_j], i<j, exact
    s = Fraction(0)
    for a in range(1,38):
        if C(a-1,i-1)==0: continue
        for b in range(a+1,38):
            w = C(a-1,i-1)*C(b-a-1,j-i-1)*C(37-b,6-j)
            if w: s += Fraction(w*a*b, OMEGA)
    return s
EIJ = {}
for i in range(1,7):
    for j in range(i+1,7): EIJ[(i,j)] = joint_E(i,j)
def Ecross(i, j):                       # E[(N_i-19)(N_j-19)] within one draw, exact
    if i==j: return EN2[i] - 38*mu[i] + 361
    a,b = min(i,j), max(i,j)
    return EIJ[(a,b)] - 19*mu[a] - 19*mu[b] + 361
tot_var = sum(Ecross(i,i) for i in range(1,7))/6
P(f"exact E[N_i], i=1..6 (= 38i/7)      : " + ", ".join(f"{float(mu[i]):.4f}" for i in range(1,7)))
P(f"pooled mean (exact)                 : {float(sum(mu[i] for i in range(1,7))/6)}  (= 19)")
P(f"pooled variance (exact)             : {float(tot_var)}  (= (37^2-1)/12 = 114, the marginal of a")
P(f"                                      random element of the draw is uniform on 1..37)")
P(f"between-position component           : {float(sum((mu[i]-19)**2 for i in range(1,7))/6):.4f}")
z = S.reshape(-1).astype(float)                       # the sorted stream, 1968 values
zc = z - z.mean(); den = float((zc*zc).sum())
P("")
P("Theoretical autocorrelation of the concatenated SORTED stream under H0 (independent draws,")
P("global-mean estimator) against the observed one. Positions i'=((i-1+L) mod 6)+1 :")
P(f"    {'lag L':>6}{'theoretical r_L':>18}{'observed r_L':>15}{'difference':>13}")
acf_rows=[]
for L in range(1,21):
    th = sum((Ecross(i, ((i-1+L)%6)+1) if (i-1+L) < 6 else (mu[i]-19)*(mu[((i-1+L)%6)+1]-19)) for i in range(1,7))/6/tot_var
    ob = float((zc[:-L]*zc[L:]).sum()/den)
    acf_rows.append((L,float(th),ob))
    P(f"    {L:>6}{float(th):>18.4f}{ob:>15.4f}{ob-float(th):>13.4f}")
P(f"    max |observed - theoretical| over lags 1..20 = {max(abs(o-t) for _,t,o in acf_rows):.4f}")
P("")
P("Verdict: the whole correlogram of the sorted stream, including r_6 = 0.7285 / r_12 / r_18")
P("reported as 'periodicity of 6' and 'all 20 lags outside the band', is produced EXACTLY by")
P("the fact that the 6 numbers are written in increasing order. It is a property of the")
P("representation (order statistics have different means per position), is fully predicted")
P("under the hypothesis of independent fair draws, and is NOT a period of any map f.")
P("The correct statement: 'the sorted stream is not stationary, so the autocorrelation")
P("estimator with a single global mean does not test independence here'.")

# --------------------------------------------- S6: lattice structure
P(""); P("="*100)
P("S6. LATTICE STRUCTURE (sriggim) AND THE TRANSITION STATISTICS THAT ARE VALID")
P("="*100)
P("Correct frame: Omega = {x in Z^6 : 1<=x_1<x_2<...<x_6<=37} with the componentwise (Gale)")
P("order. Meet and join are componentwise min and max, which preserve strict increase, so")
P("(Omega,<=) is a LATTICE, distributive, isomorphic to L(6,31) = partitions inside a 6x31 box")
P("(substitute y_i = x_i - i, then 0<=y_1<=...<=y_6<=31), i.e. an interval of Young's lattice.")
P("Rank function: r(x) = sum(x) - 21, graded from 0 to 186.")
ok = True
for n in range(0, N, 7):
    for m2 in range(0, N, 11):
        a, b = S[n], S[m2]
        if not (np.all(np.diff(np.maximum(a,b))>0) and np.all(np.diff(np.minimum(a,b))>0)): ok = False
P(f"join/meet closure verified on sampled observed pairs: {ok}")
P("NOT a sublattice of the Boolean lattice 2^[37]: the 6-element subsets form the rank-6 LEVEL,")
P("which is an antichain under inclusion and is not closed under union or intersection.")
P("So any claim phrased as 'lattice of the draws' must mean the Gale/Young lattice above.")
P("")
# exact count of ordered pairs x <= y (Gale), DP over positions; validated by brute force
def comparable_pairs(nn, kk):
    dp = [[1 if b>=a else 0 for b in range(nn+2)] for a in range(nn+2)]
    for _ in range(kk-1):
        suf = [[0]*(nn+3) for _ in range(nn+3)]
        for a in range(nn,0,-1):
            for b in range(nn,0,-1):
                suf[a][b] = dp[a][b] + suf[a+1][b] + suf[a][b+1] - suf[a+1][b+1]
        dp = [[(suf[a+1][b+1] if b>=a else 0) for b in range(nn+2)] for a in range(nn+2)]
    return sum(dp[a][b] for a in range(1,nn+1) for b in range(a,nn+1))
def brute(nn, kk):
    allx = list(itertools.combinations(range(1,nn+1), kk))
    return sum(1 for x in allx for y in allx if all(u<=v for u,v in zip(x,y)))
P("validation of the comparability DP against brute force:")
for nn,kk in [(4,2),(5,2),(6,3),(7,3),(8,4)]:
    P(f"  n={nn} k={kk}: DP {comparable_pairs(nn,kk)}  brute force {brute(nn,kk)}  match {comparable_pairs(nn,kk)==brute(nn,kk)}")
comp = comparable_pairs(37,6)
Pleq = Fraction(comp, OMEGA*OMEGA)
Pcomp = 2*Pleq - Fraction(1, OMEGA)
P(f"exact number of ordered pairs (x,y) in Omega^2 with x <= y (Gale) : {comp:,}")
P(f"P(x <= y) exact = {float(Pleq):.6e}   P(x,y comparable) = {float(Pcomp):.6e}")
le = int(sum(1 for n in range(T) if np.all(S[n] <= S[n+1])))
ge = int(sum(1 for n in range(T) if np.all(S[n] >= S[n+1])))
P(f"observed comparable consecutive pairs: {le} ascending + {ge} descending = {le+ge} out of {T}")
P(f"expected under H0 = {T}*P(comparable) = {T*float(Pcomp):.4f}")
rr = (S.sum(axis=1) - 21).astype(float)
pc = float(Pcomp)
sd_c = math.sqrt(T*pc*(1-pc))
P(f"approximate z on all {T} (overlapping, hence dependent) transitions : "
  f"({le+ge} - {T*pc:.2f})/{sd_c:.2f} = {((le+ge)-T*pc)/sd_c:+.3f}")
dj = list(range(0, T, 2))
cd = int(sum(1 for n in dj if np.all(S[n] <= S[n+1]) or np.all(S[n] >= S[n+1])))
p2 = float(2*min(binom.sf(cd-1, len(dj), pc), binom.cdf(cd, len(dj), pc)))
P(f"EXACT test on the {len(dj)} disjoint (independent) transitions: {cd} comparable, "
  f"expected {len(dj)*pc:.2f}, two-sided exact binomial p = {p2:.4f}")
dj2 = list(range(1, T, 2))
cd2 = int(sum(1 for n in dj2 if np.all(S[n] <= S[n+1]) or np.all(S[n] >= S[n+1])))
p2b = float(2*min(binom.sf(cd2-1, len(dj2), pc), binom.cdf(cd2, len(dj2), pc)))
P(f"EXACT test on the other {len(dj2)} disjoint transitions (odd index): {cd2} comparable, "
  f"expected {len(dj2)*pc:.2f}, two-sided exact binomial p = {p2b:.4f}")
# exact variance of the overlapping count: needs P(I_n=1 and I_(n+1)=1) over 3 independent draws
def triple_count(rel1, rel2, nn=37, kk=6):
    R = range(1, nn+1)
    def ok(a,b,c): return (a<=b if rel1 else b<=a) and (b<=c if rel2 else c<=b)
    dp = [[[1 if ok(a,b,c) else 0 for c in range(nn+2)] for b in range(nn+2)] for a in range(nn+2)]
    for _ in range(kk-1):
        suf = [[[0]*(nn+3) for _ in range(nn+3)] for _ in range(nn+3)]
        for a in range(nn,0,-1):
            for b in range(nn,0,-1):
                for c in range(nn,0,-1):
                    suf[a][b][c] = (dp[a][b][c] + suf[a+1][b][c] + suf[a][b+1][c] + suf[a][b][c+1]
                                    - suf[a+1][b+1][c] - suf[a+1][b][c+1] - suf[a][b+1][c+1]
                                    + suf[a+1][b+1][c+1])
        dp = [[[(suf[a+1][b+1][c+1] if ok(a,b,c) else 0) for c in range(nn+2)] for b in range(nn+2)] for a in range(nn+2)]
    return sum(dp[a][b][c] for a in R for b in R for c in R)
tc = {(r1,r2): triple_count(r1,r2) for r1 in (True,False) for r2 in (True,False)}
P11 = Fraction(sum(tc.values()), OMEGA**3)
P(f"exact P(two consecutive transitions both comparable) = {float(P11):.6f}  (vs p^2 = {pc*pc:.6f})")
var_exact = T*float(Pcomp)*(1-float(Pcomp)) + 2*(T-1)*(float(P11) - float(Pcomp)**2)
P(f"EXACT variance of the count over the {T} overlapping transitions = {var_exact:.4f} "
  f"(sd {math.sqrt(var_exact):.3f}; the independence value would be {T*pc*(1-pc):.4f})")
zc_ex = ((le+ge) - T*pc)/math.sqrt(var_exact)
P(f"  => exact-variance z = {zc_ex:+.3f}, two-sided normal p = {float(2*norm.sf(abs(zc_ex))):.4f}")
P("  => a mild excess of order-comparable consecutive pairs, not significant on the exact")
P("     independent test, and it is one statistic among the many computed here, so no")
P("     multiplicity-corrected claim can be made from it. It is NOT evidence of a monotone")
P("     (order-preserving) transition rule: such a rule would force comparability in 327/327.")
P("")
P("Exact distribution of the rank r = sum - 21 under H0 (DP over subsets, no approximation):")
cnt = [0]*(187)
dpr = {(0,0):1}
for v in range(1,38):
    nd = dict(dpr)
    for (c,s),w in dpr.items():
        if c < 6: nd[(c+1,s+v)] = nd.get((c+1,s+v),0) + w
    dpr = nd
for (c,s),w in dpr.items():
    if c == 6: cnt[s-21] += w
assert sum(cnt) == OMEGA
pr = [c/OMEGA for c in cnt]
mean_ex = sum(i*pr[i] for i in range(187)); var_ex = sum((i-mean_ex)**2*pr[i] for i in range(187))
m4_ex = sum((i-mean_ex)**4*pr[i] for i in range(187))
P(f"  exact E[r] = {mean_ex:.6f}   exact Var(r) = {var_ex:.6f}   exact 4th central moment = {m4_ex:.1f}")
s2 = float(rr.var(ddof=1))
sd_s2 = math.sqrt((m4_ex - var_ex**2)/N)
P(f"  observed sample variance = {s2:.2f}; sd of the sample variance under H0 = {sd_s2:.2f}"
  f"  -> z = {(s2-var_ex)/sd_s2:+.3f}, two-sided p = {float(2*norm.sf(abs((s2-var_ex)/sd_s2))):.4f}")
# equiprobable binning from the exact cdf: every expected count ~ N/10
cdf = np.cumsum(pr)
qs = [0.1*t for t in range(1,10)]
edges = sorted(set(int(np.searchsorted(cdf, q)) for q in qs))
bins = [0] + [e+1 for e in edges] + [187]
O, E2 = [], []
for t in range(len(bins)-1):
    O.append(int(((rr>=bins[t]) & (rr<bins[t+1])).sum()))
    E2.append(N*sum(pr[bins[t]:bins[t+1]]))
x2r = sum((o-e)**2/e for o,e in zip(O,E2)); dfr = len(O)-1
P(f"  equiprobable-cell goodness of fit (all expected >= {min(E2):.1f}): chi-square = {x2r:.4f}, "
  f"df = {dfr}, p = {float(chi2.sf(x2r,dfr)):.4f}")
P(f"    cell edges r = {bins}")
P(f"    observed/expected: " + ", ".join(f"{o}/{e:.1f}" for o,e in zip(O,E2)))
lo = int((rr <= 55).sum()); hi = int((rr >= 131).sum())
P(f"    tails: r<=55 observed {lo}, expected {N*sum(pr[:56]):.2f} ; "
  f"r>=131 observed {hi}, expected {N*sum(pr[131:]):.2f}")
P("  The rank variance excess, the equiprobable-cell fit and the comparability excess above are")
P("  NOT independent findings: over-dispersed sums put more draws near the top and bottom of the")
P("  lattice, and extreme draws are Gale-comparable to many more states. One phenomenon, counted")
P("  three times. About 25 statistics are computed in this file, so at the 0.05 level ~1.25")
P("  nominal rejections are expected by chance; nothing here survives a multiplicity correction.")
P("")
P("Exploratory (NOT a confirmatory test - this slice was chosen after seeing the split above):")
day = {}
try:
    import csv as _c
    rows8 = list(_c.DictReader(open("/home/user/Moe/independent_agent_2/full_table_8_columns.csv", encoding="utf-8")))
    days = [r["Day"] for r in rows8][:N]
    pairs = {}
    for n in range(T):
        comp = bool(np.all(S[n] <= S[n+1]) or np.all(S[n] >= S[n+1]))
        k = (days[n], days[n+1])
        pairs.setdefault(k, [0,0])
        pairs[k][0] += 1; pairs[k][1] += int(comp)
    P(f"    {'weekday(n) -> weekday(n+1)':<34}{'transitions':>12}{'comparable':>12}{'expected':>10}")
    for k,(tot,c) in sorted(pairs.items(), key=lambda kv:-kv[1][0]):
        P(f"    {k[0]+' -> '+k[1]:<34}{tot:>12}{c:>12}{tot*pc:>10.2f}")
    P("    Reading: the mild excess is spread across all the frequent weekday pairs (55/43.8,")
    P("    35/30.2, 17/13.6, 16/14.9), it is not concentrated in one direction of the weekly cycle.")
    P("    So the even/odd discrepancy is noise inside a uniform mild excess, and it is consistent")
    P("    with the rank over-dispersion above rather than with any direction-dependent rule.")
    P("    To test the excess at all it must be recomputed on draws NOT in this table.")
except Exception as e:
    P(f"    (weekday column unavailable: {e})")
P("")
P("Rank (height) dynamics on the lattice - exact null moments:")
rr = (S.sum(axis=1) - 21).astype(float)
P(f"  E[sum] = 6*38/2 = 114 -> E[r] = 93 exactly ; Var(sum) = k(n+1)(n-k)/12 = 6*38*31/12 = 589")
P(f"  observed mean r = {rr.mean():.4f} (null 93, sd of the mean {math.sqrt(589/N):.4f}, z = {(rr.mean()-93)/math.sqrt(589/N):+.3f})")
P(f"  observed var  r = {rr.var(ddof=1):.2f} (null 589, exact test below)")
d1 = rr[1:]-rr[:-1]
P(f"  increments r(n+1)-r(n): mean {d1.mean():+.4f} (null 0), sd {d1.std(ddof=1):.2f} (null sqrt(2*589)={math.sqrt(2*589):.2f})")
r1 = float(np.corrcoef(rr[:-1], rr[1:])[0,1])
P(f"  lag-1 autocorrelation of r(n) = {r1:+.4f}, 95% band +-{1.96/math.sqrt(N):.4f} -> {'inside' if abs(r1)<1.96/math.sqrt(N) else 'OUTSIDE'}")
P("  (this series IS stationary under H0, unlike the sorted stream of S5, so here the")
P("   autocorrelation estimator is the right tool.)")
P("")
P("Overlap of consecutive draws |S_n and S_(n+1)| - exact hypergeometric null:")
ov = np.array([len(set(S[n].tolist()) & set(S[n+1].tolist())) for n in range(T)])
P(f"    {'k':>3}{'observed':>10}{'P(k) exact':>14}{'expected':>10}")
obs_c = [int((ov==k).sum()) for k in range(7)]
pk = [C(6,k)*C(31,6-k)/OMEGA for k in range(7)]
for k in range(7): P(f"    {k:>3}{obs_c[k]:>10}{pk[k]:>14.6f}{T*pk[k]:>10.2f}")
P(f"  exact E|S_n and S_(n+1)| = 36/37 = {36/37:.6f} ; observed mean = {ov.mean():.6f}")
cells = [(obs_c[0],pk[0]),(obs_c[1],pk[1]),(obs_c[2],pk[2]),(sum(obs_c[3:]),sum(pk[3:]))]
x2 = sum((o-T*q)**2/(T*q) for o,q in cells); dfree = len(cells)-1
P(f"  chi-square (cells 0,1,2,>=3) = {x2:.4f}, df = {dfree}, p = {float(chi2.sf(x2,dfree)):.4f}")
P("  => the carry-over of numbers between consecutive draws matches independent draws.")
P("")
P("Strong number ABC as a Markov chain on 7 states (the one genuine transition matrix here):")
Tm = np.zeros((7,7), dtype=int)
for n in range(T): Tm[AB[n]-1, AB[n+1]-1] += 1
P("    rows = ABC(n), cols = ABC(n+1):")
for i in range(7): P("      " + " ".join(f"{Tm[i,j]:3d}" for j in range(7)) + f"   row sum {Tm[i].sum():3d}")
exp_cell = T/49
x2t = float(((Tm-exp_cell)**2/exp_cell).sum())
P(f"    chi-square vs uniform 1/49 = {x2t:.4f}, df = 48, p = {float(chi2.sf(x2t,48)):.4f} (expected cell {exp_cell:.2f})")
P(f"    chi-square of independence (rows x cols) :")
rs = Tm.sum(1); cs = Tm.sum(0); E = np.outer(rs,cs)/T
mask = E>0
x2i = float((((Tm-E)**2/np.where(mask,E,1))[mask]).sum())
P(f"      X2 = {x2i:.4f}, df = 36, p = {float(chi2.sf(x2i,36)):.4f}  -> no dependence of ABC(n+1) on ABC(n)")

# --------------------------------------------- S7: geometry
P(""); P("="*100)
P("S7. GEOMETRY")
P("="*100)
P("Omega is the set of lattice points of the cone x_1<x_2<...<x_6 intersected with the box")
P("[1,37]^6; after y_i = x_i - i it is the set of lattice points of the order polytope / the")
P("6x31 box partitions, |Omega| = C(37,6). The rank r = sum(x)-21 is a LINEAR functional, so")
P("'height' dynamics (S6) is the projection of the walk on the diagonal direction.")
# exact moments of the L1 step under H0, including the dependence of consecutive steps
V37 = np.arange(1,38)
A = np.abs(V37[:,None]-V37[None,:]).astype(float)
p1 = {j: np.array([float(pmf[j][v]) for v in range(1,38)]) for j in range(1,7)}
def Pjoint(j,k):
    if j==k: return np.diag(p1[j])
    a,b = min(j,k), max(j,k)
    M = np.zeros((37,37))
    for u in range(1,38):
        if C(u-1,a-1)==0: continue
        for w in range(u+1,38):
            ww = C(u-1,a-1)*C(w-u-1,b-a-1)*C(37-w,6-b)
            if ww: M[u-1,w-1] = ww/OMEGA
    return M if j<k else M.T
PJ = {(j,k): Pjoint(j,k) for j in range(1,7) for k in range(1,7)}
EL1 = sum(float(p1[j] @ A @ p1[j]) for j in range(1,7))
ED2 = sum(float(np.sum(PJ[(j,k)] * (A @ PJ[(j,k)] @ A))) for j in range(1,7) for k in range(1,7))
varD = ED2 - EL1**2
h = {j: A @ p1[j] for j in range(1,7)}
Eg2 = sum(float(h[j] @ PJ[(j,k)] @ h[k]) for j in range(1,7) for k in range(1,7))
rho_ex = (Eg2 - EL1**2)/varD
st = np.abs(S[1:]-S[:-1]).sum(axis=1).astype(float)
ste = np.sqrt(((S[1:]-S[:-1]).astype(float)**2).sum(axis=1))
P(f"exact E[L1 step] between two independent draws = {EL1:.4f}   exact sd = {math.sqrt(varD):.4f}")
P(f"observed mean L1 step = {st.mean():.4f}   observed sd = {st.std(ddof=1):.4f}")
P("Consecutive steps share the middle draw, so they are POSITIVELY correlated under H0 and the")
P("+-1.96/sqrt(n) band does not apply. Exact null correlation of two adjacent steps:")
P(f"  Corr(|X-Y|_1, |Y-Z|_1) for X,Y,Z iid uniform on Omega = {rho_ex:+.4f}  (exact)")
sa = float(np.corrcoef(st[:-1], st[1:])[0,1])
P(f"  observed lag-1 autocorrelation of the step length = {sa:+.4f}")
P(f"  => {abs(sa-rho_ex):.4f} away from the exact H0 value; the naive band +-{1.96/math.sqrt(T-1):.4f}")
P(f"     would have called {sa:+.4f} an anomaly. It is not one: it is the exact H0 prediction.")
sem_c = math.sqrt(varD*(1+2*rho_ex)/T)
P(f"mean step, variance-corrected for the dependence: z = ({st.mean():.3f} - {EL1:.3f})/{sem_c:.3f} "
  f"= {(st.mean()-EL1)/sem_c:+.3f}, two-sided p = {float(2*norm.sf(abs((st.mean()-EL1)/sem_c))):.4f}")
P(f"observed mean L2 step = {ste.mean():.4f}")
P("No geometric drift, contraction or attractor: the step-length distribution is what independent")
P("sampling from Omega gives, once the shared-draw dependence is accounted for exactly.")

# --------------------------------------------- S8: verdict
P(""); P("="*100)
P("S8. VERDICT TABLE")
P("="*100)
V = [
 ("transition rules (autonomous, driven by a projection of the state)",
  "FALSIFIED. Every rule next = g(N_j), g(rank), g(sum mod 37), g(#odd), g(residues mod 7), g(min,max) "
  "is refuted: the same input occurs with different outputs (S2b), for every g."),
 ("transition function f on the full state",
  "EXISTS BUT VACUOUS. The 328 states are pairwise distinct (expected: 0.023 coincidences), so any "
  "successor assignment is a function. Existence is interpolation on 328 points, not a law."),
 ("equations over the field F_37, more equations than unknowns",
  "CORRECTLY SET UP ONLY FOR d<=4 (unknowns C(d+6,6) < 327), AND INCONSISTENT AT EVERY SUCH d: "
  "no polynomial f of total degree <=4 over F_37 fits (exact rank test). For d>=5 the system is "
  "underdetermined, so a solution exists by counting alone and proves nothing. Over a finite field "
  "EVERY map is polynomial (37^6 reduced coefficients), so unbounded-degree claims are vacuous."),
 ("best non-degenerate affine rule per coordinate",
  "NOT SIGNIFICANT. Exhaustive over all 1332 rules a!=0 per coordinate, exact order-statistic null, "
  "exact binomial on 164 independent transitions: smallest p = 4.5e-4 against the Bonferroni "
  "threshold 6.3e-6 (the minimum of ~1332 p-values is expected near 7.5e-4)."),
 ("cycle structure of f",
  "UNDEFINED / UNIDENTIFIABLE. The data are one orbit SEGMENT: a directed simple path, 0 cycles, "
  "0 fixed points, 0 periodic points. Consistent extensions: (|Omega|-327)! permutations "
  "(log10 ~ 1.1e8). Every cycle type with a cycle of length >= 328 is realisable."),
 ("'periodicity of 6' in the existing report",
  "ARTIFACT, QUANTIFIED. Predicted exactly by H0 for the sorted stream (theoretical r_6 = 0.7539 vs "
  "observed 0.7285); it is a period of the representation, not of a map."),
 ("'J = 8 cycles', 'pass ratio'",
  "DIFFERENT CONCEPTS. Random-walk zero-crossing cycles; NIST proportion of passing subsequences. "
  "Neither is a cycle of f nor a transition matrix."),
 ("lattice",
  "VALID FRAME: Gale order on Omega = distributive lattice L(6,31) (6x31 box partitions), rank "
  "sum-21 in [0,186]. The level of 6-subsets inside the Boolean lattice 2^[37] is an antichain, "
  "not a lattice. Exact P(two draws Gale-comparable) = 0.339285 (DP validated against brute "
  "force); observed 132/327 vs 110.9 expected."),
 ("valid dynamical observables, consistent with independent draws",
  "rank r: mean 94.25 vs exact 93 (z=+0.93), equiprobable-cell fit to the exact pmf p=0.50, lag-1 "
  "ACF -0.003 inside the band; overlap |S_n and S_(n+1)|: mean 1.046 vs exact 36/37, chi-square "
  "p=0.21; ABC 7x7 transition matrix: independence p=0.64, uniformity p=0.42; L1 step length "
  "36.68 vs the exact 35.04 (dependence-corrected z=+1.49, p=0.138), and its lag-1 ACF +0.175 "
  "against the EXACT H0 value +0.205, not against the inapplicable band +-0.109."),
 ("the only mild deviations found, and they are one phenomenon",
  "sample variance of the rank 678 vs exact 589 (z=+2.05, p=0.041) and Gale-comparable consecutive "
  "pairs 132 vs 110.9 with the EXACT dependent variance (z=+2.20, p=0.028). Over-dispersed sums "
  "mechanically raise comparability, so this is one effect counted twice; with ~25 statistics "
  "computed here (expected ~1.25 nominal rejections) neither survives any multiplicity correction. "
  "Neither is a transition rule: they are properties of the one-step marginal, not of f."),
 ("orientation of time",
  "The table is reverse chronological (row 1 = 06/10/2026, row 328 = 19/12/2023). Any f fitted as "
  "row n -> row n+1 runs backwards in time."),
]
for k, v in V:
    P(f"* {k}")
    for ln in [v[i:i+96] for i in range(0,len(v),96)]: P(f"    {ln}")
open("/home/user/Moe/v3_dynamics/results/machine_output.txt","w",encoding="utf-8").write("\n".join(OUT)+"\n")
import csv as _csv
with open("/home/user/Moe/v3_dynamics/results/acf_sorted_stream_theory_vs_observed.csv","w",newline="",encoding="utf-8") as fh:
    w=_csv.writer(fh); w.writerow(["lag","theoretical_r_under_H0","observed_r","difference"])
    for L,t,o in acf_rows: w.writerow([L,f"{t:.6f}",f"{o:.6f}",f"{o-t:.6f}"])
with open("/home/user/Moe/v3_dynamics/results/field_equations_consistency.csv","w",newline="",encoding="utf-8") as fh:
    w=_csv.writer(fh); w.writerow(["degree_d","unknowns_per_coordinate","equations_per_coordinate","rank_M"]+[f"coord_{j+1}" for j in range(6)])
    for d,u,r,verd in cons_table: w.writerow([d,u,T,r]+verd)
print("\n[written] v3_dynamics/results/")
