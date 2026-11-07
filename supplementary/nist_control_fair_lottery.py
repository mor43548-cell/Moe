# SUPPLEMENTARY (not part of the requested test protocol): checks whether the NIST results are
# reproduced by a perfectly fair simulated 6-of-37 lottery with identical encoding. Run: python supplementary/nist_control_fair_lottery.py 60
# Control: NIST battery on simulated FAIR 6-of-37 lotteries, same layout/encoding as the real data.
import sys, json, numpy as np
sys.path.insert(0, '/home/user/Moe/scripts')
import sp800_22 as T
rng = np.random.default_rng(2026); SIMS = int(sys.argv[1])
perm = [5, 4, 1, 2, 3, 0]   # file column order Y,T,W,E,R,Q = sorted positions 6,5,2,3,4,1
res = {}
for s in range(SIMS):
    d = np.sort(np.argsort(rng.random((328, 37)), axis=1)[:, :6] + 1, axis=1)[:, perm]
    bits = np.array([int(b) for v in d.ravel() for b in format(int(v), '06b')], dtype=np.int8)
    for f in T.ALL:
        r = f(bits)
        if not r['applicable']: st = 'NA'
        else: st = 'PASS' if all(p >= 0.01 for p in r['p_values']) else 'FAIL'
        res.setdefault(r['name'], []).append(st)
for k, v in res.items():
    print(f"{k:40s} FAIL {v.count('FAIL')}/{len(v)}  PASS {v.count('PASS')}  NA {v.count('NA')}")
