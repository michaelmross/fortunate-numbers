#!/usr/bin/env python3
"""Analyze control-anchor data: T = ncand * q, q = [prod p/(p-1)]/(theta(p_k)+log R).
Usage: python3 control_analyze.py fortunate_R3.jsonl"""
import json, math, sys
import numpy as np
from scipy import stats

PATH = sys.argv[1] if len(sys.argv) > 1 else 'fortunate_R3.jsonl'
recs = {}
for line in open(PATH):
    r = json.loads(line); recs[r['k']] = r
rows = [recs[k] for k in sorted(recs)]
R = rows[0]['R']

LIM = 300000
s = bytearray([1])*LIM; s[0]=s[1]=0
for i in range(2, int(LIM**.5)+1):
    if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
primes = [i for i in range(LIM) if s[i]]
pidx = {p:i for i,p in enumerate(primes)}
theta = np.cumsum([math.log(p) for p in primes])
lA = np.cumsum([math.log(p/(p-1)) for p in primes])
assert all(r['pk1'] == primes[r['k']] for r in rows)

def T_of(r, key):
    q = math.exp(lA[r['k']-1]) / (theta[r['k']-1] + math.log(R))
    return (pidx[r[key]] - pidx[r['pk1']] + 1) * q

ks = np.array([r['k'] for r in rows for _ in (0,1)])
T = np.array([T_of(r,key) for r in rows for key in ('F','f')])
n = len(T)
D,p = stats.kstest(T, lambda v: 1-np.exp(-np.clip(v,0,None)))
print(f"anchor R={R}: n_k={len(rows)}, k in [{rows[0]['k']},{rows[-1]['k']}], pooled n={n}")
print(f"mean T={T.mean():.4f} ({(1-T.mean())*math.sqrt(n):+.2f} z-dev from 1)  KS p={p:.4f}")
Ds,ps = stats.kstest(T/T.mean(), lambda v: 1-np.exp(-np.clip(v,0,None)))
print(f"shape after rescale: KS p={ps:.3f}")
edges = np.linspace(rows[0]['k'], rows[-1]['k'], 6).astype(int)
for lo,hi in zip(edges[:-1], edges[1:]):
    sel = (ks>=lo)&(ks<=hi)
    if sel.sum()>10:
        print(f"  k in [{lo:4d},{hi:4d}]: n={sel.sum():4d}  mean={T[sel].mean():.4f}")
print("\nCOMPARISON KEY: primorial anchor showed mean ~0.85 for k>1216.")
print("Control ~1.00 here -> effect is anchor-specific (p_k# itself).")
print("Control ~0.85 here -> systematic in the shared method; hunt inward.")
