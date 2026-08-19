#!/usr/bin/env python3
"""Analyze rscan.jsonl: mean T per multiplier R, with correct q_k(R).
Usage: python3 rscan_analyze.py rscan.jsonl"""
import json, math, sys
import numpy as np
from scipy import stats

PATH = sys.argv[1] if len(sys.argv) > 1 else 'rscan.jsonl'
rows = [json.loads(l) for l in open(PATH)]

LIM = 2_000_000
s = bytearray([1])*LIM; s[0]=s[1]=0
for i in range(2, int(LIM**.5)+1):
    if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
plist = [i for i in range(LIM) if s[i]]
pidx = {p:i for i,p in enumerate(plist)}
theta = np.cumsum([math.log(p) for p in plist])
lA = np.cumsum([math.log(p/(p-1)) for p in plist])

def qval(kk, R):
    logNphi = lA[kk-1]; Rf = R
    for pp in plist:
        if Rf == 1: break
        if Rf % pp == 0:
            if pp > plist[kk-1]:
                logNphi += math.log(pp/(pp-1))
            while Rf % pp == 0: Rf //= pp
        if pp*pp > Rf and Rf > 1:
            if Rf > plist[kk-1]: logNphi += math.log(Rf/(Rf-1))
            Rf = 1
    return math.exp(logNphi) / (math.log(R) + theta[kk-1])

by_R = {}
for r in rows:
    q = qval(r['k'], r['R'])
    i0 = pidx[r['pk1']]
    Tu = (pidx[r['F']] - i0 + 1) * q
    Td = (pidx[r['f']] - i0 + 1) * q
    by_R.setdefault(r['R'], []).extend([Tu, Td])

print(f"{'R':>4s} {'mean_T':>8s} {'z_dev':>7s} {'KS_p':>7s} {'n':>5s}")
for R in sorted(by_R):
    x = np.array(by_R[R]); n = len(x)
    D, p = stats.kstest(x, lambda v: 1-np.exp(-np.clip(v, 0, None)))
    z = (1 - x.mean()) * math.sqrt(n)
    tag = "  <-- bare primorial" if R == 1 else ""
    print(f"{R:4d} {x.mean():8.4f} {z:+7.2f} {p:7.4f} {n:5d}{tag}")

if 1 in by_R:
    others = [t for R in by_R if R != 1 for t in by_R[R]]
    if others:
        t, p = stats.ttest_ind(by_R[1], others)
        print(f"\nR=1 vs all R>1: mean {np.mean(by_R[1]):.4f} vs "
              f"{np.mean(others):.4f}, t={t:.2f} p={p:.4f}")
        print("INTERPRETATION:")
        print("  R=1 alone anomalous (low), others ~1.0 -> effect is the BARE")
        print("    primorial; report as finding.")
        print("  anomaly present across R -> multiplier/method effect; the R=1")
        print("    'anomaly' is not primorial-specific.")
