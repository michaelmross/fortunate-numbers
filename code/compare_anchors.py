#!/usr/bin/env python3
"""Matched-k two-sample comparison of primorial vs control anchors.
Usage: python3 compare_anchors.py fortunate.jsonl fortunate_R3.jsonl"""
import json, math, sys
import numpy as np
from scipy import stats

def load(p):
    d={}
    for line in open(p):
        r=json.loads(line); d[r['k']]=r
    return d
prim, ctrl = load(sys.argv[1]), load(sys.argv[2])
R = next(iter(ctrl.values())).get('R',3)

LIM=300000
s=bytearray([1])*LIM; s[0]=s[1]=0
for i in range(2,int(LIM**.5)+1):
    if s[i]: s[i*i::i]=bytearray(len(s[i*i::i]))
primes=[i for i in range(LIM) if s[i]]
pidx={p:i for i,p in enumerate(primes)}
theta=np.cumsum([math.log(p) for p in primes])
lA=np.cumsum([math.log(p/(p-1)) for p in primes])

def T_of(r,key,logR=0.0):
    q=math.exp(lA[r['k']-1])/(theta[r['k']-1]+logR)
    return (pidx[r[key]]-pidx[r['pk1']]+1)*q

kk=sorted(set(prim)&set(ctrl))
Tp=np.array([T_of(prim[k],key) for k in kk for key in ('F','f')])
Tc=np.array([T_of(ctrl[k],key,math.log(R)) for k in kk for key in ('F','f')])
print(f"matched k: {len(kk)} values in [{kk[0]},{kk[-1]}]; n={len(Tp)} per anchor")
print(f"primorial mean={Tp.mean():.4f}   control(R={R}) mean={Tc.mean():.4f}")
t,pv=stats.ttest_ind(Tp,Tc); D,pk=stats.ks_2samp(Tp,Tc)
print(f"two-sample: t={t:.2f} p={pv:.4f}   KS D={D:.4f} p={pk:.4f}")
for name,x in (("primorial",Tp),(f"control R={R}",Tc)):
    D1,p1=stats.kstest(x,lambda v:1-np.exp(-np.clip(v,0,None)))
    print(f"  {name} vs Exp(1): mean-z={(1-x.mean())*math.sqrt(len(x)):+.2f}  KS p={p1:.4f}")
