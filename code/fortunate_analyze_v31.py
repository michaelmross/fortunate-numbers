#!/usr/bin/env python3
"""v3: candidate-count (geometric clock) test, the primary form of the law.
T_k = (#prime candidates up to winner) * q_k,  q_k = [prod_{p<=p_k} p/(p-1)]/theta(p_k)
Conjecture: T ~ Exp(1), zero free parameters.
Usage: python3 fortunate_analyze_v3.py [fortunate.jsonl] [oos_split_k]"""
import json, math, sys, bisect
import numpy as np
from scipy import stats

PATH = sys.argv[1] if len(sys.argv) > 1 else "fortunate.jsonl"
SPLIT = int(sys.argv[2]) if len(sys.argv) > 2 else 0

recs = {}
for line in open(PATH):
    r = json.loads(line); recs[r["k"]] = r
rows = [recs[k] for k in sorted(recs)]
kmax = rows[-1]["k"]

LIM = max(200000, 20 * max(max(r["F"], r["f"]) for r in rows))
LIM = max(LIM, 8 * rows[-1]["pk1"])
sieve = bytearray([1]) * LIM
sieve[0] = sieve[1] = 0
for i in range(2, int(LIM ** .5) + 1):
    if sieve[i]:
        sieve[i*i::i] = bytearray(len(sieve[i*i::i]))
primes = [i for i in range(LIM) if sieve[i]]
pidx = {p: i for i, p in enumerate(primes)}
theta = np.cumsum([math.log(p) for p in primes])
lA = np.cumsum([math.log(p / (p - 1)) for p in primes])

assert all(r["pk1"] == primes[r["k"]] for r in rows), "pk1/sieve mismatch"

def T_of(r, key):
    k = r["k"]
    q = math.exp(lA[k-1]) / theta[k-1]
    return (pidx[r[key]] - pidx[r["pk1"]] + 1) * q

Tup = np.array([T_of(r, "F") for r in rows])
Tdn = np.array([T_of(r, "f") for r in rows])
T = np.concatenate([Tup, Tdn]); n = len(T)
print(f"n_k={len(rows)}, k<= {kmax}, pooled n={n}.  Target: T ~ Exp(1)\n")

def rep(name, x):
    D, p = stats.kstest(x, lambda v: 1 - np.exp(-np.clip(v, 0, None)))
    se = x.std(ddof=1) / len(x) ** .5
    print(f"{name:24s} n={len(x):5d} mean={x.mean():.4f} "
          f"({abs(x.mean()-1)/se:.2f} se)  KS D={D:.4f} p={p:.3f}")

rep("pooled", T); rep("upper (Fortunate)", Tup); rep("lower (Goldbach)", Tdn)
if SPLIT:
    oos = np.array([T_of(r, s) for r in rows if r["k"] > SPLIT
                    for s in ("F", "f")])
    rep(f"out-of-sample k>{SPLIT}", oos)

print("\ntail, exact one-sided binomial:")
for t0 in (1.5, 2, 3, 4, 5):
    obs = int((T > t0).sum())
    pv = stats.binomtest(obs, n, math.exp(-t0), alternative="greater").pvalue
    print(f"  T>{t0}: obs {obs:4d}  exp {n*math.exp(-t0):7.1f}  p={pv:.3f}")

nb = max(4, len(rows) // 150)
edges = np.linspace(rows[0]["k"], kmax, nb + 1).astype(int)
ks_rows = np.array([r["k"] for r in rows]); ks_arr = np.concatenate([ks_rows, ks_rows])  # v3.1 FIX: align with T=concat(up,dn)
print("\nmean T by k-band (thy 1.0, se ~ 1/sqrt(n_band)):")
for lo, hi in zip(edges[:-1], edges[1:]):
    sel = (ks_arr >= lo) & (ks_arr <= hi)
    print(f"  k in [{lo:4d},{hi:4d}]: n={sel.sum():4d}  mean={T[sel].mean():.4f}")

i = int(T.argmax()); k_at = ks_arr[i]
print(f"\nmax T={T.max():.3f} at k={k_at}; Gumbel location ln(n)={math.log(n):.3f}; "
      f"P(max >= obs) ~ {1-math.exp(-n*math.exp(-T.max())):.3f}")

# v3.1 additions: formal trend and late-slice diagnostics
from scipy.stats import linregress
sl, ic, rv, pv, se_ = linregress(ks_arr, T)
print(f"\ntrend test: slope={sl:.3e}/k, p={pv:.3f}")
for cut in (1000, 1216, 1333):
    sel = ks_arr > cut
    if sel.sum() > 50:
        Dc, pc = stats.kstest(T[sel], lambda v: 1 - np.exp(-np.clip(v, 0, None)))
        print(f"slice k>{cut}: n={sel.sum()}  mean={T[sel].mean():.4f}  KS p={pc:.3f}")
