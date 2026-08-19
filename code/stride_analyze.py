#!/usr/bin/env python3
"""Stride & independence analyzer for the primorial first-arrival data.

Settles the one open question: are T_k independent across k, or does
autocorrelation inflate the significance of the k>1216 deficit?

Outputs:
  1. strided means (point-estimate stability under decorrelating subsamples)
  2. autocorrelation function of T_k, raw and locally-detrended
     (raw ACF conflates smooth mean-structure with noise correlation;
      residual ACF after a moving-average detrend isolates noise correlation,
      which is what determines whether naive z-scores are honest)
  3. batch-means variance inflation factor (VIF) and corrected z
  4. moving-block bootstrap z for the k>1216 mean deficit

Decision rule printed at the end.

Usage: python3 stride_analyze.py fortunate.jsonl
"""
import json, math, sys
import numpy as np

PATH = sys.argv[1] if len(sys.argv) > 1 else 'fortunate.jsonl'
CUT = int(sys.argv[2]) if len(sys.argv) > 2 else 1216   # anomaly slice boundary

recs = {}
for line in open(PATH):
    r = json.loads(line); recs[r['k']] = r
rows = [recs[k] for k in sorted(recs)]

LIM = 2_000_000
s = bytearray([1])*LIM; s[0] = s[1] = 0
for i in range(2, int(LIM**.5)+1):
    if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
plist = [i for i in range(LIM) if s[i]]
pidx = {p: i for i, p in enumerate(plist)}
theta = np.cumsum([math.log(p) for p in plist])
lA = np.cumsum([math.log(p/(p-1)) for p in plist])

def T_of(r, key):
    q = math.exp(lA[r['k']-1]) / theta[r['k']-1]
    return (pidx[r[key]] - pidx[r['pk1']] + 1) * q

ks = np.array([r['k'] for r in rows])
Tu = np.array([T_of(r, 'F') for r in rows])
Td = np.array([T_of(r, 'f') for r in rows])
Tk = (Tu + Td) / 2.0                     # per-k statistic
print(f"n_k={len(rows)}, k in [{ks[0]},{ks[-1]}]")

# ---------- 1. strided means on the anomaly slice ----------
sel = ks > CUT
print(f"\n1. strided means, k>{CUT} (pooled both sides; naive se in parens):")
for stride in (1, 2, 5, 10, 20):
    sub = np.where(sel)[0][::stride]
    x = np.concatenate([Tu[sub], Td[sub]])
    se = x.std(ddof=1)/math.sqrt(len(x))
    print(f"   stride {stride:2d}: mean={x.mean():.4f} ({se:.4f})  n={len(x)}")
print("   [mean stable across strides + se growing like sqrt(stride) = expected")
print("    for BOTH a real effect and correlated noise; see ACF below to split]")

# ---------- 2. autocorrelation ----------
def acf(x, maxlag):
    x = x - x.mean(); v = np.dot(x, x)/len(x)
    return [float(np.dot(x[:-l], x[l:])/((len(x)-l)*v)) for l in range(1, maxlag+1)]

def detrend(x, w=101):
    # centered moving average; residuals isolate noise from smooth mean drift
    pad = w//2
    xp = np.concatenate([x[:pad][::-1], x, x[-pad:][::-1]])
    ma = np.convolve(xp, np.ones(w)/w, mode='valid')
    return x - ma

MAXLAG = 20
print(f"\n2. ACF of per-k statistic T_k (se under independence ~ "
      f"{1/math.sqrt(len(Tk)):.3f}):")
raw = acf(Tk, MAXLAG)
res = acf(detrend(Tk), MAXLAG)
print("   lag :", "  ".join(f"{l:5d}" for l in (1,2,3,5,10,20)))
print("   raw :", "  ".join(f"{raw[l-1]:+.3f}" for l in (1,2,3,5,10,20)))
print("   detr:", "  ".join(f"{res[l-1]:+.3f}" for l in (1,2,3,5,10,20)))
print("   sides at same k (should be ~0 if independent):",
      f"corr(Tu,Td) = {np.corrcoef(Tu,Td)[0,1]:+.3f}")

# ---------- 3. batch means VIF on the anomaly slice ----------
x = Tk[sel]
print(f"\n3. batch-means variance inflation, k>{CUT} (n_k={len(x)}):")
naive_var = x.var(ddof=1)
for L in (8, 16, 32):
    nb = len(x)//L
    bm = x[:nb*L].reshape(nb, L).mean(axis=1)
    vif = L*bm.var(ddof=1)/naive_var
    print(f"   L={L:3d}: VIF={vif:.2f}  (1.0 = independent; "
          f"corrected z divides by sqrt(VIF))")

# ---------- 4. block bootstrap for the pooled deficit ----------
pool = np.concatenate([Tu[sel], Td[sel]])
naive_z = (1 - pool.mean())*math.sqrt(len(pool))/pool.std(ddof=1)
rng = np.random.default_rng(0)
B, blk = 4000, 20
n = len(x)
means = []
for _ in range(B):
    idx = []
    while len(idx) < n:
        st = rng.integers(0, n - blk)
        idx.extend(range(st, st + blk))
    xs = x[np.array(idx[:n])]
    means.append(xs.mean())
means = np.array(means)
boot_z = (1 - x.mean())/means.std(ddof=1)
print(f"\n4. k>{CUT} deficit: naive z={naive_z:.2f} (pooled sides), "
      f"block-bootstrap z={boot_z:.2f} (per-k, block={blk})")

# ---------- decision ----------
print("\nDECISION RULE:")
print("  detrended ACF ~ 0 at all lags AND VIF ~ 1 AND bootstrap z ~ naive z")
print("    -> T_k are independent; the k>1216 deficit's original significance")
print("       (4.4 sigma) STANDS; the 'autocorrelation dissolution' is refuted;")
print("       the anomaly is back on the table and the next step is a dense")
print("       (step-1 or step-2) control run at one multiplier, n>=690.")
print("  detrended ACF substantially positive / VIF >> 1 / bootstrap z << naive")
print("    -> correlation is real; corrected significance is the bootstrap z;")
print("       if that is < ~2, the deficit is not established and the law's")
print("       clean verification stands as the paper's result.")
