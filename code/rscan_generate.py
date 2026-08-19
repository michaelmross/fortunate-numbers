#!/usr/bin/env python3
"""R-scan: first-arrival T at anchors N = R * p_k# for several odd multipliers R,
over a fixed band of k. Discriminates 'bare primorial R=1 is uniquely anomalous'
from 'anomaly tracks a property of the multiplier R'.

Each record: {k, R, pk1, F (upper offset), f (lower offset)}.
The analyzer computes T = ncand * q_k(R) with the correct N/phi(N) for R.

Usage:
  python3 rscan_generate.py --kmin 1250 --kmax 1750 --step 10 \
          --Rset 1,5,7,11,35 --procs 8
Checkpoint: rscan.jsonl (append-only, resumable).

Cost note: each (k,R) does two first-arrival searches on a ~(k*log k)-digit
number, ~hundreds of BPSW tests each. At k~1500 budget ~1-3 s per (k,R) per core.
A 50-k x 5-R grid is ~500 cells -> minutes-to-an-hour at --procs 8.
"""
import gmpy2, json, os, sys, time
from multiprocessing import Pool

def arg(flag, default, cast=str):
    if f'--{flag}' in sys.argv:
        return cast(sys.argv[sys.argv.index(f'--{flag}')+1])
    return default

KMIN = arg('kmin', 1250, int)
KMAX = arg('kmax', 1750, int)
STEP = arg('step', 10, int)
RSET = [int(x) for x in arg('Rset', '1,5,7,11,35').split(',')]
PROCS = arg('procs', None, lambda x: int(x))
CKPT = 'rscan.jsonl'

def primorial_upto_k(kk):
    p = 0
    for _ in range(kk):
        p = int(gmpy2.next_prime(p))
    return int(gmpy2.primorial(p)), p

def one_cell(args):
    kk, R = args
    Nprim, pk = primorial_upto_k(kk)
    N = R * Nprim
    pk1 = int(gmpy2.next_prime(pk))
    m = pk1
    while not gmpy2.is_prime(N + m):
        m = int(gmpy2.next_prime(m))
    F = m
    m = pk1
    while not gmpy2.is_prime(N - m):
        m = int(gmpy2.next_prime(m))
    return dict(k=kk, R=R, pk1=pk1, F=F, f=m)

def main():
    done = set()
    if os.path.exists(CKPT):
        for line in open(CKPT):
            r = json.loads(line); done.add((r['k'], r['R']))
    grid = [(k, R) for k in range(KMIN, KMAX+1, STEP) for R in RSET
            if (k, R) not in done]
    print(f"{len(done)} done, {len(grid)} cells to compute "
          f"(k {KMIN}-{KMAX} step {STEP}, R={RSET})", flush=True)
    t0 = time.time()
    with Pool(PROCS) as pool, open(CKPT, 'a') as out:
        for i, rec in enumerate(pool.imap_unordered(one_cell, grid)):
            out.write(json.dumps(rec) + '\n'); out.flush()
            if (i+1) % 20 == 0:
                print(f"{i+1}/{len(grid)}, {time.time()-t0:.0f}s "
                      f"(latest k={rec['k']} R={rec['R']})", flush=True)
    print("complete.", flush=True)

if __name__ == '__main__':
    main()
