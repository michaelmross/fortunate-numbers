#!/usr/bin/env python3
"""expand_primes.py -- reconstruct the probable primes stored implicitly in
the checkpoint records {k, R, pk1, F, f}. Each record certifies two primes:
R*p_k# + F (upper) and R*p_k# - f (lower).

Single record:
    python expand_primes.py 1467 12281              # k, m  (R=1, upper)
    python expand_primes.py 1467 12281 --R 3 --side lower
File mode (incremental primorial; fast for whole checkpoints):
    python expand_primes.py --file fortunate_R3.jsonl            # summaries
    python expand_primes.py --file fortunate_R3.jsonl --full     # full digits
        -> writes fortunate_R3_expanded.txt
"""
import sys, json, gmpy2

sys.set_int_max_str_digits(200000)

def summary(n):
    d = str(n)
    return f"{len(d)} digits: {d[:30]}...{d[-30:]}" if len(d) > 70 else d

def single(argv):
    k, m = int(argv[0]), int(argv[1])
    R = int(argv[argv.index('--R')+1]) if '--R' in argv else 1
    side = -1 if '--side' in argv and argv[argv.index('--side')+1] == 'lower' else 1
    p = 0
    for _ in range(k):
        p = int(gmpy2.next_prime(p))
    n = R * int(gmpy2.primorial(p)) + side * m
    tag = '+' if side > 0 else '-'
    body = str(n) if '--full' in argv else summary(n)
    print(f"{R}*p_{k}# {tag} {m}  =  {body}")
    print(f"BPSW probable prime: {bool(gmpy2.is_prime(n))}")

def filemode(path, full):
    recs = sorted((json.loads(l) for l in open(path)), key=lambda r: r['k'])
    out = path.rsplit('.', 1)[0] + '_expanded.txt'
    P, p, kcur = 1, 0, 0
    with open(out, 'w') as f:
        for r in recs:
            while kcur < r['k']:                    # incremental primorial
                p = int(gmpy2.next_prime(p)); P *= p; kcur += 1
            R = r.get('R', 1)
            for key, sign, tag in (('F', 1, '+'), ('f', -1, '-')):
                n = R * P + sign * r[key]
                body = str(n) if full else summary(n)
                f.write(f"k={r['k']} R={R} {tag}{r[key]}  {body}\n")
    print(f"wrote {2*len(recs)} primes -> {out}")

if __name__ == '__main__':
    if '--file' in sys.argv:
        filemode(sys.argv[sys.argv.index('--file')+1], '--full' in sys.argv)
    else:
        single(sys.argv[1:])
