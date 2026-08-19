#!/usr/bin/env python3
"""Two-sided Fortunate numbers (primorial anchor), resource-polite v2.

Same politeness mechanisms as control_generate.py v2: workers at reduced OS
priority (default belownormal), two cores left free by default, optional
synchronized --breathe windows for Task Scheduler idle detection.

Usage:
  python fortunate_generate.py 2000
  python fortunate_generate.py 2000 --procs 8 --priority idle
  python fortunate_generate.py 2000 --breathe 90,900

Checkpoint: fortunate.jsonl (unchanged format; resumes existing data).
First 58 upper-side terms self-checked against OEIS A005235.
"""
import gmpy2, json, os, sys, time
from multiprocessing import Pool

A005235_58 = [3,5,7,13,23,17,19,23,37,61,67,61,71,47,107,59,61,109,89,103,
              79,151,197,101,103,233,223,127,223,191,163,229,643,239,157,
              167,439,239,199,191,199,383,233,751,313,773,607,313,383,293,
              443,331,283,277,271,401,307,331]

def parse():
    a = {'K': int(sys.argv[1]), 'procs': None,
         'priority': 'belownormal', 'breathe': None}
    if '--procs' in sys.argv:
        a['procs'] = int(sys.argv[sys.argv.index('--procs') + 1])
    if '--priority' in sys.argv:
        a['priority'] = sys.argv[sys.argv.index('--priority') + 1]
    if '--breathe' in sys.argv:
        S, E = sys.argv[sys.argv.index('--breathe') + 1].split(',')
        a['breathe'] = (float(S), float(E))
    if a['procs'] is None:
        a['procs'] = max(1, (os.cpu_count() or 4) - 2)
    return a

A = parse()
CKPT = 'fortunate.jsonl'

def set_priority(level):
    try:
        if sys.platform == 'win32':
            import ctypes
            classes = {'idle': 0x40, 'belownormal': 0x4000, 'normal': 0x20}
            ctypes.windll.kernel32.SetPriorityClass(
                ctypes.windll.kernel32.GetCurrentProcess(), classes[level])
        else:
            nice = {'idle': 19, 'belownormal': 10, 'normal': 0}[level]
            if nice:
                os.nice(nice)
    except Exception as e:
        print(f"[warn] could not set priority: {e}", flush=True)

def init_worker():
    set_priority(A['priority'])

def maybe_breathe():
    if A['breathe']:
        S, E = A['breathe']
        t = time.time() % E
        if t < S:
            time.sleep(S - t)

def one_k(k):
    p = 0
    for _ in range(k):
        p = int(gmpy2.next_prime(p))
    N = int(gmpy2.primorial(p))
    pk1 = int(gmpy2.next_prime(p))
    m = pk1
    while True:
        maybe_breathe()
        if gmpy2.is_prime(N + m):
            break
        m = int(gmpy2.next_prime(m))
    F = m
    m = pk1
    while True:
        maybe_breathe()
        if gmpy2.is_prime(N - m):
            break
        m = int(gmpy2.next_prime(m))
    return dict(k=k, pk1=pk1, F=F, f=m)

def main():
    set_priority(A['priority'])
    done = set()
    if os.path.exists(CKPT):
        for line in open(CKPT):
            done.add(json.loads(line)['k'])
    todo = [k for k in range(3, A['K'] + 1) if k not in done]
    print(f"{len(done)} done, {len(todo)} to compute "
          f"(procs={A['procs']}, priority={A['priority']}"
          f"{', breathe='+str(A['breathe']) if A['breathe'] else ''})",
          flush=True)
    t0 = time.time()
    with Pool(A['procs'], initializer=init_worker) as pool, \
         open(CKPT, 'a') as out:
        for i, rec in enumerate(pool.imap_unordered(one_k, todo)):
            if rec['k'] <= 58:
                assert rec['F'] == A005235_58[rec['k'] - 1], \
                    f"MISMATCH vs A005235 at k={rec['k']}"
            out.write(json.dumps(rec) + '\n'); out.flush()
            if (i + 1) % 25 == 0:
                print(f"{i+1}/{len(todo)}, {time.time()-t0:.0f}s "
                      f"(latest k={rec['k']})", flush=True)
    print("complete.", flush=True)

if __name__ == '__main__':
    main()
