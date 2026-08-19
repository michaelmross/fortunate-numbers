#!/usr/bin/env python3
"""Control-anchor experiment: two-sided first prime arrivals at N = R * p_k#.

v2 (resource-polite): workers run at reduced OS scheduling priority
(default: belownormal) and two CPU cores are left free by default, so the
run yields instantly to interactive use and normal-priority background jobs
while still consuming all genuinely idle cycles. Optional --breathe gives
Windows Task Scheduler's "only when idle" condition periodic windows in
which system CPU visibly drops to ~0.

Usage:
  python control_generate.py 1732 --kmin 1481 --step 4 --R 7
  python control_generate.py 1732 --kmin 1481 --step 4 --R 7 --procs 8 --priority idle
  python control_generate.py 1732 --kmin 1481 --step 4 --R 7 --breathe 90,900
      (all workers pause together for 90 s every 900 s of wall clock)

Defaults: procs = cpu_count - 2, priority = belownormal, no breathe.
Checkpoint: fortunate_R{R}.jsonl (append-only, resumable; format unchanged,
existing checkpoints resume as before).
"""
import gmpy2, json, os, sys, time
from multiprocessing import Pool

def parse():
    a = {'K': int(sys.argv[1]), 'R': 3, 'procs': None, 'kmin': 3, 'step': 1,
         'priority': 'belownormal', 'breathe': None}
    for flag in ('R', 'procs', 'kmin', 'step'):
        if f'--{flag}' in sys.argv:
            a[flag] = int(sys.argv[sys.argv.index(f'--{flag}') + 1])
    if '--priority' in sys.argv:
        a['priority'] = sys.argv[sys.argv.index('--priority') + 1]
    if '--breathe' in sys.argv:
        S, E = sys.argv[sys.argv.index('--breathe') + 1].split(',')
        a['breathe'] = (float(S), float(E))
    if a['procs'] is None:
        a['procs'] = max(1, (os.cpu_count() or 4) - 2)
    return a

A = parse()
CKPT = f"fortunate_R{A['R']}.jsonl"

def set_priority(level):
    """OS-level 'DoEvents': let anything at normal priority preempt us."""
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
    """Synchronized pool-wide pause: all workers share the wall clock, so
    they sleep in the same window and system CPU visibly drops."""
    if A['breathe']:
        S, E = A['breathe']
        t = time.time() % E
        if t < S:
            time.sleep(S - t)

def one_cell(k):
    p = 0
    for _ in range(k):
        p = int(gmpy2.next_prime(p))
    N = A['R'] * int(gmpy2.primorial(p))
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
    return dict(k=k, R=A['R'], pk1=pk1, F=F, f=m)

def main():
    set_priority(A['priority'])
    done = set()
    if os.path.exists(CKPT):
        for line in open(CKPT):
            done.add(json.loads(line)['k'])
    todo = [k for k in range(max(3, A['kmin']), A['K'] + 1, A['step'])
            if k not in done]
    print(f"{len(done)} done, {len(todo)} to compute "
          f"(R={A['R']}, procs={A['procs']}, priority={A['priority']}"
          f"{', breathe='+str(A['breathe']) if A['breathe'] else ''})",
          flush=True)
    t0 = time.time()
    with Pool(A['procs'], initializer=init_worker) as pool, \
         open(CKPT, 'a') as out:
        for i, rec in enumerate(pool.imap_unordered(one_cell, todo)):
            out.write(json.dumps(rec) + '\n'); out.flush()
            if (i + 1) % 25 == 0:
                print(f"{i+1}/{len(todo)}, {time.time()-t0:.0f}s "
                      f"(latest k={rec['k']})", flush=True)
    print("complete.", flush=True)

if __name__ == '__main__':
    main()
