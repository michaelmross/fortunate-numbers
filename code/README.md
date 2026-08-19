# Code for Data Generation

"Two-Sided Fortunate Numbers and Goldbach Amplification for Primorial Multiples"

Python 3, requires `gmpy2` (generation) and `numpy`, `scipy` (analysis).
All generators write append-only, resumable JSONL checkpoints; all
analyzers read them and recompute every derived quantity from
`(k, R, pk1, F, f)` alone. Reproduces every figure and table in the paper.

## scripts/  (canonical)

### Generation
- **fortunate_generate.py** — bare primorial anchors N = p_k#. Two-sided
  first prime arrivals; self-checks the first 58 upper values against
  OEIS A005235. Writes `fortunate.jsonl`.
      python fortunate_generate.py 2000 --procs 8
- **control_generate.py** — multiplier anchors N = R·p_k#. Writes
  `fortunate_R{R}.jsonl`.
      python control_generate.py 1732 --kmin 1200 --step 1 --R 3
- **rscan_generate.py** — multi-R scan at one stride over a k-band; writes
  `rscan.jsonl`. Used for the pre-registered δ_u measurement.

  Generators share resource-politeness flags: `--procs N` (default
  cpu−2), `--priority {idle,belownormal,normal}` (default belownormal),
  `--breathe S,E` (synchronized pause S s every E s, for Task Scheduler
  idle detection).

### Analysis
- **fortunate_analyze_v31.py** — primary analyzer. Candidate-count
  statistic T = (candidates to first success)·q_k, tested against
  Exp(1): pooled/upper/lower means, KS, tail binomials, k-band table,
  trend test, out-of-sample slices. (v3.1 fixes a k-label alignment bug
  in v3; see superseded/.)
- **control_analyze.py** — same statistic with q_k(R) = [∏p/(p−1)] /
  (θ(p_k)+log R) for a control anchor.
- **compare_anchors.py** — matched-k two-sample contrast of primorial vs
  a control (the δ_p statistic).
- **rscan_analyze.py** — per-R means from a scan file (the δ_u table).
- **stride_analyze.py** — independence battery: strided means,
  autocorrelation (raw and detrended), batch-means variance inflation,
  moving-block bootstrap. Settles whether naive significance is honest.

### Utility
- **expand_primes.py** — reconstruct the probable primes a record
  certifies, single or whole-file, with `--full` for complete decimals.
      python expand_primes.py 1467 12281 --full

## Reproduction
    python fortunate_generate.py 2000
    python fortunate_analyze_v31.py fortunate.jsonl      # law + deep slices
    python control_generate.py 1732 --kmin 1200 --step 1 --R 3
    python compare_anchors.py fortunate.jsonl fortunate_R3.jsonl
    python stride_analyze.py fortunate.jsonl             # independence
