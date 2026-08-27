# Two-Sided Fortunate Numbers

Data, code, and results for **[Two-Sided Fortunate Numbers and
Goldbach Amplification for Primorial Multiples](https://doi.org/10.5281/zenodo.22018781)**
(M. M. Ross, 2026).

[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.1339494240-blue.svg)](https://doi.org/10.5281/zenodo.22018881)

For an anchor `N = R · p_k#` (a multiple of the k-th primorial), every
Goldbach complement `N − p` is coprime to `p_k#`, and the least offsets
`m > 1` with `N ± m` prime — the two-sided, multiplier-generalized
**Fortunate numbers** — are themselves prime whenever they fall below
`p_{k+1}²`. This [repository](https://github.com/michaelmross/fortunate-numbers) measures the distribution of those offsets.

## Main result

Normalize each first arrival by the presieved Hardy–Littlewood rate
`q_k(R) = (N/φ(N)) / log N`. The paper conjectures, derives (under a uniform
Hardy–Littlewood hypothesis), and measures the law

    T = (prime candidates up to the first success) · q_k(R)  ⟶  Exp(1),

whose rate constant is the Goldbach amplification factor `e^γ`. Across a
fourteen-anchor dataset of **8,946 probable primes** (up to
**7,482 digits** at k ≤ 2000), the constant is confirmed to three digits
below k ≈ 1200 and again beyond k ≈ 1732 — but on the intervening band of
levels the whole anchor family departs from the law by a common
**6.4 % ± 1.8 %**, established by three pre-registered measurements,
bounded on both sides by conformity, and unexplained by any correction the
derivation supplies.

## Layout

    data/      raw append-only JSONL checkpoints, one per anchor multiplier R
    results/   fortunate_all.csv — the merged, analysis-ready dataset
    code/      generators, analyzers, and the prime expander

## Data

Each JSONL record is a ~30-byte certificate of two probable primes:

    {"k": 1467, "R": 1, "pk1": 12251, "F": 12281, "f": 13171}
      →  upper prime  R·p_k# + F,   lower prime  R·p_k# − f

`results/fortunate_all.csv` merges every anchor into one table (4,473 rows)
with all derived columns precomputed — `q_k`, candidate counts, the
normalized statistics `T_up`/`T_dn`, offsets, and anchor digit-lengths — so
the paper's tables reproduce from a spreadsheet. Column definitions and
per-anchor coverage are in `results/README.md`.

Coverage: the bare primorial (`R = 1`) at every level 3 ≤ k ≤ 2000; `R = 3`
to k = 2000; `R = 5, 7, 11` and ten further multipliers over the deep band;
plus an earlier multiplier scan (`rscan.jsonl`, R ∈ {1,5,7,11,35}). The two
runs agree exactly on all overlapping cells.

## Code

Requires Python 3 with `gmpy2` (generation) and `numpy`, `scipy` (analysis).
Generators write resumable checkpoints; analyzers recompute every quantity
from `(k, R, pk1, F, f)` alone. Full script-by-script guide in
`code/README.md`. Quick reproduction:

    python code/fortunate_generate.py 2000
    python code/fortunate_analyze_v31.py data/fortunate_R1.jsonl
    python code/control_generate.py 1732 --kmin 1200 --step 1 --R 3
    python code/compare_anchors.py data/fortunate_R1.jsonl data/fortunate_R3.jsonl
    python code/stride_analyze.py data/fortunate_R1.jsonl

To reconstruct any prime in full:

    python code/expand_primes.py 1467 12281 --full

## License

Data and code released under [choose: CC-BY-4.0 / MIT]; see LICENSE.
