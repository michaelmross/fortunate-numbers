# fortunate_all.csv

Comprehensive dataset for "The Distribution of Two-Sided Fortunate Numbers
and the Goldbach Amplification Constant for Primorials and Their Multiples."

One row per (anchor multiplier R, level k). Each row certifies two probable
primes: the upper arrival R*p_k# + F and the lower arrival R*p_k# - f.
4,290 rows; 8,580 probable primes; anchors up to ~7,500 decimal digits.

## Columns

| column     | meaning |
|------------|---------|
| `R`        | anchor multiplier; anchor is N = R * p_k#. R=1 is the bare primorial. |
| `k`        | level: number of primes in the primorial. |
| `p_k`      | k-th prime. |
| `p_k1`     | (k+1)-th prime = smallest candidate offset. |
| `digits_N` | decimal digits of the anchor N (= (theta(p_k)+log R)/log 10). |
| `F`        | least m > 1 with N + m prime (upper Fortunate number). |
| `f`        | least m > 1 with N - m prime (lower Fortunate number). |
| `u_up`     | normalized offset (F/p_k1 - 1). |
| `u_dn`     | normalized offset (f/p_k1 - 1). |
| `ncand_up` | number of prime candidates up to and including F. |
| `ncand_dn` | number of prime candidates up to and including f. |
| `q_k`      | presieved Hardy-Littlewood rate (N/phi(N))/log N. |
| `T_up`     | normalized first arrival = ncand_up * q_k  (~ Exp(1) under the law). |
| `T_dn`     | normalized first arrival = ncand_dn * q_k. |

## Reconstructing a prime

    N_upper = R * primorial(p_k) + F      # e.g. Python: gmpy2.primorial(p_k)
    N_lower = R * primorial(p_k) - f

All values are Baillie-PSW probable primes. Offsets below p_k1^2 are
provably prime given the primality of the anchor's complement (see paper).

## Coverage by anchor

    R=1   : k = 3..2000        bare primorial, full density
    R=3   : k = 600..2000      dense to 1732, stride 2 on 1732..2000
    R=5   : k = 600..1732 dense, plus stride-10 points to 1750
    R=7   : k = 600..1200 & 1481..1729 (stride 4), stride-10 points to 1750
    R=11  : k = 600..1200 & 1481..1729 (stride 4), stride-10 points to 1750
    R=35  : k = 1250..1750 (stride 10)
    R=2,13,15,21,33,39,55,77 : k = 1481..1729 (stride 4)

Coverage combines the final measurement campaigns with an earlier
multiplier scan (R in {1,5,7,11,35}, stride 10 over 1250..1750); the two
agree exactly on all 100 overlapping cells. Because strides differ by
provenance, select on (R, k) when a uniform-stride series is needed.
4,445 rows; 8,890 probable primes.

Derived columns recomputed from (R,k,F,f) alone; primes regenerate
deterministically. Note the eight stride-4 campaign anchors terminate at
k=1729 (a consequence of stepping by 4 from k=1481); this affects no
statistic in the paper.
