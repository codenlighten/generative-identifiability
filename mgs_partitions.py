#!/usr/bin/env python3
"""
EXPERIMENT 24 -- does the crossing effect survive every observation partition?

Frozen before running.

Experiment 22 found, on ONE balanced 2+2 lumping, that two resets inside a
single observational block identify nothing while two that cross blocks
identify 58.55% of the class, at equal cost. That is one lumping of one
family, and the obvious question is whether the effect is a property of the
aggregation or an accident of that particular one.

Here we enumerate ALL 15 set partitions of the four hidden states and, for
each, all 16 subsets of resets. To make that exhaustive sweep affordable the
generator family is coarser than Experiment 20's: transition rows are drawn
from {0, 1/2, 1} rather than quarters, giving 10 rows and 10^4 = 10,000
generators, still enumerated completely. Horizon T = 5.

DEFINITIONS
  Pi          a set partition of {0,1,2,3}; the observation reports which
              block the hidden state lies in
  coverage    c(A) = number of blocks of Pi that A meets
  crossings   the number of pairs in A lying in different blocks

CRITERIA (fixed before running)
  X1  Within every (Pi, |A|) cell containing at least two coverage levels,
      mean H(G|O_A) is non-increasing in coverage. Reported as the fraction
      of such cells satisfying it; we pre-commit to calling the conjecture
      SUPPORTED only at >= 90%.
  X2  The effect is not an artefact of block sizes: it survives when
      restricted to partitions whose blocks are all the same size.
  X3  Coverage beats cardinality as a predictor: within fixed |A|, the
      spread of H across coverage levels exceeds 0.10 bits somewhere.
"""
import numpy as np, math, itertools
from collections import defaultdict

NS, T = 4, 5
ROWS = np.array([r for r in itertools.product(range(3), repeat=NS) if sum(r) == 2],
                dtype=np.int64)          # denominator 2
NGEN = len(ROWS) ** NS

def partitions(items):
    if len(items) == 1:
        yield [items]; return
    first, rest = items[0], items[1:]
    for p in partitions(rest):
        for i in range(len(p)):
            yield p[:i] + [[first] + p[i]] + p[i+1:]
        yield [[first]] + p

PARTS = [sorted([sorted(b) for b in p]) for p in partitions([0, 1, 2, 3])]
PARTS.sort(key=lambda p: (len(p), [len(b) for b in p]))

def block_of(part):
    m = np.empty(NS, dtype=np.int64)
    for bi, b in enumerate(part):
        for x in b:
            m[x] = bi
    return m

def laws(pi_vec, obsv, k):
    M = np.empty((NGEN, NS, NS), dtype=np.int64)
    idx = np.arange(NGEN)
    for pos in range(NS):
        M[:, pos, :] = ROWS[idx // (len(ROWS) ** (NS - 1 - pos)) % len(ROWS)]
    out = np.empty((NGEN, k ** T), dtype=np.int32)
    col = [0]
    def rec(alpha, depth):
        if depth == T:
            out[:, col[0]] = alpha.sum(axis=1); col[0] += 1; return
        nxt = np.matmul(alpha[:, None, :], M)[:, 0, :]
        for y in range(k):
            rec(np.where(obsv[None, :] == y, nxt, 0), depth + 1)
    for y0 in range(k):
        a0 = np.broadcast_to(np.where(obsv == y0, pi_vec, 0).astype(np.int64),
                             (NGEN, NS)).copy()
        rec(a0, 1)
    return out

def ids_of(a):
    v = np.ascontiguousarray(a).view(
        np.dtype((np.void, a.dtype.itemsize * a.shape[1]))).ravel()
    _, inv, cnt = np.unique(v, return_inverse=True, return_counts=True)
    return inv.astype(np.int64), cnt

def hgo(cnt):
    return float((cnt * np.log2(cnt)).sum() / NGEN)

def main():
    print("=" * 78)
    print("EXPERIMENT 24  the crossing conjecture across all observation partitions")
    print("=" * 78)
    print(f"  |H| = {NGEN:,} generators ({math.log2(NGEN):.2f} bits), "
          f"{len(PARTS)} partitions of 4 states, 16 reset sets each, T = {T}\n")
    rows = []
    for part in PARTS:
        k = len(part)
        obsv = block_of(part)
        if k == 1:
            continue                       # observation is constant; nothing to learn
        rid = []
        for j in range(NS):
            pm = np.zeros(NS, dtype=np.int64); pm[j] = 4
            rid.append(ids_of(laws(pm, obsv, k))[0])
        lab = "|".join("".join(map(str, b)) for b in part)
        for r in range(1, NS + 1):
            for A in itertools.combinations(range(NS), r):
                comb = np.stack([rid[j] for j in A], axis=1)
                _, cnt = ids_of(comb)
                cov = len({obsv[x] for x in A})
                cross = sum(1 for a, b in itertools.combinations(A, 2)
                            if obsv[a] != obsv[b])
                rows.append((lab, k, len(A), cov, cross, hgo(cnt)))
    # ---- X1 ----
    cells = defaultdict(lambda: defaultdict(list))
    for lab, k, sz, cov, cross, H in rows:
        cells[(lab, sz)][cov].append(H)
    ok = tot = 0
    viol = []
    for (lab, sz), by in cells.items():
        if len(by) < 2:
            continue
        tot += 1
        means = [(c, sum(v)/len(v)) for c, v in sorted(by.items())]
        if all(b[1] <= a[1] + 1e-12 for a, b in zip(means, means[1:])):
            ok += 1
        else:
            viol.append((lab, sz, means))
    frac = ok/tot if tot else 0.0
    print(f"  X1  mean H non-increasing in coverage, within (partition, |A|) cells:")
    print(f"      {ok}/{tot} cells = {frac:.1%}   "
          f"{'SUPPORTED' if frac >= 0.90 else 'NOT SUPPORTED'} (threshold 90%)")
    for lab, sz, means in viol[:4]:
        print(f"        exception: partition {lab}, |A|={sz}, "
              + ", ".join(f"cov{c}:{m:.3f}" for c, m in means))

    # ---- X2: uniform-block partitions only ----
    uni = {"|".join("".join(map(str, b)) for b in p)
           for p in PARTS if len({len(b) for b in p}) == 1 and len(p) > 1}
    ok2 = tot2 = 0
    for (lab, sz), by in cells.items():
        if lab not in uni or len(by) < 2:
            continue
        tot2 += 1
        means = [(c, sum(v)/len(v)) for c, v in sorted(by.items())]
        if all(b[1] <= a[1] + 1e-12 for a, b in zip(means, means[1:])):
            ok2 += 1
    print(f"\n  X2  same test, equal-block partitions only: {ok2}/{tot2}"
          f" = {ok2/tot2:.1%}" if tot2 else "\n  X2  no equal-block cells")

    # ---- X3 ----
    best = 0.0; where = None
    for (lab, sz), by in cells.items():
        if len(by) < 2: continue
        means = [sum(v)/len(v) for v in by.values()]
        if max(means)-min(means) > best:
            best = max(means)-min(means); where = (lab, sz)
    print(f"  X3  largest H spread across coverage at fixed |A|: {best:.3f} bits"
          f" (partition {where[0]}, |A|={where[1]})   "
          f"{'PASS' if best > 0.10 else 'FAIL'}")

    # ---- the 2+2 case, for comparison with Experiment 22 ----
    print(f"\n  the balanced 2+2 partitions in this coarser family:")
    print(f"    {'partition':<12}{'A':<10}{'cov':>5}{'H(G|O_A)':>11}")
    for lab, k, sz, cov, cross, H in rows:
        if k == 2 and lab.count("|") == 1 and len(lab) == 5 and sz == 2:
            A = "?"
            print(f"    {lab:<12}{'|A|=2':<10}{cov:>5}{H:>11.3f}")
    print("\n  (H is in bits of generator identity left unresolved; lower is better.)")

main()
