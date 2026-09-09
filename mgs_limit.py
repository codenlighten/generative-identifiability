#!/usr/bin/env python3
"""
EXPERIMENT 32 -- the limiting semantic kernel.

Frozen before running.

Three of Experiment 31's findings turned out to be theorems: kernel
contraction with horizon (marginalisation), zero kernel under injective
observation (disjoint row supports), and dim ker = n - k for a strongly
lumpable k-block partition (identical within-block rows). None of them is an
empirical regularity, and this experiment does not re-measure them except as
implementation checks.

What is NOT settled by those theorems is the unconstrained case. Kernels form
a descending chain K_1 >= K_2 >= ... in an n-dimensional space, so the
dimension can strictly drop at most n times, but a single plateau does not
prove stabilisation forever. Write

    K_inf = intersection over T of ker B_{O,T}

for the structurally blind directions, and d_inf = dim K_inf.

The question this experiment asks is whether K_inf is richer than a state
partition. Let x = y mean states x and y induce identical observable laws --
classical behavioural equivalence -- and let D = span{ e_x - e_y : x = y }.
Then D is contained in K always. If K is strictly larger for some generator,
the semantic kernel detects linear observational equivalences among MIXTURES
of hidden states that no equivalence relation on states can express.

CRITERIA
  M1  IMPLEMENTATION CHECKS of the three theorems: ker B_{T+1} subset ker B_T;
      injective O gives the zero kernel; strong lumpability gives dim = n - k
      with the within-block difference basis. Zero failures.
  M2  DESCRIPTIVE. How often is dim K > dim D, i.e. how often does the kernel
      contain mixture-only directions? No direction pre-committed.
  M3  DESCRIPTIVE. Do two generators share a state-equivalence partition while
      differing in K? That would show K is not a function of the partition.

Exact rational arithmetic; kernels are computed as bases, not dimensions.
"""
import itertools, math
from fractions import Fraction as F
from collections import Counter, defaultdict

NS = 4
ROWS = [r for r in itertools.product(range(3), repeat=NS) if sum(r) == 2]
NR = len(ROWS)
PROF = [r[0] + r[1] for r in ROWS]
BYPROF = {a: [i for i in range(NR) if PROF[i] == a] for a in (0, 1, 2)}

def partitions(items):
    if len(items) == 1:
        yield [items]; return
    first, rest = items[0], items[1:]
    for p in partitions(rest):
        for i in range(len(p)):
            yield p[:i] + [[first] + p[i]] + p[i+1:]
        yield [[first]] + p
PARTS = [sorted([sorted(b) for b in p]) for p in partitions(list(range(NS)))]
PARTS = [p for p in PARTS if len(p) > 1]
PARTS.sort(key=lambda p: (len(p), [-len(b) for b in p]))

def law_row(M, x, obsv, k, steps):
    out = []
    def rec(a, d):
        if d == steps:
            out.append(sum(a)); return
        nxt = [sum(a[u] * M[u][v] for u in range(NS)) for v in range(NS)]
        for y in range(k):
            rec([nxt[v] if obsv[v] == y else 0 for v in range(NS)], d + 1)
    a0 = [0] * NS; a0[x] = 1
    for y in range(k):
        rec([a0[v] if obsv[v] == y else 0 for v in range(NS)], 1)
    return tuple(out)

def left_kernel(B):
    """Exact basis of { v : vB = 0 }, by reducing [B | I] and reading off the
    augmented rows whose B-part vanished."""
    n = len(B); m = len(B[0])
    aug = [[F(v) for v in B[i]] + [F(1 if j == i else 0) for j in range(n)]
           for i in range(n)]
    r = 0
    for c in range(m):
        piv = next((i for i in range(r, n) if aug[i][c] != 0), None)
        if piv is None: continue
        aug[r], aug[piv] = aug[piv], aug[r]
        pv = aug[r][c]
        aug[r] = [v / pv for v in aug[r]]
        for i in range(n):
            if i != r and aug[i][c] != 0:
                f = aug[i][c]
                aug[i] = [a - f * b for a, b in zip(aug[i], aug[r])]
        r += 1
        if r == n: break
    return [row[m:] for row in aug[r:]]

def span_rank(vecs):
    if not vecs: return 0
    m = [[F(v) for v in row] for row in vecs]
    R, Cn, r = len(m), len(m[0]), 0
    for c in range(Cn):
        piv = next((i for i in range(r, R) if m[i][c] != 0), None)
        if piv is None: continue
        m[r], m[piv] = m[piv], m[r]
        pv = m[r][c]; m[r] = [v / pv for v in m[r]]
        for i in range(R):
            if i != r and m[i][c] != 0:
                f = m[i][c]
                m[i] = [a - f * b for a, b in zip(m[i], m[r])]
        r += 1
    return r

def main():
    print("=" * 78)
    print("EXPERIMENT 32  the limiting semantic kernel")
    print("=" * 78)
    TMAX = 7
    free = [tuple(g) for g in itertools.product(range(NR), repeat=NS)][::311]
    lump = [(r0, r1, r2, r3) for a in (0,1,2) for b in (0,1,2)
            for r0 in BYPROF[a] for r1 in BYPROF[a]
            for r2 in BYPROF[b] for r3 in BYPROF[b]][::37]

    # ---------- M1 ----------
    bad_mono = bad_inj = bad_lump = 0
    for part in PARTS:
        k = len(part); obsv = [0]*NS
        for bi, b in enumerate(part):
            for x in b: obsv[x] = bi
        for g in free:
            M = [list(ROWS[i]) for i in g]
            dims = []
            for T in range(1, TMAX + 1):
                Ks = left_kernel([law_row(M, x, obsv, k, T) for x in range(NS)])
                dims.append(len(Ks))
            for a, b2 in zip(dims, dims[1:]):
                if b2 > a: bad_mono += 1
            if k == NS and dims[-1] != 0: bad_inj += 1
    for part in PARTS:
        k = len(part); obsv = [0]*NS
        for bi, b in enumerate(part):
            for x in b: obsv[x] = bi
        if any(len({tuple(sum(list(ROWS[g[s]])[u] for u in t) for t in part)
                    for s in blk}) > 1 for g in lump for blk in part):
            continue
        for g in lump:
            M = [list(ROWS[i]) for i in g]
            Ks = left_kernel([law_row(M, x, obsv, k, TMAX) for x in range(NS)])
            if len(Ks) != NS - k: bad_lump += 1
    print(f"\n  M1  implementation checks of the three theorems")
    print(f"      kernel monotone in horizon: {bad_mono} failures")
    print(f"      injective O gives zero kernel: {bad_inj} failures")
    print(f"      lumpable gives dim = n - k: {bad_lump} failures")
    print(f"      {'PASS' if bad_mono == bad_inj == bad_lump == 0 else 'FAIL'}")

    # ---------- M2, M3 ----------
    richer = total = 0
    plateau = Counter()
    bypart = defaultdict(set)
    for part in PARTS:
        k = len(part); obsv = [0]*NS
        for bi, b in enumerate(part):
            for x in b: obsv[x] = bi
        lab = "|".join("".join(map(str, b)) for b in part)
        for g in free:
            M = [list(ROWS[i]) for i in g]
            B = [law_row(M, x, obsv, k, TMAX) for x in range(NS)]
            K = left_kernel(B)
            # behavioural equivalence classes of states at TMAX
            groups = defaultdict(list)
            for x in range(NS):
                groups[B[x]].append(x)
            D = [[1 if u == blk[0] else (-1 if u == y else 0) for u in range(NS)]
                 for blk in groups.values() for y in blk[1:]]
            dK, dD = len(K), span_rank(D) if D else 0
            total += 1
            if dK > dD: richer += 1
            bypart[(lab, tuple(sorted(tuple(sorted(v)) for v in groups.values())))].add(dK)
            dims = [len(left_kernel([law_row(M, x, obsv, k, T) for x in range(NS)]))
                    for T in range(1, TMAX + 1)]
            stab = next((T for T in range(1, TMAX)
                         if all(d == dims[T-1] for d in dims[T-1:])), TMAX)
            plateau[stab] += 1
    print(f"\n  M2  kernel richer than the state partition (dim K > dim D):")
    print(f"      {richer}/{total} generator/partition pairs = {richer/total:.1%}")
    print("      (D is the span of within-equivalence-class state differences;")
    print("       D is contained in K always, so any excess is a mixture-only direction)")
    varying = {key: ds for key, ds in bypart.items() if len(ds) > 1}
    print(f"\n  M3  same observation partition and same state-equivalence classes,")
    print(f"      but different kernel dimension: {len(varying)} cases")
    for key, ds in list(varying.items())[:3]:
        print(f"        partition {key[0]}, classes {key[1]} -> dims {sorted(ds)}")
    print(f"\n  horizon at which dim ker stabilises (FREE class, all partitions):")
    for T in sorted(plateau):
        print(f"      T = {T}: {plateau[T]} pairs")

main()
