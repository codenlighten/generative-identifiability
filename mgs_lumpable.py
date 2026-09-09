#!/usr/bin/env python3
"""
EXPERIMENT 29 -- can the observation map CREATE a dependency?

Frozen before running.

Experiment 28 found zero observation-created dependencies: no j ever lay in
cl_Z(A) without lying in cl_G(A). That is not a theorem. For general
deterministic experiment maps it is false -- with R_0, R_1 independent and
Z_0 = Z_1 = R_0 XOR R_1, we get 1 in cl_Z({0}) while 1 is not in cl_G({0}).
The question is whether the reset-law architecture can exhibit it, and the
answer should be yes by way of lumpability: microscopic transitions that the
quotient observation erases exactly.

CONSTRUCTION. Fix the balanced partition {0,1} | {2,3}. Constrain each row so
that its BLOCK transition sums depend only on which block the source state is
in -- states 0 and 1 both send mass a to block 0 and 2-a to block 1, states 2
and 3 both send b and 2-b -- while the split WITHIN each block stays free.
Then R_0 and R_1 are not functions of one another, yet the observed process is
Markov on blocks and reset 0 and reset 1 induce the same observable law.

CLASSES (four hidden states, rows from the same 10-row set)
  FREE            unconstrained, 10^4
  TIED            R3 = R0                       (generator-induced redundancy)
  LUMPABLE        block sums tied by block       (observation-induced, target)
  TIED_LUMPABLE   lumpable and R3 = R2           (both mechanisms)
  RANDOM          uniform subset matched to LUMPABLE's size, fixed seed

DECOMPOSITION. With FD_G = {(A,j) : j in cl_G(A) \\ A} and FD_Z likewise,
  preserved  P = FD_G and FD_Z      destroyed D = FD_G minus FD_Z
  created    C = FD_Z minus FD_G
Each class is characterised by (|P|, |D|, |C|).

CRITERIA
  L1  The LUMPABLE class has |C| > 0: some reset output is determined by
      others without the corresponding generator coordinates being so.
  L2  The size-matched RANDOM class has |C| = 0, so the effect is the
      construction and not the cardinality.
"""
import numpy as np, math, itertools

NS, T = 4, 5
ROWS = [r for r in itertools.product(range(3), repeat=NS) if sum(r) == 2]
NR = len(ROWS)
ROWA = np.array(ROWS, dtype=np.int64)
PROF = [r[0] + r[1] for r in ROWS]                 # mass sent to block {0,1}
BYPROF = {a: [i for i in range(NR) if PROF[i] == a] for a in (0, 1, 2)}

def partitions(items):
    if len(items) == 1:
        yield [items]; return
    first, rest = items[0], items[1:]
    for p in partitions(rest):
        for i in range(len(p)):
            yield p[:i] + [[first] + p[i]] + p[i+1:]
        yield [[first]] + p

PARTS = [sorted([sorted(b) for b in p]) for p in partitions([0, 1, 2, 3])]
PARTS = [p for p in PARTS if len(p) > 1]
PARTS.sort(key=lambda p: (len(p), [-len(b) for b in p]))
BALANCED = [0, 1], [2, 3]

def build():
    free = np.array(list(itertools.product(range(NR), repeat=NS)), dtype=np.int64)
    tied = np.array([[a, b, c, a] for a in range(NR) for b in range(NR)
                     for c in range(NR)], dtype=np.int64)
    lump = np.array([[r0, r1, r2, r3]
                     for a in (0, 1, 2) for b in (0, 1, 2)
                     for r0 in BYPROF[a] for r1 in BYPROF[a]
                     for r2 in BYPROF[b] for r3 in BYPROF[b]], dtype=np.int64)
    tl = np.array([[r0, r1, r2, r2]
                   for a in (0, 1, 2) for b in (0, 1, 2)
                   for r0 in BYPROF[a] for r1 in BYPROF[a]
                   for r2 in BYPROF[b]], dtype=np.int64)
    rng = np.random.default_rng(29)
    rand = free[np.sort(rng.choice(len(free), size=len(lump), replace=False))]
    return [("FREE", free), ("TIED", tied), ("LUMPABLE", lump),
            ("TIED_LUMPABLE", tl), ("RANDOM", rand)]

def ncls(cols):
    if not cols: return 1
    comb = np.stack(cols, axis=1).astype(np.int64)
    v = np.ascontiguousarray(comb).view(
        np.dtype((np.void, comb.dtype.itemsize * comb.shape[1]))).ravel()
    return len(np.unique(v))

def zids(rowidx, obsv, k, j):
    n = len(rowidx)
    M = np.empty((n, NS, NS), dtype=np.int64)
    for pos in range(NS):
        M[:, pos, :] = ROWA[rowidx[:, pos]]
    out = np.empty((n, k ** T), dtype=np.int32)
    col = [0]
    def rec(alpha, depth):
        if depth == T:
            out[:, col[0]] = alpha.sum(axis=1); col[0] += 1; return
        nxt = np.matmul(alpha[:, None, :], M)[:, 0, :]
        for y in range(k):
            rec(np.where(obsv[None, :] == y, nxt, 0), depth + 1)
    pm = np.zeros(NS, dtype=np.int64); pm[j] = 4
    for y0 in range(k):
        rec(np.broadcast_to(np.where(obsv == y0, pm, 0).astype(np.int64),
                            (n, NS)).copy(), 1)
    v = np.ascontiguousarray(out).view(
        np.dtype((np.void, out.dtype.itemsize * out.shape[1]))).ravel()
    _, inv = np.unique(v, return_inverse=True)
    return inv.astype(np.int64)

SUBS = [frozenset(A) for r in range(NS + 1)
        for A in itertools.combinations(range(NS), r)]

def fdset(colfn):
    fd = set()
    for A in SUBS:
        base = ncls([colfn(j) for j in sorted(A)])
        for j in range(NS):
            if j in A: continue
            if ncls([colfn(i) for i in sorted(A | {j})]) == base:
                fd.add((A, j))
    return fd

def main():
    print("=" * 78)
    print("EXPERIMENT 29  can the observation map create a dependency?")
    print("=" * 78)
    res = {}
    for name, rowidx in build():
        n = len(rowidx)
        FDG = fdset(lambda j: rowidx[:, j])
        # the balanced partition the lumpable class was built for
        obsv = np.array([0, 0, 1, 1], dtype=np.int64)
        zid = {j: zids(rowidx, obsv, 2, j) for j in range(NS)}
        FDZ = fdset(lambda j: zid[j])
        P, D, C = FDG & FDZ, FDG - FDZ, FDZ - FDG
        res[name] = (len(P), len(D), len(C))
        print(f"\n  {name:<14} |H| = {n:>6,}    "
              f"preserved {len(P):>3}   destroyed {len(D):>3}   created {len(C):>3}")
        if C:
            ex = sorted(C, key=lambda t: (len(t[0]), sorted(t[0]), t[1]))[:3]
            for A, j in ex:
                print(f"      created: Z_{sorted(A)} determines Z_{j}, "
                      f"but R_{sorted(A)} does not determine R_{j}")
        # how the effect depends on which partition is used
        made = []
        for part in PARTS:
            k = len(part)
            ov = np.empty(NS, dtype=np.int64)
            for bi, b in enumerate(part):
                for x in b: ov[x] = bi
            zz = {j: zids(rowidx, ov, k, j) for j in range(NS)}
            c = len(fdset(lambda j: zz[j]) - FDG)
            if c: made.append(("|".join("".join(map(str, b)) for b in part), c))
        print(f"      partitions with created dependencies: {len(made)}/{len(PARTS)}"
              + (f"   {', '.join(f'{l}({c})' for l, c in made[:4])}" if made else ""))
        # reconcile with Experiment 28, which aggregated over all partitions
        tot_c = tot_d = 0
        for part in PARTS:
            k = len(part)
            ov = np.empty(NS, dtype=np.int64)
            for bi, b in enumerate(part):
                for x in b: ov[x] = bi
            zz = {j: zids(rowidx, ov, k, j) for j in range(NS)}
            fz = fdset(lambda j: zz[j])
            tot_c += len(fz - FDG); tot_d += len(FDG - fz)
        print(f"      summed over all {len(PARTS)} partitions: "
              f"destroyed {tot_d}, created {tot_c}")

    print()
    print("=" * 78)
    print("PRE-REGISTERED CRITERIA")
    print("=" * 78)
    l1 = res["LUMPABLE"][2] > 0
    l2 = res["RANDOM"][2] == 0
    print(f"  L1  LUMPABLE class creates dependencies      "
          f"{'PASS' if l1 else 'FAIL'}   (|C| = {res['LUMPABLE'][2]})")
    print(f"  L2  size-matched RANDOM class does not       "
          f"{'PASS' if l2 else 'FAIL'}   (|C| = {res['RANDOM'][2]})")
    print(f"\n  {'class':<16}{'|P|':>6}{'|D|':>6}{'|C|':>6}")
    for k2, v in res.items():
        print(f"  {k2:<16}{v[0]:>6}{v[1]:>6}{v[2]:>6}")

main()
