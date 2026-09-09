#!/usr/bin/env python3
"""
EXPERIMENT 35 -- the semantic kernel in belief space.

Rows of B_T are probability laws over words, so B_T 1 = 1 and therefore

    v B_T = 0  =>  0 = v B_T 1 = v.1,

i.e. K_T is contained in 1-perp. Every semantic-kernel direction has zero
coordinate sum and so is a probability-PRESERVING redistribution of hidden
belief. For an interior belief mu and small enough eps, mu + eps v lies in the
simplex and induces exactly the same observable law. The kernel is therefore
the nullspace of hidden-belief observability, and

    d_obs = (n - 1) - dim K_inf

counts the identifiable degrees of freedom of the hidden belief.

DRV  K_T is contained in 1-perp, for every case (forced; checks the code).
EXH  An explicit pair of distinct beliefs mu != nu in the simplex with
     mu B = nu B, in a system where no two PURE states are equivalent.
"""
import itertools
from fractions import Fraction as F

NS = 4
ROWS = [r for r in itertools.product(range(3), repeat=NS) if sum(r) == 2]
NR = len(ROWS)

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

def law_row(M, x, obsv, k, steps):
    out = []
    def rec(a, d):
        if d == steps:
            out.append(sum(a)); return
        nxt = [sum(a[u]*M[u][v] for u in range(NS)) for v in range(NS)]
        for y in range(k):
            rec([nxt[v] if obsv[v] == y else 0 for v in range(NS)], d+1)
    a0 = [F(0)]*NS; a0[x] = F(1)
    for y in range(k):
        rec([a0[v] if obsv[v] == y else F(0) for v in range(NS)], 1)
    return out

def left_kernel(B):
    n, m = len(B), len(B[0])
    aug = [[F(v) for v in B[i]] + [F(1 if j == i else 0) for j in range(n)]
           for i in range(n)]
    r = 0
    for c in range(m):
        piv = next((i for i in range(r, n) if aug[i][c] != 0), None)
        if piv is None: continue
        aug[r], aug[piv] = aug[piv], aug[r]
        pv = aug[r][c]; aug[r] = [v/pv for v in aug[r]]
        for i in range(n):
            if i != r and aug[i][c] != 0:
                f = aug[i][c]; aug[i] = [a-f*b for a,b in zip(aug[i], aug[r])]
        r += 1
        if r == n: break
    return [row[m:] for row in aug[r:]]

def main():
    print("=" * 78)
    print("EXPERIMENT 35  the semantic kernel is a belief-space nullspace")
    print("=" * 78)
    gens = [tuple(g) for g in itertools.product(range(NR), repeat=NS)][::311]
    bad = tot = 0
    for part in PARTS:
        k = len(part); obsv = [0]*NS
        for bi, b in enumerate(part):
            for x in b: obsv[x] = bi
        for g in gens:
            M = [[F(v, 2) for v in ROWS[i]] for i in g]
            B = [law_row(M, x, obsv, k, 5) for x in range(NS)]
            for v in left_kernel(B):
                tot += 1
                if sum(v) != 0: bad += 1
    print(f"\n  DRV  every kernel vector has zero coordinate sum:")
    print(f"       {tot} kernel vectors checked, {bad} with nonzero sum   "
          f"{'PASS' if bad == 0 else 'FAIL'}")

    # ---- the explicit belief witness ----
    print(f"\n  EXH  two distinct beliefs with identical observable laws,")
    print(f"       in a system where no two PURE states are equivalent")
    M = [[F(0),F(1),F(0),F(0)], [F(0),F(1),F(0),F(0)],
         [F(1),F(0),F(0),F(0)], [F(1,2),F(1,2),F(0),F(0)]]
    obsv = [0, 1, 0, 0]          # partition {0,2,3} | {1}
    B = [law_row(M, x, obsv, 2, 5) for x in range(NS)]
    pure_equal = [(x, y) for x in range(NS) for y in range(x+1, NS) if B[x] == B[y]]
    K = left_kernel(B)
    print(f"       pure states with identical laws: "
          f"{pure_equal if pure_equal else 'none'}")
    print(f"       dim K = {len(K)}, basis {[ [str(c) for c in v] for v in K ]}")
    v = K[0]
    nu = [F(1,4)]*NS
    eps = F(1,4)
    mu = [nu[i] + eps*v[i] for i in range(NS)]
    print(f"       nu = {[str(c) for c in nu]}")
    print(f"       mu = {[str(c) for c in mu]}   (= nu + (1/4) v)")
    print(f"       both in the simplex: "
          f"{all(0 <= c <= 1 for c in mu) and sum(mu) == 1}")
    lm = [sum(mu[x]*B[x][j] for x in range(NS)) for j in range(len(B[0]))]
    ln = [sum(nu[x]*B[x][j] for x in range(NS)) for j in range(len(B[0]))]
    print(f"       mu and nu are distinct: {mu != nu}")
    print(f"       observable laws identical: {lm == ln}")
    print(f"\n       d_obs = (n-1) - dim K = {NS-1} - {len(K)} = {NS-1-len(K)}")
    print("       so of three hidden-belief degrees of freedom, "
          f"{NS-1-len(K)} are observable.")

main()
