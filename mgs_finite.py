#!/usr/bin/env python3
"""
EXPERIMENT 33 -- implementation check of finite determination.

NOT an empirical experiment. Proposition 10 states that for an n-state chain
under a deterministic observation with k non-empty blocks,

    K_{T+1} = K_T  =>  K_s = K_T for all s >= T,      and   K_inf = K_{n-k+1}.

Proof sketch. With E_y the diagonal selector for symbol y and
q_w = E_{y_0} P E_{y_1} P ... P E_{y_{T-1}} 1, the columns of B_T are the q_w,
so K_T = V_T^perp for V_T = span{q_w : |w| = T}. Since q_{yw} = E_y P q_w we
have V_{T+1} = span_y E_y P V_T, so V_T = V_{T-1} forces V_{T+1} = V_T and the
chain stabilises. dim V_1 = k because the E_y 1 have disjoint supports, and
each pre-stabilisation step strictly increases dimension, leaving at most
n - k increases.

Everything below is therefore a DRV check: a failure means the code is wrong,
not that the mathematics is. What it guards is the claim now made in the
paper, that a mixture-only kernel vector found at the bounded horizon is
PERMANENTLY invisible rather than merely invisible so far.
"""
import itertools
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

def kernel_dim(M, obsv, k, T):
    B = [law_row(M, x, obsv, k, T) for x in range(NS)]
    m = [[F(v) for v in r] for r in B]
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
        if r == R: break
    return NS - r

def main():
    print("=" * 78)
    print("EXPERIMENT 33  finite determination of the belief-observability kernel (DRV check)")
    print("=" * 78)
    free = [tuple(g) for g in itertools.product(range(NR), repeat=NS)][::311]
    lump = [(r0, r1, r2, r3) for a in (0,1,2) for b in (0,1,2)
            for r0 in BYPROF[a] for r1 in BYPROF[a]
            for r2 in BYPROF[b] for r3 in BYPROF[b]][::37]
    gens = free + lump
    TMAX = 8
    over_bound = plateau_broken = 0
    tsem = defaultdict(Counter)
    for part in PARTS:
        k = len(part); obsv = [0] * NS
        for bi, b in enumerate(part):
            for x in b: obsv[x] = bi
        bound = NS - k + 1
        for g in gens:
            M = [list(ROWS[i]) for i in g]
            dims = [kernel_dim(M, obsv, k, T) for T in range(1, TMAX + 1)]
            # first T with K_T = K_{T+1}; equal dimension implies equal subspace
            t = next((T for T in range(1, TMAX) if dims[T-1] == dims[T]), None)
            if t is None or t > bound:
                over_bound += 1
            else:
                tsem[k][t] += 1
            # once equal, must stay equal
            if t is not None and any(d != dims[t-1] for d in dims[t-1:]):
                plateau_broken += 1
    print(f"\n  checked {len(gens)} generators x {len(PARTS)} partitions "
          f"= {len(gens)*len(PARTS):,} cases, horizons 1..{TMAX}\n")
    print(f"  T_sem exceeds the bound n-k+1:        {over_bound}   "
          f"{'PASS' if over_bound == 0 else 'FAIL'}")
    print(f"  a plateau later broken by a new drop: {plateau_broken}   "
          f"{'PASS' if plateau_broken == 0 else 'FAIL'}")
    print(f"\n  observed T_sem by number of observable blocks k:")
    print(f"    {'k':>3}{'bound n-k+1':>14}   distribution of T_sem")
    for k in sorted(tsem):
        d = ", ".join(f"T={t}: {c}" for t, c in sorted(tsem[k].items()))
        print(f"    {k:>3}{NS-k+1:>14}   {d}")
    print("\n  For k = n the bound is 1, agreeing with the injective-observation")
    print("  corollary. For the balanced two-block case the bound is 3, so a kernel")
    print("  vector exhibited at T = 3 is certified permanently invisible.")

main()
