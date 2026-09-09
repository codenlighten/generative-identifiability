#!/usr/bin/env python3
"""
EXPERIMENT 34 -- the emission-rank bound.

Corollary 11 states that for a finite HMM with emission matrix C,

    K_inf = K_{n - rank(C) + 1},

generalising the deterministic case by replacing the block count k with
rank(C). The plateau argument needs only sum_y E_y = I and P1 = 1, which hold
for stochastic emissions; only the starting dimension changes, since E_y 1 is
the y-th column of C so dim V_1 = rank(C).

Two parts, and they are different kinds of claim:

  DRV  The bound holds. A failure means the code is wrong, not the theorem.
  EXH  The bound is TIGHT for each attainable emission rank -- that is, some
       (transition, emission) pair actually attains T_sem = n - rank(C) + 1.
       Experiment 33 established tightness for the deterministic case; this
       asks the same question of the generalisation, and it is a measurement.

Emission matrices are constructed with prescribed rank and their rank is then
verified rather than assumed.
"""
import itertools
from fractions import Fraction as F
from collections import Counter, defaultdict

NS = 4
ROWS = [r for r in itertools.product(range(3), repeat=NS) if sum(r) == 2]
NR = len(ROWS)
Q = F(1, 4)

EMITS = {
 "rank-1 uniform":    [[F(1,4)]*4 for _ in range(NS)],
 "rank-2 paired":     [[F(1,2),F(1,2),F(0),F(0)], [F(1,2),F(1,2),F(0),F(0)],
                       [F(0),F(0),F(1,2),F(1,2)], [F(0),F(0),F(1,2),F(1,2)]],
 "rank-3 mixed":      [[F(1,2),F(1,2),F(0),F(0)], [F(1,2),F(1,2),F(0),F(0)],
                       [F(0),F(0),F(1),F(0)],     [F(0),F(0),F(0),F(1)]],
 "rank-4 stochastic": [[F(3,4) if i==j else F(1,12) for j in range(4)]
                       for i in range(NS)],
}

def rank_of(rows):
    m = [list(r) for r in rows]
    R, Cn, r = len(m), len(m[0]), 0
    for c in range(Cn):
        piv = next((i for i in range(r, R) if m[i][c] != 0), None)
        if piv is None: continue
        m[r], m[piv] = m[piv], m[r]
        pv = m[r][c]; m[r] = [v/pv for v in m[r]]
        for i in range(R):
            if i != r and m[i][c] != 0:
                f = m[i][c]; m[i] = [a-f*b for a,b in zip(m[i], m[r])]
        r += 1
    return r

def kernel_dim(P, C, T):
    """dim of the left kernel of B_T, exactly."""
    m = len(C[0])
    cols = []
    # columns q_w = E_{y0} P E_{y1} P ... E_{y_{T-1}} 1, accumulated from the right
    def build(depth, vec):
        # vec is the column being accumulated from the right
        if depth == 0:
            cols.append(list(vec)); return
        for y in range(m):
            # left-multiply by P then by E_y  ->  (E_y P vec)_x = C[x][y] * sum_v P[x][v] vec[v]
            pv = [sum(P[x][v]*vec[v] for v in range(NS)) for x in range(NS)]
            build(depth-1, [C[x][y]*pv[x] for x in range(NS)])
    for y in range(m):
        build(T-1, [C[x][y] for x in range(NS)])
    B = [[cols[j][x] for j in range(len(cols))] for x in range(NS)]
    return NS - rank_of(B)

def main():
    print("=" * 78)
    print("EXPERIMENT 34  the emission-rank bound  K_inf = K_{n - rank(C) + 1}")
    print("=" * 78)
    gens = [tuple(g) for g in itertools.product(range(NR), repeat=NS)][::677]
    print(f"\n  {len(gens)} transition matrices x {len(EMITS)} emission matrices, n = {NS}\n")
    print(f"  {'emission':<20}{'rank C':>8}{'bound':>7}   observed T_sem")
    over = 0; tight = {}
    for name, C in EMITS.items():
        rC = rank_of([list(r) for r in C])
        bound = NS - rC + 1
        seen = Counter()
        for g in gens:
            P = [[F(v, 2) for v in ROWS[i]] for i in g]
            dims = [kernel_dim(P, C, T) for T in range(1, bound + 3)]
            t = next((T for T in range(1, len(dims)) if dims[T-1] == dims[T]), None)
            if t is None or t > bound: over += 1
            else: seen[t] += 1
        tight[name] = (rC, bound, max(seen) if seen else None)
        d = ", ".join(f"T={t}: {c}" for t, c in sorted(seen.items()))
        print(f"  {name:<20}{rC:>8}{bound:>7}   {d}")
    print(f"\n  DRV  T_sem exceeds n - rank(C) + 1:  {over}   "
          f"{'PASS' if over == 0 else 'FAIL'}")
    attained = [n for n,(rC,b,mx) in tight.items() if mx == b]
    print(f"  EXH  bound attained for {len(attained)}/{len(EMITS)} emission ranks: "
          + ", ".join(attained))
    for n,(rC,b,mx) in tight.items():
        if mx != b:
            print(f"       {n}: bound {b}, largest observed {mx} -- not attained")
    print("\n  The rank-1 gap is structural, not a sampling artefact. A row-stochastic")
    print("  C of rank 1 has all rows equal, since they are proportional and each sums")
    print("  to 1; so emission is independent of state, V_1 = span{1}, and")
    print("  V_2 = span_y E_y P 1 = span_y E_y 1 = V_1. Hence T_sem = 1 always, and")
    print("  the bound n - 1 + 1 = n cannot be attained for n > 1.")
    allsame = all(C[0] == C[x] for x in range(NS)) if True else None
    C1 = EMITS["rank-1 uniform"]
    print(f"  check: rank-1 matrix has all rows equal: "
          f"{all(C1[0] == C1[x] for x in range(NS))}")

main()
