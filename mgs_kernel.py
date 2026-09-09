#!/usr/bin/env python3
"""
EXPERIMENT 31 -- the semantic kernel.

Frozen before running.

Experiments 27-30 produced three mechanisms for experimental redundancy, each
with its own criterion. They are one equation. Let b_x be the observable law
of the next T-1 symbols starting from hidden state x, and let B be the matrix
with those rows. Resetting to s and taking one transition leaves the hidden
distribution P_s, so the observable SUFFIX law is exactly

    S_s^(T) = P_s B_{O,T-1}(G),

and therefore

    S_s = S_t   <=>   (P_s - P_t) B = 0   <=>   P_s - P_t in ker B.

Identity ties are the case P_s - P_t = 0, which lies in every kernel. A
permutation tie survives when P_s(I - Pi)B = 0, which is weaker than requiring
Pi B = B pointwise -- explaining the 358 cases of Experiment 30 that agreed
without behavioural compatibility, as weighted cancellation rather than
coincidence. Lumpability makes ker B non-trivial, so distinct rows can have a
difference the observation is blind to.

CRITERIA
  K1  IMPLEMENTATION CHECK, not a discovery. For every generator, partition
      and reset pair, S_s = S_t iff (P_s - P_t)B = 0. Zero failures; the
      proposition is a one-line proof and a failure means our code is wrong.
  K2  dim ker B is 0 under the discrete partition (every state its own block)
      and positive under the balanced partition for the lumpable class.
  K3  The kernel contracts with horizon: dim ker B_{O,T} is non-increasing
      in T, for every generator and partition tested.

Exact rational arithmetic throughout; no rank tolerances.
"""
import itertools, math
from fractions import Fraction as F
from collections import Counter, defaultdict

NS = 4
ROWS = [r for r in itertools.product(range(3), repeat=NS) if sum(r) == 2]
NR = len(ROWS)
PERM = (1, 0, 3, 2)
PIDX = {r: i for i, r in enumerate(ROWS)}
PERM_OF = [PIDX[tuple(ROWS[i][p] for p in PERM)] for i in range(NR)]
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
    """Observable law of `steps` symbols from state x, as an integer vector."""
    if steps == 0:
        return (1,)
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

def rank_exact(rows):
    """Exact rank over Q of a list of integer row vectors."""
    m = [[F(v) for v in r] for r in rows]
    R = len(m); Cn = len(m[0]) if R else 0
    r = 0
    for c in range(Cn):
        piv = next((i for i in range(r, R) if m[i][c] != 0), None)
        if piv is None: continue
        m[r], m[piv] = m[piv], m[r]
        pv = m[r][c]
        m[r] = [v / pv for v in m[r]]
        for i in range(R):
            if i != r and m[i][c] != 0:
                f = m[i][c]
                m[i] = [a - f * b for a, b in zip(m[i], m[r])]
        r += 1
        if r == R: break
    return r

def build(kind):
    trip = list(itertools.product(range(NR), repeat=3))
    if kind == "FREE":
        return [tuple(g) for g in itertools.product(range(NR), repeat=NS)][::311]
    if kind == "EQUAL":
        return [(a, b, c, a) for a, b, c in trip][::31]
    if kind == "PERM":
        return [(a, b, c, PERM_OF[a]) for a, b, c in trip][::31]
    if kind == "LUMPABLE":
        L = [(r0, r1, r2, r3) for a in (0,1,2) for b in (0,1,2)
             for r0 in BYPROF[a] for r1 in BYPROF[a]
             for r2 in BYPROF[b] for r3 in BYPROF[b]]
        return L[::37]
    raise ValueError(kind)

CLASSES = ["FREE", "EQUAL", "PERM", "LUMPABLE"]

def main():
    print("=" * 78)
    print("EXPERIMENT 31  the semantic kernel: (P_s - P_t) B = 0")
    print("=" * 78)
    T = 5
    k1_bad = k1_tot = 0
    dimker = defaultdict(list)
    for kind in CLASSES:
        gens = build(kind)
        for part in PARTS:
            k = len(part)
            obsv = [0] * NS
            for bi, b in enumerate(part):
                for x in b: obsv[x] = bi
            lab = "|".join("".join(map(str, b)) for b in part)
            for g in gens:
                M = [list(ROWS[i]) for i in g]
                B = [law_row(M, x, obsv, k, T - 1) for x in range(NS)]
                dimker[(kind, lab)].append(NS - rank_exact(B))
                for s, t in itertools.combinations(range(NS), 2):
                    diff = [M[s][v] - M[t][v] for v in range(NS)]
                    prod = [sum(diff[v] * B[v][c] for v in range(NS))
                            for c in range(len(B[0]))]
                    inker = all(v == 0 for v in prod)
                    Ss = tuple(sum(M[s][v] * B[v][c] for v in range(NS))
                               for c in range(len(B[0])))
                    St = tuple(sum(M[t][v] * B[v][c] for v in range(NS))
                               for c in range(len(B[0])))
                    k1_tot += 1
                    if (Ss == St) != inker:
                        k1_bad += 1
    print(f"\n  K1  S_s = S_t iff (P_s - P_t)B = 0")
    print(f"      {k1_tot:,} reset pairs tested, {k1_bad} failures   "
          f"{'PASS' if k1_bad == 0 else 'FAIL'}")
    print("      (a one-line proof; this checks the implementation)")

    print(f"\n  dim ker B by class and observation partition (mean over generators):")
    print(f"  {'partition':<12}" + "".join(f"{c:>11}" for c in CLASSES))
    for part in PARTS:
        lab = "|".join("".join(map(str, b)) for b in part)
        row = []
        for c in CLASSES:
            v = dimker[(c, lab)]
            row.append(sum(v) / len(v) if v else float('nan'))
        print(f"  {lab:<12}" + "".join(f"{x:>11.3f}" for x in row))

    disc = "0|1|2|3"
    bal = "01|23"
    k2a = all(v == 0 for v in dimker[("LUMPABLE", disc)])
    k2b = all(v > 0 for v in dimker[("LUMPABLE", bal)])
    print(f"\n  K2  dim ker = 0 under the discrete partition (lumpable class): "
          f"{'yes' if k2a else 'no'}")
    print(f"      dim ker > 0 under the balanced partition (lumpable class): "
          f"{'yes' if k2b else 'no'}   {'PASS' if k2a and k2b else 'FAIL'}")

    print(f"\n  K3  does the kernel contract with horizon?")
    print(f"      {'partition':<12}" + "".join(f"{'T='+str(t):>8}" for t in range(2, 7)))
    bad3 = 0
    for kind in ("LUMPABLE", "FREE"):
        gens = build(kind)[:6]
        for part in PARTS[:6]:
            k = len(part)
            obsv = [0] * NS
            for bi, b in enumerate(part):
                for x in b: obsv[x] = bi
            lab = "|".join("".join(map(str, b)) for b in part)
            means = []
            for t in range(2, 7):
                ds = [NS - rank_exact([law_row([list(ROWS[i]) for i in g], x, obsv, k, t - 1)
                                       for x in range(NS)]) for g in gens]
                means.append(sum(ds) / len(ds))
            for a, b2 in zip(means, means[1:]):
                if b2 > a + 1e-12: bad3 += 1
            if kind == "LUMPABLE":
                print(f"      {lab:<12}" + "".join(f"{m:>8.2f}" for m in means))
    print(f"      non-monotone steps across both classes: {bad3}   "
          f"{'PASS' if bad3 == 0 else 'FAIL'}")

main()
