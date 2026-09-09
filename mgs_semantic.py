#!/usr/bin/env python3
"""
EXPERIMENT 30 -- when does a generator dependency survive the semantic map?

Frozen before running.

Experiments 27-29 established that observation can preserve, destroy, or
create functional dependencies among reset outputs, and sorted four constraint
families into "robust" and "fragile". That is a taxonomy, not a criterion.
This experiment tests a proposed criterion.

Write the semantic map as Phi_{O,T} : G -> (Z_0^(T), ..., Z_{n-1}^(T)).
Decomposing a reset law by its first symbol,

    Z_s^(T)  =  delta_{O(s)}  tensor  sum_x P_sx Z_x^(T-1),

so whether two resets agree after the first symbol depends on the transition
row AND on what the destination states subsequently do.

  IDENTITY  R_i = R_j gives sum_x P_ix Z_x = sum_x P_jx Z_x term by term, for
            every O and every T. Robustness should be provable, so P1 is an
            implementation check of that argument rather than a discovery.

  PERMUTATION  R_j = pi(R_i) gives sum_y P_iy Z_{pi(y)}^(T-1) against
            sum_y P_iy Z_y^(T-1). These agree if the permutation moves each
            reachable destination to a BEHAVIOURALLY EQUIVALENT one, where
            s ~ t iff Z_s^(T-1) = Z_t^(T-1) for that generator. That is the
            criterion under test, in the sufficient direction only.

  LUMPABILITY  Created dependencies should occur exactly on partitions where
            the class satisfies strong lumpability: states in a block send
            equal aggregate mass into every block.

CRITERIA (sufficient directions only; necessity is not claimed)
  P1  EQUAL class: suffix laws of resets 0 and 3 coincide for every generator,
      every observation partition. Expected to pass by construction.
  P2  PERM class: whenever pi maps every reachable destination to a
      behaviourally equivalent state, the suffix laws coincide. Zero
      counterexamples required. The converse rate is reported, not required.
  P3  LUMPABLE class: the partitions admitting created dependencies are
      exactly those satisfying strong lumpability.
"""
import numpy as np, math, itertools
from collections import defaultdict

NS, T = 4, 5
ROWS = [r for r in itertools.product(range(3), repeat=NS) if sum(r) == 2]
NR = len(ROWS); ROWA = np.array(ROWS, dtype=np.int64)
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
PARTS = [sorted([sorted(b) for b in p]) for p in partitions([0, 1, 2, 3])]
PARTS = [p for p in PARTS if len(p) > 1]
PARTS.sort(key=lambda p: (len(p), [-len(b) for b in p]))

def law_from(M, s, obsv, k, steps):
    """Exact law of (Y_0..Y_{steps-1}) starting deterministically in state s."""
    res = []
    def rec(alpha, d):
        if d == steps:
            res.append(int(alpha.sum())); return
        nxt = alpha @ M
        for y in range(k):
            rec(np.where(obsv == y, nxt, 0), d + 1)
    a0 = np.zeros(NS, dtype=np.int64); a0[s] = 2 ** (steps)
    for y0 in range(k):
        rec(np.where(obsv == y0, a0, 0), 1)
    return tuple(res)

def suffix_law(M, s, obsv, k, steps):
    """Law of (Y_1..Y_{steps-1}): one transition, then observe."""
    nxt = np.zeros(NS, dtype=np.int64); nxt[s] = 1
    nxt = nxt @ M
    res = []
    def rec(alpha, d):
        if d == steps - 1:
            res.append(int(alpha.sum())); return
        nn = alpha @ M
        for y in range(k):
            rec(np.where(obsv == y, nn, 0), d + 1)
    for y in range(k):
        rec(np.where(obsv == y, nxt, 0) * (2 ** steps), 1)
    return tuple(res)

def main():
    print("=" * 78)
    print("EXPERIMENT 30  when does a generator dependency survive observation?")
    print("=" * 78)

    trip = list(itertools.product(range(NR), repeat=3))
    equal = [(a, b, c, a) for a, b, c in trip]
    perm  = [(a, b, c, PERM_OF[a]) for a, b, c in trip]
    lump  = [(r0, r1, r2, r3) for a in (0,1,2) for b in (0,1,2)
             for r0 in BYPROF[a] for r1 in BYPROF[a]
             for r2 in BYPROF[b] for r3 in BYPROF[b]]

    # ---------- P1 ----------
    bad1 = 0
    for part in PARTS:
        k = len(part)
        obsv = np.zeros(NS, dtype=np.int64)
        for bi, b in enumerate(part):
            for x in b: obsv[x] = bi
        for g in equal[::7]:
            M = ROWA[list(g)]
            if suffix_law(M, 0, obsv, k, T) != suffix_law(M, 3, obsv, k, T):
                bad1 += 1
    print(f"\n  P1  IDENTITY (R3 = R0): suffix laws of resets 0 and 3 coincide")
    print(f"      counterexamples over {len(equal[::7])*len(PARTS)} generator/partition "
          f"pairs: {bad1}   {'PASS' if bad1 == 0 else 'FAIL'}")

    # ---------- P2 ----------
    compat_eq = compat_ne = incompat_eq = incompat_ne = 0
    for part in PARTS:
        k = len(part)
        obsv = np.zeros(NS, dtype=np.int64)
        for bi, b in enumerate(part):
            for x in b: obsv[x] = bi
        for g in perm[::7]:
            M = ROWA[list(g)]
            beh = {x: law_from(M, x, obsv, k, T - 1) for x in range(NS)}
            reach = [y for y in range(NS) if M[0, y] > 0]
            compat = all(beh[y] == beh[PERM[y]] for y in reach)
            same = suffix_law(M, 0, obsv, k, T) == suffix_law(M, 3, obsv, k, T)
            if compat and same: compat_eq += 1
            elif compat and not same: compat_ne += 1
            elif same: incompat_eq += 1
            else: incompat_ne += 1
    print(f"\n  P2  PERMUTATION (R3 = pi(R0)): does behavioural compatibility of pi")
    print(f"      imply the suffix laws coincide?")
    print(f"      {'':<22}{'suffix equal':>14}{'suffix differs':>16}")
    print(f"      {'pi compatible':<22}{compat_eq:>14}{compat_ne:>16}")
    print(f"      {'pi incompatible':<22}{incompat_eq:>14}{incompat_ne:>16}")
    print(f"      counterexamples (compatible but differing): {compat_ne}   "
          f"{'PASS' if compat_ne == 0 else 'FAIL'}")
    if incompat_eq:
        print(f"      sufficient but not necessary: {incompat_eq} cases coincide "
              f"without pi being compatible")

    # ---------- P3 ----------
    print(f"\n  P3  LUMPABILITY: which partitions admit created dependencies?")
    print(f"      {'partition':<12}{'strongly lumpable':>20}{'creates deps':>14}")
    ok3 = True
    for part in PARTS:
        k = len(part)
        obsv = np.zeros(NS, dtype=np.int64)
        for bi, b in enumerate(part):
            for x in b: obsv[x] = bi
        # strong lumpability of the CLASS: every generator lumpable under this partition
        lumpable = True
        for g in lump[::23]:
            M = ROWA[list(g)]
            for blk in part:
                sig = {tuple(int(M[s, list(t)].sum()) for t in part) for s in blk}
                if len(sig) > 1:
                    lumpable = False; break
            if not lumpable: break
        # does reset 0 equal reset 1 as a random variable over the class?
        z0 = [law_from(ROWA[list(g)], 0, obsv, k, T) for g in lump[::23]]
        z1 = [law_from(ROWA[list(g)], 1, obsv, k, T) for g in lump[::23]]
        creates = (z0 == z1)
        lab = "|".join("".join(map(str, b)) for b in part)
        print(f"      {lab:<12}{('yes' if lumpable else 'no'):>20}"
              f"{('yes' if creates else 'no'):>14}")
        if lumpable != creates:
            ok3 = False
    print(f"      lumpability matches dependency creation exactly: "
          f"{'PASS' if ok3 else 'FAIL'}")

    # ---- POST-HOC refinement, not pre-registered ----
    print("\n  POST-HOC. P3 fails on the discrete partition 0|1|2|3, which is")
    print("  strongly lumpable only vacuously: with every block a singleton no two")
    print("  states share a block, so there is nothing for the quotient to collapse.")
    print("  Refined criterion: lumpable AND some block of size > 1.")
    ok3b = True
    for part in PARTS:
        k = len(part)
        obsv = np.zeros(NS, dtype=np.int64)
        for bi, b in enumerate(part):
            for x in b: obsv[x] = bi
        lumpable = True
        for g in lump[::23]:
            M = ROWA[list(g)]
            for blk in part:
                if len({tuple(int(M[s, list(t)].sum()) for t in part) for s in blk}) > 1:
                    lumpable = False; break
            if not lumpable: break
        nontrivial = any(len(b) > 1 for b in part)
        z0 = [law_from(ROWA[list(g)], 0, obsv, k, T) for g in lump[::23]]
        z1 = [law_from(ROWA[list(g)], 1, obsv, k, T) for g in lump[::23]]
        if (lumpable and nontrivial) != (z0 == z1):
            ok3b = False
    print(f"  refined criterion matches dependency creation on all {len(PARTS)} "
          f"partitions: {'yes' if ok3b else 'no'}")

main()
