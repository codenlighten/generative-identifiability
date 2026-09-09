#!/usr/bin/env python3
"""
EXPERIMENT 27 -- is "no experimental compression" a property of interventions,
or of the hypothesis class being a product?

Frozen before running.

Experiment 26 found that across all 14 non-trivial observation partitions the
only spanning set of the intervention polymatroid is J itself: no reset is
redundant given the others. That family is a full product, H = ROWS^4, with
the four transition rows chosen independently. A product class is exactly the
kind that should resist compression, so the result may say nothing about
interventions at all.

We compare three hypothesis classes on the SAME system and the SAME 14
partitions:

  A  PRODUCT     all 10^4 = 10,000 generators, rows independent   (Exp 26)
  B  TIED        row 3 constrained equal to row 0, giving 10^3 = 1,000
  C  RANDOM      a uniformly random subset of A of size 1,000, fixed seed

B and C are matched in cardinality, so any difference between them is
structural rather than a size effect. C is the control that decides whether
compression, if it appears in B, is about the tie or merely about there being
fewer generators to tell apart.

CRITERIA
  Y1  Class B admits, for at least one partition, a spanning set smaller
      than J.
  Y2  Class C admits no such spanning set for any partition. (If Y2 fails,
      compression is a cardinality effect and the tie explains nothing.)
  Y3  Where B compresses, the omitted reset lies in the closure of the rest:
      H(Z_j | Z_{J minus j}) = 0 exactly.

Spanning is tested exactly by integer class counts, never by a tolerance.
"""
import numpy as np, math, itertools

NS, T = 4, 5
ROWS = np.array([r for r in itertools.product(range(3), repeat=NS) if sum(r) == 2],
                dtype=np.int64)
NR = len(ROWS)

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

def build_classes():
    prod = np.array(list(itertools.product(range(NR), repeat=NS)), dtype=np.int64)
    tied = np.array([[a, b, c, a] for a in range(NR) for b in range(NR)
                     for c in range(NR)], dtype=np.int64)
    rng = np.random.default_rng(27)
    pick = rng.choice(len(prod), size=len(tied), replace=False)
    rand = prod[np.sort(pick)]
    return [("PRODUCT (rows free)", prod), ("TIED (row3 = row0)", tied),
            ("RANDOM subset", rand)]

def laws(rowidx, pi_vec, obsv, k):
    n = len(rowidx)
    M = np.empty((n, NS, NS), dtype=np.int64)
    for pos in range(NS):
        M[:, pos, :] = ROWS[rowidx[:, pos]]
    out = np.empty((n, k ** T), dtype=np.int32)
    col = [0]
    def rec(alpha, depth):
        if depth == T:
            out[:, col[0]] = alpha.sum(axis=1); col[0] += 1; return
        nxt = np.matmul(alpha[:, None, :], M)[:, 0, :]
        for y in range(k):
            rec(np.where(obsv[None, :] == y, nxt, 0), depth + 1)
    for y0 in range(k):
        a0 = np.broadcast_to(np.where(obsv == y0, pi_vec, 0).astype(np.int64),
                             (n, NS)).copy()
        rec(a0, 1)
    return out

def nclasses(idlists, n):
    if not idlists:
        return 1
    comb = np.stack(idlists, axis=1)
    v = np.ascontiguousarray(comb).view(
        np.dtype((np.void, comb.dtype.itemsize * comb.shape[1]))).ravel()
    return len(np.unique(v))

def ids_for(rowidx, obsv, k, j):
    pm = np.zeros(NS, dtype=np.int64); pm[j] = 4
    L = laws(rowidx, pm, obsv, k)
    v = np.ascontiguousarray(L).view(
        np.dtype((np.void, L.dtype.itemsize * L.shape[1]))).ravel()
    _, inv = np.unique(v, return_inverse=True)
    return inv.astype(np.int64)

def main():
    print("=" * 78)
    print("EXPERIMENT 27  does the no-compression result survive a non-product class?")
    print("=" * 78)
    results = {}
    for name, rowidx in build_classes():
        n = len(rowidx)
        print(f"\n  {name}   |H| = {n:,} ({math.log2(n):.2f} bits)")
        print(f"  {'partition':<12}{'classes(J)':>12}{'identified':>12}"
              f"{'min basis':>11}   bases")
        small = 0
        for part in PARTS:
            k = len(part)
            obsv = np.empty(NS, dtype=np.int64)
            for bi, b in enumerate(part):
                for x in b:
                    obsv[x] = bi
            lab = "|".join("".join(map(str, b)) for b in part)
            ids = {j: ids_for(rowidx, obsv, k, j) for j in range(NS)}
            nJ = nclasses([ids[j] for j in range(NS)], n)
            spanning = []
            for r in range(NS + 1):
                for A in itertools.combinations(range(NS), r):
                    if nclasses([ids[j] for j in A], n) == nJ:
                        spanning.append(set(A))
            minimal = [A for A in spanning if not any(B < A for B in spanning)]
            msz = min(len(A) for A in minimal)
            if msz < NS:
                small += 1
            shown = ", ".join("{" + ",".join(map(str, sorted(A))) + "}"
                              for A in sorted(minimal, key=lambda a: (len(a), sorted(a)))[:3])
            if len(minimal) > 3: shown += ", ..."
            print(f"  {lab:<12}{nJ:>12,}{nJ/n:>12.1%}{msz:>11}   {shown}")
        results[name] = small
        print(f"    -> {small}/{len(PARTS)} partitions compress "
              f"(a spanning set smaller than J)")

    print()
    print("=" * 78)
    print("PRE-REGISTERED CRITERIA")
    print("=" * 78)
    y1 = results["TIED (row3 = row0)"] > 0
    y2 = results["RANDOM subset"] == 0
    print(f"  Y1  TIED class compresses somewhere              "
          f"{'PASS' if y1 else 'FAIL'}   "
          f"({results['TIED (row3 = row0)']}/{len(PARTS)} partitions)")
    print(f"  Y2  RANDOM class of equal size does not          "
          f"{'PASS' if y2 else 'FAIL'}   "
          f"({results['RANDOM subset']}/{len(PARTS)} partitions)")
    print(f"  PRODUCT class, for reference: "
          f"{results['PRODUCT (rows free)']}/{len(PARTS)} partitions compress")
    if y1 and y2:
        print("\n  -> compression is structural, not a cardinality effect.")
    elif y1 and not y2:
        print("\n  -> both classes compress at this size; the tie explains nothing")
        print("     that shrinking the class does not, and Experiment 26's result")
        print("     is about |H| relative to the resolving power of a subset.")
    else:
        print("\n  -> no compression anywhere; the Experiment 26 finding is not")
        print("     specific to product structure.")

main()
