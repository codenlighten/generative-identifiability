#!/usr/bin/env python3
"""
EXPERIMENT 28 -- which generator dependencies survive the observation map?

Frozen before running.

Experiment 27 showed that tying two transition rows collapses the minimal
intervention basis from four resets to three, while a size-matched random
thinning does not. That is one dependency of one kind. It establishes a
sufficient mechanism, not an equivalence, and it leaves the real question
open: a constraint among generator coordinates need not survive into a
functional dependency among reset outputs, because Z_j depends on the whole
generator and not only on row j.

The composition under test is

    generator constraints  ->  reset laws  ->  functional dependencies among Z_j

and the observation map sits in the middle, able in principle to preserve,
destroy, or manufacture a dependency.

------------------------------------------------------------------- CLASSES
All on four hidden states with rows drawn from the same 10-row set. Every
constrained class has exactly 1,000 generators, so cardinality is held fixed
and only the constraint varies.

  FREE        all 10^4, rows independent                        (reference)
  EQUAL       R3 = R0
  PERM        R3 = pi(R0), pi swapping coordinates 0<->1, 2<->3
  DERIVED     R3 = row index (idx(R0) + idx(R1)) mod 10
  ADJACENT    R2 = R1   (a different pair tied, same shape)
  RANDOM      1,000 generators drawn uniformly from FREE, fixed seed (control)

----------------------------------------------------------------- CLOSURES
Two closure operators, both tested exactly by integer class counts:

  cl_G(A) = { j : H(R_j | R_A) = 0 }    among generator coordinates
  cl_Z(A) = { j : H(Z_j | Z_A) = 0 }    among reset outputs

------------------------------------------------------- CRITERIA (Z1-Z3)
Z1  DESCRIPTIVE. For each class, report the minimal bases of cl_Z and whether
    the class compresses. No direction is pre-committed.

Z2  SEPARATION CANDIDATE. Within each class, the set of minimal bases is the
    same for all 14 observation partitions -- i.e. the observation map moves
    rank values but not basis structure. Reported per class; a single class
    failing refutes the general form.

Z3  INCLUSION. For every class, every A and every j,
    j in cl_G(A)  implies  j in cl_Z(A).  Zero violations required.
    The converse is reported descriptively: a j in cl_Z(A) but not cl_G(A) is
    a dependency the observation map created rather than inherited.
"""
import numpy as np, math, itertools
from collections import defaultdict

NS, T = 4, 5
ROWS = [r for r in itertools.product(range(3), repeat=NS) if sum(r) == 2]
NR = len(ROWS)
ROWA = np.array(ROWS, dtype=np.int64)
PERM = (1, 0, 3, 2)
PIDX = {r: i for i, r in enumerate(ROWS)}
PERM_OF = [PIDX[tuple(ROWS[i][p] for p in PERM)] for i in range(NR)]

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

def classes_def():
    free = np.array(list(itertools.product(range(NR), repeat=NS)), dtype=np.int64)
    trip = list(itertools.product(range(NR), repeat=3))
    equal    = np.array([[a, b, c, a]              for a, b, c in trip], dtype=np.int64)
    perm     = np.array([[a, b, c, PERM_OF[a]]     for a, b, c in trip], dtype=np.int64)
    derived  = np.array([[a, b, c, (a + b) % NR]   for a, b, c in trip], dtype=np.int64)
    adjacent = np.array([[a, b, b, c]              for a, b, c in trip], dtype=np.int64)
    rng = np.random.default_rng(28)
    rand = free[np.sort(rng.choice(len(free), size=len(equal), replace=False))]
    return [("FREE", free), ("EQUAL", equal), ("PERM", perm),
            ("DERIVED", derived), ("ADJACENT", adjacent), ("RANDOM", rand)]

def ncls(cols):
    if not cols:
        return 1
    comb = np.stack(cols, axis=1).astype(np.int64)
    v = np.ascontiguousarray(comb).view(
        np.dtype((np.void, comb.dtype.itemsize * comb.shape[1]))).ravel()
    return len(np.unique(v))

def laws_ids(rowidx, obsv, k, j):
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
        a0 = np.broadcast_to(np.where(obsv == y0, pm, 0).astype(np.int64),
                             (n, NS)).copy()
        rec(a0, 1)
    v = np.ascontiguousarray(out).view(
        np.dtype((np.void, out.dtype.itemsize * out.shape[1]))).ravel()
    _, inv = np.unique(v, return_inverse=True)
    return inv.astype(np.int64)

SUBSETS = [frozenset(A) for r in range(NS + 1)
           for A in itertools.combinations(range(NS), r)]

def closure(colfn):
    cl = {}
    for A in SUBSETS:
        base = ncls([colfn(j) for j in sorted(A)])
        cl[A] = frozenset(j for j in range(NS)
                          if ncls([colfn(i) for i in sorted(A | {j})]) == base)
    return cl

def bases(cl):
    J = frozenset(range(NS))
    span = [A for A in SUBSETS if cl[A] == J]
    return sorted((sorted(A) for A in span if not any(set(B) < set(A) for B in span)),
                  key=lambda a: (len(a), a))

def main():
    print("=" * 78)
    print("EXPERIMENT 28  which generator dependencies survive the observation map?")
    print("=" * 78)
    incl_viol = 0
    created = 0
    z2 = {}
    for name, rowidx in classes_def():
        n = len(rowidx)
        clG = closure(lambda j: rowidx[:, j])
        gdep = sorted(f"{sorted(A)}->{sorted(clG[A] - A)}"
                      for A in SUBSETS if clG[A] - A)
        print(f"\n  {name:<9} |H| = {n:>6,}   generator-coordinate dependencies: "
              f"{len(gdep)}")
        if gdep:
            print(f"    e.g. {gdep[0]}" + (f", {gdep[1]}" if len(gdep) > 1 else ""))
        per_part = []
        for part in PARTS:
            k = len(part)
            obsv = np.empty(NS, dtype=np.int64)
            for bi, b in enumerate(part):
                for x in b:
                    obsv[x] = bi
            zid = {j: laws_ids(rowidx, obsv, k, j) for j in range(NS)}
            clZ = closure(lambda j: zid[j])
            for A in SUBSETS:
                if not (clG[A] <= clZ[A]):
                    incl_viol += 1
                created += len(clZ[A] - clG[A])
            per_part.append((("|".join("".join(map(str, b)) for b in part)),
                             tuple(tuple(b) for b in bases(clZ))))
        distinct = {b for _, b in per_part}
        z2[name] = len(distinct) == 1
        shown = ", ".join("{" + ",".join(map(str, b)) + "}" for b in sorted(distinct)[0])
        print(f"    minimal intervention bases: {shown}"
              f"   (size {len(sorted(distinct)[0][0])})")
        print(f"    invariant across all {len(PARTS)} observation partitions: "
              f"{'yes' if z2[name] else 'NO -- ' + str(len(distinct)) + ' distinct sets'}")

    print()
    print("=" * 78)
    print("PRE-REGISTERED CRITERIA")
    print("=" * 78)
    allz2 = all(z2.values())
    print(f"  Z2  basis structure invariant across observation partitions, "
          f"every class:  {'PASS' if allz2 else 'FAIL'}")
    for nm, v in z2.items():
        print(f"        {nm:<9} {'invariant' if v else 'VARIES'}")
    print(f"  Z3  cl_G(A) subset cl_Z(A) always:  "
          f"{'PASS' if incl_viol == 0 else 'FAIL'}   ({incl_viol} violations)")
    print(f"      dependencies present in cl_Z but not cl_G (observation-created): "
          f"{created}")

main()
