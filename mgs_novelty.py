#!/usr/bin/env python3
"""
EXPERIMENT 25 -- what makes conditional intervention novelty large?

DESCRIPTIVE, not a hypothesis test. Experiment 24 refuted the crossing
conjecture; this decomposes why, rather than proposing a replacement rule to
test. No pre-registered criteria, and none should be read into the output.

Because each reset's observable law Z_j is an exact deterministic function of
the generator G, information gain is exactly an entropy of those variables:

    F(A) = I(G; Z_A) = H(Z_A),          H(G|O_A) = log2|H| - H(Z_A)

so the value of adding reset j to a set A is its conditional novelty

    Delta(j | A) = H(Z_j | Z_A)

and for a pair, H(Z_i, Z_j) = H(Z_i) + H(Z_j) - I(Z_i; Z_j). Block crossing
was only ever a proxy for this. We compute the underlying quantities over all
15 partitions of the four hidden states, on Experiment 24's family.
"""
import numpy as np, math, itertools
from collections import defaultdict

NS, T = 4, 5
ROWS = np.array([r for r in itertools.product(range(3), repeat=NS) if sum(r) == 2],
                dtype=np.int64)
NGEN = len(ROWS) ** NS
LOGN = math.log2(NGEN)

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

def ent(cnt):
    p = cnt / NGEN
    return float(-(p * np.log2(p)).sum())

def joint_ent(*idlists):
    comb = np.stack(idlists, axis=1)
    return ent(ids_of(comb)[1])

def main():
    print("=" * 78)
    print("EXPERIMENT 25  conditional novelty H(Z_j | Z_A), decomposed")
    print("=" * 78)
    print(f"  |H| = {NGEN:,} generators ({LOGN:.2f} bits), {len(PARTS)} non-trivial "
          f"partitions, T = {T}\n")

    same_I, cross_I = [], []
    by_blocksize = defaultdict(list)
    ident_checked = ident_ok = 0
    print(f"  {'partition':<12}{'state':>6}{'|B(j)|':>8}{'H(Z_j)':>9}"
          f"{'   pair redundancies I(Z_i;Z_j)'}")
    for part in PARTS:
        k = len(part)
        obsv = np.empty(NS, dtype=np.int64)
        for bi, b in enumerate(part):
            for x in b:
                obsv[x] = bi
        blk = {x: len(b) for b in part for x in b}
        lab = "|".join("".join(map(str, b)) for b in part)
        ids, H1 = {}, {}
        for j in range(NS):
            pm = np.zeros(NS, dtype=np.int64); pm[j] = 4
            ids[j] = ids_of(laws(pm, obsv, k))[0]
            H1[j] = ent(np.bincount(ids[j])[np.bincount(ids[j]) > 0])
            by_blocksize[blk[j]].append(H1[j])
        for j in range(NS):
            row = []
            for i in range(NS):
                if i == j: continue
                Hij = joint_ent(ids[i], ids[j])
                I = H1[i] + H1[j] - Hij
                row.append(f"{i}:{I:.2f}")
                if i < j:
                    (same_I if obsv[i] == obsv[j] else cross_I).append(I)
                    # exact identity check: log2|H| - H(G|O_{ij}) == H(Z_i,Z_j)
                    _, cnt = ids_of(np.stack([ids[i], ids[j]], axis=1))
                    hgo = float((cnt * np.log2(cnt)).sum() / NGEN)
                    ident_checked += 1
                    if abs((LOGN - hgo) - Hij) < 1e-9:
                        ident_ok += 1
            print(f"  {lab if j == 0 else '':<12}{j:>6}{blk[j]:>8}{H1[j]:>9.3f}"
                  f"   {'  '.join(row)}")
        print()

    print("=" * 78)
    print(f"  identity  log2|H| - H(G|O_A) == H(Z_A):  {ident_ok}/{ident_checked} pairs")
    print("  (definitional, so this checks the implementation, not the mathematics)\n")

    print("  individual intervention information by observation-block size:")
    for sz in sorted(by_blocksize):
        v = by_blocksize[sz]
        print(f"    |B(j)| = {sz}:  mean H(Z_j) = {sum(v)/len(v):6.3f} bits "
              f"over {len(v)} states   [{min(v):.3f}, {max(v):.3f}]")

    print("\n  pairwise redundancy I(Z_i;Z_j):")
    for name, v in (("same observation block", same_I), ("different blocks", cross_I)):
        if v:
            print(f"    {name:<26} mean {sum(v)/len(v):6.3f} bits over {len(v)} pairs"
                  f"   [{min(v):.3f}, {max(v):.3f}]")
    # ---- exact decomposition of the Experiment 24 reversal ----
    print("\n" + "=" * 78)
    print("  WHY THE CROSSING EFFECT REVERSES: the two cases are driven by")
    print("  different terms of the same decomposition")
    print("=" * 78)
    print(f"  {'partition':<10}{'pair':<8}{'cross?':>8}{'H(Z_i)+H(Z_j)':>15}"
          f"{'I(Z_i;Z_j)':>12}{'H(Z_i,Z_j)':>12}{'H(G|O_A)':>11}")
    for part, pairs in ((["0"], None), ):
        pass
    for lab_want, pr in (("0|123", [(1, 2), (0, 1)]), ("01|23", [(0, 1), (0, 2)])):
        for part in PARTS:
            lab = "|".join("".join(map(str, b)) for b in part)
            if lab != lab_want:
                continue
            k = len(part)
            obsv = np.empty(NS, dtype=np.int64)
            for bi, b in enumerate(part):
                for x in b:
                    obsv[x] = bi
            idm, H1 = {}, {}
            for j in range(NS):
                pm = np.zeros(NS, dtype=np.int64); pm[j] = 4
                idm[j] = ids_of(laws(pm, obsv, k))[0]
                c = np.bincount(idm[j]); H1[j] = ent(c[c > 0])
            for (i, j) in pr:
                Hij = joint_ent(idm[i], idm[j])
                I = H1[i] + H1[j] - Hij
                _, cnt = ids_of(np.stack([idm[i], idm[j]], axis=1))
                hgo = float((cnt * np.log2(cnt)).sum() / NGEN)
                print(f"  {lab:<10}{'{'+str(i)+','+str(j)+'}':<8}"
                      f"{('yes' if obsv[i] != obsv[j] else 'no'):>8}"
                      f"{H1[i]+H1[j]:>15.3f}{I:>12.3f}{Hij:>12.3f}{hgo:>11.3f}")
            # attribute the gap
            (a, b), (c, d) = pr
            Ha = H1[a] + H1[b]; Hc = H1[c] + H1[d]
            Ia = H1[a] + H1[b] - joint_ent(idm[a], idm[b])
            Ic = H1[c] + H1[d] - joint_ent(idm[c], idm[d])
            gap = (Ha - Ia) - (Hc - Ic)
            print(f"    gap {gap:+.3f} bits = individual term {Ha-Hc:+.3f} "
                  f"+ redundancy term {Ic-Ia:+.3f}")
            print()

    print("\n  A pair is valuable when both terms are large and the redundancy small:")
    print("     H(Z_i, Z_j) = H(Z_i) + H(Z_j) - I(Z_i; Z_j).")
    print("  Whether that tracks block crossing is exactly what Experiment 24 found")
    print("  it does not, in general.")

main()
