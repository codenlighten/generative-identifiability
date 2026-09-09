#!/usr/bin/env python3
"""
EXPERIMENT 20 -- pre-registered replication on a new stochastic family.

This file was written and frozen BEFORE it was run for the first time.  No
threshold, metric, or criterion below was chosen or adjusted with knowledge of
the outcome.  Experiment 19's interventional quadrant result was post-hoc; this
is the confirmatory test of it on a system that is not a parameter tweak of
Experiment 19.

--------------------------------------------------------------------- SYSTEM
Hidden state      X_t in {0,1,2,3}                     (4 states, was 3)
Observation       O(0)=O(1)=0,  O(2)=O(3)=1            (balanced 2+2 lumping)
Transition rows   (p0,p1,p2,p3), p_i in {0,1/4,1/2,3/4,1}, sum 1  -> 35 rows
Hypothesis class  ALL 35^4 = 1,500,625 generators, enumerated completely
Horizon           T = 6 observed symbols

------------------------------------------------------------------ PROTOCOLS
PASSIVE         initial hidden state uniform, P(X_0=i)=1/4
INTERVENTIONAL  the experimenter may initialise each hidden state separately;
                the observation is the 4 resulting laws pooled.  Intervention
                means this and only this in both protocols.

-------------------------------------------------------------------- METRICS
Generator uncertainty   I(G)   = log2 |{G' : law(Y_0..Y_5 | G') = law(.. | G)}|
Behavioural uncertainty P_1(G) = H(Y_5 | Y_0..Y_4, G)   -- generator KNOWN,
                                 evaluated under the uniform operating
                                 condition in BOTH protocols, so that only the
                                 identification channel differs between them.

Quadrant labels:  identifiable  <=>  I(G) = 0
                  predictable   <=>  P_1(G) < 0.10 bits/step
(Both carried over unchanged from Experiment 19.)

------------------------------------------------------- CRITERIA (C1-C5)
C1  Quadrant support.  Under the interventional protocol each of the four
    quadrants holds at least 1% of the hypothesis class.

C2  Identifiability does not determine predictability.  Among generators with
    I(G)=0, the 10th-90th percentile spread of P_1 is at least 0.50 bits.

C3  Predictability does not determine identifiability.  Among generators whose
    P_1 lies within 0.02 bits of the class median P_1, the 10th-90th percentile
    spread of I is at least 1.00 bit.

C4  Exact observational-equivalence witness.  Some G1 != G2 induce an identical
    observable law while sharing P_1 < 0.10 bits/step.

C5  Intervention changes regime.  At least 1% of generators occupy a different
    quadrant under the interventional protocol than under the passive one.

C2 and C3 are distributional rather than existence tests on purpose: over 1.5M
generators an existence claim is nearly free, whereas a percentile spread can
fail.
"""
import numpy as np, math, itertools

ROWS = np.array([r for r in itertools.product(range(5), repeat=4) if sum(r) == 4],
                dtype=np.int64)                     # 35 rows, denominator 4
NS, T = 4, 6
OBSV = np.array([0, 0, 1, 1])
DEN = 4**T
NGEN = len(ROWS)**NS
CHUNK = 100_000

def laws_for(chunk_idx, pi_vec):
    """Exact law of (Y_0..Y_{T-1}) for a chunk of generators, int16 over DEN."""
    n = len(chunk_idx)
    i0 = chunk_idx // (35**3) % 35; i1 = chunk_idx // (35**2) % 35
    i2 = chunk_idx // 35 % 35;      i3 = chunk_idx % 35
    M = np.empty((n, 4, 4), dtype=np.int64)
    M[:, 0, :] = ROWS[i0]; M[:, 1, :] = ROWS[i1]
    M[:, 2, :] = ROWS[i2]; M[:, 3, :] = ROWS[i3]
    out = np.empty((n, 1 << T), dtype=np.int16)
    col = [0]
    def rec(alpha, depth):
        if depth == T:
            out[:, col[0]] = alpha.sum(axis=1); col[0] += 1; return
        nxt = np.matmul(alpha[:, None, :], M)[:, 0, :]
        for y in (0, 1):
            rec(np.where(OBSV[None, :] == y, nxt, 0), depth + 1)
    for y0 in (0, 1):
        a0 = np.where(OBSV[None, :] == y0, pi_vec[None, :], 0).astype(np.int64)
        a0 = np.broadcast_to(a0, (n, 4)).copy()
        rec(a0, 1)
    return out

def shannon(counts):
    p = counts.astype(np.float64) / DEN
    with np.errstate(divide='ignore', invalid='ignore'):
        t = np.where(p > 0, -p * np.log2(p), 0.0)
    return t.sum(axis=1)

def ids_of(arr):
    """Exact equivalence-class ids for rows of a 2-D integer array."""
    v = np.ascontiguousarray(arr).view(
        np.dtype((np.void, arr.dtype.itemsize * arr.shape[1]))).ravel()
    _, inv, cnt = np.unique(v, return_inverse=True, return_counts=True)
    return inv, cnt

def main():
    print("=" * 78)
    print("EXPERIMENT 20  pre-registered replication, 4-state partially observed chain")
    print("=" * 78)
    print(f"  |H| = {NGEN:,} generators ({math.log2(NGEN):.2f} bits), "
          f"T = {T}, balanced 2+2 lumping\n")

    passive = np.empty((NGEN, 1 << T), dtype=np.int16)
    inter = [np.empty((NGEN, 1 << T), dtype=np.int16) for _ in range(NS)]
    P1 = np.empty(NGEN, dtype=np.float64)
    uni = np.array([1, 1, 1, 1], dtype=np.int64)

    for s in range(0, NGEN, CHUNK):
        idx = np.arange(s, min(s + CHUNK, NGEN))
        L = laws_for(idx, uni)
        passive[idx] = L
        P1[idx] = shannon(L) - shannon(L[:, ::2] + L[:, 1::2])
        for k in range(NS):
            pm = np.zeros(NS, dtype=np.int64); pm[k] = 4
            inter[k][idx] = laws_for(idx, pm)
        if (s // CHUNK) % 5 == 0:
            print(f"    ... {min(s+CHUNK, NGEN):>9,} / {NGEN:,}")

    inv_p, cnt_p = ids_of(passive)
    I_pass = np.log2(cnt_p[inv_p])
    comb = np.stack([ids_of(inter[k])[0] for k in range(NS)], axis=1).astype(np.int64)
    inv_i, cnt_i = ids_of(comb)
    I_int = np.log2(cnt_i[inv_i])

    HGO_p = (cnt_p * np.log2(cnt_p)).sum() / NGEN
    HGO_i = (cnt_i * np.log2(cnt_i)).sum() / NGEN
    print(f"\n  H(G|O) passive       : {HGO_p:.2f} bits over {len(cnt_p):,} distinct laws")
    print(f"  H(G|O) intervening   : {HGO_i:.2f} bits over {len(cnt_i):,} distinct laws")
    print(f"  mean P_1             : {P1.mean():.2f} bits/step")

    PT = 0.10
    def quad(I):
        return (I == 0.0).astype(np.int8) * 2 + (P1 < PT).astype(np.int8)
    for name, I in (("PASSIVE", I_pass), ("INTERVENING", I_int)):
        q = quad(I)
        print(f"\n  {name}")
        print(f"  {'':<20}{'predictable':>20}{'not predictable':>20}")
        for bit, lab in ((2, "identifiable"), (0, "not identifiable")):
            a = int((q == bit + 1).sum()); b = int((q == bit).sum())
            print(f"  {lab:<20}{a:>11,} ({100*a/NGEN:5.2f}%)"
                  f"{b:>11,} ({100*b/NGEN:5.2f}%)")

    # ---------------- criteria ----------------
    qi = quad(I_int); qp = quad(I_pass)
    c1 = all((qi == v).sum() >= 0.01 * NGEN for v in (0, 1, 2, 3))

    sel = I_int == 0.0
    spread_P = (np.percentile(P1[sel], 90) - np.percentile(P1[sel], 10)) if sel.any() else 0.0
    c2 = spread_P >= 0.50

    med = np.median(P1)
    sel2 = np.abs(P1 - med) < 0.02
    spread_I = (np.percentile(I_int[sel2], 90) - np.percentile(I_int[sel2], 10)) if sel2.any() else 0.0
    c3 = spread_I >= 1.00

    dup = (cnt_p[inv_p] > 1) & (P1 < PT)
    c4 = bool(dup.any())

    migr = float((qi != qp).mean())
    c5 = migr >= 0.01

    print("\n" + "=" * 78)
    print("PRE-REGISTERED CRITERIA")
    print("=" * 78)
    print(f"  C1  four quadrants each >= 1% (intervening)          "
          f"{'PASS' if c1 else 'FAIL'}")
    print(f"  C2  P_1 spread among identifiable >= 0.50 bits       "
          f"{'PASS' if c2 else 'FAIL'}   ({spread_P:.2f}, n={int(sel.sum()):,})")
    print(f"  C3  I spread at fixed P_1 >= 1.00 bit                "
          f"{'PASS' if c3 else 'FAIL'}   ({spread_I:.2f}, n={int(sel2.sum()):,})")
    print(f"  C4  exact observational-equivalence witness          "
          f"{'PASS' if c4 else 'FAIL'}   ({int(dup.sum()):,} generators)")
    print(f"  C5  >= 1% change quadrant under intervention         "
          f"{'PASS' if c5 else 'FAIL'}   ({migr:.1%})")
    print(f"\n  overall: {sum([c1,c2,c3,c4,c5])}/5 criteria met")

main()
