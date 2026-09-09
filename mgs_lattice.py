#!/usr/bin/env python3
"""
EXPERIMENTS 22 and 23 -- theorem-seeking on Experiment 20's family.

Not a new system. The same 35^4 = 1,500,625 four-state partially observed
Markov chains, interrogated for structure rather than for another data point.

E22  INTERVENTION LATTICE.  The interventional observation is the tuple of
     the four laws obtained by resetting to each hidden state; the passive
     observation is their uniform mixture. If so, passive is a deterministic
     coarsening of interventional, and preimage containment is forced rather
     than observed. We check the mixture identity exactly, then verify
     containment over every pair A subset B of reset sets, and test whether
     information gain is submodular.

E23  HORIZON REFINEMENT.  Equivalence at horizon T can split but never merge
     as T grows, so i_T(g) is non-increasing. We locate T* -- the horizon at
     which the partition stops refining -- which is what would license
     upgrading finite-horizon statements to full observational equivalence.

CLAIMS TESTED (stated before running):
  P1  passive law == uniform mixture of the four reset laws, exactly
  P2  A subset B  =>  i_B(g) <= i_A(g) for every generator (no exceptions)
  P3  quadrant motion under refinement is one-way: not-identifiable ->
      identifiable, never the reverse, because p_1 is held fixed
  P4  information gain is submodular: for A subset B, j not in B,
      H(G|O_A) - H(G|O_A+j) >= H(G|O_B) - H(G|O_B+j)
  P5  the horizon partition stabilises at some finite T* <= 8
"""
import numpy as np, math, itertools

ROWS = np.array([r for r in itertools.product(range(5), repeat=4) if sum(r) == 4],
                dtype=np.int64)
NS = 4
OBSV = np.array([0, 0, 1, 1])
NGEN = len(ROWS)**NS
CHUNK = 60_000

def laws_for(idx, pi_vec, T):
    n = len(idx)
    i0 = idx // (35**3) % 35; i1 = idx // (35**2) % 35
    i2 = idx // 35 % 35;      i3 = idx % 35
    M = np.empty((n, 4, 4), dtype=np.int64)
    M[:, 0, :] = ROWS[i0]; M[:, 1, :] = ROWS[i1]
    M[:, 2, :] = ROWS[i2]; M[:, 3, :] = ROWS[i3]
    out = np.empty((n, 1 << T), dtype=np.int32)
    col = [0]
    def rec(alpha, depth):
        if depth == T:
            out[:, col[0]] = alpha.sum(axis=1); col[0] += 1; return
        nxt = np.matmul(alpha[:, None, :], M)[:, 0, :]
        for y in (0, 1):
            rec(np.where(OBSV[None, :] == y, nxt, 0), depth + 1)
    for y0 in (0, 1):
        a0 = np.broadcast_to(np.where(OBSV == y0, pi_vec, 0).astype(np.int64),
                             (n, 4)).copy()
        rec(a0, 1)
    return out

def ids_of(arr):
    v = np.ascontiguousarray(arr).view(
        np.dtype((np.void, arr.dtype.itemsize * arr.shape[1]))).ravel()
    _, inv, cnt = np.unique(v, return_inverse=True, return_counts=True)
    return inv.astype(np.int64), cnt

def hgo(cnt):
    return float((cnt * np.log2(cnt)).sum() / NGEN)

def main():
    T_I, T_H = 6, 8
    uni = np.array([1, 1, 1, 1], dtype=np.int64)

    # ---------- per-start laws at T=6, plus the passive law ----------
    print("=" * 78); print("EXPERIMENT 22  intervention lattice"); print("=" * 78)
    reset_ids, reset_cnt = [], []
    mixture_ok = True
    passive6 = np.empty((NGEN, 1 << T_I), dtype=np.int32)
    resets6 = [np.empty((NGEN, 1 << T_I), dtype=np.int32) for _ in range(NS)]
    for s in range(0, NGEN, CHUNK):
        idx = np.arange(s, min(s + CHUNK, NGEN))
        passive6[idx] = laws_for(idx, uni, T_I)
        acc = np.zeros((len(idx), 1 << T_I), dtype=np.int64)
        for k in range(NS):
            pm = np.zeros(NS, dtype=np.int64); pm[k] = 4
            L = laws_for(idx, pm, T_I)
            resets6[k][idx] = L
            acc += L
        if not np.array_equal(acc // 4, passive6[idx].astype(np.int64)):
            mixture_ok = False
    print(f"  P1  passive law == (1/4) * sum of the four reset laws, exactly: "
          f"{'CONFIRMED' if mixture_ok else 'REFUTED'}")
    print("      so the passive experiment is a deterministic coarsening of the")
    print("      interventional one, and preimage containment is forced, not observed.")

    for k in range(NS):
        i, c = ids_of(resets6[k]); reset_ids.append(i); reset_cnt.append(c)
    pid, pcnt = ids_of(passive6)
    i_passive = np.log2(pcnt[pid])

    # p_1 under the uniform operating condition, exactly as in Experiment 20
    def shannon(a):
        p = a.astype(np.float64) / (4 ** T_I)
        with np.errstate(divide='ignore', invalid='ignore'):
            return np.where(p > 0, -p * np.log2(p), 0.0).sum(axis=1)
    p1 = shannon(passive6) - shannon(passive6[:, ::2] + passive6[:, 1::2])
    del passive6

    # ---------- every subset of resets ----------
    sub_i, sub_H = {}, {}
    print(f"\n  {'reset set A':<16}{'H(G|O_A)':>10}{'classes':>11}"
          f"{'P[i=0]':>9}{'max |P|':>10}")
    order = sorted(itertools.chain.from_iterable(
        itertools.combinations(range(NS), r) for r in range(NS + 1)), key=lambda a: (len(a), a))
    for A in order:
        if not A:
            ids = np.zeros(NGEN, dtype=np.int64); cnt = np.array([NGEN])
        else:
            comb = np.stack([reset_ids[k] for k in A], axis=1)
            ids, cnt = ids_of(comb)
        iA = np.log2(cnt[ids]); sub_i[A] = iA; sub_H[A] = hgo(cnt)
        lab = "{}" if not A else "{" + ",".join(map(str, A)) + "}"
        print(f"  {lab:<16}{sub_H[A]:>10.3f}{len(cnt):>11,}"
              f"{float((iA == 0).mean()):>9.2%}{int(cnt.max()):>10,}")
    print(f"  {'passive (mixture)':<16}{hgo(pcnt):>10.3f}{len(pcnt):>11,}"
          f"{float((i_passive == 0).mean()):>9.2%}{int(pcnt.max()):>10,}")

    # ---------- P2: containment over the lattice ----------
    viol = 0
    for A in order:
        for B in order:
            if set(A) < set(B) and np.any(sub_i[B] > sub_i[A] + 1e-12):
                viol += 1
    print(f"\n  P2  A subset B => i_B(g) <= i_A(g) for every generator: "
          f"{'CONFIRMED' if viol == 0 else f'REFUTED ({viol} pairs)'}")
    v2 = int(np.sum(sub_i[(0, 1, 2, 3)] > i_passive + 1e-12))
    print(f"      full intervention vs passive mixture: {v2} violations")

    # ---------- value of intervention ----------
    V = hgo(pcnt) - sub_H[(0, 1, 2, 3)]
    print(f"\n  value of full intervention over passive: "
          f"{hgo(pcnt):.3f} - {sub_H[(0,1,2,3)]:.3f} = {V:.3f} bits")
    print(f"  = I(G ; O_I | O_P), since O_P is a function of O_I")

    # ---------- P3: one-way quadrant motion ----------
    PT = 0.10
    qp = (i_passive == 0).astype(np.int8) * 2 + (p1 < PT).astype(np.int8)
    qi = (sub_i[(0, 1, 2, 3)] == 0).astype(np.int8) * 2 + (p1 < PT).astype(np.int8)
    lost = int(np.sum((qp >= 2) & (qi < 2)))
    gained = int(np.sum((qp < 2) & (qi >= 2)))
    pred_moved = int(np.sum((qp % 2) != (qi % 2)))
    print(f"\n  P3  quadrant motion under refinement is one-way: "
          f"{'CONFIRMED' if lost == 0 and pred_moved == 0 else 'REFUTED'}")
    print(f"      became identifiable {gained:,}; lost identifiability {lost}; "
          f"predictive coordinate moved {pred_moved}")

    # ---------- P4: submodularity ----------
    bad, tested, worst = 0, 0, 0.0
    for A in order:
        for B in order:
            if not set(A) <= set(B): continue
            for j in range(NS):
                if j in B: continue
                dA = sub_H[A] - sub_H[tuple(sorted(set(A) | {j}))]
                dB = sub_H[B] - sub_H[tuple(sorted(set(B) | {j}))]
                tested += 1
                if dA < dB - 1e-9:
                    bad += 1; worst = max(worst, dB - dA)
    print(f"\n  P4  information gain submodular (diminishing returns): "
          f"{'CONFIRMED' if bad == 0 else 'REFUTED'}")
    print(f"      {tested} (A,B,j) triples tested, {bad} violations"
          + (f", worst {worst:.4f} bits" if bad else ""))

    # ---------- greedy gain curve ----------
    print(f"\n  greedy intervention sequence (each step adds the best next reset):")
    cur, chosen = (), []
    print(f"    passive-equivalent start   H = {sub_H[()]:.3f} bits (no observation)")
    while len(cur) < NS:
        best = max((k for k in range(NS) if k not in cur),
                   key=lambda k: sub_H[cur] - sub_H[tuple(sorted(set(cur) | {k}))])
        nxt = tuple(sorted(set(cur) | {best}))
        print(f"    + reset state {best}            H = {sub_H[nxt]:.3f} bits "
              f"(gain {sub_H[cur]-sub_H[nxt]:.3f})")
        cur = nxt
    del resets6, reset_ids

    # ---------- E23: horizon refinement ----------
    print()
    print("=" * 78); print("EXPERIMENT 23  horizon refinement (passive protocol)"); print("=" * 78)
    arr = np.empty((NGEN, 1 << T_H), dtype=np.int32)
    for s in range(0, NGEN, CHUNK):
        idx = np.arange(s, min(s + CHUNK, NGEN))
        arr[idx] = laws_for(idx, uni, T_H)
    rows = []
    for T in range(T_H, 0, -1):
        ids, cnt = ids_of(arr)
        rows.append((T, hgo(cnt), len(cnt), int(cnt.max()),
                     float((cnt[ids] == 1).mean())))
        if T > 1:
            arr = arr[:, ::2] + arr[:, 1::2]
    rows.reverse()
    print(f"  {'T':>3}{'H(G|O_T)':>11}{'classes':>12}{'max |P|':>12}{'identified':>12}")
    for T, H, nc, mx, fr in rows:
        print(f"  {T:>3}{H:>11.3f}{nc:>12,}{mx:>12,}{fr:>12.2%}")
    stab = None
    for a, b in zip(rows, rows[1:]):
        if a[2] == b[2] and stab is None:
            stab = a[0]
    print(f"\n  P5  partition stabilises at a finite T* <= {T_H}: "
          + (f"CONFIRMED, T* = {stab}" if stab else "NOT OBSERVED within this range"))
    print("      (class count is non-decreasing in T by construction; stabilisation")
    print("       is what would license upgrading finite-T claims to full")
    print("       observational equivalence.)")

main()
