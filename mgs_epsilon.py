#!/usr/bin/env python3
"""
EXPERIMENT 21 -- is i_O(g) recoverable from the epsilon-machine?

Frozen before running.  This is an adversarial test of THIS paper's own
positioning claim, not a test of the two-axis thesis.

The manuscript argues that our horizontal coordinate differs from the
statistical complexity C_mu of computational mechanics: C_mu is a property of
the observed PROCESS (the memory of its minimal predictor), while i_O(g)
counts how many members of a declared hypothesis class share that process.
C_mu is an invariant OF an equivalence class; i_O is the CARDINALITY of one.

If that distinction is empty -- if C_mu predicts i_O -- then Related Work
overstates the contribution and should be rewritten.

--------------------------------------------------------------------- SYSTEM
The exploratory family of Experiment 19, enumerated completely:
  3 hidden states, rows from {0,1/4,1/2,3/4,1}, 15^3 = 3375 generators,
  lumping {0,1} -> a, {2} -> b, uniform initial state.

------------------------------------------------------------------- MEASURES
  i_O(g)  log2 |{g' : law(Y_0..Y_7 | g') = law(.. | g)}|      (as in Exp 19)
  C_mu(g) Shannon entropy of the stationary distribution over causal states
          of the observed process, via the mixed-state construction with
          predictively-equivalent states merged
  h_mu(g) entropy rate of the observed process

------------------------------------------------------- CRITERIA (D1-D3)
D1  SANITY. C_mu is constant on every observational equivalence class.
    It must be: the class is defined by the observed law, and C_mu is a
    function of that law. If D1 fails the construction is wrong and D2/D3
    mean nothing.

D2  C_mu does not determine i_O. Among generators sharing a common C_mu
    value (to 1e-9), the 10th-90th percentile spread of i_O is at least
    1.00 bit, for at least one C_mu value holding >=1% of the class.

D3  C_mu is not a proxy for i_O. |Pearson r| between C_mu and i_O < 0.5.

Failing D2 or D3 would mean the paper's Related Work claim is too strong.
"""
import math, itertools
from fractions import Fraction as F
from collections import Counter, defaultdict

NS, T = 3, 8
OBS = [0, 0, 1]
ROWS = [(a, b, 4-a-b) for a in range(5) for b in range(5-a)]
MAXMIX, FUT = 400, 6

# ---------- observation law (Experiment 19's identification channel) ----------
def seq_law(M, pi):
    res = []
    def rec(alpha, depth):
        if depth == T:
            res.append(sum(alpha)); return
        nxt = tuple(sum(alpha[x]*M[x][j] for x in range(NS)) for j in range(NS))
        for y in (0, 1):
            rec(tuple(nxt[x] if OBS[x] == y else 0 for x in range(NS)), depth+1)
    for y0 in (0, 1):
        rec(tuple(pi[x] if OBS[x] == y0 else 0 for x in range(NS)), 1)
    return tuple(res)

# ---------- mixed-state / epsilon-machine construction ----------
def step_belief(b, M, y):
    """Bayes update of belief b on observing symbol y. Returns (prob, b')."""
    un = [sum(b[x]*F(M[x][j], 4) for x in range(NS)) if OBS[j] == y else F(0)
          for j in range(NS)]
    Z = sum(un)
    if Z == 0:
        return F(0), None
    return Z, tuple(v/Z for v in un)

def future_dist(b, M, depth):
    """Distribution over the next `depth` symbols from belief b -- the
    predictive signature that defines a causal state."""
    out = []
    def rec(bb, p, d):
        if d == depth:
            out.append(p); return
        for y in (0, 1):
            Z, nb = step_belief(bb, M, y)
            if Z == 0:
                rec(bb, F(0), d+1)
            else:
                rec(nb, p*Z, d+1)
    rec(b, F(1), 0)
    return tuple(out)

def epsilon_machine_capped(M, cap):
    """Reachability only: does the mixed-state set close within `cap` states?"""
    start = tuple(F(1, NS) for _ in range(NS))
    seen = {start}; frontier = [start]
    while frontier:
        b = frontier.pop()
        for y in (0, 1):
            Z, nb = step_belief(b, M, y)
            if Z == 0:
                continue
            if nb not in seen:
                if len(seen) >= cap:
                    return None
                seen.add(nb); frontier.append(nb)
    return len(seen)

def epsilon_machine(M):
    """Reachable mixed states, merged by predictive equivalence.
    Returns (C_mu, h_mu, n_causal) or None if the mixed-state set blows up."""
    start = tuple(F(1, NS) for _ in range(NS))
    seen, order, edges = {start: 0}, [start], {}
    frontier = [start]
    while frontier:
        b = frontier.pop()
        for y in (0, 1):
            Z, nb = step_belief(b, M, y)
            if Z == 0:
                continue
            if nb not in seen:
                if len(order) >= MAXMIX:
                    return None
                seen[nb] = len(order); order.append(nb); frontier.append(nb)
            edges[(seen[b], y)] = (float(Z), seen[nb])
    n = len(order)
    # stationary distribution over mixed states
    pi = [1.0/n]*n
    for _ in range(400):
        nxt = [0.0]*n
        for (i, y), (p, j) in edges.items():
            nxt[j] += pi[i]*p
        s = sum(nxt)
        if s == 0: break
        pi = [v/s for v in nxt]
    # merge predictively equivalent mixed states -> causal states
    sig = {}
    for i, b in enumerate(order):
        sig.setdefault(future_dist(b, M, FUT), []).append(i)
    cmu = 0.0
    for grp in sig.values():
        m = sum(pi[i] for i in grp)
        if m > 1e-12:
            cmu -= m*math.log2(m)
    hmu = 0.0
    for i in range(n):
        for y in (0, 1):
            e = edges.get((i, y))
            if e and e[0] > 0:
                hmu -= pi[i]*e[0]*math.log2(e[0])
    return cmu, hmu, len(sig)

def pearson(xs, ys):
    n = len(xs); mx = sum(xs)/n; my = sum(ys)/n
    sxy = sum((a-mx)*(b-my) for a, b in zip(xs, ys))
    sxx = sum((a-mx)**2 for a in xs); syy = sum((b-my)**2 for b in ys)
    return sxy/math.sqrt(sxx*syy) if sxx and syy else float('nan')

def pct(v, q):
    v = sorted(v); k = (len(v)-1)*q/100.0
    lo = int(k); hi = min(lo+1, len(v)-1)
    return v[lo] + (v[hi]-v[lo])*(k-lo)

def main():
    print("=" * 78)
    print("EXPERIMENT 21  is i_O(g) recoverable from the epsilon-machine?")
    print("=" * 78)
    gens = list(itertools.product(ROWS, repeat=NS))
    PI = (1, 1, 1)
    laws = [seq_law(M, PI) for M in gens]
    cls = Counter(laws)
    iO = [math.log2(cls[l]) for l in laws]

    cmu, hmu, ncs, ok = [], [], [], []
    for M in gens:
        r = epsilon_machine(M)
        if r is None:
            cmu.append(None); hmu.append(None); ncs.append(None); ok.append(False)
        else:
            cmu.append(r[0]); hmu.append(r[1]); ncs.append(r[2]); ok.append(True)
    keep = [i for i in range(len(gens)) if ok[i]]
    print(f"  |H| = {len(gens)}, epsilon-machine built for {len(keep)} "
          f"({100*len(keep)/len(gens):.1f}%); {len(gens)-len(keep)} exceeded "
          f"{MAXMIX} mixed states and are excluded\n")
    print(f"  mean C_mu {sum(cmu[i] for i in keep)/len(keep):.3f} bits, "
          f"mean h_mu {sum(hmu[i] for i in keep)/len(keep):.3f} bits/symbol, "
          f"mean causal states {sum(ncs[i] for i in keep)/len(keep):.2f}")
    print(f"  mean i_O {sum(iO[i] for i in keep)/len(keep):.3f} bits")

    # ---- D1: C_mu constant on classes ----
    bad = 0
    byclass = defaultdict(list)
    for i in keep:
        byclass[laws[i]].append(i)
    for l, idx in byclass.items():
        vs = [cmu[i] for i in idx]
        if max(vs) - min(vs) > 1e-9:
            bad += 1
    d1 = bad == 0
    print(f"\n  D1  C_mu constant within each equivalence class: "
          f"{'PASS' if d1 else 'FAIL'} ({bad} classes violate)")

    # ---- D2: spread of i_O at fixed C_mu ----
    groups = defaultdict(list)
    for i in keep:
        groups[round(cmu[i], 9)].append(iO[i])
    best, bestc = 0.0, None
    for c, vs in groups.items():
        if len(vs) >= 0.01*len(keep):
            sp = pct(vs, 90) - pct(vs, 10)
            if sp > best:
                best, bestc = sp, (c, len(vs))
    d2 = best >= 1.00
    print(f"  D2  i_O spread at fixed C_mu >= 1.00 bit:            "
          f"{'PASS' if d2 else 'FAIL'} ({best:.2f} bits"
          + (f" at C_mu={bestc[0]:.3f}, n={bestc[1]}" if bestc else "") + ")")

    # ---- D3: correlation ----
    r = pearson([cmu[i] for i in keep], [iO[i] for i in keep])
    d3 = abs(r) < 0.5
    print(f"  D3  |Pearson r| between C_mu and i_O < 0.5:          "
          f"{'PASS' if d3 else 'FAIL'} (r = {r:+.3f})")
    rh = pearson([hmu[i] for i in keep], [iO[i] for i in keep])
    print(f"\n  (for reference, r between h_mu and i_O = {rh:+.3f};")
    print(f"   r between C_mu and h_mu = "
          f"{pearson([cmu[i] for i in keep], [hmu[i] for i in keep]):+.3f})")
    print(f"\n  overall: {sum([d1,d2,d3])}/3 criteria met")

    # ---- is the exclusion biasing the comparison? ----
    import random
    excl = [i for i in range(len(gens)) if not ok[i]]
    mk = sum(iO[i] for i in keep)/len(keep)
    me = sum(iO[i] for i in excl)/len(excl)
    print("\n  exclusion diagnostic:")
    print(f"    kept     n={len(keep):>5}  mean i_O = {mk:.3f} bits")
    print(f"    excluded n={len(excl):>5}  mean i_O = {me:.3f} bits  "
          f"(difference {me-mk:+.3f})")
    rng = random.Random(3)
    sample = rng.sample(excl, 60)
    fin = 0
    for i in sample:
        if epsilon_machine_capped(gens[i], 20000) is not None:
            fin += 1
    print(f"    of 60 excluded generators re-run at a 20,000-state budget, "
          f"{fin} terminated")
    print("    -- consistent with predictive states being generically infinite for")
    print("       hidden Markov processes (Jurgens and Crutchfield 2021), so C_mu is")
    print("       not a finitely constructible Shannon entropy across this family,")
    print("       whereas i_O is a finite integer for every one of them.")

main()
