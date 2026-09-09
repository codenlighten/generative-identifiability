#!/usr/bin/env python3
"""
Minimal Generative Systems VI: partially observed Markov chains.

A deliberately adversarial fourth system.  The first three are deterministic,
fully observed, and algorithmic; this one is stochastic, partially observed, and
its unpredictability is aleatoric rather than computational.  If the two-axis
separation is an artefact of algorithmic systems, it should fail here.

PRE-REGISTERED CRITERIA (fixed before running; reported honestly either way):

  C1  All four quadrants are non-empty, with at least 1% of the hypothesis
      class in each of the two off-diagonal cells.
  C2  |Pearson correlation between the two axes| < 0.3 across the class,
      i.e. the axes are not measuring the same thing.
  C3  An explicit witness exists: G1 != G2 inducing an identical observation
      law, whose shared predictive entropy is < 0.1 bits per step.

Hidden state X_t in {0,1,2}; transition rows drawn from all probability
vectors with entries in {0,1/4,1/2,3/4,1}, giving |H| = 15^3 = 3375 generators.
Observation lumps states {0,1} -> 'a' and {2} -> 'b'.

  identifiability  I(G) = log2 |{G' : law(Y_0..Y_{T-1} | G') = law(... | G)}|
  predictability   P(G) = H(Y_{T-1} | Y_0..Y_{T-2}, G)   -- generator KNOWN

I is epistemic uncertainty about the mechanism; P is aleatoric uncertainty
about the consequence given the mechanism.  They are different quantities and
the point of the experiment is whether they vary independently.
"""
import math, itertools
from collections import Counter, defaultdict

NS, T = 3, 8
OBS = [0, 0, 1]                      # the lumping
ROWS = [(a, b, 4-a-b) for a in range(5) for b in range(5-a)]   # 15 rows, /4
DEN = 3 * 4**(T-1)                   # common denominator of every sequence prob

def seq_law(M, pi):
    """Exact law of (Y_0..Y_{T-1}) as integer numerators over DEN.
    Adjacent entries differ in the LAST symbol."""
    res = []
    def rec(alpha, depth):
        if depth == T:
            res.append(alpha[0]+alpha[1]+alpha[2]); return
        n0 = alpha[0]*M[0][0] + alpha[1]*M[1][0] + alpha[2]*M[2][0]
        n1 = alpha[0]*M[0][1] + alpha[1]*M[1][1] + alpha[2]*M[2][1]
        n2 = alpha[0]*M[0][2] + alpha[1]*M[1][2] + alpha[2]*M[2][2]
        nxt = (n0, n1, n2)
        for y in (0, 1):
            rec(tuple(nxt[x] if OBS[x] == y else 0 for x in range(NS)), depth+1)
    for y0 in (0, 1):
        rec(tuple(pi[x] if OBS[x] == y0 else 0 for x in range(NS)), 1)
    return tuple(res)

def ent(counts, den):
    h = 0.0
    for p in counts:
        if p:
            q = p/den
            h -= q*math.log2(q)
    return h

def cond_last(law):
    """H(Y_{T-1} | Y_0..Y_{T-2}) in bits."""
    full = ent(law, DEN)
    marg = [law[i]+law[i+1] for i in range(0, len(law), 2)]
    return full - ent(marg, DEN)

def pearson(xs, ys):
    n = len(xs); mx = sum(xs)/n; my = sum(ys)/n
    sxy = sum((a-mx)*(b-my) for a, b in zip(xs, ys))
    sxx = sum((a-mx)**2 for a in xs); syy = sum((b-my)**2 for b in ys)
    return sxy/math.sqrt(sxx*syy) if sxx and syy else float('nan')

def main():
    PI = (1, 1, 1)                                   # uniform initial state
    gens = list(itertools.product(ROWS, repeat=NS))
    print("=" * 78)
    print("EXPERIMENT 19  partially observed Markov chains")
    print("=" * 78)
    print(f"  |H| = {len(gens)} generators ({math.log2(len(gens)):.2f} bits), "
          f"{NS} hidden states, binary observation, horizon T = {T}\n")

    laws, preds = [], []
    for M in gens:
        law = seq_law(M, PI)
        laws.append(law); preds.append(cond_last(law))
    classes = Counter(laws)
    ident = [math.log2(classes[l]) for l in laws]

    HGO = sum(v*math.log2(v) for v in classes.values())/len(gens)
    print(f"  H(G|O) over the whole class          : {HGO:.2f} bits")
    print(f"  distinct observation laws            : {len(classes)}")
    print(f"  largest indistinguishable class      : {max(classes.values())} generators")
    print(f"  mean predictive entropy H(Y|past,G)  : {sum(preds)/len(preds):.2f} bits/step")

    # ---- quadrants -------------------------------------------------------
    PT = 0.10          # bits/step below which we call the process predictable
    q = Counter()
    for i in range(len(gens)):
        q[(ident[i] == 0.0, preds[i] < PT)] += 1
    n = len(gens)
    print(f"\n  quadrants (identifiable = 0 bits unresolved; "
          f"predictable = < {PT} bits/step):\n")
    print(f"  {'':<20}{'predictable':>16}{'not predictable':>18}")
    for I, lab in ((True, "identifiable"), (False, "not identifiable")):
        print(f"  {lab:<20}{q[(I,True)]:>10} ({100*q[(I,True)]/n:4.1f}%)"
              f"{q[(I,False)]:>10} ({100*q[(I,False)]/n:4.1f}%)")

    r = pearson(ident, preds)
    print(f"\n  Pearson correlation between the two axes: r = {r:+.3f}")

    # ---- witness for the non-identifiable / predictable cell -------------
    bylaw = defaultdict(list)
    for i, l in enumerate(laws):
        bylaw[l].append(i)
    wit = None
    for l, idxs in bylaw.items():
        if len(idxs) > 1 and cond_last(l) < 0.10:
            wit = idxs; break
    print("\n  witness, distinct mechanisms with identical observable law:")
    if wit:
        a, b = gens[wit[0]], gens[wit[1]]
        fmt = lambda M: " ".join("[" + " ".join(f"{v}/4" for v in row) + "]" for row in M)
        print(f"    G1 = {fmt(a)}")
        print(f"    G2 = {fmt(b)}")
        print(f"    {len(wit)} generators share this law; predictive entropy "
              f"{cond_last(laws[wit[0]]):.3f} bits/step, so the OBSERVABLE future is")
        print(f"    pinned down while the MECHANISM is {math.log2(len(wit)):.2f} bits unresolved.")
    else:
        print("    none found")

    # ---- does intervention help? ----------------------------------------
    print("\n  intervention: reset the hidden state instead of accepting nature's:")
    pooled = {}
    for M in gens:
        pooled[M] = tuple(seq_law(M, tuple(3 if x == k else 0 for x in range(NS)))
                          for k in range(NS))
    pc = Counter(pooled.values())
    HGO_i = sum(v*math.log2(v) for v in pc.values())/len(gens)
    print(f"    passive (uniform start)            : {HGO:.2f} bits, "
          f"{len(classes)} distinguishable laws")
    print(f"    intervening (pool all three resets): {HGO_i:.2f} bits, "
          f"{len(pc)} distinguishable laws")
    print(f"    intervention recovers {HGO-HGO_i:.2f} bits")

    # ---- verdict ---------------------------------------------------------
    off1 = q[(True, False)]/n; off2 = q[(False, True)]/n
    c1 = q[(True,True)] and q[(True,False)] and q[(False,True)] and q[(False,False)] \
         and off1 >= 0.01 and off2 >= 0.01
    c2 = abs(r) < 0.3
    c3 = wit is not None
    print("\n" + "=" * 78)
    print("PRE-REGISTERED CRITERIA")
    print("=" * 78)
    print(f"  C1  all four quadrants non-empty, off-diagonals >= 1%   "
          f"{'PASS' if c1 else 'FAIL'}  ({off1:.1%}, {off2:.1%})")
    print(f"  C2  |r| < 0.3 between the axes                          "
          f"{'PASS' if c2 else 'FAIL'}  (r = {r:+.3f})")
    print(f"  C3  witness: same law, distinct mechanism, low entropy  "
          f"{'PASS' if c3 else 'FAIL'}")

    # =================================================================
    # POST-HOC.  Not pre-registered.  Reported as exploratory.
    #
    # C1 fails because under PASSIVE observation the identifiable row is
    # essentially empty: a binary lumping of a 3-state chain from a uniform
    # start destroys too much for any generator to be pinned down.  The
    # question this raises is whether the empty row is a property of the
    # SYSTEM or of the ADMISSIBLE OBSERVATION CLASS.
    # =================================================================
    print()
    print("=" * 78)
    print("POST-HOC (exploratory): does quadrant occupancy depend on the")
    print("                        admissible observation class?")
    print("=" * 78)
    identP = ident
    identI = [math.log2(pc[pooled[M]]) for M in gens]
    for name, idt in (("passive (uniform start)", identP),
                      ("intervening (reset the hidden state)", identI)):
        qq = Counter((idt[i] == 0.0, preds[i] < PT) for i in range(n))
        r2 = pearson(idt, preds)
        print(f"\n  {name}   r = {r2:+.3f}")
        print(f"  {'':<20}{'predictable':>20}{'not predictable':>20}")
        for I, lab in ((True, "identifiable"), (False, "not identifiable")):
            print(f"  {lab:<20}"
                  f"{qq[(I,True)]:>11} ({100*qq[(I,True)]/n:4.1f}%)"
                  f"{qq[(I,False)]:>11} ({100*qq[(I,False)]/n:4.1f}%)")
    print("\n  robustness: relaxing the identifiability threshold barely moves the")
    print("  identifiable row under intervention, so the effect is not an artefact")
    print("  of demanding exactly zero bits:")
    print(f"    {'I(G) <=':>9}{'id+pred':>10}{'id+unpred':>12}")
    for th in (0.0, 0.5, 1.0, 2.0):
        a = sum(1 for i in range(n) if identI[i] <= th and preds[i] < PT)
        b = sum(1 for i in range(n) if identI[i] <= th and preds[i] >= PT)
        print(f"    {th:>9.1f}{a:>10}{b:>12}")

main()
