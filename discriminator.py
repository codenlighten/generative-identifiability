#!/usr/bin/env python3
"""
ModelDiscriminator -- a worked demonstration of the framework as a tool.

The scientific content is deliberately dull: a device with three internal
states and sixteen candidate fault mechanisms. The point is the inference
machinery, not the domain.

WHAT THIS DEMONSTRATES, AND WHAT IT DOES NOT
  It answers three questions in order, and keeps them separate:
    decision    which test should I run next, per unit cost
    scientific  what remains unknowable, and is more passive data any use
    proof       the certificate a technical user can check by hand

  We do NOT claim that detecting non-identifiability is new. Structural
  identifiability analysis, Fisher-information singularity and the
  differential-algebra tools do it, and Bayesian optimal experimental design
  ranks experiments by expected information gain. The narrower combination
  demonstrated here is: given an explicit finite family of forward models and
  an observation regime, determine exactly which distinctions are permanently
  unobservable, certify that after a bounded horizon, and then rank
  interventions by information added beyond what is already known.

THE FORWARD-MODEL CONTRACT
  Every hypothesis must be able to predict the outcome distribution of every
  experiment it is ranked on. A hypothesis that cannot is not silently
  skipped: the experiment is refused. No forward model, no information
  calculation.

TWO MODES, NEVER MIXED
  EXACT       finite class, exact rational laws. May report CERTIFIED.
  APPROXIMATE simulations or estimated laws. May only report ESTIMATED, and
              must never say permanently impossible. Not implemented here;
              the mode banner exists so the distinction is structural.
"""
from fractions import Fraction as F
from collections import defaultdict
import math, itertools

MODE = "EXACT"
NS = 3                                   # device internal states
SENSOR = [0, 0, 1]                       # passive regime: states 0,1 look alike


# ----------------------------------------------------------------- hypotheses
class Hypothesis:
    """A candidate fault mechanism. It is a forward model or it is nothing."""

    def __init__(self, name, P, start, desc=""):
        self.name, self.P, self.start, self.desc = name, P, start, desc

    def law(self, horizon, sensor=None, start=None, P=None):
        """Exact distribution over sensor-reading sequences of `horizon` steps."""
        sensor = SENSOR if sensor is None else sensor
        start = self.start if start is None else start
        P = self.P if P is None else P
        k = max(sensor) + 1
        out = []

        def rec(alpha, d):
            if d == horizon:
                out.append(sum(alpha)); return
            nxt = [sum(alpha[u] * P[u][v] for u in range(NS)) for v in range(NS)]
            for y in range(k):
                rec([nxt[v] if sensor[v] == y else F(0) for v in range(NS)], d + 1)

        for y in range(k):
            rec([start[v] if sensor[v] == y else F(0) for v in range(NS)], 1)
        return tuple(out)

    def predict(self, experiment, horizon):
        """Outcome distribution under an experiment, or None if unmodelled."""
        return experiment.apply(self, horizon)


class Experiment:
    def __init__(self, name, cost, fn, unmodelled=()):
        self.name, self.cost, self.fn, self.unmodelled = name, cost, fn, set(unmodelled)

    def apply(self, h, horizon):
        if h.name in self.unmodelled:
            return None                   # the contract, enforced
        return self.fn(h, horizon)


# ------------------------------------------------------------ the fault family
def build():
    """Sixteen faults, constructed so the family contains three situations a
    naive tool would miss:
      F8, F9, F10  differ only in where the device starts, and F10's law is
                   exactly the mixture of the other two -- invisible to any
                   pairwise comparison of mechanisms
      F0, F3       are observationally identical under every experiment offered
      F11..F13     start already in the fault state, so the evidence excludes them
    """
    def M(a, b, c):
        return [[F(a, 4), F(4 - a, 4), F(0)],
                [F(0), F(b, 4), F(4 - b, 4)],
                [F(4 - c, 4), F(0), F(c, 4)]]

    d0  = [F(1), F(0), F(0)]
    d1  = [F(0), F(1), F(0)]
    d2  = [F(0), F(0), F(1)]
    mix = [F(1, 2), F(1, 2), F(0)]
    shared = M(1, 2, 3)

    H = [
        Hypothesis("F0",  M(1, 1, 1), d0, "bearing wear"),
        Hypothesis("F1",  M(2, 1, 3), d0, "seal leak"),
        Hypothesis("F2",  M(3, 2, 1), d0, "controller drift"),
        Hypothesis("F3",  M(1, 1, 1), d0, "coupling wear"),      # twin of F0
        Hypothesis("F4",  M(2, 3, 1), d0, "sensor bias"),
        Hypothesis("F5",  M(3, 1, 2), d0, "thermal cycling"),
        Hypothesis("F6",  M(1, 3, 2), d0, "supply ripple"),
        Hypothesis("F7",  M(2, 2, 2), d0, "contamination"),
        Hypothesis("F8",  shared, d0,  "stuck-at, starts nominal"),
        Hypothesis("F9",  shared, d1,  "stuck-at, starts degraded"),
        Hypothesis("F10", shared, mix, "intermittent, start uncertain"),
        Hypothesis("F11", M(3, 3, 3), d2, "hard failure"),
        Hypothesis("F12", M(1, 1, 3), d2, "latched trip"),
        Hypothesis("F13", M(3, 1, 1), d2, "open circuit"),
        Hypothesis("F14", M(2, 3, 3), d0, "vibration"),
        Hypothesis("F15", M(3, 3, 2), d0, "lubrication loss"),
    ]
    return H


def experiments():
    d0 = [F(1), F(0), F(0)]
    d2 = [F(0), F(0), F(1)]
    fine = [0, 1, 2]                      # a probe that separates states 0 and 1
    def disable(h, T):
        Q = [row[:] for row in h.P]
        Q[2] = [F(1), F(0), F(0)]         # force the fault state back to nominal
        return h.law(T, P=Q)
    return [
        Experiment("Reset subsystem A", 20, lambda h, T: h.law(T, start=d0)),
        Experiment("Probe bus B",       50, lambda h, T: h.law(T, sensor=fine)),
        Experiment("Disable module C",  40, disable),
        Experiment("Voltage test D",    10, lambda h, T: h.law(T)),
        Experiment("Inject fault E",    35, lambda h, T: h.law(T, start=d2)),
        Experiment("Clamp rail G",      15, lambda h, T: h.law(T, sensor=fine),
                   unmodelled=("F5", "F6")),
    ]


# ------------------------------------------------------------------- machinery
def entropy_of_partition(labels):
    c = defaultdict(int)
    for l in labels: c[l] += 1
    n = len(labels)
    return -sum((v / n) * math.log2(v / n) for v in c.values())


def left_kernel(rows):
    """Exact basis of { v : v B = 0 } for integer/rational B."""
    n = len(rows); m = len(rows[0])
    aug = [[F(x) for x in rows[i]] + [F(1 if j == i else 0) for j in range(n)]
           for i in range(n)]
    r = 0
    for c in range(m):
        piv = next((i for i in range(r, n) if aug[i][c] != 0), None)
        if piv is None: continue
        aug[r], aug[piv] = aug[piv], aug[r]
        pv = aug[r][c]; aug[r] = [v / pv for v in aug[r]]
        for i in range(n):
            if i != r and aug[i][c] != 0:
                f = aug[i][c]
                aug[i] = [a - f * b for a, b in zip(aug[i], aug[r])]
        r += 1
        if r == n: break
    return [row[m:] for row in aug[r:]]


def main():
    H = build()
    EXPTS = experiments()
    T = 4
    evidence = None                       # observed trace filter, if any

    print("=" * 66)
    print(f"ModelDiscriminator            mode: {MODE}")
    print("=" * 66)

    # ---- which hypotheses survive the evidence ----
    laws = {h.name: h.law(T) for h in H}
    # evidence: the first telemetry reading was nominal, which rules out any
    # mechanism that starts the device already in the fault state
    def consistent(h):
        return h.start[2] == 0
    live = [h for h in H if consistent(h)]

    by_law = defaultdict(list)
    for h in live: by_law[laws[h.name]].append(h.name)
    amb = entropy_of_partition([laws[h.name] for h in live])

    print("\nINFERENCE STATUS\n")
    print(f"  Candidate mechanisms:      {len(H)}")
    print(f"  Consistent with evidence:  {len(live)}")
    print(f"  Distinguishable classes:   {len(by_law)}")
    print(f"  Mechanism ambiguity:       {amb:.3f} bits")
    for law, names in by_law.items():
        if len(names) > 1:
            print(f"    indistinguishable passively: {', '.join(names)}")

    # ---- what passive observation can never resolve ----
    #
    # The kernel of the hypothesis-level law matrix says which MIXTURES of
    # mechanisms observation cannot separate. Its stabilisation, however, is
    # NOT certified by a plateau at this level: the plateau argument needs a
    # single shared transition operator, and each hypothesis has its own. We
    # certify instead through the joint chain, where the theorem does apply:
    # (mechanism, device state) is one hidden Markov chain whose transition is
    # block diagonal -- the fault does not change -- and whose emission depends
    # only on the device state.
    names = [h.name for h in live]
    B = [list(h.law(T)) for h in live]
    K = left_kernel(B)

    n_joint = len(live) * NS
    rank_C = max(SENSOR) + 1
    bound = n_joint - rank_C + 1

    def joint_kernel_dim(t):
        rows = []
        for h in live:
            for x in range(NS):
                st = [F(1) if u == x else F(0) for u in range(NS)]
                rows.append(list(h.law(t, start=st)))
        return len(left_kernel(rows))

    jdims = [joint_kernel_dim(t) for t in range(1, 9)]
    jstab = next((t for t in range(1, len(jdims)) if jdims[t-1] == jdims[t]), None)

    hdims = [len(left_kernel([list(h.law(t)) for h in live])) for t in range(1, 9)]

    print("\nPASSIVE OBSERVABILITY\n")
    print(f"  Mixture directions hidden, by horizon: " +
          ", ".join(f"T={t}:{d}" for t, d in enumerate(hdims, 1)))
    print(f"  Transiently hidden:  {hdims[0] - hdims[-1]}")
    print(f"  Structurally hidden: {hdims[-1]}")
    if K:
        print("\n  Conclusion (sparsest relations first):")
        for v in sorted(K, key=lambda w: sum(1 for c in w if c != 0))[:3]:
            pos = [(names[i], v[i]) for i in range(len(v)) if v[i] > 0]
            neg = [(names[i], -v[i]) for i in range(len(v)) if v[i] < 0]
            fmt = lambda t: " + ".join(f"{c}*{n}" if c != 1 else n for n, c in t)
            print(f"    no passive trace separates  {fmt(pos)}  from  {fmt(neg)}")
        if len(K) > 3:
            print(f"    ... and {len(K)-3} further relations")
    else:
        print("\n  Passive observation can in principle separate every candidate.")

    # ---- rank the interventions ----
    print("\nINTERVENTION OPTIONS\n")
    base = [laws[h.name] for h in live]
    rows, refused = [], []
    for e in EXPTS:
        preds = [h.predict(e, T) for h in live]
        missing = [h.name for h, p in zip(live, preds) if p is None]
        if missing:
            refused.append((e, missing)); continue
        joint = list(zip(base, preds))
        gain = entropy_of_partition(joint) - entropy_of_partition(base)
        rows.append((e, gain))
    rows.sort(key=lambda r: -r[1])
    print(f"  {'Test':<20}{'Cost':>7}{'New information':>18}{'bits/$100':>12}")
    print("  " + "-" * 55)
    for e, g in rows:
        print(f"  {e.name:<20}{'$'+str(e.cost):>7}{g:>13.3f} bits"
              f"{100*g/e.cost:>12.2f}")
    for e, missing in refused:
        print(f"  {e.name:<20}{'$'+str(e.cost):>7}{'REFUSED':>18}")
        print(f"      cannot rank: {len(missing)} hypotheses supply no outcome "
              f"distribution ({', '.join(missing)})")

    best = rows[0]
    dead = [e.name for e, g in rows if g == 0]
    print(f"\n  Recommended next test: {best[0].name}")
    print(f"  Reason: largest reduction in remaining mechanism ambiguity "
          f"({best[1]:.3f} bits).")
    if dead:
        print(f"  Avoid: {', '.join(dead)} -- redundant with existing evidence "
              f"(0.000 bits).")

    # ---- the proof layer ----
    print("\nCERTIFICATE\n")
    if MODE != "EXACT":
        print("  status:   ESTIMATED -- approximate mode may not certify impossibility.")
    elif not K:
        print("  status:   no kernel; nothing to certify.")
    else:
        v = K[0]
        prod = [sum(v[i] * B[i][j] for i in range(len(B))) for j in range(len(B[0]))]
        print(f"  status:    CERTIFIED (exact rational arithmetic, finite class)")
        print(f"  vector:    v = ({', '.join(str(c) for c in v)})")
        print(f"             over ({', '.join(names)})")
        print(f"  check:     v B_T = 0  ->  {all(x == 0 for x in prod)}")
        print(f"  sum:       sum(v) = {sum(v)}   (probability-preserving direction)")
        print()
        print(f"  permanence is certified through the joint chain, not by the")
        print(f"  hypothesis-level plateau -- the plateau argument requires one")
        print(f"  shared transition operator and each mechanism has its own:")
        print(f"    joint chain:      {len(live)} mechanisms x {NS} device states "
              f"= {n_joint} hidden states")
        print(f"    emission rank:    {rank_C}")
        print(f"    horizon bound:    T_obs <= n - rank(C) + 1 = {bound}")
        print(f"    joint kernel dim: " +
              ", ".join(f"T={t}:{d}" for t, d in enumerate(jdims, 1)))
        if jstab:
            print(f"    joint kernel stabilised at T = {jstab}; one plateau of the")
            print(f"    joint chain is permanent, so the hypothesis-level kernel it")
            print(f"    determines is final too.")
        else:
            print(f"    joint kernel had not plateaued by T = {len(jdims)};")
            print(f"    NOT CERTIFIED -- run to the bound of {bound} to decide.")
        print(f"  meaning:   no passive observation of any length separates these.")
    print()


main()
