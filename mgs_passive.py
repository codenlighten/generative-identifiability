#!/usr/bin/env python3
"""
Minimal Generative Systems V: identifiability without intervention.

Every identifiability win measured so far came from SETTING the system's state
-- choosing a random initial row, choosing evaluation points.  Here the observer
may not choose anything.  Nature picks the initial condition; the observer
arrives, possibly late, and watches.
"""
import random, math, itertools
from collections import Counter, defaultdict

N = 12
STATES = 1 << N

# neighbourhoods physically PRESENT in a ring state (rule-independent)
NBMASK = [0] * STATES
for s in range(STATES):
    m = 0
    for i in range(N):
        m |= 1 << (((s >> ((i-1) % N)) & 1)*4 + ((s >> i) & 1)*2 + ((s >> ((i+1) % N)) & 1))
    NBMASK[s] = m

POP = [bin(x).count("1") for x in range(256)]

def successor_table(r):
    t = [(r >> k) & 1 for k in range(8)]
    succ = [0]*STATES
    for s in range(STATES):
        o = 0
        for i in range(N):
            nb = ((s >> ((i-1) % N)) & 1)*4 + ((s >> i) & 1)*2 + ((s >> ((i+1) % N)) & 1)
            if t[nb]:
                o |= 1 << i
        succ[s] = o
    return succ

def orbit_analysis(succ):
    """For every start state: the mask of neighbourhoods its whole forward orbit
    exercises, the mask its ATTRACTOR alone exercises, and the attractor period."""
    orb = [None]*STATES        # transient + attractor
    cyc = [None]*STATES        # attractor only
    per = [0]*STATES
    for s0 in range(STATES):
        if orb[s0] is not None:
            continue
        path, pos, s = [], {}, s0
        while orb[s] is None and s not in pos:
            pos[s] = len(path); path.append(s); s = succ[s]
        if s in pos:                                  # closed a fresh cycle
            j = pos[s]; ring = path[j:]
            cm = 0
            for v in ring: cm |= NBMASK[v]
            p = len(ring)
            for v in ring: orb[v], cyc[v], per[v] = cm, cm, p
            tail = path[:j]
        else:                                          # ran into known territory
            tail = path
        for v in reversed(tail):
            nx = succ[v]
            orb[v] = NBMASK[v] | orb[nx]
            cyc[v] = cyc[nx]
            per[v] = per[nx]
    return orb, cyc, per

def main():
    print("=" * 76)
    print("EXPERIMENT 16  cellular automata you are not allowed to restart")
    print(f"              {N}-cell ring, exact over all {STATES} initial states x 256 rules")
    print("=" * 76)

    best_one_step = min(8 - POP[NBMASK[s]] for s in range(STATES))
    mean_one_step = sum(8 - POP[NBMASK[s]] for s in range(STATES)) / STATES

    tot_early = tot_late = 0.0
    zero_early = zero_late = 0
    per_hist = Counter()
    worst = []
    for r in range(256):
        succ = successor_table(r)
        orb, cyc, per = orbit_analysis(succ)
        e = sum(8 - POP[orb[s]] for s in range(STATES))
        l = sum(8 - POP[cyc[s]] for s in range(STATES))
        tot_early += e / STATES
        tot_late  += l / STATES
        zero_early += sum(1 for s in range(STATES) if POP[orb[s]] == 8)
        zero_late  += sum(1 for s in range(STATES) if POP[cyc[s]] == 8)
        for s in range(0, STATES, 8): per_hist[per[s]] += 1
        worst.append((l/STATES, r))

    tot = 256 * STATES
    print("\n  three ways to look at the same 256 rules:\n")
    print(f"  {'regime':<46}{'bits unresolved':>16}")
    print(f"  {'-'*46}{'-'*16}")
    print(f"  {'INTERVENE: choose the state, watch 1 step':<46}{best_one_step:>16.2f}")
    print(f"  {'PASSIVE: nature picks, you arrive at t=0':<46}{tot_early/256:>16.2f}")
    print(f"  {'PASSIVE: nature picks, you arrive late':<46}{tot_late/256:>16.2f}")
    print(f"\n  (a randomly chosen state watched for one step would leave "
          f"{mean_one_step:.2f} bits --\n   intervention is worth {mean_one_step-best_one_step:.2f} bits "
          f"even before the system starts moving)")

    print(f"\n  fully identified (0 bits), as a fraction of all rule/start pairs:")
    print(f"    arriving at t=0 and watching forever : {zero_early/tot:6.1%}")
    print(f"    arriving after the transient is gone : {zero_late/tot:6.1%}")
    print(f"    intervening                          : 100.0%  (guaranteed, in one step)")

    print("\n  attractor periods encountered (the observer's entire future is one period):")
    tp = sum(per_hist.values())
    for p, c in sorted(per_hist.items())[:8]:
        print(f"    period {p:>4}: {c/tp:6.1%}  {'#' * int(round(60*c/tp))}")
    print(f"    ...{len(per_hist)} distinct periods, max {max(per_hist)}")

    worst.sort(reverse=True)
    print("\n  least identifiable rules once the transient is gone:")
    for b, r in worst[:6]:
        print(f"    rule {r:>3}: {b:.2f} bits unresolved on average, "
              f"{2**b:.0f} rules indistinguishable")

    # ---------------------------------------------------------------
    print()
    print("=" * 76)
    print("EXPERIMENT 17  polynomials sampled wherever the system happens to go")
    print("=" * 76)
    DEG, C = 4, 2
    H = list(itertools.product(range(-C, C+1), repeat=DEG+1))
    DOM = list(range(-4, 9))
    col = {x: [sum(c*x**i for i, c in enumerate(h)) for h in H] for x in DOM}
    bitsH = math.log2(len(H))
    print(f"  H = degree<={DEG}, coeffs in [-{C},{C}]  ({len(H)} candidates, {bitsH:.2f} bits)")
    print(f"  The driving system has an attractor of period p, so the observer sees the")
    print(f"  polynomial at p distinct points -- and never at a {len(DOM)+1}th, however long they watch.\n")

    def bits_for(points):
        g = Counter(zip(*[col[x] for x in points]))
        return sum(v*math.log2(v) for v in g.values()) / len(H)   # H(G|O)

    rng = random.Random(11)
    print(f"  {'p':>3}{'chosen points':>16}{'nature-given points':>34}")
    print(f"  {'':>3}{'(intervention)':>16}{'mean':>10}{'best':>8}{'worst':>8}{'p(id)':>8}")
    for p in range(1, 7):
        inter = bits_for(list(range(p)))
        vals = []
        for _ in range(400):
            vals.append(bits_for(rng.sample(DOM, p)))
        ident = sum(1 for v in vals if v == 0.0) / len(vals)
        print(f"  {p:>3}{inter:>16.2f}{sum(vals)/len(vals):>10.2f}"
              f"{min(vals):>8.2f}{max(vals):>8.2f}{ident:>8.0%}")
    print("\n  With intervention, p is a threshold: past it you always win.")
    print("  Without, p only sets the odds -- identifiability becomes a lottery")
    print("  over which states the system happened to visit.")

    # ---------------------------------------------------------------
    print()
    print("=" * 76)
    print("WHY MORE TIME DOES NOT HELP")
    print("=" * 76)
    print("  A passive observer's information saturates after ONE attractor period.")
    print("  Past that the system repeats states it has already shown, and every")
    print("  further observation is a duplicate. The cap on identifiability is the")
    print("  DIVERSITY of the attractor, not the DURATION of the observation.")
    med = sorted(per_hist.elements())[tp//2]
    print(f"\n  median attractor period on this ring: {med}")
    print(f"  share of starts landing on a period-1 or period-2 attractor: "
          f"{(per_hist[1]+per_hist[2])/tp:.1%}")

main()


# =====================================================================
# EXPERIMENT 18  removing the clustered-attractor idealisation
#
# Experiment 17 modelled "nature gives you p points" as a uniform random
# p-subset of the domain.  That is generous: it silently assumes the
# observation points are spread out.  Here the points are not modelled at
# all -- they are whatever the attractor of an explicit driving map visits.
# =====================================================================
def experiment_18():
    from functools import lru_cache
    DEG, C = 4, 2
    Hp = list(itertools.product(range(-C, C+1), repeat=DEG+1))
    DOM = list(range(-4, 9))
    COL = {x: tuple(sum(c*x**i for i, c in enumerate(h)) for h in Hp) for x in DOM}
    NH = len(Hp)

    @lru_cache(maxsize=None)
    def bits(points):
        g = Counter(zip(*[COL[x] for x in points]))
        return sum(v*math.log2(v) for v in g.values()) / NH       # H(G|O)

    def attractor(phi, s):
        seen, order = {}, []
        while s not in seen:
            seen[s] = len(order); order.append(s); s = phi[s]
        return tuple(sorted(set(order[seen[s]:]))), tuple(sorted(set(order)))

    rng = random.Random(2027)
    lo, hi = min(DOM), max(DOM)

    def clip(v): return lo if v < lo else hi if v > hi else v

    models = []

    # A. structureless nature: a uniform random map on the domain
    obs_late, obs_early, ps = [], [], []
    for _ in range(3000):
        phi = {x: rng.choice(DOM) for x in DOM}
        s = rng.choice(DOM)
        cyc, orb = attractor(phi, s)
        obs_late.append(bits(cyc)); obs_early.append(bits(orb)); ps.append(len(cyc))
    models.append(("random map on the domain", obs_late, obs_early, ps))

    # B. realistic dissipative nature: integer affine contractions
    obs_late, obs_early, ps = [], [], []
    for num in range(-9, 10):
        a = num / 10.0
        for b in DOM:
            phi = {x: clip(int(round(a*x + b))) for x in DOM}
            for s in DOM:
                cyc, orb = attractor(phi, s)
                obs_late.append(bits(cyc)); obs_early.append(bits(orb)); ps.append(len(cyc))
    models.append(("affine contraction round(ax+b)", obs_late, obs_early, ps))

    # C. the idealisation used in experiment 17, at matched p
    def matched(ps_ref):
        out = []
        for p in ps_ref:
            out.append(bits(tuple(sorted(rng.sample(DOM, min(p, len(DOM)))))))
        return out

    print()
    print("=" * 76)
    print("EXPERIMENT 18  the observation points are the attractor, not a random subset")
    print("=" * 76)
    print(f"  H = deg<={DEG}, coeffs in [-{C},{C}]  ({NH} candidates, {math.log2(NH):.2f} bits)\n")
    print(f"  {'model of nature':<34}{'mean p':>8}{'late':>9}{'early':>9}{'P(id) late':>12}")
    print(f"  {'-'*34}{'-'*8}{'-'*9}{'-'*9}{'-'*12}")
    for name, late, early, ps in models:
        print(f"  {name:<34}{sum(ps)/len(ps):>8.2f}{sum(late)/len(late):>9.2f}"
              f"{sum(early)/len(early):>9.2f}{sum(1 for v in late if v == 0)/len(late):>12.0%}")
        ideal = matched(ps)
        print(f"  {'  ^ idealised: uniform p-subset':<34}{sum(ps)/len(ps):>8.2f}"
              f"{sum(ideal)/len(ideal):>9.2f}{'':>9}"
              f"{sum(1 for v in ideal if v == 0)/len(ideal):>12.0%}")
        gap = sum(late)/len(late) - sum(ideal)/len(ideal)
        print(f"  {'  ^ idealisation was optimistic by':<34}{'':>8}{gap:>9.2f} bits\n")

    name, late, early, ps = models[1]
    print("  where the affine-contraction observer ends up (late, on the attractor):")
    h = Counter(round(v, 1) for v in late)
    for v, c in sorted(h.items())[:9]:
        print(f"    {v:>5.2f} bits : {c/len(late):6.1%}  {'#'*int(round(52*c/len(late)))}")
    per = Counter(ps)
    print(f"\n  attractor sizes seen: " +
          ", ".join(f"p={k} {v/len(ps):.0%}" for k, v in sorted(per.items())[:5]))
    print("  A contraction has a single fixed point, so the late observer sees ONE value")
    print("  forever, and identifiability is decided entirely by WHERE that point sits.")

experiment_18()
