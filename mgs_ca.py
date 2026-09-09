#!/usr/bin/env python3
"""
Minimal Generative Systems II: elementary cellular automata.

Regime check.  For {1,+,*} the forward map was cheap and the inverse map
was ill-posed.  ECAs invert the situation: the inverse map is trivial,
and knowing the generator tells you almost nothing.
"""
import zlib
from collections import defaultdict

N = 12                      # ring size for exact state-space work
STATES = 1 << N

def rule_table(r):
    return [(r >> k) & 1 for k in range(8)]

def successor_table(r):
    """succ[s] for every state s of the N-ring, states packed as ints."""
    t = rule_table(r)
    succ = [0] * STATES
    for s in range(STATES):
        o = 0
        for i in range(N):
            l = (s >> ((i - 1) % N)) & 1
            c = (s >> i) & 1
            rr = (s >> ((i + 1) % N)) & 1
            if t[l * 4 + c * 2 + rr]:
                o |= 1 << i
        succ[s] = o
    return succ

def eventual_states(succ):
    """Iterate the image until it stops shrinking: the states on cycles."""
    img = set(range(STATES))
    while True:
        nxt = {succ[s] for s in img}
        if len(nxt) == len(img):
            return nxt
        img = nxt

def attractors(succ):
    """Number of distinct cycles and the longest cycle length."""
    colour, cycles, maxlen = {}, 0, 0
    for s0 in range(STATES):
        if s0 in colour:
            continue
        path, s = {}, s0
        while s not in colour and s not in path:
            path[s] = len(path)
            s = succ[s]
        if s in path:                       # found a fresh cycle
            cycles += 1
            maxlen = max(maxlen, len(path) - path[s])
        for v in path:
            colour[v] = 1
    return cycles, maxlen

# --------------------------------------------------------------------
W, T = 201, 150

def spacetime(r, seed="point"):
    t = rule_table(r)
    if seed == "point":
        row = bytearray(W); row[W // 2] = 1
    else:
        import random; random.seed(7)
        row = bytearray(random.getrandbits(1) for _ in range(W))
    out = bytearray()
    seen = set()
    for _ in range(T):
        out += row
        nxt = bytearray(W)
        for i in range(W):
            l, c, rr = row[(i - 1) % W], row[i], row[(i + 1) % W]
            seen.add(l * 4 + c * 2 + rr)
            nxt[i] = t[l * 4 + c * 2 + rr]
        row = nxt
    return bytes(out), seen

def complexity(r, seed="point"):
    st, _ = spacetime(r, seed)
    return len(zlib.compress(st, 9))

# --------------------------------------------------------------------
def mirror(r):
    t = rule_table(r); out = 0
    for l in (0, 1):
        for c in (0, 1):
            for rr in (0, 1):
                if t[rr * 4 + c * 2 + l]:
                    out |= 1 << (l * 4 + c * 2 + rr)
    return out

def complement(r):
    t = rule_table(r); out = 0
    for l in (0, 1):
        for c in (0, 1):
            for rr in (0, 1):
                if not t[(1 - l) * 4 + (1 - c) * 2 + (1 - rr)]:
                    out |= 1 << (l * 4 + c * 2 + rr)
    return out

def orbits():
    seen, classes = set(), 0
    for r in range(256):
        if r in seen:
            continue
        classes += 1
        frontier = {r}
        while frontier:
            x = frontier.pop()
            if x in seen: continue
            seen.add(x)
            frontier |= {mirror(x), complement(x), mirror(complement(x))} - seen
    return classes

# --------------------------------------------------------------------
def main():
    print("=" * 70)
    print(f"EXPERIMENT 5  exact state space of the {N}-cell ring ({STATES} states)")
    print("=" * 70)
    retention, attr = {}, {}
    for r in range(256):
        succ = successor_table(r)
        ev = eventual_states(succ)
        retention[r] = len(ev)
        attr[r] = attractors(succ)
    rev = [r for r in range(256) if retention[r] == STATES]
    dead = [r for r in range(256) if retention[r] <= 2]
    print(f"  information-preserving (reversible on the ring): {len(rev)} rules  {rev}")
    print(f"  collapse to <=2 eventual states:                  {len(dead)} rules")
    print(f"  median eventual-state count over all 256 rules:   "
          f"{sorted(retention.values())[128]}")
    print("\n  a Hamming-1 pair, same neighbourhood table but one bit apart:")
    for a, b in [(150, 151), (110, 111), (204, 205), (90, 91)]:
        print(f"    rule {a:>3}: {retention[a]:>5} eventual states, {attr[a][0]:>4} attractors"
              f"   |   rule {b:>3}: {retention[b]:>5}, {attr[b][0]:>4}")

    print()
    print("=" * 70)
    print("EXPERIMENT 6  is the generator->complexity landscape smooth?")
    print("=" * 70)
    C = {r: complexity(r) for r in range(256)}
    lo, hi = min(C.values()), max(C.values())
    span = hi - lo
    jumps = []
    for r in range(256):
        for k in range(8):
            jumps.append(abs(C[r] - C[r ^ (1 << k)]))
    jumps.sort()
    print(f"  compressed spacetime size ranges {lo} .. {hi} bytes (span {span})")
    print(f"  over all 2048 single-bit rule mutations:")
    print(f"    median |delta C| = {jumps[len(jumps)//2]}   "
          f"90th pct = {jumps[int(len(jumps)*0.9)]}   max = {jumps[-1]}")
    print(f"    fraction of mutations changing C by >25% of the full span: "
          f"{sum(1 for j in jumps if j > 0.25*span)/len(jumps):.1%}")
    print("\n  rule 110 (universal) and its eight one-bit neighbours:")
    print(f"    110 -> C={C[110]}")
    for k in range(8):
        m = 110 ^ (1 << k)
        print(f"    {m:>3} -> C={C[m]:>5}   ({C[m]/C[110]:.2f}x)")

    print()
    print("=" * 70)
    print("EXPERIMENT 7  the inverse problem: recover the rule from the output")
    print("=" * 70)
    exact = 0
    for r in range(256):
        st, _ = spacetime(r, seed="random")
        rows = [st[i*W:(i+1)*W] for i in range(T)]
        rec, ok = {}, True
        for a, b in zip(rows, rows[1:]):
            for i in range(W):
                nb = a[(i-1) % W]*4 + a[i]*2 + a[(i+1) % W]
                if rec.setdefault(nb, b[i]) != b[i]:
                    ok = False
        if ok and len(rec) == 8 and sum(v << k for k, v in rec.items()) == r:
            exact += 1
    print(f"  from a RANDOM initial condition: rule recovered exactly for "
          f"{exact}/256 rules, in one linear pass.")
    amb = defaultdict(list)
    for r in range(256):
        _, seen = spacetime(r, seed="point")
        amb[8 - len(seen)].append(r)
    print("\n  from a SINGLE-CELL seed, neighbourhoods never exercised -> ambiguity:")
    for miss in sorted(amb):
        print(f"    {miss} unexercised -> preimage of size {2**miss:>4}  "
              f"({len(amb[miss])} rules)")
    tot = sum(len(v) * 2**k for k, v in amb.items()) / 256
    print(f"  mean number of rules consistent with a single-cell-seed run: {tot:.1f}")

    print()
    print("=" * 70)
    print("EXPERIMENT 8  how many genuinely distinct generators are there?")
    print("=" * 70)
    print(f"  256 rule numbers collapse to {orbits()} equivalence classes "
          f"under reflection + colour swap")

main()
