#!/usr/bin/env python3
"""
Minimal Generative Systems IV: identifiability curves.

One axis for three systems.  For a hypothesis space H of generators and an
observation O, the preimage is  P(O) = { h in H : h consistent with O }.
We plot  log2 |P(O)|  -- the bits of generator identity still unresolved --
against how much / how well we observed.
"""
import random, itertools, math
from collections import defaultdict

def bar(v, vmax, width=34, ch="#"):
    n = 0 if vmax <= 0 else int(round(width * v / vmax))
    return ch * n + "." * (width - n)

# ====================================================================
# SYSTEM A -- elementary cellular automata.  |H| = 256, 8 bits.
# ====================================================================
def eca_observe(rule, W, T, density, seed=0):
    """Return the set of neighbourhood indices actually exercised."""
    t = [(rule >> k) & 1 for k in range(8)]
    rng = random.Random(seed)
    if density is None:                       # single-cell seed
        row = [0] * W; row[W // 2] = 1
    else:
        row = [1 if rng.random() < density else 0 for _ in range(W)]
    seen = set()
    for _ in range(T):
        nxt = [0] * W
        for i in range(W):
            nb = row[(i-1) % W]*4 + row[i]*2 + row[(i+1) % W]
            seen.add(nb); nxt[i] = t[nb]
        row = nxt
    return seen

def eca_curve():
    W = 64
    designs = [("single cell", None), ("density 0.02", 0.02), ("density 0.10", 0.10),
               ("density 0.25", 0.25), ("density 0.50", 0.50)]
    Ts = [1, 2, 3, 5, 10, 25, 100, 400]
    print("=" * 78)
    print("EXPERIMENT 13  ECA: unresolved generator bits vs observation length")
    print("              |H| = 256 rules, so 8 bits to resolve")
    print("=" * 78)
    print(f"{'initial condition':<16}" + "".join(f"{'T='+str(t):>8}" for t in Ts))
    curves = {}
    for name, d in designs:
        row = []
        for T in Ts:
            tot = 0
            for r in range(256):
                cov = len(eca_observe(r, W, T, d))
                tot += 8 - cov                      # log2 preimage
            row.append(tot / 256)
        curves[name] = row
        print(f"{name:<16}" + "".join(f"{v:>8.2f}" for v in row))
    print("\n  same data, plotted (bar length = unresolved bits, max 8):")
    for name, row in curves.items():
        print(f"    {name:<14} T=400  {bar(row[-1], 8)}  {row[-1]:.2f} bits")
    plateau = curves["single cell"]
    print(f"\n  single-cell seed: {plateau[2]:.2f} bits unresolved at T=3, "
          f"{plateau[-1]:.2f} at T=400  ->  running 130x longer gains "
          f"{plateau[2]-plateau[-1]:.2f} bits.")
    stuck = [r for r in range(256) if len(eca_observe(r, W, 400, None)) < 8]
    print(f"  {len(stuck)}/256 rules are permanently unidentifiable from a "
          f"single-cell seed, at ANY observation length.")
    return curves

# ====================================================================
# SYSTEM B -- polynomials.  H = bounded-degree, bounded-coefficient Z[x].
# ====================================================================
def poly_curve():
    print()
    print("=" * 78)
    print("EXPERIMENT 14  polynomials: preimage vs number of evaluation points")
    print("=" * 78)
    for DEG, CO in [(4, 2), (6, 2)]:
        H = list(itertools.product(range(-CO, CO+1), repeat=DEG+1))
        bits = math.log2(len(H))
        print(f"\n  H = degree<={DEG}, coeffs in [-{CO},{CO}]   |H| = {len(H)}  "
              f"({bits:.1f} bits)")
        for label, pts in [("consecutive x=0,1,2,...", list(range(0, 12))),
                           ("all-same x=1,1,1,...",    [1]*12),
                           ("two points repeated",     [0,1]*6)]:
            row = []
            for k in range(0, 13):
                use = pts[:k]
                groups = defaultdict(int)
                for h in H:
                    key = tuple(sum(c*x**i for i, c in enumerate(h)) for x in use)
                    groups[key] += 1
                # H(G|O) = (1/|H|) sum_c n_c log2 n_c -- expected log preimage,
                # NOT log of the mean preimage (Jensen: the latter is an upper bound)
                row.append(sum(v*math.log2(v) for v in groups.values()) / len(H))
            print(f"    {label:<24}" + "".join(f"{v:>6.2f}" for v in row))
        print(f"    {'k =':<24}" + "".join(f"{k:>6}" for k in range(13)))
        print(f"    Lagrange bound for unbounded coeffs would be deg+1 = {DEG+1};")
        print(f"    the measured threshold is lower because H is coefficient-bounded")
        print(f"    (see experiment 14b).")

# ====================================================================
# SYSTEM C -- hypergraph rewriting.  H = a generated rule library.
# ====================================================================
MAXE = 400
def find_matches(edges, lhs):
    idx = defaultdict(set)
    for e in edges:
        for v in e: idx[v].add(e)
    out = []
    def rec(k, bind, used):
        if len(out) > 3000: return
        if k == len(lhs):
            out.append((dict(bind), tuple(used))); return
        pat = lhs[k]
        bound = [bind[v] for v in pat if v in bind]
        cand = set.intersection(*[idx[b] for b in bound]) if bound else edges
        for e in cand:
            if len(e) != len(pat) or e in used: continue
            nb, ok = dict(bind), True
            for pv, ev in zip(pat, e):
                if nb.setdefault(pv, ev) != ev: ok = False; break
            if ok: rec(k+1, nb, used + [e])
    rec(0, {}, [])
    return out

def step(edges, lhs, rhs):
    matches = find_matches(edges, lhs)
    consumed, produced = set(), []
    for b, used in matches:
        if any(e in consumed for e in used): continue
        consumed.update(used); produced.append(b)
    if not produced: return edges
    nxt = set(e for e in edges if e not in consumed)
    fresh = max((max(e) for e in edges), default=0) + 1
    for b in produced:
        local = dict(b)
        for pat in rhs:
            for pv in pat:
                if pv not in local: local[pv] = fresh; fresh += 1
            nxt.add(tuple(local[pv] for pv in pat))
    return nxt

def wl(edges, rounds=3):
    col = {v: 0 for e in edges for v in e}
    for _ in range(rounds):
        sig = defaultdict(list)
        for e in edges:
            cs = tuple(col[v] for v in e)
            for pos, v in enumerate(e): sig[v].append((len(e), pos, cs))
        col = {v: hash((col[v], tuple(sorted(sig[v])))) for v in col}
    return hash((tuple(sorted(col.values())),
                 tuple(sorted(tuple(col[v] for v in e) for e in edges))))

def degseq(edges):
    d = defaultdict(int)
    for e in edges:
        for v in e: d[v] += 1
    return tuple(sorted(d.values()))

def rewrite_curve():
    print()
    print("=" * 78)
    print("EXPERIMENT 15  rewriting: preimage vs richness of the observation")
    print("=" * 78)
    VS = ["x", "y", "z"]
    pairs = [p for p in itertools.product(VS, repeat=2)]
    lib = []
    for size in (1, 2, 3):
        for rhs in itertools.combinations(pairs, size):
            lib.append(list(rhs))
    lhs = [("x", "y")]
    print(f"  |H| = {len(lib)} rules with LHS {{xy}} and RHS of 1-3 edges over x,y,z"
          f"  ({math.log2(len(lib)):.1f} bits)")
    traces = []
    for rhs in lib:
        st, Es, Vs, Hs, Ds = {(1,2),(2,3),(3,1)}, [], [], [], []
        for _ in range(5):
            st = step(st, lhs, rhs)
            Es.append(len(st)); Vs.append(len({v for e in st for v in e}))
            Hs.append(wl(st)); Ds.append(degseq(st))
            if len(st) > MAXE: break
        traces.append((Es, Vs, Hs, Ds))
    ladder = [
        ("final |E| only",            lambda t: (t[0][-1],)),
        ("|E| trajectory",            lambda t: tuple(t[0])),
        ("|E| and |V| trajectory",    lambda t: (tuple(t[0]), tuple(t[1]))),
        ("final degree sequence",     lambda t: (t[3][-1],)),
        ("final graph (WL hash)",     lambda t: (t[2][-1],)),
        ("full graph trajectory",     lambda t: tuple(t[2])),
    ]
    print(f"\n  {'observation':<26}{'classes':>9}{'mean preimage':>15}{'bits left':>11}")
    for label, f in ladder:
        g = defaultdict(int)
        for t in traces: g[f(t)] += 1
        hgo = sum(v*math.log2(v) for v in g.values()) / len(traces)
        mean_pre = sum(v*v for v in g.values()) / len(traces)
        print(f"  {label:<26}{len(g):>9}{mean_pre:>15.2f}{hgo:>11.2f}"
              f"   {bar(hgo, math.log2(len(lib)), 20)}")
    print("  (deterministic scheduler seed throughout, so the order-dependence")
    print("   measured in Experiment 10 is held fixed here, not eliminated.)")

# ====================================================================
def summary(eca):
    print()
    print("=" * 78)
    print("SUMMARY  the identifiability curve is a property of the OBSERVATION")
    print("=" * 78)
    rows = [
        ("ECA, single-cell seed, T=400",   eca["single cell"][-1],  8.0),
        ("ECA, density 0.5, T=3",          eca["density 0.50"][2],  8.0),
        ("polynomials, k=1 point",         None, None),
    ]
    print(f"  {'setting':<34}{'bits unresolved':>17}{'of':>6}")
    for name, v, tot in rows[:2]:
        print(f"  {name:<34}{v:>17.2f}{tot:>6.0f}")
    print("\n  Same objects. Same mathematics. The difference is what you looked at.")

if __name__ == "__main__":
    e = eca_curve()
    poly_curve()
    rewrite_curve()
    summary(e)
