#!/usr/bin/env python3
"""
Minimal Generative Systems III: hypergraph rewriting.

Tests the proposed hypothesis  C = G(P, R, I)  directly.
The claim is that complexity is a function of primitives, relations and
iteration.  A rewriting system lets us check whether G is a function at all.
"""
import random
from collections import defaultdict

MAXE, STEPS, TRIALS = 1200, 6, 12

# ---------------- matching ----------------
def index_by_vertex(edges):
    idx = defaultdict(set)
    for e in edges:
        for v in e:
            idx[v].add(e)
    return idx

def find_matches(edges, lhs, cap=4000):
    """All (binding, matched_edges) for the pattern list lhs."""
    idx = index_by_vertex(edges)
    out = []
    def rec(k, bind, used):
        if len(out) >= cap:
            return
        if k == len(lhs):
            out.append((dict(bind), tuple(used)))
            return
        pat = lhs[k]
        bound = [bind[v] for v in pat if v in bind]
        cand = set.intersection(*[idx[b] for b in bound]) if bound else edges
        for e in cand:
            if len(e) != len(pat) or e in used:
                continue
            nb = dict(bind)
            ok = True
            for pv, ev in zip(pat, e):
                if nb.setdefault(pv, ev) != ev:
                    ok = False; break
            if ok:
                rec(k + 1, nb, used + [e])
    rec(0, {}, [])
    return out

def apply_step(edges, rules, order_seed):
    """One synchronous step: a maximal non-overlapping match set,
    selected in an order determined by order_seed."""
    rng = random.Random(order_seed)
    allm = []
    for lhs, rhs in rules:
        for b, used in find_matches(edges, lhs):
            allm.append((b, used, rhs))
    rng.shuffle(allm)
    consumed, new, produced = set(), set(edges), []
    for b, used, rhs in allm:
        if any(e in consumed for e in used):
            continue
        consumed.update(used)
        produced.append((b, rhs))
    if not produced:
        return edges, 0
    nxt = set(e for e in edges if e not in consumed)
    fresh = max((max(e) for e in edges), default=0) + 1
    for b, rhs in produced:
        local = dict(b)
        for pat in rhs:
            for pv in pat:
                if pv not in local:
                    local[pv] = fresh; fresh += 1
            nxt.add(tuple(local[pv] for pv in pat))
    return nxt, len(produced)

# ---------------- canonical invariant (sound, incomplete) ----------------
def wl_hash(edges, rounds=3):
    col = {v: 0 for e in edges for v in e}
    for _ in range(rounds):
        sig = defaultdict(list)
        for e in edges:
            cs = tuple(col[v] for v in e)
            for pos, v in enumerate(e):
                sig[v].append((len(e), pos, cs))
        col = {v: hash((col[v], tuple(sorted(sig[v])))) for v in col}
    ec = sorted(tuple(col[v] for v in e) for e in edges)
    return hash((tuple(sorted(col.values())), tuple(ec)))

# ---------------- rule library ----------------
E = lambda *s: [tuple(t.split()) for t in s]
RULES = {
 "identity      {xy}->{xy}":            [(E("x y"), E("x y"))],
 "subdivide     {xy}->{xz}{zy}":        [(E("x y"), E("x z", "z y"))],
 "sprout        {xy}->{xy}{yz}":        [(E("x y"), E("x y", "y z"))],
 "triangle      {xy}->{xy}{yz}{zx}":    [(E("x y"), E("x y", "y z", "z x"))],
 "contract      {xy}{yz}->{xz}":        [(E("x y", "y z"), E("x z"))],
 "star          {xy}{yz}->{xw}{yw}{zw}":[(E("x y", "y z"), E("x w", "y w", "z w"))],
 "reverse       {xy}->{yx}":            [(E("x y"), E("y x"))],
 "fold          {xy}{yz}->{xz}{zy}":    [(E("x y", "y z"), E("x z", "z y"))],
}
SEED = {(1, 2), (2, 3), (3, 1)}

def run(rules, seed_order, steps=STEPS):
    st, sizes = set(SEED), []
    for _ in range(steps):
        if len(st) > MAXE:
            break
        st, n = apply_step(st, rules, seed_order)
        sizes.append(len(st))
        if n == 0:
            break
    return st, sizes

# ---------------- experiments ----------------
def main():
    print("=" * 74)
    print("EXPERIMENT 9  growth class of each rewrite rule (seed = triangle)")
    print("=" * 74)
    print(f"{'rule':<38} {'|E| by step':<26} class")
    growth = {}
    for name, rules in RULES.items():
        _, sizes = run(rules, 0)
        growth[name] = sizes
        if len(sizes) >= 3 and sizes[-1] > 0 and sizes[0] > 0:
            r = (sizes[-1] / sizes[0]) ** (1 / max(1, len(sizes) - 1))
            cls = "exponential" if r > 1.6 else "polynomial" if r > 1.05 else "bounded"
        else:
            cls = "bounded"
        print(f"{name:<38} {str(sizes):<26} {cls}")

    print()
    print("=" * 74)
    print("EXPERIMENT 10  does the OUTPUT depend on the update ORDER?")
    print("            (C = G(P,R,I) assumes it does not)")
    print("=" * 74)
    print(f"{'rule':<38} {'distinct outcomes / ' + str(TRIALS):<22} {'|E| spread'}")
    nonconf = 0
    for name, rules in RULES.items():
        hs, es = set(), []
        for s in range(TRIALS):
            st, sizes = run(rules, s)
            hs.add(wl_hash(st)); es.append(len(st))
        if len(hs) > 1:
            nonconf += 1
        flag = "  <-- NOT order-independent" if len(hs) > 1 else ""
        print(f"{name:<38} {len(hs):<22} {min(es)}..{max(es)}{flag}")
    print(f"\n  {nonconf}/{len(RULES)} rules provably depend on the scheduler.")
    print("  (WL hash is a sound invariant: a difference PROVES non-isomorphism.)")

    print()
    print("=" * 74)
    print("EXPERIMENT 11  the phase transition: add ONE edge to a rule's RHS")
    print("=" * 74)
    base = E("x y")
    variants = [
        ("{xy} -> {xz}",                E("x z")),
        ("{xy} -> {xz}{zy}",            E("x z", "z y")),
        ("{xy} -> {xz}{zy}{yx}",        E("x z", "z y", "y x")),
        ("{xy} -> {xz}{zy}{yx}{zw}",    E("x z", "z y", "y x", "z w")),
    ]
    for label, rhs in variants:
        _, sizes = run([(base, rhs)], 0, steps=8)
        tail = sizes[-1] if sizes else 0
        rate = (sizes[-1] / sizes[0]) ** (1/max(1, len(sizes)-1)) if len(sizes) > 1 else 1
        print(f"  {label:<32} sizes={str(sizes):<34} growth/step={rate:.2f}")

    print()
    print("=" * 74)
    print("EXPERIMENT 12  inverse problem: can the rule be read off the result?")
    print("=" * 74)
    lib = list(RULES.items())
    print(f"  observing only the final hypergraph after {STEPS} steps:")
    amb = []
    for name, rules in lib:
        target, _ = run(rules, 0)
        th = wl_hash(target)
        cons = []
        for oname, orules in lib:
            for s in range(TRIALS):
                st, _ = run(orules, s)
                if wl_hash(st) == th:
                    cons.append(oname); break
        amb.append(len(cons))
        print(f"    {name:<38} consistent rules: {len(cons)}  {cons if len(cons)>1 else ''}")
    print(f"\n  mean preimage size: {sum(amb)/len(amb):.2f} rules")

main()
