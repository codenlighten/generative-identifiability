#!/usr/bin/env python3
"""
Minimal Generative Systems: bounded-budget reachability experiments.

Question: given a primitive set P and a budget of n primitive tokens,
how many distinct objects are reachable, and where does adding ONE
primitive change the growth class?
"""
from collections import defaultdict

# ---------------------------------------------------------------
# 1. INTEGER REACHABILITY.  Leaves are all the literal 1.
#    Budget n = number of leaves.
# ---------------------------------------------------------------
ADD = ('+', lambda a, b: a + b)
MUL = ('*', lambda a, b: a * b)
SUB = ('-', lambda a, b: a - b)

def integer_reach(ops, max_leaves, cap=10**7):
    reach = {1: {1}}
    for n in range(2, max_leaves + 1):
        s = set()
        for i in range(1, n):
            for a in reach[i]:
                for b in reach[n - i]:
                    for _, f in ops:
                        v = f(a, b)
                        if abs(v) <= cap:
                            s.add(v)
        reach[n] = s
    return reach

def cumulative(reach, max_leaves):
    out, acc = [], set()
    for n in range(1, max_leaves + 1):
        acc |= reach[n]
        out.append(len(acc))
    return out

# ---------------------------------------------------------------
# 2. BOOLEAN CLOSURE.  How many of the 2^(2^k) functions of k
#    variables can a gate set actually build?
# ---------------------------------------------------------------
def bool_closure(k, gates, unary=(), seeds=()):
    rows = 1 << k
    full = (1 << rows) - 1
    varstt = []
    for i in range(k):
        tt = 0
        for r in range(rows):
            if (r >> i) & 1:
                tt |= 1 << r
        varstt.append(tt)
    S = set(varstt) | set(seeds)
    while True:
        new = set()
        for a in S:
            for u in unary:
                v = u(a) & full
                if v not in S:
                    new.add(v)
            for b in S:
                for g in gates:
                    v = g(a, b) & full
                    if v not in S:
                        new.add(v)
        if not new:
            return S
        S |= new

AND = lambda a, b: a & b
OR  = lambda a, b: a | b
XOR = lambda a, b: a ^ b
NOT = lambda a: ~a
NAND= lambda a, b: ~(a & b)

# ---------------------------------------------------------------
# 3. POLYNOMIAL REACHABILITY.  Leaves are x and 1.
# ---------------------------------------------------------------
MAXDEG, MAXCOEF = 5, 200

def ptrim(p):
    while len(p) > 1 and p[-1] == 0:
        p = p[:-1]
    return p

def padd(a, b):
    n = max(len(a), len(b))
    return ptrim(tuple((a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)))

def pmul(a, b):
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            r[i + j] += x * y
    return ptrim(tuple(r))

def ok(p):
    return len(p) - 1 <= MAXDEG and all(abs(c) <= MAXCOEF for c in p)

def poly_reach(ops, max_leaves):
    reach = {1: {(1,), (0, 1)}}          # the literal 1, and x
    for n in range(2, max_leaves + 1):
        s = set()
        for i in range(1, n):
            for a in reach[i]:
                for b in reach[n - i]:
                    for f in ops:
                        p = f(a, b)
                        if ok(p):
                            s.add(p)
        reach[n] = s
    return reach

# ---------------------------------------------------------------
def main():
    B = 18
    print("=" * 68)
    print("EXPERIMENT 1  integer reachability from the literal 1")
    print("=" * 68)
    sets = [("{1,+}", [ADD]), ("{1,+,*}", [ADD, MUL]), ("{1,+,*,-}", [ADD, MUL, SUB])]
    curves = {}
    for name, ops in sets:
        r = integer_reach(ops, B)
        curves[name] = cumulative(r, B)
        if name == "{1,+,*}":
            mincost = {}
            for n in range(1, B + 1):
                for v in r[n]:
                    mincost.setdefault(v, n)
    print(f"{'budget n':>9} " + " ".join(f"{n:>12}" for n, _ in sets))
    for i in range(B):
        print(f"{i+1:>9} " + " ".join(f"{curves[n][i]:>12,}" for n, _ in sets))
    print()
    for name, _ in sets:
        c = curves[name]
        print(f"  {name:<10} n=18 reach {c[-1]:>10,}   ratio c[18]/c[17] = {c[-1]/c[-2]:.3f}")

    print()
    print("=" * 68)
    print("EXPERIMENT 2  Boolean closure, k=3 variables (256 functions total)")
    print("=" * 68)
    gs = [
        ("{AND}",           [AND], [],      []),
        ("{AND,OR}",        [AND, OR], [],  []),
        ("{XOR}",           [XOR], [],      []),
        ("{XOR,AND}",       [XOR, AND], [], []),
        ("{XOR,AND,1}",     [XOR, AND], [], [0xFF]),
        ("{AND,NOT}",       [AND], [NOT],   []),
        ("{NAND}",          [NAND], [],     []),
    ]
    for name, g, u, s in gs:
        c = len(bool_closure(3, g, u, s))
        print(f"  {name:<14} reaches {c:>4} / 256   ({100*c/256:5.1f}%)  {'UNIVERSAL' if c == 256 else ''}")

    print()
    print("=" * 68)
    print("EXPERIMENT 3  polynomial reachability from leaves {x, 1}")
    print(f"              (capped: degree <= {MAXDEG}, |coef| <= {MAXCOEF})")
    print("=" * 68)
    P = 10
    a = cumulative(poly_reach([padd], P), P)
    b = cumulative(poly_reach([padd, pmul], P), P)
    print(f"{'budget n':>9} {'{x,1,+}':>12} {'{x,1,+,*}':>12}")
    for i in range(P):
        print(f"{i+1:>9} {a[i]:>12,} {b[i]:>12,}")
    r2 = poly_reach([padd, pmul], P)
    seen = set()
    degs = defaultdict(int)
    for n in range(1, P + 1):
        for p in r2[n]:
            if p not in seen:
                seen.add(p)
                degs[len(p) - 1] += 1
    print("\n  degrees reached by {x,1,+,*}:", dict(sorted(degs.items())))
    lin = set()
    for n in range(1, P + 1):
        lin |= poly_reach([padd], P)[n]
    print(f"  {{x,1,+}} reaches only degree <= {max(len(p)-1 for p in lin)}")

    print()
    print("=" * 68)
    print("EXPERIMENT 4  the inverse problem G^-1 for {1,+,*}")
    print("=" * 68)
    print("  min number of 1s needed to write each integer (integer complexity):")
    for lo in range(1, 61, 20):
        print("   " + "  ".join(f"{v}:{mincost[v]}" for v in range(lo, lo + 20)))
    worst = max(range(1, 200), key=lambda v: mincost.get(v, 0))
    print(f"\n  most expensive integer under 200: {worst} costs {mincost[worst]} ones")
    ratio = [(v, mincost[v] / (v.bit_length())) for v in range(2, 200)]
    print(f"  cost is NOT monotone in v: "
          f"cost(107)={mincost[107]} > cost(128)={mincost[128]}")

main()
