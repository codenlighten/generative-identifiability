#!/usr/bin/env python3
"""
Claim-by-claim numerical audit.  No number in the manuscript is trusted:
each is re-derived from the script that is supposed to produce it, or
recomputed here when it was only ever produced ad hoc.

Provenance tags:
  EXH  exhaustive computation over a finite class
  SMP  sampled computation (a fixed seed, not exhaustive)
  DRV  derivation / definitional
  AUX  recomputed inside this audit because no committed script emitted it
"""
import subprocess, re, sys, math, itertools, os
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = "/tmp/mgs_audit"
os.makedirs(CACHE, exist_ok=True)

def run(script):
    out = os.path.join(CACHE, script + ".out")
    if not os.path.exists(out) or os.path.getmtime(out) < os.path.getmtime(os.path.join(HERE, script)):
        print(f"  running {script} ...", file=sys.stderr)
        r = subprocess.run([sys.executable, os.path.join(HERE, script)],
                           capture_output=True, text=True, timeout=3600)
        open(out, "w").write(r.stdout)
    return open(out).read()

OUTS = {}
def text(script):
    if script not in OUTS:
        OUTS[script] = run(script)
    return OUTS[script]

# ---------------------------------------------------------------- AUX recomputes
def aux_integer_complexity():
    reach = {1: {1}}
    for n in range(2, 23):
        s = set()
        for i in range(1, n):
            for a in reach[i]:
                for b in reach[n-i]:
                    s.add(a+b); s.add(a*b)
        reach[n] = {v for v in s if v <= 10**7}
    cost = {}
    for n in range(1, 23):
        for v in reach[n]:
            cost.setdefault(v, n)
    oeis = [1,2,3,4,5,5,6,6,6,7,8,7,8,8,8,8,9,8,9,9,9,10,11,9,10,10,9,10,11,10,11,10]
    matches = [cost[v] for v in range(1, 33)] == oeis
    drops = sum(1 for v in range(1, 200) if cost[v+1] < cost[v])
    memo = {}
    def forms(t, b):
        if b == 1: return 1 if t == 1 else 0
        if (t, b) in memo: return memo[(t, b)]
        c = 0
        for i in range(1, b):
            j = b - i
            for a in reach[i]:
                if a and t % a == 0 and (t // a) in reach[j]:
                    c += forms(a, i) * forms(t // a, j)
                if (t - a) in reach[j]:
                    c += forms(a, i) * forms(t - a, j)
        memo[(t, b)] = c
        return c
    return cost, matches, drops, forms(107, cost[107])

def aux_single_point():
    H = list(itertools.product(range(-2, 3), repeat=5))
    def bits(x):
        g = Counter(sum(c*x**i for i, c in enumerate(h)) for h in H)
        return sum(v*math.log2(v) for v in g.values())/len(H)
    return {x: round(bits(x), 2) for x in (0, 1, 2, 3, 4, 5)}

# ---------------------------------------------------------------- claim table
def grab(script, pattern, group=1):
    m = re.search(pattern, text(script), re.M)
    return m.group(group).replace(",", "") if m else None

CLAIMS = []
def C(claim, prov, script, got, expected, tol=None):
    CLAIMS.append((claim, prov, script, got, expected, tol))

print("running scripts (cached in %s) ..." % CACHE, file=sys.stderr)

# --- E1-E4, algebra ---
C("reachable integers, {1,+,*}, budget 18", "EXH", "mgs.py",
  grab("mgs.py", r"\{1,\+,\*\}\s+n=18 reach\s+([\d,]+)"), "404")
C("growth ratio c[18]/c[17], {1,+,*}", "EXH", "mgs.py",
  grab("mgs.py", r"\{1,\+,\*\}\s+n=18 reach\s+[\d,]+\s+ratio c\[18\]/c\[17\] = ([\d.]+)"),
  "1.37", 0.0051)   # 2-dp rounding of 1.365
C("Boolean closure {AND,OR} on 3 vars", "EXH", "mgs.py",
  grab("mgs.py", r"\{AND,OR\}\s+reaches\s+(\d+) / 256"), "18")
C("Boolean closure {XOR,AND} on 3 vars", "EXH", "mgs.py",
  grab("mgs.py", r"\{XOR,AND\}\s+reaches\s+(\d+) / 256"), "128")
C("polynomials reachable, {x,1,+,*}, budget 10", "EXH", "mgs.py",
  grab("mgs.py", r"^\s+10\s+[\d,]+\s+([\d,]+)\s*$"), "12120")
C("polynomials reachable, {x,1,+}, budget 10", "EXH", "mgs.py",
  grab("mgs.py", r"^\s+10\s+([\d,]+)\s+[\d,]+\s*$"), "65")

cost, oeis_ok, drops, trees107 = aux_integer_complexity()
C("integer complexity cost(107)", "AUX", "audit.py", str(cost[107]), "16")
C("integer complexity cost(128)", "AUX", "audit.py", str(cost[128]), "14")
C("cost drops in [1,200]", "AUX", "audit.py", str(drops), "49")
C("minimal ordered expression trees for 107", "AUX", "audit.py", str(trees107), "26624")
C("integer complexity matches OEIS A005245, n=1..32", "AUX", "audit.py",
  str(oeis_ok), "True")

# --- E5-E8, cellular automata ---
C("reversible elementary rules on the ring", "EXH", "mgs_ca.py",
  grab("mgs_ca.py", r"reversible on the ring\): (\d+) rules"), "6")
C("ECA equivalence classes under reflection+swap", "EXH", "mgs_ca.py",
  grab("mgs_ca.py", r"collapse to (\d+) equivalence classes"), "88")
C("rules recovered from a random initial condition", "EXH", "mgs_ca.py",
  grab("mgs_ca.py", r"rule recovered exactly for\s+(\d+)/256"), "256")
C("compressed spacetime span", "EXH", "mgs_ca.py",
  grab("mgs_ca.py", r"ranges \d+ \.\. \d+ bytes \(span (\d+)\)"), "3083")
C("single-bit mutations moving C by >25% of span", "EXH", "mgs_ca.py",
  grab("mgs_ca.py", r"changing C by >25% of the full span:\s*([\d.]+)%"), "10.7", 0.05)
C("rule 110 compressed spacetime size", "EXH", "mgs_ca.py",
  grab("mgs_ca.py", r"^\s+110 -> C=(\d+)"), "1562")

# --- E13-E15, identifiability curves ---
C("ECA plateau, single cell, T=400", "EXH", "mgs_ident.py",
  grab("mgs_ident.py", r"^single cell\s+(?:[\d.]+\s+){7}([\d.]+)"), "1.58", 0.005)
C("ECA plateau, density 0.50, T=1", "EXH", "mgs_ident.py",
  grab("mgs_ident.py", r"^density 0.50\s+([\d.]+)"), "0.00", 0.005)
C("rules unidentifiable from single-cell seed (W=64)", "EXH", "mgs_ident.py",
  grab("mgs_ident.py", r"(\d+)/256 rules are permanently unidentifiable"), "132")
C("polynomial curve, consecutive, k=3", "EXH", "mgs_ident.py",
  grab("mgs_ident.py", r"consecutive x=0,1,2,\.\.\.\s+[\d.]+\s+[\d.]+\s+[\d.]+\s+([\d.]+)"),
  "1.06", 0.005)
C("rewriting residual, full graph trajectory", "EXH", "mgs_ident.py",
  grab("mgs_ident.py", r"full graph trajectory\s+\d+\s+[\d.]+\s+([\d.]+)"), "1.09", 0.005)
C("rewriting hypothesis class size", "EXH", "mgs_ident.py",
  grab("mgs_ident.py", r"\|H\| = (\d+) rules with LHS"), "129")

# --- E16-E18, passive ---
C("passive ECA, intervene, one step", "EXH", "mgs_passive.py",
  grab("mgs_passive.py", r"INTERVENE: choose the state, watch 1 step\s+([\d.]+)"), "0.00", 0.005)
C("passive ECA, arrive at t=0", "EXH", "mgs_passive.py",
  grab("mgs_passive.py", r"PASSIVE: nature picks, you arrive at t=0\s+([\d.]+)"), "0.51", 0.005)
C("passive ECA, arrive late", "EXH", "mgs_passive.py",
  grab("mgs_passive.py", r"PASSIVE: nature picks, you arrive late\s+([\d.]+)"), "2.96", 0.005)
C("fully identified, arriving late", "EXH", "mgs_passive.py",
  grab("mgs_passive.py", r"arriving after the transient is gone :\s*([\d.]+)%"), "20.6", 0.05)
C("value of intervention before the system moves", "EXH", "mgs_passive.py",
  grab("mgs_passive.py", r"intervention is worth ([\d.]+) bits"), "1.46", 0.005)
C("starts landing on a period-1 attractor", "EXH", "mgs_passive.py",
  grab("mgs_passive.py", r"period\s+1:\s*([\d.]+)%"), "30.2", 0.05)
C("starts on period 1 or 2", "EXH", "mgs_passive.py",
  grab("mgs_passive.py", r"period-1 or period-2 attractor:\s*([\d.]+)%"), "42.0", 0.05)
C("least identifiable rule, bits unresolved", "EXH", "mgs_passive.py",
  grab("mgs_passive.py", r"rule 255: ([\d.]+) bits"), "7.00", 0.005)
C("affine contraction, bits late", "SMP", "mgs_passive.py",
  grab("mgs_passive.py", r"(?i)affine contraction round\(ax\+b\)\s+[\d.]+\s+([\d.]+)"), "2.56", 0.005)

sp = aux_single_point()
C("single evaluation point, x=1", "AUX", "audit.py", f"{sp[1]:.2f}", "7.90", 0.005)
C("single evaluation point, x=0", "AUX", "audit.py", f"{sp[0]:.2f}", "9.29", 0.005)
C("single evaluation point, x=5 (injective)", "AUX", "audit.py", f"{sp[5]:.2f}", "0.00", 0.005)

# --- E20, pre-registered ---
S = "mgs_preregistered.py"
C("pre-registered class size", "EXH", S, grab(S, r"\|H\| = ([\d,]+) generators"), "1500625")
C("class size in bits", "DRV", S, grab(S, r"\(([\d.]+) bits\), T = 6"), "20.52", 0.005)
C("passive H(G|O)", "EXH", S, grab(S, r"H\(G\|O\) passive\s+:\s*([\d.]+)"), "3.24", 0.005)
C("passive distinct laws", "EXH", S, grab(S, r"passive\s+:\s*[\d.]+ bits over ([\d,]+)"), "303229")
C("intervening H(G|O)", "EXH", S, grab(S, r"H\(G\|O\) intervening\s+:\s*([\d.]+)"), "1.37", 0.0051)   # 2-dp rounding of 1.365
C("intervening distinct laws", "EXH", S,
  grab(S, r"intervening\s+:\s*[\d.]+ bits over ([\d,]+)"), "1235809")
C("C2 spread of p_1", "EXH", S, grab(S, r"C2 .*\(([\d.]+), n="), "0.55", 0.005)
C("C3 spread of i_O", "EXH", S, grab(S, r"C3 .*\(([\d.]+), n="), "12.00", 0.005)
C("C4 qualifying generators", "EXH", S, grab(S, r"C4 .*\(([\d,]+) generators\)"), "81482")
C("C5 quadrant migration", "EXH", S, grab(S, r"C5 .*\(([\d.]+)%\)"), "82.0", 0.05)

# ---------------------------------------------------------------- report
w = max(len(c[0]) for c in CLAIMS)
print("=" * (w + 54))
print("CLAIM-BY-CLAIM NUMERICAL AUDIT")
print("=" * (w + 54))
print(f"{'claim':<{w}}  {'prov':<4} {'expected':>12} {'recomputed':>12}  status")
print("-" * (w + 54))
bad = 0
for claim, prov, script, got, exp, tol in CLAIMS:
    if got is None:
        st = "NOT FOUND"; bad += 1
    elif tol is None:
        st = "ok" if got == exp else "MISMATCH"
    else:
        try:
            st = "ok" if abs(float(got) - float(exp)) <= tol else "MISMATCH"
        except ValueError:
            st = "MISMATCH"
    if st != "ok":
        bad += 1 if st != "NOT FOUND" else 0
    print(f"{claim:<{w}}  {prov:<4} {exp:>12} {str(got):>12}  {st}")
print("-" * (w + 54))
print(f"{len(CLAIMS)} claims audited, {len(CLAIMS)-bad} ok, {bad} needing attention")

with open(os.path.join(HERE, "RESULTS.txt"), "w") as f:
    for claim, prov, script, got, exp, tol in CLAIMS:
        key = re.sub(r"[^A-Za-z0-9]+", "_", claim).strip("_").upper()
        f.write(f"RESULT_{key}={got}\n")
print(f"\nmachine-readable values written to RESULTS.txt")
