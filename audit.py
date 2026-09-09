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
SEEN = set()
def C(key, claim, prov, script, got, expected, tol=None):
    """Register one audited claim under an explicit, stable export key.

    Keys are assigned by hand rather than derived from the claim text: a
    derived key silently changes when prose is reworded, and two claims that
    differ only in punctuation collide. Both happened."""
    if key in SEEN:
        raise SystemExit(f"duplicate result key: RESULT_{key}")
    SEEN.add(key)
    CLAIMS.append((key, claim, prov, script, got, expected, tol))

print("running scripts (cached in %s) ..." % CACHE, file=sys.stderr)

# --- E1-E4, algebra ---
C("INT_REACH_ADD_MUL_B18", "reachable integers, {1,+,*}, budget 18", "EXH", "mgs.py",
  grab("mgs.py", r"\{1,\+,\*\}\s+n=18 reach\s+([\d,]+)"), "404")
C("INT_GROWTH_RATIO_ADD_MUL_B18", "growth ratio c[18]/c[17], {1,+,*}", "EXH", "mgs.py",
  grab("mgs.py", r"\{1,\+,\*\}\s+n=18 reach\s+[\d,]+\s+ratio c\[18\]/c\[17\] = ([\d.]+)"),
  "1.37", 0.0051)   # 2-dp rounding of 1.365
C("BOOL_CLOSURE_AND_OR", "Boolean closure {AND,OR} on 3 vars", "EXH", "mgs.py",
  grab("mgs.py", r"\{AND,OR\}\s+reaches\s+(\d+) / 256"), "18")
C("BOOL_CLOSURE_XOR_AND", "Boolean closure {XOR,AND} on 3 vars", "EXH", "mgs.py",
  grab("mgs.py", r"\{XOR,AND\}\s+reaches\s+(\d+) / 256"), "128")
C("POLY_REACH_ADD_MUL_B10", "polynomials reachable, {x,1,+,*}, budget 10", "EXH", "mgs.py",
  grab("mgs.py", r"^\s+10\s+[\d,]+\s+([\d,]+)\s*$"), "12120")
C("POLY_REACH_ADD_ONLY_B10", "polynomials reachable, {x,1,+}, budget 10", "EXH", "mgs.py",
  grab("mgs.py", r"^\s+10\s+([\d,]+)\s+[\d,]+\s*$"), "65")

cost, oeis_ok, drops, trees107 = aux_integer_complexity()
C("INT_COMPLEXITY_107", "integer complexity cost(107)", "AUX", "audit.py", str(cost[107]), "16")
C("INT_COMPLEXITY_128", "integer complexity cost(128)", "AUX", "audit.py", str(cost[128]), "14")
C("INT_COMPLEXITY_DROPS_1_200", "cost drops in [1,200]", "AUX", "audit.py", str(drops), "49")
C("INT_MIN_TREES_107", "minimal ordered expression trees for 107", "AUX", "audit.py", str(trees107), "26624")
C("INT_COMPLEXITY_OEIS_MATCH", "integer complexity matches OEIS A005245, n=1..32", "AUX", "audit.py",
  str(oeis_ok), "True")

# --- E5-E8, cellular automata ---
C("ECA_REVERSIBLE_COUNT", "reversible elementary rules on the ring", "EXH", "mgs_ca.py",
  grab("mgs_ca.py", r"reversible on the ring\): (\d+) rules"), "6")
C("ECA_EQUIV_CLASSES", "ECA equivalence classes under reflection+swap", "EXH", "mgs_ca.py",
  grab("mgs_ca.py", r"collapse to (\d+) equivalence classes"), "88")
C("ECA_RECOVERED_RANDOM_INIT", "rules recovered from a random initial condition", "EXH", "mgs_ca.py",
  grab("mgs_ca.py", r"rule recovered exactly for\s+(\d+)/256"), "256")
C("ECA_COMPRESSED_SPAN", "compressed spacetime span", "EXH", "mgs_ca.py",
  grab("mgs_ca.py", r"ranges \d+ \.\. \d+ bytes \(span (\d+)\)"), "3083")
C("ECA_MUTATION_BIG_MOVE_PCT", "single-bit mutations moving C by >25% of span", "EXH", "mgs_ca.py",
  grab("mgs_ca.py", r"changing C by >25% of the full span:\s*([\d.]+)%"), "10.7", 0.05)
C("ECA_RULE110_COMPRESSED", "rule 110 compressed spacetime size", "EXH", "mgs_ca.py",
  grab("mgs_ca.py", r"^\s+110 -> C=(\d+)"), "1562")

# --- E13-E15, identifiability curves ---
C("ECA_PLATEAU_SINGLECELL_T400", "ECA plateau, single cell, T=400", "EXH", "mgs_ident.py",
  grab("mgs_ident.py", r"^single cell\s+(?:[\d.]+\s+){7}([\d.]+)"), "1.58", 0.005)
C("ECA_PLATEAU_D050_T1", "ECA plateau, density 0.50, T=1", "EXH", "mgs_ident.py",
  grab("mgs_ident.py", r"^density 0.50\s+([\d.]+)"), "0.00", 0.005)
C("ECA_UNIDENT_SINGLECELL_W64", "rules unidentifiable from single-cell seed (W=64)", "EXH", "mgs_ident.py",
  grab("mgs_ident.py", r"(\d+)/256 rules are permanently unidentifiable"), "132")
C("POLY_CURVE_CONSECUTIVE_K3", "polynomial curve, consecutive, k=3", "EXH", "mgs_ident.py",
  grab("mgs_ident.py", r"consecutive x=0,1,2,\.\.\.\s+[\d.]+\s+[\d.]+\s+[\d.]+\s+([\d.]+)"),
  "1.06", 0.005)
C("REWRITE_RESIDUAL_FULL_TRAJ", "rewriting residual, full graph trajectory", "EXH", "mgs_ident.py",
  grab("mgs_ident.py", r"full graph trajectory\s+\d+\s+[\d.]+\s+([\d.]+)"), "1.09", 0.005)
C("REWRITE_CLASS_SIZE", "rewriting hypothesis class size", "EXH", "mgs_ident.py",
  grab("mgs_ident.py", r"\|H\| = (\d+) rules with LHS"), "129")

# --- E16-E18, passive ---
C("PASSIVE_ECA_INTERVENE", "passive ECA, intervene, one step", "EXH", "mgs_passive.py",
  grab("mgs_passive.py", r"INTERVENE: choose the state, watch 1 step\s+([\d.]+)"), "0.00", 0.005)
C("PASSIVE_ECA_EARLY", "passive ECA, arrive at t=0", "EXH", "mgs_passive.py",
  grab("mgs_passive.py", r"PASSIVE: nature picks, you arrive at t=0\s+([\d.]+)"), "0.51", 0.005)
C("PASSIVE_ECA_LATE", "passive ECA, arrive late", "EXH", "mgs_passive.py",
  grab("mgs_passive.py", r"PASSIVE: nature picks, you arrive late\s+([\d.]+)"), "2.96", 0.005)
C("PASSIVE_ECA_LATE_IDENT_PCT", "fully identified, arriving late", "EXH", "mgs_passive.py",
  grab("mgs_passive.py", r"arriving after the transient is gone :\s*([\d.]+)%"), "20.6", 0.05)
C("PASSIVE_ECA_INTERVENTION_VALUE", "value of intervention before the system moves", "EXH", "mgs_passive.py",
  grab("mgs_passive.py", r"intervention is worth ([\d.]+) bits"), "1.46", 0.005)
C("PASSIVE_PERIOD1_PCT", "starts landing on a period-1 attractor", "EXH", "mgs_passive.py",
  grab("mgs_passive.py", r"period\s+1:\s*([\d.]+)%"), "30.2", 0.05)
C("PASSIVE_PERIOD12_PCT", "starts on period 1 or 2", "EXH", "mgs_passive.py",
  grab("mgs_passive.py", r"period-1 or period-2 attractor:\s*([\d.]+)%"), "42.0", 0.05)
C("PASSIVE_WORST_RULE_BITS", "least identifiable rule, bits unresolved", "EXH", "mgs_passive.py",
  grab("mgs_passive.py", r"rule 255: ([\d.]+) bits"), "7.00", 0.005)
C("PASSIVE_CONTRACTION_BITS_LATE", "affine contraction, bits late", "SMP", "mgs_passive.py",
  grab("mgs_passive.py", r"(?i)affine contraction round\(ax\+b\)\s+[\d.]+\s+([\d.]+)"), "2.56", 0.005)

sp = aux_single_point()
C("POLY_SINGLE_POINT_X1", "single evaluation point, x=1", "AUX", "audit.py", f"{sp[1]:.2f}", "7.90", 0.005)
C("POLY_SINGLE_POINT_X0", "single evaluation point, x=0", "AUX", "audit.py", f"{sp[0]:.2f}", "9.29", 0.005)
C("POLY_SINGLE_POINT_X5", "single evaluation point, x=5 (injective)", "AUX", "audit.py", f"{sp[5]:.2f}", "0.00", 0.005)

# --- E20, pre-registered ---
S = "mgs_preregistered.py"
C("E20_CLASS_SIZE", "pre-registered class size", "EXH", S, grab(S, r"\|H\| = ([\d,]+) generators"), "1500625")
C("E20_CLASS_BITS", "class size in bits", "DRV", S, grab(S, r"\(([\d.]+) bits\), T = 6"), "20.52", 0.005)
C("E20_PASSIVE_HGO", "passive H(G|O)", "EXH", S, grab(S, r"H\(G\|O\) passive\s+:\s*([\d.]+)"), "3.24", 0.005)
C("E20_PASSIVE_LAWS", "passive distinct laws", "EXH", S, grab(S, r"passive\s+:\s*[\d.]+ bits over ([\d,]+)"), "303229")
C("E20_INTERVENE_HGO", "intervening H(G|O)", "EXH", S, grab(S, r"H\(G\|O\) intervening\s+:\s*([\d.]+)"), "1.37", 0.0051)   # 2-dp rounding of 1.365
C("E20_INTERVENE_LAWS", "intervening distinct laws", "EXH", S,
  grab(S, r"intervening\s+:\s*[\d.]+ bits over ([\d,]+)"), "1235809")
C("E20_C2_SPREAD_P1", "C2 spread of p_1", "EXH", S, grab(S, r"C2 .*\(([\d.]+), n="), "0.55", 0.005)
C("E20_C3_SPREAD_IO", "C3 spread of i_O", "EXH", S, grab(S, r"C3 .*\(([\d.]+), n="), "12.00", 0.005)
C("E20_C4_GENERATORS", "C4 qualifying generators", "EXH", S, grab(S, r"C4 .*\(([\d,]+) generators\)"), "81482")
C("E20_C5_MIGRATION", "C5 quadrant migration", "EXH", S, grab(S, r"C5 .*\(([\d.]+)%\)"), "82.0", 0.05)

# --- E22/E23, lattice and horizon (new in Paper 1: T* and the exact migration) ---
L = "mgs_lattice.py"
C("E22_MIXTURE_IDENTITY", "passive law == mean of the four reset laws", "EXH", L,
  grab(L, r"reset laws, exactly: (\w+)"), "CONFIRMED")
C("E22_CONTAINMENT", "A subset B => i_B <= i_A, no violations", "EXH", L,
  grab(L, r"i_B\(g\) <= i_A\(g\) for every generator: (\w+)"), "CONFIRMED")
C("E22_H_FULL_INTERVENTION", "H(G|O_A), full reset set", "EXH", L,
  grab(L, r"^  \{0,1,2,3\}\s+([\d.]+)"), "1.371", 0.005)
C("E22_H_PASSIVE", "H(G|O), passive mixture", "EXH", L,
  grab(L, r"^  passive \(mixture\)\s+([\d.]+)"), "3.245", 0.005)
C("E22_H_SINGLETON", "H(G|O_A), a single reset", "EXH", L,
  grab(L, r"^  \{0\}\s+([\d.]+)"), "4.504", 0.005)
C("E22_H_WITHIN_LUMP", "H(G|O_A), two resets within one lump", "EXH", L,
  grab(L, r"^  \{0,1\}\s+([\d.]+)"), "3.356", 0.005)
C("E22_H_CROSS_LUMP", "H(G|O_A), two resets crossing the lump", "EXH", L,
  grab(L, r"^  \{0,2\}\s+([\d.]+)"), "2.618", 0.005)
C("E22_IDENT_CROSS_LUMP", "identified, two resets crossing the lump", "EXH", L,
  grab(L, r"^  \{0,2\}\s+[\d.]+\s+[\d,]+\s+([\d.]+)%"), "58.55", 0.05)
C("E22_VALUE_OF_INTERVENTION", "value of full intervention over passive", "EXH", L,
  grab(L, r"= ([\d.]+) bits"), "1.874", 0.005)
C("E22_BECAME_IDENTIFIABLE", "generators becoming identifiable", "EXH", L,
  grab(L, r"became identifiable ([\d,]+)"), "1230684")
C("E22_LOST_IDENTIFIABILITY", "generators losing identifiability", "EXH", L,
  grab(L, r"lost identifiability (\d+)"), "0")
C("E22_PREDICTIVE_MOVED", "generators moving in the predictive coordinate", "EXH", L,
  grab(L, r"predictive coordinate moved (\d+)"), "0")
C("E22_SUBMODULAR", "information gain submodular", "EXH", L,
  grab(L, r"diminishing returns\): (\w+)"), "CONFIRMED")
C("E23_TSTAR", "horizon at which the passive partition stabilises", "EXH", L,
  grab(L, r"CONFIRMED, T\* = (\d+)"), "6")
C("E23_CLASSES_T6", "passive equivalence classes at T=6", "EXH", L,
  grab(L, r"^    6\s+[\d.]+\s+([\d,]+)"), "303229")
C("E23_CLASSES_T8", "passive equivalence classes at T=8", "EXH", L,
  grab(L, r"^    8\s+[\d.]+\s+([\d,]+)"), "303229")

# --- E21, the C_mu comparison ---
E = "mgs_epsilon.py"
C("E21_EPSILON_BUILT_PCT", "generators with a finite mixed-state set", "EXH", E,
  grab(E, r"built for \d+ \(([\d.]+)%\)"), "40.8", 0.05)
C("E21_EXCLUDED_MEAN_IO", "mean i_O of excluded generators", "EXH", E,
  grab(E, r"excluded n=\s*\d+\s+mean i_O = ([\d.]+)"), "2.217", 0.005)
C("E21_KEPT_MEAN_IO", "mean i_O of retained generators", "EXH", E,
  grab(E, r"kept\s+n=\s*\d+\s+mean i_O = ([\d.]+)"), "2.325", 0.005)
C("E21_TERMINATED_AT_20K", "excluded generators terminating at a 20,000-state budget",
  "SMP", E, grab(E, r"20,000-state budget, (\d+) terminated"), "0")

# --- E24, the observation-partition sweep behind Paper 2's refutation ---
Q = "mgs_partitions.py"
C("E24_CLASS_SIZE", "partition-sweep hypothesis class size", "EXH", Q,
  grab(Q, r"\|H\| = ([\d,]+) generators"), "10000")
C("E24_X1_CELLS_OK", "cells where H is non-increasing in coverage", "EXH", Q,
  grab(Q, r"(\d+)/\d+ cells = "), "3")
C("E24_X1_CELLS_TOTAL", "cells tested for the crossing conjecture", "EXH", Q,
  grab(Q, r"\d+/(\d+) cells = "), "23")
C("E24_X1_FRACTION", "fraction of cells supporting the crossing conjecture", "EXH", Q,
  grab(Q, r"cells = ([\d.]+)%"), "13.0", 0.05)
C("E24_X1_VERDICT", "pre-registered verdict on the crossing conjecture", "EXH", Q,
  grab(Q, r"cells = [\d.]+%\s+(NOT SUPPORTED|SUPPORTED)"), "NOT SUPPORTED")
C("E24_X2_EQUAL_BLOCK_OK", "equal-block cells satisfying the conjecture", "EXH", Q,
  grab(Q, r"equal-block partitions only: (\d+)/"), "3")
C("E24_X2_EQUAL_BLOCK_TOTAL", "equal-block cells tested", "EXH", Q,
  grab(Q, r"equal-block partitions only: \d+/(\d+)"), "3")
C("E24_X3_MAX_SPREAD", "largest H spread across coverage at fixed |A|", "EXH", Q,
  grab(Q, r"across coverage at fixed \|A\|: ([\d.]+) bits"), "1.065", 0.005)
C("E24_REVERSAL_COV1", "unbalanced partition, coverage 1, H", "EXH", Q,
  grab(Q, r"partition 0\|123, \|A\|=2, cov1:([\d.]+)"), "4.852", 0.005)
C("E24_REVERSAL_COV2", "unbalanced partition, coverage 2, H", "EXH", Q,
  grab(Q, r"partition 0\|123, \|A\|=2, cov1:[\d.]+, cov2:([\d.]+)"), "5.917", 0.005)
C("E24_BALANCED_COV1", "balanced partition, coverage 1, H", "EXH", Q,
  grab(Q, r"^    01\|23\s+\|A\|=2\s+1\s+([\d.]+)"), "4.187", 0.005)
C("E24_BALANCED_COV2", "balanced partition, coverage 2, H", "EXH", Q,
  grab(Q, r"^    01\|23\s+\|A\|=2\s+2\s+([\d.]+)"), "3.799", 0.005)

# --- E25, decomposition of the crossing effect ---
V = "mgs_novelty.py"
C("E25_IDENTITY_PAIRS", "pairs satisfying log2|H| - H(G|O_A) = H(Z_A)", "DRV", V,
  grab(V, r"== H\(Z_A\):\s+(\d+)/\d+ pairs"), "84")
C("E25_H_SINGLETON_BLOCK", "H(Z_j) for a reset into a singleton block", "EXH", V,
  grab(V, r"^  0\|123\s+0\s+1\s+([\d.]+)"), "5.017", 0.005)
C("E25_H_THREE_BLOCK", "H(Z_j) for a reset into a three-state block", "EXH", V,
  grab(V, r"^  0\|123\s+0\s+1\s+[\d.]+.*\n\s+1\s+3\s+([\d.]+)"), "6.028", 0.005)
C("E25_REDUNDANCY_SAME_BLOCK", "mean I(Z_i;Z_j), same observation block", "EXH", V,
  grab(V, r"same observation block\s+mean\s+([\d.]+)"), "4.935", 0.005)
C("E25_REDUNDANCY_CROSS_BLOCK", "mean I(Z_i;Z_j), different blocks", "EXH", V,
  grab(V, r"different blocks\s+mean\s+([\d.]+)"), "5.970", 0.005)
C("E25_REVERSAL_INDIVIDUAL_TERM", "reversal gap, individual-entropy term", "EXH", V,
  grab(V, r"^  0\|123     \{0,1\}.*\n    gap [+-][\d.]+ bits = individual term \+?([-\d.]+)"),
  "1.011", 0.005)
C("E25_REVERSAL_REDUNDANCY_TERM", "reversal gap, redundancy term", "EXH", V,
  grab(V, r"^  0\|123     \{0,1\}.*\n    gap.*redundancy term \+?([-\d.]+)"),
  "0.054", 0.005)
C("E25_BALANCED_INDIVIDUAL_TERM", "balanced gap, individual-entropy term", "EXH", V,
  grab(V, r"^  01\|23     \{0,2\}.*\n    gap [+-][\d.]+ bits = individual term \+?([-\d.]+)"),
  "0.000", 0.005)
C("E25_BALANCED_REDUNDANCY_TERM", "balanced gap, redundancy term", "EXH", V,
  grab(V, r"^  01\|23     \{0,2\}.*\n    gap.*redundancy term \+?([-\d.]+)"),
  "-0.387", 0.005)

# ---------------------------------------------------------------- report
w = max(len(c[1]) for c in CLAIMS)
print("=" * (w + 54))
print("CLAIM-BY-CLAIM NUMERICAL AUDIT")
print("=" * (w + 54))
print(f"{'claim':<{w}}  {'prov':<4} {'expected':>12} {'recomputed':>12}  status")
print("-" * (w + 54))
bad = 0
for key, claim, prov, script, got, exp, tol in CLAIMS:
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
    f.write("# Generated by audit.py. Keys are stable; values are recomputed.\n")
    for key, claim, prov, script, got, exp, tol in CLAIMS:
        f.write(f"RESULT_{key}={got}\n")
print(f"\nmachine-readable values written to RESULTS.txt")
