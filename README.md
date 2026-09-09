# Generative Identifiability

Exact measurements of how much of a generating rule an observation actually pins down,
and whether knowing the rule tells you what the system will do.

Companion code and manuscript for *Identifiability and Predictability Are Distinct:
Measured Preimages Across Generative Systems* ([`identifiability.tex`](identifiability.tex),
[`identifiability.pdf`](identifiability.pdf)).

## The question

A recurring proposal in complexity science is to search for the smallest rule that
generates an observed structure. That question is not falsifiable as posed. The smallest
non-trivial generative system available already refuses it: for expressions over
`{1, +, ×}`, the integer 107 costs 16 ones, admits 26,624 minimal expression trees, has no
known closed-form cost function, and its cost is not even monotone — `cost(107) = 16`
while `cost(128) = 14`.

What replaces it is a question about **preimages**. For a hypothesis space `H` of candidate
generators and an observation map `O`,

```
P(O) = { h ∈ H : O(h) = O(g) }
```

with two per-generator coordinates:

```
i_O(g) = log₂|P_O(g)|                  which mechanism is it?
p_k(g) = H(Y_{t+k} | O≤t, G = g)       what happens next, given the mechanism?
```

Identifiability is `|P(O)| = 1`. Everything here is enumerated exactly, not sampled.

## Results

**1. Unresolved information splits in two, and only one part is fixable.**

```
H(G|O) = Δ_O + D_H(𝒜)
```

with `O* = argmin_{O ∈ 𝒜} H(G|O)`, `D_H(𝒜) = H(G|O*)`, `Δ_O = H(G|O) − H(G|O*)`. Both
terms are non-negative by construction; the content is their measured values and the fact
that `D_H` moves when the admissible class `𝒜` does.

| system | \|H\| (bits) | best observation | residual |
|---|---|---|---|
| elementary cellular automata | 8.00 | 0.00 | 0.00 |
| polynomials, deg ≤ 4, C = 2 | 11.61 | 0.00 | 0.00 |
| hypergraph rewriting | 7.01 | 1.09 | **1.09** |

**2. The fixable part is fixable only under intervention.**

Same 256 automata, same ring, varying only what the observer is allowed to do:

| regime | bits unresolved | fully identified |
|---|---|---|
| intervene — choose the state, watch one step | 0.00 | 100% |
| passive — nature picks, arrive at t = 0 | 0.51 | 66.9% |
| passive — nature picks, arrive late | 2.96 | 20.6% |

A passive observer's entire future is one attractor period; 42.0% of initial states land on
a period of 1 or 2. More observation time buys nothing past that.

**3. Identifiability and predictability do not determine one another.**

Tested by a pre-registration frozen before it was run
([`mgs_preregistered.py`](mgs_preregistered.py)): all `35⁴ = 1,500,625` four-state
partially observed Markov chains, enumerated completely. **5/5 criteria met.**

| protocol | id + pred | id + not pred | not id + pred | not id + not pred |
|---|---|---|---|---|
| passive | 0.00% | 0.00% | 5.43% | 94.57% |
| intervening | 1.85% | 80.16% | 3.58% | 14.41% |

Changing only the identification protocol moves **82.0%** of generators between regimes, so

```
Q = Q(S, 𝒜, Σ, O)      not      Q(S)
```

*Identifiable* and *predictable* are not intrinsic labels on a mechanism. They describe it
relative to an observation and intervention regime and a declared prediction target.

## Running it

Python 3, standard library only, except `mgs_preregistered.py` which needs NumPy.

```sh
python3 mgs.py                # reachability, Boolean clones, integer complexity
python3 mgs_ca.py             # elementary cellular automata
python3 mgs_rewrite.py        # hypergraph rewriting
python3 mgs_ident.py          # identifiability curves and the plateau
python3 mgs_passive.py        # observation without intervention
python3 mgs_stochastic.py     # exploratory Markov family (its criterion failed; kept)
python3 mgs_preregistered.py  # frozen confirmatory run — a few minutes, ~1 GB
python3 audit.py              # re-derives every number in the manuscript
```

## The audit

No number in the manuscript is manually trusted. [`audit.py`](audit.py) re-derives each
quantitative claim from the script that produces it, recomputing from scratch the few that
no committed script emitted, and tags provenance as exhaustive computation, sampled
computation, derivation, or recomputed-in-audit. Current status: **45/45**
([`AUDIT.txt`](AUDIT.txt)), with values emitted to [`RESULTS.txt`](RESULTS.txt) for direct
LaTeX import.

Validation against known results: integer complexity reproduces
[OEIS A005245](https://oeis.org/A005245) for n = 1..32; the six reversible elementary rules
and the 88 equivalence classes under reflection and colour swap both come out right;
`{∧, ∨}` generates exactly the 18 monotone functions of three variables.

## What this does not claim

- **Prior art.** Separating a structural axis from a randomness axis and plotting the pair
  over an enumerated space of systems is the complexity–entropy diagram of computational
  mechanics (Feldman, McTague & Crutchfield, *Chaos* 18:043106, 2008). Our `p_k` is their
  entropy rate. What differs is the horizontal coordinate: `C_μ` is the memory of the
  minimal predictor and a property of the *process*, whereas `i_O(g)` measures the size of
  the equivalence class relative to a declared `(H, O)`.
- **Finite horizon.** Identification in the Markov experiments is equality of the
  finite-dimensional law at `T = 6`, not of the full process law. Every preimage reported
  is an upper bound on the preimage available with unlimited passive data.
- **Structural non-identifiability.** That aggregation creates observational equivalence
  classes is established (Ito, Amari & Kobayashi, *IEEE Trans. Inf. Theory* 38:324, 1992).
  The contribution is the exhaustive quantity, not the phenomenon.
- **Graph isomorphism.** Rewriting comparisons use a Weisfeiler–Leman invariant, sound but
  incomplete, so the 1.09-bit residual is a lower bound on distinguishability.

## Status

Draft. Bibliography verified against primary sources; the numerical claims are audited; the
interpretive claims are not, and the manuscript has not had an external technical reader.
