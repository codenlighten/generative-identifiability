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

Identifiability is `|P(O)| = 1`. The principal finite hypothesis-space results are
enumerated exactly; the few sampled results are explicitly marked.

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

## Papers

- **Paper 1** — [`identifiability.tex`](identifiability.tex) /
  [`identifiability.pdf`](identifiability.pdf). *Identifiability and Predictability Are
  Distinct: Measured Preimages Across Generative Systems.* The measurement framework and
  the two-coordinate thesis, with the pre-registered Markov experiment as the
  load-bearing result.
- **Paper 2** — [`paper2.tex`](paper2.tex) / [`paper2.pdf`](paper2.pdf). *The Information
  Order of Experiments: Intervention Lattices and Generator Identifiability.* The order
  structure on observations: refinement contracts preimages, and informativeness is
  partial rather than total &mdash; a single intervention can be *less* identifying than
  passive observation. Two claims in early drafts did not survive. Submodularity of
  information gain turned out to be a theorem (Shannon's inequality), not the empirical
  finding it was first reported as. And a natural conjecture that valuable interventions
  cross observation blocks is **refuted in general** by a sweep over all 15 partitions of
  the state space; intervention value is governed by the joint information carried by
  reset-specific observable laws, `H(Z_A)`, for which block crossing is only a proxy.

## A worked demonstration

[`discriminator.py`](discriminator.py) shows the framework as a tool rather than a paper.
A device with three internal states and sixteen candidate fault mechanisms; the program
reports what remains unknowable, whether more passive data would help, and which
diagnostic test buys the most information per dollar.

It reports the decomposition rather than a single number:

```
Total mechanism uncertainty:   3.700 bits     log2 of the surviving candidates
Resolved by the passive law:   3.085 bits     H(O)
Residual ambiguity H(G|O):     0.615 bits     what is actually still unknown
```

Four things it deliberately gets right, which are the point of the demonstration:

- **It refuses to rank an experiment when any hypothesis supplies no outcome
  distribution.** No forward model, no information calculation — the contract is enforced
  rather than papered over.
- **It certifies permanence through the joint chain, not the hypothesis-level plateau.**
  A plateau only implies permanence when there is one shared transition operator, and each
  mechanism has its own. A family of models is not a dynamical system until it is lifted
  into the joint state space `(h, x)`, where the mechanism identity becomes a persistent
  hidden coordinate and the transition matrix is block diagonal. The bridge is one line:
  for hypothesis `h` with initial device distribution `μ_h`, the embedding
  `J : v ↦ w`, `w_(h,x) = v_h · μ_h(x)`, satisfies `Σ_h v_h L_h(T) = w B_T(joint)`
  identically, so `K_T(hypothesis) = J⁻¹(K_T(joint))` for every `T` — and once the joint
  kernel reaches its limit, so does the hypothesis-level kernel it determines. The demo
  verifies that embedding on its own kernel basis rather than asserting it.
- **The certificate is a detected plateau, not an exhaustive watch.** The worst-case bound
  here is `T ≤ 39 − 2 + 1 = 38`, but the joint kernel plateaus at `T = 6`. The claim is not
  "we watched to 38" — it is "we found a plateau at 6, and the theorem says a plateau
  cannot later break." That distinction is what makes the certificate cheap.
### Three levels of indistinguishability

Ordinary model discrimination has two states: identified, or not. This has three.

| | |
|---|---|
| **Level 1** | different mechanisms, different observable laws — separable |
| **Level 2** | different mechanisms, *same* observable law — pairwise comparison finds these |
| **Level 3** | different **mixtures** of mechanisms, same observable law |

Level 2 is what pairwise model discrimination catches. Level 3 is what the kernel adds, and
it matters because scientific uncertainty is rarely "definitely A or definitely B" — it is
a distribution over possibilities. If two different distributions over mechanisms imply the
same observable process, the ambiguity is real at the belief level even though every pure
candidate looks unique.

In plain terms: *the system does not only find models that look identical; it finds
different mixtures of models that are observationally identical, even when no pairwise
comparison reveals the ambiguity.*

- **It reports hidden *directions*, not hidden mechanisms, and decomposes them.** The
  certified kernel is split basis-independently into the part any pairwise comparison
  would find and the part it would not: here 4 pairwise-equivalence directions plus one
  higher-order relation, `F10` inseparable from `½F1 + ½F9`. That last line is the
  research result of Paper 2 surfacing as a product feature — no comparison of two
  mechanisms could ever produce it.
- **It separates EXACT from APPROXIMATE mode.** Only exact rational arithmetic over a
  finite class may report CERTIFIED; estimated laws may report expected information gain
  and must never claim permanent impossibility.

We claim no novelty for detecting non-identifiability, which structural identifiability
analysis and Fisher-information methods already do, nor for ranking experiments by
expected information gain, which is Bayesian optimal experimental design. The narrower
combination shown here is: given an explicit finite family of forward models and an
observation regime, determine exactly which distinctions are permanently unobservable,
certify that after a bounded horizon, and rank interventions by information added beyond
what is already known.

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
python3 discriminator.py      # worked demonstration: fault diagnosis
python3 audit.py              # re-derives every number in the manuscript
```

## The audit

No number in the manuscript is manually trusted. [`audit.py`](audit.py) re-derives each
quantitative claim from the script that produces it, recomputing from scratch the few that
no committed script emitted, and tags provenance as exhaustive computation, sampled
computation, derivation, or recomputed-in-audit. Current status: **163/163**
([`AUDIT.txt`](AUDIT.txt)), with values emitted to [`RESULTS.txt`](RESULTS.txt) for direct
LaTeX import.

A full audit re-executes all nineteen experiment scripts, several of which do exact
rational linear algebra and one of which enumerates 1,500,625 generators, so a cold run
from a fresh clone takes roughly half an hour. Outputs are cached per script and keyed on
mtime, so in practice only changed scripts re-run. It is a release check, not a
pre-commit hook.

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

## Research team

Gregory J. Ward, Bryan W. Daugherty, Shawn M. Ryan.

## Citing

Metadata is in [`CITATION.cff`](CITATION.cff); GitHub's *Cite this repository* button reads
it directly. The manuscript is an unpublished draft, so cite it as such:

> Gregory J. Ward, Bryan W. Daugherty and Shawn M. Ryan, "Identifiability and
> Predictability Are Distinct: Measured Preimages Across Generative Systems" (2026),
> unpublished draft.
> https://github.com/codenlighten/generative-identifiability

## License

MIT for the scripts, [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) for the
manuscript and documentation. See [`LICENSE`](LICENSE).

## Status

Draft. Bibliography verified against primary sources; the numerical claims are audited; the
interpretive claims are not, and the manuscript has not had an external technical reader.
