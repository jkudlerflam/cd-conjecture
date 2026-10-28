# Referee report on the species-graph counterexample to c ≤ d − 1

## Claim

The proposal under review ([PROPOSED_COUNTEREXAMPLE](../research_continuation/PROPOSED_COUNTEREXAMPLE.md), with [IR_AUDIT](../research_continuation/IR_AUDIT.md), [SCALE_HIERARCHY](../research_continuation/SCALE_HIERARCHY.md), [MICROSCOPIC_HAMILTONIAN](../research_continuation/MICROSCOPIC_HAMILTONIAN.md) and the folders they cite) asserts that for every large prime p a one-site translation-invariant, nearest-neighbour, Hermitian chain with d ≤ 192C(p²+p+1)[1+log₂(2(p²+p+1))]+48, C independent of p, flows at one fixed nonzero coupling to N+1 compact bosons, N = (p+1)(p²+p+1), so c > d − 1. The logical flavours are the edges of the PG(2,p) incidence graph, whose girth six forbids neutral flavour-changing monomials with fewer than six fields.

My conclusion is that a careful theoretical physicist should not yet accept this as a refutation of c ≤ d − 1. I found no fatal error and no relevant or marginal operator that the infrared analysis misses, but two steps that carry the conclusion are argued in prose rather than established, and no microscopic instance of the mechanism has been tested.

## Status

Referee report, with independent re-implementations of the decisive finite checks and a small calibrated numerical test.

## How it was established and checked

- Author: Claude (Opus 5.5) referee agent, 2026-10-28, for Jonah Kudler-Flam. No Codex/GPT was used for the checks or the report.
- I formed the verdict below before reading the parallel Astra review ([codex/work_C24](../../codex/work_C24/RESULT.md)). It reaches the same ranking of the two serious gaps independently.
- Every script here was written from the statements in the notes, without importing the repository's own check scripts. The FES driver reuses the VUMPS solver in `../d2_followup_220_1282/vumps.py`, unmodified.
- No human has refereed this report.

## Verdict and ranked problems

The finite algebra is right, and everything I re-implemented passes. The infrared operator list in IR_AUDIT §2 to §3 is complete given exact species and marker conservation. Every resonant non-Gaussian vertex has free dimension at least three except single and paired backscatterings, whose mismatches are π/[6(N+1)] and π/6 as stated. The order of limits (fix p, the core, the wrapper, then a small coupling s, then L → ∞) is the right one, and the scale hierarchy closes formally because the exact core has no splitting that s must beat. What is missing is the existence of the core in the claimed form and a controlled reduction of the full mobile chain.

### Fatal

None found.

### Serious

1. The initialized core is a composition argument, not a construction. The 2^N-dimensional code, its positive gap, the READY readout contract and, above all, the p-independence of C rest on [EXISTENTIAL_CORE](../torus_bottleneck/EXISTENTIAL_CORE.md) and [COMPILER_INTERFACE](../torus_bottleneck/COMPILER_INTERFACE.md). What [controller_instantiation](../controller_instantiation/README.md) actually supplies is a host-side circuit generator, a preloaded dispatcher that cannot write its own program, and a standalone length counter that uncomputes its result before any handoff. Generating a verified instruction tape from blank, synchronizing the endpoint shuttle with the clock (conditional in COMPILER_INTERFACE §4), and READY correctness at every history vertex are not built anywhere. This layer has a demonstrated error rate. The original format verifier accepted both values of a neutral bit (I reproduce this below). It was repaired, but nothing guarantees it was the last such alias. A universal-machine argument of Gottesman–Irani type makes a p-independent C plausible, and I found no resource obstruction, but it is not shown.

2. The thermodynamic step is asserted at its decisive point. Passing from the first-order free model to c = N+1 requires that the excluded machinery (history clocks whose only stated gap bound is of order 8^(−n)/n⁴, wrong-length cells, vacancies, helices, selector heads) can be integrated out into a quasi-local action whose corrections are o(s). IR_AUDIT §9 gives a correct cellwise defect-energy inequality (8), but the connected-excursion bounds (9)–(10), with marker projectors that move during an excursion and retained propagation left unexpanded, are stated rather than derived. I agree with C24 that they do not follow from the defect count alone. The project's standard does not require a theorem here, but it does require a controlled physical derivation or one tested instance, and neither exists. Not even the d = 240 one-mode fixture in `controller_instantiation/local` has been simulated.

### Fixable

3. The RG tolerance should be stated with its N-dependence. There are of order N² marginal inter-flavour backscatterings with mismatch down to π/[6(N+1)], and generically their one-loop flow is enhanced by N in collective channels, so the requirement is N g log(1/δ_pay) ≪ 1 rather than g log ≪ 1, with dressed k_a shifts below δ_pay/8 ≈ p^(−3)/10. Harmless at fixed p, but IR_AUDIT quotes only a generic tolerance.

4. Six-field conversions are irrelevant only through the full metric bound, not through positivity. With the natural edge order, 9 of the 28 hexagons at p = 2 (110 of 234 at p = 3) carry an exactly non-oscillating conversion, of free dimension three. A one-third softening of the dimension form along a cycle direction makes one marginal. IR_AUDIT does impose κ > 2/3, so this is a matter of stating it as the operative condition.

5. Bookkeeping. The visibility quoted as of order p^(−5) omits the unevaluated readout fractions p_L p_R and compensator norms of order b², so B_p ∼ b⁴p⁵/(p_L p_R). The length density is √2/12 in torus_bottleneck §7 but 1/2 in IR_AUDIT. Both meet the margins, but ν = 1/2 makes the h = 2 length umklapp commensurate, harmless only while K_ℓ > 1/2.

6. Interpretation. The infrared theory has N+1 different velocities and, at incommensurate fillings, ground states in size-dependent zero-mode sectors (C13 §7). Calling this c = N+1 is standard (it is how one assigns c = 1+1 to the Hubbard chain), but the final claim should state that this is the reading of "CFT(s)" and "total c" being used.

### Cosmetic

7. C is not evaluated, and the first violation is astronomically large. Counting only components that are already explicit (circuit record D_exact ≈ 6q log₂q, the dispatcher factor 26,488, the fifteen-letter clock, three static marks, and the FAIL and READY qubits) gives C ≥ 2.2 × 10⁷, and the first prime with N+1 > d − 1 then satisfies p ≳ 4 × 10¹¹, at d ≈ 8 × 10³⁴. The program generator, length counter, verifier families and shuttle are not in that count. A further factor of 10² to 10⁶ for them moves the first violation to d ∼ 10⁴¹ to 10⁵³. The core length grows as a fixed power of p (the explicit circuit already uses 8 × 10⁴ sites and 1.2 × 10⁹ local gates at p = 2), so the coupling s must lie below exp(−O(p^K)). None of this affects existence, but no violating chain can ever be written down or simulated, so acceptance must rest on the argument alone.

8. Several composition claims are tagged "[proved]" where "argued" is meant, and STATUS's "c > d − 1 for sufficiently large p" should give the scale of p. The claim files changed during this review (MICROSCOPIC_HAMILTONIAN and controller_instantiation appeared).

## What I re-implemented

`graph_checks.py` builds PG(2,p) from F_p³ and enumerates every multiset of fermion fields, counting creation and annihilation separately. At p = 2 (q = 14, N = 21, kernel rank 8) it covers all 12.3 million multisets of up to six fields, and at p = 3 (q = 26, N = 52, kernel rank 27) all 5.36 million of up to four. No neutral monomial with at most four fields changes a flavour. At p = 2 there are exactly 56 neutral six-field conversions, two orientations of each of the 28 hexagons, as the circulation argument predicts.

`isometry_checks.py` builds the full two-port isometry W for the six-cycle (186,624 output states), a graph with a four-cycle and an isolated vertex, and a path. Isometry, both compressions εf† and ε(1−1/2N)f†P, vanishing non-edge compressions and the channel identities hold to 10⁻¹⁵ (the dressed identity also on the Fano graph, on 2²¹ states), and no nonzero entry of V_L violates species covariance. The basis transducer of CIRCUIT §1, run on actual species words, conserves every species multiset, is injective, and reproduces V_L exactly on all 64 × 9 six-cycle branches, and it conserves species injectively on 8,400 sampled Fano branches. The load, exchange and query gates are charge-preserving involutions of the full two-site alphabet, the three-gate controlled swap is correct on every three-site word at (q, m) = (3, 2), (5, 3), and the prefix amplitude preparation is exact for all M ≤ 256.

`input_selection.py` reproduces the format alias (the original format accepts 8 words of a four-slot register instead of 4, the repaired blank format exactly 4). A small Feynman–Kitaev history (check, copy to FAIL, encode, readout, decode, uncheck) has kernel dimension 4 = 2^N with an all-time FAIL penalty and 9 with a terminal-only one, since FAIL is uncomputed. The gap obeys the acceptance lemma, and the readout compresses to (|S|/C)V†FV to 3 × 10⁻¹⁷.

`hierarchy.py` tabulates ε, η, δ_pay and the component lower bound on d against p.

## Numerical test of a reduced logical chain

The interacting p = 2 logical chain has 2²² states per cell, so I test the reduced logical chain the infrared argument relies on, a length fluid at ν = 1/2 plus n_pay payload flavours with the proposal's dispersion, k_F ∈ (π/4, π/3) and, for two flavours, the p = 2 mismatch 2|k₁ − k₂| = 0.024. I add every two-site term IR_AUDIT and C13 list as generated at finite coupling (length-dependent payload hopping, payload-dependent and correlated marker hops, all on-site and nearest-neighbour densities), with random signs and strength g up to half the hopping. `fes_effective.py` runs VUMPS finite-entanglement scaling, and `analyze_fes.py` gives c from 6 dS/d ln ξ and from dS/d ln χ, both calibrated on the free point.

| Flavours | g | χ range | ξ_max (cells) | c from 6 dS/d ln ξ | c from dS/d ln χ | exact |
| --: | --: | :-- | --: | --: | --: | --: |
| 1 | 0 | 32–96 | 29 | 1.96 | 1.99 | 2 |
| 1 | 0.2 | 48–128 | 37 | 2.07 | 1.88 | 2 |
| 1 | 0.5 | 32–96 | 30 | 1.98 | 1.97 | 2 |
| 2 | 0 | 32–96 | 10 | 3.69 | 3.04 | 3 |
| 2 | 0.2 | 24–64 | 8 | 3.33 | 2.94 | 3 |

With one payload flavour both estimators agree with the free calibration to within 0.1 at couplings up to half the hopping, with no sign of a gap. With two flavours the ξ estimator is not converged even at the free point (it overshoots by 0.7), and the χ estimator agrees with its calibration to within 0.1. These ran locally (the Typhon array was cancelled before it reached larger χ). The two-flavour runs reach only ξ ≈ 10 cells, well short of the 1/δ ≈ 40 cells at which the small-mismatch backscattering is cut off, so they cannot address item 3.

This confirms the textbook ingredient, that a multicomponent Luttinger liquid with nonresonant margins survives generic symmetric interactions, and nothing more. It does not test the encoded core or the elimination of its machinery.

## What would make it acceptable

At the project's physics standard I would accept the construction as a refutation once (i) a single initialized autonomous controller is specified as a finite transition table or a symbolic compiler with an invariant proof, composed with the verifier, the shuttle and the clock, tested on all inputs of small analogues, and counted, giving an explicit C, and (ii) the local elimination is either derived for the moving projectors or demonstrated on one complete microscopic instance, for example the d = 240 fixture or a two-mode analogue with a history-type core, showing c = n+1 at a fixed small s with two calibrated estimators and second-order couplings computed from the actual cell. The final claim should also state the tolerance of item 3, the metric condition of item 4 and the convention of item 6. An optimized d, a simulation at a violating prime, or a rigorous lattice-to-CFT theorem should not be required.

## Reproducing

Run from this folder with numpy and scipy:

```
python3 graph_checks.py        # ~30 s, writes graph_results.json
python3 isometry_checks.py     # ~1 min, writes isometry_results.json
python3 input_selection.py     # ~2 min, writes input_selection_results.json
python3 hierarchy.py           # instant, writes hierarchy_results.json
python3 fes_effective.py NAME --npay 1 --g 0.2 --seed 1 --kF 0.28 --chis 16 24 32 48 64 96 128
python3 analyze_fes.py 'data/*.json'
```

`data/` holds the FES outputs. `tasks.txt` and `fes.slurm` are the Typhon array (job 17514212, cancelled early).

## Relation to the rest of the repository

This reviews item 1 of "Where help is most useful" in [STATUS](../../STATUS.md). It supports keeping the species-graph route as the leading candidate, labelled a conditional construction rather than a refutation. If completed, it would refute c ≤ d − 1 while remaining consistent with the repository's view that c is bounded by an onset length times log d, since here the onset length grows as a power of p.
