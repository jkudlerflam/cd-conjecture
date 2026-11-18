# V4 locking operators: exact selection rules, conditional stability

## Claim

**The actual symmetries do not forbid the quartet locking term. They permit
an open c=2 region conditionally, but do not establish that V4 lies in it.**
The nearest-neighbour, bulk translation-invariant, Hermitian qubit model is

\[
H=\sum_j(9X_jY_{j+1}/16+7Y_jX_{j+1}/16+Z_jZ_{j+1})
-7\sum_jZ_j/4.
\]

Its microscopic dimension is d=2, without blocking or auxiliary states.
The target is unitary compact Virasoro central charge (or the sum of two
Gaussian contributions with unequal velocities), not c_eff. We assume a
translation- and valley-symmetry-preserving two-component Gaussian regime
only to derive its possible instabilities; its realization at V4 is unproved.

For generic incommensurate valley centres and balanced *dressed* species
densities, all nonderivative bulk vertices in the bosonic local operator
lattice that survive the symmetries are

\[
\cos[2\sqrt2(r\theta_c-s\phi_s)],\qquad(r,s)\in\mathbb Z^2.
\]

The primitive quartet (r,s)=(1,0) changes the candidate emergent charge by
four; (0,1) is relative-density backscattering. Both are allowed. Mixed
vertices cannot generally be dropped. Resonances away from generic
incommensurability are covered by the exact momentum equation in
[LEDGER.md](LEDGER.md), rather than a finite list asserted to be exhaustive.

## Status and main corrections

- **[proved, microscopic algebra]** Translation, parity, and the two
  antiunitary/unitary symmetries and their constant phases are derived for
  this exact bond. Uniform undressed even-separation pair sources are
  forbidden; odd-separation sources and the resonant quartet are not.
- **[proved within the specified Gaussian theory]** With logarithmic
  covariance G, define a=G(theta_c,theta_c), b=G(phi_s,phi_s),
  h=G(theta_c,phi_s). The dimension of the generic vertex is
  **x(r,s)=2(a r²+b s²-2h rs)**. Checking only the pure terms is insufficient
  when h is nonzero. Full stability and the resonance qualification are in
  [GAUSSIAN.md](GAUSSIAN.md).
- **[refuted] C180's covariance formula with drift.** Its oscillator
  symplectic matrix is not the spatial boson commutator metric. A free,
  symmetry-allowed opposite-drift example gives x4=2 exactly but C180 gives
  2.5. The corrected calculation agrees with C191's existing chiral code
  and an independent real-Fourier-oscillator calculation to <3e-14.
  C180's zero-drift formulas and its generic quartet selection survive.
- **[conditional, not a V4 amplitude calculation]** The leading projected
  pair source is 2 Delta³/(3V²), verified including all blockaded backgrounds
  on a ten-site ring. A quartet can be generated at second order in this
  source. Its coefficient depends on the full interacting effective model;
  neither a nonzero order-one coefficient nor a V4 locking length has been
  established. [AMPLITUDES.md](AMPLITUDES.md) separates what is derived from
  power counting and corrects a factor-two arithmetic slip in C179.

The previous finite-size slope 2.5669 remains provisional and is not used
as a scaling dimension. Nothing here certifies c=2 or c>2, excludes a flow
to c=1 or a gap, or establishes an arbitrary-c family at d=2.

## Most informative future measurement

Measure the **connected, momentum-pi inter-valley pair correlation exponent
eta_pair**, in a demonstrated bulk scaling window with bond-dimension and
length convergence. For its leading vertex exp(i sqrt(2) theta_c),
**x_quartet=2 eta_pair**, even with Gaussian mixing and unequal velocities.
Thus eta_pair<1 makes a nonzero quartet relevant; eta_pair>1 makes this
particular perturbation irrelevant. Equality needs the marginal flow.
Use a short physical pair operator such as b_j b_{j+2}, project its pi
component, and verify its overlap and contamination by other harmonics.
This avoids identifying raw density with a conserved charge.

This one measurement targets the main number-breaking threat; it cannot
certify c=2 by itself. Relative backscattering, mixed vertices, resonances,
quadratic stability and a possible symmetry-broken phase still require the
other entries of G and low-energy mode information. Existing short-window,
anchor-dependent pair correlations do not provide a defensible exponent
bound. We have not fitted them to claim relevance or irrelevance.

## How it was established and checked

Author: OpenAI Codex, GPT-6 family. Prior field theory by Claude Opus 5.5
(C179, C180, C189, C191) was read first and is credited in the detailed
notes. No independent human or second-model review is claimed. checks.py
uses one worker/BLAS thread: eight-site symmetry identities, ten-site
virtual paths, finite integer bookkeeping, and 4x4/8x8 Gaussian algebra.
It takes about **0.2 seconds**. No MPS checkpoint was loaded, no eigenstate
production or DMRG was run, and no larger chain was launched. The user's
battery-saving stop remains in force.

## Reproducing

With numpy and threadpoolctl, from this folder:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B checks.py
```

checks.json contains results, runtime, and hashes of the old sources used;
checks.txt contains a compact readable subset. The check imports only the
existing C191 Gaussian functions, not its production scripts. The algebra
checks are not finite-size evidence for the infrared phase.

## Relation to existing work

This specializes and corrects the field-theory ledger following
[the V4 interpretation audit](../v4_interpretation_audit_20261117/README.md).
Old code, outputs and contributions remain unchanged. Recommended source-
audit note: C180's drift covariance and uncoupled stability inequalities
must not be used at nonzero drift; C191 already uses the proper chiral
metric. Quartet selection is robust, while its V4 amplitude and dimension
remain unresolved. The coordinator owns shared source/status edits.
