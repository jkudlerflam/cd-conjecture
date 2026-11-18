# Results of the four chats on large qubit c and fixed small dimension

Date: 2026-11-18. Author and AI assistance: OpenAI Codex, GPT-6 family.
Status: review and synthesis of the four completed GPT-6 Astra Max
contributions. This review read the source arguments and scope statements;
it did not rerun their computations or provide independent human review.
All four chats were observed idle with their assigned turns completed.

The largest advance is an **arbitrary-n leading quadratic band construction
at physical d=9**. The corrected d=10 mobile model also gains exact finite
sector certificates and a useful charge-field modification. The two qubit
tasks identify stronger interaction obstructions, without establishing a
new microscopic central charge. No full fixed-d unbounded-c CFT theorem
follows from this round.

## Modified d=9 family: the all-n band obstruction is resolved

The [markerless all-n contribution](../markerless_all_n_20261118/README.md)
proves that every rational moment coefficient L_j in the earlier loop
criterion is strictly negative. After reversing the polynomial variable,
orthogonal-polynomial interlacing and a signed resolvent identity give
opposite strict signs for the two factors in each summand. Every summand
is negative, so cancellation cannot occur. This resolves the exact
nonvanishing question left open by the eleven-state analysis.

**[Proved in the contribution; internally checked.]** For every finite n,
an explicitly modified nearest-neighbour nine-state bond has a leading
quadratic logical symbol with exactly 2n simple Fermi points and distinct
sums for all different unordered pairs, including repeated-point pairs.
Its physical coefficients depend on n but stay fixed as chain length grows.
No new local state is added.

The proof uses a boundary-biased vacuum clock, an independently scaled
active rail and a physical L1-only field. It includes an exact invariant
grammar, a clock gap at least 1/2 and zero-port length selection, including
head-number competitors. It is not an all-n proof of the original
unbiased, uniform-charge-field certificate schedule, whose certified
members remain exactly n=9,10,11,12,14,16,20.

Finite-cluster quadratic error allowances can be chosen after the band
margins. Nonzero-port thermodynamic ground-sector selection,
volume-independent elimination and interacting CFT stability remain open.
The new field is not the conserved charge. In particular, the mobile
charge-field result below cannot simply be applied to this readout.

Simplicity has not improved: the proof introduces more independently
controlled scales and retains severe coefficient hierarchies. The stated
counts are 53 real Hermitian basis elements, or 90 literal matrix-unit
monomials in a different explicitly stated basis. Neither count is a new
physical dimension. The result is an analytic existence construction,
not a moderate-coefficient realization.

## Original d=10 mobile family: finite sector control and an exact field repair

The [mobile contribution](../markerless_mobile_control_20261118/README.md)
studies the corrected contact-weighted model at n=9. On a canonical ring
with 20 material sites and two vacancies, J/t=1/25, exact integer
certificates show that the intended length word (10,10) beats (9,11),
(20), and (6,7,7). The smallest certified gap exceeds 1.20e-4 t.
The comparisons persist across a stated selector interval and a small
J/t interval. They include both head-number directions, but are still
three specified finite competitors, not all words in the thermodynamic
limit. The occupied density here is 10/11, not half occupancy.

**[Proved finite extension.]** A charged-head domination inequality and
an explicit finite port/field allowance extend these comparisons. With
f_ext>=0 and sqrt(2)|lambda|+f_ext/2<=1/1000, the smallest remaining certified gap
is above 4.00e-5 t. This whole-word norm estimate is confined to this
finite ring; it is not a uniform defect bound with arbitrary spectators.

**[Exact algebraic improvement.]** Moving only the conserved charge field
outside the material-contact embedding makes it proportional to the exact
global charge. This removes the position-dependent bit-field error at
every J/t without adding states or off-diagonal processes. To match the
leading frozen band at occupied density rho, the normalization is
f_ext=alpha_infinity(rho) f9. Other material energies must remain contact
weighted; moving all of them outside would repeat the original failure.

Second-order vacancy corrections are unfavorable for some balanced and
shorter-cell competitors, so the earlier d=22 favorable-sign argument does
not transfer. A uniform relative energy-response bound at fixed finite
density, especially half occupancy, is still needed. The selector budget
and the dressed band's 2^-35 normalized allowance are different; there is
no certified J/t schedule for the latter. The original n=9 coefficient
spread remains above 10^407 at half occupancy in the stated basis.

## Qubit two-mode route: actual quartet mixing and failure of a leading cancellation

The [two-mode contribution](../qubit_stable_two_modes_20261118/README.md)
evaluates complete interacting intermediate-state sums for finite resonant
N to N+4 transitions in the literal four-term qubit family. A nonzero
amplitude scales as Delta^4/V^3, and microscopic finite-ring spectra
approach the predicted avoided crossing. The calculation keeps the
generated attractions and blockade; it does not replace the intermediate
theory with free particles.

**[Exact finite perturbation result.]** A special twelve-site state is
dark in the truncated leading model, but its cancellation is lifted in
the full microscopic Hamiltonian. The next off-diagonal coupling is

```text
(6 sqrt(30)/5) Delta^6/V^5 + O(Delta^8/V^7).
```

These results establish finite mixing and expose a false protection
mechanism. They do not determine the finite-density continuum quartet
coefficient or its scaling dimension at V4. The finite crossing fields
used for matching are measurement devices, not size-dependent proposed
thermodynamic Hamiltonians.

Within the original four-term family the balanced quartet momentum stays
resonant. Adding v XX detunes its bare momentum but generates another
pair source. A rational five-term example has coefficient spread 4064;
it is not an established stable c=2 region. The chat also calibrated a
connected momentum-pi pair estimator and a complete integer-lattice
stability test that detects mixed vertices missed by pure-term tests.
No inadequate V4 short-window data were fitted and no production restarted.

The next decisive physics input remains finite-density matching, at one
fixed Hamiltonian, of the quartet amplitude, pair exponent and full
Gaussian covariance. A nonzero finite-ring amplitude alone does not prove
that the thermodynamic mode gaps.

## Qubit route beyond two modes: a crystal preempts the tested third valley

The [three-mode contribution](../qubit_three_mode_deformations_20261118/README.md)
exhausts the physical normal deformations near the m=3 parent under its
specified lower-order kernel constraints. There are five real directions
normal to the full parent manifold and seven leading controls after the
complex resonance displacement is included. No independent repulsion-only
direction is hidden among the remaining bond entries at this order.
Particle-number-changing companions are retained.

**[Proved all-length inequality for the leading model, exact arithmetic.]**
For the number-conserving slice, with
gamma=|Re zeta|+|Im zeta| in (0,1/10], mu>=0.7 gamma selects exactly the
three period-three crystals on compatible rings L>=12. For strict
inequality the gap is at least mu-0.7 gamma in the stated leading energy
units. The highest bare valley minimum exceeds 0.71 gamma. Thus increasing
the chemical potential to fill that third bare valley reaches a gapped
crystal first. This is a local positivity certificate extended to every
compatible ring, not a finite-size phase fit.

A six-Pauli-term number-changing escape was also tested on bounded rings.
It retains overwhelming support on a single residue class and does not
establish three liquids. Its coefficient ratio can be as low as 4.61 at
epsilon=1/5, but the perturbative phase identification is uncontrolled
there. The theorem is restricted to the stated leading model and slice;
it is not a universal c<=2 bound or a uniform microscopic exclusion.
Mixed-residue interface energy and motion are the next discriminating
calculation for the number-changing escape.

## Verification and publication

The original contributions report passing bounded exact/numerical checks:
literal-bond and moment identities for the all-n construction; integer
energy witnesses for the mobile sectors; complete spectral sums and
rational cancellation lifting for the two-mode qubit study; and rational
positive-block certificates plus full-sector elimination for three modes.
No new microscopic central charge was measured in this round.

All four contributions are published:

- d76e315a: interacting qubit quartet amplitudes and cancellation lifting.
- b8181dda: complete three-mode normal map and crystal preemption.
- e57633f6: all-n markerless moment theorem and band construction.
- 35bb8752: finite mobile sector certificates and external charge field.

The collaborator also published [informal sharing notes](../collaborator_notes/README.md).
Their charge-locking length and phase-1.53 protection discussion remains
an estimate with an uncomputed thermodynamic coefficient. The new finite
quartet matching does not determine that coefficient. No direct assembly
of the new all-n clock/readout with the original mobile estimates has
been established.
