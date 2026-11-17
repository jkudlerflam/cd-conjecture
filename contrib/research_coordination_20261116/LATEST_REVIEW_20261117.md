# Latest review and third coordinated follow-up

Date: 2026-11-17. Author and AI assistance: OpenAI Codex, GPT-6 family.
Status: synthesis of completed contributions and new research assignments.
No independent verification of their full proofs or numerical campaigns is
claimed by this review. The small coefficient and selector-margin identities
below are elementary arithmetic, not new phase results.

This checkpoint follows [the simplicity-first assignments](FOLLOWUP_20261117.md).
All three completed contributions were committed and pushed before this
review. The new prompts were delivered to the same three project chats,
and all three were observed actively working on them.

## What the latest round established

### A simpler qubit lead, with central charge still unknown

The [simpler-qubit campaign](../simple_qubit_candidates_20261117/README.md)
completed 46 runs. Its most attractive simplicity follow-up is the exact
V=4, kappa=a=1 point. A common onsite rotation gives, in Pauli normalization,

```text
H = sum_j [(9/16) X_j Y_(j+1) + (7/16) Y_j X_(j+1)
           + Z_j Z_(j+1)] - (7/4) sum_j Z_j.
```

There are four literal monomials when the uniform onsite field is counted
once, and the largest/smallest nonzero coefficient ratio in this displayed
basis is four. The family still has four independent real couplings before
removing the overall scale. The rotation does not remove the invariant
pairing/DM ratio eight or ZZ/DM ratio sixteen; it puts the small DM term
into the difference of two cross coefficients. It is useful to report both
the shorter expression and the remaining physical hierarchy.

**[Numerical.]** Fresh full-chain L=96 calculations in both parity sectors
passed the registered local, sweep and bond-dimension gates. The final
entropy changes were about 0.0051 and 0.0056. This is an accuracy result
for those finite chains, not a thermodynamic CFT identification. In
particular, the density range was about 0.301, and there is not yet a
calibrated new size sequence or a robust field interval for this point.
Entropy estimates around 2.5--3.1 depend on the size, parity and fit window;
**they are not established values of c**. Its c is unknown.

The V=8 extension reached L=128 and chi=256, but missed the registered
bond-dimension entropy gate (changes about 0.013--0.014, versus 0.01).
Its apparent c estimates vary with window and parity, and the exact c=2
control itself has appreciable finite-size corrections at these sizes.
V=12 also remained unconverged. The older V=40 candidate remains the
stronger large-size c approximately 2 benchmark; its unresolved accuracy
and compute limitations are not repaired by the simpler-point campaign.
No new c>1 counterexample was established in this round.

### Exact exclusions on two faces, with simpler interacting faces open

The [minimal-terms classification](../qubit_minimal_terms_20261117/README.md)
studied

```text
H(a,b,c,h) = sum_j [a(XX-YY) + b(XY-YX) + c ZZ] + h sum_j Z_j.
```

Here c in the Hamiltonian is an exchange coupling, not the central charge.

| Removed structure | Established reduction or result | Limit of the result |
| --- | --- | --- |
| a=0 | Uniform-field XXZ after a DM twist; the established exact-solution infrared classification gives central charge at most one | Includes the conformal, finite-velocity qualification |
| c=0 | Exact free-fermion/BdG analysis gives central charge at most one | Singular endpoints must be treated separately |
| b=0 | Exact staggered-field XXZ reduction, with explicit rigorous gapped regions | The remaining interacting strip is unresolved |
| h=0 | Exact asymmetric Ashkin--Teller/Majorana representation and controlled c=1 neighborhoods | Comparable-coupling and slow/flat-band regions remain unresolved |

Thus pairing and ZZ interaction are necessary for c>1 within this family,
under the stated XXZ classification. Necessity of the DM term or field
has not been established. The zero-field analysis identifies an especially
simple unresolved line, after a common rotation:

```text
H(g) = sum_j [2 X_j Y_(j+1) + g Z_j Z_(j+1)].
```

The rational points g=1 and g=1/2 have coefficient spreads two and four.
They have two literal monomials and one independent ratio after scale.
They are residual tests, not demonstrated counterexamples. Counting two
microscopic Majorana chains does not prove a global interacting c bound.

### A uniform mobile vacuum-clock estimate, not yet unbounded c

The [mobile vacuum-clock contribution](../mobile_vacuum_clock_20261117/README.md)
supplies a **[proved, awaiting independent review]** positive compressed
clock-pin cost in the literal d=22 legal sector. A sufficient allowance is

```text
0 < q <= 1/4,   n f <= 1/64,   0 < J/t <= 1/80,
t/2 <= s <= 3t/4.
```

It applies to the positive-length-tilt ensemble and its minimizing canonical
vacancy blocks, retaining the gapless vacancy dynamics and actual contact
terms. Its J/t condition has neither a charged auxiliary-mass factor D nor
a factor growing with cell length. At q=1/4 the uniform pin cost is at
least 2J/4369. This is a compressed-sector energy comparison, not a full
spectral gap or a dynamical effective-Hamiltonian theorem.

The terminal selector breaks the exact common clock kernel. A separate
local excitation-weight estimate, including sufficiently small ports,
does not supply the missing dressed dynamical matching. Wrong-length
selection, illegal-sector exclusion and the claimed multichannel CFT
remain open. The proposed c=n+1 family is consequently still conditional.

The parameter audit also shows that Delta=4 exactly removes a fine
coefficient adjustment while preserving the stated finite-order mechanism.
Setting D=0 removes a physical term and preserves the selector/clock
arguments, but loses the existing large-D band proof. Neither change
certifies a simpler arbitrary-c Hamiltonian. The architecture still has
many terms and delicate payload scales.

## Follow-ups issued to the existing chats

1. **Counterexamples: diagnose the V=4 bulk phase.** Reuse existing states
   to measure density/bond modulation, connected correlations and structure
   factors. Register a bounded, accuracy-controlled sequence extending
   beyond L=96, in both parity sectors. Distinguish decaying boundary
   oscillations, ordered states and crossover effects before attempting
   another c estimate. Keep the rational bulk coefficients fixed; label
   any boundary probe separately. Use at most two CPU production workers.
   Assigned output: `contrib/v4_bulk_diagnosis_20261117/`.

2. **Bounds: analyze the two-term zero-field line.** Seek a stronger exact
   reduction, a rigorous subfamily gap/bound, or a controlled phase argument
   for 2XY+gZZ, focusing on g=1 and g=1/2. Account for physical boundaries,
   neutral excitations, order and possible emergent small velocities.
   Modest finite-chain checks should distinguish explicit alternatives;
   they should not duplicate the other chat's DMRG campaign.
   Assigned output: `contrib/two_term_qubit_phase_20261117/`.

3. **Fixed-d unboundedness: test mobile length selection.** Begin with
   lambda=f=0 in the exact legal vacuum sector. Compare different cell
   words at fixed total material length with the charge bulk retained.
   Seek an error estimate that keeps the small terminal-overlap scale,
   or a precise failure of the proposed selector chemical potential.
   Do not replace the gapless fluid by an uncontrolled contact average or
   divide by its vanishing gap. Keep a restricted comparison distinct
   from arbitrary wrong-length exclusion and the final CFT.
   Assigned output: `contrib/mobile_length_selection_20261117/`.

The third prompt includes a deliberately weak selector test. Compare
n+2 cells of n bits with n+1 cells of n+1 bits. Both have total occupied
length (n+1)(n+2), but different marker counts. With the frozen selector
choice mu_sel=-(n+2)E_n, their material-energy difference is exactly

```text
(n+1)(mu_sel+E_(n+1)) - (n+2)(mu_sel+E_n)
    = (n+1)E_(n+1).
```

This is the small margin the mobile comparison must resolve. A common
contact expectation would multiply it by J at first order; the prompt
asks the chat to audit that step and its corrections in the actual model.
The easier balanced pair (n-1,n+1) versus (n,n) alone does not address
this marker-changing competitor. No positive mobile margin is inferred
from the frozen arithmetic.

## Status and checks for this review

The completed contributions were read together with their conclusions and
limitations. Their original small checks and datasets remain in the linked
folders; they were not rerun by this coordination-only update. Coefficient
ratios and the marker-count identity above were checked exactly, and local
documentation links and the patch were checked. All three follow-up chats
were observed in progress after accepting the new prompts.

There is still no general upper bound B(d), no verified new c>1 qubit
counterexample from this round, and no established unbounded-c family at
fixed d. The new assignments target distinct gaps while giving priority
to simpler Hamiltonians.
