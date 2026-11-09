# Counterexample continuation, 2026-11-09

AI-assisted research by OpenAI Codex (GPT-6 family), with independent
parallel Codex reviews. No human referee or new counterexample is claimed.

## Starting point and priorities

The checkout was synchronized to `a45ecf5` before work. The latest project
status separates three routes:

1. **Most mature all-ring construction:** the d=13 mobile-cell bond with
   predicted total Virasoro central charge c=14. The finite construction is
   extensively checked; the physical infrared argument retains stated
   hypotheses. Reconcile its finite parameter prescription with the
   independent d=13 red-team report before treating the displayed member
   as a complete certificate.
2. **Highest-value direct test:** C173's nearest-neighbour qubit bond at
   the Fibonacci cusp, with candidate c=2 at d=2. C180 supplies the
   symmetry and field-theory analysis. Existing entropy estimates are
   unconverged. Valley miscibility, number-changing corrections and a
   converged microscopic central-charge measurement remain open.
3. **Stronger long-term statement:** unbounded c at fixed d=4464, whose
   final microscopic-to-infrared step is still conditional.

These are research priorities, not new theorem claims. The original
conjecture concerns the total c of unitary CFT sectors, including sectors
with distinct velocities: [Latorre and Sierra](https://arxiv.org/abs/2403.17242).

## Work in this continuation

- Independent exact few-body scattering for the Fibonacci candidate,
  including its generated attractions and a no-attraction calibration.
- Independent enumeration of its perturbative reduction and
  number-changing terms, with tests at the actual candidate parameters.
- A finite-ring boundary-condition and symmetry diagnostic audit.
- An explicit repaired scale ledger for the d=13 construction, preserving
  the distinction between a finite certificate and its infrared input.

The existing C173/C176/C177 cluster DMRG tasks will not be duplicated.
Each new calculation will include source code, saved results and its
scope. Contributions live in separate `contrib/*_20261109/` folders;
upstream collaborator files are left intact unless a checked correction
is warranted. GitHub synchronization happens at reviewable checkpoints.

## Initial evidence

- [Current project status](../../STATUS.md), REPORT sections 9.137-9.138.
- [C173 candidate](../../codex/work_C173/RESULT.md) and
  [C180 field theory](../../codex/work_C180/RESULT.md).
- [d=13 red-team report](../redteam_d13/README.md).
- Open GitHub issues 1-5 concern the earlier construction's matching,
  defect comparison and infrared estimates. Open PR 6 concerns a
  replacement bound, not the new qubit candidate.

## Checked checkpoint

The continuation pulled C175 (`1a4015f`) during the calculation. Its new
binding result agrees with our independent two-body calculation. The
d=13 route is therefore the strongest mature route; the qubit route
needs a finite-density pairing test or a modified microscopic bond.

- **[proved] Exact binding boundary.** For the symmetric effective model
  with nearest hop u, next-nearest hop +1, V2=-2+v and V3=-g,
  0<u<4, 0<g<2, a bound pair in the opposite-valley total-momentum
  sector exists exactly when v<g u^2/[2(2-g)]. A half-line square
  factorization proves both directions. A 24-component integer vector
  independently certifies binding greater than 1/100 at the candidate.
  See [the theorem and checks](../cusp_binding_certificate_20261109/README.md).
- **[proved algebra; numerical evaluations] Exact pair and scattering
  tools.** An algebraic bound-state energy, normalized wavefunction,
  pair mass, and a probability-flux scattering audit agree with an
  independent banded diagonalization. At u=1 the bound energy is
  -4.5195265346664 versus the free threshold -4.5. See
  [the calculation](../cusp_scattering_20261109/README.md).
- **[proved local coefficient; finite-matrix checks] Microscopic
  reduction correction.** The first omitted term is distance-three
  pair creation, 2 Delta^3/(3V^2), already relative order Delta/V.
  Direct low-energy atomic pairing is off resonance, but a separate
  elimination is required before claiming a relative-order
  (Delta/V)^2 bound on its effect. See
  [the independent audit](../cusp_validation_20261109/README.md).
- **[proved under C92's hypotheses] Stronger d=13 contact and pin
  estimates.** Contact exceeds 1/264, and pin activation is at least
  J/135, replacing 1/4096 and J/4096 under the same legal-space,
  positive-tilt assumptions. See
  [the proof and route ledger](../route_audit_20261109/README.md).

The coordinating agent reviewed the derivations and reran the checks;
two other Codex agents independently reviewed the binding theorem.
No finite-density c=2 conclusion or unconditional d=13 infrared
identification follows. The next high-value analytic task for d=13 is
the dressed contact-clock reduction and common-bulk defect comparison.
