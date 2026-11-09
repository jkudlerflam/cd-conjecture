# Independent review of the cusp binding certificate

Date: 2026-11-09. Reviewer: the route-audit Codex agent, separate from the
authoring root agent. AI-only review; no human review.

**Verdict: [checked].** The proof in
[the binding contribution](../cusp_binding_certificate_20261109/README.md)
establishes the claimed two-body criterion for 0<u<4 and 0<g<2:

\[
 v<g u^2/[2(2-g)].
\]

The following steps were rederived independently:

- The half-line adjacency square has boundary diagonal 1, giving
  `H-E_th = 2(A+u/2)^2 + v P_1 - g P_2` with the stated threshold.
- Completing the boundary square gives the exact coefficient
  `v-g*u*u/(2*(2-g))` of the first amplitude.
- Below the threshold coefficient, the recurrence has two distinct
  unit-modulus roots. Its solution is bounded for fixed u,g, and a
  linear cutoff over R sites gives squared residual O(1/R). This
  supplies a negative quadratic form for sufficiently large R.
- The finite-rank boundary change leaves the essential threshold
  unchanged. Negative spectrum is therefore a normalizable discrete
  state. At v=0, rank-one subtraction bounds its multiplicity by one.
- The robustness allowance controls the two-body perturbation and the
  free threshold separately, including complex changes of hopping
  amplitudes at the fixed momentum. It does not include number-changing
  terms or certify the full microscopic reduction.

The author's standard-library checker was rerun: all 180 rational
identity cases pass, and the finite-support witness has binding margin
`319791/30580778 > 1/100`. Output: `binding_review.log`.

A separate dense matrix rebuild at nine (u,g,v) choices also verifies
the square identity and nonbinding side. The finite-section square
requires adding +2 at the last site to represent the omitted exterior
row of the infinite adjacency operator; this term was included in the
independent check. This numerical corroboration is not needed for the
analytic theorem.

The scope in the author's report is appropriate: this excludes the
unbound dilute-atom premise in the second-order model. It leaves the
many-body phase near density 0.17 and the full qubit chain open.
