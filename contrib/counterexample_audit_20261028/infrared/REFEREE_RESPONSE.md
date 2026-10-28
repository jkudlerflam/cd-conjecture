# Response to the independent d=24 referee

Status: **[assessment]** The referee's central distinction is justified.
"Verified at the physics standard" is too strong if it means that this audit
has established the microscopic d=24 counterexample. A defensible description
is: **a credible analytic counterexample candidate, with certified finite
ingredients and a substantially checked conditional infrared argument**.
The narrower statement that the vacancy mechanism survives the identified
infrared attacks is supported. No new fatal error has been found.

Read independently: `contrib/referee_d24/README.md` as supplied in commit
8c1ab8f. Response by the same GPT-6-family worker that wrote this folder's
infrared audit, 2026-10-28. The original report's opening and verdict have
been corrected to make the matching hypothesis explicit.

## What my new calculations actually add

The report goes beyond merely repeating "if matching works, c=24":

- It derives the charge susceptibility and checks its normalization.
- It verifies that the Euclidean smooth-density kernel has no inverse slow
  velocity, and explicitly distinguishes it from the full contact response.
- It solves the dangerous slow/fast Gaussian determinant. The raw mixed term
  is O(J), but its relative stiffness vanishes when J/t is chosen last.
- It controls a connected **second-order** quantum-interface integral, with
  a finite first spatial moment for payload dimensions above 1/2.
- It identifies a finite length cutoff for near-competing cells and an
  extensive bound for long cells that avoids using a uniform internal gap.
- It uses the exact port norm to bound certain boundary changes by O(J lambda).
- Independent ED checks the phase fluid's contact and density responses.

These eliminate concrete proposed singularities and strengthen the physical
case for matching. They do not constitute a finite-J construction of the
coupled microscopic low-frequency action in every competing sector.

## The exact remaining hypothesis

There is no newly discovered allowed relevant vertex. The outstanding issue
is **uniform matching**, which my initial verdict treated too readily as
automatic once the leading calculations behaved well.

At fixed h23, one needs a local matching scale and remainder satisfying both:

1. The normalized bulk stiffness, Fermi-point displacements and all relevant
   matching couplings approach the specified 23-band-plus-charge reference
   as J/t tends to zero, while retaining all gapless oscillators.
2. For bounded-size material defects and interfaces, after common bulk
   contributions are subtracted, the correction to the energy difference
   between the selected sector and a sector with k such defects is bounded by

       |delta E_defect(J,L)| <= J k epsilon_h(J/t) + boundary rounding terms,
       epsilon_h(s) -> 0,

   with constants uniform in L, defect separation and allowed defect pattern.
   Fixed port changes may instead be budgeted explicitly by O(J lambda),
   below the selector margin.

**Same-day correction from the [long-cell attack](../matching_attack/README.md):**
the original wording allowed an arbitrarily long wrong cell to count as one
defect in this absolute bound. That version is false: removal of the selected
bulk's charge-relaxation energy gives a positive correction proportional to
J-squared times the replaced length. This favors the selected sector and does
not refute the one-sided exclusion inequality. For unbounded wrong cells,
use a length-weighted remainder controlled by their extensive bare cost, or
retain the favorable contribution in a one-sided comparison. Neither repair
has yet been established for arbitrary arrangements and interfaces.

The second-order calculation supports the second assertion around the chosen
reference. It does not establish it uniformly for the exact interacting ground
states in all defect sectors. In particular, the connected payload correlation
estimate already uses a Gaussian-phase hypothesis; a bounded Fourier response
gives an extensive bound, while a defect comparison needs the common bulk
subtraction and interface control as well. The exact port norm does not by
itself provide that subtraction for the charge/walker vacuum.

Higher smooth-charge cumulants vanish in the Gaussian reference, and the
remaining operators are irrelevant. These are sound reasons to expect a
regular Wilsonian completion. They are not a separately executed matching
calculation for the full h23 source. An adequate physical derivation need not
be a mathematical CFT theorem or provide every numerical constant, but it
must explain this uniform completion rather than invoke it as the desired
conclusion. The strongest precise statement currently supplied by this audit
is therefore **P and M imply c=24**, together with nontrivial checks supporting
M. P denotes the frozen c=23 phase; M includes the two uniformity statements.

This distinction matters independently of whether the final claim is true.
Failure to complete this argument is not evidence for a different phase.
Conversely, not finding an instability is not verification that every matching
remainder has the required bound.

## What is and is not a logical objection in the referee report

Lack of a measured central charge is not by itself a logical defect in an
analytic existence argument. Nor is the use of "sufficiently small fixed
couplings": once a nonempty admissible region is established, it defines
counterexamples even without choosing decimal values for one point. Here the
uncomputed matching constants limit how concretely that region has been
established and how the claim can be reproduced, rather than proving that it
is empty.

The extreme scale hierarchy is also not a contradiction. It makes numerical
validation of n=23 effectively unavailable and places the burden on the
analytic transfer. A demand for human or different-model review is an
evidentiary acceptance criterion, not an additional Hamiltonian hypothesis.
It is reasonable, but neither authorship nor agreement changes the logical
content of a derivation.

The referee's cell-period imprint is a real effect. It is included in the
need to bound corrections relative to the actual period-24 payload velocities
and momentum margins. Having the same period is necessary to avoid the toy's
new dimerization mechanism; it does not by itself bound the size of the
correction. This does not add a missing vertex to the current operator audit.

## Can a small-n diagnostic settle the issue?

**No, not logically for n=23.** A calibrated c=n+1 measurement for the same
microscopic bond family at n=1 or 2 would be valuable. It could expose an
incorrect embedding, an unaccounted cell-period perturbation, a bad defect
selection argument, or a mismatch between the proposed action and the
microscopic model. Jointly measuring defect energies and the J-dependence of
the effective couplings would test M more directly than c alone.

A positive result at n=1 or 2 does not establish the frozen n=23 phase, its
extreme velocity hierarchy, or a matching estimate uniform over its defect
sectors. The construction explicitly does not assume uniformity in n, and
the certificates and smallness budgets change with n. A small-n failure at
moderate J is likewise not a refutation of the sufficiently-small-J claim
until the required hierarchy and length crossover are demonstrably resolved.

Thus the referee's proposed diagnostics would improve confidence and could
justify a physicist's practical acceptance, but they would not logically
remove P and M. The analytic path to a stronger conclusion is to complete
the uniform local matching/sector comparison; the numerical path supplies
independent evidence about that mechanism. They are complementary tests.

## Recommended wording

"The finite bond, material representation and nonresonant quadratic reference
have been checked. The vacancy lift has no identified infrared obstruction;
its susceptibility, Gaussian stability and leading defect dressing agree
with the proposed 24-mode theory. The d=24 counterexample remains conditional
on the frozen interacting phase and uniform finite-J matching, so we classify
it as a credible analytic candidate rather than a verified microscopic
counterexample."
