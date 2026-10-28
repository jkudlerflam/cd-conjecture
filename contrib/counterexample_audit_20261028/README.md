# Fresh adversarial audit of the d=24 construction

Final verdict: **verified as a controlled counterexample argument at the
repository's stated theoretical-physics standard**. The microscopic
ingredients, vertex classification, and scale hierarchy pass independent
audits. No refutation of the model was found. The finite-J conclusion uses
the explicitly stated conventional local matching and Gaussian RG input.
No rigorous lattice-to-CFT theorem or microscopic c=24 measurement is
claimed. Prior AI endorsements were not used as evidence.

The user requested an independent audit through a justified verification
or refutation verdict. Four workers were used concurrently: a coordinator
and three fresh subagents, followed by complementary cross-review and
new quantum-impurity diagnostics. All work is dated 2026-10-28 and is
AI-assisted (OpenAI Codex GPT-6 family); no human referee is claimed.

## Independent results

* [Microscopic audit](microscopic/README.md): literal full-Fock two-cell
  calculations at n=1,2 independently recover the effective hopping,
  onsite terms, phases and normalization; complete 512-bit n=23 matrices
  were regenerated exactly; all 46 simple roots and quadratic/quartic
  margins recertified. Frozen-sector bounds and the primary-source local
  Schrieffer-Wolff application pass.
* [Vacancy audit](vacancy/REPORT.md): independently constructed physical
  and squeezed matrices agree, including complex oriented bonds, seams,
  parity and degenerate charge seas. Number and flux curvatures recover
  K_gap=4. A separate scan checks 1,167,480 ordered nonforward quartics
  and 56,039,040 potentially dangerous marker branches.
* [Infrared audit](infrared/README.md): derives the dynamic Euclidean
  response, quantum interface decay, scale ordering, and the distinct
  treatment of arbitrarily long competing cells. Independent phase-model
  response calculations reach 92,378 states.
* [General Gaussian bound](GAUSSIAN_BOUND.md): an operator-geometric-mean
  proof covers arbitrary current mixing and unequal chiral velocities.
  A relative stiffness error eta<5/13 keeps every reference-dimension-three
  non-Gaussian vertex above dimension two.
* [Finite-J cross-review](vacancy/FINITE_J_ADDENDUM.md): audits possible
  inverse-J constants and distinguishes charge-stiffness corrections
  from induced internal interactions.
* [Quantum-impurity diagnostic](microscopic/quantum_impurity/README.md):
  noncommuting one- and two-impurity sources with machinery gaps J Delta
  retain finite charge dressing. The pair test resolves an additional
  J^3 correction caused by induced quantum pair flips. This supports the
  mechanism but is not a simulation of the complete 24-state chain.
* [Loop-scale check](microscopic/LOOP_SCALE.md): the tiny loop's visible
  momentum shifts are consistent with the exceptionally small slopes at
  every Fermi crossing. The onsite expectation divided by the actual
  slope reproduces the displacement, and the certified centers remain
  separated. The largest hopping entry is not the relevant denominator.

## Meaning of the verdict

The reviewed statement is an existence construction at fixed n=23:
choose g sufficiently small, then nonzero lambda sufficiently small,
then nonzero J/t sufficiently small, all independent of the eventual
ring length. The resulting symmetry-preserving weak-coupling phase has
24 right and 24 left compact Gaussian modes, hence c_L=c_R=24, conventionally
reported as c=24. The repeated microscopic bond has d=24, so this exceeds
d-1=23. It is not a claim about J=0 or about a size-dependent sequence
of Hamiltonians.

The certificate fixes the internal geometry and proves a finite-g
reference. It does not provide a certified numerical interval for the
final lambda and J. The saved moderate-scale bond fixture is a structural
check, not a numerically demonstrated violating parameter point. The
weak-coupling existence inference is the same type of bosonization/RG
reasoning permitted by this repository's declared standard.

This audit verifies more than a count of free modes. The argument keeps
the gapless charge field, checks the whole low-dimension local vertex
lattice, compares corrections with the slowest retained velocity, and
addresses defect/interface energies separately from extensive bulk errors.
The charged and internal U(1) symmetries remain exact. The stable region
is asserted inside that symmetry family.

The stronger statement still absent is a uniform microscopic all-order
linked-sector inequality and a rigorous lattice-to-CFT limit for the full
24-state Hamiltonian. Neither a finite-band certificate nor our finite
numerics is that theorem. The verdict must retain this limitation; it
does not constitute an unconditional mathematical disproof or a direct
microscopic measurement of c=24.

## Clarifications uncovered

The gap coordinates require a root coordinate and its cyclic quotient
on a finite physical ring. The Gaussian response bound is for Euclidean
smooth density, not the complete lattice contact susceptibility or a
pointwise real-time response. Interacting oscillatory bilinears need not
have dimension exactly one; the weaker bound Delta_min>1/2 still makes
the quantum interface kernel's first moment finite. Pure charge stiffness
can change by O(J), which is small relative to t, while induced internal
energy corrections are O(J^2/t). None of these corrections currently
provides a refutation.

The continuum-to-lattice matching and the uniform all-order microscopic
per-defect inequality must not be confused. The reports assess the former
at the stated controlled weak-coupling physics standard; they do not claim
to have proved the latter. See each report for exact hypotheses and
reproduction commands. The original construction is
[here](../counterexample_campaign/COUNTEREXAMPLE_D24.md).
