# Fresh adversarial audit of the d=24 construction

Status at this checkpoint: **microscopic ingredients verified; independent
vacancy and infrared audits support the construction at the repository's
theoretical-physics standard; coordinating final review in progress**.
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
