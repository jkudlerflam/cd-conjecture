# C131 result

2026-11-07. The compactness reduction is proved. C(2)<infinity remains open.

[Proved] For the explicit even-ring target C_even(2), finiteness is equivalent to an ambient, length-uniform upper bound in a neighborhood of every normalized gauge class. The proof gives |Q_L(h)-Q_L(k)|<=4L^2 M([h]-[k]) and then takes a finite subcover. It also applies to all lengths, or to a fixed bond-dependent admissibility graph, if the same convention is used on both sides. See PROOF.md, Sections 1 and 2.

[Proved qualification] All even rings are an upper-bound envelope, not a universal conformal-vacuum prescription. With C107's ferromagnetic XX energies the doubling Q_L is exactly zero for odd L and tends to pi/2 for even L. The U(1) all-length target of C129 can therefore be pursued without claiming a common raw Casimir limit. The gauge sphere has dimension 11; the further compact rotation quotient has generic dimension 8 and is not asserted to be an eight-sphere.

[Open local lemmas] PROOF.md, Section 4 gives the requested table, including U(1), free, gapped, non-integrable critical, and multicritical loci, the precise local obligations, and their literature tools. The table retains a residual class because no exhaustive classification of the general spin-1/2 gapless locus is known from the supplied sources. Once all ambient lemmas and exhaustion are established, the finite cover follows.

[Scope corrections]

1. C129 addresses the U(1) class itself. C130 is only a draft and its proposed perturbations preserve fermion parity. Neither currently supplies an ambient neighborhood of a general U(1) or free center. Transverse fields, off-axis DM, and the crossover off the critical hypersurface need their own control.
2. A gapped phase needs both a ring-energy convergence theorem with uniform constants and a valid stability theorem. Hastings-Koma clustering alone is not a comparison of different ring sizes. At coexistence, uniform stability above the exact ground space can fail even for the classical ferromagnet with Q_L=0. A bound on the competing phase energies is a better formulation there.
3. Borgs-Kotecky Theorem 4.1 is a classical cubic-torus result. The checked Borgs-Kotecky-Ueltschi Theorem 2.1 assumes spatial dimension >=2. They cannot be cited as a direct theorem for a quantum spin-chain ring. Section 6 states the missing logarithmic partition-function estimate uniform in imaginary time.
4. The supplied C130 draft has no quantitative radius or ambient ball, so a numerical fraction of critical surfaces inside those balls is undefined. No coverage percentage was invented.

[Numerical] NUMERICS.md records the small-ring multistart search, separate odd-length results, and Heisenberg benchmarks through L=8. Local search values are lower bounds on the suprema. The larger-ring multistart search is packaged in submit_c131.slurm and was not submitted because Typhon is inaccessible. See RUNNING.md for reproduction and job details.

[Source discrepancy] The shared target quotes a 2M product-gluing bound, while C98 proves 4M. The reduction uses the proved 4M bound only to absorb finitely many short rings. This does not affect the compactness conclusion.

[Next steps] (i) Complete the relative C129 and C130 estimates, preserving their actual symmetry and uniformity scope. (ii) Attack the uniform gapped ring-energy comparison in a phase with known gap stability, and a phase-envelope bound at classical coexistence. (iii) Run the full multistart array before drawing conclusions about the motion of maximizers with L. No analytical attack on the non-integrable critical or multicritical rows was made here.

Only work_C131 was modified. The local STATE.md is the handoff for the project memory, respecting the task's folder-confinement rule. There was no cluster submission, shared REPORT edit, commit, or push.
