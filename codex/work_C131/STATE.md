# C131 handoff

2026-11-07. Work confined to this directory. This is a local project-memory update because the run prohibits writes outside work_C131.

[Proved] Compactness equivalence for C_even(2), C_all(2), or a fixed admissibility graph. Lipschitz estimate 4 L^2 in the quotient norm. A putative divergence accumulates at a bond lacking a uniform local upper bound, with lengths tending to infinity. Finite small lengths can be absorbed using the proved C98 bound Q_L<=4L.

[Conventions] E_N is the full periodic spin ground energy and H_2 includes both orientations. General target here admits all even L, an envelope containing any compatible even vacuum subsequence. C129 pursues all L on U(1). Evenness alone is not a universal physical-vacuum prescription. The gauge sphere is dimension 11, the further rotation quotient has generic dimension 8. It is not an eight-sphere.

[Proved corrections] C107's ferromagnetic XX example has doubling Q_L=0 exactly on odd L and Q_L->pi/2 on even L; C107's original Q denoted a different observable. Exact-ground-space gap stability on an ambient neighborhood is false at classical ferromagnetic coexistence under longitudinal fields, although Q_L=0. Use competing phase energies for coexistence. A gapped free Ising example has E_2=-sqrt(5)/2, E_4=-sqrt((5+sqrt(17))/2), positive doubling defect despite no Fermi points. Thus C130's free leg needs an M/band-norm bound, not only Fermi-point velocities if it includes gapped chains.

[Open cover] PROOF.md Section 4 is the requested table. C129 is relative U(1); C130 is drafted and relative to parity-even perturbations, not an ambient ball. Gapped energy convergence and uniform stability are separate inputs. The general frustrated case remains open. An exhaustive physical classification of the gapless set is also missing, so a residual-class obligation is retained. No analytical attack on rows (d) or (e).

[Literature checked] Hastings-Koma Theorem 2.8 proves clustering, not comparison of different ring sizes. Michalakis-Zwolak Theorem 1 needs frustration-freeness, Local-TQO and Local-Gap. Borgs-Kotecky Theorem 4.1 is classical and cubic-torus. Borgs-Kotecky-Ueltschi Theorem 2.1 explicitly assumes spatial dimension >=2 and cannot be used directly for chains. The needed fixed-circumference estimate is written as an O(beta exp(-N/xi)) bound on log Z, with constants uniform as beta->infinity. A relative O(exp(-N/xi)) error bounded independently of beta would incorrectly force exact linearity of E_N.

[Numerical] Small multistart optimization at L=2,3,4 and separate Heisenberg benchmarks through L=8. NUMERICS.md is the authoritative final numerical summary; raw outputs include termination status and numerical gauge brackets. Best even values found are 4 and 2.7912527285, at Heisenberg or a staggered-rotation image. These are empirical lower bounds on suprema, not certificates. Full multistart through 16 sites is deferred to submit_c131.slurm, 224 tasks, 32 starts for each L=2,4,6,8,3,5,7. Nothing was submitted or sent to Typhon.

[Unavailable numerical input] No C130 radius, ambient theorem, or surface sampling measure exists in the supplied draft, so a coverage fraction is undefined and was not manufactured.

[Files] RESULT.md summarizes conclusions. PROOF.md contains the lemma, table, gapped formulation and source audit. NUMERICS.md contains values and limitations. optimize_q.py implements the search; checks.py checks ring construction, gauge and exact baselines; pilot.py and refine_pilot.py reproduce small searches; write_numerics.py regenerates the report; submit_c131.slurm and summarize.py support the full cluster run. RUNNING.md explains reproduction.

[Next] Finish the relative C129/C130 estimates, prove a uniform ring-energy comparison in a stable gapped phase and a coexistence version, run the full cluster search, and preserve the unresolved transverse, multicritical and exhaustion obligations. The shared target's 2M gluing claim conflicts with C98's proved 4M; the reduction does not rely on the sharper value. C(2)<infinity is not proved.
