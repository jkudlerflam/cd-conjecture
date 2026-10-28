# C36 progress

[proved, task record] Read the project STATE.md and the five requested source notes. No prior PROGRESS.md existed. All writes and computations will stay in work_C36, with Python, NumPy and SciPy only and no network. The project STATE.md cannot be updated under this task restriction.

[assumed, inherited finite inputs] The complete initialized compiler/core remains an external input. This task studies its infrared consequences conditional on the exact finite-cell contracts and the wrapper bounds in C8.

[proved, imported result] The Gaussian certificate in contrib/digit_compilation/audit/GAUSSIAN_CERTIFICATE.md supplies the full matrix inequality kappa=sqrt((1-delta_G)/(1+delta_G)); delta_G<5/13 implies kappa>2/3. It does not bound microscopic corrections.

[heuristic, active investigation] Conserved marker order should identify the mobile legal sector with a constrained tensor product of direct-sum cell Hilbert spaces. The low projector then becomes a fixed on-site projector onto the two selected lengths and their code states. This is the route being tested for a controlled local elimination. A global spectral-band argument is excluded.

[proved, outstanding questions] Need distinguish bounds valid at finite perturbative order from a convergent all-order local expansion, classify non-oscillating hexagons with the stated edge ordering, account for the collective N*g*log(1/delta_pay) tolerance, and give fixed-p assumptions and limits in RESULT.md.

## Derived and checked

[proved, conditional on the finite core] Ordering conserved markers identifies the legal sector, on a rooted cover, with a fixed-length-constrained tensor product of direct-sum open-cell Hilbert spaces. Cell defects are fixed on-site projectors in these coordinates. The exact species and marker chemical fields commute with every microscopic term and do not enter virtual denominators within a charge sector. The remaining perturbations are bounded nearest-cell-boundary operators, including length transfers.

[proved] Second-order connected kernels have local norm at most C*s^2*B_c^2*exp(-delta_* tau); their time integrals are bounded by C*s^2*B_c^2/delta_*. Every fixed higher order has finite local constants. A positive majorant gives geometric and exponential bounds for the restricted histories in which every insertion touches a currently active defect. This does not include all linked histories of the microscopic model.

[proved, obstruction] An unrestricted interaction-picture operator-norm version of IR_AUDIT (10) is false, even for an immobile defect coupled locally to a finite-range low chain. A domain-wall Jacobi chain gives a kernel matrix entry growing at rate s*(lambda_max-lambda_min)-delta_*, while the corresponding ground-state entry decays at least as exp(-delta_* tau). The analytic example will be given in RESULT.md. It locates the missing energy restriction and linked-background cancellation, rather than an actual infrared instability.

[proved] An undressed six-cycle is non-oscillating precisely when its six distinct edge labels partition into equal-sum triples. There is at most one such unordered partition. The alternating creation/annihilation pattern and the chirality signs are independent. Full integer momentum congruence agrees with this criterion.

[numerical, p=2,3,5,7] checks.py enumerated respectively 28, 234, 3875, 26068 hexagons and found 9, 110, 2290, 17556 resonant hexagons in the stated natural ordering. Complete p=2 and p=3 lists are saved. All 64 chirality assignments were compared for every p=2 and p=3 cycle. Additional finite checks used Jacobi matrices of dimensions 5,9,17,33 and positive majorants of dimensions 8,32,128. All checks passed, in under one second. They do not test the microscopic compiler or its finite-coupling stiffness.

[heuristic, current conclusion] The RG implication can be quantified, including N*g*log(1/delta_pay), but the full low-frequency elimination remains an assumption. A final report will distinguish the proved restricted-excursion bounds from the missing full connected-background bound, and state acceptance conditional on resolving it and completing the core.

## Completion

[proved, deliverable] RESULT.md is complete. It proves the fixed cell-coordinate representation, the second-order bound, and geometric/time bounds for the explicitly restricted active-defect histories. It supplies two exact diagnostic obstructions to naive extensions: an interaction-picture all-band kernel with poles accumulating at s=0 as L grows, and independent dressed two-level systems whose global bare-Q threshold tends to zero. Neither is an infrared instability. The corrected low-external-energy bound is proved for the first example.

[proved, hexagon refinement] At nu=1/2 in the natural point-major edge ordering, there are no odd-h dressed six-cycle resonances. A six-cycle has two labels in each of three distinct point rows; this prevents one label from equaling the sum of the other five. Thus the equal-triple hexagons and their even-h dressings exhaust this class. The updated integer checks include the absence of odd-h resonances at p=2,3 and pass.

[heuristic, final status] The Gaussian certificate and RG majorants give a strict kappa>2/3 neighborhood, with N*g*log(1/delta_pay)<<1 and separate filling-shift control. One fixed nonzero s can satisfy all tolerances at fixed p if the local low-frequency estimate E holds. The full retained-background linked estimate E and the complete compiler/core remain unestablished. The candidate is acceptable at the physics standard conditional on closing those gaps and the standard compact-boson interpretation, but this task does not certify a complete refutation.

[numerical, verification] All checks pass after the odd-h refinement. Graph sizes are p=2,3,5,7; Jacobi dimensions are 5,9,17,33; positive-majorant dimensions are 8,32,128. Runtime is about one second. No full microscopic chain, dressed stiffness, or central charge was computed. source_manifest.json records hashes of the sources consulted. No source outside work_C36 was changed.
