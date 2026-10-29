# C96 progress

[Scope] All writes and Python execution are confined to work_C96. No network commands or delegated agents are used. The project memory and all four requested predecessor folders have been read, including check implementations and retained outputs. The project STATE.md is outside the allowed write scope and will remain unchanged.

[Imported] C95 controls the zero material-momentum contact response and density relaxation. Its invariant-subspace isometry does not identify the tagged kernel with the physical local density-product kernel. C94 requires a tagged kernel and a separate relative thermal estimate. C93 requires a dynamical dressed clock statement in addition to either static input.

[In progress] Deriving a full finite-ring tagged-contact spectral formula using Slater determinants with material-momentum-dependent twists. Investigating a weaker kernel norm for the product-cell comparison, based on cancellation of the equal mean in translated selected cells. The full thermal remainder will be treated with the frozen material Hamiltonian retained at its actual coupling.

[Proved, completed finite formula] Constructed and normalized the full twisted-Slater gap eigenbasis. Proved completeness by a subset-orbit count. Derived the tagged-contact form factor by Cauchy-Binet with two pinned physical positions. Its zero material mode reduces to C95; nonzero modes use different twists and cannot be restricted to two particle-hole replacements on one grid.

[Numerical] Checked the full basis through N=4 and the contact kernel through N=8 against independent sparse gap resolvents. The largest kernel discrepancy is below 3e-15. The specified density subtraction produces moments near 0.1102/t at half density. No asymptotic inference is made.

[Proved, new comparison criterion] Split the product-cell source Gram matrix into means and connected within-cell covariances. The selected mean minus its cell average has a uniformly bounded discrete primitive on each good run. Summation by parts transfers two differences to the kernel. The resulting norm is sum|k|+b^2 sum(1+dist)|Delta k|+|c|. It accommodates a smooth nonoscillatory r^-2 tail and is insensitive to a canonical constant except through its bounded absolute row sum. Resonant oscillatory tails still need stronger control.

[Proved, thermal coefficient] Derived a positive-matrix bound for imaginary-time material covariances, including complex noncommuting sources. This extends the weaker coefficient comparison to isolated product thermal cells with their full internal traces retained at the actual J/T. It does not compare the nonzero-port trace or an interacting bath response.

[Proved, remainder diagnostic] An exact two-spin bath with vanishing frozen material Hamiltonian still has a low-temperature divergent fourth-order coefficient. Its full residual is uniformly quadratic in J, but the quadratic free-energy coefficient misses its leading low-temperature value by a factor of two. Replication preserves this after the thermodynamic limit. Retaining frozen dynamics alone does not resum induced material splittings.

[Check correction] The first auxiliary norm check stopped at weaker_norm<36, with actual value 39.0379841724 at N=32. The provisional threshold omitted the absolute row sum. Telescoping gives the correct analytic bound 36+pi*coth(pi)<40. The corrected check and all other checks pass. Details and tables are in CHECKS.md.

[Open] No uniform estimate of the exact tagged spatial norm, no full W-relative quantum thermal remainder, and no C93 dressed contact-clock theorem were obtained. RESULT.md leads with this limitation. The exact formulas, the weaker sufficient criterion, and the thermal diagnostic are the new results.

[Completed deliverables] PROGRESS.md, RESULT.md, PROOF.md, CHECKS.md, local STATE.md, two reproducible Python scripts and their small outputs. The project-level memory remains unchanged because writes outside this folder are prohibited.
