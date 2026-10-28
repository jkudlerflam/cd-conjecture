# C34 progress

2026-10-28. Completed. Read RESULT.md for proofs, scope and next tests.

[proved] New domain-wall gap bound. For S=diag(a,a,c), c<0, D=(d,0,d_z), B=(0,0,b), set kappa=|c|-sqrt(a^2+d_z^2)-2|d|. At b=0,kappa>0, exactly two polarized ground states and gap >=4 kappa for every periodic L>=4. At b!=0,kappa>=0, unique polarized ground state and gap >=2|b|. The proof is the finite operator inequality H-cL>=2 kappa W+b sum Z, W=sum(1-ZZ)/2. It closes dimension-3 TR-only and reflection-only subsets and a dimension-4 symmetry-free subset, all including models outside the supplied common-axis certificate.

[proved] At B=0 every uniform product eigenvector requires a uniaxial S. Thus the simple-spectrum dimension-5 TR-only region has no such eigenvector. This is not a statement about dressed or nonuniform vacua.

[proved, conditional on the supplied theorem] Direct Bravyi-Gosset-Liu-Wong applicability is equivalent to a nonzero projection of B onto E_-=ker(J+||J|| I) intersect ker(J^T+||J|| I). Optimal resulting gap bound is 2||P_E B||. Arbitrary site-dependent product rotations cannot enlarge this criterion: average their longitudinal axes around the periodic chain. With D nonzero the only possible uniform axis is D. Thus the criterion fails on TR-only, reflection-only with D nonzero, and generic misaligned-D symmetry-free bonds.

[proved] RESULT Section 5 gives an exact algebraic envelope after those two certificates, with seven disjoint rows and maximal dimensions 5,3,5,4,8,6,5. The overall maxima remain 5,5,8 for TR-only, reflection-only and no internal symmetry. The actual unresolved critical set is not known; controlled continuum neighborhoods lack bare-coupling radii. No finite scan supplies a uniform exclusion.

[proved] Uniform-product threshold for the exact-vacuum two-fluid slice is mu_prod=-2 cos(phi)+(A-sqrt(A^2+4|g_L+g_R|^2))/2, A=V+2 cos(phi). Above it a product state beats the vacuum. This does not exclude a finite-density two-fluid phase.

[numerical] checks.py passed in under one second with NumPy/SciPy only. Nine complete spectra check the inherited certificate at L=4,6,8; eighteen complete operator inequalities and spectra check the new wall bound at L=4,6,8, including both field signs and kappa=0 with nonzero field. Worst negative remainder eigenvalue -3.6e-15. Also checked seven symmetry examples, a rotating-axis certificate, 200 transverse-block norms and 100 two-site product expectations. No new unknown c estimate, ED scan, or DMRG run.

[numerical, inherited evidence audit] Multicritical hunt uses rings up to L=20 and VUMPS chi<=64, but its tilted-field Neel line is unresolved at L<=16. It does not cover generic TR-only, nonzero-D reflection-only, or no-symmetry strata. The two-fluid hunt has converged one-mode examples through chi=48, but S1 and S2 remain unresolved and its dressed-vacuum extension is open. Read RESULT Sections 6-7 for the distinctions.

[heuristic] Full proof not currently within reach from the supplied exclusions. Needed: microscopic exhaustion of complete low-energy descriptions, including neutral sectors, non-product CFTs and singular neighborhoods. Also the unitary Z2 family remains open in C15. Next bounded tasks are the unresolved tilted-field junction, soft-wall tracking from the new bound, S1/S2 convergence, and neutral-tower response tests in the generic symmetry-free class.

[proved, workspace record] Read project memory first and all requested sources. No network, no subagents, no writes outside this folder. Parent STATE.md intentionally unchanged under the task rule. Local STATE.md is the handoff. Files are RESULT.md, PROGRESS.md, STATE.md, checks.py, checks_results.json, checks_output.txt. Reproduce with python3 checks.py. Unrefereed.
