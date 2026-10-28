# C47 result

[proved, verdict] The requested infrared theorem is not established. This work proves finite and local reductions for the ten-letter demonstrator, including its exact XXZ anisotropy, but does not prove c=1 for that microscopic chain or c=2 for the n=1 wrapper. Full derivations and the assumption ledger are in [PROOF.md](PROOF.md).

## What is now proved here

[proved] In the invariant period-three sector, each 32-dimensional cell has an exact two-dimensional code and gap 2/5. The complete second-order logical bond is

    h2 = -20505/11648
         -(9/32)(sigma+ sigma- + sigma- sigma+)
         +(243/11648) ZZ.

After the staggered rotation, its spin-1/2 XXZ normalization is

    J = (9/16)lambda^2,
    Delta_XXZ = 27/182.

The proof includes both diagonal and hopping coefficients, the vanishing of cross-interface second-order terms, and the exact U(1) and charge-conjugation symmetries. At kappa=0 only the scalar -3/2 remains, as required by the scalar-resolvent obstruction.

[proved] A range-three generator gives an exact finite-order unitary identity

    U H U^dagger = D+R,
    P D P = lambda^2 sum_j h2_(j,j+1),
    ||R||_a <= C_a |lambda|^3.

PROOF Section 3 gives explicit constants and a full convergent BCH estimate, uniform in ring size and independent of temperature. D is locally block diagonal. R still couples the code to excited sectors, so this is not the requested unitary reduction to a spin chain.

[proved] A spectral embedding into one logical fermion and 60 massive auxiliary fermions preserves charge, grading and charge conjugation. Every auxiliary energy is at least 1/5. For interaction incidence norm w<1/5, the enlarged and original frozen-sector spectra coincide below a positive ground-relative threshold, and their free-energy densities obey

    0 <= f_original-f_enlarged
       <= T log[1+((2^61-32)/2) exp(-(1/5-w)/T)].

This comparison is uniform in volume and temperature. The embedding is an isometry into a larger space, not a unitary elimination of the massive sector.

## The resummed step, with its assumption exposed

[assumed, imported theorem M] The newer local note [FINITE_CELL_MATCHING.md](../../contrib/direct_attack_20261028/infrared/FINITE_CELL_MATCHING.md) supplies a finite-cell matching theorem using Pedra-Salmhofer Theorems 2.4 and 4.5. Its proof and independent local audit were read. The primary paper was not available in the supplied local material, and no network was used. Its determinant-expansion input and effective-action norm conventions are treated explicitly as an imported assumption.

[proved conditional on M] The demonstrator meets the finite-cell, vacuum-annihilation, coefficient-conversion and smallness conditions of M after choosing lambda last. PROOF Section 5 gives the choices, including a sufficient spatial weight a=768. Its zero-temperature effective logical action differs from the XXZ action by at most 4C_768|lambda|^3. Relative to the hopping scale this is at most (128C_768/9)|lambda|, with no volume-dependent coupling choice.

[proved, limitation] This is a retarded action with general many-body vertices. It is not a finite-range spin Hamiltonian. The finite-temperature baseline also contains thermal auxiliary terms. The comparison does not supply the interacting T^2 coefficient or a complete current/stress scaling limit. It does avoid expanding the order-one excited-cell logical dependence in the raw Grassmann series that failed in C46.

## Why no published XXZ theorem is invoked

[proved, application boundary] No Benfatto-Mastropietro theorem is invoked as proving this model's asymptotics. The available sources summarize relevant results but do not supply a complete primary theorem statement for this retarded, arbitrary-degree action. PROOF Section 6 lists every checked input and missing hypothesis.

[proved] Lowering lambda does not make the reference anisotropy smaller. Its density-interaction/hopping ratio remains 27/91. A theorem valid in an unspecified sufficiently small neighborhood of the free XX chain cannot therefore certify Delta=27/182. Half filling and its umklapp also need treatment in the actual theorem. These are application gaps, not evidence against the numerical Luttinger liquid.

[assumed, remaining infrared input] A suitable constructive RG theorem must cover the matched action, its thermal terms and microscopic sources, and prove either the thermal coefficient with identified velocity or the complete compact-boson scaling limit. Neither a local O(lambda^3) norm bound nor selected correlation exponents alone supply that conclusion. The error in the available direct pressure comparison stays fixed as T tends to zero.

## Full-chain and wrapper scope

[proved] The period-three type sectors are invariant, and their selection at lambda=0 follows from an explicit shared-port bound. At nonzero lambda, the supplied run lower bounds alone do not prove a ground-relative cost of dilute defects. A positive per-site separation from every other type configuration is impossible uniformly in volume: changing one site's type costs at most 4||h||. A positive cost per defect remains possible. PROOF Section 8 keeps this separate from frozen-sector infrared matching.

[proved] The n=1 wrapper still has an unbounded microscopic cell-length space and two retained modes. A finite length cutoff is not invariant under marker motion, and the finite-cell embedding constants are not controlled uniformly as that cutoff grows. A single XXZ theorem would also not cover the interacting payload and length fluids. The new finite-cell argument therefore does not establish its c=2 conclusion.

## Checks and final status

[numerical] `python3 -B checks.py` passes. Exact Fraction fixtures verify the second-order formulas. Independent floating-point matrices of dimension at most 64 give cell gap 0.4, second-order discrepancy below 1.8e-15, and the expected cubic remainder for the explicit unitary. Results are in [checks.json](checks.json). No many-cell central charge was measured, and the conservative analytic bounds do not certify the previously simulated couplings .05, .10 or .15.

[numerical, imported evidence] The older demonstrator's c near one remains numerical evidence. It is not upgraded to a theorem by this work.

[proved, final boundary] The new unconditional results are the finite algebra, the uniform finite-order transformation, the symmetry-preserving embedding and its spectral/thermal comparison, and zero-coupling type selection. The resummed logical action estimate is conditional on M. The requested all-order unitary reduction, the published infrared theorem application, c=1 for the full ten-letter chain, and c=2 for the n=1 wrapper remain unproved.

[proved, work record] PROOF.md, RESULT.md, PROGRESS.md, STATE.md, checks.py and checks.json are confined to work_C47. No network, subagents, cluster jobs or writes outside this folder were used. The project-level STATE.md was left unchanged under the explicit write restriction.
