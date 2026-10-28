# C47 local handoff

[proved, outcome] The requested microscopic infrared theorem is not closed. Read RESULT.md, then PROOF.md. The project-level STATE.md cannot be updated because writes outside work_C47 were forbidden.

[proved, fixed data] Ten-letter demonstrator at D=1,kappa=3/5,Jp=2,Jz=-1 has cell gap 2/5, code dimension two, second-order scalar -20505/11648, hopping -9/32, ZZ coefficient 243/11648, and staggered XXZ anisotropy 27/182. Full derivation is in PROOF Sections 1-2.

[proved, local theorem] PROOF Section 3 gives a range-three anti-Hermitian generator and explicit C_a for ||R||_a<=C_a|lambda|^3, uniformly in volume. D is locally block diagonal but R still leaks. Do not describe this as an exact spin-chain elimination.

[proved, new embedding] Pair the 15 even excited cell eigenstates with their charge-conjugate odd states. Four auxiliary modes per pair give 60 massive modes plus one logical mode, preserving U(1), grading and conjugation. Auxiliary gap floor is 1/5. PROOF Section 4 gives the exact isometry and the positive invalid-sector cost and thermal comparison.

[assumed, imported input] The finite-cell massive-matching theorem in contrib/direct_attack_20261028/infrared/FINITE_CELL_MATCHING.md is used as assumption M. Its local proof and audit were read; the primary Pedra-Salmhofer theorem was not available offline. PROOF Section 5 verifies the demonstrator application conditional on M. Use field weight 2+4sqrt(60), spatial norm weight 768, temporal rate 1/10, and choose lambda after those constants. The relative zero-temperature action error is at most (128 C_768/9)|lambda|.

[assumed, still missing] A primary perturbed-XXZ theorem verified for arbitrary-degree retarded vertices, its quantitative inclusion of Delta=27/182 or a different in-family parameter choice, half-filled RG treatment, sources, thermal winding effects and a c-identifying conclusion. No such theorem is invoked here. Lowering lambda does not reduce the reference anisotropy.

[proved, scope warning] Frozen-sector invariance is exact. Full-chain selection is proved here only at zero coupling. The README's run lower bounds alone do not establish a dressed dilute-defect gap at nonzero coupling. A fixed positive per-site separation from all other type sectors is impossible, because changing one site costs O(1).

[proved, wrapper boundary] The n=1 wrapper's unbounded cell lengths prevent a direct uniform use of this finite embedding; a length cutoff is not invariant. Two coupled retained modes still require a two-mode infrared theorem and physical-ensemble transfer. No c=2 proof was obtained.

[numerical, checks] python3 -B checks.py passes in a small local run. Uses Fraction, numpy and scipy only, matrices at most 64 dimensional. The first run had an array-shape bug in the new control, repaired before the successful runs. No network, external writes or delegated agents were used.
