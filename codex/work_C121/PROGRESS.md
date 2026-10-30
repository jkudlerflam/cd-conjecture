# C121 progress

2026-10-29.

[Work record] Read the project STATE.md, all updates in the bound background, REPORT sections 9.68 through 9.80, and the C116, C118 and C120 result, proof and handoff notes. The project-level memory is not edited because this task forbids writes outside work_C121. All new Python and output stay in this folder.

[Proved] For a fixed Hamiltonian, a seam pinching inserted into positive imaginary-time evolution loses at most a factor d per insertion. A unitary-average expansion and Schatten Holder inequality give the complementary upper bound by the unmeasured partition trace. The comparison holds at arbitrary spacing, without a spectral gap or recovery assumption.

[In progress] Formulate the positive recorded histories, the Cauchy-Schwarz step and the remaining reconnection functional. Coarse measurements do not make the seam a spatial separator between measurements. Test the proposed positivity step on XX, Heisenberg and SU(3), and distinguish coarse pinching from the continuous-measurement limit.

[Open] No general ratio bound or velocity-weighted theorem has been obtained.

[Proved] Completed the arbitrary-spacing pinching theorem, the coherent and positive-history Cauchy-Schwarz inequalities, and the exact reduction R <= 3N log d+G. Stated the remaining reconnection lemma with explicit constants in the conditional bound on A. A three-site stoquastic witness explains why endpoint seam labels do not factor the two sides.

[Proved] Continuous pinching converges to the Gibbs trace for D(H). On XX its ground-energy defect is 1/2+(1/2)tan(pi/(2L)), tending to M. Derived a full finite-temperature one-pinch XX formula retaining both Jordan-Wigner sectors.

[Numerical] The no-cost reconnection inequality fails on all three requested models at tau = L/M, N = 1. Large XX checks through L8192 also violate the comparison with the unmeasured doubled ring. All proved channel bounds pass. Full details and limitations are in CHECKS.md.

[Work record] An initial summary-file serialization error was fixed by converting a numpy count to a Python integer. The final suite completes and writes both output files. Read primary literature on spatial approximate Markov chains and complex-time locality using the browsing tool only. These results do not supply the missing temporal-history reconnection estimate.

[Open] The universal ratio bound remains unproved. RESULT.md, PROOF.md, CHECKS.md and the local STATE.md record the theorem obtained, the false stronger step and the precise remaining lemma. No independent review or external file edit was performed.
