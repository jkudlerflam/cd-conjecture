# C38 progress

[Proved, scope] Read the requested background, status, 23-letter candidate, C32 results and proofs, C31, C26, C20, infrared control, and C37/C33 progress. No previous C38 progress existed. All research artifacts are in work_C38. The project STATE is outside the permitted write scope.

[Proved, inherited] C31 already proves that separated hopping eigenvalues and a common nonzero field can remove the centered reference's non-pair resonances. It explicitly excludes pair transfer. C32 proves all-n crossings, not distinct centers. C33 now reports a certified n=23 finite-loop witness. C37 is pursuing algebraic nonidentity certificates through n=23. C38 will not repeat those computations.

[In progress] Test the hierarchy implication itself. The key distinction is bandwidth separation versus momentum separation. Derive scale-weighted Schur bounds for cross-flavor mixing, a quantitative lacunary target lemma covering the complete Fermi-point list, and the extra hypotheses needed to apply either to the frozen and mobile cores.

## Located obstruction

[Proved, derivation in preparation] In the large-stub-mass short-path limit the Gram matrix is a^(-3) S C S, S_jj=a^(-2j), with C_ij=(2i+2j+2)!/[(2i+1)!(2j+1)! 2^(2i+2j+3)]. For the first two indices C_01^2/(C_00 C_11)=2/5. The cross-scale Schur correction is therefore a fixed fraction of the narrow coefficient, not exponentially small in the scale separation. Adjacent normalized squared correlations tend to one at large index. Diagonalizing the Gram matrix removes its cross terms, but gives no smallness or nonidentity result for the independently induced loop and field matrices.

[Proved, elementary] More generally h=S F(k) S has exactly the same Fermi zeros as F(k) for every positive diagonal S. Arbitrarily strong energy hierarchies can coexist with either exact resonances or nonresonance. Cross-flavor elimination requires a normalized spectral gap at the relevant momentum, not a bandwidth gap.

[In progress] Write a constructive all-n lacunary endpoint target and explicit perturbation bounds, keeping its realization by the 23/28-letter cores as a missing hypothesis. This avoids repeating C33/C37 finite certificates.

## Completed analysis and checks

[Proved] PROOFS.md gives the natural Gram obstruction, the congruence and Schur-complement argument, and the mixed loop/antisymmetric-hopping center response. It gives an explicit full-endpoint lacunary target with margin d_n=(7 pi/256)16^(-(2n-1)). A scale-weighted C1 matrix error at most d_n^2/100000 preserves exactly 2n simple roots, every pair-sum margin, and a velocity lower bound. This is a conditional target theorem, not a realization by the candidate cores.

[Proved, computer-assisted check] Fraction arithmetic checks the target at n=1,2,3,4,8,16,23,32,64,128, up to 256 endpoints and 32896 pair sums. It verifies the Gram correlation formulas through index 63 and the scale-independent two-flavor Schur correction at four ratios. No finite check substitutes for the all-n proofs.

[Numerical, n=2] Ten 8-dimensional Lyapunov solves check the frozen/mobile large-mass Gram limits. Both approach the proved squared correlation 2/5. The full microscopic chain was not diagonalized.

[Proved, mobile conditional arithmetic] For the designed endpoints, every length density in [1/5,3/10] has mixed momentum margin at least 59 pi/160. C32 residual length-pair bonds are expanded into all four occupation-polynomial coefficients, with a total triangle bound 9 max||E_pq||. Their scalar components are retained. Matching these operators to dressed continuum quantities remains assumed.

[Proved, coordination] Re-read C37/C33 progress and read C39/C40 progress once available. C37 has finished certified frozen targets for n=1,...,23. C39 is pursuing a position ruler. C40 is working on the RG and massive-sector budget. This report does not duplicate those tasks.

[Proved, scope accounting] No network, subagents, large computation, or project-file edits outside this folder were used. One accidental empty scratch file was created in /tmp and immediately removed. No research output was stored there. The project STATE is outside scope; a local STATE.md supplies the handoff.

## Final handoff

[Proved, completed] RESULT.md and local STATE.md are written. The proof constants, saved rational results, and numerical output were reviewed. The all-n arithmetic target and perturbation theorem are proved; their realization in the frozen/mobile cores is not. The hierarchy-only argument fails at the relative cross-scale estimate and at the inference from bandwidths to loop-dependent center differences.

[Proved, latest coordination] Final progress reads find C33 and C37 completed, C39 now checking a two-bus plus biased-vacuum proposal, and C40 deriving connected finite-order bounds and a velocity-dependent RG budget. No part of those ongoing constructions is claimed as proved here.
