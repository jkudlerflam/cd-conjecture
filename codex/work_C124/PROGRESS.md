# C124 progress

2026-10-30.

[Proved, inherited] Read the background with all updates, reconnection specification, C121, C123 and C118 derivations, and requested report updates. The independent parallel attempt has not been accessed.

[In progress] Prove reconnection for commuting nearest-neighbor bonds in a basis adapted to their interaction algebras. The separator is the pair of independent virtual factors at each measured site.

[In progress] Audit C123's temporal-transfer boundary factor. A diagonal three-level on-site family makes that sufficient boundary estimate grow with L although G is nonpositive. The full thermal transfer spectrum must be retained.

[Planned numerical checks] Reuse locally copied finite-grid XX determinant and full-spectrum ED kernels, calibrate them, and check all new applicable identities on XX, Heisenberg, XXZ and SU(3). All computations and outputs remain in this folder. No cluster or network shell commands.

[Environment] The command python is unavailable. Use python3.

[Proved] For commuting adjacent bonds, the finite-dimensional interaction-algebra decomposition gives an adapted on-site basis and a nonnegative matrix B on the full recorded histories such that q_s=B_ss and q_D(s,t)=B_st B_ts. Therefore G<=0 at every spacing. A separate transfer matrix on central sectors gives R_L<=log d, hence G<=(2N+1)log d in any measurement basis.

[Proved] The sufficient boundary-overlap estimate in C123 is false uniformly, even in a diagonal on-site d=3 family. With k=diag(0,L^(-2),1), h=(k tensor I+I tensor k)/2, M=1/2 and tau=2aL, its boundary factor is log z(2tau)+(2L-2)log z(tau), asymptotic to 2L log 2. The actual G=log z(2beta)-2log z(beta) is nonpositive. Tensoring with stoquastic XX preserves a nontrivial conformal infrared sector and the obstruction.

[In progress] Running small full-matrix commuting checks, fresh large-ring XX determinant checks, and full-spectrum interacting ED. The requested general noncommuting stoquastic lemma remains open.

[Completed, proved] Wrote the adapted-basis commuting theorem, the arbitrary-basis bound, the diagonal boundary-factor obstruction, its critical stoquastic d=4 extension, and the full thermal increment reduction in PROOF.md. The general noncommuting case remains open. RESULT.md states that limitation first.

[Completed, numerical] Main survey: 12 commuting quantum cases, 72 diagonal boundary cases, 16 dense/free-fermion calibrations, 46 large-ring XX cases and 90 ED cases. Additional surveys: nine direct-sum commuting cases, nine rotated-basis cases, 45 critical-spectator cases, 156 finite-interval cases, and 23 temporal-transfer/diagonal checks. Some surveys overlap. All applicable assertions pass. The attempted extension of G<=0 beyond commuting bonds fails in every requested model family at sampled parameters.

[Completed, numerical] Maximum commuting-separator relative residual 3.50e-14. Dense/free-fermion recorded-history residual 8.92e-15. Direct temporal-transfer versus charge-block ED residual in G 1.87e-13. Retained all spectral sectors. The largest charge block has dimension 924. No large interacting extrapolation or SU(3) length-12 calculation is claimed.

[Completed, audit] All Python sources parse, numerical output files load, and Markdown files contain no em-dash or en-dash. Checked the algebraic decomposition, endpoint multiplicities, reflected-history factorization, optimized M in both counterexample families, tensor factorization of the full top eigenspaces, and the ceiling's time-slab band. Updated local STATE.md only, respecting the prohibition on external writes.
