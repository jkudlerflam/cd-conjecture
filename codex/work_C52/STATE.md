# C52 local state

[proved, task status] Finished the exact-cell redesign audit. No successful fixed-d nearest-neighbor construction and no reduction of S to a published theorem. RESULT.md summarizes the outcome; PROOF.md contains all finite arguments and their scope.

[proved, inherited theorem] C8 RESULT.md Section 10 proves that an exactly invariant product-cell code with separated boundary regions has commuting adjacent connected logical bonds. Do not treat this theorem as a new C52 result.

[proved, new application] For number-conserving logical bonds, the outer-cell block of their three-cell commutator is T_left T_right. The C39 target has invertible T=-W[Gamma(w)+i eta Gamma(wv)], with sigma_min(T)>=W lambda_min(Gamma(w)). Genuine normal-ordered quartics do not change this one-particle matrix element. Exact invariant product-cell implementations of this target are impossible. The same proof applies to diagonal fixed-marker blocks of an invariant direct-sum mobile code.

[proved, C44 marker leakage] With the notation in PROOF Section 5, the directed head output has norm squared s_b, overlap A_b squared, and leakage squared s_b-A_b^4>=s_b(1-s_b)>0. The full wrapper factor multiplying it is nonzero. A boundary-local repair retaining the sine-head ground states also fails if it leaves two bulk head positions untouched. Changing only the payload parent cannot fix this.

[proved, further bounds] Single-cut compression gives 2 rank(T)<=d^2, or rank(T)<=floor(d^2/4) with integer-charge-preserving encoding. An n-orbital nearest-neighbor Bloch matrix with rank-r hopping has at most 2r isolated simple roots. Exact dressing of physical radius R remains subject to the separated-endpoint theorem for cell lengths greater than 4R+2.

[assumed, remaining construction] Cross-boundary non-product codes and dressing with physical range comparable to n are not excluded. No literal bond, wrapper replacement, exclusion proof, or theorem-level infrared application for such a design was obtained. Do not describe these restricted obstructions as a proof that every redesign is impossible or that the original conditional phase is false.

[numerical, reproduction] Run OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B checks.py. All finite algebra controls pass, largest matrix 64 by 64. checks.json contains results. No microscopic CFT computation was performed.

[proved, work scope] All writes are inside work_C52. The parent project STATE is outside the permitted write scope. No network, agents or cluster jobs were used. Neighboring S runs were consulted through their progress files; their work was not duplicated.
