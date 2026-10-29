# C71 result

[assumed, target unresolved] Lemma TP remains unproved for weak off-half-filled XXZ with C65's smooth trace-compatible regulator. No counterexample to TP or to the interacting thermal coefficient was established. The next obstruction is now explicit at the ultraviolet preparation step, and concerns the available absolute majorant rather than the actual coefficient norm.

[proved obstruction to a specified bound] Let a_Lambda be the smoothly regulated chiral propagator after a fixed smooth infrared subtraction. For every fixed delta>0,

    integral_0^delta dt integral dx |a_Lambda(t,x)|^2
                    = log(Lambda)/(4 pi A) + O(1).

This is the absolute two-line density loop at order g^2. Taking absolute values before the internal spatial integrations gives the same logarithm times ||v0||_1^2 in a generated four-field kernel. The zeroth-moment majorant already diverges. Smooth annular Gram bounds and a fixed field-weight reserve cannot make this particular ultraviolet estimate uniform. See [PROOF.md, Section 3](PROOF.md#3-a-logarithm-in-the-absolute-ultraviolet-tree-bound).

[proved cancellation and partial advance] The integrated loop is finite. At zero temperature, throughout the interaction's transfer support and for every sufficiently large Lambda,

    Pi_Lambda^0(Omega,p) = p/[2 pi(A p-i Omega)].

After spatial interaction smearing and the infrared subtraction, this coefficient is exactly cutoff independent and has every fixed absolute spacetime moment. The short-distance identity B_full=-(i/(2 pi)) partial_x c moves a derivative onto the smooth spatial potential. It is essential to retain this cancellation before taking the coefficient norm. This differs from C65's sharp-window result, where the actual kernel failed even the first moment.

[proved thermal coefficient estimate] For this prepared density-loop quartic, including two speed derivatives,

    ||partial_A^j(K_Lambda^T-Per_beta K_Lambda^0)||_(m,beta)
                                      <= C_(m,P) T^P,
    j<=2,  P>1,

uniformly in Lambda at the fixed preparation scale. Equal-image singular products cancel in the thermal difference. The remaining distinct-image products have integrable singularities and arbitrary polynomial decay. Thus m=5 and P=3, with 0<theta<1, work for this coefficient. [PROOF.md, Sections 4-5](PROOF.md#4-the-same-coefficient-after-the-loop-is-integrated) give the derivation. This is not an all-order effective-action estimate or an all-scale q_h^theta defect bound.

[assumed, precise next work] Prove the cutoff-independent analytic ultraviolet preparation (18) with its loop cancellations, full external-field kernels, fixed moment reserve, thermal derivatives and complete scalar. Then prove the common-domain lattice/reference bounds (19)-(20), retaining both q_h^theta and (T/q_h)^P. The existing four-history closure and shared-short-memory scale sum only apply after these estimates are established. Trace convergence alone supplies neither the prepared-kernel topology nor finite-regulator infrared matching. The latter logical distinction is stated without asserting an actual regulator-induced mismatch.

[numerical] `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B checks.py` passes. The measured doubling increment of the logarithm is 0.04242992310627913, versus 0.042429923106279153 from the proof. Dynamic-bubble Matsubara checks have maximum extrapolation discrepancy 2.82e-11; the image identity residual is 7.95e-17. These are coefficient and covariance controls, not an interacting XXZ RG calculation.

[proved scope] PROGRESS.md, PROOF.md, RESULT.md, SOURCE_AUDIT.md, STATE.md and the checks are inside work_C71. No files outside this folder were modified. No retarded-action or multichannel theorem, all-order vacuum estimate, or independent review is claimed.
