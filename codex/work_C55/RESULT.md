# C55 result

[proved, verdict] The single-channel retarded Luttinger theorem remains unproved. This run establishes an all-scale theorem for the Gaussian linear response to arbitrary-degree initial kernels and completes the off-half-filled finite initialization. It states a precise nonlinear endpoint/source lemma sufficient to finish the theorem. No Luttinger liquid or central charge is claimed for the port chain at fixed nonzero coupling.

## New results

[proved using imported finite-integration results] The collaborator's UV initialization extends uniformly to p_F in any compact interval away from the band edges and half filling. For e_p(k)=2(cos p-cos k), a safe UV determinant constant is 2 sqrt(1+exp(16)). A compact-patch crossover has a smooth covariance with every fixed polynomial moment. C50's finite-integration theorem then absorbs arbitrary-degree retarded vertices through the first infrared shell. The surviving local terms are the chemical-potential constant, time and velocity derivatives, and one forward quartic. See [PROOF.md](PROOF.md), Sections 2–4.

[proved, all-scale linear endpoint theorem] Let A_h V be the coefficient linear in V of log E exp[V], after integrating free shells above h. For a sufficiently large fixed input field radius H and a smaller output radius rho,

    ||(A_(h-1)-A_h)V||_(rho,s)
                  <= C gamma^(2h) ||V||_(H,s+1),       h<=0.

The bound includes every initial degree and every number of contractions inside that one vertex. Its proof uses a determinant cofactor interpolation with one marked shell contraction. The linearized shell is odd, so g_h^lin(0)=0; curvature and displacement give |g_h(z)|<=C gamma^(2h)(1+|z|). One initial kernel moment pays that displacement. The remaining determinant has a uniform cumulative Gram bound. Field radius is spent once, not at every scale. See PROOF Sections 5–6.

[proved, linear counterterm] The induced quadratic and quartic coefficients converge with tails O(epsilon gamma^(2h)). Subtracting the limiting quadratic constant gives a dimensionless relevant term O(epsilon gamma^h). The localized coefficients depend continuously on the trial Fermi point, and an intermediate-value argument gives a first-order root shift at most 2C epsilon/v_min. These are statements at the Gaussian linearized level, not a self-consistent nonlinear counterterm. See PROOF Section 7.

[proved conditional implication] A joint finite-dimensional Brouwer map constructs the marginal history, relevant history and physical Fermi point from continuous beta functions satisfying the explicit bounds (28)–(29). Diagonal compactness removes the infrared cutoff. Unlike the collaborator's contraction construction, this existence variant needs neither a history Lipschitz bound nor a unique continuous counterterm branch. It still requires the uniform absolute beta estimates and exact UV matching. See PROOF Section 10.

## The remaining lemma

[assumed, not established] Lemma N in PROOF Section 9 requires the normalized renormalized-tree majorant for arbitrary retarded initial endpoints, with no cumulative field-radius or moment loss. Trees reaching an initial irrelevant endpoint must carry a summable short-memory gain. History marks and density/pair source marks must preserve the required estimates. The precise beta inequalities, normalization condition, source-ratio estimate and two-insertion remainder bound are displayed there.

The one-endpoint theorem does not prove this lemma. With several vertices, their relative centers are integrated by contraction lines, so a marked line's displacement cannot be paid by a moment of one initial kernel. This is where nested localization and the reference beta cancellation enter. The collaborator's polynomial tree norm pays a bounded number of factors per edge, but a compatible bound on those loads through the full renormalized construction still has to be proved.

[proved conditional implication] With Lemma N, reference matching transfers the density and pairing scaling dimensions K and 1/K. Their two-point decay powers are 2K and 2/K. The free pair overlap is -2i sin p_F, nonzero in the specified interval. Source comparison, rather than a presumed instantaneous current Ward identity for the retarded action, is needed to transfer these statements. See PROOF Section 11.

## Model and scope

[imported] The nearest-neighbor choice is the collaborator's Jp/Jz tuning. C53's diagonal BPPB counterterm has physical range four. Massive matching and its strong initial norm are imported from the finite-cell theorem and the collaborator's connecting-path extension. C47's finite-order SW identity alone is not an all-order matching theorem. An order-W uniform conserved-charge field implements the off-half-filled reference without closing the finite-cell gap at sufficiently small port coupling.

[scope] No nonlinear power law, thermal T^2 coefficient, full physical-sector selection, or central-charge conclusion is proved. The linear theorem is in the bulk zero-temperature shell representation; arbitrary finite-grid last-shell uniformity remains in Lemma N. Existing C56–C58 work and the collaborator's counterterm/source notes were read. No delegated agents were used.

[numerical] `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B checks.py` passes. Shell cancellation, determinant interpolation, history sums, relevant-tail recursion and a synthetic root bracket are checked. These are lemma controls, not an interacting RG calculation. Results are in [checks.json](checks.json).

[record] All writes stayed in work_C55. [STATE.md](STATE.md) is the local handoff because the project STATE is outside the permitted write scope. The proof has not received independent review.
