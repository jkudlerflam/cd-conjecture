# C66 result

[proved, outcome] The n=1 specialization is written explicitly in [PROOF.md](PROOF.md). It has a payload fermion, a length fermion, four chiral branches and six marginal density couplings. The action retains arbitrary-degree retarded kernels, quadratic counterterms, field normalizations, scalar terms and sources. The corrected Lemma N and C50's comparison are carried through as conditional statements. Their missing global routing and marked comparison estimates remain unproved, including at n=1.

## The action and its coefficients

[proved from the imported projected representation] Let m_j be the length occupation, and let Q_ab=lambda^2 K2_ab be the four physical projected payload bonds. With

    A=Q_00, B=Q_10-Q_00, C=Q_01-Q_00,
    D=Q_11-Q_10-Q_01+Q_00,

the static projected bond is

    -j_l(d_j^dagger d_(j+1)+h.c.)
       +A_j+m_j B_j+m_(j+1) C_j+m_j m_(j+1) D_j.

The scalar length interaction is V_l=tr(D)/4. The traceless B,C,D terms cannot be removed by a marker chemical potential or by replacing length occupations with their mean. They include six- and eight-field operators already at this projected order.

[proved, coefficient dictionary] PROOF Section 3 gives all six first-order localized coefficients directly in terms of the number-conserving two-qubit matrices A,B,C,D. It includes Gaussian contractions of their higher-field terms. The local action is

    g_p rho_pR rho_pL + g_l rho_lR rho_lL
        +sum_(sigma,omega=R,L) g_sigma,omega rho_psigma rho_lomega.

For a pure nearest-neighbor interaction V_a, the within-species coefficient is 4(V_a/W)sin^2(k_a). The full length coefficient also receives the payload expectation of the traceless part of D. The exact matched six coefficients require finite matching of the retarded action; the numerical c calibration does not supply them.

## The scheduled velocities

[proved from C39/C45's reference] In time units t'=Wt,

    W=lambda^2 J0/s,   J0=2g^2 q^4 D_aux^(-6)/(1+q^2+q^4)^2,
    t_p=sqrt(1+eta^2)/256,
    v_p=sqrt(4t_p^2-mu_p^2),
    j_l/W=1,   v_l=2sin(pi nu).

The four signed reference velocities are (v_p,-v_p,v_l,-v_l). Physical cell speeds have the extra factor W. These are the scheduled reference speeds, before finite-parameter and interaction corrections. For small eta and mu_p, v_p is near 1/128 and the length speed is of order one.

[proved, two-channel changes] Mixed couplings must be controlled in units sqrt(v_p v_l), and within-payload couplings in units v_p. Propagator and Gram bounds cost inverse minimum speed but no inverse difference of speeds. Coinciding speeds do not cause a new singularity. Separate charge conservation excludes payload-length conversion, while mixed exchange terms need the momentum margins. PROOF Section 3 also explains why gauging away the payload hopping phase must be distinguished from checking an unrestricted pair-sum hypothesis.

## Lemma N, comparison and attraction

[proved, conditional implications] The four-branch routing has the same Taylor orders and paired-cut identity as C61. A routing proof based only on uniform line/Gram bounds and kernel moments would extend to four branch labels with a fixed exponential coefficient loss, not a depth-dependent loss. That premise has not been proved. The corrected source statement keeps D_O(0) separate and assigns an interaction prefactor only to D_O(u)-D_O(0). Adding a spectator length channel does not remove C60's free pair descendant.

[proved conditional on the stated RG estimates] In a velocity-weighted norm, C50's comparison becomes

    B_h <= A_eff d_for 2^(theta h)
       +b eta sum_(k>h) 2^(-kappa(k-h)) sum_(j=h+1)^k B_j.

Here eta bounds the marginal history and d_for bounds comparison forcing. They are not automatically the same parameter or the square of one another. PROOF Section 7 controls the four relevant directions conditionally, derives this comparison, and proves its summable closure with C50's explicit geometric constant. C56's continuous counterterm construction then gives the four shifted roots at fixed physical Hamiltonian. No extra microscopic chemical-potential tunings are required.

[proved, attraction constraint] For the isolated attractive length t-V chain, V_l<=-2j_l gives linear thermodynamic energy at fixed interior density, by a positive two-site decomposition and a droplet upper bound. Thus the length-fluid reference needs |V_l|<2j_l. C44/C45's complete equalization gives |V_l|<=4W d_eq/J0 with j_l=W. Lowering lambda alone does not improve this ratio. The coupled Gaussian action additionally needs its full positive stiffness reserve; the mixed-block Schur complement shows why the length-only condition is insufficient for the two-channel problem.

[numerical, imported evidence and limitation] verify_C44 finds calibrated c near 2 for its projected M1 point with V_l/|t| about -0.32, and c near 1 for C1 with V_l/|t| about -6.3. Both use j_l/|t|=1, rather than C45's j_l/(W t_p)=1/t_p. These tests support the projected two-mode regime and its collapse control, but do not test the scheduled velocity hierarchy or prove the retarded-action comparison.

## Checks and remaining work

[numerical] `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B checks.py` passes. The 18 six-coupling checks have maximum residual 1.36e-14. Bond reconstruction, slopes, the attraction matrix, coincident-speed Gaussian stability, comparison summation, relevant recurrence and free-source subtraction also pass. [CHECKS.md](CHECKS.md) records the results. These are finite algebra checks, not microscopic RG calculations.

[assumed, remaining] Prove the simultaneous routing and normalized marked comparison for the actual kernels, then the critical/source estimates. The scalar thermal estimate remains separate. No microscopic c=2, complete compact-boson limit or finite-ring spectrum is established here. All writes stayed in work_C66; local [STATE.md](STATE.md) supplies the handoff because the parent state is outside the permitted write scope.
