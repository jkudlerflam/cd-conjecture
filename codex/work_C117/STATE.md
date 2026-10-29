# C117 local state

2026-10-29.

[Open] General A<=C M f(d), without symmetry or crossover hypotheses, remains unproved. The three routes are analyzed in RESULT.md and derived in PROOF.md. Parent project memory was not edited because writes outside work_C117 are forbidden.

[Proved] Fixed-support translation-invariant marginal and local ground-state SDPs have an L-independent feasible set for L>2R. Shorter-ring ground states are feasible. For A>0 this enforces a fixed energy-density deficit and blocks Casimir resolution. Under an O(R) compatible-shorter-length assumption, R=o(L) cannot give o(1/L) energy error. This excludes only local fixed-support programs, not global or winding constraints.

[Proved] The local energy first moment is b1(cos q-1)+b2(cos 2q-1), with |b1|<=8M^3 and |b2|<=4M^3. Hence its q^4 coefficient has absolute value at most 3M^3. This coefficient is not purely infrared. The gapped optimized-gauge Ising bond -JXX-g(ZI+IZ)/2 has M=sqrt(g^2+J^2), A=0, and coefficient gJ^2/4+O(J^4/g). The perturbation derivation is in PROOF Section 6; moment_check.py checks normalization on six sites.

[Proved] Exact k=3 optima including all one-site ground positivity are -1/2 for Heisenberg, -(1+sqrt(33))/16 for XXZ at Delta=1/2, and -1 for SU3 swap. Certificates come from open three-site lower bounds and shorter-ring feasible marginals.

[Proved] All 15 local certificates and 10 full-ring intervals pass independent integer verification. bootstrap.py uses scipy dual search, with W=BB^T, and then rationalizes the coefficients. verify.py rebuilds the dual matrices using integer arithmetic. A quantized Cholesky factor plus an exact residual row bound certifies PSD. ring_checks.py adds exact rational Rayleigh upper bounds. Matrices are at most 27 for local SDPs and 81 for tiny full rings. Random d=3 k=3 searches hit iteration limits and are not claimed optimal.

[Assumed sufficient lemmas] (i) Filtered inverse moment F(q)=integral_0^(Lambda|q|) Sbar/omega: stress matching plus limsup F/q^2<=K M f implies A<=24pi K M f. (ii) Ring certificates ell_L<=E(L), u_2L>=E(2L), with limsup L(u_2L-2ell_L)<=K M f imply A<=4K M f/pi. (iii) The C106 central-charge-counting and logarithmic-support lemma, with gamma=2beta log d<1, gives A<=4e M d^(2k0)/[eta(1-gamma)]. None of the three missing inequalities is proved universally.

[Proved work record] Read the full campaign background and REPORT 9.68 to 9.78, plus C100/C103/C106 results. Consulted primary literature via the browsing tool. One Typhon SSH attempt exited 255 with "Connection timed out during banner exchange". No remote job ran. run_typhon.slurm is prepared for k up to 5 at d=2 and 4 at d=3, using the same one-site ground constraints. Extended mode refuses to run without a SLURM allocation. No install, delegation, commit, push, or edit outside this directory.

[Open next steps] Obtain Typhon access for larger blocks and genuinely circumference-dependent moment matrices. For analytic progress, prove the filtered inverse-moment budget or construct paired upper/lower ring certificates with the stated rate. Ground-localizing certificates are not automatically operator SOS and must be distinguished from marginal-only certificates in any handoff to C114.
