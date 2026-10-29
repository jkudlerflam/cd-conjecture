# C96 result

[Status] Partial resolution. The finite-N tagged contact kernel now has an exact determinant formula, and C94's first-moment hypothesis can be weakened. The uniform tagged-kernel bound and the full low-temperature relative remainder are still unproved. The d=13 mobile defect comparison and C93's contact-clock theorem remain open.

[Proved, exact tagged kernel] PROOF equations (4)-(14) give a complete free spectrum and every tagged contact form factor. A material Fourier mode theta=2 pi m/N requires physical Slater momenta

    p_a=(2 pi/L)(I_a+xi+m/N),
    sum_a I_a+N xi+m=0 modulo L,

where xi=0 or 1/2 according to particle parity. The form factor of A_0 is the coefficient of zw in an N by N determinant built from the physical positions 0,1 and 2,...,L-1. Thus

    chi_N(r)=sum_(m,I excited)
          exp(i theta_m r)|F_(m,I)|^2/(E_(m,I)-E_0),
    k_N(r)=chi_N(r)+c_dens/N.

The zero mode reduces to C95's two-particle-hole calculation. Nonzero material modes change the momentum grid. Restricting those sectors to two particle-hole excitations on the original grid is not valid. A local physical n_x n_(x+1) correlator is a different kernel.

[Open, requested asymptotic] This exact finite sum has not been bounded uniformly in N. No r^(-2) tail or bounded first absolute moment is claimed. At half density, direct checks through N=8 give first moments after the specified subtraction near 0.1102/t. Those sizes do not settle an asymptotic tail.

[Proved, replacement for C94] Define the discrete Laplacian Delta k(r)=2k(r)-k(r+1)-k(r-1) and

    M_b(k,c)=sum_r |k(r)|
           +b^2 sum_r(1+dist(r,0))|Delta k(r)|+|c|.

For the same zero-port product-cell coefficient as C94,

    |B_sigma-B_good|<=32 M^2 b W M_b(k,c).

The proof separates source means from within-cell covariances. The translated selected means have equal averages over complete cells. Their zero-mean part has a bounded discrete primitive, allowing two differences to be transferred to the bath kernel. Wrong cells require only their material length, regardless of their excitation gap.

A smooth nonoscillatory r^(-2) tail has finite M_b even though its first absolute moment diverges logarithmically. A resonant oscillatory tail can still produce logarithmic interface costs. The distinction is necessary. One may also apply the criterion directly with k=chi,c=0, since Delta kills the canonical constant. The actual tagged kernel has not yet been shown to satisfy this weaker uniform criterion.

[Proved, thermal quadratic extension] Retain J a_0 Hmat in the reference and interpolate only the centered contact interaction. The real part of a thermal material covariance satisfies

    0<=Re C_mat(tau)<=Re C_mat(0),    0<=tau<=beta.

The real free-bath covariance is positive semidefinite at each imaginary time. These facts bound the connected resummed quadratic coefficient by the static bath response contracted with the material equal-time covariance. For isolated product thermal cells, the same W-relative estimate holds, with the thermal static kernel. This keeps all states in each cell and is uniform in J/T and the wrong-cell gap. A uniform thermal bath norm and the nonzero-port matching remain missing.

[Proved, limitation of quadratic resummation] A new exact two-spin bath diagnostic has Hmat=0, so its frozen dynamics is already retained exactly. Its residual is

    R=-T log{[1+(cosh(sqrt(t^2+J^2/4)/T)/cosh(t/T))^2]/2}.

Nevertheless its fourth-order coefficient diverges as -1/(128t^2T). As T tends to zero, R=-J^2/(4t)+O(J^4/t^3), while its quadratic thermal coefficient gives -J^2/(8t). Replication preserves this discrepancy when volume tends to infinity first. The exact bound -J^2/(4t)<=R<=0 remains uniform. Bath-induced material interactions must be retained inside the thermal trace, or the full residual must be bounded directly. This diagnostic is not a C88 counterexample.

[Open, full thermal remainder] PROOF equation (42) states the exact remaining relative interpolation estimate. The free determinant formula and the product-cell quadratic bound do not control the interacting response along that interpolation. No uniform full-auxiliary, nonzero-port, wrong-length remainder is established here.

[Conditional, d=13] At b=15, a bound t M_15<=C would give |Delta B|<=480 C M^2 W/t. The full legal comparison still needs the relative thermal estimate. C92's activity sum would then yield O(T exp(-delta/T))=o(T^2), with volume first. Finite grammar penalty separately needs the relative penalty comparison. C93's dressed clock threshold and auxiliary heat bound do not follow from a static kernel estimate.

[Numerical and scope] The determinant formula agrees with independent gap resolvents to below 3e-15 through bath dimension 6435. The thermal covariance, weaker-norm, and exact diagnostic checks pass. An initially incorrect numerical threshold in the test-kernel check was corrected using its analytic bound, with the failure and correction recorded in CHECKS.md. All writes and Python execution stayed in work_C96. No network commands, delegated agents, large computations, or external edits were used.
