# C114. A ring-only bound on the Casimir coefficient

2026-10-29. Full derivations are in [PROOF.md](PROOF.md). All work is confined to this folder.

[Open] I did not prove A <= C M f(d) for general nearest-neighbor chains, for any explicit finite C and f, and did not find a counterexample. The new results are an exact coherent-reconnection inequality, an analytic obstruction to that entire trial family on XX, and a quantitative formulation of the missing positive-certificate estimate.

## A symmetry-free coherent-reconnection inequality

[Proved] Let psi_L be a translation eigenvector in the ring ground space. Write e_L=E_L/L, let sigma_L be its one-site marginal, and let rho_r be its contiguous r-site marginal. Define

    q_L = 2[tr(h (sigma_L tensor sigma_L))-e_L],
    S_L = sum_(r=0)^(L-1) tr(rho_r^2),   tr(rho_0^2)=1.

Then

    E_(2L)-2E_L <= q_L/S_L.                          (1)

This improves bare gluing by a factor S_L and requires no on-site symmetry. The trial coherently superposes all prefix-swap reconnections of two copies of psi_L. Its Gram matrix consists of interval purities. The real part of its energy-defect matrix is q_L times the identity, which makes the uniform coherent sum explicit.

[Proved] If h and psi_L are real and q_L>0, equation (1) is the minimum over the whole span of these reconnection states, including arbitrary complex coefficients. Reweighting the seam positions cannot improve it.

[Proved] On XX, this optimized cost is not O(1/L). An elementary number-variance estimate gives

    tr(rho_r^2) <= min(r,L-r)^(-alpha),   alpha=2/pi^2,
    q_L = 2/[L sin(pi/L)] >= 2/pi,
    S_L <= 3+[2/(1-alpha)](L/2)^(1-alpha).

Consequently

    L q_L/S_L -> infinity.

This is an analytic obstruction, independent of a CFT entropy assumption. The true XX doubling defect is O(1/L). The result excludes this coherent reconnection subspace, not arbitrary adapted local unitaries and not the desired Casimir bound.

[Assumed CFT calibration; numerical consistency] Standard interval-purity scaling predicts a cost proportional to L^(-3/4) in the c=1 examples and L^(-1/2) for SU(3). The finite XX correlation-matrix checks support the first power. Neither prediction is needed for the analytic obstruction above.

## The other two proposals

[Proved] An energy-preserving lift to a kL ring would give E_(kL)<=k E_L. For positive A the vacuum expansion gives the opposite strict sign,

    E_(kL)-k E_L = (pi A/6L)(k-1/k)+o(L^-1)>0.

XX supplies an exact finite-ring counterexample. A product of copies is a valid state, but each new seam replaces a correlated closing bond by a product marginal and costs O(M).

[Proved] A convex combination of larger-ring energy densities requires an error

    (pi A/6L^2)[1-sum_i w_i/k_i^2]+o(L^-2).

Thus the necessary quantitative correction is already proportional to A. Neither convexity nor subadditivity bounds its coefficient. The scalar-energy obstruction established in C107 still applies.

[Proved] An exact finite-range identity h-e0=P+tau(B)-B with P>=0 would wrap with no residual seam and imply E_L>=L e0. Adding a finite local Hamiltonian commutator does not help, since its ring-ground-state expectation vanishes. Such an identity is incompatible with a positive Casimir coefficient.

[Proved conditional reduction; open hypothesis] A support-r identity with positive part, coboundary, Hamiltonian commutator, and residual norm

    ||R_r|| <= B M f(d)/r^2+o(r^-2)

would imply

    A <= (24B/pi) M f(d).

No construction with this error rate and coefficient was obtained. Exact identities and o(r^-2) residuals are excluded when A>0. Averaging an O(M) seam over its positions does not reduce its total expectation in a translation eigenstate.

## Verification and limitations

[Numerical] Small full-ring diagonalization tested XX, Heisenberg, XXZ at Delta=-1/2 and +1/2, and SU(3). The qubit comparisons use L=4 and 8, and SU(3) uses L=3 and 6. All are vacuum-compatible. The zero-cost lift and zero-error convex comparison fail; the coherent inequality and its real-state optimum pass. The largest Hilbert space has dimension 729 and the maximum eigenpair residual is 1.33e-13. See [CHECKS.md](CHECKS.md).

[Numerical] Additional single-ring ED evaluates the exact coherent cost through L=8 for qubits and L=6 for SU(3), still with dimension at most 729. See [SINGLE_RING_CHECKS.md](SINGLE_RING_CHECKS.md). The scaled cost L q_L/(M S_L) increases across these tested sizes in every family. No asymptotic interacting fit is claimed.

[Numerical] XX correlation matrices through L=64, of matrix dimension at most 32, give L q_L/S_L increasing from 1.103776 to 1.965547. Three arbitrary complex d=2 bonds also pass the symmetry-free identities. See [PURITY_CHECKS.md](PURITY_CHECKS.md).

[Operational limitation] Three Typhon SSH attempts, separated by analytical and local checking work, each exited 255 with 'Connection timed out during banner exchange'. No job was submitted and no computation ran on the login node. The expanded ED suite is prepared in [run_ed.slurm](run_ed.slurm) and guarded against execution outside a SLURM allocation. It remains unrun. No heavy computation ran locally.

[Open next step] A successful ring argument needs either (i) a positive certificate with the quantitative r^-2 residual above, or (ii) a genuinely adapted reconnection beyond coherent translations of bare seams. Neither estimate has been proved here. The universal bound remains open.
