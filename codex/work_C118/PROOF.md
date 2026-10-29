# C118. Partition-function ratio estimates and the remaining gap

2026-10-29. Natural logarithms. All results below are finite-dimensional statements unless infrared input is explicitly stated. Work is confined to this folder.

## 1. Scope and the necessary zero-temperature bound

[Assumed] Use the setting, optimized bond scale M, and vacuum-compatible ring lengths of specs/_bound_background.md and C116. The Hamiltonian is H_N = sum_i h_(i,i+1), the local dimension is d, and Z_N(beta) = Tr exp(-beta H_N). Write

    g_N(beta) = log Z_N(beta),
    R_L(beta) = 2g_L(beta) - g_(2L)(beta),
    Delta_L = E_(2L) - 2E_L,
    W_N(beta) = Tr exp[-beta(H_N-E_N)].

No reflection symmetry, conserved charge, or stoquastic basis is assumed in the general target. Scalar bond shifts and coboundaries leave R and Delta invariant. C116 already proves |R| <= 4 beta M and the infrared extraction of A. Those settled results are used here.

[Proved] The identity

    R_L(beta) = beta Delta_L + 2 log W_L(beta) - log W_(2L)(beta)       (1)

implies that the proposed bound with C = k f(d),

    R_L(beta) <= K(d) + C M beta/L,                                  (2)

would force

    L Delta_L <= C M                                                  (3)

at every allowed L. Indeed, W_N(beta) tends to the ground degeneracy at fixed N, so dividing (1) by beta and sending beta to infinity gives Delta_L. K(d)/beta vanishes. This necessary condition is stronger than bounding the limiting Casimir coefficient. It cannot be bypassed by bounding an entropy term alone.

[Proved conditional implication] Conversely, a bound (3) together with

    2 log W_L(beta) - log W_(2L)(beta) <= K(d)                         (4)

would prove (2). Neither (3) nor (4) is proved here for arbitrary quantum chains. Equation (4) is only a sufficient additional condition. It is not claimed to follow from (2), since spare linear slope in (2) could absorb a growing thermal remainder.

[Assumed and proved implication] Under the vacuum Casimir expansion recorded in C116,

    Delta_L = pi A/(4L) + o(1/L),

(3) implies A <= (4k/pi) M f(d), without assuming convergence of full thermal traces. C116's stronger fixed-aspect-ratio trace hypothesis gives the same conclusion by taking L to infinity at beta = tau L, then tau to infinity. The two orders of limits should not be conflated.

## 2. A complete all-temperature estimate for XX

[Proved] For the periodic spin chain

    h = Sx tensor Sx + Sy tensor Sy,
    d = 2,  M = 1/2,

and every even L >= 4 and beta >= 0,

    R_L(beta) <= 2 log 2 + (5 pi/3) beta/L.                           (5)

Thus the requested form holds in this class with K(2) = 2 log 2, f(2) = 1 and k = 10 pi/3. Multiplying h by J replaces beta by |J| beta, and M by |J|/2. For negative J an alternating on-site rotation gives this equivalence on even rings. This coefficient is loose. The stronger free-chain bound on A from C101 remains the appropriate bound on the infrared coefficient itself.

### 2.1 Full spin trace and the parity factors

[Proved, Jordan-Wigner calculation] At fixed fermion-number parity, the hopping dispersion is cos k. Even fermion number uses k_j = 2 pi(j+1/2)/L, and odd fermion number uses k_j = 2 pi j/L. Let s be 0 or 1/2 and define

    P_s = product_(j=0)^(L-1) [1 + exp(-beta cos k_j)],
    Q_s = product_(j=0)^(L-1) tanh(beta cos k_j/2).

Expanding these products in occupation subsets, the sum with a plus sign counts all subsets, while the product of 1-exp(-beta cos k_j) weights each subset by its parity. Consequently the even and odd traces on a fixed grid are P_s(1+Q_s)/2 and P_s(1-Q_s)/2. The physical spin trace is

    Z_L = [P_(1/2)(1+Q_(1/2)) + P_0(1-Q_0)]/2.                     (6)

This expression retains every occupation sector. It is not a grand-canonical replacement for the spin trace.

[Proved] If L is divisible by four, the periodic grid contains zero modes, so Q_0 = 0. The antiperiodic grid has L/2 negative energies and no zero modes. L/2 is even, so Q_(1/2) >= 0. If L is two modulo four, the antiperiodic grid contains zero modes, so Q_(1/2) = 0. The periodic grid has L/2 negative energies, an odd number, so Q_0 <= 0. In both cases the nonzero parity correction in (6) has the enhancing sign and magnitude at most one. Therefore

    (P_0 + P_(1/2))/2 <= Z_L <= P_0 + P_(1/2).                     (7)

At beta = 0, both Q's vanish and (6) gives Z_L = 2^L. These bounds also hold at beta = 0 by continuity.

### 2.2 A uniform shifted-grid quadrature bound

[Proved] If f is a twice continuously differentiable 2 pi-periodic function and f_hat_n its Fourier coefficients, two integrations by parts give, for n != 0,

    |f_hat_n| <= (1/(2 pi n^2)) integral_0^(2 pi) |f''(k)| dk.

The absolutely convergent Fourier series can be summed on any shifted L-point grid. All frequencies except multiples of L cancel. Thus

    |sum_j f(2 pi(j+s)/L) - L f_hat_0|
      <= L sum_(m != 0) |f_hat_(mL)|
      <= [pi/(6L)] integral_0^(2 pi) |f''(k)| dk.                  (8)

Here sum_(m != 0) 1/m^2 = pi^2/3. The estimate is uniform in the shift s and needs no asymptotic expansion in beta/L.

[Proved] Apply (8) to

    f_beta(k) = log[1 + exp(-beta cos k)],
    n_beta(k) = 1/[1+exp(beta cos k)],
    p(beta) = (1/(2 pi)) integral_0^(2 pi) f_beta(k) dk.

Direct differentiation gives

    f'_beta(k) = beta sin k n_beta(k),
    f''_beta(k) = beta cos k n_beta(k)
                 + beta^2 sin^2 k n_beta(k)[1-n_beta(k)].          (9)

Since n_beta(k+pi) = 1-n_beta(k),

    integral |cos k| n_beta(k) dk = 2.

For the other term, set u = cos k on each half-period and let n(u) = 1/(1+exp(beta u)). Then

    integral sin^2 k n_beta(k)[1-n_beta(k)] dk
      = 2 integral_(-1)^1 sqrt(1-u^2) n(u)[1-n(u)] du
      <= 2 integral_(-1)^1 n(u)[1-n(u)] du
      = (2/beta)[n(-1)-n(1)]
      <= 2/beta.                                                   (10)

The expression has its continuous interpretation at beta = 0. By the triangle inequality in (9),

    integral |f''_beta(k)| dk <= 4 beta.

It follows that both grids obey

    |log P_s - L p(beta)| <= a_L,
    a_L = 2 pi beta/(3L).                                          (11)

### 2.3 Bulk cancellation

[Proved] Insert (11) into (7). The lower bound is

    Z_L >= exp[L p(beta)-a_L],

and the upper bound is

    Z_L <= 2 exp[L p(beta)+a_L].

Hence the finite-size correction Xi_L = g_L-Lp satisfies

    -a_L <= Xi_L <= log 2 + a_L.                                  (12)

Use the upper bound for both copies of L and the lower bound for 2L, at the same beta. The pressure cancels, giving

    R_L = 2 Xi_L-Xi_(2L)
        <= 2 log 2 + 2a_L + a_(2L)
        = 2 log 2 + (5 pi/3) beta/L,

as asserted. The step which supplies 1/L is (8), a controlled momentum quadrature error. It is not a spatial Schmidt-rank estimate.

[Proved checks] beta has inverse-energy units, so after rescaling J the right side depends on |J| beta/L = 2 M beta/L. At beta = 0, R = 0. At fixed L and beta tending to infinity, (5) bounds Delta_L. In the infrared it gives A <= 40M/3, consistent with A = 2M for XX but much weaker. No factor of two from a parity projection is missing, because both grids were included in (6).

[Numerical] checks.py verifies (6) by a stable even/odd occupation recursion against full exact diagonalization for L = 4,6,8. The maximum log-trace discrepancy is 2.14e-14. It tests (5) at L = 4,6,8,16,64,256,1024,4096 and ten aspect ratios from 0.01 through 32, including doubled length 8192. It also tests (12) against independent quadrature for sixteen choices of L,beta. Every check passes. The proof does not rely on these samples.

[Proved extension] For a tensor product of m independent XX chains with couplings J_a on each site, log Z and R add. Thus

    R_L <= 2m log 2 + (5 pi/3) beta sum_a |J_a|/L.

For this subclass M = sum_a |J_a|/2. The upper bound follows from the norm of the sum of centered independent bonds. The matching lower bound follows from tensor products of the bond's extremal Bell eigenvectors, all of whose one-site marginals are maximally mixed, so a coboundary changes neither extremal expectation. Here d = 2^m, and (5) extends with K(d) = 2 log d and the same k = 10 pi/3, f(d) = 1. This includes arbitrarily unequal velocities in this subclass. It does not extend the argument to interacting sectors or to the fixed-d constructions with unbounded sector count.

## 3. A complete classical comparison without transfer positivity

[Proved] Suppose h is diagonal in a fixed on-site product basis, h|a,b> = epsilon_ab|a,b>. It may be oriented, with epsilon_ab != epsilon_ba. Define the entrywise positive d by d matrix

    T_ab = exp(-beta epsilon_ab).

The classical ring trace is Z_L = Tr T^L. Let B = T^L. B is entrywise nonnegative even if it is neither Hermitian nor positive semidefinite. Therefore

    (Tr B)^2 <= d sum_a B_aa^2
             <= d sum_(a,b) B_ab B_ba
             = d Tr B^2.

Consequently

    R_L <= log d                                                   (13)

for every L and beta. This is (2) with K(d) = log d and zero slope. Matrix positive semidefiniteness, which C116 used for one possible transfer argument, is unnecessary in this classical subclass. Entrywise nonnegativity is enough. A general quantum bond does not supply a fixed-d spatial transfer matrix of this kind.

[Numerical] An oriented three-state classical bond was checked at four lengths and four temperatures. All sixteen tests of (13) pass. This is a calibration of the elementary matrix proof, not evidence for a quantum extension.

## 4. Seam averaging gives an exact comparison, but no small norm

[Proved] Label the two L-rings by a_i and b_i, with i modulo L. On their common Hilbert space let

    H_0 = sum_i [h_(a_i,a_(i+1)) + h_(b_i,b_(i+1))],
    V_i = h_(a_i,b_(i+1)) + h_(b_i,a_(i+1))
          - h_(a_i,a_(i+1)) - h_(b_i,b_(i+1)).                    (14)

H_0+V_i is one oriented ring of length 2L. All i give unitarily equivalent Hamiltonians. Set Vbar = (1/L) sum_i V_i. Convexity of log Tr exp gives

    g(H_0+Vbar) <= (1/L) sum_i g(H_0+V_i) = g_(2L),
    R_L <= g(H_0)-g(H_0+Vbar),                                    (15)

where g(H) = log Tr exp(-beta H). This is a valid coupling-convexity comparison with the correct direction.

[Proved obstruction to a norm estimate] The 1/L in Vbar does not imply ||Vbar|| = O(M/L). For h = J Z tensor Z, evaluate (14) in the product state with every a-spin up and every b-spin down. Every original bond has expectation J and every crossed bond -J. Hence <V_i> = -4J for all i and

    ||Vbar|| >= 4|J| = 4M.

The reverse bound follows from the four bond terms, so equality holds. This is true for every L. Averaging positions does not make the seam operators cancel on that state.

[Proved obstruction in the requested critical models] For XX, take every a-spin in the +x state and every b-spin in the -x state. Each original bond has expectation 1/4 and each crossed bond -1/4, so <Vbar> = -1 and ||Vbar|| >= 1 = 2M. The same states give the same value for Heisenberg since their y and z expectations vanish. For SU(3), use orthogonal colors |0> on all a sites and |1> on all b sites. Original swap bonds have expectation one and crossed swap bonds zero, so <Vbar> = -2 and ||Vbar|| >= 2 = 2M. These statements hold for all L, including vacuum-compatible lengths. Thus the attempted bound ||Vbar|| <= const M/L is false on the calibration models themselves.

[Proved scope] This defeats the norm estimate after the specific averaging (14). It does not prove that (15) must lose a factor L for every model, because its Gibbs expectation may be much smaller than its norm. Bounding that expectation sharply is additional work, not a consequence of averaging.

[Numerical] For XX and Heisenberg at L = 4 and SU(3) at L = 3, the averaged seam norms are 1.172603940, 1.322875656 and 2.774851773. The corresponding zero-temperature upper bounds in (15), divided by beta, are 0.404587196, 0.770621461 and 2.464816242. The actual Delta_L values are 0.215301195, 0.348906591 and 1.394448725. These small rings demonstrate slack, not an asymptotic lower bound on the slack.

### 4.1 A constant partial swap is not a symmetry

[Proved necessary condition for a different spreading construction] A common twist-spreading argument starts with a continuous on-site symmetry, then makes its rotation angle vary by O(1/L) between neighboring sites. The two-copy Hamiltonian has the discrete symmetry that swaps every a site with its b partner. That does not imply a continuous symmetry generated by the sum of on-site swaps.

On two neighboring sites define S_i as the swap of a_i and b_i, and

    h_0 = h_(a_i,a_j) + h_(b_i,b_j),
    G = S_i+S_j.

If [h_0,G] != 0, then exp(i theta G) does not preserve h_0 even at a constant angle. Any operator estimate which bounds the entire resulting bond change solely by a power of theta_i-theta_j would give zero at constant angle and would therefore be false. A bulk change must be accounted for.

[Proved explicit witness] Order the four factors as (a_i,a_j,b_i,b_j). For XX,

    H_0 |1010> = (|0110>+|1001>)/2,
    G |1010> = 2 |1010>,
    G H_0 |1010> = |1100>+|0011>.

Thus <1100|[H_0,G]|1010> = -1. For Heisenberg the extra diagonal term on |1010> is -|1010>/2 and cancels out of this matrix element, so the same value -1 holds. For SU(3), use colors 0 and 1. The two swaps in H_0 have amplitude one instead of 1/2, giving the value -2. These nonzero matrix elements establish the absence of the proposed continuous symmetry in all three models.

[Numerical] The commutator norm and the constant-angle mismatch at theta = pi/4 are evaluated for XX, Heisenberg and SU(3), and are nonzero in each case. This check rejects that continuous-symmetry premise for the natural partial-swap interpolation. It does not rule out another auxiliary symmetry or a different gluing construction.

### 4.2 Rank does not control the insertion's condition number

[Proved] A rank-one insertion can have arbitrarily large trace effect. For Q = c |0><0| and X = |0><0|, Tr(XQ)/Tr X = c for every c > 0, although rank Q = 1. Even restricting 0 <= Q <= I only gives an upper bound on that ratio, while its inverse is unbounded as c tends to zero. A two-sided seam comparison needs an appropriate norm and lower-overlap bound, not rank alone.

[Proved limitation inherited from C116] A Trotter spatial cut crosses the entire imaginary-time history. The rank at one time slice is at most d^2, but the history rank grows with the number of slices. A factor exp(const beta M log d) still loses L at beta = tau L. Replacing beta M by beta M/L requires a new estimate on propagation or compression across a spatial width L. Neither the local Schmidt rank nor (14) supplies it. An error bound for imaginary-time evolution also cannot be inferred directly from a real-time Lieb-Robinson bound.

## 5. Convexity and a single dimensionless temperature

[Proved] For a finite spectrum,

    g'_N(beta) = -<H_N>_beta,
    g''_N(beta) = <(H_N-<H_N>)^2>_beta >= 0.

For the ratio,

    R''_L(beta) = 2 Var_(L,beta)(H_L) - Var_(2L,beta)(H_(2L)).       (16)

Both variances are nonnegative, but their difference has no determined sign. Convexity in a coupling constant has the same issue when taking a signed combination of log traces. Also, an upper bound on a convex function at one point does not upper-bound its slope at larger arguments.

[Proved explicit failure of both global signs] For XX at L = 4, the full spin spectrum has a unique ground energy -sqrt(2), and first gap

    delta_4 = sqrt(2)-1.

For L = 8, the unique ground energy is -csc(pi/8), and the first gap is

    delta_8 = csc(pi/8)-1-sqrt(2).

Each first excited level has multiplicity two. These statements follow by filling the periodic odd and antiperiodic even grids from Section 2.1; the finite enumeration is also checked by the script. Since 0 < delta_8 < delta_4, the low-temperature expansion of (16) is

    R''_4(beta) = -2 delta_8^2 exp(-beta delta_8)
                 + o(exp(-beta delta_8)) < 0                     (17)

for sufficiently large beta. On the other hand R'_4(0) = 0, because the infinite-temperature energies are extensive, while

    R'_4(infinity) = 2 sqrt(2)-csc(pi/8) > 0.

A differentiable globally concave function could not have its derivative increase from zero to that positive limit. Thus R_4 is neither globally convex nor globally concave. This is an analytic obstruction on a vacuum-compatible ring, not a twist artifact.

[Numerical] On 121 aspect ratios from 0.025 to 32, (16) has both signs for XX L = 4, Heisenberg L = 4 and SU(3) L = 3. The observed minimum and maximum values are respectively (-0.0243785,0.1574501), (-0.0268226,0.2832090), and (-0.5142742,2.8067394). R_L-beta Delta_L is nonpositive on all these samples. That last observation is not a proof of (4), and would still leave (3) if established.

### 5.1 The CFT scaling form alone cannot make the reduction

[Proved within the stated continuum assumptions] Choose any fixed compact unitary modular-invariant c = 1 CFT with unique vacuum and smallest positive scaling dimension x_* > 0. Let F(y) be its torus log trace, with the convention of C116. Vacuum dominance and modular invariance give

    F(y) = pi/(6y) + O(exp(-2 pi x_*/y))   as y -> 0.

For n independent copies, all with velocity v_n = n^(-1/2), the equal-temperature scaling ratio is

    r_n(tau) = n [2 F(v_n tau)-F(v_n tau/2)].

The 1/y terms cancel, so at every fixed tau_0 > 0,

    |r_n(tau_0)| = O(n exp[-2 pi x_* sqrt(n)/tau_0]) -> 0.          (18)

But the large-tau slope is

    lim_(tau -> infinity) r_n(tau)/tau = (pi/4) n v_n
                                     = (pi/4) sqrt(n) -> infinity. (19)

Every member of this continuum family has finitely many sectors and finite torus traces, and the velocities have a common upper bound. Thus a bound at one tau_0, the exact CFT scaling form, and a per-sector velocity ceiling do not by themselves bound A. Information linking the microscopic d,M to the allowed collection of sectors is indispensable.

[Proved scope] Equations (18)-(19) are not a fixed-d lattice counterexample. No fixed-d microscopic realization with these specified velocities is asserted. They refute only a proposed deduction using continuum scaling data and a single-temperature bound without further microscopic input. The known slow-sector constructions make discarding this distinction especially problematic.

## 6. Exact relative-entropy identities

[Proved] Use H_0 and H_1 = H_0+V from (14), and define their normalized Gibbs states

    sigma = exp(-beta H_0)/Z_L^2,
    rho = exp(-beta H_1)/Z_(2L).

Both are full rank at finite beta. Their logarithms are known explicitly. Substituting them in D(a||b) = Tr a(log a-log b) gives

    D(sigma||rho) = beta <V>_sigma - R_L,
    D(rho||sigma) = R_L - beta <V>_rho.                            (20)

Equivalently,

    R_L = beta <V>_sigma - D(sigma||rho),
    R_L = beta <V>_rho + D(rho||sigma).                            (21)

These signs matter. Dropping the first relative entropy gives an upper bound by the seam energy in the split state. Dropping the second gives a lower bound, not the desired upper bound.

[Proved] Let A and B be the two halves in the H_0 tensor factorization, and let gamma_L be the ring Gibbs state. Since sigma = gamma_L tensor gamma_L,

    D(rho||sigma) = I(A:B)_rho
                   + D(rho_A||gamma_L) + D(rho_B||gamma_L).        (22)

To verify (22), use log sigma = log gamma_L tensor I + I tensor log gamma_L, and add and subtract S(rho_A)+S(rho_B). This yields

    R_L = beta <V>_rho + I(A:B)_rho
           + D(rho_A||gamma_L) + D(rho_B||gamma_L).                (23)

The two halves are equivalent by translation around the doubled ring, but their marginals need not equal the thermal state of a closed L-ring. Dropping those two positive marginal relative entropies cannot give an upper bound on R.

[Proved] Finite dimension bounds mutual information, I(A:B) <= 2 min(log dim A,log dim B), but not relative entropy. For instance, on one qubit let H = Delta |1><1| and gamma_beta = exp(-beta H)/Tr exp(-beta H). Then

    D(|1><1| || gamma_beta) = beta Delta + log(1+exp(-beta Delta)),

which is unbounded as beta tends to infinity although log dim = log 2. Thus a seam region containing w sites has an information-capacity ceiling of order w log d for its entropy or mutual information, but no such ceiling for its relative entropy to a cold reference state.

[Proved] For any fixed finite ring with Delta_L > 0, (1) makes R_L grow linearly in beta, while every mutual information between its fixed finite subsystems stays bounded. Therefore R_L <= I(A:B), or R_L <= I(A:B)+a fixed seam-capacity constant, cannot hold at all beta. Allowing a beta/L term avoids this contradiction, but its coefficient is then the quantity still to be bounded.

### 6.1 The thermal mutual-information area law retains beta M

[Proved] Write an arbitrary bipartition Hamiltonian as H_A+H_B+H_boundary. The Gibbs state rho minimizes Tr rho H - S(rho)/beta. Compare it with rho_A tensor rho_B, which has the same separate-region energies. The variational inequality gives

    I(A:B)_rho <= beta Tr[H_boundary(rho_A tensor rho_B-rho)]
                <= 2 beta ||H_boundary||.                       (24)

For the two halves of a ring there are two crossing bonds, so in an optimized gauge ||H_boundary|| <= 2M and

    I(A:B)_rho <= 4 beta M.

The dimensional ceiling 2L log d is also valid. Neither bound has the needed beta/L dependence. A strip of width proportional to L gives a dimensional ceiling proportional to L, and taking L to infinity at fixed tau does not remove it.

[Assumed published context] This is the elementary Gibbs-variational mechanism behind the thermal mutual-information area law in Wolf, Verstraete, Hastings and Cirac, “Area laws in quantum systems: mutual information and correlations,” Phys. Rev. Lett. 100, 070502 (2008), https://arxiv.org/abs/0704.3906. The primary abstract explicitly distinguishes the temperature-independent classical prefactor from the quantum prefactor which may grow with inverse temperature. The finite-dimensional derivation (24) is supplied above rather than relying on a quoted theorem to infer a stronger bound.

[Numerical] Direct density-matrix calculations for XX and Heisenberg on eight total sites, and SU(3) on six total sites, verify (20)-(23) and Jensen (15). The largest identity residual is 2.60e-12. At beta = 16, the values (R_L, I(A:B), D(rho_A||gamma_L)+D(rho_B||gamma_L)) are

    XX, L = 4:          (3.370382150, 1.832322103, 5.657854099),
    Heisenberg, L = 4:  (5.581805924, 2.096886450, 12.043550829),
    SU(3), L = 3:       (22.311179590, 3.322852827, 40.864071625).

In each case R exceeds the mutual information. The seam term beta<V>_rho is negative and cancels a substantial part of the relative entropy. Treating the relative entropy as a capacity-limited correction misses this cancellation.

## 7. What is established and what remains open

[Proved] There are explicit all-temperature ratio bounds for the XX chain, its independent tensor products, and all diagonal classical bonds. The XX proof gets beta/L from a momentum quadrature estimate with integrated second derivative bounded by beta. The classical proof gets a constant from the finite entrywise nonnegative transfer matrix.

[Proved or numerical with tags above] The direct extensions examined here fail at identified steps. Uniform seam averaging does not reduce the operator norm. A discrete replica swap does not automatically generate a continuous symmetry. Rank does not bound the condition number of an insertion. R is not convex or concave even on XX. Continuum scaling plus a bound at one tau cannot account for slow sectors. Mutual information omits both the seam energy and the marginal relative entropies.

[Open] No finite K(d), f(d), and k valid for all chains in the question have been obtained. No counterexample to that estimate or to A <= C M f(d) has been constructed. A useful new quantum statement would have to control the full finite-size free-energy correction, or the cancellation in (21), with an error proportional to beta/L. The computations in this folder do not establish such control.

[Work record] All local checks complete in about two seconds with one BLAS thread. They include every eigenvalue of the small interacting rings and linear-cost free-fermion traces at large L. The provided SLURM script extends full-spectrum checks to Heisenberg length 12 and SU(3) length 12. The latter has a largest occupation block of dimension 34650 and requires substantial memory and runtime. It was not run locally or submitted to Typhon. No network shell command, external write, or independent agent/referee review was used.
