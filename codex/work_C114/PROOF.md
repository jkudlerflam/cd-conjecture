# C114. Ring energies and coherent reconnection

2026-10-29. The universal bound remains open. Tags apply to the arguments following them. This work has not received independent review.

## 1. Scope and inherited input

[Assumed] Let H_L=sum_j h_(j,j+1) be a one-site translation-invariant nearest-neighbor spin ring. The local Hilbert space has dimension d. On a specified vacuum-compatible sequence, including its fixed integer multiples, assume

    E_L = L e0 - K/L + o(L^-1),
    K = pi A/6,             A = sum_a c_a v_a >= 0.

The coefficient concerns the fixed Hamiltonian followed by L tending to infinity. No uniform onset length across Hamiltonians is assumed. The target is A <= C M f(d), where

    M = inf_(alpha,o) ||h + alpha I + o tensor I - I tensor o||.

[Proved, inherited from C98] Gauge shifts cancel from E_L-L e0 and from topology changes. One can work with a representative of norm at most M+epsilon and then let epsilon decrease to zero. Define

    D_L = E_(2L)-2E_L = pi A/(4L)+o(L^-1).             (1)

C98 already proves |E_L-L e0|<=2M, |D_L|<=4M, and the dyadic identity

    L e0-E_L = sum_(j>=0) 2^(-j-1) D_(2^j L).         (2)

Those are inputs here, not new results. C107 already supplies countermodels to deriving a coefficient bound from scalar almost additivity, concavity, and finite-size continuity. Its countermodels are not local Hamiltonians. Neither their construction nor the earlier gradient and modulation arguments is repeated.

[Assumed] A decoupled massive sector with ring correction bounded by C L^p exp(-L/xi), at fixed finite xi, has zero contribution to lim L D_L. This follows from L^(p+1) exp(-L/xi) tending to zero. The arguments below use ring energies and trial states on rings, not local deformations intended to isolate the anomaly. Suppression is taken after fixing the massive sector, not uniformly as its gap closes.

## 2. An exact lift has the wrong consequence

[Proved] Suppose an L-ring ground state could be lifted to a state sigma on the kL ring with energy k E_L. The variational principle would imply

    E_(kL) <= tr(H_(kL) sigma) = k E_L.                (3)

But the assumed vacuum expansion gives, for fixed integer k>1,

    E_(kL)-k E_L = K(k-1/k)/L + o(L^-1).              (4)

If A>0, the right side is positive for sufficiently large compatible L. Thus the proposed energy-preserving lift does not exist. This contradicts the premise of route (i), not the desired inequality.

[Proved explicit example] In the inherited XX normalization h=SxSx+SySy, even rings have E_L=-csc(pi/L). For even L>=4, put x=pi/(2L). Then

    D_L = 2 csc(2x)-csc(x)
        = [sec(x)-1]/sin(x) > 0.                     (5)

This is an exact finite-ring counterexample to (3). No CFT limit is needed.

[Proved] A product of k copies is a legitimate lift as a state, but not an energy-preserving one. Choose a translation eigenvector psi_L in the ground space, let sigma_L be its one-site marginal, and set

    e_L=E_L/L,
    b_L=tr[h(sigma_L tensor sigma_L)].

Replacing the k internal closing bonds by k intercopy bonds gives

    trial energy = k E_L + k(b_L-e_L),
    b_L-e_L <= 2M.                                   (6)

Every copy retains L-1 bonds with their old energy. The new bonds see independent one-site marginals. A spatially periodic quantum state need not preserve the closing-bond correlations of each finite ring. If an entire L-site block marginal were the pure state psi_L, that block would factor from its complement, which makes the loss of closing correlations explicit.

## 3. Convex comparisons of longer rings

[Proved] Write e_n=E_n/n. Let w_i>=0, sum_i w_i=1, with finitely many fixed k_i>=2. The infinite-chain energy density may also be included by setting k_i=infinity. The vacuum expansion gives

    sum_i w_i e_(k_i L)-e_L
      = K[1-sum_i w_i/k_i^2]/L^2 + o(L^-2).           (7)

The coefficient in brackets is at least 3/4. Therefore the zero-error lower comparison e_L>=sum_i w_i e_(k_i L) fails for every positive A. With an error B M f(d)/L^2, equation (7) instead requires

    K[1-sum_i w_i/k_i^2] <= B M f(d).                (8)

An independently proved comparison of this strength would solve the problem. Convexity and the exact-lift premise do not supply it. The same observation applies to a comparison with e0 alone.

[Proved] Averaging bare seam positions as a density matrix cannot reduce the energy cost. Translation invariance makes each translated trial have the same energy, and energy is linear in the density matrix. Coherent superposition is a different operation and is analyzed in Section 5.

## 4. Finite-range positive certificates cannot be exact

[Proved] Suppose there were an identity in the infinite-chain local algebra

    h-e0 I = P + tau(B)-B,              P >= 0,        (9)

where P and B have finite support and tau translates one site. For a ring larger than these supports, periodization and summation give

    H_L-L e0 I = sum_j tau^j(P) >= 0.                (10)

The coboundary telescopes on the ring. There is no residual seam. Thus E_L>=L e0 for every sufficiently large L, contradicting K>0. Calling an infinite ground state a minimizer of a translation-invariant functional does not establish an identity such as (9).

[Proved] The same conclusion holds if a finite local Hamiltonian commutator is included. Let delta(Q)=i[H,Q], with the commutator evaluated as the finite sum of interaction terms meeting the support of Q. If

    h-e0 I = P + tau(B)-B + delta(Q),                (11)

then periodization gives

    H_L-L e0 I = sum_j tau^j(P) + i[H_L,Q_L],
    Q_L=sum_j tau^j(Q).

The expectation of the last commutator in a ring ground eigenvector is zero, since both occurrences of H_L can be replaced by E_L. Equation (10) follows at the level of its ground-state expectation. This argument concerns spatial coboundaries and Hamiltonian commutators with the displayed finite supports. It does not dismiss all nonlocal or state-dependent certificates.

[Proved conditional estimate] An approximate certificate makes the missing rate explicit. Suppose for every integer r there is

    h-e0 I = P_r + tau(B_r)-B_r + delta(Q_r) + R_r,
    P_r >= 0,                 ||R_r|| <= epsilon_r,  (12)

with all written local terms supported in intervals of length at most r, so periodization is valid for L>2r. The ground-state expectation gives

    E_L-L e0 >= -L epsilon_r.                        (13)

Choose r=floor((L-1)/2) along vacuum-compatible lengths. If an independent construction proves

    epsilon_r <= B M f(d)/r^2 + o(r^-2),             (14)

then (13) gives K<=4 B M f(d), or

    A <= (24 B/pi) M f(d).                          (15)

Conversely, (13) and the Casimir expansion require

    liminf r^2 epsilon_r >= K/4                     (16)

along this associated sequence of r. In particular an exact certificate or an o(r^-2) remainder is impossible for positive A. Finite-range positivity has not supplied (14). Equation (15) is a conditional reduction, not a proof of the target.

[Proved] If a proposed decomposition leaves a bounded term B_s at seam position s, translation averaging gives (1/L)sum_s B_s. In a translation eigenstate each term has the same expectation, so the averaged expectation equals the original one. The factor 1/L does not make the expectation smaller. Likewise the norm estimate gives only ||(1/L)sum_s B_s||<=max_s ||B_s||. A cancellation or a separately controlled tail is required.

## 5. A coherent reconnection inequality with no on-site symmetry

[Proved] There is a useful intermediate bound beyond bare gluing. Let psi be a normalized translation eigenvector of H_L with energy E_L. Translation and H_L commute, so such a vector exists even if the ground space is degenerate. Define the interval purities

    p_r = tr(rho_[0,r-1]^2),   p_0=1,   p_(L-r)=p_r,
    S_L = sum_(r=0)^(L-1) p_r,
    q_L = 2(b_L-e_L).                                (17)

Then

    E_(2L)-2E_L <= q_L/S_L.                          (18)

No internal symmetry, reflection symmetry, real matrix elements, or gap is used in (18). It is a variational inequality between ring ground energies. The purities describe the trial-state construction, not an additional estimator for A.

[Proved construction] Put two copies on rungs j=0,...,L-1, with layers A and B. Let H0=H_L^A+H_L^B and Omega=psi tensor psi. At a cut j define

    V_j = h_(jA,(j+1)B)+h_(jB,(j+1)A)
          -h_(jA,(j+1)A)-h_(jB,(j+1)B).             (19)

All indices are modulo L and all bond orientations are retained. The operator Hstar=H0+V_(L-1) is the ordinary 2L ring ordered as A0,...,A_(L-1),B0,...,B_(L-1).

Let S_j swap the two spins on rung j and define

    U_j = product_(a=0)^(j-1) S_a,
    phi_j = U_j Omega,                  0<=j<L.       (20)

The empty product is the identity. U_L is the full replica swap and fixes Omega. Each U_j creates a second reconnection opposite the fixed reference cut. These are finite unitaries of growing support; the final trial is a coherent sum of the resulting states.

[Proved Gram identity] The swaps commute, and U_i U_j swaps the replicas on the interval between i and j. For any bipartition with Schmidt decomposition psi=sum_a s_a |a>_R |a>_outside, contraction of the two-copy partial swap gives sum_a s_a^4=tr(rho_R^2). Hence

    G_ij = <phi_i|phi_j> = p_(|i-j|).                (21)

Translation invariance and purity of psi make G circulant, with nonnegative real entries and common row sum S_L. Consequently the vector of all ones has eigenvalue S_L. Every eigenvalue has absolute value at most the maximum absolute row sum S_L, so S_L is the largest eigenvalue. In particular the coherent sum is nonzero and has squared norm L S_L.

[Proved energy identity] Define

    K_ij = <phi_i|(Hstar-2E_L)|phi_j>.

For j>0,

    U_j H0 U_j = Hstar+V_(j-1),
    (Hstar-2E_L) phi_j = -V_(j-1) phi_j.

Conjugating V_(j-1) by U_j changes its sign. It follows that

    K_ij = <Omega|U_i U_j V_(j-1)|Omega>,             (22)

with j-1 interpreted as L-1 at j=0. For j=0, (22) follows directly from Hstar-H0=V_(L-1).

For i=j the two U operators cancel. The two crossed bonds have expectation b_L each and the two removed bonds have expectation e_L each. Thus K_jj=q_L.

For i!=j, U_i U_j swaps precisely one endpoint of the cut j-1. Swapping that endpoint exchanges the straight and crossed bonds in (19), giving

    (U_i U_j) V_(j-1) (U_i U_j) = -V_(j-1).

The product (U_i U_j)V_(j-1) is therefore anti-Hermitian. Its expectation is purely imaginary. We have proved

    Re K = q_L I.                                   (23)

K is itself Hermitian, so its imaginary part is antisymmetric. For the real coefficient vector (1,...,1), all imaginary off-diagonal terms cancel. Therefore

    <sum_j phi_j|(Hstar-2E_L)|sum_j phi_j>
       / ||sum_j phi_j||^2 = L q_L/(L S_L)=q_L/S_L.

The variational principle proves (18).

[Proved optimality in a stated class] If h and psi are real in a product basis, then K is real. Equation (23) becomes K=q_L I. When q_L>0, the minimum Rayleigh quotient over all complex coefficient vectors in the span of the phi_j is

    min_(a!=0) q_L a^dagger a/(a^dagger G a)
        = q_L/lambda_max(G) = q_L/S_L.               (24)

The same formula covers every normalized state in that span. Thus no different weighting or momentum projection in this coherent seam family can improve (18) for real XX, XXZ, Heisenberg, or SU(3) ground vectors. This optimality assertion is not made for arbitrary complex h; only (18) is symmetry-free.

[Proved elementary dimension bound] The Schmidt rank across an r-site interval is at most d^min(r,L-r). For a unit-trace positive matrix of rank R, Cauchy-Schwarz gives tr(rho^2)>=1/R. Thus

    S_L >= 1 + 2 sum_(r=1)^floor((L-1)/2) d^(-r)
                 + 1_(L even) d^(-L/2).              (25)

The difference b_L-e_L is gauge invariant. Scalar shifts add the same scalar to both terms, while a coboundary has zero expectation in both a translation eigenstate and sigma_L tensor sigma_L. In a representative with norm at most M+epsilon, b_L-e_L<=2(M+epsilon), hence q_L<=4M after taking epsilon to zero. Equations (18) and (25) improve the bare constant 4M by the denominator in (25). For d>1 that denominator tends to (d+1)/(d-1). It stays of order one. The improvement still does not bound lim L D_L. For d=1 every ring is scalar and A=0.

## 6. The entire coherent seam family fails to achieve 1/L on XX

[Proved] The failure can be shown without assuming CFT interval entropy. For even XX rings L>=4,

    q_L = 2/[L sin(pi/L)] >= 2/pi.                   (26)

The half-filled fermion correlation matrix is a rank-L/2 projector. With the parity-compatible periodic or antiperiodic grid, a geometric sum gives for nonzero displacement s

    |C_(x,x+s)| = |sin(pi s/2)|/[L |sin(pi s/L)|].    (27)

The phase of this matrix element will not be needed. For a contiguous interval of r sites, the spin reduced density matrix has the same spectrum as the Gaussian fermion reduced density matrix. If nu_a are the restricted correlation eigenvalues, diagonalizing its one-particle modes gives

    p_r = product_a [nu_a^2+(1-nu_a)^2],
    V_r = tr[C_R(1-C_R)]
        = sum_(x in R,y outside R) |C_xy|^2.          (28)

The last identity follows from C^2=C. V_r is the particle-number variance, used only as an algebraic bound on p_r.

[Proved] By complementarity it suffices to take m=r<=L/2. For each odd displacement 1<=s<=m-1, there are s pairs crossing the right boundary and s crossing the left boundary. These displacement classes are distinct since s<L/2. Keeping only those terms in the positive sum in (28), and using sin(pi s/L)<=pi s/L, gives

    V_r >= (2/pi^2) sum_(1<=s<=m-1, s odd) 1/s.

Write H_n=sum_(j=1)^n 1/j. The odd harmonic sum is H_n-(1/2)H_floor(n/2)>=H_n/2, and H_(m-1)>=log m. Hence

    V_r >= log(m)/pi^2.                              (29)

For m=1 this states only V_r>=0 and remains valid. The strict upper cutoff m-1 avoids double counting the displacement L/2 when r=L/2.

[Proved] The elementary inequality 1-x<=exp(-x), applied with x=2nu_a(1-nu_a), gives

    p_r <= exp(-2V_r) <= m^(-alpha),
    alpha=2/pi^2,              m=min(r,L-r).          (30)

Therefore

    S_L <= 1+2 sum_(m=1)^(L/2) m^(-alpha)
        <= 3 + [2/(1-alpha)](L/2)^(1-alpha).          (31)

The second inequality follows by comparison with the integral of x^(-alpha) from 1 to L/2 and then dropping a negative constant. The bound intentionally overcounts the midpoint purity once, which is harmless for an upper bound.

Combining the exact variational optimum (24), (26), and (31) proves

    optimized coherent cost
      >= (2/pi)/{3+[2/(1-alpha)](L/2)^(1-alpha)},

    L times optimized coherent cost -> infinity.     (32)

Thus no fixed constant B makes this entire trial subspace satisfy cost<=B M f(2)/L at all sufficiently large even L, for any fixed finite f(2). The actual doubling defect is O(1/L) by (5). Equation (32) disproves this variational construction as a route to the desired scaling, not the Casimir bound itself and not every adapted local-unitary construction.

[Assumed CFT calibration] The sharper interval-purity scaling is p_r asymptotic to a model-dependent constant times [(L/pi)sin(pi r/L)]^(-c_total/4), away from microscopic endpoint regions. Under the associated endpoint control, S_L grows as L^(1-c_total/4) for 0<c_total<4, as log L for c_total=4, and stays bounded for c_total>4. With q_L tending to a positive constant, (24) then gives L^(-3/4) for c_total=1 and L^(-1/2) for c_total=2. This is conditional calibration, not used to prove (32). The XX Renyi expansion and its corrections are studied in [Calabrese and Essler, arXiv:1006.3420](https://arxiv.org/abs/1006.3420). In particular this construction depends on unweighted interval entanglement, so it cannot assume a bound on the number of slow sectors.

## 7. Checks, literature, and remaining gap

[Numerical] checks.py tests all candidate finite-ring comparisons on XX, XXZ at Delta=-1/2 and +1/2, Heisenberg, and SU(3). The qubit tests compare L=4 and 8; SU(3) compares L=3 and 6. These are vacuum-compatible lengths. Every zero-cost lift and zero-error longer-ring comparison fails with the predicted positive defect. The coherent bound and its real-state optimality pass, as do |D_L|<=4M and the exact XX energy formula. The sharp finite-L Heisenberg bound without a remainder fails as already recorded in C98. The maximum eigenpair residual is 1.33e-13, and the largest Hilbert space has dimension 729.

[Numerical] purity_check.py evaluates exact XX correlation matrices through L=64, using matrices of dimension at most 32. The minimum coherent cost has L times cost equal to 1.103776 at L=4 and 1.965547 at L=64. The quantity L^(3/4) times cost approaches about 0.695 on these sizes. The assertions in (29)-(30) pass. Three small arbitrary complex bonds additionally verify (21), (23), and the symmetry-free upper bound (18). These computations establish finite-size checks, not the universal coefficient bound. The failure theorem (32) is analytic.

[Numerical] single_ring_check.py also evaluates q_L/S_L from ground vectors through qubit L=8 and SU(3) L=6, without constructing their doubled Hilbert spaces. All norm and dimension bounds pass. No approximate positive certificate was constructed, so there is no numerical test of the additional hypothesis (14).

[Assumed literature context] There are rigorous local-marginal hierarchies for ground-energy lower bounds. [Kull, Schuch, Dive, and Navascues, arXiv:2212.03014](https://arxiv.org/abs/2212.03014) give a renormalization-based hierarchy whose effectiveness depends on the coarse-graining choice. [Eisert, arXiv:2301.06142](https://arxiv.org/abs/2301.06142) discusses Anderson and semidefinite lower bounds with energy-density guarantees. The present work does not import an r^-2, crossover-free certificate from either source. The periodic joining geometry also appears in [Morin-Duchesne, Parez, and Lienardy, arXiv:2008.08952](https://arxiv.org/abs/2008.08952); a ground-state overlap by itself is not the variational energy estimate required here.

[Open] None of (i) convex comparisons, (ii) finite-range positive identities, or (iii) bare or coherent seam reconnection has supplied an explicit universal C and f(d). A viable positive-certificate proof must construct the quantitative remainder in (14), or control a genuinely nonlocal wrapping term. A viable variational proof must modify more than the position and coherent weighting of a bare reconnection. General adapted reconnection unitaries remain unexcluded and unconstructed. C98's integrated two-seam response condition is another unproved sufficient input, not an established estimate.
