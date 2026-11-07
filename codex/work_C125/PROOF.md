# C125 proof and audit

2026-11-06. Lattice spacing and hbar are one. Natural logarithms are used. A tag at the start of a paragraph applies to its displayed calculations. The results concern the input stated in C125. They do not settle the proposed bound for actual fixed-d nearest-neighbor chains.

## 1. What positivity gives

[Proved] Let T be positive and let Z(t)=Tr T^t be finite on an interval containing t and 2t. Restrict to the support of T and write epsilon=-log T. Scalar shifts of epsilon cancel from R. Then

    rho_t = T^t/Z(t),
    Tr rho_t^2 = Z(2t)/Z(t)^2,
    R(t)=2 log Z(t)-log Z(2t)=S_2(rho_t).

Where the derivatives exist,

    (log Z)'(t)=-<epsilon>_t,
    (log Z)''(t)=Var_t(epsilon)>=0,
    R'(t)=2[<epsilon>_(2t)-<epsilon>_t]
         =-2 integral_t^(2t) Var_u(epsilon) du <=0.                 (1)

Approximation by finite spectral projections gives the monotonicity without assuming differentiability at endpoints. Also R>=0. Positivity of matrix entries is not enough for this argument.

[Proved] Suppose a positive transfer representation, or positive finite-Trotter approximations, reproduces the same physical ring family on even lengths. Taking the Trotter limit in (1), followed by beta tending to infinity at each fixed ring length, gives

    D(L)=E_(2L)-2E_L>=0,
    D(L+2)<=D(L).                                                  (2)

Degeneracies contribute only O(1) to R and disappear after division by beta. This establishes the corrected items 1 and 2 on the allowed sublattice. It does not establish that L D(L) decreases. Multiplication by L is the entire missing issue.

[Proved] A positive transfer operator supplies a real-length partition function only in its domain of Schatten summability. It need not have finite Tr T^t at every t>0. The exact XX example below is finite for t>1 and diverges for t<=1. Moreover, the existence of the physical beta limit at integer lengths alone does not prove its existence at every interpolated real length. If that limit exists, E(t)=-lim_beta beta^(-1) log Z(t) is concave in t, since log Z is convex. Concavity of E(t) is weaker than convexity of -t E(t) in t^2.

[Proved] There is a useful finite-Trotter formulation that avoids claiming convergence of transfer operators acting on changing Hilbert spaces. Prove positivity and (1) at each Trotter number, then pass to the physical traces at allowed integer lengths. These limits are ordinary finite-dimensional Trotter limits in the original spin Hilbert space. The statements in (2) need no infinite-column reconstruction.

[Proved spectral reconstruction for even L>=4] For fixed beta, order the eigenvalues lambda_(m,j) of T_m decreasingly, padding with zeros. Since Tr T_m^2=Z_m(2,beta) converges, it is bounded by some C_beta. Thus lambda_(m,j)<=sqrt(C_beta/j). A diagonal subsequence converges at every j. For each real t>2 the summable envelope (C_beta/j)^(t/2) permits passage through the sum. The limiting diagonal positive operator T_beta therefore has Tr T_beta^L=Z(L,beta) for every even L>=4, and finite traces for every real t>2. This constructs a common positive spectrum rather than assuming one. At L=2, spectral mass could escape to zero in this argument; no equality at that endpoint is claimed from compactness alone. It is true explicitly for XX below. Reconstructing local insertions and continuous time rotations requires further structure beyond this eigenvalue compactness argument.

## 2. Reflection positivity and the models in the question

### 2.1 An explicit two-site transfer construction

[Proved] For an even ring, let

    Z_m(L,beta)=Tr[(exp(-beta H_even/m) exp(-beta H_odd/m))^m].

A column is a periodic history s=(s_0,...,s_(2m-1)), with s_(2m)=s_0. Put g=exp(-beta h/m), and define matrices on columns by

    K_e(s,t)=product_(j even) g_(s_(j+1),t_(j+1);s_j,t_j),
    K_o(s,t)=product_(j odd)  g_(s_(j+1),t_(j+1);s_j,t_j).

Contracting the checkerboard tensor network, with no approximation beyond the displayed Trotter product, gives

    Z_m(L,beta)=Tr[(K_e K_o)^(L/2)].                               (3)

Suppose the two-site gate has a positive realignment in a compatible staggered basis,

    g=sum_alpha w_alpha A_alpha tensor conjugate(A_alpha),
    w_alpha>=0.

Realignment means that row and column indices are (a,c) and (b,d) for g_(a,b;c,d). Its matrix is sum_alpha w_alpha |A_alpha><A_alpha|. Both K_e and K_o are tensor products of this matrix, with permutations of the column coordinates, so both are positive. Therefore

    U_m=K_e^(1/2) K_o K_e^(1/2)>=0,
    Tr U_m^n=Tr(K_e K_o)^n,
    T_m=U_m^(1/2),
    Z_m(L,beta)=Tr T_m^L                                         (4)

for even L. No assertion that the usual one-site QTM is itself positive is needed. Squaring an arbitrary non-Hermitian QTM would not prove (4).

[Proved] If -h=sum_alpha J_alpha A_alpha tensor conjugate(A_alpha), J_alpha>=0, up to a scalar, then the gate condition holds for every positive Trotter step. Expand exp(-t h) in powers. Each word is B tensor conjugate(B), its coefficient is nonnegative, and the convergent sum retains positive realignment. This is the direct local reflection-positivity certificate used here.

[Assumed published context] The necessary and sufficient across-plane coupling-matrix criterion is developed by Jaffe and Janssens, [Characterization of Reflection Positivity: Majoranas and Spins](https://arxiv.org/abs/1506.04197), and, for general matrix algebras, [Reflection Positive Doubles](https://dspace.library.uu.nl/bitstream/handle/1874/350519/doubles.pdf?isAllowed=y&sequence=1), Theorem 7.9. The gate calculation above supplies the sufficiency needed here. Necessity can also be seen by differentiating the reflection quadratic form at beta=0 on a traceless test operator. A negative eigenvalue of the across-plane coupling form gives a negative first derivative.

### 2.2 XX, XXZ and Heisenberg

[Proved] Use h=Sx Sx+Sy Sy+Delta Sz Sz, without a magnetic field. Rotate alternate sites through pi about y. On an even ring this is a consistent on-site unitary transformation and gives

    h'=-Sx Sx+Sy Sy-Delta Sz Sz,
    -h'=Sx tensor conjugate(Sx)+Sy tensor conjugate(Sy)
         +Delta Sz tensor conjugate(Sz).                         (5)

For Delta>=0, (5) proves (4). This includes XX, the antiferromagnetic spin-1/2 Heisenberg chain, both the critical interval 0<=Delta<=1 and the Ising-like region Delta>1. The transfer statement does not require a critical infrared. The same argument covers the antiferromagnetic bilinear Heisenberg interaction at every spin. Changing the sign of the XY exchange is an alternating rotation about z on even rings.

[Proved] For spin 1/2, no on-site antiunitary reflection gives this across-bond positivity for Delta<0 with nonzero XY exchange. On the three-dimensional space of traceless Hermitian operators, an antiunitary acts by an orthogonal transformation O of determinant -1. Reflection invariance and positivity require -diag(1,1,Delta) O to be symmetric positive semidefinite. If Delta is nonzero this matrix is invertible; its determinant is Delta, up to a positive normalization. A positive matrix has positive determinant. Thus Delta<0 is excluded. At Delta=0 the explicit certificate (5) applies.

[Numerical] The actual two-column matrices in (3), not just a local sign test, have negative and complex eigenvalues for the tested negative anisotropies. For beta=2, m=4, the unnormalized minimum real parts are approximately -0.06703 at Delta=-0.9 and -0.03413 at Delta=-0.5. These finite-Trotter findings do not by themselves constitute an all-beta theorem about every possible representation of exact ring traces. They do refute positivity of those Trotter transfer matrices. There is also an exact-ring numerical obstruction, without the two-site ring convention: at Delta=-0.9 and beta=21.928721377, R_4=1.073477335 and R_6=1.137473641. This violates (1), so no common positive T can reproduce those exact even-ring traces at that beta. No such R witness was found for Delta=-0.5 in the tested sample; failure of its local RP criterion remains the analytic conclusion.

[Assumed published numerical evidence] Complex QTM eigenvalues and incommensurate finite-temperature correlations in the attractive XXZ regime were studied by Fabricius, Klumper and McCoy, [Temperature dependent spatial oscillations in the correlations of the XXZ spin chain](https://arxiv.org/abs/cond-mat/9812012). Their result is additional evidence against extending the positive spatial propagation construction to that regime, not a substitute for the local coupling test.

[Proved] Odd rings are excluded from the staggered construction. This is a physical restriction, not an innocuous interpolation convention. On XX, E_3=-1/2 whereas E_6=-2, so D(3)=-1. A positive T reproducing all integer spin-ring traces would have R(3)>=0 and could not produce this limit. Similar odd-ring failures occur for the Heisenberg chain. Vacuum compatibility and even bipartiteness must be stated separately for other models.

### 2.3 SU(N) swap and symmetric-representation interactions

[Proved] The SU(2) swap is P=2 S.S+I/2 and is covered on even rings. The antiferromagnetic SU(N) fundamental swap for N>2 does not have the same bond-centered reflection positivity. Write

    P=I/N+2 sum_a t^a tensor t^a,
    tr(t^a t^b)=delta_ab/2.

Positivity of the across-plane form would force the antiunitary reflection to act as minus the identity on all traceless generators. This is impossible for N>2. Choose a diagonal traceless Hermitian X with tr X^3 nonzero. Antiunitary conjugation preserves the real number tr X^3, while X -> -X changes its sign. Alternating fundamental and conjugate representations is a different Hamiltonian and does not repair the same-representation swap chain.

[Proved] The same obstruction applies to the integrable interaction on Sym^k(C^N), N>2. Here is the extra step needed beyond naming a complex representation. The Hermitian operator space decomposes without multiplicity into SU(N) irreducibles labelled l=0,...,k; l=1 is the adjoint. In a Hilbert-Schmidt orthonormal Hermitian basis let C be the real symmetric coefficient matrix of h, with the scalar component removed. It is scalar on each of those irreducibles. Its coefficient on the adjoint is strictly positive. Indeed let

    X=sum_a t_R^a tensor t_R^a.

Both X and h are scalar on the two-site irreducibles [2k-j,j], j=0,...,k. Their eigenvalues are strictly ordered in the same direction: h_j=H_(k-j), while X_j is one half of the total quadratic Casimir minus the two fixed one-site Casimirs. The latter decreases strictly as j increases. Since tr X=0,

    tr(h X)=sum_j dim_j (h_j-mean(h)) X_j >0.

The strict inequality follows by writing the covariance as a sum of products (h_i-h_j)(X_i-X_j) with positive dimension weights. Thus C is positive on the adjoint subspace.

[Proved] Let O denote the orthogonal action of the proposed on-site antiunitary on Hermitian operators. Positivity requires -C O to be symmetric positive semidefinite, with the harmless transpose convention absorbed in O. The polar decomposition of the symmetric matrix -C forces O=-sign(C) on its nonzero spectral subspace. One can also square the positive factor to obtain C^2 and then solve on its support. On every positive-eigenvalue subspace of C this forces O=-I. In particular Theta(t_R^a)=-t_R^a. A diagonal Cartan generator on Sym^k(C^N) has nonzero cubic trace for a suitable choice when N>2. One can take the fundamental Cartan eigenvalues (N-1,-1,...,-1); the weights on Sym^k are N n_1-k with multiplicity binomial(k-n_1+N-2,N-2), and their cubic sum divided by dim Sym^k is k(N-1)(N-2)(N+k)(N+2k)/[(N+1)(N+2)], which is positive for N>2 and k>=1. Antiunitary conjugation again cannot reverse them all. This excludes the on-site DLS-type reflection for these same-representation chains, including every Sym^2 family in the supplied data with N>2.

[Proved] For SU(2), Sym^2 is the spin-1 Takhtajan-Babujian bond. In the normalization of the question, with X=S.S,

    h=(X-X^2)/4+3I/2,
    h|_(J=0,1,2)=(0,1,3/2).

Take Hilbert-Schmidt orthonormal Cartesian spin and quadrupole operators A_(1,m), A_(2,m). The realignment of -h after the alternating spin time-reversal rotation has eigenvalues

    -7/2 on the scalar operator,
     3/4 on the three spin operators,
     1/4 on the five quadrupoles.

This follows by inserting the displayed quadratic polynomial and resolving the scalar, antisymmetric and traceless symmetric parts of S_a S_b. The scalar can be shifted freely. All remaining coefficients are positive, so (4) applies on every even ring.

[Proved general-k extension] The same positivity holds for every SU(2) Sym^k integrable bond h=H_J on total-spin J=0,...,k. Here is a finite recoupling proof, including the sign of every coefficient. Put d=k+1, w_J=2J+1, and use Hilbert-Schmidt normalized rank-l tensor operators. Their realignment coefficients for -h are

    lambda_l=-sum_(J=0)^k w_J r_l(J) H_J,
    r_l(J)=(1/d) sum_(r=0)^min(l,J)
       [(-l)_r(l+1)_r(-J)_r(J+1)_r]
       /[(1)_r(-k)_r(k+2)_r r!].                                 (5a)

Here (a)_r is the rising factorial. This is the finite SU(2) recoupling coefficient, equivalently (-1)^(k+J+l) times the six-j symbol {k/2,k/2,J;k/2,k/2,l}. The displayed terminating polynomial fixes its normalization and signs without a six-j convention ambiguity.

[Proved] Define a_J=(J+1)[d^2-(J+1)^2]/2, with a_(-1)=a_k=0, and

    (D f)_J=[a_J(f_J-f_(J+1))+a_(J-1)(f_J-f_(J-1))]/w_J.

Substitution of the terminating sum in (5a) gives D r_l=l(l+1) r_l. This can be checked coefficient by coefficient in the polynomial variable J(J+1); the ratio of successive coefficients is

    [(-l+r)(l+1+r)] / [(r+1)^2(-k+r)(k+2+r)],

with the J-dependent factors retained as (-J)_r(J+1)_r. D is self-adjoint with weights w_J by cancellation of neighboring edge terms. Also r_l(0)=1/d. For l>=1, summing D r_l with weights shows sum_J w_J r_l(J)=0. Since H_(J+1)-H_J=1/(J+1), direct subtraction gives

    (D H)_J=1/2-(d^2/2) delta_(J,0).

Therefore

    l(l+1) sum_J w_J r_l(J) H_J
      =sum_J w_J r_l(J)(D H)_J=-d/2,
    lambda_l=d/[2l(l+1)]>0, l=1,...,k.                           (5b)

For l=0, r_0=1/d and interchanging the finite harmonic sums gives lambda_0=-d H_k+k/2. A scalar bond shift changes only lambda_0, so all the required non-scalar coefficients are positive. This proves the reflection-positive transfer construction on even rings for all SU(2) symmetric integrable representations. spin_certificates.py checks the recoupling spectrum directly for k=1,...,8; the proof does not extrapolate those samples.

[Numerical] SU(3) has complex eigenvalues already at beta=2, m=2 in (3), with maximum imaginary part about 5.5214 before normalization. Its exact ring energies also refute an ordinary positive T reproducing all even ring traces: E_4=-2 and the computed E_8<-4 give D(4)<0. The numerical margin is large and can equally be used as a variational upper bound from the computed normalized eigenvector. This does not exclude a different representation restricted only to vacuum-compatible lengths 3n or 6n. No positive transfer construction on those restricted subsequences is proved here. Likewise, absence of the DLS certificate for N>2 symmetric chains is not a theorem excluding every conceivable auxiliary spectral representation on a selected subsequence.

## 3. Exact crossed spectrum of XX

[Proved] This provides both a positive infinite-column spectrum and a test of the proposed Matsubara assumptions. Keep h=Sx Sx+Sy Sy, so M=1/2 and the fermion dispersion is cos k. Define

    p_c(beta)=(1/(2pi)) integral_0^(2pi) log[2 cosh(beta cos k/2)] dk,
    p_s(beta)=(1/(2pi)) integral_0^(2pi) log|2 sinh(beta cos k/2)| dk,
    delta_beta=p_c-p_s>0,
    kappa_n^o=asinh((2n+1)pi/beta), n>=0,
    kappa_n^e=asinh(2n pi/beta), n>=1.

For every even L>=2 the physical spin partition function is

    Z_L=exp(L p_c) { (1/2)[product_(n>=0)(1+exp(-L kappa_n^o))^4
                          +product_(n>=0)(1-exp(-L kappa_n^o))^4]
              +2 exp(-L delta_beta)
                    product_(n>=1)(1+exp(-L kappa_n^e))^4 }.       (6)

[Proved derivation] Start from the two Jordan-Wigner parity sectors as in C118. On either periodic or antiperiodic momentum grid,

    P_s=product_j [1+exp(-beta cos k_j)],
    Q_s=product_j [1-exp(-beta cos k_j)],
    Z_L=(P_0+P_(1/2)+|Q_nonzero|)/2.

Exactly one Q grid vanishes for even L, and the other has the enhancing sign after the physical parity projection. Use

    cosh z=product_(n>=0)[1+4z^2/(pi^2(2n+1)^2)],
    sinh z=z product_(n>=1)[1+z^2/(pi^2 n^2)].

For nu>0 set kappa=asinh nu. Factoring cos k-i nu in z=exp(i k) and taking the product over an L-point grid gives

    product_j(nu^2+cos^2 k_j)
      =(exp(kappa)/2)^(2L)
        [1-(-1)^(L/2) exp(2pi i s) exp(-L kappa)]^4.

The two s choices interchange the plus and minus signs. The constant factors, integrated over k, yield exp(L p_c) or exp(L p_s). On the nonzero sinh grid, product_j |cos k_j|=2^(2-L), supplying the factor four in |Q_nonzero|. Combining the parity sectors gives (6). Convergence follows from exp(-L kappa_n)=O(n^(-L)) for L>1; the original L=2 product is also convergent.

[Proved] Formula (6) is a positive spectral expansion. The first sector consists of four fermionic levels of each energy kappa_n^o, restricted to an even total number of excitations. The second consists of two identical copies of a fermionic Fock spectrum with four levels of each energy kappa_n^e and a vacuum offset delta_beta. The common vacuum energy is -p_c. Thus all multiplicities are nonnegative integers. It supplies a positive T_beta independent of L and reproduces the even spin-ring traces. Its real-length trace is finite for L>1. At L<=1, already the two-excitation part of the first sector diverges.

[Proved] The momentum labels can be assigned by the Fourier modes on the time circle. At each positive frequency nu the four levels have temporal momenta +nu,+nu,-nu,-nu. In the first sector nu=(2n+1)pi/beta; the even occupation constraint makes the total momentum an integer multiple of 2pi/beta. In the second sector nu=2n pi/beta and the two vacua have momentum zero. This is the Fourier diagonalization underlying the same determinant factors, with two Fermi points and two signs of frequency. In particular the spectrum contains states

    p_m=2pi m/beta,
    epsilon_m-epsilon_0
      =asinh(pi/beta)+asinh((2m-1)pi/beta), m>=1.                  (7)

They use one level at +pi/beta and one at +(2m-1)pi/beta. These are not invisible extra sectors: the local density Sz has a nonzero two-fermion form factor in this sector, as also follows from the bubble calculation below.

[Proved] At fixed beta, (7) grows logarithmically with m. Therefore for any finite V,

    epsilon-epsilon_0 >= |p|/V                                  (8)

is false for the full lattice transfer spectrum, even for locally accessed states in XX. More strongly, for any real a!=0, one sign of p_m gives

    exp[-L(epsilon_m-epsilon_0)-a p_m] -> infinity.

Consequently Tr[T_beta^L exp(-a P_tilde)] diverges for every nonzero a. A unitary family exp(i theta P_tilde) commuting with T does not imply a finite positive real-tilted trace. It supplies a real-angle trace, with modulus bounded by the untilted trace, but not an analytic continuation to a neighborhood of a=0. Items (ii) and (iii) in the proposed lattice input thus require correction, not just an improved constant.

[Proved independent local-correlator check] In the infinite XX chain the one-fermion Matsubara resolvent is

    G_x(i nu)=integral_(0)^(2pi) (dk/2pi) exp(i k x)/(i nu-cos k).

Its pole inside the unit circle is z=-i sign(nu) exp(-asinh|nu|), so

    |G_x(i nu)|=exp(-|x| asinh|nu|)/sqrt(1+nu^2).

Wick contraction expresses the Sz correlation at omega_m=2pi m/beta as a sum over fermionic nu=(2j+1)pi/beta of products G_x(i nu)G_(-x)(i nu-i omega_m), with an overall sign immaterial here. The least decay exponent is attained at nu=pi/beta and nu=omega_m-pi/beta, and is (7). Both contributions have the same sign on the even spatial subsequence and a nonzero coefficient. All other exponents are larger, except this endpoint degeneracy. Thus the rate in (7) is detected by an original local spin operator. This also verifies its momentum assignment without relying on a choice of labels for the full QTM spectrum.

[Numerical] The ordinary fermion trace agrees with full spin ED through L=8 to 2.28e-13. The crossed product (6) agrees with it through (L,beta)=(1024,32768); its largest absolute log-trace discrepancy is 5.59e-9 on a log trace of order 10^7. Doubling the product cutoff from 20000 to 40000 modes gives no resolved change in the sampled cases. The limiting gap beta delta_beta tends to pi/2, with value 1.570796328 at beta=32768. These are numerical evaluations of the proved product, not a finite-Trotter extrapolation.

## 4. What Lieb-Robinson actually proves at Matsubara frequency

[Proved, using the C109 bound] Let A and B have disjoint single-site supports at separation x, and consider an infinite-volume beta-KMS state. For positive bosonic Matsubara frequency omega=2pi k/beta, the Lehmann representation, or its KMS spectral version, gives

    G_AB(i omega)=integral_0^beta d tau exp(i omega tau)<A(tau)B>
                =i integral_0^infinity dt exp(-omega t)<[A(t),B]>.
                                                                    (9)

To verify the sign, a matrix element with delta=E_n-E_m contributes (w_m-w_n)/(delta+i omega) on the first side. The commutator has coefficient (w_n-w_m)exp(i delta t), whose Laplace integral multiplied by i is the same expression. At nonzero omega the disconnected part vanishes.

[Proved] C109 supplies, in the optimized gauge and its epsilon limit,

    ||[A(t),B]|| <= C min{1, exp(-mu x+4M exp(mu)t)},
    C=2||A||||B||, mu>0.                                        (10)

Set V=4eM and first use mu=1. Split (9) at t=x/V. For 0<omega<V, direct integration gives

    |G_AB(i omega)|
      <= C[(exp(-omega x/V)-exp(-x))/(V-omega)
                         +exp(-omega x/V)/omega]
      <= C V/[omega(V-omega)] exp(-omega x/V).                   (11)

For omega=V the bound is C(x+1)exp(-x)/V. For omega>V it is bounded by

    C[(omega-V)^(-1)+omega^(-1)] exp(-x).                         (12)

There is no all-frequency linear rate in these expressions. The early-time tail of the commutator, absent in a relativistic theory with exact microcausality, is what prevents it.

[Proved] Optimizing mu gives a convenient asymptotic rate

    alpha(omega)=omega/(4eM),              0<omega<=4eM,
                 log(omega/(4M)),         omega>=4eM.             (13)

More explicitly, at fixed mu the exponential rate is min(mu,mu omega/(4M exp(mu))). If omega<=4eM it is maximized at mu=1; otherwise the maximum is at mu=log(omega/(4M)). At the equality point a factor linear in x occurs. One may either retain that factor or replace alpha by any slightly smaller exponent. Equations (11)-(13) give constants and the precise frequency range. They are uniform in beta apart from the allowed discrete frequencies.

[Proved conditional spectral implication] If a spatially reflection-positive representation expresses the reflected A,A Matsubara correlator as a positive spectral measure of transfer gaps, (11) forbids support below omega/V for 0<omega<V. Otherwise a positive mass below that threshold would contradict the exponential upper bound at large even x. The same statement with alpha(omega) applies at high frequency. For a general A,B correlator cancellations prevent this inference; for sectors not generated by these insertions there is no inference at all. Thus a low-frequency locally detected spectrum condition is available, but the global condition required for the real-tilted trace is false.

[Numerical] In XX, beta=1 and k=100 give omega=628.3185, the exact local rate in (7) is about 8.99, while omega/(4eM) is about 115.57. The nonlinear bound (13) is compatible with this rate. Frequencies beta=1,10,100 and k=1,2,10,100,1000 are checked in checks.json.

## 5. A positive, finite-dimensional spectral countermodel even granting (i)-(iii)

### 5.1 Definition

[Proved construction] Fix the labels d=2 and M>0, choose v=1 after setting M=1 if desired, and require v<=V=4eM. Let n be any positive integer. Use n copies of the compact self-dual free-boson CFT, each with c=1 and physical velocity v. A basis of its crossed circle Hilbert space has

    Delta=sum_(s=1)^n [(m_s^2+w_s^2)/2+N_s^++N_s^-],
    j=sum_(s=1)^n [m_s w_s+N_s^+-N_s^-],
    m_s,w_s integers, N_s^+,N_s^- nonnegative oscillator levels.

Each oscillator level has the usual integer partition multiplicity. The vacuum is unique and |j|<=Delta term by term, because |m w|<=(m^2+w^2)/2.

[Proved construction] Fix b=M/2. For each beta>0 retain the largest complete set of lowest Delta shells whose total dimension D_beta obeys

    D_beta<=floor(exp(b beta)).                                  (14)

There is always at least the vacuum. Ties are retained as complete shells. This defines a finite-dimensional Hilbert space depending on beta, with an ordinary counting trace, not fractional spectral multiplicities. On it set

    epsilon=epsilon_0(beta)+(2pi/(beta v))Delta,
    p=(2pi/beta)j,
    T_beta=exp(-epsilon),
    mu_beta=sum_retained_states delta_(epsilon,p).                (15)

For definiteness choose the scalar vacuum energy

    epsilon_0(beta)=
      -[pi n/(6 beta v)] beta^2/[beta^2+pi n/(3vM)]
      -(log 2) exp(-beta M).                                    (16)

Only its large-beta limit is relevant for the CFT; its choice never affects R. Equation (16) is asymptotic to the CFT vacuum energy -pi n/(6 beta v). The smooth adjustment at small beta is included to satisfy elementary original-channel norm bounds as well.

### 5.2 Verification of every stated input

[Proved] T_beta is positive, and all its powers have finite traces. P_tilde is diagonal, commutes with T_beta, and exp(i beta P_tilde)=I. W(L,a) is finite for every real (L,a), and its logarithm is jointly convex: its Hessian is the covariance matrix of (epsilon,p) in the normalized finite spectral trace. Also

    epsilon-epsilon_0=(2pi/(beta v))Delta
                    >=|p|/v>=|p|/V.                             (17)

Thus the model satisfies a stronger version of (iii), in every sector, and a stronger version of the claimed tilt analyticity.

[Proved] For an ordinary D-dimensional density matrix, its squared trace is at least 1/D by Cauchy-Schwarz. Therefore

    0<=R_L(beta)<=log D_beta<=b beta                              (18)

for every L>0. In particular it obeys the requested trivial bound at any L_0>=1,

    R_(L_0)<=2L_0 log d+4 beta L_0 M.

It even satisfies the stronger inherited reconnection ceiling |R_L|<=4 beta M for every L. Its R decreases with L by (1). None of these estimates depends on n.

[Proved optional stronger check] For L>=1, (16) gives

    0<=-epsilon_0(beta)<=M beta/2+log 2,
    log Z_L<=-L epsilon_0+log D_beta
             <=L log 2+M beta L.

Moreover -epsilon_0>=log 2 exp(-beta M)>=log 2-beta M, so

    log Z_L>=L log 2-M beta L.

At beta tending to zero, D_beta=1 and Z_L tends to 2^L. Thus the elementary dimension and bond-norm trace bounds do not exclude the example either. It is not the partition function of a single finite-d nearest-neighbor Hamiltonian for all beta. In particular its shell cutoff can jump with beta. Compatibility with one microscopic Hamiltonian is additional information absent from (i)-(iii), not a property claimed for this model.

[Proved infrared CFT limit] Fix n and an aspect ratio beta/L>0, and let L,beta tend to infinity. The retained conformal-level cutoff tends to infinity, while

    L[epsilon_0(beta)+pi n/(6 beta v)] ->0.

The trace is bounded by the full compact-CFT excitation trace at that aspect ratio. Dominated convergence of the positive convergent series therefore gives the full modular-invariant n-boson CFT torus partition function, not just its vacuum Casimir term. The same argument works for real tilts in the interior |a|<L/v at fixed scaled aspect ratio. The infrared theory has A=n v and a unique vacuum.

### 5.3 Its real-length ground energy and the failed convexity

[Proved] The compact-boson counting asymptotic is

    log N(Delta)=2pi sqrt(n Delta/3)+o(sqrt(Delta)).                (19)

For completeness, two chiral oscillator partition products have

    log sum_Delta mult(Delta) exp(-s Delta)
       =n pi^2/(3s)+O(log(1/s)), s down to zero.

This follows by integrating -2n log(1-exp(-s k)); the zero-mode Gaussian sum changes only the logarithmic remainder. Minimizing s Delta+n pi^2/(3s) gives the exponent in (19). A matching lower bound is obtained from the oscillator partition asymptotic, or by a saddle-point coefficient estimate of the same positive product; adding the Gaussian zero modes does not change its leading square-root exponent. Thus (19) is the usual explicit free-boson counting asymptotic, not an assumed density of hidden states.

[Proved] Put

    B=pi n v/6,
    e=(epsilon-epsilon_0)/beta.

Equation (19) becomes log N(beta e)/beta ->2 sqrt(B e), where N on the left counts crossed excitation energy. The dimension cutoff (14) therefore restricts e to

    0<=e<=e_max=b^2/(4B).

The last included complete shell has the same limiting endpoint, since adjacent conformal levels are separated by O(beta^(-2)) in e and (19) has a continuous leading exponent. The Laplace principle, which here follows directly by upper and lower Riemann-sum bounds on positive level sums, gives

    g(L)=lim_(beta->infinity) beta^(-1) log Z_L
        =sup_(0<=e<=e_max) [2 sqrt(B e)-L e]
        = b-L b^2/(4B),     L<=ell=2B/b,
          B/L,              L>=ell.                            (20)

The scalar (16) contributes zero to this limit. Hence E(L)=-g(L), e_0=0, and

    Theta(L^2)=L g(L)
             =bL-b^2 L^2/(4B), L<=ell,
              B,               L>=ell.                         (21)

For 0<L<ell,

    d^2 Theta/d(L^2)^2=-b/(4L^3)<0.                              (22)

Theta is increasing up to ell and constant afterward. Its first derivative with respect to L^2 joins continuously to zero. It is concave, not convex, before the crossover. E itself is concave, as positivity requires: it is linear below ell, equals -B/L above ell, and has matching first derivatives.

[Proved] Its doubling defect is

    D(L)=2g(L)-g(2L).

For L<=ell/2, D(L)=b. For L>=ell,

    D(L)=3B/(2L)=pi n v/(4L),
    lim_(beta->infinity) L R_L(beta)/beta=pi n v/4.                (23)

Thus at fixed d=2 and M, the family has unbounded L R_L/beta as n increases. Each member has a finite number n of CFT sectors; no infinite-central-charge infrared is being used. Choosing n arbitrarily large also places the interval of failed Theta convexity beyond any proposed L_0<=g(d). Positivity, the entire real-tilt trace, the stronger spectrum condition (17), the trivial bounds, and the full infrared CFT form all hold. They imply neither proposed conclusion.

[Numerical] checks.py independently constructs compact-boson level multiplicities from two integer zero modes and two oscillator partition products, then takes n=1,2,4,8 copies and the finite-rank truncation. It checks R<=log D, positivity, and the approach toward the thermodynamic mechanism. The finite enumeration uses M=1 and b=1/2, as in (14). Momentum-resolved checks in audit_checks.py also test the covariance Hessian, the full cone inequality and the sign of R'. The analytic family (20)-(23) is sampled up to n=1024. No enormous Hilbert space was allocated.

## 6. Why the null-edge argument does not transfer

[Proved] For a sector of physical speed v_s, crossed conformal excitations satisfy epsilon-epsilon_0=|p|/v_s on chiral states. Its unrestricted chiral trace has its null edge at |a|=L/v_s. A common spectrum ceiling with V>=v_s guarantees at most |a|<L/V. If v_s<V, that region stops strictly before its null edge. Letting a approach L/V therefore cannot implement that sector's null limit. No inference should replace all v_s by V in a boost identity that has not been proved.

[Proved] There is a second, independent obstruction. Even if all speeds equal v and W is entire, its joint convexity only constrains derivatives of its own tilted free energy. It does not identify them with derivatives of the untilted pressure at a different radius. The finite-rank model has all speeds equal and W entire, yet (22) holds. It therefore isolates the missing rotation relation from the issue of multiple velocities.

[Assumed sufficient structural addition] An exact rotation or boost Ward identity determining the tilted thermodynamic free energy from its untilted slice, with convergence up to the appropriate sector null limits, would restore the continuum argument supplied in the question. It is not a lattice consequence of the ingredients audited here. With multiple velocities one needs the corresponding sectorwise identities and access to their separate edges, or a direct aggregate substitute. Merely adding more analyticity or the full spectral cone does not suffice, by the countermodel.

[Proved sufficient ground-state addition] A weaker scalar addition than a full Ward identity would be

    L D(L)<=L_0 D(L_0) for all allowed L>=L_0.                    (24)

Then the trivial bound at L_0 and beta tending to infinity yield D(L_0)<=4L_0 M, and the infrared Casimir limit gives

    A <= (4/pi) L_0 D(L_0) <= (16/pi) L_0^2 M.                   (25)

Using the sharper inherited gluing bound D(L_0)<=4M improves the last expression to (16/pi)L_0 M. Equation (24) is an additional scale comparison, not a result of ordinary Renyi monotonicity. It rules out the constant-D interval in (20). There is no uniquely defined mathematically weakest extra axiom; (24) states a sufficient missing scalar comparison without imposing an artificial bound on fast crossover lengths.

[Proved sufficient Theta addition] If E(L) has the infrared expansion with bulk e_0 and Theta(q)=-L E(L) is convex for q>=L_0^2, then

    C(q)=Theta(q)+e_0 q -> pi A/6.

A convex function with a finite limit at infinity is nonincreasing and at least its limiting value. Consequently

    A<=c_A(L_0),
    c_A(L)=(6/pi)C(L^2),
    c_A(L_0)<=(12/pi)L_0^2 M.                                   (26)

The last bound follows by choosing a centered optimized bond, so |e_0|<=M and |E(L_0)|<=L_0 M. The infinitesimal condition is

    Theta''(q)=-(L^2 E''(L)+L E'(L)-E(L))/(4L^3)>=0.              (27)

Positivity gives only E''<=0, not (27). Equations (21)-(22) separate them explicitly.

[Proved scope correction] The condition in (26) controls ground energies and the beta>>L limit of the ratio. It is not by itself an all-temperature ratio theorem. Indeed

    R_L=beta D(L)+2 log Tr exp[-beta(H_L-E_L)]
                       -log Tr exp[-beta(H_(2L)-E_(2L))].         (28)

A bound on the thermal remainder, or a finite-temperature scale comparison, is still needed. Also the C116 estimate R_L<=K+k M f(d) beta/L becomes L R_L/beta<=K L/beta+k M f(d), not K+k M f(d) at arbitrary aspect ratio. These two formulations must not be identified silently. The countermodel refutes either uniform conclusion because (23) already violates it in the low-temperature limit.

## 7. Checks on physical finite-size energies

### 7.1 XX

[Proved] On even spin rings E_L=-csc(pi/L), e_0=-1/pi, A=1. For L>1 the convergent expansion of x/sin x has positive coefficients in even powers of x. Thus

    Theta(q)=q/pi+pi/6+7pi^3/(360q)+31pi^5/(15120q^2)+...,
    c_A(L)=1+7pi^2/(60L^2)+31pi^4/(2520L^4)+... .

Each inverse power of q is convex and decreasing. This proves Theta convexity and c_A monotonicity for this real interpolation on L>1. It is also the interpolation from (6): its beta-limit gives g(L)=L/pi+(2/pi) integral_0^infinity cosh(u) log(1+exp(-L u)) du=csc(pi/L), where expansion of the logarithm and summation of its alternating partial fractions proves the last equality. It is stronger information specific to the free dispersion, not a consequence of (i)-(iii). The coefficient 31pi^4/(2520 L^4) follows from multiplying the csc expansion by 6/pi.

[Numerical] The complete spin trace is tested on even rings through paired length 16384, with beta/L from 0.01 to 128. The inherited C118 bound R<=2 log 2+(5pi/3) beta/L passes every sample. The large-beta QTM spectrum is checked via (6), not by claiming that m=4 approximates a low-temperature continuum Trotter limit.

### 7.2 A genuine lattice failure of Theta convexity

[Numerical] For XXZ Delta=-0.9, the Bethe bulk energy and ED ground energies give

    L       4             6             8             10            12
    c_A     0.190911175   0.214541263   0.229533760   0.238419817   0.243662250.

The Theta secant slopes in q=L^2 decrease. Thus both proposed signs fail on actual spin-1/2 rings. This model lies outside the positive-transfer class established in Section 2.

[Proved small-anisotropy mechanism] The failure can be pushed beyond any fixed length cutoff by approaching Delta=-1 from above. After an alternating pi rotation about z, write Delta=-1+u, u>0:

    H=-sum_i S_i.S_(i+1)+u sum_i Sz_i Sz_(i+1).

At u=0 the ferromagnetic ground multiplet has total spin L/2. In its zero-magnetization state, permutation symmetry and (sum Sz)^2=0 give

    <Sz_i Sz_j>=-1/[4(L-1)], i!=j.

The projected perturbation is minimized at magnetization zero for even L. Finite-dimensional degenerate perturbation theory at fixed L therefore gives

    E_L=-L/4-u L/[4(L-1)]+O_L(u^2),
    Theta(L^2)=L^2/4+u L^2/[4(L-1)]+O_L(u^2).                   (29)

The smooth function f(q)=q/(sqrt(q)-1) obeys

    f''(q)=(3-L)/[4L(L-1)^3]<0 for L>3.

Take any fixed even triple of lengths above a proposed cutoff. The leading strict discrete concavity from (29) dominates the O_L(u^2) remainders for sufficiently small u. Hence no dimension-only onset of Theta convexity exists for the unrestricted XXZ family. This statement concerns the ground-energy inequality, not a fixed-d counterexample to the Casimir bound itself.

### 7.3 Heisenberg, other XXZ, SU(3), and Sym^2 Bethe data

[Numerical] The tested Delta=-0.5,0,0.5,1,2 rings have nonnegative discrete Theta curvature on lengths 4,6,8,10,12. c_A decreases in the critical cases among those samples. SU(3) passes both signs on the vacuum-compatible triple 3,6,9. These are finite samples. In particular their success does not prove a positive transfer operator for SU(3).

[Proved exact point added to the supplied table] The supplied Sym^2 JSON starts at L=2N, not L=N. For L=N, the SU(N) singlet constructed by taking two epsilon tensors and symmetrizing their two indices at each site has every pair in the [2,2] channel. Its bond energy H_0 is zero. Every bond in the stated normalization is positive semidefinite, so E_N=0 is exact. This supplies the missing first point; it is not a Bethe data entry read from the file.

[Numerical] With that point added, all requested N=2,3,4,5,7,8 families have c_A nonincreasing and Theta discretely convex on the supplied lengths through 120 (or the largest available multiple of N). The nonuniform-coordinate check is

    [Theta_(i+1)-Theta_i]/[L_(i+1)^2-L_i^2]
      >=[Theta_i-Theta_(i-1)]/[L_i^2-L_(i-1)^2].                  (30)

For N=8, c_A(8)/M=7.125295634 and A/M=6.597344573. Using only the JSON's first point, L=16, would instead give 6.768553904. The smallest positive difference of consecutive slopes in the N=8 table is 5.07e-8. Input energies are floating-point Bethe solutions; this is a numerical sign test, not an interval-certified proof of their exact sign.

## 8. What was and was not done

[Proved] The abstract sufficiency question has a negative answer even with stronger hypotheses than requested, by the finite-dimensional construction (14)-(23). The proposed global crossed spectrum condition and real tilt fail already for XX, by (6)-(13). Positive-transfer Renyi monotonicity remains valid in the reflection-positive class and does imply the unscaled monotonicity (2).

[Numerical] Small full-spectrum ED, small finite-Trotter QTM spectra, exact large-ring free-fermion traces, and the specified Bethe-data sign tests were run locally. See CHECKS.md and checks.json for the scope. No large interacting computation was run locally. run_typhon.slurm is prepared but was not submitted, since the sandbox cannot access Typhon.

[Assumed or unresolved] The supplied continuum rotation argument is accepted as conditional input, not independently refereed here. The spectral compactness construction gives Schatten summability for real length greater than two; an operator reconstruction with local insertions and continuous time rotations in full generality is not supplied. Nor is a positive auxiliary representation on every selected SU(N) vacuum subsequence. These gaps are not used in the countermodel or the XX refutation. No actual fixed-d Hamiltonian violating A<=C M f(d) has been constructed. No independent referee or sub-agent review was used.
