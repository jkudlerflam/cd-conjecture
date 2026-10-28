# C33 finite-dimensional proofs

All statements in Sections 1 through 7 are proved. Numerical spectral and momentum results are identified separately in RESULT.md. The proof concerns the frozen legal-cell sector of the supplied 23-letter chain. It does not assert a thermodynamic CFT limit.

## 1. Parameters and the matrix being studied

[Definition] Keep all three charged/neutral hopping magnitudes, the vacuum hopping, and the junction switch equal to one. Set the loop to zero for now. Write

    E_v = -2 cos(pi/(n+2)),       Delta = 4 + E_v,
    D_x = D_M = D,              g_+ = g_- = g,
    delta = 0,                 phi_x,0 = pi/(2n),
    phi_x,+ = phi_x,- = 0,      phi_port = pi/2.

All stub and switch phases vanish. These are parameters of the original 23-letter bond, with no extra letters. The vacuum code and its endpoint visibility are inherited from the supplied construction,

    alpha = sqrt(2/(n+2)) sin(pi/(n+2)).

The one-cell excited Hamiltonian, minus E_v, has scalar diagonal 4 on the n stub vertices and 4+D on the n+1 express/marker vertices. At g=0 it is a path with 2n+1 vertices. Removing the neutral rail phase makes its scalar matrix real. The right neutral source acquires exp(-i n phi_x,0) P = -i P, where P is logical parity. All conversions remain real in this gauge.

[Proved] The excited gap is at least 2-|g|. Indeed the scalar hopping matrix has norm at most two, and the conversion matrix has norm one before multiplication by g. At a fixed token position its nonzero conversion blocks are disjoint 2 by 2 blocks. Thus, for |g|<=1, the full cell gap is at least

    min(1, 2[cos(pi/(n+2))-cos(2pi/(n+2))]).

The vacuum code is exactly 2^n dimensional. The inherited selector, U(1), parity and zero first-order assertions are unchanged at zero logical field.

[Definition and distinction] H_2(g) is the exact coefficient of lambda^2 on two cells. It is not exactly quadratic at finite g. There are two hopping matrices:

* T_1(g), its off-diagonal one-particle block, as in hermpart.py;
* T_q(g), the hopping in its Hilbert-Schmidt projection onto quadratic fermion operators, as in logical_chain.py.

Their g^2 coefficients agree. The proof below controls both at finite g. The bands mentioned below are the quadratic reference bands, not the spectrum of the full interacting H_2.

## 2. Exact scalar reduction at order g^2

[Definition] Let H_D be the real path matrix just described, p the stub junction and M the marker. For a data vertex i define

    K(t) = exp(-t H_D),
    f_s,i(t) = integral_0^t K_s,i(t-u) K_i,s(u) du,
    F_ij = integral_0^infinity f_M,i(t) f_p,j(t) dt.

The g^2 hopping is

    T^(2) = -2 alpha^4 F,                  T(g)=g^2 T^(2)+O(g^4).       (1)

[Proved] One conversion on each cell changes the neutral-pair source into a charged-pair source. The two source orderings have opposite fermion-parity signs. Their phase difference gives 2i exp(i phi_port) sin(n phi_x,0). It is -2 at the parameters above. No other two-conversion process transfers a fermion between cells. This gives (1), including the endpoint factor alpha^4. Same-cell two-conversion terms are scalar on the real tree. Thus the g^2 one-cell quadratic matrix vanishes.

An independent finite-matrix expression fixes the normalization and supplies the computation. Put

    A = H_D tensor I + I tensor H_D,
    e = |M,p>,       x = A^(-1)e,
    U_i = (|i><i| tensor I)x,
    V_j = (I tensor |j><j|)x.

Then

    F_ij = 2 U_i^T A^(-1) V_j.                                      (2)

To prove this identity, write each inverse as a Laplace integral. The two orderings of the insertions |i><i| and |j><j| divide the positive time quadrant into two regions. The associated matrix elements are transposes of one another because A is real symmetric. Their sum is the integral of the two single-cell convolutions. Equation (2) needs a (2n+1)^2 matrix, one source solve and n further right-hand sides. It never builds a logical Fock space.

[Proved] For finite g, only isolated-cell charge pairs (0,1) and (1,0) can connect |B_j> to |A_i>. For cell charge q the basis size is

    (2n+1) [binom(n,q)+binom(n,q-1)+binom(n,q+1)].

Terms with invalid binomial indices are zero. Species 0,+,- carry respectively q,q-1,q+1 data particles. This gives the polynomial-size exact finite-g algorithm in finite_g.py. In particular the largest required two-cell sector has dimension

    (2n+1)^2 (n+1) [n+1+binom(n,2)].                                 (3)

The other charge sectors appearing in a source cannot overlap the desired final source. Omitting them is exact for T_1, not a low-particle approximation to the resolvent.

## 3. The express track is a positive delay kernel

[Proved] Let C be the express block, of size m=n+1, with diagonal D+4 and off-diagonal -1. Schur complementation and the endpoint cofactor of this tridiagonal matrix give

    [(H_D+z)^(-1)]_M,i = beta(z) [(H_D+z)^(-1)]_p,i,
    beta(z) = 1/det(C+z).                                           (4)

The Laplace transform of f_s,i is the square of the scalar resolvent entry. Therefore

    f_M,i = beta(0)^2 (b * f_p,i),
    Laplace(b)(z) = product_(r=1)^m [c_r/(c_r+z)]^2,                 (5)

where c_r are the eigenvalues of C. This is an exact finite-D identity. The function b is the convolution of two normalized exponential kernels for each c_r. Consequently

    integral b = 1,       integral t b(t) dt = 2 sum_r 1/c_r
                                             <= 2m/(D+2).           (6)

No spectral-window approximation is used. Write beta=beta(0)>0 and

    F_hat = F/beta^2.

This normalization removes the express attenuation, not the poor conditioning of the data filters.

## 4. A finite-D norm bound

[Proved] Let A_s=4I-A_path,n be the isolated stub, and let f_i^s be its junction filters. Define its Gram matrix

    Gamma_ij = integral f_i^s(t) f_j^s(t) dt.

The following deliberately loose estimate is uniform in D>0:

    ||F_hat-Gamma|| <= n(n+2)/D.                                   (7)

Here and below ||.|| is the spectral norm. The proof has two parts.

First, H_D>=2I and A_s>=2I. The high block obeys C>=(D+2)I. Solving the coupled heat equations with initial state in the stub gives

    ||K_D,stub(t)-exp(-t A_s)|| <= (t/D) exp(-2t).

For completeness, the high component is bounded by

    integral_0^t exp[-(D+2)(t-u)] exp(-2u) du
      <= exp(-2t)/D.

Inserting it once in the low equation gives the stated t/D bound. Subtracting the two products in the filter convolution now gives

    |f_p,i(t)-f_i^s(t)| <= (t^2/D) exp(-2t).

Also |f_p,i(t)| and |f_i^s(t)| are at most t exp(-2t). Their L2 norms are less than 1/4, and their difference has L2 norm less than 1/(4D).

Second, ||H_D e_p||<=sqrt(18)<5 independently of D. Differentiating the convolution at the source endpoint gives

    |f_p,i'(t)| <= (1+5t) exp(-2t),       ||f_p,i'||_2 < 2.

Extend each filter by zero to negative time. Its value at zero is zero, so the extension has the same weak derivative bound. Translation in L2 gives

    ||b*f-f||_2 <= (integral t b(t) dt) ||f'||_2.

Apply these inequalities to the n-column maps from C^n to L2. Their norms are at most sqrt(n)/4, their difference at most sqrt(n)/(4D), and the derivative map at most 2sqrt(n). The delay contribution is at most nm/D. Replacing the finite-D Gram matrix by the stub Gram matrix costs at most n/(8D). Their sum is less than the right-hand side of (7).

## 5. An explicit Gram lower bound and all-n schedule

[Proved] Number data vertices i=0,...,n-1, and put m_i=2i+1. The shortest path is unique, so

    f_i^s(t) = exp(-4t) t^(m_i)/m_i! [1+r_i(t)].

For 0<=t<=1/4, walk counting gives

    0 <= r_i(t) <= 16*4^i*t^2 <= B*t^2,     B=32 n 4^n.             (8)

Indeed a walk from the endpoint to i has length i+2u. At total excess length 2r there are r+1 choices of the two walk lengths and at most 4^(i+r) walks. After dividing by the leading coefficient, the factorial ratio is at most one. Summing (r+1)(4t^2)^r for r>=1 bounds it by 16t^2. This proves (8).

Let

    Q_ij = 1/(m_i+m_j+1),
    q = det(Q)/(tr Q)^(n-1),
    det(Q) = product_(i<j) [2(j-i)]^2 / product_(i,j)(2i+2j+3).

The Cauchy determinant is positive, and lambda_min(Q)>=q. Choose the smallest integer ell>=2 satisfying

    4 n B^2 2^(-4 ell) <= q,
    t0 = 2^(-ell),
    L_n = (t0 q/32) [t0^(2n-1)/(2n-1)!]^2.                         (9)

These are explicit positive rational numbers. On t=t0 x, the normalized filter columns differ from x^(m_i) in L2[0,1] operator norm by at most sqrt(n) B t0^2<=sqrt(q)/2. Their least singular value is therefore at least sqrt(q)/2. Restoring the factors exp(-4t), dt and t0^(m_i)/m_i!, and using exp(-8t0)>=1/8, proves

    Gamma >= L_n I.                                                 (10)

Now set

    D(n) = ceil[4 n(n+2)/L_n],
    s_n = L_n/2,
    b_n = [D(n)+6]^(-(n+1)).                                       (11)

Equation (7) gives Herm(F_hat)>=3L_n I/4, hence at least s_n I. Since every eigenvalue of C is at most D+6, beta>=b_n. This proves finite-D positivity for every n. It is not an extrapolation from the integer certificates.

[Proved scale bound] One also has ||F_hat||<=n/32. Each filter has squared L2 norm at most 1/32, and convolution with b is an L2 contraction. Thus the leading Hermitian hopping margin is at least 2 alpha^4 g^2 b_n^2 s_n, while the leading hopping norm is at most alpha^4 g^2 beta^2 n/16. These bounds are very conservative.

## 6. Finite nonzero conversion and the common field

[Proved] On two excited cells write the exact inverse as (A+gC_2)^(-1), with A>=4I and ||C_2||<=2. The source norm is sqrt(3) alpha^2. Odd conversion orders vanish after source compression because each conversion changes total data number by one, while the initial and final code states have the same total charge. The geometric series gives, for |g|<=1,

    ||H_2(g)-H_2(0)-g^2 H_2^(2)||
       <= 3 alpha^4 g^4/[64(1-g^2/4)]
       <= alpha^4 g^4/16.                                         (12)

For the one-particle hopping, compression does not increase this bound. For the quadratic projection, each off-diagonal coefficient is four times a normalized trace against c_j^dagger c_i. Its absolute error is at most four times the remainder norm. Hence

    ||delta T_q|| <= n alpha^4 g^4/4,
    ||delta mu_q|| <= n alpha^4 g^4/2.                              (13)

Here mu_q sums the two one-cell matrices. These estimates are independent of the dimension 2^(2n).

An explicit sufficient choice is

    g(n) = b_n s_n/[16(n+1)],
    J_*(n) = alpha^4 g(n)^2 b_n^2 s_n,
    epsilon(n) = J_*(n)/2.                                        (14)

It satisfies g^2<=b_n^2 s_n/[16(n+1)] and g<=1. Equations (11)-(13) imply

    Herm(-T_1(g)) >= J_* I,
    Herm(-T_q(g)) >= J_* I,
    ||mu_q|| <= J_*/32.

The last constant follows from the same inequalities with room to spare. Therefore

    h(k)=mu_q+epsilon I+T_q exp(ik)+T_q^dagger exp(-ik)

is negative definite at k=0 and positive definite at k=pi. It is periodic and continuous. Each ordered eigenvalue changes sign on both semicircles. This proves crossings for all n at finite parameters, counting multiplicity. It does not by itself prove simple roots, separate centers, or momentum mismatch margins.

The common data field is inserted as a physical perturbation lambda^2 epsilon N_data. It is not part of the exactly degenerate unperturbed code Hamiltonian. Changing the excitation denominators by this field changes higher orders in lambda, not the displayed second-order coefficient.

[Proved limitation] These choices do not establish c=n. The residual H_2 is interacting, and finite-cluster estimates are not a thermodynamic integration theorem.

## 7. Integer certificates at practical D

[Proved, algorithm implemented in certify_tree.py] Take integer D and a=4. Order the path from M through the express track to the far stub endpoint. Its two-carrier inverse has diagonal entries at least eight and at most four off-diagonal neighbors, all -1. Jacobi iteration is a contraction in the infinity norm with ratio at most 1/2.

Compute N_C=det C exactly by

    d_0=1, d_1=D+4, d_(j+1)=(D+4)d_j-d_(j-1).

Let S=2^bits. Floor-rounded integer Jacobi iteration for the source N_C |M,p> gives X/S approximating x/beta from below. After

    N=bit_length(N_C)+bits+10

steps its componentwise error is at most 3/S. The truncation contribution is at most 2^(-N) N_C/4, and accumulated rounding at most 2/S. Applying the same iteration to the restricted columns of X gives Y/S approximating R V/beta. The total componentwise error, including the source error, is at most 4/S. All entries are nonnegative.

Consequently every exact normalized F entry is enclosed by

    2 sum_s X_i,s Y_i,s,j / S^2
    <= F_hat_ij
    <= 2 sum_s (X_i,s+3)(Y_i,s,j+4) / S^2.                          (15)

Symmetrize these integer numerators. Call the lower matrix N and the radius matrix E, and put e=max_i sum_j E_ij. Fraction-free Bareiss elimination checks that all leading principal minors of N are positive. Let d=det N and m be the sum of its principal minors of order n-1. Then

    lambda_min(N/S^2) >= d/(m S^2),
    ||actual Herm(F_hat)-N/S^2|| <= e/S^2.

The exact integer test d>e m proves positivity. It also gives

    (d-e m)/(m S^2) <= lambda_min(Herm(F_hat))
                              <= (n d+e m)/(m S^2).                (16)

The upper bound uses the harmonic mean of eigenvalues. No floating-point eigenvalue is used in either positivity decision. Decimal logarithms in the tables only display these rational bounds.

[Proved, computer-assisted] This test passes at D=1000 for n=2,4,8,12,16,23,32, using 256 bits through n=23 and 384 bits at n=32. At every such n the rational lower endpoint from (16) can replace s_n in (14), and beta=1/N_C can replace b_n. This gives substantially larger, still very small, explicit finite-g couplings.

## 8. The remaining center problem

[Proved] On the tree, all free and conversion matrices can be made real, while the three right-port source coefficients at our phases are purely imaginary. Thus the quadratic hopping is real. Particle-hole conjugation c_i <-> c_i^dagger, together with swapping the + and - species, makes the same-cell quadratic term antisymmetric under transposition. Hermiticity makes it purely imaginary. On the real tree it is consequently zero. The quadratic Bloch spectrum is reflection symmetric. Every occupied interval has the same center modulo pi, and interflavor Cooper transfer has zero momentum mismatch for n>=2. A common field does not remove this pair-transfer resonance.

[Proved] With a loop, the same-cell g^2 matrix is still purely imaginary by the particle-hole argument, but it need not vanish. It must be computed as an imaginary expression directly, not as a small difference of large scalar energies. For example, in the gauge with right neutral source -iP, let x=R_00|M,p>, y=R_cc|M,p>, and let U,V be their data-row and data-column restrictions, while W,Z are the restrictions of y. Then

    T = -2 alpha^4 Re[(Z^dagger R_c0 U)^T + W^dagger R_0c V],
    mu = -2i alpha^4 Im[U^dagger R_c0 U + W^dagger R_0c W
                        + V^dagger R_0c V + Z^dagger R_c0 Z].       (17)

Charged + and - propagators are identical here. These formulas follow by listing the four matching one-conversion source pairs for T and pairing creation and annihilation processes for mu. They are implemented in fast_loop.py and agree with the independently enumerated source blocks.

[Open] Positivity is an open condition, so a sufficiently small loop preserves the crossings. It does not follow that every difference of distinct pair sums of Fermi points becomes nonzero. To prove this for all n, one must exclude an identically vanishing resonance function of the loop parameters and bound the smallest nonzero mismatch. Neither (4), Gram positivity, nor nonzero individual center shifts proves that assertion. Numerical results and the precisely retained limitations are in RESULT.md.

## 9. A finite-loop, finite-g certificate at n=23

[Proved, computer-assisted] The center obstruction in Section 8 is removed at n=23 without enlarging the alphabet. Use

    n=23, D_x=D_M=1000, Delta=4-2 cos(pi/25),
    t=t_x=t_switch=1, delta=0,
    phi_x,0=pi/46, phi_x,+ = phi_x,- = 0,
    phi_port=pi/2,
    phi_d,0=pi, phi_d,+ = phi_d,- = pi/2.

Let d_0=1, d_1=1004, d_(j+1)=1004 d_j-d_(j-1), and beta=1/d_24. Set

    t_d = beta 2^(-136),
    g_+=g_-=2^(-365),
    E_unit = 2 alpha^4 beta^2 g^2,
    epsilon = E_unit 2^(-214).                                    (18)

The physical common field is lambda^2 epsilon. Energies and velocities quoted for the second-order logical Hamiltonian acquire a factor lambda^2. The certificate is independent of that common positive factor.

After removing the express phases, all three species now have the same scalar path-plus-loop matrix. The loop entry from the far stub to M is -i t_d. The three source phases are still -iP, i, -i. Thus (17) reduces to

    T/g^2 = -4 alpha^4 Re(U^dagger R V),
    mu/g^2 = -4i alpha^4 Im(U^dagger R U + V^dagger R V),            (19)

at leading conversion order, with one common complex two-carrier inverse R. This uses a phase choice within the existing weak-link terms, not an extra local state or a new interaction range.

### Entry enclosures

[Proved] full_loop_integer.py evaluates (19) by fixed-point integer arithmetic at 512 bits. Normalize source vectors by beta, as in Section 7. The complex Jacobi contraction in the infinity norm is at most

    rho <= 1/2 + t_d/4.

Each real and imaginary component is rounded by integer division. There are two divisions per step, one for the rational loop coefficient and one for the positive diagonal. The complex rounding error per step is less than 2/S, where S=2^512. The accumulated error is less than 5/S. With N=bit_length(d_24)+512+20 steps, the tail is less than 1/S. The restricted-source error contributes less than another 2/S in the next inverse. The code uses the conservative common bound 10/S for both source and resolved restricted vectors. The assertion d_24 2^136>100N controls the small change from the tree contraction ratio.

If x,y denote scaled complex approximations with componentwise errors at most 10/S, then the real or imaginary part of their product has error at most

    [10(|x_R|+|x_I|+|y_R|+|y_I|)+100]/S^2.

Here the x and y on the right are integer numerators. Summing this bound encloses every entry of T and mu. Antisymmetrizing the onsite numerator preserves the exact particle-hole identity and averages its two valid enclosures. The resulting rational matrices and entry radii are saved in full_loop_n23.npz. No finite-loop Taylor truncation is made.

### Exact conditioning of the narrow channels

[Proved] Direct floating-point diagonalization is inadequate at these scales. Fraction-free symmetric elimination of the integer tree Gram enclosure gives an integer matrix R_0 with

    R_0 N R_0^T = diag(d_i),       d_i>0.

The identity is checked in integer arithmetic. Choose b_i=floor(log2(d_i)/2) and define

    W_0 = 2^256 diag(2^(-b_i)) R_0.

Its transformed tree metric is diagonal with entries in [1,4). All congruences of the full-loop matrices are performed on integer numerators before any conversion to floating point. Their exact entry-error row bound in this metric is less than 1.697e-50. A rational Gershgorin test proves a positive Hermitian hopping margin greater than 0.5028 in this metric.

### Root existence, simplicity and pair-sum margins

[Proved, computer-assisted] Floating-point generalized eigenvalues only propose 46 disjoint brackets. Each endpoint is then replaced by a dyadic rational t=tan(k/2). At such a point,

    cos(k)=(1-t^2)/(1+t^2),       sin(k)=2t/(1+t^2),

so the Bloch matrix and its propagated entry enclosure are rational. Its real and imaginary parts are enclosed by fixed-point matrices at 80 bits. Numerical eigenvectors propose a further 46-bit dyadic complex congruence. Integer diagonal-dominance checks prove that this congruence is invertible and prove the inertia of the enclosed exact Bloch matrix. Numerical eigenvalue signs are never used for acceptance.

For each of the 46 brackets the endpoint inertias differ by one. The counts run from 0 to 23 along the negative-momentum roots and back to 0 along the positive-momentum roots. The smallest remaining endpoint congruence margin is greater than 1.9305e-8. The determinant of an n by n nearest-neighbor Bloch symbol, multiplied by z^n, is a polynomial in z=exp(ik) of degree at most 2n. Thus 46 distinct bracketed roots exhaust the degree bound. Each root is simple, and there are no additional roots.

The quartic test includes all unordered pairs with replacement, 1081 pairs at n=23. At a rational bracket midpoint, exp(ik) is a rational point on the unit circle. Form every pair product, sort their phases by exact rational orientation tests, and test every neighboring pair on the circle. The coordinate difference gives a rational lower bound on chord distance. The root uncertainty is subtracted using

    |2 atan(t)-2 atan(t_mid)| <= 2|t-t_mid|.

All remaining bounds are strictly positive. Since angular distance is at least chord distance, the exact rational tests prove

    delta_2 > 0.0140193,
    delta_4 > 3.1558e-7.                                          (20)

The decimal inequalities are rounded outward from saved rational lower bounds. Pair transfer between two bands is among these tests. Their sea centers are therefore distinct modulo pi. The exact brackets and rational lower bounds are in root_certificate_exact_n23.json.

### Extension from g^2 to the exact finite-g quadratic projection

[Proved] On the looped graph, the free two-carrier inverse has A>=3I. For g<=1/2, the same even-conversion expansion used in (12) gives

    ||remainder|| <= alpha^4 g^4/4.

The quadratic Bloch error is at most 4n alpha^4 g^4. In E_unit units, after congruence by W_0, this is at most

    2n g^2 beta^(-2) tr(W_0 W_0^dagger).                            (21)

The trace is computed rationally. The code verifies, with integers, that g=2^(-365) makes (21) at most 2^(-33). Every root-endpoint inertia certificate remains strict under this additional error, and the positive Hermitian hopping certificate remains strict as well. The degree argument still applies to the exact quadratic projection at finite g. Therefore (20), all 46 simple crossings, and distinct band centers hold for that finite-g matrix, not just for its leading coefficient. The one-particle hopping is also positive because its remainder bound is smaller.

[Proved limitation] The nonquadratic part of H_2 is not discarded by this statement. Its norm is controlled by the same conversion expansion, but an interacting lattice-to-CFT conclusion also requires the infrared and length-sector arguments. This run proves the finite band conditions at n=23 and all-n tree positivity. It does not prove distinct centers for every n.
