# C37 proofs and certificate specification

[Proved, scope] The result below concerns the frozen 23-letter candidate, for every integer 1 <= n <= 23. The finite-field calculations and the inherited C33 positivity certificates prove statements about specified finite matrices over the rationals. They do not prove a microscopic interacting CFT theorem. The local elimination and bosonization inputs remain assumed where stated.

## 1. Parameters and normalization

[Proved construction] Keep the original alphabet and every original letter rule. Set

    t = t_x = t_sw = 1,
    D_x = D_M = 1000,
    Delta_n = 4 - 2 cos(pi/(n+2)),   delta = 0,
    g_+ = g_- = g,
    phi^x_0 = pi/(2n),   phi^x_+ = phi^x_- = 0,
    phi^d_0 = pi,        phi^d_+ = phi^d_- = 0,
    phi_port = pi/2,
    t_d = 10^(-2000).

All stub and switch phases are zero. The loop phase pi is a permitted coefficient of the existing bond. No letter or rule has been added. The marker selector and illegal-sector penalty are separate inherited requirements.

[Proved convention] Write

    alpha_n = sqrt(2/(n+2)) sin(pi/(n+2)),
    a4_n = alpha_n^4,
    e = 10^(-1000),
    W_n = lambda^2 g^2 a4_n.

The physical uniform data field is epsilon_physical = W_n e. The normalized leading quadratic symbol is

    h(z) = T z + M + e I + T^dagger / z,       z = exp(i k).       (1)

T and M in this note have neither g^2 nor a4_n included. The leading physical symbol is W_n h. The supplied structured.py includes a4_n, whereas modular.py omits it. This normalization matters in reading the C33 bounds.

[Proved] The vacuum energy without the marker chemical potential is -2 cos(pi/(n+2)). Subtracting it makes the excited stub diagonal 4 and the express/marker diagonal 1004. The vacuum visibility contributes a common a4_n. Thus all entries of T and M are in Q(i), even though the individual express phase and the vacuum visibility need not lie in Q(i).

[Proved] Gauge the neutral carrier at data site j by the Jordan-Wigner string Z_0 ... Z_(j-1). At the marker this is the full data parity. The departure-Z factors on both rails then disappear. Gauge away the express phases on the tree. The total phase is n phi^x_0 = pi/2, so the neutral marker source has phase -i. The loop phase in this gauge is phi^d_0 - n phi^x_0 = pi/2. Neutral loop hopping is therefore -i t_d in one orientation and +i t_d in the other. Charged loop hopping is real. The data conversions remain local Jordan-Wigner creation or annihilation operators. These are the matrices used by modular.py.

[Proved] The folded graph is a path with a closing edge of weight t_d <= 1. Every row has total hopping magnitude at most 2. Each shifted one-carrier Hamiltonian A_s is therefore at least 2 I, independently of n. The two-carrier Hamiltonian A_st = A_s tensor I + I tensor A_t is at least 4 I. Conversion has one-cell norm at most one, since its nonzero blocks are disjoint two-state blocks. The two-cell conversion C has norm at most 2.

[Proved] In the target cell the excited-token gap is at least 2-|g|. The full legal-cell gap is bounded below by

    min(2-|g|, 2[cos(pi/(n+2))-cos(2pi/(n+2))]).             (2)

In particular Delta_n >= 2 > |g| for the schedule below. This is a target-cell bound, not a new unrestricted-ring selector theorem.

## 2. Leading quadratic coefficient without a logical Fock-space calculation

[Proved] Let Z be the pair of endpoint source maps, with a4_n removed. Its norm is at most 3. Let R_0 be the inverse of the two-carrier Hamiltonian at g=0. Expanding the pair resolvent gives the g^2 coefficient

    K = - Z^dagger R_0 C R_0 C R_0 Z.                       (3)

The sign is negative. The g=0 coefficient is scalar on the logical register, and odd conversion counts vanish. A conversion changes neutral versus charged carrier status, whereas both source and sink have an even number of charged carriers.

[Proved] At order g^2 there are only two possibilities. One conversion in each cell gives an intercell fermion bilinear, with the marker parity supplying its Jordan-Wigner grading. Two conversions in one cell give a scalar plus a same-cell bilinear. The loop has scalar holonomy after the gauge above, so it does not add logical operators. Consequently (3) is exactly quadratic plus a scalar, for arbitrary n. This assertion concerns the Taylor coefficient, not the complete finite-g bond.

[Proved] It is therefore enough to evaluate (3) between the vacuum and the 2n one-particle logical states. The single-carrier dimension is s=2n+1 and every inverse block has dimension s^2. There are only the species-pair blocks 00, 0c, c0 and cc, with c denoting either charged species. The c0 and 0c blocks are related by swapping the cells.

[Proved] The implementation first computes R_0 Z for each source species. Applying one conversion restricts a carrier coordinate to data position j. The resulting vectors are grouped by their output species pair. A block of their bilinear contractions is

    - V_minus^T A_st^(-1) V_plus.                           (4)

The minus subscript means substitution i -> -i. It is ordinary complex conjugation over Q(i), not a conjugation operation inside the residue field. The logical bit dependence enters only through conversion signs, occupation admissibility and the marker parity. These signs combine the blocks (4) into the vacuum and one-particle matrix elements. Subtracting the vacuum diagonal from the one-particle diagonal gives the quadratic coefficients. The two onsite blocks add to M.

[Proved, implementation size] modular.py performs these operations in a finite field, caching the conversion sources. It never builds a 2^n logical space. At n=23 its largest inverse is 2209 by 2209 and the final determinant uses a 23 by 23 matrix.

[Numerical, independent literal check] literal_coefficient.py instead differentiates the supplied literal cell by applying (3) directly in its complete excited-pair space at n=1,2,3. The pair dimensions are 324,3600,28224. It agrees with the structured hopping to relative errors 5.8e-15,1.9e-13,1.4e-12. Onsite errors divided by the hopping scale are 6.5e-14,9.8e-12,2.3e-10. The non-Gaussian residual of this coefficient is at most 1.2e-9 of the hopping at n=3. These are floating-point consistency checks of the derivation, not its proof.

[Numerical limitation] An earlier finite-g Richardson subtraction lost accuracy in the onsite coefficient, with errors up to 4.5e-5 of the hopping at n=3 and no clean decrease on reducing g. literal_results.json retains that failed precision check. The direct coefficient calculation above replaces it as the numerical validation. Neither subtraction is used by the exact certificate.

## 3. Crossing, inherited from C33

[Proved, inherited C33 certificate method] At t_d=0, M=0 and T is real. C33's normalized product-path inverse and integer enclosures certify positivity of S=-(T+T^T)/2. This work does not rederive its positivity theorem. Existing certificates at n=2,4,8,12,16,23 were reused. The same C33 certifier was evaluated only at the missing sizes 1,3,5,6,7,9,10,11,13,14,15,17,18,19,20,21,22. Every size passed at D=1000.

[Proved by integer comparisons] In C33's saved notation, let L be its lower integer matrix, E its entry-radius matrix, d_C=det C_express, b its fixed-point precision, d_L=det L, and m_L the sum of principal codimension-one minors. Put r_E=max_i sum_j E_ij. Its certified bound in the present normalization is

    lambda_min(S) >= 2(d_L-r_E m_L)/(m_L 2^(2b) d_C^2).      (5)

verify_bounds.py checks directly, with integers, that the right side is greater than 10^(-250), for every n=1,...,23. No comparison of floating-point eigenvalues or displayed logarithms is used for this decision. crossing_bounds.json records the exact comparisons and their source certificates.

[Proved] The loop changes the pair Hamiltonian by norm at most 2 t_d. Since both inverses have norm at most 1/4,

    ||R_loop-R_tree|| <= t_d/8.

Telescoping the three resolvents in (3) then bounds the change of the complete g^2 coefficient by

    3 * 9 * 2^2 * (1/4)^2 * t_d/8 < t_d.                  (6)

Extracting the quadratic symbol increases this by at most 32n on 1/2 <= |z| <= 2. One conservative derivation bounds each coefficient of a traceless fermion bilinear by four times the input operator norm and uses row sums for the two onsite and two hopping blocks.

[Proved] For r=10^(-250), e+32n t_d < r/2 throughout the certified range. At k=0 every eigenvalue of (1) is negative, at most -3r/2. At k=pi every eigenvalue is positive, at least 3r/2. Each ordered band has a zero on each open half-circle. This gives at least 2n zeros, counted with multiplicity. The determinant polynomial below has degree 2n, so there are exactly 2n once its simplicity is established.

## 4. The pair-product certificate

[Proved] Define

    f(z) = det[T z^2 + (M+e I)z + T^dagger].                (7)

Its coefficients lie in Q(i). Let z_1,...,z_m be its roots, m=2n, and let

    N = m(m+1)/2 = n(2n+1),
    Q_pair(w) = product_(1 <= i <= j <= m) (w-z_i z_j).      (8)

Repeated indices are included. The polynomial is symmetric in the roots, so its coefficients are in Q(i). If its discriminant is nonzero, all unordered pair products are distinct. This also forces the original roots to be distinct and nonzero when f has nonzero constant and leading coefficients.

[Proved] On the unit circle, z_i z_j = exp(i(p_i+p_j)). Thus a repeated pair product is exactly a quartic resonance modulo 2pi. Distinct unordered pairs are compared; the identical-pair density/current processes are not excluded. Taking pairs (i,i) and (i,j) also excludes coincident Fermi points. Taking the two endpoints of two seas excludes coincident centers modulo pi.

[Proved] Exact total U(1) forbids number-changing terms. Therefore the pair-product condition covers every potentially relevant number-conserving quartic momentum, including pair transfer, exchange and umklapp. Including repeated indices is conservative, since some corresponding local fermion representatives require derivatives. Pairing and all-plus four-field momenta need not be tested in this U(1)-preserving candidate.

[Proved, finite-field test] Use the prime p=65521 and the map

    Z[i] -> F_p,       i -> 24297,
    24297^2 = -1 mod p.

Trial division verifies primality in audit.py. The rational inputs 10^(-2000) and 10^(-1000) are evaluated by modular inversion, not by underflowed floating-point numbers. All required carrier inverses exist modulo p, as shown by the nonzero pivots and the exact residual tests. The leading and constant coefficients of f are nonzero modulo p.

[Proved, polynomial algorithm] Evaluate the n by n determinant (7) at 2n+1 distinct field points and interpolate. To form (8), let s_k=sum_i z_i^k, computed from the monic f by Newton identities. The power sums of the unordered pair products are

    t_k = (s_k^2+s_(2k))/2.                                 (9)

Indeed s_k^2 counts every unequal unordered pair twice and every equal pair once. Adding s_(2k) gives two copies of every pair. Newton identities applied to t_k give Q_pair. Here N<=1081<p, so every required division by k and by 2 is valid. Compute gcd(Q_pair,Q_pair') over F_p.

[Proved by exact computation, n=1,...,23] In every case this gcd is 1. The complete determinant and pair-product polynomial, their gcd, field parameters and three extra determinant evaluations are saved in pairs_n*.json. Therefore the characteristic-zero discriminant is nonzero. An identically zero discriminant could not have a nonzero reduction at a place where all its denominators and the leading coefficient are invertible.

[Proved] Combining this test with Section 3, all 2n roots lie on the unit circle and are simple. Every ordered band has exactly two Fermi points. Every nontrivial pair sum differs modulo 2pi, including the pair sums that determine the sea centers.

[Proved, arithmetic implementation] Sparse elimination, interpolation, Newton identities and polynomial gcd use integer residues. The only floating-point operations used for finite-field arithmetic are NumPy dot products of nonnegative integers. Before each such product, the code checks

    inner_dimension * (p-1)^2 < 2^53.

Every product and every possible partial sum is then an exactly representable integer in binary64. Reducing the result modulo p and converting to int64 is exact. Sparse ordering uses SciPy only to choose a permutation; no result of its floating-point factorization enters a finite-field solve. Every solve used in the certificate has an exact modular residual check.

[Proved, independent finite audits] audit.py compares sparse solves against a separate dense modular inverse at product dimensions 9,25,49,81,121, compares the pair-product polynomial against direct multiplication for up to 46 prescribed roots and 1081 pairs, and compares the cached coefficient extraction with a direct uncached contraction at n=1,2,3. All checks agree exactly. These checks supplement the algebraic derivations and the stored certificates.

## 5. Explicit momentum and velocity margins

[Proved] The modular certificate establishes nonvanishing without estimating a tiny root difference. The following height bound converts it into a stated positive margin. It is deliberately conservative.

[Proved] Set q=10^2000 and d=(2n+1)^2. Each integer Gaussian matrix q A_st has norm at most q*2014. Its determinant is a positive integer, bounded by (q*2014)^d. There are only three distinct determinants, for 00,0c,cc. Their product, denoted d_*, clears the denominators of every inverse block. The three resolvents in (3) have common denominator d_*^3. A common denominator for T,M and e is therefore

    Q_0 = 10^1000 d_*^3,
    Q_0 <= Qbar_n := 10^1000 (2014*10^2000)^(9(2n+1)^2).    (10)

The exact logical extraction adds and subtracts entries and introduces no new denominator.

[Proved] A bound for the norm of h(z), throughout 1/2 <= |z| <= 2, is B_n=32n. For example ||K||<=9*2^2/4^3=9/16 in (3); subtracting the vacuum scalar and adding the two onsite blocks gives entry bounds far smaller than B_n after row sums. Also e<1. Define

    H_n = (B_n Qbar_n)^n,
    N_n = n(2n+1),
    delta_n = (2 H_n)^(-N_n(N_n-1)).                        (11)

[Proved] The integer Gaussian polynomial

    F(z)=det[Q_0 T z^2+Q_0(M+eI)z+Q_0 T^dagger]

has nonzero leading coefficient a_m=det(Q_0 T), with 1<=|a_m|<=H_n. Each w_i=a_m z_i is an algebraic integer: substituting it in F gives the monic equation with coefficients a_(m-1), a_(m-2)a_m,...,a_0 a_m^(m-1). The monic polynomial with roots w_i w_j has Gaussian-integer coefficients, since these coefficients are both algebraic integers and elements of Q(i). Its nonzero discriminant therefore has modulus at least one.

[Proved] Scaling every pair root by a_m^2 gives

    |Disc(Q_pair)| >= |a_m|^(-2N_n(N_n-1)).                 (12)

There are N_n(N_n-1)/2 pair-root differences in the square root of the discriminant, and each has modulus at most 2. If d_pair is the smallest such difference, (12) implies

    d_pair >= H_n^(-N_n(N_n-1))
              * 2^(-(N_n(N_n-1)/2-1)) >= delta_n.          (13)

[Proved] Circular angular distance is at least chordal distance. Consequently the leading quadratic symbol obeys

    delta_2 >= delta_n,
    delta_4 >= delta_n,
    min_(a != b) dist(theta_a-theta_b, pi Z) >= delta_n/2.   (14)

For n=1 the center comparison is empty; the quadratic and conservative quartic tests still apply. The result includes umklapp modulo 2pi, rather than only mismatches as real numbers.

[Proved] At a root, differentiating det h(e^(ik)) gives modulus

    |a_m| Q_0^(-n) product_(j != i) |z_i-z_j|
       >= Qbar_n^(-n) delta_n^(2n-1).

The product of the other n-1 eigenvalues of h has modulus at most B_n^(n-1). The zero eigenvalue is simple. Its speed therefore satisfies

    v_min >= v_n := Qbar_n^(-n) delta_n^(2n-1)/B_n^(n-1).  (15)

This is a normalized speed; multiply by W_n for the physical cell-coordinate speed.

## 6. Explicit robustness and finite conversion

[Proved] Put

    r_n = delta_n/64,
    tau_n = Qbar_n^(-n) r_n^(2n)
            /[2^(n+4) n (B_n+1)^(n-1)].                   (16)

Suppose a perturbation E(z) is analytic on 1/2<|z|<2, obeys E(1/conj(z))=E(z)^dagger, and has norm at most tau_n on that annulus. This is a sufficient, not necessary, error criterion.

[Proved] On the circle |z-z_i|=r_n, every other root has distance at least delta_n-r_n>=delta_n/2. Since |z|<=2, the factorization of det h gives

    |det h(z)| >= 2^(-n) Qbar_n^(-n)
                  r_n (delta_n/2)^(2n-1).                 (17)

For matrices of norm at most B_n and errors at most one, multilinearity of the determinant gives

    |det(h+E)-det h| <= n ||E|| (B_n+1)^(n-1).              (18)

Equation (16) is smaller than needed in (17)-(18). Rouche's theorem gives one simple zero in each disk for the full homotopy h+tE, 0<=t<=1. On the rest of the unit circle, every distance to a root is at least r_n; the weaker bound Qbar_n^(-n) r_n^(2n) and (18) exclude additional zeros there.

[Proved] Hermitian Laurent symmetry pairs a zero z with 1/conj(z). A simple unit-circle zero cannot leave the circle under a continuous perturbation preserving this symmetry without colliding with its partner. The disks each contain one simple zero throughout the homotopy, so no such collision occurs. Each Fermi point moves by angular distance at most 2r_n. A quartic mismatch moves by at most 8r_n=delta_n/8. In particular all perturbed quadratic and quartic margins remain greater than delta_n/2, and the center separation remains greater than delta_n/4 modulo pi.

[Proved] For |g|<=1/2, the part of the complete second-order bond beyond g^2 has norm, with a4_n removed, at most

    (9/4) (g/2)^4/[1-(g/2)^2] < g^4.                     (19)

This is the norm-convergent resolvent series, with only even powers retained. Its normalized quadratic symbol after dividing by g^2 differs from (1) by at most 32n g^2 on the annulus. The non-Gaussian part of the two-cell bond has normalized norm at most

    C_NG(n) g^2,      C_NG(n)=2+16n^2.                     (20)

One conservative bound sums the at most (2n)^2 bilinear coefficients, each bounded by four times the remainder norm, and includes the scalar subtraction.

[Proved explicit schedule] Choose

    g_n^2 = tau_n/(256n).                                  (21)

Then g_n<1/2 and the quadratic error is at most tau_n/8. The exact finite-g, second-order quadratic projection therefore retains all margins above. Any smaller positive g also works.

[Proved normalized interaction bound] Combining (15),(16),(20),(21), and C_NG(n)<=18n^2, gives

    epsilon_NG/(v_n delta_n)
       <= 18/[256 * 2^(n+4) * 64^(2n)] < 10^(-6).           (22)

Also epsilon_NG*log(1/delta_n)/v_n < 10^(-6), since delta log(1/delta)<=1 for 0<delta<1. These are bounds on the normalized lattice coefficients. They are not a derivation of continuum matching constants or an interacting root-dressing theorem. A further reduction of g accommodates any specified finite matching constants without changing the leading roots.

## 7. Higher orders in the port coupling and the infrared scope

[Proved, finite clusters] At fixed n a finite cluster has an isolated unperturbed code band and analytic perturbation theory at sufficiently small lambda. Odd powers vanish on open cell chains. The data field epsilon_physical=W_n e is order lambda^2. Its effect inside virtual denominators, as well as the other next corrections, begins beyond the retained order and is included in the order-lambda^4 remainder. A finite-cluster constant can be bounded from a contour resolvent; it is independent of lambda but not uniform in cluster size by this argument.

[Assumed, existing local-elimination input] To discuss the infinite chain, assume the massive-sector elimination supplies a normalized quadratic error with an analytic annulus norm bounded by

    C_n lambda^2/(a4_n g_n^2),                              (23)

with finite C_n independent of ring length. An exponentially weighted local bound that controls the annulus is sufficient. The same assumption must bound the retained non-Gaussian coefficients and the interaction-induced self-energy, with the necessary matching constants. This is a quantitative version of the existing leakage/local-elimination input, not a theorem proved by the finite matrices here.

[Proved conditional schedule] Given such C_n>=1, take

    0 < lambda^2 <= a4_n g_n^2 tau_n/(4 C_n),                (24)

and decrease it if necessary for the inherited selector/exclusion hierarchy. Equations (21),(24) put the normalized quadratic corrections below tau_n. Thus the explicit momentum margins survive those corrections. The tiny physical field still removes the bare particle-hole umklapp, since the permitted root shifts were budgeted against the complete delta_n, not against the crossing amplitude alone.

[Assumed/heuristic, infrared conclusion] With (23), controlled dressing, positive Gaussian stiffness and the standard bosonization input of ir_control and C20, the weak non-Gaussian couplings yield the n-mode Luttinger description. This remains a physics-standard conditional conclusion. The present work does not establish the local-elimination hypothesis, a thermodynamic central charge, or a mobile length-fluid gluing rule.

## 8. Band scales and limits of the result

[Proved, inherited tree identity] Let beta_n=1/det C_express, where C_express has n+1 vertices, diagonal 1004 and off-diagonal -1. Then

    1006^(-(n+1)) <= beta_n <= 1002^(-(n+1)).                (25)

C33's transfer identity expresses the marker filter as beta_n^2 times a convolution of the junction filter with a unit-mass positive kernel. Since the shifted carrier matrix is at least 2, each conversion filter has modulus at most t exp(-2t), whose squared L2 norm is 1/32. Cauchy-Schwarz and the convolution contraction then give

    ||T_tree|| <= n beta_n^2/16.

Together with (6), a safe loop bound is

    J_n := W_n ||T||
      <= W_n [n beta_n^2/16 + 32n*10^(-2000)].              (26)

[Proved] The crossing excursion at k=0,pi is at least (3/2) W_n*10^(-250), and the leading minimum speed is at least W_n v_n. The microscopic vacuum gap is of order n^(-2), while a4_n is of order n^(-6). The express attenuation in the leading hopping is beta_n^2, of order D^(-2n-2). The conversion and port schedules shrink the physical scale much further. The bounds are existence bounds, not predictions of attainable simulation scales.

[Numerical evaluation of proved bounds] scales.py evaluates only logarithms of (10)-(26). At n=23, log10(delta_23) is approximately -1.06948e15 and log10(g_23) approximately -2.45981e16. The exact definitions (11),(16),(21), not these rounded displays, specify the margin and coupling. The tree hopping upper bound before multiplying by lambda^2 g^2 is about 10^(-149.73). C33's much sharper actual scale estimate is compatible with this upper bound.

[Proved limitation] This establishes the finite-dimensional conclusion for every n<=23. It is not an all-n theorem. The positivity requirement needs an increasing D beyond a fixed certified range. The modular pair-product noncollision also needs an all-n nonidentity argument or more finite certificates.

[Proved extension mechanism] At any fixed n with the present modular certificate, replacing D by 1000+65521 j, and replacing either decimal exponent E,L by E+65520 k,L+65520 l, preserves the residue-field calculation. The inherited positive Gram limit allows sufficiently large D, then arbitrarily small field and loop preserve crossing. Thus a successful noncollision certificate is compatible with the scale hierarchy, rather than being confined to a real parameter point outside the crossing region. This argument does not assert that the same modular certificate succeeds for untested n.

[Proved scope of C20] There is no length density in this frozen candidate. If a separately established mobile implementation has the same certified payload symbol, its mixed low-dimension resonances require the additional finite C20 density test. The C32 core has different matrices, so these certificates cannot be assigned to it without another calculation.
