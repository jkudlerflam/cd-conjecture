# C42. All-n nonresonance of the two-bus payload

[proved, scope] For each positive integer n there is a choice of finite parameters in the 46-letter initialized core of C39 with exactly 2n simple Fermi points in the quadratic projection of its exact second-order bond. Every two different unordered pairs of endpoints have different sums modulo 2pi, including pairs with repeated endpoints. In particular the sea centers are distinct modulo pi. The parameter prescription below is a terminating rational search followed by explicit error inequalities. Its termination is proved for every n. It is not an extrapolation of the n=1,2,3 checks.

[proved, provenance] The local two-bus construction and its two geometric scales are C39's construction, read from ../work_C39/PROOFS.md. This note supplies a detailed proof of the singular-pencil eigenvector limit, an exact acceptance test avoiding asymptotic thresholds, and explicit finite-D and finite-g budgets. The implementation files are attributed snapshots of C39, with C42 modifications identified in RESULT.md. The original 23-letter all-n problem remains unresolved here.

[assumed, infrared scope] An interacting central charge at fixed nonzero microscopic couplings requires local elimination and Gaussian stability, including the mobile mixed vertices. Those are not consequences of the finite-band theorem. The conclusion here is the payload input requested in C42, with the vacancy compatibility described in Section 8.

## 1. Local definition and exact code

[proved, definition] The letters are L_a,R_a (four), d_(v,a) (two), A_empty,B_empty,A_v,B_v (four), and, for each bus r=1,2, d_(r,s,a),x_(r,s,a),A_(r,s),B_(r,s), where a=0,1 and s=0,+,- (eighteen per bus). Thus d_core=46 for every n. A legal initialized word has A first, B last and n intervening data sites, with exactly one token. Data to its left are L, data to its right are R. Bulk arrows are L->L, R->R, and L->R for a token. A_empty exits L, A with a token exits R; B_empty enters R, B with a token enters L. Internal pairs may not have A on the right or B on the left. Mismatched arrows and wrong boundary letters cost K>0.

[proved, definition] Vacuum motion follows A,d_0,...,d_(n-1),B. Each directed forward vacuum edge is the positive operator

    (q|initial>-|final>)(q<initial|-<final|),  0<q<1.

An excited token of bus r follows the tree

    A -- d_0 -- d_1 -- ... -- d_(n-1)
         |
         x_0 -- x_1 -- ... -- x_(n-1) -- B.

The d-to-d hopping magnitude is w_r, where w_1=w and w_2=wv. All other edge magnitudes are one. The junction switch d_0<->x_0 is conditioned on A on its left. Every rightward neutral departure from a data site carries Z of that site's bit. Each of the n neutral express hops, including x_(n-1)->B, also carries exp(i*pi/(2n)). Charged hops and the A-to-d_0 link are real. There is no loop link.

[proved, definition] The excited diagonal is a=4 on d and a+D on A,x,B. Conversion at a data token is

    g |r,+,0><r,0,1| + g |r,-,1><r,0,0| + h.c.

The integer onsite charge is the data bit plus q_s, with q_+=1,q_-=-1 and q_0=q_v=0. All moves conserve it. In particular the total U(1) and its parity subgroup are exact.

[proved, positive local implementation] Put Q_left+2Q_right on each internal bond, where Q counts an excited token, and add the excited off-diagonal moves on that bond. Each left token has at most one move of magnitude at most one; each right token has at most two. This matrix is diagonally dominant with nonnegative diagonal, hence positive. This remains true on illegal pair states and with two buses, because a token has only its own bus's moves. Add the onsite operator (a-3)Q+D Q_high+g C. Its excited eigenvalues are at least 1-|g|. Split it into halves on incident bonds, and restore the missing halves at the two boundaries. Also restore 2Q at A and Q at B. These boundary terms are positive. The resulting initialized Hamiltonian has the stated diagonals and is a sum of positive local terms for |g|<1. The vacuum edges and grammar projectors are separately positive.

[proved] Kinetic moves preserve the outer arrows and marker positions, so violation sectors are invariant and cost at least K. In a legal zero state the strictly positive excited onsite term excludes every excited token. For each data string the remaining connected vacuum path has kernel

    G|b> = S^(-1/2) sum_(ell=0)^(n+1) q^ell |ell,v;b>,
    S=sum_(ell=0)^(n+1) q^(2ell).

There are exactly 2^n independent zero vectors. The path's nonzero eigenvalues are 1+q^2-2q cos(j*pi/(n+2)), j=1,...,n+1. To see this, its incidence matrix has rows q e_j-e_(j+1); the row Gram is the constant tridiagonal matrix with diagonal 1+q^2 and off-diagonal -q. Consequently the complete initialized core gap is at least

    gamma=min(K,1-|g|,1+q^2-2q cos(pi/(n+2)))>0.

Exact degeneracy refers to this unperturbed core. The weak field below is a retained perturbation, as in the original candidate.

## 2. Ports and the effective bond

[proved, definition] On endpoint e=A or B and bus r put

    F_(r,0)=|r,0><v|+|v><r,0|,
    F_(r,+)=|r,+><v|+|v><r,-|,
    F_(r,-)=F_(r,+)^dagger.

The port perturbation lambda V uses

    V=sum_r b_r [F_(r,0)^B F_(r,0)^A
           +exp(i phi_r) F_(r,+)^B F_(r,-)^A
           +exp(-i phi_r) F_(r,-)^B F_(r,+)^A],
    b_1=1, b_2=sqrt(eta), phi_1=pi/2, phi_2=pi.

This conserves total charge and has zero first-order compression. Two different spatial bonds cannot contribute at second order: an outer cell touched once is left excited. Different buses cannot interfere at second order either, because unperturbed propagation preserves the bus label and the final source must remove the same label.

[proved, tree gauge] Multiplying a neutral token at data position j by Z_0...Z_(j-1) removes departure strings. Removing express phases then makes every species' scalar graph the same real matrix H_D(w_r). Conversion becomes g c_j or g c_j^dagger, with the usual Jordan-Wigner c_j. At A the neutral source is I and at B it is -iP, P=product Z_j. Charged sources are I. These identities follow on every occupation vector, so they are full Fock-space identities, not one-particle fits.

[proved] Let H_r include the three excited species and g on one cell. Let B_r be the two-cell source map formed from endpoint insertion and coefficients 1,exp(i phi_r),exp(-i phi_r). Then the exact coefficient of lambda^2 is

    K2=-V0 sum_r b_r^2 B_r^dagger
                        (H_r tensor I+I tensor H_r)^(-1) B_r,
    V0=q^(2n+2)/S^2.

Here ||B_r||=sqrt(3). This formula retains all conversion orders and the complete code. On three cells the coefficient is K2_12+K2_23.

[proved, leading conversion] At order g^2, two conversions in the same cell give a scalar. The real scalar kernel is symmetric in their sites, and the + and - channels combine into {c_i,c_j^dagger}=delta_ij. One conversion in each cell gives hopping. The two allowed charge paths and the right neutral parity factor give 2i exp(i phi_r), so the relative hopping phase of bus 2 is i. This is a relative virtual-channel phase, which cannot be removed by one data-basis gauge.

[proved, Gram limit] Put A(w)=4I-w adjacency(P_n), K(t)=exp(-t A(w)), and

    f_j(t;w)=integral_0^t K_(0,j)(s) K_(j,0)(t-s) ds,
    Gamma(w)_ij=integral_0^infinity f_i(t;w) f_j(t;w) dt.

Writing p=2n+4, the hopping divided by J0=2 V0 g^2 D^(-p) tends, at order g^2, to

    T0=-Gamma(w)-i eta Gamma(wv).

A left endpoint propagator reaches the stub through one high vertex and a right endpoint propagator through n+1 high vertices. Each filter contains two such propagators. More explicitly their Laplace transforms contain the squares of 1/(D+4+z) and 1/det(C_D+z), respectively, times the squared junction-to-data resolvent. C_D is the n+1 vertex express block. The junction Schur correction tends to zero. The two high transfer kernels, after multiplication by D^2 and D^(2n+2), converge to unit point masses. Their normalized versions are convolutions of positive exponential densities. The stub heat kernels are bounded by a fixed decaying exponential at w<=1/4. Dominated convolution convergence proves the displayed Gram limit. Equivalently the finite pair inverse is rational in 1/D, and its first nonzero hopping coefficient is -2 Gamma at degree p. Section 6 bounds its remainder.

[proved, exact computability] Gamma is rational for rational w. For example let R=(A tensor I+I tensor A)^(-1), z=R|0,0>. Then

    Gamma_ij=2 sum_(x,y) z_(i,y) R_((i,y),(x,j)) z_(x,j).

This follows by the two possible time orderings of the insertions on the two tensor factors. It supplies a finite rational definition and the implementation in certify.py.

[proved, positivity and analytic normalization] The first nonzero small-t term of f_j is w^(2j)t^(2j+1)/(2j+1)!. Distinct lowest powers prove linear independence and Gamma(w)>0. Put D_w=diag(w^(2j)) and C(w)=D_w^(-1) Gamma(w) D_w^(-1). Walk expansions, integrated against exp(-8t), show that C is analytic at w=0 and

    C(0)_ij=(2i+2j+2)!/
         [(2i+1)!(2j+1)! 8^(2i+2j+3)].

This is a positive monomial Gram matrix. Analyticity also follows from the rational inverse above after canceling the powers w^(2i+2j).

## 3. The second geometric scale

[proved] Let R_read count a data bit only when its own token is the vacuum token. It is a repeated onsite diagonal operator. Its compression is

    G^dagger R_read G=(q^2/S) sum_j q^(2j) n_j.

Choose q=(wu)^2. Add lambda^2 f R_read with f=J0 mu S/q^2. Its compressed one-particle field, divided by lambda^2 J0, is mu D_w E_u D_w, E_u=diag(u^(4j)). The readout need not preserve the code. Its leakage affects higher orders in lambda, not this coefficient or the exact kernel of the unperturbed core.

[proved] The limiting reference symbol, after congruence by D_w^(-1), is

    h(k)=-2 A cos k+2 eta B(v) sin k+mu E_u,
    A=C(w), B(v)=D_v C(wv) D_v, D_v=diag(v^(2j)).       (1)

Congruence preserves zeros and their multiplicities. It is not being used to identify energies or to discard cross entries. Both full positive matrices remain in (1).

## 4. All-n separation of centers and widths

[proved, eigenvalue scales] Fix n and w=1/4. Let r_0(v)>=...>=r_(n-1)(v)>0 be the generalized eigenvalues B(v)x=r A x. For sufficiently small v, positive constants c,C bound C(wv) between cI and CI. Bounds on A and the min-max principle therefore give

    (c/||A||) v^(4j) <= r_j(v)
        <= (C/lambda_min(A)) v^(4j).

In particular all r_j are distinct and r_(j+1)/r_j tends to zero. This argument retains the off-diagonal Gram correlations.

[proved, leading coefficients] Expanding det(B-tA), the least power of v in the coefficient of (-t)^(n-k) selects rows and columns {0,...,k-1} from B. Every other selection has strictly larger power. Thus the elementary symmetric polynomial of order k is

    c_k v^(2k(k-1))(1+o(1)),
    c_k=det C(0)[0:k,0:k] det A[k:n,k:n]/det A, c_0=1.

The scale bounds above imply that the product of the k largest eigenvalues dominates this polynomial. Taking successive quotients proves

    r_j(v)=kappa_j v^(4j)(1+o(1)), kappa_j=c_(j+1)/c_j>0. (2)

[proved, eigenvectors with all tail entries retained] Normalize x_j by its j-th coordinate. Then

    x_j(v) -> x_j^0,
    (x_j^0)_i=0 (i<j), (x_j^0)_j=1,
    (A x_j^0)_i=0 (i>j).                                (3)

Here is a block proof of this otherwise singular limit. Split P={0,...,j-1}, T={j,...,n-1}. In the P equations factor D_P on both sides. Since r_j D_P^(-1) A_PP D_P^(-1)=O(v^4), the resulting P block has a bounded inverse. Eliminating it gives x_P=O(v^2)||x_T||. Its Schur complement on T, divided by v^(4j), tends to

    b_j e_0 e_0^T-kappa_j A_TT,
    b_j=det C(0)[0:j+1]/det C(0)[0:j]>0.

To check the discarded terms, the leading B complement is D_T times the Schur complement of C(wv) times D_T. After division by v^(4j), only its (j,j) entry survives. The effective A differs from A_TT by O(v^2), since D_P^(-1) C_PP^(-1) C_PT D_T=O(v^2). Higher terms contain another r_j D_P^(-2)=O(v^4). The limiting kernel is one-dimensional and proportional to A_TT^(-1)e_0, whose first coordinate is positive. Normalizing that coordinate gives (3). This also proves that normalization is valid for small v.

[proved, centers] At mu=0 the zeros are

    p_(j,+)=theta_j+pi/2, p_(j,-)=theta_j-pi/2,
    theta_j=-arctan(eta r_j).

Choose positive eta sufficiently small, and then v sufficiently small, so that |theta_0|<1/8 and |theta_(j+1)|<|theta_j|/32. Every root is simple. In the generalized eigenbasis its scalar factor is -2 cos k+2 eta r_j sin k, with derivative magnitude 2 sqrt(1+eta^2 r_j^2) at the zero.

[proved, width derivatives] Differentiating the kernel equation at mu=0 gives

    dp_(j,+)/dmu=-s_j, dp_(j,-)/dmu=+s_j,
    s_j=(x_j^T E_u x_j)/
        [2 sqrt(1+eta^2 r_j^2) x_j^T A x_j]>0.           (4)

At v=0, (3) makes the numerator u^(4j)(1+O(u^4)); its denominator has a finite strictly positive limit, with r_j->0 for j>0 and r_0->kappa_0. First choose u sufficiently small that these limiting positive slopes obey s_(j+1)^0<s_j^0/64. Then choose v sufficiently small that the actual slopes obey s_(j+1)<s_j/32, while retaining the center bounds. Choosing v before u without controlling its leakage would not justify this step.

[proved, all quartic comparisons] A difference of two unordered endpoint pairs is described by integer counts h_(j,s), sum h=0 and sum |h|<=4. At mu=0 its momentum is

    sum_j m_j theta_j+(pi/2) sum_(j,s) s h_(j,s),
    m_j=h_(j,+)+h_(j,-).

The second term is an integer multiple of pi. The first has magnitude at most 4|theta_0|<1/2. Odd multiples of pi stay away from zero modulo 2pi. For even multiples and m!=0, the earliest nonzero m_j dominates: all later terms together are at most 4|theta_(j+1)|, so the absolute value is at least (7/8)|theta_j| and hence (7/8)|theta_(n-1)|. If m=0, the comparison has the same flavor multiset. Different endpoint pairs then have nonzero coefficients -sum_s s h_(j,s) in the mu derivative, with total absolute value at most four. The same dominance bound applied to (4) is at least (7/8)s_(n-1). This includes the n=1 umklapp comparison of two identical right points with two identical left points. Repeated labels are deliberately included, even when Fermi statistics would remove their local vertex.

[proved, finite mu and margins] Each of the finitely many roots is analytic in mu near zero. Let L bound its first derivative and R its second derivative on a sufficiently small common interval. Choose mu>0 so that 4L mu<|theta_(n-1)|/8 and 2R mu<s_(n-1)/8. Nonzero constant terms retain a margin; zero constant terms have a linear term larger than the remainder. A conservative margin is

    delta4 >= min(|theta_(n-1)|/4, mu s_(n-1)/4)>0.       (5)

Endpoint differences are covered by comparing pairs with one common label. Sea centers are covered by comparing opposite-endpoint pairs. This proves all-n nonidentity. No independence of the individual flavor parameters is assumed.

## 5. A terminating rational schedule and its quantitative certificate

[proved, algorithm] The following is an explicit schedule in n, although not a practical complexity bound. Fix w=1/4. Enumerate dyadic eta,u,v,mu in (0,1) by increasing sum of their positive exponents. Dovetail this with increasing rational interval precision. For each tuple calculate Gamma(w),Gamma(wv) exactly and use the physical, uncongruent reference

    H0(k)=-2 Gamma(w) cos k+2 eta Gamma(wv) sin k
                    +mu D_w E_u D_w.                  (6)

Also require mu<t_*, where t_*=det Gamma(w)/(tr Gamma(w))^(n-1). This makes H0(0) negative definite and H0(pi) positive definite. Use x=tan(k/2). Multiplying (6) by 1+x^2 gives a real symmetric rational quadratic matrix polynomial. Its determinant P(x) is a real rational polynomial of degree at most 2n. Accept only if it has 2n simple real roots in (-2,2), each isolated in a disjoint rational interval I_i=[l_i,r_i], and the pair tests below are strict. Sturm division and gcd(P,P') supply exact root counts and simplicity. Dyadic roots are handled as exact singleton roots, or by rational endpoints avoiding them. The concrete small-n program uses witnesses with no dyadic roots.

[proved, pair acceptance without floating-point angles] For each interval midpoint m_i put z_i=(1+i m_i)/(1-i m_i), a rational point of the unit circle. A root in I_i moves its phase by at most r_i-l_i, since |d(2 atan x)/dx|<=2. Use the larger safe error e_i=2(r_i-l_i). For every two different unordered pairs (i,j),(k,l), compute

    d_(ij,kl)=max(|Re(z_i z_j-z_k z_l)|,
                     |Im(z_i z_j-z_k z_l)|)-e_i-e_j-e_k-e_l.

Accept if all d_(ij,kl)>0. Let Delta be their minimum. The chord lower bound is also an angular lower bound. The same list bounds endpoint separation by comparing (i,i) with (i,j). It proves delta2>=Delta, delta4>=Delta and center separation modulo pi>=Delta/2. The checks use integers and rational arithmetic only. An angle-interval implementation with rigorously bounded arctangents, such as certify.py, is an equivalent acceptance test.

[proved, termination] Section 4 supplies an open set with simple roots and all strict inequalities. Its roots are within a small neighborhood of +/-pi/2, hence x is in (-2,2). Sequentially sufficiently small dyadic eta,u,v,mu enter that set. Rational root isolation and pair enclosures eventually resolve every nonzero separation. Dovetailing therefore accepts a tuple for every n. There is no claim that fixed formulas v=1/1024,u=1/64 work at arbitrary n; those are only the verified small-n witnesses.

[proved, explicit perturbation budget from the accepted certificate] Use the physical matrices (6). At each of the 4n rational endpoints x_e calculate H_e exactly. Put

    M_e=max(1, n max_(i,j)(|Re H_e,ij|+|Im H_e,ij|)),
    beta_e=|det H_e|/M_e^(n-1)>0,
    t_*=det Gamma(w)/(tr Gamma(w))^(n-1)>0,
    tau=min(1,t_*/8,min_e beta_e/4).                    (7)

All are positive rationals. Since ||H_e||<=M_e, its smallest singular value is at least beta_e. Any Hermitian nearest-neighbor symbol within tau of (6) at these endpoints preserves their determinant signs. Each bracket still has a zero, since P changes sign there. A nearest-neighbor n-by-n determinant has degree at most 2n in exp(ik) after clearing negative powers. The 2n distinct bracketed zeros therefore exhaust that degree, are simple, and have no additional competitors. They retain the same pair and center margins Delta and Delta/2. Since mu<t_* and tau<=t_*/8, the global error budgets below also preserve negative definiteness at k=0 and positive definiteness at k=pi. Each semicircle has n crossings in the same inertia direction. In particular these are n pairs of left and right Fermi endpoints. This argument does not extrapolate a momentum grid or assume derivatives of a small perturbation.

[proved, velocity bound] Also impose ||T-T0||<t_*/2, ensured by the budgets below. Then -Herm T>=t_* I/2, so |det T|>=(t_*/2)^n. Set

    B=1+2 tr Gamma(w)+2 eta tr Gamma(wv)+mu.

The corrected symbol has norm at most B+1. The unit-circle root polynomial has leading coefficient det T. At a root its derivative has magnitude |det T| times the product of the other 2n-1 chord distances, each at least Delta. Dividing by the product of the n-1 nonzero eigenvalues of H gives

    |v_F| >= (t_*/2)^n Delta^(2n-1)/(B+1)^(n-1).         (8)

This is in J0 units. Multiplication by lambda^2 J0 gives physical speeds. The argument includes n=1 with the empty denominator product.

## 6. Explicit finite-express and finite-conversion errors

[proved, analytic disk] For a single bus at g=0 write the two-carrier Hamiltonian as z^(-1)N+H_f, z=1/D. N counts high vertices and is at least one on the high subspace. Each one-cell finite part has norm at most 7 by its row sums, so ||H_f||<=14<16. The low-low block is at least 7I at w<=1/4. Put

    r=1/(64*16^2)=1/16384, p=2n+4, B0=128.

On |z|<=r, the high inverse has norm <=2r. Its correction to the low block is <=2r*16^2=1/32. The block inverse formula bounds the complete inverse by

    (1+32r)^2/(7-512r)+2r < 1.

It extends analytically to z=0. Since ||C_pair||<=2 and ||B_r||=sqrt(3), each entry of the g^2 source coefficient has modulus <=12<B0 on that circle.

[proved, Cauchy remainder] The hopping coefficient vanishes below degree p and its degree-p term is -2 Gamma(w), with the appropriate bus phase, by Section 2. Cauchy's estimate therefore gives, for 0<z<=r/2,

    ||D^p T2_r/2 - T0_r|| <= n B0/(D r^(p+1)).

T2_r here has its endpoint visibility removed; its phase is included. The two buses, with weights 1 and eta, contribute a Bloch-symbol error at most

    E_D=2n B0(1+eta)/(D r^(p+1)).

There is no nonscalar onsite term at g^2. Choose an integer

    D >= max(2/r, 8n B0(1+eta)/(tau r^(p+1))).            (9)

Then E_D<=tau/4. This explicitly pays for all narrow scales through tau.

[proved, nonzero conversion] Each scalar excited-cell matrix is at least I by row sums at a=4. Its pair matrix is at least 2I and ||C_pair||<=2. Source compression eliminates odd conversion orders, since one conversion changes the parity of the number of charged tokens whereas both source states have an even number. For 0<g<=1/4 the full remainder beginning with g^4 is bounded conservatively by

    ||K2-K2(0)-g^2 K2^(2)|| <=6(1+eta)V0 g^4.

For example summing the inverse Neumann series from order four gives 3(2g)^4/[2^5(1-g)]<=2g^4 for each bus, which is stronger. Dividing by J0 gives epsilon_bond<=3(1+eta)g^2 D^p. In the Hilbert-Schmidt quadratic projection each bilinear coefficient is bounded by four times the operator-norm error. Summing onsite left/right and both hopping row bounds gives a Bloch error at most 16n epsilon_bond. Thus choose a positive dyadic g with

    g<=1/4,
    g^2<=tau/[192n(1+eta)D^p].                          (10)

Then the finite-g symbol error is <=tau/4. Equations (9),(10) together preserve all brackets, simple roots, centers, margins and the velocity bound. The exact second-order quadratic projection is nearest-neighbor, so Section 5's degree argument applies at finite g.

[proved, non-Gaussian remainder] The same full-code norm bound controls all non-Gaussian terms. There are at most (2n)^2 number-conserving bilinears. Bounding their coefficient norms as above, including the centered scalar, gives the safe estimate

    ||non-Gaussian part||/J0
          <=(2+64n^2) 3(1+eta)g^2 D^p.                 (11)

One may reduce g further so this is below any prescribed positive interaction tolerance, including one defined using (8) and Delta. Reducing g keeps the reference field tied to J0 and does not reduce the reference momentum margins. This fixes the loop-versus-conversion ordering problem identified in C38.

[proved, finite clusters] All microscopic parameters are now finite. The readout coefficient f is finite, although it can be very large. On a fixed cluster with unperturbed gap gamma, the resolvent on the circle of radius gamma/2 has a convergent Neumann series when |lambda| ||V_cluster||+lambda^2 |f| ||R_cluster||<gamma/4. The alternate-cell excited-token parity sends lambda to -lambda and acts trivially on the code, while leaving the readout invariant. The canonical effective Hamiltonian thus has no odd powers and its remainder after the displayed lambda^2 coefficient is O(lambda^4). Its constants are finite and can be bounded by that resolvent series. This is not a uniform-in-volume block-diagonalization claim.

## 7. Verification and limitations of the calculations

[numerical, n=1,2,3] The literal 46-letter onsite/two-site matrices are compared with the tree-gauge matrices in literal.py. The pair matrix has dimension 2116. Code dimensions are 2,4,8 and complete legal-cell dimensions are 54,160,424. Matrix agreement is at floating-point precision. Positivity and exact charge conservation are proved in Sections 1 and 2, not inferred from these decimals.

[proved, computer-assisted at n=1,2,3] certify.py constructs exact rational Gram matrices and uses a rational Sturm chain to isolate every tangent root. Alternating-series arctangent bounds and a rational Machin bound for pi certify every pair mismatch, including modular wraparound. These are certificates for the limiting references. Equations (9),(10), not a large-D floating-point scan, transfer them to finite microscopic parameters. The all-n statement instead follows from Section 4 and the terminating search in Section 5.

[numerical, specified sectors] ed.py diagonalizes the untruncated invariant legal sectors of two and three cells at n=1,2,3 and compares their low-energy matrices to the full-in-g coefficient K2, including the readout. RESULT.md lists dimensions, residuals and errors. Numerical spectra are not rigorous finite-dimensional proofs. The finite-dimensional identities they test were proved above. At the moderate test parameters the narrowest hoppings are much smaller than scalar shifts, so absolute ED agreement does not certify relative accuracy for every tiny hopping entry. The direct g^2 source contraction is a separate check of the Gram normalization and phase.

## 8. Compatibility with the vacancy completion

[proved, initialized core input] This core has a complete positive initialized kernel, positive finite gap, local endpoints, exact charge and silent ports. These are the finite hypotheses of the existing C8/C32 fixed-alphabet wrapper. Its endpoint vacuum weights are unequal but positive and known, q^0/sqrt(S) and q^(n+1)/sqrt(S). Endpoint occupancy compensation acts separately at A and B; it multiplies these bare amplitudes rather than requiring them to be equal. The readout moves with occupied labels and is independent of their physical spacing.

[proved, using the inherited wrapper theorem] Substituting 46 for C32's 28-letter core gives d_wrap=96*46+48=4464. No new n-dependent symbol is introduced. C32's full second-order resolvent estimate applies to the sum of the two bus sources, with their norms included. It makes every length-pair coefficient arbitrarily close to K2, scalar part included. One can take its marker-motion coefficient zero and detune its selector to a single length before vacancy doping. Explicitly, at the inherited b/(b+1) coexistence point let G_sel=t_sel eta_b-U>0 be the wrong-length bound and lower the selector marker chemical potential by epsilon with 0<epsilon<min(mu_sel/2,b G_sel/2). Relative to the new b-cell energy density, an l-cell cost changes by epsilon(l-b)/b. The b+1 cost becomes epsilon/b>0; the inherited G_sel(1+|l-b|) bound keeps every other length positive. The marker chemical potential stays positive, and the original fragment and marker-free exclusions remain valid. On rings divisible by b the unperturbed selected code therefore consists only of b-cells. This assertion invokes the existing finite wrapper and selector theorem, not a new implementation or numerical audit of its 4464-letter bond. The bare initialized dimension 46 does not include a spontaneous selector, and a 47-letter full selected chain is not claimed.

[proved, vacancy algebra, independent of the wrapper] For any resulting repeated internal bond h on C^d, form

    h_vac=-t sum_a (|vac,a><a,vac|+h.c.)
          +J (iota tensor iota) h (iota tensor iota)^dagger
          -(mu_occ/2)(n_occ tensor I+I tensor n_occ).

This has dimension d+1, exact occupied number, exact internal U(1), and parity. The whole h, including its half-onsite terms and the readout, is lifted into occupied pairs. Ordered occupied coordinates preserve internal label order; crossing the ring cut gives the cyclic label permutation. Thus the exact squeezing identity of the specified vacancy completion applies without any special property of the original 23 letters. Its leading thermodynamic internal coefficient is the same positive scalar alpha(rho)=rho-sin(pi rho)^2/(pi^2 rho) multiplying every term of h. Scalar multiplication leaves all payload momenta and their mismatch margins unchanged. Using the inherited selected internal wrapper gives a fixed, conservative alphabet d=4465 after adding the vacancy.

[assumed, finite-coupling mobile conclusion] The gapless charge sector cannot be eliminated using its finite-ring level spacing as a uniform denominator. Finite-J corrections, fractional charge harmonics tied to marker order, physical momentum gluing, and interacting dressing still need the mobile infrared analysis. This note proves compatibility of the local construction and leading payload target with that completion. It does not infer c=n+1, uniqueness on every ring, or a counterexample to a central-charge bound merely from squeezing.
