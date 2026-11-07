# C135 proofs

## 1. Conventions and the LTI dual

[Proved, conventions] Trace distance here means the unhalved norm. Write

    D_n(X) = X tensor I - I tensor X,
    M(h) = inf_(a,o) ||h+aI+o tensor I-I tensor o||.

The scalar changes both energy densities equally. The divergence has zero expectation in an LTI state and sums to zero on a ring. C132 section 1 proves that the dual of M on differences with zero trace and equal partial traces is the trace norm. No vacuum-compatible length restriction is used in this note.

[Proved, F1] If a ring state is invariant under cyclic translation, its n-site consecutive marginal has equal first and last (n-1)-site marginals, for every 2 <= n <= L, including n=L. The same statement follows from translation invariance on Z. Averaging a ring ground state over translations preserves its energy. Consequently

    E_LTI(n) <= e0,
    E_LTI(n) <= E_L/L  (L >= n).

Tracing an end site of an LTI state gives an LTI state. Indeed, trace another end in Tr_1 rho_n=Tr_n rho_n; commutation of the two endpoint traces gives the new equality. Iterating shows that all translated consecutive k-site marginals agree wherever they fit inside the n-site interval. Thus E_LTI(n) is nondecreasing.

[Proved, F2] For the primal, use the equality A(rho)=Tr_n rho-Tr_1 rho=0. Its adjoint is A*(X)=D_n(X). The Lagrangian for minimizing Tr(H rho), with H=h tensor I_(n-2), is

    Tr[(H+D_n(X)-tI)rho]+t.

Its infimum over rho>=0 is t if the coefficient is positive, and minus infinity otherwise. The maximally mixed state is strictly feasible in the primal. A sufficiently negative t with X=0 is strictly feasible in the dual. Finite-dimensional SDP duality therefore gives attainment and

    E_LTI(n)=max_(X=X*,t) {t: H+D_n(X)>=tI}.

Embed each inequality on the ordered cyclic windows of a ring with L>=n. Every translated X occurs once with each sign, including when L=n. The sum is H_L>=LtI. This proves the telescoping certificate without an assumption about saturation.

[Proved, F3] Put n=L in F1 and in an assumed uniform bound e0-E_LTI(n)<=KM/n². Then E_L/L>=e0-KM/L². The marginal inclusion (K) follows from C132's norm identification and convex separation. The implication preserves K. A directional lower bound E_LTI(n)-E_L/L is valid for every n, but can be positive only when n>L. The condition n>L in the computational brief is necessary for a useful positive witness, not for validity of the inequality.

[Proved, terminology] A ground-space mixture exposes a face of Omega_L. Its marginal need not itself be an exposed extreme point. This distinction does not affect the directional energy identity. Finite increases in a numerical sequence cannot alone refute existence of a finite universal K. Unbounded growth requires an infinite family with controlled bounds.

## 2. Periodization and supports

[Proved, F4] Let rho_n be LTI and b=rho_(n-1) its common endpoint marginal. Put P=supp(b). Positivity implies

    rho_n=(I tensor P)rho_n(I tensor P)
          =(P tensor I)rho_n(P tensor I).

For example, Tr[(I tensor (I-P))rho_n]=0. For positive operators, zero expectation of a projection implies that the projection annihilates the state on both sides. This proves the support statements, including vanishing cross blocks.

Identify rho_n as the target on the shifted memory plus a new site. Its Petz map on P is

    R(Y)=sqrt(rho_n)(b^(-1/2)Y b^(-1/2) tensor I)sqrt(rho_n),

with inverse on the support. It is CP, and

    Tr R(Y)=Tr[Y b^(-1/2)(Tr_last rho_n)b^(-1/2)]=Tr(PY),
    R(b)=sqrt(rho_n)(P tensor I)sqrt(rho_n)=rho_n.

For any output state tau, the completion

    R_full(Y)=R(PYP)+Tr[(I-P)Y]tau

is CP and trace preserving. It agrees with the support map on every input encountered below, including when that input is entangled with earlier sites. No full-rank assumption is required. The same construction works for any CPTP R with R(b)=rho_n. A replacement channel proves that the feasible set is nonempty.

[Proved] Start with rho_n on sites 1 through n. Repeatedly apply R to the last n-1 sites, appending one site. The last n sites always have state rho_n, since their input marginal is b. At the final update touching the bond (j,j+1), the enclosing n-site input is rho_n, and its frozen marginal is

    sigma2 = Tr_(B2...B_(n-1),C) (id_A tensor R)(rho_n).

Later updates act strictly to its right. Every fixed finite initial window eventually stabilizes, so a consistent half-line state exists. Translate arbitrarily long interior windows to Z and average translates. Compactness and a diagonal subsequence give a TI state, and every bond in it has marginal sigma2. Therefore sigma2 belongs to Omega_inf and e0<=Tr(h sigma2).

[Proved] Its two one-site marginals equal rho_1, by trace preservation on A and R(b)=rho_n on the output. The gauge drops out of its energy difference from rho_2. The finite states before stabilization are not a projectively consistent sequence, which is why the stabilization argument is needed. This proof uses LTI alone, not a ring extension.

## 3. Reflection and twists

[Proved] Let J_n reverse the order of n tensor factors. Then

    Tr_1(J_n rho J_n)=J_(n-1)(Tr_n rho)J_(n-1),
    Tr_n(J_n rho J_n)=J_(n-1)(Tr_1 rho)J_(n-1).

Thus reversal preserves LTI. Its bond energy is Tr[(FhF)rho_2], using equality of all internal bond marginals. For FhF=h, averaging rho with its reversal preserves feasibility and objective. If h is real, averaging with complex conjugation also preserves them. These operations commute, so a real reversal-invariant optimizer exists. On-site symmetry twirling likewise preserves LTI and the objective. For U(1), one may use fixed-magnetization blocks. These are symmetries of the full primal, not assumptions about the optimizer.

[Proved, definition] Suppose g is unitary and [h,g tensor g]=0. Let T send the tensor factors (1,...,L) to (2,...,L,1), and let U=g_1 T, with the convention that the factor arriving at site 1 is conjugated by g. A symmetry-twisted ring is

    H_L(g)=sum_(j=1)^(L-1) h_(j,j+1)
             +(I_L tensor g_1)h_(L,1)(I_L tensor g_1*) .

Here h_(L,1) is ordered with site L first. Direct conjugation shows that U cycles all L terms. For the step taking the seam to bond (1,2), use the common on-site symmetry. Hence [H_L(g),U]=0 and U^L=g^(tensor L).

[Proved] Choose a U-invariant ground state, for example the normalized ground-space projector, or average any ground state under conjugation by U. The compact closure of the powers of U supplies the averaging if g has infinite order. Every bulk window in the fixed linear order 1,...,L is LTI. To see this, shift an observable supported on sites 1,...,L-1 one site right. Its shifted support misses site 1, so the extra conjugation by g has no effect. U-invariance equates its expectations on the two consecutive windows. All bulk h terms and the twisted seam have the same expectation, so

    Tr(h rho_2)=E_L(g)/L,
    E_LTI(n)<=E_L(g)/L  for n<=L.

[Proved, qualification] It is false that every ordinary cyclic window of every twisted-ring ground state is LTI. First, a ground state in a degenerate ground space need not be U-invariant. Second, a window crossing the seam needs the corresponding g transport to compare it with a bulk window. Ordinary cyclic permutation does not include this transport.

An explicit example is XX with field b=1/2, L=3 and g=Z. The unique ground vector in the sector with one up spin is

    |w>=(|100>-|010>+|001>)/sqrt(3),

where 1 denotes the up spin in this displayed occupation convention. Its energy is -5/4 and it is the unique full ground state, as follows from the two magnetization-sector spectra. In the linear order its neighboring coherences are both -1/3, so it is LTI_3. A cyclic permutation gives amplitudes proportional to (1,1,-1); its two adjacent coherences are +1/3 and -1/3, so it is not LTI_3. Twisted transport repairs this. The unqualified window interpretation in MAP section 7 needs this correction.

## 4. Ring saturation is false, even with symmetry twists

[Proved] Use the real reflection-symmetric bond

    h_b=(XX+YY)/4+b(ZI+IZ)/4,   b=1/2.

Its four eigenvalues are -1/2,-b/2,b/2,1/2. Its gauge-optimized scale is M=1/2. Indeed, the spread of h is 1, and swapping the two sites changes the sign of every divergence. Convexity of the centered operator norm, followed by this swap average, shows that a divergence cannot improve the scale of a swap-symmetric bond.

[Proved, exact n=3] The staggered W vector of section 3 has both two-site marginals equal. Its hopping energy per bond is -1/3 and its magnetization is -1/6. Thus its energy is -5/12. A matching dual is

    X=-(XX+YY)/8-ZI/12-IZ/24.

Substitution gives

    h_12+D_3(X)=(h_XX,12+h_XX,23)/2+(Z1+4Z2+Z3)/24.

Its spectrum is

    {-5/12,-1/3,-1/4,-1/6,1/6,1/4,1/3,5/12}.

For verification, the one-up sector is the tridiagonal matrix with diagonal (-1/6,1/12,-1/6) and adjacent off-diagonal entries 1/4. Its antisymmetric eigenvalue is -1/6; the symmetric 2 by 2 block has eigenvalues -5/12,1/3. Spin reversal gives the opposite sector's eigenvalues, and the two polarized states give +/-1/4. Therefore E_LTI(3)=-5/12 exactly.

[Proved, all twisted rings L>=4] Every on-site symmetry of this bond preserves Z, since the nonzero partial trace of h is proportional to Z. Up to scalar phase, g is therefore a z rotation. These twists are arbitrary fermionic fluxes. At fixed occupation m, the most negative possible hopping sum over L equally spaced momenta is

    -sin(pi m/L)/sin(pi/L).

This follows by summing m consecutive unit-circle points and centering their sum on the negative real axis. The minimizing occupied momenta are consecutive; a disconnected occupied set can be improved by replacing a point farther from pi by an unoccupied closer point. Flux can realize the centering for either Jordan-Wigner parity. Hence the exact minimum over all symmetry twists is

    min_g E_L(g)/L
      =min_(0<=m<=L) {-sin(pi m/L)/(L sin(pi/L))+b(m/L-1/2)}.

For b=1/2 and L=4,5,6,7 the minima are, respectively,

    -3/8,
    -(1+sqrt(5))/10-1/20,
    -sqrt(3)/6-1/12,
    -2 cos(pi/7)/7-3/28.

The last three exceed -3/8. These finite comparisons can also be checked directly in the formula over the L+1 values of m. Particle-hole reflection excludes m>L/2 for b>0.

[Proved, tail bound] Filling the infinite-line dispersion cos(k)+b gives

    e0= -[sqrt(1-b²)+b arcsin(b)]/pi
       =-sqrt(3)/(2pi)-1/12  at b=1/2.

The momentum-box proof of C132 section 8 applies to an arbitrary flux grid in this number-conserving chain. It preserves the on-site occupation and multiplies hopping by sinc(pi/L), giving

    E_L(g)/L >= e0-pi²M/(6L²).

For L>=8 its right side is at least

    -sqrt(3)/(2pi)-1/12-pi²/768 > -3/8.

For example, substituting the rational bounds 3.14159<pi<3.14160 and 1.73205<sqrt(3)<1.73206 proves the strict inequality with a margin greater than 0.003. This is a bound on the entire tail, not an extrapolation. Thus the minimum over e0 and every symmetry-twisted ring with L>=4 is exactly -3/8. The untwisted L=4 ring attains it. The untwisted L=3 ring has minimum -1/4. It follows that the proposed untwisted saturation at n=3 fails by exactly 1/24, or 1/12 in units of M.

[Proved, exact n=4 feasible point] Consider the real Toeplitz occupation covariance

    C=(1/10000) [[3187,-2870,2007,-756],
                [-2870,3187,-2870,2007],
                [2007,-2870,3187,-2870],
                [-756,2007,-2870,3187]].

The leading principal minors of 10000C are

    3187, 1920069, 93784040, 2524309676.

Those of 10000(I-C) are

    6813, 38180069, 143496795960, 7082071749676.

All are positive integers, so Sylvester's criterion gives 0<C<I. There is a number-conserving Gaussian state with this covariance. One construction diagonalizes C and puts independent occupation weights equal to its eigenvalues in the corresponding fermionic modes. This defines a positive density operator on four spins through Jordan-Wigner. Use the convention c_j=(product_(k<j)(-Z_k))S_j^-, so c_j* c_j=(1+Z_j)/2 and the hopping has the positive sign above.

The two three-site restrictions of C are identical. Gaussian reduced states are determined by their covariance. Their corresponding spin reduced states are also identical: every even spin monomial on an interval is a local even CAR monomial under the interval's Jordan-Wigner identification; odd expectations vanish by parity. Therefore the four-spin state is LTI. Its energy is the exact rational number

    C_12+(1/2)(C_11-1/2)=-7553/20000=-0.37765 < -3/8.

This refutes saturation by rings with symmetry twists as well, with a certified feasible gap at least 53/20000=0.00265 in energy, or 0.0053 in units of M. It also refutes the untwisted statement at n=4. The state is real and reflection invariant as a spin state: alternatively average with spin reversal of site order, which preserves its energy and LTI by section 3, without needing any convention for reflection of fermions.

[Numerical] The optimized n=4 value is -0.37837433475, corresponding to a gap about 0.00337433475 from the best twisted ring, or 0.00674866950 in units of M. The rational witness, rather than this decimal optimum, proves failure. free_covariance.py reproduces the small covariance SDP and checks all principal minors with integer arithmetic.

[Proved, closed finite-dimensional formulation] For the whole XXZ family h_Delta+b(ZI+IZ)/4, the n=3 and n=4 SDPs can be written without a full unconstrained density matrix. Use blocks R_k indexed by subsets A of {1,...,n} with |A|=k. Their sizes are (1,3,3,1) and (1,4,6,4,1). Set rho_(A,B)=R_|A|(A,B) when |A|=|B| and zero otherwise. Impose R_k>=0, sum_k Tr R_k=1 and, for all equal-cardinality subsets A,B of {1,...,n-1},

    sum_(s=0,1) rho_((s,A),(s,B))
       =sum_(s=0,1) rho_((A,s),(B,s)).

Reversal averaging permits R_k(A,B)=R_k(reverse A,reverse B), and reality permits real symmetric blocks. The objective is the sum of the diagonal coefficients Delta z1 z2/4+b(z1+z2)/4, plus the real matrix entries connecting configurations 01 and 10 on the first two sites, with coefficient 1/2 in each direction. This is a closed linear SDP for arbitrary Delta,b. Its dual is the blockwise inequality of section 1, with X on 2 or 3 sites. No closed formula for its optimum at all Delta,b is asserted. The n=3 example and the rational n=4 witness settle the universal conjecture without that formula.

[Assumed, remaining conjecture] Exact saturation, with or without twists, should be discarded. A defensible weaker question is whether the difference between the best twisted-ring energy density at lengths >=n and E_LTI(n) is uniformly O(M/n²), together with a uniform twisted-ring bound. Neither estimate is proved here for interacting bonds. The direct LTI rate remains a cleaner sufficient conjecture; generic bonds may have only scalar on-site symmetries, so arbitrary adjustable fluxes cannot be inserted without changing the bond.

## 5. A uniform LTI rate for quadratic qubit bonds

[Proved, theorem] Suppose, after a common on-site change of basis and a gauge change, h is a nearest-neighbor even quadratic Jordan-Wigner bond. This includes arbitrary XY anisotropy, transverse field and z-directed DM, with no ZZ interaction. For every n>=2,

    0 <= e0(h)-E_LTI(n)
       <= M(h)[1-cos(pi/(n+1))]
       <= pi² M(h)/(2(n+1)²).

In particular (L-LTI) holds on this entire class with K=pi²/2, uniformly in gaps, fields and crossover scales. This is an LTI statement, stronger than the previously established ring statement on this class. It does not extend to arbitrary interacting U(1)-symmetric bonds.

[Proved, covariance construction] Average an LTI state under total fermion parity. This preserves its energy and LTI. Its generalized one-particle covariance G_n, in sitewise Nambu order, satisfies 0<=G_n<=I and the CAR particle-hole identity. These inequalities follow by applying positivity to arbitrary linear combinations of creation and annihilation operators and using their anticommutators. Its 2 by 2 site blocks are Toeplitz. Indeed, the corresponding bilinears are even Jordan-Wigner strings confined between their endpoints, and all translated strings with that support have equal expectations by LTI.

Set

    v_j=sqrt(2/(n+1)) sin(pi j/(n+1)),  1<=j<=n,
    s=sum_(j=1)^(n-1) v_j v_(j+1)=cos(pi/(n+1)).

Normalization and the value of s follow from the normalized lowest sine eigenvector of the path adjacency matrix with off-diagonal entries 1/2. For each k, let V(k) have site blocks v_j exp(ikj) I_2. Then V(k)*V(k)=I_2, and

    G_infinity(k)=V(k)*G_n V(k)

satisfies 0<=G_infinity(k)<=I_2. Because the v_j are real, the particle-hole identity of G_n gives the required particle-hole relation between G_infinity(k) and G_infinity(-k). Thus this is the symbol of a TI even Gaussian CAR state. It induces a TI spin state: define local even spin moments by the interval Jordan-Wigner map and set odd moments to zero. Positivity and compatibility follow from the positive Gaussian reduced states on each interval, and translation invariance follows from the Toeplitz covariance.

[Proved, existence detail] The Gaussian existence step can be checked on finite intervals. Convert a physical generalized covariance to its real antisymmetric Majorana covariance. An orthogonal change of Majorana basis puts it into 2 by 2 blocks with parameters nu_j in [-1,1]; this interval is equivalent to 0<=G<=I. In the associated canonical modes take a product of diagonal one-mode density matrices with eigenvalues (1+nu_j)/2 and (1-nu_j)/2, choosing the order to match the covariance convention. Transforming back by the finite Clifford representation gives a positive even state with the prescribed covariance. Its higher even moments obey Wick's rule, which makes the state unique within the Gaussian class and makes restrictions Gaussian. For the symbol above, Fourier transformation gives a bounded positive contraction on the infinite one-particle space, so every finite interval restriction is physical. The finite Gaussian states are consistent because they have the same restricted covariance and Wick moments. Their consistent limit gives the required infinite CAR state. The interval spin states are also consistent: Jordan-Wigner strings outside an interval cancel in even products, while odd expectations vanish. This justifies the spin extension without an infinite Jordan-Wigner string operator.

[Proved] If G_r denotes the r-th Toeplitz block of G_n, Fourier integration gives

    (G_infinity)_r = G_r sum_(j=1)^(n-|r|) v_j v_(j+|r|),   |r|<n.

In particular G_0 is unchanged and G_1 is multiplied by s. Therefore the energy of the constructed infinite state is

    e_trial=e_on+s e_cross,
    e=Tr(h rho_2)=e_on+e_cross.

Scalar terms are included in e_on. Since e0<=e_trial and e<=e0 for an LTI optimizer,

    0<=e0-e <= (s-1)e_cross <= (1-s)||h_cross||.

To compare with M, parity-average a minimizing centered gauge representative. The one-site gauge can be chosen even. Let P be local parity. Cross-site quadratic terms are odd on each site, whereas the on-site and gauge terms are even. Thus

    h_cross=[h_tilde-(P tensor I)h_tilde(P tensor I)]/2,
    ||h_cross||<=||h_tilde||=M.

This proves the theorem. Attainment of the finite-dimensional gauge minimum can alternatively be replaced by an arbitrarily close minimizing sequence. The inequality 1-cos x<=x²/2 proves the displayed constant.

[Proved, checks] A pure on-site bond has e_cross=0 and zero error. Multiplication by a positive scalar rescales both sides equally. The construction is valid for n=2, where s=1/2. No massless approximation is made, so near-critical TFIM at g=1+c/n is included uniformly. Pairing is covered by the generalized covariance and its particle-hole condition, not by treating anomalous covariances as independent matrices. The rational four-site counterexample is consistent with the theorem: exact ring saturation is unnecessary for the rate.

## 6. An explicit XX telescoping certificate

[Proved] Let h_j=(X_j X_(j+1)+Y_j Y_(j+1))/4 on an n-site interval. Keep v and s from section 5 and set

    a_j=v_j v_(j+1)/s,  sum_(j=1)^(n-1) a_j=1,
    x_j=sum_(k=1)^j a_k-1,  1<=j<=n-2,
    X_n=sum_(j=1)^(n-2) x_j h_j  on n-1 sites.

The empty sum is allowed at n=2. Coefficient comparison gives

    h_1+D_n(X_n)=sum_(j=1)^(n-1) a_j h_j.

For the first bond the coefficient is 1+x_1=a_1, for an interior bond it is x_j-x_(j-1)=a_j, and for the last it is -x_(n-2)=a_(n-1). This verifies telescoping explicitly.

[Proved, positivity] Define normalized annihilation modes A(k)=sum_j v_j exp(ikj)c_j. CAR gives {A(k),A(k)*}=1, so 0<=A(k)*A(k)<=I. Fourier integration gives

    sum_j v_j v_(j+1) h_j
      = integral_(-pi)^pi cos(k) A(k)*A(k) dk/(2pi).

For cos(k)>=0 the integrand is positive. For cos(k)<0 it is bounded below by cos(k)I. Consequently

    sum_j v_j v_(j+1) h_j >= -I/pi,
    h_1+D_n(X_n) >= -I/[pi cos(pi/(n+1))].

This is an explicit positive decomposition: subtract -I/pi before division by s, and the integrand is cos(k)A*A on the positive half of the dispersion and |cos(k)|(I-A*A) on the negative half. It is a continuous sum of positive quadratic operators. The same weighted XX inequality already appears in C105 section 3, proved there using D_v |T| D_v. The new use here is the normalized sine choice and conversion into the LTI telescoping certificate. No historical-priority claim is made.

[Proved, explicit constant] Since e0=-1/pi and M=1/2, the deficit of this certificate is (sec(theta)-1)/pi, theta=pi/(n+1). Put m=n+1>=3. The inequality

    cos(pi/m)>=1-pi²/(2m²)>=(1-1/m)²

follows from 2m-1>=5>pi²/2. Hence

    n²(sec(theta)-1)
      <= n² theta²/[2 cos(theta)] <= pi²/2.

It follows that

    h_1+D_n(X_n) >= [e0-pi M/n²]I.

Thus c=pi is explicit for XX at every n>=2. The numerical optimum is generally better. This certificate consists only of translated nearest-neighbor bonds. Its proof of positivity uses extended fermionic modes, which are multi-site Jordan-Wigner strings as spin operators. Multi-site terms are therefore not necessary in X for a Casimir-rate XX certificate, even if a solver returns them in an optimal X.

[Proved, relation to lesson 2] This is a quadratic-chain inequality controlling the full finite-size deficit, not an extraction of the conformal anomaly from a local deformation. It makes no claim that the same weighted positivity estimate survives interactions. C105 already gives an interacting obstruction to a naive replacement of its XX coefficient by an interacting energy density. The proof uses CAR positivity and the quadratic dispersion, not a CFT local response expansion.

## 7. Complementary slackness does not supply the recovery estimate

[Proved] Let t=E_LTI(n), S=h_1+D_n(X)-tI>=0 and rho be an optimal primal. Strong duality gives Tr(S rho)=0. Positivity implies S rho=rho S=0, by applying the vanishing trace to sqrt(rho) S sqrt(rho). This is the stated identity.

For any n-site state eta obtained after recovery and discarding the appended site,

    Tr[h_1(eta-rho)]=Tr(S eta)-Tr[D_n(X)eta].

The first term is nonnegative but has no available small upper bound. The second is a boundary functional of eta and need not vanish because eta need not be LTI. If instead eta is an n-site marginal of the final TI periodized state, its divergence expectation vanishes and

    Tr(h sigma2)-t=Tr(S eta)>=0.

There is still no upper estimate. Slackness controls S on the support of rho, while the recovered state may have weight outside that kernel. Neither the magnitude of that weight nor its S-energy is controlled by slackness.

[Proved obstruction to a channel-independent inference] The admissible replacement map R(Y)=Tr(Y)rho_n freezes the product rho_1 tensor rho_1. For antiferromagnetic Heisenberg its energy is zero in a spin-invariant optimizer, whereas t approaches 1/4-log(2)<0. This yields an order-one defect despite exact complementary slackness. Thus the identity plus admissibility cannot bound every channel. This is not a disproof for the Petz map or the energy-minimizing channel.

[Proved subclasses] For a diagonal classical LTI state, use the order-(n-1) classical conditional extension. It leaves the memory symbols unchanged, appends with the conditional weights of rho_n, and preserves the old bond. Its defect is zero. If a bondwise minimizing infinite state exists, then E_LTI(n)=e0=lambda_min(h) for every n, so the energy rate has K=0, irrespective of recovery quality. For the symmetric ferromagnetic Heisenberg mixture rho_n=P_sym,n/(n+1), C132's support calculation applies at ell=n and gives optimal and Petz defect 1/[2n(n+1)], hence n² delta/M=n/(n+1). This particular ansatz need not have zero defect even on a frustration-free chain.

[Proved, limitation] Section 5 closes the LTI energy rate for the quadratic class using covariance smoothing. It does not prove the same estimate for Petz or for a channel constrained to send rho_(n-1) to rho_n. The smoothing state need not have that prescribed n-site marginal. No new recovery estimate is proved for general interacting U(1) bonds.

[Proved audit] C132 section 5 correctly notes that the conditional-entropy continuity argument retains a logarithm of n. It does not rule out a uniform trace-norm bound when crossover scales are unrestricted. Moreover, full-window recovery error and two-site marginal distance are different quantities. The obstruction established here is the absence of an upper estimate from the slackness identity, not a theorem that trace-norm methods cannot work.

## 8. Literature and the remaining lemma

[Numerical, literature] Kull, Schuch, Dive and Navascues report approximately n^(-2) LTI error in their critical benchmarks, with parity dependence and faster behavior for gapped examples. Their compressed relaxations reach larger effective sizes; this is not a proved uniform rate in h. Sections II.6 and V.2, Figures 2 and 5 were opened and read. [Source](https://arxiv.org/html/2212.03014v3).

[Proved, literature and elementary specialization] Fawzi, Fawzi and Scalet Appendix B gives convergence by tiling. In one dimension, translate-average a product of copies of an LTI_n state. Of every n bonds, n-1 have energy e and one has the product-seam energy. Centering the optimal gauge bounds their difference by 2M, giving e0-e<=2M/n. This does not give the rate sought here. [Source](https://arxiv.org/html/2305.06855v2#A2).

[Proved, literature] Blakaj and Wolf prove that the infinite-chain two-site marginal set is not semialgebraic for d>=2. A fixed finite LTI level is a spectrahedral shadow and is semialgebraic, so it cannot equal Omega_inf. This excludes finite-level exactness for all bonds, not an n^(-2) approximation rate or the particular free-chain construction above. Theorem 2 and its XY-family argument were opened and read. [Source](https://arxiv.org/html/2308.14585v2).

[Assumed, single remaining sufficient lemma] There exists a finite K, independent of the interacting qubit bond and n, such that every h admits a Hermitian X on n-1 sites satisfying

    h tensor I+D_n(X)>=[e0(h)-K M(h)/n²]I.

This is proved here for quadratic bonds, and trivially for bondwise frustration-free energy problems and classical diagonal bonds. It remains open for general interacting qubit bonds. Ring saturation cannot replace it, because section 4 refutes both the untwisted and symmetry-twisted versions at d=2. No finite numerical scan proves or disproves this remaining uniform statement.
