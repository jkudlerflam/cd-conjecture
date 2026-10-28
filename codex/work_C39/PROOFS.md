# C39 construction and proofs

[proved, scope] The main finite result below is a two-bus modification of C32. It uses a biased vacuum walk to read out a geometric density profile. Every assertion labeled proved refers to finite matrices or their analytic limits. The mobile wrapper is inherited explicitly. A thermodynamic central charge is not proved here.

## 1. Local core and alphabet

[proved, definition] Start with C32's initialized word A, n data sites, B. Retain L_a,R_a, d_(v,a), A_empty,B_empty,A_v,B_v. For each r=1,2 add d_(r,s,a),x_(r,s,a),A_(r,s),B_(r,s), where s=0,+,- and a=0,1. There are 10+2*18=46 letters. The arrows, endpoint conditions and penalties are those of C32 PROOFS Section 1, with each excited species duplicated. The charge is data occupation plus q_s, with q_0=q_v=0 and q_+=1,q_-=-1. The token type is vacuum or one of the two buses.

[proved, definition] There is no weak loop edge. The excited graph of either bus is

    A -- d_0 -- d_1 -- ... -- d_(n-1)
         |
         x_0 -- x_1 -- ... -- x_(n-1) -- B.

Its edges have magnitude one except for d_j--d_(j+1), whose magnitude is w_r. Set w_1=w, w_2=wv, with 0<w,v<1. The physical moves are C32's local moves. Every neutral departure from a data site carries Z on its data bit. Neutral express hops, including x_(n-1)--B, carry exp(i*pi/(2n)); there are n of these hops. All other phases vanish. The vacuum uses the path A,d_0,...,d_(n-1),B and never the express rail.

[proved, definition] Replace each forward vacuum edge by the positive matrix

    (q|initial>-|final>)(q<initial|-<final|),  0<q<1.       (1)

Choose uniform excited diagonal a=4 on all data-rail vertices, and a+D on A, B and x. The two conversions at each data site in bus r are

    g |d_(r,+,0)><d_(r,0,1)|
      +g |d_(r,-,1)><d_(r,0,0)| + h.c.                    (2)

[proved, local positive decomposition] These diagonal choices have a repeated positive implementation. On each ordered internal bond put Q on its left site and 2Q on its right, plus all excited off-diagonal moves on that bond. Here Q counts any excited token. The resulting matrix is positive by diagonal dominance. A token on the left of a bond has at most one allowed move on that bond. A token on the right has at most two, the A--d and d--x moves at the junction. Each magnitude is at most one. This also holds on illegal pair states and for the two-bus duplication because only the token's own bus acts. More than one token only increases the allocated diagonal. Equivalently subtract each edge's endpoint diagonals from Q_left+2Q_right, leaving a nonnegative diagonal, and add its positive two-state edge matrix.

[proved] Add the positive onsite block (a-3)Q+D Q_high+g C, where C is (2) with g removed. Its least excited eigenvalue is at least a-3-|g|. Split this onsite term over incident bonds. At A supply the missing 2Q and onsite half as a positive boundary term, and at B supply Q and the other half. This gives diagonal a or a+D at every excited vertex. Add positive grammar and boundary penalties K. Thus the entire initialized core is a sum of positive one- and two-site terms when |g|<1. No negative compensation term is being hidden in the core.

[proved, complete code] Positivity of the onsite blocks forces every zero vector to have no excited tokens. The grammar then leaves exactly one vacuum token on the connected path. Equation (1) forces its amplitudes to have ratio q along that path. Consequently, with k=n+2 and S=sum_(ell=0)^(n+1) q^(2ell),

    G|b> = S^(-1/2) sum_(ell=0)^(n+1) q^ell |ell,v;b>,
    dim ker H_core = 2^n.                                (3)

This proves completeness, including illegal words, not only a code embedding in a selected block. The weighted path has nonzero spectrum
1+q^2-2q cos(r*pi/k), r=1,...,k-1. This follows from the nonzero singular values of the incidence matrix with rows q e_j-e_(j+1), whose row Gram is the constant tridiagonal matrix with diagonal 1+q^2 and off-diagonal -q. Therefore

    gap_core >= min(K,1-|g|,1+q^2-2q cos(pi/k))>0.         (4)

[proved] Every move and conversion conserves the onsite charge. The code charge is the number of occupied data bits. This is exact U(1), with its exact fermion-parity subgroup. It is independent of the perturbative approximations below.

## 2. Sources, phase control and the exact second-order bond

[proved, definition] At endpoint e=A,B, use the original source maps separately for each bus,

    F_(r,0)=|r,0><v|+|v><r,0|,
    F_(r,+)=|r,+><v|+|v><r,-|, F_(r,-)=F_(r,+)^dagger.

The two-cell perturbation is

    V=sum_r b_r [F_(r,0)^B F_(r,0)^A
       +exp(i phi_r) F_(r,+)^B F_(r,-)^A
       +exp(-i phi_r) F_(r,-)^B F_(r,+)^A],               (5)

with b_1=1,b_2=sqrt(eta), phi_1=pi/2,phi_2=pi and eta>0. This is Hermitian and charge conserving. Its first-order code compression vanishes.

[proved, gauge] Remove the neutral departure strings and express phases on the excited tree. The data-rail conversion then becomes g c_j in the +,0 block and g c_j^dagger in the -,0 block. Here c_j=(product_(i<j) Z_i)|0><1|_j. Every species has the same scalar graph Hamiltonian K_r. The neutral left source is identity, the neutral right source is -iP, with P=product_j Z_j. The charged sources are identity on the data. The phase -i follows from n express hops. This gauge is only a calculation of the physical local rules above.

[proved, exact formula] Let H_r be this complete excited-cell Hamiltonian, including g. Its dimension is 3(2n+2)2^n. Define S_(L,s),S_(R,s) by inserting the token at A,B, with the preceding data matrices, and

    B_r=S_(R,0) tensor S_(L,0)
       +exp(i phi_r) S_(R,+) tensor S_(L,-)
       +exp(-i phi_r) S_(R,-) tensor S_(L,+).

The full second-order coefficient, before any expansion in g or D, is

    K2= - q^(2n+2)/S^2 sum_r b_r^2
          B_r^dagger (H_r tensor I+I tensor H_r)^(-1) B_r. (6)

Different buses have no cross term: unperturbed propagation preserves the bus of each excited token, and a different bus's endpoint operator annihilates it. Two different spatial bonds also have no cross term, since an outer cell remains excited after its sole port insertion. Thus the exact three-cell coefficient is K2_12+K2_23. This assertion holds on the full code, including all particle sectors and all non-Gaussian terms in (6).

[proved] The g^2 coefficient is quadratic plus a scalar. Two conversions in one cell sum into the fermion anticommutator, since + and - have identical scalar propagation. They leave no nonscalar onsite term. One conversion in each cell gives hopping. In the real tree gauge the two interfering charge paths have the same port phase exp(i phi_r), and their neutral-source factors give 2i. The hopping is therefore a real matrix times i exp(i phi_r). In particular the second bus multiplies the first bus's phase by i. This phase is a physical difference between virtual channels, not a change of basis on the data.

[proved, large-D limit] Put A(w)=aI-w adjacency(P_n), K(t)=exp(-t A(w)), and

    f_j(t;w)=int_0^t K_(0,j)(s) K_(j,0)(t-s) ds,
    Gamma(w)_ij=int_0^infinity f_i(t;w)f_j(t;w) dt.       (7)

Eliminating the high vertices gives the leading hopping

    T= -J [Gamma(w)+i eta Gamma(wv)]+o(J),
    J=2 lambda^2 g^2 [q^(2n+2)/S^2] D^(-2n-4).           (8)

The little-o can be taken first as D tends to infinity at the g^2 coefficient and then as g tends to zero relative to the finite D scale. A left source must traverse A and return, giving two high factors. A right source must traverse B and n express vertices and return, giving 2n+2 high factors. The remaining low propagators are (7). This is the same Schur expansion as C32, now with uniform stub diagonal and variable stub hopping. No return path of smaller high order reaches the logical data. The phase computation above fixes the sign and the factor two.

[proved, computability] Gamma(w) is a rational matrix when w,a are rational. One definition requiring only finite linear algebra is C32's Lyapunov construction

    Bcal=[[I tensor A,-I tensor |0><0|],[0,I tensor A]],
    z=(0,vec I), Bcal Y+Y Bcal^T=z z^T.

Take the diagonal-coordinate entries of Y in its first block. The real parts of the eigenvalues of Bcal are positive, so this equation has a unique solution. The resulting matrix is (7).

[proved] The shortest walk gives

    f_j(t;w)=w^(2j)t^(2j+1)/(2j+1)!+O(t^(2j+2)).

The first nonzero powers differ. Thus the functions are linearly independent and Gamma(w) is positive definite for every w>0. Define D_w=diag(w^(2j)) and C(w)=D_w^(-1) Gamma(w) D_w^(-1). It extends analytically and positively to w=0, where

    C(0)_ij=(2i+2j+2)! /
      [(2i+1)!(2j+1)!(2a)^(2i+2j+3)].                   (9)

This is the Gram matrix of exp(-at)t^(2j+1)/(2j+1)!. Analytic extension follows either from the convergent walk expansion at |w|<a/2 or the finite rational Lyapunov equation after canceling the displayed powers.

## 3. A second geometric scale supplied by the code

[proved] The repeated onsite readout R counts a data bit only when its site's token is the vacuum token. Its exact compression is

    G^dagger R G = (q^2/S) sum_(j=0)^(n-1) q^(2j) n_j.   (10)

Choose an additional rational u in (0,1) and q=(wu)^2. Add the retained perturbation

    H_field= J mu S/q^2 sum_sites R_site.                 (11)

Then its compressed single-particle matrix is J mu D_w E_u D_w, where E_u=diag(u^(4j)). The unperturbed code is still the frustration-free code (3). The readout generally leaks out of that code; it is an order-lambda^2 perturbation. Its leakage contributes only at higher orders in lambda and is not silently included in the zero-kernel hypothesis.

[proved] After dividing by J and making a positive congruence by D_w^(-1), the limiting Bloch matrix is

    h(k)= -2A cos k+2 eta B(v) sin k+mu E_u,
    A=C(w)>0,
    B(v)=D_v C(wv) D_v>0.                               (12)

A congruence preserves the zeros and their inertia, not the energies or velocities. The two different Gram matrices remain in (12). Thus the argument does not infer Fermi-point separation from a bandwidth hierarchy.

## 4. The generalized eigenvalue lemma

[proved] Fix n, rational w in (0,1/2), and A=C(w). Number the generalized eigenvalues of B(v)x=r A x in decreasing order r_0,...,r_(n-1). Set

    c_k=det C(0)[0:k,0:k] det A[k:n,k:n]/det A,
    c_0=1, k=0,...,n,
    kappa_j=c_(j+1)/c_j>0.

Empty determinants are one. Then, as v tends to zero through positive values,

    r_j(v)=kappa_j v^(4j)(1+o(1)).                       (13)

[proved, derivation] Expand det(B-tA). In the coefficient of t^(n-k), choose k rows and columns from B. Its entries carry factors v^(2i+2j). The unique least exponent is attained by using row and column sets {0,...,k-1}; that exponent is 2k(k-1). Its coefficient is det C(0)[0:k] det A[k:n], with the determinant sign already accounted for by (-t)^(n-k). Every other minor has a strictly larger exponent. Dividing by det A proves that the k-th elementary symmetric polynomial of the generalized eigenvalues is c_k v^(2k(k-1))(1+o(1)). Positivity and the Newton polygon, or induction on these elementary symmetric polynomials, give the n different valuations 0,4,...,4n-4 and coefficients (13). For the inductive version, the first k largest eigenvalues dominate the k-th symmetric polynomial once the preceding valuation is separated; the next coefficient gives their successive quotient. Equivalently substitute t=v^(4j)y into the characteristic polynomial. The exponent of its k-th term is 2k(k-1)+4j(n-k), whose two minima are k=j and k=j+1. After division by the common power the limiting polynomial is (-y)^(n-j-1)[-c_j y+c_(j+1)], up to a common sign. Its nonzero root kappa_j is simple. The implicit-function theorem supplies that root. There are n such disjoint scales and n roots in total.

[proved, eigenvectors] Normalize x_j(v) so that its j-th coordinate is one. It tends to the vector x_j^0 with

    (x_j^0)_i=0 for i<j,
    (A x_j^0)_i=0 for i>j,
    (x_j^0)_j=1.                                        (14)

Thus its tail is A[j:n,j:n]^(-1)e_0 divided by that inverse's first entry. Here is the scaled elimination in detail. Partition the first j coordinates as I and the tail as T. For j=0 there is no elimination. For j>0,

    B_II=D_I C_II(wv) D_I,
    ||B_II^(-1)||=O(v^(-4j+4)),
    ||B_II^(-1)B_IT||=O(v^2),
    r_j ||B_II^(-1)||=O(v^4).

The middle estimate follows entrywise from D_I^(-1) C_II^(-1) C_IT D_T; every exponent is at least two. Hence the first block equation gives x_I=O(v^2)||x_T||, and its Schur complement is

    D_T [C_TT-C_TI C_II^(-1)C_IT] D_T
       -r_j A_TT+o(r_j).

The terms omitted here contain either B_II^(-1)B_IT, its transpose, or r_j B_II^(-1), so the displayed estimates bound them in norm. Divide by v^(4j). The first matrix tends to c_* e_0 e_0^T, where c_*=det C(0)[0:j+1]/det C(0)[0:j]>0. Thus the limiting tail equation is

    c_* e_0 e_0^T x_T=kappa_j A_TT x_T.

It forces x_T to be proportional to A_TT^(-1)e_0. Its first entry is nonzero, so normalization by that entry gives (14). The identity kappa_j=c_* (A_TT^(-1))_00 follows from the displayed determinant ratios. All normalized subsequences have this same limit. Invertibility on the complementary directions supplies the continuous, in fact analytic after the scaling, eigenvector branch. This proves the eigenvector limit without setting v=0 in the unscaled singular pencil.

## 5. Centers, widths and all pair sums

[proved] Fix eta>0 so eta*kappa_0<1/8. At mu=0 the roots are

    p_(j,+)=theta_j+pi/2, p_(j,-)=theta_j-pi/2,
    theta_j=-arctan(eta r_j).                            (15)

They are simple for small positive v, since the r_j are distinct and the derivative on x_j is nonzero. Shrink v so all |theta_j|<1/4 and |theta_(j+1)|<=|theta_j|/32. Equation (13) proves this finite choice exists.

[proved] The derivatives with respect to mu are

    dp_(j,+)/dmu=-s_j, dp_(j,-)/dmu=+s_j,
    s_j=(x_j^T E_u x_j)/
           [2 sqrt(1+eta^2 r_j^2) x_j^T A x_j]>0.        (16)

This is the differentiated zero-eigenvalue equation, using cos p=eta r_j sin p. At v=0 the numerator on x_j^0 is u^(4j)(1+O(u^4)), and the denominator tends to a strictly positive constant. Hence one can first choose a sufficiently small rational u so the limiting slopes obey s_(j+1)^0<s_j^0/64, and then choose v sufficiently small so the actual slopes obey s_(j+1)<s_j/32. This order is essential. The q in (3) is fixed after u, and before D,g,lambda. The case n=1 has no ratio constraints.

[proved, elementary lacunarity bound] If positive numbers z_j obey z_(j+1)<=z_j/32, then any nonzero integer combination with sum_j |m_j|<=4 has magnitude at least (7/8) z_last. Indeed the least index carrying a nonzero coefficient contributes at least z_j, and all remaining coefficients together contribute at most 4z_(j+1)<=z_j/8. This bound applies to the absolute centers, which all have one sign, and to the slopes.

[proved, pair-sum classification] Compare any two different unordered pairs of Fermi-point labels, allowing repeated labels. Their difference at mu=0 is

    sum_j m_j theta_j + (pi/2) sum_(j,s) s h_(j,s),

where the signed label counts h have total zero and l1 norm at most four, and m_j=h_(j,+)+h_(j,-). The second term is an integer multiple of pi. The center term has magnitude below one, so an odd multiple of pi cannot resonate modulo 2pi. If m is nonzero and the multiple is even, the lacunarity bound is at least (7/8)|theta_(n-1)|. If m=0, the comparison uses the same flavor multiset. Different label pairs then have a nonzero slope combination

    -sum_(j,s) s h_(j,s) s_j,

whose coefficient l1 norm is at most four. Its magnitude is at least (7/8)s_(n-1). Thus the only zero-order coincidences split linearly at nonzero mu. This covers interflavor pair transfer, equal-width backscattering, repeated-point comparisons and umklapp. Quadratic root differences follow either by the same proof or by adding one common label to both pairs.

[proved, velocity check] At mu=0 the physical, J-normalized Bloch matrix is D_w h(k) D_w. Its zero vector is D_w^(-1)x_j. The magnitude of its energy derivative is

    2 sqrt(1+eta^2 r_j^2) x_j^T A x_j /
                         ||D_w^(-1)x_j||^2 >0.

The finite minimum over j is positive. Small mu and the later small matrix perturbations preserve a positive minimum. Congruence has therefore not been used to assign unphysical velocities.

[proved, finite nonzero field] All 2n roots are analytic in mu on some interval around zero. Let L bound their first derivatives and R bound their second derivatives there. Both bounds are finite and computable from the finite polynomial pencil. Choose a rational mu>0 small enough that the roots remain simple, 4L mu<|theta_(n-1)|/8, and 2R mu<s_(n-1)/8. Comparisons with nonzero constant term retain a positive margin. Those with zero constant term have a linear term exceeding their quadratic remainder. A conservative lower bound is

    delta4 >= min(|theta_(n-1)|/4, mu s_(n-1)/4),
    delta2 >= delta4.                                   (17)

The small interval can also ensure the sea centers differ by at least one quarter of the smallest zero-field center separation. Their first-order shift vanishes by (16), but even a first-derivative bound suffices. No claim of a uniform-in-n margin is made.

[proved, effective parameter algorithm] Fix w=1/4. Compute A and C(0) by rational linear algebra. Choose rational eta<min(1,1/(8kappa_0)). Compute the vectors (14), then search dyadic u for the strict limiting slope inequalities above. Their positive asymptotic coefficients prove termination. Next search dyadic v for separated generalized eigenvalues, the center ratios and slope ratios. Equation (13) and (14) prove termination. Finally search dyadic mu for simple roots and strictly distinct pair sums. The preceding classification proves termination for all sufficiently small mu. Polynomial root isolation and rational interval arithmetic can perform every test. Replacing exp(ik) by z makes det h(k) a degree-2n Laurent polynomial with rational complex coefficients. Simplicity, unit-circle root counts, pair products and their separation therefore have finite algebraic certificates. This is a constructive terminating search, not an inference from small-n numerical scans. Its running time is not bounded by minutes for arbitrary n; no such computation is required to state the all-n family.

## 6. Returning to finite microscopic parameters

[proved] At fixed n,w,u,v,eta,mu the reference (12) has nonzero velocities and strictly positive margins. Its root and pair-sum conditions persist in a sufficiently small C1 neighborhood of the Bloch matrix. The scale D_w is invertible, so this is also an open condition in the original coordinates. A tolerance can be computed from isolated root intervals, minimum slopes and the minimum absolute determinant outside those intervals. The field and hopping are held fixed while that tolerance is chosen.

[proved] The finite-D coefficient is rational in D and analytic in z=1/D near zero after multiplying its hopping by D^(2n+4). Schur inversion gives its limit (8). The scalar onsite terms are removed only after computing the full bond; no nonscalar onsite g^2 term remains. Choose an integer D large enough to put the normalized g^2 hopping within one quarter of the preceding tolerance. The convergent finite-dimensional rational expansion proves termination. One may certify the tail using a complex resolvent disk as in C32 PROOFS Section 6, or compute the rational coefficient with interval enclosures. Both use a strictly positive low block a-2w and a high block D+a with degree at most three.

[proved] At this finite D the Neumann expansion in g converges for |g|<1. The source denominator is at least 2(1-|g|). Each conversion changes a neutral to a charged token, so only even total conversion counts return to the source. The non-Gaussian remainder starts at g^4. A crude uniform bound for one bus is 6g^4/(1-|g|) after the source visibility is removed, for |g|<=1/4; using a larger constant is harmless. Therefore, after dividing by the nonzero scale 2g^2 D^(-2n-4), its remainder tends to zero with g. Choose rational g sufficiently small for the normalized quadratic remainder and full non-Gaussian remainder to lie in the desired finite tolerance. This can cost powers of D depending on n, but is finite at every fixed n.

[proved, finite-cluster qualification] The readout (11), including its prefactor S/q^2, is finite at fixed n. Take lambda last. The exact local code gap (4) and bounded port/readout operators give a convergent finite-cluster effective expansion at sufficiently small lambda. On an open chain, alternate-cell excitation parity sends lambda to -lambda while fixing the readout. Thus odd powers of lambda vanish. Readout leakage, mixed readout-port processes and higher port terms are bounded higher-order errors, which can be made arbitrarily small relative to J by this final choice. This is a finite-cluster statement and does not prove a uniform thermodynamic block diagonalization.

## 7. Frozen and mobile wrappers

[assumed inherited finite construction] Apply C8's complete initialized-core wrapper, occupation filter, selector and separate endpoint compensators, as used in C32. The required hypotheses are positivity, an initialized complete code, local endpoint ports and a finite core gap, all established in Sections 1 and 2. The charge of every extra wrapper track is zero. The port amplitudes at A and B are unequal but known from (3). The wrapper compensates physical endpoint occupation separately; it need not equate the two bare core amplitudes.

[proved alphabet arithmetic under that inherited construction] The full onsite space is

    C^4_selector tensor
      [(C^3_tag tensor C^2_filter tensor C^46_core) direct-sum C^3_hole]
      tensor C^2_helixA tensor C^2_helixB,
    d0=96*46+48=4464.                                   (18)

The frozen implementation can use the same alphabet with marker hopping zero and the selector in a single-length window. The mobile implementation selects physical lengths n+3,n+4 with the same occupied core of n+2 sites and the same n-mode payload. The number 46 is the initialized-core alphabet, not a claim that the spontaneous selector has already been included in 46.

[proved conditional on the inherited wrapper bounds] C32's noncommuting second-order estimate applies to the entire two-bus source map, including its scalar part. Choose the vacancy ratio sufficiently small, then the core scale divided by the actual fast gap sufficiently small, to approach (6) in every length pair. Constants now include both ports and the biased core gap. The field readout moves with occupied letters; its sum in occupied coordinates is (10) for every vacancy pattern. It needs no new flavor addressing. At finite wrapper ratios the residual bond error is retained and bounded, not claimed to vanish exactly. The known positive core allows the same proof of identical payload transport by marker motion.

[proved, mixed margins conditional on cell-coordinate translation] Once the payload points in (17) are fixed, only finitely many length densities are excluded by the norm-four and norm-eight local vertices of C20, with sum of payload charges zero imposed by exact U(1). Choose an interior density away from those finitely many points and pure length commensurabilities. An interval with positive mixed margin exists. C8/C27's identity-transport gluing is the inherited microscopic translation input. The finite bounds can be made smaller than the combined payload and mixed margins in the same ordered hierarchy as C32.

[assumed thermodynamic input] Local elimination of massive sectors and weak-coupling Gaussian stability are still needed to conclude c=n in the frozen sector or c=n+1 in the mobile model at fixed nonzero coupling. The finite construction does not prove these statements, uniqueness on every finite ring, or a common Lorentz velocity. They are separate from the all-n finite-band construction.


## 8. Exact small-n certificates and numerical checks

[proved, rational certificate] certify.py computes Gamma using a second formula independent of the Lyapunov implementation. Put K=A(w) tensor I+I tensor A(w), R=K^(-1), z=R e_(0,0). Then

    Gamma_ij=2 z^T(P_i tensor I) R(I tensor P_j)z.        (19)

To derive it, expand each convolution in (7), partition the four time variables by the order of the two conversions, and integrate the three positive time intervals. Each integral is a copy of R. Transposition identifies the two orders, giving the factor two. This is also the two-conversion pair-resolvent expansion. All matrices in (19) are rational at the displayed parameters.

[proved, certificate algorithm] Substitute x=tan(k/2) in (12). Its zero condition becomes

    det[2A(x^2-1)+4 eta B x+mu E_u(1+x^2)]=0.            (20)

The degree is at most 2n. Rational interpolation at 2n+1 integers computes this polynomial exactly. Sturm sequences isolate all 2n roots in (-2,2), each to interval width less than 2^(-220). There can be no additional root outside that interval because the degree and root count agree. The leading coefficient det(2A+mu E_u) is positive, so k=pi is not a missing root. The Sturm isolation and positive discriminant implied by the root count certify simplicity.

[proved, angular certification] The intervals are mapped to k=2 arctan x using arctan x=sign(x) pi/4+arctan((x-sign(x))/(1+sign(x)x)). An alternating series with its next-term bound encloses the second term. Machin's formula pi=16 arctan(1/5)-4 arctan(1/239), with the same bounds, encloses pi. Every bound is an exact Fraction. Outward rounding to a denominator 2^250 reduces arithmetic size without weakening validity. Comparing all unordered pair sums, including repeated indices, to all relevant multiples of 2pi yields certified lower bounds. Floating-point numbers in the summaries only display those rational bounds.

[proved, certified scope] The certificates concern the limiting Gaussian reference (12), with w=1/4, v=1/1024, u=1/64, eta=1/100, and mu=10^(-(8n+10)), for n=1,2,3. They do not claim these same fixed parameters work for every n, and do not certify the finite-D diagnostic used for ED. The all-n result uses the terminating search in Sections 4 through 6.

[numerical, independent local rules] literal.py constructs the repeated 2116 by 2116 physical core bond from the 46 letters and its positive boundary operators. Restriction to the legal spaces of dimensions 54,160,424 agrees with model.py after the explicit tree gauge, with maximum row-sum differences 0,5.9e-16,2.1e-16. The physical onsite-charge commutator vanishes entrywise. The excited kinetic blocks pass the diagonal-dominance check. An initial test incorrectly applied that check also to the weighted vacuum block; it failed because its diagonal q^2 can be smaller than q. That block is positive by its rank-one form (1), and the corrected test checks diagonal dominance only on the excited blocks.

[numerical, complete bond source ED] checks.py diagonalizes every excited-cell state, dimensions 24,72,192 per bus at n=1,2,3, and evaluates (6) on all 4,16,64 two-cell code states. No excitation or charge cutoff enters that bond. The saved three-cell coefficients are its two embeddings on code dimensions 8,64,512, as proved by the outer-cell rule following (6). These are full finite-g bonds, not the g^2 coefficients or their free fits. The direct g^2 pair inverse at D=10,100,1000 approaches (7) with relative matrix errors decreasing from .877 to .0223, .947 to .0302 and .979 to .0380 respectively. This checks the high-energy power and normalization, not an all-n extrapolation.

[numerical, finite-coupling physical ED] ed.py diagonalizes the actual legal-cell cluster Hamiltonian at n=1,2,3 and two/three cells in total charges Q=0,1 and even parity for each excited bus. These are exact invariant sectors containing every code state of those charges. No vacuum-walker or carrier excitation is truncated inside them. A polar code identification gives the effective finite-coupling matrix, compared against lambda^2 times (6) and the compressed readout. The other total charges are not included in this cluster ED. They are included in the complete source bond just described. The wrapper, wrong-marker sectors and mobile charge/helix excitations are not diagonalized.

[numerical, eigenvalue scope] The n=3 three-cell block solver returns the code-sector low vectors with residual below 2.8e-12. It does not return an additional excited level, so no next-level gap is reported for that run. ARPACK supplies that extra level for the smaller runs. The fourth-order error law is evidence for the finite-cluster reduction, not a rigorous eigenvalue enclosure or a central-charge measurement.

[proved, infrared scope] Finite-dimensional coefficient bounds do not imply a convergent extensive Schrieffer-Wolff band. C40's completed result separately requires an irreducible low-frequency massive-sector remainder assumption and matching constants for the slowest velocity. The present construction supplies the missing bare momentum witness; it does not prove that additional assumption. Reduce g to control the lambda^2 g^4 interaction relative to the lambda^2 g^2 hopping before choosing lambda. Reducing lambda alone does not reduce this ratio.
