# C96. Tagged contact, a weaker comparison norm, and thermal resummation

[Status] This is a partial resolution. The finite-volume tagged contact kernel is computed by an exact determinant spectral sum. Its uniform first absolute moment is not established. A weaker sufficient norm for C94 is proved, including a resummed quadratic version for product thermal cells. The full relative quantum thermal remainder and the C93 dressed contact-clock theorem remain open. No full d=13 mobile defect theorem is claimed.

[Convention] Proved tags derivations below. Imported retains a predecessor result with its scope. Numerical tags finite floating-point checks. Conditional tags implications with unproved hypotheses. Open identifies an estimate not obtained. A tag governs the following argument until the next tag. The gap ring has N material labels and v vacancies, L=N+v, t>0, and N>=3. All distances without a qualifier are material-label distances. No files outside work_C96 are changed.

## 1. Which kernel is required

[Imported, C94 section 7 and C95 sections 1 and 5] The required zero-temperature canonical kernel is

    chi_N(r) = <0|(A_0-a) Q(H_gap-E_0)^(-1)Q(A_r-a)|0>,
    A_j = 1_(g_j=0),
    H_gap = -t sum_j(E_j^+ E_(j+1)^- + h.c.).             (1)

C94 uses chi_N(r)=k_N(r)-c_N/N and a first absolute moment of k_N. C95 proves a bound on sum_r chi_N(r), and computes the thermodynamic coefficient c_dens. Neither determines the nonzero material-momentum eigenvalues of (1).

A physical contact D_x=n_x n_(x+1) and a tagged contact A_j are different operators. Only their sums are identified by the cyclic-invariant isometry. The exact answer below keeps every material-momentum sector. In particular, the two-particle-hole truncation used for the summed contact does not compute (1).

[Convention] When making a definite subtraction below, set

    c_N = c_dens(rho),     rho=N/L,
    k_N(r) = chi_N(r)+c_dens(rho)/N.                     (2)

This is a definition, not a proof of a uniform norm. The canonical response is the density-relaxed candidate minus c_dens/N. Adding a spatial constant changes only its zero Fourier mode. It cannot cancel a nonconstant algebraic tail.

[Imported, C95] With q=v/N,

    e''(q)=2 pi t rho^3 sin(pi rho),
    a_q=-rho^2+(rho/pi)sin(2 pi rho)-sin^2(pi rho)/pi^2,
    c_dens=a_q^2/(2e'').                                (3)

These constants and their derivations are not recomputed here. At rho=1/2, t c_dens=0.0785757974922801. This is the one-sided resolvent convention.

## 2. A complete finite-ring free spectrum with the tags retained

[Proved] Let xi=0 for odd N and xi=1/2 for even N. For each m=0,...,N-1 put theta_m=2 pi m/N. Choose an N-element subset I of {0,...,L-1} satisfying

    sum_(a in I) a + N xi + m = 0 modulo L.              (4)

The momenta, energy, and wavefunction are

    p_a = (2 pi/L)(a+xi+m/N),
    E_(m,I) = -2t sum_(a in I) cos p_a,
    x_1=0,  x_j=j-1+sum_(i=0)^(j-2) g_i,
    psi_(m,I)(g) = det_(j,a)[exp(i p_a x_j)]
                              /sqrt(N L^(N-1)).         (5)

Here g_j>=0 and sum g_j=v. The ordered x_j occupy a ring with a particle fixed at the origin. The common boundary phase and total physical momentum obey

    exp(i L p_a)=(-1)^(N-1) exp(i theta_m),
    exp(i sum_a p_a)=1.                                 (6)

To verify the eigenvalue equation, first apply the free hopping recurrence to the determinant in unrestricted ordered coordinates. Each column is a one-particle hopping eigenfunction. A forbidden collision has two identical rows and gives zero. A collision across the physical seam also gives zero because the rows at 0 and L are proportional by (6). Moving the particle at the root requires translating all relative coordinates back to x_1=0. Its determinant multiplier is exp(+/- i sum p_a)=1. Every remaining allowed move is precisely a vacancy transfer between two adjacent gaps. Thus the recurrence is H_gap psi=E psi, including moves at the gap seam.

Under the left cyclic rotation Rg=(g_1,...,g_(N-1),g_0), the relative coordinates become

    0, x_3-x_2, ..., x_N-x_2, L-x_2.

Translating all columns gives exp(-i x_2 sum p_a)=1. Moving the old first row to the last gives (-1)^(N-1), canceled by the first factor in (6). Therefore

    psi_(m,I)(Rg)=exp(i theta_m) psi_(m,I)(g).            (7)

Different m are orthogonal. At fixed m the physical one-particle plane waves on the common twisted grid are orthogonal. To transfer their Slater normalization, sum over occupation sets and over their N occupied roots. The root phase cancels in a product of two wavefunctions with the same m. The N choices of occupied root and the L translations of a rooted gap tuple give

    sum_g det_p(g)^* det_k(g)
        = (N/L) sum_(X, |X|=N) det_p(X)^* det_k(X)
        = N L^(N-1) delta_(I,J).                        (8)

This proves the normalization in (5).

For completeness, count the states. The gap space has dimension binom(L-1,N-1). Partition the N-element subsets I into orbits under I -> I+1 modulo L. If such an orbit has length d, then dN is divisible by L. Write h=gcd(N,L). Along the orbit, sum I runs through a residue class modulo h, with each member repeated d/(L/h) times. Among the N consecutive values m=0,...,N-1, exactly N/h meet the required residue class in (4). Thus the orbit contributes

    [d/(L/h)](N/h) = dN/L

admissible pairs (m,I). Summing over the orbits gives

    number of states = (N/L) binom(L,N)
                     = binom(L-1,N-1).                 (9)

The normalized eigenvectors (5) therefore form a complete basis. No assumed Bethe-state completeness is needed. The ground state is the centered sea k_a=2 pi[a-(N-1)/2]/L, a=0,...,N-1, in sector m=0. Its positive form and uniqueness agree with the inherited gap ground state.

## 3. Exact tagged contact form factors and kernel

[Proved] For two momentum sets p,k from (5), define N by N matrices

    S_(ab) = sum_(x=2)^(L-1) exp[i(k_b-p_a)x],
    U_(ab) = 1,
    V_(ab) = exp[i(k_b-p_a)].                            (10)

The exact tagged contact matrix element is

    F(p,k) = <p|A_0|k>
      = {det(S+U+V)-det(S+U)-det(S+V)+det S}
                                             /(N L^(N-1)). (11)

Every entry in S is a finite geometric sum. Equations (4),(5),(10),(11) give a finite expression in roots of unity, without a many-body matrix inversion.

To derive (11), A_0=1 fixes the first two particles at 0 and 1. The other N-2 occupied positions are an arbitrary subset of {2,...,L-1}. Cauchy-Binet applied to

    det[S+z U+w V]

expands it into products of the two Slater determinants, weighted by z when 0 is occupied and w when 1 is occupied. Since U and V each have rank one, the determinant has degree at most one in each variable. Its zw coefficient is the numerator of (11). That coefficient sums precisely the gap configurations with g_0=0. Dividing by the two normalization factors from (5) proves the formula. This derivation also covers singular S or singular full overlap matrices.

Let k denote the ground sea and let F_(m,I)=F(p_(m,I),k). Changing variables by R^r in the gap sum and using (7) gives

    <m,I|A_r|0> = exp(i theta_m r) F_(m,I).              (12)

Consequently the exact answer for (1), and for the specified k_N, is

    chi_N(r) = sum_((m,I) != ground)
           exp(i theta_m r) |F_(m,I)|^2/(E_(m,I)-E_0),
    k_N(r) = c_dens(rho)/N + chi_N(r).                  (13)

Reflection symmetry makes these quantities real and even. Equivalently, if the Fourier convention is chi_hat(theta)=sum_r exp(-i theta r)chi_N(r),

    chi_hat(theta_m) = N sum_(I in sector m, excited)
                               |F_(m,I)|^2/(E_(m,I)-E_0). (14)

Each eigenvalue is nonnegative. The subtraction in (2) adds c_dens to the m=0 eigenvalue and leaves every other eigenvalue unchanged.

There are two immediate normalization checks. Completeness and A_0^2=A_0 give sum_(m,I)|F_(m,I)|^2=a. At m=0, sum_j A_j has matrix element N F_(0,I). It becomes the physical total contact, so only two-particle-hole states contribute away from the ground sea, and

    N F_(0,I) = (2/L)[cos(p_1-k_1)-cos(p_2-k_1)]         (15)

up to the Slater basis sign, with p_1+p_2=k_1+k_2 modulo 2 pi. Thus (14) at m=0 reproduces C95, including its normalization.

For m!=0, the initial and final momentum grids differ by theta_m/L. The orthogonality cancellation that gave (15) no longer applies. Formula (11) is a determinant of overlaps of all occupied orbitals. It is not a sum restricted to two replacements of the original sea.

[Numerical] The finite check N=v=6 includes m=5 and I=(1,3,7,8,10,11), with |F|=0.0301981882773300. The ground set on the untwisted index grid is (0,1,2,9,10,11). Three indices are replaced, in addition to the change of twist. This illustrates the distinction in (11), rather than defining a particle-hole convention between different grids.

[Proved, physical-coordinate comparison] One can compute the different kernel of D_x by ordinary one- and two-particle-hole form factors. In a centered physical sea put g_1=L^(-1)sum_(h occupied)cos h and rho=N/L. Normal ordering D_0 gives the one-particle-hole amplitude

    V_(p,h) = [rho(1+exp(i(h-p)))
                         -g_1(exp(i h)+exp(-i p))]/L.   (16)

The two-particle-hole amplitude, for distinct ordered pairs, is

    W_(p1,p2;h1,h2)
       = (exp(-i p2)-exp(-i p1))
                         (exp(i h2)-exp(i h1))/L^2,     (17)

up to the basis sign. These follow either by normal ordering

    D_0 = c_0^dagger c_1^dagger c_1 c_0

or by evaluating the two-site Slater minors. The physical local resolvent is the sum of |V|^2/Delta_ph and |W|^2/Delta_2ph, each multiplied by the corresponding physical momentum phase at separation x. A spatial sum kills (16) and imposes zero momentum transfer in (17). A tagged Fourier mode instead changes the twist as in (5). Substituting the physical local kernel for (13) would compute the wrong observable.

[Open] Formula (13) computes every finite-N kernel, but an all-N estimate of its spatial tail has not been proved. In particular this work does not decide whether sup_N t(K_N+|c_N|) is finite in the requested density window. It establishes neither an r^(-2) tail nor its absence. Fixed-volume agreement of (13) with a resolvent is not a uniform asymptotic estimate. The sums over twisted seas must be estimated before taking a first absolute moment.

## 4. A sufficient norm that does not require the first moment of k

[Proved] There is a weaker replacement for C94's product-cell theorem. For a real even circulant kernel define

    S(k)=sum_r |k(r)|,
    Delta k(r)=2k(r)-k(r+1)-k(r-1),
    D(k)=sum_r [1+dist_N(r,0)] |Delta k(r)|,
    M_b(k,c)=S(k)+b^2 D(k)+|c|.                          (18)

Let chi=k-c/N be positive semidefinite. In C94's zero-port product-cell setup, with ||h_i||<=M and the same W, one has

    |B_sigma-B_good| <= 32 M^2 b W M_b(k,c).             (19)

It is enough to take k=chi,c=0 in this statement. Thus the weaker theorem does not require a separate identification of the canonical constant term. When a decomposition is known, Delta chi=Delta k and S(chi)<=S(k)+|c|.

Here is the proof, including the cancellation that C94's entrywise estimate did not use. Write

    m_i=<h_i>,
    G_ij=m_i m_j+C_ij,
    B=sum_ij chi(j-i)G_ij.                              (20)

Product cells make C block diagonal in the cells. Let L_bad be the material length of wrong cells and pinned selected cells, so L_bad<=bW. Each complete good cell has the same mean profile of period b. Let e be the average of that profile, |e|<=M.

On every good run choose a discrete primitive u, reset to zero at each cell boundary, such that

    m_i=e+(u_i-u_(i-1)),    |u_i|<=2bM.                 (21)

The primitive exists because the mean of m_i-e over a whole selected cell is zero. Set u=0 throughout the bad set. Then globally

    m=e+d+D_1 u,
    supp d subset bad,   |d_i|<=2M,   |D_1 u_i|<=2M,
    D_1 u_i=u_i-u_(i-1).                                (22)

The reference has d=0 and a b-periodic primitive u_good bounded by 2bM. On every good run u is a translate of u_good.

In the quadratic form of k, the common constant term cancels. Since sum_i D_1u_i=0, its cross term with e vanishes. The other terms involving d obey

    |2e (sum k)(sum d)| <=4M^2 L_bad S(k),
    |d^T k d| <=4M^2 L_bad S(k),
    |2d^T k D_1u| <=8M^2 L_bad S(k).                   (23)

For the remaining term, discrete summation by parts on the ring gives

    (D_1u)^T k(D_1u)=u^T(Delta k)u.                    (24)

This is an exact identity, without boundary terms. Apply the translated-good-run count from C94 to u instead of m. If there are k_run good runs, at most 2k_run dist_N(r,0) origins have their shortest arc leave a good run. Bad origins and the omitted reference fraction cost at most 2(2bM)^2 L_bad. Thus

    |sum_i u_i u_(i+r)-sum_i u_good,i u_good,(i+r)|
       <=2(2bM)^2 L_bad
                         +4(2bM)^2 k_run dist_N(r,0).  (25)

Since k_run<=W when W>0, summing against |Delta k(r)| bounds (24)'s difference by 16 b^3 M^2 W D(k). Equations (23)-(25) bound the regular mean part by

    16M^2 bW [S(k)+b^2D(k)].                            (26)

For the constant part, sum_i m_sigma,i differs from Ne by at most 2M L_bad, while each total mean has magnitude at most MN. Therefore

    (|c|/N)| (sum m_sigma)^2-(Ne)^2 |
                                  <=4 |c| M^2 L_bad.   (27)

Finally, C on each cell is a positive semidefinite covariance matrix with trace at most M^2 times its length. Put Lambda=||chi||_op<=S(k)+|c|. Each cell contributes between zero and Lambda M^2 ell to Tr(chi C). Every identical selected cell has the same contribution, regardless of its absolute position. Their contributions cancel in the comparison. The wrong cells and the missing reference cells together cost at most

    2 Lambda M^2 L_bad.                                (28)

Combining (26)-(28) gives coefficients 18 for S, 16 for b^2D, and 6 for |c|, all below the conservative common constant 32 in (19). This proves the claim. It uses no gap in a wrong cell and retains the connected quantum source covariance.

[Proved, scope of the replacement] A smooth nonoscillatory r^(-2) tail need not spoil (19). For example take kappa(r)=A/(1+r^2), A>0, on the infinite line and periodize it,

    k_N(r)=sum_(z in Z) kappa(r+zN).                     (29)

Its first absolute moment grows at least logarithmically in N, since the terms 1<=r<=N/2 already contain A r/(1+r^2). But sum|kappa| is finite. For r>=2, the second-derivative formula gives

    |Delta kappa(r)| <=8A/(r-1)^4.

Indeed Delta is minus the average of kappa''(r+u) with triangular weight 1-|u| on [-1,1], and |kappa''(x)|<=8A/x^4 for x>=1. Hence sum_r(1+|r|)|Delta kappa(r)| is finite. Periodization commutes with Delta, and dist_N(r,0)<=|r+zN|, proving uniform S(k_N) and D(k_N). Subtracting c/N adds at most |c| to S and nothing to D.

For the test kernel (29), an explicit bound is available. Its second difference is positive only at r=0 and is nonpositive elsewhere. Summing Delta kappa gives zero, while Delta kappa(0)=A, so sum|Delta kappa|=2A. The first-moment telescoping identity is

    sum_(r=1)^R r[-Delta kappa(r)]
        =kappa(0)-(R+1)kappa(R)+R kappa(R+1) -> A.

Thus sum|r||Delta kappa|=2A and D(kappa)=4A. The Fourier transform of A/(1+x^2) is A pi exp(-|q|), obtained by closing the contour on its pole at +i or -i. Summing its Fourier series at integer positions gives S(kappa)=A pi[1+2sum_(j>=1)exp(-2pi j)]=A pi coth(pi). Periodization therefore gives M_3(k_N,0)<=A[36+pi coth(pi)]<40A.

This conclusion requires cancellation or regularity of the tail, not just an upper bound C/(1+r^2). For example kappa(r)=(-1)^r A/(1+r^2) resonates with a period-two mean profile. After multiplying m_i by (-1)^i, a phase-shifted selected run becomes a sign domain for the positive kernel A/(1+r^2). Its cross-boundary quadratic cost contains

    sum_(i in domain, j outside) 1/[1+(i-j)^2].

A domain of length R has a number of crossing pairs of separation d proportional to d for 1<=d<=R. The sum grows as log R. Two wrong odd-length cells can create and restore that phase shift while their total wrong length stays fixed. Bounded changes inside those cells contribute only O(S(kappa)). This is an admissible positive semidefinite test kernel: the unmodulated sequence has a positive Fourier series by the preceding periodization formula, modulation by (-1)^r shifts that positive symbol, and periodization on an even ring preserves positivity. Thus an arbitrary oscillatory r^(-2) tail cannot be handled by an unweighted absolute sum alone. The norm (18) correctly remains divergent for this example.

[Open] Neither S(chi_N) nor D(chi_N) has been bounded uniformly for the exact tagged kernel (13). Equation (19) is a proved replacement criterion. It is not a claim that its hypotheses have now been verified for C88.

## 5. Thermal quadratic control with the actual frozen dynamics

[Proved] Fix J,T and the free bath ensemble, canonical or tilted, with contact mean a_0. Introduce an independent interpolation parameter z and set

    Href_sigma = H_bath+J a_0 Hmat_sigma,
    Y_sigma = sum_i(A_i-a_0)h_sigma,i,
    H_sigma(z)=Href_sigma+z JY_sigma,   0<=z<=1.         (30)

The frozen material Hamiltonian is held at its actual J, not Taylor expanded at J=0. Define the free-energy residual at z=1 as in C95. The product reference gives <Y>_(z=0)=0. Differentiating in z twice and integrating yields the exact finite-volume identity

    R_sigma=-J^2 integral_0^1(1-z) Xi_sigma(z) dz,
    Xi_sigma(z)=integral_0^beta
                 <Y_sigma(tau);Y_sigma>_(H_sigma(z))d tau. (31)

Here the semicolon means the connected thermal correlator. Its spectral weights are nonnegative, so Xi_sigma(z)>=0. The product reference trial also gives R_sigma<=0 directly.

At z=0 put

    C_b,ij(tau)=<A_i(tau)A_j>_b-a_0^2,
    G_sigma,ij(tau)=<h_i(tau)h_j>_(J a_0 Hmat_sigma),
    B_sigma,res=(1/2)sum_ij integral_0^beta
                               C_b,ij(tau)G_sigma,ij(tau)d tau.
                                                               (32)

Then Xi_sigma(0)=2B_sigma,res. No replacement of h_i(tau) by its expectation or its normalized trace has occurred.

[Proved, a uniform quadratic inequality] Define the one-sided static thermal bath matrix

    chi_beta=(1/2) integral_0^beta C_b(tau)d tau.          (33)

The same determinant formula also gives this matrix exactly at fixed N,v,T. If a,b label all states in (5), set w_a=exp(-beta E_a)/Z and

    L_ab=(w_a-w_b)/(E_b-E_a) if E_a!=E_b,
    L_ab=beta w_a if E_a=E_b.

With F_ab computed by (11), its row is

    chi_beta(r)=(1/2)sum_ab L_ab |F_ab|^2
                 exp[i(theta_b-theta_a)r]-beta a_0^2/2. (33a)

This follows by integrating each spectral exponential. The equal-energy prescription includes all degeneracies. As T tends to zero at fixed canonical volume, the two ground-to-excited orientations cancel the factor 1/2 and reproduce (13). A length-tilted trace also sums the fixed-v blocks and adds the between-block contact covariance from C95. None of these finite-volume identities justifies interchanging that temperature limit with the thermodynamic limit.

For the real gap Hamiltonian and the real diagonal A_i, C_b(tau) is a real symmetric positive semidefinite matrix at every tau. In a real bath eigenbasis its off-diagonal-state terms are positive multiples of v v^T, with v_i=<a|A_i|b>. Its diagonal-state connected part is a covariance matrix and is also positive semidefinite.

For arbitrary Hermitian material sources, let m_i=<h_i> and let C_m(tau) be their connected correlator. As a real quadratic form,

    0 <= Re C_m(tau) <= C_m(0)_real,   0<=tau<=beta.     (34)

To see this, take any real vector u and O=sum u_i h_i. In a material eigenbasis pair energies E_a<=E_b, write delta=E_b-E_a>=0 and w_a=exp(-beta E_a)/Z. The pair's contribution is

    w_a[exp(-tau delta)+exp(-(beta-tau)delta)]|O_ab|^2.

The bracket is positive and convex in tau, with equal endpoint values 1+exp(-beta delta). It is bounded above by that endpoint value. The diagonal part after subtracting <O>^2 is the nonnegative variance of the diagonal entries and is independent of tau. This proves (34) for every real u, including degeneracies.

Taking the trace against the positive bath matrix in (32) therefore gives

    0 <= B_sigma,res-m^T chi_beta m
       <= Tr[chi_beta C_m(0)_real]
       <= ||chi_beta||_op sum_i Var(h_i)
       <= ||chi_beta||_op M^2 N.                        (35)

There is no inverse material gap, no bound on imaginary-time operator norms, and no factor exp(J/T) in this inequality. A full spectral proof works also for complex material Hamiltonians, since the real bath matrix contracts only the real symmetric part of C_m.

[Proved, relative product-thermal-cell version] Suppose the zero-port material Hamiltonian is a sum of isolated cell Hamiltonians and each h_i acts within its assigned cell. This includes isolated legal cells with their auxiliary states retained whenever they have that factorization, not just a vacuum-clock vector. The frozen thermal state is then a product over cells, and C_m(tau) is block diagonal. Each identical good cell makes the same contribution to the left side of (35), independent of its position in the bath. The cancellation used in (28) therefore works for B_res too. The mean profile is b-periodic on every good run and obeys |m_i|<=M at every J/T. Applying (21)-(27) to those thermal means proves

    |B_sigma,res-B_good,res|
         <=32 M^2 bW [S(k_beta)+b^2D(k_beta)+|c_beta|],
    chi_beta=k_beta-c_beta/N.                           (36)

This is uniform in J/T, the internal dimensions, and the wrong-cell energy gaps, at any fixed bath kernel. It extends the conditional coefficient comparison to product thermal cells without the frozen-subtraction Taylor error. A uniform bound on the displayed thermal bath norm, with N tending to infinity first, remains an additional hypothesis.

The argument does not extend the good-cell cancellation through the nonzero port. There the frozen material correlator need not be block diagonal. Inequality (35) still bounds its total-volume coefficient, but a total-volume bound cannot replace the W-relative estimate.

## 6. Even the resummed quadratic coefficient is not a uniform full remainder

[Proved, exact diagnostic] Keeping J a_0 Hmat inside the reference fixes C95's specific first-order material-weight expansion. It does not by itself give a temperature-uniform remainder in (31). Consider two independent two-level bath spins and one material projector h,

    H_b=t(tau_z,1+tau_z,2),
    A_1=(1+tau_x,1)/2,  A_2=(1+tau_x,2)/2,
    h_1=h, h_2=-h, h=diag(0,1).                         (37)

Every contact is a projector, each source has norm one, a_0=1/2, and Hmat=h_1+h_2=0. Thus the frozen material reference is already retained exactly at every J/T. The interpolation (30) becomes

    H(z)=H_b+(zJ/2)(tau_x,1-tau_x,2)h.

With d_z=sqrt(t^2+z^2J^2/4), its exact residual is

    R(z)=-T log{[1+(cosh(d_z/T)/cosh(t/T))^2]/2}.        (38)

The h=0 block has the original bath trace and the h=1 block has two spins of energy scale d_z. This proves (38) without a perturbative expansion.

For fixed T>0, expansion in z gives

    R(z)=-z^2 J^2 tanh(t/T)/(8t)
       +z^4 J^4[tanh(t/T)-t/T]/(128t^3)+O_T(z^6J^6).   (39)

For the algebra, use d_z-t=z^2J^2/(8t)-z^4J^4/(128t^3)+..., expand twice the log of the cosh ratio, and then use

    log[(1+exp x)/2]=x/2+x^2/8+O(x^4).

The two fourth-order terms containing sech^2(t/T) and tanh^2(t/T) sum to one. The fourth-order coefficient in (39) diverges as -J^4/(128t^2T).

At fixed J>0 and T->0,

    R(1)->-2[sqrt(t^2+J^2/4)-t]
          =-J^2/(4t)+O(J^4/t^3),
    -J^2 B_res -> -J^2/(8t).                            (40)

Thus the error in the resummed quadratic free-energy coefficient is of order J^2/t, not o(J^2/t), uniformly at low T. The bath is gapped and there is no vacancy crossing. Replicate the independent units and divide the free energy by their number. The same identities hold after taking the thermodynamic limit at every fixed T>0. The thermodynamic order of limits does not remove this diagnostic.

The exact integrated bound is nevertheless

    -J^2/(4t) <= R(1) <= 0.                             (41)

To prove it, put

    delta_eff(J,T)=2T log[cosh(d_1/T)/cosh(t/T)].

Since the derivative of T log cosh(x/T) is tanh(x/T) between zero and one for x>=t>0,

    0<=delta_eff<=2(d_1-t)<=J^2/(4t).

Equation (38) is the free energy of the material operator -delta_eff h minus its zero-Hamiltonian free energy. Monotonicity and the operator-norm Lipschitz bound for free energies give (41). Retaining this induced operator inside the material exponential performs the required resummation in this diagnostic.

[Scope] This is not a C88 counterexample. It disproves a proposed general inference from bounded contacts, a free-bath quadratic bound, and retention of the frozen material dynamics to a uniform smaller remainder around that quadratic free-energy coefficient. Bath-induced material splittings on the scale J^2/t must themselves be allowed to reorganize the low-temperature trace. A direct relative bound on the full residual can still hold, just as (41) does.

[Open, precise full-chain target] In the notation of (31), the needed one-sided estimate is

    integral_0^1(1-z)[Xi_sigma(z)-Xi_good(z)] dz
           <= [epsilon(J/t)/J] W+B_N/J^2,
    epsilon(x)->0,
    B_N/N->0 at every fixed T>0.                        (42)

A direct bound on a relatively matched bath-induced operator or influence functional, retained inside the material trace, could prove (42). Neither (13), (19), nor (36) controls Xi(z) at z>0. No bound on (42) uniform in arbitrary wrong lengths, the full auxiliary trace, and the nonzero port has been obtained here. It would be incorrect to report (36) as control of this full remainder.

## 7. What follows for d=13 and for the contact-clock problem

[Imported] C92's positive contact expectation, legal selected-pin thermal replacement, and whole-sector exclusions remain valid with their original scope. C95's signed contact response and density-relaxation calculation also remain available. No predecessor certificate is weakened by the new determinant representation.

[Proved conditional] For b=15, a uniform bound t M_15(k,c)<=C in (18) would replace the first-moment hypothesis of C94 by

    |Delta B|<=480 C M^2 W/t.                            (43)

The same numerical prefactor applies to (36) under its product-thermal-cell and bath-norm hypotheses. It is deliberately conservative. A smooth nonoscillatory r^(-2) tail would therefore not by itself prevent this coefficient comparison.

[Conditional] If the full legal estimate (42) is proved with error below C92's selector margin, the already established length-weighted activity sum gives a pressure difference O(T exp(-delta/T))=o(T^2), with N or L tending to infinity first. Finite grammar penalty separately requires C92's relative penalty correction. Nothing in the present static calculation proves that correction.

[Open, C93] The tagged static resolvent is not a dressed clock gap theorem. It integrates over bath excitation energies and supplies no lower bound on a massive material complement, no product auxiliary heat trace, and no mixed-bath multi-time influence estimate. C93's order-J contact-clock leakage remains in the microscopic Hamiltonian. Neither boundedness of a static kernel nor the product-cell comparison proves a dressed threshold or certifies thirteen mobile payload channels.

## 8. Checks and final ledger

[Numerical] tagged_kernel.py evaluates (11)-(14) and independently solves the gap resolvent. At N<=4 it also constructs the complete determinant eigenbasis and verifies normalization, orthogonality, the eigenvalue equation and each A_0 form factor. Half-density checks through N=8 have maximum bath dimension 6435. The largest row discrepancy is below 3e-15 and the largest independent resolvent residual is below 4e-13. The full basis checks have errors below 7e-15. These are small finite controls, not a thermodynamic computation.

With c_N fixed to (3), the half-density first absolute moments at N=3,...,8 are approximately

    0.1087344, 0.1099512, 0.1100844,
    0.1102580, 0.1102139, 0.1102238                       (t=1).

They do not resolve a tail or prove a uniform bound. No logarithmic asymptotic is inferred from them. CHECKS.md records the additional thermal and norm controls.

[Proved, consistency] The dimensions in (9) agree with the gap compositions. Equation (15) reproduces the inherited zero-momentum form factor. The equal-time spectral sum equals a because A_0 is a projector. At v=0 the gap space has one state and the connected canonical kernel vanishes. All response kernels and their norms have units 1/t. Delta annihilates a constant kernel. Equation (34) includes arbitrary degeneracies, so no hidden wrong-cell gap enters (36). Equations (38)-(41) exhibit a bounded exact residual despite divergent Taylor coefficients, with the thermodynamic limit taken first by replication.

[Final ledger] The exact finite-N tagged spectral formula, the weaker comparison criterion, the product-thermal quadratic bound, and the resummed-remainder diagnostic are proved. Uniform tagged spatial norms, their low-temperature thermodynamic counterparts, the nonzero-port relative remainder, and the dressed contact-clock theorem remain unproved. The full d=13 mobile comparison is therefore still conditional.
