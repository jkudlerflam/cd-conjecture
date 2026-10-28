# C44 mobile-wrapper audit

[proved, scope] For every positive integer n, C39's initialized core admits the C8 mobile wrapper on 4464 letters. This note proves its finite-dimensional hypotheses, exclusion bounds, exact projected maps and second-order source formula, including bounds on the length-dependent scalar and payload terms. It supplies finite parameter choices at each n. No bound uniform in n is claimed.

[assumed, infrared scope] A thermodynamic local elimination theorem with errors small relative to the slowest retained velocity, followed by stability of the resulting multicomponent Gaussian theory, is still required for c=n+1 at fixed nonzero microscopic couplings. A finite-cluster band expansion does not supply that theorem. Statements about suppression of symmetry-allowed continuum operators below refer to their finite bare coefficients and momentum margins, not a proof of their infrared flow.

## 1. Notation and the actual repeated bond

[proved, definitions] Use k=n+2 occupied sites, b=k+1, physical cell lengths N_p=k+p with p=1,2, and logical dimension d_log=2^n. To avoid collisions of notation write rho=w_vac/U for the wrapper vacancy ratio, t_sel for the selector head hopping, u_mov for marker hopping, s for the overall core scale, and J_h for the endpoint-helix scale. C39's parameters are w,v,u,eta,mu,D,g,q, with q=(wu)^2. None depends on the eventual physical ring size L.

[proved, construction] Let h_c,a_c,b_c be the positive repeated bond and positive endpoint terms in C39 PROOFS Section 1. On occupied sites tensor the qubit occupation filter onto the 46-letter core. Its pitch is phi_f=2pi/(k+2), its left pin is |+0>, and its right pin is |+(k-1)phi_f>. The occupied bond and endpoints are

    h_occ=s h_c tensor I + J_f I tensor h_phi_f,
    a_occ=s a_c tensor I + J_f I tensor a_f,
    b_occ=s b_c tensor I + J_f I tensor b_f.

The wrapper letters are occupied tags L,B,R, each carrying this 92-dimensional space, and holes H_0,H_1,H_2. Their arrows are respectively 0->1,1->1,1->2 and j->j. On matching occupied pairs lift h_occ; add oriented endpoint terms (a_occ-2U) on a left-tag letter in the left slot and (b_occ-2U) on a right-tag letter in the right slot. Add 6U on each mismatched wrapper bond. The charge bond has diagonal U per occupied-hole interface and hopping -w_vac on

    H_0 L_alpha <-> L_alpha H_1,
    H_1 B_alpha <-> B_alpha H_1,
    H_1 R_alpha <-> R_alpha H_2.

[proved, construction] Tensor the four-letter selector #,L,X,R onto this wrapper. Its repeated bond is C8 selector_wrapper (4): illegal-selector penalty K_sel, half marker chemical potentials on each slot, head hopping -t_sel(LX<->XR), and the wrapper bond on legal pairs whose right selector is not #. This last condition cuts every wrapper bond entering a marker and provides the arrow reset. Add the two endpoint-helix bonds on the same internal selector projector. On the cut projector X# + R# pin track B on the left and track A on the right. All pins and helix bonds are positive. This specifies one uniform two-site operator, without an imposed marker sector.

[proved, count] The onsite dimension is

    4 [3*2*46+3] 2*2 = 4464.

The entire construction is nearest neighbor. The n-dependent pitches, phases, masses and coefficients are entries of that same fixed-size bond. A constant alphabet does not mean a bond independent of n.

[proved, perturbations] Lift C39's two-bus endpoint operators to the appropriate occupied endpoint tag and multiply by the opposite endpoint-helix compensator specified in Section 4. Enable their product on X# + R#. Marker hopping is C8's R#<->#L move, with H_0->H_2 on the transferred site and the two local helix rotations. The conditional data readout is a repeated onsite operator on occupied letters, split in halves between physical bonds if desired. It vanishes on holes. A marker chemical-potential correction is also a repeated onsite term. Thus no flavor-addressed physical letter or nonlocal transport operator is needed.

## 2. Core hypotheses, including biased initialization

[proved] C39's excited kinetic part allocates Q_left+2Q_right on each internal bond. Its allowed excited off-diagonal moves have row sums bounded by that allocation, also on illegal words and with both buses present. Subtracting positive two-state edge blocks leaves a nonnegative diagonal. The remaining onsite block is (1+D P_high)Q+g C and has least excited eigenvalue at least m=1-|g|>0. Restore the missing onsite halves and 2Q at the left endpoint and Q at the right endpoint. These endpoint terms are positive. The vacuum edge is (q|initial>-|final>)(q<initial|-<final|), a positive rank-one matrix. Positive grammar penalties complete the definition.

[proved] The outer arrows of every kinetic rewrite are unchanged, and conversion preserves all arrows. Hence core violations are conserved. They cost at least K_c. Positivity first excludes excited tokens from a zero vector. The legal initialized word then has precisely one vacuum token, and its connected path imposes amplitude ratio q. The complete code and a lower gap bound are

    G|a> = S^(-1/2) sum_(ell=0)^(k-1) q^ell |ell,v;a>,
    S=sum_(ell=0)^(k-1) q^(2ell),
    gamma_c=min(K_c,1-|g|,1+q^2-2q cos(pi/k))>0.

The last eigenvalue formula follows by taking the row Gram of the incidence matrix with rows q e_j-e_(j+1). It has diagonal 1+q^2 and off-diagonal -q. There are precisely 2^n code vectors. This is an initialized, uniform local Hamiltonian. C8 does not require equal history amplitudes, only termwise annihilation by positive repeated terms.

[proved] Every wrapper track is neutral under C39's exact integer charge. Vacancy swaps move the complete occupied label. Marker motion transfers a hole, never a charged occupied label. Both bus port pairs conserve total charge. Therefore U(1) and its parity subgroup are exact on the unrestricted physical space.

[proved] Every endpoint F changes a vacuum token to an excited token or conversely, so G^dag F G=0. The two buses are separately conserved by the unperturbed Hamiltonian. These facts, rather than the Bell-syndrome rule of C25's different example, establish the selection rules here.

## 3. Occupation, arrows, arbitrary lengths and the selector

[proved] Untwisting the occupation-filter helix gives ferromagnetic singlet projectors. Their kernel is the symmetric subspace. The left pin forces a product state, and the right pin accepts it iff the occupied length r equals k modulo k+2. In 2<=r<=k+2 only r=k is accepted. Moving the first-site pin to site j with adjacent swaps gives

    <sum_j n_j> <= r(2r-1)<H_left>.

On the complement of the product vector sum_j n_j>=1, so the left-pinned gap is at least kappa_r=1/[r(2r-1)]. For the rejected right pin let beta_r=sin^2((r-k)phi_f/2). Writing psi=c g+z and using ||z||^2<=E/kappa_r and sqrt(beta_r)|c|<=sqrt(E)+||z|| gives

    F_r >= kappa_r beta_r / [(1+sqrt(kappa_r))^2+beta_r].

Consequently set

    kappa_*=1/[(k+2)(2k+3)], beta_*=sin^2(pi/(k+2)),
    epsilon_*=kappa_* beta_*/[(1+sqrt(kappa_*))^2+beta_*],
    J_f=2U/epsilon_*,  0<rho<1/(k+1).

[proved] At zero vacancy hopping, a wrong occupied number in a legal two-endpoint word either has at least two interfaces, costing 2U, or a contiguous occupied interval where the initialized filter costs at least 2U. Switching on hopping lowers any energy by at most w_vac(N_p-1). A target end-hole droplet has trial energy U. Thus every wrong occupation lies at least

    U-w_vac(N_p-1) >= U[1-rho(k+1)] > 0

above the selected code. This proof never uses the spectrum of C39 at other occupied lengths.

[proved] A legal wrapper word has at most one left and one right endpoint tag. With both present, their reward is -4U and the ordered-coordinate Hamiltonian is

    H_ch tensor I + sum_j D_j tensor (s h_c,j+J_f h_f,j)
       + I tensor [s(a_c+b_c)+J_f(a_f+b_f)] -4U.

Here D_j indicates adjacency of occupied labels j,j+1. Every vacancy configuration appears with a strictly positive amplitude in the unique charge ground vector chi_p. Subtract its ground energy e_p. Positivity forces any zero vector into chi_p tensor the intersection of all core and filter kernels. Thus the complete selected code is chi_p tensor G tensor g_f, with energy -4U+e_p, 0<=e_p<=U. Missing one endpoint costs at least U, missing both costs at least 3U.

[proved] Every wrapper violation projector commutes with the Hamiltonian because the three hole rewrites preserve the pair's outer arrows. With z bad bonds on an open chain there are at most 2(z+1) endpoint tags. Dropping positive terms gives

    H_wrap(l)>=6Uz-4U(z+1)>=-4U

for every l, including arbitrary illegal core sectors. On a periodic wrapper chain the endpoint count is 2z and H_wrap,ring>=2Uz>=0. These bounds are independent of D, J_f and the norms of all positive core terms. This independence is necessary for C39's lacunary schedule.

[proved] Define a_l=2(1-cos(pi/l)) and

    delta_b=(b+1)a_b-b a_(b+1), mu_c=2-delta_b,
    e_c=a_(b+1)-a_b<0,
    eta_b=min((a_(b-1)-2a_b+a_(b+1))/2,
              (a_b-2a_(b+1)+a_(b+2))/3)>0.

Convexity of a(x) for x>=2 gives eta_b>0, of order b^-4. Also mu_c>0 for b>=4. Choose

    t_sel=2U[1/eta_b+b/mu_c+1/|e_c|],
    mu_sel=t_sel mu_c+4U-(b+1)e_1+b e_2,
    e=t_sel e_c+e_2-e_1,
    K_sel=4t_sel+9U.

Then mu_sel>0, e<0 and E_b=eb,E_(b+1)=e(b+1). For l>=2 not selected, convex increasing secant slopes give

    a_l-a_b-(l-b)(a_(b+1)-a_b)>=eta_b(1+|l-b|).

Using H_wrap(l)>=-4U and 0<=e_p<=U now gives

    E_l-el >= (t_sel eta_b-U)(1+|l-b|)>0.             (1)

[proved] Legal marker-free sectors are periodic nonnegative wrapper chains and lose to eL. Cutting at z selector violations gives legal strings, each with at most two incomplete fragments. Each fragment costs at worst -2t_sel-4U, while complete cells obey H>=e times length. Therefore

    H_sector>=eL+(K_sel-4t_sel-8U)z=eL+Uz.           (2)

Positive endpoint helices cannot weaken either bound. For every L>=b(b+1), a b/(b+1) tiling exists by choosing M in [L/(b+1),L/b]. These tilings with the codes described above exhaust the unrestricted unperturbed ring ground space. No claim of uniqueness or of this same tiling description on all shorter rings is needed.

[proved, explicit gap] Let h_* =K_c+D+12 and a_* =K_c+D+5. Row sums of C39's local matrices give ||h_c||<=h_* and ||a_c||,||b_c||<=a_*. Put

    M_c=(k-1)h_*+2a_*,
    M_occ=s M_c+J_f(k+1),
    gamma_occ=min(s gamma_c,J_f/[k(2k-1)]).

Let delta_ch,p be the charge gap and d_p=min_j <D_j>. Decomposing a vector into the charge ground space and its complement, then applying the triangle inequality to the square root of the positive occupied interaction, gives

    Delta_occ,p >= delta_ch,p d_p gamma_occ /
       [(sqrt(delta_ch,p)+sqrt(M_occ))^2+d_p gamma_occ]. (3)

Indeed the complementary norm squared is <=E/delta_ch,p, the compressed occupied energy is >=d_p gamma_occ times the ground-component norm squared, and its square root is <=sqrt(E)(1+sqrt(M_occ/delta_ch,p)). Rearranging proves (3). Minimize (3) with U[1-rho(k+1)] to obtain Delta_p. The full base gap is at least the minimum of Delta_1,Delta_2,

    J_h/[(b+1)(2b+1)],
    2t_sel[cos(pi/b)-cos(2pi/b)],
    2t_sel[cos(pi/(b+1))-cos(2pi/(b+1))],
    t_sel eta_b-U, U, |e|b(b+1).                       (4)

It is volume independent above the entire unperturbed code manifold. It is not a gap inside the interacting length fluid.

## 4. Charge scale, endpoint compensation and exact transport

[proved] H_ch has dimension m_p=binomial(k+p,p), diagonal U times the interface count, and off-diagonal -w_vac on adjacent exchanges. Its graph is connected with diameter kp. Perron positivity proves uniqueness. Reflection gives equal endpoint weights r_p at the two ends even though C39's vacuum is biased.

[proved, fully explicit charge-gap bound] Put a_p=rho/(N_p-1). Since the charge bond is positive, e_p>=0. The ground eigenvector equation gives chi_i>=a_p chi_j on every neighboring configuration. Thus min chi>=a_p^(kp)/sqrt(m_p). The ground-state transform gives the quadratic form

    w_vac sum_edges chi_i chi_j |f_i-f_j|^2,
    ||chi f||^2=sum_i chi_i^2 |f_i|^2.

For mean-zero f, the variance identity and a path of at most kp edges between every pair bound the variance by kp sum_edges |f_i-f_j|^2. Therefore the conservative positive bound

    delta_ch,p >= w_vac a_p^(2kp)/(m_p kp)              (5)

is sufficient in every parameter inequality below. This avoids relying on an unquantified asymptotic coefficient to choose s.

[proved, sharper asymptotic] The two end-hole droplets are the only charge minima at w_vac=0. Moving p holes from one end to the other takes kp swaps. Every shortest-path intermediate has more than one interface. Schur expansion onto the two end droplets gives

    delta_ch,p = C_(k,p) U rho^(kp)+O(U rho^(kp+1)),
    C_(k,p)=2 sum_shortest_paths product_intermediates 1/(I(S)-1)>0. (6)

Reflection makes the diagonal entries equal and the shortest off-diagonal paths have a common sign. The p=2 gap thus scales as U rho^(2n+4). Equation (5), or an enclosed finite-matrix gap, should be used for a literal finite budget. Equation (6) explains its physical origin.

[proved] For rho(k+1)<1, projecting the charge eigenvector equation outside the two end droplets gives

    zeta=rho(k+1)/[1-rho(k+1)],
    0<=1-d_(p,j)<=zeta^2,
    r_p >= (1-zeta^2)/2.                                (7)

The first bound follows because the complementary diagonal energy is at least 2U, hopping norm is at most w_vac(N_p-1), and e_p<=U. Every occupied bond is present at either end droplet. Reflection divides the end-droplet weight equally and proves the last inequality.

[proved] Choose phi_h=1/b and let u_p=|+(-(N_p-1)phi_h)>. A positive left compensator is

    M=[r_1^-1(I-|u_2><u_2|)+r_2^-1(I-|u_1><u_1|)]/
          sin^2(phi_h/2).

It obeys <u_p|M|u_p>=1/r_p. The right compensator is its conjugate. If zeta<=1/2,

    ||M|| <= M_h=4/[(1-zeta^2)sin^2(1/(2b))].           (8)

A smaller Hermitian interpolant, C8 (8a), has the same expectations and is used in the numerical diagnostic. Positivity of a perturbing port observable is not required.

[proved] Each port-source projection onto the charge and helix ground factors is r_p times 1/r_p times the bare core source. This compensates physical endpoint occupation, not the unequal core amplitudes 1/sqrt(S) and q^(k-1)/sqrt(S). First-order port compression remains zero at both lengths. No sine-head factor occurs because the cut projector includes both X# and R#.

[proved] The physical marker move transfers only H_0 from the right cell to H_2 in the left. Ordered occupied labels, their biased vacuum code, and the occupation filter are unchanged. Its wrapper overlap is beta_v I, where

    beta_v=[sum_(j=0)^k chi_1({j}) chi_2({j,k+1})]^2>0.

Rotating the two helix qubits on the transferred site by diag(1,exp(i b phi_h)) and multiplying by exp(-i b phi_h) leaves overlap cos(phi_h/2)^(2b) on the untouched sites. The head overlap is A_b^2, with

    A_b=sum_(j=1)^(b-1) sqrt(2/b)sin(pi j/b)
                            sqrt(2/(b+1))sin(pi j/(b+1)).

Hence the exact first-order length hopping is

    J_len=u_mov beta_n,
    beta_n=A_b^2 beta_v cos(1/(2b))^(2b)>0.              (9)

The operator is the identity on the two complete 2^n-dimensional payloads. Omitting the typed-hole conversion makes this overlap zero.

## 5. Exact second-order bond at all four length pairs

[proved] In the legal source sector of bus r and hole number p define

    H_(p,r)=(H_ch,p-e_p) tensor I
             +s[sum_j D_j tensor h_(c,j,r)+I tensor(a_c,r+b_c,r)].

Its entire dimension is m_p*3(2n+2)2^n. The occupation-filter ground state is invariant under this operator, the source and charge motion: every gated filter term annihilates that same ordered product state separately. The selector head is also unchanged by a source. These factors can therefore be removed exactly from the source calculation, without an excitation cutoff.

[proved] Let H_(p,r)|alpha>=E_(p,r,alpha)|alpha>, retaining every eigenvector. Let A_(p,r,e,sigma),alpha,i be the amplitudes of the source

    (O_e chi_p) tensor F_(r,sigma)^e G|i>

in this eigenbasis. O_e is physical endpoint occupation. The reached endpoint-helix space consists of its ground vector and all one-flip states after untwisting and rotating the product state to |0>. On one flip its matrix is

    T_N=(1/2)L_path,N+|N-1><N-1|.

Write its eigenpairs as h_l,v_l. The helix source measure at the opposite endpoint has energies 0,J_h h_l and weights

    omega_0=r_p^-2,
    omega_l=|<u_p^perp|M|u_p>|^2 |v_l(0)|^2.             (10)

The other helix track stays in its ground state. No higher-flip state is reached, since one local one-qubit operator acts on the product vacuum and the parent conserves flip number.

[proved, exact formula] For each bus form

    Z_(alpha,beta;ij)=sum_(sigma,tau) c_(r,sigma,tau)
       A_(p,r,R,sigma),alpha,i A_(q,r,L,tau),beta,j,

where c has entries 1 for 00, exp(i phi_r) for +-, and exp(-i phi_r) for -+, with phi_1=pi/2,phi_2=pi. Then

    K2_(p,q) = -sum_r b_r^2 sum_(alpha,beta,l,m)
       omega_(p,l) omega_(q,m) Z_(alpha,beta)^dag Z_(alpha,beta) /
       [E_(p,r,alpha)+E_(q,r,beta)+J_h h_(p,l)+J_h h_(q,m)], (11)

with h_0=0, b_1=1,b_2=sqrt(eta). This is the coefficient of lambda^2. The denominator is positive. Its complete single-cell source sizes, including the reached helix, are

    binomial(n+2+p,p) * 3(2n+2)2^n * (n+3+p).

Different buses have no cross term. Two different cuts leave an outer cell excited, hence the second-order coefficient on three cells is the sum of adjacent bonds. Equation (11) retains every scalar, one-cell, hopping and non-Gaussian contribution at finite g.

[proved, unwrapped comparison] The bare bond is

    K_bare=-V0 sum_r b_r^2 B_r^dag(H_r tensor I+I tensor H_r)^-1 B_r,
    V0=q^(2n+2)/S^2.

This is C39's complete coefficient before multiplication by lambda^2. Its Gaussian hopping tends to J0[-Gamma(w)-i eta Gamma(wv)], where

    J0=2g^2 V0 D^(-2n-4).                                (12)

All four wrapped bonds obey s K2_(p,q)=K_bare+E_(p,q), with the finite bound below. In particular the factor 1/s cannot be omitted from the physical payload scale.

## 6. Finite equalization bound and actual mixed terms

[proved] In a pair-source sector write the full denominator A+sW, where A is the shifted charge plus helix operator and P its ground projector. Put Q=1-P, X=PZ,Y=QZ. Uniform bounds over the four length pairs are

    delta=min(delta_ch,1,delta_ch,2,J_h/[(b+1)(2b+1)]),
    L_c=2 M_c, L_b=2(k-1)h_*, a=1-|g|,
    x=sqrt(3(1+eta)V0), y=x M_h^2,
    t=s/delta.

Here ||X||=x, ||Y||<=y. Orthogonality of the three source charges and the two bus labels gives the value of x. Before projection the charge endpoint factors have norm at most one and each helix source norm is at most M_h. With zeta<=1/2 the compressed pair core is at least 2(1-zeta^2)(1-|g|)>=a. The bare pair has the same lower bound. Also ||W||<=L_c and the difference between compressed and bare pairs is at most zeta^2 L_b.

[proved] Let C=QWP, D_Q=Q(A+sW)Q, B=PWP. Block inversion gives

    s Z^dag(A+sW)^-1 Z
      =(X-s C^dag D_Q^-1 Y)^dag
        (B-s C^dag D_Q^-1 C)^-1
        (X-s C^dag D_Q^-1 Y)+s Y^dag D_Q^-1 Y.

Since W>=0, ||D_Q^-1||<=1/delta. If t L_c^2<=a/2, the Schur inverse has norm at most 2/a and differs from B^-1 by at most 2t L_c^2/a^2. Applying the triangle inequality to the displayed quadratic form, then the resolvent identity to B versus its bare value, yields

    ||E_(p,q)|| <= epsilon_eq,
    epsilon_eq = x^2 zeta^2 L_b/a^2
       +2t L_c^2 x^2/a^2
       +(2/a)[2t L_c x y+t^2 L_c^2 y^2]+t y^2.         (13)

This is an operator-norm bound, so noncommuting logical sources cause no extra assumption. Every constant is explicit in n and the chosen C39 parameters. The filter's large J_f does not enter (13), because its ground factor is exactly invariant in these source sectors. It does enter the full base-gap bound (3).

[proved] Equation (13) shows why finite source compensation is insufficient for exact equalization. The mean adjacency factors d_(p,j) change the core propagator at order rho^2; fast charge and helix source components add the remaining error. The ordered choices are rho first, then s/delta. Neither ratio may be replaced by a fixed n-independent small number.

[proved, length and payload decomposition] Write l_i=p_i-1 in {0,1}. The full scaled bond has the unique decomposition

    s K2_(p_i,p_(i+1))=K_bare+E_0+l_i E_L+l_(i+1) E_R
                                      +l_i l_(i+1) E_LR,
    E_0=E_11,
    E_L=E_21-E_11, E_R=E_12-E_11,
    E_LR=E_22-E_21-E_12+E_11.                            (14)

Their norms are at most epsilon_eq,2epsilon_eq,2epsilon_eq,4epsilon_eq. Normalized trace gives scalar length terms with the same bounds. Their traceless parts have at most twice those bounds. In particular the scalar density interaction c_22-c_21-c_12+c_11 is bounded by 4epsilon_eq and is generally nonzero. Constant and linear length terms can be retuned. This density interaction cannot. The nonscalar parts of E_L,E_R,E_LR are the leading length-payload couplings. They all carry physical prefactor lambda^2/s.

[proved, mixed marker-port selection rule] A marker move preserves the occupied sequence of each cell, including whether its token is vacuum or in bus r. The unperturbed resolvent has the same property. A single V creates an excited token on two cells and cannot return to the code through any sequence of marker moves. Thus the mixed u_mov lambda term vanishes separately in both orderings, and more generally every word with exactly one V vanishes. Pure marker excursions stay in the vacuum sector, whose core matrices act as the identity on data strings. Their effective operators are scalar on the payload, even if they leave the selected lengths temporarily.

[proved, first potentially nonzero higher mixed coefficient] On a fixed finite open cluster let R=QH0^-1Q after the common ground-energy shift, P the full code, T the unit marker perturbation, and V the unit port perturbation. The canonical third-order coefficient of lambda^2 u_mov is

    P(V R V R T+V R T R V+T R V R V)P
           -(1/2){P V R^2 V P,P T P}.                  (15)

Terms in the anticommutator with one V vanish by the same token rule. If ||R||<=1/Delta, ||V||<=v_c, ||T||<=t_c, its norm is at most 4v_c^2 t_c/Delta^2. The u_mov^2 coefficient is -P T R T P and has norm at most t_c^2/Delta; it is scalar on payloads. Neither (15) nor the length-dependent parts of (14) are asserted to vanish.

[proved, readout] Lifting the repeated conditional readout to occupied labels gives, at both lengths,

    P sum R_site P=(q^2/S)sum_(j=0)^(n-1) q^(2j) n_j.

It has no first-order length dependence. Its physical coefficient is (lambda^2/s) f_n with f_n=J0 mu S/q^2. Its unit coefficient has norm at most one per legal initialized cell, and at most k on the full occupied space. One V and one readout cannot return to the code. Marker-readout terms may be nonzero and are bounded on a cluster by 2|u_mov|(lambda^2/s)||T|| ||f_n R_read||/Delta. Readout self-terms are bounded by (lambda^2/s)^2||f_n R_read||^2/Delta. Two-V/readout terms have the bound (15) with T replaced by f_n R_read and u_mov replaced by lambda^2/s. These terms belong in the higher-order budget.

[proved, finite-ring qualification] Alternating the sign of excited tokens on alternate cells removes odd port orders on an open or even cell chain. This transformation is unavailable on an odd periodic cell ring. A winding word with M ports can visit every cut, so a blanket assertion that every finite ring has only even powers would be wrong. The local second-order results and the vanishing one-V mixed terms do not use that assertion.

## 7. Closing the all-n schedule without suppressing the narrow scale

[proved, inherited finite payload input] Fix a C39/C42 terminating rational witness. The nonresonance proof and root certificate are inputs to this wrapper audit, not rederived here. At zero readout the centers are theta_j=-atan(eta r_j), with r_j asymptotic to kappa_j v^(4j). Their readout slopes obey s_j asymptotic to d_j u^(4j) after choosing u and then v. The schedule enforces successive center and slope ratios below 1/32, and then chooses positive mu small enough that

    delta_payload >= min(|theta_(n-1)|,mu s_(n-1))/4 >0. (16)

A dovetailed root-isolation search avoids waiting forever on a trial tuple with an exact coincidence. D and g are then finite. Use an actual positive margin Delta from the accepted brackets after these choices; it can be smaller than (16). The fixed tuple u=1/64,v=1/1024 is only a small-n witness, not an all-n formula.

[proved, explicit inherited finite budgets] One available C42 budget uses r0=1/16384, p0=2n+4, B0=128 and a positive normalized symbol tolerance tau0 from rational root brackets:

    D>=max(2/r0,8n B0(1+eta)/(tau0 r0^(p0+1))),
    0<g<=1/4, g^2<=tau0/[192n(1+eta)D^p0].

A further decrease of g controls the non-Gaussian conversion remainder relative to J0; decreasing lambda alone does not do that. Obtain a remaining positive normalized band tolerance tau_band and minimum velocity v_* from the finite-g quadratic bond. These are finite computable numbers. For example C42's conservative normalized velocity bound is

    v_* >= (t_*/2)^n Delta^(2n-1)/(B+1)^(n-1),
    t_*=det Gamma(w)/(tr Gamma(w))^(n-1),
    B=1+2 tr Gamma(w)+2eta tr Gamma(wv)+mu,

provided the stated hopping perturbation budget is retained. The core's physical bare velocity is J0 v_*. Keeping the factor J0 in the wrapper budget is essential.

[proved, mixed momentum choice] A conservative finite list suffices at the free reference point. Let h range over integer vectors on the 2n payload endpoints with sum h=0 and sum |h|<=4. Include every equation

    2pi m nu + sum_a h_a p_a = 0 mod 2pi, m=1,2,       (17)

on I=[1/4,3/4]. This list contains all mixed length-bilinear/payload-bilinear resonances, pure length umklapp, and extra relations that need not actually be relevant. There are at most H_n=(4n+1)^4 such h: represent each by four choices from zero and the 4n signed endpoint labels, padding by zero. For each h and m there are at most m+1 roots in I. Hence F_n=5H_n bounds the number of forbidden points.

[proved] Include the two interval endpoints, sort the forbidden set and choose the midpoint of a longest component. Its width is at least 1/[2(F_n+1)]. An irrational displacement smaller than one eighth of that width still leaves distance at least 1/[8(F_n+1)] from every forbidden point and the endpoints. Since (17) has slope 2pi m, this gives the explicit margin

    delta_mix >= pi/[4(F_n+1)] >0.                     (18)

Roots just outside I cause no problem because the chosen point is at least the same distance from its endpoints. The actual density can be chosen by finite interval enclosures using the certified payload points. Payload endpoint errors must also be budgeted, since up to four enter (17). Exact marker number forbids length pairing, and exact payload U(1) forbids payload number-changing vertices. Classifying higher-dimension vertices as irrelevant and retaining that conclusion under interactions belongs to the assumed Gaussian stability step.

[proved, sufficient wrapper tolerances] Choose a fixed dimensionless length hopping j_l>0, for example j_l=1, and define

    delta0=min(1,Delta,delta_mix),
    epsilon=min(J0 tau_band/(128n),
                J0 v_* delta0/(512n),
                j_l delta_mix^2/[10^7(b+1)^2], j_l/512).     (19)

Shrinking this positive number further is allowed. A norm-epsilon two-cell operator changes the projected nearest-neighbor quadratic symbol by at most 16n epsilon: expand in the trace-orthogonal fermion bilinears, bound each coefficient by four times the operator norm, and bound each row sum. Thus (19) leaves room for the finite root brackets and the narrow velocity scale. The nonscalar length couplings in (14) and the scalar density interaction are small compared with both retained velocity scales. This last statement is a finite coefficient bound. It is not a proof that interacting momenta or stiffnesses stay fixed.

[proved, explicit ordered wrapper choices] Set x,a,L_c,L_b as in Section 6, and choose

    z=min(1/2,sqrt(epsilon a^2/(2x^2 L_b))),
    rho=z/[2(k+1)(1+z)], U=1,
    M_h0=16/[3sin^2(1/(2b))], y0=x M_h0^2,
    F=2L_c^2 x^2/a^2+4L_c x y0/a+2L_c^2 y0^2/a+y0^2.

Compute delta from (5) and the helix bound, and choose

    0<s<=delta min(1,a/(2L_c^2),epsilon/(2F)).           (20)

Then zeta<z, the adjacency part of (13) is <=epsilon/2, and its remaining part is <=t F<=epsilon/2. Thus ||E_(p,q)||<=epsilon at all four length pairs. All quantities are explicit functions of the accepted lacunary witness and n. Equations (5),(19),(20) also show how much smaller s must become when the narrowest momentum margin collapses. No unproved uniformity of the constants in (6) is used.

[proved, nonzero transport bound] If a_p=rho/(N_p-1), the minimum-amplitude bounds used for (5) imply

    beta_v >= (k+1)^2 a_1^(2k) a_2^(4k)/(m_1 m_2)>0.

Together with (9) this supplies an explicit positive lower bound on beta_n if an enclosed finite diagonalization is inconvenient. Choose

    u_mov=(lambda^2/s)j_l/beta_n.                       (21)

The positive beta_n may be small, but lambda is chosen after it. No coefficient varies with L.

[proved, filling and scalar retuning in the reference model] At fixed cell number, let e_pay be the per-cell ground energy of the common reference payload Hamiltonian, including its length-independent scalar and readout. The free length energy per cell is e_len(nu)=-2j_l sin(pi nu)/pi. The physical energy density has the form

    [kappa+e_pay+e_len(nu)]/(b+nu).

Its derivative vanishes at the chosen nu if

    kappa=(b+nu)e_len'(nu)-e_len(nu)-e_pay.              (22)

Since e_len''=2pi j_l sin(pi nu)>0 in the selected interval, this is a strict interior minimum. The scalar e_pay can be bracketed to any prescribed accuracy from finite open payload chains, with boundary-energy error bounded by a fixed multiple of the bond norm divided by the chain length. Thus (22) does not require an oracle for a thermodynamic energy. The marker chemical potential receives the correction (lambda^2/s)kappa. A finite ring chooses an integer cell count; its density approaches this target with an error of order (b+1)^2/L. This is a filling statement for the reference model. The tangent-subtracted numerator in the energy difference is at least (sqrt(2)pi j_l/2)(nu-nu_target)^2 on I. A perturbation with local norm epsilon changes the physical energy density by at most epsilon/b. The squared delta_mix restriction in (19) therefore also leaves a finite energy barrier against shifting the minimizing density by delta_mix/(32pi), allowing several such error terms. The reference density stays interior under the bounded second-order corrections. Establishing the required bound for the microscopic thermodynamic remainder remains part of the infrared input. The huge common scalar in K_bare must be included in e_pay, not ignored.

[proved, finite-cluster choice of lambda] For a fixed finite open cluster let Delta be the positive gap obtained from (4), v_c bound the sum of unit port terms, and A_c bound all retained order-lambda^2/s perturbations, including j_l T/beta_n, readout and (22). For instance each cut port has norm at most 3(1+sqrt(eta))M_h0^2, each unit marker bond norm at most one, and each readout cell norm at most k|f_n|. Cell counts multiply these bounds. Then

    q_c=||W|| <= |lambda|v_c+(lambda^2/s)A_c.

When q_c<=Delta/8 there is a canonical low-space effective matrix whose error beyond PWP-PWRWP is at most 32 q_c^3/Delta^2. One finite proof is to write its invariant subspace as a graph X over P. With D0=QH0Q the equation is

    X=D0^-1[-QWP-QWQ X+X PWP+X PWQ X].

On ||X||<=2q_c/Delta this is a contraction: its Lipschitz constant is at most 2q_c/Delta+4(q_c/Delta)^2<1. It gives ||X+D0^-1QWP||<=5q_c^2/Delta^2. Conjugating PWP+PWQ X by (I+X^dag X)^(1/2) makes it Hermitian and adds at most another 12q_c^3/Delta^2. The displayed 32 bound is conservative. This proves the finite-matrix statement without an extensive perturbation assumption.

[proved] One sufficient nonzero lambda satisfies

    lambda <= min(s v_c/A_c, Delta/(16v_c),
                  epsilon Delta^2/(512s v_c^3)),
    lambda^2 <= epsilon s Delta/(2 A_c^2).              (23)

The first two inequalities imply q_c<=2lambda v_c<=Delta/8; the third controls the cubic remainder relative to lambda^2/s. The last controls the second-order terms containing two order-lambda^2/s insertions. Mixed second-order terms containing just one V vanish by Section 6, including a marker chemical potential. The bounds allocate at most epsilon to these discarded pieces in scaled units. On an odd ring one may instead keep its winding terms in the same finite-cluster bound.

[assumed, precise boundary of this result] The constants v_c,A_c grow with cluster size. Equation (23) is not to be imposed with the whole ring as the cluster and then extrapolated. The missing thermodynamic statement must replace it by local elimination with volume-independent control in the appropriate interaction norm, preserve the auxiliary threshold and physical momentum gluing, and establish Gaussian stability relative to J0 v_* and the mixed margins. These are the remaining infrared assumptions. They are substantial, especially for the lacunary schedule.

## 8. Physical ring gluing and numerical checks

[proved] Identity transport gives the reference cover (x,l,z) with x in Z_L, M cell bits l_i in {0,1}, sum l_i=r=L-bM, and M payload states z_i. Its physical quotient is

    (x,l,z) ~ (x+b+l_0,S l,S z).

The dimension is (L/M) binomial(M,r) 2^(nM). Let the common payload ring have energy E_a and translation label m_a, with wavefunction chi_a(Sz)=exp(2pi i m_a/M)chi_a(z). Fourier transforming x and removing the gap-coordinate gauge gives length hopping with flux 2pi ell/L. Jordan-Wigner diagonalization of that hopping yields

    E=kappa M+E_a-2J_len sum_(j=1)^r cos(q_j+2pi ell/(LM)),
    q_j=2pi[I_j+(r-1)/2]/M,
    sum_j I_j+r(r-1)/2+m_a+ell=0 mod M.                 (24)

The half-integer shift is the hard-core ring boundary condition. The final congruence is the quotient constraint, so both the gap and payload translation labels enter. This unitary finite decomposition proves completeness. It applies to any common translation-invariant payload bond, including the full finite-g bond, not just an XX payload. At finite wrapper ratios (14) spoils exact separation, so (24) is the common-bond reference spectrum rather than the exact microscopic spectrum.

[numerical, source calculation] audit_numeric.py reconstructs each positive C39 occupied bond from the physical 46-letter rules, using attributed local copies of C39's model.py and literal.py. C32's mobile.py is copied as wrapper_tools.py for charge diagonalization, stable source block diagonalization and helix measures. Parameters are a=4,D=.2,g=.4,q=.8,w=.8,v=.65,eta=.6. These diagnostic points do not implement the lacunary witness or (19). Every source level in (11) is retained. The eight runs take about 15 seconds with one BLAS thread. The table gives full operator norms in units s K2; no scalar subtraction is used for these errors.

| n | rho | s/min(delta_ch) | max ||sK_pq-sK_11|| | max ||sK_pq-Kbar_pq|| | scalar E_LR |
| :-- | :-- | :-- | :-- | :-- | :-- |
| 1 | .16 | .03 | .02824377 | .02863682 | -.00619062 |
| 1 | .16 | .003 | .00306282 | .00342529 | -.00069091 |
| 1 | .10 | .03 | .02853931 | .02867629 | -.00622421 |
| 1 | .10 | .003 | .00330111 | .00343035 | -.00069266 |
| 2 | .16 | .03 | .01427237 | .01440654 | -.00312597 |
| 2 | .16 | .003 | .00159029 | .00172320 | -.00034838 |
| 2 | .10 | .03 | .01437501 | .01442649 | -.00313429 |
| 2 | .10 | .003 | .00167447 | .00172575 | -.00034870 |

[numerical] Kbar is the exact charge-ground compressed-core resolvent with its compensator source equalized. Its maximum distance from K_bare is .00045153,.00016330 at n=1 and .00017811,.00006848 at n=2 as rho decreases .16 to .10. These approach the expected quadratic missing-adjacency correction. The finite-s scalar correction dominates the displayed runs. The pair-source dimensions range from 480^2 to 7560^2 per bus, contracted by spectral sums without forming a dense pair matrix. The largest single-cell eigensystem before the helix convolution is 1080. The block-graph residual is below 5.8e-13.

[numerical, resolution limitation] The bare traceless norms are only 1.97e-8 at n=1 and 4.30e-10 at n=2. Some computed nonscalar differences of the full source matrices are affected by floating-point cancellation and ill-conditioned high blocks. The saved tiny mixed-payload norms are diagnostics, not certified relative hopping measurements. The large scalar length interaction is resolved. The all-n small-error conclusion follows from (13),(19),(20), not from an extrapolation of these numerical ratios.

[numerical, cell ED] With n=1, rho=.20,s=.2 and J_f=433.952933, the complete legal initialized core, every vacancy configuration and every occupation-filter state give matrices of dimensions 1728 and 4320 at lengths 4 and 5. Sparse shift-invert diagonalization finds the two-dimensional zero code in both, with next levels .01475806347 and .00003363095421. Code residuals are below 3.4e-16, code-subspace errors below 1.6e-13, and eigenpair residuals below 1.1e-13. The endpoint helices factor exactly; their complete matrices at lengths 4,5,6 have dimensions 16,32,64 and unique kernels, with gaps .07612046749,.04894348370,.03407417371, above their analytic bounds. The selector head is a separate sine path. Illegal sectors and wrong occupations are excluded by the proved operator inequalities, not by this ED.

[numerical, correction recorded] An initial cell subspace assertion failed because the degenerate ARPACK eigenvectors were not orthogonal and V V^dag was used as a projector. QR orthogonalization of their span fixes the check; the original code residual was already below 4e-16. The final reported subspace comparison uses that orthonormal basis.

[numerical, first-order maps] Direct sparse codeword contractions at n=1,2 give typed-hole transport beta_v times I_4,I_16 with errors below 1.4e-16, and zero overlap without hole conversion. The explicit helix product-state contraction agrees with cos(1/(2b))^(2b) to 4.2e-15. Every computed PVP is zero. These use the actual biased C39 vacuum code, not a replacement uniform payload.

[numerical, length-fluid and gluing ED] transport_ring.py constructs the physical marker graph, including cyclic payload relabeling when a marker crosses the origin. It adds the actual common C39 finite-g second-order bond, with its tiny traceless part rescaled to norm .3 for numerical visibility and its floating-point Hermitian part taken. Every eigenvalue agrees with the independent momentum decomposition (24) at

| n | M | r | L | complete projected dimension | maximum spectral error |
| :-- | :-- | :-- | :-- | :-- | :-- |
| 1 | 3 | 1 | 13 | 104 | 5.2e-15 |
| 1 | 4 | 2 | 18 | 432 | 1.7e-14 |
| 2 | 3 | 1 | 16 | 1024 | 2.4e-14 |

[numerical, scope] These are full spectra of the projected physical marker graph and single-cell ED of the stated invariant microscopic sectors. No full many-cell tensor product of 4464-state sites, no finite-coupling microscopic mobile-chain central charge, and no lacunary-scale spectrum was computed. The free length fluid is verified as the leading projected operator. Its survival at fixed nonzero coupling is assumed in the infrared step.

[proved and numerical, selector checks] bounds_checks.py computes the shortest-path coefficients by exact Fraction arithmetic. At n=1,2,3,4 the p=1 coefficients are all 2; the p=2 coefficients are 1/2,1/4,1/8,1/16. No extrapolation of that pattern is used. A floating-point check of the selector secant inequalities through length 500 at each of those n agrees with (1). The wrong-length coefficients t_sel eta_b-U are 1.69357336,1.39840043,1.30289697,1.25125262. The illegal-selector margin is U=1. The proof of (1) covers every length, not just that scan.

[numerical, independent source consistency] consistency_checks.py compares the bare source contraction with the separate C39 pair-eigenbasis implementation at n=1,2, with errors below 6.3e-17. The physical bond restrictions agree with the gauge matrices to 2.0e-16; charge commutators vanish. Each individual legal core term annihilates the code to 4.8e-17 and has no eigenvalue below -1.8e-15, the roundoff scale. The repeated readout compression agrees with its geometric formula. Direct inverses of the 96- and 240-dimensional charge-core source blocks at n=1 agree with the stabilized spectral calculation to 8.1e-14. These checks do not rely on resolving a tiny many-cell hopping splitting.
