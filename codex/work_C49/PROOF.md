# C49. Quasi-local dressing, a finite-volume spectral window, and the missing infrared estimate

[proved, status convention] A tag covers the following paragraphs and displays until the next tag. Imported results are tagged assumed here even where the cited files give proofs. The theorems below have full proofs. Numerical checks are separate. No assertion about the microscopic central charge is tagged proved.

[proved, outcome] The product-cell reference permits an explicit exponentially accurate local dressing. Its error is uniform in the number of cells in an interaction norm. A second, independent estimate makes the low-energy spectral error quadratic in the remaining leakage. Neither estimate is an o(1/M) result as M tends to infinity at fixed couplings. A sufficient ground-relative spectral comparison is stated and proved in Section 7. Its hypotheses are not established for this chain. Input S remains open.

## 1. Imported chain and conventions

[assumed, finite construction] Use C39/C42's fixed-n core and C44's mobile wrapper. The physical dimension is 4464, the selected lengths are b=n+3 and b+1, and the complete selected code per ordered cell has dimension 2^(n+1). On the rooted tensor product of M cell spaces,

    A = sum_j A_j,    A_j P_j=0,    A_j >= Delta Q_j,
    Q_j=1-P_j,        P=product_j P_j,    N_Q=sum_j Q_j.       (1)

The cell space is a direct sum over lengths. Delta is the complete positive wrapper gap, including the small charge tunneling gap. It is not the large penalty U. The bounded-interaction and invariant malformed-sector qualifications in C45 Section 3 remain in force.

[assumed, C45 reduction] C45 PROOF Sections 3 and 4 give an exact unitary U0 and an interaction decomposition

    U0 H U0* = A+B+R,       [B_X,P_X]=0,
    P_X=product_(j in X) P_j,
    ||R||_a <= epsilon = C_SW lambda^3,
    ||B||_a <= d = 2 lambda v0+2 lambda^2 G,
    G=f0+48v0^2/Delta,
    C_SW=2^16(v0 G/Delta+v0^3/Delta^2),                     (2)

for lambda<=min(v0/G,Delta/(192v0)). Here

    ||Phi||_a = sup_j sum_(X contains j) exp(a|X|)||Phi_X||.

Supports are connected intervals in ordered-cell coordinates. Constants and the exponent a are fixed before increasing M. The instantaneous retained bandwidth is

    W=lambda^2 c_n,     c_n=J0/s>0.                          (3)

The reference slow speed is W v_* with v_*>0 supplied by C42 and the length band. This is a reference speed, not a proved dressed microscopic velocity.

[assumed, physical locality contract] The rooted representation identifies an interval of r consecutive cells with the physical string between its outer markers. Local cell interactions preserve the total length on that interval and fix its outer markers. C44's port and marker-transfer rules have this property. They preserve exact payload U(1), total physical length, marker number and cyclic gluing. The cell inverse and every operation below preserve the same symmetries.

[proved, scope] Estimates for the rooted tensor product do not automatically compare different marker numbers on one fixed physical ring. Section 6 gives the precise additional energy offset needed when restricting to a physical sector. No stability assumption about the optimizing marker density is inserted into the proof.

## 2. Elementary interaction estimates

[proved, local pinching and inverse] For a self-adjoint term C_X write

    C_X^d=P_X C_X P_X+Q_X C_X Q_X,
    C_X^o=P_X C_X Q_X+Q_X C_X P_X,    Q_X=1-P_X.

Both maps are contractions in operator norm. For pinching this follows by averaging C_X with (P_X-Q_X)C_X(P_X-Q_X). For the off-diagonal part it follows from the block matrix with entries 0,T*,T,0, whose norm is ||T||<=||C_X||. Let

    R_X=Q_X (A_X|Q_X)^(-1) Q_X,
    L_X(C_X)=R_X C_X P_X-P_X C_X R_X.

Then

    L_X(C_X)*=-L_X(C_X),
    [L_X(C_X),A_X]=-C_X^o,
    ||L_X(C_X)||<=2||C_X||/Delta.                            (4)

The commutator follows from A_X P_X=0 and A_X R_X=Q_X. The reduced inverse is bounded by 1/Delta. Only this lower bound is used, not an upper bound on A_X.

[proved, loss-of-weight commutator bound] For 0<h<a,

    ||[C,D]||_(a-h) <= (4/h)||C||_a||D||_a.                 (5)

Terms on disjoint supports commute. For a fixed marked site j, split the overlapping pairs into those with j in X and those with j in Y. The first contribution is bounded by

    2 sum_(X contains j) e^((a-h)|X|)||C_X||
         sum_(Y intersects X) e^((a-h)|Y|)||D_Y||
    <= 2||D||_a sum_(X contains j) |X|e^(-h|X|)e^(a|X|)||C_X||
    <= (2/h)||C||_a||D||_a.

The second is the same with C,D interchanged. We used |X union Y|<=|X|+|Y| and r e^(-hr)<=1/h. Multiple terms with the same support may be kept as declared summands. Their absolute sum defines the bound.

[proved, convergent conjugation bound] Spend h/k of weight at each of k commutators. Equation (5), monotonicity in the weight and k!>=(k/e)^k give

    ||ad_S^k C||_(a-h)/k!
       <= (4e||S||_a/h)^k ||C||_a.                          (6)

Thus, if zeta=4e||S||_a/h<1, the BCH series of bounded interactions converges absolutely in the weaker norm. For A, use the bounded commutator [S,A] rather than ||A||. The identity

    e^(ad S) A=A+[S,A]+sum_(k>=1) ad_S^k[S,A]/(k+1)!

is valid on the domain of A. At finite volume a bounded S with bounded [S,A] preserves that domain, as follows by integrating the commutator of e^(tS) with A. The formulas can also be checked at finite cell-length cutoff and passed to the limit using the bounds (4)-(6). No upper spectral cutoff occurs in a constant.

## 3. Explicit finite iteration

[proved, theorem D] Suppose (1)-(2) hold. Define

    Lambda=d+2 epsilon,       h=64e Lambda/Delta,
    m=floor[a Delta/(128e Lambda)],
    a_j=a-jh,                epsilon_j=2^(-j)epsilon.        (7)

If Lambda=0 the conclusion is immediate with the identity. Otherwise there is an explicit product U_m of m quasi-local unitaries such that

    U_m U0 H U0* U_m* = A+D_m+O_m,
    [D_(m,X),P_X]=0,
    O_(m,X)=O_(m,X)^o,
    ||O_m||_(a/2) <= epsilon_m,
    ||D_m||_(a/2)+||O_m||_(a/2) <= Lambda.                  (8)

For m=0 these statements mean local pinching of B+R and no additional transformation. The retained operator K_m=P(A+D_m)P is quasi-local on the logical chain. In the common norm a/2,

    ||K_m-P(A+B)P|| <= 2 epsilon.                           (9)

The norm in (9) is an interaction norm after code compression, not an extensive operator norm.

[proved, construction and induction] Start with

    D_0=B+R^d,    O_0=R^o.

Their combined a-norm is at most Lambda. At step j choose S_j=L(O_j) term by term, and set

    T_j=sum_(k>=1) ad_(S_j)^k D_j/k!
             +sum_(k>=1) [k/(k+1)!] ad_(S_j)^k O_j,
    D_(j+1)=D_j+T_j^d,     O_(j+1)=T_j^o.                 (10)

The exact transformed Hamiltonian is A+D_(j+1)+O_(j+1). Indeed [S_j,A]=-O_j, and the coefficient of ad_(S_j)^k O_j after cancellation is 1/k!-1/(k+1)!=k/(k+1)!.

Write o_j=||O_j||_(a_j) and d_j=||D_j||_(a_j). If d_j+o_j<=Lambda, (4),(6) imply

    zeta_j <= 8e o_j/(Delta h)=o_j/(8 Lambda)<=1/8,
    ||T_j||_(a_(j+1))
       <= zeta_j(d_j+o_j)/(1-zeta_j) <= o_j/7.              (11)

Consequently o_(j+1)<=o_j/2 and

    d_(j+1)+o_(j+1)<=d_j+2o_j/7<=d_j+o_j<=Lambda.

This proves the induction. The product is U_m=e^(S_(m-1))...e^(S_0), and a_m>=a/2. The correction to K starts with PR^dP and then PT_j^dP. Its norm is at most epsilon+sum_j o_j/2<=2epsilon, proving (9). A useful accompanying bound is

    sum_(j=0)^(m-1) ||S_j||_(a/2) <= 4epsilon/Delta.         (12)

[proved, exponential error] Since floor(x)>=x-1,

    epsilon_m <= 2epsilon exp[-(log 2)a Delta/(128e Lambda)]. (13)

This is an explicitly constructed prethermal remainder. No optimality among all possible dressings is claimed. At fixed n and fixed wrapper parameters, Lambda=O(lambda), so the exponent grows as a positive constant times 1/lambda. If (7) gives m=0, only the original finite-order estimate is obtained. Taking lambda smaller makes m arbitrarily large.

[proved, what was special] Only transitions from the local all-code space to its complement are inverted. Excited-to-excited resonances are left in D_j. In particular [D_j,N_Q]=0 is neither assumed nor concluded. There is no denominator involving differences between two excited energies, or between two logical eigenvalues. This is why (4) remains available throughout the finite iteration. Extending the iteration indefinitely at fixed Lambda would exhaust the locality weight; (11) is not a convergence proof at infinite order.

## 4. Lieb-Robinson locality and the mobile coordinates

[proved, finite-time locality lemma] For a self-adjoint interaction G define

    J_mu(G)=sup_x sum_(X contains x) |X| e^(mu diam X)||G_X||.

For disjoint supports X,Y the Heisenberg evolution generated by G obeys

    ||[tau_t(A_X),B_Y]||
       <= 2||A_X||||B_Y|| |X| exp[-mu dist(X,Y)]
                  [exp(2|t|J_mu(G))-1].                   (14)

For a time-dependent G replace |t|J_mu by integral J_mu(G(s))ds. Here is the path estimate. Remove by an interaction-picture rotation the terms acting internally to the current support. The integrated commutator inequality iterates over chains of overlapping interaction supports Z_1,...,Z_k, with Z_1 meeting X and Z_k meeting Y. Each commutator costs 2||G_Z||. Choose the first marked site in X, then an overlap site at each step. Summing each next support and its possible next site costs at most J_mu. The sum of the diameters along a chain reaching Y is at least dist(X,Y), giving exp[-mu dist(X,Y)]. Ordered time integrals give |t|^k/k!. The initial commutator bound is 2||A_X||||B_Y||. Summing k>=1 gives (14). This argument also proves the time-dependent version with the time-ordered integrals bounded by the kth power of the integrated norm divided by k!.

[proved, ordered-cell locality] Choose 0<mu<a/2. Since r e^(-(a/2-mu)r) is bounded, J_mu(iS_j) is bounded by a fixed multiple of ||S_j||_(a/2). Equation (12) therefore gives a volume-independent Lieb-Robinson bound for the entire product U_m, regarded as a time-dependent evolution. C45's U0 is generated by range-three S_loc and obeys the same conclusion. This proves quasi-locality in ordered-cell distance, without a gap inside the retained chain.

[proved, physical support lemma under the contract in Section 1] Generic cell-local operators need not be physical-local when cells can be arbitrarily long. The present generators have an extra property that repairs this problem. Each S_(j,X) has a local P_X on one side in both of its two terms. If |X|=r, that side has total physical length ell in [br,(b+1)r]. Since the generator preserves total length, both sides have that same ell. Thus S_(j,X) vanishes in all other length sectors. Its physical realization is supported between the two outer markers, an interval with at most (b+1)r+1 sites. The same argument applies to the two generators in U0. The inverse A_X^(-1) does not change total length.

[proved, quantitative lift] For a translation-covariant interval decomposition write s_r for the norm of the term on r ordered cells, keeping extra term labels if necessary. The cell interaction norm is bounded using sum_r r e^(alpha r)s_r. Split the physical operator by ell, combining all internal marker patterns at fixed ell as a single operator. This is a direct-sum compression of the cell operator, so its norm is at most s_r. For a marked physical site there are at most ell+1 intervals of that length containing it, and at most r+1 allowed ell. Hence, if mu(b+1)<alpha,

    ||S^physical||_mu
       <= C_(b,alpha,mu) ||S^cell||_alpha,
    C_(b,alpha,mu)
       =2(b+2)e^mu sup_(r>=1) r exp[-(alpha-mu(b+1))r].      (15)

Indeed (r+1)((b+1)r+1)<=2(b+2)r^2. The remaining factor r is absorbed by the displayed exponential reserve. The direct sum over marker patterns avoids an exponential multiplicity from enumerating those patterns. On rings one uses cyclic intervals and the same bound on their number; terms covering the whole ring satisfy the same estimate. The invariant sectors on which the generators vanish may be extended by zero. No marker-free generator is introduced.

[proved, conclusion and scope] Apply (15) with alpha=a/2 and, for example, mu=a/[8(b+1)], retaining a further reserve when applying (14). Equations (12),(14),(15) prove a physical quasi-local dressing on the marked sectors covered by the imported representation, with constants independent of physical ring length at fixed n. This support argument concerns the generators. It does not assert that every arbitrary Q_X D_X Q_X term has bounded physical support merely because it is cell local. It also does not prove ground-state selection among distinct physical sectors.

## 5. Leakage is quadratic in the low-energy spectral estimate

[proved, theorem F] Consider any decomposition

    H'=A+D+O,
    [D_X,P_X]=0,    O_X=O_X^o,
    J=||D||_0,    epsilon'=||O||_0,
    delta=Delta-2J>0,
    K=PH'P,        E=min spec K.                            (16)

For a finite number M of rooted cells, assume

    M(epsilon')^2 <= delta^2/8.                             (17)

Let k_i be the increasingly ordered eigenvalues of K, including multiplicities. For every k_i<=E+delta/8, the correspondingly ordered eigenvalue h_i of H' satisfies

    k_i-t_M <= h_i <= k_i,
    t_M=8M(epsilon')^2/delta.                               (18)

In particular

    |(h_i-h_0)-(k_i-E)| <= t_M.                             (19)

For an eigenvector psi=p+q at one of these energies,

    ||q|| <= [8 sqrt(M)epsilon'/delta] ||p||.                (20)

Finite cell spaces suffice. For direct-sum cells the same variational argument applies when their resolvents are compact; C44's linear wrong-length cost and finite multiplicity at each length give that property at finite M. A fixed physical length is finite-dimensional.

[proved, reference threshold, attributed to C45] The completely positive extension argument of C45 Section 4 yields

    A+D >= E I+delta N_Q.                                  (21)

For clarity, extend each code operator to an excited cell by its expectation in a fixed code density matrix. The product extension is positive and unital, hence sends K>=E to an operator >=E. On each support the difference from D_X is bounded below by -2||D_X||Q_X. Since Q_X<=sum_(j in X)Q_j, summing gives (21). This cites the settled reference comparison, not a new stability theorem for H'.

[proved, square-root-volume leakage] Put T=QOP. If X and Y are disjoint, then

    P O_X O_Y P=0.

O_Y creates at least one defect in Y from the all-code state. O_X cannot remove that defect when its support is disjoint. Thus

    ||T||^2
       <= sum_(X intersects Y) ||O_X||||O_Y||
       <= sum_j [sum_(X contains j)||O_X||]^2
       <= M(epsilon')^2.                                  (22)

The possible overcounting of a pair by several sites only enlarges this upper bound. The argument holds for entangled code states and for moving defects. No independence of logical states is assumed.

[proved, form bound and the complementary threshold] For any normalized psi, writing r_X=||O_X||,

    |<O>_psi|
       <=2 sum_X r_X ||P_X psi|| ||Q_X psi||
       <=2 [(sum_X r_X)(sum_X r_X<Q_X>_psi)]^(1/2)
       <=2epsilon' sqrt(M<N_Q>_psi).

Young's inequality gives, for every u>0,

    O >= -u N_Q-M(epsilon')^2/u.                            (23)

Take u=delta/2 in (21). With Q=1-P and B_Q=QH'Q,

    B_Q >= [E+delta/2-2M(epsilon')^2/delta]Q
         >= [E+delta/4]Q.                                  (24)

This is only a finite-volume statement under (17). It must not be promoted to a uniform dressed-defect gap.

[proved, Feshbach comparison and eigenvalue count] For every real z<=E+delta/8 the inverse of B_Q-z exists and

    F(z)=K-z-T*(B_Q-z)^(-1)T,
    0<=T*(B_Q-z)^(-1)T<=t_M I_P.                           (25)

The Schur factorization of H'-z by invertible triangular matrices has diagonal blocks F(z) and B_Q-z. Its negative index therefore equals the negative index of F(z). This can also be seen by completing the square in q. Restricting trial spaces to P gives h_i<=k_i by min-max. If z<k_i-t_M, then F(z)>=K-z-t_M has at most i eigenvalues below zero when indexing from i=0. Hence H' has at most i eigenvalues below z, proving h_i>=k_i-t_M by continuity. This controls multiplicities and excludes extra levels below the relevant ordered levels; it is more than an approximate-eigenvector assertion. Since h_i-k_i and h_0-E both lie in [-t_M,0], their difference lies in [-t_M,t_M], proving (19). The Q-component equation gives q=-(B_Q-h_i)^(-1)Tp and (20).

[proved, uniform ground-energy estimate] Minimizing delta x-2epsilon'sqrt(Mx) over x>=0 in the stronger form preceding (23) gives, without (17),

    E-M(epsilon')^2/delta <= min spec H' <= E.               (26)

It is an energy-density estimate of order (epsilon')^2/delta. At fixed nonzero epsilon' it is not a Casimir or central-charge estimate. Subtracting the two ground energies does not remove the unknown state-dependent shifts in (18).

## 6. Application and the finite spectral window

[proved conditional on the imported reduction] Apply theorem F to (8). If Lambda<=Delta/4, one can take

    delta=Delta-2Lambda>=Delta/2,
    epsilon'=epsilon_m=2^(-m) C_SW lambda^3.                (27)

Using a smaller delta than the actual threshold only weakens the bounds. The exact retained Hamiltonian K_m differs in weighted interaction norm from C45's complete second-order code operator by at most 2C_SW lambda^3. This includes scalar length interactions, forward terms and every instantaneous many-body correction generated by the specified finite iteration. It does not say that K_m is a free Hamiltonian.

[proved, spacing criterion] Compare (19) to a reference oscillator scale v/M, where v=W v_*=lambda^2 c_n v_*. A relative tolerance eta>0 follows from

    M <= sqrt(eta delta v/8)/epsilon_m,
    M epsilon_m^2 <= delta^2/8,                             (28)

and the reference levels in question must lie within delta/8 of E. The first inequality is the useful spectral restriction. For fixed n, the construction certifies a window growing at least as a constant times

    lambda^(-2) exp[(log 2)a Delta/(128e Lambda)]            (29)

as lambda decreases, provided the other finite-sector and reference-spectrum hypotheses are met. The word window refers to this error budget, not to a proved Gaussian spectrum of K_m. The positivity of the reference speed v does not prove that the interacting K_m has that speed. At fixed nonzero lambda the right-hand side of (28) is a finite certified size. The upper bounds do not show failure beyond it, but cease to certify the required accuracy there.

[proved, physical-sector qualification] Let C be a conserved total-length, charge and cyclic-root sector with P_C nonzero, and let E_C be the ground energy of K restricted to C. Define d_C=E_C-E>=0. The unrestricted inequality (21) remains true after restriction. If d_C<=delta/16, then under (17)

    Q_C H' Q_C >= E_C+3delta/16.

For k_i^C<=E_C+delta/16, the same proof has denominator at least delta/8, so (18)-(20) hold unchanged in C. More generally the denominator is the explicitly reduced threshold from (24), minus d_C and the chosen energy window. This is the necessary finite-sector hypothesis. Establishing d_C=O(W/M) for the interacting K_m at the desired physical density, and comparing competing marker numbers, is part of the unresolved infrared and ensemble problem. The free-reference rounding estimate cannot silently be substituted for it.

[proved, symmetry] Every generator, compressed correction and residual respects the microscopic U(1), physical translation and marker/length constraints whenever these are represented on the chosen sector. The iteration does not introduce symmetry-violating pairing. Momentum nonresonance of the quadratic reference alone does not provide bounds on the renormalized matrix elements needed as M increases.

## 7. A weaker statement that would suffice at the logical spacing

[proved, theorem S_low, conditional comparison] The following is a sufficient replacement for exact block diagonalization of the entire code band. It isolates the missing estimate rather than assumes a separated extensive band.

Let K_M have ground energy E_M and a known reference low-energy spectrum. Let v>0 be fixed and q_M=v/M. Suppose a volume-dependent unitary puts the microscopic Hamiltonian into blocks A_M,T_M,B_M relative to an embedding of the full logical space. A scalar shift s_M is allowed. Suppose

    B_M >= (s_M+E_M+gamma)I,       gamma>0 independent of M. (30)

For z in a fixed neighborhood of s_M+E_M, define

    F_M(z)=A_M-z-T_M*(B_M-z)^(-1)T_M.

Assume there are eta_M tending to zero such that, as quadratic forms on the logical space,

    |<u,[F_M(z)-(K_M+s_M-z)]u>|
       <= eta_M <u,(K_M-E_M+q_M)u>                         (31)

uniformly for the z used below. For each fixed R>0 and every ordered reference level k_i-E_M<=R q_M, the corresponding microscopic level obeys

    |h_i-(s_M+k_i)| <= eta_M(R+1)q_M,
    |(h_i-h_0)-(k_i-E_M)| <= eta_M(R+2)q_M.                (32)

These errors are o(1/M). The result can be imposed separately in each conserved momentum/charge sector.

[proved, proof] Put x=z-s_M-E_M and G_M=K_M-E_M>=0. Equation (31) gives

    (1-eta_M)G_M-x-eta_M q_M
       <= F_M(z) <= (1+eta_M)G_M-x+eta_M q_M.              (33)

For sufficiently large M, eta_M<1 and the comparison energies below lie below (30). The negative-index identity used in (25) sandwiches h_i-s_M-E_M between

    (1-eta_M)(k_i-E_M)-eta_M q_M
    and
    (1+eta_M)(k_i-E_M)+eta_M q_M.

This proves the first inequality in (32). For i=0 the bound is eta_M q_M, which gives the second. The needed z window may be taken to contain |x|<2(R+1)q_M for each fixed R. Increasing R handles every fixed scaling-energy window without controlling the entire high-energy code band.

[proved, why this is weaker and why an unweighted projection is inadequate] Equation (31) permits errors growing with logical excitation energy. It requires neither an inverse on the overlapping full code spectrum nor a unitary eliminating all transitions at all energies. But the weight is on the entire logical quadratic form. A small norm of Pi_low R Pi_low alone is insufficient because Pi_low R(1-Pi_low) can shift the low levels through a Schur complement. For the two-by-two matrix [[0,g],[g,delta]], the low-low remainder is zero while the lowest level shifts by (delta-sqrt(delta^2+4g^2))/2. This is why completeness or an integrated complementary self-energy must be part of any projected criterion.

[assumed, connection to c] If K_M is independently proved to have the n+1 Gaussian modes with their compact charge lattice, and the sectorwise scaling spectrum and required local observables transfer, (32) preserves that spectrum. Under an independently established unitary CFT interpretation, this supplies the same mode central charges, c=n+1. A spectral comparison alone is not a construction of a conformal operator algebra. If one instead extracts c from the vacuum Casimir energy, the scalar s_M must also be controlled through order 1/M, with the dressed velocities retained. An arbitrary s_M can hide that coefficient. C46 Section 8 and C51 address the current/thermal identification separately.

[proved, the missing implication] The bounds already proved give a Schur self-energy no larger than 8M epsilon_m^2/delta. On the logical ground state the right side of (31) is eta_M v/M. Thus these bounds yield at best an error parameter of order

    M^2 epsilon_m^2/(delta v),                             (34)

rather than eta_M tending to zero. A scalar subtraction may remove a common bulk shift, but none of (4)-(25) bounds the remaining state-dependent part by logical excitation energy with a coefficient tending to zero. There is no proved implication from the cell gap and Lieb-Robinson locality to (31). This is the unresolved estimate for this route. It is not a counterexample to the chain's proposed phase.

## 8. Why the customary shortcuts do not close the gap

[proved, a bare complementary compression need not stay massive] Take one auxiliary two-level system per cell,

    h=[[0,t],[t,Delta]],    P_j=|0><0|,    Delta>0, t!=0,
    e_-=(Delta-sqrt(Delta^2+4t^2))/2<0,
    alpha=(1+Delta/sqrt(Delta^2+4t^2))/2 in (0,1).

For H=sum_j h_j, let psi be its product ground vector, E_H=M e_-, and p=||P psi||^2=alpha^M. The normalized vector Qpsi/sqrt(1-p) has energy

    <Qpsi,H Qpsi>/(1-p)
       = E_H+(-M e_-)p/(1-p).                             (35)

To verify this, expand Q=1-P, use Hpsi=E_H psi and PHP=0. The numerator is E_H(1-2p), so subtracting E_H after division gives -E_H p/(1-p). By the variational principle,

    0<=min spec(QHQ)-E_H<=(-M e_-)alpha^M/(1-alpha^M).       (36)

This tends to zero exponentially. Nevertheless the physical excitation gap is sqrt(Delta^2+4t^2), and a product of onsite rotations diagonalizes H exactly. Tensoring with any spectator logical chain preserves the example. Therefore the bare Q resolvent near the true ground need not be uniformly massive, and exponentially small overlap of two vacua is not by itself an obstruction to a local unitary. A fixed finite-order approximate dressing leaves this same issue with a smaller nonzero residual rotation.

[proved, no isolated full logical band from a cell gap] The retained code has an extensive many-body bandwidth at nonzero W. The auxiliary threshold is ground-relative and O(Delta). These facts do not give a separating contour enclosing the complete code spectrum with a volume-independent gap. Excited logical energies can coincide with auxiliary energies. The local inverse (4) avoids those denominators; ordinary spectral flow of the whole code band would need an additional separation hypothesis not present here.

[proved, standard low-energy spectral-flow estimate loses uniform locality] If one instead follows a ground-state or fixed scaling-window projector of a gapless reference, the available separating gap is of order v/M, when such a projector is separated at all. An exact inverse-frequency filter for a generic gapped spectral-flow argument has Fourier transform proportional to 1/omega above that gap gamma_M. Its time-kernel L1 norm is at least a constant times 1/gamma_M because |Fhat(omega)|<=||F||_1. Wherever it agrees smoothly with 1/omega, its first absolute time moment is at least a constant times 1/gamma_M^2 by differentiating the Fourier transform. Hence the generic locality bound based on that filter and (14) is not uniform as gamma_M tends to zero. This does not exclude special gapless paths with extra cancellations. Those cancellations would be new input, not a consequence of the cell gap.

[proved, local indistinguishability is absent] For any retained logical occupation n_a there is an operator supported on one selected finite cell whose compression is n_a. It has different expectation values in two code states. Thus the product code does not obey a condition of the form P O_X P=c_X P for every bounded local region X. Gapped stability arguments requiring that condition, or a uniformly gapped infrared subspace along the path, cannot be applied by checking only (1). This statement is a direct hypothesis test, not an attribution of a new theorem to Bravyi-Hastings-Michalakis.

[assumed, literature boundary] The supplied local review of Chen-Lucas-Yin 2303.07386 and Yin-Lucas 2209.11242 describes a general obstruction to demanding exact quasi-local block diagonalization, and prethermal alternatives. Their primary statements are not stored in the supplied text sources and were not retrieved under the no-network instruction. No theorem-specific exponent or no-go for this particular chain is claimed here. Equations (4)-(15) are an independent constructive proof under the imported product-cell hypotheses. An exceptional exactly decouplable chain is not ruled out by the general warning, as (35) illustrates.

## 9. Checks, provenance and remaining work

[numerical, checks passed] checks.py checks synthetic cells with code dimension two, total dimensions 9,27,81,243, and all 4,8,16,32 retained levels respectively. Every leakage, form, threshold, eigenvalue and eigenvector bound passes. The one-step BCH formula agrees with exact matrix exponentiation to 1.24e-16. The product-rotation identity agrees through eight cells to 3.34e-16. The scalar recurrence and weight budget pass. These are controls of the new operator inequalities, not simulations of the 4464-state chain. Analytic proofs establish the bounds for arbitrary cell count; finite floating-point checks do not replace those proofs. Reproduce with `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B checks.py`. The results are in checks.json.

[proved, provenance] The local inverse and ground-relative reference threshold are C45's ingredients. The loss-of-weight iteration, the quadratic-leakage spectral window and the weighted Schur comparison are derived above. C46's already established thermal-reference theorem is not re-derived. C48 pursues resummed Euclidean matching, C50 the logical RG, C51 central-charge extraction, and C52 exact redesign. Their progress was re-read while this work was in progress; no result from those unfinished tasks is assumed proved.

[assumed, unresolved chain statements] To finish the requested result one still needs (i) ground-relative, energy-weighted self-energy control such as (31), or a different all-scale matching theorem, (ii) the interacting logical RG and central-charge identification, and (iii) the fixed-physical-length sector/marker-density comparison for the dressed operator. The finite iteration proves a prethermal dressing and the finite-volume estimates (18)-(29). It does not prove S or c=n+1 at fixed nonzero couplings.
