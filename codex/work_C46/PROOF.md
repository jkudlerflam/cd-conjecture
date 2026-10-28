# C46. Rigorous infrared reductions and the remaining gaps

[proved, status convention] A statement tagged proved below is either an abstract theorem with its proof here, or an explicitly conditional implication with the conditions stated. Imported construction results are tagged assumed, even when the cited earlier file supplies a proof. Numerical statements describe only the controls in checks.py. Heuristic statements are proposals, not conclusions of an estimate. Definitions carry the tag of the argument in which they are used. Each tag applies to its displayed formulas and continuation prose until the next tag or heading.

[proved, outcome] The results below do not prove the requested microscopic decoupling or central charge theorem. They prove a uniform thermal comparison for the locally block-diagonal reference, identify an obstruction to a raw Grassmann expansion even at zero temperature, give explicit second-order time bounds, strengthen the finite momentum test to include the length fermion, and prove the thermodynamic coefficient of the quadratic reference. Sections 3, 6 and 8 identify the additional statements needed for the original chain.

## 1. Inputs, scales and scope

[assumed, imported finite construction] Use C39 PROOFS Sections 1-7, C42 PROOF Sections 1-6, and C44 PROOFS Sections 1-8. At fixed n the physical alphabet has size 4464, the selected physical cell lengths are b=n+3 and b+1, and their complete retained space has dimension d_P=2^(n+1). On rooted ordered cells, the unperturbed Hamiltonian is A=sum_j A_j with A_j P_j=0 and A_j>=Delta Q_j, Q_j=1-P_j. Delta is the complete wrapper gap. It includes the charge tunneling splitting, not just the charge penalty U.

[assumed, imported explicit bounds] C42 defines

    p=2n+4, q=(wu)^2, w=1/4,
    S_q=sum_(j=0)^(n+1) q^(2j), V0=q^(2n+2)/S_q^2,
    J0=2 V0 g^2 D^(-p), W=lambda^2 J0/s.

Its rational root certificate gives delta_pay>0, tau>0 and

    v_cert=(t_*/2)^n delta_pay^(2n-1)/(B+1)^(n-1)>0.

These are speeds in J0-normalized cell units before the wrapper factor. With r=1/16384 and B0=128, its symbol errors obey

    E_D <= 2n B0(1+eta)/(D r^(p+1)),
    E_g <= 48n(1+eta)g^2 D^p,
    epsilon_ng <= 3(1+eta)(2+64n^2)g^2 D^p.                (1)

The last bound is on the complete second-order non-Gaussian bond in J0 units. These estimates can be made small in the order D, then g. Lowering lambda does not improve epsilon_ng.

[assumed, imported wrapper constants] From C44 Section 6, with k=n+2,

    a_c=1-|g|, x^2=3(1+eta)V0,
    zeta=rho(k+1)/(1-rho(k+1)),
    M_h=4/[(1-zeta^2) sin^2(1/(2b))], y<=x M_h^2,
    delta_fast=min(delta_ch,1,delta_ch,2,J_h/[(b+1)(2b+1)]).

Here delta_ch,p has the explicit lower bound

    delta_ch,p >= w_c [rho/(k+p-1)]^(2kp)/[binomial(k+p,p) kp].

With L_c=2M_c, L_b=2(k-1)(K_c+D+12) and t=s/delta_fast,

    epsilon_eq = x^2 zeta^2 L_b/a_c^2
       +2t L_c^2 x^2/a_c^2
       +(2/a_c)(2t L_c x y+t^2 L_c^2 y^2)+t y^2            (2)

bounds every full scaled pair difference sK2_(p,q)-K_bare, including the scalar. The complete Delta is the minimum specified in C44 (3)-(4), using the charge bound above. The source and conversion constants in this note can therefore be expressed in finite cell parameters. None is uniform in n.

[assumed, imported local reduction] C45 Section 3 gives an exact identity on its rooted cell representation,

    U H U^*=D+R, D=A+sum_X B_X,
    [B_X,P_X]=0, P_X=product_(j in X)P_j,
    H_P=PDP=lambda^2(PFP-PVQA^(-1)QVP),
    ||R||_a <= epsilon_R=C_SW lambda^3.                    (3)

The interaction norm is sup_j sum_(X containing j) exp(a|X|)||B_X||. Taking

    v0=max(1,||V||_a), f0=||F||_a,
    G=f0+48v0^2/Delta,
    C_SW=2^16(v0 G/Delta+v0^3/Delta^2),
    lambda<=min(v0/G,Delta/(192v0))                        (4)

gives J=sup_j sum_(X containing j)||B_X||<=4lambda v0<=Delta/48. The proof of (3)-(4) is in C45, and is not repeated or upgraded to an infinite-order transformation here. In particular [D,Q_j]=0 is not assumed. The bounded-perturbation and invariant malformed-sector qualifications in that proof remain in force.

[proved, scope of subsequent applications] The pressure theorem in Section 2 applies to a tensor product of rooted cells before the total physical length and cyclic quotient are imposed. A theorem about this tensor product is not automatically a theorem about the full physical Gibbs ensemble. A separately stated cell-coordinate Gaussian reference is used in Section 7. The physical quotient, equilibrium marker density and illegal-selector sectors require the transfers described in Sections 7-8. This distinction is necessary because a product partial trace need not preserve a fixed total-length subspace.

## 2. A uniform thermal decoupling theorem for the reference D

[proved, theorem T] Let each cell have a finite-dimensional code P_j of rank d_j and an orthogonal complement Q_j. Allow Q_j to be infinite-dimensional. Let A_j vanish on P_j, preserve Q_j and obey A_j>=Delta_j Q_j. Assume the interaction B_X is bounded, supported on X, and block diagonal for P_X versus 1-P_X. Define

    J_j=sum_(X containing j)||B_X||,
    t_j(beta)=Tr_(Q_j) exp[-beta(A_j-2J_j)],
    Z_P=Tr_P exp(-beta H_P), Z_D=Tr exp(-beta D).

If all t_j(beta) are finite, then for every finite number M of cells,

    1 <= Z_D/Z_P <= product_j [1+t_j(beta)/d_j].            (5)

For identical finite cells with complement dimension d_Q and Delta>2J,

    0 <= f_P(beta,M)-f_D(beta,M)
      <= beta^(-1) log[1+(d_Q/d_P)e^(-beta(Delta-2J))].      (6)

Here f=-(beta M)^(-1)log Z. The constants are independent of M and the finite-dimensional estimate holds for every beta>0. No gap within the code Hamiltonian is assumed.

[proved, product extension] On one code matrix define a map into the full cell algebra by

    E_j(O)=O on P_j, (Tr O/d_j) I on Q_j.

It is unital and positive, and its tensor product is positive. One direct verification of the latter is to decompose the full space by the set S of cells in Q. In this block E=product E_j is the normalized partial trace over the code factors in S, tensored with identities on their Q factors. Partial trace is positive since, on a vector v, its expectation is the sum of the nonnegative expectations on v tensor an orthonormal basis. Thus no tensor-positivity assertion is being inferred merely from positivity of an arbitrary map.

[proved, comparison operator] E_X(P_X B_X P_X) agrees with B_X in the P_X block and has norm at most ||B_X||. The norm claim follows either from the preceding partial-trace formula, or from positivity and the inequalities -||B_X||I<=P_X B_X P_X<=||B_X||I. Both operators are block diagonal for P_X. Consequently

    B_X-E_X(P_X B_X P_X) >= -2||B_X||(1-P_X).

The commuting projections obey 1-P_X<=sum_(j in X)Q_j. Summing gives

    D >= E(H_P)+sum_j (A_j-2J_j Q_j) =: K.                (7)

Also [D,P]=0, so the full trace contains the positive summand Z_P. This proves the lower inequality in (5).

[proved, trace exponential order] If self-adjoint finite matrices H>=K, their increasingly ordered eigenvalues obey lambda_r(H)>=lambda_r(K). To see this, in the minimum over r-dimensional subspaces of the maximum Rayleigh quotient, every Rayleigh quotient for H is at least the corresponding one for K. Therefore Tr exp(-beta H)<=Tr exp(-beta K). The same argument with the variational eigenvalues applies to semibounded operators with compact resolvent and trace-class heat kernels. The operator exponential itself is not being asserted to be order monotone.

[proved, partial-trace inequality] Let S contain r erased code factors and let d_S=product_(j in S)d_j. Put H_S=d_S^(-1)Tr_S H_P. Choose an orthonormal eigenbasis v_i of H_S and any orthonormal basis a of the erased factors. Twice using scalar convexity of exp gives

    exp(-beta <v_i,H_S v_i>)
      <= d_S^(-1) sum_a exp(-beta <v_i,a|H_P|v_i,a>)
      <= d_S^(-1) sum_a <v_i,a|exp(-beta H_P)|v_i,a>.

The second inequality follows by expanding v_i tensor a in an eigenbasis of H_P. Summing i yields

    Tr exp(-beta H_S) <= Z_P/d_S.                          (8)

[proved, completion of theorem T] In the block with Q precisely on S, K is H_S on the remaining code factors plus the mutually commuting A_j-2J_j on the Q factors. Equations (7)-(8) imply

    Z_D <= Tr exp(-beta K)
      <= Z_P sum_(S subset {1,...,M}) product_(j in S)t_j(beta)/d_j.

The subset sum factorizes and proves (5). If A_j>=Delta Q_j, t_j<=d_Q exp[-beta(Delta-2J)], giving (6). This proof includes defect hopping and transitions changing the number of defects inside a local Q_X block.

[proved, ground and thermal corollaries] When Delta>2J and cells are finite, E(H_P)>=E_P I by positivity, so (7) implies D>=E_P+(Delta-2J)sum_j Q_j. The code attains E_P, hence D and H_P have the same ground energy. Moreover the right side of (6) is o(beta^(-m)) for every fixed m>0. Indeed log(1+x)<=x, and exp(beta delta) grows faster than every power by its power series. Thus any already-established low-temperature power coefficient of H_P passes unchanged to D. This is a statement about D, not D+R.

[proved, direct-sum extension] The same conclusions hold whenever t_j(beta)<=C exp(-beta delta) for beta>=beta_0, with d_Q exp(-beta delta) in (6) replaced by C exp(-beta delta). The proof above only used the heat trace, so an infinite complement is allowed. The bounds can first be proved with spectral cutoffs on A_j and then passed to the full trace, or proved directly using (7) and the variational eigenvalues.

[assumed, cell counting for the wrapper application] In the well-formed rooted representation of C45, use C44's all-length cost gamma(1+|ell-b|), gamma=t_sel eta_b-U>0, for legal selector cells of length ell>=2 outside {b,b+1}. Treat ell=1 and the two selected lengths as exceptional finite spaces with excluded gap at least Delta, choosing Delta<=gamma. The number of states of length ell is at most d0^ell, d0=4464. The ell=1 treatment uses C45's penalized ## convention. It is not a new assertion about all possible illegal selector strings.

[proved, explicit thermal constant under that input] Put

    D_fin=d0+d0^b+d0^(b+1), beta_0=log(2d0)/gamma,
    C_th=[D_fin+b d0^(b-1)+2d0^(b+2)]/d_P.

For beta>=beta_0,

    t(beta)/d_P <= C_th exp[-beta(Delta-2J)].               (9)

To prove it, the exceptional spaces contribute at most D_fin exp(-beta Delta). For 2<=ell<=b-1, sum their dimensions by b d0^(b-1), and their cost is at least 2gamma. For ell=b+2+r, r>=0, the sum is

    d0^(b+2)e^(-3beta gamma)/(1-d0 e^(-beta gamma)).

Its denominator is at least 1/2. Multiply by exp(2beta J)/d_P and use Delta<=gamma. This proves (9). Under (4), Delta-2J>=23Delta/24, so Delta/2 is a valid weaker exponent.

[proved, comparison with the original rooted Hamiltonian] Suppose (3) holds on this same trace-class space. Since sum_X||R_X||<=M||R||_0, -M epsilon_R I<=R<=M epsilon_R I. The eigenvalue comparison used above gives |f_(D+R)-f_D|<=epsilon_R. Unitary invariance identifies f_H=f_(D+R). Thus

    |f_H-f_P| <= C_SW lambda^3
          + beta^(-1)log[1+C_th e^(-beta Delta/2)]         (10)

for beta>=beta_0 and the specified rooted ensemble. In finite cells replace C_th by d_Q/d_P and use (6) for every beta>0.

[proved, why (10) does not identify a central charge] At fixed lambda>0, beta^2 C_SW lambda^3 is unbounded as beta tends to infinity. Subtracting the two ground energies only changes the estimate to 2C_SW lambda^3 plus the thermal term. Thus (10) controls energy and pressure densities but does not control the coefficient of beta^(-2). This conclusion follows directly from the inequality; it does not require a proposed gapping mechanism.

[proved, range where the error is small] The quadratic reference's coefficient is A_n/W with A_n>0 from Section 7. Ignoring the separately bounded exponential thermal term, (10) gives a small relative thermal error if

    C_SW lambda^3 W << A_n beta^(-2).

With W=c_n lambda^2, this requires T=beta^(-1) much larger than a constant times lambda^(5/2). Linearization also requires T<<W e_lin, a constant times lambda^2. These inequalities have a nonempty interval for sufficiently small lambda at fixed n. They describe a shrinking comparison window as lambda decreases. They are not the T->0 limit at one fixed Hamiltonian, nor a proof that H_P itself has the Gaussian coefficient.

## 3. Why the requested general Grassmann theorem needs more hypotheses

[proved, first obstruction, thermal logical dependence] Consider one cell with two logical fermions, an even two-state auxiliary with q=0,1, and

    H_c=q(Delta+u n_1 n_2), Delta>0, u>0.                 (11)

Its complete code is q=0, with four degenerate states, and its excluded gap is Delta. Exact logical U(1) is preserved. At zero intercell coupling, the logical weights after tracing q are A=1+e^(-beta Delta) on 00,10,01 and B=1+e^(-beta(Delta+u)) on 11. Directly summing four weights gives

    <n_1 n_2>-<n_1><n_2>=A(B-A)/(3A+B)^2 !=0.             (12)

All flavor-off-diagonal and anomalous two-point functions vanish by the independent phase symmetries. Wick's identity for a Gaussian number-conserving state would make the left side zero. Thus no scalar auxiliary partition factor times a Gaussian logical state reproduces these tensor-factor observables at lambda=0. A theorem with non-Gaussian kernels vanishing as a positive power of lambda, uniformly for every beta>0, is false for this class of cells. It needs thermal terms or stronger factorization assumptions.

[proved, definition of the stronger norm issue] A raw external-field influence functional is obtained by integrating only the auxiliary at fixed logical coherent-state fields. At different times, the even Grassmann monomials

    N(t)=bar(psi_1)(t)psi_1(t)bar(psi_2)(t)psi_2(t)

are independent generators. Define the absolute coefficient norm by summing their ordered-monomial coefficients, integrating all their time arguments per unit time, and giving each external fermion leg a weight a>=1. This is the unweighted-time restriction of the type of fully summed kernel norm in C40 Section 3. Antisymmetric kernels with a 1/p! convention give the same ordered coefficient sum, since their permutations cancel that factorial. Adding positive temporal or spatial weights can only increase the norm.

[proved, zero-temperature raw-excursion obstruction] Add the local auxiliary flip lambda sigma_x to (11), and use auxiliary vacuum boundary conditions. In the vacuum-to-vacuum Duhamel expansion, the order-lambda^2 contribution consists of creating q=1 at s and removing it at t>s. It is

    lambda^2 integral_(s<t) ds dt e^(-Delta(t-s))
              exp[-u integral_s^t N(v)dv].                (13)

The sign of this expression in minus the action is immaterial for its absolute norm. At order lambda there is no vacuum term, so the logarithm has the same order-lambda^2 coefficient. This formula can be obtained without a continuous coherent-state assumption. With time step h, an occupied auxiliary slice has the exact logical matrix element multiplier 1+(e^(-hu)-1)N_l relative to its free coherent-state overlap. The N_l commute and square to zero, so their product is exp[(e^(-hu)-1)sum_l N_l]. Taking h->0 at fixed field degree gives (13). Two auxiliary flips supply its endpoints.

[proved, divergence calculation] Let T=t-s and consider exactly m external quartic insertions in (13). Its absolute coefficient, integrated over these insertions and over the excursion duration per unit starting time, is

    lambda^2 (u^m/m!) integral_0^infinity e^(-Delta T)T^m dT
       = lambda^2 u^m/Delta^(m+1).                        (14)

Integration by parts, starting from integral e^(-Delta T)dT=1/Delta, proves the factorial integral. The m insertion times are unordered in the exponential; equivalently ordered times occupy a simplex of volume T^m/m!. Endpoints can be integrated out first and give the same positive absolute integral. The contribution to the leg-weighted norm is therefore

    (lambda^2/Delta) sum_(m>=1) (a^4 u/Delta)^m.            (15)

It diverges whenever a^4 u>=Delta. The code, gap and exact U(1) have not changed. An added temporal exponential weight cannot repair divergence already present in (15). On a finite time interval of length L, the per-time expression has an additional factor (1-T/L)_+; monotone convergence as L increases gives (14). Thus no bound uniform in observation time is hidden by the infinite-line notation.

[proved, precise force and limit of this counterexample] The order-lambda^2 coefficient of this raw influence action is not an element of the requested summed-coefficient space when u>=Delta. Hence there is no analytic, or norm-asymptotic-through-order-two, decoupling theorem in that space for arbitrary logical-dependent gapped cells using this raw expansion. This does not prove nonexistence of every effective action, of every low-frequency resummation, or of the action for the particular C39/C44 cell. In fact the operator second-order effective Hamiltonian of this example is the finite polynomial

    -lambda^2/(Delta+u n_1n_2)
       = -lambda^2/Delta
           +lambda^2 u/[Delta(Delta+u)] n_1n_2.            (16)

The identity follows from (n_1n_2)^2=n_1n_2. It explicitly exhibits the cancellations/resummation that an uncontracted-time coefficient norm fails to capture. The one-cell gap and the operator norm bound on the resolvent are both consistent with (15).

[proved, an additional finite-temperature norm obstruction] For the model (11), integrating q with the thermal trace gives, in the same raw fields,

    F[N]=log[(1+r exp(-u X))/(1+r)],
    r=e^(-beta Delta), X=integral_0^beta N(t)dt.            (17)

Consider the ordinary complex function F_beta(z)=log[(1+r exp(-u beta z))/(1+r)] at z=0. Its closest singularities are

    z=-Delta/u +/- i pi/(u beta),
    R_beta=sqrt(Delta^2+pi^2/beta^2)/u.                    (18)

Indeed the logarithm's argument first vanishes at these points; inside that disk it is nonzero and has an analytic logarithm. Its derivative has a pole at each displayed point, so the radius is not larger. If F_beta(z)=sum c_m z^m, the absolute integrated coefficient at m quartic insertions is |c_m|, before the harmless division by beta. With leg weight a it is |c_m|a^(4m). This follows by expanding X^m and integrating m times; it can also be seen from the discrete formula where h^m m! binomial(beta/h,m) tends to beta^m. Therefore the summed norm diverges whenever a^4>R_beta. Taking u=2Delta makes this happen for a>=1 at every beta>pi/(sqrt(3)Delta).

[proved, interpretation of the thermal result] Every fixed field-degree coefficient in (17) tends to zero exponentially as beta grows, times a polynomial in beta. The sum over all degrees is nevertheless divergent in the stated norm at those beta. Thus adding a small scalar thermal error to a raw all-degree action estimate does not solve the problem. The order of these limits matters. This statement concerns the raw action, not its already-contracted physical pressure, which is finite and satisfies theorem T.

[proved, source-observable qualification] The tensor logical observable O and the code-projected observable P O P differ on excited states. Even when the auxiliary is completely independent of the logical Hamiltonian, repeated P O P insertions impose auxiliary projectors at their insertion times. Their auxiliary trace is not a single source-independent scalar multiplying the trace of the tensor logical O's. This follows by writing the tensor trace and leaving the intervening P's in its auxiliary factor. Any exact decoupling theorem for the code-projected observables in the literature question must specify source-dependent insertion kernels or dressed observables. Astra B's free-auxiliary theorem explicitly uses tensor-factor observables.

[heuristic, implication for the wrapped cell] C39's excited carrier propagation acts nontrivially on the logical register, and C44 adds charge, helix and length excursions. Its cell gap does not impose an external-field analyticity radius larger than one. Equations (14)-(16) explain why a gap-only polymer proof cannot be inferred from C42/C44 norm constants. A viable proof should resum same-cell logical operator products and retained contact contractions before estimating external legs, and separately control winding histories. No such resummation theorem for this wrapped cell is proved here.

## 4. What can be bounded in imaginary time for the actual source

[proved, theorem E, one vacuum excursion] Let A_Q>=Delta>0 on an excited space and let Z map a retained pair to it. Define

    K(t)=-lambda^2 Z^* exp(-t A_Q) Z, t>=0.

For any integer m>=0 and 0<=kappa<Delta,

    integral_0^infinity t^m e^(kappa t)||K(t)||dt
       <= lambda^2 ||Z||^2 m!/(Delta-kappa)^(m+1).         (19)

Also integral K(t)dt=-lambda^2 Z^*A_Q^(-1)Z, and, for real omega,

    ||lambda^2 Z^*[(A_Q+i omega)^(-1)-A_Q^(-1)]Z||
        <= lambda^2 ||Z||^2 |omega|/Delta^2.              (20)

[proved, proof] The spectral decomposition gives ||exp(-t A_Q)||<=exp(-Delta t). The gamma integral proved in Section 3 gives (19). Integrating each spectral eigenvalue proves the resolvent identity for the static term. Finally

    (A_Q+i omega)^(-1)-A_Q^(-1)
        =-i omega (A_Q+i omega)^(-1)A_Q^(-1),

and each inverse on the right has norm at most 1/Delta. This proves (20), with no restriction on the real frequency. For complex z inside |z|<Delta the corresponding denominator is Delta(Delta-|z|).

[proved, application conditional on C44 source data] Treat the two buses as an orthogonal direct sum in Z. C44's ground and complementary fast-source bounds give

    ||Z||^2 <= x^2+y^2 <= B_src,
    B_src=3(1+eta)V0(1+M_h^4).                            (21)

The first inequality follows since the two source components have orthogonal ranges. Use the complete pair excitation gap at least Delta in (19). There are at most two bonds incident on one ordered cell. Thus a spatial norm weighted by exp(mu|X|) satisfies

    ||K_2||_(mu,kappa,m)
      <= 2e^(2mu) lambda^2 B_src m!/(Delta-kappa)^(m+1).    (22)

For kappa=Delta/2 this is explicit in C42/C44's parameters. Its time integral is the complete finite-g second-order source formula in C44 (11), including every scalar and mixed length term.

[proved, scope of (22)] K(t) is an operator-valued vacuum transfer kernel, not a full raw Grassmann action at independently specified intermediate times. Passing its operator norm to arbitrary external-time coefficients is precisely the invalid step exposed by (14)-(16). Thermal winding contributions and connected multi-excursion histories are not included in (22). Equations (19)-(22) therefore prove fixed-order retardation and its moments, not the requested all-order decoupling theorem.

[proved, why the available coefficient estimates do not close it] C40's finite-order estimate contains a factor [16(r+1)]^r at perturbative order r. Its r-th root grows without bound. Therefore the sum of those upper bounds times any fixed positive coupling to the r-th power diverges. Divergence of an upper bound does not imply divergence of the actual coefficients, but it supplies no convergent majorant. C45's estimate concerns a fixed second-order generator and its convergent BCH exponential; it supplies (3), not a summation of all subsequent auxiliary histories. No geometric polymer activity bound follows from either statement.

[heuristic, missing decoupling lemma in usable form] The shortest sufficient replacement is an exact low-external-frequency action, including matched source insertions, whose norm is finite after operator/contact resummation and subtraction of retained propagation. It should have an error epsilon_dec(lambda)=o(W) in a coefficient norm with two time and space moments, uniformly in cell count, with thermal winding terms controlled separately in the pressure or in an appropriate source norm. Exponential decay at rate of order Delta would be sufficient but is stronger than the required moments. One must specify the external-field radius, the treatment of the infinite length sum and the stability of excluded physical sectors. Calling these operations a cluster expansion does not prove their convergence.


## 5. Exact charges and a strengthened finite momentum certificate

[proved, exact payload U(1)] The C39 on-site charge is the data bit plus the excited-token charge, with token charges 0,+1,-1. Neutral and charged hopping preserve both their species and the transported data label. Their Z factors and phases are diagonal in charge. Conversion |+,0><0,1| preserves total charge 1, and |-,1><0,0| preserves total charge 0. Their adjoints do the same. F_+ raises charge by one, F_- lowers it by one and F_0 is neutral. Each port product is F_+ F_-, F_- F_+ or F_0 F_0, so every port term is neutral. Every wrapper track is neutral; vacancy swaps transport the entire occupied label, marker moves transport a neutral hole, and the readout is diagonal. Thus each physical local term commutes with the sum of the on-site charges on its support. Summing proves [H,Q_pay]=0 on the unrestricted physical chain. The subgroup exp(i pi Q_pay) is its exact fermionic parity. This removes the lack of U(1) discussed for earlier candidates in the literature notes.

[proved, retained length charge] In the rooted direct sum define Q_len=sum_j(ell_j-b). Ports do not move boundaries and marker hopping changes two adjacent lengths by +1 and -1. All remaining terms preserve lengths. Thus [H,Q_len]=0 in this enlarged representation. In its code Q_len=sum_j l_j, l_j=0,1. After the usual occupation representation of the two-state length variable, it is conserved length-fermion number. This is a kinematic symmetry of ordered cells before fixing total physical length. It is not a new on-site U(1) charge of the original fixed-L physical alphabet.

[proved, preservation by covariant elimination] Spectral projections and resolvents of a Hamiltonian commuting with a charge commute with that charge. Products and commutators of charge-preserving operators do also. Hence the local generators of C45 preserve both rooted charges. An exact auxiliary trace with source fields transformed along with the charge is invariant as well: conjugation inside the traced auxiliary space leaves the trace unchanged. Expanding this invariant polynomial or formal logarithm in charge eigenmonomials kills every term with nonzero total charge. This proves the selection rule if such an action is constructed. It does not prove existence of its norm limit.

[proved, finite marginal momentum classification] Let p_i, i=1,...,2n, be the payload Fermi points and let the length points be +pi nu and -pi nu. A local quadratic vertex conserving each of the two numbers is either branch diagonal or has a momentum difference within one species. A quartic vertex with two creations and two annihilations has momentum equal to the outgoing pair sum minus the incoming pair sum. If both pairs are identical as multisets, normal ordering reduces it to densities and density products, with possible lower-degree contact terms. Every other such vertex has a nonzero pair mismatch when the full pair-sum condition below holds. This is algebraic momentum conservation, not a bosonization claim. A repeated Grassmann field at the same point vanishes; allowing repeated labels in the finite test is a stronger condition and avoids having to remove those cases.

[proved, an extra finite check if using an unrestricted multichannel theorem] C42 certifies different payload pair sums, but that test alone does not compare a pair sum with zero. For the enlarged list of Fermi points, the length pair (+pi nu,-pi nu) always sums to zero. Therefore a sufficient full-pair theorem also needs

    delta_zero=min_(i<=j) dist(p_i+p_j,2pi Z)>0.            (23)

This is not implied by abstract pair-sum separation. For example, the two points -alpha,+alpha have distinct repeated and mixed pair sums for generic alpha, but their mixed sum is zero. If one uses the separate length/payload conservation law, conversion of a length pair into a payload pair is forbidden and (23) is unnecessary. It is useful for fitting a stronger theorem statement with only one total U(1).

[proved, compatibility of (23) with C39's all-n construction, conditional on its asymptotic roots] At mu=0 C39's roots are theta_j+/-pi/2, with all theta_j<0, |theta_j|<1/8 and min_j|theta_j|=theta_min>0. A pair sum is theta_i+theta_j plus 0 or +/-pi. The first kind is at distance at least 2theta_min from 2pi Z; the other kind is at distance at least pi-1/4. Root shifts of less than theta_min/4 preserve a lower bound theta_min for every pair sum. The root continuity and arbitrarily small positive readout supplied by C39/C42 permit this choice simultaneously with their positive pair-difference margins. An additional rational acceptance test on |z_i z_j-1|, where z_i encloses exp(i p_i), therefore has terminating witnesses. Tighten the finite-D and finite-g tolerance to preserve that test. This adds a finite certificate to the schedule and no physical letter or interaction. It is not a claim that every previously accepted C42 tuple already included this extra test.

[proved, all-channel density selection] Suppose the payload list has both delta_pay and delta_zero. Set N=2n+2, Pairs=N(N+1)/2 and

    F_n=3 binomial(Pairs,2).

For each distinct pair of unordered pairs from the enlarged list, its sum difference is a*pi nu+c modulo 2pi, with integer a in [-4,4]. A comparison with a=0 is independent of nu. It is nonzero by delta_pay, delta_zero or the payload root separation, except for comparisons of the identical pair, which are not in the test. If a!=0, the zeros in I=[1/4,3/4] are separated by 2/|a|>=1/2. There are at most two, and the factor three in F_n is a safe upper bound. Thus at most F_n forbidden points occur in I.

[proved, explicit mixed margin] Insert the two endpoints of I and choose the midpoint of a longest complementary interval. Its length is at least 1/[2(F_n+1)]. Move this midpoint by less than one quarter of that interval's length if needed to choose a rational nu. Its distance from every forbidden point and from both endpoints is at least

    d_n=1/[8(F_n+1)].

A zero outside I is also at distance at least d_n because an endpoint lies in between. For a nonzero slope, the distance of a*pi nu+c from 2pi Z equals |a|pi times the distance to its nearest zero on the real line. Consequently the enlarged pair margin is at least

    delta_all=min(delta_pay,delta_zero,pi d_n)>0.           (24)

Root separation follows by comparing pairs with one common label. The pure length umklapp at nu=1/2 is among these exclusions. This proves finite nonresonance for all 2(n+1) branches, with an explicit conservative bound. If only symmetry-allowed quartics are tested, many exclusions can be removed.

[proved, robustness of that finite test] If every Fermi point changes by at most r, a pair difference changes in circular distance by at most 4r by the triangle inequality. Thus r<delta_all/8 retains at least half the margin. Likewise a nonzero slope persists under a C1 perturbation that changes it by less than half its magnitude. These facts describe the budget needed for a dressed Fermi surface. They do not bound the interacting self-energy or the dressed root displacement.

[assumed, available root perturbation bounds] C40 Section 4 and C45 (23) supply explicit C1 symbol tolerances using isolated patches, v_cert, a band separation d_*, an off-patch gap e_off and a derivative bound M1. C42's finite rational certificate provides the initial roots, while these patch estimates extend the finite argument beyond nearest-neighbor quadratic symbols. Their application to an interacting frequency-dependent self-energy still requires a constructive estimate.

[proved, velocities, and what is not certified] Transversality means a nonzero signed derivative at each root. Distinct roots do not logically imply different derivative magnitudes. Already for the n=1 reference

    epsilon(k)=-2A cos k+2B sin k+C,

the two root derivatives have equal magnitudes sqrt(4(A^2+B^2)-C^2): write the first two terms as a single shifted cosine and differentiate at its two crossings. Thus a demand that every absolute chiral speed be different is not generally satisfied, even in the simplest payload. C42 provides a strictly positive minimum, not a separation certificate for all speed magnitudes. Its singular congruence is used to locate roots and cannot supply physical speeds.

[proved, an optional repair for length versus payload speeds] At the chosen nu, length hopping j_l W has speed 2j_l W sin(pi nu). If avoidance of the finitely many payload speed magnitudes is desired, choose j_l in [1,2] outside the finite set v_i/[2sin(pi nu)] in normalized units. Such a choice exists since finitely many points cannot fill an interval. Retune the length chemical potential to keep nu. This proves a length-speed repair, not an all-n distinct-payload-speed theorem. If distinct velocities means only different nonchiral channels, the n=1 equality is not an objection; an explicit all-channel speed-separation certificate is still not supplied by C42. Equal speeds by themselves do not change the finite momentum test or the quadratic thermodynamic calculation.

## 6. Audit against the cited constructive RG frameworks

[assumed, precise source boundary] The local reviews ../../contrib/ir_control/README.md and ../../contrib/lit_ir_step/README.md describe Mastropietro-Porta 2108.13708 and Porta-Scola-Singh 2509.08665 as multichannel constructive fermion frameworks with simple Fermi points, number conservation, momentum nonresonance and weak short-range interactions. The first review attributes the multichannel reference beta-function result to Theorem 8.5 of 2108.13708 and its microscopic assumptions to Assumption 2.1. It describes a one-dimensional multiorbital model and transversal Fermi points for 2509.08665. These are imported descriptions, not primary passages independently checked here. No local primary text was located in the supplied directories, and the user forbade network access. No exact theorem-number application or unverified hypothesis list is used as a proof below.

[proved, logical limitation of a source audit from those descriptions] Matching selected items in a summary does not establish all hypotheses of a theorem whose full statement has not been checked. Moreover, none of the new estimates in Sections 2-5 constructs the complete action or controls its RG flow. Consequently this note cannot assert that either published theorem applies to the microscopic chain. The following table identifies the inputs actually supplied and the gaps, rather than treating an extension of a theorem as part of its statement.

| Input | Status for the wrapped chain |
| :-- | :-- |
| Exact charge conservation and parity | [proved] Section 5 proves payload U(1) on the full physical chain and length number on rooted cells. The earlier missing-U(1) objection is repaired. |
| Finite set of transverse free Fermi roots | [assumed] C42 supplies 2n payload roots and v_cert>0; the explicit free length dispersion supplies two more. |
| Quartic momentum nonresonance | [proved conditional] Section 5 upgrades the C39/C42 witness by the optional zero-pair test and a finite density choice, giving (24). |
| Different absolute velocities | [proved] The supplied certificates do not assert this, and n=1 has an equal left/right pair. The local reviews describe arbitrary velocities; an all-distinct condition is not assumed to be required. |
| Small complete interaction norm | [assumed for the needed action] C42/C44 and C45 give finite operator and local remainder norms, but not an exact, summed retarded action norm. Sections 3-4 show the distinction. |
| General many-body and retarded vertices | [assumed extension] The local summaries start from number-conserving two-body Hamiltonians. Their direct applicability to the resulting arbitrary-degree, time-nonlocal action has not been proved here. |
| Dressed Fermi points and counterterms | [assumed] A uniform self-energy/root bound at fixed physical couplings is missing. A theorem proved with independently tuned channel counterterms needs a matching argument for this restricted bond family. |
| Fixed physical length and cyclic gluing | [assumed transfer] C44 supplies the exact common-bond reference quotient. Passing an interacting rooted action and its observables to the physical Gibbs ensemble remains part of the application. |
| Low-temperature pressure or stress tensor | [assumed extra conclusion] The local reviews report correlation/response conclusions and beta-function control, not the thermal coefficient required here. No such coefficient is deduced from them in this note. |

[proved, distinction between a counterterm method and microscopic tuning] A running diagonal quadratic term is allowed by total U(1) and translation. Those symmetries alone do not set its coefficient to zero. To apply a theorem stated at prescribed Fermi momenta one must either realize its counterterms by allowed microscopic parameters or prove an inverse map to shifted momenta for the fixed microscopic model. The latter option can avoid individually tunable chemical potentials, but requires a self-energy and susceptibility bound. It is not supplied by the finite-matrix root certificate. This explains the counterterm row without asserting that either cited paper requires extra physical tunability.

[proved, what momentum support would give once roots are controlled] If a vertex involves at most four branch momenta and has mismatch at least delta_all/2, four external momenta each within delta_all/16 of their root cannot sum to a reciprocal lattice vector. Translation invariance therefore sets that low-patch kernel to zero. This follows because their possible displacement sum is at most delta_all/4. Thus detuning removes the non-density quartic support after the patches have shrunk sufficiently. It does not control the flow before that scale, the accumulated density interactions, or the complete sum of higher vertices.

[assumed, Gaussian vertex dictionary and stiffness control] C40 Sections 5-6 and C45 Sections 6-7 describe the compact-boson charge lattice, the lower free dimension three for the remaining non-Gaussian primaries, and the stiffness estimate giving D>=3sqrt(3/5)>2 at relative stiffness error at most 1/4. These are conditional Gaussian statements there. They require an established continuum dictionary and a bound on the total drift into that neighborhood. This note does not promote them to a microscopic constructive RG proof.

[heuristic, possible repair of the RG step] First construct a resummed massive-sector action in a norm that can be inserted into the multiscale determinant expansion. Then check the full primary statements and prove the needed ultraviolet initialization for arbitrary-degree retarded vertices, preserving Ward identities and the two charges. Use either a fixed-bare-parameter Fermi-surface theorem or an inverse counterterm argument with constants depending on the smallest speed. The finite detuning scale makes the preliminary RG interval finite at each n, but its constants can be arbitrarily large. The step is plausible and still unproved here.

## 7. A rigorous thermal coefficient and its central-charge meaning

[proved, theorem F, quadratic lattice reference] Let h(k) be a finite Hermitian Bloch matrix, continuous and piecewise C2 on the reciprocal circle. Suppose all zeros of its eigenvalues are isolated simple crossings, with no other zero eigenvalues. Denote their signed derivatives by u_i!=0, i=1,...,N_ch. In the thermodynamic limit of the free number-conserving fermion chain,

    f_0(beta)=e_0-(pi/12) sum_i |u_i|^(-1) beta^(-2)
                    +o(beta^(-2)).                        (25)

Twists of a bounded angle on the ring do not change this thermodynamic formula. They only shift the Riemann sum mesh by order 1/M.

[proved, free partition function] Diagonalizing each finite-ring one-particle matrix gives independent occupation numbers with energies epsilon_alpha(k). The trace for one occupation is 1+exp(-beta epsilon), so dividing log Z by M and taking its Riemann-sum limit yields

    f_0(beta)=-beta^(-1) integral_(BZ) dk/(2pi)
                           sum_alpha log(1+e^(-beta epsilon_alpha(k))),
    e_0=integral_(BZ) dk/(2pi) sum_alpha min(epsilon_alpha(k),0).

The identity log(1+e^(-beta e))+beta min(e,0)=log(1+e^(-beta|e|)) gives

    f_0(beta)-e_0=-beta^(-1) integral dk/(2pi)
                           sum_alpha log(1+e^(-beta|epsilon_alpha(k)|)). (26)

[proved, asymptotic evaluation] Choose disjoint monotone patches around the simple roots. Compactness bounds |epsilon| away from zero outside the patches, so that contribution to (26) is exponentially small. On patch i use e=epsilon_i(k), whose Jacobian j_i(e)=|epsilon_i'(k(e))|^(-1) is continuous and tends to 1/|u_i| as e->0. Set x=beta e. The integral on this patch, multiplied by beta, tends to

    |u_i|^(-1) integral_(-infinity)^infinity log(1+e^(-|x|))dx
       = pi^2/(6|u_i|).

For domination, j_i is bounded on a sufficiently small patch and log(1+e^(-|x|))<=e^(-|x|); use its integrability. The one-sided integral is sum_(m>=1)(-1)^(m+1)/m^2=pi^2/12, obtained by integrating the logarithm series. Absolute integrability of that series follows from sum m^(-2)<infinity. The sum's value follows from the Fourier series of x^2 on [-pi,pi], whose constant term is pi^2/3 and whose m-th cosine coefficient is 4(-1)^m/m^2 by two integrations by parts; evaluating at pi gives sum m^(-2)=pi^2/6, and subtracting twice the even terms gives the alternating sum. Combining patches with the prefactor in (26) proves (25).

[proved, application conditional on the finite-band input] Apply theorem F to the C42 quadratic symbol in cell coordinates, with its physical scale W, and append the free length dispersion -2W j_l cos k minus its chemical potential. There are 2n+2 simple crossings. Writing their normalized signed speeds as bar_u_i,

    f_G(beta)=e_G-(pi/(12W))
                     sum_(i=1)^(2n+2)|bar_u_i|^(-1) beta^(-2)
                        +o(beta^(-2)).                    (27)

For the two length crossings |bar_u|=2j_l sin(pi nu). In the finite list at each fixed n, every coefficient is finite. This is a theorem about the selected quadratic reference. The complete second-order Hamiltonian H_P contains non-Gaussian terms bounded by (1)-(2), so H_P is not identified with this quadratic model at fixed nonzero parameters.

[proved, Gaussian boson counterpart] Suppose an infrared theory has N independent nonchiral compact Gaussian modes, with positive right and left speeds v_(a,R),v_(a,L), a=1,...,N. Its nonzero oscillators of one chirality on circumference M have frequencies 2pi v m/M, m=1,2,... . Their thermal log partition function is -sum_(m>=1)log(1-exp(-2pi beta v m/M)). Dividing the free energy by M and taking the Riemann-sum limit gives

    (1/(2pi beta)) integral_0^infinity dk log(1-e^(-beta v k))
       = -pi/(12v beta^2).

To evaluate the integral, expand -log(1-e^(-x))=sum_(m>=1)e^(-mx)/m, integrate its nonnegative terms, and use the same sum of inverse squares. The logarithmic singularity at zero is integrable, which justifies the improper Riemann-sum limit. A positive definite compact zero-mode quadratic form contributes only O(log M) to log Z at fixed beta: each integer Gaussian sum is bounded above by a constant times 1+sqrt(M/beta) and below by one. A fixed finite number of such sums therefore vanishes after division by M. Thus

    f_G(beta)=e_G-(pi/12) sum_a
                 [v_(a,R)^(-1)+v_(a,L)^(-1)] beta^(-2).    (28)

Cutoff corrections, if present, require their own small-temperature estimate; (28) is the continuum Gaussian oscillator result.

[proved, normalization correction] Equation (28) becomes the question's -(pi/6)sum_a v_a^(-1) beta^(-2) only when v_(a,R)=v_(a,L)=v_a for each pair. Exact U(1) does not imply reflection or time-reversal symmetry. The two-bus reference has no imposed symmetry giving that equality at every n. One could define v_a as the harmonic mean of the two chiral speeds to rewrite the scalar coefficient, but that number would not in general be either propagation speed. The safe thermodynamic target is (28) with N=n+1.

[proved, canonical free-length ensemble] The thermodynamic free energy at fixed length-fermion density nu in (0,1) is the Legendre transform of the free grand-canonical pressure. Here is a proof avoiding a finite-size ensemble assumption. Divide momentum space into K bins and replace the continuous dispersion in each by a constant e_j. The energy of any occupation configuration changes by at most M times the maximum binwise error, so its free energy per cell changes by at most that error. Let M_j be the number of momenta in bin j. For the piecewise constant problem the canonical trace is exactly

    sum_(sum r_j=r) product_j binomial(M_j,r_j) exp(-beta e_j r_j).

There are at most (M+1)^K summands. Its log differs from the largest summand's log by at most K log(M+1). Integral comparison for sum_(l=1)^m log l gives log(m!)=m log m-m+O(log(m+1)); applying it to the three factorials yields, uniformly in 0<=r_j<=M_j,

    log binomial(M_j,r_j)=M_j h(r_j/M_j)+O(log(M+1)),
    h(x)=-x log x-(1-x)log(1-x), h(0)=h(1)=0.

Thus at fixed K, M->infinity gives the minimum of sum_j w_j[e_j x_j-T h(x_j)] subject to sum_j w_j x_j=nu, with w_j the bin fractions. The discrete constraints approach this compact feasible polytope: round K-1 occupation counts and adjust the last, or distribute the bounded rounding error over bins with available capacity. Continuity gives the same limiting minimum. Its strict convexity for T>0 and its diverging endpoint derivatives give a unique chemical potential mu with x_j=(1+exp(beta(e_j-mu)))^(-1). Substitution proves the finite-bin Legendre formula. Refining bins removes their uniform dispersion error and proves the continuum formula. The argument is uniform for nu in compact subintervals of (0,1).

[proved, canonical low-temperature coefficient] For the cosine length band and nu in a compact subinterval of (0,1), the zero-temperature chemical potential mu_0(nu) is unique, is away from the band edges, and has nonzero compressibility: differentiating the filled interval gives a positive sum of inverse Fermi speeds divided by 2pi. The finite-temperature occupation equation from the preceding minimization has a unique solution mu_T. Its monotonicity and convergence to the zero-temperature filling at any two chemical potentials bracketing mu_0 force mu_T->mu_0, uniformly on a smaller compact density interval. The simple-root proof of (25) is uniform in the corresponding compact chemical-potential interval because its patches, Jacobian bounds and off-patch gap can be chosen uniformly. Applying the elementary maximum version of lemma (30) below to f_grand(T,mu)+mu nu therefore proves the same T^2 coefficient at fixed density. In particular the canonical length contribution is -pi T^2/[6(2J_len sin(pi nu))]. No finite-size level-spacing comparison is used.

[assumed, imported exact quotient spectrum] For the common-bond reference, C44 (24) gives, at L=bM+r and 0<=r<=M,

    E=kappa M+E_a-2J_len sum_(j=1)^r cos(q_j+2pi ell/(LM)),
    q_j=2pi[I_j+(r-1)/2]/M,
    sum_j I_j+r(r-1)/2+m_a+ell=0 mod M,
    ell=0,...,L-1.                                        (29a)

Here a labels a simultaneous eigenbasis of the common payload Hamiltonian and cell translation, m_a is its translation label, and the I_j are distinct length-fermion momentum labels. This formula includes the hard-core boundary twist. It is not the finite-wrapper Hamiltonian with its length-dependent error terms.

[proved, theorem Q, quotient pressure at fixed M and L] Let Z_prod(M,r) be the product payload/length partition function with the same kappa and q_j but with the extra flux in (29a) set to zero. Then

    b e^(-4pi beta |J_len|) Z_prod(M,r)
       <= Z_phys,ref(L,M)
       <= (b+1)e^(4pi beta |J_len|) Z_prod(M,r).             (29b)

For fixed a and {I_j}, the congruence in (29a) has between b and b+1 representatives ell in {0,...,L-1}, since L=bM+r. For any such ell, the additional angle has magnitude less than 2pi/M. The cosine is 1-Lipschitz, so its total energy shift is bounded by 2|J_len|r(2pi/M)<=4pi|J_len|. Bound each Boltzmann factor above and below and sum over the exact complete spectrum to obtain (29b). This proof does not require a quadratic payload, only the common-bond spectrum (29a).

[proved, thermodynamic consequence of theorem Q] Taking logarithms of (29b) and dividing by beta M proves that the quotient and product reference have the same thermodynamic free energy per cell at any fixed beta and density limit r/M->nu, whenever the product limit exists. The error is bounded by 4pi|J_len|/M plus log(b+1)/(beta M). Dividing instead by L gives the same conclusion per physical site. This closes the cyclic-quotient part of the common-bond reference pressure transfer. It does not remove the finite-wrapper length-dependent interactions or prove their Gaussian stability. The thermodynamic limit is taken before the low-temperature limit.

[proved, physical length conversion, conditional on a common reference ensemble] At fixed length density nu, physical circumference is L=ellbar M with ellbar=b+nu. A wave of cell wavenumber k has physical wavenumber k/ellbar, so its physical speed is ellbar times its cell speed. Free energy per physical site is free energy per cell divided by ellbar. For the common-bond reference, theorem Q justifies using the product pressure in this conversion. Equations (27)-(28) are therefore covariant:

    f_phys(beta)=e_phys-(pi/12)sum_i |u_i,phys|^(-1)beta^(-2)
                          +o(beta^(-2)),
    u_i,phys=ellbar u_i,cell.                              (29)

This is a coordinate calculation for a specified reference state, not a proof of equivalence of the full microscopic ensembles.

[proved, equilibrium-density lemma] A useful sufficient ensemble transfer is the following. On a compact density interval suppose

    g_T(nu)=g_0(nu)-A(nu)T^2+r_T(nu),
    sup_nu |r_T(nu)|/T^2 ->0,

where g_0 and A are continuous, g_0 has a unique minimum nu_0, and every minimizing nu_T stays in this interval. Then

    min_nu g_T(nu)=g_0(nu_0)-A(nu_0)T^2+o(T^2).            (30)

Indeed uniform convergence implies nu_T->nu_0: otherwise compactness yields a subsequence a fixed distance away and a strictly positive ground-energy difference, contradicting minimality. Comparing g_T(nu_T) with g_T(nu_0) gives

    0<=g_0(nu_T)-g_0(nu_0)
      <=T^2[A(nu_T)-A(nu_0)]+o(T^2)=o(T^2).

Substitute back to prove (30). This result needs no asserted differentiability of an interacting ground energy. It does need the uniform thermal expansion and the minimizing-density control. C44's free density curvature supplies the ground minimum for its common-bond reference, not these assumptions for the full microscopic chain.

[proved, physical Gaussian reference with variable marker count, conditional on the imported reference spectrum] The preceding results determine the full selected Gaussian reference as well. Write its zero-temperature physical energy as

    g_0(nu)=[kappa+e_pay-(2J_len/pi)sin(pi nu)]/(b+nu).

Choose kappa so that a chosen nu_0 in (1/4,3/4) is stationary, as in C44. The numerator of g_0' is (b+nu)e_len'(nu)-e_len(nu)-kappa-e_pay; its derivative is (b+nu)e_len''(nu)>0 on (0,1). Therefore the stationary point is the unique global minimum on [0,1]. At finite T, the free canonical and grand-canonical formulas converge uniformly to their ground energies, since each fermionic occupation contributes at most T log 2 of entropy. Thus minimizing densities eventually lie in a compact interior neighborhood of nu_0. On that neighborhood, the uniform expansion just proved and the payload expansion (27) give the hypotheses of (30). The equilibrium coefficient is consequently the cell coefficient divided by b+nu_0, as in (29).

[proved, summing the reference marker sectors] At fixed L, the permitted M lie in [L/(b+1),L/b] and number at most L+1. The sum of their positive partition functions lies between its largest term and (L+1) times that term. Hence its free energy per physical site differs from the smallest sector free energy by at most T log(L+1)/L. The allowed density mesh becomes dense as L grows. The uniform canonical bin argument, the payload Riemann sums and (29b) give convergence to min_nu g_T(nu). This proves the thermodynamic and low-temperature reference coefficient with the physical quotient and variable marker count, not just the rooted grand-canonical formula. It still concerns the exactly common, quadratic-payload reference, not the finite-coupling wrapped Hamiltonian or its full non-Gaussian second-order block.

[proved, conditional passage of a coefficient] If an interacting model and a Gaussian comparison satisfy

    f_H(beta)-e_H-[f_G(beta)-e_G]=o(beta^(-2)),              (31)

then their beta^(-2) coefficients coincide, by subtracting their expansions and multiplying by beta^2. Theorem T establishes (31) between D and H_P whenever H_P already has a thermal expansion. It does not establish (31) between H_P and the Gaussian model, or between H and D. A constructive pressure theorem for the resummed action is a sufficient way to prove these two remaining comparisons.

[proved, why normalized correlations do not by themselves prove (31)] For any generating functional Z_beta[J], replacing it by exp[M a(beta)] Z_beta[J], with a independent of J, leaves every normalized source derivative unchanged and shifts the free energy by -a(beta)/beta. Taking a(beta)=kappa/beta changes the beta^(-2) coefficient. This is an algebraic observation about information lost when a proof normalizes away the vacuum functional. It is not a claim that an arbitrary temperature-dependent multiplier is generated by a temperature-independent microscopic Hamiltonian. It shows why a correlation theorem requires additional control of the scalar vacuum terms to imply (31).

## 8. Central charge, a conditional current bound, and the unresolved theorem

[proved, Gaussian central charge] One chiral compact-boson oscillator algebra has modes J_m with [J_m,J_n]=m delta_(m+n,0) and J_m annihilating the vacuum for m>=0. Its stress generators are L_m=(1/2)sum_l :J_(m-l)J_l:. Commuting with one current gives [L_m,J_n]=-n J_(m+n), by using the current commutator in each of the two factors. Therefore [L_m,L_n]-(m-n)L_(m+n) commutes with every current and is scalar on the oscillator representation; it vanishes unless m+n=0 by the L_0 grading. For m>0 the norm calculation

    L_(-m)|0>=(1/2)sum_(l=1)^(m-1)J_(-l)J_(-(m-l))|0>,
    ||L_(-m)|0>||^2=(1/2)sum_(l=1)^(m-1)l(m-l)
                      =(m^3-m)/12

fixes that scalar. The finite sum follows from the formulas for sum l and sum l^2; the two contractions give the factor 1/2, including the coincident case. Thus the Virasoro central charge is one. Independent chiral oscillators add their central terms. N right and N left Gaussian branches have c_R=c_L=N, and the convention assigning c=1 to a nonchiral boson gives c=N. Velocities multiply the Hamiltonian, not the current commutators. If the microscopic chain is proved to have exactly the N=n+1 Gaussian modes and no additional massless sector, this establishes c=n+1 in the conjecture's multi-velocity sense.

[proved, conditional lower bound from currents] A weaker route assumes a unitary CFT limit with r independent dimension-one chiral U(1) currents, positive level matrix, and the usual vacuum representation. Change basis to normalize that matrix to the identity. Let L^J be their Sugawara stress from the preceding calculation, so its central charge is r. Assume the full stress acts on each current as [L_m,J_n^a]=-n J_(m+n)^a, the dimension-one primary relation. Then L_m-L_m^J commutes with all currents and hence with L_n^J. Subtracting the two Virasoro algebras yields a commuting stress with central charge c-r. Its vacuum norm is

    ||(L_(-2)-L_(-2)^J)|0>||^2=(c-r)/2 >=0.

To justify the equality, use the resulting commutator at m=2,n=-2, the annihilation of the vacuum by positive modes and L_0|0>=0. Unitarity proves c>=r. Apply this separately to both chiralities when both sets of currents are established. This is a conditional CFT theorem, not an identification of lattice currents by a free-root count.

[proved, limitation of the current route] The exact microscopic payload U(1) supplies one conserved charge, not n independent chiral current algebras. The rooted length number supplies another kinematic conservation law, not the remaining emergent chiral charges. To use the preceding bound with r=n+1 one must construct all those chiral currents in the scaling limit and prove their positive, nondegenerate level matrix and stress transformation. The finite nonresonance certificate alone does not do so. A lower bound also would not exclude a further neutral massless sector, so it would not establish equality without additional information.

[assumed, outstanding microscopic statement D*] For the fixed C42/C44 schedule, there exists a symmetry-preserving, source-matched, low-frequency effective action after eliminating the massive complement, with a well-defined summed kernel norm and two time/space moments, whose error relative to the complete known second-order retained Hamiltonian is o(W) as lambda->0. The constants are independent of M, have the required low-temperature control, and retain a positive dressed cost of the excluded physical sectors. Thermal windings and the total-length/quotient constraints are included. Section 3 shows why an unresummed expansion is not a proof of D*.

[assumed, outstanding RG statement R*] The action in D* admits a constructive multichannel RG at fixed bare parameters, with dressed momenta remaining within delta_all/8, positive speeds and a controlled Gaussian stiffness. Its arbitrary-degree retarded interactions and all source insertions are allowed, and the resulting infrared limit has n+1 Gaussian modes. The full primary assumptions of the proposed cited theorem, or a proved extension to this class, are needed. C45's standard input S includes this step; it is not proved by its coupling schedule.

[assumed, outstanding thermodynamic statement T*] The microscopic scalar vacuum functional and the physical ensemble transfer satisfy (31) relative to the dressed Gaussian fixed point, or a scaling-limit theorem identifies the complete current/stress algebra and excludes additional massless degrees of freedom. These are alternative sufficient routes to a central-charge conclusion. They are not consequences of the local norm bound on R.

[proved, conditional closure] D*, R* and the pressure version of T*, together with the imported finite construction, imply (28)-(29) with n+1 right/left pairs and c=n+1 if the Gaussian scaling-limit identification in R* is included. The proof is the coefficient transfer (31), the oscillator calculation (28), and the current-algebra computation in this section. Theorem T removes the massive thermal contribution of the reference D from this list of assumptions. It leaves leakage and the interacting logical pressure as separate issues. A beta^(-2) coefficient alone, without measured speeds or a mode/stress identification, is velocity weighted and does not by itself determine an unweighted central charge.

[heuristic, shortest path to closure] First target D* at low temperature with polynomial time moments, rather than the stronger raw all-beta exponential action theorem. Use operator-valued excursions to preserve identities such as (16), subtract retained propagation, and prove a connected, resummed polymer norm that includes source insertions. Theorem T can control the reference thermal defect trace, but does not yet control winding terms after leakage is restored. Second check and extend the multichannel RG initialization to that norm, retaining the exact U(1)'s and (24). Third keep the scalar vacuum terms through the RG and prove their thermal remainder is o(T^2), then transfer the ensemble with a uniform version of (30). If only a violation under the conjecture's CFT premise is sought, constructing n+1 positive-level chiral currents could replace the equality calculation, but it still requires a rigorous scaling limit and matching of those currents.

## 9. Checks, files and theorem boundary

[numerical] checks.py verifies (7) and (5) on two, three and four cells of local dimension three and code dimension two, using local Q-block interactions that move and change the defect count. The commutator with total defect number has norm 0.0321094939. All nine partition-ratio inequalities pass. The smallest operator slack is -2.40e-16, consistent with roundoff in the code block. This checks the proof's distinction between local P-block preservation and individual defect conservation.

[numerical] For Delta=1,u=2,beta=4, (12) is -0.00113393226 and the scalar Taylor radius (18) is 0.635777138. Absolute coefficient partial sums through degrees 10,20,40,80 are approximately 47.3,1952,5.34e6,2.70e14. Standard-library Decimal arithmetic uses 120 digits for this diagnostic. The divergence itself is proved by (18), not inferred from the finite list.

[numerical] A two-band Hermitian nearest-neighbor reference with nonsymmetric chiral speeds has four signed speeds approximately -2.08160529,-0.609131854,1.93582270,0.622857788. Direct quadrature of (26) gives -(f-e)/T^2=1.111135346 at T=0.0025 versus the predicted 1.111118049, with errors decreasing as T decreases. This is a control of the chiral-speed coefficient, not a test of C42's tiny lacunary parameters or of interacting central charge.

[proved, completed theorem boundary] Theorem T is a uniform thermal decoupling result for locally block-diagonal references, valid even with a gapless code. Theorem Q controls the common-bond physical cyclic quotient in the thermodynamic pressure. The raw-action counterexamples in Section 3, the second-order excursion bounds (19)-(22), the augmented finite momentum test (23)-(24), and the quadratic/continuum-Gaussian thermodynamic formulas (25)-(30) have proofs here. Their application to the supplied microscopic parameters is conditional on the explicitly imported finite construction and C45 reduction. The exact charge proof uses the actual local rules. No all-order massive-sector action, microscopic RG flow, interacting thermal coefficient, or c=n+1 theorem has been proved for the full 4464-letter chain.

[proved, work record] New files are PROOF.md, RESULT.md, PROGRESS.md, STATE.md, checks.py and checks.json. The local controls use only Python's standard library, numpy and scipy and take less than a second in the recorded run. No network, delegated agents, large computation or external file writes were used. The project STATE was not modified because the user forbids writes outside this directory. STATE.md here is the local handoff.
