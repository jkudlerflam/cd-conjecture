# C62. Thermal XXZ comparison and a trace-defined reference

## 1. Verdict and scope

[unproved target] Lemma TP of C59 Section 7 is not proved here. In particular, there is no unconditional conclusion about the thermal coefficient of XXZ at fixed nonzero anisotropy. The obstruction is an unestablished uniform estimate, not a counterexample to the predicted coefficient. Sections 5 and 7 specify the estimate and distinguish it from an elementary scale summation.

[proved partial results] This note supplies a trace-defined Hamiltonian reference with a convergent finite-fermion regulator, an explicit conversion between reference parameter conventions, a finite-cutoff current-algebra calculation, and a joint-mark scale summation. The trace construction fixes the scalar of this chosen reference. It does not prove that its effective actions satisfy TP, or identify its scalar with the different Euclidean reference used in BM by an unproved universality argument.

[read, imported inputs] C57's four-history identity and bootstrap, C59's curvature and Gaussian-symbol bounds, and C51/C59's oscillator pressure evaluation are retained. Their proofs are not repeated. The collaborator's WINDING_ESTIMATE proves the fixed-data single-shell comparison; THERMAL_LOCALIZATION supplies p=3; MOMENT_RESERVE Section 6 allocates a fixed thermal mark conditional on ordinary distance routing; ensemble/THERMAL_HISTORY proves the nonlinear closure conditional on uniform derivatives. None of these is an unconditional all-depth XXZ estimate. C60's source-prefactor counterexample does not refute TP, whose defect amplitude need not vanish at zero interaction.

[convention] Use C59's Hamiltonian, with lattice spacing and k_B equal to one,

    H_X = -t sum_x(c_x^dagger c_(x+1)+h.c.)
          +V sum_x(n_x-1/2)(n_(x+1)-1/2)-mu sum_x n_x,
    Delta=V/(2t).

The Fermi point lies in a fixed compact subset of (0,pi) away from pi/2. Physical counterterms are fixed at zero temperature. The ultraviolet energy unit is one, q_h=2^(-h), and T<=q_H<2T. All constants may depend on the fixed filling and t. The thermodynamic limit precedes T tending to zero. No bound uniform near half filling or a band edge is claimed.

## 2. Two reference conventions must be distinguished

[read, primary-source input] BM [0907.2837, (29)-(33)](https://arxiv.org/pdf/0907.2837) uses a Grassmann reference with a spacetime-smoothed interaction and a frequency-momentum cutoff. Its (52) and Appendix D distinguish this reference's parameter tau=lambda_infinity/(4 pi v) from the Hamiltonian Luttinger parameters. The reference matching in Lemma 2.1 concerns infrared correlations. Its scalar is not identified with a Hamiltonian trace there. BFM [1303.3684, Section 2.1](https://arxiv.org/pdf/1303.3684) explicitly describes its related effective model as non-Hamiltonian.

[proved parameter conversion] Consider canonically normalized right and left densities, with low-momentum Hamiltonian stiffness

    H = pi integral dx [A(rho_+^2+rho_-^2)+2B rho_+rho_-],
    A>|B|.

The current Hamiltonian formulas, already established in C51, give

    u=sqrt(A^2-B^2),       K=sqrt((A-B)/(A+B)).             (1)

Given u=v>0 and K>0, solving these equations gives

    A=(v/2)(K+K^(-1)),     B=(v/2)(K^(-1)-K).             (2)

Indeed A-B=vK and A+B=v/K. Their product is v^2 and their ratio is K^2. Thus (2) is positive and stable. Using the BM convention K=(1-tau)/(1+tau), |tau|<1, gives

    A=v(1+tau^2)/(1-tau^2),
    B=2v tau/(1-tau^2),
    sigma=B/A=2tau/(1+tau^2).                            (3)

If g_H denotes the coefficient of integral rho_+rho_-, then

    g_H=2pi B=lambda_infinity/(1-tau^2).                 (4)

This is an algebraic conversion of infrared parameters, not an equality of the full regulated actions. Field and source amplitudes still need their own matching.

[proved consistency check] Setting instead A=v and g_H=lambda_infinity would give

    u_naive=v sqrt(1-4tau^2),

for |tau|<1/2. Its thermal coefficient would differ from -pi/(6v) at order tau^2. Smallness of the discrepancy in a bare parameter is not enough to make the pressure error o(T^2). Equations (2)-(4) avoid this error. At tau=0 they give A=v and B=0. Reversing tau changes the sign of B and replaces K by K^(-1), while leaving A and u unchanged.

## 3. A specific Hamiltonian reference and its fermion regulator

### 3.1 Definition and imported current representation

[read, imported representation] The filled-sea fermion Fock space, its current algebra, its complete charged oscillator basis, and the unitary diagonalization of a stable spatially smoothed Luttinger Hamiltonian are standard operator constructions. A primary account is [Langmann, Lebowitz, Mastropietro and Moosavi, 1511.01884, Section 3 and Appendix B](https://arxiv.org/pdf/1511.01884). This representation is used here as an imported theorem. It is not deduced from lattice density correlations. The following choice has only finitely many interacting oscillator pairs at each fixed L, which makes the regulator argument particularly direct.

[definition] Choose an even real C-infinity function chi(k), equal to one near k=0, between zero and one, and zero for |k|>=kappa. Keep A from (2) constant at all momenta and put B(k)=B chi(k). In a spatial circle of length L, take antiperiodic fermion momenta and transfer momenta k_m=2pi m/L. Relative to the filled sea define

    H_ref = H_0(A) + (2pi/L) sum_(m in Z)
                       B(k_m) rho_+(k_m)rho_-(-k_m).     (5)

Both densities are normal ordered relative to that sea. At m=0 their values are the integer charges Q_+,Q_-. The two chiral algebras commute. Thus there is no same-chirality density square in the interaction in (5). The spatial potential is Schwartz; no exponential spatial bound is asserted for this compact Fourier cutoff.

[proved using the imported current representation] With canonical annihilators a_m,b_m for m>0, (5) is

    H_ref = H_Q + sum_(m>0) k_m [A(N_a,m+N_b,m)
                                +B(k_m)(a_m^dagger b_m^dagger+a_m b_m)],
    H_Q=(pi/L)[A(Q_+^2+Q_-^2)+2B Q_+Q_-].              (6)

The normalization follows from rho_+(k_m)=sqrt(m) a_m and rho_-(-k_m)=sqrt(m) b_m. Combining the m and -m terms in (5) gives the pair creation and annihilation terms in (6). The free current identity supplies H_0(A).

For each interacting pair a finite two-mode Bogoliubov transformation gives energy u(k_m)k_m per excitation, where

    u(k)=sqrt(A^2-B(k)^2)>=sqrt(A^2-B^2)=v.

Only finitely many pairs have B(k_m)!=0 at fixed L. Their product transformation is therefore unitary on the Fock space. The total vacuum shift is the finite sum

    E_ref(L)=sum_(m>0) k_m[u(k_m)-A].                    (7)

The charge form is positive because A>|B|. The diagonalized operator is a positive charge form plus positive oscillators and (7). Its Gibbs trace is finite for every L,beta>0: the charge sum is a convergent two-dimensional Gaussian sum and the oscillator logarithm converges since k_m grows linearly with m. This proves that the scalar trace in (5) is well defined.

[read, imported thermal evaluation] C59 Section 6 now applies to this reference. Since u(k)->v at zero, the thermodynamic-first trace obeys

    f_ref(T)-e_ref=-pi T^2/(6v)+o(T^2).                 (8)

The evaluation and the subextensive charge contribution are already proved in C51/C59 and are not rederived. Equation (8) concerns (5), with its full trace normalization.

### 3.2 Finite fermion windows converge in the Gibbs trace

[proved regulator lemma] Let P_N keep every fermionic occupation inside a symmetric window of 2N momenta per chirality arbitrary, fixing the occupations outside the window to their filled-sea values. Take nested windows exhausting all momenta. The finite-dimensional space F_N=Ran P_N has dimension 2^(4N). Define H_N as the quadratic-form compression of H_ref to F_N. Then for each fixed L,beta>0,

    Tr_(F_N) exp(-beta H_N) increases to
                                      Tr_F exp(-beta H_ref),
    min spec H_N decreases to E_ref(L).                (9)

Here and below the word increases permits equality.

Proof. First the finite-excitation space is a form core. For a two-mode vector psi,

    |<psi,(a^dagger b^dagger+ab)psi>|
                 <=<psi,(N_a+N_b+1)psi>.              (10)

This follows by writing the expectation of a^dagger b^dagger as the inner product of a psi and b^dagger psi, then applying 2|<x,y>|<=||x||^2+||y||^2. The charge cross term satisfies 2|Q_+Q_-|<=Q_+^2+Q_-^2. Consequently, with b=|B|/A<1,

    (1-b)H_0(A)-C_L <= H_ref <= (1+b)H_0(A)+C_L,
    C_L=sum_(m>0) k_m |B(k_m)|<infinity,              (11)

as forms. The interaction has relative form bound less than one, so the closed form domains agree and their shifted form norms are equivalent. Finite free oscillator occupation vectors in finitely many charge sectors form a core for the diagonal H_0. Each such vector is a finite product of current modes acting on a charge-shifted sea. Each current mode on a finite fermion excitation has only finitely many nonzero summands: deep filled pairs are Pauli blocked and high empty pairs annihilate the vector. Hence every such oscillator vector lies in some F_N. This proves that the union of F_N is a form core for H_ref.

Let E_j be its ordered eigenvalues, with multiplicity, starting at j=0, and let E_j^(N) be the ordered Ritz eigenvalues in F_N. Compact resolvent follows from the diagonalization above. The min-max principle gives E_j^(N)>=E_j and monotonic decrease as the trial space grows. Approximate the first j+1 exact eigenvectors in the form norm by vectors in a common F_N. Their Gram and form matrices converge to those of the exact vectors, so the maximal Rayleigh quotient on their span tends to E_j. The opposite min-max inequality then gives E_j^(N)->E_j for each j. This proves the second assertion in (9).

For the trace, extend exp(-beta E_j^(N)) by zero when j>=dim F_N. These terms increase pointwise to exp(-beta E_j). Monotone convergence for the counting measure gives the first assertion in (9). Finiteness was proved above. This proof uses spectral min-max, not the false inference that exponentiation preserves arbitrary operator inequalities. QED.

[proved cutoff identification] For the cross-chirality interaction (5), H_N is a finite-mode fermion Hamiltonian of degree at most four. In fact P_N=P_(N,+) tensor P_(N,-), and

    P_N rho_+(k)rho_-(-k)P_N
      = [P_(N,+)rho_+(k)P_(N,+)] tensor
                            [P_(N,-)rho_-(-k)P_(N,-)]. (12)

Thus it is the product of the two truncated normal-ordered currents. The free Hamiltonian is diagonal in occupations and compresses directly. Zero-mode normal ordering produces the corresponding finite constants and quadratic terms when written in un-ordered occupations. Those terms are included, not discarded. The factorization (12) would not be valid for an arbitrary product of two currents in the same chirality. That modification is unnecessary for our reference.

### 3.3 Exact Grassmann normalization at finite regulator

[proved finite-regulator identity] Write H_N=H_(0,N)+V_N and Z_(0,N)=Tr exp(-beta H_(0,N)). The Duhamel identity gives

    Z_N/Z_(0,N)
      = sum_(r>=0) (-1)^r integral_(0<t_r<...<t_1<beta)
          <V_N(t_1)...V_N(t_r)>_(0,N) dt_1...dt_r.     (13)

To justify it, solve the interaction-picture differential equation for exp(-beta H_N) by iteration and then take the finite-dimensional trace. At fixed N,L,beta, all operators are bounded and the ordered series converges, bounded for example by exp[beta sup_(0<=s<=beta)||V_N(s)||]. The free thermal expectation is quasifree, so every coefficient equals its time-ordered Grassmann Wick expression, including its equal-time prescription. The right side is consequently the normalized fermionic functional integral, with Z_(0,N) as its Gaussian partition factor. Equivalently, a finite coherent-state time slicing converges to (13) by the finite-dimensional Trotter formula.

Any exact covariance splitting and Grassmann Gaussian integration preserve this full expression. At finite regulators the connected logarithm is its analytic germ at zero interaction; every covariance-change determinant and constant term must be retained. The total real logarithm is fixed by the positive Hamiltonian trace. This argument supplies no cutoff-independent disk of convergence for the individual connected RG series.

Combining (9) and (13), with the time regulator removed first, proves

    lim_(N->infinity) [-T/L log Z_N - E_0^(N)/L]
         = -T/L log Tr exp(-H_ref/T)-E_ref(L)/L         (14)

at every fixed L,T>0. Taking the thermodynamic limit afterwards identifies its scalar with the oscillator reference (8). No uniform-in-T trace convergence is needed for this specified order of exact limits.

[proved boundary] Equations (9)-(14) remove a scalar normalization ambiguity for a Hamiltonian-derived reference. They do not give an all-depth analytic RG theorem for its finite-window kernels. In particular, neither the sharp window's position-space bounds nor a smooth ultraviolet decomposition uniform in N have been verified here in the common norm of TP. Nor do (13)-(14) identify the different spacetime-smoothed BM reference with (5). One must either prove TP using this Hamiltonian-derived reference, or compare the other reference's ground-subtracted scalar with it to o(T^2).

## 4. Why finite-cutoff bosonization cannot be substituted for (9)

[proved finite-CAR identity] In one chiral window label the modes j=-N,...,N-1, and for 1<=m<=N set

    J_m^(N)=sum_(j=-N)^(N-1-m) c_j^dagger c_(j+m),
    J_(-m)^(N)=(J_m^(N))^dagger.

Then

    [J_m^(N),J_(-m)^(N)]
      =sum_(j=-N)^(-N+m-1) n_j
                      -sum_(j=N-m)^(N-1) n_j.         (15)

Proof. The CAR give

    [c_i^dagger c_j,c_k^dagger c_l]
                   =delta_(j,k)c_i^dagger c_l
                                  -delta_(l,i)c_k^dagger c_j.

In the two current sums both Kronecker conditions set their summation indices equal. The result is sum_(j=-N)^(N-1-m)(n_j-n_(j+m)). Cancelling overlapping ranges leaves (15).

On the filled sea (15) has eigenvalue m. On the full finite Fock space it is not m times the identity. Indeed

    m I-[J_m^(N),J_(-m)^(N)]
       =sum_(j=-N)^(-N+m-1)(1-n_j)
                         +sum_(j=N-m)^(N-1)n_j.       (16)

Its operator norm is 2m for every N>=m. The normalized finite-dimensional trace of (15) is zero, as it must be for a commutator. Thus the anomalous infinite-sea current relation has no operator-norm approximation of this type. On every vector with excitations in a fixed inner window the right side of (16) does vanish for N sufficiently large. This is convergence on a core, not norm convergence.

[proved free check] For free energies epsilon_j=(2pi A/L)(j+1/2), occupations factor in the finite Gibbs state. Every edge hole or particle in (16) has cost at least (2pi A/L)(N-m+1/2), so

    <m I-[J_m^(N),J_(-m)^(N)]>_0
        <=2m exp[-(2pi A/(LT))(N-m+1/2)].              (17)

This follows from (1+exp(x))^(-1)<=exp(-x) for x>=0. It is not an interacting estimate. Equation (9) uses a different argument and does establish trace convergence for the interacting reference.

[proved implication for the attempted shortcut] An exact finite-dimensional current algebra with central commutator m I is impossible. Therefore replacing truncated fermion currents by canonical oscillators inside each regulated Grassmann integral is invalid. Removing the sea cutoff with a justified trace argument, as in Section 3, is a legitimate alternative. This is an obstruction to that shortcut, not to TP.

## 5. There is no conflict between the two scale powers

[proved summation lemma] Suppose a completed renormalized-tree majorant with a defect attached at scale k<=h and a thermal mark at scale l<=h has the weight

    2^(-sigma max(h-k,h-l)) q_k^alpha (T/q_l)^p,
    sigma>0, alpha>0, p>0.                              (18)

There is only one short-memory factor here. For any

    0<theta<min(alpha,sigma),

the sum over both attachment scales obeys

    sum_(k,l=0)^h (18)
       <= q_h^theta (T/q_h)^p
            /[(1-2^(-(sigma-theta)))(1-2^(-p))].        (19)

Proof. Since q_k<=1, replace q_k^alpha by q_k^theta. Set r=h-k, s=h-l. The summand divided by q_h^theta(T/q_h)^p is

    2^(theta r-p s-sigma max(r,s))
                      <=2^(-(sigma-theta)r)2^(-p s).

Extend both finite sums to infinity and evaluate their geometric series. QED. Polynomial factors in r,s are also summable with strict inequalities, by differentiating the convergent geometric series a fixed number of times.

[proved conditional consequence] If a GN construction supplies (18), a positive short-memory reserve and the collaborator's fixed-mark moment allocation, then p=3 and any theta<min(alpha,sigma,1) satisfy p>2+theta. The thermal mark does not require a second independent copy of the short-memory factor on a shared ancestry path. Taking theta=sigma can instead produce an h factor in the r sum, so a strict reserve matters.

[unproved application] Equation (18) is not established here for the actual joint XXZ/reference expansion in TP's complete kernel norm, with scalar roots and its required derivatives. Equation (19) only performs the remaining scale sum if that majorant has first been obtained. The collaborator's MOMENT_RESERVE already proves how a fixed extra thermal moment is paid after ordinary routing; it does not establish (18) for this comparison.

## 6. A global Cauchy argument needs matched critical interpolation

[proved obstruction to an inference] Uniform analyticity along two matched flows does not imply uniform analyticity on their affine interpolation in bare coordinates. This matters if one tries to promote the single-shell Cauchy argument to an all-depth argument without controlling the critical manifold.

An exact example is the analytic two-coordinate map

    F(nu,r)=(2nu+r^2,r/2).                              (20)

Its stable curve is nu=-(4/7)r^2: substitution gives

    2[-(4/7)r^2]+r^2=-(1/7)r^2
                                  =-(4/7)(r/2)^2.

The points x=(-4a^2/7,a) and y=(0,0) remain near zero under every iterate. Their affine interpolation at 0<t<1 is x_t=(-4ta^2/7,ta). Its displacement from the stable curve at r=ta is

    e_0=-(4/7)a^2 t(1-t).

For e=nu+(4/7)r^2, direct substitution into (20) gives e'=2e. Therefore

    r_h=ta 2^(-h),
    nu_h=-(4/7)t^2 a^2 2^(-2h)
                              -(4/7)a^2t(1-t)2^h.     (21)

Every fixed interior interpolation leaves any fixed small domain at sufficient depth, although both endpoint trajectories stay inside it.

[proved scope] This example is not an XXZ flow or a counterexample to the local four-history bootstrap. Interpolations in a single rescaled state domain can still work. It shows why a bare interpolation used for Cauchy differentiation must either remain on the matched critical manifold, or be controlled in appropriate scale-dependent neighborhoods. The zero-temperature existence of the two endpoints does not supply that control. TP2 and TP3 explicitly ask for the requisite common domains and derivatives; this note has not constructed them.

## 7. The remaining obstruction for XXZ

[unproved estimate, exact stopping point] Let Y be the complete effective action of the trace-defined model (5), including all generated kernels, not just A,B and a local quartic. Let X be the XXZ action in the same slow-field coordinates, with the zero-temperature counterterms and field amplitude matched. What is still needed is a normalized GN comparison, uniform in infrared depth and the removed regulators, giving in C59's norm

    ||z_h||=||X_h^0-Y_h^0||<=C q_h^theta,
    ||a_h||=||Y_h^T-Y_h^0||<=C(T/q_h)^p,
    ||d_h||=||(X_h^T-Y_h^T)-(X_h^0-Y_h^0)||
                              <=C q_h^theta(T/q_h)^p,  (22)

and, for the complete vacuum map,

    ||DE_(h,T)||+||D^2E_(h,T)||<=C q_h^2,
    ||D(E_(h,T)-E_(h,0))||<=C q_h^2(T/q_h)^p,          (23)

with p>2+theta, in compatible domains. Equivalent weighted history estimates are allowed. For example, the collaborator's THERMAL_HISTORY requires weighted sums of history Hessians and thermal derivatives such as

    sum_(j,k>=0) M_(j,k) 2^(theta j-p k)<infinity,
    sum_(j>=0) G_j 2^(theta j)<infinity,               (24)

uniformly at every root scale. Bounded localized couplings do not prove (24).

At the tree level, the missing operation is an absolute summation with both (i) a curvature, matched-history or full-kernel reference difference and (ii) a winding, thermal-localization or thermal-history replacement. Ordinary Taylor distances, derivative half-lines, determinant cofactors, and normalized Z ratios must be bounded in that same summation. The output includes vacuum roots. One fixed reserve of field weight and moments must suffice at all depths. Proving an individual-scale estimate while losing that reserve at every scale is insufficient.

The reference construction above does not solve this problem. It adds the requirement that the ultraviolet preparation of the finite-fermion reference (or an equivalent trace-compatible regulator) admit the same uniform marked bounds. The scalar convergence (9) has no rate or kernel norm and cannot supply them. Conversely, using BM's spacetime-smoothed reference without changing its construction leaves an additional ground-subtracted scalar comparison to (5) unproved.

[read, available ingredients and boundary] The ordinary BM theory proves weak-coupling flow and correlation comparisons for XXZ. This note does not assert that its ordinary instantaneous-model tree estimates are missing. The unverified extension is their simultaneous finite-temperature, reference-subtracted, differentiated vacuum version, with the chosen common reference regulator. The supplied retarded-endpoint routing problem is broader and should not be confused with failure of the ordinary XXZ construction.

[proved conditional conclusion, imported from C59] If (22)-(23), their ultraviolet and terminal counterparts, and an exact scalar expansion are established for the reference (5), C57/C59 imply

    [f_X(T)-e_X]-[f_ref(T)-e_ref]=O(T^(2+theta)).        (25)

Section 3 fixes the scalar on the reference side, so (8) then gives the desired thermal c=1. Equation (25) is not proved for XXZ in this note. No temperature-dependent retuning, source-normalized partition function, or independently assumed lattice CFT is used to replace it.

## 8. Checks and work boundary

[proved checks] The Hamiltonian parameter map preserves positive stiffness and has the correct free limit. Its finite-L trace is finite, its fermion projections are nested, and min-max fixes the direction of both limits in (9). The central current term is recovered on every fixed finite-excitation vector but not in operator norm, consistent with finite-dimensional trace cyclicity. The thermal weight in (19) has the same exponent p as at its attachment scale, with no conflict with p>2+theta. C57's settled first interaction coefficient is not recomputed.

[numerical/exact finite controls] checks.py tests the current commutator using independently assembled finite CAR matrices, the parameter conversion with rational arithmetic, the two-mark scale sum, and the unstable interpolation example. It also constructs the two-chirality reference Hamiltonian at N=1,2, with only transfer modes 0,+1,-1 interacting, and compares its ordered eigenvalues and two Gibbs traces with the full charged oscillator trace. The min-max inequalities, increasing traces and decreasing ground energies pass; the largest matrix is 256 by 256. These controls do not compute an interacting XXZ flow or test (22)-(23). Their output is checks.json.

[work boundary] All changes are confined to work_C62. No subagents, messages to authors, shell network commands, or large computations were used. The local STATE.md is the handoff under the explicit folder restriction. Neither TP nor a rigorous interacting XXZ thermal coefficient has been established here. The new arguments have not had independent or human review.
