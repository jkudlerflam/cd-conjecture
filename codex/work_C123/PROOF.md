# C123. Bounded reconnection cost and temporal transfer operators

2026-10-29. Natural logarithms. Every status tag applies to the paragraph and displayed formulas that follow it.

## 1. Definitions and theorem obtained

[Assumed conventions] Use the background Hamiltonians H_L, physical on-site dimension d, and optimized bond half-spread M. Set Z_L(beta)=Tr exp(-beta H_L) and R_L=2 log Z_L-log Z_(2L). For N>=1 let tau=beta/N, F=exp(-tau H_L/4), E=F^2, and

    W_s = F P_(s_N) E ... E P_(s_1) F,
    q_s = Tr(W_s W_s^dagger),
    S_N = sum_s q_s^2.

On the 2L-ring use the same tau,N and pinch both sites 0,L. Its channel trace is Q_D. The quantity under investigation is G=log(S_N/Q_D). All histories and both Jordan-Wigner sectors refer to these definitions. In particular, a coherent Euclidean history is not substituted for q_s.

[Proved, inherited] C121 establishes

    d^(-N) Z_L <= Q_L <= Z_L,
    d^(-2N) Z_(2L) <= Q_D <= Z_(2L),
    Q_L^2/d^N <= S_N <= Q_L^2.

Consequently,

    R_L-3N log d <= G <= R_L+2N log d.                  (1)

These are inputs. For beta=0, direct counting gives S_N=d^(2L-1), Q_D=d^(2L), and G=-log d for any N. The repeated records must then be constant. For M=0 the scalar ring evolution cancels from this ratio, giving the same value.

[Proved theorem, XX] For H_L=sum(Sx_i Sx_(i+1)+Sy_i Sy_(i+1)), d=2, M=1/2, even L>=4, every fixed a>0 and

    N=max(1,ceil(beta M/(aL))),

one has

    G <= 2 log 2 + [2 log 2+(10pi/3)a] N.               (2)

This is a bound uniform in L,beta,N, with explicit constants. The proof is short. C118 proves, including both physical parity sectors,

    R_L <= 2 log 2+(5pi/3) beta/L.

Insert it in the upper side of (1), then use beta/L<=2aN. This proves (2). The Gaussian formulas below compute the actual G and give its scaling limit, but are not needed to infer (2). This proof uses the already settled XX partition-ratio bound rather than claiming a new proof of its infrared content.

[Proved extension, inherited subclass] For independent XX stacks, d=2^m, C118 gives R_L<=2 log d+(10pi/3)M beta/L, so (2) holds with log 2 replaced by log d. C119's two-band subclass instead gives the same slope and the constant 3 log d, which carries over unchanged as b_0. No all-temperature ratio theorem for every paired or multi-band quadratic Hamiltonian is inferred from the separate known Casimir theorem A<=2M.

[Proved conditional consequence] Any bound G<=b_0+b_1 N implies

    A <= [4/(pi a)] [3 log d+b_1] M,

under the inherited vacuum Casimir hypothesis. For (2) this gives

    A <= [40/3+20 log 2/(pi a)] M.

It is weaker than the settled free-fermion bound A<=2M. The purpose of (2) is to establish the requested reconnection lemma on XX, not to improve that bound.

## 2. There is no finite maximal spacing under this criterion

[Proved] Suppose the general lemma holds at one a_0>0 with constants b_0,b_1. Write c=3 log d+b_1 and N_a=max(1,ceil(beta M/(aL))). Equation (1) gives

    R_L <= b_0+c N_(a_0),
    N_(a_0) <= 1+(a/a_0)N_a,
    G_(N_a) <= R_L+2N_a log d.

Combining these three inequalities proves

    G_(N_a) <= b_0+c+[(a/a_0)c+2 log d]N_a.             (3)

Thus existence of finite constants at one positive spacing implies existence at every finite positive spacing. Here a is fixed before defining the constants, so dependence on that chosen numerical value is permitted. The set of admissible spacings is either empty or all positive spacings. A finite threshold can arise only after imposing an additional restriction, such as a prescribed cost b_1. No such restriction is in the target lemma.

[Proved distinction] Constants uniform over all a simultaneously are impossible even on XX. The large-spacing result in Section 5 gives G/N asymptotic to pi a/2-log 2. The meaningful optimization is of (3 log d+b_1(a))/a, not the largest a with some finite b_1(a).

## 3. Exact XX multi-pinch determinants

### 3.1 The contour and parity sectors

[Proved] Cyclicity in Tr(W_s W_s^dagger) leaves 2N projectors separated by exp(-tau H/2). Their labels, in one valid cyclic order, are

    s_1,...,s_N,s_N,...,s_1.

Their times are t_j=j beta/(2N), j=0,...,2N-1. For Q_D there are two projectors, at sites 0,L, at every such time. Each fixed doubled-ring history therefore has 4N occupation insertions. Histories have 2^N or 2^(2N) independent labels, respectively.

[Proved, inherited convention] On an even n-site XX ring, the two one-particle grids are

    k_j=2pi(j+s)/n,       epsilon_j=cos k_j,
    (s,p)=(1/2,0) or (0,1),

where p is the required fermion parity. The spin trace of any number-conserving insertion string O is

    Tr_spin O = (1/2) sum_(s,p)
                  [Tr_Fock O_s+(-1)^p Tr_Fock((-1)^F O_s)].  (4)

The occupation projectors commute with total parity, so the same projection applies separately to every history. There is no selection of only the lower-energy grid.

### 3.2 Occupation determinants

[Proved] List the k insertions by their time t_i, position r_i and label b_i in {0,1}. Put e_i=1-b_i. On one nonsingular Gaussian trace, let f_j be its mode occupation contraction, and define

    K_ij = (1/n) sum_m exp[ik_m(r_j-r_i)]
             * exp[(t_i-t_j)epsilon_m]
             * { f_m,      t_i>=t_j,
                 f_m-1,    t_i<t_j } .                    (5)

For ordinary thermal traces f_m=1/(1+exp(beta epsilon_m)). For parity-inserted traces without zero modes, f_m=-1/(exp(beta epsilon_m)-1). For distinct sites at the same time either order in (5) agrees, since their anticommutator vanishes. Coincident insertions on the same site do not occur at beta>0 in the present contour.

[Proved] Wick expansion of the ordered product of k densities gives det K. To see the sign, pair each creation operator with one annihilation operator. A pairing is indexed by a permutation of the k annihilators; its fermion reordering sign is the permutation sign, and its contractions are the matrix entries (5). This is the determinant expansion. Replacing the factor n_i by 1-n_i replaces its row's diagonal entry by K_ii-1 and introduces an overall minus sign. Multilinearity therefore gives

    <ordered product_i P_(b_i)(r_i,t_i)>
       = (-1)^(sum_i e_i) det[K-diag(e_i)].                 (6)

This can also be derived by introducing factors 1+z_i n_i, applying the fermionic trace determinant, and taking coefficients. The trace determinant is the standard identity Tr Gamma(A)=det(1+A), proved by diagonalization and extension by continuity. A primary reference is [Klich, Full Counting Statistics](https://arxiv.org/abs/cond-mat/0209642). No literature assumption beyond this finite-dimensional identity is needed here.

[Proved] Every empty label appears twice on the reflected contour, so the sign outside (6) is positive for our histories. Multiply (6) by the corresponding unnormalized thermal trace, apply (4), and obtain each q_s. Finally,

    S_N=sum_(2^N histories) q_s^2,
    Q_D=sum_(2^(2N) histories) q_(D,s).

These are finite sums of determinants, with determinants of order 2N or 4N, rather than determinants of order L. The kernel formation costs O(LN^2) for each grid. Enumerating Q_D costs exponentially in N, which is why the large-ring survey stops at N=6.

[Proved limitation] Gaussian individual histories do not make the nonselective channel Gaussian. For example, pinch one site of the Slater state (|10>+|01>)/sqrt(2). The result has contractions <n_1>=<n_2>=1/2 and <c_1^dagger c_2>=0, but <n_1 n_2>=0, whereas a Gaussian state with those contractions would give 1/4. Likewise, summing identical labels on two replicas introduces the projector 1-n_A-n_B+2n_A n_B, which is quartic in fermion operators. A single Gaussian covariance update or an unexplained single determinant for the summed channel would therefore be incorrect. The finite sums above are exact.

### 3.3 Zero-mode parity traces without a regulator

[Proved] One grid has g=2 zero modes. Its parity-inserted partition trace vanishes, but its inserted traces need not vanish. Substituting a zero partition trace times a singular normalized contraction in (6) is not valid. The following bordered determinant gives the finite answer.

[Proved] Give only the zero modes a temporary fugacity z. Their contraction is z/(1+z). The kernel can be written

    K(z)=A+[z/(1+z)]UV,
    U_(i,alpha)=exp(-ik_alpha r_i)/sqrt(n),
    V=U^dagger.

Here A includes all nonzero-mode parity contractions and the zero-mode term -(UV)_ij when t_i<t_j, with no zero-mode term when t_i>=t_j. Thus A remains finite as z tends to -1. The empty-projector shifts are included in A by replacing A with A-diag(e).

[Proved] The zero-mode partition factor is (1+z)^g. The determinant is a polynomial of degree at most g in x=z/(1+z). Its degree-g coefficient is

    det [ A-diag(e)    U ]
        [    -V       0 ].                                  (7)

For invertible A-diag(e), the Schur complement proves this identity by comparing det(A-diag(e)+xUV) with det(I+xV(A-diag(e))^(-1)U). Polynomial continuity proves it without the invertibility restriction. Therefore

    lim_(z->-1) (1+z)^g det[K(z)-diag(e)]
       = (-1)^g times (7).                                 (8)

If g exceeds the insertion rank the bordered determinant vanishes automatically, as it should.

[Proved] Write E_s=sum_(epsilon<0) epsilon and n_-=#{epsilon<0}. Factoring out the filled sea, the ordinary trace prefactor is

    exp(-beta E_s) product_m(1+exp(-beta|epsilon_m|)),

and the nonzero-mode parity prefactor is

    exp(-beta E_s)(-1)^(n_-)
       product_(epsilon_m!=0)(1-exp(-beta|epsilon_m|)).       (9)

Combine (7)-(9) with the factor (-1)^p/2 from (4). This is the formula implemented in checks.py. There is no tiny mass, complex fugacity extrapolation, eigenvalue floor, or dropped parity sector.

[Proved, numerical implementation] To avoid bulk cancellation, every trace is multiplied by exp(beta E_n), where E_n=-csc(pi/n). The difference E_s-E_n is either zero or tan(pi/(2n)); which grid has the smaller sea switches at n modulo four. The final energy factor is evaluated as

    Delta_L=2 sin(x/2)^2/[sin x cos x],    x=pi/(2L),
    G=beta Delta_L+log(sum_s qhat_s^2)-log Qhat_D.            (10)

This identity is inherited from C120/C121. Exponential factors in (5) are evaluated using exp(-c|epsilon|), with c>=0, and expm1 for the parity denominator. No growing Boltzmann factor is formed.

## 4. Uniform XX theorem versus fixed-N asymptotics

[Proved] Fix a>0 and a finite integer N, set tau=2aL and beta=2aNL, and let L run through even lengths. The determinants imply

    q_s/Z_L = 2^(-2N)[1+O_(a,N)(L^(-2))],
    q_(D,s)/Z_(2L) = 2^(-4N)[1+O_(a,N)(L^(-2))].            (11)

The constants in this statement may depend on N. This is not used as a uniform-in-N proof of (2).

[Proved derivation] Distinct measurement times differ by at least tau/2=aL, and their cyclic complements also differ by at least aL. Every off-diagonal time kernel in (5) is therefore bounded by n^(-1) times a sum of exp(-cL|cos k|). The zeros of cos k are simple, with |cos k| bounded below by a constant times distance to those zeros. On either equally spaced grid this sum is bounded by a convergent geometric series independent of L. Hence those entries are O_(a,N)(1/L). In the parity kernel, nonzero energies obey L|epsilon|>=c_*>0, so the extra denominators 1-exp(-beta|epsilon|) stay bounded away from zero at fixed a,N. The same estimate applies.

[Proved derivation] Energy pairing epsilon_(k+pi)=-epsilon_k makes the ordinary diagonal K_ii=1/2. On the zero-free parity grid the diagonal is also 1/2. On the doubled ring the two sites at equal times have even separation L. Pairing k with k+pi makes their equal-time kernel one-half of the complete plane-wave sum, which vanishes for distinct sites. With zero modes excluded its residual is O(1/L). Thus, on each nonsingular grid, K-diag(e) is a diagonal matrix with entries plus or minus 1/2 plus off-diagonal O(1/L) corrections. The determinant has no linear correction because those corrections have zero diagonal. Its relative correction is O_(a,N)(L^(-2)).

[Proved derivation] On the zero-mode parity grid, (7) has two U columns and two V rows, each of size O(L^(-1/2)), and its remaining entries are bounded. Its determinant is O_(a,N)(L^(-2)). The ordinary shifted trace prefactors and the nonzero-mode parity prefactors in (9) are bounded at fixed beta/n, while the full physical shifted Z_n is at least one. Thus the parity projection preserves the relative O(L^(-2)) error. Counting the k=2N or 4N factors 1/2 proves (11).

[Proved] Summing the finitely many records in (11) yields

    S_N/Z_L^2 = 2^(-3N)[1+O_(a,N)(L^(-2))],
    Q_D/Z_(2L) = 2^(-2N)[1+O_(a,N)(L^(-2))],
    G = R_L-N log 2+O_(a,N)(L^(-2)).                        (12)

This proves the observed mechanism. At macroscopic time separations the local recorded weights become equal to leading order, while the full ring free-energy difference remains.

## 5. Dependence on spacing

[Proved, using inherited XX scaling function] C119 computes the full spin scaling function r_XX(y)=lim_(L->infinity)R_L(yL), with v=1. Therefore (12) gives

    lim_(L->infinity) G/N = r_XX(2aN)/N-log 2.               (13)

Its inherited asymptotics are

    r_XX(y)=4 exp[-pi/(2y)]-6 exp[-pi/y]+... ,     y->0,
    r_XX(y)=pi y/4-2 exp[-pi y/4]+... ,           y->infinity.

Thus the small-a limit at fixed N is -log 2, and the large-a limit at fixed N is pi a/2-log 2 plus exponentially small corrections. Taking L to infinity first and then N to infinity also gives pi a/2-log 2. This is an iterated limit; no interchange with N growing with L is asserted.

[Proved distinction] The sign of G is a separate question from bounded cost. In the last iterated limit it changes at a=2 log 2/pi. That sign change is not a maximal admissible spacing. Positive G/N remains finite for every fixed a.

[Numerical] The survey evaluates (10) for L=4,8,32,128,512,2048,8192; a=0.125,0.25,0.5,1,2,4,16; N=1,2,4; and additional N=6 cases. At a=1,L=8192 the values G/N for N=1,2,4,6 are 0.685705292, 0.839919367, 0.876720683, 0.877622285. The iterated-limit value is pi/2-log 2=0.877649146. CHECKS.md records the remaining data and calibration.

## 6. An exact general temporal transfer reduction

### 6.1 Positive operators on Hilbert-Schmidt space

[Proved] Let K_S=H_L tensor I+I tensor H_L on two independent copies, and K_D=H_(2L) on the reconnected ring. These spaces have the same physical dimension d^(2L). Define

    B_S=exp(-tau K_S/2),       B_D=exp(-tau K_D/2),
    J_S(X)=sum_s (P_s tensor P_s) X (P_s tensor P_s),
    J_D(X)=sum_(s,t) P_(s,t) X P_(s,t),
    A_j(X)=B_j X B_j,
    T_j=J_j A_j J_j,          u_j=J_j(B_j),    j=S,D.         (14)

The Hilbert-Schmidt inner product is <X,Y>=Tr(X^dagger Y). Both J maps are orthogonal projections for this inner product. J_D is a complete pinching. J_S keeps only the equal-label blocks, so it is not unital. Treating it as a complete pinching would change S_N.

[Proved] Each T_j is self-adjoint and positive semidefinite as an operator on Hilbert-Schmidt space. Indeed,

    <X,A_j X>=Tr(X^dagger B_j X B_j)
              =||B_j^(1/2) X B_j^(1/2)||_2^2>=0,

and compression by J_j preserves both properties. T_j also preserves the cone of positive matrices, because it is a composition of completely positive maps.

[Proved] Expanding the Kraus labels gives the exact identities

    S_N=<u_S,T_S^(N-1)u_S>,
    Q_D=<u_D,T_D^(N-1)u_D>.                                (15)

To check the endpoints, the initial half interval sends I to B_j, the last half interval is the linear functional X->Tr(B_j X), and each interior interval is A_j. Move the first and last J_j onto B_j using self-adjointness and use J_j^2=J_j between intervals. This gives (15). On the numerator, every shared Kraus label factors into the same label on the two copies, producing q_s^2. On the denominator, the two labels at each step are independent, giving Q_D.

### 6.2 Spectral rate and boundary factor

[Proved] Let lambda_j be the largest eigenvalue of T_j, let Pi_j be its full top eigenspace projection, and set

    w_j=||Pi_j u_j||_2^2,
    C_L(tau)=log[||u_S||_2^2/w_D],
    kappa_L(tau)=log(lambda_S/lambda_D).                     (16)

The full top eigenspace is used, so degeneracies cause no basis dependence. Its weight is nonzero. One elementary justification is to maximize the Rayleigh quotient on Hermitian matrices. Replacing X by |X| does not decrease it, because positivity of T makes Tr[X_+ T(X_-)] nonnegative, while ||X||_2=|||X|||_2. A positive maximizing eigenmatrix therefore exists. B_D is strictly positive, and J_D(B_D) is strictly positive on every block. Their inner product is positive. The same reasoning works on the selected support for J_S.

[Proved] Spectral decomposition with nonnegative eigenvalues and weights gives

    S_N <= ||u_S||_2^2 lambda_S^(N-1),
    Q_D >= w_D lambda_D^(N-1),
    G <= C_L(tau)+(N-1)kappa_L(tau).                        (17)

Furthermore, at fixed finite L,tau,

    lim_(N->infinity) G/N = kappa_L(tau).                   (18)

Both identities retain the boundary vectors. Positivity of a temporal map alone controls neither its spectral radius relative to the other geometry nor its boundary overlap.

[Open precise sufficient residual estimate] One sufficient completion of the general proof is, for tau M/L<=a,

    lambda_S(tau)/lambda_D(tau) <= exp[b_1(d)],
    ||J_S(B_S)||_2^2 / ||Pi_D J_D(B_D)||_2^2
        <= exp[b_0(d)+b_1(d)].                             (19)

Then (17) proves G<=b_0+b_1 N. For the prescribed ceiling choice of N, the first estimate is needed in the band a/2<=tau M/L<=a whenever N>=2. The shorter N=1 intervals need only the direct bound log(S_1/Q_D)<=b_0+b_1. This distinction avoids imposing an unnecessarily strong rate estimate outside the used band.

[Proved necessity of part of the residual] If the target lemma holds for every N at fixed tau=aL/M, then (18) forces kappa_L(tau)<=b_1+b_0/N in the limit N->infinity, hence kappa_L(tau)<=b_1. The boundary estimate in (19) is sufficient and need not be necessary, since excited transfer eigenvalues can improve Q_D. A top-eigenvalue bound by itself is not a proof for all N.

[Open scope] T_S and T_D act on the full remaining quantum state, not on a d- or d^2-dimensional seam register. Their projectors have few labels but large ranks. A construction reducing (19) to a seam-only estimate with controlled constants is not supplied. Equation (19) identifies the missing spectral and overlap estimates without replacing them by a claim of locality.

## 7. Why the proposed locality and gap inputs do not close (19)

[Proved, inherited from C122] A full operator-order comparison of the two Gibbs operators is too strong. C122 shows on XX

    log c_order >= (2 beta/L)csc(pi/L)-2 log 2

for the minimum c_order in its doubled-ring comparison. At beta proportional to L this grows linearly in L, despite the scalar XX bound (2). Therefore the full order constant cannot supply the O(1) comparison per macroscopic slab by this direct route. This does not exclude a weaker comparison only on the relevant transfer boundary vectors.

[Proved, inherited from C122] Its explicit imaginary-time estimate has constant

    exp[6b+2(exp(2b)-1)],      b=beta M,

and the guaranteed localization width grows exponentially in b. Substituting tau M of order L therefore does not yield a bounded seam column. This is a limitation of the available estimate, not a lower bound on every possible trace-specific localization estimate. The primary reference [Perez-Garcia and Perez-Hernandez](https://pmc.ncbi.nlm.nih.gov/articles/PMC10113381/) develops the applicable complex-time locality setting; it does not provide (19).

[Proved gap-only limitation] Suppose a ring has ground multiplicity g, Hilbert dimension D and first gap delta. Then

    0 <= log[exp(beta E_0)Z/g]
       <= log[1+(D-g)exp(-beta delta)/g].                   (20)

This follows by bounding every excited Boltzmann factor by exp(-beta delta). With D=d^(2L), delta of order v/L and beta of order L/M, the right side need not be bounded. An O(v/L) gap alone gives no control of the number of low-energy states. There is also no lower bound v>=cM in the background assumptions. Slow sectors are allowed.

[Proved] Even a justified replacement by ground projectors would retain exp[beta(E_(2L)-2E_L)], together with seam matrix elements. It therefore retains the unknown Casimir coefficient. A low-energy/high-energy decomposition would have to prove a multiplicative estimate for the relevant contracted spectral weights, such as (19), rather than only a gap or an operator-norm tail bound.

## 8. Status

[Proved] The bounded-cost lemma holds for XX and the inherited free subclasses with explicit constants, for every fixed finite spacing. Multi-pinch history determinants, including zero-mode parity contributions, give (12). The general positive temporal-transfer formulation gives the exact identities (15), the rate formula (18), and the sufficient residual (19).

[Numerical] All completed finite-grid and full-spectrum checks pass. Interacting ED extends to 12 total qubit sites and 8 total SU(3) sites. The SU(3) L=4 comparison is outside its vacuum-compatible subsequence and is kept only as a finite-channel check. It cannot replace a length-6 versus length-12 vacuum comparison.

[Open] No proof of (19) for general interacting chains, no uniform interacting finite-size extrapolation, and no new general bound on A is obtained. Large interacting ED was not run locally or on Typhon. The task forbids network commands, and no cluster connection or submission was attempted. No file outside this folder was changed.
