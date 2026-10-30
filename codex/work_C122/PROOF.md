# C122. Imaginary-time locality and the seam normalization

2026-10-29. Natural logarithms. Claims below distinguish finite-dimensional proofs, published input, infrared assumptions and numerical calculations.

## 1. Setting and the bound that is actually available

[Assumed, inherited] Write H0 = H_L tensor I + I tensor H_L, H1 = H0 + V = H_(2L), Y = exp(-beta H0), X = exp(-beta H1), sigma = Y/Tr Y, rho = X/Tr X. Use the optimized representative with bond norm at most M. The four reconnection terms give ||V|| <= 4M. Let C contain w physical sites, q = d^w, and let D be its complement. C can be disconnected. For the minimal physical seam used here C = {0,L-1,L,2L-1}.

[Proved, inherited from C120 Section 4] The functional

    mu(W) = min Tr Omega,   Omega >= 0,   I_C tensor Omega >= W

is homogeneous, monotone and subadditive on positive matrices. It satisfies

    Tr W/q <= mu(W) <= q Tr W,
    Q_C := log[mu(Y)/mu(X)] = R_L - Hmin(sigma) + Hmin(rho),
    |Q_C-R_L| <= 2 log q.                                      (1)

The conditional min-entropies refer to C|D. C118/C116 give |R_L| <= beta ||V|| <= 4 beta M. Consequently the safe general bound from these inputs is

    |Q_C| <= 4 beta M + 2 w log d.                            (2)

[Proved distinction] A Hamiltonian order H0 >= H1-cI does not imply exp(-beta H0) <= exp(beta c) exp(-beta H1). The exponential is not operator monotone. The often suggested bound Q_C <= beta ||V|| cannot be justified by that argument. Equation (2) is sufficient to state the trivial scale without using the invalid order step. No claim is made here that every possible sharper scalar bound on Q_C is false.

[Proved] For fixed w, a bound on Q_C of the requested form is equivalent, up to 2w log d in the intercept, to the original bound on R_L. Passing to conditional min-entropy has not weakened the unknown scalar estimate. For XX, the already proved C118 bound gives

    Q_C <= 2(w+1) log 2 + (5 pi/3) beta/L                     (3)

on every even L >= 4. For classical diagonal bonds it gives Q_C <= (2w+1) log d. These are inherited scalar consequences, not new interacting estimates from locality.

## 2. The correct operator-order comparison and its obstruction

[Proved] Set t = beta/2 and

    F(t) = exp(t H1) exp(-t H0).

Then

    X^(-1/2) Y X^(-1/2) = F(t) F(t)^dagger,
    c_ord := min{c : Y <= c X} = ||F(t)||^2,
    Q_C <= log c_ord.                                         (4)

Monotonicity and homogeneity of mu prove the last step. Bounds on the half-time expansional and its inverse therefore suffice for an operator comparison. Using the norm of the nonsymmetric full-time product without this congruence does not directly give (4).

[Proved] Normalization yields

    log c_ord = Dmax(sigma || rho) + R_L.

Operator monotonicity of log gives Dmax >= D. The exact Gibbs identity in C118 is D(sigma||rho) = beta Tr(sigma V)-R_L. Hence

    log c_ord >= beta Tr(sigma V).                            (5)

This is a lower bound on every possible full operator-order comparison, independent of which locality theorem is used to estimate it.

[Proved on XX, including the scaling point] Use h = Sx Sx + Sy Sy, M = 1/2, even L >= 4. Spin rotation about z and global spin flip make each one-site Gibbs marginal I/2. The crossed bonds have zero expectation in the split state. Translation invariance on either L-ring gives

    Tr(sigma V) = -2 U_L(beta)/L,
    U_L(beta) = Tr[H_L exp(-beta H_L)]/Z_L.

The settled XX ground energy is E_L = -csc(pi/L). The Gibbs variational principle, compared with a ground state, gives U_L-S_L/beta <= E_L. Since S_L <= L log 2,

    U_L <= E_L + L log 2/beta.

Combining this with (5) proves, for beta > 0,

    log c_ord >= (2 beta/L) csc(pi/L) - 2 log 2.               (6)

At beta = tau L,

    log c_ord >= 2 tau csc(pi/L)-2 log 2
              = (2/pi) tau L - 2 log 2 + O(tau/L).            (7)

The proposed K(d)+k M f(d) tau is bounded as L tends to infinity at fixed tau, whereas (7) diverges. Thus full Gibbs operator order cannot establish the desired scale through (4), even on XX. This is a counterexample to the stronger order lemma, not to the scalar mu-ratio target. The latter is already true for XX by (3).

[Proved checks] At fixed L and beta tending to infinity, (6) has slope 2 csc(pi/L)/L, tending to 2/pi. The true XX doubling defect has slope pi/(4L)+o(1/L). At beta = 0, c_ord = 1; the negative constant in (6) is harmless. Rescaling h by a positive J replaces beta by J beta in all these formulas, so the obstruction has energy scale M, as required.

## 3. What the literature supplies

[Assumed published context] The original source is Huzihiro Araki, [Gibbs states of a one dimensional quantum lattice](https://doi.org/10.1007/BF01645134), Communications in Mathematical Physics 14, 120-157 (1969). The publisher metadata and abstract were checked. The quantitative locality statements used for comparison below were checked in the accessible primary papers of the later authors, rather than inferred from the original abstract.

[Assumed published results] Araki's finite-range 1D analyticity and locality results are formulated explicitly in Perez-Garcia and Perez-Hernandez, Theorems 2.3 and 3.1 of [arXiv:2004.10516v1](https://arxiv.org/pdf/2004.10516). The observable estimate is a positive generating series with coefficients Omega_k^*(4|t|). For nearest-neighbor interactions these coefficients reduce to (4|t| Omega_1)^k/k!. The resulting locality scale grows exponentially with |t| times the interaction strength. Uniformity in volume at fixed t does not give uniformity when t grows with the volume.

[Assumed published results] Bluhm, Capel and Perez-Hernandez, [Proposition 3.2 and Remark 3.3](https://arxiv.org/html/2104.04419v2#S3.SS2), state bounded expansionals and inverse expansionals with factorial spatial tails. Their constant G depends on beta, interaction strength and range. Their later convention absorbs beta into the Hamiltonian. It does not remove this dependence.

[Assumed published results] Kuwahara and Saito, [Lemma 1 and equations (18)-(19)](https://arxiv.org/pdf/1807.08424), obtain an imaginary-time tail controlled by zeta_l = 6 e k t/[gamma log ceil(l/k)] in their unit-strength convention. The bound requires zeta_l < 1. They choose a length exponential in beta for this step. Their belief-propagation insertion has a different, polynomial locality length, but their comparison still requires the imaginary-time step. Neither estimate is a beta/L bound on mu.

[Proved scope] These are sufficient estimates. Their exponential length scales do not prove that the actual support radius grows exponentially in every chain. The obstruction (7), in contrast, is a lower bound on the actual order constant.

## 4. An explicit nearest-neighbor expansional bound

[Proved] A direct counting argument gives an explicit Araki-type constant for this reconnection without importing unspecified constants. Fold the two chains into L rungs (a_i,b_i), each of dimension d^2. Both H0 and H1 are nearest-neighbor Hamiltonians on this rung circle, with each rung bond of norm at most

    J = 2M.

Their difference V is supported on the two seam rungs and has norm v <= 4M. Scalar shifts cancel and the optimized coboundaries telescope on both rings. The following counting uses L >= 3; the length-two numerical calibration is not needed for this theorem.

[Proved] More generally let H be a nearest-neighbor circle Hamiltonian with bond norms <= J. Let V be supported on a connected interval of s sites, with ||V|| <= v. Put

    F_V(t) = exp[t(H+V)] exp(-tH),    t >= 0.

Direct differentiation gives

    F_V'(t) = [H,F_V(t)] + V F_V(t),    F_V(0) = I.             (8)

Expand the finite-dimensional superoperator exp[t(ad_H+left_V)] applied to I. Regard I as supported in the initial s-site interval. An internal commutator does not grow the support. After k growth events, the support lies in an interval of length at most s+k. There are at most s+k internal bonds, including the full-circle endpoint case. Their total superoperator norm is at most 2J(s+k). Multiplication by V adds at most v. At most two boundary bonds can enlarge the interval, each by one site, with combined norm at most 4J. This bound remains valid when opposite ends meet, because the actual number of growing bonds then decreases.

[Proved] Grouping Taylor words by the number k of boundary growth events is dominated coefficient by coefficient by

    a_0(0)=1,   a_k(0)=0 for k>0,
    a_k' = [v+2J(s+k)] a_k + 4J a_(k-1),    a_(-1)=0.

Solving this triangular system gives

    a_k(t) = exp[(v+2Js)t] [2(exp(2Jt)-1)]^k/k!.

For k=0 this is immediate. Differentiating the displayed expression for k>0 gives the preceding recurrence, so it holds at every k. Summing the nonnegative majorant yields

    ||F_V(t)|| <= exp[(v+2Js)t + 2(exp(2Jt)-1)].               (9)

The sum converges for every finite t and independently of the circle length. The same argument applies with -t, since changing signs does not affect the norm counting. It also applies with H+V as reference when that Hamiltonian has the same bond norm bound, as it does here.

[Proved] For the ring seam, s=2, J=2M, v=4M, t=beta/2. With b = beta M, define

    G(b) = exp[6b + 2(exp(2b)-1)],
    B(b) = 2 log G(b) = 12b + 4(exp(2b)-1).

Both half-time products needed for the two order directions have norm <= G(b). Therefore

    |Q_C| <= B(beta M),
    |Q_C| <= min{4 beta M+2w log d, B(beta M)}.                (10)

G is a double exponential upper bound in b; its logarithm is exponential. No optimal growth law is claimed. At small b, B(b)=20b+8b^2+O(b^3). At large b, (2) is the better scalar bound. Neither term in (10) stays bounded at beta=tau L. The scalar partition ratio itself retains the stronger settled estimate |R_L|<=4 beta M, so combining (10) with the entropy ceiling does not improve the general bound on R_L.

## 5. Enlarging the seam and optimizing its width

[Proved] The same counting controls spatial truncation. Let I_r contain the original two seam rungs and r extra rungs on each side, before they cover the circle. Drop Hamiltonian bonds outside I_r, keeping V. Denote the truncated half-time product by F_r. Taylor terms that never reach the omitted bonds agree. A term that differs must have at least r growth events. Thus a safe tail bound is

    ||F-F_r|| <= 2 exp(6b) sum_(k>=r) a^k/k!,
    a = 2(exp(2b)-1).                                        (11)

The factor two bounds the two discarded tails. This is a sufficient bound, with deliberately nonoptimal endpoint counting.

[Proved] If r >= 2 e a, Stirling's elementary bound k! >= (k/e)^k gives

    sum_(k>=r) a^k/k! <= 2^(1-r).

For a fixed absolute error epsilon it is enough that, in addition,

    r >= [6b+log(4/epsilon)]/log 2.

The resulting guarantee has r=O(exp(2b)+log(1/epsilon)). It cannot provide a small fraction of the circle when b=M tau L. This conclusion concerns what (11) certifies, not a necessary spatial width for every model.

[Proved] There is also a direct way to see the cost of replacing an expansional by a local congruence in mu. The SDP dual is

    mu(W) = max Tr(TW),    T>=0,    Tr_C T <= I_D.             (12)

Strict primal feasibility proves strong duality. For fixed W>=0 the function sqrt(mu(AWA^dagger)) is a seminorm in A, because it is the supremum over feasible T of ||sqrt(T) A sqrt(W)||_2. If A acts only on C, primal feasibility gives

    mu(AWA^dagger) <= ||A||^2 mu(W).

For arbitrary A the dimension bounds in (1) instead give

    mu(AWA^dagger) <= q^2 ||A||^2 mu(W).

Let B = exp(-tH0) exp(tH1) = F^dagger, so Y=B X B^dagger, and let B_r=F_r^dagger be supported on C=I_r, containing w<=4+4r physical sites. With epsilon_r=||B-B_r||, the triangle inequality proves

    Q_C <= 2 log(||B_r||+q epsilon_r).                        (13)

The q multiplying the approximation error must not be dropped. In the range where I_r has its full size, q=d^(4+4r). From (11), a sufficient choice r>=2 e a d^4 makes q epsilon_r <= 4 d^4 exp(6b) 2^(-r). This remains an exponential-in-b width guarantee. The entropy cost is 2w log d, while ||B_r|| is still bounded by G(b), not by a beta/L quantity.

[Proved limitation] Merely enlarging C does not replace V by a term of norm M/w. For the exact order route, c_ord is independent of C altogether. Even if one assumed a more favorable but unproved cost c beta M/w, balancing it against 2w log d would give

    inf_(w>0) [2w log d+c beta M/w]
      = 2 sqrt(2 c beta M log d),

which diverges as sqrt(L) at beta=tau L. A fixed-width entropy ceiling and a local spreading estimate do not by themselves produce the desired finite scaling limit.

## 6. Low-energy restriction and what a ring gap does not control

[Proved, quantitative extension of C120] Let E_i, P_i and g_i=Tr P_i be the ground energies, projectors and degeneracies of H_i. Write

    exp[-beta(H_i-E_i)] = P_i + T_i(beta),
    eta_i(beta) = log[mu(P_i+T_i)/mu(P_i)] >= 0.

Then the exact identity is

    Q_C = beta(E1-E0) + log[mu(P0)/mu(P1)] + eta_0-eta_1.      (14)

If delta_i is the first positive gap and the total Hilbert dimension is N=d^(2L), then 0<=T_i<=exp(-beta delta_i)I. Since mu(I)=N/q and mu(P_i)>=g_i/q, monotonicity and subadditivity give

    0 <= eta_i <= log[1+(N/g_i)exp(-beta delta_i)].            (15)

Another useful bound is eta_i <= log[1+q^2 Tr(T_i)/g_i]. Both require control of the spectral tail for uniform use. At fixed L they prove eta_i -> 0 as beta tends to infinity.

[Proved obstruction to the gap-only argument] If delta_i is of order v/L and beta=tau L, the gap suppression is exp(-const tau v), while N grows exponentially in L. Equation (15) is then not small. This is not merely a defect of a proof based on the smallest eigenvalue: finite spectra with one level at zero and N-1 levels at delta have tail trace (N-1)exp(-beta delta). Thus the gap alone contains no bound on the multiplicity of the tail. This spectral example is not claimed to be a fixed-bond nearest-neighbor lattice counterexample. Additional microscopic input would have to supply the missing spectral information.

[Proved] With unique ground states, mu(P_i)=(Tr sqrt((P_i)_C))^2. Its value lies between 1 and min(dim C,dim D), so the ground offset in (14) is bounded. Nevertheless, even an exact replacement of both Gibbs operators by their ground projectors leaves beta(E1-E0). Low-energy restriction cannot determine that energy difference from the gap. In particular, a bound on excited-state counting is not a bound on the relative vacuum normalization.

[Assumed infrared input, proved consequence] For vacuum-compatible lengths, E1-E0=pi A/(4L)+o(1/L). Fixed C and (14) reproduce the C120 limit

    lim_(beta->infinity) Q_C/beta = E1-E0.

If the fixed-aspect-ratio thermal limit assumed in C116 exists, (1) also shows that Q_C and R_L have the same large-tau slope after taking L to infinity, because their difference divided by tau is bounded by 2w log d/tau. Neither statement allows the order of limits to be exchanged from a gap estimate alone. The permitted models have no dimension-only lower bound on sector velocities, and no dimension-only bound on sector count is used here.

## 7. The remaining lemma

[Open, sufficient formulation] Fix a seam of w0 sites, independent of L. Let F_C be the feasible dual set in (12). A genuinely weaker replacement for operator order would construct, for each dual maximizer T0 of mu(Y), a T1 in F_C such that

    Tr(T0 Y) <= exp[K'(d)+k M f(d) beta/L] Tr(T1 X).           (16)

Then (12) and the entropy bound prove R_L <= K'(d)+2w0 log d+k M f(d) beta/L and hence A <= (4k/pi) M f(d), under the inherited infrared assumption. It is enough for the A conclusion to prove an asymptotic version at beta=tau L with an additional exponent epsilon_L(tau)+u(tau), where epsilon_L(tau) tends to zero as L tends to infinity at each fixed tau and u(tau)/tau tends to zero afterward.

[Open] Equation (16) is not proved by the cited locality results. Allowing the dual witness to change is crucial: the stronger comparison against every positive test operator is operator order and fails by (7). A useful new lemma must exploit the partial-trace constraint on the witness and the cancellation of the two vacuum normalizations. Writing (16) is a precise specification of the missing work, not an independent solution of it.

[Proved conditional alternative] Under unique ground states, a bound L(E1-E0)<=k M f(d), together with a uniform upper bound on eta_0 in (14), would suffice, because the ground offset is bounded by log q and -eta_1<=0. The first of those hypotheses is the original finite-size Casimir estimate. Restricting to a low-energy subspace has not removed it.

## 8. Numerical method and checks

[Proved method] All tested Hamiltonians are real and invariant under independent global rotations of each conserved color. Twirling a primal optimizer by these rotations, then averaging with its complex conjugate, preserves the objective and feasibility. Therefore Omega can be restricted to real symmetric blocks at fixed complement color counts without changing the full SDP optimum. The physical positive-semidefinite constraint still includes every total-color sector. No thermal sector or excitation is discarded.

[Numerical] sdp.py minimizes Tr Omega-t sum_blocks log det(I_C tensor Omega-W) along a decreasing barrier path with Newton steps. For each reported result it checks the minimum slack eigenvalue and adds a scalar repair. The positive dual matrices t times the inverse slack are rescaled until Tr_C T<=I. The resulting primal and dual objective values bracket the optimum to the reported floating-point accuracy. These are numerical feasibility checks, not directed-rounding certificates.

[Numerical] checks.py uses W_i=exp[-beta(H_i-E_i)] to avoid overflow and restores log mu(exp(-beta H_i))=-beta E_i+log mu(W_i). The full spectra are retained. It independently checks beta=0, where mu(I)=dim D, and every unique-ground-state SDP against the Schmidt-value formula. It checks (1), the entropy interval and the conditioned full-order comparison (4)-(5). Results and limitations are in CHECKS.md.

[Open conclusion] The new obstruction excludes a full-norm order proof, and the explicit locality bound explains the nonuniformity in beta. No universal improvement over the O(beta M) scale for the scalar ratio has been proved. There is no counterexample here to Q_C<=K'(d)+k M f(d) beta/L or to the proposed bound on A.
