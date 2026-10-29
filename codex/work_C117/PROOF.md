# C117 proof record

2026-10-29. Claims are tagged. Lattice spacing and hbar are one.

## 1. Setting and scope

[Assumed] H_L = sum_i h_i is a one-site translation-invariant nearest-neighbor ring Hamiltonian on local dimension d. Choose a centered gauge representative with norm M, where M is the gauge-optimized half spectral width. On vacuum-compatible lengths the infrared assumption is

    E_0(L) = e_infinity L - pi A/(6L) + o(L^(-1)),
    A = sum_a c_a v_a.

For unequal chiral speeds, v_a denotes their average. Compatible lengths and boundary sectors must be fixed before using any finite-size estimator. The campaign's many slow sectors are allowed throughout. Neither their number nor their total central charge is bounded here.

[Open] No general inequality A <= C M f(d) is proved in this work. The results below are a local-bootstrap obstruction, certified small energy bounds, an exact local moment bound with a massive contamination example, and conditional reductions specifying the missing inequalities.

## 2. The marginal and ground-state SDPs

[Definition] A connected-block relaxation at level k uses a density matrix rho on k sites and the constraints

    rho >= 0, Tr rho = 1,
    Tr_first rho = Tr_last rho.

The last equality identifies the ordered (k-1)-site density matrices. Repeated tracing imposes consistency on all connected overlaps. Assign this same rho to all cyclic translates of a k-site block on the ring. When L >= 2k, intersections of these blocks do not split into two components, and this imposes all their nonempty overlap constraints. For smaller rings there can be extra disconnected-overlap constraints. Omitting them still gives a valid relaxation, but is not full ring extendibility. A global ring density matrix is imposed in the separate k=L reference calculation.

[Proved] Set C = h_12 tensor I_(k-2). Minimizing Tr(C rho) gives an energy-density lower bound ell_k for every ring for which these local constraints are valid. In particular E_0(L) >= L ell_k. A translation average of an actual ring ground state is a feasible point, so no symmetry assumption on the ground vector is needed.

[Definition] For a basis O_a supported strictly inside the k-site block, put

    D_a = [H,O_a],
    G_ab = (O_a^dagger D_b + D_a^dagger O_b)/2.

Here [H,O_a] includes both boundary bonds touching its support. All terms fit inside the block. Impose stationarity for the needed products and impose the matrix inequality

    G(rho)_ab = Tr(rho G_ab) >= 0.

[Proved] In a stationary state,

    <O_a^dagger[H,O_b]> - <[O_a^dagger,H]O_b>
        = <[H,O_a^dagger O_b]> = 0.

Thus G(rho)_ab equals <O_a^dagger[H,O_b]>. For a ground eigenvector Omega and O = sum_a z_a O_a,

    z^dagger G(rho) z
      = <O Omega,(H-E_0) O Omega> >= 0.

A ground-state mixture has the same property. The ground-state conditions are necessary convex constraints; arbitrary feasible truncated marginals need not extend to a state on the line or the ring.

[Proved dual certificate] Let X be Hermitian on k-1 sites and let S_j be Hermitian operators whose expectation vanishes by stationarity. Define

    D(X) = X tensor I - I tensor X,
    W >= 0,
    K = C + D(X) + sum_j y_j S_j - sum_ab W_ba G_ab.

If K >= b I, then for every feasible rho,

    Tr(C rho) = Tr(K rho) + Tr(W G(rho)) >= b.

Hence E_0(L) >= L b. This implication does not require optimizer convergence, strict feasibility, or strong duality.

[Proved qualification] With W=0, the local positive operator C+D(X)+sum yS-bI gives an operator sum-of-squares certificate after translation summation if the stationarity terms also sum to a zero operator, or if they are absent. In the marginal-only calculations y=W=0, so this cancellation is automatic for D(X). With W nonzero, the localizing term G is positive in ground-state expectation, not necessarily as an operator. Such a certificate is a ground-state certificate and is not automatically a literal SOS decomposition of H_L-LbI. This distinction matters for C114.

## 3. Why fixed-support local programs cannot resolve the Casimir term

[Proved, precise class] Consider a translation-invariant hierarchy whose variables and constraints are translates of operators and products fitting in intervals of at most R sites. Include any local density positivity, local stationarity, and local ground-state Gram matrices of this kind. Do not include a Gram matrix coupling arbitrarily distant translates, winding operators, or explicit circumference-dependent constraints. For L > 2R its feasible set and its objective per site are independent of L. Denote the optimum ell_R.

The reason is algebraic. Every overlap, commutator, product, and partial trace occurring in the program can be evaluated on a line segment without identifying distinct sites. Translation averaging removes the location label. The matrices specifying the SDP then contain h, d, and the finite support patterns, but contain no L.

[Proved shorter-ring bound] For any n > 2R, take a translation-invariant ground-state density matrix on an n-site ring. Its local data satisfy the same program. Therefore

    ell_R <= E_0(n)/n.                                      (1)

The statement applies to the strengthened ground-state program as well. The shorter ring satisfies all the stated ground-state inequalities, not only marginal positivity.

[Proved under the assumed Casimir expansion] If A>0, choose a sufficiently large, finite, vacuum-compatible n > 2R. Then

    e_infinity - ell_R >= e_infinity-E_0(n)/n
                          = pi A/(6n^2) + o(n^(-2)) > 0.    (2)

Consequently the direct upper estimator supplied by this lower bound is

    (6/pi) L [e_infinity L - L ell_R]
      = (6/pi) L^2 (e_infinity-ell_R),                       (3)

which diverges with L. Knowing e_infinity exactly does not repair this. A fixed-support relaxation may fit a few short-ring values, but cannot isolate the asymptotic 1/L correction.

[Proved conditional extension] Suppose R=R(L) tends to infinity with R=o(L), and vacuum-compatible shorter lengths n(L)>2R can be chosen with n=O(R). Then

    E_0(L)-L ell_R >= E_0(L)-L E_0(n)/n
      = (pi A/6) [L/n^2-1/L] + o(L/n^2)+o(1/L).             (4)

Multiplying by L shows that the error cannot be o(1/L). The main positive term grows as (L/n)^2. The compatible-length hypothesis is explicit; it holds for the usual even or fixed-modulus vacuum sequences. For bounded R, (2) already proves the obstruction. This is an obstruction to these local programs, not to SDPs with long-distance or winding information.

[Proved sufficient certificate lemma, not established universally] Let ell_L be a certified lower bound for E_0(L), and u_2L a certified upper bound for E_0(2L). Assume both L and 2L belong to the vacuum sequence, and suppose one can prove

    limsup_(L->infinity) L [u_2L-2 ell_L] <= K M f(d).        (5)

Since

    E_0(2L)-2E_0(L) <= u_2L-2 ell_L,
    E_0(2L)-2E_0(L) = pi A/(4L)+o(1/L),

(5) gives the actual target inequality

    A <= (4K/pi) M f(d).                                    (6)

[Open] The missing lemma is (5) with a universal finite K and f. Two lower bounds cannot be subtracted to obtain the required upper bound on the doubling defect. A convergent hierarchy without a suitable rate does not prove (5). This is a concrete ring-certificate target for C114, rather than a claim that the existing local certificates have the needed rate.

## 4. Exact small-block examples and dual structure

[Proved] For spin-1/2 Heisenberg h=S.S, the two-site relaxation has optimum -3/4. The singlet marginal has equal partial traces and attains the smallest bond eigenvalue.

For k=3 choose X=-h/2. Then

    C+D(X) = (h_12+h_23)/2.

The three-spin open Hamiltonian has eigenvalues -1, 0, 1/2. One derivation writes it as S_2.(S_1+S_3). Coupling the outer spins to spin zero gives zero. Coupling them to spin one and then to total spin 1/2 or 3/2 gives -1 or 1/2. Therefore

    (h_12+h_23)/2 + I/2 >= 0,
    E_0(L) >= -L/2.                                       (7)

The four-site ring has energy -2, also obtained from the two-dimensional block below at Delta=1. Its translation-invariant three-site marginal is feasible and has energy density -1/2, including all one-site ground-state Gram constraints. Thus both the basic and strengthened three-site programs have optimum -1/2. In particular, the failure to improve this example is not an optimizer failure.

[Proved] For h=Sx.Sx+Sy.Sy+Delta Sz.Sz with Delta=1/2, the open three-site Hamiltonian has a fully polarized eigenvalue Delta/2. In each one-down or one-up sector its matrix is

    [[0, 1/2, 0], [1/2, -Delta/2, 1/2], [0, 1/2, 0]].

The antisymmetric endpoint vector has eigenvalue zero; the other eigenvalues are

    (-Delta +/- sqrt(Delta^2+8))/4.

The same choice X=-h/2 yields

    ell_3 >= (-Delta-sqrt(Delta^2+8))/8
           = -(1+sqrt(33))/16
           = -0.4215351654086268... .                      (8)

For completeness the four-site zero-magnetization symmetric subspace has the matrix

    [[-Delta, sqrt(2)], [sqrt(2), 0]],

in the normalized sums of the two alternating and four adjacent-pair configurations. Its lower eigenvalue is (-Delta-sqrt(Delta^2+8))/2. Equation (8) is already a universal lower bound on the four-site energy divided by four, so this eigenvector must attain the ground energy. Its three-site marginal also satisfies the ground constraints. Both three-site optima therefore equal (8).

[Proved] For the SU(3) swap h=P, P>=-I. The completely antisymmetric three-qutrit state is an eigenvector with eigenvalue -1 for all three swaps. Its marginals attain ell_2=ell_3=-1, and its three-site density matrix satisfies the one-site ground constraints. Thus these small levels do not improve the trivial swap bound. No statement about the large-L Casimir coefficient follows.

[Proved SOS form] Summing the translates of the positive three-spin operator in (7) gives

    H_L + (L/2)I = sum_i [(h_i+h_(i+1))/2 + I/2].

Taking a positive square root of each summand gives a local SOS. The coefficient is independent of L, which is the structural problem. The analogous XXZ certificate follows from (8). Full-ring certificates can depend on L, but their factors need not have bounded support.

## 5. Implementation and exact verification

[Numerical procedure] bootstrap.py searches the dual with scipy L-BFGS-B. It smooths lambda_min(K) to -tau log Tr exp(-K/tau), and uses W=B B^T. This parametrization makes the search nonconvex, even though the underlying SDP in W is convex. Search outputs are not asserted optimal unless Section 4 supplies a matching feasible point. All tested bonds are real. Taking a real ground-state density matrix by conjugation averaging is sufficient for these bonds; this is a computational restriction on the examples, not an assumption in the general results.

[Numerical scope] The local tests use k=2 and 3, all one-site matrix units for O at k=3, and Heisenberg, XXZ at Delta=1/2, SU(3) swap, and one seeded real symmetric random bond at each d=2,3. No internal spin symmetry is imposed in the solver. Local PSD matrices have dimension at most 27. Full-ring reference SDPs have k=L, with L=4,6 for d=2 and L=3,4 for d=3, so their matrices have dimension at most 81. The SU(3) L=4 ring is a finite-energy diagnostic, not a vacuum-sector Casimir datum. No infrared theory is assumed for the random bonds.

[Proved exact-arithmetic verification] Bond entries are dyadic. Round X and stationarity coefficients to multiples of 2^(-16), and B to multiples of 2^(-12). Reconstruct K exactly as an integer matrix divided by 2^32. A floating eigensolver only proposes a trial scalar e and a Cholesky factor Q. Quantize Q to denominator 2^24. Compute the residual exactly,

    R = K-eI-Q Q^T,
    delta = max_i sum_j |R_ij|.

For real symmetric R, every eigenvalue is at least -delta. Therefore

    K >= (e-delta) I.                                    (9)

All entries of the residual have denominator 2^48, and the stored lower endpoint is a rational number with this denominator. verify.py independently reconstructs K from h, X, y, and B using integer arithmetic and checks (9). It does not trust a floating eigenvalue or a solver status. W remains positive because it is reconstructed as B B^T.

[Proved explicit SOS conversion] R+delta I is symmetric diagonally dominant with nonnegative diagonal. It decomposes into terms

    |R_ij| (e_i+sign(R_ij)e_j)(e_i+sign(R_ij)e_j)^T, i<j,

and diagonal remainders R_ii+delta-sum_(j!=i)|R_ij| >= 0. Together with Q Q^T, this gives a directly checkable positive decomposition of K-(e-delta)I. For the ground-localizing certificates, the qualification in Section 2 still applies.

[Proved reference intervals] ring_checks.py applies the same method to the entire finite-ring H_L. It also quantizes a trial ground vector v. The upper endpoint v^T H_L v/(v^T v) is evaluated as an exact rational quotient. Thus both endpoints of each of the ten stored finite-ring intervals are certified. These full-marginal SDPs are equivalent to diagonalizing H_L; they are reference calculations, not an efficient large-ring algorithm.

[Proved verification outcome] All fifteen local certificates and all ten ring intervals pass the integer checks. Their decimal summaries and the small floating residual checks appear in CHECKS.md. The largest reference interval width is below 3.83e-6. Random d=3 k=3 optimizations reached the iteration cap at the last three smoothing parameters, both with and without ground constraints. The certificates remain valid; optimality is not established. The tiny changes between ground and marginal searches do not demonstrate a useful improvement from the ground constraints.

## 6. An exact moment bound, and why its q^4 coefficient is not infrared data

[Definition] Use the C100 convention

    h_q = sum_j exp(-iqj) h_j,
    m_1(q) = <[h_-q,[H,h_q]]>/(2L).

At zero temperature this is the first moment of the positive, momentum-averaged absorption spectrum. C100 proves 0<=m_1(q)<=8 M^3 q^2 for allowed momenta. The following local expansion refines a different part of that information.

[Proved] Write b_r=<[h_0,[H,h_r]]>. If |r|>2 the operator is zero. Jacobi and stationarity imply

    <[h_0,[H,h_r]]> = <[h_r,[H,h_0]]>,

and translation invariance gives b_r=b_-r. Also sum_r b_r=0, since sum_r h_r=H. Consequently on rings large enough to avoid overlapping images, or on the infinite chain,

    m_1(q) = b_1(cos q-1)+b_2(cos(2q)-1).                 (10)

This identity on a ring is asserted at allowed q. Its trigonometric polynomial extension need not be nonnegative at arbitrary real q on a fixed finite ring.

There are at most two nested commutators in b_1 and one in b_2. Each has norm at most 4 M^3. Thus

    |b_1| <= 8 M^3,    |b_2| <= 4 M^3.

Expanding (10),

    m_1(q) = a_2 q^2+a_4 q^4+O(q^6),
    a_2 = -b_1/2-2b_2,
    a_4 = b_1/24+2b_2/3,
    |a_4| <= 3 M^3.                                      (11)

This is an unconditional local coefficient inequality, with explicit constant 3. It is not a bound on A or on sum c_a v_a^3. In particular large-separation nested energy commutators vanish as operators; they cannot be used to remove a massive contribution by themselves.

[Proved within the continuum stress-tensor dictionary, inherited from C100] At T=0 the stress-tensor contribution is

    Sbar_h(q,omega) = sum_a c_a v_a^2 |q|^3/(24pi)
                     delta(omega-v_a|q|),
    m_1,CFT(q) = |q|^4 sum_a c_a v_a^3/(24pi).              (12)

The total lattice moment is not just (12). Taking four momentum derivatives also probes analytic ultraviolet and massive contributions.

### A massive example in an optimized bond gauge

[Proved model facts] Consider the transverse-field Ising bond

    h_i = -J X_i X_(i+1) - (g/2)(Z_i+Z_(i+1)), 0<J<g.

Its two-site eigenvalues are +/-sqrt(g^2+J^2) and +/-J. The extremal eigenvectors have identical left and right one-site marginals because they are symmetric under interchange. Therefore their energy difference is unchanged by every coboundary o tensor I-I tensor o, and

    M = sqrt(g^2+J^2).

The standard fermion solution has dispersion 2 sqrt(g^2+J^2-2gJ cos k), hence a nonzero gap 2(g-J), and A=0. The exact solution also ensures analyticity of local ground-state data near J=0. The spectral and phase statement is the conventional Ising solution of [Pfeuty](https://www.math.ucdavis.edu/~bxn/pfeuty1970.pdf), with Pauli normalization used here.

[Proved leading perturbative coefficient] At J=0 the ground state is all Z spins up. Denote by |i,i+1> the state with two neighboring spins down. First-order perturbation gives

    |Omega> = |0> + [J/(4g)] sum_i |i,i+1> + O(J^2).

This expression is used for fixed-order matrix elements, not as a uniformly normalized infinite-volume vector. The matrix element of h_q to |i,i+1>, including the first-order correction of the excited vector needed for orthogonality, is

    exp(-iqi) [-J + (J/(4g)) g(2+2cos q)]
      = (J/2)(cos q-1) exp(-iqi) + O(J^2).                 (13)

To see the second term directly, the field part of h_q changes by g(1+exp(iq))[exp(-iqi)+exp(-iq(i+1))] when those two spins are flipped. This is g(2+2cos q) exp(-iqi). At nonzero q the vacuum expectation of h_q vanishes, so subtracting the vacuum component of the excited state does not add another term.

The states |i,i+1> are orthonormal and cost energy 4g at this order. In the normalization of m_1, summing their squared matrix elements gives

    m_1(q) = g J^2 (1-cos q)^2 + O(J^4/g),
    a_4 = g J^2/4 + O(J^4/g) > 0                          (14)

for sufficiently small nonzero J/g. The even expansion follows from the staggered Z rotation taking J to -J on the infinite chain or an even ring. Analyticity of the finitely supported b_r permits expanding their q^4 coefficient term by term. Since (10) is local, the leading coefficient survives the thermodynamic limit.

[Proved implication] A strictly gapped nearest-neighbor chain, in an optimized bond gauge, contributes at order q^4 while its Casimir coefficient vanishes. Thus matching a_4 to (12) is false, even after focusing on the order where the vacuum CFT term occurs. This example does not disprove any upper bound on A. It identifies a failed identification, not a general no-go for every combination of moments.

[Numerical check] moment_check.py on a six-site ring finds a_4/(gJ^2)=0.2499414, 0.2497658, 0.2490645 at J/g=0.025, 0.05, 0.1. The polynomial (10) agrees with direct double commutators at allowed q to below 3e-18. Small negative a_2 values on these finite rings are not positivity violations because nonallowed momenta are not physical ring Fourier modes. These checks verify the leading term; the analytic argument establishes the massive contamination.

[Proved additional obstruction] Even a hypothetical upper bound on C_3=sum c_a v_a^3 would not bound A without more input. Algebraically, n sectors with c_a=1 and v_a=V/sqrt(n) have C_3=V^3/sqrt(n) but A=V sqrt(n). This is a counterexample to that inference, not a construction of a fixed-d microscopic chain.

## 7. A positive infrared filter that does isolate the desired weight

[Definition] Take the infinite-volume limit first. For a positive number Lambda exceeding every chiral propagation speed define

    F(q) = integral_(0, Lambda|q|] Sbar_h(q,omega) d omega/omega.

At finite volume this uses the positive spectral projector

    P_(0,Lambda|q|)(H-E_0) (H-E_0)^(-1),

with the inverse restricted to strictly positive energies, inserted between the two h_q operators and averaged over +/-q with factor 1/(2L). Degenerate zero-energy states are excluded. Then take compatible thermodynamic limits before q tends to zero.

[Proved spectral property] Any spectral contribution whose threshold is a fixed Delta>0 vanishes from F(q) once Lambda|q|<Delta. This statement concerns excitations with a positive threshold. Massive virtual processes can still change the infrared parameters; they are not an additional positive-threshold spectral contribution. No uniform lower bound on Delta is used.

[Assumed lattice-to-continuum identification] Assume the low-energy energy-density spectral weight contains the stress-tensor contribution with the normalization (12), up to o(q^2) in the filtered inverse moment. An asymptotic ring energy formula alone does not prove this operator matching. Density improvements can add other contributions. For the sufficient conclusion it is enough that the stress contribution gives

    liminf_(q->0) F(q)/q^2 >= A/(24pi).                     (15)

For a pure stress-tensor scaling limit there is equality. With unequal chiral speeds the two weights are c_a v_R^2/(48pi) and c_a v_L^2/(48pi); dividing by their respective frequencies gives c_a(v_R+v_L)/(48pi)=c_a v_a/(24pi). Thus the same A appears without parity symmetry.

[Assumed filtered inverse-moment lemma; proved implication] If the microscopic theory satisfies

    limsup_(q->0) F(q)/q^2 <= K M f(d),                    (16)

then (15) proves

    A <= 24pi K M f(d).                                   (17)

For equal chiral speeds, the inherited Lieb-Robinson input v<=4eM permits Lambda>4eM. If only the recorded bound on the average speed is used in an unequal-speed sector, each positive chiral speed is at most twice the average, so Lambda>8eM suffices. A propagation bound applied separately to each chirality can improve this choice. The cutoff does not specify a matching temperature or finite crossover length.

[Open] Equation (16) is not proved. Positivity and the known upper bounds on nonnegative spectral moments do not supply an upper bound on this inverse moment. For a positive measure with a weight w at energy epsilon, its inverse moment is w/epsilon while its positive moments are w epsilon^r. Moving weight toward zero keeps the latter bounded while increasing the former. More formally, a fixed polynomial in nonnegative powers of omega cannot majorize 1/omega on an interval reaching zero. Cauchy-Schwarz gives the opposite direction, m_(-1)>=m_0^2/m_1.

[Proved scope qualification] This moment argument rules out the inference from those moment upper bounds, not a future theorem using the entire family of local ground-state constraints. The full, untruncated family characterizes ground states in the relevant setting. A useful new ingredient must control the inverse Liouvillian on the shrinking infrared window, or produce ring certificates such as (5). Local norm estimates alone have not done so.

## 8. Conformal and modular bootstrap with lattice inputs

[Proved information content] The space of k-site operators has dimension d^(2k). It bounds the rank of a specified linear map from those operators into infrared fields. It does not bound how many independent fields have nonzero overlap with one image vector. In particular a single bond can contain a linear combination of many stress tensors. To count all of them one needs a surjectivity or independence statement with a quantitative support bound.

[Proved representation-theory distinction] One stress tensor spans one dimension irrespective of c. Its level-two vacuum vector has norm squared c/2. The generic Virasoro vacuum descendant counting at fixed level for c>1 does not grow linearly with c. Counting stress tensors therefore counts factors, not the central-charge weight of general factors.

[Verified literature statement] The basic Hellerman modular bound, under its hypotheses of unitarity, modular invariance, c_L,c_R>1 and no extended chiral algebra, is

    Delta_1 < (c_L+c_R)/12 + 0.473695.

It is an upper bound on the first primary dimension given the central charges. Rearranging a supplied lower gap gives a lower constraint on central charge, not the needed upper constraint. See [Hellerman](https://arxiv.org/abs/0902.2790). Numerical modular gap bounds likewise constrain spectra at specified c, as in [Collier, Lin, and Yin](https://arxiv.org/abs/1608.06241).

[Proved modular limitation] Modular invariance of an infrared theory contains no M or lattice support size. Changing an overall velocity does not change its dimensionless modular spectrum. Tensor products remain consistent continuum theories while adding central charge. Sector-by-sector modular bootstrap requires each factor to be a full modular-invariant theory, which is an additional assumption if the sectors are only chiral subalgebras or components subject to selection rules. Even granting that assumption, the microscopic velocity-weighted budget is still missing. Factorization of the infrared theory also does not imply factorization of the on-site Hilbert space, so one cannot allocate independent microscopic dimensions to its factors without another theorem.

[Assumed sufficient operator lemma, inherited structure from C106] Let V=4eM and

    N(u)=sum_(v_a>=u) c_a.

For every 0<u<=V, suppose there is a field space W_u with

    dim W_u >= eta N(u),

and a surjective linear map onto W_u from operators on k(u) sites, where

    k(u) <= k_0(d) + beta(d) log(V/u),
    gamma(d)=2 beta(d) log d < 1.

The first inequality must count central charge, not merely the number of factors. The second is a quantitative support statement for the fast sectors. Neither follows here from modular bootstrap or Lieb-Robinson.

[Proved conditional consequence] Dimension counting gives

    N(u) <= eta^(-1) d^(2k_0) (V/u)^gamma.

Tonelli, or a finite sum over the assumed sectors, gives

    A = integral_0^V N(u) du
      <= 4e M d^(2k_0) / [eta(1-gamma)].                   (18)

This permits unbounded N(u) as u approaches zero and does not impose a forbidden bound on the number of slow sectors. Replacing 4e by the earlier conservative 6e reproduces the C106 constant. The improvement here is only the inherited C109 velocity constant, not a proof of the operator lemma.

[Open] Modular bootstrap has not established either required part of that lemma. The fixed-support operator space and the propagation cone alone do not justify it. Applying a CFT formula at a scale chosen from v or M would also require the crossover control already disproved in C106.

## 9. Literature and limitations of the completed work

[Verified literature context] [Lawrence, Bootstrapping Lattice Vacua](https://arxiv.org/abs/2111.13007) develops positivity-based ground-energy lower bounds, including infinite-volume formulations. [Berenstein, Hulsey, and Lloyd, Numerical exploration of the bootstrap in spin chain models](https://arxiv.org/abs/2406.17844) studies finite-size energy and conformal-data extraction and notes periodic solutions arising in a translation-invariant relaxation. [Berenstein and Lloyd, Successes and challenges of using Semidefinite Programming for the study of Spin Chain Hamiltonians](https://arxiv.org/abs/2606.02013) studies additional constraints and approximate Virasoro data and identifies scaling limitations away from the free case. These works motivate circumference-sensitive constraints; none is being cited here as a universal bound on A.

[Proved work record] Only work_C117 was modified. No package was installed. The Typhon SSH attempt exited 255 with "Connection timed out during banner exchange". Thus no large-block or large-ring SLURM experiment was run. run_typhon.slurm prepares the larger block search and bootstrap.py rejects its extended mode outside a SLURM allocation. The local computations were small validation calculations, with dimensions at most 81, and completed in seconds. The parent STATE.md is outside the authorized write area; STATE.md here supplies the handoff instead.

## 10. Consistency checks

[Proved] Rescaling h by a positive lambda rescales M, velocities, A, F, and energy certificates by lambda. It rescales m_1 and its q^4 coefficient by lambda^3. The dimensions of (6), (11), (17), and (18) agree. At J=0 the Ising ground state is a product, and the massive q^4 contribution in (14) vanishes. The expansion in (14) is used only near J/g=0, not at the critical point J=g.

[Proved] The conditional counting bound integrates an allowed divergent N(u) near u=0, so it does not contradict unboundedly many slow sectors. The spectral filter removes a fixed positive excitation threshold only after q is taken sufficiently small for that Hamiltonian. It makes no uniform crossover claim. The fixed-block obstruction requires A>0 and therefore makes no false claim that product or frustration-free gapped states cannot have exact finite-block energy certificates.

[Proved] The numerical ring bounds are invariant statements about the chosen Hamiltonians, with no CFT fit assumed. The Heisenberg and XXZ exact short-ring saturators independently match the local certificate constants. The SU(3) nonvacuum L=4 ring is excluded from Casimir interpretation. These checks separate finite-volume certification from an unproved infrared extrapolation.
