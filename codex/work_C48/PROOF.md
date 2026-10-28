# C48. Resummed massive matching and its application boundary

## 1. Verdict and provenance

[proved, scope] A resummed cluster expansion gives the massive-matching part of S for an explicitly specified enlarged, rooted cell model at sufficiently low temperature. The bounds below are uniform in the number of cells and in beta above a fixed threshold. They include every external field degree, two time moments, scalar terms and bounded sources. The matching constant and the time-decay rate in units of the complete cell gap can be chosen independently of the wrapper ratios, at each fixed n.

[assumed, remaining application] This is not a proof of the complete physical input S. The construction needs a bounded local interaction representation in the legal rooted sector. Its auxiliary embedding introduces invalid states. Their low-energy threshold and pressure contribution can be bounded, but its thermal source functional is not identical to the original constrained physical trace. Extending the comparison to all physical sectors, fixed total length, the cyclic quotient, filling derivatives and the subsequent multichannel RG remains necessary. No central-charge theorem is asserted.

[proved, attribution] The operator-valued embedding and coarse-time expansion are already in ../../contrib/direct_attack_20261028/infrared/INFINITE_CELL_PROGRESS.md, with an independent audit in its sibling alternatives folder. C48 uses that route, supplies a self-contained conservative version of the cluster estimate, and specializes its dimensionless constants to C42/C44/C46. C42/C44's finite spectral and source identities and C45's finite-order local transformation are imported, not rederived. No primary literature was accessed over the network. The combinatorial convergence argument needed here is given below.

[proved, tagging convention] A tag applies to its following paragraphs and displays until the next tag or section. Mathematical statements about the wrapper are conditional on the imported finite construction wherever that is stated.

## 2. Precise rooted model and finite coefficient conversion

[assumed, theorem hypotheses] A cell is a separable graded Hilbert space K with positive self-adjoint A, a finite code P K isomorphic to m fermion modes, and

    ker A = P K,        A >= Delta Q,       Q=1-P,
    Tr_Q exp(-t A) <= C0 exp(-t Delta),     t >= t0.
                                                        (1)

On a cell ring or interval, suppose the Hamiltonian is A_tot+B+R. Each bounded interval term B_X reduces the local all-code projection P_X. Both interactions are even. For a>=0 use

    ||B||_a = sup_j sum_(X contains j) exp(a |X|) ||B_X||.
                                                        (2)

A local fermionic representation with these interval supports is part of the hypothesis. A dimension count alone does not assert that a particular fermionization preserves the desired Gaussian reference. The fields retained below are those of that representation.

[proved, exact embedding] Let E = C Omega direct_sum Q K, set a_aux Omega=0 and a_aux|QK=A|QK, and define the even isometry

    J |u,code> = |u>_log tensor Omega,
    J |s,exc>  = |0>_log tensor |s>_aux.                 (3)

The grading of |s> is retained. This definition uses a vector-space identification of the code, not a tensor factorization of the original excited space. It intertwines A and I_log tensor a_aux. For each interval define embedded operators with J_X and J_X^dagger. Put

    H_log,X = P_X B_X P_X, extended on logical factors,
    K_X = J_X B_X J_X^dagger - H_log,X tensor I_aux,
    Rtilde_X = J_X R_X J_X^dagger.
                                                        (4)

Then the enlarged Hamiltonian is sum a_aux + H_log + K + Rtilde. Each auxiliary coefficient of K_X kills Omega_X from both sides. Indeed (4) vanishes on every |u> tensor Omega_X, its adjoint also does, and logical matrix units are linearly independent. The spectrum or bandwidth of H_log is irrelevant to this assertion.

[proved, coefficient conversion] On r cells, expand a bounded operator O in logical matrix units. There are 4^(mr) blocks, each with auxiliary operator norm at most ||O||. Every logical matrix unit can be written with creation operators, the vacuum projector product(1-n_i), and annihilation operators. Expansion of that projector uses at most 2^(mr) terms. CAR reduction cannot increase their number or degree, which is at most 2mr. Therefore for h>=1,

    sum_I h^degree(I) ||O_I|| <= (8 h^2)^(mr) ||O||.
                                                        (5)

The estimate is independent of the auxiliary dimension. Odd coefficients are interpreted in the graded tensor product. No auxiliary vacuum projector is expanded in auxiliary fermion modes.

Fix h=4. For an auxiliary-operator-valued external Grassmann interaction V, define

    v_L(V) = sup_j sum_(X contains j) exp(L |X|)
                  sum_I 4^degree(I) ||V_(X,I)||.
                                                        (6)

This is a density at each insertion time, not an integral over beta. With a=L+m log(128)+1, (5) gives

    k := v_L(K) <= 2 ||B||_a,
    r := v_L(Rtilde) <= ||R||_a.                        (7)

Scalar coefficients are included. Time-local coefficient kernels may be measures, including delta measures, rather than functions.

## 3. Quantitative cluster theorem

[proved, statement] Fix an output spatial exponent alpha>=0. Choose

    L = 8(alpha + log(32) + 10),
    T >= max(2L, Delta t0, log C0 + 2),
    tau = T/Delta,       kappa = L/(16 tau),
    v_* = 1/(16 tau),    B_* = 2^16 exp(L)/tau.
                                                        (8)

The logarithm of the full auxiliary trace with all logical fields left external has convergent connected kernels when v_L(V)<=v_*. The vacuum-boundary version has the same bound. For beta>=2tau, their absolute coefficient norm, with field weight 4 and weight

    exp(alpha spatial_diameter + kappa temporal_diameter),
                                                        (9)

is at most B_*, uniformly in the number of cells, beta and auxiliary cutoffs. At finite beta use the circular time metric. For the spatial weight one may use the number of cells in the smallest connected interval covering the cluster support, including intermediate auxiliary sites; on a ring use the shortest covering arc, or the whole ring if necessary. This stronger support weight dominates the external spatial diameter and matches the interval convention in (2). At zero temperature use the line time metric. Scalar coefficients have norm per cell per unit time. For a nonconstant translation-invariant kernel, divide the absolute fully integrated coefficient measure by spacetime volume, or equivalently pin one external coordinate and sum its species. The latter dominates a supremum over its species. This is the convention for the norm below. Additional local support or external-leg marking is paid by the reserves in the proof and, for leg marking, by a fixed loss from field weight 4 to 2.

If k+r<=v_*/4 and K has the two-sided vacuum-annihilation property in (4), then

    ||W_infinity[K+Rtilde]|| <= A0 r,
    ||W_beta[K+Rtilde]|| <= A0 r + B_* exp[-beta/(4tau)],
    A0 = 4 B_*/v_* = 2^22 exp(L).                      (10)

Here W is log of the auxiliary influence functional; changing the action sign leaves every norm unchanged. The free auxiliary scalar may be retained, as in (10), or removed in differences. All estimates include the scalar density.

[proved, block bound] Divide the time circle into N_t=floor(beta/tau) blocks, with common width beta/N_t in [tau,2tau]. Within a block expand by Duhamel's formula using the contraction semigroups exp(-s sum a_aux), s>=0. Let S be the set of distinct active interval labels. Write

    v_X = sum_I 4^degree(I) ||V_(X,I)||,
    z_X = exp(2tau v_X)-1.

Summing all words which use every X in S at least once gives an operator/coefficient norm at most product_(X in S) z_X. To see the factorial explicitly, multiplicities m_X have total p, the simplex volume is at most (2tau)^p/p!, and there are p!/product_X m_X! label words. Summing each m_X>=1 gives z_X. This bound is valid for noncommuting coefficients. It never uses exp(sa_aux) V exp(-sa_aux).

At inactive sites split a full-block free propagator into P_Omega and exp(-tau_t a_aux)Q_Omega. A node of the spacetime square grid is occupied if it lies in an active interval or uses the passive Q term. Connected components, with nearest-neighbor spatial and temporal adjacency, are decorated polymers. The active support fills the whole interval, including sites without a field insertion.

[proved, factorization] A vacant node carries the rank-one P_Omega. Cutting this wire separates the auxiliary contractions. Disjoint spatial supports in a block factor after Duhamel interleavings are summed, by the tensor-product semigroup identity. Every contraction connecting different times on a site crosses consecutive occupied nodes. Hence different grid components factor. Coarse adjacency may merge independent contractions but cannot split a connected contraction.

A completed polymer has even external parity, since an odd auxiliary vacuum matrix element or trace is zero and the total vertex parity is even. Its weight therefore belongs to the commuting even part of the external Grassmann algebra. This yields the exact hard-core gas identity

    Z = sum_(pairwise nonadjacent polymers) product_gamma w_gamma.
                                                        (11)

This identity can first be read on a finite time lattice, with finite spatial and spectral cutoffs.

[proved, thermal trace bound] The heat hypothesis and (8) imply

    Tr_E exp(-tau a_aux) <= 1+C0 exp(-T) < e.

For one completed component with spatial support S, insert P_Omega at its absent nodes. Split every free semigroup into two equal halves. At a passive Q site the first halves over the whole block commute with every active insertion in that block and cost exp(-Delta tau_t/2). Keep the second halves of all segments. If their original durations are s_i with sum_i s_i=beta, Schatten Holder with exponents beta/s_i gives

    |Tr product_i [exp(-s_i A_S/2) O_i]|
       <= Tr exp(-beta A_S/2) product_i ||O_i||
       <= exp(|S|) product_i ||O_i||.                  (12)

All remaining first halves and projectors are contractions. This proves a passive cost at most exp(-T/2)<=exp(-L) and an additional factor at most exp(N) for a polymer with N nodes. The calculation is dimension independent and includes winding components. Vacuum boundary conditions need no trace factor and give a stronger bound.

[proved, summing arbitrary interval decorations] Set epsilon=exp(-L) and zhat_X=exp(L|X|) z_X. Since 2tau v_* =1/8,

    eta := sup_j sum_(X contains j) zhat_X
       <= 2tau exp(2tau v_*) v_* = exp(1/8)/8 < 1.
                                                        (13)

For a specified occupied support Y of N nodes, the passive set in each layer is determined by its active intervals. Their total support cardinality is at least the size of their union. Thus each decoration, before the trace factor, is bounded by epsilon^N product zhat_X. Summing all active sets contained in Y, even ones not producing exactly Y, gives

    sum_(decorations of Y) ||w||
       <= epsilon^N exp(eta N) exp(N)
       <= exp[-(L-2)N].                                (14)

This sum uses the activities, not a count of all interval subsets. In particular no finite maximum interval length, or factor growing as 2^(N^2), is needed.

[proved, connected graph bound] A connected square-grid support of N nodes containing a specified node has at most 16^N possibilities: fix an ordering of the four neighbors, choose the corresponding deterministic spanning tree, and encode its depth-first traversal by at most 2N steps. The same upper bound holds on finite rings and time circles, whose degree is at most four.

For completeness, the hard-core logarithm has coefficients

    phi(gamma_1,...,gamma_p)
      = sum_(connected graphs G on p vertices)
             product_({i,j} in E(G)) [-1_(gamma_i incompatible gamma_j)],

and factorial 1/p!. Its absolute value is bounded by the sum over spanning trees of incompatibility indicators. One proof partitions connected subgraphs of a fixed incompatibility graph into intervals [tree, closure(tree)]. Root at vertex 1, assign breadth-first distances, choose the lowest-labeled parent in the preceding level, and put in the closure all allowed same-level edges and preceding-level edges whose lower-level label is no smaller than the chosen parent. These are exactly the edges that preserve the breadth-first tree and parent choices. The signed sum on an interval is (-1)^(p-1)(1-1)^k, with k the number of optional edges. Its magnitude is at most one. Summing over trees proves the bound.

The rooted-tree recursion is also explicit. For nonnegative activities z_gamma, let F_gamma be the sum of trees attached to a prescribed root, with its own activity excluded. Truncating tree depth and attaching an unordered collection of children gives

    F_gamma <= exp[sum_(eta incompatible gamma) z_eta F_eta].

If sum_(eta incompatible gamma) z_eta exp(|eta|)<=|gamma|, induction in depth proves F_gamma<=exp(|gamma|). Monotone convergence then proves absolute convergence of the rooted connected expansion. This is the only cluster convergence criterion used here.

[proved, weights and reserves] A connected cluster of total N occupied nodes, including multiplicities, has spatial diameter at most N and time diameter at most 2tau N+4tau. The extra 4tau safely accommodates endpoints and a shifted block origin. Its weight (9) is bounded by

    exp(L/4) exp[(alpha+L/8)N].

For the convergence test reserve exp(4N), one factor for the tree recursion and three for roots, tails and boundary losses. Equations (14) and the animal bound give the per-node sum

    sum_(N>=1) q^N,    q=16 exp[-7L/8 + alpha + 6] < 1/100.
                                                        (15)

An N-node polymer has at most 5N possible incompatible anchor nodes. Hence its incompatible sum is at most 5N q/(1-q)<N, which proves the recursion criterion with substantial slack. A factor N for marking a support site, or N^2 for a second marking, is absorbed by the unused exponential reserve. More explicitly, the rooted-tree bound gives at most q/(1-q)<1/99 for clusters attached to a specified grid node after the displayed reserves. Counting a cluster once per possible root only overcounts its contribution to the density. The endpoint factor is exp(L/4), and division by the block width costs at most 1/tau. Support marking consumes a reserve because N<=exp(N) and N^2<=2exp(N). Boundary-layer shifts by at most two layers cost at most exp(4) when using the remaining tail reserve. Even allowing a factor 64 for these markings and the two line/circle tails gives 64 exp(L/4+4)/tau < B_*. Thus the stated B_* bounds the density, scalar terms and boundary tails. No positivity of quantum activities is required.

[proved, continuous fields and kernels] Normal-order each finite logical polynomial before replacing it by external fields. At a vertex put its fields at ordered split times inside an interval of length epsilon_time, with later creation times. The length of that interval is independent of the number of fields. Use delta measures or regularized deltas of total variation one. Fixed-order Duhamel coefficients converge as the split goes to zero. Bounds (13)-(15) are uniform under sufficiently small splits and dominate the coefficient measures, so the connected expansion converges in the stated norm.

Equivalently one can start with coherent-state resolutions on a finite logical time lattice, keeping auxiliary operators unexpanded. The logical trace gives antiperiodic boundary conditions. The graded tensor product supplies the signs of odd auxiliary coefficients. The block bounds control the time-lattice limit. Ordinary cyclic permutations of odd Grassmann-valued auxiliary operators are not assumed.

Translation invariance in continuous time follows from partition independence of the exact functional. For vacuum slabs use arbitrary block origins and partial first/last blocks. A passive Q term touching a vacuum boundary vanishes; active estimates only improve in a shorter block. Clusters reaching a boundary have exponentially suppressed weight by (15). Thus interior kernels have a common, continuously translation-invariant limit. Their total absolute integral divided by spacetime volume equals the relative-coordinate integral with one external coordinate fixed, by a change of variables of unit Jacobian. Scalar terms have no external coordinate and retain the density convention. This also explains why a general inhomogeneous pinned theorem would require a marked-source formulation.

[proved, analyticity and the leakage factor] The block activities are entire in bounded vertices and their convergent majorant is uniform on v_L<=v_*. Hence W is analytic in that ball as a map to the weighted kernel space. Along K+s Rtilde, 0<=s<=1, the complex disk of radius v_*/(4r) in the direction Rtilde stays inside v_L<=v_*/2. Cauchy's estimate bounds its derivative by 4B_* r/v_*. Integrating over s gives

    ||W[K+Rtilde]-W[K]|| <= A0 r.                       (16)

For r=0 this is immediate. At zero temperature every K-only vacuum word is zero: its first insertion on Omega is zero. Thus W_infinity[K]=0, including its scalar density.

At finite temperature any nonzero K-only polymer must meet every time layer. Otherwise cut at an empty layer, obtain a vacuum amplitude, and use the same cancellation. Such a polymer has at least N_t nodes. The spare exponential reserve in (15) bounds all such clusters, including their external coefficient sums, by B_* exp[-beta/(4tau)], since N_t>=beta/(2tau). This proves (10). An exponentially small pressure alone would not justify this coefficient bound; the layer argument does.

[proved, temperature-independent endpoint] The same reserve compares W_beta[V] with the periodization of W_infinity[V]. A cluster whose occupied layers leave a cut can be lifted to the line, with the fermionic antiperiodic signs. Lifts are counted modulo common translation by beta. Clusters that wrap or whose line lifts span a full period have at least N_t-O(1) nodes. Their weighted sum is exponentially suppressed by the remaining reserve. Enlarging the numerical factor within B_* if needed, the difference is bounded by B_* exp[-beta/(4tau)]. Periodization does not increase the norm: its absolute integral unfolds into the line integral and circular distances do not exceed distances of a lift. This statement concerns the enlarged influence action.

## 4. Time moments, frequencies and sources

[proved] Give the output fields weight 2 instead of 4. For p external fields and temporal diameter d, the sum of distances to an anchored field is at most p d. Using p^2/2^p<=4 and x^2<=2 exp(x),

    2^p (1+Delta sum_i |t_i-t_anchor|)^2
      <= [2+16(Delta/kappa)^2] 4^p exp(kappa d).
                                                        (17)

The same estimate holds on the circle with circular distances. Reserving a spatial exponent supplies spatial moments in cell coordinates. For frequencies with max_i |omega_i|<=omega,

    |exp(i sum_i omega_i(t_i-t_anchor))-1|
      <= omega p d.

Since p/2^p<=1 and d exp(-kappa d)<=1/kappa, the induced frequency difference is at most A0 r omega/kappa. A single explicit constant sufficient for both conclusions is

    A_dec = 2^22 exp(L) [2+16(16T/L)^2].                 (18)

Consequently the zero-temperature action satisfies the S(a)-type estimates

    ||delta S_R||_LE <= A_dec ||R||_a,
    ||delta S_R(omega)-delta S_R(0)||_LE
        <= A_dec ||R||_a omega/Delta.                  (19)

These are bounds on the massive influence functional before integrating any retained propagator. There is no factor 1/W.

[proved, bounded sources] Add a bounded local source as another coefficient vertex and differentiate within a smaller analytic ball. Cauchy's bound controls any fixed number of such derivatives and preserves exponential diameter weights. For a code-projected observable, embed the actual projector-containing operator, not a tensor logical observable silently extended through the excited sector. The source-free scalar log functional must be kept. Neither this construction nor (19) computes its later gapless T^2 coefficient.

[assumed, unbounded sources] The physical cell length is unbounded away from the code and does not belong to (6). Its chemical-potential derivatives and the fixed-length projection require a relative-bound or ensemble argument. They are not consequences of the bounded-source statement.

## 5. Substitution of the C42/C44 constants

[assumed, imported finite construction] Fix n, the accepted C42 certificate and C44 wrapper parameters. Write b=n+3, k=n+2, d0=4464, and m=n+1 retained modes. All these data are fixed before choosing lambda. This section applies to C45's well-formed initialized rooted sector, where the readout has norm at most one because there is exactly one token. It does not extend that boundedness statement to arbitrary illegal strings.

[proved conditional on those inputs, complete gap] The following is an explicit lower-gap prescription using C44 (3)-(5). Use its notation rho=w_vac/U and impose zeta=rho(k+1)/(1-rho(k+1))<=1/2. Put

    gamma_c=min(K_c,1-|g|,1+q^2-2q cos(pi/k)),
    h_*=K_c+D+12,       a_*=K_c+D+5,
    M_c=(k-1)h_*+2a_*,
    M_occ=s M_c+J_f(k+1),
    gamma_occ=min(s gamma_c,J_f/[k(2k-1)]),
    d_*=1-zeta^2.

For p=1,2 let N_p=k+p, m_p=binomial(k+p,p), and

    delta_ch,p = w_vac [rho/(N_p-1)]^(2kp)/(m_p kp),
    delta_occ,p = delta_ch,p d_* gamma_occ /
        [(sqrt(delta_ch,p)+sqrt(M_occ))^2+d_* gamma_occ].

Here the right sides are lower bounds, not the actual eigenvalues. Take Delta positive and no larger than the minimum of delta_occ,1, delta_occ,2 and

    U[1-rho(k+1)], J_h/[(b+1)(2b+1)],
    2t_sel[cos(pi/b)-cos(2pi/b)],
    2t_sel[cos(pi/(b+1))-cos(2pi/(b+1))],
    gamma=t_sel eta_b-U, U, |e|b(b+1).
                                                        (20)

The length-one penalty is also included if smaller. The positive eta_b, e and filter parameters are those in C44 Section 3. In particular Delta<=gamma. Replacing this Delta by U or by the bare core gap would invalidate the estimate.

[proved conditional on the legal counting input, heat bound] C46 Section 2 supplies the following bound for the same legal rooted space. Exceptional lengths 1,b,b+1 have total dimension at most

    D_fin=d0+d0^b+d0^(b+1).

The shorter nonselected lengths contribute at most b d0^(b-1) exp(-beta Delta). For lengths b+2+r, r>=0, C44's linear cost gives the geometric sum

    sum_(r>=0) d0^(b+2+r) exp[-beta gamma(3+r)]
      = d0^(b+2) exp(-3beta gamma)/(1-d0 exp(-beta gamma)).

For beta>=t0=log(2d0)/gamma this is at most 2 d0^(b+2) exp(-beta Delta). Hence (1) holds with

    C0 = D_fin + b d0^(b-1) + 2 d0^(b+2).
                                                        (21)

This is a trace over the entire unbounded legal length space, not a fixed length cutoff. It uses the stated state-count upper bound d0^ell and the legal wrong-length energy bound. It is not a heat-trace proof for all illegal selector words.

[proved, dimensionless simplification] Since Delta<=gamma,

    Delta t0 <= log(2d0).

Choose once and for all

    T_n = max(2L, log(2d0), log C0+2).                   (22)

Equations (8)-(19) then hold with

    beta >= 2T_n/Delta,
    kappa = [L/(16T_n)] Delta,
    A_dec,n = 2^22 exp(L) [2+16(16T_n/L)^2].             (23)

At fixed alpha and n these constants depend on neither the core scale nor rho, g, D, the charge tunneling splitting or the logical bandwidth. The input support exponent a=L+(n+1)log128+1 also has that property. This improves the bookkeeping in the existing abstract theorem: taking its heat time as an unrelated parameter would obscure the independence from wrapper ratios. Dependence on n is allowed and can be large.

[proved conditional on C45, explicit final coupling] Recompute the C45 finite-order transformation at the strengthened exponent a. Set

    v0=||V||_a, f0=||F||_a,
    G=f0+48v0^2/Delta,
    C_SW=2^16 (v0 G/Delta+v0^3/Delta^2).                 (24)

C45 gives B=lambda V^d+lambda^2 G2^d, and

    ||B||_a <= 2lambda v0+2lambda^2 G,
    ||R||_a <= C_SW lambda^3,

provided lambda<=min(v0/G,Delta/(192v0)). Combined with (7),

    k <= 4lambda v0+4lambda^2 G <= 8lambda v0,
    r <= C_SW lambda^3.                                (25)

Thus the following sufficient, explicit choices imply k+r<=v_*/4:

    lambda <= min(v0/G, Delta/(192v0),
                  Delta/(1024 T_n v0),
                  [Delta/(128 T_n C_SW)]^(1/3)).        (26)

The harmless case of zero v0 is treated by choosing a positive upper bound for it. All norms here are full physical cell interaction norms, including compensators, marker overlap inverses, readout and scalar adjustments. A symbolic specification of these norms is a finite bound from the given local operators, not a numerical certificate at the extreme C42 witness values.

The retained payload scale from C42/C45 is

    J0=2 g^2 q^(2n+2) D^(-2n-4)/S_q^2,
    S_q=sum_(ell=0)^(n+1) q^(2ell),
    c_n=J0/s,             W=c_n lambda^2.               (27)

Consequently the massive matching error relative to W obeys

    ||delta S_R||_LE/W <= A_dec,n C_SW lambda/c_n.       (28)

For any prescribed epsilon>0, also require

    lambda <= epsilon c_n/(A_dec,n C_SW).               (29)

This gives a fixed positive lambda independent of cell count. There is no division by a global many-body spectral gap or by an individual logical velocity. Slow velocities enter the later logical tolerance and RG, which are separate.

[proved, static errors stay separate] The complete second-order H_log is retained, including scalar and length-dependent terms. C42 bounds its conversion remainder in J0 units by

    (2+64n^2) 3(1+eta) g^2 D^(2n+4).

C44 bounds the four unequal length-pair bonds, including their scalar table. These must already be made small by the earlier D,g,rho,s choices. Decreasing lambda in (29) cannot remove an error proportional to W. Conversely, the constants in (23) are fixed before those choices, so this decoupling estimate introduces no new circular dependence into their tolerance. This does not prove that the unknown RG constants have the same independence.

[proved, temperature qualification] Formula (10) supplies a thermal baseline of size B_* exp[-beta Delta/(4T_n)]. To make the full initial action, rather than its change under R, epsilon-close in W units, beta must additionally be large enough that this term is <=epsilon W, with the moment factor from (17) included. For each fixed lambda this is a finite low-temperature threshold. It is not an estimate vanishing with lambda uniformly at a fixed beta. The beta-independent endpoint is W_infinity; thermal matching approaches it exponentially.

## 6. The C46 counterexample and what the resummation changes

[proved, obstruction retained] C46 considers

    A=q_aux(Delta+u N),        N=n1 n2,
    V=lambda sigma_x,         u>=Delta>0.

In the original tensor logical coordinates the order-lambda^2 raw auxiliary excursion has absolute coefficient sum

    (lambda^2/Delta) sum_(j>=1) (h^4 u/Delta)^j,

which diverges for h>=1. Reordering an absolutely divergent nonnegative majorant cannot cure it. A Mayer tree inequality applied after that expansion also cannot cure it. Our theorem imposes smallness on (6); putting u q_aux N into V would violate that condition, so the theorem makes no false assertion about this raw influence functional.

[proved, explicit spectral resummation] Instead keep the whole A in the cell spectral problem. Its code is the four states with q_aux=0. Introduce four orthogonal auxiliary excited states e_s labelled by s=00,10,01,11, of energies

    E_s = Delta, Delta, Delta, Delta+u.

Map |s,0> to |s> tensor Omega and |s,1> to |0> tensor e_s. The gap is Delta and Tr_Q exp(-tA)<=4 exp(-tDelta), independently of u. The flip becomes

    Vtilde/lambda = sum_s [ |0><s|_log tensor |e_s><Omega|
                         + |s><0|_log tensor |Omega><e_s| ].
                                                        (30)

Its operator norm is one and (5) bounds its coefficient norm by (8h^2)^2. Crucially u is now in the auxiliary semigroup, not an external quartic vertex. At order lambda^2 there are only finitely many external monomials, at the two endpoints; their duration factors are exp(-E_s t). Their coefficient norm is bounded by a fixed multiple of lambda^2/(Delta-kappa). No series in u/Delta occurs. The all-order cluster theorem converges for sufficiently small lambda with constants independent of u.

[proved, operator consistency check] Integrating the operator excursion gives

    -lambda^2 sum_s |s><s|/E_s
      = -lambda^2/Delta
         +lambda^2 u/[Delta(Delta+u)] N.                (31)

This agrees with C46. Equation (31) is an operator identity. One must not replace two endpoint Grassmann polynomials by their coincident product and use N(t)^2=N(t), which is false in the external Grassmann algebra. The isometry changes how retained variables extend off the code; logical contact contractions and source matching change with it. For the wrapper, performing C45's local transformation first fixes the complete static second-order H_log and obtains (28) from the exact cancellation W_infinity[K]=0. It avoids an unsupported identification of a raw retarded vertex with its static operator product.

[proved, source distinction] The raw tensor observable n_i has auxiliary-sector-dependent representation under (30). Code-projected P n_i P instead becomes n_i tensor |Omega><Omega| on the valid subspace. These are different sources. The resummed theorem bounds the actual embedded sources. It does not claim equality of raw action coefficients in the original off-code coordinates, so C46's counterexample remains valid.

## 7. Original versus enlarged model and the unsolved parts of S

[proved, exact low-energy comparison] In the enlarged model the valid projection JJ^dagger on every cell is conserved. Embedded interactions touching an invalid cell vanish. Let

    w=sup_j sum_(X contains j) ||(B+R)_X|| < Delta.

In an invariant sector with invalid-site set I, cut all interactions touching I. Invalid sites cost at least Delta each. The remaining Hamiltonian on the other sites has some ground energy E_rest. An all-valid trial state consisting of that ground state and code vectors at the cut sites has energy at most E_rest+w|I| after restoring the cut bonds. Therefore

    E_(invalid I) >= E_original + (Delta-w)|I|.         (32)

The spectra of original and enlarged models below E_original+Delta-w coincide, including multiplicities. This includes energies of order 1/M for sufficiently large M, without proving that such energies have a CFT organization. Conditions (26), together with a stricter lambda if needed, ensure w<Delta/2. This comparison is between the original rooted tensor model and its embedding, not between different physical marker counts.

[proved, pressure comparison] The one-cell invalid and valid heat traces are

    z_inv(beta)=(2^m-1) Tr_Q exp(-beta A),
    z_val(beta)=2^m+Tr_Q exp(-beta A).

In an invalid sector, compare its cut-system trace with the all-valid trace after restoring bonds. The change in log partition function is at most beta w |I|, by the variational eigenvalue bound for a bounded perturbation. Summing invalid subsets gives

    0 <= f_original-f_enlarged
       <= beta^(-1) log[1+exp(beta w) z_inv/z_val]
       <= beta^(-1) log[1+C0 exp(-beta(Delta-w))].       (33)

This is uniform in M at beta>=t0. It is exponentially small at fixed couplings and hence smaller than every thermal power. An analogous exact identity for arbitrary thermal source functionals has not been proved by (33).

[proved, why all beta needs qualification] A gap and (1) need not imply a trace for every beta>0. For example, a direct sum with 2^ell states of energy Delta+gamma ell has a finite heat trace only for beta gamma>log2. Thus no theorem for a grand-canonical infinite cell space can remove a heat-time restriction using those hypotheses alone. A fixed physical-length spin chain has a finite trace at every positive beta, but obtaining it by a total-length projection of this low-temperature rooted expansion is a different argument. Our uniformity means all M and all beta>=2T_n/Delta, which includes the infrared limit at fixed microscopic coefficients.

[assumed, unresolved physical application] The missing obligations are explicit.

1. Verify the bounded interval coefficient representation, with the required logical operator lattice and local fermionization, for the complete rooted physical sector. C45 proves bounded readout in initialized well-formed cells. Positive readouts can be unbounded on illegal strings and cannot be placed in (6) without a separate argument.
2. Transfer the retarded action and its sources to fixed physical length and the cyclic root quotient, including the unbounded length source and smooth filling response. C46's quotient result treats the common-bond reference, not all the resummed interactions.
3. Establish dressed exclusion for the remaining physical sectors under the full coupling. The invalid states in (32) are artifacts of the embedding and are not physical illegal cells, so (32) does not settle this question.
4. Prove the gapless multichannel RG for the arbitrary-degree retarded kernels obtained here, with source and scalar control. The cluster expansion integrates no retained line and therefore cannot establish that RG or identify n+1 currents and the stress tensor.

[proved, final theorem boundary] The proof establishes an all-degree, exponentially weighted, leakage-linear massive matching theorem for the stated enlarged rooted model. It also supplies a uniform low-energy and pressure comparison with the valid rooted model, and an explicit C42/C44 conditional schedule. It does not establish full S, a microscopic central charge, or an unconditional counterexample to the c-d conjecture. Every claim of uniformity is at fixed n; neither the constants nor the small-coupling threshold are claimed uniform as n grows.
