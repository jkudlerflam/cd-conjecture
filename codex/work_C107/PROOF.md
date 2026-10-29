# C107. Compactness and semicontinuity of the Casimir coefficient

2026-10-29. Tags apply to the statements and derivations following them. [Proved] denotes a mathematical argument given here, [assumed] an additional physical or mathematical hypothesis, [numerical] a finite computation, and [heuristic] an interpretation without a proof. [Open] marks an unresolved statement. This investigation has not had an independent referee review.

## 1. The compact parameter space

[Proved] Let V be the real vector space of Hermitian operators on C^d tensor C^d, and set

    G = {alpha I + o tensor I - I tensor o : alpha real, o Hermitian}.

The scale in the task is the quotient norm

    M([h]) = inf_(g in G) ||h+g||_op.                         (1)

Minimizing over alpha centers the extreme eigenvalues, giving half their difference. The subspace G is closed, so the infimum is attained. It vanishes precisely on G. Positive scaling multiplies M and the physical A by the same factor. Negative scaling is not an equivalence of ground-state problems.

The literal set {h in V : M(h)=1} is not compact, since h+nI remains in it and is unbounded. The correct compact set is

    K_d = {[h] in V/G : M([h])=1}.                           (2)

For d>=2 the quotient has real dimension d^4-d^2. It can be represented in the Hilbert-Schmidt orthogonal complement G^perp. If P is that orthogonal projection and h_min is an optimal representative, then

    ||P h||_op <= ||P h_min||_HS <= ||h_min||_HS <= d M(h).

Thus this gauge-fixed sphere is closed and bounded in finite dimension. For d=1 the quotient is zero and there is no nonzero normalized bond.

On a periodic ring the coboundary telescopes, while alpha shifts E_L by L alpha and e by alpha. Consequently

    Q_L(h) = (6/pi) L [L e(h)-E_L(h)]                       (3)

is gauge invariant, as are its limits when defined. Here E_L is the full spin-ring ground energy, with no sector imposed unless stated.

## 2. What finite-volume continuity proves

[Proved] If delta=||h-g||_op, the variational principle gives

    |E_L(h)-E_L(g)| <= L delta,
    |e(h)-e(g)| <= delta,
    |Q_L(h)-Q_L(g)| <= (12/pi) L^2 delta.                  (4)

The second inequality follows by taking the thermodynamic limit of E_L/L. Existence of this limit is also obtained below. The last inequality remains true with delta=M([h]-[g]), by choosing an optimal gauge for the difference and canceling its scalar shift in (3). Each Q_L is continuous on K_d. The estimate is not equicontinuity.

[Proved] Write B_L for the open-chain ground energy in a representative with ||h||<=M. Joining two open chains adds one bond, so

    |B_(n+m)-B_n-B_m| <= M.

Thus B_L+M is subadditive and B_L-M is superadditive. Fekete's lemma gives a common limiting energy density e, since the two sequences differ by 2M, and

    |B_L-L e| <= M.

Closing a chain adds another bond. Therefore

    |E_L-L e| <= 2M,
    |Q_L| <= (12/pi) M L.                                  (5)

This recovers the established C98 estimate. The energy density is controlled uniformly at order 1/L, while the Casimir limit asks for its order 1/L^2 term.

[Proved] Neither a pointwise limit of continuous functions nor their limsup or liminf is automatically upper semicontinuous. In particular,

    limsup_L Q_L = inf_N sup_(L>=N) Q_L,
    liminf_L Q_L = sup_N inf_(L>=N) Q_L.                    (6)

A supremum of continuous functions is lower semicontinuous. An infimum is upper semicontinuous. The outer operations in (6) do not preserve the property needed here. Section 6 supplies a compact parameter countermodel where both limits exist, are finite everywhere, and are unbounded.

If continuous R_L decrease with L on a common domain and converge to A, then A=inf_L R_L is upper semicontinuous. Eventual monotonicity with a starting length depending on the point is not sufficient for this argument on a neighborhood.

## 3. The raw full-length limit is not always the physical A

[Proved] Consider the ferromagnetic XX spin chain

    H = -t sum_j (sigma_j^+ sigma_(j+1)^-
                         + sigma_j^- sigma_(j+1)^+), t>0.   (7)

The bond eigenvalues are -t,0,0,t. Its swap symmetry prevents any coboundary from decreasing its half-spread. Swapping reverses a coboundary, and averaging the two representatives returns the original bond. Convexity and unitary invariance of the half-spread give M=t.

In the N-particle sector the Jordan-Wigner boundary condition is periodic for odd N and antiperiodic for even N. In either case the N occupied momenta can be centered at zero. A geometric sum gives

    E_L(N) = -2t sin(pi N/L)/sin(pi/L).

Minimizing takes N nearest L/2, so

    E_L = -2t csc(pi/L),          L even,
    E_L = -t csc(pi/(2L)),        L odd,
    e = -2t/pi.                                           (8)

Expanding csc x=1/x+x/6+O(x^3),

    lim_(L even) Q_L = 2t,
    lim_(L odd) Q_L = t/2.                                 (9)

The two Fermi velocities have magnitude 2t. The massless Dirac theory has c=1, so physical A=2t. Here limsup Q_L recovers A, the ordinary full-length limit fails to exist, and liminf Q_L=A/4. This issue occurs at M=1 and d=2.

[Proved] A gapped example is h=Z tensor Z, with M=1. Its ring energies are -L for even L and -L+2 for odd L, while e=-1. Thus Q_L=0 on even lengths and Q_L=-12L/pi on odd lengths. Its raw liminf is minus infinity, not the stipulated zero. The frustrated-ring cost is not an infrared central charge.

[Assumed] Identifying A=sum_a c_a v_a with a Casimir limit requires a vacuum-compatible sequence or a justified prescription for sector and commensurability effects. Let C denote the subset where a physical A is defined and finite. Equation (9) does not establish that limsup always supplies the required prescription. A subsequence chosen separately for each bond also does not supply a common continuous sequence on parameter space.

These sector issues agree with the periodic-chain discussion in Nakayama and Suzuki, [arXiv:2601.18867, Sections II.3 and III.3](https://arxiv.org/html/2601.18867v1). Equations (8)-(9) were derived directly in the stated normalization.

## 4. A physical obstruction to uniform convergence

[Proved] In Pauli normalization take

    H_gamma(g) = -sum_j [(1+gamma) X_j X_(j+1)/2
                        +(1-gamma) Y_j Y_(j+1)/2 + g Z_j]. (10)

Put half the field at each endpoint of the bond. Its eigenvalues are +/-sqrt(g^2+gamma^2) and +/-1. The swap-averaging argument of Section 3 gives

    M(h_gamma(g)) = max(1,sqrt(g^2+gamma^2)).                (11)

The Jordan-Wigner quasiparticle dispersion is

    epsilon(k) = 2 sqrt((g-cos k)^2+gamma^2 sin^2 k).        (12)

For gamma=1 and g>1 the dispersion is gapped and analytic in a nonzero complex strip about the real axis. The ground energy is the antiperiodic sum of -epsilon/2. Its sum-minus-integral difference is exponentially small at fixed g>1, by the Fourier coefficient estimate for an analytic periodic function. Hence Q_L tends to zero. At g=1 the cusp gives A=1, as calculated in Section 5.

Therefore A is discontinuous at g=1, in a direction compatible with upper semicontinuity. Uniform convergence of continuous Q_L on a neighborhood containing g=1 and g>1 would imply continuity. Such convergence is impossible.

An explicit diagonal sequence is g_L=1+L^-3. The bond difference from g=1 has operator norm L^-3. Equation (4) gives

    |Q_L(h_1(g_L))-Q_L(h_1(1))| <= 12/(pi L).

Thus Q_L(h_1(g_L)) tends to 1 although A(h_1(g_L))=0 for every L. After division by the continuous norm (11), the same obstruction holds on K_2 with nonzero limit 1/sqrt(2).

This disproves uniform two-sided convergence. It does not disprove upper semicontinuity or exclude a one-sided upper approximation to A by suitable finite-scale quantities.

## 5. No universal monotonicity for the raw Casimir quantity

[Proved] Fix g=1 and 0<gamma<=1 in (10), and write

    f_gamma(k) = sqrt((1-cos k)^2+gamma^2 sin^2 k).

For even L, the physical antiperiodic even-parity vacuum has energy

    E_NS(L) = -2 sum_(n=0)^(L/2-1) f_gamma((2n+1)pi/L),
    e = -(1/pi) integral_0^pi f_gamma(k) dk.                (13)

In the periodic sector the zero mode at k=0 allows the parity constraint at zero cost. Euler-Maclaurin gives E_R=Le+pi gamma/(3L)+O(L^-3). Equation (15) below gives E_NS=Le-pi gamma/(6L)+O(L^-3). Therefore E_R-E_NS=pi gamma/(2L)+O(L^-3)>0 at every fixed gamma>0 for all sufficiently large even L. The asymptotic NS expression is the full spin ground energy.

Expanding at zero gives

    f_gamma(k) = gamma k + b_gamma k^3 + O(k^5),
    b_gamma = 1/(8gamma)-gamma/6
            = (3-4gamma^2)/(24gamma).                     (14)

Odd derivatives vanish at pi, since f_gamma(pi+q) is even in q. With a=2pi/L, the midpoint Euler-Maclaurin formula reads

    sum_n f((n+1/2)a)
       = (1/a) integral_0^pi f
         -(a/24)[f'(pi)-f'(0)]
         +(7a^3/5760)[f'''(pi)-f'''(0)] + O(a^5).

Here f'(0)=gamma and f'''(0)=6b_gamma. Hence

    E_L = L e - pi gamma/(6L)
                  + 7 pi^3 b_gamma/(60L^3) + O(L^-5),     (15)

    Q_L = gamma + a_gamma/L^2 + O(L^-4),
    a_gamma = 7 pi^2 (4gamma^2-3)/(240gamma).               (16)

The remainder constants may depend on gamma. This is a fixed-positive-gamma asymptotic expansion.

For gamma=1, a_gamma>0, and Q_L approaches A=1 from above. For gamma=1/2, a_gamma<0, and Q_L approaches A=1/2 from below. In particular,

    Q_(2L)-Q_L = -3a_gamma/(4L^2)+O(L^-4)

has opposite signs. Neither universal decreasing nor universal increasing monotonicity holds, even along common even lengths. The full Euler-Maclaurin expansion also determines the corresponding eventual sign for consecutive even lengths.

[Proved] No finite constant C makes Q_L+C M/L^2 eventually decreasing for every bond. At g=1, M_gamma=sqrt(1+gamma^2), so

    Q_L+C M_gamma/L^2
       = gamma+(a_gamma+C M_gamma)/L^2+O(L^-4).

As gamma tends to zero, a_gamma tends to minus infinity while M_gamma tends to one. Choose a fixed positive gamma with a_gamma+C M_gamma<0. Then the corrected sequence increases under all sufficiently large doublings. This remains a counterexample if the monotonicity onset is allowed to depend on the bond. A correction o(L^-2) at each fixed gamma also cannot repair the increasing example.

This does not exclude slower-decaying corrections, different finite-scale observables, or model-dependent coefficients diverging near gamma=0. Those possibilities require separate estimates.

[Proved checks] At gamma=1, f_1(k)=2 sin(k/2), and (13) gives E_L=-2 csc(pi/(2L)). Expanding,

    E_L = -4L/pi-pi/(6L)-7pi^3/(1440L^3)+O(L^-5),

in agreement with (15). At gamma=0 the polarized ground state has E_L=-L and e=-1, so Q_L=0. Physical A=gamma tends to zero at this nonconformal quadratic endpoint. The divergence of a_gamma is a nonuniform crossover, not a divergent A.

The sign-changing correction agrees with the irrelevant-field factor 3/4-gamma^2 in Campostrini, Pelissetto, and Vicari, [arXiv:1401.0788, Section III.1 and Eqs. (53)-(56)](https://arxiv.org/html/1401.0788). The normalization and coefficients used here follow from (13)-(16).

[Numerical] `checks.py` compares (13) with nine spin matrices of dimension at most 256, with maximum energy discrepancy 1.95e-14. Quadrature and finite sums verify the sign and coefficient in (16) at gamma=0.1,0.25,0.5,1 and L=64,128,256. CHECKS.md records the outputs. These small computations check the analytic mechanism, not a universal bound.

## 6. A compact countermodel to the elementary argument

[Proved logical counterexample] Let t range over [0,1], and for L>=1 define

    e(t) = -3t^2,
    d_L(t) = t/(1+Lt^2),
    E_L(t) = L e(t)-d_L(t).                                (17)

This is an abstract energy family, not the spectrum of a fixed-d nearest-neighbor bond. It isolates the insufficiency of the elementary hypotheses.

Each E_L is real analytic near [0,1], and

    E_L'(t) = -6Lt-(1-Lt^2)/(1+Lt^2)^2,
    E_L''(t) = -6L+2Lt(3-Lt^2)/(1+Lt^2)^3.                (18)

The positive part of the second term in E_L'' is at most 6Lt<=6L. Thus E_L is concave. Also |E_L'|<=6L+1<=7L, and |e'|<=6. This gives parameter concavity and the order-L perturbation estimate.

The boundary error obeys 0<=Le-E_L=d_L<=1. Its maximum over [0,1] is actually 1/(2sqrt(L)), attained at t=1/sqrt(L). Thus E_L/L approaches e uniformly faster than the basic boundary estimate.

For real L>0,

    partial_L^2 E_L(t) = -2t^5/(1+Lt^2)^3 <=0.

The energy is concave in length, its density increases with length, and its doubling defect is positive. Because d_L decreases with L,

    0<=E_(n+m)-E_n-E_m
       =d_n+d_m-d_(n+m)<=2.                               (19)

An auxiliary open sequence B_L=E_L+1 satisfies

    |B_(n+m)-B_n-B_m|<=1,
    |E_L-B_L|=1.

Thus it obeys the same almost-additivity and closure estimates used in (5), with boundary scale one. In particular B_L+1 is subadditive and B_L-1 is superadditive.

Nevertheless,

    Q_L(t) = (6/pi) Lt/(1+Lt^2),
    A_model(t) = lim_L Q_L(t) = 6/(pi t),    t>0,
    A_model(0) = 0.                                        (20)

Both limsup and liminf equal these ordinary limits. The function is finite everywhere, unbounded on a compact interval, and not upper semicontinuous at zero. Taking t_k=1/k gives A_model(t_k)=6k/pi while E_L(t_k) tends to E_L(0)=0 for every fixed L.

The crossover length is of order t^-2. Beyond it,

    Le-E_L = 1/(Lt)+O(1/(L^2 t^3)),

whereas at fixed L this defect tends to zero as t tends to zero.

[Proved scope] Continuity, finite-L analyticity, coupling concavity, order-L Lipschitz control, uniform convergence of energy densities, almost additivity, length concavity, and positive doubling defects together do not prove boundedness or upper semicontinuity of the Casimir limit. The finite-L functions in (17) are even rational in t. Semialgebraicity of finite-L data alone does not fix the inference.

Positive degree-one homogeneity can be imposed too. For s>0 and 0<=u<=s set Ehat_L(s,u)=s E_L(u/s). This concave perspective is homogeneous, with (17) as its s=1 slice. Homogeneity alone does not remove the obstruction. No local Hamiltonian realization is claimed.

## 7. The exact missing compactness hypothesis

[Proved] Let C be any subset of K_d on which physical A is finite and nonnegative. The following are equivalent:

1. A is bounded on C.
2. Every h in K_d has a neighborhood U_h and finite B_h such that A(g)<=B_h for g in U_h intersect C.
3. A has a finite nonnegative upper semicontinuous majorant on K_d.

A constant majorant proves 1 implies 3. Upper semicontinuity gives local upper bounds, proving 3 implies 2. A finite subcover of the neighborhoods in 2 gives a global bound max_j B_j, proving 2 implies 1.

The canonical nonnegative extended-real majorant is

    A^*(h)=inf_(r>0) sup {A(g): g in C, dist(g,h)<r},         (21)

with empty supremum defined as zero. If A^*(h)<b, some ball has supremum below b. Every sufficiently close h' has a smaller ball contained in that ball, so A^*(h')<b. This proves upper semicontinuity. If U is any nonnegative upper semicontinuous majorant, its local upper bound at h yields A^*(h)<=U(h). Thus (21) is the smallest such majorant.

An exact extension also requires A^*(h)=A(h) for h in C, which demands relative upper semicontinuity of A on C. A majorant is sufficient for the bound, so exact extension is stronger than necessary.

[Open] The substantive missing statement is finiteness of (21), equivalently local boundedness at every singular bond. Assigning zero or another finite number at nonconformal bonds does not prove it. Defining the value there by (21) makes semicontinuity automatic but may assign infinity. Compactness of an extended-real upper semicontinuous function does not yield a finite bound.

[Proved] If A is unbounded on C, choose h_k with A(h_k)>=k. A subsequence converges to h in K_d, and every neighborhood of h contains arbitrarily large A. Hence A^*(h)=infinity. Conversely, an infinite value of (21) produces a divergent sequence approaching h. This is the precise blow-up alternative. An upward value at a critical point is allowed by upper semicontinuity. Unbounded nearby values at a point assigned a finite value are not.

[Assumed stability] An open neighborhood containing only phases assigned A=0 satisfies the local bound. A gap above an arbitrarily degenerate ground space is not itself such an open-neighborhood hypothesis. This distinction is already recorded in C16.

[Heuristic] Near a finite-c fixed point, one could seek a common continuum description for nearby lattice theories and a bound on their velocity-weighted coefficients. This requires uniform matching and exclusion of additional sectors below the matching scale. It does not follow from finite-volume continuity or the usual Lorentz-invariant c-theorem for one specified flow.

## 8. What a divergent sequence would require

[Assumed infrared identification] A finite-range Lieb-Robinson bound supplies a finite propagation speed controlled by the bounded bond norm. If the proposed conformal sectors are genuine local scaling-limit modes, their velocities satisfy a corresponding bound v_a<=V_d at M=1. The lattice theorem is standard, see Nachtergaele and Sims, [arXiv:math-ph/0506030](https://arxiv.org/abs/math-ph/0506030). Identifying the causal bound with every proposed infrared sector is additional scaling-limit input here.

[Proved conditional implication] Then A<=V_d sum_a c_a, so divergent A requires divergent total central charge. It need not require an increasing number of modes whose velocities stay away from zero. An abstract list of k sectors with c_a=1 and v_a=k^-1/2 has A=sqrt(k), even though its fastest velocity tends to zero. No fixed-d lattice realization of this list is asserted.

Define the weighted mode count

    N_h(u)=sum_(a: v_a>u) c_a.

Integrating the indicator for each sector gives

    A(h)=integral_0^V_d N_h(u) du.                         (22)

An h-uniform integrable envelope N_h(u)<=G_d(u) would give an explicit bound. Finite bounds at each fixed u>0 do not suffice unless their growth near zero is integrable. In the k-sector example, eventually no modes exceed any fixed positive threshold, but the integral diverges. Equation (22) restates the needed velocity tradeoff, not a proved mode-counting estimate.

[Open] No actual fixed-d normalized bond sequence with A tending to infinity was constructed, and no general exclusion theorem was proved. The inherited quasi-free bound A<=2M and the recorded XXZ formulas exclude it in those classes only.

## 9. Making a compactness bound explicit

[Proved conditional implication] A finite cover {U_j} with certified local bounds B_j yields F(d)=max_j B_j. This is explicit only if the cover and its bounds are supplied quantitatively. Existence of a finite real-valued upper semicontinuous extension gives no numerical value or effective maximum-finding procedure by itself.

A finite-scale sufficient condition makes the missing data concrete. If a common length L_0 and an explicit r_d(L_0) obey

    A(h)<=Q_(L_0)(h)+r_d(L_0) for all h in C,                (23)

then (5) gives, at M=1,

    F(d)<=12L_0/pi+r_d(L_0).                                (24)

A continuous decreasing sequence Q_L+r_d(L) converging to A on a common domain would imply (23). Section 5 rules out the simple universal choice r_d(L)=C/L^2 already at d=2. No more general suitable correction was found.

For an epsilon-net {h_j} of K_d in its quotient metric, (4) also gives

    sup_(h in K_d) Q_L(h)
       <=max_j Q_L(h_j)+(12/pi)L^2 epsilon.                (25)

At each net point the energy density can be bracketed using (5):

    |e(h_j)-E_N(h_j)/N|<=2/N.

Replacing e by that finite-volume approximation adds at most 12L^2/(pi N) to Q_L. Thus finite-L suprema admit finite matrix bounds at a certified resolution. What is missing is (23), or another controlled tail bound. A bare finite-size search cannot certify it.

No algorithmic impossibility theorem is asserted. Compactness supplies neither the crucial local boundedness hypothesis nor the quantitative information needed for an explicit F(d).

## 10. Outcome

[Proved] The gauge quotient is compact. Raw limits have sector problems. Uniform two-sided convergence fails near ordinary critical points. The Casimir estimator has no universal monotonicity, and no universal C M/L^2 correction repairs it. The elementary analytic and energy-comparison properties admit the unbounded compact countermodel (17).

[Open] A finite upper semicontinuous physical majorant, local boundedness at every singular bond, and a general finite F(d) remain unproved. These results rule out specific elementary implementations of the compactness route, not the velocity-weighted bound itself.
