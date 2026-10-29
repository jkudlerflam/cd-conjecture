# C69. Velocity-weighted entropy rates and crossover bounds

2026-10-28. All logarithms are natural unless a base is displayed. Lattice spacing and hbar are one. This is a closed-book investigation using the four requested earlier reports. No network access or new many-body simulation was used.

## 1. Result and unresolved question

[Proved conditional result] There is a CFT entropy rate with the requested velocity weighting. Prepare each sector in a regularized conformal boundary state with the same spatial extrapolation length ell. The late-time half-line entropy rate is

    r_CFT(ell) = pi A/(12 ell),     A = sum_a c_a v_a.

The lattice small-incremental-entangling bound is r <= 18 M log d per cut, where M is the optimized bond half-spread. If the lattice realizes this rate with relative lower accuracy alpha at ell = ell_0, then

    A <= [216/(pi alpha)] ell_0 M log d.                 (1)

Thus a uniform ell_0 would prove a universal bound, with an explicit but nonsharp constant. The infrared CFT assumption does not supply a uniform ell_0.

[Proved restricted obstruction] A scale-covariant entropy difference for a geometry of size ell has a physical time derivative proportional to 1/ell. Multiplying the derivative by ell recovers a coefficient proportional to A, but multiplies its SIE bound by ell as well. This excludes a scale-free inference from this class of observables and the SIE inequality alone. It is not a theorem excluding every possible CFT observable.

[Proved crossover result under the stated thermal assumption] Let B = sum_a c_a/v_a and assume at one temperature T_* > 0 that

    s(T_*) >= alpha pi T_* B/3,       0 < alpha <= 1.

Then, for c = sum_a c_a and v_H = c/B,

    c <= [3/(pi alpha)] xi_H log d,
    xi_H = v_H/T_* <= xi_max = v_max/T_*.

In particular,

    c <= [3/(pi alpha)] xi_max log d.                  (2)

No LR theorem or common sector velocity is required. Section 7 gives the factor-of-two conversion to the entropy convention in the question, a conservative LR corollary, and a separate spatial crossover theorem.

[Assumed/open] I have neither proved A <= C M log d with C independent of crossover scales nor disproved it. I also have not proved that no such observable exists. Those two possibilities are not exhausted by the failure of an entropy-rate argument. The sharp proposal A <= pi M log_2 d remains a conjecture. The statements below identify a concrete missing lower bound on entangling capacity, rather than treating the obstruction as a solution of that conjecture.

## 2. Inputs, conventions, and local bounds

[Assumed setting] Consider a nearest-neighbor chain with on-site dimension d >= 2 and a Hermitian translation-invariant bond h. Its infrared description consists of finitely many nonchiral unitary CFT sectors with c_a > 0 and physical propagation speeds v_a > 0. Thus

    c = sum_a c_a,       A = sum_a c_a v_a,
    E_0(L) = e_0 L - pi A/(6L) + o(1/L)

on vacuum-compatible rings. A sector's left and right Virasoro charges both equal c_a. A is the oscillator Casimir coefficient, not a coefficient contaminated by zero-mode rounding on incompatible rings.

[Definition] Write

    M = inf_(o=o^dagger, b real) ||h + o tensor I - I tensor o - b I||.

The minimization over b is the spectral half-spread. The o term telescopes on rings. We work with arbitrarily near-minimizing representatives, so no attainment assumption is needed.

[Assumed input, published SIE theorem] For a bipartition aA | Bb, where A and B are the two finite systems on which a Hamiltonian g acts, the ancillary SIE theorem gives

    |d S(aA)/dt| <= 18 ||g|| log min(dim A, dim B).       (3)

The ancillary systems a and b can be arbitrarily large. This is the full ancillary bound of Van Acoleyen, Marien, and Verstraete, building on Bravyi's entangling-rate work. The numerical constant 18 is used throughout, not an unspecified order-one constant. Entropy and the logarithm in (3) use the same base. We invoke this theorem and do not reprove it here. References are listed in Section 10; they were not downloaded or rechecked in this task.

[Proved from (3)] In a finite open chain, only the bond crossing a fixed cut changes the entropy of one side. Hamiltonians supported wholly on either side generate local unitaries and leave that entropy invariant. Taking A and B to be the adjacent physical sites in (3) yields

    |S_L(t_2)-S_L(t_1)| <= 18 M log d |t_2-t_1|.         (4)

A mixed state is covered by a purification with the purifiers assigned to the appropriate sides and left stationary. For an interval with two boundary bonds the right side doubles. For a cut on an infinite chain, apply (4) in finite volumes and take a limit of finite entropy differences, when that limit exists. We do not subtract two undefined infinite entropies.

[Proved gauge compatibility] Replacing the crossing bond h by h + o tensor I - I tensor o changes it only by operators local to the two sides. Those terms contribute zero to the entropy derivative in every state, since local unitaries preserve the reduced-state spectrum. A scalar also contributes zero. Thus the instantaneous entangling rate of the physical Hamiltonian can be bounded using every gauge representative of h. Taking the infimum gives M in (4), even if the optimizing o has large norm. At the distant ends of an open chain the telescoping terms are also local to one side. This argument is specific to physical entropy rates; it does not assert that all local definitions of energy current are gauge invariant.

[Proved] Let

    E(h) = sup |dS_L/dt|

where the supremum allows finite ancillary systems and all initial pure states across a single bond. It is invariant under the same additions of local terms, and E(h) <= 18 M log d. A sufficient new theorem for a universal Casimir bound would be

    E(h) >= kappa A                                     (5)

for a dimension-independent kappa > 0. It would imply A <= (18/kappa) M log d. None of the long-distance calculations below proves (5).

## 3. A CFT observable that does measure A

### 3.1 Equal spatial regularization of conformal boundary states

[Assumed preparation] In each sector choose a conformal boundary condition B_a admitting the usual strip one-point function of the replica twist field and its continuation to n = 1. Prepare the normalized state

    |psi_a> proportional to exp(-tau_a H_a)|B_a>,
    ell_a = v_a tau_a.

Take the product over sectors and evolve by the physical Hamiltonian sum H_a. These are pure states. Availability of such a product preparation in the microscopic chain, with controlled errors, is an additional assumption. The existence of a ground-state CFT alone does not provide it.

[Proved CFT calculation under this input] Use the Euclidean coordinate w = x + i v_a tau on the strip 0 < Im w < 2 ell_a. The map z = exp(pi w/(2 ell_a)) sends the strip to the upper half-plane. The full scaling dimension of the replica twist field is

    x_n = (c_a/12)(n-1/n).

Its one-point function transforms with the ratio of |dz/dw| to 2 Im z. At x = 0 this ratio is

    pi/[4 ell_a sin(pi v_a tau/(2 ell_a))].

The ket and bra preparation place the insertion at the strip midpoint. Continuing tau = tau_a + i t turns the sine into cosh(pi v_a t/(2 ell_a)). With a short-distance regulator a_uv,

    Tr rho_half^n = C_n [ (4 ell_a/(pi a_uv))
                         cosh(pi v_a t/(2 ell_a)) ]^(-x_n).

Taking -partial_n at n = 1 uses partial_n x_n|_1 = c_a/6 and gives

    S_a(t)-S_a(0) = (c_a/6) log cosh(pi v_a t/(2 ell_a)),
    dS_a/dt = pi c_a v_a/(12 ell_a)
              tanh(pi v_a t/(2 ell_a)).                (6)

The constants C_n, boundary entropy, and ultraviolet logarithm cancel in the difference. These formulas are statements in the specified continuum preparation, not exact identities for a lattice quench.

[Proved] Set every ell_a = ell. Additivity of entropy for the product state gives

    lim_(t->infty) [S(t)-S(0)]/t = pi A/(12 ell).        (7)

All v_a are positive and their number is finite, so the limit follows term by term from log cosh x = x - log 2 + o(1). The finite-time derivative increases from zero to that value. This is the requested CFT velocity weighting. The size of the preparation is indispensable.

[Proved normalization check] Equal temporal preparation, tau_a = tau_0, instead gives the rate pi c/(12 tau_0). One must set equal spatial lengths, not equal times, to obtain A. The energy-matched effective temperature of sector a in this preparation is 1/(4 tau_a), and its thermal entropy density is pi c_a/(12 ell_a). The interval growth slope before saturation is twice the half-line slope. These checks reproduce the two boundary cuts in the lattice bound.

### 3.2 What is needed to transfer the result to a lattice

[Proved finite-time estimate] An infinite-time plateau is not necessary for a conditional theorem. Let v_min = min_a v_a, choose q > 2 log 2/pi, and set T_obs = q ell_0/v_min. Since log cosh x >= x-log 2 for x >= 0, (6) implies

    Delta S_CFT/T_obs >= pi A/(12 ell_0) - c log 2/(6 T_obs)
      >= [1-2 log 2/(pi q)] pi A/(12 ell_0),

where the last step uses c v_min <= A. If the lattice entropy increment over this one interval is below the CFT increment by at most delta [pi A/(12 ell_0)] T_obs, set alpha = 1-2 log 2/(pi q)-delta > 0. Applying (4) already proves (1). For example q = 4 gives alpha = 1-log 2/(2 pi)-delta. The finite-volume boundaries must be far enough away for the assumed matching to hold. This requires a longer observation time when a sector is slow, not a change in the coefficient A.

[Assumed matching, finite-increment form] Suppose a family of finite-chain states at a fixed preparation length ell_0 has time intervals of duration Delta t, before boundary returns, on which

    Delta S_L >= alpha [pi A/(12 ell_0)] Delta t - o(Delta t),
    0 < alpha <= 1,

with a limit in which the error divided by Delta t vanishes. Equivalently, one may assume a matching entropy-growth slope in a controlled thermodynamic limit. A lower estimate suffices. Ground-state scaling, local correlator convergence, or convergence of entropy values at isolated times is weaker than this hypothesis. In particular, the derivative of a small entropy error need not be small.

[Proved] Apply (4), divide by Delta t, and take the stated limit. This gives (1). There is no exchange of a derivative with a lattice scaling limit. A finite interval gives the same coefficient, because both its CFT growth slope and its SIE upper bound double. If only ell >= ell_cross is controlled, the conclusion is A <= [216/(pi alpha)] ell_cross M log d, when matching holds at ell_cross. If matching is only asymptotic as ell grows, no bound at the endpoint ell_cross is justified.

[Assumed/open] At fixed finite lattice preparation length, irrelevant interactions can eventually change a quench state or equilibrate sectors prepared at different effective temperatures. An arbitrarily long conformal plateau is not automatic. The finite-increment hypothesis can instead be imposed on a controlled interval with an explicitly bounded error, in which case that error must be retained. This is a second reason why a statement about the infrared vacuum cannot silently be upgraded to the matching hypothesis above.

[Proved conditional implication] If all chains in the proposed class admitted the matching hypothesis with common alpha > 0 and ell_0 <= L_0, then

    A <= C M log d,       C = 216 L_0/(pi alpha).

This proves the claimed implication, not the uniform-preparation premise. It also supplies no argument for the sharp coefficient pi/log 2.

## 4. What the usual current and vacuum observables measure

[Proved CFT identities under the nonchiral convention] The right-moving thermal energy and entropy densities of a sector are

    e_R(T) = pi c_a T^2/(12 v_a),
    s_R(T) = pi c_a T/(6 v_a).

Multiplication by the speed cancels the inverse velocity. Perfectly transmitting chiral reservoirs at common temperatures therefore give

    J_E = (pi/12) c (T_L^2-T_R^2),
    J_S = (pi/6) c (T_L-T_R).                           (8)

The second expression is the net thermodynamic entropy flux, not automatically the derivative of a chosen subsystem's von Neumann entropy. Any use of SIE for it needs a purification and a proved identification of that subsystem entropy rate with the flux. Thus even before the velocity issue, an entropy-flux formula alone is insufficient.

[Proved conditional observation] Sector-dependent hot temperatures T_(a,L) = v_a/ell and cold right movers give J_S = pi A/(6 ell). This again measures A/ell. Preparing such reservoirs and identifying their flux with a lattice entropy rate are extra assumptions. They do not supply a scale-free alternative to (7).

[Proved, using the already derived C35 identity] The CFT total energy-current variance per length obeys

    Var(J_E)/L = (pi/3) A T^3.

C35 RESULT Section 7 derives this identity and explains the ultraviolet spectral-weight problem. For a centered bond representative g, the local estimate ||i[g_12,g_23]|| <= 2 ||g||^2 does not bound the coefficient after division by T^3. Near-minimizing representatives have ||g|| arbitrarily close to M, but changing the density gauge also changes the local current. SIE controls a first derivative of entropy under a crossing Hamiltonian, not this current susceptibility. No new susceptibility estimate is proved here.

[Proved] The ground state of H is stationary under H. Its fixed-region entropies and mutual informations therefore have zero physical time derivative, even though its Casimir coefficient is nonzero. Moving an interval endpoint is a change of subsystem, not time evolution by the crossing bond. For example, S(L) = (c/3) log(L/a_uv) + constant gives dS/dt = c dot L/(3L) under a prescribed endpoint motion. The factor 1/L remains. Boosted or modular time derivatives use a different generator, whose coupling strengths must be bounded separately before applying SIE with physical M.

## 5. A restricted no-go statement and its limits

### 5.1 Scale-covariant entropy differences

[Assumed class of observables] Let F_ell(t) be a fixed finite linear combination of region entropies in a CFT preparation whose distances and smoothing lengths all scale with ell. Keep all dimensionless shapes, sector velocity ratios, and dimensionless state parameters fixed. Suppose scale covariance gives

    F_(lambda ell)(lambda t) = F_ell(t) + b log lambda,

where b is independent of t. The logarithmic anomaly is absent for many mutual-information combinations and cancels from entropy differences in time. This class includes (6) and UV-finite fixed-shape mutual-information or conditional-mutual-information protocols with the stated covariance.

[Proved] Differentiating the covariance relation gives

    partial_t F_(lambda ell)(lambda t)
        = (1/lambda) partial_t F_ell(t).                (9)

Thus any nonzero rate in this class scales as a velocity divided by a length. A rate k A with nonzero k independent of ell cannot hold throughout this scale family. The conformal limit at increasing ell makes the unrescaled rate smaller, not larger.

[Proved consequence of SIE] Suppose the lattice expression is sum_j q_j S(R_j), and region R_j has n_j crossing bonds. Define W = sum_j |q_j| n_j. Equation (4) and the triangle inequality give

    |partial_t F_ell| <= 18 W M log d.

For fixed W, a coefficient extracted as ell partial_t F_ell is bounded only by 18 W ell M log d. Choosing weights proportional to ell makes W grow by the same factor. Increasing the number of cuts also increases the upper bound. Cancellations might permit a new, stronger estimate, but the triangle-inequality application of SIE does not provide one.

[Proved restricted conclusion] Within the assumed scale family and these SIE estimates, neither the continuum limit, a multiplication by ell, nor replication of cuts produces a universal bound on A. This is an obstruction to an inference. It is not a classification of all possible information-theoretic functionals, nor a prohibition on an independent ultraviolet-to-infrared comparison theorem.

### 5.2 The inequalities are compatible with arbitrarily large formal A

[Proved logical nonimplication] Fix d, M > 0, and a speed v_0 below a chosen velocity cap. Let K = 18 M log d. For any formal value c > 0, set A = c v_0 and choose

    ell_0 >= max{1, pi A/(12K), pi c/(12 log d)}.

For all ell >= ell_0 use the CFT entropy-increment profile

    Delta S_ell(t) = (c/6) log cosh(pi v_0 t/(2 ell)).

Its derivative is at most pi A/(12 ell) <= K. Its associated thermal entropy density pi c/(12 ell) is at most log d. Its propagation speed remains v_0. Hence arbitrarily large A is compatible with these rate, entropy-capacity, and speed inequalities if the preparation length grows. These functions are a logical test of the inference, not a constructed finite-dimensional Hamiltonian or a counterexample to the conjecture. C35 Section 7 gives the analogous test for thermal entropy and energy capacities.

[Proved scope limitation] This construction does not show that no other observable exists. Such a conclusion would require a separate theorem or a microscopic counterexample. In particular, a new estimate on stress-tensor susceptibility, a uniform lower bound like (5), or a monotone comparison between microscopic entangling capacity and the infrared coefficient could escape the obstruction.

## 6. A conservative LR bound, with its continuum assumption exposed

[Proved lattice estimate] For a centered nearest-neighbor representative g with ||g|| <= M + delta, let X and Y be single-site observables at distance r >= 1. The standard iterated commutator LR inequality gives a sum over chains of overlapping interaction bonds. The first bond has at most two choices and each subsequent bond at most three, including repetition. At least r bonds are needed to connect the sites. Consequently

    ||[X(t),Y]||
      <= 2 ||X|| ||Y|| sum_(n>=r) (2(M+delta)|t|)^n
                            [2*3^(n-1)]/n!
      <= (4/3) ||X|| ||Y||
         exp[-mu r + 6(M+delta) exp(mu)|t|],            (10)

for every mu > 0. The last step uses 1_(n>=r) <= exp(mu(n-r)) and sums the exponential series. Terms supported on only one side can also be absorbed into an interaction picture; retaining them in this count only weakens the estimate. Coboundaries do not change bulk dynamics. Taking delta down to zero and minimizing exp(mu)/mu at mu = 1 gives the conservative cone speed

    v_LR = 6 e M.                                      (11)

This count deliberately does not use the sharper 2 e M quoted in the earlier reports. Nothing below requires that sharper constant.

[Assumed real-time identification; proved implication] To infer v_a <= v_LR, assume each sector is visible in suitable local lattice observables whose real-time commutators near its propagation cone have a nonzero algebraic scaling limit. A cone with v_a > v_LR would then have a nonzero algebraic signal where (10) is exponentially small, a contradiction. Low-energy eigenvalue fits alone do not state this dynamical scaling hypothesis. Under it,

    A <= c v_LR <= 6 e c M.

This still contains c. LR bounds signal propagation, not the number of conformal sectors. No factor log d enters (10).

## 7. Crossover theorems

### 7.1 Thermal theorem with the weakest necessary matching assumption

[Assumed] The thermodynamic entropy density exists at T_* > 0 and obeys the single lower estimate

    s(T_*) >= alpha (pi T_*/3) B,
    B = sum_a c_a/v_a,       0 < alpha <= 1.            (12)

For a two-sided relative-error window one can take alpha = 1-epsilon. We require neither the upper error estimate nor matching at every smaller temperature. If a window endpoint is defined only as a supremum, use temperatures approaching it from below and pass to the limit, rather than assume the endpoint is attained.

[Proved finite-dimensional input] For L sites, positivity of relative entropy against I/d^L gives

    0 <= D(rho_L || I/d^L) = L log d - S(rho_L).

Divide by L and take the thermodynamic limit to obtain s(T_*) <= log d. Combining with (12) gives

    B <= 3 log d/(pi alpha T_*).                       (13)

Define v_H = c/B. Multiplying (13) by v_H proves

    c <= 3 [v_H/T_*] log d/(pi alpha).                 (14)

Since 1/v_H = sum_a (c_a/c)/v_a, v_H lies between the minimum and maximum sector speeds. Replacing v_H by v_max proves (2). This is an algebraic consequence of the stated thermal matching and finite local dimension, with no entanglement-growth hypothesis.

[Definition and caveat] The operational thermal crossover length here is xi_cross = v_max/T_*. It is measured in sites. It is a thermal wavelength at the temperature where conformal thermodynamics is quantitatively valid. It is not automatically the longest correlation length of an arbitrary operator. Such correlation lengths also depend on scaling dimensions. The stronger weighted length xi_H = v_H/T_* is useful when the sector data are known.

[Proved numerical constants] At relative accuracy epsilon = 0.1, (2) has C' = 3/(0.9 pi), approximately 1.06103. With exact matching alpha = 1, C' = 3/pi. These are constants of the stated wavelength convention, not claims of sharpness over microscopic chains.

[Proved conditional LR corollary] Under the real-time identification in Section 6,

    xi_max <= 6 e M/T_*,
    c <= [18 e/(pi alpha)] (M/T_*) log d.              (15)

This is explicit and conservative. It depends on a crossover temperature, as it must in this argument.

### 7.2 The pi/6 convention in the question

[Proved normalization conversion] If the entropy input is instead written literally as

    s(T_*) >= alpha (pi T_*/6) sum_a c_a/u_a,

then the same proof, with u_H = c/(sum_a c_a/u_a), gives

    c <= [6/(pi alpha)] (u_H/T_*) log d
      <= [6/(pi alpha)] (u_max/T_*) log d.              (16)

This is the requested C' = 6/[pi(1-epsilon)] if xi_cross is defined as u_max/T_* and alpha = 1-epsilon.

[Proved convention audit] There are two usual ways to reconcile this with (14). If c_a retain their ordinary nonchiral meaning and the Casimir convention in Section 2, then u_a = v_a/2 makes the two entropy formulas identical; the length in (16) is half the length in (14). Alternatively pi/6 is the thermal entropy coefficient when a labels individual chiral sectors with their physical speeds. Then sum_a c_a counts both chiralities and equals twice the ordinary nonchiral central charge in a parity-symmetric theory. One cannot keep physical u_a = v_a, nonchiral c_a, and the pi/6 coefficient simultaneously. No physical inequality depends on which consistent convention is chosen.

### 7.3 A spatial crossover theorem without thermal assumptions

[Assumed spatial matching] Let S(n) be the von Neumann entropy of n consecutive sites in an infinite-chain state. For some integers n,m >= 1, suppose its entropy increment satisfies

    S(n+m)-S(n) >= alpha (c/3) log(1+m/n).              (17)

The c/3 coefficient is the ground-state interval coefficient of the nonchiral CFT convention. This assumption controls an increment, so the additive nonuniversal entropy constant cancels. It is stronger than merely asserting the asymptotic coefficient without specifying an onset scale and error.

[Proved] Subadditivity, with the extra m sites called B, gives

    S(n+m)-S(n) = S(B|A) <= S(B) <= m log d.

Therefore

    c <= [3m/(alpha log(1+m/n))] log d.                (18)

At a doubling scale m = n = xi_cross this becomes

    c <= [3/(alpha log 2)] xi_cross log d.             (19)

If the one-site increment is already controlled at n = xi_cross, then

    c <= 3 log d/[alpha log(1+1/xi_cross)]
      <= [3/alpha](xi_cross+1) log d,                 (20)

where log(1+x) >= x/(1+x). This version contains no velocity and no identification of thermal and spatial crossover lengths. The constant in (19) is about 4.32809/alpha. It is weaker than the thermal coefficient when the two lengths happen to coincide, but the assumptions are different and that coincidence is not automatic.

[Proved scope] The bare asymptotic statement S(n) = (c/3) log n + constant + o(1) ensures (17) for large enough doubling scales at any fixed alpha < 1. It places no uniform bound on the first such n. Likewise s(T)/T -> pi B/3 as T -> 0 supplies some sufficiently small T_*, but no uniform lower bound on it. This is why the crossover theorems are meaningful and do not prove c <= constant log d without an onset assumption.

## 8. Checks and consistency

[Proved normalization checks] The conformal half-line rate vanishes at t = 0, approaches pi A/(12 ell) at late times, and doubles for two independent interval boundaries. Rescaling H by a positive factor lambda sends M and all v_a to lambda times their original values. At fixed spatial ell the rate and both sides of (1) scale by lambda. If T_* also scales by lambda, xi_max and both sides of (2) are unchanged. The formulas are additive in independent CFT sectors before imposing the common microscopic norm bound.

[Proved limiting checks] With one sector, v_H = v_max = v and (14) is just pi c T_*/(3v) <= log d divided by alpha. With widely separated speeds, v_H <= v_max keeps the inequality in the correct direction; replacing v_max by v_min would not be justified. In the limit T_* -> 0 the crossover bound becomes weak. Taking ell -> infinity also weakens (1). These outcomes agree with the locations of the missing uniform scale assumptions.

[Assumed inherited model input] The Heisenberg chain in h = P normalization has d = 2, M = 1, c = 1, v = pi. It saturates the proposed sharp inequality A <= pi M log_2 d, so any universal natural-log coefficient C must be at least pi/log 2, approximately 4.53236. Our conditional SIE coefficient 216/pi, approximately 68.7549, is compatible with that requirement and is not competitive with it. For the SU(N) symmetric-rank-k family, use C35's established conditional inequality and the supplied survey's numerical results. They were not recomputed here.

[Numerical, elementary formula checks only] checks.py evaluates the strip derivative, its late-time limit, the geometric LR tail majorant, the velocity-weighted harmonic-mean inequality, and the constants in Sections 3 and 7. It writes CHECKS.md. These checks test factors and algebra, not the CFT identification of a microscopic Hamiltonian, the SIE theorem, or the existence of a uniform preparation scale. No many-body numerical claim is added by these checks.

## 9. What remains to be proved

[Assumed/open] The concrete unresolved target is a dimension-independent lower bound on microscopic entangling capacity in terms of A, such as (5), or a different positive lattice response with a controlled infrared coefficient and no inverse preparation scale. A proof must bridge the lattice to the infrared without evaluating a continuum formula at one site by assumption. The previously identified current susceptibility is another possible target, with its ultraviolet and gauge problems still present.

[Proved conclusion about this investigation] The candidate entropy observable exists, but it measures A/ell. SIE proves (1), and finite-dimensional entropy capacity proves the crossover inequalities (14), (16), and (18). The restricted scaling obstruction does not settle the universal Casimir conjecture and does not justify a claim that every possible observable is excluded.

## 10. Sources and local handoff

[Assumed source inputs] The supplied reports used here are ../../contrib/refined_conjecture/README.md, ../../contrib/casimir_survey/README.md, ../../contrib/lit_casimir_d2/README.md, and ../work_C35/RESULT.md with its PROGRESS.md. In particular, the failed dimension-independent Casimir bound and the all-symmetric-family logarithmic comparison are inherited conclusions. The source reports distinguish their analytic continuum inputs from finite-size numerical evidence.

[Assumed external theorem and CFT framework] The relevant references are S. Bravyi, "Upper bounds on entangling rates of bipartite Hamiltonians," Physical Review A 76, 052319 (2007); K. Van Acoleyen, M. Marien, and F. Verstraete, "Entanglement Rates and Area Laws," Physical Review Letters 111, 170501 (2013); and the conformal boundary-state quench framework of Calabrese and Cardy. The full ancillary constant used here is from the 2013 SIE theorem. The strip calculation needed for this report is given explicitly in Section 3, so it does not rely on an unverified recalled quench coefficient. No fresh literature search was performed.

[Proved work record] All created files are in work_C69. The project-level STATE.md is outside the permitted writable folder and was not edited. PROGRESS.md and the local STATE.md provide the resumable handoff. Nothing was committed, pushed, or submitted to the cluster.
