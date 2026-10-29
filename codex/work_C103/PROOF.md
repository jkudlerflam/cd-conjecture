# C103. Physics argument and its missing input

2026-10-29.

## 1. Scope and conventions

[Assumed] The microscopic Hamiltonian is a one-site translation-invariant nearest-neighbor chain with finite site dimension d. Its complete infrared theory consists of finitely many unitary CFT sectors with positive velocities in physical lattice units. The coefficient is

    A = sum_a c_a v_a,
    E_0(L) = e_0 L - pi A/(6L) + o(1/L).

Use a compatible vacuum sequence of ring lengths and boundary sectors. For different right and left speeds, v_a means their average. Thermodynamic and dynamical CFT identifications used below are additional parts of the standard continuum description, not consequences of the finite-size formula alone. All explicit examples below have equal left and right speeds.

[Definition] The bond scale is

    M = inf_o (lambda_max(h + o tensor I - I tensor o)
               - lambda_min(h + o tensor I - I tensor o))/2.

A scalar shift permits a representative g with ||g||=M, or an arbitrarily close approximant. Both A and M scale linearly when H is multiplied by a positive constant. Logarithms in entropies are natural, whereas log_2 d = ln d/ln 2 in the proposed inequality.

[Open] No general derivation of A <= pi M log_2 d is obtained here. The argument below identifies a genuine scale-free mechanism for current-complete theories and tests why it does not yet extend to all chains. Two explicit free-chain calculations also disprove proposed monotonicity shortcuts. Neither is a counterexample to the inequality.

## 2. The strongest positive mechanism available here

### 2.1 An optical budget

[Assumed Q1] There are commuting conserved on-site charges Q^alpha, with bond conservation [h,Q^alpha tensor I + I tensor Q^alpha]=0. A uniform phase phi_alpha is introduced across each bond by conjugation on its right site. On a finite ring this defines H(phi). Choose a smooth ground-state branch at phi=0 and take derivatives before the thermodynamic limit. Suitable lengths or a background twist can remove shell degeneracies.

[Proved gauge qualification] Average an arbitrary coboundary o over the on-site charge rotations. The bond is invariant and spectral width is convex, so this cannot increase the width. The optimization of M can therefore be restricted to o commuting with the charges. Such a coboundary still telescopes after introducing the uniform phases. Thus this construction does not require replacing M by a larger, symmetry-restricted bond scale.

[Proved] Write J_alpha=partial_alpha H and T_alpha_beta=partial_alpha partial_beta H at zero phase. Second-order perturbation theory gives

    D_alpha_beta = lim_L (1/L) partial_alpha partial_beta E_0(phi)
                = K_alpha_beta - R_alpha_beta,
    K_alpha_beta = lim_L <T_alpha_beta>/L,
    R_alpha_beta = lim_L (2/L) sum_(n>0)
       Re[<0|J_alpha|n><n|J_beta|0>]/(E_n-E_0).

For every real vector x, x^T R x is a sum of nonnegative squared matrix elements divided by positive gaps. Thus R >= 0 and D <= K in matrix order. This is an all-frequency optical inequality and needs no conformal matching frequency. The stiffness is not the curvature of the infinite-volume energy after minimizing over winding sectors, which would erase the response.

[Assumed Q2] The complete infrared stiffness has the form

    D = (1/pi) sum_a v_a G_a,

where the matrices G_a are positive semidefinite and are fixed by the charge normalization and infrared current algebra or Luttinger parameters. There is no additional neutral CFT sector omitted from this formula. For a single boson with unit-normalized conserved particle charge, G=K_LL and c=1.

[Assumed Q3] A positive matrix W, calibrated from independent current-algebra or Luttinger-parameter information, satisfies

    tr(W G_a) >= c_a for every sector a.

[Proved conditional implication] Multiplying Q2 by W and summing gives

    A <= pi tr(W D) = pi tr[W(K-R)] <= pi tr(W K).       (1)

This sums all sectors before applying the microscopic budget. It is the valid version of a channel-capacity argument when the channels are physical conserved currents.

[Assumed Q4, microscopic budget] If an independent lattice estimate proves

    tr[W(K-R_lower)] <= M log_2 d,
    0 <= R_lower <= R,                                 (2)

then (1) proves

    A <= pi M log_2 d.                                 (3)

[Proved scope] Equations (1)-(3) contain no unknown crossover scale. Q1 and the spectral positivity step are microscopic statements. Q2 and Q3 are standard inputs in specified solvable Luttinger liquids, but not in arbitrary critical chains. Q4 is an additional inequality, not an established consequence of unitarity or locality. Setting R_lower=0 recovers the mechanism used in C101. Supplying R_lower independently is important, since setting R_lower=R and using only infrared information about D would just recast the desired bound.

### 2.2 A nontrivial lower estimate on the regular optical weight

[Proved] At finite L and fixed W >= 0 define the positive measure

    dmu_W(omega) = (1/L) sum_(n>0) sum_alpha_beta
       W_alpha_beta Re[<0|J_alpha|n><n|J_beta|0>]
       delta(omega-E_n+E_0) d omega.

Let s_j = integral omega^j dmu_W. Then tr(W R)=2s_(-1). Cauchy-Schwarz applied to sqrt(omega) and 1/sqrt(omega) in this measure gives

    s_0^2 <= s_1 s_(-1),
    tr(W R) >= 2 s_0^2/s_1.                            (4)

If the measure vanishes, its contribution is zero and the quotient is defined as zero. The statement can be used at finite size and passed to a limit when the moments converge. Here

    s_0 = (1/L) sum_alpha_beta W_alpha_beta
                    <delta J_alpha delta J_beta>_real,
    s_1 = (1/(2L)) sum_alpha_beta W_alpha_beta
                    <[J_alpha,[H,J_beta]]>_real.

The second identity follows by inserting energy eigenstates into the double commutator. Consequently a sufficient scalar version of Q4 is

    tr(W K) - 2 s_0^2/s_1 <= M log_2 d.                (5)

[Proved conditional conclusion] Under Q1-Q3 and the thermodynamic limits just stated,

    A <= pi [tr(W K) - 2 s_0^2/s_1].                   (6)

This retains a certified piece of the regular optical weight and is stronger than discarding it. It does not evaluate the left side of (5) from d and M. In particular, s_0 is a spatially integrated ground-state current covariance, not a single-bond norm.

[Proved limitation] Equality in (4) requires the positive measure to be supported at one nonzero frequency. Thus a single-moment estimate need not be sharp enough for a saturating chain with a broad regular spectrum. Higher moments or additional operator identities would be needed. Equation (6) is a useful restricted-class inequality, not the general theorem requested.

[Numerical] On six- and eight-site XXZ rings, checks.py verifies the double-commutator expression for s_1, the inequality (4), and D=K-R against an independent finite difference of the twisted ground-state energy. These finite-size tests check the spectral algebra and its factors, not the thermodynamic continuum matching.

### 2.3 Heisenberg tests the budget, not just its dimensions

[Assumed standard exact data] For h=S_i dot S_(i+1), M=1/2, e_0=1/4-ln 2, v=pi/2, c=1 and K_LL=1/2. These are the inherited Heisenberg normalization and continuum data. For a twist of S^z, Q3 requires W>=2.

[Proved from these inputs] The phase multiplies S_i^+ S_(i+1)^- by exp(i phi), so the second derivative of the bond is minus its XY part. SU(2) invariance gives

    K = -<Sx Sx + Sy Sy> = -(2/3)e_0
      = (2/3)(ln 2 - 1/4) = 0.2954314537...,
    D = v K_LL/pi = 1/4,
    R = K-D = 0.0454314537....

For the optimal scalar W=2,

    W K = (4/3)(ln 2 - 1/4) = 0.5908629074... > M.

[Proved objection] The simplest diamagnetic budget Q4 with R_lower=0 fails on the saturator itself. It gives A <= 1.85625... rather than A=pi/2. Recovering the sharp bound by this route requires all of the displayed regular weight at Heisenberg. This is a substantive missing step. The optical argument is not a universal proof with an innocuous choice of normalization.

[Proved obstruction to current completeness] A neutral sector has G_a=0, so tr(WG_a)>=c_a>0 is impossible. For example a critical Ising sector can contribute c=1/2 without any continuous on-site charge current. Non-Abelian WZW theories can also contain a nontrivial coset relative to their Cartan currents. Q2-Q3 cannot be assumed for an arbitrary list of unitary sectors. Adding a conserved energy current does not repair this at zero temperature, since its thermal coefficient has extra powers of temperature, as derived below.

## 3. What a velocity-weighted c-theorem would have to establish

[Assumed restricted relativistic setting; proved implication] If a flow has one fixed physical velocity v at every scale and obeys the usual unitary relativistic c-theorem, then

    c_IR <= c_UV  implies  A_IR=v c_IR <= v c_UV=A_UV.

The same conclusion holds for a tensor product of independent such flows with fixed velocities. This is a legitimate scale-free continuum argument. To prove the lattice claim one would still need an ultraviolet embedding with A_UV <= pi M log_2 d in the original microscopic units.

[Proved distinction] A finite-dimensional lattice site is not a UV CFT, and ln d is not its ultraviolet Virasoro charge. Standard blocking generates longer-range interactions and changes the block dimension to d^b. If physical time is kept fixed, the velocity measured in block-site units is v/b. Setting it back to one by rescaling time cannot be used to compare it with the original M. These bookkeeping changes do not supply a monotone normalized by M log_2 d.

[Heuristic, not a theorem] At a perturbative multivelocity fixed point one might try a running expression sum_a v_a(l)c_a(l). Even if each c_a had a signed beta function, differentiation would contain

    d/dl sum_a v_a c_a = sum_a [v_a dc_a/dl + c_a dv_a/dl].

There is no signed formula for the second term supplied by the usual relativistic c-theorem. Away from a product fixed point the individual c_a(l) need not even have a canonical meaning. This identifies missing structure, not a counterexample to every possible weighted theorem.

[Literature scope] A monotone RG potential in a weak-coupling multivelocity model is not automatically a Casimir or current coefficient. [Chang, Chen and Lin](https://arxiv.org/abs/cond-mat/0508660) construct such a potential at one loop. [Hartman and Mathys](https://arxiv.org/abs/2310.15217) construct spectral and null-energy versions of the two-dimensional c-theorem within relativistic QFT. Neither result supplies the required lattice ultraviolet budget. These sources were checked through the browsing tool, with no network shell commands.

[Open] I have neither proved nor disproved endpoint ordering A_IR<=A_UV for the full class of local flows between products of unitary CFTs with physical time fixed. The next two examples rule out identifying a weighted monotone with two natural observables along a general lattice flow. They do not rule out a differently constructed function.

## 4. The Casimir estimator need not decrease toward the infrared

### 4.1 A nearest-neighbor realization

[Proved] Start with a complex-fermion chain in an auxiliary site coordinate r,

    H = -sum_r [t_1 f_r^dag f_(r+1) + t_3 f_r^dag f_(r+3) + h.c.],
    epsilon(k) = -2t_1 cos k - 2t_3 cos 3k,
    0 < t_3 < t_1/3.

Factorization gives

    epsilon(k) = -2 cos k [t_1-3t_3+4t_3 cos^2 k].

The bracket is strictly positive. There are only two Fermi points, at +/-pi/2, and c=1. Expanding at pi/2 gives

    epsilon(pi/2+p) = v_0 p + b p^3 + O(p^5),
    v_0=2(t_1-3t_3),  b=-t_1/3+9t_3.

[Proved] Group three consecutive auxiliary sites into each physical site. The physical dimension is d=8 and the Hamiltonian is now one-site translation-invariant and nearest neighbor. With three-component f_j, its on-site and intersite matrices are

    A_0 = -t_1 (|1><2|+|2><1|+|2><3|+|3><2|),
    B = -t_3 I_3 - t_1 |3><1|,
    H = sum_j [f_j^dag A_0 f_j + f_j^dag B f_(j+1) + h.c.].

Split A_0 equally between adjacent bonds. A Jordan-Wigner representation turns each even fermionic term on one or two adjacent physical sites into an ordinary operator on those sites. Thus this is within the stated finite-dimensional nearest-neighbor class, not a range-three exception to it. Boundary parity selects the usual fermionic spin structure on a ring.

[Proved units check] If the auxiliary chain has N=3L sites, the velocity in physical site units is

    v=v_0/3,  A=2(t_1-3t_3)/3.

Failing to divide by three would make the comparison with the physical bond incorrect.

### 4.2 Exact finite-ring energy

[Proved] Choose N divisible by four and the antiperiodic fermion vacuum sector, so no level lies at a Fermi point. The occupied interval is (-pi/2,pi/2). For odd p its discrete cosine sum is

    sum_(occupied k) cos(pk) = sin(p pi/2)/sin(p pi/N).

It follows that

    E_0(N) = -2t_1 csc(pi/N) + 2t_3 csc(3pi/N),
    e_infty = -2t_1/pi + 2t_3/(3pi).

The antiperiodic half-filled sea contains N/2 fermions, an even number. In the periodic fermionic sector the negative levels contain N/2-1 fermions, an odd number, and the two endpoint modes have zero energy. The minimum periodic-sector energy is

    E_P(N)=-2t_1 cot(pi/N)+2t_3 cot(3pi/N).

It follows directly that

    E_P(N)-E_0(N)
       =2t_1 tan(pi/(2N))-2t_3 tan(3pi/(2N))
       =pi(t_1-3t_3)/N+O(N^-3)>0

for sufficiently large N. These are the two Jordan-Wigner parity sectors of the physical spin ring. Filling any positive-energy level or emptying a negative-energy level cannot lower either sector's energy. Thus E_0(N) is the actual vacuum energy on the stated large-length sequence, not an oscillating non-vacuum estimator.

[Proved] Using csc x=1/x+x/6+7x^3/360+O(x^5),

    E_0(N) = e_infty N
             - pi(t_1-3t_3)/(3N)
             - 7pi^3(t_1-27t_3)/(180N^3) + O(N^-5).

Define the physical-site estimator

    A_ring(L) = (6L/pi)[3 e_infty L-E_0(3L)].

Then

    A_ring(L) = A + [7pi^2/(810L^2)](t_1-27t_3) + O(L^-4).    (7)

For t_1/27<t_3<t_1/3, the coefficient is negative. Hence A_ring(L)<A at sufficiently large L and increases with L along any fixed-ratio compatible sequence. This contradicts nonincreasing Casimir monotonicity even arbitrarily close to the infrared limit.

[Proved free-theorem check] The inherited quasi-free theorem gives M>=||B||_1 and A<=2M. Here ||B||_1>=|B_31|=t_1, and already

    A=2(t_1-3t_3)/3 < 2t_1/3 <= 2M.

Thus the example has ample slack in A<=pi M log_2 8. It invalidates the proposed monotonic observable, not the desired bound. Taking t_3=0 reverses the sign of the correction in (7), as expected for the ordinary XX dispersion. The sign can be changed by an irrelevant dispersion term without changing the infrared central charge.

[Numerical] checks.py verifies the exact energy by an independent occupied-level sum, verifies the blocked Bloch spectrum against the unfolded dispersion, and evaluates (7) on L=8 through 96 at t_1=1,t_3=0.1. The analytic argument is the proof; these are normalization checks.

## 5. The thermal current coefficient is not a monotone either

### 5.1 Why it measures A in the infrared

[Assumed standard CFT thermodynamics; proved calculation] In one sector,

    <H_R>/L = <H_L>/L = pi c T^2/(12v).

Canonical energy fluctuations obey Var(H_R)/L=T^2 partial_T(<H_R>/L)=pi c T^3/(6v), and similarly on the left. The two chiral sectors are independent in the thermodynamic CFT. Since J_E=v(H_R-H_L),

    Var(J_E)/L = pi c v T^3/3.

Adding independent sectors gives

    R_J(T) := 3 Var(J_E)/(pi L T^3) -> A.               (8)

This is an infrared normalization. It is not a positivity theorem for its temperature derivative. A lattice thermal Drude identification requires additional transport assumptions; the next example has a conserved current and avoids that complication.

### 5.2 An explicit XX calculation

[Proved] Consider H=-t sum_j(f_j^dag f_(j+1)+h.c.) at half filling, with t>0. Its dispersion is epsilon(k)=-2t cos k and the energy current from the nearest-neighbor continuity equation is

    J_E=sum_k epsilon(k) epsilon'(k) n_k

up to an irrelevant overall sign convention. It commutes with H. In the grand-canonical state at chemical potential zero,

    Var(J_E)/L = integral_(-pi)^pi dk/(2pi)
         epsilon(k)^2 epsilon'(k)^2 f_T(epsilon)[1-f_T(epsilon)],
    f_T(epsilon)=1/(1+exp(epsilon/T)).

There are two momenta at each energy in the band. Changing variables gives

    Var(J_E)/L = (1/pi) integral_(-2t)^(2t) d epsilon
         epsilon^2 sqrt(4t^2-epsilon^2) f_T(epsilon)[1-f_T(epsilon)].

[Proved] Expand sqrt(4t^2-epsilon^2)=2t-epsilon^2/(4t)+O(epsilon^4/t^3). The required integrals are

    integral_R epsilon^2 f_T(1-f_T) d epsilon = pi^2 T^3/3,
    integral_R epsilon^4 f_T(1-f_T) d epsilon = 7pi^4 T^5/15.

The band-edge errors are exponentially small at low T. Inserting these into (8) yields

    R_J(T) = 2t - 7pi^2 T^2/(20t) + O(T^4/t^3).        (9)

Consequently R_J increases as T is lowered toward zero. A nonincreasing function under coarse graining toward the infrared would have the opposite behavior.

[Proved ultraviolet check] At high temperature f_T(1-f_T)->1/4. The Brillouin-zone average of epsilon^2 epsilon'^2 is 2t^4, so

    Var(J_E)/L -> t^4/2,
    R_J(T) ~ 3t^4/(2pi T^3) -> 0.

Its lattice ultraviolet limit is zero, although A=2t>0. One cannot start this normalized observable at the lattice ultraviolet endpoint and apply downward monotonicity to obtain the conjectured bound.

[Proved consistency] The XX bond has eigenvalues {-t,0,0,t} and M=t, so A=2M. This saturates the free-fermion theorem and has Q=2/pi<1. Equation (9) is consistent with the theorem. It also shows that current conservation and spectral positivity are insufficient for the proposed monotonicity.

[Numerical] Independent integrations over energy and momentum agree in checks.py. The low-temperature expansion and high-temperature limit are checked separately. These computations use no fitted many-body spectrum.

## 6. Why entanglement capacity does not close the argument

[Inherited proved bound] Small incremental entangling gives an instantaneous entropy rate no greater than 18 M ln d across a cut, with arbitrary ancillas. On-site coboundaries can be absorbed into terms acting on only one side, which do not change the cut entropy. C69 and C99 establish the application and its qualifications.

[Assumed continuum preparation] For the equal-spatial-extrapolation-length boundary-state protocol in the supplied background,

    dot S = pi A/(12 ell).

At a controlled ell, comparison with the microscopic rate gives

    A <= (216/pi) ell M ln d,

before matching errors. The length cannot be put equal to one just because a cut intersects one bond. Infrared sector wave packets may require many sites for their continuum description. The physical entropy transported by one sector per time is proportional to c_a v_a divided by its preparation length, not to c_a v_a alone.

[Proved logical distinction] SSA constrains entropies of microscopic tensor factors. A decomposition into infrared sectors is not a decomposition of each d-dimensional site into independent factors. Treating it as one would incorrectly reintroduce a bound on sum c_a, contradicted by the campaign's slow-sector constructions. LR controls propagation speed, but has no central-charge multiplicity factor. Combining a bound on max v_a with a nonexistent bound on sum c_a does not control A.

[Inherited proved limitation] C99 treats signed combinations, repetitions and integrations of SIE-controlled rates. Those arguments retain the smallest length at which the conformal formula is controlled. No new joint entropy inequality is established here. A claim of a direct microscopic transport budget for all c_a v_a would be an additional assumption, not standard SSA.

## 7. Why spectral weight at frequency of order M does not close it

[Inherited proved bound] C100 establishes for the optimized local energy density

    0 <= m_1(q,T) <= 8M^3 q^2,
    m_1,CFT(q,T) = pi A T^2 q^2/6
                      + q^4 sum_a c_a v_a^3/(24pi).

Matching at T_* gives a bound proportional to M^3/T_*^2. Infrared matching alone cannot set T_*=const times M. Making this replacement assumes the missing uniform window.

[Proved scaling obstruction for this comparison] At zero temperature the stress-tensor contribution to the positive first moment is of order q^4, while the available lattice estimate is of order q^2. Dividing by q^4 therefore diverges as q tends to zero. At finite temperature the coefficient carrying A has the additional factor T^2. Neither limit extracts A from the given upper bound with a uniform constant.

[Inherited CFT distinction] The local energy-current spectral density at zero temperature is proportional to c omega^3, with the velocities canceling in the local temporal correlator. The spatially integrated thermal current variance measures A, as in (8). These are different observables. Positivity of one cannot be used with the coefficient of the other.

[Proved logical distinction] Low-energy positive spectral weight can be bounded by a full-frequency moment, but this is useful only if the weight defining A has a finite normalization without a matching scale. The charge-stiffness route in Section 2 has such a Ward identity. The energy-density positive moment above does not. The inverse-frequency static susceptibility has a CFT q^2 term A/(12pi), but also allows ultraviolet analytic terms and density-improvement terms on the lattice. C100 already records this qualification. No universal bound on that inverse moment is assumed here.

## 8. Checks against C101 and the sharp free theorem

[Assumed inherited continuum identifications; proved comparisons in C35 and C101] The table records the tests without re-deriving settled Bethe constants. Q denotes A/(pi M log_2 d).

| Family | Q or an established upper bound | Implication for the proposed argument |
| :-- | :-- | :-- |
| XXZ, Delta=cos gamma, -1<Delta<=1 | sin(gamma)/gamma | Heisenberg saturates. Changing Delta is a path through Hamiltonians, not an RG trajectory. |
| Zero-field critical XYZ equivalents | Same XXZ ratio | Consistent after the allowed unitary transformations and vacuum-sector choice. |
| SU(N) fundamental | 2(N-1)/(N log_2 N) | Equality at N=2. No dimension-independent bound on c follows. |
| SU(N), symmetric rank k | 2k(N^2-1)/[N(N+k)H_k log_2 binom(N+k-1,k)] | At most one; below 3/4 for k>=2. |
| Spin-S Takhtajan-Babujian | N=2 symmetric-family ratio | Included with the inherited normalization. |
| SO(N) vector | 1/log_2 N | Consistent. |
| Sp(2r), r>=2 | [2(r-1)/r+1/(r+1)]/log_2(2r) | Uses the assumed complete multivelocity assignment, not a single speed. |
| Balanced repulsive Hubbard | <=3/(2pi) | Complete charge and spin stiffnesses with K_c>=1/2 give a scale-free proof. |
| Balanced attractive Hubbard | <=1/(2pi) | Requires the stated spin-gap and K_c>=1 inputs. |
| Symmetric supersymmetric t-J, q=2 | <=3/[2 log_2 3] | C101 uses a four-site operator certificate in addition to stiffness. |
| Symmetric supersymmetric t-J, q>=3 | <=2/log_2(q+1) | Complete currents and K_c>=1/q provide the budget. |
| Any quasi-free chain, d=2^m | <=2/(pi m) | The sharp theorem A<=2M is stronger than the conjecture. |

[Inherited qualifications] Arbitrary graded integrable chains and arbitrary multisea dressed-charge matrices are not covered by C101. The Sp assignment has the stated analytical and finite-size support, not an all-rank rigorous completeness theorem. No missing cases are promoted to established results here.

[Proved check using inherited formulas] Symmetric SU(8), k=2 has d=36, M=3/2 and A/pi=3.15. Hence A/(2pi M)=1.05, while Q=0.406195448.... A universal argument eliminating the dimension dependence in favor of A<=2pi M would be false. The free theorem A<=2M applies only to quasi-free Hamiltonians, not to this interacting chain.

[Proved normalization check] For m decoupled half-filled XX species, A=2 sum_i t_i and M=sum_i t_i, so A=2M while Q=2/(pi m). For m decoupled reflection-symmetric Heisenberg species, A=pi M and Q=1/m. Tensor products of saturating qubit chains do not saturate the logarithmic bound at their larger site dimension.

[Proved logical check] The XXZ increase of A/M from 2 at the free point to pi at Heisenberg disproves extrapolating the sharp free constant to interactions. It does not disprove endpoint RG monotonicity, since changing an exactly marginal infrared parameter is not an RG flow from the free lattice Hamiltonian to the interacting one. Likewise C101's Hubbard point with v_c>2t invalidates a sectorwise free-speed ceiling, not the summed bound.

[Numerical] CHECKS.md evaluates representative inherited formulas, including the SU(8), k=2 test, and independently verifies the two new counterexamples to monotonicity. The small checks all pass. No new Bethe ansatz calculation or large diagonalization was performed.

## 9. What remains to be shown

[Heuristic] The viable physical picture is an additive budget for the infrared stress response carried through a microscopic bond. The successful charge-stiffness examples and the free-fermion trace-norm theorem implement versions of this picture. They work because an exact microscopic structure relates the infrared coefficients to a positive finite budget before the infrared limit is taken.

[Open] For general interacting chains, the needed structure must include neutral sectors, respect the gauge optimization of M, retain the bond-energy pairing rather than bound ||grad_h A||, and allow arbitrarily delayed conformal onset. A new Ward identity or an independently controlled susceptibility could do this. The round-one gradient divergence rules out replacing this task by a universal gradient-norm estimate.

[Open] A weighted c-theorem with a lattice ultraviolet normalization would solve the problem, but both that normalization and the relevant monotonicity remain additional claims. A complete-current optical theorem with a sharp bound on its remaining spectral weight would also solve a substantial class. Neither follows from the stated infrared assumptions alone in the arguments developed here. The final conjectured inequality remains scale-free and unproved, not replaced by an estimate with an unnamed crossover scale.
