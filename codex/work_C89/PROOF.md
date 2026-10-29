# C89. Pin defects and the mobile sector comparison

[Status] The requested full mobile ground-state selection is not established. This note proves finite logical-model pin exclusion, derives the finite-volume dressed cost through second order in J/t with quantum clock sources, and gives a nonperturbative large one-marker exclusion in the legal projected model. It formulates a sufficient length-weighted comparison with the correct physical scales. The missing low-frequency bound is identified explicitly. None of these results discharges C88's S_pin for the finite-K microscopic chain.

[Tags] Proved means an argument given here or an explicitly inherited result. Numerical means floating-point calculations at stated finite sizes. Assumed identifies the input still required. Heuristic identifies an infrared extrapolation. Definitions and scope statements are labeled separately. All energies below are canonical unless a chemical potential is explicitly included.

## 1. Inputs and the actual schedule

[Inherited proved] C88 PROOF Sections 1 to 4 give the pinless eleven-letter material bond, its two-color vacancy lift, exact pin conservation, and the legal projection. With N occupied sites and V vacancies, L=N+V, the rooted total-momentum-zero representation is

    H = H_gap + J sum_j A_j h_j,
    H_gap = -t sum_j (E_j^+ E_(j+1)^- + h.c.),
    A_j = 1_(g_j=0), sum_j g_j=V.                         (1)

The shifts have unit matrix elements. The material bond h_j includes both distributed onsite halves. Moving those halves to physical onsite energies changes the model. The physical root and cyclic gluing must be retained when interpreting other momentum sectors.

[Inherited proved] C88, C85, and C83 distinguish the frozen selector certificate from mobile matching. For C88 the selected material length is b=15, with fourteen data sites and thirteen active bits. Write

    E_l = mu_sel - 2 cos(pi/(l+1)), e_b=E_15/15,
    mu_sel = 1.8891880765550189...,
    e_b = -0.004825498950096...,
    min_(l != 15) (E_l/l-e_b) > 2.8116274e-5.            (2)

These are the previously certified constants, not new certificates. The finite-K all-string correction is at most (l-2)_+/(K-4). A pinned length-15 cell has energy at least mu_sel-epsilon_pin, where epsilon_pin=2^32/(K-2^17). The pin cost over a vacuum cell is therefore at least

    r_pin = 2 cos(pi/16)-epsilon_pin > 1.                 (3)

[Definition, inherited schedule] The material parameters are g=2^-272, lambda=2^-608, K>=2^2560, D=1000, Delta=4-2cos(pi/16), t_d=beta 2^-66, and f=lambda^2 E_unit 2^-120. Here beta=1/d_15, d_0=1, d_1=1004, d_(j+1)=1004d_j-d_(j-1), and

    alpha_end=sqrt(2/16)sin(pi/16),
    E_unit=alpha_end^4 beta^2 g^2.                        (4)

[Scope] C88 MEMBER.json specifies no numerical J/t. Its K_mobile is an implicit inequality containing an unspecified local remainder constant and infrared tolerance. Thus there is no scheduled, numerically certified mobile point at which an unknown dressing coefficient can be evaluated. The material minimum K=2^2560 does not resolve this omission. No new power of two for J is certified here.

## 2. Reference costs with the common bulk subtracted

[Inherited proved] The free vacancy ground energy and mean contact in (1) are

    E_gap=-2t sin(pi N/L)/sin(pi/L),
    a_L=N/L-sin^2(pi N/L)/(N L sin^2(pi/L)),
    a(rho)=rho-sin^2(pi rho)/(pi^2 rho).                  (5)

At half occupancy, a=1/2-2/pi^2=0.29735763271532445. The factor a multiplies all material terms in the first-order degenerate perturbation operator. It does not multiply the bare kinetic energy t.

[Proved] At lambda=f=0, K=infinity, the common frozen bulk is e_b per material site. A wrong cell has excess

    d_l=E_l-l e_b,                                       (6)

and a pinned b-cell has excess d_pin=2cos(pi/16). The leading mobile costs are J a_L d_l and J a_L d_pin. Comparing E_l with E_b instead would compare different occupied lengths and would give the wrong selector test. A single wrong cell can be assigned (6) as an island excess; a finite compatible ring must compensate its length elsewhere.

[Numerical arithmetic from inherited constants] At half occupancy,

| Defect | Frozen excess | Leading mobile cost divided by J |
| :-- | --: | --: |
| Pin in a length-15 cell | 1.96157056080646 | 0.583287978365481 |
| Length 14 | 0.000449860388753 | 0.000133769420252 |
| Length 16 | 0.000449860388753 | 0.000133769420252 |
| Length 30 | 0.0442143982741 | 0.0131474888027 |
| Length 60 | 0.181369833791 | 0.0539317044221 |

The equality d_14=d_16 follows from the midpoint selector. A (14,16) replacement of (15,15) therefore costs 0.000267538840504 J to first order at half occupancy.

[Proved] For long cells,

    d_l=|e_b| l+(mu_sel-2)+pi^2/(l+1)^2+O(l^-4).         (7)

The bare excess grows with length. This is the correct resource for absorbing a length-weighted error, including any extensive difference in charge relaxation.

[Definition] Count a wrong-length pinned cell once as a wrong-length cell. Let k_p count pinned cells of length 15 and define

    W=sum_(wrong-length cells) l + k_p.                   (8)

[Inherited proved, applied] C88's frozen cutting proof gives, on compatible rings,

    E_sigma^mat-E_good^mat
      >= gamma sum_wrong l+r_pin k_p-5 B_port k,
    gamma>2^-17, B_port=sqrt(2)|lambda|, k<=W.            (9)

Here E_good^mat is the exact common selected bulk on that ring, including its port energy. Hence the weaker but uniform bound is (gamma-5B_port)W. The all-string finite-K correction is already included in gamma. No two-sided error proportional just to the number of unbounded wrong cells is asserted.

## 3. An exact pin inequality in the logical clock model

[Definition] Restrict (1) to K=infinity, lambda=f=0 and the zero-auxiliary vacuum-clock block. The internal bits are spectators in this block. Cell c of length b_c has clock positions q_c=0,...,b_c-1. Its marker position in material coordinates is m_c. Its diagonal source is

    (mu_sel/2)(A_(m_c-1)+A_(m_c)).

The link q_c -> q_c+1 has amplitude -J A_(s_c+q_c), where s_c is the first data position. A pin freezes the cell at its terminal all-L position and removes these clock links. The marker remains E. This is an exact invariant block at zero port, not a simulation of the full finite-K Hamiltonian.

[Proved, finite-system pin exclusion] Let psi_p be a ground state with any specified nonempty set P of pins, with all other cell lengths unchanged. In the unpinned sector replace each frozen internal vector by the normalized sum of its terminal and penultimate clock vectors. Use the same bath and other-cell state psi_p. These two clock positions have identical diagonal terms. Their one mutual link contributes -J A_(m_c-1). All links to other clock positions have zero expectation. Terms of other cells are unchanged, including the marker energy split across its two contacts. The resulting trial energy gives

    E_pinned-E_unpinned
       >= J sum_(c in P) <psi_p,A_(m_c-1) psi_p>.         (10)

At finite N,V, with J,t>0, the hopping graph of this clock block is connected. Vacancy transfers connect all gap compositions. For every clock link some composition has the required zero gap, so every clock position is reachable. All nonzero off-diagonal entries are negative. Perron-Frobenius therefore gives a strictly positive ground vector on this graph. Each A_(m_c-1) has nonempty support. Every term on the right of (10) is strictly positive.

[Scope] Equation (10) proves strict pin exclusion at every finite size in this specified model. It neither bounds the contact uniformly below in the thermodynamic limit nor treats nonzero port, auxiliary sectors, or grammar excursions. An order-J thermodynamic pin gap would follow from a uniform positive contact lower bound in the pinned ground states; that additional bound has not been proved here. A tiny finite-size gap cannot be promoted to a uniform gap by Perron-Frobenius alone.

[Scope] Removing the terminal L_1 changes its charge by two. This comparison minimizes over Q, as required in C88's definition of E_good. It is not a theorem at artificially fixed Q. For identical cells the enlarged rooted clock ground state is positive and invariant under cyclic cell relabeling, so the physical zero-momentum subspace retains its minimum. The artificial distinction of a root does not create a lower unphysical minimum in this sign-free model.

## 4. Dressed cost through second order, with quantum sources

[Proved, finite volume] Set t=1 temporarily. Let |0> be the free bath ground state, P=|0><0| tensor 1, Q=1-P, and V_mat=sum_j A_j h_j. Degenerate perturbation theory in J gives

    H_eff=E_gap+J a_L H_mat-J^2 B+O_(N,V)(J^3),
    B=sum_ij chi_ij h_i h_j,
    chi_ij=<0|(A_i-a_L)Q(H_gap-E_gap)^-1 Q
                    (A_j-a_L)|0>.                       (11)

Indeed P V_mat P=a_L H_mat and the second term is -P V_mat Q(H_gap-E_gap)^-1 Q V_mat P. The contact matrix chi is positive semidefinite. No replacement of h_i by its expectation is made. In a real clock ground vector phi, its coefficient is

    B_phi=sum_ij chi_ij <h_i phi,h_j phi> >=0.            (12)

The clock ground vector is nondegenerate within each spectator-bit block, so its energy expansion follows by diagonalizing the first-order operator before evaluating (12). Restoring t, define chi using the unit-t bath and obtain

    E_sigma-E_good = J a_L (E_sigma^mat-E_good^mat)
                    -(J^2/t)(B_sigma-B_good)
                    +O_(N,V)(J^3/t^2).                 (13)

[Scope] These are fixed-volume expansions. The bath excitation gap vanishes with size. The constants in the last term are not volume-uniform, and (13) does not itself prove an isolated thermodynamic defect energy. Individual B's are nonnegative, while their difference has either sign. An extensive bound on either B separately does not control that difference for a single local defect.

[Numerical] second_order.py solves the projected bath resolvent and retains the full Gram matrix of quantum clock sources. The largest linear-solve residual is 2.49e-13. Representative coefficients are

| Cells compared | N,V | B_defect-B_good | Measured residual divided by J^2, smallest J |
| :-- | :-- | --: | --: |
| One pin in (3,3) | 6,6 | -0.0701028652 | +0.0700498682 |
| (2,4) against (3,3) | 6,6 | -0.0184402179 | +0.0184257622 |
| One pin in (15) | 15,2 | -0.416266463 | +0.412345896 |
| One pin in (15,15) | 30,1 | -0.123230510 | +0.127850190 |
| (14,16) against (15,15) | 30,1 | -0.000310860430 | +0.000330653639 |
| (30) against (15,15) | 30,1 | +0.00904218905 | -0.0106919894 |

The smallest J/t is .01 for the first three rows and .005 for the last three. The finite-J discrepancies are consistent with higher-order terms; no remainder certificate is inferred from these values. The length-30 row explicitly has an unfavorable dressing correction, while its total defect cost stays positive.

## 5. What the vacancy infrared argument supplies

[Inherited derived continuum input] The collaborator's DEFECT_TRANSFER.md obtains K_gap=rho^-2 and v_gap=2t rho sin(pi rho) for the phase model. At half occupancy K_gap=4 and v_gap=t. A tagged material contact has a smooth density term of dimension one, dimension-two descendants, and a first harmonic of dimension four. The noncontact part of the static response of the smooth term vanishes because

    integral_0^infinity (x^2-v^2 tau^2)/(x^2+v^2 tau^2)^2
                         d tau=0, x!=0.                 (14)

Its primitive is tau/(x^2+v^2 tau^2), which vanishes at both endpoints. Taking absolute values before the time integral would destroy this cancellation. The remaining descendant and harmonic static tails scale at worst as |x|^-3 and |x|^-7. Their first spatial moments converge.

[Inherited derived finite-density constraint] Canonical fixed vacancy number adds a rank-one term. For a scalar source s_j, the smooth response has the form C(delta_ij-1/N)/t, with C=0.07857579749... at half occupancy. The rank-one energy contribution is +(J^2/t)C(sum_j s_j)^2/N. Changing a bounded source on W sites changes it by O(J^2 W/t), since

    [(S+delta S)^2-S^2]/N=2(S/N)delta S+(delta S)^2/N,
    |S|=O(N), |delta S|=O(W), W<=N.                     (15)

For a bounded scalar source, summability of the remaining kernel's first spatial moment also controls interfaces between translated selected patterns. At distance r only O(r) pairs cross each interface. This justifies the weight W at quadratic order for that scalar comparison.

[Proved, finite-cell cancellation] If h_(c,j) are a cell's distributed bond contributions, sum_j h_(c,j)=H_c. Define f_j=h_(c,j)-H_c/m and B_j=sum_(r<=j) f_r. Since B_0=B_m=0, summation by parts gives

    sum_j A_j h_(c,j)=Abar_c H_c
                    +sum_(j=1)^(m-1)(A_j-A_(j+1))B_j.   (16)

Off-diagonal matrix elements between different isolated cell energies therefore couple to a contact difference. The smooth low-momentum matrix element acquires a factor q. This improves its infrared power counting. It is not a bound on the response for arbitrary frequencies, arbitrary quantum histories, or arbitrarily long cells.

[Inherited proved UV result] The collaborator's SINGLE_EXPONENTIAL.md provides a local Gaussian-filter normal form with first-spatial-moment remainder

    ||R_UV||_1 <= C t s^2 delta^-10
                   +C J delta^-5 exp(-c/delta^2), s=J/t. (17)

Taking delta=s^a, 0<a<1/10, makes this o(J). This controls an interaction remainder. The retained low-frequency operator is a filter in energy transfer, not a projection onto low absolute bath energy. Its quantum-cell ground-state comparison remains open in that work.

[Assumed, missing input] To use (9) in the microscopic mobile chain one needs, after common-bulk subtraction and with the cyclic constraint accounted for,

    R_sigma(J)-R_good(J) >= -J epsilon_ch(s) W,
    epsilon_ch(s) -> 0 as s -> 0,                        (18)

uniformly over defect patterns and unbounded wrong lengths. Here R is the difference between mobile energy and E_gap+J a_L E^mat, with the finite-K terms treated consistently. A length-weighted two-sided version is sufficient too. Equation (18) must apply to operator-valued sources, retain the slow payload dynamics, and have a finite-defect version without an uncontrolled negative boundary constant. No derivation of (18) follows from (14), (16), or (17). In particular, replacing the operator sources by static scalars is not a proof of it.

[Heuristic] The finite static contact response suggests epsilon_ch(s)=C s[1+|log s|], with C depending on fixed material data. The long-cell weight prevents the known false absolute per-cell estimate. This is a plausible scaling form, not a computed or bounded C for C88. The supplied analyses and the finite ED do not make it a controlled all-pattern estimate at a specified J.

## 6. The lambda-squared band and finite-K scales

[Proved in the isolated selected-cell reduction] Write V_port for the port with lambda removed. Its norm on a pair is sqrt(2). The vacuum code has P_code V_port P_code=0. With a uniformly averaged contact, the leading port operator is

    -J a_L lambda^2 P_code V_port Q
          (H_aux-E_code)^-1 Q V_port P_code.             (19)

There are two powers of J a_L in the vertices and one in the energy denominator, leaving J a_L lambda^2, not J^2 lambda^2/t. The selected carrier has diagonal 4 on the D rail and 1004 on the X and M rail. Its row-sum bound gives a single-excitation gap at least 2-|g|, and a port pair has gap at least 4-2|g|>3 at the schedule. Consequently this isolated second-order bond has norm below (2/3)J a_L lambda^2. Endpoint projections can improve this bound but are not needed here.

[Inherited proved, scope] The nontrivial thirteen-flavor dispersion has the much smaller scale J a_L lambda^2 E_unit. The finite material certificate bounds the normalized g and K errors and the isolated lambda^4 remainder. It does not prove a mobile many-cell effective interaction with the same norm. In (9) the exact frozen common bulk already includes the band ground energy. On unmatched boundaries the conservative full port norm B_port covers its effects without assuming a uniform auxiliary gap in every long wrong cell. Thus (9) never discards the lambda-squared energy, although it uses the weaker O(|lambda|) bound there.

[Numerical arithmetic, not floating-point simulation of tiny couplings] The schedule has

    log2(E_unit) = -858.577410080...,
    log2(lambda^2 E_unit) = -2074.577410080...,
    log2(f) = -2194.577410080....                         (20)

These logarithms avoid underflow. All explicit port budgets are negligible relative to gamma>2^-17. That scale separation says nothing about an unspecified C in (18).

[Inherited proved] C88 gives the complete leading grammar correction

    D_2=-(t^2 C_col+J^2 A_tail+2J^2 lambda^2 Z)/(JK),
    ||D_2||_local <= (2t^2+3J^2+6J^2 lambda^2)/(JK),
    ||R_2||_local,exponential <= C_3 nu^3/(JK)^2.         (21)

C_col counts material R_0 or E followed by a vacancy, A_tail counts the specified wrong-predecessor vacuum move, and Z counts wrong-end port channels. The color counter depends on bits and clock position. It is not a common scalar and has no lambda^2 factor. The mixed tJ/(JK) term vanishes. One may take nu<=2^13 J+4t+2|mu|.

[Proved bookkeeping; assumed use in a sector comparison] Dividing the local strengths by J gives

    eta_K,local <= 2/(s^2 K)+3/K+6lambda^2/K
                   +C_3 (nu/t)^3/(s^3 K^2).             (22)

A comparison coefficient C_geo is needed to convert this local interaction control into a common-bulk defect estimate. For explicitly matched localized interaction differences it comes from the support count and first spatial moment. Such a matching has not been constructed for the entire low-frequency microscopic theory. It is not legitimate to insert N ||D_2|| in place of W ||D_2|| when excluding finitely many defects.

[Proved scale check] Even at K=2^2560 the color term divided by J is 2^(1-2560)/s^2. Making J arbitrarily smaller at fixed K eventually makes this larger than the selector margin. For example s=2^-2000 makes it 2^1441. Increasing K after fixing J is mandatory. The band tolerance from C88 is still smaller, being proportional to lambda^2 E_unit times its root margin and matrix normalization. A defect-only K bound cannot certify that band.

## 7. A sufficient comparison and why it is still conditional

[Proved conditional algebra] Suppose (18) holds, and suppose finite-K matching gives a defect remainder bounded by J eta_K W. Use the fully frozen comparison (9), so finite-K terms already included in the frozen reference are subtracted once. Bound any additional matched terms conservatively by (22) times C_geo. On compatible half-density rings, with no unresolved gluing or boundary term, one obtains

    E_sigma-E_good >= J W
       [a_L(gamma-5sqrt(2)|lambda|)-epsilon_ch(s)-eta_K]. (23)

In the thermodynamic half-density limit the conservative value a gamma is 2.26865869686e-6. Sufficient conditions are

    5sqrt(2)|lambda| < gamma/4,
    epsilon_ch(s) < a gamma/4,
    eta_K < a gamma/4.                                  (24)

They leave at least J a gamma W/4. The first inequality holds by the explicit schedule. If (18) and the finite-K defect matching were known with finite constants, J could be chosen first and K increased afterward to meet the other two. This would establish the requested positive comparison. Those constants and the low-frequency estimate have not been supplied or obtained here.

[Scope] A residual negative boundary term O(J), or a cyclic flux allowance O(t/L), cannot be dropped in a claim about every finite-number defect. One must either cancel it between the compared sectors or check that it is smaller than the positive margin. An inequality with such an uncontrolled term proves neither the requested finite-defect statement nor positivity at all compatible sizes. The pin theorem (10) avoids this issue in its stated logical model only.

[Heuristic density check] In the projected selected reference the energy per physical site is

    e(rho)=-(2t/pi)sin(pi rho)
              +J e_b [rho^2-sin^2(pi rho)/pi^2]+... .    (25)

Thus tuning half occupancy gives mu=J e_b+... and e''(1/2)=2pi t+4J e_b+..., positive at sufficiently small J/t. The half-filled grand potential is -2t/pi+O(J), below the empty sector. The omitted terms include mobile matching and the band. These reference thermodynamic checks do not prove global density selection for the microscopic chain.

## 8. Nonperturbative exclusion of a whole giant cell

[Proved, marker-free legal sector] A legal cyclic material word without markers consists entirely of L letters or entirely of R letters. A head D or X would change the automaton state from zero to one without any marker to reset it. The same obstruction forbids an L-to-R transition in a marker-free cycle. Every vacuum rewrite on an all-R word would create an illegal predecessor junction, so its projected matrix element vanishes. The remaining material energy is the nonnegative field. Taking all bits zero attains zero material energy. Therefore E_marker-free=E_gap in the legal projected model. The selected trial below has lower energy by J a_L N |e_b| for every compatible N>=15 with a_L>0. This conclusion allows nonzero port because no port can act without a marker.

[Proved, legal projection] In the one-marker legal sector a port cannot act: it would require both ends of the same cell to be at their incompatible vacuum endpoints. All diagonal contributions are nonnegative. The remaining material transition matrix has absolute row sum at most 3+|t_d|+|g|<4 at the scheduled parameters, by counting the two rail steps, the switch and the possible conversion or direct endpoint link. Multiplication by contacts A_j cannot increase this row sum. Therefore, at K=infinity,

    E_one-marker >= E_gap-4J.                            (26)

This bound includes pin configurations and auxiliary states in that legal sector. It uses no cell excitation gap.

[Proved variational upper bound] On compatible N=15M, use the product of identical empty-bit vacuum clocks and the free bath state. Its port expectation is zero. The trial can be symmetrized over equivalent marker roots, since both the identical clock product and the free bath amplitudes respect their cyclic relabeling. It is an allowed zero-momentum state of the physical ring. Its expectation is E_gap+J a_L N e_b. Consequently

    E_one-marker-E_good >= J[a_L |e_b| N-4].              (27)

At half occupancy, a_L>.29 for N>=15, as follows from sin x>=x-x^3/6 at x=pi/(2N). With |e_b|>.0048254, every compatible N>=2865 makes (27) positive. This is an exclusion of an entire giant-cell ring in the legal projected model at all J,t>0. It is not a proof for a long cell embedded in selected bulk, nor for finite-K illegal excursions. The latter situations need the common-bulk matching described above.

## 9. ED checks and reproducibility

[Numerical] check_defects.py builds the exact gap compositions and the clock matrices of Section 3. It independently checks the free energy and mean contact in (5). The checks comprise 24 quantum-clock comparisons and 24 scalar contact comparisons. Every defect cost is positive. The maximum clock dimension is 6750, the maximum scalar dimension is 6435, and the largest eigenpair residual is 8.89e-12. No DMRG or large computation was needed.

| Quantum-clock comparison | N,V | J/t | Cost | Cost divided by first order |
| :-- | :-- | --: | --: | --: |
| One pin in (3,3) | 6,6 | .01 | .004145896761 | 1.00169248 |
| (2,4) against (3,3) | 6,6 | .01 | .000617587216 | 1.00299244 |
| One pin in (15) | 15,2 | .01 | .017051901841 | 1.00242404 |
| One pin in (15,15) | 30,1 | .005 | .009484120632 | 1.00033712 |
| (14,16) against (15,15) | 30,1 | .005 | .000004356917 | 1.00190090 |
| (30) against (15,15) | 30,1 | .005 | .000213435625 | .99874920 |

[Scope] The length-15 quantum-clock tests have high occupied density, not half density. The half-density quantum tests use b=3 and its own midpoint selector. Both retain the distributed marker energies and quantum contact-dependent clock motion. The half-density scalar tests at N=6,8 use e_b sum_j A_j as a common background and insert the actual b=15 bare defect costs as scalar contact sources. They test bath relaxation only. They do not reproduce the spatial profile, quantum clock, or payload band of a length-15 defect. In particular a scalar source labeled wrong_60 does not contain a sixty-site cell.

[Numerical] At N=8,L=16,J/t=.01 the scalar pin cost is .005747523498, which is .99414028 times first order. The scalar (14,16) cost is .000002651795, which is 1.00000699 times first order. This difference from the clock corrections is another reason not to identify scalar and quantum matching coefficients. The scalar potential is nonnegative, so its positive sign also follows from operator monotonicity; the useful numerical information is its resolved relaxation magnitude.

[Scope] The numerical calculations set lambda=f=0 and K=infinity. They do not resolve the scheduled lambda^2 band or finite-K terms. Those scales are retained analytically in Sections 6 and 7. No central charge or microscopic thermodynamic pin gap is measured.

[Reproduction] Run, from this folder,

    python3 -B check_defects.py
    python3 -B summarize_checks.py
    python3 -B second_order.py

The files ed_results.json, checks_summary.json and second_order_results.json retain parameters, energies, dimensions and residuals. All scripts use only Python, NumPy and SciPy. All outputs remain in this folder.

## 10. Conclusion and provenance

[Result] Pin exclusion is proved for the finite zero-port logical vacuum-clock model. The computed finite-volume quantum dressing, the length-weighted reference margins, and the legal giant-cell bound support the proposed selected sector. The full C88 S_pin remains open because a one-sided low-frequency operator-valued defect comparison is missing. No claim of a completed d=13 mobile ground-state theorem is made.

[Provenance] The starting inputs were ../../specs/_simplest_background.md, C88, C85, C83, the collaborator's counterexample campaign and DEFECT_TRANSFER.md, frozen/DEFECT_PRESSURE.md, direct_attack_20261028/vacancy, and matching_attack/README.md. The latter two explicitly retain the gap that appears as (18). No network commands, subagents, or writes outside work_C89 were used. The project-level memory remains untouched under the user's folder restriction.
