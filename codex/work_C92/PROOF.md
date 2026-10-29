# C92. Mobile defect traces

[Status] The full finite-penalty thermal comparison is not proved. There is a new uniform contact bound in a positive length-tilt ensemble and a resulting activated pin comparison for the legally projected d=13 chain. The remaining wrong-length comparison is stated as a specific, common-bulk-subtracted thermal response estimate. A separate relative penalty estimate is needed to pass to the full finite-penalty ring. Neither is established by the supplied frozen certificates.

[Convention] A tag governs the following argument until the next tag. Proved means derived here. Imported, proved means a result used with its stated scope from the supplied sources. Assumed means a missing hypothesis, including hypotheses of conditional theorems. Numerical means a finite matrix check. No subjective likelihood estimates are assigned.

## 1. The trace that must be compared

[Imported, proved] C88 and C90 give the legal projected Hamiltonian in ordered material coordinates. With N material sites and gaps g_j in the nonnegative integers, it is

    H_sigma(s) = H_gap + s sum_j g_j + J sum_j A_j h_sigma,j,
    H_gap = -t sum_j (E_j^+ E_(j+1)^- + E_j^- E_(j+1)^+),
    A_j = |0><0|_j,   E_j^+|g> = |g+1>,   E_j^- = (E_j^+)^*.
                                                               (1)

Here sigma records ordered cell lengths and pin flags. The internal space is finite at each N, even when cell lengths vary. All internal states in that sector are retained. The material operator includes the distributed onsite halves. In particular the marker energy is a contact term. The shifts have norm one. The projection onto legal words is independent of the gaps after the unique-color isometry.

[Definition] We work first with independently labeled periodic gaps, N>=3, and sum over all vacancy numbers with s>0. A term sN or a material chemical potential is a common scalar at fixed N and can be restored later. Set

    Z_sigma(T,s) = Tr exp[-H_sigma(s)/T],
    F_sigma(T,s) = -T log Z_sigma(T,s).

All couplings, including s in a compact positive interval, are fixed before N tends to infinity and then T tends to zero. The gap trace exists because the hopping and the material interaction are bounded at fixed N, whereas s sum g_j has a geometric heat trace. The comparison is not a fixed-vacancy-number variational argument.

[Imported, proved] C58's root compression and bounded-seam argument, as transferred in C90, identifies the rooted thermodynamic pressures with physical ones at positive density. It does not identify finite-ring spectra. The arguments below may cut one seam for the whole ring. Its energy is O(t+J M), independent of the number of defects, and vanishes per site after the thermodynamic limit. Here M is a uniform bound on the compressed material bond norm. Cutting a charge seam at every defect would instead cost O(t) per defect and would invalidate a small-J selector estimate.

At finite K the corresponding physical seam constant can also depend on JK. It is still a fixed global constant. Throughout, bounds stated for independently labeled gaps are finite-volume statements in that representation. Their physical interpretation below is at the level of pressure after removing this seam and the root multiplicity.

[Definition] At finite K the selected sector means the exact invariant sector with all material cell lengths b and all pin flags zero, retaining its entire Hamiltonian and heat trace. It includes virtual grammar excursions and thermally excited illegal words with those same conserved labels. It is not the compression P H P onto legal words. The latter is not invariant at finite K. Comparison to a separately matched low-energy trace would require an additional, exact massive thermal matching result.

For C88 the conserved labels are independent of whether a word is legal. Thus every full word with markers belongs to one ordered-length and pin sector, including all its illegal internal configurations. Words with no markers and words with no material sites are separate classes. This observation permits a complete sector sum without pretending that instantaneous illegal-bond locations are conserved.

## 2. Uniform contact expectation from vacancy removal

[Proved, contact lemma] Let (1) act on a gap tensor product and an arbitrary finite internal space. The operators h_sigma,j may fail to commute with one another. Suppose

    ||h_sigma,j|| <= M,    s in [s_min,s_max],
    c = s_min - J M > 0.

For its Gibbs state rho, at every site j, every N>=3, every internal sector and every

    0 < T <= c/(4 log 2),

one has

    Tr rho A_j >= a_* = min(1/2, c^2/(256 t^2)).              (2)

Take a_*=1/2 if t=0. The same bound holds in zero-temperature limits of these Gibbs states. It is independent of cell lengths, pin density, internal excitation gaps and N.

Proof. Fix j, write P=A_j, Q=1-P, S=E_j^+, and P_1=|1><1|_j. The identities S^*S=1 and SS^*=Q imply that

    Phi(rho) = P rho P + S^* rho S                            (3)

is a trace-preserving channel. It removes one vacancy from the nonempty part of the chosen gap and leaves the empty part unchanged. Its Stinespring isometry has a two-dimensional environment. Subadditivity of entropy gives

    S(rho) - S(Phi(rho)) <= log 2.                            (4)

This bound allows arbitrary entanglement with all other gaps and the material system. There is no replacement of a quantum material source by a scalar.

On the range of S, the compression identities are

    S^* H_sigma(s) S = H_sigma(s) + s - J P h_sigma,j,
    S H_sigma(s) S^* = Q H_sigma(s) Q - sQ + J P_1 h_sigma,j.
                                                               (5)

Indeed S^* E_j^+ S=E_j^+ and S^* E_j^- S=E_j^-, while S^*g_j S=g_j+1 and S^*P S=0. All other material contact terms commute with S. Only the two gap-hopping bonds meeting j connect P and Q, so

    ||P H_sigma(s) Q|| <= 2t.                               (6)

Put a=Tr rho P. Positivity of rho and the Hilbert-Schmidt Cauchy-Schwarz inequality give

    |Tr rho(PHQ+QHP)| <= 2||PHQ|| sqrt(a(1-a)).

Using P_1<=Q in (5), the energy change under (3) obeys

    Tr H(Phi(rho)-rho)
       <= -(s-JM)(1-a) + 4t sqrt(a(1-a)).                    (7)

The Gibbs variational principle, (4) and (7) imply

    0 <= F_H(Phi(rho))-F_H(rho)
      <= -c(1-a) + 4t sqrt(a(1-a)) + T log 2,                (8)

where F_H(omega)=Tr H omega-T S(omega). If a<=1/2, then c/2<=4t sqrt(a)+c/4, so a>=c^2/(256t^2). If a>1/2 the other branch of (2) applies. For t=0, (8) excludes a<=1/2 under the stated temperature restriction. This proves (2).

All entropies and energies used above are finite. One can first restrict to sum g_j<=V_max, which the channel maps into itself, apply the variational argument, and then remove the cutoff. The hopping preserves total gap occupation, and the linear tilt supplies convergence of both the trace and its first moment. No cutoff depending on N or T is retained in (2).

[Proved, concrete interval] On

    t/2 <= s <= 3t/4,    JM <= t/4,

one can take

    a_*=1/4096,    T<=t/(16 log 2).                          (9)

The interval contains the free half-density tilt 2t/pi. The bound is deliberately crude and does not assert proximity to the free contact value 1/2-2/pi^2. For the legal C88 model M is finite and independent of K, N and the cell lengths, so the condition on J can be met before any finite-penalty choice.

For an explicit conservative C88 choice, M=2^11 is sufficient. The scheduled diagonal u has maximum below 1006. There are 24 directed rewrite entries before their Hermitian conjugates, each with magnitude at most one, so ||w||<=48 by the triangle inequality. Thus ||b_11||<1054<2^11, and compression cannot increase the norm. Taking J/t<=2^-13 verifies JM<=t/4. This is a sufficient range for the contact and projected pin theorems only. It does not certify the missing mobile selector or infrared estimates at that coupling.

[Scope] Removing a vacancy changes total physical length. Therefore (3) is admissible in the tilted ensemble and is not admissible in a fixed-L, fixed-N variational problem. A uniform canonical contact bound has not been proved here. Uniform pressure bounds on a tilt interval can nevertheless be transferred to canonical free energies by the C58 duality argument, with interior optimizing tilts. These are distinct statements.

## 3. A thermal pin replacement

[Imported, proved local input] The C88 rule table and Part IV of the C88 referee report imply the following replacement property at K=infinity. A pinned b-cell has only L letters followed by E, M_A or M_B, and has no intracell kinetic move. Replace its entire internal state by a zero-bit vacuum clock state on b positions. The marker chemical energy is unchanged, the nonnegative bit and charged-marker energies do not increase, and the clock hopping expectation added is

    -2J sum_(r=0)^(b-2) phi_r phi_(r+1) Tr rho A_(s_c+r).
                                                               (10)

The replacement preserves the complete charge-gap marginal and all internal marginals outside the replaced cells. The positive charged-marker terms remain contact-weighted before their removal. This is why replacing them is an energy decrease, not an unproved statement about a dressed gap.

[Proved, entropy extension] Fix any ordered cell lengths and any specified set P of k pinned b-cells. Other cells, including wrong-length cells, may be arbitrary. Let rho be the full Gibbs state of this sector. Trace out the k pinned internal spaces and tensor in a product of clock vectors. If D_pin bounds the dimension of one pinned cell, subadditivity gives

    S(rho)-S(rho_trial) <= k log D_pin.                       (11)

The bound is log D_pin per cell, including correlations between pinned cells and the bath. No claim of an isospectral replacement is needed. For C88,

    D_pin = 3*2^(b-2),    b=15,    D_pin=24576.              (12)

At zero port, optimizing each clock vector in (10) gives the thermal version of the referee bound,

    F_P - F_unpinned
       >= J sum_(c in P) lambda_max(Path_c) - T k log D_pin.
                                                               (13)

Path_c is the b-vertex path with edge weights Tr rho A_(s_c+r). The expectations in this expression are taken in the pinned Gibbs state, not in a free bath or an unpinned ground state.

To prove (13), apply the Gibbs variational principle in the target sector to rho_trial, use (10) for its energy, and use (11) for its entropy. Every target state lies in the same length sector. Payload charge is summed, as it is in the stated full-pressure problem. The map need not preserve an artificially fixed payload charge.

[Proved, nonzero port allowance] Write B_port=sqrt(2)|lambda| for the material port norm. At most 2k port bonds meet the reset cells. Bounding the difference of their expectations in the two states by twice their norm costs at most 4J B_port k. This bound does not require a sign of the port or a positivity property of its matrix elements.

Choose the normalized positive sine vector on each b-vertex path. By (2), its energy gain is at least 2J a_* cos(pi/(b+1)). Thus

    F_P-F_unpinned >= k[delta_pin - T log D_pin],
    delta_pin = J[2a_* cos(pi/(b+1))-4B_port].                (14)

This is uniform in the other cell lengths and all other quantum states. At the C88 port value and with (9),

    delta_pin >= J/4096 > 0.                                (15)

The inequality follows from 2cos(pi/16)>1.96 and the stated lambda=2^-608. It is a conservative analytic estimate, not a claim to reproduce the reference pin cost near 0.58J.

[Proved, pin trace sum] For an ordered length pattern with m_b cells of length b, sum (14) over the subsets of those cells that are pinned. Pins in wrong-length cells can remain part of the wrong-cell class. The same target trace, with every b-cell unpinned, is used for every subset. Consequently

    1 <= Z_all_b_pin_flags/Z_no_b_pins
       <= [1+D_pin exp(-delta_pin/T)]^m_b.                  (16)

The pressure error per material site is at most

    (T/b) log[1+D_pin exp(-delta_pin/T)],                    (17)

since m_b<=N/b. This is o(T^2), uniform in the other length pattern. It also remains true after any positive sum over those patterns. A single global seam changes only a subextensive logarithm. It must not be charged once per pin.

[Scope] Equations (14)-(17) prove the pin part of the tilted thermal comparison for the exact legal projected C88 chain, including its auxiliary states and nonzero port. They do not prove a finite-K statement. They are stronger than a fixed-volume ground-state inequality, but have a different ensemble scope from a canonical contact assertion.

## 4. Entire marker-free and one-marker sectors

[Proved, legal marker-free trace] In C88 a legal cyclic material word without markers is all L or all R. Every projected material transition vanishes and the remaining field is nonnegative. Its internal dimension is at most 2^(N+1). If F_gap is the free energy of H_gap+s sum g_j, min-max comparison gives

    Z_mf <= 2^(N+1) exp(-F_gap/T).                           (18)

For N divisible by b use the Gibbs state of the free gap Hamiltonian and a product of the selected zero-bit vacuum clocks as a trial state in the selected sector. The free gap Gibbs state is invariant under material-index translations. Put a_0(T,s)=Tr rho_gap A_j. The contact lemma with M=0, or the weaker common bound (9), gives a_0>=a_* at the temperatures under consideration. If E_b=b e_b<0 is the inherited frozen vacuum-cell energy, the trial gives

    F_good <= F_gap + J a_0 N e_b,
    Z_mf/Z_good <= 2 exp[-N(J a_* |e_b|-T log 2)/T].        (19)

The port expectation in this trial is zero, so nonzero lambda does not alter the argument. The C88 certificate supplies e_b<-0.0048254. Thus the entire legal marker-free sector has strictly smaller tilted pressure at sufficiently small T. No estimate of the interacting selected pressure or of its T^2 coefficient was used.

[Imported, proved, then thermal consequence] C89 proves H_one_marker>=H_gap-4J on the legal one-marker sector, including auxiliaries and pins, with no assumed long-cell excitation gap. Its internal dimension is at most 11^N. Hence

    Z_one_marker/Z_good
      <= exp[4J/T - N(J a_* |e_b|-T log 11)/T].              (20)

It too has smaller pressure at sufficiently small T. Its vanishing clock level spacing is irrelevant to this estimate.

[Scope] Neither (19) nor (20) is an estimate for a marker-free fragment or a giant cell embedded in an infinite selected background. In an island comparison the untouched charge fluid must remain common. Cutting it out at every island loses an O(t) boundary energy. Those islands are included in the weight and the missing comparison below.

[Proved, all-vacancy class] At infinite penalty there are q constant-color all-vacancy words, with q=2 for C88, and no hopping or material energy. At finite penalty the vacancy-only Hamiltonian is the classical cost JK for each change of color. Its ring trace satisfies

    Z_vac <= q[1+(q-1)exp(-JK/T)]^L.                        (21)

This follows from the q by q transfer matrix with diagonal entries one and off-diagonal entries exp(-JK/T), whose largest eigenvalue is the bracket in (21). For the one-vacancy d=24 lift, q=1. At a chemical potential near the half-density value the selected free-charge trial has physical energy density -2t/pi+O(J)<0, while (21) has pressure tending to zero. Thus the all-vacancy class cannot win the physical pressure at sufficiently small T. This conclusion is uniform for, for example, |mu|<=t/4 and sufficiently small JM/t. It does not settle nonempty marker-free sectors at finite K.

The same argument disposes of sequences with N=o(L). At fixed N the hopping and negative material terms have norm at most C N, with C finite at fixed parameters. Keep only penalties on vacancy-vacancy bonds. For any occupancy pattern the vacancy intervals contribute at most q^(N+1) lambda_vac^(L-N), where lambda_vac=1+(q-1)exp(-JK/T). Internal material labels contribute at most d_mat^N and occupancy patterns at most binom(L,N). Min-max therefore bounds the fixed-N trace by these factors times exp(CN/T), including the material chemical potential in C. For N/L tending to zero their logarithms per L, except log lambda_vac, vanish. These sectors have at most the pressure of (21), strictly below the selected half-density trial at small T.

## 5. The required wrong-length charge estimate

[Imported, proved] The C88 material selector certificate gives

    E_ell = mu_sel - 2cos(pi/(ell+1)),
    e_b=E_15/15<0,
    E_ell-ell e_b >= gamma ell,  ell!=15,
    gamma>2^-17.                                           (22)

Its conservative finite-penalty certificate can be used for gamma; the legal projection does not worsen that isolated-cell lower bound. At zero port E_pin,b-E_b=2cos(pi/16). Long-cell cost grows as |e_b|ell+O(1). The subtraction is ell e_b, not E_b.

[Definition] For a sector with markers set

    W(sigma) = sum_(ell_c!=b) ell_c
                  + number of pinned b-cells.               (23)

A wrong-length pinned cell is counted only in the first sum. Short illegal fragments within a cell are not given a fictitious conserved label. At finite K all their states belong to the full trace of that cell-length sector. If a separate defect-resolved representation is used, its fragments must receive their material lengths in addition to any wall weight.

[Proved, frozen thermal comparison] For the legal frozen material chain use the same selected-sector thermal bulk on both sides. Let f_mat,g(T;u a_0) denote its exact free energy per material site when the material Hamiltonian is multiplied by u a_0. A selected open run differs from its limiting bulk free energy by at most one port norm. A ring differs by at most two. This follows by concatenation and min-max, exactly as in the supplied frozen DEFECT_PRESSURE proof.

For a length pattern with k wrong cells, delete at most 2k ports touching them. The remaining selected runs contribute at most k further port norms. Bound each wrong cell trace by its full dimension, at most 11^ell, times its ground-energy exponential. The selected product-clock trial gives f_mat,g<=u a_0 e_b. Consequently, for b-compatible N and no b-pins,

    F_sigma^mat(T;u a_0)-F_good^mat(T;u a_0)
       >= u a_0(gamma-3B_port) W
                         -T(log 11)W - B_boundary,          (24)

where B_boundary can be chosen 2u a_0 B_port for the periodic comparison. Boundary conventions and a finite number of root offsets can instead be included as fixed entropy factors. Only one global boundary allowance remains. No cell excitation gap was used in this argument. Pinned wrong-length states obey the same cell lower bound and are included in the full dimension bound.

For later summation one may replace 3B_port by the more conservative 5B_port. At the C88 schedule this change is immaterial. The proof at this stage concerns the legal material chain; it is not an assertion that forbidden-arrow positions remain conserved in the full colored chain.

[Definition, exact common-bulk subtraction] For fixed T,s let a_0 be the free gap contact expectation above and define, for 0<=u<=J,

    H_sigma(u)=H_gap+s sum g_j+u V_sigma,
    V_sigma=sum_j A_j h_sigma,j,
    H_sigma^mat=sum_j h_sigma,j,
    R_sigma(u,T,s)
       =F_sigma^mob(u,T,s)-F_gap(T,s)
                            -F_sigma^mat(T;u a_0).          (25)

R_sigma(0,T,s)=0 for every internal sector, including its entire internal entropy. The comparison uses R_sigma-R_good, not separate extensive upper bounds on the two terms. In particular the extensive selected charge-relaxation energy is retained.

[Assumed, C_T, sufficient integrated form] A sufficient legal wrong-length input is

    R_sigma(J,T,s)-R_good(J,T,s)
       >= -J epsilon(J/t) W(sigma) - B_N(T,s),
    epsilon(x)->0 as x->0,    B_N(T,s)/N->0,                (26)

uniformly over all ordered length patterns, all N, s in the chosen source interval and all sufficiently small T. The subextensive allowance is global, not per defect. For pressure comparison it may be O(1) at fixed T. No absolute value is required. A two-sided bound with the same length weight is sufficient but stronger.

Since section 3 already removes b-pins, it is enough to prove (26) on sectors with no b-pins. Pinned wrong-length cells remain included among wrong cells. The version displayed with the full W is a uniform extension, not an additional pin hypothesis needed for the projected argument.

[Proved, exact response formulation] Here is an explicit quantum estimate that would imply (26). For a Hamiltonian H define

    chi_H(X)=integral_0^(1/T)
                <(X(tau)-<X>)(X-<X>)>_H d tau,
    X(tau)=exp(tau H) X exp(-tau H).

At finite N and T, differentiation of the trace gives

    dF_sigma^mob/du = <V_sigma>_mob,
    dF_sigma^mat(T;u a_0)/du = a_0 <H_sigma^mat>_mat,
    d^2R_sigma/du^2
       = -chi_mob(V_sigma)+a_0^2 chi_mat(H_sigma^mat).        (27)

These are Duhamel response identities, not equal-time variances. Translation symmetry of the free gap Gibbs state implies R_sigma'(0)=0 even for an inhomogeneous material pattern. Thus, writing

    D_sigma(u)=chi_mob,sigma(V_sigma)
                        -a_0^2 chi_mat,sigma(H_sigma^mat),

one has the exact identity

    R_sigma(J)-R_good(J)
       = -integral_0^J (J-u)[D_sigma(u)-D_good(u)] du.       (28)

The one-sided estimate

    D_sigma(u)-D_good(u) <= C W(sigma)/t                    (29)

uniformly in N,T,s,u and all wrong lengths would prove (26) with epsilon(J/t)=CJ/(2t) and zero boundary term. Only the weighted integral in (28) needs to be bounded; pointwise (29) is a convenient stronger sufficient condition. An integrated bound with logarithmic J/t dependence also suffices if epsilon tends to zero.

Every source in (27)-(29) is the full quantum operator, including clocks, auxiliaries and ports. The covariance subtraction includes the full frozen material response. A lower bound on <A_j> does not bound these connected, imaginary-time responses and does not prove (29).

[Status] Equation (29), or its weaker integrated version (26), has not been proved here. The contact lemma proves a lower bound, not a source-uniform closeness estimate to a_0 and not the comparison of these two coupled response functions. The gradient endpoint cancellation in matching_attack and C89 improves power counting for off-diagonal cell transitions. It does not give this uniform thermal estimate for arbitrary quantum histories and unbounded cells.

[Imported, proved obstruction] The matching_attack note exhibits a one-cell sector whose charge residual differs from the selected residual by a favorable extensive term of order J^2 N/t. It rules out an absolute bound proportional only to the number of long wrong cells. It does not contradict (26), (28) or (29), since (26) is one-sided and W=N for that cell. Subtracting another identical scalar bulk from both sectors cannot repair the false per-cell absolute estimate.

## 6. Illegal words at finite penalty

[Proved, exact accounting] Let F_sigma^K denote the full finite-K tilted free energy in an invariant ordered-length and pin sector, with all color words and all illegal material words retained. Define F_sigma^infinity by (1) and

    Q_sigma^K(T,s)=F_sigma^K(T,s)-F_sigma^infinity(T,s).     (30)

These quantities are exact, not truncated effective Hamiltonians. The comparison can sum over sigma without any classical decomposition by illegal-bond positions. Illegal configurations within the selected labels already occur in F_good^K on the denominator side. They are not counted again as independent defects.

[Assumed, K_T, relative penalty comparison] A sufficient finite-penalty estimate is

    Q_sigma^K-Q_good^K
        >= -J eta_K W(sigma)-B_N^K(T,s),
    B_N^K(T,s)/N->0,                                     (31)

uniformly in the same thermal and source windows. The comparison must also hold for the nonempty marker-free class with weight N, and for the whole one-marker sector when (20) is used. The small constant eta_K must fit below the selector margin after (26). For pins alone a much weaker tolerance fits below (15).

[Imported, proved local budget; unproved relative use] C89 and C90 record the leading C88 correction

    -(t^2 C_col+J^2 A_tail+2J^2 lambda^2 Z)/(JK).

Divided by J, a local bound including the next remainder is

    2/[(J/t)^2 K]+3/K+6lambda^2/K
                     +C_3(nu/t)^3/[(J/t)^3 K^2].           (32)

The term t^2/(JK) is not port-suppressed. A local estimate such as (32) bounds a total-volume error; it does not prove the relative W-weighted estimate (31). Converting it requires an exact common-bulk matching or another relative trace argument. Choose J before enlarging K, and verify any K-dependence of matching constants. Holding K fixed while making J arbitrarily smaller does not control (32).

[Proved, why an extra grammar activity is not automatic] In the colored chain N_bad does not commute with H. Thus summing instantaneous illegal positions as if they were invariant sectors is invalid. A finite-order unitary elimination retains a nonzero remainder. Dropping it gives an error independent of T, which cannot be called o(T^2) at fixed couplings. Keeping the full F_good^K and proving (31) avoids that mistake. Alternatively an exact defect-resolved trace representation must retain virtual processes in the common good trace and prove an activated bound for its remaining real excitations. Such a representation is not supplied by the frozen DEFECT_PRESSURE argument or by C58's invariant wrapper flags.

[Scope] With (31) all illegal words in marker-containing sectors are included, not discarded. Illegal fragments inside a long cell are paid by its length in W. Illegal fragments inside a selected-length unpinned cell occur on both sides through the exact selected sector. Whole marker-free words require their own N-weight comparison, and N=0 is handled by (21). If "defect-free" is intended to exclude even the virtual illegal components of a selected eigenstate, it does not define an invariant finite-K sector and is not the denominator used here.

## 7. Conditional full trace comparison and the length sum

[Proved, conditional on C_T and K_T] Combine (24), (26) and (31), or first remove pins with (14). At fixed parameters satisfying

    delta_W = J[a_*(gamma-5B_port)-epsilon(J/t)-eta_K] > 0,
    delta_P = delta_pin-J eta_K > 0,                        (33)

the full relative free energy of a sector with markers has a lower bound of the form

    F_sigma^K-F_good^K
       >= delta W(sigma) - T c_ent W(sigma)-B_N,
    delta=min(delta_W,delta_P)>0.                           (34)

Here c_ent is a finite constant depending on the fixed alphabet and b. It includes full internal dimensions of wrong cells and D_pin for pins. B_N/N tends to zero at each fixed T. One can apply (31) directly to each pin pattern and use (14) in its legal comparison; no assumption about the sign of a finite-K pin correction is made.

The sector count in material coordinates is exponential in the length weight, with constants independent of N. To see this, give each wrong cell its starting material index and its positive length. Give each b-pin its starting index and weight one. All remaining cells have length b and zero pin flag. Their marker offsets have at most b choices per remaining connected run. With k defects there are at most k runs, up to a single global offset. Thus the offsets and the finitely many defect types contribute at most A_0^W times a fixed global factor. The internal state multiplicities have already been included in c_ent. A configuration consisting entirely of wrong cells is also covered. Root rotations contribute at most a polynomial factor in N, not an extensive new entropy.

Let A absorb A_0 and exp(c_ent). For each possible anchor, sum all positive weights. Dropping overlap and total-length constraints increases the positive sum. Therefore a valid upper bound is

    z(T) = sum_(w>=1) A^w exp(-delta w/T)
         = A exp(-delta/T)/[1-A exp(-delta/T)],
    Z_markers^K/Z_good^K <= exp(B_N/T) poly(N) [1+z(T)]^N.
                                                               (35)

This is the required length-weighted remainder summation. In particular an arbitrarily long cell costs an exponentially decaying factor in its material length. A count-only factor would not dominate its internal entropy.

The displayed comparison uses N divisible by b so that its denominator exists. This restriction does not discard the other material-count residues in the full finite-K physical pressure. Given a physical ring with material count N, append r occupied sites with 0<=r<b and N+r divisible by b, cutting the one original joining bond. Compress the enlarged ring to fixed letters on those r sites. Only O(b) literal bonds differ from the original operator. Peierls compression and a norm bound give Z_(L,N)<=exp(C_b/T) Z_(L+r,N+r), with C_b finite at fixed K. The target trace here is the full trace at that compatible material count, to which (35) applies. Summing over the finitely many r and the O(L) counts supplies only subextensive factors. This argument is used at finite K, where an inserted word is allowed even if it is illegal. The independent projected pin estimate (16) needs no divisibility assumption.

If (31) also transfers (19) with a positive remaining density margin, the nonempty marker-free class contributes an additional term exp[-kappa N/T+O(N)] to the right side of (35). The all-vacancy class is separately dominated as in (21). Hence, at fixed s,

    0 <= p_full^K(T,s)-p_good^K(T,s)
       <= T log[1+z(T)] = O(T exp(-delta/T)) = o(T^2).       (36)

Here p is pressure per material site in the tilted representation. The global boundary and polynomial factors disappear as N tends to infinity at fixed T. This proof never requires a gap in the selected sector or differentiability of its thermal pressure.

[Proved, thermodynamic transfer under the stated interval condition] The uniform version of (36) passes to fixed interior physical densities using the C58 canonical duality and bounded seams. If q(T,s) is pressure per material site with s multiplying total physical length, its lengths are at least one, so q(T,s+h)<=q(T,s)-h for h>0. A uniform pressure error r(T) therefore moves a pressure zero by at most r(T). A supremum of two functions differing by at most r(T) differs by at most r(T). These two elementary bounds respectively give the physical summed-length and canonical transfers. For pressure per selected cell the lower length is b and the zero displacement is at most r(T)/b.

The optimizing tilts or pressure zeros must lie inside the interval on which the comparison is uniform. The positive-curvature free reference brackets half density, and C90 explains how a sufficiently small uniform norm bound preserves a finite-difference bracket for the selected sector. This placement is an explicit hypothesis of the use of (36), not a conclusion about arbitrary fixed sources. No derivative of an uncontrolled thermal remainder is taken.

All common scalar terms remain in both exact selected and full free energies. Since the error in (36) tends to zero as T tends to zero, their ground pressures coincide under the hypotheses. Subtracting that common ground pressure then preserves the same o(T^2) bound. There is no discarded O(J^2) bulk term.

## 8. Dimensions 13, 14 and 24

[Imported, proved data] The appropriate selected lengths and color numbers are

| Construction | d | b | Vacancy colors |
| :-- | --: | --: | --: |
| C88 pinless lift | 13 | 15 | 2 |
| C83 lift | 14 | 14 | 2 |
| C85 lift | 14 | 15 | 3 |
| Original route-A vacancy lift | 24 | 24 | 1 |

[Proved transfer of the new lemma] The contact proof uses only norm-one gap shifts, a positive linear length tilt and uniformly bounded internal contact sources. It therefore applies unchanged to the legal projected version of every row, with that row's material bond norm. For C88 the local pin replacement and D_pin are verified above. C85 and C83 have different pin penalties and alphabets, so their thermal pin statements must use their own replacement map if such flags are kept. No new d=14 pin spectrum is asserted here.

[Proved transfer of the conditional theorem] For each row, replace (22) by its inherited selector and fragment margins and replace the alphabet dimension in the entropy factors. The proof from (26), (31) and the all-sector margins to (36) is unchanged. The d=24 chain has no color-ambiguity correction, but it still needs the quantum thermal comparison (26) and exact finite-penalty relative control. Its frozen thermal theorem is not a proof of these mobile inputs. C83 and C85 need their own finite-penalty constants. No d=13 constant is silently assigned to them.

## 9. Checks and final ledger

[Numerical] Running python3 -B checks.py in this folder performs 48 contact-channel checks with noncommuting Pauli material sources, across two tilts, two couplings, two temperatures and two total-vacancy cutoffs. It checks the trial free-energy inequality, the entropy loss and the contact lower bound. The largest matrix dimension is 112. Sixteen vacuum-clock pin trace checks have largest dimension 336. All assertions pass. CHECKS.md records the residual summary. These checks do not establish the thermodynamic limits; the channel and variational proofs do.

[Proved consistency checks] At J=0 the contact bound stays positive on a fixed positive tilt interval. As s_min-JM decreases to zero its lower bound vanishes, so it makes no claim in a dilute or zero-tilt limit. At full occupancy every A_j=1 and the local pin energy gain reduces to 2J cos(pi/(b+1)), agreeing with the frozen clock. The entropy terms vanish as T tends to zero. Long-cell costs and errors both carry material length. All exponentials have energy divided by T. A single physical seam disappears before the T limit. The finite-K correction grows relative to J if J decreases at fixed K, as required by (32).

| Claim | Status | Scope |
| :-- | :-- | :-- |
| Uniform positive contact expectation | Proved | Legal projected, positive tilt, JM<s_min, all patterns and lengths |
| Activated removal of b=15 pins | Proved | C88 legal projected thermal trace, including auxiliaries and nonzero port |
| Entire legal marker-free and giant-cell traces | Proved | Tilted pressure, compatible N, sufficiently small T |
| All-vacancy sector | Proved | Dominated at the half-density source, including finite-K vacancy color walls |
| Wrong-cell thermal charge comparison | Assumed | Equation (26), sufficiently implied by the response bound (29) |
| Finite-K illegal-word comparison | Assumed | Exact relative estimate (31), not a discarded local remainder |
| Full-ring o(T^2) comparison | Proved conditional | Equations (33)-(36), source interval and exact full selected trace |

[Conclusion] The contact input needed by the pin replacement can be proved in the length-tilt ensemble. This closes the projected thermal pin step. The full mobile defect theorem still requires the common-bulk quantum response estimate (26) and the relative finite-penalty estimate (31). The supplied references and the new contact argument do not prove either. No unconditional full-pressure equality or central-charge conclusion is claimed.

[Provenance] Read the project memory, C89, C90, C58, the C88 referee report, the frozen DEFECT_PRESSURE proof and matching_attack/README.md. The local C88 rule table was consulted to check the pin dimension and conserved-sector definitions. All files and Python execution stay in work_C92. No network commands, delegated agents, external edits or large computations were used. The project STATE.md remains unchanged under the user's explicit restriction.
