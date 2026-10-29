# C94. The legal charge comparison

[Status] The literal pointwise estimate in C92 equation (29) is false with its stated uniformity in N, T, the length tilt s and 0 <= u <= J. This can be proved for the full legal C88 trace, including the auxiliary states and the scheduled port. The integrated comparison in C92 equation (26), which allows a global subextensive term, is not disproved and is not proved here.

[Scope] Every statement concerns the independently labeled periodic gap representation used in C92. A finite-ring response statement is not automatically a statement about the physical ring with a root constraint. The obstruction below applies directly to C92's stated representation. The numerical interacting models are the exact zero-port vacuum-clock blocks, with distributed marker energies. They do not simulate the full auxiliary trace or the scheduled payload band.

[Convention] Proved labels an argument supplied here, or an imported result with its source identified. Numerical labels floating-point evidence. Assumed identifies a hypothesis not established here. Heuristic labels an extrapolation. A tag governs the following argument until the next tag.

## 1. Definitions and the order of limits

[Proved, imported from C92] Write

    H_sigma(u,s) = H_gap + s V + u X_sigma,
    V = sum_i g_i,
    X_sigma = sum_i A_i h_sigma,i,
    A_i = |0><0|_i,
    Hmat_sigma = sum_i h_sigma,i.

Here V is the vacancy number, not the interaction. The free contact is a_0(T,s). Define

    R_sigma(u) = Fmob_sigma(u)-Fgap-Fmat_sigma(u a_0),
    D_sigma(u) = chi_mob_sigma(X_sigma)
                    - a_0^2 chi_mat_sigma(Hmat_sigma).

The susceptibility is the connected imaginary-time integral over 0 <= tau <= 1/T. C92 proves

    R_sigma(J)-R_good(J)
       = - integral_0^J (J-u)[D_sigma(u)-D_good(u)] du.       (1)

The requested sufficient pointwise estimate is

    D_sigma(u)-D_good(u) <= C W/t,                           (2)

with a single C independent of N, T, s and the pattern. The weaker target is

    R_sigma(J)-R_good(J)
       >= -J epsilon(J/t) W - B_N(T,s),
    epsilon(x) -> 0,  B_N(T,s)/N -> 0.                      (3)

The thermodynamic limit in (3) is taken before T tends to zero. Fixed-volume perturbation theory interchanges these limits and does not prove (3).

## 2. The response at u=0 is an exact bath-kernel contraction

[Proved] Let tau_sigma denote normalized trace on the entire finite legal internal space of sector sigma. At u=0 this space has zero Hamiltonian. The joint Gibbs state is the free gap Gibbs state tensor the normalized identity, regardless of the auxiliary energies at nonzero coupling. Define

    G_sigma,ij = tau_sigma(h_sigma,i h_sigma,j),
    K_ij(T,s) = integral_0^(1/T)
                  [<A_i(tau) A_j>_gap-a_0^2] d tau.

Direct expansion of the two susceptibilities gives

    chi_mob(X_sigma)
      = sum_ij integral <A_i(tau) A_j> G_sigma,ij d tau
                        - (a_0^2/T) tau_sigma(Hmat_sigma)^2,
    a_0^2 chi_mat(Hmat_sigma)
      = (a_0^2/T)[tau_sigma(Hmat_sigma^2)
                                      -tau_sigma(Hmat_sigma)^2].

Consequently

    D_sigma(0) = sum_ij K_ij(T,s) G_sigma,ij.                (4)

The material subtraction removes the response of a bath with constant contact a_0. It does not remove fluctuations of the mean contact between conserved vacancy sectors.

## 3. A conserved-sector crossing gives a 1/T term

[Proved] Suppose the free tilted gap Hamiltonian has two ground states, one in vacancy sector v and one in sector v+1, with equal energy, and has a positive gap to all other eigenvalues. Within either fixed-vacancy sector the hopping graph is connected. Perron-Frobenius gives a unique, translation-invariant ground vector. Denote the two contact expectations by a_v and a_(v+1).

The A_i conserve vacancy number, so their matrix elements between the two ground states vanish. Their restrictions to the two ground lines are respectively a_v and a_(v+1), independently of i. The zero-energy part of the spectral representation of K is therefore

    K_ij(T,s_*) = (a_v-a_(v+1))^2/(4T) + O_N(1).           (5)

The two ground lines have equal multiplicity. Contributions from nonzero energy denominators remain bounded as T tends to zero. Thermally occupied excited lines contribute exponentially small corrections to (5), including the terms divided by T. At fixed N the positive tilt makes the infinite vacancy sum convergent.

Combining (4) and (5) yields

    D_sigma(0)-D_good(0)
      = (a_v-a_(v+1))^2/(4T)
          [tau_sigma(Hmat_sigma^2)-tau_good(Hmat_good^2)]
        + O_N(1).                                         (6)

A positive bracket rules out (2). This is a statement about second moments, not just the square of the mean material energy.

[Proved] Excluding the single endpoint u=0 does not fix the uniform estimate. At each positive T the response is continuous in u. For every large value reached by (6), an interval of positive u has a response of the same sign and comparable magnitude. Equivalently, one can take u sufficiently smaller than T as T tends to zero. Such u remain in (0,J) for every fixed J>0.

## 4. A full C88 counterexample to the pointwise bound

[Proved, imported finite spectrum] C89 gives the free phase-model ground energy and contact at fixed N,v, with t=1 and L=N+v,

    e_N(v) = -2 sin(pi N/L)/sin(pi/L),
    a_N(v) = N/L - sin^2(pi N/L)/(N L sin^2(pi/L)).         (7)

Take N=15 and v=15. Set

    s_* = e_15(15)-e_15(16).

The rational interval calculation in crossing_certificate.py establishes the following facts using (7). The tilt lies in (1/2,3/4). Precisely the vacancy sectors 15 and 16 minimize e_15(w)+s_* w. Their contact expectations differ. All other vacancy sectors are separated in energy.

The calculation encloses pi using Machin's identity and alternating rational arctangent sums. It evaluates sine by rational interval Taylor bounds. It checks w=0,...,39 individually. For w>=40, the operator bound H_gap>=-2N and s_*>0 give a strictly larger energy than the common ground energy. Thus the check does not truncate an uncontrolled vacancy tail. The certificate output reports decimal summaries of exact rational inequalities. It gives s_*=0.610116148316937...t, a_15-a_16=0.0223207131683484..., a gap of more than 0.04743t to the other vacancy-sector ground energies, and a lower bound greater than 534.73 for the coefficient T(D_sigma-D_good) derived below.

[Proved, full internal spaces] Compare the pin-free sectors

    good = (15),       sigma = (3,3,3,3,3),       W=15.

All legal auxiliary states are retained. C88 PROOF section 1 gives the dimension of a pin-free cell with p data sites,

    d_p = (5p+3) 2^(p-1).

The X and M states carrying the additional diagonal energy D occupy

    (2p+2) 2^(p-1)

of these states. This counts the two X rails, with p positions each, and the two M endpoints. The other diagonal terms are nonnegative. Every hopping and port term has zero trace in the word basis.

For the five p=2 cells, using D=1000,

    tau_sigma(Hmat_sigma) >= 5 (6/13) 1000 = 30000/13,
    tau_sigma(Hmat_sigma^2) >= (30000/13)^2.                (8)

The second inequality is the Hilbert-Schmidt Cauchy-Schwarz inequality. It holds with the nonzero scheduled port, whose trace vanishes.

There is one marker in the good sector. The port annihilates this entire sector, as proved in C89 section 8. Its diagonal energy is bounded above by

    mu_sel + Delta + D + 14 f < 1007.

The remaining hopping and conversion row sum is less than 4, again the C89 one-marker bound. Gershgorin's bound or the row-sum operator norm therefore gives the conservative estimate

    ||Hmat_good|| < 1016,
    tau_good(Hmat_good^2) < 1016^2.                        (9)

The positive field at the schedule obeys 14 f<1. No floating-point representation of that field or of the port is needed for (8)-(9).

It follows that the bracket in (6) is strictly greater than

    (30000/13)^2 - 1016^2 > 4.29 * 10^6.                   (10)

Equations (6)-(10) prove that D_sigma(0)-D_good(0) tends to positive infinity as T tends to zero, at a fixed allowed N, W and s. This disproves (2) for the full legal C88 trace. The small-J condition in C92's contact lemma does not help because u=0 belongs to every such coupling interval.

[Scope] This obstruction does not use the vacuum-clock truncation, a scalar replacement, a finite vacancy cutoff or finite grammar penalty. It addresses the stated finite-N uniformity. It does not prove divergence of a thermodynamic response density with the thermodynamic limit taken first.

## 5. Why the integrated problem survives

[Proved] A finite-volume conserved-sector crossing also warns against deleting the global boundary term from (3). In the limit T->0 at the free crossing, a_0 tends to

    a_bar = (a_v+a_(v+1))/2.

Let E_sigma denote the frozen material ground energy. At fixed N, first-order degenerate perturbation theory within the two free ground lines gives

    Fmob_sigma(J,T=0,s_*)
      = e_* + J min(a_v E_sigma,a_(v+1) E_sigma)+O_N(J^2),
    R_sigma(J,0,s_*)
      = -J |a_v-a_(v+1)| |E_sigma|/2 + O_N(J^2).           (11)

Here R at T=0 uses the limiting free Gibbs contact a_bar. The finite internal Hilbert space and the isolated two-line bath eigenspace suffice for this expansion. An internal degeneracy is handled by diagonalizing Hmat first.

For the sectors in section 4, the imported isolated-cell certificate gives E_good=E_15<0, with |E_good|<0.073. The wrong sector obeys

    E_sigma >= 5 E_3 - 5 sqrt(2)|lambda| > 2.37.

Deleting the five ports costs at most their total norm. Thus (11) gives a negative residual difference of order J at this finite-volume tilt crossing. A uniform bound with epsilon->0 and B_N=0 is false here as well.

The allowance B_N in (3) changes this conclusion. A finite-volume density-rounding contribution can be subextensive, and a counterexample at N=15 does not contradict a bound allowing such a contribution. Nothing here shows that the crossing term remains extensive at fixed interior density as N tends to infinity. C92 explicitly allows a global subextensive error and takes that limit first.

[Assumed, corrected analytic target] The remaining useful estimate is a bound on the weighted integral in (1), with the finite-volume density-rounding term kept, or a response estimate established after taking the thermodynamic limit. A pointwise bound uniform in T before that limit is unnecessary and, as stated in (2), impossible.

## 6. Quantum second order at fixed vacancy number

[Proved, imported and restated] At fixed N,v the free bath has a unique ground state |0>. Set Q=1-|0><0| and write

    chi_ij = <0|(A_i-a) Q(H_gap-e_0)^(-1)Q(A_j-a)|0>.

This is the time integral from zero to infinity of the ground-state connected contact function. It has units 1/t. It is a positive semidefinite matrix and is circulant. Degenerate perturbation theory gives the operator on the internal space

    Heff = e_0 + J a Hmat - J^2 sum_ij chi_ij h_i h_j
                                      + O_(N,v)(J^3).     (12)

For a nondegenerate frozen ground vector phi,

    B_sigma = sum_ij chi_ij <h_i phi,h_j phi>,
    R_sigma(J) = -J^2 B_sigma+O_(N,v,sigma)(J^3).           (13)

The Gram matrix in (13) retains every quantum clock transition. Replacing it by <h_i><h_j> drops the connected internal part and changes the answer. The finite-volume zero-temperature Duhamel curvature equals 2 B in this ground-state limit. The factor of two distinguishes it from the one-sided ground-state time integral defining chi.

The error in (13) is not uniform in N or the cell lengths. Bath gaps close with N, and long-cell clock gaps close with cell length. An ordinary resolvent expansion can introduce inverse powers of both. The following result improves the second-order coefficient without asserting an improvement of that remainder.

[Proved, nonperturbative finite-volume bound] There is also an all-order ground-energy estimate that does not use a material excitation gap. Let Delta be the free bath gap at fixed N,v, let ||X_sigma||<=NM, and assume 2JNM<Delta. With P=|0><0| tensor 1, the lower-right block of H-e in the Q space is bounded below by Delta-2JNM at a ground eigenvalue e. This follows from H_QQ>=e_0+Delta-JNM and the product trial e<=e_0+JNM. The Schur complement then gives

    -J^2 ||P X_sigma Q||^2/(Delta-2JNM)
       <= R_sigma(J) <= 0.                                (13a)

The upper bound is the frozen ground state tensor the free bath trial. The lower bound follows by bounding the positive Schur-complement correction in operator norm, while retaining the entire P-space material operator. If a ground vector has no P component, its energy lies above the product trial under the stated gap condition, so it cannot invalidate the argument.

In particular, for 4JNM<=Delta,

    R_sigma-R_good >= -2J^2 N^2 M^2/Delta.                 (13b)

This holds for arbitrary bounded quantum sources and arbitrary internal degeneracies. For W>=1 it implies a finite-volume estimate of the form (3), with B_N=0 and epsilon_N=2J N^2 M^2/Delta. It proves epsilon_N->0 at each fixed N,v, including the full legal internal spaces. Its constant grows with N and it does not imply a thermodynamic statement. It also does not apply to the tilted bath at a vacancy crossing, where Delta=0 for a one-line P.

## 7. A length-weighted quantum second-order comparison

[Proved, conditional theorem] Consider the zero-port vacuum-clock model with cell lengths at least two, pins permitted, b-compatible N and a product of isolated clock ground states. Each distributed bond source has norm at most M independent of length. For the vacuum model one may take M=1+mu_sel. Suppose the free canonical kernel has a decomposition

    chi_N(r) = k_N(r) - c_N/N,
    K_N = sum_(r mod N) [1+dist_N(r,0)] |k_N(r)|.          (14)

No boundedness of K_N or c_N is assumed for the finite-N statement. Then

    |B_sigma-B_good| <= 4 M^2 b W (K_N+|c_N|).             (15)

The right side is uniform over the locations and lengths of wrong cells at this fixed bath kernel. If sup_N t(K_N+|c_N|) is finite in an interior-density window, (15) is a uniform O(W/t) bound on the full quantum second-order coefficient. It requires no lower bound on a wrong cell's clock gap.

Proof. Call the union of wrong cells and pinned b-cells the bad set. Its material length is

    L_bad = sum_wrong ell + b k_pin <= b W.

The remaining maximal runs consist of complete unpinned b-cells. Their number k is at most W, unless W=0, in which case the comparison is an equality.

For the product frozen ground vector write G_ij=<h_i phi,h_j phi>. Cauchy-Schwarz gives |G_ij|<=M^2. If i and i+r lie in the same good run and the chosen shortest arc between them stays in that run, their Gram entry agrees with the corresponding translated selected Gram entry. The latter is periodic in i with period b. Summing its value over a complete good run therefore gives the fraction length(run)/N of its full-ring sum, for each fixed r.

There are at most 2 k dist_N(r,0) good origins whose arc leaves the run. For these pairs use the bound 2 M^2 on the difference. The bad origins and the omitted fraction of the selected ring together cost at most 2 M^2 L_bad. Consequently, at each displacement,

    |sum_i G_sigma,i,i+r - sum_i G_good,i,i+r|
       <= 2 M^2 L_bad + 4 M^2 k dist_N(r,0).              (16)

This argument matches whole good runs by translation. It does not assume that a length defect leaves the marker phase on its far side unchanged. It includes correlations between different cells through their mean sources.

Multiplying (16) by |k_N(r)| and summing bounds the regular-kernel contribution by 4 M^2 b W K_N.

For the constant kernel, phi is an eigenvector of Hmat, so

    sum_ij G_sigma,ij = E_sigma^2,
    B_sigma,constant = -c_N E_sigma^2/N.

Let e_b=E_b/b. Additivity of the zero-port frozen energy gives

    |E_sigma-N e_b| <= 2 M L_bad,
    |E_sigma|, |N e_b| <= M N.

It follows that

    |E_sigma^2-E_good^2|/N <= 4 M^2 L_bad.

This bounds the constant-kernel difference by 4 |c_N| M^2 b W, completing (15).

[Scope] The proof uses a product of isolated ground states and no port. It is not a theorem about the full thermal trace, a coupled payload ground state, or a nonzero-port remainder. The vacuum marker halves are included in h_i. In this block both halves belong to the two final bond indices of their cell, so each good run in (16) consists of complete copies of the same source pattern.

[Assumed] The uniform free-kernel hypothesis sup_N t(K_N+|c_N|)<infinity is not established here. It is stronger than checking a finite set of static susceptibility matrices. C89's continuum cancellation motivates subtracting the canonical constant kernel before taking a spatial first moment. Using the raw first moment can falsely manufacture growth from that constant term.

[Heuristic] The smooth contact field has a static local response plus the canonical constant term; its noncontact time integral cancels. Descendants and higher harmonics then suggest a summable first spatial moment for k_N at interior density. This is the power counting recorded in C89 and matching_attack. It is not an exact, all-distance phase-model correlation bound. The exact finite-cell gradient identity does not prove a uniform all-order expansion either, since its partial-sum operators can grow with cell length.

[Assumed, still missing even if (14) is bounded] To obtain (3), one still needs control of the common-bulk-subtracted terms beyond second order, uniformly in W, all wrong lengths, temperature, and volume, with the density-rounding contribution retained. A bound on the free two-point function alone does not control the coupled multi-time quantum sources.

## 8. Exact diagonalization tests

[Numerical] checks.py uses the exact phase-model composition basis and the C89 clock construction, copied locally as inherited_ed.py. The quantum-source resolvent calculation is copied locally as inherited_second_order.py. It checks the exact free energy and mean contact formulas before each bath calculation. The new driver computes 50 canonical comparisons at J/t=0.04, 0.02, 0.01, 0.005 and 0.0025. No scalar surrogate is included.

The systems are (N,v)=(4,4),(6,6),(9,4),(30,1), with selected b=2,3,3,15 respectively. They include pins, local compensated length defects, phase shifts across wrong cells, and a whole long-cell replacement. The largest matrix has dimension 13365. The largest eigenpair residual is 1.29e-11. The largest bath resolvent residual is 4.28e-13. All 50 total defect costs are positive.

Define the observed one-sided coefficient and its two-sided diagnostic by

    epsilon_observed = max(0,-Delta R)/(J W),
    epsilon_absolute = |Delta R|/(J W).

These are finite-sample quantities, not a supremum over all patterns or sizes. In particular, a zero observed one-sided value means that the residual has favorable sign, not that the response vanishes.

For the unfavorable length-30 replacement of (15,15), at N=30,v=1, epsilon_observed decreases from 3.15549e-5 at J/t=0.04 to 8.21655e-7 at J/t=0.0025. The second-order limit is

    (B_sigma-B_good)/W = 0.000301406301669...

and epsilon_observed/(J/t) approaches it. This is a high-occupancy b=15 test, not a half-density b=15 calculation.

At half density, the b=3 wrong pattern (2,4) has favorable Delta R. Its epsilon_absolute decreases from 1.22549e-4 to 7.68192e-6 over the same J range. A single pin in (3,3) gives 2.79568e-3 and 1.75224e-4. These values resolve the linear decrease of the normalized residual rather than relying on its favorable sign alone.

[Numerical, thermal traces] Ten full logical-clock thermal comparisons use selected (4), wrong (2,2), N=4, s=0.65t and total vacancy cutoff 9. Each conserved-vacancy block is diagonalized completely, with maximum dimension 880. The tests take T/J=0.1 and 1 at the five couplings above. They subtract the complete frozen material free energy at the same free contact a_0(T,s), including its entropy.

At T=J the residual is unfavorable. Its epsilon_observed decreases from 4.68637e-4 at J/t=0.04 to 9.47357e-6 at J/t=0.0025. The last three values of Delta R/(J^2 W) are approximately -0.00383056, -0.00378039 and -0.00378943. At T=0.1J the sign is favorable and the absolute normalized residual decreases from 6.67822e-4 to 4.15113e-5.

[Proved cutoff estimate, numerically evaluated] For each omitted vacancy block,

    dim_gap(v) = binom(v+N-1,N-1),
    Hmob(v) >= s v - 2tN - J sum_i ||h_i||.

The ratio of successive upper-bound terms beyond v=v_0 is at most

    q = [(v_0+N)/(v_0+1)] exp(-s/T) < 1.

The first omitted term divided by 1-q bounds the missing partition sum. The same bound with J=0 controls the free contact error, since 0<=A_i<=1. If r_sigma,r_good,r_free are the resulting relative tail bounds, the residual-difference error is at most

    T(r_sigma+r_good)
      + J(||Hmat_sigma||+||Hmat_good||) r_free.             (17)

The output retains these bounds. The largest evaluated residual-difference cutoff error bound is 9.11e-10, on the J=T=0.04 row, compared with a residual magnitude 7.50e-5. Zero decimal tail outputs mean numerical underflow of a positive exponentially small bound, not an assertion that the tail is absent.

[Numerical, response obstruction] At u=0 the separate b=4 logical trace test uses s_*=e_4(4)-e_4(5)=0.532518623638...t. It evaluates (4) from complete spectral decompositions of the bath blocks with v<=7. This cutoff contains both true ground sectors and hence reproduces the exact 1/T coefficient. The finite-temperature values are cutoff diagnostics. As T/t decreases from 0.01 to 0.0001, D_sigma-D_good grows from 0.511428 to 44.640966. The product T(D_sigma-D_good) approaches

    (a_4-a_5)^2 [3 mu_sel^2+1/2]/4
      = 0.00445752917643...

The observed value at the smallest temperature is 0.00446409655920. The factor 3 mu_sel^2+1/2 follows from normalized traces of the two clock Hamiltonians. This is a separate small-matrix confirmation of the crossing mechanism, not the proof of the full C88 obstruction in section 4.

[Numerical, size dependence at second order] coefficient_scaling.py computes the half-density kernel at N=v=3,6,9, with maximum bath dimension 24310. For one b=3 pin, Delta B/W is -0.0932690, -0.0701029 and -0.0616633. For a compensated (2,4) replacement at N=6,9 it is -0.00307337 and -0.00304739. The finite-size diagnostic choice c_N=-N chi_N(floor(N/2)) gives regular first moments 0.0931373, 0.106396 and 0.107905. These values do not establish the uniform hypothesis in section 7.

## 9. Consistency checks and remaining work

[Proved] At full occupancy every A_i=1, so chi=0 and R=0. The response crossing mechanism requires two distinct vacancy sectors and is absent in a fixed-vacancy ensemble. At fixed vacancy number the second-order B is nonnegative, while its difference between patterns has no fixed sign. The numerical long-cell row has the unfavorable sign, and the pin rows have the favorable sign. Both occur in the formulas above.

The rank-one term in (14) scales as E_sigma^2/N. Its difference is O(W), not O(N), after the common selected bulk is matched. All coefficients in (12)-(15) have the required units, with chi, K_N and c_N proportional to 1/t. An extensive long-cell replacement carries W equal to its material length, in agreement with matching_attack. No charge seam is cut at every defect.

[Status] The pointwise estimate (2) is disproved with its stated uniformity. A length-weighted quantum second-order comparison is proved under an explicit free-kernel hypothesis, without a wrong-cell gap assumption. Finite canonical and thermal ED show the normalized integrated residual decreasing toward zero on the tested models. A uniform thermodynamic remainder estimate, the full auxiliary thermal comparison away from the endpoint obstruction, and the finite-K relative correction remain open. No mobile sector-selection or central-charge theorem follows from these finite tests.

[Provenance] The project STATE.md, C92, C89, the C88 referee README and matching_attack README were read. C88's literal grammar was additionally consulted for the full legal trace count in (8). All scripts, outputs and notes are in work_C94. No network commands, delegated agents, external edits or cluster jobs were used. Project STATE.md was left unchanged under the explicit folder restriction.
