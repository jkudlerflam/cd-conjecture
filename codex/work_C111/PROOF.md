# C111. The proposed smooth-modulation lemma is false as stated

## 1. Statement and scope

[Proved, main result] There is a family of one-site translation-invariant nearest-neighbor bonds on local dimension d = 8, each attaining its gauge-optimized half-spread M = 3/2, for which the coefficient of the quadratic smooth-modulation deficit is unbounded. The uniform infrared theory is the ordinary XX conformal sector, with c = v = A = 1, together with a strictly gapped spectator. The conformal identification of XX is the input already recorded in the background, not a new lattice-to-CFT theorem here. All statements about the spectator and the energy inequalities below are elementary finite-dimensional results.

[Proved] More explicitly, fix lambda = 1/10 and let 0 < delta <= 1/20. For smooth positive periodic profiles w_i = W(i/L),

    liminf_(L -> infinity, L even) Def(w)/Q(w) >= lambda^2/(2 delta),
    Q(w) = sum_i (w_(i+1)-w_i)^2/w_i.

For the quadratic open-interval weights w_j = j(N-j)/N, 1 <= j <= N-1,

    Def(w) >= (lambda^2/delta) log N + O_delta(1),
    Q(w) = 2 log N + O(1).

Here an O_delta(1) in a lower bound denotes a signed remainder bounded in absolute value at fixed delta. Therefore no finite universal C f(8) bounds Def/(M Q). The limit in length is taken first for each fixed positive delta. No uniform crossover assumption is used.

[Assumed, interpretation of the proposed remainder] For bounded smooth periodic profiles, a remainder must be o(Q) if the inequality is to bound a gradient coefficient. For the quadratic modular profile it must be o(log N), as in the usual O(1) remainder. Allowing another arbitrary multiple of Q or log N in the remainder would not imply the desired bound on A. The counterexample rules out the stated lemma with these useful interpretations.

[Proved, limitation] This refutes the proposed smooth-modulation lemma, not A <= C M f(d). In the counterexample A = 1 while M = 3/2. It also refutes identifying the logarithmic weighted deficit with A/(12 pi) for an arbitrary energy-density representative, including every representative attaining M. A separately specified representative and a subtraction of local response terms would be a different problem.

## 2. The basic two-level response

[Proved] Let X and Z be Pauli matrices, with eigenvalues +/-1, and set

    b_delta = (delta/2)(X tensor I + I tensor X)
                  + lambda(Z tensor I - I tensor Z).

On a uniform ring the second term telescopes and

    H_1 = delta sum_i X_i,       e0 = -delta.

For positive periodic weights define

    u_i = (w_i+w_(i-1))/2,      v_i = w_i-w_(i-1).

The weighted Hamiltonian is a sum of independent on-site matrices,

    H_w = sum_i [delta u_i X_i + lambda v_i Z_i].

Their negative eigenvalues give the exact deficit

    D_delta(w) = sum_i [sqrt(delta^2 u_i^2 + lambda^2 v_i^2) - delta u_i]
               = sum_i lambda^2 v_i^2 /
                    [sqrt(delta^2 u_i^2 + lambda^2 v_i^2) + delta u_i].       (2.1)

The identity sum_i u_i = sum_i w_i was used. For an open chain use virtual weights w_0 = w_N = 0, with u_i = (w_i+w_(i-1))/2 and v_i = w_i-w_(i-1) for 1 <= i <= N. The same identity and formula hold.

[Proved] For W positive, nonconstant, smooth and periodic, w_i = W(i/L), the lower bound on W gives u_i bounded away from zero and v_i = O(L^-1). Rationalization or Taylor expansion yields, at fixed delta,

    D_delta(w) = (lambda^2/(2 delta)) sum_i v_i^2/u_i + O_delta(L^-3).

Also

    sum_i v_i^2/u_i = Q(w) + O(L^-2),
    Q(w) = L^-1 integral_0^1 W'(s)^2/W(s) ds + O(L^-2).

For the first estimate compare 1/u_i with 1/w_(i-1). Their difference is O(L^-1), multiplied by v_i^2 = O(L^-2) and summed over L sites. The second is a Riemann sum. Thus

    D_delta(w)/Q(w) -> lambda^2/(2 delta).                       (2.2)

[Proved, checks] Formula (2.1) vanishes at constant w on a ring and at lambda = 0. Both sides are homogeneous of degree one in w. Under simultaneous rescaling of delta and lambda by an energy factor, D scales by that factor, as does lambda^2/delta. The coefficient diverges as delta -> 0 although the gap is strictly positive at every fixed member of the family.

[Proved] The simple bond already shows why a commutator-norm argument cannot control this response. Adjacent bonds have

    [b_(12), b_(23)] = delta lambda [X_2,Z_2],
    norm([b_(12),b_(23)]) = 2 delta |lambda|.

This norm tends to zero while the quadratic response coefficient tends to infinity. The operator norm of b_delta is sqrt(delta^2+4 lambda^2), so using the raw bond norm instead of M does not repair a universal bound either. The gauge-optimized M of this simple bond is delta, attained after removing the divergence. The next construction makes the offending representative itself an optimizer with M fixed.

## 3. A bond attaining M with the same response

[Proved, definition] Take a four-dimensional site space

    V = V_A direct-sum V_B,       dim V_A = dim V_B = 2.

The label A or B at each site is conserved by the bond. Set a = 3/4 and define the block-diagonal bond s_delta by

    AA block:  -a I + b_delta,
    BB block:  P_B, the swap of the two B qubits,
    AB and BA blocks:  +I.

The parameter range remains lambda = 1/10 and 0 < delta <= 1/20.

[Proved] The AA eigenvalues are

    -a - 2r,  -a,  -a,  -a + 2r,
    r = sqrt((delta/2)^2 + lambda^2).

Indeed b_delta is the sum of (delta X/2 + lambda Z) on the first site and (delta X/2 - lambda Z) on the second. Their spectra are +/-r. In the stated range,

    a + 2r <= 3/4 + sqrt(17)/20 < 1.

All AA eigenvalues lie strictly between -1 and 1. The BB swap has eigenvalues -1 and 1, and both mixed blocks equal 1. Hence the raw half-spread of s_delta is 1.

[Proved] This is also the gauge-optimized half-spread. Let eta_- be the B singlet and eta_+ the B product vector |0_B 0_B>. Each has identical left and right reduced density matrices. For every Hermitian on-site o,

    <eta_+|s_delta + o tensor I - I tensor o|eta_+> = 1,
    <eta_-|s_delta + o tensor I - I tensor o|eta_-> = -1.

The variational characterization of the largest and smallest eigenvalues forces a spread at least 2 in every gauge. Scalar shifts leave the spread unchanged. Since s_delta has spread 2,

    M(s_delta) = 1, and s_delta attains the minimum.             (3.1)

[Proved, stronger gauge observation] For this four-state spectator, the minimizing representative is unique up to a scalar shift. To see this, omit the scalar shift and suppose s_delta + o tensor I - I tensor o still has spread 2. The two witnesses force its spectrum to lie in [-1,1]. Every normalized mixed product vector |a_A b_B> then has expectation

    1 + <a_A|o_AA|a_A> - <b_B|o_BB|b_B> <= 1.

The reversed product gives the opposite inequality. Thus all these expectations of o agree, so o_AA and o_BB are the same scalar on their respective blocks. Subtract that scalar. Each mixed vector now has expectation 1 and hence must be an eigenvector of the bounded-above operator with eigenvalue 1. The AA and BB components of its image under the gauge term are

    -|a_A> tensor o_AB|b_B>,       o_BA|a_A> tensor |b_B>.

They belong to orthogonal label sectors and must both vanish for all a_A,b_B. Consequently o_AB = o_BA = 0, and o was scalar. This uniqueness claim is for s_delta; uniqueness is not needed or asserted for the eight-state bond below.

## 4. The uniform ground state and gap

[Proved] In the all-A sector of a periodic chain,

    H_A = -a L + delta sum_i X_i,
    E_A = (-a-delta)L.

Its unique ground vector is the product of the X = -1 vectors. Its gap is 2 delta.

[Proved] No other label sector has lower energy. First consider a uniform open A segment with k sites. For k >= 2 its energy is

    E_A^open(k) = -a(k-1) - delta(k-2) - 2r
               = k(-a-delta) + a + 2 delta - 2r.

The interior sites have on-site coefficient delta X; the two end sites have delta X/2 +/- lambda Z, giving -r each. The excess a + 2 delta - 2r is positive throughout the parameter range. For k = 1 the energy is zero, also above k(-a-delta).

[Proved] For three B qubits the operator P_12 + P_23 has minimum eigenvalue -1. One derivation decomposes by the number of qubits in state 1. On the one-excitation basis |100>, |010>, |001> the matrix is

    [1 1 0]
    [1 0 1]
    [0 1 1],

whose eigenvalues are 2, 1, -1. The two-excitation block is identical; the zero- and three-excitation blocks have value 2. Pair the consecutive bonds of an open B segment into such two-bond groups. An unpaired bond is bounded below by -1. Thus, for k sites,

    E_B^open(k) >= -ceil((k-1)/2) >= -k/2.                     (4.1)

For an all-B ring with L >= 3, summing the translated three-site inequality gives 2 H_B >= -L and hence E_B >= -L/2.

[Proved] In a mixed label configuration, the Hamiltonian is the sum of independent open A and B segment Hamiltonians and +1 at every interface. If n_A and n_B are the numbers of sites of each type and m is the number of interfaces, its energy satisfies

    E >= (-a-delta)n_A - n_B/2 + m,
    E - (-a-delta)L >= (a+delta-1/2)n_B + m > 0.

The all-B ring is also above the all-A ground energy by at least (a+delta-1/2)L. Thus the uniform four-state chain has the unique all-A product ground state and, for L >= 3, gap 2 delta. Label defects have a larger positive cost. In particular the thermodynamic ground energy per bond is

    e_s = -a-delta.                                           (4.2)

[Proved] In the all-A sector of a weighted chain the exact energy is

    E_A(w) = -a sum w_i - sum_i sqrt(delta^2 u_i^2 + lambda^2 v_i^2).

The unrestricted weighted ground energy can only be lower. Therefore

    e_s sum w_i - E_0(H_w^s) >= D_delta(w).                    (4.3)

No claim that the weighted ground state remains in the all-A sector is needed. A trial state is being used in the valid direction to disprove an upper deficit bound.

## 5. Adding a conformal infrared sector at fixed dimension

[Proved, construction] At each site take C^2 tensor V, of dimension d = 8, and let

    h_delta = h_XX + s_delta,
    h_XX = (X tensor X + Y tensor Y)/4.

The two summands act on the independent XX and spectator factors of adjacent sites. Uniform and weighted Hamiltonians split into these two factors. The uniform energy density is

    e0 = -1/pi - a - delta.

[Assumed, supplied infrared input] The XX factor has c = 1 and v = 1 in this normalization, as recorded in the background and C105. Section 4 proves that the spectator has a strictly positive gap. Thus at each fixed delta > 0 the infrared has the same single conformal sector and A = 1. There are no additional gapless label modes or zero-energy local spectator modes.

[Proved] The raw extremal eigenvalues of h_delta are +/-3/2, since the eigenvalues of its independent bond factors add. To show gauge optimality, take the tensor product of the XX antisymmetric one-particle bond vector with the B singlet for the lower witness, and the XX symmetric one-particle bond vector with |0_B 0_B> for the upper witness. These vectors have energies -3/2 and 3/2. Each has identical left and right one-site marginals on the full eight-dimensional site space. Every on-site gauge term therefore has zero expectation in each witness. As in Section 3,

    M(h_delta) = 3/2, and h_delta itself attains it.             (5.1)

[Proved] Deficits add for the two independent factors. On an open interval the XX deficit relative to -1/pi is nonnegative, by using the restriction of its infinite-chain ground state as a trial density matrix. For even periodic lengths it is also nonnegative. In that case the uniform XX energy is

    E_XX(L) = -csc(pi/L) <= -L/pi.

For completeness, the Jordan-Wigner fermions at half filling have occupied consecutive momenta separated by 2 pi/L, centered on pi. The boundary condition is periodic or antiperiodic according to the fermion parity, placing those L/2 momenta symmetrically around pi. The sum of their cosines is

    -sin((L/2)pi/L)/sin(pi/L) = -csc(pi/L).

A translation average of a uniform finite-ring ground state has each bond expectation E_XX(L)/L. Its weighted expectation bounds E_0(H_w^XX) above by E_XX(L) sum w_i/L. The inequality sin t <= t gives the asserted sign of the deficit.

[Proved] Combining this with (4.3), on open intervals and on even rings,

    Def_h(w) = Def_XX(w) + Def_s(w) >= D_delta(w).              (5.2)

Equations (2.2) and (5.1) now imply

    liminf Def_h(w)/(M Q(w)) >= 1/(300 delta).

Given any proposed finite K = C f(8), choose a fixed positive delta <= 1/20 with 1/(300 delta) > K, then take the smooth limit. This contradicts the proposed inequality with an o(Q) remainder. The uniform infrared sector and the optimized bond scale stay fixed throughout the family.

## 6. The modular profile gives the same obstruction at logarithmic order

[Proved] Let w_j = j(N-j)/N for 1 <= j <= N-1, and set w_0 = w_N = 0. At site i put y_i = i-1/2. Then

    v_i = 1 - 2 y_i/N,
    u_i = y_i(1-y_i/N) - 1/(4N).

Write u_i^0 = y_i(1-y_i/N). A direct partial-fraction identity gives

    v_i^2/u_i^0 = 1/y_i + 1/(N-y_i) - 4/N.

The replacement of u_i^0 by u_i changes the summed expression by O(1/N). Indeed u_i is uniformly comparable to u_i^0, and on the left half u_i^0 is comparable to y_i, so the absolute error is at most a constant times (N y_i^2)^-1. Sum this convergent square tail and reflect the right half. Hence

    sum_(i=1)^N v_i^2/u_i
       = 2 sum_(i=1)^N 1/(i-1/2) - 4 + O(1/N)
       = 2 log N + O(1).                                      (6.1)

[Proved] For t >= 0, the elementary bounds

    1 + t/2 - t^2/8 <= sqrt(1+t) <= 1+t/2

give

    0 <= lambda^2 v_i^2/(2 delta u_i)
              - [sqrt(delta^2 u_i^2+lambda^2 v_i^2)-delta u_i]
       <= lambda^4 v_i^4/(8 delta^3 u_i^3).

Here |v_i| <= 1 and u_i is bounded below by a fixed positive multiple of min(i-1/2,N-i+1/2). Thus sum_i v_i^4/u_i^3 is uniformly bounded in N. Combining with (6.1), at fixed delta,

    D_delta(w) = (lambda^2/delta) log N + O_delta(1).            (6.2)

The constant in the remainder need not remain bounded as delta -> 0. No such uniformity is required for the contradiction.

[Proved] For the Q in the proposed lemma, summed over adjacent positive bonds,

    Q_N = sum_(j=1)^(N-2) (w_(j+1)-w_j)^2/w_j
        = 2 log N + O(1).                                    (6.3)

To check this directly, w_(j+1)-w_j = 1-(2j+1)/N. Replacing its square by (1-2j/N)^2 changes Q_N by O(log N/N), since the numerator changes by O(1/N) and sum_j 1/w_j = O(log N). The resulting summand is

    1/j + 1/(N-j) - 4/N.

Summing proves (6.3). Adding or omitting a bounded number of endpoint terms with nonzero denominators changes only O(1).

[Proved] Equations (5.2), (6.2) and (6.3) imply

    liminf_(N -> infinity) Def_h(w)/(M Q_N) >= 1/(300 delta).

Thus the obstruction persists in the specific modular geometry used to extract A. It is not a finite-ring Casimir offset, an unsmoothed staggered profile, or a boundary term of bounded size. It contributes to the same logarithm as the proposed conformal anomaly.

## 7. What happens to the three proposed proof routes

[Proved] Route (i), a local dilation of the uniform ground state, gives E_0(H_w) <= <H_w>_trial and hence a lower deficit bound. To prove the requested upper bound it would require an inverse comparison valid for every state. The counterexample shows why a universal quadratic cost cannot hold for such a comparison. On the all-A subspace, rotating the local spin into the direction delta u_i X + lambda v_i Z lowers its energy by (2.1). The rotation angle satisfies tan theta_i = lambda v_i/(delta u_i). The inverse gap is unavoidable even though the Hamiltonian is local and its bond norms are bounded. Section 2 also gives an adjacent-commutator norm tending to zero while this response diverges. Lieb-Robinson propagation estimates do not remove this explicit static susceptibility.

[Proved, using C105] Route (ii), positive block-Hamiltonian bounds on variable windows, is already covered by C105's obstruction. If w_i = sum_I c_I 1_(i in I), c_I >= 0, with intervals of any lengths and overlaps, then sum_I c_I >= max_i w_i. For XX, C105 proved

    |I| e0 - E_0(H_I^open) >= g_* = 1/sqrt(2)-2/pi > 0.

Therefore the sum of the scalar block allowances is at least g_* max w. This remains true if every interval length is chosen from a local gradient scale. The premise that C105 excluded only equal-length windows is narrower than its actual theorem. Signed corrections or additional boundary operators are outside that obstruction, but they cannot establish the universal inequality disproved here.

[Proved] Route (iii), curvature in the weights, exposes the same failure directly. For the on-site energy

    e(u,v) = -sqrt(delta^2 u^2+lambda^2 v^2),

one has

    partial_v^2 e(u,0) = -lambda^2/(delta u).

This is unbounded as delta -> 0 at fixed u. Positivity, concavity and bounded local dimension do not supply a uniform quadratic curvature constant. The d = 8 construction embeds this exact curvature in a model with M fixed and a nonempty conformal infrared theory. A Lieb-Thirring-type bound that includes inverse-gap or susceptibility data may still be possible, but it is not the proposed crossover-free lemma.

## 8. Consequences for the infrared matching

[Proved] A discrete divergence produces a field proportional to the first difference of w. Its linear expectation can vanish and its quadratic energy response can nevertheless be of order Q(w). In this construction the response comes entirely from a gapped degree of freedom. Integrating it out gives the local term -lambda^2 (w')^2/(2 delta w), plus higher gradients. Equation (2.1) proves this statement without an effective-field-theory assumption.

[Proved] Consequently, the assertion in the modular README that total derivatives are subleading is not valid for every bond representative in the requested class. They can be subleading in the operator expansion and still contribute at the same gradient order as the ground-energy anomaly. Norm minimization of the bond does not exclude them, as (3.1) and (5.1) show. The special free and reflection-symmetric examples checked previously do not establish the general matching.

[Assumed, open corrected target] A possible repair is to define an energy-density representative by an additional physical prescription and separate the nonuniversal local gradient response from the conformal stress contribution. No theorem establishing such a prescription, bounding the subtracted coefficient, or treating every infrared sector follows from this work. The velocity-weighted conjecture remains open beyond the classes already settled in the background.

## 9. Numerical checks and provenance

[Numerical] checks.py verifies the 16 by 16 spectator spectrum, the 64 by 64 stacked bond spectrum, and the three-B-site operator inequality. Full, unrestricted spectator diagonalizations at L = 3 and 4 reproduce E0 = (-a-delta)L and gap 2 delta for delta = .05, .01, .001. Smooth scalar evaluations agree with (2.2) and (6.2). The scalar sums use at most 65536 sites and require no many-body diagonalization.

[Numerical] The same script checks smooth finite-ring deficits on XX, Heisenberg, XXZ(1/2), SU(3), and four seeded random Hermitian two-qubit bonds. The maximum Hilbert-space dimension is 729. The reference is E0(H_1) sum w/L, removing the finite-ring uniform energy offset. For random bonds the reported scale is the raw half-spread, not a claimed gauge optimum. These small checks are diagnostics, not infrared limits or certified eigenvalue bounds.

[Proved, infrastructure record] The Typhon connection command returned "Connection timed out during banner exchange". No job was submitted. large_gaussian.py and cluster_job.sbatch prepare the requested large Gaussian XX calculation, including the stacked spectator lower bound. The Python script rejects execution without a SLURM compute allocation and rejects login hostnames. No large calculation ran locally or on a login node.

[Proved, sources and file boundary] Inputs are the task specifications, REPORT.md Sections 9.68 through 9.77, contrib/fable_bound_modular/README.md, and work_C105. No new literature theorem is required for the counterexample. No files outside work_C111 were modified. The project-level memory was not edited because the task explicitly restricts writes to this folder; STATE.md here records the result for integration.
