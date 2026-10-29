# C106 result. Counting fast modes

2026-10-29.

[Proved] The proposed uniform thermal crossover lemma is false, even for a nearest-neighbor free chain with d = 4. The ratio v/M can stay fixed and positive while the temperature interval of fixed-fraction CFT entropy matching collapses to zero.

[Open] This does not disprove the velocity-weighted bound or an integrable bound on N(u). No general interacting counting bound is proved here. The operator route still needs two explicit additional lemmas, stated below.

## 1. A microscopic obstruction to the thermal route

[Proved] Put two complex fermion modes on each site and take

    H_epsilon = sum_x [f_x^dagger (-i t sigma_z) f_(x+1) + h.c.
                       + f_x^dagger (epsilon sigma_x-a epsilon 1) f_x],
    t > 0,       a > 1.

The Jordan-Wigner representation is a one-site translation-invariant nearest-neighbor spin chain with d = 4. The bands are

    E_+(k) = -a epsilon + sqrt(4t^2 sin^2 k+epsilon^2),
    E_-(k) = -a epsilon - sqrt(4t^2 sin^2 k+epsilon^2).

For small epsilon there are four simple Fermi points in E_+, giving two c = 1 sectors with the common speed

    v = 2t sqrt(1-a^(-2)) sqrt[1-epsilon^2(a^2-1)/(4t^2)].

The gauge-optimized bond scale is exactly

    M = sqrt(4t^2+epsilon^2).

The proof of the gauge optimization uses the two extremal bond eigenstates. Each has equal left and right one-site marginals, so their energy difference cannot be changed by a bond coboundary.

[Proved] For any fixed r in (0,1), a smooth choice a(epsilon) keeps v/M = r exactly. The entropy has two different limiting slopes:

    s(T)/T -> 2pi/(3v)       as T -> 0 at fixed epsilon,
    s(T)/T -> pi/(3t)        when epsilon/T -> 0 and T/t -> 0.

Consequently,

    s(T) / [2pi T/(3v)] -> r

in the second regime. Fix a matching fraction alpha in (0,1), choose r < alpha, and define T_alpha as the largest endpoint for which

    s(T) >= alpha 2pi T/(3v)

holds at every smaller positive temperature. Then

    T_alpha = O(epsilon),       T_alpha/M -> 0,

although d = 4 and v/M = r are fixed. Rescaling each Hamiltonian also fixes M = 1 and v = r. Thus no positive lower bound T_alpha >= M F_alpha(d,v/M) holds universally. This includes T_alpha >= kappa v^2/M and every fixed-power variant.

[Proved mechanism] The Bloch matrix has uniformly bounded first and second derivatives, but the upper-band curvature at a Fermi point diverges as 4t^2/(a^3 epsilon). The avoided-crossing gap appears in the denominator when differentiating individual band eigenvalues. A bound on the local interaction norm or on propagation speed does not control this gap.

[Proved scope] This is not a counterexample to A <= C M f(d). Here

    N(u) = 2 for u <= v,       A/M = 2r < 2,

so the already proved free-chain bound is respected. The example rules out a uniform window with a fixed matching fraction. It does not rule out every weaker entropy inequality, for example one with a sufficiently small velocity-dependent fraction.

## 2. The counting reductions that do hold

[Proved] With V = C_LR M and the inherited local real-time identification, C_LR = 6 e is sufficient. Then

    A = integral_0^V N(u) du.

An estimate N(u) <= g(d)(M/u)^(1-delta), 0 < delta <= 1, gives

    A <= [C_LR^delta/delta] M g(d).

A uniform N(u) <= g(d) needs no logarithm, but also bounds the total central charge. The borderline N(u) <= g(d) M/u leaves log(V/v_min), which is not uniform. A logarithmic tail in V/u is integrable.

[Assumed thermal window lemma; proved implication] Suppose that for each velocity window w <= v_a < 2w, at T_w = kappa M(w/M)^p,

    s(T_w) >= (alpha pi T_w/3) sum_window c_a/v_a.

Entropy capacity gives

    sum_window c_a <= K (M/w)^(p-1),
    K = 6 log d/(alpha pi kappa).

For 1 < p < 2 this gives

    N(u) <= K (M/u)^(p-1)/(1-2^(1-p)).

At p = 1 it gives N(u) <= K[1+log_2(V/u)]. A direct dyadic sum gives, for 1 <= p < 2,

    A <= [2^(p-1) K C_LR^(2-p)/(1-2^(p-2))] M.

[Assumed entropy-increment lemma; proved implication] The quadratic case can also suffice if the matching controls increments. Let w_j = V/2^(j+1), T_j = kappa w_j^2/M, and I_j = {a: w_j < v_a <= 2w_j}. If

    s(T_j)-s(T_j/2) >= (alpha pi T_j/6) sum_(a in I_j) c_a/v_a,

then the disjoint temperature intervals spend at most the total entropy log d, giving

    A <= 24 M log d/(alpha pi kappa).

Two entropy lower bounds cannot be subtracted to obtain this increment lemma. PROOF.md Section 4.7 also shows that the universal versions of both thermal matching assumptions above fail in the explicit family, for any fixed positive alpha and kappa. Their conditional implications remain valid for restricted classes.

## 3. The operator lemma still needed

[Proved] A k-site operator space has dimension d^(2k), but one local observable can overlap with many independent infrared fields. Nonzero overlap with each stress tensor does not establish a surjective scaling map. Also, one stress tensor counts one sector, irrespective of c. Its vacuum descendant has norm squared c/2, while its span remains one-dimensional. General unitary CFT sectors need not have spin-one currents.

[Assumed sufficient operator lemma] For each u, require an infrared field space W_u and a scaling map onto it from operators on a fixed k-site block, with

    dim W_u >= eta N(u),
    k <= K(d)+beta(d) log(V/u),
    gamma(d) = 2 beta(d) log d < 1.

[Proved conditional implication] Dimension counting and integration then give

    N(u) <= eta^(-1) d^[2K(d)] (V/u)^gamma(d),
    A <= 6 e M d^[2K(d)]/[eta(1-gamma(d))].

[Open] Neither the required number of independent fields per unit central charge nor the quantitative block-support bound follows from the stated assumptions. A support bound proportional to M/u gives only an exponential tail through this argument, which is not integrable. This identifies the unresolved operator input without treating it as proved.

## 4. Verification and handoff

[Numerical] The small bond-matrix and entropy-quadrature checks in checks.py pass. They verify the norm formula, equal extremal marginals, fixed v/M, both thermal regimes, and failures of the one-temperature and increment matching inequalities. CHECKS.md records the results. The proof itself is analytic.

[Proved work record] PROOF.md gives all derivations and qualifications. PROGRESS.md and local STATE.md give the resumable handoff. All files are in work_C106. No heavy numerical job, cluster submission, network command, commit, or push was run. The project-level STATE.md was not edited because it is outside the authorized folder.
