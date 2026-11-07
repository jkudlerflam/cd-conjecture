# C132 derivations

All statements carry a status tag. The bond convention is H_L = sum_x h_(x,x+1), including both oriented bonds when L = 2. The gauge scale is M. No general finite K(2) is claimed. When M = 0, h is a scalar plus a telescoping gauge term, so all energy defects here vanish. The structural results below hold for every finite d, including d = 2.

## 1. Marginal reformulation [Proved]

Let G consist of scalar operators and o tensor 1 minus 1 tensor o. Centering the spectrum shows

    M(h) = inf_(a,o) ||h + a 1 + o tensor 1 - 1 tensor o||.

This is the quotient norm on Hermitian bonds modulo G. All two-site marginals in Omega_L and Omega_infinity have unit trace and equal left and right one-site marginals. Their differences annihilate G. The dual norm of such a difference D is

    ||D||_* = sup_(M(h)<=1) Tr(h D).

In fact ||D||_* = ||D||_1, not merely the inequality in the source. One inequality follows by centering and optimizing the gauge. For the reverse inequality, every ||h|| <= 1 has M(h) <= 1, and ordinary trace-norm duality applies. Thus on this annihilator the purported dual-M distance is ordinary trace distance, without the conventional factor 1/2.

Translation averaging gives E_L/L = min_(rho in Omega_L) Tr(h rho). Compactness of the state space and the variational principle give e_0 = min_(sigma in Omega_infinity) Tr(h sigma). Taking support functions in direction -h proves, in both directions,

    E_L >= L e_0 - K M/L for all h
    iff Omega_L is contained in Omega_infinity + (K/L^2) B_*.

The sets and the ball live in the common affine space. Separation takes place in its finite-dimensional translation space. This verifies claim 1, and strengthens its norm identification.

A qualification matters. The equivalence holds at a fixed length when the energy inequality is required for every h. An h-dependent restriction to vacuum-compatible lengths does not directly give an inclusion for all of Omega_L. In this note the proposed sufficient lemma is stated for all lengths, which is stronger and avoids that quantifier problem. Its CFT consequence is taken only along a sequence with the vacuum Casimir expansion.

## 2. Arbitrary CPTP recovery periodizes [Proved]

Fix 2 <= l <= L-1. Write rho_l = rho_AB with A one site and B comprising l-1 sites. Identify the shifted l-site target with rho_BC = rho_l, and rho_B = rho_(l-1). Let R be any CPTP map with R(rho_B) = rho_BC.

Start with rho_l on sites 1 through l. At step n append a site by applying R to sites n-l+2 through n and keeping all earlier sites. The last l sites remain rho_l. Indeed, tracing the first site of the last l+1 sites gives R(rho_(l-1)) = rho_l. This proves the invariant by induction.

The bond (k,k+1) is last touched at n = k+l-1. Immediately before that step its enclosing l-site window has state rho_l, so immediately afterward its marginal is

    sigma_2(R,l) = Tr_(B_2 ... B_(l-1) C) [(id_A tensor R)(rho_AB)].

Later steps act to its right and cannot change it. More generally every fixed initial finite window stabilizes once the moving channel has passed it. This defines a consistent half-line state. Embed successively longer interior windows into Z and average translates. Compactness gives a translation-invariant limit, and every bond in that limit is sigma_2(R,l). Therefore

    e_0(h) <= Tr[h sigma_2(R,l)].

Neither Petz structure nor h-independence was used. A Petz map initially defined on supp(rho_B) is completed by

    R_full(X) = R(P X P) + Tr[(1-P)X] tau,

where tau is any output state. Inputs encountered in the construction are supported on P, so this completion changes no result.

Both one-site marginals of sigma_2 are rho_1. The A marginal is preserved by trace preservation, and the B_1 marginal follows from R(rho_B) = rho_l. Consequently, for X = (id tensor R)(rho_AB),

    delta(R,l) = Tr[(h_tilde_(A B_1) tensor 1)(X-rho_ABC)]
              = Tr[h sigma_2(R,l)] - E_L/L,
    |delta(R,l)| <= M ||X-rho_ABC||_1.

This verifies claims 2 and 3 and their extension. The finite sequence of whole states need not be consistent before the moving channel has passed a window, which is why stabilization is used.

The source's l range is empty at L = 2. Define the L = 2 trial state separately as the product of rho_1 on Z. Its defect is at most 2M, so K >= 8 absorbs this length. Equivalently permit l = 1 with empty memory. No limiting argument at L = 2 is needed.

## 3. Lemma P prime and its SDP [Proved reduction; assumed estimate]

Define

    p_L(h,rho) = min_(2<=l<=L-1) min_(R CPTP, R(rho_(l-1))=rho_l)
                    Tr[h sigma_2(R,l)] - E_L/L.

Use the product convention above at L = 2. The minimum exists. Each feasible channel set is nonempty, compact and convex in its Choi matrix, and l takes finitely many values. It includes Petz.

[Assumed, open] Lemma P prime. For the normalized projector onto the ground space of every d-dimensional nearest-neighbor ring, p_L <= K(d) M/L^2 with a finite constant independent of h and L. A version for every translation-invariant ground state is stronger; no equivalence between these ground-state quantifiers is assumed.

[Proved] Lemma P prime implies variant (a) and, on a vacuum-compatible CFT sequence, A <= 6 K(d) M/pi. It is no stronger than Lemma P with the same choices of states and lengths. The unqualified marginal inclusion of section 1 does not imply this particular recovery ansatz is adequate.

Here is an explicit SDP. Let b = rank(rho_B), q = rank(rho_BC). Work in eigenbases on these supports, with input eigenvalues w_i and target D_out. Write the compressed input rho_AB as blocks X_(aa') on B, and let

    H_(a'a) = V^dagger (h_(a'a) tensor 1_(remaining l-1 sites)) V,
    C = sum_(a,a') X_(aa')^T tensor H_(a'a).

For J = sum_(i,j) |i><j| tensor R(|i><j|), minimize Tr(CJ) subject to

    J >= 0,
    Tr_out J = 1_b,
    sum_i w_i J_(ii) = D_out.

The objective is linear and real. To verify its indices, R(X) = sum_(ij) X_(ij) J_(ij), and Tr[(X^T tensor H)J] = sum_(ij) X_(ij) Tr[H J_(ij)]. Sum over A indices to recover Tr[h sigma_2].

Support reduction is exact before numerical eigenvalue thresholds. For any positive input Y supported in supp(rho_B), Y <= c rho_B for some c. Hence 0 <= R(Y) <= c rho_BC, so its range is in supp(rho_BC). Linear extension covers all inputs. Conversely a reduced CPTP channel extends by the completion in section 2. The replacement channel has J = 1_b tensor D_out, which is strictly positive on the reduced support and feasible. Slater's condition gives strong SDP duality.

The dual maximizes Tr Y + Tr(Z D_out), subject to

    C - Y tensor 1_q - diag(w) tensor Z >= 0,

with Y and Z Hermitian. This is also a concrete place to search for analytic optimality certificates. A lower dual bound alone does not prove the desired upper bound on p_L, for which a feasible channel is required.

[Proved] If rho_AB, rho_BC and h share an on-site symmetry, averaging R over it preserves feasibility and objective. An optimal channel may therefore be chosen covariant. This is a structural restriction on the optimizer, not an estimate on its energy.

## 4. Rotated and twirled maps [Proved]

With support powers and P the Petz map, define

    R_t(X) = rho_BC^(-it) P(rho_B^(it) X rho_B^(-it)) rho_BC^(it).

Each is CPTP on the input support and maps rho_B to rho_BC. Completion as in section 2 makes it a full channel. Integrating R_(t/2) against beta_0(t) = pi/[2(cosh(pi t)+1)] preserves these properties. Thus rotated maps, their convex mixtures, and the universal twirled map are feasible in the SDP. This is the convention of Junge, Renner, Sutter, Wilde and Winter, equations (18)-(21), https://arxiv.org/html/1509.07127v3. Admissibility here also follows directly from the displayed formula; their recovery inequality is not needed.

## 5. What the entropy argument does and does not prove [Proved audit]

For epsilon = ||X-rho_ABC||_1/2, data processing and conditional-entropy continuity imply

    I(A:C|B)_rho <= 2 epsilon log d + (1+epsilon) h_2(epsilon/(1+epsilon)).

This part of the source is correct. But substituting epsilon <= f(d)/L^2 yields only c <= O(f(d) log L) after inserting the CFT entropy expression. For every fixed finite c this inequality becomes weaker at large L. A family with growing c can have an arbitrarily growing onset length. The source therefore has no contradiction, even accepting its conditional fixed-d constructions. The sentence that a trace-norm proof is impossible is not established.

There is a second distinction. Recovery error on l+1 sites and distance between two-site marginals are different quantities. Section 1 identifies the optimal two-site distance with trace distance. A successful dimension-only marginal theorem would itself be a trace-norm statement. It does not require small error on the full recovered window.

[Heuristic] Energy weighting remains a sensible strategy because it respects slow-sector scales. It is not logically mandatory by the argument above. The CFT expression for the discrete entropy curvature remains a CFT input, not a lattice theorem uniform in h and L.

## 6. Exact recovery optimum for the ferromagnetic Heisenberg mixture [Proved, d = 2]

Take h = -S dot S, with S = sigma/2, M = 1/2. For every L >= 2 the ground space is the totally symmetric spin-L/2 representation. To see this, h = 1/4 - P_swap/2 and each bond is minimized at swap eigenvalue +1. Adjacent swaps generate the permutation group. The common +1 eigenspace is the symmetric subspace and attains the bondwise lower bound. Its normalized projector has every l-site marginal

    rho_l = P_sym,l/(l+1).

For 2 <= l <= L-1, the input support B has spin j = (l-1)/2 and dimension l, while the output support BC has spin k = l/2 and dimension l+1. The target and input are maximally mixed on these irreducible spaces. Within the output support, the first site's spin is J_out/l by permutation symmetry. The input A B state obeys

    rho_AB = 1_(2l)/(2l) + [2/(l(l+1))] sum_alpha S_A^alpha tensor J_in^alpha.

Indeed the projector onto total spin j+1/2 is [2 S_A dot J_in + (j+1)1]/(2j+1), and division by l+1 gives this formula. Tr(S^alpha S^beta) = delta_(alpha,beta)/2 then gives the energy functional

    e(R) = -1/[l^2(l+1)] sum_alpha Tr[J_out^alpha R(J_in^alpha)].

Equivalently the Choi cost is -sum_alpha (J_in^alpha)^T tensor J_out^alpha/[l^2(l+1)]. Spin time reversal unitarily maps (J_in^alpha)^T to -J_in^alpha. Thus the cost is unitarily equivalent to

    J_in dot J_out/[l^2(l+1)].

Addition of angular momenta gives eigenvalues [s(s+1)-j(j+1)-k(k+1)]/2, for s from |k-j| = 1/2 to j+k. The smallest is -j(k+1) = -(l-1)(l+2)/4. Since Tr J = l, every feasible channel satisfies

    e(R) >= -(l-1)(l+2)/[4l(l+1)]
          = -1/4 + 1/[2l(l+1)].

This lower bound is attained. In the time-reversed representation choose J = (l/2) P_(s=1/2). Its partial trace over the output is 1_l by irreducibility and total trace l. Since the input is 1_l/l, its recovered state is the output partial trace divided by l, which is 1_(l+1)/(l+1), again by irreducibility and trace. It is feasible.

It is also the Petz channel. On these supports the formula reduces to

    R_P(X) = [l/(l+1)] P_sym,l (X tensor 1_2) P_sym,l.

Its two Kraus operators transform as a spin-1/2 doublet under the input-conjugate times output representation. Their Choi vectors span the minimal-spin subspace, and trace preservation fixes their sum of projectors to (l/2)P_(s=1/2) after time reversal. Hence its Choi matrix is the one just constructed. Alternatively covariance and the displayed Kraus formula verify its output and energy directly.

We have proved, for every allowed l,

    delta_P(l) = delta_opt(l) = 1/[2l(l+1)],
    min_l L^2 delta/M = L/(L-1),  L >= 3.

The product convention at L = 2 gives delta = 1/4 and L^2 delta/M = 2. Thus both Lemma P and P prime hold on this class for all L with K = 2. For L >= 3 the sharp uniform constant for this ansatz is 3/2. The actual Casimir lower bound on this frustration-free class has K = 0, since E_L = L e_0 = -L/4. The positive recovery defect is an artifact of this ansatz.

Consistency checks. At l = 2 the defect is 1/12. At l = L-1 and L = 4,6,8,10 the scaled values are 4/3,6/5,8/7,10/9, agreeing with the reproduced run. At l approximately L/2 the scaled value tends to 4 instead of 1. Multiplying h by a positive scalar multiplies delta and M equally. The statement holds after a uniform on-site unitary, a positive scaling, a scalar shift or a gauge change. For even L, the staggered pi rotation about z maps XXZ at Delta = -1 to this ferromagnet. Its normalized ground-space projector and the recovery optimization transform covariantly, so the same formulas hold there. Shifting the staggered pattern is a uniform spin rotation, under which the symmetric mixture is invariant. No all-d symmetric-representation extension is claimed here. The smallest nontrivial local dimension is 2; d = 1 is trivial.

## 7. Candidate energy mechanisms [Proved identities; open estimates]

(a) Ground-state stationarity. If rho is any ground-space state and O is Hermitian, unitary conjugation gives

    Tr rho i[H_L,O] = 0,
    Tr rho [O,[H_L,O]] = 2 Tr rho O(H_L-E_L)O >= 0.

For a CPTP map T supported on an interval, the finite change is

    Tr H_L[(T tensor id)(rho)-rho] >= 0.

Only bonds meeting the interval enter its adjoint expression. These are lower bounds on the total energy cost. They do not upper-bound the one boundary bond in delta. Even an on-site rotation can raise a product ground state's local energy by order M, with no L suppression. Therefore stationarity alone supplies neither a small parameter nor an L^-2 upper bound. A proof using it must add a channel-specific cancellation or a scale-dependent estimate. Averaging over translated boundaries does not improve a translation-invariant scalar cost.

Applying the recovery on a ring and appending a site also changes the other boundary of its interval. Internal output bonds have the prescribed rho_2, but two boundary terms remain. The variational principle bounds E_(L+1) from above by this trial energy; it does not upper-bound either of those terms individually.

(b) Modular flow. The exact formula contains rho_BC^(1/2) and rho_B^(-1/2), or exponentials of K_BC and K_B separately. Noncommutativity prevents replacement of their product by exp[-(K_BC-K_B)/2]. A Baker-Campbell-Hausdorff expansion would require control of commutators and domains. Reduced states of product or frustration-free ground states already have kernels, so their full modular Hamiltonians are unbounded even at d = 2. Restriction to supports repairs the map but does not give local bounds on the inverse eigenvalues or commutators. Thus a uniformly bounded full-space modular mismatch is false under the stated hypotheses. A weaker energy-pairing estimate on the supported product is not refuted, but remains the missing estimate. CFT locality and its inverse-velocity weights are heuristic input here. Near-critical massive response prevents replacing that input by the excluded local energy functionals.

(c) Longer-ring or open-chain dynamics. The replacement channel R(X) = Tr(X)rho_l is always feasible, but its frozen bond is rho_1 tensor rho_1. For antiferromagnetic Heisenberg its energy is zero, and its defect is -E_L/L, of order M. So a prescription that discards the old correlations fails even at d = 2. A channel prepared from another ring generally maps its own shorter marginal to its own longer marginal, not the required rho_(l-1) to rho_l. Equality of a ground energy does not specify these density matrices. An additional channel-interpolation theorem is needed to make such a proposal feasible. Energy variational comparisons alone only supply the familiar order-M seam cost. This obstructs these specific arguments, not all possible uses of ring dynamics.

(d) Toolbox. The supplied toolbox has no uniform energy-weighted recovery theorem. Thermal recovery estimates have constants depending badly on inverse temperature, and their zero-temperature limit is uncontrolled here. LSM twists compare states on the same ring, without the prescribed marginal constraint. Bethe asymptotics at fixed parameters do not provide uniformity at Delta -> -1. Perturbative RG does not cover all d = 2 bonds. The frustration-free energy bound is applicable and yields K = 0 for that energy problem, but it does not make Petz defect vanish, as section 6 proves. Free momentum smoothing proves a different trial-state bound, not Lemma P prime automatically.

[Proved] Classical diagonal ground-space mixtures satisfy P prime with zero defect for L >= 3. At l = 2 the channel measures the input symbol b, leaves b unchanged and appends c with weights rho_2(b,c)/rho_1(b) on nonzero diagonal entries. It maps rho_1 to rho_2 and preserves the old A,b marginal. Translation invariance makes rho_2 a stationary classical edge distribution. This extends to every d. For d = 2 and L = 2 a nonproduct classical extension can also be used directly for the energy theorem, even though our product convention for P prime need not have zero defect.

[Assumed, single remaining sufficient lemma] The general dimension-only upper bound on p_L in section 3 remains open, including at d = 2. None of the identities in this note supplies that bound.


## 8. Uniform Casimir bound for the quadratic class [Proved, smallest d = 2]

The source's momentum-box construction gives a useful explicit M-normalized consequence. Suppose the bond is quadratic in m on-site fermionic modes, with on-site even terms h_on and cross-site hopping or pairing h_cross. Let P be on-site fermion parity. Global parity averaging of a minimizing gauge representative preserves h and cannot increase its norm. Thus one may take the one-site gauge operator even. After centering this representative h_tilde has norm M and

    h_cross = [h_tilde - (P tensor 1) h_tilde (P tensor 1)]/2.

The on-site and gauge terms cancel in this difference, while every cross-site quadratic term changes sign. It follows that ||h_cross|| <= M and |e_cross| <= M in every state.

For a Gaussian eigenstate in a physical Jordan-Wigner parity sector, extend the discrete covariance symbol as a constant on momentum boxes of width 2pi/L centered on each allowed momentum. Positivity of the symbol is preserved. In the paired case the inversion-symmetric momentum grid and conjugate boxes also preserve its particle-hole constraint. Direct Fourier integration gives

    G_infinity(0) = G_L(0),
    G_infinity(plus or minus 1) = sinc(pi/L) G_L(plus or minus 1).

Therefore e_trial = e_on + sinc(pi/L)e_cross and

    e_0 - E_L/L <= [sinc(pi/L)-1] e_cross
                <= [1-sinc(pi/L)] M <= pi^2 M/(6L^2).

The physical spin ground space can be resolved into Gaussian occupation eigenstates in its parity sectors. Mixtures are treated by applying this construction to each component and averaging the resulting infinite-chain states. Thus variant (a) holds with K = pi^2/6 for this quadratic class, at d = 2 and for its m-mode extensions d = 2^m, for all L >= 2 with both physical Jordan-Wigner parity sectors included. This bound is uniform in gaps, velocities and crossover scales. It covers XX, TFIM and C106. It is a trial-state proof of variant (a), not a proof that the recovery-restricted optimum p_L obeys this constant.

Checks. For XX, e_on = 0 and e_0 = -1/pi, so e_0 <= sinc(pi/L) E_L/L rearranges to E_L >= -csc(pi/L), with equality on even rings in the reproduced checks. The weaker M-normalized inequality still holds at L = 2. A pure on-site quadratic bond has e_cross = 0 and no defect. For independent quadratic stacks the same parity argument uses the total local parity and the full optimized M, without counting modes. The induced asymptotic bound A <= pi M is weaker than the previously proved sharp free-family bound A <= 2M in REPORT section 9.71; its purpose here is the uniform finite-L statement.
