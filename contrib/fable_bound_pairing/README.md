# The stiffness pairing: Kohn stiffnesses are bounded by the bond, and what that gives for A

Author: Claude (Anthropic, Fable 5.1), 2026-10-29, as a sub-agent of the coordinating session. Closed book apart from the files named in Section 8. Numerics are numpy and scipy only, in this folder, all run locally in minutes (nothing here needed the cluster). No Codex run: the lattice part is elementary and the obstruction found at the end is an open question rather than a derivation. Every claim is tagged proved, numerical, heuristic, assumed or recalled. Natural logarithms unless a base is written.

## Claim

Setting as in the background note. H = sum_x h_{x,x+1} is a one-site translation-invariant nearest-neighbour chain on (C^d)^L, with infrared a finite set of unitary CFT sectors (c_a, v_a), Casimir coefficient A = sum_a c_a v_a, and M the gauge-optimized half-spread of the bond. The target is A <= pi M log_2 d.

The free-fermion proof of A <= 2M (previous folder) rested on the f-sum step sum_F |v_F| <= 2 pi |<h_hop>|: the Fermi velocities are controlled by the ground-state expectation of the bond itself. This note works out the interacting version of that step. The lattice half generalizes completely and without any crossover scale; the CFT half pays for the generality by weighting the velocities with Luttinger parameters (dressed charges) rather than central charges.

1. [Proved] Stiffness-bond inequality. Let q be a Hermitian one-site operator whose total sum_x q_x is conserved, with the bond chosen in the charge-conserving gauge class (this costs nothing: the gauge optimum can be taken to commute with q (x) 1 + 1 (x) q). Let D_q = lim_L (1/L) d^2 E_0/d theta^2 at theta = 0 be the Kohn stiffness for a uniform twist theta per bond. Then for the centred representative h~ (spectrum in [-M, M]),

       D_q  <=  <h''>_0 := -< [q, [q, h~]] >_0   <=   || [q, [q, h~]] ||   <=   sigma_q^2 M ,

   where q acts on the first site of the bond, the expectation is in the infinite-chain ground state, and sigma_q is the largest charge transferred across the bond by h (at most the spread of the spectrum of q). The first inequality is the variational (f-sum) step, the last is Bernstein's inequality for the entire function theta -> e^{i theta q} h~ e^{-i theta q}, which is of exponential type sigma_q and bounded by M on the real axis.

2. [Proved] Torus budget. If an orthonormal set of commuting charges q^alpha (traceless, tr q^alpha q^beta = delta on C^d) is conserved, then

       sum_alpha D_alpha  <=  < Gamma h~ + h~ Gamma - 2 sum_alpha q^alpha h~ q^alpha >_0  <=  4 M ,

   with Gamma = sum_alpha (q^alpha)^2 (x) 1, 0 <= Gamma <= 1. For the maximal torus (all d - 1 diagonal charges conserved) the middle operator is 2 (h~ - P_1 h~) with P_1 the pinching onto the site-1 diagonal, and the budget is 2 || h~ - P_1 h~ ||, which is at most 4M and at most 2M at d = 2. The budget does not grow with the number of charges or with d.

3. [Proved] Non-abelian budget. If a compact Lie algebra G acts on the site and the bond is G-invariant, then with C_G h~ = sum_alpha [t^alpha, [t^alpha, h~]] the Casimir superoperator (generators orthonormal),

       sum_alpha D_alpha  <=  -< C_G h~ >_0  <=  || C_G h~ || .

4. [Proved given the standard Luttinger and WZW dictionary, Section 2] Consequences for A.
   (a) A single Gaussian sector with Luttinger parameter K (defined by D_q = v K/pi): v <= (pi/K) min( <h''>_0 , sigma_q^2 M ). On the XXZ line this is A <= 2 (pi - gamma) M, within a factor 2 of the exact A = pi sin(gamma)/(2 gamma) for all |Delta| < 1, and it is an equality in the f-sum step at the free point.
   (b) A complete Gaussian infrared with dressed-charge rows Z_a in an orthonormal conserved torus: sum_a v_a |Z_a|^2 <= 4 pi M, hence A <= 4 pi M / min_a |Z_a|^2. There is no d in this bound.
   (c) An SU(N)_k sector realized by an SU(N)-invariant bond: c v <= pi || C_SU(N) h~ || / (k + N).
   (d) Explicitly, for the Heisenberg chain and for every SU(N) swap chain h = P, all three routes give A <= 2 pi M. For N >= 4 this is at least as strong as the conjectured pi M log_2 N, and strictly stronger for N >= 5. For the Hubbard chain the same mechanism is C101's stiffness proof (Q <= 3/(2 pi)).

5. [Proved, within the CFT] The energy-current boost H -> H + theta J_E carries no Casimir information: the CFT ground state is annihilated by J_E, the first-order change of the chiral velocities v_a -> v_a (1 +- theta v_a) cancels in A, and the second-order response is a non-universal lattice susceptibility.

6. [Proved consequence of 2 and the dictionary] Since the stiffness budget is d-independent, a chain can approach A = pi M log_2 d at large d only through sectors whose Luttinger weight |Z_a|^2 is of order 4/log_2 d or smaller, through non-abelian sectors with large level (where c v / D = 2 pi/(k+2) for SU(2)_k), or through neutral sectors that no conserved charge sees. The SU(N) swap chains, the only known family with A/M of order one at large d, stay at A <= 2 pi M.

No bound of the form A <= C M f(d) is proved for general chains. The strongest unconditional inequality here is item 2, sum_alpha D_alpha <= 4M, and the strongest conditional one is item 4(b), sum_a v_a |Z_a|^2 <= 4 pi M. What separates these from the conjecture is a lower bound on Luttinger parameters in terms of d, which I did not find and which may not exist.

## Status

Items 1 to 3 are proved at the level of a variational inequality plus Bernstein's inequality and elementary operator algebra, for finite rings and then in the limit; the only assumptions are conservation of the charge and existence of the stiffness limit. Item 4 uses the textbook Luttinger-liquid and WZW identifications of Kohn stiffnesses (Section 2) as inputs; these are recalled, not proved, and checked numerically here at the SU(2)_1, SU(3)_1 and SU(2)_2 points. Item 5 is a CFT computation. Item 6 is a corollary. All inequalities are checked by exact diagonalization on rings up to L = 18 (d = 2) and L = 12 (d = 3), Section 5.

## How it was established and checked

- Author: Claude (Anthropic, Fable 5.1), 2026-10-29.
- AI assistance: none beyond the author. Nothing was sent to Codex.
- Checks: `stiffness_ed.py` computes D_L by a symmetric finite difference of the ring energy with the twist on one bond (two step sizes, Richardson-combined; the two raw values agree to 1e-6), the bond expectation <h''>_L in the same ground state, the operator norms, and the Casimir fit for A. The chain of inequalities D <= <h''> <= ||h''|| <= sigma^2 M holds at every size for every chain run, and the f-sum step is an equality to 1e-6 on the XX chain where the proof predicts it. `casimir_budget.py` checks the operator inequalities of items 2 and 3 on 900 random bonds and evaluates the budgets for the reference bonds. `pinch_constant.py` probes the sharp constant of the pinching inequality. Nobody else has rerun these.

## 0. Setting and conventions

[Definition] A one-site charge is a Hermitian q on C^d. It is conserved if [Q, H] = 0 for Q = sum_x q_x. For a nearest-neighbour bond this is implied by [q (x) 1 + 1 (x) q, h] = 0, and conversely if Q is conserved there is a gauge representative of h that commutes with the bond charge q (x) 1 + 1 (x) q (the commutator of the bond charge with h is a telescoping one-site difference, which a gauge change removes). All statements below use such a representative.

[Proved, gauge compatibility] The gauge optimization in M can be restricted to the charge-conserving class without loss. If h_* = h + o (x) 1 - 1 (x) o - b has half-spread M, average o over the U(1) generated by q: o_bar = integral dphi/(2 pi) e^{i phi q} o e^{-i phi q}. Since h commutes with the bond charge, h + o_bar (x) 1 - 1 (x) o_bar - b = integral dphi/(2 pi) e^{i phi(q (x) 1 + 1 (x) q)} h_* e^{-i phi (...)} has norm at most M, and o_bar commutes with q. The same averaging over a compact group G works for G-invariant bonds.

[Definition, twist] For a conserved q and a ring of length L, H(theta) = sum_x e^{i theta q_x} h_{x,x+1} e^{-i theta q_x}. Conjugating by exp(i theta sum_x x q_x) moves the whole twist to the bond (L-1, 0), where it becomes Phi = L theta, so E_0(theta) is the ground energy with a Peierls phase Phi on one bond. The Kohn stiffness is D_q = lim_L (1/L) E_0''(0) = lim_L L d^2 E_0/dPhi^2. In the ring sequences used here theta = 0 is the minimum of E_0(theta) and D_q >= 0; the upper bounds below do not need this.

[Recalled, dictionary; Section 2] For a Gaussian sector with velocity v coupled to q with dressed charge z, D_q = v z^2/pi. For a single sector one writes z^2 = K: K = 1 for free fermions with q = n, and K = 1/2 for the Heisenberg chain with q = S^z. For SU(N)_k with generators t^alpha orthonormal in tr(t t) = 1, D_alpha = k v/pi for every alpha.

## 1. The lattice inequalities

[Proved, Lemma 1 (variational step)] Fix L. The ground state |0> of H(0) is non-degenerate in the vacuum sector, and E_0(theta) <= <0| H(theta) |0> for all theta with equality at theta = 0. Both sides are smooth in theta (the right side is entire; the left is analytic near 0 by non-degeneracy). Two functions that agree at a point, one dominating the other, have E_0''(0) <= d^2/dtheta^2 <0|H(theta)|0> at 0. The right side is L times <0| h''_{x,x+1} |0> with h'' = d^2/dtheta^2 [e^{i theta q} h e^{-i theta q}]_{theta=0} = -[q, [q, h]], q on the first site of the bond, by translation invariance. Dividing by L and taking the limit, D_q <= <h''>_0 = -< [q,[q,h]] >_0. In the charge-conserving gauge class the operator [q,[q,h]] is gauge independent: a gauge term o (x) 1 - 1 (x) o with [o, q] = 0 is killed by ad_q. So h can be replaced by h~.

The same computation in second-order perturbation theory gives (1/L) E_0'' = <h''> - (2/L) sum_{m>0} |<m|H'|0>|^2/(E_m - E_0), so the gap between D_q and <h''>_0 is the regular part of the conductivity, nonnegative. This is the f-sum rule; for a single free band the ground state is an eigenstate of the current H' = sum_k v_k n_k, so H' has no off-diagonal matrix elements out of |0> and the regular part vanishes, which is why the XX chain saturates the step.

[Proved, Lemma 2 (Bernstein step)] Let phi(theta) = e^{i theta q} h~ e^{-i theta q}, q acting on site 1. Decompose h~ = sum_m h_m into eigen-operators of ad_q, ad_q h_m = m h_m, where m runs over differences lambda - lambda' of eigenvalues of q connected by h~; let sigma_q = max |m|. Then phi(theta) = sum_m e^{i m theta} h_m extends to an entire function of exponential type sigma_q, and ||phi(theta)|| = ||h~|| = M for real theta because conjugation is unitary. For unit vectors u, w the scalar function f(theta) = <u| phi(theta) |w> is entire of exponential type sigma_q and bounded by M on the real axis, so Bernstein's inequality gives |f''(0)| <= sigma_q^2 M. Taking the supremum over u, w: ||phi''(0)|| = ||[q,[q,h~]]|| <= sigma_q^2 M. For integer-spaced charges phi is 2 pi-periodic and this is the trigonometric Bernstein inequality.

Sharper, for the d = 2 chains: the double commutator kills every charge-diagonal component of the bond, so for XXZ at any Delta, including |Delta| > 1 at finite magnetization, ||[S^z,[S^z,h~]]|| = ||S^x S^x + S^y S^y|| = 1/2 regardless of Delta, whereas M grows with |Delta|. The stiffness is controlled by the charge-moving part of the bond alone.

[Proved, Lemma 3 (torus budget)] Let q^alpha = diag(lambda^alpha) (x) 1 with sum_a lambda^alpha_a = 0 and sum_a lambda^alpha_a lambda^beta_a = delta_{alpha beta}, and let G = sum_alpha lambda^alpha (lambda^alpha)^T, an orthogonal projector on R^d. For any operator X, sum_alpha [q^alpha,[q^alpha, X]] = sum_alpha ((q^alpha)^2 X + X (q^alpha)^2 - 2 q^alpha X q^alpha) = Gamma X + X Gamma - 2 sum_alpha q^alpha X q^alpha with Gamma = sum_alpha (q^alpha)^2 = diag(G_aa) (x) 1, so 0 <= Gamma <= 1. The map X -> sum_alpha q^alpha X q^alpha is R^dag (1 (x) X) R with R = sum_alpha q^alpha (x) e_alpha, and ||R||^2 = ||Gamma|| <= 1, so its norm is at most 1. Hence || sum_alpha [q^alpha,[q^alpha, h~]] || <= 4 ||h~|| = 4M, and with Lemma 1 applied to each conserved q^alpha, sum_alpha D_alpha <= -< sum_alpha [q^alpha,[q^alpha,h~]] >_0 <= 4M. In components, the superoperator acts on the site-1 matrix unit |a><b| (x) Y by multiplication with w_ab = G_aa + G_bb - 2 G_ab = |G e_a - G e_b|^2 <= 2, and w_aa = 0. For the maximal torus G = 1 - (1/d) J, w_ab = 2 (1 - delta_ab), and the superoperator is 2 (1 - P_1) with P_1 the pinching that keeps the components diagonal on site 1. Then sum_alpha D_alpha <= 2 < P_1 h~ - h~ >_0 <= 2 ||h~ - P_1 h~||. Since the pinching has norm 1, ||h~ - P_1 h~|| <= 2M; at d = 2 it is the off-diagonal block of a 2 x 2 block matrix and ||h~ - P_1 h~|| <= M. [Proved, lower bound] The bond (J - (d/2) 1) (x) 1 gives ||h~ - P_1 h~|| / ||h~|| = 2 (1 - 1/d), so the constant 2 cannot be improved uniformly in d. [Numerical] Random search finds nothing above 2(1 - 1/d) (`pinch_constant.py`, values 1.000, 1.311, 1.097 at d = 2, 3, 4 against 1, 4/3, 3/2; the optimizer is weak at d = 4). If 2(1 - 1/d) is the sharp constant, the maximal-torus budget is 4 (1 - 1/d) M.

[Proved, Lemma 4 (non-abelian budget)] If the bond is invariant under a compact G with orthonormal generators t^alpha acting on the site, each t^alpha is conserved and Lemma 1 applies to each: sum_alpha D_alpha <= -< C_G h~ >_0 <= ||C_G h~|| with C_G = sum_alpha ad_{t^alpha}^2 the (positive) Casimir superoperator on site 1. C_G acts on the component of h~ transforming in the irreducible representation R of the site-1 adjoint action as multiplication by C_2(R). For spin-1/2 sites End(C^2) = 1 + adjoint and C_G h~ = 2 (h~ - h~^{(0)}) with h~^{(0)} the site-1 scalar part; for fundamental SU(N) sites C_G h~ = 2N (h~ - h~^{(0)}) in the tr(t t) = 1 normalization. The Casimir does not preserve centring, so ||C_G h~|| is C_2(R) times ||h~ - h~^{(0)}||, which lies between M and 2M.

## 2. The CFT dictionary (inputs, recalled)

[Recalled, standard] For a Luttinger liquid of r bosons with velocities v_a and dressed-charge matrix Z (row a gives the coupling of boson a to the charge directions alpha), the stiffness matrix is D_{alpha beta} = (1/pi) sum_a v_a Z_{a alpha} Z_{a beta}, and the compressibility matrix is kappa = (1/pi) Z^{-1} V^{-1} Z^{-T}. For one boson and one charge, D = v K/pi and kappa = K/(pi v), so v^2 = D/kappa and K^2 = pi^2 D kappa. The normalization is fixed by free fermions, D = v/pi with q = n (K = 1), and by the Heisenberg chain, D = v/(2 pi) = 1/4 with q = S^z (K = 1/2). Checked here: the XX chain gives D_L = <h''>_L = 1/pi to 1e-4 at L = 18 and the Heisenberg chain D_L = 0.2815 at L = 18 decreasing logarithmically toward 1/4 (Section 5).

[Recalled] For an SU(N)_k WZW sector with the site carrying the fundamental (k = 1) or, for SU(2), spin S = k/2, the stiffness tensor on the Lie algebra is an invariant form, hence proportional to the Killing form: D_alpha = k v/pi for tr-orthonormal generators, c = k (N^2 - 1)/(k + N). In SU(2) spin normalization (generators S^alpha), D_{S^z} = k v/(2 pi). Checked here: SU(3)_1 swap chain D_L -> 0.667 = v/pi with v = 2 pi/3 (0.750 at L = 12, decreasing; both Cartan directions give the same D_L to all digits, as invariance requires), and the Takhtajan-Babujian spin-1 chain D_L = 2.53 at L = 12 decreasing toward k v/(2 pi) = 2 with v = 2 pi (Section 5). Convergence at these two points is slow because of the marginal current-current operator.

## 3. Consequences for A

[Proved given Section 2] (a) Single Gaussian sector: A_a = v_a = pi D_q/K <= (pi/K) <h''>_0 <= (pi/K) ||[q,[q,h~]]|| <= pi sigma_q^2 M/K. XXZ, Delta = cos gamma, K = pi/(2(pi - gamma)): A <= 2(pi - gamma) ||S^x S^x + S^y S^y|| = (pi - gamma) = 2(pi - gamma) M. Against the exact A = pi sin(gamma)/(2 gamma): the ratio A/bound is 1/2 at Heisenberg, 2/pi at XX, and tends to 1/2 at the ferromagnetic point; the bound is never off by more than 2 on the whole critical line. (Table in Section 5.)

(b) Complete Gaussian infrared on a conserved orthonormal torus: sum_a v_a |Z_a|^2 = pi tr D <= 4 pi M (Lemma 3), hence A = sum_a v_a <= 4 pi M / min_a |Z_a|^2. For SU(N)_1 realized on fundamental sites, |Z_a|^2 = 1 for the N - 1 Cartan bosons, and the sharper maximal-torus form gives A <= 2 pi ||P - P_1 P|| = 2 pi for the swap chain (||P - P_1 P|| = 1 since P - P_1 P = sum_{a != b} |ab><ba| has spectrum {1, -1, 0}). So A <= 2 pi M for every SU(N) swap chain, with the exact value 2 pi (1 - 1/N).

(c) SU(N)_k: sum over all N^2 - 1 generators gives (N^2 - 1) k v/pi <= ||C_SU(N) h~||, and c v = k (N^2 - 1) v/(k + N) = pi (sum_alpha D_alpha)/(k + N) <= pi ||C_SU(N) h~||/(k + N). For the swap chain ||C h~|| = 2N ||P - 1/N|| = 2(N + 1), so c v <= 2 pi again. For the Heisenberg chain ||C h~|| = 2 ||S.S|| = 3/2 and c v <= (2 pi/3)(3/2) = pi = 2 pi M. For the Takhtajan-Babujian spin-1 chain (M = 3, ||C h~|| = 16, k = 2) the bound is c v <= 8 pi = 25.1 against the exact c v = (3/2)(2 pi) = 9.42 and the conjectured pi M log_2 3 = 14.9; here the route is weaker than the conjecture because C_2 of the rank-2 tensor components of the bond is 6 while the stiffness only sees 2S = 2 units of charge.

(d) Where the route beats the conjecture [proved given the dictionary]: for any chain with SU(N)-invariant bond, N >= 4, whose infrared is SU(N)_1 (equivalently any complete Gaussian infrared with all |Z_a|^2 >= 1 on a maximal torus), A <= 2 pi ||h~ - P_1 h~|| <= 4M (1 - 1/d)-type constants times pi, which is below pi M log_2 d from d = 16 on with the proved constant 4 and from d = 5 on for the swap chains. The conjecture's logarithm is therefore not the right growth for Gaussian chains with order-one Luttinger parameters; it can only be approached by small |Z_a|^2 (item 6 of the Claim).

[Proved, exact powers] All lattice quantities here are first powers of M: the stiffness is a ground-state response at zero temperature, so unlike the thermal moments of C100 (m_1 <= 8 M^3 q^2, needing a matching temperature) no crossover scale enters. The entire price of the route is the factor 1/K.

## 4. Why the pairing route does not reach the conjecture

[Proved] The velocities enter every stiffness weighted by dressed charges. The conjecture needs them weighted by central charges. For Gaussian sectors the conversion factor is 1/|Z_a|^2, for SU(2)_k it is 2 pi/(k+2) per unit stiffness, and for neutral sectors (Ising, parafermions, cosets with no conserved on-site charge coupling to them) there is no stiffness at all. The lattice gives no lower bound on |Z_a|^2: the only exact statement available is the compressibility, kappa = K/(pi v), and a lower bound on kappa would be an upper bound on the curvature of the energy density e''(n), which is not controlled by any norm (e is convex with range 2M, but its second derivative is local). I did not find a lattice argument bounding e'' by M, and the relation v^2 = D/kappa shows that this is the whole gap: D <= O(M) is proved here, and v <= O(M) would follow from 1/kappa <= O(M).

[Proved, CFT] The energy-current boost. With J_E = sum_x j_x, j_x = i[h_{x-1,x}, h_{x,x+1}], the CFT limit of J_E is sum_a v_a^2 P_a with P_a the sector momentum, so H + theta J_E has chiral velocities v_a (1 +- theta v_a), whose sum per sector is unchanged. The ring ground energy E_0(L, theta) = e(theta) L - (pi/(6L)) sum_a c_a v_a + o(1/L) + O(theta^2) has no first-order dependence on theta in the Casimir term, the ground state being annihilated by the CFT current, and the bulk term e(theta) = e_0 - (theta^2/2) chi_J + ... involves the lattice susceptibility chi_J of J_E, which is not universal (the lattice J_E differs from the CFT current by irrelevant operators that do not annihilate the ground state). The CFT instability at theta = 1/v_max (left-moving modes of the fastest sector becoming static) is a statement about excited states and leaves no trace in the ground-state energy at finite L. Nothing in this family measures A.

[Heuristic] The same conclusion applies to any lattice deformation whose CFT limit is a conserved current: conserved-current deformations of a CFT move it along the conformal manifold or rescale velocities symmetrically, and the Casimir coefficient changes only at second order through non-universal mixings. Casimir information is first order in the stress tensor, which is why the two-site pairing Tr[h grad_h A] = A of the previous folders, and not a twist, is the exact identity, and why the only lattice routes that measure A directly are the ring-to-line marginal difference and the energy-current Drude weight in the thermal channel (C100), the latter with its matching temperature.

[Proved given the dictionary, structural consequence] Items 2 and 4(b) of the Claim say that the total stiffness over all conserved charges is at most 4M whatever d is. On the slow-mode constructions of the campaign (c unbounded at fixed d, each extra sector exponentially slow) this is automatically satisfied, since the stiffness of a sector is v_a |Z_a|^2 and v_a is tiny. On the SU(N) chains it caps A at 2 pi M. So within everything the repository knows, the velocity-weighted bound has slack that grows with d, and the stiffness budget is the quantitative reason.

## 5. Numerical checks

[Numerical] `stiffness_ed.py` (output `stiffness_ed.out`), sparse Lanczos in the zero-charge sector, D_L from E_0(+-Phi) with Phi = 0.02 and 0.04 on one bond, Richardson-combined (the two raw values agree to 1e-6 in every row). <h''>_L is the ring expectation of -[q,[q,h~]] on one bond. Infinite-volume D for the integrable chains from the dictionary; A from the fit E_0/L = e_0 - a/L^2 + b/L^4 over the four largest sizes.

| chain | d | M | ||[q,[q,h~]]|| | D_L (largest L) | D exact | <h''>_L (largest L) | D_L/<h''>_L | A (fit) | A exact | bound (pi/K) ||h''|| |
|---|---|---|---|---|---|---|---|---|---|---|
| Heisenberg, L = 18 | 2 | 1/2 | 1/2 | 0.2815 (log decrease) | 1/4 | 0.2971 | 0.948 | 1.578 | 1.5708 | 3.1416 = 2 pi M |
| XXZ Delta = 0.5, L = 18 | 2 | 1/2 | 1/2 | 0.3121 | 0.3101 | 0.3141 | 0.993 | 1.2994 | 1.2990 | 2.0944 |
| XX, L = 18 | 2 | 1/2 | 1/2 | 0.3199 | 1/pi | 0.3199 | 1.0000 | 0.9999 | 1 | 1.5708 |
| XXZ Delta = -0.5, L = 18 | 2 | 1/2 | 1/2 | 0.3119 | 0.3101 | 0.3124 | 0.998 | 0.6496 | 0.6495 | 1.0472 |
| SU(3) swap, Cartan, L = 12 | 3 | 1 | 2 | 0.7500 (decreasing) | 2/3 | 0.7890 | 0.951 | 4.215 | 4.1888 | 2 pi (torus) |
| spin-1 TB, S^z, L = 12 | 3 | 3 | 5.464 | 2.528 (decreasing) | 2 | 2.815 | 0.898 | 9.535 | 9.4248 | 25.13 (Casimir) |
| spin-1 ULS, S^z, L = 12 | 3 | 1 | 4 | 1.5001 | 4/3 | 1.5781 | 0.951 | 4.215 | 4.1888 | same as SU(3) |

Every row obeys D_L <= <h''>_L <= ||[q,[q,h~]]|| <= sigma_q^2 M. The XX row is an equality in the first step at every L (the regular part vanishes for free fermions), as the proof predicts. The ULS row is the SU(3) swap chain in spin-1 language (h = P + 1): D and <h''> are exactly twice the Cartan values because S^z = sqrt(2) times the orthonormal Cartan charge, a check of the normalization in Lemma 3. The two SU(3) Cartan directions give identical D_L to all digits, as SU(3) invariance requires. The TB ground energy e_0 = -4.0000 is the exact Takhtajan-Babujian value in these units.

[Numerical] `casimir_budget.py` (output `casimir_budget.out`): on 900 random complex Hermitian bonds at d = 2, 3, 4 the ratio ||h~ - P_1 h~||/||h~|| never exceeds 1.061 and the maximal-torus Casimir ratio never exceeds 2.122 (proved bounds 2 and 4); the explicit bond (J - 1) (x) 1 at d = 3 gives 4/3. Budgets for the reference bonds: ||h~ - P_1 h~|| = M for Heisenberg, all XXZ, all SU(N) swap chains and ULS, and 1.124 M for TB; ||C_SU(N) P|| = 2(N + 1) for N = 2, 3, 4, 5, 8; along the spin-1 bilinear-biquadratic arc cos(theta) S.S + sin(theta)(S.S)^2 the single-twist norm ||[S^z,[S^z,h~]]||/M runs from 1.82 (TB) through 0.94 (theta = 0) to 4 (ULS) while ||h~ - P_1 h~||/M stays within [0.94, 1.12].

## 6. What remains

1. The missing lemma is a lower bound on Luttinger weights. In the single-sector language: is there f(d) > 0 with K >= f(d) for every nearest-neighbour chain with local dimension d whose infrared is a single Luttinger liquid at finite velocity? Equivalently, is the curvature of the ground energy density, e''(n) = pi v/K, bounded by C(d) M? Known examples do not decide it: the XXZ line has K >= 1/2, the t-V chain reaches K = 1/4 only as M diverges. A proof would give A <= 4 pi M/f(d) for all complete Gaussian chains, and a counterexample with K -> 0 at fixed d, M and v would show that the stiffness pairing cannot reach any d-only bound even in the Gaussian class.
2. Neutral sectors (Ising, Z_k parafermions, cosets) carry no stiffness. The free-fermion proof handled Majorana chains because there the Fermi data is directly available; no interacting analogue is in sight. This is the same gap named in REPORT section 8 (a locality budget for neutral stress-tensor sectors).
3. For non-abelian sectors the route loses the ratio C_2(R_max)/k between the Casimir of the highest tensor rank in the bond and the level. A refinement that uses only the current-algebra part of the bond (projecting h~ onto the components of rank 1 before taking the norm) is not valid as written, because Lemma 1 needs the full h''.
4. The sharp constant in ||h~ - P_1 h~|| <= c_d ||h~|| (conjecturally 2(1 - 1/d)) is a finite-dimensional question that could be settled; it would sharpen Lemma 3 to 4(1 - 1/d) M but changes nothing structural.

## 7. Reproducing

Python 3 with numpy and scipy, run from this folder.

- `stiffness_ed.py [quick]`: the table of Section 5. Full run (L <= 18 at d = 2, L <= 12 at d = 3) about 4 minutes on a laptop; `quick` stops at L = 14 and 10 and takes 10 seconds. Output `stiffness_ed.out`.
- `casimir_budget.py`: operator-norm budgets and random-bond checks of Lemmas 3 and 4, seconds. Imports bond constructors from `stiffness_ed.py`. Output `casimir_budget.out`.
- `pinch_constant.py`: Nelder-Mead search for the sharp pinching constant, about a minute. Output `pinch_constant.out`.

No cluster job was needed.

## 8. Sources used

The background note for this round, `contrib/fable_bound_proof/README.md` (the Euler pairing, Proposition 3 on M as a dual norm, the energy-comparison obstruction), `contrib/fable_kappa/README.md` (D_2 unbounded, so gradient-norm arguments fail), `contrib/fable_bound_freefermion/README.md` (the Drude step sum |v_F| <= 2 pi |<h_hop>| and the lifting lemma), `codex/work_C100/RESULT.md` (thermal moments and the energy-current Drude weight, with their matching scales), and `codex/work_C101/PROOF.md` Section 4 (the stiffness matrix inequality D <= K_dia used for Hubbard and t-J, which is Lemma 1 here in matrix form; this note adds the Bernstein, torus and Casimir norm bounds that tie it to M, and the assessment of what the route can and cannot reach). Bethe-ansatz velocities and Luttinger parameters for XXZ (des Cloizeaux-Gaudin, Luther-Peschel), the Kohn stiffness formula and the WZW stiffness normalization are recalled from memory and checked numerically here.

## Relation to the rest of the repository

This closes the "interacting Drude route" named in the background note. Outcome: the lattice half of the free-fermion argument generalizes (Lemmas 1 to 4, no crossover scale, first power of M), the CFT half replaces c_a by Luttinger weights, and the resulting inequality sum_a v_a |Z_a|^2 <= 4 pi M is a different, d-independent constraint which is stronger than the conjecture on the SU(N) swap chains for N >= 5 and silent on neutral sectors. The conjecture A <= pi M log_2 d is untouched for general chains; the precise missing input is a lower bound on Luttinger weights in terms of d. "Where things stand" in the top-level README is unchanged.
