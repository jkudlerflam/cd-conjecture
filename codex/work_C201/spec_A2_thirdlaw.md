# Task: a uniform decay theorem for the ring energy defect ("crossed third law"), or its refutation

## Common definitions (complete; you have no other context)

Chain. Fix d >= 2 and a Hermitian bond h on C^d ⊗ C^d. The ring Hamiltonian on (C^d)^{⊗L}, L >= 2, is H_L = sum_{j in Z_L} h_{j,j+1}. E_L is its lowest eigenvalue. e_0 = lim_{L} E_L / L exists and equals the infimum of tr(h rho_2) over two-site marginals rho_2 of translation-invariant states on the infinite chain Z.

Gauge and M. Bonds h and h' = h + X ⊗ 1 - 1 ⊗ X + c (X Hermitian on C^d, c real) have the same H_L up to the constant cL. Define M(h) = min_X (lambda_max - lambda_min)/2 of h + X ⊗ 1 - 1 ⊗ X. M is a seminorm on bonds, vanishing exactly on gauge-trivial bonds, and it is the dual norm for differences of two-site states with equal one-site marginals: for two-site states rho, rho' whose left marginals agree and whose right marginals agree, |tr h (rho - rho')| <= M ||rho - rho'||_1.

Marginal sets. Omega_L = set of two-site (nearest-neighbour) reduced states of translation-invariant states on the L-ring. Omega_inf = same for the infinite chain. Both are compact convex. E_L / L = min over Omega_L of tr(h rho), e_0 = min over Omega_inf.

Infrared assumption (physics level, used only to interpret): when the low-energy theory is a sum of decoupled unitary CFT sectors a with central charges c_a and velocities v_a, then on vacuum-compatible lengths E_L = L e_0 - pi A / (6 L) + o(1/L) with A = sum_a c_a v_a.

Target (open). Lemma P: there is K(d) with E_L >= L e_0 - K(d) M / L for all h and all L. Equivalent forms: Omega_L is contained in the (K/L^2)-neighbourhood of Omega_inf in the norm dual to M; and (on doubling) C(d) = sup_{h,L} L (E_{2L} - 2 E_L)/M < infinity. It would give A <= (6/pi) K(d) M.

Known facts you may use without proof.
(F1) Trivial bound: L e_0 - E_L <= 2M for all h, L (cut the ring state at one bond and take tensor powers: the state with two-site marginal ((L-1) rho_2 + rho_1 ⊗ rho_1)/L lies in Omega_inf).
(F2) Lemma P holds for every quasi-free (quadratic fermion, any number of modes, with pairing) chain with K = pi^2/6, for classical (diagonal) and commuting bonds with K = 0, for frustration-free bonds with K = 0, and for the spin-1/2 XXZ chain with Dzyaloshinskii-Moriya term and field for anisotropy in [-1 + delta_0, 1] with an explicit K(delta_0).
(F3) Lieb-Robinson: every CFT sector has v_a <= 4 e M (e = 2.718...).
(F4) For any conserved on-site charge q (H_L commutes with sum_j q_j), the Kohn stiffness obeys D_q <= ||[q,[q,h]]|| <= sigma_q^2 M with sigma_q the spectral spread of q, proved variationally by spreading the twist uniformly. For commuting orthonormal charges sum_alpha D_alpha <= 4M. So sectors carrying on-site charge with order-one weight cannot form a fast tower.
(F5) Free fermions with m orbitals per site (d = 2^m): A <= 2M whatever m is.
(F6) Blocking b sites into one multiplies M and log d by about b and divides A (in block units) by b.
(F7) Known examples: Heisenberg d = 2: L (L e_0 - E_L) -> pi^2/12 with M = 1/2 in the normalization h = S.S. SU(N) symmetric-representation chains show that sup A/M grows slowly with d, consistent with A <= pi M log_2 d.
(F8) Constructions with unboundedly many CFT sectors at fixed d exist (conditionally), but their sector velocities are suppressed like a high power of a small parameter, so A/M stays bounded.

Recorded failed routes and their obstructions (do not repeat without a new ingredient).
(O1) Any local energy functional (smooth modulations sum w_j h_j, gradient norms of A in h, energy-current sum rules, fixed-support marginal bootstraps) receives contributions from nearly critical massive modes at the same order as the conformal term, or loses O(M) per cut.
(O2) Positivity of the spatial transfer operator alone gives E_{2L} - 2E_L >= 0 and nonincreasing, and concavity of E_L in L, but there is an abstract positive countermodel with A unbounded, so crossed-channel positivity data do not suffice.
(O3) Recovery-map (Petz, Fawzi-Renner) periodization of the ring ground state controls the two-site error only through the square root of a conditional mutual information, which gives O(M) for L e_0 - E_L at best, the trivial order. In detail: building a state on Z by appending blocks C of n sites with a recovery channel acting on the last k sites B disturbs one bond per n, with trace-norm error at most sqrt(2 ln 2 I(a:C|B)); pigeonhole gives I <= 2 n log d / L for general states and the CFT value is of order c (n/k)^2, and in both cases the resulting bound on L e_0 - E_L is of order M, independent of n.
(O4) The doubling defect E_{NL} - N E_L is the energy of a Z_N permutation twist in N copies; spreading it by on-site fractional permutations fails for interacting chains because fractional powers of the copy permutation are not symmetries.
(O5) Counting low-lying ring states fails (slow sectors flood fixed windows). The number of infrared sectors is not bounded in terms of d. Crossover lengths can be arbitrarily long.

Style of the answer. Complete arguments with every step shown. Tag each step [rigorous], [rigorous given cited theorem X] with X stated precisely, [numerical], or [physics-level]. A theorem whose hypothesis secretly assumes the conclusion is not acceptable. A sharp, proved obstruction or a counterexample is a good outcome. No probability estimates. Plain declarative prose, no em-dashes.

## The specific question

Define g_d(L) = sup over bonds h on C^d ⊗ C^d with M(h) = 1 of (L e_0(h) - E_L(h)). Equivalently g_d(L) = L^2 times the one-sided Hausdorff distance from Omega_L to Omega_inf in the M-dual norm (h may depend on L). (F1) gives g_d(L) <= 2. Lemma P is g_d(L) <= K(d)/L. Nothing between these is known.

Interpretation. If a single chain had N(L) infrared sectors with velocity of order M and crossover length at most L, one expects L e_0 - E_L of order M N(L)/L. So g_d(L) -> 0 is the statement that the number of fast sectors visible at scale L is o(L), and g_d(L) bounded away from 0 would be an extensive "fast tower". In the crossed (space-as-time) channel, for chains with a positive spatial transfer operator, s_L := L E_{L+1} - (L+1) E_L is the entropy per unit Euclidean time of the crossed Gibbs state at crossed temperature 1/L, it is nonincreasing in L and bounded by 2M, the CFT value is (pi/3) A / L, and lim s_L = lim (L e_0 - E_L). So the question is a third law for the crossed channel.

Prove as much as you can of the following chain, in order of strength, and say exactly where the proof stops.
(T0) Pointwise: for every fixed h, limsup_{L} (L e_0 - E_L) <= 0. (Is even this known or provable in general? Note E_L - L e_0 can be positive and order M on vacuum-incompatible lengths, that direction is not at issue.)
(T1) Uniform: g_d(L) -> 0 as L -> infinity at fixed d, with any explicit rate.
(T2) A power: g_d(L) <= K(d) L^{-p} for some p > 0.
(T3) Lemma P itself, p = 1.
Alternatively refute: exhibit, for some fixed d, bonds h_L with M(h_L) = 1 and L e_0(h_L) - E_L(h_L) >= c > 0 for infinitely many L (equivalently translation-invariant ring states whose two-site marginal stays at distance at least c/L from Omega_inf), or a single h violating (T0).

Ingredients that have not been tried and may help (use, improve or discard them).
(I1) Pure-state second-order identity. Let psi be the ring ground state, H psi = E psi, and let a channel with Kraus operators K_i act on a region R. Writing K_i psi = c_i psi + chi_i with chi_i orthogonal to psi, the energy change is exactly sum_i <chi_i|(H - E)|chi_i>, and sum_i ||chi_i||^2 = 1 - <psi| channel(psi) |psi>. So energy errors are second order in the amplitude error whenever the channel nearly fixes the global pure state, in contrast to the first-order (square-root) loss in (O3). The difficulty is to bound the excitation energies <chi_i|(H-E)|chi_i>/||chi_i||^2 and to use this for a change of topology (ring L to ring 2L, or to Z).
(I2) Site insertion. Lemma P follows from the variational statement E_{L+1} <= (L+1) E_L / L + K M / L^2 for all L >= L_0 together with control of lim (E_L - L e_0) (explain precisely what is needed; note that summing increments alone leaves an additive constant undetermined). The CFT value of the excess is pi A (2L+1) / (6 L^2 (L+1)). (F1)-type arguments give only K M / L. One needs an isometry from the L-ring to the (L+1)-ring that stretches the state uniformly.
(I3) Kinematic freedom. Since the statement is about marginal sets, h is only a witness. One may try to prove the inclusion Omega_L ⊂ Omega_inf + eps(L) L^{-1} B directly for arbitrary translation-invariant ring states, pure or mixed, with eps -> 0. Decide whether this purely kinematic inclusion can hold. The extreme points of Omega_L are exposed by ground states, so no generality is lost, but a kinematic proof cannot use (I1).
(I4) For mixed or for symmetric (permutation-invariant) ring states the de Finetti gap is Theta(1/L) in some norms. Decide whether permutation-invariant states, Dicke states, W-type momentum states, or ground states of L-dependent bonds near a ferromagnetic or other z = 2 point give two-site marginals at distance of order 1/L from Omega_inf, or whether they are within O(1/L^2).

Give complete proofs for whatever you establish. If you prove (T0) or (T1), state all constants. If you find the kinematic inclusion false at rate o(1/L), give the explicit states and the explicit separating bonds h_L with the computation of e_0(h_L) to the needed accuracy. Small exact computations (for example d = 2, L <= 6) are welcome as checks, with the code or the closed forms shown.
