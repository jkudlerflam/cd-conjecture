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
