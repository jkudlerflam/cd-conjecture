# The velocity-weighted bound in the quasi-free class: A <= 2 M, sharp, with no log d

## Claim

Setting. A translation-invariant nearest-neighbour quasi-free chain with m fermion modes per site, local dimension d = 2^m, written in Majorana form with n = 2m Majoranas gamma_{x,a} per site,

    H = (i/4) sum_{x,y} gamma_x^T A_{x-y} gamma_y,   A_0 real antisymmetric,  A_1 = B real,  A_{-1} = -B^T,

so that the single-particle Bloch matrix is h(k) = i(A_0 + B e^{ik} - B^T e^{-ik}), an n x n Hermitian matrix with bands eps_j(k), j = 1..n. The number-conserving subclass is H = sum_x [f_x^dag A_0 f_x + f_x^dag B f_{x+1} + h.c.] with A_0 Hermitian and B an arbitrary complex m x m matrix, Bloch matrix h(k) = A_0 + B e^{ik} + B^dag e^{-ik}. The same symbol B is used for the hopping matrix in whichever language is in force; ||B||_1 is its trace norm (sum of singular values) and ||B||_op its largest singular value.

The two quantities compared are the ones of the background note. A = sum_a c_a v_a is the Casimir coefficient, E_0(L) = e_0 L - pi A/(6L) + o(1/L). For a quasi-free chain the low-energy theory is a collection of chiral fermions, one for each zero crossing of a band, and the standard dictionary gives

    A = (1/2) sum over Fermi points of |d eps/dk|          (number-conserving, m complex bands),
    A = (1/4) sum over all zero crossings of the n Majorana bands of |d eps/dk|,

where crossings are counted with multiplicity (an orthonormal basis of ker h(k_F) diagonalizing the compressed velocity matrix P_F h'(k_F) P_F). The two formulas agree on number-conserving chains, where each Fermi point appears twice in the Majorana count. M is the gauge-optimized half-spread of the two-site bond h_{12} on C^d x C^d, M = min over alpha and Hermitian one-site o of (lambda_max - lambda_min)/2 of h_{12} + alpha + o x 1 - 1 x o.

Theorem 1 (bond scale) [proved]. For every quasi-free nearest-neighbour chain, whatever the on-site term A_0,

    M >= ||B||_1          (number-conserving, B the m x m hopping matrix),
    M >= ||B||_1 / 2      (Majorana form, B the n x n real hopping matrix),

with equality when A_0 = 0. The general dual formula is M = max of tr(Delta h_{12})/2 over Hermitian Delta with tr Delta = 0, ||Delta||_1 <= 2 and equal one-site marginals tr_1 Delta = tr_2 Delta.

Theorem 2 (Drude bound) [proved]. For every quasi-free nearest-neighbour chain, in either language,

    A <= pi |<h_hop>_GS| <= pi ||h_hop||_op = pi x (lower bound of Theorem 1) <= pi M,

where h_hop is the hopping part of one bond and <h_hop>_GS its ground-state expectation. This already gives A <= pi M <= pi M log_2 d for every d = 2^m, the conjectured bound with the factor log_2 d removed, and pi is the value reached by the interacting Heisenberg chain at d = 2.

Theorem 3 (sharp constant at d = 2 without pairing) [proved]. For m = 1 number-conserving chains, A <= 2 ||B|| = 2 M, with equality exactly for the XX chain at half filling (any phase of the hopping). More generally A <= 2 rank(B) ||B||_op for any m.

Theorem 4 (sharp constant) [proved; proof found by Codex, checked by me]. For all quasi-free nearest-neighbour chains, in either language,

    sum over crossings of |v| <= 4 ||B||_1,   hence   A <= 2 ||B||_1 <= 2 M   (number-conserving),   A <= ||B||_1 = 2 (||B||_1/2) <= 2 M   (Majorana),

with equality in the first inequality iff B is normal and A_0 B = 0 (so A_0 = 0 when B is invertible; A_0 may act on ker B, where it must be invertible). Equality in A = 2 M needs in addition M = ||B||_1 (or ||B||_1/2), which an on-site term on ker B spoils (it raises the bond spread), so the chains with A = 2 M are exactly A_0 = 0 with B normal: the stacks of XX chains at half filling with arbitrary complex hopping phases, in an arbitrary common on-site unitary basis, together with (Majorana form) decoupled single Majorana chains. So

    A <= 2 M <= (2/pi) pi M log_2 d   for every d = 2^m,

the XX chain is optimal among free chains, and the optimal constant on the quasi-free class is 2, independent of d. Interactions are needed to go from 2 M to pi M.

Interactions [exact Bethe ansatz, recalled; first order, proved]. For the XXZ chain at half filling, M = J/2 for |Delta| <= 1 and A/M = pi sin(gamma)/gamma with Delta = cos gamma. Repulsion pushes A/M above the free value 2, monotonically up to pi at the isotropic point, where both the free-fermion bound A <= pi M and the conjectured pi M log_2 d are saturated. First-order Hartree-Fock reproduces the slope, A/M = 2 + 2 V/(pi t) + O(V^2) with V = J Delta, but extrapolated to Delta = 1 it would overshoot pi (3.27 against 3.14); the exact velocity stays at or below pi M. Nothing at first order in an interaction can push A/M above pi m starting from a free chain, since every free chain has A/M <= 2.

## Status

Theorems 1 to 4 are proved below with all steps. Theorem 4 was found as a conjecture by random sampling and adversarial optimization (with an independent brute-force recomputation of A for every candidate; m <= 3 complex modes, n <= 6 Majoranas), then proved by gpt-6-astra from a self-contained spec; I checked each step of its proof by hand and the two non-obvious identities numerically (`lift_check.py`). The interacting statements are first order plus the exactly solvable XXZ line; no claim is made about strongly interacting chains with m >= 2.

## How it was established and checked

- Author: Claude (Anthropic, Fable 5.1), 2026-10-29, as a sub-agent of the coordinating session.
- AI assistance: the sharp-constant question (Theorem 4) was posed to gpt-6-astra (effort high) through the shim, closed book, with the self-contained spec `spec_sharp_constant.md`; its output is `codex_sharp_constant.out`. The lifting trick and the trace lemma of Section 5 are its. I had derived the J-orthogonality relation (a) of Section 5 and the A_0 = 0 reduction before, and I verified its proof line by line and numerically. Theorems 1 to 3 and Section 6 are mine.
- Checks: the f-sum chain of inequalities behind Theorem 2 is verified numerically on random chains in both languages (`fsum_check.py`), with equality reproduced where the proof predicts it (decoupled bands). The gauge-optimized M is computed by a convex solver (`gauge_M.py`) and agrees with the dual lower bound in every case where the proof says they coincide. Each adversarial optimum is recomputed by a fine k-grid band scan independent of the companion-matrix root finder (`verify_optimum.py`); two early candidates above 2 were rejected by this check and traced to the root finder counting a near-unimodular pair twice.

## Proof

### 1. The free-fermion dictionary for A [standard input]

Near a simple zero crossing k_F of a Majorana band with slope v = eps'(k_F), the modes gamma_k with k near k_F form a chiral Majorana fermion with velocity |v| and chiral central charge 1/2. The particle-hole redundancy gamma_k^dag = gamma_{-k} pairs the crossing (k_F, j) with (-k_F, j') of the band eps_{j'}(k) = -eps_j(-k); a pair is one chiral complex fermion (chiral central charge 1), a self-paired crossing (k_F = 0 or pi) is one chiral Majorana. In both cases the Casimir energy of the ring is E_0 = e_0 L - (pi/12 L) sum over chiral sectors of c_chiral v, so A = (1/2) sum_chiral c_chiral v = (1/4) sum over all crossings of |v|, uniformly. In the number-conserving language each Fermi point of a complex band is one chiral complex fermion, A = (1/2) sum_F |v_F|. At a degenerate crossing (dim ker h(k_F) = p > 1) the p analytic branches through zero have slopes equal to the eigenvalues of P_F h'(k_F) P_F (first-order degenerate perturbation theory, exact for analytic Hermitian families), and each branch is one chiral mode. Checks: XX at half filling (t = 1) gives A = 2, critical Ising (J = g = 1) gives A = 1/2, in agreement with c v.

### 2. The bond scale M

Lemma 2.1 (duality) [proved]. For Hermitian X on C^d x C^d, (lambda_max - lambda_min)/2 = max over density matrices rho, sigma of tr((rho - sigma) X)/2. The gauge term tr((rho - sigma)(o x 1 - 1 x o)) = tr((rho_1 - rho_2 - sigma_1 + sigma_2) o) with rho_1 = tr_2 rho and so on, and the alpha term is tr(rho - sigma) alpha = 0. The function (o, alpha; rho, sigma) -> tr((rho - sigma)(h + alpha + o x 1 - 1 x o))/2 is bilinear, the set of pairs (rho, sigma) is compact convex and the set of (o, alpha) is convex, so Sion's minimax theorem gives

    M = max over (rho, sigma) of inf over (o, alpha) of ... = max over (rho, sigma) with rho_1 - rho_2 = sigma_1 - sigma_2 of tr((rho - sigma) h)/2,

the infimum being -infinity unless the marginal constraint holds. Writing Delta = rho - sigma, which ranges exactly over Hermitian traceless Delta with ||Delta||_1 <= 2, gives the formula in Theorem 1. In particular any single two-site state psi with equal one-site marginals, paired with any other such state phi, gives M >= (<psi|h|psi> - <phi|h|phi>)/2, and any pair of translation-invariant chain states gives M >= (e_max - e_min)/2 in terms of energy densities.

Lemma 2.2 (the hopping bond) [proved]. Number-conserving case: with the singular value decomposition B = sum_s sigma_s |a_s><b_s|, define new modes g_{1s} = a_s^dag f_1 on site 1 and g_{2s} = b_s^dag f_2 on site 2 (unitary changes of basis on each site separately). Then h_hop = f_1^dag B f_2 + h.c. = sum_s sigma_s (g_{1s}^dag g_{2s} + h.c.) is a sum of m commuting dimers, each with one-particle energies +-sigma_s and zero- and two-particle energies 0. So spec h_hop = {sum_s e_s sigma_s : e_s in {-1, 0, 1}} and ||h_hop||_op = ||B||_1. Majorana case: with the real singular value decomposition B = O_1 Sigma O_2^T, the Majoranas gamma~_1 = O_1^T gamma_1, gamma~_2 = O_2^T gamma_2 give h_hop = (i/2) gamma_1^T B gamma_2 = (i/2) sum_s sigma_s gamma~_{1s} gamma~_{2s}, a sum of commuting terms with eigenvalues +-sigma_s/2, so spec h_hop = {sum_s (+-) sigma_s/2} and ||h_hop||_op = ||B||_1/2.

Lemma 2.3 (lower bound, Theorem 1) [proved]. Let |Psi_+> and |Psi_-> be the two-site states with every dimer in its top, respectively bottom, one-particle level (number-conserving: Psi_+- = prod_s (g_{1s}^dag +- g_{2s}^dag)/sqrt2 |0>; Majorana: the joint eigenstates of all i gamma~_{1s} gamma~_{2s} with eigenvalue +1, respectively -1). Then <Psi_+|h_hop|Psi_+> - <Psi_-|h_hop|Psi_-> = 2 ||h_hop||_op. Both are fermionic Gaussian states of definite parity (the product of all i gamma~_{1s} gamma~_{2s} is, up to sign, the total parity). Their one-site reduced states: the correlation matrix restricted to site 1 is <g_{1s}^dag g_{1t}> = delta_{st}/2 with no anomalous part (number-conserving), or vanishing restricted Majorana covariance <i gamma_{1a} gamma_{1b}> = 0 for a != b (Majorana), and the same on site 2; a Gaussian state with vanishing restricted covariance is the maximally mixed state, so rho_1 = rho_2 = 1/d for both states. The reduced state in the spin (Jordan-Wigner) sense coincides with the fermionic one for definite-parity states: even operators on a site agree in the two languages, and odd ones have zero expectation in both. So Delta = |Psi_+><Psi_+| - |Psi_-><Psi_-| is admissible in Lemma 2.1 with vanishing marginals, and Lemma 2.1 gives M >= tr(Delta h_{12})/2 = tr(Delta h_hop)/2 = ||h_hop||_op, because on-site terms of h_{12} are of the form X x 1 and 1 x X and have zero expectation against a Delta with vanishing marginals. This is ||B||_1 or ||B||_1/2 by Lemma 2.2. When A_0 = 0 the bare bond is h_hop itself, whose half-spread is the same number, so M = ||h_hop||_op exactly. For reference, with the on-site term split evenly between the two sites, M <= ||B||_1 + ||A_0||_1/2 in the number-conserving case (half-spread of f^dag A_0 f is ||A_0||_1/2).

### 3. Counting and velocity lemmas

Lemma 3.1 (counting) [proved]. z^n det h(k) = det(B^dag + z A_0 + z^2 B) (number-conserving: replace n by m; Majorana: det(i(-B^T + z A_0 + z^2 B))) is a polynomial in z whose monomials have degree between n - r and n + r, r = rank B, since a term with more than r factors of z^2 B or of B^dag vanishes. If det h is not identically zero (no exactly flat zero band; such a band carries no velocity and is outside the conformal hypothesis) there are at most 2r roots on the unit circle counted with multiplicity, and the multiplicity of a root k_F is at least dim ker h(k_F), because each analytic eigenvalue branch through zero contributes at least one factor of (k - k_F). So the number of crossings, counted as in Section 1, is at most 2 rank B <= 2n.

Lemma 3.2 (velocity) [proved]. For u in ker h(k_F), f(k) = <u|h(k)|u> = alpha + 2 Re(e^{ik} beta) with alpha = <u|A_0|u> (times i in the Majorana case) real and beta = <u|B|u> (times i in the Majorana case). f(k_F) = 0 and f'(k_F) = <u|h'(k_F)|u> = v, so v^2 = 4|beta|^2 - alpha^2 <= 4 |<u|B|u>|^2 <= 4 ||B||_op^2.

Theorem 3 follows: A <= (1/2)(2 rank B)(2 ||B||_op) = 2 rank(B) ||B||_op in the number-conserving case, and for m = 1 this is A <= 2|B| = 2 ||B||_1 <= 2 M. Equality needs alpha = 0 and |sin k_F| = 1 at both Fermi points, i.e. A_0 = 0 (half filling), which is the XX chain with hopping B = t e^{i phi}; the phase is removed by a momentum shift and does not change A or M. The same counting gives A <= rank(B) ||B||_op <= 2 rank(B) M in the Majorana language, which is weaker than Theorem 2 as soon as rank B >= 2; Theorem 2 is the statement that does not degrade with the number of bands.

### 4. The Drude identity and Theorem 2

Let P(k) be the spectral projector of h(k) onto its negative eigenvalues and E(k) = tr(P(k) h(k)) = sum of the negative eigenvalues. Away from Fermi points 0 is not in the spectrum, P is analytic in k (resolvent integral over a contour enclosing the negative spectrum), and so is E; band crossings between two negative or two positive bands do not matter.

Step 1 (derivatives). E' = tr(P h'), because tr(P' h) = 0: from P^2 = P one has P' = P P' + P' P hence P P' P = 0 and (1-P) P' (1-P) = 0, while h = P h P + (1-P) h (1-P). Then E'' = tr(P h'') + tr(P' h'). With P' = (1/2 pi i) oint R h' R dz, R = (z - h)^{-1}, and eigenpairs (eps_j, u_j),

    tr(P' h') = sum over j occupied, l unoccupied of 2 |<u_l|h'|u_j>|^2 / (eps_j - eps_l) <= 0,

since eps_j < 0 < eps_l (residues of 1/((z - eps_j)(z - eps_l)) at both poles cancel when both are inside or both outside the contour).

Step 2 (jumps). E' = tr(P h') is periodic and piecewise smooth. At a Fermi point k_F, P gains or loses the kernel projector P_F (a branch with slope nu < 0 enters the occupied set, one with nu > 0 leaves it), so E' jumps by - sum over the kernel branches of |nu|. Periodicity gives

    integral over the smooth part of E'' dk = sum over all crossings of |v|   (with multiplicity, exactly as in A).

Step 3 (dropping the regular part). Since tr(P' h') <= 0,

    sum over crossings of |v| <= integral over the circle of tr(P(k) h''(k)) dk.

Step 4 (the hopping energy). Number-conserving: h'' = -(e^{ik} B + e^{-ik} B^dag) = -(h - A_0) is minus the Bloch symbol of h_hop, and integral dk/(2 pi) tr(P(k)(e^{ik} B + e^{-ik} B^dag)) = <f_x^dag B f_{x+1} + h.c.>_GS = <h_hop>_GS. So sum |v_F| <= 2 pi (-<h_hop>_GS) <= 2 pi ||h_hop||_op = 2 pi ||B||_1, hence A = (1/2) sum |v_F| <= pi ||B||_1 <= pi M by Theorem 1. Majorana: h'' = -i(e^{ik} B - e^{-ik} B^T) = -(h - i A_0) is minus the symbol of the hopping part, and the ground-state hopping energy per bond is <h_hop>_GS = (1/2) integral dk/(2 pi) tr(P(k) i(e^{ik} B - e^{-ik} B^T)) (the factor 1/2 is the particle-hole redundancy; checked on the single Majorana chain H = i t sum gamma_x gamma_{x+1}, where e_0 = -2t/pi). So sum |v| <= 4 pi (-<h_hop>_GS) <= 4 pi ||h_hop||_op = 2 pi ||B||_1 and A = (1/4) sum |v| <= (pi/2) ||B||_1 = pi (||B||_1/2) <= pi M. This proves Theorem 2.

Where the constant comes from. Step 3 is an equality iff h' never couples occupied to unoccupied bands, which holds when A_0 and B are simultaneously diagonalizable (decoupled bands). Then -<h_hop>_GS = sum_s (2/pi) |b_s| sin k_{F,s} <= (2/pi) ||B||_1 and the chain of inequalities gives A <= 2 ||B||_1, the value of Theorem 4, attained at half filling. For coupled bands -<h_hop>_GS can reach ||B||_1 (the dimerized rank-one chain with A_0 = 0 has flat bands at +-sigma and hopping energy -sigma per bond), but then Step 3 is maximally lossy (A = 0). So pi is the constant that the two inequalities give when each is saturated separately; and Theorem 4 shows they are never saturated together.

Remark (the first inequality is an exact Drude weight). The quantity sum |v_F| / (4 pi) is Kohn's Drude weight D = (L/2) d^2 E_0/d phi^2 at zero twist, and the chain of inequalities is the standard f-sum statement D <= -<h_hop>_GS / 2 (the regular part of the conductivity is nonnegative). The novelty here is only the pairing of this with the gauge-invariant bond scale through Theorem 1, and its extension to the Majorana form.

### 5. The sharp constant

Theorem 3 settles m = 1 without pairing by counting. The general case is Theorem 4, proved at the end of this section after the evidence that led to it.

What the numerics say [numerical]. (i) Random sampling of 2000 chains per size (complex m = 1..4, Majorana n = 1..6) never exceeds A/||B||_1 = 2 (complex) or A/(||B||_1/2) = 2 (Majorana). (ii) Direct maximization of these ratios by Nelder-Mead from 30 to 60 random starts, for complex m = 2, 3 and Majorana n = 3, 4, 6, with every candidate recomputed on a fine k-grid, ends at 2.00000 for all cases that converged, and at the optimum A_0 vanishes and B is normal to 10^-5 or better (files `best_*.npz`, logs `log_search_*.txt`). The n = 6 Majorana search did not converge within the budget (best 1.655) but the explicit stacked chains give 2 there too. (iii) Since M >= ||B||_1 (or ||B||_1/2) with equality at A_0 = 0, and the optima have A_0 = 0, the same numbers hold for A/M.

Why the natural lemmas fail [proved by counterexample]. One would like to prove sum |v_F| <= 4 ||B||_1 from Lemma 3.2 by bounding the frame operator S = sum_F |u_F><u_F| of the Fermi vectors by 2. This is false: for B = diag(J_1, J_2), A_0 = t sigma_x with t^2 < 4 J_1 J_2 there are four Fermi points whose kernel vectors are (1, -+ sqrt(J_1/J_2)) (two each), and S has largest eigenvalue 4 max(J_1, J_2)/(J_1 + J_2), up to 4. The Fermi vectors pile up on the weakly hopping orbital, but the velocities there are 4 |sin k_F| J_1 J_2/(J_1 + J_2), and A/||B||_1 = 8 J_1 J_2 |sin k_F|/(J_1 + J_2)^2 <= 2. The velocity weighting is essential; any proof must use it.

Structure available for a proof [proved]. (a) For distinct Fermi points F, G, subtracting <u_G|h(k_F)|u_F> = 0 from <u_G|h(k_G)|u_F> = 0 gives <u_G|B^dag|u_F> = z_F z_G <u_G|B|u_F>, so Gamma = Z U^dag B U (U the matrix of Fermi vectors, Z = diag(e^{i k_F})) is Hermitian off the diagonal with diagonal -alpha_F/2 - i v_F/2, and sum |v_F| = -2 Im tr(C B) with C = sum_F sgn(v_F) e^{i k_F} |u_F><u_F|. (b) When A_0 = 0 the Fermi points come in pairs (z, -z) sharing the kernel vector u_s, distinct u_s are B-orthogonal, B^dag u_s = lambda_s B u_s with |lambda_s| = 1, and the claim becomes sum_s |<u_s|B|u_s>| <= ||B||_1 for this B-orthogonal family; for B > 0 it follows from Kantorovich's inequality. (c) Convex-geometric form, A_0 = 0: with W = {tr(P B) : 0 <= P <= 1}, a planar convex set, E(k) = 2 min over w in W of Re(e^{ik} w), Fermi points are the directions in which W has a flat edge, |v_F| is twice the edge length, and the claim reads: the total edge length of W is at most 2 ||B||_1. For normal B, W is the zonotope sum_s [0, b_s], all of whose perimeter 2 sum |b_s| is edges, which is the equality case. For B = sigma |a><b| with a orthogonal to b, W is a disc of radius sigma/2 with no edges. With A_0 the same holds in R x C for K = {(tr(P A_0), tr(P B))} probed along the directions (1, 2 e^{ik}).

Proof of Theorem 4 [proved]. Work in the number-conserving notation h(k) = A_0 + z B + conj(z) B^dag with B complex; the Majorana case is the same with A_0 -> i A_0 and B -> i B, which does not change ||B||_1. Index the Fermi terms by a: z_a = e^{i k_a}, u_a a unit vector in ker h(k_a) from the velocity-diagonalizing basis, velocity v_a = <u_a|h'(k_a)|u_a>.

Step 1 (lift). Define the 2n x 2n Hermitian matrix J = [[0, iB], [-iB^dag, 0]] and the unit vectors x_a = (u_a, z_a u_a)/sqrt2 in C^{2n}. Then x_a^dag J x_a = (i/2) u_a^dag (z_a B - conj(z_a) B^dag) u_a = v_a/2, and J^2 = diag(B B^dag, B^dag B) gives tr|J| = 2 ||B||_1.

Step 2 (J-orthogonality). For z_a != z_b, relation (a) above, u_b^dag B^dag u_a = z_a z_b u_b^dag B u_a, gives x_b^dag J x_a = (i/2)(z_a u_b^dag B u_a - conj(z_b) u_b^dag B^dag u_a) = (i/2) z_a (1 - |z_b|^2) u_b^dag B u_a = 0. For z_a = z_b, x_b^dag J x_a = (1/2) u_b^dag h'(k_a) u_a = (v_a/2) delta_ab by the choice of basis. So all Fermi terms are mutually J-orthogonal, with J-norms d_a = v_a/2.

Step 3 (trace lemma). Let J be Hermitian and x_1, ..., x_p Euclidean unit vectors with x_i^dag J x_j = 0 for i != j and d_j = x_j^dag J x_j. Then sum_j |d_j| <= tr|J|. Proof: keep the vectors with d_j > 0, X = [x_j], so X^dag J X = D = diag(d_j) > 0. They are linearly independent and J is positive definite on their span W (w = Xc gives w^dag J w = c^dag D c). Write X = V C with V an orthonormal basis of W and C invertible, H = V^dag J V, G = C^dag C = X^dag X (so G_jj = 1). From C^dag H C = D, tr H = tr(D G^{-1}) = sum_j d_j (G^{-1})_jj, and Cauchy-Schwarz 1 = |(G^{1/2} e_j)^dag (G^{-1/2} e_j)|^2 <= G_jj (G^{-1})_jj gives (G^{-1})_jj >= 1. Hence sum_{d_j > 0} d_j <= tr H = tr(P_W J) <= tr J_+, the last step because tr J_+ - tr(P_W J) = tr((1 - P_W) J_+) + tr(P_W J_-) >= 0. The same for -J and the vectors with d_j < 0 gives sum_{d_j<0} |d_j| <= tr J_-, and adding, sum |d_j| <= tr|J|. Equality forces (G^{-1})_jj = 1 for all j, hence (equality in Cauchy-Schwarz) e_j eigenvectors of G with G = 1, so the vectors of each sign are orthonormal, and W equal to the corresponding spectral subspace of J; together with J-orthogonality this makes the vectors with d_j != 0 an orthonormal eigenbasis of ran J.

Step 4 (conclusion). sum_a |v_a| = 2 sum_a |d_a| <= 2 tr|J| = 4 ||B||_1. With Theorem 1 this is A <= 2 M in both languages.

Step 5 (equality). If sum |v_a| = 4 ||B||_1 and B != 0, the lifted vectors with v_a != 0 are an orthonormal eigenbasis of ran J: J x_a = d_a x_a, whose block equations read i z_a B u_a = d_a u_a and -i B^dag u_a = d_a z_a u_a, i.e. B u_a = b_a u_a and B^dag u_a = conj(b_a) u_a with b_a = -i d_a conj(z_a). Then z_a b_a + conj(z_a b_a) = 0, so the Fermi equation h(k_a) u_a = 0 reduces to A_0 u_a = 0. Since the x_a span ran J = ran B (+) ran B^dag, the u_a span ran B = ran B^dag, on which B and B^dag commute by the eigen-relations, and both vanish on the complement, so B is normal; and A_0 annihilates ran B, so A_0 B = 0. Conversely, if B is normal with A_0 B = 0, h(k) splits into the scalar blocks 2 Re(e^{ik} b_s) on the eigenvectors of B with b_s != 0, plus A_0 restricted to ker B (which must be invertible, or det h would vanish identically); each scalar block has two Fermi points with velocities -+ 2|b_s|, and sum |v| = 4 sum_s |b_s| = 4 ||B||_1. Scalar check: n = 1, A_0 = a, B = b: sum |v| = 2 sqrt(4|b|^2 - a^2) <= 4|b| with equality iff a = 0.

Both identities of Steps 1 and 2 and the lemma of Step 3 are checked numerically on 300 random chains per size in `lift_check.py` (off-diagonal J-overlaps below 10^-13, ratio sum|d|/tr|J| at most 1.00000).

Why the Drude bound is weaker. Theorem 2 drops the regular part of the conductivity and then bounds the ground-state hopping energy by its operator norm; the two steps are never tight together. Theorem 4 instead works directly with the Fermi data, and the lift to C^{2n} is what makes all Fermi terms (including the two at +-z sharing a kernel vector when A_0 = 0) simultaneously orthogonal for one Hermitian form whose trace norm is 2 ||B||_1. In the convex-geometric picture (c) above it says that the total edge length of {tr(P B)} is at most 2 ||B||_1, with the zonotopes of normal B as the only equality cases.

Saturating chains [proved]. A_0 = 0 and B normal. Number-conserving: B = sum_s b_s |s><s| in an orthonormal on-site basis, each band 2 Re(e^{ik} b_s) = 2|b_s| cos(k + arg b_s) at half filling, A = 2 sum |b_s| = 2 ||B||_1 = 2M. These are m decoupled XX chains with complex hoppings, rotated by a common on-site unitary (which preserves both A and M). Majorana form: a real normal B is orthogonally block-diagonal with 2 x 2 blocks sigma R(theta) and 1 x 1 blocks +-sigma; a 2 x 2 block is the XX chain with complex hopping sigma e^{i theta}, a 1 x 1 block is the single Majorana chain i sigma sum gamma_x gamma_{x+1} (two self-paired crossings at k = 0, pi with |v| = 2 sigma, A = sigma, M = sigma/2). Both have A/M = 2. Adding any on-site term A_0 that commutes with B and acts on ran B lowers A (velocities sqrt(4|b_s|^2 - mu_s^2)) without lowering M; an A_0 acting only on ker B adds gapped on-site modes, keeps sum |v| = 4 ||B||_1, but raises M by half its spread, so A/M drops below 2.

### 6. Weakly interacting perturbations

XXZ line [exact, recalled; M proved]. The spinless chain -t sum (c^dag c + h.c.) + V sum (n_x - 1/2)(n_{x+1} - 1/2) is the XXZ chain with J = 2t, Delta = V/(2t), at half filling. The bond spectrum is {V/4, V/4, -t - V/4, t - V/4}, so the bare half-spread is t for -2t <= V <= 2t, and the dual states (|10> +- |01>)/sqrt2 (equal marginals 1/2) give M >= t; hence M = t for |Delta| <= 1. The Bethe-ansatz velocity is v = (pi J/2) sin(gamma)/gamma, Delta = cos gamma, c = 1, so A/M = pi sin(gamma)/gamma: 2 at the free point, pi at Delta = 1, and smaller than 2 for attraction. `interactions.py` tabulates it.

First order [proved]. Hartree-Fock for the spinless NN interaction at filling k_F: the Hartree term is k-independent and is absorbed in the chemical potential at fixed density; the Fock term renormalizes the hopping to t_eff = t + V <c_x^dag c_{x+1}> = t + V sin(k_F)/pi, so

    v = 2 t sin k_F + (V/pi)(1 - cos 2 k_F) + O(V^2),   delta A = (V/pi)(1 - cos 2k_F),

which at half filling is 2V/pi and matches the expansion of the Bethe result. M is unchanged to first order (the extremal bond states have <(n_1 - 1/2)(n_2 - 1/2)> = -1/4 on both one-particle states, a common shift that the alpha gauge removes; exactly, M = t). So A/M = 2 + 2V/(pi t) + O(V^2): repulsive interactions push the ratio above the free value 2. Extrapolating the first-order line to Delta = 1 gives 3.27 > pi, while the exact value is pi; first-order perturbation theory overestimates the velocity at finite coupling and is not a reliable tool for the question of saturation.

General multi-band first order [proved, elementary]. For a two-site interaction V_int added to a quasi-free chain with simple Fermi points, the first-order velocity shifts are the slopes at k_F of the Hartree-Fock self-energy, delta v_F = d/dk <u_F| Sigma_HF(k) |u_F>, and delta A is their sum with weight 1/2 (or 1/4). M changes at first order by at most ||V_int||_op, and by (<psi_max|V_int|psi_max> - <psi_min|V_int|psi_min>)/2 when the gauge-optimized bond has nondegenerate extremal eigenvectors (an envelope argument gives this as an upper bound on delta M even when the optimal gauge moves). Both shifts are O(V), so A/M moves continuously from its free value. Since every free chain has A/M <= 2 (Theorem 4), no first-order perturbation of a free chain with m >= 2 can approach the conjectured pi m, and for m = 1 the only way to reach pi is at finite coupling, as XXZ does at the isotropic point. The Hubbard chain (m = 2, d = 4) is a case in point [recalled]: away from half filling its charge and spin velocities shift by equal and opposite amounts at first order in U, so delta A = 0, while M >= 2t; at half filling the charge sector gaps and A decreases.

## Numerical checks

Python 3 with numpy and scipy only. All scripts run from this folder.

- `ff_bound.py`: Bloch matrices, Fermi points (unimodular eigenvalues of the companion pencil of B^dag + z A_0 + z^2 B), velocities (degenerate Hellmann-Feynman), A, the many-body two-site bond via Jordan-Wigner, the complex-to-Majorana map, and the first random scans. Running it reproduces the XX and Ising sanity checks and the random-scan maxima (all 2.00 or below). About 2 minutes.
- `gauge_M.py`: gauge-optimized M by an annealed log-sum-exp smoothing with Hellmann-Feynman gradients, polished by Nelder-Mead; reference values XX (M = 1), XX with chemical potential (M = 1), Ising (M = 0.7071, above the lower bound 0.5), two XX chains (M = 2), SSH (M = 1.4142, above the lower bound 1). Seconds.
- `fsum_check.py`: the two inequalities of Section 4 on 300 random chains per size, with the hopping energy integrated piecewise between Fermi momenta by Gauss quadrature. Output: max of sum|v| / (2 pi (-<h_hop>)) = 1.00000 (complex) and of sum|v| / (4 pi (-<h_hop>)) = 1.00000 (Majorana), and max of sum|v| / (2 pi ||B||_1) = 0.6366 = 2/pi in both languages. About 2 minutes.
- `search_ratio.py kind n trials seed` with `params.py`, `verify_optimum.py`: adversarial maximization of A/||B||_1 (complex, `nc`) or A/(||B||_1/2) (Majorana, `mj`); each candidate is recomputed by `verify_optimum.brute_A` (20001-point grid, sign changes, finite-difference slopes) and only brute-force-confirmed improvements are recorded. Logs in `log_search_*.txt`, optima in `best_*.npz`. Results: 2.00000 for nc n = 2, 3 (two seeds) and mj n = 3, 4; 1.655 for mj n = 6 (not converged). Ten to forty minutes per job.
- `table.py`: A, the proved lower bound on M, the bare half-spread, the gauge-optimized M and A/M for eleven reference chains (XX at two fillings, Ising, critical Kitaev, Majorana stacks, two XX chains with and without on-site mixing, SSH, complex hopping, a random m = 2 chain). The random chain has M = 4.151 strictly between the lower bound 4.128 and the bare value 4.164, so the gauge optimization is active there. All ratios are at most 2.
- `interactions.py`: the XXZ table of Section 6 and the Hartree-Fock check at three fillings.
- `lift_check.py`: the lifted vectors x_a and the form J of Theorem 4 on 300 random chains per size (n <= 4, both languages): J-orthogonality of distinct Fermi terms to 10^-13 and sum |d_a| <= tr|J|. Under a minute.

Tolerances: Fermi points are accepted when the companion eigenvalue has modulus within 10^-6 of 1 and the kernel eigenvalue is below 10^-6 of the spectral radius; the brute-force recomputation protects against double counting of near-unimodular pairs, which was the source of two spurious ratios (2.15 and 3.0) in a first run that is not kept.

## What remains

1. The gap between 2 (free) and pi (free bound, and Heisenberg) at d = 2 is where interactions act. A quantitative version, A/M as a functional of the Luttinger parameter along other integrable deformations (XYZ, dimerized XXZ), would say whether pi M log_2 d can be exceeded at d = 2 by any interaction; XXZ never exceeds it.
2. For m >= 2 the free bound is 2 M, far below pi m M. Whether interacting chains with m >= 2 (Hubbard away from half filling, SU(4) at d = 4 with A = 3 pi M/2) can approach pi m M is the real content of the conjecture and is untouched here.
3. Theorems 2 and 4 extend verbatim to longer-range quasi-free chains with the bond redefined, and to period-p chains with the inter-cell hopping matrix in place of B; the bound then compares A to the inter-cell scale, which is the natural one for the blocked chain. Not written out.

## Relation to the rest of the repository

This addresses the velocity-weighted bound A <= pi M log_2 d of `specs/_bound_background.md` on the quasi-free class, where C1 proved c <= log_2 d. It sharpens the remark in `contrib/refined_conjecture/` that A <= pi ||h|| for number-conserving free fermions: the constant is now the sharp 2, the bound uses the gauge-optimized M, covers pairing, and comes with the exact saturation structure (free chains reach 2 M only for XX-type stacks, the Heisenberg point reaches pi M). It does not change "Where things stand": the conjecture holds in this class with room log_2 d, and the first example that reaches the conjectured value, at d = 2, is interacting.

## Appendix: the Codex run

`spec_sharp_constant.md` is the self-contained statement given to gpt-6-astra (it contains the known facts F1 to F5 and the reformulations R1 to R4 as hints); `codex_sharp_constant.out` is its answer, which is the proof reproduced in Section 5 with the equality characterization B normal, A_0 B = 0. My check: the J-orthogonality (Step 2) is the relation (a) I had derived independently; the lemma (Step 3) is elementary and I re-derived it; the equality analysis I verified against the explicit saturating family and against the one-band formula. `lift_check.py` confirms the two identities numerically.
