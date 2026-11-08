This report was written by a Claude Opus agent (2026-11-08), not by Codex. Astra was not used (out of credits).

# C161: referee report on C160 (c <= 1 near the ferromagnetic Heisenberg point)

Tags: [proved], [derived], [numerical], [heuristic], [assumed], as in C160. Scripts and console captures are in this folder (twisted_hF.py/.out, soliton_valley.py/.out, soliton_pair.py/.out, csl_phase.py/.out, bloch_hessian.py/.out). Notation and equation numbers are those of codex/work_C160/PROOF.md. All local runs took seconds to about a minute. No Typhon job was needed.

## Summary of grades

| Step | Grade |
|---|---|
| 1. Blow-up and hbar_eff | sound with repair (the singular set is a curve, and hbar_eff is not the IR semiclassical parameter at onsets) |
| 2. Semiclassical correspondence (SC) | fails as stated (not uniform). It is repairable at magnon onsets, and open at soliton onsets |
| 3. M1 (Mather graph plus Poincare-Bendixson) | sound (continuum), heuristic on the lattice as C160 says |
| 4. M2 (no Bloch zero mode at k != 0) | sound, with a rephrasing of the scalar reduction |
| 5. M3 (pairing with the phason) | sound for dense textures. The C160 numerics claim about "floating + SSB" is family-specific |
| 6. Lemma 3' (magnon minima at 0 or pi) | sound for magnons. The a_s != 0 case needs one more line. It does not extend to solitons |
| 7. Two-magnon, McGuire collapse, lattice step | sound for the U(1) part, heuristic otherwise |
| Overall | C160 controls the corner eps -> 0 of h_F only. It does not remove the C159 Proposition C window |

## 1. The blow-up (PROOF section 1, 2)

1a. [proved] The gradient expansion (1.2) and the gauge form (1.3) are correct. The pure-DM check is correct too. J = -1 + eps A(D) gives J J^T = 1 + eps^2 (|D|^2 - D D^T), so Delta_AF = -1/sqrt(1 + eps^2 |D|^2).

1b. [derived] hbar_eff is uniformly small as an upper bound. Every term of the bond is O(eps) per site while the stiffness at h_F is 1. So no classical texture varies on the lattice scale, and the coherence length obeys l >= c eps^{-1/2} in every direction u. Hence 1/(s l) <= C eps^{1/2} uniformly. The directions with only a field (polarized, gapped), only a traceless anisotropy (zero-field XYZ near -1: Theta-broken pair, or XXZ easy-plane c = 1, matching Baxter) and the pure-DM wedge all fall under this bound. In none of them does the coherence length stay O(1).

1c. [proved, numerical check twisted_hF.out] Repair: the singular set through h_F is a curve, not a point. For any rotation R(phi) about an axis m, J = -R(phi) (B = 0) is h_F after sigma_j -> R(j phi) sigma_j. Open chains at L = 8 have the h_F spectrum to 5e-14 at phi = 0.3, 1, 2, pi. Modulo common rotations this is a curve gamma = { -R(phi), 0 <= phi <= pi } in B_2. Its far end is -R(pi) = XX + YY - ZZ in Pauli form, the staggered image. The curve enters every punctured ball around h_F. Expanding, -R(phi) = -1 - sin(phi) A(m) + (1 - cos phi)(1 - m m^T), so it lies in C160's pure-DM wedge W_Z at the point where W in (2.2) cancels the DM-induced easy plane exactly. There the rescaled functional is gauge-equivalent to |m'|^2/2, scale-free, with a degenerate classical ground state. So "punctured neighbourhood of h_F" must read "U_eps0 minus gamma". The classical part of C160 (section 3) still applies near gamma, because W in (2.2) is arbitrary. The correct semiclassical parameter there is again 1/(s l) with a diverging l. This repair costs nothing inside U_eps0. Its consequence for the overall verdict is in section 8.

1d. [derived] hbar_eff = 1/(s l) is the semiclassical parameter of the order-parameter field. It is not the semiclassical parameter of the infrared boson near onsets. At the edge of the cone in W_C (theta_0 -> 0) the conjugate density is s(1 - cos theta_0), approximately s theta_0^2/2, which goes to zero. The magnons form a dilute Lieb-Liniger gas with gamma_LL of order eps/density, which diverges, so K -> 1 (Tonks). The same holds at every Pokrovsky-Talapov edge of a soliton lattice. There the solitons are dilute, their exponential interactions vanish, and K -> 1. So the claim "hbar_eff uniformly small on U_eps0" is true for the classical field but false for the low-energy boson in boundary layers next to every z = 2 surface. That is where step 2 fails.

## 2. The semiclassical correspondence (SC)

2a. [derived] Fails as stated. SC asserts that the quantum infrared follows from the classical soft-mode count uniformly as hbar_eff -> 0. By 1d this cannot hold in the dilute layers next to onsets, which are open sets of every U_eps0. There the species count is fixed by the quasiparticle band at the onset, and that is a quantum statement. For magnons C160 supplies it (Lemma 3' plus the two-magnon problem, sections 4 and 5), so the magnon onsets are repaired. For solitons C160 supplies nothing. The soliton analogue of Lemma 3' is false, as shown by the following explicit bond.

2b. Explicit counterexample to "one soliton species at onset" [numerical, soliton_valley.out]. Take the K-stratum bond

    v = (-1, -1 + mu eps^2/2, -1 + kappa eps^2/2, 0, 0, -eps, -h eps^2, 0, 0)   (est.py coordinates),

with kappa = 0.05 and mu = 2.3, small eps and h near 0.0265. This is C160's CSL family with an in-plane anisotropy mu added. Its blow-up is E = int [ |n'|^2/2 - (n x n')_z - h n_x + (kappa/2) n_z^2 + (mu/2) n_y^2 ]. The Z2 n_z -> -n_z is Theta R_pi(z), exact. The background +x is the global minimum of W, since h + kappa > 0 and h + mu > 0, and it is Z2-symmetric. Results:
- The planar 2 pi soliton reaches zero energy at h = 0.024 to 0.03. There its out-of-plane Hessian has lowest eigenvalue -0.067 to -0.068. This is consistent with -0.117 at kappa = 0 from the onset scan, since kappa shifts the eigenvalue by +kappa. So the soliton is tilt-unstable over a stable background.
- Full minimization gives two degenerate tilted solitons (max |n_z| = 0.55), with energy 0.02 below the planar one at h = 0.02 to 0.06. Their oriented solid angles are 2 pi (1 -+ 0.443) mod 4 pi.
- [derived, Haldane PRL 57, 1488 (1986)] A ferromagnetic texture carries crystal momentum P = s Omega mod 2 pi. With s = 1/2 a planar 2 pi soliton has P = pi (this is the s = 1/2 Berry phase at the soliton scale), and the two tilted solitons sit at P = pi -+ 0.443 pi. In the collective coordinates (X, tilt) the Berry term is s Omega(tilt) dX/dt. So the tilt is the soliton momentum, and E_sol(tilt) is the soliton dispersion E(P). A double well in the tilt is a soliton band with two degenerate minima at incommensurate momenta, not related by any DM shift.

This is the J1-J2 two-minimum mechanism, the one C160 Task 2 was asked to exclude, but now for solitons rather than magnons. It is not a Hessian zero mode at k != 0, so M2 does not see it. M3 sees it only as a tilt conjugate to the phason in the dense regime.

2c. Does it open a c > 1 window? [numerical, soliton_pair.out, csl_phase.out] In this example it does not. Like tilts attract and opposite tilts repel, with |U| proportional to exp(-T sqrt(h + kappa)). The decay rate fitted from T = 18 and 22 is 0.27, against sqrt(0.0765) = 0.28. So E_sol + U_same(T = 10) < 0 while E_sol > 0, and the onset is first order. The dilute two-valley gas is preempted, and the uniform state goes directly into a same-tilt (conical-like) CSL, <n_z> = 0.44. That CSL is the classical ground state for 0 <= h <= 0.0265, below both the planar CSL and the uniform state. Its infrared is one phason plus a spontaneously broken uniform Theta R_pi(z), c = 1. [heuristic] The general reason is as follows. The Z2 forces D . n_bg = 0, so the out-of-plane tails decouple from the in-plane ones and are monotone. Monotone same-sign tails of a massive scalar attract, so equal valleys attract more than unequal ones (g_intra < g_inter). This is the immiscible case of the dilute two-component Bose gas (Kolezhuk et al., PRB 85, 064420 (2012)), which gives one valley and Z2 SSB, not c = 2. A c = 2 miscible soliton gas would need g_intra > g_inter at the onset. That needs a valley-independent repulsion that decays more slowly than the tilt tails, together with a tilt-unstable core, and the two pull against each other (slow in-plane tails want mu close to -h, the tilt instability wants mu - kappa large). I did not find such a bond. I also did not exclude it.

2d. The other quantum effects named in the spec [derived or heuristic]:
- Order by disorder acts within O(hbar_eff) of a classical degeneracy and selects among classical states. It cannot create a second canonical pair. C160's treatment is sound.
- The s = 1/2 Berry phase at the soliton scale is not a topological term (a ferromagnet has none). It does give every planar 2 pi soliton momentum pi, as in 2b. That shifts the momentum assignment of phason vertex operators in a CSL by pi per soliton. C160 never uses it, but C159's commensurability bookkeeping (Lemma A, "half-floating") assumed U(1)-parent momenta. Whether the extra pi per soliton changes which harmonics are commensurate in a CSL far from the U(1) class is unchecked [heuristic remark, not a finding against C160's count].
- Quantum melting of a soliton lattice: away from onsets K is of order 1/hbar_eff and the lattice is a single boson. At onsets see 2a to 2c.
- Soliton bound states: the like-tilt attraction of 2c binds soliton pairs classically. That is the mechanism of the first-order onset, not a second fluid.

## 3. Lemma M1

[proved, given Mather] L = |v|^2/2 + A(n).v + W(n) with A = D x n is Tonelli on TS^2. A is a global 1-form, so the magnetic 2-form is exact, and H^1(S^2) = 0 makes the minimizing class unique. Mather's graph theorem applies to time-independent Tonelli Lagrangians with magnetic terms. The relevant set is the Mather set, the supports of minimizing measures, and that is what the bulk ground state needs, since its mean energy is the minimal one. Aubry and Mane orbits add only isolated solitons or domain walls at zero density. Poincare-Bendixson holds for any continuous flow on S^2, so recurrent points are fixed or periodic. C160's Lipschitz extension is harmless because the set K is invariant. Coexisting ergodic minimizers project to disjoint closed curves and fixed points. They are first-order coexistence, or nongeneric continuous families (for example the twisted Heisenberg points on gamma, where every cone in the rotating frame is minimizing). Grade: sound.

Lattice [heuristic, as C160 says]. The discrete graph theorem needs a twist condition. P_n J P_n' degenerates when n' is perpendicular to n, so it holds only for nearly parallel neighbours, which is the case near h_F. Poincare-Bendixson has no analogue for maps, so the single-phason statement on the lattice rests on persistence of the continuum invariant circle. Near homoclinic limits (PT edges) the circle breaks into cantori, which are pinned and c = 0. This is acceptable.

## 4. Lemma M2

[derived] Sound. The worry in the spec is that DM makes the linearization a two-component gyroscopic problem that is not a scalar Hill equation. The reduction is nonetheless valid where it is needed. The 4 x 4 monodromy M of the linearized Euler-Lagrange flow preserves the energy differential dE (dE . M = dE). So an eigenvector with multiplier mu != 1 satisfies dE(v) = mu dE(v), which forces dE(v) = 0, and it lies in the energy-preserving 3-dimensional subspace. Modulo the flow direction this subspace is a 2-dimensional symplectic space with multipliers lambda, 1/lambda. Hence every Bloch zero mode at k != 0 is an eigenvector of the reduced 2 x 2 linearized Poincare flow, and no off-shell reduction is ever used.

What needs rephrasing is the "scalar Hill equation of the Jacobi-Maupertuis-Randers metric". At the Mather energy c(L) the Randers metric can be degenerate, because c(L) is the Mane critical value. The clean statement uses Lagrangian subspaces. Orbits in the Mather set are semistatic, so they have no conjugate points (Contreras and Iturriaga). Then the image of the vertical line in the reduced flow never returns to vertical, so the rotation number is zero. That excludes elliptic multipliers (rotation number theta/2 pi) and negative ones (half-turn per period). So lambda > 0, and e^{ikT} is a multiplier only at k = 0.

[numerical, bloch_hessian.out] Test on a non-planar texture with DM along the moment, as the spec asked. The tilted CSL of 2c at h = 0.015 (kappa = 0.05, mu = 2.3) has T = 8.85 and D.n up to 0.60, and its two tangent components are gyroscopically coupled. Its Bloch Hessian bands: lambda_0 rises monotonically from 0 at k = 0 to 0.0149 at k = pi/T, and lambda_1 runs from 0.066 to 0.049, never near zero. On the planar CSL at kappa = 0.15 the tilt band has its minimum, 0.015, at k = 0. At kappa = 0.25 the second band has its minimum at k = pi/T (0.0348 against 0.0407 at k = 0). This is allowed: the bottom of a higher band can sit at k != 0, and only zeros are forbidden. So the answer to the spec's question is this [derived plus numerical]. The bottom of the lowest band is at k = 0 with value 0 (the phason). A higher band can have its minimum at nonzero Bloch momentum, but it cannot reach zero there at a global minimizer.

Limit of M2. It constrains classical Hessian zero modes only. The two-valley soliton band of 2b has its minima at P = pi -+ P0 without any Hessian zero mode at k != 0, so M2 says nothing about it.

## 5. Lemma M3

[derived] Sound for dense textures. Equation (3.1) is correct: Omega(n_0', psi) = s int g |n_0'| dx, and a periodic solution of a disconjugate equation has no zeros. In the parabolic case the reduced monodromy must be a Jordan block, because the identity case would make every solution periodic, including one with a zero. Case (b) is also correct. Along the family of orbits e'(T) = -(E + e)/T and e''(T*) = -E'(T*)/T*, so at T'(E) = 0 the mean energy is cubic in E - E_0 and the orbit is not a minimizer (given T'' != 0). One wording issue: at T'(E) -> 0^- one has e'' -> +infinity, so the soft point is e'' = 0 (r kappa = h^2), as C160's last sentence says, not the point T'(E) = 0. The "alternative is not a minimizer" claim stands.

The numerics claim does not stand as a general statement. C160 RESULT Task 4 says the "floating + SSB" phase does not exist. That holds in the one-parameter family checked (mu = 0). With in-plane anisotropy the tilted CSL is the classical ground state on an open set of the K stratum near h_F (csl_phase.out). This does not contradict M3. The tilt is the nodeless partner of the phason, so the planar-to-tilted CSL transition is z = 2 (classical) or first order. Quantum mechanically this is the coupling sigma d_x theta (Ising order field times the phason's conjugate, dimension 9/8, relevant). C159 counts it as one extra tuning for "Ising x boson with uniform Theta-odd order", so the c = 3/2 point there has codimension 2. The period-2 version (alternating tilts, P1 proper) is a 2T-periodic orbit. M2 forbids its continuous onset from a T-periodic minimizer, and in the example opposite tilts repel. So P1 is not reopened. C160 should still drop the sentence "the floating + SSB phase does not exist" or restrict it to mu = 0.

## 6. Lemma 3'

[proved] Sound. The BdG form (4.1) is general for a one-site cell and a nearest-neighbour bond. Pair terms come only from a_j a_{j+1}, which gives cos k after symmetrization. For s = 1/2 there is no single-ion anisotropy, so no on-site pair term. Linear terms vanish at a stationary reference. The proof for a_s = 0 is correct.

For a_s != 0 (DM along the moment) C160 cites C154 without proof. Here is the one-line version. Zeros of omega = a_s sin k + sqrt(F(cos k)) satisfy G(x) := F(x) - a_s^2 (1 - x^2) = 0 with a_s sin k <= 0. Stability means G >= 0 on [-1, 1], and G is quadratic. So the soft set is a double root (one k, the sign of sin k fixed), or x = +-1 (k = 0 or pi, both only at codimension 2), or G identically zero. The last case requires a_c = 0, A0 = |a_s| and 2|g| = |a_s|, a flat half-band. Near h_F, a_c is about -4, so it cannot occur. So one soft species at every magnon onset near h_F [proved].

Lemma 3' concerns magnons around uniform states. Its soliton analogue is false (2b).

## 7. Two magnons, McGuire collapse, continuum-to-lattice

[proved for the U(1) part] The two-magnon bound state (5.2) is correct (at h_F it gives E_b = 2 - 2 cos K). The string energies per magnon J sinh(gamma) tanh(N gamma/2)/N decrease in N, which gives first-order reversal on the easy-axis side. [heuristic] "Pair creation sources the pair channel, so there is no independent order parameter" is a plausibility argument. It is adequate at the physics standard because the U(1)-breaking scale eps exceeds the binding eps^2.

[heuristic] Lattice step. On the lattice the map is not autonomous, so the multiplier structure {1, 1, lambda, 1/lambda} holds only up to discreteness corrections. These are exponentially small in l (adiabatic invariance), except near lambda = 1 where they can matter. Near PT edges hyperbolic circles break into cantori. Both effects give pinning (c = 0), not extra modes. Acceptable as stated.

## 8. Overall verdict

C160's classical analysis (M1, M2, M3, Lemma 3', the two-magnon problem) is correct at the physics standard, and I found no wrong step in it. M2 in particular survives the two-component gyroscopic test. But C160 does not remove the window left by C159 Proposition C, for three reasons.

(i) Scope in K [numerical, twisted_hF.out]. Proposition C's window is K > 15/8 for Ising x boson and K > 1 for two bosons. For XXZ, K = 15/8 is at Delta_AF = -0.669. On C160's pure-DM line that is eps = 1.11, hbar_eff = 2.2. K = 1 is the XX point, nowhere near h_F. C160's theorem holds for 0 < eps < eps0 with eps0 unquantified and (SC) valid only for hbar_eff << 1. So it covers only the K -> infinity end of the window, not the band 15/8 < K < K(eps0) where the P1 count of C159 is live.

(ii) Scope in the twist angle [proved for the curve, heuristic for its neighbourhood]. The "DM-twisted neighbourhood of h_F" is the neighbourhood of the whole curve gamma = {-R(phi)}. C160 covers only phi of order eps0. Near -R(phi) with phi of order one, the spiral-frame parent is XXZ near Delta = -1, and the non-axial perturbations carry lattice momenta phi and 2 phi. Those vary on the lattice scale, so the continuum functional (1.2) does not describe them. Heuristically the IR there is still one boson: charge-m operators carry momentum m phi, so an emergent U(1) appears for irrational phi/2 pi, with Z_m clock pinning at phi = 2 pi/m. At phi = pi the m = 2 harmonic becomes uniform, and the problem reduces to C160's blow-up with staggered fields entering at second order. But none of this is in C160.

(iii) The onset boundary layers [derived, numerical]. (SC) fails where K -> 1. At soliton onsets the soliton band can have two degenerate minima (explicit K-stratum bond in 2b, P = pi -+ 0.443 pi). In that example the like-valley attraction makes the onset first order into a tilted CSL with c = 1. A miscible (c = 2) case is not excluded by any argument in C160. My heuristic (monotone tails attract when D . n_bg = 0) points against it.

So the honest status of C160 is this. It proves (derived, given SC away from onsets) that no c > 1 phase lives in U_eps0 minus gamma outside the onset layers. At the magnon onsets it is repaired by Lemma 3'. The soliton onsets need a new lemma ("at a soliton onset over a Z2-symmetric background, valleys are immiscible"), which I support only heuristically. No explicit c > 1 bond was found.

Concrete repairs for C160:
- (a) Replace "punctured neighbourhood" by "U_eps0 minus gamma".
- (b) Restrict (SC) to hbar_IR = 1/K << 1, and add an onset lemma for magnons (already present) and for solitons (missing).
- (c) Rephrase M2 via Lagrangian Green bundles.
- (d) Drop or restrict the sentence that "floating + SSB does not exist".
- (e) State plainly that the theorem covers the K -> infinity corner of the C159 window and the phi -> 0 end of gamma only.
