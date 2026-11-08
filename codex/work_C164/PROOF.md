# C164 derivations and numerics: searching for accidental two-valley onsets

This run was done by a Claude Opus agent (2026-11-08), not by Codex. Astra was not used (out of credits).

Tags: [proved], [derived] (physics standard), [numerical], [heuristic], [assumed], [literature]. A tag applies to the paragraph or displayed claim it opens. Conventions are those of est.py (v = (Jx, Jy, Jz, Dx, Dy, Dz, fx, fy, fz), Pauli matrices, field f . sigma_j per site, ring momentum k = 2 pi m / L) and of C160 to C162. Scripts and console captures are in this folder.

## 1. What a counterexample has to look like

[derived, from C162 section 1] A c = 2 phase from this route needs, over a uniform gapped background, a one-particle band with two isolated minima k_A, k_B, both soft at the same point (one tuning for the onset, one for the degeneracy), a dilute gas past onset that does not bind, and a pseudospin chain (1.3 of C162) that is not an Ising ferromagnet, J_z > -|J_perp|.

Two kinematic conditions are needed for the two valley numbers to be emergent conserved charges [derived]. The vacuum can emit a pair of particles in valleys (v, v') only with total momentum zero, so 2 k_A, 2 k_B and k_A + k_B must avoid 0 mod 2 pi, and k_A, k_B must avoid 0 and pi (a single valley at 0 or pi with pair creation allowed gives an Ising-type onset, c = 1/2, not a gas). Otherwise the charges are conserved up to oscillating, irrelevant terms at generic density, as in C162 point (3).

A remark on K [heuristic]. The C161/C162 valleys sit at P = pi -+ delta, so P_+ + P_- = 2 pi = 0 mod 2 pi. A (+, -) pair is therefore emitted from the vacuum without any momentum cost. In the blow-up this is suppressed by the topology of the 2 pi winding and is exponentially small, but at eps of order one only N_+ - N_- is protected. This does not weaken C162's no-go (an immiscible gas stays immiscible), but it does mean that a K-type two-valley gas at eps of order one would not even be a two-charge system.

## 2. Symmetry and the T class

[literature, proved] Lieb-Schultz-Mattis with time reversal (Watanabe, Po, Vishwanath and Zaletel, PNAS 112, 14551 (2015); Ogata, Tachikawa and Tasaki, CMP 385, 79 (2021); cited from memory, not rechecked in this run): a translation-invariant spin-1/2 chain with time reversal Theta (Theta^2 = -1 per site) has no unique gapped symmetric ground state. So in T (B = 0) every uniform gapped background breaks Theta, and its one-particle band has no Theta protection (Theta maps it to the band over the reversed background). Two-valley onsets in T are accidental, as in N. [numerical, local/q_uT.jsonl] In the uT ring screen, 24 of 44 apparently gapped k = 0 ground states at L = 12 had exactly symmetric valleys +-k0. Of the 12 rechecked at L = 16, eleven moved their ground state to k = pi. The twelfth kept k = 0, but its gap fell from 0.241 to 0.146 between L = 12 and 16, faster than 1/L. None is a gapped uniform background, consistent with LSM.

[proved, algebra; numerical kpik_check.out] The linear-spin-wave band of C160 (4.1), omega(k) = a_s sin k + sqrt((A0 + a_c cos k)^2 - pi_0^2 cos^2 k), is invariant under k -> pi - k when a_c = 0. This is an exact symmetry of the full quantum chain in one case. For the Pi-symmetric bond a(XX - YY) + s(XY + YX) + D(XY - YX) + Jz ZZ + f Z, the unitary U = I R_{pi/2}(z) prod_{j odd} R_pi(z)_j maps H to itself. It sends momentum k to pi - k in the odd sector of Pi = prod sigma^z (the one-particle states) and k to -k in the even sector. The staggered rotation flips the transverse block and shifts the momentum by pi on Pi-odd states. R_{pi/2}(z) restores the traceless symmetric transverse part and keeps the DM. Inversion restores the DM and reverses k. Ring spectra at L = 12 agree to 7e-14 in both parity sectors, and a field along x breaks the symmetry (deviations of order one). This is a protected two-valley band (k*, pi - k*) at incommensurate momenta, a stratum C162 did not list because it uses a staggered (non-on-site) operation. It lies inside C152's parity class on its cos(phi) = 0 surface (pure imaginary hopping w = 2iD). There the two valleys are the two Majorana combinations A, B of C152 result 3, which form one Dirac fermion (c = 1) or two Ising lines. So the protected case is already covered by C152, c <= 1 [derived there]. The valley sum is pi, not 0, so pair emission is forbidden by momentum, as it should be for this two-fermion structure.

[numerical, local/q_I.jsonl, q_IhF.jsonl] Inversion-symmetric bonds (D = 0) have E(k) = E(-k) exactly. In 93 + 102 uniform gapped samples at L = 12 no band had isolated incommensurate minima (0 flagged). This is the quantum counterpart of C160 Lemma 3' (a_s = 0, minima at 0 or pi only), which is proved at the LSW level.

## 3. Classical lattice screen (Task 1)

### 3.1 Code and calibration

clat.py is the lattice classical energy E = sum_j [n_j^T J n_{j+1} + f . n_j] with no gradient expansion. It contains the uniform minimum, LSW, periodic-texture stability on rings of length 2 to 24, open chains with clamped background ends, the solid angle Omega with analytic gradient (P = Omega/2 mod 2 pi), constrained relaxation and the tangent Hessian.

- [numerical, lsw_check.out] The classical LSW frequencies equal the exact one-magnon band of U(1) rings (XXZ + D_z + field along z) to 1e-14 for stable references, which also fixes the momentum convention to est.py's.
- [numerical, calib_part1.out] C161 K bond (kappa = 0.05, mu = 2.3) on the lattice, planar soliton at the lattice onset. The lowest tangent-Hessian eigenvalue divided by eps^2 is -0.0544, -0.0645, -0.0671, -0.0677 at eps = 0.2, 0.1, 0.05, 0.025, out-of-plane weight 1.00. C161 continuum: -0.067. Lattice planar onsets h_c = 0.0266, 0.0244, 0.0238, 0.0237 against 0.0238.
- [numerical, calib_part2.out] C162 off-axis K bond (psi = 0.3, h = 0.02886), rings of length 2T. U_same/eps = -0.0082, -0.0027, -0.00041 (eps = 0.1) and -0.0084, -0.0028, -0.00042 (eps = 0.05) at T = 10, 15, 22, against C162's -0.0088, -0.0029, -0.00044. U_alt = +0.0152 and +0.0158 against +0.0169 at T = 10. Like-valley binding (E_sol + U_same < 0) is reproduced, and the approach to the continuum is linear in eps.

### 3.2 Screen

solscreen.py, for each bond: stable uniform background (global minimum on S^2, LSW positive definite, no lower periodic texture on rings of length 2 to 16), then 72 starts of 2 pi loops about six axes perpendicular to the background, two widths, three tilts, both windings, on an open chain of 40 free sites. Each relaxed texture is kept only if it is a genuine local minimum (tangent Hessian >= -1e-7 after polishing), has a large excursion (max |n - n_bg| > 1), and is localized. A two-valley flag needs two such minima with |P_A - P_B| > 0.1. A candidate also needs 0 < E_A < omega_min (the soliton is the lowest excitation) and E_B - E_A < 0.3 E_A. Pair interactions at fixed core separation (one core spin of each soliton clamped, like and unlike, both orders) are coded in pair_interactions().

[numerical, local/cl/] Pilot coverage is 188 bonds: 60 each from the pure-DM wedge on the lattice (J = -1 + eps A(z) + eps^2 s S, f = eps^2 s b, eps in [0.3, 1], s in [0.1, 5]), from near h_F (eps in [0.3, 1.5], N) and from S^8, plus 8 from near h_F in T. 121 had stable uniform backgrounds and 39 had genuine soliton minima, 28 of them below the magnon gap. None had a second soliton valley. A first version without the Hessian test flagged 20 of 51 bonds. Every one of those was a saddle or an unconverged lump, which shows why the test is needed.

[numerical, t_ck.py output in PROGRESS] For the K bond at eps = 0.3 the classical solitons are 18 to 28 sites wide (deviation > 0.3). At eps = 0.5 they are 10 to 12 sites wide, and the tilt instability survives only at larger h (P = pi -+ 0.137 pi at h = 0.15, E = 0.31 against omega_min = 0.35).

[numerical, t_qk.py] The quantum band of the same bonds on a 14-site ring has its minimum at k = 0 (the magnon, gap 0.026 to 0.34). At eps = 0.5, h = 0.15 the band at k = pi -+ 0.14 pi is 0.72, against the classical 0.31. The soliton either does not fit in the ring (eps = 0.3) or is strongly raised by quantum and finite-size effects (eps = 0.5). [heuristic] At eps of order one the classical soliton band is a poor guide to the quantum band, and the quantum screen of section 4 is the decisive one.

## 4. Quantum one-particle bands (Task 2)

### 4.1 Tools

qband.py computes momentum-resolved ring ED (est.py tables, exact to 1e-14), the one-particle band E1(k) (the lowest level above E0 in each sector, the second one at the ground-state momentum) and its trigonometric interpolation. For an isolated band whose particle is small compared with L this interpolation is exact up to aliasing of order e^{-L/(2 xi)}.

qopt.py is the adversarial objective, minimized by Nelder-Mead restarts on S^8 or inside an ensemble. It requires a uniform background (unique ground state at k = 0 at L = 12 and 14) and two valleys with a barrier of at least 5 percent of the bandwidth. The valleys must agree between L = 12 and 14 (positions within 0.3, energies penalized), avoid 0 and pi by 0.35, and satisfy k_A + k_B != 0 mod 2 pi. They must lie below the discrete two-particle floor. The band must be smooth, meaning the weight of the two highest harmonics stays below 2 percent of the bandwidth at both sizes. A gap that closes between L = 12 and 14 like 1/L is penalized as gapless. The objective is |E_B - E_A|/bw + 0.25 E_A/bw plus penalties.

qtune.py continues a two-valley bond toward its onset. It takes minimal-norm Newton steps on the tracked (E_A, E_B) toward a common target gap that shrinks by a factor 0.6 per step, accepting only steps that keep the quality conditions.

### 4.2 Findings

- [numerical] Without the smoothness condition the optimizer reached O < 0.03 easily. Rechecks at L = 14, 16, 18 showed those bands to be jagged: the lowest level alternates between even and odd m (period-2 structure or level crossings), and the "valleys" moved by up to 0.7 between sizes. They were artifacts. The smoothness condition removes them (tail weight 0.19 against the 0.02 threshold for the best artifact).
- [numerical, local/qopt_all_pilot2.jsonl] With the full conditions, 8 restarts of 150 evaluations gave smooth, isolated, near-degenerate two-valley bands in N, for example k_A = 0.43, k_B = 2.15, E = 0.273 and 0.275 with bandwidth 1.04 (dm wedge, eps about 0.6), and k_A = 5.78, k_B = 3.80, E = 0.580 and 0.581. Accidental two-valley bands therefore exist at moderate gaps.
- [numerical, local/qtune_p*.jsonl] Continuing toward onset, the barrier shrinks and the valleys move together. The best endpoint reached gap/bw = 0.066 with k_A = 0.41, k_B = 1.14. Its L = 14 band is flat to 15 percent of the gap from k = 0.45 to 1.8 (E1 = 0.081, 0.084, 0.078, 0.092 at m = 1 to 4), one broad flat valley with ripples rather than two separated valleys. The other three endpoints stalled at gap/bw of 0.38 to 0.87.
- [numerical, twobody.py at that endpoint] The lowest two-particle levels on the 14-site ring lie above the non-interacting floor by +0.20 (AA), +0.09 (BB) and +0.27 (AB), with no bound state below the floor in these sectors. A small binding of -0.012 appears only at total momentum 0 from k = 0 particles, which are not valley states. At this ring size the shifts are mostly confinement energy of point-like hard-core particles, so no scattering length or J_z can be extracted. [heuristic] So the like-valley channel does not bind, but the sign of J_z is not determined locally.

[heuristic] Why the valleys merge. In a one-site nearest-neighbour cell, the LSW band has only first harmonics (C160 4.2), and two minima need a second harmonic generated by quantum dressing. Dressing is strong only where the first harmonics are small, that is near a flat band. Lowering the gap at fixed degeneracy seems to drive the band back toward a single broad minimum. This matches C160's "one species per site" heuristic, but it is an observation over a few continuations, not a theorem.

## 5. Mass imbalance (the correction C162 left undetermined)

[numerical, masscheck.py, masscheck.out] Impenetrable two-species lattice gas on a ring, no passing, hopping t_A = 1 + eta and t_B = 1 - eta, and nothing else. This is the leading dilute description of two valleys with masses 1/(2 t).

    E(ABAB) - E(AABB) = +3.4e-4, +3.2e-3, +1.4e-2   at eta = 0.1, 0.3, 0.6   (Lr = 24, N = 4).

The difference scales as eta^2 (ratio 9.2 between eta = 0.1 and 0.3) and as n^2 (ratio 2.25 between n = 1/6 and 1/9). With six particles the energies are ordered by the number of A-B interfaces: -10.0198 (AAABBB, 2 interfaces), -10.0058 (AABABB, 4), -9.9927 (ABABAB, 6), roughly additive. So each A-B interface costs about 0.3 eta^2 hbar^2 n^2 / M (lattice values at n <= 1/6).

[derived] Mass imbalance favours segregation, with a ferromagnetic pseudospin coupling of order eta^2 n^2/M per bond. The scattering-length coupling of C162 (1.3) is of order n^3 Delta a / M. So in the strict dilute limit every mass-imbalanced accidental two-valley gas is immiscible, Delta -> -infinity, whatever the scattering lengths. A miscible (or Neel) pseudospin sector is possible only where n Delta a exceeds about 0.1 eta^2 with the miscible sign of Delta a, that is at finite density or with masses tuned close to equal. This settles C162's heuristic in its favour, and it shrinks (A_sol) to a density window, without closing it. [heuristic] The mechanism is the positive surface (Dirichlet) energy of a fermionized gas at a heavy-light interface. The Born-Oppenheimer limit gives the same sign (two light particles in one box cost 5 units, against 8 when split into two half boxes).

## 6. Direct test past the onset (Task 3)

onset_dmrg.py (version 2, uploaded to Typhon). At a tuned endpoint it computes the directions u_mu (both valleys down) and u_lam (valley splitting) from the L = 12/14 Jacobian. It then runs a coarse scan at depth 2 delta (delta = 0.05 bandwidth) over three splittings, bisects 7 times in lam on the difference of helicity-resolved valley weights, and finishes with L = 64 and 128 at the switching point lam* plus one deeper point. The observables are entanglement fits (C155 profile_fits), connected correlations C^{ab}(r) from L/4 (checked against ED on 10 sites to 6e-14, test_corr.out), and the helicity-resolved transverse structure factor S_+(q). That last one separates q from -q. The eigenvalues of S^{ab}(q) do not, because S(-q) is the complex conjugate of S(q) for real correlators. Valley weights are the maxima of S_+ within 0.3 of k_A and k_B.

Criteria, fixed before any data. A c = 2 candidate needs entanglement fits above 1.5 that do not fall from L = 128 to 256, truncation below 1e-7, two incommensurate S_+ peaks at k_A and k_B with comparable weights at lam*, and an FES slope above 1.5 at chi >= 256. Points with c_prof(64) > 1.2 are promoted to L = 128 automatically. L = 256 and FES are follow-ups for any survivor.

[numerical, smoke tests only, L = 16 to 48] The pipeline runs. At L = 48 on the flat-bottom endpoint it gives c_prof = 0.71 with a single valley peak near |k_A| on one side of the splitting, and a gapped state (c = 0) on the other. The linear extrapolation of the ring band is a poor guide: the deepest point (4 delta) came out gapped. This is why the job scans and bisects rather than trusting one point.

## 7. Scope (Task 4)

[derived estimates]
- Soliton length. Classical K solitons are 10 to 28 sites at eps = 0.3 to 0.5. The quantum two-valley particles found in section 4 are lattice scale, since band smoothness at L = 12 to 14 needs a particle size of at most about 3 sites.
- Pseudospin scale. With M of order 1 (band curvature of order the bandwidth), depth mu = 2 delta = 0.1 bandwidth, n = sqrt(2 M mu)/pi is about 0.1 to 0.2. The scattering-length coupling is J of order (pi^2/3) n^3 Delta a / M, about 3e-3 to 3e-2 in bond units for Delta a of one site. The mass term is about 0.3 eta^2 n^2/M, about 6e-4 at eta = 0.3, n = 0.15. Valley exchange is not tunnelling-suppressed for point-like particles, so |J_perp| is of the same order as J.
- System size. L = 64 holds 6 to 13 particles and L = 128 holds 13 to 26, which is enough for a two-component entanglement signal if the pseudospin sector is critical (its contribution goes like (1/6) ln N_p, independent of its velocity).
- Visibility. A miscible window in the splitting has width of order J, about 3e-3 or more. The bisection resolves 4 delta / 2^7, about 1.6e-3 of the bandwidth, so such a window is resolvable but only marginally at its lower end. A window much below 1e-3 would be missed, and the null would then bound the window width by about 1e-3 of the bandwidth at the points tested.
- The blow-up regime (eps <= 0.1) remains out of reach, as C162 estimated (pseudospin scale about 1e-7).

## 8. Theorem-style summary

- [proved, literature] In T every uniform gapped background breaks Theta (LSM). Two-valley bands in T are accidental.
- [derived] Emergent two-charge structure needs k_A, k_B not in {0, pi} and k_A + k_B != 0 mod 2 pi. The K valleys of C161 (sum 2 pi) fail the second condition at eps of order one [heuristic consequence].
- [proved, numerical check] The Pi-symmetric, a_c = 0 surface (C152's cos phi = 0) carries an exact symmetry k -> pi - k on one-particle states, a protected incommensurate two-valley band. It is the two-Majorana case of C152, c <= 1 there.
- [numerical] D = 0 bonds show no incommensurate minima in 195 uniform samples (quantum analogue of Lemma 3').
- [numerical] The classical lattice code reproduces C161's tilt eigenvalue (-0.0671 at eps = 0.05) and C162's K binding (U_same(10) = -0.0084 against -0.0088). The classical pilot screen of 188 bonds found genuine soliton minima but no second valley.
- [numerical] Accidental, smooth, isolated two-valley quantum bands exist in N at gap/bw of about 0.25. Continuation toward onset merges or flattens them (best gap/bw 0.066, flat bottom). No two-valley onset with separated valleys was reached locally.
- [numerical, derived] Mass imbalance favours segregation, with an interface cost of about 0.3 eta^2 n^2/M. It dominates the dilute limit, so accidental two-valley gases with unequal masses are immiscible as n -> 0.
- No explicit c = 2 bond. The direct test is in Typhon job 17601039, with pre-fixed criteria.
- Residual (A_sol') [assumed]: no accidental two-valley onset in N or T at eps of order one, at finite density n with n Delta a > 0.1 eta^2 (or near-equal valley masses), has separated, non-binding valleys and J_z > -|J_perp|. This is narrower than C162's (A_sol): the mass-imbalance sign is now fixed, T is reduced to N-type accidents by LSM, and the protected staggered case is identified and covered by C152.

## 9. What is not done

- The large classical and quantum screens, the continuation of every optimizer result, the past-onset DMRG and the controls at L >= 64 are in the Typhon job and had not run at write-up.
- The quantum soliton band at eps of 0.2 to 0.4, where solitons are 10 to 30 sites, needs ring sizes beyond ED (L >= 30) or an MPS excitation ansatz. It was not computed.
- Scattering lengths at a two-valley onset (and thus J_z directly) were not extracted. Rings of 14 to 16 sites are too small, and the DMRG test measures the outcome (miscible or not) instead.
- Continuations toward onset are few (4 locally). That the valleys merge near onset is an observation, not a result.
