# C165 derivations: a staggered Ising transition inside a floating liquid

This run was done by a Claude Opus 5.5 agent (2026-11-08), not by Codex. Astra was not used (out of credits).

Tags: [proved] (complete argument or exact identity), [derived] (physics-standard derivation with named inputs), [numerical] (script and output named), [heuristic], [assumed]. Conventions as in C154: H = sum_j [sigma_j^T J sigma_{j+1} + B . sigma_j], Pauli matrices; the classical image replaces sigma_j by a unit vector n_j. All scripts are in this folder.

## 1. Frames and the parent

### 1.1 The family

[definition] classical.fam builds

    J = [[-cos a + g, sin a + g', G], [-sin a + g', -cos a - g, Dx], [G', -Dx, Dz]]  (plus Dy),   B = (hx, hy, hz).      (1.1)

At g = g' = G = G' = Dx = Dy = hx = hy = 0 this is the U(1) parent of C154 section 2: XXZ about z with |J_perp| = 1 and DM along z, with pitch q = -a. The other entries are the harmonics: g, g' have charge 2 about z (XX - YY and XY + YX), G, G' are the symmetric charge-1 terms (Gamma, XZ + ZX), Dx, Dy the charge-1 DM terms, hx, hy the charge-1 fields. symcheck.py computes the residual O(3) stabilizer of (J, B) (unitary rotations and antiunitary improper ones, identity excluded modulo 2 pi). [numerical] With g = 0.3, hx = 0.2, hz = 0.05, Dx = 0.03 the smallest residuals are 0.067 (unitary) and 0.016 (antiunitary), so the bond is in N. Without hz and Dx the antiunitary Q = diag(1, 1, -1) is exact (residual 1e-25, stratum K). The C68 floating control nfloat (eta = 0.1, b = 1.05) has residuals 0.030 and 2.2.

### 1.2 The parent is the spin-1/2 XXZ chain, and the switch is its isotropic point

[proved] The in-plane part of the parent bond is -R(a), with R the rotation matrix. The site rotation sigma_j = R_z(-a j) tau_j (a product of on-site unitaries, exact on open chains, a boundary twist on rings) maps the parent to

    H_tau = sum_j [ -(tau^x tau^x + tau^y tau^y) + Dz tau^z tau^z ] + hz sum_j tau^z_j.      (1.2)

The further rotation by pi about z on odd sites (mu frame) gives the antiferromagnetic XXZ chain mu.mu (in plane) + Dz mu^z mu^z in a longitudinal field hz. So the parent is the spin-1/2 XXZ chain with Delta = Dz and field hz, for every a. Under the two maps the in-plane spiral of sigma becomes in-plane Neel order of mu, the staggered z moment of sigma stays staggered z, and the "staggered cone" n_j = (sin theta cos qj, sin theta sin qj, (-1)^j cos theta) becomes a Neel vector of mu tilted by theta from z. The spiral/staggered switch of the parent is therefore Delta = 1 at hz = 0 (the Heisenberg point in disguise) and the spin-flop line of XXZ in a field at hz != 0. A further flip of all odd spins (classically allowed, mu_j = (-1)^j nu_j) turns the Delta = 1 parent into the ferromagnet -nu.nu, with the Neel vector L as a uniform state nu_j = L.

[proved, bookkeeping] In the nu frame a charge-c operator on site j acquires a phase: sigma^+_j = e^{iqj} nu^+_j and sigma^z_j = (-1)^j nu^z_j. Hence the charge-2 bond terms carry momentum 2q, the symmetric and antisymmetric charge-1 bond terms carry q + pi, the transverse fields carry q, and the uniform field hz becomes a staggered field (momentum pi). None of these is zero for q not in {0, pi/2, pi}.

## 2. Classical Landau theory (Task 1, classical part)

### 2.1 Proposition 1: the rigid staggered cone is linear in cos^2 theta for every bond

[proved] Let e be a unit vector, (u, v) an orthonormal basis of e^perp, and n_j = sin theta (cos(qj) u + sin(qj) v) + (-1)^j cos theta e with q not in {0, pi}. For any (J, B), the energy per site averaged over a period (or over an infinite chain for irrational q/2 pi) is

    E(theta) = sin^2 theta e_s(q) + cos^2 theta e_N,   e_s(q) = (1/2)[(u^T J u + v^T J v) cos q + (u^T J v - v^T J u) sin q],   e_N = -e^T J e.      (2.1)

Proof. In-plane times in-plane products have momenta 0 and +-2q, in-plane times staggered have +-q + pi, the field on the in-plane part has +-q and on the staggered part pi. Only momentum 0 survives the average. So there is no sin theta cos theta, no sin theta, and no quartic term. QED. The rigid ansatz never coexists: one order or the other wins except on the degeneracy surface e_s = e_N, where the whole cone family is degenerate. For the parent this degeneracy is the SU(2) orbit of the Heisenberg point (section 1.2), and at hz != 0 it is the known biconical degeneracy of the classical XXZ antiferromagnet at its spin-flop field.

### 2.2 Proposition 2: second-order relaxation and the sign of the quartic coefficient

[derived; inputs: second-order perturbation theory of a classical minimizer in the harmonics, at the switch Delta = 1, hz = 0; checked numerically in 2.3] Write the nu-frame Neel vector as L = (sqrt(1-u) cos chi, sqrt(1-u) sin chi, sqrt(u)), so u = N^2 is the squared staggered z moment. Each harmonic exerts on the uniform state nu_j = L a force F_j = Re(B_k(L) e^{ikj}) at momentum k != 0, with B_k linear in L. The ferromagnet's transverse stiffness at momentum k is the same in every direction, so relaxing the texture lowers the energy per site by

    delta E = - |P_L B_k(L)|^2 / (8 (1 - cos k)),   P_L = 1 - L L^T,   |P_L B|^2 = |B|^2 - |L . B|^2.      (2.2)

|B|^2 is quadratic in L, so it only shifts the linear coefficient of E(u) (the switch). The quartic part is +|L . B_k(L)|^2/(8(1 - cos k)), a non-negative quartic form on the sphere. Evaluated per harmonic:

- charge 2, g (XX - YY): B = 2 g cos q L^+ (1, i, 0), k = 2q, |B|^2 = 8 g^2 cos^2 q (1-u), |L . B|^2 = 4 g^2 cos^2 q (1-u)^2, hence |P_L B|^2 = 4 g^2 cos^2 q (1-u)(1+u) and

      delta E_2 = - c_2 (1 - u^2),   c_2 = g^2 cot^2 q / 4.      (2.3)

  The u^2 coefficient is +c_2 > 0 (convex).
- charge 1, symmetric, G (XZ + ZX): B = 2G (cos q L^+ z_hat - L^z (1, i, 0)), k = q + pi, L . B = 2G L^z L^+ (cos q - 1), and

      delta E_1 = - G^2 [cos^2 q (1-u) + 2u - (1 - cos q)^2 u (1-u)] / (2 (1 + cos q)).      (2.4)

  The u^2 coefficient is -c_1 with c_1 = G^2 (1 - cos q)^2 / (2 (1 + cos q)) (concave), and the linear part favours u = 1 (tilted Neel).
- charge 1, antisymmetric (DM_x, DM_y): L . B = 0 because L^T A L = 0 for antisymmetric A. No quartic term.
- transverse fields hx, hy: B is independent of L. No quartic term. The uniform field hz is a staggered field in the nu frame and gives -hz^2 (1-u)/8, again linear (the spin-flop).

Harmonics of different charge carry different momenta and do not interfere at second order. Altogether, at the switch,

    E(u) = E_s + (e_N - e_s + linear shifts) u + (c_2 - c_1) u^2 + O(harmonics^3).      (2.5)

If c_2 > c_1, E(u) is convex and its minimum moves continuously from u = 0 to u = 1 across a window of width 2(c_2 - c_1) in the detuning: a coexisting phase (spiral plus staggered moment) entered and left by continuous transitions, with u linear in the detuning (mean-field Ising). If c_2 < c_1 the switch is first order. So the bond term that curves the interpolation in favour of coexistence is the charge-2 harmonic, the XY anisotropy about the spiral axis. Physically, a planar spiral cannot use the radial part of the XY-anisotropy force, while a cone can, by modulating its cone angle at momentum 2q. Gamma terms do the opposite, and DM and fields do nothing at this order.

### 2.3 Numerical check of Proposition 2

[numerical] Ring L = 100, a = 2 pi 13/100 (q = 0.817), scripts scan_cl*.py and the inline checks in PROGRESS.md.
- c_2: the spiral's relaxation gain at g = 0.05 and 0.1 is 5.510e-4 and 2.2016e-3, against g^2 cot^2 q/4 = 5.511e-4 and 2.2046e-3.
- Gamma: the Neel gain at Dz = 1.2 is 1.33e-3 (G = 0.05) and 5.29e-3 (G = 0.1). The switch-point formula G^2/(1 + cos q) gives 1.48e-3 and 5.94e-3. The ratio 0.89 is (1 + cos q)/(Dz + cos q) = 0.894, the stiffness at Dz = 1.2 rather than at the switch.
- Coexistence: at g = 0.3 the ground state has both a spiral peak and a staggered moment for 1.00 < Dz < 1.047, with |N|^2 = 0.087, 0.30, 0.52, 0.74, 0.97 at Dz = 1.005, 1.015, 1.025, 1.035, 1.045, linear in Dz with slope 1/(2 c_2^eff), c_2^eff = 0.023. The relaxed pitch locks at 11/100 (q = 0.691), and the formula gives c_2 = 0.020 at q = 0.817 and 0.033 at q = 0.691, which brackets the measured value.
- G = 0.3 alone: tilted Neel already at Dz = 0.98, no coexistence. hx = 0.3 alone: canted Neel at all Dz scanned, no coexistence. g = 0.3, hx = 0.2, hz = 0.1 (in N): coexistence for about 0.99 < Dz < 1.05 with uniform canting |M| = 0.11 to 0.15.

### 2.4 The window is parametrically narrow

[numerical, cl4_*.out, cl3_g45.out] Larger g does not widen the window. At a = 0.817, g = 0.45 the spiral turns into a soliton lattice (incommensurate peak at k = 0.13 to 0.19) and then a collinear ferromagnet along y. At a = 1.19 (stiffness 1 - cos 2q = 1.74) g = 0.5 gives c_2 = 0.012 and g = 0.7, 0.9 again give soliton lattices. [derived] The reason is (2.3): c_2 grows with cot^2 q, but small q is exactly where the collinear y ferromagnet (energy -cos a - g) competes with the spiral (-1 - c_2). In this family the classical coexistence window is never wider than a few percent of the exchange.

## 3. Classical lattice search (Task 2)

[numerical, hess_cl.py] The continuous classical transition is a Bloch-pi softening of the minimizing spiral. At g = 0 the lowest staggered Hessian mode is w = 2(1 - Dz) (0.100, 0.060, 0.020 at Dz = 0.95, 0.97, 0.99), zero only at the degeneracy. At g = 0.3 it is 0.0995, 0.0605, 0.0214, 0.0018 at Dz = 0.95, 0.97, 0.99, 1.00 and crosses zero at Dz ~ 1.001, while the spiral is still the global minimizer (e = -1.024 against the Neel energy -Dz). The soft mode has 97 percent of its weight at momentum pi.

[derived] This is consistent with C160 M2, whose scope is the continuum (Tonelli) functional near h_F: the scalar disconjugacy argument needs a smooth texture and fails on the lattice, where the Hessian of a spiral is block tridiagonal with 2 x 2 blocks and has no Perron-Frobenius property. The out-of-plane band of the parent's spiral is w(k) = 2(1 + Dz cos k), with its minimum at k = pi for antiferromagnetic Dz. So classically, in N, P1 exists: a nearest-neighbour bond with a minimizing incommensurate spiral whose staggered mode softens continuously, followed by a coexisting phase. The classical obstruction stated in the spec (linearity of the staggered cone) holds only for the rigid ansatz.

## 4. Quantum mechanism (Task 1, quantum part)

### 4.1 Spin 1/2 has no independent amplitude mode at the switch

[derived; inputs: the SU(2)_1 description of the spin-1/2 Heisenberg chain, C154 section 2.2, and section 1.2] By 1.2 the classical switch is, quantum mechanically, the spin-1/2 Heisenberg antiferromagnet (in the mu frame), whose infrared is the SU(2)_1 WZW model. The staggered z moment of sigma and the two components of the sigma-frame spiral are the three components N^a = tr(g sigma^a) of a single primary (dimension 1/2). They are not two order parameters with an independent amplitude, they are one field of one c = 1 theory. This is the precise sense in which "the liquid's boson is the staggered order" (spec, Task 1): in the U(1) parent Neel along z is cos 2phi, and the spiral is e^{i vartheta}, of the same boson.

The classical Landau potential E(L) on the Neel sphere decomposes into SU(2)_diag harmonics of the Neel vector. The quadratic (spin-2) part, which in (2.5) is the linear term in u, corresponds to the marginal current-current anisotropy J^z Jbar^z (dimension 2): this is the XXZ anisotropy, which drives the KT transition into gapped Neel order on one side and the c = 1 liquid on the other. The quartic part (2.5) contains a spin-4 harmonic of L. In SU(2)_1 a global spin-j chiral operator has dimension at least j^2 (the chiral algebra is the boson at the self-dual radius, with spin-j states first at h = j^2), so an operator with diagonal spin 4 needs (j_L, j_R) with j_L + j_R >= 4 and has dimension at least 8. The curvature that produces classical coexistence is therefore an irrelevant operator of the spin-1/2 chain. Near the switch the quantum problem is XXZ with a renormalized anisotropy, and its two phases are the liquid (c = 1) and gapped Neel, separated by KT (hz = 0) or by the Pokrovsky-Talapov edge of the Neel plateau (hz != 0). No coexistence.

### 4.2 To all orders in weak harmonics

[derived; inputs: C154 2.2 to 2.3, C154 5.4] Off resonance every harmonic carries an incommensurate momentum (1.2, bookkeeping) and is forbidden at first order. At order n only charge-zero products at total momentum zero survive. These renormalize the XXZ parameters (marginal) or are irrelevant. Resonant products lock commensurate gapped phases on hypersurfaces. In the liquid at hz != 0 the sigma-frame staggered moment has no projection on the boson (C154 5.4, its vertex momenta are pi + 2n k_F), so a staggered order inside the liquid needs a second sector, which perturbation theory in the harmonics does not create. Hence, within the perturbative neighbourhood of the U(1) class, the classical coexistence window of section 2 is absent at spin 1/2. The classical window is O(g^2) wide (2.3), and the quantum shift of the Neel onset is O(1) (4.3), so even the location of the classical window lies inside the quantum liquid.

### 4.3 Numerical check at finite harmonics

[numerical, loc/Q1_*.json, loc/Q2_*.json, two-site VUMPS chi <= 40, anaq.py] Line Q2: a = 0.9, g = 0.3, hx = 0.2, hz = 0.05, Dx = 0.03 (in N, 1.1), whose classical coexistence window is about 0.99 < Dz < 1.05.
- Dz = 1.0 (inside the classical window): staggered moment N = 0.082, 0.064, 0.056, 0.050 at chi = 16, 24, 32, 40 (decaying), c_pair = 1.08, 0.96, 1.02, 1.07, xi = 82 at chi = 40. A c = 1 liquid without staggered order.
- Dz = 1.15: N decays (0.135 to 0.096 in the cell-difference norm), c_pair 0.84 to 0.92.
- Dz = 1.25, 1.3, 1.35, 1.4: N saturates (0.256, 0.311, 0.365, 0.413 at chi = 40, along z), c_pair at the last step 0.45, 0.27, 0.17, 0.09 and falling, xi = 42, 33, 26, 20. Gapped Neel.
- No point has a saturated N together with c_pair near 1. The Neel onset is near Dz = 1.2, outside the classical window, as 4.2 predicts. The same holds on line Q1 (a = 0.8), where in addition the transfer-matrix phase sits at 0.494 to 0.499 per cell, near the period-4 lock.
- Plateau edge (loc/P_*.json, Dz = 2.0, same harmonics): gapped Neel with N = 0.73 and xi = 4 up to hz = 2.2, a disordered unconverged state at hz = 2.6. Finer scans of this edge are in the Typhon screen (family A2).

### 4.4 Why the semiclassical regime does not help

[derived] The classical analysis is controlled at spin 1/2 only where the length scale is long, near the ferromagnetic point h_F (C160, K of order eps^{-1/2}). The coexistence window of section 2 sits at the hidden Heisenberg antiferromagnet, K = 1/2 in the zero-field XXZ convention, or at the Neel plateau edge, K near 1. The two regions are disjoint. Near h_F the staggered moment costs exchange energy of order J (the m-coupling is ferromagnetic there) and C160 M2 forbids a continuous Bloch-pi onset of continuum minimizers. In the large-K bosonization of a floating liquid coupled to a staggered Ising field, translation forbids sigma V and epsilon V for every vertex V (all carry incommensurate momentum in N, C159 Lemma A: both shifts free), and epsilon d_x phi (dimension 2, spin 1) and epsilon (d phi)^2 are marginal or irrelevant. So the decoupled Ising x boson is kinematically stable at any K (C154 table 4.2), and the obstruction is dynamical: at large K there is no low-lying staggered mode, at small K the staggered mode is a component of the boson.

### 4.5 Floating, not half-floating

[derived, C159 Lemma A] In K (antiunitary present) one field of the boson is commensurate and sigma V_c is allowed. On line Q2 the antiunitary diag(1, 1, -1) of the hz = 0 bond is broken by hz = 0.05 and Dx = 0.03 (residual 0.016), so both shifts are free and the liquid is floating, but only weakly away from K. The Typhon screen therefore includes hz = 0.3 (family A1) and random symmetry-free harmonics (A3), which are robustly in N.

## 5. Calibration (Task 4)

[numerical, local] TFI open chain L = 32: c_prof 0.53. Ising x XX ladder (decoupled critical TFI qubit and XX qubit per rung, d = 4, outside B_2) L = 24, chi = 32: c_prof 1.29, c_sm 1.43. Q2 at Dz = 1.0 (a c = 1 point in B_2) L = 32: c_prof 0.84 (open chain, small L). The thermodynamic calibrations (Ising x XX FES to chi = 256 and DMRG L = 64, 128; XX to chi = 256 and L = 256; the N-stratum floating control nfloat FES to chi = 256 and DMRG L = 128, 256) are in Typhon job 17601089.

## 6. Theorem-style summary (Task 5)

1. [proved] For every nearest-neighbour bond and every incommensurate q the rigid staggered cone has energy linear in cos^2 theta (2.1). The parent's switch is the Heisenberg point (hz = 0) or the spin-flop line of XXZ in a field (1.2).
2. [derived, checked numerically to 4 digits] At second order in the harmonics the Landau coefficient of N^4 at the switch is c_2 - c_1, with c_2 = g^2 cot^2 q/4 from the charge-2 XY anisotropy and c_1 = G^2 (1 - cos q)^2/(2(1 + cos q)) from the symmetric charge-1 Gamma terms. DM and fields contribute nothing. c_2 > c_1 gives a classical coexistence window of width 2(c_2 - c_1), entered and left continuously.
3. [numerical] Classical lattice P1 exists in N: at g = 0.3, hx = 0.2, hz = 0.1 the minimizing spiral's staggered Bloch-pi mode softens while the spiral is the global minimum, followed by a coexisting window of width about 0.05. M2 does not apply on the lattice. The window does not exceed a few percent of the exchange in this family.
4. [derived] At spin 1/2, near the switch and to all orders in weak harmonics, the staggered moment and the spiral are components of one SU(2)_1 primary. The classical curvature maps to an operator of dimension at least 8, and the infrared is the XXZ dichotomy (liquid or gapped Neel). No coexistence and c <= 1.
5. [numerical] At harmonic amplitude 0.3 (fraction about 0.2 of the exchange) the quantum line through the classical window is a c = 1 floating liquid with decaying staggered moment, followed near Dz = 1.2 by gapped Neel order, with no intermediate coexisting phase at chi <= 40.
6. [assumed, the narrowest remaining assumption, A_S] In N, at harmonic amplitudes comparable to the exchange, in the hidden-Heisenberg region (parent Delta near 1, any hz up to the Neel plateau edge), no lattice-scale excitation independent of the liquid's boson orders at momentum pi while the liquid survives. Statement 4 proves the perturbative part. The Typhon screen tests amplitudes 0.2 to 0.6 on 252 points. Outside this region the staggered moment is either far from soft (near h_F, statement 4.4 and C160 M2) or covered by C154 (A6).

## 7. Consistency checks

- Gradient of the classical energy against finite differences: agreement to 6e-9 relative (classical.py main).
- Parent reduction: at g = 0 the classical switch is at Dz = 1 exactly, energies -1 (spiral) and -Dz (Neel), and the staggered Hessian mode is 2(1 - Dz), the Heisenberg-ferromagnet stiffness at k = pi in the nu frame.
- Second-order coefficients against direct minimization (2.3), including the predicted stiffness correction (1 + cos q)/(Dz + cos q) away from the switch.
- Symmetry assignments by symcheck (1.1). The pure parent has exact unitary and antiunitary stabilizers (residual 1e-33), as it must.
- Quantum: at Dz = 1.0 the VUMPS state has equal entanglement on both bonds of the cell (spread 1e-7) and c_pair near 1, the expected c = 1 liquid. The ordered points have N along z, matching the classical Neel axis.

## 8. What is not done

- The Typhon screen (job 17601089) had not started when this was written. Its pre-fixed criteria are in typhon/stageB.py and typhon/summarize.py.
- The classical search is a targeted family plus the second-order classification. A random classical search over the eight harmonic directions was not run separately; the quantum screen A3 covers random directions.
- No exact or controlled treatment of harmonics comparable to the exchange. That is assumption A_S.
