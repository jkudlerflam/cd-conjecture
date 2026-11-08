# C160 derivations: c <= 1 near the ferromagnetic Heisenberg point

This run was done by a Claude Opus agent (2026-11-08), not by Codex. Astra was not used (out of credits).

Standard: theoretical physics. Tags: [proved] (complete argument, or an exact identity checked to machine precision), [derived] (physics-standard derivation with named inputs), [numerical] (sizes, scripts, outputs named), [heuristic], [assumed]. A tag applies to the paragraph or displayed claim it opens.

Conventions as in C154 and C159. H = sum_j [ sigma_j^T J sigma_{j+1} + B . sigma_j ], Pauli matrices, J = S + A(D) with A(D)_ab = eps_abc D_c, so that sigma_j^T A(D) sigma_{j+1} = D . (sigma_j x sigma_{j+1}). est.py coordinates v = (J_x, J_y, J_z, D_x, D_y, D_z, f_x, f_y, f_z). The spin is s = 1/2 per site and sigma = 2S. h_F = -(XX + YY + ZZ), i.e. J = -1, B = 0.

## 1. Blow-up at h_F

### 1.1 Coordinates

[proved] Every bond near h_F is, up to a positive scale and a common rotation,

    J = -1 + eps (S0 + A(D)),   B = eps b,   S0 symmetric traceless, |(S0, D, b)| = 1,      (1.1)

because the trace part of J only rescales to first order. Rotations fix h_F and act on u = (S0, D, b) by conjugation (S0 -> R^T S0 R, D -> R^T D, b -> R^T b). So the directions form S^10 / SO(3), of dimension 7, and with eps the neighbourhood has dimension 8, as B_2 does.

### 1.2 Classical functional

[derived; input: spin coherent states] For a product state with unit vectors n_j, <sigma_j^a sigma_{j+1}^b> = n_j^a n_{j+1}^b, so the classical energy is E = sum_j [ n_j^T J n_{j+1} + B . n_j ]. With -n_j . n_{j+1} = -1 + |n_{j+1} - n_j|^2 / 2, n_j^T S0 n_{j+1} = n^T S0 n + (total derivative) + O(eps |n'|^2), and D . (n_j x n_{j+1}) = D . (n x n') + O(n'^2), the gradient expansion is

    E[n] = int dx [ |n'|^2 / 2 + eps D . (n x n') + eps V(n) ],   V(n) = n^T S0 n + b . n,      (1.2)

up to a constant and corrections of relative order eps (stiffness renormalization) and k^2 (lattice dispersion). The DM term is the Lifshitz invariant. It is a connection: |n' + eps D x n|^2 / 2 = |n'|^2/2 + eps D . (n x n') + eps^2 |D x n|^2 / 2, so

    E[n] = int dx [ |n' + eps D x n|^2 / 2 + eps V(n) - (eps^2/2) |D x n|^2 ].      (1.3)

The rotating frame n = R_D(eps |D| x) m (rotation about D-hat) removes the connection, turns -(eps^2/2)|D x n|^2 into an easy-plane anisotropy perpendicular to D of order eps^2, and makes the parts of V that are not axial about D rotate with wave vector eps|D| and 2 eps|D|. This is the continuum form of the C154 spiral frame. Near h_F the parent of C154 is a ferromagnetic XXZ chain with anisotropy of order eps (from V) plus eps^2 (from the DM), and the harmonics have amplitude eps and momentum eps|D|, 2 eps|D|.

[proved] Check at the pure DM direction, u = (0, D, 0). J = -1 + eps A(D) has J J^T = 1 + eps^2 (|D|^2 - D D^T), so in the C154 spiral frame the parent is XXZ with J_perp = -sqrt(1 + eps^2 |D|^2) and J_par = -1, i.e. Delta_AF = -1/sqrt(1 + eps^2 |D|^2) > -1. This is the U(1) class, c = 1, with K -> infinity as eps -> 0, as (1.3) predicts (easy plane of strength eps^2 |D|^2 / 2 perpendicular to D).

### 1.3 Quantum theory and the semiclassical parameter

[derived] The coherent-state path integral is

    Z = int Dn exp{ i int dt [ s sum_j omega(n_j) . dn_j/dt - E[n] ] },      (1.4)

with omega the monopole connection and s = 1/2. The equation of motion is the Landau-Lifshitz equation s dn/dt = n x dE/dn (up to orientation). Linearizing (1.2) around n = e_z with no anisotropy gives s omega = k^2, i.e. omega = 2 k^2, which matches the exact one-magnon band 4(1 - cos k) of h_F in Pauli units (twomag.py, c0 + 2|a| cos(k - pi) with c0 = 4, |a| = 2).

Rescale x = l xi with l the shortest length set by (1.2). Every term of (1.2) has the same size eps per site at l = eps^{-1/2} when V carries the gap, and eps^2 at l = 1/eps in the pure-DM wedge (section 2.4). The action becomes

    S = (s l) int dtau dxi [ omega(n) . d_tau n - E_resc[n] ],     hbar_eff = 1/(s l) = 2 eps^{1/2}  (generic),  2 eps  (pure-DM wedge),      (1.5)

with E_resc of order one. So the blow-up is semiclassical in the sense of a large spin per coherence length, s l = eps^{-1/2}/2 or 1/(2 eps), and the loop expansion of (1.4) is an expansion in hbar_eff. This is the same mechanism that gives the XXZ Luttinger parameter K = pi/(2 sqrt(2 eta)) at Delta_AF = -1 + eta: K is proportional to 1/hbar_eff, because eta is of order eps (easy plane from V) or eps^2 (from the DM).

In the generic scaling the DM term of (1.2) enters E_resc with weight eps^{1/2}, the same order as hbar_eff, and V enters at order one. In the pure-DM wedge everything enters at order one.

## 2. Phases of the blow-up

### 2.1 Leading order: minimizers of V on S^2

[proved] Let a_1 <= a_2 <= a_3 be the eigenvalues of S0, e_i the eigenvectors, b_i = b . e_i. Minimizing V(n) = n^T S0 n + b . n on S^2 is the trust-region problem. The Lagrange condition is 2 (S0 - mu) n = -b, the global minimizer has mu <= a_1, and it is unique unless the hard case occurs (b_1 = 0 and sum_{i >= 2} b_i^2 / (4 (a_i - a_1)^2) <= 1). In the hard case with a_1 simple there are two minimizers n_0 +- t e_1, exchanged by the reflection n_1 -> -n_1. In the hard case with a_1 = a_2 the minimizers form a circle (a cone about e_3).

[proved] A circle of minimizers of a quadratic-plus-linear function on S^2 requires the m = 1 and m = 2 harmonics of V along that circle to vanish. If the circle is a cone about e_3 this holds when S0 is axial about e_3 and b is parallel to e_3 (the U(1)-symmetric case), and also, accidentally, when the m = 1 harmonic of the off-diagonal S0_{13}, S0_{23} cancels that of (b_1, b_2) at the cone angle. Both are of codimension 4 in (S0, b).

### 2.2 Phase list as a function of u at leading order

[derived; inputs: 2.1 and section 3]
- P (open set). Unique nondegenerate minimizer of V. Uniformly polarized or canted ferromagnet, magnon gap of order eps, c = 0. This is the generic case.
- F (codimension 1). Hard case with a_1 simple, strictly inside the trust region. Two canted minima, exchanged by n_1 -> -n_1. This is a first-order surface between two gapped ferromagnets (c = 0). On spins the reflection is Theta R_pi(e_1), exact when D is parallel to e_1 (the K stratum, and T at B = 0), broken at relative order eps^{1/2} otherwise.
- I (codimension 2, codimension 1 inside K). The edge of F, sum b_i^2/(4(a_i - a_1)^2) = 1, where the two minima merge (pitchfork). One real soft mode at k = 0 with a gapped conjugate, so a quantum Ising transition, c = 1/2. With D not along e_1 the reflection is only approximate and the edge is an Ising critical endpoint with emergent Z_2 (still codimension 2).
- C (codimension 4 in (S0, b), any D). Circle of minimizers. At leading order a U(1) degeneracy. The subleading terms decide (section 2.3).
- Z (the D-sphere, u = (0, D, 0)). V = 0. The U(1) class XXZ with Delta_AF -> -1^+, c = 1 (1.2).
- h_F itself is the excluded centre (SU(2), z = 2, degenerate multiplet).

[derived] In T (B = 0) near h_F one has b = 0, so the hard case holds with sum = 0 < 1. The ground states are the Theta-broken pair +-e_1 (c = 0) unless a_1 = a_2. There is no pitchfork I in T near h_F. This agrees with C159 Lemma B (no symmetric Ising point in T).

### 2.3 The wedge around C

[derived] Near C write V = V_ax + eta V_1, with V_ax having the circle of minima at polar angle theta_0 about e_3 and eta V_1 the perturbation. The polar mode has a gap of order eps (curvature of V_ax). Below that scale the theory is a single compact boson, the azimuth vartheta with conjugate density s(1 - cos theta), with

    H_eff = int dx [ rho_s (d_x vartheta - q)^2 / 2 + (chi^{-1}/2) delta m^2 + eps eta U(vartheta) + eps^2 U_D(vartheta) ],      (2.1)

where rho_s = sin^2 theta_0, q = eps D . e_3 (DM along the axis, the Lifshitz twist), U has harmonics 1 and 2 only (V is quadratic in n), and U_D (harmonic 2) comes from the transverse DM through (1.3). The Luttinger parameter is of order s l rho_s, large. The vertex operators e^{i p vartheta} have dimension p^2/(4K), of order hbar_eff, so they are strongly relevant and act classically. The phases are therefore those of the classical sine-Gordon problem with twist and harmonics 1, 2:
- pinned (commensurate, c = 0) when eps eta |U| is larger than about rho_s q^2,
- floating (chiral soliton lattice, an incommensurate c = 1 liquid with K of order 1/hbar_eff) otherwise,
- separated by Pokrovsky-Talapov transitions (z = 2), and inside the pinned region by double-sine-Gordon Ising lines (the minimum of U splitting into two, c = 1/2, vartheta gapped),
- at the edge of the cone (theta_0 -> 0) the polar mode softens and rho_s -> 0: z = 2 saturation onset, one species (a dilute Bose gas of magnons).
So the floating phases occupy a wedge of relative width eps around C, and inside it there is one scalar field. [derived] A scalar field cannot carry two soft modes in a periodic texture: the Hessian is a scalar Hill operator, its lowest band is minimal at Bloch momentum 0 with the nodeless phason as ground state (Floquet and Sturm theory), so no second gapless branch exists. Section 3 shows the same without the scalar reduction.

### 2.4 The wedge around Z

[derived] Near Z, i.e. for |(S0, b)| of order eps relative to |D|, the twist energy eps^2 |D|^2 and the potential eps V are comparable. Rescaling x = xi/(eps|D|) gives (1.2) with all three terms of order one: the general one-dimensional chiral ferromagnet

    E_resc[n] = int dxi [ |n'|^2/2 + D-hat . (n x n') + W(n) ],   W quadratic plus linear on S^2, arbitrary,      (2.2)

with hbar_eff = 2 eps. Its classical phase diagram contains uniform states, helices, cones, fans and soliton lattices, with first-order and continuous transitions between them. There is no scalar reduction. Section 3 handles (2.2) in general, and it contains all the previous cases as limits.

## 3. Soft modes of the classical ground states

Throughout, L(n, v) = |v|^2/2 + A(n) . v + W(n) on the tangent bundle of S^2, with A(n) = D x n (any D, including 0) and W any smooth function. The "time" of this mechanical problem is the chain coordinate x. The functional (1.2) and (2.2) are both of this form. The classical ground state minimizes the mean energy per length.

### 3.1 Lemma M1: classical ground states are uniform or singly periodic

[proved, given Mather's graph theorem] L is a Tonelli Lagrangian (fibrewise strictly convex and superlinear) on the compact manifold S^2. Every invariant probability measure that minimizes int L dmu (a Mather measure, whose value is the minimal mean energy) is supported on a fixed point of the Euler-Lagrange flow (a uniform state) or on one periodic orbit (a texture with a single period).

Proof. By Mather's graph theorem (J. Mather, Math. Z. 207, 169 (1991)), the union of the supports of minimizing measures projects injectively to S^2, and the inverse of the projection is Lipschitz. So on the compact projected set K the flow is generated by a Lipschitz vector field X. Extend X to a Lipschitz vector field on all of S^2 (McShane extension in R^3, then project to the tangent plane). Orbits starting in K are orbits of the extension, by uniqueness. Let mu be an ergodic minimizing measure. By Poincare recurrence, mu-almost every point is recurrent. For a flow on S^2 a recurrent point is a fixed point or lies on a periodic orbit (Poincare-Bendixson: a non-closed orbit of a planar flow cannot return arbitrarily close to its starting point, by the Jordan curve theorem). Hence mu is carried by a fixed point or by a periodic orbit. QED.

Consequences. [derived] (i) No classical ground texture of (1.2) or (2.2) is quasi-periodic with two incommensurate wave vectors, because the configuration space S^2 carries no invariant 2-torus that projects injectively. So no classical ground state has two independent phasons. (ii) The lattice corrections of relative order eps turn the flow into a nearby symplectic map. A hyperbolic minimizing periodic orbit (3.2) persists as an invariant circle, locked or incommensurate, with one phason. Elliptic islands, which could support two-frequency tori, do not exist around minimizers [heuristic for the lattice, persistence of normally hyperbolic invariant circles].

### 3.2 Lemma M2: no Bloch zero mode at k != 0

[derived; inputs: absence of conjugate points along globally minimizing orbits (Jacobi, Morse; for Mather orbits, Contreras and Iturriaga, "Global minimizers of autonomous Lagrangians", 1999), and the instability of minimizing closed orbits on surfaces (Poincare)] Let n_0(x) be a periodic ground texture of period T. The Hessian H of the energy at n_0 acts on tangent fields eta(x). H eta = 0 is the linearized Euler-Lagrange (Jacobi) equation. Bloch zero modes at quasi-momentum k correspond to Floquet multipliers e^{ikT} of the 4-dimensional linearized flow. These multipliers are {1, 1, lambda, 1/lambda} (the double 1 from the flow direction and the energy). At fixed energy the Jacobi fields transverse to the orbit obey a scalar Hill equation y'' + K(x) y = 0 (the 2-degree-of-freedom reduction, K the curvature of the Jacobi-Maupertuis-Randers metric, including the magnetic term). A minimizing orbit has no conjugate points, so this equation is disconjugate on the whole line, and a disconjugate periodic Hill equation has real positive multipliers (an elliptic or negative multiplier forces oscillating solutions). Therefore lambda > 0, and e^{ikT} with k != 0 is never a multiplier. So the Hessian has no zero mode at nonzero Bloch momentum: no soft mode at an incommensurate or period-doubled momentum relative to the texture.

### 3.3 Lemma M3: any extra zero mode at k = 0 is the canonical partner of the phason

[derived] At k = 0 the periodic zero modes are the phason n_0' and possibly one more, psi. There are two cases.

(a) psi comes from lambda = 1 (a bifurcation of the orbit). Then psi has zero linearized energy (it is the limit of eigenvectors with lambda != 1, which are symplectically orthogonal to the flow), so its normal component g solves the scalar Hill equation of 3.2 and is periodic. A periodic solution of a disconjugate equation has no zeros, so g > 0 up to sign.

(b) psi = d n_E / dE along the family of periodic orbits, which is periodic only if T'(E) = 0. Then psi changes the energy, and a two-mode Bloch analysis shows that the state is not a stable minimizer. In span{e^{ikx} n_0', e^{ikx} psi} the Hessian is [[kappa k^2, -i h k], [i h k, r]] to leading order, with kappa the phason stiffness, r the curvature along psi at fixed period and h the coupling between psi and compression. The family of stationary points at period T + delta T exists only if h != 0 (r psi + a psi^3 = h delta T), and at r = 0 the determinant k^2 (r kappa - h^2) < 0. Equivalently, with e(T) the energy density along the family, e''(T) = -E'(T)/T, so e'' diverges at T'(E) = 0 and the energy along the family is cubic there (e = e_0 - (T''/(3T)) (E - E_0)^3 + ...). So case (b) never occurs at a ground state. The instability that does occur first is r kappa = h^2, i.e. e''(T) = 0, a divergent compressibility of the soliton lattice (z = 2).

In case (a) the spin symplectic form of (1.4), Omega(u, v) = s int n_0 . (u x v) dx, pairs psi with the phason:

    Omega(n_0', psi) = s int_0^T psi . (n_0 x n_0') dx = s int_0^T g(x) |n_0'(x)| dx != 0,      (3.1)

because g has one sign and |n_0'| >= 0 is not identically zero. Geometrically (3.1) is the change of the solid angle enclosed by the orbit on S^2 under the deformation psi, and a nodeless normal deformation always changes it. So psi is the momentum conjugate to the phason. The two zero modes form one canonical pair with no restoring force at k = 0, and the mode has omega^2 proportional to (kappa k^2)(r + O(k^2)), which at r = 0 gives omega proportional to k^2 (z = 2). Away from the instability (r > 0) the phason is a linear mode, a single compact boson, c = 1.

[numerical, csl_tilt.out] On the K-stratum soliton lattice (planar CSL with field perpendicular to D) the out-of-plane Hessian is a scalar Hill operator. Its band minimum is at Bloch k = 0 for h = 0.05 to 0.55, its ground state is nodeless (min psi_0 from 0.031 to 0.059), and the normalized overlap with the phason is 0.995 to 1.000. These are Lemmas M2 and M3 in the planar case.

### 3.4 Uniform ground states: continuum form of Lemma 3

[proved, algebra] Around a uniform minimizer n* with tangent Hessian eigenvalues eps alpha, eps beta >= 0 (axes e_u, e_v) and DM component D_n = D . n*, the quadratic part of (1.2) is (u'^2 + v'^2)/2 + eps D_n (u v' - v u') + eps (alpha u^2 + beta v^2)/2 (the tangential DM components give a total derivative). The spin-wave frequency is

    s omega(k) = 2 eps D_n k + sqrt( (k^2 + eps alpha)(k^2 + eps beta) ).      (3.2)

Zeros satisfy y^2 + (eps(alpha + beta) - 4 eps^2 D_n^2) y + eps^2 alpha beta = 0 with y = k^2, a quadratic in y. At an onset (stable state, a zero) there is either one double root y_0 > 0, giving one incommensurate complex mode at k_0 = (eps^2 alpha beta)^{1/4} (the sign condition keeps one of +-k_0), or y = 0 with alpha beta = 0, giving one real mode at k = 0 (Ising, c = 1/2, when the other eigenvalue is positive) or, with alpha = beta = 0, one quadratic mode (z = 2). This is C154 Lemma 3 at scale eps, and it never gives two independent soft modes.

### 3.5 Counting

[derived from 3.1 to 3.4] At every classical ground state of (1.2) or (2.2) the soft sector is one of:
- none (gapped, c = 0);
- one linear canonical pair, the phason of a periodic texture or the Goldstone mode of a uniform state with an exact U(1) (one compact boson, c = 1);
- one real linear mode at k = 0 at a uniform pitchfork (c = 1/2);
- one quadratic canonical pair (z = 2: Pokrovsky-Talapov, saturation onset, divergent compressibility, or the case (a) pairing of 3.3), not conformal.
At first-order surfaces two macroscopically distinct classical states coexist, each with one of the above. Their critical endpoints are pitchforks of uniform states (Ising) or of textures (paired by 3.3, z = 2). So at most one gapless canonical pair is ever present, and two linear sectors (c = 3/2 or 2) never are.

The pattern P1 of C154 and C159 (Ising x half-floating boson, with a "floating + SSB" neighbour) needs a soft real mode independent of the phason. In a texture that mode would be either at k != 0 (excluded by M2) or at k = 0 (paired with the phason by M3, or not a minimizer). The pattern P3 (two floating bosons) needs two phasons (excluded by M1). P2 (Ising x boson with a uniform even order parameter) needs a soft internal mode at fixed density, which is case (b) of M3 (not a minimizer).

## 4. Task 2: magnon kinematics of polarized states

### 4.1 Lemma 3 refined

Setting of C154 Lemma 3: a uniform polarized reference along n, local frame (u, v, n), J' = R^T J R. Linear spin waves give H = sum_k [ A_k a_k^+ a_k + (g cos k a_k a_{-k} + h.c.) ] with

    A_k = A0 + a_c cos k + a_s sin k,  A0 = -4 J'_nn - 2 B'_n,  t = J'_uu + J'_vv + i (J'_vu - J'_uv),  a_c = 2 Re t,  a_s = -2 Im t,
    g = J'_uu - J'_vv - i (J'_uv + J'_vu),   omega(k) = a_s sin k + sqrt( Abar_k^2 - 4 |g|^2 cos^2 k ),  Abar_k = A0 + a_c cos k.      (4.1)

[numerical] lsw_minima.py checks (4.1) against direct diagonalization of the bosonic BdG problem on a 24-site ring to 7e-14.

Lemma 3'. [proved] If a_s = 0 (no DM component along n) and the reference is energetically stable (Abar_k >= 2|g||cos k| for all k), then every minimum of omega(k) lies at k = 0 or k = pi, unless omega vanishes identically. In particular the item "pair {k_1, -k_1} with a_s = 0" of C154 Lemma 3 never occurs for a stable reference, and no stable reference has two degenerate minima +-k_0 with k_0 not in {0, pi}.

Proof. Write x = cos k and F(x) = omega^2 = (A0 + a_c x)^2 - pi_0^2 x^2 with pi_0 = 2|g|. If a_c^2 <= pi_0^2, F is concave or linear in x and its minimum on [-1, 1] is at an endpoint. If a_c^2 > pi_0^2, F is convex with unconstrained minimum value -A0^2 pi_0^2 / (a_c^2 - pi_0^2) <= 0 at x* = -A0 a_c/(a_c^2 - pi_0^2). If x* is interior and A0 pi_0 != 0 then F(x*) < 0, contradicting stability. If pi_0 = 0 then x* is a zero of Abar, which changes sign there, contradicting Abar >= 0. If A0 = 0 then x* = 0 and Abar = a_c x changes sign at 0, again a contradiction (unless a_c = 0, where stability forces pi_0 = 0 and omega = 0 identically). So the minimum is at x = +-1. QED.

[numerical, lsw_minima.out] 30000 random references with a_s = 0, of which 5916 energetically stable: no interior global minimum and no interior soft point. Near h_F (J' = -1 + eps M, eps = 0.3, 0.1, 0.03, 3000 references each) every stable reference has its minimum at k = 0.

With a_s != 0, Lemma 3 of C154 stands: at onset the soft set is one momentum (the DM-shifted minimum), so again one species. Near h_F the bandwidth is 8 and the gap is of order eps, so the only low-energy region is around k = 0 (shifted by about eps D_n), matching (3.2).

### 4.2 Why two minima need a second harmonic

[derived] In a one-site cell with a nearest-neighbour bond, A_k and Pi_k contain only the first harmonic. Two minima at +-k_0 need either a second harmonic (cos 2k, which a next-nearest-neighbour coupling provides; in the frustrated ferromagnetic J1-J2 chain the k^2 coefficient changes sign at J2/|J1| = 1/4 and the minima split) or a k^4-stabilized Lifshitz point, i.e. a vanishing stiffness. Near h_F the stiffness is 1 - O(eps) and never vanishes. So the J1-J2 mechanism for multi-component and nematic liquids near saturation is absent. Lemma 3' and (3.2) make this precise for every reference state near h_F.

Answer to the spec's question in Task 2. [proved for uniform references, derived for textures] No: a one-site nearest-neighbour bond cannot have a stable polarized state whose magnon band has two degenerate minima at +-k_0 not related by the DM shift. There is therefore no dilute two-component gas to derive. C154 Lemma 3 covers this case, and its pair item is empty for stable references. For textured references the analogue is Lemma M2.

## 5. Task 3: bound states

### 5.1 Exact two-magnon problem of the number-conserving part

[proved, standard; numerical check twomag.out] In the frame of the moment n = e_z, split the bond into its part that conserves magnon number (XXZ about n with J_perp = (J'_uu + J'_vv)/2, a DM component D_n, J_z = J'_nn and the field B'_n) and the rest (pair creation g of (4.1), Delta N = +-2, and the cubic terms from J'_un, J'_nu, Delta N = +-1). For the conserving part, a magnon (down spin above the all-up state) has

    e(k) = c0 + 2|a| cos(k - phi),   c0 = -4 J_z - 2 B'_n,   |a| = 2 sqrt(J_perp^2 + D_n^2),   tan phi = D_n / J_perp (up to orientation),      (5.1)

and two magnons at total momentum K reduce, in the relative coordinate r >= 1 (hard core), to a half-line hopping problem with amplitude tau(K) = 2|a| cos(K/2 - phi) and a contact energy V = 4 J_z at r = 1. The ansatz f(r) = z^r gives z = -tau/V and

    E_b(K) = 2 c0 + V + tau(K)^2 / V,   bound iff |V| > |tau(K)|,   size xi = 1/ln|V/tau|.      (5.2)

The DM only shifts K by 2 phi (the Bethe 2-string with a twist). For h_F, E_b(K) = 2(1 - cos K) is the Heisenberg 2-string, with binding K^4/16 relative to the bottom of the continuum 16 sin^2(K/4), so the binding vanishes at K = 0 (the pair at K = 0 is an SU(2) descendant). twomag.py checks (5.2) against exact diagonalization of the two-magnon sector of 32-site rings for Heisenberg, easy-axis, easy-axis with DM and field, easy-plane, and easy-plane with DM. The agreement is 1e-4 to 3e-4 for bound states smaller than L/8, an exponentially small finite-size tail.

Near h_F, at the band bottom (K = 2 phi) the pair is bound iff |J_z| > |J'_perp| := sqrt(J_perp^2 + D_n^2) (easy axis in the twisted frame), with binding -4(|J_z| - |J'_perp|)^2/|J_z| = O(eps^2). [derived] The zero-energy scattering solution f(r) = r - a gives a = |J_z|/(|J_z| - |J'_perp|), so the dilute magnons form a 1D Bose gas with mass m = 1/(2|a|) and contact coupling

    g_eff = 8 |J'_perp| (|J'_perp| - |J_z|) / |J_z|,      (5.3)

repulsive on the easy-plane side (the Lieb-Liniger gas, one Luttinger liquid, which is XXZ with Delta_AF > -1), zero at SU(2) (free bosons), attractive on the easy-axis side.

### 5.2 Pairs and single magnons are never soft together at finite density

[proved for the number-conserving part, literature] On the attractive side the exact N-strings of the ferromagnetic XXZ chain (Delta = cosh gamma > 1) have energy E_N(K) = J sinh gamma (cosh N gamma - cos K)/sinh N gamma. Per magnon at K = 0 this is J sinh gamma tanh(N gamma/2)/N, strictly decreasing in N. So the lowest energy per flipped spin belongs to the infinite cluster (a reversed domain). In the weakly attractive continuum limit this is McGuire's result E_N = -m g^2 N(N^2 - 1)/24 for the attractive delta Bose gas. Pairs never form a stable fluid. When the field makes magnons favourable the system jumps to the reversed or canted state: a first-order transition, c = 0 on both sides. This is the classical spin flop of section 2.2 (F and P).

[derived] On the repulsive side there is no bound state near the band bottom (|tau| > |V| at K = 2 phi), so the gas has one species, a single Luttinger liquid with K of order 1/sqrt(gamma_LL), large. This is the known field-induced XXZ result near saturation: one species and a z = 2 onset. It is the limit theta_0 -> 0 of the cone of section 2.3.

[derived] The terms that break magnon number near h_F are of order eps, while the pair binding is of order eps^2 at the band bottom. So outside the wedge of section 2.3, where U(1) about n holds only up to relative order eps, the two-magnon problem is not the right description: bound states are the classical breathers of the Landau-Lifshitz equation (1.4) and are part of the classical analysis of section 3. The pair-creation term g a_k a_{-k} is a direct source for the pair field at K = 0. So the pair channel has no independent order parameter. A "pair fluid plus magnon fluid" would require a classical ground state with two independent soft canonical pairs, excluded by 3.5. Inside the wedge, the bound states belong to the single boson (2.1).

Answer to Task 3. [derived] No regime near h_F has bound pairs and single magnons simultaneously soft at finite density. The attractive side collapses (first order), the repulsive side has one species, and the U(1)-breaking terms source the pair channel. The comparison case (field-induced XXZ near saturation, one species) is reproduced. The difference from the frustrated J1-J2 chain, where a nematic pair liquid and a c = 2 two-fluid phase exist, is that there the next-nearest-neighbour coupling makes N >= 3 clusters unbound while pairs stay bound. With a nearest-neighbour qubit bond the cluster energy per particle is monotone (5.2).

## 6. From classical to quantum

[assumed, the semiclassical correspondence (SC)] As hbar_eff -> 0 uniformly on the neighbourhood, the quantum infrared of (1.4) is fixed by the classical ground state and its soft modes. Each linear canonical pair (phason, U(1) Goldstone) becomes one compact boson with K of order 1/hbar_eff (c = 1). A real linear mode at a symmetric or tuned pitchfork becomes an Ising critical point (c = 1/2). A quadratic pair stays non-conformal (z = 2). Gapped classical modes stay gapped. This is the standard large-spin (here large spin per coherence length) correspondence. It reproduces K = pi/(2 sqrt(2 eta)) for XXZ at Delta_AF = -1 + eta, and the dilute Bose gas at saturation.

[derived, with SC] Quantum corrections are of relative order hbar_eff. They can close a classical gap only within O(hbar_eff) of a classical instability, and there the soft content is that of the instability (section 3.5). Non-perturbative effects are of order exp(-c/hbar_eff). They split degenerate classical vacua (tunnelling), which can only produce an Ising transition at the end of a Z_2 coexistence, and they cannot create a second boson. Order-by-disorder lifts accidental classical degeneracies (the accidental circles of 2.1) and leaves at most one mode.

## 7. Numerics

- [numerical] lsw_minima.out, twomag.out, csl_tilt.out, csl_fine.out, csl_firstorder.out (local, 1 to 13 s each), described above and in PROGRESS.md.
- [numerical] K-stratum CSL family (the most dangerous kinematic candidate for P1, a planar soliton lattice that could tilt spontaneously, giving floating + SSB with K large). Continuum E = int [ |n'|^2/2 - (n x n')_z + (kappa/2) n_z^2 - h n_x ], lattice bond v = (-1, -1, -1 + kappa eps^2/2, 0, 0, -eps, -h eps^2, 0, 0). The CSL period is T = 6.29, 6.36, 6.61, 7.22, 8.29 at h = 0.05, 0.15, 0.3, 0.45, 0.55 (diverging toward the Dzyaloshinskii value pi^2/16 = 0.617). The tilt instability is at kappa_c = -0.999, -0.989, -0.955, -0.898, -0.846. A first-order transition to the canted uniform state (n_z = +-0.86 to +-0.99, Theta R_pi(z) broken, c = 0) always comes first, at kappa_1 = kappa_c + 0.0000, 0.0002, 0.0031, 0.0194, 0.0571 (roughly h^4/2). The two-component minimization at kappa_c +- 0.05, +- 0.15 and the fine scan at h = 0.45 find no tilted soliton lattice below both. So in this family the P1 neighbour (floating + SSB) does not exist at all, and the CSL ends first order. The limit h -> 0 is the rotating-frame Heisenberg point, where all three states are degenerate (the U(1) class at Delta = -1).
- [numerical, pilot, inconclusive] texture_periodic.py and texture_refine.py: lattice ground textures of random near-h_F bonds at N = 400 to 600. Most bulk spectra are harmonics of one wave vector, but for long periods (100 to 250 sites) the bulk window and the local minima found by L-BFGS do not allow a clean test. Not used as evidence. The N = 2000 version is in the Typhon job.
- [pending] Typhon job 17599583 (one node, 48 cores, 8 h, PENDING at the end of this run). Stage D: open-chain DMRG of the CSL family at eps = 0.5 (period about 13 to 17 sites), points (h, kappa) = (0.3; -0.40, -0.70, -0.85, -0.95, -1.05, -1.30), (0; -0.60) as a U(1) calibration, (0.45; -0.60, -0.85), at L = 64, 128, 192, 256 with chi up to 160. It probes lengths up to 256 sites, i.e. 15 to 20 soliton periods, and should read c near 1 on the CSL side and near 0 on the canted side, with no window above 1.2. The crossover length of a liquid with K of order 1/hbar_eff is long, so the expected finite-size signature is a slow approach to 1 from below, as at XXZ near Delta = -1. Stage X: 48 random bonds of the pure-DM wedge, N = 2000, 16 starts, spectral test of single periodicity.

## 8. Theorem-style summary

Notation. U_eps0 = { (1.1) with 0 < eps < eps0 }, a punctured neighbourhood of h_F in B_2. W_C and W_Z are the wedges of relative width O(eps) around the circle locus C and the pure-DM sphere Z (sections 2.3, 2.4).

- Lemma 1 [proved]. Coordinates (1.1), directions S^10/SO(3) of dimension 7.
- Lemma 2 [derived]. The continuum functional (1.2), its gauge form (1.3) and the semiclassical parameter hbar_eff = 2 eps^{1/2} (generic) or 2 eps (pure-DM wedge), uniformly small on U_eps0.
- Lemma 3' [proved]. Refinement of C154 Lemma 3: with no DM along the moment, a stable polarized state has its magnon minima at k in {0, pi} only. No two-minimum dispersion exists for a one-site nearest-neighbour qubit bond.
- Lemma M1 [proved, given Mather's graph theorem]. Classical ground states of the continuum functional are uniform or singly periodic.
- Lemma M2 [derived; Jacobi/Morse, Poincare, Contreras-Iturriaga]. A periodic ground texture has no Hessian zero mode at nonzero Bloch momentum.
- Lemma M3 [derived]. An extra zero mode at Bloch momentum 0 is the canonical partner of the phason (z = 2), or the state is not a minimizer.
- Lemma 4 [proved for the number-conserving part (Bethe strings), derived for the full bond]. Two-magnon bound states (5.2) and the effective coupling (5.3). The attractive side collapses (first order), the repulsive side has one species, pair channels are sourced by pair creation. No two-fluid regime.
- Proposition 5 [derived, from M1 to M3, 3.4 and Lemma 3']. Every classical ground state in U_eps0 has at most one gapless canonical pair. Linear pairs (c = 1) occur only in W_C and W_Z, real linear modes (c = 1/2) only on the pitchfork loci I (codimension 2, codimension 1 in K) and on double-sine-Gordon lines inside W_C, quadratic pairs (z = 2) on Pokrovsky-Talapov, saturation-onset and compressibility surfaces.
- Theorem [derived, conditional on (SC) of section 6]. There is eps0 > 0 such that every conformal bond of U_eps0 has c in {1/2, 1}. In particular c <= 1 on a full punctured neighbourhood of h_F, every attractive (K > 1) liquid there is a single boson, and the assumption (A6_TK) of C159 holds on U_eps0. This removes the window that C159 Proposition C left open around h_F, in T, in K and in N alike (the argument uses no symmetry).

Narrowest remaining assumption. (SC), the standard semiclassical correspondence at large spin per coherence length, plus the use of Mather's graph theorem and the non-conjugacy of minimizers for the continuum functional (lattice corrections handled by persistence, heuristic). The size eps0 is not quantified. The Typhon DMRG at eps = 0.5 tests the finite-eps end in the most dangerous family.

## 9. Consistency checks

- XXZ near Delta_AF = -1 (pure-DM direction, and the axial easy-plane case): c = 1 with K proportional to 1/hbar_eff (1.2, 1.5).
- Heisenberg two-magnon problem: E_b(K) = 2(1 - cos K) and binding of order K^4 reproduced by (5.2) and by ED.
- Field-induced XXZ near saturation: one species (5.2), z = 2 onset (3.4 with alpha = beta = 0).
- Easy-axis XXZ ferromagnet in a longitudinal field: first-order reversal, from the monotone string energies (5.2).
- J1-J2 ferromagnetic chain: the two-minimum and nematic mechanisms need the second harmonic, absent here (4.2).
- Dzyaloshinskii's CSL: h_c = pi^2/16 for unit twist, consistent with the diverging T(h) of csl_tilt.out.
- C159 Lemma B: no Ising point in T near h_F (2.2), since b = 0 puts T in the hard case with the Theta-broken pair +-e_1.
- C159 Proposition C: its attractive windows (x > 1 for two bosons, x > 15/8 for Ising x boson) are exactly the large-K liquids of W_C and W_Z, where 3.5 leaves one boson.
- C154 Lemma 3: covered and sharpened (Lemma 3').

## 10. What is not done

- (SC) is assumed, not derived from the lattice. A rigorous version would need a uniform large-spin expansion of the Landau-Lifshitz path integral in 1+1 dimensions.
- The non-conjugacy input for M2 and the Mather/Poincare-Bendixson step are for the continuum functional. For the lattice (a 4-dimensional symplectic map) persistence of hyperbolic minimizing circles is a heuristic step.
- The phase diagram of (2.2) is not drawn in full. Only its soft-mode content is controlled.
- The Typhon job has not run. Its texture stage would test M1 on the lattice, its DMRG stage the finite-eps end of the CSL family.
