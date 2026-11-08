# C162 derivations: two-valley soliton onsets are immiscible

This run was done by a Claude Opus agent (2026-11-08), not by Codex. Astra was not used (out of credits).

Tags as in C160 and C161: [proved], [derived], [numerical], [heuristic], [assumed]. A tag applies to the paragraph or displayed claim it opens. Notation is that of codex/work_C160/PROOF.md (blow-up at h_F) and codex/work_C161/REFEREE.md (section 2b). Scripts and console captures are in this folder.

## 0. Setting and scales

The pure-DM wedge of C160 (eq. 2.2), in the K stratum (D along the Z2 axis z, field in the plane), is

    E[n] = int dx [ |n'|^2/2 - (n x n')_z + W(n) ],   W(n) = W(n_x, n_y, n_z^2),        (0.1)

with W quadratic plus linear and even in n_z. The Z2 is n_z -> -n_z (on spins Theta R_pi(z)). Lattice length l = 1/eps, lattice energy eps E, and the quantum theory is semiclassical with hbar_eff = 1/(s l) = 2 eps (C160 1.3).

[derived, C161 2b] A soliton is a 2 pi winding of the in-plane angle over a Z2-symmetric background n_bg (n_z = 0). Its collective coordinates are the position X and the oriented solid angle Omega, which are canonically conjugate with [X, Omega] = i hbar_eff in rescaled units. The lattice crystal momentum is P = s Omega (Haldane). So the soliton band is E(P) = min{ E[n] : Omega[n] = P/s }. A tilt-unstable soliton has two degenerate minima P_pm = pi -+ delta, exchanged by Z2 (which is antiunitary and sends P -> -P).

Scales [derived]. In lattice units the soliton energy is eps e(Omega) with e of order one, so the valley mass is M = 1/E''(P) = 1/(4 eps e''(Omega_pm)), of order 1/eps. The interaction potential between two solitons is eps U(R/l) with U of order one and range of order l. The two-body problem is therefore controlled by k2 = M_resc/hbar_eff^2, of order eps^{-2} in rescaled units, which is large. Every two-body statement below is about the regime k2 -> infinity, with eps small but finite.

## 1. Effective theory of the dilute two-valley soliton gas

### 1.1 Tails and channels

[proved, algebra] Linearize (0.1) around n_bg with tangent coordinates u (in-plane) and w (out-of-plane, n_z). Because D = z is perpendicular to n_bg, the DM term has no quadratic part around n_bg (the linear part u' is a total derivative). By Z2, W has no u w cross term at n_bg. So the static tail equations decouple,

    -u'' + m_u^2 u = 0,   -w'' + m_w^2 w = 0,   m_u^2 = d^2_phi W(n_bg),   m_w^2 = d^2_theta W(n_bg),      (1.1)

and a soliton of valley v has tails u ~ r_u^v e^{-m_u x}, w ~ r_w^v e^{-m_w x} on the right and l_u^v e^{m_u x}, l_w^v e^{m_w x} on the left. Z2 gives r_u^- = r_u^+, l_u^- = l_u^+, r_w^- = -r_w^+, l_w^- = -l_w^+.

[derived; input: Manton's tail formula] Two well separated solitons, valley i on the left and j on the right at distance R, interact through

    V_ij(R) = - sum_k c_k r_k^i l_k^j e^{-m_k R} + (faster terms),   c_k = 2 m_k > 0,      (1.2)

the standard "like charges attract" rule of a massive scalar. For the in-plane channel the 2 pi winding forces r_u l_u < 0 (the angle leaves the background upward and returns from below), so the u channel is repulsive and valley independent. The w channel is valley odd.

### 1.2 Two-body problem and the hard-core gas

[derived] Near the valley minima each soliton has kinetic energy (P - P_v)^2/(2M), equal masses by Z2. The relative motion of a pair (i, j) at zero energy obeys psi'' = k2 V_ij psi. Two 2 pi solitons can pass through each other only by tunnelling through a 4 pi core, and two opposite-valley solitons can exchange valleys (both tilts flip, total momentum conserved) only by tunnelling over the tilt barriers. Both amplitudes are of order exp(-S/hbar_eff) with S of order one [heuristic for the value of S, derived for the exponential form]. So at low density the gas is a 1D mixture of hard-core particles whose order along the line is fixed up to exponentially small exchange.

[derived; inputs: strong-coupling mapping of 1D mixtures (Ogata-Shiba, Matveev, Volosniev et al. 2014, Deuretzbacher et al. 2014)] For hard-core particles with fixed order the ground-state energy at density n is the Tonks value with the length reduced by the excluded lengths, E = (pi^2 hbar^2 N^3)/(6 M (L - sum_bonds a_b)^2). So each neighbouring pair contributes c a_{s s'} with c = pi^2 hbar^2 n^3/(3M) > 0 and a_{ss'} the zero-energy scattering length (the zero of the outer solution psi ~ R - a) of the ordered pair (s left of s'). Writing the valley as a pseudospin S^z = s/2 and adding the exchange gives the squeezed-lattice chain

    H_ps = sum_k [ J_perp (S^x_k S^x_{k+1} + S^y_k S^y_{k+1}) + J_z S^z_k S^z_{k+1} ] - h_ps sum_k S^z_k,
    J_z = c (a_{++} + a_{--} - a_{+-} - a_{-+}),   |J_perp| ~ c e^{-S/hbar_eff},      (1.3)

with h_ps = 0 in K by Z2. The S^z_tot conservation is emergent: a single valley flip changes the crystal momentum by 2 delta, which is incommensurate, so it is oscillating and irrelevant away from isolated commensurate densities. A nearest-neighbour chain can be sublattice-rotated, so the sign of J_perp is immaterial. With Delta = J_z/|J_perp| the phases at h_ps = 0 are
- Delta < -1: Ising ferromagnet, all solitons in one valley. Immiscible, Z2 spontaneously broken (uniform n_z), c = 1 (the charge boson only).
- -1 < Delta <= 1: pseudospin Luttinger liquid, c = 2. This is the miscible two-component gas.
- Delta > 1: pseudospin Neel, valleys alternate along the line, pseudospin gapped, c = 1. With h_ps != 0 (possible only without Z2) a partially polarized window of width of order |J_perp| near h_ps = J_z reopens c = 2.

Check against the delta-interacting two-component Bose gas [derived]. With g_11 = g_22 = g and g_12 = g', one has a_11 = -2/(Mg), a_12^even = -2/(Mg'), a_12^odd = 0, which gives Delta = 1 - 2g'/g. So g' = g gives the SU(2) ferromagnet (Eisenberg-Lieb), g -> infinity gives the antiferromagnetic Heisenberg chain (the Gaudin-Yang equivalent), and g' > g gives the immiscible Ising ferromagnet, as in the literature.

### 1.3 Statistics and Eisenberg-Lieb

[heuristic] The 2 pi twist operator prod_j exp(i theta_j S^z_j), with theta_j ramping from 0 to 2 pi, is local up to a c-number, because exp(2 pi i S^z) = -1 for s = 1/2. So the solitons are created by local bosonic operators, and the Haldane momentum pi of the planar soliton is this -1 per site. For bosons the like-valley pair is in the spatially even channel, and the opposite-valley pair is in the even channel for the pseudospin-symmetric combination. The exact borderline then reads [derived]: Delta < -1 if and only if a_{++} < a_{+-}^even (for a_{+-}^even <= a_{+-}^odd, which holds for repulsive potentials because the odd channel carries an extra Dirichlet condition). For fermionic statistics the comparison would be with odd-channel data instead. In the semiclassical regime even and odd scattering lengths differ only by exp(-S/hbar_eff), so statistics affects nothing below.

Eisenberg-Lieb (bosons, valley-independent interactions) puts the gas exactly at Delta = -1, the SU(2) ferromagnet. Three things break it here. The interactions are valley dependent (V_{++} != V_{+-}). There is no SU(2) relating the valleys, only Z2 times an emergent U(1). And valley exchange occurs by flip-flop as well as by passage. Section 3 shows that the breaking always goes in the ferromagnetic (Ising, immiscible) direction in K, and by a margin that is a power of eps while |J_perp| is exponentially small.

## 2. Scattering-length calculus

Lemma S1 (monotonicity) [proved]. Let V_lambda be a one-parameter family of potentials on R > 0 with a hard core at R = 0 (or a repulsion that makes the zero-energy solution vanish there), and psi_lambda the zero-energy solution normalized to R - a_lambda at infinity. Then

    d a_lambda / d lambda = k2 int_0^infinity (d V_lambda/d lambda) psi_lambda^2 dR.      (2.1)

Proof. chi = d psi/d lambda solves chi'' = k2 V chi + k2 (dV) psi. The Wronskian (psi chi' - chi psi')' = k2 (dV) psi^2. At R = 0 both functions vanish. At infinity psi = R - a and chi = -da/dlambda, so the Wronskian tends to da/dlambda. QED. So V_1 <= V_2 pointwise, with no bound state in V_1, implies a_1 <= a_2, strictly if V_1 < V_2 on a set where psi does not vanish. [numerical, scatlen.out part 3] finite differences agree in sign.

Lemma S2 (one channel is additive) [proved]. For V = g e^{-lam R}, psi'' = k2 g e^{-lam R} psi is solved by psi = K_0(z), z = (2/lam) sqrt(k2 g) e^{-lam R/2}, and

    a = (2/lam) [ ln( sqrt(k2 g)/lam ) + gamma_E ].      (2.2)

[numerical, scatlen.out part 1] the solver reproduces (2.2) to 1e-6 at k2 = 1e2 to 1e6. So for V_ij = g r^i l^j e^{-lam R} one gets a_ij = const + (ln|r^i| + ln|l^j|)/lam. This is a mixture of hard rods with additive lengths, and J_z of (1.3) vanishes identically. The pseudospin coupling therefore never comes from the slowest channel alone.

Lemma S3 (sign rule) [derived]. Let the slowest channel (rate lam_1) be repulsive and let a second channel (rate lam_2 > lam_1) contribute -c_2 x^i y^j e^{-lam_2 R}, with x, y its right and left tail amplitudes. To first order in the second channel, which is exact as k2 -> infinity because psi lives where e^{-lam_2 R}/e^{-lam_1 R} -> 0, Lemma S2 shifts psi_ij by (ln|r^i l^j|)/lam_1 and

    J_z / c = a_{++} + a_{--} - a_{+-} - a_{-+} = - c_2 I (X^+ - X^-)(Y^+ - Y^-),
    X^i = x^i |r^i|^{-lam_2/lam_1},  Y^j = y^j |l^j|^{-lam_2/lam_1},  I = k2 int e^{-lam_2 R} psi_0^2 dR > 0.      (2.3)

So J_z scales as k2^{-(lam_2/lam_1 - 1)}, a power of eps, while |J_perp| is exponentially small. Hence Delta = J_z/|J_perp| -> -infinity or +infinity according to the sign of (2.3), and the c = 2 window |Delta| < 1 is reached only where (2.3) vanishes. [numerical, scatlen.out part 2] the sign of (2.3) and the scaling k2^{-1.5} (lam_2/lam_1 = 2.5) are reproduced in three cases with a common hard core added.

## 3. The K stratum

### 3.1 Rearrangement lemma

Lemma K1 [proved, continuum]. For (0.1), let n be any minimizer of E in a class defined by the in-plane data (boundary conditions at infinity, winding of phi = arg(n_x + i n_y), periodicity, or constraints on phi at fixed points). Then n_z does not change sign. If n_z is not identically zero, it has no zeros, and its tail coefficients on both sides are nonzero and of the same sign.

Proof. Replace n_z by |n_z|. The gradient term is unchanged (|d|w|| = |dw| almost everywhere), W is even in n_z, and the DM density (n x n')_z = n_x n_y' - n_y n_x' does not involve n_z. Constraints on phi are untouched. So |n_z| is also a minimizer, hence a smooth solution of the Euler-Lagrange equation. Constraint forces act on phi only and are in-plane at planar configurations. At a zero x_0 of n_z the function |n_z| has a minimum, so n_z(x_0) = n_z'(x_0) = 0. The equator with tangent velocities is invariant under the Euler-Lagrange flow (Z2 fixed set), so by uniqueness n_z vanishes identically. Hence n_z has no zeros if it is not zero. In the tails w solves -w'' + (m_w^2 + q(x)) w = 0 with q decaying, whose decaying solution is unique up to scale and is asymptotic to a nonzero multiple of e^{-m_w |x|}. Since w > 0 (say), both coefficients are positive. QED.

Lattice form [proved]. The classical lattice energy is sum_j [n_j^T J n_{j+1} + B . n_j] with J = diag(J_x, J_y, J_z) + D_z A(z) and B in the plane. n_z enters only through J_z n_z^j n_z^{j+1}. For J_z <= 0 (near h_F, J_z = -1 + O(eps^2)) the replacement n_z^j -> |n_z^j| does not raise the energy and strictly lowers it across any bond with a sign change. At a site with n_z^j = 0 stationarity forces n_z^{j-1} + n_z^{j+1} = 0, so |n_z| = 0 propagates. The conclusion is the same.

Corollary K2 (classical inter-valley potential is always higher) [proved]. Define V_{+-}(R) as the energy (relative to two isolated solitons) of the local minimizer with solitons at fixed separation R (phi = pi at the two centres) and n_z > 0 near the first, n_z < 0 near the second, whenever it exists, and V_{++}(R) as the minimum in the same constraint class with one sign. Then

    V_{++}(R) < V_{+-}(R)   for every R.      (3.1)

The |n_z| image of the (+-) minimizer has the same energy and lies in the (++) class, and it is not smooth, so it is not stationary and can be lowered. The same argument shows that in every sector (2 pi N winding, periodic textures of any period) the classical minimizer has all solitons in one valley.

[numerical, soliton_pair.out of C161 and lemma_check.out] (3.1) holds at T = 10 to 22 in C161's family and in a family with the field rotated by psi = 0.3 off the principal axes, where the extra classical symmetry x -> -x, n_y -> -n_y is broken and the tilted soliton is strongly asymmetric (fitted tail coefficients B_R = 0.36, B_L = 1.03, decay rate 0.2785 = m_w, min n_z/max n_z = 5e-4).

### 3.2 Consequence for the dilute gas

Proposition K3 (immiscibility) [derived; inputs K2, S1, the hard-core mapping 1.2]. At a soliton onset in K whose two valleys form a dilute gas (no two-body bound state in any channel), J_z < 0, and

    for bosons:  a_{++} < a_{+-}^even, hence Delta < -1 at every eps;
    for any statistics, semiclassically:  J_z = -(positive power of eps) c, |J_perp| = c e^{-S/hbar_eff}, Delta -> -infinity.      (3.2)

So one valley wins, Z2 breaks spontaneously, and c = 1. The first line uses S1 along V_lambda = V_{++} + lambda (V_{+-} - V_{++}), which is pointwise increasing by (3.1) and has no bound state if V_{++} has none. The second line uses S3. In the tail language, with the u channel valley even, (2.3) gives J_z/c = -c_w I (2 X_w)(2 Y_w) with X_w Y_w > 0 by K1. Up to the flip-flop amplitude, which shifts the even-channel potential by a term of order exp(-S/hbar_eff), this is the classical Eisenberg-Lieb mechanism, strengthened from Delta = -1 to Delta < -1.

### 3.3 Which tail is slow, and why the onset is first order

[proved, algebra] For a planar soliton the DM energy is -2 pi (topological), so its shape does not depend on the DM. It obeys the first integral phi'^2/2 = V_eq(phi) - V_eq(phi_bg) with V_eq(phi) = W on the equator, and its gradient-plus-potential energy scales as sqrt(s) under W -> sW. So the planar onset is at s_c = (2 pi/E_gp(1))^2 exactly. Writing W = Q(phi) cos^2 theta - L(phi) cos theta on S^2 (gauge S_zz = 0, Q quadratic, L = B . n) one finds the exact operator identity

    H_perp = H_par + (m_w^2 - m_u^2) + 2 phi'(x) + 2 [ L(phi(x)) - L(phi_bg) ],      (3.3)

where H_par = -d^2 + V_eq''(phi) is the in-plane Hessian (ground state phi', eigenvalue 0) and H_perp = -d^2 + U, U = phi'(2 - phi') + d^2_theta W, is the tilt Hessian. Hence

    lam0(H_perp) >= m_w^2 - m_u^2 + min_x [ 2 phi' + 2 (L(phi) - L_bg) ].      (3.4)

[numerical, bound_scan.out] Over 4719 random admissible K shapes (background the Z2-symmetric global minimum on S^2) the bracket in (3.4) was negative in only 2 shapes, and lam0 - (m_w^2 - m_u^2) >= +0.088 in all of them. 21 shapes were tilt-unstable at onset, every one with m_u^2/m_w^2 >= 4.6. Targeted searches (opt_caseB.py, probe_caseB.py) under the constraint m_w >= m_u found lam0 >= +0.62. [derived from this] At a K two-valley onset the tilt tail is the slow one, m_w < m_u. By K1 its coefficients have one sign, so by (1.2) like valleys attract at long range with strength of order B^2, which exceeds the zero-point scale (of order eps^2 in rescaled units) unless the tilt is within O(eps) of the bifurcation. So like-valley pairs bind, and by K1 every minimizing cluster is single valley with same-sign tilt tails, so clusters attract clusters at long range, E_N/N decreases, and the onset is first order into the same-tilt soliton lattice. That is C161's finding (lemma_check.out reproduces it off-axis: U_same(10) = -0.0088 with E_sol = +0.0116). Proposition K3 covers the remaining (non-generic) continuous case. Either way c <= 1.

Remark (scaling direction) [derived]. Under W -> sW the relative DM weight in H_perp is 2/sqrt(s), so tilt instability strengthens as s grows, while the soliton energy increases with s (W >= W_bg on the minimizer). So "tilt-unstable at s_c" is the correct onset criterion, up to subcritical tilted branches, to which K1 to K3 apply unchanged.

## 4. The other protected case: U = I R_pi(n_bg)

[derived] A symmetry protects two degenerate soliton valleys over a fixed background if it fixes n_bg, maps the 2 pi winding sector to itself and reverses P. Antiunitary candidates are Theta R with R n_bg = -n_bg, i.e. Theta R_pi(a) with a perpendicular to n_bg. Classically that is the reflection through the plane perpendicular to a, which is the K stratum (D along a, B perpendicular to a). Unitary candidates must contain site inversion I (which reverses P). I alone maps a soliton to an antisoliton. I R_pi(n_bg) maps the winding sector to itself (inversion and the pi rotation both reverse the in-plane angle). Invariance under U = I R_pi(x) requires S block diagonal (x | yz), D in the yz plane, B along x. After a rotation in the yz plane this is a 5-dimensional class, not contained in K unless D lies on a principal axis of the yz block. Theta I R maps P to P and protects nothing. So K and U are the only strata with symmetry-protected two-valley soliton bands.

[derived, tail level] In U, D is perpendicular to n_bg = x, so the tails decouple into the two principal directions of the yz block, both reversed by R_pi(x). U maps the right tail of valley + to minus the left tail of valley -, x^- = -y^+, y^- = -x^+ in every channel. The slowest channel is additive (S2). The next one gives, by (2.3),

    J_z / c = - c_2 I (X^+ + Y^+)^2 <= 0,      (4.1)

so the U class is immiscible at onset too. [numerical, scatlen.out] confirmed for a U-like amplitude set. There is no analogue of K1 here (W is not even in the tilt coordinate n . D-hat when D is off the principal axes), so (4.1) is a statement about the dilute, semiclassical limit, not about all separations. If the slow channel is self-attractive the onset is first order, as in K.

## 5. N and T

[derived] T (B = 0, Kramers Theta only) has no Theta-fixed background (Theta sends n -> -n), and N has no symmetry. So a two-valley soliton band over a fixed background is accidental there, codimension 1. The dilute theory (1.3) still applies, now with h_ps from the valley splitting and with the sign rule (2.3) unconstrained by symmetry. Two situations differ.
- Near K or U the strict negativity of J_z persists. A Z2-breaking parameter zeta changes X^+ - X^- only at O(zeta^2) (Z2 sends zeta -> -zeta and exchanges the valleys), so J_z stays negative for |zeta| below an eps-independent bound, and the immiscible Ising ferromagnet is fully polarized for any h_ps. No c = 2 there.
- Far from K and U nothing fixes the sign. scatlen.out exhibits the miscible sign for a self-repulsive, valley-dependent fast channel. In N the valley masses also differ, and mass imbalance makes the hard-core energy sequence dependent already at order n^2 (second order in the imbalance), with a sign I have not determined [heuristic: it is expected to favour segregation]. So an accidental N or T two-valley onset with non-binding interactions, monotone tails and the miscible sign is not excluded. If it existed, the c = 2 region would be an open but exponentially thin sliver (|Delta| < 1, or the partially polarized window near h_ps = J_z when Delta > 1). Oscillating tails, which occur in N whenever |D . n_bg| > |m_u - m_w|/2, bind solitons semiclassically and make such onsets first order [heuristic].

## 6. Beyond the dilute limit

[proved, classical] By K1 every minimizing periodic texture in K, with any number of solitons per period, has single-signed n_z. So alternating (+-) and mixed soliton lattices are never classical ground states, and the classical ground state at every density is uniform, planar CSL or single-valley tilted CSL.

[derived] A one-valley lattice with a gapless second band would need a soft valley flip. A single flipped soliton costs a finite energy, strictly positive by (3.1) at every lattice spacing, and a Z2 domain wall costs a finite energy too. Both are gapped classically and exceed every quantum correction (O(hbar_eff) shifts and exponentially small exchange). The C161 Bloch spectrum of the tilted CSL (lambda_1 >= 0.049) is this gap in linear form. M2 (no Bloch zero mode at k != 0) and M3 (a k = 0 zero mode pairs with the phason) of C160 see the classical side. What they cannot see is the pseudospin sector, whose couplings are the scattering-length differences of section 2 and the exponentially small exchange. That sector is controlled by (3.1), which holds uniformly in R, so there is no intermediate density at which the sign changes. At the dilute end the sector is the Ising ferromagnet of K3, at the dense end it is the classical flip gap, and the two match. Hence c = 1 (phason) plus a spontaneously broken uniform Z2 throughout the two-valley region of K.

## 7. Theorem-style summary

- Lemma K1 [proved, continuum and lattice for J_z <= 0]: in K every minimizer in any in-plane constraint class has single-signed n_z, with nonzero same-sign tilt tails.
- Corollary K2 [proved, classical]: V_{++}(R) < V_{+-}(R) at every separation, and every N-soliton or periodic minimizer is single valley.
- Lemma S1 [proved]: da/dlambda = k2 int psi^2 dV/dlambda. Lemma S2 [proved]: single-exponential scattering lengths are additive (no pseudospin coupling). Lemma S3 [derived, k2 -> infinity]: J_z/c = -c_2 I (X^+ - X^-)(Y^+ - Y^-).
- Mapping (1.3) [derived; standard strong-coupling theory of 1D mixtures, exchange exponentially small in 1/hbar_eff, heuristic for its action].
- Soliton onset lemma, K stratum [derived, given (1.3)]: at a two-valley soliton onset over a Z2-symmetric background the gas is immiscible (Delta < -1 for bosons at any eps, Delta -> -infinity semiclassically for any statistics). The phase above onset is one valley with uniform Z2 breaking, c = 1.
- Tail ordering [proved inequality (3.4), numerical for its bracket]: tilt instability at onset occurs only with m_w < m_u, so like valleys attract at long range and the onset is first order into the same-tilt CSL (c = 1), as in C161.
- U class [derived, tail level]: immiscible by (4.1).
- Beyond the dilute limit [proved classically, derived quantum]: no two-valley mixed phase at any density in K.
- Remaining assumption (A_sol) [assumed]: no accidental (unprotected, codimension 1) two-valley soliton onset in N or T far from K and U has non-binding interactions together with the miscible sign of (2.3) (including the mass-imbalance correction of section 5).

## 8. Consistency checks

- C161's numbers are reproduced by kfun.py: lam0 = -0.0670 at kappa = 0.05, mu = 2.3, h = 0.0238 (C161 -0.067), and the off-axis ring energies (lemma_check.out) match C161's soliton_pair.out to about 1 percent at psi = 0.3.
- Delta-interaction limit: (1.3) gives Delta = 1 - 2 g'/g, so Eisenberg-Lieb at g' = g, Gaudin-Yang (AF Heisenberg) at g = infinity, immiscible for g' > g.
- K_0 formula (2.2) checked to 1e-6, and the scaling J_z ~ k2^{-(lam_2/lam_1 - 1)} reproduced (factor 1e3 between k2 = 1e3 and 1e5 for lam_2/lam_1 = 2.5).
- J1-J2 literature [literature, from memory, not rechecked in this run]: miscible two-component (c = 2) regions near saturation are reported in the antiferromagnetic J1-J2 chain, and the immiscible vector-chiral phase in other ranges (Kolezhuk, Heidrich-Meisner, Greschner and Vekua, PRB 85, 064420 (2012), and earlier work by Okunishi and by Hikihara et al.). This does not contradict K3. Magnons are point-like quantum particles whose valley dependence comes from the momentum dependence of a lattice-scale vertex. No rearrangement inequality like (3.1) exists for them. K3 uses the extended semiclassical structure of solitons (tails, the |n_z| replacement) that magnons lack.
- The tilted soliton's tail decay rate in lemma_check.out is 0.2785 = sqrt(h + kappa) at h = 0.0289, as (1.1) predicts.

## 9. What is not done

- The exchange action S (pass-through and flip-flop) is not computed. Only its exponential form is used.
- The bracket in (3.4) is not proved non-negative. Tail ordering at onset is numerical over random shapes.
- U class: no rearrangement lemma valid at all separations, and no numerics on an actual U bond.
- N and T accidental two-valley onsets: not searched. The sign of the mass-imbalance contribution is not derived.
- No DMRG. With no miscible region in K or U, Task 4 does not apply, and the onset energy scales (exponentially small pseudospin gaps) are beyond lattice numerics at the eps where the blow-up is valid. No Typhon job was submitted.

## 10. Construction side (coordinator's priority shift, 2026-11-08)

The coordinator asked to weight the search for a miscible region as heavily as the no-go, with four pointers. Here is what each one gives within the theory of sections 1 and 2.

(1) Accidental degeneracy in N suffices [agreed, derived]. A codimension-1 c = 2 locus would refute c <= 1, and in N the two band minima can be made degenerate by one tuning. Section 5 already allows this. What decides it is the sign of J_z in (1.3) on that locus. Two refinements matter for the construction. (a) Near K or U the sign is fixed: J_z < 0 with a margin that does not shrink with eps (X^+ - X^- changes only at O(zeta^2)), so the c = 2 locus, if it exists, is a finite distance in the blow-up parameters from every protected stratum. (b) In K every two-valley onset found was Case A (tilt tail slow, like valleys bind, section 3.3). So the N candidate cannot be obtained by deforming the K mechanism. It needs a different way of making two soliton valleys, one in which the slow channel is self-repulsive in both valleys and the fast channel is valley dependent with the miscible sign of (2.3).

(2) Oscillatory tails [derived for the sign change, heuristic for binding]. With D . n_bg = d != 0 the static tail equations are -u'' + 2 d w' + m_u^2 u = 0, -w'' - 2 d u' + m_w^2 w = 0, with decay rates lambda solving (m_u^2 - lambda^2)(m_w^2 - lambda^2) + 4 d^2 lambda^2 = 0. They are complex (oscillating tails) for |m_u - m_w| < 2|d| < m_u + m_w, an open condition in N, and the upper end is the magnon-softening (Lifshitz) line of the background. Then V_ij(R) contains Re[A_ij e^{-(kappa + i q) R}], and the valley-odd part of J_z in (2.3) oscillates with R, so its sign can indeed flip as the parameters move the zero-energy wavefunction across a node. The obstruction is that the valley-even part oscillates too. Each V_ii then has attractive wells of depth of order |A| e^{-kappa R}, and the bound-state criterion is k2 |V_well| R^2 of order one with k2 of order eps^{-2}. So in the controlled regime eps -> 0 every channel binds, the onset is a first-order transition or a gas of bound clusters (Bak's result for commensurate-incommensurate transitions with oscillating soliton interactions), and the lowest cluster species generically wins alone. Oscillatory tails can give the miscible sign only where binding fails, that is at eps of order one, outside the blow-up.

(3) Emergent pseudospin conservation [derived, confirmed]. In K, P_pm = pi -+ delta with delta = 0.443 pi. A single valley flip changes the crystal momentum by 2 delta, a pair conversion (++ -> --) by 4 delta = 1.77 pi, neither 0 mod 2 pi. The flip-flop (+- -> -+) conserves momentum and is the J_perp term of (1.3), which conserves S^z. At order m a process changing S^z by m must absorb m (2 delta) mod 2 pi from the charge sector, whose low-energy momenta are 2 pi n j. This is possible only on the lines 2 m delta = 2 pi (k + n j), a measure-zero set of densities, so S^z_tot is an emergent U(1) at every order. The relative-phase vertex (single flip) carries momentum 2 delta (the coordinator's 2 P_0). So the miscible side is indeed a gapless relative-phase boson with c = 2, and Eisenberg-Lieb (Delta = -1) is its ferromagnetic boundary. The quantitative problem is that "easy plane" means |J_z| < |J_perp|, and J_perp (valley exchange, by pass-through or flip-flop) is of order exp(-S/hbar_eff) while J_z is a power of eps (Lemma S3). So the easy-plane window is exponentially thin around the surface J_z = 0. It is open, and it exists if and only if J_z changes sign inside a non-binding two-valley region. That sign change is the precise construction target.

(4) J1-J2 templates [derived comparison]. In the J1-J2 chain near saturation the two components are magnons at plus or minus k_0, point-like quantum particles with lattice-scale contact vertices. Their exchange (the inter-valley vertex) is of the same order as the intra-valley vertex, so the pseudospin anisotropy Delta is of order one and the miscible window is a finite range of J2/J1. For blow-up solitons the analogue of J_perp is exponentially small and Delta is generically huge. So matching to Kolezhuk et al. works structurally (same chain (1.3), same Eisenberg-Lieb boundary) but not quantitatively. The soliton route reaches J1-J2-like values of Delta only when hbar_eff is of order one, that is eps of order one, where the soliton is a few sites wide and the blow-up gives no control.

Feasibility of the lattice test [derived estimate]. With eps = 0.1 the soliton width is l of order 10 sites. A dilute gas needs spacing much larger than l, say 30 sites, so n <= 0.03 and L = 256 holds about 8 solitons. The pseudospin scale is c = pi^2 n^3/(3M) with M of order 1/(4 eps e''), about 4e-5 e'' in units of the bond, times exp(-S/(2 eps)) = exp(-5 S) for J_perp. A gap or velocity of order 1e-7 cannot be resolved by DMRG or FES at L <= 256. So a lattice test is meaningful only for eps of order one, which is again outside the controlled theory.

Conclusion of the construction side. No candidate bond was found, so no Typhon job was queued. The most promising remaining region is the one points (1) to (4) all lead to: N bonds at eps of order one, with a two-minimum soliton band produced by a mechanism other than the K tilt (which always binds), self-repulsive slow tails in both valleys, and a valley-dependent fast channel with the sign J_z > 0 of (2.3). That region is not reachable by the blow-up and would have to be searched on the lattice directly, by an adversarial scan for two-component entanglement scaling near soliton onsets (C155-style, but targeted at bonds whose classical soliton band has two near-degenerate minima).
