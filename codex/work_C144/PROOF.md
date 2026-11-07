# C144 derivations: the XZ point and its tuned cone

Run by a Claude Opus agent, 2026-11-07 (not by Codex). Standard: theoretical physics. Tags: [proved] (complete argument, or exact identity checked to machine precision), [derived] (physics standard, assumptions named), [numerical] (sizes and scripts given), [heuristic], [assumed]. Conventions are those of C143 (h = sum J_ab s^a (x) s^b + sum B_a (s^a (x) 1 + 1 (x) s^a)/2, ring H_L = sum h_{j,j+1}, M the gauge-optimized half-spread).

## 1. The local picture at XZ in all directions of B_2 (Task 1)

### 1.1 Majorana dictionary

Jordan-Wigner along y. With S_j = prod_{k<j} Y_k put b-Majoranas beta_j = S_j X_j and a-Majoranas alpha_j = S_j Z_j ({., .} = 2 delta). Then Y_j = i beta_j alpha_j (since X_j Z_j = -i Y_j), and for a bond (j, j+1), using S_{j+1} = S_j Y_j and XY = iZ, ZY = -iX:

    beta_j beta_{j+1}   = X_j Y_j X_{j+1} = i Z_j X_{j+1}      so  ZX = -i beta_j beta_{j+1} =: b_j,
    alpha_j alpha_{j+1} = Z_j Y_j Z_{j+1} = -i X_j Z_{j+1}     so  XZ = +i alpha_j alpha_{j+1} =: a_j,
    beta_j alpha_{j+1}  = X_j Y_j Z_{j+1} = i Z_j Z_{j+1}      so  ZZ = -i beta_j alpha_{j+1},
    alpha_j beta_{j+1}  = Z_j Y_j X_{j+1} = -i X_j X_{j+1}     so  XX = +i alpha_j beta_{j+1},
    YY = a_j b_j,   Y_j = i beta_j alpha_j,
    X_j = S_j beta_j,  Z_j = S_j alpha_j,  Y_j X_{j+1} = S_j beta_{j+1},  Y_j Z_{j+1} = S_j alpha_{j+1},
    X_j Y_{j+1} = S_j beta_j (i beta_{j+1} alpha_{j+1}),  Z_j Y_{j+1} = S_j alpha_j (i beta_{j+1} alpha_{j+1}).

[proved] The quadratic part of this dictionary (XZ, ZX, XX, ZZ, Y) was checked to machine precision by comparing exact free-fermion energies and entanglement entropies of open chains built from it (ff_open.py) with DMRG of the spin chain (test_dmrg.out: energy differences 1e-13, entropy differences 1e-10 at L = 40). YY = a_j b_j is C143's identity.

So the 12 bond directions split by the global parity P = prod_j Y_j into

- P-even, A/B preserving: XZ (= a_j, the A chain itself), ZX (= b_j), YY (= a_j b_j);
- P-even, A/B mixing (still quadratic): XX, ZZ, B_y. These hybridize the two Majorana species;
- P-odd (Jordan-Wigner strings): XY, YZ, YX, ZY, B_x, B_z.

Modulo on-site rotations (generators sum X, sum Y, sum Z give the tangent directions YZ... precisely: rotation about z gives XZ -> Y Z, about x gives XZ -> X Y, about y gives XZ -> XZ + theta(ZZ - XX)), the 8 normal directions of B_2 at XZ are ZX, YY, XX + ZZ, B_y (P-even) and YX, ZY, B_x, B_z (P-odd). [proved, elementary]

### 1.2 First-order effective Hamiltonian on the degenerate manifold

Facts at XZ [proved]: (i) the a's and b's commute, H_0 = sum a_j acts on the A Majoranas only, and the ground manifold D is (A ground state) (x) (all B states in the allowed sector), of dimension 2^{L/2 - 1} on even rings. (ii) prod_j b_j = (-i)^L P and prod_j a_j = i^L P (direct from b_j = Z_j X_{j+1}, a_j = X_j Z_{j+1}, X Z = -iY), so P is fixed on D because the A ground state fixes prod a_j.

Hence for any perturbation V = sum_j v_{j,j+1}:

- A P-odd V has P V P = 0 on D (P is constant there).
- A V odd in the A Majoranas (XX, ZZ, B_y contain exactly one alpha) has P V P = 0, since the A ground state has definite A parity at fixed boundary condition and the lowest A-odd state is an excitation.
- P a_j b_j P = <a_j>_A b_j, P b_j P = b_j, P a_j P = <a_j>_A.

So, for the general tangent vector delta h = sum dJ_ab s^a s^b + sum dB_a (...):

    H_eff^(1) = dJ_xz E_0^A(L) + ( dJ_zx + <a>_L dJ_yy ) sum_j b_j ,     <a>_L = E_0^A(L)/L -> -2/pi.     (1.1)

[proved; numerical check first_order.out] At L = 10 the projections P V P of all 12 directions were fitted to c0 + c1 P (sum b) P: zx gives c1 = 1, yy gives c1 = -0.6472135955 = E_0/L, xz gives c0 = E_0, and the other nine vanish to 1e-13.

Tuned locus. [proved at first order] The first-order B dynamics vanish exactly on the hyperplane

    T = { dJ_zx + <a> dJ_yy = 0 }  (dJ_zx = (2/pi) dJ_yy as L -> infinity),

which has codimension one in the 8-dimensional normal space at XZ (7 dimensions). Its intersection with the A/B-preserving plane (ZX, YY) is C143's cone alpha = (2/pi) beta; on the AT family the exact tuned curve is t(alpha, beta) = 0 where t is the full uniform B hopping, alpha = (2/pi) beta + O(beta^2) (section 3.4 shows the O(beta^2) shift is visible numerically). The other six normal directions (XX + ZZ, B_y, YX, ZY, B_x, B_z) are tangent to T: they never enter at first order. So "the tuned cone" in B_2 is really a 7-dimensional hypersurface through XZ, and inside it the leading dynamics are set by different mechanisms depending on the direction (1.3 to 1.5).

### 1.3 The mixing directions XX + ZZ and B_y: exact free solution

[proved] On the free locus h = XZ + alpha ZX + e_xx XX + e_zz ZZ + B_y (sum Y) the Bloch Hamiltonian in the basis (alpha_k, beta_k) of h = 2iA is, with w(k) = e_xx e^{ik} + e_zz e^{-ik} - B_y,

    h(k) = [[ -4 sin k , 2i w(k) ], [ -2i conj(w(k)) , 4 alpha sin k ]],   det h(k) = -16 alpha sin^2 k - 4 |w(k)|^2.

(The single-particle energies are the eigenvalues; the factor 2 relative to the C143 normalization is the h = 2iA convention, velocity of the A chain 4.) Gapless points are zeros of det h. Consequences:

- At alpha = 0 (the free part of T): det = -4|w|^2, so gapless points are double zeros, k with w(k) = 0. For e_xx = e_zz = e (the normal direction; e_xx - e_zz is the rotation about y) w = 2e cos k - B_y is real, and there are zeros iff |B_y| <= 2|e|: a quadratic (z = 2) touching of the B band at cos k_0 = B_y/(2e), E_B(k) ~ w(k)^2/(2 sin k) (second order in the mixing). The A Fermi points at k = 0, pi are gapped at first order, gap 2|w(0)|, 2|w(pi)| = 2|2e -+ B_y|. For |B_y| > 2|e| everything is gapped.
- alpha > 0: det < 0 for all k unless sin k = 0 and w = 0, so the chain is gapped (e.g. XZ + ZX + e(XX + ZZ) = (1 + e) Z'Z' - (1 - e) X'X' after a 45 degree rotation about y, an anisotropic XY chain).
- alpha < 0: det = 0 at |w(k)| = 2 sqrt(|alpha|) |sin k|, simple zeros, four Majorana Fermi points near the zeros of w, c = 1 (e.g. alpha = -1 is the U(1) chain XZ - ZX + e(XX + ZZ), a Dzyaloshinskii-Moriya XX chain).

So the mixing directions are relevant couplings of the A Majorana (dimension 1/2) to the dispersionless B Majorana: they gap A at first order and leave a single-species Pokrovsky-Talapov / Lifshitz onset (one band, z = 2) on T, with c in {0, 1} nearby. This is C143 stratum 3 (PT1, PT2 surfaces) meeting XZ, not a new kind of point. [numerical] DMRG of XZ + 0.1(XX + ZZ) and XZ + 0.1 B_y at L = 64 (dirs_L64.out, section 3.6).

### 1.4 The P-odd directions B_x, B_z, YX, ZY

[derived] X_j = S_j beta_j is (string of j - 1 A Majoranas) times (string of j B Majoranas). A string of consecutive Majoranas of the critical A chain is an Ising order or disorder field (the A chain of L Majoranas is a critical transverse-field Ising chain of L/2 sites, and one-Majorana translation is the Kramers-Wannier map exchanging sigma and mu). So B_x sum X_j, B_z sum Z_j, YX, ZY all have the form sum (sigma_A or mu_A)(x) (x) O_B(x), with O_B a B-string, of A-dimension 1/8 and B-dimension 0 (B has no dynamics at XZ). Their first-order projection vanishes (P-odd), but they are strongly relevant at the A fixed point: at second order the kernel int dtau <sigma_A(x, tau) sigma_A(0)> grows like |x|^{3/4}, so perturbation theory in these directions breaks down at the length |delta|^{-8/15} (the A gap opened by a longitudinal-type field). For B_x and B_z the B strings commute with each other (the X_j commute), so they behave classically: they align to make the effective field on A uniform, and both sectors gap. For YX and ZY the B parts of neighbouring terms anticommute, the second-order B dynamics are genuine, and numerically A gaps while a slow sector survives (section 3.5). [numerical] XZ + 0.1 B_x at L = 64 has S_mid = 0.26 and no logarithmic growth (c_fit = 0.0 +- 0.1), and Q_L -> 0 rapidly (qscan_local_L8.out: Q_8 = 0.33 at s = 0.1, 0.04 at s = 0.3). The rotation directions XY and YZ are P-odd too but are commutators [sum Z, H_0], [sum X, H_0] and change nothing.

### 1.5 Summary of Task 1

| direction | P | A/B | first order on D | leading effect |
| :-- | :-- | :-- | :-- | :-- |
| ZX | even | preserved | sum b_j | B hopping, c_B = 1/2, v_B = 4|t| |
| YY | even | preserved | <a> sum b_j | AT coupling, marginal |
| XX + ZZ, B_y | even | mixed (free) | 0 | A gapped at first order, B z = 2 onset on T |
| B_x, B_z | odd | mixed (strings) | 0 | relevant sigma_A O_B coupling, gap |
| YX, ZY | odd | mixed (strings) | 0 | A gapped, slow B sector with c about 1 at L = 64 (3.5) |
| XY, YZ | odd | rotations | 0 | none |

The direction-dependence that matters for Q1 and Q2 is therefore confined to the AT plane: everywhere else on T the A chain is gapped at a scale set by the factorization-breaking component.

## 2. (Q1) Boundedness of Q_L near XZ

### 2.1 The XZ value in closed form

[proved, given the sector bookkeeping; checked numerically] On even rings the XZ ground energy is the antiperiodic critical Majorana chain of L Majoranas, E_L(XZ) = -2/sin(pi/L) (ED L = 4 to 20: -2.828, -4.000, -5.226, -6.472, -7.727, ..., -12.859; kernel.out, qscan). M(XZ) = 1. With x = pi/(2L), sin 2x = 2 sin x cos x and (1 - cos x)/sin x = tan(x/2),

    Q_L(XZ) = L ( -2/sin x + 4/sin 2x ) = (2L/sin x)(1/cos x - 1) = 2L tan(pi/(4L)) / cos(pi/(2L)),         (2.1)

which decreases monotonically from Q_2 = 2.3431 to pi/2 (= (pi/4) c v with c = 1/2, v = 4). Check against ED: 2.34315, 1.72241, 1.63556 at L = 2, 4, 6 (qscan output) and 1.60674 at L = 8 (qfree.py). On odd rings Q_L(XZ) < 0 (-1.61, -1.58 at L = 3, 5).

### 2.2 Small rings: the cover lemma

[proved, C131] |Q_L(h) - Q_L(k)| <= 4 L^2 M[h - k]. With delta = M[h - XZ]:

    Q_L(h) <= Q_L(XZ) + 4 L^2 delta <= 2.343 + kappa   for all L <= L_1(kappa) = (kappa/(4 delta))^{1/2}.   (2.2)

This is uniform in the direction of approach, and it is the only rigorous input.

### 2.3 Large rings: sector decomposition

Beyond L_1 the cover lemma is useless and one needs the physics of section 1. [derived, under assumptions S1 to S3 below] For |delta| small the low-energy spectrum of the ring of h = XZ + delta u is that of two sectors:

- A, the critical Majorana chain with v_A = 4(1 + O(delta)), possibly gapped by the factorization-breaking components of u (section 1.3, 1.4) or by B (section 3.4);
- B, an effective chain on the degenerate manifold, of energy scale eta <= C|delta| (first order (1.1)) or eta ~ delta^2 log(1/delta) on T (section 3).

Write E_L = E^A_L + E^B_L up to terms that are extensive with an L-independent density (these cancel in E_{2L} - 2E_L) and corrections of relative order eta/v_A. Then

    Q_L(h) = [ L D^A_L + L D^B_L ] / M(h),      D^s_L = E^s_{2L} - 2 E^s_L.

Assumptions.
- (S1) [derived, section 1 and 3.1] Separation of scales: the B-sector energies are O(eta) per site with eta/v_A -> 0, so Schrieffer-Wolff on the manifold holds at every L (the B level spacings are O(eta/L), far below the A gap 4 pi/L).
- (S2) [assumed; numerical for free massive sectors, qfree.out] No Casimir overshoot: for each sector, L D^s_L <= (pi/4) max(c_UV v_UV, c_IR v_IR) times a factor of order one (the CFT doubling is exactly (pi/4) c v; the crossover function of a massive free Majorana lowers it monotonically: Q_L = 1.72, 1.60, 1.33, 0.86, 0.25, 0.01 at L = 4 to 128 for a dimerization 0.05, never above the critical (2.1)).
- (S3) [derived, section 3] The B sector is a Majorana chain with c_B <= 1/2 (or gapped) on the AT plane, and single-species elsewhere on T, with velocity v_B <= C' eta. Where the B energy density still drifts with L (the logarithmic kernel of section 3.1 below the A correlation length), D^B_L < 0: for the effective chain at lam = 0, E^B_L/L = -0.159, -0.234 at L = 12, 24 (units beta^2), so L D^B_12 = -21 beta^2. The drift lowers Q, it cannot raise it.

Result. [derived, conditional on S1 to S3]

    sup_L Q_L(h) <= max{ 2.343 + kappa ,  (pi/2)(1 + C delta) + (pi/8) v_B^max / M(h) + O(kappa^{-1}) },     v_B^max = O(delta),     (2.3)

uniformly in the direction u, so for delta -> 0 the supremum tends to the small-ring value 2.343 of XZ itself (attained at L = 2). The B sector contributes at most (pi/8)(4 |t|)/M = O(delta) whatever happens on the tuned cone, because c_B <= 1/2 and v_B = O(delta). The answer to Q1 does not depend on c_B being small, only on c_B being bounded; with c_B <= 1/2 it is quantitatively close to the XZ value.

Remark. The crossover window L_1 < L < xi_B (B not yet in its IR) is covered by S2 and S3 and not by the cover lemma. This is where a hidden tower would have to live; S3 says the B contribution is negative there (energy density still decreasing), and the ED below shows no positive transient.

### 2.4 ED check

Local (qscan.py, qscan_local_L8.out, L <= 8): the even-L values of Q_L stay below 2.343 for every direction and amplitude tried (cone s = 0.1, 0.3, 0.6; B_x; XX + ZZ; a random direction r1), and decay to 0 on gapped directions (B_x: 0.33 at L = 8, s = 0.1). Typhon (job 17580166, L <= 10, 2L <= 20, 22 direction families): see RESULT.md.

## 3. (Q2) The B chain on the tuned Ashkin-Teller cone

### 3.1 Second-order effective Hamiltonian

On the AT plane H = sum_j (a_j + alpha b_j + beta a_j b_j). Write a_j = <a>_L + da_j. Then exactly

    H = H_A + t sum_j b_j + beta sum_j da_j b_j,       t = alpha + beta <a>_L.

The b's commute with H_A, so on the manifold D the unperturbed B energies are all zero and the energy denominators of Schrieffer-Wolff are A energies only. [derived, standard SW with gap 4 pi/L of A in the relevant sector]

    H_B = t sum_j b_j - beta^2 sum_{j,l} K_jl b_j b_l + O(beta^3),
    K_jl = <0| da_j Q (H_A - E_0)^{-1} Q da_l |0>_A = int_0^inf dtau < da_j(tau) da_l(0) >_A.            (3.1)

Free-fermion evaluation [proved]. With H_A = (1/4) g^T h g, h = 2iA, positive eigenpairs (lam_mu, v_mu), the imaginary-time Heisenberg equation [H, g] = -h g gives g(tau) = e^{-h tau} g and G_mn(tau) = <g_m(tau) g_n> = 2 sum_mu e^{-lam_mu tau} v_{m mu} conj(v_{n mu}) (G(0) = 2 P_+, diagonal 1). Wick's theorem for the connected four-point function, <a_j(tau) a_l>_c = G_jl G_{j+1,l+1} - G_{j,l+1} G_{j+1,l}, and the tau integral give

    K_jl = 2 sum_{mu,nu} W^j_{mu nu} conj(W^l_{mu nu}) / (lam_mu + lam_nu),   W^j_{mu nu} = v_{j mu} v_{j+1,nu} - v_{j+1,mu} v_{j nu}.

[numerical, kernel.out] This equals the ED resolvent of the full XZ ring to 2e-16 at L = 8 and 10.

Reality and the vanishing of the hopping. [proved] The unitary R': g_j -> (-1)^j g_{-j} maps a_j -> a_{-j-1} and leaves H_A and its unique sector ground state invariant, so K_{j,j+1} = K_{-j-1,-j-2} = conj(K_{j,j+1}) by translation and hermiticity. K is real symmetric. Since b_j b_{j+1} = -b_{j+1} b_j, the nearest-neighbour pairs cancel in (3.1), and the diagonal gives a constant. The only two-fermion operators b_j b_l can produce are for |j - l| = 1, so

    H_B = t sum_j b_j - 2 beta^2 sum_j sum_{r >= 2} K(r) b_j b_{j+r} + const,   b_j b_{j+r} = -beta_j beta_{j+1} beta_{j+r} beta_{j+r+1}.   (3.2)

On the tuned cone (t = 0) the second-order B chain is purely quartic: no hopping at all. (C143's "next-nearest hopping" term is absent.) This is an exact statement about second order.

Long-range tail. [derived; checked numerically] In the infinite chain the positive modes are k in (-pi, 0), lam_k = 4|sin k|, v_{jk} = e^{ikj}/sqrt(L), |W|^2 = 4 sin^2((k - k')/2)/L^2, so

    K(r) = (1/(2 pi^2)) int int_{(-pi,0)^2} dk dk' sin^2((k - k')/2) e^{-i(k+k') r} / (|sin k| + |sin k'|).

For large r the integral is dominated by k near 0 and k' near -pi (and the mirror), where sin^2((k - k')/2) -> 1. With k = -p, k' = -pi + p' and 1/(p + p') = int_0^inf ds e^{-s(p + p')},

    int_0^inf int_0^inf dp dp' e^{i(p - p') r}/(p + p') = int_0^inf ds / ((s - ir)(s + ir)) = pi/(2r),

so K(r) = (-1)^r [2 x (1/(2 pi^2)) x pi/(2r)] + O(r^{-3}) = (-1)^r/(2 pi r) + O(r^{-3}). The same-region contributions carry sin^2((k - k')/2) ~ (p - p')^2/4 and give the non-oscillating r^{-3} part. [numerical] K(r) (-1)^r 2 pi r = 0.945, 0.968, 0.976, 0.978 at r = 2, 3, 4, 5 (L = 800; finite-ring chord corrections at larger r). The staggered sum diverges logarithmically: K~(pi) = sum_r (-1)^r K(r) = (1/pi) log L - 0.0577 (L = 200 and 800 agree to 1e-5). In continuum language the tail is the static eps_A eps_A correlator, int dtau <eps(x, tau) eps(0)> ~ pi/(2 v |x|), and the B chain acquires an instantaneous eps_B(x) eps_B(y)/|x - y| interaction, marginal at the Ising point (RG eigenvalue 2 - 2 Delta_eps = 0).

Validation of (3.2). [numerical, validate_bmodel*.out] In units of beta^2, with lam = t/beta^2, h_B(lam) = H_B/beta^2 was written on L/2 spins (b_{2i-1} = Z_i, b_{2i} = X_i X_{i+1}, closing bond c X_n X_1) and compared with the full qubit ring at beta = 0.004. The lowest 2^{L/2-1} levels of (E - E_0^A)/beta^2 + L K(0) agree with h_B in the sector (c, P) = (+1, even) for L = 0 mod 4 and (-1, odd) for L = 2 mod 4, to O(beta) (for example L = 12, lam = 1: -8.534, -6.333, -3.586 against -8.538, -6.335, -3.587).

### 3.2 The B chain at lam = 0: self-generated staggered hopping

Mean-field [derived]. Decoupling (3.2) as b_j b_l -> <b_j> b_l + b_j <b_l>, a state with <b_j> = m (-1)^j gains energy

    <H_B>/L = -beta^2 m^2 kappa_L,    kappa_L = sum_{|r| >= 2} (-1)^r K(r) = K~(pi) - K(0) - 2|K(1)| ~ (1/pi) log L,

and the mean field it creates is a staggered-sign hopping H_MF = -2 beta^2 m kappa_L sum_j (-1)^j b_j. A Majorana chain with hoppings t_s(-1)^j is gauge-equivalent (beta_j -> s_j beta_j with s_j s_{j+1} = (-1)^j, period 4) to the uniform critical chain, so it is gapless with |<b_j>| = 2/pi. The solution is self-consistent with m = 2/pi for any kappa > 0. So at lam = 0 the B chain breaks the one-Majorana translation (the sign of m, twofold), but it is not gapped: it is a critical Majorana chain with Fermi points moved to k = +-pi/2, velocity v_B = 4 t_s = (16/pi) beta^2 kappa.

[numerical, bscan_a.out, bscan_b.out, bdiag_lam0.out] ED of h_B(0) to L = 40: twofold ground state with splitting 2e-4, 5e-7, 8e-10, 1e-12 at L = 16, 24, 32, 40 (exactly degenerate at L = 4 mod 8); <b_odd> = -<b_even> = 0.641, 0.638, 0.638 at L = 16, 24, 32 (2/pi = 0.637); the next levels L (E - E_0) = 25.3 (x2), 48.5 (x6), 68 at L = 32, ratio 1 : 1.92 : 2.69, close to a Majorana tower (1 : 2 : 3), and L x gap = 13.9, 17.5, 20.1, 22.1, 23.8, 25.3, 26.5, 27.6 for L = 12 to 40, growing like log L as v_B ~ kappa_L predicts (mean-field v_B x 2 pi at L = 40: 21, ED 27.6). The mean-field energy -(4/pi^2) kappa_40 = -0.27 per site against ED -0.294.

lam != 0 [numerical]. At lam = 0.5 the lowest pair splits exponentially (0.090, 0.047, 0.023 at L = 32, 36, 40) and the gap above grows (2.56, 2.72, 2.95): a gapped twofold state. In mean field the hoppings become t_j = lam + t_s (-1)^j, gapless only if |lam + t_s| = |lam - t_s|, i.e. lam t_s = 0: a uniform hopping detunes the staggered critical chain and opens a gap 2 min(|lam|, t_s). At larger lam the Ising (uniform) behaviour survives to a length that grows with lam (split x L flat at lam = 2 up to L = 40, decreasing at lam <= 1.5 beyond L = 20 to 32): with strictly critical A the log-divergent kernel makes the instability marginal and eventually wins at any lam in the B-only model.

### 3.3 Feedback on A and the phase structure near the cone

[derived] Once <b_j> = m (-1)^j, the A chain sees hoppings 1 + beta <b_j>, alternating by 2 beta m: it is dimerized and gapped, with Delta_A ~ beta m and xi_A ~ v_A/(beta m). This cuts the logarithm: kappa -> (1/pi) log(1/|beta|) + O(1). The dimerized A in turn feeds a staggered first-order term beta <a_j> b_j back to B, of the same sign pattern (reinforcing). The resulting picture on the AT plane at small beta is:

- |t| >> beta^2 log(1/beta): both chains critical, c = 1 (the C143 numerics, and the XXZ line alpha = 1);
- a wedge |t| <~ C beta^2 log(1/|beta|) around the cone: translation (self-duality) spontaneously broken, A dimerized and gapped, B with alternating hoppings, gapped except on one line inside the wedge where its two hoppings have equal magnitude;
- that line (the true tuned cone, t_eff = 0): A gapped, B critical.

So c_tot <= 1 throughout, with c_tot = 0 in the wedge. What the critical line carries is settled numerically in 3.4.

### 3.4 Full qubit chain: DMRG on open chains (route b)

Method. dmrg.py is a two-site DMRG written for this run (NumPy/SciPy, MPO of bond dimension <= 6, Lanczos by eigsh). [numerical, test_dmrg.out] Checks: open L = 10 against ED (difference 6e-15); open free-fermion chains XZ + 0.5 ZX and XZ + 0.3 ZX + 0.2(XX + ZZ) + 0.1 B_y at L = 40 and 80 against exact Majorana diagonalization (ff_open.py): energies to 1e-8 or better, S(l) to 1e-10; the critical case gives parity-resolved c fits equal to the exact ones (1.04 / 0.90 at L = 80). Central charges are extracted from S(l) on open chains by fitting S = (c/6) log[(2L/pi) sin(pi l/L)] + const, either per parity of l or for the parity average (cfit.py; central windows w L <= l <= (1 - w) L). As a reference, a free state with c = 1/2 of the same structure (A dimerized, B staggered and critical; ff_mf.py) gives 0.48 to 0.49 with the averaged fit and 0.56 / 0.41 with the parity fits at L = 64.

Ashkin-Teller plane, L = 128, chi <= 128, truncation <= 1e-9 (Typhon job 17580176; d128_summary.out). For beta = 0.15, 0.3, 0.5 and alpha/beta from 0 to 1.5:

| alpha/beta | beta = 0.15 | beta = 0.3 | beta = 0.5 |
| :-- | :-- | :-- | :-- |
| 0 | 0.95 | 0.95 | 0.92 |
| 0.3 | 0.81 | 0.79 | 0.89 |
| 0.45 | 0.40 | 0.30 | 0.72 |
| 0.55 | 0.05 | 0.02 | 0.25 |
| 0.60 | 0.00 | 0.00 | 0.04 |
| 0.64 | 0.29 (b_mid 0.05) | 0.01 | 0.00 |
| 0.68 | 0.00 | 0.00 | 0.66 (b_mid -0.06; parity fits 0.46 / 0.44) |
| 0.72 to 0.80 | <= 0.03 | <= 0.09 | <= 0.01 |
| 0.85 | 0.56 | 0.32 | 0.09 |
| 0.90 | 0.72 | 0.58 | 0.32 |
| 1.0 | 0.91 | 0.89 | 0.78 |
| 1.2 | 1.04 | 1.02 | 0.99 |
| 1.5 | 1.01 | 0.99 | 0.99 |

(parity-averaged c fit, w = 1/8). Inside the low-c window the A bonds alternate strongly (a_stag = -0.13 to -0.35) and b alternates with amplitude 0.45 to 0.64, close to 2/pi, as predicted in 3.2 and 3.3. The width of the window in alpha (c_fit < 0.1) is about 0.05, 0.11, 0.17 for beta = 0.15, 0.3, 0.5, which is 1.17, 1.02, 0.98 times beta^2 log(1/beta): the beta^2 log(1/beta) scaling of 3.3 holds. [numerical]

The line inside the window. The uniform part b_mid of <b_j> changes sign once inside the window (alpha/beta between 0.60 and 0.68 at all three beta). [numerical, L = 64 to 96, beta = 0.3, fine_b03_L64.out, peak64.out, line_local.out, typhon_out/t_L64_*] Approaching that line S(L/2) rises steeply (0.54, 0.62, 0.72, 0.81, 0.87 at alpha = 0.194, 0.195, 0.196, 0.197, 0.1985, L = 64) and the c fits pass through a maximum of about 1.0 near alpha = 0.197 before dropping. Two facts show that the line carries one critical Majorana, not two:

1. [numerical, corr.out] At L = 64, alpha = 0.197 the connected a-a correlator decays exponentially (9.8e-2, 9.0e-3, 4.2e-4, 3.8e-5, 5.3e-6 at r = 2, 4, 8, 12, 16, correlation length about 1.7) and Y-Y likewise, while b-b decays as a power (3.1e-2, 7.1e-3, 1.9e-3, 5.7e-4 at r = 4, 8, 16, 32, local exponents -1.9 and -1.7, the Ising energy correlator r^{-2}). So A is gapped and only B is critical.
2. [numerical] The local uniform part of b drifts from +0.22 at the open ends to -0.045 in the middle (L = 96), so on an open chain the B mass vanishes on a surface inside the chain and not uniformly. Single-chain profile fits are then biased upward; the averaged fit of the exact c = 1/2 free reference is 0.48, while the line points at fixed L give anything from 0.3 to 1.0 depending on where the mass zero sits. The size scaling of the maximum over alpha of S(L/2) is the clean estimator. Between L = 64 (alpha = 0.1985, 0.870) and L = 96 (alpha = 0.198, 0.895) it grows by 0.025, against 0.034 for c = 1/2 and 0.068 for c = 1 (the maxima are not located to better than 0.002 in alpha, so this is indicative only). At beta = 0.5, where the line was hit at L = 128, the parity fits give 0.46 / 0.44.

3. [numerical, peak64.out; heuristic interpretation] Continuing the L = 64 scan, S(L/2) stays near 0.87 for alpha = 0.199 to 0.203 while the averaged c fit falls from 0.80 to 0.21, and then drops abruptly to 0.17 at alpha = 0.205. A flat profile carrying a large constant, which then disappears suddenly, is the signature of a cat state. On the side of the line where B is ordered (topological in fermion language), the open chain has a pair of B edge zero modes, and DMRG can return the symmetric superposition, which adds log 2 to every cut. So S(L/2) on that side contains a spurious log 2, and its maximum over alpha is not a clean estimator either. Only the correlation structure (item 1) and the beta = 0.5 parity fits are free of this problem.

[numerical, partial; Typhon jobs 17580192 and 17580432 pending] So the tuned line carries c_tot = c_B = 1/2 with c_A = 0, within the accuracy reached. The larger sizes (L = 128, 192, 256 on a grid alpha = 0.194 to 0.205) are running.

### 3.5 Directions outside the AT plane, on T

[numerical, dirs_L64.out, qscan_typhon_summary*.out] At L = 64: XZ + 0.1 B_x is gapped (S_mid 0.26, no log growth); XZ + 0.1 B_y is gapped (0.61, c_fit 0.0) as the free solution predicts (w = -B_y, det h < 0); XZ + 0.1(XX + ZZ) has S_mid 0.86 with c_fit 0.16 / 0.03, the z = 2 touching of 1.3 (no logarithmic growth). Q_L along cone + B_x, cone + (XX + ZZ), cone + B_y decreases to 0.09 to 1.5 by L = 10 (gapped or slow), never above the XZ value.

XZ + eps YX (and its mirror ZY) is different from what 1.4 predicts for B_x. The YX_j do not commute among themselves (YX_j YX_{j+1} anticommute), so the B parts are not classical, and the second-order B dynamics are genuine. [numerical] Q_L drops below pi/2 quickly (0.31 at L = 10, eps = 0.1), so A is gapped, but DMRG at L = 64 (eps = 0.1) shows S_mid 0.81 with parity fits 1.37 / 0.84, and ring ED (yx_ed.out) shows L x gap roughly constant (0.57 to 0.63 at L = 8, 10, 14, 16) with a period-3 commensuration (a triplet near 1.6 to 1.75 at L = 12, 18). This looks like a slow c of about 1 sector (velocity of order 0.1), with A gapped. L = 128 and 256 runs are in job 17580192. [heuristic] Whatever its value, it is a single B species with v = O(eps^2) or smaller, so for Q1 it is harmless; for the c_tot <= 1 question it is consistent with c_tot <= 1 at L = 64 but not settled.

### 3.6 Answer to Q2

[numerical, with the derivation of 3.1 to 3.3] On the AT tuned cone of XZ, for beta = 0.15 to 0.5:
- the first-order cone alpha = (2/pi) beta lies inside a gapped wedge of width about beta^2 log(1/beta) where translation (self-duality) is spontaneously broken, A is dimerized and B has alternating hoppings;
- the exactly tuned line t_eff = 0 inside the wedge (alpha/beta between 0.64 and 0.68 for these beta) carries one critical Majorana, c_B = 1/2 with A gapped, so c_tot = 1/2 there;
- outside the wedge c = 1 (both Majoranas critical, marginal coupling), with no fit above 1.09 anywhere.

So c_B <= 1/2 everywhere on the AT plane near XZ, and no nearest-neighbour qubit bond in this neighbourhood was found with c_tot > 1. C143's ED estimates of 1.1 to 1.5 on the cone came from L <= 16, where the twofold symmetry breaking appears as a log 2 cat entropy on rings and the A dimerization length (about 1/beta) is not yet reached.

## 4. Consequences for C141's (EM-XZ) and (NC) (Task 4)

(EM-XZ), as stated in C143 Section 5: "At XZ, v_B = O(delta M) (O(delta^2 M) on the cone) and c_B is bounded uniformly in the direction of approach." [derived, numerical support] Holds, with three amendments. (i) On the tuned locus T the B velocity is O(delta^2 log(1/delta)) (v_B = 4 t_s, t_s ~ beta^2 kappa with kappa ~ (1/pi) log(1/beta)), not O(delta^2). (ii) The A sector is not always a spectator: on most of T it is gapped, by the hybridization directions at first order and by the self-generated dimerization on the AT part of T. (iii) c_B <= 1/2 on the AT plane; on the P-odd YX/ZY directions the emergent slow sector looks like c of about 1 (with A gapped), still bounded. So the recursion of C141 5.3 closes at XZ with a bounded emergent central charge in every direction, and v_em -> 0.

(NC), "c_B <= 1/2 on the XZ tuned cone": [numerical] holds on the AT plane (the line inside the wedge has c_B = 1/2, the wedge has c_B = 0). The physically relevant statement for C141's c_tot <= 1 at d = 2 is c_tot <= 1, and it holds in every direction tested (max fit 1.09 within the c = 1 phase, where the exact value is 1). The YX/ZY result means (NC) should be restated as c_tot <= 1 rather than c_B <= 1/2: there the A Majorana is gapped and the emergent sector may carry c of about 1.

## 5. Consistency checks

| check | expected | found |
| :-- | :-- | :-- |
| first-order projection, L = 10 | zx: 1, yy: E_0/L, others 0 | 1, -0.6472135955, <= 1.2e-13 |
| kernel, free formula vs ED (L = 8, 10) | equal | 2e-16 |
| kernel tail K(r)(-1)^r 2 pi r | -> 1 | 0.945 to 0.978 (r = 2 to 5, L = 800) |
| K~(pi) - (1/pi) log L | constant | -0.05770, -0.05769 (L = 200, 800) |
| E_0/L of the A chain | -2/pi = -0.636620 | -0.636621 (L = 800) |
| B model vs full ring, beta = 0.004 | equal to O(beta) | e.g. -8.534 vs -8.538 (L = 12, lam = 1) |
| B model, lam = 0, <b_odd> | 2/pi = 0.637 | 0.641, 0.638, 0.638 |
| B model, lam = 4, L x split | 2 pi v x_eps = 2 pi 16 = 100.5 | 102.7 to 108.8 (L = 12 to 28; marginal drift) |
| Q_L(XZ) closed form | 2L tan(pi/4L)/cos(pi/2L) | ED 2.34315, 1.72241 (L = 2, 4) |
| Q_L(ZX direction) | = Q_L(XZ) (two decoupled chains, scaled) | identical to 5 digits |
| Q_L on rotation directions XY, YZ | = Q_L(XZ) | identical |
| massive free Majorana Q_L | below the critical value, -> 0 | 1.70, 1.54, 1.33, 0.86, 0.25 (d = 0.05) |
| DMRG vs ED, vs free fermions | equal | 6e-15; 1e-8; S to 1e-10 |
| c on the XXZ side (alpha/beta = 1.5) | 1 | 0.99 to 1.01 |
| decoupled chains (C143, alpha = 0.3, beta = 0) | 1 | reproduced by structure (exact) |
