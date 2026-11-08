# C163 derivations: the attractive liquids at intermediate K

This run was done by a Claude Opus agent (2026-11-08), not by Codex. Astra was not used (out of credits).

Standard: theoretical physics. Tags: [proved] (complete argument, or an exact identity checked to machine precision), [derived] (physics-standard derivation with named inputs), [numerical] (sizes, scripts and outputs named), [literature], [heuristic], [assumed]. A tag applies to the paragraph or displayed claim it opens.

Conventions as in C159 PROOF. H = sum_j [ sigma_j^T J sigma_{j+1} + B . sigma_j ], Pauli matrices, J = S O (polar, O a rotation by theta about the spiral axis m). The spiral parent is XXZ about m with J_perp = s_perp, J_par = s_m, field m . B, and Delta = s_m/|s_perp| in the antiferromagnetic convention (Delta = cos gamma, K = pi/(2(pi - gamma)), K > 1 for Delta < 0, K = 15/8 at Delta = cos(11 pi/15) = -0.669). Bosonization sigma^m ~ (2/pi) d_x phi + a (-1)^j cos 2phi, sigma^+ ~ e^{i vartheta}[ b (-1)^j + c cos 2phi ] in the AF frame, x_{p,n} = p^2/(4K) + n^2 K for V_{p,n} = e^{i(p vartheta + 2n phi)}. "Harmonics" are the parts of the bond with U(1)_m charge q = +-1, +-2, which in the spiral frame carry lattice momentum q theta (C154 section 2).

The window W. By C159 Proposition C, a c > 1 candidate occupying an open set or a hypersurface of T or K needs a half-floating boson whose commensurate vertex V_c (= e^{2 i phi}, dimension x = K) obeys K > 15/8 (Ising x boson, sigma at momentum pi) or x1 + x2 > 2 (two bosons, the second with V_c at pi, or at momentum 0 with x2 > 2). C160 controls the K -> infinity corner next to h_F away from the twisted curve gamma = { -R(phi) } (C161 1c). W is the rest of the attractive half-floating liquids: K between 1 and the semiclassical range, plus the neighbourhood of gamma at phi of order one.

## 1. Mode inventory in the exact parents (Task 1, first half)

### 1.1 Zero-field parents (all of T and the sheet K_I)

[proved, C159 Lemma A and Lemma 1] The spiral parent of every T bond and of every K_I bond is at zero field (m . B = 0), and a half-floating liquid of T or K_I connected to it has zero magnetization along m (phi shift pi/2).

[literature, derived] The zero-field XXZ chain in the critical regime -1 < Delta < 1 has, below the lattice scale, only the states of one compact boson. The finite-size spectrum computed from the Bethe ansatz (de Vega and Woynarovich 1985, Woynarovich 1989, Klumper, Batchelor and Pearce 1991) is exactly the Gaussian tower x = p^2/(4K) + n^2 K + integers with c = 1, at every gamma, and its corrections are irrelevant operators. In the attractive regime Delta < 0 (gamma > pi/2) the Bethe states include strings (bound states of flipped spins). At zero magnetization these are the massless breathers of the massless sine-Gordon scattering description of the same boson (the attractive regime of the vartheta-type perturbations, whose vertex e^{i vartheta} has dimension 1/(4K) < 1/4). They are part of the c = 1 theory and carry no extra central charge. So the zero-field parent has no gapped branch below the multi-boson continuum, and nothing that is "almost soft". A second sector in T or K_I would have to come from lattice-scale states (energy of order J at the relevant momentum) brought down by harmonics of order J. This is the first structural fact about W. It is non-perturbative by construction.

### 1.2 Finite-field parents (N and the sheet K_II)

[literature] In a field the Bethe ground state is a Fermi sea of 1-strings, and the n-strings (n >= 2) are bound states with positive dressed energy at T = 0 (Takahashi, Thermodynamics of one-dimensional solvable models, ch. 9, gapless regime). They form gapped branches that become massless as h -> 0 and are bound magnon pairs (and n-tuples) near saturation.

[proved, gauge covariance] In the U(1) class the DM vector along m is the gauge twist sigma_j^+ -> e^{i j alpha} sigma_j^+, so every state of U(1)_m charge q has its momentum shifted by q alpha. Above the polarized state the n-magnon bound states of XXZ have a single-harmonic band E_n(P) = A_n - B_n cos(P - n alpha) (Bethe; for n = 2 this is C160's exact E_b(K) = 2 c0 + V + tau(K)^2/V with tau = 2|a| cos(K/2 - phi), checked there by ED to 1e-4). Hence the band minima sit at n alpha or n alpha + pi. Every soft wavevector of the U(1) parent lies in Lambda = Z alpha + Z pi.

[literature, established in the programme] In the U(1) class c = 1 at every gapless point, so no string branch ever reaches zero energy while the liquid survives (the saturation transition of -1 < Delta < 1 is single-magnon Pokrovsky-Talapov, and for Delta < -1 the reversal is first order by McGuire collapse, C160 Task 3).

### 1.3 Conversion locking: a bound-state branch cannot become an independent boson

Lemma L (conversion locking). [derived; inputs: two-component bosonization, the momentum bookkeeping of 1.2] Suppose that inside a liquid with phase vartheta_1 (vertex e^{i vartheta_1} at momentum Q_1) a branch of charge-n composites softens and forms a second boson vartheta_2 with e^{i vartheta_2} at momentum Q_2 = n Q_1 mod (2 pi, pi Z). Then the conversion vertex

    O_conv = e^{i (vartheta_2 - n vartheta_1)}                                                     (1.1)

is translation invariant, and under Theta (or the K-stratum antiunitary) it maps to itself up to a phase that can be absorbed by a constant shift, so a real combination is allowed. Its dimension in the decoupled two-boson theory is

    x_conv = n^2/(4 K_1) + 1/(4 K_2).                                                             (1.2)

For an attractive host (K_1 > 1) and n = 2 this is below 1 + 1/(4 K_2), so O_conv is relevant unless the composite fluid is extremely repulsive (K_2 < 1/4, beyond the crystallization point of a commensurate composite fluid). A relevant O_conv locks vartheta_2 = n vartheta_1 and leaves one boson (c = 1). So a second boson made from bound states of the host's own quanta is excluded in attractive liquids. The microscopic origin is the ordinary pair-single exchange psi_2^dagger psi_1 psi_1, allowed by every symmetry, and attraction makes the phase vertices strongly relevant.

Proof sketch. Translation multiplies e^{i vartheta_2} by e^{i Q_2} and e^{-i n vartheta_1} by e^{-i n Q_1}, so O_conv is invariant when Q_2 = n Q_1 mod 2 pi. If Q_2 = n Q_1 + pi, use O_conv cos 2phi_1 (dimension x_conv + K_1, irrelevant) at first order and O_conv^2 (dimension n^2/K_1 + 1/K_2) at second order, again relevant for attractive hosts unless K_2 is small. Theta acts on both bosons with phi reversed (Lemma A), hence vartheta_i -> vartheta_i + c_i as fields, and a coupling g O_conv + h.c. is invariant if and only if c_2 = n c_1 mod 2 pi. For a branch made of n quanta of the host this holds, because the lattice operator that creates the composite is the product of n host creation operators and Theta acts on it as on that product. So g cos(vartheta_2 - n vartheta_1 + arg g) is allowed with generic g. The dimension (1.2) holds at zero inter-species coupling, and a density-density coupling changes it continuously, which does not affect relevance for K_1 well above 1. QED.

Consequence. [derived] Within the U(1) class and to every finite order of perturbation theory around it, every soft wavevector lies in Lambda (1.2 and section 2.1), so by Lemma L no second boson can form. In T and K_I there is in addition no field, hence no finite-density composite fluid to start from.

### 1.4 The Ising modes

Two kinds of real order parameter are available in T and K_I (C159 table 4.1).

(a) Uniform and odd under the antiunitary (a moment along m). [derived] The Theta-odd uniform density of the liquid is sigma^m ~ (2/pi) d_x phi itself. An Ising field sigma with these quantum numbers couples through sigma d_x phi (dimension 9/8, spin +-1 but allowed in pairs with sigma d_x vartheta as in C159 4.1), so Ising x boson needs two tunings (C159). On the attractive side this is the natural instability. The m-coupling of an attractive parent is ferromagnetic (s_m < 0 in the original frame for either sign of s_perp, because Delta = s_m/|s_perp| < 0), so a uniform tilt out of the spiral plane lowers the m-energy. The classical conical tilt becomes unstable only at Delta = -1, and by C160 M3 the tilt is the canonical partner of the phason (solid-angle pairing), not an independent mode. So the soft "Ising" direction that attraction favours belongs to the boson and gives z = 2 or first-order transitions (C161 2c found exactly this, a tilted CSL with c = 1).

(b) At momentum pi, Theta-odd (staggered moment along m) or Theta-even (dimerization). [derived, C159 4.2] The coupling sigma V_c with V_c = cos 2phi or sin 2phi has dimension 1/8 + K and is irrelevant only for K > 15/8. For 15/8 < K < 2 the ordered side is pinned (<sigma> V_c has dimension K < 2), so a c = 3/2 hypersurface there would separate the liquid from gapped Neel (or dimer) order. [heuristic, energetics] In the window the parent's m-coupling is ferromagnetic with |s_m| > 0.67 |s_perp|. Classically a staggered m component costs |s_m| per bond relative to the in-plane spiral, and dimerization needs frustration, which a nearest-neighbour bond does not have (its only route to dimerization is the umklapp cos 4phi, relevant for K < 1/2). So an independent staggered Ising mode has to be created by harmonics of the size of the exchange, which is the regime of 1.1. This is a heuristic statement, not a derivation, and it is the first place where the narrowest assumption bites (section 4).

[derived, continuity remark] A c = 3/2 segment of the liquid-Neel boundary at K > 15/8 cannot join continuously a Kosterlitz-Thouless segment of the same boundary (K = 1/2 there), because the boson's K is continuous along a line of fixed points. So such a segment would be an island bounded by first-order or multicritical points, which is the kind of object assumption (H) of C154 already isolates.

### 1.5 The two-boson option

[derived] By Lemma L, two half-floating bosons in W need two independent incommensurate wavevectors, beta_2 not in Q beta_1 + Q pi. The microscopic bond supplies only theta (harmonics carry q theta), and the parent bands supply only Lambda (1.2). An independent wavevector must therefore be generated dynamically, by a composite whose band has at least two harmonics and an interior minimum. C161 2b gives the only known example in this programme, a soliton band with minima at P = pi -+ 0.443 pi on a Z2-symmetric K background near h_F. There the like-valley attraction makes the onset first order into a tilted CSL with c = 1, and C162 works on soliton onsets in general. For nearest-neighbour magnons around polarized references such interior minima are excluded (C160 Lemma 3', proved). In T the second boson must in addition be anomaly free (C159 4.1 (iv)), and if its V_c sits at momentum 0 it needs x_2 > 2.

## 2. Perturbation theory around the attractive parents

### 2.1 All orders in the harmonic amplitude

[derived; inputs: C154 section 2 (spiral frame), momentum conservation, standard multiscale RG] Write the bond as parent + a sum_q H_q, with H_q of charge q carrying momentum q theta in the spiral frame. At order k in a the generated operators have charge Q = q_1 + ... + q_k and momentum Q theta (mod 2 pi), and in the parent's boson they are vertices e^{i(Q vartheta + 2 n phi)} times polynomials in derivatives, with momentum Q(theta + pi) + n pi (AF frame). Such a term survives the oscillation average only if

    Q (theta + pi) + n pi = 0  (mod 2 pi),                                                         (2.1)

which for fixed Q singles out theta in (pi/Q) Z. Away from these resonances every generated term is an oscillating vertex, irrelevant, or a renormalization of the charge-0 part (K, v, and the DM twist, which shifts theta). At a resonance the relevant term e^{i Q vartheta} (dimension Q^2/(4K), relevant for Q^2 < 8K, which in an attractive liquid includes Q = 1, 2 and at K > 9/8 also Q = 3) locks vartheta, with a tongue width of order a^{k(Q)/(2 - Q^2/(4K))} where k(Q) = ceil(|Q|/2) is the lowest order that builds charge Q from harmonics of charge +-1, +-2. The phase diagram near the parent is therefore a devil's staircase: c = 1 liquid outside the tongues, gapped clock-locked phases inside (breaking Theta or translation, as LSM requires in T), Pokrovsky-Talapov edges and KT points. Every operator in this expansion acts in the one-boson Hilbert space, so no second sector appears at any finite order. Combined with 1.1, a second sector needs a level crossing of order J, i.e. harmonic amplitudes of order the exchange.

### 2.2 Numerical check of the window's population

[numerical, kest.py, kest_c159.out] A Luttinger estimator independent of momentum assignments is K_S(L) = (L/4) S_mm(2 pi/L), from s^m ~ (2/pi) d_x phi, which gives S_mm(q) -> 2K|q|/pi. On XXZ rings it is exactly 1 at XX for every L and approaches K from below (L = 16 against exact: 1.47/1.50 at Delta = -0.5, 1.87/1.98 at -0.7, 2.82/3.48 at -0.9). Applied to C159's stage-B points (job 17597033, completed 2026-11-08), the attractive TFM and T points (parent Delta from -0.64 to -1.18, harmonic fraction 0.10 to 0.58) have K_S(16) from 1.9 to 7.5, inside the Ising x boson window, and their open-chain DMRG reads c_prof(128) between 0.98 and 1.03 with two-size slopes 6 dS/ln 2 between 0.99 and 1.02 at 64 -> 128 (hi_TFM_60, 62 to 66, dist_TFM_77, 78, hi_T_34, 36, 39, ising_T_41, 43, 46). hi_T_37 (K_S about 7) has slopes 1.42 then 1.15, falling. These sizes probe chains of 128 sites. For large K the healing length of the liquid grows (about 1/sqrt(1 + Delta) lattice spacings), so for hi_T_34 to 37 the probed range is 10 to 30 healing lengths. Note that the small-ring level at relative momentum pi is not a usable K estimator (xpi_c159.py: off by -0.5 in the median even at harmonic fraction below 0.15, because the spiral twist L theta places e^{i vartheta} states near pi).

## 3. The twisted curve J = -R(phi) at phi of order one (Task 2)

### 3.1 Spiral frame

[proved, algebra; numerical identity in c163lib.py] With sigma_j = R(-j phi) tau_j (R = R_m, rotation about m),

    sigma_j^T (-R(phi)) sigma_{j+1} = - tau_j . tau_{j+1},                                         (3.1)

so -R(phi) is h_F after a site-dependent rotation (open 6-site chains agree with h_F to 2e-14 at phi = 1). A perturbation J = -R(phi) + eps X, B = eps b becomes

    H = - sum_j tau_j . tau_{j+1} + eps sum_j sum_q e^{i q j phi} [ tau_j^T X'_q tau_{j+1} + b_q . tau_j ],      (3.2)

with X' = X R(-phi), and X'_q, b_q the U(1)_m charge-q parts (q in {0, +-1, +-2} for X', q in {0, +-1} for b), because conjugation by R(j phi) multiplies the charge-q part by e^{i q j phi}.

### 3.2 The charge-zero part is the U(1) class near h_F

[proved, algebra] The q = 0 part of (3.2) is XXZ + DM along m + field along m, a member of the U(1) class at distance eps from h_F. The DM part is a further gauge twist delta = O(eps), which shifts phi to phi + delta, so phi below means the total twist. By the U(1) results (Bethe ansatz) this part has c <= 1 everywhere: a gapped polarized ferromagnet (easy-axis or field-dominated), an easy-plane liquid with K of order eps^{-1/2} (Delta_AF = -1 + O(eps)), and z = 2 boundaries.

### 3.3 The non-axial parts average out except at commensurate twists

[derived; inputs: momentum conservation, Schrieffer-Wolff with the excitation gap of the q = 0 model at momentum q phi] The low-energy states of the q = 0 model have momenta in windows of width O(1/l) (l the healing length, of order eps^{-1/2}) around 0 and pi, shifted by charge times delta. The term of charge q in (3.2) changes the momentum by q phi. Its first-order matrix elements between low-energy states therefore vanish unless q phi = 0 or pi mod 2 pi up to O(eps^{1/2}).

Intermediate states at momentum q phi cost at least E_*(q phi), the lower edge of the q = 0 spectrum at that momentum. In the easy-plane liquid near h_F this lower edge is the out-of-plane soliton branch with E_* of order eps^{1/2} J (its momentum is s Omega by Haldane's Berry phase, so every P in (0, 2 pi) is reached), and in the polarized ferromagnet it is the magnon or q-string branch, of order J (1 - cos q phi)/q. So the second-order terms are of order eps^2/E_* <= eps^{3/2}, small compared with the eps-size q = 0 couplings. They are of two kinds: charge-0 renormalizations (q + q' = 0), and charge Q = q + q' terms at momentum Q phi, which are again oscillating unless Q phi = 0 mod 2 pi. The same holds at every order. The resonant terms are the clock operators

    eps^{k(m)} cos(m vartheta + const)  at  phi = 2 pi p/m,                                        (3.3)

with k(m) = ceil(m/2) (charges +-2 at first order, so Z_3 and Z_4 appear already at second order), plus their n = 1 relatives e^{i m vartheta} cos 2phi at phi = (2p + 1) pi/m, of dimension m^2/(4K) + K, irrelevant at large K. For phi/2 pi irrational no resonant term exists and the infrared has an emergent U(1)_m. Near phi = 2 pi p/m with small detuning the clock term oscillates slowly, which is the commensurate-incommensurate problem of one boson.

[derived] Endpoints. At phi -> 0 the q = +-1, +-2 terms become uniform at first order and the problem is C160's blow-up. At phi = pi the q = +-2 part is uniform at first order (in-plane anisotropy), the q = +-1 part is staggered, and -R(pi) = XX + YY - ZZ is the staggered image of h_F, so the staggered rotation maps this end of gamma onto C160's blow-up as well (C161 1c).

### 3.4 No second sector near gamma

Proposition G. [derived, conditional only on the asymptotic validity of the Schrieffer-Wolff expansion in 3.3, the standard input for spatially oscillating perturbations] For every closed interval I of phi inside (0, pi) there is eps_I > 0 such that every conformal bond J = -R(phi) + eps X, B = eps b with phi in I, 0 < eps < eps_I has c <= 1. Its infrared is one of: a polarized ferromagnet (c = 0), the easy-plane liquid of the q = 0 model (c = 1, K of order eps^{-1/2}), a Z_m clock-locked phase inside a tongue around phi = 2 pi p/m (c = 0, broken Theta or translation), or the KT, Pokrovsky-Talapov and z = 2 boundaries between these. Together with C160 at the two ends of gamma, the whole curve gamma has a neighbourhood with c <= 1, away from C161's soliton-onset layers (C162).

Proof. By 3.2 and 3.3 the effective low-energy Hamiltonian is the q = 0 model (c <= 1 by Bethe ansatz) plus single-boson vertices. A second sector needs an excitation outside the q = 0 model's low-energy sector to reach zero energy. The candidates are the q-string and magnon branches at momentum q phi (gap of order J (1 - cos q phi)/q) and the soliton branch (gap of order eps^{1/2} J, and part of the q = 0 model's own spectrum). Corrections of order eps^{3/2} cannot close these gaps for small eps. At a resonance the softening "mode" is the uniform charge-m state, whose operator is the boson's own e^{i m vartheta} (3.3), so it locks the boson instead of adding a sector. The clock terms are relevant (m^2/(4K) < 2 since K is of order eps^{-1/2}), so the tongues are gapped. QED.

Consistency checks. [derived] (i) At D = 0 and B = 0, J = -R(phi) + eps X with X in the XYZ form reproduces Baxter (the q = +-2 part at phi = pi gives the zero-field XYZ near -1, critical only when two couplings coincide). (ii) For irrational phi the emergent U(1) agrees with C161's heuristic. (iii) The Z_3 resonance at phi = 2 pi/3 appears at order eps^2, consistent with the cyclic-stratum picture of C154 (period-3 structure, A_4 symmetry) in the special case where X has the cyclic form. [pending numerical] twisted.py runs 16 bonds (phi = 1, 2, 2 pi/3, pi/2; eps = 0.15, 0.3; types T and N) at L = 64, 128 in job 17600999. Prediction: c_prof and c_two at or below 1, with slow approach for the easy-plane liquids (K of order eps^{-1/2}), and possibly c read near 1 at L = 128 inside the Z_3 and Z_4 tongues, whose gaps scale like eps^{2/(2 - m^2/(4K))} and can exceed these sizes.

## 4. Numerics on the critical path (Task 3)

### 4.1 The completed C159 job 17597033

[numerical, c159_summary.txt, c159_refine/] Stage A covered 64,000 ring points. Stage B refined 78 points by DMRG at L = 32, 64, 128 (chi 128, truncation 1e-20 to 5e-9). Against C159's pre-fixed criterion (c_two and the profile fits above 1.2 at L = 128 with truncation below 1e-7):
- dist_T_50 meets it (c_prof 1.266, c_two 1.266, slopes 1.20 then 1.34, rising). It is within about 0.2 of the XZ point (singular values of J in ratio 1 : 0.20 : 0.04, leading singular vectors orthogonal to 0.04) and K_S(16) = 0.65, so it is a repulsive liquid near the KT value inside the XZ window that C152 and C158 diagnosed as a slow crossover. It is outside W.
- hi_T_40 (slopes 1.08 then 1.21, c_prof(128) 1.15) does not meet it yet. Its K_S is 0.42 to 0.47, again near the KT value, also outside W.
- ising_TFM_74 meets it formally (c_two 1.76, c_prof 2.72), but S(1) = ln 2 at every L and S_mid is flat at L = 32, 64 (0.951, 0.970) and jumps at L = 128. This is a Theta-broken ferromagnetic cat with a domain wall in the L = 128 run, an artifact, not a critical state. ising_K_14 and 15 have c_prof up to 3.6 but S_mid slopes of 0.66 to 0.78 with large incommensurate oscillations, so their profile fits are unreliable and their entropies grow slower than c = 1.
- Everything else reads c <= 1 at L = 128 or falls with L, and every point in W reads 1 (section 2.2).
So C159's job found no c > 1 in T, K or the ferromagnetic and antiferromagnetic neighbourhoods, and none in W.

### 4.2 Job 17600999 (submitted by this run, pending)

[pending] One node, 48 cores, 8 h, /scratch/lustre/jkudlerflam/c163/. Tasks (typhon/tasks.txt):
- big/: DMRG at L = 256 (chi 192) for hi_T_40, dist_T_50, hi_T_37, hi_K_1, at L = 192 (chi 160) for hi_T_40, dist_T_50 and two 0.03 displacements of hi_T_40, and L = 128 for those displacements. The two-size slopes are printed by summarize.py together with C159's L = 32, 64, 128.
- fes/T40.json, fes/T50.json: two-site-cell damped VUMPS, chi 8 to 192.
- twisted/: section 3.4.
- paths/: 40 paths from attractive parents into T (16, Delta = -0.95 to -0.67, parent K from 4.9 to 1.9), K_I (12, same Delta, with in-plane anisotropy and in-plane field) and N (12, Delta = -0.85 to -0.3, parent field 0.15 to 0.85 of saturation). Gap-first targeting in the form appropriate to an Ising sector inside a liquid (whose total gap stays 1/L): an Ising mode that becomes critical adds (1/12) ln min(xi_I, L) to S_mid, a bump in S_mid(lambda) at fixed L = 48 of width of order 1/L. Each path is scanned on 16 values of lambda in [0, 1.2] (harmonic norm relative to the parent), the two most prominent bumps (excess over linear interpolation of neighbours) and the c_prof maximum are zoomed twice, and the best point gets DMRG at L = 64, 128 with chi 128 and K_S(16).
- Criteria fixed in advance. A c > 1 candidate needs c_two and c_prof above 1.2 at L = 128 with truncation below 1e-7 and K_S(16) > 15/8 (Ising window) or > 1 (two-boson window), then FES at chi >= 256. For hi_T_40 and dist_T_50 the decisive numbers are the slopes 128 -> 192 -> 256: falling towards 1 means crossover, saturating above 1.2 together with an FES slope above 1.2 means a c > 1 candidate outside W, which would contradict C159 Proposition C at K near 1/2 and would have to be examined first for a hidden symmetry.
- Length scales. L = 48 at the scan stage resolves bumps from Ising lengths below about 50 sites. The final DMRG reaches 128 sites (256 for the four refined points), that is 10 to 30 healing lengths for K between 2 and 7.

C155's job 17594766 (running, results folders empty at the time of writing) covers adversarial searches in N, and C160's job 17599583 (running) covers the CSL family at the K -> infinity end. This run does not duplicate them.

## 5. Theorem-style summary (Task 4)

- Fact 1 [proved, C159]: in T and K_I the spiral parent is at zero field, and a connected half-floating liquid has zero magnetization along m.
- Fact 2 [literature, derived]: the zero-field attractive XXZ parent has only the c = 1 Gaussian tower below the lattice scale. Its strings are massless breathers of the same boson.
- Fact 3 [proved for the parent bands, literature for c = 1 in the U(1) class]: in the U(1) class all soft wavevectors lie in Z alpha + Z pi (gauge covariance, single-harmonic n-string bands), and no string branch softens while the liquid survives.
- Lemma L [derived]: a second boson whose wavevector is an integer multiple of the host's (mod pi) is locked by the conversion vertex, relevant in attractive liquids. So bound-state second fluids are excluded in W.
- Proposition P [derived]: to every order in the harmonic amplitude around an attractive parent, the infrared is one boson with a devil's staircase of clock-locked tongues at theta in pi Q. A second sector needs a level crossing of order the exchange.
- Lemma U [derived, C159 table and C160 M3]: the uniform Theta-odd Ising mode, which attraction favours, is the boson's own density channel, and Ising x boson through it has codimension 2.
- Proposition G [derived]: c <= 1 on a neighbourhood of every compact piece of the twisted curve gamma at phi of order one. With C160 at the ends, all of gamma is covered away from soliton-onset layers.
- Numerics [numerical]: 16 liquids with measured K_S between 1.9 and 7.5 (inside the Ising x boson window) read c = 0.98 to 1.03 at L = 128. The only C159 point that meets the pre-fixed c > 1 criterion is a repulsive (K near 1/2) point in the XZ window, outside W, now refined to L = 256 and FES in job 17600999 together with 40 targeted paths through W.

Verdict. c <= 1 is not derived on all of W. It is derived near the U(1) class (all orders, Proposition P), near gamma (Proposition G) and in the K -> infinity corner (C160), and the candidate mechanisms that rely on bound states or on the uniform tilt are excluded (Lemma L, Lemma U). No c > 1 bond was found. The narrowest remaining assumption is

    (A_W) [assumed] In an attractive half-floating liquid of T or K_I (K > 1), with harmonic amplitude of the order of the exchange, no lattice-scale excitation outside the boson's charge lattice reaches zero energy while the liquid survives. Concretely, (a) no staggered (momentum pi) Ising mode independent of e^{+-2 i phi} softens when K > 15/8, and (b) no composite develops a band minimum at a wavevector outside Z theta + Z pi that reaches zero energy when K > 1.

Where it bites. (a) Parents with -1 < Delta < -0.67 and harmonic fraction above roughly 0.3 (below that, Proposition P applies at the physics standard), away from the semiclassical neighbourhoods of h_F and gamma. (b) The same with -1 < Delta < 0, and in particular at soliton onsets on Z_2-symmetric backgrounds (C161 2b), which C162 treats. In N the Ising x floating boson option is codimension 1 at every K (C154), so A_W is not the relevant assumption there, and N remains under C154's (A6).

## 6. What is not done

- The numerical tests of section 4.2 are queued (job 17600999). Until they report, the numerical support for W is the 16 C159 points of section 2.2.
- Lemma L assumes that the composite fluid, if it formed, would be a Luttinger liquid with K_2 above about 1/4. A strongly repulsive composite crystal is not covered, but it would be gapped.
- The energetic argument against a staggered Ising mode (1.4 (b)) is heuristic.
- The Schrieffer-Wolff step of Proposition G is asymptotic. eps_I is not quantified and shrinks as I approaches 0 or pi, where C160 takes over with its own unquantified eps0.
