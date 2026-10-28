# The explicit route-A candidate: one bond on ℂ²³ with c = n

## Claim

For every integer n ≥ 1 there is a Hermitian nearest-neighbour bond h_n on ℂ²³ ⊗ ℂ²³, built from the same 23 letters and the same list of terms for every n, with only the coefficients depending on n. The chain H_n = Σ_i h_n(i, i+1) on a ring is one-site translation invariant. Its ground states break translation spontaneously into frozen cells of length b = n + 1, each holding n data qubits and one marker. Each cell carries an exactly degenerate 2ⁿ-dimensional code, neighbouring cells interact only at second order in a port coupling λ, and the resulting logical chain is n bands of complex fermions with exact U(1) and fermion parity, every band crossing a common logical chemical potential, with distinct Fermi-sea centres. The claim is that on rings whose length is a multiple of n + 1 the low-energy theory has c = n. Since d − 1 = 22, the first member that would beat c ≤ d − 1 is n = 23 (c = 23, cell length 24). As it stands, the bond and its finite-dimensional properties are explicit and checked, and the logical bands with distinct centres are verified numerically through n = 3. Positivity of the leading hopping matrix, which makes every band cross a common field, is checked here through n = 6, certified in integer interval arithmetic at n = 23 for one parameter point by the running Astra C33, and proved for every n with an explicit schedule for the mobile variant of Astra C32. Distinct centres with non-resonant momentum margins are established only through n = 3, and the passage from bands to c = n is standard bosonization input. The cells are frozen, so the chain has c = n on top of a spontaneously broken translation, not one CFT with a unique ground state on every ring.

## Status

Heuristic as a whole. Every finite-dimensional statement about the bond is proved or checked by exact computation, and the bond has been re-implemented here from the text below and matched entry by entry to the repository code. The band structure with distinct centres is numerical for n ≤ 3. Positivity of the leading hopping is numerical for n ≤ 6 here and certified at n = 23 by C33 (in progress, unreviewed). No central charge has been measured for this chain. The upgrades that matter are, in order, distinct centres with non-resonant margins for every n up to 23, a calibrated measurement of c at n = 2, and either an argument that frozen cells are admissible or the mobile version (Section 6).

## How it was established and checked

- Author: Claude (Opus 5.5) subagent for Jonah Kudler-Flam, 2026-10-28. This folder assembles the 23-letter chain of `contrib/routeA_spatial/` (a Claude subagent building on Astra C22) with the selector of `codex/work_C22/`, the stability inputs of `contrib/ir_control/`, `codex/work_C13/` and `codex/work_C20/`, and the mobile evidence of `contrib/mobile_cells/`. No Codex call was made for this folder.
- New checks in this folder. `bond23.py` builds the literal 529 × 529 bond from Section 1 of this text, independently of `ycell.py`, and compares it with the repository. `extend_T.py` and `gap_delta.py` compute the second-order hopping by a method with no scalar background, extend the positivity check to n = 5, 6, and regenerate the published logical-chain numbers.
- Astra results folded in: C30 and C31 (finished), C32 (finished at 07:50 EDT), C33 (PROGRESS and certificate files as of 08:00 EDT, run still going). Their claims are quoted with their own tags and were not re-derived here.
- Checked by: nobody other than the AI authors. The ED of `contrib/routeA_spatial/` is an internal check by its author. Not refereed.

## 1. The bond

### 1.1 Alphabet

Twenty-three letters. Four carry only a data bit and a position relative to the token, fourteen carry the token on a data site together with that site's bit, and five are markers.

| Letters | Count | Meaning | Arrows (in, out) | Charge q |
| :-- | --: | :-- | :-- | :-- |
| L_a, a = 0, 1 | 2 | data site left of the token, data bit a | (L, L) | a |
| R_a | 2 | data site right of the token | (R, R) | a |
| Xd_{s,a}, s ∈ {v, 0, +, −} | 8 | token of species s on the data track, sitting on a data site with bit a | (L, R) | a + q_s |
| Xx_{s,a}, s ∈ {0, +, −} | 6 | token on the express track, on a data site with bit a | (L, R) | a + q_s |
| M_∅ | 1 | empty marker | (R, L) | 0 |
| M_s, s ∈ {v, 0, +, −} | 4 | marker holding the token | (L, L) | q_s |

The species charges are q_v = q_0 = 0, q_± = ±1. The species v is the vacuum token, 0, + and − are excited. The vacuum token never enters the express track. The basis order used in `bond23.py` (and `ycell.py`) is the order of the table, L_0, L_1, R_0, R_1, Xd_{v,0}, Xd_{v,1}, Xd_{0,0}, …, Xd_{−,1}, Xx_{0,0}, …, Xx_{−,1}, M_∅, M_v, M_0, M_+, M_−.

### 1.2 The terms

The bond acts on (site i, site i+1), left factor first. Write z_s(a) = (−1)^a for s = 0 and z_s(a) = 1 for s = +, −, so the neutral token carries a Pauli Z of the data bit it departs from. In every line a, b ∈ {0, 1}, s runs over the excited species {0, +, −}, and every term is added together with its Hermitian conjugate.

```
h = h_kin + K Π_bad + λ V + ½ (o ⊗ 1 + 1 ⊗ o)

h_kin =
 (K1)  −t                     |L_a Xd_{v,b}⟩⟨Xd_{v,a} R_b|        vacuum token, data track
 (K2)  −t   e^{iφ_s}   z_s(a) |L_a Xd_{s,b}⟩⟨Xd_{s,a} R_b|        excited token, data track ("stub")
 (K3)  −t_x e^{iφ^x_s} z_s(a) |L_a Xx_{s,b}⟩⟨Xx_{s,a} R_b|        excited token, express track
 (K4)  −t                     |L_a M_v⟩⟨Xd_{v,a} M_∅|            vacuum token into the marker
 (K5)  −t_x e^{iφ^x_s} z_s(a) |L_a M_s⟩⟨Xx_{s,a} M_∅|            express track into the marker
 (K6)  −t_d e^{iφ^d_s} z_s(a) |L_a M_s⟩⟨Xd_{s,a} M_∅|            weak loop link, stub into the marker
 (K7)  −t_sw e^{iφ^sw_s}      |m Xx_{s,a}⟩⟨m Xd_{s,a}|            junction switch, for each m ∈ {M_∅, M_v, M_0, M_+, M_−}

K Π_bad:  K on every ordered pair (l, r) whose out-arrow of l differs from the in-arrow of r

λ V = λ [ F_0^M ⊗ F_0^X + e^{iφ} F_+^M ⊗ F_−^X + e^{−iφ} F_−^M ⊗ F_+^X ]
      F_0^M = |M_0⟩⟨M_v| + |M_v⟩⟨M_0|,          F_+^M = |M_+⟩⟨M_v| + |M_v⟩⟨M_−|,          F_−^M = (F_+^M)†
      F_0^X = Σ_a (|Xd_{0,a}⟩⟨Xd_{v,a}| + h.c.), F_+^X = Σ_a (|Xd_{+,a}⟩⟨Xd_{v,a}| + |Xd_{v,a}⟩⟨Xd_{−,a}|), F_−^X = (F_+^X)†
      (F^M vanishes on M_∅ and on non-marker letters, F^X vanishes on every letter other than Xd)

o = ε Σ_{letters ℓ with data bit a = 1} |ℓ⟩⟨ℓ|
  + Σ_{s ∈ 0,+,−} (Δ + δ_s) [ Σ_a |Xd_{s,a}⟩⟨Xd_{s,a}| + Σ_a |Xx_{s,a}⟩⟨Xx_{s,a}| + |M_s⟩⟨M_s| ],   δ_0 = 0, δ_± = ±δ
  + D_x Σ_{s,a} |Xx_{s,a}⟩⟨Xx_{s,a}| + D_M Σ_{s ∈ 0,+,−} |M_s⟩⟨M_s|
  + μ_sel Σ_{m ∈ M_∅, M_v, M_0, M_+, M_−} |m⟩⟨m|
  + g_+ (|Xd_{+,0}⟩⟨Xd_{0,1}| + h.c.) + g_− (|Xd_{−,1}⟩⟨Xd_{0,0}| + h.c.)
```

In words: the token moves right on the data track by (K1, K2) and into the marker by (K4) if it is the vacuum species. Excited tokens cannot enter the marker from the stub except through the weak link (K6). They reach the marker by switching to the express track at the first site of a cell (K7, the only move whose left letter is a marker) and running along it (K3, K5). The conversions g_± turn a neutral token into a charged one while flipping the data bit under it, conserving q. The port V acts on (marker of one cell, first site of the next), turning both vacuum tokens into excited ones. The Hamiltonian of a ring of L sites is H = Σ_{i=1}^{L} h(i, i+1), so the on-site operator o is counted once per site.

### 1.3 Coefficients

All couplings are real and the phases are angles. The values used in the checks are

| Parameter | Role | Checked value | Schedule in n (Section 3) |
| :-- | :-- | :-- | :-- |
| t, t_x, t_sw | hopping on data track, express track, switch | 1, 1, 1 | unit of energy |
| K | penalty on illegal pairs | not needed inside the legal sector (10 in `bond23.py`) | K ≥ 8t + margin |
| μ_sel = 2t − δ_sel | marker chemical potential, selects b = n + 1 | 0 inside fixed-marker computations | δ_sel at the midpoint of (δ_{n+1}, δ_n), Section 2.3 |
| Δ | excited energy on every excited letter | 0 | Δ ≥ \|g\| with margin, for example 0.1 (C33 uses Δ = 4 − 2cos(π/(n+2)) ≈ 2) |
| δ | ± splitting (particle-hole asymmetry) | 0 | 0 |
| D_x, D_M | extra energy on express and excited-marker letters | 4 (n ≤ 3), 10 (n = 4) | increasing, ≥ 20 at n = 5, 6, and 1000 at n = 23 with Δ ≈ 2 (C33) |
| g_+ = g_− = g | conversion | 0.05 (bands), 0.01 (infrared estimate at n = 2) | decreasing, Section 3 |
| φ_s, φ^sw_s | stub and switch phases | 0 | 0 |
| φ^x_0, φ^x_± | express phases | 0.4, 0 | φ^x_0 = π/(2n), φ^x_± = 0 |
| t_d, φ^d_0, φ^d_± | loop link and its phases | 10⁻³ (n = 1, 2), 10⁻⁵ (n = 3); 0.9; 0 | decreasing |
| φ | port phase | π/2 | π/2 |
| ε | data field = logical chemical potential | λ² · ½ λ_min(R_s) (Section 2.4) | same rule |
| λ | port coupling | 0.05 to 0.2 in ED | last and smallest, Section 3 |

The express phase is the one coefficient where the checked value should not be extrapolated. The leading hopping is proportional to sin(nφ^x_0) (Section 2.4), and with φ^x_0 = 0.4 this nearly vanishes at n = 8 (|sin 3.2| = 0.058) and is 0.22 at n = 23.

### 1.4 Cross-check of the text against the repository

[checked, `bond23.py`, `bond23_results.json`] The 529 × 529 matrix built from the list above is Hermitian, commutes with q ⊗ 1 + 1 ⊗ q and with Π_bad, at the working point and at a generic point with all phases, t_d, δ, ε and μ_sel nonzero. It has 795 nonzero entries at the working point with t_d = 0, λ = 1, K = 10. Restricted to a frozen M_∅ followed by one legal cell it equals `ycell.build_cell(n).H` entry by entry (difference 0 in floating point, n = 1 to 4, generic point with g_+ ≠ g_−), and the port operators agree exactly. On the complete open-chain spaces of 23³ = 12167 and 23⁴ = 279841 states (n = 1, 2) the legal-cell subspace is invariant, with zero leakage, and the restricted matrix agrees with the cell builder. The legal-cell dimensions are (7n + 4)2ⁿ = 22, 72, 200, 512.

## 2. What the bond does

### 2.1 Legal cells and the code

[proved, C22 §2 and `contrib/routeA_spatial/` §2] Every move preserves the outer arrows of its pair, so Π_bad commutes with H, and marker positions are exactly conserved. Between two markers a zero-penalty word with n data sites is one of

```
L^j Xd_s R^(n−j−1) M_∅   (s = v, 0, +, −),   L^j Xx_s R^(n−j−1) M_∅   (s = 0, +, −),   L^n M_s   (s = v, 0, +, −),
```

with j = 0, …, n − 1 and an arbitrary data string, so each cell holds one token and the legal-cell space has dimension (7n + 4)2ⁿ.

[proved] H0 (everything but λV) never changes vacuum into excited species. In the vacuum sector the token walks on a path of b = n + 1 vertices (data sites 0, …, n − 1, then the marker) with hopping −t, acting as the identity on the data. At ε = 0 its ground energy is E0 = μ_sel − 2t cos(π/(n + 2)) for every data string, and its wavefunction is

```
G|a⟩ = Σ_{j=0}^{n} v_j |vacuum token at position j⟩ ⊗ |a⟩,   v_j = √(2/(n+2)) sin(π(j+1)/(n+2)),
```

for every data string a. The endpoint amplitude is α_n = v_0 = v_n = √(2/(n+2)) sin(π/(n+2)). The field ε adds ε N_data on the code, an exact logical chemical potential.

[numerical, `gap_delta_results.json`] The code is the complete ground space of the cell at every tested point. The vacuum-walker gap equals ω_b = 2t[cos(π/(b+1)) − cos(2π/(b+1))] to machine precision. The excited-token gap at Δ = 0, g = 0.05, D_x = 4 is 0.70, 0.27, 0.12, 0.058, 0.027, 0.011 for n = 1 to 6. It falls roughly like n⁻³ and is lowered by about g/2 by the conversions (n = 6, D_x = 20: 0.019 at g = 0.05, 0.039 at g = 0.01). With Δ = 0.1 it is 0.12 at n = 6. [heuristic] This is why the schedule needs Δ ≥ |g| with a margin: at Δ = 0 and fixed g, long cells would put an excited token below the code. C22's proof that excited states lie Δ − |g| above the code does not carry over, because the excited graph here is the folded path of 2n + 1 vertices, not the vacuum path.

### 2.2 Symmetries

[proved] Every term conserves Q = Σ_sites q, which on the code is the logical fermion number Σ_i n_i. The on-site unitary Π_sites (−1)^{q} is the logical parity P = Π_i Z_i on the code. Exact U(1) forbids pairing terms at every order, which removes the instability seen in `contrib/routeA_numerics/`. The neutral hops carry Z strings, so after C22's gauge U_j = Z_0 ⋯ Z_{j−1} the conversions become Jordan–Wigner fermions c_j and the right port sees the string P through a local wire, not through a nonlocal term.

[proved] Every F changes vacuum to excited, so G†F G = 0 and first order vanishes. An outer cell touched once by λV is left excited and cannot return under H0, so there is no three-cell term at second order. Conjugation by (−1)^{N_exc} on alternate cells sends λ → −λ and is the identity on the code, so odd orders vanish on open chains and on rings with an even number of cells.

### 2.3 Selection of the cell length

[proved for the vacuum sector, C22 §3] With μ_sel = 2t − δ_sel, a_b = 2t[1 − cos(π/(b+1))] and δ_b = (b+1)a_b − b a_{b+1}, the energy density e_b = (a_b − δ_sel)/b of a legal ring is minimized uniquely at length b when δ_b < δ_sel < δ_{b−1}, with margin γ_b = min_{k≠b} k(e_k − e_b) = Θ(t b⁻⁴) and window width Θ(t b⁻³). Midpoint values, from these formulas:

| n | b | window (δ_b, δ_{b−1}) | μ_sel | γ_b | ω_b | α_n |
| --: | --: | :-- | --: | --: | --: | --: |
| 1 | 2 | (1.828, 3.0), use 1.9 | 0.100 | 3.6e-2 | 2.00 | 0.707 |
| 2 | 3 | (1.197, 1.828) | 0.487 | 0.105 | 1.41 | 0.500 |
| 3 | 4 | (0.838, 1.197) | 0.982 | 4.5e-2 | 1.00 | 0.372 |
| 4 | 5 | (0.617, 0.838) | 1.272 | 2.2e-2 | 0.73 | 0.289 |
| 23 | 24 | (0.0443, 0.0479) | 1.954 | 7.6e-5 | 0.047 | 0.035 |

[heuristic for this alphabet] The C22 statements that excited sectors cannot undercut the vacuum cell at any length, and that illegal sectors cost at least K − 4t, were proved for the 17-letter alphabet. Here they need the excited gap at every length b′ (Section 2.1, hence Δ ≥ |g|) and a hopping bound that includes the express track, the switch and t_d. Neither has been re-derived. The second-order energies, about ½λ²α_n⁴ per bond, must also stay below γ_b, which is automatic for the λ used below.

### 2.4 Second order and the logical chain

[proved] The exact λ² coefficient of one bond on code ⊗ code is

```
H2 = − P V Q (H0 − 2E0)^(−1) Q V P,
```

and the logical chain to this order is

```
H_log = λ² Σ_j [ c_j† T c_{j+1} + h.c. + c_j† μ c_j ] + ε Σ_j c_j† c_j + λ² (non-Gaussian part, O(g⁴)) + O(λ⁴),
```

with c_j the n-component fermion of cell j, Jordan–Wigner ordered by cell and then by data site, and T_xy = ⟨A_x|H2|B_y⟩ the one-particle matrix element moving a fermion from data site y of the right cell to data site x of the left cell.

[proved structure, `contrib/routeA_spatial/` §3] At t_d = 0 the excited token graph of a cell is a path M, x(n−1), …, x(0), d(0), …, d(n−1). The right port of a cell sits at M and the left port at d(0), so both ports see data site j at graph distance j beyond d(0). On a tree the resolvent from M to d(j) factorizes through d(0) with a flavour-independent transfer factor. As D_x → ∞ that factor becomes energy independent, the right-port filters become a constant times the left-port filters f_j(τ), and T → e^{iΦ} c Γ with Γ_ij = ∫ f_i f_j dτ positive definite (f_j starts at order τ^{2j+1}). Then h(k) = μ + T e^{ik} + T† e^{−ik} equals μ ± 2R_s at k = −Φ and π − Φ, where R_s is the Hermitian part of e^{−iΦ}T. If R_s > 0 and ‖μ + ε‖ < 2λ_min(R_s), every band changes sign on both arcs and there are at least 2n Fermi points. This is the folding that repairs C22, whose two ports saw the flavours in opposite orders.

[numerical, this folder] At t_d = 0, T = g_+g_− T⁽²⁾ + O(g⁴). The coefficient T⁽²⁾ is a phase times a real matrix (imaginary residue after removing the phase below 6 × 10⁻¹⁶ of max|T| for n ≤ 6), and it is proportional to sin(nφ^x_0) when φ^x_± = 0 (ratio of |T⁽²⁾| at φ^x_0 = 0.4 and π/(2n) equals sin(0.4n) to four digits at n = 2, 3, 5, 6). Real express rails give no hopping, which is C22's sine cancellation. Hence the schedule φ^x_0 = π/(2n).

[numerical] Smallest eigenvalue of R_s divided by max|T|, from T⁽²⁾ at Δ = 0, t_d = 0 (`extend_T_results.json`; the n ≤ 4 entries at D_x = 4, 10 agree with `hermpart.json` within the O(g²) difference between T⁽²⁾ and the full-in-g T). The largest eigenvalue is between 1.2 and 2.9 in these units.

| n | D_x = 4 | D_x = 10 | D_x = 20 | D_x = 40 |
| --: | --: | --: | --: | --: |
| 2 | 1.8e-2 | 3.0e-2 | 3.5e-2 | 3.8e-2 |
| 3 | 2.2e-4 | 9.2e-4 | 1.4e-3 | 1.6e-3 |
| 4 | −6.9e-6 | 1.4e-5 | 3.3e-5 | 4.3e-5 |
| 5 | −5.3e-5 (two negative) | −5.6e-8 | 8.3e-7 | 1.3e-6 |
| 6 | −1.2e-4 (two negative) | −1.5e-6 | 1.6e-8 | 3.3e-8 |

So R_s is positive definite at D_x = 4 for n ≤ 3, at D_x = 10 for n = 4, and at D_x = 20 for n = 5, 6. The ratio λ_min/λ_max falls by a factor of about 40 per flavour at large D_x, and |T⁽²⁾| by about 10³ per flavour at D_x = 20 (7.2e-10, 4.2e-13, 3.7e-16, 3.2e-19, 3.5e-22 for n = 2 to 6). A small Δ does not change the picture (Δ = 0.1 at D_x = 20 gives smallest eigenvalues 3.4e-2, 1.3e-3, 4.1e-5, 1.1e-6 for n = 2 to 5).

[numerical, precision] The full-in-g extraction used in `contrib/routeA_spatial/` sums an O(1) scalar background, so its absolute floor is about 10⁻¹⁷ in units where |T| is 10⁻¹² to 10⁻¹⁸. Two independent full-in-g methods (the Cauchy form of `ycell.py` and a sparse solve in the one-particle charge sectors here) agree to 3e-12, 1e-11, 3e-10, 1e-8 (relative) for n = 1 to 4 at D_x = 4 and 8e-6 at n = 4, D_x = 10, but only 7e-3 at n = 4, D_x = 20, where they give smallest eigenvalues 4.5e-5 and 3.0e-5. The D_x = 20 row of `hermpart.json` is therefore at the double-precision floor. The coefficient T⁽²⁾ computed as −v_out† A0⁻¹ C A0⁻¹ C A0⁻¹ v_in has no scalar background, agrees with a Richardson extrapolation of the full-in-g T to 6e-4 (n = 3) and 4e-3 (n = 4), and is the source of the table above. With λ_min/λ_max near 40^{−(n−1)}, even this method loses the smallest eigenvalue in double precision around n ≈ 10, so n = 23 needs a theorem or exact arithmetic.

[numerical, this folder] At large D_x the coefficient falls as |T⁽²⁾| ∝ D_x^{−(2n+2)} (doubling D_x from 20 to 40 divides it by 2^{5.9}, 2^{7.7} and 2^{13.4} at n = 2, 3, 6), the square of the transfer through the (n + 1)-vertex express block.

[Astra C33, running, not reviewed; `codex/work_C33/`, read 2026-10-28 08:00 EDT] C33 derives a tree transfer identity, Laplace(f_M,i) = β(z)² Laplace(f_0,i) with β(z) = 1/det(C + z) for the express block C, and certifies with 256-bit integer interval arithmetic (`certify_tree.py`, `certificate_n*.json`) that the Hermitian part of the leading coefficient is positive definite at n = 2, 4, 8, 12, 16 and 23. Its point is stub energy a = 4 above the vacuum (Δ = 4 − 2cos(π/(n+2)) ≈ 2), D_x = D_M = 1000, φ^x_0 = π/(2n), t_d = 0. At n = 23 the certified smallest eigenvalue lies between 10^{−63.6} and 10^{−62.3} against a largest of 10^{−2.35}, and D_x = 100 is not enough (floating-point minimum −3 × 10⁻⁶ relative). Its floating-point scan gives |T⁽²⁾| ≈ 10^{−152} at n = 23, D_x = 1000, in its normalization. A finite-g transfer with explicit g(n) is reported in progress. This point differs from the working point of §1.3 in Δ and D_x.

[Astra C32, finished, not reviewed; `codex/work_C32/RESULT.md` §3] For its 28-letter positive variant (Section 5), C32 proves T⁽²⁾(D) = −c_n D^{−(2n+4)} Γ_n + O(D^{−(2n+5)}), with Γ_n the positive Gram matrix and c_n = 2/(n+2)², and an explicit sufficient schedule D_n, g, t_d, ε under which every band crosses on both arcs, for every n. The bounds are very conservative (log₁₀ D_n from 30 to 79 for n = 1 to 5) and say nothing about simple zeros, distinct centres or pair-sum margins.

### 2.5 Distinct centres

[numerical, `contrib/routeA_spatial/` §4] At t_d = 0 the Hilbert–Schmidt-projected T is a phase times a real matrix to 4e-12 at n = 2 (checked here) and μ vanishes, so all Fermi seas share one centre and inter-flavour Cooper transfer is resonant. Breaking particle-hole symmetry instead (δ ≠ 0) gives |μ|/|T| ≈ 2.4 × 10³ at n = 2, far too large. The weak link t_d closes the loop M → express → d(0) → stub → M with a species-dependent flux (φ^d_0 − nφ^x_0 for the neutral token in this orientation, zero for the charged ones). Then μ = iA with A real antisymmetric and proportional to t_d, the k → −k symmetry of h(k) is broken, and the centres separate linearly in t_d. The coefficient grows with n (|μ|/|T| per unit t_d is 0.57 at n = 2 and 8.2 at n = 3), because the direct link bypasses the express suppression, and too large a t_d makes R_s indefinite (n = 2, t_d = 10⁻² at D_x = 10, Δ = 0.5).

### 2.6 The logical chain at the checked points

[numerical, `contrib/routeA_spatial/logical_n*.json`; the n = 1, 2 rows regenerated here with the repository code to every printed digit, `extend_T_results.json`] Working point of §1.3 with D_x = 4, g = 0.05 and ε = ½λ_min(R_s) (in λ = 1 units). The quadratic part is the Hilbert–Schmidt projection referred to half filling. δ₂ is the smallest Fermi-point separation, δ₄ the smallest |p₁ + p₂ − p₃ − p₄| mod 2π over distinct pairs (it includes inter-flavour Cooper transfer, 2(centre_a − centre_b)), v_min the smallest Fermi velocity.

| n | t_d | cell gap | \|T\| (λ = 1) | λ_min/λ_max of R_s | Fermi points | centres | δ₂ | δ₄ | v_min/\|T\| | non-Gaussian/\|T\| |
| --: | --: | --: | --: | --: | --: | :-- | --: | --: | --: | --: |
| 1 | 10⁻³ | 0.70 | 3.6e-7 | 1 | 2 | 0 | 2.64 | 1.01 | 1.94 | 2.1e-3 |
| 2 | 10⁻³ | 0.27 | 1.1e-8 | 1.4e-2 | 2 + 2 | ±0.0017 | 0.54 | 6.8e-3 | 0.047 | 5.9e-3 |
| 3 | 10⁻⁵ | 0.12 | 1.3e-10 | 1.0e-4 | 2 + 2 + 2 | 0.0013, 0, −0.0014 | 0.68 | 2.6e-3 | 9.5e-4 | 2.2e-2 |

[numerical] The non-Gaussian part scales as g² (2.3e-4 at n = 2 and 8.7e-4 at n = 3 for g = 0.01, band data unchanged). The g⁴ correction to the one-particle hopping is 0.5, 1.5, 5 and 14 percent of T at n = 1 to 4 for g = 0.05 (`bond23_results.json`), so g has to shrink with n.

[numerical, `contrib/routeA_spatial/ed_*.json`] ED of open chains of two and three complete legal cells (dimensions 484, 5184, 10648, 373248) at λ = 0.2, 0.1, 0.05 reproduces M E0 + λ² Σ_bonds H2 with error/λ⁴ constant (0.005 to 0.027), confirming the absence of odd orders and of three-cell terms. The ED hopping matrix converges to T as λ² (relative error between 0.01 λ² and 0.33 λ² at the two parameter points). Eigenpair residuals are below 6 × 10⁻¹⁴ and the retained band is at least 0.27 below the next level. No central charge is measured.

### 2.7 From bands to c = n

[heuristic, verdict (b) of `contrib/ir_control/`] With exact U(1) and parity, n bands with 2n simple Fermi points, distinct centres, δ₂ and δ₄ > 0, positive Gaussian stiffness, and corrections ε_NG with ε_NG ≪ v_min δ₂ and (ε_NG/v_min) log(1/δ₄) ≪ 1, standard bosonization gives n compact bosons, so c = n. The multi-channel Luttinger theory is close to the Mastropietro–Porta framework but no theorem covers it, and the local elimination of the gapped cell sectors (leakage, H6 of ir_control) is assumed. At n = 2 with g = 0.01 the ratio ε_NG/(v_min δ₂) is 0.009. At n = 3 it is 1.35 at g = 0.01, so g ≈ 10⁻³ is needed (ratio about 0.013 and (ε_NG/v_min) log(1/δ₄) ≈ 0.05 by the measured g² scaling, not run). The statement in `contrib/routeA_spatial/` §7 that all conditions are met at n = 3 therefore holds only at this smaller g.

## 3. Parameter schedule in n and the scale hierarchy

[heuristic, fixed n] The couplings are chosen one after another, each as an upper bound in terms of the ones before. No step bounds a smaller scale from below by a larger one, so the choice closes at every fixed n. Nothing needs to be uniform in n, since each h_n is a separate chain.

1. b = n + 1, δ_sel at the midpoint of the selector window (§2.3), K ≥ 8t plus margin.
2. D_x = D_M = D_x(n) large enough that R_s⁽²⁾ > 0 at t_d = 0, with margin r_n = λ_min(R_s)/|T|. Checked values: D_x = 4 (n ≤ 3), 10 (n = 4), 20 (n = 5, 6) at Δ = 0, and D_x = 1000 at n = 23 with Δ ≈ 2 (C33 certificate). C32 gives an explicit, very conservative D_n for its variant, and C33 reports one for this cell in progress.
3. φ^x_0 = π/(2n), φ^x_± = 0, φ = π/2, all other phases zero.
4. Δ ≥ |g| with margin (Δ = 0.1 is checked to leave R_s positive for n ≤ 5), so the code stays the cell ground space at every length.
5. t_d(n) small enough that iA + ε stays below 2λ_min(R_s), yet nonzero, which separates the centres (δ₄ ∝ t_d). Checked: 10⁻³ (n = 2), 10⁻⁵ (n = 3).
6. g(n) small enough that the g⁴ corrections move T and the centres by much less than r_n and δ₄, and that ε_NG ≪ v_min δ₂ and (ε_NG/v_min) log(1/δ₄) ≪ 1. Checked: g = 0.01 at n = 2, about 10⁻³ needed at n = 3.
7. ε = ½λ_min(R_s) in λ = 1 units, so the physical data field is λ²ε.
8. λ last: λ ≪ cell gap and λ²α_n⁴ ≪ γ_b, and the relative λ⁴ correction to T (up to about 0.3 λ² in the ED) ≪ r_n and v_min δ₂/|T|. Estimated: λ ≲ 10⁻² at n = 2.

The resulting energy scales, in decreasing order, with t = 1:

| Scale | Size | n = 2 | n = 3 | n = 23 (C33 point) |
| :-- | :-- | --: | --: | --: |
| t, t_x, t_sw, K, D_x | top | 1 to 10 | 1 to 10 | D_x = 1000 |
| vacuum-walker gap ω_b | ≈ 3π²t/b² | 1.41 | 1.00 | 0.047 |
| excited cell gap | ≈ c n⁻³ − g/2 + Δ | 0.27 | 0.12 | ≈ Δ ≈ 2 |
| selector margin γ_b | Θ(t b⁻⁴) | 0.105 | 0.045 | 7.6e-5 |
| second-order scalar energy per bond | ≈ ½λ²α_n⁴ (measured −0.13, −0.029, −0.0084 λ² at n = 1, 2, 3) | 0.029 λ² | 0.0084 λ² | about 10⁻⁶ λ² |
| logical bandwidth J_n = λ² g² \|T⁽²⁾\| | ∝ D_x^{−(2n+2)}, ~10³ per flavour at D_x = 20 | λ² × 1.1e-8 (g = 0.05) | λ² × 1.3e-10 | about λ²g² × 10⁻¹⁵² |
| narrowest band λ² λ_min(R_s) | ~40 per flavour at Δ = 0, more at Δ ≈ 2 | J × 1.4e-2 | J × 1.0e-4 | about J × 10⁻⁶¹ |
| loop splitting, δ₄ | ∝ t_d | 6.8e-3 | 2.6e-3 | tiny, set by t_d(23) |
| non-Gaussian couplings | ∝ λ²g⁴ | ≪ v_min δ₂ at g = 0.01 | needs g ≈ 10⁻³ | needs g(23) tiny |

[numerical, C33] At the certified n = 23 point the logical bandwidth is about 10⁻¹⁵² g²λ² t, and the narrowest band another 10⁻⁶¹ smaller. That is allowed for existence at fixed n, but it means no simulation will ever see the n = 23 member directly. The violating claim has to rest on an all-n argument plus the infrared theory.

## 4. Status of each step

| Step | Statement | Tag | Evidence, sizes | Source |
| :-- | :-- | :-- | :-- | :-- |
| 1 | Bond Hermitian, exact U(1), arrows preserved, one 529 × 529 matrix from §1.2 | proved; checked | entry-by-entry match with ycell, n = 1–4; complete spaces 12167, 279841 | `bond23.py`, `ycell.py` |
| 2 | Legal cells (7n+4)2ⁿ, one token per cell, markers frozen | proved | structural; dimensions 22, 72, 200, 512 checked | routeA_spatial §2, C22 §2 |
| 3 | Code = 2ⁿ vacuum-walker states, exactly degenerate | proved | residual 2e-16 | ycell.code, `bond23.py` |
| 4 | Code is the complete cell ground space | numerical | n = 1–6, several D_x, g, Δ | `gap_delta_results.json` |
| 5 | Selector picks b = n + 1 with margin Θ(b⁻⁴) | proved for the vacuum sector (C22); heuristic for this alphabet | needs Δ ≥ \|g\| | C22 §3 |
| 6 | Illegal sectors cost ≥ K − O(t) | proved for 17 letters; assumed here | none | C22 §3 |
| 7 | Zero first order, no three-cell λ² term, no odd orders on open chains | proved | ED error ∝ λ⁴ | routeA_spatial §2, §6 |
| 8 | Exact λ² bond H2 | proved; checked by ED | M = 2, 3 cells, n = 1, 2, dims up to 373248 | `ed_*.json` |
| 9 | T⁽²⁾ = phase × sin(nφ^x_0) × real matrix at t_d = 0 | numerical; structural reason (tree) | n = 2–6 | `extend_T_results.json` |
| 10 | Gram limit T → e^{iΦ}cΓ, Γ > 0, and D_x^{−(2n+2)} scaling | proved structure; explicit error bound proved for C32's variant, in progress for this cell | scaling checked n = 2–6 | routeA_spatial §3, C32 §3, C33 |
| 11 | R_s > 0 at finite D_x | numerical; certified at one point (Astra, unreviewed) | n ≤ 3 (D_x = 4), n = 4 (10), n = 5, 6 (20), two full-in-g methods for n ≤ 4; certified n = 2, 4, 8, 12, 16, 23 at Δ ≈ 2, D_x = 1000 | `hermpart.json`, `bond23_results.json`, `extend_T_results.json`, C33 `certificate_n*.json` |
| 12 | Distinct centres from the loop link | numerical | n = 2, 3 | `logical_n*.json` |
| 13 | 2n Fermi points, δ₂, δ₄ > 0 | numerical | n ≤ 3 | `logical_n*.json` |
| 14 | Non-Gaussian part ∝ g² | numerical | n = 2, 3, g = 0.05, 0.01 | routeA_spatial §5 |
| 15 | Hierarchy closes at fixed n | heuristic | explicit at n = 2; n = 3 needs g ≈ 10⁻³ (not run) | §3 here |
| 16 | Bands → c = n | heuristic (standard bosonization) | no measured c | ir_control, C20 |
| 17 | Leakage only renormalizes the logical couplings | assumed | none | ir_control H6 |
| 18 | Steps 12–13 for every n, in particular n = 23 | open | crossing for every n proved for C32's variant; centres and margins only n ≤ 3 | C32 §3, C33 |
| 19 | c > d − 1 at n = 23 | heuristic, frozen cells only | rests on 16–18 | this folder |

## 5. Frozen cells, and what the mobile version would add

[proved] Marker positions are exactly conserved. On a ring of length L divisible by b = n + 1, the ground space contains the b translates of the marker pattern, related by the one-site translation, and each carries the same logical chain. The low-energy spectrum is b copies of the c = n spectrum and the ground state is b-fold degenerate. On rings with b ∤ L there are defect cells of the wrong length, costing about γ_b each, with degeneracy growing with L (C7, C22, STATUS step 6). So the chain is c = n on top of spontaneously broken translation. If the conjecture's "unitary CFT(s)" is read to allow identical sectors related by a broken translation, on compatible rings, the frozen chain counts. Under the stricter reading (one CFT, unique ground state on every ring), it does not.

[heuristic, with toy evidence] The mobile version lets markers move so that cells of two adjacent lengths b and b + 1 coexist as a free c = 1 length fluid, with the payload riding on it as in spin–charge separation (Ogata–Shiba gluing). Then every ring has a unique ground state and c = n + 1. C13 and the blind C20 find no unavoidable mixed instability, given exact marker conservation, exact payload parity, Fermi points jointly non-resonant with the length density (margins of order n⁻⁴), positive stiffness and controlled local elimination. The d = 6 toy of `contrib/mobile_cells/` has a unique ground state on every ring L = 8 to 18 (primes included) and c = 2 from three estimators, where the frozen version is degenerate.

What the 23-letter cell lacks for this, and how Astra C32 (finished 2026-10-28, `codex/work_C32/RESULT.md`, not reviewed) supplies it:

1. Equal payloads. Here n = b − 1, so lengths b and b + 1 carry different numbers of flavours, and C20's equal-payload pair does not exist. C32 adds a second (left) marker with five letters A_∅, A_v, A_0, A_+, A_−, so a core is n data sites between two markers, and wraps it in C8's padding. The two physical lengths n + 3 and n + 4 both carry n modes [proved].
2. A frustration-free core. C8's wrapper theorem needs a positive initialized core, and the vacuum walker here has a sine ground state. C32 replaces every hop by its positive two-state edge matrix, so the vacuum is flat and the code exact, with gap at least min(K, m − |g|, 2[1 − cos(π/(n+2))]), m the excited on-site mass [proved].
3. Second-order equalization. C25 showed that endpoint compensation alone does not equalize the two lengths' second-order bonds at finite scales. C32 extends C25's bound to the folded bus, ‖E_(p,q)‖ ≤ C_n (w/U)² + O(s/δ), so the bonds, scalar parts included, agree to any prescribed accuracy at finite nonzero parameters chosen in order. The leading length Hamiltonian is a free hard-core chain, the mixed uλ term vanishes, and there is no second-order three-cell term [proved].
4. Cost. The alphabet is 4(3 · 2 · 28 + 3) · 4 = 2736, independent of n, and the extra marker costs two more powers of D_x in the logical exchange. The mobile chain would beat c ≤ d − 1 only at c = n + 1 > 2735, that is n ≥ 2735. Any cheaper mobility mechanism would change this.

C32 leaves open the same items as the frozen chain (distinct centres and non-resonant margins for every n, local elimination, bosonization stability). Its ED is of the reduced fast-wrapper model at n = 1, 2 with all legal core states, not of the 2736-letter chain.

## 6. What is missing, ranked

1. Distinct centres and non-resonant margins for every n. Band crossing at leading order is now largely in hand: C32 proves an explicit all-n schedule for its variant, and C33 certifies positivity of the leading Hermitian part at n = 23 for this cell (one parameter point) and reports an all-n D(n) and g(n) in progress. No one yet has an argument that the loop link separates all n centres and keeps δ₂, δ₄ and C20's mixed margins positive at every n, in particular n = 23, and today this rests on n ≤ 3. Positivity does not imply it, since a definite diagonal T with a common field puts every opposite-branch pair at the same momentum sum (C32 §3). Double precision cannot reach n = 23 (λ_min/λ_max is 10⁻⁶¹ at the C33 point), so the argument must be analytic or use exact arithmetic of the kind C33 uses. The finite-g transfer for this cell, with explicit g(n), also needs to be finished and checked.
2. A measured central charge. No estimator has been run on this chain. The cheapest decisive check is the n = 2 logical chain (local dimension 4, with the complete 16 × 16 two-cell bond H2 from `ycell.second_order`, so the non-Gaussian remainder is kept) at g = 0.01, with two calibrated estimators (entanglement and Casimir with measured velocities), as `contrib/routeA_incommensurate/` does for the 64-dimensional cell. Then n = 3 at g ≈ 10⁻³.
3. The infrared transfer at fixed coupling. Show at n = 2, 3 with explicit numbers that the λ⁴, g⁴ and loop-induced terms keep the Fermi points non-resonant (shifts below δ₄/16 in C20's sense), and state the local elimination of the gapped cell sectors as an explicit physics input (ir_control H6). At the physics standard the latter is accepted as standard, the former is a computation.
4. The selector and illegal-sector bounds for this alphabet. Re-derive C22's proofs with the express track, the switch and t_d, including the excited gap at every length (Δ ≥ |g|) and λ²α_n⁴ ≪ γ_b. Bookkeeping, but not done.
5. Scope or mobility. Either argue that the conjecture's "CFT(s)" admits spontaneously broken translation on compatible rings, or accept the mobile version (C32, d = 2736), which then needs item 1 and the C13/C20 margins including the length density, and violates the bound only at n ≥ 2735.
6. Independent review, human or otherwise, of the bond in Section 1, of steps 9–13, and of the C32 proofs and C33 certificates.

## 7. Other fixed-alphabet cells

These are alternatives, not parts of the candidate above.

- `contrib/routeA_compile/` central-switch cell. Spinless fermions with exact fermion number, frustration-free Bell-dimer ports, frozen alphabet 48 without a selector. T is Hermitian and definite, giving c = n by Gaussian count at n = 1 to 4, with displaced seas shown at n = 2. Astra C30 (finished) adds a nine-letter symmetric selector, giving d = 120 frozen and 8112 mobile with compensators. Astra C31 (finished) finds the centred-sea pair transfer irrelevant at one loop at n = 2, 3, 4 and leaves general n open. A violation here would need n ≥ 120.
- C28 (17 letters). Frustration-free charged bus with exact U(1), a good doped n = 2 chain with irrelevant pair transfer at one loop, and a proof that diagonal readouts give at most four independent fields, so its field route fails for n ≥ 5.
- C29 (26 letters). Two buses give normal hopping with independent phases at n = 2. The same rules are non-normal at n = 3, with a rational certificate.

## Reproducing

numpy and scipy only, run from this folder. `bond23.py` and `extend_T.py` import `../routeA_spatial/ycell.py` (and `logical_chain.py`) read-only for the comparisons.

- `python3 -B bond23.py` (about a minute) builds the bond from §1.2, runs the Hermiticity, charge, grammar, cell and full-space checks, and computes T at n = 1 to 4 by three methods. Output `bond23_results.json`, `bond23_stdout.txt`.
- `python3 -B extend_T.py` (a few minutes) runs the Richardson test, T⁽²⁾ at n = 2 to 6 for D_x = 4, 10, 20, 40 and two express phases, the reality of the projected T at t_d = 0, and regenerates `logical_chain.analyse` at n = 1, 2. Output `extend_T_results.json`, `extend_T_stdout.txt`.
- `python3 -B gap_delta.py` (about a minute) gives the cell gaps for n = 1 to 6 and the Δ dependence of R_s. Output `gap_delta_results.json`, `gap_delta_stdout.txt`.

## Relation to the rest of the repository

This is the route-A item of STATUS "Where help is most useful" (realize route A's excited-spectrum construction spatially), written as one checkable candidate. It does not change "Where things stand". It suggests three corrections to `contrib/routeA_spatial/README.md`. The D_x = 20 row of its Hermitian-part table is at the double-precision floor. The express phase must scale as π/(2n), because T ∝ sin(nφ^x_0). The infrared conditions at n = 3 need g ≈ 10⁻³, not 0.01. It also records that the cell gap at Δ = 0 falls like n⁻³ and is lowered by the conversion, so Δ ≥ |g| belongs in the schedule. The C33 statements in §2.4, §3, §4 and §6 are from a run still in progress at the time of writing and should be revisited when its RESULT.md appears.
