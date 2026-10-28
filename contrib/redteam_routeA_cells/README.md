# Red team: the 23-letter folded-bus chain and the switch cell

## Claim

Under review is `contrib/routeA_spatial/`. It claims a one-site translation-invariant, nearest-neighbour, Hermitian bond on C^23 with frozen cells of length b = n + 1, an exactly degenerate 2^n code per cell, exact U(1) and parity, zero first order, no three-cell term at second order, and a second-order logical chain of n free bands with distinct centres. With the weak-coupling input of `contrib/ir_control/`, it concludes c = n at d0 = 23 on compatible rings for every n, hence c > d0 − 1 once n ≥ 23. Secondarily, `contrib/routeA_compile/` claims a number-conserving, frustration-free switch cell with c = n at n ≤ 4.

My verdict. The finite-dimensional claims survive an independent reimplementation, and I found nothing that breaks the logical chain at fixed n. The scale hierarchy closes from the top down at fixed nonzero couplings. Three things stop this from being a counterexample today. The chain is frozen, so it fails Latorre and Sierra's hypothesis on most ring lengths (§5). The stated working point is not uniform in n (at Δ = 0 the code stops being the cell ground space at n = 8), and the n = 3 point is outside its own weak-coupling criterion. Both are fixable. And n ≥ 23 is an extrapolation of a structure I could verify only through n = 10.

## Status

Refereed. Numerical and heuristic. Three things would upgrade it: a mobile version with a unique ground state on every ring, a proof of positivity at finite express detuning D_x, and one explicit parameter set written down for general n.

## How it was established and checked

- Author: Claude (Opus 5.5) red-team subagent for Jonah Kudler-Flam, 2026-10-28.
- AI assistance: none beyond this session. No Codex, since every check is a short script whose specification would be as long as the code.
- Checked by: me. I read `ycell.py` for conventions, but built the bond from the README's rules as one explicit 529 × 529 matrix, generated every cell and chain from that matrix and the arrow grammar, and computed the second-order bond two new ways. One is a sparse direct solve in the two-cell excited space, the other a data-free single-token formula derived here (`rt_token.py`). Neither uses the Cauchy eigen-sum of the original.

## Findings

### 1. What survives

- Bond. Hermitian, [h, Q⊗1 + 1⊗Q] = 0, marker number conserved on each site, outer arrows preserved by every move, all 23 letters used. The count is 2 (L) + 2 (R) + 8 (stub token, species v,0,+,−) + 6 (express token) + 5 (markers) = 23. The selector is included, since μ_sel sits on all five marker letters and the vacuum walker is C22's path.
- Literal rings. Two-cell rings built from the split bond (n = 1, 2) have the same spectrum as cells plus port, to the last bit.
- Selector. On the full 23^4 ring at n = 1 (all 279,841 words, illegal ones included) with μ_sel = 0.1, the eight lowest states sit at −1.8 = 2(a_2 − δ), two marker patterns times 2^2 code states, and the next level is the single b = 4 cell at −1.518. At λ = 0.1 the same eight form a band 0.28 below the rest. No published check sets μ_sel ≠ 0.
- Code. 2^n-fold degenerate to 8 × 10⁻¹⁵ at n = 1, 2, 3, gaps 0.704, 0.272, 0.121, each vector the walker ground state times physical data. First order ≤ 3 × 10⁻¹⁵. H2 commutes with fermion number to 4 × 10⁻¹⁷.
- No three-cell term. In the full three-cell space (10,648 states at n = 1, 373,248 at n = 2) P V R V P equals the sum of the two-cell bonds to 6 × 10⁻¹⁶.
- T and R_s. The sparse solve reproduces |T| = 3.56 × 10⁻⁷, 1.0996 × 10⁻⁸, 1.2721 × 10⁻¹⁰ and the eigenvalues of R_s/|T|, (0.01758, 1.2219) at n = 2 and (2.1111 × 10⁻⁴, 0.06416, 2.0163) at n = 3, to at least five digits. The single-token formula agrees up to a relative O(g²) (2 × 10⁻⁵ at g = 10⁻³, n = 3).
- Fourth order. In the full two-cell legal space the data-dependent λ⁴ terms are suppressed like T itself, T4/T2 = 0.20 (n = 2) and 0.26 (n = 3), with one-cell terms 10⁻⁴ and 2 × 10⁻⁵ of |T2| and a real part zero to noise. Particle-hole protection survives, and two-cell ED at n = 2 matches λ²H2 + λ⁴H4 to 6 × 10⁻¹⁵. Fourth-order one-body terms stay below a tenth of the narrowest window for λ ≲ 0.13 (n = 2) and 0.013 (n = 3).

### 2. The Gram limit, checked through n = 10

The single-token formula exposes the structure. At t_d = 0, T = (phase) × sin(nφ^x_0) × M_n(D_x), with M_n real to 10⁻¹⁶. Two channels carry conjugate express phases and cancel except for the sine. As D_x grows, M_n/|M_n| approaches the bare-stub Gram matrix Γ_ij = ∫ f_i f_j dτ with deviation 0.6/D_x, so the limit claimed in the README holds. The smallest eigenvalue of R_s/|T| at large D_x is

| n | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| λ_min(R_s)/\|T\| | 4.0e-2 | 1.8e-3 | 5.1e-5 | 1.6e-6 | 4.4e-8 | 1.4e-9 | 4.1e-11 | 1.3e-12 | 4.1e-14 |
| D_x needed (Δ = 0) | ≤ 2 | ≤ 2 | 4–10 | 10–30 | 10–30 | 10–30 | 20–30 | 20–30 | 30–50 |

That is about a factor 31 per flavour, while D_x grows slowly, far better than the naive D_x ~ |Γ|/λ_min(Γ) but still numerical. Since |T| scales like D_x^{−2(n+1)}, an extrapolation to n = 23 puts the narrowest logical band below 10⁻¹²⁰ in units of t, times λ²g².

### 3. Problems, ranked

| # | Problem | Rank |
|:--|:--|:--|
| 1 | Frozen cells. Not conformal on rings with L ∉ bZ, b-fold translation breaking on L ∈ bZ (§5) | Fatal for the claim that this chain refutes Latorre–Sierra as stated. Not fatal for the program |
| 2 | No mobile version. The payload n = b − 1 grows with the cell, so C20's equal-payload pair of lengths does not exist. The code is the ground state of a frustrated hopping walker (the selector itself), while C8's wrapper theorem needs a positive frustration-free core | Serious |
| 3 | n ≥ 23 rests on finite-D_x positivity, verified only to n = 10, with margins near 30^{−n} | Serious, probably fixable |
| 4 | At the working point (Δ = 0, g = 0.05) the excited-token sector drops below the code at n = 8. The cell gap is −0.0044, −0.0080, −0.0103 at n = 8, 9, 10, so the code is not the ground space and long segments undercut the selector. Δ shifts excited levels rigidly, so Δ ≳ g restores both, and R_s stays positive at Δ = 0.3 through n = 8 (D_x = 100) | Fixable |
| 5 | The n = 3 point is outside its own criterion. The non-Gaussian part over v_min δ₂ is about 33 at g = 0.05 and 1.3 at g = 0.01. Self-consistent Hartree–Fock keeps 2 + 2 + 2 Fermi points, but at g = 0.05 it reverses the centre splitting (±0.0013 becomes ∓0.0027), which is the quantity cutting off pair transfer. At g = 0.003 the shifts are below 10⁻⁴ \|T\|. The README's "all conditions are met at n = 2 and 3" holds at n = 2 only | Fixable (smaller g) |
| 6 | The margin-setting process is inter-band opposite-branch pair transfer, momentum 2(θ_a − θ_b). With exact U(1) and t_d = 0 it is resonant (δ₄ = 2.6 × 10⁻¹⁰ at n = 2) because T is real after its phase. The loop cures this, but the largest t_d keeping R_s > 0 falls with the express suppression (n = 3: 3 × 10⁻⁴, 6 × 10⁻⁵, 3 × 10⁻⁶ at D_x = 4, 10, 30), so t_d must follow D_x and precede g | Fixable |
| 7 | φ^x_0 = 0.4 is n-specific. At t_d = 0, |T| ∝ \|sin(nφ^x_0)\| exactly, so T vanishes at nφ^x_0 = π and is suppressed 17-fold at n = 8. Use φ^x_0 = π/(2n) | Cosmetic |
| 8 | "μ vanishes" holds in half-filling normal order. Referred to the empty state, Re μ/\|T\| = 0.014 (n = 2) and 0.037 (n = 3), a Hartree term, not a symmetry failure. The ED's "error/λ⁴" is almost all scalar (traceless H4 2.8 × 10⁻⁹ against scalar 5.0 × 10⁻³ at n = 2) | Cosmetic |

Nothing else I tried leaked. Umklapp is not near resonance at the computed fillings (0.71, 0.50, 0.23 at n = 3, and same-flavour umklapp is irrelevant near K = 1). Illegal sectors cost of order K, defect cells cost γ_b ~ t b⁻⁴, and the cell gap is O(Δ − g). All three are fixed by n and Δ, while every logical scale carries λ², so the tiny band does not compete with them.

### 4. Satisfiability of the hierarchy

At fixed n the ordering is Δ (O(0.1)), then D_x ≥ D*(n, Δ) and φ^x_0 = π/(2n), then t_d below t_d*(n, D_x), then g below Δ and small enough for the non-Gaussian criterion, then ε inside the window, then λ, with λ² e0 below the selector margin. Each choice is an upper bound in terms of earlier ones. The only lower bound, D_x ≥ D*, sits on a large scale. I found no circular constraint, so at fixed n a fixed nonzero coupling set exists, granting the finite-D_x positivity of §2.

### 5. The frozen-cell issue

Latorre and Sierra assume periodic boundary conditions, one-site translation invariance and "a collection of unitary CFTs" describing the low-energy spectrum (arXiv:2403.17242, Eq. 3). They do not discuss spontaneous translation breaking. They do say that blocking k sites raises the local dimension to d^k.

The frozen chain has an exact b-fold degenerate ground state on L ∈ bZ, the translates of one marker pattern, each carrying the n-band theory. Read as CFT ⊗ (Z_b broken), c_T = n. On every other L the selector forces defect cells, the ground state is at least L-fold degenerate with immobile defects, and the spectrum is not conformal (C8 §12). So the hypothesis holds only along the subsequence L ∈ bZ, or for the infinite chain's symmetry-broken sectors. Under that reading the chain would be a counterexample at the physics standard. Under the all-ring reading it lies outside the conjecture rather than against it. Physically, period-b breaking turns the chain into one with a b-site unit cell (here 23^b ≥ 2^n states), which is the enlargement Latorre and Sierra's own blocking remark accounts for. A referee would most likely ask for the mobile version, with a unique ground state on every ring and c = n + 1. For this cell that is problem 2, not a formality.

### 6. The switch cell (secondary)

The same-schedule lemma is correct, since h(k) is congruent to N + 2cos k and Sylvester inertia moves one band through zero at each crossing. Exact fermion number and a frustration-free code make this the better candidate for a mobile wrapper. Three problems remain.

- The arm automorphism that makes T Hermitian and definite also pins every sea centre to 0 or π, so 2(θ_a − θ_b) ≡ 0 mod 2π and every inter-flavour pair transfer is resonant. Unequal arms without a phase give centres 0 and −π (`switch_break_m2.json`), still resonant. Displacement needs asymmetry and a phase, shown only at n = 2, the same tension as item 6 of §3. Fixable.
- No selector. Letter types are frozen and every pattern (P_L L^a M R^c P_R)* is zero-penalty, so a = c and the length are not selected at d = 48. Serious until built.
- Frozen, as in §5.

### 7. Verdict

The 23-letter chain is a sound finite construction. Every decisive finite claim reproduces, generated terms through fourth order are harmless, and the fixed-n hierarchy is satisfiable. It is not yet a counterexample to c ≤ d − 1 as Latorre and Sierra state it, because it is frozen, and its large-n regime rests on numerically supported but unproven positivity at margins near 30^{−n}. The highest-value next steps are (i) a mobile, equal-payload version of the folded bus, which will cost letters, (ii) a bound on the finite-D_x correction to the Gram limit, and (iii) one parameter schedule written for general n with Δ > g, φ^x_0 = π/(2n) and g chosen against the pair-transfer margin.

## Reproducing

numpy and scipy only, run from this folder with `python3 -B`.

- `rt_bond.py`: the literal bond and its checks (seconds).
- `rt_cell.py`: grammar-generated cells, code, sparse second-order solve. `rt_quick.py` prints §1's code, T and R_s numbers (`rt_quick_stdout.txt`).
- `rt_selector.py ring | full | comp nmax Δ g`: literal rings (`ring_check_stdout.txt`, minutes), the 23^4 ring (`full_ring_L4.json`, 25 min on 8 cores), competing lengths (`competitors_*_stdout.txt`, minutes to n = 10).
- `rt_threecell.py 1 2` (`threecell.json`, about 20 min at n = 2).
- `rt_token.py` (single-token formula), with `scan_gram.py`, `scan_phase.py`, `scan_td.py`, `scan_delta.py` for §2 and items 4, 6, 7 (minutes).
- `rt_fourth.py n td g` (fourth order and ED), `rt_ir.py n td g ...` (Hartree–Fock; the n = 3, g = 10⁻³ row is at the double-precision floor), `rt_margins.py` (margins).
- `slurm_sel.sh`: the Typhon array used for `full_ring_L4.json`.

## Relation to the rest of the repository

This reviews `contrib/routeA_spatial/` and, more briefly, `contrib/routeA_compile/`, in the spirit of `contrib/redteam_construction/`. It does not change "Where things stand". STATUS's description of the 23-letter chain is accurate, but should add that the working point needs Δ > g for n ≥ 8, that the n = 3 numbers need smaller g, and that "frozen" means the chain does not meet the conjecture's hypothesis on most ring lengths.
