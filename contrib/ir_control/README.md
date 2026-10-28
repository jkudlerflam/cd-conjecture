# Infrared control of the encoded multi-flavour fermion target

## Claim

Setting. After the frozen-marker reduction (rings with b | L, one marker sector), the candidate of `contrib/programmable_cell/` is a ring of M identical cells. Each cell has a finite space K, an exactly degenerate code G : F_n → K (F_n the Fock space of n complex fermions) with gap Δ above it, and neighbouring cells couple through λV. Schrieffer–Wolff (SW) gives the logical chain W(H_* + δH) plus leakage into excited cells. Here H_* is the incommensurate n-flavour target of C8 §6 (`codex/work_C8/RESULT.md`), and δH is quasi-local with relative size ε.

The statement needed, IR_n: at fixed small λ, the low-energy spectrum on rings of M cells is that of a multi-velocity Gaussian theory of n complex fermions, so c = n.

1. IR_n is covered by no theorem. It is standard but not rigorous, verdict (b). I found no genuine obstruction for the incommensurate target, provided h is even under an on-site Z₂ that acts on the code as logical fermion parity.
2. The fermionic core falls inside the multi-channel Benfatto–Mastropietro framework (Mastropietro–Porta 2022, Porta–Scola–Singh 2025). It differs from their hypotheses in two ways that look routine, no U(1) and a general even interaction. Those theorems control correlations, not c.
3. The SO(2n)-symmetric first-order point and SU(N)₁ sit on a separatrix of their marginal flow. Incommensurate detuning removes every non-Gaussian marginal channel. C8 §11's mass-channel repair is then unnecessary.
4. C8's n = 40 example has Fermi-point margin δ₂ = 4.5 × 10⁻⁵. A Sidon design attains the optimal orders δ₂ ≥ π/(6n) and δ₄ ≥ π/(576n²).
5. Local SW holds to any order k with an exponentially local remainder of order C^k k!(λ/Δ)^{k+1}Δ per cell. That controls energy densities but not level spacings of order 1/M. Exact integration of excited cells is proved only for a special cell.

## Status

Heuristic, verdict (b). §7 lists the three steps to a theorem. The last, a spectral or thermodynamic c for a non-integrable Luttinger liquid, has not been done for any model.

## How it was established and checked

- Author: Claude (Anthropic, Opus 5.5), sub-session of the coordinating session, 2026-10-27.
- AI assistance: gpt-6-astra (high), closed book via `gpt_solve.sh`. Spec A (operator classification, margins, symmetric-point flow) produced `astra_A/`. Spec B (local SW, leakage, imaginary time) produced `astra_B/`.
- Checked by me. I reran `astra_A/margins.py`. `checks_ir.py` recomputes δ₂ for C8's example and both margins of the Sidon scaffold, and integrates the full one-loop so(4) flow. I checked by hand the 2 × 2 annihilation criterion (A §A) and the counterexamples of B §A–§D. I did not re-derive the constants of B's SW theorem. I read the literature myself.

## Argument

### 1. What is rigorously known

- Spinless lattice fermions: Benfatto–Gallavotti–Procacci–Scoppola (CMP 160, 1994) and Benfatto–Mastropietro (J. Stat. Phys. 115, 2004; CMP 258, 2005) [recalled]. Schwinger functions have anomalous exponents analytic in the coupling, and the beta function vanishes. Benfatto–Mastropietro prove D/κ = v² for any weakly interacting chain, integrable or not [read, 1012.0403 abstract].
- Spinful: Benfatto–Falco–Mastropietro, CMP 330 (2014) [read, 1303.3681 §1.1–1.2]. Extended Hubbard, SU(2) and U(1), |v(x)| ≤ Ce^{−κ|x|}, positivity λv̂(2 arccos μ̄) ≥ 0, away from half filling. Response functions (L → ∞, then β → ∞) are Borel summable with logarithmic corrections. This is the only controlled marginally irrelevant flow I found.
- Multi-channel. Mastropietro–Porta, 2108.13708, CMP 395 (2022) [read, Assumption 2.1, Thm 4.1, Rem. 4.2, Eqs. (6.44)–(6.50), Thm 8.5]. Any number of modes with arbitrary velocities, counter-propagating allowed, a number-conserving finite-range two-body interaction, and k_{ω1} − k_{ω2} = k_{ω3} − k_{ω4} only in pairs, with margin γ. Thm 8.5 proves the beta function vanishes for the multi-channel Luttinger reference model. Porta–Scola–Singh, 2509.08665 [read, §1, Assumption 2.1], treat a one-dimensional chain with M internal states and N_f transversal Fermi points.
- Central charge: only the non-solvable 2D Ising model, c = 1/2 from finite-size free-energy corrections with a ferromagnetic inequality (Giuliani–Mastropietro, CMP 324, 2013) [read, 1211.3030 abstract]. Nothing exists for Luttinger liquids.
- Gapless stability outside constructive RG: none found. Fröhlich–Pizzo (1812.02457) and Del Vecchio–Fröhlich–Pizzo–Rossi (2007.07667) need a gap [read, abstracts].
- SW. Bravyi–DiVincenzo–Loss [read, 1105.0675 §4.1–4.5] prove the linked-cluster theorem (Thm 2). With the local SW of Datta–Fernández–Fröhlich–Rey-Bellet [recalled] they get a ground-energy error O(Nε^{n+1}) (Thms 1, 3). That error is extensive, and §4.1 notes that degenerate local ground spaces leave no global gap. Abanin–De Roeck–Ho–Huveneers (1509.05386) and Yin–Lucas (2209.11242) make the per-site remainder (stretched-)exponentially small [read, abstracts].

### 2. The statement and its hypotheses

- H0, selector [heuristic, elementary]. Conserved markers decouple sectors exactly. Defect sectors must cost δ_sel ≫ W uniformly in M. C8's mobile selector would instead add a c = 1 length fluid coupled to the payload, which I have not analyzed.
- H1, exact grading [design]. An on-site Z₂ commutes with every term of h and restricts to P_log on the code. Jordan–Wigner then makes the sector a local fermionic chain. C8 §11 builds graded joint codes, and C8 Thm 3 forbids only large simple symmetries.
- H2, local SW [proved modulo the Cauchy coefficient bound, Astra B §A]. For λ/Δ ≤ [32C(k+1)]⁻¹ there is a generator of range ≤ k + 1 with remainder ‖R_k‖_μ ≤ C_*^k k! λ^{k+1}/Δ^k in an exponentially weighted norm. This corrects my spec. Block diagonality holds for the all-code projection P only, not for the number of excited cells, since degenerate excited configurations mix at first order.
- H3, target [design]. 4n simple zeros of the Bogoliubov symbol, minimum velocity v_min, and margins δ₂ = min d(p_i ± p_j, 2p_i) and δ₄ = min d(Σ_{i∈I} ±p_i) over four distinct indices. These are the momenta of every non-neutral parity-even operator of dimension ≤ 2 (A §B–§C) [derived, bosonization input].
- H4, smallness. A quadratic perturbation of C¹-norm below 10⁻³v_min δ₂ keeps exactly 4n simple zeros, and this is sharp up to the constant, since opposite chiralities annihilate once |h| > vδ/2 (A §A) [proved]. Interacting: ε ≪ v_min δ₂ and ε log(1/δ₄) ≪ 1 [heuristic]. Dressed Fermi points must also be non-resonant. That fails on a parameter set of measure at most 8·C(2n,4)(1 − b_*)^{−2n} r/π for margin r (A §D) [proved, given Lipschitz dressing].
- H5, fermionic IR theorem [open, routine-looking], §3.
- H6, leakage [open], §5.
- H7, correlations to spectrum [open for every non-integrable gapless chain]. Under the c–d premise that the spectrum is a collection of unitary CFTs, H5 supplies 2n independent chiral U(1) currents, so c ≥ n by Sugawara and unitarity [standard CFT; identifying lattice currents with CFT currents is heuristic]. A counterexample needs only this lower bound.

### 3. Is the target within the Benfatto–Mastropietro framework?

Structurally, yes. In Mastropietro–Porta, separated Fermi points make quadratic kernels diagonal, and non-resonance makes quartic kernels density–density (their Eqs. (6.44), (6.47)). The flow reduces to velocities, wave-function renormalizations, per-mode counterterms and density couplings (Eq. (6.49)), which Thm 8.5 controls. Our target differs in three ways.

1. No U(1). Pairing and number-changing quartics carry sum momenta, which δ₂ and δ₄ cover. Below the mismatch scale they vanish by support, as differences do in their Assumption 2.1(d), and Nambu doubling handles the Grassmann representation. The shifts θ_a ≠ 0 of C8 §6 make p_{a,+} + p_{a,−} = 2θ_a non-resonant. With inversion-symmetric Fermi points pairing would gap the chain (B §E).
2. A general even many-body interaction. The construction needs only exponentially decaying kernels. The two-body form serves the lattice Ward identities for the conductance.
3. The output is correlation asymptotics, not c (H7).

So the fermionic core is an unwritten but low-risk extension.

### 4. Targets that do not work as stated, parity, and margins

At the symmetric point every Majorana flavour sits at k = 0 and π, and every current coupling g_{AB}J_R^A J_L^B is allowed. Under ġ_{EF} = ½f^{ACE}f^{BDF}g_{AB}g_{CD}, the isotropic irrelevant coupling lies on the separatrix for m ≥ 4, and an arbitrarily small SO(2) × SO(m − 2) anisotropy runs away (A §G) [proved at one loop; checked: for so(4) at g = −0.05, a +0.01 anisotropy blows up at ℓ ≈ 61 while −0.01 flows to the Gaussian line]. The same Kosterlitz–Thouless structure makes SU(N)₁ need an easy-plane sign on its O(λ³/Δ²) anisotropies. It also has no weakly interacting fermion reference, even at N = 2, where the Heisenberg point lies far from the XX point that constructive RG controls. After detuning, every non-neutral quartic component oscillates, so C8 §11's attractive Gross–Neveu sign is harmless once ε log(1/δ₄) ≪ 1. Equal k_a = k_b must be avoided, since 2(k_a − k_b) is a four-index momentum.

Without exact Z₂, twist fields have dimension n/4. They are irrelevant for n ≥ 9 (A §F) [derived, bosonization input], but they are non-local in the fermions, so no constructive method reaches them. H1 is essential to any proof.

C8's n = 40 example crowds its Fermi points near ±π/2. It has δ₂ = 4.515 × 10⁻⁵ (from p_{31,−} − p_{34,−}) [checked] and δ₄ = 2.1 × 10⁻⁸. For any design, packing gives δ₂ ≤ π/(2n) and δ₄ ≤ 4π/(n(n − 1)). Astra's Sidon points, A_i = 2qi + (i² mod q) with q prime, attain δ₂ ≥ π/(6n) and δ₄ ≥ π/(576n²) [proved]. At n = 40 they give δ₂ = 0.02351, δ₄ = 7.915 × 10⁻⁵ and v_min = 1.41 [checked]. The c–d route has n ≈ q₊q₋, which is huge, so ε must be ≲ v/n.

### 5. Leakage

Per-cell remainders fix the ground-energy density and the free energy per cell uniformly in M and β (B §B) [proved]. They cannot fix the infrared spectrum, since εΣb†b shifts levels by more than 2πv/M. Exact quasi-local removal of leakage would make the dressed code space invariant, while a cell excitation can decay into about Δ/W logical quanta. So exact removal is not expected, though not disproved (B §C). Imaginary time suffices for the infrared, and there excursions decay as e^{−Δs}.

Astra B §D–§E [checked by hand, not re-derived] shows three things. A graded factorization K = K_aux ⊗̂ F_n exists if and only if D₊ = D₋ = Ld/2, after padding. A naive Gaussian-reference expansion fails when H_c acts on the logical factor of excited cells, already at order e^{−βΔ}. Integration is exact when one free auxiliary fermion per cell (H_c = Δa†a) couples linearly, V = h₁ + g(a†B + B†a). That gives a retarded even interaction of weighted kernel norm ≤ 2g²b_n²Δ(λ/Δ)^{3/2}, with b_n = 8^{2n}. For general cells an excursion-polymer expansion is the missing new mathematics.

### 6. Numerics

[RESULTS PENDING]

### 7. Verdict

(b), standard but not rigorous. Only Fermi-point shifts survive as relevant operators. The encoding forces no resonance, because the first-order bond Σt_{ab}R_aΓ_b and on-site bilinears reach any quadratic target. A theorem needs three steps.

1. Extend Mastropietro–Porta and Porta–Scola–Singh to U(1)-free, general even interactions with sum-type non-resonance (H5). This looks routine.
2. Prove a leakage lemma (H6), either by building cells of the B §E form or by an excursion-polymer expansion. This is new but standard in spirit.
3. Bridge to c (H7). Either use the current-algebra bound c ≥ n under the c–d premise, or prove a spectral or thermodynamic c for a Luttinger liquid, which would be a first.

Without step 2, steps 1 and 3 plus a dilute-gas bound on thermally excited cells would still give a rigorous prethermal statement. The free energy matches the logical chain's to O(ε_k), so c = n would be visible above T ≈ √(vε_k), which is exponentially small at optimal k.

## Reproducing

- `spec_A_operator_classification.md` and `spec_B_leakage_expansion.md` are the Astra specs. `astra_A/` holds RESULT.md, `margins.py` (seconds) and `numerical_results.json`. `astra_B/RESULT.md` holds the other output. The `.out` and `.err` files are the raw runs.
- `checks_ir.py` holds my checks (numpy and scipy, seconds). Its output is `checks_ir_stdout.txt`.
- `toy_ir.py` imports `contrib/counterexample_search/idmrg.py`. `tasks.txt` and `slurm_toy.sh` are the Typhon array, and the results go in `toy_ir_results.jsonl`.

## Relation to the rest of the repository

This addresses open item (ii) of `contrib/programmable_cell/` §5 and REPORT §9.6. "Where things stand" could say the following. Infrared control of the incommensurate encoded target is standard physics, not a theorem. Its fermionic core is a routine extension of the multi-channel constructive RG. Leakage and spectral c are the open mathematical steps, and neither looks like an obstruction. Exact on-site Z₂ and designed Fermi-point margins are hard requirements, and the symmetric and SU(N)₁ targets should be dropped.
