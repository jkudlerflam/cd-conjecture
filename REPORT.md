# c–d conjecture REPORT (2026-10-21): can c ≤ d − 1 be proved?

Paper: J. I. Latorre and G. Sierra, "The c-d conjecture", arXiv:2403.17242 (v1 March 2024, v2 Nov 2024), J. Stat. Mech. (2024) 113103.
Inputs: both arXiv versions, read in full. Literature map `lit.md` (Opus agent, every claim tagged [read] or [recalled]).
Codex: C1 is the free-fermion theorem and C2 is the lattice-Virasoro numerics, specs in `specs/` and transcripts in `codex/`, all closed-book on gpt-6-astra high (§6).
Nothing here is refereed. §6 lists what I checked.

## Bottom line
The conjecture is open and probably true as stated. No citing paper proves a case, and no counterexample is known within its scope.
A general proof is blocked, and the block can be located.
No known lattice quantity bounds c from above and grows linearly in d.
Entanglement is capped at O(log d) per site. Kinematic data (‖h‖, the Lieb–Robinson velocity, sum rules) only bound c times a power of the width of the conformal window. Anomalies bound c from below.
The modular c-theorem toolkit cannot prove it, and §3.1 shows that no lattice-scale entropic argument can.
Proofs exist now in restricted classes: quasi-free chains (c ≤ log₂ d, which is stronger), toroidal IR theories whose bosons come from on-site U(1)s, and WZW IR theories obeying a "parton" level bound (§4).

## 1. The precise statement [read, v2 Eqs. (2)–(3)]
The setting is H_L = Σ_{i=1}^{L} h_{i,i+1} on (ℂ^d)^{⊗L} with periodic boundary conditions.
h is one Hermitian d²×d² matrix, the same on every bond, so the chain is translation invariant by one site.
If the low-energy spectrum is "a collection of unitary CFTs" with central charges c_a, possibly with different velocities (as for Hubbard spin and charge), then Σ_a c_a ≤ d − 1.
- Scope. "Nearest-neighbour" means h acts on two sites, and longer range is reduced to it by blocking (a term spanning r + 1 sites gives d^r). "Critical" means a conformal, z = 1, low-energy theory, and the authors stress that conformality is essential.
- Not covered: non-unitary chains (c_eff), open boundary conditions, anyonic chains, and the 2D chiral analog. The authors list all four as open.
- Chirality is not discussed. It need not be, since Kapustin–Spodyneiko (PRL 123, 060601) forbid equilibrium energy currents on such lattices, which gives c_L = c_R [read].
- Support is a table of about a dozen solvable or well-studied chains, of which only Uimin–Lai–Sutherland (ULS, SU(N)_1) saturates. The paper also gives the coarse bound (c/3) log L < L log d, superadditivity of c_max(d) under stacking, and a heuristic in which decoupled levels are dropped.
- v1 had an appendix proposing a ground-state Koo–Saleur ratio as a proof route. v2 removed it, and §3.3 shows it fails.
- Since ULS attains c = d − 1, the conjecture is equivalent to the equality c_max(d) = d − 1.

## 2. Literature status [lit.md]
INSPIRE lists seven citing papers. None proves a case or gives an in-scope counterexample.
Li–Lin–McGreevy (2509.04596) find agreement in small-system numerics. Prembabu–Verresen (2509.20431) show that the qutrit XX chain (d = 3, c = 2) saturates both c–d and the Lieb–Schultz–Mattis (LSM) lower bound. Kumar–Pujari–D'Emidio (2506.16424) report a "non-unitary violation" at non-integer N, which is really a tensor-product effect (§3.6).
Each closest miss drops one hypothesis:
- Range: Rahmani et al. (PRB 1505.03966) have d = 2 and c = 3/2 with a σˣ_jσˣ_{j+2} term, and Σ S_j·S_{j+2} trivially gives c = 2.
- Translation invariance: aperiodic d = 2 chains reach c_eff = ln 4/ln 3, but they are not CFTs.
- Conformality: Motzkin chains.
- Integer d: anyonic chains (§3.6).

LSM gives c ≥ d − 1 in the SU(N)-fundamental and Z_d × Z_d classes, so with c–d it would pin c = d − 1 there.

## 3. Routes

### 3.1 The entropic / modular toolkit (C_± of the note): no-go [checked]
There are two obstructions.

(a) Structure. C_± rests on half-sided modular inclusions from positive-energy null translations, the ant formula and nesting, all of which need Lorentz covariance.
A lattice has none at the lattice scale. Its only continuum object is the IR CFT, a fixed point, so there is no QFT flow whose UV is "the lattice".

(b) Quantitative. Suppose a lattice c-function C(ℓ) tends to c in the conformal regime, is non-increasing, and obeys C ≤ K log d at one lattice spacing. Then c ≤ K log d.
- For the discrete Casini–Huerta function 3ℓ[S(ℓ+1) − S(ℓ)], subadditivity gives K = 3.
- For the lattice analog of C_±, −6a ∂_a I(half-line : interval), the chain rule I(L:A) − I(L′:A) = I(x:A|L′) ≤ 2 log d gives K = 12.
- ULS has c = d − 1 > K ln d once d ≥ 7 (K = 3) or d ≥ 48 (K = 12).

For those chains every such function must increase somewhere, so monotonicity fails on the lattice.
No entropic argument starting at the lattice scale can prove c–d, because it would prove the false bound c = O(log d).
Equivalently, the SU(d) conformal regime can only begin at ℓ* ≳ (d − 1)/(3 log d), up to logarithms.
A proof needs an input that counts something linear in d.

### 3.2 Thermodynamics, Lieb–Robinson, sum rules [checked]
In the conformal window T ≤ T*, s = πcT/(3v). With s ≤ log d this gives cT*/v ≤ (3/π) log d.
Lieb–Robinson, v ≤ v_LR = O(‖h‖), then gives cT* ≤ (3/π) v_LR log d.
The f-sum rule for an on-site current algebra at level k gives k·v ≤ const.
Each bound constrains c times a positive power of the window, never c alone.
ULS makes this concrete: h = P has spectrum {±1} and ‖h‖ = 1 for every d, and v_LR is bounded uniformly, while c = d − 1.
So no bound in terms of ‖h‖, v_LR or the eigenvalues of h can exist. d must enter through multiplicities, and the SU(d) window must close at T*/v ≤ 3 log d/(π(d − 1)).

### 3.3 Lattice Virasoro, the v1 Koo–Saleur ratio: fails [C2 + my analysis]
v1 proposed c₂(N) = 8m₀³/m₁², with m₀ = ⟨h₂h₋₂⟩ and m₁ = ⟨h₂(H − E₀)h₋₂⟩.
High-energy states contribute q²·(spectral weight of the total energy current J_tot)/Δ² ∝ 1/N to m₀, while the CFT signal is ∝ 1/N².
That 1/N weight is absent only when J_tot is conserved, the lattice shadow of "energy current = momentum density", that is, of Lorentz invariance. So generically c₂ ~ K/N → 0.

Even conservation is not enough. For critical Ising with the centred density, J_tot is conserved, but a q⁴ form factor still feeds m₁ at order N⁻³, the same order as the CFT term. The result is c₂ → 9/98, not 1/2 [Astra, analytic; arithmetic rechecked].
With the literal density −X_jX_{j+1} − Z_j the same chain gives c₂ ≈ 9π⁴/(32N), so c₂ depends on how H is split into local terms.

Exact diagonalization up to N = 22 agrees. For the non-integrable O'Brien–Fendley chain (λ = 0.2), c₂ falls from 0.83 to 0.26, still preasymptotic (current moments give K = 0.005 to 0.33 across the O'Brien–Fendley and J₁–J₂ chains).
The Milsted–Vidal estimator, which projects onto the stress-tensor state, extrapolates to 0.49996 (Ising), 0.498 (O'Brien–Fendley) and 1.0005 (J₁–J₂). Heisenberg comes close even unprojected (c₂ = 1.043 at N = 22), which is why v1's figure looked fine.

The lesson: extracting c needs a projection onto the low-energy subspace, and lattice locality says nothing about that projected weight. This route computes c but does not bound it.

### 3.4 Representation theory: a reformulation [checked, script]
Suppose the IR is a G_k WZW model with G on-site, acting on ℂ^d through a representation R, so c = k dim G/(k + h^∨), increasing in k.
Gutzwiller-projected partons (one fermion per site in d flavours) give SU(d)_1 ⊃ G_{ind R}, where ind is the embedding index.
GKO coset unitarity gives c(G_{ind R}) ≤ d − 1, so the level bound k ≤ ind(R) implies c ≤ d − 1.
The table lists k_max, the largest k with c(G_k) ≤ d − 1:

| G | R | d | ind R | k_max |
|---|---|---|---|---|
| SU(N) | fundamental | N | 1 | 1 |
| SO(N) | vector | N | 2 | 2 |
| Sp(2n) | fundamental | 2n | 1 | 1 |
| SU(2) | spin 1 | 3 | 4 | 4 |
| SU(2) | spin ≥ 3/2 | 2s+1 | – | none needed |
| SU(N) | antisymmetric² | N(N−1)/2 | N−2 | N−2 |
| E₆ | 27 | 27 | 6 | 6 |
| G₂ | 7 | 7 | 2 | 3 |

For conformally embedded R (the rows with k_max = ind), the conjecture restricted to WZW IR theories is equivalent to the level bound.
For example, a d = 4 Sp(4) chain must not reach Sp(4)_2 (c = 4), and an SO(N)-vector chain must not reach SO(N)_3.
The level bound is a parton heuristic, not a theorem. Lattice data fix k·v (Drude weight) and k/v (susceptibility), whose product ∝ k² is computable but has no lattice bound.

Toroidal class. Suppose the IR is c compact bosons whose momentum U(1)^c is the image of an on-site torus acting faithfully on ℂ^d modulo phases. That torus has rank ≤ d − 1, so c ≤ d − 1.
The saturating ULS chain is of this type, through the Cartan U(1)^{d−1}. So is t–J (d = 3, c = 2 [recalled]), through charge and S^z.
In both, d − 1 counts the independent on-site densities.

### 3.5 Anomalies and frustration-free chains [lit.md]
LSM results bound c from below. Yao–Hsieh–Oshikawa plus Sugawara give c ≥ N − 1 for SU(N) fundamentals, and there is a Z_d × Z_d analog. None gives an upper bound.
Gapless frustration-free chains have gap O(1/L²) [Gosset–Mozgunov, recalled], so they are never conformal.

### 3.6 Integrality of d is essential [own arithmetic, checked]
SU(2)_k spin-1/2 anyon chains are unitary, translation invariant and nearest-neighbour in the anyon sense. They have c = 1 − 6/((k+1)(k+2)), which exceeds d_q − 1 = 2cos(π/(k+2)) − 1 for every k ≥ 2. For the golden chain that is 0.7 > 0.618.
The one-Majorana-per-site chain has d = √2 and c = 1/2 > √2 − 1.
Kumar et al.'s non-integer-N "violation" is the same Temperley–Lieb algebra at N = d_q, g = (k+1)/(k+2), where it is unitary. So it comes from the missing tensor product, not from non-unitarity (lit.md reached the same conclusion).
A proof must use the integer tensor-product structure. Arguments that see only Hilbert-space growth (entropy, Cardy counting, loop algebras) cannot work.

### 3.7 Lin–McGreevy's entropy function (arXiv:2303.05444, Jonah's suggestion) [read in full; bounds checked by script]
LM work with three consecutive intervals A, B, C with cross ratio η. They show that a 1+1D CFT ground state satisfies the vector fixed-point equation K_Δ|ψ⟩ = (c/3)h(η)|ψ⟩, with K_Δ = K_AB + K_BC − η(K_A + K_C) − (1 − η)(K_B + K_ABC) and h the binary entropy. Equivalently, the ground state is a stationary point of S_Δ = η[S(B|A) + S(B|C)] + (1 − η) I(A:C|B), which equals (c/3)h(η) there.
They conjecture the converse, that this characterizes CFT ground states. From one- and two-site entanglement Hamiltonians they reconstruct H, P and the other global conformal generators, with intrinsic normalization (no velocity factor).
Li–Lin–McGreevy (2509.04596, lit.md §1) is the follow-up that searches for such fixed points on 4 qubits and 4 qutrits and finds agreement with c–d.

What it buys for c–d.
- The velocity ambiguity of §3.2 is gone, since H_rec is normalized by the modular flow.
- c is read off at a single scale as the value of an entropy function, which turns c–d into a statement about vector fixed points on chains of d-dimensional sites.
- The (1 − η) piece, the CMI of adjacent intervals, is the project's κ object and is cutoff-free in AQFT. The η piece is a sum of conditional entropies and has no AQFT meaning.
- For an exact fixed point whose 3-site marginal ρ_ABC has full rank, the vector equation is equivalent to the operator identity K_Δ = (c/3)h(η)·1 on ABC, because (K_Δ − λ)ρ_ABC^{1/2} = 0 [my observation, elementary]. This is a weighted quantum-Markov condition on finite-dimensional marginals, a finite algebraic problem that could plausibly be classified for d = 2, 3.

Why it is not yet a proof. S_Δ is entropic, so the §3.1 ceiling applies.
- With blocks of ℓ sites, S(B|A) ≤ S_B and I(A:C|B) ≤ 2 min(S_A, S_C), so S_Δ ≤ 2ℓ log d. At η = 1/4 this gives c ≤ (6/h(1/4)) ℓ ln d ≈ 10.7 ℓ ln d.
- If the state is an exact fixed point at every scale, the (1,1,k) equations leave a two-parameter solution family, which forces S(n) = (c/3) ln n + S(1) (I checked that this family solves every triple). Subadditivity S(2) ≤ 2S(1) then gives c ≤ 3 log₂ d.
- SU(d)_1 exceeds these bounds once d ≥ 41 and d ≥ 12 respectively. So the ULS ground state becomes a vector fixed point only at block size ℓ₁ ≳ (d − 1)/(10.7 ln d).

Result of C3 (d = 2 exact single-site fixed points; Astra high, closed-book, with my independent root-finder agreeing).
- Classical case, proved: the only full-support solution in a common on-site basis is the tracial state, with c_Δ = 3(½ ln 2)/h(1/4) = 1.8489.
- Quantum case, numerical (about 6,600 starts, and my 300 independent ones). There are eight classes modulo U ⊗ U, some isolated and some one-parameter families, all with τ = I/2. They take six values of c_Δ: 1.1288, 1.1579, 1.4617, 1.4889, 1.5019 and 1.8489 (tracial). So every nontrivial exact solution has c_Δ > d − 1 = 1, and none has 0 < c_Δ ≤ 1. The global supremum is unproved, and the rigorous bound is 6 ln 2/h(1/4) = 7.40.
- One solution is fully explicit: σ = (I + r XY)/4 and ρ₃ = (I + r(XYI + IXY))/8, with r = tanh j₀ and tanh j₀ = tanh(4√2 j₀/3)/√2. I checked it to 10⁻¹⁶ and it gives c_Δ = 1.4889.
- Purity does not remove them. Pure translation-invariant states realize the tracial marginals (a graph state with stabilizers X_j Z_{j−2} Z_{j+2}). The XX class with c_Δ = 1.1288 is realized by a free-fermion Fermi sea with a finite-range dispersion.
- None is scale consistent. The (1,1,2) equation fixes ρ₄ uniquely from ρ₃, and the (2,1,1) and (1,2,1) equations then follow automatically, because h(1/9) = (4/3)h(1/3) − (8/9)h(1/4) (checked). For every nontrivial class the marginal of ρ₄ misses ρ₃ (Frobenius 0.005 to 0.011). The tracial state passes the marginals but has λ₁₁₂ off by 0.070 (my check gives 0.0698). In the classical case this proves no full-support scale-consistent solution exists. In the quantum case none was found.
- For comparison, true critical chains are not exact fixed points at this scale. Their single-site values are c_Δ = 0.5246 for critical Ising (I reproduced this) and 1.164 for Heisenberg, which already exceeds d − 1 at the lattice scale.

So the single-scale LM condition at the lattice scale is too weak. It admits exact "fixed points" with c_Δ > d − 1, and none of them survives one more scale.
LM therefore locates where conformality starts, but cannot provide the linear-in-d input. A proof along these lines would need a bound on c that is linear in d at the onset scale ℓ₁, and nothing entropic gives one.

## 4. Verdict
- General: blocked. The conjecture is a "rank-like" bound, and no known lattice functional is rank-like.
- Provable now:
  - (i) Quasi-free NN chains of any period. c = (number of roots of det h(z) on |z| = 1)/4 ≤ rank A(1)/2 ≤ n/2 = log₂ d, with n Majoranas per site, c_L = c_R, and sharp examples. This is a theorem. Its proof is a Laurent-degree count, independent of the period because only one bond crosses each cell boundary. It is strictly stronger than d − 1 for d ≥ 4, and interactions are needed to reach d − 1 (SU(3)_1 has c = 2 > log₂ 3).
  - (ii) Toroidal IR theories with on-site momentum U(1)s (rank argument).
  - (iii) WZW IR theories with k ≤ ind(R) (coset). For conformally embedded R this class is equivalent to the conjecture.
  - (iv) Trivially, SU(d)-symmetric fundamental NN chains, since the only invariant h is a + bP.
- Modular c-theorem project: no proof route. The toolkit lives in the continuum, and §3.1(b) shows that even a lattice version would give only O(log d).

## 5. Next calculations, ranked
1. A falsification scan at d = 2. NN qubit chains form an 8-parameter family (12 couplings minus SO(3) minus scale). Search the multicritical points, where Ising lines meet c = 1 regions, for c > 1 with iDMRG entanglement scaling and ED spectra. It is cheap and the most decisive test available.
2. Level tests where they would bite. The Sp(4)-symmetric d = 4 chain is a one-angle family (projectors on 1 ⊕ 5 ⊕ 10): look for Sp(4)_2. Do the same for SO(N)-vector chains reaching SO(N)_3.
3. Settle the parton level bound k ≤ ind(R). Either prove it, for instance through Drude weight × susceptibility ∝ k², or find a G-symmetric NN chain that violates it.
4. Look for an interacting analog of channel counting: a lattice functional linear in d (for example the number of independent local conserved densities) that equals d − 1 on ULS and dominates c on the lit.md examples.
5. The Lin–McGreevy problem (§3.7). d = 2 single-site is done (C3): exact solutions exist, all have c_Δ > 1, and none is scale consistent. What remains is a quantum no-go for the two-scale system, (1,1,1) plus (1,1,2), at d = 2. The interesting version is the approximate one: maximize the c_Δ that is stable across scales over approximate fixed points, extending Li–Lin–McGreevy's 4-site search. The onset scale is still the crux.
6. Make the shrinking window explicit for ULS: T*/v and ℓ* versus d from Bethe-ansatz thermodynamics, and the non-monotone lattice C(ℓ) at d = 7 or 8 by iDMRG.

## 6. Provenance and what was checked
- Checked by me:
  - the paper's statement and scope (both versions);
  - §3.1 (bounds and thresholds by script) and §3.2 (standard formulas);
  - the §3.4 table (scratchpad script) and the coset argument;
  - the §3.6 arithmetic;
  - §3.7: LM v2 read in full, plus the S_Δ bounds, thresholds and the two-parameter solution family (by script). The equivalence of the vector and operator equations at full rank I checked by hand. C3 (`specs/C3_LM_fixedpoints_d2.md`, `codex/C3_LMfixed_astra.out`, code and data in `codex/work_C3_astra/`, my scripts in `codex/opus_check_C3/`): my own root-finder found the same six c_Δ values, and I checked the four-site failure, the explicit E₁₂ solution, the h(1/9) identity and the Ising value. I did not check Astra's classical proof line by line (it is short and looks right), its pure-state realizations, or the Heisenberg value;
  - the C1 proof, line by line;
  - the C2 1/N estimate (my own derivation agrees);
  - the arithmetic of Astra's 9/98, including its two k-integrals. I did not rederive its pair form factor.
- Taken from lit.md without reopening: everything in §2 and §3.5 marked [read] there. The Gosset–Mozgunov and t–J statements are [recalled].
- Codex: both specs first ran on gpt-5.6-sol, which was the shim default at the time, and that was my mistake. Both were rerun on gpt-6-astra high, closed-book in empty directories.
  - C1 (`codex/C1_freefermion_{sol,astra}.out`): both models give the same theorem, proof structure, sharp examples and non-Hermitian caveat. The literature agent found the same bound independently.
  - C2 (`codex/C2_kooSaleur_{sol,astra}.out`; code and data in `codex/work_C2_{astra,sol}/`): the numbers agree on every common entry I compared. Sol read the data as already ∝ 1/N and assumed conserved-current models were safe. Astra corrected both points, and §3.3 follows Astra.
- `routes_draft.md` holds my working derivations for §3.1, §3.2, §3.4 and §3.6.

## 7. Astra campaign, 2026-10-21 (five parallel runs, gpt-6-astra high, closed-book)
Specs are `specs/C4_{P1,P2,D1,D2,D3}.md` (shared `_background.md`). Outputs are in `codex/C4_*_astra.out`, with full reports in `codex/work_C4_*/`.
Verdict: neither a proof nor a counterexample. Nothing below has been refereed or rechecked by me beyond reading it, except where stated.

- **P1, general proof.** Not proved.
  - Proved: SU(d)-symmetric bonds reduce to a + bP, so ULS is the only conformal case. Chains whose IR stress tensor is exhausted by the Sugawara tensor of an on-site reductive algebra, embedded simultaneously in SU(d), satisfy c ≤ d − 1 given level ceilings k_i ≤ I_i.
  - New diagnostic: the inverse moment F_h(q) = ∫S_h/ω² dω of the energy-density spectral function tends to C|q|/(24π), with velocities cancelling (I checked the normalization). It is not invariant under h → h + λ(A_j − A_{j+1}), which raises it without bound at the same order, so the raw version cannot be bounded.
  - Missing lemma: a bound of (d − 1)/(24π) on an improvement-independent version. P1 notes this is as strong as the conjecture.
- **P2, structured classes.** Not proved, no violation.
  - Settled: SU(2)-invariant spin-½ chains.
  - New sufficient condition: ind(R) ≥ s(R), the sum of Dynkin labels, so Σ k_a ≤ s(R) implies c ≤ d − 1. This covers rectangular SU(N) fusion chains if their IR is SU(N)_k.
  - Yang–Baxter constraint: invariant orthogonal/symplectic bonds bP + gK must satisfy bg[2b + (d−2)g] = 0 or bg[2b − (d+2)g] = 0, leaving only the permutation, Temperley–Lieb and rational branches.
  - Graded-permutation chains have quadratic branches, so they are not CFTs.
  - Missing: an NN-specific upper bound on IR levels, and control of all massless TBA sectors.
- **D1, integrable and non-compact candidates.** No counterexample for c or c_eff.
  - The reconstructed Izergin–Korepin (d = 3) and alternating six-vertex Hamiltonians have bulk non-Hermitian terms and complex finite-size eigenvalues, so they lie outside the conjecture's unitary scope. The black-hole loophole therefore does not apply at the tested parameters.
- **D2, all NN qubit chains (8 parameters).** No reliable c > 1.
  - Run: 1,338 global samples plus targeted and optimization runs. The raw maximum was c_eff = 2.70 at L = 10, a finite-size artifact.
  - Of 16 candidates followed to L = 22, the best entropy slope was 1.70, falling to 1.35 once short blocks are excluded. None is a clean CFT with c > 1, but not every above-one point reached L = 22.
  - The prescribed estimators also misbehave on Ising and XX controls at these sizes.
  - A SLURM continuation is prepared (`work_C4_D2/cluster_array.slurm`) but was not submitted.
- **D3, d = 3 and d = 4 symmetric families.** No verified violation in 3,526 calculations.
  - An Sp(4) energy estimate of c ≈ 4.72 is consistent with c = 3 split across two velocities (the energy estimator measures Σ c_a v_a).
  - Small-angle Sp(4) has steep entanglement growth with roughly 1/L² gaps, which looks non-conformal but is unresolved.
  - A 154-job Typhon continuation is prepared, not submitted.

Net: the two proof attempts both end at the same missing ingredient found in §3–§4, a microscopic invariant linear in d (an improvement-independent stress response, or an NN-specific level ceiling). The searches leave the conjecture intact, with two loose ends that only larger sizes can settle: the D2 qubit candidates with entropy slope 1.35 to 1.7, and small-angle Sp(4).

Closed 2026-10-22 at Jonah's request: further numerics are unlikely to settle the question either way, so the D2/D3 Typhon continuations were not submitted.

## 8. Interacting generalizations of the free-fermion bound (2026-10-22, Astra high, closed-book)
Specs `specs/C5_F1.md` (partons) and `specs/C5_F2.md` (Green's functions, RG, Bethe). Reports `codex/work_C5_F1/f1_results.md` and `codex/work_C5_F2/f2.md`. I read both and checked the items marked [checked].
log₂ d cannot survive interactions (SU(3)_1 has c = 2 > log₂ 3), so the target was an interacting invariant that reduces to rank A(1)/2 and is ≤ d − 1.

Two invariants come out, both ≤ d − 1 by on-site linear algebra alone, both equal to c only when the IR is "complete" in their modes.
- Affine dimension of the local occupation code (F1, Theorem 4). Write the d allowed states as occupation vectors of partons and gauge every linear constraint that is constant on them. The coset then has c equal to the affine dimension of the code, which is ≤ d − 1. For the minimal code (one parton of d flavours per site) this is d − 1, saturated by ULS. It assumes the IR is the torus coset.
- Rank of the susceptibility matrix of commuting on-site charges (F2, Theorem 6). rank χ ≤ dim span{q_α mod 1} ≤ d − 1 holds unconditionally [checked: identity directions only shift e(μ) linearly]. If the full IR is a Gaussian (multi-component Luttinger) theory whose modes those charges resolve, then c = rank χ. It covers spinful Luttinger liquids (2), Mott Hubbard (1) and SU(d)_1 ULS (d − 1). It fails for SU(2)_2 (rank 1, c = 3/2) and Ising (rank 0), which carry neutral non-Gaussian sectors.

What fails, with proofs.
- Gauging one U(1) is not hard projection (F1, Prop. 3) [checked]. Four partons with charges (1, 2, 4, 7) constrained to total charge 7 leave d = 2 states, {7} and {1,2,4}. The formal coset U(4)_1/U(1)_70 has c = 3 and keeps the current 2j₁ − 3j₂ + j₃, whose density vanishes on both states. So "c_free − c_gauge" is not a bound without completing the constraints.
- The bond rank does not survive projection. Inter-site hopping projects to zero on the minimal code, so the physical exchange is quartic and the determinant argument has nothing to act on.
- Green's-function zero counts are not complete (F2, Theorem 1). A large-U spinful fermion chain with only spin exchange (d = 4, NN) has charge-gapped, k-independent electron Green functions. Depending on the couplings, its spin sector is a c = 1 CFT or gapped, with the same Green-function data. Fixed-Matsubara windings vanish identically by spectral positivity. I followed the construction but did not rerun its numerics.

Conditional theorems (F2).
- RG continuity: chains in a "continuum-matching" neighbourhood of a quasi-free critical chain have c ≤ c_free ≤ log₂ d. The neighbourhood is not shown to contain an open ball, and ULS-type chains cannot be relativistic descendants of a complete free parent.
- Regular TBA: c equals the number of occupied rapidity intervals. So c ≤ d − 1 for nested gl(d) chains with one interval per level and no massless strings. Takhtajan–Babujian violates the no-string hypothesis.

Net. The interacting generalization of "channels through one bond" is "independent on-site densities", and every version is bounded by d − 1 (or log₂ d in the quasi-free case) through the same rank of commuting on-site observables. What no route controls is a neutral sector: Majorana, parafermion or coset modes not resolved by any on-site U(1). Both reports independently name a locality-based budget for neutral stress-tensor sectors as the missing lemma. That is the general conjecture in disguise.

## 9. Campaign of 2026-10-27: bootstrap tools, d = 2, and whether c is bounded at all
Jonah's framing: any upper bound c ≤ f(d), however weak, would be interesting.

- **Tools** (`tools/`, three Claude reading agents). The entanglement-bootstrap program and McGreevy's notes give velocity-free normalizations, exact multi-scale relations and a finite-dimension no-go template. Every quantity in them that equals c is entropic or controlled by entropies, so each tool gives only c ≲ (onset length of conformal behaviour) × log d. The unexplored lever is building Kac–Moody currents from the wavefunction (2403.18410 §8.3, 2206.01656).
- **d = 2 classification** (C6, Astra). c ≤ 1 is proved on the free-fermion slice and holds, given standard input, on the U(1) and XYZ slices (codimensions 5, 5 and 6 in the 8-parameter family). The generic residual set (a single Z₂, only antiunitary, or no symmetry) is open. Exact diagonalization up to L = 22 found no consistent c > 1, and two apparent values above 1 were estimator artifacts.
- **Any bound** (C7, Astra; items marked [checked] were verified by me):
  - Modular invariance forces at least 1 + (πe/12) c e^{πc/6} states with Δ ≤ c/4 + 2 [checked].
  - If those states appear on a ring of length at most K(d), then c ≤ 1 + (6/π) K(d) log d [checked]. This is the onset-length problem in spectral form.
  - c < d² − 1, with no level ceiling, whenever the IR stress tensor is built entirely from currents of on-site symmetries [checked].
  - An explicit d = 4 translation-invariant NN chain selects any spontaneous period b ≥ 5, with L-independent couplings and a gap for each fixed b. The cell energy is a quantum-walk energy −2t cos(π/b), which evades the classical bound period ≤ d [checked]. With a spectator qubit (d = 8) each cell carries 2^b degenerate states.
  - First-order coupling between blocks through one microscopic bond cannot produce an SU(N) swap with N > d, by operator Schmidt rank [checked].
- **Unbounded c** (Claude agent, still running when this was written; draft in `contrib/unbounded_c/`, unchecked).
  - No construction and no bound exist in the literature. Fixed-coupling translation-invariant constructions give either one global clock with gap 1/poly(L) or blocks that never talk to each other.
  - Stability of a large-c target is not the obstruction: Narain CFTs with no relevant operators exist for every c ≥ 31 (Siegel–Weil count, reproduced by the agent).
  - The obstruction is bond rank. First order cannot exceed the connected operator Schmidt rank r(h) ≤ d² − 1. Higher orders can exceed it, but the non-commuting terms come later than dominating diagonal ones. The agent conjectures c ≤ r(h).
  - That obstruction still leaves a window d − 1 < c ≲ (d² − 1)/2 reachable at FIRST order, because fermion-like channels cost rank 2 per unit of c while ULS costs about c². Exact unit-norm Majorana compressions are limited to log₂ d, but scaled ones are not: for a qubit, three anticommuting Majoranas at scale 1/3 compress completely positively [checked].
  - The agent has been redirected to build an explicit (d, h) with c > d − 1 in this window. The hurdles are dimension accounting for the selector (about q ≥ 8 spectators per site, so d ≈ 32), a block code whose two boundary sites compress to many anticommuting Majoranas, tuning to criticality, and commensurability on rings.
- **Still running when this was written:** d = 2 numerics (DMRG at L = 96 and Typhon jobs on the C4-D2 candidates, `contrib/d2_classification/`), the d ≥ 3 counterexample search (`contrib/counterexample_search/`) and the unbounded-c agent. Their folders will be committed when they finish.

Net. Whether c is bounded at fixed d is now a sharp dichotomy. A bound needs c-many conformal levels to appear on rings of size K(d). A counterexample to c–d needs a first-order, rank-efficient coupling between spontaneously formed blocks. The second route currently looks the more promising.

### 9.1 Update, 2026-10-27 evening: the unbounded-c investigation finished (`contrib/unbounded_c/`)
- **G1, second-order exchange through a tiny port** [checked: I reran the script from scratch]. A cell is a port C^d tensored with an auxiliary register, with code G|i⟩ = d^{-1/2} Σ_a |a⟩|a,i⟩. Every traceless port operator then compresses to zero, and each one sends the code into its own excited "syndrome" space. Engineering the cell's excited spectrum makes the second-order effective coupling between two cells exactly aI + bP_N, for every N and every d ≥ 2. The examples (2,3), (2,4), (3,4) were confirmed by full diagonalization: the swap coefficient over λ² converges to the prediction, the splitting scales as λ², and the residual as λ³. So neither bond rank nor perturbative order obstructs large c, once a cell's internal Hamiltonian is unrestricted.
- **The first-order Majorana route** [unchecked]. It needs a fermion-parity grading on the code. Anticommuting compressions alone give a Heisenberg chain (AKLT edges yield c = 1, not 3/2). The economical candidate is q = 7 at d = 22, with c = 24 > 21, but it needs about 48 tuned masses.
- **Agent's lean** [opinion]: modest confidence that c is unbounded at fixed d and that c–d fails at some finite d₀.
- **The remaining crux.** It is realizing a prescribed cell (its code, its port compressions and its excited spectrum) as the low-energy structure of a translation-invariant NN chain that uses the same h on every bond. Astra C9 (`specs/C9_cell_realization.md`) is working on exactly this. Astra C8 (autonomous, resumable) is also running.

### 9.2 Update: Astra C8, autonomous run 1 (`codex/work_C8/`, gpt-6-astra high; items marked [checked] were verified by me)
- **Joint boundary compression.** For every m ≤ q² − 1 there is a completely positive map from a cell's two q-dimensional boundary sites onto one code, sending the left port to αΓ_a and the right port to αiΓ_aP (α = 1/(4m√q), P the parity) [checked: the right images are Hermitian because P anticommutes with Γ_a; Choi positivity is checked numerically by Astra]. The port coupling ΣF_a⊗F_a compresses to free Majorana hopping, one Majorana per cell per flavour. That is automatically critical, with c = m/2 (Astra's numerics give c/n = 1.000000). With the C7 selector, q = 9 would give d = 36 and c = 40 > 35. The missing ingredient is a block Hamiltonian made of the same repeated h whose ground space is this code.
- **G1 extends to identical cells** with distinct left and right ports. That is the form translation invariance needs. Astra's explicit N = 3, q = 2 cell has dimension 48, and the error is 8 × 10⁻¹⁶.
- **A stable target with no symmetry.** A multi-flavour Luttinger liquid with incommensurate Fermi points has no relevant translation-invariant perturbations for n ≥ 9, so it needs no exact large symmetry. This uses standard bosonization input, not a lattice theorem.
- **No-gos.** A one-pass history clock cannot realize the code. Exact SU(N) protection of a swap through a q-dimensional port is impossible at every order for N > q.
- **Rank conjecture likely false.** Declaring each G1 cell a single physical site gives a fixed bond rank with c = N − 1, unbounded. This is conditional on the effective theory controlling the IR, and it does not affect c–d.
- **Status.** Both Astra lines (C8, and C9 on cell realization) now point to the same missing step: a local, repeated-h realization of a prescribed cell code, plus IR control. C8 has been relaunched from its progress file.

### 9.3 Update: Astra C9, cell realization (`codex/work_C9/`, gpt-6-astra high, unchecked except as noted)
- **Any logical bond through two qubit ports at second order.** With an exactly degenerate N-dimensional code, vanishing first order and a finite cell gap, any Hermitian logical nearest-neighbour bond can be synthesized, including an antiferromagnetic swap with positive coefficient, in a cell of dimension D = 4⌈(6N³ + 13N)/4⌉. This generalizes G1. It still uses an arbitrary matrix cell Hamiltonian.
- **SU(N) protection no-go** (independently of C8). With on-cell SU(N) symmetry, fundamental codes and q-dimensional ports, N > q forces every logical interaction to be scalar, at all orders.
- **Majorana route.** Explicit joint two-port encodings exist, but the route needs exact physical parity, an even number of channels and protection beyond second order. A maximal ungraded three-channel example already violates logical fermion parity at second order.
- **The bottleneck, sharpened.** A programmable cell must reproduce the gadget's boundary matrix elements across the whole relevant excited band, not just its spectrum. A Feynman–Kitaev history-state cell dilutes them by the clock length. IR stability is a separate open problem: either protect the target beyond perturbative order, or simulate a lattice target that is stable without symmetry (§9.2).

### 9.4 Update: d ≥ 3 counterexample search finished (`contrib/counterexample_search/`)
- **No counterexample** in the spin-orbital (302 points), spin-3/2 (400 plus a 129-point chart), SO(5)-vector and small-angle Sp(4) scans. Every apparent excess was traced to rotor towers, dimer doublets, ferromagnetic multiplets or incommensurate sizes. The calibrated half-chain estimator is within 0.07 on every control. The small-angle Sp(4) region is dimerized. G2 at d = 7 is inconclusive.
- **Sp(4)_k excluded** [checked]. Pure Sp(4)_k cannot be the IR of an Sp(4) chain on ℂ⁴ at any level. Under the principal SU(2) the 4 is spin 3/2, with Dynkin index 5 against 1/2 for the fundamental, so Sp(4)_k ⊃ SU(2)_{10k}, an even level that cannot carry the half-integer-spin LSM anomaly.
- **Sharpest target so far.** The c = 3 phase of spin-3/2 chains is SO(5)₁ × Ising, with a singlet Majorana whose mass translation forbids. By the Grover–Sheng–Vishwanath and Rahmani et al. mechanism, such a protected Majorana can pass through a tricritical (emergent-supersymmetry) point on the way to dimerization. If the SO(5)₁ sector stays gapless there, the total is c = 5/2 + 7/10 = 16/5 > 3 = d − 1 at d = 4, which would be a counterexample. The scanned part of the phase boundary shows dimerization instead (iDMRG D ≈ 0.107), and the rest is unexplored. A dedicated search is running.
- **Direct sums** give c_max(d₁ + d₂) ≥ c_max(d₁) + c_max(d₂) + 1, which saturates d − 1 but cannot exceed it. Any violation at some d₀ would propagate to every larger d [heuristic].

### 9.5 Update: Astra C8, autonomous run 2 (`codex/work_C8/`, unchecked)
- **Exact codes ruled out.** Exact product codes cannot produce the parity-dressed hopping target when cell endpoints are disjoint. A leakage bound and a tradeoff between gap and cell length explain why the perturbative route behaves differently and remains open.
- **Second-order stability.** The full flat-parent second-order calculation has an attractive singlet mass-channel coefficient, which would gap the target. An explicit parity-preserving term repairs it.
- **Ring lengths.** Conserved markers force at least L ground states on large prime rings, which breaks the CFT description there. A nearest-neighbour marker hop −u(|#L⟩⟨R#| + h.c.) at the crossing of two preferred cell lengths makes the cell lengths a free XX chain at first order instead: a gapless c = 1 "length fluid", with L·gap ≈ 2.39 to 2.31 for L = 33 to 66, matching the free formula to 3 × 10⁻¹⁴. This offers a way around the frozen-defect degeneracy, provided the logical state can travel coherently with the moving markers, which is not yet constructed.

### 9.6 Update: programmable cell (`contrib/programmable_cell/`)
- **Ports need not be diluted** [checked]. If the cell's history-state computation never gates the port, and the port enters only through energy terms, every port operator commutes with the circuit and maps the history subspace into itself. The code and its port compression are then exact, and the two-cell exchange converges to the target as 1/J_c in the clock coupling (Astra P1, gpt-6-astra high, with the agent's own independent script). A single gate on the port leaves an error that survives J_c → ∞. The published 1D translation-invariant universal Hamiltonians (KPBC, 2003.13753) likewise keep inputs in place and remove the dilution by idling [read, by the agent].
- **An exactly degenerate first-order Majorana cell** [checked]. H = Δ(1 − |Ω⟩⟨Ω|) + g Σ_a (|Ω_a⟩⟨Ω| + h.c.)⊗γ_a closes, for each logical vector v, on the 2×2 block [[0, g√m], [g√m, Δ]], independently of v. Its 2^n-dimensional ground space is therefore exactly degenerate for every g. With channels odd under a port grading, both ends compress to C8's joint form αγ_a and αiγ_aP (the agent checked this up to q = 4).
- **Resulting outline** [unchecked]. A translation-invariant chain at the level of tracks and terms, with d₀ = 4 k_c q, where k_c is the alphabet of a segment-confined computation track (guessed at 10² to 10³). The second-order route gives c = N − 1 at d₀ = 8 k_c, independent of N, which would make c unbounded at fixed d₀. The first-order route gives c = q²/4 > d₀ − 1 once q > 16 k_c.
- **No obstruction theorem of the naive kind.** Boundary-visible complexity cannot be bounded by a function of d alone (Astra gives a d = 4 example, checked by the agent). The only true bounds are tradeoffs between gap, scale and size, and these bound c only through the conformal onset scale, as in the C7 dichotomy.
- **Open.** (i) Writing out the confined computation track, which looks routine but long. (ii) Infrared control of the logical chain at finite coupling, which is decisive. (iii) Ring lengths, where C8's mobile selector is a candidate fix.

### 9.7 Update: Astra C8, autonomous run 3 (`codex/work_C8/` §15–19, unchecked)
- **Exact mobile selector.** An exact d = 4, c = 1 mobile selector works on every ring length (equal-hopping selector). Entropy and Casimir estimates give c ≈ 1.004 and 0.9999997.
- **Padding transport.** Padding a two-port code with vacuum symbols lets the logical state move coherently with the markers (checked by Astra to 5 × 10⁻¹⁶ on a small code). An exact cat-padded code cannot be the sole ground space of a nearest-neighbour parent. A vacancy-hopping repair works for any frustration-free core with a uniform positive bond (checked on a three-qubit core).
- **Finite-witness theorem.** By the Hilbert basis theorem, for each q there is a finite K_FF(q) such that a uniform positive nearest-neighbour bond with a zero-energy open-chain state at length K_FF(q) has one at every length. Selecting arbitrarily long cells as the longest frustration-free segment is therefore impossible at fixed q, which gives c_target < K_FF(q) log₂ q for that particular selection route. It does not cover C7's kinetic-energy selection or boundary computations.
- **Still missing.** A high-channel frustration-free uniform core, selection of the right occupation-number sectors on a full ring, and control of the coupled length–payload theory at finite coupling.

### 9.8 Update: Astra C8, autonomous run 4 (`codex/work_C8/` §20–24, unchecked)
- **Boundary Schmidt-rank obstruction.** If a cell's code contains even one product vector across the boundary, the first-order Clifford construction is capped at c_target ≤ q − 1. Viable codes must therefore be entangled throughout. The earlier Choi code has full boundary rank, checked analytically.
- **Sharp witness length.** K_FF(2) = 4 for the frustration-free selection route (§9.7).
- **Conditional bound.** c ≤ ⌊r(h)/2⌋ for Gaussian modes resolved by commuting microscopic charges.
- An explicit graded-swap construction reaches only the existing channel budget. The verdict remains nonconclusive.

### 9.9 Update, 2026-10-27 evening: d = 2 settled, the d = 4 target negative, IR control at the physics standard
Jonah's standard for this question is theoretical physics (standard RG, controlled perturbation theory, numerics calibrated with at least two estimators), not mathematical rigor.
- **d = 2: no c > 1 anywhere we looked** (`contrib/d2_followup_220_1282/`) [checked: I recomputed the finite-entanglement slopes from the raw data]. The two flagged points are incommensurate "floating" spirals with c = 1. Finite-entanglement scaling S against log ξ gives 1.00 (point 220, ξ up to 730) and 1.02 (point 1282, ξ up to 90). The same pipeline gives Ising 0.50, XX 1.00, XXZ 0.99, cluster-Ising 1.5 and two decoupled XX chains 2.08. The earlier doubling estimates near 2 were open-boundary artifacts.
- **d = 4 tricritical target: not found** (`contrib/tci_target_d4/`, numerical, unchecked except for the operator count). The glued Sp(4)₁ × tricritical Ising point (c = 16/5) has exactly two relevant operators invariant under the lattice symmetries, ε′ and Mε, both of dimension 6/5 [checked: ε and the Sp(4)-singlet mass M each carry momentum π, so only ε′ and their product are translation invariant]. It is therefore codimension 2, in principle reachable in the two-parameter spin-3/2 family. The full map of the c = 3 phase shows no tricritical fingerprint: the soft singlet's partners sit at free-Majorana ratios of 2 to 4, never 6. Every c_half above 3.15 sits next to a ferromagnetic boundary or a forming dimer doublet. Estimate now about 3%.
- **IR control: standard physics, no obstruction found** (`contrib/ir_control/`, unchecked). The incommensurate multi-flavour target needs an on-site Z₂ acting as logical fermion parity. It falls within the multi-channel constructive-RG framework of Mastropietro–Porta (2022) and Porta–Scola–Singh (2025), up to routine differences (no U(1), general even interactions). Incommensurate detuning removes every non-Gaussian marginal channel. The SO(2n)-symmetric first-order point and SU(N)₁ instead sit on separatrices. Local Schrieffer–Wolff holds to any finite order, with exponentially local remainders. A rigorous spectral central charge is out of reach for any non-integrable Luttinger liquid, which is a limitation of the field rather than of this construction.
- **Running:** C10 (explicit computation track and d₀), C11 (route-B Gross–Neveu sign), encoded-chain numerics at finite coupling, a red-team review, and C8 run 5.

### 9.10 In progress, 2026-10-27 ~20:00 (Astra drafts, unchecked)
- **C11, route-B infrared** (`codex/work_C11/RESULT.md`, finished; the m = 7 sign flip is unchecked by me; the full second-order calculation passed Astra's two- and three-cell spectral checks with O(λ³) residuals; r = 2 has a stable weak-coupling c = 1 region; for q ≥ 8 no choice of g/Δ repairs the sign). The four-Majorana interaction starts at order λ²g⁴/Δ⁵. In the flavour-symmetric weak-dressing family q₊ = 1, q₋ = m, the non-abelian Gross–Neveu coupling flips to the marginally relevant sign at m = 7, and a common chemical potential does not cure it. So the symmetric Majorana route does not give stable large c. A conditional fix keeps exact U(1)^m flavour charges with distinct small chemical potentials, which makes the exchange and umklapp processes oscillate and leaves an open multi-component Luttinger region with c = m. That is consistent with the incommensurate prescription of `contrib/ir_control/`. The second-order answer also depends on the cell's excited-state resolvent, not only on the port compressions.
- **C10, computation track** (`codex/work_C10/RESULT.md`, draft). It gives an explicit translation-invariant NN track for any hard-coded finite gate program, with path-clock spectrum, initialization and illegal-sector bounds, and preservation of the selector ordering. But the alphabet grows as k_c = 2m + T + 2 with cell size and program length. A fixed d₀ therefore needs a fixed-alphabet decoder of the Gottesman–Irani or KPBC type, which is not supplied. Energy terms at chosen clock times are not exact ground states at finite clock coupling, consistent with the O(1/J_c) convergence of §9.6.
- **Prior work** (`contrib/prior_work/`, finished; literature claims tagged [read] or [recalled] by the agent, not reopened by me). The claim that c is unbounded at fixed d appears to be new, and nothing found contradicts it. Translation-invariant universal-Hamiltonian constructions are all finite-size with one global clock, and Kohler et al. (2003.13753) list boundary conformal symmetry as open. Zhou–Aharonov (2102.02991) list fully translation-invariant strong universality as open. González-Guillén–Cubitt (1810.06528) show history-state gaps close at least as fast as 1/n, so the logical bandwidth shrinks with N, consistent with the hierarchy. SU(N) t–J (U → ∞ Hubbard) chains give c = d − 1 with mobile carriers, another saturating family. The quantum Wielandt inequality is the published analogue of the finite-witness theorem (§9.7). There are no new citations of Latorre–Sierra.
