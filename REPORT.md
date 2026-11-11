# c–d conjecture REPORT (2026-10-21): can c ≤ d − 1 be proved?

Paper: J. I. Latorre and G. Sierra, "The c-d conjecture", arXiv:2403.17242 (v1 March 2024, v2 Nov 2024), J. Stat. Mech. (2024) 113103.
Inputs: both arXiv versions, read in full. Literature map `lit.md` (Opus agent, every claim tagged [read] or [recalled]).
Codex: C1 is the free-fermion theorem and C2 is the lattice-Virasoro numerics, specs in `specs/` and transcripts in `codex/`, all closed-book on gpt-6-astra high (§6).
Nothing here is refereed. §6 lists what I checked.

## Bottom line (as of 2026-10-21, superseded)
This section, §4 and §5 predate the cell construction of §9. The current state is in [STATUS.md](STATUS.md).

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
- **d = 4 tricritical target: not found** (`contrib/tci_target_d4/`, numerical, unchecked except for the operator count). The glued Sp(4)₁ × tricritical Ising point (c = 16/5) has exactly two relevant operators invariant under the lattice symmetries, ε′ and Mε, both of dimension 6/5 [checked: ε and the Sp(4)-singlet mass M each carry momentum π, so only ε′ and their product are translation invariant]. It is therefore codimension 2, in principle reachable in the two-parameter spin-3/2 family. The full map of the c = 3 phase shows no tricritical fingerprint: the soft singlet's partners sit at free-Majorana ratios of 2 to 4, never 6. Every c_half above 3.15 sits next to a ferromagnetic boundary or a forming dimer doublet.
- **IR control: standard physics, no obstruction found** (`contrib/ir_control/`, unchecked). The incommensurate multi-flavour target needs an on-site Z₂ acting as logical fermion parity. It falls within the multi-channel constructive-RG framework of Mastropietro–Porta (2022) and Porta–Scola–Singh (2025), up to routine differences (no U(1), general even interactions). Incommensurate detuning removes every non-Gaussian marginal channel. The SO(2n)-symmetric first-order point and SU(N)₁ instead sit on separatrices. Local Schrieffer–Wolff holds to any finite order, with exponentially local remainders. A rigorous spectral central charge is out of reach for any non-integrable Luttinger liquid, which is a limitation of the field rather than of this construction.
- **Running:** C10 (explicit computation track and d₀), C11 (route-B Gross–Neveu sign), encoded-chain numerics at finite coupling, a red-team review, and C8 run 5.

### 9.10 In progress, 2026-10-27 ~20:00 (Astra drafts, unchecked)
- **C11, route-B infrared** (`codex/work_C11/RESULT.md`, finished; the m = 7 sign flip is unchecked by me; the full second-order calculation passed Astra's two- and three-cell spectral checks with O(λ³) residuals; r = 2 has a stable weak-coupling c = 1 region; for q ≥ 8 no choice of g/Δ repairs the sign). The four-Majorana interaction starts at order λ²g⁴/Δ⁵. In the flavour-symmetric weak-dressing family q₊ = 1, q₋ = m, the non-abelian Gross–Neveu coupling flips to the marginally relevant sign at m = 7, and a common chemical potential does not cure it. So the symmetric Majorana route does not give stable large c. A conditional fix keeps exact U(1)^m flavour charges with distinct small chemical potentials, which makes the exchange and umklapp processes oscillate and leaves an open multi-component Luttinger region with c = m. That is consistent with the incommensurate prescription of `contrib/ir_control/`. The second-order answer also depends on the cell's excited-state resolvent, not only on the port compressions.
- **C10, computation track** (`codex/work_C10/RESULT.md`, draft). It gives an explicit translation-invariant NN track for any hard-coded finite gate program, with path-clock spectrum, initialization and illegal-sector bounds, and preservation of the selector ordering. But the alphabet grows as k_c = 2m + T + 2 with cell size and program length. A fixed d₀ therefore needs a fixed-alphabet decoder of the Gottesman–Irani or KPBC type, which is not supplied. Energy terms at chosen clock times are not exact ground states at finite clock coupling, consistent with the O(1/J_c) convergence of §9.6.
- **Prior work** (`contrib/prior_work/`, finished; literature claims tagged [read] or [recalled] by the agent, not reopened by me). The claim that c is unbounded at fixed d appears to be new, and nothing found contradicts it. Translation-invariant universal-Hamiltonian constructions are all finite-size with one global clock, and Kohler et al. (2003.13753) list boundary conformal symmetry as open. Zhou–Aharonov (2102.02991) list fully translation-invariant strong universality as open. González-Guillén–Cubitt (1810.06528) show history-state gaps close at least as fast as 1/n, so the logical bandwidth shrinks with N, consistent with the hierarchy. SU(N) t–J (U → ∞ Hubbard) chains give c = d − 1 with mobile carriers, another saturating family. The quantum Wielandt inequality is the published analogue of the finite-witness theorem (§9.7). There are no new citations of Latorre–Sierra.
- **C8 run 5** (`codex/work_C8/` §25–29, unchecked). A nearest-neighbour construction at local dimension 6q + 3 turns any initialized frustration-free core into complete ground codes at two adjacent cell lengths, with explicit sector and gap bounds, which is the mobile-cell ingredient with a fixed alphabet. Nonnegative joint Clifford codes support quadratically many channels and independently tunable flavour densities, which is the cell-level ingredient the U(1)^m incommensurate repair of C11 needs. Still open: a high-channel spatial core, the periodic implementation and infrared control.
- **C10, computation track** (finished, unchecked). The explicit partial construction for a small hard-coded program has k_c = 24, q = 2 and d₀ = 192, with proved illegal-sector and selector-margin bounds. Two gaps remain: a fixed-alphabet compiler (so d₀ does not grow with the program) and exact history ground states at finite clock coupling, which the explicit route-B schedule fails, as shown algebraically and numerically. At the physics standard, the O(1/J_c) deviation is acceptable if it only renormalizes the effective logical couplings. That is to be checked against the infrared analysis.
- **C12, a simple-counterexample search at moderate d** (in progress, unchecked). A qutrit (d = 3) candidate: a three-fluid mixture of "atoms" and bound pairs in two momentum valleys ±Q, with an explicit 9×9 bond (`codex/work_C12/DERIVATION.md`). Periodic L = 21 entropy fits give c ≈ 3.05 (von Neumann) and 3.03 (Rényi-2) at δ = 0.9, which would exceed d − 1 = 2. But the small-momentum gaps are equally consistent with z = 2, the open chains are still converging and show less entropy, and valley miscibility is unsettled. An independent verification on Typhon is running (`contrib/qutrit_candidate/`).
- **d = 4 tricritical target, final** (`contrib/tci_target_d4/`, numerical and negative). The glued Sp(4)₁ × tricritical Ising point needs two tunings (ε′ and Mε, both Δ = 6/5; field theory by Astra, rederived and checked by Claude), so on the two-parameter spin-3/2 SU(2) family it can only be isolated. ED at L = 8, 12 over the whole two-sphere and 292 boundary points finds no tricritical-singlet fingerprint: the soft momentum-π singlet along the c = 3/dimer boundary has x_eff ≈ 0.5 to 0.7 against 0.18 in a tricritical control at the same sizes, and every c_half > 3.15 sits next to a ferromagnetic boundary or a level crossing. L = 16 ED on three boundary points is still running on Typhon.
- **IR control of the encoded chain, final** (`contrib/ir_control/`; Astra on the operator classification and leakage expansion, checked by Claude). At the physics standard the incommensurate multi-flavour target is under control, provided every term of h is even under an on-site Z₂ that acts on the code as logical fermion parity. Incommensurate detuning removes every non-Gaussian marginal channel, so no mass-channel repair is needed, while the SO(2n)-symmetric point sits on a one-loop separatrix and runs away under small anisotropy (checked for so(4)); a small n = 2 iDMRG toy agrees (the symmetric chain dimerizes, the detuned one does not). A Sidon-set choice of Fermi momenta gives margins δ₂ ≥ π/(6n), δ₄ ≥ π/(576n²). The closest rigorous framework is multi-channel Benfatto–Mastropietro (Mastropietro–Porta 2108.13708, Porta–Scola–Singh 2509.08665), which assumes U(1) and proves correlations rather than c. A Typhon leakage array (17511966) is still running.
- **Minimal d₀, final** (`contrib/minimal_d0/`; Astra on the frame span, checked by Claude, who also reran `frame_checks.py`). The old count d₀ = 4k_c q overcounts: only two letters (port and partner) need the q-dimensional register, and one static Bell-check term on the (port, partner) bond is the only term that must touch the port, since every port coupling either route needs is a positive time-average of it conjugated by gates on the partner side [proved, Astra; checked at q = 2, 4]. On a Bausch–Cubitt–Ozols-class history machine (42 dimensions, the smallest translation-invariant one in the literature) this gives d₀ ≈ 130 for route A with the C7 walker tensored on, and d₀ ≈ 45 if the machine's own clock selects the cell length, which needs a machine whose legal history is the longest penalty-free path in every cell (not yet built). Route B gives d₀ ≈ 76 to 88 (merged) or 174 to 208 (tensored), with c ≈ 81 to 210. These counts are heuristic, rest on adapting BCO's machine to marker-bounded cells, and are for frozen cells; the mobile version adds the equalizer and padding registers.
- **Red team of the construction, final** (`contrib/redteam_construction/`, a Claude agent that reimplemented the decisive checks rather than rerunning them). No fatal flaw in the finite-dimensional ingredients: G1 (residual 2.5 × 10⁻¹⁵), the exact degeneracy of the dressed Bell-pair cell, the port-idle identity, and the second-order formula with the three-cell term (error falling as λ³ against exact diagonalization of three cells) all pass. Ranked problems: (1) frozen cells give at least L ground states on prime rings, fatal under the reading that every ring must be conformal (the mobile-cell repair is in progress, see STATUS step 6); (2) the SU(N)₁ target of route A sits on a separatrix and should be replaced by the incommensurate target, which still gives c = n ≈ log₂ N, unbounded; (3) the computation track is unwritten; (4) no rigorous thermodynamic Schrieffer–Wolff reduction, which matters for a theorem but not at the physics standard. New: the second-order chain of the route-B cell, computed directly, forbids every mass and has a repulsive singlet Gross–Neveu sign, but it is never untuned. Its first-order hopping is purely imaginary, so every flavour's Fermi points sit at δk and π − δk, the inter-flavour pair-tunnelling "quartets" stay at zero momentum, and the flow runs away under small flavour noise. My reading [checked against `contrib/ir_control/astra_A/`]: this is the degenerate case that the incommensurate prescription excludes. With independent hopping phases per flavour (available at first order from C8 run 5) the Fermi-point centres differ, the quartets carry momentum 2(θ_a − θ_b), and the Sidon margin δ₄ covers them.
- **C8 run 6, final** (`codex/work_C8/` §30 to §35, Astra; the five check scripts were rerun by Claude and pass). New explicit pieces: a complete composition of the two-length selector with the cell wrapper, the exact all-ring spectrum of the projected length fluid, an equalizer that makes the first-order logical bond independent of adjacent cell lengths (local dimension 96q + 48), the complete projected spectrum of length fluid plus any payload with one momentum-gluing constraint (c = 2.0000 from the Casimir energy for an XX payload), and error bounds for approximate positive cores. Still open, in Astra's ranking: a uniform high-channel core with a constant alphabet, and the infrared behaviour of the microscopic moving-cell model at fixed coupling (now C13).
- **C12, final** (`codex/work_C12/RESULT.md`, Astra). No counterexample. The qutrit model has an analytically controlled atom band and two bound-pair minima, but those minima disappear in a saddle-node at δ* ≈ 0.570, so the δ = 0.75 and 0.9 points that gave the large entropy estimates do not have the structure the three-fluid argument needs, and their estimators fail a size-scaling check. Cell, ladder, gauge and squeezed-space screens found nothing either.
- **Refined conjecture, final** (`contrib/refined_conjecture/`; f-sum constants and the SU(d) thermodynamic Bethe ansatz by Astra, rederived and checked against exact Schur–Weyl diagonalization by Claude). Bounds of the form c ≤ F(d, window) at the physics level: c ≤ 3ℓ* log d/(π(1−ε)), with ℓ* the length where conformal behaviour sets in; c ≤ 24R²/π with R = ‖h‖/T*, independent of d; and for the Casimir coefficient A = Σ c_a v_a, A ≤ π‖h‖ for number-conserving free fermions and A ≤ 8πC₂(R)‖h‖/(k + h^∨) for current algebras of on-site generators. The SU(d) chain has A = 2π(1 − 1/d)‖h‖ [checked: v = 2π/d for h = P] and a conformal window that closes at T* ≈ 6.1‖h‖/d², with ℓ* ≈ d. Proposed refinement of c–d: A = Σ c_a v_a ≤ 2π‖h‖ for every chain in the class, with no d in it, which together with the entropy bound gives c² ≤ 6R log d/(1−ε). The cell construction is consistent with it, since its velocities shrink with the cell length, so it trades velocity for central charge rather than violating the refined bound [heuristic].
- **C13, mobile cells at finite coupling, final** (`codex/work_C13/RESULT.md`, Astra; the soluble controls rerun by Claude reproduce c = 1.00003, 2.00006, 3.00009 from the Casimir energy and 1.00008, 2.0003, 3.0007 from cell entropy for zero, one and two payload flavours plus the length mode). Conditional positive result at the physics standard: at fixed n and cell length, with exact marker conservation, exact payload parity, displaced and jointly non-resonant Fermi points (including the length fluid's), controlled elimination of gapped cell sectors and a positive stiffness matrix, the mobile construction is an (n + 1)-boson Gaussian phase, and the extra length mode forces no relevant mixed resonance; an explicit open weak-coupling region exists in standard bosonization. Finite-ring uniqueness is proved for a positive-hopping subclass and checked over the stated ranges. Qualifications: (i) at irrational densities the zero modes drift with L, so the spectrum converges to CFT towers in L-dependent charge and flux sectors rather than to one fixed vacuum sector, exactly as for the XXZ chain in a generic field, which is universally counted as c = 1; (ii) the local-boundedness use of the c-theorem near a conformal point needs a uniform crossover hypothesis; (iii) the high-channel microscopic core and its coefficients are still not specified, and generic interacting numerics at accessible sizes do not yet meet the two-estimator standard.
- **Mobile cells, toy model** (`contrib/mobile_cells/`, in progress, Claude agent with one Astra run). A d = 6 chain whose cells carry a spin payload (c_p = 1): the frozen chain has degenerate ground states (3 on rings with 3 | L, up to L otherwise), while the mobile chain has a unique ground state on every ring tested including primes, z = 1 with two velocities, and c = 2 from three estimators. Squeezed-space factorization (Ogata–Shiba) is exact and checked against brute force to 10⁻¹⁴.
- **Encoded route-B chain at finite coupling** (`contrib/encoded_chain_numerics/`, in progress). For q = 2 (r = 2), c = 1.00 ± 0.03 for λ/Δ from 0.05 to 0.8, with z = 1 and no translation breaking, so the first-order prediction c = r/2 survives at finite coupling. The q = 3 case is still running.

### 9.11 Chat material imported and status reconciled, 2026-10-27

Codex (GPT-6) compared the available text of "Prove c d conjecture" against
commit `58ef8a9` and added [contrib/chat_followup/](contrib/chat_followup/README.md)
at the user's request. The import adds the species-torus charge-incidence
argument excluding marginal Cooper transfer, a direct-sum realization of
finite cell excitations by massive auxiliary fermions, explicit closure and
symmetry arguments for controlled classes, and a sourced qubit Hamiltonian
catalog. The original chat's model identity was unavailable. Two construction
replies and the catalog were truncated in retrieval; the contribution states
this and does not reconstruct missing text or references.

The finite-dimensional arguments are separated from the proposed
c = r², d = O(r log r) microscopic family [heuristic]. In particular, the
compiler, initialization/illegal-sector checks and thermodynamic RG transfer
remain to be supplied. The source chat's stronger suggestion that only the
last RG extension remains is not adopted. The earlier C1/C6 results and
generic boundary-compression lemmas were already recorded and are linked
rather than treated as new work. No new numerical evidence for a violation
of c-d is added.

Separately, the already committed [C14 progress log](codex/work_C14/PROGRESS.md)
supersedes the favorable general Casimir-conjecture summary above. It reports
the symmetric-rank-2 SU(8) chain with d = 36, c = 12.6, v = π/4 and
‖h‖ = 3/2, so cv/(2π‖h‖) = 1.05 [Astra; algebraic and TBA checks, standard
Bethe/string and WZW inputs; full write-up and independent review pending].
It also reports a family inconsistent with any dimension-independent
constant in that bound. These statements concern the velocity-weighted
proposal, not the original c-d conjecture. STATUS and the refinement note
now point to this newer result. The on-site-current f-sum argument's actual
constant is 8πC₂(R)/(k+h^∨); it did not establish 2π in general.

### 9.12 C14 final, 2026-10-27
- **C14, the Casimir bound is false** (`codex/work_C14/RESULT.md`, Astra; the bond spectrum and the velocity independently checked by Claude with `contrib/refined_conjecture/check_C14_velocity.py`). The integrable SU(N) chain with the rank-k symmetric representation on each site, h = Σ_j 2(H_{k−j} − H_k)Π_j (for k = 2, h = −X²/8 + 3X/4 − 1 with X = Σ E_ab ⊗ E_ba and spectrum {0, −1, −3}), has ‖h‖ = H_k, c = k(N² − 1)/(N + k) (SU(N)_k) and v = 2π/N independent of k, so A/(2π‖h‖) = k(N² − 1)/(N(N + k)H_k). This exceeds 1 for k = 2 and N ≥ 7 (d = 28, ratio 64/63) and grows like N/(2 log N) for k = N, so no d-independent constant works. Claude's ED gives v = 3.142 at N = 2, L = 10 (target π) and 2.088 at N = 3, L = 9 (target 2.094). Positive by-product: every even quadratic Majorana or BdG chain, any number of orbitals per site, obeys A ≤ π‖h‖ [Astra, physics level]. The family does not threaten c ≤ d − 1 (c < N² while d = C(N + k − 1, k)). For k = N, A/‖h‖ grows roughly like log d/log log d, so a bound of the form A ≤ C‖h‖ log d is still open.


### 9.13 Sparse-graph continuation: a candidate physics-level refutation, 2026-10-28 UTC

At Pace's request, Codex (GPT-6 family) pursued the most promising route
with three parallel research workers, starting from collaborator commit
`5360803`. The result is a [candidate counterexample argument](contrib/research_continuation/PROPOSED_COUNTEREXAMPLE.md)
at this project's stated theoretical-physics standard. Its finite
construction claims are exact; the conformal conclusion uses explicitly
stated local massive-sector integration and stable Gaussian RG inputs.
No rigorous thermodynamic CFT theorem, human review, or numerical
measurement of a violating chain is claimed.

- **Sparse simultaneous ports [proved, independently checked].** For any
  simple bipartite graph with u upper vertices and N edges, an explicit
  isometry gives left ε f† and right ε(1−1/(2N))f†P, with
  ε=1/(2√(uN)). The physical edge-swap bond has norm equal to its largest
  hopping amplitude. The map is species preserving and has a sparse
  reversible basis description. Full small-case isometries reproduce
  both compressions to below 7×10⁻¹⁶. See [exact_route](contrib/exact_route/README.md).
- **Graph protection [proved].** The projective-plane incidence graph has
  q=2(p²+p+1) species and N=(p+1)(p²+p+1) logical fermions. Girth six
  excludes every nonzero neutral flavour transfer with fewer than six
  fields. Thus all flavour-changing quartics are forbidden exactly.
  Enumeration checked 143,528 short charge vectors at p=2; graph and
  momentum checks also cover p=3,5,7. See [torus_bottleneck](contrib/torus_bottleneck/README.md).
- **Local circuit and compiler [proved by constructive composition with
  classical reversible computation].** Local species queries, loaded
  exchanges and neutral reversible operations prepare the code with
  O(q log q) local states. Exact interval amplitudes need only O(log q)
  indexed algebraic rotations, so approximate gate synthesis is optional.
  A classical oracle scheduler has fixed finite control and an input-
  independent trajectory. The moving head overlays one species record,
  avoiding a hidden q² on-site factor. The existing fifteen-letter C8
  clock executes 2^n ordered full sweeps. Independent exhaustive checks
  cover all 50,625 length-four clock words and the quantum sweep blocks.
- **Initialized positive core [proved finite theorem].** Static endpoint
  marks, trusted neutral scratch and a single boundary-pinned query helper
  initialize the computation. A reversible full-format checker and
  all-time FAIL penalties select exactly 2^N inputs. A noncommuting
  acceptance estimate gives gap at least aγ/(a+B+γ). Cyclic preparation,
  flagged readout and exact inverse preparation give exact endpoint
  matrices times positive computable history weights. Compact scans and
  SWAP walks avoid storing a history proportional to n log n in n sites.
  See [EXISTENTIAL_CORE](contrib/torus_bottleneck/EXISTENTIAL_CORE.md).
- **Mobile composition and infrared inference [exact composition;
  conditional physics conclusion].** The existing wrapper preserves the
  species torus and costs d=96d_core+48. Choose
  k_a=π/4+aπ/[12(N+1)] and length density in [5/12,7/12]. Payload
  mismatches are at least π/[6(N+1)], mixed mismatches at least π/6,
  and all other non-Gaussian resonances start at dimension three.
  With both retained energy scales first order, fixed sufficiently small
  coupling puts the higher-order local corrections inside a strict
  Gaussian stability neighborhood. The [infrared audit](contrib/research_continuation/IR_AUDIT.md)
  derives an additive local defect-energy bound and gives the controlled
  physics excursion argument, while explicitly withholding an all-order
  mathematical convergence claim. Under these standard infrared inputs,
  c=N+1 grows as p³ while d=O(p² log p), hence c/d→∞.
- **Route-A repairs [proved finite algebra].** The original diagonal
  syndrome prescription cannot encode odd hopping factors while keeping
  the desired parity. Off-diagonal graded weights repair it. Independent
  left/right dilations can share the logical register with a complete
  gapped product code even when the gadget Hamiltonians do not commute.
  Full resolvent checks agree below 9×10⁻¹⁵, with relative errors reported
  against the small exchange scale. See [graded_route_a](contrib/graded_route_a/README.md)
  and [separated_port_cells](contrib/separated_port_cells/README.md).

All ten recorded check scripts were rerun by the coordinating agent;
independent workers also rederived the key port, acceptance, sweep and
infrared-selection arguments. The exact-amplitude checks cover 256
interval states and 1,024 two-branch preparations, with error below
2×10⁻¹⁵. Another 29,359 verifier-helper round trips, including 18,364
malformed helper-like target records, pass exactly. These checks concern finite ingredients, not an expanded universal
machine or a thermodynamic central-charge fit.

The candidate's status is **complete at the construction-argument level,
conditional on the stated standard infrared inputs, and in need of
independent end-to-end review**. Its compiler constant is not numerically
expanded or optimized. Unlike the older route-A goal, this family does
not prove unbounded central charge at one fixed microscopic dimension.

### 9.14 Route A, route B and d = 2, evening of 2026-10-27 (coordinator; check scripts of every Astra run below rerun by Claude and passing)
- **C19, a uniform route-A cell** (`codex/work_C19/RESULT.md`, Astra). A uniform family of finite matrix cells with two qubit ports, n logical fermion modes, an exactly degenerate 2^n code and exact parity, whose second-order exchange is independent flavour hopping with arbitrary phases and chemical potentials. It needs only n + 1 port-accessible excitation energies, with Jacobi-polynomial coefficients from bounded loop bodies; as a matrix the cell has dimension 16(n + 2)2^n. Two obstructions are proved: a pole-count bound 2n ≤ 4(k − 1) for k accessible energies, and (from C18) the failure of the cheapest frustration-free interface, whose source resolvents are scalar. No local controller and no d₀ yet.
- **C18, clock-selected cell length** (`codex/work_C18/RESULT.md`, Astra). Found the frustration-free interface obstruction above; the clock and interface checks pass. No unbounded-c d₀.
- **C8 run 7** (`codex/work_C8/RESULT.md` §36 to §39). A 15-letter initialized nearest-neighbour clock with one legal history and explicit exclusion of unwanted sectors, checked through 100,000-state ring spectra. Missing: the clean-readout program that produces both high-channel boundary maps.
- **C20, blind check of the mobile-cell theory** (`codex/work_C20/RESULT.md`, Astra, closed book). Without seeing C13, it finds an open stable region with c = n + 1 and momentum margins of order n⁻⁴, independently confirming C13 at the level of standard bosonization.
- **Route A numerics at finite coupling** (`contrib/routeA_numerics/`, Claude agent). G1 identical cells with the SU(3)₁ target look like c ≈ 2 at λ ≤ 0.1 (two estimators, ξ up to about 30 cells), but SU(3) breaking starts at λ³ and the spectrum drifts at λ = 0.2, consistent with SU(3)₁ sitting on a separatrix. An identical cell built for an incommensurate XX target generated a relevant pairing term (about 0.4λ) and gapped the chain, because the Fermi points sat at ±k_F. Route A cells must keep an exact U(1) or displace each flavour's Fermi-sea centre.
- **Route B numerics, final** (`contrib/encoded_chain_numerics/`). r = 2: c = 1 confirmed (0.98 to 1.05) at λ/Δ from 0.05 to 0.8. r = 4: the derived logical chain gives c = 1.97 ± 0.06, with the doping predicted by one loop matching VUMPS (0.061 against 0.062 per cell); the full chain of big cells could not be pushed to long distances.
- **d = 2 program, first round** (`codex/work_C15/` to `codex/work_C17/`, Astra). No proof and no counterexample. Established: exclusions for conserved one-site markers, literal sublattice decoupling and exact U(1) chains; the K ≥ 1/2 lore is not universal (strong repulsion gives K = (1 − ρ)², which can make atom–pair conversion irrelevant), but nearest-neighbour attraction favours larger clusters over a stable adjacent-pair fluid; compactness cannot exclude an isolated finite-c > 1 point. The leading open target is an isolated neutral critical sector (for example Ising) alongside a c = 1 mode at intermediate coupling.
- **Now running on route A**: C21 (compile C19's cell into a constant-alphabet chain with a non-frustration-free interface), C22 (a computation-free spatial cell), C23 (end-to-end assembly with the full scale hierarchy), C8 run 8 (clean-readout table), and a Claude agent building a microscopic end-to-end demonstrator (`contrib/routeA_demonstrator/`). Reviewing the collaborator's species-graph candidate: Astra C24 and a Claude referee (`contrib/review_species_graph/`).

### 9.15 Executable microscopic components and initialization repair, 2026-10-28 UTC

At Pace's request to instantiate the proposed Hamiltonian, Codex and three
parallel workers added [controller_instantiation](contrib/controller_instantiation/README.md).
This is concrete finite implementation progress, not a newly verified
counterexample or thermodynamic central-charge measurement.

- **Encoding circuit [implemented; finite exact/numerical checks].** The
  complete two-port transducer now generates explicit Boolean gates,
  charged routing, algebraic branch preparation and inverse cleanup,
  lowered to 17 two-site families. The p=2 manifest has 21 logical modes,
  317 local record states, 81,342 circuit sites and 287,798 abstract gates.
  Its 1,213,761,392 localized gates were counted, not all simulated. A
  one-edge fixture executes all 36,140 actual local gates coherently and
  restores every auxiliary to blank; amplitude error is below 1.2e-16.
- **Controller tables [implemented; proved local permutations].** A
  preloaded dispatcher has 88 neutral symbols, 150 direction-fixed
  controls and 13,200 collision entries. Its head-overlay factor is
  26,488 times the quantum-record dimension. A separate blank-to-length
  primitive has 12 neutral symbols and 16 controls, computes n-1 and
  reverses, in exactly twice n^2+3n-4-2 popcount(n-1) steps. The latter
  passed complete cycles at all lengths four through 128, totaling
  1,461,216 elementary steps. These are component counts, not the entire
  initialized core's constant C.
- **Input-format repair [proved and checked].** TYPE/species queries
  could not distinguish neutral bit zero from one. The accepted untrusted
  neutral format is now explicitly the unique blank; a full-space
  blank-to-bit0 involution initializes circuit scratch after verification.
  The fixed helper and verifier scratch remain separately pinned. The
  construction notes have been corrected. An independent all-sector
  d=150 two-site history fixture also checks the initialization mechanism
  on all 22,500 basis states and 14,800 propagation edges.
- **Actual physical bond [implemented; independently projected].** The
  uniform bond assembler accepts explicit core/port matrices and includes
  all selector, wrapper, helix, hopping, marker and field terms. A complete
  d=240 one-logical-mode fixture is saved as a 57,600-square sparse bond
  with 122,192 entries. Every column passes Hermiticity and charge checks.
  Independent contraction of the actual hopping entries through cell
  endpoint reductions reproduces all four short/long length combinations
  to relative error below 5.8e-13. Negative controls detect omission of
  X# cuts or endpoint compensators. This fixture is not a c-d violation.

Source-loop accounting also supplies explicit circuit/dispatcher resource
bounds: for A=(N+4q)^2 times the maximum bit width, at most 32A neutral
records, 128A abstract gates, 70,000A^2 local gates, and 2,310,035A^3
preloaded tape sites. For projective planes the latter is below
2^50 p^21; this excludes the autonomous generator and is not a final
selected core length. The coordinating agent reran all nine new check
scripts successfully.
The full autonomous program generator, verified program initialization,
endpoint transport and elementary-clock READY protocol remain to be
composed. The dispatcher demonstrably conserves its program symbols and
therefore cannot generate its own instruction tape from blank. The new
length reader is a separate primitive and needs a reversible handoff.
The full controller count, a numerical violating choice, and the
fixed-positive-coupling infrared conclusion are still outstanding. No
human review is claimed.

- **Casimir survey and literature, after C14** (`contrib/casimir_survey/`, Claude agent with two Astra Bethe-solver runs, in progress; `contrib/lit_casimir_d2/`, Sonnet literature agent, final). Exact diagonalization of the Sym^k SU(N) chains at (N, k) = (2, 1 to 4), (3, 1 to 3), (4, 1 and 2), (5, 1) reproduces the Bethe/WZW values of e₀, v, c and A to about 1%, and at N = 3, 4 confirms that rank two beats the fundamental chain by the predicted margin. The central charge c = k(N² − 1)/(N + k) is in the literature (Führinger et al., 0806.2563); no paper quotes v for general (N, k), but the string density equations are k-independent, giving v = 2π/N. No paper bounds Σ c v or the Casimir coefficient by a local norm. Outside the symmetric family nothing examined exceeds the fundamental SU(d) value of A/(2π‖h‖); across everything computed, A ≤ π‖h‖ log₂ d holds, with equality only for the Heisenberg chain [numerical], and the symmetric family forces any bound A ≤ K(d)‖h‖ to grow at least like log d/log log d. For qubit chains, every published c > 1 example (c = 3/2, 2) needs range three or more; Chepiga–Mila's nearest-neighbour interacting Kitaev chain has c ≤ 1.

### 9.16 Autonomous program/readout cycle and a fixed generator, 2026-10-28 UTC

At Pace's request to continue the most promising implementation direction,
Codex and three parallel workers added
[autonomous_controller](contrib/autonomous_controller/README.md).
The high-channel program remains unfinished; no new c-d violation or
central-charge measurement is claimed.

- **Generated program, execution and readout [implemented; finite checks].**
  One completed local collision permutation now starts with blank dynamic
  tracks, counts the sites, copies n-1 into a binary instruction tape,
  executes noncommuting CV/CNOT gates, exposes two physical endpoint
  records, and uncomputes the generated program and counter. It has 80
  neutral symbols, 40 controls, 3,200 table entries and 497 prescribed
  rules. All 125 lengths 4–128 pass complete cycles, totaling 1,590,466
  moves; 535,782 elementary gate vertices pass endpoint-safety checks.
  This is a demonstration program, not the desired sparse encoding.
- **Initialized microscopic core [finite construction and checks].**
  Full-space local work gates are combined with the fifteen-letter C8
  clock in explicit sparse bulk and boundary matrix-column evaluators.
  A LEFT-controlled seed on the first initialization departure creates
  exactly one head; the dynamic tracks are initially pinned during the I sweep.
  The demonstration's core dimension is 194,400, independent of length.
  Its complete ground space is the history encoding of n free input
  qubits, dimension 2^n. With C=(2n+4)2^n+n-4 and
  gamma=1-cos(pi/C), the gap is at least gamma/(1+C+C gamma).
  An exponential code dimension alone is not a high-channel port theorem.
  The actual clock has positive readout weights from n=7 in the checked
  examples, and a full controller cycle fits for every n>=9; n=4,5,6
  have zero readout weight. At n=9 the endpoint weights are exactly
  658/11269 and 656/11269, independent of quantum input.
- **Full circuit generator [implemented word-machine intermediate].**
  A frozen 3,115-instruction program, with 34 functions and 27 opcode
  families, constructs the graph, circuit, local gates and dispatcher
  tape. Exact abstract comparisons cover all 287,798 operations at p=2
  and all 1,588,624 at p=3. The complete one-edge local/dispatcher streams
  also match. A journal inverse clears its tested computations, including
  adversarial array aliasing and nested calls. Source-loop accounting
  bounds word instructions by 18463(128A)^6 for the earlier parameter A.
  Arbitrary-length arithmetic, addressed arrays and journal operations
  still need compilation to finite local tape rules; the word instruction
  count is not the final local alphabet count.
- **Parameter-input shortcut [explicit gates and circuit checks].**
  Two neutral permutation families expose bit_j(p) and the last index
  through the existing helper H_(j,c). Since the family already has
  p-dependent local matrices, this avoids decoding p from cell length
  without adding a species register or a growing number of controls.
  A 10m-4-gate extraction/cleanup word restores all untrusted records;
  its compact walk still needs finite-table compilation.

Independent AI audit checks the actual initialized clock path, including
carry and idle vertices, and verifies the finite-core initialization and
gap argument. The decisive remaining finite task is compiling and joining
the real generator, format verifier, charged dispatcher and readout into
one counted high-channel core. The standard fixed-coupling infrared inputs
remain separate obligations. All artifacts identify their proof and test
scope; no human review is claimed.

The coordinating agent reran all nine new verification scripts successfully;
the combined result is saved in the contribution's verification summary.
- **Casimir survey, final** (`contrib/casimir_survey/`). The violation of A ≤ 2π‖h‖ is now seen directly on the lattice, without WZW input: finite-size Bethe roots for the rank-two SU(N) chain (solver written by Astra, validated against 11 exact-diagonalization energies to 4 × 10⁻¹¹, rerun and fitted by the Claude agent on rings up to L = 105 and 120) give A/(2π‖h‖) between 1.014 and 1.018 at N = 7 (d = 28) for every fit form, and 1.05 at N = 8 (d = 36). This assumes the continued root configuration remains the ground state. Elsewhere nothing exceeds the fundamental SU(d) value, which a twirling argument shows is a local maximum to first order; the best ratio at each d grows like 0.38 ln d/ln ln d, and A ≤ π‖h‖ log₂ d holds for everything computed.
- **Mobile cells, final** (`contrib/mobile_cells/`, Claude agent with one Astra run). In a d = 6 toy whose cells carry a two-qubit payload (c_p = 1), frozen markers give ground-state degeneracy 3 when 3 | L and 5 to 14 otherwise (L = 9 to 16), while mobile markers give a unique ground state on L = 8 to 18 including primes, z = 1 with both mode velocities measured, Casimir energy over the sum of velocities 1.02 and 0.98, ring-entropy c = 1.90 to 2.11 and finite-entanglement c = 2.04 to 2.12, so c = c_p + 1. Breaking the reflection gaps the payload and leaves c = 1.00 from the length mode alone. The squeezed-space (Ogata–Shiba) reduction is exact and reproduces every eigenvalue at equal hopping to 4 × 10⁻¹⁴; power counting leaves no relevant cross-coupling at generic long-cell density, with isolated commensurate exceptions (a spin-Peierls point at ν = ½, near C8's natural density, so μ must be detuned). Moving a multi-site payload costs one vacancy symbol, d₀ = 4(k_c q + 1). The c_p = 2 extension (c = 3) has a unique ground state but c = 3 is not yet established; Typhon jobs are running.

### 9.16a Route A assembled, and the first review of the species-graph candidate (2026-10-28)
- **Route A assembly** (`contrib/routeA_assembly/`, Claude agent; Astra C23 in `codex/work_C23/` reached the same verdict independently; the coordinator reran `scales.py` and `third_order.py`, reproducing the visibility ρ² = 5.79 × 10⁻⁶ and 9.49 × 10⁻⁸, the clock gaps 1.1 × 10⁻⁶J and 4.9 × 10⁻⁸J, and the third-order coefficients 2.26, 4.47 and 6.52 at n = 1, 2, 3). Nothing found kills route A. Every condition in the scale hierarchy bounds a smaller scale by larger ones, so at each fixed n the scales can be chosen from the top down; nothing has to be uniform in n. Third-order corrections are only a few times ρ²λ³/Δ², allowing λ/Δ around 10⁻², but they include same-flavour pairing of 0.14 to 0.20 λ/Δ (the mechanism behind the pairing gap in `contrib/routeA_numerics/`), and C19's geometric Fermi points are unusable from n = 3 on, so the target must use a Sidon design. The accumulation point is gapless and not conformal (z = 2 wall motion, exponentially large degeneracy). Two architectural constraints: the clock-selected length works only for frozen cells, so mobile cells need C8's wrapper and d₀ above about 2000; and C8's wrapper theorem needs a frustration-free core, so the port interface must keep the code frustration-free (syndrome-mixing positive checks do). Ranked gaps: G1 a fixed-alphabet port interface giving parity-odd, logical-operator-valued second-order exchange with an exact frustration-free code (no-gos: bare Bell checks, unitary-decorated paths and trees, distant product layouts; live candidate: a syndrome-mixing check, being tested in `contrib/routeA_demonstrator/`); G2 the whole-cell compiler; G3 second-order composition with the wrapper; G4 finite-coupling numerics with the Sidon design; G5 C19 with a Sidon target; G6 the second-order infrared transfer; G7 gluing data; G8 the alphabet count; G9 scope.
- **C24, Astra review of the species-graph candidate** (`codex/work_C24/RESULT.md`). No fatal algebraic error, unavoidable relevant resonance or contradiction in the fixed-coupling hierarchy. Accepted: the girth-six selection rule (independently checked for p = 2, 3, 5, 7), the sparse port construction, the charge-preserving circuit and the conditional Gaussian stability argument. Not yet accepted as an established refutation: two serious verification gaps remain, the complete autonomous sparse-encoding controller meeting its initialization and readout contract, and a local low-frequency elimination for the mobile microscopic model with corrections small relative to its first-order Hamiltonian. A Claude referee (`contrib/review_species_graph/`) is reviewing independently.

### 9.17 Exact qubit subclasses and a staggered field reduction, 2026-10-28

At the user's request to prioritize a simplified proof, Codex and three
parallel workers added [d2_gap_reduction](contrib/d2_gap_reduction/README.md).
The general d=2 conjecture remains open. The contribution combines exact
algebra with a recently available rigorous gap theorem; the finite spectra
check conventions rather than establish an unknown CFT numerically.

- **Complete ferromagnetic rank-one class [proved using a cited theorem].**
  For J=−κnnᵀ and arbitrary uniform field B, nonzero B·n gives a unique
  ground state and Δ_L≥2|B·n|. At B·n=0 this is transverse-field Ising,
  with its sole conformal locus |B|=κ>0 and c=1/2. The gap input is
  Bravyi–Gosset–Liu–Wong, arXiv:2607.14401v1, Definition 1 and Theorem 1,
  read directly and independently checked for applicability. More generally
  the same exclusion holds when Jn=Jᵀn=λn and
  λ≤−||(I−nnᵀ)J(I−nnᵀ)||op. Arbitrary transverse fields are allowed.
  This adds genuinely interacting points outside the previous free,
  uniform-U(1), and zero-field-XYZ sets.
- **Exhaustive additive-charge reduction [proved algebraically].**
  For any exact Q=Σ n_j·σ_j, local commutation is equivalent to
  Ω(n_j)J=JΩ(n_(j+1)) and n_j×B=0. Classifying these equations at every
  rank of J shows that zero-field models reduce to XXZ or commuting Ising.
  At nonzero field the only additional interacting class is
  a(XX−YY)+uZZ+bZ, whose charge is staggered and whose alternating-spin
  rotation gives XXZ in a staggered field. This strengthens C6's example
  into an exhaustive classification of this microscopic charge route.
  It does not constrain emergent or multispin conserved currents.
- **Exact gapped regions of that remaining class [proved].**
  On even rings with a≠0, put r=sqrt(a²+b²/4). For u>r the only ground
  states are the two Néel products and the exact gap above them is 4(u−r).
  A positive-bond inequality supplies the lower bound and a closed
  one-defect sector attains it. At u=r the L+1 ground states and quadratic
  mode recover C15's nonconformal factorization line. For u≤−|a| and b≠0
  the source gap theorem gives Δ_L≥2|b|. Thus, after a=1 and b>0,
  the unresolved candidates in this charge class lie in
  −1<u<sqrt(1+b²/4), u≠0. This is not an assertion that this whole region
  is critical. Zhang–Song, arXiv:2609.29728v1, provides a matching exact
  line and gapped side; its broader phase diagram uses DMRG.
- **Exact limits of the rotation route [proved].**
  Arbitrary product rotations followed by standard Jordan–Wigner do not
  enlarge the known quadratic locus. Product rotations can turn a uniform
  zero-field bond into homogeneous XYZ in the bulk if and only if
  JJᵀ=JᵀJ, which is already the XYZ/axial-U(1) union. These statements
  do not classify entangling circuits or nonlocal dualities.

Saved NumPy checks cover 15 gap-theorem spectra at L=4,6,8; nine exact
staggered-family gaps and nine boundary degeneracies; 14 charge-kernel cases
and rotated copies; and independent full-matrix commutators. The largest
charge-basis commutator residual is below 1.2×10⁻¹³. Critical Ising matrices
also match their exact finite-size energies. Two workers independently
derived and cross-reviewed the charge classification; the coordinator
reran the saved checks. No human review is claimed.

### 9.18 Local digit components and a quantitative Gaussian certificate, 2026-10-28

At the user's request to advance the strongest counterexample route,
Codex and three workers added
[digit_compilation](contrib/digit_compilation/README.md). The sparse
species-graph first-order route remains the focus; no new counterexample,
final violating dimension, or microscopic central-charge measurement is
claimed.

- **Parameter reader [explicit finite construction; checked].** A literal
  96-entry collision permutation, eight controls and twelve neutral symbols
  implement the parameter-oracle walk. Its standalone head/record dimension
  is 204 D(q). It returns the parameter bits and restores all quantum
  records after 2m sweeps, then autonomously resets in another 2m. Full-space
  involution/charge checks include malformed H/T records, 289 exhaustive
  short payloads and 222 larger cycles. Actual primitive events equal the
  previous 10m-4 extraction word and its inverse. An independent worker
  checked another 504 cases against a separately derived schedule.
- **Arithmetic [explicit finite local rules; checked].** A frozen 4,428-rule
  relation implements eleven aligned binary operations, with inputs and
  scratch restored after 2b+3 transitions at width b. It uses 900 cell
  symbols, independent of b; that is not a completed sweeping-controller
  unitary dimension. Tests exhaust 60,060 signed cases, add 1,056 cases
  through 4,096 bits, and check inverse intermediate states. Memory,
  operand placement and the unitary sweep lift remain separate tasks.
- **Real generator integration [checked finite execution].** Eight opcodes
  execute through the local rules inside the frozen generator. Short
  workloads reproduce p=2,3 incidence graphs and the complete one-edge
  circuit. The full p=2 run executes 1,175,426 local arithmetic calls at
  128 bits (304,435,334 local transitions), matches every word-machine
  stack/instruction state and all 287,798 abstract circuit gates, and checks
  the final complete state. The larger run retains no whole-VM inverse
  journal; the short cases check exact journal reversal. This is hybrid
  execution, not a fully local compiled generator.
- **Gaussian certificate [proved matrix inequality; conditional infrared
  application].** For delta=||G_ref^(-1/2)(G-G_ref)G_ref^(-1/2)||,
  the complete vertex dimension form is bounded below by
  sqrt((1-delta)/(1+delta)) times its reference form. Hence delta<5/13
  keeps every reference-dimension-three-or-higher vertex irrelevant.
  The proof includes arbitrary density-current blocks and noncommuting
  velocities, and is sharp already for one mode. Independent numerical
  tests cover 53 modes, velocity ratios through 10,000 and large-cell
  coordinates; exact resonance checks extend through N=186. Actual
  dressed fillings and nonresonance margins must also be controlled.

The Gaussian bound is now linked from the candidate's infrared audit.
The all-order local connected-excursion estimate with mobile cell
projectors remains a stated physics input; this continuation does not
upgrade it to a theorem. All four new check suites and the larger circuit
integration passed coordinator reruns. Code, finite rule tables, proof,
resource scope and reproduction reports are retained in the contribution.
No human review is claimed.

### 9.19 Reversible single-head digit compilation, 2026-10-28

Continuing the finite compiler work, Codex and three parallel workers added
[digit_compiler](contrib/digit_compiler/README.md). This supplies a complete
spatial reversible compilation layer for a supplied finite source table;
the whole high-channel generator is still not such a joined table. No new
counterexample, final violating dimension or central-charge estimate is
claimed. The parameter-reader and Gaussian work in §9.18 is complementary.

- **Finite arithmetic and serial memory [implemented and checked].**
  A fixed eight-tape, six-symbol ALU has 4,749 active controls and 28,494
  digit transitions. It implements all signed arithmetic operations needed
  by the frozen generator, including multiplication, binary long division,
  modulo and variable shifts. A separate five-tape machine implements
  GET, variable-width functional SET, APPEND and LEN with 97 controls and
  207 rules. A 20-rule copy machine preserves results before exact cleanup.
  These are actual current-cell transition tables, with finite alphabets
  independent of operand width. Canonical input and blank-scratch framing
  are hypotheses, not yet a complete Hamiltonian format theorem.
- **Reversible spatial compilation [proved and independently checked].**
  Interleaved tagged cells simulate the digit machines with one physical
  head and 98 arithmetic or 62 memory symbols. An independent audit gives
  an explicit noninjectivity witness for this deterministic intermediate.
  A subsequent construction repairs it: source rule labels are carried
  to a right-hand history frontier, appended, and the head returns to its
  marked work cell. Fixed incoming directions and unique target/written
  pairs prove partial reversibility on the complete declared domain.
  Its time and space overhead are explicitly polynomial. Actual lowered
  ALU computations 0==1 and 1+1 take 648,266 and 3,055,626 reversible
  physical microsteps; complete inverse execution restores all work and
  history without a host journal.
- **Full-space gates and autonomous fixtures [proved gates; checked cycles].**
  Stationary transitions become reversible two-move excursions. Fresh and
  updated head sectors complete collisions to onsite involutions, and a
  two-site involution moves the head. All malformed collision inputs and
  zero/two-head bond inputs have specified action. Independent exhaustive
  checks cover 407 collision inputs and 165,649 movement pairs. Repeating
  the same gate layers gives 63 complete erasure/reset cycles, with 14,202
  sweeps, correct intermediate output and exact restoration. No host
  direction change drives these cycles. They are finite initialized
  fixtures, not a whole-generator readout or illegal-sector proof.
- **Real-generator interfaces [checked hybrid execution].**
  All 25,234 arithmetic events in the p=2,3,5 graph, width and complete
  one-edge circuit runs agree with finite-table arithmetic; 1,712 distinct
  operand pairs are replayed. Memory checks cover all 8,165 events in the
  p=2,3 graphs and 906-gate one-edge circuit, using 4,436 distinct replays
  and 75,446,718 digit transitions, including every SET event. The host
  still selects heap blocks and schedules VM instructions. Source-loop
  bounds give an arithmetic-only cost below 2^80 A^8 under the existing
  generator assumptions; no word-step count is reused as a digit-history
  count.

All eleven new verification scripts pass coordinator reruns. Separate
workers implemented arithmetic, memory and the compiler audit; independent
review caught and repaired a missing declared return-control family. The
remaining finite task is the complete interpreter's stack, heap aliases,
allocation, dispatch and output table, then initialization, parameter input,
charged execution and endpoint readout in one counted core. The full
history kernel/gap and fixed-coupling infrared claims remain separate. No
human review is claimed.

### 9.19 Route A gaps, round one (2026-10-28; check scripts of every run below rerun by Claude and passing)
- **C21, a fixed-alphabet cell** (`codex/work_C21/`, `LOCAL_CURSOR.md`). A 16-state site supports an explicit nearest-neighbour cell with an exactly degenerate 2^n code, exact parity and a nonscalar source resolvent; its second-order exchange is exactly Gaussian with hopping rank 2n and between 2n + 1 and 6n accessible energies. So local repeated terms do overcome the scalar-resolvent obstruction. Limitation: its hopping matrix is J = −2RG_f (a reversal times a positive matrix), so at zero field the Fermi points sit only at 0 or π and the one-flavour chain is gapped; independent flavour fields and phases are not yet built. Also proved: C19's original weights cannot be realized on a product code in a cell of length four or more (left and right first moments fail to commute), repaired by delayed Jacobi filters with n + s poles.
- **C22, a computation-free cell** (`codex/work_C22/`). A 17-state cell with an exact 2^n code, exact parity, vanishing first order and enough spatial poles, checked by full ED up to dimension 110,592; but its n = 2 second-order hopping matrix is provably gapped, and a uniform flux does not give independent phases. Tree buses with unitary return paths are excluded at second order.
- **C25, cell inside the mobile wrapper (G3)**. The length mismatch of the second-order bond is tunable, O((w/U)²) + O(Δ_src/Δ_wrap), but Δ_wrap includes a vacancy-tunnelling gap of order w⁸/U⁷, so equalization is not exact at finite parameters; the mixed uλ term vanishes for this cell, and the proposed interface keeps the core frustration-free.
- **C26, Sidon target and infrared transfer (G5, G6)**. The Sidon target preserves C19's exact construction, every traceless effective coefficient carries ρ² at every order, and the fourth-order norms over ρ² are 1.24 and 4.38 (one bond) and 0.86 and 1.61 (adjacent bonds) at n = 1, 2, so the second-order infrared transfer stays viable.
- **C27, gluing data (G7)**. Done for C8's projected identity-transport wrapper: the common Jordan–Wigner boundary sign, the translation character and the winding lattice; a physical harmonic carries momentum ρ(q_pay + 2πD) and its length-field dressing recovers C20's resonance and dimension test (finite-ring statements proved, full-spectrum checks to dimension 512 at 2.7 × 10⁻¹⁴).
- **C8 run 8**. An explicit autonomous clean-readout core of dimension 446400(q + 1) with q²/8 logical fermions, exact U(1), both Clifford ports and coloured density readouts (a first-order, route-B-type candidate; finite-coupling infrared control still needed).
- **Net for route A.** Fixed-alphabet cells with exact codes, parity and growing hopping rank now exist (C21 at 16 states, C22 at 17). The single missing ingredient is independent flavour fields and hopping phases (or an exact U(1)) at fixed alphabet, so that the Fermi points are not pinned to 0 and π. C28 is working on exactly this.
- **d = 2 multicritical hunt** (`contrib/d2_multicritical/`, Claude agent). No c > 1 anywhere: where Ising lines meet c = 1 regions in XYZ in a field (including the transverse-field XXZ and dual ANNNI cases), the junction is gapped or non-conformal (VUMPS finds a gapped strip of ξ ≈ 7 between them), adding a Dzyaloshinskii–Moriya term moves the Ising line only to a Lifshitz point, and XXZ in a tilted field stays at c_S ≤ 1.21. Calibrated on two c = 3/2 controls (1.506, 1.42) and two XX chains (2.008).
- **d = 2 two-fluid hunt** (`contrib/d2_two_fluid/`, Claude agent). No c = 2 phase. With a U(1) the only bonds are twisted XXZ, the pair-band minimum always sits at twice the single-particle minimum (so the conversion vertex never oscillates), and binding needs V < −2t, where the fully polarized state already undercuts the vacuum; this matches C15 and C17. Without U(1), terms that create a particle next to an existing one unpin the pair minimum (57 of 10,100 random bonds pass every two-body filter), but the same terms pull a dense canted state below the vacuum before the second Fermi sea fills (0 of 270 product-state checks, 39 of 53 survivors jump to a dense state at L = 12). On the best candidate, finite-entanglement scaling gives c = 1.01 to 1.03 in the surviving window, against 2.02 and 2.14 on two c = 2 controls. Two candidates are unresolved at small bond dimension.

### 9.20 A fixed-alphabet route-A interface with exact fermion number, and the referee report on the species-graph candidate (2026-10-28)
- **Route A, the port interface (G1)** (`contrib/routeA_compile/`, Claude agent; the switch-cell ED rerun by the coordinator: exact 2^n code, zero number-changing terms, zero first order). Ports are fermionic Bell dimers whose excitations are mobile carriers that trade a fermion with data sites only when present. The cell stays frustration-free with an exactly degenerate 2^n code and vanishing first order, and fermion number and parity are exact at every order, so the pairing that gapped the earlier XX target cannot be generated. C18's obstruction does not apply because the port excitation moves into the cell rather than sitting under a static check. The second-order bond is Gaussian at leading order, with non-Gaussian parts falling as τ² (0.215, 0.055, 0.014 at τ = 0.3, 0.15, 0.075). The plain cell fails as C21 and C22 did (the two ports see the data in opposite orders). A same-schedule lemma shows that if the junction hopping is Hermitian and definite, all n bands are partially filled. A central-switch cell P_L L^m M R^m P_R meets it: the junction matrix is symmetric to 10⁻¹⁶ and definite, and the effective chain has c = n at n = 1 to 4 once τ is small enough; at n = 2 a 5% arm asymmetry plus a carrier phase moves the Fermi seas off 0 and π, giving pair tunnelling a nonzero momentum. Local dimension 48 for the frozen pattern, about 120 with C7 walker tracks (heuristic). Open: a selector that forces equal arm lengths, the scales (λ_min/λ_max drops by about 10⁻³ per flavour), displaced Fermi seas beyond n = 2, and the infrared fate of the τ² and λ³ non-Gaussian terms.
- **Referee report on the species-graph candidate** (`contrib/review_species_graph/`, Claude agent, independent of C24 and reaching the same two serious gaps). Reimplemented checks pass: the graph rule by exhaustive enumeration at p = 2 and 3 (exactly 56 six-field conversions at p = 2, two per hexagon), the two-port isometry to 10⁻¹⁵ with exact species covariance, the circuit components, and the repaired input format accepting exactly 2^N words. Serious: the autonomous controller does not yet exist as a whole, and the step from the first-order free model to c = N + 1 at fixed coupling rests on infrared bounds that are stated, not derived; no microscopic instance has been simulated. Fixable: the tolerance should read N·g·log(1/δ) ≪ 1, and six-field conversions are irrelevant only through the stiffness bound κ > 2/3. Counting only components already explicit gives C ≥ 2.2 × 10⁷, so the first violation would be at p ≳ 4 × 10¹¹ and d ≈ 10³⁵ or larger. One payload flavour plus the length mode stays at c ≈ 2 up to interaction 0.5 by two estimators.
- **Route A, an explicit 23-letter computation-free chain** (`contrib/routeA_spatial/`, Claude agent; the coordinator reran the band check and the two-cell ED at n = 1, whose error over λ⁴ is constant at 0.027). Extending C22's 17-letter walker cell with six express-track letters makes both ports see the flavours in the same order, and a weak loop link separates the Fermi-sea centres. Proved: an exactly degenerate 2^n code, exact U(1) charge and parity, zero first order, no three-cell term at second order, and the C7-type length selector inherited from C22, so the whole bond, selector included, has d₀ = 23. Numerical: the Hermitian part of the second-order hopping is positive definite at n = 2, 3, 4 (n = 4 needs a larger express energy), every band crosses a small common field with distinct centres (n = 2: centres ±0.0017, margins δ₂ = 0.54, δ₄ = 0.0068), and full ED of two and three cells up to dimension 373,248 matches the exact second-order bond. Heuristic: at each fixed n, with n-dependent small couplings, the logical chain sits in the stable incommensurate region with pairing forbidden by U(1). Limits: the narrowest band shrinks by a factor of about 50 to 100 per flavour (the band scale at n = 3 is about 10⁻¹⁰), the cells are frozen and their payload grows with length, so the equal-payload mobile pair of C20 is not yet available, and the infrared step is taken from `contrib/ir_control/`.
- **C28, fields and phases for the bus cell** (`codex/work_C28/`, Astra; checks and the algebra verification rerun by Claude). A frustration-free, repeated 17-letter initialized core with an exact 2^n code, exact total U(1) and a nonscalar source resolvent (the vacuum walk replaced by a Laplacian, so C8's wrapper applies). A uniform density gives a common chemical potential and endpoint densities give two independent fields at n = 2; doping the bus repairs the zero-field band gap, giving four simple Fermi points with two distinct Fermi momenta, and the generated inter-flavour pair tunnelling is on the irrelevant side of the one-loop flow (pairing is forbidden by U(1) to all orders). Two- and three-cell ED up to dimension 110,592 confirms the fields and the full second-order bond with λ⁴ remainders. Obstruction of restricted scope: diagonal two-site readouts on this code give at most four independent field profiles, so they cannot set n independent fields for n ≥ 5. The n = 2 logical chain is saved for numerics.
- **C29, a two-bus computation-free cell** (`codex/work_C29/`, Astra; checks and the rational certificate rerun). A 26-letter repeated bond with two spectrally distinct virtual buses gives normal, diagonal exchange with independently adjustable phases and chemical potentials at n = 2 (at leading order λ²g²) and four distinct displaced Fermi points; arbitrary position-dependent flux alone cannot do this in C22's single-bus family. The same rules are provably non-normal at n = 3, so the arbitrary-n mechanism remains open.
- **Net for G1 after this round.** At n = 2 there are now several explicit fixed-alphabet cells with exact U(1) and a good two-flavour logical chain (C28, C29, `contrib/routeA_compile/`, `contrib/routeA_spatial/`). For general n, `contrib/routeA_spatial/` (23 letters) has all bands crossing with distinct centres through n = 3 and `contrib/routeA_compile/` reaches c = n through n = 4, while the simplest extensions of C28 and C29 fail at n ≥ 5 and n ≥ 3. The open question is a mechanism that provably works for every n; C31 and C32 are on it.
- **Route A with the incommensurate target at n = 2** (`contrib/routeA_incommensurate/`, Claude agent; the coordinator reran the cell construction, reproducing J₀ = 0.0222). A D = 64 cell (two qubit ports Bell-paired with partners, plus ℂ⁴) with an exact on-site Z₂ equal to logical parity; C9's recipe as written breaks that grading, so the hopping is built from excited states that mix odd and even syndromes. The second-order bond equals the two-flavour incommensurate target times J₀ = 0.022 up to a constant, to 10⁻¹⁵. Two- and three-cell ED gives a residual scaling as λ³, a three-cell term of about 1.5λ²J, and generated same-flavour pairing of 0.34λ and 0.42λ (in units of J = J₀λ²); the Fermi points move by less than 0.06 up to λ = 0.2 and stay at least 0.54 from every resonance, so none of these terms is resonant. On three-cell rings the logical chain reproduces the full chain's 64 lowest levels to 10⁻⁴J at λ = 0.05. Logical-chain rings up to L = 9 track the free target at λ = 0.1 and 0.2, while the symmetric contrast develops size-independent splittings (masses). A calibrated long-distance c was not obtained: VUMPS fails even on the free target at these fillings, and full-chain VUMPS did not converge. Open-chain and four-cell runs are pending on Typhon.
- **C30, a selector for symmetric switch cells** (`codex/work_C30/`, Astra; checks rerun by Claude, all assertions pass). Symmetric frozen cells P_L L^m M R^m P_R are selected by one chemical potential at d₀ = 120, with explicit margins, checked by walker spectra through arm length 8 and complete selector rings through 59,049 states. At the crossing of adjacent lengths asymmetric cells also appear, so the mobile version keeps a symmetric occupied core and varies only neutral padding, giving d₀ = 8112 (2028 without endpoint compensators). The all-n carrier band structure and the finite-coupling infrared step remain open.
- **C8 run 9** (`codex/work_C8/`, Astra). A complete periodic-bond implementation, a simpler uniform-vacancy wrapper, and rational Fermi points with uniformly positive velocities in place of the exponentially slow target; finite checks through legal-sector dimension 2,903. Remaining: connected local elimination of the massive cell excitations and the infrared identification.
- **Route A works end to end in one explicit Hamiltonian** (`contrib/routeA_demonstrator/`, Claude agent; the coordinator reran the selector margins, the second-order bond and its λ³ residual). A single ten-letter nearest-neighbour bond (letters P, a port qubit, and B, eight states holding two port partners and one logical qubit) with five terms: K n_B n_B, −μ n_B, a port check on (P, B) and (B, P) pairs that annihilates the Bell state Φ⁺ with a κ term trading a charged syndrome for Φ⁻ and moving the charge onto the logical qubit, and λ(XX + YY − ZZ) on (P, P) pairs. It has an exact U(1), charge conjugation and reflection. Proved: site types are conserved and the period-3 cell pattern wins with a margin of at least 0.021 per site; each cell has an exactly degenerate two-dimensional code with gap 0.4 and every port compression vanishes; the second-order bond is exactly a logical XXZ chain with J_xy = J⊥J_zκ²/[4D(D² − κ²)], the two-cell residual scaling as λ³; at κ = 0 the second-order bond is scalar (the C18 obstruction), which the κ term escapes. Numerical, on the full microscopic chain: finite-entanglement scaling gives c = 1.01 ± 0.04, 1.01 ± 0.03 and 1.00 ± 0.03 at λ = 0.05, 0.10 and 0.15, with correlation lengths up to 30 cells, against 0.99 for XX and for the effective chain; microscopic rings up to 18 sites match the effective chain's gaps to 0.1 to 0.6% and its Casimir estimate to 3 × 10⁻⁴; z = 1. This is the first end-to-end microscopic confirmation that spontaneously formed cells coupled only at second order through idle ports give a conformal chain. Limits: frozen cells, one flavour, λ below about 0.23.
- **C31, displaced seas and infrared of the switch cell** (`codex/work_C31/`, Astra; `verify.py` rerun by Claude, matching the full computation to 10⁻¹⁵). The switch cell's leading hopping is definite for every n, with an explicit positive lower bound (proved), and non-Gaussian corrections start at λ²τ⁴ and λ³τ². At n = 2, 3 and 4 two routes give an open weak-coupling c = n region: arm asymmetry, a neutral-carrier phase and a uniform logical chemical potential remove every non-density quartic resonance; or, with centred seas, the interactions the cell actually generates put pair transfer on the irrelevant side of its one-loop flow [numerical with heuristic RG]. A correction: at zero field, staggered particle-hole symmetry leaves umklapp resonances even with distinct centres, and a small uniform chemical potential is needed. The all-n separation of centres and the all-n sign of the relative stiffness remain unproved.
- **The explicit candidate, written out** (`contrib/explicit_candidate/README.md`, Claude agent; the coordinator reran its three scripts). One fixed bond on ℂ²³ for every n (the selector, the cell rules K1 to K7, a penalty, the port coupling λ and an on-site operator, with a table of coefficients and a schedule of couplings in n) whose low-energy theory is c = n on frozen cells of length n + 1, so n = 23 would beat c ≤ 22. Re-coded from the text alone, the bond matches the source code entry by entry for n = 1 to 4 with zero leakage on the complete 12,167- and 279,841-state spaces. New findings: the hopping matrix is proportional to sin(nφ), so the express phase must scale as π/(2n) (the value 0.4 nearly kills the hopping at n = 8); positivity of its Hermitian part needs the express energy D_x ≥ 10 at n = 4 and ≥ 20 at n = 5, 6 (now checked to n = 6); the cell gap falls like n⁻³ at zero detuning, so the schedule needs Δ ≥ |g|; and the n = 3 infrared conditions need g ≈ 10⁻³. Ranked gaps: distinct centres and margins for every n up to 23, a measured c at n = 2, explicit infrared-transfer numbers, the selector and illegal-sector bounds for 23 letters, the frozen-cell scope, and review.
- **C32, mobile version and all-n crossing** (`codex/work_C32/RESULT.md`, `PROOFS.md`, Astra). A mobile construction at d₀ = 2736 with controlled second-order equalization, and finite all-n bounds proving that every band crosses (two- and three-cell checks at n = 1, 2 with λ⁴ errors). Still open: distinct centres and Sidon margins for every n. C33, still running, has a 256-bit certificate that the leading hopping's Hermitian part is positive at n = 23 at one parameter point.
- **C34, d = 2 round two** (`codex/work_C34/`, Astra, checks not yet rerun). No full physics-level proof. The open strata keep their maximal dimensions (5, 5, 8) as search regions, and the actual unresolved critical subsets cannot yet be assigned dimensions; the missing step is a dynamical classification that controls neutral sectors, isolated interacting fixed points and singular neighbourhoods, which the numerical hunts do not supply.
- **C35, the log-Casimir bound** (`codex/work_C35/`, Astra, checks not yet rerun). A ≤ π‖h‖ log₂ d survives every solvable family checked (small-ring Sp fits were inconsistent and rejected). Entropy and current sum rules leave one susceptibility bound missing for a proof. If the bound holds, c² ≤ 3R(log d)²/[(1 − ε) log 2], and the route-A hierarchies respect it through their suppressed velocities, including the mobile length mode.
- **C36, infrared bounds of the species-graph candidate** (`codex/work_C36/`, Astra, checks not yet rerun). Derives the irrelevance of six-field conversions for stiffness κ > 2/3 (enumeration checks at p = 2, 3, 5, 7) and gives the physics-standard route to c = N + 1 with the corrected tolerance N g log(1/δ) ≪ 1 at fixed p and fixed coupling, conditional on the complete core and a local low-frequency resummation that controls the retained background.
- **C8 run 10** (`codex/work_C8/`, Astra). Uniform local dressing, ground-energy and defect bounds, length-tail estimates and low-frequency source bounds for the C8 chain, with five check scripts; two soluble controls show these do not by themselves establish a CFT at fixed coupling. The remaining gap is a connected local effective action followed by infrared identification.
- **Red team of the 23-letter chain** (`contrib/redteam_routeA_cells/`, Claude agent, independent reimplementation as one 529 × 529 bond). Every decisive finite claim reproduces: d₀ = 23 with the selector included (the full 23⁴ ring at n = 1, illegal words included, selects b = 2), the exact 2^n code, U(1), parity, zero first order, no three-cell λ² term (to dimension 373,248), T and the positivity of its Hermitian part at n = 2, 3, and fourth-order terms suppressed like T. Ranked problems: (1) frozen cells, fatal for the conjecture as Latorre and Sierra state it (rings with b ∤ L are degenerate and non-conformal, and spontaneous breaking on compatible rings is akin to the blocking d → d^k they set aside), so the mobile version is required; (2) the red team did not find a mobile version (C32's, with a positive core and two cell lengths carrying the same payload, needs to be checked against its objections); (3) n ≥ 23 is extrapolation: positivity holds through n = 10 with an express detuning growing to about 30 to 50, and margins scale like 30⁻ⁿ; (4) the stated working point breaks at n = 8 (an excited-token state drops below the code) unless Δ ≳ g; (5) the n = 3 point needs g ≈ 0.003 rather than 0.05 to pass its weak-coupling test; (6) pair transfer sets the margin and needs the loop link; (7) the express phase must be retuned with n. No circular constraint in the ordering of scales at fixed n. The switch cell is frozen, has no selector, and its symmetric arms pin the Fermi-sea centres.
- **C33 and C37: the frozen 23-letter chain certified through n = 23** (`codex/work_C33/`, `codex/work_C37/`, Astra; the coordinator reran C33's positivity certificate at n = 8, C37's crossing bound for n = 1 to 23 and its pair certificates at n = 3, 8, 12). With the unchanged 23-letter alphabet and an explicit parameter schedule, every n from 1 to 23 has distinct band centres and excludes every nontrivial number-conserving quartic resonance, pair transfer and umklapp included, with explicit (very small) margins and finite-g correction bounds; at n = 23 C33 certifies 46 simple Fermi points with δ₂ > 0.014 and δ₄ > 3.2 × 10⁻⁷ at g = 2⁻³⁶⁵, and common-field crossings are proved for every n. The finite-dimensional input for the frozen chain at d = 23 is therefore complete. Open: an all-n theorem (needed for unbounded c and for the mobile chain, which must reach n of order its local dimension) and the infrared step.
- **C38, all-n margins from the hierarchy** (`codex/work_C38/`). The scale hierarchy alone is insufficient (cross-scale mixing changes the narrow coefficient by a fraction 2/5 even as the energy separation diverges), but an explicit all-n lacunary endpoint target has quantitative resonance margins covering every quadratic cross term, checked exactly through n = 128. Realizing those inequalities in the actual frozen or mobile core is open.
- **C40, infrared transfer** (`codex/work_C40/`). The velocity hierarchy imposes stricter coupling requirements but no new obstruction at fixed n; connected bounds at every finite perturbative order and the RG through all band scales are given. A careful physicist can accept c = n (frozen) or n + 1 (mobile) at fixed n, conditional on an exact connected low-frequency remainder bound and on non-resonant witnesses for the actual core.
- **c at n = 3 and 4 on the logical chains** (`contrib/routeA_n34_chains/`, Claude agent, not rerun by the coordinator). Switch cell: all 2n Fermi points survive self-consistent Hartree–Fock with every generated term and the exact third-order bond, the dimensionless couplings at the Fermi points are at most 1.1 × 10⁻³ (n = 3) and 1.3 × 10⁻⁴ (n = 4), and a rescaled proxy chain with the same Fermi points gives c = 3.00 and 4.00 from the Casimir energy, with entropies matching the c = n calibrations to 0.004. 23-letter chain at n = 3: c = 3 at g = 0.05 (estimators 3.41 and 2.98), unresolved at g = 0.1, and c = 2 at g = 0.2, where Hartree–Fock empties one band, consistent with the red team's requirement of small g. At n ≥ 3 the phase must shrink (10⁻⁵ at n = 3, 10⁻⁷ at n = 4) and a small field is needed, without which particle-hole symmetry leaves a mixed umklapp resonant and pair transfer is exactly resonant. The true chains reach their conformal regime only beyond about 2000 cells, so c was measured on the proxy chains; n = 4 of the 23-letter chain needs extended precision.
- **C39, all-n margins in a 46-letter core** (`codex/work_C39/RESULT.md`, `PROOFS.md`, Astra, unrefereed; the coordinator reran `checks.py` and the exact rational certificate at n = 2). A fixed-alphabet core of 46 letters with two independently phased buses and a biased frustration-free vacuum walker. Its proof gives, for every n, distinct band centres and positive pair-sum margins through a terminating parameter search with consecutive centre and slope ratios below 1/32 (lacunary), with exact rational certificates (Fraction arithmetic, Sturm root isolation) at n = 1, 2, 3 and full second-order bonds checked by two- and three-cell ED with fourth-order remainders. With the inherited mobile wrapper, d₀ = 4464. If it survives review, this is the all-n finite-dimensional input that unbounded c at fixed d₀ needs; the thermodynamic elimination of massive sectors remains the stated assumption.
- **C8 run 11** (`codex/work_C8/` §55 to §58). Factorial length tails, control of the ground-state marker density, and exclusion of illegal selector sectors after strengthening one penalty, checked through dimension 7,007; the controlled local infrared reduction remains.

### 9.21 One-vacancy d=24 counterexample at the physics standard (2026-10-28)

[proved local ingredients; controlled physics inference; independently
reviewed by AI] The
[counterexample campaign](contrib/counterexample_campaign/COUNTEREXAMPLE_D24.md)
adds one vacancy state to the explicit 23-letter bond. Species-independent
vacancy hopping preserves the ordered occupied labels, and the original
bond acts only on adjacent occupied pairs. This is a literal Hermitian,
one-site invariant nearest-neighbour bond on d=24 states. At half occupied
density its weak-coupling infrared theory contains 23 payload compact
bosons and one charge boson, giving c=24>d-1=23. The conclusion uses
conventional Gaussian RG stability; it is not a mathematical
lattice-to-CFT theorem or a microscopic numerical c=24 determination.

[proved and checked] The campaign independently replayed C33's
n=23 matrix/inertia certificate, completed all-length and illegal-sector
exclusion for compatible material rings at fixed nonzero port coupling,
and checked squeezing, cyclic twists and the exact charge-sea projection.
The 576-by-576 bond fixture was independently reconstructed with zero
entry error. Exact momentum optimization through L=256 checks a positive
free-charge threshold for restoring nonforward payload momentum, even
with arbitrarily many particle-hole excitations and bounded number
rounding. On arbitrary physical rings, rounding occupied number to a
multiple of 24 costs O(1/L) and avoids fixed-energy material defects.

[derived continuum analysis; numerical support] Material gap coordinates
turn vacancy motion into the local q=0 phase model. At half physical
occupancy K_gap=4. Physical translation makes every nonzero charge
umklapp dressing cost at least four in scaling dimension, including the
fractional marker harmonics. The other non-Gaussian resonances have
dimension at least three. Local defect dressing is O(J^2/t) per defect,
with the canonical zero-mode term explicitly separated; tagged-defect
checks reach L=20. Choose g, then lambda, then J/t sufficiently small,
all fixed before the thermodynamic limit and compared with the smallest
payload speed and momentum margins.

[read; proved application] Bravyi-DiVincenzo-Loss's fixed-order local
Schrieffer-Wolff expansion supplies an exact local-unitary normal form
with a full-space interaction-strength remainder O(lambda^(m+1)). Its
block-diagonal stability lemma, applied after subtracting half the cell
defect penalty, gives a uniform massive auxiliary threshold above the
reference ground. These statements avoid an uncontrolled extensive-band
projection and make C40's stronger all-order exact kernel unnecessary
for the conventional phase-stability inference. Two independent
primary-source audits are linked in the synthesis.

Four Codex GPT-6-family workers (coordinator plus three repeatedly
redeployed agents) produced and cross-checked this campaign. All three
reviewers accept the final argument at the repository's physics standard;
no human review is claimed. Fresh collaborator C41 progress independently
agrees with the harmonic and number-rounding analysis, but its own final
matching assessment remains pending. The new C39 all-n core and C42 work
concern the stronger fixed-d unbounded-c program and are separate from
this finite-d violation.

- **c at n = 2 on four explicit cells** (`contrib/routeA_n2_chains/`, Claude agent, not rerun by the coordinator). No cell flows away through interactions: at every bond dimension each interacting chain matches its own Hartree–Fock free chain (entropies to 2 × 10⁻³, 0.02 for the 23-letter chain), with interactions at the Fermi points of 10⁻⁴ to 0.06. The flows observed come from generated quadratic terms: C28 at the source's fixed fields fills its narrow band for g ≳ 0.15, and C29 at its source fields sits 0.03 × 2π from the ν = 1 umklapp and becomes insulating at g = 0.5; the switch cell's exact particle–hole symmetry puts it on that umklapp at zero field. Retuning one uniform field restores c = 2: finite-entanglement scaling gives 1.89 to 2.06 for C28 and C29 (calibrations: XX 0.99, free chains 2.02 and 1.98), with no translation breaking. The second estimator (ring Casimir) failed its own calibration at the reachable sizes, so these are single-estimator results; for the 23-letter and switch cells, c = 2 rests on the free reference plus the resonance margins (the 23-letter chain's closest resonance, pair transfer, is detuned by 0.011).
- **C41, independent theory of the d = 24 vacancy chain** (`codex/work_C41/`, Astra, checks not yet rerun). Agrees with the collaborator's harmonic classification: including fractional marker harmonics, an open density window around half filling has no unavoidable relevant mixed vertex, so c = 24 is a supported prediction. Qualifications it keeps: local matching at finite J/t and positive dressed defect costs are assumed, and exact uniqueness on every large ring is not established, so it is a conditional rather than a completed counterexample.
- **C42, all-n finite-band result** (`codex/work_C42/`, Astra). Proves, in C39's 46-letter core, the all-n band statement with a terminating rational parameter schedule and explicit finite-D, finite-g error bounds, with two- and three-cell ED at n = 1, 2, 3; a conservative vacancy completion has fixed dimension 4465. The 23-letter all-n problem and the interacting mobile infrared transfer remain open.
- **C43, vacancy completion of the all-n core** (`codex/work_C43/`, Astra). Full microscopic ED at n = 1, 2 (103,823 states each) finds a sixfold ground-state degeneracy at those sizes, and larger invariant sectors confirm the squeezing, twists and first-order projection; uniqueness on every large ring is not established, so the fallback remains the wrapper version at d₀ = 4464.
- **Referee of C39's all-n proof** (`contrib/referee_C39/`, Claude agent). No fatal error; the finite all-n band claim holds. Checked by hand: the second-order bond formula and normalization, the Newton-polygon lemma, the slope formula and the pair-sum classification, which with exact U(1) is complete (it covers pair transfer, exchange and every umklapp; the smallest margin is the slowest flavour's 4k_F umklapp). The bond rebuilt from the rule table matches the stored bonds to 7 × 10⁻¹⁷ at n = 1, 2, 3, and margins stay positive in 900- to 1300-digit arithmetic at n = 4 to 8. Serious, as C39 itself flags: the mobile wrapper is inherited from C8 and C32 without audit (its hypotheses do hold for this core), and the thermodynamic elimination is assumed. Fixable: the parameter search lacks an explicit depth bound, and the cluster ED residual exceeds the logical hopping, so the ED does not test the bands. The d₀ = 4464 count is correct. The one-vacancy completion does not apply as it stands, because this core is frustration-free at every cell length and a selector would have to be added.

### 9.22 Independent vacancy continuation and microscopic control (2026-10-28)

The [independent continuation](contrib/vacancy_continuation/README.md)
adds four complementary results, using a coordinator and three parallel
GPT-6-family Codex workers. No human review is claimed. Existing
collaborators' construction files were left unchanged.

- **Density without exact tuning [proved; independently checked].** For
  an arbitrary internal Hermitian bond h with one vacancy added, zero
  occupied-particle chemical potential and J||h||/t<=1/1000, every
  thermodynamic ground-state density lies in (19/40,21/40). The exact
  finite-L bound is |N/L-1/2|^2<=3pi J||h||/(16t)+pi^2/(32L^2).
  This entire window is below 1/sqrt(3); the existing physical-harmonic
  analysis therefore retains its dimension-three stability budget.
  The inherited C33 reference gives physical mismatch >6.24585e-9
  throughout the window. This removes exact half-filling tuning under
  the same matching/RG assumptions; the density estimate itself is
  independent of those assumptions and does not prove criticality.
- **Dynamical charge response [derived continuum bound].** The Gaussian
  gap-density kernel is bounded by 1/e_gap'' at every momentum and
  Euclidean frequency, equal to 4/(pi t) at half filling. Smooth
  operator-valued internal sources therefore have an O(J^2/t) influence
  coefficient. This strengthens the matching argument without replacing
  it by a lattice operator-norm theorem. A fresh primary-source review
  accepts the finite-order local SW construction and auxiliary threshold.
- **All-n finite input [checked analytical argument; exact checks].**
  An independent positive-matrix-pencil derivation validates C39's
  lacunary eigenvalues, field slopes, repeated-label pair comparisons
  and parameter order. A new exact n=4 reference certificate has eight
  simple roots and 630 distinct-pair comparisons, independently bounded
  by delta4>1.0437e-54. The Gram evaluator is inherited; a separate
  determinant and rational Cayley verifier checks the resulting data.
  This does not implement a selector/mobile completion. C43's final
  report, received during synchronization, proves that the naive d=47
  vacancy lift lacks a selector and has exponential ground degeneracy;
  its fallback is the inherited d=4464 wrapper. The new finite-band
  audit does not repair that separate obstruction.
- **Microscopic mechanism test [numerical].** A literal d=3 vacancy/XX
  chain on L=8,12,16, at J/t=.2,.1,.05, gives half-chain size slopes
  c=2.026 to 2.038, or 2.011 to 2.019 on the two larger rings. Independent
  flux and number-sector responses resolve the charge and spin velocities;
  their sum agrees with the Casimir coefficient to within 0.2 percent
  in the three-size fit. That fit has three coefficients for three data
  points, so the agreement is a consistency diagnostic. Single-size
  entropy shapes still drift, and no precision thermodynamic c is
  claimed. An exact free-charge calibration validates the physical
  partial trace. The largest response sector has 2,018,016 states.
  This d=3,c=2 control saturates the conjectured bound and is not the
  proposed d=24 counterexample.

The new density proof, diagnostics, numerical scripts and data are
reproducible from the linked folder. The main d=24 conclusion remains a
controlled theoretical-physics inference with explicit continuum inputs,
not a rigorous lattice-to-CFT theorem or a measured microscopic c=24.
- **Margins audit, n = 1 to 23** (`contrib/margins_audit/`, Claude agent, high-precision numerics at 115 to 313 digits, not a certificate). At one explicit schedule inside the 23 letters, every n has exactly 2n simple Fermi points, distinct centres, and every quartic resonance (umklapp included) missed by a stated margin: δ₂ ≥ 0.014, umklapp ≥ 3.8 × 10⁻³, δ₄ and pair transfer ≥ 4.2 × 10⁻⁸ (the worst cases, n = 21 and 15, rise above 3.6 × 10⁻⁷ when the field power moves by one). The minimum velocity falls by about 8.7 decades per flavour, to about 10⁻²⁰⁶ at n = 23. It reproduces C33's n = 23 point to every printed digit, and finds C37's certificate consistent but its literal parameter point practically degenerate (margins of order 10⁻¹⁹⁰⁰ for pair transfer).
- **Independent referee of the d = 24 vacancy counterexample** (`contrib/referee_d24/`, Claude agent). No fatal error. The 576 × 576 bond rebuilt from the prose matches the collaborator's fixture entry for entry and conserves its charges; C33's n = 23 computation reproduces. The vacancy step removes the frozen-cell objection in the letter of the conjecture: the 24 crystal offsets become charge-flux sectors split by O(1/L), marker order decays as a power law, and every long ring can choose an occupied number divisible by 24 (unstated consequence: the charge sector is a charge-24 liquid, and defect-free ground states are guaranteed only on very long rings); uniqueness on every ring is not proved (a toy gives a ±P doublet, which does not change c). The infrared operator bookkeeping is correct (K_gap = 4, charge umklapp dimension 4h², other non-Gaussian vertices at least 3, consistent order of limits), but the two decisive inputs are assumed rather than derived: c = 23 for the frozen chain at finite coupling, and the matching of the gapless charge fluid. The hierarchy needs λ ≲ 10⁻³¹ and J/t ≲ 10⁻⁵¹⁰, beyond any numerical test. The couplings t, J, μ and λ are not fixed numerically, so the claim defines a family. Suggested decisive test: a calibrated c = n + 1 for the vacancy lift of this bond family at n = 1 or 2 (in progress in `contrib/vacancy_mobile_numerics/`).


### 9.23 Fresh adversarial audit of the d=24 argument (2026-10-28)

[proved finite ingredients; checked calculations; controlled physics verdict]
The [new audit](contrib/counterexample_audit_20261028/README.md) used a
coordinator and three fresh parallel workers, followed by cross-review.
It verifies the construction at the repository's physics standard. No
refutation was found. The verdict retains conventional local matching
and Gaussian RG, and does not claim an all-order microscopic per-defect
inequality, a certified final lambda/J interval, or a measured c=24.

New checks include a literal full-Fock two-cell calculation of the g-squared
coefficient at n=1,2 (hopping agreement within 5.2e-16 relative error), exact
regeneration of the full 512-bit n=23 certificate, and an independent scan
of 1,167,480 ordered nonforward quartics and 56,039,040 dangerous marker
branches. Independent phase-model number/flux curvatures recover K=4.
An operator-geometric-mean proof extends the stiffness bound to arbitrary
current mixing and unequal chiral velocities. A direct loop-scale check
explains why a tiny onsite loop term moves the very slow Fermi crossings
by the certified amount.

The audit sharpened several statements: finite rings need the full rooted
coordinate quotient; the bounded smooth response is Euclidean and does
not bound the full contact operator; interacting bilinears need only have
dimension above one half for the connected interface estimate; charge
stiffness may shift at order J while induced internal energies shift at
order J-squared/t. Long wrong cells have a separate volume selector bound,
avoiding use of a nonexistent uniform walker gap. These clarifications
are incorporated in the linked synthesis and reports.

A new exact diagnostic couples the vacancy fluid to noncommuting quantum
impurities whose gaps shrink as J. At 92,378 bath states, the single-impurity
energy shift at J=.01 agrees with its retained-gap second-order resolvent
to 1.42e-6 relative error. Adjacent impurities resolve the additional
J-cubed correction from induced pair flips. This supports the matching
mechanism and is explicitly not a simulation of the entire counterexample.
All reports, scripts and data are saved; the reviews are AI-assisted,
not human referee reports.
- **C44, wrapper audit on the all-n core** (`codex/work_C44/`, Astra). The finite wrapper audit passes for every n at d₀ = 4464, with explicit parameter budgets: equalization is approximate at finite parameters, the scalar length interaction and the mixed payload terms are retained and bounded against C39's lacunary scales, and the mixed marker-port term at order uλ vanishes. Small-cell ED confirms the code, and transport and projected-ring spectra confirm identity transport and momentum gluing. No full microscopic many-cell ED or central charge.
- **C45, infrared step for the wrapped all-n chain** (`codex/work_C45/`, Astra). Concludes yes at the theoretical-physics standard, using one explicit standard Wilsonian RG input: with C44's audit and a volume-independent local Schrieffer–Wolff remainder bound, the coupling hierarchy closes for every fixed n (the schedule includes equalization errors, mixed terms and the smallest wrapper gap), giving c = n + 1 at fixed d₀ = 4464, hence unbounded c at fixed local dimension. Exact all-order infrared kernels and a rigorous microscopic CFT theorem remain unproved. An independent referee of C44 and C45 is running (`contrib/referee_C44_C45/`).


**Correction to §§9.21 and 9.23, same-day referee cross-review:** the
verified-counterexample headline was too strong and is withdrawn. The
finite identities, certificate, normal-form and conditional RG checks
remain valid. They do not complete the finite-J, common-bulk-subtracted
sector comparison uniformly over ring lengths and defect patterns. The
current classification is a credible conditional candidate, with no
refutation of the model. See the
[infrared reassessment](contrib/counterexample_audit_20261028/infrared/REFEREE_RESPONSE.md)
and [referee comparison](contrib/counterexample_audit_20261028/vacancy/REFEREE_COMPARISON.md).
New workers are pursuing the specific missing bound and possible
countermechanisms. A human review or small-n simulation would improve
confidence, but neither substitutes logically for that comparison.
- **Literature on the infrared step** (`contrib/lit_ir_step/`, Sonnet literature agent). No published theorem covers eliminating a gapped sector at fixed coupling and concluding a gapless CFT. The Hamiltonian-simulation framework (Cubitt–Montanaro–Piddock, Kohler et al.) needs couplings growing with system size; perturbation gadgets at fixed coupling (Bravyi–DiVincenzo–Loss–Terhal) leave extensive energy errors; a fully block-diagonalizing quasi-local unitary does not exist in general (Chen–Lucas–Yin), and rigorous dressing results need a gapped infrared. The clean physical confirmation is the large-U Hubbard chain (Ogata–Shiba, Lieb–Wu, two flavours with velocity ratio growing as U/t), with SU(N) Hubbard checks to N = 4; no confirmation with many flavours or strongly hierarchical velocities was found. Known failures all involve marginal flows or symmetry-allowed relevant operators generated at higher order (for example a relevant term appearing only at third order in Schrieffer–Wolff, Prembabu–Verresen 2509.20431); none involves a robustly critical target with all non-Gaussian operators detuned. The sharpest theorem that would close the step is an imaginary-time decoupling statement (exponentially decaying retarded kernels for the logical Grassmann action, uniformly in system size and temperature), feeding an extension of Mastropietro–Porta.

### 9.24 Targeted attacks on the remaining d=24 matching step (2026-10-28)

[proved partial statements; numerical diagnostics; unresolved conclusion]
Three parallel audit workers isolated the outstanding comparison more
precisely. The [long-cell attack](contrib/counterexample_audit_20261028/matching_attack/README.md)
shows that a two-sided remainder bound cannot count one arbitrarily long
wrong cell as one bounded defect. The discrepancy is favorable to selection,
and an exact variational comparison excludes the whole-ring giant-cell
sector for sufficiently large length at each fixed J. This is neither a
refutation of the candidate nor a proof for islands and arbitrary patterns.

The [quantum-selector diagnostic](contrib/counterexample_audit_20261028/selector_numerics/README.md)
performs 264 exact diagonalizations of the zero-port microscopic restriction
at selected lengths b=2 and b=3, through L=20 and 2,217,072 quotient states.
Independent physical-word and quotient constructions agree within 1.3e-14.
No wrong sector undercuts the selected one. Noncommuting b=3 walkers give a
substantial quantum correction beyond fixed scalar-source calculations.
Scripts, complete energies, residuals, and reproduction instructions are saved.

The [uniform-matching note](contrib/counterexample_audit_20261028/uniform_matching/README.md)
derives an exact cell-gradient cancellation, an all-order Gaussian energy
inequality, and gap-free spectral-filter algebra. It states the missing
spatially weighted RG bound explicitly and explains why available
single-impurity or homogeneous sine-Gordon theorems do not supply it directly.
The corrected conditional verdict remains in force; the counterexample has
not been verified or refuted by these additions.
- **Referee of C44 and C45** (`contrib/referee_C44_C45/`, Claude agent). No fatal error. The wrapper pieces that do not need the 46-letter core were rebuilt and hold for every n (charge gaps reproduced to every quoted digit, the gap bound and its asymptotics, the occupation filter, the selector bounds out to n = 4467, the equalization inequality), and the exclusion bounds do not depend on D or the core norms, as the lacunary schedule requires. The local Schrieffer–Wolff remainder bound is correct and volume-independent (reimplemented: λ³ scaling, linear in chain length, 10⁷ below C45's constant), and the coupling schedule is not circular (every inequality holds at n = 1, 2, 3; the payload scale is near 10⁻⁶⁵⁰⁰⁰ at n = 1). Verdict: a careful theoretical physicist should accept the argument given its one stated input, a standard Wilsonian elimination of the kind used from large-U Hubbard to Heisenberg, but here protected by momentum margins of 10⁻¹⁶ to 10⁻⁴¹ rather than a symmetry and never tested on this family. c should be read from conformal towers or interval entanglement, since the Casimir coefficient has no limit at irrational density. What would settle it: an independent rebuild of C44's source formula and transport overlap (checked so far only by the authors' code at n ≤ 2), a calibrated c for the n = 1 member, and a human read of C39 to C45.
- **One-vacancy mobile completion, numerics** (`contrib/vacancy_mobile_numerics/`, Claude agent; proposal credited to the collaborator). For the ten-letter demonstrator bond plus one vacancy (d = 11), exact diagonalization matched to the literal bond shows a unique ground state on L = 5 to 11, primes included (an exact Q = ±½ doublet when the occupied number is odd), against 3- to 22-fold degeneracy for the frozen chain; the frozen triplet's partners become charge-twist levels scaling as 1/L; the low spectrum matches free charge plus J times the frozen spectrum with error ∝ J²; the charge sector has z = 1 and calibrated Casimir c = 0.998. The payload is about 360 times slower than the charge, so full-chain finite-entanglement scaling sees only the charge boson (c ≈ 1.0 to 1.1) and a direct two-estimator c = 2 was out of reach; it rests on the demonstrated separation plus the demonstrator's c_payload = 1. Short rings keep a pattern defect (a slow rotor band, ground state still unique); the picture holds for J from 0.25 to 4 and densities 0 < ρ < 1 except ρ = 3/4. For d = 24 at n = 1, the frozen two- and fivefold ground states become unique on L = 5 to 8.
- **C46, toward a rigorous infrared step** (`codex/work_C46/`, Astra). Proved: uniform thermal suppression of defects in the block-diagonal reference, explicit second-order time bounds, and the Gaussian reference's thermal coefficient with the physical gluing. A counterexample shows that a cell gap alone cannot justify the raw summed Grassmann expansion, so the decoupling theorem needs a resummed treatment. The full chain's c = n + 1 remains a conjecture at the rigorous level; the shortest path is a resummed matching of the massive sector, constructive RG for that action, and control of the thermal vacuum terms to o(T²).
- **Final assembly of route A in mobile form** (`contrib/routeA_final/README.md`, Claude agent). One refereeable document covering the claim, every bond (written out or referenced), step tables with status tags, the order of limits, ranked gaps with closing calculations, and the likely referee attacks, in two forms: Form V (the 23-letter bond plus one vacancy, d = 24, c = 24 at n = 23, the primary claim against c ≤ d − 1) and Form W (C39's core in the C8 wrapper, d₀ = 4464, the route to unbounded c). New in the folder: C37's finite-field certificate ported to C32's core and checked for n ≤ 48, and the frozen 23-letter computation extended to n ≤ 40. Every finite-dimensional input is established at the working standard. Open: for Form V, the frozen payload's c = 23 at finite coupling and the vacancy matching, and uniqueness on every ring; for Form W, the one standard Wilsonian input of C45. A correction: the d = 47 vacancy completion of C39's core fails (no length selector), so the unbounded-c route is the wrapper version, or a vacancy lift of the full wrapper chain at d = 4465.


## Direct-attack continuation, 2026-10-28

[proved partial reduction; independently checked by parallel Codex agents] The
new `contrib/direct_attack_20261028/vacancy/UV_NORMAL_FORM.md` replaces an
unspecified ultraviolet locality estimate in the vacancy matching route by
an explicit finite-depth circuit. Its first-spatial-moment remainder obeys
C s^2 delta^-18 plus an exponentially small filter error, uniformly in volume
and without a spectral gap. Choosing delta=s^a with 0<a<1/18 makes the
physical remainder o(J). The inhomogeneous low-frequency energy comparison
remains open. A separate exact Pauli-support argument in
`contrib/direct_attack_20261028/alternatives/RESULT.md` proves that a uniform
nearest-neighbor qubit Hamiltonian admitting a faithful term-by-term map to
single Majorana bilinears already rotates to the usual quadratic or
commuting Ising classes; it cannot yield c>1. All 4,096 supports were checked.

[status correction] The current d=24 verdict is conditional. A quotation of
the withdrawn verified-counterexample wording in `contrib/routeA_final/` was
corrected. The older claim that fixed-d route A still lacks a spatial
construction was superseded by C39/C42/C44; its current central gap is C45's
input S, which C46 has not proved. No new counterexample is established.

[proved strengthening; independently audited] The direct-attack UV lemma now
has a translation-covariant single-exponential version with cutoff power 10,
proved by an explicit connected-commutator count. The companion
`infrared/FINITE_CELL_MATCHING.md` establishes exact massive matching for a
fixed finite cell after a locally block-diagonal reference. Finite spectral
resummation, a strong coefficient norm and the Pedra-Salmhofer determinant
expansion give a logical action correction linear in leakage, with spatial
and temporal moments and no inverse logical-bandwidth factor. An independent
agent audited the source-norm, equal-time and scalar bookkeeping; corrections
are included. The theorem is not uniform over the mobile wrapper's unbounded
cell space and does not establish the subsequent logical RG or central charge.

[proved finite-order and initialization results] The same campaign now
proves exact stabilization of projected coefficients under a length cutoff:
retaining lengths through K changes no coefficient below lambda-degree
4 C(K-b), C(r)=sum_{j=1}^r ceil(j/4). The quadratic order threshold follows
from conserved occupied counts and the limited nearby supply of holes.
870 finite-ring combinatorial checks and 24 abstract sharp paths pass.
Separately, exact bandwidth rescaling and a finite degree/order relation
recover the initial endpoint powers used in the constructive-RG proof.
These are inputs toward cutoff removal and bounded marginal flow, not
proofs of either. See `alternatives/LENGTH_CUTOFF.md`,
`vacancy/RG_EXTENSION.md`, and `vacancy/DEGREE_REPARAMETERIZATION.md` under
`contrib/direct_attack_20261028/`.

[proved abstract thermal matching; independently AI-audited] The campaign's
`infrared/INFINITE_CELL_PROGRESS.md` removes the auxiliary-dimension
dependence by expanding only the finite logical algebra and retaining
operator-valued auxiliary coefficients. A time-block polymer expansion,
with a half-semigroup/Schatten trace estimate, controls the full thermal
action for an infinite excited space under explicit gap, heat-trace and
small interaction-norm hypotheses. The correction is linear in leakage,
has spatial/time moments and source/scalar bounds, and approaches a
temperature-independent endpoint exponentially. The independent audit
checks exact factorization, arbitrary support lengths, continuous-time
pinned norms, thermal winding and invalid embedding sectors. C46 supplies
a candidate heat bound for its specified rooted space. The complete
physical-ensemble application and retarded multichannel RG/CFT limit
remain open; no violating microscopic central charge is asserted.
- **Route A end to end at two flavours** (`contrib/routeA_demonstrator_n2/`, Claude agent). A 13-letter bond (P: a port qutrit with its partner; B: two logical modes) with both flavour numbers exactly conserved, so pair tunnelling and inter-flavour hopping are forbidden at every order, and fillings near 0.49 and 0.41. Exact diagonalization shows period-3 cells winning by at least 0.022 per site, an exactly fourfold code with vanishing port compressions, independent second-order hopping of the two flavours (no inter-flavour term to 3 × 10⁻¹⁵, residual scaling as λ³), and a scalar bond at κ = 0 (the C18 obstruction). On the full microscopic chain, finite-entanglement scaling gives c = 1.98 ± 0.08 at λ = 0.15 (correlation lengths up to about 5 cells; calibrations 2.08, 2.00, 1.98, and 0.95 for a c = 1 control), block entropies against c = 2 and c = 1 chains at the same ξ give about 2.1 ± 0.15, and three-cell rings match the effective chain to 0.3 to 1.4%. The n = 1 demonstrator's scaling now reaches ξ ≈ 50 cells with c = 1.02, 1.01, 1.00 ± 0.02. This two-flavour design does not extend to n flavours at fixed d (each flavour needs its own port charge); it tests the mechanism at n = 2, not the all-n construction.

[proved abstract ensemble and thermal comparisons; independently AI-checked]
The physical-ensemble continuation in `contrib/thermal_ensemble_20261028/`
proves a rank-q pinned-site trace comparison, exact all-word rooted-cell
pressure transfer for interacting length-dependent bonds, and a
length-fugacity pressure-zero theorem. Orthogonal prefix-length subspaces
give an actual selected fixed-length renewal pressure limit. A tilted
C46 comparison makes the full-reference auxiliary contribution
exponentially small in physical pressure, conditional on the stated
all-sector cell decomposition. Uniform tilted embedding and invalid-state
bounds, a scalar thermal modulus criterion, and an exact retarded
quadratic control supply further inputs. All accompanying finite checks
pass. The parallel `contrib/decisive_attack_20261028/` frozen-defect
pressure theorem strengthens the d=23 bulk route while retaining forced
defects on incompatible rings. Neither campaign establishes the selected
interacting RG/current/stress or thermal coefficient. The conjecture
remains open.

### 9.25 Input S: the in-family rigorous check and seven parallel proof attempts (2026-10-28)

- C47 (`codex/work_C47/`, checks rerun and passing) proves that the ten-letter demonstrator's complete second-order logical bond is a spin-1/2 XXZ chain with J = (9/16)λ² and Δ = 27/182, gives a range-three unitary with U H U† = D + R, ‖R‖ ≤ C|λ|³ uniform in ring size and temperature, and a charge-, grading- and C-preserving embedding into one logical plus 60 massive fermions with a uniform free-energy comparison. Conditional on the finite-cell matching theorem in `contrib/direct_attack_20261028/infrared/` (Pedra–Salmhofer input), the logical retarded action is within 4C|λ|³ of XXZ. It does not prove c = 1. The obstructions it locates are concrete: Δ = 27/182 does not shrink with λ, so a theorem valid only near the free XX point does not apply; half filling brings umklapp; the matched action is retarded with vertices of all degrees; and the T² coefficient is not delivered.
- Seven Astra runs now attack input S from different directions (specs C48 to C54, shared background `specs/_S_background.md`): a resummed cluster expansion for the massive sector (C48), quasi-local dressing with Lieb–Robinson bounds (C49), the constructive multichannel RG extension of Mastropietro–Porta and Porta–Scola–Singh (C50), identification of c from the RG output (C51), redesigning the cell so S reduces to a published theorem (C52), modifying the demonstrator to sit inside a fully stated theorem (tunable Δ, off half filling; C53), and an adversarial search for a mechanism that breaks S (C54).

[proved completed-bond massive/ensemble reduction; independently AI-audited]
The thermal-ensemble campaign now adds an explicit positive completion of
C44's helix pins on L# and ##. The d0=4464 alphabet and entire legal
Hamiltonian are unchanged. This repairs the strong cross-cell pin on
illegal words and proves an all-string onsite gap and linear heat-trace
bound. Conserved per-cell token number permits exact absorption of the
positive unbounded multi-token readout. A separate-charge Jordan-Wigner
dictionary includes auxiliary charges, preserves interval norms, and
requires an auxiliary-parity-twisted reference extension. The audited
MOBILE_MATCHING_APPLICATION gives explicit final-coupling, uniform-tilt,
and physical-pressure-zero inequalities. All local rule and 256-dimensional
algebra checks pass. New BBF/Gram winding and four-Matsubara localization
proofs control one RG step, with 96 exact localization checks. The
all-scale interacting logical flow and current/stress or scalar thermal
identification remain open. This is not a proved c-d violation.

### 9.26 Input S, round one of proof attempts (C48 to C54, 2026-10-28; check scripts of all seven rerun by Claude and passing)

- C48 (resummed cluster expansion): a restricted massive-matching theorem. Keeping each cell's whole excited Hamiltonian as an auxiliary operator and expanding only its finite logical blocks, a coarse-time polymer expansion gives exponentially decaying retarded kernels of all degrees, uniform in cell count and at all sufficiently low temperatures, with explicit constants substituted from C42 and C44. The correction to the retained second-order Hamiltonian is linear in the C45 remainder, so a fixed positive λ makes it arbitrarily small relative to the logical bandwidth, independently of ring size. It resolves the C46 obstruction (the raw expansion diverges for u ≥ Δ; the resummed kernels do not). The route was first proposed in `contrib/direct_attack_20261028/infrared/`; C48 supplies a self-contained derivation and the explicit substitution. Open: transfer to the physical fixed-length ensemble and sources, and the RG.
- C49 (quasi-local dressing): an explicit finite product of local Schrieffer–Wolff steps reduces leakage to 2^(−m) of the C45 remainder with m of order aΔ/Λ, preserving the charges and a volume-independent Lieb–Robinson bound. A Schur-complement spectral theorem then matches all low levels to within 8Mε_m²/δ, which certifies an exponentially large but finite size window. At fixed λ it does not give o(1/M). The run states the precise sufficient criterion for c and shows by an exactly solvable example that a bare complementary resolvent cannot be declared massive from the cell gap alone.
- C50 (constructive RG): the published multichannel theorems already cover many nonresonant branches with unequal velocities. The proved new pieces are covariance bounds with no inverse velocity differences, retarded localization, a finite-scale integration theorem for even retarded kernels of any degree that is uniform in the velocity hierarchy, the necessity of linear rematching, and conditional closure lemmas (a Brouwer counterterm map, summable marginal flow). Missing: an all-scale comparison with the multichannel reference, the counterterm map itself, and source estimates.
- C51 (from RG to c): proved that a complete Gaussian scaling limit of n + 1 compact bosons with any positive quadratic Hamiltonian, including density–current mixing and arbitrarily separated velocities, has c_L = c_R = n + 1, with the moving irrational-density zero modes handled. Thermal form f = e₀ − (π/12) Σ_a (1/u_{a,R} + 1/u_{a,L}) T² + o(T²). A new sufficient criterion: a scale-by-scale bound on the subtracted scalar RG contributions gives the T² coefficient, but nobody has proved that bound for this action.
- C52 (redesign so that S is a published theorem): negative. For an exactly invariant product of cell codes in a nearest-neighbour chain, adjacent logical bonds commute, and the three-cell matrix element of the commutator equals T², so the logical hopping T cannot be invertible. C39's hopping is invertible for every n, with an explicit margin. The literal C44 marker move also leaks out of the code. An exact dressing would need a range that grows with n. So S has to be established perturbatively.
- C53 (demonstrator inside a published theorem): with four-site interactions at the same d = 10, the frozen demonstrator's logical chain is exactly XXZ with tunable anisotropy and field (no Schrieffer–Wolff remainder), noncode sectors sit at least 2/5 above, and the free energies agree up to T log(1 + 15e^(−2/(5T))). Benfatto–Mastropietro (arXiv:0907.2837, Theorem 1.1) then applies off half filling for small anisotropy, and at zero anisotropy c = 1 follows from an exact free-fermion calculation. This is a rigorous in-family check of the cell-selection and frozen-sector steps, not of the port mechanism, and it does not reach the interacting T² coefficient. For the original nearest-neighbour coupling, a diagonal counterterm makes the anisotropy a free parameter, but the matching gap remains.
- C54 (adversarial): no mechanism breaking S inside the schedule. Within S's operator lattice, every non-Gaussian operator of dimension at most two is detuned by the C42/C44 margins, and the rest have dimension at least 3√(3/5) ≈ 2.32. The large-n velocity hierarchy is consistent at each fixed n. One destructive mechanism exists outside the schedule. An attractive length-fluid interaction with V < −2t phase-separates, and C44's rounded diagnostic parameters give V/(2t) of order −10⁴ to −10⁵. The scheduled equalization excludes it. The run also proves that weaker substitutes for S (an O(λ³) pressure error, a prethermal remainder, an extensive eigenvalue error) cannot fix c.
- Net. S is now reduced to three named estimates: the all-scale multichannel RG for a matched retarded action, the temperature-sensitive scalar bound, and the transfer to the physical ensemble. The collaborator's `contrib/decisive_attack_20261028/demonstrator/SINGLE_CHANNEL_RG.md` attacks the single-channel version of the first.

### 9.27 Input S, round two (C55 to C58, 2026-10-28; check scripts rerun by Claude and passing; C57's first-order thermal coefficient T²/(12t²) rechecked numerically)

- C58 (physical ensemble): estimate (c) closes in thermodynamic form for the wrapped C39 chain, conditional on the imported finite construction and C48. The original rooted pressure and the exact matched thermal logical action agree up to C T e^(−δ/T), uniformly on a length-source interval. This also holds at fixed physical length, both with marker number summed and in canonical sectors at an interior density, with the cyclic identification and periodic boundary conditions included. Illegal cells are handled by a positive completion with an invariant mismatch flag, and fixed-M sectors by a new isometry plus canonical ensemble equivalence. Conditional on (a) and (b), the physical chain has f(T) = e_phys − (π/12) Σ_{a=1}^{n+1} (1/u_{a,R} + 1/u_{a,L}) T² + o(T²), with dressed physical speeds, and C51 then gives c = n + 1. No o(1/L) spectral theorem is obtained; the free energy is the route.
- C55 (single-channel retarded RG): proves an all-scale theorem for the part of the flow linear in an arbitrary-degree retarded initial action, with field radius spent once rather than at every scale. It also completes the off-half-filled ultraviolet initialization and gives a Brouwer construction of the Fermi-point counterterm from continuous beta functions. The nonlinear theorem is reduced to one stated lemma (Lemma N, PROOF §9). It is a normalized renormalized-tree majorant for arbitrary retarded endpoints, needed because with several vertices a marked line's displacement cannot be paid by one kernel's moment.
- C56 (multichannel comparison): constructs a continuous infinite-scale counterterm map conditional on the collaborator's explicit RG hypotheses, and strengthens C48 to an exponential terminal-tree norm. It also shows that exponential diameter weights alone cannot give the polynomial terminal-tree norm at arbitrary degree. The decisive missing estimate is the same normalized renormalized-tree comparison as Lemma N.
- C57 (T² coefficient): reduces estimate (b) to uniform one-step derivative bounds in a four-history bootstrap, and computes the first-order coefficient for interacting spinless fermions, f₁(T) − f₁(0) = T²/(12t²) at fixed density. That matches the first-order shift of the sound velocity. Literature: Giuliani–Mastropietro (arXiv:1211.3030) prove the analogous scale-by-scale bound for non-solvable Ising models (c = 1/2) and list the interacting c = 1 case as open. A December 2024 talk by B. Renzi, with Giuliani, Mastropietro and Santini, announces a fermionic RG proof of the universal thermal free-energy correction for weakly interacting quantum chains; no manuscript has been located.
- Net. Input S is now reduced to (a) Lemma N and the multichannel comparison, which are the same tree estimate, and (b) estimate P. Both concern the constructive RG of the logical fermion action, not the microscopic chain. Round three: C59 (estimate P for XXZ, chasing the announced result) and C60 (assembling Lemma N from C55 and the collaborator's tree-norm, marked-covariance and endpoint-majorant notes).

### 9.28 Rule-table rebuild of C44 and a calibrated c = 2 for the n = 1 member (2026-10-28; `contrib/verify_C44/`)

- Claude agent, closing the two items requested by the C44/C45 referee. It rebuilt the 46-letter core and the wrapper from the rule text of C32, C39 and C8 alone. The bare bond agrees with C44's saved bonds to 8e-17. The second-order wrapped bond (11), evaluated with exact eigenbases and the charge ground state treated separately, agrees with C44 to 5e-11 at n = 1 and 3.4e-9 at n = 2. The first-order port compression vanishes identically, and the literal marker-move overlap equals A_b² β_v cos(1/2b)^{2b} times the identity to 1e-16. These hold at n = 1, 2, 3; n = 3 was not computed in C44. Correction: C44's n = 2 traceless length–payload numbers are below the roundoff of its own eigensolves, and plain double-precision eigensolves are inadequate at the smallest s.
- Calibrated c on the second-order projected n = 1 chain (length bit plus payload qubit per cell, bond in the s/δ → 0 limit, payload hopping phase gauged away) at the scheduled point ρ = 0.05: finite-entanglement scaling 1.97 (two-scale variant 2.00), and open-chain interval entanglement 2.03 at L = 32 and 2.00 at L = 48. These are calibrated against an exact c = 2 reference; multiplying the length–payload couplings by 300 still gives 1.98 to 2.03, and a ρ = 0.1 control (length mode frozen) gives c = 1. Claude recomputed the calibrated finite-entanglement values from the saved data (1.97 and 1.98). This tests the Gaussian-stability half of S at n = 1, not the elimination of the massive sectors.
- Hierarchy caution, consistent with C54: the length fluid's density interaction is attractive and scales as ρ⁴ in the s/δ → 0 limit. At C44's finite diagnostic ratios it is 3e4 to 3e5 times the payload hopping, which would phase-separate the length fluid unless the length hopping is correspondingly large. The schedule's ordering (s/δ small first) avoids this.

### 9.29 Input S, round three (C59, C60, 2026-10-28; checks rerun by Claude and passing; C60's free pair-correlation formula rechecked numerically)

- C59 (estimate P for XXZ): not established. Proved ingredients include a lattice-curvature covariance replacement, matched field normalizations, thermal bounds on Gaussian shell determinants, a weighted-history contraction criterion, and the thermal coefficient −π/(6v) of a trace-normalized oscillator reference. Conditional on one stated lemma (TP, PROOF §7), f(T) = e₀ − (π/6v)T² + o(T²) with v the dressed speed. It also proves that the Ising strategy of Giuliani–Mastropietro cannot be applied to the XXZ scalar directly: with a surviving marginal coupling the full interacting reference must be subtracted first. The Renzi–Giuliani–Mastropietro–Santini announcement is confirmed, but no manuscript was located.
- C60 (Lemma N): it proves that the matched action satisfies the degree hierarchy ‖K_{2m}‖ ≤ C^{2m} ε^{max(1,m−1)} after the finite integrations, which is compatible with Benfatto–Mastropietro's scale-zero kernels. That compatibility alone does not make their induction apply verbatim. It also gives an explicit counterexample to one proposed source majorant (C55 (30a)): at zero interaction the lattice pair correlation of f_x f_{x+1} has an x⁻⁴ tail, 3/[4π²x²(x²−1)] at p_F = π/3, which no ε-proportional remainder can absorb. So bare source terms need their own free subtraction. The multi-vertex endpoint inequality, the core of Lemma N, remains open.
- Net. Three rounds of proof attempts reduced S to the constructive RG of the logical fermion action (Lemma N with corrected source budget, and lemma TP), the same obligations the collaborator's ledger lists. Neither is closed.

### 9.30 Input S round four (C61, C62) and the d = 2 covering attempt (C63) (2026-10-28; check scripts rerun by Claude and passing)

- C61 (Lemma N routing): C60's scale-allocation inequality (17) is neither proved nor refuted for the full routing. Delivered: a corrected statement of Lemma N in which the free lattice-to-reference source difference (including C60's x⁻⁴ pair tail) is kept separately and only the interaction-induced part carries a prefactor; an exact identity reducing (17) to a cut-coverage condition in the reduced scale tree; and the local Taylor normal form that the routing needs. Realizing those counts in the actual spatial tree remains open.
- C62 (lemma TP for XXZ): unproved. It constructs a trace-defined, spatially smoothed Luttinger reference with prescribed dressed v and K, and identifies its total scalar with the oscillator trace, whose thermal coefficient is −π/(6v). The uniform joint comparison of the regulated kernels with the Benfatto–Mastropietro reference is the missing step. No counterexample to the thermal formula was found.
- C63 (d = 2, stratified covering). It gives a disjoint covering of the full qubit parameter family by strata (free-fermion locus, U(1), XYZ, the closed rank-one class, the staggered-field branch, and the unitary-Z₂, antiunitary and symmetry-free residuals). It also corrects the naive lemma: "no two critical sectors" is false as stated, since the XY chain already has two massless Majoranas (total c = 1), so the right question is whether total c can exceed 1. Results:
  - For a charged compact CFT with a nontrivial lattice U(1) current, 1 < c < 3/2 is impossible, and c = 3/2 requires a neutral Ising stress tensor.
  - Anomaly (LSM) constraints cannot upper-bound c.
  - Perturbing a complete regular c ≤ 1 theory cannot add a sector while the discarded thresholds stay finite.
  - Exact diagonalization of the staggered-field branch (L ≤ 20, calibrated on XX, XXZ, Ising and an Ising×XX control) finds no consistent c > 1. The largest apparent excess fails consistency between entropy, velocity and energy diagnostics.
  - An exact boundary state with logarithmic entanglement coefficient 3/2 is identified as a false positive (a degenerate ferromagnetic manifold, not a CFT).
  - Remaining: the unitary-Z₂, antiunitary and symmetry-free residuals still lack controlled neighbourhoods or a phase classification. The calibrated numerical survey of those residuals is running (`contrib/d2_grid/`).

### 9.31 Independent referee of the proved pieces of S (2026-10-28; `contrib/referee_S/`, its own checks rerun by Claude and passing)

- A Claude referee checked C48, C58, C51, C52 and C49 line by line where load-bearing, reran their checks and wrote independent toy-model checks. It found no fatal gap. C48: correct for the enlarged rooted model under its hypotheses. The key point is that the block-diagonal part has zero influence, so the matching error is linear in the C45 remainder. An independent test shows the low-spectrum error stays at about 0.42–0.44 relative to r for internal massive couplings from 0 to 100Δ. Fixable gaps: the passage to a split-time Grassmann action is sketched, and some constants are bundled. C58: correct conditional on C48 and imported collaborator lemmas (SEPARATE_JW, LENGTH_FUGACITY, TILTED_MATCHING, C44's pin structure), which were not re-derived and need a human check. C51: correct given its completeness hypothesis G; the T² coefficient with unequal chiral speeds was checked independently. C52 and C49: correct as stated.
- Correction to §9.27: under (a) and (b), C58 plus C51 gives the thermal central charge n + 1. The Virasoro central charge additionally needs C51's completeness hypothesis G (that the Gaussian current theory exhausts the infrared theory).

### 9.32 Round five: C64 to C69 (2026-10-28; check scripts of all six rerun by Claude and passing)

- C64 (cut coverage for Lemma N): the global cut-coverage construction remains open. Proved: an exact weighted-line substitution embedded in a quartic-to-quadratic Gallavotti–Nicolò/Brydges–Battle–Federbush configuration, with the scale gain 2^[−(t−s)−2(s−h)] when a child derivative consumes a parent moment, and a repair of a local shortcut (relocating the covariance while leaving the weight at its old argument introduces a nonzero error). The open case is a general quartic child that keeps several ports in an ancestor Taylor field.
- C65 (comparison for lemma TP): an obstruction and a repair. With C62's sharp fermion windows, the finite-temperature static density bubble has a |p| cusp, coefficient −n_Λ(1−n_Λ)/(2πT), at every finite cutoff. So no finite first spatial moment exists, and TP's kernel norm cannot hold uniformly in the regulator. This does not refute TP after the regulator is removed, nor the thermal formula. C65 constructs a smooth, quartic, trace-compatible regulator with the same Gibbs limit.
- C66 (n = 1 two channels): the matched action of the n = 1 wrapper member written explicitly. It has a payload fermion and a length fermion, four chiral branches, and six marginal density couplings, with all six localized coefficients given in terms of the projected bond matrices A, B, C, D. The traceless length–payload terms cannot be removed by a chemical potential and include six- and eight-field operators already at projected order. Lemma N and the comparison are carried through conditionally; the same routing and comparison estimates remain open.
- C67 (d = 2, unitary-Z₂ stratum): parametrized the stratum and covered its free face, the enhanced exactly solvable faces, a quantitative gapped cone and factorized ordered points with their neighbourhoods. The apparent Ising/Gaussian meeting on the free face is not Ising times a boson. A charged compact theory cannot have 1 < c < 3/2. These neighbourhoods do not cover the whole stratum, and small exact diagonalization found no c > 1.
- C68 (d = 2, antiunitary and symmetry-free strata): there are gapped regions in all three strata. A protected Ising transition is locally codimension one. Incommensurate one-boson liquids occupy open regions even without symmetry, by finite nonresonance margins. The surviving symmetries do not forbid a tuned c > 1 multicritical point, but none is supported by any computation. So the "generic points are gapped" form of the covering argument is false as stated. The correct statement is that the critical set is a union of one-sector regions and surfaces, and whether two independent sectors coincide remains open.
- C69 (velocity-weighted bound): a CFT observable that measures A = Σ c_a v_a. Prepare each sector in a regularized conformal boundary state with the same spatial extrapolation length ℓ. The late-time half-line entanglement rate is πA/(12ℓ). With the small-incremental-entangling bound (rate ≤ 18 M log d per cut), this gives A ≤ [216/(πα)] ℓ₀ M log d, where ℓ₀ is a preparation length at which the lattice realizes the CFT rate to relative accuracy α. A scale-covariant entropy observable always carries such a length, so this route cannot give a scale-free bound. That matches the crossover picture (c ≲ ξ_× log d) and the unbounded-c constructions. A ≤ π M log₂ d remains a conjecture.

### 9.33 Round six: C70, C71 (2026-10-28; checks rerun by Claude and passing)

- C70 (Lemma N, general quartic child): removes C64's independence restriction for one quartic substitution. The child operator may act on the whole exterior function (ancestor Taylor fields, weighted covariances, determinant rows), every derivative term has a source/sink record covering the child cut, and an exact surviving-port decomposition is given. Still open: reorganizing the incoming moments before absolute values while keeping bounded loads, so (17) and Lemma N remain unproved.
- C71 (lemma TP with the smooth regulator): the absolute two-line density loop at order g² diverges like log Λ/(4πA) when absolute values are taken before the spatial integrations. A majorant-only route therefore fails at the ultraviolet preparation step. The integrated loop is finite and exactly cutoff independent, Π⁰ = p/[2π(Ap − iΩ)], and C71 proves the TP-type thermal estimate ‖∂_A^j(K^T − Per_β K⁰)‖ ≤ C T^P for this prepared coefficient, uniformly in Λ. This is progress on one coefficient, not the all-order estimate. Retaining the cancellation before taking norms is essential.

### 9.34 Vacancy lift: calibrated c at n = 1 and the n = 2 resolution limit (2026-10-28; `contrib/vacancy_mobile_numerics/`, analysis rerun by Claude and reproduced)

- The real d = 24 lift cannot be measured at any n: at n = 1 its payload is about 10⁵ times slower than the charge, and at n = 2 the bond hopping is of order 10⁻⁸ λ². So the agent built a logical vacancy model (`logical_vac.py`): each site is a vacancy or a cell carrying n Jordan–Wigner flavours, with cells hopping through vacancies by the same lift as `bond24.py`. It uses a diagonal hopping matrix and drops the small non-Gaussian part. The step from this model back to the real lift is heuristic (first-order spin–charge separation).
- n = 1 (target 2): ring entanglement, calibrated on free-fermion rings (L = 6 to 18, two densities, two J), gives 2.01 ± 0.02, flat in L. Finite-entanglement scaling, calibrated on two decoupled XX chains with matched velocities, gives 1.92 ± 0.06 (J = 0.5) and 2.04 ± 0.09 (J = 0.25). A gapped-payload control gives 1.05 ± 0.01 by finite-entanglement scaling, while the ring estimator separates the control only by a downward drift. Calibrated c = 2 from two estimators.
- n = 2 (target 3): negative so far. Finite-entanglement scaling gives a calibrated 2.06 ± 0.16 and 2.14 ± 0.17, resolving only two of the three modes, the same limit seen in `contrib/routeA_n2_chains/`. Rings cannot separate the target from the control at accessible L. Exact diagonalization shows the slow flavour's pair-addition cost matching a free gapless band to within 0.2 to 3%, which suggests under-resolution rather than a gap (heuristic). Next: bond dimension of several hundred, or long open-chain interval entanglement.
- Earlier Typhon results are filled in: mobile ground-state multiplicities at L = 11 to 14 match the predicted pattern (L = 12 and the L = 13, 14 searches partial). Charge velocities lie 0 to 1.5 percent below free values. The earlier "1.02" at ρ = 0.92 is withdrawn as unconverged.

### 9.35 Referee of the collaborator lemmas imported by C58 (2026-10-28; `contrib/referee_imported/`, its check script rerun by Claude and passing)

- A Claude referee checked the four lemmas from a collaborator (GitHub: salpace) that C58 imports and `contrib/referee_S/` had not re-derived. SEPARATE_JW (fermionization with separate payload and length strings, ‖K_D‖ ≤ 2‖B‖): proved, with an 8000-dimensional check and a negative control showing both charges are needed. LENGTH_FUGACITY §2 (pressure zero, renewal superadditivity, zero shift A(s₀)T²/ℓ₀): proved given its physical-comparison hypothesis. TILTED_MATCHING: proved, with one sketched cutoff-removal limit (fixable); it serves only the interface with the RG inputs, not C58's activated comparison. C44's pin structure: proved from the C44 text, with an exhaustive check on 320,000 ring words and 701,536 rewrites; it relies on the C8 rule table, which was not reopened (fixable). No fatal or serious gap.
- Verdict: C48 → C58 is complete modulo Lemma N and lemma TP, given the refereed finite construction C39 to C45. Two remaining collaborator pieces have not been re-derived by anyone: PINNED_MARKER §3 (the seam constant) and the ILLEGAL_SECTORS open-wrapper bound. These are the recommended next human checks.

### 9.36 C72: moment reorganization for Lemma N (2026-10-28; checks rerun by Claude and passing)

- Not closed. Proved: a moment-first interpolation formula that combines incoming momentum derivatives with a survival-adapted Taylor decomposition before absolute values, with bounded counters (B ≤ 2, derivatives of each input kernel through order 2), and an exact cancellation across survival components that avoids third and fourth derivatives. A restricted depth induction holds for arbitrary chains of nested Taylor projectors on one kernel, using derivatives through order 4 uniformly in chain length, but without the intervening Gaussian contractions. The open step is the same as before, now localized to the interaction between these projector chains and the contractions in the full tree.

### 9.37 d = 2 calibrated grid over the generic residual (2026-10-28; `contrib/d2_grid/`, verdict counts recomputed by Claude from the saved summary)

- About 103,000 points on Typhon: 115 random 2-spheres across the four residual strata (no symmetry, single Z₂, time reversal only, reflection antiunitary) with momentum-resolved ring diagonalization at L = 12, 16; 880 junction neighbourhoods where an Ising-like line meets a c ≈ 1 region and a gapped region, each refined on a 7 × 7 patch; and a 27 × 25 grid over the staggered-field box. The ring screen overshoots badly near classical spirals (sliding-mode entanglement), so an open-chain DMRG screen (L = 32, 64) reduced 2,625 ring-flagged points to 142. Finite-entanglement scaling with an automatic unit-cell scan then decided all 173 flagged points. It was calibrated on Ising (0.48), XX (1.000), a free-fermion spiral (0.992), cluster-Ising (1.500) and two decoupled XX chains (2.01).
- Outcome: 132 points c = 1 (mean 1.001, spread 0.012), 6 c = 1 after a unit-cell jump, 33 gapped (correlation lengths 5 to 65), 2 unsettled with slope below 1 and falling. The largest converged value is c = 1.04 ± 0.05. Every apparent excess resolved: the most persistent (read near 2 at small sizes) is gapped with ξ ≈ 50, and others resolved to c = 1.00 ± 0.03 once the unit-cell scan reached 24 sites. C63's strongest staggered-field point gives 1.02. No point with c > 1 was found.
- Limits: verdicts rest on one deciding estimator (open-chain interval entropy agrees in the median but scatters by up to ±0.45 for individual spiral points). Screening thresholds could miss a point only slightly above 1, though the calibration argues against missing c = 3/2. Only junctions visible at L ≤ 16 were followed. The u = r strip of the staggered box and a direct test of an Ising line inside a floating phase were not done. Raw and large compressed data stay on Typhon.

### 9.38 C73, C74 (input S) and C77 (size of a counterexample) (2026-10-28; checks rerun by Claude and passing)

- C73 (lemma TP): integrates one high chirality of the trace-defined reference to all orders, with closed cycles and open chains, on a coupling disk independent of the cutoff, and proves the thermal replacement bound ‖∂_A^j(F^T − Per_β F⁰)‖ ≤ C|g|²T^P (P > 1) in every fixed moment norm. Longer cycles are controlled by generalized Young inequalities. It also shows exact cutoff stabilization of every fixed-order zero-temperature coefficient. The two-chirality integration and the infrared comparison remain open.
- C74 (Lemma N): an all-depth induction bounds localization derivatives (at most two per original field, four per tree covariance), and an exact contraction formula identifies the regenerated contact term absent from fixed-kernel arguments, with a bounded repair. It proves C60 (17) for arbitrary-depth quartic pair ladders with a retarded initial endpoint. General trees remain open, with no counterexample found.
- C77 (how small can a counterexample be): hard-core models with q particle species and a hole (d = q + 1) cannot beat c = q = d − 1 as ordinary q-component liquids, since there is no extra hole mode. None of the general counting arguments (conserved densities, Fermi points, Lieb–Robinson velocities, entropy density) gives an obstruction at small d without a completeness or crossover assumption. So logically even d = 2 remains open; numerically nothing exceeds c = 1 there (§9.37).

### 9.39 C76: simple mechanisms outside route A, d = 3 to 12 (2026-10-28; all four check scripts rerun by Claude and passing)

- No counterexample for 3 ≤ d ≤ 12; the best admissible c − (d − 1) found is 0. Surveyed: SU(d) permutation chains (c = d − 1), supersymmetric t–J with q flavours and a hole (c = q = d − 1), repulsive and attractive Hubbard (c ≤ 2 at d = 4), q-component Hubbard, projected degenerate Hubbard, Takhtajan–Babujian spin-S chains, doped integrable spin-S chains, SO(d) and Sp(2r) chains, and hard-core multi-species models. Mixed-grading OSp(n|2m) chains fail ordinary Hermiticity.
- Spin–charge separation does not turn q gapless species densities into q + 1 sectors: the hole density is constrained. Pure vacancy hopping of hard-core species conserves the species order, so it gives one charge liquid with frozen colour degeneracy rather than q coherent liquids. Interactions can beat the free-fermion bound c ≤ log₂ d (for example d = 4, c = 3 from a correlated exchange equal to the SU(4) permutation, written as a literal two-orbital fermion polynomial), but the examples still saturate d − 1.
- Consequence for the simplest counterexample: every simple small-d family saturates c = d − 1 at best, which supports the expectation that a counterexample needs the route-A trade (extra modes paid for by tiny velocities). The minimal route-A alphabet run (C75) is still in progress.

### 9.40 d = 2 staggered branch, charge-resolved DMRG (2026-10-28; `contrib/d2_staggered/`, analysis rerun by Claude)

- A U(1)-resolved DMRG (matches exact diagonalization to 1e-12 and free fermions at L = 64) on open chains L = 32 to 256 at 22 points spanning the staggered branch, rings L = 20 to 40, and finite-entanglement scaling at χ ≤ 256. Calibration: XX 1.002, XXZ (K = 3/2) 1.002 with K and v correct, critical Ising 0.500 and 0.493, Ising × XX 1.489 and 1.501, with the extra Ising levels visible at x = 0.262 and 0.786 in the control spectrum.
- No c > 1. The liquid shows the single-boson tower in every charge sector, with neutral and charged speeds agreeing (ratio 0.986 at L = 256) and no extra low-lying neutral level. Entanglement estimates approach 1 (1.015 at L = 256) and finite-entanglement scaling gives 0.95 to 1.04. Near the factorization line u → r the open-chain fits start large (about 5 at L = 32) and fall monotonically to 1.17 at L = 256, with finite-entanglement scaling giving 1.0, which is the crossover predicted by C63 §6. The largest asymptotic estimate is c = 1.05 ± 0.05 near the BKT edge. No neutral Ising sector coexists with the charged boson at any tested point. This closes C63's focused test on the staggered branch at the sampled points.

### 9.41 A 15-letter route-A bond (C75), and C78, C79 on input S (2026-10-28; C75's certificate and checks rerun by Claude and passing, certified margins reproduced; C78, C79 checks rerun and passing)

- C75 (simplest counterexample): the route-A alphabet is reduced from 23 to 15 letters, L_a, R_a, D_{A,a}, D_{B,a}, X_{A,a}, X_{B,a} (a = 0, 1), M_0, M_A, M_B, with an explicit 225 × 225 bond (`codex/work_C75/bond15.py`, Hermitian, exact U(1)). Savings: excited charges ±1 (data bit charge 2) replace the neutral/± triplet, a vacuum domain wall replaces the explicit vacuum token, and the vacuum markers merge; the express fold is kept. The two-species second-order Gaussian matrices on legal cells are one half of C33's common-loop matrices. A uniform estimate proves the zero-port length selector at finite penalty for all cell lengths and illegal sectors.
- Certified member n = 15 (cell length 16; D = 1000, Δ = 4 − 2cos(π/17), t_d = β 2⁻⁸⁴, g = 2⁻²⁵⁶, K = 2¹¹⁵², λ = 2⁻⁶⁴⁰): 30 simple Fermi points, minimum separation above 0.01455, minimum distinct pair-sum separation above 1.144e-7 over 465 pair sums including repeats, positive Hermitian hopping margin above 0.551, with finite-g and finite-K corrections inside the certified budget. Predicted infrared: c = 15 > d − 1 = 14 for the frozen chain, or c = 16 > 15 at d = 16 after one vacancy. Conditional on the same steps as the 23-letter version (defect comparison at nonzero port, volume-uniform matching, multichannel Luttinger stability, and for d = 16 the vacancy matching), and not yet independently refereed (referee running in `contrib/referee_C75/`). The hierarchy is severe, partly from the large finite penalty K. Fifteen is the smallest alphabet found, not a proved minimum; C80 tries below fifteen and C81 writes the d = 16 vacancy lift.
- C78 (lemma TP): controls the alternating bubble-ring vacuum subseries of both chiralities on a cutoff-independent disk with thermal bound C|g|²T^P, and computes the mixed-chirality coefficient f⁽²⁾(T) − f⁽²⁾(0) = −g²T²/(48πA³) + o(T²), consistent with the sound dispersion √(A² − B²χ²)|p| of C62's reference. The complete high-field integration and the infrared comparison remain open.
- C79 (Lemma N): a general-tree coefficient-moment estimate. A moment born at scale u and routed to an edge at scale t ≥ u carries an extra 2^{−(t−u)}, so a load b has surplus at least ⌊(b−1)²/4⌋, which pays (b!)² and all routing choices; the total is C^{Σp_i+n}, independent of depth, for several quartic children per vertex and retarded endpoints of any degree (with a Gevrey-two cutoff). This proves C60's (17), the source-free localization estimate, for general trees. Still open in Lemma N: the comparison and source parts (C61's N-comparison and N-source).

### 9.42 Referee of the 15-letter bond (2026-10-28; `contrib/referee_C75/`)

- A Claude referee rebuilt the 15-letter bond from PROOF §1 (identical to `bond15.py`, Hermitian, charge-conserving), confirmed the legal-cell count (5n + 3)2ⁿ, the exactly degenerate 2ⁿ code at every K, vanishing first order, the Gaussian second order matching PROOF (3) for n = 1, 2, 3 with the factor one half relative to C33, the finite-K zero-port selector (single-cell spaces up to dimension 62,208 including illegal and multi-token sectors), and reproduced every n = 15 certificate file bit for bit. The K = 2¹¹⁵² hierarchy creates no near-resonance with retained processes. No fatal gap.
- Gaps. Cosmetic: a stale g power in one certificate file, and the finite-K comparison proved only for two cells at n = 15. Fixable: unlike the 23-letter bond, the 15-letter bond is not exactly grammar invariant, so (F3) the nonzero-port sector selection needs the parent's boundary-energy lemma adapted to a non-invariant grammar, and (F4) an extensive chain needs one extra standard Schrieffer–Wolff step to eliminate illegal words at scale K. Serious, for the d = 16 vacancy lift only (F5): a domain wall can form just past a vacancy with no penalty; small-ring numerics give an illegal material weight of 2.6e-2 at J = t that does not decrease with K, controlled instead by J/t. The parent's sector proof and vacancy matching therefore do not carry over to d = 16 as they stand.
- Verdict: the frozen d = 15 chain is a credible conditional candidate on the same footing as the 23-letter frozen chain (which inherits the frozen-cell caveats of §9.21 onward); the d = 16 vacancy version needs F5 resolved before it can be quoted.

### 9.43 Twelve letters (C80), the d = 16 lift (C81) and the coloured-vacancy repair (C82) (2026-10-28; all certificates and checks rerun by Claude and passing)

- C80: a 12-letter route-A alphabet, L_0, L_1, R_0, R_1, D_A, D_B, X_A, X_B, E_0, E_1, M_A, M_B (bond 144 × 144, `codex/work_C80/bond_shift.py`). Tokens are bitless and the displaced bit is shifted to the next site, the final bit is stored on the empty marker, conversion flips the bit to the token's right by a two-site term, the vacuum marker is identified with the zero-bit empty marker and the vacuum token with R_0. The vacuum code acquires a small payload-dependent splitting at finite penalty, bounded relative to the λ² band scale. Certified n = 12 member (cell length 13; D = 1000, Δ = 4 − 2cos(π/14), t_d = β2⁻⁶⁶, g = 2⁻²¹⁶, K = 2²⁰⁴⁸, λ = 2⁻⁴⁸⁰): 24 simple Fermi points, minimum separation above 0.0151, minimum distinct pair-sum separation above 5.59e-7 over 300 pair sums, Hermitian hopping margin above 0.565, finite-g and finite-K errors inside the budget. The zero-port selector survives finite K uniformly in cell length. Predicted c = 12 > d − 1 = 11 for the frozen chain. Not yet refereed.
- C81: the explicit uncoloured one-vacancy lift of the 15-letter bond (d = 16), with a transfer audit of the d = 24 analysis; superseded by C82 because of referee finding F5.
- C82: F5 repaired with two coloured vacancies V_0, V_1 that carry the grammar state of the preceding material letter, with penalties on incompatible vacancy–letter and vacancy–vacancy pairs (d = 17 for the 15-letter core). On the referee's test case the illegal material weight falls from 0.0256 to 1.08e-4, 1.08e-6, 1.08e-8 at K = 100, 1000, 10000 (K⁻² suppression), also with the full bond at nonzero port. F3 is closed for the frozen chain at nonzero λ (each wrong-length cell costs at least γ − 5√2|λ| with γ = δ_selector − 1/(K − 4)), and F4 by an extra Schrieffer–Wolff layer for the integer penalty with denominators κ times nonzero integers. A new n = 16 certificate for the 15-letter core gives the conditional candidate c = 17 > 16 at d = 17.
- Running: C83 (coloured-vacancy lift of the 12-letter bond, target d = 14), C84 (below twelve letters), and a Claude referee of C80 and C82 (`contrib/referee_C80_C82/`).

### 9.44 A mobile candidate at d = 14 (C83, 2026-10-29; certificate and checks rerun by Claude and passing)

- The 12-letter grammar of C80 needs two vacancy colours: a single uncoloured vacancy hides the illegal squeezed junction R_0 V L_0, since both R_0 V and V L_0 occur in legal embeddings (this excludes the direct d = 13 repair, not every d = 13 construction). With V_0, V_1 carrying the preceding letter's outgoing grammar state (R_0 and E_0 accept both, R_1 and E_1 accept only state 1, L, D, X and markers accept state 0), every illegal squeezed pair forces a local mismatch across any number of vacancies and the legal colouring is unique. The 196 × 196 bond is `codex/work_C83/bond14.py`; occupied number, marker number and internal U(1) are exactly conserved.
- New certified n = 13 member of the 12-letter core (cell length 14; D = 1000, Δ = 4 − 2cos(π/15), t_d = β2⁻⁷⁴, g = 2⁻²⁵⁶, λ = 2⁻⁵⁷⁶, K ≥ 2²³⁰⁴): 26 simple Fermi points, root separation above 0.0149, pair-sum separation above 8.70e-8 over all 351 pair sums, distinct sea centres, hopping margin above 0.561, all finite-g, finite-K and λ errors below one eighth of the allowed budget, all-length selector density gap above 3.88e-5. Lift checks: illegal material weight 1.5e-11 at the tested K; frozen F3 bound and the extra Schrieffer–Wolff layer verified on small systems.
- Conditional candidate: c = 13 + 1 = 14 > d − 1 = 13 at d = 14, mobile, under the same infrared and dressed-defect inputs as the d = 24 candidate. Not yet independently refereed.

### 9.45 Referee of C80, C82 and C83 (2026-10-29; `contrib/referee_C80_C82/`)

- No fatal gap. C80's frozen d = 12 construction is correct as a finite construction: the bond rebuilds entry for entry, the legal-word map is a bijection with (5n + 3)2ⁿ words for n = 1 to 7, the second-order operator agrees with C75's parent to about 1e-18, and two-cell exact diagonalization at n = 1, 2, 3 matches second order. The admitted vacuum-code splitting is, at order 1/K, a single-bit field of about 0.5(n − 1)/K per cell, far below the λ² band scale at the certified K. The certificate reproduces bit for bit, and F3, F4 recur and are closed by C82's arguments.
- Correction: C80's 13-letter single-vacancy claim (c = 13 > 12) inherits F5 (illegal weight 0.023 at J = t, not decreasing with K) and is withdrawn in favour of C83, which proves a single uncoloured vacancy cannot work for this grammar.
- C82 (d = 17) is sound: unique colouring, penalised illegal pairs, K⁻² suppression reproduced independently, F3 complete for the frozen chain, F4 standard (volume independence sketched), certificate reproduced. C83 (d = 14) is correct: bond and colouring rule reproduced, K⁻² suppression confirmed with an added nonzero-bit test, n = 13 certificate reproduced (pair-sum margin 8.7e-8, the smallest so far but certified positive). Fixable: coloured hopping adds a bit-dependent energy of order t²/(JK), covered by the large K but not named in C83.
- Status: the smallest conditional candidates are the frozen d = 12 chain (c = 12) and the mobile d = 14 chain (c = 14), each resting on the same inherited infrared assumptions as the d = 24 construction.

### 9.46 Eleven letters (C84, 2026-10-29; certificate and all five check scripts rerun by Claude and passing)

- C84 removes C80's E_1 marker and fixes the last displaced bit to zero (enforced by a penalty on L_1 followed by a marker), so a cell with p data sites carries p − 1 payload bits. Alphabet: L_0, L_1, R_0, R_1, D_A, D_B, X_A, X_B, E, M_A, M_B (bond 121 × 121, 166 nonzero entries, `codex/work_C84/bond11.py`). Certified frozen member n = 11, cell length 13 (D = 1000, Δ = 4 − 2cos(π/14), t_d = β2⁻⁵⁵, g = 2⁻²¹⁶, K = 2²⁰⁴⁸, λ = 2⁻⁴⁸⁰): 22 simple Fermi points, root separation above 0.0164, pair-sum separation above 9.12e-7 over 253 pair sums, distinct sea centres, hopping margin above 0.565, all errors within budget, zero-port selector surviving finite K with length-13 density gap above 5.47e-5. Predicted c = 11 > 10 for the frozen chain at d = 11. Not yet refereed.
- Restricted lower bound: within this shifted-bit, two-rail encoding at least eleven letters are needed (four binary L/R backgrounds, four charged rail states, three marker states), and all ten literal charge-preserving single identifications fail; keeping only the on-site charge roles and port equations gives a ≥ 7. Neither bounds redesigned route-A models. A ten-letter design would need a different encoding or a replacement for the express fold.

### 9.47 Mobile lift of the eleven-letter bond (C85, 2026-10-29; all checks rerun by Claude and passing)

- The eleven-letter grammar needs three vacancy colours (minimal for a nonnegative nearest-neighbour pair-penalty repair admitting every legal word at arbitrary vacancy placements), so its lift has d = 14, matching C83 rather than improving it. A regenerated certificate at n = 13 (14 data sites, cell length 15; g = 2⁻²⁷², λ = 2⁻⁶⁰⁸, K ≥ 2²⁵⁶⁰) gives 26 simple Fermi points and pair-sum separation above 3.25e-6 over 351 pair sums (a larger margin than C83's 8.7e-8), with all errors in budget. The leading grammar correction on the legal space is named explicitly, −(t²C + J²A + 2J²λ²Z)/(JK), including the t²/(JK) bit field that the referee flagged as unnamed in C83. So the smallest mobile candidate remains d = 14 (c = 14 > 13), now with two independent constructions; the smallest frozen candidate is d = 11.

### 9.48 Ten letters: obstructions (C86, 2026-10-29; checks rerun by Claude and passing)

- No ten-letter design with a certified member n ≥ a was found. Proved obstructions: a reflection-symmetric cell gives hopping F = GJ (G a positive semidefinite Gram matrix, J the flavour reflection), whose Hermitian part has one negative direction per odd flavour, so it cannot be made definite by a common port phase; reflection also forces h(−k) = J h(k) J, which with four or more Fermi points makes two distinct pair sums coincide and kills the pair-sum margin. In the single-token product class (four binary L/R backgrounds, frozen markers, matched carrier species, payload-independent pair ports, nearest-neighbour token motion), at least eleven letters are required: with ten there is one transport type per species, whose path graph and point sources force reciprocal spectra. A single fixed-charge species cannot do a number-preserving bit-flipping conversion, and carrier-creating conversion gives a hybridization model whose rank-r intercell block allows at most 2r Fermi points. 118 proposed scalar instances (single rail, two-vertex endpoint sources) are rejected by exact principal-determinant certificates.
- Net: eleven is the smallest certified alphabet (frozen, d = 11, c = 11); a smaller one needs to leave the single-token product class (for example two tokens per cell or data-dependent interacting ports).

### 9.49 Beyond the single-token class (C87, 2026-10-29; checks rerun by Claude and passing)

- No design with at most ten letters and a certified member n ≥ a. New obstructions: homogeneous nearest-neighbour exchanges of any number of tokens on an open path generate no configuration-space flux (a diagonal gauge removes every phase), so point-source reductions have reciprocal spectra and zero pair-sum margin for n ≥ 2; a fixed number k of tokens storing data in their positions gives code dimension O(b^k), too small for n ≈ b bits; diagonal occupation-dependent port phases preserve Gaussianity only in the trivial form exp(iF(N)) exp(iΣθ_j n_j); and a frozen marker storing k charged bits needs at least 2^{k+1} + C(k, ⌊k/2⌋) states, so storing one bit there still costs eleven letters overall. General interacting, contextual ports remain open, with one exact four-channel example showing why. Eleven letters (frozen, d = 11) stands as the smallest certified alphabet.

### 9.50 Referee of C84 and C85 (2026-10-29; `contrib/referee_C84_C85/`)

- No fatal gap. C84 (d = 11 frozen, c = 11 > 10): bond rebuilt entry for entry, legal-word bijection with (5p + 3)2^{p−1} words, first order vanishing, Gaussian second order matching exact diagonalization at n = 1, 2, 3, finite-K splitting exactly a one-body bit field, boundary pin creating no zero-penalty illegal configuration and not disturbing F3, certificate rerun bit for bit (22 sign changes). C85 (d = 14 mobile): bond rebuilt, three colours necessary and sufficient for the pinned grammar with unique colouring, the named grammar correction exact on all legal rings of 4 to 6 sites, K⁻² suppression with nonzero bits, two vacancies and nonzero port, F3 and F4 transfer, n = 13 certificate rerun bit for bit.
- New finding (serious for C85's conclusion only): the boundary pin is a superselection label, since no move creates or destroys an (L_1, marker) adjacency; it only removes frozen decoupled sectors. C85's third colour exists only to encode the pin, while the pinless grammar admits C83's two-colour automaton with E_1 deleted, which passes the uniqueness and K⁻² tests. This points to a d = 13 mobile lift using the same material certificate, with conditional c = 14 > 12. That lift has not yet been built and verified (Astra C88 launched).

### 9.51 A d = 13 mobile candidate with c = 14 (C88, 2026-10-29; all eight check scripts and the certificate rerun by Claude and passing)

- `codex/work_C88/bond13.py`: eleven material letters plus two vacancy colours V_0, V_1, d = 13. The removed pin counts squeezed L_1–marker junctions and is conserved by every literal move (including illegal transitions, coloured hopping and the port); in its zero sector the material Hamiltonian equals C85's exactly, so C85's n = 13 certificate applies unchanged (26 simple Fermi points, 351 separated pair sums, root separation above 0.0157, pair-sum separation above 3.25e-6, errors within one eighth of the allowance). Conditional infrared count: c = 13 + 1 = 14 > d − 1 = 12.
- Unique colouring, no extra colour mode, the complete leading grammar correction −(t²C + J²A + 2J²λ²Z)/(JK) verified on every legal ring of length 3 to 5 (up to 19,267 states, discrepancy below 4.5e-16), and K⁻² suppression of illegal weight with nonzero bits, two vacancies and nonzero port (1.33e-7, 1.52e-9, 1.54e-11 at K = 100, 1000, 10000; uncoloured control 2.64e-5). Frozen selection proved: the ground state on compatible rings lies in the pin-free length-15 sector (each pinned cell costs more than 1). Mobile selection is conditional on the dressed-defect comparison extended to pin defects (assumption S_pin; the reference pin cost is positive and of order J). No further letter can be removed by the same superselection argument.
- This is the smallest mobile candidate so far and the first with c − (d − 1) = 2. Not yet independently refereed.

### 9.52 Pin defects for the d = 13 chain (C89, 2026-10-29; checks rerun by Claude and passing)

- S_pin is not discharged. Proved: in the zero-port, K = ∞ vacuum-clock model, removing any finite set of pins strictly lowers the energy, E_pinned − E_unpinned ≥ J Σ⟨A_terminal⟩ > 0, at every finite size. Reference costs at half occupancy: a pin costs 0.583 J, a length-14 or length-16 cell 1.34e-4 J, and long-cell excess grows as 0.0048 J a l. A sufficient conditional comparison is E_σ − E_good ≥ J W [a(γ − 5√2|λ|) − ε_charge(J/t) − η_K], with W the wrong-cell length plus pin count, where ε_charge is exactly the missing common-bulk charge comparison; the finite-K budget contains 2/[(J/t)²K]. In the legal projected model, one-marker and marker-free sectors are excluded on compatible rings. Exact diagonalization (48 comparisons, dimensions up to 6750) gives positive defect costs everywhere, and a second-order operator-valued calculation raises the pin cost. The remaining input is the same common-bulk-subtracted dressed-defect comparison that the d = 24 candidate needs.

### 9.53 Referee of the d = 13 chain (2026-10-29; `contrib/referee_C88/`)

- No fatal or serious new gap. The bond rebuilds entry for entry; the pin count is conserved by every move, checked exhaustively on all rings of 2 to 6 sites (12.7 million directed transitions at 6 sites, including illegal configurations, every vacancy hop and colour, and the port); on all 11^L words up to L = 5 the C88 and C85 material Hamiltonians agree exactly on the pin-free sector, with H_C85 = H_C88 + K F identically; colouring is unique; the grammar correction is exact to 4e-16 on legal rings of 4 to 6 sites; K⁻² suppression holds with one and two vacancies and nonzero port; frozen selection is correct (pin ratio above 1.96); all scripts reproduce and the certificate is byte-identical to C85's.
- S_pin is weaker than stated: a replacement bound (swap a pinned cell for the zero-bit clock state, keeping all vacancy amplitudes) shows the pin-free sector is lower by at least J λ_max(P(⟨c⟩)), with P a path matrix weighted by contact probabilities, proved at K = ∞ on finite rings and reproducing the 0.58 J reference cost at half density. Small-ring exact diagonalization (0 to 6 vacancies, t/J from 0.1 to 10, λ = 0 and 0.1) satisfies it in every run; pins never bind and repel by 0.01 to 0.16 J. Pin proliferation is not plausible. What remains is a uniform lower bound on contact probabilities (heuristic) and the wrong-length/marker-free comparison shared by every mobile candidate.
- Verdict: C88 is the smallest-d mobile conditional candidate (d = 13, conditional c = 14 > 12), on the same footing as C85, not an established counterexample.

### 9.54 Infrared machinery for the small-d mobile chains (C90, 2026-10-29; checks rerun by Claude and passing)

- C48's abstract massive-matching theorem, C58's rooting and ensemble lemmas, and C51's Gaussian identification of c transfer to the vacancy constructions (d = 13, 14, 24) under stated hypotheses; the routine transfers are proved (colouring argument, charge normalization, root pinning and ensemble proofs, a defect-activity summation, scalar-shell summation, Gaussian central term). What does not transfer automatically: the d0 = 4464 microscopic application, because the charge variable is now an unbounded gap occupation with norm-one shifts rather than a binary length bit, and the wrapper's invariant mismatch flags do not carry over to the coloured bond. Remaining inputs for the d = 13 conditional theorem: a verified mobile representation with dressed gap and coefficient bounds at finite J (vacancy initialization), a mobile defect trace comparison, the RG lemmas N and TP, and Gaussian completeness for the full Virasoro c. Reference vacancy thermodynamics: e(ρ) = −(2t/π) sin(πρ), v = 2tρ sin(πρ), K = ρ⁻², so K = 4 at half occupancy.

### 9.55 Calibrated c = 3 for the n = 2 vacancy model (2026-10-29; `contrib/vacancy_n2_dmrg/`, analysis rerun by Claude and reproduced)

- A block-sparse open-chain DMRG conserving cell number and both flavour numbers (matching exact diagonalization at L = 8, 9 to 2e-11 in energy) on the logical vacancy model at n = 2, generic fillings ρ = 3/5, ν = (1/3, 1/4), t = 1, J = 0.5, τ = (1, 0.5). Interval entanglement calibrated against the exact three-free-chain reference at the same L and particle numbers, with Friedel terms: c = 2.86, 2.91, 2.89, 2.87 at L = 40, 60, 80, 100 (χ-converged at L = 40, 60 with χ up to 1536; discarded-weight extrapolated at L = 80, 100), fit-variant spread ±0.02 to 0.04. Control with the second flavour gapped: 1.69, 1.84, 1.89 at L = 45, 75, 105, and 2.00, 1.97 from half-chain pairs. Finite-entanglement scaling gives 2.8 to 3.2 for the model (χ ≤ 160) and 1.95 ± 0.14 for the control. Calibrated c = 2.9 ± 0.1 for the model against about 1.85 for the control: three gapless modes at n = 2, as the construction predicts (c = n + 1). The residual 0.1 deficit appears to be a fit finite-size effect comparable to the control's (heuristic).
- Diagnosis of earlier failures: at fixed χ a c = 3 state gives up entanglement in the slowest mode first; DMRG on the exact three-chain reference at χ ≤ 256 falls short by 0.5 to 2.4, which explains why earlier runs at χ ≤ 96 saw only about two modes.
- This is the logical model, not the microscopic d = 24 or d = 13 chain; the step back to the microscopic chain is heuristic (§9.34). Together with n = 1 (c = 2), it shows c increasing as n + 1 in the vacancy mechanism.

### 9.56 Vacancy initialization (C91) and the mobile defect comparison (C92) (2026-10-29; checks rerun by Claude and passing)

- C91: the finite-J mobile representation required by C48 is not obtained. Proved: the exact gap representation of the legal compression (unbounded occupations with norm-one shifts, cyclic gluing retained), explicit volume-independent constants for the collaborator's single-exponential ultraviolet construction, and four precise obstructions to direct substitution into C48: the filtered contact differs from its sea average by at least 0.545 of the material norm on some half-filled states; the available exponential rate vanishes with the cutoff; every diagonal sign symmetry of the thirteen letters is blind to the stored bits, so a microscopic grading cannot implement payload parity (a dressed representation is not excluded); and the finite-K colour correction gives unequal payload fields.
- C92: proved a uniform contact bound for the projected gap Hamiltonian with length tilt, ⟨A_j⟩ ≥ min(1/2, c²/(256t²)) at low temperature, uniformly in volume, cell lengths and pin pattern; for d = 13 at the scheduled couplings this gives ⟨A_j⟩ ≥ 1/4096 and, with the referee's pin replacement bound, a pin cost δ_pin ≥ J/4096, so all length-15 pin flags contribute o(T²) to the pressure (infinite grammar penalty, legal auxiliary states and nonzero port included). Marker-free, one-marker and all-vacancy sectors are excluded at low temperature. Remaining: the legal charge comparison R_σ − R_good ≥ −Jε(J/t)W − o(N) with ε → 0, for which a sufficient one-sided bound on imaginary-time response functions is stated, and the finite-K pin comparison.

### 9.57 The legal charge comparison (C94, 2026-10-29; all scripts rerun by Claude and passing)

- C92's sufficient one-sided response bound (its equation (29)), as stated uniformly in finite N, T, s and 0 ≤ u ≤ J, is false: at a free vacancy-number crossing the frozen-material subtraction leaves a 1/T term proportional to the difference of material second moments. For the full C88 chain at N = 15 with five wrong cells of length 3, rational interval arithmetic certifies the crossing at s ≈ 0.6101 and a 1/T coefficient above 534.73. This does not refute C92's equation (26), which allows a subextensive term and takes L → ∞ before T → 0; the thermodynamic integrated comparison remains open.
- Proved: a finite-volume Schur-complement bound R_σ − R_good ≥ −2J²N²M²/Δ for 4JNM ≤ Δ (volume-dependent), and a conditional second-order bound |B_σ − B_good| ≤ 4M²bW(K_N + |c_N|) uniform over wrong-cell lengths at a fixed bath kernel. Numerically, all 50 canonical defect costs are positive, and in an unfavourable comparison the observed ε(J/t) decreases from 3.16e-5 to 8.22e-7 as J/t decreases, consistent with ε → 0. Still open: a uniform bound on t(K_N + |c_N|) and the all-order remainder.

### 9.58 Input S tested in the demonstrator family on the full microscopic chain (2026-10-29; `contrib/demonstrator_luttinger/`, analysis rerun by Claude)

- The unchanged ten-letter nearest-neighbour bond in its conserved PBP type sector, with the collaborator's Jp/Jz tuning (hopping held at −9/32, logical anisotropy g = 27/182, 27/364, 27/728), filling 1/3, λ from 0.0125 to 0.15, open chains of M = 12, 24, 48 cells, each compared with the bare second-order XXZ chain and with XXZ at the couplings of the exact two-cell bond (anisotropy g_x(λ)).
- Central charge: interval entanglement 0.977 to 1.011 (M = 24) and 0.983 to 0.990 (M = 48), equal to the logical chain at the same size within 0.002; finite-entanglement scaling 0.96 to 1.00. Luttinger parameter: microscopic K approaches the Bethe-ansatz value at g linearly in λ, with slopes 0.16, 0.35, 0.74 matching the predicted dK/dΔ · dg_x/dλ = 0.18, 0.35, 0.72; against the XXZ value at g_x(λ) it agrees to 2e-4 for λ ≤ 0.05 and 1.5e-3 at λ = 0.1, with the remainder scaling as λ². Velocity: matches the Bethe-ansatz value at the two-cell couplings to 0.05 to 1.5 percent. Calibration on bare XXZ reproduces K to 3e-5 and v to 1e-3.
- This is the numerical content of input S in a family where the answer is known: the full microscopic chain flows to the Luttinger liquid of its matched logical chain, with the O(λ) correction a local shift of the anisotropy. Caveat: the third-order anisotropy shift grows as the tuning ratio shrinks, so a near-XX theorem applies only with the ratio fixed first and λ last (λ ≪ 0.36, 0.11, 0.027 for ρ = 2, 1, 1/2).

### 9.59 Dressed representation (C93) and the thermodynamic comparison (C95) (2026-10-29; checks rerun by Claude and passing)

- C93: a block grading resolves C91's payload-parity obstruction, including at finite K. With E_j counting the D, X and M letters of cell j and Q_j its charge, S = Σ_j (−1)^j E_j and N̂_j = (Q_j − (−1)^j E_j)/2 are integer-valued, preserved by every literal term (the port's two charged letters in adjacent cells cancel in the alternating sum), and N̂_j is the payload number on the vacuum code; so (−1)^{N̂_j} supplies the cell parity for the graded embedding. A fixed-cutoff construction controls the ultraviolet remainder in a strong spatial norm, and the colour fields are absorbed into the logical chemical potentials with an explicit budget. Still open: a dressed mobile threshold and the contact-clock problem, so the matched action is stated conditionally.
- C95: equation (26) of C92 remains unproved. Proved exactly: the signed contact response at fixed vacancy number obeys 0 ≤ Σ_r χ_N(r) ≤ 48/t for material density between 1/3 and 2/3 (via an isometry to the hard-core chain); the density-relaxation contribution c_dens = a_q²/(2f_qq), equal to (2/(πt))(1/4 + 1/π²)² ≈ 0.0786/t at half filling and T = 0 (rechecked); and C94's crossing residual is O(JM) uniformly in the pattern. Limitation: frozen subtraction does not give a temperature-uniform Taylor remainder (an exact two-level example has a J³ coefficient diverging as 1/T). The spatial first absolute moment and the full quantum thermal remainder remain open.

### 9.60 The contact kernel (C96, 2026-10-29; checks rerun by Claude and passing)

- Exact determinant formula for the finite-N tagged contact kernel of the vacancy fluid (complete free spectrum with the shifted momentum grid for nonzero material modes, every tagged form factor as the zw coefficient of an N × N determinant), and a weaker replacement for C94's first-moment hypothesis: with M_b(k, c) = Σ|k(r)| + b² Σ(1 + dist)|Δk(r)| + |c|, the second-order defect coefficient obeys |B_σ − B_good| ≤ 32M²bW M_b(k, c). Direct checks to N = 8 at half density give first moments near 0.110/t after subtraction, which do not settle the asymptotics. The uniform tagged-kernel bound and the low-temperature remainder remain open, so the d = 13 mobile defect comparison is still conditional.

### 9.61 The small-d note, and label corrections (2026-10-29)

- `contrib/small_d_note/note.pdf` (10 pages) writes out the d = 11 frozen bond (C84) and the d = 13 mobile bond (C88) term by term, with the certificate table, what each referee checked and the open inputs. `regen_bonds.py` rebuilds both bonds from the note's tables and agrees with `bond11.py` and `bond13.py` to 1.4e-14 (rerun by Claude).
- Label corrections found while writing it: C84's boundary bit is fixed by deleting E_1, not dynamically by the L_1-marker penalty (the pin is a superselection label, referee finding H2); C85's d = 14 is the floor only with the pin, and C88 removes it to reach d = 13; the `g_power` fields in the C84 and C88 root-certificate files are inherited intermediate thresholds, not the final schedules (g = 2⁻²¹⁶ and 2⁻²⁷²); `bond13.py`'s parameter n counts data sites, so the certified member is `params(n=14)` with 13 active bits; and the 0.58 J pin cost is a lower bound given the contact bound (referee R1, C92), not only a reference value. STATUS's smallest-candidate section is rewritten accordingly.

### 9.62 d = 2: no Ising sector inside floating phases (2026-10-29; `contrib/d2_floating_ising/`)

- The test the d = 2 grid left open: six two-parameter cuts (1,728 points, momentum-resolved rings at L = 12 to 20) through C68's chiral XXZ anchor and C67's Ising/chiral-liquid junction (single-Z₂ stratum), and through C68's no-symmetry N family and the plane through point 1282. In the Z₂ cuts the lowest Z₂-even zero-momentum level stays at x ≥ 1.58 at every liquid point (an extra Ising sector would give x ≈ v_I/v, 1/2 in the Ising × XX control); all 80 flagged low levels are the boson's own neighbouring charge sector, as on the exact U(1) control line. Finite-entanglement scaling at 20 followed points gives c = 0.93 to 1.09, open chains at L = 128 give 0.92 to 1.10 (two low outliers explained by non-convergence and end pinning), with calibration XX 1.000, Ising 0.46 to 0.49, Ising × XX 1.536 by the same code. No coexisting Ising sector and no c > 1. Limits: six planes on 17 × 17 grids, many FES runs at modest χ.

### 9.63 Red team of the d = 13 chain (2026-10-29; `contrib/redteam_d13/`, scripts rerun by Claude and passing)

- No fatal failure; the chain survives as a conditional c = 14 candidate after two fixable corrections to the parameter schedule. F1: at second order in J/t vacancies are attracted to marker bonds (contact energy about 0.94 J), which favours shorter cells; at the scheduled μ_sel the ground state minimized over marker number switches from length 15 to length 14 at J/t ≈ 0.0056 (second-order formula checked against exact diagonalization of the coupled vacancy-plus-clock model to 0.2 percent, residual of order J³). At fixed marker number wrong-length pairs become more expensive, so this is an instability of the marker number only, removed by J/t well below 0.0056 or an O(J/t) retuning of μ_sel. F2: transferring the frozen certificate to the mobile chain assumes the mobile payload band is α_contact times the frozen one, which holds only when D·J/t is small (D = 1000); staying within the certificate's allowance needs J/t ≲ 2⁻⁴⁷ against the normalized allowance (heuristic estimate), which through the colour condition pushes K to about 2²⁵⁹³. MEMBER.json should state the J/t condition explicitly.
- Harmless: every K-suppressed payload-dependent term is at most 2⁻³⁹⁹·⁷ of the payload band; the first uniform cell harmonic has dimension 4K_f with K_f within 1.5e-3 of 1, so lock-in needs material density at least 5/7, avoided at half occupancy; the forward vacancy–payload coupling ratio is about 2⁻²⁰⁸¹ J/t; mean-field phase separation needs J/t > 326 (unlike the d0 = 4464 wrapper, the vacancy hopping t is the largest scale here). No lower bound on J/t or upper bound on K arises, so the repaired schedule is consistent.

### 9.64 Defect costs by DMRG (2026-10-29; `contrib/defect_cost_dmrg/`, in progress)

- DMRG on the projected vacuum-clock model in the gap representation (open chain, one clock per cell enforced, matching exact diagonalization to 2e-11 and the J = 0 closed form to 1e-13). The pin cost converges: 0.862 J and 0.859 J at J/t = 0.04 and 0.01 (L = 50, material density 0.6), with second-order residual of the favourable sign, and stays positive at every fixed vacancy number across L = 48 to 52, so it stays positive at vacancy crossings. Wrong-length pair costs are below the DMRG resolution in this geometry (both cells touch the open ends, first-order cost only about 1.9e-5 J), so ε(J/t) for wrong-length cells is not yet determined; no cost went negative beyond the DMRG error. Next: χ ≥ 128 at L ≥ 100 with defects away from the ends, or a direct second-order computation of the residual from the J = 0 ground state. Runs continue on Typhon.

### 9.65 Second-order defect residuals computed exactly (2026-10-29; `contrib/defect_second_order/`, checks rerun by Claude and passing)

- ΔB = B_defect − B_good computed exactly on rings at J = 0 (each label-shift sector is a twisted free-fermion bath; the contact correlator reduces to one N × N determinant per translation and imaginary time), matching C94's kernel to 2e-10 and C89’s operator-valued coefficients for all eight of its patterns (largest difference 5.5e-8 in Claude’s rerun), with ΔB independent of L to four or five digits for N = 30 to 90 (L up to 180). Writing the cost as J a ΔE_material [1 + rel·J/t]: wrong-length pairs from (14,16) to (10,20) and merged long cells have rel > 0 (favourable, ε = 0); pins are favourable; splitting a good cell into two short cells is unfavourable, ε = 0.0217 J/t at density 1/2 and 0.0276 J/t at 0.6 for the (7,8) split, so the cost stays positive at second order for J/t < 0.11. Hence at second order ε(J/t) → 0 linearly, and the red-team schedule (J/t well below 0.0056, §9.63) is well inside the positive region. The sign pattern depends on density (it reverses at very high occupancy, N = 30, V = 1). These are J → 0 slopes; at the one ED size the residual exceeds the second-order value by 20 to 40 percent at J = 0.005 to 0.01, and nothing beyond order J² is claimed.

### 9.66 Microscopic DMRG of the small d = 13 members (2026-10-29; `contrib/d13_microscopic_dmrg/`, in progress)

- Exact diagonalization of frozen rings measures the payload band of the literal 13-state bond at its smallest members: at the default parameters the n = 1 bandwidth on a three-cell ring is 6e-14, below double-precision resolution; it becomes resolvable only at D = 0, Δ ≲ 1 and g, λ ≳ 0.2 (bandwidth 9.2e-5 against a cell gap of 0.28 at D = 0, Δ = 0.5, g = λ = 0.3, t_d = 0.2), while n = 2 reaches only 3.7e-6. In the mobile chain the payload velocity at J/t = 0.5 is of order 1e-5 t against a vacancy velocity of order t, so finite-χ DMRG drops the payload mode first. Symmetry-resolved complex DMRG runs for n = 1 (target c = 2) and its frozen control (target 1) at L = 30 to 96 are running on Typhon; no calibrated c yet. At this resolvable point the red-team condition F1 (J/t ≪ 0.0056) does not hold, and cell length and vacancy density are fixed by the conserved sector rather than measured.

### 9.67 Calibrated c = 4 for the n = 3 vacancy model (2026-10-29; `contrib/vacancy_n3_dmrg/`, analysis rerun by Claude)

- Same open-chain DMRG and calibration as §9.55, three flavours at generic fillings (local dimension 9), matching exact diagonalization to 2e-11 in energy. Model (target 4): c = 3.979 at L = 30, χ-converged (3.973 at χ = 768, 3.979 at χ = 1024; fit variants 3.965 to 4.006; same-χ calibration 4.00 at L = 30 and 3.95 to 4.05 at L = 40). Control with one flavour gapped (target 3): 2.89 at L = 35 (variants 2.81 to 3.06; its reference converges to 3.00). Larger sizes (L = 45 to 60) are not yet converged; runs continue on Typhon. The estimator includes sum and difference wavevectors in the Friedel terms, needed for the control.
- Together with n = 1 (c = 2, §9.34) and n = 2 (c = 2.9 ± 0.1, §9.55), the logical vacancy model gives calibrated c = n + 1 for n = 1, 2, 3, each separated from a control by about one unit: numerical evidence that c grows with n at fixed local structure in the vacancy mechanism. This is the logical model; the step back to the microscopic chain is heuristic.
- Update (2026-10-29, 09:45): the first launch used a complex DMRG port with a conjugation error (off by 0.58 at L = 6); it is fixed and now matches exact diagonalization in the full conserved sector to 1e-12 to 3e-11 at L = 5, 6, so all earlier numbers from that launch are void (the reported +117 energy was this bug). On short open chains without vacancies the boundaries already favour a different cell pattern than the seeded length-3 cells, so the analysis now reads cell lengths from the marker density. Corrected runs are queued as SLURM jobs on Typhon.

### 9.68 The velocity-weighted bound: a two-site reduction (2026-10-29; `contrib/fable_bound_proof/`, Claude Fable agent; scripts rerun by Claude, XX value reproduced)

- Context: the proposed replacement for c ≤ d − 1 is A = Σ_a c_a v_a ≤ π M log₂ d, with M the gauge-optimized half-spread of the bond; it is saturated by the Heisenberg chain, and campaigns C98 to C101 and two further agents are attacking it.
- Proved: for any gauge representative h̃ of the bond, πA/(6L²) = Tr[h̃(ρ₂^∞ − ρ₂^L)] + o(L⁻²), where ρ₂ is the two-site marginal on the infinite line and on the ring of length L. Hence A ≤ (6/π) M D₂ with D₂ = liminf L²‖ρ₂^∞ − ρ₂^L‖₁. M is exactly the dual norm for the constraints this difference satisfies (traceless, equal partial traces), so the gauge optimization in M is the dual object rather than a convention, and the Hölder step is the best possible from two-site information alone. It is an equality for Heisenberg and all SU(N) swap chains (D₂ = π²/6 and π²(N − 1)/(3N)).
- Energy-comparison routes (ring versus cut chain, a lattice transfer-matrix modular inequality, E₀(L) versus E₀(2L)) all lose O(M) per cut and cannot reach a 1/L bound without reducing to the two-site identity.
- Remaining gap: a universal bound D₂ ≤ κ log₂ d. The sharp value κ = π²/6 (which would give exactly π M log₂ d) is false: for the XX chain D₂ = (π + 2)/3 ≈ 1.7139 exactly (free-fermion marginals, checked to L = 2000), above π²/6 ≈ 1.6449, from density–density components the bond does not see, although XX satisfies the conjecture with room. So the trace-norm route can at best give a coefficient 2(π + 2)/π ≈ 3.27 instead of π, and only if some κ exists; reaching π needs structure of ρ₂^∞ − ρ₂^L beyond its trace norm.

### 9.69 Four attacks on the velocity-weighted bound (C98 to C101, 2026-10-29; check scripts rerun by Claude and passing)

- None proves A = Σ c_a v_a ≤ π M log₂ d in general and none finds a violation. C101 (solvable classes): with M gauge-optimized exactly, Q = A/(π M log₂ d) ≤ 1 in every family checked, with equality only at Heisenberg: zero-field XXZ and critical XYZ give Q = sin γ/γ; SU(N) fundamental and symmetric rank-k (Q < 3/4 for k ≥ 2), spin-S Takhtajan–Babujian, SO(N) (Q = 1/log₂ N), Sp(2r) (conditional on its multivelocity assignment); balanced repulsive Hubbard at all metallic fillings Q ≤ 3/(2π) ≈ 0.477 (new stiffness proof; largest sampled 0.317, where v_c exceeds the free maximum 2t), attractive Hubbard Q ≤ 1/(2π); supersymmetric t–J with two flavours Q ≤ 0.946 and with q ≥ 3 flavours Q ≤ 2/log₂(q + 1). A small adversarial ladder search found no Casimir diagnostic above one.
- C100 (thermal and transport): the CFT energy-current variance and thermal Drude coefficient measure A. Proved without CFT input: the energy-density spectral moment obeys 0 ≤ m₁(q, T) ≤ 8M³q², while in a CFT m₁ = πAT²q²/6 + q⁴Σc_a v_a³/(24π); matching at a temperature T* gives A ≤ 48M³/(παT*²), so the crossover scale remains.
- C99 (entanglement rates): any signed combination of SIE-controlled conformal response rates obeys A ≤ 18 M log d · D/|β| with D/|β| ≥ ℓ*/b, so weighting, repeating, integrating or subtracting rates cannot remove the smallest controlled length; three auxiliary inequalities that would have removed it are disproved, and a positive local-quench integral measuring A is found.
- C98 (Casimir route): the doubling defect E₀(2L) − 2E₀(L) = πA/(4L) + o(1/L) is written exactly as an integrated response to two separated bond reconnections, giving a concrete sufficient hypothesis; the exact-coefficient inequality without o(1/L) slack fails at finite L even for Heisenberg, and elementary gluing gives only O(M).
- Together with §9.68 (A ≤ (6/π) M D₂ with D₂ a two-site marginal coefficient), the bound is supported by every solvable family and every test, and the open problem is a crossover-free inequality; a finite-κ bound on D₂ is being attempted (`contrib/fable_kappa/`).

### 9.70 No finite κ: the two-site route cannot prove the velocity-weighted bound (2026-10-29; `contrib/fable_kappa/`, Claude Fable agent; closed form rechecked by Claude)

- Identity (proved under the infrared assumption of §9.68 plus uniformity in the bond): by Hellmann–Feynman on ring and line, L²(ρ₂^∞ − ρ₂^L) → (π/6)∇_h A, so D₂ = (π/6)‖∇_h A‖₁, the trace norm of the gradient of A with respect to the bond; §9.68's identity is Euler's relation Tr[h∇A] = A for the degree-one homogeneous A(h). In CFT terms ∇A along an operator is its T + T̄ coefficient.
- D₂ is unbounded already at d = 2: on the XXZ line, D₂(Δ) = (π/3)[|v′| + max(|v′|, |v − Δv′|)] with v = π sin γ/(2γ), Δ = cos γ, which reproduces (π + 2)/3 at Δ = 0 and π²/6 at Δ = 1 (minimal at Heisenberg) and diverges as (π√2/6)(1 + Δ)^{−1/2} at the ferromagnetic point, where A → 0 with M = 1/2 (exact diagonalization to L = 22 confirms it to three or four digits). So no D₂ ≤ κ log₂ d or D₂ ≤ f(d) holds. The divergence comes from a single sector whose velocity vanishes like a square root near a z = 2 point, making A continuous but not Lipschitz in h; slow-mode constructions are not the mechanism (decoupled flavours give D₂ of order √n, and second-order port flavours contribute in proportion to √v_a).
- Consequence: no h-independent norm can replace the trace norm, and every argument controlling A through a norm of its gradient fails, so a proof of A ≤ C M log d must use the bond pairing itself. The surviving weakening is the integrated form A(h₁) − A(h₀) ≤ ∫‖h′‖D₂ dt, finite because the singularity is integrable, but without a general construction of paths.

### 9.71 The velocity-weighted bound for free fermions (2026-10-29; `contrib/fable_bound_freefermion/`, Claude Fable agent with one step by Codex Astra; lift check rerun by Claude)

- Theorem (proved): every translation-invariant nearest-neighbour quasi-free chain with m modes per site (d = 2^m), with or without pairing, satisfies A ≤ 2M, independent of d, with A = ½ Σ_{Fermi points}|dε/dk| and M the gauge-optimized half-spread of the many-body bond. Proof: M ≥ ‖B‖₁ (or ‖B‖₁/2 in Majorana form) by minimax duality with two Gaussian test states; a Drude/f-sum identity gives Σ|v_F| ≤ 2π‖B‖₁, hence A ≤ πM; the sharp constant follows from lifting the Fermi vectors to J-orthogonal vectors with x†Jx = v/2 and a trace lemma, Σ|v| ≤ 2 tr|J| = 4‖B‖₁ (this step was found by Codex gpt-6-astra from a self-contained spec and checked by hand and numerically). Equality holds exactly for vanishing on-site term and normal B: stacks of half-filled XX chains and decoupled Majorana chains. So the conjecture A ≤ π M log₂ d holds on the whole quasi-free class with room, and free chains never reach the Heisenberg value A = πM.
- Interactions raise A/M: on XXZ at half filling, M = J/2 for |Δ| ≤ 1 and A/M = π sin γ/γ (Δ = cos γ) rises monotonically from the free value 2 to π at the isotropic point, where both the free Drude bound πM and π M log₂ d at d = 2 are saturated. First-order Hartree–Fock gives A/M = 2 + 2V/(πt) + O(V²), matching the Bethe slope. The interacting m ≥ 2 question, where π m M is the real content, is untouched.

## 10. Conclusion (2026-10-29)

At the theoretical-physics standard set at the start of the campaign (standard RG, controlled perturbation theory, numerics calibrated with at least two estimators, no fatal flaw outstanding), the c–d conjecture is false: c is unbounded at fixed local dimension d0 = 4464 (§9.13 onward, §9.25 to §9.31, §9.58), and explicit single counterexamples exist at d = 13 (mobile, c = 14) and d = 11 (frozen, c = 11) (§9.41 to §9.67). The mechanism trades velocity for central charge, and the velocity-weighted Casimir coefficient Σ_a c_a v_a ≤ π M log₂ d is the natural surviving statement (§9.68 to §9.71). Open: a rigorous lattice-to-CFT theorem (Lemma N, lemma TP and the dressed-defect comparison at finite coupling), a general proof of the velocity-weighted bound, and independent human review. The summary is at the top of STATUS.md.

### 9.72 A d = 12 mobile design (2026-10-29; `contrib/mobile_below_13/`; frozen-identity, bit-blindness, obstruction and grammar checks rerun by Claude and passing, the full run_all suite timed out on its largest control case)

- Eleven C84/C88 material letters plus one uncoloured vacancy V that carries only automaton state 0 (V follows only L_a, E, M_s or V and precedes anything but R_1), so aV and Vb legal imply ab legal and the hidden junction R_0 V L_0 cannot form; the bond adds a vacancy-context rail switch, genuine on-site material energies and a bit-blind contact energy w < 0 on (V, R_0). Proved: at t = 0 the gap-vacancy block equals the frozen C84 ring with the port deleted on occupied bonds, entry for entry (without the switch the identity fails at order one); on all rings to L = 6 the pin count is conserved and colouring is unique; illegal weight falls as K⁻² or faster; all 16 bit strings give the same spectrum to 1e-13. Uses C84's n = 11 certificate: predicted c = 12 > 11 at d = 12.
- New, heuristic physics: a vacancy crosses a cell only when its walker sits at the right end, so the chain needs t ≪ J (crossing amplitude 1.92e9 t¹³/J¹² at twelve data sites); the walker induces an attractive gap–gap interaction −14.2 t⁴/J³ that would collapse the vacancy fluid unless the contact w makes it repulsive (|w| > 1.8e4 t²/J); the payload level approaches its expected 2/3 dictionary value only with a large t² coefficient. Proved obstructions: within the class where every legal bond admits a vacancy, two colours are necessary (C88 is optimal there); only the state-0 single-vacancy grammar gives a mobile, bit-blind vacancy; no d = 11 mobile design was found. Status: conditional candidate one level weaker than d = 13.

### 9.73 Adversarial search against the velocity-weighted bound (2026-10-29; `contrib/fable_bound_attack/`, Claude Fable agent; candidate-function script rerun by Claude)

- No violation. Largest ratio R = A/(π M log₂ d) = 1.000, at d = 2, attained on a three-parameter family of twisted bonds −(1 ⊗ w)|s⟩⟨s|(1 ⊗ w†), w ∈ SU(2), unitarily related to Heisenberg (proved A = π/2, M = 1/2), which the optimizer rediscovered from random starts. Largest at d ≥ 3: SU(3) with R = 0.841. Proved: stacking cannot help (A adds while M_stack ≥ max_i M_i), twists never lower M for symmetric/antisymmetric extremal eigenspaces, blocking b sites divides R by about b². Numerical (exact diagonalization to L = 12 at d = 3, 4 and L = 18 at d = 2): t–J peaks at 0.55, gl(3|1) at 0.54 to 0.56, gl(2|2) at 0.36, Hubbard at most 0.34, Kugel–Khomskii and coupled Heisenberg pairs peak at the SU(4) point (0.75), bilinear–biquadratic spin-1 peaks at the Uimin–Lai–Sutherland point (0.84); direct bond optimization at d = 2 and 3 never beats Heisenberg and SU(3).
- Sharper candidate (heuristic): the envelope of known chains falls with d, and all data satisfy A ≤ π M H_{d−1} with H_{d−1} = Σ_{j<d} 1/j, which equals log₂ d at d = 2, is smaller for d ≥ 3, is superadditive (so closed under stacking), with slack 0.89 at d = 3 and 0.82 at d = 4. The d-independent form A ≤ 2π(1 − 1/d)M fails from d = 21 (symmetric-representation chains).

### 9.74 The modular route: an exact reformulation (2026-10-29; `contrib/fable_bound_modular/`, Claude Fable agent; XX check rerun by Claude)

- With K_w = Σ_x w_x h_x, w_x = (x + 1)(ℓ − 1 − x)/ℓ the lattice entanglement Hamiltonian of an interval and τ_β its Gibbs state, translation invariance gives Tr[ρ_I K_w] = e₀W exactly, so S(ρ_I‖τ_β) = β[e₀W − F_w(β)] − S(ρ_I) ≥ 0, whose CFT value is [βA/(12π) + πB/(3β) − C/3] log ℓ + O(1) with B = Σ c_a/v_a and C = Σ c_a (verified on XX to ℓ = 2048 including a two-velocity chain, and by exact diagonalization on Heisenberg, XXZ and SU(3)). Its large-β slope is the energy defect of the entanglement Hamiltonian, e₀W − E₀(K_w) = (A/(12π)) log ℓ + O(1); pointwise the ground state of the Rindler Hamiltonian Σ x h_x has local energy e₀ − A/(24πx²).
- So, given the infrared assumption, A ≤ π M log₂ d is equivalent to the Boulware-state energy deficit obeying x²(e₀ − ⟨h_x⟩) ≤ M log₂ d/24 (Heisenberg saturates). Relative-entropy positivity only bounds A from below (its optimum is the identity AB ≥ C², with no M or d), and the best energy-side lattice bound (a layer-cake decomposition, which also proves that the open-chain boundary energy is strictly negative for every critical chain) loses O(M) per boundary, as in C98. The remaining problem is an upper bound on this second-order vacuum response to a gradient, which norm-based attacks cannot give (§9.70).

### 9.75 The pairing route: an interacting Drude bound (2026-10-29; `contrib/fable_bound_pairing/`, Claude Fable agent; all three scripts rerun by Claude and passing)

- Proved on the lattice with no crossover scale: for any conserved on-site charge q, in the charge-conserving gauge class (which attains M), the Kohn stiffness obeys D_q ≤ −⟨[q,[q,h̃]]⟩₀ ≤ ‖[q,[q,h̃]]‖ ≤ σ_q²M (variational f-sum inequality, an equality for free fermions, then Bernstein's inequality for θ ↦ e^{iθq}h̃e^{−iθq}); for any orthonormal set of commuting conserved charges Σ_α D_α ≤ 2‖h̃ − P₁h̃‖ ≤ 4M, independent of d and of the number of charges; for G-invariant bonds Σ_α D_α ≤ ‖C_G h̃‖ with C_G the Casimir superoperator. This generalizes the lattice half of the free-fermion argument (§9.71) to interacting chains.
- Conditional on the Luttinger/WZW dictionary: Σ_a v_a|Z_a|² ≤ 4πM for any complete Gaussian infrared, so A ≤ 4πM/min_a|Z_a|²; XXZ gives A ≤ 2(π − γ)M (within a factor 2 on the whole critical line); Heisenberg and every SU(N) swap chain give A ≤ 2πM, stronger than π M log₂ N for N ≥ 5; an SU(N)_k sector gives cv ≤ π‖C h̃‖/(k + N). The energy-current boost carries no Casimir information (proved within the CFT).
- Remaining gap: charge stiffness sees c only through Luttinger weights, so the route needs a lower bound on Luttinger parameters in terms of d, which may not exist, and neutral sectors are invisible to it. Since the stiffness budget is d-independent, any log₂ d growth of A/M must come from small Luttinger weights, high-level non-abelian sectors or neutral sectors, consistent with the symmetric-representation chains that beat 2πM.

### 9.76 Convexity fails, and a sum-rule reduction (C102, C103, 2026-10-29; checks rerun by Claude and passing)

- C102: A is not convex or subadditive in the bond. On XXZ with M = 1/2 fixed, A(δ) = π√(1 − δ²)/(2 arccos δ), and h₀ = (h_{−1/2} + h_{1/2})/2 has A(h₀) = 1 > 9√3/16 = [A(h_{−1/2}) + A(h_{1/2})]/2 (rechecked). On the critical XXZ cone A is instead concave and superadditive, so a reduction to extreme bonds by convexity is unavailable.
- C103: no general physics-standard proof. Best positive argument: if conserved on-site currents exhaust the infrared and their stiffness is D = (1/π)Σ_a v_a G_a, then for any positive W with tr(WG_a) ≥ c_a the optical sum rule gives A ≤ π tr[W(K − R)], with K the diamagnetic response and R ≥ 0 the regular optical weight, and a Cauchy–Schwarz bound tr(WR) ≥ 2s₀²/s₁ from its moments gives the scale-free A ≤ π[tr(WK) − 2s₀²/s₁]. Discarding R fails at Heisenberg (tr(WK) = 0.591 > M), so the regular weight must be kept, and neutral sectors need not carry on-site currents.
- Round three on A ≤ C M f(d) with any f is in progress (C104 to C107; snapshots committed): C104 has corrected a representation-normalization factor in the stiffness route and notes that arbitrary compact-boson radii defeat the purely continuum inference; C106 reports a possible obstruction to the proposed crossover lemma (a free band crossing whose thermal matching window collapses while its velocities stay finite); C107 has shown finite-L Casimir quantities converge non-uniformly near critical points.

### 9.77 Round three on A ≤ C M f(d) (C104 to C107, 2026-10-29; all checks rerun by Claude and passing)

- C104 (charges and stiffness): proved, conditional on the WZW stiffness dictionary, a crossover-free polynomial bound for full simple WZW factors coupled to conserved on-site semisimple symmetries, A_WZW ≤ (4π/3) d(d² − 1) M (coefficient 2π/3 for simply laced factors), allowing any number of copies with different levels and velocities and with no bound on levels assumed. It also corrected a representation-normalization factor in the stiffness route (D = kv/(2πI_R) with I_R the on-site Dynkin index). Not proved: the charged U(1) case, since charge quantization and vertex locality do not force a lower bound on Luttinger weights; coset stress tensors and neutral sectors need separate accounting.
- C105 (Boulware deficit): no general bound. Proved for XX a nonlocal positive certificate, e₀Σw − E₀(K_w) ≤ (1/2π)Σ(a_{i+1} − a_i)² for w_i = a_i a_{i+1}, giving defect ≤ (1/4π) log N + O(1) for tent weights; proved that every positive decomposition into uniform open intervals loses at least g* = 1/√2 − 2/π per unit of max w, so that method cannot work; and found that the literal unsmoothed lattice Boulware statement needs smoothing before it is equivalent to the bound.
- C106 (counting fast modes): the proposed thermal crossover lemma is false: a d = 4 nearest-neighbour free chain with two c = 1 sectors keeps v/M fixed while the temperature window of CFT entropy matching collapses. This does not disprove the bound; the operator-counting route needs two further lemmas, stated there.
- C107 (compactness): the unit sphere is compact only in the quotient by gauge, which repairs the premise; the finite-L Casimir quantities Q_L are continuous with Lipschitz constant of order L², so their limit is not controlled; the raw limit need not equal the physical A (ferromagnetic XX gives 2t on even and t/2 on odd lengths), so a vacuum-sector prescription is essential. Upper semicontinuity and a finite F(d) remain open; no fixed-d sequence with A/M → ∞ was constructed.

### 9.78 Rounds four and five on the bound (C108 to C113, 2026-10-29; all checks rerun by Claude and passing)

- C108 (U(1) sectors): proved that the chemical potential ranges over an interval of length at most 8M for an integer charge (replacing non-maximal charge states by a local channel costs at most 4M per site), so the curvature measure of e(n) has mass at most 8M and ∫ v(n)² dn ≤ 4σ_q²M²(1/δ₋ + 1/δ₊) over smooth Luttinger branches: an averaged velocity bound. A pointwise bound is not obtained.
- C109 (neutral sectors): Lieb–Robinson with the optimized bond gives v_a ≤ 4eM, so each Ising factor has cv ≤ 2eM, each minimal model cv < 4eM, each Z_N parafermion cv < 8eM; but no bound on the number of such factors in terms of d is proved.
- C110 (vacuum-sector semicontinuity): twist infima have the wrong sign (ferromagnetic XX gives −4t against A = 2t), tail infima over vacuum-compatible lengths approach A from below; no upper semicontinuous extension or finite F(d) obtained.
- C111, C112, C113 (the smooth-modulation lemma of `specs/_smooth_lemma.md`): within CFT the deficit is exactly Def = (A/(24π))∫w′²/w, with κ = 1/(24π) for any sectors, charged or neutral, and any velocities (three independent derivations agree, including the exact XX response). But the lemma as a lattice inequality is false: a gapped spectator of gap 2δ tensored onto XX (d = 8, M = 3/2, A = 1) gives Def/(M Σ(∇w)²/w) ≥ 1/(300δ), unbounded as δ → 0, because massive modes near criticality respond to the modulation with the same local quadratic form; a z = 2 ferromagnetic chain has a positive deficit under arbitrarily slow modulation; and the representative dependence of the deficit gives further counterexamples with A = 0. None of these violates A ≤ C M f(d).
- Lesson: every local response to bond deformations sees nearly critical massive modes at the same order as the conformal anomaly, so a proof must use a genuinely infrared quantity (the ring's 1/L Casimir term, in which massive modes are exponentially suppressed) rather than a local energy functional.

### 9.79 Round six on the bound: ring, counting, partition-function and bootstrap routes (C114 to C117, 2026-10-29; checks rerun by Claude and passing)

- C114 (ring-only): a symmetry-free coherent-reconnection inequality, E_{2L} − 2E_L ≤ q_L/S_L, with q_L = 2[tr(h σ_L⊗σ_L) − e_L] and S_L the sum of interval purities of the ring ground state, improving bare gluing by the factor S_L; the whole trial family is analytically obstructed on XX, and the missing positive-certificate estimate is stated quantitatively.
- C115 (counting fast sectors): the proposed tower inclusion runs the wrong way (fast sectors cost more energy, slow sectors flood fixed O(1/L) windows); two decoupled XX chains with hopping ratio ε have fixed d, M and one fast sector but a divergent number of low-lying ring states as ε → 0, so unfiltered finite-size counting cannot work.
- C116 (ring partition functions): the ratio R_L(β) = log[Z(L,β)²/Z(2L,β)] isolates A, A = (4/π) lim_τ lim_L R_L(τL)/τ, so a crossover-free estimate R_L(β) ≤ K(d) + k M f(d) β/L would prove A ≤ (4k/π) M f(d); only |R_L| ≤ 4βM is proved (loses a factor L). Other natural ratios retain the reciprocal-velocity sum Σc_a/v_a, which slow sectors make diverge.
- C117 (bootstrap): fixed-support translation-invariant marginal relaxations give an energy density independent of the circumference once supports fit, so they cannot resolve the Casimir term; the q⁴ coefficient of the exact energy-density moment obeys |[q⁴]m₁| ≤ 3M³ but receives massive contributions (gapped transverse-field Ising gives gJ²/4 with A = 0); certified small SDP bounds are recorded. Bootstrap constraints must grow with L or act at momentum of order 1/L.
- Status: every route now reduces to one missing crossover-free lattice estimate, cleanest in C116's form R_L(β) ≤ K(d) + k M f(d) β/L or C114's positive-certificate form. No counterexample to A ≤ C M f(d) has been found.

### 9.80 The partition-function ratio estimate (C118, 2026-10-29; scripts rerun by Claude and passing)

- Proved for the XX chain (d = 2, M = 1/2), every even L ≥ 4 and every β: R_L(β) ≤ 2 log 2 + (5π/3)β/L, the target form with K = 2 log 2, k = 10π/3, f = 1 (both Jordan–Wigner parity sectors, via a quadrature error bound), extending to independent XX stacks with K = 2 log d. Proved for every diagonal classical nearest-neighbour bond: R_L(β) ≤ log d (Cauchy–Schwarz on the nonnegative spatial transfer matrix).
- Proved structure: the target at all β forces L(E_{2L} − 2E_L) ≤ kM f(d), so the finite-size ground-energy problem is unavoidable; averaging a reconnection over seam positions leaves a norm of at least 2M (XX, Heisenberg, SU(3)); R_L is neither convex nor concave in β; and an exact information identity R_L = β⟨V⟩_ρ + I(A:B)_ρ + D(ρ_A‖γ_L) + D(ρ_B‖γ_L) holds, whose marginal relative entropies are not bounded by dimension alone. The general quantum estimate remains open; C119 and a ring-aware SDP bootstrap on Typhon are running.

### 9.81 The ratio estimate in free families, and dimension-bounded divergences (C119, C120, 2026-10-29; scripts rerun by Claude and passing)

- C119: proved uniformly for XX, arbitrary independent XX stacks with different velocities, and the C106 two-band family (including its massive comparison), even L: R_L(β) ≤ 3 log d + (10π/3) Mβ/L, the target form with f = 1, by a matrix trace-convexity and quadrature estimate. Exact finite-size XX ratio in the scaling variable y = βv/L: r_XX(y) = 4e^{−π/(2y)} − 6e^{−π/y} + … at small y and πy/4 − 2e^{−πy/4} + … at large y; a slow sector is negligible at fixed β/L but contributes r_XX(s) at β = sL/v.
- C120: max-relative entropy, Petz and sandwiched Rényi divergences of every order 0 < α < 1, and negative log fidelity are not bounded by dimension (a faithful qubit pair makes each arbitrarily large); support min-relative entropy loses the ratio information. Conditional min-entropy does give a dimension bound and an exact identity R = H_min(C|D)_σ − H_min(C|D)_ρ + log[μ(Y)/μ(X)] ≤ 2 log q + log[μ(Y)/μ(X)] for a seam region of dimension q, but the whole unknown slope sits in the normalization ratio μ(Y)/μ(X), which entropy capacity does not bound. The seam-localized channel route (C121) is running.

### 9.82 The seam-localized quantum-to-classical channel (C121, 2026-10-29; checks rerun by Claude and passing)

- Proved: N seam pinchings in positive imaginary-time evolution obey d^{−N}Z_L ≤ Q_L ≤ Z_L at arbitrary spacing (operator pinching inequality and Schatten–Hölder after a unitary-average expansion), with no gap or Markov assumption; with N = max(1, ⌈βM/(aL)⌉) the log-trace cost is at most (1 + βM/(aL)) log d per seam. Both the coherent and the positive recorded-history Cauchy–Schwarz steps hold, giving the exact reduction R_L ≤ 3N log d + G with G = log[Σ_s q_s²/Q_D].
- The zero-cost reconnection step G ≤ 0 is false: at one insertion with τ = L/M, G = 0.763 (XX, L = 4), 1.888 (Heisenberg, L = 4), 2.293 (SU(3), L = 3), and at L = 8192, β = 4L on XX log[Σq²/Z_{2L}] = 0.987, because a coarse record leaves virtual couplings between the two sides of the seam (cross terms absent classically; stoquasticity does not remove them). Continuous pinching instead converges to Tr exp[−βD(H)] and on XX shifts the ground energy by 1/2 + (1/2)tan(π/(2L)) → M, so dense measurement costs order βM.
- Remaining lemma: Σ_s q_{L,s}² ≤ exp[b₀(d) + b₁(d)N] Q_D with N = max(1, ⌈βM/(aL)⌉), which would give A ≤ [4/(πa)][3 log d + b₁(d)]M. The observed G values are of order one per insertion, including at L = 8192, which is what the lemma requires; it is comparable in strength to the target and not yet proved.

### 9.83 Ring-aware SDP bootstrap (2026-10-29/30; `contrib/bootstrap_ring/`, Claude agent, SLURM jobs on Typhon)

- Translation-invariant SDPs over expectation values of operator strings on the ring, with contiguous k-site marginals, local stationarity and positivity, two-block marginals at every separation, and ⟨O_q†O_q⟩ ≥ 0, ⟨O_q†[H,O_q]⟩ ≥ 0 at every q = 2πm/L; the exact ground state satisfies all constraints to 1e-13 in eight representative cases. With one-site momentum modes the bound per site plateaus below e₀ and cannot see the Casimir term (C117's obstruction). Adding the energy-density mode h_q at the lowest momenta, with positive two-block marginals, the Casimir term appears: the gap to exact diagonalization is 2 percent of πA/(6L) at L = 12 and 20 percent at L = 20, and the bound's own 1/L coefficient is 1.46 for Heisenberg (exact 1.571), 1.20 for XXZ (1.299), 1.10 for SU(3) (4.19), and about zero for a gapped random bond.
- Obstruction: a residual per-site offset (about 5e-4 to 8e-4 for critical d = 2 chains, 0.037 for SU(3)) remains at fixed support, so the certified A still grows like L²; certifying A would need the offset to fall faster than 1/L², hence supports growing with L. No certified bound on A; the bootstrap 1/L coefficients all lie below π M log₂ d. Lesson: the ring's energy-density mode at q = 2π/L is the right object for a lattice proof.

### 9.84 Calibrated c = 3 for the n = 2 member of the d0 = 4464 chain (2026-10-30; `contrib/verify_C44_n2/`)

- Second-order projected n = 2 chain of the C39/C44 construction (two payload modes plus the length fluid), at a resolvable core point (w = 1.8, a = 3.5, η = 2, v = 0.3, D = 1, g = 0.22, ρ = 0.008; core gap 0.509, four-dimensional code; payload Fermi velocities 0.857 and 1.144 and length velocity 1 in units of the mean payload velocity; payload interaction 0.175; V/j_l = −0.157; length density 3/8). A new two-charge block-sparse DMRG (matching exact diagonalization at L = 4 to 5e-11) gives open-chain interval entropies calibrated against the exact three-free-mode reference: c = 2.98 to 3.01 at L = 32 for χ = 512 to 1536 and 2.99 to 3.00 at L = 48, also with the length–payload coupling multiplied by 1000; the length-frozen control gives 2.003 and 2.002 against exact 2. Finite-entanglement scaling (χ = 120 to 288, correlation lengths 6 to 11 sites) is the weaker estimator: the agent reports a calibrated 3.00 on the χ set common to model and reference; Claude's recomputation from the saved points gives 2.99 for the stress-test model and 3.2 for the scheduled model when its additional saved points are included, so the profile estimator carries the result. The raw reference entropies rise from 1.77 to 2.65 with χ against an exact 2.83, so calibration is essential, but the calibrated values stay within 0.04 of 3 while raw values move by 0.9.
- Together with n = 1 (c = 2, §9.28), the projected d0 = 4464 chain gives c = n + 1 at n = 1, 2. Not tested: C39's flavour-margin readout field, whose effect (a splitting below 4e-5) is invisible at these sizes, so the pair-transfer resonance at equal pocket widths is not probed.

### 9.85 The μ ratio and the reconnection lemma (C122, C123, 2026-10-29/30; all checks rerun by Claude and passing)

- C122: the operator-order version of the locality route is ruled out on XX. The comparison constant c_ord = min{c : Y ≤ cX} = ‖exp(βH₁/2)exp(−βH₀/2)‖² obeys log c_ord ≥ (2β/L)csc(π/L) − 2 log 2 on XX, which is (2/π)τL + O(1) at β = τL, a factor L too large. A volume-independent half-time expansional bound G(b) = exp[6b + 2(e^{2b} − 1)], b = βM, is proved but grows too fast for β of order L; seam smearing balances only at order √(βM log d); and projecting onto low energies leaves the Casimir slope in the vacuum normalization. The remaining lemma must compare the optimized SDP witnesses of the two Gibbs operators, not the operators themselves. Seventy full finite-temperature and fourteen ground-projector SDPs were solved on XX, Heisenberg and SU(3).
- C123: the reconnection lemma G ≤ b₀ + b₁N is proved for XX at every spacing a > 0, uniformly in even L and β: G ≤ 2 log 2 + [2 log 2 + (10π/3)a]N (from C118 and C121), extending to independent XX stacks. Exact multi-pinch XX computations to L = 8192 (156 cases) give G/N bounded and saturating in N: at a = 1, G/N = 0.686, 0.840, 0.877, 0.878 for N = 1, 2, 4, 6; for each fixed a and N, G = R_L − N log 2 + O(L⁻²). Correction: at L = 8192, β = 4L, N = 1 the actual G is 2.373 (C121's 0.987 was a lower bound). Interacting exact diagonalization (Heisenberg, XXZ at Δ = −0.9, −0.5, 0.5 on L = 4, 6; SU(3) on L = 3, 4; 192 cases) gives G/N of order one throughout, with no growth from the smaller to the larger ring except a sign change for SU(3) at these tiny sizes.
- Status: for the bound A ≤ C M f(d), the remaining lemma is the bounded per-interval reconnection cost for interacting chains, which free fermions satisfy and small interacting rings are consistent with; it remains unproved.

### 9.86 The reconnection lemma is equivalent to the target (independent Opus 5.5 attempt, 2026-10-30; `contrib/opus_reconnection/`; two of its scripts rerun by Claude and passing)

- Exact decomposition for every finite chain: G = R_L − 2ε_L + ε_D − H₂(π), with ε_L, ε_D the single- and doubled-ring pinching losses and H₂(π) the Rényi-2 entropy of the record distribution, each bounded by a multiple of N log d (verified to about 1e-15 in 99 interacting exact-diagonalization cases and 54 XX cases to L = 8192). Consequently the lemma G ≤ b₀ + b₁N is equivalent, with explicit constants, to the ratio bound R_L ≤ K + κMβ/L itself: at fixed L and slab length τ, G/N lies within O(log d) of τ(E_{2L} − 2E_L), so the surviving cross terms cost the reconnection energy per interval, and the ring level spacing does not enter. Under the infrared CFT plus clustering, G = R_L + N log[Σp⁴/(Σp²)²] + o(1), so the best b₁ is (πa/4)(A/M) + O(log d): a universal b₁(d) exists exactly when A/M is bounded.
- Also proved: continuous pinching gives exact gluing (G → a nonpositive limit) for every chain, at a pinching loss of order βM; for stoquastic chains the zero-cost inequality is equivalent to a comparison X ≤ Y of cross terms Σ B_xx B_yy against Σ B_xy B_yx (numerically X/Y can exceed one, for example 1.055 for Heisenberg at β = 6, L = 4); for commuting chains (assuming the Bravyi–Vyalyi structure theorem) R_L ≤ log d. Larger exact-diagonalization runs are queued on Typhon (jobs 17530900, 17530901).
- Consequence: the coarse-seam formulation is a reformulation, not a reduction. The open problem is the ratio bound itself, whose zero-temperature content is L(E_{2L} − 2E_L) ≤ C(d)M. The Astra attempt (C124) is still running.

### 9.87 Astra attempt at the reconnection lemma (C124, 2026-10-30; scripts rerun by Claude and passing)

- Independent of the Opus attempt (§9.86). Proved: the lemma holds with b₀ = b₁ = 0 (G ≤ 0 at every spacing, N and L) for every chain with commuting adjacent bonds, including complex bonds with entangled eigenstates, using an on-site basis adapted to the two commuting interaction algebras so that recorded histories form a nonnegative matrix with q_{L,s} = B_ss and q_{2L,s,t} = B_st B_ts; in an arbitrary basis R_L ≤ log r ≤ log d and G ≤ log r + 2N log d, with r the number of central sectors. These chains have A = 0 under the infrared assumption. A sufficient boundary-overlap estimate proposed in C123 cannot hold uniformly even for diagonal chains (explicit three-level counterexample). The lemma remains open for noncommuting stoquastic and general interacting chains, consistent with §9.86's finding that it is equivalent to the ratio bound.
- Launched C125 (spec from Jonah): a crossed-channel c-theorem for the quantum transfer matrix, rewriting R_L as the Rényi-2 entropy of a crossed Gibbs state on the Trotter circle and asking whether positivity, Trotter-circle rotations and a Lieb–Robinson crossed spectrum condition give convexity of Θ(q) = −L E(L) in q = L².

### 9.88 The crossed-channel c-theorem does not close the bound (C125, Astra, 2026-11-06; `codex/work_C125/`; scripts not yet rerun by Claude)

- Proved: for a positive spatial transfer operator T, R_L = S₂(T^L/Tr T^L) and D(L) = E_{2L} − 2E_L is nonnegative and nonincreasing on the allowed length subsequence, but neither monotonicity of L D(L) nor convexity of Θ in L² follows. Positive even-ring constructions (checkerboard Trotterization, two sites per step) exist for XX, XXZ with Δ ≥ 0, the antiferromagnetic Heisenberg chain at any spin and integrable SU(2) Sym^k. The on-site reflection criterion fails for the SU(N) swap with N > 2 (SU(3) has D(4) < 0) and for SU(N) Sym^k with N > 2. Odd rings are excluded in general (XX has D(3) = −1).
- Proved: the global crossed spectrum condition ε − ε₀ ≥ |p̃|/V is false. XX has crossed states reached by a local S^z insertion whose crossed energy grows only like log m at crossed momentum 2πm/β, so the real crossed tilt diverges for XX. A valid low-frequency replacement from C109 is |G(iω, x)| ≤ 2‖A‖‖B‖V e^{−ωx/V}/[ω(V − ω)] for 0 < ω < V = 4eM.
- Proved: a positive spectral countermodel (n compact c = 1 bosons in the crossed channel, shells truncated at rank exp(β/2)) satisfies positivity, an entire tilt, the cone condition and the trivial bounds, yet has A = nv unbounded at fixed d = 2, M = 1. Crossed-channel positivity data therefore cannot bound A. The missing input is a relation that fixes the tilted free energy from its untilted slice, as the continuum rotation identity does, or a direct scale comparison of the finite-size anomaly. Conditional closure: Θ convex for L ≥ L₀ gives A ≤ c_A(L₀) ≤ (12/π)L₀²M.
- Numerical: XXZ at Δ = −0.9, outside the positive class, violates R monotonicity and Θ convexity; Heisenberg, XXZ at Δ = −0.5, 0, 0.5, 2 and SU(3) on vacuum-compatible lengths pass discrete Θ convexity at small L; the Sym² Bethe families pass c_A monotonicity (N = 8 gives c_A(8)/M = 7.125 against A/M = 6.597).
- Consequence: lesson for the next round, recorded in `specs/_Cd_target.md`, that a proof of a finite C(d) = sup L(E_{2L} − 2E_L)/M must use fixed-d locality beyond positivity and crossed-channel structure.

### 9.89 Three Astra routes to a finite C(d) (C126 to C128, Astra, 2026-11-06; `codex/work_C126/` to `work_C128/`; scripts not yet rerun by Claude)

- C126 (permutation twist). Proved: E_{NL} − N E_L is the energy of a Z_N permutation twist in N decoupled copies, so N = 2 is the doubling defect and N → ∞ the Casimir variant. Quadratic chains give A ≤ πM through the spread U(N) copy twist. The on-site path from the identity to the permutation is obstructed for interacting chains: the Heisenberg stabilizer is disconnected, fractional swaps cost at least 1/3, and the large-N Ising cat costs of order JN. The remaining lemma is an integrated curvature condition equivalent to the target.
- C127 (crossed majorization). Proved: a scalar bosonic majorant of the crossed spectrum is equivalent to the target, the number of directions reached by one insertion is not a spectral multiplicity, and K(d) ≥ log d is forced. For XX, log Z̃ ≤ log 3 + (2 − 4/π)β/L. All 240 scalar checks pass with S = 6, V = 4eM, K = log 3.
- C128 (world-line majorization). World-line Fock domination is open, because paired jumps give a singular measure. For XX, R_L ≤ 2 log 3 + 2.256 Mβ/L, and the logarithmic high-frequency decay is summable.

### 9.90 Two Claude routes: periodization and the crossed Fock space (2026-11-06; `contrib/cd_free_route/`, `contrib/cd_crossed_fock/`, Claude agents; checks rerun by Claude and passing)

- Periodization (cd_free_route). Proved: variant (a) at fixed d is equivalent to an inclusion of marginal sets, Ω_L ⊂ Ω_∞ + (K/L²)B; iterating the Petz map of a ring marginal gives a translation-invariant state on Z with e_0 ≤ Tr[hσ₂*]. Lemma P, δ_L = Tr[hσ₂*] − E_L/L ≤ KM/L², would give A ≤ (6/π)KM. Numerical: L²δ_L/M is bounded and decreasing on every test, and tends to πA/3 on critical chains, twice the optimal constant. The claimed trace-norm obstruction applies to the full recovered window only, not to the two-site statement (C132 correction, §9.92).
- Crossed Fock space (cd_crossed_fock). Proved: XX majorization with (V′/M, K) = (4.16, 2.04) or (5.54, 1.44); the route reduces to a hypothesis (H2) at least as strong as the target, which would give A ≤ (d² − 1)(12.95 + 1.91/γ₀)M. Adjoint selection rule: Heisenberg's spin-0 dimer sits below the two-particle Fock levels, so operator inequalities built from exact single-insertion levels cannot prove (H2). Numerical: (H2) holds on eight positive-transfer d = 2 chains with S = 2 and V′/M from 1.2 to 6.9.

### 9.91 The d = 2 program and its cover (C129, C131, Astra, 2026-11-07; `codex/work_C129/`, `work_C131/`, `contrib/rigorous_toolbox_map/`)

- Literature map (Claude agent): no general rigorous theorem gives the 1/L term; the best tools are Bethe ansatz for the U(1) class, constructive RG near free points, and gap stability. Six arXiv identifiers were checked against their abstracts.
- C129. Proved: on the U(1) spin-1/2 class the field drops out sector by sector, the DM term becomes a twist, and a twist costs at most Jπ²/(4L); the optimized M has a closed form (checked by Claude against gauge minimization on 40 points); for J_z = 0, L(E_{2L} − 2E_L) ≤ (8 − 4√2 + π²/4)M (checked by Claude by ED on 40 random points, worst 2.23); zero defect for J_z ≤ −J; elementary dilute and near-ferromagnetic bands. Kozlowski arXiv:1508.05741 gives no finite-size rate at half filling.
- C131. Proved: for each L, Q_L is Lipschitz on the compact gauge quotient, so C(2) < ∞ is equivalent to a finite list of local, length-uniform lemmas; a stratification table records each lemma and its tool. Numerical: multistart maximization at L = 2, 4 finds Heisenberg (Q = 4, 2.79); the best odd value at L = 3 is 2.04. The larger search is Typhon job 17578807, queued.

### 9.92 Lemma P with an optimal recovery channel (C132, Astra, 2026-11-07; `codex/work_C132/`)

- Proved: periodization works with any channel R satisfying R(ρ_{l−1}) = ρ_l, so the optimal channel is a semidefinite program; for the ferromagnetic Heisenberg ground-space mixture Petz is optimal with δ = 1/[2l(l+1)] (K = 2); for every quadratic chain E_L ≥ e_0 L − π²M/(6L). The dual-M norm equals the trace norm on the gauge-annihilated slice, which corrects the cd_free_route obstruction.
- Numerical: the optimal channel coincides with Petz for XX and Heisenberg at L ≤ 6, so it does not remove the factor 2 there; it gains substantially only for TFIM. No d = 2 bond shows a constant growing with L.

### 9.93 The fast-tower adversary (C133, Astra, 2026-11-07; `codex/work_C133/`)

- No fast tower was built and no general no-go was proved. Proved: unbounded A does not require many sectors of velocity of order M, since c_k = 1, v_k = M/k gives A = M H_N; the right missing lemma is a crossover-independent bound Σ_{v_a ≥ u} c_a ≤ B(d)(M/u)^p with p < 1, which gives A ≤ B(d)(4e)^{1−p}M/(1 − p). The c = n + 1 encoder is slow for structural reasons (payload velocities suppressed by D^{−2n−4} and an unbiased path gap of order s/n²), and an exact first-order quadratic compression carries at most ⌊(d² − 1)/2⌋ gapless sectors. Open candidates: hierarchical emergence, many weakly charged bosons, several neutral factors. Numerical: 98 ED ground states at d = 2 to 6 show no growth.

### 9.94 The attractive U(1) class, rigorously (C134, C136, Astra, 2026-11-07; `codex/work_C134/`, `work_C136/`; each refereed by an independent Claude agent, verdict correct as stated)

- C134. Proved: for the half-filled XXZ ring at minimizing flux and Δ ∈ [−1 + δ, 0], F_L ≥ L e_0 − 1200J/(δ³L) and L(F_{2L} − 2F_L) ≤ 3000J/δ³ for every L ≥ 2. The mechanism is the negative scattering kernel (L¹ norm 2γ/π − 1 < 1), which makes the counting function monotone, places the roots at midpoint nodes with O(1/L) cumulative error, and, with an exact dressed-energy quadrature identity and clipping of the tails, gives an O(1/L) total error with no log L. The ground state is identified without assuming completeness of the Bethe ansatz. Also proved: Δ = 1/2 exactly with C = 2π²/9 from the supersymmetric twisted identity E = −3JL/8, which Claude verified by ED on even rings L = 4 to 12.
- C136. Proved: every sector at every filling satisfies |F_L(n) − J L e_bulk| ≤ 10¹⁵ δ^{−25/3} J/L, so L(E_{2L} − 2E_L) ≤ 6·10¹⁵ δ^{−25/3} M on the whole U(1) class with J_z/J ∈ [−1 + δ, 0], every field, DM angle and L. The referee confirmed uniformity in the Fermi boundary and found no hidden log L.
- Open: the ferromagnetic crossover. The estimates start at L of order (1 + Δ)^{−1} while the elementary band reaches only (1 + Δ)^{−1/2}; the missing input is an O(J/L) error once (1 + Δ)L² ≥ A (C138).

### 9.95 The marginal side: translation-invariant relaxation and the window lemma (2026-11-07; `contrib/marginal_toolbox_map/` Claude agent, C135 Astra, C137 Claude Opus agent)

- Literature map: no 1/L² rate is known for any interacting class; the best proved rate comparing ring and infinite-chain marginal sets is O(1/n) (Fawzi, Fawzi, Scalet; Blakaj and Wolf). By SDP duality the locally translation-invariant relaxation at level n equals the best t with h ⊗ 1 + X ⊗ 1 − 1 ⊗ X ≥ t, and summing around a ring removes X, so a rate e_0 − t ≤ KM/n² proves the target.
- C135. Proved: this rate holds with K = π²/2 for every quadratic qubit bond; an explicit XX certificate with K = π, X a sine-weighted sum of bonds (Claude checked its positivity at n = 4 to 10); ring saturation is false, with E_LTI(3) = −5/12 against −3/8 for h = (XX + YY)/4 + (ZI + IZ)/8 (reproduced by Claude). Numerical: restricting X to weighted bonds loses nothing for XX through n = 8 and for Heisenberg through n = 6.
- C137. Proved: the window class is the translated-bond restriction of the relaxation's dual; the fixed sine window is a certificate for every quadratic qubit bond; K = 0 for frustration-free bonds. Numerical: the optimal window is the sine-square window at even n and the sine window at odd n, and on the chains tested it reproduces E_n/n, so there the window lemma is equivalent to the target rather than weaker. Lesson 2 does not bite: on TFIM at g = 1 + c/n to n = 1024 the constant peaks at criticality (0.37) and falls like 1/|c|. Constants n²(e_0 − λ_min)/M at n = 16: Heisenberg 1.66, XX 1.05, XXZ (0.5) 1.37, XXZ (2) 0.78. The remaining directions are ZZ, transverse DM and transverse fields.

### 9.96 The ferromagnetic crossover, reduced (C138, Claude Opus agent, 2026-11-07; `codex/work_C138/`; not yet refereed)

- Proved: no rescaling of rapidities is needed near Δ = −1; dividing by η = π − γ gives a kernel symbol bounded below by 2/π uniformly, and ηL is the crossover parameter. At fixed L the η → 0 Bethe equations become a strictly convex log-gas whose exact sum rule reproduces C129's first-order energy. Theorem A: |F_L(n) − J L e_bulk| ≤ (5π/4)ηJ for every L, n and Δ ∈ (−1, 0], which gives the missing input on ηL ≤ A₀.
- Conditional (Theorem B): the missing input holds at every L with K* ≤ 125, closing the attractive U(1) class with C ≤ 628, given a pointwise density estimate |σ_L − r_Q| ≤ C_D/(ηL) and four continuum bounds. C136's cumulative root error grows near a Fermi boundary as η → 0, but the pointwise density error does not.
- Numerical: ηL·sup|σ_L − r| stays between 0.23 and 0.32 for L = 16 to 8192 and 1 + Δ from 10⁻² to 10⁻⁶ (120 further cases at most 0.317, against C136's bound of 11 to 1061); 0 < LΛ/M₀ ≤ πv/2 in every case.

### 9.97 The repulsive U(1) class, conditionally (C139, Claude Opus agent, 2026-11-07; `codex/work_C139/`; not yet refereed)

- Conditional on hypothesis (I), that a real, ordered, reflection-symmetric Bethe solution has energy F_L(n): for Δ ∈ [0, 1 − δ], every filling and both parities, |F_L(n) − J L e_bulk| ≤ (1086/δ) J/L, so L(E_{2L} − 2E_L) ≤ (5400/δ) M on the repulsive U(1) class with every field, DM angle and L. Hypothesis (I) is Yang–Yang's result at even L (as restated by Kozlowski) and is unproved at odd L. Yang–Yang convexity fails for Δ > 0; coercivity survives.
- Proved: all 21 uses of the kernel sign in C134 and C136 have replacements, built on the hole kernel (I + K)^{−1}K, which is positive with norm below 1/2. The continuum side is unconditional and uniform in Δ (r_0 ≤ r_Q ≤ (1 + 4/π) r_0, restricted inverse of norm at most 3). The two compact theorems join at Δ = 0.
- Heisenberg, conditional on (M): if the counting density σ_L = a₁ − (1/L)Σ_k K(x − x_k) is nonnegative on the Fermi interval for the ground-state roots, every constant is uniform and L(E_{2L} − 2E_L) ≤ (2215 + π²/4) M on J_z/J ∈ [0, 1], with Δ = 1 by continuity; (M) also gives (I) at both parities. Numerical: (M) holds in every case computed (Δ up to 1, L up to 4096, all fillings), thinnest at the outermost root; the half-filled defect decreases in L toward πv/2, with maximum 4 at L = 2, Δ = 1.

### 9.98 A physics-level reduction of C(d) < ∞ (C141, Claude Opus agent, 2026-11-07; `codex/work_C141/`)

- Derived (physics standard): by compactness of the normalized bond space, A/M is bounded if and only if it is locally bounded, so a hidden tower must accumulate at one limit bond h*. Bonds with an ordinary CFT or gapped infrared are harmless whatever their crossover length, given perturbative RG and a c-theorem for weak perturbations; C106 is of this kind. Near bonds with a gapped degenerate manifold, Schrieffer–Wolff and Lieb–Robinson give emergent velocities v ≤ C|δ|M, and emergence from z > 1 bands is also slow, which explains why the c = n + 1 sectors are slow. A late fast sector therefore needs a z = 1, multi-velocity, non-Lorentz-invariant flow that increases the central charge above a velocity threshold; none was excluded.
- d = 2: the singular qubit bonds found are the Fibonacci classical point, the ferromagnetic frustration-free points and the quadratic Lifshitz points, each with a single emergent species of vanishing velocity (ED at the Fibonacci point confirms v/δ approaching its projected value). C(2) < ∞ in asymptotic form follows given a finite infrared at each bond, the weak c-theorem and completeness of that list; non-coexistence of sectors would add c_tot ≤ 1 and A ≤ 4eM.
- Budgets: the ground energy and the energy-current f-sum both fail (a sector enters with weight c v/ξ² or c v³/ξ², and only dimension-one conserved charges give crossover-free sum rules); the entanglement rate gives p = 1 at best. The partial positive is operator rank at second order, A_em ≤ ν η² n a [1 + (4/π²) log(4b/a)] with ν ≤ d² − 1 port operators, bounded for gapped cells and growing at most like log n for critical ones.
- Remaining: a c-theorem resolved by velocity, and the emergence budget at third order and beyond or for interacting payloads.

### 9.99 The repulsive counting density (C140, Claude Opus agent, 2026-11-07; `codex/work_C140/`; C139 and C140 refereed together by an independent Claude agent, verdict correct after minor repairs)

- Proved: for 0 ≤ Δ ≤ 63/64, every L ≥ 2 and 1 ≤ n ≤ L/2, the counting density satisfies σ_L ≥ θ r_0 between the outermost roots, with θ = 1 − 8(2 − ln 2)(π − γ)/π³ > 0, and its negative mass beyond them is at most 1/L. The method writes σ_L through the positive unimodal hole kernel, pairs each gap between roots with a neighbouring root, bounds the outer correction by positive Neumann series, locates the outermost root from the exact tail mass, and closes by continuity in γ from the free point. This proves C139's identification hypothesis at both parities on this range, so C139's repulsive theorem is unconditional there, with L(E_{2L} − 2E_L) ≤ 345600 M on 0 ≤ J_z/J ≤ 63/64.
- Remaining repulsive estimate (L*): L r_0(x_n) ≥ ((2 − ln 2)/π²) b ≈ 0.132 b at the outermost root; the proved bound tends to 0.125 b at Heisenberg, about 6% short. Numerical: L r_0(x_n)/b ≥ 0.442 in every case (L from 4 to 1024, Δ up to 0.9999), and σ_L ≥ 0.952 r_0 on the whole line for Δ from 0.01 to 1 and L to 4096. Given (L*), the repulsive class including Heisenberg has C ≈ 19938 (C139's constant rescaled, not rederived line by line).
- Not proved: the ferromagnetic estimate (D) and C138's continuum bounds; for Δ < 0 the hole kernel is negative and the pairing argument does not carry over.
- Referee (C139 with C140): the repulsive theorem holds, L(E_{2L} − 2E_L) ≤ 345600 M on 0 ≤ J_z/J ≤ 63/64 (the sharp supremum of the constant on that range is 325954). The ground-state identification works at both parities through Perron–Frobenius at flux Lπ and analyticity in γ, without Yang–Yang convexity. The 63/64 cutoff is tight but valid. ED at L ≤ 12 matches Bethe to 4.8e-14, and with field and DM at L ≤ 6 the ratio L(E_{2L} − 2E_L)/M peaks at 3.97. The sharper constant 19935 + π²/4 is only a sketch and must not be cited as proved.

### 9.100 Heisenberg included, and (D) up to a logarithm (C142, Claude Opus agent, 2026-11-07; `codex/work_C142/`; not yet refereed)

- Proved: the 6% shortfall in (L*) came from one tail integral in C140 Lemma 3.3. The sharp form is ∫_Q^∞ r_0 ≤ (π/(2b)) r_0(Q), because atan(t)(1 + t²)/t is increasing on (0, 1] with value π/2 at t = 1 (C140 used the constant 2). This gives L r_0(x_n) ≥ b/(2(π − γ)) ≥ 0.159 b against the needed 0.1324 b, and the continuation runs on all of 0 ≤ Δ < 1 with margin θ₁ ≥ 0.168. I checked the lemma and the constant chain by hand.
- Proved: C139 Section 8.2 redone line by line gives |F_L(n) − L e_bulk| ≤ 815/L for L ≥ 29, and L(E_{2L} − 2E_L) ≤ (5·815 + π²/4) M ≤ 4078 M on 0 ≤ J_z/J ≤ 1, every field, DM angle and length, Heisenberg included by continuity at fixed L. This replaces 345600 on [0, 63/64] and the unchecked 19935 + π²/4. The exact γ-dependent constant on a grid is 4059.
- Numerical: 160 Bethe cases pass with large margins (L|F − L e_bulk| ≤ 0.91 against 815), and Bethe matches ED to 2.3e-14 at L ≤ 12.
- Proved: near Δ = −1, ηL |σ_L − r_Q| ≤ (4/π)(log(1/η) + π/2) at every point of the occupied interval and every filling, with a matching bound for (D1). This is (D) up to a logarithm, and improves C136's 1/(η²L).
- Open: the logarithm-free (D) on the shell C₀ ≤ ηL r ≤ C log(1/η) (true value at most 0.317, proved 9.97 at 1 + Δ = 1e-6), and C138's (W1)–(W3). Given (W1)–(W3) alone, the strip constant is about 1235 log(1/η), length-uniform at each η but divergent as Δ → −1.

### 9.101 The list of singular qubit bonds is incomplete (C143, Claude Opus agent, 2026-11-07; `codex/work_C143/`; not yet refereed)

- New singular point: the XZ bond h = X⊗Z, which up to on-site rotations is a single point of B_2 (Z⊗X is its mirror). After a Hadamard on odd sites it becomes the isotropic quantum compass chain. Its spectrum is a critical Majorana chain (c = 1/2) times an exactly degenerate space of dimension 2^(L/2−1), so the degenerate manifold is not gapped from a conformal part, and C141's emergence lemma (5.3) and non-coexistence claim (6.2) do not apply as written. I checked this by ED at L = 8, 10, 12: the ground degeneracy is 8, 16, 32 = 2^(L/2−1), the lowest gap times L is 3.18, 3.17, 3.16, and the Casimir term gives c v ≈ 2.03, consistent with c = 1/2 at v = 4 (a Majorana chain of L Majoranas, i.e. critical Ising on L/2 sites).
- Corrections to C141's list: the frustration-free stratum is a two-parameter family (Bravyi–Gosset), including XYZ bonds with J_x = −J_y at the Kurmann–Thomas–Müller field, which lie in a staggered-U(1) class outside C129–C140. Lifshitz-type band-onset transitions are codimension-one surfaces in every symmetry class, not isolated points, and include the XXZ lower critical field at Δ > 1.
- Proved: a local conserved operator of any range forces rank ≤ 1 for the bond's coupling matrix, and at range two such operators exist only at classical bonds and at XZ (range 3 checked numerically on the rank-one family). The amended list is complete on classical bonds, on the free-fermion locus (z ≤ 2 except a one-sided z = 3 curve, no quartic onsets), and on all frustration-free bonds. In the uniform U(1) class it is complete given the standard Bethe-ansatz completeness.
- Numerical: a screen of 3,780 random bonds with refinement found nothing outside the amended list. One refined minimizer converged onto XZ, another sits next to a band-onset surface. Typhon job 17579868 (2,880-bond screen done, 48 refinements running) writes to /scratch/lustre/jkudlerflam/c143/v3_refined.jsonl.
- Near XZ: XZ + αZX + βYY is a self-dual Ashkin–Teller chain. Off a tuned cone the emergent theory has c = 1 and velocity O(δ). On the cone α ≈ (2/π)β the emergent Majorana chain is strongly interacting and its central charge is not settled (ED to L = 16 gives estimates falling from about 2.5 to 1.1–1.5, not converged). This is the one candidate for a qubit bond with c > 1.
- Conditional derivation of C(2)_asymptotic < ∞, with assumptions (IR) finitely many conformal sectors, (WCT) the weak-perturbation c-theorem at regular bonds, (COMP) completeness of the amended list (proved on the classes above, heuristic in strata without symmetry), (EM-G), (EM-Z) and (EM-XZ) the emergent budgets at Fibonacci, at onset and frustration-free points, and at XZ (the tuned cone is open), and (NEST) finite nesting depth. C141's sharper c_tot ≤ 1 additionally needs c ≤ 1/2 on the tuned cone.

### 9.102 The XZ point: Q_L bounded nearby, no c > 1 on the tuned cone (C144, Claude Opus agent, 2026-11-07; `codex/work_C144/`; not yet refereed; large-L Typhon runs pending)

- Exact: at XZ on even rings, Q_L = L(E_{2L} − 2E_L)/M = 2L tan(π/(4L))/cos(π/(2L)), falling from 2.343 at L = 2 to π/2. I checked this by ED at L = 4, 6, 8 (agreement to six digits, M = 1). On odd rings the compass mapping fails and the defect is negative (about −1.6 at L = 3, 5, 7).
- (Q1), derived at physics level: near XZ, uniformly in the direction of approach, Q_L tends to (π/2)(1 + O(δ)) plus at most O(δ) from the B sector, with the cover lemma at small L. Assumptions: the A/B separation of scales, no Casimir overshoot within a sector (true for massive free Majoranas), and c_B bounded (not small). Numerical: 86 ED runs at L ≤ 10 (2L ≤ 20) over all directions and the cone give a maximum of 2.3476 at L = 2 and 1.847 over L ≥ 4.
- Tuned locus: the first-order B dynamics vanish on the hyperplane dJ_zx + ⟨a⟩ dJ_yy = 0 (⟨a⟩ → 2/π), 7-dimensional in the 8-dimensional normal space. XX + ZZ and the y-field stay free (solved exactly), gap A, and leave a quadratic B touching. The x- and z-fields gap everything. XY and YZ are rotations.
- (Q2): on the cone the second-order B chain has no hopping and is purely quartic, with kernel tail (−1)^r/(2πr). Full-chain DMRG (own NumPy code, checked against ED and free fermions) at L = 128 and β = 0.15, 0.3, 0.5 shows the first-order cone sitting inside a gapped wedge of width about β² log(1/β), with A dimerized and translation spontaneously broken. Inside the wedge one critical line at α/β ≈ 0.64–0.68 has A gapped and one critical Majorana (b–b correlations decaying as about r⁻², entropy fits 0.46 and 0.44 at β = 0.5), so c_tot = 1/2 there. Fits at β = 0.3 are unreliable (0.3 to 1.0). Outside the wedge c ≈ 1, and no fit anywhere exceeds 1.09.
- Not converged: the YX and ZY directions gap A and leave a slow sector that looks like c ≈ 1 at L = 64. Jobs 17580192 and 17580432 (L = 128 to 256) are pending, with outputs in /scratch/lustre/jkudlerflam/c144/.
- Consequence for C141 and C143: (EM-XZ) holds with v_B = O(δ² log(1/δ)) on the tuned locus. The non-coexistence condition should be restated as c_tot ≤ 1 rather than c_B ≤ 1/2, and it holds on the Ashkin–Teller plane, pending YX/ZY.

### 9.103 The ferromagnetic strip reduced to three inequalities (C145, Claude Opus agent, 2026-11-07; `codex/work_C145/`; not yet refereed)

- The logarithm in (D) is not removed. Inside C138 Theorem B the core term is already log-free (η·Var(r/σ) ≤ 13.19 there, proved). The log enters only through the first-order cost of the shell 2 ≤ ηL r < 4C_D(η), which holds about log(1/η)/η roots, so it cannot be counted away. Routes (a) and (b) of the spec both fail at one derivative term, leaving a residual C_D1 log(1/η)/λ².
- Remaining hypothesis (G), on the shell: σ_L ≥ (3/4) r_Q and η·Var_S(r_Q/σ_L) ≤ C_GS. Numerical: in 36 cases up to L = 16384 and 1 + Δ = 1e-6, min σ/r = 0.995 and η·Var ≤ 0.0018, against 0.75 and C_GS = 1.
- Proved: (W2) uniformly in η and Q by a reflection supersolution, 0 ≤ r_0 − r_Q ≤ r_0(2Q − x) + r_0(2Q + x), so r_Q ≥ (1 − 4e^(−2d)) r_0 at distance d from the Fermi point (κ = 0.459 at d = 1). No Wiener–Hopf factorization is needed. I checked the argument line by line.
- Reduced, not proved: (W1) follows from monotonicity of v_Q on [0, Q] (then C_W = 4π; holds numerically in 55 of 56 cases, the exception a 3.4e-6 dip in an underresolved case). (W3) follows from η Z_Q(Q)² ≥ c_Z tanh(2Q) (then C_e = π/c_Z; numerical minimum c_Z = 0.808, constant uses 3/4).
- Refereed: C138 Theorem B is correct as stated, every line and the arithmetic (124.87 ≤ 125, 627.5 ≤ 628) re-derived. With the proved κ its constant becomes 647.4.
- Conditional on (G) and the two (W) reductions, the whole attractive side −1 < J_z/J ≤ 0 has |F_L − L e_bulk| ≤ 136.61/L and L(E_{2L} − 2E_L) ≤ 686 M (661.6 with sampled values), so δ₀ = 1, and the U(1) class on [−1, 1] has the single constant 4078.

### 9.104 Referee of C142 and C145 (C147, independent Claude Opus agent, 2026-11-07; `codex/work_C147/REFEREE.md`)

- C142 Part I: correct after minor repairs, none touching a constant. The margin θ₁ ≥ 0.16803 (ratio 1.202 to the needed value), the majorants 46.96 + 765.50 + 2.47 = 814.93 ≤ 815, the grid supremum 811.30, 5·815 + π²/4 = 4077.47 and the seam bound 4·28 = 112 are reproduced. Repairs: the γ* wording in the continuation (the step shows γ* = γ₀), one written-out sentence for the Q = ∞ case of Lemma I.5, and the location of the 43.42 maximum (near γ = 0.12). Clause (i) of the continuation at n = 1 is not covered by the C140 citation but is unused downstream.
- C142 Part II (Theorem II.6, log-bounded (D) and (D1)): correct, with 3.61 to read 3.6113.
- C145: the reduction of the strip to (G), (W1'), (W3') is correct, the reflection proof of (W2) is correct, Lemma B.2 and Proposition B.3 are correct, and C138 Theorem B is correct (124.87, 627.47, 647.41, 685.5 reproduced). The analytic C_G1 is 13.196, so write C_G1 ≤ 13.20; the strip constant becomes 685.52 and 686 stands. Part D.1 obstruction computations should be tagged heuristic.
- Numerics rerun with independent code: Bethe matches ED to 5e-13 at L ≤ 12 on both sides, flux π per bond is the sector minimizer (twisted ED, L ≤ 10), the ED doubling defect peaks at 4.000 at Heisenberg L = 2, all Part I inequalities hold in 128 repulsive cases to L = 256, and (D), (D1) hold up to the Fermi point in 72 attractive cases at no more than 0.061 of the proved bound.
- Overall: L(E_{2L} − 2E_L) ≤ 4078 M on 0 ≤ J_z/J ≤ 1 is refereed correct. The C129/C136 reductions (field, DM, flux, M ≥ M_0) were read, not re-refereed, and rest on their earlier referee. The strip theorem is correct as a conditional statement on (G), (W1'), (W3').

### 9.105 Two of the three strip inequalities proved (C146, Claude Opus agent, 2026-11-07; `codex/work_C146/`; not yet refereed)

- Proved, (W3'): η Z_Q(Q)² ≥ πQ/(4(1 + 16Q/(3π²))) for every η and Q, so c_Z = 0.373 and e_bulk'' ≤ 8.43 η². The argument has three steps. W(Q) = ∫ Z_Q satisfies W'(Q) = 2 Z_Q(Q)² (a Hadamard-type identity, from symmetry of the resolvent). Z_Q(Q) is nondecreasing in Q (shift the larger interval so the right edges coincide, then use positivity of the resolvent), so W is convex and Z_Q(Q)² ≥ W(Q)/(2Q). A semicircle trial function in the variational principle for W gives the lower bound. I checked all three steps. No Wiener–Hopf factorization is needed.
- Proved, replacing (W1'): Theorem A.2 needs only Var(v_Q) bounded, and Var(v_Q) ≤ 92.9 for every η and Q (59.2 as η → 0), from an exact first-order identity for v_Q' whose edge terms integrate to Z(Q) − 1. Monotonicity of v_Q itself is still unproved and no longer needed.
- Partly proved, (G): its derivative half follows from its comparability half at the cost of log log(1/η). What remains is one one-sided inequality, σ_L ≥ (3/4) r_Q on the shell 2 ≤ ηL r_Q < 4C_D(η). Numerical margin: min σ/r = 0.995 against 0.75 (C145, 36 cases, L to 16384, 1 + Δ to 1e-6).
- Strip constant: conditional on the comparability inequality alone, C_strip runs from 838 at η = 0.1 to 867 at η = 1e-12, length-uniform at each η but growing like log log(1/η). Given full (G) with C_GS = 1, C_strip = 773. The constant on the whole U(1) class [−1, 1] would stay 4078.

### 9.106 The easy-axis side, and the U(1) class assembled (C148, Claude Opus agent, 2026-11-07; `codex/work_C148/`; not yet refereed)

- Proved: for every J_z/J > 1, every field, DM angle and L ≥ 2, L(E_{2L} − 2E_L) ≤ 10781 M. The sector theorem |F_L(n) − L e_bulk(n/L)| ≤ K/L holds with K ≤ 4311 M_0 (M_0 = J(1 + Δ)/4), matched doubling gives (5/2)K ≤ 10778 M_0, and the C129 twist adds π²/4. The exact λ-dependent constant on a grid is 10386, worst as Δ → 1+. K/M_0 ≤ 462 for Δ ≥ 2.08 and tends to zero as Δ → ∞.
- Method: on Δ = cosh λ the Bethe roots live on a circle, and the C139/C140/C142 machinery transfers with constants uniform in the Fermi point, half filling and the dilute regime near h_c1 included. Positivity of the counting density uses a circle version of C140's pairing argument and an elliptic-function bound on the outermost root, L r_0(x_n) ≥ 4 sin(π/8) k K_e/π², which give σ_L ≥ (1 − Φ) r_0(x_n) with Φ ≤ 0.804 for λ ≤ 1 (analytic) and Φ ≤ 0.939 on 1 ≤ λ ≤ 1.36 (cell bracketing of explicit monotone functions). For λ ≥ 1.36, σ_L > 0 for every configuration since min a_1 > (1/2) max K. I redid the cell bracketing in interval arithmetic (pos_interval.py, mpmath iv at 30 digits) after checking the monotonicity direction of each factor: worst cell 0.93895 < 1, direct margin 0.0210 > 0. So the step is rigorous.
- Half filling: the odd-ring spinon gives a positive defect decaying in L, and on even rings the defect is exponentially small beyond the correlation length.
- Glue at Δ = 1: C142's Bethe branch continues analytically through Δ = 1 in XXX rapidities, which identifies the ground state on Δ > 1 at both parities, and both sector bounds are uniform up to Δ = 1 from their sides with the same C129 reductions.
- Numerical: the largest computed defect is 3.9993 M_0 (ED at L = 2, Δ = 1.001, against Heisenberg's 4), and the largest sector error is L|F − L e_bulk| = 1.82 M_0 against the proved 4311 M_0.

The U(1) class, assembled. For the bond h = J_⊥(SˣSˣ + SʸSʸ) + D(SˣSʸ − SʸSˣ) + J_z SᶻSᶻ − (B/2)(Sᶻ⊗1 + 1⊗Sᶻ), with J = √(J_⊥² + D²) and Δ = J_z/J, at every field B, DM angle and length L ≥ 2:

| Δ = J_z/J | bound on L(E_{2L} − 2E_L)/M | source | status |
|---|---|---|---|
| Δ ≤ −1 | 0 | C129 | proved |
| −1 < Δ < −1 + δ | 773, given σ_L ≥ (3/4) r_Q on the shell | C145, C146 | conditional |
| −1 + δ ≤ Δ ≤ 0 | 6·10¹⁵ δ^(−25/3) | C134, C136 | refereed |
| 0 ≤ Δ ≤ 1 | 4078 | C139, C140, C142 | refereed (C140 referee, C147) |
| Δ > 1 | 10781 | C148 | proved, not refereed |

So the bound is proved on all of the class except a strip −1 < Δ < −1 + δ next to the ferromagnet, for every δ > 0, and given one inequality on the shell (numerical margin 0.995 against 0.75) the whole class has the single constant C = 10781.

### 9.107 The strip by comparing L with 2L (C150, Claude Opus agent, 2026-11-07; `codex/work_C150/`; not yet refereed)

- Proved, unconditional on the strip (η₀ = 1/2, every L ≥ 2, 0 ≤ n ≤ L/2): L|F_{2L}(2n) − 2F_L(n)| ≤ (85.8 log(1/η) + 180) J, and |D| ≤ π sin(η) n/L for every Δ ∈ (−1, 0]. This is the best unconditional strip bound so far, still logarithmic.
- Proved, conditional on (H_I): ∫ over S′ = {2 ≤ ηL r_Q < 2C_D(η)} of dy/(ηL σ_L) ≤ C_I, the bound is L·D ≤ (112.2 + 4.87 C_I) J. Any lower bound σ_L ≥ c_* r_Q on S′ with c_* > 0 implies (H_I) with C_I = 3.18/c_*, so C145's comparability with c_* = 3/4 gives 133. No derivative or variation bound is needed, and there is no log log.
- Decomposition: with the quantum numbers matched exactly (t′_{2j∓1} = t_j ∓ 1/(4L)), the 2L trial is the refined grid z_L⁻¹(t′) on a central block and doubled roots outside. D = A + B with A a Richardson difference of midpoint rules (zero off the block) and B the self-consistency correction.
- Key lemma, proved (checked by me): for any trial configuration with residual ρ in counting units, |E(true) − E(trial)| ≤ 2a Σ|ρ_k|, a = π sin η. Along the homotopy F(y(s)) = (1 − s)ρ, the Gaudin Jacobian is a symmetric M-matrix for every configuration with J·1 = 2πL a₁, so a₁ᵀJ⁻¹ = 1ᵀ/(2πL), and |e′| ≤ 2a a₁ cancels the 1/η amplification exactly.
- Checks: 216 cases (1 + Δ from 0.1 to 1e-6, L to 2048, eight fillings). Δ = 0 is exact, the CFT ratio is 0.99997 at 1 + Δ = 1e-2, L = 2048, and the η → 0 log-gas limit agrees to 1e-3. Numerical margin for (H_I): λ_L/λ_Q ≥ 0.9993 on S′ (any positive c_* suffices), I₁ ≤ 0.97 against the 6.45 used.

### 9.108 The strip by Destri–de Vega (C151, Claude Opus agent, 2026-11-07; `codex/work_C151/`; not yet refereed)

- Proved, conditional on σ_L ≥ c r_Q on the shell S = {M₀ ≤ λ < 4C_D(η)} for some fixed c > 0 (nothing else: no (D1), no (W1), no variation of r/σ): for η ∈ (0, π/4], all L ≥ 2 and 0 ≤ n ≤ L/2, |F_L(n) − L e_bulk| ≤ K/L, with K = 76.7 at c = 3/4, M₀ = 2 (so C = 5K + π²/4 = 386), and K = 107, 203, 492 at c = 1/2, 1/4, 1/10. Any positive c closes the strip.
- Proved, unconditional: K₀(η) ≈ 31.4 log(1/η) + 140, a logarithm with about one eighth of C142's slope.
- Key estimate: Re σ_L > 0 on the strip |Im z| < min(η, π/4), which gives an exact contour identity at every filling with vertical sides at the counting boundary. Harnack's inequality gives 2π Im θ(x + iη/2) ≥ 1.954 ηL σ_L(x), so the decay of the nonlinear terms above a point needs only the density at that point. The tail is charged with the sawtooth, and the dense set sums as a geometric series.
- Why the shell bound does not follow from DdV: a bootstrap through the NLIE leaves log(1/η) from the scales between η and 1 (it closes only for λ ≳ log log(1/η)), and the sparse tail reintroduces C145's boundary obstruction.
- Relation to C150: σ ≥ c r implies C150's (H_I) with C_I ≤ 6.07/c, so both routes now need the same single input, a positive lower bound on σ_L/r_Q on the shell (numerically ≥ 0.99, 34 cases, L to 16384, 1 + Δ to 1e-6).
- Checks: the contour identity holds at every height and filling tested, the free point gives −π/6, the Fermi-edge pieces reproduce −πv_F/(6L), and the NLIE solved without Bethe input matches the Bethe roots and energy down to η = 0.01.

### 9.109 The strip by a discrete comparison principle (C149, Claude Opus agent, 2026-11-07; `codex/work_C149/`; not yet refereed)

- Proved, unconditional: the ground-state Bethe map at flux Lπ on −1 < Δ ≤ 0 is an M-function, so F(y) ≥ F(x) implies y ≥ x for all real configurations and any quantum numbers (uniqueness of the real solution follows). Exact residual identity F_j(y)/(2πL) = [Z₀(y_j) − I_j/L] − T*(N₀ − N_y)(y_j) for every configuration. The energy identity taken at the continuum Fermi point Q_ν (m(Q_ν) = ν) has no density-mismatch term, so (W3) is not needed.
- The strip is reduced to (R_κ), an inequality about explicit trial configurations with no Bethe roots in it: the quantiles y⁺ of the continuum counting function Z₀ at the shifted levels I_j/L + κ₀/L satisfy F(y⁺) ≥ 0 (the subsolution follows by reflection).
- Conditional on (R_κ) with κ₀ = 0.4: L sup|N_x − N₀| ≤ κ₀ + 1/2, and with the factorization T = T_A * T_B into positive unimodal kernels, log-free (D) and (D1) at every point, ηL|σ_L − s_ν| ≤ 0.573 and η²L|σ_L' − s_ν'| ≤ 1.459. C138 Theorem B rerun at Q_ν gives K = 136.4 and L(E_{2L} − 2E_L) ≤ 685 M on the whole attractive side (η₀ = π/2), so the U(1) class has C = 10781 overall (4078 on J_z/J ≤ 1).
- Numerical: (R_κ) at κ₀ = 0.4 holds in all 173 cases (1 + Δ from 0.5 to 1e-6, L from 16 to 8192, nine fillings), worst margin 0.088 η/L at half filling against a twist gain 0.255 η/L. On the true roots: cumulative error ≤ 0.758 (bound 0.9), ηL sup|σ − s| ≤ 0.317 (bound 0.573), η²L sup|σ' − s'| ≤ 0.217 (bound 1.459).
- Numerical caveat for earlier sections: the graded Nyström continuum solver used in C138 to C146 is mesh-limited near Δ = −1 (r_Q(0) off by 5 to 10 percent at 1 + Δ = 1e-6, Q = 8, interpolation error amplified by B = π/(2η)). C149 replaced it by an FFT solver (uniform grid h = η/10, order-7 Gregory end weights, m(Q) converged to 1e-12). With the new solver the (D) value 0.317 on true roots is confirmed, which gives σ/r ≥ 1 − 0.317/λ ≥ 0.84 on λ ≥ 2, so the earlier shell-comparability numerics stand with a corrected reference.
- What remains for (R_κ): a sharp sawtooth-quadrature lemma with sparse-limit constant b/(2π) (a crude bound loses a factor 2, the whole window b/4 < κ₀ < 1/2), and control of r_Q on the kernel scale near the Fermi edge.

### 9.110 Literature: no known nearest-neighbour one-site qubit chain with c > 1 (C156, Claude Opus agent, 2026-11-07; `codex/work_C156/`)

- No credible translation-invariant, nearest-neighbour, one-site-cell qubit chain with c > 1 was found, exact or numerical. Chains with a two-site cell whose neighbouring bonds and fields are related by one fixed rotation are in the class after a site-dependent rotation (spectra checked equal to 1e-14), and all of these report c ≤ 1: Kitaev–Γ, Kitaev–Heisenberg–Γ, Kitaev–Γ–DM (c = 1 with emergent SU(2)₁, c = 1/2 near the antiferromagnetic Kitaev point), compass and 90° compass, Kitaev–XX–Γ. Also c ≤ 1: XYZ in fields, the interacting Kitaev chain, XYZ with DM, XXZ with DM and transverse field, Zhang–Song's staggered-field XXZ.
- Near-misses outside the class: Rahmani et al.'s interacting Majorana chain (c = 3/2) needs a next-nearest X_jX_{j+2} term; O'Brien–Fendley and Titvinidze–Japaridze need three-site terms; the J1–J2 chain in a field (c = 2) and Fitouchi–Chepiga (c = 3/2, 2) need next-nearest neighbours; Kitaev zigzag chains and Kitaev–Γ ladders are not one-site chains.
- Open case in the class: the antiferromagnetic Kitaev chain in a [110] field (Sørensen et al., arXiv:2209.06221), which is the XZ point with a uniform field along the bisector of its two axes. It is reported critical for 0 < h < 0.511K with no central charge given, and its low-energy scale is about 1e-3 K. C144 did not test this direction. Typhon job 17594698 (DMRG to L = 128) was pending.
- Rigorous: de Leeuw–Pribytok–Ryan classify integrable nearest-neighbour qubit chains (Hermitian ones are XYZ-, XXZ- or free-fermion type, c ≤ 1); frustration-free gapless chains have z ≥ 2 (Bravyi–Gosset, Gosset–Mozgunov, Masaoka–Soejima–Watanabe); Benfatto–Falco–Mastropietro control perturbations of free fermions; no general bound on c from the local dimension exists.

### 9.111 Effort narrowed to the XZ point (2026-11-07)

On Jonah's decision, the parity-class (C152) and staggered-box (C153) agents were stopped mid-run so the effort concentrates on a possible counterexample at the XZ point. Their partial work is committed as is in `codex/work_C152/` and `codex/work_C153/` (no RESULT.md; PROGRESS.md records what was done). Typhon jobs 17594700 (C153 grid) and 17594696 (C152 random scan) were cancelled; job 17594701 (C152 DMRG, which includes a near-XZ family) was kept. C153 left a VUMPS / finite-entanglement-scaling code (`vumps.py`, `fes.py`) that the XZ work can reuse.

### 9.112 Symmetry-free bonds and the XZ directions (C154 and C158, Claude Opus agents, 2026-11-07; `codex/work_C154/`, `codex/work_C158/`; not yet refereed; Typhon runs pending)

- C154, proved: beyond the uniform symmetric strata there is one new site-dependent-symmetry stratum, αXZ + βYX + γZY at B = 0 (the cyclic stratum). A period-3 cyclic relabelling of the axes on successive sites maps it exactly to a period-3 XYZ chain, and at α = β = γ to the Heisenberg chain. The YX and ZY directions at XZ lie on it, so there the low-energy theory is either emergent SU(2)₁ (c = 1) or dimerized (c = 0). Ring ED to L = 24: c = 1 above ε_c ≈ 0.1 to 0.15 (Casimir 1.08 at ε = 0.2), dimerized below. The full stratum list is N(8), T(5), K(5), parity(4), U(1)(3), staggered(2), cyclic(2), XYZ(2).
- C154, derived: the codimension count does not exclude c > 1 (Ising × floating boson at codimension 1, two floating bosons at codimension 0 are kinematically allowed); every bond is a U(1)-class parent plus harmonics carrying lattice momentum; a one-species lemma gives at most one complex incommensurate soft mode at an instability of a uniform polarized phase. c ≤ 1 on all of B_2 follows conditionally, through a closure lemma, from the weak c-theorem, the controlled classes and strata, and (H): no critical component is an island disconnected from them.
- C158, numerical (VUMPS with finite-entanglement scaling, χ ≤ 56, validated on Ising 0.50 and XX/Heisenberg/XZ+ZX 1.00 to 1.02): on the bisector-field line (Kitaev chain in a [110] field) c = 1.01 ± 0.03 at h = 0.45 and 1.07 ± 0.05 at h = 0.3, with the Hadamard × reflection symmetry unbroken; h = 0.6 is gapped with that symmetry broken. The structure is a gapped k = π sector plus one critical incommensurate compact boson; transfer-matrix ratios 1 : 2.05 : 4.28 : 5.19 : 5.54 at momenta Q, 0, Q, 0, 2Q fit one boson, not Ising × boson. Small-χ readings of c = 2 to 6 are an artifact of the saturating k = π length. YX at 0.2, 0.3 gives c = 1.00 ± 0.05.
- Pending: C158 Typhon job 17595611 (χ to 384, open chains L = 192, 256; estimated start 2026-11-12), C154 jobs 17594697, 17594764, 17594765.

### 9.113 The staggered-charge box (C153, Claude Opus agents, resumed, 2026-11-07; `codex/work_C153/`; not yet refereed; Typhon job 17596442 pending)

- Phase diagram inside the box (b > 0, −1 < u < r = √(1 + b²/4)): a gapped polarized phase (c = 0) for u < u_BKT(b), a BKT line u_BKT = 1/√2 + κb + O(b²) with κ ≈ 0.25, a Luttinger liquid of the staggered charge Q (c = 1) for u_BKT < u < r, and the factorization line u = r (z = 2, L + 1 ground states, not conformal; proved). At large b the liquid strip has width g_c/b³ with g_c ≈ 0.6 to 1, and g = b³(r − u) is the controlling scaling variable (derived; the collapse onto one continuum problem is heuristic).
- c ≤ 1 everywhere in the box: c = 1 carried by one boson for Q, with the neutral total-density mode gapped. Derived near b = 0 and near u = r at fixed b; heuristic at large b (the conversion vertex has dimension well below 2 and locks the neutral mode); numerical in between (c = 0.95 to 1.05 for b ≤ 4, 120 small-ring ED points with no third phase, low-χ checks at b = 6, 8 giving 0.9 to 1.1). No c > 1 candidate. One feature to watch: at b = 6, 8 on rings L ≤ 20 a zero-momentum neutral level sits next to the boson doublet, read as a finite-size effect; open chains at L = 64, 128 in job 17596442 will decide.

### 9.114 The parity class (C152, Claude Opus agents, resumed, 2026-11-07; `codex/work_C152/`; not yet refereed; Typhon job 17597010 pending)

- Gapless loci of the interacting Kitaev chain with hopping phase: a chiral Luttinger phase (c = 1, open and 4-dimensional, containing the U(1) liquid and the D = 0 floating phase); Ising surfaces (c = 1/2: topological/trivial, CDW/disordered, and the attractive-side ferromagnetic pair of Peotta et al.); Pokrovsky–Talapov edges (z = 2) and KT edges (c = 1) of the c = 1 phase; an exact asymmetric self-dual Ashkin–Teller locus at cos φ = 0, h_z = 0 (c = 1), with an exact duality (t, Δ, V) → (Δ, t, −V) checked by ED to 2e-14, so crossing Ising surfaces give one boson rather than Ising × boson; zero-field XYZ critical lines (c = 1); Ising–chiral Lifshitz junctions (c ≤ 1); possible tricritical Ising points (7/10, not located); and the XZ window, controlled by C144, C154 and C158.
- Standard of c ≤ 1: proved on the free locus and for every period ≤ 2 mean-field state; exact on the U(1) and zero-field XYZ loci; derived near these, on and near the Ashkin–Teller plane, at the junctions away from XZ and at strong repulsion; heuristic at strong attraction with Δ ≠ 0; numerical only at intermediate couplings far from all controlled regions.
- No c > 1 found. Open readings of 1.13 to 1.55 at L = 64 in the XZ window (φ = 1.45, small h_z, |V| = 0.4 to 0.8) are an interaction effect (the free chain reads 0.99 to 1.01 there) of the same crossover type C158 showed converges to 1; ring-ED spikes of 2 to 6 at φ = 0.9 became c = 1 under DMRG at L = 96. Job 17594701 completed (283 runs, `typhon_out/`); job 17597010 (finite-size scaling to L = 192 at the worst points, attractive transitions, strong attraction) is pending.

### 9.115 The time-reversal strata T and K (C159, Claude Opus agent, 2026-11-07; `codex/work_C159/`; not yet refereed; Typhon job 17597033 pending)

- T (B = 0, only Kramers time reversal): symmetric gapless phases must carry the LSM anomaly. The liquid is "half-floating" (one field locked, zero magnetization along the spiral axis, k_F = π/2), stable for Luttinger K > 1/2, ending by KT into Néel or dimerized order with a c = 1 line between them. No symmetric critical point has c < 1 (derived for Ising, tricritical Ising, Potts; assumed for other minimal models), so c = 1/2 occurs only inside a broken background, as at C144's line near XZ. For c > 1: SU(2)₂ cannot carry the anomaly, SU(2)₃ needs 9 tunings in a 5-dimensional stratum, Ising × boson needs 2 tunings unless K > 15/8, two bosons need none only if K₁ + K₂ > 2. So c > 1 is not excluded by symmetry, but on large loci it can only sit next to attractive liquids (K > 1).
- K: exactly XYZ plus a DM vector along one principal axis plus a perpendicular field, two sheets by the sign of J_xJ_y + D_z². No anomaly; the liquids are also half-floating; the c > 1 counting matches T.
- Standard: c ≤ 1 derived near the U(1) and XYZ strata, the antiferromagnetic Heisenberg point, the cyclic SU(2)₁ region and ferromagnetic onsets, proved on the free locus in T, and supported by the literature (Kitaev–Γ, Kitaev–Heisenberg–Γ mapped into T, Garate–Affleck, tilted-field XXZ, Kitaev chain in a [110] field). Globally conditional on one assumption narrower than C154's one-species statement: no point hosts two independent gapless sectors one of which is an attractive (K > 1) half-floating boson.
- Remaining window for c > 1: the attractive liquids of T and K, especially the easy-plane neighbourhood of the ferromagnetic Heisenberg point. No unexplained critical phase found; one pilot flag (1.48 at rings 12/16) fell to 0.71 at rings 14/18.

### 9.116 Adversarial search for c > 1, local pilot (C155, Claude Opus agent, 2026-11-07; `codex/work_C155/`; main Typhon search job 17594766 pending)

- Estimator: ring ED at L = 12, 16 (about 1.5 s per point), the minimum of a half-ring-entropy and a Casimir estimate. Calibration: Ising 0.50 (Casimir 0.45), c = 1 points 0.75 to 1.05, range-3 controls outside B_2 1.90 (c = 2) and 1.39 (c = 3/2), so it reads low by about 0.1 there; flat-band, ferromagnetic and saturation points read 0. Open-chain DMRG accurate to about ±0.04 at L = 64. A gap-first stage, neighbourhood ball tests, XZ seeds and the Kitaev [110]-field line were added at the coordinator's request.
- Local pilot coverage: 320 uniform points (5.6 percent read gapless), 11 c-maximization restarts, 3 gap-first restarts, 5 XZ directions to L = 128, 6 fields on the Kitaev line. Screen maxima of 4.0 near the Kitaev line and 2.40 in an unpenalized restart were artifacts that collapse at larger rings or under DMRG. Largest refined value: 1.13 at L = 64 for one uniform point, falling from 1.18 at L = 32 (two-size estimate 0.80), read as c = 1. XZ neighbourhood (YX, ZY, YX + ZY): 0.93 to 1.00 at L = 128. Kitaev [110] line: undecided by this method (ring gaps about 1e-3; C158's finite-entanglement scaling gives c = 1.01 ± 0.03 there).
- Pending: job 17594766 (127 XZ and Kitaev points with DMRG to L = 256, 12,000 uniform points, 4 h each of gap-first and c-maximization search, refinement to L = 128), outputs in /scratch/lustre/jkudlerflam/c155/, digest in summary.txt.

### 9.117 The XZ loophole closes at second order (C157, Claude Opus agent, 2026-11-07; `codex/work_C157/`; not yet refereed; Typhon job 17595396 pending)

- Proved: every direction at XZ even under the parity ΠY is exactly a two-chain Majorana model whose only four-fermion term is the YY plaquette. On the tuned hyperplane the effective B chain is λ Σ b_j − 2 Σ_r K(r) b_j b_{j+r} with λ = t/β², K(2) = 0.150 β², K(3) = 0.102 β², tail (−1)^r/(2πr). K is positive semidefinite and peaks at q = π, so the induced interaction is attractive in every channel: matched to Rahmani–Zhu–Franz–Affleck it is g = 2β²K(2) > 0, while their c = 3/2 phase needs −2.86 < g/t < −0.285 (the repulsive sign).
- Derived: no direction gives B more than c = 1/2 at second order. The Lifshitz route needs h3/h1 < −1/3 or > 1, and mean field gives at most 0.94 (code reproduces RZFA's Lifshitz point −0.29); ED of the B chain at small λ is twofold gapped and dimerized; the tricritical route (c_tot = 6/5) is preempted by the log-divergent tail, which dimerizes B at g/t ≈ 0.1; XX + ZZ and B_y gap A and the plaquette cannot add Fermi points, so c_tot ≤ 1.
- Bisector field (Kitaev chain in a [110] field): exact two-qubit-cell rewriting (checked to 1e-14); on L = 48 to 128 (C156's job) the pair estimates are 0.94, 0.98 at h = 0.4 and 0.89, 0.97 at h = 0.45, against 1.24 to 1.69 at L = 128 for the same pipeline in RZFA's c = 3/2 phase. Consistent with C158's c = 1.01 ± 0.03.

### 9.118 The ferromagnetic Heisenberg neighbourhood (C160, Claude Opus agent, 2026-11-08; `codex/work_C160/`; not yet refereed; Typhon job 17599583 pending)

- Blow-up: writing J = −1 + ε(S₀ + A(D)), B = εb, the chain near h_F is a quantum Landau–Lifshitz ferromagnet E[n] = ∫ [|n′|²/2 + εD·(n × n′) + ε(nᵀS₀n + b·n)] with effective ħ = 2ε^{1/2} (2ε in the wedge around pure DM), so the infrared is fixed by classical ground states and their soft modes. Phases: gapped ferromagnets (open), first-order surfaces between canted states, Ising pitchforks (c = 1/2) at their edges, and in wedges of width ε chiral soliton lattices (c = 1, K ~ 1/ħ_eff), pinned phases, Pokrovsky–Talapov transitions and saturation onsets.
- One soft pair at most (derived, conditional on the semiclassical correspondence): with the chain coordinate as time the energy is a Tonelli Lagrangian on S²; Mather's graph theorem plus Poincaré–Bendixson make ground states uniform or singly periodic (no second phason); non-conjugacy of minimizers excludes Hessian zero modes at nonzero Bloch momentum; a Bloch-zero mode pairs canonically with the phason (z = 2). This rules out two bosons, Ising × boson and floating + broken symmetry near h_F in T, K and N alike, removing the window left by C159 Proposition C.
- Proved: a stable polarized state with no DM along the moment has its magnon minima only at k = 0 or π, so the J1–J2 two-minimum mechanism is absent for one-site nearest-neighbour bonds (5916 random references; checked independently by the coordinator on 80,000 random stable BdG forms A = a₀ + a₁cos k, B = b₀ + b₁cos k, no interior minimum). Exact two-magnon bound state: the Bethe 2-string shifted by the DM twist (ED 1e-4 at L = 32); easy-plane side repulsive (one liquid), easy-axis side McGuire collapse (first order). Numerical: on the K-stratum chiral soliton lattice a first-order transition to the canted ferromagnet precedes the tilt instability, so floating + SSB does not occur there.
- Open: the semiclassical correspondence is assumed, the continuum-to-lattice step is heuristic, ε₀ is not quantified. Job 17599583 (DMRG of the soliton-lattice family at ε = 0.5 to L = 256, lattice texture periodicity at N = 2000) is pending in /scratch/lustre/jkudlerflam/c160/.

### 9.119 Referee of C160: classical lemmas sound, scope much narrower (C161, Claude Opus agent, 2026-11-08; `codex/work_C161/`)

- Sound: M1 (Tonelli with exact magnetic term, Mather set, Poincaré–Bendixson on S²), M2 (Bloch zero modes at k ≠ 0 are energy-preserving, so the reduced 2 × 2 Poincaré flow applies; rephrase via conjugate-point-free Lagrangian subspaces; checked on a non-planar DM-coupled tilted chiral soliton lattice, lowest band rises from 0 at k = 0), M3, Lemma 3′ (the case with DM along the moment needs one line, supplied), the two-magnon problem.
- Scope corrections. (i) Proposition C's window is K > 15/8 (Δ_AF < −0.669, ε = 1.11 and ħ_eff ≈ 2.2 on the pure-DM line) or K > 1 (from the XX point), while C160 holds only for ε < ε₀ unquantified, so it covers the K → ∞ corner only. (ii) The singular set through h_F is the curve J = −R(φ), h_F after a site-dependent rotation (open-chain spectra agree to 5e-14, checked independently by the coordinator); C160 covers only φ of order ε₀. (iii) The semiclassical correspondence is not uniform: at onsets the infrared boson has K → 1 even though ħ_eff is small. Magnon onsets are repaired by Lemma 3′. Soliton onsets are not.
- Soliton onsets: for the K-stratum bond v = (−1, −1 + με²/2, −1 + κε²/2, 0, 0, −ε, −hε², 0, 0) with κ = 0.05, μ = 2.3, h ≈ 0.02 to 0.03, the planar 2π soliton is tilt-unstable (lowest out-of-plane eigenvalue −0.067) over a stable Z₂-symmetric background, and minimization gives two degenerate tilted solitons with solid angles 2π(1 ∓ 0.443), hence crystal momenta π ∓ 0.443π (Haldane, P = sΩ). Reproduced by the coordinator. This is the two-minimum mechanism for solitons. Here like tilts attract (decay rate 0.27 against √(h + κ) = 0.28), so the onset is first order into a tilted chiral soliton lattice with c = 1. A miscible c = 2 onset is not excluded, only argued against heuristically.
- C160's statement that floating + SSB does not exist holds only for its μ = 0 family. With μ ≠ 0 the tilted lattice is the classical ground state on an open set, which C159 already counts as needing two tunings for c = 3/2.
- Verdict: C160 does not remove the C159 Proposition C window. No c > 1 bond was found.

### 9.120 The attractive liquids at intermediate K (C163, Claude Opus agent, 2026-11-08; `codex/work_C163/`; not yet refereed; Typhon job 17600999 pending)

- Derived: the zero-field attractive XXZ parent of every T and K_I liquid has only the Gaussian c = 1 tower below the lattice scale (its strings are massless breathers of the same boson), so no mode is close to softening. In the U(1) class every soft wavevector lies in Zα + Zπ. Lemma L: a second boson built from n-quantum composites is locked to the host by the conversion vertex e^{i(ϑ₂ − nϑ₁)}, of dimension n²/4K₁ + 1/4K₂, relevant for attractive hosts unless K₂ < 1/4. Proposition P: to every order in the U(1)-breaking amplitude the infrared is one boson with clock-locked tongues. Lemma U: the uniform Θ-odd tilt that attraction favours is the boson's own density channel (C160 M3), so Ising × boson through it has codimension 2.
- Proposition G (derived, Schrieffer–Wolff asymptotic): near the twisted curve −R(φ) at φ of order one, non-resonant charge-q terms average out and resonant Z_m clock terms gap the boson, so c ≤ 1 there. With C160 at the ends, the whole curve is covered away from soliton-onset layers.
- Numerical (C159's completed job 17597033): 16 attractive liquids with ring K_S(16) = 1.9 to 7.5 read c_prof(128) = 0.98 to 1.03. Two further points in the window are still above 1 but falling, hi_T_37 (K_S ≈ 7, c_prof 1.28, 1.33, 1.08 at L = 32, 64, 128) and hi_K_1 (K_S ≈ 2.3, 1.40, 1.37, 1.26, c_two 0.80). The only point meeting C159's pre-fixed c > 1 criterion is dist_T_50 (c_prof 1.14, 1.18, 1.27, rising), a repulsive liquid (K_S ≈ 0.65) within 0.2 of XZ, outside the attractive window but in the XZ crossover region; hi_T_40 (K_S ≈ 0.45) also rises (1.02, 1.06, 1.15). If either saturated above 1.2 it would contradict C159 Proposition C and would first need a hidden-symmetry check. Job 17600999 runs these four to L = 256 with FES, 16 twisted-curve bonds and 40 gap-first paths through the window.
- Remaining assumption (A_W): in an attractive liquid of T or K_I with U(1)-breaking amplitude comparable to the exchange, no excitation outside the boson's charge lattice softens while the liquid survives, i.e. (a) no independent staggered Ising mode for K > 15/8 (heuristic energetics: no frustration in a nearest-neighbour bond) and (b) no composite band minimum outside Zθ + Zπ for K > 1 (the soliton two-valley band of C161 is the known instance, C162). N stays under C154's (A6).

### 9.121 Two-valley soliton onsets are immiscible where protected (C162, Claude Opus agent, 2026-11-08; `codex/work_C162/`; not yet refereed)

- Effective theory [derived]: the dilute hard-core soliton gas with two valleys at π ± δ is a pseudospin XXZ chain, J_z = c(a₊₊ + a₋₋ − a₊₋ − a₋₊) with c = π²ħ²n³/3M from 1D scattering-length differences (a power of ε), J_⊥ from valley exchange (exponentially small in 1/ħ_eff). Valley number is conserved except on measure-zero density lines, since a flip carries the incommensurate momentum 2δ. c = 2 needs |J_z| < |J_⊥|.
- Rearrangement lemma [proved, continuum; lattice for J_z ≤ 0]: in K the DM density (n × n′)_z does not involve n_z and W is even in n_z, so n_z → |n_z| preserves the energy. Minimizers have single-signed n_z and V₊₊(R) < V₊₋(R) at every separation, a classical form of Eisenberg–Lieb pushed past the SU(2) point. Checked by the coordinator. Scattering lengths are monotone in the potential, a single exponential channel gives additive hard rods, and the subleading channel gives J_z < 0 in K and in U (site inversion times a π rotation about the background spin, the only other stratum with a protected two-valley band) [derived]. In 4719 random K shapes every tilt-unstable onset had the tilt tail slowest (m_u²/m_w² ≥ 4.6), so like valleys bind and the onset is first order into a same-tilt chiral soliton lattice (c = 1), as in C161. No mixed or alternating soliton lattice is a classical minimizer at any density.
- Construction side (at the coordinator's request after Jonah flagged the c = 2 route): an accidental two-valley onset in N must sit a finite distance from K and U and cannot come from deforming the K tilt mechanism; oscillatory tails flip the valley-odd sign but bind every channel semiclassically; the easy-plane pseudospin window is exponentially thin (J1–J2 magnons reach Δ of order one only because they are point-like); at ε = 0.1, L ≤ 256 the pseudospin scale is about 1e-7, below any DMRG resolution.
- Remaining assumption (A_sol): no accidental two-valley soliton onset in N or T far from K and U, at ε of order one, combines non-binding interactions with J_z > −|J_⊥| (the order-n² mass-imbalance correction is not derived). No c = 2 bond found, no Typhon job.

### 9.122 Lattice search for accidental two-valley soliton onsets (C164, Claude Opus agent, 2026-11-08; `codex/work_C164/`; Typhon job 17601039 pending), and C160's Typhon output

- Calibrations reproduced: C161's tilt eigenvalue (−0.0671 against −0.067), C162's K-stratum like-valley binding (−0.0084 against −0.0088 at separation 10), DMRG c = 0.977 for XX at L = 32.
- Local coverage [numerical]: 188 classical lattice bonds (39 with genuine Hessian-tested soliton minima, none with a second valley); at ε 0.3 to 0.5 the solitons are 10 to 28 sites wide and the quantum band sits far above the classical energy, so the classical screen is a poor guide there. 610 quantum ring screens and 16 optimizer restarts in N found smooth near-degenerate two-valley bands (best: valleys at k = 0.43 and 2.15, energies 0.273 and 0.275, bandwidth 1.04), but driving toward onset merged the valleys into one flat bottom (best residual gap 0.066 of the bandwidth). Like-valley pairs did not bind; J_z could not be extracted on 14-site rings.
- Mass imbalance [numerical, derived], the term C162 left open: each A–B interface costs about 0.3η²n²/M (six-particle energies ordered by interface number), which dominates the n³ scattering-length coupling as n → 0, so unequal-mass two-valley gases are immiscible in the dilute limit. T reduces to accidental cases by LSM. A new exact k → π − k symmetry on C152's cos φ = 0 surface (checked to 7e-14) is the two-Majorana case already covered by C152.
- Resolution: particles about 3 sites, density 0.1 to 0.2 past onset, pseudospin scale 3e-3 to 3e-2, valley-splitting bisection resolving about 1.6e-3 of the bandwidth. A miscible window narrower than about 1e-3 would be missed. Residual (A_sol′): no accidental two-valley onset in N or T at ε of order one has separated non-binding valleys, J_z > −|J_⊥|, and n·Δa above about 0.1η² (or near-equal masses). No c = 2 bond found so far. Job 17601039 runs the large screen, controls (XX to L = 256, two decoupled XX chains as a c = 2 control) and the past-onset DMRG.
- C160's job 17599583 finished (`codex/work_C160/typhon_out/summary.txt`). DMRG of the K-stratum chiral soliton lattice family at ε = 0.5 to L = 256: the U(1) calibration reads 0.997, the CSL points oscillate with soliton commensuration between about 0.7 and 1.10 with no trend above 1, and the canted side reads 0. The lattice texture test failed its own pre-fixed criterion: 46 of 48 bonds have unexplained (non-harmonic) peak weight above 1e-2, with clusters of nearby wavevectors that look like unconverged multi-domain minimizations rather than quasi-periodic textures. Per the criterion they must be re-minimized before anything is read from them, so the lattice version of M1 remains untested.

### 9.123 C160 lattice texture test redone: all 48 bonds uniform or single-period (coordinator, Claude Opus 5.5, 2026-11-08; `codex/work_C160/remin/`; Typhon job 17601046)

- The original test (texture_task.py) failed because of its analysis, not the textures. Its Hanning window has sidelobes at ±2.5 bins at −31 dB, its peak finder counted them as peaks, and its ±2-bin weights double-counted main-lobe power, so a clean helix scored about 0.5 of "non-harmonic" weight. On near-uniform bonds the open-chain DM edge twist also reached the bulk window.
- Redone (remin.py): best uniform cone helix from the exact per-site energy; ring cells up to P = 300 near its pitch; the N = 2000 open chain re-minimized from the tiled best ring and 8 fresh starts; bulk energy on the middle half compared with the best ring cell; Blackman-window spectrum (sidelobes below −58 dB) with every peak above 1e-3 tested as a harmonic of the strongest; single-period cells matched to the fitted pitch up to P = 1200.
- Result [numerical]: 34 bonds uniform (bulk variance below 1e-6 and bulk energy within 1e-7 of the best ring cell, mostly 1e-8), 14 singly periodic spirals (q between 0.17 and 0.19, peaks only at q, 2q, 3q), no candidate. The lowest open-chain bulk energies lie within 1e-4 of the best single-period cell, the residual differences coming from pitches not representable by P ≤ 300 cells and slow phason relaxation. One bond (seed 10, t 2) crashed in least squares on Typhon and was rerun locally at N = 1200, uniform. Scope: the C160 ensemble, random bonds in the doubly rescaled pure-DM wedge at ε = 0.2.

### 9.124 A staggered Ising transition inside a floating liquid does not occur at spin 1/2 (C165, Claude Opus agent, 2026-11-08; `codex/work_C165/`; not yet refereed; Typhon job 17601089 pending)

- Classical [proved, derived, numerical]: the rigid staggered cone is linear in cos²θ for every bond, but at second order in the U(1)-breaking harmonics the N⁴ coefficient at the switch is c₂ − c₁, with c₂ = g²cot²q/4 from the charge-2 XY anisotropy (favours coexistence) and c₁ = G²(1 − cos q)²/(2(1 + cos q)) from Γ terms (first order); DM and fields do not enter. Checked against direct minimization to 4 digits. On the classical lattice, at g = 0.3, hx = 0.2, hz = 0.1 (in N), the spiral's momentum-π Hessian mode softens while the spiral is the global minimum, followed by a coexistence window of width about 0.05. So the construction exists classically.
- Spin 1/2 [proved mapping, derived]: the parent is XXZ about z with DM along z, which a site-dependent rotation and a π rotation on odd sites map exactly to the antiferromagnetic XXZ chain with Δ = D_z in a longitudinal field. The spiral/staggered switch is its Heisenberg point (zero field) or spin-flop line (field). There the staggered moment and the two spiral components are the three components of one SU(2)₁ primary, and the classical quartic curvature is a spin-4 harmonic of dimension at least 8. To all orders in weak harmonics the infrared is the XXZ dichotomy, a c = 1 liquid or gapped Néel, with no coexistence.
- Numerical (two-site VUMPS, χ ≤ 40), line a = 0.9, g = 0.3, hx = 0.2, hz = 0.05, Dx = 0.03 (in N): at D_z = 1.0, inside the classical window, the staggered moment decays with χ (0.082 to 0.050) and c reads 0.96 to 1.07; from D_z = 1.25 the moment saturates (0.26 to 0.41) and c falls to 0.09, gapped Néel setting in near D_z = 1.2. No point is both ordered and gapless.
- Remaining assumption (A_S): with harmonics comparable to the exchange near the hidden Heisenberg point, no lattice-scale excitation independent of the liquid's boson orders at momentum π while the liquid survives. Job 17601089 runs a 252-point screen at amplitudes 0.2 to 0.6 with follow-up to χ = 256 and L = 256, and an Ising × XX (c = 3/2) calibration.

### 9.125 Typhon output of C163 and C155: open-chain and infinite-chain numerics disagree near XZ (coordinator, Claude Opus 5.5, 2026-11-08; `codex/work_C163/typhon_out/`, `codex/work_C155/typhon_out/`)

- C163, job 17600999 (completed). Open-chain DMRG at L = 192 and 256 (χ 192, truncation below 1e-9). hi_T_40 reads c_prof 1.23 and 1.27, with two-size slopes 6 dS/d ln L of 1.08, 1.21, 1.36, 1.44 over 32 → 64 → 128 → 192 → 256. dist_T_50 reads 1.32 and 1.31, with slopes 1.20, 1.34, 1.46, 1.41. Neither falls towards 1. At fixed distance x from the edge S(x) still grows with L, and the profile is CFT-shaped with c ≈ 1.25 in the chord variable. hi_T_37 falls to 0.35 (gapped). hi_K_1 reads 1.45 at L = 256 after 1.37 and 1.26 at L = 64 and 128, so it is non-monotone.
- Two-site VUMPS on the same two bonds (χ 8 to 192, `fes/`) gives the opposite picture. S saturates at 0.7007 (hi_T_40) and 0.6834 (dist_T_50) from χ ≈ 96, while ξ grows only from 39 to 46 and from 30 to 35. The FES slope is about 0.01 at χ ≥ 96, and the energy converges to 1e-9 over χ = 96 → 192, the behaviour of a gapped state. The VUMPS energy per site (−0.95930920, −0.83712373) is variational and lies below every open-chain two-size energy slope. Those slopes approach it from above (gaps 2.7e-4, 9.7e-5, 4.5e-5, 2.6e-5 for hi_T_40), so both methods are consistent with the same e₀. The leading transfer-matrix eigenvalues are a complex pair at arg ±0.38π and ±0.60π (incommensurate correlations). The open chains have S(1) = ln 2 exactly, so the DMRG state is time-reversal symmetric.
- Verdict [numerical]: by C163's pre-fixed criterion (an FES slope above 1.2 is required) neither point is a c > 1 candidate. The conflict is not resolved, however. A gapped bulk with ξ ≈ 46 should saturate an open chain of 256 sites, and it does not. Candidate explanations, none tested: (i) the bulk breaks time reversal weakly and the symmetric open-chain state is crossing over to a cat (S_∞ + ln 2 ≈ 1.39, although the per-doubling growth of S_mid is still increasing at L = 256); (ii) a two-site cell or χ ≤ 192 traps VUMPS in a gapped state that is not the ground state (its energy is lower than any open-chain estimate, which argues against this but does not exclude it); (iii) a crossover length well above 256 sites.
- C163 twisted curve: 16 bonds near γ read c_prof(128) ≤ 1.10, consistent with Proposition G. Paths into the attractive window: 40 paths, of which three meet the first-stage criterion (c_two and c_prof(128) above 1.2, K_S(16) > 15/8). They are KI22 (Δ = −0.75, c_prof 1.31 at L = 64 and 128, c_two 1.76, K_S 3.2), KI26 (1.64 at 128 but −0.35 at 64) and N34 (1.68 at 128, −0.24 at 64). All three sit on a jump in S_mid(λ), that is at a level crossing, and the latter two are noisy. FES at χ ≥ 256, the second stage of the criterion, has not been run.
- C155, job 17594766, hit its 11 h limit before writing summary.txt. Its summarize.py was run locally on the copied outputs.
  - Stage A: 12,000 uniform points; 636 gapless by the ring criterion; 28 with c_obj > 1.1 and gapless.
  - Stages B and C: 334 gap-closing and 367 c-maximization restarts, 150,328 evaluations. The ring objective saturates at its cap (4.0), so the screen alone carries no weight.
  - Stage D: 78 candidates refined by DMRG to L = 128 (truncation below 1e-8). Most read 0.8 to 1.1. Those still rising or above 1.3 at L = 128 are cands 60, 61, 63 (L64 → 128: 1.28 → 1.70, 1.03 → 1.91, 1.19 → 1.71), 16 (0.98 → 1.45), 59 (1.08 → 1.37), 22 (1.04 → 1.74), and 1, 2 (about 1.4 at both sizes).
  - Every one of these except 16 and 22 (which carry a field 0.44 along z) has nearly the same bond as dist_T_50: diagonal about (−0.58, 0, 0.55) with DM along y about −0.6, i.e. the XZ neighbourhood.
  - Stage E (104 XZ directions and Kitaev-line points): the U(1) directions read 1.00 at L = 256. T8 at d = 0.3 reads 1.24 and 1.27 at L = 128 and 256, but its two-size slope falls to 0.58. No point of stage E rises at L = 256.
- Summary: no c > 1 bond survives a pre-fixed criterion. Every reading that keeps rising sits in one region, the XZ neighbourhood, where open-chain DMRG and infinite-chain VUMPS disagree. Deciding it needs the same test on all of them: VUMPS with cells commensurate with the incommensurate wavevector and χ up to 512, open chains to L = 512, bond-energy comparison in the middle of the chain, and a pinning-field test for time-reversal breaking.

### 9.126 The XZ conflict is a Néel cat, and the C164 two-valley search is null (coordinator, Claude Opus 5.5, 2026-11-08; `codex/work_C166/`, `codex/work_C164/typhon_out/`)

- C166 (Claude Opus 5.5 agent, spec `specs/C166_xz_reconcile.md`), local runs. The pass criterion (C1 to C3 in its RESULT.md) was written before any new data.
  - Codes agree [numerical]. Open-chain DMRG matches ED to 7e-13 at L = 8, 12, 16 on all seven bonds, the VUMPS bond operator equals the DMRG one to 1e-16, and C163's χ = 16 VUMPS energies are reproduced to 1e-10. Hypothesis (iv) of §9.125 is excluded.
  - VUMPS cells n = 1 to 8, 10 to 14, 20, 21 (two seeds, χ ≤ 32). Every even cell converges to the same state as the two-site cell, equal in energy to 1e-12. Odd cells give non-injective cats at higher energy. No evidence for hypothesis (ii) at this χ.
  - The bulk is Néel ordered along z [numerical]. The staggered moment is |m| = 0.345 (hi_T_40) and 0.359 (dist_T_50), bond energies are uniform, and the state is not time-reversal invariant. C163's χ = 192 VUMPS files already carry this moment.
  - Unpinned open chains are the symmetric superposition of the two Néel states. The staggered correlation at r = L/2 is 0.076 and 0.073 at L = 64 and 128 for hi_T_40 (m² = 0.119), nearly flat, where a critical state would decay. The middle bond energy sits about 1.5e-5 above e₀ and does not fall like 1/L².
  - Pinning test [numerical]. Edge fields at both ends along the VUMPS order turn the open chain into the VUMPS state. For hi_T_40 at L = 64 and 128, S_mid is 0.671 and 0.690 (VUMPS 0.7007), c_prof is 0.48 and 0.29, the central staggered moment matches VUMPS to 1e-3, and the middle bond energy is within 4e-6 of e₀. dist_T_50 at L = 64 behaves the same way (S_mid 0.653, c_prof 0.48). A field perpendicular to the order leaves the profile nearly unchanged and S_mid still growing (0.918 → 1.039).
  - Reading [numerical, with the length scale heuristic]. This is hypothesis (i) of §9.125, with strong rather than weak order. The unpinned chain's extra entanglement relaxes towards ln 2 over a length of about 130 to 190 sites, much longer than ξ ≈ 46, and the rising c ≈ 1.3 is that crossover. Why the length is so long (cheap domain walls emitted at a free edge) is a guess.
  - Registered prediction before the Typhon data: S_mid(512) ≈ 1.35 (bounded by S_VUMPS + ln 2 = 1.394) for hi_T_40 and ≈ 1.36 (bounded by 1.377) for dist_T_50, against about 1.38 and 1.43 for a critical bulk with c ≈ 1.3. Exceeding a bound by 0.02 falsifies the reading.
  - Typhon job 17602270 (pending at write-up) runs open chains to L = 512 on hi_T_40, dist_T_50, cand_60 and cand_61, pinned chains at L = 128 and 256, VUMPS to χ = 512, and the FES second stage on the C163 path points KI22, KI26, N34. At χ = 32 N34 is a near-product state (S = 0.013), consistent with a level-crossing artefact. KI22 and KI26 wait on the large-χ runs. Collection: scp `oc`, `vc`, `logs`, `summary.txt` into `codex/work_C166/typhon_out/`, then `python3 summarize.py typhon_out`.
- C164, job 17601039 (completed, 7 h 24 min), output in `codex/work_C164/typhon_out/`.
  - Classical screen: about 48,000 lattice bonds in N, T, the DM class and near h_F. 10 classical candidates have an isolated soliton below the magnon gap with a second valley within 30 percent in energy.
  - Quantum optimizer: 5992 restarts of the two-valley objective, 2424 tuned to onset, 35 endpoints with a small two-valley gap taken past onset by DMRG.
  - Past-onset DMRG: the largest c_prof at L ≥ 64 among the 35 endpoints is 1.06, then 0.84, 0.75. None reaches the 1.2 promotion threshold, so no point went on to L = 128 for c, and none meets the c = 2 criterion. Most endpoints sit in a gapped phase on one side of λ* (c ≈ 0, structure-factor peaks at 0 and π).
  - Controls: XX gives c_prof 0.975, 0.983, 0.990 at L = 64, 128, 256 and an FES slope of 0.98 (χ 128 → 256). The two-XX c = 2 control gives an FES slope of 1.91 (χ 96 → 128), but its open-chain DMRG is under-resolved (1.73 at L = 64 with truncation 1.3e-5, and L = 32 failed with an SVD error). The pipeline sees c = 2 by FES, and its open-chain c = 2 reading at L = 64 runs low.
  - Verdict [numerical]: null within the coverage above. No accidental two-valley onset with a miscible gas was found, and (A_sol) stands as the residual assumption in the form C164's RESULT.md gives.

### 9.127 Red team of the conditional theorem (C174, Claude Opus 5.5 agent, 2026-11-08; `codex/work_C174/`)

Part of the swarm C169 to C174 (`specs/C169_C174_swarm.md`), launched alongside C167 (A_W) and C168 ((H)/(A6)). C174 attacked the logic of C154 PROOF 7, which turns c ≤ 1 on every symmetric class plus (H) into c ≤ 1 on B_2. The implication is valid. It is not sound modulo only (A_W), (H)/(A6) and A_S.

- R1 [proved, given (A2) and gap stability; checked by the coordinator]. The closure of a c > 1 component contains no gapped bond and no controlled conformal point. So on c > 1 components (H) is equivalent to (H*), which says every such component accumulates on a controlled singular point. With (A5), that is the conclusion on N itself. Upper bounds only travel from higher c to lower c, so the theorem reduces B_2 to N and proves nothing inside N. A connectivity argument cannot establish (H).
- R2, R3 [high]. (A5) has to be restricted to an explicit list S_ctrl of controlled singular sets, with control on full punctured neighbourhoods rather than direction by direction. Three neighbourhoods are not yet covered in that sense.
  - XZ, along curved approaches; the Néel-cat region of §9.126 sits here.
  - h_F at soliton onsets, where ε → 0 and L → ∞ do not commute (C161).
  - The onset layers along γ.
- R4 [medium-high]. The stratum results hold at generic points. A c > 1 component could meet a stratum only at an isolated multicritical point, and those are not inventoried: the tricritical points in the parity class, the large-b corner of the staggered-charge class, and the SU(2)₁/ordered boundary of the cyclic stratum.
- R5 [medium-high]. The candidate list misses irreducible c > 1 theories with one species and only gapped or c ≤ 1 neighbours. Examples are Z_N parafermions with N ≥ 5 (c = 8/7 at N = 5, codimension 2), N = 1 superconformal minimal models, W₃ minimal models, and SU(2)₃ off the cyclic stratum. A natural home is the melting of Z_m-locked tongues with m ≥ 5 at order-one harmonics. This is new pattern P4, passed to C168 and C169.
- R6 [medium]. The K_II sheet is under no assumption. C159's counting leaves a repulsive c > 1 window there.
- R7 [medium]. The chain also uses (A1), (A2), C160's semiclassical correspondence, (A_sol′), the Schrieffer–Wolff step in C163 Proposition G, Lemma 3 at spin 1/2, and jobs that were never collected (17597010, 17596442, 17594765).
- R8 [partial fix, numerical to 2e-14 for p ≤ 7]. Only one bond crosses each cell boundary, so det h(K) is a degree-2 trigonometric polynomial for any period p. In the parity class this gives c ≤ 1 at mean-field level for every period, and at most two soft momenta at any onset. The central charge c must be defined per superselection sector, or cats and first-order surfaces read c = 2.
- R9 to R12 [low]. The strata list is complete for on-site product symmetries, and spatial-reflection strata are unanalysed. Gap stability, constancy of c on a component, and the torus count are sound. No connectivity test in N has been completed (C68's point 1282).
- Reading. What c ≤ 1 at d = 2 needs is a dynamical exclusion of c > 1 fixed points in N away from S_ctrl, P4 included, plus junction inventories on each stratum. The critical path is longer than the three named assumptions.

### 9.128 Referee of C165 (C170, Claude Opus 5.5 agent, 2026-11-08; `codex/work_C170/`)

The main conclusion survives. For weak harmonics there is no staggered Ising transition inside the liquid at spin 1/2, so c ≤ 1 near the switch. Three statements are wrong, and one regime needed a new argument.

- §1.2, the exact map of the parent to the spin-1/2 AF XXZ chain in a field: correct. Spectra agree to 5.5e-14, and the full bond commutes with T_μ exp(−i(a−π)S^z) to 4e-15. The ν frame is classical bookkeeping only.
- §2.1, rigid-cone linearity: correct for every q ∉ {0, π}.
- §2.2, the classical N⁴ coefficient: wrong for DM. "DM gives no quartic term" fails because in the ν frame the bond matrix is not antisymmetric.
  - Corrected coefficient c_2 − c_1, with c_2 = (g² + g_xy²) cot²q/4 and c_1 = sin²(q/2)|tan(q/2)Γ_c − iD_c|², where Γ_c = Γ_xz − iΓ_yz and D_c = D_y + iD_x.
  - It matches exact linear response on eight random harmonic vectors. In a full minimization, D_x = 0.25 turns the g = 0.2 coexistence window into a first-order jump.
  - Charge-2 terms remain the only convex source. The sentence "DM and fields do not enter" in §9.124 is corrected here for DM. The fields part stands.
- §4.1 and §4.2, the SU(2)_1 multiplet with the quartic of dimension ≥ 8, to all orders in weak harmonics: correct at Δ_eff = 1, h_eff = 0, and stronger than claimed.
  - The selection rule P − Qα ∈ 2πZ is exact, so no relevant momentum-π operator is ever generated while the infrared stays in the SU(2)_1 basin.
  - A uniform field is generated at second order (ED scaling as amplitude²), so the point sits at h_eff = 0.
  - Charged relevant operators appear only in resonance tongues near q ∈ {0, π}.
- The spin-flop line: a gap, now closed [derived]. At hz ≠ 0 the switch is the Néel plateau edge, a Pokrovsky–Talapov transition, so the dimension-8 argument does not reach it. The replacement:
  - off resonance the chain is U(1) with one boson;
  - coexistence needs a magnon bound below the two-spinon threshold;
  - weak harmonics cannot create one.
- "The quantum shift of the Néel onset is O(1)": wrong at hz = 0, where both shifts are O(harmonics²).
- "Gapped Néel up to hz = 2.2 at Dz = 2": wrong. The plateau ends at hz_c1 = 0.78. The VUMPS energy at hz = 2.2 is about 0.1 above two ED estimates, so those states are metastable.
- A_S restated. Along the quantum switch set (the KT point and the plateau edge h_eff = hz_c1(Δ_eff) for Δ_eff > 1), at order-one harmonic amplitude, there is no second low-energy sector, in particular no bound magnon. The old wording "fields up to the plateau edge" covered an exponentially small range (hz_c1 = 0.0004 at Δ = 1.1).
- Caveat for C165's Typhon job 17601089, family A2. Its fields mostly lie above the parent plateau edges, so any A2 point that reads gapped Néel above the edge needs an energy comparison against a magnetized start before it is used. The restated A_S has been passed to C172.
- C165 §4.4 and §4.5 were checked for consistency only, not refereed.

### 9.129 (A_W) narrowed to its charge-1 part (C167, Claude Opus 5.5 agent, 2026-11-08; `codex/work_C167/`)

- Lemma R0 [proved; checked to 1.2e-13 on 200 random bonds]. The U(1)_m-breaking harmonics split by charge. The charge-1 parts are the parts odd under the π rotation about the spiral axis m, so a bond with no charge-1 part is a parity-class bond (XYZ + D_m). In T the charge-1 part is the exchange block S_{m⊥}. In K_I it is the in-plane field and nothing else. (The coordinator checked the identity.)
- Proposition R [derived; inputs are C152's G1 liquid, π-parity and momentum bookkeeping]. Expand to every order in the charge-1 part, at any charge-2 amplitude. The infrared is one boson, with liquid, KT-ordered period-2 and clock-locked phases, all with c ≤ 1. The large-charge-2 part of (A_W) therefore rests on A4 (c ≤ 1 on the parity class), already an input of C159. C163's Proposition P is the special case of small charge-2.
- Lemma C [derived, constant not computed]. The charge-1 harmonic couples to the boson's own e^{iϑ} at the twist θ_F. The liquid survives only for g < C₀ v (K θ̃)^{2−1/(4K)}, so a lattice-size charge-1 harmonic needs θ̃ of order one.
- Lemma Θ2, case (b) [derived]. The antiunitary pairs each composite minimum P* with −P*. A c = 2 phase then has one of two forms.
  - With the antiunitary unbroken, the second boson needs K₂ > 2, so a local operator at ±P* outside the charge lattice has quasi-long-range order with η₂ < 1/4.
  - With the antiunitary broken, there is spontaneous uniform magnetization along m at zero field.
- Lemma O, case (a) [derived]. The diagnostic is x_π, the dimension of the leading momentum-π operator. In the liquid x_π = K, and the boson orders by KT at x_π = 1/2. Next to an independent Ising surface x_π drops to 1/8.
- Negative results.
  - The antiunitary forces no mixing.
  - Perron–Frobenius covers only θ ∈ {0, π} with no charge-1 part.
  - Concavity of E₀ gives upper bounds only.
  - C165's SU(2)₁ mechanism has no analogue here.
- hi_K_1 [numerical; local VUMPS on two- and four-site cells, χ ≤ 96]. This was the one window point whose open-chain c read above 1 at L = 256.
  - The bulk is period-2 ordered, with a staggered moment of 0.183 converged in χ and c_pair = 0.08.
  - Pinning the edges along the order lowers c_prof(64) from 1.37 to 0.62. Like the XZ points of §9.126, the open-chain reading is a crossover in an ordered bulk.
  - Along its harmonic path, x_π starts at the parent K = 1.06 and falls continuously through 1/2 (0.98 → 0.47 for λ from 0.90 to 0.98), and c_pair never exceeds 1.06. The order is the boson's own KT, so the liquid before it is repulsive and hi_K_1 is not in the window. The ring K_S is not a reliable window test at large harmonics.
- Two attractive T points on the infinite chain [numerical, χ ≤ 48]. ising_T_43 reads c_pair 1.00. hi_T_34 reads 0.91 to 0.98, and its slow in-plane moment is unidentified.
- Residual (A_W1) [assumed]. It applies when the charge-1 harmonic is about 0.3 of the exchange or more and θ̃ is of order one, with charge-2 arbitrary.
  - Case (a), K > 15/8: x_π stays equal to K.
  - Case (b), K > 1: no c = 2 phase of either form in Lemma Θ2.
  - It bites for −1 < Δ < −0.67 in case (a) and −1 < Δ < 0 in case (b), away from h_F and γ. There is no small parameter, extra symmetry or hidden Heisenberg point to continue with.
- Typhon job 17602938 (pending) runs:
  - the hi_K_1 path to χ = 192 (256 at λ = 1) and past λ = 1;
  - a four-site cell, and the charge-separated bonds;
  - pinned open chains at L = 128 and 256.
  It writes `summary.txt`. Collection: scp `summary.txt`, `vc`, `oc` into `codex/work_C167/typhon_out/`. Prediction registered before the output: x_π crosses 1/2 between λ = 0.98 and 0.985, c_pair stays at or below about 1.05, and pinned S_mid stays within 0.05 of the VUMPS entropy.

### 9.130 Referee of C159 and C163: the T and K_I reduction to (A_W) is incomplete (C171, Claude Opus 5.5 agent, 2026-11-08; `codex/work_C171/`)

Most of the individual lemmas hold. They were re-derived by hand, with an operator-counting script (`opcount.py`) that reproduces every row of C159's table 4.1. But five routes to c > 1 on T and K_I end neither in (A_W) nor in a derived exclusion.

- C159. Correct as written:
  - Lemma 1, the normal form of K and its sheets;
  - Ogata–Tasaki on T;
  - the liquid-stability conditions;
  - Lemma B (the D and E series are assumed, as tagged);
  - the SU(2)₂ and SU(2)₃ exclusions;
  - the rows of table 4.1.

  Gaps and errors:
  - Lemma A is correct for one boson on which translation acts as a shift. It does not cover translation acting by charge conjugation, nor two bosons swapped by time reversal.
  - Table 4.1 as a complete list, and Proposition C as used: wrong.
  - Dismissing codimension ≥ 2 candidates, and replacing H_T by A6_TK: wrong.
- Missing routes.
  - R1. Three-state Potts order at momentum 2π/3 × a half-floating boson, c = 9/5. It needs one tuning at every K > 1/2 and no attraction. C154's table for N lacks it too, though (A6) covers it as worded. Whether a chiral Huse–Fisher term makes the transition non-conformal is not settled.
  - R2. Codimension-2 loci are generic, so "they must lie near symmetric strata" is a non-sequitur. In T, uniform Θ-odd Ising × boson (P2) has codimension 2 at any K and only c = 1 neighbours.
  - R3. C163's Lemma L locks Q₂ = nQ₁ only for n = 1, 2. For n ≥ 3 it holds only where n²/(4K₁) + 1/(4K₂) < 2. The nQ₁ + π branch is wrong, and n = 0 is untreated. Two-boson phases with wavevectors inside Zθ + Zπ survive at codimension 0.
  - R4. A second boson whose θ is the field time reversal flips has x₂ = 1/(4K₂). This gives a codimension-0 c = 2 candidate with a repulsive host (K₁ = 0.9, K₂ = 0.2).
  - R5. A sector on which translation acts by charge conjugation reaches codimension 1 for K₁ > 3/2.
- C163.
  - Fact 1: correct, tag [derived].
  - Fact 2: correct.
  - "Every soft wavevector lies in Zα + Zπ": wrong at finite field, harmless on T and K_I.
  - Lemma U: the count is correct.
  - "Dimerization needs frustration": wrong, but not load-bearing.
  - Proposition P: correct, with the 0.3 threshold heuristic.
  - Proposition G: correct as conditionally tagged. It covers only the large-K corner, not the window.
  - KI22 is still an open first-stage pass, and its FES is in C166's job 17602270.
- Corrected residual (A_W′) [assumed]. In T and K_I, outside the controlled neighbourhoods, no gapped sector becomes critical or forms an independent liquid while a half-floating liquid is present, at any K and codimension. This is C154's (A6) restricted to T ∪ K_I. The "harmonics comparable to the exchange" restriction survives. The "attractive" restriction survives only for the codimension-1 Ising-at-π row (K > 15/8), the charge-conjugation row (K₁ > 3/2) and the composites Lemma L locks. C167's charge-1 reduction (§9.129) applies to the (A_W) part. K_II is covered by none of this.
- Two lemmas would restore part of the narrowing:
  - an exclusion of Z_N order (N ≥ 3) independent of the boson, the period-3 analogue of C165;
  - a dynamical argument that a codimension-2 c > 1 point forces a codimension ≤ 1 one nearby. Kinematically this is false for P2.
- The new rows were passed to C168, C169 and C173. The Potts × boson route is now a construction target for C169, alongside C174's parafermion route.

### 9.131 Typhon output of C158 (FES near XZ) and C165 (staggered Ising screen) (coordinator, Claude Opus 5.5, 2026-11-09; `codex/work_C158/typhon_out/`, `codex/work_C165/typhon_out/`)

- C158, job 17595611, completed in 18 h 39 min. Digest: `fesfit_summary.txt` (C158's fesfit.py). Per C158's rule, the bisector chain is read only at bond dimensions where the incommensurate correlation length ξ_ic leads and ξ_π has saturated.
  - Bisector field [numerical]. At hp = 0.2, 0.3, 0.4 and 0.45, c = 1.00 ± 0.02 out to χ = 384. The three-point fit c_ic3 is 1.001, 0.997, 1.002 and 0.985 (two-site cell). Four-site cells, a second seed and a one-site cell agree (0.99 to 1.06). This is one boson, consistent with C158's 1 : 2 : 4 tower and with C157.
    - hp = 0.05 and 0.1 stay in the crossover: ξ_π leads, VUMPS errors of about 1e-4 are unconverged at 800 iterations, and the pair slopes swing in sign. They are undetermined, not c > 1.
    - hp = 0.5 is gapped up to χ = 160, then jumps at χ = 192 to an unconverged state (error 5e-2). It is undetermined. hp = 0.6 is gapped, as the control.
  - The five tuned-hyperplane points are all gapped, with c fits of 0.000 to 0.006 and S saturated.
  - YX at ε = 0.1 to 0.5 and ZY at 0.2 give c = 1.00 (six-site cell). YX at 0.05 reads 0.46, in C154's dimerization window.
  - Ashkin–Teller plane at β = 0.3: α = 0.18 and 0.198 are gapped. α = 0.09 and 0.3 have falling slopes (0.37 and 0.74 at χ = 256), so nothing is above 1.
  - Controls. XX and XZZX give 1.000, Heisenberg 1.007, and the XZ ε points 0.98 to 0.99. The critical TFI reads 0.50 at χ ≤ 48 but drifts to 0.34 at larger χ, where VUMPS stops converging (800 iterations). So late-χ fits on very long ξ are not reliable.
  - Open-chain DMRG of the bisector chain: S_mid = 1.41 (hp = 0.1, L = 192), 1.25 (0.2, L = 256), 1.21 (0.4, L = 256). There is only one size each, so these give no c.
  - Verdict: no point near XZ has c above 1 where FES is past its crossover. The small-field end hp ≤ 0.1 of the bisector line is still unresolved.
- C165, job 17601089, hit its 14 h limit. No summary.txt was written, and the DMRG directory `dm/` is empty, so the long open-chain tasks did not finish. C165's summarize.py was run locally on the copied VUMPS output (`summary_local.txt`).
  - Calibrations. The Ising × XX ladder (c = 3/2, outside B_2) gives an FES slope of 1.49 at χ = 256, so the pipeline sees 3/2. XX gives 0.997. TFI again drifts low (0.37) at large χ.
  - The N-stratum floating control (`nfloat`) is not a c = 1 control. It reads 0.03 at χ = 256 with a staggered moment of 0.14 to 0.16, so it is ordered and gapped.
  - Screen (252 points, χ ≤ 64). 169 points are ordered and 0 are flagged as coexisting (floating + staggered order). On the 24 stage-B points continued to χ = 256, every c fit is at most 1.05, and the ordered ones go to 0.
  - Q2 at Dz = 1.2, 1.25, 1.3 is gapped Néel (c → 0.12, 0.02, 0.00, with N saturating at 0.18, 0.25, 0.31). This agrees with C165's local reading of the onset near Dz = 1.2.
  - Six A2 points read c_last of 1.3 to 3.0 at χ = 64. Their fields (hz = 1.6 to 2.7) lie far above the parent plateau edges (0.27 at Dz = 1.6, 0.78 at Dz = 2), where C170 showed the gapped-Néel VUMPS states are metastable. They sit in the magnetized liquid regime. The one continued to χ = 256 (A2_Dz1.6_h0_hz2.0) falls to 0.23. The other five were not continued, and a χ = 64 last-pair slope is noise at this level. Per C170's caveat, none is usable without an energy check against a magnetized start.
  - Verdict by the pre-fixed criteria: no c = 3/2 candidate. Two things are incomplete: the DMRG stage, lost to the time limit, and the five uncontinued A2 points, which need χ ≥ 256 from a magnetized start before they say anything.

### 9.132 A_S narrowed to a two-kink binding statement (C172, Claude Opus 5.5 agent, 2026-11-09; `codex/work_C172/`)

The target was C170's restated A_S: no bound magnon along the quantum switch set at order-one harmonic amplitude.

- [proved] An exact on-site SU(2) exists only on the twisted Heisenberg curve (J a multiple of a rotation, B = 0). In N, ring momentum is the only exact quantum number, so the SU(2)_1 towers have no non-crossing protection. Continuity of the towers and a hidden exact symmetry both fail as routes.
- Locking lemma [derived]. An independent momentum-π Ising sector locks to any relevant liquid vertex whose momentum lies within the locking scale of π. On T this recovers C159's second tuning. Near T, any c = 3/2 hypersurface of N is pinched onto that two-tuning locus.
- Kink-parity dichotomy [derived]. Kink parity is exact on rings. At the edge of a period-2 phase, a route to c = 3/2 needs a kink-even state to soften before the kink, which means a bound state below the two-kink continuum. This is C170's bound magnon, now without needing a U(1). Any binding opens such a window, so the exclusion has to be exact.
- Lemma W, the Ising end [proved at first order, derived for the chain]. Write the transverse bond as A σ⁺σ⁻ + P σ⁺σ⁺ + h.c.
  - For any charge-0 and charge-2 harmonics and any longitudinal field, no bound magnon exists when |P| ≤ |A|, with margin 2(|A| − |P|).
  - The threshold is sharp. At a = π a bound state appears for every |P| > |A|, matching −(w + t²/w) to four digits at L = 30.
  - For C165's family, |P| ≤ |A| means g² + g_xy² ≤ 1, the whole polar-parent range.
  - Charge-1 terms are not covered by the proof. Numerically, in 192 bonds they never move the threshold.
- The KT end [derived, conditional]. If the KT point stays in the SU(2)_1 basin at that amplitude, the massive spinons of the SU(2) Thirring model have no bound states. The width of that neighbourhood is not computed.
- Ring ED, L = 10 to 18, at Dz = 1.6 and 2 [numerical]. Harmonics of amplitude up to 1.4 shift the binding measure (kink-even gap minus twice the kink gap) by at most 0.105, against twice the kink gap of 0.9 to 1.9. Along field paths the kink gap closes first. The resolution is about 0.1, so weak binding is not excluded.
- Residual A_S′ [assumed]. On the plateau-edge part of the switch set at intermediate anisotropy, between the KT crossover and the Ising regime, at order-one amplitude with |P| < |A|, the two-kink channel has no bound state below the two-kink continuum. This is a two-body statement on a gapped background, testable directly by odd and even sector DMRG at resolution well below 0.1. The liquid interior stays under (A6), with the locking lemma's detuning condition added.

### 9.133 The symmetry-free stratum N: residual (R_N), and a corrected sampling bound (C168, Claude Opus 5.5 agent, 2026-11-09; `codex/work_C168/`)

(H)/(A6) is not derived. It is reduced to a strictly narrower residual, and no c > 1 bond was found.

- Kinematics cannot do it [proved/derived]. A symmetry-free qubit chain with three-site terms has c = 2 (two floating bosons) on an open set. Translation alone carries no LSM anomaly. So any exclusion in N must use the nearest-neighbour range and the dynamics.
- Lemma O [derived]. This replaces C154's Lemma 3 in N, as C174 R9 asks, and works for any gapped periodic reference. At a generic onset one species softens, so c ≤ 1. Translation acting by charge conjugation, and swapped boson pairs, come only from a pair of real modes at 0 and π, which gives one boson.
- Lemma D [derived, with numerical checks]. Two independent species softening together out of a gapped phase is a codimension-2 event. Dilute particles in one dimension are impenetrable, so the mass imbalance costs 0.31 η² n²/M per unlike neighbour pair (C164 found 0.3). Passing costs only order n³. So the two-species gas separates or alternates, and c ≤ 1 on a full neighbourhood. A c = 2 region can come out of an onset only at codimension-3 exceptions: equal masses, a zero-energy transmission resonance, or a bound state.
- Lemma P [derived]. This addresses C174's P4. At weak locking, clock tongues of every order melt through c ≤ 1. A Z_m parafermion with m ≥ 5 needs clock amplitudes of order the exchange, at codimension at least 2.
- Potts × floating boson in N [derived]. It has codimension 2, not 1, because nothing in N forbids the chiral Huse–Fisher operator. On T and K_I the antiunitary may forbid it, which is C171's case, and that remains open.
- Lemma H [derived]. c ≤ 1 holds on a full neighbourhood of the AF Heisenberg point in N, uniformly in direction.
- Residual (R_N) [assumed], on N outside the controlled set:
  - (i) inside a c = 1 liquid, no independent excitation softens, at momentum 0, π, 2π/3 or incommensurate, composites included (A_S and A_W restated for N);
  - (ii) there is no P1 corner, where a c ≤ 1 transition of a gapped phase meets the onset of an independent branch;
  - (iii) none of Lemma D's codimension-3 exceptions occurs with an easy-plane pseudospin;
  - (iv) there is no P4 fixed point at amplitudes of order the exchange;
  - (v) there is no c > 1 region whose boundary is entirely first order.

  The method stops at (i) and (ii). They need a lattice-scale branch to come down at harmonics of order the exchange, where nothing is small. A Wentzel–Bardeen-type instability argument for a dilute guest inside a dense liquid was tried and rejected, since the Hubbard chain near saturation is a stable counterexample. K_II is outside scope.
- Sampling, a correction to §9.115 and §9.125 [numerical, job 17602943].
  - C155's stage-A flags, 28 gapless points with c_obj > 1.1, were never refined. C155's refine stage reads only the search and gap-closing outputs.
  - On 1600 range-3 qubit chains with known c, C155's rule flags only 1.8 percent of near-free c = 2 chains and 0.6 percent of interacting ones. A two-stage ring screen (c_S > 1.5 at rings 12/16 and again at 16/20) flags 34 and 18 percent.
  - All 416 candidates went to rings 16/20, and 45 went to DMRG at L = 32, 64, 128. Of those 45, 43 read c_prof(128) ≤ 1.16.
  - Two points pass the first stage. idx 10193 reads 0.51, 0.61, 1.94, but its L = 128 energy lies 0.06 above the trend of the smaller sizes, so the state is probably unconverged. idx 450 reads 1.62 at L = 64 and 1.21 at 128, falling.
  - If both fail the second stage, the normalized volume of any open c > 1 set is below 7.3e-4 to 1.4e-3 at 95 percent confidence. With C155's rule alone the bound would have been about 4e-2.
  - Caveats: the proxies have range 3 and are close to free, the bound uses C155's coordinate measure, and it says nothing about small-volume sets.
  - Job 17603367 (pending) runs FES on a four-site cell to χ = 256, and DMRG reseeds at L = 128 and 192, on 10193, 450 and the control 10540.

### 9.134 Counterexample hunt by structure: no candidate, with an obstruction for each route (C169, Claude Opus 5.5 agent, 2026-11-09; `codex/work_C169/`)

- Exact maps cannot create a counterexample. c is invariant, and the preimage of B_2 is B_2 relabelled, so a map can only recognise a known c > 1 model inside B_2. None does, under any class tested.
  - Single-site unitaries, any period [proved]. They preserve the support of every Pauli term. A period-p chain lies in B_2 exactly when it has a twisted translation, which generalizes C156's screw lemma.
  - Clifford dualities: KW, JW, cluster entanglers, any Clifford QCA [proved].
    - For translation-commuting maps, the spread of anticommutation offsets is invariant. It is at most 2 in B_2 and 3 to 4 in every target model.
    - For any Clifford map, the number of independent infinite pieces of the commutation graph is invariant. It is at most 2 in B_2, 3 for cluster-Ising and Ising × XX, and 4 for two XX legs.
    - The part of B_2 that stays nearest-neighbour under KW lies on the free locus. The cluster entangler acts on B_2 as the site reflection.
  - Blocking, any block size [proved]. The coupling across a block boundary must factor as one qubit-to-qubit bond, and every target fails.
  - General QCAs carry no invariant beyond the spectrum itself. Integrable c > 1 models are excluded for symmetric J [literature, Yamaguchi–Chiba–Shiraishi]. With DM this step is assumed.
- In every controlled setting the obstruction is the nearest-neighbour range, and the bounds hold at every period. Translation by one site is never what blocks c > 1.
  - Free fermions [proved]. c ≤ 1 at any period, since two Majoranas per site with nearest-neighbour coupling allow at most four zero crossings. At range r the bound is c ≤ r, and Titvinidze–Japaridze reaches it at r = 2. This is consistent with C174's R8.
  - Hidden free fermions [proved by enumeration]. Among Pauli-aligned bonds, every one with a free-fermion-solvable commutation graph is an ordinary JW-free chain.
  - Spin waves [algebra proved, physics derived]. There is at most one soft species about any periodic classical state. The J1–J2 two-minimum route to c = 2 needs range 2.
  - Exact decoupling [proved]. At any period the only bond that splits into two commuting sectors is XZ + ZX, which is two Majorana chains with c = 1. This proves C157's uniqueness claim.
  - Ising × XX through XZ. The coupling needed is a four-site term, generated only virtually, and with the wrong sign (C157).
- The three routes from the referees.
  - Z₅ parafermion at a melting tongue (C174 R5) [derived]. B_2 has no on-site Z_m with m ≥ 3 other than U(1), so Z₅ can only be translation acting on a period-5 tongue.
    - No lattice symmetry in any stratum forbids the chiral term, so the parafermion point needs three tunings.
    - For m ≥ 5 the clock and disordering couplings cannot both be relevant in a one-boson theory, since their dimensions multiply to at least m²/4. So controlled tongues melt through a c = 1 floating phase.
    - At spin 1/2 the Z₅ clock term first appears at third order, and it is relevant only near the ferromagnet, where C160 excludes disordered melting.
    - Period-5 Rydberg chains also melt in two steps [literature]. What remains open is order-one amplitude only.
  - Potts × half-floating boson, c = 9/5 (C171 R1) [derived]. The chiral term is allowed on T and K_I as well, so R1 needs two tunings there and at least three in N. C168 independently finds codimension 2 in N. A period-3 order on one qubit per site is either the host's own locked spiral, which gaps the host, or a period-3 density wave, which no nearest-neighbour U(1) liquid has (all have K > 1/4). So the Potts sector needs a second species, and it falls under (A_W′) or (R_N).
  - Repulsive-host c = 2 (C171 R4) [derived near U(1), assumed beyond]. It needs a fluid of bound pairs with K₂ < 1/4, the n ≥ 2 gap of C163's Lemma L. In a repulsive host that branch never softens near the U(1) class, and beyond that it is (A6).
- Open: interacting nearest-neighbour chains of period p > 1, where a structured period-2 search would test whether one-site translation matters at all; bonds not aligned with a Pauli frame, in the hidden-free-fermion enumeration; the DM integrability step; and order-one harmonic amplitude, where all three routes would have to live (PROOF §9 names a Z₅ test family).

### 9.135 Typhon output of C157: the bisector line, the Ashkin–Teller wedge and the mixing directions near XZ (coordinator, Claude Opus 5.5, 2026-11-09; `codex/work_C157/typhon_out/`)

Job 17595396 completed in 22 h 33 min. The estimates below use C157's own `cscale.py` and `cfit.py`, plus central pair estimates 6 ΔS_mid / Δ ln L computed by the coordinator. Truncation is 1e-8 or below throughout.

- Bisector field, open chains at L = 192 and 256.
  - At h = 0.3 the smoothed profile fits are 0.83 and 0.81, and c_pair(192, 256) = 0.89.
  - At h = 0.45 the fits are 1.02 and 1.00, and c_pair = 1.01.
  - Window fits at h = 0.45 lie between 0.97 and 1.05 at both sizes.

  So the line is one boson out to L = 256, which agrees with the C158 FES result of §9.131 at the same fields. On the C157 estimator calibration, the known c = 3/2 RZFA phase read 1.5 to 1.7 at L = 128, so a c = 3/2 sector would have been visible.
- Ashkin–Teller wedge boundary at β = 0.5, α = 0.21 to 0.27. The pair estimates are all below 0.83 and fall with L. At α = 0.225 they go 0.74 → 0.68 over 64 → 128 → 192, and at α = 0.25 they go 0.52 → 0.38. Nothing reaches 1.
- Mixing plus plaquette, XZ + α ZX + 0.3 YY + 0.1(XX + ZZ). At α = −0.15, −0.05, 0.05, c_pair(64, 128) = 1.04, 1.06, 1.09, consistent with one boson at these sizes. α = 0.19 is gapped (S_mid = 0.068 at both sizes).
- Verdict [numerical]: no c > 1 reading near XZ in C157's directions. This closes the numerical side of C157's argument that the XZ loophole has the wrong sign. The mixing points were taken only to L = 128.

### 9.136 Typhon output of C167: the hi_K_1 harmonic path to χ = 192 (coordinator, Claude Opus 5.5, 2026-11-09; `codex/work_C167/typhon_out/`)

The job ran the hi_K_1 path (harmonics scaled by λ) on the infinite chain with two-site VUMPS to χ = 192 (256 at λ = 1), a four-site cell at λ = 1, the charge-1-only and charge-2-only bonds to χ = 128, and open chains of hi_K_1 at L = 128 unpinned and at L = 128 and 256 with edge fields along the VUMPS order (truncation below 5e-10). The digest is `summary.txt`. C167 registered three predictions before the output (§9.129).

- x_π (last two χ pairs) falls continuously along the path, 0.97, 0.80, 0.73, 0.64, 0.58 and 0.49 at λ = 0.90, 0.95, 0.96, 0.97, 0.975 and 0.98, then 0.30, 0.08 and 0.01 at λ = 0.985, 0.99 and 0.995. Interpolating, it crosses 1/2 at λ ≈ 0.979. The registered window was 0.98 to 0.985, so this part of the prediction misses by about 0.001 in λ. A finite-χ estimate of x_π at a KT point carries logarithmic drift of at least this size, and the qualitative claim (a continuous fall through 1/2, then order) holds.
- c_pair at the last χ pair is at most 1.03 everywhere on the path (λ = 0.98), and drops to 0.89, 0.36, 0.06 and 0.00 as the order sets in. The prediction c_pair ≤ about 1.05 holds. The staggered moment converges in χ from λ ≈ 0.99 on (0.183 at λ = 1, 0.243 and 0.277 at λ = 1.01 and 1.02).
- At λ = 1 the four-site cell reaches the same e₀ = −0.65620519 as the two-site cell, so there is no lower period-4 state. The charge-1-only and charge-2-only bonds at λ = 1 are c = 1 liquids (c_pair 1.00 at χ = 128 for both, ξ still growing).
- Pinned open chains. S_mid = 0.786 at L = 128 and 0.796 at L = 256, against the VUMPS bond entropies 0.776 and 0.819 (mid-chain pair at L = 256, 0.818 and 0.775). c_prof drops to 0.34 and 0.13, and the mid-chain bond energy at L = 256 is −0.656206, within 1e-6 of the VUMPS e₀. The prediction (within 0.05 of VUMPS) holds. The unpinned L = 128 chain reads c_prof 1.26 and S_mid 1.01, the slowly relaxing cat already seen at XZ (§9.126).
- Verdict [numerical]. hi_K_1 is a gapped period-2 state reached through the boson's own KT, as C167 concluded from χ ≤ 96, and the window point that read c above 1 on open chains is closed. The residual (A_W1) is unchanged. Two of the three registered predictions hold, and the third misses its window narrowly.

### 9.137 C173: the zigzag-ladder route, an obstruction in most regimes and a candidate c = 2 corner near the Fibonacci point (Claude Opus 5.5 agent, 2026-11-08 to 09; `codex/work_C173/`)

A nearest-neighbour qubit chain is a zigzag ladder (even and odd sites as legs) with rung couplings only, and two floating bosons would mean one band with two occupied valleys. C173 asks when the generated intra-leg dynamics can make two valleys.

- Obstruction [derived]. One valley, hence no second single-particle species, (i) whenever the first-order hop dominates the second-order ones it generates (Lemma S [proved], |t₁| > 8|t₂| + 2Σr²|t_r|), (ii) at large field (gapped), (iii) in the classical spiral, where the rungs lock the relative leg phase at first order, and (iv) in the dilute part of the Fibonacci corner, where inter-valley scattering is impenetrable and a second valley gains no kinetic energy, except in a wedge |δ| ≲ n³ around the J_x + J_y = 0 sub-stratum. The controlled weak-rung points are large field, XZ (§9.135) and the Fibonacci point.
- The Fibonacci cusp [derived, checked]. At J ≈ J_z ZZ, B ≈ 2J_z ẑ, with the first-order hop J_x + J_y + 2iD tuned to order Δ²/V (Δ = J_x − J_y, V = 4J_z), the pairing acts only virtually across the blockade and generates a leg hop t₂ = −Δ²/(2V) with leg interactions V₂ = 2t₂, V₃ = (2/3)t₂. The effective model matches the full L = 12 chain to 1.2e-2, 1.5e-3, 1.9e-4 at V = 20, 40, 80 (error ∝ V⁻³). At w = 0 it is the staggered-charge U(1) case with c = 1 (C153). C152 §7 assumed |t| ≫ Δ²/V, and the single-species Fibonacci results of C141 and C143 treat PXP and U(1) hopping, so neither covers this corner.
- Numerics [numerical, not converged]. Effective-model DMRG, two-size estimator on the mean of the three middle cuts, controls 0.90 (c = 1) and 1.88 (c = 2) for L = 48 → 96. Away from φ ≈ π/2 every point is c ≤ 1. At φ = π/2, n = 0.17, the mid-chain entropy is 1.646, 1.790, 1.869 at L = 48, 72, 96, giving 2.13 (48 → 72) and 1.93 (48 → 96). The same reading holds at φ = 1.56, it is high at φ = 1.54 (L = 48 only), and one valley by φ = 1.5. The truncation at L = 72 and 96 is 5e-7 and 4e-6, which fails the first stage of C163 §4.2, and the c = 2 control at L = 96 is equally unconverged. I re-read the jsonl outputs. The estimator moves between 1.8 and 2.1 with the choice of cut window, because the incommensurate Friedel oscillations are large at these sizes.
- Candidate bond. h = 0.5 XX − 0.5 YY + 0.00625 (XY − YX) + 10 ZZ − 9.988125 (Z⊗1 + 1⊗Z), energy scale |t₂| = 0.0125, on the J_x + J_y = 0 sub-stratum of the parity class. There an exact symmetry (staggered rotation × inversion, C173 PROOF 3.3) maps one valley onto the other. C173 argues the region is open in B₂ if the reading holds, since the parity-odd terms couple only to off-resonant vertices [derived, conditional].
- Coordinator's note. Physically this is the question familiar from frustrated J1–J2 chains near saturation, where two-valley magnon gases either choose one valley (a vector-chiral phase, c = 1) or stay two-component (c = 2), depending on the sign of the inter-valley interaction. That sign is the same one as C162's immiscibility argument. C162 found immiscibility in its corner, and ring ED here shows a valley-balanced ground state. At an exactly symmetric point that is also what a symmetric superposition of a valley-polarized pair looks like, so the analytic sign (a three-body computation in the blockaded model) is the deciding input alongside the Typhon run.
- On C171's routes. R3 (composites) is not covered in any regime. R4 is covered where the second boson would be a single-particle valley, in the strong-rung regime and in the dilute cusp off the wedge. The commensurate cusp (φ = 0) reads c = 1.15 at n = 0.25, one boson.
- Still open. A converged reading (L to 192 at χ = 192 for φ = π/2, 1.54, 1.52, both controls at L = 144, and the full qubit chain at V/Δ = 40 and 80) is on Typhon and pending. Then FES at χ ≥ 256, the pseudospin gap on rings against 1/L, and the analytic pseudospin sign. If the reading holds, C152 Theorem C(e), the Fibonacci single-species claims and (H)/(A6) of C154 fail in this corner.

### 9.138 C180: field theory and literature for the Fibonacci-cusp candidate (Claude Opus 5.5 agent, 2026-11-09; `codex/work_C180/`)

C180 bosonized the two-valley gas of C173's effective model and surveyed the literature. Verdict: no symmetry or kinematic obstruction to c = 2. The question reduces to whether the valley pseudospin is easy-plane once the generated attractions V₂ and V₃ are included, which is C175's computation.

- Symmetries [proved]. At φ = π/2 the gauge c_j = i^j d_j turns the effective model into a real J1–J2 chain with J1/J2 = |w|, plus the blockade, an attractive J2^z = −J2 and a third-neighbour attraction. The valleys sit at π/2 ∓ k₀ with cos k₀ = −|w|/4, so k_a + k_b = π. Inversion × complex conjugation is an exact antiunitary of every parity-class bond and keeps each valley. The valley exchange Σ maps (J_x, J_y, D) → (−J_y, −J_x, D), so it is exact only on J_x + J_y = 0. I checked both by hand.
- Operator table [derived, enumeration checked by `optable.py`]. At the symmetric point the only momentum-zero vertices are inter-valley backscattering cos 2√2φ₋ (dimension 2K₋, relevant if K₋ < 1) and, in the parity class only, the two-pair vertex e^{2√2iθ₊} at momentum 2(k_a + k_b) = 2π (dimension 2/K₊, relevant if K₊ > 1). Off the symmetric point the valley offset is a pseudospin field and no vertex survives. All other vertices are resonant only on codimension-one lines and cut c ≤ 1 tongues, one of which (a dressed pair transfer, dimension 2/K₋ + 2K₊) falls at n ≈ 0.161 for |w| = 1, near the candidate's n ≈ 0.17.
- Stability of c = 2 [derived]. (S1) a positive inverse valley susceptibility including the drift terms allowed by Σ, whose failure is the vector-chiral instability, (S2) K₋ ≥ 1 (marginal on the irrelevant side), (S3) K₊ < 1 in the parity class, (S4) off the resonance tongues. Off the symmetric point c = 2 survives for |δ| of order n³, which reproduces C173's wedge. Off the parity class the region is open in B₂ when (S1) to (S4) hold with margin, protected by incommensurability.
- Dilute limit [derived]. The impenetrable gas maps to an XXZ pseudospin chain with anisotropy set by the same-valley excluded length a_s and the inter-valley length matrix (A_d, A_o). c = 2 if A_d − |A_o| < a_s < A_d + |A_o|, vector chiral below, valley crystal above. A δ-interacting gas reproduces the standard exchange 4π²n³/(3m²g).
- Competing c = 1 phases. Vector chiral (Σ broken), valley-crystal liquid, ab or aa bound-state liquids, commensurate charge locking, and ΔN = 4 locking of θ₊.
- Literature (arXiv numbers opened by C180). Mishra, Greschner and Santos (arXiv:1411.7270) study the hard-core frustrated triangular ladder, which at V → ∞ and t′ = −t is the effective model at |w| = 1 without V₂ and V₃. They find the infinite-blockade limit on the two-component side, no chiral region in the squeezed model, and c = 2 for ρ < 0.4 at V = W = 8. Kolezhuk et al. (arXiv:1112.5720) cover 0 ≤ J2^z/J2 ≤ 1 only, so the candidate's J2^z = −J2 and third-neighbour attraction are outside every published calculation. A Born-level estimate, not reliable at hard core, says both attractions push toward vector chiral. Verresen, Vishwanath and Pollmann (arXiv:1903.09179) confirm t₂ = −Δ²/(2V) and V₃ = −Δ²/(3V) independently in the Rydberg-constrained Kitaev chain, which has real hopping and so no valley degeneracy. No integrable match was found (Corcoran, de Leeuw and Pozsgay, arXiv:2405.15848, partial list).
- Known or excluded? Neither. Every published c = 2 spin-1/2 phase found has range at least 2 or is a ladder, and no theorem reaches the Fibonacci corner. If confirmed, the candidate gives c = 2 at d = 2, against Latorre and Sierra's c ≤ d − 1 by one unit. The replacement bound is not affected.
- Registered predictions for the other agents (before their output). Two velocities with u₊ ≈ 4|t₂| and u₋ well below, a valley-imbalance gap scaling as 1/L, K₊ between 0.3 and 0.5, two equal n(k) peaks, a possible c = 1 dip near n ≈ 0.161 at |w| = 1, and the V₂ = V₃ = 0 limit two-component.
- Open. The pseudospin sign with V₂ and V₃ (C175). A possible near-threshold resonance at φ = π/2 (C173 saw two-body levels 0.002 to 0.008 above threshold at φ = 1.3), which would break the dilute expansion. The size of the drift correction at n = 0.17, which C180 estimates as comparable to the susceptibility itself.

### 9.139 C175: the two valleys bind into pairs at the symmetric point (Claude Opus 5.5 agent, 2026-11-09; `codex/work_C175/`)

C175 computed the valley-pseudospin sign of the Fibonacci-cusp gas (§9.137) from the exact two-body problem of C173's effective model. The dilute limit is not two-component. One particle from each valley binds into a pair, the dilute phase is a pair liquid with a pseudospin gap, and c = 1 there. At the candidate density n ≈ 0.17 the question stays open.

- Criterion [derived]. The dilute impenetrable gas maps to an XXZ pseudospin chain, and it is c = 2 iff the same-valley length a_s lies inside (A_d − |A_o|, A_d + |A_o|), with A the inter-valley length matrix. This agrees with C180's criterion, reached independently. The lattice matrix has a non-Hermitian part from the band's cubic term, independent of V₂ and V₃, which drops out.
- Closed form [proved, matches numerics]. At φ = π/2 with u = |w|, A_d = 5/2 + 8/u², with eigenvalues a_∓ = [5u² + 16 ∓ √(u⁴ + 160u² + 256)]/(2u²).
- Exact resonance [proved]. The mixed-valley scattering length diverges on 2V₂V₃ + 4V₂ + V₃u² + 4V₃ + 8 = 0. At V₃ = 0 this is V₂ = 2t₂ for every |w|, which is the Schrieffer–Wolff value. At the Schrieffer–Wolff point (V₂, V₃) = (2t₂, (2/3)t₂) the left side is −(2/3)u², on the bound side, and higher orders (relative size (Δ/V)²) cannot undo it [derived]. I checked the substitution by hand.
- Bound pair [numerical, two independent codes in C175]. It sits at total momentum π, with E_b = 0.0017, 0.0195 and 0.1225 at |w| = 0.5, 1 and 2 (E_b ≈ |w|⁴/32 at small |w|), and is about 21 sites across at |w| = 1. Trimers bind at |w| = 1 and 2, and larger clusters at |w| = 2.
- Coordinator check [numerical, a separate two-particle exact diagonalization on rings of 160 and 240 sites, written without reference to C175's code]. E₂ − 2ε_min at |w| = 1 is −0.0217 and −0.0202 (L = 160, 240) with both attractions, against +0.0007 with V₂ = V₃ = 0 and a value that changes sign with L (−0.0007, +0.0003) for V₂ alone. At |w| = 2 the numbers are −0.1236 and −0.1226. Both match C175's binding energies, and the V₂-only channel behaves as a zero-energy resonance should.
- Calibration [numerical]. With V₂ = V₃ = 0 the gas is two-component (at |w| = 1, a_s = 0.629 inside (−0.016, 1.016)), which reproduces Mishra, Greschner and Santos. With both attractions a_s = −3.17, below the window, and the scattering states alone would be vector chiral. The bound pair dominates both.
- Dilute window [derived at leading order in n]. Pairs win while the valley offset is below E_b (|φ − π/2| ≲ 0.005 at |w| = 1), and a single valley beyond. The dimer-plus-atom mixture between is not two-component because trimers bind (this step also uses C164's mass-imbalance result). The c = 2 window at small n is empty, replacing C173's wedge.
- Candidate density [numerical and heuristic]. Ring ED (L = 12 to 36, momentum resolved) puts every even-N ground state at n ≥ 0.13 in the valley-balanced pair sector at K = Nπ/2, with valley-polarized sectors 0.46 to 0.66 higher, so vector chiral is ruled out there. At n = 0.17 and |w| = 1 the pairs overlap strongly (n·a₊ ≈ 3.5), and a Yang–Gaudin-type estimate gives a pair gap of about 1e-4 |t₂| and a correlation length of order 1e4 sites. That would make C173's c ≈ 2 at L ≤ 96 a crossover. The attractive same-valley channel pushes the other way, so a gapless window at intermediate density is not excluded.
- Against its registration. P1 fails through its no-bound-state condition, so by C175's pre-fixed criterion the dilute gas is not c = 2. P2(ii) and P2(iii) pass. P2(i) failed as stated (the non-Hermitian part above). P3 was not testable as registered.
- Open. The phase at |w| = 1, n ≈ 0.17. The tests are the parity gap E(N+1) + E(N−1) − 2E(N) against L, pair correlations at momentum π against single-particle ones, and c at lower density or larger |w|, where the pair gap is large. A Typhon run of these is pending, with the registered prediction c(96 → 192) < 1.4 at the pairing points and > 1.6 for the V₂ = V₃ = 0 control. Also open are clusters beyond trimers at |w| = 1, and whether the full qubit chain carries the bound pair (C177).

### 9.140 Independent continuation: exact cusp binding, a local reduction correction, and a stronger d=13 pin bound (OpenAI Codex, GPT-6 family, 2026-11-09)

The [continuation](contrib/counterexample_push_20261109/README.md) began at
`a45ecf5` and pulled C175 (`1a4015f`) while independent calculations were
running. It credits C175 for the binding discovery and adds the following.

- [proved] **Exact binding boundary** in the symmetric second-order
  blockaded model. In the real gauge, at total momentum zero, take hops
  u,+1 and potentials V₂=-2+v, V₃=-g, with 0<u<4 and 0<g<2. There is a
  bound state below the global free-pair threshold exactly when
  v<g u²/[2(2-g)]. The identity H-E_th=2(A+u/2)²+vP₁-gP₂ for half-line
  adjacency A, followed by completing the boundary square, proves both
  directions. At v=0 rank one gives uniqueness. At u=1,g=2/3 the
  threshold is v=1/4. An exact integer-vector variational certificate
  gives binding >1/100. [Proof and reproduction](contrib/cusp_binding_certificate_20261109/README.md).
- [proved algebra; numerical evaluation] **Explicit pair state.** The
  unique p in (0,1) solving (1-p)(3-p)(1+p)²=u²p³ gives
  E_pair=-2[p+1/p+u²p/(1+p)²], a normalized wavefunction and its momentum
  curvature. At u=1, E_pair=-4.5195265346664 and inverse pair mass is
  3.373597375632. Independent cutoff diagonalizations at seven u values
  agree to 7e-11 or better. A unit-probability-flux scattering calculation
  independently checks C175's Hermitian length-matrix convention at
  three parameter choices. [Derivation and data](contrib/cusp_scattering_20261109/README.md).
- [proved local coefficient; numerical matrix checks] **Scope correction
  to §9.139's robustness estimate.** In H=V D+Delta T+(Delta²/V)W,
  the canonical effective Hamiltonian contains
  (2Delta³/3V²)P sum_j(b†_j b†_(j+3)+h.c.)P. Its coefficient relative to
  |t₂| is 4Delta/(3V), or 1/30 at V/Delta=40. Thus the full local
  remainder is not only relative order (Delta/V)². Direct pairing into
  the valley minima is off resonance, so this does not overturn the
  leading binding result; bounding its higher-order effect needs an
  explicit elimination. Forty-five spin-basis coefficient checks
  through L=12 pass. The short-ring spectral check is vacuum-side even
  at the candidate couplings and is not finite-density validation.
  [Proof, checks and limitations](contrib/cusp_validation_20261109/README.md).
- [proved under C92's stated hypotheses] **Stronger mobile contact and
  pin bounds.** The vacancy-removal channel has entropy loss at most
  h₂(a), where a is the contact probability, improving C92's log 2.
  In its unchanged legal-space positive-tilt window this gives a>1/264
  (and a≥1/257 at zero temperature), and pin activation at least J/135
  instead of J/4096. Exact arithmetic, 72 quantum-channel examples and
  seven saturating examples pass. The d=11/d=13 bond reconstruction
  also reruns to roundoff. [Proof and route ledger](contrib/route_audit_20261109/README.md).

The coordinating Codex agent reviewed the proofs and reran the checks;
two parallel Codex agents independently reviewed the binding theorem.
No human review is claimed. The new two-body results do not imply a
particular many-body paired phase or settle c at density 0.17. The
stronger pin estimate is still a legal projected tilted-ensemble result,
not the missing finite-K or wrong-length comparison. The d=13 mobile
construction remains the main analytic priority of this continuation.

### 9.141 A new transverse-DM qubit family with a binding-reducing direction (OpenAI Codex, GPT-6 family, 2026-11-09)

The [constructive follow-up](contrib/counterexample_push_20261109/FACILITATED_ROUTE.md)
uses the exact binding boundary of §9.140 to seek a microscopic escape
from C175's obstruction. This is a new test family, not a counterexample.

- [proved local reduction; independent exact enumeration] Add
  lambda sum_j(n_j X_(j+1)-X_j n_(j+1)), which on a periodic ring is
  exactly the nearest-neighbour transverse DM term
  (lambda/2) sum_j(X_j Z_(j+1)-Z_j X_(j+1)). It vanishes on compression
  to the blockade. At second order it generates a distance-two
  repulsion 2lambda²/V, a real nearest hop lambda²/V and a number
  shift. Crucially, its mixed term with the original pairing is
  (Delta lambda/V)P sum_j(n_(j+2)-n_(j-2))X_j P. Three independent
  implementations agree on the full operator. The repulsion-only
  shortcut omits a term of the same effective scale and is invalid.
- [derived and numerical] With r=lambda/Delta, compensate the direct
  real hop and retain the complete leading model. The coefficient of r
  maps an atom at momentum K into a distance-two pair with amplitude
  4i sin K. At the original pair momentum pi its two-to-one channel
  vanishes; its two-to-three channel remains. Exact pair wavefunctions
  and interacting two-/three-body resolvents give, at mu_eff=-2.27,
  atomic r² shifts -0.825224871396 and -2.663400053460 and pair shift
  -1.163523967762. The three-body gap-coordinate calculation reaches
  dimension 57,121; cutoff changes in the pair coefficient are below
  4.1e-7. These are isolated-branch perturbation coefficients, not a
  proof of the vacuum or its thermodynamic stability.
- [derived] The perturbation breaks valley exchange. An additional
  effective real hop q r², q=-0.474614791170, balances the dressed
  valley energies through this order. The pair at pi has zero linear
  response to that hop. The resulting formal binding expansion is
  B(r)=0.0195265346664-2.32510095709 r²+O(r⁴).
  The apparent zero r≈0.0916 is an uncontrolled extrapolation; its
  fourth-order remainder and branch isolation have not been bounded.
- [numerical, finite rings] Exact constrained-spin diagonalization
  retains every particle number and periodic spin boundary conditions.
  At L=24, mu_eff=-2.30 and q=-0.36413 (the retuning for that source,
  rounded), the pi-pair energy minus two tracked valley-like energies
  is -0.041485743 at r=0, -0.019439073 at .10 and +0.002211339 at .15.
  A nearby tuning gave a much larger positive value on L=20, so size
  drift is substantial. These discrete levels do not certify an
  infinite-volume unbinding transition or c=2.
- [checked finite algebra] The contribution includes a generator and
  16 literal Hermitian 4-by-4 qubit bonds at V=40,80, with the real-hop
  and field tunings spelled out. Their local occupation-basis identity
  is checked up to a scalar and a telescoping one-site term. They are
  test points with no assigned central charge.

Independent reviews and reproducible data are in
[the reduction audit](contrib/cusp_validation_20261109/FACILITATED_REPAIR.md),
[the pair calculation](contrib/cusp_scattering_20261109/FACILITATED_EXTENSION.md),
[the atom calculation](contrib/cusp_binding_certificate_20261109/facilitated_atom.py),
and [the independent review](contrib/route_audit_20261109/FACILITATED_REVIEW.md).
The next direct search should track the dressed bound-state pole and
continuum in this complete model, retune the valley offset, and check
larger clusters before attempting a central-charge claim. The d=13
construction remains the mature analytic route; this new qubit family
is the continuation's preferred next numerical test.

### 9.142 C176: a blockade-exact DMRG for the cusp model, and preliminary readings (Claude Opus 5.5 agent, 2026-11-09; `codex/work_C176/`)

C176 built and validated an open-chain DMRG for C173's effective model in which the blockade and the particle number are exact. Its production runs (L to 256, bond dimension to 1024, five densities, φ = π/2 to 1.50, both controls, N ± 1 for the parity gap) are queued on Typhon, and nothing beyond L = 64 has been read.

- Construction [proved]. Left-canonical tensors carry (N_left, n_{i−1}) and right-canonical tensors (N_right, n_i), so both bases are orthonormal and the constraints are structural zeros, with no penalty. Labelling both sides by the occupation to the left gives non-orthogonal right bases and wrong Schmidt values.
- Validation [numerical]. Against C173's open-chain ED at L = 12 to 24 over ten cases: energies to 3e-12, S(l) to 6e-11, correlators to 1e-9, ⟨n_j n_{j+1}⟩ = 0 exactly. The open-chain bulk energy at n = 1/6 agrees with ring ED within the ring's shell oscillation.
- Densities [numerical, derived]. The band bottom is −9/4, and open-chain edges hold about two fewer particles, so the bulk density is about N/L + 2/L. At fixed N/L this drift pulls two-size estimates down by about 0.1 for 128 → 256, which makes a pass conservative. C176 added a secondary fixed-bulk-density estimator, not registered.
- Preliminary [numerical, L = 48 → 64, bare Schrieffer–Wolff couplings, truncation 1e-8 to 4e-8, about ±0.3]. Two-size c is 0.98 for the c = 1 control, 1.92 for the c = 2 control and 1.79 at n = 1/8 (2.26 at a fixed bulk density of 0.17). The valley-polarization correlation decays. L times the parity gap is 2.70 and 2.40 at L = 48 and 64. These are at the couplings where the exact two-body problem has a bound pair (§9.139), so they show that sizes this small cannot see the pairing. The deciding outputs are at n = 1/8 and L = 256 (parity gap and pair against single-particle correlations).
- Registered before production. The first-stage criterion of C173 PROOF 5.4 at n = 3/16 and L = 128, control bands of ±0.15 and ±0.25, predictions for 128 → 256, and, after C175, predictions for the parity gap and pair correlations.
- Open. The production job. Pair gaps of order 1e-4 are invisible at L ≤ 256. Open chains only, so no FES stage.

### 9.143 C184: multicritical points with c > 1, an inventory (Claude Opus 5.5 agent, 2026-11-09; `codex/work_C184/`)

A fine-tuned point counts as a counterexample, so C184 inventoried the ways a nearest-neighbour qubit bond could reach a multicritical point with c > 1. It found none, and no candidate bond. Its registered open-chain test of the Peotta meeting region is on Typhon and pending.

- Closure lemma [derived, agrees with C174]. It does not exclude multicritical c > 1, because bounds pass only from higher to lower c. Only (H*)/(A6) excludes such points.
- Lemma M [proved for quadratic chains, derived for weak interactions; inputs C152 Prop. F, Lemma P2, C174 R8]. Near any quadratic reference a nearest-neighbour qubit chain has at most four soft Majorana–Weyl modes at any crossing, while every c > 1 candidate needs at least six or a strongly interacting species. So a c > 1 crossing can sit only where a non-quadratic species exists: the XZ flat band, the Fibonacci cusp, strong attraction, or near the ferromagnet. Crossing two transitions changes the number of tunings, not the requirement of two independent soft sectors.
- Inventory [derived]:
  - Ising × Gaussian needs three tunings at a commensurate Gaussian surface, because the σ·V couplings are relevant. The exception is K > 15/8 with σ even under the π rotation, which is C167's (A_W1)(a).
  - TCI × Ising needs four tunings perturbatively. In the parity class the only inter-species term is also the TCI-driving self-interaction, so the crossing flows away (heuristic at V ~ t).
  - Ising² × boson needs eight Majorana–Weyl modes.
  - Potts × boson and Z₅ reduce to C168, C169 and C171, and SU(2)₂ reduces to P1.
  - Two Gaussian surfaces give c = 1 when they belong to one boson. Only the cusp supplies two bosons.
  - Lifshitz points are not conformal.
- Two old jobs, never collected until now [numerical]. In C152's job, the worst XZ-window point reads c_pair 1.46, 1.24, 1.13 over L = 32 to 192, falling toward 1, the Peotta meeting region reads at most 1.15, and at V = −1.5t everything is gapped. In C154's job, the SU(2)₁/ordered boundary of the cyclic stratum looks like a level crossing, with the Casimir c falling from 0.93 to 0.69 and no SU(2)₃ tower.
- Ring ED of the Peotta region [numerical, φ = 1, h_z = 0, L ≤ 20]. The DM-split Ising band reads c_S(16, 20) = 0.96 to 0.99, and its meeting with the Luttinger region is a level crossing.
- Estimator calibration [numerical, exact free fermions]. c_pair at (48, 96) and (96, 192) reads 0.99 for XX, 0.51 for TFI, 1.98 to 1.99 for XX × XX and 1.50 for TFI × XX.
- Open. The registered DMRG test. The K > 15/8 window. A fine scan of the cyclic SU(2)₁/ordered boundary for SU(2)₃, the one place where c > 1 needs a single tuning inside a stratum and nothing excludes it dynamically, which was not run. TCI × Ising at V ~ t.

### 9.144 C185: which bonds have an exponentially degenerate ground space, and the transverse-DM route at the Fibonacci point (Claude Opus 5.5 agent, 2026-11-09; `codex/work_C185/`)

The Fibonacci cusp works because a dominant bond has an extensively degenerate ground space in which second-order dynamics can win. C185 looked for every such point. It found no closure and no confirmed c > 1 candidate, plus one new route whose test is queued.

- Classification. Only two orbits in B₂ have an exponentially degenerate ground space: the Fibonacci orbit (J ZZ + bZ at |b| = 2J, J > 0, the rank-1 projector onto |aa⟩) and XZ.
  - For classical bonds this is [proved].
  - For frustration-free bonds it is [proved for |aa⟩, numerical otherwise]: every other rank-1 kernel gives exactly L + 1 states, and ranks 2 and 3 give at most 2, checked exactly to L = 14.
  - Among bonds with local integrals of motion, only XZ qualifies (C143).
  - Points with polynomial degeneracy, such as the ferromagnetic multiplet, cannot host the mechanism [derived], because their degenerate space has zero entropy density and a single-minimum magnon band.
  - XZ stays closed, because its induced interaction is attractive in every channel (C144, C157).
- First order at the Fibonacci point [proved]. Modulo rotations and scale, the directions invisible at first order are the pairing Δ ∈ ℂ and the transverse DM terms λ_y(XZ − ZX) and λ_x(YZ − ZY). In the parity class only Δ survives, so at second order every parity-class bond reduces to C173's model, where the pair is bound (§9.139, §9.140). Unbinding inside the parity class therefore needs higher orders, which agrees with the V/Δ = 8 fit of C177.
- Transverse DM at second order [derived, matches the exact full-chain extraction to 1e-4]. In units |t₂| the term gives:
  - a repulsive V₂ = −2 + v with v = 16λ²/Δ²;
  - a real NN hop 4λ²/V;
  - nothing at distance 3;
  - a conditioned ΔN = 1 creation κ Σ σ⁻_s(n_{s+2} − n_{s−2}) with κ = 4λ/Δ. The uniform part cancels exactly, and v = κ².
  This agrees with the collaborator's transverse-DM family (§9.141) and with C186's general form.
- Number-conserving part [numerical]. For v above the exact unbinding boundary u²/4 the pair is unbound. The dilute pseudospin is easy-plane (Δ_σ = 0.27 to 0.95 at u = 0.5 to 2), and no trimer binds at u = 1. This is the first setting near the cusp that is two-component by C175's criterion.
- The catch [numerical, ring ED at L = 20 and 24]. κ is of order one wherever the pair is unbound. At κ = 1 it breaks N and the valley symmetry, moves the valley-balance point to φ* ≈ 1.64, and raises the density at fixed μ. Local open chains at L = 36 carry no information on c.
- Candidate bond (not established), at V/Δ = 40 in the convention of §9.137 plus λ(XZ − ZX): J_z = 10, J_x − J_y = 1, λ = 0.3, J_x + J_y = −0.009 − 0.0125 sin δφ, D = 0.00625 cos δφ with δφ ≈ 0.07 to 0.1, h_z = −19.98525.
- Pending tests. A small backfill job runs the κ = 1 φ scan around φ* at L = 48 and 96, a κ = 0, v = 1 reference, and two full-chain runs at V/Δ = 8. A larger job runs the number-conserving part, more full-chain φ shifts, and an unbiased search for further degenerate bonds.

### 9.145 C181: the partially polarized pair route off the symmetric line (Claude Opus 5.5 agent, 2026-11-09; `codex/work_C181/`)

Off φ = π/2 the valley offset acts as a field on the bound pairs. In the attractive Gaudin–Yang gas the corresponding phase, partially polarized, is a two-component liquid with c = 2. C181 tested whether that survives here. The route is neither closed nor established, and the first numerics lean against it at the sizes reached. The registered test is queued.

- Mechanism [derived]. In any imbalanced two-valley liquid the backscattering term that opens the pair gap oscillates (the 1D FFLO mechanism), so c = 2 is generic where such a liquid is homogeneous. The exceptions are composite plateaus (m = 0, n/3, n), resonance tongues and first-order jumps. The FFLO pair momentum is K_p ± π(ρ_a − ρ_b).
- Few-body data.
  - The offset is δ = 4u√(1 − u²/16)(π/2 − φ) [proved], with a mass imbalance of only δ/(16 − u²).
  - Pair binding is constant in φ to 1.5% down to φ = 1.48 [numerical].
  - The aab trimer stays bound at every φ [numerical]. At |w| = 2 its binding falls from 0.0189 to 0.0109 over φ = π/2 to 1.50. At |w| = 1 it is about 6e-4, about 95 sites across.
  - Which species appears first from the vacuum (pairs, trimers or atoms) follows from E_b and the trimer binding B₃ [derived], and dilute mixtures phase-separate by C164's mass-imbalance coefficient.
  - So there is no c = 2 window at small density, as C175 found, but this holds only for n below about 0.01 at |w| = 2 and a few times 0.01 at |w| = 1.
- |w| = 2 at its liquid density (bulk n ≈ 0.35) [numerical]. The open chain is a self-bound cluster liquid with short-ranged single-particle and pair correlations, and the 28- and 32-site rings prefer the trimer composition. Composite locking seems to pre-empt c = 2 there, and a partial 36-site ring makes this less clear.
- |w| = 1, n = 3/16 [numerical].
  - The valley imbalance grows continuously with the offset, and there is no parity gap at L = 64.
  - At φ = 1.53 the two-size estimate for L = 64 → 96 reads 0.81, 0.99 and 1.15 over the three cut windows.
  - This is not decisive. These chains hold only one to three unpaired particles.
  - At these sizes the fixed-window estimator is uncalibrated for c = 2: the c = 2 control reads 1.54 to 3.46. A window proportional to L, chosen after seeing the controls and so not registered, reads 1.00 and 1.79 on the controls and 1.06 at φ = 1.53.
- Vertex check off the stratum [derived]. No vertex is resonant at generic points. The two-pair vertex O₄ is detuned only by q₄ = 2u²ϑ/(8 − u²/2), with ϑ = π/2 − φ, which is 0.004 to 0.008 at |w| = 1. Openness needs q₄ξ₄ ≫ 1, which holds at V/Δ = 40 and fails at V/Δ = 8.
- Best bond to test, in the parity class off J_x + J_y = 0 (|w| = 1, φ = 1.555, V/Δ = 40, μ_eff ≈ −2.143): J_x = 0.500098723, J_y = −0.499901277, D = 0.006249220, J_z = 10, h_z = −19.97589125. This is not a c > 1 claim.
- Open. The registered test M1 to M6 (L = 64 → 128 with N ± 1, both |w|). Larger offsets at L ≥ 128. FES, ring spectroscopy and the full chain at this bond. The reality of the |w| = 2 trimer plateau. The FFLO momentum shift. The locking length of O₄ in the imbalanced liquid.

### 9.146 C186: unbinding the cusp pair by changing the bond, and where the V/Δ = 8 unbinding comes from (Claude Opus 5.5 agent, 2026-11-09; `codex/work_C186/`)

The collaborator's exact boundary (§9.140) says the cusp pair is bound unless v = V₂ + 2 exceeds g u²/[2(2 − g)] in units |t₂|. C186 asked whether any nearest-neighbour bond moves the generated couplings across it.

- Second order, every NN bond [derived, full-chain check to 3.6e-15]. Only the pairing generates the NNN hop and V₃, so g = 2/3 is fixed. Symmetric conditioned-creation terms (XZ + ZX, YZ + ZY, transverse fields) only renormalize the NN hop, the on-site energy and a uniform ΔN = 1 term, which the field cancels. The antisymmetric part, a transverse DM term Γ(XZ − ZX), gives v = |d|²/|p|² ≥ 0. In the parity class v = 0. This agrees with C185 (§9.144).
- The price [derived]. v > 0 always comes with an uncancellable vertex √v |t₂| Σ P b_j†(n_{j+2} − n_{j−2}) P + h.c., which also breaks the valley symmetry Σ, so φ must be retuned.
- What the vertex does [numerical, exact ring ED at L = 48 converged in particle number to 6, calibrated on the collaborator's E_b]. It weakens the binding. At u = 1.035 the pair is unbound within resolution for v ≥ 0.1, while the model without the vertex has its boundary at 0.268 and is still bound by 0.0088 and 0.0015 at v = 0.1 and 0.2. C186's registered prediction, that the vertex rebinds the pair, failed. After retuning φ two degenerate dressed valleys survive, with no pair or cluster below the single-quasiparticle onset (v = 0.1, 0.3). Without the vertex, past its boundary, the dilute pseudospin is easy-plane (Δ_σ = 0.60 to 0.95), with no trimer. The pseudospin with the vertex was not computed.
- Candidate bond (an unbound two-valley gas, not a demonstrated c = 2 point), V/Δ = 40, u = 1.035, v = 0.3, φ = 1.5881, in h = J_x XX + J_y YY + D(XY − YX) + J_z ZZ + (h_z/2)(Z⊗1 + 1⊗Z) + Γ_xz XZ + Γ_zx ZX: J_x = 0.49895031, J_y = −0.50104969, D = 0.0064695034, J_z = 10, h_z = −19.978903, Γ_xz = −Γ_zx = 0.13693064 (= Δ√v/4). For u = 2 the vertex removes the two valleys, so the window is u ≲ 1.5.
- The V/Δ = 8 conflict, resolved [numerical, derived]. C177 fitted V₂ = −1.733 at V/Δ = 8 (v ≈ 20(Δ/V)², pair unbound), while C186's des Cloizeaux extraction gives v ≈ 3.5(Δ/V)². The difference is the chemical potential. C177's fit was taken at μ_eff = −12, and C186 reproduces it there (V₂ = −1.722, −1.917, −1.986 at V/Δ = 8, 16, 40 against C177's −1.723, −1.917, −1.986). Every virtual pair costs kV − 2μ_bare, so v(μ_eff) = v(−6) − 2(μ_eff + 6)(Δ/V)², with measured slopes −1.53, −1.85, −1.98. At the dilute onset (μ_eff ≈ −2.2), which is where the gas lives, v ≈ 3(Δ/V)², and the resonance function at V/Δ = 8 is −0.49, on the bound side. ΔN = 2 mixing contributes only about 1e-3|t₂| at μ_eff = −12. So the parity-class unbinding at V/Δ = 8 reported from C177's fit does not carry over to the gas. Near the onset the ΔN = 2 vertex is close to resonance and its effect on binding is not computed; the deciding test is a full-chain two-magnon calculation at the onset h_z with the four-particle sector kept.
- Pending. Block DMRG of the model with the vertex at v = 0.1 to 0.3, L = 96 and 192, with the c = 1 and c = 2 controls, and ring runs at higher particle cutoff.

### 9.147 C177: the full qubit chain at the cusp, first readings (Claude Opus 5.5 agent, 2026-11-09; `codex/work_C177/`)

C177 ran the actual nearest-neighbour qubit chain, not the effective model, on C173's bond family at V/Δ = 8 to 24. It added a Davidson DMRG with an exact fermion-parity selector (matches ED to 1e-11). The registered criterion (P1 to P3: L = 256 chains, FES at χ ≥ 256, energy comparison) awaits its Typhon job.

- Open chains on the symmetric stratum J_x + J_y = 0 [numerical]:
  - At V/Δ = 8, S̄ = 1.45, 1.64 and 1.763 at L = 48, 96 and 128, with c_prof 2.2 to 2.4 and truncation about 1e-8.
  - At V/Δ = 12 and L = 128, S̄ = 1.791 with c_prof 2.20.
  - The two-size estimate from L = 48 to 96 is 2.12, against 1.93 for the effective model at the same cuts. The full chain's S̄ is about 0.2 lower at equal L, with the same slope.
  - A one-valley control on the same chain (J_x + J_y = 0.15|t₂|) reads c_prof 0.98 to 1.11.
  - The valley-polarization correlation ⟨PP⟩ (P = XX + YY) decays to about 1e-3 at L = 96, against a plateau of 0.23 in the control.
  - VUMPS at V/Δ = 8 so far reaches χ = 160 and ξ = 420, with pairwise FES slopes 2.09, 1.56, 1.41, 1.38, 1.47.
- Parity gap [numerical]. L·(E_odd − E_even) on open chains at L = 24, 32, 48 stays bounded (+0.20, −0.07, +0.13 |t₂| at V/Δ = 8), so no pair gap is resolved at these sizes.
- Pair binding in the full chain [numerical]. Excitations at the physical μ are heavily dressed by virtual pairs (one-excitation ⟨N⟩ = 1.3 to 1.8 at V/Δ = 8), and number-conserving fits are poor there. Measured directly on rings, the lowest K = π two-excitation level falls below twice the one-excitation minimum by L = 18 to 20 at V/Δ = 8 and 12. The margin is smaller than in the second-order model. C177 withdrew its earlier "unbound at V/Δ = 8", which was measured at μ_eff = −12, in agreement with C186 (§9.146). So the pairing obstruction is present but weaker at V/Δ = 8 and 12, and the c ≈ 2 readings at L ≤ 128 may be a crossover.
- Parity-breaking perturbations [numerical].
  - A generic point off the parity class (J_x + J_y = 0.02, a symmetric XZ coupling 0.25, h_x = 0.25, all in units |t₂|) is gapped: S̄ = 0.371 at L = 128 and 256, VUMPS ξ = 16.8, density locked at 0.2000, and valley polarization ⟨P⟩ = 0.23.
  - At L = 48, an h_x or symmetric XZ coupling of only 0.02|t₂| drops S̄ from 1.45 to about 0.7 to 0.8, while J_x + J_y = 0.02|t₂| alone changes nothing.
  - These terms produce a uniform ΔN = 1 vertex (§9.146), and at these amplitudes it appears to lock the density and gap the chain. So the c ≈ 2 reading is not shown open in B₂ [heuristic: the parity-breaking terms look relevant].
  - For the transverse-DM route, this means the uniform ΔN = 1 part generated by symmetric terms has to be tuned away, as C186 noted. A fine-tuned c = 2 point would still be a counterexample.
- Tuning caveat. D was tuned at μ_eff = −12, so at the physical μ the V/Δ = 8 bond has |w|/|t₂| ≈ 0.93 to 0.95. It stays exactly on the symmetric stratum.
- Open. P1 to P3 from the Typhon job. Whether a pair gap below about 0.2/L |t₂| makes the reading a crossover. Whether any c = 2 window survives very small parity-breaking terms.

### 9.148 C168's second stage: neither uniform-sample point survives (coordinator, Claude Opus 5.5, 2026-11-09; `codex/work_C168/typhon_out/b/`)

C168's recalibrated screen of the uniform sample (§9.133) left two stage-A points, idx 10193 and 450, for the pre-fixed second stage: FES slope above 1.2 at χ ≥ 256, and reseeded DMRG at L = 128 and 192 keeping c_prof and the two-size slope above 1.2.

- FES [numerical]. Two- and four-site cells to χ = 256. For 10193, 450 and the control 10540 the entropy saturates and every pair slope at χ ≥ 128 is between −2.4 and 0.2. The fitted slopes (c_ls3) are 0.15 and −1.56 for 10193 (n = 2, 4), −1.26 for 450 (n = 2), and 0.015 for the control. The VUMPS iterations hit their 600-step cap with residuals 0.2 to 0.8, and a few χ steps returned garbage energies, so the FES numbers are rough. They show no growth.
- DMRG reseeds [numerical, truncation below 1e-11]. Seeds 1 and 2 at L = 128 give identical energies and profiles for all three points. From L = 128 to 192, c_prof falls from 1.94 to 0.83 at 10193, from 1.21 to 1.04 at 450, and from 1.16 to 1.03 at the control. S_mid falls at 10193 (0.955 → 0.797) and 450 (1.004 → 0.891), and the energy per site at L = 192 is lower than at L = 128.
- Verdict, by the pre-fixed rule. Both points fail the second stage. C168's statistical bound stands: zero detections give p < 7.3e-4 (95%) with the near-free efficiency and p < 1.4e-3 (95%) with the interacting one, for the normalized volume of any c > 1 set in C155's coordinates. It says nothing about small-volume sets such as the cusp corner.

### 9.149 Coordinator quick test of the cusp effective model at fitted couplings, and a partial desktop run (Claude Opus 5.5, 2026-11-09; `codex/work_C187q/`)

C176's blockade-exact DMRG on C173's effective model at |w| = 1, φ = π/2. "fit8" uses the couplings C177 fitted at V/Δ = 8 (V₂ = −1.7334, V₃ = −0.4930). §9.146 later showed these were taken at μ_eff = −12 and do not describe the chain at the gas density, so fit8 is a model point with an unbound pair, not a real bond. "sw" uses the bare couplings, where the pair is bound.

- Typhon, N/L = 1/6, L = 48, 96, 144, χ to 512, truncation 1e-8 to 3e-8 [numerical, `typhon_out/summary.txt`]:
  - Two-size c over windows of 3, 7 and 13 cuts:
    - fit8 1.57, 1.46, 1.40 for 48 → 96 and 1.77, 1.74, 1.69 for 96 → 144;
    - sw 1.71, 1.67, 1.70 and 1.77, 1.75, 1.74;
    - c = 1 control 0.98 to 1.00;
    - c = 2 control 1.92 to 2.06.
  - Parity gap at fit8: L·gap = 3.57, 2.82, 2.72 at L = 48, 96, 144.
  - Both points read about 1.75, between the controls and nearer 2, and the readings rise slightly with L. As expected, the bound and unbound points do not separate at these sizes and this density.
- Jonah's desktop, N/L = 1/12, partial (L = 48 and 96 complete; L = 144 still running, so its rows are not used) [numerical, `desktop_run2_summary_partial.txt`]:
  - fit8 reads 2.07, 1.98, 1.78 for 48 → 96, sw reads 1.62, 1.57, 1.49, and the controls read 0.91 to 1.15 and 1.95 to 2.35.
  - L·gap goes 2.74 → 2.34 for fit8 and 2.51 → 2.25 for sw.
  - At this lower density sw sits below fit8, which is the direction pairing predicts. The sw gap at L = 96 (0.023) is already close to the binding energy 0.0195, so L = 144 should show whether it levels off.

### 9.150 C181's backfill run: the partially polarized liquid passes the first stage at L = 128 (coordinator, Claude Opus 5.5, 2026-11-09; `codex/work_C181/typhon_out_bf/`)

C181's registered M1 (§9.145) predicted that at |w| = 1, n = 3/16, φ = 1.565 to 1.54, the effective model at bare couplings (where the opposite-valley pair is bound) is a partially polarized two-component liquid. Off φ = π/2 the valley offset acts as a field on the pairs, and the pair-gap term oscillates (the 1D FFLO mechanism), so pairing should not gap the liquid. The backfill job ran φ = 1.555 (in M1) and 1.53 (first-stage criterion only) at L = 64 and 128 with N ± 1, plus both controls, at χ ≤ 384.

- Controls [numerical]: c_two(64 → 128) = 1.06, 1.01, 1.05 for c = 1 and 2.15, 1.86, 2.08 for c = 2 (windows W = 1, 6, 12). Calibration passes.
- φ = 1.555 [numerical, truncation 7e-8 to 9e-8 at L = 128]:
  - c_two(64 → 128) = 2.18, 2.05, 2.18;
  - c_prof 2.12 at L = 128;
  - valley imbalance p = 0.009 at L = 64 and 0.065 at L = 128;
  - profiles homogeneous to 1 percent;
  - L·Δ_P = −0.22 → 1.50 (the L = 64 value is negative, a shell effect).
  - The first-stage criterion of C173 PROOF 5.4 passes, and M1 passes at this point.
- φ = 1.53 [numerical]:
  - c_two(64 → 128) = 1.66, 1.83, 1.91;
  - c_prof 2.14;
  - p = 0.12 to 0.13;
  - L·Δ_P = 2.24 → 1.53.
  - The first stage passes. The earlier 64 → 96 reading of about 1 (§9.145) came from too few unpaired particles.
- p grows with π/2 − φ, as M1 predicts, and no parity gap is resolved at L = 128.
- Reading. In the effective model of the parity-class cusp, an imbalanced two-valley liquid reads c ≈ 2 at L = 128, level with the c = 2 control, with no sign of the pair gap. The FFLO mechanism makes this c = 2 even though the pairs are bound at φ = π/2, provided the liquid stays homogeneous and off composite plateaus. Not yet tested: the 128 → 256 stage, FES, the FFLO momentum shift (M6), the full qubit chain at this bond, and the C179 two-pair vertex O₄, whose detuning protects openness at V/Δ = 40 but not at V/Δ = 8.

### 9.151 C178: ring spectroscopy of the cusp model at φ = π/2, interim (Claude Opus 5.5 agent, 2026-11-09; `codex/work_C178/`)

C178 wrote a momentum- and N-resolved ring ED of C173's effective model in the blockaded space (matches C173's ED to 1e-13). It labels levels by valley, scans the flux, and relaxes the charge zero mode exactly by flux. What remains is the pseudospin spectrum, labelled by its momentum κ. The "twist-free" series retune |w| so the bare valleys do not twist the pseudospin chain. Predictions P1 to P4 were registered before any new spectrum. Verdict by that rule: inconclusive, leaning toward pairing at the symmetric point.

- P1, pseudospin gap ∝ 1/L: passes, with log-slopes −1.01 to −1.05 and a valley-balanced ground state. It also passes in the pair-liquid calibration (|w| = 2), so it does not discriminate at these sizes.
- P2, R = D_Néel/D_ps above 1 and not decreasing: fails. R falls with N in every series, without slowing:
  - 1.72, 1.42, 1.12 at n = 1/4 (N = 4, 6, 8, twist-free);
  - 1.67, 1.51, 1.35 at n = 1/5;
  - 1.68, 1.57 at n = 1/6;
  - 0.86, 0.64 at n = 1/12;
  - down to 0.02 or below in the |w| = 2 pair calibration.
- P3, two velocities: passes. At n = 1/5 the velocities are v₁ ≈ 0.9 and v₂ ≈ 2.05, with combination levels matching to 1%. v₁ softens with N.
- P4, Casimir: passes. The coefficient is 2.83 against v₁ + v₂ = 2.96 at n = 1/5, and agrees to 3 to 6% at n = 1/6 and 1/4. The controls reproduce their known values to 1 to 3%.
- Parity gap. L·Δ_pb is flat (2.75 → 2.70 at n = 1/6, |w| = 1), which bounds the pair gap only at about 0.07.
- Reading [heuristic]. R falling without deceleration, already below 1 at n = 1/12 and collapsing at |w| = 2, looks like the gapped side of a marginal flow. That fits C175's bound inter-valley pair, a pair liquid with c = 1 whose crossover scale lies beyond these rings at n ≥ 1/5. Three sizes cannot exclude a flow that stops at an easy-plane fixed point.
- The charge sector is very soft, with √(v_J/v_N) ≈ 3.2. In C179's convention that would make the two-pair vertex relevant on the symmetric line. At |w| = 2 the compressibility turns negative, consistent with C175's clustering.
- Off φ = π/2, small L only: at φ = 1.54 the valley offset partly polarizes the pseudospin, and at 1.50 the ground state is already imbalanced. This is the regime of C181's partially polarized route (§9.150), which these φ = π/2 diagnostics do not test.
- Corrections. C175's claim 9 (a growing L·E for the valley flip) is not twist-relaxed, and relaxed it is flat at 3.75 and 3.72. C173's valley_ed.py puts n(k) on the wrong grid at even N.
- Pending, registered before the output. If the flow continues, R < 1.0 at n = 1/4, L = 40, N = 10, and R ≤ 1.50 at n = 1/6, L = 48, N = 8. If it stops, R ≥ 1.1 and R ≥ 1.55.

### 9.152 C183: the full qubit chain near the cusp, interim (Claude Opus 5.5 agent, 2026-11-09; `codex/work_C183/`)

C183 searches the full qubit chain around the cusp at moderate V/Δ. Its many-body tests are on Typhon. The local results:

- Unbinding depends on the field [numerical]. C183 extended C177's two-body fit to V/Δ = 4 to 16 and |w| = 0.5 to 3 (ring ED at L = 14 to 20, C175's pseudospin criterion).
  - At μ_eff = −12 the parity-class pair is unbound for roughly V/Δ < 10/|w|, with Δ_σ = 0.36 to 0.96, inside the two-component window.
  - At |w| = 1 the resonance moves from V/Δ ≈ 10 (μ_eff = −12) to about 8 (−8), about 6 (−5) and below 4 (−4).
  - So near the onset the pair is bound at every V/Δ ≥ 4 tested. This independently confirms §9.146. The fit degrades above μ_eff ≈ −4, so the onset itself is reached only by the trend.
- One resonance, two sides [derived, all 35 map points]. Δ_σ approaches +1 from below on the unbound side and −1 from below on the bound side.
- No trimers at seven unbound far-field points [numerical].
- The gas window [numerical, V/Δ = 8, |w| = 1]. Below a sharp onset the chain is a gapped dressed vacuum, with up to 0.1 virtual particles per site. Above the gas window there is a first-order jump to a packed blockaded phase at n ≈ 0.45. The gas lives in between, at densities 0.05 to 0.3. A registered local pilot at V/Δ = 5 landed in the vacuum (S̄ flat at 0.46) and tested nothing, and the density targeting of the Typhon job was fixed before it started.
- The transverse-DM bond needs no extra tuning [numerical, derived]. The uniform ΔN = 1 element vanishes exactly, because Γ(XZ − ZX) together with the rest of the bond has the antiunitary symmetry inversion × Π_j Z_j × complex conjugation. The vacuum ⟨X⟩ is 2e-16, against −0.0115 when a symmetric XZ coupling of 0.001 is added. So the DM bond sits on a symmetric stratum and needs no h_x tuning.
- Retuning in the full chain. At V/Δ = 16, v = 0.3, degenerate valleys need J_x + J_y = −0.00400 and D = 0.015212, against the second-order values −0.00525 and 0.016174. The bond is J_x = 0.4979983, J_y = −0.5020017, D = 0.0152122, J_z = 4, Γ_xz = −Γ_zx = 0.13693064, with h_z ≈ −7.9463 at the onset, set by a density ladder.
- Pending. A small job runs the two V/Δ = 8 parity-class points (on and off the symmetric line) to L = 256, with fixed-parity gaps and controls. A larger job screens 180 points of the DM-bond family at V/Δ = 12 and 16 and of C181's off-symmetric point at V/Δ = 12, 16 and 24. It confirms up to ten of them at L = 128 and 256, with parity gaps and FES.

### 9.153 C182: cluster bands and the phase diagram of the cusp model, interim (Claude Opus 5.5 agent, 2026-11-09; `codex/work_C182/`)

C182 asked whether bound clusters reopen the two-valley question one level up. It found no c > 1 opportunity on that side. Its systematic DMRG screen is queued.

- Valley symmetry of cluster bands [proved]. At φ = π/2 every N-cluster band satisfies E_N(K) = E_N(Nπ − K), so odd clusters (trimers, pentamers) have two degenerate minima.
- The ab pair band has a single minimum at K = π for every |w| from 0.3 to 3.9 [numerical]. The collaborator's exact inverse pair mass stays positive on (0, 4) [proved algebra, evaluated], falling from 4 to 0.643.
- aa and bb pairs are unbound for |w| ≤ 3.6 [numerical], and bound above that, 0.06 to 0.13 above the ab pair.
- Pre-emption [proved, elementary]. If two cluster species with equal energy per particle bind to each other, their bound state takes over. With the valley offset δ as a field, the pair/atom crossing at δ = E_b is pre-empted by the aab trimer wherever B₃ > 0, which holds at every |w| checked. That opens a window E_b − 2B₃ < δ < E_b + B₃. At φ = π/2 trimers never lead, because B₃ (0.0006, 0.019, 0.068 at |w| = 1, 2, 3) is below E_b/2.
- Equation of state [numerical, box bounds and L = 40]:
  - at |w| = 3 the dilute gas is self-bound;
  - at |w| = 2, clusters of 4 to 6 form droplets of density about 0.35;
  - at |w| = 1 the chemical potential rises monotonically to n = 0.5, so the candidate region is a uniform phase;
  - at |w| = 0.5 a flat chemical potential for n ≈ 0.2 to 0.3 and a dense droplet near 0.47 indicate a coexistence region.
- Uniform liquids at |w| = 2 and 3 (n = 0.3 to 0.4) read c(40 → 80) = 0.85 to 0.96 with balanced valleys and truncation below 1e-10, against controls of 1.10 and 2.17.
- High readings of 2.4 to 3.9 at |w| ≤ 1 and inside the |w| = 2 droplet regime are inhomogeneous, with mid-chain density 0.03 to 0.12 above N/L that drifts with L. The |w| = 2, δ = E_b lead fails C182's registered homogeneity rule. A homogeneity filter was added to the stage-2 selection before the job started. C181's off-symmetric points of §9.150 pass homogeneity (maximum deviation 0.008 at L = 128).
- Open. The screen and its stage-2 runs, the |w| ≲ 1 near-symmetric question (C176, C178, C181), the small-|w| coexistence region, and the full chain.

### 9.154 Desktop run at density 1/12: the parity gap does not resolve pairing at L ≤ 144 (coordinator, Claude Opus 5.5, 2026-11-09; `codex/work_C187q/desktop_run2_summary_L144.txt`)

This is the low-density version of the quick test in §9.149. It uses the effective model at |w| = 1, φ = π/2 and N = L/12, and compares the unbound model point fit8 with the bare couplings sw, where the pair is bound with E_b = 0.0195. It was run on Jonah's desktop with C176's DMRG, with χ up to 384. At L = 144 the truncation is 1.5e-9 for fit8 and 5.7e-9 for sw. The L = 192 runs are unconverged and not used. [numerical]

- Two-size c, windows of 3, 7 and 13 cuts:

  | | 48 → 96 | 96 → 144 |
  | :-- | :-- | :-- |
  | fit8 | 2.07, 1.98, 1.78 | 2.01, 1.96, 1.85 |
  | sw | 1.62, 1.57, 1.49 | 1.72, 1.69, 1.62 |
  | c = 1 control | 1.15, 1.04, 0.91 | 1.11, 1.04, 0.94 |
  | c = 2 control | 2.35, 2.23, 1.95 | 2.27, 2.20, 2.04 |

- L × parity gap E(N+1) + E(N−1) − 2E(N), at L = 48, 96 and 144:
  - fit8: 2.74, 2.34, 2.19;
  - sw: 2.51, 2.25, 2.16.
- At L = 144 both gaps are 0.015, below the two-body binding energy, and both are still closing roughly as 1/L.
- Reading. The unbound point reads about 0.9 times the c = 2 control at both size pairs. The bound point reads lower but rises with L. At this density (n times the pair size ≈ 1.75) the many-body pair gap of the bound model is well below 0.015, and neither the gap nor the entropy resolves it at L ≤ 144. So the parity gap cannot settle pairing on the symmetric line at accessible sizes. Routes where the pair gap is irrelevant (C181's imbalanced liquid) or absent (the transverse-DM bond) are the ones a finite-size test can decide.

### 9.155 Independent full-qubit DMRG, a rational transverse-DM family, and a uniform-Y correction (OpenAI Codex, GPT-6 family, 2026-11-09)

The d=2 route was prioritized explicitly. Fresh full-spin calculations,
independent of the effective-model blockade numerics, now test the literal
microscopic Hamiltonian. This is AI-assisted research without human review.

- [checked numerical] TeNPy MPOs and ground states agree with independent
  Kronecker-product ED at L=8 (matrix error below 3e-14, entropy below
  3e-12), and C183 L=12 also agrees with the separately run repository
  Davidson solver and sparse ED. No number or parity projection is used
  for the transverse-DM family. Pauli fields are uniform on open chains.
- [numerical] C183's V=16 retuned bond at hz=-7.9463 has three-central-cut
  entropy .32228975 at L=32 and .44301798 at L=64. The two-size c estimate
  is 1.04504. The simpler rational bond below gives 1.02180. Controls
  measured by the same method give .94672 (XX, c=1) and 1.94135 (two XX
  legs, c=2). The latter has discarded weight 1.9e-6, enough for this
  coarse distinction, not a high-precision calibration. Candidate final
  discarded weights are below 1e-10. The high central-half profile
  fits fall from 2.25 to 1.15 (C183) and from 1.99 to 1.01 (rational).
  These readings do not support c>1 at those fields; they do not settle
  the infrared limit or exclude a nearby retuned window.
- [numerical] Independent full-chain field scans and matched L=24/48
  controls are recorded in `contrib/d2_dmrg_priority_20261109/`. Density
  drifts and fit-window dependence must be retained. The original V=40,
  r=.15 point quoted in the chat has Sbar=.04685 at L=32 and .01748 at
  L=64, with middle-half density falling .0292 to .0082. It is not a
  demonstrated finite-density critical point.
- [explicit construction; numerical comparison] A simpler one-field
  test family fixes Jx=255/512, Jy=-257/512, Jz=4, Dz=1/64, Dy=9/64 and
  varies hz. Full-ring ED through L=16 reproduces the reference's nearby
  finite-size behavior but shifts its onset by about 5e-4 in field.
  Removing XX+YY entirely creates a much larger valley imbalance and
  is not a faithful dilute simplification. The rational family is not
  a counterexample claim.
- [proved symmetry; exact finite-ring algebra] C183's antiunitary PIK
  forbids uniform X but permits uniform Y. Its X-only check does not
  prove every uniform charge-one source vanishes. In the canonical
  effective Hamiltonian with real pairing p, facilitated amplitude f,
  w=(a+ib)/V and mu=m/V, exact Gaussian-integer perturbation on L=6,7,8
  gives the first vacuum-to-one source i*2*p^3*f*b/V^5. Lower orders
  vanish exactly. This is a finite-ring polynomial identity, not an
  infinite-volume remainder theorem. The leading Pauli hy coefficient
  at C183's V=16 point is only 2.54e-7, so it does not by itself establish
  a relevant perturbation or a gap. Pure open-chain transverse DM also
  differs from the facilitated form by transverse edge fields.
- [independent data audit] C181's phi=1.53 effective-model c>1 signal at
  L=128 survives multiple fitting windows and visible chi convergence;
  it is the cleaner simpler backup (no transverse DM). Incomplete L=192
  data and the old mismatched c=1 control cannot be used. The microscopic
  parity-family tests remain separate from those effective-model data.
- [new numerical microscopic evidence] The simpler parity-family bond
  w=exp(1.53i)/80, Delta=1, V=40, hz=-19.97589125, transverse DM=0
  gives Sbar=1.35323214 at L=32 (odd parity) and 1.60105887 at L=64
  (even parity), selecting the lower of the two independently optimized
  sectors at each size. Thus c(32 to 64)=2.14523. The L=64 profile
  fits over four preselected windows range from 2.05 to 2.40.
  This was encouraging short-size evidence, superseded as a persistence
  test by the negative L=128 follow-up below. Bulk density
  drifts .25494 to .21826, parity switches, and last sweep energy changes
  remain around 5e-7 (odd sector 8e-7), despite discarded weights below
  2e-9. Local solver convergence remains unresolved.
- [negative larger-size follow-up, incompletely optimized] Both parity
  sectors at L=128 return Sbar=1.68631521 (even) and 1.69191583 (odd).
  Nominal lower-energy selection gives c(64→128)=.78647; every parity
  pairing lies below one (.738–.931). Central-half profile fits fall to
  .749 and 1.076. Interior density stays smooth with modest mean drift.
  Thus the microscopic c≈2 signal is not sustained by the larger test.
  Both chi=128 runs hit the sweep cap. The even branch's last energy
  change 5.77e-6 is comparable to the 9.00e-6 parity splitting, and its
  mean entropy changes .00405 in the final sweep. Recovering c=2 would
  require about .1402 additional central entropy. The data distinguish
  neither a physical crossover from optimization bias nor c=1 from a
  small gap. Small discarded weights are not a convergence certificate.
- [explicit simplification] A nearby parity-preserving rational candidate
  is 10 ZZ + (XX-YY)/2 + (XX+YY)/4000 + (XY-YX)/160, with uniform
  field -199759 Z/10000. It removes transverse DM; d=2, nearest-neighbour
  range and one-site translation invariance are exact. Both parity
  sectors pass independent L=10 ED checks. Independent larger-size
  Davidson results and residual limitations are retained in
  `contrib/d2_dmrg_priority_20261109/SIMPLE_PARITY_CANDIDATE.md`.
- [independent intermediate-size replication, incomplete convergence]
  Both parity sectors of that rational bond are complete at L=32 and 64.
  The selected three-cut entropies are 1.35324023 and 1.59917937, giving
  c_eff=2.12889; L=64 profile fits range 2.05–2.37. This reproduces the
  shorter-size signal with an independent Davidson solver. Its registered
  local residual target 1e-8 fails: final maxima are approximately 1.1e-6
  at L=32 and 1.18e-5/2.16e-5 at L=64. The result is not a converged c
  determination and does not undo the nearby unrounded L=128 downturn.
  An exact integer 4×4 bond representation, all sector data and source
  hashes are provided. All local runs in this contribution are complete.
- [numerical follow-up] C183's transverse-DM point at L=128 reads
  c(64 to 128)=1.10032 and central-half c_profile=.98583. A different
  L=64 initialization reproduces its energy within 6.8e-8 and entropy
  within .000195. Its apparent short-chain excess therefore does not
  become a two-mode signal in this extension. The exact XX controls
  give .946721 and 2.047422 at 32 to 64, and .972275 and 2.025433 at
  64 to 128, exposing a downward finite-chi bias of .106 in the recorded
  two-leg DMRG control. No strict thermodynamic bound is inferred.
- [provisional infinite-chain diagnostic] Corrected environments and a
  coherent excluded-gas initialization find a C183 liquid branch whose
  chi=32→64 entropy/correlation-length proxy is 1.04635; an independent
  ARPACK local solve changes it to 1.04667. Global convergence is not
  established, and the chi=64 local residual remains 1.48e-5. Earlier
  failed infinite runs are retained and excluded. This diagnostic does
  not recover a c≈2 signal; details are in `INFINITE_LIQUID_BRANCH.md`.
- [bounded original-point field scan] At hz=-19.9771875, the original
  V=40 transverse-DM chain has rapidly growing central entropy from
  L=32 to 64 but density nearly halves (.1344→.0708), while L=64
  profile fits vary from 4.05 to 14.58. These inhomogeneous, unconverged
  onset data cannot support a CFT interpretation. Two other registered
  L=32 field points are retained in the same raw-results folder.

Reproducible scripts, full profiles, convergence histories and limitations:
[TeNPy runs](contrib/d2_tenpy_20261109/RESULT.md),
[independent Davidson runs](contrib/d2_dmrg_priority_20261109/RESULT.md),
[simplification](contrib/d2_simplification_20261109/README.md), and
[source/data audit](contrib/d2_priority_audit_20261109/README.md).

### 9.156 C181's main job: the |w| = 1 effective model reads c ≈ 2 over a range of offsets and densities, but without the FFLO momentum shift; calibrated controls for the L = 192 and 256 stages (coordinator, Claude Opus 5.5, 2026-11-09; `codex/work_C181/typhon_out/`, `codex/work_C187q/`)

C181's main Typhon job ran all 124 tasks of the registered grid (§9.145), with the |w| = 1 core already read in §9.150. Every L = 128 run ends with truncation below 1e-7. The controls of §9.150 calibrate it. All c values are c_two(64 → 128) for windows W = 1, 6, 12.

- M1, |w| = 1, n = 3/16 [numerical]:
  - φ = 1.54: 1.57, 1.74, 1.81, c_prof(128) = 2.09, p = 0.13, L·Δ_P = 1.23 → 0.00, homogeneous. The first stage passes.
  - φ = 1.50: 1.84, 1.80, 1.93, c_prof = 2.20, p = 0.26, L·Δ_P = 1.42 → 0.85. The first stage passes. M1's registered prediction for this point (one valley, p > 0.9, c_two < 1.3) is falsified, the liquid stays partially polarized.
  - φ = 1.565: 1.80, 1.76, 1.88, but p = 0.001. The criterion fails on the polarization window only.
- M2, |w| = 1, n = 1/8 [numerical]: φ = 1.55 passes (1.63, 1.80, 1.93, c_prof 2.19, p = 0.11, L·Δ_P = −0.03 → 0.16). φ = 1.565 reads 1.74, 1.63, 1.80 with p = 0.002 and fails on p only.
- The reading does not track p. At N = 24, moving one particle between valleys changes p by 1/12, so p ≈ 0 at φ ≥ 1.555 is either a pseudospin gap below a critical offset or finite-size quantization. Either way, the unpolarized points read c ≈ 1.8 at L = 128, as the polarized ones do. At L ≤ 128 the entropy therefore does not separate the FFLO regime from the symmetric-line regime, where the bound pair says the reading should eventually fall to 1 (§9.139).
- M6 fails [numerical]. The pair correlation D(r) peaks within 0.003 to 0.03 of K_p at every passing point (K_p − π ≈ −0.005 to −0.009, C181 PROOF), against the registered windows 0.5 to 1.5 × πpn, which are 0.04 to 0.12 at φ = 1.53 and 1.54 and 0.08 to 0.23 at φ = 1.50. The unwindowed transform over the 62-point correlation has width 0.13 at φ = 1.565, 1.54 and 1.50 alike, so a doublet at K_p ± 0.16 at φ = 1.50 would have been resolved. This counts against the FFLO reading of the c ≈ 2 signal, though not against c ≈ 2 itself. Either the imbalance from the n(k) weights is not the Fermi-sea imbalance that sets q, or the polarized state is not the two-Fermi-sea liquid the mechanism assumes. These data do not say which.
- |w| = 2 [numerical]:
  - At n = 9/32, c_two reads 4.2 to 5.3 at every φ, and the V₂ = V₃ = 0 control reads 2.54, 1.73, 1.41. S_mid jumps from about 1.5 at L = 64 to 2.1 at L = 128, with c_prof(64) between 0.15 and 1.9. These are not central charges. The L = 64 chain is a self-bound droplet with depleted edges (C181 RESULT), so the two sizes are different states and the readings are uninterpretable.
  - M3 fails as registered (c_two(W = 6) = 4.29, and L·Δ_P grows by 42 percent rather than more than 50). M4's one prediction, saturation at φ = 1.48, is falsified (p = 0.19).
  - At n = 3/8, c_two = 0.87, 0.93, 0.95 at all three φ, with p < 0.04 and L·Δ_P = 41.5 → 45. This is a one-mode liquid with no sign of the second mode.
- M5 [numerical]. At φ = π/2, E/N has its minimum at n = 0.156 (N = 10 at L = 64), inside the registered range, so M5 passes. The minimum is shallow, E/N changing by 5e-4 between n = 1/8 and 3/16. At φ = 1.535, E/N rises monotonically from N = 4, so there is no self-bound liquid there.
- Calibration of the next stage [numerical]. The corrected c = 1 control (default t₂) reads 1.04, 1.00, 1.01 from 128 to 192 (desktop) and 1.04, 1.00, 1.01 from 128 to 256 (Typhon), with truncation 4e-13. The c = 2 control reads 2.09, 1.91, 2.00 from 128 to 192, truncation 2.5e-9. Both stages are calibrated under their registered rules (`codex/work_C187q/desktop_run3/README.txt`, `codex/work_C187q/L256_job/README.md`). The candidate runs have not finished. At L = 192 they are at χ = 192 with truncation 1.4e-6. At L = 256 the φ = 1.53, N = 48 run had reached χ = 768 at sweep 23, with truncation 3e-8. The configuration-error run ctl1_L192_N36 (t₂ = 0) reads 0.43 to 0.45 and is not used.

Reading. At |w| = 1 the effective model keeps reading c ≈ 2 at L = 128 for φ from 1.50 to 1.565 and at n = 1/8 and 3/16. The reading is the same where the imbalance vanishes, and the pair correlations show no FFLO shift. Together these suggest L = 128 is not yet asymptotic, rather than confirming the FFLO mechanism. The full-qubit test of §9.155 points the same way: the phase-1.53 parity-family bond reads 2.15 from 32 to 64 but below 1 from 64 to 128 (incompletely optimized). The L = 192 and 256 stages, now calibrated, decide the effective model.

### 9.157 The L = 192 stage of C181's partially polarized point passes its registered rule (coordinator, Claude Opus 5.5, 2026-11-09; `codex/work_C187q/desktop_run3/`)

Jonah's desktop ran the L = 192 stage registered in `codex/work_C187q/desktop_run3/README.txt` (C176's fibdmrg.py, unmodified, χ to 512). The model is C173's effective model at bare couplings, |w| = 1, n = 3/16, φ = 1.53 and 1.555, each at N = 35, 36 and 37. The L = 128 references are C181's backfill runs (§9.150). Windows below are summarize.py's 3, 7, 13 middle cuts, with W = 1, 6, 12 in brackets.

- Calibration [numerical]. The corrected c = 1 control (ctl1fix, default t₂) reads 1.04, 1.00, 1.01 from 128 to 192, and the c = 2 control reads 2.07, 1.99, 1.92 [2.09, 1.91, 2.00]. Both are within the registered tolerances, so the stage is calibrated. The config-error run ctl1 (t₂ = 0) reads 0.43 to 0.45 and is not used.
- Convergence [numerical]. All six candidate runs completed the full schedule, ending at χ = 512 with truncation 7.6e-8 to 9.0e-8, just inside the registered 1e-7.
- φ = 1.53 [numerical]: c_two(128 → 192) = 1.89, 1.81, 1.78 [1.84, 1.76, 1.83], c_prof(192) = 2.03. L·Δ_P = 2.24, 1.53, 0.67 at L = 64, 128, 192.
- φ = 1.555 [numerical]: c_two(128 → 192) = 1.79, 1.77, 1.81 [1.72, 1.79, 1.83], c_prof(192) = 2.02. L·Δ_P = −0.22, 1.50, 0.73.
- Verdict under the registered rule: consistent with c = 2 at both φ. c_two is at least 1.6 in all three windows, the truncation is below 1e-7, and L·Δ_P falls rather than growing (the rule allowed up to 1.5 times the L = 128 value). No point meets the crossover criterion.

Reading. The effective model's c ≈ 2 now holds over a factor of three in length, at 1.8 against 1.9 to 2.1 for the control, and the parity gap closes faster than 1/L, with no sign of a pair gap. The caveats of §9.156 stand: the pair correlations showed no FFLO shift at L = 128, and an unpolarized point read the same c there. This is the effective model at second order, not the qubit chain. The full-qubit test of the nearby phase-1.53 bond reads below 1 from 64 to 128 (§9.155), at a different density and with incompletely optimized states. The L = 256 stage on Typhon is still running.

### 9.158 C176's production job: on the symmetric line the effective model reads c ≈ 1.9 from L = 128 to 256, with a flat L·Δ_P and power-law single-particle correlations (coordinator, Claude Opus 5.5, 2026-11-09; `codex/work_C176/typhon_out/`)

C176's job ran all 72 tasks of its registered grid (C176 RESULT, registered before any production output). The model is C173's effective model at bare couplings, on the symmetric line φ = π/2, |w| = 1, at fixed even N, with χ up to 1024 (768 for the N ± 1 runs). This is the line where the two-body pair is bound (§9.139). c values are c_two for windows W = 1, 6, 12.

- Controls [numerical]. c = 1 (hard rods): 1.04, 0.98, 1.05 for 64 → 128, and 1.03, 0.98, 1.01 for 128 → 256. c = 2 (two XX legs): 2.10, 1.92, 2.03 and 2.05, 1.95, 1.98. Calibrated at every size pair.
- Registered first stage, n = 3/16 [numerical]: truncation 7.7e-10 at L = 128, c_two(64 → 128) = 1.86, 1.82, 1.94, c_prof(128) = 2.14, valley weights 0.999, C_P(L/4) = 4.7e-4, below threshold and falling. It passes, and the same five tests pass at n = 5/32.
- Registered large-size predictions, n = 3/16 [numerical]: c_two(128 → 256) = 1.88, 1.87, 1.90, against 1.82 at W = 6 for 64 → 128, with c_prof(256) = 2.04 and truncation 1.3e-8 at χ = 1024. C_P(L/4) falls by a factor of 4 from 128 to 256, faster than 1/r. The c = 2 prediction holds. The pseudospin-gap prediction (c_two falling, below 1.4) and the valley-polarized prediction both fail. Over every size pair from 64 to 256, c_two stays between 1.82 and 1.94 at W = 6 and 12.
- Other densities [numerical]. n = 5/32: 1.86, 1.83, 1.87 for 128 → 256. n = 1/8: 1.86, 1.81, 1.83. n = 1/6 (C173's sizes): 1.87, 1.84, 1.90 for 96 → 192. At n = 1/4 the reading falls from 2.04 to 1.63, 1.65, 1.69 between the two size pairs. That density is commensurate with period 4 (C181's amendment), and the drop is what locking would give.
- Parity gap [numerical, same χ for the three energies]. L·Δ_P = 2.44, 2.19, 2.38, 2.37, 2.37 at L = 64, 96, 128, 192, 256 for n = 3/16, falling 2.44 → 2.24 at n = 5/32 and 2.70 → 2.18 at n = 1/8. By the registered reading all three series are gapless at these sizes. The χ = 512 and 768 values agree to 0.3 percent.
- Pair against single-particle correlations [numerical]. At L = 256 and n = 3/16, |G_f| falls by a factor of 17.6 between r = 16 and 63, with a steady envelope exponent near −2.0 at L = 128, 192 and 256, and |D| by 2.9 with exponent −0.78. Both are power laws. The registered pairing reading (|G_f| falling by more than 100 while |D| falls by less than 10) is not met at any density.
- Off the symmetric line, n = 3/16 [numerical]: φ = 1.50 and 1.54 give c_two(128 → 256) = 1.81, 1.79, 1.82 and 1.98, 1.94, 1.98. φ = 1.52 falls to 1.26, 1.30, 1.32, with c_prof(256) = 1.71 and lower S_mid than both neighbours. I have no explanation for that one point. Its energy extrapolation residual (2.5e-5) is the same as its neighbours', so a stuck state is possible but not shown.

Reading. This is the strongest effective-model evidence yet, and it is on the symmetric line itself. The bare-coupling model reads c ≈ 1.9 out to L = 256 at three densities, against 1.95 to 2.05 for the control. The parity gap shows no pair gap, and the correlations show no pairing. What it cannot exclude is a pseudospin gap far below the finite-size spacing. L·Δ_P ≈ 2.4 means level spacing near 0.01 at L = 256, and C175's estimate of the many-body gap at n ≈ 0.17 is about 1e-4, which open chains would not see below L ≈ 1e4. In one dimension an attractive two-component gas has a spin gap that is exponentially small in the coupling, so a crossover at very large L is the natural c = 1 alternative. Telling them apart needs the pseudospin gap directly (C175's scattering data at finite density, or infinite-chain methods), not larger open chains. As before, this is the effective model at second order, not the qubit chain.

### 9.159 C184's open-chain test of the Peotta junction: P1 and P2 fail narrowly near the first-order boundary (coordinator, Claude Opus 5.5, 2026-11-09; `codex/work_C184/typhon_out/`)

C184's job ran the parity-class slice at φ = 1 and 1.2, h_z = 0 (C152's fermion units), on open chains at L = 48 and 96, and at L = 192 on a subgrid, with the XX and critical transverse-field Ising controls. c_pair is C184's estimator (mean of the four central cuts, pairs (48, 96) and (96, 192)).

- Controls [numerical]: XX 0.989, 0.993, and TFI 0.515, 0.507. P4 holds for these two. The Ising × XX control did not run in the job, so the c = 3/2 calibration rests on the exact free-fermion values (1.504, 1.501, `ff_controls.out`).
- P2 fails at one point [numerical]. At φ = 1, Δ = 0.6, V = −2.2, c_pair = 1.23 and 1.26, both above the registered 1.2. Truncation is below 1e-11.
- P1 fails [numerical]. The highest (96, 192) readings are 1.31 at Δ = 0.7, V = −2.4 (1.14 at (48, 96)), 1.26 at Δ = 0.5, V = −2.2, and 1.26 at the P2 point. The registered bound for the best-tuned junction was 1.1.
- Context [numerical]. These points sit next to the first-order boundary. At Δ = 0.3 the entanglement collapses from 0.85 to 0.16 between V = −2.1 and −2.2, and similar jumps appear across the grid, which is the cat-state signature of P3. Ring ED at L ≤ 20 reads c ≈ 1 on the band (C184 RESULT). The readings rise between the size pairs, but they stay well below 3/2, and the ring confirmation that a c > 1 candidate needs is missing.

Reading. C184's registered predictions fail narrowly, and the pass criterion for a candidate is not met either, because the readings are near 1.3 and unconfirmed by rings. The likeliest reading is a c = 1 liquid close to a first-order line, where the open-chain entropy overshoots. A short follow-up (L = 384 at the three points, and the Ising × XX control) would settle it.

### 9.160 C191: the pseudospin RG from the scattering data puts the candidate densities deep on the gapless side, and fixes the frame (Claude Opus 5.5 subagent, 2026-11-09; `codex/work_C191/`)

C191 built the low-energy theory of the two-valley gas of C173's effective model (3.2) on the symmetric line from the two-body S-matrix, with each valley treated in the Bethe ansatz and the inter-valley terms at first order. Its prediction was registered before C189 had numbers.

- The inter-valley backscattering is weak. The reflection amplitude at k_F is 0.02 to 0.10 and passes through zero near n ≈ 0.167. The shallow bound state makes the even channel attractive only for q ≲ 0.15, so the zero-range Yang–Gaudin template of §9.139 does not apply at k_F ≈ 0.3.
- Predicted K_s = 1.34, 1.53, 1.67 at n = 1/8, 5/32, 3/16 (bands of about ±0.2), so cos(√8 φ_s) is irrelevant and there is no pseudospin gap. The separatrix sits at n ≈ 0.07 to 0.085.
- Frame. The valleys are weakly coupled only as mutually bosonic fields. A 600-site two-body ring matches the boson-frame levels to 2.7e-6, and the fermion frame misses by half a level spacing. So the inter-valley pair exponent is 1/K_c, the lattice boson exponent is 1/(4K_c) + 1/(4K_s), and the Jordan–Wigner fermion carries a string over the other valley.
- Checks: Gaudin–Yang gap exponent to 0.8%, hard rods K = (1 − ρ)² to 1e-4, free two-leg XX exponents 1, 1/2, 3/2.
- It warns that C189's registered relations are in the fermion frame and that the open end reflects mostly within a valley (|R_intra| ≈ 0.95 to 0.98).

### 9.161 On the symmetric line the effective model has K_s > 1, from two calibrated exponents (coordinator, Claude Opus 5.5, 2026-11-09; `codex/work_C192/`)

A direct test of §9.160 on C176's converged data (L = 128, 192, 256, χ up to 1024). The frame test was stated before the boson data were opened, with a predicted boson exponent of about 0.32. Scripts are `onebody.py`, `pairconf.py` and `sq2.py`, and `CHECKPOINT.md` holds the numbers.

- Boson one-body exponent, from site pairs symmetric about the centre and the open-chain form tan(πr/2L)^(−η). The c = 2 control (two XX legs) reads 0.500 to 0.502 against the exact 1/2. The candidate reads η_b = 0.326, 0.335, 0.351 at n = 3/16, 5/32, 1/8 at L = 256, and moves by at most 0.004 from L = 128. A value below 1 is impossible in the fermion frame, so this confirms C191's frame.
- Inter-valley pair exponent at momentum π, with the open-chain conformal form. The control reads 0.995 to 1.007 against the exact 1. The candidate reads 0.61 to 0.64 at all three densities and sizes, so 1/K_c ≈ 0.62. The density structure factor (controls within 3%) gives K_c ≈ 1.8, so 1/K_c ≈ 0.55.
- K_s from 1/K_s = 4η_b − 1/K_c.

| n | 4η_b | K_s (pair exponent) | K_s (structure factor) |
|---|---|---|---|
| 3/16 | 1.30 | 1.47 | 1.34 |
| 5/32 | 1.34 | 1.39 | 1.27 |
| 1/8 | 1.40 | 1.27 | 1.17 |

- K_s > 1 needs only K_c < 3.3 at n = 3/16 and K_c < 2.5 at n = 1/8, against a measured 1.6 to 1.8. The backscattering has dimension 2K_s ≈ 2.3 to 2.9.
- Independent support. The odd-even energy excess L·(Δ_P − E'') is 1.52, 1.48, 1.45 at L = 128, 192, 256 and extrapolates to no pair-breaking gap above about 3e-4. The valley imbalance responds linearly to the valley offset (C181-type runs), 1.16 to 1.33 times the free value at L = 256, with no threshold. The valley structure factor slope is 1.2 to 1.55 times free.

Reading. On the symmetric line (3.2) is a c = 2 two-component liquid with K_c ≈ 1.6 to 1.8 and K_s ≈ 1.2 to 1.5, and the hidden pseudospin gap of §9.139 is not there. The values sit at the low edge of C191's bands, with the same rise in density. Two limits. First, the two estimates of 1/K_c differ by about 10%, which moves K_s but not its side of 1. Second, this is the effective model. For the qubit chain C179's two-pair vertex O_4 (dimension 2/K_c ≈ 1.1 to 1.25) is still resonant and relevant on the symmetric line itself, so the qubit candidate has to sit off the line, where O_4 and the backscattering both carry momentum. What this section removes is the threshold. With no pseudospin gap, any nonzero valley offset gives the partially polarized two-valley liquid. The converged full-qubit test is C190.

### 9.162 Off the symmetric line the exponents do not move, and the two-pair mismatch is measured (coordinator, Claude Fable 5, 2026-11-09; `codex/work_C192/`)

§9.161 was on the symmetric line. The qubit candidate has to sit off it, so the same two fits were run on the off-line data (C176 at φ = 1.50, 1.52, 1.54 with L = 128 and 256, and C181 at φ = 1.53, 1.555, 1.565 with L = 128). The script is `offline.py`, and `CHECKPOINT.md` items 14 to 20 hold the relations and the full table.

Off the line no symmetry exchanges the valleys, so the charge and pseudospin labels are dropped. Write ⟨(θ_α(r) − θ_α(0))(θ_β(r) − θ_β(0))⟩ = M_αβ ln r for the two valley phases, and N = M⁻¹ for the dual fields. A vertex e^{i m·θ} has two-point exponent (1/2) mᵀ M m. The two boson exponents give M_aa and M_bb, and the inter-valley pair exponent gives M_ab = η_pair − η_a − η_b. Then the backscattering has dimension N_aa + N_bb − 2N_ab (2K_s on the line), and the two-pair vertex O_4 has dimension 2η_pair and momentum q_4 = 2(Q − π), where Q is the measured pair momentum.

| φ | L | majority fraction | η (majority, minority) | η_pair | backscattering dim | O_4 dim | Q − π |
|---|---|---|---|---|---|---|---|
| π/2 | 256 | 0.50 | 0.326, 0.326 | 0.620 | 2.92 | 1.24 | 0.0000 |
| 1.54 | 256 | 0.55 | 0.329, 0.325 | 0.623 | 2.92 | 1.25 | 0.0053 |
| 1.52 | 256 | 0.59 | 0.329, 0.321 | 0.629 | 2.98 | 1.26 | 0.0095 |
| 1.50 | 256 | 0.64 | 0.327, 0.331 | 0.624 | 2.89 | 1.25 | 0.0168 |
| 1.555 | 128 | 0.53 | 0.323, 0.321 | 0.624 | 3.01 | 1.25 | 0.0049 |
| 1.53 | 128 | 0.57 | 0.325, 0.323 | 0.621 | 2.96 | 1.24 | 0.0075 |

All at n = 3/16 and |w| = 1. The L = 128 rows of the C176 runs agree with their L = 256 rows to 0.01 in every exponent. At n = 1/8 (φ = 1.55 and 1.565, L = 128) the backscattering dimension is 2.6, as on the line.

- Checks. The sum of the two fitted valley momenta equals Q to about 1e-4 in every run, so the fitted pair is the inter-valley one. The c = 2 control reads η = 0.500 at both valleys and 1.00 for the inter-leg pair, which gives M_ab = 0 and dimension 2 for both operators. The momentum systematic is about ±0.002 at L = 128.
- Nothing moves. Between φ = π/2 and φ = 1.50, with the valley polarization rising to 0.28, both valley exponents stay at 0.32 to 0.33 and the pair exponent at 0.62 to 0.63. The backscattering is irrelevant throughout (dimension 2.9 to 3.0), and off the line it also carries momentum. The φ = 1.53 point of the collaborator's bond reads the same as the line.
- O_4. Its dimension is 1.23 to 1.26, so it is relevant when resonant, with RG eigenvalue about 0.75. Off the line its mismatch is q_4 = 0.011, 0.019, 0.034 at φ = 1.54, 1.52, 1.50, and about 0.015 at φ = 1.53. Among parity-even vertices of dimension below 2, O_4 is the only one with a small momentum. The others sit at momenta of order 0.5 or near π.
- Whether the mismatch protects the liquid is a commensurate-incommensurate question. The bare O_4 coupling is fourth order in the pairing, of relative size about 2(Δ/V)² ≈ 1.3e-3 at V/Δ = 40, which puts its locking length at a few thousand sites. The mismatch period at φ = 1.53 is about 420 sites. So the mismatch wins by a factor of order 10 in this estimate, and by more at φ = 1.50, where the exponents are the same. This is an estimate of orders of magnitude, not a calculation.
- The 10% disagreement between the structure-factor slope (1.75 to 1.80) and ΣN/2 from the exponents (1.60) is the same off the line as on it. It is still unexplained. A fit systematic looks more likely than a physical mixing term.

Reading. In the effective model the two-valley liquid is the same liquid on and off the line, with no sign that polarization brings any operator closer to marginal. Limits. This is the effective model. The C181 points are L = 128 only. And at L ≤ 256 the O_4 mismatch is below or close to 2π/L, so these sizes cannot show the cutoff at work. The full-qubit test C190 is the missing piece, and at accessible sizes it should track the effective model whether or not O_4 locks at much longer scales.

### 9.163 C189: further K_s estimators and an infinite-chain calculation both read a stable c = 2 liquid in the effective model, after the registered frame failed (Claude Opus 5.5 subagent, checked by the coordinator, Claude Fable 5, 2026-11-09; `codex/work_C189/`)

C189 asked whether the pseudospin backscattering is relevant in C173's effective model on the symmetric line (|w| = 1, n = 3/16, with 5/32 and 1/8 as checks). It used C176's converged open-chain data and a new infinite-chain code. `RESULT.md` holds the registration and its three amendments, and `PROOF.md` the derivations. The coordinator reran `kest_cp.py` and `summarize_inf.py` and compared the other numbers below with `kest_run.out`, `kest_boson.out` and `validate_inf.out`. They reproduce.

- The registered estimators fail. The registration assumed the fermion frame. Its pair-exponent estimator gives K_s = −0.21 ± 0.18 at n = 3/16 and L = 256 (−0.24 and −0.32 at 5/32 and 1/8). Its Friedel estimator gives −0.52 ± 0.21 and is unstable in L. The one-body estimator failed its controls and was excluded. Under the registered rule a value below 1 reads "relevant", but a negative K_s is impossible, so these numbers say that the frame is wrong and say nothing about K_s. The registered structure-factor window (0.05 ≤ q ≤ 0.3) reads K_c = 1.02 ± 0.17 with passing controls, but πS(q)/q is not flat in the candidate. It rises from 1.0 at q = 0.4 to about 1.55 at q = 0.06 to 0.1.
- Boson-frame estimators (amendment A3). These are post hoc. They were added after the registered ones failed and after the frame check of §9.161, so they do not carry the weight of a registered test. At L = 256:

| n | K_c (pair exponent) | η_B | K_s (from η_B) | 2σ range | K_s (C_P tail) |
|---|---|---|---|---|---|
| 3/16 | 1.61 ± 0.05 | 0.318 ± 0.029 | 1.54 | 1.13 to 2.42 | 1.27 ± 0.11 |
| 5/32 | 1.61 ± 0.05 | 0.327 ± 0.029 | 1.45 | 1.08 to 2.22 | 1.28 ± 0.09 |
| 1/8 | 1.65 ± 0.06 | 0.340 ± 0.030 | 1.33 | 1.01 to 1.96 | 1.25 ± 0.09 |

  The controls read 0.504 for the boson exponent of two XX legs (exact 1/2), 1.020 for their inter-leg pair (exact 1), and 0.918 for hard rods (exact 0.889). The last column is K_s = π²A/Z_P², where −A/r² is the tail of the valley-polarization correlation and Z_P = 1.87 at n = 3/16 (1.92 at 1/8) comes from the off-line runs. It uses densities only. The same procedure on free two-valley fermions reads 0.94 to 1.04 at n = 3/16 and 1/8, against the exact 1, and up to 1.21 at one size at 5/32. It is flat in L (1.256, 1.271, 1.266 at L = 128, 192, 256 for n = 3/16), with no drift toward 1.
- Consistency in the boson frame. K_c = 1.6 agrees with the long-wavelength structure factor (1.5 to 1.65). The predicted Jordan–Wigner fermion exponent K_c + 1/(4K_c) + 1/(4K_s) is 1.93 at (K_c, K_s) = (1.6, 1.5), against a measured 1.9 to 2.0. The slowly decaying 2k₀ wall modulation that spoiled the Friedel fit is the inter-valley coherence, with bulk exponent 1/K_s ≈ 0.7.
- Infinite chain. `idmrg.py` is blockade-exact and number-conserving, and matches the finite code to 1e-11 on 8 to 32 sites. The registered quantity is c_fit, the slope of S against (1/6) log ξ over the largest four χ, together with the correlation lengths by sector.

| run | χ | largest ξ (sites) | c_fit | growth of ξ₁, ξ₂, ξ₀ over the last doubling of χ |
|---|---|---|---|---|
| hard rods (c = 1 control) | 16 to 96 | 281 | 1.015 | 2.47, 2.48, 2.50 |
| free two-valley fermions (c = 2 control) | 16 to 96 | 95 | 1.989 | 2.00, 2.01, 2.00 |
| candidate, n = 1/8 | 16 to 256 | 236 | 2.011 | 1.81, 1.78, 1.85 |
| candidate, n = 3/16 | 16 to 256 | 205 | 2.011 | 1.81, 1.79, 1.83 |

  Both controls are inside their registered tolerances (0.15 and 0.25). In the candidate no sector's correlation length saturates, the ratio ξ₁/ξ₂ changes by a factor of at most 1.02 over the last doubling, and the slope between consecutive χ is 1.99 to 2.02 from χ = 64 on. The registered reading at both densities is a c = 2 liquid with no pseudospin gap up to the ξ reached. The n = 3/16 run stood at χ = 96 (c_fit = 1.980, ξ = 89) when `RESULT.md` was written. The same run has since reached χ = 256, and the row above is the refreshed `local/summary_local.txt`. The two-leg XX control cannot be run in the infinite code, since its legs conserve number separately, so the c = 2 calibration rests on the free two-valley run alone. An unregistered side check is the growth of ξ itself. Finite-entanglement scaling gives ξ ∝ χ^κ with κ = 6/(c(√(12/c) + 1)), so a factor 1.83 per doubling at c = 2 and 2.54 at c = 1. The candidate reads 1.8 and the c = 1 control 2.5, but the c = 2 control at χ ≤ 96 reads 2.0, so this is a rough check only.
- Limits.
  - The open-chain parity gap is not accounted for. The compressibility gives v_c/K_c = 0.50 and the valley susceptibility (χ_s = 0.40 ± 0.04) gives v_s/K_s = 1.59. For a valley-preserving wall these predict L·Δ_P = (π/2)(v_c/K_c + v_s/K_s) = 3.3, against the measured 2.38, which is about 30% lower. This is unexplained.
  - The C_P estimator assumes Z_N = 1. A smaller Z_N would raise K_s.
  - The infinite-chain test resolves a gap only up to the ξ reached, about 200 sites. A gap of the size C175 estimated has a correlation length of order 1e4 sites and would not show. What speaks against it is K_s above 1, and the K_s estimators that work are post hoc.
  - The bipartite fluctuations of the valley imbalance, asked for in the spec, were not done.
  - This is the effective model. Nothing here tests the qubit chain.
- The 10% gap of §9.162. C189's pair-exponent K_c = 1.61 is an independent fit and agrees with the exponent-based 1.60. Its long-wavelength structure factor (1.5 to 1.65) is also consistent with 1.6, while the tapered fit of §9.161 gave 1.75 to 1.80. With πS(q)/q running from 1.0 to 1.55 across the fitted range, the slope depends strongly on the window, so the gap looks like a systematic of the structure-factor fit and not a physical mixing term. No single controlled S(q) analysis has closed it yet.

Reading. Two further routes agree with §9.160 and §9.161. K_s is 1.25 to 1.55 at all three densities, so the backscattering has dimension 2K_s ≈ 2.5 to 3.1, and the infinite chain reads c = 2.01 at n = 1/8 and 3/16 out to correlation lengths of about 200 sites. In the effective model the two-valley liquid on the symmetric line is stable. What remains open is the parity gap above and the qubit chain itself, which is C190's test. Larger χ and the remaining controls are running on the cluster, and `desktop_run/` holds a package for the desktop.

### 9.164 C190 interim: the full-qubit solver is validated and the test is registered, production runs are in progress and there is no c reading yet (Claude Opus 5.5 subagent, checked by the coordinator, Claude Fable 5, 2026-11-09; `codex/work_C190/`)

This is an interim entry. No c_two exists yet and no verdict is drawn.

What C190 is. It is the converged test of the nearest-neighbour qubit chain itself at the candidate (`specs/C189_C191_decisive_tests.md`). The bonds are p153 (the collaborator's φ = 1.53 bond, w = e^{1.53i}/80, Δ = 1, V = 40), their rational simplification, and the symmetric-line bond φ = π/2 at V/Δ = 40 and 80. Chains are open, with L = 64, 128, 256, both parity sectors, and χ up to 384, 512, 768.

Registration (RESULT.md, written before any production output). There are convergence targets C1 to C6 (local residual ≤ 1e-8, stationarity between final sweeps, discarded weight ≤ 1e-7, a χ-stage change in S̄₆ of at most 0.01, agreement between initializations, density within ±0.006 of target) and verdict rules V1 to V6.
- The chain "tracks the effective model" if the ground-branch c_two(64 → 128) and c_two(128 → 256) are both ≥ 1.6 at W = 6 and 12, c_prof(256) ≥ 1.6, and the controls are calibrated.
- It "falls toward 1" if c_two ≤ 1.3 at W = 6 and 12 for either size pair.
- The collaborator's downturn reproduces if p153 at their field gives c_two(64 → 128, W = 1) ≤ 1.2, and does not reproduce if that value is ≥ 1.6.

What is established so far.
- Solver. z2dmrg.py is a two-site DMRG with exact Z₂ parity and a preconditioner built from the exact block Hamiltonians at fixed physical indices (PROOF eq. 5.1). In the p153 pilot and in the Typhon sweeps checked so far, every local solve reaches its residual target (1e-8 or 1e-9) in 50 to 90 matrix-vector products. The diagonal preconditioner used in C177, C183 and by a collaborator (GitHub: salpace) stalls at 1e-7 to 1e-6.
- Validation. At L = 10, in both parity sectors, energies match exact diagonalization to 7e-12 and entropies to 6e-10 for every candidate bond. At L = 32 on the rational bond, both sectors come out below the collaborator's energies, by 6e-8 and 9e-8, with the same parity splitting 7.24e-5.
- Pilot at L = 64, the collaborator's field, even sector (validation only, not a registered reading). I checked these against `local/out/`. S̄₁ = 1.5993, 1.6125, 1.6140 at χ = 128, 192, 256, with discarded weight 9e-9, 1.5e-9, 3e-10. The collaborator's χ = 128 value is 1.60106. So the entropy still moves by 0.015 after χ = 128, although the discarded weight there is already 9e-9.
- The pilot also gives L·Δ_GC/|t₂| = +0.25 (the even sector is lower), below the bound of about 1.2 for an unpaired liquid. The sectors differ in n_mid (0.218 against 0.202), which moves S̄, so a density-corrected estimator was added to the registration, still before production output. S̄₆ sits 0.23 below the effective model at equal n_mid, close to C177's offset at V/Δ = 8.
- Field tuning has finished for one key only, the symmetric-line bond at n = 3/16 (μ_eff = −2.1676). At equal μ_eff the full chain is denser than the effective model by about 0.009.

How far the pilot lead goes. With c_two = 6 ΔS̄/ln 2, a deficit of 0.015 in S̄ at one size is 0.13 in c. The collaborator's downturn (c from 64 to 128 of 0.74 to 0.93, against about 2 at the smaller pair) needs a deficit of about 0.13 in S̄ at L = 128. So under-converged χ = 128 states explain the downturn only if the deficit grows nearly tenfold from L = 64 to 128. That is plausible, since the entropy to be captured grows with L, but it is not shown. The L = 128 production runs decide it under rule V3.

A small correction to RESULT.md's R1. Its χ = 256 pilot values (S̄₁ = 1.6142, E = −650.0414852592) are from the first sweep at that χ. The last sweep gives 1.6140 and −650.0414852455. Nothing above changes.

State of the production runs. They run on Typhon as two payloads, the decisive runs and the controls (`codex/work_C192/TYPHON_OVERLAP.md`). Half an hour in, the three p153 runs at the collaborator's field (L = 256 in both sectors, L = 128 even) were at χ = 128, with sweep times of 500 to 800 s and local residuals at target. The agent estimates about 3 h per sweep at χ = 768 for L = 256. So the 64 → 128 readings come first, and the 128 → 256 reading may need the resume-from-checkpoint path.

Limits. No production run has finished. That the collaborator's convergence problem sits in the local eigensolver is the agent's reading. The lower energies at L = 32 support it, and it is untested at L = 128.

### 9.165 C188: the SU(2)_1/ordered boundary of the cyclic stratum has no c > 1 point, and it is a z = 2 ferromagnetic instability, not the registered first-order line (Claude Opus 5.5 subagent, checked by the coordinator, Claude Fable 5, 2026-11-09; `codex/work_C188/`)

Question. In the cyclic stratum, h = XZ + β YX − tβ ZY with β in (0, 1] and t in [0, 1], an SU(2)_1 liquid meets an ordered phase. A continuous transition there with z = 1 would need a CFT with c > 1, so the boundary was scanned as a possible counterexample.

Verdict. There is no c > 1 point on this boundary and no candidate bond. The boundary is the ferromagnetic instability of the SU(2)_1 liquid, of the same type as the XXZ chain at Δ = −1. The velocity falls linearly to zero, the gaps close like 1/L², and the point is not conformal. The SU(2)_1 side has c = 1 all the way to it.

Evidence (rings to L = 24 and 30, open-chain DMRG to L = 192).
- Location. t_c = 0.565, 0.62, 0.66, 0.72, 0.79, 0.82, 0.87, 0.92 at β = 1, 0.8, 0.7, 0.5, 0.4, 0.3, 0.2, 0.1.
- Velocity. At β = 1, v = 1.66, 1.37, 1.08, 0.78, 0.48, 0.18 at t = 0.30 to 0.55, a constant decrement that extrapolates to zero at t = 0.581. The other cuts behave the same way.
- Not first order. The slope jump per site at the ring crossing falls with size, 0.114, 0.048, 0.024 at L = 12, 18, 24 (β = 1), and 0.054, 0.017, 0.0035 at β = 0.4. L² times the gap at the crossing is 14.5, 12.2, 11.5 at β = 1, so z = 2.
- Central charge on the liquid side. c_cas(18, 24) = 0.92, 0.83, 0.43 at β = 1, t = 0.45, 0.50, 0.55, against 1.004 for the Heisenberg control. Open chains give c_pair at (48, 96) and (96, 192) of 0.928 and 0.931 at β = 1, t = 0.5, 0.827 and 0.862 at t = 0.54, 0.862 and 0.884 at β = 0.7, t = 0.62, against 0.916 and 0.934 for the control. No point away from the XZ corner reads above 1.01.
- The ring entropy slope rises to 1.37 within 0.05 of the boundary. The c = 1 XXZ chain near Δ = −1 does the same (1.04, 1.11), so this is the vanishing velocity and not a larger c.
- A lemma, proved. A uniform π rotation about σ^x maps (α, β, γ) to (−α, −β, γ), so spec H(1, β, −γ) = −spec H(1, β, γ), checked to 7e-14.

Against the registration. Q2 (a first-order crossing) is wrong, the boundary is continuous with z = 2. Q1's direction is wrong, t_c increases toward the XZ corner. Q3, Q4 and Q5 hold. Q6's argument stands with its second case realized. Q7 holds away from the corner.

What I checked. The c_pair values above agree with `dmrg_summary.out` row by row, and the ordered side shows up there as it should (c_pair = 0.52, 0.13 and negative values at β = 1, t = 0.56, 0.58, 0.60, with a nonzero order parameter from t = 0.58).

Limits and open points.
- Near the XZ corner the readings exceed 1 at small sizes and fall with L (β = 0.2: c_pair 1.322 → 0.939 at t = 0.8 and 1.144 → 1.057 at t = 0.6, from (48, 96) to (96, 192)). This is the slow A-Majorana crossover of C154. One point is not settled by the data. At β = 0.2, t = 0.85, next to the boundary, c_pair(48, 96) = 1.666, the L = 96 run there is not converged (dE = 5e-6), and there is no L = 192 run. The pattern at t = 0.8 suggests the same crossover, but that is not shown at this point.
- The junction where the dimer line meets the v = 0 line near the XZ corner is not resolved at these sizes. It has v = 0, so it is not a CFT, but its structure is not pinned down.
- A narrow window on the ordered side has an L-dependent ground-state momentum before the ferromagnet takes over. It does not bear on c at the boundary.
- The field-theory step that excludes a continuous z = 1 boundary (PROOF §2.3) is heuristic. The numerics show that none occurs along the cuts tested.
- The Typhon job was still running at write-up (L = 30 at two more points and the L = 24 grid). It refines numbers and is not needed for the verdict.

### 9.166 The L = 256 stage of the partially polarized points passes its registered rule, at c ≈ 1.72 to 1.85 and not yet converged in bond dimension (coordinator, Claude Fable 5, 2026-11-09; `codex/work_C187q/L256_job/`)

The Typhon run registered in `codex/work_C187q/L256_job/README.md` completed its whole schedule. It is C173's effective model at bare couplings, |w| = 1, n = 3/16, φ = 1.53 and 1.555, each at N = 47, 48 and 49, with C176's fibdmrg.py unmodified and χ up to 768. The L = 128 references are C181's backfill runs (§9.150). The windows are the registered W = 1, 6, 12 of C181's summarize.py. All numbers are in `pair256.out`, written by `pair256.py`.

- Calibration [numerical]. The corrected c = 1 control (ctl1fix) reads 1.040, 1.002, 1.007 from 128 to 256, and the c = 2 control reads 2.085, 1.920, 2.001. Both are within the registered tolerances. The config-error run ctl1 (t₂ = 0) reads 0.67 to 0.68 and is not used.
- Convergence [numerical]. All six candidate runs reached sweep 24 at χ = 768, with truncation 3.5e-8 to 4.0e-8, inside the registered 1e-7. They are not converged in χ. In the last sweep the energy still moves by 3e-7 and S̄₆ by 1e-3, and S̄₆ of the two main runs rises by 0.018 from χ = 384 to 512 and by 0.011 from 512 to 768. The L = 128 references end at χ = 384 with a last increment of 0.007, and the c = 2 control at L = 256 with 0.002.
- φ = 1.53 [numerical]: c_two(128 → 256) = 1.743, 1.720, 1.751, and c_prof(256) = 1.96. L·Δ_P = 1.53 at L = 128 and −0.18 at L = 256.
- φ = 1.555 [numerical]: c_two(128 → 256) = 1.773, 1.842, 1.851, and c_prof(256) = 2.00. L·Δ_P = 1.50 at L = 128 and 0.49 at L = 256.
- Verdict under the registered rule: consistent with c = 2 at both φ. c_two is at least 1.6 in all three windows, the truncation is below 1e-7, and L·Δ_P does not grow. Neither point meets the crossover criterion.

What the pass does not show.

- The size of the reading. The candidate reads 1.72 to 1.85 where the c = 2 control reads 1.92 to 2.09 in the same protocol. The remaining χ drift is of the right size to account for the shortfall, since a further rise of 0.015 in S̄(256) raises c_two by 0.13, but that is an extrapolation and not a measurement. At matched χ = 384 on both sizes the candidates read 1.47 to 1.50 and 1.53 to 1.61 against 1.86 to 2.03 for the control. That comparison is biased low, because the longer chain is less converged at a fixed χ.
- The trend with size. From the backfill files alone, c_two(64 → 128) is 1.66, 1.83, 1.91 at φ = 1.53 and 2.18, 2.05, 2.18 at φ = 1.555. With §9.157's 128 → 192 values the implied 192 → 256 readings are 1.61 to 1.66 at φ = 1.53 and 1.85 to 1.92 at φ = 1.555. Over that short lever arm an error of 0.01 in S̄ moves c by 0.2, so these do not establish a drift, and they do not exclude a slow one at φ = 1.53.
- The parity gap. L·Δ_P is 2.24, 1.53, 0.67, −0.18 at L = 64, 128, 192, 256 for φ = 1.53, and −0.22, 1.50, 0.73, 0.49 for φ = 1.555. The L = 256 values are stable to 0.002 from χ = 384 on. This is not the flat L·Δ_P ≈ 2.4 of the symmetric line (§9.158). A second difference in N that changes sign excludes a growing pair gap, which is what the rule tests. One possible origin is a shell effect of the valley imbalance, whose preferred value off the line need not be an integer of the right parity. That has not been checked.

Exponents [numerical, `codex/work_C192/offline_C187b.out`]. The two-valley estimators of §9.162 give the same numbers at L = 256 as at L = 128. The backscattering dimension is 2.93 at φ = 1.53 and 2.92 at φ = 1.555 (2.96 and 3.01 at L = 128), the two-pair vertex dimension is 1.25 at both, the structure-factor slope is 1.60 and 1.61, and the pair momentum mismatch Q − π is 0.0072 and 0.0046 (0.0074 and 0.0048 at L = 128). So the backscattering stays irrelevant by a wide margin and nothing flows between the two sizes.

Reading. Off the line the effective model stays consistent with c = 2 out to L = 256, with size-independent exponents. Beyond §9.157 this adds a factor of 4/3 in length at the same reading of about 1.8, not a sharper number. A continuation of the six runs to χ = 1024 from their checkpoints would settle how much of the shortfall against the control is truncation. This is still the second-order effective model and not the qubit chain, for which C190 is the test (§9.164).

### 9.167 C190 interim: at χ = 256 the qubit chain reads c ≈ 2.0 from L = 64 to 128 at the collaborator's bond and field, and the earlier downturn does not reproduce (Claude Opus 5.5 subagent, checked by the coordinator, Claude Fable 5, 2026-11-09; `codex/work_C190/`)

This is an interim reading and not the registered verdict. No p153c run has reached its registered final χ (384 at L = 64, 512 at L = 128, 768 at L = 256), so V1 to V6 of `codex/work_C190/RESULT.md` are not evaluated. Details are in RESULT.md, R4.

The bond is the φ = 1.53 parity-family bond at V/Δ = 40 with the collaborator's field h_z = −19.97589125 (§9.155). The solver is the parity-resolved z2dmrg.py of §9.164. All numbers below are from the last sweep at χ = 256 with the mixer off, with every local solve at its residual target (1e-8 at L = 128, 1e-9 at L = 64) and discarded weight 3.5e-9 at L = 128.

- States [numerical]. At L = 128 the even sector has E = −1290.1134121766, S̄₁ = 1.8476, S̄₆ = 1.8386, S̄₁₂ = 1.8314 and n_mid = 0.2211. The odd sector has E = −1290.1134075210 and S̄₁ = 1.8403. A product start with 28 particles reaches the same even state as the random start, to 7e-10 in energy and four digits in S̄. The L = 64 partner is the local pilot of §9.164.
- c_two(64 → 128) [numerical]. On the ground branch (even at both sizes) it is 2.022, 1.991, 2.064 for W = 1, 6, 12. The four parity pairings lie between 1.93 and 2.24. The collaborator's values for the same quantity are 0.74 to 0.93.
- Why the two computations differ [numerical]. At L = 128 the collaborator's states (χ = 128) lie 4.4e-5 (even) and 3.1e-5 (odd) above these in energy, which is 7 to 9 times the parity splitting, and their branch order is reversed. At L = 64 the difference is 2.7e-6. Their S̄₁ at L = 128 is 1.686, against 1.817 at χ = 192 (a tuning run at a neighbouring field) and 1.848 at χ = 256 here. At L = 64 the same sequence is 1.599, 1.613, 1.614. A χ = 128 truncation therefore costs about ten times more entropy at L = 128 than at L = 64, which lowers c_two by about 1.3. The downturn appears to be that truncation.
- Parity [numerical]. The even sector is lower at both sizes, by 4.90e-5 at L = 64 and 4.66e-6 at L = 128, so L·Δ/|t₂| falls from 0.25 to 0.05. There is no pair-gap signature in this size pair.
- Against the effective model at equal n_mid [numerical]. S̄₆ of the qubit chain is below C176's by 0.232 at L = 64 and 0.246 at L = 128. The effective model's c_two between the same two densities is 2.11, against 1.99 here.

What I checked. I recomputed S̄ and all four c_two pairings from the json files without the agent's scripts and they agree to the digits quoted. The χ = 256 stage at L = 128 has two sweeps with the mixer off, between which the energy moves by 4.5e-7 and S̄₁ by 0.0012.

Limits.

- S̄(128) is not converged in χ at the registered level. The step from χ = 192 to 256 is 0.031 (taken from a neighbouring-field run, since sweeps with the mixer on record no entropy), against the target 0.01. A further rise would raise c_two, by 0.09 per 0.01.
- The L = 64 partner is the local pilot at χ = 256 and not the registered χ = 384 run.
- The mid-chain density is 0.221, above the 0.196 of the matched effective-model point, as in §9.155.
- One size pair. The registered question is the pair 128 → 256, which does not exist yet. The two L = 256 runs are in their χ = 192 stage, and their ground branch may be odd.

Timing, revised. Measured sweep times put a complex L = 256 run at 7 to 40 h per χ = 768 sweep, so 1.5 to 6 days for its whole schedule, well above the estimate in §9.164. The L = 128 runs need another 7 to 27 h for their χ = 384 and 512 stages. The host allocations of the two payloads end on 2026-11-10 and 2026-11-11, after which the payloads have to be moved by the rule in `codex/work_C192/TYPHON_OVERLAP.md`. The GPU backend (specs/C190g_gpu_port.md) is being built for the L = 256 runs.

### 9.168 A GPU backend for the C190 solver passes its four gates, runs 48 to 60 times faster per bond at χ = 768, and is in production (Claude Opus 5.5 subagent, checked by the coordinator, Claude Fable 5, 2026-11-09; `codex/work_C190/gpu/`)

The port (`z2dmrg_xp.py`, spec `specs/C190g_gpu_port.md`) is the solver of §9.164 with an array-module switch between numpy and cupy. The algorithm, the precision (complex128) and every tolerance are unchanged. Gates, tables and job details are in `codex/work_C190/gpu/RESULT.md`.

Gates [numerical].

- G1, numpy backend against the original code. Bit-identical at L = 10 in every case. At L = 32, sweep by sweep, within 1.8e-11 in energy, which is the run-to-run noise of the original code with threaded BLAS (4.3e-11).
- G2, GPU against exact diagonalization at L = 10. Within 6.5e-12 in energy and 6.0e-10 in entropy at the candidate bonds, with the SVD on the host or on the GPU.
- G3, GPU against the L = 64 pilot at the collaborator's field. Same 13 sweeps, energies within 2.5e-10 (even) and 1.8e-10 (odd), S̄₁ and S̄₆ within 4e-6. Tolerances were 2e-8 and 5e-4.
- G4, benchmark on the complex L = 256 problem (A100, 40 GB). Per bond update with 80 matvecs, GPU against the original code with 6 threads on the same node.

| χ | GPU (s) | CPU (s) | speed-up per bond |
| :-- | :-- | :-- | :-- |
| 128 | 0.57 | 1.48 | 2.6 |
| 256 | 0.96 | 6.31 | 6.6 |
| 384 | 0.73 | 15.5 | 21 |
| 512 | 0.94 | 29.2 | 31 |
| 768 | 1.57 | 74.9 | 48 (60 with the SVD on the GPU) |

The GPU matvec costs about 3.7 ms at every χ, so the gain grows with χ and there is none below χ = 128, where production stays on the host. Device memory at χ = 768, L = 256 peaks at 27 GB.

Registration amendment, written before any GPU production output. The GPU runs are additional runs of the registered configurations, started from scratch. For each key, size and parity sector the run with the lowest final energy among those passing C1 to C4 is used, CPU or GPU, and the C5 comparison applies between them. V1 to V6 are unchanged.

Production [numerical]. Six A100 GPUs on apollo, one run each, in the order p153c at L = 256, then L = 128 and 64, then sym40_n316, the controls and the other tuned keys. With the mixer off the GPU and CPU runs of the same configuration agree to 1e-9 to 2e-8 in energy at the same sweep, and S̄₆ agrees within 2e-5 at L = 128, χ = 256. A χ = 256 sweep at L = 128 takes 129 to 143 s against 1207 to 1354 s on the CPU.

What I checked. The G3 differences and the χ = 768 benchmark row (speed-up 47.8 per bond, 111 per matvec, 27.2 GB) from the json files. The live state at 22:16 EST: the four p153c L = 128 runs are in their χ = 384 and 512 stages with the mixer off and residuals at target, and the two L = 256 runs are at χ = 192.

Limits.

- G3 ran with the SVD on the host. The production path (host below χ = 128, then GPU with the GPU SVD) is checked at L = 10 and L = 32 and by the agreement with the CPU production runs above, not by a rerun of G3.
- No production sweep at χ ≥ 512 at L = 256 has been timed yet. The estimate that the two p153c L = 256 runs finish around midnight EST on 2026-11-09 is an extrapolation, uncertain by a factor of about 1.5. It replaces the CPU timing of §9.167 for the decisive pair.
- The script that applies the amendment's selection across the CPU and GPU folders is not written. `analyze190.py` reads one folder at a time.
- cupy's import hung twice on apollo01, so the jobs are pinned to apollo02.

### 9.169 C190 at its registered final χ: the qubit chain reads c ≈ 2.2 from L = 64 to 128 at the collaborator's bond and field, and the earlier downturn does not reproduce (coordinator, Claude Fable 5, 2026-11-09; `codex/work_C190/gpu/`)

The GPU runs of §9.168 finished the L = 64 and L = 128 chains of the key p153c (w = e^{1.53i}/80, Δ = 1, V = 40, h_z = −19.97589125) at the registered final bond dimensions, χ = 384 and 512. This replaces the χ = 256 interim reading of §9.167 for this size pair. The L = 256 runs are still in progress, so the rules V1 and V2 are not evaluated here.

Convergence. There are eight runs, two initializations (random and product) in each parity sector at each size. All eight pass C1 to C4: the local residuals are 1e-9, the energy changes over the last two sweeps are at most 1e-9, the discarded weights are 3e-11 (L = 64) and 2.6e-10 (L = 128), and S̄₆ moves by 0.0005 and 0.0046 from the previous χ stage, against the limit of 0.01. The two initializations agree to 1e-9 in energy and to four decimals in the entropies (C5). The one CPU run that has reached its final χ (L = 64, even sector) agrees with its GPU partner to 1e-10 in energy and 1e-4 in S̄₁.

| L | sector | E | S̄₁ | S̄₆ | S̄₁₂ | n_mid |
|---|---|---|---|---|---|---|
| 64 | even | −650.0414853399 | 1.6146 | 1.6090 | 1.5933 | 0.2183 |
| 64 | odd | −650.0414362973 | 1.5941 | 1.5885 | 1.5727 | 0.2018 |
| 128 | even | −1290.1134157022 | 1.8682 | 1.8591 | 1.8515 | 0.2211 |
| 128 | odd | −1290.1134107394 | 1.8598 | 1.8510 | 1.8436 | 0.2130 |

Reading. The even sector is the ground branch at both sizes, and c_two(64 → 128) = 6 [S̄_W(128) − S̄_W(64)]/ln 2 on it is 2.196, 2.165 and 2.234 for W = 1, 6 and 12. Over all four parity pairings the W = 1 value lies between 2.12 and 2.37. The residual χ dependence is small on this scale. The L = 128 entropy rose by 0.016 from χ = 256 to 384 and by 0.005 from 384 to 512, and a geometric continuation leaves about 0.002, or 0.02 in c.

Verdict on V3. The registered rule says the collaborator's downturn (0.74 to 0.93 from 64 to 128) reproduces if the converged ground-branch value at W = 1 is at most 1.2, and does not if it is at least 1.6. The value is 2.20, so the downturn does not reproduce. It was a property of the under-converged L = 128 states of §9.155, as §9.167 suggested.

What this does not show.

- The reading sits above 2 by about 0.2, as the collaborator's 32 → 64 reading (2.15) did. Drifts of this size are common at these lengths, and the two-leg XX control at the same sizes, which calibrates this pair under the registration, has not finished. Until it has, the number says "not falling toward 1" and no more.
- The valley weight ratio is 0.78 (even) and 0.76 (odd) at L = 128, a little under the 0.8 that V6 asks for at L = 128 and 256. Both valleys are occupied, with the partial polarization expected off the symmetric line.
- L·Δ_GC/|t₂| is +0.251 at L = 64 and +0.051 at L = 128, inside the bound of V5 with no growth so far. V5 needs L = 256.
- The decisive pair is 128 → 256. At collection the two L = 256 runs were at χ = 384 (even) and 512 (odd) of 768, and no reading is taken from them.

Also started, with Jonah's approval: a continuation of the effective-model L = 256 stage of §9.166 to χ = 1024, with the c = 2 control and the L = 128 partners carried to the same χ, so that the reading is taken at matched bond dimension. The rule is registered in `codex/work_C187q/L256_job/chi1024/README.md`.

Addendum, the same size pair in the effective model (C190 subagent, Claude Opus 5.5, arithmetic checked by the coordinator; `codex/work_C190/RESULT.md` R5). Interpolating C176's effective-model entropies to the full chain's mid-chain densities gives S̄₆ = 1.8405 at L = 64 and 2.0842 at L = 128, so the effective model itself reads c_two(64 → 128) = 2.11 at W = 6, against 2.165 for the qubit chain. The offset S̄₆(qubit) − S̄₆(effective) is −0.2315 at L = 64 and −0.2251 at L = 128 with the final χ = 512 value, a change of +0.006, inside the ±0.05 window of V4 for this pair. The reference carries an interpolation error of about 0.01 in S̄, or 0.1 in c. So the excess over 2 at these sizes is shared by the effective model and reads as a finite-size drift, and over 64 → 128 the qubit chain follows the effective model's size dependence with a constant offset. V4 itself needs L = 256.

### 9.170 The corner point of C188 that read 1.67 crosses over, c_pair(96, 192) = 0.82 once converged (coordinator, Claude Fable 5, 2026-11-09; `codex/work_C188/corner085/`)

Why. C188 (§9.165) left one loose end. At (α, β, γ) = (1, 0.2, −0.17), t = 0.85, which is 0.02 in t from the ferromagnetic boundary near the XZ corner, the open-chain pair (48, 96) read c_pair = 1.67, with the L = 96 state not converged and no L = 192 run. A rule was registered before the new runs (README.md in the folder): convergence items K1 to K5 at L = 96 and 192 with χ up to 320, then "crossover to c = 1" if c_pair(96, 192) ≤ 1.15, "stays elevated" if ≥ 1.40.

Result. Every convergence item passes at both sizes. Energy steps are 2e-9 and 1e-8, the discarded weights 8e-10 and 3e-9, and ⟨X⟩ and ⟨Z⟩ vanish to 1e-5, so there is no order. The entropy step from χ = 240 to 320 is 0.0006 at L = 96 and 0.0099 at L = 192, the second marginally inside its 0.01 bound.

| pair | χ = 160 | χ = 240 | χ = 320 |
|---|---|---|---|
| c_pair(48, 96) | 1.666 | 1.706 | 1.711 |
| c_pair(96, 192) | 0.320 | 0.737 | 0.817 |

By the registered rule this is the crossover to c = 1. The (48, 96) value was right as a number (the converged L = 96 state gives 1.71), and it is a small-size transient. The Heisenberg control reads 0.93 for the same pair, and C188's other points near the boundary read 0.86 to 0.88, so 0.82 is in family.

Limits. The L = 192 entropy still rises with χ. The last step is a fifth of the one before it, and one more step of the same size would give 0.90. Reaching 1.15 would take four of them. The profile fits fall with size (0.57 at L = 96, 0.24 at L = 192), which says the entropy profile is flattening in the bulk. The mid-chain bond alternation keeps falling (0.0034 to 0.0024 from L = 96 to 192), so these sizes do not show a dimerized state, and the data do not separate a c = 1 liquid with a small velocity from a very weakly gapped one. Either way the point is not a c > 1 candidate, and C188's verdict stands without the exception.

### 9.171 C190 at L = 256: the qubit chain reads c ≈ 1.94 from 128 to 256 at the collaborator's bond and 2.07 on the symmetric line, the fail condition is excluded, and the L = 256 entropy is not yet converged in χ (coordinator, Claude Fable 5, 2026-11-10; `codex/work_C190/gpu/`)

This is the decisive size pair of C190, from the GPU runs at the registered final bond dimension χ = 768. It is a partial collection. Finished: the odd sector of p153c at L = 256 (the collaborator's φ = 1.53 bond and field), and all four L = 256 runs of sym40_n316 (the symmetric-line bond at V/Δ = 40, tuned to n = 3/16) with one L = 128 partner. Not finished: the even sector of p153c at L = 256 (running at χ = 768), the remaining sym40_n316 runs at L = 64 and 128, and the controls at L = 256. Tables are in `codex/work_C190/gpu/RESULT.md`, R6.

Readings.

| key | pair | sectors | c_two, W = 1 | W = 6 | W = 12 | c_prof(256) |
|---|---|---|---|---|---|---|
| p153c | 64 → 128 | ground branch | 2.196 | 2.165 | 2.234 | |
| p153c | 128 → 256 | ground branch (even, odd) | 1.902 | 1.936 | 1.971 | 2.028 |
| sym40_n316 | 128 → 256 | even, even | 2.049 | 2.074 | 2.108 | 2.032 |

- Convergence. Every finished L = 256 run passes C1 (residual 1e-9), C2 (energy steps near 1e-9, entropy steps 1e-5 or less) and C3 (discarded weight 3e-10). Two initializations of sym40_n316 agree to 1.5e-9 in energy and 1e-4 in S̄₆ in each sector (C5), and its density is on target (C6).
- C4 fails at L = 256 in every run. S̄₆ rises by 0.024 to 0.027 from χ = 512 to 768, against a bound of 0.01. For the p153c odd sector the sequence at χ = 384, 512, 768 is 2.0203, 2.0559, 2.0828.
- Faithfulness. The offset between the qubit chain and the effective model at the same mid-chain density is −0.2315, −0.2251, −0.2129 at L = 64, 128, 256 for p153c. V4's difference is +0.019, inside its ±0.05. The qubit chain tracks the effective model with a constant offset.
- Parity. L·Δ_GC/|t₂| = +0.251, +0.051, −0.374 for p153c at L = 64, 128, 256 (the last uses the unfinished even run), and −0.129 for sym40_n316 at L = 256. All are far inside V5's 1.5. No pair gap is resolved.
- Valleys. W_min/W_max = 0.999 on the symmetric line and 0.772 at the collaborator's bond, the second slightly under V6's 0.8.

Against the registration (`codex/work_C190/RESULT.md`).
- V2, "falls toward 1", needs c_two ≤ 1.3 at W = 6 and 12 for either pair. It does not hold. The lowest reading is 1.90.
- V1, "tracks the effective model", has its thresholds met (c_two ≥ 1.6 at W = 6 and 12 for both pairs, c_prof(256) ≥ 1.6). It also requires C1 to C4 in every run and the two-leg control within 2 ± 0.25 at the same pairs. C4 fails at L = 256 and the control has not finished, so V1 is not declared.
- As registered, the outcome so far is "unresolved", and the reason is convergence in χ at L = 256, not a low reading.

Direction of the bias. At fixed size the entropy rises with χ, and the L = 128 partner passes C4 (step 0.0046). An unconverged L = 256 state lowers c_two, so 1.94 errs on the low side. Outside the registration, a linear extrapolation in the discarded weight through the last two stages gives S̄₆ = 2.094 at L = 256, 1.861 at L = 128 and c_two = 2.01 at W = 6.

What is left. The even-sector p153c run and the L = 256 controls (in progress on the GPU), the CPU runs of the same configurations as a cross-check, and a further χ stage at L = 256 to close C4. A stage at χ = 1024 needs about 48 GB of GPU memory, more than the 40 GB cards used so far.

Addendum (2026-11-10, 00:50 EST). The even-sector p153c run at L = 256 has finished its χ = 768 stage. E = −2570.2571001943, S̄₁ = 2.0892, S̄₆ = 2.0830, S̄₁₂ = 2.0793, n_mid = 0.2192, c_prof = 2.011. C1 to C3 pass (residual 1e-9, last two sweeps |ΔE| = 1e-9 and |ΔS̄₆| = 1.3e-4, discarded weight 3.7e-10). C4 fails by more than in the odd sector (step 0.087 from χ = 512, where this run sat at 1.996). The two sectors now agree in S̄₆ to 2e-4 (2.0830 and 2.0828), although they came up through χ = 384 and 512 along different paths (1.912, 1.996 against 2.020, 2.056). The ground sector at L = 256 stays odd, L·Δ_GC/|t₂| = −0.374, and the ground-branch readings above are unchanged (1.902, 1.936, 1.971). Over all parity pairings c_two(128 → 256) at W = 1 runs from 1.90 to 1.99. The XX control on the GPU reads 0.939 from 64 to 128 at W = 6, equal to the exact free-fermion value to 1e-8 in entropy. The two-leg control and the XX control at L = 256 are still running, so the verdict stands as written above.

### 9.172 Four CPU jobs collected: the ring anisotropy ratio keeps drifting (C178), a V/Δ = 8 qubit point reads 1.76 from 128 to 256 at χ = 256 (C183), and two backfill screens ended before their L = 192 and L = 96 runs (C185, C186) (coordinator, Claude Fable 5, 2026-11-10; `codex/work_C178/typhon_out/`, `codex/work_C183/typhon_out_f/`, `codex/work_C185/typhon_out_bf/`, `codex/work_C186/typhon_out_bf/`)

C178, ring spectroscopy of the effective model (registered item 1 of its RESULT.md, written before these outputs existed). The question was whether the pseudospin ratio R = D_Néel/D_ps keeps falling with size or stops.
- Twist-free point at n = 1/4, |w| = 1.5307. R = 1.724, 1.424, 1.122, 0.819 at (L, N) = (16, 4), (24, 6), (32, 8), (40, 10). Registered: the flow continues if R < 1.0 at L = 40, it has stopped if R ≥ 1.1. It continues, and R has crossed 1.
- Twist-free point at n = 1/6, |w| = 1.0353. R = 1.681, 1.573, 1.463 at (24, 4), (36, 6), (48, 8). Registered: the flow continues if R ≤ 1.50, it has stopped if R ≥ 1.55. It continues.
- The effective anisotropy label Δ_R (one ratio, not a fit, see C178 claim 11) is 0.50, 0.72, 0.92, 1.13 in the first series and 0.52, 0.64, 0.72 in the second. D_ps keeps its 1/L law (slope −1.01 in both).
- The parity gap stays flat at |w| = 1, n = 1/6. L·Δ_pb = 2.75, 2.70, 2.68 at L = 24, 36, 48.
- Reading. No easy-plane fixed point is visible on rings up to L = 48. At n = 1/4, |w| = 1.53 the ratio has passed the isotropic value. Near the candidate density (n = 1/6, |w| ≈ 1) it is still on the easy-plane side and its drift slows (steps of 0.12 and 0.08 in Δ_R). This is the ring-level form of the open alternative in STATUS item 1, a pseudospin gap that open chains at L ≤ 256 would not resolve. It does not decide it, because a flow that ends at the isotropic point is also gapless.

C183 stage F, the full qubit chain at V/Δ = 8 (job ended at its 12 h limit, so the set is incomplete).
- S1, on the stratum. c₂(64 → 128) = 2.68 and c₂(128 → 256) = 1.76, with χ = 192, 256, 256 and discarded weight 5.5e-10, 3.8e-9, 1.0e-8 (bound 1e-7). c_prof = 1.49, 3.00, 2.05. Bulk density 0.097 at all three sizes. E_odd − E_even = −2.1e-5 and −7.9e-7 at L = 64 and 128, so the parity sectors are degenerate at this resolution.
- Against the registration. F1 (both c₂ at or above 1.6, parity gap not growing) holds at S1. F2 does not. P1 fails on its clause c_prof(256) ≥ c_prof(128) − 0.03, because c_prof drops from 3.00 to 2.05. P4 passes.
- S2, at φ = 1.555. c₂(64 → 128) = 1.40. The L = 256 run did not finish, so F1 and F2 are not evaluated. F1 already fails on the first pair.
- Controls. XX reads 0.86 and 0.98 for the two pairs. The two-leg control reads 1.82 from 64 to 128, and its L = 256 run did not finish, so the 128 → 256 pair has no c = 2 calibration here.
- Reading. A second nearest-neighbour qubit point, at a much weaker blockade than C190's V/Δ = 40, reads close to 2 at the largest pair. It is a first-stage number. It is one initialization at χ = 256, where C190 needed χ = 768 at L = 256 and still saw the entropy rise.

C186 backfill, rule (B) of its RESULT.md. Each task ran under an 11.5 h limit.
- The two decisive v = 0.3 runs at L = 192 did not finish (their L = 96 partners took 8.3 h and 5.2 h). At L = 96, S̄ = 1.138 at μ offset 0.03 (n_mid 0.166, discarded weight 1.6e-8) and 0.953 at offset 0.06 (n_mid 0.199, 5.7e-9). No two-size reading exists, so rule (B) gives no verdict at v = 0.3.
- Estimator controls from 96 to 192 at χ ≤ 128. Hard rods read c = 1.07 (discarded weight 6e-11). Two XX legs read 1.37 with discarded weight 1.6e-5 at L = 192, above the rule's 1e-6 bound. So at χ ≤ 128 the screen would read a converged two-component state near 1.4 at this size pair, between its two thresholds. Rule (B) needs χ of at least 256 at L = 192 before it can separate the cases.
- The v = 0.4 control and the v = 0.1, 0.2 points are in the main job, which is still running.

C185 backfill, R6′ and R7 of its RESULT.md.
- The four blockk runs at L = 96 did not finish, so R6′ has no two-size reading from this job. The same runs are in the main job, still running. At L = 48 (χ = 192, discarded weight below 2e-10) the three-cut mean entropy is 1.153 for the κ = 0, v = 1 reference (n_mid 0.096), and 0.959, 0.996, 0.496 for κ = 1 at φ = 1.625, 1.640, 1.655 (n_mid 0.230, 0.234, 0.188).
- Full qubit chain with the transverse-DM term, V/Δ = 8, L = 64, χ = 256. With λ = 0.3 at φ shift 0.06, S_mid = 0.692 and c_prof = 1.19 (n_mid 0.144). The λ = 0 reference reads S_mid = 1.477 and c_prof = 2.05 (n_mid 0.178). R7 counts support if the DM run exceeds the reference by 0.1. It is lower by 0.78, so there is no support at this φ shift. The other shifts (0, 0.10, 0.14) are in the main job.

What is left. The main jobs of C185 and C186. For rule (B) of C186, a rerun of the L = 192 points and the two-leg control with a longer limit and χ ≥ 256. For C183, the L = 256 runs of S2 and the two-leg control, and a χ sequence at S1.

### 9.173 The χ = 1024 continuation of the partially polarized points: the reading is converged at 1.73 to 1.86, the registered verdict is E1, and the shortfall against the control is of the size of a shell effect in N (coordinator, Claude Fable 5, 2026-11-10; `codex/work_C187q/L256_job/chi1024/`)

The continuation registered in `codex/work_C187q/L256_job/chi1024/README.md` finished all ten runs. The six candidate runs at L = 256 (φ = 1.53 and 1.555, N = 47, 48, 49) and the c = 2 control got three more sweeps at χ = 1024, and the three L = 128 partners went from χ = 384 through 512 and 768 to 1024. All numbers are in `typhon_out/analyze1024.txt`, written by `analyze1024.py`.

- Convergence [numerical]. From χ = 768 to 1024, S̄₆ moves by 0.0025 to 0.0027 in the six candidate runs at L = 256 (bound 0.005), by 0.0001 in the L = 128 partners (bound 0.002) and by 0.0003 in the control at L = 256. Every run counts as converged under the registered rule. At χ = 1024 the discarded weight is 1.2e-8 to 1.3e-8 at L = 256 and 7e-10 at L = 128, and the last sweep moves the energy by 7e-9 and S̄₆ by 7e-6.
- R1, matched χ = 1024 [numerical]. c_two(128 → 256) = 1.748, 1.725, 1.756 at φ = 1.53 and 1.777, 1.847, 1.856 at φ = 1.555, for W = 1, 6, 12. The control reads 2.086, 1.921, 2.002. At matched χ = 512 the candidates read 1.61 to 1.64 and 1.67 to 1.75, and at χ = 768 they read 1.70 to 1.73 and 1.76 to 1.83.
- R2, extrapolated in discarded weight [numerical]. 1.757, 1.734, 1.765 at φ = 1.53 and 1.786, 1.856, 1.865 at φ = 1.555, with an uncertainty of 0.003 from the two-point extrapolation. The control gives 2.088, 1.923, 2.004, so the procedure is calibrated.
- Parity gap [numerical]. L·Δ_P at L = 256 is −0.177 at φ = 1.53 and +0.488 at φ = 1.555, the same as at χ = 384 to 0.003. Neither is above 1.5 times its L = 128 value.
- Other readings at χ = 1024 [numerical]. c_prof goes from 2.16 at L = 128 to 1.98 at L = 256 for φ = 1.53 and from 2.14 to 2.02 for φ = 1.555, against 1.93 to 1.97 for the control. The valley polarization is 0.143 and 0.069 at L = 256 (0.132 and 0.065 at L = 128).

Verdict under the registered rule: E1, the reading stays at c ≈ 2. R1 is at least 1.6 in all three windows at both φ, and R2 is at least 1.75 in two of the three windows at both φ. The margin at φ = 1.53 is small, 1.757 and 1.765 against the 1.75 threshold, with 1.734 at W = 6.

What this settles. §9.166 left open how much of the shortfall against the control was truncation. Almost none of it is. The extrapolation adds 0.009 to the matched reading, and the earlier unmatched number (L = 256 at χ = 768 against L = 128 at χ = 384) was already within 0.005 of it, because the remaining rise at the two sizes cancelled. So with the bond dimension converged the registered reading sits below the control by 0.34, 0.20, 0.25 at φ = 1.53 and by 0.31, 0.07, 0.15 at φ = 1.555 (W = 1, 6, 12).

A shell effect in N of the same size (not registered, post hoc). The registered reading uses one particle number at each size, N = 24 and N = 48. At fixed L the entropy is not a smooth function of N. Let d₂ = S̄_W(N) − [S̄_W(N − 1) + S̄_W(N + 1)]/2, which removes a smooth density dependence to first order.

- At L = 256 and χ = 1024, S̄₆ is 2.2554, 2.2480, 2.2727 for N = 47, 48, 49 at φ = 1.53, so d₂ = −0.016, and 2.2574, 2.2630, 2.2762 at φ = 1.555, so d₂ = −0.004.
- At L = 128, where the N = 23 and 25 runs exist only at χ = 384, S̄₆ is 2.0230, 2.0466, 2.0404 at φ = 1.53, so d₂ = +0.015, and 2.0213, 2.0475, 2.0570 at φ = 1.555, so d₂ = +0.008.
- The sign of d₂ follows the parity gap. The N = 24 run at L = 128 is favoured in energy (L·Δ_P = 1.5) and has more entropy than its neighbours, and the N = 48 run at L = 256 is disfavoured or nearly neutral (L·Δ_P = −0.18 and 0.49) and has less.
- Read against the same N = 24 partner, the three L = 256 runs give c_two(W = 6) = 1.79, 1.73, 1.94 at φ = 1.53 and 1.80, 1.85, 1.96 at φ = 1.555.
- Removing d₂ at both sizes shifts c_two by +0.23, +0.27, +0.27 at φ = 1.53 and by +0.16, +0.11, +0.13 at φ = 1.555, to 1.98, 1.99, 2.03 and 1.94, 1.95, 1.98. The control reads 2.09, 1.92, 2.00.

So the shortfall of the registered reading against the control is matched, in size and in sign, by the dependence of the entropy on the particle number at fixed length. This is the same structure as the sign-changing L·Δ_P of §9.166, most likely a shell effect of the valley occupations at partial polarization. The limits of this statement are that the correction was chosen after seeing the data, that the L = 128 second difference is taken at χ = 384, and that the control has no N ± 1 runs, so its own shell effect is not measured.

Reading. Off the symmetric line the effective model stays consistent with c = 2 out to L = 256 with the bond dimension converged, and the registered rule returns E1. The single-N reading of 1.73 to 1.86 is a converged number and is not raised by more χ. Its distance from the control is of the size of the N-dependence of the entropy, so a reading averaged over neighbouring N, or taken in a grand-canonical setup at fixed chemical potential, is the cleaner estimator for partially polarized points. This is still the second-order effective model. The qubit chain itself is the subject of §9.171.

What is left for these points. The N = 23 and 25 partners at L = 128 continued to χ = 1024, and N ± 1 runs of the control at both sizes, would turn the shell correction from a post hoc observation into a calibrated reading.

### 9.174 The full-qubit test: both controls calibrate both size pairs, and a χ = 1024 stage at L = 256 is registered and queued (coordinator, Claude Fable 5, 2026-11-10; `codex/work_C190/gpu/`)

- [numerical] Controls with the same code and schedules, collected at their registered final χ (`codex/work_C190/gpu/RESULT.md` R9). At W = 6 the XX chain reads c_two = 0.939 (64 → 128) and 0.974 (128 → 256), and the two-leg XX chain 1.917 and 1.896. The registered calibration (1 ± 0.15 and 2 ± 0.25 at W = 6) holds for both size pairs. The candidate's ground branch reads 2.165 and 1.936 at the same window (§9.171).
- [numerical] The two-leg control at L = 256 agrees with the exact free-fermion entropies to 4e-4. Its last energy step (3.2e-8 after six sweeps) is above the 1e-8 of C2, and every other criterion passes.
- [numerical] C4 is the criterion that fails for the candidate at L = 256. The stages at χ = 384 and 512 have two sweeps each and are not converged in sweeps (the first χ = 768 sweep of the even sector moves S̄₆ by 0.06), so the registered step from χ = 512 mixes sweep convergence with the dependence on χ.
- A χ = 1024 stage is registered (R8, written before any χ = 1024 output) and queued on the h200 GPUs. It continues the finished χ = 768 runs of the collaborator's bond and of the symmetric-line bond in both parity sectors, and the two-leg control. C4 is then read between χ = 768 and 1024, with a two-point extrapolation in the discarded weight reported next to it. V1 and V2 are unchanged.

### 9.175 Two more tuned qubit bonds read c ≈ 2 from 128 to 256 at χ = 768, with the same criterion open (coordinator, Claude Fable 5, 2026-11-10; `codex/work_C190/gpu/`)

- [numerical] The collaborator's rational bond tuned to n = 3/16 (h_z = −19.975988986), GPU runs at L = 64, 128, 256 in both parity sectors (`codex/work_C190/gpu/RESULT.md` R10). The odd sector is the ground sector at all three sizes. Ground branch at W = 1, 6, 12: c_two(64 → 128) = 2.119, 2.102, 2.180 and c_two(128 → 256) = 2.004, 2.033, 2.065. The even sector reads 1.878, 1.917, 1.952 from 128 to 256.
- [numerical] The symmetric-line bond tuned to n = 1/8 (h_z = −19.976333694). The ground sector alternates (odd, even, odd, with splittings 4.7e-5, 1.6e-6, 5.1e-6). Ground branch: c_two(64 → 128) = 2.353, 2.309, 2.401 and c_two(128 → 256) = 1.930, 1.978, 2.018. Sector by sector from 128 to 256 the even sector reads 1.840, 1.868, 1.902 and the odd sector 2.141, 2.164, 2.193.
- [numerical] Every run meets the residual, sweep-to-sweep and discarded-weight criteria (C1 to C3), and the L = 256 mid-chain density is within the registered window of its target for both keys (0.1995 against 0.1957, 0.1359 against 0.1326). The step in S̄₆ from the previous χ stage (C4, bound 0.01) passes at L = 64 and 128 and fails at L = 256, with 0.025 for the rational bond and 0.014 to 0.015 at n = 1/8.
- Against the registration. The fail condition (c_two ≤ 1.3) is excluded for both keys. The values meet the pass thresholds, but the pass is not declared while C4 is open at L = 256, as for the two keys of §9.171.
- Four tuned keys now read about 2 from 128 to 256 at χ = 768, at two densities and three bonds, with both controls calibrated at the same sizes (§9.174). These two keys are not part of the queued χ = 1024 stage.

### 9.176 C182's Typhon screen of the effective model: the cluster ladder comes out as registered, and four screen points pass the stage-2 reading, all in families already under study (coordinator, Claude Fable 5, 2026-11-10; `codex/work_C182/`)

- [numerical] Cluster ladder at φ = π/2. The onset species among N ≤ 5 is the pair at |w| = 1 and 1.5, the tetramer at |w| = 2, 2.5 and 3, and the pentamer at 3.5. The hexamer per particle (R = 16, 18) lies below both for |w| ≥ 2.5 (−0.10178, −0.11959, −0.10130 at |w| = 2.5, 3, 3.5) and is level with the tetramer at |w| = 2 (−0.06334 against −0.06345, still moving with R). So the registered prediction (pairs at |w| ≤ 1.5, a ladder that keeps falling at |w| ≥ 2) holds, marginally at |w| = 2. Off the line the trimer is the onset species at one point only (|w| = 1.5, δ = E_b), and for |w| ≥ 2 the lowest cluster among N ≤ 5 is the tetramer or the pentamer up to δ = E_b.
- [numerical] DMRG screen, 210 points in (|w|, φ, n) at L = 48 and 96. 101 are homogeneous (mid-chain density within 0.02 of N/L). Of these, 81 read c_two(48 → 96, W = 6) between 0.8 and 1.2, 7 below 0.8, none between 1.2 and 1.4, and 13 above 1.4. Another 85 points read above 1.4 but are self-bound or phase-separated. The controls read 0.98 and 1.02 (c = 1) and 1.91 and 1.81 (c = 2) at W = 6.
- [numerical] Stage 2 (L = 96, 144, 192, χ to 512, with N ± 1 partners) on the six highest, with calibrated controls (c_two(96 → 192, W = 6) = 0.979 and 1.923). Four pass the registered reading: |w| = 0.5 on the symmetric line at n = 1/12 (1.661, from 1.531 at the screen), and three off-line points, |w| = 1.5 at φ = 1.30, n = 1/4 (2.414), |w| = 1.5 at φ = 1.45, n = 1/6 (1.748) and |w| = 1 at φ = 1.45, n = 1/6 (1.566). Two fail (one on the parity-gap rule, and the point that read 7.0 at the screen reads 1.03).
- Against the registration. The prediction that no screen point passes the converged stage is not borne out, and the prediction that φ = 1.45 and 1.30 read below 1.3 everywhere is wrong.
- Assessment. The |w| = 0.5 point is the crossover the registration anticipated (balanced valleys, three consistent size pairs at 1.63 to 1.71, a slowly falling L·Δ_P). The registration counts it only with an independent pair-gap estimate, which does not exist yet. It poses the same question as the main candidate, at weaker binding. The three off-line points are valley-imbalanced (weight ratio 0.17 to 0.31), their two size pairs disagree by more than 1, and their profile fits read 2.2 to 3.7, so they are not converged readings. They lie in the partially polarized two-valley regime of §9.156, §9.157 and §9.173, farther from the line. No new mechanism appears, and the task's verdict on the cluster-band side stands.
- The job's own stage-2 summary mixes the N + 1 runs into the size series and is not valid. The reading here is from `stage2_read.py`, which selects runs by tag (`codex/work_C182/RESULT.md`, last section).

### 9.177 C188's final outputs: the verdict stands on the full L = 24 grid, and at β = 1 the ring crossing at L = 30 lies below t = 0.55 (coordinator, Claude Fable 5, 2026-11-10; `codex/work_C188/`)

- [numerical] The L = 24 ring grid is complete (96 points, 13 momenta each). On the SU(2)_1 side, c_cas(18, 24) is at most 1.06 for β ≥ 0.3, and at β = 1 it is 0.989, 0.917, 0.832, 0.431 at t = 0.30, 0.45, 0.50, 0.55. Readings above 1.3 at β ≥ 0.3 occur only past the boundary or deep in the ordered region, with velocities below 0.1.
- [numerical] The XZ corner keeps falling with size. At β = 0.2 the reading is 1.11 to 1.19 for t = 0.30 to 0.75, against 1.25 to 3.32 at (12, 18). At β = 0.1 it is 1.38 to 6.99, against 1.97 to 17.1. No point shows an SU(2)_3 signature.
- [numerical] L = 30 at β = 1. At t = 0.45 and 0.50, c_cas(24, 30) = 0.934 and 0.872 with the L = 24 velocity (0.946 and 0.897 with the L = 30 one). At t = 0.55 the lowest levels at magnetization m = 10 and 14 lie below the m = 15 singlet by 0.033 and 0.011, so that point is past the ring crossing at L = 30. The crossing sat at t = 0.566, 0.563, 0.566 for L = 12, 18, 24, so its stability in L holds for L ≤ 24 only. The other momentum sectors were not run at that point.
- [numerical] Two more L = 96 open chains (β = 0.4, t = 0.79 and 0.85) are ordered, with c_pair(48, 96) = −1.6 and −1.5.
- The verdict of C188 is unchanged. No c > 1 point appears on the boundary between SU(2)_1 and the ordered phase of the cyclic stratum, which is a ferromagnetic instability.

### 9.178 The full-qubit test is complete at χ = 768: six tuned bonds read 1.94 to 2.03 from 128 to 256, and the CPU and GPU codes agree (coordinator, Claude Fable 5, 2026-11-10; `codex/work_C190/gpu/`)

- [numerical] The GPU production finished, and all 51 runs are done (`codex/work_C190/gpu/RESULT.md` R11). The last two keys are the collaborator's φ = 1.53 bond tuned to n = 3/16 (h_z = −19.975988793) and the symmetric-line bond at V/Δ = 80 tuned to n = 3/16 (h_z = −39.988023246). Ground branch at W = 1, 6, 12: the first reads c_two(64 → 128) = 2.120, 2.104, 2.182 and c_two(128 → 256) = 1.966, 1.996, 2.029, and the V/Δ = 80 bond reads 2.286, 2.253, 2.344 and 1.915, 1.942, 1.975. Sector by sector from 128 to 256 the values lie between 1.88 and 2.07.
- [numerical] Every run meets C1 to C3 (residual 1e-9, last energy step at most 1.9e-9, discarded weight at most 3.1e-10), and the L = 256 density is on target (0.1955 and 0.1952 against 0.1957). The step in S̄₆ from the previous χ stage (C4, bound 0.01) passes at L = 64 and 128 and fails at L = 256, with 0.025 for the first key and 0.034 to 0.035 at V/Δ = 80.
- [numerical] Two initializations of the first key at L = 128 (product and random) end at the same state in both sectors, with energies equal to 1e-10.
- [numerical] CPU against GPU, rule C5. Three CPU runs of p153c with the same schedule have finished (L = 64 even, L = 128 even and odd). Their energies agree with the GPU runs to 4e-12, 3.3e-10 and 1.8e-10, and S̄₆ agrees to 1e-6 against the allowed 0.01. The two codes land on the same states. The CPU runs at L = 256 are still going.
- Against the registration. The fail condition (c_two ≤ 1.3) is excluded for all six keys. On the ground branch at W = 6 they read 1.936, 1.946, 2.033, 1.978, 1.996 and 1.942 from 128 to 256, at two densities, three bonds on the φ line and two values of V/Δ, with both controls calibrated at the same sizes (§9.174). The reading does not move between V/Δ = 40 and 80, where the corrections to the effective model differ by a factor of four. The pass is still not declared, because C4 is open at L = 256 for every key. The queued χ = 1024 stage addresses that for two of them.

### 9.179 The effective model on the infinite chain reads c = 2.0 out to ξ ≈ 500 to 600 sites at three densities, with the c = 1 control narrowly outside its tolerance (coordinator, Claude Fable 5, 2026-11-10; `codex/work_C189/`)

- [numerical] C189's infinite-chain runs on Typhon have reached χ = 768 (`codex/work_C189/RESULT.md`, last section). The hard-rod control is finished. The three candidate densities and the free two-valley control are still running χ = 1024, so their numbers are interim. This extends §9.163, which stopped at χ = 256 and ξ of 205 to 236 sites.
- [numerical] Candidates on the symmetric line at |w| = 1. The registered c_fit over χ = 256 to 768 is 1.995 at n = 3/16, 1.997 at n = 1/8 and 1.993 at n = 1/6, with the local slope between 1.98 and 2.01 on every step. The correlation lengths at χ = 768 are 511, 591 and 532 sites. Over the last doubling of χ the single-particle, pair and neutral lengths all grow by 1.76 to 1.80, and ξ₁/ξ₂ changes by less than one percent. The registered reading is a c = 2 liquid with no pseudospin gap at all three densities, up to the largest ξ reached.
- [numerical] Calibration. The free two-valley control reads 1.977, inside its tolerance of 0.25. The hard-rod control reads 0.847, which misses its tolerance of 0.15 by 0.003, so by the registration the candidate readings are uncalibrated on the c = 1 side. The control hit its iteration cap at every χ ≥ 64 and has ξ = 2900 sites at χ = 768, five times the candidates'. Its local slope is 0.93 to 1.07 up to ξ = 940 and scatters beyond. Fitted over the ξ range the candidates cover it reads 0.96 (post hoc). The candidates stop on tolerance well below the cap.
- Stated limit. The test resolves a pseudospin correlation length up to about 500 to 600 sites. A gap as small as C175's estimate (ξ_s of order 1e4) is still invisible to it.

### 9.180 C177's Typhon job: the registered criterion fails at V/Δ = 8 to 24, on open chains that are not converged at L = 256, and the V/Δ = 8 infinite chain reads about 1.4 (coordinator, Claude Fable 5, 2026-11-10; `codex/work_C177/`)

- [numerical] C177's job is collected (`codex/work_C177/RESULT.md`, last section, `typhon_out/summary_local.txt`). It ran the full qubit chain at band-tuned bonds with V/Δ = 8, 12, 16 and 24, a single-valley c = 1 control and a generic point off the parity class, with the older solver (open chains to L = 256 and χ ≤ 384, VUMPS to χ ≤ 384).
- [numerical] Registered criterion. P1 (open chains) fails at every bond, so by the registration the full chain at these V/Δ does not support the c = 2 candidate at the sizes reached. The two-size slopes from 128 to 256 are 2.420, 0.595, 1.715 and −0.105 at V/Δ = 8, 12, 16, 24. V/Δ = 8 reads 2.3 to 2.5 on every size pair and fails only on the secondary clauses (the profile fit drops from 2.42 to 2.29, and the discarded weight is 1.0e-8 to 1.1e-8 against the bound 1e-8).
- [numerical] The outcome does not discriminate. The L = 256 states have last energy steps of −5e-7 to −1.1e-6 (about 1e-8 at L = 128), S̄ is not monotonic in L at V/Δ = 16, and the well-converged c = 1 control reads 0.74 against 1 with the same estimator. C190 supersedes this half of the job at V/Δ = 40 and 80, with residual 1e-9 and both parity sectors (§9.178).
- [numerical] VUMPS. The slope of S against ln ξ over the last three χ is 1.17 at V/Δ = 8 (1.33 in the four-site cell), 1.56 at 12, 1.98 at 16 and 2.26 at 24. Most runs hit the iteration cap with gradient errors of 1e-4 to 5e-3. The one reasonably converged stretch is V/Δ = 8 at χ = 64 to 192 (ξ = 170 to 500 sites, errors of order 1e-5), and it reads 1.33 to 1.49 in both cells, below the open chains at the same bond and between the two candidates. A crossover length that shortens as V/Δ decreases would give this pattern, and so would an under-converged VUMPS. The job cannot tell them apart.
- [numerical] The generic point off the parity class is gapped at both sizes (S̄ = 0.3707), as in the local runs.
- Open. A converged rerun at V/Δ = 8 to 24 with the C190 solver, and a symmetric infinite-chain code on the full chain, to settle the V/Δ = 8 reading.

### 9.181 K_s on the qubit chain from the boson exponent: 1.4 to 1.6 at L = 256, with a fall from L = 128 that triggers the registered gap-sign clause by the letter (coordinator, Claude Fable 5, 2026-11-10; `codex/work_C193/`)

- Scope. C193 (registered in `specs/C193_C194_sign_tests.md`) reads the pseudospin Luttinger parameter from the stored correlators of the six converged C190 keys (§9.178), with K_s = 1/(4η_b − 1/K_c) as in §9.161. No new runs.
- [numerical] What the stored data hold. The stored C_P is the connected correlator of the bond operator XX + YY, not the inter-valley pair function that the spec assumed. On the candidates it is dominated by q = 0 with exponent 2.04 to 2.09 and carries no information on K_s. So 1/K_c could not be measured on the qubit chain. It is taken from the effective model at matched density (0.62, band 0.55 to 0.64, §9.161), by an addendum written before any candidate correlator was opened. The boson one-body matrix covers a 48-site window at L = 128 and 256.
- [numerical] Calibration, with fit ranges [8, 28] and [12, 34] frozen on the controls. The two-leg control gives η_b = 0.502 and 0.501 and a C_P exponent of 0.999 and 1.004 at L = 128 and 256 (exact values 1/2 and 1). The XX chain gives η_b = 0.504 to 0.507.
- [numerical] Ground branch.

| key | η_b(128) | K_s(128) | η_b(256) | K_s(256) |
| :-- | :-- | :-- | :-- | :-- |
| p153c | 0.293 | 1.81 | 0.309 | 1.63 |
| sym40_n316 | 0.299 | 1.73 | 0.315 | 1.56 |
| rat_n316 | 0.298 | 1.75 | 0.312 | 1.59 |
| sym40_n18 | 0.319 | 1.53 | 0.333 | 1.41 |
| p153_n316 | 0.298 | 1.75 | 0.313 | 1.58 |
| sym80_n316 | 0.297 | 1.76 | 0.313 | 1.58 |

  Errors on K_s are 0.05 to 0.11, and the band on 1/K_c moves K_s by −0.17 to +0.06. The lowest K_s in either parity sector is 1.36 (sym40_n18 at L = 256). In the same window the qubit chain's η_b lies below the effective model's by 0.033, 0.017 and 0.004 at L = 64, 128 and 256.
- [numerical] Change with size. K_s falls by 0.12 to 0.19 from L = 128 to 256 at every key in both sectors, because η_b at fixed separations rises with L (0.265, 0.298, 0.313 at L = 64, 128, 256 at n = 3/16). The same estimator on the effective model gives 1.47 → 1.46 at n = 3/16 and 1.30 → 1.28 at n = 1/8, so the effective model shows no comparable change. The mid-chain density differs by at most 0.012 between sizes and not monotonically, so it does not account for the change. The cause is not identified. Truncation at L = 256, where criterion C4 fails, is one candidate, and the χ = 1024 stage tests it.
- [numerical] Kosterlitz–Thouless bound. Reading the change as dy_∥/dl = −y_1² gives |y_1| = 0.38 to 0.48 (upper limit at most 0.58) against y_∥ = 0.72 to 1.25 (lower limit at least 0.62). So y_∥ > |y_1| at every key, most tightly at sym40_n18 (0.57 against 0.62). The window samples the same separations at both sizes, so the change is not a flow over a doubling of scale and the bound is indicative only. If K_s kept falling at the same rate per doubling it would reach 1 at L of order 1000 to 3000, which these data cannot exclude.
- Against the registration. Q0 holds. Q1 holds in its η_b clause, and its pair clause cannot be evaluated. Q2 fails on its drift clause alone. Q3 holds by the letter, through its drift clause alone (a fall of more than 0.1). Its other two clauses do not hold, since K_s(256) ≥ 1.36 and the bound keeps y_∥ above |y_1|. So the registered outcome is the gap sign, on a clause whose premise (a change of scale) the stored window does not supply.
- Open. A second stage (C193b, registered in the same spec before any new measurement) measures the inter-valley pair function, the density correlator and the one-body function across the whole chain on the stored C190 states, with fit ranges that scale with L.

### 9.182 The pair function and whole-chain correlators on the stored qubit-chain states: the pair exponent and K_c match the effective model, K_s reads 1.3 to 2.0 at L = 256, and the registered letter is the gap sign through the structure-factor estimate alone (coordinator, Claude Fable 5, 2026-11-10; `codex/work_C193/ckpt_meas/`)

- Scope. C193b (registered in `specs/C193_C194_sign_tests.md` before any measurement) measured, on all 42 final C190 states (six keys and two controls, both parity sectors, L = 64, 128, 256, one GPU, no DMRG rerun), the boson one-body function for pairs symmetric about the centre out to r = L − 1, the connected inter-valley pair function at separations 2 and 3, and the connected density correlator. The loaded states reproduce the stored energy to 2.7e-12, S̄₆ to 2e-14 and the stored one-body window to 1.2e-14, and the contraction code matches exact diagonalization on small chains to 3e-15.
- [numerical] Calibration. On the two-leg control at L = 128 and 256, η_b = 0.504 and 0.502, the pair exponent is 1.014 and 1.001, and K_c from the structure factor is 1.026 and 1.039. Recorded before any candidate output was opened: the structure-factor estimator fails at L = 64 (1.48 for an exact value of 1).
- [numerical] Ground branch at L = 256, fit ranges fixed in r/L. The pair exponent is 0.634 to 0.643 (effective model 0.623 at n = 3/16 and 0.617 at n = 1/8), and K_c from the structure factor is 1.70 to 1.73 (effective model 1.79 and 1.81). So 1/K_c on the qubit chain agrees with the effective model to 2 to 4% by two estimators. η_b is 0.282 to 0.286 at n = 3/16 and 0.297 at n = 1/8, lower than the effective model's 0.326 and 0.352 by 0.04 to 0.055. K_s = 1/(4η_b − 1/K_c) is 1.84 to 2.05 with the pair exponent (errors 0.2 to 0.3) and 1.64 to 1.85 with the structure factor.
- [numerical] Change over a doubling of scale. With the pair estimate K_s changes by −0.32 to −0.62 from L = 64 to 128 and by −0.08 to −0.09 from 128 to 256 (errors 0.3 to 0.7), and the Kosterlitz–Thouless bound gives |y_1| ≤ 0.53 to 0.61 against y_∥ ≥ 1.27 to 1.59. With the structure-factor estimate K_s is negative at L = 64, 3 to 4 at L = 128 (errors 1.2 to 2.3) and 1.6 to 1.9 at L = 256. The same estimator gives the same sequence on the effective model (K_c = 0.41, 1.27, 1.79), so that sequence is the finite-size behaviour of the K_c estimator and not a pseudospin flow.
- Against the registration. B0 holds. B1 fails at one key (sym40_n18, η_b = 0.297 against 0.351, a difference of 0.054 against the tolerance 0.05), and every pair exponent lies in the registered range. B2 holds with the pair estimate at every key and sector and fails with the structure-factor estimate. B3 holds by the letter through the structure-factor estimate alone (the bound allows |y_1| ≥ y_∥ at sym40_n18 and p153_n316, and the acceleration clause holds at every key). The rule fixed before the candidates were opened was B3 if either estimate meets it, so the registered outcome is the gap sign. The acceleration clause uses L = 64, where the estimator had already failed its control. The bound clause uses the change of the structure-factor K_s from L = 128 to 256, which comes from K_c moving from 1.2 to 1.7 while the pair exponent, a second measure of 1/K_c, stays at 0.63 to 0.64. So both clauses are carried by the finite-size error of one K_c estimator. The control has K_c = 1, where that estimator is good to 3% at L = 128, and this did not carry over to K_c ≈ 1.6. The registration should have required a calibration that covers the sizes and the K_c a clause uses, and that flaw is the coordinator's.
- [numerical, heuristic, unregistered] The open ends. On the qubit chain the local η_b falls with r/L (0.322, 0.296, 0.274, 0.219, 0.131 over ranges from (0.05, 0.10) L to (0.30, 0.80) L at L = 256, and nearly the same function of r/L at L = 128, 0.292, 0.263, 0.208, 0.125 on the last four ranges), while both controls are flat at 0.500 ± 0.003 and the effective model rises slightly. A proposed cause: the chain conserves only parity, so an open end carries a pair field of boundary dimension 1/K_c ≈ 0.6 < 1, which is relevant and pins the charge phase. A two-exponent fit of the corresponding form gives η = 0.325 to 0.330 at n = 3/16 and 0.353 at n = 1/8, equal to the effective model's values, with a charge part of 0.16 to 0.17 (0.21 at n = 1/8) against 1/(4K_c) ≈ 0.15, and it returns a charge part of −0.002 to +0.006 on both controls and on the effective model. With this η, K_s at L = 128 and 256 is 1.54 to 1.63 and 1.46 to 1.52 at n = 3/16, and 1.32 and 1.30 at n = 1/8, a change of −0.02 to −0.12 (error about 0.11). The bound gives |y_1| ≤ 0.47 to 0.53 against y_∥ ≥ 0.78 to 0.88 at n = 3/16, and a tie at n = 1/8 (0.477 against 0.478). The fit residual stays at 0.015, so the form captures most but not all of the end effect.
- Reading. Every estimate of K_s that is calibrated at the sizes it uses lies between 1.3 and 2.05 at L = 256, and changes by −0.02 to −0.12 from L = 128 on the ground branch, against −0.12 to −0.19 at fixed separations in §9.181. The spread between estimates is a systematic error from the treatment of the open ends. n = 1/8 is the marginal density, in line with C191's separatrix at n ≈ 0.07 to 0.085 (§9.160).
- Open. The χ = 1024 states can be measured with the same job (`codex/work_C193/ckpt_meas/RESULT.md`, last section). A measurement of η_b without open ends (the infinite chain, C196) would remove the end effect.

### 9.183 The infinite qubit chain: neither iDMRG nor VUMPS converges at the candidate field, so C196's registered test was not run (coordinator, Claude Fable 5, 2026-11-10; `codex/work_C196/`)

- Scope. C196 (registered in `specs/C195_C196_third_size_and_infinite_chain.md`) asked for entanglement scaling on the infinite nearest-neighbour qubit chain at the p153c and sym40_n316 fields. An agent built two infinite-chain codes with exact parity labels on C190's solver (two-site iDMRG and one-site VUMPS), validated them locally, and found that neither converges at p153c at bond dimension 32. No Typhon job was submitted, and none of I0 to I4 is evaluated. I1 was not relaxed.
- [numerical] What works. Untruncated iDMRG reproduces open-chain exact diagonalization for L ≤ 10 to 4e-14 at both candidate bonds and on the Ising chain (rerun by the coordinator). The two-leg control matches to L = 8 and is 0.2 too high at L = 10, which `PROGRESS.md` attributes to zero-weight Schmidt states in that degenerate free model. The remaining numbers in this section are the agent's and were not rerun. On a gapped Ising chain (h = 1.05) iDMRG reaches the exact energy density to 1e-13 and VUMPS reaches a gradient of 3e-11. VUMPS on the XX chain converges at D = 16 to 64 and reads c = 1.015, 1.037, 0.975, 0.985 between consecutive D. Its finite-D state breaks particle-number symmetry, so the two parity sectors carry the same transfer spectrum and the per-sector correlation lengths of I2 and I3 would not be informative with this code.
- [numerical] What fails. iDMRG at χ = 32 has no fixed point at p153c. The energy density drifts from −10.00105 to −10.00095 (C190's bulk value is −10.0011226), and the two cell sites alternate between densities 0.08 and 0.34. The XX chain shows the same behaviour at incommensurate filling. The stated cause is that the particle number is conserved up to virtual pairs, so at density 0.215 each two-site step must add 0.43 particles, while exact parity allows only even increments. VUMPS at D = 32 from a random start stalls at gradients of 3e-4 to 1.6e-2 in every stable variant (plain, C solve preconditioned, damped, annealed in the field from μ_eff = −1.0), with every local eigenproblem converged, and the state stays near the dressed vacuum (density 0.002 to 0.003 against 0.215). The annealed run does not reach a liquid even at μ_eff = −1.0, deep in the dense regime.
- Reading. This is a failure of the method and carries no information on c. The effective model's infinite-chain result (§9.179) held the density with exact particle-number conservation, which the qubit chain does not have.
- Against the registration. Not evaluated.
- Open. Three routes for a second attempt, none tried: a uniform-MPS optimizer or VUMPS in a larger cell started from the bulk of a converged C190 open chain, at a bond dimension where the liquid is well represented, or an approximate particle-number label after a unitary that removes the virtual pairs. All candidate-field tests so far were at D = 32 from random or annealed starts. Without the infinite chain, the end-free check of the boson exponent asked for in §9.182 falls to the L = 512 run, for which a prediction of the open-end form is registered in the same spec before any L = 512 output exists.

### 9.184 L = 512 at the collaborator's field: the runs are started, with the registration, gates and first sweeps recorded (coordinator, Claude Fable 5, 2026-11-10; `codex/work_C195/`)

No result yet. This section records what was fixed before production, so that the L = 512 reading can be judged against it.

- Registration. T1 to T4 of `specs/C195_C196_third_size_and_infinite_chain.md` were copied into `codex/work_C195/RESULT.md` before any code was changed. The pass is c_two(256 → 512) ≥ 1.6 at W = 6 and 12 on the ground branch with T1 holding, the fail is c_two ≤ 1.3. The open-end clauses O0 to O3 of the same spec apply to the final states.
- Code. The C190 GPU solver, with only the two environments and two tensors of the current bond on the device and everything else on the host. Only `dmrg()` changed, and no tolerance moved.
- [numerical] Gates. On a GPU at L = 128 the modified code reproduces the unmodified one bit for bit in all 9 sweeps and in the final measurements (S̄₆ = 1.85434606 in both). The local numpy check at L = 20 agrees to 4e-14 in E and 4e-12 in S. A run stopped and resumed differs by 6.9e-11 in S after the resume.
- [numerical] Benchmark at L = 512 (100 matvecs per bond): 1.42, 2.33 and 5.83 s per bond at χ = 768, 1024 and 1536, with 33, 58 and 128 GB of host memory. χ = 2048 would need about 230 GB, above the 160 GB allowed.
- Schedule, fixed before production. Both parity sectors of p153c from the registered random start. Stages to χ = 512 as at L = 256, then χ = 768 and 1024 each run to a stop rule (3 to 6 sweeps, |ΔE| < 4e-9 and |ΔS̄₆| < 5e-4 between sweeps), then a final χ = 1536 (3 to 8 sweeps, same rule). Running the stages to a stop rule separates sweep convergence from the χ dependence, which the two-sweep stages at L = 256 mixed (§9.178).
- [numerical] First five sweeps, checked by the coordinator from the collected logs. Every local solve met its residual target. The bulk energy per site (E(512) − E(256))/256 agrees with (E(256) − E(128))/128 of the L = 256 run at the same stage to 1e-7 to 1.1e-6, and n_mid agrees with L = 256 to 0.0012 (T4 asks for 0.006 at the end).
- Estimate. 14 to 33 h per run, the two sectors in parallel on two GPUs, with a chained job that resumes from the checkpoint. The final-state checkpoints stay on the cluster in the C190 format, and their paths for the §9.182 measurement are in `codex/work_C195/RESULT.md`.

Open. The runs, then T1 to T4 and O0 to O3. S̄(256) for T2 comes from χ = 768 until the χ = 1024 stage at L = 256 has run.


### 9.185 C175's Typhon job: the pairing test in the effective model does not discriminate (coordinator, Claude Fable 5, 2026-11-10; `codex/work_C175/typhon_out/`)

The job registered in §9.139 ended at its time limit. Eleven of fourteen runs finished. The three that did not are the L = 192 runs of both |w| = 1 pairing points and of the |w| = 1 control (each hit the 1400 minute cap per run and left an empty file), so the |w| = 1 part of the registration has no two-size reading. Numbers below are from `summarize_typhon.py`, with its output in `typhon_out/summary.txt`. The estimator is the registered one, c = 6[S̄(192) − S̄(96)]/ln 2 with S̄ the mean of the three middle block cuts, at χ = 128.

| run | n (mid) | c(96 → 192), 3 cuts | 1, 5, 7 cuts | middle half | final truncation (96, 192) | registered |
|---|---|---|---|---|---|---|
| W2_290 (\|w\| = 2, pairing point) | 0.36 | 1.98 | 2.00, 1.92, 1.85 | 1.40 | 4e-11, 2e-9 | below 1.4 |
| W15_255 (\|w\| = 1.5, pairing point) | 0.22 to 0.23 | −0.19 | | | 3e-6, 4e-6 | below 1.4, not reportable |
| W2M_250 (V₂ = V₃ = 0 control) | 0.14 | 1.19 | 0.74, 1.11, 1.12 | 1.20 | 2e-8, 2e-7 | above 1.6 |
| W1R_220 (V₃ = 0, resonance) | 0.07 to 0.08 | 1.75 | 1.65, 1.91, 2.08 | 2.24 | 7e-10, 3e-7 | none |

Against the registration, by the letter.
- W2_290 reads 1.98 with truncation below 1e-7 at both sizes. The prediction (below 1.4) fails, and the registered failure clause applies: the pairing picture fails at this size and density.
- The V₂ = V₃ = 0 control reads 1.19 against a prediction above 1.6. The control fails as well.
- W15_255 has truncation above 1e-6 at both sizes and is not reportable under the registered rule. Its energy still moved by 5e-4 between the last two sweeps at L = 192.
- The clause comparing with C173's reading at n = 0.17 for the same |w| is not evaluable, since C173's runs are at |w| = 1 (and one at w = 0).

What the numbers support [numerical, with the caveats stated]. Less than the letter suggests, in both directions.
- The W2_290 reading is not a clean c = 2. The entropy profile is not of Calabrese–Cardy form (the stored profile fit gives c_cc = −0.46 and −0.33). It is nearly flat with a modulation of 0.2 to 0.3 on a scale of about 12 to 15 blocks, the middle of the L = 96 chain sits at a minimum of that modulation, and the two-size number moves from 2.00 to 1.40 as the window grows from one cut to the middle half of the chain. So this point contradicts the registered pairing prediction but does not establish a two-component liquid at |w| = 2, where §9.139 found that the gas clusters.
- The control cannot calibrate the test. Its L = 192 profile is asymmetric at the 1e-3 level and its energy moved by 3e-5 between the last two sweeps, so the run is not converged at χ = 128, and a two-component point that reads 1.2 means the estimator at these sizes and this χ does not separate c = 1 from c = 2.
- The single-size points at L = 96 (W1_225 at n = 0.06, W1_224 at n = 0.10, W1M_215 at n = 0.06) are converged but carry no reading without their partners.

Conclusion. The test as registered is inconclusive: one pairing point fails its prediction on a profile that is not conformal, the control fails on an unconverged run, and the |w| = 1 points are missing at L = 192. Nothing here changes the case for or against c = 2 on the qubit chain, which now rests on §9.178, §9.181 and §9.182 and on the gap test in progress (C194), which follows the pseudospin gap directly with a converged solver.

Open. A rerun of the three unfinished L = 192 states and of the control with a longer limit and larger χ would complete the registration. It is low priority next to C194, which asks the same question on the model that matters.

### 9.186 Current ground-sector selection audit (OpenAI Codex, GPT-6 family, 2026-11-10)

[numerical arithmetic audit; no new simulation] Recomputing the completed C190g files finds odd ground branches at both L=128 and 256 for sym40_n316. Their W=6 entropy slope is 1.9463, whereas 2.074 is the even-sector comparison retained from the earlier partial collection. The six ground-branch values now span 1.94–2.03. STATUS and §9.178 are corrected; historical partial snapshots remain as recorded. The L=256 bond-dimension test C4 still fails for every candidate, and the committed L=512 snapshots remain unfinished. The stronger c≈2 candidate assessment is unchanged.

The independent coordinating-team reviewer confirmed the sector-selection correction. Reproduction, source hashes, both final-sweep and canonical-state entropy conventions, and all six results: [audit](contrib/d2_status_audit_20261110/README.md). This is AI-only review.

### 9.187 The infinite qubit chain, second attempt: the iDMRG state breaks particle number at the candidate, so the registered test still cannot be run (coordinator, Claude Fable 5, 2026-11-10; `codex/work_C196/b/`)

C196b retried the infinite-chain test of §9.183 with a commensurate cell (32 sites holding 6 particles at sym40_n316, 16 sites holding 2 at sym40_n18) and a start embedded from the effective model, as `specs/C196b_infinite_chain_second_attempt.md` set out. Every route failed the gates, so no production ran and I0 to I4 are not evaluated.

- G1 passes [numerical]. The q-site-cell iDMRG matches exact diagonalization to 1.3e-11 without truncation (q = 2, 3, 4, chains up to 16 sites, at sym40_n316, p153c and the transverse-field Ising chain).
- G2 passes locally at χ ≤ 96 only [numerical]. The XX chain reads c = 0.99, the two-leg XX control 2.15 (2.02 with the even-sector ξ), and the Ising control saturates at ξ = 5.24. The controls were not run on the production χ sequence.
- G3 fails at both candidate bonds [numerical]. At sym40_n316 the energy per site scatters by 5e-10 to 3e-8 between iterations at χ = 64 and S̄ moves by up to 1.6e-3, so I1 is never reached. The cell does not pin the number of particles, which settles at ⟨N_cell⟩ = 6.12, the open-chain bulk density times 32. A field scan moves ⟨N_cell⟩ continuously and leaves the scatter unchanged, so the allowed retune of the field would not help. sym40_n18 and a 64-site cell behave the same way.
- Cause [numerical, with a heuristic reading]. At finite χ the state breaks the approximate U(1) of the real particles into a pair condensate. The anomalous correlator |⟨b b⟩| stays near 0.08 out to 50 sites (1e-6 on the XX control), and the particle-number variance of a block grows linearly with its size (1.60 at 32 sites, 2.77 at 64) where on XX it grows logarithmically. Only the weak third-order pair term pins the condensate's phase, so the iteration wanders. `codex/work_C196/b/PROOF.md` §5 argues that no exact label prevents this. This is a property of the finite-χ variational state, not evidence about the ground state. The open chains of §9.178 conserve particle number up to virtual pairs and do not show it.
- Route 2 (VUMPS with the same cell) collapses to the dressed vacuum, as in the first attempt. Route 3 was not run, since the embedding is not the problem (the embedded start has 12 particles in two cells and reaches the bulk energy within ten iterations).
- Energy reference, corrected. The value −10.0011861 quoted in the spec for (E(256) − E(128))/128 at sym40_n316 is wrong, and it is presumably the even-sector value, while the ground branch at sym40_n316 is odd at both sizes (the collaborator's audit, §9.186, found the same sector issue in the entropy reading). From the C190 GPU energies (−13.1890724885 and −16.4089104718, plus the field offset −9.976047251 per site) the ground branch gives −10.0012022, which I checked by hand. The χ = 64 infinite-chain energy, −10.0011998, sits 2.4e-6 above it.

What this means for the paper. The infinite-chain reading on the qubit chain (referee item 5 in the ordered list) is not available by these routes, and the case rests on the open chains (§9.178, §9.181, §9.182, and L = 512 in progress) and on the effective model's infinite chain (§9.179). A third attempt would need an ansatz that keeps the real particle number symmetric (a U(1)-projected or number-conserving construction in the effective basis), or a different convergence criterion, which would need a new registration. A three-GPU job of C196b runs χ = 64 to 128 at both bonds and the two-leg control. It changes the verdict only if the scatter falls below 1e-10 at χ = 96 or 128, and its output is in `/scratch/lustre/jkudlerflam/c19x/c196b/out_g3/`.


### 9.188 C194, the pseudospin gap followed in density: a gap is seen at low density and closes near n ≈ 0.05 to 0.10, below the candidate, but the registered verdict is formally G3 through an unconstrained fit (coordinator, Claude Fable 5, 2026-11-10; `codex/work_C194/`)

C194 followed the pair-binding gap of §9.139 from the dilute limit up to the candidate density in the effective model, on open chains (parity gap Δ_pb against L, χ = 128) and on the infinite chain (iDMRG, correlation lengths ξ₁ of the pseudospin sector and ξ₂ of the charge sector), with the registration G0 to G4 of `specs/C193_C194_sign_tests.md`. All numbers are from local runs and C176's stored runs. The Typhon payload and the desktop copy, which take the gaps to χ = 384 and the infinite chain to χ = 512, have not run. I reran `gap_table.py` and it reproduces the agent's table to the last digit.

- G0, partly met [numerical]. On both controls the registered gap fit gives |Δ_∞| ≤ 8e-5, consistent with zero, and DMRG matches the exact free-fermion energies (3e-10 for hard rods, gap to 1.5e-5 for the two-leg XX chain at L = 128). The K_s estimator reads 0.500 and 0.998 against exact 1/2 and 1. The ξ₁/ξ₂ part rests on C189's free two-valley run, and the hard-rod infinite-chain control is in the unrun Typhon payload.
- The gap, open chains [numerical]. Δ_pb(L → ∞) = 0.0129 ± 0.0028 at n = 1/48, 0.0085 ± 0.0029 at 1/32, 0.0049 ± 0.0034 at 1/24, 0.0028 ± 0.0008 at 1/20, 0.0013 ± 0.0013 at 1/16, and consistent with zero at 1/8, 5/32 and 3/16. The dilute end extrapolates to about 0.020, against the binding energy 0.018 to 0.020 derived in §9.139. At the largest L the gap still falls with χ, so χ = 128 slightly overstates it.
- The gap, infinite chain [numerical]. ξ₁/ξ₂ falls with χ at n = 1/32 (0.70 to 0.15) and 1/20 (0.83 to 0.47), so the pseudospin correlation length saturates while the charge one grows, the signature of a pseudospin gap. At 1/16 it is flat and then falls to 0.84 by χ = 384, and at 3/40 and 1/12 it is flat (ξ to 266 at 1/12), with entanglement fits c = 1.08, 1.29, 1.72, 2.01 and 2.01 at 1/32, 1/20, 1/16, 3/40 and 1/12. This brackets the boundary at 1/16 < n_c ≲ 1/12 (not a registered estimator).
- K_s [numerical]. At L ≈ 128 it rises with n, 0.89, 0.95, 0.99, 1.29, 1.41 and 1.50 at n = 1/12, 3/32, 1/10, 1/8, 5/32 and 3/16 (the first three at N = 12, uncertain by about 0.25). From L = 128 to 256 it drifts by −0.02 at 1/8 and −0.04 at 3/16.
- n_c [numerical]. From the gapless side, K_s crosses 1 at 0.05 to 0.09 and 1 + |y₁|/2 at 0.08 to 0.10, depending on the fit form. From the gapped side, the registered Kosterlitz–Thouless fit of Δ_∞(n) gives 0.138 with an interval of 0.066 to 0.20, which covers the whole grid and constrains nothing. The 1/ξ₁ fit is not possible, since only n = 1/32 saturates.

Against the registration, by the letter. G1 is met only across densities and provisionally (Δ_∞ above three errors at 1/48 and 1/20, and the infinite-chain condition met only in an unconverged snapshot at 1/32), so G2 is not met. G3 is triggered formally, because the central value of the gapped-side n_c, 0.138, lies above 0.11. Every substantive G3 condition fails: no density at or above 0.11 shows a saturating ξ₁ or a nonzero gap, and K_s at 1/8 and 3/16 stays well above 1. G4, C191's separatrix, fails formally on the same unconstrained estimate, while the gapless-side estimates and the infinite-chain bracket lie inside its 0.05 to 0.11. C191's predicted ξ_s is too long by a factor of two at 1/32 and by more than ten at 1/16.

Reading. The registered outcome is G3 by the letter, and in substance the test is undecided between G2 and G3, as the spec's last clause allows. What it adds is a positive control the earlier tests lacked: the estimators do see the pair gap where theory says it is open, at n ≲ 1/20, and see it close as the density rises, with the boundary near 0.05 to 0.10 on every estimator that constrains it. The candidate densities 3/16 and 0.19 sit 0.09 to 0.14 above that, and n = 1/8 about 0.03 to 0.07 above it, consistent with n = 1/8 reading marginal in §9.182. The χ = 384 gaps at L = 240 and 256 and the infinite chain to χ = 512 (queued on Typhon, with a copy packaged for the desktop) are what would turn the gapped-side estimate from unconstrained into a number.

### 9.189 Initial comparable-coupling search: null preferred-ratio tests and a wider-ratio finite-size lead (OpenAI Codex, GPT-6 family, 2026-11-10)

[numerical; registered local screen, not a counterexample] The new search prefers a largest-to-smallest nonzero Pauli coefficient ratio R at most 3, labels 3<R<=10 exploratory, and keeps wider bridges separate. It records exchange-only ratios, term counts, field resolution, cancellations and invariant exchange scales. An overall rescaling cannot improve these ratios. There are 139 new Hamiltonians, two references and four controls on periodic rings of 10 and 12 sites; selected momentum-resolved checks reach 16 sites, and controls reach 20. The r=40 reference is missed by tiny rings, and an odd/even-leg two-XX control fails their size-pair calibration. Both limitations are retained, not used to exclude a narrow critical region.

[numerical] Full microscopic, parity-resolved C190 DMRG completed 51 MPS runs, including four control states and an independent start. All final sweep gates pass. Three preferred-ratio bridge seeds (R=7/4, 2 and 8/3) and a rational attractive R=2 point do not give a stable entropy excess on the tested sizes. One early flag has a transverse ordered-cat signature. The R=8/3 seed has smoothed entropy slopes 1.002 and 1.004 through 96 sites, with strong density oscillations and much less stable narrow-cut estimates. A lower-density R=14/3 anchor has saturated entropy and parity gap about 0.18. The exploratory R=28/3 bridge reads about 0.90 from 24 to 48 sites. The existing attractive Peotta reference (R=7.77) reads 1.13 from 48 to 96 and remains unresolved; a c slightly above one is not excluded.

[numerical follow-up lead] The simple wider bridge

H = sum_j [(XX-YY)/2 + (XY-YX)/16 + ZZ]_(j,j+1) - (7/4) sum_j Z_j

has R=28, exchange spread 16, and four grouped terms. Its two registered entropy windows give slopes 2.780/3.024 from 24 to 48 sites and 2.707/2.683 from 48 to 96. At 48, chi96->192 changes the ground-sector entropy by 0.00181, and a second start agrees to 1.5e-7 in entropy. At 96, chi192 is converged in sweeps but has no converged larger-chi comparison; discarded weight is 4.6e-8 and energy variance per site is 5.0e-7. Density remains strongly oscillatory, its bulk mean drifts 0.334->0.301->0.296, profile fits depend on the window, and no conformal dispersion or infinite-chain result is available. Thus it is a finite-size anomaly, not an identified CFT. Its ratio is about 114 times smaller than the old reference's but still outside both preference bands. A nearer R=14 bridge has an even larger, inconsistent short-chain slope and is not promoted to a central-charge measurement.

The physical microscopic sites and bonds remain qubits and nearest neighbour throughout; the two-XX range-two chain is a control only. Checkpoints are local and preserved; no cluster job was submitted and no run remains active. The next useful work is a larger-chi/size test of the wider bridge, fine field/density targeting toward the exploratory points, and rational real-hopping offsets; the proposed further grid has not been run. Code, protocols, raw records, source audits, plots and all null results: [contribution](contrib/d2_order_one_20261110/README.md). Independent AI checks, no human review.

### 9.190 A simpler exploratory candidate with entropy estimates near 1.45 through 192 sites (OpenAI Codex, GPT-6 family, 2026-11-10)

[numerical; promising finite-size candidate, not an identified CFT] The continued comparable-coupling search found the fixed physical-qubit nearest-neighbour Hamiltonian

H = sum_j [(XX-YY)/2 + (XY-YX)/10 + (5/8) ZZ]_(j,j+1) - (7/8) sum_j Z_j.

Here X,Y,Z are Pauli matrices. Multiplying H by 40 gives integer coefficients (20,-20,4,-4,25,-35). There are four grouped terms and six Pauli strings. The nonzero coefficient spread is R=8.75, within the exploratory 3<R<=10 band but outside the preferred R<=3 band. Exchange singular values are (0.625,0.600,0.400), spread 1.5625. The exact Jx+Jy=0 relation and all tuning/cancellation scales are disclosed.

[numerical] Full microscopic DMRG completed both parities at L=24,48,64,96,192. The lower variational branch is even throughout. Central-three/fractional entropy size estimates are 1.6826/1.7038 at 24→48, 1.4185/1.4363 at 48→96, and 1.4475/1.4530 at 96→192. The intermediate L64 check gives a smooth trend near 1.43. Final L192 profile fits are 1.4539/1.4630 with RMS residuals 0.00201/0.00181. All selected sweep gates pass. At L192, chi192→256 changes the two entropies by only 0.000261/0.000247; an independent L48 seed agrees and its chi128→192 entropy change is below 1e-6.

[numerical diagnostics] Central-window transverse order m² decreases 0.4474→0.3623→0.2890→0.2231 across the doubling sizes; the far transverse RMS and density modulation also decrease. Changing binary parity weights contributes -0.00470 to the last fractional slope, so it does not explain the excess. Both fixed-parity sequences retain excess growth. A post-hoc energy combination has correction power 1.0361, close to inverse-size behavior (same-size exact XX control: 0.9497); this is not a velocity, dynamical exponent or independent c measurement. Renyi order-2/order-3 size estimates are 1.4809/1.3702 at 96→192.

Limits. Periodic exact diagonalization at L=14,16,20 remains close to c=1, without independent confirmation of an extra bulk mode. Larger-scale crossover, boundary effects and nonrelativistic behavior remain possible. A nearby exact weighted-Dicke control demonstrates that large, clean entropy scaling alone need not identify a CFT. No conformal spectrum, thermodynamic c, or counterexample is established. The preferred R<=3 attractive alternatives weakened at larger size; no preferred-band candidate is promoted.

All local workers launched for this continuation are stopped. Under a preregistered allocation rule, the wider R=28 L192 runs were deferred at saved checkpoints after this candidate passed validation; incomplete values are excluded. Twelve winning-candidate MPS states and the two unfinished wider states have hash-verified, Git-ignored workspace copies. The unavailable pre-upgrade L96 even chi192 MPS is explicitly disclosed; its raw JSON and prior diagnostics survive. Protocols, complete null/deferred accounting, replayable state manifests, independent AI review, code, raw records and plots: [continuation contribution](contrib/d2_order_one_continuation_20261110/README.md). No cluster job or remote publication is claimed.


### 9.191 C197: the c ≈ 1.45 reading of §9.190 does not survive small changes of the field or the DM term (coordinator, Claude Fable 5, 2026-11-10; `codex/work_C197/`)

The chain of §9.190, (J_x, J_y, J_z, D, h_z) = (1/2, −1/2, 5/8, 1/10, −7/8), sits close to a special point of its family. At D = 0 and h_z = −3/4 = −√(4J_z² − 1) the chain is frustration-free, with one exact ground state for every particle number, the weighted Dicke states of the collaborator's own control (an 8-site ring has a 9-fold degenerate ground level there, checked by exact diagonalization). Such states have entanglement growing as (1/2) ln L, which an open-chain estimator reads as c ≈ 3. This suggested testing whether 1.45 is a crossover rather than a phase, by moving the field and D at fixed sizes.

Method [numerical]. The collaborator's wrapper around the C190 solver, unchanged, at L = 48 and 96, χ = 128, both parity sectors, with c_eff = 6[S(96) − S(48)]/ln 2 from their central-three and fractional-window entropies (`run197.py`, `analyze197.py`, output in `summary.txt`). The even sector is the lower one at every point except D = 0.2, and the ground-branch readings are:

| D | h_z | c_eff, central three | c_eff, fractional window |
|---|---|---|---|
| 0.10 | −25/32 | 0.43 | 0.51 |
| 0.10 | −13/16 | 2.06 | 2.02 |
| 0.10 | −27/32 | −0.33 | −0.26 |
| 0.10 | −7/8 (§9.190) | 1.418 | 1.435 |
| 0.10 | −29/32 | −0.33 | −0.34 |
| 0.10 | −15/16, −1 | 0.00 | 0.00 |
| 0.05 | −7/8 | 0.14 | 0.16 |
| 0.15 | −7/8 | 0.80 | 0.94 |
| 0.20 | −7/8 | 0.97 | 0.92 |

At §9.190's point this reproduces their 48 → 96 values (1.4185 and 1.4363) to 0.002 at the lower χ.

Reading [numerical, with a heuristic interpretation]. The reading is not stable. A field step of 1/32 moves it from 1.42 to −0.33 on both sides, while the bulk density changes smoothly (0.29, 0.27, 0.24). A negative two-size reading is not a central charge, so at these sizes the estimator is dominated by a non-universal finite-size effect, most likely the fit of the density-modulated state to the chain length. Doubling L keeps that fit, which would explain why §9.190's value repeats from 48 → 96 to 96 → 192. A c = 3/2 phase would give a plateau, and none appears. For h_z ≤ −15/16 the parity gap is size independent (0.029 and 0.088) and c_eff = 0, a gapped phase. At D = 0.2, away from both the gapped region and the frustration-free point, c_eff ≈ 1, an ordinary one-component liquid, consistent with the periodic rings of §9.190 (about 1.04) and with the single slow species found earlier at frustration-free ferromagnetic points (§9.98).

Conclusion. The §9.190 chain is most naturally read as c = 1 with large finite-size corrections near a gapped phase and a frustration-free point. This does not rule out a narrow c > 1 region elsewhere in the family, but it removes the reason to look for one at this point.


### 9.192 The χ = 1024 stage at L = 256: truncation converges, and the registered pass V1 is declared at the collaborator's bond and on the symmetric line (coordinator, Claude Fable 5, 2026-11-10; `codex/work_C190/gpu/` R12)

The registered χ = 1024 continuation of §9.178 (R8 of `codex/work_C190/gpu/RESULT.md`) finished on an H200 node in under an hour.
- [numerical] Truncation. The entropy step from χ = 768 to 1024 is 0.006 to 0.008 for the four candidate runs at L = 256, inside the registered tolerance of 0.01, so C4 passes. Residual, energy step and discarded weight pass C1 to C3 (1.0e-9, at most 1.6e-9, about 1e-10).
- [numerical] Readings on the ground branch. c_two(128 → 256) at W = 1, 6, 12 is 1.964, 1.998, 2.032 at the collaborator's bond and field (p153c), and 1.968, 1.999, 2.032 on the symmetric line at n = 3/16 (sym40_n316), up from 1.94 and 1.95 at W = 6 and χ = 768. c_prof(256) reads 2.08 and 2.09. From 64 → 128 the readings are 2.17 to 2.36.
- [numerical] Controls. The two-leg XX chain reads 1.917 and 1.898 at W = 6 over the two size pairs, and the XX chain calibrates both.
- Against the registration. Every condition of V1 holds for both keys, so by the pre-registered rule of `codex/work_C190/RESULT.md` the full qubit chain tracks the effective model's c = 2 liquid at both bonds, through L = 256. This is the first declared pass of the registered full-chain test. The other four tuned bonds stay undeclared, because their L = 256 runs stop at χ = 768.

Open. L = 512 (C195, running) and the gap test at larger χ (C194, running). The rise of c_two with χ means the finite-χ readings were lower bounds, as expected from truncation.


### 9.193 C166's L = 512 open chains: the XZ-neighbourhood reading of c ≈ 1.3 is a Néel cat, as registered (coordinator, Claude Fable 5, 2026-11-10; `codex/work_C166/typhon_out/`)

The registered prediction (`codex/work_C166/RESULT.md`) was that the unpinned open-chain S_mid approaches S_VUMPS + ln 2 from below, with S_mid(384) ≈ 1.31 and 1.34, S_mid(512) ≈ 1.35 and 1.36, bounded by 1.394 (hi_T_40) and 1.377 (dist_T_50), and falsified above the bound by 0.02.
- [numerical] Measured, `summarize.py` on the collected outputs: hi_T_40 reads S_mid(384) = 1.318 and S_mid(512) = 1.363, dist_T_50 reads 1.346 and 1.369, all with truncation below 4e-11. Both are below their bounds, close to the predicted values, and the prediction passes.
- [numerical] The two-size slope collapses at the largest sizes: 1.34 → 0.94 (hi_T_40) and 1.03 → 0.46 (dist_T_50) over 256 → 384 → 512, against the 1.3 to 1.4 that suggested a critical bulk at L ≤ 256. The staggered order correlation at r = L/2 rises toward the VUMPS value of m² (0.107 against 0.117, and 0.128 against 0.133), and edge fields along the order parameter remove most of the entropy (S_mid ≈ 0.70 at L = 256, against 1.23 to 1.28 unpinned).
- Pre-fixed criterion. C1 fails for both XZ bonds and for the C155 candidates cand_60 and cand_61, so none is a c > 1 candidate. The FES second stage fails for KI22, KI26 and N34.

Conclusion. The c ≈ 1.3 open-chain reading in the XZ neighbourhood is a gapped Néel bulk seen through a slowly relaxing cat state, as hypothesis (i) of C166 proposed. This closes the last open item of that neighbourhood.


### 9.194 C185 and C186 main jobs: the bond-changing routes to unbinding give no support beyond the effective model's U(1) part (coordinator, Claude Fable 5, 2026-11-10; `codex/work_C185/typhon_out/`, `codex/work_C186/typhon_out/`)

These extend the backfill readings of §9.172. Readouts are `summarize_typhon.py` in each folder, with output in `typhon_out/summary_coord.txt` and, for C185, the job's own `summary_eff.txt`.

C185, the transverse-DM route at the Fibonacci point.
- R6, U(1) part [numerical]. The V₂ = V₃ = 0 calibration reads c_two(64 → 128) = 2.03 at W = 6, inside its tolerance. dmv1 (|w| = 1, v = 1) passes the first-stage criterion with c_two = 1.63, 1.72, 1.86 at W = 1, 6, 12, c_prof(128) = 1.72 and truncation 3e-9, and its parity gap is flat (L Δ_P = 4.818 and 4.832 at L = 64 and 128), as predicted for a gapless pseudospin. dmv2 (v = 2) reads 1.23 at W = 6, so the registered prediction (≥ 1.6) fails for it.
- R6′, with the ΔN = 1 term [numerical]. The code calibration reproduces C173's penalty-model reading (S̄ = 1.646 and 1.886 against 1.646 and 1.869, c₂ = 2.08). At κ = 1, v = 1 and φ = π/2, c₂(48 → 96) = 1.57. The φ scan at 1.625 to 1.655 has only L = 48, since its L = 96 runs did not finish, so the amended R6′ is not evaluable.
- R7, full chain at V/Δ = 8 [numerical]. With λ = 0.3 the L = 64 runs read S_mid = 0.54, 0.69, 0.74 and 0.56 at φ shifts 0, 0.06, 0.10 and 0.14, against 1.48 for the λ = 0 reference. No shift comes within 0.7 of the reference, so there is no support.

C186, the vertex route.
- (B), block DMRG at u = 1.035 [numerical]. Only the μ offset 0.06 runs at v = 0.1 and 0.2 finished at L = 192. Both read c(96 → 192) = −0.27 with truncation below 1e-7, which is "one mode" by the registered rule. The decreasing entropy suggests a gap. The other points, including the v = 0.4 control, have L = 96 only.
- (A), ring binding at L = 48 [numerical, grid threshold]. Where N_max = 4 and 5 agree in sign, the pair is bound only at v/v_c = 0.25 (u = 0.522 and 1.035) and unbound from v/v_c = 0.5 up. The re-read against the off-grid threshold (true_threshold.py) was not done, because the Typhon records lack the field that script reads.

Conclusion. Neither bond change gives a c = 2 reading beyond the effective model's U(1) part, and the full-chain transverse-DM runs lose entropy rather than gain it. These were alternative routes to the candidate, which does not depend on them.


### 9.195 C199: two readings of the falling size sequence, with predictions for 256 → 512 fixed before the converged L = 512 result (coordinator, Claude Fable 5, 2026-11-10; `codex/work_C199/`)

The qubit-chain readings fall with size, c_two = 2.17 to 2.36 (64 → 128) and 2.00 (128 → 256), and the interim 256 → 512 reading of C195 is 1.58 at matched χ = 1024, not converged. Two pictures fit the converged points, and they predict different things at the next size. `flowfit.py`, output in `flowfit.out`.

- (A) c = 2 with a correction from an irrelevant operator [heuristic]. The deviation from 2 shrinks by a factor of 5 to 100 per doubling (+0.17 to +0.36 at 64 → 128, −0.04 to +0.03 at 128 → 256). The prediction for 256 → 512 is 1.99 to 2.01, within ±0.04 at every window.
- (B) a slow flow toward a pair liquid [heuristic]. The decline per doubling (−0.17 to −0.34) continues, giving 1.73 to 1.83 at W = 6 for 256 → 512 and 1.46 to 1.66 for 512 → 1024.
- K_s, a consistency check [heuristic, from §9.182]. With the pair estimate, K_s changes by −0.32 to −0.62 from 64 → 128 and by −0.08 to −0.09 from 128 → 256, a slowdown by a factor of 0.13 to 0.28 per doubling. In a Kosterlitz–Thouless flow on the irrelevant side the change per doubling scales as 2^(−4(K−1)), and the observed factor corresponds to K ≈ 1.46 to 1.74, which matches the measured K_s (1.4 to 2.0 at L = 256). A flow toward a gap would make the change per doubling grow as K_s approaches 1, and the data show it shrinking. If the slowdown continues, K_s moves by a further −0.01 to −0.04 in total.

Fixed now, before the converged L = 512 result. With S̄(256) at χ = 1536 (R13 of `codex/work_C190/gpu/RESULT.md`, which may lower the reading by up to about 0.03), a converged c_two(256 → 512) at W = 6 of 1.93 or more is read as favouring (A), 1.85 or less as favouring (B), and anything in between as undecided. The same reading is applied to the symmetric-line bond (C195b) when it finishes. C198 tests whether the open-chain sequence of the effective model, whose infinite chain is c = 2, has the same shape.
