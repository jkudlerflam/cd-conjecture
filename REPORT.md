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
