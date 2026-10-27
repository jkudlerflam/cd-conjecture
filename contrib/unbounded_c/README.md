# Large c at fixed d: block constructions and where they stop

## Claim

Class: one-site translation-invariant, nearest-neighbour, Hermitian chains H = Σ h_{i,i+1} on (ℂ^d)^⊗L with periodic boundary conditions, whose low-energy spectrum on compatible ring lengths is a unitary CFT (possibly several, with different velocities) of total Virasoro central charge c. Questions: is c ≤ f(d) for some f, and is there a chain with c > d − 1?

Neither is settled here. The findings are these.

1. The literature has no construction and no bound. All fixed-coupling translation-invariant constructions produce decoupled blocks, and their gapless part, when present, is a decoupled factor that is not a large-c CFT (§1).
2. Blocks are free at fixed d. C7 Theorem 4 selects any period b at d = 4 with L-independent couplings. Letting the walker carry no spectator gives cells with a q^{b−1}-dimensional register at d = 3q + 1 (§2).
3. If the internal Hamiltonian of a cell is unrestricted, every bond-rank obstruction disappears. G1 (Astra) proves that a swap P_N between N-dimensional cell codes arises at second order through a port of any dimension d ≥ 2, even through a rank-one bond, for every N, with vanishing first order and no one-cell terms. At first order, scaled Clifford compressions give r ≤ q² − 1 Majorana channels per bond and c = r/2 (§3).
4. Combining 2 and 3 gives a candidate counterexample to c–d at d = 22 (q = 7) with c = 24, and, through G1, a candidate family with unbounded c at fixed d. Both rest on one unproved step, realizing the prescribed cell (code, port compressions, and for G1 the port spectral weights) as the ground space of a nearest-neighbour Hamiltonian built from the same h. I call this local realizability (§4, §5).
5. Protecting the infrared theory is not the obstacle at large c. Narain CFTs with no primary of dimension ≤ 2 besides the vacuum exist for every c ≥ 31, so they need no symmetry (§5).

## Status

Heuristic. Items 2 and 3 are proved (C7, G1, and the Clifford dilation argument, which is elementary), and I checked them. Item 4 is a construction outline with the open step stated. It would be upgraded to a counterexample by an explicit nearest-neighbour cell with the Clifford port compressions of §4, even at small q as a proof of principle.

## How it was established and checked

- Author: Claude (Anthropic, Opus 5.5), as a sub-session of the coordinating session, 2026-10-27.
- AI assistance. Two Sonnet 5.5 literature agents read the papers in §1, tagged [read, agent]. C7 is the coordinator's closed-book gpt-6-astra run (`codex/work_C7/C7_result.md`). G1 is a closed-book gpt-6-astra (high) run I set up (`spec_G1_port_exchange.md`, output in `G1/`).
- Checked by me. The selector by my own exact diagonalization at a coupling different from C7's. G1's second-order formula by hand (Eqs. 17 to 19 of `G1/G1_result.md`) and by rerunning its dense check for d = 2, N = 3. The Siegel–Weil count by my own script, which agrees with the literature agent's numbers and with the ball-volume asymptotics. I spot-checked the BCLP and Harley et al. abstracts but did not reread the other papers.

## Argument

### 1. What the existing constructions guarantee

- Gottesman–Irani (0905.2419) [read, agent, §5.8, §7]. A finite chain holds one history state over the whole chain with gap Ω(1/N³), a global quadratic clock mode. On the infinite chain energetics selects independent segments of length N, but the couplings carry Θ(log N) bits and a wrong length costs only Ω(1/N⁹) per site.
- AGIK (0705.4077, d = 12) and Hallgren–Nagaj–Narayanaswami (1312.1469, d = 8) are not translation invariant and hold one global history state [read, agent].
- Bausch–Cubitt–Lucia–Perez-Garcia (1810.01858) [read, agent, Thm 1, Lemma 21, Thms 22 and 25]. Couplings are fixed and the marker Hamiltonian is block-diagonal in the marker pattern. Segments of the halting length are selected and the ground state is a product over them. A leftover segment costs O(1), and periodic boundary conditions at fixed d are left open. The gapless part is a decoupled factor H_dense, which in the 2D papers is XX + YY + 2Z [read, agent, 1502.04573 §2.1.2], with dispersion 4 + 4cos k (my check), hence z = 2 and not a CFT.
- Kohler–Piddock–Bausch–Cubitt (2003.13753) [read, agent, §3.3, §4]. Translation-invariant 1D universality runs one history state over the chain, with energies exp(poly n) and global operator-norm errors. Zhou–Aharonov (2102.02991) reach poly(n) only without translation invariance [read, agent].
- Harley et al. (2306.13739) prove no-gos for superlinear encodings and for locality reduction at global precision η → 0 [read, agent, Prop. 10, Thm 26]. Cells with fixed local error avoid both.
- Entanglement. Gottesman–Hastings (0901.1108, d = 9, not translation invariant) reach S ≳ (Δ log 1/Δ)^{-1/4}, and Irani (0901.1107, d = 21) reaches S = Ω(min(m, n)) with gap ~ n^{-12} [read, agent]. At fixed d the area-law constant therefore grows polynomially in the correlation length, far beyond (c(d)/6) log ξ, but none of these chains is conformal. Colored Motzkin and Fredkin chains have z ≥ 2 [read, agent, 1408.1657, 1707.02317], and Gosset–Mozgunov (1512.00088) make every gapless frustration-free chain non-conformal.

Nobody has used these tools on a central charge. The encodings are global (one clock, gaps 1/poly(L)), or local with blocks that never interact.

### 2. Cells at fixed d

I wrote down C7's selector independently, with letters (M, L, W, R) for (#, L, X, R). Exact diagonalization on a ring of L = 10 at t = 1, K = 8, μ = 0.9 gives a five-fold ground level −1.436068 = 2(μ − 2cos(π/5)), so period 5 > d = 4, and the next level −1.346264 from tilings 4 + 6, both as predicted (`checks.py`). Over μ = 0.9 to 1.999 the period runs from 5 to 172, close to π(3t/(2t − μ))^{1/2}, with gaps t/p⁴ to other cell lengths and 3π²t/p² inside a cell.

For a register, attach ℂ^q to M, L and R but not to W, so d = 3q + 1. The hop LW → WR then carries the spectator from one site to the next without reading it. The selector analysis is unchanged, and each cell has ground space (walker ground state) ⊗ (ℂ^q)^{⊗(b−1)}, with the spectators of the last R and of the next M as ports. This improves the 4q of C7 Corollary 5.

### 3. Unrestricted cells defeat every rank obstruction

First-order compression cannot raise the connected rank of a bond above d² − 1 (C7 Theorem 6). Two ways around it exist once the internal cell Hamiltonian is unrestricted.

Second order (G1). Let the port operators F_a map the code into orthonormal syndrome spaces and choose the cell's excited spectrum so that the matrix-valued port spectral weights, paired by the Cauchy kernel 1/(E_e + E_f), whiten into a Hilbert–Schmidt basis. Then H^(1) = 0 and H^(2) = a·1 − sτ² P_N exactly, with either sign, for every N and every d ≥ 2, even with a rank-one bond (G1 §5). The cell dimension grows like N³ and the cell Hamiltonian is a generic operator on it. G1 verified (d, N) = (2, 3), (2, 4), (3, 4) by full diagonalization, with splittings ∝ λ² and non-swap terms starting at λ³. So there is no lower bound on the order in terms of log N, and the lower-order terms need not be configuration dependent. My simpler example (`checks.py` part d) reaches only a diagonal coupling of rank N − 1.

First order. A unital completely positive map M_q → M_D with F_α ↦ s γ_α, where γ_1, …, γ_r anticommute, exists for s small enough (Choi positivity, for example s ≤ 1/3 for the three Pauli matrices), and Stinespring turns it into a code. A chain of such cells gives a free Majorana chain if three things hold. Each cell code is a graded Clifford module V_L ⊗ V_R. The left port compresses to s γ^L_α. The right port compresses to s iγ^R_α P, with P the cell parity. Then Jordan–Wigner makes the inter-cell coupling Σ_α F_α ⊗ F_α into J Σ_α iγ^R_α(j)γ^L_α(j+1), and an intra-cell term κ Σ_α iγ^L_α γ^R_α completes r Majorana chains. At |κ| = |J| every flavour is critical, and c = r/2 ≤ (q² − 1)/2.

The grading matters. An AKLT edge also compresses S^a to anticommuting Pauli matrices, but there is no P, and the chain of edges is a spin-½ Heisenberg chain with c = 1, not SO(3)_1 with c = 3/2.

### 4. A candidate counterexample to c–d at d = 22

The construction takes d = 3q + 1 with q = 7, so d − 1 = 21 and the target is c = (q² − 1)/2 = 24. The ingredients follow.

(i) The selector of §2 with at least 18 spectators per cell (18 log₂7 > 48), enough to hold a 2^48-dimensional code. [proved, C7 and §2]

(ii) A port bond F_α ⊗ F_α = qP_q − 1 on the letter pair (R, M), which is isotropic in the channel index. [explicit]

(iii) A cell code of §3 on the spectator register, with r = 48 channels at each port. It must be the ground space of the intra-cell spectator Hamiltonian, which is built from h on the pairs (M, L), (L, L), (R, R) and weighted by the walker profile. [open]

(iv) Criticality. A generic Majorana flavour needs one tuning, so r = 48 tunings in all, or one if an exact symmetry acts irreducibly on the 48 channels. Parameter counting is no obstacle, since h has d⁴ ≈ 2.3 × 10⁵ real parameters, and c–d does not require stability. Higher orders add quartic, Gross–Neveu-type marginal terms, which must lie in the marginally irrelevant cone. [open]

(v) Commensurability. On rings with b dividing L the spectrum is b degenerate copies of the CFT, one per translate of the cell pattern. On other rings an immobile defect cell is forced, and the defect's position gives an exact degeneracy of order L/b, which is not a CFT spectrum. So the claim would be that H_L violates c–d along L ∈ bℤ. Whether that is in scope depends on reading "the low-energy spectrum is a collection of unitary CFTs" for all large L or along a subsequence. ULS itself needs N | L for a singlet ground state but has conformal spectra on every ring, so this reading is an extra hypothesis Jonah should rule on. Making defects mobile gives non-conformal defect bands instead.

A symmetry that protects isotropy runs into rigidity. If an on-site group acts on ℂ^q as a unitary 2-design (SU(q), or the qudit Clifford group), its only invariant bonds are a + bP or a + bΦ. That forces exact SU(q), and then the code must be the spinor of SO(q² − 1) restricted to SU(q). The infrared theory would be SO(q² − 1)_1 ⊃ SU(q)_q, a conformal embedding with c = (q² − 1)/2, and it would also break the level ceiling k ≤ ind(R) = 3 of REPORT §3.4. I see no way to get that spinor edge from an inhomogeneous ULS or Temperley–Lieb cell, so the realistic version of (iv) tunes 48 masses without symmetry.

The second-order route of G1 is cheaper in d. With a qubit spectator (d = 7), c = N − 1 > 6 needs only N ≥ 8, but the target is SU(N)_1. Its N² − 1 currents give marginal SU(N)-breaking current-current terms, arising here at third order, which can flow away over a length exp(O(1/λ)).

### 5. Unbounded c, protection, and local realizability

For unbounded c at fixed d, G1 cells give ULS on N-dimensional logical sites for every N, and a robust target removes the protection problem. A Narain CFT whose lattice has no nonzero vector with Δ = (p_L² + p_R²)/2 ≤ 2 has no relevant operator. Its only marginal operators are the J^iJ̄^j moduli and the spin-two bilinears that change velocities, and both keep a multi-velocity Gaussian theory of the same c. By the Siegel–Weil average (2006.04839, Thm 4.1 [read, agent]) the expected number of such vectors is 4.11, 1.69 and 0.67 at c = 30, 31 and 32 (`checks.py`). They come in ± pairs, so such lattices exist for every c ≥ 31. A Hermite-constant bound requires c ≥ 13 [read, agent]. The remaining caveats are that this stability is RG lore, not a lattice theorem, and that the leading-order logical model must already sit in that strongly interacting basin.

Everything therefore hinges on local realizability. Can a gapped nearest-neighbour Hamiltonian on a segment of fixed-dimensional sites, built from the one translation-invariant h, have a prescribed code as ground space, with prescribed compressions (and, for G1, prescribed spectral weights) at its end sites, up to small errors? Abstractly the data exist (G1 Lemma 1, Stinespring). A natural route is to run a translation-invariant universal simulation (KPBC type) inside each selected cell, since a cell is a finite target and its energies and errors stay fixed as L → ∞. The cells would then talk to each other only through their port sites. Nobody has shown that such a simulation preserves end-site compressions, but I see no reason it would not. The cost is a very large d0 for the c–d counterexample, because q must outgrow the number of simulator tracks.

### 6. Verdict

Neither verdict is reached, but the obstruction has moved. Arbitrary periods, large cell manifolds, rank raising across a single bond (first order with Clifford ports, or second order for arbitrary N), and infrared protection at c ≥ 31 are all available or proved. The single unproved step is local realizability of a prescribed finite cell inside a translation-invariant chain. If it holds, c–d fails at some finite d0 (22 in the most economical version, much larger through generic simulation), and c is unbounded at fixed d. I now lean towards that outcome, with modest confidence, because no principle I can find forbids it and every specific obstruction proposed so far (classical period, bond rank, perturbative order, symmetry protection) has fallen.

## Reproducing

- `checks.py` (numpy, scipy, fractions): (a) exact diagonalization of the selector, L = 10; (b) period and gaps against μ; (c) Siegel–Weil counts; (d) the second-order diagonal rank example. About four minutes, almost all in (a).
- `spec_G1_port_exchange.md` is the closed-book spec. `G1/G1_result.md`, `G1/gadget.py` (runs in under a minute) and the JSON files are what gpt-6-astra (high) returned. `python3 G1/gadget.py --cases '2,3' --lambdas 0.01,0.005 --method dense` reproduces my rerun.

## Relation to the rest of the repository

This builds on C7 (`codex/work_C7/`) and addresses Jonah's "any f(d)" question and the refocused question of a c–d counterexample. It does not change the top-level status line, but if the coordinator agrees "Where things stand" could note that block selection and rank raising are no longer obstacles, that a candidate counterexample at d = 22 reduces to realizing one finite cell locally, and that the commensurability reading in §4(v) needs a decision.
