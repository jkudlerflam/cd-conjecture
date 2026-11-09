This run was done by Claude Opus 5.5 (claude-opus-5-5), not by Codex. Astra was not used (out of credits).

# C169 result: constructive counterexample hunt

Derivations are in PROOF.md, the log in PROGRESS.md. Scripts (all local, a few seconds each): frustration.py, claw.py, interface.py, freecut.py, with outputs in the matching .out files. No Typhon job and no DMRG were run.

## Verdict

No candidate. Every construction fails, and the reason is stated precisely below. Exact maps cannot create a counterexample, since c is invariant and α^{-1}(B_2) is B_2 relabelled. They can only recognize a known c > 1 model, and none of the known ones is equivalent to a B_2 bond under on-site unitaries (any period), Clifford dualities (KW, JW, cluster entanglers, any Clifford QCA) or blocking. Each class has a computable invariant that separates every target. Asked whether the obstruction is range, translation by one site or dynamics, the answer is range in every free, spin-wave and exact-decoupling setting, and these bounds hold at every period, so one-site translation is never the obstruction. Dynamics is the obstruction in the interacting settings: the sign of the induced XZ interaction, the SU(2)_1 multiplet at spin 1/2, and one species per site. The three irreducible-CFT targets added by the coordinator (Z_5 parafermion, Potts × half-floating boson, repulsive-host c = 2) are all blocked by derived dynamical statements in the controlled regions and fall back on (A6)/(A_W) outside them.

## Key claims

Exact maps.
- [proved] Products of single-site unitaries, site-dependent with any period, preserve the support set of the Pauli decomposition. No target, all of which have range ≥ 2, becomes nearest-neighbour. Equivalence to a one-site chain holds iff the chain has a twisted translation (C156's screw lemma, generalized).
- [proved] Clifford invariants. For translation-covariant maps the diameter of the anticommutation offset set D(F, G) is invariant. In B_2 it is at most 2 for all 144 family pairs, against 3 or 4 for cluster-Ising, RZFA, O'Brien–Fendley, Titvinidze–Japaridze, J1–J2 and zigzag two-XX. For any Clifford map the number of infinite frustration components is invariant. In B_2 it is at most 2 over all 4095 Pauli-aligned subsets, against 3 (cluster-Ising, Ising × XX) and 4 (two XX legs).
- [proved] The KW-stable part of B_2 lies inside the free locus. The cluster entangler acts on its nearest-neighbour-preserving subspace as the site reflection, so it gives nothing new.
- [proved] Blocking (any block size, arbitrary block unitaries). The cut must be a qubit interface: dim V_L ≤ 3, Clifford relations, and triple product ±i. Every target fails (interface.out). Cluster-Ising fails only the triple product, since its cut algebra is M_2 ⊕ M_2.
- [proved, elementary] General non-Clifford QCAs carry no invariant beyond the spectrum. [literature, derived] Integrable c > 1 targets are still excluded for bonds with symmetric J, because their images would carry local charges, and Yamaguchi–Chiba–Shiraishi (arXiv:2411.02162) prove that only the known c ≤ 1 integrable chains do. [assumed] For bonds with DM this step is not covered.

Range against translation.
- [proved] Free fermions at any period p: z² det h(z) has degree ≤ 4, so c ≤ 1 (≤ 1/2 if one bond in the cell has rank-1 in-plane block). Range r gives c ≤ r, saturated by Titvinidze–Japaridze. This extends C143's free-locus bound to bonds that are free only after a periodic site-dependent rotation (for example XZ as the Kitaev chain). freecut.py: span ≤ 4 and ≤ 4 unit roots in 3800 random chains.
- [proved, enumeration] No free fermion in disguise in Pauli-aligned B_2. All 20 claw-free classes are ordinary JW-free.
- [proved algebra, derived physics] Linear spin waves about any periodic product state give at most two soft momenta and a single magnon minimum in the number-conserving case. So the J1–J2 two-minimum c = 2 mechanism needs range 2 at every period. This extends C154 Lemma 3 and C160.
- [proved] Exact term-wise decoupling H = H_A + H_B, at any period, occurs only for α XZ + α′ ZX up to rotation: two Majorana chains, c = 1. This is C157's "XZ is the unique decoupled second sector", proved in strong form.
- [proved / derived] Ising × XX needs three commuting Majorana chains and B_2 allows two. Through XZ, the needed B self-coupling b_jb_{j+2} is a four-site Pauli string, so it exists only as a virtual-A effect, whose second-order sign is attractive (C157). The obstruction is range at first order and dynamics at second.

Coordinator targets.
- Z_5 parafermion at a melting tongue (C174 R5).
  - [proved] B_2 has no on-site Z_m (m ≥ 3) beyond U(1), so Z_5 is translation acting on a locked tongue.
  - [derived] No lattice symmetry, in any stratum, maps kinks to antikinks. Reflection acts as CP and time reversal as K, and both leave the chiral clock invariant. So the chiral term is never symmetry-forbidden and the FZ locus has codimension 3.
  - [derived] Lagrange's identity gives Δ(O^m) Δ(W) ≥ m²/4 for the clock vertex and any disordering vertex, so for m ≥ 5 they cannot both be relevant. In every Gaussian-controlled region the tongue melts through a c = 1 floating phase.
  - [derived] At spin 1/2 the Z_5 clock vertex is third order in charge-≤2 harmonics and relevant only for K > 25/8, the semiclassical h_F regime where quantum order–disorder melting is excluded (C160). The φ-side period-5 vertex has Δ = 25K > 25/4 because K > 1/4 for every nearest-neighbour U(1) liquid.
  - [literature] Period-p Rydberg chains (p ≥ 5) melt in two steps.
  - An FZ point could exist only at order-one harmonic amplitude, inside the (A6) region. It is a frontier point of the c = 1 floating component and is not excluded by the closure lemma alone (item for C174).
- Potts × half-floating boson, c = 9/5 (C171 R1).
  - [derived] Chirality is allowed in T and K_I too, since Θ acts as K, so the Potts locus has codimension 2 there, and at least 3 in N (where a relevant ε∂_xφ coupling is also allowed). R1's codimension 1 should be corrected.
  - [derived] With one qubit per site a period-3 order is either the host's own locked spiral angle, which gaps the host, or a φ-side up-up-down wave with Δ = 9K > 9/4, which does not exist to all orders in weak harmonics. Classical nearest-neighbour Ising order has period ≤ 2 [proved]. So the Potts sector needs a second species, which is (A_W)(b) or (A6).
- Repulsive-host c = 2, K_1 = 0.9, K_2 = 0.2 (C171 R4).
  - [derived] Kinematically fine. K_2 < 1/4 cannot be a single-species qubit fluid, so it must be n ≥ 2 composites, which is exactly the gap in C163 Lemma L. The conversion vertex has dimension 1.53 at n = 1 (locked) and 2.36 at n = 2 (free).
  - [derived near U(1), assumed beyond] In a repulsive host the 2-string branch lies above the continuum and never softens (C163 Fact 3). Beyond the U(1) neighbourhood this is (A6).

## What is left open

1. Interacting nearest-neighbour qubit chains with a period p > 1 (a superset of B_2). No known example has c > 1, and every controlled bound above holds at any period. If a period-2 chain with c > 1 were found, any proof for B_2 would have to use one-site translation. A structured period-2 search was not run (the spec forbids repeating sampling in B_2, and this is a different space), and it is the natural next numerical test of whether translation matters.
2. Non-Pauli-aligned bonds in the free-in-disguise enumeration (§5.2), and decompositions into commuting sectors whose terms are longer than nearest-neighbour (§7).
3. The integrability step of §4 for bonds with DM.
4. Order-one harmonic amplitudes, the common home of the FZ, Potts × boson and pair-fluid routes, are the same territory as (A6), (A_W) and A_S. No construction there was found, and none was excluded. A concrete test family for the Z_5 route, if one is wanted: attractive parent Δ′ ≈ −0.8, spiral pitch 2π/5, and charge-1 (B_x) and charge-2 (XX − YY) harmonics of order one, with the width of the floating window on the commensurate line measured as the amplitude grows. A pinch would be the signal.

Sources: [Yamaguchi, Chiba, Shiraishi, arXiv:2411.02162](https://arxiv.org/abs/2411.02162); [Samajdar, Chepiga, Mila, PRR 3, 023049](https://doi.org/10.1103/PhysRevResearch.3.023049); [Chepiga, arXiv:2308.07144](https://arxiv.org/pdf/2308.07144); [arXiv:2606.16128](https://arxiv.org/pdf/2606.16128); [arXiv:2401.08087](https://arxiv.org/pdf/2401.08087); [arXiv:2412.11201](https://arxiv.org/pdf/2412.11201).
