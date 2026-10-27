# C9: Realize a programmable cell inside a translation-invariant nearest-neighbour chain of fixed local dimension

Closed book: do not read or list any files you did not create. Work only inside the current working directory.
Python (numpy/scipy) is allowed; keep total runtime under about 90 minutes. Work hard and be honest: a precisely
located obstruction is valuable, an overclaimed construction is not. Flag every step you assume rather than prove.
Write PROGRESS.md as you go and RESULT.md at the end.

## Why this matters
The c–d conjecture (Latorre–Sierra) says a critical, translation-invariant, nearest-neighbour Hermitian chain on
(C^d)^{⊗L} has total central charge c ≤ d − 1. Whether c is bounded by any f(d) is also open. Two ingredients for a
counterexample, possibly with unbounded c at fixed d, are now established, and one is missing.

(1) Spontaneous cells at fixed d. A d = 4 translation-invariant NN chain (letters #, L, X, R; allowed pairs
#L, #X, LL, LX, XR, X#, RR, R#; Hamiltonian K Σ Q_i + μ Σ |#><#|_i − t Σ (|LX><XR| + h.c.), with Q_i the projector
on forbidden pairs, K > 4t, μ = 2t − δ) has exactly conserved marker positions. Between markers a single walker X
hops on a path, giving cell energy μ − 2t cos(π/b), and choosing δ between crossing values selects any period
b ≥ 5 uniquely, with couplings independent of L and a gap for each fixed b. Tensoring each site with a spectator
space gives cells carrying large degenerate registers (a spectator qubit gives 2^b states at d = 8; letting only
the walker carry nothing gives d = 3q + 1 for a q-dimensional spectator).

(2) Exchange through a tiny port at second order (proved and numerically verified). Let a cell be C^D = port C^d ⊗
auxiliary K with code G: C^N → C^D, G|i> = d^{-1/2} Σ_{a=1}^d |a>_port |a,i>_K. Then G^†(O⊗1)G = (tr O/d) 1_N, so
every traceless port operator compresses to zero (no first-order term), and S_a = (F_a ⊗ 1)G (F_a a traceless
Hermitian basis, tr F_a F_b = d δ_ab) are orthogonal isometries into excited "syndrome" spaces. A lemma lets the
cell Hamiltonian H_A have kernel exactly ran G and arbitrary prescribed excited spectral weights S^† P_e S = M_e
(positive, summing to the identity) at energies E_e, provided D − N ≥ (number of levels) × (d² − 1) N. With the
port interaction h = Σ_a F_a ⊗ F_a = d P_d − 1 between the boundary ports of neighbouring cells and suitably
chosen M_e, E_e (a Cauchy-matrix argument), the second-order effective Hamiltonian between two cells is exactly
a I + b P_N, b ≠ 0, for EVERY N ≥ 2 and every d ≥ 2; first order vanishes. Explicit examples (d, N) = (2,3), (2,4),
(3,4) need cell dimensions D = 22, 40, 69 and were checked by full diagonalization (splitting ∝ λ², residual ∝ λ³).
A chain of such cells is therefore, at second order, the SU(N) Uimin–Lai–Sutherland chain (if b > 0), whose IR is
SU(N)_1 with c = N − 1.

(3) MISSING: the G1 cell Hamiltonian H_A is an arbitrary operator on C^D, not something built from a fixed
nearest-neighbour h. The task is to close this gap.

## Task
Construct, or prove impossible, a family of translation-invariant NN Hamiltonians h_N on (C^{d0})^{⊗L} with FIXED
d0 (independent of N), whose low-energy physics is that of cells (of length b(N), spontaneously selected) each
realizing the G1 gadget for N (or another construction with the same effect), coupled through the boundary bond.
Couplings may depend on N (the conjecture allows any fixed h, of any norm), but not on L. Address each requirement:
(a) Inside a cell, the sites see the same h, but the cell's internal configuration (marker, walker position,
    letters) can make terms effectively position-dependent. Use this, or a Feynman–Kitaev clock, history states,
    perturbative gadgets, or translation-invariant universal-simulation ideas (finite cells, so energy scales
    may be large but finite for each N), to simulate the prescribed H_A on the cell's low-energy subspace to
    accuracy ε, including its excited spectral structure as seen from the port, which the second-order result
    uses.
(b) The cell's code must be degenerate (or split far below the second-order coupling), port compressions must
    vanish at first order, and all cell machinery (walker, clock, gadget excitations) must be gapped above the
    inter-cell energy scale, so that no extra gapless modes appear.
(c) Error control. Corrections of order ε and higher orders in the port coupling produce SU(N)-breaking and
    longer-range terms in the effective cell chain. SU(N)_1 has relevant SU(N)-breaking perturbations in
    general. Either show the corrections preserve enough symmetry (for example by building SU(N) or a subgroup
    into the cells exactly), or switch to a target CFT with no relevant operators (Narain CFTs with no primary of
    dimension ≤ 2 other than the vacuum exist for every c ≥ 31), and say what is then needed.
(d) Commensurability: on rings with L not a multiple of the period, an immobile defect appears. State the
    result for L in bZ and what happens otherwise.
(e) Give d0 explicitly, and state what this implies: a counterexample to c ≤ d0 − 1 once N > d0, and c
    unbounded at fixed d0 if N is arbitrary.
Also consider the alternative first-order route: cells whose two boundary sites compress to many mutually
anticommuting SCALED Majoranas (scaling escapes the multiplicative-domain limit log₂ d; for a qubit, three
anticommuting Majoranas at scale 1/3 are completely positive), with a fermion-parity grading on the code, giving
an effective critical Majorana ladder with c up to about (q² − 1)/2 at d ≈ 3q + 1 to 4q, which beats d − 1 for
q ≳ 7. Compare the two routes.

## Deliverable
RESULT.md: the verdict (construction with d0, or obstruction), proofs or detailed constructions with all energy
scales, small-size numerical checks of any cell you build (exact diagonalization of one or two cells), and the
precise remaining gaps. Save code in this folder.
