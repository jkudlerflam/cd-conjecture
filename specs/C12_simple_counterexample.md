# C12: Find a SIMPLE explicit counterexample to c ≤ d − 1 at moderate d

Closed book: do not read or list any files you did not create. Work only inside the current working directory.
Python (numpy/scipy) allowed; keep total runtime under about 100 minutes. The standard is theoretical physics
(standard RG, bosonization, controlled perturbation theory, numerics calibrated with at least two estimators), not
mathematical rigor. Be honest: report negative results clearly. Write PROGRESS.md as you go and RESULT.md at the end.

## The conjecture
Latorre–Sierra: a critical, one-site translation-invariant, nearest-neighbour Hermitian chain H = Σ h_{i,i+1} on
(C^d)^{⊗L} (periodic boundary conditions) whose low-energy spectrum is a collection of unitary CFTs has total central
charge c_T = Σ c_a ≤ d − 1, saturated by the SU(d) swap chain (SU(d)_1).

## What is known (use it)
- No counterexample is known at d ≤ 7 in symmetric families. At d = 2 extensive searches find c ≤ 1. At d = 4 the
  candidate Sp(4)_1 × tricritical Ising (c = 16/5) on spin 3/2 has codimension 2 and was not found numerically.
- A complicated construction suggests c is unbounded at fixed d0 ~ 10² to 10⁴. Sites spontaneously form long cells
  (a d = 4 period selector: letters #, L, X, R with a quantum walker whose kinetic energy −2t cos(π/b) selects the
  cell length b). Cells run internal history-state computations preparing logical codes, and couple through one
  small bond. The key reason c can exceed d − 1 is that a cell's logical space grows with its length, while the
  interaction through a single bond only needs a modest operator rank. Fermion-like channels cost operator rank 2
  per unit of c at first order, while SU(d)_1 costs rank d² − 1 for c = d − 1.
- Free-fermion NN chains satisfy c ≤ log₂ d. Directly projecting one bond between block codes cannot produce an SU(N)
  swap with N > d (rank argument). At second order, any logical bond can be produced through qubit ports, but the
  cell Hamiltonian is then generic.
- Spin–charge separation (t–J, Ogata–Shiba) gives c_max(d₁ + d₂) ≥ c_max(d₁) + c_max(d₂) + 1. It saturates but does not
  exceed d − 1.

## Task
Find the SIMPLEST explicit NN translation-invariant Hermitian chain with c > d − 1, ideally at d ≤ 64 and better
d ≤ 16, by any mechanism. Ideas to explore (and invent others):
(a) cells with FRUSTRATION-FREE cores whose ground space has many boundary-accessible states (no computation track),
    for example AKLT-like or MPS parent Hamiltonians whose edge spaces grow with the cell length, or constrained
    (Rydberg, Fibonacci) cells;
(b) squeezed-space constructions, where mobile objects carry large internal spaces assembled from several sites;
(c) emergent gauge constraints or a non-on-site symmetry producing many gapless flavours from few states per site;
(d) multi-leg ladders encoded in a chain by spontaneous period formation;
(e) a direct numerical search (optimization over h at small d, say d = 3 or 4, targeting c above d − 1 with
    calibrated estimators), if you think a structured family is promising.
For each idea, decide quickly whether it can beat d − 1. For any candidate, give h explicitly, estimate c with two
calibrated estimators at the largest sizes you can reach (exact diagonalization, DMRG), and check z = 1. Report the
best candidate, or explain why simple mechanisms fail and what minimal complexity seems necessary.
