# C10: Write out the segment-confined computation track explicitly
Closed book: do not read or list any files you did not create. Work only inside the current working directory.
Python (numpy/scipy) allowed; keep total runtime under about 90 minutes. Work hard and be honest: a precisely located
obstruction is valuable, an overclaimed construction is not. Flag every step you assume rather than prove. Write
PROGRESS.md as you go and RESULT.md at the end.

## Context (given; do not re-derive)
The c–d conjecture (Latorre–Sierra) says a critical, translation-invariant (one-site), nearest-neighbour Hermitian chain
on (C^d)^{⊗L} has total central charge c ≤ d − 1. A construction is being assembled to show c can exceed d − 1, and even
be unbounded, at a fixed local dimension d0:
(1) Cell selector (d = 4, letters #, L, X, R; allowed pairs #L, #X, LL, LX, XR, X#, RR, R#; H = K Σ Q_i + μ Σ |#><#|_i
    − t Σ (|LX><XR| + h.c.), Q_i the projector on forbidden pairs, K > 4t, μ = 2t − δ). Markers # are conserved; a
    single walker X hops between markers, so a cell of length b has energy μ − 2t cos(π/b), and δ selects any period
    b ≥ 5 uniquely, with gaps of order t/b² (walker) and t/b⁴ (wrong lengths). Couplings are independent of L.
(2) Port-idle history cells. Each site also carries a computation-track letter and a q-dimensional register. Inside
    a cell a standard-form history-state computation (clock t = 0..T, gates W_t on workspace registers only, energy
    terms h_t at chosen clock times) prepares and maintains a logical code. If no gate ever acts on the PORT register
    (the register on the boundary sites next to #), every port operator commutes with the circuit and maps the
    history subspace into itself exactly, so port compressions are exact and the effective coupling between cells
    converges to the target as 1/J_c (clock coupling). Every cell has a clock tower J_c[1 − cos(kπ/L)].
(3) Route B cell (first order). Port p and partner s (each q-dimensional, q = q₊ + q₋, with a grading under which
    the traceless port operators F_a, a = 1..r, r = 2 q₊ q₋, are odd), |Ω> a fixed maximally entangled (p, s) state,
    |Ω_a> = (F_a ⊗ 1)|Ω>, a logical register with Clifford generators γ_a, and B = Σ_a |Ω_a><Ω| ⊗ γ_a. The cell term
    H_L = Δ(1 − |Ω><Ω|) + g(B + B†) has an exactly 2^n-fold degenerate ground space for every g (each logical ψ spans
    with r^{-1/2}Bψ a block [[0, g√r], [g√r, Δ]]). Both ends (with a mirror term using η_a = iγ_a P, P the logical
    parity) compress exactly: G†F_a G = α γ_a on the left port and α η_a on the right, with tan 2θ = 2g√r/Δ and
    α = −sin 2θ/√r. The inter-cell bond λ Σ_a F_a ⊗ F_a (right port of cell j, left port of cell j+1) then gives at
    first order the free Majorana chain λα² Σ_{j,a} i γ_{j,a} γ_{j+1,a} with ONE Majorana per cell per flavour, which
    is critical with c = r/2 (right movers at k = 0, left movers at k = π).
(4) Route A (second order, G1/C9) gives an exact SU(N) swap between cells for any N through qubit ports.

## Task
Write out explicitly, as a translation-invariant nearest-neighbour Hamiltonian on a finite alphabet, a computation track that lives inside each selector cell (between consecutive markers #), is blind to the walker letters L, X, R (it must commute with the selector), and implements a port-idle history-state computation for a prescribed finite program (gate list on workspace registers plus energy terms h_t). The program for route B must prepare the Bell pair |Ω> on (port, partner) as an energy term (never a legality constraint) and implement the coupling B through a finite gate set acting on partner, ancilla and logical registers only, keeping the port idle.
Requirements:
(a) Give the alphabet (its size k_c), the allowed pairs and the transition rules, and the full bond h, explicitly enough to code.
(b) Prove that the legal, correctly initialized sector's ground space is the history space of the program, and bound the gap of every illegal sector (no head, several heads, cycles, histories cut by a marker, wrong cell lengths) below by an explicit function of the couplings. The selection margin of the selector (of order t/b⁴) must not be destroyed.
(c) Explain how the cell length b, or a coupling phase, determines the program (for example a unary or binary counter decoding N from b), and how k_c stays independent of N.
(d) Verify a small instance numerically by exact diagonalization (one or two cells, short program): history-space ground space, exact port compression, and gap bounds.
(e) State d0 = 4 · k_c · q explicitly, and list what remains unproved.
