C7 session memory, October 27, 2026.

Scope. Closed book, confined to work_C7. No pre-existing project files were read or listed. The root project STATE.md was not accessed or changed because the user explicitly prohibited this. This file records the session locally.

Files created.

- C7_result.md is the theorem-style deliverable, including all assumptions, proofs, counterexample attempts, and limitations.
- c7.py contains small numpy/scipy consistency checks. It reads no input files. Run with python3 c7.py. The shell has no python alias.
- C7_STATE.md is this memory.

Verdict. No unconditional f(d) bound and no fixed-d conformal counterexample were established. Do not cite the long-period model as a CFT or as evidence that a conformal family with unbounded central charge has been constructed.

Proved conditional bounds.

1. If the full infrared stress tensor is Sugawara for independent microscopic on-site compact symmetry currents, c_T <= dim(g) <= d^2 - 1, with no affine-level bound. The uncontrolled remainder without completeness is c_coset.
2. For a unitary compact modular-invariant CFT with c >= 1, the number of states below Delta = c/4 + 2 is at least 1 + (pi e/12)c exp(pi c/6). A modular linear functional proves this without a uniform Cardy assumption. If these states are represented on a ring L <= K(d), then c_T <= 1 + (6/pi)K(d) log d. The bounded-ring realization is an extra hypothesis, not a result.

Explicit partial construction.

Alphabet #, L, X, R, d = 4. Allowed pairs are #L, #X, LL, LX, XR, X#, RR, R#. Let Q_i project onto forbidden pairs and M_i onto #. The Hamiltonian is K sum Q_i + mu sum M_i - t sum (|LX><XR| + h.c.), with K > 4t and mu = 2t - delta. Every Q_i and M_i is conserved.

A cell from a marker to the next has length b including the first marker, b - 1 head positions, and lowest energy a_b - delta, a_b = 2t(1 - cos(pi/b)). Define delta_b = (b+1)a_b - b a_(b+1). The delta_b decrease. Choosing delta_b < delta < delta_(b-1) selects the unique energy-density minimum at b. For every b >= 5 this choice has 0 < delta < 2t.

For rings divisible by b there are b ground states, one for each translated marker pattern. Illegal sectors have energy >= e_b L + (K - 4t)q, where q is the number of forbidden bonds. The gap above the ground space is positive uniformly in compatible ring length for each fixed b. The interaction norm is bounded by K + 3t independently of b. At delta = 0 the gap instead behaves as 2t(1 - cos(pi/L)), of order L^-2, so the accumulation point is nonconformal.

Adding a free qubit per site gives d = 8 and N = 2^b ground states per block. This has extensive zero-energy degeneracy and is not a CFT. Producing the ULS_N exchange with an L-independent, locally controlled, symmetry-protected effective Hamiltonian is unresolved.

Boundary obstruction. Direct compression of one microscopic boundary bond to product block codes has connected operator Schmidt rank <= d^2 - 1. The connected swap P_N - I/N has rank N^2 - 1. Thus a nonzero ULS_N bond cannot arise by direct projection or first-order perturbation theory for N > d. This does not bound arbitrary higher-order resolvents or dressed encodings.

Remaining construction issues. Uniform-in-volume local perturbation control, critical stability including marginal operators, encoded symmetry and correct translation action, incompatible periodic ring lengths, and removal of extra gapless clocks or defects. Finite-size Hamiltonian universality alone settles none of these.

Checks. All final c7.py checks pass in under one second. They cover local commutators, period selection, small rings, modular-functional signs, and connected swap rank. An initial unshifted sparse eigensolver omitted isolated zero-energy product states. The final code applies a scalar shift and explicitly checks those states. With t = 1, K = 8, midpoint delta: L = b = 5 has E0 = -0.794593112989, degeneracy 5, next gap 0.794593112989; L = b = 6 has E0 = -0.551711241165, degeneracy 6, next gap 0.551711241165. These are consistency checks, not the proofs.
