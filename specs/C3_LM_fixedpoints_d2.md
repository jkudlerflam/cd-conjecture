# C3: Exact single-site Lin–McGreevy vector fixed points for qubit chains (d = 2)

Closed book: do not read or list any files you did not create. Work only inside the current working
directory. You may write and run Python (numpy, scipy). Keep total runtime under about 30 minutes.

## Background and definitions
Infinite chain Z of qubits (local dimension d = 2). Let omega be a translation-invariant (TI) state, with
marginals rho_X on finite sets X of consecutive sites. Write tau := rho_{site} (2x2), sigma := rho_{2 sites}
(4x4, with Tr_1 sigma = Tr_2 sigma = tau), rho_3 := rho_{3 sites} (8x8) with Tr_3 rho_3 = Tr_1 rho_3 = sigma.
Entanglement Hamiltonians K_X := -log rho_X (assume full rank so they are finite). Entropies S_X.

For three consecutive intervals A, B, C with lengths a, b, c, the vacuum cross ratio is
  eta = a c / ((a+b)(b+c)),
and the Lin–McGreevy (arXiv:2303.05444) combination is
  K_Delta = (K_AB + K_BC) - eta (K_A + K_C) - (1 - eta)(K_B + K_ABC),     (each K_X tensored with identity)
  S_Delta = <K_Delta> = eta [S(B|A) + S(B|C)] + (1 - eta) I(A:C|B) >= 0.
A 1+1D CFT vacuum satisfies the vector equation K_Delta |psi> = (c/3) h(eta) |psi> for every consecutive
triple, with h(x) = -x log x - (1-x) log(1-x) and one universal c. For a purification |psi> of omega,
the vector equation for the triple ABC is equivalent to (K_Delta - lambda) rho_ABC^{1/2} = 0; when rho_ABC
has full rank this is the operator identity
  K_Delta = lambda * 1   on ABC.                                                    (*)

Single-site triple: a = b = c = 1, eta = 1/4, h(1/4) = 0.5623351446. Then (*) reads
  (3/4) log rho_3 = log sigma_12 + log sigma_23 - (1/4)(log tau_1 + log tau_3) - (3/4) log tau_2 + lambda,
so rho_3 is determined by sigma, and (*) becomes a fixed-point condition on sigma:
  rho_3 = exp( (4/3)[ log sigma_12 + log sigma_23 - (1/4) log tau_1 - (1/4) log tau_3 - (3/4) log tau_2 ] + (4/3) lambda ),
  with Tr rho_3 = 1 (fixes lambda), and Tr_3 rho_3 = sigma, Tr_1 rho_3 = sigma.                (**)
Define c_Delta := 3 lambda / h(1/4) (= 3 S_Delta / h(1/4)).

Known example (check it): the tracial state (all marginals maximally mixed) solves (**) with
lambda = (1/2) log 2, so c_Delta = 1.849 > d - 1 = 1. More generally for the tracial state
K_Delta = 2 b eta log 2 for every triple (a, b, c), so its "c" is triple dependent, i.e. it is not scale
consistent.

## Tasks
1. Classify all solutions sigma of (**) with full-rank marginals (up to the local unitary symmetry
   sigma -> (U x U) sigma (U x U)^dagger). Do both:
   (a) analytically, the commuting/classical case (sigma diagonal in a product basis: a stationary
       two-site probability p(x,y) with equal marginals q), including all solutions, not only products;
   (b) the general quantum case, numerically: solve (**) by damped fixed-point iteration and by
       root finding from many random full-rank starts (and starts near pure, near Bell-diagonal,
       near Werner states, near known critical-chain marginals such as the critical transverse-field
       Ising and Heisenberg 2-site reduced states). Report the solution set (isolated points or
       families, their dimension), and lambda and c_Delta on each. Is c_Delta bounded on the solution
       set? What is its supremum? Is any solution with 0 < c_Delta <= 1 found?
2. Scale consistency. A genuine lattice CFT state must satisfy (*) for larger triples with the SAME c.
   For each solution of task 1, ask whether it extends to a TI 4-site state rho_4 (marginals rho_3 on
   sites 123 and 234) satisfying (*) for the triples (1,1,2) and (2,1,1) (eta = 1/3) and (1,2,1)
   (eta = 1/9), with eigenvalues (c/3) h(1/3) and (c/3) h(1/9) for the same c = c_Delta. Note that (*)
   for (1,1,2) determines rho_4 from rho_3, sigma and tau in the same way as (**). Find all solutions of
   the combined system (1,1,1)+(1,1,2)+(2,1,1)+(1,2,1), analytically in the classical case and
   numerically in general. Do any nontrivial (c > 0) scale-consistent solutions exist? If yes, give c.
3. Purity. LM's states are pure (ground states on a ring or line). Say whether the 3-site marginals
   you found can arise from a pure TI state on Z (or on a ring of L qubits for some L), at least in
   simple cases; the tracial example in particular.
4. Also report, for comparison, S_Delta(1,1,1) and c_Delta for the exact 3-site marginals of the critical
   transverse-field Ising chain and the Heisenberg chain on the infinite line (free-fermion formulas or
   DMRG/ED extrapolation, whichever you can do reliably), and how far they are from satisfying (*).

## Deliverable
Theorem-style statements for whatever you can prove (especially the classical case), a clear account of the
numerics (method, number of starts, convergence tolerances, residuals of (**)), the list of solutions with
lambda and c_Delta, and a direct answer: do nontrivial exact single-site fixed points exist for d = 2, which
are scale consistent, and what are their c? Save code as c3.py in the working directory. Flag anything
assumed rather than proved.
