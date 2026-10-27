# G1: Exchanging large logical spins through one small bond

Closed book: do not read or list any files you did not create. Work only inside the current working directory.
Python with numpy and scipy is allowed. Keep the total runtime of your scripts under 40 minutes. Be honest: a clean
negative result or a precisely located obstruction is as valuable as a construction. Flag every step you assume
rather than prove.

## Setting

Two blocks A and B. Block A has Hilbert space C^d (the "port") tensored with a finite-dimensional space K_A of any
size; block B likewise, C^d ⊗ K_B. The unperturbed Hamiltonian is H0 = H_A ⊗ 1 + 1 ⊗ H_B with Hermitian H_A ≥ 0,
H_B ≥ 0. Each has an exactly degenerate ground space at energy 0, G_A and G_B, of dimension N, with all other
eigenvalues ≥ Δ > 0. Fix isometries identifying G_A ≅ C^N ≅ G_B. The perturbation is V = λ h ⊗ 1_{K_A} ⊗ 1_{K_B},
where h is a Hermitian operator on port_A ⊗ port_B = C^d ⊗ C^d only. Write Π for the projector on G_A ⊗ G_B and
Q = 1 − Π. The effective Hamiltonian on G_A ⊗ G_B (Schrieffer–Wolff / Bloch, any standard convention) is
H_eff = Σ_n λ^n H^(n). For example H^(1) = Π h Π and H^(2) = −Π h Q H0^{-1} Q h Π.

Call an operator on C^N ⊗ C^N "connected" after subtracting its one-block parts (X ⊗ 1 + 1 ⊗ Y). Its connected rank is
the rank of its coefficient matrix in traceless operator bases. The swap P_N on C^N ⊗ C^N has connected rank N² − 1.
Known fact (you may use it): the connected rank of H^(1) is at most d² − 1, because it is a compression of h. So
H^(1) cannot contain a nonzero multiple of P_N when N > d.

## Questions

1. Construction. For N > d, find explicit (H_A, H_B, h) such that the first nontrivial connected term of H_eff is
   proportional to the swap. Precisely: for some order k, every H^(n) with n < k is of the form a_n·1 (one-block terms
   proportional to the identity, no connected part), and H^(k) = a_k·1 + b_k·P_N with b_k ≠ 0. Start with d = 2 and
   N = 3 or 4, then d = 3 and N = 4 or larger. Report the smallest k you achieve and the sign of b_k (b_k > 0 is the
   antiferromagnetic sign we prefer, but report either). Verify numerically: diagonalize H0 + λV exactly for several
   small λ, take the N² lowest levels, and check that they split into the symmetric and antisymmetric subspaces of
   C^N ⊗ C^N (degeneracies N(N+1)/2 and N(N−1)/2), with splitting ∝ λ^k. Also check that the effective Hamiltonian
   extracted from the exact low-energy eigenvectors (for example the des Cloizeaux effective Hamiltonian) agrees with
   a·1 + b·P_N to the expected order.
2. Second order. Is a nonzero swap coupling possible at k = 2 when N > d? If yes, how large can N be for given d (give
   the best bound you can prove, and an example near it). If no, prove it.
3. Lower bound on the order. Is there a bound of the form k ≥ g(N, d) (for example k ≳ 2 log_d N, an information
   counting through the port) for the first order at which a swap coupling can appear with all lower orders trivial?
   Prove it or explain where a proof fails.
4. Locality. Your H_A may be any Hermitian operator. Say whether the construction could plausibly be realized with H_A
   a nearest-neighbour Hamiltonian on a chain of qudits of fixed dimension (the port being the last site), and what
   would grow with N (block length, energy ratios, k).

## Deliverable

A report gadget.md with the constructions or proofs, theorem-style where possible, and code gadget.py that reproduces
every numerical check (numpy and scipy only).
