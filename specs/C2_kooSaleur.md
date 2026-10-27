# C2: Does a ground-state-only Koo-Saleur ratio converge to c for non-integrable critical chains?

Closed book: do not read any files you did not create. Work only inside the current working directory;
write and run Python (numpy, scipy) there. Keep total runtime under about 15 minutes.

## Background
Periodic chain of N sites, H = sum_j h_j with h_j a local term (two or three sites). Unscaled Fourier modes
h_n = sum_j exp(2 pi i n j / N) h_j. Ground state |I>, energy E_0. Define
  m_0(N) = <I| h_2 h_{-2} |I>,        m_1(N) = <I| h_2 [H, h_{-2}] |I> = <I| h_2 (H - E_0) h_{-2} |I>,
  c_2(N) = 8 m_0^3 / m_1^2.
(Here h_{-2} = h_2^dagger. c_2 is invariant under H -> lambda H and under constant shifts of h_j.)
In a CFT, h_n -> (2 pi v/N)(L_n + Lbar_{-n}) on CFT states and H - E_0 -> (2 pi v / N)(L_0 + Lbar_0), which
gives m_0 = (c/2)(2 pi v/N)^2, m_1 = c (2 pi v/N)^3 and c_2 = c, PROVIDED only CFT (low-energy) intermediate
states contribute. This formula was proposed in Latorre-Sierra arXiv:2403.17242v1 (appendix) as a route to the
central charge from the ground state alone.

## Tasks
A. Analytic. Estimate the contribution of high-energy intermediate states (energies O(1) in lattice units) to
   m_0 and m_1 as N -> infinity. Use [H, h_q] = O(q) x (Fourier mode of the lattice energy current) and
   express the leading high-energy contribution through the ground-state fluctuations of the total energy
   current J_tot (defined by continuity, d h_j/dt = -(J_j - J_{j-1})). Show that generically m_0 and m_1 pick up
   terms of order 1/N (compared with the CFT terms of order 1/N^2 and 1/N^3), so that c_2(N) -> 0 like 1/N,
   unless the high-energy weight of J_tot|I> vanishes, as happens when J_tot commutes with H (XXZ chain,
   transverse-field Ising chain). Be explicit about powers of N and the coefficient.
B. Numerics (sparse exact diagonalization, periodic boundary conditions, N as large as feasible, at least
   N = 8..20 in steps of 2 for spin-1/2; use momentum-space tricks only if needed). For each model report
   c_2(N), m_0 N^2, m_1 N^3, and the ground-state variance <J_tot^2> - <J_tot>^2 divided by N.
   Also report the Milsted-Vidal low-energy estimate c_MV(N) = 2 |<T| h_{-2} |I>|^2 / (E_T - E_0)^2 * (2 pi/N)^2 *
   4 ... -- derive the correct normalization yourself so that c_MV -> c in a CFT, where |T> is the eigenstate in
   the momentum sector 2 x (2 pi/N) with the largest overlap with h_{-2}|I> among the lowest ~10 states there.
   Models:
   1. Critical transverse-field Ising H = -sum_j (X_j X_{j+1} + Z_j)  (integrable, c = 1/2).
   2. Heisenberg H = sum_j S_j . S_{j+1}  (integrable, c = 1).
   3. Non-integrable Ising-critical chain (O'Brien-Fendley): H = sum_j [ -X_j X_{j+1} - Z_j
      + lam (X_{j-1} X_j Z_{j+1} + Z_{j-1} X_j X_{j+1}) ], Kramers-Wannier self-dual and Ising-critical
      (c = 1/2) for 0 <= lam < about 0.42. Use lam = 0.2 and lam = 0.35. Take h_j centred on site j.
   4. Non-integrable c = 1 chain: H = sum_j [ S_j . S_{j+1} + J2 S_j . S_{j+2} ] with J2 = 0.2
      (critical, c = 1 for J2 < 0.2411).
   Check translation invariance of the computed ground state and that it is unique (for N where the
   ground state is degenerate, say so and skip).
C. Conclusion: does c_2(N) converge to c for the non-integrable models? Does c_MV(N)? Fit the N-dependence.

Deliverable: save the code as ks_test.py in the working directory and give the tables and fits in your answer,
together with the analytic estimate of part A.
