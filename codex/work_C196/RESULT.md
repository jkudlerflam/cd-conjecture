This run was done by a Claude Opus 5.5 agent (C196, 2026-11-10), not by Codex. Astra was not used (out of credits).

# C196 result: entanglement scaling on the infinite qubit chain

Spec: `specs/C195_C196_third_size_and_infinite_chain.md`. Derivations and the method are in PROOF.md, the log in PROGRESS.md.

## Registration (copied verbatim from the spec before any code was written)

- (I0) Controls with the same code, χ sequence and fit window: the XX chain reads c within 0.1 of 1, and the two-leg XX control reads c within 0.15 of 2. A gapped control (any bond with a known finite ξ) shows ξ saturating.
- (I1) Convergence, per point: energy density change over the last iterations below 1e-10 per site and gradient or residual norm below 1e-7. Points that miss this are flagged and excluded from the fit.
- (I2) Pass. At each candidate bond the slope of S against ln ξ over the last three converged χ gives c ≥ 1.7, ξ grows by at least a factor 1.5 over the last χ doubling in both parity sectors, and I0 holds.
- (I3) Fail. c ≤ 1.3 from converged points at a candidate bond, or the leading correlation length of one sector grows by less than 1.15 over the last χ doubling while the other keeps growing.
- (I4) Consistency. The bulk energy density agrees with (E(256) − E(128))/128 of the C190 ground branch to 1e-5, and the density with the C190 mid-chain value to 0.006. A miss is reported, not fatal.
- Anything between I2 and I3 is undecided and reported with the numbers.

## Status (2026-11-10, written after the local tests; no Typhon production exists)

Outcome: not evaluated. Neither infinite-chain method I built converges at the candidate field p153c, even at D = 32. So no Typhon controls or production were submitted, and none of I0 to I4 is evaluated. Following the spec, I report the failure with its numbers rather than relax I1. Details are in PROGRESS.md and PROOF.md.

What works [numerical, local, numpy].
- iDMRG on the C190 solver (idmrg196.py). Untruncated, it matches open-chain ED for L ≤ 10 to 5e-14. The TFIM (h = 1.05, χ = 32) reaches the exact energy density to 1e-13, with λ₀ = 1 to 1e-15 and ξ_odd = 19.2 (exact asymptotic value 20.5).
- VUMPS with the one-site version of C190's preconditioner (vumps196.py). TFIM at D = 16 reaches ε = 3e-11. On the XX chain at D = 16 to 64, ε ≤ 1e-9 and the energy error is 5e-5 to 1e-6. c from S against ln ξ_even between consecutive D reads 1.015, 1.037, 0.975, 0.985. The finite-D state is a Z2 cat of two U(1)-broken components, so its two parity sectors carry the same transfer spectrum.

What fails at p153c (h_z = −19.97589125, C190 bulk e = −10.0011226, n_mid = 0.215).
- iDMRG, χ = 32, 400 steps. No fixed point. e drifts from −10.00105 to −10.00095, the cell sites alternate between densities 0.08 and 0.34, and the prediction residual stays near 3e-2. The XX control shows the same thing (e wanders by ±5e-3). Cause (PROOF 2): at incommensurate filling, in a chain that conserves N up to virtual pairs, each step must add 0.43 particles while parity forces even increments.
- VUMPS, D = 32. Plain: ε = 3e-3 to 1.6e-2. With the C solve preconditioned through the A_C preconditioner, every local solve converges but ε stays near 1e-2. Damped (0.3): ε = 3e-4 to 1e-3. Annealed from μ_eff = −1.0: ε = 9e-4 at the target. In every variant the state stays near the dressed vacuum (n = 0.002 to 0.003, against 0.215). Projected-C and two-site-centre variants are unstable (ε ≈ 1). This is C177's failure, now with the local eigenproblems converged, so the local eigensolver is not the cause.

What a working attempt would need. A cure at the level of the outer iteration rather than the local solver. Candidates are (i) an exact U(1)-like label for the real-particle number, after a unitary that removes the virtual pairs to the order needed (C189's labels then apply at a commensurate density), (ii) a gradient-based uniform-MPS optimizer with a line search started from a uniformized C190 bulk state, or (iii) VUMPS in a q-site cell started from a converged C190 open chain.
