This log was written by a Claude Opus 5.5 agent (C196, 2026-11-10), not by Codex.

# C196 progress log

1. Registration I0 to I4 copied into RESULT.md before any code was written.
2. Copied z2dmrg.py, models.py (from work_C190) and z2dmrg_xp.py (from work_C190/gpu), unchanged.
3. idmrg196.py: two-site McCulloch iDMRG on the C190 solver (exact Z2 labels, kron4-preconditioned Davidson, residual 1e-9, at most 300 iterations), with a two-site cell transfer matrix resolved by parity and a thick-restart Arnoldi.
   - validate196.py ed: untruncated iDMRG reproduces open-chain ED (bulk MPO, P = 0) for L = 2 to 10 to 5e-14 at p153c, sym40_n316 and the TFIM. The two-leg XX control matches to L = 8 and is 0.2 too high at L = 10. Cause: zero-weight Schmidt states are dropped, and the separately conserved leg numbers then cannot change (the same trap as C189's amendment A1).
   - validate196.py tm: Arnoldi agrees with dense eig of the cell transfer matrix to 2e-13 in both sectors.
   - From the empty state the XX chain never fills (a localized particle costs energy), so a finite C190 DMRG of 32 sites (random start, mixer) was added as the initialization, plus White's mixer for the first 40 iterations and a 1e-6 random admixture to each Davidson start vector.
   - TFIM h = 1.05, χ = 32: e exact to 1e-13, λ0 = 1 to 1e-15, ξ_odd = 19.2 and ξ_even = 8.7 (exact asymptotic values 20.5 and 10.2). The fixed point is reached in about 60 iterations. The code is correct.
   - XX chain, χ = 32, 400 iterations: no fixed point. e wanders by ±5e-3 around the exact value, the density drifts between 0.184 and 0.199, and the prediction residual stays at 1e-2 to 0.5.
   - p153c, χ = 32, 400 iterations: no fixed point. e drifts upward from −10.00105 to −10.00095 (C190 bulk −10.0011226), the two cell sites develop densities 0.08 and 0.34 that alternate between iterations, and the prediction residual stays near 3e-2.
   - Diagnosis (PROOF 2). In a chain whose particle number is conserved (XX) or conserved up to virtual pairs (the candidate), the growing window must add 2n = 0.39 to 0.43 particles per step at incommensurate filling. Exact parity forces an even number. There is no fixed point, only a quasi-periodic orbit.
4. vumps196.py: one-site VUMPS with exact Z2 labels, MPO environments by GMRES with the fixed-point projector, A_C from C190's Davidson with the one-site analogue of the kron4 preconditioner, and C from Davidson.
   - TFIM, D = 16: ε = 3e-11, e 4e-10 above exact, λ0 = 1.
   - XX chain, D = 16 to 64: ε below 1e-9 at every D, e error 5.4e-5 down to 1.2e-6. The state is a Z2 cat of two U(1)-broken components, so the even and odd transfer spectra coincide and the odd sector has |λ| = 1. c from S against ln ξ_even between consecutive D: 1.015, 1.037, 0.975, 0.985.
   - p153c, D = 32, random start. Plain VUMPS: ε stays at 3e-3 to 1.6e-2, e jumps by up to 5e-3, the C solve hits its 300-iteration cap, and the state sits near the dressed vacuum (n ≈ 0.002). C from the projection A_L† A_C instead: unstable (ε about 1). A two-site centre with C190's TwoSite: unstable (ε about 1). The C solve preconditioned through the A_C preconditioner (H_C = A_L† H_AC A_L): the C solves converge (35 to 117 matvecs), but ε stays near 1e-2. Damping 0.3: ε 3e-4 to 1e-3, e within ±3e-4, still near the vacuum (n = 0.0027).
   - Annealing from a dense liquid (μ_eff = −1.0 down to the target −2.1426) at D = 32: see below.
   - Anneal (local/vanneal.py, D = 32, damping 0.5, 200 iterations per field, μ_eff = −1.0, −1.5, −1.8, −2.0, −2.1, −2.1426): ε = 0.43, 0.86, 3.2e-2, 1.1e-2, 5.2e-3, 9.2e-4, with n = 0.074, 0.041, 0.0038, 0.0023, 0.0027, 0.0028. Even at μ_eff = −1.0 the state does not reach a dense liquid, and no field converges.
5. Stopped development here. No Typhon job was submitted (see RESULT.md, status).
