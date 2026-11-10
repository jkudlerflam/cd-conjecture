# Bulk central-charge follow-up, 2026-11-10

OpenAI Codex, GPT-6 family. Registered before new numerical results.
Input: local commit 7014ef41 and its fixed candidate
`(Jx,Jy,Jz,Dz,hz)=(1/2,-1/2,5/8,1/10,-7/8)` in Pauli units.
Previous raw results and states are immutable. No field retuning.

The question is whether the open-chain effective entropy coefficient near
1.45 is a bulk central charge. First discriminants:

1. Periodic full-spin DMRG using independent TeNPy 1.1.1, with a finite
   open MPS representation of the ring Hamiltonian. Validate energies and
   entropies against independent full-Hilbert-space ED at L8 in both
   parities, and at L16 against saved ED if available. Then L32,48,64 in
   both parities, initially chi128, with chi192/256 checks on promising
   completed states. Increase sizes only after checking accuracy. The
   periodic estimator is `3 Delta S/log(L2/L1)` and profile slope times 3,
   using the central three and fractional [3L/8,5L/8] cuts as before.
   Fixed-parity comparisons and both-parity availability remain explicit.
2. Infinite-chain DMRG with 2-site computational cell, physical one-site
   uniform Hamiltonian, initially chi16,24,32,48,64. Compare independent
   initializations and 4-site cells before promoting a finite-entanglement
   fit `S=(c/6) log(xi)+const`. Cell size does not change the Hamiltonian.
   Contracted finite-range environments bypass the previously documented
   imaginary-scalar iterative-environment initialization issue. Require
   independently contracted local energy, returned/local agreement, stable
   late energy/entropy, small pre-canonical norm defect, and checked local
   solver residuals. Preserve failed runs; do not fit unconverged points.
3. Controls: exact small rings and XX finite-entanglement scaling. Calibrate
   against c=1 without using the calibration to force a target result.
   Additional sizes, bond dimensions, or diagnostics will be dated adaptive
   amendments, not retroactively treated as initial registration.

Do not infer c from an inverse-size energy term without a velocity, or from
a single entropy fit. Ground-state selection remains variational. Several
velocities, metastability, symmetry breaking, oscillations, and nonrelativistic
critical behavior may invalidate a naive single-scale interpretation.

Primary references for the conditional estimators:
- Calabrese and Cardy, https://arxiv.org/abs/hep-th/0405152
- Calabrese and Cardy review, https://arxiv.org/abs/0905.4013
- Pollmann et al., https://arxiv.org/abs/0812.2903

Every process uses one BLAS thread. Checkpoints are local and Git-ignored;
raw JSON records, source hashes, logs, and explicit status are retained.
No cluster or remote publication is authorized by this protocol.
