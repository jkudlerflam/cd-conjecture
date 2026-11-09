Quick test of the Fibonacci-cusp candidate (coordinator, Claude Opus 5.5, 2026-11-09).
Effective model of C173 (PROOF eq. 3.2) with C176's blockade-exact, N-conserving DMRG (fibdmrg.py, unmodified).
  fit8*  : couplings fitted from the full chain at V/Delta = 8 (V2 = -1.7334, V3 = -0.4930), |w| = 1, phi = pi/2, n = 1/6,
           with N - 1 and N + 1 for the parity gap. The pair is unbound at these couplings.
  sw     : bare Schrieffer-Wolff couplings (V2 = -2, V3 = -2/3), where the pair is bound.
  c1, c2 : known c = 1 (hard rods) and c = 2 (two decoupled XX legs) controls at the same sizes.
L = 48, 96, 144, bond dimension up to 384 (512 at L = 144).
Needs python3 with numpy and scipy. Run:  bash run_desktop.sh      then  python3 summarize.py
This is the effective model with fitted couplings, not the full qubit chain.
