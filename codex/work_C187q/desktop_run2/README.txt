Discriminating test of the Fibonacci-cusp candidate at low density (coordinator, Claude Opus 5.5, 2026-11-09).
Effective model of C173 (PROOF eq. 3.2) with C176's blockade-exact, N-conserving DMRG (fibdmrg.py, unmodified), at |w| = 1, phi = pi/2, N = L/12.
  fit8*  : couplings fitted from the full chain at V/Delta = 8 (V2 = -1.7334, V3 = -0.4930), where the opposite-valley pair is UNBOUND.
  sw*    : bare Schrieffer-Wolff couplings (V2 = -2, V3 = -2/3), where the pair is BOUND (binding 0.0195, size about 21 sites).
  Both with N - 1 and N + 1 for the parity gap. c1, c2: c = 1 and c = 2 controls at the same N/L.
L = 48, 96, 144, 192. The decisive output is the parity gap: for sw it should level off near 0.02 (pairing, c = 1);
for fit8 it should keep closing like 1/L if the unbound gas is two-component.
Needs python3 with numpy and scipy. Run: bash run_desktop.sh   then  python3 summarize.py
This is the effective model with fitted couplings, not the full qubit chain.
