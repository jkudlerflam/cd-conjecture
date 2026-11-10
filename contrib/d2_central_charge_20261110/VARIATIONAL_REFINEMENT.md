# Fixed-chi refinement — 2026-11-10

The complex-frame chi32 two-site infinite run has stable energy and entropy,
but its final sweep norm defect 1.39e-5 exceeds the registered 1e-5 gate.
It is retained and excluded from accepted fits; post-run canonicalization
does not erase this failure. Add one-site DMRG at fixed chi, initialized
from the saved two-site state, with no mixer and the same local-residual,
energy, entropy, and pre-canonicalization sweep-norm acceptance checks.
This removes the recurring two-site truncation from the variational
optimization and tests whether a stationary state exists at the same chi.
One-site and two-site sequences are reported separately. A new convergence
pass must be earned; the threshold is not relaxed.

## Periodic extension — 2026-11-10 18:44 UTC

The L48 chi128 ring has stable energy/entropy but its final two-site norm
defect is 6.62e-5, above the same 1e-5 gate. It is excluded. Apply the same
fixed-chi one-site procedure to this completed state and, if needed, to the
completed L48/L64 chi192 states. Wait for the source worker's final record;
do not overwrite or weaken any two-site result. New one-site runs have
4 minimum and 12 maximum sweeps, the same independently checked local
residuals and all original acceptance gates. Failed refinements stay failed.

At 19:07 UTC, also queue the same one-site procedure for completed L48
odd chi128 and both L64 chi128 states if their original gates fail. This
lets the chi128→192 entropy/energy changes be assessed using stationary
states at both dimensions, without relaxing the recurring two-site norm
failure. No extra chain size or coupling is introduced.
