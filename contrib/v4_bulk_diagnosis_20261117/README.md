# V=4 bulk-phase diagnosis

**In progress; central charge UNKNOWN.** This bounded round tests the
large density modulation of the fixed nearest-neighbour qubit Hamiltonian
`sum[(9/16)XY+(7/16)YX+ZZ]-(7/4)sum Z`, with full endpoint fields.
No entropy estimate near 2.5--3.1 is treated as a measurement of c.

Author: OpenAI Codex, GPT-6 family, 2026-11-17. Numerical work reuses the
unchanged C190 parity solver (Claude Opus 5.5 agents). New density, bond
and pair contractions are independently checked against literal L=8 ED;
there is no independent human review or second DMRG implementation.

[Registration](REGISTRATION.md) fixes the size/accuracy sequence,
observable definitions, controls, conditional extensions and resource
cap. [Sequential decisions](DECISIONS.md) preserve adaptive choices.
Saved L=96 states are reused; new production is restricted to two workers
with one BLAS thread each. Outputs are JSON; state checkpoints stay local
and ignored. The final recommendation will be continue, reject, or a
specific missing test, not an unsupported central-charge assignment.

## Checkpoint: completed 48/96 comparisons and interim 144

Both L=48 sectors pass the existing accuracy gates and the new observable
stability checks. The reused L=96 states also retain their prior local
gate passes. Their central-half density metrics are:

| L | Parity | Density RMS | Profile amplitude | q/pi | Connected peak / M |
|---:|---|---:|---:|---:|---:|
| 48 | even | 0.11211 | 0.14892 | 0.88770 | 0.01678 |
| 48 | odd | 0.11982 | 0.15932 | 0.88721 | 0.01786 |
| 96 | even | 0.10851 | 0.15767 | 0.88428 | 0.01111 |
| 96 | odd | 0.10588 | 0.15385 | 0.88428 | 0.01084 |

Connected density RMS at separations M/3 to M/2 falls from about
0.0115--0.0124 to 0.0064--0.0066. This is distinct from the large
one-point modulation. The wavevector's fine grid is not evidence for
commensurability at that precision; the physical resolution is limited
by the finite observation window.

Both L=144,chi=256 runs are now complete. Density RMS is about 0.0923
and profile amplitude about 0.128 in both sectors, with q near 0.880 pi.
The entropy chi-step gate fails (0.01894/0.01840), while the other local
gates pass. The registered chi=320 continuations pass the resource
benchmark and are pending/in progress, as is the independent L=96 start.
The current boundary-probe amplitude trigger is not met.

At the user's subsequent request, a **post-hoc** central-charge estimator
has been added in `entropy_analysis.py`. The lower-variational-energy
branch gives c_ent=2.818 for L=48 to 96 and 2.454 for L=96 to 144
using the established W=6 average and 6 delta S/log(size ratio).
The three-size fit gives 2.698; its geometry-corrected version gives
2.631. These are finite-size entropy slopes, not an identified
thermodynamic c or the non-unitary quantity c_eff. The large-size chi
gate failure and size/window/sector dependence remain explicit.

The same estimator on fixed-parameter exact free controls gives
0.995 and 1.954 for known c=1 and c=2 at L=96 to 144. The new valley
channel XX+YY has decreasing uniform normalized weight, about
0.064--0.065, 0.026, and 0.0157 at L=48,96,144. Small periodic-ring
spectra and correlations are saved separately; their entropies have
large size/sector effects and do not determine c.

`audit_results.py` checks covariance identities and positivity, pair
subtraction/Cauchy bounds and agreement of total energy reconstructed
from density/bond profiles with the independent energy contraction.
The saved stage/baseline observations pass at this checkpoint. Final
counts and all completed/running tags are recorded in `audit_results.json`.
