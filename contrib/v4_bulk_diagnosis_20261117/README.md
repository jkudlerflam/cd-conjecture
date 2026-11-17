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

At **interim** L=144,chi=192, density RMS is about 0.093 and profile
amplitude about 0.129 in both sectors, with q near 0.880 pi. Those states
are still undergoing the registered chi=256 comparison. They are marked
as interim in `analysis.json`; no bulk-phase conclusion is drawn from
this stage. The independent L=96 start, periodic spectra, exact controls,
and secondary XX+YY channel are pending.

`audit_results.py` checks covariance identities and positivity, pair
subtraction/Cauchy bounds and agreement of total energy reconstructed
from density/bond profiles with the independent energy contraction.
All eight saved stage/baseline observations pass at this checkpoint.
