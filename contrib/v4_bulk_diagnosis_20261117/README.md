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
