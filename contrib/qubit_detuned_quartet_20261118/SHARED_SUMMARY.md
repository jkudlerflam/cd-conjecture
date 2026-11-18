# Handoff: connected p153c quartet and detuning

2026-11-18, OpenAI Codex, GPT-6 Astra (Max). All new work is confined to
this contribution; no shared overview or earlier contribution was edited.

**Verdict: undecided for thermodynamic two-boson stability / Virasoro c=2.**
The new result is a connected four-single-leg finite-density response and
an exact source-projected Cauchy bound, evaluated at fixed p153c bulk field.
The full leading blockade model retains both attractions and both pair
sources: gamma2=-0.00203925056208, gamma3=1/30 in units t=1/80.

The reference is the grand ground state found over every N and momentum.
A two-mode Legendre transform at external kA=-pi+2pi/L, kB=-2pi/L makes
every anomalous Wick, folded-energy and external-leg counterterm vanish
by charge/momentum. Its source-squared coefficient is a sum of 90
interacting-resolvent words, with five denominators and four single legs.

| L | selected N | `g4,L/t` | Cauchy bound |
|---|---|---|---|
|20|3|-0.000309319086|0.001169701615|
|24|5|-0.001233862897|0.002907939935|

Here `g4,L=L c6/(chiA^2 chiB^2)` is the coefficient of the formal
source parameter squared, not the full lambda=1 response or a continuum
cosine. The L=16 independent 2,207-state real-space Rayleigh recursion
agrees with the momentum-resolvent implementation to 8.5e-15 relative.
The off-symmetry B2 interference reduces the two magnitudes relative to
their B3-only values by 1.93% and 0.74%. No free intermediate denominator
is used.

Two explicit obstructions were resolved or exposed:

- Amputating a connected **composite pair** response is insufficient:
  an exactly quadratic squeezed oscillator has a nonzero pair inverse
  response but identically zero four-single-leg vertex. The true
  four-leg calculation is supplied separately.
- A global hard spectral projection removes almost all soft external-leg
  weight and changes the normalized response drastically. Its exact
  finite response is not a local Wilson coefficient. The size dependence,
  normal-order/fusion map, microscopic error, and matching scale remain
  unresolved; no order-one C4 is assumed.

RG.md fixes `y4=4 pi a0^2 |lambda4|/u` and proves a coupled comparison
including bounded dimension drift and OPE mixing. For a maintained
mismatch lower bound, let `l*=log(1/(a0 q_min))` (clipped at zero),
`A=max(0,2-x_min+D rho^2 l*)`. The sufficient total initial budget is
`Y0<=rho/[exp(A l*)+(C rho/A)(exp(A l*)-1)]`, with the continuous A=0
limit. The full asymmetric parity-even operator lattice and generated
near-resonances must be included. Failure is not proof of instability.

For the deliberately hypothetical calibration a0=4, q_min=0.01,
x_min=1.5, C=D=1, rho=0.1, this requires
`|lambda4|/t <= 7.7712e-5 (u/t)` if other vertices use no budget.
The actual amplitude, matching map/error, covariance, velocities, a0,
OPE bounds and positive thermodynamic mismatch lower bound are unknown.
The companion interface's finite-window q4~+0.0166 and bare q4=-0.0108673
are explicitly not interchangeable with a certified cutoff.

The saved C195 `CP` was verified to be a neutral XX+YY bond correlator,
not a creation-pair or dynamic quartet measurement. The dated p153c
protocol in REGISTRATION.md identifies the minimum additional data;
the historical sym40 protocol is unchanged and no production was launched.
The separate chirality-Ising result is a different Hamiltonian and does
not identify the operator content here.

Reproduction: `medium.py` (7.3 s), `check_matching.py` (0.2 s),
`flow.py` and `audit_inputs.py` (under 2 s combined), all bounded and
single-worker. JSON files retain parameters, individual word terms,
source-channel interference, bounds, checks and input hashes.
