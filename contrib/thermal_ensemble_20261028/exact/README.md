# Interacting all-length ensemble transfer

Status: **proved abstract reductions; no solution of the c-d conjecture**.

[LENGTH_FUGACITY.md](LENGTH_FUGACITY.md) supplies three connected results.

1. The physical pressure is the unique interior zero of a rooted cell
   pressure with physical-length fugacity, conditional on the exact
   physical/rooted trace comparison. Interacting length-dependent bonds
   and arbitrarily long cells are allowed.
2. C46's massive-reference thermal comparison survives this fugacity
   with the same interaction bound. The unbounded length tilt stays
   inside the on-cell heat trace. The reference's physical pressure
   correction is exponentially small at low temperature when its
   length-tilted heat trace is uniformly exponentially small.
3. The selected b,b+1 renewal ensemble has an actual fixed-L pressure
   limit, proved by orthogonal prefix compressions. A uniform T^2
   expansion of the rooted scalar pressure passes to its zero using
   only one ground-pressure derivative.

These results address the ensemble obstacle left in C46. They do not
control leakage or prove the retained interacting pressure is Gaussian.
[COMPLETION_AUDIT.md](COMPLETION_AUDIT.md) independently accepts the
sibling's positive modification that retains malformed physical sectors
with a complete gap and linear heat-trace cost. The full physical pin
theorem belongs to the sibling [ensemble folder](../ensemble/README.md).

[SEPARATE_JW.md](SEPARATE_JW.md) gives separate charge gradings for the
payload and length modes, including arbitrary charged auxiliary states.
It preserves neutral interval supports and operator norms and gives an
isometric extension whose retained reference is quadratic in n+1 CAR
modes. The auxiliary parity factors outside the vacuum are essential.
`jw_checks.py` verifies all algebra identities on a 256-dimensional
two-cell control, with zero observed matrix error.

`checks.py` uses NumPy only. It tests 18 noncommuting length-tilted thermal
inequalities and four all-length renewal controls. The random interactions
preserve total length and the local all-code block, while moving defects
and changing their number. All checks passed. Results are in
[checks.json](checks.json).

The system Python lacks NumPy. In this workspace the successful command was

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  /Users/pace/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 \
  -B contrib/thermal_ensemble_20261028/exact/checks.py
```

Provenance: OpenAI GPT-6-family Codex agent `exact_escape`, 2026-10-28.
The coordinator checked the pressure conversion and requested the
orthogonal-prefix and one-derivative refinements. The sibling agent's
pinning proof was independently rederived here and its 13 matrix cases
and 39 thermal checks rerun. No human review is claimed.
