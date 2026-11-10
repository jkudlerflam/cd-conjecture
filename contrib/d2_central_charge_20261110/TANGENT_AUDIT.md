# Independent infinite-environment audit — 2026-11-10 19:18 UTC

The critical Ising one-site chi24 refinement became unstable at its final
iteration and fails several original gates. It remains excluded. Other
accepted Ising and candidate records are unchanged, but this motivates
an independent stationarity check rather than assuming stable sweep
averages exhaust the infinite-environment error.

On immutable accepted checkpoints, solve the infinite MPO environments
with TeNPy's MPOTransferMatrix eigenvector method (tol=1e-12), independently
of the finite-range initialization/environment sweeps used by DMRG.
Compare the two environment energies with the directly contracted energy
and evaluate the left/right one-site tangent projections used by VUMPS.
Check accepted real-frame candidate one-site chi64,96,128 and the two-site
chi192 state, plus accepted critical-control one-site chi16 and two-site
chi24 retry. No checkpoint is modified and no optimization is performed.

For a fixed-chi one-site stationary state, aim for tangent norms below
1e-7. A two-site truncated state may retain a larger tangent residual;
report its value rather than treating a pre-truncation local Ritz test as
a guarantee about the final truncated state. Large residuals would block
a stronger physical interpretation and motivate an exact-environment
variational refinement. This audit has no assumed central-charge outcome.
