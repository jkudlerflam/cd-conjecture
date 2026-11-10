# Exact-environment variational refinement — 2026-11-10 19:20 UTC

The independent environment audit agrees with the recorded energies within
2.2e-10. Candidate one-site tangent norms are 8.8e-7, 9.6e-7 and 1.76e-7
at chi64,96,128, above the new 1e-7 diagnostic target. The critical Ising
one-site chi16 reference has tangent norm 1.0e-9, while a truncated
two-site state has the expected larger residual (4.4e-5 at Ising chi24).

Refine the accepted candidate one-site chi64,96,128 states and two-site
chi192 and chi256 states with TeNPy single-site VUMPS at the same chi.
This solves infinite environments anew at every update. Also check the
accepted critical Ising chi24 two-site retry. These are separate
variational records; do not merge them silently with either DMRG family.
Use at most 40 sweeps, at least 6, no mixer, and the same local Krylov
settings (limit500). Check norm and split errors, four-record energy and
entropy spans, independently contracted energy, and the final exact
environment tangent norm (<1e-7). TeNPy marks its UniformMPS support
experimental; retain failures and require the independent checks.

This is a stationarity audit at already requested bond dimensions, not a
new coupling search or an assumed gapped result. A material lowering of
energy or change of entropy would change the physical assessment.

At 19:23 UTC, the first candidate chi64 check misses only the tangent gate
(1.0334e-7 versus 1e-7); its entropy changes by 1.03e-6 and its energy by
5.4e-12 relative to the source. Keep that record failed. Permit a separately
labelled warm VUMPS retry for a completed check which passes every other
gate and has tangent norm below 1e-6. The same 1e-7 gate remains mandatory.

At 19:42 UTC, the original chi256 check has finished with tangent norm
4.5422e-7 and all other gates passing. Its entropy differs from its DMRG
source by 5.2233e-6. The qualifying warm retry started at 19:41 UTC.
For clarity, the continuation budget is one warm retry per original
state, with the same 40-sweep limit and unchanged gates. A failed warm
retry terminates that stationarity check; it does not trigger an
unbounded sequence of retries or a relaxed threshold.
