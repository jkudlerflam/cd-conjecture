# Higher-chi discrimination — 2026-11-10 18:41 UTC

The real-frame two-site states at chi48,64,96 all pass the original gates.
Their (S,xi) are (1.05449634,84.53352), (1.06902549,100.51955), and
(1.08364724,119.63520). The successive conditional slopes are 0.503307
and 0.503924. The nearly symmetric branch is consistent with a half-unit
slope over this limited range, but xi grows slowly with chi and this does
not yet distinguish a critical mode from a massive crossover.

Extend the completed chi96 state to chi128 with the same Hamiltonian,
solver and acceptance gates. Add one-site chi96 refinement, retaining it
as a separate sequence. Consider chi192 after the completed chi128 result.

For the periodic chi64-size run, take an immutable copy of the current
L64 chi128 even checkpoint to initialize chi192. This is an initialization,
not an accepted observation. The original chi128 worker continues, and
the chi192 worker must independently pass all gates. Both parity energies
will be reported using an explicitly separate energy-only gate if the
odd-state entropy continues to move within its near-degenerate subspace;
such a state will never enter an entropy fit without passing the entropy
gate too.

Extend the independently initialized, original-complex-frame four-site
chi24 state directly to chi96. Agreement at chi24 alone does not test
the later nearly symmetric branch. This run checks whether restricting
the computational cell to two sites missed a lower-energy state.
