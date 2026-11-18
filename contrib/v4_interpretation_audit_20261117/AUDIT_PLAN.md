# Saved-state interpretation audit, 2026-11-17

Author: OpenAI Codex, GPT-6 family. This is a post-hoc audit of known data,
not a blind prediction. The exact d=2 Hamiltonian and full endpoint fields
are fixed. The preceding completed/interrupted archive is read-only.
No DMRG optimization, larger chain, field scan or automatic resume occurs.

## Bounded work

1. Reproduce the saved entropy arithmetic, including both fixed parities,
   the lower-energy branch, W=0,1,6,12, central-quarter/half fractions,
   and averaged-log-chord geometry. Retain complete and interrupted states
   separately. Decompose observed changes; no chi difference is an error bound.
2. Compare density, wavevector and parity-energy/density changes with the
   entropy changes. No three-size fit with an arbitrary correction law.
3. Regress cut entropy on a constant and log chord on the central half
   and middle three quarters. Compare its residual with the cut-crossing
   physical bond, transverse bond, and the adjacent-site averaged density,
   each residualized against the same smooth basis. Retain raw correlations,
   cross-validated left/right-half prediction, and slope sensitivity to one
   nuisance profile at a time. This is a correlation diagnostic, not a
   validated subtraction formula for the interacting model.
4. Apply identical diagnostics to exact fixed-parameter one-/two-leg XX
   controls. The one-leg bond is available; for two legs retain the
   crossing-leg kinetic energy with its explicit support convention.
   Do not subtract their finite-size errors from V=4.
5. Audit the exact parity/translation and valley-exchange symmetry, the
   full leading projected Hamiltonian (blockade, t2, V2, V3, chemical
   shift), and the allowed number-breaking operators. Contrast one- and
   two-component momentum interpretations without treating raw density
   as an exactly conserved quasiparticle density.
6. If needed, read preserved checkpoints for one-body matrices only,
   with no tensor optimization. Use a central-quarter window at every L.
   First benchmark the old L=96 even state and reproduce its archived
   matrices. Proceed only if a conservative chi^3 window^2 projection
   for L=144,chi=320 is <60 seconds and measured RSS <2 GB. At most five
   further reads (48 both parities, 144 both at chi256, interrupted odd
   chi320); one CPU worker, one BLAS thread; total measurement cap 180 s.
   Existing L=96 odd data are reused. Skip on missing checkpoints or cost
   failure. Report Fourier resolution, gauge and matrix normalization.
7. Any small exact algebra/control validation uses the same single worker.
   No additional eigenstate production is authorized. All read-only JSON
   reductions are small and timed. Finish with a concrete missing test,
   its cost from existing sweep timings, and predeclared stop criteria.

## Interpretation constraints

No c assignment from an entropy slope, balanced valley weights, a static
peak, or a finite-size parity gap alone. No effective-model K or velocity
measured at V=40 or at density 1/8--3/16 is transplanted to this V=4 point
near raw density 0.29. The target is compact unitary Virasoro c; c_ent is
only an entropy slope and is not the non-unitary quantity c_eff.

Cross-project recommendations stay in this contribution; AGENTS.md,
the coordinator's folder and top-level status documents are not edited.
