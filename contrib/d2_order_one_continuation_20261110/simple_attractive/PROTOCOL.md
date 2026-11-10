# Small-integer attractive comparison

Registered before new output, 2026-11-10. This is a new rational comparison
to the existing `attr_ref` (phi=1, Delta=.7, V=-2.4), whose old L48→96
entropy slope is about1.13 and archived L96→192 slope about1.31. The
continuation independently checks that192-site reference. This new point
has not previously been tested in this contribution:

`H=sum[7XX-YY+5(XY-YX)-7ZZ]/12`.

Its field is exactly zero, raw coefficient spread is7, there are five Pauli
strings and four grouped terms with three distinct magnitudes. In the
fermion convention, `t=sqrt(136)/12`, `phi=atan(5/3)`,
`Delta/t=8/sqrt(136)`, `V/t=-28/sqrt(136)`. This approximates the reference
without trigonometric couplings or a fine field adjustment. Record the
exchange singular values and cancellation scales as well as raw spread.

Run L48 and96, both parity sectors, seed17, finalchi128, with intermediate
stages16,32,64,128, at most8 clean sweeps per stage and the unchanged local
residual/stopping gates. Use the same two entropy windows, full density,
long correlations and global variance. One sequential worker, isolated
records/checkpoints, no changes to earlier records. Select the lower
variational energy and retain both sector results. Do not promote based on
one parity or on a single high entropy. Only if both window slopes retain
above-one growth beyond numerical uncertainty, register192 and a chi192
comparison; the latter are not yet running or claimed complete.

This is a bounded rational comparison, not a parameter optimization. No
unlisted offsets or extra points are allowed without a new registration.
