# Comparable parity precision — 2026-11-10 18:47 UTC

The complete ED scan through L24 finds the even k=0 ground state and an
odd k=0/pi doublet. The minimum odd-even gaps are 0.02130686147,
0.01722003298, and 0.01449702172 at L16,20,24; low-level solves within
the ground momentum/parity block find no smaller same-block excitation.
No gap extrapolation is registered as a proof of criticality or a CFT.

To compare larger-ring parity energies at matched chi, initialize real
chi192 odd L48/L64 workers from immutable snapshots of their current
chi128 states. Run four two-site sweeps to expand the basis and then,
if any gate fails, up to twelve one-site sweeps. Every accepted energy
must pass the explicitly separate energy gate defined in analyze.py;
entropy still requires all gates. Retain all chi128 and chi192 results
and disclose their differences rather than treating a fitted intercept
as a thermodynamic gap. No coupling changes or cherry-picked sizes.

The old original-frame L32 odd chi128 worker uses the pre-fix serializer.
After it finishes, explicitly record its expected postprocessing failure,
verify its saved final checkpoint, and initialize an odd chi192 stage
(four two-site sweeps, then one-site refinement if required). This gives
a fresh accepted record and matches the even L32 chi192 precision.
