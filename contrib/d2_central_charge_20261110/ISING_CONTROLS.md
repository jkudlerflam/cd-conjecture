# Critical versus massive calibration — 2026-11-10 19:02 UTC

After the candidate's chi96→128 ordinary-entropy coefficient decreased
to 0.308 (and the higher Rényi orders to about 0.22), add two fixed controls
to distinguish the numerical behavior of a true half-unit critical theory
from a small-gap crossover. No candidate parameter is changed.

Use H=-sum XX-h sum Z, with h=1 and h=1.01. The Jordan–Wigner dispersion
is epsilon(k)=2 sqrt(1+h^2-2h cos k), so the first control has the Ising
critical point (c=1/2), while the second has exact fermion gap 0.02.
The exact energy density is minus the integral of the square root over
[0,pi], divided by pi; at h=1 it is -4/pi. Compare final energies with
this independent integral. This is also the Ising critical benchmark
in the finite-entanglement literature (Pollmann et al., arXiv:0812.2903).

Run both controls in the original arithmetic implementation, two-site
cell, independently initialized seed23, chi16,24,32,48, up to 120 sweeps
per chi. Keep the same convergence gates. Exclude any failed state rather
than using the expected control result to excuse it. An entropy/xi ratio
formed from nearly saturated differences can be unstable even for a known
massive control, so report entropy and xi themselves as well as slopes.

At 19:06 UTC, critical chi24 and massive chi32 have finished with stable
energies but fail the unchanged four-record entropy gate. Retain these
failures. Permit one separately labelled warm retry of each failed control
point, up to 400 sweeps, with all other settings and gates unchanged.
The retry must earn acceptance and is not presumed successful because the
model is a control. No additional control coupling is introduced.

At 19:10 UTC, the higher-chi critical control continues to show entropy
drift under two-site updates. Add the same fixed-chi one-site variational
refinement already used for the candidate, separately at chi16,24,32,48,
with 32 minimum and 600 maximum sweeps. Its sequence is reported separately
from the two-site sequence and must pass all unchanged gates. This checks
the optimization procedure rather than accepting an expected c=1/2 value.
