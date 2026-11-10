# Test of the apparent plateau — 2026-11-10 18:48 UTC

The completed chi128 two-site state passes all gates: E=-0.676213102281,
S=1.089661802, xi=134.506026. The ordinary-entropy slope falls from about
0.504 at chi64->96 to about 0.308 at chi96->128. Thus the apparent half-unit
slope is not stable under this extension. This is a useful discrimination,
not a reason to change the fitting window or acceptance threshold.

Extend the identical real-frame state to chi192 and refine chi128 with
the one-site method. Set the new chi192 local Krylov limit to 500 because
many chi128 local solves reached the 250-iteration cap even while passing
the explicitly measured residual threshold. The residual gate remains
1e-8; no numerical or physical acceptance threshold is relaxed.

Assess whether entropy and transfer lengths approach saturation. A chi256
extension will be considered after chi192 completes. A gapped crossover
is possible; neither c=1/2 nor c=3/2 is assigned by these calculations.
