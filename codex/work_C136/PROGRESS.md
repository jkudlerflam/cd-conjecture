# C136 progress

2026-11-07.

[Proved, inherited] Read the C129 reductions and the C134 real-axis proof and referee report. C134 is independently refereed as correct; its earlier status files are stale. Its half-filled doubling at odd L is not the matched-sector doubling required here.

[Proved] The negative kernel gives a globally increasing finite counting function at every filling. The restricted integral operator has the same Neumann bound, independent of its interval. The full-line discrepancy used in C134 does not stay small below half filling because the missing counting strips have width (1-gamma/pi)(1-2n/L).

[Work record] Investigating a clipped cumulative distribution equation to compare finite roots to the finite-interval density without endpoint denominators. All writes and Python execution stay in work_C136. Large scans will be prepared for SLURM. No cluster access is available.

[Assumed, unresolved] No constant uniform through Delta=-1 has yet been proved. The dilute estimate n^3/L does not cover all sequences with n/L tending to zero.

## Completed work

[Proved] Obtained an explicit finite-boundary theorem uniform down to zero density, then closed its half-filled density limit using exponentially decaying resolvent bounds, comparison to the full-line sech density, and variation bounds for the dressed velocity. Task 1 now holds at every filling and both length parities. The compact physical theorem has C(delta)=6 * 10^15 delta^(-25/3), for every field and DM angle. The proof and smaller constants are in PROOF Sections 2-6 and 9. No independent referee pass was performed.

[Proved, limitation] Tracking the C134 constants gives a cutoff of order (1+Delta)^(-1) and a remaining large-size coefficient of order (1+Delta)^(-3/2). Neither this nor the new constants becomes uniform beyond the crossover length. The uniform ferromagnetic estimate remains missing. PROOF Section 8 states the exact additional estimate needed.

[Numerical] All 1,032 paired Bethe checks pass, with ED at every sector through L=12, and sampled Bethe pairs through 514 sites. Fields, twists, kernel moments, and constant arithmetic pass. Of 54 continuum diagnostics, 43 pass after refinement; eleven narrow-kernel cases remain underresolved. Failures are retained and quoted in CHECKS.md. No full 98,304-case grid or large interacting computation ran locally.

[Work record] RESULT.md, PROOF.md, CHECKS.md, LITERATURE.md, STATE.md, scripts, and outputs are complete in this folder. Two SLURM scripts prepare the full grid and continuum refinements. Neither was submitted. Task 2 and the full uniform assembly are unfinished; the optional positive-anisotropy task was not pursued.
