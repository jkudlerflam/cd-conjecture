# Registered bridge continuation

OpenAI Codex, 2026-11-10. Written before this folder's numerical output.
This is a full microscopic nearest-neighbour qubit search. The previous
contribution remains read-only. No new ED screen or projected model enters
selection. No cluster job or publication is attempted.

Use Delta=Jx-Jy=1, q=Jx+Jy, Jz=r/4. Initial exchange choices and fields:

| r | D | Fields hz |
|---|---|---|
| 2 | 3/16 | -3/4, -11/16, -5/8, -9/16 |
| 2 | 1/4 | -7/8, -13/16, -3/4, -11/16 |
| 5/2 | 1/10 | -15/16, -7/8, -27/32, -13/16 |

For each use q=0,1/16,1/8 and Jx=(1+q)/2, Jy=(q-1)/2: 36 Hamiltonians.
These span the preferred literal coefficient spread <=3 and exploratory
spread <=10. Exact absent terms are allowed. The q offsets and cancellation
fraction |q|/(|Jx|+|Jy|) are reported; no basis rotation is used to advertise
a smaller hierarchy. Store singular values and symmetric/antisymmetric norms.
There are six expanded Pauli strings; q=0 has four grouped terms and nonzero
q has five. Field resolution is initially1/16, except the last family also
has1/32 intervals; these are sampling resolutions, not critical tolerances.

Stage A: L=24, both exact parities, chi64, seed17; C190 solver and previous
validated wrapper are imported without changing updates. Sweep convergence
targets remain local residual<=1e-8, delta E<=1e-8, delta S<=1e-3; final
local solves target1e-10. Save every sweep checkpoint in a distinct directory.
Record the lowest variational-energy branch only after both parities finish.
Report both sectors and any failed gate. Reuse the identical-code XX/twoXX
calibrations in the previous contribution; repeat the exact same r2,D3/16,
hz=-1/2,q=0 point and XX L24 as plumbing controls in this new destination.

Density refinement: for each of the nine (r,D,q) choices, find adjacent
fields whose selected central densities straddle0.275. Insert their exact
rational midpoint, up to two rounds, stopping if any point is within0.015
of0.275. At most18 extra fields. If there is no bracket, report that fact
instead of silently widening it. This is a density anchor, not evidence
that raw density represents mobile particles: virtual dressing matters.
Shell or first-order jumps must remain visible.

Stage B: select at most two points per (r,D) family, six total, after Stage A
and any registered refinement. One q=0 point minimizes |n_mid-0.275|; the
nonzero-q point is chosen from those within0.08 of0.275 by the score
clip(min(c_profile_1/4,c_profile_3/8),0,3) - |c_profile_1/4-c_profile_3/8|/2
-5|n_mid-0.275| -2*density_bulk_range. If none meets the density interval,
choose the closest density and explicitly call it an unflagged control.
Write the selected list before L48 output. L48 uses chi128 and both parities.
Two-size slopes use the already fixed central-three and fractional windows.
A high single-size fit is not a pass; a plausible promotion requires size
growth in both windows, controlled chi and sweep errors, calibrated controls,
and no obvious ordered-cat or phase-separation signature. Any c>1 can qualify;
the screening flag1.2 is not an exclusion threshold.

Coordinate L96 or larger checks with the parent before launching them. This
initial child task uses at most two concurrent single-BLAS-thread workers.
