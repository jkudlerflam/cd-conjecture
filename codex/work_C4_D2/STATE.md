# D2 local handoff, completed 2026-10-21

This directory is a closed-book numerical attack on the d=2 c-d conjecture.
No pre-existing project files were read or listed. All created files are in
this directory. The parent project STATE.md was not accessed or modified.

The verdict is no established counterexample, not a proof of c_T<=1. Large
finite-size entropy slopes persist in some candidates, but none has a
validated conformal interpretation with consistent spectra and identified
stress states. Several candidates remain unresolved beyond L=22.

## Timing and completed work

The local run began at 20:43:24 UTC. All ED was finished by 21:28:32 UTC;
final numerical verification was complete at 21:29:50 UTC. That is 46 minutes
26 seconds including development and analysis to that checkpoint, below the
75-minute wall-time budget. The final handoff files were written immediately
afterwards. Final artifacts were inspected at 21:32:20 UTC, 48 minutes
56 seconds after the start. Multiple local processes ran concurrently; this is elapsed wall
time, not an aggregate CPU-time accounting. No local ED job remains running.
No cluster job was submitted and no network sources were consulted.

There are 1,338 global samples (138 structured and 1,200 random canonical),
159 additional boundary samples, 340 initial optimization evaluations and
130 larger-size optimization evaluations. The 238-point pilot is included
in the 1,338 and must not be counted again. There are 212 L=14,16 follow-ups,
254 velocity-screen samples, and 16 selected parameter points with all-momentum
ground searches and the three diagnostics through L=22. Four calibration
chains, two near-ferromagnetic controls, and two finite-density XX controls
also reach L=22.

The literal all-above-one follow-up was not completed. Of the global samples,
797 already have c_eff(L=12)>1. Screening and targeted larger-L calculations
were used instead. The cluster scripts extend every point if run there.
Boundary point b43 was completed through L=18; its slow, nearly degenerate
L=20 run was interrupted, with no L=20 result used. The other interrupted
trials were superseded ferromagnetic-control runs described in README.md.

## Main results, do not reinterpret raw slopes as c

The largest raw c_eff across sampled sizes is 2.7033363661 at point 830,
L=10. Its regression SE is 0.1185474643, not a finite-size error bar. The
same fixed couplings give c_eff=1.4373281715 at L=22.

The largest raw L=12 slope is 2.2687097506 at point 854, falling to
1.4165058809 at L=22. The largest raw L=22 slope among the selected points
is 1.7023424011 at point 1256. Excluding short blocks changes its L=22 slope
to 1.3495580639. Its half-chain scaling and excitation velocities are strongly
window dependent. Its selected matrix-element ratios are about 59.10 and
4.079, not identified stress central charges.

The initial optimized winner opt68 falls from 1.6184387733 at L=12 to
1.2567135941 at L=22. The second winner opt16_121 falls from 1.8793533156
at L=16 to 1.6634797428 at L=22; its half-chain scaling slope is
0.9054 with a 0.2719 sensitivity scale, and its matrix-element ratios are
0.1949 and 5.9208. It is not a consistent CFT determination.

The velocity-screen leaders 1004 and 1057 have half-chain scaling slopes
1.0104 +/- 0.037 and 0.9972 +/- 0.022, respectively, consistent with one.
The JSON metrics contain the full-precision sensitivities. No generic point
is assigned a central charge above one.

## Estimator failures that must remain visible

The literal left-site-field density with largest overlap does not select the
Ising stress state at these sizes. At L=22 it yields 0.8736784780. Splitting
the field symmetrically leaves H unchanged and selects the stress branch,
giving 0.4848743524 and extrapolating toward 1/2. An energy-density improvement
can compete at the same scaling order as the stress matrix element.

Even the c=1 XX chain at nonzero field has zero-mode/shell corrections that
invalidate the naive polynomial Casimir fit. Exact ring energies at L=8,...,22
give c_E=-0.3866994658 for h_z=0.3 and 1.7122577963 for h_z=1.1, using the exact
velocity. Largest individual Fourier overlaps also oscillate because of
finite-size splitting and charge-sector effects. Do not use disagreement of
these raw formulas alone as proof that a candidate is nonconformal.

The exact ferromagnetic boundary has raw c_eff=1.3600315087 at L=22 in the
half-filled symmetric state, but E0=-L, nonzero-momentum h_n|0>=0, and
one-magnon gap 4(1-cos(2*pi/L)), hence z=2. This gives an explicit false
positive for an entropy-only search. The derivation is in methods.md.

## File map

- d2.py is the ED toolkit and CLI. It supports generic J and field, full
  momentum searches, entropies, energy fits, both Fourier directions, scans,
  and local optimization. J is in Pauli normalization. Default field density
  split is 1/2; --field-split 1 restores the literal left-site density.
- report.md is the verdict and main table. methods.md derives conventions,
  formulas, and their limitations. calibration.md, candidates.md, and
  scan_summary.md are regenerated by summarize.py and write_report.py.
- results/best_candidates.json consolidates all 16 selected L=22 trajectories,
  parameters, and metrics. results/best_candidate_points.json is the minimal
  point list for a cluster rerun. results/candidate_metrics.json contains the
  fit-window sensitivity estimates, not rigorous confidence intervals.
- results/calibration_fixed.json is the canonical calibration. Do not use
  results/calibration.json, which contains a superseded degeneracy failure.
- results/ferromagnet_half.json is the canonical near-ferromagnetic control.
  Do not use results/ferromagnet.json, whose fixed momentum was wrong for
  L=2 mod 4. results/ferromagnet_correct.json is a valid but interrupted
  slower run through L=20, superseded by the charge-resolved result.
- results/global_scan.json, boundary_scan.json, optimization.json,
  optimization16.json, shortlist.json, global_shortlist.json, and
  velocity_screen.json retain all screening data and parameters.
- results/diagnostics.png and diagnostics.svg are the inspected summary plot.
- exact_controls.py, free_fermion_control.py, density_check.py, and
  verify_results.py contain independent checks. make_boundary_points.py and
  screen_velocity.py reproduce the additional search stages.
- cluster_array.slurm and submit_cluster.sh are prepared, not submitted.
  The wrapper batches arrays and supports D2_POINTS, D2_BATCH_SIZE,
  D2_CONCURRENCY, and a user-specified authorized partition. It is CPU code.

## Verification and next steps

The direct L=6 spectrum, projection, and Fourier-action checks are at about
1e-14. The final 68 comparisons of full-space and momentum-space runs have
maximum energy difference 3.322e-13 and unique-ground-state entropy difference
1.783e-11. The largest stored ground residual is 1.042e-10; the largest stored
target residual is 2.600e-9. Sixteen exact XX-ring energies agree with ED to
7.106e-14. On-site rotation, density-gauge invariance of H, and charge-block
entropy were checked. Final Python compilation and bash syntax checks passed.

For continuation, prioritize zero-mode-aware spectral analysis and independent
larger-system methods for 1256, opt16_121, and 1282, rather than treating the
raw Casimir or largest-overlap formulas as calibrated for every generic chain.
The coarse all-above-one list can be extended using the cluster array. No
claim about the supremum over the continuous eight-parameter family was proved.
