# Quartet study progress

2026-11-25, 11:46 EST. **Active execution; estimator calibration failure found.**

- Chat: `01a21370-d364-7c60-a006-ee9c4074ddf3`.
- AI: OpenAI Codex, **gpt-6-astra, max**, verified from this chat's
  local `turn_context` metadata. No subagents or additional chats.
- Initial upstream: `389eadff7d05752b49711f238142fd9cb15c4e84`;
  clean `main`, pulled and equal to `origin/main`.
- Scope: campaign 1 of the November 25 plan. Fixed physical qubit
  Hamiltonian, delta=3/16, h=1.783; h=1.772 and two-boson estimator
  controls. No dynamics or d=10 campaign.
- IAS inventory: user queue empty; project run area contains only the
  two four-site validation runs. No available quartet warm starts there.
- Six-file historical table reproduced exactly. The 48-site filter biases
  a pure p=4 signal to p=5.376, and an exact XX x XX quartet to p=5.381.
  Forward-filtering repairs this pure-power bias, but the target's operator
  and harmonic-envelope systematics remain unresolved. No replacement
  charge-four dimension is claimed. See README and analysis files.
- Dense Hamiltonian/field/gauge, complex correlator and save/load checks
  passed; production-offset matrix checks passed on an independent L=24
  parity cat (704 entries, error below 7e-16).
- Pilot submitted: `quartet-pilot-20261125a`, Slurm **17764244**, running
  on typhon-node15 (4 CPUs, 16 GB, one hour). Current low-chi sweep timing
  is about 16 seconds with process RSS below 200 MB; not a high-chi estimate.
- Next: complete pilot/checkpoints, compare numerical thread counts,
  measure enlarged operator matrix and assess justified follow-up. Full
  L=256/384 chi=2048 expansion is not yet justified by estimator calibration.
- Target chi64 saved, warm restart into chi128 successful; chi64 is explicitly
  unconverged (41 sweeps, E=-133.56333751015325). Roundtrip errors are zero
  for E/Z/S and below 2e-14 for overlap. A separate forced-deadline L8 test
  reloads and converges to the sector ED energy within 2e-15.
- Thread test **17764268** completed, exit 0:0, 77 s, 211648 KB Slurm
  MaxRSS. MKL reports one thread after each requested 1/4/8 setting, so the
  equal timings do not establish parallel scaling. Raw evidence retained.
- Bounded c=1 control **17764277**, `quartet-c1control-20261125a`, submitted:
  L128 even, h=1.772, chi64/128/256, 1 CPU, 4 GB, one hour. This adds
  matching enlarged probes to test neutral mixing, not high-chi production.
- GitHub checkpoint published: `afc5716f3324d2a239f034cfbfa25560fd639c6e`.
- Both chi128 states converged by sweep tolerances and full quartet/matrix
  outputs were fetched. Target E=-133.56423306491507, Smid=1.8315498;
  this differs from the old chi128 Smid=1.67048 and is a convergence/basin
  warning, not portability certification. Control Smid=.75452675 agrees
  with the old chi128 .75453. Both ladders have advanced to chi256.
- Initial matrix diagnostics show order-one sensitivity to a fixed bond
  probe subtraction, also in the control. Spatial origins of different
  shapes differ; an aligned follow-up centres all Hermitian probes at
  i+1/2. Its 648 additional independent matrix checks pass below 3e-16.
  The aligned recipe is prepared, not yet submitted. No charge-rank claim.
- Future wrapper revision separates optimization and measurements into
  bounded processes, checks remaining time before stage launch, and
  preserves the undefined first energy change as null. A one-sweep L24,
  chi4 workflow test passes and correctly remains labelled unconverged.
- The initial summary's Sbar label denoted a central-half average;
  `analysis/run_summary.json` recomputes the correct central-quarter
  average from the raw profiles. Raw files are preserved.
- Control 17764277 completed through chi256: exit 0:0, 14:18, MaxRSS
  568092 KB. E=-133.2253186480565, Smid=.75457451 (same-chi reference
  .75457). The small Q(2,5) signal changes 4.6% from chi128 to 256.
- Aligned measurement 17764302 completed, exit 0:0 in 110 s, MaxRSS
  204512 KB. Target chi128 small-quadrature signals remain strongly
  sensitive to bond-probe subtraction and reference position, also
  qualitatively present in the control. No rank/charge assignment certified.
- Target pilot remains at chi256, still moving appreciably in E and S.
  The one-hour deadline will preserve a partial stage if necessary. The
  estimator and wavefunction convergence gates do not support full scale-up.
- Published second checkpoint: `6235e09c9f25e77bc48549cc8fe30d03897b75e1`.
- Scale-up gate: correctness, recoverable checkpoints, measured resources,
  convergence of correlators, and calibrated operator extraction. Energy
  convergence alone is insufficient. Higher-chi/longer-chain work will
  be selected from these checks, not scheduler exit status alone.

See `RUNS.md` for runs and `VALIDATION.md` for validation outcomes. Large
states remain on IAS outside fetched results. Historical files are intact.
