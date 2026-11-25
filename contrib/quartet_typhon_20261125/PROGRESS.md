# Quartet study progress

2026-11-25, 12:09 EST. **Bounded pilot and audit completed; full scale-up gate not passed.**

- Chat: `01a21370-d364-7c60-a006-ee9c4074ddf3`. OpenAI Codex,
  **gpt-6-astra, max**, verified from this chat's local turn metadata.
  No subagents or additional chats.
- Initial upstream: `389eadff7d05752b49711f238142fd9cb15c4e84`;
  clean main, safely pulled. Scope is campaign 1 only: physical d=2,
  delta=3/16, h=1.783, with h=1.772 and exact d=4 XX x XX controls.
- All five IAS jobs finished with allocation exit 0:0. The final queue
  check was empty. No L256/384 or L512 production was launched.
- Pilot **17764244**, `quartet-pilot-20261125a`: 4 CPUs, 16 GB,
  one-hour cap; elapsed 50:58, TotalCPU 3:21:35, MaxRSS 509872 KB.
  Retained chi64/128/256 states. Chi256 was shelved after 24 sweeps at
  E=-133.56448160561783, abs(Delta E)=1.473e-6 and
  abs(Delta S)=.006754, above the 1e-9 and 1e-5 stopping tolerances.
  Chi384/512 were not reached. Allocation success is not convergence.
- Control **17764277** completed through chi256 in 14:18, MaxRSS
  568092 KB. Same-chi midpoint entropies agree with old rounded values.
  The target misses its same-chi reference even at sweep-converged chi128.
- Timing run **17764268** completed in 77 s. MKL stayed at one thread
  in its requested 1/4/8 tests, so equal timings do not establish scaling.
- Centred-matrix run **17764302** completed in 110 s. Final read-only
  run **17764463** completed in 6:32, MaxRSS 238904 KB, with full
  pair/quartet quadratures and centred matrices from both chi256 states.
  Repeated measurements of the unchanged control agree exactly.
- Dense Hamiltonian/full-field/gauge and independent correlator checks
  pass. Forced-deadline small-L restart reaches sector ED within 4e-15.
  Pilot roundtrips preserve energy, local observables and entropy;
  final saved-state norm errors are below 6e-14 and parity is +1.
- Six-file historical table reproduces exactly. The original 48-site
  procedure maps a pure p=4 power to p=5.376, and the exact two-boson
  quartet to p=5.381 over r=32–56. Forward filtering repairs that
  pure-power bias; operator/envelope/reference systematics remain large.
- Final target chi128-to-256 quartet changes are 280% for Q(2,5) and
  93% for Q(3,6), relative L2 with the new-state norm, reference47,
  width32, r24–44. Holding momentum fixed leaves these essentially
  unchanged. Pair signals and fitted momentum also move.
- Centred enlarged matrices remain sensitive to candidate-neutral
  subtraction and reference position. This does not certify charge four
  or a neutral rank. The old x4≈2.05±.2 is uncalibrated; no replacement
  x4, locking verdict, asymptotic c=2 theorem or .05 uncertainty is claimed.
- Published checkpoints: `afc5716f3324d2a239f034cfbfa25560fd639c6e`,
  `6235e09c9f25e77bc48549cc8fe30d03897b75e1`,
  `9586405de2090431a9fe299bb29cb721b44ae365`.
  This final results checkpoint adds terminal accounting and measurements.

See [README.md](README.md) for conclusions and reproduction,
[RUNS.md](RUNS.md) for recipes/resources, [VALIDATION.md](VALIDATION.md)
for checks, and [analysis/final_summary.json](analysis/final_summary.json)
for the compact measured record. Exact submitted sources and small raw
outputs are retained. Large restart MPS files remain on IAS outside fetched
results; no historical raw data or prior scientific sections were overwritten.
