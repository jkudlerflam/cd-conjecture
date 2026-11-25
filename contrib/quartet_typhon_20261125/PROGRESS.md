# Quartet study progress

2026-11-25, 11:15 EST. **Active execution; estimator calibration failure found.**

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
- Scale-up gate: correctness, recoverable checkpoints, measured resources,
  convergence of correlators, and calibrated operator extraction. Energy
  convergence alone is insufficient. Higher-chi/longer-chain work will
  be selected from these checks, not scheduler exit status alone.

See `RUNS.md` for runs and `VALIDATION.md` for validation outcomes. Large
states remain on IAS outside fetched results. Historical files are intact.
