# Counterexample continuation, 2026-11-09

AI-assisted research by OpenAI Codex (GPT-6 family), with independent
parallel Codex reviews. No human referee or new counterexample is claimed.

## Starting point and priorities

The checkout was synchronized to `a45ecf5` before work. The latest project
status separates three routes:

1. **Most mature all-ring construction:** the d=13 mobile-cell bond with
   predicted total Virasoro central charge c=14. The finite construction is
   extensively checked; the physical infrared argument retains stated
   hypotheses. Reconcile its finite parameter prescription with the
   independent d=13 red-team report before treating the displayed member
   as a complete certificate.
2. **Highest-value direct test:** C173's nearest-neighbour qubit bond at
   the Fibonacci cusp, with candidate c=2 at d=2. C180 supplies the
   symmetry and field-theory analysis. Existing entropy estimates are
   unconverged. Valley miscibility, number-changing corrections and a
   converged microscopic central-charge measurement remain open.
3. **Stronger long-term statement:** unbounded c at fixed d=4464, whose
   final microscopic-to-infrared step is still conditional.

These are research priorities, not new theorem claims. The original
conjecture concerns the total c of unitary CFT sectors, including sectors
with distinct velocities: [Latorre and Sierra](https://arxiv.org/abs/2403.17242).

## Work in this continuation

- Independent exact few-body scattering for the Fibonacci candidate,
  including its generated attractions and a no-attraction calibration.
- Independent enumeration of its perturbative reduction and
  number-changing terms, with tests at the actual candidate parameters.
- A finite-ring boundary-condition and symmetry diagnostic audit.
- An explicit repaired scale ledger for the d=13 construction, preserving
  the distinction between a finite certificate and its infrared input.

The existing C173/C176/C177 cluster DMRG tasks will not be duplicated.
Each new calculation will include source code, saved results and its
scope. Contributions live in separate `contrib/*_20261109/` folders;
upstream collaborator files are left intact unless a checked correction
is warranted. GitHub synchronization happens at reviewable checkpoints.

## Initial evidence

- [Current project status](../../STATUS.md), REPORT sections 9.137-9.138.
- [C173 candidate](../../codex/work_C173/RESULT.md) and
  [C180 field theory](../../codex/work_C180/RESULT.md).
- [d=13 red-team report](../redteam_d13/README.md).
- Open GitHub issues 1-5 concern the earlier construction's matching,
  defect comparison and infrared estimates. Open PR 6 concerns a
  replacement bound, not the new qubit candidate.

This is the initial checkpoint; results will be appended as independently
checked calculations finish.
