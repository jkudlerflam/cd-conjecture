# Coordinator handoff and proposed shared status text

The team completed an audited first investigation of the unchanged
physical-qubit Hamiltonian
`H=sum[(19/32)XY+(13/32)YX+ZZ]-(1783/1000)sum Z`, with full on-site fields.
Its verdict is **finite-scale near-two numerical candidate, conditional
two-boson effective description, asymptotic microscopic CFT unresolved**.

## Proposed brief status entry

The [fixed-target team investigation](README.md) independently reproduces
the saved iMPS/OBC entropy slopes 1.97265 +/- .00986 and 1.95401 +/- .07500.
It corrects the identification of staggered chirality with smooth valley
imbalance: their translation characters differ. Under a matched generic
two-boson dictionary, the quoted pair/one-body exponents make the quartet
the only relevant generic harmonic even with mixing (x4=1.8386), while
relative/mixed harmonics obey x>=2/.5083. But the pair momentum/raw target
correlators and quartet coupling are unavailable. An explicit alternative
operator assignment reproduces the two scalar exponents with every
generic lock irrelevant. The raw-density and pure-pair assignments are
also mutually inconsistent at unit density normalization, even with
mixing. The two mode count, compactification, physical speeds and
asymptotic Virasoro c therefore remain open.

## Concrete new evidence and audit

- Original-frame periodic exact diagonalization at L=8,10,12 finds the
  same low pi transition in pair, bond and chirality probes; the strongest
  separation-three static pair peak need not be at pi. This is finite-ring
  operator evidence, not an independent mode count or CFT tower.
- A separate auditor built dense literal Pauli matrices and reproduced
  all L=8 parity/momentum sector minima and selected complex residues
  within 1e-14, including field, gauge and boundary checks.
- Independent Gaussian algebra incorporates all symmetry-allowed mixing,
  gives a causal microscopic speed ceiling 153e/16 (about 26, not a speed
  estimate), and derives a conditional normal/anomalous pair ratio that
  could measure the matched quartet coupling.
- [CLAIMS.md](CLAIMS.md) and [audit/AUDIT.md](audit/AUDIT.md) separate exact,
  numerical, conditional, rejected and unresolved statements. The audit
  is independent AI checking, not human review.

## User-requested prioritized to-do list

[TODO.md](TODO.md) is reviewed and ready to share on GitHub. Its five
necessary identification inputs are:

1. Recover the fitted pair carrier and raw complex rows, then obtain the
   anomalous/cross-separation quadrature data needed to match the primary.
2. Identify independent small-q bulk current oscillators and velocities,
   treating each momentum sign and the real-gauge folding correctly.
3. Match the full response/covariance, including density normalization,
   phase gradients, mixing and analytic high-energy contributions.
4. Establish local charge/winding labels, zero-mode constraints and any
   twisted sectors before assigning a torus, orbifold or rational CFT.
5. Test allowed locking amplitudes and longer-scale convergence, rather
   than extending entropy fits without discriminating observables.

Each entry states precise conventions, competing outcomes, the minimum
measurement, saved-data versus production requirements, and a resource
estimate or an explicit unknown. It does not assume the target is
asymptotically critical. No new long production is authorized or running.

## Completion record

Only `contrib/qubit_ir_cft_team_20261123/` was edited. The four reproducible
checks pass; see [verification.json](verification.json). The numerical
work was bounded, serialized and single-threaded. No previous contribution
or shared top-level summary was modified by this team.

Actual worker IDs and models are in [PROGRESS.md](PROGRESS.md). The manager
selected `gpt-6.2-astra/max` for the workers instead of the specifically
requested available `gpt-6-astra/max`; that selection error is disclosed,
and the coordinator instructed preservation of useful existing work.
The Git synchronization checkpoints and final commit identifiers are
reported in PROGRESS and the manager's completion message.
The scientific checkpoint is `1f8c2051625fc81ecc2f29eb3a5dfd224989c80d`,
successfully pushed to `origin/main` and verified there at 22:41 UTC.
