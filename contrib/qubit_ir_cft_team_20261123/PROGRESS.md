# Fixed-qubit IR investigation: progress

All timestamps are UTC. This is a first investigation of the fixed Pauli
Hamiltonian in `README.md`; no parameters may be moved. Original MPS states
are unavailable in Git. Long DMRG and background production remain prohibited.

## Team registered — 2026-11-23 22:22 UTC

The collaboration tool returned the following actual canonical agent IDs (it
did not expose a separate UUID). The three workers were explicitly launched
with `model=gpt-6.2-astra`, `reasoning_effort=max`, with no nested delegation.

| Agent ID | Role | Current task | Owned files |
|---|---|---|---|
| `/root` | Manager and synthesis | Integrate evidence; track claims; serialize numerical work and Git | This directory's top-level files |
| `/root/theory` | Symmetries, continuum dictionary, stability | Independently derive exact symmetries and conditional allowed vertices | `theory/` |
| `/root/evidence` | Numerical evidence and mode identification | Inventory saved observables, reproduce fits, propose bounded discriminating calculation | `evidence/` |
| `/root/auditor` | Independent AI audit | Check literal model and normalizations; independently reproduce a decisive calculation | `audit/` |

Achieved: read repository contribution instructions; verified clean `main`
at entry; launched the requested three workers. Git synchronization uses the
shared atomic `/tmp/cd-conjecture-git-sync.lock` directory and an owner marker.
Only the manager may mutate Git. No numerical production has been launched.

Unresolved: microscopic CFT identification, two-mode status, charge lattice,
Gaussian metric/mixing, velocities, and leading locking perturbations.
Auditor verdict: pending. Next: read source evidence, establish a claim ledger,
and authorize at most one bounded numerical worker at a time.

## First independent findings — 2026-11-23 22:26 UTC

- Startup pull completed: already up to date after the coordinator's
  deployment commit `a553796e`. Initial restricted Git write failed; the
  authorized elevated retry succeeded. The shared lock was acquired and
  released only by its owner.
- Evidence worker has the sole numerical slot: original-frame periodic
  exact diagonalization at L=8,10,12, one thread, 90-second cap. No other
  substantial computation is authorized concurrently.
- Inventory: target correlators are saved as fitted summaries, without the
  fitted pair momentum or raw complex distance-dependent arrays. The
  export's three zero check fields are placeholders, not measurements.
  Target odd-even data list L=64,96,128, all with positive differences;
  sign-changing examples belong to other parameter points.
- Theory and auditor independently agree that staggered distance-two
  chirality is translation-odd, whereas smooth valley imbalance is
  translation-even. Their proposed identification fails in the stated
  unbroken two-valley dictionary.
- Theory and auditor independently derive conditional Gaussian bounds:
  a matched momentum-pi pair exponent p below one implies a relevant
  doubled allowed vertex of dimension 2p. If it is the undressed
  inter-valley pair, the vertex is the primitive quartet. The target
  table does not establish that momentum/primary assignment.
- With the additionally matched one-body exponent, all generic mixed
  backscattering harmonics can be bounded even with Gaussian mixing.
  Details and exact hypotheses are being written and independently checked.
- The unit-normalized raw-density and pair assignments are inconsistent
  at their quoted central values even allowing Gaussian mixing. This
  rejects that joint dictionary, not the two-mode candidate itself.

Auditor verdict: accepts the exact translation-parity mismatch and the
conditional algebra; target IR identification remains unresolved. Next:
finish bounded ED, independently reproduce decisive data, and integrate
the conditional charge lattice with a quantitative measurement protocol.

Model provenance correction: the coordinator clarified that the requested
worker setting was specifically `gpt-6-astra / max`. The manager actually
launched the workers with `gpt-6.2-astra / max`, as recorded above, despite
`gpt-6-astra` being available. This was a selection error, not an availability
constraint. Existing work is preserved as instructed; no agent was relaunched
or result relabeled. The manager's requested setting remains GPT-6 Astra Max.

## Scientific investigation complete — 2026-11-23 22:40 UTC

All three workers have completed their assigned research and integration
follow-ups. `/root/evidence` reviewed the final synthesis and TODO for data
availability, exact windows, quoted numbers and zero-mode constraints.
`/root/theory` supplied complete symmetry/charge-lattice/mixing derivations
and an explicit same-exponent alternative operator assignment.
`/root/auditor` independently checked the decisive claims and required
corrections to pair quadrature recovery, operator rank, matching and
all-scale inference; the revised protocol is accepted.

Achieved:

- Reproduced target saved entropy regressions and audited the measurement
  code and accessible-data inventory.
- Computed complete L=8,10,12 periodic spectra and all-momentum static
  probe sums, retaining low-lying operator-resolved residues. The first
  run used 0.151 seconds with a 90-second cap.
- Independently reproduced L=8 sector minima and four pi probe residues
  with dense tensor-product Pauli matrices; largest discrepancy 9.8e-15.
  The first independent run used 0.30 seconds with a 30-second cap.
- Independently derived the conditional Gaussian mixing square, density
  consistency bound, neutral counterexample, and causal speed ceiling.
- Added the user-requested prioritized `TODO.md`, the claim ledger and a
  separate `SHARED_SUMMARY.md` for the coordinator. All new work stays in
  the team's contribution directory.
- Manager reran all four scripts sequentially at 22:39 UTC. All passed;
  the two independent audit outputs exactly matched saved JSON. The
  system Python lacked NumPy, so verification used the existing bundled
  Python 3.12.14/NumPy 2.3.5 environment; no installation was performed.
  Numerical calculations are finished and no job is running.

Final audit verdict: accepted exact and numerical checks, accepted
conditional theory with explicit hypotheses, rejected the identified
operator/mode/normalization overclaims. The microscopic CFT, number of
asymptotic modes, compactification, metric/mixing, actual speeds and
quartet coupling remain unresolved. No CFT identification is manufactured.

Next required action: manager's serialized commit, rebase/push and upstream
verification. Scientific work is complete; Git completion is not yet
claimed at this checkpoint.
