# Large central charge for qubits and unbounded central charge at small dimension

Date: 2026-11-18. Author and AI assistance: OpenAI Codex, GPT-6 family.
Status: research priorities and executable assignments, not new physics results.

The user prioritizes two objectives: large unitary compact Virasoro central
charge at microscopic d=2, and unbounded central charge at fixed small d.
Four new project chats are authorized, all using **GPT-6 Astra, Max reasoning**
(`gpt-6-astra`, `max`). The three earlier chats have completed their rounds
and remain idle. The [latest synthesis](ROUND5_RESULTS_20261118.md) is the
starting checkpoint.

Two chats address each objective. One qubit task pursues the best supported
two-mode mechanism while improving its stability and simplicity; the second
explicitly targets c>2. For small fixed d, the markerless family replaces
further dimension-specific d=20/22 work. Its two immediate gaps are an
arbitrary-n band construction and microscopic mobile sector selection.

## 1. Stabilize two modes in a simple qubit chain

Owned output: `contrib/qubit_stable_two_modes_20261118/`.

Start with the literal nearest-neighbour family

```text
H = sum_j [a X_j Y_(j+1) + b Y_j X_(j+1) + u Z_j Z_(j+1) - h Z_j],
(a,b,u,h) = (9/16, 7/16, 1, 7/4) at V4.
```

The [locking audit](../v4_locking_operators_20261118/README.md) gives the
correct symmetry-allowed vertex lattice and Gaussian covariance, but neither
the actual quartet amplitude nor a defensible scaling exponent at V4.
Finding a stable two-mode region requires these quantities from the same
physical Hamiltonian, not a freely chosen Luttinger matrix.

**First calculation:** identify the first nonzero quartet-generating process
in a controlled member of this family, keeping the full interacting
intermediate states and the dressed-field convention. Combine it with
microscopically matched forward/relative couplings to test the complete
allowed integer vertex lattice. Decide whether a low-complexity detuning can
make the dangerous terms irrelevant or cancel their leading source without
introducing a worse allowed instability. Test the next surviving order if a
leading cancellation is proposed.

Use V4 as the simple target and the controlled blockade regime as an anchor;
do not extrapolate the latter's conclusion to V4 without control. A useful
result is an actual stable parameter region, a physically constrained
exclusion within this family, or a new explicit amplitude/sign calculation
that sharply discriminates the two. Prefer a region with rational interior
points over an isolated long-decimal cancellation. If theory cannot fix the
phase, implement a validated pair-observable/stability diagnostic for a later
converged calculation; do not fit an exponent from the known inadequate data.

Relevant sources: the locking audit and its `AMPLITUDES.md`, `GAUSSIAN.md`,
`LEDGER.md`; `contrib/v4_interpretation_audit_20261117/`;
`contrib/simple_qubit_candidates_20261117/`; `codex/work_C179/`,
`codex/work_C180/`, `codex/work_C189/`, `codex/work_C191/`.

## 2. Find an interacting route beyond c=2 in qubits

Owned output: `contrib/qubit_three_mode_deformations_20261118/`.

The [root-of-unity parents](../qubit_resonance_hierarchy_20261116/README.md)
provide exact resonances for every m, but the known undeformed effective
models have ferromagnetic threshold kernels and quadratic excitations.
They are not m-component conformal liquids. The
[detuning study](../qubit_resonant_liquids_20261116/DETUNINGS.md) found useful
two-valley slices and explains why a putative repulsion-only change is
incomplete: number-changing companion terms occur at the same order.

**First calculation:** at the m=3 parent, classify general physical
Hermitian two-site perturbations transverse to the parent family, modulo
constants, telescoping terms and a common energy scale. Map the allowed
directions into the complete leading constrained many-body operator,
including number-changing terms. Determine whether one can obtain three
occupied low-energy valleys with positive compressibility and relative-mode
stiffness, while avoiding clustering, valley polarization and relevant
pair-transfer/number-changing locking. Start with exact local elimination
and scattering; use bounded finite-size calculations as checks, not as a
central-charge determination.

The novelty is the full physical deformation space rather than another
motion along the known projector curve. Keep the microscopic coefficients
independent of chain length. A restricted obstruction is valuable if it
identifies exactly which physical control is missing; do not promote it to
a universal d=2 bound. If m=3 is structurally excluded, identify and test one
specific escape within the same nearest-neighbour qubit class before merely
raising m. Count literal Pauli terms as well as independent real controls.

Relevant sources: `contrib/qubit_resonance_hierarchy_20261116/`,
`contrib/qubit_resonant_liquids_20261116/`, and
`contrib/qubit_cancellation_classification_20261116/`.

## 3. Extend the markerless nine-state band construction to arbitrary n

Owned output: `contrib/markerless_all_n_20261118/`.

The [markerless bond](../markerless_d9_d10/README.md) has nine physical
material states and a 2n+2-vertex carrier. Its seven certified logical-band
members are n=9,10,11,12,14,16,20. More sampled n would strengthen evidence,
but would not address the arbitrary-n gap directly.

**First calculation:** derive the all-n small-loop response of the literal
markerless carrier, preserving finite express delay and its actual physical
field profile. Transfer the methods of the
[eleven-state loop analysis](../eleven_state_loop_splitting_20261118/README.md)
only after checking the deleted-vertex geometry, source locations, fluxes and
available physical controls. Seek an analytic nonvanishing/interlacing/sign
argument that gives 2n simple roots and separates every dangerous pair sum
for each finite n. An explicit recurrence with a proved sign or a certified
induction is a suitable deliverable; numerical signs alone are not.

No new physical states are allowed for this task. If a biased clock, extra
field or altered port is helpful, first write its literal nearest-neighbour
operator in the nine-letter alphabet and verify the invariant grammar.
Give a coefficient schedule depending on n but not system size, or isolate
the remaining exact identity if a complete theorem is out of reach. Track
whether stronger margins can reduce the severe coefficient hierarchy.
Do not infer many-body CFT stability solely from band nonresonance.

Relevant sources: `contrib/markerless_d9_d10/README.md`, `REFEREE.md`,
`scalar9.py`, `certificate/`; `contrib/eleven_state_loop_splitting_20261118/`;
and `contrib/minimal_d_unbounded_audit_20261117/` for the earlier eleven-state
reduction and source inventory.

## 4. Control the ground sector of the corrected ten-state mobile model

Owned output: `contrib/markerless_mobile_control_20261118/`.

The corrected bond is

```text
h10 = J embed(b9) - t sum_a (|a V><V a| + h.c.)
      - (mu/2)(n_occ tensor 1 + 1 tensor n_occ).
```

Every material energy in `embed(b9)` is weighted by material contact. The
first version with genuine material on-site energies failed its selector.
The corrected version has exact legal-sector invariance and encouraging
one-cell checks, but not a multi-cell selection theorem or a quantified
J/t allowance. The selector's n=9 density gap is about 1.82e-4; the band
certificate allows a normalized perturbation of only 2^-35. These are
different margins and both require the correct physical normalization.

**First calculation:** compare competing cell-length words at the same
occupied/vacancy numbers using the literal corrected model. Include changes
in head number, not just balanced (b-1,b+1) defects or one cell versus a
headless ring. Derive a bulk-subtracted response bound or an explicit
reduced scalar comparison that retains the gapless vacancy fluid. Start at
n=9 with a declared density/ensemble, then identify which estimates extend
to arbitrary n. Avoid dividing by a fictitious volume-independent vacancy
gap or bounding an extensive bulk correction by a local defect cost.

Also test the specific simplification suggested by the referee: move only
the conserved charge field out of the material-contact embedding. Verify
the charge algebra, retune the physical field to the intended logical
field, and determine whether this removes the unwanted position-dependent
bit energy without spoiling the selector or adding states. This is a
proposed modification, not an established repair of all mobile errors.

Produce a quantitative sufficient J/t budget if justified, or a precise
remaining inequality with its first nontrivial computation. Distinguish
zero-port selector control from nonzero-port, charged-sector and CFT
claims. Reuse the earlier mobile response/defect techniques only after
checking their hypotheses for the markerless geometry.

Relevant sources: `contrib/markerless_d9_d10/README.md`, `REFEREE.md`
(including its corrected-bond addendum), `bond9.py`, `bond10.py`,
`mobile_selector.py`, `selector.py`; `contrib/mobile_length_selection_20261117/`,
`contrib/mobile_vacuum_clock_20261117/`, and
`contrib/small_d_defect_stability_20261116/`.

## Common scope and collaboration

- Every Hamiltonian must remain Hermitian, nearest-neighbour and one-site
  translation invariant at the microscopic dimension stated. Count every
  physical marker/vacancy state. Keep parameters fixed as L grows.
- Report finite algebra, numerical evidence and infrared assumptions
  separately. A root count, parent degeneracy, entropy slope or encoded
  dimension alone is not a conformal central charge.
- Simplicity is a scientific objective: report independent real controls,
  literal operator terms and normalized nonzero coefficient spread. Remove
  constants, telescoping terms and common positive scale. Distinguish small
  emergent velocities from unequal microscopic coefficients.
- The previous battery-saving stop on long DMRG remains in force. Use
  analytic work, exact/symbolic algebra, saved results and bounded one-worker
  checks. Do not launch background production or large scans. New chats are
  research capacity, not authorization to restart stopped production.
- Read `AGENTS.md` and `CONTRIBUTING.md`. Each chat owns only its new output
  folder; the coordinator owns this file and shared status pages. Do not
  spawn additional chats or subagents in this round.
- Serialize Git mutation in this shared checkout using the atomic directory
  lock `/tmp/cd-conjecture-git-sync.lock`. Acquire with `mkdir`, record an
  owner, hold it for the complete stage/commit/pull/push transaction, and
  release only your own lock. If busy, continue independent work or wait;
  do not delete another owner's lock. Never stage another chat's files.
  Preserve other uncommitted work; if it prevents rebase, fetch and state
  the pending integration rather than stashing it or claiming a clean sync.
- Finish a concrete derivation or reproducible calculation, not only a
  literature summary or future plan. End with the strongest established
  result and the next calculation that could change the conclusion.

## Dispatch record

The coordinator will record the four returned chat IDs and observed startup
states here after creation. The requested configuration is GPT-6 Astra Max
for all four chats; it must not be silently replaced by another model.
