# Research coordination and a residue-sector check

Date: 2026-11-16. OpenAI Codex, GPT-6 family. This records the first
coordinated follow-up to three project chats and an independent check of
their newly proposed sixth-order qubit models. No human review is claimed.

## Reviewed checkpoint and assigned follow-ups

The [chirality-chain field scan](CHIRALITY_FIELD_DIAGRAM_20261118.md)
assigns one new GPT-6 Astra Max chat to the full h-axis phase diagram of
the simple four-term qubit Hamiltonian. It prioritizes the proposed
c=3/2 transition and nearby apparent c about 2 region while checking
the low- and high-field phases. Startup is confirmed; no new phase
claim or long production run is implied.

The [completed stability results](STABILITY_RESULTS_20261118.md) give
volume-independent finite-order local interaction bounds and exact
nonzero-port frozen sector selection for modified d=9; a connected qubit
quartet response and a calibrated protocol for the missing dressed medium;
and a uniform half-occupancy response channel for original mobile d=10.
All four GPT-6 Astra Max chats completed the
[long-distance stability round](STABILITY_ROUND_20261118.md). No new
microscopic CFT or improved coefficient hierarchy is established.

The [completed focused-round results](FOCUSED_RESULTS_20261118.md) give an
arbitrary-n leading-band theorem on a modified d=9 alphabet, finite d=10
mobile sector certificates and an exact charge-field improvement, finite
qubit quartet amplitudes, and a leading-model crystal-preemption theorem
for the tested three-valley route. All four chats have completed; no new
microscopic CFT central charge is established.

The [focused four-chat round](FOCUSED_ROUND_20261118.md) prioritizes large
c at d=2 and unbounded c at fixed small d. It assigns two GPT-6 Astra Max
chats to qubit stability and physical three-mode deformations, and two to
arbitrary-n markerless bands and the corrected d=10 mobile ground sector.
These are research assignments, not new phase claims.

The [latest completed results](ROUND5_RESULTS_20261118.md) include the
newly published markerless d=9/d=10 family, the separate single-rail
nine-state three-flavour instability, the eleven-state analytic loop
splitting, and the corrected qubit locking criteria. All three chats
completed that round; no new production run is implied by this summary.

The [2026-11-18 GitHub source check](GITHUB_D10_AUDIT_20261118.md) records
the expanded search following the user's clarification, the new nine- and
eleven-state results, and their next calculations. The d=10 source was
unidentified at that checkpoint; publication of
[the markerless family](../markerless_d9_d10/README.md) later that day
supersedes that source-search result.

The [smaller-d fourth-round plan](SMALLEST_D_FOLLOWUP_20261117.md) records
the user's new priority, the source audit of the reported d=10 family,
and the active tasks on qubits, <=10-state mechanisms and the existing
eleven-state construction before the markerless source was published.

The [third-round results](ROUND3_RESULTS_20261117.md) summarize the latest
completed work on the V=4 bulk phase, the two-term zero-field qubit line,
and mobile length selection. The [review and assignments](LATEST_REVIEW_20261117.md)
record why those tasks were chosen, following the earlier
[simplicity-first assignments](FOLLOWUP_20261117.md).
The checkpoint below retains the preceding round.

| Direction | Completed checkpoint | Next assigned task |
| --- | --- | --- |
| Establish a counterexample | [Microscopic qubit audit](../qubit_convergence_20261116/README.md): the next large-chi test remains undecided; accessible checkpoints and sufficient production memory are absent | Resonant qubit liquids: study the actual many-body physics and admissible detunings of the new three-/four-valley effective branches |
| Bound c at fixed d | [Cancellation classification](../qubit_cancellation_classification_20261116/README.md): complete stated low-order normal forms and exceptional longer-range hopping; no upper bound | Resonance hierarchy: derive exact root-of-unity parent curves and study general resonance order, beginning with m=5,6 |
| Unbounded c at small fixed d | [Frozen defect theorem](../small_d_defect_stability_20261116/README.md): positive nonzero-port defect costs for d=20 on compatible rings | Mobile matching: derive the literal d=22 contact representation and control a bulk-subtracted defect or auxiliary-excitation cost with the gapless vacancy sector retained |

These are work assignments, not claims that the next results are proved.
The original qubit convergence registration remains intact; its compute
blocker is not repaired by investigating another candidate.
The assigned output directories are respectively
`contrib/qubit_resonant_liquids_20261116/`,
`contrib/qubit_resonance_hierarchy_20261116/`, and
`contrib/small_d_mobile_matching_20261116/`; their existence alone does not
indicate completion.

The frozen proof was read in full by this coordinating chat. Its common
interacting bulk reference, bounded-port cutting estimate and absorption
of the finite-ring remainder address the previously identified extensive
error problem. This is a mathematical review of that argument, not a
rerun of its literal-bond implementation or a new CFT theorem.

## The exceptional bands have a ferromagnetic threshold

Use the notation and complete leading operators of
[DILUTE_SIXTH.md](../qubit_cancellation_classification_20261116/DILUTE_SIXTH.md).
Remove the stated scalar, cancel the nearest-neighbour hop with the
available sixth-order counterterm, and write the remaining chemical term
as `-mu6 N`. We inspect the particular value `mu6=-2t`.

For q=1, let m=4, t=1, and forbid separations 1,2,3. Then

```
H4 = -(T4+T4†) - 2D4 + (10/3)D5 + (14/3)D6 + 2N.
```

For q=3, let m=3, t=128/3, and forbid separations 1,2,4. Then

```
H3 = -t(T3+T3†) - 2t D3 + 8G8 + 2t N,
G8 = sum_j n_j n_(j+8)(1-n_(j+3))(1-n_(j+5)).
```

Every operator is compressed to its entire stated constrained space.
Throughout this section take a periodic ring with L divisible by m and
large enough to distinguish the local distances. These formulas concern
the sixth-order effective models, not the full finite-epsilon qubit bond.

**[Proved from the displayed effective operators.]** Let A be the
configuration adjacency matrix for allowed hops by plus/minus m, counting
their multiplicities. Let L_graph=diag(degree)-A. A hop cannot pass another
particle because the minimal particle spacing is m. The cyclic ordered
word of particle residues modulo m is conserved.

For m=4 a hop toward the next particle is forbidden exactly when the
consecutive gap is 4,5,6 or 7. For m=3 it is forbidden exactly for gaps
3,5 or 7; a gap of four is absent from the space already. Increasing the
gap on the other side by m cannot violate a constraint. Counting both
orientations gives the identities

```
H4 = L_graph + (16/3)D5 + (20/3)D6 + 2D7,
H3 = t L_graph + 2t(D5+D7) + 8G8.
```

The relevant short-distance D terms count consecutive gaps: no allowed
intermediate particle can fit in them. Both Hamiltonians are therefore
positive graph Laplacians plus nonnegative diagonal potentials.

In a single residue class, the extra potentials vanish. On the resulting
coarse ring of M=L/m sites, the exact operator is

```
H_single = t sum_j (I - SWAP_(j,j+1)).
```

Equivalently it is the hard-core t-V chain at V=-2t with the chemical
term chosen above. The kernel consists of the completely symmetric
spin states, one for every particle number 0,...,M. Its one-particle
dispersion is `2t[1-cos(k)]`, hence quadratic at its minimum. This is
the ferromagnetic threshold, not an m-component conformal liquid.

## Kernel of the full effective model

**[Proved.]** Each connected configuration component contributes a zero
mode of the graph Laplacian only as a constant-amplitude vector. If any
vertex in that component has positive diagonal potential, that component
contributes no zero mode of H_m.

Consider a nonempty ordered residue word with a change between two
neighboring particle labels. Move the earlier particle toward the next
one in steps of m. Every step before reaching the minimal admissible gap
is allowed: the gap behind increases, and all other gaps are unchanged.
For m=4 the final gap is 5,6 or 7; for m=3 it is 5 or 7. The reached
configuration has strictly positive potential. Thus no mixed-residue
component has a zero mode. This argument includes jammed components,
which already have a positive-potential gap.

The remaining components have a single residue label. At fixed particle
number their coarse nearest-neighbour hard-core hopping graph is
connected, including the empty and fully occupied limits. There are m
such components for each N=1,...,M, and one common vacuum. Consequently

```
dim ker H_m = m M + 1 = L + 1.
```

This matches the kind of enlarged compatible-ring ground space suggested
by the parent transfer resonance. It is a result about the two explicitly
displayed leading models only. It supplies neither an all-m identity nor
a theorem for their finite-epsilon microscopic continuations.

The original three-/four-minimum observation is correct as a band
statement. It does not supply independent positive-velocity modes: the
leading models above instead have conserved ordered labels, ferromagnetic
single-residue sectors and quadratic excitations. Other chemical
potentials, admissible microscopic detunings and higher-order residue
motion require separate analysis. They are the counterexample chat's
follow-up task, not exclusions established here.

## Independent checks

`check_residue_laplacians.py` constructs the complete constrained occupation
spaces directly, without importing the sibling implementation. It checks
the identities after multiplying by three, so all matrix and graph tests
use exact integer arithmetic. It separately constructs the coarse SWAP
operator in its full spin basis and checks every residue embedding.
Connected-component counting verifies L+1 zero modes and conservation of
the complete cyclic ordered residue word.

Run from the project root with NumPy and SciPy:

```sh
OPENBLAS_NUM_THREADS=1 python3 -B contrib/research_coordination_20261116/check_residue_laplacians.py
```

The checks cover m=4 at L=12,16,20 and m=3 at L=9,12,15,18. The saved
`residue_laplacian_checks.json` records all dimensions and counts. These
are independent checks of the published effective formulas, not of the
microscopic Schur calculation and not measurements of conformal c.
