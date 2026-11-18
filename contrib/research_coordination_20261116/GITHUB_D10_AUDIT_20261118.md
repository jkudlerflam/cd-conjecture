# GitHub source check and the next smallest-d calculations

Date: 2026-11-18. OpenAI Codex, GPT-6 family.
Status: source audit, review of completed contributions and follow-up
assignments. No new Hamiltonian or independent proof verification is claimed.

The user clarified that the reported d=10 unbounded-c example was probably
on GitHub. This extends the [previous source audit](SMALLEST_D_FOLLOWUP_20261117.md)
to the repository's live branches, issues, comments and relevant history.

## Search scope and result

**No matching d=10 arbitrary-c construction was identified.** This remains
a search result, not a nonexistence theorem or a conclusion that the
collaborator's report is wrong.

The GitHub API was read directly after synchronizing the local checkout:

- `main` was at `b80b1c35f7b6443d5f5a29dbc374abb5ac752b1e`.
- The only other returned branch was `claude/loving-franklin-kj8ay2`,
  at `10980ccbaf985c6a2f8ee0cee135046c308db29c`.
- The repository returned five issues (#1--5) and one pull request (#6).
  All descriptions and all four issue comments were searched. The PR
  returned no inline review comments.
- Discussions were disabled. No forks were returned to the current
  account; this does not exclude unavailable or unrelated repositories.
- Targeted searches of the fetched Git history covered ten-state wording,
  d=10, unboundedness and the history of the identified demonstrators.
- Additional indexed GitHub code queries used central charge with
  unboundedness, d=10 and ten-state wording. The inspected results did
  not identify the reported construction. Search indexing and the chosen
  wording limit any negative inference.

No issue, comment, PR or external message was posted. The live API
inventory is a snapshot, not a permanent assertion about future GitHub
contents. Raw API responses stayed outside the repository.

## The identifiable ten-state references

1. **[The route-A demonstrator](../routeA_demonstrator/README.md).** Its
   physical alphabet is two port states plus eight block states. It has
   one logical flavour, with numerical c approximately one in the stated
   regime. Its original commit is
   [fdc2e7e7](https://github.com/jkudlerflam/cd-conjecture/commit/fdc2e7e7).
   That record explicitly describes c=1, not arbitrary c.
2. **[The mobile two-flavour toy](../mobile_cells/README.md).** Its ten
   states support two payload flavours and a proposed total c=3. Even
   that microscopic c is not established by the recorded tests. Replacing
   the two payload bits by n bits gives d=2+2^(n+1), not fixed d=10.
3. **[C53's exact completion](../../codex/work_C53/RESULT.md).** This
   retains ten states but uses four-site physical interactions and one
   logical flavour. Its exact free c=1 member does not meet the requested
   nearest-neighbour arbitrary-c specification.

[GitHub issue #2](https://github.com/jkudlerflam/cd-conjecture/issues/2)
and its comments concern matching for the ten-letter demonstrator. They
explicitly separate that finite-cell problem from the mobile wrapper's
unbounded cell lengths and from the final CFT. An unbounded cell length
is not itself an unbounded central charge. We have not identified this
issue as the collaborator's intended reference.

The more detailed [state and scope ledger](../minimal_d_unbounded_audit_20261117/SOURCE_LEDGER.md)
counts every marker, carrier and vacancy and records the known d=9--14
proposals. None should be relabeled as a verified d=10 unbounded-c model.

## What the completed smaller-d round actually added

**[Finite algebra, numerical matching and perturbative RG interpretation.]**
The [nine-state resonant seed](../ten_state_resonance_audit_20261117/README.md)
has a two-flavour leading symbol with four simple Fermi points. Its
generated pair-transfer coupling is nonzero, but its leading weak-coupling
flow lies in the marginally irrelevant basin at the audited point.
Thus the strict non-resonance certificate is sufficient, not necessary
at this perturbative level. The selected-sector prediction is c=2,
well below d-1=8. It is neither a microscopic CFT theorem nor an
arbitrary-flavour construction. A one-vacancy d=10 lift is not automatic.

**[Proved in the contribution; not independently reviewed here.]**
The [eleven-state transfer](../minimal_d_unbounded_audit_20261117/README.md)
retains all physical states while extending the biased clock code,
zero-port full length selection and leading simple-band arguments to
every finite n. The tree still has exactly coincident opposite-momentum
pair sums for n>=2. Its non-Cooper margin does not supply the missing
full pair-sum margin. The actual C84 loop may split those resonances,
but its all-n nonidentity and margin are open. Ground-sector selection
at nonzero port, volume-uniform matching and interacting infrared
stability also remain open.

**[Saved-data numerical audit.]** The [qubit interpretation study](../v4_interpretation_audit_20261117/README.md)
shows that the provisional 2.5669 entropy slope has a large contribution
from parity-branch selection and another from interrupted accuracy
refinement. The density oscillation is compatible with interference
between the same two valleys. This supports no c>2 claim and settles
neither c=2 nor a possible eventual c=1 or gap. No DMRG was restarted.

## Follow-ups delivered to the same three chats

1. **Nine physical states:** implement the actual five-control extension
   proposed in the seed's LIMITS.md at three logical flavours. Seek six
   simple Fermi roots and evaluate the complete coupled resonant
   interaction flow. Do not replace the physically available tridiagonal
   onsite family by an arbitrary logical matrix. A three-flavour success
   is a scaling test, not an all-n result.
   Output: `contrib/nine_state_three_flavours_20261118/`.
2. **Eleven physical states:** calculate the first nonzero splitting of
   opposite-Fermi-point pair sums produced by the genuine C84 loop, with
   its inactive endpoint and physical control restrictions retained.
   Seek an all-n nonidentity/margin argument or an exact obstruction;
   compare errors with the new splitting scale without adding states.
   Output: `contrib/eleven_state_loop_splitting_20261118/`.
3. **Qubits:** derive the exact symmetry selection rules and conditional
   scaling-dimension criteria for operators that can lock the two valley
   modes. Reuse existing field theory and data, preserving the distinction
   between valley centres and Fermi endpoints. No long numerical
   calculation or interrupted DMRG run is authorized by this assignment.
   Output: `contrib/v4_locking_operators_20261118/`.

These are research assignments, not claims of their outcomes. The d=20/22
extensions remain deferred. All source and inference qualifications above
remain even if the next finite calculations pass.
