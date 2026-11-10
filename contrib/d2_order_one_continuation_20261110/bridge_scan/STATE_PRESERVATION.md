# Checkpoint preservation audit

The first L96 even-sector chi192-to256 continuation preserved its completed
chi192 JSON only. Its chi192 MPS NPZ was overwritten and is unavailable.
The raw JSON and completed pre-upgrade diagnostics (including spectra and
stable before/after hashes) remain, but the old chi192 MPS cannot be loaded
for new diagnostics. Do not claim all old-chi states are replayable.

Before either L48 chi128-to192 continuation began, both chi128 JSON records
and NPZ checkpoints were preserved and verified byte-for-byte. The current
L96 even chi256 and odd chi192 records/checkpoints were preserved at the
same time. `state_preservation.json` lists hashes and exact snapshot paths.
Checkpoint snapshots live under `/private/tmp/d2-order-one-continuation-bridge/validation_snapshots`.

The L192 chi192-to256 validation uses its separately registered preservation
procedure and manifest to save both the raw JSON and NPZ before overwrite.
The raw saved JSON retains its original live checkpoint path for provenance;
use the preservation manifest's checkpoint path when loading an old state.

At completion, `final_state_preservation.json` hash-verifies immutable raw
JSON and NPZ copies of all eight final lead states: both parities at L24,
L48, L96 and L192, using each sector's highest completed chi. Its `preserved`
array has the same per-state schema as `state_preservation.json`. This does
not recover the unavailable L96 even chi192 NPZ; that exception remains
explicit in the final manifest.

The parent also copied the winning-candidate checkpoints to the durable,
Git-ignored workspace directory `../checkpoints_local`. Its
`../local_checkpoint_manifest.json` verifies those copies by SHA256 and
includes independent L64/seed23 checks and the pre-upgrade L192 chi192
state. Prefer that manifest for durable replay after `/private/tmp` cleanup.
