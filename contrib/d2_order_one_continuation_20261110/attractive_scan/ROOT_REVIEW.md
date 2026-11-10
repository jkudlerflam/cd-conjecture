# Skeptical read-only review of the R=8.75 bridge

OpenAI Codex, GPT-6 family, 2026-11-10. This review reads the root README,
analyze.py, audit.py, relevant raw records, preservation manifests, variance
and Renyi summaries. It performs no new optimization and edits no file
outside attractive_scan. The largest-size chi 256 continuation is still
pending at review time; final claims must use its completed result.

## Model and selection checks

No model-geometry or normalization defect was found. The displayed model
has Jx=.5, Jy=-.5, Jz=.625, Dz=.1, hz=-.875. It has six expanded Pauli
strings, four grouped terms, four nonzero magnitudes, literal spread 8.75,
exchange-only spread 6.25, and singular values (.625,.6,.4). Its exact
Jx+Jy=0 cancellation and the smaller DM scale are correctly disclosed.
The shared wrapper applies the physical field once on every site, including
the ends. Its staggered computational gauge is a product basis change,
not blocking or a change in microscopic geometry. All inspected selected
records use the unchanged solver SHA256
8b55ff8a0c8ada4d3da36cc208a92566eb0459c045ea93b9facdc74a2980e847.

Every currently selected row in size_comparisons.json passes its sweep gate.
The independent 48 seed and 64 intermediate size agree with the stated smooth
entropy trend. The literal model, entropy conventions and small periodic
comparison use the same physical parameters. No missing-parity substitution
was found in the currently cited 24/48/64/96 comparisons.

Two latent guard weaknesses were identified and corrected by the parent:

- Cross-seed ground-branch selection could admit a lower-energy unconverged
  record even though within-seed version ranking preferred convergence.
  analyze.py now filters sweep-converged rows before cross-seed/sector
  selection and labels pair convergence explicitly.
- audit.py previously exited successfully when failed_sweep_gates was
  nonempty but its generic errors list was empty. It now exits nonzero
  for either condition. Current numerical selections did not change.

For a final quiescent audit, no checkpoint check should remain marked
'retry after completion'. Such race handling is appropriate during live
updates but is not a final integrity pass. Historical snapshots must use
their preservation manifests: a copied JSON can legitimately retain its
original live-checkpoint path, which may now name a later state.

The Renyi slope builder has an analogous provisional-selection issue:
its caches originally omitted sweep convergence and its slope selection
could run with only one completed parity at a new size. The currently
quoted 48-to-96 values agree with valid selected records. Final 192 Renyi
claims should require both completed converged sectors or explicitly carry
provisional metadata. This issue was reported to the parent.

The README's common checkpoint-root sentence should point to exact paths
inside each raw record: the bridge and attractive branches use different
/private/tmp directories. The final preferred-attractive paragraph also
needs the now-completed 192 crossover outcome instead of 'under examination'.

## Physics assessment

The bridge is more credible than the attractive entropy flags: its parity
weight changes contribute little to the excess, transverse order and
density modulation decrease across the tested sizes, the parity splitting
is approximately inverse length, and an independent intermediate size
reproduces the same entropy growth. These facts support further pursuit,
but do not turn a finite-size coefficient into an asymptotic central charge.

A cheap post-hoc energy consistency check is encouraging. Using selected
energies at 24/48/96 and the newly completed 192/P0 chi 192 energy
-129.9212708835584 gives
D(24)=-.00176176447 and D(48)=-.000858745594, where
D(L)=E(4L)-3E(2L)+2E(L) cancels bulk and boundary constants.
The effective correction power -log2[D(48)/D(24)] is 1.03672, versus .94970
for exact XX on these sizes. This is compatible with a leading 1/L Casimir
correction, unlike the attractive reference's .38987. It is post-hoc,
assumes the selected finite-size branch, and does not independently give
c without velocity information. It must be refreshed after chi 256.

The main unresolved alternatives remain material:

- A c approximately 1 bulk theory with a long boundary or dilute-sector
  crossover. Small periodic rings currently look near 1 and do not identify
  a second bulk mode; an inverse-length parity splitting alone can also
  involve boundary states.
- A nonrelativistic or multicritical crossover rather than a unitary CFT.
  No conformal spectrum or independent velocities have been established.
- An ordered sector whose normalized order has not yet reached its limiting
  value. Binary Schmidt-parity subtraction rejects one narrow cat mechanism,
  not all symmetry-breaking cats; the attractive Z-order reserve illustrates
  that distinction directly.
- Remaining entropy finite-size corrections: Renyi alpha 2/3 estimates are
  lower than von Neumann, and half-filled free controls do not universally
  bound interacting or incommensurate corrections at density about .266.

The README largely states these limitations appropriately. Subject to the
finished largest-size chi check, both parity states, and the final strict
audit, 'promising finite-size candidate' is defensible. An established
c>1 unitary CFT or counterexample is not supported by these checks alone.

## Final resolution — 2026-11-10

The selected L192 even state has now completed chi256, and both L192
parities are completed and sweep-converged. The final audit reports 126
completed record versions and two explicitly interrupted explorations,
with no failed sweep gates, pending checkpoint checks, retry markers, or
errors. The cross-seed, audit-exit and Renyi selection guards are resolved;
all reported main Renyi pairs explicitly have converged selections and
both converged parities. The final README also resolves the checkpoint-path
and attractive-route status corrections. No final packaging blocker was found.

The final 96-to-192 von Neumann slopes are 1.44746 (central three) and
1.45303 (fractional window). The corresponding fractional Renyi alpha2
and alpha3 estimates are 1.48087 and 1.37017. Thus the earlier statement
that both higher-alpha estimates are below von Neumann applies to the
48-to-96 data, not the final pair. A common asymptotic coefficient remains
unestablished. The largest-size chi192-to256 entropy change is about
0.000247, much smaller than the observed excess.

The updated energy check uses E192=-129.92127124133077 and gives an
effective inverse-size correction power of 1.03612; the three-smaller-size
inverse-length fit predicts E192 within 2.18e-5. This supports consistency
with inverse-size energy corrections and does not independently determine
c, velocities, or a dynamical exponent.

The final evidence supports the README's restrained description,
'promising finite-size candidate'. Boundary or c=1 crossover behavior,
nonrelativistic criticality, an eventual order plateau, and remaining
Renyi/size corrections remain unresolved; small periodic rings still do
not identify a second bulk mode. No established c>1 unitary CFT or
counterexample is claimed. This addendum uses completed stored artifacts
only; no new numerical calculation was run.
