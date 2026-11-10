# Momentum-resolved cross-check — 2026-11-10, before spectrum output

Add exact translation-orbit ED for even-parity rings at L=8,12,16,20,24
in momenta 0 and ±2pi/L, with the identical XX control. Validate the
orbit representation against full-space ED at L8 and explicit Hermiticity.
The lowest spin-±1 excitation supplies a *candidate sound velocity*
`v±(L)=L[E(±2pi/L)-E0]/(2pi)`. Compare its size trend with the periodic
ground-state Casimir coefficient. This extraction is conditional on these
levels being conformal descendants of the same bulk sector. It does not
independently rule out an additional slow sector or establish a full tower.
Raw spectra, residuals, sector sizes and assumptions will be retained.

Operational correction: an initial completed XX chi24 run failed only at
JSON serialization because TeNPy mutated the options dictionary. Its final
checkpoint is retained and the record marked postprocessing_failed. The
runner now copies options before handing them to the engine, and retries
have distinct attempt labels. Earlier source versions are hash-archived.
No incomplete output enters an entropy fit. Inputs are pinned at 7014ef41
while another operation updates the shared Git checkout.
