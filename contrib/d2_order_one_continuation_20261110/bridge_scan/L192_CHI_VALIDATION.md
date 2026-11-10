# Largest-size bond-dimension validation

Registered before any L192 chi256 output. After both L192 chi192 parity
sectors have completed and passed the sweep gates, choose the lower-energy
sector and continue it at chi256 using the unchanged microscopic wrapper.
Keep the exact Hamiltonian and seed17 fixed. Preserve its completed raw
chi192 JSON byte-for-byte and copy its completed MPS checkpoint before
starting the continuation. The preservation manifest records both hashes
and the separate checkpoint snapshot path under `/private/tmp`; the raw
JSON is not edited to redirect its original checkpoint provenance.

This check estimates the remaining finite-chi error in both entropy
windows, energy, density and correlation diagnostics at the largest size.
It does not establish infinite-size convergence or a conformal spectrum.
The independent small periodic results near c=1 and all size drift remain
explicit. Run at most two local single-thread numerical workers; do not
launch field scans before the L192 assessment.
