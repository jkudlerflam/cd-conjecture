# p153c connected-amplitude measurement protocol

Prepared 2026-11-18. This protocol is **not launched**. It is specific to
the fixed asymmetric p153c Hamiltonian in README.md and does not modify
the historical sym40 registration. Auxiliary derivatives are extrapolated
to zero source at fixed physical coefficients; no L-dependent field tuning
is permitted. The one-worker finite checks in this folder are complete.

## Saved-record audit

[source_audit.json](source_audit.json) hashes the actual L=512, parity-even
C195 p153c records at maximum bond dimensions 1536 and 2048, their model
and measurement code, and the companion medium interface as read. Both
records have a 48-site one-body matrix beginning at physical site 232.
The stored `CP` is the connected **neutral bond** correlator of
`P_j=XX+YY`, anchored at site 128. It is not a separation-two creation
pair or a dynamic quartet cumulant. Its long distance range cannot replace
the missing time integrations or external propagators.

Neither record stores the six-insertion response, dynamic single-leg
kernel, a source derivative of it, or the fused vertex normalization.
Interacting higher correlations cannot be reconstructed by imposing Wick
factorization on these equal-time arrays. The companion analysis audits
additional long-pair records, but they likewise do not determine the
amputated source-squared kernel. Their apparent carrier shift is useful
for designing windows, not a certified thermodynamic cutoff.

The Pauli coefficients have spread 3198.80418 after constants,
telescoping terms and a common positive scale are removed (none of the
first two occur in this displayed basis). There are five independent
real controls before scale removal, four after, and six literal monomials.
The small effective t and possible small velocities are separate from
this microscopic hierarchy. No simplicity improvement is claimed.

## Minimum amplitude data and conventions

1. **Reference and projection.** Establish the correct thermodynamic
   state branch and both dressed modes at this fixed point, with the
   companion response protocol. Fix a0, the momentum patches, the
   ultraviolet regulator, dynamic/static order of limits, and the
   canonical compactification. Do not assume that the smallest finite
   leading-model grand sector is the physical large-chain branch.
2. **Single-field response.** Measure the full normal and anomalous
   frequency/momentum two-point kernels for the selected microscopic or
   SW-dressed probes. Record residues, harmonic admixtures and analytic
   UV contacts. Check conditioning before inverting. At general external
   momenta the full Nambu matrix is required; it cannot be replaced by
   an equal-time occupation or a fitted single exponent.
3. **Connected quartet response.** At the *physical* source strength,
   measure the four-single-leg connected kernel, subtract all three
   two-point products and amputate using the measured kernel from item 2.
   This includes higher source orders and has no need for an exact
   microscopic U(1). Alternatively, in a controlled leading-model
   calculation extract the coefficient of two X dagger insertions and
   four single legs by a six-insertion cumulant or mixed energy response,
   as in (6)--(7) of MATCHING.md. This latter calculation isolates a
   perturbative term; it is not the full physical response at lambda=1.
4. **Local and composite matching.** Vary external residual momenta and
   frequencies inside the patches, compare matching regulators and
   scales, and separate infrared running from the local coefficient.
   Determine the fusion/normal-order map to
   `:exp(2i theta_A+2i theta_B):_(a0)`, including operator mixing. The
   product of measured single-leg residues alone does not fix this map.
   Bound errors from finite size, source differentiation, state branch,
   regulator, omitted source orders and microscopic elimination.
5. **Cross channels.** Retain both B2 and B3 and their interference.
   Keep other near-resonant pair, conversion, density-assisted and
   relative-sector vertices in the same convention; quantify their
   generated contributions and total norm for RG.md (6). Test the
   complete parity-even lattice using covariance error bounds. A vanishing
   measured quartet term alone would not certify c=2.

For an energy-derivative implementation use several decreasing source
strengths, both real and imaginary source phases, and all needed mixed
sign combinations. Compare with a spectral/correlation derivative on a
bounded control. Extrapolate the derivative before claiming a coefficient;
source-induced branch crossings invalidate a single analytic expansion.
Finite ground-energy accuracy must resolve a sixth-order difference, not
merely the total energy. The real-space Rayleigh check in this contribution
calibrates the algebra without relying on such a delicate subtraction.

The necessary momentum accuracy is set by the **lower** mismatch bound
and the running-amplitude budget in RG.md, not by bare band minima.
The separate four-channel covariance/velocity measurement remains required.
No production restart, additional microscopic deformation, or extrapolation
from the separate boson-times-Ising candidate is implied by this protocol.
