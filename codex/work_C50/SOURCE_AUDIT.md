# C50 primary-source audit

[read] Sources consulted on 2026-10-28. The MP PDF identifies itself as v3, 15 September 2022. The PSS PDF identifies itself as v1, September 2025. The HTML rendition of PSS displayed a different manuscript date, so equation references below were cross-checked against the PDF. This is a targeted reading of hypotheses, initialization, counterterms, and beta-function statements, not an independent verification of every proof in either paper.

## Mastropietro and Porta

[read] [2108.13708v3](https://arxiv.org/pdf/2108.13708v3), Assumption 2.1(a)-(d), requires an edge-only energy window, smooth transverse nondegenerate edge branches, localized eigenvectors, thermodynamic convergence, and elastic Fermi-momentum scattering. The microscopic model is a two-dimensional cylinder, not this chain. Section 6.2 produces nonlocal, arbitrary-degree effective kernels. Theorem 8.5, equation (8.61), proves asymptotic vanishing of the multichannel reference velocity and interaction beta functions on constant coupling sequences. Section 9 proves that theorem. Proposition 6.2 and Section 10 transfer control to the specified lattice model. Proposition 10.2 compares fermion, vertex and current correlations with the reference theory. Theorem 4.1 concerns edge conductance. None is stated as stability under an arbitrary externally supplied retarded action.

## Porta, Scola and Singh

[read] [2509.08665v1](https://arxiv.org/pdf/2509.08665v1), equation (2.18), uses translation-invariant finite-range hopping and density-density interaction, with orbital-independent potential. Assumption 2.1(i)-(iv) requires finitely many separated smooth nondegenerate transverse branches and elastic quartic scattering. It does not require equal velocities or separated velocity magnitudes. Equations (4.48)-(4.53) add and subtract counterterms and shift Fermi points, without changing the physical Hamiltonian. Equation (4.52) bounds the generated kernels, followed by a field-degree-dependent coupling estimate. Proposition 4.5 controls running couplings uniformly in volume and temperature. Theorem 2.6 is a dynamical-response theorem. Sections 4-6 use a multichannel reference theory. Arbitrary-degree kernels occur internally, but a theorem for arbitrary supplied kernel families is not stated.

## Application ledger

[assumed] The source facts above are used only within their published classes. The following entries are our comparison with the assumed C48/C49 output, not additional claims about either publication.

| Item | Logical-chain status | Consequence |
| :-- | :-- | :-- |
| Finite number of orbitals, translation invariance | Assumed in the logical output | Fits the one-dimensional setting. |
| Isolated transverse free roots | Imported C39/C42 data, plus the length band | Fits the free spectral hypotheses once common patches are chosen. |
| Elastic quartic scattering | Payload certificate alone is insufficient for the enlarged branch list | Use C46 Section 5's full certificate, or prove the symmetry-restricted version using both exact charges. |
| Arbitrarily unequal nonzero velocities | Allowed at fixed reference data | No new flavour or velocity-comparability theorem is required. Constants may depend on inverse minimum speed. |
| Instantaneous density-density microscopic interaction | Generally absent after elimination | The published microscopic theorem is not a direct application. |
| Summed arbitrary-degree retarded kernel norm | Must be specified by the assumed output | Individual degree bounds or an operator norm do not substitute for it. |
| Relevant quadratic counterterms | Necessary in the analysis | They are bookkeeping variables. Their fixed-Hamiltonian construction must survive the extension. |
| Separate branch gauge invariance | Absent microscopically in general | It belongs to the reference marginal theory and must emerge through a comparison estimate. |
| Sources for desired observables | Not implied by an action norm | Needed to identify nonzero leading amplitudes and anomalous exponents. |
| Uniformity in auxiliary gap and wrapper ratios | Not automatic from fixed-model constants | Normalize kernels in retained-bandwidth units first, as in PROOF Section 2. |
| Physical ensemble and central charge | Outside the present RG reduction | C51 addresses the additional implications. |

[proved, scope] The published multichannel reference result can be used unchanged after the marginal theory has been matched. What remains is its stability under the supplied perturbations. Calling the latter an extension does not prove it.
