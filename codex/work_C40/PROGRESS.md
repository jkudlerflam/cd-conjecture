# C40 progress

[proved, scope] Read the project memory, route-A background, STATUS, explicit candidate, C26, C31, C36, infrared-control report, and C8 sections 40 onward. C32 and C20 are being completed in bounded reads after earlier output truncation. C33 and C37 progress were read. C38 and C39 progress were absent at the initial read. All writes stay in work_C40; project STATE is outside the authorized scope.

[proved, inherited distinction] C33 reports a rational certificate for 46 simple Fermi roots and finite quadratic/quartic margins at n=23, including a finite conversion bound. C37 is auditing pair-product certificates through n=23. No all-n nonresonance theorem is supplied by these reports. C40 will not duplicate their finite-band calculations.

[proved, scope of available locality] C8 sections 46,48 give fixed cell coordinates and a volume-independent one-step local dressing. Finite-order connected coefficients use only massive denominators. C8 sections 49,51,52 and C36 rule out replacing an irreducible low-frequency action by finite-time control, bare-source correlations, or a global bare-Q inverse. Their remaining connected-action assumption is not proved here.

[in progress] Derive a conservative multiscale RG budget incorporating the slowest velocity, band formation scales, finite momentum cutoffs, and collective interaction norms. Derive fixed-order connected spatial/time bounds explicitly, and formulate the precise additional low-frequency remainder needed for fixed nonzero coupling. Track conversion and port orders separately, including the distinction between C31 and the 23-letter bond.

## Hierarchical matching and coordination

[proved, scope] Completed the requested C32/C20 reads and the omitted portions of earlier bounded reads. C38 and C39 progress now exist and were read. C38 treats scale-weighted quadratic mixing and a lacunary momentum target; C39 treats an auxiliary ruler. C40 confines its new work to connected finite-order elimination, unequal-scale RG and the infrared acceptance conditions.

[proved, finite-order strategy] On a connected interval of at most r+1 cells, a Riesz-projector/direct-rotation expansion at radius proportional to Delta/[(r+1)J] gives explicit finite coefficient bounds. Inclusion-exclusion over interval bonds removes disconnected terms. The resulting bound grows at least like the majorant (C r)^r; it does not imply a convergent thermodynamic series.

[heuristic, RG result being written] Use the full free lattice propagator for each flavor until its individual linearization scale, retaining every slower band. A final cutoff below both the smallest linearization scale and v_min times every low-degree momentum mismatch gives only Gaussian marginal terms. Finite matching constants may contain powers of bandwidth ratios; the marginal running interval also contains their logarithms. Smallness is measured against the slowest mode and the full interval, not merely against the largest hopping.

[proved, important distinction] Reducing lambda suppresses higher port orders but does not reduce the ratio of the lambda^2 g^4 non-Gaussian terms to lambda^2 g^2 hopping. Conversion g (tau in C31) must be chosen inside the RG tolerance separately. The mobile length-dependent scalar at order lambda^2 also requires separate equalization, before lambda is chosen.

## Completed derivation and exact checks

[proved] PROOFS.md now contains an explicit connected coefficient majorant C_r(mu) J^r/Delta^(r-1), obtained from finite-cluster Riesz projections and canonical direct rotations followed by bond inclusion-exclusion. The support is at most r+1 cells. Fast time moments are bounded at every fixed order by reduced-resolvent integrals. No convergent all-order expansion is claimed.

[proved] A conservative legal-cell gap for the original bond follows from the excited token graph row sum: Delta_mass>=2+|g|+gamma ensures a gap at least gamma at every legal cell length when the rail hoppings are one and t_d<=1. This does not replace the missing illegal-sector comparison.

[proved by exact computation] checks.py passes exact neutral-vertex enumeration for n=1 through 6, up to 3102 oriented vectors. It also checks a two-component hierarchy/positivity example, the rational Gaussian threshold, and the integrated RG majorant. No floating-point signs, microscopic central-charge measurements, or numerical extrapolations enter these checks.

[heuristic, conditional infrared result] The sufficient budget includes full lattice-band matching constants, the collective marginal norm including forward terms, the logarithm down to the slowest momentum-detuning energy, relative Gaussian stiffness, and dressed filling margins. M2 states the remaining irreducible low-frequency action assumption without identifying it with a bare-source correlation or global bare-Q inverse.

## Final handoff

[proved, completed] RESULT.md, PROOFS.md and local STATE.md are written. The Python syntax and prose checks pass. The exact integer/rational checks pass. An audit corrected the isolated-root projector estimate and included forward marginal couplings and possible linear irrelevant feedback in the RG majorant.

[proved, latest coordination] Read completed C33 and C37 results and C38's result. C37 supplies certified frozen targets for n=1,...,23, not the C32 mobile core or an all-n result. C38 independently confirms that a crossing-only loop/conversion schedule does not close the nonresonance budget. C39's last progress includes a proposed all-n two-bus/biased-vacuum derivation whose local construction is still being audited; C40 does not claim or duplicate that result.

[heuristic, final assessment] The conditional fixed-n infrared argument is complete under explicitly stated F, B, M2 and R. The new connected estimates prove every fixed perturbative order but leave the exact irreducible low-frequency remainder M2 open. The original 23/28-letter all-n nonresonance premise remains open. No microscopic central charge or unconditional all-n counterexample is claimed. Only work_C40 was modified.
