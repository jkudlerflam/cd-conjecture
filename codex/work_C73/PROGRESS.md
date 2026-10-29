# C73 progress

2026-10-28.

[proved work record] Read the project memory, S2 and S backgrounds, and the complete C59, C62, C65 and C71 materials, including scripts and saved checks. Their trace construction, parameter matching and second-order bubble bounds are retained as inputs. All writes are confined to work_C73.

[assumed target] Prove TP with cancellations performed before coefficient norms, or identify an actual coefficient obstruction. A divergent graph majorant alone is not such an obstruction.

[heuristic under investigation] The chiral contractions of density vertices decompose into paths and cycles. Two-cycles require C71's distributional bubble preparation. Higher cycles have different local integrability and cancellation properties. The full chiral loop identity must be distinguished from its high-covariance version, where low-cutoff insertion terms remain. I am checking both the fixed-order ultraviolet power counting and the factorial cost of the resulting organization before making an analyticity claim.

[assumed remaining] Even an analytic fixed-scale ultraviolet preparation would still require the infrared joint defect/thermal comparison and compatible vacuum derivatives in C59 TP. No all-order claim has yet been established.

[proved partial result] PROOF Sections 2-4 now establish an analytic integration of one high chirality, uniformly in the sea cutoff. A generalized Young inequality bounds every cycle of length at least three with one fixed moment reserve. C71 supplies the exceptional two-cycle. Open chains have finite moments and periodize exactly; cycle windings give arbitrary fixed thermal power with only a polynomial cost in cycle length.

[proved partial result] PROOF Section 6 proves that every fixed-order zero-temperature coefficient of the full reference high integration stabilizes once Lambda exceeds max(K,K_e)+(n+1)kappa+1. The proof integrates closed-loop frequencies before taking norms and uses bounded momentum transfer along open chains. It includes vacuum coefficients. The finite-order bound does not control growth with perturbative order.

[proved caution] The high-covariance Ward identity contains low-cutoff insertion terms. Odd cycles cancel, but an explicit high four-cycle is nonzero. It is finite and thermally controlled, so it is not a coefficient counterexample to TP. The complete two-chirality analytic bound and the infrared joint comparison remain unresolved.

[numerical, complete] The initial cycle image check failed with residual 1.30e-11 at a 2e-12 tolerance because its finite image range was too short. Increasing the range without changing tolerance gives a final maximum residual 1.08e-16. The independent frequency contour and cutoff plateau checks pass. checks.json records the final controls.

[proved, final write-up] RESULT.md and STATE.md record the partial theorems and their boundaries. PROOF distinguishes infinite-spatial-volume coefficient bounds from the full regulator-uniform TP requirement, and fixed-order stabilization from an analytic all-order sum. SOURCE_AUDIT records the primary-source context. No complete TP proof or failed TP coefficient is claimed.
