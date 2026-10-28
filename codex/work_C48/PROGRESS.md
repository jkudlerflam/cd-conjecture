# C48 progress

[proved, scope] Work is confined to work_C48, with no network, subagents or large computations. The requested project memory, background, lead claim, formal S, C39/C42/C44 finite inputs, C45 reduction, C46 obstructions, literature summary and C47 progress have been read. No C49-C52 progress files were present at the initial check.

[proved, existing work] The newer contrib/direct_attack_20261028/infrared/INFINITE_CELL_PROGRESS.md already proposes an operator-valued coarse-time polymer theorem, with a separate analytical audit. C48 will attribute that route and supply a self-contained conservative bound and explicit substitution of the C44/C46 constants, rather than claim discovery of the representation.

[proved, new parameter simplification] The legal rooted heat bound has C0 depending only on n and d0=4464, with Delta no larger than the wrong-length slope gamma. Thus a time block tau=T_n/Delta can satisfy both the heat condition and the polymer condition with T_n independent of wrapper ratios. The resulting exponential rate is a positive n-dependent multiple of Delta, and the Lipschitz matching constant can be dimensionless and independent of Delta.

[assumed, application boundary] A bounded, locally fermionized interaction in the well-formed rooted tensor product is required. The infinite-cell embedding enlarges the thermal Hilbert space. Its invalid states have a positive threshold and exponentially small pressure contribution, but this is not an exact equality of finite-temperature source functionals for the original physical fixed-length ensemble. The all-beta wording also requires qualification when the rooted cell trace exists only above a positive beta threshold. RG and the CFT limit remain separate.

[heuristic, work in progress] Write the block/polymer proof with explicit constants, a uniform low-temperature theorem, a leakage-linear estimate after C45, and the exact resolution of the C46 toy model. Independently check finite-dimensional semigroup/embedding identities and the numerical margins in the chosen cluster constants.

[proved, completed derivation] PROOF.md now gives the operator-valued block construction, thermal Schatten bound, arbitrary-interval activity sum, hard-core tree identity and rooted recursion, analytic leakage bound, source qualifications, complete gap and heat substitution, and explicit fixed-lambda schedule. The dimensionless A_dec,n is independent of wrapper ratios, with kappa a positive n-dependent multiple of Delta. RESULT.md and STATE.md state the remaining physical and RG obligations.

[numerical, checks complete] checks.py passes. The integrated C46 spectral kernel has error 8.7e-15; 729 periodic polymer histories factor to 6.8e-21; the full noncommuting partition sum agrees to 6.7e-16. Spectral embedding and pressure checks pass on dimensions 8 and 20. The conservative cluster ratio is below 1e-39. These are controls, not a microscopic CFT test.

[proved, final scope] Full S is not established. The theorem is for the enlarged rooted model, uniformly in cell count and beta>=2T_n/Delta. Low-energy and pressure comparisons remove its artificial invalid states in those senses only. Physical length/source/quotient transfer, complete illegal-sector control and constructive gapless RG remain open. All writes stayed in work_C48.
