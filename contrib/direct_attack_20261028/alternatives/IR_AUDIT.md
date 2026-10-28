# Independent audit of finite-cell auxiliary matching

**Verdict: checked under the stated finite-cell and published determinant-expansion hypotheses.** I found no remaining fatal error in the revised [FINITE_CELL_MATCHING.md](../infrared/FINITE_CELL_MATCHING.md). Its conclusion is a summable auxiliary-induced logical action with an error proportional to the leakage interaction. It is not a proof of the gapless RG, a thermal central charge, or the all-length mobile construction.

This review covered [RESULT.md](../infrared/RESULT.md), the earlier [spectral embedding](../../chat_followup/MASSIVE_AUXILIARY.md), and the revised matching proof. I independently opened [Pedra–Salmhofer, Theorem 4.5 and its norm definitions](https://www.wias-berlin.de/people/koenig/www/FG/PedraSalm.pdf#page=18). The review is algebraic and analytical; there is no claim of formal verification or human refereeing.

## Points that pass

1. **Finite embedding.** A dedicated auxiliary occupation state for each excited eigenvector gives an even isometry and exactly reproduces its energy. Each auxiliary energy is at least half the original cell gap. All logical states with auxiliary vacuum are valid, so every invalid state has positive auxiliary occupation. Extending each complete-support interaction by its isometry makes the valid-state projections commute with the enlarged Hamiltonian. The cut-interaction energy and heat-trace comparisons remain valid for summable supports with their stated incidence norm.

2. **Coefficient conversion.** The bound \((8h^2)^N\|O\|\) follows from the number of occupation-basis matrix units, the vacuum-projector expansion, and the degree ceiling \(2N\). It is deliberately large but valid. It bounds a sum over every species and monomial, which is stronger than a supremum over one pinned field label. The required spatial support weight may therefore be enormous. It is still finite for a fixed finite cell and must be chosen before the final small coupling.

3. **Massive covariance.** With auxiliary energies bounded below by \(g>0\), each fixed row or column has weighted circular-time integral at most \((g-\kappa)^{-1}\). Species diagonality avoids an extra species sum in this integral. The revised \(\alpha=2/(g-\kappa)\) safely covers the source's larger time-difference convention. The determinant bound is independent of the logical bandwidth. Logical covariance is zero during this integration, and external auxiliary fields are set to zero afterwards.

4. **Spatial and temporal weights.** In a connected expansion, on-site auxiliary propagators join vertex supports at a common cell. The output spatial span is bounded by the sum of the input support lengths. The external time diameter is bounded by the spanning-tree time lengths, including internal point-split diameters. Thus weights can be assigned to input vertices and absolute tree edges without altering the determinant estimate for the remaining contractions.

5. **Linear dependence on leakage.** If the input norms obey the stated smallness condition, the auxiliary effective-action map is analytic on a ball around the block-diagonal reference. The complex Cauchy disk in Section 5 remains inside that ball. Its radius and bound give \(4r\), with no logical propagator and no inverse slow bandwidth. A purely logical action factors out in the Grassmann integral. The revised derivation uses a flat logical Gaussian reference, so it does not invoke Wick's theorem for an interacting logical Hamiltonian.

6. **Vacuum cancellation.** The required condition is the two-sided operator identity \(K_D P_{\rm vac}=P_{\rm vac}K_D=0\), not just zero vacuum expectation. The local block-diagonal hypothesis supplies the stronger identity. Every surviving normal monomial then has both an auxiliary annihilator and creator. The rightmost vertex kills the auxiliary vacuum in every nonempty zero-temperature time-ordered product. This justifies the vanishing baseline action with arbitrary external logical fields. It is not a finite-temperature cancellation.

## Corrections made during the review

The initial draft had three norm/regularization gaps, now repaired by its author:

- A scalar perturbation has zero nonconstant kernel seminorm. Its scalar density is now retained separately as \(r_0\), including in the pressure bound.
- A point splitting whose span grows with vertex degree is unsuitable for unbounded quasi-local support sizes at fixed norm. The revised prescription puts all fields of each vertex inside an interval of length \(\epsilon\), independently of degree. Creation fields are later than annihilation fields, giving the normal-order equal-time convention.
- A pinned-field norm alone does not bound a new locally marked density after summing species. For example, pinning a field in \(\sum_s a_s^\dagger a_s\) fixes \(s\), whereas pinning a new marker sums over \(s\). The revised input norm also bounds the complete anchored density coefficient sum. This makes the local source/marker estimate valid.

The trivial case of zero auxiliary modes is also now separated, avoiding division by a zero determinant constant.

## Limits to retain in subsequent claims

The zero-temperature baseline cancellation is not a bound on all finite-temperature source kernels. The separate block-diagonal pressure comparison does not supply such a bound. A derivation of a thermal central-charge coefficient must keep this distinction.

The finite coefficient conversion is not uniform in a growing length cutoff. The auxiliary vacuum projector alone has coefficient norm \((1+h^2)^m\). Since the mobile construction has an unbounded cell-length space and marker motion crosses a finite cutoff, this theorem cannot yet be applied to the complete mobile chain by taking larger cutoffs.

Finally, the output contains general retarded, high-degree logical interactions. The current constructive-RG source theorems do not automatically cover that entire input class or identify its Virasoro central charge. Those are further mathematical obligations, even after the finite-cell massive matching step is accepted.

Reviewed 2026-10-28 by a separate GPT-6-based Codex agent after exchanging concrete corrections with the author. The conclusions above refer to the revised source with the separate scalar density, strengthened density norm, degree-independent point splitting, and flat logical Gaussian reference.
