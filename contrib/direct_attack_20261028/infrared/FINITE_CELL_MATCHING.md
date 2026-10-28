# Finite-cell massive matching after an exactly block-diagonal reference

Status: a theorem for a fixed finite cell, proved using the published determinant expansion stated below. It is **not** an all-length theorem for the C44 mobile wrapper, and not a theorem about the subsequent gapless RG or central charge.

Provenance: Codex `infrared_theorem`, 2026-10-28. Independently checked by the separate Codex `alternative_attack` agent after correcting scalar bookkeeping, source norms and point splitting; see [the audit](../alternatives/IR_AUDIT.md). The parent also checked the original determinant theorem. This is AI review, not human review. The constants are deliberately inefficient.

## 1. Setting and precise output

Fix a finite parity-graded cell K, an even positive operator A, and a code isometry G from the Fock space of n logical complex fermions. Suppose

    ker A = ran G,       A >= Delta (1-G G^dagger), Delta>0.

On a ring of M cells write

    H = sum_j A_j + B + R.

The interactions B and R are even, translation invariant and summable with a sufficiently strong exponential support weight. Every term B_X is block diagonal between the local all-code space P_X and its complement. The code chain H_P=P B P can be interacting and gapless. No assumption about its spectrum is made.

Use the finite spectral embedding defined in
[the existing auxiliary note](../../chat_followup/MASSIVE_AUXILIARY.md): one dedicated auxiliary mode represents an odd excited eigenstate, two represent an even one. All such modes have energies e_s>=g=Delta/2. Let m be their number, F=n+m, and J the even isometry. The enlarged Hamiltonian is

    H_tilde = H_mass + H_log + K_D + R_tilde,
    H_mass = sum_(j,s) e_s a^dagger_(j,s) a_(j,s),
    H_log = the code interaction acting on logical modes in every aux sector,
    K_D = J B J^dagger - H_log,
    R_tilde = J R J^dagger.                              (1)

Here the embeddings of interactions act on their complete supports. K_D is reducing for the auxiliary vacuum on each support and vanishes there:

    K_D,X |logical,Omega_X> = 0,
    <logical,Omega_X| K_D,X = 0.                         (2)

If m=0 there is no auxiliary integration, K_D=0, and the matching statement is the direct bound on the logical residual. In the remainder assume m>=1.

Take 0<kappa<g and a>=0. For nonconstant Grassmann kernels define a pinned-coordinate coefficient norm by taking the supremum over one external field coordinate, integrating all the others, multiplying by H to the total external degree, and summing absolute coefficients. Include the weight

    exp[a diam_x + kappa diam_t],                       (3)

where diam_t uses circular time distance at finite beta and ordinary distance at zero temperature. Diam_x is the smallest connected interval length in cell coordinates; use cyclic distance on rings. For an instantaneous input interaction its time diameter is zero. Scalar kernels will be handled separately.

Choose H>=2 and H'=H+2 delta, where the safe determinant bound delta=2 sqrt(m) may be used. Put

    alpha = 2/(g-kappa),      omega = 2 alpha/delta^2.

Separate the scalar density r0 of R_tilde before defining kernel norms; its action contribution is exactly beta M r0, with the Hamiltonian/action sign convention understood. Let k and r dominate both the input pinned H'-norms and the stronger instantaneous density norm

    sup_j sum_(X containing j) exp(a |X|)
          sum_(monomials I in V_X) |c_I| (H')^(degree I).

Use these bounds for K_D and the nonconstant part of R_tilde respectively. Summing all monomials and species in this stronger norm is needed for the local-marker/source argument in Section 6; a pinned-field norm alone would not control that sum. Formula (8) controls this stronger norm directly. K_D has no scalar coefficient by (2). Assume

    k+r <= 1/(4 omega).                                  (4)

[proved] Integrating only the auxiliary fermions gives an exact effective logical action. In the notation W used below, apply the published effective-action map and then set every *external auxiliary* field to zero, retaining all external logical fields. This restriction cannot increase a coefficient norm. Its change on adding R_tilde has norm

    || W(K_D+R_tilde)-W(K_D) ||_(H;a,kappa) <= 4r,        (5)

uniformly in M and beta. The logical action H_log factors out of this integration exactly; no inverse logical bandwidth occurs in (5).

[proved] In the auxiliary-vacuum, zero-temperature limit,

    W(K_D)=0,                                           (6)

including its scalar part. Consequently the zero-temperature effective logical action is the action of H_log plus a connected, spatially and temporally summable correction bounded by 4r. This is an exact all-order statement inside (4), not truncation at a fixed perturbative order.

The theorem gives convergence and matching, not stability of the logical infrared phase. At finite beta W(K_D) can contain thermal terms. Equation (5) does not erase those terms or identify their precise asymptotics. The block-diagonal thermal comparison in C46 controls their partition-function contribution separately.

## 2. Finite coefficient conversion

[proved] For an operator O on N fermion modes, the sum of absolute normal-monomial coefficients with field weight h>=1 is at most

    (8h^2)^N ||O||.                                     (7)

Expand in 4^N occupation matrix units. Every matrix entry is bounded by ||O||. Write each matrix unit with creation operators, the vacuum projector product_i(1-c_i^dagger c_i), and annihilation operators. Expanding the projector gives at most 2^N terms. CAR reduction leaves degree at most 2N. These counts prove (7).

An interaction on r cells therefore costs at most

    exp[r F log(8h^2)] ||O_X||.                         (8)

For K_D use the sum of the embedded bond norm and its logical extension norm. Thus a physical interaction norm with exponential weight larger than

    a + F log(8(H')^2)

bounds k and r, up to the fixed support/anchoring convention. This is an explicit finite requirement. It can be imposed when taking a finite-order local SW transform: choose its weight first, recompute its constants, and choose its coupling last.

This conversion is extremely expensive, but finite. It is not uniform when the local cell dimension tends to infinity.

## 3. Covariance and the published convergence input

For positive auxiliary energy e and 0<t<beta the absolute thermal covariance is

    |C_e(t)| = exp(-et)/(1+exp(-beta e)).

It is diagonal in both cell and species. The row and column integrals, with circular time-distance weight, are bounded by

    integral_0^beta exp(kappa d_beta(t,0)) |C_e(t)| dt
       <= 1/(e-kappa) <= alpha/2.                       (9)

Indeed d_beta(t,0)<=t and one may extend the positive integral to infinity. Translated rows and columns have the same integral. The chosen alpha has an extra factor two to cover the published decay-norm convention integrating time differences over [-beta,beta]. At zero temperature this reduces to the one-sided exponential covariance.

[read] Pedra and Salmhofer, *Determinant Bounds and the Matsubara UV Problem of Many-Fermion Systems*, [Theorems 2.4 and 4.5](https://www.wias-berlin.de/people/koenig/www/FG/PedraSalm.pdf#page=18), give a determinant bound for thermal fermion covariances and analyticity of the effective action with arbitrary external fields. Their input is a pinned-coordinate field-weighted interaction norm, not an operator norm. The effective-action remainder after order P is bounded by omega^P v^(P+1)/(1-omega v), with omega=2 alpha/delta^2 and a field-weight loss 2 delta.

[proved application of the determinant bound] Take the direct sum of m flat one-particle bands e_s, with normalized Brillouin-zone measure on each band. The total measure is m. Theorem 2.4, or its Gram representation with orthogonal species vectors, permits delta=2 sqrt(m), independent of volume and beta. Extending the covariance by zero on all logical field labels preserves this determinant bound: minors containing a logical row or column vanish. Thus the logical fields can be retained as external variables throughout. The sharper species-resolved constant is unnecessary.

[proved weighted extension of the cited estimate] In the connected tree expansion used for Theorem 4.5, each covariance belonging to the spanning tree is estimated in absolute value and the remaining contractions form a determinant. For instantaneous vertices the external time diameter is at most the sum of circular time distances on the tree: connect any two external vertices by the tree path and use the triangle inequality. The external spatial diameter is at most the sum of the vertex-support interval lengths, since a covariance is on-site. Thus (3) is bounded by the product of the vertex spatial weights and the tree-edge time weights. Replace each tree-edge decay integral by (9); retain the original determinant bound for the remaining contractions. Every other combinatorial step of Theorem 4.5 is unchanged. This proves the same estimate in the weighted norm:

    ||W(V)-W_1(V)||_(H;a,kappa)
         <= omega v^2/(1-omega v),                     (10)
    ||W_1(V)||_(H;a,kappa) <= v.

Here v is the input H'-norm with vertex spatial weights. More generally its order-P remainder is omega^P v^(P+1)/(1-omega v).

The spatial statement uses connected supports; it does not require a bound on the physical system length. For an exponentially decaying interaction, approximate by finite range and use the absolutely summable bound to pass to the limit.

## 4. Equal-time vertices and the exact action identity

[proved] Normal-order every finite operator polynomial first. At a single vertex put all creation and annihilation fields at point-split times consistent with that normal ordering. Put creation fields later than annihilation fields, so an on-vertex occupation contracts to the Fermi factor, not one minus that factor. For a vertex of degree l, fit every split time inside an interval of length epsilon, with spacing at most epsilon/(l+1). This bounds the extra weight by exp(kappa epsilon) uniformly even for quasi-local interactions with unbounded support degrees. Equivalently use regularized delta kernels of total mass one at these times. Their pinned norms converge to the instantaneous coefficient norm, and their extra time-diameter weights tend to one. The covariance is bounded, so each fixed-order integral converges to the normal-ordered, time-ordered Wick expression. Coincidences between distinct integrated vertex times have measure zero. The geometric bound (10), with a fixed reserve in (4), is uniform for sufficiently small splits, so the entire connected series converges to that limit.

For finite M and beta, use H_mass together with a flat, zero logical Hamiltonian as the Gaussian reference. Expand the arbitrary even H_log, K_D and R_tilde in the finite operator Duhamel formula. Wick's rule for the auxiliary factors alone gives exactly the coefficients just constructed, while the logical fields stay external. Purely logical vertices commute with auxiliary integration in the Grassmann algebra and factor as exp(V_log). The connected auxiliary series is W. A chosen quadratic part of H_log can subsequently be moved into the logical Gaussian covariance; no Wick rule for an interacting Hamiltonian was used. At finite volume, this is an identity of the finite-temperature operator expansion, rather than an asymptotic correspondence. Uniform convergence gives the infinite-volume local kernel limit. One need not integrate a logical propagator, or assume the logical Gibbs state has a gap, to obtain W.

Signs depend on whether one writes V=-integral H_int or absorbs the minus sign into its definition. All estimates above concern absolute coefficients and are unchanged.

## 5. Lipschitz estimate and the vacuum cancellation

[proved] In this paragraph R_tilde denotes its nonconstant part after the exact scalar extraction above. Let U_s=K_D+s R_tilde for 0<=s<=1. In the complex disk |z|<=1/(4 omega r), the interaction U_s+z R_tilde has input norm at most

    k+r+1/(4omega) <= 1/(2omega).

By (10), its effective action has output norm at most 1/omega. Cauchy's estimate therefore bounds its derivative in z at zero by 4r. Integrating over s proves (5). The case r=0 is immediate.

This argument gives a bound proportional to the residual r even when k is much larger than the logical bandwidth. Bounding the entire induced action by O(k^2/g) would lose precisely this advantage.

[proved] Because of (2), every normal monomial in K_D contains at least one auxiliary annihilation and at least one auxiliary creation operator. The coefficients of pure-creation, pure-annihilation and auxiliary-free monomials vanish; this follows by acting on the vacuum and by taking the adjoint. In a time-ordered zero-temperature product of K_D vertices, the rightmost vertex annihilates the auxiliary vacuum. Consequently every auxiliary expectation with only K_D vertices is zero, with logical variables left as external coefficients. All its connected cumulants vanish too. Absolute convergence permits taking this equality through the summed expansion and proves (6).

For finite beta this cancellation fails because the thermal auxiliary density matrix includes occupied states. No assertion that W(K_D)=0 at finite temperature is being made.

## 6. Time moments, scalar terms and sources

[proved] Two time moments follow from the exponential estimate with a loss in external field weight. If there are l external fields, the sum of their distances from an anchored time is at most l diam_t. Use H=2h. Then

    h^l (1+g sum_i |t_i-t_anchor|)^2
       <= C_(g/kappa) H^l exp(kappa diam_t),

since l^2/2^l is bounded and x^2 is bounded by 2 exp(x). A safe constant is 2+16(g/kappa)^2. Spatial moments are obtained in the same way by reserving part of a>0. Thus (5) supplies the finite time moments needed for low-frequency matching after this fixed norm loss.

[proved] Scalar terms can be bounded without pretending that Theorem 4.5's nonconstant kernel norm already contains them. Choose one spatial anchor for each translated interaction density, for example its leftmost site, and mark each residual density anchored at (j,t) by its own local external pair eta_bar_j(t) eta_j(t). Give all these fields zero covariance; a single global Grassmann pair would not implement the local rooted bound. Differentiating the resulting effective action at zero marker gives the rooted connected expectation of that density. Marking a residual vertex costs at most (H')^2 r in the input norm. Apply the same Cauchy estimate as in Section 5, extract the eta_bar eta coefficient by dividing the output norm by H^2, and integrate the interpolation parameter s. This bounds the scalar change per beta M by

    |r0| + 4 (H'/H)^2 r,                               (11)

provided the same smallness reserve holds. The |r0| term is necessary: a scalar perturbation has zero nonconstant kernel seminorm but shifts the pressure. The marker construction also controls local source insertions with finite coefficient norms. Additional derivatives use a further fixed Cauchy-radius reserve. Equation (11) controls the matching scalar; it does not yet identify a T^2 coefficient.

## 7. The embedding does not introduce low-energy competitors

[proved] The valid-state projection J J^dagger on each cell commutes with the enlarged Hamiltonian. Every interaction touching an invalid cell vanishes. Let

    w = sup_j sum_(X containing j) ||(B+R)_X||.

For an invalid set B_bad, cutting all interactions that meet B_bad costs at most w |B_bad|. Every invalid cell has auxiliary energy at least g. Reinstating the cut bonds in an all-valid trial state therefore proves an energy cost at least (g-w)|B_bad| above the all-valid ground. If w<g, the spectra below this positive threshold coincide with the original finite-cell model.

The corresponding partition-function argument factors the cut cells. Their valid one-cell heat trace is at least 2^n; the invalid trace is at most D_inv exp(-beta g), where D_inv is the finite number of invalid Fock states. Summing invalid subsets yields

    0 <= f_original-f_enlarged
       <= beta^(-1) log[1+(D_inv/2^n) exp(-beta(g-w))].

This extends the nearest-neighbor comparison in the imported auxiliary note to summable support interactions. It supplies a thermal comparison and does not claim uniform convergence of arbitrary nonlocal sources.

## 8. Application boundary for the current conjecture campaign

[proved conditional on the imported finite-order SW estimates] For a **fixed finite cell space**, C45's transformed Hamiltonian has k tending to zero with lambda and residual r bounded by C_emb C_SW lambda^3 in the strengthened norm. Choose lambda after the finite embedding constants, so (4) holds. Its logical scale is W=c_n lambda^2 with c_n>0 fixed. Equations (5)-(6) then give a relative zero-temperature matching error at most

    4 C_emb C_SW lambda/c_n -> 0.

No many-body gap in the code and no isolation of its extensive excited spectrum was used. This closes the finite-cell massive-matching implication under the stated norm and invariance hypotheses, rather than the whole infrared conclusion.

[unresolved for C44/C45] The actual mobile rooted representation includes arbitrarily long cells. A finite cutoff is not invariant under marker motion. Its spectral embedding has m increasing with the cutoff, and even the vacuum projector has coefficient norm (1+h^2)^m. The finite-cell theorem therefore cannot be asserted uniformly in that cutoff. A length-weighted or separately resummed infinite-cell theorem is still required, or a different exactly finite representation.

[unresolved after matching] The resulting logical action has general high-degree retarded vertices. Its weak-coupling multichannel RG and current/stress or thermal scaling limit must still be proved. The precise original hypotheses of the cited constructive-RG papers do not make those steps automatic. The c-d conjecture remains open on the evidence established here.
