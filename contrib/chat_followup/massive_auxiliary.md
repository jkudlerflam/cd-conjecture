# Direct-sum embedding into free massive auxiliary fermions

## Claim and status

An arbitrary finite, parity-graded cell with an exactly degenerate fermionic code can be embedded into a larger fermion Fock space whose **on-site Hamiltonian is quadratic**. The embedding need not factor the original cell into logical and auxiliary tensor factors. The additional states are excluded from the low-energy spectrum by a volume-independent energy bound at sufficiently weak intercell coupling; their contribution to the free-energy density is exponentially small at low temperature. These statements are [proved] below.

This is new relative to [the existing infrared-control analysis](../ir_control/README.md), whose general-cell reduction encountered an obstruction to tensor factorization and whose exactly integrated example assumed a special auxiliary Hamiltonian. It gives a different representation of a general cell. It does **not**, by itself, prove a uniformly small effective action after auxiliary integration or the claimed multicomponent Luttinger-liquid infrared limit. Those further statements remain conditional here.

**Provenance and checking.** Imported on 2026-10-27 from the user's earlier research discussion, with AI assistance from Codex; the source discussion's model was not recorded. The finite-dimensional and partition-function arguments below were rederived during import, with explicit constants and hypotheses added. They have not received an independent review. The source's blanket claim that the massive integration was already controlled uniformly in volume is downgraded to a proposed application of a specific convergence theorem. No numerical calculation is asserted or required to reproduce the algebra below.

## 1. Precise finite-cell setting

Let \(\mathcal K\) be finite-dimensional, with grading \(\Pi\), and let

\[
H_0\ge0,\qquad [H_0,\Pi]=0,\qquad
\ker H_0=G\mathcal F_n,\qquad
\Pi G=G(-1)^{N_{\rm log}},
\]

where \(G\) is an isometry and \(\mathcal F_n\) is the Fock space of \(n\) logical complex fermions. Assume

\[
H_0\ge\Delta(1-GG^\dagger),\qquad \Delta>0.
\]

Choose an orthonormal excited eigenbasis \(|s\rangle\), with energies \(E_s\ge\Delta\) and parity \(\pi_s=\pm1\). There are finitely many such states, including multiplicities.

An optional torus symmetry \(U(1)^k\) can also be preserved. Assume it commutes with \(H_0\) and \(\Pi\), and choose the excited basis with charges \(C_s\in\mathbb Z^k\). On the code require

\[
C(\mathbf n)=C_0+\sum_{a=1}^n n_a w_a.
\]

For the species-torus version, assume there is a nonzero parity homomorphism \(\chi:\mathbb Z^k\to\mathbb Z_2\) such that \(\chi(w_a)=1\) and the grading on a charge-\(C\) state is \((-1)^{\chi(C-C_0)}\). This uses parity relative to the code vacuum. If the initially chosen physical parity assigns that vacuum a minus sign, its overall on-site sign must first be adjusted; the induced even/odd operator grading is unchanged. This baseline convention is necessary for a grading-preserving map to the usual even Fock vacuum.

## 2. Exact embedding and charge assignment [proved]

For each excited basis vector introduce a dedicated set of auxiliary modes:

\[
k_s=\begin{cases}1,&\pi_s=-1,\\2,&\pi_s=+1.\end{cases}
\]

Different states use disjoint modes. Give each of the \(k_s\) modes energy \(E_s/k_s\), put \(m=\sum_s k_s\), and write

\[
\widehat{\mathcal K}=\mathcal F_n\widehat\otimes\mathcal F_m,
\qquad
H_{\rm mass}=\sum_s\sum_{\nu=1}^{k_s}
\frac{E_s}{k_s}a_{s,\nu}^\dagger a_{s,\nu}.
\]

Here the tensor product is graded. With a fixed ordering of auxiliary creation operators, define

\[
J G|\mathbf n\rangle=|\mathbf n\rangle_{\rm log}\otimes|\Omega\rangle,
\qquad
J|s\rangle=|0\rangle_{\rm log}\otimes
\left(\prod_{\nu=1}^{k_s}a_{s,\nu}^\dagger\right)|\Omega\rangle.
\]

The image vectors are orthonormal, and their parities agree with those of their preimages. Consequently \(J\) is an even isometry. Each image vector has precisely the original on-site energy, so

\[
H_{\rm mass}J=JH_0.
\]

For an independent torus symmetry, assign the auxiliary charge vectors any values summing to \(C_s-C_0\), and give the enlarged Fock vacuum baseline charge \(C_0\). If the torus itself implements parity as in Section 1, every auxiliary mode can be assigned odd \(\chi\)-charge: for odd \(C_s-C_0\), use that vector for the single mode; for even \(C_s-C_0\), choose any integer vector \(u\) with \(\chi(u)=1\) and use \(u\) and \(C_s-C_0-u\). Both are odd. Logical modes carry charges \(w_a\). These assignments make \(J\) intertwine the torus as well as the grading. Auxiliary charges may be negative or large; these modes are an analytical representation, not new microscopic sites.

Let \(P=JJ^\dagger\) denote the valid-state projection and \(Q=1-P\). The valid subspace is

\[
\bigl(\mathcal F_n\otimes|\Omega\rangle\bigr)
\oplus
\operatorname{span}\left\{
|0\rangle\otimes\prod_\nu a_{s,\nu}^\dagger|\Omega\rangle
\right\}_s.
\]

It is invariant under \(H_{\rm mass}\). Every occupation state in its orthogonal complement contains an auxiliary particle, because all logical states with auxiliary vacuum are already valid. Every auxiliary energy is at least

\[
g=\Delta/2,
\qquad QH_{\rm mass}Q\ge gQ.
\]

For an operator \(O\) on \(\mathcal K\), its extension \(JOJ^\dagger\) has the same norm, vanishes on \(Q\widehat{\mathcal K}\), and preserves the valid subspace. It is a finite polynomial in the enlarged creation and annihilation operators: matrix units between occupation states are products of creation operators, annihilation operators, and vacuum projectors \(\prod_\alpha(1-n_\alpha)\). Its polynomial degree and coefficient norm can be very large, but are finite for a fixed cell. Even operators have even polynomials, and symmetry invariance is preserved.

This is an embedding, not a unitary factorization of \(\mathcal K\). In particular, it does not contradict the parity-multiplicity obstruction in [the previous leakage analysis, Section D](../ir_control/astra_B/RESULT.md).

## 3. Invalid sectors have a uniform energy cost [proved]

Consider \(M\) cells on a ring, already expressed as a graded fermionic chain:

\[
H_M=\sum_j H_0(j)+\lambda\sum_j V_{j,j+1},
\qquad \|V_{j,j+1}\|\le w,
\]

where each bond is even. For an ordinary tensor-product spin chain, first perform the appropriate Jordan–Wigner identification; the usual parity-dependent boundary condition on a ring must be retained. Even nearest-neighbor bonds stay local in this representation.

Extend a bond by the two-cell even isometry:

\[
\widetilde V_{j,j+1}=(J_j\widehat\otimes J_{j+1})
V_{j,j+1}(J_j\widehat\otimes J_{j+1})^\dagger,
\quad
\widetilde H_M=\sum_jH_{\rm mass}(j)
+\lambda\sum_j\widetilde V_{j,j+1}.
\]

All the \(P_j\) commute with \(\widetilde H_M\). The sector with every site valid is exactly isospectral to \(H_M\). In a sector with an invalid-site set \(B\), every bond incident on \(B\) is zero, so at most \(2|B|\) bonds have been removed.

Let \(E_{\overline B}\) be the ground energy of the valid cells outside \(B\), with those incident bonds cut. The cut all-valid problem has the same energy, since the isolated valid cells in \(B\) can each occupy the zero-energy code. Reinstating the cut bonds gives

\[
E_{\rm valid}\le E_{\overline B}+2|B|\,|\lambda|w.
\]

The invalid-sector energy obeys

\[
E_B\ge E_{\overline B}+|B|g
\ge E_{\rm valid}+|B|g_{\rm eff},
\qquad g_{\rm eff}=g-2|\lambda|w.
\]

Thus, whenever \(2|\lambda|w<g\), no invalid sector contributes below \(E_{\rm valid}+g_{\rm eff}\), uniformly in \(M\). The enlarged and original spectra agree exactly in this energy window, including multiplicities. This statement requires no gap within the logical chain and no perturbative block diagonalization.

## 4. Free-energy comparison [proved]

Write \(D=\dim\mathcal K\), \(D_{\rm inv}=2^{n+m}-D\), and

\[
z_{\rm v}(\beta)=\operatorname{Tr}_{\mathcal K}e^{-\beta H_0}
\ge2^n,\qquad
z_{\rm i}(\beta)=\operatorname{Tr}_{Q\widehat{\mathcal K}}
e^{-\beta H_{\rm mass}}\le D_{\rm inv}e^{-\beta g}.
\]

For a fixed \(B\), factoring the isolated on-site traces gives

\[
Z_B=Z_{\overline B}z_{\rm i}^{|B|},
\qquad Z_{\rm valid}^{\rm cut}=Z_{\overline B}z_{\rm v}^{|B|}.
\]

The elementary finite-dimensional bound
\( |\log\operatorname{Tr}e^{-\beta(A+K)}-
\log\operatorname{Tr}e^{-\beta A}|\le\beta\|K\| \)
follows by differentiating along \(A+tK\). Applying it to the removed bonds proves

\[
\frac{Z_B}{Z_{\rm valid}}
\le\left[A_{\rm inv}e^{-\beta g_{\rm eff}}\right]^{|B|},
\qquad A_{\rm inv}=D_{\rm inv}/2^n.
\]

Summing over invalid-site subsets yields, for free energy per cell,

\[
0\le f_{\rm physical}(T)-f_{\rm enlarged}(T)
\le T\log\left(1+A_{\rm inv}e^{-g_{\rm eff}/T}\right).
\]

For no invalid states, take \(A_{\rm inv}=0\). The constants depend on the fixed cell, not the ring length. Per original physical site, divide the bound by the fixed cell length. Hence the two models have the same coefficient of a low-temperature power-law free-energy expansion, if such an expansion exists. This does not establish the existence or the CFT interpretation of that expansion.

The enlargement leaves the original microscopic dimension \(d\) unchanged: it is a comparison representation of the blocked problem, with an explicitly controlled extra sector.

## 5. What the Gaussian representation does and does not establish

Let \(h\) be the compression of a physical bond onto two code spaces and \(\widehat h\) its extension acting on the logical tensor factors in every auxiliary sector. Then

\[
\widetilde H_M=H_{\rm mass}+\lambda\sum_j\widehat h_{j,j+1}
+\lambda\sum_jK_{j,j+1},\qquad
K=\widetilde V-\widehat h,
\]

with \(\|K\|\le2w\) and

\[
\langle\Omega,\Omega|K|\Omega,\Omega\rangle_{\rm aux}=0.
\]

These are [proved] identities. If the compression is quadratic, the first two terms are a Gaussian reference Hamiltonian on independent logical and auxiliary modes. An approximate quadratic compression has an additional purely logical error which must be estimated separately. In contrast to the original-cell factorization attempt, \(H_{\rm mass}\) now acts only on auxiliary modes.

For an auxiliary energy \(e\ge g\), its antiperiodic thermal covariance is spatially on-site and, for \(0<\tau<\beta\), has magnitude

\[
|C_e(\tau)|=\frac{e^{-e\tau}}{1+e^{-\beta e}},
\qquad
\int_0^\beta|C_e(\tau)|\,d\tau
=\frac{1-e^{-\beta e}}{e(1+e^{-\beta e})}\le\frac1g.
\]

Thus its summable time-decay bound is uniform in temperature and volume. The vacuum expectation of the residual interaction is zero. With normal ordering and the corresponding point-split convention, its first auxiliary cumulant vanishes at zero temperature. At positive temperature the first cumulant can be nonzero; for a fixed finite polynomial its coefficients are bounded by a cell-dependent constant times \(|\lambda|e^{-\beta g}\). Thermal auxiliary occupations cannot simply be dropped.

**Relevant literature [read], with an incomplete application.** Pedra and Salmhofer, *Determinant Bounds and the Matsubara UV Problem of Many-Fermion Systems*, Theorem 4.5, gives analyticity and a remainder bound for general interaction vertices in a pinned-coordinate kernel norm. It requires finite determinant and decay bounds and smallness in a stronger field-weighted interaction norm. The theorem retains external fields and is not restricted to quartic vertices. [Primary paper, Section 4.3](https://www.wias-berlin.de/people/koenig/www/FG/PedraSalm.pdf#page=18).

An intended application would use covariance zero on the external logical fields and the massive covariance on auxiliary fields. To claim the desired theorem for the present Hamiltonians, still supply the following details:

1. Specify a coefficient norm for the normally ordered, possibly very high-degree bond polynomials. Operator norm \(\|K\|\le2w\) alone is not the needed kernel-norm bound; the conversion constants depend on the fixed cell and number of modes.
2. Verify uniform determinant bounds, the field-weighted smallness condition, and the treatment of equal-time vertices in a time discretization with a controlled continuum limit.
3. Prove the weighted spatial and temporal kernel estimates needed by the later infrared argument, with the logical fields left external. The effective action's constant term must also be tracked if the observable is free energy.
4. Keep the positive-temperature first cumulant and thermal normalization, and verify the order of the thermodynamic, time-continuum, and zero-temperature limits.

The proposed result is a bound of the form

\[
\|\mathcal V_{\rm induced}\|_{\rm loc}
\le C_{\rm cell}\frac{\lambda^2}{\Delta}
+C'_{\rm cell}|\lambda|e^{-\beta g}
\]

for small enough \(|\lambda|/\Delta\), in a stated local kernel norm, with constants independent of \(M\) and \(\beta\). **This bound is conditional in this note.** Reading the convergence theorem makes the route plausible; it is not a substitute for checking its hypotheses and the limiting argument. The direct-sum embedding and the free-energy estimate above remain valid independently of this unfinished step.

This qualification also preserves the finite-temperature counterexample in [the previous leakage analysis, Section D](../ir_control/astra_B/RESULT.md): for general original cells, a coupling-only error bound that vanishes at \(\lambda=0\) cannot identify their finite-temperature logical state with a Gaussian state. The enlarged comparison introduces exponentially small thermal effects and changes the definition of logical observables away from the code; it does not erase that counterexample.

## 6. Scaling and the remaining infrared obligation

The source discussion's first-order proposal has logical scale

\[
\eta=|\lambda|p_Lp_R\epsilon^2,
\]

where \(p_L,p_R>0\) are endpoint clock weights and \(\epsilon>0\) is the compression scale. If the purely logical encoding-error kernel is bounded by \(C|\lambda|p_Lp_R\delta_{\rm enc}\), and if the conditional massive-integration estimate above holds, rescaling time by \(\eta\) gives at zero temperature

\[
\rho\le C_1\frac{\delta_{\rm enc}}{\epsilon^2}
+C_2\frac{|\lambda|}{\Delta p_Lp_R\epsilon^2}.
\]

This is a **conditional scaling estimate**, not an established uniform theorem for the microscopic construction. At nonzero temperature there is also the rescaled thermal term. The constants and nonresonance/velocity thresholds must be controlled at the chosen fixed cell. Reducing \(\lambda\) suppresses the second term but cannot repair a too-large encoding error in the first; any dependence of the first constant on the requested encoding accuracy must be resolved before declaring that accuracy arbitrarily attainable.

Even after the integration estimate, one must establish the low-energy behavior of the resulting one-dimensional fermions with high-degree, retarded interactions. The proposed inputs are simple Fermi points, an exact species torus, the relevant momentum nonresonances, small interactions relative to all logical velocities, and a marginal flow that remains in the positive-stiffness multiboson regime. A symmetry-based list of density/current marginal operators does not alone prove this flow or its conformal thermodynamics. In particular, a complete comparison/RG argument and the inference of total **Virasoro** central charge \(c=n\) remain missing.

The retrieved source's proposed RG statement ends before a complete theorem is given; no unseen proof is imported here. The earlier repository's separate questions about microscopic compilation, frozen versus mobile cells, and finite-coupling infrared control also remain in force. This contribution improves the representation and separation of excited-cell states; it is not a completed counterexample to the c–d conjecture.
