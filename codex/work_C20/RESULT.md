# C20. Stability of the length mode and payload

## Status and scope

[Derived with standard bosonization input] In the usual product compactification and cell-translation description, there is no unavoidable non-Gaussian resonance. For each finite number of flavors, an open region of positive Gaussian stiffnesses and nonresonant momenta supports an infrared theory of `n + 1` compact bosons. The exclusion must include pairing and flavor conversion. Total parity alone does not forbid either process.

[Assumed] Two additional pieces of information are needed to apply this conclusion literally to the underlying lattice. They are the action of physical translation on the compact fields and the compactification and boundary-sector identifications imposed by the Ogata-Shiba gluing. The phrase “translation together with a shift tied to the length fluid” does not specify these data. Below, the cell-coordinate answer is explicit, and the dependence on these missing data is kept visible.

[Proved] Positivity and a finite nonresonance condition define an open set. Irrationality itself does not define an open set in the real density axis. Thus an “open set at irrational density” means either an open region intersected with the irrational-density subspace, or an open region containing irrational densities and also sufficiently high-denominator rational densities.

[Assumed] This analysis is closed book. The inputs are the setting in the question and the standard bosonization conventions stated below. There is no microscopic derivation of the gluing or a nonperturbative proof of the lattice phase in this report. Tags distinguish algebraic proofs, consequences of standard bosonization, assumptions, and perturbative physical reasoning. Numerical checks are identified as checks within those categories.

## 1. Complete local vertex lattice

[Assumed] First use a coordinate `x` measured in cells, a unit cell spacing, and the product compactification of one length liquid and `n` ordinary complex-fermion liquids. Write `vartheta_alpha` for the payload dual field, reserving `theta_alpha` for the Fermi-sea center specified in the question. The length fields remain `phi_l, theta_l`. The commutator convention is

\[
[\phi_i(x),\partial_y\vartheta_j(y)]=i\pi\delta_{ij}\delta(x-y),
\qquad \vartheta_0=\theta_\ell,
\]

[Assumed] A payload annihilation field has the long-distance form

\[
\psi_{\alpha s}(x)\sim
e^{i(\theta_\alpha+s k_\alpha)x}
e^{-i(s\phi_\alpha-\vartheta_\alpha)},\qquad s=\pm1.
\]

[Derived with standard bosonization input] A general charge-neutral-in-the-length-sector vertex is

\[
V_{m,a,b}=
\exp i\left[2m\phi_\ell+
\sum_{\alpha=1}^n(a_\alpha\phi_\alpha+b_\alpha\vartheta_\alpha)\right].
\tag{1}
\]

[Derived with standard bosonization input] Its complete integer constraints in this compactification are

\[
m\in\mathbb Z,\qquad a_\alpha,b_\alpha\in\mathbb Z,
\qquad a_\alpha\equiv b_\alpha\pmod2,
\qquad \sum_\alpha b_\alpha\equiv0\pmod2.
\tag{2}
\]

[Proved] To obtain (2), let `h_alpha,R` and `h_alpha,L` be the signed numbers of annihilation fields in a local monomial, counting creation fields negatively. Then

\[
a_\alpha=-h_{\alpha R}+h_{\alpha L},\qquad
b_\alpha=h_{\alpha R}+h_{\alpha L},
\qquad
h_{\alpha R}=\frac{b_\alpha-a_\alpha}{2},\quad
h_{\alpha L}=\frac{b_\alpha+a_\alpha}{2}.
\]

[Proved] Integer chiral counts give the first congruences. Total fermion parity gives the last congruence. Conversely these formulas construct integer chiral counts for every allowed `(a,b)`. In particular, `sum b` need not vanish. Imposing its vanishing would incorrectly add total payload U(1) conservation.

[Derived with standard bosonization input] Conservation of `r` forbids every vertex containing a nonzero coefficient of the length dual field. Length density harmonics are `exp[2im(phi_l - pi nu x)]`. This gives the remaining label `m`. Repeated identical fermions require point splitting or derivatives in their microscopic representative. The dimension of their bosonized primary already includes the resulting Pauli suppression. One must not add those derivatives again to that dimension.

[Derived with standard bosonization input] In the stated cell convention, the oscillatory factor multiplying (1) is `exp(i q x)`, with

\[
q(m,a,b)=-2\pi m\nu+
\sum_\alpha(b_\alpha\theta_\alpha-a_\alpha k_\alpha)
\quad (\bmod\ 2\pi).
\tag{3}
\]

[Derived with standard bosonization input] The conjugate vertex has momentum `-q`. A translation-invariant perturbation survives coarse graining without oscillation if and only if `q = 0 mod 2 pi`. Reciprocal-lattice momentum, rather than only the real number zero, is the correct umklapp test.

[Proved] The bosonic Hamiltonian locality check is also consistent with (2). The conformal spin of (1), up to reversal of the spatial orientation, is `sum a_alpha b_alpha / 2`. Since `a_alpha b_alpha = b_alpha mod 2`, total parity makes this an integer.

### Dimensions and the quadratic metric

[Derived with standard bosonization input] At a decoupled Gaussian point, the dimension is

\[
\Delta_0(m,a,b)=m^2K_\ell+
\frac14\sum_\alpha
\left(a_\alpha^2K_\alpha+\frac{b_\alpha^2}{K_\alpha}\right).
\tag{4}
\]

[Derived with standard bosonization input] Velocities do not enter (4). Formula (4) is not the dimension formula after arbitrary intercomponent quadratic couplings have been introduced.

[Assumed] To include all such couplings, let `N = n + 1`, order the fields as

\[
X=(\phi_\ell,\phi_1,\ldots,\phi_n;
\theta_\ell,\vartheta_1,\ldots,\vartheta_n)^T,
\]

[Assumed] Use the real symmetric stiffness matrix

\[
H_0=\frac1{2\pi}\int dx\, (\partial_xX)^T\mathcal H(\partial_xX),
\qquad
\mathcal H=\begin{pmatrix}A&W\\W^T&B\end{pmatrix},\qquad
\mathsf J=\begin{pmatrix}0&I_N\\I_N&0\end{pmatrix}.
\tag{5}
\]

[Assumed] Here `A` contains density stiffnesses, `B` current stiffnesses, and `W` density-current terms. The stated symmetries do not forbid `W`. If the question's quadratic Hamiltonian is restricted to density-density and current-current terms, set `W = 0`.

[Derived with standard bosonization input] For positive `H`, define the scaling metric

\[
\mathcal G=\mathcal H^{-1/2}
\left|\mathcal H^{1/2}\mathsf J\mathcal H^{1/2}\right|
\mathcal H^{-1/2},
\qquad |C|=(C^2)^{1/2}\text{ for real symmetric }C.
\tag{6}
\]

[Derived with standard bosonization input] For

\[
\lambda=(2m,a_1,\ldots,a_n;0,b_1,\ldots,b_n)^T,
\]

[Derived with standard bosonization input] The full dimension is

\[
\boxed{\Delta(\lambda)=\frac14\lambda^T\mathcal G\lambda.}
\tag{7}
\]

[Proved] One way to check (6) is to diagonalize `C = H^(1/2) J H^(1/2) = O diag(c_j) O^T`. The matrix `S = H^(-1/2) O diag(sqrt(|c_j|))` obeys

\[
S^T\mathcal H S=\operatorname{diag}(|c_j|),\qquad
S\operatorname{diag}(\operatorname{sgn}c_j)S^T=\mathsf J,
\qquad SS^T=\mathcal G.
\]

[Proved] Thus the transformation to normalized chiral fields gives (7) from their standard dimension formula. It also gives `G J G = J`.

[Derived with standard bosonization input] When `W = 0`, the result takes the more familiar form

\[
\mathcal G=\begin{pmatrix}\mathcal K&0\\0&\mathcal K^{-1}\end{pmatrix},
\qquad
\mathcal K A\mathcal K=B,
\qquad
\mathcal K=A^{-1/2}(A^{1/2}BA^{1/2})^{1/2}A^{-1/2}.
\tag{8}
\]

[Derived with standard bosonization input] For diagonal `A_ii = v_i/K_i` and `B_ii = v_i K_i`, this reduces to (4). Consequently the individual `K_alpha` alone are insufficient to specify dimensions in the coupled theory.

### All relevant and marginal operators

[Derived with standard bosonization input] The exhaustive criterion for a nonconstant vertex perturbation is (2), together with

\[
\lambda\ne0,\qquad q(\lambda)=0\pmod{2\pi},\qquad
\Delta(\lambda)\le2.
\tag{9}
\]

[Derived with standard bosonization input] For each such vertex one can add `g V + g* V†`. It is relevant for dimension below two and marginal at dimension two. The latter statement is a scaling classification, not a statement about its beta function.

[Derived with standard bosonization input] Descendants are also included by this prescription. A vertex times derivatives with total derivative level `d` has dimension `Delta + d` and the same momentum. Keep every independent descendant with `Delta + d <= 2`, identifying total derivatives and equation-of-motion redundancies. In particular, a current times a dimension-at-most-one vertex can be marginal or relevant. Excluding all primaries in (9) excludes all such non-Gaussian descendants as well.

[Derived with standard bosonization input] In the zero-vertex sector, the complete bulk list through dimension two consists of the identity, linear currents `partial_x X_i` of dimension one, and quadratic current products of dimension two. Linear currents select background densities or winding sectors. Quadratic products give (5). Single second derivatives are total derivatives on the ring. There are no further non-Gaussian operators in this sector through dimension two.

[Derived with standard bosonization input] The following table displays the processes most easily missed. Its dimensions use the decoupled metric (4). For coupled liquids, substitute their vertex labels into (7). Every row is subject to the momentum and dimension tests (9). Here `s,t = +/-1`, `alpha != beta`, `D_alpha = psi_alpha,R psi_alpha,L`, and `F_alpha = psi_alpha,R† psi_alpha,L`.

| Process | Oscillatory momentum `q` | Decoupled dimension |
|---|---|---|
| `[Derived with standard bosonization input]` Length harmonic `exp(2im phi_l)` | `-2 pi m nu` | `m^2 K_l` |
| `[Derived with standard bosonization input]` Intraflavor backscattering `F_alpha` | `-2 k_alpha` | `K_alpha` |
| `[Derived with standard bosonization input]` Intraflavor opposite-branch pairing `D_alpha` | `2 theta_alpha` | `1/K_alpha` |
| `[Derived with standard bosonization input]` Same-branch pairing, point split | `2(theta_alpha + s k_alpha)` | `K_alpha + 1/K_alpha` |
| `[Derived with standard bosonization input]` Flavour conversion `psi_alpha,s† psi_beta,t` | `p_beta,t - p_alpha,s` | `(K_alpha + 1/K_alpha + K_beta + 1/K_beta)/4` |
| `[Derived with standard bosonization input]` Interflavor pairing `psi_alpha,s psi_beta,t` | `p_alpha,s + p_beta,t` | `(K_alpha + 1/K_alpha + K_beta + 1/K_beta)/4` |
| `[Derived with standard bosonization input]` Opposite-branch pair transfer `D_alpha† D_beta` | `2(theta_beta - theta_alpha)` | `1/K_alpha + 1/K_beta` |
| `[Derived with standard bosonization input]` Density-wave exchange `F_alpha F_beta†` | `2(k_beta - k_alpha)` | `K_alpha + K_beta` |
| `[Derived with standard bosonization input]` Two-flavor umklapp `F_alpha F_beta` | `-2(k_alpha + k_beta)` | `K_alpha + K_beta` |
| `[Derived with standard bosonization input]` Length harmonic times any payload vertex | `-2 pi m nu + q_payload` | `m^2 K_l + Delta_payload` |
| `[Derived with standard bosonization input]` General higher umklapp or multipair process | Equation (3) with (2) | Equation (4) |

[Proved] For example, the same-branch pair has `a_alpha = -2s, b_alpha = 2`. Therefore its dimension is `K_alpha + 1/K_alpha >= 2`, with equality at `K_alpha = 1`. Treating it as an opposite-branch pair of dimension `1/K_alpha` would miss the Pauli constraint.

[Derived with standard bosonization input] All orientations of pair tunnelling, four-fermion processes involving three or four distinct flavors, and mixed length-payload umklapps are covered by (1)-(9), whether or not they appear by name in the table. For instance,

\[
e^{2i(m\phi_\ell+\sum_\alpha u_\alpha\phi_\alpha)}:
\quad q=-2\pi m\nu-2\sum_\alpha u_\alpha k_\alpha,
\quad \Delta_0=m^2K_\ell+\sum_\alpha u_\alpha^2K_\alpha.
\]

[Derived with standard bosonization input] An `h`-fold opposite-branch pair of one flavor has `q = 2h theta_alpha` and dimension `h^2/K_alpha`. Such higher processes cannot be omitted if the stiffnesses are unrestricted.

[Proved] Irrational `nu` removes every pure length umklapp, but does not remove mixed umklapp. A relation involving `m nu`, payload centers, and payload widths can still hold. Conversely, no nonzero vertex in (2) has identically vanishing cell momentum as a function of all independent `nu, theta_alpha, k_alpha`. An identity in these variables would force `m = a = b = 0`.

## 2. Exact nonresonance conditions and margins

[Proved] At a specified positive Gaussian point, the exact condition for absence of non-Gaussian relevant and marginal operators in this cell description is

\[
\boxed{
\operatorname{dist}\left(
-m\nu+\sum_\alpha\frac{b_\alpha\theta_\alpha-a_\alpha k_\alpha}{2\pi},
\mathbb Z\right)>0
\quad\text{for every allowed }\lambda\ne0
\text{ with }\lambda^T\mathcal G\lambda\le8.}
\tag{10}
\]

[Proved] This is both necessary and sufficient for absence of the operators, not necessary and sufficient for gaplessness. A surviving marginal operator may be marginally irrelevant, and some resonant perturbations need not gap a mode.

[Proved] The set tested in (10) is finite. If `G >= gamma I` with `gamma > 0`, then every candidate obeys

\[
4m^2+\sum_\alpha(a_\alpha^2+b_\alpha^2)\le\frac8\gamma.
\tag{11}
\]

[Proved] For `R = sqrt(8/gamma)`, a conservative upper bound on the number of candidates modulo conjugation is

\[
M_\gamma\le\frac12\left[
(2\lfloor R/2\rfloor+1)(2\lfloor R\rfloor+1)^{2n}-1\right].
\tag{12}
\]

[Proved] The bound deliberately ignores the parity restrictions and the spherical cutoff. A uniform lower bound on `G` gives a finite set that works throughout a stiffness neighborhood. Excluding the finite set of resonance surfaces for these candidates gives an open momentum region.

[Proved] A stronger all-orders condition is rational independence of

\[
1,\ \nu,\ \theta_1/(2\pi),\ldots,\theta_n/(2\pi),
\ k_1/(2\pi),\ldots,k_n/(2\pi).
\tag{13}
\]

[Proved] It is sufficient to eliminate every nonzero resonant vertex, regardless of dimension. It is also necessary for that stronger statement: multiply any integer relation by two to make every payload vertex coefficient even, which satisfies (2). Condition (13) is stronger than needed for stability and has no open interior. Arbitrarily high-order resonance surfaces are dense. There is no need to exclude them when their dimensions are uniformly above two.

### A finite list and an explicit bound near the free point

[Assumed] To quantify the size of the problem in terms of `n`, take a neighborhood of the free scaling metric satisfying

\[
\mathcal G\succeq(1-\epsilon)I,\qquad 0\le\epsilon<\frac13.
\tag{14}
\]

[Assumed] A strict operator-norm neighborhood `||G - I|| < epsilon` is one sufficient way to arrange (14). Different positive velocities at the decoupled free point are allowed.

[Proved] For every vector satisfying (2), `||lambda||^2` is divisible by four. An even `(a,b)` contributes zero modulo four. An odd `(a,b)` contributes two, and the number of such flavors is even by total parity. The length contribution is `4m^2`. Thus the only potentially relevant or marginal nonzero vectors under (14) have

\[
\|\lambda\|^2=4\quad\text{or}\quad8.
\tag{15}
\]

[Proved] Every other nonzero vector has dimension at least `3(1-epsilon) > 2`. Its uniform irrelevance margin is at least `1 - 3 epsilon`.

[Proved] These shells can be counted without assuming any flavor U(1). For one payload flavor, the possibilities through squared norm eight are a zero vector, four odd vectors `(+-1,+-1)` of norm two, four even vectors `(+-2,0),(0,+-2)` of norm four, and four even vectors `(+-2,+-2)` of norm eight. The number of odd-flavor vectors must be even. Including `m = +/-1` gives

\[
\begin{aligned}
N_4&=2+4n+16\binom n2,\\
N_8&=12n+48\binom n2+192\binom n3+256\binom n4.
\end{aligned}
\]

[Proved] For example, the term `192 binomial(n,3)` chooses one norm-four even vector and two norm-two odd vectors among three flavors. The term `256 binomial(n,4)` chooses four odd vectors. Length-payload combinations contribute twice the payload norm-four count to `N_8`.

[Proved] Conjugation has no fixed nonzero vector, so a sufficient candidate count is

\[
\boxed{M(n)=\frac{N_4+N_8}{2}
=1+8n+32\binom n2+96\binom n3+128\binom n4.}
\tag{16}
\]

[Proved] Some candidates yield the same resonance surface, so this is an upper bound on the number of distinct exclusions needed. It counts vertex pairs, not independent beta functions. For large `n`, its leading term is `(16/3)n^4`.

[Proved] Here is a quantitative existence argument. Set `tau_alpha = theta_alpha/(2 pi)` and `kappa_alpha = k_alpha/(2 pi)`, and use the parameter domain

\[
0<\nu<1,\qquad 0\le\tau_\alpha<1,\qquad0<\kappa_\alpha<\tfrac12.
\]

[Proved] Its Lebesgue volume is `2^(-n)`. For one nonzero candidate, the slab where the distance in (10) is at most `delta < 1/2` has volume `2 delta` times the volume of this domain. To see this, integrate first over a `tau_alpha` with nonzero `b_alpha`, or over `nu` when `m != 0`. Either coordinate traverses an integer number of periods. If both are absent, all `a_alpha` are even by (2), and integrating over a nonzero `kappa_alpha` again traverses an integer number of periods. This proves the slab-volume statement in every case.

[Proved] Subadditivity of volume bounds the excluded union by `2 M delta` times the domain volume. Consequently the strict margins

\[
\operatorname{dist}(q_\lambda/(2\pi),\mathbb Z)>\delta,
\qquad \delta=\frac1{4M(n)},
\tag{17}
\]

[Proved] These margins can be satisfied simultaneously on a nonempty open set. This is a sufficient margin of order `n^(-4)`, not an optimal bound or a necessary separation between every pair of Fermi points.

[Proved] For any prescribed irrational `nu`, the same argument applies to the payload variables with the pure length candidate omitted. A sufficient choice is

\[
\delta(\nu,n)=\min\left\{
\frac{\operatorname{dist}(\nu,\mathbb Z)}2,\frac1{4M(n)}\right\}.
\tag{18}
\]

[Proved] There is no uniform positive detuning bound depending only on `n` for every irrational `nu`, since an irrational density can lie arbitrarily close to zero or one. Away from the free neighborhood, other pure length harmonics can matter, and the corresponding minimum over `dist(m nu,Z)` must replace the first entry of (18).

[Proved] The margins are robust against small momentum shifts. In the shells (15), the sum of absolute coefficients in `q/(2 pi)` is

\[
|m|+\sum_\alpha(|a_\alpha|+|b_\alpha|)\le8.
\]

[Proved] Thus changing every normalized coordinate `(nu,tau,kappa)` by less than `delta/16` preserves a detuning greater than `delta/2`. This also quantifies how small interaction-induced shifts of the Fermi points must remain for the same exclusions to apply.

[Heuristic] A nonzero momentum cuts off the growth of an oscillating perturbation only after the length scale `1/|q|`, where `q` is reduced to its distance from a reciprocal vector. If `y` is its dimensionless bare coupling and `a_uv` a short-distance cutoff, perturbative control for a relevant primary requires

\[
|y|\,(a_{\rm uv}|q|)^{-(2-\Delta)}\ll1
\quad\text{when }a_{\rm uv}|q|\ll1.
\tag{19}
\]

[Heuristic] For a marginal perturbation with a quadratic beta function, a sufficient perturbative condition has the form `C |y| log[1/(a_uv |q|)] << 1`, with `C` determined by its operator products and the coupled flow. These are weak-coupling estimates, not universal exact thresholds. They explain why “small coupling” must be measured against the detuning margin, not only against a bandwidth. The statement of the problem does not determine a universal microscopic coupling bound as a function of `n`.

## 3. Quadratic stability

[Proved] The full quadratic theory has a positive energy for every nonzero field gradient if and only if

\[
\boxed{\mathcal H\succ0.}
\tag{20}
\]

[Proved] In the blocks of (5), the equivalent Schur-complement conditions are

\[
A\succ0,\qquad B-W^TA^{-1}W\succ0.
\tag{21}
\]

[Proved] For `W = 0`, both `A` and `B` must be positive definite. Positivity of each isolated `v_i,K_i` is not sufficient if the off-diagonal couplings are unrestricted.

[Proved] To prove (21), complete the square in the two gradient vectors `u,z`:

\[
u^TAu+2u^TWz+z^TBz
=(u+A^{-1}Wz)^TA(u+A^{-1}Wz)
+z^T(B-W^TA^{-1}W)z.
\]

[Proved] The necessity and sufficiency follow by independently varying the two vectors.

[Derived with standard bosonization input] The propagation eigenvalues are those of `J H`, up to a convention for which direction is called right moving. This matrix is similar to the real symmetric matrix `H^(1/2) J H^(1/2)`. Positive `H` makes these velocities real and nonzero, with `N` of each sign. Thus the positive quadratic theory has `N` right-moving and `N` left-moving gapless branches.

[Proved] A useful stiffness margin is the following. If `H = H_* + E`, `lambda_min(H_*) = h_* > 0`, then `||E||_op < h_*` ensures positivity. If every entry of the `2N` by `2N` perturbation is bounded in magnitude by `e`, the sufficient bound `2N e < h_*` follows from the row-sum norm. Entrywise smallness therefore needs an `O(1/n)` restriction in this conservative bound. This energy bound is separate from the scaling-metric condition (14).

[Derived with standard bosonization input] If an eigenvalue of `H` vanishes, a gradient direction has no quadratic energy and the nonzero-velocity Luttinger description fails. Higher spatial derivatives or nonlinear terms then become necessary. A negative eigenvalue makes the proposed homogeneous Gaussian state energetically unstable.

[Heuristic] Depending on the microscopic system, that instability can lead to phase separation, a changed current or density background, or another ordered or gapless state. Failure of positivity does not by itself determine which outcome occurs, and it is not equivalent to proving that one cosine gaps one mode. Conversely, positivity alone does not protect against a permitted relevant cosine.

## 4. Translation gluing and the scope of the answer

[Proved] Equations (2) and (3) use two logically distinct inputs: a lattice of local operators and a translation character on that lattice. To state the corresponding exact test for an explicitly glued theory, denote its local lattice by `Lambda_loc` and the translation phase of a vertex by `exp(i Q_T(lambda))`. If translation acts diagonally on these vertices, the test is

\[
\lambda\in\Lambda_{\rm loc},\quad \lambda\ne0,\quad
Q_T(\lambda)=0\pmod{2\pi},\quad
\tfrac14\lambda^T\mathcal G\lambda\le2.
\tag{22}
\]

[Proved] If translation permutes vertices or sectors, one must instead take its invariant linear combinations. A vertex phase cannot then be assigned without first resolving that action.

[Assumed] Equating (22) to the explicit cell test (9) requires that the gluing preserve the local lattice (2) and the zero-momentum criterion (3), possibly after an explicitly known change of variables. A gluing that only changes allowed ring sectors is compatible with all the preceding bulk conclusions. A gluing that attaches a length-dependent vertex or changes the translation shift must be included before using the cell momenta or the decoupled dimensions.

[Proved] For example, if the translation action in a fixed compactification is a field shift `X -> X + t`, its vertex character is `Q_T(lambda) = lambda^T t mod 2 pi`. Specifying cell Fermi momenta does not specify an additional gluing contribution to `t`. Likewise, a change of field coordinates changes both vertex coefficients and stiffnesses; changing just a momentum formula leaves an inconsistent dimension calculation.

[Proved] Therefore the stated symmetries and the informal reference to Ogata-Shiba gluing do not supply an exact physical momentum formula or an exact finite-ring selection rule. Equations (10), (16)-(18), and the absence of an unavoidable resonance are proved for the explicit cell implementation above. They become a conclusion about the microscopic system when its gluing is shown to implement that lattice and character. There is no evidence here of an unavoidable cell resonance, but an unspecified gluing is not a proof that none exists in the physical translation representation.

## 5. Infrared theory and ring spectrum

[Derived with standard bosonization input] Assuming the cell implementation, strict positivity, the uniform high-dimension margin, and the finite nonresonance margins, the Gaussian fixed theory has `N = n + 1` nonchiral bosons. In the usual convention its central charge is `c = n + 1`, or `c_R = c_L = n + 1` when the two chiral algebras are specified separately. Generic unequal velocities do not give a single common Lorentz velocity. They do not change the number of bosons or their contribution to the ground-state entanglement logarithm.

[Heuristic] For local couplings sufficiently small relative to the stiffness and momentum margins, the oscillatory terms renormalize the Gaussian parameters and the nonoscillatory higher-dimension terms flow to zero. The infrared theory is then the Gaussian theory just described, with emergent payload charge conservation. This is the perturbative Luttinger-liquid stability conclusion. The exact symmetries do not make the separate flavor charges exact microscopic quantum numbers.

### General positive Gaussian spectrum

[Assumed] Let `L_c` be the circumference in the cell coordinate. Let `w` denote the dimensionless winding vector, defined by `X(x + L_c) - X(x) = pi w` after subtraction of the chosen thermodynamic background. Let `W_glue` be its allowed affine lattice or union of affine sectors, including boundary twists and fixed parity. The microscopic gluing is needed to specify `W_glue` and the relation of `L_c` to the number of underlying lattice sites.

[Derived with standard bosonization input] For a given allowed winding sector, write the positive branch speeds as `v_j,R` and `v_j,L`, `j = 1,...,N`. The Gaussian energies are

\[
\begin{aligned}
E-e_\infty L_c-E_{\rm linear}
={}&\frac{\pi}{2L_c}w^T\mathcal H w
-\frac{\pi}{12L_c}\sum_{j=1}^{N}(v_{j,R}+v_{j,L})\\
&+\frac{2\pi}{L_c}\sum_{j=1}^{N}
\big(v_{j,R}\mathcal N_{j,R}+v_{j,L}\mathcal N_{j,L}\big)
+o(L_c^{-1}),
\qquad w\in\mathcal W_{\rm glue}.
\end{aligned}
\tag{23}
\]

[Derived with standard bosonization input] Here `mathcal N_j,s = sum_{p>=1} p n_j,s,p`, with nonnegative integer occupations, and `E_linear` accounts for the chemical potentials or other linear sources used to choose the background. The little-o statement applies in the stable infrared regime; it is absent for the strictly quadratic continuum theory itself. State-dependent boundary twists are evaluated in their corresponding sector.

[Derived with standard bosonization input] The first term follows directly by substituting the constant gradient `partial_x X = pi w/L_c` into (5). Each right- or left-moving bosonic branch contributes `-pi v/(12 L_c)` and oscillator spacing `2 pi v/L_c` by standard bosonic mode quantization. These give the remaining displayed terms in (23).

[Derived with standard bosonization input] In the unglued, decoupled, periodic-fermion example, write `w = (-delta N_i, delta J_i)`. Up to the chosen reference sea and twists, `N_i` and `J_i` are integer number and current labels with `N_i = J_i mod 2`, since `N_R = (N+J)/2` and `N_L = (N-J)/2` must be integers. Fixing total payload parity restricts `sum_alpha N_alpha mod 2`, rather than fixing each flavor number. The length gas may have its own boundary twist from its microscopic statistics. These example congruences are not a substitute for `W_glue`.

[Derived with standard bosonization input] In this decoupled example, (23) reduces to

\[
E-e_\infty L_c-E_{\rm linear}
=\frac{\pi}{2L_c}\sum_{i=0}^{n}v_i
\left(\frac{\delta N_i^2}{K_i}+K_i\delta J_i^2\right)
+\frac{2\pi}{L_c}\sum_{i=0}^{n}v_i
\left(\mathcal N_{i,R}+\mathcal N_{i,L}-\frac1{12}\right)
+o(L_c^{-1}).
\tag{24}
\]

[Derived with standard bosonization input] The corresponding zero-mode conformal weights can be written

\[
h_i=\frac18\left(\frac{\delta N_i}{\sqrt{K_i}}+
\sqrt{K_i}\delta J_i\right)^2,\qquad
\bar h_i=\frac18\left(\frac{\delta N_i}{\sqrt{K_i}}-
\sqrt{K_i}\delta J_i\right)^2,
\]

[Derived with standard bosonization input] Here right and left are interchanged if the orientation convention is reversed. Oscillators add their integer levels. The slowly varying continuum momentum contains `2 pi (h - bar h + N_R - N_L)/L_c`; Fermi-point phases and the gluing are also required for the full underlying-lattice translation eigenvalue.

### Irrational length density

[Proved] A finite ring with an integer number of cells and integer conserved charge cannot have `r/L_c = nu` for irrational `nu`. The thermodynamic density is implemented by a sequence of integer sectors `r_c` with `r_c/L_c -> nu`. In a background convention referenced to the infinite-system density, the length number zero mode is

\[
\delta N_\ell=r_c-\nu L_c.
\tag{25}
\]

[Proved] It is an offset integer lattice, not a continuous zero mode and not an additional gapless field.

[Derived with standard bosonization input] For a decoupled length liquid, its number and current contribution is

\[
E_{\ell,\rm zm}
=\frac{\pi v_\ell}{2L_c}
\left[\frac{(r_c-\nu L_c)^2}{K_\ell}
+K_\ell(J_\ell-\zeta_\ell)^2\right],
\tag{26}
\]

[Derived with standard bosonization input] This expression is subject to the applicable charge-current selection rule and boundary twist `zeta_l`. With coupling, the full quadratic form in (23) replaces (26), including cross terms and possible state-dependent twists.

[Proved] If the number sector can be minimized independently and a nearest integer is allowed, then `|r_c - nu L_c| <= 1/2`. Its contribution is of order `1/L_c`. The offset generally does not approach a fixed value as the ring grows, so scaled finite-size energies can show shell-dependent shifts. Restrictions from gluing or parity can change which nearest allowed sector is used. Conservation of `r` also permits working in one chosen sequence of fixed sectors instead of minimizing over them.

[Proved] In the independent-number example the extra ground-state term is

\[
\frac{\pi v_\ell}{2K_\ell L_c}
\operatorname{dist}(\nu L_c,\mathbb Z)^2.
\]

[Proved] It must be accounted for before identifying the Casimir coefficient from ground-state energies. Taking a rational finite-size approximation does not create a fixed low-order thermodynamic umklapp. If `r_c/L_c` is reduced to lowest terms, its denominator must become unbounded along any sequence converging to an irrational limit. The order and, at fixed positive scaling metric, the dimension of its pure length umklapp therefore become unbounded.

## 6. Checks and conclusion

[Proved] The `n = 0` limit leaves only the length vertices. Irrational `nu` excludes all of their nonzero harmonics, and the positive Gaussian theory has one boson. Formula (16) gives one candidate conjugate pair near the free point, namely `m = +/-1`.

[Derived with standard bosonization input] As a resonance check, setting every `theta_alpha = 0` makes intraflavor opposite-branch pairing translation invariant. At the free point it has dimension one and can gap payload modes. Pair transfer between flavors is also then nonoscillatory and has free dimension two. This special case demonstrates why generic centers cannot be silently replaced by inversion-symmetric Fermi seas.

[Derived with standard bosonization input] As a mixed-umklapp check, when `k_alpha = pi nu`, the vertex with `m = 1, a_alpha = -2, b = 0` has zero cell momentum and dimension `K_l + K_alpha` in the decoupled theory. It is marginal when both parameters equal one. Irrational length density alone therefore does not establish stability.

[Proved] The independent integer enumeration in `check.py` reproduces (16) for `n = 1,...,8`, giving respectively `9, 49, 217, 737, 1961, 4369, 8569, 15297` conjugate vertex pairs. This checks the finite combinatorics rather than an interacting phase diagram.

[Derived with standard bosonization input] Numerical algebra checks of (6) give a maximum residual of `4.98e-15` in `G J G = J` for a fixed positive coupled matrix, and `6.66e-16` in its decoupled reduction to `diag(K,K^(-1))`. They support the formulas but do not replace their derivation. The initial command named `python` was unavailable in this environment. The complete check passed using `python3 check.py` and NumPy.

[Derived with standard bosonization input] The conditional answer is yes. With the explicit local lattice and cell translation character above, positive stiffness, finite nonresonance margins, and sufficiently weak perturbations give an infrared `n + 1`-boson Gaussian theory on an open parameter region. No individual flavor U(1) is needed, and no unavoidable resonance appears in this implementation.

[Assumed] The remaining input for an unconditional statement about the system in the question is its explicit gluing map. It must determine (i) the physical translation character and any dressing of local vertices, (ii) the allowed winding and parity sectors, and (iii) the conversion between cell circumference and underlying-lattice circumference. Without these data, the full physical momentum classification and the exact ring-sector spectrum are undetermined. The cell result and the Gaussian energy formula are complete under their stated assumptions.
