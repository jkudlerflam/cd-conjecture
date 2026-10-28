# Localization on the fermionic Matsubara grid

Status: **[proved finite-temperature localization lemma]**. This supplies
an explicit one-step localization map whose difference from line Taylor
localization has a thermal power greater than two. It does not construct
the complete RG flow.

Provenance: OpenAI GPT-6-family Codex thermal worker, 2026-10-28.
AI-assisted; no human review.

## 1. A four-point projection for a quadratic kernel

Let a=pi T. For a C^4 function f on [-3a,3a], define

    A_T f = [9(f(a)+f(-a))-(f(3a)+f(-3a))]/16,
    B_T f = [27(f(a)-f(-a))-(f(3a)-f(-3a))]/(48a),
    (L_T f)(omega)=A_T f+omega B_T f.                  (1)

Every evaluation in (1) is at an allowed fermionic Matsubara frequency.
The maps extract the constant and linear coefficients of the cubic
interpolating polynomial at those four nodes. In particular,

    A_T(1)=1,  A_T(omega)=A_T(omega^2)=A_T(omega^3)=0,
    B_T(omega)=1,
    B_T(1)=B_T(omega^2)=B_T(omega^3)=0.                 (2)

Thus L_T is a projection onto affine polynomials: L_T^2=L_T. It
reproduces the same relevant and marginal polynomial space as
L_infinity f=f(0)+omega f'(0). No new quadratic coupling is introduced.
The symmetric real weights preserve any conjugation relation of the form
f(-omega)=conjugate(f(omega)), with the corresponding parity of the
linear coefficient. The construction acts separately on the allowed
branch and charge labels, so it does not alter their selection rules.

**[proved error estimate]** With M4=sup|f''''| on that interval,

    |A_T f-f(0)| <= (15/32) M4 a^4,
    |B_T f-f'(0)| <= (3/16) M4 a^3.                   (3)

Proof. Subtract the Taylor polynomial of degree three at zero, which
(2) reproduces in the two coefficients. Its remainder at x has absolute
value at most M4|x|^4/24. The sum of the absolute A weights times |x|^4
is (18+162)a^4/16. For B it is (54+162)a^3/48. Divide by 24 to obtain
(3). This proof works for complex f by the integral Taylor remainder.

## 2. Application to a scaled, retarded quadratic endpoint

Let K_2 be an integrable line kernel, or finite measure, with

    s^(-1) integral |K_2(t,x)|(1+s(|t|+|x|))^4 <= u,
    0<T<=s.                                          (4)

For its Fourier transform f(omega,k), differentiation under the measure
gives sup|partial_omega^4 f|<=u s^(-3). Consequently

    |A_T f-f(0)|/s <= (15 pi^4/32)u(T/s)^4,
    |B_T f-partial_omega f(0)|
                            <=(3 pi^3/16)u(T/s)^3.    (5)

These are exactly the normalized relevant and temporal-kinetic
coefficients. In particular, their comparison has a common thermal
power p=3. If the quadratic kernel contains an already-localized affine
piece, it cancels identically by (2), so (4) need only control its
integrable remainder. There is no attempt to treat a derivative delta
as a finite measure.

The spatial kinetic coefficient can be extracted by A_T applied to
partial_k f at the Fermi point. Requiring a fifth spacetime moment in
(4) bounds partial_omega^4 partial_k f by u s^(-4), and gives

    |A_T partial_k f-partial_k f(0)|
                              <= C u(T/s)^4.          (6)

Here momentum is continuous after the spatial thermodynamic limit.
Applying (1) in frequency and first-order Taylor extraction in momentum
therefore handles the complete local quadratic polynomial with at most
five moments. Finite-volume momentum interpolation is a separate issue
and is not needed for the thermodynamic-first pressure calculation.

The decomposition using L_T is exact: f=L_T f+(I-L_T)f. Moreover

    (I-L_T)f=(I-L_infinity)f+(L_infinity-L_T)f.         (7)

The first term has the ordinary quadratic Taylor remainder. The second
is the explicitly bounded affine thermal correction (5)-(6). It is
retained in the RG estimate, not declared zero. This is how one avoids
assuming a uniform zero at unavailable frequency zero on the thermal
Matsubara grid.

## 3. Quartic localization at legal frequencies

For a number-conserving quartic kernel choose three independent fermion
frequencies omega_1,omega_2,omega_3 and set

    omega_4=omega_1-omega_2+omega_3.

When each independent frequency is an odd multiple of pi T, so is the
fourth. Evaluate the tensor product A_T^(1)A_T^(2)A_T^(3) to extract the
local quartic constant. This is a projection onto constant quartic
kernels and uses only legal conserving Matsubara assignments. Other
sign conventions give an equivalent odd sum with the same property.

Suppose each pure fourth partial derivative in the three independent
frequency variables is bounded by M4 on the node cube. The absolute
weight norm of A_T is 5/4. Telescope the three tensor factors against
evaluation at zero and apply (3) at each step. This yields

    |(A_T tensor A_T tensor A_T)f-f(0,0,0)|
       <= (15/32)[1+5/4+(5/4)^2] M4 a^4.              (8)

For a quartic line kernel of normalized fourth-moment norm u,
M4<=C u s^(-4); the constant accounts for the chosen independent
frequency convention. Thus (8) is O(u(T/s)^4). The projected quartic
term has the same charge/branch content as the original localization.
Nonresonant branch combinations can be excluded before applying this
map, exactly as in the zero-temperature construction.

The displayed tensor formula distinguished one leg when eliminating the
conservation constraint. Before using it as a Grassmann-kernel localization,
restore the fermionic permutation symmetries. Let Alt be normalized
antisymmetrization over the four complete field labels, including their
coordinates and creation/annihilation labels. Use

    L_(T,4)^ferm = Alt L_(T,4) Alt.                    (8a)

Equivalently, average the formula over distinguished-leg/permutation
choices with their fermionic signs. Its image consists of antisymmetric
constant quartic tensors. It is the identity on that space, hence remains
a projection. On antisymmetric input, subtracting the line localization
commutes with the outer Alt, so the estimate follows by averaging the
permuted versions of (8). The diameter coefficient norm is invariant under
these permutations; any different anchored norm convention costs at most
a fixed finite permutation factor. The thermal power is unchanged.

Charge labels move with the fields, so this operation preserves number
conservation. If the chosen coefficient convention has a separate
Hermitian-conjugation involution, average also with its conjugate map.
Both maps have the same constant image, act as its identity, and have
symmetric real Matsubara weights. This retains the projection and the
conjugation relation. The derivative bound is required for each permuted
frequency assignment; the global Fourier moment bound above supplies it.

## 4. What periodization and this lemma establish together

For an exactly periodized line kernel, its Fourier coefficients at allowed
Matsubara assignments equal the line Fourier transform at those same
assignments. Therefore (5)-(8) are a bound on the commutator between
periodization and localization. The maps are linear. Applied to a
reference-subtracted kernel of normalized size omega(s), their error is

    O(omega(s)(T/s)^3),                               (9)

with the appropriate quadratic/marginal dimensional factors restored.
The single-scale determinant estimate in WINDING_ESTIMATE.md supplies
the other forcing term, from integrating a shell. Together they give
the required one-step finite-temperature comparison under a common
five-moment norm.

Kernels must be represented before restricting their external fields to
the next shell. The effective-action map with unrestricted external
Grassmann sources supplies that representation, so the four sampled
frequencies remain defined even near the last thermal scale. Their
values are not guessed from a polynomial known only on a smaller
external Fourier support.

The actual finite-temperature running kernels also differ because all
higher RG steps were performed at temperature T. A global history
estimate must propagate both forcing terms, retain the normalization
ratios and scalar prefactors, and prove a common moment/field-weight
reserve. This lemma does not imply that bootstrap by itself. Its finite
temperature coefficients are outputs of the flow with one fixed bare
Hamiltonian; no physical counterterm is retuned as a function of T.
