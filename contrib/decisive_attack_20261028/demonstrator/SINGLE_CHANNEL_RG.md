# Single-channel retarded RG: exact reductions and the remaining theorem

Status: **proved finite-cell symmetry, norm and ultraviolet reductions;
proved conditional marginal-flow lemma; audited published inputs.
Not a completed RG/CFT theorem.**

This note concerns only the one logical fermion per PBP cell demonstrator.
It supplements [the exact tuning calculation](README.md). The physical local
dimension remains ten. The enlarged fermion spaces below are proof devices.
No statement here settles selection of the PBP sector among all type sectors.

The charge-conjugation implementer and equal-time convention in Sections
3--4 have a separate read-only AI audit by the `frozen_scope` agent. Seventeen
32-by-32 checks of the physical cell and of one encoded excited pair give
zero numerical residual, including all possible even charge offsets
`q_s=-2,0,2`. The proofs below do not rely on those numerical checks.
The reproducible script [verify_single_channel.py](verify_single_channel.py)
also checks twenty-four scalar heat-kernel identities, twenty-seven matrix
cutoff identities, finite-ring filter norms and the flow-sum constant.
Run it with a Python environment containing NumPy; the saved output is
[single_channel_checks.json](single_channel_checks.json). No infrared-flow
or central-charge assertion is tested by this script.
The scalar ultraviolet argument has two separate AI reviews, including
the explicit [UV_AUDIT.md](../exact/UV_AUDIT.md). Its general row-filter
lemma also explains the determinant step of the matrix extension. None
of these reviews is a human review.

## 1. The regime that must be proved

Take fixed `D>0`, `0<|kap|<D` and nonzero `Jp,Jz`. The second-order code
Hamiltonian, with an irrelevant scalar suppressed, is

    lambda^2 sum_j [Jxy (s+_j s-_(j+1) + h.c.) + Jzz Z_j Z_(j+1)],

    Jxy = Jp Jz kap^2 / [4D(D^2-kap^2)],
    Jzz = 3 Jp^2 kap^4 / [32D(D^2-kap^2)(4D^2-kap^2)].

Thus `2 Jzz/|Jxy|` can be made arbitrarily small by choosing the fixed ratio
`|Jp/Jz|` first. The hopping bandwidth scale is `W=|Jxy| lambda^2`.
After finite-cell massive matching, division by W makes the additional
zero-temperature action `O(lambda/|Jxy|)` in an exponential kernel norm.
Choose lambda only after the hopping ratio and all finite-cell constants.
Both choices precede the thermodynamic and zero-temperature limits.

The extra action has vertices of unbounded external degree and is generally
nonlocal in imaginary time. A theorem for a specified instantaneous
density-density Hamiltonian cannot simply be applied to that action.

## 2. A closer published result, and its exact scope

**[read]** Benfatto, Gallavotti and Lebowitz,
[*Disorder in the 1D spinless Holstein model*](https://arxiv.org/pdf/chao-dyn/9505010),
Sections 2--3, integrate quantum optical phonons. Equation (2.7) is a
time-retarded quartic density interaction, with the exponential decay in
(2.4). After ultraviolet integration, (2.10)--(2.11) give all-degree kernels
with exponential connecting-path bounds and degree-dependent coupling powers.
The paper applies the BGPS infrared construction to these kernels and proves
the anomalous fermion two-point asymptotics (1.6). Its microscopic interaction
is the specified Holstein interaction, not arbitrary even kernels. Section 1
explicitly does not prove ground-state uniqueness. It supplies neither a
thermal central-charge theorem nor an abstract initial-kernel theorem.

**[read]** Benfatto, Gallavotti, Procacci and Scoppola,
[*Beta function and Schwinger functions for a many fermions system in one
dimension. Anomaly of the Fermi surface*](https://doi.org/10.1007/BF02099791)
([author-hosted source](https://ipparco.roma1.infn.it/pagine/deposito/1980-1993/bb.ps.gz)),
Theorem 1 and (3.18), bound their ultraviolet kernels. Theorem 2 and (5.52)
give convergent renormalized trees conditional on small running couplings.
Section 7 compares the marginal beta function with the Luttinger reference;
(7.28) controls dependence on the higher-scale history. The construction
starts from its stated two-body Hamiltonian. Treating arbitrary input kernels
with comparable bounds still requires checking the endpoint estimates in
that proof, rather than citing its conclusion as an abstract theorem.

**[read]** Gentile and Mastropietro,
[*Renormalization Group for One-Dimensional Fermions. A Review on Mathematical
Results*](https://www.mat.uniroma3.it/users/gentile/ricerca/lavori/gm5.pdf),
Theorem 1 of Section 8 and (8.76)--(8.80), state the conditional tree bounds.
Section 9.1 identifies contributions containing initial irrelevant endpoints
as suppressed at low scales. Section 10, especially (10.8)--(10.11), separates
the Luttinger beta function from lattice corrections and describes short
memory. Section 18 treats the retarded quantum-phonon interaction. These
passages support a precise extension strategy; they do not state the desired
arbitrary-kernel theorem or identify a universal thermal coefficient.

Consequently, exponential retardation itself is already treated in a primary
rigorous construction. The unanswered question is the uniform renormalized
tree estimate for the demonstrator's more general initial endpoints.

## 3. An exact charge-conjugation-equivariant spectral embedding

**[proved]** Write A for the five-qubit cell Hamiltonian and

    Gamma = Z_pL Z_sL Z_c Z_sR Z_pR,
    Q = (Z_pL-Z_sL+Z_c-Z_sR+Z_pR)/2,
    C = X_pL X_sL X_c X_sR X_pR.

Direct substitution in each port check gives

    [A,Gamma]=[A,Q]=[A,C]=0,
    C Q C = -Q,     C Gamma C = -Gamma,     C^2=1.

For example C sends `Phi-` to `-Phi-`, interchanges `01,10`, and
interchanges the logical raising and lowering operators, so the two signed
off-diagonal terms in the check are interchanged. On the code, C interchanges
the logical states `G|0>,G|1>`. The gap is `Delta=D-|kap|>0`.

The even parity space has dimension sixteen, with one code vector. Choose
fifteen orthonormal simultaneous A,Q eigenvectors v_s in its excited part,
with energies `E_s>=Delta`. Their partners `C v_s` form the odd excited
space. Put `q_s=Q(v_s)-1/2`. This is an even integer: if the five-qubit
occupation numbers are n_i, then

    Q-1/2 = -n_pL+n_sL-n_c+n_sR-n_pR,

whose parity equals the total occupation parity.

For each s introduce four distinct auxiliary modes `a_s1,a_s2,b_s1,b_s2`,
all of energy `E_s/2`. With one logical mode f, define an isometry J by

    J G|0> = |0>_f Omega,
    J G|1> = f^dagger Omega,
    J v_s = a_s1^dagger a_s2^dagger Omega,
    J C v_s = f^dagger b_s1^dagger b_s2^dagger Omega.

There are sixty massive modes and one logical mode. All displayed vectors
are orthonormal. J is even and intertwines A with the auxiliary number
Hamiltonian H_mass. Every auxiliary energy is at least `Delta/2`.

Assign auxiliary U(1) charges

    (q_a1,q_a2,q_b1,q_b2)=(1,q_s-1,-1,1-q_s),

and define

    Q_tilde=1/2-n_f+sum_aux q_aux n_aux.

These charges are odd integers. The two encoded excited vectors have charges
`1/2+q_s` and `-1/2-q_s`, respectively. Thus `Q_tilde J=J Q`.

There is a canonical CAR automorphism

    C_tilde: f <-> f^dagger,    a_sr <-> b_sr.

To verify its implementation and phases, let U_swap be second quantization
of the permutation interchanging a_sr and b_sr, chosen to fix the vacuum,
and let Gamma_aux be total auxiliary parity. Then

    U_C=(f+f^dagger) Gamma_aux U_swap

is a self-adjoint unitary with square one. Conjugation sends f to f^dagger.
The two signs from conjugating an auxiliary annihilator by Gamma_aux and
by `f+f^dagger` cancel, leaving precisely its partner. On the encoded pair
states, U_C sends `a1^dagger a2^dagger Omega` to
`f^dagger b1^dagger b2^dagger Omega` with positive sign. It interchanges
the code vectors and reverses the same steps on their partners. Therefore

    U_C J=J C,    [U_C,H_mass]=0,    U_C Q_tilde U_C=-Q_tilde.

The valid subspace and its complement are invariant. Embedded interactions,
including their zero extension on the invalid subspace, preserve these
symmetries whenever the physical interactions do. A symmetry-equivariant
local SW construction can be used: its spectral projections, reduced
resolvents and commutator operations all commute with the physical symmetry
actions. This qualification must be checked for a different chosen transform.

Reflection can also be retained. The cell reflection swapping the two ports
and their partners commutes with A,Q,Gamma,C. Choose the v_s to diagonalize
it as well, with signs r_s. Assign reflection signs `(r_s,1)` to each
auxiliary pair a and b. Their product is r_s on both encoded excited
vectors, and the code reflection is the identity.

After Jordan--Wigner ordering by cells, physical charge conjugation acts on
the logical lattice fermion by

    f_j -> (-1)^(j-1) f_j^dagger.

Each auxiliary annihilator is sent to the paired auxiliary annihilator with
the same staggering. The sign comes from flipping each preceding cell's
parity. The on-site massive covariance is invariant under this operation.
Thus integrating out auxiliary variables preserves a *linear* particle-hole
action on the logical variables. The original forty-five-mode encoding need
not have this useful property; replacing it by this finite encoding costs
constants, not physical dimension or an auxiliary gap.

## 4. Equal-time ordering: the relevant counterterm requires care

**[proved algebraic qualification]** Particle-hole invariance of an operator
does not imply that its normal-ordered polynomial has no one-body term:

    Z_j Z_(j+1)=4 n_j n_(j+1)-2 n_j-2 n_(j+1)+1.

Ignoring those contact terms would give a false kernel symmetry.

A symmetry-preserving scheme can instead be formulated in antisymmetrized
CAR, or fermionic Weyl, coordinates. In a finite Majorana algebra, map an
exterior monomial to the fully antisymmetrized operator product. This is an
invertible linear map: the ordered Majorana products form a basis of the
matrix algebra. It intertwines every linear canonical transformation,
because antisymmetrization commutes with applying the same linear map to
each factor. In particular the symbol of `n_f-1/2` is the corresponding
quadratic Grassmann monomial, and particle-hole conjugation acts linearly
without adding a constant.

This symbol convention has an exact Euclidean realization. At equal-time
insertions use fully antisymmetrized time ordering, and for a Gaussian state
use half the two one-sided covariance limits. Wick's identity for such
insertions is the Pfaffian of the antisymmetrized two-point matrix: it follows
by antisymmetrizing the ordinary Gaussian Wick identity, term by term.
Expand the finite-dimensional Duhamel formula with this ordering at every
instantaneous vertex. It reproduces the operator Weyl polynomial, including
all contact corrections. Unequal-time insertions retain ordinary time
ordering. Hence it gives the same partition function and correlators, not a
new Hamiltonian. Conversion between this convention and normal ordering has
a coefficient cost at most geometric in the number of modes: expand the
finitely many local factors `n_i=(n_i-1/2)+1/2`. This is absorbed by increasing
the support weight in the finite-cell theorem. The zero-temperature
auxiliary-vacuum cancellation is an operator identity and survives the change
of symbols; it must not be re-proved by asserting each Weyl monomial separately
annihilates that vacuum.

Choose logical momentum/frequency cutoffs invariant under
`(k0,k)->(-k0,pi-k)`, and under spatial reflection. For the finite spatial
regulator, use even logical ring lengths with periodic or antiperiodic
fermion boundary condition, or another twist whose grid is explicitly
invariant under this transformation. A generic odd ring is not suitable
for this symmetry argument. The thermodynamic limit can be taken through
these even rings; the ultraviolet norm estimates below hold on all rings.
The finite-mode Grassmann integral preserves the linear symmetry by change
of variables. Its quadratic kernel on the allowed momenta satisfies

    K_h(k0,k)=-K_h(-k0,pi-k).

On a finite torus choose the constant localization as a symmetric average
over allowed frequency/momentum pairs related by this transformation;
the average is zero. In the zero-temperature, infinite-volume Fourier
description at a fixed infrared cutoff, analyticity near either
`k=+pi/2` or `-pi/2` gives `K_h(0,k)=0` directly. Thus the localized
chemical-potential constant is identically zero in this scheme. This is
a statement about each regularized Wilsonian kernel, not a proof of
convergence as h tends to minus infinity. A normal-ordering scheme may
instead carry a nonzero Hartree subtraction.

**[limitation]** The exact Green function obeys its corresponding
particle-hole relation in the symmetric Gibbs state. That identity alone
does not imply gaplessness. A gapped phase can have a Green-function zero
and a singular inverse Green function at the symmetric momentum. Nor does
symmetry alone prove uniqueness or exclude spontaneous symmetry breaking.

## 5. Strengthening the matching norm to a connecting-path norm

**[proved weighted extension of the existing matching argument]** The
[finite-cell theorem](../../direct_attack_20261028/infrared/FINITE_CELL_MATCHING.md)
uses exponential spatial and temporal diameters. Its same spanning-tree
proof gives the stronger external connecting-path weight used in the RG
papers, after a fixed reduction of decay rates. Independently,
[TREE_NORM_EXTENSION.md](../../direct_attack_20261028/infrared/TREE_NORM_EXTENSION.md)
obtains a terminal spanning-tree norm from the occupied network in the
operator-polymer construction, including its thermal periodization estimates.
Here the network instead comes directly from the finite-auxiliary Gaussian
contraction tree. The factor-two geometry is the same; the convergence
majorants are distinct.

For a connected auxiliary contraction tree, connect all legs of each input
vertex along its spatial support interval, at its insertion time. The total
horizontal length is at most the sum of the support interval lengths. Each
tree contraction joins two legs at the same cell, so join them vertically;
the total vertical length is the sum of the tree-edge time distances. The
union is a connected network containing every external point. Doubling its
edges gives a tour visiting every external point. Triangle-inequality
shortcuts remove internal points. Therefore the shortest external visiting
path d_(a,kappa) obeys in the metric
`a|x-y|+kappa|t-s|`,

    d_(a,kappa) <= 2a sum_vertices |X_v|
                  + 2kappa sum_tree_edges |t_e-t'_e|.

The same proof holds with circle distances at finite volume and temperature.
Thus an output weight `exp(d_(a,kappa))` is paid for by input spatial
weights with rate `2a` and tree covariance time weights with rate `2kappa`.
Take `2kappa<g`, where g is the auxiliary gap, and use the weighted covariance
integral `2/(g-2kappa)` in the determinant theorem. All remaining combinatorics
are exactly those of the existing connected expansion. The result is uniform
in external degree: it has no factor proportional to that degree in the
exponent. For point splitting within a time interval of width epsilon,
include a vertical interval of that width at every cell in a vertex's
support. The extra network length is at most
`epsilon sum_vertices |X_v|`, independently of field degree; a small reserve
in the input spatial weight absorbs it uniformly before taking epsilon to
zero. A shortest connecting tree is no longer than this path, so its
weight is controlled as well.

This supplies a stronger norm, not a new infrared estimate. Together with
the finite-order support/degree lemma in
[DEGREE_REPARAMETERIZATION.md](../../direct_attack_20261028/vacancy/DEGREE_REPARAMETERIZATION.md),
it supplies the desired initial degree powers and exponential moments after
rescaling time. It does not claim norm analyticity in the physical lambda
through zero after that singular time rescaling.

### 5.1 An explicit uniform ultraviolet initialization lemma

**[proved using the published determinant theorem]** In hopping units take
the free dispersion `e(k)=2 cos k`, with an inessential overall sign change
if needed. Define

    g(k0,k)=1/(-i k0+e(k)),
    chi(k0,k)=exp[-k0^2-e(k)^2],
    g_IR=chi g,       g_UV=(1-chi)g.

This is the Gaussian decomposition used in the Holstein paper's (2.8),
with the present dispersion normalization. The argument below supplies
explicit uniform bounds suited to arbitrary small initial endpoints.

The covariance g_UV has a volume- and temperature-independent determinant
bound, and its space-time kernel has an integrable exponential weight at
any fixed spatial and temporal rate. Thus it can integrate arbitrary small
retarded input vertices by the effective-action theorem, not just
density-density vertices.

For decay, use the exact cancellation

    g_UV(k0,k)=(i k0+e(k)) integral_0^1
                             exp[-s(k0^2+e(k)^2)] ds.                 (3)

Let `q_s(t)=(4 pi s)^(-1/2) exp[-t^2/(4s)]`. With the Fourier transform
convention fixed, inverse transformation replaces i k0 by either plus or
minus the time derivative; only its absolute value is needed. For
`a>a'>=0`, shifting the spatial contour by `+/- i a` gives

    |g_UV(t,x)| <= exp[-a|x|] integral_0^1 exp[s E_a^2]
                         (|q'_s(t)|+E_a q_s(t)) ds,
    E_a=2 cosh a.

Indeed e is entire and bounded by E_a on that contour, and
`|exp(-s e^2)|<=exp(s E_a^2)`. The cancellation in (3) removes the apparent
Fermi denominator before performing the contour shift. The elementary
inequality `kappa |t|<=2 kappa^2 s+t^2/(8s)` yields

    integral_R exp[kappa|t|] |q'_s(t)| dt
       <= 2 exp[2 kappa^2 s]/sqrt(pi s),
    integral_R exp[kappa|t|] q_s(t) dt
       <= sqrt(2) exp[2 kappa^2 s].

Consequently the full weighted row integral is at most

    [1+2/(exp(a-a')-1)] exp[E_a^2+2 kappa^2]
                             [4/sqrt(pi)+sqrt(2) E_a],                (4)

a finite constant. If the decay norm integrates a doubled time interval,
twice this constant is safe. Spatial and antiperiodic temporal
periodization preserve this bound with circular distances: for every lift
the circular distance is no larger, so the triangle inequality bounds the
periodized absolute sum by the corresponding line sum. This proves the
bound uniformly on finite rings and at every inverse temperature.

For the determinant bound, the full lattice covariance has bound delta=2
with normalized Brillouin-zone measure, by Pedra--Salmhofer Theorem 2.4.
The IR covariance is its convolution in time by q_1, and in space by the
kernel b with Fourier multiplier `exp[-e(k)^2]`. This spatial kernel has
absolute row sum at most one. To see this without a special-function
identity, let S be the unit lattice shift. Since `e=S+S^(-1)`,

    exp[-e^2]=exp[-2] exp[-S^2] exp[-S^(-2)].

The two exponential series have absolute convolution norms at most exp(1)
each. Their product with exp(-2) therefore has norm at most one, on both
the line and a ring. The time convolution has absolute integral one;
its antiperiodic periodization has absolute integral at most one.

Expand an n-by-n determinant of `g-g_IR` by multilinearity in its rows.
There are 2^n choices. A row chosen from g_IR is an integral of a row of g
with shifted space and time, against a kernel of total variation at most
one. After interchanging these integrals with the finite determinant, each
integrand is a determinant of the original covariance, bounded by `2^(2n)`.
This also applies with the auxiliary Gram factors in the definition of a
determinant bound, since those factors are independent of the coordinate
shifts. Hence the determinant bound

    delta_UV=2 sqrt(2)                                               (5)

is valid uniformly. The half-jump equal-time convention is the limit of
averaging the two one-sided row times and preserves the same bound.

**[read]** Pedra and Salmhofer,
[*Determinant Bounds and the Matsubara UV Problem of Many-Fermion
Systems*](https://www.wias-berlin.de/people/koenig/www/FG/PedraSalm.pdf),
Theorem 2.4 supplies the full-covariance determinant bound. Theorem 4.5
accepts arbitrary coefficient-summable Grassmann vertices once determinant
and row/column decay bounds are finite. It supplies an analytic
effective-action map with field-weight loss `2 delta_UV` and geometric
remainder estimates. Their Section 5 also treats ultraviolet integration;
the explicit analytic cutoff and exponential estimate (3)--(5) above are
the application used here.

To retain a connecting-path weight for retarded input vertices, connect
each vertex's leg coordinates by its input path, and connect vertices by
the contraction spanning tree. The network contains all output terminals.
Doubling it bounds the output path by twice the sum of input path lengths
and tree-edge lengths. Reserve this factor two in the input and covariance
rates, apply (4) to tree edges, and keep (5) for the remaining determinant.
Theorem 4.5's other estimates are unchanged. Instantaneous kernel measures
are obtained by point-split approximation with the same norm bound.

More explicitly, let alpha_UV be twice the row-integral bound (4), using
spatial and temporal rates twice the desired output rates, and set
`omega_UV=2 alpha_UV/delta_UV^2`. Measure the input nonconstant vertices
with the doubled connecting-path weight and field weight
`H'=H+2 delta_UV`. If this norm is v and `omega_UV v<=1/2`, then the
output with the desired weaker path weight obeys

    ||W(V,g_UV)||_H <= v+omega_UV v^2/(1-omega_UV v) <= 2v.

Scalar coefficients require the separate rooted density convention already
used in the matching theorem. They are not included in this seminorm.

For completeness the endpoint degree hierarchy also survives this step.
Suppose each initial degree-p vertex, for p>=2 even, vanishes below order
`d_p=max(2,p/2-1)` in an artificial analytic parameter z, with a convergent
field-weighted norm on a fixed disk. A connected term with n vertices of
degrees p_i and l output fields has at least n-1 contractions, so

    l <= sum_i p_i-2(n-1),
    sum_i d_(p_i) >= max(2n,l/2-1) >= d_l.

Thus the output degree-l kernel has no Taylor coefficient below d_l.
The analytic bound on a smaller fixed z-disk and Cauchy's formula then
give `C^l |z|^(d_l)` coefficient bounds after a further fixed field-weight
loss. The degree reparameterization cited above constructs such a family
at any fixed microscopic lambda. Linear particle-hole symmetry preserves
field degree in the Weyl convention, so that reparameterization respects
it. No analyticity in the physically rescaled lambda at zero is asserted.

Both covariance pieces preserve particle-hole and reflection symmetry,
since chi is invariant under `(k0,k)->(-k0,pi-k)`. This finishes a uniform
ultraviolet initialization for the stated class of small initial kernels.
It does not perform the subsequent infinitely many infrared integrations.

### 5.2 Finite-range matrix dispersions: the same initialization

**[proved extension of Section 5.1]** Let there be F finite internal fermion
components and let the free one-particle operator be the finite-range
Hermitian convolution

    h(k)=sum_(r=-R)^R h_r exp(i k r),    h_(-r)=h_r^dagger.

The argument requires no nondegenerate bands, no separated eigenvalues,
and no gap at the Fermi level. Define the fixed finite constants

    A_0=max_i sum_(r,j) |(h_r)_(ij)|,
    E_a=sum_r ||h_r||_op exp(a |r|),
    C_B=exp(A_0^2).

Hermiticity gives the same bound A_0 for column sums. The convolution
row-sum norm is submultiplicative. Therefore the matrix convolution kernel
`B=exp(-h^2)` obeys

    max_a sum_(r,b) |B_(ab)(r)| <= exp(A_0^2)=C_B.                     (6)

This follows directly by its convergent exponential series and
`||h^2||_row<=A_0^2`; unlike the special scalar cosine dispersion, no bound
by one is asserted. The same statement holds on every finite ring, with
either periodic or twisted boundary condition, by periodizing coefficients
and using absolute values.

Take matrix covariances

    g=(-i k0 I+h(k))^(-1),
    g_IR=exp[-k0^2 I-h(k)^2] g,
    g_UV=g-g_IR.

The full covariance has the safe determinant bound `delta_full=2 sqrt(F)`.
For an explicit application of the published proof, measurably diagonalize
the Hermitian h(k) and use the spectral measure consisting of normalized
Brillouin measure and counting measure on its F bands. The components of
each unitary eigenvector have absolute value at most one. The Gram vectors
in the chronological representation of Pedra--Salmhofer's Theorem 2.4
therefore have squared norms at most F, and the same chronological
determinant argument gives this bound, including the additional Gram
factors in its definition. Finite-volume Fourier sums give the identical
estimate. No smooth eigenvector choice is used or needed.

The IR covariance is time convolution by q_1 followed by left spatial and
component convolution with B. Expanding an n-row determinant of `g-g_IR`
as before, every replaced row has a coordinate-and-component integral with
total variation at most C_B. Summing which rows were replaced gives
`(1+C_B)^n delta_full^(2n)`. Thus

    delta_UV = 2 sqrt[F(1+exp(A_0^2))]                                (7)

is a uniform determinant bound. This argument uses the original full
covariance after row shifts and changes of component; it introduces no
volume factor.

The heat-kernel cancellation is a matrix identity, since all its factors
are functions of the same matrix:

    g_UV=(i k0 I+h(k)) integral_0^1
                            exp[-s(k0^2 I+h(k)^2)] ds.                (8)

On `|Im k|<=a`, `||h(k)||_op<=E_a`, and the exponential series gives
`||exp[-s h(k)^2]||_op<=exp(s E_a^2)`. Repeating the contour and Gaussian
integrals in Section 5.1 bounds the operator norm of each spatial kernel.
Summing its at most F entries in a fixed row or column is safely bounded
by F times that operator norm. With `a>a'>=0`, a safe decay constant is

    alpha_UV = 2F [1+2/(exp(a-a')-1)] exp[E_a^2+2 kappa^2]
                                 [4/sqrt(pi)+sqrt(2) E_a].            (9)

The factor two covers the doubled time-interval convention. Periodization
again preserves the bound with circular distances, uniformly in ring
length and inverse temperature. Matrix band crossings do not affect this
proof: the entire matrix exponential, not analytically continued
individual eigenvectors or eigenvalues, is used for the contour estimate.

With delta_UV and alpha_UV as in (7)--(9), use
`omega_UV=2 alpha_UV/delta_UV^2`, field-weight loss `2 delta_UV`, doubled
input path rates, and `omega_UV v<=1/2` exactly as in Section 5.1. The same
effective-action estimate, analyticity and endpoint-degree hierarchy
follow for arbitrary small retarded kernels with these finitely many
components. Constants can grow very badly with F,R and h_r; all are fixed
before choosing the interaction strength. The result is consequently
available for a fixed finite multichannel construction such as n=23,
once its specified free covariance is put in the displayed finite-range
matrix form. It is not uniform over an unbounded number of internal
components. The band-localization and infrared GN estimates remain separate.

**[cutoff qualification]** The Gaussian g_IR still has noncompact momentum
support. It is not already an exactly supported sum of Fermi-branch patches.
In particular, a nonresonant vertex cannot be declared identically zero
solely from this Gaussian cutoff. After the actual finitely many Fermi
points, band projectors and patch separation have been checked, one may
perform a further fixed crossover to compact branch patches. The
complementary covariance has Gaussian frequency decay and a uniform band
gap on its support. Its Fourier absolute integral is finite, giving a
Gram bound by the usual square-root factorization; smoothness on that
gapped support gives summable space-time decay with any prescribed finite
number of moments by integration by parts. The finite-volume Gram bound
uses the same gap and a uniformly bounded Gaussian Matsubara sum; the
decay bounds pass to finite volume and temperature by periodization.
The effective-action
theorem therefore permits this additional finite integration for small
input vertices. Its constants depend on the chosen patch gap and cutoff.

A compact smooth cutoff need not preserve exponential spatial locality:
one must use the resulting fixed-moment norms, or a specified Gevrey
cutoff's subexponential norm, consistently in the subsequent GN proof.
The exact support exclusions in Section 6 refer to a branch decomposition
with these properties established, not to the Gaussian tails before
this crossover. This note does not substitute such a crossover for
checking the actual multichannel Fermi geometry.

### 5.3 A polynomial tree norm compatible with compact patches

**[proved norm conversion and finite-step estimate; conditional endpoint
application]** A fixed global diameter moment is not a substitute for the
connecting-tree norm at unbounded degree. There is, however, a stronger
polynomial *tree* norm that resolves the smooth-cutoff mismatch without
requiring exponential decay of a compactly supported Fourier cutoff.

Let `r(z,z')=|x-x'|+|t-t'|`, with circular distances on a torus, and set

    d_q(z,z')=q log[1+r(z,z')],
    T_q(P)=min_(trees T on P) sum_(ij in T) d_q(z_i,z_j),
    w_q(P)=exp[T_q(P)]
          =min_T product_(ij in T) [1+r(z_i,z_j)]^q.                 (10)

The function d_q is a metric for q>0: the triangle inequality for r and
`1+u+v<=(1+u)(1+v)` give its triangle inequality. Repeated external points
cost no distance. Use w_q in the same pinned coefficient norm as before.

The connected-network argument works in any metric. Take a minimizing
terminal tree inside each input vertex, add the covariance spanning-tree
edges, double the resulting connected network and shortcut it to an
external visiting walk. Removing cycles gives an external terminal tree.
It follows that

    w_q(P_external)
       <= product_vertices w_(2q)(P_vertex)
          product_covariance_tree_edges [1+r(z_e,z'_e)]^(2q).        (11)

Thus input w_(2q) norms and the covariance's fixed 2q moment suffice for
the same finite-step determinant theorem with output w_q norm. No degree
factor multiplies q. A complementary compact-patch covariance with the
gap and smoothness described above has this moment for any fixed q.
The crossover therefore produces small retarded initial kernels in a
polynomial tree norm, not merely a global polynomial-diameter norm.

The exponentially weighted matching output supplies the required input
norm. For every eta>0,

    C_(q,eta)=sup_(r>=0) (1+r)^q exp(-eta r)<infinity.

Choose a terminal tree minimizing the available positive exponential
space-time metric. If its smallest coordinate rate is eta, then on each
edge `(1+r)^q<=C_(q,eta) exp[eta r]`, and eta r is no greater than that
edge's exponential metric length. For p external fields,

    w_q(P)<=C_(q,eta)^(p-1) exp[T_exp(P)].                            (12)

A fixed loss in field weight absorbs this geometric degree factor. This
works at any positive available exponential rates; the temporal rate need
not exceed q. If the available output controls a visiting-path length
rather than a terminal-tree length, it also controls the latter.

Finally suppose the GN endpoint routing lemma being used has the following
explicit bounded-load property: for a chosen terminal tree on p endpoint
legs, its localization differences expand into at most C^p terms; in each
term no edge carries more than b coordinate factors or Taylor derivatives.
Choose a minimizing tree in (10). Each term's endpoint moment factor is
bounded by

    product_edges r_e^(m_e)
       <=product_edges (1+r_e)^b<=w_q(P),    m_e<=b, q>=b.             (13)

Hence the same endpoint integration bound follows from the w_q norm, with
only the routing factor C^p and a fixed field-weight loss. In particular,
loads at most two or four are paid by any fixed q>=4. If another portion
of the GN proof needs additional moments, choose q larger once, before
fixing the ultraviolet constants. This statement assumes the bounded-load
and geometric-count routing lemma; it does not prove that lemma or the
full renormalized-tree induction. It does identify a norm in which the
smooth compact crossover and such an endpoint lemma are compatible.

## 6. What the symmetries leave in the localized action

**[proved algebraic classification]** Write the slowly varying fields near
`+pi/2,-pi/2` as psi_+,psi_-. Number conservation gives equal numbers of
creation and annihilation fields. Translation by one cell excludes the
off-diagonal quadratic mass, whose momentum transfer is pi. It does permit
half-filled umklapp, whose transfer is 2pi. Particle-hole symmetry, with
the ordering qualification above, kills the diagonal constant quadratic
term. Reflection relates the two velocities and wavefunction coefficients.
The remaining marginal quadratic terms are the common time derivative and
the opposite-sign spatial derivatives on the two branches.

At one space-time point, antisymmetry leaves exactly one nonzero
number-conserving quartic monomial,

    psi_+^+ psi_-^+ psi_-^- psi_+^-.

A same-branch contact quartic vanishes. The apparently resonant umklapp
monomial `psi_+^+ psi_+^+ psi_-^- psi_-^-` also vanishes. For its nonlocal
kernel, antisymmetrizing the two identical creation legs forces one
difference; antisymmetrizing the two identical annihilation legs forces a
second. Taylor expansion turns these into two field derivatives, with
coefficients bounded by the corresponding two kernel moments. Its free
scaling dimension is therefore at least four. General derivative quartics
have dimension at least three, and six or more fermion fields have dimension
at least three even before Pauli improvements. Retardation only introduces
time derivatives into this classification; the exponential moments bound
their coefficients.

Consequently the weak one-channel problem has one forward quartic coupling
and a velocity correction, besides wavefunction renormalization. There is
no relevant-counterterm fixed point to solve in a symmetry-preserving scheme.
This classification alone does not control their infinite RG evolution.

## 7. A complete conditional marginal-flow lemma

**[proved implication, with explicit unproved input for this model]** Fix
`gamma>1` and `0<alpha<theta`. Let mu_h be the finite vector of marginal
couplings, `h=0,-1,...`, satisfying `mu_(h-1)=mu_h+b_h`. Suppose that whenever
all already constructed couplings have norm at most epsilon, the normalized
beta function has the bound

    |b_h| <= A epsilon^2 gamma^(theta h)
       + B epsilon sum_(k=h+1)^0 gamma^(theta(h-k)) |mu_k-mu_h|.       (1)

The first term is the summed initial-irrelevant and propagator-curvature
correction plus any finite-ultraviolet remainder of the reference beta
function. The second is the history dependence after subtracting its
constant-history value. In the published construction these bounds come
from renormalized-tree short memory and the vanishing Luttinger beta
function. Equation (1) must be proved for the actual general initial
endpoints before applying this lemma to the demonstrator.
Here epsilon is a provisional smallness parameter, chosen as in the degree
reparameterization so that every initial irrelevant endpoint is O(epsilon^2).
It is not automatically the physical quartic coupling. A linear velocity
correction in the actual endpoint size is then O(epsilon^2); no cancellation
of all first-order corrections is assumed. The beta functions must be
defined on the stated small-history domain, as part of the conditional
tree theorem.

Put

    C_(alpha,theta) = 1 / [(1-gamma^(-theta))
                           (gamma^(theta-alpha)-1)],
    K=2A.

Assume `|mu_0|<=epsilon/2` and

    B epsilon C_(alpha,theta)<=1/2,
    K epsilon/(1-gamma^(-alpha))<=1/2.

Then all scales exist within this small-coupling domain and

    |mu_(h-1)-mu_h| <= K epsilon^2 gamma^(alpha h).                    (2)

Proof: proceed downward in h. At h=0 the history sum is empty. At a later
step, previously proved increments give

    |mu_(h+n)-mu_h|
       <=K epsilon^2 gamma^(alpha h) sum_(r=1)^n gamma^(alpha r).

Multiplying by `gamma^(-theta n)` and summing n yields

    sum_(n>=1) gamma^(-theta n) sum_(r=1)^n gamma^(alpha r)
      = C_(alpha,theta).

Since `gamma^(theta h)<=gamma^(alpha h)` for h<=0, (1) is at most
`(A+B K epsilon C_(alpha,theta)) epsilon^2 gamma^(alpha h)`, which is
bounded by the right side of (2). The sum of all previous increments is
at most `K epsilon^2/(1-gamma^(-alpha))<=epsilon/2`, so the running vector
stays in the required domain. This also makes the continuation induction
noncircular. Summing (2) proves existence of mu_infinity and

    |mu_h-mu_infinity|
        <= K epsilon^2 gamma^(alpha h)/(1-gamma^(-alpha)).

This elementary implication is useful here because exact half filling
removes the separate relevant-direction fixed-point construction. It is
not a replacement for the renormalized-tree estimate that establishes (1).

## 8. The remaining decisive work

The demonstrated reductions supply a finite symmetric embedding, a genuine
connecting-path norm, uniform ultraviolet initialization, controllable
endpoint degree powers, no relevant
constant in a symmetry-preserving scheme, and exactly the single-channel
marginal operator list. The missing RG statement is now specific: prove
the conditional tree majorant and its history-continuity estimate for those
initial kernels, uniformly in volume, inverse temperature and infrared
scale. The source papers prove them for kernels arising from their stated
microscopic interactions. Their proofs strongly indicate the extension,
but this note does not pass a proposed extension off as a verified theorem.

Even a completed marginal-flow estimate and fermion/current asymptotics
would not by themselves prove that the physical chain has central charge
one in the conjecture's required sense. A stress-tensor or thermal scalar
identification, and the relevant source matching, would remain. Sector
selection is a separate issue as well. No numerical anisotropy, power-law
two-point function, or formal bosonization is used here as a substitute.
