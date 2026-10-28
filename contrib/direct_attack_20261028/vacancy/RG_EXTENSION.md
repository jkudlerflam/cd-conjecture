# From finite-cell matching to multichannel RG: an initialization lemma

Status: **[proved initialization estimates; full RG and CFT extension
unproved]**. This note applies to the fixed finite-cell setting of
[FINITE_CELL_MATCHING.md](../infrared/FINITE_CELL_MATCHING.md), conditional
on that theorem's hypotheses and norm estimates. It does not remove the
unbounded-cell issue in the mobile wrapper or prove the vacancy comparison.

Author: OpenAI GPT-6-family Codex vacancy worker, 2026-10-28. AI-assisted;
no human review. All estimates below are derived here. The source audit
was performed by opening the primary paper.

The follow-up [DEGREE_REPARAMETERIZATION.md](DEGREE_REPARAMETERIZATION.md)
uses the actual microscopic relation between perturbative order and field
degree to recover the stronger endpoint powers appearing in the published
GN estimates. This narrows the extension problem beyond the generic
kernel-norm assumptions used below; it still does not claim a completed
multichannel RG or CFT theorem.

## 1. Precisely what the primary source supplies

**[read]** Porta, Scola and Singh,
[arXiv:2509.08665v1](https://arxiv.org/html/2509.08665v1), assume finite-range
quadratic hopping, weak density interactions, distinct transverse Fermi
crossings and elastic quartic scattering (Assumption 2.1, equation 2.18).
Their ultraviolet integration already generates general nonlocal kernels
(4.27–4.32). The GN estimates (4.111–4.119) use Proposition 4.5's running
coupling bounds; the proposition depends on the reference model's vanishing
beta function. Theorem 2.6 concerns density/current response. Proposition
6.2 compares specified density-current correlations with reference-model
correlations. These statements are not a theorem for arbitrary initial
retarded vertices or a central-charge theorem.

**[inference from the proof organization]** General kernel shape is not
the decisive new difficulty. The required extension is stability of the
bounded marginal flow and the counterterm construction under a small
analytic family of additional endpoints. Using the GN bounds to assert
that flow would reverse their logical dependence.

## 2. Assumed output of finite-cell matching

Work at zero temperature in cell coordinates. Let the quadratic logical
bandwidth be `W>0`; its Fermi velocities divided by `W` are fixed and
nonzero. Suppose the induced, translation-invariant, even action is

    V_K(psi)=sum_(p even, p>=2) integral dZ_1 ... dZ_p
                          K_p(Z_1,...,Z_p) psi(Z_1)...psi(Z_p).

Species and creation/annihilation labels are included in each coordinate.
Fix any consistent factorial convention for antisymmetric kernels; the
same convention is used in all norms. Fix one field coordinate, integrate
the other `p-1` coordinates, sum species, and suppose

    sum_(p>=2) H^p ||K_p||_(a,kappa) <= r,                   (1)

where the norm includes the weight

    exp[a diam_x + kappa diam_t],        a,kappa>0.

The kernel can be a signed measure with finite total variation, allowing
instantaneous delta kernels; the integrations in (1) then mean total
variation integrations. Scalar action densities are recorded separately.

No quartic-degree truncation is assumed. The hypotheses are stronger than
smallness of an operator norm: absolute coefficient summability and the
diameter weights are essential.

## 3. Exact bandwidth rescaling

**[proved]** Set `u=W t` and leave the lattice fields unrescaled. This
preserves the kinetic term `integral dt psi^dagger partial_t psi`. The
transformed `p`-field kernel is

    K'_p(x_1,u_1;...;x_p,u_p)
       =W^-p K_p(x_1,u_1/W;...;x_p,u_p/W).                  (2)

When `p-1` relative times are integrated, the Jacobian is `W^(p-1)`.
Therefore the pinned norm, at every degree, picks up exactly `W^-1`:

    sum_p H^p ||K'_p||_(a,kappa/W) <= r/W =: epsilon.         (3)

The same computation holds for measure kernels by pushforward. In
particular, high field degree does not produce extra inverse powers of
the small bandwidth. A scalar energy density `r0` becomes `r0/W`, since
`beta M r0=(W beta)M(r0/W)`.

For the finite-cell schedule `r<=C lambda^3` and `W=c_n lambda^2`, (3)
gives `epsilon<=C lambda/c_n`. The rescaled memory length is `W/kappa`,
which shrinks with the logical bandwidth. The constants may be very large
for a fixed cell, but no dependence on the number of cells enters (3).

## 4. A fixed field-weight reserve controls all Taylor remainders

**[proved]** Choose `0<h<H`, and write `theta=h/H<1`. Let the last
coordinate be the anchor. For external frequency and momentum offsets
bounded in absolute value by `Lambda`, the Fourier phase of a kernel
has magnitude bounded by

    |phase| <= p Lambda [diam_x + diam_u].

Using separate anisotropic coordinates would improve constants but is
unnecessary. Define the rescaled range

    ell=a^-1+W/kappa.

The exponential moment in (3) implies, for every nonnegative integer `m`,

    integral (diam_x+diam_u)^m |K'_p|
       <=m! ell^m ||K'_p||_(a,kappa/W).                    (4)

One proof expands the power with the binomial theorem and bounds each
factor by its exponential series; the resulting sum is
`m! (a^-1+W/kappa)^m` times the weighted norm.

The integral remainder for `exp(i z)` satisfies

    |exp(i z)-sum_(j=0)^(m-1)(i z)^j/j!|<=|z|^m/m!        (5)

for real `z`. Expand about any fixed collection of Fermi momenta, whose
phases have modulus one. Let `T_(m-1)` denote the Fourier Taylor polynomial
in the external offsets through total degree `m-1`. Equations (3)–(5) give

    sum_p h^p sup_(offsets<=Lambda)
        |Khat'_p-T_(m-1)Khat'_p|
       <= C_(theta,m) epsilon (ell Lambda)^m,              (6)

where

    C_(theta,m)=sup_(p>=2) p^m theta^p < infinity.

All species sums are already present in (3). There is no truncation in
`p`, no field-degree dependent choice of cutoff, and no factorial
divergence. Purely temporal localization gives the sharper replacement
`ell -> W/kappa`.

In particular, quadratic Taylor remainders have two powers of external
scale after their constant and linear terms are extracted. Quartic
remainders have one power after their constant terms are extracted.
These are genuine bounds on the initial kernels, not a claim about
kernels after arbitrarily many RG steps.

## 5. Linear rescaling of the irrelevant initial interaction

**[proved linear-rescaling statement]** In one space and one time
dimension, a local term with `p` fermion fields and `m` derivatives has
engineering dimension

    D=2-p/2-m.

This follows directly by scaling both coordinates by `Lambda^-1` and
each massless fermion field by `Lambda^1/2` in its action. Its dimensionless
coefficient at scale `Lambda` consequently contains the factor `Lambda^-D`.

For the degree-`p>=6` part, `Lambda^(p/2-2)<=Lambda` when `0<Lambda<=1`.
For the quartic remainder from (6), one derivative gives `D=-1`.
For the quadratic remainder, two derivatives also give `D=-1`.
Combining (3) and (6), the normalized initial irrelevant interaction is
bounded by

    C epsilon (1+ell+ell^2) Lambda.                        (7)

The constant can depend on the finite number of Fermi branches and on
their fixed nonzero dimensionless velocities. At a lower fixed starting
scale one can ensure that all branch neighborhoods are disjoint and the
nonresonant quadratic/quartic momentum combinations are absent.

Projection from the original finite orbital labels to finitely many
nondegenerate band branches can multiply a degree-`p` estimate by
`C_band^p`. Reserve field weight for this fixed conversion by requiring
`C_band h<H` instead of just `h<H`. The first two derivatives of the band
eigenvectors are bounded on the fixed disjoint neighborhoods; their
polynomial-in-`p` factors are absorbed by the same strict weight reserve.
No uniformity in an increasing number of bands is asserted.

This proves that arbitrary degree and exponentially short retardation
do not, by themselves, spoil the initial irrelevant-kernel power count.
It is explicitly a statement before further Gaussian contractions.
A high-degree vertex can contract to a quadratic or quartic vertex and
shift its initial local coupling by order `epsilon`. Such contributions
must be retained in the counterterms and running couplings; (7) does
not declare them nonexistent or automatically summable over all scales.

## 6. The theorem that is still required

**[unproved target, stated independently of any citation]** Start with
finitely many simple Fermi crossings, a positive elastic-scattering
margin, small forward quartic couplings, and an additional even,
number-conserving initial interaction satisfying (3). Establish an
infinite RG iteration with the following properties:

1. The relevant quadratic counterterms can be solved at fixed microscopic
   parameters. The resulting Fermi-point shifts stay within the certified
   momentum margin, and every renormalized velocity remains nonzero.
2. The beta function for the forward quartic couplings is summable, with
   arbitrary-degree retarded endpoints included. This requires comparison
   with the vanishing-beta reference theory and a short-memory estimate,
   not just negative dimensions of irrelevant endpoints.
3. A fixed analytic coefficient norm controls the degree sum throughout
   the iteration. Repeatedly spending a fixed fraction of field weight at
   every scale would exhaust the radius; (6) spends it only once during
   initialization.
4. Sources are integrated with the action. Physical density/current
   sources acquire induced kernel corrections during massive matching.
   Exact number conservation must be expressed for the complete source
   functional, including the resulting contact terms.

The initialization estimates here supply useful input for such a theorem.
They do not establish any of its four conclusions.

## 7. Why response does not yet identify c

**[proved logical limitation]** Normalized correlators are invariant if
the generating functional is multiplied by a source-independent factor.
That factor can change the free energy and its `T^2` term. Consequently
matching normalized density/current responses does not, by itself,
determine the thermal central-charge coefficient.

There is a second distinction: a response formula is a statement about
specified observables, whereas a unitary CFT scaling limit includes a
local operator algebra and stress tensor. The existence and normalization
of the full stress sector cannot be inferred merely from the number of
microscopic charge labels.

Two sufficient continuations are available. One may retain all scalar
vacuum terms and prove the pressure comparison to accuracy `o(T^2)` after
subtracting the ground energy. Alternatively one may construct the full
positive-level chiral current algebra in an independently established
unitary conformal scaling limit. The latter gives a Sugawara lower bound
on `c`, as proved separately in C46; a large enough lower bound could
already refute the conjecture without identifying every neutral sector.
Neither continuation is supplied by (3), (6), or (7).

The present concrete result is therefore an initial-kernel norm theorem
and a precise remaining RG statement. It strengthens the route from
finite-cell matching to RG while keeping both the marginal-flow and
central-charge steps open.
