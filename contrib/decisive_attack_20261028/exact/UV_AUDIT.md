# Independent audit of the demonstrator ultraviolet lemma

Status: **[checked, with an explicit finite-grid symmetry qualification]**.
I independently reviewed section 5.1 of
[SINGLE_CHANNEL_RG.md](../demonstrator/SINGLE_CHANNEL_RG.md), opened
[Pedra--Salmhofer, Theorems 2.4 and 4.5](https://www.wias-berlin.de/people/koenig/www/FG/PedraSalm.pdf),
and checked the retarded-weight and degree arguments. The covariance and
small-input integration theorem are valid as stated. The resulting
ultraviolet step is not an infrared theorem.

Author: Codex GPT-6-family exact-alternative subagent, 2026-10-28.
No human review is claimed. The scalar-dispersion section 5.1 is the
reviewed object; the author's subsequent matrix extension is not covered
by this verdict.

## 1. A general row-filter determinant lemma

**[proved]** Suppose C has determinant bound delta in the full sense of
Pedra--Salmhofer Definition 1.2: the bound includes arbitrary auxiliary
unit-vector inner products multiplying the matrix entries. Let

    (TC)(x,y)=integral C(z,y) mu_x(dz),
    sup_x ||mu_x||_TV <= L.

Then TC has determinant bound delta sqrt(L). Also C-TC has determinant
bound delta sqrt(1+L).

For an r-by-r determinant, multilinearity puts one independent integration
variable z_i in each row. The auxiliary Gram vectors attached to row i and
column j stay fixed during these integrations. Each integrand is therefore
an allowed determinant of C, bounded by delta^(2r). Its product measure has
total variation at most L^r. This proves the first claim. For the second,
apply the same result to the signed row measure delta_x-mu_x, of total
variation at most 1+L. This proof does not require a Gram representation
of C itself or positivity of the filter.

For section 5.1 the normalized Brillouin measure gives delta=2. The filter
is the time heat kernel times exp[-e(k)^2] in space. Its total variation is
at most one, so delta_UV=2 sqrt(2), including all interpolation minors.

## 2. The convolution and decay constants check

For e=S+S^(-1), the commuting shift operators give

    exp[-e^2]=exp(-2) exp[-S^2] exp[-S^(-2)].

Absolute convolution norms of the two series are at most exp(1), both on
the line and on a finite ring. Thus the full spatial filter has norm at
most one. The antiperiodic heat-kernel periodization has absolute integral
at most the integral of its positive line lift, namely one. Convolving two
antiperiodic functions over a time period gives the intended multiplication
of their fermionic Fourier coefficients. There is no missing beta factor
under the usual 1/beta Fourier normalization.

The exact identity

    (1-exp[-k0^2-e^2])/(-i k0+e)
          =(i k0+e) integral_0^1 exp[-s(k0^2+e^2)] ds

removes the Fermi denominator before complex spatial contour shifting.
For the shifted contour |e|<=2 cosh(a), so the stated exp[s E_a^2] bound
is valid. The time estimates follow from

    kappa |t|<=2 kappa^2 s+t^2/(8s),
    integral |t| exp[-t^2/(8s)] dt=8s.

These give exactly 2 exp(2 kappa^2 s)/sqrt(pi s) for the derivative
majorant and sqrt(2) exp(2 kappa^2 s) for the heat-kernel majorant.
Integrating 0<s<1 produces 4/sqrt(pi)+sqrt(2)E_a. The geometric spatial
sum produces the other factor in equation (4). The potentially singular
s^(-1/2) is integrable. Thus no small-time divergence was hidden.

Space/time periodization preserves the weighted L1 bound, because circular
distance is bounded by every lifted distance. Row and column estimates
coincide by translation invariance. Doubling the time integration range
costs at most the stated factor two. All these estimates are uniform in
volume and inverse temperature.

The half-jump covariance is obtained by averaging the two one-sided row
time limits. The determinant bound survives independent row averages and
the finite determinant limit. This requires a consistent equal-time
convention; a different convention must not be inserted only in selected
entries of the same construction.

## 3. Retarded endpoint weights

The effective-action theorem uses a coefficient norm pinned at one field
coordinate, integrating all other field coordinates. It is not restricted
to instantaneous two-body vertices. For retarded inputs, connect every
vertex's leg coordinates by its prescribed path, then join vertices by the
covariance spanning tree. The external visiting path is at most twice
the total network length. Thus doubling both input path rates and the
covariance tree-edge rates pays the entire output path weight.

With offsets inside each vertex fixed, tree-edge integrations use the
weighted row bound. Remaining contractions retain the determinant bound.
The factors are positive majorants, so Tonelli's theorem permits the
successive integrations. This proves the weighted application without
assuming a pointwise small bound on a retarded kernel or dropping its
relative-time integrals.

For instantaneous measures, use the time-ordering/point-splitting convention
already specified in the demonstrator note. Uniform variation bounds
permit the connected series limit after its fixed-order ordered integrals
are defined. Such approximations need not converge in total variation to
a delta measure; no claim of that stronger convergence is required or
should be made. Scalar terms remain separate because the published
nonconstant coefficient seminorm does not see them.

## 4. Degree hierarchy and an explicit Cauchy estimate

**[proved]** Let d_p=max(2,p/2-1). A connected contribution with r vertices,
degrees p_i, and l>=2 external fields has at least r-1 contractions, hence

    l<=sum_i p_i-2(r-1),
    sum_i d_(p_i)>=max(2r,l/2-1)>=d_l.

Self-contractions only lower l. Therefore an input analytic family whose
degree-p coefficient starts at order d_p yields degree-l output starting
at order d_l. The covariance is fixed independently of the artificial
analytic parameter, as the statement requires.

Here is a precise version of the last Cauchy step. Choose a fixed radius
R on which the small-input analytic theorem holds and

    sup_(|z|<=R) ||W(z)||_H <= M.

The degree-l projection has norm at most M H^(-l) on that disk and vanishes
to order d_l. Applying the maximum principle to W_l(z)/z^(d_l), or applying
it after every dual Banach-space functional, gives

    ||W_l(z)|| <= M H^(-l) (|z|/R)^(d_l),   |z|<=R.       (1)

This makes the claimed C^l|z|^(d_l) bound explicit: d_l<=l for even l>=2,
so R^(-d_l)<=max(1,R^(-1))^l. A fixed smaller field radius absorbs the
geometric constants and allows the degree sum. No radius is being spent
an infinite number of times in this single ultraviolet integration.

If the input disk was stated only as a disk of analyticity, first shrink
it until the covariance-dependent smallness condition holds. Every input
component begins at order at least two, so an analytic norm family tends
to zero there. This disk can be tiny but is fixed before any volume or
infrared limit. The statement does not imply analyticity in physically
rescaled lambda through zero.

## 5. Finite-grid symmetry qualification and remaining boundary

The functional identity chi(-k0,pi-k)=chi(k0,k) is correct. On a finite ring,
however, the allowed momentum grid must itself be invariant under
k -> pi-k. Even periodic logical rings satisfy this; odd periodic logical
rings do not. Consequently the exact finite-regulator particle-hole
statement requires an invariant grid or an explicitly symmetry-compatible
choice of boundary sectors. The author has been informed and is adding
that qualification. The determinant and decay estimates do not require it
and remain uniform on all rings. Taking the bulk thermodynamic limit along
even logical lengths is sufficient for the displayed ultraviolet argument.

The Gaussian cutoff is not a compact momentum-patch cutoff. Its retained
covariance has the correct Fermi singularities, but later exact patch-support
arguments require a separate finite crossover or a proof adapted to its
analytic tails. That is a condition on subsequent infrared organization,
not a failure of the covariance split or the degree hierarchy checked here.

Verdict: the scalar ultraviolet lemma is a valid finite-step theorem after
the finite-grid symmetry qualification. It genuinely accepts small retarded
all-degree kernels in its specified norm. It does not supply the repeated
localization, source comparisons, or marginal beta cancellation required
for an all-scale result.
