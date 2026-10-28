# A single-scale thermal bound with determinant contractions

Status: **[proved single-scale theorem using the published BBF identity;
full multiscale application still open]**. This supplies the fixed-history
vacuum-root comparison left as an obligation in THERMAL_CRITERION.md.
The theorem allows arbitrary even field degree and retarded endpoint
kernels. It never expands a Gram determinant into Wick pairings.

Provenance: OpenAI GPT-6-family Codex thermal worker, 2026-10-28.
AI-assisted; no human review or formal verification.

## 1. Precise single-scale setting

Space is Z and imaginary time is R, or its circle of circumference beta.
The number of fermion labels is finite and fixed. The only nonzero Gaussian
pairing is between a minus and a plus field; individual even endpoints
need not separately conserve their number. Any unbalanced total monomial
has zero expectation and can be discarded.

Fix an integer p>=3 and a scale 0<s<=1. A translation-invariant line
covariance g_s and its antiperiodization

    g_(s,beta)(t,x)=sum_(m in Z)(-1)^m g_s(t+m beta,x)  (1)

obey, with constants independent of s and beta whenever beta s>=1,

    sup_(t,x) |g_s(t,x)| <= D s,
    integral_R dt sum_x |g_s(t,x)|(1+s|t|)^p <= D_p/s,
    |g_s(t,x)| <= D'_p s [1+s(|t|+|x|)]^(-p-3).       (2)

Label sums can be included in these bounds. Both covariances have Gram
representations with vector norms at most sqrt(D_G s). A tensor product
with any BBF interpolation Gram matrix therefore obeys the same bound.
These are hypotheses on a smooth **single scale**, not a false Gram
claim about an uncut ultraviolet covariance.

These hypotheses hold for finitely many smooth frequency-momentum shells
of transversal Fermi crossings with fixed nonzero velocities, after the
usual field normalization. For example, on a shell of radius comparable
to s, the absolute Fourier integral of
chi_s(k0,k)/(-ik0+v k) is O(s): its support has area O(s^2), its
denominator is bounded below by a constant times s. Taking square roots
of its absolute Fourier amplitude gives the Gram vectors. The Matsubara
version has the same bound when beta s>=1. Smooth cutoff derivatives and
integration by parts give (2). Constants may depend on the fixed cutoff,
number of branches, and smallest velocity.

An endpoint V is specified by even kernels K_l, l>=2, with one field
coordinate anchored. The coefficient convention is fixed once; any 1/l!
can be included in K_l. Define the line coefficient norm

    ||K_l||_p = sum_labels integral |K_l(0,X_2,...,X_l)|
                                      (1+s diam_t X)^p,
    u_s(V)=sum_(l even>=2) H^l s^(l/2-2) ||K_l||_p.     (3)

Spatial relative coordinates are also integrated/summed in (3). Kernels
may be finite signed measures, including instantaneous delta kernels.
The smooth scale covariance removes equal-time ambiguities; a fixed
normal-order point split can also be retained throughout.

The thermal endpoint is exactly the periodization of this **same** line
kernel: integrate its anchor time in [0,beta), integrate all relative
coordinates on the line, and evaluate every field as an antiperiodic
function. This is an unambiguous definition of periodization even for
retarded vertices. A vertex translated by beta is unchanged because l
is even. Scalar endpoints are treated separately and retained exactly;
an exactly periodized scalar density has no thermal difference.

Let

    F_(s,beta)(V)=-(beta |Lambda|)^(-1)
                          log E_(g_(s,beta)) exp(-V)   (4)

after the spatial thermodynamic limit. Let F_(s,infinity) be its line
counterpart, defined by its connected expansion per spacetime volume.
The Gaussian expectations in (4) are normalized. Thus F is the full
**scalar produced by this integration**, not a normalized correlation
function. Any Gaussian partition prefactor or field-rescaling Jacobian
outside (4) must remain in the overall RG bookkeeping.

## 2. The theorem

Choose H sufficiently large, depending only on the covariance constants
and p. There are u_*>0 and C_p<infinity, independent of beta and s, such
that for u_s(V)<=u_* the connected scalar expansion converges absolutely,
is analytic in V, and

    |F_(s,beta)(V)|+|F_(s,infinity)(V)| <= C_p s^2 u_s(V),
    |F_(s,beta)(V)-F_(s,infinity)(V)|
                    <= C_p s^2 u_s(V)(beta s)^(-p).    (5)

If V=V_*+R and u_s(V_*)+u_s(R)<=u_*/4, then

    |[F_(s,beta)(V)-F_(s,beta)(V_*)]
       -[F_(s,infinity)(V)-F_(s,infinity)(V_*)]|
                    <= C'_p s^2 u_s(R)(beta s)^(-p),  (6)
    |F_(s,infinity)(V)-F_(s,infinity)(V_*)|
                    <= C'_p s^2 u_s(R).              (7)

All kernels and covariances in these comparisons use fixed running data.
The theorem does not assert that the kernels actually produced at
temperature T are exact periodizations of the zero-temperature kernels.
Their separate history difference is an application step, stated below.

## 3. Connected expansion without factorial Wick sums

**[read identity]** The BBF formula writes a truncated expectation of m
even field clusters as a sum over spanning trees between those clusters,
with one propagator per tree edge and one determinant for every remaining
contraction. Its interpolation matrix is positive with diagonal entries
one; its interpolation measure has total mass one. This exact identity is
stated in [Porta, Scola and Singh, (4.109)-(4.110)](https://arxiv.org/html/2509.08665v1).
The Gram interpolation property is also explained in
[Pedra and Salmhofer, Section 1, (4)-(7)](https://www.wias-berlin.de/people/koenig/www/FG/PedraSalm.pdf).

Expand log E exp(-V) as its connected cumulant series. At order m, let
l_i be the degree of endpoint i and P=sum_i l_i. For a nonzero term P
is even. After m-1 tree contractions, the remaining determinant has

    r=P/2-m+1 >= 1                                    (8)

rows. Fix all relative endpoint coordinates and use the tree edge
differences as coordinates for the m-1 relative vertex anchors. This
change of variables has unit Jacobian, also on lattice space. Integrating
the absolute tree propagators costs at most (D_p/s)^(m-1), and the
determinant is bounded by (D_G s)^r. Endpoint factors have dimensions
s^(2-l_i/2) according to (3). The combined scale power is exactly

    s^[-(m-1)+r+sum_i(2-l_i/2)] = s^2.                 (9)

The tree sums and field attachments cost exponentially in degree, rather
than factorially. For m>=2 the weighted Cayley identity gives

    sum_(labelled trees T) product_i l_i^(deg_T i)
                   = (product_i l_i) P^(m-2).         (10)

It follows from the Prufer-word bijection: each of the m-2 letters in the
word contributes the weight of its label, and each vertex contributes one
additional l_i. The left side bounds all choices of fields for tree
attachments, even if some choices are inadmissible. Using

    product_i l_i <= 2^P,
    P^(m-2)/(m!) <= exp(P)/[m(m-1)],                   (11)

the 1/m! in the cumulant series pays the tree count. The second inequality
is the term of order m-2 in exp(P). Choose H to absorb (2e)^P, the Gram
degree constants, and a further fixed reserve. Summation over degrees is
then bounded by a convergent geometric series in C u_s(V). For m=1 there
is no tree and the same estimate follows immediately from its determinant.
This proves convergence and the first bound in (5), initially at finite
cutoffs and then by the absolute majorant.

## 4. Exact unfolding along a tree

Apply the same BBF formula on the time circle and periodize the endpoint
kernels as specified above. Expand **only the m-1 tree propagators** using
(1). Root the tree at endpoint 1. The integer image index on each tree
edge determines a unique integer lift of each other vertex anchor.
Summing these image indices and integrating the anchors over [0,beta)
is exactly integration of their lifted anchors over R.

The signs can be checked without guessing a thermal loop convention.
Moving a vertex anchor by n_i beta multiplies each of its l_i fields by
(-1)^(n_i). The tree image signs are absorbed by the complementary
row/column signs of the residual determinant. The total sign is
product_i (-1)^(n_i l_i)=1, because all l_i are even. Internal relative
endpoint times stay on the line throughout. Consequently the thermal
amplitude after unfolding has exactly the same endpoint kernels and line
tree propagators as its zero-temperature counterpart. Only its residual
determinant still uses g_(s,beta) at the lifted field coordinates.

This is why one must retain a determinant: the thermal comparison can be
made after unfolding the small spanning tree, without expanding all
fermionic loops or paying their factorial count.

For m=1 there is no tree to unfold. The endpoint's anchor is fixed by
translation invariance and all its relative coordinates are already on
the line by the endpoint definition. Its thermal expectation is exactly
the residual determinant of g_(s,beta) at those coordinates. Thus the
following comparison applies unchanged, with R_t equal to that single
endpoint's diameter. There is no missing one-vertex loop sign.

## 5. Short configurations: telescope the determinant

For a fixed lifted configuration define

    R_t=sum_(tree edges e) |t_e|
                          +sum_(vertices i) diam_t(endpoint i).       (12)

The time separation of any two fields is at most R_t: follow the unique
tree path between their vertices and include the endpoint displacements
at its vertices. If R_t<=beta/2, (1)-(2) give uniformly in space

    |g_(s,beta)(t,x)-g_s(t,x)|
                         <= C_p s (beta s)^(-p).       (13)

Indeed every nonzero image has |t+n beta|>=(|n|-1/2)beta,
and the image sum converges with the p+3 decay in (2).

Telescope the residual determinant by its columns, replacing thermal
columns one at a time by line columns. Expand only the one difference
column. This produces at most r^2 terms, each with one entry bounded by
(13) and an (r-1)-by-(r-1) mixed cofactor.

The mixed cofactor has a Gram bound. Use the direct sum of the line and
thermal Gram spaces for every row vector; each column selects its own
summand. Row norms grow by at most sqrt(2), column norms do not grow.
The BBF interpolation vectors can be tensored into both representations.
Hence every such cofactor is bounded by

    (sqrt(2) D_G s)^(r-1).                             (14)

Combining (13)-(14) bounds the determinant difference by
C_p r^2 (C s)^r (beta s)^(-p). The polynomial r^2 is paid by the fixed
field-weight reserve in H. No expansion into determinant permutations is
used.

## 6. Long configurations: use one moment

On R_t>beta/2, bound both determinants separately by their Gram bounds.
Use

    1_(R_t>beta/2) <= (2R_t/beta)^p,
    (sum_(j=1)^k a_j)^p <= k^(p-1)sum_(j=1)^k a_j^p.  (15)

There are at most 2m-1 summands in (12). Assign the single p-th moment
in (15) to one tree line or one endpoint. A tree-line moment is paid by
(2); an endpoint moment is paid by (3). Thus the integrated bound is the
ordinary s^2 tree majorant times

    C_p (2m-1)^p (beta s)^(-p).                        (16)

For m=1 this is just the moment of its sole endpoint. The polynomial in
m is summable after reducing the fixed analytic coupling radius, since
sum_m m^p z^m converges for |z|<1. The already-reserved degree weight
handles the additional r^2 factor in the short-region estimate.

The short and long regions exhaust the same lifted integration domain.
Their bounds prove the second line of (5), uniformly before removing
spatial and auxiliary degree cutoffs. Absolute convergence permits their
removal and the spatial thermodynamic limit.

## 7. Extracting a marked reference defect

Let D_beta(V)=F_(s,beta)(V)-F_(s,infinity)(V). The preceding estimates
apply to complex endpoint coefficients and give analyticity on a ball
of radius u_*. For every z0 in [0,1], the disk

    V_*+z0 R+zR,       |z|<=u_*/[4u_s(R)]

lies in a smaller fixed analytic ball under the hypotheses of (6).
Apply Cauchy's derivative bound to D_beta in this disk and integrate
z0 from zero to one. Equation (5) gives (6); the same argument without
the thermal subtraction gives (7). When R=0 both claims are immediate.

This argument has the same effect as requiring a marked defect endpoint
in every contributing tree, but it avoids an additional combinatorial
enumeration. Its bound is linear in the defect even when the reference
interaction V_* is larger. In particular, if u_s(R)<=A omega(s), then
(6) has exactly the form

    C A s^2 omega(s)(T/s)^p                            (17)

required by THERMAL_CRITERION.md. The one-step scalar issue is therefore
closed under the stated coefficient and covariance hypotheses.

## 8. Derivative endpoints and one localization step

**[proved extension]** Smooth-shell localization produces derivatives of
fields, which can be represented directly instead of treating derivative
deltas as finite measures. Give each field a fixed-order time derivative
or lattice difference label a_f, with order at most d0, and put
A=sum_f |a_f| at an endpoint. Suppose the covariance derivative bounds
corresponding to (2) hold: each differentiated field brings one power of
s per derivative. On a smooth shell these follow from its Fourier
multipliers, including the lattice difference multiplier exp(ik)-1 after
the Fermi oscillation is removed. The Gram vector for a differentiated
field has norm at most C s^(|a_f|+1/2). The constant can depend on d0.

Replace the dimension factor in (3) by

    H^l s^(l/2+A-2),                                  (18)

and sum over the finitely many derivative labels per field. The tree and
determinant contractions consume each field exactly once, so their total
extra scale factor is s^A. This cancels the extra s^(-A) of its coefficient
in the dimensionless bound. Equations (9), (13)-(17) therefore hold without
change. Image differentiation commutes with the absolutely convergent
periodization sum. The mixed Gram argument uses the differentiated row
and column vectors. This proves the derivative-endpoint version of the
single-scale theorem.

For one temporal Taylor subtraction, the identity

    psi(t+r)-psi(t)=r integral_0^1 dtheta
                                         partial_t psi(t+theta r)    (19)

represents the remainder by such derivative endpoints. Higher fixed-order
remainders follow by iterating the fundamental theorem of calculus.
Intermediate times stay in the convex hull of the original times. The
kernel is a pushforward of its original measure times at most d powers of
its diameter for a remainder of order d. Consequently its p-th-moment
bound follows from the original (p+d)-th-moment bound, since

    (s diam)^d(1+s diam)^p <= (1+s diam)^(p+d).          (20)

This inequality already includes the s^d in (18). It also handles products
of subtractions involving different fields. For lattice spatial
subtractions, use a telescoping sum of nearest-neighbor differences along
an interval between the original sites; the total variation costs at most
the corresponding spatial separation. Replace diam_t by a spacetime
diameter in the input moments when both sorts of subtraction occur.

Thus one localization step preserves the thermal power p if d additional
endpoint moments are available. This does not spend a new field-weight
factor for every Taylor path: the derivative order bound d0 and the
coefficient pushforward are included in the norm once. It does not prove
that an arbitrarily deep GN tree has a common (p+d)-moment budget; the
standard redistribution of Taylor distances over descendant lines must
still be verified in the complete multiscale construction.

## 9. The same bound for the whole one-step effective action

**[proved extension]** Fix an output moment order q>=0 and field weight
h>=1. Require the input and covariance bounds above through moment p+q,
and choose the input field weight H>C_(p,q)h. Define

    E_(s,beta)(V;psi)=-log E_(g_(s,beta)) exp[-V(psi+chi)].             (21)

Here chi is the integrated field and psi remains external. Keep the
scalar density and all external even kernels in E. Equip the output with
the norm

    ||E||_(s,h,q)=s^(-2)|scalar density|
       +sum_(E even>=2) h^E s^(E/2-2) ||kernel_E||_q,                (22)

using circular relative-time moments for the thermal kernels. Then

    ||E_(s,beta)(V)-Per_beta E_(s,infinity)(V)||_(s,h,q)
                  <= C_(p,q) u_s(V)(beta s)^(-p).                  (23)

The marked-defect version is linear in u_s(R), as in (6), and the
undifferenced effective action has norm at most C_(p,q)u_s(V). Periodization
in (23) includes each external field's antiperiodic sign. The scalar
component is interpreted per spacetime volume, as before.

Proof. At a fixed term with P input fields and E retained external fields,
the residual determinant has size

    r=(P-E)/2-m+1 >= 0.

Choosing which fields remain external costs at most 2^P. This is absorbed
by the field-weight reserve. If E>0, pin its first external coordinate
and root the BBF tree at the endpoint containing that field. The root
vertex anchor is then determined by that endpoint's relative coordinates.
The same unit-Jacobian tree unfolding integrates all other anchors.
When a lifted vertex carries external fields, their antiperiodic
periodization signs supply the additional factors absent from the
contraction determinant. Together with the internal fields they again
give (-1)^(n_i l_i)=1 at every even input endpoint.
The normalized scale factor is now

    s^[-(m-1)+r+sum_i(2-l_i/2)] = s^(2-E/2),            (24)

which cancels (22). If r=0, no residual contraction remains after the
tree unfolding. That term periodizes exactly and contributes zero to
the difference (23). If r>0, the determinant comparison is unchanged.

Each relative output field coordinate is connected to the pinned field
by a tree path and endpoint displacements. Its temporal diameter is
therefore bounded by R_t from (12). The output weight costs at most
(1+s R_t)^q. In the long-configuration estimate, multiplying this by
(s R_t)^p requires at most a (p+q)-th moment of one tree line or endpoint,
by the same finite-sum inequality (15). In the short region only the
q-th moment is needed in addition to the pointwise image bound. Output
circular diameter is at most the diameter of any chosen lift, so these
are valid estimates of the thermal kernel norm. Triangle inequality
allows summing the lifts after taking coefficient absolute values.
The weighted Cayley and analytic-radius bounds then prove (23). Cauchy's
estimate proves its marked version exactly as in Section 7.

Derivative endpoints from Section 8 are included by giving each external
field its derivative label and adding its total derivative order to the
output exponent E/2-2. The same power cancellation holds because each
internal field is contracted once and each external field is retained
once.

This provides a finite-temperature forcing bound for an entire single
RG map. It does not justify spending H/h>1 anew at infinitely many scales.
A full iteration needs a common analytic endpoint hierarchy or another
proof that this reserve is available uniformly, as already emphasized
in the repository's degree-reparameterization work.

## 10. A complete shell partition has no extra terminal measure

**[proved conditional on the stated decomposition]** The separate final
term q(T) in THERMAL_CRITERION.md can be zero by convention. Suppose the
entire infrared covariance is decomposed into smooth shells h<=0 with
scale s_h=gamma^h and support

    |omega|<=B s_h,                                    (25)

for a fixed B>0. The shells form an exact partition on the infrared
Fourier support. Fermion frequencies are omega=(2m+1)pi T, so every shell
with B s_h<pi T has identically zero covariance. Put

    H_T=ceil(log_gamma(pi T/B)).                       (26)

Only h>=H_T can contribute. Integrate every one of those scale measures,
retaining the constant term produced at each step. After the last such
integration, the remaining covariance is exactly zero. Gaussian
expectation with zero covariance evaluates the remaining Grassmann
polynomial at zero. If its constant was recorded in the last scalar
increment, there is no remaining integral and q(T)=0 exactly.

The endpoint case B s_(H_T)=pi T may make that shell zero as well if its
smooth cutoff vanishes at its support boundary. Including this zero
integration is harmless. The ratio s_(H_T)/T is bounded above and below
by fixed positive constants. The estimates in this note extend from
beta s>=1 to beta s>=c>0, with constants depending on c: alternatively
rescale the shell variable so B<=pi. Thus this terminal convention is
compatible with the thermal winding estimate and with the geometric
summation in THERMAL_CRITERION.md.

At zero temperature there are additional shells h<H_T. If their
reference-subtracted scalar bounds are (7) with defect omega(s_h), their
sum is bounded by

    C s_(H_T)^2 sup_(0<s<=s_(H_T)) omega(s)=o(T^2).      (27)

No separate terminal comparison is required under this exact complete
partition convention. This does **not** remove Gaussian partition
prefactors or field-renormalization Jacobians: they must be included
among the recorded scalar increments or kept as explicit known factors.
Nor does it compare the thermal and line running kernels; that remains
an independent obligation. The observation was suggested by the parent
coordinator and checked here from the Matsubara support condition.

## 11. What this proves for a GN tree, and what it does not

The theorem applies directly to a single GN integration vertex, including
a vacuum root, once its child kernels obey (3) with one common field-weight
reserve and its normalized scale covariance obeys (1)-(2). It includes
all scalar contractions produced at that vertex. It also applies to a
complete single-scale analytic interaction without choosing a particular
finite number of endpoint fields.

Four additional assertions are needed before (17) becomes the full
microscopic scalar estimate.

1. **Reference-subtracted kernels.** At every scale one must construct
   V_h and V_*h in a common representation and prove
   u_(s_h)(V_h-V_*h)<=A omega(s_h), including the relevant counterterm,
   the forward interactions, and all irrelevant kernels. Bounded marginal
   couplings alone do not say that this difference tends to zero.
2. **Thermal history.** The actual finite-temperature child kernels are
   not exactly the periodizations of their line values. Their difference
   must be propagated through localization and the RG history. The causal
   lemma in THERMAL_CRITERION.md gives a sufficient algebraic closure
   after the forcing and history Lipschitz bounds are supplied. No
   temperature-dependent bare tuning is allowed.
3. **Localization and normalization.** Taylor subtraction can place time
   moments on descendants and introduce derivatives or interpolated
   coordinates. A complete GN proof must show that these operations
   preserve one common moment and field-weight reserve. Section 8 proves
   the one-step statement, including derivative endpoints, with d extra
   moments for a subtraction of order d. The complete proof must retain the
   Gaussian normalization constants and anomalous field-rescaling
   Jacobians. If different reference covariances are used, their exact
   determinant prefactors cannot be dropped when rewriting them with
   one common covariance.
4. **Reference evaluation.** The dressed interacting reference's thermal
   coefficient must be evaluated with its actual left and right
   velocities. The free energy of a bare free-fermion reference is not
   generally the correct coefficient. Section 10 removes an additional
   terminal remainder if the exact shell partition is exhausted; an
   alternative convention stopping earlier must bound its remaining
   integral relative to the reference by T^2 omega(T).

The massive matching theorem's exponentially accurate periodization is a
useful initial input for item 2. The multichannel beta-function work
addresses parts of item 1. This note proves the determinant winding step
needed at fixed running data; it does not silently promote those other
obligations to established hypotheses of the actual wrapped chain.
