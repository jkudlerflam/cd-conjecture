# Fixed moments from the matched tree norm, and localization accounting

Status: **[proved initial norm transfer and allocation lemmas; ordinary
replacement-endpoint GN routing remains open]**. The temporal radius need
not shrink with field degree. The remaining issue is repeated *internal*
Taylor-distance allocation, not a failure of the matched initial norm.

Provenance: OpenAI GPT-6-family Codex thermal worker, 2026-10-28.
AI-assisted; no human review.

## 1. The actual input being audited

The direct-attack [tree-norm extension](../../direct_attack_20261028/infrared/TREE_NORM_EXTENSION.md)
controls, for a fixed field weight H,

    sum_(l even>=2) H^l integral |K_l(X)| exp[mathcal T_(a,kappa)(X)]
                                                        <= r,       (1)

with a pinned external coordinate. Here mathcal T is the minimum
spanning-tree length for metric a|dx|+kappa|dt|. Its construction charges
this length to the occupied polymer network; it never uses l times the
external diameter. This note uses that theorem as input and checks its
consequences for the thermal RG norms.

Choose **a>0**, not merely a>=0. With a=0, (1) alone cannot control a
spatial derivative of a Fourier kernel: a kernel supported on two distant
sites can have fixed temporal norm and unbounded first spatial moment.
The polymer theorem allows a positive a to be chosen before its input
support exponent and the final microscopic coupling.

## 2. Rescaled time and a fixed reserve

Set u=W t with fixed W>0. The transformed pinned kernels have the common
factor 1/W in their integrated coefficient norm, independently of l.
Their tree metric is

    d_W=a|dx|+(kappa/W)|du|,
    mathcal T_W=mathcal T_(a,kappa)(x,u/W).              (2)

In particular their bound is r/W in (1) with mathcal T_W. Every pair of
external points is connected by the minimum tree, so

    a diam_x <= mathcal T_W,
    (kappa/W) diam_u <= mathcal T_W.                    (3)

No field-count factor occurs in (3). Fix 0<delta<1 and reserve the
remaining exponential weight exp[(1-delta)mathcal T_W]. For integers
q_x,q_u>=0, q=q_x+q_u, and 0<s<=1,

    (s diam_x)^q_x (s diam_u)^q_u
                    exp[(1-delta)mathcal T_W]
       <= q! delta^(-q) s^q a^(-q_x)(W/kappa)^q_u
                                                 exp[mathcal T_W].  (4)

Indeed (3) bounds the polynomial by its displayed constant times
mathcal T_W^q, and y^q exp(-delta y)<=q!delta^(-q).
Consequently five spacetime moments and a nonzero tree-weight reserve
are available with one fixed choice, for example delta=1/2. The 5!
is a fixed constant, not a factorial in field degree or perturbative
order. A bound for (1+s diam_x+s diam_u)^5 follows by its finite binomial
expansion. The particular fourth-temporal/first-spatial moment needed
for the velocity coefficient in THERMAL_LOCALIZATION.md costs at most

    5! 2^5 s^5 a^(-1)(W/kappa)^4.                     (5)

The moment bound becomes better as the normalized memory W/kappa
shrinks. It does not demand a degree-dependent temporal radius.

At the initial normalized scale s=1, (4) proves exactly the fixed-moment
input of WINDING_ESTIMATE.md. Its factors s^(l/2-2) are harmless for
l>=4 and 0<s<=1. At l=2 a factor 1/s remains; that is the usual relevant
quadratic dimension and must be handled by localization/counterterms.
It is not justified to claim the entire unlocalized initial quadratic
kernel is uniformly small after arbitrary further scale reductions.

## 3. One field reserve covers degree polynomials

Suppose changing from orbital fields to the fixed finite branch basis
costs at most C_b^l, and a fixed number of marks or momentum derivatives
costs l^m. Choose h with q0=C_b h/H<1. Then

    sum_l h^l C_b^l l^m ||K_l||
          <= [sup_(l>=2) l^m q0^l] sum_l H^l ||K_l||.  (6)

The supremum is finite and independent of l. For example, writing
c=-log q0>0, l^m exp(-cl)<=m!c^(-m), by the exponential series.
For a *fixed* m this is one finite loss. It can be combined with (4)
because the field and tree weights are separate factors in (1).

The output matching theorem permits any fixed H, at the cost of a larger
input support exponent. Thus H can be selected at the beginning to pay
the branch conversion, the five moments, and a fixed number of history
or thermal marks. This observation does not authorize spending H/h at
every one of infinitely many scales. The degree hierarchy in the direct
attack is intended to supply the global GN sum with one reserve.

The branch conversion cost can be checked directly. On each fixed
nondegenerate band patch, choose a smooth eigenvector frame and a smooth
cutoff. Its inverse spatial Fourier transform B(x) has finite weighted
sum

    C_5=sum_x (1+2|x|)^5 ||B(x)|| < infinity,           (6a)

by integration by parts more than six times. There are finitely many
patches. Each external field is convolved with one such B. If the
old coordinates are Y_i and the new ones are X_i=Y_i+delta_i, then

    1+diam(X) <= (1+diam(Y)) product_i(1+2|delta_i|).

Thus their fifth-moment kernel norm grows by at most C_5^l, including
the finite branch-label sum. Equation (6) pays this with one field-weight
reserve. The Fermi oscillations have modulus one. A smooth cutoff need
not preserve an exponential tree norm after this convolution; only the
fixed polynomial moment norm is claimed after branch conversion. This
is sufficient for the initial winding and localization estimates.

## 4. Two exact simplifications of repeated localization

Let L_2 be Taylor projection through total degree one on a quadratic
slow-momentum kernel, and let L_4 be projection to a quartic constant.
Write R_i=I-L_i.

**[proved]** Consecutive identical localization operations on a line with
no intervening contraction reduce to one, because R_i^2=R_i. This also
holds for the fixed-temperature projections of THERMAL_LOCALIZATION.md.
Thus the scale length of a unary string does not by itself create new
Taylor derivatives or moments. Scalar field-normalization factors commute
with these projections. A step that changes external variables or performs
new contractions is not such a unary string and is not erased by this
observation.

There is a useful derivative routing statement for zero-temperature
Taylor localization. Represent derivatives on explicit external field
half-lines. Suppose a term already contains a homogeneous slow-momentum
polynomial P_d of total degree d from those derivatives. For a projector
removing Taylor degrees below z, where z=2 for a quadratic cluster and
z=1 for a quartic cluster,

    R_z(P_d F) = P_d R_(z-d) F       when d<z,
    R_z(P_d F) = P_d F              when d>=z.          (7)

Here R_0=I. This follows by multiplying the Taylor series of F by P_d;
only degrees below z-d can survive in the local part. It remains true
termwise for a sum of homogeneous derivative polynomials.

**[proved allocation consequence]** In a termwise derivative-field
implementation, apply only the deficit z-d in (7). A new localization
then never raises an individual external half-line derivative order
above two. Surviving derivatives are carried forward. When a differentiated
half-line is contracted, its derivative is assigned to that propagator
and is no longer an external derivative of a larger cluster. The two ends
of a contracted line therefore carry at most four derivatives in total.
This proves a bounded derivative order for this representation; it does
not prove the full coefficient/moment bound for the nested expansion.

The exact factorization (7) is a zero-temperature Taylor statement.
For the finite-temperature projection, separate it as
R_T=R_infinity+(L_infinity-L_T). The first piece uses (7); the second is
the thermal mark bounded in THERMAL_LOCALIZATION.md. One must not assume
that arbitrary non-affine derivative polynomials commute exactly with
the four-point interpolation map.

## 5. The distinct distance-allocation issue

Taylor localization creates both differentiated external fields and
coordinate differences inside the kernel. The previous derivative count
does not automatically bound how many coordinate differences a naive
tree-path expansion can assign to one **internal** propagator or one
initial endpoint.

There is a concrete obstruction to an arbitrary allocation. For t>=0,
take K(t)=exp(-2t), whose norm with weight exp(t) is finite. Then

    integral_0^infinity t^(2r) K(t)dt = (2r)!/2^(2r+1). (8)

If r second-order Taylor distances are all charged to that same t,
even their factors 1/(2!)^r do not make the result bounded by C^r for
any fixed C. A field weight of the form H^l cannot cure this at fixed l.
This is a counterexample to unrestricted distance allocation, **not**
to a correctly reorganized GN expansion or to the conjecture. Applying
R_2 repeatedly to this one kernel would actually be idempotent, so such
a repetition cannot be used as a physical counterexample either.

An ordinary replacement-endpoint proof must give its zero/path routing
rule and show how it avoids (8), for example by bounded load on each
line and endpoint-tree segment, or by an explicit summation that pays
the high-load alternatives. Citing fixed-order moments alone does not
supply this statement. The tree norm's exponential reserve does supply
every fixed moment; it does not supply factorially many moments with
only an exponential constant in their order.

## 6. A fixed thermal mark adds no new accumulation problem

Here is an allocation lemma that separates the extra thermal work from
that ordinary routing problem. Consider a completed connected BBF/GN
term whose *ordinary* Taylor-distance allocation has already been bounded.
Flatten its nested spanning trees into a tree joining its N initial/local
endpoints. Add a minimum terminal tree inside each nonlocal endpoint.
Let E_net be the number of network segments. It is at most P+N, where
P is the total initial/local field degree. Every covariance segment is
at a scale s_k>=s_h, where h is the root scale.

Any external temporal diameter is bounded by the sum of the temporal
lengths of these segments. For a fixed integer p,

    (s_h sum_e length_t(e))^p
       <= E_net^(p-1) sum_e (s_h length_t(e))^p.        (9)

On a covariance segment of scale k,

    (s_h |t_e|)^p <= (s_k |t_e|)^p.                   (10)

Hence a fixed root moment can be assigned to just **one** line or one
initial endpoint per summand. If the ordinary routing uses line moments
through order m0, this requires only order m0+p there. For an initial
endpoint, (4) supplies the additional moment with one fixed reserve.
The number of choices contributes at most (P+N)^p. Equation (6) and a
single reduction of an analytic endpoint-count radius pay that polynomial:
sum_N N^p z^N converges for |z|<1. There is no loss per scale.

Similarly, telescoping a product of residual determinants creates only
a polynomial field count. Every original field is either external or
contracted once, so the sum of their sizes is at most P/2. Thus

    sum_v (size of determinant v)^2 <= (P/2)^2.        (11)

The mixed-cofactor estimate of WINDING_ESTIMATE.md can be applied at the
one marked determinant. Its image factor is

    (T/s_k)^p=(T/s_h)^p gamma^[-p(k-h)],               (12)

which improves the root-to-mark scale sum. It consumes no dimensional
short-memory reserve. A localization thermal mark has the same property
with p=3. A separate reference-defect or history mark is handled by the
ordinary marked-tree reserve.

Equations (9)-(12) prove that **once ordinary zero routing is controlled**,
the fixed thermal comparison does not introduce a new unbounded derivative
or moment order. They do not replace the missing ordinary routing proof.
In particular, it would be circular to use a repeated single-step H-to-h
loss as that proof and then assert a common field radius at all scales.

## 7. Audit conclusion

The strengthened massive-matching norm passes the initial thermal norm
test: choose a>0, reserve a fixed part of its tree exponential, and make
one field-weight reduction. This supplies all five moments, branch-label
conversions and fixed marked-degree factors required here, uniformly in
field degree. The normalization/covariance lemmas in the direct attack
also retain any fixed derivative/moment order.

The currently read direct-attack files still label the ordinary
replacement-endpoint zero routing and global field sum as unfinished.
They do not yet provide a theorem to which this note could attach a final
microscopic thermal conclusion. The precise remaining accounting issue
is the distance load in Section 5. The elementary derivative routing and
the fixed-mark allocation above narrow it; they do not conceal it behind
the phrase that localization is routine.
