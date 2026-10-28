# Fermi-shift existence without differentiated beta functions

Status: **[proved conditional existence theorem]**. The conditions are the
ordinary uniform beta-function bounds in
[COUNTERTERM_INDUCTION.md](COUNTERTERM_INDUCTION.md), continuity at each
fixed infrared cutoff, and the exact ultraviolet matching relation below.
No derivative with respect to a Fermi point is needed. This does not prove
those beta bounds for the entire retarded-endpoint class or identify a CFT.

Author: OpenAI GPT-6-family Codex vacancy worker, 2026-10-28. AI-assisted;
no human review. The coordinating agent proposed the Brouwer/compactness
replacement; this note gives its complete finite-cutoff proof.

## 1. The obstruction that existence avoids

A branch-changing microscopic endpoint, expressed in fields centered at
Fermi points, carries a phase depending on its absolute anchor. Changing
a Fermi point differentiates that phase and produces the anchor position.
Exponential moments in relative coordinates do not control that factor
in the supremum over anchors. Consequently an initial relative-coordinate
norm bound alone does not imply a uniformly small derivative of the
infinite-scale counterterm map.

Existence of suitable Fermi points is enough here. A derivative estimate
would supply a stronger uniqueness or smooth-dependence statement, but
it is unnecessary for existence. The finite-scale argument below never
differentiates that phase and never requires continuity in the supremum
over anchors of a charged endpoint.

## 2. Precise hypotheses, including physical matching

Let there be a fixed finite number m of simple Fermi branches. Write p in
R^m for candidate shifts of their reference centers. Choose a fixed
compact neighborhood P of zero on which all velocities are nonzero,
branches stay separated, and the quadratic/quartic nonresonance margins
used for localization persist. Constants may depend on this fixed cell
and on these margins. There is no claim of uniformity as m grows.

After division by the logical bandwidth, the free dispersion evaluated
at the displaced centers has the form

    d(p)=V p+Q(p),
    V=diag(v_1,...,v_m),  v_*:=min_a |v_a|>0,
    ||Q(p)||_infinity <= C_Q ||p||_infinity^2.              (1)

The same proof allows any invertible V, replacing 1/v_* by its induced
norm. Assume the following on P, or on a fixed sufficiently small
subneighborhood independent of epsilon.

**H1: ordinary RG bounds.** The causal reduced beta functions obey
(G1) and (G2) of COUNTERTERM_INDUCTION, with uniform constants, for every
p and every infrared cutoff. Here (G2) compares running histories at
the same p and the same initial marginal vector g_0(p). It does not
compare different Fermi shifts. The initial forward couplings, initial
irrelevant data, and velocity displacements have the smallness required
there; the reference velocities remain in a fixed compact admissible set.

**H2: fixed-cutoff continuity.** For each fixed scale h, its causal beta
function is continuous in p and its finite history from h to zero.
The initial marginal vector and fixed initial irrelevant kernels enter
through a continuous finite-cutoff construction. All single-step
normalizations have already been solved continuously. No continuity
uniform in h is required.

**H3: exact initial matching.** Before the artificial local relevant
counterterm is varied, the physical ultraviolet action has initial
relevant coefficient

    nu_0^phys(p)=A(p)d(p)+u(p),                            (2)

where A is continuous and invertible, A and A^-1 are uniformly bounded,
and u is continuous with ||u(p)||<=C_u epsilon. The beta functions in H1
are formed by keeping this physical initial action, including its
initial marginal and irrelevant kernels, fixed at p, and varying only
the localized relevant coefficient nu_0. Thus equality of that one
coefficient recovers the exact physical action.

The simplest convention is A=I and u=0. Equation (2) also permits a
small finite ultraviolet self-energy and finite field normalization.
It must be checked in an application: an O(epsilon) counterterm theorem
does not license changing the physical Hamiltonian. In particular,
microscopic chemical potentials need not be independently adjustable.
All dependence of the physical ultraviolet data on d(p) is retained in
the fixed p-dependent data, not silently discarded when nu_0 is varied.

## 3. The continuous finite-cutoff counterterm

Fix N>=0 and retain scales h=-N,...,0. Put nu_(-N-1)=0 as the terminal
boundary condition. The finite version of the relevant fixed-point map is

    (T_(N,p) nu)_h
      =-sum_(j=-N)^h 2^(j-h-1) r_j(p,nu),  -N<=h<=0.     (3)

The causal marginal trajectory is obtained from the same finite history
by

    g_h=g_0(p)+sum_(j=h+1)^0 b_j.                          (4)

Every estimate in Sections 4--6 of COUNTERTERM_INDUCTION remains valid
when these sums are finite. In particular, for epsilon sufficiently
small, uniformly in p and N, (3) is a contraction of the finite weighted
ball with a contraction factor q<1 independent of N,p. It has a unique
fixed point and associated marginal increments satisfying

    |nu_h^(N)(p)| <= A_0 epsilon 2^(theta h),
    |b_h^(N)(p)|  <= B_0 epsilon 2^(theta h).              (5)

H2 implies continuity of this fixed point for each fixed N. One direct
proof starts Picard iteration at zero: every iterate is continuous and
the convergence is uniform on P at geometric rate q. The causal marginal
recursion preserves continuity. This argument uses no p derivative.

Set

    c_N(p)=A(p)^-1 [nu_0^(N)(p)-u(p)].                    (6)

This map is continuous and has the uniform bound

    ||c_N(p)||_infinity <= C_c epsilon.                   (7)

All finite ultraviolet self-energy contributions are included in u and
therefore in (6). No unproved cancellation of those contributions is
used. Constants can be enlarged so C_c>0.

## 4. Brouwer chooses physical centers at each cutoff

Take R=2 C_c/v_*. For sufficiently small epsilon the cube

    K_epsilon={p: ||p||_infinity<=R epsilon}              (8)

lies in P. Define

    F_N(p)=V^-1[c_N(p)-Q(p)].                             (9)

It is a continuous map of the cube into itself because

    ||F_N(p)||_infinity
      <=(C_c epsilon+C_Q R^2 epsilon^2)/v_*
      <=R epsilon,                                       (10)

where the last inequality holds when C_Q R^2 epsilon<=C_c.
Brouwer's fixed-point theorem gives a p_N in K_epsilon with

    d(p_N)=c_N(p_N),
    nu_0^phys(p_N)=nu_0^(N)(p_N).                         (11)

Thus the counterterm needed at this cutoff is exactly provided by
rewriting the given physical action around the shifted centers. It is
not an additional physical coupling. This argument proves neither
uniqueness of p_N nor its differentiability in microscopic parameters.

## 5. Remove the infrared cutoff by diagonal compactness

The centers p_N lie in the compact finite-dimensional cube (8), so pass
to a subsequence with p_N->p_*. For each fixed h, (5) bounds
nu_h^(N)(p_N), b_h^(N)(p_N) and g_h^(N)(p_N) in compact sets whenever
N>=-h. A diagonal subsequence gives simultaneous limits at all integer
scales. If normalized beta functions retain additional bounded
single-step normalization ratios as explicit coordinates, include these
finitely many coordinates per scale in the same diagonal extraction.
No convergence of the anomalously running Z_h as h->-infinity is used.

For any fixed h, the beta equations depend only on the finite history
from h to zero. By H2 their limit is the same causal beta equation at
p_*. Consequently

    b_h=B_h^g(p_*,{g_k,nu_k}_(k=h)^0),
    r_h=B_h^nu(p_*,{g_k,nu_k}_(k=h)^0),
    g_(h-1)=g_h+b_h,
    nu_(h-1)=2 nu_h+r_h.                                 (12)

The uniform bounds (5) pass to the limit. They imply nu_h->0 and
absolute summability of b_h. Backward iteration of the relevant
recursion gives, for M sufficiently large that -M<h,

    nu_h=2^(-M-h-1)nu_(-M-1)
           -sum_(j=-M)^h 2^(j-h-1)r_j.                   (13)

The terminal term tends to zero by (5), and the r_j bound from the
abstract induction makes the series absolutely convergent. Hence the
infinite-scale tail equation holds exactly. Also

    g_*=g_0(p_*)+sum_(j<=0)b_j,
    |g_h-g_*|<=C epsilon 2^(theta h).                    (14)

The limiting velocities remain nonzero. Finally, pass to the limit in
the second equality of (11), using continuity of (2) and convergence
of nu_0^(N)(p_N), to obtain

    nu_0^phys(p_*)=nu_0.                                 (15)

This is the exact initial matching condition for the infinite-scale
trajectory of the given physical model. The shifts obey |p_*|<=R epsilon
and stay within all fixed certified momentum margins.

The proof does not assert convergence of the finite-N sequences in the
weighted norm sup_h 2^(-theta h)|nu_h|: the moving terminal scale need
not converge in that norm. Coordinatewise compactness and the uniform
decay suffice. It also does not require an a priori continuous
infinite-IR map c(p). If such a map is proved separately, Brouwer can
instead be applied directly to it.

## 6. How fixed-cutoff continuity avoids the anchor problem

Here is a sufficient concrete criterion for H2. Use smooth covariance
cutoffs and centered dispersions

    e_a(k_F,a+p_a+q)-e_a(k_F,a+p_a),                      (16)

keeping the subtracted constant in (2). At a fixed infrared cutoff,
the finitely many scale covariances have uniform integrable moment
bounds on P and are continuous in the corresponding fixed-scale
weighted L1 norms. Assume each fixed-degree initial pinned kernel is
continuous in its weighted total-variation norm as p varies, and the
exponential-moment bounds and all-degree hierarchy hold uniformly on P.
An alternative sufficient condition is a common positive dominating
measure with an integrable exponential moment, with p-dependent
densities continuous almost everywhere and bounded by a common
integrable majorant. This includes singular time-local kernels.
Pointwise continuity together with bounded L1 norms alone would not
be sufficient: it would allow concentrating spikes. Let the
renormalized finite-scale tree expansion be absolutely convergent with
uniform majorants on P, as is required in the construction behind H1.

Each localized quadratic or forward quartic beta root is branch neutral.
Pin that root at the origin. Propagators conserve the branch label, so
the sum of all endpoint branch charges in a contributing tree vanishes.
The potentially troublesome phases can therefore be expressed entirely
in terms of positions relative to the pinned root. For fixed positions
they are continuous in p and have modulus one. When the p dependence
is such a phase multiplying a fixed exponentially integrable kernel,
that kernel's total variation is the required common domination. No
anchor-supremum norm continuity is needed for this claim. The stated
weighted continuity or domination criterion justifies continuity of
each tree integral, and the uniform convergent tree/degree
majorants permit passage through the sums. Localization and fixed-scale
normalization remain continuous because the retained local channels
are fixed by the positive margins. This proves H2 under the stated
criterion.

This argument allows charged high-degree endpoints and their translated
phases. It does not differentiate them. It requires neither an infinite
Diophantine condition on higher-degree resonances nor a uniform
anchor-norm continuity statement. The branch neutrality used here is
a property of the beta root and total charge balance of a diagram;
it is not an assumption of separate branch conservation at every
microscopic endpoint.

## 7. What is, and is not, now settled

Together with the abstract recursion theorem, H1--H3 prove existence of
a physical infinite-scale trajectory with O(epsilon) Fermi shifts,
summable marginal increments and nonzero velocities. Only history
differences at fixed p enter the contraction estimates. The previously
listed differentiated-Fermi-point bound is unnecessary for this result.

The outstanding analytic task is still to derive the ordinary reduced,
normalized (G1)--(G2) bounds for all allowed retarded endpoints, including
the extra small factor in the history-difference terms. Fixed-cutoff
continuity does not supply that factor. Source-dependent current
correlations, pressure/stress control and a conformal scaling limit are
further results; none follows solely from this existence theorem.
