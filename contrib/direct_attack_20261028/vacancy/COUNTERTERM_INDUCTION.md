# A closed marginal-flow and counterterm induction

Status: **[proved abstract recursion theorem; application to every
retarded endpoint is not yet proved]**.

This note closes the nonlinear sequence argument once explicitly stated
GN beta-function estimates are available. It separates that theorem
from the still-needed derivation of those estimates for the entire
replacement-kernel class. No CFT or central charge is claimed.

Author: OpenAI GPT-6-family Codex vacancy worker, 2026-10-28. AI-assisted;
no human review. The elementary recursion proof below is new to this
folder. The construction is independent of the collaborator's C47
single-mode campaign.

## 1. Why the anomalous wave-function flow must remain normalized

The multichannel calculation cannot copy a single-chirality proof by
setting all wave-function ratios to one. Those ratios generally have
an anomalous, nonsummable logarithmic flow. Work instead with the reduced
beta functions in normalized fields: single-scale field normalization
is solved inside the RG step, and its small ratio is retained in the
propagator and determinant bounds.

The sequence variables below comprise only forward quartic couplings
and velocities. The reference-theory comparison concerns these reduced
beta functions. A separate estimate `|log(Z_(h-1)/Z_h)|<=C eta` suffices
for the dimensional bounds; no convergence of `Z_h` is assumed.

## 2. The explicit beta-function hypotheses

Let scales be integers `h<=0`. The finite-dimensional vector `g_h`
contains the marginal couplings and velocities; `nu_h` is the vector
of normalized relevant quadratic terms. All vector norms are fixed
finite-dimensional norms. Fix `0<theta<1/4` and constants `C>0`.

The initial forward couplings have magnitude at most `epsilon` and the
initial velocities lie in a fixed compact region bounded away from zero.
Let `eta=sqrt(epsilon)`. Throughout a provisional domain, forward
couplings and velocity displacements remain less than `eta`, while
the velocities themselves stay in that fixed compact region. The
normalized relevant coefficients also have magnitude less than `eta`.

Assume finite-cutoff RG defines causal reduced beta functions

    b_h = B_h^g({g_k,nu_k}_{k=h}^0),
    r_h = B_h^nu({g_k,nu_k}_{k=h}^0),
    g_(h-1)=g_h+b_h,
    nu_(h-1)=2 nu_h+r_h.                                  (1)

Any implicit dependence from single-step normalization is already solved
in these functions. Their definitions and the constants are uniform in
the infrared cutoff and in volume after the usual thermodynamic limit.

**(G1), absolute bound.** For either component `r=g,nu`,

    |B_h^r| <= C 2^(2 theta h) [epsilon+|g_h-g_0|]
       + C eta sum_(k=h)^0 2^[-4 theta(k-h)]
                              [|g_k-g_h|+|nu_k|].          (2)

Define the nonconstant-history size

    Q_h(g,nu)=sum_(k=h)^0 2^[-4 theta(k-h)]
                               [|g_k-g_h|+|nu_k|].

**(G2), difference bound with fixed `g_0`.** For two histories in the
same domain, put `delta g=g-tilde g`, `delta nu=nu-tilde nu`. Then

    |delta B_h^r| <= C 2^(2 theta h) |delta g_h|
       + C eta sum_(k=h)^0 2^[-4 theta(k-h)]
                [|delta g_k-delta g_h|+|delta nu_k|]
       + C |delta g_h|[Q_h(g,nu)+Q_h(tilde g,tilde nu)].  (3)

These are sufficient bounds, not asserted consequences of smallness
alone. The first term permits a first-order but scale-decaying velocity
correction. The history terms carry the extra small factor `eta`.
The exact linear relevant evolution `2 nu_h` is excluded from the
beta function. Counting it again in (3) would invalidate the contraction.
The last term in (3) is needed when differentiating the coefficient of
a history deviation. For example, `g_h(g_(h+1)-g_h)` produces precisely
such a term. Omitting it would demand a stronger property not implied
by a generic analytic marked-tree estimate. The alternatives worker
identified this issue in the first version of the note.

For finite-volume applications one can add a cutoff error tending to
zero before the infrared limit, provided it obeys the same summable
weighted bounds. It is omitted here to keep the theorem precise.

## 3. Abstract theorem

**[proved]** Under (G1)–(G2), for sufficiently small `epsilon`, there is
a unique relevant sequence in a ball

    ||nu||_theta:=sup_(h<=0) 2^(-theta h)|nu_h|
                     <=A epsilon,                         (4)

with its associated causal marginal trajectory, satisfying (1) and the
infrared boundary condition specified by the convergent series (16)
below. Constants depend on `C,theta` and the finite-dimensional domain,
not on the number of RG scales. The marginal increments obey

    ||b||_theta<=B epsilon,
    |g_h-g_*|<=B' epsilon 2^(theta h).                     (5)

The bare relevant counterterm is therefore `nu_0=O(epsilon)`, all
forward couplings stay `O(epsilon)`, and all velocity shifts are
`O(epsilon)`. In particular they remain in the provisional domain
`O(epsilon)<eta` when `epsilon` is small enough.

The constants can be inefficient. This theorem supplies a closed
induction and fixed point, rather than assuming bounded running
couplings as its conclusion.

## 4. The triangular kernel, whose norm need not be small

For a fixed relevant sequence, the marginal flow is computed causally:

    g_h=g_0+sum_(j=h+1)^0 b_j.                             (6)

Use (6) inside (2). Interchanging the finite positive sums in the
history term yields constants `C0,C1,C2` such that

    |b_h| <= C0 epsilon 2^(2 theta h)
            +(K0 |b|)_h+(K1 |b|)_h
            +C2 eta ||nu||_theta 2^(theta h),              (7)

where

    (K0 f)_h=C0 2^(2 theta h) sum_(j=h+1)^0 f_j,
    (K1 f)_h=C1 eta sum_(j=h+1)^0
                            2^[-4 theta(j-h)] f_j.        (8)

For example,

    sum_(k=h)^0 2^[-4 theta(k-h)] |g_k-g_h|
       <=[1/(1-2^(-4 theta))]
          sum_(j=h+1)^0 2^[-4 theta(j-h)] |b_j|.

On the weighted sequence norm, `||K1||<=C3 eta`, since
`sum_(n>=1)2^(-3 theta n)` converges. The operator `K0` need not have
norm smaller than one, so a direct supremum-norm contraction would
be unjustified.

Nevertheless `K0` is strictly triangular. Let

    S_theta=sum_(j<=0)2^(theta j)=1/(1-2^(-theta)).

Repeated ordered summation gives

    ||K0^n||_theta <= (C0 S_theta)^n/n!.                    (9)

To see this directly, write `f_j=2^(theta j) tilde f_j`.
The normalized kernel is `C0 2^(theta h)2^(theta j)` for `j>h`.
An `n`-fold path has strictly increasing indices
`h<j_1<...<j_n<=0`. Its intermediate factors
`2^(2 theta j_i)` are bounded by `2^(theta j_i)`; the initial
factor is at most one. The sum over strictly ordered products is
bounded by the full product sum divided by `n!`. This proves (9).

Thus the positive Volterra resolvent exists with

    ||(I-K0)^-1||_theta <= E_theta:=exp(C0 S_theta).         (10)

If `E_theta C3 eta<=1/2`, (7) implies

    ||b||_theta <=C4 [epsilon+eta ||nu||_theta].             (11)

These estimates hold first on every finite set of scales, where the
causal equations are ordinary finite recursions. The bounds are uniform;
each fixed scale depends on only finitely many higher scales. Therefore
they define the infinite causal trajectory without a cutoff-dependent
selection. A continuation/bootstrap argument keeps it inside the
provisional domain, because (11) is `O(epsilon)` on (4).

## 5. Lipschitz dependence of the marginal trajectory

For two relevant sequences with the same `g_0`, use (3) and

    delta g_h=sum_(j=h+1)^0 delta b_j.

The absolute bootstrap (11), applied to each history, first gives

    Q_h(g,nu)+Q_h(tilde g,tilde nu)
                                      <=C epsilon 2^(theta h)

when both relevant sequences are in (4). Indeed each marginal history
difference is a partial sum of the bounded increments, and convolution
against `2^[-4 theta(k-h)]` preserves the `2^(theta h)` decay.
The positive kernels therefore give

    |delta b_h| <= (K0 |delta b|)_h+(K1 |delta b|)_h
                    +(K2 |delta b|)_h
                    +C2 eta ||delta nu||_theta 2^(theta h).

Here

    (K2 f)_h=C epsilon 2^(theta h) sum_(j=h+1)^0 f_j,
    ||K2||_theta<=C epsilon S_theta.

After reducing epsilon so
`E_theta (C3 eta+C epsilon S_theta)<=1/2`, apply the same resolvent:

    ||delta b||_theta <=C5 eta ||delta nu||_theta.           (12)

Consequently

    sup_h |delta g_h|<=S_theta C5 eta ||delta nu||_theta.    (13)

This is the small factor needed for the relevant fixed point. It did
not require the full marginal recursion to be a contraction in the
unweighted running-coupling norm.

## 6. The relevant fixed point and its bare counterterm

Let `r_h(nu)` be the relevant beta function evaluated on the causal
marginal trajectory associated with `nu`. Equations (2), (6), and (11)
imply

    ||r(nu)||_theta <=C6 epsilon+C7 eta ||nu||_theta.        (14)

Likewise (3), (12), (13), and the just-proved `Q_h` bound imply

    ||r(nu)-r(tilde nu)||_theta
                      <=C8 eta ||nu-tilde nu||_theta.      (15)

The solution of the relevant recursion with a bounded decaying infrared
tail must satisfy

    (T nu)_h=-sum_(j=-infinity)^h 2^(j-h-1) r_j(nu).       (16)

This series converges absolutely in the weighted norm. The linear tail
operator has norm

    D_theta=1/[2(1-2^[-(1+theta)])].                       (17)

Choose `A=2D_theta C6`, then take `epsilon` small enough that
`D_theta C7 eta<=1/2` and `D_theta C8 eta<1`. Equations (14)–(17)
show that `T` maps (4) into itself and is a strict contraction.
Banach's fixed-point theorem proves existence and uniqueness in that
ball. Rearranging (16) verifies

    nu_(h-1)=2 nu_h+r_h.

It also identifies the bare counterterm explicitly:

    nu_0=-sum_(j<=0)2^(j-1)r_j,       |nu_0|<=A epsilon.   (18)

Finally, the absolute summability of `b_h` gives
`g_*=g_0+sum_(j<=0)b_j`; summing its tail proves (5).
This completes the abstract induction.

## 7. Counterterm versus changing the physical Hamiltonian

The theorem above fixes the quadratic reference and solves for a bare
quadratic counterterm. Application to a given physical Hamiltonian must
still absorb this counterterm into its Fermi-point locations; it must
not silently add independently tunable physical terms.

A sufficient existence argument is proved in
[FERMI_SHIFT_EXISTENCE.md](FERMI_SHIFT_EXISTENCE.md). Require (G1)--(G2)
uniformly over a fixed small neighborhood of reference Fermi centers,
continuity of every finite-cutoff beta function in those centers, and
the exact initial matching relation stated there. For each infrared
cutoff N the finite version of this theorem supplies a continuous
counterterm map `c_N(p)` with `|c_N(p)|<=C epsilon`. If the free dispersion
displacement is `V p+Q(p)`, with invertible V and `Q(p)=O(|p|^2)`, then

    p -> V^-1[c_N(p)-Q(p)]                               (19)

maps a fixed `O(epsilon)` cube into itself. Brouwer provides a solution.
A compact subsequence of centers, together with diagonal limits of
every fixed-scale running history, removes the infrared cutoff and
preserves the exact physical matching relation.

This proves existence of the needed Fermi-point adjustment without
derivatives in p and without independent microscopic chemical-potential
parameters. The shifts stay inside every fixed certified momentum
margin. A small derivative of the counterterm map would yield a stronger
uniqueness statement, but is not required for existence. In particular,
relative-coordinate moments alone need not control derivatives of
branch-changing endpoint phases in an anchor-supremum norm; the new
argument never makes that inference.

## 8. Progress toward the GN hypotheses, and the unresolved part

For the replacement endpoints,
DEGREE_REPARAMETERIZATION proves the needed perturbative powers.
The label/degree sums themselves can be controlled: for even `p>=2`,

    max(1,p/2-1)>=1/2+p/8.

Thus an endpoint weight can reserve a factor `epsilon^(p/8)` to absorb
`C^p` field-label and Gram constants and a factor `sqrt(epsilon)` for
the endpoint count. The usual renormalized internal dimensions have
one unit of gain before the small anomalous corrections; a fixed
fraction of that gain is available for short memory. These elementary
facts remove the combinatorial degree obstacle.

They do not by themselves establish (2)–(3). The required additional
verification is the *reduced normalized beta-function comparison*,
including its history-difference bounds at fixed Fermi centers. In
particular:

- The reference contribution at constant marginal history must have
  the stated scale decay for velocities and forward couplings.
- Trees with an initial irrelevant endpoint or a curvature insertion
  must retain that gain after localization and differentiation.
- History differences must carry the extra factor `eta` in (2)–(3),
  after all wave-function normalization terms are included; coefficient
  differences may also produce the explicitly allowed `Q_h` term.
- These bounds must be uniform in a fixed neighborhood of Fermi centers,
  with continuous fixed-cutoff beta functions and correct initial
  physical matching. No Fermi-center derivative estimate is needed.

The first item is the published multichannel reference input. The
endpoint hierarchy and Taylor bounds address pieces of the second.
I have not written a complete estimate for the last two items for
arbitrary retarded endpoints. Therefore this note does **not** tag an
endpoint-only multichannel stability theorem as proved.

What is proved is that these explicit beta estimates, if established,
close the full marginal/counterterm induction with uniform constants.
There is no remaining sequence-level small-denominator problem hidden
behind a claim that a beta function is small.

## 9. Current correlations and c remain distinct outputs

Even the completed induction would first construct bounded forward
couplings, nonzero limiting velocities, and anomalous fermion
renormalizations. Gaussian chiral current correlations additionally
require source-dependent GN comparison and the reference Ward identities.
The physical density/current source kernels must undergo the same
massive matching as the action.

A full unitary CFT limit or an `o(T^2)` pressure comparison is stronger
still. Normalized current correlations do not retain scalar vacuum
terms. Nothing in the recursion theorem identifies the stress tensor
or the thermal central-charge coefficient. Those outputs remain
separate obligations rather than consequences of bounded `g_h` alone.
