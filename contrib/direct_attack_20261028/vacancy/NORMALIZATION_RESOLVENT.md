# Eliminating dependent normalization ratios without losing short memory

Status: **[proved sequence and analytic-mark lemmas]**. The hypotheses
are raw GN estimates in which the small normalization increments are
temporarily independent variables. This note solves their dependence
on the running history; it does not assume that dependence is absent.
It complements infrared/MARKED_COVARIANCE.md and GN_MARK_CALCULUS.md.

Author: OpenAI GPT-6-family Codex vacancy worker, 2026-10-28. AI-assisted;
no human review.

## 1. An explicit triangular convolution inverse

For h<=k<=0 set

    W_a(h,k)=2^[-a(k-h)],           a>0.                 (1)

The same definitions apply first to any finite scale interval [-N,0].
For sequences with the supremum norm,

    ||W_a||<=1/(1-2^-a).                                 (2)

Let q>=0 and q<1-2^-a. The positive Neumann series of I-qW_a converges.
Writing n=k-h, its entries are exactly

    [(I-qW_a)^-1-I](h,k)
        =q (1-q)^[-(n+1)] 2^(-a n),                     (3)

and

    [(I-qW_a)^-1 W_a](h,k)
        =(1-q)^-1 2^(-a' n),
    a'=a+log_2(1-q)>0.                                 (4)

**Proof.** For m>=1 the m-fold convolution has coefficient
`binom(n+m-1,m-1)2^(-a n)`, obtained by distributing the nonnegative
distance n among m factors. Summing q^m uses
`sum_(r>=0) binom(n+r,r)q^r=(1-q)^(-n-1)`. This proves (3);
the same sum with one final W_a proves (4). Finite-interval matrix
entries have the same formula whenever h,k lie in that interval.
All constants are independent of N. QED.

Thus a small diagonal-and-triangular feedback loses only O(q) in the
short-memory exponent. It need not be handled by an uncontrolled sum
over the number of scales.

## 2. Raw normalization equations and their analytic solution

Let z_h be a finite vector of single-step normalization increments,
such as the coefficient that produces `Z_(h-1)=Z_h(1+z_h)`.
The vector may include both temporal and spatial self-energy derivatives.
Then the next velocity is computed from
`v_(h-1)=(v_h+z_(h,1))/(1+z_(h,0))`; a raw covariance written using that
next velocity is a function of the present velocity and z_h, not an
unexplained dependence on a future running history.
Let x denote the independent running interaction and velocity histories.
They range over a complex neighborhood, with interaction amplitudes
of magnitude O(eta), and velocities in a fixed complex neighborhood
of nonzero real values. The raw normalized GN construction gives

    z_h=Phi_h(x,z),                                      (5)

where Phi_h depends only on coordinates at scales >=h. We allow z_h
itself on its right side. Assume on a complex ball |z_h|<=R eta^2:

    |Phi_h(x,0)|<=C_0 eta^2,
    ||partial_(z_k) Phi_h||<=q W_a(h,k),                 (6)

with q=C_1 eta or any smaller bound. All functions are analytic and
these bounds hold uniformly on a slightly larger history domain than
the one ultimately used. Initial irrelevant contributions of size
epsilon^2=eta^4 are included in the first bound.

In the raw covariance construction it is useful to prove the bounds on
a larger fixed ratio domain `|z_h|<r_0`, where r_0 is a sufficiently
small constant independent of eta. This retains a fixed dimensional
reserve and permits raw z marks without an eta^-2 Cauchy loss. A z mark
changes a covariance or normalization factor but removes no interaction
endpoint, so the normalization self-energy often gives the stronger
q=O(eta^2). Restricting afterwards to the invariant O(eta^2) ball in
(6) recovers the small physical anomalous exponent.

Choose R=2C_0, and then eta so that
`q/(1-2^-a)<=1/2`. Equations (2) and (6) imply that Phi maps the z ball
to itself and is a contraction. It has a unique fixed point z(x),
uniformly bounded by R eta^2. The Picard iterates are analytic and
converge uniformly on each smaller complex history domain, so z(x)
is analytic at every finite cutoff. No limit in the anomalous product
Z_h itself is asserted.

The same construction can be made directly on bounded infinite
sequences when Phi is defined there. Alternatively the finite-cutoff
solutions are consistent on each fixed higher-scale block, by
causality and uniqueness, and define the infinite sequence. Constants
are independent of the number of scales.

## 3. One scale-specific history mark

Let x_l be a local-amplitude history coordinate at scale l. Suppose,
with z held independent,

    ||partial_(x_l) Phi_h||<=C_2 eta W_a(h,l),  h<=l,     (7)

and the derivative is zero for h>l. The norm also covers a fixed
finite number of branch components. Differentiate (5), or differentiate
its uniformly convergent Picard iterates. Taking absolute values gives

    |D_(x_l) z_h|<=C_2 eta W_a(h,l)
                +q sum_(k=h)^l W_a(h,k)|D_(x_l)z_k|.

By (4),

    ||D_(x_l) z_h||
        <=C_2 eta (1-q)^-1 2^[-a'(l-h)].                 (8)

This preserves both the extra eta factor and exponential short memory.
The independent-variable bound (7) must be proved from the raw tree
count. For the normalization self-energy, the purely local contribution
has at least two interaction endpoints: a one-local-quartic tadpole
is momentum independent and has no wave-function derivative. An
initial irrelevant endpoint supplies its eta^4 budget instead.

Velocity history marks can obey the same or a stronger estimate, since
marking a velocity does not remove an interaction endpoint. The
canonical relevant term is not part of Phi.

## 4. A second common-history mark by one-variable Cauchy

Use the common-history coordinate u of GN_MARK_CALCULUS, formally
independent of the deviations d_k. Keep all physical initial irrelevant
endpoint kernels fixed. Suppose (6)--(8) hold uniformly when u is
shifted by a complex scalar t times any unit vector, |t|<=c eta, on a
larger domain. Cauchy's estimate applied to the analytic function
`D_(x_l) z_h` gives on the inner domain

    ||D_u D_(x_l) z_h||
          <=C_3 2^[-a'(l-h)].                            (9)

Only one auxiliary complex parameter is used. There is no factor equal
to the number of RG scales. If desired the bound can be written with
any fixed exponent below a', so that q is chosen uniformly over the
whole complex neighborhood. Direct two-variable differentiation of
Phi leads to the same result but has more terms.

The fixed initial endpoint data are essential: this is a derivative
of the formal running-history dependence, not of the microscopic
Hamiltonian or the Fermi-point anchors.

### Real C2 alternative for smooth moving cutoffs

Complex analyticity is unnecessary if the raw construction supplies
the corresponding first and mixed real derivatives. This is useful
when a smooth cutoff depends on a real running velocity. Write
`A=D_z Phi`, `P_l=partial_(x_l)Phi`, and `a_l=D_(x_l)z`.
For a common physical-history direction u, first solve

    (I-A)D_u z=partial_u Phi.

The raw bound `||partial_u Phi||_infinity<=C eta` and (2) give
`||D_u z||_infinity<=C' eta`. Hence differentiation along the solved
history is a bounded common direction U on the combined raw (x,z)
coordinates. Suppose the raw mixed estimates, valid for such U, are

    |D_U P_(h,l)|<=C W_a(h,l),
    |D_U A_(h,k)|<=C eta W_a(h,k).                       (9a)

The second estimate retains one interaction factor: a ratio mark
removes none and a common amplitude mark removes at most one from the
two-endpoint normalization self-energy. Differentiating
`(I-A)a_l=P_l` gives the exact equation

    (I-A)D_u a_l=D_U P_l+(D_U A)a_l.                    (9b)

Use (8) on its second term, the elementary convolution estimate in
Section 5, and then the positive inverse (3). For any fixed exponent
below a', this proves

    |D_u D_(x_l) z_h|<=C 2^[-a''(l-h)],
    0<a''<a',                                          (9c)

after choosing a fixed gap of exponents and reducing eta if needed.
All constants are uniform in the scale interval. For example one may
use successive reserves a, 3a/4, a/2, a/4. This is a real C2 implicit
differentiation argument, not a holomorphic continuation of a smooth
cutoff. It yields the same mixed-mark conclusion as (9), with a smaller
fixed exponent.

## 5. Passing the marks to the eliminated beta function

Let a raw marginal or relevant beta be B_h(x,z), and put
`B_h^red(x)=B_h(x,z(x))`. Suppose its raw derivatives satisfy

    ||partial_(x_l) B_h||<=C_B eta W_a(h,l),
    ||partial_(z_k) B_h||<=C_B W_a(h,k).                 (10)

The first estimate may include the additional UV-decaying marginal
history term allowed in Section 2 of GN_MARK_CALCULUS. Treat that term
separately as there. These raw beta bounds are needed on the solved
O(eta^2) z domain and its analytic history neighborhood. Unlike the
raw covariance construction they need not hold with the eta factor
throughout a fixed z ball. The chain rule and (8) give the bound

    ||D_(x_l) B_h^red||
      <=C_B eta W_a(h,l)
       +[C_B C_2 eta/(1-q)]
          sum_(k=h)^l W_a(h,k) W_(a')(k,l).              (11)

Take eta smaller so a'>=3a/4. For n=l-h the sum is at most
`(n+1)2^(-3a n/4)`, hence at most C_a 2^(-a n/2).
The constant C_a is independent of q and eta. Consequently

    ||D_(x_l) B_h^red||<=C'_B eta 2^[-(a/2)(l-h)].       (12)

Cauchy in the same common-history parameter, on an eta-sized complex
disk, gives

    ||D_u D_(x_l) B_h^red||<=C''_B 2^[-(a/2)(l-h)].      (13)

In the real C2 alternative, directly differentiate the chain rule. The
raw mixed bounds `|D_U partial_(x_l)B_h|<=C W_a(h,l)` and
`|D_U partial_(z_k)B_h|<=C W_a(h,k)`, together with (8), (9c) and (10),
give the same O(1) mixed-mark estimate after another fixed reserve of
the exponent. Only convolutions of exponentially decaying kernels
occur. Thus a smooth velocity-dependent cutoff can be treated entirely
on real histories when its direct one- and two-mark covariance bounds
are supplied.

Thus eliminating the history-dependent normalization ratios preserves
the exact one- and two-mark forms needed in GN_MARK_CALCULUS. Choose
theta so the final reserved exponent is at least 4theta (a/2 in the
analytic argument, or a smaller fixed fraction in the real C2 argument).
This requires only a smaller
fixed theta; no optimization of the exponent is needed.

The z derivative in (10) deliberately does not require an extra eta.
If z includes the spatial derivative, the velocity beta
`(z_1-v_h z_0)/(1+z_0)` has an O(1) z derivative. Equation (8) supplies
the needed eta after eliminating z. A quartic normalization factor
instead differentiates a coupling of size eta, and normalizing the
relevant canonical term supplies nu_h=O(eta), so those components obey
stronger bounds. These distinctions prevent silently dropping the
direct velocity dependence on the solved self-energy derivatives.

## 6. Scope of the result

Equations (6), (7), and (10) are raw marked-tree hypotheses, not
consequences of a bound on |z_h| alone. The direct covariance proof
must provide the Gram and moment estimates under independent ratio
marks; the tree count must then supply the stated endpoint factors
and path weights. Retarded initial kernels may enter Phi but are
fixed data with their own eta^4 and field-degree budgets.

Subject to those explicit raw estimates, the dependence of the
normalization increments on the running history is now completely
solved. It does not add another open fixed-point or differentiated
normalization problem. Combining (12)--(13) with the constant-history
reference comparison gives the corrected G1/G2 through
GN_MARK_CALCULUS, followed by COUNTERTERM_INDUCTION and
FERMI_SHIFT_EXISTENCE. Source correlations, pressure and conformal
identification are outside this statement.
