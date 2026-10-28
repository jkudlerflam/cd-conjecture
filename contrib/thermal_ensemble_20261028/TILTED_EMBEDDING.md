# Length tilts in the auxiliary representation

Status: **[proved abstract bounds; microscopic application conditional on the
stated cell data and bounded-interaction representation]**. No infrared or
central-charge conclusion is asserted.

Author: coordinating OpenAI GPT-6-family Codex agent, 2026-10-28. The
length-preserving embedding was independently suggested by the coordinator
of `direct_attack_20261028`. No human review.

Independent AI review: the thermal worker and the direct-attack coordinator
rederived (6)-(11), including the partial-trace direction. Finite-matrix
controls in `check_tilted_embedding.py` pass 18 thermal comparisons and
30 derivative bounds; output is saved in `tilted_embedding_checks.json`.

## 1. Hypotheses

Let a cell space be the orthogonal sum `K = direct_sum_(ell>=1) K_ell`, with
length operator `l`. Let a positive self-adjoint `A` commute with length and
grading. Its kernel is a code `P K` of dimension `d_P=2^(n+1)`, identified
with `n` payload fermions and one length fermion. In the code,

    l_P = b + N_len,             b>=1,

and the two length eigenspaces have equal rank. Write `Q=1-P`. Assume

    A >= Delta Q,                         Delta>0,
    |l-b| Q <= b Q + A/gamma,              gamma>0,       (1)
    Tr_Q exp(-t0 A) < infinity

for some `t0>0`. All inequalities in (1) are joint spectral inequalities,
since `A` and `l` commute. A linear cost for wrong cell lengths, together
with finitely many exceptional lengths, supplies (1). In C44 the selected
lengths are `b,b+1`; the length-one exception has `|l-b|<=b`; the stated
wrong-length cost supplies the other cases. This is an application of the
imported all-length bound, not a new proof of that microscopic bound.

The full cell Hamiltonian is `sum A_j + sum B_X`, where the interactions
are bounded and even, conserve total physical length, and have finite
interaction strength. Unbounded positive readout terms in invariant illegal
blocks must be included in `A`, if this decomposition is used there. The
bounded-interaction hypothesis is not inferred from a bounded physical bond
after cells of arbitrary length are grouped.

## 2. An embedding that preserves physical length

Let `E=C Omega direct_sum QK`. Embed the original space in
`F_(n+1) graded_tensor E` by

    J G|u> = |u> tensor Omega,
    J|s> = |0> tensor |s>,             |s> in QK.         (2)

Choose logical vacuum `|0>` in the code's length-`b` subspace. Define

    a Omega=0,          a|QK=A|QK,
    q Omega=0,          q|QK=(l-b)|QK,
    l_ext=b+N_len+q.                                     (3)

These definitions give `l_ext J=J l` and `(I tensor a)J=JA` exactly.
Thus length-conserving operators stay length conserving when extended as
`J_X B_X J_X^dagger`. Total length conservation is not replaced by an
unrelated logical particle-number conservation law.

For a real length tilt `s`, define `a_s=a+s q`. On the valid subspace,

    (a_s+s(b+N_len)) J = J(A+s l).                       (4)

The unbounded part of the tilt belongs to the auxiliary semigroup; the
bounded logical term `s N_len` belongs to the logical reference Hamiltonian.
The scalar `sb` is retained in the partition function. Expanding `s l` as
a bounded perturbation on the all-length cell space would be invalid.

## 3. A uniform auxiliary gap and heat trace

Put

    s0 = min(gamma/4, Delta/(4b)).                        (5)

For `|s|<=s0`, (1) gives on the excited auxiliary space

    a_s >= A-|s|(b+A/gamma)
        >= (3/4) A-Delta/4 >= A/2 >= Delta/2.            (6)

In particular

    Tr_(E minus Omega) exp(-t a_s)
        <= Tr_Q exp(-t A/2),                             (7)

so a common heat-trace time is `2t0`. This proves the needed uniform gap
and trace constants without any dependence on the logical bandwidth.
The same estimates hold for a complex tilt in a disk of this radius with
`a_s` replaced by its real part; simultaneous spectral calculus is available.

For real `s` in the interval and integers `k>=0`, one also has

    Tr_Q |d_s^k exp(-t a_s)|
       <= k! (tb+4/gamma)^k Tr_Q exp(-t A/4).             (8)

Indeed the derivative is `(-t q)^k exp(-t a_s)`. Reserve one factor
`exp(-tA/4)` from (6), expand `(tb+tA/gamma)^k`, and use
`x^j exp(-x)<=j!` with `x=tA/4` on the remaining factor. The binomial sum
is bounded by `k!(tb+4/gamma)^k`. Equation (8) applies when the trace on
its right is finite, for example `t>=4t0`. It gives uniform fixed-order
derivative control of coarse time blocks. It does not prove differentiated
infrared RG estimates.

Any auxiliary matching theorem with constants depending only on a lower
auxiliary gap, a coarse-time heat trace, and the bounded coefficient norm
can therefore use common constants on `|s|<=s0`. Applying a particular
matching theorem still requires all its other hypotheses and its proof.

## 4. Invalid embedding states at nonzero tilt

The complement of `J K` at one site is

    I = (logical states orthogonal to |0>) tensor QK.

Extend each interaction by `J_X B_X J_X^dagger`, so it vanishes if any
site in its support is invalid. The set of invalid sites is then invariant.
Take

    w = sup_j sum_(X containing j) ||B_X||.

Let `Z_valid(s)` denote the original full-cell partition function and
`Z_enl(s)` the enlarged-space trace with onsite terms from (4). Define

    z_I(beta,s) = Tr_I exp[-beta (a+s(b+N_len+q))],
    zhat_I(beta,s) = exp[beta s(b+1/2)] z_I(beta,s).

For any finite number `M` of cells with finite traces,

    1 <= Z_enl(s)/Z_valid(s)
      <= [1+ exp(beta w) zhat_I(beta,s)/d_P]^M.          (9)

Proof. Fix a set `S` of invalid cells. All interactions touching `S` vanish,
so that block's trace is `z_I^|S|` times the remaining valid-cell trace
with those interactions cut. In the original valid Hilbert space, compress
the cells in `S` to their code and take their normalized partial trace.
Scalar convexity of the exponential in an eigenbasis, twice applied, gives

    Z_valid(s) >= d_P^|S| Tr_rest exp[-beta H_av].

This is the same elementary partial-trace inequality as C46 (8), used in
the reverse direction to compare with a compression. The averaged onsite
tilt on each erased code factor is `s(b+1/2)`. The averaged touched
interactions have total norm at most `w|S|`. Ordered-eigenvalue comparison
therefore bounds the remaining cut trace by

    Z_valid(s) d_P^(-|S|) exp[beta(s(b+1/2)+w)|S|].

Multiplication by `z_I^|S|` and summation over all subsets proves (9).
The all-valid block supplies the lower bound. This proof permits arbitrary
logical interactions and does not assume a logical spectral gap.

On an invalid product state, the energy inside the centered trace is

    A+s(l-b+N_len-1/2).

Equations (5)-(6), `N_len` in `{0,1}`, and `b>=1` show this is at least
`A/2-|s|/2 >= A/4`. Consequently

    zhat_I(beta,s) <= (d_P-1) Tr_Q exp(-beta A/4).        (10)

For `beta>=4t0`, its right side is at most

    (d_P-1) Tr_Q exp(-t0 A)
                       exp[-(beta/4-t0) Delta].        (11)

If `w<Delta/4`, (9)-(11) bound the enlarged/valid free-energy difference
per cell by an exponentially small function of beta, uniformly for
`|s|<=s0`. Thus this embedding artifact does not alter an already-proved
low-temperature power coefficient. The condition on `w` is sufficient
and deliberately inefficient; it can be met by choosing the final small
coupling after the fixed cell data. This is an embedding comparison, not
a proof that arbitrary physical illegal sectors have been included in
the chosen cell representation.

## 5. Role in the conjecture program

The physical-pressure transfer uses a zero of the rooted pressure as a
function of the length tilt. At fixed `n`, (5) is a positive interval.
After the coexistence energy is subtracted, the retained energy scale is
`W=O(lambda^2)`; a proof placing the pressure zero in this interval can
choose lambda last. The accompanying pressure-root theorem must supply
that placement rather than assume it.

Together, (4)-(11) remove a specific possible obstruction: the physical
length chemical potential need not enter an auxiliary expansion as an
unbounded source, and its derivative/heat-trace constants can be uniform
in a neighborhood of the pressure zero. They do not prove the thermal
vacuum-root RG estimate, Gaussian scaling, or any value of microscopic
Virasoro central charge.
