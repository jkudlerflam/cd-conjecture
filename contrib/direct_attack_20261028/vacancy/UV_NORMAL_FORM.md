# A finite-depth ultraviolet normal form with a first spatial moment

## 1. Statement and scope

**[proved, given the standard finite-range Lieb-Robinson bound]** Consider
a one-dimensional finite interval or ring with an arbitrary on-site
Hilbert space. Let

    H0=sum_x a_x,       V=sum_x v_x,
    ||a_x||<=M0,        ||v_x||<=M1,

where every term is Hermitian and supported within a fixed distance `r0`
of its anchor `x`. A fixed finite number of terms per anchor can be
collected into one term. Neither translation invariance nor a spectral
gap is assumed. The constants below can depend on `r0,M0,M1` but not on
the ring length, boundary conditions, cell pattern, or any low-energy gap.

For a specified decomposition into interactions define

    ||Phi||_1 = sup_y sum_(X contains y)(1+diam X)||Phi(X)||.       (2)

The diameter uses graph distance; using the size of a supporting interval
only increases this norm. We construct the decomposition, so there is no
implicit infimum over alternative decompositions.

For every `0<delta<=delta0` and every real `s`, there is a finite-depth
unitary `U_delta(s)` and a local Gaussian-filtered interaction `W_delta`
such that

    U_delta(s)(H0+sV)U_delta(s)* = H0+s W_delta + R_delta(s),

    ||R_delta(s)||_1
       <= C s^2 delta^-18 + C |s| delta^-5 exp(-c/delta^2).       (3)

`W_delta` is defined in (8), not by an unproved low-energy projection.
The statement holds uniformly on small rings as well by taking local
regions equal to the entire ring when necessary. Large constants and the
nonoptimal exponent 18 are intentional.

For the physical normalization `H=t H0+J V`, choose

    s=J/t,             delta=s^a,             0<a<1/18.

Then

    ||t R_delta(s)||_1/J
       <= C s^(1-18a)+C s^(-5a) exp[-c s^(-2a)] -> 0.           (4)

This supplies a controlled Gaussian-filter version of the polynomial
ultraviolet estimate, replacing the unspecified exponent in
`uniform_matching` (10). Gaussian energy-transfer tails remain, so this is
not an exact compact spectral cutoff. It is not a low-frequency matching
theorem and makes no statement about the central charge.

## 2. Filter, locality input, and a strictly local generator

Let `tau_u(B)=exp(iuH0) B exp(-iuH0)` and put

    g_delta(u) = delta/(2 sqrt(pi)) exp(-delta^2 u^2/4),
    F_delta(u) = -(i/2) sign(u) erfc(delta |u|/2).                (5)

The convention is `integral g_delta(u) exp(i omega u)du
=exp[-(omega/delta)^2]`. Direct integration by parts gives

    integral F_delta(u) exp(i omega u)du
       = (1-exp[-(omega/delta)^2])/omega,

with zero value at `omega=0`, and in distributions

    F_delta'=-i delta_Dirac+i g_delta.

Thus the formal full generator `S_delta=integral F_delta(u)tau_u(V)du`
is anti-Hermitian and satisfies

    [S_delta,H0]=-V+integral g_delta(u)tau_u(V)du.                (6)

We will not conjugate by this extensive quasi-local generator. Instead,
set `T=delta^-2`, and choose an integer

    R >= r0+2+(2 v_LR+2)T,                R<=C_R delta^-2,

where `v_LR` is a finite-range Lieb-Robinson velocity, increased to at
least one if necessary. Let `I_x` be the radius-`R` interval about `x`;
let `H0,I_x` contain precisely the `a_y` whose supports lie inside `I_x`.
Define

    s_x = integral_(-T)^T F_delta(u) tau_u^(I_x)(v_x)du,
    S_R = sum_x s_x.                                           (7)

Every `s_x` is anti-Hermitian, has support `I_x`, and obeys
`||s_x||<=C M1/delta`, because `integral |F_delta|=C/delta`.

The local filtered interaction in the theorem is

    W_delta = sum_x w_x,
    w_x = integral_(-T)^T g_delta(u) tau_u^(I_x)(v_x)du.          (8)

In particular, `||w_x||<=M1` and its support is `I_x`. This bound does not
assert that `w_x` is close to its vacuum expectation.

The only dynamical estimate used below is the standard finite-range
Lieb-Robinson estimate. In the form needed here it says that a bounded
observable initially near the center of `I_x` has a commutator with the
bounded boundary-crossing terms of `H0` at most

    C M1 exp[-mu(R-r0-v_LR |u|)]                              (9)

up to a fixed change of constants. The same estimate compares full and
restricted dynamics at times shorter than the distance to the boundary.
This estimate follows from the usual iterated-commutator/path bound for
bounded finite-range interactions. Its constants do not involve the
on-site dimension. The phase-model hopping and contact terms are bounded,
so infinite-dimensional gap coordinates cause no difficulty here.

For literature context, the locality machinery, including
infinite-dimensional local spaces with bounded interactions, is discussed
in Bachmann, Michalakis, Nachtergaele and Sims,
[arXiv:1102.0842](https://arxiv.org/abs/1102.0842). The abstract was opened
for this contribution. We use the standard local-dynamics estimate (9),
not that paper's gapped ground-state equivalence theorem.

## 3. First-order cancellation with an exponentially small local error

Write `H0=H0,I_x+H_out+H_boundary`, where `H_out` commutes with `s_x`
and `H_boundary` is a fixed number of terms crossing the two boundaries.
Since

    [tau_u^(I_x)(v_x),H0,I_x]=i d tau_u^(I_x)(v_x)/du,

integration by parts gives the exact identity

    v_x+[s_x,H0]
       = w_x + e_x,

    e_x = i F_delta(T)tau_T^(I_x)(v_x)
          -i F_delta(-T)tau_(-T)^(I_x)(v_x)
          +integral_(-T)^T F_delta(u)
                      [tau_u^(I_x)(v_x),H_boundary]du.          (10)

The error has support in `I_x` enlarged by the original interaction range.
The two endpoint terms are bounded by
`C M1 exp[-1/(4 delta^2)]`, using the elementary `erfc` tail bound.
By (9), the boundary term is bounded by
`C delta^-1 exp[-c delta^-2]`. Consequently

    ||e_x||<=C delta^-1 exp[-c/delta^2].

A site belongs to at most `O(R)` supporting intervals, and their diameter
is `O(R)`. Assigning `e_x` to this interval therefore gives

    ||sum_x e_x||_1
       <= C R^2 delta^-1 exp[-c/delta^2]
       <= C delta^-5 exp[-c/delta^2].                           (11)

This argument directly bounds the actual local error. It does not infer a
local interaction bound from a finite-volume operator norm.

For clarity about the word "filter," the ideal filtered anchored operator
is `w_x,ideal=integral_R g_delta(u)tau_u(v_x)du`. The omitted time tail
and the comparison of restricted to full dynamics show

    ||w_x-w_x,ideal||<=C exp[-c/delta^2]                        (12)

after changing constants. Equation (12) is an anchored operator statement;
it is not used as an interaction-norm bound. The theorem uses the strictly
local definition (8) and the stronger, explicitly decomposed error (11).
In finite volume (6) shows the ideal filter multiplies a matrix element
with energy transfer `omega` by `exp[-(omega/delta)^2]`. The local version
implements that filter up to (12), with a volume factor if anchored terms
are summed in ordinary operator norm.

## 4. A circuit whose Taylor remainder can be counted directly

Color the intervals `I_x` so intervals in the same color are disjoint.
The overlap graph has maximum degree at most `4R+O(r0)`; a greedy coloring
therefore uses

    q<=C R

colors. On an interval an explicit residue coloring also suffices. Define

    U_c(s)=product_(x of color c) exp(s s_x),
    U_delta(s)=U_(q-1)(s)...U_0(s).                            (13)

Each factor is unitary and each color layer is a product of disjoint local
gates. At `s=0`,

    d[U_delta(s) B U_delta(s)*]/ds = [S_R,B].                   (14)

No Trotter approximation or limit is used: (13) is the chosen exact
change of basis.

For comparisons between different material patterns, fix the coloring and
the layer order using the lattice geometry alone and use that same choice
in every pattern. Independently recoloring two patterns could change the
description far from the defect and invalidate an attempted cancellation
of their local interactions. A common coloring avoids that problem.
It guarantees equality of the transformed terms only where the original
coefficient arrays agree at the same anchored positions throughout the
relevant causal cone. A geometrically fixed coloring is not generally
translation covariant. A translated or rephased selected bulk requires a
separate comparison; common coloring alone does not identify its terms
with those of the original selected phase.

For a term `B_x` of original range `r0`, its causal cone through the
`q` layers lies in an interval of diameter

    D<=C R q+C r0<=C R^2.                                    (15)

Within a layer the relevant centers are separated by more than `2R`.
An interval of length `D` consequently meets at most `C(D/R+1)` gates of
that color. The whole causal cone contains at most

    n_gate<=C q(D/R+1)<=C R^2                                (16)

gates. Gates outside the causal cone cancel from the conjugation for all
`s`, so the conjugation is exactly a product of at most this many unitary
conjugations. This argument works also on a ring when the cone covers the
whole ring: use `min(N,C R^2)` for its length. If `N<=2R+O(r0)`, use one
gate per color, giving an even smaller bound.

For anti-Hermitian generators `b_1,...,b_n` and any bounded `B`, repeated
product differentiation of a product of unitary conjugations gives,
uniformly at every real `s`,

    ||d tau_circuit,s(B)/ds||
       <= 2(sum_j ||b_j||)||B||,

    ||d^2 tau_circuit,s(B)/ds^2||
       <= 4(sum_j ||b_j||)^2 ||B||.                            (17)

One way to check the second inequality is to view conjugation by each gate
as `exp(s ad_(b_j))`. Each such map is an isometry in operator norm;
each derivative inserts `ad_(b_j)` of norm at most `2||b_j||`. The sum
over all ordered pairs of derivative locations, including two derivatives
at the same gate, is bounded by the square in (17).

Equations (7) and (16) imply

    sum_(causal cone) ||s_x|| <= C R^2/delta.                  (18)

Taylor's theorem with an integral remainder, (14), (17), and (18) now give

    U_delta(s)a_x U_delta(s)*
       = a_x+s[S_R,a_x]+r_x^(0),
    ||r_x^(0)|| <= C s^2 R^4 delta^-2,

    s U_delta(s)v_x U_delta(s)*
       = s v_x+r_x^(1),
    ||r_x^(1)|| <= C s^2 R^2 delta^-1.                        (19)

Each remainder in (19) is supported in the explicit interval from (15).
Subtracting the first derivative does not enlarge that interval.
The estimates are uniform in `s`: no analytic radius in volume is hidden
in the argument.

## 5. Interaction norm and the exponent 18

Assign each remainder in (19) to its anchored interval of length `O(D)`.
Any site belongs to at most `O(D)` such intervals, and the first diameter
weight is at most `O(D)`. Thus

    ||sum_x(r_x^(0)+r_x^(1))||_1
       <= C D^2 [s^2 R^4 delta^-2+s^2 R^2 delta^-1]
       <= C s^2 [R^8 delta^-2+R^6 delta^-1]
       <= C s^2 delta^-18,                                   (20)

using `D=O(R^2)`, `R=O(delta^-2)`, and `delta<=1`.
Summing (19) and substituting the exact first-order identity (10) yields

    U_delta(s)(H0+sV)U_delta(s)*
       =H0+s W_delta
          +sum_x(r_x^(0)+r_x^(1))+s sum_x e_x.

Equations (11) and (20) prove (3). The weak condition `a<1/18` then proves
(4). Improving the cutoff choice, the circuit layout, or the norm estimates
can greatly improve the exponent; none is needed for the existence of an
`o(J)` ultraviolet error in the weighted norm.

## 6. What this contributes to vacancy matching

The phase-model vacancy bath has bounded nearest-neighbor hopping. Its
contact projector is bounded by one. With fixed bounded microscopic
material interactions the assumptions of Section 1 hold in squeezed
coordinates, uniformly over material sector patterns. The proof uses no
walker gap and is therefore not invalidated by a long wrong cell's small
level spacing. Applied there, (4) controls the ultraviolet remainder in a
norm that can pay for finitely many interfaces after identical bulk
interactions have actually been canceled.
For this particular circuit, cancellation follows directly only for
coefficient arrays equal at the same anchored positions away from the
modified region. Translated selected domains are not covered by that
observation without an additional argument.

There are still two essential restrictions.

1. The local filtered interaction `W_delta` is not small in operator norm
   relative to its vacuum expectation. An energy-transfer filter does not
   project onto low absolute energies. Controlling its ground-state
   influence still requires a relative-form estimate or a ground-state RG
   calculation. The endpoint/contact-gradient identity from
   `uniform_matching` supplies extra structure but is not used or proved
   sufficient here.
2. A weighted bound on the change-of-basis remainder does not establish a
   weighted bound on the connected *ground-state energy* of `H0+sW_delta`.
   One must still prove the spatially uniform, inhomogeneous low-frequency
   estimate proposed in `uniform_matching` (12), or a different one-sided
   sector comparison. Subtracting bulk terms is not automatic.

The statement therefore removes one technical ultraviolet task from the
list of missing inputs. It does not close the low-frequency matching,
the dressed-defect comparison, or the CFT inference for the `d=24` chain.
