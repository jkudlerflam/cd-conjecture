# Applying massive matching to the positively completed mobile chain

Status: **[proved parameter reduction using the recorded C44/C45 finite
construction and the proved two-charge fermionic representation]**. This note
instantiates the strong norm, all-length, tilt and pressure-zero hypotheses.
The logical infrared theorem remains open. The separate Jordan-Wigner
dictionary is supplied by [SEPARATE_JW.md](exact/SEPARATE_JW.md), including
auxiliary charge exchange; it is not inferred from an auxiliary gap.

Author: coordinating OpenAI GPT-6-family Codex agent, 2026-10-28.
AI-assisted research, without human review.

## 1. Fixed microscopic data and the repaired bond

Fix one integer n and every parameter preceding lambda in C44/C45's
ordered schedule. Set N=n+1 and b=n+3. Use the positive completion
specified in [ILLEGAL_SECTORS.md](ensemble/ILLEGAL_SECTORS.md):

    h_completed = h_C44
       + (|L#><L#|+|##><##|) tensor (a_B tensor I+I tensor a_A).

It is one fixed Hermitian nearest-neighbor bond on 4464 letters.
It agrees with the original on the complete selector-legal invariant
space. The full rooted Hilbert space includes every word between markers.
After subtracting the fixed coexistence energy e times physical length,
the completed unperturbed cell operator has

    ker A=P,     dim P=2^N,     l_P=b+N_len,
    A>=Delta_f Q,
    A|K_ell >= kappa(ell+1) Q_ell,                      (1)

where Delta_f and kappa are the explicit positive constants in that
proof. Its heat trace obeys

    Tr_Q exp(-t A)
       <=d exp(-2t kappa)/(1-d exp(-t kappa)), d=4464.  (2)

The positive readout on cells with two or more tokens is absorbed as
`A_lambda=A+lambda^2 C`, with C>=0. Its invariant-block proof shows
that the code-connected SW calculation is unchanged. Equations (1)-(2)
hold uniformly with A_lambda in place of A. All remaining interactions
are bounded on the full, infinite direct-sum cell space.

## 2. The norm exponent is chosen before lambda

Fix an output spatial exponent a_out, a logical field weight H>=2, and
an input support exponent

    L >= 4(a_out+log 256+4),
    a = L+N log(8 H^2).                                (3)

Use this a in C45's local SW construction from the outset. At fixed n
all its constants are finite. That construction supplies a positive
lambda_SW and constants K_B,K_R,K_P such that, for 0<lambda<=lambda_SW,

    U (H_completed-e L_M) U^dagger = D_lambda+R_lambda,
    D_lambda=A_lambda+B_lambda,
    ||B_lambda||_a <= K_B lambda,
    ||R_lambda||_a <= K_R lambda^3,
    ||H_P||_0 <= K_P lambda^2,
    H_P=P_all D_lambda P_all.                          (4)

Here all scalar port shifts, readout terms and weak marker chemical
adjustments are included. K_P is not inferred from the tiny nonscalar
bandwidth. The source's explicit finite-order expressions give (4)
after increasing the constants and restricting lambda<=1. These are
local interaction norms, not extensive operator norms independent of
volume. At M cells they imply `||H_P||<=M K_P lambda^2`.

No uniformity in n is required. Increasing a in (3) can drastically
decrease the allowed lambda; it does not invalidate an existence
argument at fixed n. Every quantity in (3)-(4) is fixed before the
eventual thermodynamic limit.

## 3. Dimension-independent coefficient norms

Use the finite logical algebra and the operator-valued auxiliary
embedding of the [massive theorem](../direct_attack_20261028/infrared/INFINITE_CELL_PROGRESS.md).
Its local matrix-unit estimate on r cells is

    sum_I H^(degree I) ||O_I|| <= (8 H^2)^(Nr) ||O||.   (5)

The representation preserves interval support and both exact charges.
SEPARATE_JW identifies the payload-plus-length quadratic reference and
preserves operator norms, including vertices exchanging charge with
auxiliaries. Equation (5) concerns the
finite logical dimension, not the dimension of the excited space.

Extend H_P independently of the auxiliary algebra in CAR variables,
using the isometric extension E_JW from that dictionary, and put
`K_D=J B_lambda J^dagger-E_JW(H_P)`. In the original tensor variables
this extension includes the required bounded auxiliary parity factors;
a naive identity extension there is incorrect. Since D_lambda is locally block
diagonal, every logical coefficient of K_D annihilates the auxiliary
vacuum from both sides. Combining (3)-(5) gives

    k:=||K_D||_(L,H) <=2 K_B lambda,
    r:=||J R_lambda J^dagger||_(L,H) <=K_R lambda^3.     (6)

Thus the large code dimension is paid by the fixed support exponent,
and does not require a finite cutoff on excited cell lengths.

## 4. One uniform interval of physical length tilts

The uniform coercivity implies `|l-b|Q<=b Q+A_lambda/kappa`.
Take

    s0=min(kappa/4,Delta_f/(4b)),
    t0=log(2d)/kappa,
    tau>=max(4L/Delta_f, 2t0).                          (7)

The construction in [TILTED_EMBEDDING.md](TILTED_EMBEDDING.md) puts the
unbounded tilt in `a_(lambda,s)=a_lambda+s q_aux`, where
`q_aux=l-b` on excluded states. For |s|<=s0,

    a_(lambda,s)>=A_lambda/2>=Delta_f Q/2,
    Tr_Q exp[-tau a_(lambda,s)] <=1/(2d)<e-1.           (8)

The last bound follows from (2) at tau/2>=t0. The vacuum contributes
one. Thus a single tau works for every small lambda and every tilt in
this interval. The scalar sb and bounded logical term sN_len remain
in the logical reference exactly. Since U commutes with total length,
the same B_lambda and R_lambda in (4) apply at every s.

Use the massive theorem's constants

    v_*=1/(8tau),   B_*=100 exp(L/2)/tau,
    C_match=4B_*/v_*.                                  (9)

Any fixed harmless enlargement required by the periodization bound can
be included in B_*. These constants have no inverse logical bandwidth.
The isometric two-charge dictionary has exactly the coefficient bound
(5), so introduces no additional support or auxiliary-dimension loss.

## 5. A noncircular final coupling choice

Let W=c_n lambda^2 be the retained bandwidth convention in C45, with
fixed c_n>0. Given a desired positive normalized endpoint tolerance eta,
it is sufficient to choose lambda below all of

    1, lambda_SW,
    v_*/[16(2K_B+K_R)],
    Delta_f/[8(K_B+K_R)],
    eta c_n/(C_match K_R),                             (10)

omitting a fraction whose denominator is zero. Decrease lambda further
to satisfy the original finite C45 inequalities and

    K_P lambda^2+K_R lambda^3 <= b s0/4.                (11)

All right sides are positive quantities determined before lambda. The
first norm restriction yields k+r<=v_*/16; the second gives the
interaction threshold needed for the invalid-state comparison.
The third gives

    C_match r/W <= eta.                                (12)

The massive theorem therefore supplies a temperature-independent logical
influence action whose correction to the full second-order H_P is at
most C_match K_R lambda^3 in the stated kernel norm. Its finite-temperature
endpoint differs from the appropriate periodization by exponentially
small terms, including scalar densities. Time moments follow with a
fixed loss of the exponential time weight. This is an exact matching
statement using that representation, not a claim that the
logical action has already been integrated.

Invalid embedding sectors contribute exponentially little to the real
tilted pressure by TILTED_EMBEDDING (9)-(11), using (10). They are distinct
from physical illegal strings, which are already included in (1).

## 6. The pressure zero is inside the controlled interval

Let p_P(T,s), p_D(T,s), and p_H(T,s) be rooted pressures per cell, with
the coexistence term subtracted and all marker shifts included in (4).
Since b<=l_P<=b+1,

    p_P(T,s0)<=K_P lambda^2+T log(2^N)-b s0,
    p_P(T,-s0)>=-K_P lambda^2+b s0.                     (13)

There is an explicit uniform thermal bound on this very interval.
Set `J=||B_lambda||_0<=Delta_f/8`, by (10). Equations (7)-(8) give

    A_lambda+s(l-b-1/2)-2J
       >=A_lambda/2-Delta_f/(8b)-Delta_f/4
       >=A_lambda/8.                                   (14a)

For beta>=16t0, split `exp(-beta A_lambda/8)` into two equal
semigroup factors and apply (2) to one. The excited trace is bounded
by `(1/(2d)) exp(-beta Delta_f/16)`. The tilted comparison therefore gives

    0<=p_D(T,s)-p_P(T,s)<=delta_T,
    delta_T <=[T/(2d 2^N)] exp[-Delta_f/(16T)],          (14)

uniformly on the original fixed interval |s|<=s0. No interval reduction
or lambda-dependent heat constant is required.

Bounded perturbation and unitary trace invariance give

    |p_H(T,s)-p_D(T,s)|<=K_R lambda^3.                  (15)

At each fixed lambda satisfying (11), take T small enough that
`T log(2^N)+delta_T<=b s0/4`. Equations (11),(13)-(15) make
p_H(T,s0)<0 and p_H(T,-s0)>0. The interacting pressure-zero theorem
then applies to the full repaired physical chain, by its exact pin/seam
comparison. Its physical pressure, after subtracting eL, is the unique
zero in this same controlled tilt interval. The original pressure is
that zero minus e. This proves the required placement; it is not an
extra compressibility assumption.

## 7. Remaining logical obligations

Equations (1)-(15) provide full-space coercivity, uniform thermal matching
constants, a small logical endpoint, and the correct physical ensemble
for the explicitly completed bond. They do not turn the fixed norm
bound in (15) into `o(T^2)` at fixed lambda. That requires the interacting
logical scalar RG estimates, or a current/stress scaling construction.

The coupling hierarchy, exact charges and finite momentum margins remain
the C42/C44/C45 inputs to that logical problem. The final question is
whether the matched retarded theory retains all N=n+1 compact modes at
arbitrarily low energy with the required CFT interpretation. No answer
to that question, and no c-d violation, is asserted here.
