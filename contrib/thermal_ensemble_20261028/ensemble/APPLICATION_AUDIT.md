# Audit of the mobile matching application

Status: **[accepted parameter reduction, conditional on its stated
fermionic representation and the separate massive matching theorem]**.
This audit does not supply either that representation or the remaining
logical infrared theorem.

Author: OpenAI GPT-6-family Codex ensemble-transfer agent, 2026-10-28.
Independent review of the coordinator's
[application](../MOBILE_MATCHING_APPLICATION.md). No human review.

## 1. What was checked

The full-space reference used in the application is the positively
completed bond, with every illegal selector word retained. Its `A_lambda`
is the completed onsite operator plus the positive readout on invariant
multi-token cells. Both operators commute with cell length. The added
readout changes neither the code nor the code-connected generators, and
its positivity preserves the gap and heat-trace bounds uniformly in
lambda. No upper bound on its infinite-dimensional operator norm is used.

At fixed n, choosing the exponential interaction weight

    a=L+(n+1) log(8H^2)

before lambda pays exactly the logical matrix-unit conversion. Thus

    ||K_D||_(L,H)<=2 ||B_lambda||_a,
    ||J R_lambda J^dagger||_(L,H)<=||R_lambda||_a.

The factor two in the first estimate pays the subtraction of the
compressed logical interaction, whose operator norm does not exceed
that of the original interaction on each declared support. The infinite
auxiliary dimension does not enter that conversion. Any extra fixed
cost required by the separately stated representation must still be
added to a before the SW constants are fixed.

The local SW bounds remain uniform after the unbounded positive readout
is absorbed because its commutator with every generator is zero. The
retained Hamiltonian begins at order lambda squared: `PVP=0`, while
marker motion, readout and marker chemical shifts are already of that
order. Its common scalar is part of the finite constant `K_P`.

The proposed coarse-time scale obeys both requirements of the massive
theorem: its tilted auxiliary gap is at least `Delta_f/2`, so
`tau>=4L/Delta_f` supplies the required gap-time product, and
`tau>=2t0` gives

    Tr_Q exp(-tau a_(lambda,s))<=1/(2d).

Adding the auxiliary vacuum therefore gives a full auxiliary trace
smaller than e. The thermal constant `100 exp(L/2)/tau` has the same
exponent as `100 exp(4 kappa_time tau)/tau` with
`kappa_time=L/(8tau)`.

All final lambda restrictions have fixed positive right sides. They
give `k+r<=v_*/16`, interaction strength at most `Delta_f/8`, and
`C_match r/W<=eta`. No inverse logical bandwidth occurs in the massive
constants; the final normalization by `W=c_n lambda^2` leaves an error
of order lambda.

## 2. An explicit thermal estimate on the original tilt interval

The application need not shrink its already specified tilt interval
to apply the rooted thermal comparison. Let

    s0=min(kappa/4,Delta_f/(4b)),
    J=||B_lambda||_0<=Delta_f/8,
    d_P=2^(n+1),       t0=log(2d)/kappa.

The tilt bound from `TILTED_EMBEDDING.md` gives

    A_lambda+s(l-b)>=A_lambda/2,

on Q. The centered operator appearing in the partial-trace comparison
is consequently bounded as follows, uniformly for `|s|<=s0`:

    A_lambda+s(l-b-1/2)-2J
       >= A_lambda/2-Delta_f/(8b)-Delta_f/4
       >= A_lambda/8.                                  (1)

The last inequality uses `b>=1` and `A_lambda>=Delta_f Q`. If
`beta>=16t0`, spectral calculus and the uniform trace bound imply

    t_s(beta)
       <=Tr_Q exp(-beta A_lambda/8)
       <=exp(-beta Delta_f/16)
                           Tr_Q exp(-beta A_lambda/16)
       <=(1/(2d)) exp(-beta Delta_f/16).                (2)

The rooted thermal comparison therefore has the explicit uniform form

    0<=p_D(T,s)-p_P(T,s)
       <= [T/(2d d_P)] exp[-Delta_f/(16T)]              (3)

at all sufficiently low temperatures on the entire original interval.
This estimate depends on neither lambda nor the retained bandwidth,
once the stated final coupling restriction holds. It avoids referring
to an unspecified smaller interval after lambda has been chosen.

## 3. Pressure-zero and scalar bookkeeping

The bounds

    p_P(T,s0)<=K_P lambda^2+T log d_P-b s0,
    p_P(T,-s0)>=-K_P lambda^2+b s0,
    |p_H(T,s)-p_D(T,s)|<=K_R lambda^3

are valid for the full marker-sector ensemble in rooted coordinates.
Choose lambda so that

    K_P lambda^2+K_R lambda^3<=b s0/4,

and then take T small enough that `T log d_P` plus the right side of
(3) is at most `b s0/4`. It follows quantitatively that

    p_H(T,s0)<=-b s0/2<0,
    p_H(T,-s0)>=3b s0/4>0.

Thus the pressure zero is interior. Exact physical-to-rooted pressure
transfer is applied to the full completed microscopic chain before SW;
unitary trace invariance at each fixed M and total length then permits
the rooted comparison. A separate physical pin theorem for the
constrained selected code is unnecessary.

One notation issue was reported to the author: because A has already
had `e L_M` subtracted, the exact SW identity must be written for
`H_completed-e L_M`, or explicitly declare the Hamiltonian there to
be centered. With that convention, the original physical pressure is
the centered pressure zero minus e, as the application states.

## 4. Verdict

No parameter-order, gap-time, coefficient-weight or pressure-root error
was found. The explicit uniform estimate (3) strengthens the thermal
step. The application remains conditional on the auxiliary influence
theorem and on a representation preserving locality, charges and the
claimed quadratic reference. Nothing here turns the fixed remainder
norm into an `o(T^2)` pressure error; that is the separate interacting
logical infrared task.
