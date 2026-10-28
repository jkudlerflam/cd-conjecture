# Four trajectories and the marked thermal-history bound

Status: **[proved abstract nonlinear recursion and scalar comparison]**.
The required derivatives of the actual normalized multichannel RG maps
are hypotheses. This note does not assert their constructive-GN proof.

Provenance: OpenAI GPT-6-family Codex ensemble-transfer agent, 2026-10-28.
Derived for the reference-subtracted thermal criterion in the sibling
`thermal/` folder. No human review or microscopic CFT claim.

## 1. Exact maps in one common representation

Fix `gamma>1`, `theta>0`, `p>theta`, scales `s_h=gamma^h`, and
`H=ceil(log_gamma T)`, so `T<=s_H<gamma T`. Write

    t_h=T/s_h,        rho=gamma^(p-theta).

Let `F_h` be the exact zero-temperature normalized step, and let
`F_h^T=F_h+G_h,T` be the exact thermal step, for `H<h<=0`. Both maps act
in a common Banach-space representation of their input and output states.
Scale-dependent spaces may be used if the displayed norms already include
the rescalings. For circle kernels this includes the choice and control of
a compatible lift; one cannot subtract an unspecified circle kernel from
a line kernel and declare the difference to be small.

The same maps are used for target trajectories x and dressed-reference
trajectories y:

    x_(h-1)^0=F_h(x_h^0),    y_(h-1)^0=F_h(y_h^0),
    x_(h-1)^T=F_h^T(x_h^T),  y_(h-1)^T=F_h^T(y_h^T).    (1)

All states and the line segments and parallelograms used below lie in
a common domain on which the derivative bounds hold. A convex domain
with the necessary room for the translated points suffices. Existence
inside that domain is an assumption or a separate bootstrap obligation.

Assume throughout that

    ||D F_h||<=L<rho,
    ||D^2 F_h||<=M,
    ||D G_h,T||<=g t_h^p.                               (2)

All constants are independent of T, h and the number of steps. It is not
necessary to assume `L+g<rho`. Thermal derivative corrections will be
summed rather than replaced by their worst value at every step.

Assume the line target approaches the dressed reference and the reference
thermal trajectory obeys

    ||x_h^0-y_h^0||<=A epsilon s_h^theta,
    ||y_h^T-y_h^0||<=B t_h^p.                           (3)

The second bound follows from a one-trajectory lemma whenever the value
bound `||G_h,T||<=g0 t_h^p` holds and the initial thermal deviation is at
most `B t_0^p`: take `B>=g0/(gamma^p-L)`. Indeed the relevant difference
satisfies `b_(h-1)<=L b_h+g0 t_h^p`. The derivative restriction in (2)
is stronger than `L<gamma^p`, so this denominator is positive.

The initial four-trajectory difference must be marked as well:

    delta_h=(x_h^T-x_h^0)-(y_h^T-y_h^0),
    ||delta_0||<=B0 epsilon T^p.                        (4)

Using exactly the same temperature-independent bare input for each
model's thermal and line flows gives `delta_0=0`. The target and dressed
reference may have different bare actions; no counterterm is retuned as
T changes. If a thermal preparation step is included before h=0, (4) is
its additional marked estimate. Separate unmarked bounds on the two
initial deviations do not imply (4).

## 2. The mixed derivative preserves the reference defect

**[proved]** Under (1)-(4),

    ||delta_h||<=D epsilon s_h^theta t_h^p,
                                     H<=h<=0,           (5)

where one sufficient constant is

    F=A(MB+g),
    D=max(B0,F/(rho-L))
                exp[g/(rho(1-gamma^(-p)))].             (6)

Proof. At one scale abbreviate

    x=x_h^0, y=y_h^0, a=x_h^T-x_h^0, b=y_h^T-y_h^0.

The line contribution splits exactly as

    F_h(x+a)-F_h(x)-F_h(y+b)+F_h(y)
      = F_h(x+a)-F_h(x+b)
          +F_h(x+b)-F_h(x)-F_h(y+b)+F_h(y).             (7)

The first term has norm at most `L||a-b||`. Twice applying the fundamental
theorem of calculus expresses the second as

    integral_0^1 du integral_0^1 dv
        D^2 F_h(y+u(x-y)+v b)[x-y,b],                 (8)

so its norm is at most `M||x-y|| ||b||`. The thermal correction obeys

    ||G_h,T(x+a)-G_h,T(y+b)||
         <=g t_h^p (||x-y||+||a-b||).                  (9)

Equations (3),(7)-(9) give

    ||delta_(h-1)||
       <=(L+g t_h^p)||delta_h||
                    +F epsilon s_h^theta t_h^p.       (10)

Normalize by `epsilon s_h^theta t_h^p`, whose value at h-1 is rho
times its value at h. Let `q_h` be the resulting nonnegative ratio.
If epsilon=0, the marked initial condition and (10) give delta=0
directly, so assume epsilon>0. Then

    q_(h-1)<=(L/rho+g t_h^p/rho)q_h+F/rho.            (11)

Let `Q_h` be the maximum of all previously encountered q values and
`F/(rho-L)`. Since `(L/rho)Q_h+F/rho<=Q_h`, (11) implies

    Q_(h-1)<=(1+g t_h^p/rho)Q_h.

The thermal corrections are summable:

    sum_(h=H+1)^0 t_h^p <=1/(1-gamma^(-p)).             (12)

Use `1+z<=exp(z)` in the finite product and `Q_0<=max(B0,F/(rho-L))`.
This proves (5)-(6), uniformly as H tends to minus infinity.

Only the reference thermal bound in (3) was needed for the mixed term.
The target bound then follows from (5):
`||x_h^T-x_h^0||<=(B+D epsilon)t_h^p` because `s_h<=1`.

## 3. A weighted-memory version

The same argument allows a genuinely history-dependent step without
pretending that its entire memory is one small running coupling.
Let the input to the step h be the finite history

    X_h=(x_h,x_(h+1),...,x_0),

and similarly for Y. Let derivatives in input slot j refer to scale
h+j. Suppose nonnegative coefficient arrays, independent of T,h, bound

    ||D_j F_h||<=a_j,
    ||D_j D_k F_h||<=m_jk,
    ||D_j G_h,T||<=g_j t_h^p.                           (13)

Only slots `0<=j,k<=-h` are present. Assume the following weighted sums
are finite:

    L_eff=sum_(j>=0) a_j gamma^[-(p-theta)j] <rho,
    G_delta=sum_(j>=0) g_j gamma^[-(p-theta)j],
    G_line=sum_(j>=0) g_j gamma^(theta j),
    M_mix=sum_(j,k>=0) m_jk gamma^(theta j-p k).        (14)

Use (3)-(4) componentwise on the trajectories. Then (5) holds with

    F_eff=A(B M_mix+G_line),
    D=max(B0,F_eff/(rho-L_eff))
              exp[G_delta/(rho(1-gamma^(-p)))].        (15)

Proof. Apply the exact difference decomposition (7) to the history
vectors. Its first difference is bounded by
`sum_j a_j ||delta_(h+j)||`. The Hessian term is bounded by

    sum_jk m_jk ||x_(h+j)^0-y_(h+j)^0||
                       ||y_(h+k)^T-y_(h+k)^0||
      <=A B epsilon s_h^theta t_h^p M_mix.             (16)

The marked thermal-map difference contributes the two sums

    t_h^p sum_j g_j ||delta_(h+j)||,
    A epsilon s_h^theta t_h^p G_line.

The normalized past ratios carry the factor
`gamma^[-(p-theta)j]`. Taking their running maximum gives exactly the
last induction of Section 2 with `L_eff,G_delta,F_eff`. This proves (15).
There is no factor proportional to the number of past scales.

Any fixed finite memory length is covered. A sufficient infinite-memory
condition is exponential decay of the derivative arrays, for example
`a_j<=C_a r^j`, `g_j<=C_g r^j`, and `m_jk<=C_m r^(j+k)` with
`r gamma^theta<1`, together with the numerical restriction on L_eff.
The weighted conditions (14), not that particular exponential example,
are the actual hypotheses. A stronger thermal mark at an older scale
can be included directly by improving the corresponding g_j.

The one-trajectory input B can also be proved in this setting: if
`||G_h,T||<=g0 t_h^p`, put

    L_p=sum_j a_j gamma^(-pj)<=L_eff<gamma^p,
    B>=max(B_initial,g0/(gamma^p-L_p)).                 (17)

The same causal induction gives both unmarked thermal deviations,
provided the trajectories remain in the domain. Thus the additional
four-trajectory information is precisely the mixed second derivative
and the marked derivative of G, not another temperature-dependent
choice of the relevant counterterm.

## 4. Scalar increments and normalization factors

Let the exact scalar density produced at scale h be

    f_h^0(z)=s_h^2 Phi_h(z),
    f_h^T(z)=s_h^2[Phi_h(z)+Psi_h,T(z)].                (18)

This scalar includes every Gaussian prefactor and field-rescaling
Jacobian assigned to that step. If some factor is kept separately,
it must satisfy the same comparison as a separate scalar map and then
be added back. A normalized correlation functional is not a substitute
for (18).

For the memoryless case assume

    ||D Phi_h||<=c1,       ||D^2 Phi_h||<=c2,
    ||D Psi_h,T||<=cG t_h^p.                           (19)

Define the target-minus-reference increments by

    r_h(0)=f_h^0(x_h^0)-f_h^0(y_h^0),
    r_h(T)=f_h^T(x_h^T)-f_h^T(y_h^T).                  (20)

**[proved]** Equations (3),(5),(19) imply

    |r_h(0)|<=c1 A epsilon s_h^(2+theta),
    |r_h(T)-r_h(0)|
       <=C epsilon s_h^(2+theta)(T/s_h)^p,             (21)

where `C=(c1+cG)D+c2 A B+cG A` suffices. Apply (7)-(9) with Phi and
Psi in place of F and G. The term involving `delta_h` has the required
gain by (5); the mixed Hessian term contains the product of the line
defect and reference thermal deviation; and the derivative of Psi
marks its thermal difference with the line defect. Since `t_h<=1`,
the extra factor `t_h^p` on Psi's delta term can be discarded.

For scalar maps with memory, sufficient bounds are

    ||D_j Phi_h||<=c_j,
    ||D_j D_k Phi_h||<=d_jk,
    ||D_j Psi_h,T||<=e_j t_h^p,                         (22)

with finite sums

    C_line=sum_j c_j gamma^(theta j),
    C_delta=sum_j c_j gamma^[-(p-theta)j],
    D_mix=sum_jk d_jk gamma^(theta j-p k),
    E_line=sum_j e_j gamma^(theta j),
    E_delta=sum_j e_j gamma^[-(p-theta)j].              (23)

Then replace the first constant in (21) by `A C_line`, and the second by

    D(C_delta+E_delta)+A B D_mix+A E_line.              (24)

The proof uses the same componentwise estimates as (16). It is valid
for source-independent normalization maps when their complete derivatives
obey (22)-(23); no cancellation of such factors is assumed automatically.

## 5. Consequence and exact remaining hypotheses

For p>2, (21) supplies the active-scale estimates in
`thermal/THERMAL_CRITERION.md` with `omega(s)=epsilon s^theta`. If the
exact scalar expansion includes the ultraviolet preparation and final
integration as required there, its theorem gives an `o(T^2)` target-minus-
dressed-reference free-energy difference. When `p>2+theta`, the summed
scale contribution is `O(epsilon T^(2+theta))`. A complete shell
partition can use the zero-terminal-remainder convention proved in
`thermal/WINDING_ESTIMATE.md`; otherwise its terminal remainder still
needs its own bound.

The present result is a closure of nonlinear causal estimates. It does
not prove the following microscopic inputs:

- A common normalized Banach space and a uniform reserve for the
  arbitrary-degree, retarded effective kernels at every RG scale.
- The target's line estimate (3) relative to its correctly dressed
  reference, including the fixed relevant counterterm.
- The exact map derivatives (2), or the weighted history derivatives
  (13)-(14), through localization, covariance changes and anomalous
  wave-function normalization.
- The complete scalar derivative estimates (19), or (22)-(23), with
  every source-independent normalization factor included.

The fixed-history determinant winding theorem motivates the marked G
and Psi estimates, but does not establish this whole list. Individual
thermal histories of size `O((T/s_h)^p)` alone lose the factor
`epsilon s_h^theta` when subtracted. Uniform second derivatives and a
marked thermal derivative are the precise extra information that keeps
it in the four-trajectory comparison.

## 6. Finite controls

`python3 contrib/thermal_ensemble_20261028/ensemble/history_checks.py`
uses only the standard library. It checks 21 nonlinear examples with
240-digit decimal arithmetic, at depths through 120 scales. The cases
include two-slot memory, the mixed second-derivative term, scalar-map
comparisons, and thermal derivative coefficients for which `L+g>rho`.
All bounds pass; results are in `history_checks.json`. These are controls
of the abstract estimates, not tests of the microscopic GN hypotheses.
