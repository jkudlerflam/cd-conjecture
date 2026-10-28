# What the scalar RG calculation must prove

Status: **[proved abstract thermal comparison and history lemmas;
microscopic application not established]**. These results isolate a sufficient
finite-temperature estimate. They do not assert that the retarded action of
C44 or the vacancy chain satisfies that estimate.

Provenance: OpenAI GPT-6-family Codex thermal worker, 2026-10-28. AI-assisted;
no human review. The estimates below are derived here.

## 1. The reference must include the marginal interactions

The comparison model is the multichannel Luttinger model with the **limiting
forward couplings and limiting velocities** of the microscopic flow. It is
not generally the bare free-fermion Hamiltonian. A forward interaction can
change the sound velocities at first order and therefore change the
coefficient of the free energy at order T^2 without changing the number of
bosons.

Write f(T) for free energy per cell after the thermodynamic limit, e=f(0),
and use a star for the dressed reference. The desired scalar comparison is

    [f(T)-e]-[f_*(T)-e_*] = o(T^2).                     (1)

Every field-rescaling Jacobian and every source-independent factor is included
in f and f_*. Normalized density correlations omit these factors. Equation
(1) cannot be obtained merely by setting external fields to zero in a
normalized correlation comparison.

## 2. A thermal modulus theorem

Fix gamma>1, let s_h=gamma^h for integers h<=0, and take 0<T<1.
Put H=ceil(log_gamma T), so T<=s_H<gamma T. Suppose exact convergent
scalar expansions give

    f(T)-f_*(T)=a(T)+sum_(h=H)^0 r_h(T)+q(T),
    e-e_*=a(0)+sum_(h=-infinity)^0 r_h(0).               (2)

The remainder q includes the complete last integration, including a possible
shift in the convention for the terminal scale. It is a **difference from
the reference**. These expansions are hypotheses, not a definition that
discards unknown normalization factors.

Let omega:(0,1]->[0,infinity) be bounded, with omega(s)->0 as s->0. Suppose,
for constants C and p>2 independent of T and h,

    a(T)-a(0)=o(T^2),
    |r_h(0)| <= C s_h^2 omega(s_h),
    |r_h(T)-r_h(0)| <= C s_h^2 omega(s_h) (T/s_h)^p,
                            H<=h<=0,
    |q(T)| <= C T^2 sup_(0<s<=gamma T) omega(s).        (3)

Only values s<=1 are used; for gamma T>1 replace that upper endpoint by one.

**[proved]** Equations (2)-(3) imply (1). No power-law rate of approach
to the fixed point is required.

Proof. Subtract the two expansions. The missing low scales obey

    sum_(h<H) |r_h(0)|
       <= C s_H^2/(gamma^2-1)
                           sup_(0<s<=s_H) omega(s)
       = o(T^2).                                      (4)

The active-scale difference divided by T^2 is bounded by

    C sum_(j=0)^(-H) (T/s_H)^(p-2)
                       gamma^[-(p-2)j] omega(s_H gamma^j).
                                                        (5)

The prefactor is at most one. Extend each sum to j>=0 by setting terms
beyond -H to zero. For each fixed j its omega factor tends to zero, while
the boundedness of omega gives a summable geometric majorant. Dominated
convergence for this countable sum makes (5) tend to zero. The a and q
terms obey the required estimate directly. This proves (1).

This criterion permits logarithmic convergence, such as
omega(s)=1/(1+|log s|), which is relevant to marginally irrelevant flows.
The theorem does not establish that any particular marginal flow has that
rate or the required finite-temperature comparison.

**[proved rate]** If omega(s)<=A s^theta with theta>0 and p>2+theta, then
the non-ultraviolet error in (1) is at most C' A T^(2+theta). Indeed (5)
is bounded by

    C A T^(p-2) s_H^(2+theta-p)
                      /[1-gamma^(-(p-2-theta))]
       <= C A T^theta/[1-gamma^(-(p-2-theta))].          (6)

Equation (4) and the terminal bound have the same rate. At p=2+theta
the active sum acquires log(1/T); if 2<p<2+theta it is O(T^p). All three
cases still give o(T^2). The ultraviolet term must separately have the
claimed stronger rate if that rate is desired for the full answer.

## 3. Why a zero-temperature beta-function bound is insufficient

The zero-temperature estimates alone do not control (5). For example, set
r_h(0)=0 for all h and, on the thermal terminal scale, set

    r_H(T)=a T^2,       r_h(T)=0 for h!=H.

All zero-temperature remainders and all zero-temperature running-coupling
differences vanish, but the thermal coefficient shifts by a. This is a
counterexample to that **inference between scale estimates**, not a
microscopic Hamiltonian construction. The thermal difference and final
integration in (3) are genuine additional obligations.

A second failure is comparing with the bare free reference. For a compact
boson with speed v+delta v, the thermal free energy is
-pi T^2/[6(v+delta v)] rather than -pi T^2/(6v). Its difference is of
order delta v T^2 at fixed delta v. No small fixed norm of delta v makes
that difference o(T^2). The dressed comparison removes this error.

## 4. Finite-temperature histories can use the zero-temperature counterterm

The relevant counterterm must be selected once for the physical Hamiltonian.
One must not retune it for each T. The following elementary lemma shows
how a thermal forcing estimate can control a relevant direction with that
fixed initial datum.

Let x_h(T) and x_h(0) be two causal RG states for H<=h<=0, with
d_h=|x_h(T)-x_h(0)|. The state may include normalized relevant terms;
the Lipschitz amplification need not be less than one. Assume

    d_0 <= B A T^p,
    d_(h-1) <= L d_h + B A T^p s_h^(theta-p),
                                      H<h<=0,          (7)

where L>=0, A,B>0, theta>=0 and p>theta with

    L < gamma^(p-theta).                               (8)

**[proved]** There is a T-independent D such that

    d_h <= D A s_h^theta (T/s_h)^p,
                                             H<=h<=0. (9)

Choose D>=B and D>=B/[gamma^(p-theta)-L]. The estimate holds at zero.
If it holds at h, (7) yields

    d_(h-1) <= (LD+B) A T^p s_h^(theta-p)
              <= D A T^p s_(h-1)^(theta-p),

which proves the induction. In particular d_H<=D A (gamma T)^theta.
For the usual linear relevant amplification gamma, a finite-temperature
forcing with p>1+theta accommodates sufficiently small nonlinear losses.

This lemma is applicable to a complete RG state only after its causal norm
and Lipschitz constant have been established. It does not replace the
nonlocal history estimates in COUNTERTERM_INDUCTION.md, nor prove the
finite-temperature forcing in (7). A finite window of memory can be
included as part of x; an infinite history requires its own weighted norm.

## 5. A direct single-scale periodization estimate

There is a concrete source for the factor (T/s_h)^p at fixed running data.
Let a translation-invariant, infinite-time, single-scale fermion covariance
satisfy for every integer N>=2

    |g_h(t,x)| <= C_N s_h
                       [1+s_h(|t|+|x|)]^(-N).          (10)

Suppose the finite-temperature covariance is its antiperiodization,

    g_(h,beta)(t,x)=sum_(m in Z)(-1)^m
                                      g_h(t+m beta,x).

For |t|<=beta/2 and beta s_h>=1, absolute convergence and
|t+m beta|>=(|m|-1/2)beta give

    |g_(h,beta)(t,x)-g_h(t,x)|
       <= C'_N s_h (beta s_h)^(-N).                    (11)

This is an elementary consequence of summing
sum_(m!=0)(|m|-1/2)^(-N). It is uniform in x. A spatially integrable
version follows by starting from (10) with a further decay power.

Equation (11) by itself is **not** a bound on a connected vacuum diagram:
relative times integrated over a cylinder may cross its seam, and products
can wind. A complete proof must lift each connected cluster, classify
nontrivial windings or long paths, preserve determinant bounds rather than
sum factorially many Wick pairings, and sum the resulting size factors.
It must also compare the T-dependent running data, as in Section 4.
These are the precise steps left when using (11) in (3).

Exponential time memory of the initial retarded kernels gives the analogous
image bound before RG. It does not assert exponential decay for the
massless, already-renormalized kernels at every scale.

## 6. Source audit and actual remaining application

**[read]** Giuliani and Mastropietro,
[arXiv:1211.3030, Section 3.4, (3.83)-(3.91)](https://arxiv.org/pdf/1211.3030),
retain scalar scale contributions and compare finite- and infinite-volume
histories to establish a finite-size free-energy coefficient for a
non-solvable Ising model. The quantum multichannel application here is not
their theorem. Its limiting forward interaction is generally nonzero, so
one must compare to an interacting Luttinger reference.

**[read]** Mastropietro and Porta,
[arXiv:2108.13708v3, Section 8.2 and Theorem 8.5](https://arxiv.org/html/2108.13708v3),
give a convergent free-energy construction with bounded running couplings
and establish a reference beta-function cancellation. Analyticity of free
energy in coupling does not identify its T^2 coefficient.

**[read]** Porta, Scola and Singh,
[arXiv:2509.08665v1, Theorem 2.6 and Section 6](https://arxiv.org/html/2509.08665v1),
establish density/current response for their specified lattice model. That
result does not supply (2)-(3) for the encoded chain's arbitrary-degree
retarded endpoints.

For the current repository the new target is explicit: build the scalar
**difference** expansion (2), with its terminal factor and all normalization
Jacobians; prove its periodization and history estimates (3); and evaluate
the dressed reference pressure. The theorem above then converts those
estimates to o(T^2), even if approach to the fixed point is only logarithmic.
The normalized marginal-flow induction is valuable input but does not
already contain these statements.

Even after (1), a coefficient

    -(pi/12) sum_a [1/v_(a,R)+1/v_(a,L)]

is velocity weighted. Identifying it with N compact bosons requires the
reference current/stress construction or equivalent mode identification.
C46 supplies the Gaussian oscillator and Sugawara calculations. None of
the lemmas here proves a microscopic c=N or refutes the c-d conjecture.
