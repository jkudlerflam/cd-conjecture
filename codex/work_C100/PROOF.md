# C100. Thermal and transport routes

2026-10-29. All logarithms are natural unless a base is displayed. Lattice
spacing, Planck's constant, and Boltzmann's constant are one.

## 1. Scope, conventions, and the resulting inequalities

[Assumed setting.] The Hamiltonian is a one-site translation-invariant
nearest-neighbor chain on local dimension d. Its infrared theory consists of
finitely many unitary nonchiral CFT sectors, with central charges c_a and
positive velocities v_a. Write

    A = sum_a c_a v_a,
    B = sum_a c_a/v_a,
    C3 = sum_a c_a v_a^3,
    c = sum_a c_a.

The proposed inequality is A <= (pi/log 2) M log d. No proof of that
inequality is obtained here, and no microscopic counterexample is claimed.

[Proved lattice theorem.] For any Gibbs state at T > 0, or a stationary
ground state at T = 0, define h_q = sum_j exp(-iqj) g_j using a centered
optimal bond representative ||g|| = M. For allowed nonzero ring momenta,
with |q| <= pi, put

    m1(q,T) = <[h_-q,[H,h_q]]>_T/(2L).

Then

    0 <= m1(q,T) <= 8 M^3 q^2.                         (1)

The same bound passes to thermodynamic limits when the correlators exist.
It bounds a positive spectral moment. It does not assume current
conservation or a CFT.

[Proved CFT calculation.] For the continuum stress-tensor density,

    m1_CFT(q,T) = (pi/6) A T^2 q^2 + C3 q^4/(24 pi).  (2)

[Assumed quantitative matching; proved implication.] Suppose at T_* > 0
and some nonzero q_* a positive part of the lattice spectral moment is at
least alpha times (2), where 0 < alpha <= 1. Then

    A + C3 q_*^2/(4 pi^2 T_*^2)
        <= 48 M^3/(pi alpha T_*^2),                  (3)

and in particular A <= 48 M^3/(pi alpha T_*^2). A corresponding lower
bound on the full moment suffices formally, but matching a specified
low-energy spectral window is the physically useful hypothesis. One must
not infer that hypothesis from a large ultraviolet background.

[Assumed thermal matching; proved implication, inherited entropy input
from C69.] If at some T_s > 0,

    s(T_s) >= alpha_s (pi/3) B T_s,

then

    B <= 3 log d/(pi alpha_s T_s),
    A <= 3 v_max^2 log d/(pi alpha_s T_s).            (4)

The second line follows from A <= v_max^2 B. Under C69's local real-time
identification of CFT speeds with propagation cones, its conservative
Lieb-Robinson estimate v_max <= 6 e M gives

    A <= 108 e^2 M^2 log d/(pi alpha_s T_s).          (5)

The Lieb-Robinson theorem and its identification assumption are inherited
from C69, not rederived here. Equations (3) and (4) are complementary.
Their minimum is the strongest of these two conclusions when both matching
hypotheses hold. They have different hypotheses and neither uniformly
dominates the other. Their constants are explicit, not claimed optimal.

[Proved limitation of these estimates.] As T_* or T_s tends to zero,
these upper bounds diverge. An infrared asymptotic formula supplies
sufficiently small matching temperatures, not a uniform positive lower
bound on those temperatures in units of M. No choice such as T_* = M is
justified by the hypotheses.

## 2. Gauge and elementary microscopic bounds

[Proved.] The definition of M is equivalent to minimizing ||g|| over

    g = h + alpha I + o tensor I - I tensor o.

For each choice of o, the minimizing scalar centers the spectral interval.
The permitted additions form a finite-dimensional closed linear space.
A minimizing sequence of the matrices g has a bounded subsequence, whose
limit remains in the affine space. Thus a minimizing g exists, with

    -M I <= g <= M I.

Use this representative throughout. On a ring the coboundary telescopes;
the scalar shifts the energy origin. Neither changes the dynamics or the
CFT velocities.

[Proved.] With g_i acting on sites i,i+1, define

    j_i = i[g_(i-1),g_i],       J = sum_i j_i.

Only neighboring bonds fail to commute. Consequently

    d g_i/dt = j_i - j_(i+1),
    ||j_i|| <= 2 M^2,
    ||J|| <= 2 L M^2.                                (6)

The variance of a Hermitian operator is bounded by a quarter of its
spectral width squared. Hence

    Var_T(j_i) <= 4 M^4,
    Var_T(J)/L <= 4 L M^4.                            (7)

The second bound grows with L. Replacing it by 4 M^4 is an incorrect
inference. Indeed Var(J)/L is a spatial sum of connected correlations,
not a single-site variance.

[Proved.] At infinite temperature, trace j_i = 0. Disjoint supports have
factorized traces. There are at most five overlapping three-site supports
for each j_i, and each covariance has magnitude at most 4 M^4. Thus

    <J^2>_(infinity)/L <= 20 M^4.                      (8)

This conservative estimate needs no optimized overlap count. At finite
temperature the product-state factorization is lost. A correlation-length
estimate could replace the factor L in (7), but its temperature dependence
would have to be retained.

[Proved.] Imaginary-time separation suppresses high frequencies without
requiring a contact prescription. At the thermal midpoint,

    C_i = Tr(rho^(1/2) j_i rho^(1/2) j_i) - <j_i>^2
        <= Var_T(j_i) <= 4 M^4,
    C_i >= 0.                                        (9)

For the upper bound, in the energy basis use
sqrt(p_n p_m) <= (p_n+p_m)/2. For positivity, off-diagonal terms are
nonnegative and the diagonal connected terms form an ordinary variance.

[Proved gauge qualification.] The local energy density and regular
finite-frequency current spectra are gauge dependent. A density change
g_i -> g_i + o_i - o_(i+1) changes the continuity current by d o_i/dt,
and the total current by d(sum_i o_i)/dt. In an energy basis this addition
has matrix elements i(E_n-E_m) O_nm, so its equal-energy matrix elements
vanish. In particular the finite-volume exactly zero-frequency Drude
projection is unchanged. General finite-frequency correlators are not
unchanged. Matching the optimized representative to the continuum stress
tensor must therefore be specified; a ground-state Casimir fit alone
does not establish such matching.

## 3. CFT thermal quantities and the current correlator

[Proved within the assumed CFT.] For one nonchiral sector,

    f(T) = e0 - pi c T^2/(6v),
    u(T)-e0 = pi c T^2/(6v),
    s(T) = pi c T/(3v),
    C_V/L = pi c T/(3v).                              (10)

The two chiral energy densities each equal pi c T^2/(12v). If H_R and
H_L are their integrated energies, the continuum total energy current is

    J = v(H_R-H_L).

Chiral factorization and thermal differentiation give

    Var(H_R)/L = -d[<H_R>/L]/d beta
               = pi c T^3/(6v),
    Var(H_L)/L = pi c T^3/(6v),
    Cov(H_R,H_L)/L = 0.

Therefore

    Var(H)/L = (pi/3) B T^3,
    Var(J)/L = (pi/3) A T^3.                          (11)

These are continuum identities. On the lattice the equal-time current
variance can also contain ultraviolet weight, including ground-state
fluctuations if J is not conserved.

[Proved CFT correlator.] Normalize the chiral stress tensor by
<T(z)T(0)>_conn = c/(2 z^4) in the plane, and set h_R = v T/(2pi).
The thermal cylinder map then gives

    <h_R(x,t) h_R(0,0)>_conn
      = c v^2/(8pi^2) (pi T/v)^4
        / sinh^4[(pi T/v)(x-vt+i0)],                  (12)

with the corresponding opposite-chirality prescription for h_L. At
imaginary-time separation beta/2 the current correlator is regular,

    <j(x,-i beta/2)j(0,0)>_conn
      = sum_a [c_a pi^2 T^4/4]
        sech^4(pi T x/v_a).                          (13)

Using integral sech^4 y dy = 4/3 over the real line,

    integral dx <j(x,-i beta/2)j(0,0)>_conn
       = (pi/3) A T^3.                               (14)

At one spatial point the answer is instead

    <j(0,-i beta/2)j(0,0)>_conn = (pi^2/4) c T^4.    (15)

Thus the local norm in (9) naturally measures c, with no velocity
weight, whereas the spatial integral measures A and requires spatial
correlation control. If (15) has a lattice lower matching factor alpha
at T_*, (9) gives c <= 16 M^4/(pi^2 alpha T_*^4).
This is another crossover inequality, not a dimension-only bound on c.

[Proved local spectral sum rule.] With the same delta-function convention
as Section 4, the positive zero-temperature local current spectrum has
total weight <j_i^2>-<j_i>^2 <= 4 M^4 on the lattice. In the continuum,
the zero-temperature limit of (12) at x=0 and imaginary time tau>0 is

    <j(0,-i tau)j(0,0)>_conn = c/(4pi^2 tau^4).

Since integral_0^infinity omega^3 exp(-omega tau) d omega = 6/tau^4,
its local spectral density is c omega^3/(24pi^2), independent of the
velocities. Matching this density from zero to a specified cutoff
Omega_* with lower factor alpha gives

    alpha c Omega_*^4/(96pi^2) <= 4 M^4.

The cutoff is again indispensable. A local energy-density spectrum
instead has coefficient sum_a c_a/v_a^2 in place of c and total
connected lattice weight at most M^2. Neither local zeroth-moment sum
rule directly measures A.

[Proved CFT transport coefficient; assumed transport convention.] Use

    Re kappa(omega) = pi D_th(T) delta(omega)
                     + Re kappa_reg(omega),

with the Green-Kubo prefactor beta times the imaginary-time integral.
For a conserved continuum J the imaginary-time integral is
beta Var(J); the real part of the one-sided time integral is pi delta.
It follows that

    D_th_CFT = beta^2 Var(J)/L = (pi/3) A T.          (16)

In the convention Re kappa = 2pi D delta + ..., the quantity called D
is pi A T/6. If the delta-function coefficient itself is called the
Drude weight, it is pi^2 A T/3. These factors are conventions, not
different predictions. The convention in (16) agrees with the coefficient
reported in Gawedzki, Langmann, and Moosavi's CFT calculation.

[Proved logical qualification.] The fixed-point Hamiltonian conserves J.
A lattice Hamiltonian need not. Its infrared CFT by itself gives no
interchange theorem between T -> 0 and the infinite-time transport limit.
If an irrelevant interaction produces a nonzero relaxation rate at every
T > 0, a narrow peak can replace the delta function. This is a possible
mechanism, not a claim that every nonintegrable chain has zero Drude
weight. Positive integrated spectral weight, used below, avoids assuming
an exact lattice delta function.

[Proved CFT comparison.] For perfectly transmitting thermal reservoirs,
the steady current is

    j_NESS = (pi/12) sum_a c_a (T_L^2-T_R^2).

The velocity cancels between each chiral density and its speed. The
linear thermal conductance is pi c T/6. These quantities measure c,
not A. The lattice norm gives |j_NESS| <= 2 M^2, which again produces
only a temperature-dependent comparison if a CFT NESS is matched.

## 4. Spectral normalizations and continuum moments

[Definition.] Let the unsymmetrized dynamical structure factor have no
2pi in its delta-function normalization,

    S_h(q,omega) = (1/L) sum_(n,m) p_n
       |<m|h_q|n>|^2 delta(omega-E_m+E_n).

For nonzero q the disconnected uniform expectation does not contribute.
To include parity-breaking chains, use the momentum average

    Sbar_h(q,omega) = [S_h(q,omega)+S_h(-q,omega)]/2.

This is positive and satisfies detailed balance
Sbar_h(q,-omega) = exp(-beta omega) Sbar_h(q,omega). Expanding the
double commutator in the energy basis proves

    m1(q,T) = integral_all d omega omega Sbar_h(q,omega)
            = integral_0^infinity d omega
                omega(1-exp(-beta omega)) Sbar_h(q,omega).             (17)

Every term in the second expression is nonnegative. At zero temperature
use the positive absorption spectrum directly. A low-frequency window
can be discarded or retained without changing this ordering.

[Proved CFT spectrum, equal left and right speeds.] For q > 0 in one
sector, define

    K(q,T) = c v^2/(24pi) q [q^2+(2pi T/v)^2].

Then the continuum stress density has

    S_h(q,omega) = K(q,T) {
       delta(omega-vq)/(1-exp(-beta vq))
       + delta(omega+vq)/(exp(beta vq)-1) } .          (18)

One derivation, fixing all normalization factors, starts with the
midpoint version of (12). The Fourier identity

    integral dy exp(-iky) sech^4 y
       = pi k(k^2+4)/[6 sinh(pi k/2)]

follows by t=exp(2y): the integral is
8 integral_0^infinity t^(1-ik/2)/(1+t)^4 dt
= 8 B(2-ik/2,2+ik/2) = (4/3)|Gamma(2+ik/2)|^2.
Using |Gamma(1+ix)|^2 = pi x/sinh(pi x) gives the displayed result.
This yields the density midpoint Fourier transform K/sinh(beta vq/2).
Chiral propagation restricts the frequencies to +/- vq. Detailed balance
fixes the ratio of the two weights. These two conditions determine (18).
The current midpoint Fourier transform has an additional v^2, since
the current is +/- v times the chiral energy density.

[Proved.] The equal-time structure factor, static imaginary-time
susceptibility, and first moment are consequently

    S_h(q) = K(q,T) coth(beta vq/2),
    chi_hh(q,T) = integral_0^beta d tau
                  <h_-q(-i tau)h_q>/L
                = 2K(q,T)/(vq)
                = pi c T^2/(3v) + c v q^2/(12pi),
    m1(q,T) = vq K(q,T)
            = pi c v T^2 q^2/6 + c v^3 q^4/(24pi).  (19)

Add sectors to obtain (2), as well as

    chi_hh_CFT(q,T) = pi B T^2/3 + A q^2/(12pi).     (20)

At T = 0, S_h(q) = sum c_a v_a^2 |q|^3/(24pi).
These are the requested c and v weights. The current spectrum for each
sector is v_a^2 times the corresponding density spectrum; its zero-T
equal-time coefficient is sum c_a v_a^4 |q|^3/(24pi).

[Proved limitation.] Equation (20) is a promising quantity with A in a
q^2 coefficient and no explicit T. Its positive spectral representation
uses an inverse frequency,

    chi_hh(q,T) = 2 integral_0^infinity d omega
            (1-exp(-beta omega)) Sbar_h(q,omega)/omega.               (21)

The norm and double-commutator estimates control zeroth or positive
moments, not this inverse moment. A small positive spectral weight can
be moved to arbitrarily small frequency while respecting those bounds
and making (21) large. An assumed gap would give a comparison, but the
chain is gapless. A generic lattice susceptibility also has ultraviolet
analytic terms at order q^2, and density improvements affect them.
No upper bound on the coefficient of q^2 follows from (1).

## 5. Proof of the lattice spectral bound and its transport form

[Proved.] Set C_i = [g_i,g_(i+1)]. Direct expansion gives

    [H,h_q] = (exp(-iq)-1) sum_i exp(-iqi) C_i.       (22)

For a stationary state, <[H,C_i]> = 0. Only
j = i-1,i,i+1,i+2 have [g_j,C_i] nonzero. Hence

    <[h_-q,[H,h_q]]>
      = (exp(-iq)-1) sum_i exp(-iqi)
          sum_(j=i-1)^(i+2)
          [exp(iqj)-exp(iq(i+1/2))] <[g_j,C_i]> .    (23)

For a ring, use unwrapped local indices in this expression. Allowed
momenta make exp(iqj) periodic, so the subtraction is valid across the
seam. Rings larger than the displayed supports suffice; the result is
used in the thermodynamic limit.

The relevant four distances from i+1/2 are 3/2,1/2,1/2,3/2. Their sum
is 4. Also

    ||C_i|| <= 2 M^2,
    ||[g_j,C_i]|| <= 4 M^3,
    |exp(-iq)-1| <= |q|,
    |exp(iqj)-exp(iq(i+1/2))| <= |q| |j-i-1/2|.

Taking absolute values in (23) therefore gives

    |<[h_-q,[H,h_q]]>|/L <= 16 M^3 q^2.

Positivity from (17) supplies the lower bound, proving (1). Keeping the
first difference exact also yields m1 <= 8 M^3 |q| |2sin(q/2)|, a
slightly sharper finite-q version. Neither changes the small-q constant.

[Proved conditional theorem.] If W(q_*,T_*) is the moment (17) restricted
to a positive-frequency window and

    W(q_*,T_*) >= alpha [pi A T_*^2 q_*^2/6
                        + C3 q_*^4/(24pi)],

then W <= m1 <= 8 M^3 q_*^2. Dividing by
alpha pi T_*^2 q_*^2/6 proves (3). This formulation tolerates peak
broadening. It requires quantitative matching of the integrated weight
in the chosen window, rather than the unsupported assertion that the
lattice has the CFT delta functions at every positive temperature.

[Proved transport sum rule, with limits specified.] Let
k_q = 2sin(q/2). Continuity fixes, for nonzero frequency,

    S_j(q,omega) = omega^2 S_h(q,omega)/k_q^2.

The dissipative energy-transport kernel at nonzero q, momentum averaged
as above, has the Lehmann form

    Re kappa(q,omega)
      = [pi/T] [(1-exp(-beta omega))/omega]
           Sbar_j(q,omega)
      = [pi/(T k_q^2)] omega(1-exp(-beta omega))
           Sbar_h(q,omega).                          (24)

The expression is nonnegative. Integrating over all frequencies gives

    integral_all Re kappa(q,omega) d omega
       = 2pi m1(q,T)/(T k_q^2).                      (25)

At finite q a strictly static conserved-current component is not
determined by density continuity; (24) concerns the longitudinal
dissipative kernel. For the homogeneous transport sum rule take the
thermodynamic limit first and then the q -> 0 limit of this kernel,
retaining any peak that collapses to zero frequency. If this limit and
its integrated weight exist, define

    Theta(T) = lim_(q->0) 2m1(q,T)/k_q^2.

Then

    0 <= Theta(T) <= 16 M^3,
    integral_all Re kappa(omega) d omega = pi Theta(T)/T,
    0 <= D_th(T) <= Theta(T)/T <= 16 M^3/T.          (26)

If weight is lost in the limiting procedure, positivity still gives the
corresponding upper inequality rather than an asserted equality. The
finite-q theorem (1) requires none of these transport-limit assumptions.
This is the energy-density version of a thermal conductivity sum rule.

[Proved CFT check.] Equation (2) gives

    Theta_CFT(T) = pi A T^2/3,
    integral_all Re kappa_CFT = pi^2 A T/3.

It agrees with (16). If an actual lattice Drude coefficient obeys
D_th(T_*) >= alpha pi A T_*/3, (26) yields the same weaker form of (3).
That lower bound is an extra assumption. A ground-state CFT does not
supply it.

[Proved obstruction to a zero-temperature subtraction.] At T=0 the
continuum term in m1 is order q^4, while the lattice bound is order
q^2. The missing q^2 weight may come from ultraviolet excitations.
Even if Theta(T)-Theta(0) had a CFT T^2 term, subtracting spectra at
different temperatures destroys positivity. The norm bound on Theta
does not bound [Theta(T)-Theta(0)]/T^2 uniformly.

## 6. Free energy, high temperature, convexity, and entropy

[Proved lattice thermodynamics.] The spectrum of H lies in [-LM,LM].
For p_L(beta) = L^(-1) log Tr exp(-beta H),

    log d - beta M <= p_L(beta) <= log d + beta M,
    p_L''(beta) = Var_beta(H)/L >= 0.

In a Gibbs state s_L = p_L-beta p_L', so

    0 <= s_L(T) <= log d,
    ds_L/dT = C_V/(LT) >= 0.                         (27)

In particular, finite entropy bounds the area under C_V/T, not its
initial slope. Assuming zero residual entropy,

    integral_0^infinity [C_V(T)/L] dT/T = log d.      (28)

With residual entropy, the right side is log d-s(0). For a CFT,
C_V/L = pi B T/3 near zero. Equation (28) controls the product of B
with the width of that low-temperature regime.

[Proved high-temperature coefficient bound.] At beta = 0 put
e_infinity = Tr(g)/d^2 and
sigma_infinity^2 = lim Var_0(H)/L. A bond has spectral width 2M,
so its variance in the tracial state is at most M^2. Only the same or
adjacent bonds have nonzero covariance at infinite temperature.
Cauchy-Schwarz bounds each such covariance by M^2, yielding

    0 <= sigma_infinity^2 <= 3 M^2,
    p(beta) = log d-beta e_infinity
              + beta^2 sigma_infinity^2/2 + O(beta^3),
    f(T) = -T log d+e_infinity-sigma_infinity^2/(2T)
              + O(T^(-2)),
    s(T) = log d-sigma_infinity^2/(2T^2)+O(T^(-3)).  (29)

These displayed coefficients can be obtained at finite volume before
taking the local cumulant limit. No radius of convergence connecting
high and low temperatures is assumed. In contrast the low-temperature
coefficient is B, as in (10), not A.

[Proved counterexample to a thermodynamic inference, not a lattice
construction.] For any S > 0 and tau > 0 consider the free-energy density

    f_tau(T) = -S sqrt(T^2+tau^2),
    p_tau(beta) = S sqrt(1+tau^2 beta^2).

It obeys

    s_tau(T) = S T/sqrt(T^2+tau^2),
    C_tau(T) = S T tau^2/(T^2+tau^2)^(3/2) >= 0,
    u_tau(T) = -S tau^2/sqrt(T^2+tau^2),
    e0 = -S tau,
    f_tau(T) = e0-S T^2/(2tau)+O(T^4),
    p_tau(beta) = S + S tau^2 beta^2/2 + O(beta^4).

It is real analytic for T > 0 and for finite real beta, has p'' > 0,
has entropy between zero and S, has bounded energy, and has the correct
even leading high-temperature entropy correction. With S = log d and
S tau <= M it also obeys the pressure and energy bounds above. Nevertheless
its coefficient identified with B would be 3S/(pi tau), arbitrarily
large as tau decreases.

This is not a partition function of a constructed spin chain. Its role
is narrower. The thermodynamic inequalities and real analyticity just
listed do not imply a uniform bound on the T^2 coefficient. Choosing
v=1 and a sequence of allowed integer CFT central charges, then adjusting
tau to their desired coefficient, likewise shows that these scalar
thermodynamic constraints alone cannot bound a putative A. No microscopic
CFT realization of that assignment is asserted.

Even a nonzero fixed high-temperature variance does not repair the
argument. Sum two such pressures with entropies S1,S2, S1+S2=log d,
small scale tau and a scale Delta. Choose Delta(tau) so that
S1 tau^2+S2 Delta(tau)^2 stays fixed. The first sector's low-temperature
coefficient diverges, while the high-temperature variance and an upper
energy bound stay fixed. Constants can be chosen so the same M bounds
all energies. This disproves the inference from capacity and the leading
high-temperature variance. Knowledge of every high-temperature
coefficient and a controlled analytic continuation would specify more
information; such control has not been derived from M and d.

[Proved conditional entropy bound.] Apply s <= log d to the lower
matching assumption stated in Section 1. Divide by
alpha_s pi T_s/3 to obtain the bound on B in (4). Since
c_a v_a <= v_max^2 c_a/v_a term by term, the bound on A follows.
This is the thermal conclusion from C69 specialized to the requested
velocity-weighted quantity. It cannot be strengthened by evaluating
the CFT entropy at an unspecified temperature of order M.

## 7. Ring ground energy and modular invariance

[Proved CFT statement, with the vacuum sector assumed.] For one sector,
the rectangular torus parameter is tau_torus = i beta v/L. The CFT
partition function is invariant under tau_torus -> -1/tau_torus, with
the corresponding sum over sectors or spin structures when required.
In the vacuum-dominated regimes this relates

    E0(L) = e0 L-pi c v/(6L)+o(L^(-1))

to

    log Z(L,beta)+beta e0 L
        = pi c L/(6 beta v)+o(L/beta).

The first coefficient is c v; the second is c/v. With several speeds,
there are several moduli i beta v_a/L. Sector by sector the relation
holds, but there is no single rotation into a microscopic chain at one
common rescaled temperature that replaces B by A. With unequal chiral
speeds one uses the corresponding left and right moduli.

[Proved lattice energy bound.] Let a_L be the ground energy of the
open chain of L sites, with L-1 bonds of norm at most M. Joining two
open chains adds one bond, so

    |a_(L+K)-a_L-a_K| <= M.

Tiling n blocks of length L and then dividing by nL proves that the
infinite-chain ground energy density e0 satisfies

    |a_L-L e0| <= M.

Closing one open chain into a ring adds one more bond. Hence for its
ground energy b_L,

    |b_L-L e0| <= 2M.                                (30)

If at a specified size L_* the lower Casimir estimate

    L_* e0-b_(L_*) >= alpha_E pi A/(6L_*)

holds, then

    A <= 12 M L_*/(pi alpha_E).                      (31)

The asymptotic Casimir formula alone ensures such a comparison at
sufficiently large sizes, not at a uniformly bounded size. Thus (31)
contains a spatial crossover scale.

[Proved finite-temperature capacity.] For any finite ring,

    b_L-T L log d <= F_L(T) <= b_L,
    1 <= exp(beta b_L) Z_L(beta) <= d^L.              (32)

Inserting controlled continuum asymptotics in (32) yields entropy or
finite-size inequalities of the preceding type. It does not remove
their scale requirements. The continuum modular identity holds after
both spatial and temporal lengths have entered the continuum regime.
A microscopic Trotter lattice is anisotropic and its rotated transfer
matrix is not the original finite-dimensional Hamiltonian. Setting
L=1 or beta M=1 in the continuum identity supplies an unproved crossover
assumption. Modular invariance is a relation among continuum coefficients,
not an exact ultraviolet duality of the given spin chain.

## 8. Checks, status of claims, and open step

[Proved dimensional check.] Under H -> lambda H, M,v,T,omega each
scale by lambda. The total current scales by lambda^2. Thus its
variance scales by lambda^4, in agreement with A T^3 in (11).
m1 has dimension energy cubed, agreeing with A T^2 and C3 in (2).
Both sides of (3), (4), and (31) have dimension energy. The midpoint
correlator has dimension energy to the fourth power. The susceptibility
in (20) has dimension energy.

[Proved low-q and zero-T checks.] As q -> 0,
S_h(q) -> pi B T^3/3, the thermal energy variance. Its imaginary-time
integral tends to beta times this value, pi B T^2/3. At T=0 the
density spectrum scales as |q|^3 and m1 as q^4, as required for an
operator of scaling dimension two. At fixed T > 0 the first moment
instead starts at T^2 q^2. At beta=0 every commutator has zero tracial
expectation, so the lattice m1 is zero, consistent with (1).

[Proved unequal-speed extension.] If a nonchiral sector has the same
Virasoro charge c_a in both chiralities but speeds v_(a,R),v_(a,L),
replace each c_a v_a^n by

    (c_a/2)[v_(a,R)^n+v_(a,L)^n],

for n = 1,-1,3. Thus A uses the arithmetic mean specified in the
background, but B uses the average of reciprocal speeds, not the
reciprocal of that arithmetic mean. The momentum-averaged spectral
formulas and bounds above then hold without change. The local midpoint
current correlator still counts c_a, with no speeds. The same reasoning
extends to separately listed chiral sectors using half their charge in
these nonchiral normalizations.

[Numerical checks only.] checks.py verifies the continuity equation,
the spectral representation of m1, positivity, and the norm estimates
on four small random nearest-neighbor chains. It also checks the
midpoint spatial integral and Fourier normalization by quadrature.
All final assertions pass; residuals are in CHECKS.md. An initial
Fourier check omitted the v^2 converting density to current, and was
corrected before the final run. These calculations do not test a
lattice-to-CFT matching hypothesis or a thermodynamic Drude limit.

[Conclusion about this route, not a no-go theorem.] Transport supplies
observables proportional to A, and the spectral moment gives the
unconditional microscopic inequality (1). Comparing them gives (3),
not the conjectured bound. Entropy gives (4); finite-size energy gives
(31). Convexity, the displayed high-temperature controls, and continuum
modular invariance do not supply the missing uniform onset scale.

A possible way beyond this analysis would be a new microscopic bound
on the q^2 coefficient of the inverse-frequency susceptibility (20),
or a monotone relation that compares its infrared stress contribution
with M log d without assuming a matching scale. Neither is proved
here. The obstruction is to the specific estimates used, not a proof
that every thermal or dynamical route must fail.

## 9. Sources

The supplied background and C69 fix the conjecture, gauge, entropy
capacity argument, and conservative propagation-speed constant. Their
settled model data were not recomputed.

Primary literature checked for conventions and context:

- K. Gawedzki, E. Langmann, P. Moosavi, [Finite-Time Universality in
  Nonequilibrium CFT](https://link.springer.com/article/10.1007/s10955-018-2025-x),
  J. Stat. Phys. 172, 353 (2018), particularly its thermal Drude convention.
- B. S. Shastry, [A Sum Rule for Thermal Conductivity and Dynamical Thermal
  Transport Coefficients](https://arxiv.org/abs/cond-mat/0508711),
  Phys. Rev. B 73, 085117 (2006). Equations (22)-(26) here are derived
  directly; no model-dependent thermal-operator formula is imported.
- D. Bernard and B. Doyon, [Energy flow in non-equilibrium conformal field
  theory](https://arxiv.org/abs/1202.0239), J. Phys. A 45, 362001 (2012),
  for the reservoir-current setting.

The derivations in this note give the normalization and assumptions
needed for its inequalities. No literature claim that all nonintegrable
thermal Drude weights vanish is used.
