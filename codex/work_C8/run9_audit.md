# Ninth-run audit, explicit periodic matrix and polynomial scales

2026-10-28. Closed-book continuation by Codex, GPT-6. No independent
referee, network, cluster job or external write was used.

## 1. Verdict

Nonconclusive under the task's definition. The eighth run's initialized
core survives the present source audit and new composition checks. The
complete periodic bond now has an executable sparse column oracle. Its
reference target can have payload velocities bounded away from zero,
detunings of order m^(-3), and perturbation norm bounded by a polynomial
in m. Its finite clock gap is still very small.

The remaining thermodynamic step has not been proved. A new soluble
example shows why an inverse on the global bare excited-space projector
cannot provide that step, even when the physical massive sector is a
collection of independent, uniformly gapped spins. Connected local
elimination is still a separate requirement. This example does not
disprove local effective field theory or the candidate.

The new results are (i) a finite-matrix implementation and audit,
(ii) an elementary number-conserving momentum design,
(iii) analytic normalization bounds, and (iv) a precise obstruction to
a proposed elimination shortcut. None is a verified microscopic CFT.

## 2. The complete microscopic bond

[proved finite specification] Use the core of clean_readout_core.md,
with q divisible by four,

    m=q^2/8, k=m+6, b=k+1, Q=446400(q+1).

The on-site basis of periodic_clean_bond.py is

    selector x wrapper x equalizer A x equalizer B,
    dim = 4(6Q+3)4 = 96Q+48.

Selector labels are #,L,X,R. Wrapper labels are H0,H1,H2 and
L,B,R, each occupied tag carrying a Q-state core record and one filter
qubit. The two equalizer bits are present also on vacant sites. The
encode/decode functions specify the order of every basis vector.

The functions CoreBond.bulk and CoreBond.boundary supply the complete
positive core matrices. PeriodicBond.wrapper implements arrow_wrapper
equation (4), including its oriented boundary penalties on all sectors.
PeriodicBond.column implements selector_wrapper equation (4), both
equalizer helix bonds on internal selector bonds, both pins on cuts,
and all perturbations described below. It returns a dictionary of row
pairs and matrix entries for an arbitrary column pair. Neither input
nor output is restricted to a legal clock, role, marker or head sector.

Choose U=1, vacancy hopping w=0.9/(k+1), filter scale as in
length_filter.md, and equalizer pitch phi=1/b. Selector coefficients
are the explicit functions in selector_wrapper equation (7). Charge
ground vectors in the one- and two-vacancy sectors determine the
endpoint factors r_p, the first-two-sites occupation factors s_p,
and the cut overlap beta_wrap. They are unique positive ground vectors
of the finite charge matrices specified in vacancy_parent.py and,
independently, in positive_charge_ground in the new oracle.

The new numerical routine works in the even reflection sector and
uses positive inverse iteration. Ordinary dense diagonalization at
k=14 lost positivity in very small wavefunction tails because the
even/odd tunneling doublet is nearly degenerate. The even-sector
Stieltjes inverse avoids that problem. No negative amplitudes were
clipped. This changes coefficient evaluation, not the mathematical
definition of the charge ground state.

For definiteness use the target of Section 3 with t_f=J0=1 and nu=1/2.
Let C_uf be the raw charge-lowering record map

    C_uf = sum_(p=0,1) (-1)^p |e_a><o_(b+p)|,
    a=f div(q/4), b=2[f mod(q/4)],

with the local READY, no-head and endpoint-role projectors. Lift it
onto the appropriate wrapper endpoint, apply the helix compensator
and divide by alpha_uf. Write the resulting operator as Chat_uf.
The Hermitian cut interaction is

    v_cut = -P_cut sum_f [exp(i theta_f) Chat_Rf† tensor Chat_Lf
                        +exp(-i theta_f) Chat_Rf tensor Chat_Lf†].

This equals the J_f contraction in run8_microscopic_reduction.md.
To check, put X=C+C† and Y=i(C†-C). The right pair is (Y,-X), and

    sum_uv J_f[uv] B_u tensor A_v
      =-(cos(theta_f)/2)(X tensor X+Y tensor Y)
       -(sin(theta_f)/2)(Y tensor X-X tensor Y)
      =-exp(i theta_f) C† tensor C-exp(-i theta_f) C tensor C†.

The direct charge form is used in the code. Floating-point cancellation
between the four Hermitian components initially left forbidden entries
of size 7.45e-9 in a matrix whose allowed entries are much larger.
The charge form makes every forbidden entry absent by construction.
This was an implementation issue, not a failure of the exact symmetry.

The field is -sum_f(mu_f/2) Dhat_f on the first internal bond, with
the separate density count beta_Lf and compensation 1/s_p. Both #L
and #X are enabled. The marker move is the typed-vacancy conversion
with its local equalizer phase, divided by

    zeta=A_b^2 beta_wrap cos(phi/2)^(2b).

It therefore has retained hopping J0=1. Its matrix entries also
conserve total microscopic charge and the number of markers.

The remaining marker chemical potential is explicit. For the target
widths k_f, the payload energy per squeezed cell is

    e_pay=sum_f[-2 sin(k_f)/pi+mu_f(k_f/pi-1/2)].

Set kappa=2/pi-e_pay and use the split on-site marker density. At
fixed length density nu, the projected energy per physical site is

    e_1(nu)=[kappa+e_pay-2 sin(pi nu)/pi]/(b+nu).

Its derivative vanishes at nu=1/2. More generally its numerator after
differentiation is

    -2(b+nu)cos(pi nu)-kappa-e_pay+2sin(pi nu)/pi,

whose derivative is 2pi(b+nu)sin(pi nu)>0 in (0,1). Thus this is the
unique projected thermodynamic minimum, with positive curvature.
This calculation uses the complete field constant -mu_f/2. Omitting
it would select the wrong marker density.

All microscopic coefficients are fixed before the ring length is
varied. The constructor evaluates them numerically only for modest q.
At large q their finite spectral definitions still specify matrices;
the code is not a claim to evaluate an astronomical selected cell on
one machine. In particular no finite nonzero coupling has been
certified to produce the predicted CFT.

## 3. Rational Fermi points with uniformly nonzero payload velocities

[proved algebraic theorem] For every m>=1 define N=2m, B=2N+1,

    a_i=i+B i^2, 1<=i<=N,    A=a_N,
    p_i=pi/4+pi a_i/(8A).

Assign p_(f,L)=p_f and p_(f,R)=p_(m+f), for f=1,...,m. Set

    k_f=(p_(f,R)-p_(f,L))/2,
    theta_f=-(p_(f,R)+p_(f,L))/2,
    mu_f=2 cos(k_f).

Then the dispersion -2cos(p+theta_f)+mu_f has these two roots and

    pi/64 < k_f < pi/16,
    2sin(pi/64) < |v_f| < 2sin(pi/16).

The lower bound is independent of m. In physical coordinates the
payload speed is multiplied by b+nu, as in the inherited squeezed
coordinate dictionary. The length speed for J0=1 and nu=1/2 is
2(b+1/2), so it also remains positive.

Proof of the width bounds. The smallest numerator difference is

    a_(m+1)-a_1=m+B(m^2+2m)>A/4,
    A=2m+4Bm^2.

Every such difference is less than A. Divide by 16A and multiply
by pi. The velocity statement follows by monotonicity of sine in
the displayed interval.

[proved] All unordered pair sums a_i+a_j, including i=j, are distinct.
If two sums agree, reducing modulo B gives i+j=r+s, since all these
sums lie between 2 and 2N<B. The remaining equality gives
i^2+j^2=r^2+s^2. The two elementary symmetric polynomials agree,
so {i,j}={r,s}. Consequently any nontrivial two-branch difference
or two-minus-two branch sum obeys

    dist(momentum,2pi Z)>=delta_pay=pi/(8A).

There is no modular wrap: all p_i lie in (pi/4,3pi/8], so branch
differences have magnitude less than pi/8 and pair-sum differences
have magnitude less than pi/4. Since A=16m^3+4m^2+2m, the stated
margin is polynomial, of order m^(-3).

[proved conditional on the inherited bosonization dictionary] Exact
total U(1) restricts a local payload vertex with integer chiral charges
z_i to sum_i z_i=0. At the free point its dimension is

    Delta_pay=(1/2) sum_i z_i^2.

The nonzero vertices with Delta_pay<=2 are neutral bilinears and
neutral quartics. The pair-sum theorem excludes their uniform
non-Gaussian resonances. Equal creation/annihilation multisets instead
give densities and current products, which deform the Gaussian action.
Pauli-forced derivatives increase the dimension, so they do not add
another dangerous case. Integrality and sum z_i=0 make Delta_pay an
integer. Every remaining nonzero resonant vertex has Delta_pay>=3.

The mobile locality dictionary of C13 gives, at a uniform resonance,

    Delta_ref=h^2+Delta_pay,
    sum_i z_i p_i-2pi h nu=0 modulo 2pi.

At nu=1/2 the only new dimension-two possibility is h=+/-1 times a
neutral payload bilinear. Since |p_i-p_j|<pi/8, its mismatch is
strictly greater than 7pi/8. Pure h=+/-1 is nonresonant, and |h|>=2
has dimension at least four. Thus all non-Gaussian resonant vertices
have reference dimension at least three, including the mixed sector.

This replaces the geometric momentum prescription of RESULT Section
41 for this total-number-conserving family. It does not require the
sparse species torus of the separate collaborator construction.
Rational higher-order resonances are allowed. Their free dimension is
at least three, which is sufficient for the conditional Gaussian
stability argument. No all-order irrationality assertion is needed.

[assumed infrared input] Apply the positive Gaussian stability
criterion to the actual dressed parameters, including mixed
density-current terms and the inherited ring quotient. The matrix
certificate in contrib/digit_compilation/audit/GAUSSIAN_CERTIFICATE.md
is available as a source result: a relative Gaussian perturbation
delta<5/13 keeps every reference-dimension-three vertex irrelevant.
That is a continuum matrix statement, not a bound on the corrections
generated by this microscopic bond. Dressed momenta must separately
retain their nonzero margins. This run does not rederive that source
result or use it to assume the missing elimination.

## 4. Polynomial lower bounds for every readout visibility

[proved] For the autonomous core, with m>=2,

    alpha_uf >= 1/(80m),    beta_uf >= 1/(80m)

for both ports u and every flavor f. The exact values remain those
computed by clean_readout_counts.py; these are uniform lower bounds.

Proof. Let n=m+6, T=11m^2+105m and S=2^n. The saved machine proof
gives T<=S. Every history therefore contains floor(S/T) complete work
cycles. During a readout walk, a move between two sites strictly
inside the physical endpoints leaves the port active at all 2n+1
vertices of its complete clock sweep. There are 2(n-3) such moves
per flavor and port. For density, exclude the adjacent color site
as well, leaving 2(n-4) moves per flavor and port. None of these
moves changes READY or a record.

The total clock vertex count obeys

    C_n=(2n+4)S+n-4 <= (2n+5)S.

Using floor(S/T)>=S/(2T), and ignoring additional active carry
and overflow vertices, gives

    alpha_uf >= (n-3)(2n+1)/[T(2n+5)],
    beta_uf  >= (n-4)(2n+1)/[T(2n+5)].

For m>=2, T<=64m^2, n-4=m+2>=m, and
(2n+1)/(2n+5)>=4/5. Both bounds imply 1/(80m).
The counts are independent of the logical input by the all-time
clean-readout theorem. Truncated later cycles only add nonnegative
contributions.

This removes an apparent exponential normalization cost from the
history length. It does not remove the exponentially small clock gap.

## 5. Wrapper occupation weights, transport, and a complete norm bound

[proved] For k>=8 and 0<w<U/(k+1), both vacancy sectors p=1,2 obey

    r_p>1/4, s_p>1/4, beta_wrap>1/16.

Here r_p is the endpoint occupation weight and s_p is the weight
of configurations in which the first two sites are occupied.

Proof. In the p-hole basis write H_charge=U D-w A, where D counts
vacancy/particle domain walls and A is the hole hopping adjacency
matrix. Its degree is at most 2p, so ||A||<=2p. The only configurations
with D=1 are the left and right contiguous endpoint droplets. Every
other configuration has D>=2. A single endpoint droplet is a trial
state of energy U. Therefore the ground vector chi satisfies

    U >= <H_charge> >= U(1+W_non-droplet)-2pw,
    W_non-droplet<=2pw/U.

Perron positivity and reflection invariance give equal weights for
the two endpoint droplets. Each has weight at least
(1-2pw/U)/2>1/4, since p<=2 and k>=8. The right droplet leaves the
first two sites occupied, proving the r_p,s_p bounds. In the positive
overlap sum for beta_wrap, retain just the term pairing the last-site
one-hole droplet with the last-two-sites two-hole droplet. Squaring
gives the product of their weights, greater than 1/16.

[proved] At phi=1/b the positive endpoint and density compensators
satisfy

    ||M_L||, ||M_R||, ||M_D|| <= 8pi^2 b^2.

Indeed their norms are at most the sum of two requested inverse
occupation weights divided by sin^2(phi/2). Each inverse weight is
less than four, and sin(1/(2b))>=1/(pi b).

[proved] The complete marker visibility obeys zeta>1/64, so the
coefficient needed for retained hopping J0 is less than 64J0.
For t=b/(b+1), concavity of sine on [0,pi] gives

    sin(t x)>=t sin(x),
    A_b>=t^(3/2),    A_b^2>1/2.

Also Bernoulli's inequality and cos x>=1-x^2/2 give

    cos(1/(2b))^(2b)>=1-1/(4b)>1/2.

Combine these two inequalities with beta_wrap>1/16.

[proved conservative bound] Let v be the full physical perturbation
bond at J0=t_f=1, before multiplying by epsilon. Its norm is at most

    B_m =129+2m+1638400 pi^4 b^4 m^3+640 pi^2 b^2 m.

This bounds the actual full-space matrix, including illegal sectors.
To see this, each raw C_f has norm sqrt(2), and hence

    ||v_cut|| <=4 ||M_R|| ||M_L|| sum_f 1/(alpha_Rf alpha_Lf).

The colored density supports for different flavors are disjoint, so

    ||v_field|| <=||M_D|| max_f |mu_f|/(2 beta_Lf).

The marker terms cost at most 2u+|kappa|. With |mu_f|<=2,
|e_pay|<2m, |kappa|<=2/pi+2m, u<64, and the preceding bounds,
the displayed polynomial follows. No small norm was assumed from
the visibility of a projected matrix element.

Thus the physical perturbation normalization is at worst O(m^7).
The much smaller unperturbed gap, and the unknown all-order local
elimination constants, remain in any finite-coupling hierarchy.

## 6. Why the global bare-Q resolvent is insufficient

[proved obstruction to a shortcut] Consider L independent two-level
ancillas,

    H0=Delta sum_j |1><1|_j,
    V=g sum_j X_j,    H=H0+V,
    P=|0...0><0...0|, Q=I-P,

with Delta>0 and any fixed nonzero real g. Their exact massive gap is
Omega=sqrt(Delta^2+4g^2), independent of L. Define

    e_-=(Delta-Omega)/2,
    p=(1+Delta/Omega)/2 in (0,1),
    E0=L e_-.

Then

    0<lambda_min[Q(H-E0)Q on ran Q]
      <=(-L e_-) p^L/(1-p^L),

which tends to zero exponentially with L.

Proof. The exact ground state is the product psi=|g>^tensor L and
||P psi||^2=p^L. Since (H-E0)psi=0, one has

    <Qpsi|(H-E0)|Qpsi>=<Ppsi|(H-E0)|Ppsi>
                       =(-L e_-)p^L.

Divide by ||Qpsi||^2=1-p^L. Strict positivity at finite L follows
from the unique full ground state having nonzero overlap with P.
The same example may be tensored with any independent critical
retained chain; neither the bound nor the massive gap changes.

This model meets additive defect-energy and bounded local-insertion
hypotheses in their strongest form. Nevertheless propagation in the
global bare Q space, after resumming QVQ and shifting to the true
ground energy, does not decay at a rate bounded below independently
of L. Thus that propagation cannot be inserted as a uniform fast
kernel in a Feshbach argument.

The example is also a control showing what a correct integration must
do. An on-site rotation diagonalizes each ancilla, giving the scalar
energy density e_- and a gap Omega. Its local expansion converges
near g=0. The failure is in the global projector, not in local
integration or in the existence of an effective theory.

In fact the connected energy per site is explicitly

    e_-=sum_(r>=1) (-1)^r Cat_(r-1) g^(2r)/Delta^(2r-1),
    Cat_j=binomial(2j,j)/(j+1).

The binomial series for the square root proves the identity. The bound
Cat_(r-1)<=4^(r-1) proves absolute convergence for |g|<Delta/2,
uniformly in volume after division by L. Thus this same example has
both a vanishing global bare-Q gap and a well-controlled local
connected expansion. A global-resolvent obstruction must not be
reported as a no-go for the required local theorem.

Consequently an excursion defined only by "the system has not yet
returned to P" is not necessarily a connected local excursion.
It includes newly created defects arbitrarily far from an existing
one. Their cancellation requires an actual connected-cumulant or
local dressing construction. Merely retaining propagation in P
unexpanded does not demonstrate that cancellation. The geometric
connected-excursion estimates in the collaborator IR_AUDIT remain
stated inputs until that construction and its bound are established.
This is consistent with that file's own qualification that it is
not a convergence theorem.

## 7. A finite-order bound that does follow

[proved, finite order] For the new periodic bond, second-order
elimination around H0 is local uniformly in L. Write E_base=eL,
P for its complete code and R=Q(H0-E_base)^(-1)Q. The inherited
uniform gap gives ||R||<=1/Delta_base. For local v_i,

    H2=-sum_ij P v_i R v_j P,
    ||P v_i R v_j P||<=B_m^2/Delta_base.

There is also a support restriction. Fix an initial legal selected
marker pattern and order its cells cyclically. A single v_j changes
only one cell or two neighboring cells. A marker move shifts their
boundary by one site and leaves all other cells untouched. The
unperturbed Hamiltonian fixes the new marker positions and is a sum
of independent cell Hamiltonians. Every untouched cell remains in
its zero code during the resolvent. Thus the inverse in a single
excursion acts only on the affected cell patch, with a denominator
independent of spectator code states.

If v_i touches a disjoint cell patch, it cannot erase the nonzero
excitation left by Qv_jP. Final projection onto P therefore kills
that term. Nonzero pairs have overlapping patches. In particular
they involve at most three consecutive initial cells, or a single
boundary displacement and the cells meeting it. The conservative
physical radius R0=4(b+2) around the first insertion contains every
such process, including the one-site displacements.

Summing over first-bond anchors and second insertions in that radius
gives a safe interaction-norm bound

    ||H2||_local <= (2R0+3)^2 B_m^2/Delta_base.

Here the local norm sums norms of the represented patch interactions
whose support contains a specified physical site. The marker-pattern
projectors needed to specify a selected patch lie within the same
radius; direct summation over their orthogonal sectors does not
introduce a factor depending on the number of tilings of the ring.
For rings shorter than twice this radius the whole finite ring is
the patch, which obeys the same conservative counting bound.

This proves a finite-order statement, not convergence of the
connected expansion. At higher orders a patch can grow, and the
necessary linked subtractions and critical propagation must be
controlled. A fixed-order error in an extensive Hamiltonian does
not resolve conformal levels of order 1/L. No claim that this last
bound proves M is made.

## 8. Numerical checks and reproduction

[numerical finite checks] Run with single-thread BLAS:

    python3 rational_target.py
    python3 periodic_bond_checks.py
    python3 normalization_bounds_checks.py
    python3 massive_projection_checks.py

The rational arithmetic tests use m=1,2,3,8,18,32,64,128,256.
The largest test has 512 branches and 131328 unordered pair sums,
all distinct by integer arithmetic. Velocity checks use the stated
dispersion, not a central-charge estimator.

At q=4 the complete periodic site dimension is 214272048. The oracle
is checked on 3232 input columns separately for the base and
perturbation matrices. The samples include random full-space inputs
and exhaustive active endpoint record/equalizer combinations, both
cut selector letters, both first-internal-bond letters, and both
marker-move orientations. Every nonzero tested entry conserves total
charge and marker number. All selector changes preserve the three
overlapping legality indicators for every external selector pair.

An independent dense wrapper builder agrees on all 225 columns of a
two-state payload fixture within 9.1e-13. All four q=4 histories,
each of 5124 vertices, are contracted using the new oracle's actual
port and field definitions. The resulting 16-dimensional two-cell
Hamiltonian agrees with an independently assembled fermionic Fock
matrix with zero measured error. Helix compensation errors are below
7e-14. The projected marker-field minimum is nu=1/2.

Normalization bounds are checked at q=4,8. The q=8 calculation was
the test that exposed loss of Perron-vector positivity in the older
unrestricted diagonalization. The replacement solver has residuals
below 5e-13 by assertion and keeps every component positive.

The massive-projection example uses Delta=1, g=1/2 and
L=4,8,16,32,64,96,128. Permutation symmetry reduces the bare-Q
calculation to an L-dimensional tridiagonal matrix. At L=128 its
gap is about 3.98e-8, below the proved 4.18e-8 variational bound,
while the full massive gap is sqrt(2) at every size. Perron
uniqueness makes the lowest bare-Q vector permutation symmetric,
so this reduction includes its actual lowest eigenvalue.

No microscopic CFT estimator, thermodynamic simulation, external
source verification or independent referee was performed.

## 9. Ranked continuation

1. Prove a connected local massive-sector reduction for the concrete
   sparse periodic matrix. Use local dressed degrees of freedom or
   connected cumulants. The global bare-Q inverse is insufficient,
   as Section 6 demonstrates even for independent massive spins.
2. Carry the resulting interaction bound into the total-U(1) Gaussian
   theory, retaining actual dressed fillings and the joint ring
   quotient. Sections 3-5 remove exponential velocity and port-norm
   costs, but do not replace this step.
3. Independently audit the full periodic oracle and its core. The new
   finite tests are a second implementation path by this agent,
   not an independent referee report.
4. For numerical central charges, first find a substantially smaller
   equivalent microscopic implementation. Even q=4 has over two
   hundred million on-site states and is not a feasible direct local
   many-body job. A cluster alone does not make the violating
   example practical.
5. Keep C19's fixed-d spectral problem separate. This first-order
   family has growing d and does not establish unbounded c at fixed d.
