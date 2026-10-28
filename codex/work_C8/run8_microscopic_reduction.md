# What the new core gives for the periodic chain

Status: exact finite construction and first-order reduction. Conditional
infrared consequence under the two stated continuum inputs. No microscopic
CFT measurement or unconditional counterexample is claimed here.

## 1. A concrete microscopic family

The new core meets the hypotheses of arrow_wrapper.md,
selector_wrapper.md and helix_equalizer.md. Those earlier finite results
are used as proved inputs, not rederived. With K=446400, put

    Q=K(q+1),  m=q^2/8,  k=m+6,  b=k+1,
    d=96Q+48.

Q denotes the core's full site dimension here, not its control count.
The complete on-site alphabet and every matrix element of the positive
core are specified by clean_readout_machine.py and clean_core_bond.py.
The latter evaluates an arbitrary column of its full bulk or boundary
matrix, including clock-illegal, role-illegal and multiple-head sectors.
It does not enumerate an exponentially large block Hamiltonian first.

Use arrow_wrapper equations (2)-(5), with this h_pay,a_pay,b_pay and
selected occupied length k. Then use selector_wrapper equations (4),(7)
and the two helix tracks in helix_equalizer equations (2)-(6). This
specifies the unperturbed periodic bond h_base and the complete base
code on cells of lengths b,b+1. All its coefficients are fixed once q
and k are fixed. They do not depend on the later ring length L.

Let Abar_a=A_a/alpha_Lf and Bbar_a=B_a/alpha_Rf for a=2f-1,2f, and
Dbar_f=D_Lf/beta_Lf. Lift them through the wrapper and apply the endpoint
helix compensators. Denote the resulting physical operators by Ahat,
Bhat,Dhat. The density lift is on the first internal bond, including
both possible selector letters after #, as in helix_equalizer Section 4.
It compensates the joint endpoint/neighbor occupation factor separately
from the single-site port occupation. Then their compressions at BOTH
cell lengths are Gamma_a, i Gamma_a P and Z_f, respectively.

For a desired quadratic payload, take

    H_pay = sum_(j,f) [-t_f(exp(i theta_f)c†_(j,f)c_(j+1,f)+h.c.)
                       +mu_f(n_(j,f)-1/2)],

and its real two-by-two block in the R tensor Gamma convention,

    J_f=-(t_f/2) [[sin theta_f,cos theta_f],
                  [-cos theta_f,sin theta_f]].

The following is a short specification of the full repeated microscopic
bond, using only the explicitly defined two-site terms above:

    h(epsilon)=h_base
        +epsilon [h_marker + kappa h_marker_number
                   +P_cut sum_(f,u,v) J_f[u,v]
                           Bhat_(2f-2+u) tensor Ahat_(2f-2+v)
                   -(1/2) sum_f mu_f Dhat_f].

Here h_marker is the compensated local H_0 to H_2 marker conversion of
helix_equalizer equation (10), divided by its positive projected overlap
so its retained hopping is the chosen J_0>0. h_marker_number is the
usual split on-site marker density. kappa selects an interior length
fluid density. P_cut includes X# and R#, not only R#. The field term
has its stated first-internal-bond projector. Every displayed term is
Hermitian and nearest neighbor. There is no explicit spatial period in h.

The actual positive core bulk is not small as a matrix, even at q=4 its
site dimension is 2232000, but its sparse column oracle is executable.
The full ring bond is equally a finite-matrix formula through the cited
wrapper definitions. A full violating instance has not been materialized.

## 2. A stronger exact symmetry

The machine preserves the integer charge d+1_(v odd)+x_head on every
operation. The code charge is sum_f n_f. Give all wrapper vacancies,
selector letters and helix tracks zero charge, and give occupied wrapper
letters the charge of their core record. Vacancy hopping transports that
charge unchanged, as do the marker conversion and all positive filters.
The compensators act only on neutral helix factors.

Each physical X_f+iY_f lowers the record charge by one. Although B uses
(Y,-X), the block J_f above combines the physical operators into a
linear combination of XX+YY and YX-XY on the two records. Each conserves
the sum of charges. Density readouts are diagonal in charge. Hence the
entire microscopic h(epsilon) has exact total U(1), in addition to exact
marker conservation. No exact U(1)^m is claimed. This improves the
parity-only hypothesis used in C13 and removes uniform pairing from the
allowed corrections without relying on momentum detuning.

## 3. The exact retained theory and explicit nonresonance

The existing projection and ring quotient give a free length fluid plus
the m-flavor quadratic payload at first order. It has m+1 oscillator
families. The quotient retains the common fermion spin structure and
combined marker/payload momentum condition. It applies to every
sufficiently large integer L because both b and b+1 cells are mobile.
There is no restriction to b dividing L. As in C13, finite-size zero-mode
backgrounds can drift with L; this is the usual incommensurate-chain
interpretation, not convergence to one fixed untwisted vacuum sector.

A fully explicit, deliberately inefficient choice avoids a genericity
assumption about Fermi momenta. Use the C19 geometric design with
alpha=exp(-1), simplifying the hopping amplitudes to t_f=1:

    p_(f,R)=pi alpha^(2f-1), p_(f,L)=pi alpha^(2f),
    k_f=(p_(f,R)-p_(f,L))/2,
    theta_f=-(p_(f,R)+p_(f,L))/2,
    mu_f=2cos(k_f),
    nu=1/2+alpha^(2m+1)/10.

The payload dispersion is -2cos(p+theta_f)+mu_f. Its two roots are the
displayed p's, with speeds 2sin(k_f)>0. The length reference speed is
2J_0(b+nu)sin(pi nu)>0 in physical coordinates. Thus the first-order
spectrum is z=1, with positive speeds, for each one fixed finite m.
No speed lower bound uniform in m is asserted.

The following elementary bounds make the finite momentum margins explicit.
For distinct exponents i in {1,...,2m},

    |sum_i sigma_i alpha^i| >= alpha^(2m)(1-2alpha)/(1-alpha)

whenever sigma_i=+1 or -1 and the sum is nonempty. The smallest exponent
term dominates the entire possible tail because alpha<1/2. Also the
absolute sum of all powers is alpha/(1-alpha)<1. Therefore the same
lower bound, times pi, applies to distance from 2pi Z for every distinct
signed two- or four-branch sum. Doubled branches have distance at least
2pi alpha^(2m). In particular one may use

    delta_pay=pi alpha^(2m)(1-2alpha)/(1-alpha)>0.

For a mixed length/payload bilinear,
|p_i +/-p_j|<=pi(alpha+alpha^2). Its distance from either sign of
2pi nu, modulo 2pi, is at least

    delta_mix=pi[1-alpha-alpha^2-alpha^(2m+1)/5]>0.

The potentially dangerous pure length harmonic h=1 is excluded by
0<nu<1. Higher pure harmonics have free dimension at least four on
resonance, as derived in C13. More generally a purported exact mixed
relation sum z_i p_i-2pi h nu=2pi l would give a rational polynomial
in alpha with coefficient -h/5 at alpha^(2m+1). Its transcendence forces
h=0 and then all z_i=l=0. Only the finite low-dimension margins are
needed for weak-coupling stability; an all-order Diophantine estimate is
not being assumed.

The C13 operator classification therefore has no non-Gaussian uniform
resonance through free dimension two. Exact total number conservation
further restricts that list. Neutral currents and their quadratic
products shift fillings, speeds and the positive Gaussian stiffness.
All remaining allowed resonant vertices have free dimension at least
three. At each fixed finite m there is consequently a nonzero Gaussian
stability neighborhood, conditional on the standard RG input below.

## 4. What a finite-coupling conclusion still uses

There are two distinct inputs.

M, massive-sector integration. In the mobile-cell model, after keeping
propagation in the entire retained manifold unexpanded, integrating out
local machinery defects produces exponentially local spatial kernels
and short-time kernels. Their norm after dividing by epsilon is at most
C epsilon B^2/Delta, and their frequency dependence is suppressed by
omega/Delta. B includes the FULL physical norms of normalized ports,
fields, marker conversion and compensators. C and Delta may depend on
the chosen finite q,k, but not on L. Defect costs stay positive at
sufficiently small epsilon.

G, Gaussian infrared stability and identification. A sufficiently small
number-conserving, quasi-local deformation of the mobile m+1-mode
reference, with the stated momentum margins and positive stiffness,
flows to a compact Gaussian theory with m+1 oscillator sectors and its
specified ring zero-mode quotient. It therefore has c=m+1 and z=1.

G is the standard bosonization/RG conclusion used in C13 and ir_control.
M is the usual massive-mode elimination used in the later
research_continuation/IR_AUDIT, but that file explicitly does not prove
the connected-excursion convergence estimate for moving cells. It is
stronger than the finite-order local SW theorem. Neither M nor G was
proved for this new microscopic bond or numerically verified here.

There is substantial microscopic support for M, beyond an isolated
finite-volume gap. In a fixed legal marker pattern the established
wrapper theorem gives the additive operator lower bound

    H_base-eL >= Delta_* sum_(selected cells j) Q_j
                    +gamma sum_(wrong cells j)(1+|l_j-b|),
    gamma=t eta_b-U>0.

The selected-cell projectors are on disjoint finite cells. A wrong cell's
cost also obeys gamma(1+|l-b|)>=gamma l/(b+1), so arbitrarily long wrong
cells pay for their volume. A local perturbation creates or moves defects
in a bounded number of adjacent cells. These facts are the microscopic
inputs to a connected-excursion estimate. They do not by themselves
prove its resummation in a critical retained sector.

Assuming M and G, choose a positive tolerance tau inside the Gaussian
neighborhood, including the smallest velocity and all momentum margins.
At fixed q,k and all physical normalizations, choose one fixed epsilon>0
with epsilon B<<Delta and C epsilon B^2/Delta<tau. Then the normalized
retained action lies in that neighborhood for all sufficiently large L.
This gives the conditional conclusion c=m+1. Sending epsilon to zero
with L would not establish it and is not part of the prescription.

This conditional reduction is complete at the finite-dimensional end.
Whether M can be accepted without a separate model-specific argument is
now the decisive issue. The report does not promote a finite-order SW
error proportional to L into control of conformal levels proportional
to 1/L. No additional compiler or clean-readout assumption is hidden in M.

## 5. A concrete resource inequality, not a measured violation

For K=446400, take q=1536K. Exact integer arithmetic gives

    q = 685670400,
    m = 58767987179520000,
    k = 58767987179520006,
    Q = 306083267006400,
    d = 29383993632614448.

The conditional central charge m+1 exceeds d-1 by
29383993546905554. These integers are in clean_readout_resources.json.
They are an existence-scale example, not a practical simulation proposal.
The construction gives d=O(q), c_projected=O(q^2), so optimizing this
large controller constant is not needed for the resource inequality.

The conclusion is conditional on M and G, with no microscopic central
charge estimator. Under the task's verified-result criterion, this run
is recorded as nonconclusive. Its new result is the explicit initialized
core and the resulting finite-matrix candidate family. No full-class
upper bound and no unbounded-c construction at fixed d follows.

## 6. Next work, ranked

1. Independently audit the 58-control machine and its microscopic bond
   oracle, particularly counter restoration, reversed collision order,
   the I-pins, and every partial-time READY convention. The new proof
   and tests are by the same author; no independent referee is claimed.
2. Establish M for this concrete number-conserving first-order family,
   or provide a precise standard theorem whose hypotheses it meets.
   The additive defect-energy bound above is available. A proof of only
   finite-order SW or prethermal behavior is insufficient.
3. Apply the number-conserving version of C13's Gaussian analysis to
   that effective action. Keep the common ring quotient, relative error
   scale, velocity minimum and drifting zero modes explicit. A numerical
   check would require a much cheaper equivalent core or cluster work;
   a direct violating-chain calculation is not a reasonable local job.
4. Separately, use the completed port-touching controller interface to
   revisit C19's source filters. Ground-code readout alone does not
   implement their resolvent, so fixed-d unbounded c remains open.
