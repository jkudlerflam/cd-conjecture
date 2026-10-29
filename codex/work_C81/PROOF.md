# C81. The 16-state one-vacancy lift

[Heuristic, conditional conclusion] This family is a candidate with
c=15+1=16>d-1=15. The bond and finite identities below are explicit.
The microscopic conformal conclusion still requires the assumptions in
Section 9. No central charge is measured here.

The tags have the following meaning. [Proved] denotes an algebraic result
or a derivation with its hypotheses stated. An inherited computer-assisted
result is identified as such. [Numerical] denotes a resolved finite-matrix
check. [Heuristic] denotes a continuum mechanism. [Assumed] denotes an
input still needed to connect that mechanism to this microscopic model.
Definitions are identified separately and are not empirical claims.

## 1. One fixed bond and its charges

[Definition] There are sixteen orthogonal local states

    vac;
    L_a, R_a, D_(A,a), D_(B,a), X_(A,a), X_(B,a), a=0,1;
    M_0, M_A, M_B.

The fifteen occupied states span H_occ. Write iota for its inclusion,
n=I-|vac><vac|, and m for the projector onto the three marker states.
The integer on-site charge q is zero on vac and is

    q(L_a)=q(R_a)=2a,
    q(D_(A,a))=q(X_(A,a))=2a-1,
    q(D_(B,a))=q(X_(B,a))=2a+1,
    q(M_0)=0, q(M_A)=-1, q(M_B)=1.

[Definition] Here is the complete occupied bond h15. Allowed ordered pairs
have a marker or L on the left, or have R or M_0 on the right. Every other
pair pays K_pen. Put z_A(a)=(-1)^a and z_B(a)=1. Each displayed move includes
its Hermitian conjugate. Sum all bit and species indices shown.

| Input | Output | Amplitude |
| :-- | :-- | :-- |
| R_a R_b | L_a R_b | -1 |
| R_a M_0 | L_a M_0 | -1 |
| D_(s,a) R_b | L_a D_(s,b) | -z_s(a) |
| X_(s,a) R_b | L_a X_(s,b) | -exp(i phi_s) z_s(a) |
| X_(s,a) M_0 | L_a M_s | -exp(i phi_s) z_s(a) |
| D_(s,a) M_0 | L_a M_s | -t_d exp(i psi_s) z_s(a) |
| m D_(s,a), each marker m | m X_(s,a) | -1 |
| M_0 R_a | M_A D_(B,a) | lambda exp(-i pi/4) |
| M_0 R_a | M_B D_(A,a) | lambda exp(+i pi/4) |

Let h_moves contain the first seven rows and V the last two with lambda
removed. The on-site operator o assigns Delta to D,X and excited markers,
adds D to X and excited markers, assigns mu_sel to all markers, and assigns
f a to every data letter. It also contains

    g (|D_(B,0)><D_(A,1)| + h.c.).

Then

    h15 = h_moves + K_pen Pi_bad + lambda V
          + (o tensor I + I tensor o)/2,                         (1)
    h16 = -t sum_(a in H_occ) (|vac,a><a,vac| + h.c.)
          + J (iota tensor iota) h15 (iota tensor iota)^dagger
          - mu (n tensor I + I tensor n)/2.                      (2)

[Proved] Equations (1)-(2) are Hermitian nearest-neighbor bonds, repeated
unchanged around every physical ring. The dimension is 16 and the bond
matrix is 256 by 256. Every listed transition conserves Q=sum q, occupied
number N=sum n, and marker number M=sum m. These are independent exact
global U(1) charges. The relation N=16M is not an operator identity. It is
the proposed low-energy length locking. The other two charges do not imply
three independent charge fluids.

[Proved] Each half-site term in (1) stays inside the occupied-pair block
in (2). For an isolated occupied particle it contributes nothing. Moving
o to a physical on-site term would change the model and invalidate the
uniform first-order rescaling below. The new mu couples to N; mu_sel still
couples to markers only when inherited through the occupied contacts.

[Proved] Physical marker positions move under vacancy hopping. Their order
relative to all other occupied particles cannot change. In squeezed
coordinates the marker positions are conserved up to cyclic relabeling.
The full ordered word is not conserved, since h15 changes its letters.

## 2. Parameter schedule and inherited certificate

[Definition] Fix n=15, b=n+1=16, t=1, unit data/express/switch hoppings, and

    D=1000, Delta=4-2cos(pi/17),
    phi_A=pi/30, phi_B=0, psi_A=pi, psi_B=pi/2,
    d_0=1, d_1=1004, d_(j+1)=1004 d_j-d_(j-1),
    beta=1/d_16, a0=sqrt(2/17)sin(pi/17), t_d=beta 2^-84,
    a_l=2[1-cos(pi/(l+1))],
    delta_l=(l+1)a_l-l a_(l+1),
    mu_sel=2-(delta_16+delta_15)/2.

Here a0 is a code endpoint amplitude, unrelated to the contact factor
alpha used later. A literal starting member of C75's certified family is

    g=2^-256, K_pen=2^1152, lambda=2^-640,
    E_unit=a0^4 beta^2 g^2, f=lambda^2 E_unit 2^-140.       (3)

[Proved, inherited computer-assisted input C75] The exact second-order
quadratic projection at finite g,K_pen has 30 simple Fermi crossings,
with delta_2>0.0145540 and delta_4>1.1444e-7. The latter compares all 465
unordered pairs with replacement. The normalized finite-g and finite-K
error bounds are respectively below 1.137e-15 and 1.293e-28, against
an allowed Bloch error 2^-32. These results are taken from
[C75 PROOF, Sections 4-6](../work_C75/PROOF.md) and
[C75 STATE](../work_C75/STATE.md), without rerunning the certificate.
They certify the quadratic reference, not a many-body infrared limit.

[Definition; proved monotonicity of the inherited error bounds] To leave
room for unknown finite matching constants, use the nested family

    g_r=2^(-256-r), K_r=2^(1152+4r),
    lambda_(r,u)=2^(-640-u), f=lambda_(r,u)^2 a0^4 beta^2 g_r^2 2^-140,
    J_v=2^-v,  r,u,v nonnegative integers.                  (4)

Choose r first, then u at that fixed r, then v after the complete h15 is
fixed. The finite-g error scales as g^2 and the normalized finite-K error
as 1/(K g^2), so both improve by 2^-2r. The unnormalized nonquadratic
finite-K bound likewise improves relative to E_unit. Decreasing g at
fixed K indefinitely would not have this property. This is a necessary
change from the d24 schedule. Decreasing lambda alone also does not
suppress the order-g^4 interaction relative to the order-g^2 band.

[Assumed] Choose r,u sufficiently large to put the complete normalized
static and retarded payload interactions inside the multichannel stable
region, not merely inside the quadratic certificate. Neither C75 nor this
report gives the numerical matching constants needed to certify r=u=0
for that purpose. Section 8 specifies a sufficient conditional choice of v.
After these choices all coefficients are finite, fixed, and nonzero before
the thermodynamic limit. The powers in (3)-(4) must be retained symbolically;
K_r exceeds ordinary double precision even at r=0.

[Assumed, density tuning] Choose mu so the interacting infinite-volume
occupied density is 1/2. This uses the positive compressibility of the
surviving charge fluid. At J=0, mu=0. If e_h is the exact selected material
ground energy per occupied site, first-order thermodynamics gives

    e_phys(rho)=-2t sin(pi rho)/pi + J D_c(rho)e_h+o(J),
    D_c(rho)=rho^2-sin^2(pi rho)/pi^2,
    mu=J e_h+o(J) at rho=1/2.                              (5)

The derivative in (5) uses D_c'(1/2)=1. It is a conditional leading tuning
formula, not a computed finite-J chemical potential. mu does not depend
on the physical ring length.

## 3. Exact squeezing and the cyclic projection

[Proved] Fix 2<=N<L. Order occupied sites x_1<...<x_N and keep their labels
in that order. This is a basis bijection with

    C^(choose(L,N)) tensor (C^15)^tensor N.

A hop through the closing physical bond rotates the word by
R|a_1,...,a_N>=|a_N,a_1,...,a_(N-1)>; interior hops leave it unchanged.
Writing A_j=1[x_(j+1)=x_j+1] for j<N and
A_N=1[x_1=0,x_N=L-1], the interaction is

    V_occ=sum_j A_j tensor h15_(j,j+1).                    (6)

Bond orientation is preserved, including j=N. In an R=e^(i theta_r)
sector, theta_r=2pi r/N, the free charge fermions have

    k_a=[2pi a-theta_r-pi(N-1)]/L,
    epsilon(k)=-2t cos k.                                 (7)

The parity shift pi(N-1) is required. Occupying N consecutive momenta
a=ell,...,ell+N-1 gives the sea. Minimize with ell=0 for r<N/2 and ell=1
for r>N/2. At r=N/2 retain both branches.

[Proved] The physical translation phases of the retained seas are
exp[2pi i(r-N ell)/L]. They are distinct because N<L, including the two
half-flux branches. Translation invariance therefore removes cross-branch
matrix elements of (6). Within one branch,

    Pi_r h15_(j,j+1) Pi_r = Pi_r H15 Pi_r/N.

The free sea has

    <n_x n_(x+1)>= (N/L)^2-[sin(pi N/L)/(L sin(pi/L))]^2.

Summing contacts and dividing by N gives the exact finite-ring coefficient

    alpha_(L,N)=N/L-sin^2(pi N/L)/(N L sin^2(pi/L)),
    P_sea V_occ P_sea = direct_sum_(r,ell)
                         alpha_(L,N) Pi_r H15 Pi_r.       (8)

No equality of individual rooted contact expectations was used. This is
the same proof as [doping_check](../../contrib/counterexample_campaign/doping_check/README.md),
valid for arbitrary oriented Hermitian h15, including its large penalty.

[Proved] At fixed 0<rho<1, alpha(rho)=rho-sin^2(pi rho)/(pi^2 rho)>0.
In particular alpha(1/2)=1/2-2/pi^2. As rho approaches 1 it approaches 1,
and at small rho it is pi^2 rho^3/3+O(rho^5), consistent with rare contacts.
The leading material scale is J alpha, not J rho or J rho^2.

[Proved limitation] P_sea is a collection of seas with O(t/L) kinetic
splittings, not one exactly degenerate eigenspace. Its complement contains
charge excitations with vanishing gaps. Equation (8) cannot be used as a
volume-uniform Schrieffer-Wolff projection at fixed J.

## 4. A replacement for the frozen sector proof

[Proved, using C75's all-length inequality] In the frozen h15 chain,
fix marker positions and assign each marker to the data interval on its
left. A complete length-l cell has the full space

    K_l=(C^12)^(tensor (l-1)) tensor C^3,
    dim K_l=3*12^(l-1).                                   (9)

All off-grammar words are included. The switch on marker-to-data pairs
acts identically on each marker state and belongs to the cell on its
right. The only genuine coupling between cells is lambda V, of norm

    B=sqrt(2)|lambda|.                                    (10)

All remaining terms are cell terms. In particular, the grammar penalty
across a left marker boundary is zero. This gives an exact tensor-product
cell decomposition, without using conservation of individual bad pairs.

[Proved, inherited input] With epsilon_K=1/(K_pen-4), C75 Section 8 gives
the complete zero-port cell energy lower bound

    E_l^K >= E_l^inf-(l-2)_+ epsilon_K,
    E_l^inf=mu_sel-2cos(pi/(l+1)).                         (11)

Adding the nonnegative field preserves this lower bound. The legal empty
data vacuum remains a trial state with energy E_l^inf. The marker-free
ring has energy density at least -epsilon_K. C75 certifies

    e_b^inf=E_b^inf/b < -2^-8,
    E_l^inf/l-e_b^inf > 2^-16 for every integer l!=b.      (12)

Its monotone selector lemma supplies the all-length statement from the
two neighboring bounds. We use that established result, not a finite scan.
Put gamma=2^-16-epsilon_K>2^-17 at (3)-(4). Equations (11)-(12) imply

    E_l^K >= l e_b^inf+l gamma, l!=b.                    (13)

[Proved] Let A_m be the ground energy of an open chain of m full b-cells,
with the exact field and port, and u its infinite-chain energy per cell.
Concatenating q copies adds at most q-1 bonds of norm B. The variational
principle gives

    q A_m-(q-1)B <= A_(qm) <= q A_m+(q-1)B,
    |A_m-m*u|<=B, S_m<=m*u+2B, u<=E_b^inf.               (14)

The last bound follows from a product of empty-data legal cell vacua,
whose port expectation vanishes. The almost-additivity bounds establish
the limit defining u.

[Proved] On a compatible material ring N=mb with at least two cells,
suppose k>=1 cells have wrong lengths l_1,...,l_k and let
D_wrong=sum_i l_i. Remove at most 2k ports incident to them. The selected
cells form at most k open runs. Apply (13)-(14) to get

    E_wrong >= (u/b)N+gamma D_wrong-3kB,
    E_wrong-S_m >= gamma D_wrong-3kB-2B
                  >= (gamma-5B)D_wrong > 0,             (15)

provided 5B<gamma. This holds for the displayed lambda schedule. The last
step uses D_wrong>=k>=1. Arbitrarily long wrong cells pay a linear cost;
their shrinking internal walker gaps never enter a denominator.

[Proved] Marker-free rings cannot compete. Their energy is at least
-epsilon_K N, while a selected trial has energy N e_b^inf. All other
marker sectors are covered by (15). Thus the full frozen model selects
length 16 at nonzero lambda on compatible rings. This improves the
zero-port statement in C75, conditional only on its proved all-length
bound, which has not received an independent referee here.

[Proved distinction] Illegal words inside the selected full cells are
virtual components of their eigenstates. They are not separate conserved
sectors. For example RR R can move to R L R with amplitude -1, changing
the bad-pair count. The d24 proof based on commuting illegal-bond
projectors therefore does not transfer verbatim. Equation (15) replaces
it for the frozen model. It does not prove a mobile defect inequality.

### Thermal corollary for the frozen model

[Proved, using the same inputs] The strategy of the decisive campaign's
[defect-pressure theorem](../../contrib/decisive_attack_20261028/frozen/DEFECT_PRESSURE.md)
can also be repaired using full cells. Let f_cell(T) be the free energy
per selected full cell. Trace almost-additivity gives the analog of (14),
with |log Z_m^open+beta m f_cell|<=beta B. The variational principle gives
f_cell/b<=e_b^inf. Each wrong cell then has normalized weight at most

    w_l <= 15^l exp[-beta(l gamma-3B)]
         <= 15^l exp[-beta l(gamma-3B)].                  (16)

There are at most k selected runs after k wrong cells, so b^k is a safe
bound on their offsets. Sum lengths and possible starting positions to
obtain (1+W)^N, with W=b sum_(l!=b) w_l. At low enough T,

    W<=30b exp[-(gamma-3B)/T].                            (17)

Indeed q=15 exp[-(gamma-3B)/T]<=1/2 gives sum_(l>=1)q^l<=2q.
This overcount contains all marker sectors with wrong lengths. A completely
selected sector contributes only a subextensive prefactor. The marker-free
trace is bounded after normalization by
[12 exp(-(2^-8-epsilon_K)/T)]^N, which vanishes in the pressure limit.
No classical bad-bond sector sum is needed. Existence of the bounded-bond
thermodynamic limit and the compatible-ring selected-sector inclusion give

    0<=f_cell(T)/b-f15(T)
       <=480 T exp[-(gamma-3B)/T], b=16.                 (18)

All cell excitations, including illegal configurations, are already inside
f_cell. This result transfers an assumed frozen selected-cell thermal
coefficient. It does not establish that coefficient and does not transfer
it through the vacancy bath.

## 5. The vacancy fluid and its two field normalizations

[Proved, kinematics] Let g_j>=0 be the number of vacancies after material
particle j. Then sum g_j=L-N and the contact is A_j=1[g_j=0]. At zero
physical momentum, retaining a material root, the free bath is

    H_gap=-t sum_j(E_j^+ E_(j+1)^-+h.c.),
    E^-|g>=|g-1> for g>0, E^-|0>=0, E^+|g>=|g+1>.       (19)

There are no square-root occupation factors. Nonzero momentum and cyclic
root identifications require the usual boundary gluing. An unrestricted
rooted tensor product should not be mistaken for extra physical states.

[Proved thermodynamics; assumed Gaussian field identification] Put
r_gap=(1-rho)/rho. From the free fermion energy and coordinate change,

    e_gap(r_gap)=-(2t/pi)(1+r_gap)sin[pi/(1+r_gap)],
    v_gap=2t rho sin(pi rho),
    d^2 e_gap/dr_gap^2=2pi t rho^3 sin(pi rho).

In the convention Delta(exp(2i phi_gap))=K_gap, the Gaussian relation
e_gap''=pi v_gap/K_gap gives K_gap=rho^-2. At half filling

    K_gap=4, v_gap=t;
    K_particle=1, v_particle=2t.                         (20)

The distinction is a coordinate and compactification change, not an
interaction correction. K_pen is a third, unrelated parameter. All these
reference formulas are independent of the occupied alphabet and survive
the replacement 23 -> 15 without alteration.

[Heuristic] In the local phase-model contact expansion the density term
has dimension one and the first gap harmonic has dimension K_gap=4.
The zero-frequency noncontact density response cancels because

    integral_0^infinity (x^2-v^2 tau^2)/(x^2+v^2 tau^2)^2 d tau=0
    for x!=0.

Dimension-two descendants give at worst |x|^-3 static tails and the first
harmonic gives |x|^-7. These suggest a finite first spatial moment for the
local static kernel. On a fixed-N ring the forward kernel also has a
rank-one canonical contribution proportional to delta_ij-1/N. Its uniform
part must be separated before taking an absolute first moment. A change
on k bounded source sites changes (sum s_j)^2/N by O(k), so it can still
be paid for per defect. This is the d24 mechanism, unchanged at the free
reference, not a proof for the operator-valued finite-J bond.

## 6. Harmonics, thresholds and the conditional 16-mode reference

[Assumed] Below a positive dressed length-defect cost the material order
has b=16 occupied sites per cell. Let the particle counting field be
N(x)=rho x-phi_p(x)/pi, with Delta(exp(2ia phi_p))=a^2 at the reference.
Then y(x)=N(x)/b counts cells. The marker field is phi_p/b.

[Proved within this dictionary] A payload vertex with cell momentum q
and dimension Delta_pay has all integer harmonic branches

    Q_h=q+2pi ell,
    k_phys=rho Q_h/b mod 2pi,
    Delta_ref=Delta_pay+[Q_h/(2pi b)]^2.                  (21)

This follows by expanding the cell comb in exp[2pi i ell y] and inserting
y(x). Fractional marker harmonics are retained, not discarded. The first
marker harmonic has dimension 1/256 and correlation exponent 1/128.
It is an oscillating observable with wavevector pi/16 at half filling.

[Proved within the dictionary] A uniform Hamiltonian term obeys
rho Q_h/b=2pi h, h integer. At rho=1/2 this implies

    Q_h=4pi b h, q=0 mod 2pi,
    Delta_ref,res=Delta_pay+4h^2.                         (22)

For h!=0 the dimension is at least four. For h=0 the inherited 30-root and
465-pair certificate excludes every nonforward neutral quadratic or
quartic term. Forward terms are Gaussian. Higher neutral fermion vertices
have reference dimension at least three. Exact Q forbids anomalous
payload pairing and exact N forbids number-changing charge vertices.
Thus every resonant non-Gaussian term has reference dimension at least
three, with no rational-independence claim about higher combinations.

[Proved from the inherited margins and (21)] For dimension below three,
|Q_h|/(2pi b)<sqrt(3). Nonzero reciprocal vectors are separated by at least
pi(2-sqrt(3)); nonforward internal momenta by at least delta_int/(2b),
where delta_int=min(delta_2,delta_4). Hence

    delta_phys >= min(delta_int/32,pi(2-sqrt(3)))
                 >3.57625e-9.                           (23)

The free-charge fixed-momentum threshold also transfers. For any
number-preserving particle-hole excitation,

    E_exc=2t(sum_particles |cos k|+sum_holes |cos k|)
          >=(4t/pi) dist(P,pi Z).                        (24)

To see this, choose the nearest edge +/-pi/2 to each momentum. Its energy
is at least (4t/pi) times the distance to that edge, and the signed edge
sum is a multiple of pi. The triangle inequality gives (24). A twist
transfer Q_h/b requires compensating physical momentum Q_h/(2b) modulo
pi, so the limiting threshold for the certified nonforward set is

    omega_q >=2t delta_int/(pi b)>4.5534e-9 t.             (25)

The finite-ring correction at L=2N is -2t tan(pi/(2L)). Bounded rounding
of N changes it by O(t/L), as proved in the cited charge-threshold note.
This is a free-charge threshold, not a gap for all charge oscillators.

[Proved Gaussian criterion, inherited] For the full unequal-velocity
energy matrix G_ref, require

    eta=||G_ref^(-1/2)(G-G_ref)G_ref^(-1/2)||<1/4.

The dimension-metric inequality in the d24 infrared audit gives
Delta_G>=sqrt((1-eta)/(1+eta)) Delta_ref. At eta=1/4 the factor is
sqrt(3/5), so dimension-three terms stay strictly above two. This includes
charge-payload forward mixing and compares with the slowest velocity.

[Heuristic, conditional] Fifteen surviving payload Dirac pairs and the
charge fluid give sixteen compact Gaussian oscillators. In physical
coordinates the leading payload speeds are (b/rho)J alpha v_a^cell,
and the charge speed is 2t at half occupancy. Unequal speeds are allowed
in the proposed collection of conformal sectors. Mode counting is not
by itself a microscopic conformal-limit proof.

## 7. Compatible rings, rounding and defects

[Proved kinematics] With no length mismatch, N=bM. For any physical L,
choose the nearest N in b Z to L/2, so |N-L/2|<=8. The free charge energy
has second derivative 2pi t at half filling, giving a rounding cost
pi t (N-L/2)^2/L+O(t/L^2), with an additional O(t/L) twist cost. The
O(t b^2/L) conclusion survives a smooth interacting compressibility.

[Proved conditional on dressed defect stability] A positive fixed
defect cost of order J eventually exceeds that rounding cost. Thus an
incompatible physical L need not force a wrong material cell. This applies
when the ground state can minimize over N. If N is externally fixed at
an incompatible value, the obstruction has not been removed.

[Proved kinematics] For payload cell momentum P, the one-particle material
translations have branches theta=(P+2pi r)/b, r=0,...,15. These are the
charge boundary twists, not sixteen independent frozen physical vacua.
No oscillator is removed by the cyclic gluing. This does not prove a
unique ground state on every finite physical ring.

[Proved limitation] The decisive campaign's all-ring obstruction used
conserved marker and illegal-bond patterns. Its illegal-bond version is
not literally applicable here. Conserved frozen marker patterns can
still have large translation orbits on incompatible lengths. The mobile
kinematics permits compatible occupied counts; it does not prove that
the full h16 ground state chooses them without the next assumption.

[Assumed, required mobile comparison] After common bulk subtraction one
needs a one-sided bound, for a measure D of changed material length and
interfaces, of the form

    E_def(J)-E_sel(J)
       >=J alpha c_def D - r_h(J/t)D - C_h t/L,
    c_def>0, r_h(s)/J ->0 at fixed h15.                   (26)

Long wrong cells must be weighted by their length or handled with their
linear bare cost in (15). A two-sided correction proportional only to
the number of arbitrarily long wrong cells is not justified. A simple
operator norm O(J^2 L) also cannot establish (26). The dynamic contact
operator multiplies the entire h15, including K_pen, and is not a scalar
impurity. Static phase-model susceptibility tests do not prove (26).

## 8. Large finite penalty and the ultraviolet normal form

[Proved application of the supplied local lemma] In rooted squeezed
coordinates the phase shifts E^+/- and contacts are bounded operators,
even though each gap space is infinite dimensional. For a fixed h15,
M_h=K_pen+2^16 bounds its local norm at the schedule above. The
[single-exponential normal form](../../contrib/direct_attack_20261028/vacancy/SINGLE_EXPONENTIAL.md)
therefore applies with constants depending on this fixed M_h. It gives

    U H U*=t H0+J V_low+R,
    ||R||_1 <=C_h t s^2 delta^-10
              +C_h J delta^-5 exp(-c_h/delta^2), s=J/t,  (27)

with C_h s delta^-3 small. The norm carries one spatial moment. Choose
delta=s^(1/20). Then ||R||_1/J tends to zero as s ->0. Translation
covariance of this normal form is preserved. Constants are not uniform
as K_pen tends to infinity. One must fix r,u before taking v large.

[Proved limitation] V_low filters energy transfers, not absolute bath
energies. Equation (27) does not bound its inhomogeneous ground-state
energy, does not establish (26), and cannot eliminate the charge boson.
The same limitation is stated in the direct-attack source itself.

[Assumed sufficient final schedule] Let Delta_* be a positive massive-cell
gap, v_*>0 the smallest payload speed in cell units, and c_def>0 a material
margin, all for the already fixed h15. Given a fixed small tolerance eta,
choose an integer v so s=2^-v satisfies

    s M_h < omega_q/(100t),
    C_h s M_h^2(1+|log s|)
       <eta alpha min(Delta_*,c_def,v_*,v_* delta_phys).   (28)

For definiteness take eta=1/8 and choose the smallest integer v satisfying
these inequalities and the smallness condition of (27), after choosing a
valid fixed upper bound C_h. This is an explicit conditional selection
rule for the last dyadic parameter, not a numerical value for that unknown
constant. The same bound controls a pessimistic mixed-current coefficient of size
O(J M_h), whose squared ratio to charge and payload stiffnesses scales
as s M_h^2/(alpha v_*). The local ultraviolet conditions in (27) add
another finite upper bound on s. Every inequality has positive solutions
if the matching assumptions supply finite C_h. No resolved integer v
is certified because that low-frequency comparison remains unproved.

[Proved finite-cell facts; remaining representation work] At zero field,
the finite-K cell ground code is G_K|a>, with a common dressed orientation
vector and dimension 2^15, rather than the sine-wave path alone. Its full
space is 3*12^15, rather than the legal (5n+3)2^n space. C75's Schur bound
protects its positive auxiliary gap. A full-cell heat-trace bound is finite,
Tr_Q exp(-beta A)<= (3*12^15-2^15)exp(-beta Delta_*).
Thus fixed-cell spectral and state-counting arguments are available, but
the d24 matching constants and its legal-cell representation cannot be
copied numerically or verbatim.

[Proved charge distinction; assumed matching repair] The d24
[N23 matching application](../../contrib/decisive_attack_20261028/frozen/N23_MATCHING.md)
uses a microscopic Z2 grading that restricts to logical fermion parity.
Here Q=2N_log on the code, so logical parity is exp(i pi Q/2) there.
On odd-Q excited states this operator has eigenvalues +/-i and is not
a Z2 grading. The microscopic parity exp(i pi Q) is instead trivial on
the code. Consequently that particular full-space fermionic embedding
requires a new grading or an intermediate elimination of the half-charge
auxiliaries. C75 proves the quadratic CAR structure after projection;
it does not by itself verify every hypothesis of the N23 full-space
matching theorem. We retain the needed strong-norm massive matching as
an explicit assumption. This observation is not a refutation of a
finite-cell reduction by another representation.

## 9. Complete dependency list for c=16

| Input | Status here | Role |
| :-- | :-- | :-- |
| Literal sixteen-state bond, N,Q,M, oriented squeezing | Proved; numerical implementation checked | Defines the microscopic family |
| C75 n=15 finite-g, finite-K root and pair certificate | Proved, inherited computer-assisted result | Fifteen transverse reference bands |
| C75 all-length lower bound and fixed-cell gap | Proved, inherited; not independently refereed here | Controls full finite-K cells |
| Nonzero-port frozen selection and frozen pressure transfer | Proved here from that bound, (15),(18) | Repairs the noncommuting-grammar issue |
| Charge sea projection and free-charge threshold | Proved | Fixes the first-order model and momentum costs |
| K_gap=4 and K_particle=1 | Proved reference thermodynamics plus assumed Gaussian dictionary | Correct impurity normalization |
| Finite-K massive matching in suitable strong coefficient/source norms | Assumed for this 15-letter application | Controls nonquadratic and retarded payload errors |
| Spatially weighted ultraviolet normal form | Proved supplied lemma, applicable at fixed finite K | Only the high-frequency remainder |
| Low-frequency vacancy matching, preserving harmonic lattice and charges | Assumed | Retains charge and all payload modes with small normalized corrections |
| Uniform mobile dressed-defect bound (26) | Assumed | Locks N=16M and rules out competing material patterns |
| Half-density tuning with positive compressibility | Assumed after matching | Makes the resonance and rounding arguments apply |
| Harmonic resonance exclusion and Gaussian dimension criterion | Proved conditional on the dictionary and matching | Protects the proposed reference phase |
| Multichannel RG including moving Fermi points and arbitrary retarded vertices | Assumed; continuum mechanism heuristic | Connects weak endpoints to a stable Gaussian limit |
| Stress/current algebra or scalar thermal identification sufficient for sixteen conformal modes | Assumed | Converts that limit into the claimed central charge |

[Proved limitation] The direct-attack initialization, degree-rescaling,
counterterm, normalization and marked-tree lemmas are conditional on their
stated endpoint and beta-function estimates. They do not supply all those
estimates for h16. The decisive campaign likewise distinguishes source
response from a conformal limit. Normalized two-point responses alone
cannot fix the source-independent scalar free energy or its T^2 term.
None of these notes closes the last two rows of the table for C81.

[Heuristic, conditional conclusion] If the assumed inputs hold for one
fixed member of (4), with mu tuned as above, the predicted low-energy
theory has c=16>15 at d=16. This is a smaller conditional candidate than
d24. It is not an established microscopic counterexample.

## 10. Finite checks and reproducibility

[Numerical] `bond16.py` constructs (1) independently from the move table.
`checks.py` compares it against C75 after permuting the occupied basis.
The entry discrepancy, Hermiticity discrepancy and all three additive
charge commutators are zero at the resolved fixture. Its h16 has 437
nonzero entries. The port norm is sqrt(2).

[Numerical] Fixture parameters are n=1, D=4, Delta=3, g=.07, t_d=.003,
K_pen=30, mu_sel=1, f=.017, lambda=.13, with the n=1 phases of Section 2.
The bond-only charge check uses t=1.2,J=.013,mu=.037. Ring checks use t=1,
mu=0 and the J values below. These are diagnostic parameters, not (3).

[Numerical] On L=4,N=2 the entire 1350-dimensional Hilbert space is used.
Literal and independently squeezed hopping/interaction matrices agree
entrywise. All cyclic sectors and both half-flux seas give a 330-dimensional
P_sea. Equation (8) has maximum entry error 1.83e-14, with alpha=1/4.
Dense diagonalization retains all 120 levels originating in the untwisted
lowest sea, avoiding the missed-multiplicity issue reported in C75.

| J | Maximum error from projected first order | Error/J^2 |
| :-- | :-- | :-- |
| .0008 | 3.90450664e-5 | 61.0079163 |
| .0004 | 9.78131648e-6 | 61.1332280 |
| .0002 | 2.44783848e-6 | 61.1959621 |

Eigenpair residuals are below 1.20e-14. This checks finite-ring perturbation
theory and the lift, not a bound uniform in L.

[Numerical] Two n=1 cells with N=4, M=2, Q=2 have 464 material words when
both alternating marker origins and every allowed local letter are kept.
The exact invariant physical sectors have dimensions 2320 at L=5 and
6960 at L=6. Sparse ED at J=.04,.02,.01 gives six low eigenpairs in each.
All residuals are below 1.49e-11; literal/squeezed discrepancies are below
1.12e-16. At J=.01 the ground energies are -1.9998927219148006 and
-3.4640305942519087 respectively. These sectors contain illegal words and
mobile vacancies, but are not the entire physical Hilbert space. They do
not test global selection or fifteen payload modes. Full spectra and
parameter records are in `checks.json`.

[Numerical] `auxiliary_checks.py` independently diagonalizes (19) at
(L,N)=(8,4),(10,5),(12,6), dimensions 35,126,462. Free energy and contact
formulas agree to 1.2e-14 and 7.8e-16. A unit tagged contact has
chi_00=.0774768757,.0819044448,.0848338254. Its energy shift obeys
J alpha-J^2 chi_00+O_L(J^3), as checked at J=.04,.02,.01. The corresponding
relaxation/J^2 at J=.01 is .0773609501,.0817775696,.0846988344. These
scalar checks do not prove the quantum-source defect bound (26).

[Numerical] The same script compares full frozen two-cell sectors at
N=4, K_pen=1000, f=0 and the b=2 midpoint selector mu_sel=1-sqrt(2).
Both the selected lengths (2,2) and wrong lengths (1,3) have dimension
1296, without a charge or grammar truncation. Wrong-minus-selected
ground energies are .5855366805,.5855377041,.5855407747 for lambda=0,
.005,.01. Eigenpair residuals are below 1.7e-12. This is a small check
of the full-cell replacement, not the n=15 certificate.

Run from this folder with Python, NumPy and SciPy.

    python3 -B checks.py
    python3 -B auxiliary_checks.py

The scripts set one BLAS thread. Read-only imports from C75 suppress bytecode
writes. All new programs and outputs are in work_C81. No cluster run,
network command, external edit, subagent or referee was used.
