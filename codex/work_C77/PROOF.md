# C77 derivations

2026-10-28. Tags distinguish proved algebra, assumed continuum input,
heuristic search proposals, and numerical checks. Set the lattice spacing,
hbar, and Boltzmann constant to one. The central charge c is the usual
nonchiral value, so a compact boson has c=1 and Ising has c=1/2.
Every counted conformal sector must have a positive limiting velocity.

## 1. Hard-core particles, holes, and the mode count

[proved, local algebra] The site basis is |0>,|1>,...,|q>, with 0 a hole.
Write X^{ab}=|a><b| and n_a=X^{aa}. Then

    n_0 + sum_(a=1)^q n_a = I,
    n_a n_b = delta_ab n_a.

The hole density is not an independent density. A variation obeys

    delta rho_0 = -sum_(a=1)^q delta rho_a.

There are q independent density coordinates, not q+1. One particle per
site does not imply sum_a delta rho_a=0 when holes are mobile. Only at
strict unit filling is the total particle-density coordinate frozen.
Nor does fixing the global particle numbers remove local sound modes.

[assumed, complete ordinary liquid] Suppose the entire low-energy theory
is obtained by linearizing one occupied sea for each of q species, with
no other soft sector. There are q right movers and q left movers.
Bosonization supplies q canonical pairs (phi_a,theta_a). A positive
quadratic Hamiltonian has the form

    H_G = (1/2pi) int dx [ (partial theta)^T V_theta (partial theta)
                          + (partial phi)^T V_phi (partial phi) ],

with possible mixed gradient terms when parity is broken. A nonsingular
canonical change of basis diagonalizes the stable quadratic dynamics.
Each massless canonical pair contributes c=1. Forward scattering changes
stiffnesses and velocities without adding a pair. Relevant pinning terms
can remove modes. Thus within this complete liquid description

    c_total <= q = d-1,                                           (1)

with equality if all q modes remain massless. This is a physics-level
bound in a specified continuum regime. It is not an algebraic theorem
about every strongly interacting Hamiltonian on these site states.

[proved within that description] Define

    phi_c = q^(-1/2) sum_a phi_a,
    phi_s^alpha = sum_a M_(alpha a) phi_a,   alpha=1,...,q-1,

where the rows of M are orthonormal and orthogonal to (1,...,1).
Transform theta by the same orthogonal matrix. This is an invertible
change of q canonical pairs, not an enlargement. Spin-charge separation
therefore gives one charge mode plus q-1 relative-species modes. At an
SU(q)-symmetric level-one spin fixed point,

    c_spin = dim SU(q)/(1+h_dual) = (q^2-1)/(q+1) = q-1,
    c_charge = 1,
    c_total = q.                                                 (2)

For q=1 the spin sector is absent. In a charge-gapped state only the
q-1 relative modes can survive within this description. Spin gaps or
depleted species reduce the count further. A hole mode added to (2)
would count the charge coordinate twice. The familiar spin-charge
current decomposition is explicitly derived in Assaraf et al.,
[cond-mat/9903057, Eqs. (19)-(26)](https://arxiv.org/pdf/cond-mat/9903057).
That paper's unprojected Hubbard site has dimension 2^q. Its continuum
algebra is an input here, not a demonstration that projection preserves
finite spin velocities.

### 1.1 Fermi seas and the meaning of the projection

[proved, bare one-body count] For an independent range-one hopping
channel with complex amplitude t_a,

    epsilon_a(k) = -2 |t_a| cos(k+arg t_a) - mu_a.

A regular partial filling has two simple Fermi points and one connected
occupied arc on the Brillouin circle. The two endpoints form one Dirac
sector with c=1. They are not two nonchiral sectors. A phase shifts the
arc without splitting it. A band-edge double zero has zero velocity and
is a z=2 threshold, not another CFT. An empty or full band has no regular
Fermi points.

[proved, scope] In the hard-core Hilbert space X^{a0} is not an independent
canonical fermion for each a. For example,

    {X^{0a},X^{a0}} = X^{00}+X^{aa},

which is not I for q>1. Thus the q-sea argument is a continuum hypothesis,
not a diagonalization of the projected many-body hopping. An unprojected
q-flavor fermion Fock site has 2^q states, whereas the hard-core site has
q+1. The no-double-occupancy constraint neither provides an extra hole
sea nor proves that all q ordinary seas actually survive dynamically.

### 1.2 Literal hopping plus density interactions has frozen color order

[proved] The literal Hamiltonian specified in the question is

    H_hc = -sum_(j,a) [t_a X_j^{a0} X_(j+1)^{0a} + h.c.]
           + sum_(j,a,b) V_ab n_(j,a) n_(j+1,b)
           - sum_(j,a) mu_a n_(j,a).                              (3)

Terms involving n_0 can be eliminated using completeness. V_ab is real
and need not equal V_ba. A basis state on an open chain can be labeled by
ordered positions x_1<...<x_N and the ordered color word
w=(a_1,...,a_N). Every allowed hop moves a particle into an adjacent hole.
No two particles can pass. A density term is diagonal. Consequently

    <x',w'|H_hc|x,w> = 0 for w' != w.                             (4)

The corresponding word projectors commute with H_hc. On a ring only the
cyclic word is invariant, since a particle can cross the chosen cut.
Cyclic rotation is not arbitrary local spin exchange. No sequence of the
allowed moves exchanges two adjacent material colors, so perturbation
theory in (3) cannot generate an off-diagonal exchange between word sectors.

[proved, symmetric example] If t_a=t, mu_a=mu, and V_ab=V independently of
color, the open-chain restriction to every word is the same spinless
hard-core hopping and density-interaction Hamiltonian. At fixed N,

    H_hc |_(N) is unitarily equivalent to H_charge |_(N) tensor I_(q^N).
                                                                    (5)

This is a decomposition in ordered coordinates, not a sitewise tensor
factorization. It proves that the dynamical spectrum is one spinless
spectrum repeated q^N times. In its ordinary liquid regime the charge
sector has c=1. The colors have zero exchange energy, not q-1 finite-speed
CFT modes. Summing all words in the fixed-density thermal ensemble gives
an entropy term (N/L) log q even at T=0. The full model therefore is not
an ordinary finite-degeneracy CFT plus q-1 propagating spin sectors.
Kumar's [infinite-U solution](https://arxiv.org/html/0812.2977v1) provides
the q=2 free-charge example and its spin degeneracy explicitly.

[assumed, regular exchange extension] Adding

    H_ex = J sum_j sum_(a,b=1)^q X_j^{ab} X_(j+1)^{ba}              (6)

allows adjacent occupied particles to exchange colors. Antiferromagnetic
exchange can give the usual SU(q)_1 spin sector, alongside the charge
liquid, when the filling and couplings place the system in that phase.
This is a t-J-type model. A large but finite Hubbard repulsion generates
exchange of order t^2/U, but strict U=infinity sets it to zero. The limits
T->0 at fixed J>0 and J->0 before T->0 do not commute. Also, finite-U
Hubbard is a 2^q-state microscopic chain, not a q+1-state chain.

[assumed, established exact-model continuum; proved bond decomposition]
An explicit saturation benchmark in the exchange-extended class is the
antiferromagnetic permutation chain

    h=J P,   J>0,
    P=sum_(a,b=0)^q X^{ab} tensor X^{ba}.

Its a=0 or b=0 off-diagonal terms hop particles into holes. Its a,b>0
off-diagonal terms exchange colors. Its diagonal terms are densities.
The standard SU(q+1)_1 infrared identification gives

    c=[(q+1)^2-1]/(q+2)=q.

The SU(N) saturation is recorded in
[Latorre and Sierra](https://arxiv.org/html/2403.17242v2).
This benchmark has d=q+1 and one fixed nearest-neighbor bond. It belongs
to the enlarged class including (6), not to (3) for q>1.

[proved limitation; heuristic assessment] Equation (4) holds also for
species-dependent hopping and interactions, but (5) does not. The charge
Hamiltonian then depends on the frozen word. Ordinary homogeneous fluids
still supply no mechanism for an extra independent spin CFT. However,
(4) alone is not an exhaustive classification of every critical point of
these position Hamiltonians. It does not prove a global maximum c=1 or
c=q for all choices in (3). An assertion of such a maximum would require
an additional dynamical theorem. The defensible maximum q in (1) is the
standard complete q-mode count, saturated when color exchange is allowed.
The task's literal density-only wording must not be silently enlarged.

## 2. Counting conserved densities

[proved, linear algebra] A family of commuting Hermitian one-site charges
can be simultaneously diagonalized. The diagonal matrices form a
real d-dimensional vector space. Quotienting the identity leaves at most
d-1 linearly independent charges. Hence

    r_onsite <= d-1.                                             (7)

This concerns independent commuting additive on-site charges. It does not
bound integrable charges of arbitrary support, the energy density, or
all noncommuting symmetry generators. Independent densities must also be
distinguished from their composite operators.

[assumed CFT input; proved implication] If r independent nontrivial
abelian currents survive in a unitary CFT, their abelian Sugawara tensor
has central charge r. If there is a commuting neutral tensor,

    T = T_U(1)^r + T_neutral,
    c = r + c_neutral,   c_neutral >= 0.                          (8)

Thus current existence gives a lower bound c>=r, not the upper bound
c<=r. To obtain c<=d-1 from (7), one must additionally establish that all
propagating sectors are Gaussian and are resolved by those charges.
Emergent currents and neutral sectors invalidate that inference. Even a
critical Ising chain has c=1/2 with no continuous on-site charge. Rank of
a nonabelian symmetry is not a central-charge bound either, since
SU(2)_k has rank one and c=3k/(k+2). Lattice representation and anomaly
constraints may restrict a particular realization, but rank alone does
not establish those restrictions.

[proved, conditional refinement] If there is a nontrivial U(1) current
and the remaining compact unitary Virasoro sector has 0<c_neutral<1,
unitarity permits c_neutral=1-6/[m(m+1)], m=3,4,... . Its first positive
value is 1/2. In this setting a violation of the qubit bound must satisfy
c>=3/2. This is a bound on the first possible excess once a charged
sector is established, not a proof that qubits cannot realize it.
Without the current hypothesis this argument does not apply. These CFT
inputs and their scope were already recorded in C63 Sections 3-4.

## 3. Counting Fermi points beyond independent bands

[proved, quadratic number-conserving class] Let there be r canonical
fermion orbitals per primitive cell, with a range-one Bloch matrix

    H(k)=A+B exp(ik)+B^dagger exp(-ik).

The determinant det[H(k)-mu] is a Laurent polynomial with powers from
-r to r. Multiplication by z^r, z=exp(ik), gives a polynomial of degree at
most 2r. Unless identically zero, there are at most 2r simple Fermi
crossings, counted with multiplicity. Each paired right/left crossing
contributes c=1, so c<=r. The determinant statement is safer than claiming
that every mixed band separately crosses only twice. An identically
zero determinant signals a flat band at mu, whose zero-velocity extensive
degeneracy is not an additional regular CFT. The regular remaining
channels are handled by the corresponding reduced problem.

[assumed, inherited free theorem] The project C1 free-fermion analysis
also controls pairing and gives c<=r for regular range-one quasi-free
chains. Particle-hole redundant BdG roots must not be counted as
independent fermions. The unconstrained local Fock dimension is d=2^r,
so this class obeys the stronger c<=log_2 d<=d-1.

[proved limitation] A general interacting chain has no Bloch determinant
whose degree bounds all low-energy poles or collective sectors. Products
of neighboring bond operators already have longer support. For example,

    (X_1 X_2)(Y_2 Y_3)=i X_1 Z_2 Y_3.

Eliminating virtual states can generate effective longer-range motion and
multi-site terms even with a range-one microscopic bond. A block of ell
sites has dimension d^ell and may carry many emergent modes. Merely
folding a fixed free dispersion into an enlarged cell does not create
modes, but genuinely interacting collective dynamics is not subject to
the bare scalar cosine count. Nor is a path through lattice couplings an
RG flow to which the continuum c-theorem can automatically be applied.

## 4. Entropy, velocities, and the missing crossover scale

[proved, Hilbert-space bound] On L sites, S(rho)<=log(dim H)=L log d.
The thermodynamic entropy per site therefore obeys s(T)<=log d. For an
ordinary critical system with finitely many ground sectors, s(0)=0.
Extensive frozen-word degeneracy, as in (5), violates the ordinary CFT
premise and must be kept separate.

[assumed, conformal thermodynamics] Let a label a nonchiral sector with
central charge c_a and positive right/left speeds v_(a,R),v_(a,L). Then

    f(T)=e_0 - (pi T^2/12) sum_a c_a[1/v_(a,R)+1/v_(a,L)] + o(T^2),
    s(T)= (pi T/6) sum_a c_a[1/v_(a,R)+1/v_(a,L)] + o(T).           (9)

The entropy expansion is assumed in its own right, or follows with
appropriate differentiability control of the free-energy remainder.
For v_(a,R)=v_(a,L)=v_a this becomes s=(pi T/3) sum_a c_a/v_a+o(T).
The formula pi/6 times a single sum in the task is correct if that sum
runs over chiral branches. In that convention the sum of branch central
charges is 2c_total, where c_total=sum_a c_a. A factor of two cannot be
changed without also changing the central-charge convention.

[assumed, standard locality input] For a finite-range bounded bond there
is a Lieb-Robinson speed v_LR, of order the bond norm in lattice units,
bounding propagation of locally visible modes. It is an upper speed,
not a positive lower speed. See
[Nachtergaele and Sims](https://arxiv.org/abs/math-ph/0506030).
Applying this to all conformal branches gives

    A := (pi/6) sum_a c_a[1/v_(a,R)+1/v_(a,L)]
       >= pi c_total/(3v_LR).                                  (10)

[proved conditional inequality] Suppose a temperature T_0 lies in a
common conformal window, with a specified fractional error epsilon<1,

    |s(T_0)-A T_0| <= epsilon A T_0.

Combining this with (10) and the Hilbert-space bound gives

    (1-epsilon) pi c_total T_0/(3v_LR) <= s(T_0) <= log d,
    c_total <= 3 v_LR log d/[pi(1-epsilon)T_0].                   (11)

For physical spacing a_0, the entropy per length bound is log d/a_0,
and the last expression contains v_LR/(a_0 T_0). Thus the dimensions
match. A common rescaling of all energies scales v_LR and T_0 together
and cannot improve (11).

[proved limitation] The asymptotic expansion supplies some sufficiently
small T_0 for each model and each epsilon. It supplies no uniform lower
bound on T_0/v_LR across the class. Letting T_0->0 removes any upper
bound on c. Small velocities strengthen the entropy constraint at a
fixed temperature, but can also lower the scale where their linear
spectrum applies. Neither a velocity upper bound alone nor a lower bound
on velocities alone supplies the missing common conformal window.
A postulate T_0>=eta v_LR with eta>0 would yield an O(log d) bound, but
eta is precisely the additional scale information absent here.

[proved, logical countertest for the counting inequalities] For arbitrary
positive c and v, set

    A=pi c/(3v),   tau=log d/A,
    s_trial(T)=log d (T/tau)/sqrt[1+(T/tau)^2].                    (12)

Then s_trial=A T+O(T^3), 0<=s_trial<=log d, and

    ds_trial/dT=(log d/tau)[1+(T/tau)^2]^(-3/2)>0.

The specific heat T ds_trial/dT is positive. The compatible free energy
is e_0-log d tau(sqrt[1+(T/tau)^2]-1). It has the correct leading
low-temperature form and high-temperature entropy limit. Taking
v<=v_LR satisfies the velocity input for any c. The crossover tau
shrinks as 1/c. These are thermodynamically consistent test functions,
not claimed lattice realizations. They prove that these entropy,
positivity, and speed inequalities alone cannot exclude arbitrary c.

[heuristic, spatial version] If conformal block entropy is valid only
beyond ell_*, the bound S(ell)<=ell log d does not justify substituting
ell=1 in S(ell)=(c/3)log ell+constant. Under controlled errors at ell_*
and 2ell_*, their entropy difference is (c/3)log 2, whereas the elementary
Hilbert-space estimate gives at best O(ell_* log d). The unknown crossover
length remains. Large emergent cells are costly in scale, not excluded
by the on-site entropy bound.

## 5. Smallest dimension and a concrete search family

[proved logical conclusion from the audited inputs] At d=1 there is no
nontrivial Hilbert space or positive central charge. The counting
arguments above, and C63/C67/C68, do not exclude every nearest-neighbor
qubit Hamiltonian. Therefore the smallest dimension not excluded by
these arguments is d=2. This does not assert existence there. In
particular, the supplied background sentence that nothing exceeds c=1
at d=2 reports the absence of an established example, not a global
exclusion theorem. C63 Section 6, C67 Section 6, and C68 Sections 6-8
explicitly leave the missing dynamical classification open.

[heuristic, concrete candidate family] A short candidate bond is

    h(gamma,u,D,b) = (1+gamma) X tensor X + (1-gamma) Y tensor Y
                    + u Z tensor Z + D(X tensor Y-Y tensor X)
                    + (b/2)(Z tensor I+I tensor Z).              (13)

It has d=2, one-site translation invariance, Hermiticity, and range one.
For gamma D b !=0 it has parity product_j Z_j but no continuous additive
on-site spin symmetry. C67 identifies this same family as the unresolved
strict parity stratum, and already rejects several finite-size excess-c
signals in it. It is a candidate search space, not a newly discovered
high-c point.

[assumed, inherited continuum control] The point

    (gamma,u,D,b)=(0,0.4,0.5,0)

is C68's chiral XXZ anchor, with one boson and positive speed. Its
incommensurate weak-perturbation analysis gives controlled one-boson
neighborhoods. Their known count is c=1, so small perturbation theory
about the anchor cannot manufacture an extra sector. A possible violation
requires continuation beyond that complete weak-coupling neighborhood.

[heuristic, target] Continue (13) at nonzero gamma,D,b and ask whether a
separate neutral Ising mass can vanish while the boson remains regular.
The target is c=1+1/2=3/2>1. In a proposed boson-Ising description, enumerate
all local invariant relevant products before claiming that one mass
parameter is sufficient. For a boson vertex V with dimension x_V,

    x(sigma V)=1/8+x_V,   x(epsilon V)=1+x_V.

Their relevance thresholds are x_V<15/8 and x_V<1. A parity-odd sigma
has no uniform magnetic source, but a product with a parity-odd boson
vertex can be invariant. Translation momenta can forbid or allow that
product independently. Sigma times a current has dimension 9/8 when
allowed, and epsilon times a current is marginal. These perturbations
can require additional tunings or destroy the candidate product theory.
No existence or full-rank matching map is inferred from having four
parameters. General qubit exchange matrices supply further directions
if a branch is found, but they supply no existence proof either.

[assumed, success criteria] A successful point needs a converged entropy
coefficient and independent spectral identification of both sectors,
with finite positive velocities and controlled relevant/marginal mixed
couplings. With unequal speeds, the periodic Casimir term is

    E_0(L)=e_inf L - pi/(12L) sum_a c_a[v_(a,R)+v_(a,L)] + o(1/L),

subject to the correct vacuum and zero-mode sector. It measures a
velocity-weighted sum, not c_total divided by an arbitrary fitted speed.
Two correlation peaks can be harmonics of one boson. A zero-velocity
threshold or extensive ground degeneracy is not a second CFT. C63's
factorization boundary already gives an entropy-size slope approaching
3/2 with z=2 excitations, and C67 rejects a short-ring excess through
L=20. These inherited failures are reasons to require the two-sector
check, not proofs that (13) cannot succeed elsewhere.

## 6. Checks performed in C77

[numerical, finite algebra only] Run `python3 checks.py`. For the symmetric
literal hard-core model at q=2,L=6,N=3 with open boundaries and V=0.23,
the full 160-dimensional spectrum equals the 20-dimensional spinless
charge spectrum repeated eight times, to 1.20e-14. Every hopping matrix
element preserves the ordered color word. This checks (4)-(5).

[numerical, finite algebra only] At
(gamma,u,D,b)=(0.17,0.4,0.5,0.23), the bond (13) is Hermitian and commutes
with Z tensor Z, with zero numerical residual. The uniform additive-charge
commutator map has three positive singular values, approximately 5.2341,
3.7868,1.9233. This checks the absence of a uniform continuous spin charge
at that test point. It does not assign the point an infrared theory.

[numerical, formula check] Equation (12) passes the entropy bound,
monotonicity, and low-temperature slope tests at formal c=1,1.5,10,100,
d=2,v=1 on 1000 temperatures from 10^(-6)tau to 10^6 tau. No lattice
realization is inferred. No large local computation, cluster job, or
thermodynamic central-charge measurement was performed in C77.
