# C17. Mechanisms for c > 1 in a nearest-neighbor qubit chain

2026-10-27. Analytic audit, with small consistency checks. No counterexample
and no proof of the full c <= 1 claim were obtained.

The strongest exclusions are for exact microscopic constructions. A conserved
one-site qubit marker makes the chain classical. Literal interleaved subchains
cannot propagate with range one. Exact U(1) gives the XXZ family, including its
bound strings, with at most one conformal mode under the supplied XXZ input.
The free-fermion and zero-field XYZ exclusions remain as recorded in C6.

The remaining gap is dynamical. No argument here excludes a strongly interacting
neutral sector, an emergent higher-level current algebra, or a quantum-selected
cell with a nontrivial low-energy code. Symmetry counting and relevance counting
constrain such possibilities but do not prove they cannot occur. There is no
known exhaustion theorem that turns a list of familiar mechanisms into a
classification of every possible qubit infrared theory.

## 1. Scope, tags and evidence

We use the background's Hamiltonian and normalization,

    H = sum_i [sum_ab J_ab sigma_i^a sigma_(i+1)^b + B.sigma_i],

with one fixed real J, one fixed B, and periodic boundary conditions. Local
constants and telescoping one-site terms have been removed. A Majorana CFT has
c=1/2, a compact boson has c=1, and simultaneous gapless sectors contribute
their sum even when their velocities differ. A direct sum of alternative
vacua does not contribute the sum of their central charges.

Tags used throughout:

- [proved] An algebraic argument is given here, or a specific repository proof
  is cited. A conditional result identifies its hypotheses explicitly.
- [assumed] Standard continuum or integrable-model input supplied by the
  repository, or an explicitly stated effective description. No new literature
  verification was performed, in accordance with the no-network instruction.
- [heuristic] A physical expectation, candidate exclusion, or proposed test
  that does not establish the conclusion for the entire microscopic class.
- [numerical] A computation actually performed here, with sizes and estimators.
  Proposed computations are labeled as tests, not numerical evidence.

The principal sources are [STATUS.md](../../STATUS.md), [REPORT sections 3.4
and 8](../../REPORT.md), [C6](../work_C6_d2/c6.tex), [C7](../work_C7/C7_result.md),
[C8's qubit witness](../work_C8/qubit_ff_witness.md), [exact encoding
obstruction](../work_C8/exact_encoding.md), [boundary-rank
obstruction](../work_C8/boundary_rank_obstruction.md), and [C5's parton
analysis](../work_C5_F1/f1_results.md). This report does not promote their
conditional infrared statements to unconditional theorems.

### Mechanism summary

| Mechanism | Decision at d=2 | Status and limitation |
|---|---|---|
| Longest frustration-free segment | Cannot select a finite length above three | [proved], K_FF(2)=4, open chain without boundary penalties |
| C7-style conserved one-site marker and walker | Impossible | [proved], marker conservation forces a binary classical chain |
| Classical binary period p>2 | Cannot be a strictly preferred bulk period | [proved], minimum-cycle argument; degenerate longer words are possible |
| General quantum spontaneous period p | Not excluded | [heuristic/open], neither selection nor a useful payload has been constructed |
| Microscopic SU(2) with k>=2 | Excluded for this NN qubit bond | [proved reduction; assumed Heisenberg IR] |
| Emergent SU(2)_k, SO(3)_1 or a neutral product | Not excluded | Level ceiling requires extra assumptions; anomalies do not upper-bound c |
| Atom and cluster fluids with exact U(1) | No c>1 in the supplied XXZ classification | [proved reduction; assumed XXZ IR] |
| Atom and cluster fluids without U(1) | Generic conversions often remove modes, but no universal exclusion | [assumed bosonization; heuristic reachability] |
| Completed Abelian parton coset | c<=1 | [proved conditional count; assumed identification of physical IR] |
| General paired/Majorana partons | Not excluded by the count | Hard projection is not automatically a continuum coset |
| Literal even/odd or residue-class chains | Impossible with a nonzero NN connected bond | [proved] |
| Emergent or spontaneous decoupling | Not excluded in general | Free and certain exact-code versions are excluded |
| A tuned product CFT, e.g. Ising times boson | Not excluded | Relevant couplings need tuning or selection rules |
| Coexisting vacua, rotor towers, z!=1, nonunitary fits | Do not establish an in-scope c>1 CFT | [proved distinction; assumed CFT diagnostics] |

## 2. Spontaneous enlargement of the unit cell

### 2.1 Blocking is not a construction

[proved] Regrouping p qubits produces a block Hilbert space of dimension 2^p,
but the microscopic Hamiltonian still has its original constraints. In a fixed
blocking it contains p-1 internal bonds and one interblock bond, all inherited
from the same h. They cannot be chosen independently. Translation by one site
also remains a symmetry of H, although it is no longer on-site in block
coordinates. An arbitrary chain on C^(2^p) is not thereby realized.

[assumed] In a phase with spontaneous period p, choose one extremal vacuum.
There are p translated vacua if no further degeneracy is present. Their finite
ring superpositions can give a cat-state entropy contribution of at most log p.
For fixed p this is a constant, not another logarithmic slope. The CFT charge
is that of the fluctuations within a vacuum. Domain walls between vacua can
be massive. A p-fold ground multiplet is not evidence for p massless modes.

[heuristic/open] A sector with broken translation and propagating critical
degrees of freedom is worth testing even if one adopts a stricter definition
of a single CFT on every ring. This report does not exclude that sector by
definition. On incompatible rings its domain-wall energy and twisted sectors
must be separated from a Casimir term. A mobile-cell proposal must additionally
show how its length fluctuations and payload coexist at finite coupling.

### 2.2 What the frustration-free witness excludes

[proved, repository] For a uniform positive two-qubit bond with no boundary
terms, a nonzero zero-energy space on four open sites implies one at every
length. Thus K_FF(2)=4. A finite longest zero-energy segment has at most three
sites. We use the proof in C8 without rederiving it.

[proved, scope] This does not bound the period of a frustrated ground state,
the length selected by a nonzero ground energy, or the degeneracy of a
frustration-free manifold. In particular it does not exclude exponentially
many allowed configurations, later split by a kinetic perturbation. C16 is
analyzing that last possibility. Small nonzero energies are not exact zeros,
and boundary penalties change the theorem's hypotheses.

Test. For this specific selector, compute the nullities of open H_N at
N=2,3,4. A nonzero nullity at N=4 ends a claimed finite-longest-FF construction.
It does not settle a kinetic construction. This is a proposed test; the
repository already supplies the witness calculation.

### 2.3 A conserved one-site qubit marker cannot host C7's walker

[proved] Suppose M_i is a translated nontrivial rank-one one-site projector
and [H,M_i]=0 for every i. A global basis rotation makes M_i=(1+Z_i)/2.
Commuting with every Z_i says H is diagonal in the common computational basis.
Equivalently, expand H in Pauli strings. Distinct centered two-site supports
and single-site supports are linearly independent on rings of length at least
four. Every coefficient containing X or Y must vanish. Therefore

    J_ab = 0 unless a=b=z,       B_x=B_y=0,
    H = sum_i [J_zz Z_i Z_(i+1) + B_z Z_i] + constant.

Telescoping terms do not alter the conclusion about H. The same argument
applies if only one [H,M_i]=0 is initially specified, since translation gives
the others.

[proved] A diagonal binary NN model can be assigned directed edge costs
w_ab, a,b in {0,1}. Any cyclic word is a closed walk on a graph with two
vertices. Decompose its edges into simple directed cycles. Its mean energy
is a convex combination of their mean energies. The only simple cycles are
0, 1, and 01. Hence its minimal possible bulk energy density is

    min{w_00, w_11, (w_01+w_10)/2}.

There is always a minimizing period-one or period-two pattern. Longer
minimizing words occur at ties but are not strictly selected periods. An
odd ring can require a defect in a period-two phase, without changing this
thermodynamic statement. The graph argument is about orthogonal binary
letters, not arbitrary noncollinear coherent-spin product ansatzes.

[proved, application] C7's d=4 selector conserves M_i=|#><#|_i and uses the
three-dimensional complement {L,X,R} for the move LX <-> XR. On a segment,
the walker has energy -2t cos(pi/b), so E_b=mu-2t cos(pi/b) can select long b.
At d=2 a rank-one marker has a one-dimensional complement. Conservation
therefore eliminates the walk, as the stronger Pauli argument shows. A
literal qubit implementation of that frozen-marker mechanism is impossible.

[heuristic/open] A marker encoded over several sites, a mobile marker, or a
marker conserved only within an emergent low-energy subspace escapes this
proof. The effective interaction can depend nonlinearly on interval length.
There is no proof here that such dependence cannot stabilize a long quantum
cell. Saying that a binary alphabet has only two letters is insufficient
once the legal configurations and dynamics are quantum.

[numerical] Four-site linear algebra using all twelve J,B coefficients
leaves exactly the two directions J_zz,B_z when [H,Z_0]=0 is imposed. The
null-space projector error is 4.45e-16. Forty random directed binary cost
matrices, tested on all words at lengths 2 through 10, obey the minimum-cycle
bound to 4.45e-16. These check the algebra, not arbitrary kinetic selectors.

Test for an alternative selector. Demonstrate a period-p ground manifold,
finite wall and internal-cell gaps at fixed p, and an extensive low-energy
payload before fitting c. Measure bond and spin structure factors at 2pi/p,
the near-ground momenta 2pi j/p, and gaps on both compatible and incompatible
rings. Use a weak pinning field to select a vacuum, then extrapolate it away.
A peak that moves with L or an exponentially large zero-space is not this
demonstration. No such candidate was tested in this run.

### 2.4 Period two, period three, and modes carried by a cell

[proved] A period-two block can contain a spin-one triplet, and a period-three
block can contain several doublets or a spin-3/2 multiplet. These are
representation-theoretic possibilities. They do not provide an independently
selectable block Hamiltonian or arbitrary interactions between blocks.

[assumed] In the elementary strong-dimer construction, two ferromagnetically
bound qubits give a triplet. If S is its spin-one operator, projection gives

    P sigma_left^a P = P sigma_right^a P = S^a.

Thus a weak isotropic inter-dimer exchange projects at first order to a
bilinear spin-one exchange. It does not directly give the
Takhtajan-Babujian bilinear-biquadratic critical point with c=3/2. The
biquadratic term arises at higher order and is suppressed in that controlled
limit. Moreover, independently strong and weak bonds explicitly double the
cell and are not allowed here. A spontaneous construction would have to
derive both the hierarchy and the desired effective couplings from one h.

[assumed] An antiferromagnetic Heisenberg trimer with a separated lowest
doublet yields one effective spin-1/2 per trimer under weak exchange. This
familiar truncation gives one spin mode, not several. Keeping two low doublets
would change the retained theory and requires its own degeneracy and coupling
analysis. Neither this example nor the dimer example classifies all codes.

[proved, narrow first-order obstruction] Consider a two-qubit block intended
to carry three on-site Clifford generators Gamma_a and a parity P anticommuting
with all three, with the two ports realizing Gamma_a and i Gamma_a P by
compression. The four anticommuting involutions require a code dimension
divisible by four. A two-qubit block then has to retain its entire space.
Compression is a unitary identification, so its disjoint endpoint algebras
commute. The desired maps do not, since

    [i Gamma_a P, Gamma_a] = -2i P.

This excludes that three-Majorana first-order target for a two-site cell.
It does not exclude a spin-one target, higher-order encodings, or larger cells.

[proved, warning] The connected bond rank rank J<=3 does not by itself prove
c<=1 for a general block code. Three real Majorana channels would have c=3/2,
and higher-order effective couplings can have a larger operator span. Abstract
port compressions and a local parent realizing them are different problems.
The C8 exact-product-code theorem excludes the parity-dressed hopping target
when the full tensor product code is exactly invariant with separated ports.
It expressly leaves perturbatively dressed codes open.

[assumed, examples] If a period-two order parameter alone becomes critical in
the ordinary Ising universality class, it contributes 1/2. If period-three
order has an ordinary three-state Potts transition, it contributes 4/5. Those
are specific universality assumptions, not consequences of the period. A
commensurate-incommensurate floating phase often has one phason boson, c=1.
The finite number of translated vacua never counts these modes by itself.
An additional payload or a coincident order/disorder transition could in
principle add a neutral sector. That is one of the unresolved mechanisms.

[heuristic/open, answer for general p] There is no established implication
"period-p order in a range-one qubit chain => at most one gapless mode".
It holds in the supplied quadratic class even with a period, and in the
controlled examples just described. It has not been established for an
arbitrary interacting, spontaneously enlarged cell.

Test. For p=2,3, measure the cell multiplet and its endpoint matrix elements
on open segments, then compare the retained-chain spectrum with the microscopic
one. A credible payload requires at least two independently dispersing towers
with gaps proportional to 1/L, entropy giving their summed c, and a consistent
velocity-weighted Casimir coefficient. Compare sizes divisible by p and, when
needed, by the effective chain's additional commensurability. This is a proposed
targeted test, not evidence that such a cell exists.

## 3. Emergent non-Abelian currents and anomalies

### 3.1 Exact microscopic SU(2) is already closed

[proved] Two spin-1/2 sites decompose as 0 direct-sum 1 without multiplicity.
An SU(2)-invariant bond is a scalar on each summand, hence

    h = a I + b Swap = a' I + b' sigma_i.sigma_(i+1).

There is no invariant field and no independent biquadratic coupling for
two qubits. For H, the same reduction holds after removing telescoping terms.

[assumed, supplied Heisenberg IR] Antiferromagnetic exchange gives SU(2)_1,
c=1. Ferromagnetic exchange has quadratic magnons and a large ground multiplet,
not a finite-velocity CFT. Zero exchange is not a propagating CFT. This closes
microscopic SU(2), without needing a conjectural level bound.

### 3.2 What the embedding index does and does not say

[proved, conditional] With the fundamental qubit representation R of SU(2),
the embedding index in SU(2)_1 is one. If the physical currents embed in that
parton current algebra and the entire stress tensor is exhausted by it, coset
unitarity gives c<=1. More generally the argument in REPORT section 3.4 assumes
k<=ind(R), then uses

    c[SU(2)_k] = 3k/(k+2),      ind(fundamental)=1.

At k=1,2,3 the charges are 1,3/2,9/5. The ceiling k<=1 is the substantive
dynamical hypothesis. It does not follow from the dimension of the single-site
operator algebra. A continuum current can be a scaling limit of operators on
many sites, and one microscopic operator can contain several scaling fields.

[proved, logical scope] Without a microscopic SU(2) action there is no
microscopic representation R to which the proposed emergent SU(2) currents
must be assigned. One cannot impose the fundamental embedding index on them
by labeling the ultraviolet qubit a spin-1/2. A product of emergent current
algebras is equally outside that inference. Even when one current algebra is
identified, Sugawara gives c_total=c_current+c_coset with c_coset>=0. This is
not an upper bound on a neutral coset.

[assumed, normalization] SO(3)_1 in the three-Majorana convention is the same
local affine theory as SU(2)_2, with c=3/2. Global quotients can alter the
allowed primary sectors but do not supply a new central-charge exclusion.
The level labels should not be compared without this normalization.

### 3.3 LSM matching constrains representations, not the size of c

[assumed, standard anomaly input] A half-integer spin per unbroken unit cell
has the mixed SO(3)-translation LSM obstruction. For a pure SU(2)_k WZW theory
with the usual antiferromagnetic translation action g -> -g, its matching
invariant is k modulo two. In that setting odd k matches a qubit cell and even
k does not. This permits k=3,5,...; it never implies k<=1.

[proved given that input] The qualifier about the translation action and the
rest of the infrared theory is necessary. A further sector or broken
translation can carry the missing anomaly. After period p is selected, a
blocked cell has projective class p modulo two. For even p the original
one-cell anomaly no longer obstructs a symmetric gapped state under the
surviving translation. For odd p it remains. None of these statements limits
the number of simultaneous massless sectors.

[assumed, surviving-symmetry cases]

- With U(1), flux insertion constrains filling and the momentum shift,
  Delta P=2pi nu modulo 2pi. In a unique, gapped, period-p charge plateau it
  requires p nu integer. It does not bound the number of neutral modes.
- A single unitary Z2 on a qubit has no nontrivial projective representation
  obstruction of this type, since its generator can be rephased to square to
  one. Two pi spin rotations forming Z2 x Z2 can act projectively and restore
  an LSM obstruction, but still no upper bound on c.
- Physical spin time reversal squares to -1 on a qubit, giving a Kramers-per-cell
  LSM obstruction with translation. The plane-reflection antiunitary in C6
  squares to +1 and has no such on-site Kramers obstruction. Neither produces
  a level ceiling.
- With no nontrivial internal symmetry, translation alone admits a unique
  gapped product ground state. There is no spin or filling anomaly to match.
  Emergent symmetries need not be restrictions of microscopic on-site ones.

[heuristic] Breaking the microscopic symmetry usually permits relevant
anisotropies and masses at a putative high-level fixed point. For three
Majoranas, each diagonal mass i xi_R^a xi_L^a has dimension one. With fewer
symmetries, off-diagonal masses and other local fields can also be allowed.
Their momenta and locality in the physical spin theory must be determined.
This makes an emergent SU(2)_2 point less robust, but an unstable fixed point
still lies within the question. The eight-parameter bond family can in
principle tune several relevant couplings. Parameter counting alone does not
prove it intersects, or avoids, the fixed point's critical manifold.

Test. Search for the full affine spectrum, not just an entropy value 3/2.
For a diagonal SU(2)_k theory the primary of spin j has
x_j=2j(j+1)/(k+2), j=0,1/2,...,k/2, with symmetry and momentum assignments
to be checked in the actual modular invariant. Identify dimension-one currents
and their multiplets, measure a common current-sector velocity, and compare
entropy and Casimir estimates. Perturb along symmetry-breaking bond directions
to count which masses open. If microscopic spin is not conserved, a spin twist
or uniform-spin susceptibility is not automatically a measurement of the
emergent current level. No such candidate is identified here.

## 4. Multiple fluids, atoms, pairs and longer strings

### 4.1 Exact U(1) exhausts the bond, not just a convenient subfamily

[proved] Rotate the conserved axis to z. Invariance under simultaneous
rotations about z forces

    J = [[a,d,0],[-d,a,0],[0,0,b]],       B=(0,0,g).

The transverse exchange has one complex hopping amplitude. A rotation whose
phase advances uniformly along the chain removes its phase in the bulk,
leaving XXZ with transverse coupling sqrt(a^2+d^2). The ring has a boundary
twist. It changes finite-size zero-mode sectors, not the local bulk CFT.
This rotation does not remove generic DM terms in a transverse field or
with transverse pairing anisotropy, since those acquire spatial phases.

[proved] In hard-core particle notation, after this gauge choice and a
possible staggered sign change, the complete U(1) class has the form

    H = -t sum_i (b_i^dagger b_(i+1) + h.c.)
        + V sum_i n_i n_(i+1) - mu sum_i n_i,       t>=0,
    n_i in {0,1}.

With Pauli coefficients above, t=2 sqrt(a^2+d^2), V=4b, and the chemical
potential includes the longitudinal exchange's linear term and g. Thus
Delta=V/(2t) is the conventional XXZ anisotropy after the staggered rotation.
There is no independently tunable pair hopping or longer-range repulsion
within exact U(1) and one physical bond.

[assumed, supplied XXZ input] Its finite-velocity gapless phases have one
compact boson. Other regimes are gapped, ordered, fully polarized, or at
nonconformal endpoints. This includes the magnetized regimes and the effects
of Bethe strings. Consequently exact U(1) closes the multi-fluid proposal at
the physics standard requested here. Counting the one microscopic charge
alone would not close it, since one charge can coexist with a neutral sector.

### 4.2 Two-body binding and why it does not imply a pair fluid

[proved] Set the hopping phase to zero and first omit -mu N. In the two-body
sector use total momentum K and separation r>=1, with wavefunction f_r. Define
C=2t cos(K/2), choosing K in [-pi,pi]. The lattice Schrödinger equations are

    E f_r = -C (f_(r-1)+f_(r+1)),                  r>=2,
    E f_1 = V f_1-C f_2.

For f_r=lambda^(r-1), the bulk and boundary equations give

    E=-C(lambda+lambda^-1)=V-C lambda,
    lambda=-C/V,
    E_pair(K)=V+C^2/V.

A normalizable attractive bound state requires V<0 and |V|>|C|. The scattering
band at that K has lower edge -2|C|. Thus a pair can bind at large |K| even for
weak attraction. This is only a fixed-momentum statement.

[proved] The global two-atom threshold is -4t. A pair minimum at K=0 is
normalizable only for V<-2t. Writing u=|V|, its energy is

    E_pair(0)=-u-4t^2/u < -4t,      u>2t.

For u<2t the allowed momenta have C<u and their pair energies obey
E_pair>-2u>-4t. Such bound states do not beat two atoms at their common
single-particle minimum in the dilute ground-state problem. At u=2t the
K=0 state is at threshold and has lambda=1. Restoring chemical potential
subtracts 2mu from both competing energies and changes none of these tests.
A hopping phase merely shifts the pair momentum by twice the atom shift;
it does not tune their minima independently.

[assumed, XXZ input supported by the following limit] V<-2t is also the
attractive phase-separation regime. In the limit |V|>>t, r adjacent particles
have interaction energy (r-1)V. Joining two clusters gains another negative
bond energy V; it does not stabilize a preferred finite cluster size. At
fixed density the low-energy state aggregates, rather than producing stable
independent atom, pair and trimer gases. The exact threshold and the absence
of an intervening multi-fluid ground phase use XXZ input, not the two-body
equation alone. Repulsively bound pairs are excited states and do not change
that ground-state conclusion. The t=0 limit is diagonal and nonpropagating.

[numerical] The relative-coordinate Hamiltonian truncated at r=200 gives
the following lowest bound energies, with t=1 and mu=0:

| V | K/pi | E_pair | Global two-atom threshold | Error against formula |
|---|---:|---:|---:|---:|
| -3 | 0 | -4.333333333333 | -4 | below 1e-14 |
| -1.5 | 1/2 | -2.833333333332 | -4 | 1.1e-12 |
| -1 | 4/5 | -1.381966011250 | -4 | 2.3e-16 |

The last two are bound relative to their momentum-resolved continua and
still fail the global dilute-fluid test.

Test for longer strings. Compare E_L(N+1)+E_L(N-1)-2E_L(N), the low-energy
minima in the few-particle momentum sectors, and the convex envelope of
E_L(N) versus N. Convex-envelope linear segments and macroscopic clusters
signal coexistence. To claim two fluids one needs two finite-velocity towers
and two independently varying low-energy densities, not merely several
bound-state branches above the ground state. C15 is pursuing the detailed
few-particle and Z2 analysis; no additional scan is performed here.

### 4.3 Conversion operators with an assumed two-fluid infrared theory

[assumed, bosonization setup] Suppose, independently of microscopic
reachability, there are atom and pair liquids with phases theta_a,theta_p,
density fields phi_a,phi_p, and positive stiffnesses. For independent liquids
choose conventions

    H_0 = sum_s v_s/(2pi) integral dx
          [K_s (partial_x theta_s)^2 + K_s^-1 (partial_x phi_s)^2],
    rho_s(x)=rho_s-(partial_x phi_s)/pi+harmonics,
    dim exp(i l theta_s+2i m phi_s)=l^2/(4K_s)+m^2 K_s.

Here the annihilation field's leading bosonic harmonic is
a_s(x) ~ exp(i Q_s x-i theta_s). Fields and harmonics must be local physical
operators in the proposed theory. A hard-core atom cannot be squared on one
site, but its adjacent-site pair is nonzero, so hard-core exclusion does not
in general remove the continuum two-atom conversion vertex.

[proved given this setup] Local pair-to-two-atom conversion gives

    p^dagger a a + h.c.
      ~ cos(theta_p-2theta_a+(2Q_a-Q_p)x),
    Delta_conv = 1/K_a + 1/(4K_p).

It preserves total charge N_a+2N_p, so exact microscopic U(1) does not forbid
it. Translation permits the leading term only if 2Q_a-Q_p=0 modulo 2pi.
Density harmonics add momentum 2pi sum_s m_s rho_s, with integer m_s and
signs fixed by the chosen harmonic. They also add sum_s m_s^2 K_s to its
dimension in the independent-fluid limit. For an already selected period-p
vacuum the reciprocal cell momentum is 2pi/p, making more resonances possible.
One must include the order parameter that supplies that momentum in the
translation-symmetric description.

[assumed, standard RG] A nonoscillating scalar vertex of dimension Delta<2
is relevant, Delta>2 is irrelevant at weak coupling, and Delta=2 requires
its beta function. A relevant conversion locks one relative phase and
usually leaves one charge boson. For K_a=K_p=1 the dimension is 5/4, so this
common dilute-liquid estimate excludes a stable two-boson fixed point with
resonant nonzero conversion. It does not exclude a tuned zero of its
coefficient or a neutral critical point at a phase-locking transition.

[proved given Gaussian input] With density interactions the two K values
must be replaced by a matrix. In a Gaussian theory without density-current
mixing, define the positive K matrix by the logarithmic covariances of phi
and theta. Then

    Delta(l,m) = (1/4) l^T K^-1 l + m^T K m,
    l=(-2,1) for atom-pair conversion.

Off-diagonal entries change the relevance condition. More general chiral
quadratic terms require the corresponding full covariance matrix. A
microscopic DM term is not permission to use the independent-fluid formula
without checking this structure. Derivative density couplings have dimension
two and can renormalize the Gaussian matrix without gapping a mode. Density
locking, for example cos[2(phi_a-phi_p)], has dimension K_a+K_p in the diagonal
case and is allowed only at its corresponding density/momentum resonance.

[proved given the same setup] Conversion of r atoms into a charge-r cluster
has leading dimension

    Delta_r = r^2/(4K_a)+1/(4K_r),
    momentum = r Q_a-Q_r modulo 2pi.

Between charge-r and charge-s clusters the minimal conversion uses s/g
clusters of the first and r/g of the second, where g=gcd(r,s). Its phase
vector has those integer coefficients. More generally, every integer phase
vector l satisfying sum_r r l_r=0 is charge-neutral. Charge counting is only
one condition. Locality, statistics, momentum and the full covariance matrix
are also needed. For fermionic fields, same-branch coincident products can
require derivatives; the right-left pair can be allowed without that
extra suppression. One cannot assign the bosonic leading harmonic to every
parton or soliton operator without identifying its physical local lift.

### 4.4 The available K range does not give a universal relevance proof

[assumed, XXZ input] At half filling in -1<Delta<=1, the conventional single
charge mode has

    K = pi/[2(pi-arccos Delta)].

It runs from 1/2 at Delta=1 to infinity as Delta approaches -1 from above;
K=1 at Delta=0. The endpoint with divergent K has vanishing velocity and
is not an additional finite-velocity CFT. This formula is not the formula
at all fillings or for two hypothetical emergent fluids.

[proved within a controlled projected limit; assumed LL interpretation]
For strong repulsive V/t, below density 1/2, project to particles with no
adjacent occupancies. Ordered coordinates map to y_j=x_j-(j-1), a free
hard-core gas on length L-N, up to a boundary twist. Put q=pi rho/(1-rho).
Its thermodynamic energy and physical sound velocity are

    e(rho)=-(2t/pi)(1-rho) sin q,
    v=2t sin q/(1-rho).

The squeezed Fermi velocity is 2t sin q. A low particle-hole energy spacing
2pi(2t sin q)/(L-N), expressed as 2pi v/L on the physical ring, gives the
stated v. Differentiating twice gives e''(rho)=2pi t sin q/(1-rho)^3. The Gaussian
compressibility identity e''=pi v/K therefore gives

    K=(1-rho)^2 -> 1/4 as rho -> 1/2 from below.

The same value is understood as the charge-1/2 soliton threshold of the
period-two plateau: a dilute soliton gas has K_s=1, while the physical
density field is phi=phi_s/2 and its conjugate phase theta=2theta_s, giving
K=K_s/4. The threshold itself has vanishing velocity. Taking large but finite
V and an interior density gives K<1/2 with a nonzero velocity, so the failure
of a universal K>=1/2 claim is not confined to a singular endpoint. C15 gives
the more detailed strong-repulsion analysis.

[proved, consequence] Even the independent-fluid arithmetic does not force
conversion relevance on the basis of a claimed lower bound K>=1/2. If
K_a=K_p=1/2 its dimension would already be 5/2, not less than two. At
K_a=K_p=1/4 it would be five. These values are kinematic illustrations,
not realizations of two independent low-K liquids by a qubit NN bond.
Conversely the existence of a large K in the attractive XXZ liquid does not
produce an independently stable pair fluid. The exact U(1) conclusion rests
on the full model classification, not on a universal scaling inequality.

[heuristic/open] There is no established range for every eigenvalue and
orientation of the hypothetical multi-fluid K matrix in the generic
symmetry-free qubit problem. Supplying one would require new microscopic
control. Using the known single-mode XXZ K range twice is not justified.

### 4.5 Removing U(1) adds locking terms but does not complete the proof

[proved given operator identification] Without U(1), charge-changing
vertices must be considered as well. In the two-fluid notation, examples are

| Candidate vertex | Dimension for independent bosonic fluids | Leading momentum condition |
|---|---:|---|
| a+h.c. | 1/(4K_a) | Q_a=0 mod 2pi |
| aa+h.c. | 1/K_a | 2Q_a=0 mod 2pi |
| p+h.c. | 1/(4K_p) | Q_p=0 mod 2pi |
| p^dagger aa+h.c. | 1/K_a+1/(4K_p) | 2Q_a-Q_p=0 mod 2pi |

If the microscopic Z2 is atom-number parity, the first is forbidden and the
other three are even. If no internal symmetry survives, all can be allowed
when they are local and satisfy momentum selection. A uniform microscopic
source cannot supply arbitrary nonzero momentum. Additional density harmonics
and their dimensions must be checked. Pairing or a transverse field may also
invalidate the initial identification of independently conserved fluids.

[heuristic] Resonant low-order source and conversion operators commonly pin
phases and destroy extra bosons. This is a reason to expect c<=1 near the
controlled one-mode liquids. It is not a proof at an arbitrary strongly
interacting fixed point. Incommensurate momenta can eliminate low-order
resonances, and isolated tunings can set relevant coefficients to zero.
The NN U(1) chain cannot independently detune atom and pair minima, as shown
above. The corresponding obstruction has not been proved for all dressed
excitations of a bond without U(1).

[assumed continuum possibility; heuristic microscopic status] A locked
atom-pair theory can have an Ising transition of a residual relative sign
while its total charge mode remains gapless. Such a transition would give
Ising times a boson, c=3/2, rather than two independent conserved fluids.
Source terms, allowed spin fields, and couplings to the charge mode can
destroy or change it. No allowed qubit NN realization is supplied here,
but the conversion argument alone does not exclude this tuned possibility.

Test. At a candidate point, determine wavevectors in single-spin, adjacent-pair
and longer-string correlations, their decay exponents, and momentum-resolved
low gaps. Fit a positive Gaussian covariance matrix only if the spectrum
supports that description. Enumerate all local symmetry-allowed vertices
with dimension at most two and their momentum mismatches. If a small allowed
conversion perturbation gaps a relative mode, measure its gap scaling
Delta_gap proportional to |lambda|^(1/(2-Delta_conv)), allowing for logs
and other relevant couplings. Cross-check total c by entropy and
velocity-resolved finite-size energy. These are proposed tests, not a claim
that a second fluid was found.

## 5. Fractionalization, spinon-holon separation and partons

### 5.1 Physical charges and the minimal representation

[proved] Commuting nontrivial Hermitian on-site charges on C^2 can all be
diagonalized together. Modulo the identity, their span has dimension at most
one. Consequently their susceptibility matrix has rank at most one.
If these densities resolve every Gaussian mode, there is at most one boson.
The qualifier is substantive: a neutral Ising sector has no
such density and is not counted by this rank.

[proved, kinematics] Ordinary spinful particles with an independently
variable hole and two spin states require at least the three on-site states
|0>, |up>, |down>. A qubit supplies two. In the usual spin-chain parton
representation

    |up> = f_up^dagger |vac>,
    |down> = f_down^dagger |vac>,
    n_up+n_down=1,

the charge fluctuation is constrained away. There is no independent physical
holon occupation in that representation. This excludes the literal
three-state t-J mechanism, not every possible emergent neutral degree of freedom.

[proved conditional; assumed coset identification] A number-conserving
two-flavor free parton theory with two active Dirac channels has c_free=2.
Removing its active diagonal U(1) current subtracts one. The current algebra
decomposition is

    U(2)_1 / U(1)_2 = SU(2)_1,       c=2-1=1.

The physical chain must actually realize this coset, with no further gapless
sector, for this to bound its charge. Projecting a variational state or writing
the on-site constraint does not establish that dynamical identification.

### 5.2 The complete occupation-code bound applies conditionally at d=2

[proved conditional, C5 theorem] Suppose a qubit code consists of two distinct
occupation vectors s^(0),s^(1) in {0,1}^n. Their affine span has dimension one.
There are n-1 independent linear constraints constant on that span. If every
parton supplies one active massless number current, all these constraints
are gauged, and the physical stress tensor is the resulting torus coset, then

    c = n-(n-1)=1.

Frozen orbitals are first removed. A one-state code gives zero. Thus the
REPORT section 8 bound is strong enough for the completed Abelian subclass
at d=2. Its assumptions fail for an arbitrary coherent non-Abelian code,
a paired ansatz without those currents, or an unidentified neutral sector.

[proved, explicit failure of incomplete subtraction] With four fermionic
occupations and charges (1,2,4,7), the local constraint

    n_1+2n_2+4n_3+7n_4=7

leaves only (0,0,0,1) and (1,1,1,0), a qubit. Subtracting just one free U(1)
would formally give c=4-1=3. But the density 2n_1-3n_2+n_3 is zero on both
states, even though its formal current is orthogonal to the gauged charge
vector and has nonzero free norm. That current cannot survive the hard
projection as the claimed physical density. There are three independent
linear constraints, so their complete torus subtraction gives c=1.

[numerical] Enumerating the sixteen local occupation states confirms the
two allowed states, affine rank one, and the identically vanishing density.
This reproduces the small local check, not a continuum projection theorem.

### 5.3 Majoranas and gauge fields do not supply an automatic subtraction

[proved, kinematics] Four on-site Majoranas have a four-dimensional Fock
space. Fixing their local parity leaves a qubit. A continuum theory of four
independently massless Majoranas would have c=2, but the local parity
restriction alone does not justify that continuum theory. A single
inter-site Majorana hopping flips the constrained parity on both endpoints,
so its hard projection onto fixed parity vanishes. Physical two-site spin
interactions are even at each endpoint and are quartic in those Majoranas.

[assumed, continuum gauge fact] Gauging a finite group, when the result is a
consistent finite orbifold of a CFT, preserves the local stress tensor and
its central charge. There is no continuous Sugawara charge to subtract.
Thus "four Majoranas minus a Z2 constraint" is neither a numerical c bound
nor a derivation of a c=2 qubit chain. The crucial missing step is to show
that the hard microscopic projection has that orbifold as its infrared
limit with all four channels still present.

[assumed, physical gauge description] A pure gauge field in one spatial
dimension does not supply an independent transverse photon CFT. Gauss-law
constraints and possible global flux sectors must be distinguished from
propagating matter. Confinement or Higgs mechanisms often remove parton
modes, but a universal confinement assertion would also be insufficient
here: critical gauge-invariant composites can remain gapless. Neither
"the gauge field adds c=1" nor "all fractionalization confines" is a
valid general argument for this task.

[heuristic/open] Multiple spinon minima, paired partons, and neutral
parafermion or coset sectors are alternate descriptions of the unresolved
strongly interacting mode problem. The parton representation is not a
new independent microscopic construction. A proposed ansatz must reproduce
the physical local operator algebra and the actual Hamiltonian's spectrum.
The NN mean-field hopping range is not inherited automatically by physical
exchange after projection, and conversely the free physical-fermion bound
does not directly constrain every parton ansatz.

Test. Enumerate the local code and every independent linear constraint.
Check which proposed currents vanish after compression, and whether the
purported hopping has P H_hop P=0. Then compare the projected or microscopic
chain's operator spectrum, boundary-condition sectors, entropy and Casimir
coefficient with the proposed coset/orbifold. A mean-field count of Fermi
points or a good variational energy alone does not pass this test.

## 6. Decoupled sublattices and emergent decoupling

### 6.1 Literal microscopic residue classes

[proved] Let L be even and at least four, and suppose
H=H_even tensor I_odd+I_even tensor H_odd on the physical even/odd tensor
factorization. Every NN connected Pauli term sigma_i^a sigma_(i+1)^b acts
nontrivially on both factors. Its Hilbert-Schmidt inner product with either
summand is zero because its Pauli factor on the other sublattice is traceless.
Distinct such Pauli strings are independent. Taking their inner products
with H therefore gives J_ab=0 for every a,b. Only one-site fields and a scalar
remain, so neither factor propagates.

[proved] The same argument works for residue classes modulo R>1 on compatible
rings: a range-one edge joins different residues. Deleting every other edge
by hand breaks one-site translation and leaves dimers, not two extended
chains. An on-site basis rotation does not change the supports or their
connected Schmidt rank and cannot evade the argument.

[assumed, standard realizations] With range R one can instead use

    H_R=sum_i h_XXZ(i,i+R).

This consists of R decoupled spin chains on compatible rings. Taking each in
its Luttinger phase gives c_total=R, while translation by one site permutes
the chains. Range two therefore already permits c=2. It is the missing
spatial edge, not translation invariance itself, that closes this literal
mechanism at range one.

[numerical] For L=4, projecting the twelve-parameter Hamiltonian onto operators
that act on both the even and odd factors leaves a null space of dimension
three, exactly the uniform fields. Its null-space projector agrees with that
subspace to machine precision. This checks the support argument.

### 6.2 Dual sublattices are possible, but do not imply an excess

[proved mapping; assumed Ising IR] Consider the allowed NN Hamiltonian

    H=-J sum_i X_i X_(i+1)-g sum_i Z_i Z_(i+1).

Using tau_i^z=X_i X_(i+1), Z_i=tau_(i-1)^x tau_i^x gives

    H_dual=-J sum_i tau_i^z-g sum_i tau_(i-1)^x tau_(i+1)^x.

In dual variables the even and odd chains decouple. At J=g>0 each is
critical Ising, so the total is 1/2+1/2=1. The duality includes global
constraints and twisted boundary sectors, which do not alter this bulk
central term. The original Hamiltonian is also an XY chain after a global
spin rotation. Thus a claim that range-one qubits can never have any
description as decoupled sublattices would be false. The support theorem is
about the physical tensor factors. It does not prohibit nonlocal dualities.

[proved, supplied free theorem] Merely folding a free dispersion into a
larger Brillouin-zone cell does not add low-energy states. The repository's
quadratic NN bound applies at any finite period and gives c<=log_2(2)=1.
Multiple crossings in a folded plot must be unfolded and counted with their
Majorana/Dirac multiplicities. The dual example saturates this bound.

### 6.3 What a spontaneous version would require

[proved, conditional low-energy logic] If one sublattice is gapped and used
as a virtual mediator, eliminating it can generate effective interactions
between the remaining sites at second or higher order. Those remaining
sites form one active chain. The mediator cannot also be counted as an
independent gapless chain in the same controlled elimination.

[heuristic/open] To obtain two gapless emergent factors, a proposed spontaneous
construction must retain both and show that all their dangerous cross-couplings
are absent, irrelevant, oscillatory, or tuned away. This is not literal
microscopic decoupling. Blocking, perturbative generation of longer effective
range, and fractionalization can all lead to such a proposed description;
none is excluded by the Pauli-support proof alone. A two-sector strong-coupling
fixed point is part of the same unresolved problem as the cell payload and
neutral-sector mechanisms.

Test. First compute the microscopic cross-sublattice operator component, which
must vanish for literal decoupling. For emergent decoupling, measure separate
low-energy dispersions and identify the leading continuum cross-operators.
Test their response to the corresponding allowed microscopic perturbations.
Check whether apparent two-mode behavior survives larger L after the
smallest gap and the smallest velocity are resolved. No interacting
spontaneous decoupling candidate is established here.

## 7. Other routes and false positives

### 7.1 Coincident critical sectors and isolated multicritical points

[proved arithmetic; assumed candidate field theories] Two Ising sectors give
c=1 and do not violate the claim. Three Ising sectors give 3/2, Ising times
a boson gives 3/2, two bosons give 2, and two tricritical Ising sectors give
7/5. More general neutral cosets and parafermion sectors also deserve the
same audit. These sums describe simultaneous local sectors, not alternative
ground-state branches at a first-order crossing.

[assumed, perturbative analysis of a candidate] Near Ising times a boson,
the Ising energy epsilon has dimension one and its order field sigma has
dimension 1/8. A local, zero-momentum epsilon term is a mass that must be
tuned. If sigma is a physical local, symmetry-even, zero-momentum field,
its coefficient must also vanish. It may instead be forbidden by a Z2,
translation, or locality, depending on the actual embedding. The continuum
disorder field is not automatically a second local lattice perturbation.

[proved within the normalization of section 4] Candidate mixed vertices have
dimensions, when local,

    dim[epsilon cos(2m phi+n theta)]
        =1+m^2 K+n^2/(4K),
    dim[sigma cos(2m phi+n theta)]
        =1/8+m^2 K+n^2/(4K).

Each must satisfy the momentum and symmetry constraints of both factors.
The derivative coupling epsilon partial_x phi has dimension two and can
have nontrivial marginal flow, including velocity effects. If the boson
has no exact U(1), the phase vertices of section 4 must also be included.
Thus tuning one Ising mass beside a known c=1 phase is not automatically
sufficient, but there is no universal count of tunings before its operator
embedding is known. C16 is examining these counts in greater detail.

[heuristic/open] A special point can set several relevant coefficients to
zero without being a stable phase. Generic instability is therefore not
an exclusion. Conversely, the existence of eight microscopic parameters
does not show that the desired coefficient map has a zero or full rank.
The missing data are a candidate bond, its local scaling operators and
their responses to the allowed perturbations. An isolated neutral-sector
multicritical point is the most economical remaining violation mechanism.

Test. Follow a candidate c=1 boundary, identify an additional gap closing with
independent quantum numbers or momentum, and fit the product spectrum.
Perturb in independent J,B directions to determine the masses and mixed
couplings. Require entropy consistent with the sum, two towers with finite
velocities, and their weighted Casimir coefficient. A level crossing between
two phases, each with one mode, is insufficient.

### 7.2 Additional Fermi pockets and weak-coupling arguments

[proved, supplied free theorem] A physical quadratic NN Jordan-Wigner chain
has c<=1. Multiple folded bands, accidental BdG degeneracies and left/right
asymmetry do not evade the theorem's correct zero counting. At a band edge
the velocity vanishes, so a root of higher order cannot be counted as an
extra ordinary massless CFT channel.

[assumed, conditional RG statement] A unitary relativistic continuum flow
whose complete ultraviolet theory is the known one-boson or two-Majorana
theory cannot end at larger c. Applying the c-theorem here requires that
complete continuum matching, including all modes and an appropriate
relativistic regime. A connected region of lattice parameter space is not
a renormalization-group trajectory. New soft sectors at strong coupling or
vanishing-velocity points invalidate an argument that follows c continuously
over that region. Different velocities also require care before invoking a
single Lorentz-invariant c-theorem.

[heuristic/open] An interaction-induced reconstruction into several low-energy
pockets, or a neutral branch not visible in the single-particle Green
function, remains possible in principle outside the controlled neighborhoods.
At exact U(1) the NN XXZ solution closes it. Outside U(1), a single-particle
dispersion is not a complete accounting of low-energy degrees of freedom.
This is another description of the unresolved multi-mode mechanism, not an
additional established qubit construction.

### 7.3 Constrained Hilbert spaces and unusual critical theories

[proved, scope] A large-penalty NN bond can leave a constrained low-energy
configuration space. Although the physical tensor product remains qubits,
the retained configurations need not factorize site by site. Projected
one-site and NN perturbations can consequently have conditional multi-site
actions. K_FF(2)=4 does not exclude this degeneracy, and the single-site
parton charge rank does not classify its neutral modes. C16's blockade
analysis is the relevant complementary work.

[heuristic/open] Any proposed clock, height, parafermion, or effective anyon
description must be derived from those retained configurations and their
actual projected dynamics. Familiar unitary minimal-model order transitions
have c<1, but coupling one to another gapless sector can exceed one. An anyon
chain's fusion-space dimension is not the local dimension of the original
qubit tensor product. Its nearest-neighbor notion must also be translated
back into physical range before it can count as an in-scope example.

Test. Determine the retained space and the first nonzero effective Hamiltonian,
including virtual corrections. Resolve global sectors and boundary conditions.
Then apply the same entropy and spectral diagnostics to the microscopic chain
at finite penalty. A projected critical point at infinite penalty alone does
not demonstrate a finite-coupling counterexample.

### 7.4 Nonconformal, noncompact and nonunitary limits

[proved distinction] Saturation edges, Lifshitz points and ferromagnetic
Heisenberg magnons can have gaps proportional to L^-2 rather than L^-1.
Their degeneracies and entanglement logarithms do not define the c in this
question. The fully symmetric half-filled ferromagnetic Dicke state has
block entropy growing as (1/2) log L for fixed block fraction, from the width
of its hypergeometric magnetization distribution. Fitting that slope with
the periodic CFT coefficient c/3 would suggest 3/2, even though its dynamics
is quadratic and its choice within a degenerate ground multiplet matters.
This is an analytic cautionary example, not a numerical claim from this run.

[proved algebra, supplied free dispersion] Quadratic singularities are not
the only ones. C16 identifies the free branch

    E_-(k)=sin k-sqrt[(cos k-1)^2+sin^2 k]
          =sin k-2sin(k/2)=-k^3/8+O(k^5),       k>0.

Its partner is linear near zero. This mixed linear/cubic touching is also
nonconformal despite being on the free locus. Its free neighbors remain
subject to c<=1. An interacting neighborhood needs additional analysis;
a list consisting only of ordinary quadratic saturation points is incomplete.

[assumed, compact-boson check] Letting the radius of a single compact boson
diverge does not change its central charge from one. It produces a dense
zero-mode spectrum and difficult finite-size crossovers. A genuinely
noncompact unitary CFT with additional propagating fields and c>1 would
remain in scope if it had the required finite-velocity local scaling limit.
Noncompactness alone is not an exclusion theorem. No such qubit NN
realization is identified here.

[proved, unitarity distinction] The finite-dimensional Hamiltonian is
Hermitian and has a positive physical Hilbert-space inner product. A faithful
unitary scaling limit cannot acquire negative-norm states. A nonunitary CFT
or a nonunitary transfer-matrix effective central charge does not meet the
stated hypothesis. Excited-state fits, twisted sectors and constrained
statistical ensembles can use different quantities. They must not be
reported as the unitary vacuum central charge of this H.

[proved distinction; assumed finite-size formulas] A finite number of
coexisting vacua contributes at most an O(1) entropy term. Large rotor
multiplets, a cat state, a near-first-order level crossing, a small unresolved
velocity, or a drifting incommensurate zero-mode sector can distort an entropy
or ground-energy fit. If sectors have velocities v_a, the periodic Casimir
term is

    E_0(L)=e_infinity L-(pi/6L) sum_a c_a v_a + o(L^-1)

only in the appropriate vacuum sector. Twists and nonvacuum conformal weights
add their own 1/L terms. The Casimir coefficient by itself measures neither
sum_a c_a nor a universal single velocity. A half-chain fit alone is equally
insufficient when the purported CFT window has not been demonstrated.

Test. Compare spectral gaps at several momenta and sizes to determine z and
each velocity. Check stability of entropy fits to the fit window, symmetry
sector, and compatible ring subsequence. At a putative noncompact point,
separate the zero-mode tower from propagating oscillator towers and show
whether velocities remain finite. Compare boundary twists and pinning-field
limits. These tests distinguish false positives from a new infrared theory;
they do not assume every large fitted c is an artifact.

### 7.5 A concrete companion warning about effective sublattices

[assumed companion result, not independently rederived here] C15's updated
[PROGRESS.md](../work_C15/PROGRESS.md) studies a pairing-only, strongly
repulsive corner with vanishing ordinary hopping. Virtual pairing generates
two-site hopping, with apparent minima at 0 and pi. Its reported second-order
Hamiltonian reduces to a positive graph Laplacian plus a nonnegative
separation penalty, with ferromagnetic degeneracies in one-sublattice sectors.
It also identifies an exact factorization line with a quadratic defect
branch. Thus the displayed special line is not two stable conformal fluids.
The surrounding perturbations remain unresolved in that report.

[proved logical consequence] The example illustrates why the physical
sublattice support theorem and the free single-cosine dispersion are not
general theorems about effective Hamiltonians. Virtual processes can generate
longer effective hopping already from range-one qubit interactions. To decide
whether the extra minima represent independent CFT sectors one must keep
their interactions, conversion terms and singular energy scales. This run
does not duplicate C15's few-particle calculations or count their results as
independent numerical evidence.

## 8. Reproduction and the scope of the numerical evidence

[numerical] Run from this directory:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 checks.py > checks_output.txt
```

The script uses only numpy and scipy, reads no inputs, and makes no network
calls. All computations are small local checks. Results are in
[checks_output.txt](checks_output.txt); code is in [checks.py](checks.py).
There are no random-bond phase scans, DMRG runs, or new c>1 candidate results.

[numerical] The algebra and pair checks are specified in sections 2,4,5,6.
For estimator calibration the script diagonalizes the full periodic spin
Hamiltonian on L=6,8,10,12,14 for

    XX:     H=-(1/2) sum_i (X_i X_(i+1)+Y_i Y_(i+1)),
    Ising:  H=-sum_i X_i X_(i+1)-sum_i Z_i.

[proved free-dispersion normalization; assumed standard CFT fitting forms]
Both have v=2. XX has epsilon(k)=-2 cos k with half-filled Fermi points
at +/-pi/2. Ising has quasiparticle energy 4|sin(k/2)| at its critical point.
The exact ground energies for the chosen even rings are -2 csc(pi/L) and
-2 csc(pi/(2L)), respectively. Those values are used only to check the ED.

The two fitted estimators are

    S(ell) = s_0+(c_S/3) log[(L/pi) sin(pi ell/L)],
    E_0(L)/L = e_infinity+a/L^2+b/L^4,
    c_Cas = -6a/(pi v).

At each size, the entropy fit uses every integer ell from 2 to L/2. The
Casimir fit uses all five sizes and includes the leading L^-4 correction
in the energy density. The largest matrix dimension is 16384.

| Control | Expected c | c_S at L=6 | c_S at L=14 | c_Cas, L=6..14 |
|---|---:|---:|---:|---:|
| XX | 1 | 1.077768 | 1.027712 | 0.999204 |
| Ising | 1/2 | 0.505047 | 0.502011 | 0.499976 |

[numerical] All ground-energy checks pass at absolute error below 1e-10.
The largest reported eigenvector residual is below 6e-12. The script's
assertions on the algebra, pair formula, code enumeration and calibration
all pass. The entropy overshoot in XX demonstrates that a few-percent excess
from these small rings cannot distinguish a new CFT from ordinary finite-size
corrections. It does not establish a universal error bar for other models.

[assumed, protocol for a real candidate] With several sectors, first infer
their velocities from their dispersing towers. Then compare the energy
coefficient with sum_a c_a v_a, rather than dividing by whichever gap is
smallest. If left and right velocities differ, replace c_a v_a by
c_a(v_(a,R)+v_(a,L))/2 for a nonchiral sector with equal left/right central
charges. Track twist-induced conformal weights. The entanglement coefficient
measures the summed charge once conformal ground-state scaling holds.
Finite-temperature entropy or finite-entanglement scaling can provide
additional checks, but were not run here.

## 9. Remaining mechanisms, ranked and tied to next tests

This is a qualitative physics ranking [heuristic], not a probability
assignment. The entries overlap. An emergent SU(2)_2 description, three
Majoranas and an Ising-plus-boson critical point can describe related
central-charge budgets. A parton rewrite alone is not another construction.

1. An isolated strongly interacting multicritical point with an additional
   neutral sector. Ising times a boson is the first target; an emergent
   SU(2)_2 or three-Majorana point is a more structured variant. This asks for
   a finite number of coincident gap closings and does not require a long
   cell. Neither LSM nor the on-site charge rank excludes it. The first test
   is a neutral gap closing along a calibrated c=1 locus, followed by the
   mixed-operator and two-velocity audit in section 7.1. Higher odd SU(2)
   levels allowed by a surviving spin anomaly are less economical variants,
   not ruled out by that anomaly. No candidate was found in this analytic run.

2. A strongly interacting, incommensurate two-fluid state whose dangerous
   conversions are oscillatory, irrelevant or tuned away. This is excluded
   on the exact U(1) NN locus but not on the entire residual family. It needs
   independently located soft branches and an actual covariance matrix,
   not the substitution of two XXZ K values. Test momenta, pair/string gaps,
   convexity and all low-dimension conversions as in section 4. C15 is the
   relevant companion analysis for the axial Z2 class.

3. A small spontaneous cell with both an active payload and a neutral
   transition or a second active internal degree of freedom. Period two
   or three costs less microscopic structure than a programmable long cell.
   The conventional dimer/trimer reductions and the exact-code obstruction
   eliminate several naive versions, but not arbitrary frustrated dressed
   codes. First demonstrate spontaneous order, the retained multiplets and
   their effective couplings. Then test the CFT inside a selected vacuum.
   A p-fold cat multiplet or an arbitrary Hamiltonian written on 2^p states
   does not suffice.

4. A critical constrained manifold with an unfamiliar neutral/coset theory.
   Large NN penalties can leave extensive degenerate spaces even at d=2.
   The witness theorem does not classify their kinetic splitting. This is
   a broader and less directed version of the previous mechanisms. Compute
   the projected dynamics and its finite-penalty corrections before assigning
   a clock, anyon or parton theory. C16 addresses the simplest blockade corner.

5. Long quantum-selected cells with extended or mobile markers and a
   high-dimensional low-energy payload. The frozen one-site marker and
   longest-FF-segment routes are closed at d=2, so this requires a different
   selector as well as a useful payload and controlled intercell dynamics.
   The first decisive test is a strict length-dependent bulk energy minimum
   and stable auxiliary gaps from one repeated qubit bond. No such selector
   is supplied. An abstract walk energy or an arbitrary block code is not it.

6. A genuinely unitary noncompact fixed point with more than one propagating
   field. It is not logically excluded by finite on-site dimension alone,
   but no construction or special reason to expect it in this NN family
   was identified. First establish finite velocities and oscillator towers,
   separating a dense zero-mode or z=2 limit. Decompactifying the known single
   boson cannot produce c>1.

The audit therefore closes several specific mechanisms and identifies the
common missing statement: a NN qubit bond must not generate a second
independent low-energy stress-tensor sector outside the controlled free,
XXZ/XYZ and completed-coset descriptions. Establishing that statement would
require new dynamical input. It does not follow from the number of qubit
states, the surviving anomaly, generic conversion relevance, or a change
of blocking.
