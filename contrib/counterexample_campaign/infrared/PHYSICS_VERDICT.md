# Does the vacancy lift work at the project's physics standard?

**Verdict: yes, the infrared completion is a conventional controlled
weak-coupling construction, provided the explicit frozen h23 has the
claimed selected-cell sector and the C33 certified 23-band geometry.**
The vacancy step does not require an additional novel stability
principle. Its exact material-coordinate representation is a local
gap fluid coupled weakly to the internal chain. At half occupancy the
gap fluid has K=4, which leaves a strict margin against every local
charge locking term. Its forward coupling to the payload is Gaussian;
the remaining terms are irrelevant or have the original nonzero
payload mismatch. Choosing J/t last controls their effect even on
the very slowest payload band.

This is a theoretical-physics conclusion using bosonization,
Wilsonian matching and weak-coupling RG. It is not a proof of a
microscopic CFT limit in the mathematical sense, nor a numerical
measurement of c=24. The original fixed-bond implementation and
sector estimates remain separate finite tasks; an incorrect bond or
a lower competing sector would invalidate this verdict. Merely
asking for an all-order mathematical CFT theorem after the following
calculation would impose a stronger standard than the repository's.

## 1. Exact local material-coordinate Hamiltonian

[proved] Fix N occupied sites and read them cyclically. Let g_j be the
number of vacancies after occupied site j and before occupied site
j+1. Then g_j is a nonnegative integer and

    sum_j g_j=L-N.

The internal state a_j belongs to C^23. If B|g>=|g-1> for g>=1,
B|0>=0, and B^dagger|g>=|g+1>, a vacancy hop transfers one unit
between adjacent g's with unit amplitude. Thus, up to the rooted-ring
quotient and its boundary twist, the exact bulk Hamiltonian is

    H_gap=-t sum_j (B_j^dagger B_(j+1)+h.c.),
    H_int=J sum_j A_j h_(j,j+1),       A_j=|g_j=0><g_j=0|.           (1)

This is local in the occupied-particle label j. B is a bounded shift,
not a canonical boson operator with square-root occupation factors.
It is the q=0 limit of the q-boson chain, also called the phase model;
the operator representation and Hamiltonian appear in equations
(1), (6), (7) and the q-oscillator matrix elements of
[Bortz and Sergeev, *The q-deformed Bose gas: Integrability and
thermodynamics*](https://arxiv.org/html/cond-mat/0603093).
The infinite gap alphabet does not enlarge the physical local
dimension, which remains 24. At every finite ring its total
occupation is bounded by L-N. Equation (1) is a representation used
to analyze the physical finite-alphabet chain.

The cyclic quotient supplies boundary twists and zero-mode gluing;
it does not change the local operators in (1). A local perturbative
expansion of (1) generates local material operators. In particular,
fractional gap-number vertices are not generated. The physical
fractional *marker* harmonics arise from changing coordinates, as
explained in the other audit file, and are fully consistent with (1).

## 2. Exact gap-fluid thermodynamics fixes its Luttinger parameter

[proved thermodynamics; standard Gaussian identification] Set

    r=(L-N)/N,       rho=1/(1+r).

The free hard-core charge energy, expressed per material site, is

    e_gap(r)=-(2t/pi)(1+r) sin[pi/(1+r)].                           (2)

Physical and material wavevectors differ by rho, hence

    v_gap=2t rho sin(pi rho),
    e_gap''(r)=2pi t rho^3 sin(pi rho).

Using e''=pi v/K for a boson whose integer number is total gap
occupation gives

    K_gap=1/rho^2.

The current stiffness gives the same answer: a uniform hopping phase
theta changes (2) by a factor cos(theta), so its second derivative
is v_gap K_gap/pi. At rho=1/2, r=1,

    v_gap=t,       K_gap=4.                                       (3)

[read, independent primary-source check] Bortz and Sergeev's section
4.1, equations (42)-(45), gives the q=0 velocity and susceptibility,
with dressed charge Z=1+r. In the conventional normalization used
here, K=Z^2=(1+r)^2. Their hopping amplitude is 1/2, so their velocity
is multiplied by 2t to obtain (3). Their symbol K includes an extra
factor 1/pi; it should not be copied without converting conventions.

Its local vacuum-contact expectation is

    alpha=<A>=rho-sin(pi rho)^2/(pi^2 rho),
    alpha(1/2)=1/2-2/pi^2>0.

Here <A> means the cyclic bulk expectation. Equal individual rooted
expectations are not asserted on a finite ring.

## 3. The operator expansion has a strict charge margin

[standard bosonization] A local, number-conserving material operator
such as A_j has the low-energy form

    A_j=alpha+a_1 partial_j phi_gap
        + terms quadratic in currents or with more derivatives
        + sum_(m!=0) a_m exp[2im(pi r j-phi_gap)] + ... .            (4)

Its first harmonic has scaling dimension K_gap=4, and the m-th
harmonic has dimension 4m^2 at the reference. At r=1 these harmonics
are commensurate in the material lattice, but they are irrelevant,
with a strict dimension margin. The only nonconstant charge operator
below dimension two is the smooth conserved current. Charge-current
quadratics are Gaussian marginal terms. Exact gap-number conservation
forbids vertices of the dual phase.

The internal chain breaks material translation into cells of b=24
sites. Its order harmonics have momenta 2pi ell/b. They multiply (4)
but do not create a fractional local gap vertex. A nonzero period-b
wavevector stays oscillatory. A zero-wavevector order component only
changes alpha or the chemical potential. The physical-coordinate
calculation independently yields

    Delta_uniform=Delta_payload+4h^2,

including fractional physical marker harmonics. The two coordinate
descriptions therefore give the same absence of a relevant locking
operator.

## 4. Matching and weak coupling

[controlled physics argument] At J=0 the internal degrees of freedom
are static. Integrate charge fluctuations above a fixed material
cutoff, leaving the charge boson itself in the theory. The identity
part of (4) gives J alpha H23. The cutoff is fixed before decreasing
J/t, so higher matching coefficients have finite denominators and
start at order J^2/t times finite cutoff-dependent constants.

There is no need to project onto one charge ground state or divide by
its O(1/L) gap. In nonforward channels the momentum-restoration bound
in `CHARGE_THRESHOLD.md` gives a strictly positive free threshold.
The zero-momentum smooth sector remains explicitly in (4).

In the latter sector, the density-density propagator has Fourier
form proportional to

    v_gap K_gap q^2/(omega^2+v_gap^2 q^2).                          (5)

The nonzero-momentum static limit is a finite constant, hence a local
susceptibility, not an infrared-divergent attraction. In a canonical
gap sector the zero mode is absent, so on N material sites the
forward energy-response kernel is

    chi_forward(i,j)=C[delta_(i,j)-1/N],
    C=[alpha'(r)]^2/[2 e_gap''(r)],                                (6)

in the convention E^(2)=-J^2 sum_(i,j) s_i chi(i,j) s_j. At r=1,
C is approximately 0.07858/t. The uniform rank-one pedestal must
be kept separately: its absolute first moment is not uniformly
bounded in N. Its energy is +J^2 C(sum_j s_j)^2/N. If k localized
defects or phase interfaces change the source sum by O(k), relative
to a bounded periodic background with sum O(N), its energy change
is O(J^2[k+k^2/N])=O(J^2 k), since k<=N. A translated periodic
domain has the same bulk source sum; only its endpoint fragments
contribute to this change.

After this rank-one part is separated, the offsite static integral
of the leading smooth density correlator vanishes away from contact.
The first oscillating contact term has dimension
four, so its static large-separation tail is O(|j|^-7). Dimension-two
smooth composites and descendants can have the slower O(|j|^-3)
static tail and must also be included. A dimension-Delta two-point
function has static envelope O(|j|^(1-2Delta)); the surviving
nonlocal terms here have Delta>=2. Mixed dimension-one/two Gaussian
correlators vanish by odd Gaussian field number or reduce to
derivative/contact terms. Consequently the nonlocal static kernel
has a finite **first spatial moment** after the rank-one subtraction,

    sum_j (1+|j|) |chi_remainder(j)| < infinity.                    (7)

The first moment, stronger than ordinary summability, bounds the
number of induced bonds crossing a defect interface. Thus
charge-mediated static changes of the internal local energy are
O(J^2/t) per site, or per localized defect/interface after vacuum
pieces are subtracted. Finite constants depend on the fixed h.

One must retain (5), rather than use its static value, when counting
the charge oscillator. In the coupled Gaussian problem it generates
the ordinary density-current mixing. Its normalized coefficient is
small as J/t tends to zero at fixed internal parameters. For example,
an unnormalized mixed stiffness of order J has relative magnitude
O(sqrt(J/t)/sqrt(v_min_internal)); v_min_internal is fixed and positive
before J is chosen. The joint quadratic form stays positive for
sufficiently small J/t. This prevents phase separation as well as a
spurious loss of the fast charge mode.

After projection onto the internal low-energy fields, a charge
current times a payload current has dimension two and is Gaussian.
A charge current times payload energy has dimension three. Every
contact harmonic has dimension at least four even before its payload
factor. All remaining neutral payload non-Gaussian vertices are
either momentum-detuned or of reference dimension at least three.
The independently checked mismatch is >6.5745e-9 in physical units.

This is sufficient at the ordinary RG standard. After the finite
matching interval, only a positive Gaussian matrix and strictly
irrelevant non-Gaussian vertices remain. No marginally relevant
coefficient with an unknown sign has been omitted.

## 5. Three scales must be chosen in the right order

[controlled physics argument] The complete finite-parameter choice
can be made in this order:

1. Fix the n=23 finite loop and shifted field ratios for which C33
   certifies all 46 simple Fermi roots and every nontrivial quartic
   momentum mismatch. They give a fixed positive minimum velocity
   and a finite, potentially long, matching interval.
2. Decrease g, if necessary, below the C33 certified upper value.
   The normalized band geometry tends to its certified limit, while
   the nonquadratic g^4 contribution divided by the g^2 hopping
   tends to zero. Fix a nonzero g that keeps all RG and filling
   margins, including multiplicity factors and logarithms.
3. Choose a fixed nonzero lambda small enough that connected
   corrections beyond the intended lambda^2 logical Hamiltonian are
   inside those margins and the frozen cell defects stay massive.
   This is the usual local strong-coupling elimination of finite,
   gapped internal machinery; it is a low-energy operation, not a
   global isolation of an extensive logical spectral band.
4. With this complete h23 fixed, choose J/t positive and small enough
   that all O(J^2/t) induced internal terms are small relative to J
   times the minimum payload velocity, the sector gap, and the
   finite mismatch tolerance. At the same time the normalized
   charge-payload Gaussian mixing tends to zero. Tune mu to the
   dressed density 1/2. Its reference compressibility is strictly
   positive, so this tuning is regular.

Each requirement is a strict inequality in finite constants. None
requires a uniform bound as n grows. Every microscopic coefficient
is fixed before L tends to infinity. The hierarchy may be far beyond
any feasible experiment or simulation; that does not affect the
existence claim at the stated theoretical-physics standard.

The final theory has the 23 payload bosons plus the gap-fluid boson,
so total Virasoro c=24 at physical local dimension d=24. The
stronger fixed-d unbounded-c conclusion is not needed: this one
finite n already gives c>d-1. Its central charge has not been
numerically measured, and exact finite-ring uniqueness is not used
to infer it.

## 6. What this verdict does and does not settle

[assessment] I find no remaining *specific physical mechanism* that
would obstruct the vacancy completion once the frozen finite bond,
selected sectors, and nonresonant 23-band phase are accepted. The
material gap fluid makes its matching an ordinary local RG problem.
It is reasonable to call the resulting chain a physics-level
counterexample after those finite checks are independently accepted.

A completely rigorous mathematical result would still need to
construct the local low-frequency action and prove its limiting CFT
correlations. That general theorem is not supplied. Conversely, a
new lower-energy defect sector, an error in the explicit bond, or a
failure of the C33 certificate would be substantive failures at the
physics standard too. Those are the concrete checks to continue.

[scope after reviewing C40] The newly synchronized C40 report
explicitly leaves the frozen core's exact connected low-frequency
remainder as an input, called M2 there. The present verdict does not
claim to prove M2 by the gap-fluid calculation. The gap-fluid
calculation supplies a conventional local matching argument for the
*additional vacancy step*. The frozen chain still uses the usual
Wilsonian treatment of gapped finite-cell excitations. A rigorous
version of that input remains distinct from finite-order SW and from
the finite band certificates.

## 7. Is C40's M2 an unresolved physical mechanism?

[assessment] No specific such mechanism is identified. M2 is the
standard local strong-coupling/Wilsonian matching input for gapped
internal auxiliaries, applied here at one fixed finite n. The finite
cell has a positive gap; a local auxiliary defect changes only
finitely many interactions. Its energy correction is small when
lambda is chosen below that fixed gap, and logical modes at energies
of order lambda^2 are retained throughout. The finite n can make
every constant extremely unfavorable without changing this scale
separation. This is the same physical use of integrating out gapped
charge or orbital states that underlies ordinary strong-coupling
spin Hamiltonians.

C36 correctly disproves a stronger operator-norm claim made on
**all** retained energies, where extensive logical excitations can
resonate with the machinery. That counterexample does not exhibit
an instability of the low-frequency ground-state theory. Likewise,
vacuum dressing invalidates a bare global-Q threshold, but it does
not make a locally dressed auxiliary gapless. Neither invalid
argument is used as the physical justification here.

A nonperturbatively small generated term would still respect the
exact U(1) and locality and enter the same stable nonresonant
Gaussian basin, provided its coefficient tends to zero in the
strong-coupling limit. There is no allowed relevant resonance whose
unknown exponentially small sign would invalidate the construction.
The actual failure that would matter would be a new low-energy
defect phase or a finite normalized correction persisting as
lambda tends to zero. No route to either appears in this finite,
locally gapped core with independently chosen small couplings.

Accordingly, I would retain M2 as an explicit conventional physics
input, and accept the fixed-n infrared conclusion at the project's
standard once the finite sector and bond checks are accepted. I
would not label M2 itself a proved convergent expansion or present
the resulting claim as a rigorous mathematical lattice-CFT theorem.

## 8. Final adversarial operator audit

[proved, local algebra] The physical occupied-number U(1) and the
original internal-charge U(1) are independent exact symmetries.
Within the zero-defect code the latter is total logical fermion
number. An anomalous payload pairing term therefore cannot be
made neutral by a vacancy-sector vertex: changing vacancy charge
would violate the other U(1), and a charge current carries no net
internal charge. Virtual auxiliary excitations return to the same
code representation and cannot supply a missing final charge.
Dual gap-field vertices are forbidden by exact gap-number
conservation in material coordinates, equivalently occupied-number
conservation at fixed physical L.

[proved, locality/parity check] A finite-support internal-neutral
operator in the selected material sector is an even polynomial in
the logical Jordan-Wigner fermions. Any exterior parity strings
cancel; strings between its finitely many endpoints remain within
its finite support. A physical finite-range operator intersects
only finitely many occupied labels, regardless of vacancies.
Quasi-local elimination preserves this local operator algebra.
Consequently neutral half-integer payload charge vectors or
standalone Jordan-Wigner disorder strings are not additional local
Hamiltonian vertices. The fractional **physical marker** harmonics
are different: they are local observables and are precisely the
integer ell dressings already retained in equation (1) of the
companion README. Their resonance dimension is already included.

[proved within the free-boson operator dictionary] For an integer
payload chiral vector z with sum_i z_i=0, its free primary dimension
is sum_i z_i^2/2. Dimensions one and two give exactly the neutral
bilinear and quartic momentum tests. Higher-degree fermion
monomials either contract to these lower classes, become current
descendants, or have primary dimension at least three. Normal
ordering and point splitting cannot create an omitted lower
dimension primary. Repeated-branch quartics are covered by the
pair-with-replacement certificate; their surviving descendants
have additional derivative dimension.

[controlled physics conclusion] At half occupancy every uniform
non-Gaussian charge dressing with nonzero reciprocal index h has
charge contribution 4h^2. The h=0 cases reduce to the certified
payload selection. Gaussian charge-current, payload-current and
cross-current couplings must all be retained in the full positive
stiffness matrix; they do not remove an oscillator when its
relative perturbation is sufficiently small. The momentum and
filling-shift margins are separate requirements and are reserved
along with the stiffness margin. This audit finds no additional
allowed relevant or marginally relevant vertex.

The resulting stable basin is open **within the exact
number-conserving symmetry family**, with all couplings nonzero
and fixed before the thermodynamic limit. No claim of openness
under arbitrary symmetry-breaking perturbations is needed or
made. This qualification does not weaken the existence of a
counterexample bond. Possible finite-ring zero-mode crossings and
unequal mode velocities are likewise not additional local gapping
operators.
