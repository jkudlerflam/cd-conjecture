# Independent physics referee: the d=24 vacancy completion

**Verdict: endorse the c=24 conclusion at the repository's stated
theoretical-physics standard.** The current construction has a literal
24-state nearest-neighbor bond, independently checked finite band and
sector data, and a conventional weak-coupling route to 24 compact bosons.
I found no remaining specific instability or missing microscopic mechanism
in the vacancy completion after its exact local gap representation was
supplied. This is an analytic physics result based on controlled scale
separation and Gaussian RG, not a rigorous theorem about lattice CFT
correlations or a numerical measurement of c=24.

This review read `infrared/PHYSICS_VERDICT.md`, the route-A global-sector
proof and C33 independent audit, and checked the primary phase-model
source. Its author independently implemented the finite-ring squeezing,
charge-sea projection, cyclic-color gluing, and all-occupation charge
threshold tests in this folder. It is not a second human referee report.

## 1. The main objection was concrete, and is resolved

The naive argument averaged adjacency over a gapless charge ground state.
That does not give a uniform perturbative elimination because the charge
gap vanishes with ring length. A related concern was that changing material
momentum changes a charge boundary twist with only O(1/L) sea-energy cost.
Neither argument is sufficient on its own.

The revised argument keeps the charge mode and uses two additional facts.
First, the exact material variables are the integer vacancy counts g_j
between consecutive occupied sites. Their hopping is local,

    H_gap=-t sum_j(B_j^dagger B_(j+1)+h.c.),
    H_coupling=J sum_j |g_j=0><g_j=0| h_(j,j+1).

The shift B has unit matrix elements. The contact projector is a bounded
local operator in these coordinates. The rooted-ring quotient changes
boundary gluing, not bulk locality. Thus the required matching is a local
coupled-chain problem rather than elimination of a nonlocal adjacency
string by an unjustified global gap.

Second, physical momentum conservation forces a finite charge threshold
for fixed nonforward material transfer. The independent proof in
[CHARGE_THRESHOLD.md](CHARGE_THRESHOLD.md) includes arbitrary numbers of
particle-hole excitations and bounded density rounding on arbitrary ring
lengths. Its limiting lower bound for cell momentum q is

    2t dist(q,2 pi Z)/(pi b).

It does not eliminate forward charge excitations. Those remain in the
coupled Gaussian theory.

## 2. The charge parameter and local operator list are consistent

The gap chain is the q=0 phase-model limit of the q-boson chain. The
primary paper's operator algebra gives unit-amplitude lowering, and its
exact velocity and dressed charge agree with the proposed normalization.
The paper's hopping is 1/2, compared with t here; its symbol K also differs
by a factor of pi from the conventional bosonization parameter.
[Bortz and Sergeev, sections 2.1 and 4.1](https://arxiv.org/html/cond-mat/0603093v1).

Independently, with gap density sigma=(1-rho)/rho,

    e_gap(sigma)=-(2t/pi)(1+sigma)sin[pi/(1+sigma)],
    v_gap=2t rho sin(pi rho),
    e_gap''=2pi t rho^3 sin(pi rho).

Therefore K_gap=pi v_gap/e_gap''=1/rho^2. At half occupancy this is
K_gap=4 and v_gap=t. Any stiffness derivative with a uniform hopping phase
is understood in a fixed adiabatic current branch; minimizing over all
windings first would gauge away that phase in the thermodynamic limit.

The contact projector has a constant part alpha, a smooth density current,
current composites and derivatives, and integer gap-number harmonics.
Its first nontrivial harmonic has dimension four. Exact gap-number
conservation forbids dual-number vertices. This supplies a strict margin
against the charge locking that would otherwise be a serious concern at
integer gap density.

The physical marker harmonics are fractional and cannot be omitted.
The independent cyclic-color control and locality review include them.
At half occupancy, a uniform physical vertex obeys

    Delta_ref = Delta_payload + 4 h^2,  h in Z.

Thus the enlarged physical harmonic lattice produces no new relevant
uniform operator. The material and physical coordinate descriptions
agree on this conclusion.

## 3. Static defect response does not expose an infrared divergence

The leading smooth-density static susceptibility is local in the
Gaussian theory. Its gapless dynamical propagator must be retained when
counting modes, but it does not yield an infrared-divergent static
attraction between defects. Dimension-four contact harmonics give a
static tail proportional to |j|^-7.

A minor refinement to the verdict's wording is useful: uniform
dimension-two composites can also contribute static |j|^-3 tails.
Those are harmless too. For an interface-energy comparison one needs a
finite **first spatial moment**, not just a summable kernel. Both the
|j|^-3 and |j|^-7 tails have that property; the leading smooth term is
local. There is one finite-volume qualification: fixed total gap number
removes the zero mode, giving the canonical smooth kernel
`C(delta_ij-1/N)`. Its uniform pedestal must be treated separately; it
does not have a uniformly bounded absolute first moment. The associated
energy term is proportional to `+J^2 C (sum_j s_j)^2/N`. Changing k bounded
local sources changes it by O(J^2 k), including the quadratic correction
because k<=N. Translated bulk origins have the same mean source. After
separating this rank-one term, the remainder has the finite first moment
used above. This is the conventional reason that boundary/defect energies,
after subtracting identical bulk pieces, receive finite O(J^2/t)
corrections. The statement is a local low-frequency estimate, not a
comparison of an O(J^2 L) total-energy bound against one defect gap.

The frozen material proof supplies a strictly positive defect cost gamma
at fixed nonzero lambda. In the lifted chain its leading cost is
J alpha gamma. Finite corrections smaller by J/t cannot overturn it
when J is chosen sufficiently small after h is fixed. Relevant impurity
flows in the slow payload, if present at a deliberately inserted defect,
do not create an extensive negative correction; the original material
boundary energies are already included in the sector comparison.

## 4. The very small internal velocities are an engineering cost

The n=23 hierarchy is extreme, but the relevant constants are finite.
The independent C33 audit certifies simple roots and pair-momentum
margins at that one fixed n. Decreasing g suppresses the nonquadratic
second-order term relative to the quadratic hopping. Decreasing lambda
then suppresses higher connected internal corrections. These are local
strong-coupling operations on finite gapped cells, followed by the
ordinary nonresonant multichannel RG analysis.

Only after fixing that h is J/t decreased. A mixed stiffness of order J
has normalized size of order

    sqrt(J/t)/sqrt(v_min(h)),

up to fixed coefficients. It therefore tends to zero even for the very
slowest certified band. Charge compressibility remains positive, and all
24 Gaussian modes survive. A claim based only on J/t being numerically
small, without the minimum-velocity normalization, would fail; the
revised construction uses the correct hierarchy.

All remaining neutral non-Gaussian payload vertices are either
momentum-detuned or have reference dimension at least three. Gaussian
forward scattering changes the stiffness matrix continuously; it does
not conceal an unknown marginally relevant non-Gaussian coupling. The
stated strict dimension and momentum margins give a conventional open
weak-coupling neighborhood at fixed n.

## 5. What the endorsement means

The resulting physical chain has d=24 and c=23+1=24, exceeding d-1=23
at the declared physics standard. The finite arbitrary-ring gluing
checks support the interpretation of material origin branches as charge
zero modes, not independent frozen physical crystals. Bounded particle
rounding can satisfy the selected material cell length on every
sufficiently large physical ring. Isolated finite-size sector crossings
do not change the central charge and need not be excluded to make this
statement.

The frozen-cell low-frequency remainder called M2 in the collaborators'
C40 report is **not** proved as an all-order mathematical statement by
this review or by the gap-fluid argument. The physics endorsement uses
ordinary local elimination of gapped finite-cell excitations as an
effective-theory input, together with its scale hierarchy. The new local
gap representation addresses the additional vacancy matching step; it
must not be described as a proof of C40's separate M2 assertion.

The remaining work appropriate to a publication is independent expert
refereeing, a compact self-contained presentation of the explicit bond
and parameter hierarchy, and additional accessible-size numerical tests
of the mechanism. A mathematical lattice-to-CFT theorem would be a
stronger result. Its absence is not a newly identified physical
obstruction, and this report does not demand it as a prerequisite for
the repository's stated standard.

Provenance: independent GPT-6-family research worker, 2026-10-28.
The coordination thread supplied the local phase-model observation and
requested this adversarial review; the checks and reasoning above were
performed independently of the infrared verdict's author. No human
review or measured c=24 is claimed.
