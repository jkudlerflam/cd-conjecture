# A 24-state conditional counterexample candidate

Status: **credible conditional candidate; full microscopic c=24 not
established**. A fresh audit and referee comparison withdrew the earlier
verified-counterexample wording. The finite bond, frozen material-sector
bounds, and 23-band reference are explicit and independently checked.
The infrared conclusion still requires finite-J local matching and a
uniform dressed-defect comparison, in addition to conventional Gaussian RG.
It is not a mathematical lattice-to-CFT theorem, a measured microscopic
c=24, or an assertion of exact uniqueness on every finite ring.

The construction adds one mobile vacancy to the explicit
23-state chain. Conditional on those premises, its infrared theory has
23 internal compact bosons and one charge boson, giving total Virasoro
central charge 24 and exceeding the local-dimension bound d-1=23. The construction avoids
the earlier expensive mobile wrapper and the species-graph interpreter.

## 1. The physical Hamiltonian is a fixed nearest-neighbour bond

Let h_23(g,lambda) be the exact bond specified in
[the explicit candidate](../explicit_candidate/README.md) and implemented
in [bond23.py](../explicit_candidate/bond23.py). Its alphabet contains 23
letters. Add one orthogonal state |vac>, let iota embed the occupied
letters, and let n=I-|vac><vac|. On C^24 tensor C^24 set

    h_24 = -t sum_(a=1)^23 (|vac,a><a,vac| + |a,vac><vac,a|)
           + J (iota tensor iota) h_23 (iota tensor iota)^dagger
           - (mu/2)(n tensor I + I tensor n).

Use the same h_24 on every bond of every periodic chain. This is Hermitian
and one-site translation invariant, with exactly d=24. Occupied number,
original internal charge, and original marker number are conserved.
All parameters, including positive g, lambda and J/t, are fixed before
L tends to infinity.

The entire original h_23, including its distributed on-site terms, acts
only on occupied pairs. Moving those terms outside the occupied block
would change the reference internal model. The new occupied-number
chemical potential mu is distinct from the original marker potential.

[bond check] `doping/bond24.py` constructs this lift. The saved 576-by-576
moderate-scale fixture has 923 nonzero entries; Hermiticity, occupied
block identity and the three charge commutators vanish in the check.
It is a structural fixture, not a numerically resolved n=23 violation.

## 2. Fix the certified internal geometry

Use n=23, material cell length b=24, unit internal data/express/switch
hopping, D_x=D_M=1000, Delta=4-2cos(pi/25), and zero species splitting.
The neutral express phase is pi/46 and the charged express phases vanish.
The port phase is pi/2. The weak-link phases are pi for the neutral
species and pi/2 for both charged species; other rail/switch phases vanish.

Define d_0=1, d_1=1004, d_(j+1)=1004 d_j-d_(j-1), beta=1/d_24, and

    weak-link amplitude = beta 2^(-136),
    0 < g_+=g_-=g <= 2^(-365),
    alpha_23 = sqrt(2/25) sin(pi/25),
    E_unit = 2 alpha_23^4 beta^2 g^2,
    physical internal data field = lambda^2 E_unit 2^(-214).

Choose the usual midpoint selector for b=24. Explicitly, with

    a_l=2[1-cos(pi/(l+1))],
    delta_l=(l+1)a_l-l a_(l+1),

its marker chemical potential is 2-(delta_24+delta_23)/2. Set the illegal
pair penalty K=10 and initially restrict 0<lambda<=2^(-20).

[proved and independently replayed] The C33 finite-loop, finite-g
certificate gives exactly 46 simple Fermi points and bounds

    delta_2 > 0.0140193,       delta_4 > 3.1558e-7.

The quartic check includes every unordered pair with replacement, so it
includes pair transfer between distinct bands. Its exact interval error
only decreases when g is reduced with the normalized field fixed.
The audit regenerated the 512-bit matrices with exact agreement and
reran the rational inertia and pair-product certificates. See
[the independent audit](route_a/C33_AUDIT.md).

[proved] On compatible material rings, every competing cell length and every illegal word has
positive energy cost at fixed nonzero lambda in the stated interval.
The all-length argument compares exact interacting open-cell energies,
not an extensive error with a finite defect energy. The certified selector
margin is gamma_24=7.566715709... times 10^(-5); the sufficient inequalities
are gamma_b>5 sqrt(3)|lambda| and K>4+7 sqrt(3)|lambda|. See
[the sector proof](route_a/README.md).

## 3. The vacancy completion is a local phase-model problem

[proved] Order the occupied sites and retain their internal labels. Let
g_j>=0 be the vacancies after occupied site j. On the material chain,

    sum_j g_j=L-N,
    H_gap=-t sum_j (B_j^dagger B_(j+1)+h.c.),
    H_coupling=J sum_j |g_j=0><g_j=0| h_23(j,j+1),

where B lowers a nonzero g by one with unit amplitude. Thus the apparent
nonlocality of squeezing is removed in material coordinates. The only
global qualification is the cyclic identification, retained as a boundary
twist. The infinite auxiliary gap alphabet is a change of coordinates;
it does not change the physical local dimension 24.

More precisely, the full finite ring also retains a root position x in
Z_L and the equivalence
(x;g_1,...;sigma_1,...) ~ (x+1+g_1;rotate(g);rotate(sigma)).
There are N equivalent roots, giving dimension
L binomial(L-1,N-1) 23^N/N = binomial(L,N) 23^N.
At fixed physical momentum the root becomes seam phases and a cyclic
consistency condition. Gap compositions alone are not the full physical
Hilbert space. The [fresh vacancy audit](../counterexample_audit_20261028/vacancy/REPORT.md)
derives this quotient and checks the direct physical matrix independently.

The gap chain is the q=0 phase model. At physical occupancy rho=1/2 its
hole density per occupied site is r=1. Exact thermodynamics gives

    e_gap(r)=-(2t/pi)(1+r)sin[pi/(1+r)],
    v_gap=t,       K_gap=4,
    <|g_j=0><g_j=0|>_bulk = alpha = 1/2-2/pi^2 > 0.

The leading internal Hamiltonian is therefore J alpha H_23. The exact
finite-ring charge-sea projection, including cyclic twists and their
possible degeneracies, is proved and checked independently in
[doping_check](doping_check/README.md). Individual rooted adjacency
probabilities need not be equal; cyclic projection supplies alpha.

[read] The unit-amplitude phase model and its exact thermodynamics are
standard: [Bortz and Sergeev, sections 2 and 4.1](https://arxiv.org/html/cond-mat/0603093)
give its q=0 velocity and dressed charge Z=1+r. In our conventional
boson normalization K=Z^2=4; the paper's symbol K differs by a factor pi.
Their hopping is 1/2 and must be multiplied by 2t.

## 4. Why the extra mode does not gap the payload

[standard bosonization, independently audited] The local contact projector
has a constant part, a smooth density current, higher derivatives/current
composites, and integer gap harmonics. At r=1 the first harmonic has
scaling dimension four; all higher harmonics have dimension 4m^2.
A contact current times a payload current is Gaussian marginal. A contact
current times the payload energy density has dimension three. The
nonforward material terms retain the certified momentum mismatches.

In physical coordinates, fractional period-24 marker harmonics are also
included. A payload vertex of cell momentum q and dimension Delta_pay,
dressed by an integer cell harmonic ell, has

    k_phys = rho(q+2pi ell)/b mod 2pi,
    Delta_ref = Delta_pay + [(q+2pi ell)/(2pi b)]^2.

At rho=1/2 a uniform resonance requires

    Delta_ref,res = Delta_pay + 4h^2,       h integer.

For h=0, only the original zero-momentum payload selection remains.
For h!=0 the charge contribution is at least four. Every non-Gaussian
uniform vertex consequently has reference dimension at least three.
There is a strict physical momentum margin >6.5745e-9 for the remaining
potentially dangerous oscillatory vertices. A full Gaussian relative
stiffness correction below 1/4 keeps these dimensions above two.
See [the infrared analysis](infrared/README.md).

[proved free-charge check] Boundary-twist changes alone have O(1/L) sea
energy costs, but conservation of physical momentum requires a charge
excitation. At half filling, arbitrary particle-hole excitations obey
E_exc >= (4t/pi)dist(P,pi Z). The nonforward internal vertices therefore
have threshold at least 2t delta_4/(pi b)>8.371e-9 t. The inequality
extends to bounded particle-number rounding on arbitrary rings; exact
momentum optimization through L=256 checks it. See
[the threshold proof](doping_check/CHARGE_THRESHOLD.md).

## 5. Nonzero couplings and defect stability

[controlled physics argument] Fix the loop and normalized field geometry.
Reduce g until the g^4 non-Gaussian contribution is small relative to the
g^2 logical hopping, including the slowest band and the complete finite
matching interval. Then choose nonzero lambda small enough for local
elimination of the gapped cell excitations to preserve those margins.
Finally hold this complete h_23 fixed and choose J/t positive and small.
Local phase-model matching then changes internal terms by O(J^2/t),
whereas their leading scale is J. Its Gaussian mixing, measured relative
to the smallest retained velocity, tends to zero with J/t. The exact
charge oscillator is retained throughout; it is never globally projected
out using a finite-ring gap.

The order of choices is essential. Decreasing lambda alone does not
suppress the internal g^4/g^2 ratio, and the fast vacancy scale must be
compared with the smallest internal speed, not merely the largest band.
Every relevant tolerance is positive at this fixed n. This supplies a
nested weak-coupling regime at fixed, nonzero microscopic coefficients,
using conventional strong-coupling matching and Gaussian RG. The stable
neighborhood is within the exact number-conserving symmetry family;
stability under arbitrary symmetry-breaking changes is not required.
This is a controlled existence regime, not a certified numerical interval
for the final lambda and J. Pure charge stiffness can receive O(J)
corrections, small relative to t; the O(J^2/t) estimate above concerns
induced internal terms. Raw mixed current coefficients can be O(J),
with normalized stiffness vanishing as sqrt(J/t).

Defect energies require a local comparison. The leading static gap-density
susceptibility is local, with a separately retained canonical zero-mode
term. Contact descendants give summable tails with a finite first spatial
moment; the first harmonic gives a faster |j|^-7 tail. Thus local defects
and interfaces acquire O(J^2/t) dressing per defect, while their positive
bare costs are O(J). The independent tagged-gap test reaches L=20 and
checks single, paired and finite-density defects. It supports this
matching argument; it does not replace it. See
[defect transfer](route_a/DEFECT_TRANSFER.md).

The [fresh infrared audit](../counterexample_audit_20261028/infrared/README.md)
checks quantum sources beyond static scalar defects. Connected local
material operators with dimension Delta_min>1/2 give a charge-mediated
kernel whose first spatial moment converges; the quantitative Gaussian
bound supplies Delta_min>=sqrt(3/5). Ports contribute at most
O(J|lambda|) per affected boundary. Arbitrarily long wrong cells are
handled with a separate linear-in-length selector margin, without using
their shrinking walker gaps as uniform denominators. These are controlled
physics estimates; a full all-order microscopic per-defect inequality
is not asserted.

[proved local normal form; standard infrared inference] In the selected
material sector, block each complete 24-site cell. Its finite-dimensional
ground space has rank 2^23 and a positive auxiliary gap Delta_0. The local
Schrieffer-Wolff construction of
[Bravyi, DiVincenzo and Loss, section 4.4 and Lemma 4.2](https://arxiv.org/pdf/1105.0675)
gives an exact full-space identity at any fixed order m,

    U_m H_23 U_m^dagger = D_m + R_m,
    ||R_m||_local,exponential <= C_m |lambda|^(m+1).

D_m is locally block diagonal. The bound is on interaction strength,
independent of chain length; it is not merely an extensive energy-error
estimate. The constants include the fixed cell norm/gap ratio. At m=2,
the remainder divided by the smallest payload scale tends to zero as
lambda tends to zero at fixed g.

There is also a uniform auxiliary threshold above the reference ground.
Apply the source's block-diagonal stability lemma after subtracting
(Delta_0/2) sum_j Q_j, where Q_j projects out of the cell code. This does
not change the projected Hamiltonian or its ground energy E_eff, and gives

    D_m >= E_eff I + (Delta_0/2) sum_j Q_j.

The final internal step therefore uses conventional stability of a
Gaussian phase with massive auxiliary modes under a sufficiently small
symmetry-preserving quasi-local perturbation. It does not require C40's
stronger exact all-order low-frequency kernel. No all-order block
diagonalization or mathematical CFT theorem is claimed. See the
[local normal-form audit](route_a/M2_REVIEW.md),
[independent primary-source review](infrared/SOURCE_REVIEW.md), and
[physics assessment](infrared/PHYSICS_VERDICT.md). This argument applies
after selected-sector blocking; the vacancy charge mode remains gapless
and is handled by the separate phase-model analysis above.

## 6. Rings and the conformal conclusion

The selected internal crystal uses N=24M occupied sites. On any sufficiently
large physical ring, including prime lengths, N can be rounded to the
nearest such value at density 1/2, costing O(t b^2/L). A wrong material
cell defect has a fixed positive energy of order J, so this rounding
wins asymptotically. Tune the physical chemical potential to half density.

The b possible internal origins become charge boundary twists. They are
part of the compact-boson zero-mode gluing, not b independent frozen
physical crystals. Physical marker order has a positive power-law exponent
2/b^2 at the reference. A separately solved cyclic-color control checks
this gluing on prime and composite rings. Accidental finite-size charge
crossings do not change the oscillator count.

The original conjecture explicitly permits a collection of unitary CFTs
with different propagation speeds. It does not demand exact uniqueness
on every finite ring. The compact-boson gluing here preserves 24 oscillator
modes; see the independent [scope check](doping_check/CONJECTURE_SCOPE.md).

Conditional on finite-J matching, dressed-sector stability, and the
stated Gaussian RG inputs, the phase contains 23 payload bosons and one
vacancy boson:

    c_total=24 > d-1=23,       d=24.

This is a finite-d construction; proving arbitrarily large c at fixed d
would require more than the n=23 certificate and is not claimed here.
A numerical microscopic central-charge measurement and a fully rigorous
thermodynamic CFT theorem remain distinct stronger validations.

## Independent checks and remaining scope

Three parallel AI workers reviewed the full synthesis and accepted the
conclusion at the stated physics standard. Their reports retain the
continuum assumptions rather than treating finite numerical tests as a
measurement of c=24. The independent bond reconstruction has zero entry
error; see [the synthesis audit](doping_check/synthesis_audit.json) and
[referee report](doping_check/REFEREE.md). The collaborator's separate
[C41 progress report](../../codex/work_C41/PROGRESS.md) independently
recovers the physical harmonic and number-rounding rules; that report
was still in progress when this synthesis was completed.

A subsequent [fresh adversarial audit](../counterexample_audit_20261028/README.md)
independently rebuilds the microscopic certificate, performs new literal
two-cell comparisons and complete low-order vertex scans, and tests
noncommuting quantum-impurity dressing. Its referee cross-review withdrew the
earlier verified-counterexample wording: the finite ingredients and
conditional RG argument pass, but uniform finite-J defect matching is
not yet established. No actual refutation of the model was found. The finite-ring root quotient
and the refined quantum-interface bound above incorporate its corrections.

## Provenance

Coordinated by a GPT-6-family Codex assistant with three parallel workers,
reused across independent subtasks, on 2026-10-28. The explicit internal
bond and original C33 certificate are collaborators' work; the vacancy
lift, sector completion, local phase-model argument and independent
replays are this campaign's additions. All finite checks and limitations
are retained in the linked folders. No human review is claimed.
