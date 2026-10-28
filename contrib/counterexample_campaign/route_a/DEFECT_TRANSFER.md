# Defect energies after adding a mobile vacancy

## Outcome

[physics-level inference, with explicit matching input] Weak coupling to
the vacancy fluid does not supply an infrared mechanism that cancels the
positive material-defect energies of the selected 23-letter chain.
In material coordinates the free vacancy fluid is the phase hopping model
with Luttinger parameter K_gap=rho^-2. Near rho=1/2 this is near four.
The static response to a tagged material-contact defect is therefore
infrared finite. Its induced interaction between well-separated defects
is summable, including its first spatial moment. This supports a correction
of order J^2/t per defect, possibly with harmless logarithmic factors,
against the leading positive cost of order J.

[proved finite statements; numerical test] The exact gap-coordinate toy
below verifies the leading coefficient and computes charge relaxation
around one defect, adjacent defects, and a positive density of defects.
It is not a simulation of the complete d=24 Hamiltonian. The distinction
between a static scalar contact potential and the operator-valued
23-letter interaction is retained throughout.

The missing mathematical statement is a uniform microscopic bound on the
linked defect correction for that operator-valued interaction. The
continuum argument identifies why the relevant response is finite; it
does not silently replace this statement by an extensive O(J^2 L) bound.

## 1. Exact charge coordinates and the tagged-contact toy

[definition] There are N ordered material particles on L physical sites.
Let g_j>=0 count the vacancies after particle j, with

    sum_j g_j = L-N,        rho=N/L.

At total physical momentum zero the vacancy Hamiltonian is exactly

    H_gap = -t sum_j (E_j^+ E_(j+1)^- + h.c.),             (1)

where E^-|g>=|g-1> for g>0 and is zero at g=0, and
E^+|g>=|g+1>. There are no square-root occupation factors. A hop of a
physical particle transfers one vacancy between its two adjacent gaps.
A distinguished material defect roots the cyclic material word, so
there is no quotient by cyclic permutations of the g_j.

[proved] A material bond j is active precisely when its adjacent particles
touch. Thus its charge operator is

    A_j = 1[g_j=0].                                      (2)

The Hilbert-space dimension in (1) is binomial(L-1,N-1). A scalar defect
cost is represented exactly by H_gap+J sum_j d_j A_j. This keeps the ring
geometry and lets the charge state adapt to each defect. It does not
freeze the charge state at its free expectation.

[proved] The ground energy and contact probability at J=0 are

    E_0(L,N) = -2t sin(pi N/L)/sin(pi/L),
    alpha_L = N/L - sin^2(pi N/L)/(N L sin^2(pi/L)).       (3)

Equation (3) follows from the spinless-fermion filled sea, with its usual
parity boundary twist. The sum of all A_j is the number of physically
adjacent occupied pairs. Wick's theorem gives its expectation, and
material translation symmetry divides it equally among N bonds.
The thermodynamic value is

    alpha(rho) = rho-sin^2(pi rho)/(pi^2 rho)>0.           (4)

[proved, finite volume] If psi_0 is the normalized free ground state, the
static susceptibility matrix is

    chi_ij = <(A_i-alpha_L)psi_0,
              (H_gap-E_0)^(-1)_perp (A_j-alpha_L)psi_0>.

For fixed finite L, ordinary nondegenerate perturbation theory gives

    E(Jd)-E_0 = J alpha_L sum_j d_j
                -J^2 sum_ij d_i chi_ij d_j + O_L(J^3).   (5)

No claim of uniformity follows from finite-volume perturbation theory
alone: its spectral gap vanishes with L. Sections 2-4 address that issue.

## 2. Why the tagged-contact response is finite

[derived using standard Luttinger-liquid thermodynamics] Put
r=(L-N)/N=(1-rho)/rho, the vacancy density per material gap. The exact
energy per gap and the sound velocity in material coordinates are

    e_gap(r) = -(2t/pi)(1+r) sin[pi/(1+r)],
    v_gap = 2t rho sin(pi rho),
    e_gap''(r) = 2pi t rho^3 sin(pi rho).

Using e''=pi v/K in the convention where exp(2i phi) has dimension K,

    K_gap = rho^-2.                                      (6)

At half material filling, v_gap=t and K_gap=4. This is a coordinate and
compactification change from the physical free-fermion charge field; one
must not assign K=1 to a local **material-gap** impurity. The phase model's
relation to free fermions is also described by Bogoliubov and Timonen,
[*Correlation functions for a strongly coupled boson system and plane
partitions*](https://doi.org/10.1098/rsta.2010.0322).

[standard bosonization input] A local neutral contact observable has
the expansion

    A_j = alpha(r) + alpha'(r) delta r_j + descendants
          + sum_(m>=1) a_m cos[2m(pi r j-phi_gap)] ,
    delta r = -partial_j phi_gap/pi.                     (7)

The smooth term has dimension one; the first oscillatory term has
dimension K_gap=4 at rho=1/2. A point impurity in the oscillatory term is
irrelevant for K_gap>1; even a coherent bulk modulation is irrelevant
for K_gap>2. These strict inequalities survive a small change of density
or weak coupling to the material sector.

[derived continuum response] For the Gaussian density term, the
zero-frequency response in momentum space is a constant compressibility.
Equivalently its noncontact real-space contribution cancels:

    integral_0^infinity (x^2-v^2 tau^2)
                        /(x^2+v^2 tau^2)^2 d tau = 0
    for x!=0.                                           (8)

It must not be bounded by the integral of the absolute value, which would
create an artificial 1/|x| tail. The UV completion gives a local,
finite static response. Density descendants give faster decays; a
dimension-two contribution produces at worst |x|^-3 after the time
integral. The leading harmonic contributes |x|^(1-2K_gap), hence |x|^-7
at half filling. The thermodynamic static kernel and its first spatial
moment are summable. The canonical finite-ring zero mode must be separated
as described next. A single-impurity energy has no infrared-divergent
second-order coefficient.

The forward density coupling may be completed to a square exactly in
the Gaussian action. It changes the local density and produces a finite
quadratic energy gain. It does not flow to a reflecting impurity. The
orthogonality exponent of the corresponding wavefunctions is not an
energy divergence.

[derived; necessary finite-ring refinement] At fixed total vacancy number,
minimizing

    sum_j [ e_gap'' (delta r_j)^2/2 + J s_j alpha'(r) delta r_j ]
    subject to sum_j delta r_j=0

gives the Gaussian susceptibility

    chi_forward,ij = C [delta_ij-1/N],
    C = [alpha'(r)]^2/(2 e_gap''(r)).                     (8a)

Thus the **full finite-ring** susceptibility has a uniform negative
pedestal. Its raw absolute first moment need not be bounded as N grows.
At rho=1/2, C t=(1/4+1/pi^2)^2/(pi/2)=0.07857579749....
The rank-one part contributes +J^2 C (sum_j s_j)^2/N to the energy.
It must be handled explicitly, rather than included in an absolute
first-moment estimate for the local response.

## 3. The linked energy comparison

[proved conditional algebra] Separate the canonical rank-one term (8a)
from the static kernel. Suppose the remaining local kernel obeys
sum_r(1+|r|)|chi_local(r)|<C_0/t. For bounded scalar material sources s_j,
the quadratic charge energy is -J^2 sum_ij s_i chi(i-j)s_j. Changing s_j
on k sites changes this energy by at most

    C' (J^2/t) k,

where C' depends on the source bound, not L. This follows by expanding
the difference into its linear and quadratic terms and summing the
absolute kernel row bound.

For the removed rank-one term, a change of k bounded sources changes
sum_j s_j by O(k). The identity

    [(S+delta S)^2-S^2]/N = 2(S/N)delta S+(delta S)^2/N

and |S|=O(N), |delta S|=O(k), k<=N give an O(k) bound directly.
There is no extensive leftover from this canonical zero mode.

The same estimate, with a larger finite constant, applies to k interfaces
between translated copies of a fixed period-b source pattern. The bulk
energy on each side is identical by translation. Only linked pairs
crossing an interface remain; at separation r there are at most |r|
such pairs per interface. This is why the first spatial moment of chi,
after the rank-one zero mode is separated, is the relevant bound. A
translated period-b source has the same mean; unmatched segment endpoints
and interface cores change its total source by only O(bk), so the separate
rank-one estimate remains proportional to k at fixed b.

[standard continuum inference] After the smooth Gaussian response is
resummed, the remaining contact harmonics and descendants are irrelevant
near K_gap=4. Their connected contributions have finite static interface
energies. Ultraviolet coefficients are finite at fixed microscopic
parameters. Thus the continuum linked expansion predicts an error per
defect r(J/t)=o(J), typically C J^2/t; allowing a factor
1+|log(J/t)| does not change the conclusion. The same expansion controls
a finite density of defects by an error proportional to their number.

If Delta_def is the positive material-sector margin established in
[`README.md`](README.md), the expected dressed comparison is

    E_defect-E_selected
      >= k [J alpha(rho) Delta_def
             -C_b (J^2/t)(1+|log(J/t)|)]                 (9)

up to the usual boundary and zero-mode terms treated in the same linked
comparison. At fixed b=24 every C_b is allowed to be enormous. Choosing
one sufficiently small fixed J/t makes the bracket positive. Nothing
requires a lower bound uniform in b.

## 4. What must be added for the full d=24 chain

[scope] Its source is the operator h_j of the material Hamiltonian,
not a fixed scalar s_j. The charge response produces a retarded
operator-valued influence action. Replacing it by its static value at
all frequencies would be unjustified, since both charge and payload
are gapless.

[physics-level matching requirement] The needed statement is that the
linked UV matching and the remaining low-energy charge-payload theory
preserve the per-defect estimate (9). In implementing it:

1. The massive material excitations and interfaces have fixed positive
   gaps in units of J at the chosen parameters. Their charge dressing
   is controlled by the finite contact response, not a global charge
   spectral gap.
2. The smooth charge-density times payload-density couplings are retained
   as Gaussian forward couplings. Their positivity requires their square
   to be small compared with v_charge v_payload. Even the smallest
   payload velocity is a fixed positive number times J, so sufficiently
   small J/t satisfies this condition.
3. Momentum-allowed non-Gaussian mixed operators must be checked in the
   correct material-coordinate compactification. The large value
   K_gap near four protects the pure contact harmonics; it does not
   replace the separate payload and gluing audit.
4. Sector comparisons subtract the common linked bulk contribution.
   Neither an O(J^2 L) norm bound nor a many-body Schrieffer-Wolff gap
   argument establishes a finite defect energy.

These are standard effective-field-theory steps, and the scalar response
calculation gives positive evidence for the mechanism. A fully uniform
microscopic proof for arbitrary operator-valued h_j is **not** supplied
here. Equation (9) should be labeled a physics-level transfer conclusion
when this matching input and the separate mixed-vertex analysis are
accepted, rather than a proved lattice inequality.

### Simultaneous hierarchy for the full cell Hamiltonian

[physics-level scale audit] The fluctuating contact multiplies the entire
microscopic h_23, including its large machinery terms. This does not
require its norm to be comparable to the tiny payload bandwidth. Fix g
first, then lambda, and let M bound the resulting material interaction
strength. At this stage the machinery gap Delta_*, selector margin
gamma_*, smallest payload velocity v_*, and all momentum margins are
strictly positive fixed numbers, however small. They are measured in
the unscaled internal Hamiltonian H_23.

The leading vacancy-dressed scales are respectively J alpha Delta_*,
J alpha gamma_* and J alpha v_*. A conservative charge matching error
is C J^2 M^2/t, allowing the factor 1+|log(J/t)| if needed. Therefore
one may impose, with a chosen small tolerance eta>0,

    s M^2 (1+|log s|)
        < eta alpha min(Delta_*, gamma_*, v_*, v_* delta_k),
    s=J/t.                                              (10)

All constants on the right are fixed before s is selected. Since
s(1+|log s|) tends to zero, a nonempty positive interval satisfies
(10). There is no choice of g or lambda that depends back on J.
The factor v_* delta_k protects Fermi-point displacement as well as
velocity and stiffness changes. If an ultraviolet charge-momentum
threshold omega_q/t>0 is used in matching, also choose s M much smaller
than omega_q/t. This adds another positive upper bound on s.

Even a deliberately pessimistic forward mixing coefficient of size
C J M has squared ratio to v_charge v_payload bounded by a constant
times s M^2/(alpha v_*). Thus the same hierarchy keeps the full Gaussian
matrix positive. Actual low-energy projections often make this mixing
smaller, but no such improvement is needed for existence.

[standard local spectral estimate] A dimension-one smooth contact
fluctuation has low-frequency spectral density S_A(omega) of order
omega/t^2. A virtual machinery transition of energy J Delta_* therefore
has a ground-state shift bounded at the scaling level by

    J^2 M^2 integral_0^O(t) S_A(omega)/(J Delta_*+omega) d omega
        = O(J^2 M^2/t).

There is no 1/(J Delta_*) infrared divergence. Its on-shell decay width
is of order J^3 M^2 Delta_*/t^2, small compared with J Delta_* when J/t
is small. Higher contact harmonics are less singular. This explains
why the contact's bare operator norm JM is not the appropriate measure
of its low-frequency ability to destroy a machinery gap of order J.
The estimate is part of continuum matching, not a new rigorous spectral
theorem for the full chain.

Exact occupied-number conservation excludes charge-creating vertices;
the original internal U(1) excludes anomalous payload pairing. Remaining
mixed operators are precisely those in the separate material-coordinate
and physical-harmonic audit. Choosing J last changes their coefficients
and stiffnesses continuously inside the already strict stability region.
The extremely small payload velocity makes the allowed J/t extremely
small, but creates no conflict between the inequalities.

## 5. Reproducible finite tests

[numerical] Run

```
python3 contrib/counterexample_campaign/route_a/defect_charge_toy.py --max-L 20
```

The script constructs (1), verifies its exact positive Slater ground
state, solves the projected susceptibility equation, and independently
diagonalizes every perturbed Hamiltonian. There are no fitted central
charges. It tests L=8,10,...,20 at N=L/2, single-defect couplings
J/t=.04,.02,.01, and adjacent pairs, alternating defects and all-bond
defects at J/t=.02. All data are in
[`defect_charge_results.json`](defect_charge_results.json).

| L | Gap Hilbert dimension | chi_00 t | Absolute susceptibility row sum times t | Single-defect cost / (J alpha_L), J/t=.01 |
|---|---:|---:|---:|---:|
| 8 | 35 | .0774769 | .1304168 | .9973008 |
| 12 | 462 | .0848338 | .1446951 | .9971059 |
| 16 | 6435 | .0884536 | .1517832 | .9970038 |
| 20 | 92378 | .0905947 | .1559954 | .9969413 |

At L=20 the single-defect relaxation divided by J^2 is .0899803,
.0902880 and .0904415 as J/t=.04,.02,.01, approaching the independent
susceptibility .0905947. At J/t=.02 the cost divided by its first-order
value is .9942633 for two adjacent defects, .9959832 for alternating
defects, and .9982980 when every material bond is penalized. The charge
state adapts, but its gain is proportional to J^2 per defect in these
tests.

The finite-ring zero-mode refinement is visible directly: N chi_(0,N/2)
approaches -0.0785758, reaching -0.0780119 at N=10. The raw absolute first
moment increases from .0717373 to .1900695 between N=4 and 10. After adding
the predicted C/N pedestal back to the kernel, its first absolute moment
instead decreases from .00683854 to .00636999. These values corroborate
the separation (8a); they do not establish an infinite-volume bound from
small sizes.

Author: Codex GPT-6-family route-A worker, with the coordinating worker
independently deriving the same material-coordinate stiffness, 2026-10-28.
No human review is claimed.
