# Prioritized missing inputs for the fixed qubit infrared theory

This list is part of the requested first investigation, not authorization
to launch production. Keep the literal Hamiltonian in `README.md` fixed.
Its asymptotic criticality is **not established**. The priorities below
distinguish identification inputs from optional refinements. An apparent
entropy coefficient near two is evidence to explain, not a prerequisite
that subsequent measurements are required to reproduce.

All continuum normalizations below are conditional on a two-component
Gaussian description. Use physical site spacing and hbar equal to one.
In the uniform pairing-plus-DM frame, define
`b_j=(X_j+iY_j)/2`, `n_j=(1-Z_j)/2`,
`P_d(j)=b_j b_(j+d)`, and
`kappa_2(j)=X_j Y_(j+2)-Y_j X_(j+2)`.
The uniform rotation does not change momenta. The additional site-dependent
real computational gauge shifts a pair momentum by pi modulo 2pi.
Undo its phases explicitly before comparing different operators. Never
import momenta or fitted parameters from p153c or delta=1/8.

## Necessary identification inputs, in priority order

### 1. Recover and identify the pair carrier before interpreting its exponent

**Quantity and limit.** Save the full complex connected functions
`D_dd'(i,j)=<P_d(i)^dagger P_d'(j)>-<P_d(i)^dagger><P_d'(j)>`,
with d,d'=2,3, the anomalous connected functions
`<P_d(i) P_d'(j)>_c`, their adjoints and one-point functions, one-body
`<b_i^dagger b_j>`, and the physical-frame Fourier momenta. Equivalently
measure the complete matrix of both Hermitian pair quadratures directly.
The normal pair function alone cannot separate those quadratures when
microscopic U(1) is absent. For the leading identified pi primary extract its bulk
equal-time exponent `D_pi(r) ~ r^(-eta_pi)`, taking bond-dimension
convergence and then `1 << r << distance to either boundary`, with windows
moved outward as L grows. A maximum found by optimizing q is not a
measurement at a predetermined q=pi. A finite-window maximum is not a
proof of asymptotic momentum or of an undressed primary.

**Known.** At L=384 the saved separation-three fit is 0.9193(225);
separation two gives 0.9972 without a saved error. The estimator freely
selects q, but the printed target table omits it. Target complex arrays
and states are absent from the accessible checkout. Small rings show
substantial pi weight in both pair operators, while the global
separation-three peak shifts away from pi at L=10,12. These observations
do not identify the large-L fitted carrier.

**Discrimination.** In the generic balanced bosonic charge lattice, a
genuine pi primary with exponent below one has an allowed doubled
zero-momentum vertex of dimension `x_double=2 eta_pi<2`.
It is specifically the primitive charge quartet only after matching the
primary to `exp[i(theta_A+theta_B)]`. An exponent above one removes this
particular weak-coupling threat, but not every other lock. A different
carrier/descendant invalidates this test. Coinciding exponents alone do
not prove matching: check momenta, symmetry, cross-correlations and
compatible amplitudes across the two separations/quadratures.

**Smallest next action and resources.** Recover the already-generated
measurement JSON and estimator's fitted q values, if their original
owner still has them. C193 fits one carrier per pair separation across
the input rows, shared by its two fit windows; it does not already
provide independent per-window momenta. Refit moving windows separately.
At L=384 its existing windows span r=39–115 and 58–145, with four
central anchors; keep both the freely fitted q and a predeclared
physical pi fit (q_rot=0), and track resolution and residuals.
Reanalysis should take seconds to minutes and no
new many-body optimization. If only the old MPS states survive, contract
the selected pair/one-body matrix on the highest available chi=1024
L=192,256,384 states, checking lower/higher-chi versions if recoverable;
the L=384 entropy still moved by 0.017 and the printed exponents have
only one chi per length.
cost is unknown without states and the contraction layout, and should be
timed on one window first. If neither survives, new state production is
needed and requires a new instruction. No summary-table operation can
recover the missing complex data.

The old target JSON contains same-separation connected normal pairs,
one-body and density functions. It does not contain cross-separation
or anomalous pair rows; its `mm` and `pp` entries are dummy zeros, not
recoverable pair means. Recovering that JSON enables the carrier refit,
while the complete new quadrature matrix still needs states and fresh
contractions.

### 2. Establish two bulk acoustic modes and their velocities

**Quantity and limit.** At physical momenta q near zero measure a matrix
of connected spectral functions of local Hermitian probes, using a
consistent convention, for example

`O_a(q)=L^(-1/2) sum_j exp(-iq j) [O_a(j)-<O_a(j)>]`,

`S_ab(q,omega)=sum_(n!=0) <0|O_a(q)^dagger|n><n|O_b(q)|0>`
`                 * delta[omega-(E_n-E_0)]`.

Candidate probes include n, the real and imaginary nearest-neighbour
hopping quadratures, and the uniform distance-two current/chirality.
The **uniform** distance-two chirality may overlap a relative-density
coordinate; its **staggered** counterpart has different translation
quantum numbers. Use both q signs, a bulk thermodynamic limit before
q tends to zero, and spectral weight rather than just the lowest level.
If using open chains, normalize the central-window envelope and fit
several wavelengths within it. Identify derivative/current channels
before counting their linear branches or taking
`v_i = lim_(q->0) [E_i(q)-E_0]/|q|`.

**Known.** No large-chain target dynamic response, matched bulk current
spectrum or thermodynamic velocity is saved. The new L=8,10,12 exact spectra have
changing ground-state sectors, unresolved one-body valleys, and several
neutral poles exchanging weight. They supply operator checks, not a
velocity extrapolation. Antiunitary symmetry can duplicate momenta
without adding a mode. Pair, bond and chirality reaching the same pi
pole likewise does not establish three modes.

**Discrimination.** Two independently normalized Gaussian current
oscillator families with linear dispersions support two bosons. One
surviving family plus a positive threshold in an independently isolated
second channel supports locking/crossover to one mode. Subtract the
surviving mode and its harmonics first: it still has gapless neutral
excitations. Other critical operator content remains possible if there
is no two-current algebra. Two established Gaussian sectors with distinct
physical speeds give two Gaussian oscillator sectors with separate speeds;
their zero modes may still obey shared constraints. This does not establish
one common-velocity conformal theory. Equal speeds need a separate convergence test.
Multiple poles, multiparticle harmonics, vacuum doublets, or two peaks
in one structure factor are not by themselves independent modes.
Test each sign of q separately: summing left/right residues first can
give a false extra rank. At coincident speeds a higher-rank residue is
not enough either; nonabelian currents can also share a critical sector.
Two independent abelian oscillator/current algebras and consistent
zero modes are required for the two-boson conclusion.

**Smallest next action and resources.** First retrieve any existing
operator-resolved transfer or excitation data with original momenta and
normalizations. Unfold the real computational frame's two-site cell
using physical translation and both sublattice form factors.
Transfer-matrix correlation lengths alone do not give
velocities. Otherwise a targeted excitation or real/imaginary-time
calculation on converged target states is needed. Start with the above
four probes and at least three decreasing resolved |q| of both signs,
at two or more converged chi values with |q| xi_min much greater than one.
Use a pilot to estimate cost. Stop if residues or speed separation are
below frequency resolution or bond-dimension uncertainty. Production
cost is currently unknown. Enlarging
tiny rings alone will not resolve scales of hundreds of sites.

### 3. Determine the Gaussian metric through matched responses, not two fitted K values

**Quantity and normalization.** In an independently established bosonic
frame use `Phi=(phi_A,phi_B,theta_A,theta_B)`,
`[Phi_i(x),Phi_j(y)]=(i pi/2) eta_ij sgn(y-x)`,
`eta=[[0,I],[I,0]]`, and
`H_G=(1/(2pi)) integral (partial_x Phi)^T H (partial_x Phi) dx`.
The symmetric spatial log covariance is
`G=H^(-1/2) |H^(1/2) eta H^(1/2)| H^(-1/2)`;
`x_l=l^T G l/4` and `G eta G=eta`.
For densities and phase gradients
`X=diag(-1,-1,1,1) partial_x Phi/pi`, the matched static response is
`chi_X=(1/pi) D H^(-1) D`, with `D=diag(-1,-1,1,1)`.
Phase gradients are not automatically physical currents.

Measure the bulk, two-sided imaginary-time kernel
`chi^O_ab(q)=integral d tau <T O_a(q,tau) O_b(-q,0)>_c`
and independently match `O=Z X+...`, removing or bounding analytic
high-energy backgrounds before inverting. Include both momentum signs
and time orderings. Do not set the matching matrix Z to identity by
definition, invert a finite-ring conserved-number zero response, or
call an arbitrary spin twist an emergent-U(1) flux. Raw N is not conserved
by this microscopic Hamiltonian. Do not use H itself as an independent
uniform response direction: its zero-temperature response is null.

**Known and discriminating outcomes.** If the fitted pair and one-body
functions are the proposed primaries, two combinations of G are
constrained, leaving an unknown mixing parameter even after symmetries.
The simultaneous unit-normalized density assignment fails the necessary
inequality at the quoted central values:
`K_density * eta_pair = .675 * .9193 = .62053 < 1`.
Mixing alone cannot repair this. Distinguish density dressing/background,
wrong carrier, finite-window error, and failure of the Gaussian ansatz.
A full matched positive response identifies H, its signed propagation
speeds, and G, and tests all dimensions without redefining a parameter
from the same correlation used to validate it. Rank and normalization
checks are needed before declaring a missing mode.

**Smallest next action and resources.** Joint spectral/response data from
item 2 can be reused. First derive or calibrate low-energy overlaps and
conditioning of the probe matrix; four probes are a starting set, not
a guaranteed full-rank in-medium basis. Static or imaginary-time
contractions need states/new production; no such target kernel is
available. Small 4x4 inverse/covariance algebra is negligible once a
well-conditioned matched kernel exists. Its many-body measurement cost
is unknown and requires a pilot, not a p153c-derived cost or parameter.

### 4. Establish the compactification and local charge lattice independently

**Quantity and limit.** After identifying branches/velocities, organize
periodic bulk levels by physical translation, exact parity and the
antiunitary/unitary symmetry actions, with matrix elements of b, both
pair quadratures, density, bond and chirality. Test zero-mode dimensions,
momenta and descendant multiplicities in several sizes and commensuration
classes. In the conditional bosonic convention the proposed local
vertices are `exp[i(n dot theta-2m dot phi)]`, integer n,m; their momentum
is `n dot k+2pi m dot rho` modulo 2pi. This is a hypothesis to test, not
a lattice proven by observing c close to two. A Jordan-Wigner fermion
includes a density string and cannot be assigned this local boson
lattice without accounting for that string.

**Known and discrimination.** Exact microscopic parity/translation
actions are known; emergent separate charges, dressed densities and
allowed winding sectors are not. A torus-type Gaussian spectrum requires
a compatible local vertex lattice and two current algebras, not just
two fitted exponents. Twisted sectors, missing vertices or additional
primaries can distinguish a different compactification or an orbifold;
any named rational theory needs its specific spectrum/fusion or current
enhancement, not total c. Avoid treating OBC edge/vacuum splittings as
bulk primary levels. A lattice symmetry twist is meaningful only after
its seam and the symmetry being twisted are specified microscopically.

**Smallest next action and resources.** Reuse item 2's operator-resolved
periodic spectra if large enough. The current tiny-ring data are a
convention control and cannot supply an asymptotic tower. With unequal
velocities, excitation energies weight the different chiral dimensions
by their respective speeds; do not convert every gap with one velocity.
The exact parity defect supplies one concrete quantized seam: reverse
the XY and YX terms across that seam and keep ZZ unchanged. Its spectrum
is a diagnostic, not a K estimator by itself. Targeted
low-lying excitations or symmetry-seam spectra at substantially larger
sizes require new production; cost is unknown. Defer flux/winding
curvatures until the associated conserved emergent charge and matching
are independently established.

### 5. Test all-scale stability and the allowed locking amplitude

**Quantity and limit.** Measure independently identified pi pair and
staggered bond/chirality correlations, order-parameter distributions
and connected subtraction, alongside both current-mode gaps and
entropy, over increasing converged length scales. In a declared
two-boson matching scheme determine or bound the renormalized
coefficient of `cos[2(theta_A+theta_B)]`; also test special momentum
resonances using actual carrier momenta and dressed species densities.
A matched amplitude requires a stated cutoff and field normalization;
the presence of an allowed microscopic monomial does not set its
continuum coefficient. Preserve the fixed Hamiltonian: probe insertions
or response derivatives are diagnostics, not a parameter deformation
whose phase is substituted for the target.

**Known and discrimination.** Under the two undressed primary assignments,
`p=eta_pair=.9193`, `d=4eta_b-eta_pair=.5083`. For arbitrary
symmetry-allowed mixing t the generic dimensions obey
`x_(r,s)=2[p(r-t s)^2+s^2/d]`.
Thus the primitive quartet has x=1.8386 at central values, while
all s!=0 have x>=2/d, approximately 3.9347; special resonances remain outside this generic
test. The quartet is not symmetry-forbidden, but its target continuum
coefficient is unknown. Relevant nonzero coupling can give a very long
crossover; a true stable two-boson point would require irrelevance,
vanishing coupling on an invariant critical manifold with an explanation,
or failure of the asserted
primary/momentum assignment. One locked combination can leave c=1;
additional independent locking can gap both. Pi order must be analyzed
as translation breaking, not relabeled smooth valley polarization.
A single bare or finite-order matching coefficient equal to zero does
not establish all-scale stability: nonlinear generation and other
allowed vertices must also be controlled.

**A quantitative conditional amplitude test.** Once an independently
normalized *pure emergent annihilation primary* O at pi has been
identified, measure
`R(r)=<O(r)O(0)>_c / <O^dagger(r)O(0)>_c`.
It vanishes in the unperturbed charge-conserving Gaussian theory.
At first order in a dimensionless matched quartet coupling at cutoff a,
`R(r) ~ C g4 (r/a)^(2-x4)`, with C fixed by the normalized Gaussian
three-point integral and velocities. At the proposed exponent this
power is about 0.1614. The simple integral converges for
`1/4<x_O<1/2`, with `x4=4 x_O`, finite nonzero velocities and a window
where the perturbation is still small. The auditor derives the
homogeneity and convergence conditions in [AUDIT.md](audit/AUDIT.md).
A physical P generally mixes O and O^dagger; its uncalibrated anomalous
correlator can therefore be nonzero even at g4=0. Other matching
couplings and running parameters must also be separated. This protocol
does not extract g4 from the existing normal pair data, and a null
signal does not prove an all-scale zero coupling. It can reuse the
quadrature contractions from item 1 after charge/mixing identification.

**Smallest next action and stopping criteria.** Resolve item 1 before
allocating a long entropy run. Reuse converged states first, compare at
least two moving bulk windows and anchors, and retain the actual raw
data and convergence changes. If new production is authorized later,
start with a single L=384 convergence pilot; proceed to L=512 only when
window entropy changes <0.005 and diagnostic exponents/momenta are
stable under a further bond-dimension step to stated tolerance (for
example 0.02 for the pair exponent, a target precision not guaranteed
accuracy). Stop a pre-budgeted pilot on failure to converge or changing
state branch; do not count its apparent plateau as infrared evidence.
No finite maximum length excludes an arbitrarily long locking scale
when the magnitude of its coupling has no positive lower bound.
Near the quoted dimension,
the linearized coupling grows by only about 1.25 from L=96 to 384.

**Resources.** Existing documentation quotes 1–6 CPU-hours for a
different, chi<=256 iMPS ladder; it is not a cost estimate for converged
L=512 target states or dynamics. A useful production budget requires
pilot timing, peak memory and contraction timing. No production is
launched by this investigation.

## Optional refinements after those identification gates

- Improve entropy/profile systematics, including moving windows and
  oscillatory corrections. This helps establish scaling range, but
  cannot substitute for mode and charge-lattice information.
- Add finite-temperature entropy/energy and an independently normalized
  Casimir coefficient as cross-checks of the already measured velocities
  and central charges. Do not infer both c and v from one product.
- Test equality of the two speeds or enhanced current symmetry only if
  supported by the mode data; neither is required for two bosonic
  contributions to the entropy coefficient.
- Investigate nearby parameters or lower-cost controls only as separately
  labeled calibration/robustness work. Such results cannot identify this
  fixed target by continuity without demonstrating the intervening phase.

## Audit status

The independent AI auditor reviewed the list and accepted it after
corrections to anomalous-pair requirements, the estimator's shared
momentum, matching/current-rank tests, and all-scale stability limits.
See [AUDIT.md](audit/AUDIT.md) and [PROGRESS.md](PROGRESS.md). This is
internal AI review, not independent human review.
