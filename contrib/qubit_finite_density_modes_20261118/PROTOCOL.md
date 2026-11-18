# Dated protocol: two independent dressed bosons at p153c

Prepared 2026-11-18. **Not launched.** No production run, checkpoint transfer
or cluster submission is authorized or performed by this contribution.
This is a new protocol for the fixed asymmetric point, not an amendment of
the older sym40 production registration.

## 1. Fix the physical point and record the missing normalizations

Use exactly the Hamiltonian in README.md, including every end-site field.
Do not use the nearby `p153_n316` field records as derivatives at p153c:
their different field, finite-parity staircases and smaller chi do not
provide a response matrix at this point. Auxiliary probe derivatives, if
used, are extrapolated to **zero source at the same fixed Hamiltonian**;
they do not authorize retuning its field with length.

The bosonic convention is

\[
b_j\sim\sum_{\nu=A,B}e^{ik_\nu j}\Psi_\nu(j),\quad
\Psi_\nu\sim e^{i\theta_\nu}
       \sum_{m\in\mathbb Z}c_{\nu m}
          e^{2im(\pi\rho_\nu x-\phi_\nu)}.
\]

Different valleys have mutually bosonic fields. A physical Jordan--Wigner
fermion includes the **total** density string and cannot be treated as a
single-valley fermion without changing its charges and momenta.
The dressed densities rho_A,rho_B, centers k_A,k_B, operator residues,
probe matrix Z and matching length a0 all need a measured or controlled
definition. A momentum filter on the raw boson does not by itself normalize
a conserved valley charge. Fix compact charges by consistent density
harmonics and charge/winding insertions or by a controlled matching to the
leading projected theory; verify the emergent continuity/commutator
normalization. Microscopic U(1) must not be imposed on the full chain.

Energy units: the actual Pauli Hamiltonian has Delta=1, V=40 and
|t2|=1/80. Momenta and distances always use **physical sites**, not a blocked
two-site lattice. An effective-unit velocity must be multiplied by 1/80.
The matching energy cutoff and a0 are currently **unknown**; neither one
physical site nor the largest fitted separation is silently assigned to a0.
Vary a proposed matching cutoff within a demonstrated separation of scales
and include that variation in the error. Boundary fixed points are a separate
measurement, not assumed from a parity-gap formula.

## 2. Minimum new measurements

The existing saved records cannot supply these by further algebra. Two
complementary data groups are needed at the same point:

1. **Long complex correlations and harmonic identification.** Save full
   one-body, connected density and connected physical pair correlations for
   d=2 and d=3, both pair means, and several central anchors. Resolve both
   single-boson carriers, their sum in the pair channel, and the independently
   identified density harmonics at 2pi rho_A and 2pi rho_B. Test alternate
   harmonic assignments and unequal populations. The raw spectral weights
   and `n_mid/2` are not occupations. Positive dressed rho_A and rho_B must
   be separated from zero within their errors; two peaks alone do not pass.
2. **Four independent matched response channels.** Measure the low-energy
   density/current spectral-response tensor, including its cross terms.
   A concrete starting observable set is the SW-dressed versions of
   `{N,T1_real,T1_imag,J2}` from RESPONSE.md, using central smooth probes.
   Determine its in-medium Z, high-energy analytic contribution and smallest
   singular value. The four old leading parameter directions fail the rank
   test. If this replacement also becomes singular in the medium, add a
   distinct probe and identify which field remains invisible; do not invert
   a pseudoinverse and declare four normalized coordinates.

Four independent source coordinates and their ten real symmetric static
response entries determine the ten entries of H after matching. The dynamic
matrix provides an essential independent check: resolve two positive and
two negative linearly dispersing Gaussian branches, with appropriate
residues in both independent canonical pairs, and compare their speeds to
eta H. One physical probe may miss a mode through a small residue. A large
entropy slope, two one-body carriers, or two parity sectors cannot replace
this test. A possible boson-plus-Ising alternative requires its own operator
content and cannot be counted as a second canonical boson.

If dynamic response is inaccessible, a verified thermodynamic
density/winding energy Hessian can replace the static kernel. Independent
twists and number-sector changes must refer to **matched emergent charges**;
raw parity gaps or nonexistent exact microscopic number sectors do not
supply them. Such a Hessian is presently missing.

## 3. Bulk windows and convergence gates

Keep the Hamiltonian, operator, filter and physical separations fixed when
comparing chi; separately grow the distance window with L. Save enough
central rows that anchor changes can be tested without moving toward an
end. As initial *estimator* gates require at least 24 distances, a distance
span of at least two, and distance from every operator support to an end
at least twice r_max. These are necessary practical gates, not proofs of
asymptotic scaling. For example a centered 32--96 window can fit this
geometry at L=512; a 48--160 window needs a larger chain. Neither example
is a request to launch those runs.

Compare at least three increasing chi values on the same state branch and
windows; include both parities when identifying the finite-chain ground
branch. Report entropy, energy, density and the **actual measured response
and phase** changes. Require the measured correlation-length cutoff to lie
well beyond the window and test its effect explicitly. A small local
eigensolver residual is not a convergence test of long correlations.

For carrier phases, vary anchors, pair separation, filter width, seeds,
distance windows and subleading-harmonic models. The XX calibration in
DATA.md quantifies a known short-window bias; the nearby-harmonic test shows
why formal least-squares errors are not enough. Include finite-chi, finite-L,
filter, operator-matching and harmonic-assignment uncertainty separately.
Initially resolving quartet momentum to substantially better than 0.001
radian/site would distinguish the observed size drift from an apparent
precision artifact; the final required accuracy is set by the locking
comparison, not by that preliminary number.

For static response use several small nonzero q values in a large central
window, with several wavelengths per window. Remove independently measured
analytic high-energy terms, compare the static spectral integral to source
derivatives if both are available, and extrapolate source strength to zero.
Propagate uncertainties in Z and the contact before inverting. Then require
the inferred positivity margin of H to exceed its error. Report both
fixed-winding and relaxed Schur complements, including their relative-sector
projections. Evaluate G and signed speeds with the **symmetric chiral metric**.
Check `G eta G=eta` and independent dynamic speeds; do not enforce this
identity on an inadequate fit merely to make it look Gaussian.

An alternative equal-time determination of G needs enough independent
known vertex overlaps or four normalized smooth-probe slopes. For instance,
the ten quadratic forms for the basis
`{theta_A,theta_B,2phi_A,2phi_B}` and its six pairwise sums recover a general
symmetric G. Observables of odd parity can be correlated even though
Hamiltonian perturbations have even parity. Each harmonic assignment must
be verified. Equal-time data alone do not determine velocity scales.

## 4. Rederive the operator lattice off the symmetry line

The microscopic actions T (translation), P (parity) and Theta=I K survive.
The valley-exchanging antiunitary A and unitary Sigma of the pure-imaginary
hopping line do **not** survive Re w nonzero. In a compatible field origin,

```
T: phi_nu(x) -> phi_nu(x+1)-pi rho_nu,
   theta_nu(x) -> theta_nu(x+1)+k_nu
P: theta_nu -> theta_nu+pi
Theta: phi_nu(x) -> -phi_nu(-x), theta_nu(x) -> -theta_nu(-x)
        (antiunitary).
```

Therefore there is no forced balance, no exact `k_A+k_B=pi`, and no Sigma
orbit cancellation. Uniform even-separation pair sources, forbidden by
Sigma on the symmetric line, are now permitted. Their finite-density
amplitudes and perturbative orders belong to assignment A; symmetry
permission is not a measurement of an amplitude.

For **every** integer n=(n_A,n_B), m=(m_A,m_B), the local bosonic vertex is

\[
V_{n,m}=e^{i(n\cdot\theta-2m\cdot\phi)},\quad
\ell=(-2m_A,-2m_B,n_A,n_B),
\quad q=n\cdot k+2\pi m\cdot\rho\pmod {2\pi}.                 \tag{1}
\]

Hamiltonian vertices require `n_A+n_B` even. A strictly uniform resonance
has q=0; a detuned near-resonance must still be checked for a locking tongue.
Theta allows real cosine coefficients in this field convention, with no
additional A/Sigma projection. The dimension is `x=ell^T G ell/4`.

| Process | n | m | Actual momentum condition |
|---|---|---|---|
| Intra-valley pair A or B | (2,0) or (0,2) | (0,0) | 2k_A or 2k_B |
| Inter-valley pair | (1,1) | (0,0) | k_A+k_B |
| Quartet | (2,2) | (0,0) | 2(k_A+k_B) |
| Relative density backscattering | (0,0) | (1,-1) | 2pi(rho_A-rho_B) |
| Valley conversion | (-1,1) | (0,0) | k_B-k_A |
| Mixed quartet/backscattering | (2r,2r) | (s,-s) | 2r(k_A+k_B)+2pi s(rho_A-rho_B) |
| All pair-, conversion- or charge-harmonic dressings | any parity-even n | any integer m | full equation (1) |

There is generically no exactly resonant nonzero integer vertex when all
four dressed momenta/densities are independent incommensurate numbers.
That generic mathematical fact does not bound the detuning tongues of this
specific Hamiltonian. In particular an imbalance can detune backscattering,
but interactions could lock the densities; a raw valley peak does not
resolve which case applies.

The inter-valley pair exponent, *if that leading overlap is established*,
still satisfies `x_quartet=2 eta_pair` without balance or equal velocities.
This single relation determines neither the mixed dimensions nor G.

For ordinary simultaneous pinning also impose
`ell_i^T eta ell_j=0`, including i=j. Equivalently
`m_i.n_j+n_i.m_j=0`. Non-null relevant perturbations are retained in the
lattice test but flagged: weak relevance is not by itself proof that such
an operator simply gaps a boson. Derivative descendants raise dimension;
linear gradients change the equilibrium rather than pin fields.

## 5. Exhaustive conditional test and handoff

If a measured covariance has `||delta G||<=e_G` and
`lambda_lower=lambda_min(G)-e_G>0`, every potentially relevant vertex
obeys

\[
|n|^2+4|m|^2\le 8/\lambda_{\rm lower}.                         \tag{2}
\]

`response.dangerous_vertices` uses this proven finite radius, not an assumed
low-charge truncation. It retains parity-even representatives up to
Hermitian conjugation, mixed and density-assisted vertices, and self-null
flags. For each it reports

```
delta x <= e_G*(|n|^2+4|m|^2)/4
delta q <= sum |n_nu| delta k_nu + 2pi sum |m_nu| delta rho_nu
|q|_lower = max(0, |wrap(q)|-delta q).
```

At x=2 use the marginal flow, not a strict irrelevance label. Even with a
positive mismatch lower bound, compare the running amplitude and dimension
to detuning: schematically `a0 |y|^(-1/(2-x))` versus `1/|q|`.
The dimensionless y, velocities and a0 must be those of assignment A's
**same** matched finite-density convention. Their values are not supplied
by a vacuum S-matrix or an order-one amplitude guess.

The current handoff in interface.json deliberately contains no fitted G,
rho_A/B, velocities, matched momentum uncertainty or a0. It supplies the
fixed physical point, conventions, audited finite-window diagnostics and
the exact missing-response definition. A stable region or an instability
can be returned only after the uncertainty-aware inequalities decide it.
At present the scientifically justified result is **undecided**, with a
new density-dressing calculation, a rank obstruction, and a calibrated
route to the missing quantities.
