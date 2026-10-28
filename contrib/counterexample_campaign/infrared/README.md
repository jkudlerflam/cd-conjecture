# Vacancy completion of the 23-state candidate: infrared audit

**Status.** [proved, conditional Gaussian statement] A species-blind vacancy
lift has a simpler resonance test than the existing two-length wrapper.
At occupied density 1/2, every uniform non-Gaussian perturbation of the
reference theory has scaling dimension at least three, provided the
squeezed 23-state payload has C33's simple roots and noncollision margins.
This statement includes the fractional harmonics of its period-24 order.
It does not assume that the mobile charge is an independently compactified
free fermion.

[heuristic, microscopic conclusion] The proposed 24-state chain should
have c=24 at sufficiently separated, fixed, nonzero couplings, using
ordinary squeezed-space matching and weak-coupling Gaussian RG. This
audit does **not** establish the missing all-order, low-frequency
microscopic matching theorem. It also does not prove unique ground states
on every ring. The exact finite facts, continuum inference, and remaining
inputs are separated below. No measured c=24 is claimed.

Author: OpenAI Codex, GPT-6 family, 2026-10-28. Independent continuation
within `counterexample_campaign`. The project files read include
`explicit_candidate`, C13, C33, C36, and C37. The C33 root intervals are
independently processed here using exact rational arithmetic; their
matrix-inertia proofs are inherited, not rerun here. No human review.

## 1. Definition and exact squeezed coordinates

Let h be the 23-state bond in `contrib/explicit_candidate/bond23.py` at
the n=23 parameters certified in C33, and b=n+1=24. Add one orthogonal
vacancy state. Write n for the occupied-site projector and lift h to
zero whenever either site is vacant. Consider

    H = -t sum_(x,a) (|0,a><a,0| + |a,0><0,a|)
        + J sum_x h_lift(x,x+1) - mu sum_x n_x,
    a=1,...,23; t>0, J>0.

Every parameter is fixed before taking the ring length L to infinity.
This is one Hermitian nearest-neighbor bond on C^24. Occupied number,
the internal U(1) charge, and marker number are exactly conserved.

[proved] For fixed occupied number N, order the occupied sites
x_1<...<x_N and retain their cyclic internal word. A vacancy hop moves
an occupied site without changing this cyclic word. A lifted h acts on
consecutive internal labels only when their physical sites are adjacent.
Away from the chosen ring cut the exact representation is

    H = H_charge + J sum_j A_j(x) h_(j,j+1),
    A_j(x) = 1_(x_(j+1)=x_j+1).

Crossing the ring cut cyclically rotates the internal word. In a sector
of internal translation eigenvalue exp(iK), the charge wavefunction has
the corresponding boundary twist, together with the conventional
spinless-fermion parity factor. There are not b independent copies of
the physical Hilbert space labelled by an arbitrary choice of first
occupied site.

[proved at J=0] The charge motion is the hard-core hopping chain. At
occupied density rho its thermodynamic Fermi-sea values are

    v_charge=2t sin(pi rho),       K_particle=1,
    D(rho)=<n_x n_(x+1)>=rho^2-[sin(pi rho)/pi]^2.

The cyclic average of the N contact indicators gives

    alpha(rho)=(1/N)<sum_j A_j>=D(rho)/rho
              =rho-sin(pi rho)^2/(pi^2 rho).

After projection in a fixed internal cyclic-translation sector, the
first derivative in J is J alpha(rho) H_h. Individual rooted A_j need
not have the same expectation on a finite ring; the cyclic-sector
projection is essential. At rho=1/2,

    alpha=1/2-2/pi^2 > 0.

This finite-order projection identity does not by itself justify
eliminating the gapless charge sector at fixed J in the thermodynamic
limit. Section 6 states the matching input that is actually used.

## 2. The fractional compactification must be retained

[assumed, internal phase] Below its marker/length defect gap the squeezed
internal chain has one marker per b occupied sites, together with n=23
gapless fermionic bands. Let y count cells, and let phi_p count occupied
particles in physical coordinates:

    N(x)=rho x-phi_p(x)/pi,       y(x)=N(x)/b.

Consequently the marker counting field is phi_m=phi_p/b and its
reference parameter is K_marker=1/b^2. Marker order harmonics such as
exp(2i phi_p/b) are local physical observables. They cannot be removed
from the operator lattice merely because occupied number is conserved.

[proved under this compact-boson dictionary] A neutral local payload
vertex with cell momentum q and free payload dimension Delta_pay can
be dressed by any integer cell harmonic ell. The marker comb and its
payload phase give

    sum_j delta(y-j) exp(i q y)
       = sum_(ell in Z) exp[i(q+2pi ell)y].

After substituting y(x), its physical wavevector and scaling dimension
at the decoupled reference are

    k_phys = rho(q+2pi ell)/b mod 2pi,
    Delta_ref = Delta_pay + [(q+2pi ell)/(2pi b)]^2.       (1)

The second term is the charge contribution, not an independent
integer-harmonic charge contribution. Replacing q by q+2pi r simply
relabels ell and leaves (1) unchanged. Derivative descendants add their
usual positive integer dimensions. Exact occupied-number conservation
forbids dual-charge vertices; exact internal U(1) forbids anomalous
payload pairing. Gapped relative N-bM defects are not new gapless modes.

## 3. Half occupancy removes all new relevant resonances

[proved under (1)] A uniform physical Hamiltonian vertex must satisfy

    rho(q+2pi ell)/b=2pi h,       h in Z.

At such a resonance, (1) becomes

    Delta_ref,res = Delta_pay + h^2/rho^2.               (2)

At rho=1/2 the charge contribution is 4h^2. Therefore:

* If h is nonzero, the vertex has dimension at least four even when
  its payload factor is the period-b order parameter, of dimension zero.
* If h=0, then q=0 modulo 2pi. The selection problem is exactly the
  original squeezed payload momentum problem.

[proved from the inherited root certificate and usual fermion
bosonization] In the n=23 payload, every nontrivial neutral quadratic
operator has nonzero momentum; every nontrivial neutral quartic operator
has nonzero momentum, except forward-scattering currents and their
products. C33 checks all unordered pairs with replacement among the 46
Fermi roots, so Cooper transfer and repeated-branch possibilities are
included. Forward terms are Gaussian. Higher neutral fermionic vertices
have free dimension at least three; derivatives cannot reduce it.

Consequently every **resonant non-Gaussian** vertex has reference
dimension at least three. This is the desired stability budget. In
particular, half occupancy does not create a relevant charge umklapp:
its lowest allowed uniform charge harmonic has dimension four, even
though individual, oscillating marker harmonics have dimension 1/b^2.

The same resonance argument works whenever 0<rho<1/sqrt(2), with a
strict margin. For the especially simple uniform dimension-three
classification used next, rho<1/sqrt(3) is sufficient. Half occupancy
satisfies both inequalities.

## 4. Quantitative momentum margin, including mixed harmonics

Let delta_int be the smaller of the distinct-root and distinct-pair-sum
margins. C33's rational certificate gives

    delta_2>0.0140193,       delta_4>3.1558e-7.

[proved under (1)] For every nonconstant vertex with Delta_ref<3,

    |q+2pi ell|/(2pi b)<sqrt(3),
    |k_phys|<2pi rho sqrt(3).

At rho=1/2 the latter is strictly less than 2pi. The distance to a
nonzero reciprocal vector is at least pi(2-sqrt(3)). The distance to
zero is at least (rho/b) delta_int for a nonforward quadratic or
quartic payload, or 2pi rho/b for a nonzero pure/order harmonic. Thus
all non-Gaussian vertices of reference dimension below three have

    delta_phys >= min(delta_int/(2b), pi(2-sqrt(3)))
               > 6.5745e-9,                 b=24.       (3)

Descendants obey the same bound when they belong below dimension three.
For a forward payload current, ell=0 is a Gaussian term, while ell!=0
has the larger pure-harmonic mismatch. This exhausts the potentially
dangerous classes. There is no need for a new Diophantine search or for
an irrational density to exclude the low-order mixed terms.

The test remains finite even though ell runs over all integers: a
dimension cutoff bounds |ell+q/(2pi)| by b sqrt(3-Delta_pay). For a
fixed b,n only finitely many integer chiral charges and harmonics occur
below that cutoff. Very high harmonics stay irrelevant by positivity
of the full dimension form.

## 5. Gaussian mixing and the hierarchy of small parameters

[proved Gaussian statement] In the independent particle/payload
coordinates, let G_ref be the full positive Gaussian energy matrix,
including all unequal velocities. If

    eta=||G_ref^(-1/2)(G-G_ref)G_ref^(-1/2)||<5/13,

the geometric-mean dimension-metric inequality already proved in
`digit_compilation/audit/GAUSSIAN_CERTIFICATE.md` gives

    Delta_G >= sqrt((1-eta)/(1+eta)) Delta_ref
            > (2/3) Delta_ref.

All resonant non-Gaussian vertices therefore stay irrelevant. This
condition includes density-current terms and particle-payload forward
scattering. Comparing a correction only with the fastest velocity is
insufficient; the normalized matrix condition is the relevant one.

[heuristic, RG with fixed finite n] The many payload velocities can be
very different. Choose a matching energy below the smallest linear
window and below each oscillatory cutoff. The RG interval from the
largest bandwidth to this matching energy is long but finite, since
all 46 roots are simple and all mismatches in (3) are positive. For
each fixed n, a sufficiently small initial interaction remains weak
over that finite interval. Below it, the only surviving marginal terms
are Gaussian forward scattering, and the remaining non-Gaussian
vertices have the strict dimension margin above. Slow bands do not
require a uniform-in-n crossover bound.

There are three independently reducible ratios:

1. Keep C33's finite D and loop fixed, and decrease g from its certified
   upper value if needed. Its normalized leading band geometry remains
   fixed while the nonquadratic order-lambda^2 term, relative to the
   leading lambda^2 g^2 hopping, tends to zero as g^2. The root brackets
   persist because the finite-g error bound only improves.
2. Fix that nonzero g, then decrease lambda. The next connected frozen
   core corrections start beyond order lambda^2 and can be made small
   relative to the slowest payload scale, subject to local low-frequency
   elimination of the gapped internal machinery.
3. Fix h with those finite parameters, then decrease J/t. High-momentum
   vacancy corrections are small relative even to J times the slowest
   internal scale. Keep the low-momentum charge boson explicitly; it
   supplies an extra Gaussian mode and its forward couplings. Tune mu
   so the dressed density is 1/2, or remains in a small interval around it.

The smallness conditions contain the finite channel multiplicities,
inverse minimum velocity, and logarithms of the mismatches. Sending
lambda to zero alone does not suppress the g^4 interaction relative to
the g^2 payload, and sending J to zero alone does not repair an unstable
internal h. The above order of choices avoids both mistakes. No
numerical positive J or lambda satisfying every microscopic matching
constant is certified by this audit.

## 6. What the fixed-coupling microscopic argument still needs

[proved limitation] A finite-ring projection onto one free charge
ground state is not a uniform Schrieffer-Wolff construction: its charge
gap vanishes as 1/L. It is therefore incorrect to conclude that the
whole remainder is O(J^2/t) merely from second-order denominators of
a fixed ring. C36 also disproves an unrestricted operator-norm decay
claim for all retained energies in the internal massive-sector step.

[assumed, sufficient physics input] One needs a Wilsonian
squeezed-space matching in which high-energy charge fluctuations and
the gapped internal machinery are eliminated locally, but the
low-momentum charge boson and all payload bands remain. At each fixed
h, corrections to the normalized low-frequency action must tend to
zero with J/t; separately, the internal corrections must have the
small-g and small-lambda behavior described above. The matching must
preserve exact charges, the harmonic lattice (1), and the positive
defect cost. A logarithm of J/t multiplying J/t is harmless; a
nonvanishing normalized correction or new soft defect would not be.

[proved free-charge input] Different internal twists have charge ground
energies separated only by O(1/L), so their energy difference cannot
serve as the cutoff. Physical momentum conservation is essential. At
half occupancy, a free charge excitation of physical momentum P obeys
the thermodynamic lower bound

    E_exc >= (4t/pi) dist(P,pi Z).

An internal twist change Q=(q+2pi ell)/b shifts the charge sea's
physical momentum by Q/2 modulo pi. A uniform interaction must restore
this momentum. For a nonforward internal q, the required charge
excitation therefore has threshold at least

    2t delta_int/(pi b)>0.

The proof and finite-ring qualification are in `CHARGE_THRESHOLD.md`.
This includes arbitrarily many particle-hole pairs, not just one pair.
It does not assert a gap for the whole charge sector.

[heuristic, Wilsonian matching] Choose a cutoff below that fixed
threshold, then choose J||h|| below the cutoff. Integrate only the
charge channels that require the finite momentum-restoration energy;
keep the low-momentum charge boson and every internal degree of freedom.
The resulting high-energy matching coefficients start at order J^2
and have finite cutoff-dependent denominators. In low-energy channels
the contact operator has a smooth expansion into the charge density
and higher derivatives. Its constant term supplies alpha H_h. A
density times a payload current is a Gaussian marginal coupling;
density times payload energy is dimension three. After normalization
by the charge and payload stiffnesses, the former tends to zero at
fixed h as J/t tends to zero. Half-filled charge umklapp has dimension
four and remains irrelevant. Nonforward residual channels are subject
to the finite mismatch analysis above. This is a standard Wilsonian
separation of fast and slow channels, not a projection out of the
gapless charge boson. A full local all-order construction is still
not supplied here; in particular individual finite-order coefficients
alone do not prove convergence or defect-sector stability.

With that matching input and the frozen-chain sector and band results,
the ordinary Gaussian RG conclusion is n+1=24 compact bosons. This
would violate c<=d-1=23 at d=24. The present file establishes the
resonance part of that implication and identifies the microscopic
input; it does not silently promote the implication to a measured or
rigorously established lattice counterexample.

## 7. Rings, rounding, and the crystal origin

[proved conditional on a positive dressed marker-defect gap] A legal
zero-defect squeezed crystal has N=bM. In a grand-canonical ground
state, rounding the preferred N to the nearest multiple of b changes
the smooth charge energy by O(t b^2/L), whereas a wrong-length defect
costs a fixed positive energy of order J times its internal gap.
Thus every sufficiently large L can use a compatible occupied number;
L itself need not be divisible by b. This removes the frozen chain's
obstruction on incompatible physical ring lengths.

[proved kinematic statement] If P is the payload translation momentum
per cell, cyclic translation by one occupied site has branches

    K=(P+2pi r)/b,       r=0,...,b-1.

These branches become charge boundary twists. They are precisely the
zero-mode gluing behind the fractional harmonic formula (1), not b
independent physical choices of an origin. They do not remove an
oscillator family. The marker-density order correlation has reference
exponent 2/b^2 and therefore decays at finite distance asymptotics;
fixed positional order is replaced by fluctuating density order.

[limitation] These statements do not establish a unique ground state
for every finite L. Exact charge-sector crossings or time-reversal
degeneracies may occur. Generic chemical potentials or a uniform
vacancy-hopping phase can avoid accidental crossings, but an all-L
uniqueness claim needs its own argument. Size-dependent zero-mode
rounding is compatible with the conventional compact-Gaussian CFT
interpretation used in this repository.

## Reproduction

Run `python3 contrib/counterexample_campaign/infrared/checks.py` from
the repository root. It uses only the Python standard library. It
reconstructs all 46 rational Fermi-root midpoint phases and all 1,081
unordered pair products, checks their cyclic separation with their
interval errors, and derives the vacancy margin (3). It also checks
the contact coefficient and harmonic thresholds. This rechecks the
root-to-resonance implication, not the original matrix inertia or the
interacting central charge.
