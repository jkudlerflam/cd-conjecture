# A candidate refutation at the project's theoretical-physics standard

## Precise claim and qualification

The construction below gives a **candidate physics-level counterexample
family**, with a computably specified repeated two-site Hermitian bond.
The finite encoding, charge selection, initialized history core and
unperturbed mobile wrapper are exact construction claims. The conclusion
about the conformal limit uses standard local elimination of massive
degrees of freedom and stability of a nonresonant multicomponent
Luttinger liquid. These infrared inputs are stated below; they are not
claimed as new mathematical theorems.

Under those inputs there are finite-dimensional, one-site translation
invariant, nearest-neighbour chains with

\[
d_p\le 192C(p^2+p+1)
       [1+\log_2(2(p^2+p+1))]+48,
\qquad
c_p=(p+1)(p^2+p+1)+1,
\tag{1}
\]

for a fixed finite constant C and all sufficiently large primes p.
Consequently c_p/d_p tends to infinity, and in particular c_p>d_p-1
for sufficiently large p. This would refute c-d. It would **not** show
unbounded c at one fixed d: the local dimension in (1) grows with p.

This is an analytic candidate, not a numerically measured violation.
The constant C is given by a finite compiler construction but has not
been expanded into an optimized numerical alphabet count. The compiler
and the complete physical argument deserve independent end-to-end
review. Author: Codex, GPT-6 family, with three parallel research workers,
2026-10-28 UTC. No human review is claimed.

For the actual local operator assembly, see the
[microscopic bond specification](MICROSCOPIC_HAMILTONIAN.md). It makes
explicit which matrices are written down and which still require
expansion of the finite controller. No particular violating prime,
numerical local dimension, full bond matrix, or verified positive
coupling is currently instantiated.

## 1. More logical flavours than microscopic species

Take the point-line incidence graph of the projective plane over F_p.
It has u=v=p^2+p+1 vertices in its two parts, q=2u total vertices and
N=(p+1)u edges. Label one logical complex fermion by each edge. A fermion
on edge (alpha,beta) has the species charge

    w_(alpha,beta)=e_(upper,alpha)-e_(lower,beta).

The graph has girth six. If W is its oriented incidence matrix, any
nonzero integer vector z satisfying Wz=0 has l1 norm at least six.
It follows that every neutral monomial of at most four fermion fields
preserves each flavour separately. Neutral flavour conversions first
occur at six fields. This statement follows from a graph circulation
decomposition and is independent of the microscopic perturbation order.

See [the graph proof](../torus_bottleneck/README.md).

## 2. Exact ports and an initialized uniform core

The explicit two-port isometry W_p satisfies

\[
W_p^\dagger F_{L,a}W_p=\epsilon f_a^\dagger,
\qquad
W_p^\dagger F_{R,a}W_p=
\epsilon\left(1-\frac1{2N}\right) f_a^\dagger P,
\qquad \epsilon=\frac1{2\sqrt{uN}}.
\tag{2}
\]

It has a uniform polynomial circuit using charge-preserving local
queries, neutral Boolean operations and routing of existing species
symbols. Its amplitudes can be prepared exactly: a uniform interval
of M labels needs at most log_2 M exceptional rotations. An indexed
neutral rotation opcode contains these algebraic matrices without
requiring one controller state per label. The full data alphabet
remains O(q log q).

A reversible classical controller schedules this fixed-cardinality
gate menu without reading the quantum data. A single head overlays
only one data record at a time; it does not multiply two independent
q-dimensional charged records. The explicit fresh/updated head rule
is a unitary on the full two-site space and executes one machine step
per ordered sweep. The existing fifteen-letter clock supplies 2^n
sweeps with a complete legal path and positive illegal-sector penalties.

The input program, neutral scratch and one trusted boundary query helper
are pinned during initialization. A
coherent format checker selects exactly N occupation qubits and the
fixed charged reference words. An all-time FAIL penalty rejects every
other initialized data word; its noncommuting acceptance-gap estimate
also controls arbitrary invalid controller states. The controller
repeats preparation, flagged readout, inverse preparation and reset.
At every active endpoint time the same matrix (2) is read out.

The resulting finite core has positive local terms, exact zero-energy
code dimension 2^N, exact species conservation, and a positive gap at
one selected finite length. Its endpoint matrices are (2) times positive,
computable scalar time weights. Its local dimension is at most

    d_core <= C q(1+log_2 q),

where C is independent of p and cell length. The construction is
specified by finite algorithms for the local matrices, not by an
external execution of a circuit on the physical chain.

See [the sparse isometry](../exact_route/README.md),
[exact amplitudes](../exact_route/EXACT_AMPLITUDES.md), and
[the initialized-core proof](../torus_bottleneck/EXISTENTIAL_CORE.md).

## 3. Spontaneous mobile cells on every sufficiently large ring

Apply C8's positive-core wrapper, two-length selector and helix equalizer.
These apply to a complete initialized positive core of the preceding
form, independently of its internal implementation. They yield one
repeated bond on

    d_p = 96 d_core + 48

states. Its full unperturbed ground space consists exactly of cells of
two adjacent selected lengths, each carrying the core code. Every
sufficiently large ring has such a tiling. All other sectors have the
stated positive excitation cost; no marker or program sector is imposed
by hand. Marker hopping makes the cell lengths dynamical.

Give all added vacancy, selector and helix labels zero species charge.
Every wrapper move either preserves an occupied charged record or moves
it intact. This proves exact species conservation for the complete
microscopic bond. Marker number is conserved independently.

The equalized first-order projection is the direct length-fluid plus
payload Hamiltonian with the known momentum-gluing constraint. This
constraint affects zero modes and boundary twists; it does not remove
an oscillator family. The scalar marker motion intertwines the two
length codes. The already constructed wrapper therefore supplies the
mobile completion of this core, without a new copying operation on
its encoded quantum state.

See [selector and gap theorem](../../codex/work_C8/selector_wrapper.md),
[equalization](../../codex/work_C8/helix_equalizer.md),
[joint spectrum](../../codex/work_C8/joint_fluid.md), and the
[termwise symmetry proof](README.md#3-symmetry-inheritance-through-the-existing-mobile-wrapper).

## 4. A reference phase with strict stability margins

Choose the logical target

\[
H_{\rm pay}=\sum_{j,a}
[-t_a(c_{j,a}^\dagger c_{j+1,a}+\mathrm{h.c.})+n_{j,a}],
\quad
k_a=\frac\pi4+\frac{a\pi}{12(N+1)},
\quad t_a=\frac1{2\cos k_a}.
\tag{3}
\]

The positive hopping magnitudes lie in a bounded interval, and the
cell-coordinate velocities are tan k_a in (1,sqrt(3)). The physical
inter-port bond can assign each edge its own t_a. A common species
chemical potential supplies the density term, with its constant charge
offset per cell subtracted from the marker chemical potential.

Choose the length-fluid density nu in [5/12,7/12]. At first order its
velocity can be made comparable with all payload velocities by choosing
the marker amplitude. The projected model has N free complex fermion
channels and one free length-fluid boson, hence total Virasoro charge N+1.

For payload vertices write integer chiral charges u,v, z=u+v and m=u-v.
At the free point

    Delta_pay=(||m||^2+||z||^2)/4.

Species neutrality requires Wz=0. If z is nonzero, graph girth and
integrality give Delta_pay>=3. If z=0, only single backscatterings and
pairs of them can have dimension at most two. Their nonzero momentum
mismatches are bounded below by pi/[6(N+1)]. Exact total payload U(1),
contained in the species torus, forbids anomalous pairing.

C13's physical-coordinate locality dictionary gives the uniform
mixed-resonance condition and its dimension as

    payload momentum - 2pi h nu = 0 modulo 2pi,
    Delta_res = h^2 + Delta_pay,        h integer.

The only additional dangerous class is one length backscattering times
one same-flavour payload backscattering. Its mismatch is at least pi/6
throughout the chosen nu interval. Pure length harmonics at interior
commensurate densities have dimension at least four. Every other
resonant non-Gaussian vertex has free dimension at least three.

Thus the reference point has strict nonresonance margins and a strictly
positive full Gaussian matrix. A neighborhood in which the dimension
form stays above kappa times its reference value, kappa>2/3, has no
remaining relevant or marginal non-Gaussian vertex. Forward scattering
only deforms its Gaussian matrices. The neighborhood may shrink with
p and cell length; no uniform lower bound on its size is needed.

See [the infrared audit](IR_AUDIT.md) for the geometric charge matrix,
the operator enumeration, and the difference between a field-theory
statement and a thermodynamic lattice theorem.

## 5. Nonzero coupling, and the exact infrared input

Fix p, the full core, all selector and equalizer parameters, and the
desired first-order scale ratios. Normalize a finite microscopic
perturbation V_p so that P V_p P is the reference length/payload target
in common units. All small port weights are included in its finite
local norm B_p. Consider H_base+s V_p.

At any fixed connected perturbative order, virtual machinery excitations
produce local corrections starting at order s^2, whereas all retained
velocities are of order s. The unperturbed cell and wrong-length penalties
give positive local excitation costs. Taking one sufficiently small,
fixed s>0 gives the scale hierarchy described in
[SCALE_HIERARCHY.md](SCALE_HIERARCHY.md).

The physics-level conclusion uses these two inputs:

1. Local massive-sector elimination: for sufficiently small s, the
   excluded machinery remains massive and can be integrated out into
   a quasi-local low-frequency action whose correction, divided by s,
   tends to zero. Its exact microscopic symmetries are retained.
2. Gaussian infrared stability: a finite collection of compact bosons
   with positive stiffness, strict low-order nonresonance margins and
   only irrelevant non-Gaussian vertices remains a Luttinger liquid
   under sufficiently small such perturbations.

These are the standard effective-field-theory inputs at the project's
stated physics standard. There is no unknown marginally relevant sign
left to determine for this target. Numerical values of all second-order
coefficients are unnecessary for the existential small-s argument if
the first input holds: a strict neighborhood survives perturbations
of either sign. A rigorous uniform excursion expansion for this complete
mobile many-mode Hamiltonian has not been proved here.

In particular, a global spectral-band Schrieffer-Wolff argument is
inappropriate: the retained band has extensive width. No claim here
uses an O(L epsilon) eigenvalue error to control O(1/L) conformal gaps.
The two inputs concern local low-energy dynamics. All microscopic
parameters, including s, are fixed before L tends to infinity.

With those inputs the full chain has c_p=N+1. Combining this with the
exact alphabet count gives (1), and

    c_p/d_p grows at least as a positive constant times p/log p.

Different mode velocities are allowed in the conjecture's total c.
Density rounding and momentum gluing give the usual size-dependent
zero-mode sectors of a multicomponent Luttinger liquid. They do not
constitute a proof that one identical vacuum tower occurs on every L.

## 6. What should be challenged next

The highest-value next check is an independent review of the complete
finite-controller composition and the local massive-sector elimination.
An explicit alphabet expansion would make the construction easier to
audit and yield a concrete, though probably enormous, first violating d.
A useful numerical test would use a small compiled core and check
massive-defect stability as well as two independent central-charge
estimators; it could test the mechanism but not substitute for the
large-p dimension argument. The current scripts verify finite algebra,
clock completeness and local unitary rules, not the full thermodynamic
claim.
