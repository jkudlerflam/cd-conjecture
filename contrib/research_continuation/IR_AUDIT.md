# Infrared audit of the protected mobile first-order route

## Verdict

No unavoidable infrared instability was found for the sparse species
graph target. Its exact symmetry considerably simplifies C13's mobile
cell analysis: the additional mixed mismatch can be bounded below by
pi/6 independently of the number of payload modes. At each fixed finite
target size, a consistent weak-coupling scale hierarchy exists formally.

This is a **physics-standard implication to c=N+1 once the complete
microscopic core is supplied**, with the usual local massive-sector
integration and Gaussian RG inputs stated explicitly. At that standard,
no additional unknown sign of an effective coefficient has to be
calculated: after fixing the finite core, every correction is higher
order than the intended first-order Hamiltonian, and a strict Gaussian
neighborhood survives sufficiently small fixed coupling. Section 9
explains why the special cell structure supports this integration.

This is not a rigorous microscopic CFT theorem. A finite-order
Schrieffer-Wolff expansion or a bare global unperturbed gap alone would
not establish the thermodynamic conclusion. The unresolved mathematical
work is an all-order local excursion estimate and the standard
interacting infrared identification, rather than an identified
coefficient obstruction to the construction.

Author: OpenAI Codex, GPT-6 family, 2026-10-28. Sources inspected:
C13 RESULT; C8 selector_wrapper, helix_equalizer and joint_fluid;
ir_control README; the new torus_bottleneck, exact_route and
research_continuation notes. Statements below distinguish exact algebra,
standard bosonization input and missing microscopic estimates. No human
review or new thermodynamic simulation is claimed.

## 1. What the hypothetical core must supply

Assume a complete positive uniform core with exactly the intended
2^N-dimensional code, a positive finite-cell gap, exact species
U(1)^q, and the simultaneous endpoint maps in exact_route. Its gates,
initialization and discarded input sectors must all respect the charges.
Assume the C8 wrapper composition applies at a selected finite length b.

[proved, conditional on those inputs] The wrapper can carry the same
torus by giving vacancies, selector states and helix qubits zero species
charge. Vacancy exchange transports the occupied charge unchanged; the
marker conversion |H_2><H_0| is neutral. The two selected lengths contain
the same occupied core, so their representation offsets are identical.
Marker number M is conserved independently. Consequently a possible
background species charge C times M does not weaken the selection rule:
in a marker-conserving local process, species neutrality still requires
Wz=0 for the logical flavour-change vector z.

[proved] At zero perturbation the wrapper theorem gives a complete
ground manifold and an L-independent gap above that entire manifold.
It includes all allowed marker sectors, rather than an externally chosen
tiling. Helix equalization makes first-order payload hopping independent
of the two neighbouring lengths, and marker motion is a scalar on both
payloads. Its finite compensator norms must be retained in perturbation
bounds; they can grow polynomially in b or more if port visibility is
small.

[proved for the projected model] The joint physical-ring quotient ties
payload momentum to length-fluid momentum, but removes neither oscillator
family. The free projected model has N payload bosons and one length
boson, with central charge N+1 under the standard free-CFT identification.
The length mode is marker number, not an extra independent charge r:
at fixed L, r=L-bM.

## 2. The graph removes most of C13's dangerous operators

Let the species graph have girth six and oriented incidence matrix W.
Every nonzero integer z in ker W has ||z||_1>=6. Thus every neutral
quadratic or quartic fermion monomial preserves every flavour number.
This is exact symmetry selection, not emergent flavour symmetry assumed
at finite lattice coupling. Higher-order cycle conversions remain allowed.

In the Gaussian vertex conventions of torus_bottleneck, write chiral
integer charges u,v, and z=u+v, m=u-v. At the free point

    Delta_pay=(||m||^2+||z||^2)/4.

[proved] For Wz=0 and z!=0, integrality and m=z modulo 2 give
Delta_pay>=||z||_1/2>=3. For z=0, m=2u and Delta_pay=||u||^2.
The only nonzero payload vertices of dimension at most two are therefore
same-flavour backscattering and pairs of such backscatterings.

The species torus also contains total payload U(1): the sum of charges
on one bipartition counts logical particles up to its constant offset.
It forbids anomalous pairing exactly. Hence the parity-only construction's
need to displace all Fermi-sea centres does not apply here. Bands centred
at zero are permitted; differing widths are still needed to exclude
inter-flavour backscattering resonances.

## 3. A size-independent mixed-resonance margin

Keep the explicit payload target

    k_a=pi/4+a*pi/[12(N+1)], a=1,...,N,
    delta_pay=pi/[6(N+1)].                              (1)

All k_a lie in (pi/4,pi/3). The smallest dangerous payload mismatch
is 2|k_a-k_b|>=delta_pay. The single backscattering and same-sign
two-flavour mismatches stay of order one.

[derived using C13's locality/bosonization dictionary] For a physical
uniform resonance, label the marker harmonic by integer h. Its free
dimension and resonance condition are

    Delta_res=h^2+Delta_pay,
    payload_momentum-2pi*h*nu=0 modulo 2pi.              (2)

The physical marker Luttinger parameter is 1/(b+nu)^2. Equation (2)
already includes that factor and the squeezed-coordinate gluing; using
an unrelated free length fermion would give the wrong charge lattice.

The graph leaves just one new dimension-two mixed class:
length backscattering times a same-flavour payload backscattering.
Its mismatches are +/-2k_a +/-2pi nu. Choose

    nu in [5/12,7/12].                                 (3)

Then, uniformly in a and N,

    dist(+/-2k_a +/-2pi nu,2pi Z) >= pi/6.               (4)

Proof: 2k_a is in (pi/2,2pi/3), while 2pi nu is in
[5pi/6,7pi/6]. Their differences are separated from zero by at least
pi/6; their sums lie in (4pi/3,11pi/6), at least pi/6 from 2pi.
Other signs give the same distances. In particular nu=1/2 is safely
inside the allowed interval. A small irrational shift is optional.

Pure length h=+/-1 requires an endpoint density nu=0 or 1, excluded
by (3). Pure h>=2 terms have free dimension at least four, even when
commensurate. A quartic payload dressed by h=+/-1 has dimension at
least three. Nonzero graph-cycle conversions already have payload
dimension at least three, and length dressing cannot lower (2).

Thus C13's O(N^-2) mixed Sidon construction is unnecessary for this
target. Only the payload spacing (1) shrinks with N. Equations (1)-(4)
are checked numerically for the four recorded projective-plane graphs
in `ir_audit_checks.py`; their proofs do not rely on numerics.

## 4. The precise Gaussian neighborhood

[derived within bosonization] Require the full scaling-dimension form
to satisfy Delta>=kappa Delta_ref with kappa>2/3, and preserve the finite
mismatches (1),(4). Every remaining resonant non-Gaussian vertex then
has dimension strictly above two. Currents and their quadratic products
remain and renormalize the Gaussian theory.

The matrix condition must be checked in the correct joint coordinates.
C13's reference charge matrix is

    A_ref=T^T diag(v_0,...,v_N) T,
    B_ref=T^(-1) diag(v_0,...,v_N) T^(-T),
    T=[[b+nu,0],[-f,I]], f_a=k_a/pi.

The raw physical-coordinate Luttinger matrix is not near identity,
even at the free reference. In particular the marker parameter is
small at large b. An assertion that bare interaction coefficients are
small compared with the fastest mode does not imply stability.

For the actual quadratic theory require

    G=[[A,C],[C^T,B]]>0,

or equivalently A>0 and B-C^T A^(-1) C>0. The new
[explicit Gaussian certificate](../digit_compilation/audit/GAUSSIAN_CERTIFICATE.md)
replaces a continuity-only statement by the proved matrix bound

    delta = ||G_ref^(-1/2)(G-G_ref)G_ref^(-1/2)|| < 5/13,
    Delta_G >= sqrt((1-delta)/(1+delta)) Delta_ref > (2/3) Delta_ref.

This includes arbitrary density-current blocks C and noncommuting velocity
matrices. For example delta<=1/4 leaves every reference-dimension-three
vertex above 3 sqrt(3/5)>2.32. The relative threshold is independent of
N and b; converting a microscopic perturbation into this relative norm
still depends on the fixed core, geometry and smallest energy scales.
The reference geometry must use the actual dressed fillings, or the
fillings must be retuned; the stiffness bound alone does not certify the
momentum margins. This is a theorem about the Gaussian theory, not a
new microscopic elimination theorem.

A failure of positivity would be a real physical instability, such as
phase separation or a vanishing mode velocity. This is an independent
check, not repaired by momentum detuning. The global central-charge
claim also presupposes the usual interpretation as a compact Gaussian
collection when the velocities differ; a Casimir coefficient alone
measures sum c_a v_a, not the unweighted total c.

## 5. A consistent scale hierarchy exists at each fixed target size

Let eta>0 be the *complete* projected visibility of the payload link,
including history readout weights and any wrapper factors left after
equalization. Let chi>0 be the projected marker-hopping factor. Choose
the microscopic marker amplitude u proportional to lambda eta/chi,
so the two retained scales are comparable:

    W=lambda eta, J=u chi approximately W.

Let B absorb physical link norms, compensators, sums over channels,
the ratio eta/chi and chemical-potential counterterms. The microscopic
perturbation scale is then at most lambda B. Let Delta>0 be the actual
smallest local excluded-sector gap.

[finite-order estimate; conditional for a complete local elimination]
An induced local correction has schematic size

    epsilon_eff <= C(N,b) lambda^2 B^2/Delta + epsilon_core.

Consequently its size relative to W is controlled by

    epsilon_eff/W <= C(N,b) lambda B^2/(eta Delta)
                       + epsilon_core/(lambda eta).    (5)

To protect the spaced backscatterings, impose the stronger comparison
epsilon_eff/W much smaller than delta_pay, as well as the normalized
stiffness tolerance from Section 4. Equivalently, for an exact core a
sufficient schematic hierarchy is

    lambda/Delta much smaller than eta delta_pay/[C(N,b)B^2]. (6)

All factors are finite for a fixed constructed cell, so (6) has nonzero
solutions regardless of how small eta or Delta is. This is an existence
of compatible scales, not a proof of (5) to all infrared orders.
Approximate cores must additionally obey epsilon_core much smaller
than lambda eta delta_pay; simply sending lambda to zero with a fixed
core splitting makes the approximation worse relative to the target.

Independent physical hopping amplitudes allow the dressed k_a to be
retuned; the marker chemical potential controls nu. Retuning is optional
for existence, because sufficiently small continuous shifts already
preserve all strict margins. Their reference responses are nonsingular
in the displayed interior intervals. All choices must be fixed before
L tends to infinity.

## 6. The thermodynamic step is not supplied by global band SW

[proved limitation] The retained space has exponentially many states
and an extensive bandwidth O(LW) after coupling. Its excited logical
levels overlap energies of states containing a machinery excitation.
There is generally no globally isolated retained band at fixed lambda
as L grows. Requiring lambda L much smaller than Delta takes the wrong
order of limits. Nor can an O(L epsilon) truncation error resolve
conformal gaps of order 1/L.

This is not itself a reason the infrared construction fails: microscopic
spin chains routinely have massive degrees of freedom coupled to
gapless ones. The appropriate statement concerns low-frequency local
correlations or an effective action, not exact removal of an entire
extensive spectral band.

[needed local integration statement] The excluded selector, wrong-length,
workspace and clock excitations retain a positive dressed local cost.
Their virtual excursions produce a quasi-local low-frequency action,
with decaying temporal kernels and coefficients bounded in the relative
norm required by (5). Section 9 gives the cellwise energetic and locality
inputs for the standard connected-excursion argument. A bare global gap
would not supply those inputs for an arbitrary extensive perturbation.

At the theoretical-physics standard it is reasonable to test this by a
controlled local massive-sector integration and RG around the specified
Gaussian reference. It cannot be replaced by a finite-order SW formula
or a finite-size fit over a large preasymptotic crossover. At a strict
mathematical standard, the required general theorem for this interacting
mobile many-mode system is absent in the repository.

## 7. Locality, ring sectors, and the actual stopping point

Marker conservation keeps the number of squeezed sites fixed within
each sector. Finite-range even payload operators then map to finite-range
fermionic polynomials in squeezed order, with the usual parity strings.
The ring crossing gives the already derived global twist. No extra
local parity-twist perturbation is forced by the exact wrapper.
Extending this dictionary to the integrated quasi-local action is a
stated locality assumption of C13; it must remain part of the claim.

Generic density rounding gives size-dependent zero-mode offsets.
Avoiding all finite-L charge-sector crossings is a countable genericity
condition, not an open stability theorem. The conventional Luttinger
interpretation includes those offsets. If a stronger conjecture requires
one fixed untwisted vacuum tower along every integer L, C13 already
shows that generic incommensurate fillings do not meet it. The present
audit uses the conventional compact-Gaussian meaning.

The remaining implication is therefore explicit:

    correct positive core and ports
    + exact species and marker symmetries
    + a local massive-sector elimination satisfying (5)
    + dressed stiffness and momentum margins
    + standard compact-Gaussian infrared identification
    => c=N+1 at fixed nonzero microscopic couplings.

No fatal graph-induced mixed resonance was found. At the physics
standard used by C13, local infrared transfer follows from the usual
massive-sector integration and stable Gaussian RG argument once the
complete core meets these hypotheses. Individual higher-order
coefficient values are not necessary for the existential conclusion:
their finite bounds and extra power of coupling suffice. At a strict
mathematical standard, neither a finite-order projection nor the
following physics argument is a proof of the full CFT limit.

## 8. Review of the noncommuting acceptance-pin gap lemma

The independent acceptance lemma communicated by the torus worker is
valid. Let H_init>=gamma(I-P), P be its exact initialized-history kernel,
P_good<=P, and A>=0 satisfy ||A||<=B, A P_good=0 and
PAP>=a(P-P_good), with a>0. Then

    gap(H_init+A)>=a gamma/(a+B+gamma).                 (7)

For a vector orthogonal to P_good write x=P psi, y=(I-P)psi.
Its energy is at least

    gamma||y||^2+(sqrt(a)||x||-sqrt(B)||y||)_+^2.

Where the bracket is positive, the corresponding two-by-two quadratic
form has determinant a gamma and trace a+B+gamma, giving (7). Where
it vanishes, ||y||^2>=a/(a+B), which gives the stronger bound
a gamma/(a+B). This also proves that P_good is the complete kernel.

Thus all pulled-back pins need not commute on malformed controller
inputs. For a classical reversible verifier and all-time local FAIL
penalties, valid inputs must never raise FAIL and every invalid
initialized input must raise it at least once; then a>=eta/C_n on the
bad initialized-history space. B must bound the *sum* of the physical
FAIL terms, for example eta times their number. A terminal-only pin can
miss an error that was subsequently uncomputed. The lemma does not
construct or verify the underlying format-checking algorithm.

## 9. Why the local massive-sector integration is controlled in physics

### 9.1 A local defect-energy bound, not just an isolated global level

[proved from the C8 decomposition] Fix a legal classical marker pattern.
At zero inter-cell perturbation the cells are independent. Let Q_j be
the complement of the complete selected code in a cell of length b or
b+1, including its selector head and helix states. Let Delta_* be the
minimum of the two finite-cell gaps, and gamma=t eta_b-U>0 the
wrong-length coefficient of selector_wrapper. On this marker sector,

    H_base-eL >= Delta_* sum_{selected cells j} Q_j
        + gamma sum_{wrong-length cells j}(1+|l_j-b|).   (8)

The selected-cell projectors act on disjoint cells. Their sum commutes
with the unperturbed Hamiltonian. The second sum is a scalar on each
length pattern, so it commutes as well. This counts multiple machinery
excitations with an energy proportional to their number. Arbitrarily
long wrong cells pay an energy proportional to their volume, rather
than an O(1) domain-wall cost.

The selector-violation projectors are conserved by the specified marker
and port perturbations. They therefore do not occur as virtual
excitations of a legal state. C8 gives their separate sector cost at
least U per violation. Long marker-free fragments also have a bulk
penalty relative to eL<0: their unshifted lower bound is independent of
fragment length, whereas subtracting eL adds |e| times that length.
Small local perturbations can change fragment boundary energies by a
finite constant per boundary and bulk energies by O(s B_p) per site.
This supports retaining the legal sector when s is sufficiently small;
it does not use the norm of an extensive perturbation as a denominator.

Set delta_* = min(Delta_*,gamma) for legal virtual excursions. A local
microscopic perturbation insertion touches at most two adjacent cells.
For selected cells their physical size is bounded by b+1. For a wrong
cell of length l, the number of places at which an insertion can touch
it is O(l), while (8) pays at least a fixed multiple of l at large l.
Consequently the number of active local insertion positions is bounded
by a fixed constant times the defect weight in (8). Constants here may
depend arbitrarily badly on the one fixed cell, but not on total L.

### 9.2 Irreducible excursions retain all critical propagation

[controlled physics expansion] Normalize the target as
H(s)=H_base+s V_p, with local norm B_p and a first-order retained
Hamiltonian s H_ref. Keep propagation within the entire retained
manifold unexpanded. In particular do not try to invert its small
logical level spacings. An irreducible virtual excursion starts when a
local insertion creates a machinery defect and ends when the last such
defect disappears; intervals spent entirely in the retained manifold
are compositions of effective vertices, not part of one irreducible
vertex.

For a prescribed excursion with r insertions, unperturbed high-sector
propagation over a time interval tau with defect weight n has norm at
most exp(-delta_* n tau), by (8). The ordered time integrals thus supply
inverse defect energies. The sum of local insertions touching n active
defects is O(n B_p) (including the weighted long-cell case), so the
factor n in the branching count is compensated by the n in that energy.
Purely retained insertions away from the active defects cancel against,
or are absorbed into, the retained propagator. This is the linked
excursion reorganization; expanding those insertions independently
would restore the spurious extensive bandwidth problem.

The standard connected-cluster estimate resulting from this
reorganization has the form

    ||K_r||_local <= C_p B_p^2 s^2/delta_*
                      (C_p |s| B_p/delta_*)^(r-2),     (9)

where the local norm sums connected spatial supports containing a
specified cell. Its combinatorial constant incorporates the finite
cell dimension, finite range, and tree/ordering factors. In a
time-resolved description the same resummation gives the schematic
majorant

    ||K(tau)||_local <= C_p s^2 B_p^2
                         exp[-(delta_*-C_p|s|B_p)tau]. (10)

Each new spatial step requires an insertion. Thus the irreducible
vertices also have exponentially decaying spatial tails for
C_p|s|B_p/delta_*<1. In physical coordinates their range unit is the
fixed cell size, not one microscopic site. Integrating (10) gives
O(s^2 B_p^2/delta_*) corrections, and its frequency expansion has
relative derivative suppression by omega/delta_*.

Equations (9)-(10) are the usual massive-mode linked-cluster estimate,
not a proved convergence theorem in this note. Equation (8), bounded
support per insertion, and the weighted long-cell cost supply its
concrete microscopic inputs. A fully mathematical proof must construct
the connected expansion with moving marker projectors and establish
its combinatorial bound; simply asserting a gap of H_base would omit
that work. At the physics standard the expansion is controlled by
|s|B_p/delta_*, and gives a dressed machinery cost
delta_*-O(|s|B_p). Critical propagators are not hidden inside its fast
kernels: they remain in the retained theory.

### 9.3 Consequence for the existential claim

Fix p, the exact finite core, its chosen length and all compensator
norms first. Every constant in (8)-(10) is then finite. Choose a
positive tolerance epsilon_p lying inside the Gaussian positivity,
dimension and momentum margins of Sections 3-4. Choosing one fixed
nonzero s such that

    C_p |s| B_p^2/delta_* << epsilon_p

makes the normalized retained Hamiltonian H_ref plus an
epsilon_p-small quasi-local deformation. Exact species and marker
symmetries survive the integration. Every allowed resonant
non-Gaussian vertex has free dimension at least three; the finitely
many lower-dimension oscillatory vertices keep nonzero mismatch, and
forward scattering only deforms the positive Gaussian form.

[physics-standard conclusion] With the usual compact-boson RG
identification, this proves the existence of an N+1-mode phase for
sufficiently small fixed nonzero s, conditional on the complete exact
core/compiler and its wrapper interfaces. No unknown sign of a
second-order coefficient can overturn this sufficiently small-coupling
conclusion. It is unnecessary to compute each such coefficient before
making an existential claim at this standard. The construction still
needs its algebraic/compiler hypotheses verified; and the statement
is not upgraded here to a rigorous microscopic CFT theorem.
