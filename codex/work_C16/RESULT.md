# C16. Compactness, singular points, and the missing exhaustion theorem

Completed 2026-10-27. No counterexample to c <= 1 was found. The requested
finite, checkable reduction does not follow from compactness and the
c-theorem. There are two distinct missing inputs: an exhaustive description
of the singular neighborhoods, and an exclusion of finite-c > 1 fixed points
elsewhere in the family. Operator counting supplies neither input.

There are useful partial results. This note gives an exhaustive classification
of the dimensions of frustration-free ground spaces, exact tests for a
polarized vacuum and its magnon sector, a cubic band-touching example on the
free locus, controlled single-fluid neighborhoods, and conditional operator
counts for the proposed product CFTs. A finite list of research targets is
possible. A finite list whose tests certify the whole conjecture is not yet
available.

Tags apply to the derivations or tables immediately following them.
[proved] means derived here or explicitly attributed to a supplied repository
proof. [assumed] means standard continuum input with its hypotheses stated.
[heuristic] means a proposed scenario or an unproved extrapolation.
[numerical] means the computation, sizes, and estimator are specified. None
of this work has had an independent referee check. No external literature or
network was consulted.

## 1. What compactness actually gives

[proved] Write the physical Hamiltonian as

    H(J,B) = sum_i [sigma_i^T J sigma_(i+1) + B.sigma_i].

After removing a bond coboundary, an additive constant, and positive scale,
fix ||J||_F^2 + |B|^2 = 1. The sphere S^11 modulo the compact SO(3) action is
compact, with generic dimension eight. Singular symmetry strata do not change
this conclusion. The zero Hamiltonian is excluded. Negative scale is not an
equivalence of ground-state problems.

Let C be the set of parameters with a conformal infrared limit, and define
c only on C. Let R_reg consist of points p admitting an open neighborhood
U_p and a finite number M_p such that

    c(q) <= M_p for every q in U_p intersect C.                 (1)

If c(q_n) tends to infinity, compactness gives a subsequence q_n -> p. That
p cannot belong to R_reg, since eventually q_n would lie in U_p. Conversely,
if every p belongs to R_reg, a finite subcover gives a finite global bound
max_j M_(p_j). Neither statement gives the number one.

An isolated conformal point p with c(p)=3/2 is consistent with (1), for
example with M_p=2, and with all nearby critical points having c <= 1. It
need not lie near a nonconformal point. The continuum RG picture of a finite-c
multicritical point flowing down to lower-c theories is consistent with the
c-theorem. It is not evidence for its realizability by this qubit bond, but
it disproves the proposed logical inference about isolated points.

[assumed] The usual argument for (1) near a finite-c fixed point requires
that nearby lattice theories enter a common continuum RG neighborhood, that
no additional sector drops below the matching scale, and that a monotone
c-function applies. One chooses a finite crossover scale at which that
function is bounded and then uses monotonicity toward the infrared.
Continuity at each fixed finite lattice size alone does not justify this
uniform-in-volume matching. Competing vacua and diverging crossover lengths
must be controlled. If the continuum has several velocities, the ordinary
Lorentz-invariant c-theorem does not automatically apply to their coupled
flow. Exact mode counting or a separate continuum argument is needed.

A genuine perturbation of a UV CFT with c_UV <= 1, within the hypotheses of
the c-theorem, has c_IR <= 1. A lattice Hamiltonian whose *infrared* is that
CFT is not automatically a common UV completion for every nearby lattice
Hamiltonian. This distinction matters when trying to cover the entire
parameter space.

[proved] A spectral gap above an arbitrarily degenerate ground space is not
the same as a stable gapped phase. Section 4 gives a normalized, nonzero
Hamiltonian with a gap of one above exponentially many zero modes and
arbitrarily small critical perturbations. Thus even the gapped part of a
compactness argument needs a stability hypothesis.

## 2. Product vacua, ferromagnets, and saturation

### 2.1 Exact product-vacuum and one-magnon tests

[proved] Rotate a candidate uniform product state to |0>^L, with Z|0>=|0>.
For rings of length at least four, the only off-diagonal actions on this
state create one spin flip or a pair of adjacent flips. Their amplitudes are

    A_2 = J_xx - J_yy + i(J_xy + J_yx),
    A_1 = B_x + J_xz + J_zx + i(B_y + J_yz + J_zy).             (2)

Therefore the product is an eigenvector if and only if A_1=A_2=0. Its energy
is L(J_zz+B_z). These equations test an eigenstate, not a ground state.
Ground-state status still needs a lower bound or a spectrum calculation.

Set J_xx=J_yy=a, J_xy=-J_yx=D, J_zz=b. The block of H on states with one
flipped spin has on-site energy -4b-2B_z relative to the vacuum and hopping
amplitudes 2(a+iD) to the right and its conjugate to the left. With
|k>=sum_j exp(ikj)|j>, its dispersion is

    epsilon(k) = -4b - 2B_z + 4a cos(k) + 4D sin(k).           (3)

This is an exact excitation dispersion only if the one-magnon sector is
invariant. Equations (2) alone do not guarantee that.

For example, write Z=1-2n and impose the single-flip cancellation in (2).
The remaining mixed xz terms include

    -2 sum_i [J_xz X_i n_(i+1) + J_zx n_i X_(i+1)],          (4)

and likewise with Y. These terms can turn a single magnon into two adjacent
magnons. Invariance of the entire one-magnon subspace, for a generic ring of
length at least four, requires J_xz=J_zx=J_yz=J_zy=0. Together with (2),
this is the axial U(1) locus. For example J_xz=g, J_zx=-g, B_x=0 preserves
the vacuum but violates magnon-number invariance.

[numerical] checks.py verifies (3) at L=7 with maximum spectral error
5.33e-15. Adding g=0.13 in (4) leaves the vacuum residual zero but gives a
one-magnon leakage norm 0.97283. This is an explicit failure of using bare
spin-wave counting at an arbitrary factorized point.

### 2.2 The controlled single-species saturation transition

[proved] On the U(1) locus, put r=sqrt(a^2+D^2). For r>0, (3) has one
minimum modulo 2pi, at k_0, and

    epsilon(k_0+q) = delta + 2r q^2 + O(q^4),
    delta = -4b - 2B_z - 4r,       m = 1/(4r).                (5)

The curvature cannot vanish unless r=0, when the whole one-magnon band is
flat. A nonflat quartic one-magnon minimum and two inequivalent one-magnon
minima are therefore excluded in this particular sector. A two-particle
bound state or a competing finite-density phase can preempt delta=0. The
one-magnon calculation is not a proof that it is the actual saturation line.

[assumed, controlled dilute gas] Suppose delta=0 is the first instability,
there is only this one low-energy species, its effective mass is positive,
its short-range interaction is repulsive with finite scattering length, and
all other excitations remain separated. Set mu=-delta. The low-energy action
is

    S = integral d tau dx [psi* (d_tau - d_x^2/(2m) - mu) psi
                           + g (psi*psi)^2/2 + ...].          (6)

At mu=0 the vacuum transition has z=2, not a conformal central charge. For
mu>0, arbitrarily low density takes a repulsive one-dimensional gas to the
impenetrable limit. Equivalently the ordered particles have the leading
thermodynamics of spinless fermions:

    k_F = pi rho,
    e(rho) = pi^2 rho^3/(6m) - mu rho + O(rho^4),
    rho = sqrt(2m mu)/pi + O(mu),
    v = pi rho/m + O(rho^2),       K -> 1.                    (7)

There are two Fermi edges and one density/phase canonical pair. Bosonization
gives one compact boson, c=1. The vacuum side is gapped. Thus the largest
nearby c is one under the hypotheses of (6). Attractive collapse, threshold
bound states, an infinite scattering length, and an independent soft species
are explicit exceptions to the derivation, not exceptions already ruled out
for the full family.

A source for a nonconserved magnon number can remove the vacuum transition.
Pair creation can leave an Ising transition. Neither permits treating the
U(1) dilute-gas formula as a theorem on a generic product-state neighborhood.

### 2.3 Ferromagnetic and spiral points

[proved] For H=-J sum sigma_i.sigma_(i+1), J>0, every ground state belongs
to the fully symmetric spin-L/2 multiplet, of dimension L+1. Equation (3)
gives

    epsilon(k)=4J(1-cos k)=2J k^2+O(k^4).                    (8)

This is z=2. A pure polarized ground state is a product; other vectors in
the multiplet can have substantial entanglement. In a Dicke state with fixed
extensive magnetization, the Schmidt weights are a hypergeometric
distribution with variance proportional to L. Its entropy is
(1/2) log L+O(1), which a naive periodic CFT fit would label c=3/2. The
quadratic excitation spectrum rejects that interpretation.

[assumed] Smooth fluctuations of a ferromagnetic unit vector have the
Landau-Lifshitz action

    S = integral dt dx [s A(n).d_t n - rho_s (d_x n)^2/2
                       - V(n) + ...].                       (9)

Expanding around an ordered direction gives a Berry term
s(pi_x d_t pi_y-pi_y d_t pi_x)/2. The two transverse coordinates are a
single canonical pair. At V=0 their equations give omega=(rho_s/s)k^2.
If easy-plane anisotropy gives one transverse coordinate a positive local
quadratic cost, integrating it out gives a term proportional to
(d_t phi)^2 together with (d_x phi)^2, hence one boson, c=1. Easy-axis
anisotropy or a field can gap the fluctuations. This proves the count inside
the smooth one-field continuum description. Its validity throughout every
direction of the full microscopic neighborhood is an additional assumption.
In particular, the ferromagnetic point is not a repulsive gas with a
uniformly finite scattering length.

[proved] A site-dependent rotation about z with angle jQ maps an axial
XXZ+DM bond to an XXZ bond with transverse coupling magnitude r. At the
rotated isotropic ferromagnetic value, polarized directions become product
helices and (8) is shifted in momentum. On a ring there is a twist LQ, so
nonpolarized helices require compatible twist to be exact finite-ring ground
states. Bulk mode counting is unchanged. This covers the spiral points
unitarily equivalent to the ferromagnet, not all factorized spiral states.

[assumed, supplied input] The full axial U(1) family and zero-field XYZ
family obey c<=1 using the exact-solution input recorded in C6. This settles
those slices, including their singular limits. It does not settle arbitrary
transverse or mixed-bond perturbations out of them.

## 3. Band edges, Lifshitz points, and domain walls

### 3.1 A cubic edge already occurs on the quadratic locus

[proved] Consider a quadratic JW Hamiltonian with BdG eigenvalues

    E_+(k), E_-(k) = u sin k +/- sqrt[(m+t cos k)^2
                                      + Delta^2 sin^2 k].    (10)

The odd term is a tilt from uniform DM. Set m=-t and u=Delta>0. For small
positive k,

    E_+(k) = 2 Delta k + O(k^3),
    E_-(k) = -t^2 k^3/(8 Delta) + O(k^5).                   (11)

To obtain the second equation without cancellation, rationalize it:

    |E_-(k)| = t^2(1-cos k)^2 /
        [sqrt(t^2(1-cos k)^2+Delta^2 sin^2 k)+Delta sin k].

The numerator is t^2 k^4/4+O(k^6), and the denominator is
2 Delta k+O(k^3). There is a cubic slow branch as well as a linear branch.
A concrete Pauli bond has J_xx=1/2, J_yy=J_zz=0,
J_xy=-J_yx=1/4, and B_z=1/2. Momentum reversal may reverse the tilt
under a different JW convention, without changing the conclusion.

[numerical] With t=Delta=1, halving k from 0.1 to 0.0125 gives exponents
2.999324, 2.999831, 2.999958. |E_-|/k^3 tends to 0.125.

[proved, supplied free theorem] All conformal points in the quadratic
neighborhood still have c<=1. The example shows that a list containing
only quadratic saturation and ferromagnetic singularities is incomplete.
It does not establish a bound after arbitrary interactions are added.
At a cubic edge, finite-size velocities vanish, so fitting a common nonzero
velocity to the spectrum is invalid.

### 3.2 Pokrovsky-Talapov under its actual hypotheses

[assumed] Suppose a commensurate phase has a gapped amplitude mode and one
species of the favored domain wall, with gap Delta_dw, mass m_dw, short-range
repulsion, and no internal zero modes. A tilt h lowers its energy by Q_dw h.
Near h_c=Delta_dw/Q_dw,

    epsilon_dw(k) = Delta_dw-Q_dw h+k^2/(2m_dw)+... .          (12)

The resulting dilute soliton gas is (6). Its density is proportional to
sqrt(h-h_c), its incommensurate wavevector shift is proportional to that
density, the transition has z=2, and the doped side has c=1. The wall of
opposite topological charge remains massive in this ordinary transition.
A commensurate-incommensurate boundary meeting a BKT endpoint can instead
have a relativistic c=1 description. Neither has two independent density
fields under these hypotheses.

Alternating kink labels alone do not give independent species: if a wall
and antiwall must alternate around the ring, their densities are constrained.
Conversely, if a second independent wall species becomes soft, or a wall
carries a critical internal sector, (12) is insufficient. Checking the wall's
internal spectrum and the conversion vertices is essential. The finite FF
witness theorem does not bound the internal spectrum of frustrated walls.

### 3.3 What a genuine two-minimum Lifshitz neighborhood would require

[heuristic, conditional effective theory] A dressed low-energy excitation
can have

    epsilon(k) = delta + alpha k^2 + beta k^4 + ...,
    beta>0.                                                  (13)

At alpha=delta=0 it has z=4. For alpha<0 it has minima at
+/-sqrt(-alpha/(2 beta)). Equation (5) excludes this mechanism for an exact
single U(1) magnon band, but not for a dressed composite or a general
interacting excitation. Integrating out modes generates longer-range
processes even when the microscopic bond is nearest neighbor.

Expand the field about two minima Q_1,Q_2. If both species have finite
density and remain stable, the candidate action has two kinetic terms and
an interaction matrix g_ab rho_a rho_b/2. A two-fluid phase requires a
positive compressibility matrix, two nonzero stiffness eigenvalues, and
absence or irrelevance of every resonant conversion and number-changing
vertex. It then has c=2. A negative compressibility eigenvalue gives
phase separation or an instability; a relevant independent locking term
generically removes a mode. Bare weak-coupling miscibility inequalities are
only a first test, since one-dimensional dilution is strongly interacting.

For symmetric minima +/-Q, pair conversion psi_+^dagger^2 psi_-^2 has
momentum 4Q and usually oscillates. But without a particle-number symmetry,
psi_+ psi_- has zero momentum. For ordinary local bosonic fields its
bosonized dimension, at diagonal stiffness, is

    Delta_pair = 1/(4K_+) + 1/(4K_-).                         (14)

At dilute repulsive K_+=K_-=1 it is relevant. Thus the disappearance of
inter-valley scattering by incommensurability does not by itself protect
a c=2 phase. For local fermionic fields the permitted charge lattice and
statistics must be used instead. Two generic shifted centers need not sum
to zero, but all other resonances must then be checked. No rank-J argument
here excludes their microscopic emergence.

## 4. Frustration-free ground spaces and their neighborhoods

### 4.1 An exhaustive ground-space dimension classification

[proved] Fix a nonzero positive semidefinite two-qubit bond h. Define
S=ker h. The open-chain ground space is the intersection of the constraints
that each neighboring pair lies in S. The following is a classification of
zero-space dimensions, not of excitation gaps or perturbative phases.
It also applies when a Hamiltonian first has to be shifted by a scalar and
a bond coboundary to obtain its positive FF representative.

| rank h | Forbidden/allowed space | Open-chain zero space at long length | Periodic zero space |
| --- | --- | --- | --- |
| 1 | Forbidden vector entangled | Dimension L+1 | At most L+1 |
| 1 | Forbidden vector a tensor b, a,b independent | Dimension L+1 | Dimension 2 |
| 1 | Forbidden vector a tensor a | Dimension F_(L+2) | Dimension F_(L-1)+F_(L+1) |
| 2 | Allowed space dimension 2 | At most 2, or no extension | At most 2, or none |
| 3 | Allowed space dimension 1 | One repeated product if it extends to three sites | One, or none |
| 4 | No allowed vector | None | None |

Here F_0=0, F_1=1, and the periodic statements are for L>=3. Rank zero
would be the excluded zero Hamiltonian. The rank-two proof and the rank-three
extension criterion use the supplied C8 allowed-space classification.
The details needed for the table are as follows.

1. For a rank-one entangled forbidden vector with coefficient matrix M,
   the bilinear constraint on a product a tensor b is a^T M* b=0. Let
   epsilon=[[0,-1],[1,0]] and T=epsilon M^dagger. Then

       M* T = det(M*) epsilon.

   T is invertible. Apply the invertible site-dependent map
   G=1 tensor T tensor T^2 tensor ... tensor T^(L-1).
   On bond i,i+1 the constraint transforms to a nonzero multiple of
   epsilon, because A^T epsilon A=det(A) epsilon. Thus the transformed
   open zero space has no singlet on any adjacent pair. Each adjacent swap
   acts as the identity, and adjacent swaps generate all permutations.
   The space is Sym^L(C^2), of dimension L+1. Closing the ring adds a
   constraint and can only reduce this dimension. If T is proportional to
   a unitary, this is the twisted-ferromagnet construction after unitary
   normalization. General T is not unitary, so no spectral equivalence is
   claimed.

2. If the forbidden vector is a tensor b with independent a,b, choose a
   single invertible local map A taking the two bras a^dagger A and
   b^dagger A to the coordinate bras <0| and <1|. The forbidden substring
   becomes 01. Open binary strings without 01 are 11...1100...00,
   of which there are L+1. On a ring only the two constant strings survive.
   Invertibility preserves the kernel dimensions. Again, it need not
   preserve gaps or conditioning in the large-L limit.

3. If a and b are proportional, a uniform unitary makes the forbidden
   substring 11. The open count satisfies D_L=D_(L-1)+D_(L-2), with
   D_1=2,D_2=3, so D_L=F_(L+2). Closing the ring gives the trace of the
   two-state adjacency matrix to the Lth power, namely
   F_(L-1)+F_(L+1). This grows as phi^L, where phi=(1+sqrt(5))/2.

4. For rank two, C8 classifies the determinant on the projective line P(S).
   If every vector on the line is a product, S=a tensor C^2 or C^2 tensor b,
   and only one endpoint can remain free on an open chain. If there are two
   distinct product directions, allowed extensions follow a directed graph
   on two vertices with at most one incoming and one outgoing edge per
   vertex. There are at most two long paths, so at most two ground vectors.
   If the line is tangent to the product quadric, it has the form

       span(a tensor b, a tensor d+c tensor b).

   When b is not parallel to a there is no extension to three sites.
   Otherwise a uniform invertible basis and subtraction of the first
   vector put it in the form span(|00>, alpha|01>+beta|10>), with
   alpha beta nonzero. There is the vacuum and one weighted W state.
   Any amplitude with two excitations can be moved, using the local
   01/10 relation, until the excitations collide, where the forbidden 11
   condition sets it to zero. Hence the dimension is at most two.

5. For rank three, S=span(v). Compatibility of v_12 tensor x with
   y tensor v_23 forces v to have Schmidt rank one; compatibility of the
   shared site then forces v=a tensor a. There is a single repeated product.

Consequently an exponential FF ground-space degeneracy can occur only at a
rank-one product projector with identical factors, up to unitary rotation
and positive scale. Polynomial degeneracy remains possible in the other
rank-one cases. This is more restrictive than a generic "degenerate FF"
label, but it is not a theorem about c in their neighborhoods.

[numerical] ff_checks.py verifies representative cases on L=3,...,7,
with both open and periodic chains. A generic entangled forbidden vector
has open dimensions 4,...,8 and periodic dimension 2. The singlet parent
has dimension L+1 for both boundaries. The aligned product has open
5,8,13,21,34 and periodic 4,7,11,18,29. A rank-two tangent example has
dimension 2 throughout. These are checks of the algebraic classification,
not fitted central charges.

[assumed here, supplied by the concurrent C15 result] On final review,
C15 had completed an explicit nonunitary-filter example in the rank-one
entangled row. Its pairing-only fermion Hamiltonian has nearest-neighbor
pair amplitude delta, repulsion V>2|delta|, zero hopping, and
mu_f=V-sqrt(V^2-4delta^2). After an additive shift its bond is
V|v><v| with |v>=A|00>+B|11>, A^2=mu_f/(2V), B^2=1-A^2.
On even rings the ground dimension is L+1. A closed defect sector has

    E_-(k)=V-sqrt[V^2-4delta^2 sin^2(k/2)]
           = delta^2 k^2/(2V)+O(k^4),

where k is momentum per two-site cell. This is a concrete nonconformal
FF family beyond a unitary image of the isotropic ferromagnet. C15 also
shows how effective two-site hops arise from virtual pairing near the
blockade corner and why they do not constitute two independent stable
fluids on that exact line. Those computations were not rerun here; see
../work_C15/RESULT.md, equations (10)-(14). Its perturbed neighborhood
remains unresolved. This fits, rather than extends, the kernel table.

### 4.2 The exponentially degenerate parent has a concrete effective family

[proved] Choose the representative

    H_0 = U sum_i n_i n_(i+1),   n_i=(1-Z_i)/2,   U>0.         (15)

Let P project onto strings with no adjacent occupied sites. H_0 has gap U
above P, but dim P grows exponentially. For any perturbing uniform bond
V specified by J,B, define, within P,

    N = sum_i n_i,
    T_hop = sum_i (b_(i+1)^dagger b_i + b_i^dagger b_(i+1)),
    C_hop = i sum_i (b_(i+1)^dagger b_i - b_i^dagger b_(i+1)),
    X_c = P sum_i X_i P,      Y_c = P sum_i Y_i P.

Here each hopping term is also projected by P and b^dagger creates a down
spin. Directly on the allowed binary strings,

    P sum XX P = P sum YY P = T_hop,
    P sum XY P = -P sum YX P = C_hop,
    P sum ZZ P = L-4N,
    P sum XZ P = P sum ZX P = X_c,
    P sum YZ P = P sum ZY P = Y_c,
    P sum Z P = L-2N.                                      (16)

For instance, the pair-creation pieces of XX and YY produce adjacent
occupancy and are projected out. Whenever a constrained spin flip survives,
its adjacent spin is empty, so multiplying by that neighbor's Z leaves it
unchanged. Therefore

    P V P = L(J_zz+B_z) + (-4J_zz-2B_z)N
          + (J_xx+J_yy)T_hop + (J_xy-J_yx)C_hop
          + (B_x+J_xz+J_zx)X_c + (B_y+J_yz+J_zy)Y_c.        (17)

This is the entire first-order tangent family at (15): a complex hopping,
a chemical potential, and a complex source, with exclusion. A global z
rotation can remove the phase of a nonzero source; it cannot remove the
hopping phase without making the source spatially dependent.

[numerical] At L=8, dim P=47, a random full J,B gives a maximum entrywise
error 1.78e-15 in (17).

[assumed, perturbation theory] With a perturbation lambda V small compared
with U, the next term is

    -lambda^2 P V Q (Q H_0 Q)^(-1) Q V P,    Q=1-P.           (18)

It generates further-range and multibody operators inside the constrained
space, even though the original perturbation is two-site. An extensive
thermodynamic expansion should be interpreted locally; a bound on the
full extensive norm ||lambda V|| is not uniform in L. At a fine-tuned
critical point of (17), a small relevant term in (18) can control the
ultimate infrared. One cannot discard it when claiming an entire
neighborhood has c<=1.

The witness K_FF(2)=4 says that a positive bond admitting a zero vector on
four open sites admits one at every length. It does not bound dim P, remove
(17), or constrain the spectrum generated by (18).

### 4.3 One solvable direction, with two independent estimators

[proved] Set the source and current to zero and take H_eff=-t T_hop at
fixed particle number N<L/2. On an open ordered configuration
x_1<...<x_N with x_(j+1)>=x_j+2, the squeezed coordinates y_j=x_j-j
remove the excluded sites. The hopping becomes free hard-core hopping.
On a ring the squeezed length is M=L-N, with a boundary twist tied to total
momentum and particle parity. The zero-momentum ground state has occupied
momenta spaced by 2pi/M, symmetric about zero. Summing their cosines gives

    E_0(L,N) = -2t sin(pi N/(L-N)) / sin(pi/(L-N)).           (19)

At fixed rho=N/L, define k_F=pi rho/(1-rho). Expanding (19),

    e_infty(rho) = -2t(1-rho) sin(k_F)/pi,
    E_0 = L e_infty - pi t sin(k_F)/[3(1-rho)L] + O(L^-3).

A particle-hole excitation at the Fermi edge has energy
2t sin(k_F) times 2pi/(L-N), so

    v = 2t sin(k_F)/(1-rho).

The Casimir term is -pi v/(6L), hence c=1. This determination of v from
excitations is separate from the ground-state Casimir coefficient.
Differentiating the displayed bulk energy twice gives

    e_infty'' = 2pi t sin(k_F)/(1-rho)^3 = pi v/K,
    K = (1-rho)^2.                                        (20)

Thus one component can have K below 1/2 under nearest-neighbor exclusion;
it approaches 1/4 on the dilute-hole side of close packing. This is a
constrained strong-coupling limit of the controlled U(1) class, not a claim
that a second fluid can coexist with it. The hopping vanishes at the fully
packed endpoint, and the conformal formula is used only for 0<rho<1/2.

[numerical] checks.py diagonalizes the constrained ring at t=1, rho=1/3,
where e_infty=-4/(3pi), v=3, and K=4/9. The ground-state eigenvector is
embedded in the original spin tensor product when computing entropies.
No squeezed-coordinate entropy is substituted for physical-spin entropy.

| L | N | Hilbert-space dimension | c from Casimir, using exact e_infty and v | c from within-ring entropy profile |
| --- | --- | --- | --- | --- |
| 12 | 4 | 105 | 1.018289 | 1.176465 |
| 15 | 5 | 378 | 1.011636 | 1.132326 |
| 18 | 6 | 1386 | 1.008054 | 1.103162 |
| 21 | 7 | 5148 | 1.005906 | 1.087472 |
| 24 | 8 | 19305 | 1.004516 | 1.073984 |

The largest energy error against (19) is 2.5e-14, and the largest eigenvector
residual is 1.6e-13. A fit of physical half-chain entropy to
S=(c/3) log[(L/pi)sin(pi floor(L/2)/L)]+constant across these five sizes
gives c=1.010480. Adding a 1/L^2 correction gives c=0.997497. The
within-ring fits use all cuts 3<=ell<=floor(L/2); their larger values show
the effect of short-distance corrections and are not evidence of c>1.

These checks establish a solvable c=1 direction out of (15). They do not
bound the entire source-and-current family (17), much less all higher-order
corrections. The exponentially degenerate parent is a concrete unresolved
neighborhood, even though its controlled number-conserving slice is known.

## 5. What symmetry and operator counting do and do not decide

### 5.1 Locality, translation, and the stiffness matrix must be specified

[assumed, bosonization] For n bosons with no density-current mixing, use

    H_LL = (1/2pi) integral dx [d_x theta^T V_theta d_x theta
                               + d_x phi^T V_phi d_x phi],
    [phi_a(x),d_y theta_b(y)] = i pi delta_ab delta(x-y).

Both V matrices are positive. Define the positive K matrix by
K V_phi K=V_theta. A local vertex with charge vectors l,m has dimension

    Delta(l,m) = [l^T K l + m^T K^(-1) m]/4.                (21)

In a diagonal basis this gives dim(exp(2i phi_a))=K_a and
dim(exp(2i theta_a))=1/K_a. Generic chiral theories also have
density-current terms and require the full Gaussian covariance matrix
in place of this restricted formula. Positivity and the number of
independent canonical pairs still decide the number of Gaussian modes.

A potential perturbation has the form

    g_O exp(i q_O x) O(x) + Hermitian conjugate.             (22)

It survives coarse graining only if q_O is zero modulo a reciprocal lattice
vector of the appropriate translation-invariant sector, belongs to the
local physical operator algebra, and is invariant under all exact symmetries.
Its dimension then decides relevance. These data are not fixed by naming
the on-site group Z2, an antiunitary, or the trivial group. In a broken
translation phase the enlarged cell and domain sectors also have to be
specified.

An antiunitary is not a second unitary parity. Its action conjugates
coefficients and can exchange vertices. For a real Ising order field sigma,
A sigma A^-1=+sigma or -sigma decides whether its real coupling is allowed.
For complex bosonic vertices one must transform the full Hermitian
combination in (22). The labels "time reversal" and "plane reflection"
alone do not supply that emergent representation.

### 5.2 Ising times Luttinger

[proved, conditional on the candidate operator algebra] An Ising factor has
local primaries 1, sigma, epsilon with dimensions 0,1/8,1. A product with
one boson has total c=3/2. The Ising thermal mass epsilon is relevant and
symmetry even. If sigma is local and uniform, it is a second relevant
perturbation unless protected by an exact spin flip, antiunitary parity,
or translation momentum pi. A disorder field is not assumed independently
local; if the microscopic embedding permits additional local fields, they
must be counted too.

| Residual stratum | Specified action on sigma | Necessary relevant tunings from the Ising primary sector alone |
| --- | --- | --- |
| One unitary Z2 | sigma odd | At least the thermal mass, one |
| One unitary Z2 | sigma even and momentum zero | Thermal and magnetic fields, two |
| Antiunitary only | sigma odd | One |
| Antiunitary only | sigma even and momentum zero | Two |
| No internal symmetry | sigma has momentum pi | One |
| No internal symmetry | sigma has momentum zero | Two |

These are counts of independent continuum coefficients that must vanish,
not established codimensions inside the microscopic qubit family. Extra
mixed perturbations can increase the count. For a bosonic local vertex V
of dimension d_V, epsilon V has dimension 1+d_V and sigma V has
1/8+d_V. A product can be invariant even when its factors separately are
not, for example sigma and V both odd under the only Z2, or both carrying
momentum pi. One must test all such products.

Even with no nontrivial resonant boson vertices, currents have dimension one.
A permitted sigma d_x phi or sigma d_x theta has dimension 9/8 and is
relevant. Tuning the coefficient of sigma to zero does not in general tune
these derivative couplings to zero. The permitted epsilon d_x phi and
epsilon d_x theta have dimension two and require marginal RG analysis.
Velocity mismatch, which is allowed here, is not an additional massive
sector and does not by itself remove c=3/2.

[heuristic] With a protecting action forbidding sigma and all dangerous
mixed operators, an Ising mass crossing inside a stable one-boson phase is
a codimension-one candidate. Without protection and with no additional
selection rule, at least two tunings are needed and often more. Neither
scenario is excluded by a 4-, 5-, or 8-parameter lattice family. Neither is
proved reachable merely because enough parameters exist.

### 5.3 Two Ising factors times Luttinger

[proved, conditional] The proposed product has c=2. Before involving the
boson, its relevant tensor-product primary fields are

    epsilon_1, epsilon_2                      dimension 1,
    sigma_1, sigma_2                          dimension 1/8,
    sigma_1 sigma_2                           dimension 1/4,
    sigma_1 epsilon_2, epsilon_1 sigma_2       dimension 9/8.  (23)

The product epsilon_1 epsilon_2 has dimension two and is marginal.
Assume all listed fields are in the local product operator algebra. With
both sigma fields at zero momentum, the following necessary counts follow
by multiplying their symmetry signs.

| Exact action on sigma_1,sigma_2 | Number of invariant relevant primaries in (23) |
| --- | --- |
| One common Z2 flips both | 3: the two masses and sigma_1 sigma_2 |
| One Z2 flips only one | 4: two masses, the even sigma, and its product with the other epsilon |
| Neither is flipped | 7 |
| Antiunitary sign action on real primaries | Same counts for the corresponding signs |

With no internal symmetry and both sigmas uniform, the primary count is
seven, within the eight-dimensional generic family before other constraints
are imposed. If both have translation momentum pi, the count falls to three.
If only one does, it is four. A unitary Z2 together with translation can
act as two independent parities, leaving just the two masses from (23).
Thus "one on-site Z2" does not imply "only one selection rule."

Derivatives and mixed boson vertices can add further relevant couplings;
(23) is a lower count, not an exhaustive stability analysis. In particular
there is no universal codimension for this CFT in any of the three symmetry
strata. A symmetry can also act by exchanging factors, which requires a
separate invariant-combination calculation.

### 5.4 Two Luttinger liquids and phases without internal symmetry

[assumed] A two-boson phase has c=2 when both stiffness eigenvalues are
positive and no allowed relevant perturbation pins an independent
combination. No Ising mass is present. It can therefore have codimension
zero if every local resonant vertex is irrelevant or marginally harmless.
Incommensurate momenta can remove number-changing terms and conversions
without an exact on-site U(1), but this must be checked term by term.

For example, atom-pair conversion exp[i(2theta_1-theta_2)] at diagonal K
has dimension

    Delta_convert = 1/K_1 + 1/(4K_2),                       (24)

and momentum 2Q_1-Q_2. If it is local, resonant, and relevant, it locks the
relative phase and leaves at most one of these two Gaussian modes. If it
oscillates or is irrelevant, that one obstruction is absent, not all of
them. A pair of symmetric valleys has the separate danger (14).

Once a positive K matrix and a discrete local charge lattice are specified,
only finitely many charge vectors can have Delta<2. If lambda_* is the
smallest eigenvalue of diag(K,K^-1), (21) implies

    Delta >= lambda_* (|l|^2+|m|^2)/4.

Thus |l|^2+|m|^2<8/lambda_* suffices for a finite relevant-vertex search
at that candidate. There is no uniform search cutoff over parameter space
unless lambda_* has a positive lower bound. Soft or infinitely stiff
limits are among the singular neighborhoods at issue.

A microscopic realization must supply (i) two independent pairs of low-energy
edges or two independent hydrodynamic fields, (ii) positive susceptibility
and stiffness matrices of rank two, (iii) a local operator charge lattice
and momenta excluding every relevant locking process, and (iv) a common
scaling window in which both velocities are nonzero. Distinct peaks in a
static structure factor are insufficient. One boson produces many harmonics.

[proved limitation of the proposed rank argument] rank J<=3 is an automatic
algebraic fact about a qubit bond. It is not a count of continuum species.
The same on-site observable can have several low-energy harmonics, and
composite operators extend across sites. At second order, even a product
of two adjacent microscopic bonds produces terms of the form

    (sigma_i^a sigma_(i+1)^b)(sigma_(i+1)^c sigma_(i+2)^d)
      = delta_bc sigma_i^a sigma_(i+2)^d
        + i epsilon_bce sigma_i^a sigma_(i+1)^e sigma_(i+2)^d. (25)

Actual effective coefficients include denominators and Hermitian sums, but
the increased support is already explicit. Equation (18) gives a concrete
qubit example. Hence bare nearest-neighbor Fourier degree does not bound
the Fourier degree of a dressed band. A theorem relating rank J to a
current-complete Gaussian theory would need that extra hypothesis, which
the generic residual set does not have. This does not prove c=2 is possible
at d=2; it identifies why the rank argument fails to exclude it.

### 5.5 Why these counts cannot become a CFT classification

[proved, logical] Let g_alpha(p) be the relevant continuum coefficients
induced by lattice parameters near a proposed fixed point. Reachability
requires simultaneous equations g_alpha(p)=0, together with appropriate
marginal flow. For a regular constant-rank zero set the local codimension
is the rank of the constraints. Symmetry only removes some coefficients;
it neither proves there is a zero nor computes the remaining Jacobian.
A family may have identically vanishing or correlated couplings, intersect
a stable manifold nongenerically, or miss it altogether. A count exceeding
the number of lattice parameters is a genericity obstruction, not an
exclusion theorem.

Moreover, product Ising and Gaussian CFTs do not exhaust unitary CFTs with
c>1. A small on-site group need not equal the emergent continuum symmetry.
Compactness does not give a finite catalog of possible operator algebras.
Any proposed list containing only the examples in Sections 5.2-5.4 leaves
an additional, unresolved non-product-CFT case. C17 is responsible for the
broader mechanism and anomaly investigation; no such exclusion is assumed
here.

## 6. The reduction that is justified, and the finite list that is not

### 6.1 A precise conditional covering statement

[proved] Let G be the union of open neighborhoods on which c<=1 has actually
been established, by an exact mapping or a controlled common continuum
limit with the required stability hypotheses. Let S=P\G, where P is the
compact parameter space. Since G is open, S is compact. To prove the full
conjecture it suffices to construct, for every p in S, an open neighborhood
with the same bound. A finite subcover would then exist.

This is a covering theorem, not a procedure for finding the neighborhoods.
In particular:

- A controlled *slice* such as the U(1) or free locus does not itself give
  an open neighborhood in the eight-dimensional space.
- A bound c<=M_p with M_p>1 cannot serve as one of the c<=1 neighborhoods.
- Compactness does not bound the crossover length, number of disconnected
  critical sets, number of possible emergent CFTs, or required numerical size.
- An isolated finite-c > 1 CFT is an independent obstruction even if all
  nonconformal singular neighborhoods have been controlled.

[assumed, useful conditional local result] A neighborhood whose only soft
variables are one repulsive dilute species or one canonical phase-density
pair, with all conversions and competing sectors controlled, has conformal
limits c<=1. A relativistic scaling limit obtained by relevant or marginal
perturbation of a unitary c<=1 CFT also has c<=1 under the c-theorem's
hypotheses. Sections 2 and 3 exhibit where these statements apply and where
their premises can fail.

No argument here exhausts the full family by these local descriptions.
Consequently the requested statement "c<=1 unless one of a finite set of
specified, finitely testable structures exists" remains unproved. Adding
"some other c>1 CFT" to a list would make it logically exhaustive while
removing its value as a reduction. That is the unresolved issue, not an
exception that can be dismissed by parameter counting.

### 6.2 A finite research list, with concrete tests

[heuristic, research organization] The following six items are a useful
triage list. They are not a finite classification theorem. In particular,
items 2 and 6 still contain unbounded families of effective descriptions.
Each proposed numerical test can diagnose a candidate; a finite negative
test does not exclude an entire neighborhood without scale control.

| Item | Structure that could escape the controlled single-mode arguments | Concrete calculation | Status here |
| --- | --- | --- | --- |
| 1 | A nonstandard band edge or velocity-collapse point | Compute the lowest dispersions and local exponents; test whether an exact few-particle sector is invariant | Exact U(1) one-magnon band has one quadratic minimum; mixed-term leakage demonstrated; a free cubic edge exhibited |
| 2 | Two independent dilute species, a soft composite, or a wall with internal critical states | Locate all one-, two-, and low-cluster thresholds; compute scattering/conversion amplitudes, momentum resonances, and the susceptibility matrix | Single repulsive species controlled; no general composite or wall exhaustion theorem |
| 3 | A degenerate FF parent with a nontrivial perturbative theory | Classify ker h; compute PVP and the first relevant term in PVQ(H_0)^(-1)QVP; inspect independent modes | Kernel dimensions classified; all first-order perturbations of the exponentially degenerate parent derived; one c=1 direction verified |
| 4 | Finite-c product multicriticality, such as Ising times LL or two Isings times LL | Resolve translation and symmetry of primaries; construct the relevant-operator response matrix; tune its independent coefficients and test mixed RG couplings | Necessary primary counts and missing mixed couplings given; no reachability claim |
| 5 | A stable incommensurate multi-boson phase | Extract a rank-two stiffness/susceptibility matrix, independent gapless towers and momenta, and enumerate all relevant local resonant vertices | Required data specified; rank J<=3 does not close this case |
| 6 | An emergent non-product CFT, an unclassified singular limit, or failure of uniform RG matching | Identify conformal towers and fusion/selection data; track onset scales and all velocities toward the accumulation point | Not excluded; prevents calling the list a finite checkable reduction |

Items 1 and 2 must distinguish an actual excitation from a projected trial
band. A direct numerical check is ||Q_m H P_m||, where P_m projects onto
the proposed fixed-magnon sector. Nonzero leakage invalidates an exact
few-particle dispersion calculation unless a controlled elimination replaces
it. The example in Section 2.1 gives a cheap calibration.

For item 3, exact kernel dimensions on a few open chains can establish an
FF extension when the bond is known positive, by the supplied K_FF(2)=4
proof. They cannot establish its low-energy neighborhood. At the blockade
parent, (17) makes the next test concrete: include its complex source and
hopping phase, and retain symmetry-allowed relevant terms from (18). This
is a much smaller problem than scanning all of P, but it is still an
interacting phase-diagram problem rather than one finite matrix test.

For item 4, use the complete continuum selection rules before assigning a
codimension. Lattice perturbations integrated over the ring couple low-energy
states through matrix elements with leading scaling L^(1-Delta). Their
symmetry-resolved matrix elements can estimate the rank of the induced
relevant-coupling map. This requires identifying scaling operators and
checking size dependence, not merely counting available Pauli coefficients.

For item 5, when an effective charge is conserved, twists and particle-number
responses give stiffness and compressibility. For an emergent charge without
an exact lattice conservation law, first identify its low-energy density and
current operators; arbitrary twists of physical spin axes are not equivalent
to independent twists of every emergent fluid. A positive rank-two response
at a finite size is useful only if it persists as the scaling window grows.

[numerical protocol, proposed rather than run] For an actual c>1 candidate,
use both physical-spin entanglement scaling and the momentum-resolved
finite-size spectrum. With several independent nonchiral sectors,

    E_0(L)=e_infty L - (pi/6L) sum_a c_a v_a + o(1/L),
    S(ell)=(sum_a c_a/3) log[(L/pi)sin(pi ell/L)]+s_0+... .   (26)

The velocities in the first equation must come from separate excitation
towers, not from the same Casimir fit. With chiral drift use instead
-pi sum_a c_a(v_(a,R)+v_(a,L))/(12L). Boundary twists can shift zero-mode
energies; account for them before using (26). Vary sizes compatible with
incommensurate momenta and include a third diagnostic, such as low-temperature
entropy or the counting of independent low-energy towers, when two velocities
are too close to resolve. Small-ring checks can reject a false candidate;
they cannot rule out all criticality below a much longer crossover scale.

A first-order coexistence of two phases is a direct sum of sectors, not a
tensor product of gapless fluids. Their central charges are not added.
A cat state contributes at most the logarithm of the number of macroscopically
distinct components to the entropy, which is an additive constant for a fixed
number of phases. Extensive degeneracy and rotor towers require separate
scaling tests, as the ferromagnetic example illustrates.

## 7. Reproduction, scope, and next steps

[numerical] From this folder run

    python3 checks.py
    python3 ff_checks.py

Only NumPy, SciPy, and Python's standard library are used. BLAS thread counts
are set to one. Both checks run locally in seconds. The scripts write
checks_results.json and ff_checks_results.json. The largest constrained
basis diagonalized has dimension 19305 (L=24,N=8); the ground-kernel checks
use L<=7. All data were generated here, not read from previous scans.

The following parts of the request were completed: exact controlled
single-species descriptions and their failure conditions; a cubic free
singularity; the FF ground-space dimension classification; the full blockade
projection; conditional relevant-operator counts across the requested
symmetry strata; explicit criteria for multiple incommensurate modes; and
two-estimator checks of one singular neighborhood direction. The requested
exhaustive classification of all z!=1 points, actual multicritical
reachability and codimension, and finite certification of the full c<=1
claim remain open. No positive c>1 example was found or asserted.

Ranked next steps are:

1. Test the source-and-current constrained family (17), including the leading
   relevant terms in (18), near its nontrivial critical intersections. It is
   an explicit perturbative neighborhood of the only exponentially
   degenerate positive FF qubit parent. Coordinate its overlap with C15's
   strong-coupling analysis before running another scan.
2. At a concrete Ising-line/floating-phase intersection supplied by C15 or
   existing numerics, identify translation quantum numbers and the mixed
   Ising-boson operators. Determine whether it is one boson with a vanishing
   velocity or an independent Ising sector. The cubic example is a useful
   countercheck against a false c=3/2 interpretation.
3. Establish a microscopic single-fluid theorem beyond invariant U(1)
   sectors, or exhibit its failure through a dressed two-minimum band with
   all allowed locking vertices controlled. Bare rank J and bare Fourier
   degree do not provide this theorem.
4. Seek uniform RG-onset or singular-exhaustion input before using
   compactness to claim a finite numerical certificate. More random
   finite-size points cannot replace that missing step.

Local sources used: ../../STATUS.md; ../../specs/_d2_background.md;
../work_C6_d2/c6.tex for the parametrization and controlled slices;
../work_C8/qubit_ff_witness.md for the exact FF extension classification.
C15 and C17 had no progress files at the initial read. Both were read at
final review. C15's factorization line is attributed explicitly in Section
4.1; no new C15 or C17 computation was rerun or used to infer a global bound. The project-wide STATE.md was read
for orientation but was not modified, in accordance with the instruction
to write only inside this folder. PROGRESS.md is the resumable memory.
