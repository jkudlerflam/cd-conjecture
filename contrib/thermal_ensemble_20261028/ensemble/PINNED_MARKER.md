# Physical pressure from all-length rooted cells

Status: **[proved finite-volume trace comparison]**, including arbitrary
length-dependent interactions. The comparison uses the sum over marker
numbers. It is not a fixed-marker-number comparison or an infrared theorem.

Provenance: OpenAI GPT-6-family Codex ensemble-transfer agent, 2026-10-28.
The proof was derived here. The parent agent independently suggested the
pinned-marker method. The matrix and counting checks are in `checks.py`.
No human review or microscopic central-charge measurement is claimed.

## 1. A local pin costs only a boundary free energy

Let `H_L=sum_i h_(i,i+1)` be a Hermitian nearest-neighbor periodic chain on
`(C^d)^tensor L`, with `L>=3` and `||h||<=J`. Let `p` be any rank-`q`
one-site projection, `1<=q<=d`, and let `P=p_1 tensor I`. Define

    Z_L(beta)=Tr exp(-beta H_L),
    H_pin=P H_L P restricted to ran P,
    Z_pin=Tr_(ran P) exp(-beta H_pin).

**[proved]** At every `beta>0`,

    (d/q) exp(-4 beta J) Z_pin
      <= Z_L
      <= (d/q) exp(4 beta J) Z_pin.                     (1)

No sign condition, spectral gap, translation symmetry, conserved charge,
or low-energy approximation is needed. The same proof permits unequal
bonds, replacing `2J` below by the sum of the two incident bond norms.

Proof. Delete the two bonds incident on site 1 and write

    H_cut=I_1 tensor H_rest,       ||H_L-H_cut||<=2J.

On `ran P`, the compression differs from `I_q tensor H_rest` by at most
`2J`. If two finite-dimensional Hermitian matrices differ in norm by at
most `c`, the min-max principle gives eigenvalues differing by at most
`c`, and hence their partition functions differ by factors between
`exp(-beta c)` and `exp(beta c)`. Apply this once to each space:

    exp(-2 beta J) d Z_rest <= Z_L <= exp(2 beta J) d Z_rest,
    exp(-2 beta J) q Z_rest <= Z_pin <= exp(2 beta J) q Z_rest.

Eliminating `Z_rest` proves (1). This is a trace inequality; no operator
monotonicity of the exponential is being asserted.

In particular, if the thermodynamic free energy exists, it equals the
pinned free energy. The error per physical site is at most

    [4J+T log(d/q)]/L.                                  (2)

The thermodynamic limit is taken at each fixed `T>0` before `T` tends to
zero. A very large but fixed `J` therefore does not spoil this conclusion.

## 2. The exact all-length coordinate space

Now suppose `p` is the marker projection and

    [h, p tensor I+I tensor p]=0.                       (3)

Choose orthonormal bases of the marker and nonmarker spaces. Every basis
word in `ran P` has a marker at physical site 1 and a unique decomposition
into ordered cells, each starting with a marker and continuing up to,
but not including, the next marker. Thus define

    K_ell = C^q tensor (C^(d-q))^tensor(ell-1), ell>=1,
    K = direct_sum_(ell>=1) K_ell,
    K_(M,L) = {vectors in K^tensor M with sum_j ell_j=L}.

There is an exact basis permutation, hence a unitary,

    ran P = direct_sum_(M=1)^L K_(M,L).                 (4)

This is the full physical alphabet. Illegal selector words, excited
internal states, consecutive markers, and arbitrarily long cells are
included. Marker internal states must be retained: replacing `C^q` by a
single state generally makes (4) false.

The dimension check is

    dim K_(M,L)=binom(L-1,M-1) q^M (d-q)^(L-M),
    sum_M dim K_(M,L)=q d^(L-1).                        (5)

Equation (3) implies that `H_pin` preserves every `M` summand. Denote its
transport to (4) by `H_pin,M,L`.

## 3. Adding the rooted seam

Suppose a rooted periodic cell Hamiltonian `H_root,M,L` on the same space
satisfies the uniform seam estimate

    ||H_root,M,L-H_pin,M,L||<=B,                        (6)

where `B` is independent of `M,L` and the individual cell lengths. It is
enough that the difference is a bounded interaction at the root seam.
Set

    R_L(beta)=sum_(M=1)^L Tr_(K_(M,L)) exp(-beta H_root,M,L).

**[proved]** Combining (1) with the spectral estimate in each summand gives

    (d/q) exp[-beta(4J+B)] R_L
      <= Z_L
      <= (d/q) exp[beta(4J+B)] R_L.                   (7)

Thus the full physical pressure equals the all-length rooted pressure,
with all marker sectors summed, whenever either thermodynamic limit
exists. No common payload bond, free-fermion diagonalization, momentum
congruence, or factorization of payload and length is required.

The constants do not depend on the largest retained cell length. The
theorem therefore applies directly to an infinite direct-sum cell space:
at fixed `M,L` only finitely many lengths and states occur, so every trace
used to prove (7) is finite. A trace-class hypothesis is needed later for
the length-fugacity ensemble, not for (7).

### Why (6) is a local condition rather than another infrared premise

For the literal finite-alphabet physical chain satisfying (3), one can
construct such a rooted Hamiltonian without any low-energy reduction.
Write `h` in its local matrix units. Matrix elements on two nonmarkers act
inside one cell. A matrix element involving one marker either preserves
its position or moves it by one physical site, transferring a boundary
letter between two consecutive cells. A matrix element on two markers
acts between consecutive length-one boundaries. Ordered markers never
cross; their internal labels may change.

Assign each such term to its moving marker, and for a pair of adjacent
markers assign the marker-preserving term to the marker on the right.
Each boundary rewrite is a coefficient times a partial isometry on two
cell factors. Its norm is at most the absolute coefficient, independent
of their lengths. Hermitian conjugates use the same boundary assignment.
There are at most two orientations for each of the `d^4` physical matrix
units. A deliberately loose bound for the boundary interaction is

    ||W_boundary|| <= 2 d^4 J.                         (8)

The sum of the nonmarker internal terms and all these boundary terms is
the rooted periodic operator. Under (4), all terms except root motion
agree with `H_pin`: moving the pinned marker is killed by compression.
The rooted operator retains those seam moves by transferring the end
letter between the last and first cells. Deleting the seam-moving matrix
units therefore gives `H_pin` exactly, and their partial-isometry bound
gives (6), for example with `B=2 d^4 J`. At `M=1` these seam maps are
partial permutations of the single cell's internal word; the same norm
bound holds. No denominator involving cell length occurs.

This construction proves an exact coordinate-transfer result. It does
not assert that this generic decomposition already has the positive
onsite operator and small intercell interaction needed for massive-mode
integration. Those stronger properties must be checked for the chosen
cell decomposition.

For a specified cell Hamiltonian whose intercell terms are bounded, a
much smaller value of `B` follows by deleting the root-crossing terms.
For a one-dimensional translation-invariant interaction `Phi_X`, an
absolute first spatial moment suffices:

    B_cut <= sum_(r>=2) (r-1) ||Phi_[1,r]|| < infinity.  (9)

There are at most `r-1` translates of an interval of size `r` crossing a
fixed cut. If the pinned and periodic conventions each differ from the
same open convention, use the sum of their two cut bounds. Onsite terms,
including unbounded positive cell energies, do not cross this cut and
are not included in `B_cut`.

## 4. What this says for C44/C45

**[proved application of the coordinate theorem]** C44's physical
alphabet is the four-state selector tensored with 1116 other states.
The marker projection is therefore

    p=|#><#| tensor I_1116,       d=4464, q=1116,
    d/q=4.                                               (10)

The specified diagonal terms, selector head moves, ports, readout and
marker hopping preserve the total number of `#` selectors. Thus (3)
holds for the complete repeated physical bond, including illegal words.
Its operator norm is finite for every fixed choice of its parameters.
Equations (4)-(8) give an exact rooted full-space pressure representation
for this bond, with all finite-wrapper length dependence retained.

This extends C46's common-bond quotient pressure result to an interacting
all-length coordinate representation. It avoids taking a trace with a
cyclic permutation insertion, whose sign would be uncontrolled for a
general quantum Hamiltonian. It also avoids multiplying a CFT central
charge by the number of translated marker origins.

**[application boundary]** The bounded-perturbation decomposition in C45
is stated inside a well-formed invariant representation. In malformed
multi-token sectors its readout can be unbounded as cell length grows.
The theorem above does not discard those sectors. To apply a massive
matching theorem proved only in that representation, one must separately
show that its rooted trace agrees with the full rooted trace to the
claimed low-temperature accuracy, or extend its onsite and interaction
bounds to all the sectors in (4). A ground-state exclusion by itself is
not the missing finite-temperature trace estimate.

Any unbounded readout term placed into the onsite `A_j` is harmless for
the seam comparison, provided the retained intercell terms are bounded.
This is not permission to move a genuinely intercell unbounded operator
into `A_j`. The actual decomposition must retain that distinction.

## 5. Limits of the result

- Equation (7) is for the full sum over marker numbers. It does not prove
  the same estimate in every fixed-`M` physical sector. The local cut
  factors over the free physical site only before fixing `M`.
- It preserves a thermodynamic free-energy function, including its
  low-temperature coefficients if those are established in the rooted
  ensemble. It does not prove that the coefficient is Gaussian or give
  a velocity-independent Virasoro central charge by itself.
- It proves no per-defect dressed energy inequality and no multichannel
  RG statement. Those are distinct remaining tasks.
- The length-fugacity conversion requires a trace-class neighborhood
  containing an interior zero of the cell pressure. A marker-free phase
  may instead put the physical pressure at a trace-class boundary. The
  exact finite-volume identity (7) remains valid in that case.

The new conclusion is an exact physical-to-rooted pressure transfer for
the full interacting chain. It removes cyclic gluing as an additional
thermodynamic assumption when the full rooted ensemble is controlled.
