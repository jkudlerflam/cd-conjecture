# Complete finite cells with initialized payloads and mobile boundaries

Status: proved open-chain theorem. It supplies a single repeated
nearest-neighbor bond with the complete desired ground code at two
adjacent lengths. It is conditional on an initialized frustration-free
payload core. The resulting bond by itself does not form these cells
on a periodic ring.

## 1. Statement

Assume a q-state uniform positive payload bond h_pay and positive
one-site boundary terms a_pay,b_pay have a frustration-free k-site
open core with complete code V:C^D -> (C^q)^tensor k, k>=2.

There is an explicit repeated nearest-neighbor Hermitian bond on
dimension

    d_wrap=6q+3                                           (1)

whose complete ground code, on each of the open lengths N=k+1,k+2,
is the initialized payload code tensored with a unique helix and a
vacancy charge state. No particle-number or tag sector is imposed
externally. Every other sector has a positive explicit energy
separation. Both endpoint maps are scalar multiples of the original
payload maps, and a local vacancy-label conversion gives scalar
transport when the boundary between two cells moves.

The proof combines the two-state length filter with a finite arrow
alphabet. Details below also specify the bond and the competing
sector bounds.

## 2. Alphabet and bond

Tensor the filter of length_filter.md onto each occupied payload
site. Its dimension is Q=2q, and define

    h0=h_pay tensor I + J I tensor h_phi,
    a=a_pay tensor I + J I tensor a_filter,
    b=b_pay tensor I + J I tensor b_filter,                (2)

with phi=2pi/(k+2) and J>=2U/epsilon_* as in that note.

Use three arrow vertices 0,1,2 and the following site letters.

| Letter | Arrow | Internal dimension |
|---|---|---|
| H_0 | 0 to 0 | 1 |
| H_1 | 1 to 1 | 1 |
| H_2 | 2 to 2 | 1 |
| L_alpha | 0 to 1 | Q |
| B_alpha | 1 to 1 | Q |
| R_alpha | 1 to 2 | Q |

The H letters are vacancies. A bond is legal when the outgoing
vertex of its left letter equals the incoming vertex of its right
letter. Let Pi_bad be the projector onto all mismatched pairs.

Let n project onto the three vacancy letters. The charge part has
diagonal U(n tensor I+I tensor n-2n tensor n) and the following
off-diagonal terms, with amplitude -w for each rewrite and its
reverse, for every alpha:

    H_0 L_alpha <-> L_alpha H_1,
    H_1 B_alpha <-> B_alpha H_1,
    H_1 R_alpha <-> R_alpha H_2.                           (3)

Every rewrite preserves the incoming and outgoing vertices of the
pair. Both sides are legal. Therefore each Pi_bad on every physical
bond commutes with every hopping term, including overlapping ones.
The other terms below preserve each letter's arrow, so the full
Hamiltonian commutes with all these violation projectors.

Let h0_tilde act as h0 on internal labels of the legal occupied
pair types LB, LR, BB, BR, preserving the types, and vanish otherwise.
Let a_L=|L><L| tensor a and b_R=|R><R| tensor b. Let P_L,P_R denote
the complete corresponding letter subspaces. The repeated bond is

    h_wrap = K Pi_bad + h_charge + h0_tilde
            +(a_L-mu P_L) tensor I
            +I tensor(b_R-mu P_R).                       (4)

Choose

    0<w<U/(k+1),   mu=2U,   K=2mu+2U=6U.                (5)

All coefficients are fixed for the selected k and independent of
the length of any later ring. The bond is Hermitian. It need not
be reflection invariant, owing to its oriented boundary terms.

The charge bond is positive when U>=w. Its only nontrivial blocks
are the disjoint two-by-two hopping blocks with diagonal U and
off-diagonal -w. All other vacancy-occupied pairs have diagonal U,
and all remaining charge entries vanish. Condition (5) implies
this positivity. The filter, payload, and boundary penalties in
(4) are also positive before the endpoint rewards are included.

## 3. Selection of every tag and occupation sector

First consider a legal open word. Its arrows form a path on the
directed graph above. Such a path contains at most one L and one R,
with L before R whenever both are present.

If both are present, the occupied sequence is L B^(r-2) R for some
r>=2. Vacancies before L must be H_0, between L and R must be H_1,
and after R must be H_2. For every set of particle positions the
vacancy labels are therefore unique. The rewrites (3) implement
all usual vacancy-particle exchanges on those positions. The
ordered-coordinate Hamiltonian is exactly the one in
tagged_vacancy.md, with a constant reward -2mu. Neither reward is
lost at an open endpoint, because L is not the last site and R is
not the first site when both exist.

At r=k the complete code is

    W_p=chi_p tensor V tensor |g_k>,    p=N-k in {1,2},    (6)

with energy E_p=e_p-2mu<=U-2mu. The length-filter theorem excludes
all r!=k sectors in this legal two-endpoint class by at least
U-w(N-1).

If a legal word contains just one endpoint tag, its reward is at
most mu. All other terms are positive, so its energy is at least
-mu. A word with no endpoint tag has nonnegative energy. Relative
to (6), these sectors lie at least mu-U=U and 2mu-U=3U higher,
respectively. Missing oriented rewards at physical endpoints only
increase these bounds.

Now consider an arbitrary sector with v>=1 bad bonds. Cutting at
these bonds produces v+1 legal open paths. Each contains at most
one L and one R, hence the total number of endpoint tags is at
most 2(v+1). Keeping the positive terms and bounding the rewards
gives

    H_wrap >= Kv-2mu(v+1).

The violations commute with the Hamiltonian, so this is a bound
on the entire corresponding invariant sector, not just on basis
configurations. Comparing with E_p<=U-2mu gives

    E_bad-E_p >= (K-2mu)v-U >= U.                       (7)

Thus all illegal sectors are excluded as well. The full ground
code on each length N=k+1,k+2 is exactly (6), with no imposed
sector or external boundary term. QED.

## 4. Gap and endpoint maps

Let gamma be the gap of the initialized k-site core in (2), delta_p
the charge gap, p_min the smallest adjacent-particle probability,
and

    M=(k-1)||h0||+||a||+||b||.

The complete open-chain gap obeys

    gap(H_wrap,N) >= min {
       delta_p p_min gamma /
         [(sqrt(delta_p)+sqrt(M))^2+p_min gamma],
       U-w(N-1) }.                                      (8)

The first term is the within-code-sector bound proved in
tagged_vacancy.md. The second is the smallest of the occupation,
missing-endpoint, and illegal-sector bounds above. This may be
very small at fixed k; it is not a thermodynamic gap theorem for
the eventual array of cells.

For a payload endpoint observable F, define its physical lift on
the L letters as F tensor I_filter, zero on other letters. The
first-site compression is

    W_p^dag F_L(first) W_p = r_p V^dag F(first) V,
    r_p=1-Prob_chi_p(first site vacant).                  (9)

The analogous R lift on the last physical site gives the same
factor r_p times the original last-site compression. Core
reflection symmetry is unnecessary.

## 5. Moving a boundary requires a vacancy conversion

Typed vacancies change the earlier transport rule. Let the initial
left and right cells have lengths k+1,k+2 and shift their shared
cut one site to the right, giving lengths k+2,k+1. The site crossing
the cut is initially the first site of the right cell. For a
nonzero projected overlap it must be vacant, since both cells
still have k occupied labels.

Before the shift that vacancy is H_0, since it precedes the right
cell's first occupied site. Afterwards it must be H_2, since it
follows the left cell's last occupied site. These are orthogonal
physical letters. The identity operator on the crossing site
therefore gives zero overlap, not the scalar overlap of the
untyped construction.

Use instead the local partial isometry

    C_vac=|H_2><H_0|                                    (10)

on the crossing site, with its adjoint for the reverse move.
All remaining vacancy labels and ordered payload labels agree.
Consequently

    (W_2 tensor W_1)^dag C_vac(crossing site)
       (W_1 tensor W_2) = beta I_(D^2),
    beta=[sum_{s=0}^{k} chi_1({s}) chi_2({s,k+1})]^2>0.   (11)

The filter states also overlap to one because their occupied
coordinate sequences are unchanged. The same reasoning proves
the identity for every payload V, not just for product payloads.
If a two-site marker hop moves the cut, (10) can be included on
the transferred site within that same two-site operator. This
does not yet prove the infrared behavior of the resulting model.

## 6. Why this bond alone does not solve the ring problem

The arrow graph was deliberately acyclic apart from its loops.
A legal closed path cannot contain L or R. On a periodic ring
with v>=1 mismatches there are v legal segments, so the stronger
cyclic bound #L+#R<=2v holds. Equation (4) then gives

    H_ring >= (K-2mu)v > 0

in every violating sector. All-H_0 and all-H_2 states are legal
zero-energy ring states. Thus the uniform periodic model (4) has
no ground-state cells of the desired type. A separate marker or
selector must provide the arrow resets at cell boundaries and
the intended boundary couplings. Applying (4) unchanged on a
ring is a failed construction, not an untested conjecture.

The theorem supplies a local finite-cell wrapper to insert into
such a selector. What it removes is the need to assume particle
number and endpoint-tag sectors inside each selected cell. The
missing high-channel initialized core, consistent mobile selector,
and finite-coupling CFT control remain separate open problems.

## 7. Verification

Run

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 arrow_wrapper_checks.py

The example uses k=3, a trivial one-dimensional payload, the qubit
filter, and d_wrap=9. It constructs the complete bond (4) and
diagonalizes the entire legal space at N=4,5, including absent
and incomplete endpoint tags. Their dimensions are 315 and 961,
inside full Hilbert spaces of dimensions 6561 and 59049. The
ground energies agree with e_p-2mu within 7e-13. The ground-code
residuals are below 2.6e-13. The measured gaps are .00196022849
and 5.06281047e-7, illustrating the small tunneling scale.

The program verifies the overlapping-bond violation commutators
exactly to numerical precision, positivity of the charge bond,
and #L+#R<=2(v+1) for every word of length five. This validates the
algebra used to bound the illegal sectors. Those sectors are
excluded by (7), not by a full diagonalization of dimension 59049.
The unconverted cut overlap is zero; after (10) it is
beta=.24951049536 with error below 3e-16. Results are written to
arrow_wrapper_checks.json. There are no CFT estimates here.
