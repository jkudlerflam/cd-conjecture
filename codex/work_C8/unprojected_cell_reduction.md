# An exact cell-coordinate reduction before projection

Status: proved representation of the legal marker sectors. This reduces
the local elimination problem to an atomic model in cell coordinates.
It is not an all-order elimination theorem or a CFT identification.

## 1. The exact local Hilbert space

Fix the microscopic parameters, including the selected finite b, and
let s=6Q+3 be the wrapper alphabet. In a legal marker sector with M>0
markers, every cell has a selector word

    # L^(a-1) X R^(ell-1-a),    1<=a<=ell-1, ell>=2.

No low-energy projection is used here. The cell retains every wrapper,
filter, clock and work state. Including the two equalizer qubits,
its full Hilbert space is

    K_ell = C^(ell-1) tensor (C^(4s))^tensor ell,
    K = direct_sum_(ell>=2) K_ell.

The first factor stores the selector head position. Define the on-cell
length operator ellhat by ellhat|K_ell=ell I. For a rooted ordered list
of M cells, the physical length constraint is sum_j ellhat_j=L.
This is a global conserved charge constraint, not a local truncation
of K to its two low lengths.

[proved] There is a basis bijection between legal physical words with
a chosen reference marker and an absolute reference position x in
Z/LZ together with the ordered cell tensor word satisfying that length
constraint. Changing the reference marker identifies

    (x;u_1,...,u_M) ~ (x+ell_1;u_2,...,u_M,u_1).

Quotient by this free reference-label action, with normalized orbit
sums, recovers the physical legal Hilbert space. The action is free
because the reference marker changes, including when cells have
identical internal states. Physical translation shifts x. This is the
same ring-coordinate issue as joint_fluid.md, now with full cell spaces.

## 2. The unperturbed Hamiltonian is on-cell

[proved] On the rooted space, after subtracting e_base L, the full
unperturbed Hamiltonian is

    H0-e_base L = sum_j A_j,

where on K_ell

    A_ell = mu_selector I + H_head(ell)+H_wrap(ell)
              +H_equalizer(ell)-e_base ell I.

H_head is the open (ell-1)-vertex selector path. H_wrap is the complete
open wrapper, not its selected subspace. H_equalizer contains both
uniform helix bonds and the two endpoint pins. Pins on a physical
cut separate into a last-site B pin for the left cell and a first-site
A pin for the right cell, so no inter-cell base interaction remains.
Each marker's chemical potential is counted once.

The existing selector and wrapper proofs imply A_ell>=0. At ell=b,b+1
its kernel is the complete D-dimensional selected code, D=2^m, and
the smallest positive eigenvalue has the stated positive cell lower
bound. At every other ell,

    A_ell >= gamma_wrong(1+|ell-b|) I,
    gamma_wrong=t eta_b-1>0

for U=1, using the conservative inherited selector bound. The flat
wrapper actually permits a stronger convexity estimate, which is
not needed here. Thus A on K has a finite ground space

    ker A = K_b^0 direct_sum K_(b+1)^0
           isomorphic to C^2 tensor C^D,

and an on-cell gap Delta_at>0 independent of M and L. It also
penalizes sufficiently long cells proportionally to their length.

The ground-space projectors in this representation are fixed on-site
projectors. Their motion in physical coordinates is a coordinate
effect, not a failure of an on-cell decomposition of the full H0.

## 3. The perturbation is bounded and local in cell order

[proved] The complete microscopic perturbation becomes

    V = kappa M + sum_j F_j + sum_j (T_(j,j+1)+T_(j,j+1)†+B_(j,j+1)).

F is the colored field at the first internal bond. B is the physical
port bond across a cut, acting on the last site of one cell and the
first site of the next. Both preserve cell lengths.

T transfers the first physical site of the right cell to the end of
the left cell, converting H0 to H2 on that site's wrapper and applying
the specified two equalizer phases. Its selector condition is #L to
R#. The inverse is its adjoint. Thus

    T: K_ell tensor K_r -> K_(ell+1) tensor K_(r-1).

The condition on the right selector word makes T zero if the move
would destroy a legal cell of minimum length. Its inverse conditions
are those of the original physical bond. This is a partial isometry
times its fixed amplitude, on the whole direct sum, not just on two
selected codes. Hence ||T||<=u independently of ell,r.

The endpoint and field operators also have norms bounded independently
of ell. Their compensators and visibility factors were fixed using b,
not the actual length of an excited cell. Therefore

    ||F|| <= B_field,
    ||B|| <= B_cut,
    ||T+T†|| <= 2u,

with the full-space constants proved in run9_audit.md and
flat_vacancies.md. There is no factor proportional to the length of
a wrong cell in these perturbation norms. Although V is a sum over
physical bonds, its nonzero perturbation terms occur only at a marker
or the first internal bond following a marker. Internal clock and
wrapper propagation belong to A, not to V.

All terms preserve total length sum ellhat, microscopic total charge,
and M. They are translation invariant in the cell index away from
the rooted-ring seam. At the seam the reference coordinate changes
when the reference marker moves. Keeping x and taking the quotient
above gives the exact operator. Working at fixed physical momentum
instead gives the corresponding boundary twist, which must not be
discarded in spectral comparisons. Local norm bounds do not depend
on that twist.

Consequently the legal-sector bulk problem is an ordinary nearest-
neighbor perturbation of a sum of independent on-cell Hamiltonians,
with a finite-dimensional degenerate ground space at each cell and
uniformly bounded perturbation terms. The total-length constraint
and the global reference quotient remain part of the model.

## 4. What this resolves and what it leaves

[proved] There is no need to control moving cell projectors directly
in physical coordinates to formulate the local elimination. The
exact change of variables supplies a fixed tensor product, an on-cell
gap, and bounded neighboring-cell interactions. This is stronger
than the first-order projected cell description alone.

There are still three distinct tasks.

1. A connected low-frequency elimination theorem is needed for this
   atomic model with its finite on-cell ground degeneracy. The bare-Q
   counterexample in run9_audit applies already to product on-cell
   Hamiltonians, so the atomic reformulation does not validate a global
   inverse or an all-band Schrieffer-Wolff transformation.
2. K is an infinite direct sum before imposing fixed L. A theorem whose
   bounds depend only on the on-cell gap and bounded interaction norms
   could cover it, or one must prove a uniform high-length truncation
   estimate using A_ell's linear growth. An arbitrary finite alphabet
   theorem cannot be invoked without checking this point.
3. The resulting local effective action still needs the Gaussian
   infrared input and the inherited physical ring quotient. Illegal
   selector sectors are invariant and do not occur virtually, but
   their positive dressed defect cost must also be retained when
   identifying the full physical ground state.

The reduction is an exact finite-dimensional statement on every fixed
ring and a well-defined direct-sum statement before fixing its length.
No new numerical calculation is claimed for it. The separate periodic
column tests verify the local transfer and selector rules it uses.
It narrows M to a more familiar atomic problem without claiming to
have solved M or the original conjecture.
