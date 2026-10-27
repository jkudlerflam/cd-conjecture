# Local vacancy extension and an obstruction to exact cat padding

Status: finite-dimensional theorems, proved below. The explicit local
example uses a ferromagnetic code and has no large-c claim. The final
spin-charge discussion is conditional.

## 1. The exact cat padding cannot be an isolated local ground code

There is a necessary qualification to `padding_transport.md`. Its exact
two-placement code has no strictly local parent with only that code as
ground space, once the core is longer than the interaction range.

Let A_p=V tensor |0>^p and B_p=|0>^p tensor V. These isometries have
orthogonal occupation supports. Suppose every Hamiltonian term omits at
least one site where the two supports have different occupation. Then

    A_p^dagger H B_p = 0.                                  (1)

For the proposed p=1,2 padding and nearest-neighbor terms, k>=2 suffices:
already the first and last sites disagree, and no bond contains both.
Matrix elements vanish because at an untouched disagreeing site one
branch is vacant and the other occupied.

Suppose every vector W_p psi=(A_p psi+B_p psi)/sqrt(2) is a ground vector
of H with ground energy E. Positivity of H-E gives

    0 = <W_p psi|(H-E)|W_p psi>
      = ( <A_p psi|(H-E)|A_p psi>
          +<B_p psi|(H-E)|B_p psi> )/2.                    (2)

Each summand is nonnegative and must vanish. A positive matrix with zero
quadratic form annihilates the vector. Hence every A_p psi and every
B_p psi separately belongs to the ground space. The degeneracy is at
least 2D, not D. The antisymmetric combinations are also ground vectors.

This argument does not assume translation invariance, frustration
freeness, or commuting terms. The pure two-placement construction gives
coherent projected transport algebraically, but a local realization must
include intermediate vacancy positions or some other dressing. The next
construction does so.

## 2. A local extension theorem for any frustration-free core

Assume the k-site core has a uniform open-chain parent

    H_core = sum_(j=1)^(k-1) h0_(j,j+1),  h0>=0,
    ker(H_core)=range(V).                                 (3)

The positivity of each summand implies h0_j V=0. This is the extra
hypothesis. A general frustrated parent does not satisfy it.

Enlarge each site from C^q to C|0> direct-sum C^q. Write n=|0><0| for
vacancy number and extend h0 by zero outside the two-occupied-site sector.
On N=k+p sites, for p=1 or 2, use the SAME repeated bond

    h_ext = U(n tensor I + I tensor n - 2n tensor n)
            -w sum_(a=1)^q (|0,a><a,0|+|a,0><0,a|)
            +g h0,          U,w,g>0.                     (4)

Restrict to exactly k occupied sites. The core's q-state labels cannot
pass each other under vacancy hopping on an open chain. Therefore this
sector has the exact ordered-coordinate decomposition

    C^(binomial(N,p)) tensor (C^q)^(tensor k),             (5)

where the first factor specifies hole positions and the second lists the
occupied labels in their physical order. This is a unitary change of
basis, not an approximation or a low-energy restriction.

The first two terms of (4) become H_charge tensor I. The charge matrix
has diagonal energy U times the number of occupied-vacant interfaces
and off-diagonal -w for every allowed nearest-neighbor hole hop. Its
configuration graph is connected. Perron-Frobenius gives a unique,
strictly positive normalized ground vector chi_p and energy e_p. Reflection
invariance makes chi_p reflection symmetric.

Let D_j be the diagonal charge projector onto configurations in which
the j-th and (j+1)-th occupied labels are physically adjacent. The full
Hamiltonian in (5) is

    H_ext = H_charge tensor I + g sum_j D_j tensor h0_j.   (6)

Since all terms in the second sum are positive, H_ext>=e_p. The states

    V_p psi = chi_p tensor V psi                          (7)

attain e_p, because every h0_j annihilates V.

Conversely, a ground state must minimize the first term and annihilate
the second, so it has the form chi_p tensor phi. Since chi_p is strictly
positive, <chi_p|D_j|chi_p> is strictly positive for every j. Its energy
above e_p is zero only if every <phi|h0_j|phi>=0, hence phi is in the
core kernel. Thus (7) is the COMPLETE ground space in the fixed-k sector.

For every fixed k,p,U,w,g>0 the gap above it is positive, by finite
dimension and the exact kernel characterization. This is not a bound
uniform in k. Its small-w splitting can be very small, which matters
when coupling the cells. The theorem makes no claim about the ordering
of different particle-number sectors; a full selector must control them.

The construction preserves an arbitrary frustration-free core code, not
just a classical or product code. It uses a single bulk bond at both
lengths, with no length-dependent endpoint term.

## 3. Endpoint maps and scalar boundary transport

Extend a core one-site observable F by zero on the vacancy. If physical
site 1 is occupied, it contains the first ordered label, regardless of
all hole positions. A number-preserving endpoint observable leaves those
positions unchanged. Consequently

    V_p^dagger F_first V_p = r_p A,
    V_p^dagger F_last V_p = r_p B,
    r_p = 1-Prob_(chi_p)(first site vacant),               (8)

where A,B are the original core endpoint compressions. Equality of the
two coefficients uses reflection of chi_p, not reflection of the core.

The leading boundary-shift overlap is also scalar for arbitrary V.
On its input, the left cell has p=1 and the right p=2; after shifting
their cut right, they have p=2 and p=1. To keep k occupied labels in
each cell, the transferred site must be vacant. No ordered logical label
crosses the cut. Every surviving core overlap is V^dagger V=I. Thus

    (V_2 tensor V_1)^dagger (V_1 tensor V_2)
        = beta I_(D tensor D),                           (9)

with beta=a_left a_right>0. For N1=k+1, N2=k+2,

    a_left = sum_(s=0)^(N1-1) chi_1({s}) chi_2({s,N1}),
    a_right = sum_(s=0)^(N1-1) chi_1({s}) chi_2({0,s+1}). (10)

All amplitudes are real and positive. Reflection gives a_left=a_right.
These formulas are in the physical spectator coordinates after the cut
has moved. They do not assert that the microscopic hop has no leakage.

For fixed k, as w/U tends to zero the charge states approach the equally
weighted superposition of the two endpoint hole droplets. To check this,
at w=0 these are the only configurations with one interface, with energy
U; all others have at least two interfaces. In the reflection-even sector
the ground vector at w=0 is unique, and the gap there is at least U.
Finite-dimensional perturbation theory gives an excited-configuration
amplitude O(w/U). Reflection fixes the two endpoint weights to be equal.
It follows that

    r_1,r_2 = 1/2 + O((w/U)^2),
    beta -> 1/4.                                         (11)

Constants can depend on k. The odd reflection partner is separated only
by a high-order tunneling splitting; it is not controlled by the U gap
in the even sector. The coupling between cells must be below this actual
gap, not merely below U.

Unlike exact cat padding, the local model generally has r_1 != r_2 at
finite w. Hence its payload bond is r_p r_p' times the core bond. The
length/payload coupling can be made small at fixed k, but is not zero.
Its infrared effect requires analysis. No tuning that makes r_1=r_2
with this same bond has been proved.

## 4. Explicit check with a quantum core

`vacancy_parent.py` takes q=2, k=3, and h0=(I-swap)/2. The core kernel is
the four-dimensional symmetric spin-3/2 representation. It is quantum:
the compressed X,Y,Z endpoint operators do not commute. This example
does not have an extensive Clifford family or a large target charge.

For w=0.15,0.3,0.5, U=1, g=0.7 the script builds the repeated physical
bond (4) on both lengths, diagonalizes each fixed-k sector, and checks
the entire ground code (7), its degeneracy, the endpoint maps (8), and
the two-cell transition (9). The data are in `vacancy_parent_checks.json`.

Run `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 vacancy_parent.py`.

## 5. Why this still does not finish the conjecture

The theorem reduces two-length local compatibility to the one-length
problem in the frustration-free subclass. No uniform nearest-neighbor
frustration-free core with the required many-channel endpoint maps has
been constructed here. The earlier Choi code has an arbitrary parent
projector and cannot be inserted into (3) on that basis alone. Its six
tested local embeddings already failed the parent test in the prior run.

The contribution on programmable cells proves a useful finite-history
interface, but does not provide the microscopic repeated-bond core
needed by (3), with its required endpoint compressions. Neither that
interface nor this extension establishes a stable large-c chain.

Vacancy doping of an entire hypothetical counterexample is another
conditional idea. On C|0> direct-sum C^d, nearest-neighbor particle-vacancy
hopping preserves the order of all occupied labels. At large hopping
compared with occupied-label interactions, first-order charge projection
gives a squeezed label chain. If spin-charge separation and infrared
stability held, it would add one charge mode, taking (d,c) to (d+1,c+1)
and preserving c-(d-1). Its particle number could adjust to a multiple
of a preferred cell period on every large physical ring.

This is not an all-ring completion theorem. The charge sector is gapless,
so finite-volume perturbation theory has denominators O(1/L) and cannot
be used uniformly as L grows. Exact factorization (5) at zero label
interaction does not make the interacting squeezed Hamiltonian exact.
The earlier direct-sum contribution already uses this mechanism at the
RG level. Without an additional thermodynamic stability input, it does
not close the present construction, and no further numerical claim is
made for it.
