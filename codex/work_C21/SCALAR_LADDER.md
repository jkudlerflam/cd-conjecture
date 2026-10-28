# A local pole ladder whose return interface is scalar

## Rules and code

[proved] An eight-state site consists of a four-state Bell-syndrome
register 0,X,Y,Z and a data qubit. Give each nonzero syndrome energy
Delta and hop it with amplitude -t to a vacant neighbor, preserving
the data. In formula form,

    Q=I-|0><0|,
    T_i=sum_(a=X,Y,Z) |0,a><a,0| tensor I_data,
    H=Delta sum_i Q_i-t sum_i(T_i+T_i^dagger).              (1)

These are all terms. On an open cell the repeated bond is
Delta(Q tensor I+I tensor Q)/2-t(T+T^dagger), supplemented by
Delta Q/2 at each endpoint. Pin the first data qubit to 0 with a
positive unit-energy projector. A cell of m=n+1 sites then has an
exactly 2^n-dimensional vacuum code.

For Delta>2|t|,

    H >= (Delta-2|t|) sum_i Q_i.

The hopping inequality follows by bounding each partial isometry and
its adjoint by its initial and final occupation projectors. This excludes
every nonvacuum sector, not just the single-excitation sector. With the
pin, the complete gap is at least min(1,Delta-2|t|).

Site parity is diag(1,-1,-1,1) tensor Z_data. Physical port operators
are the Pauli matrices on the first qubit of the Bell pair, expressed
in the displayed syndrome basis. They have zero compression and the
standard odd,odd,even grading. The intercell bond is their Pauli exchange.

## Poles and their weights

[proved] The one-excitation Hamiltonian is the m-vertex scalar path,
tensored with the entire syndrome/data space. Its energies and endpoint
weights are

    E_e=Delta-2t cos(pi e/(m+1)),
    w_e=2 sin^2(pi e/(m+1))/(m+1), e=1,...,m.               (2)

All m=n+1 poles are accessible. None of the weights is a logical
operator. The second-order bond is only

    a0 I,  a0=-3 sum_(e,f) w_e w_f/(E_e+E_f).             (3)

Thus a fixed alphabet carrying enough poles does not solve the interface.

## Unitary decorations do not repair it

[proved] More generally consider a rooted tree with scalar on-site
energies, scalar edge magnitudes, and unitary fiber operators U_xy on
its edges. Set W_root=I and define W_y=U_yx W_x recursively away from
the root. There is one path to each vertex. The block-diagonal unitary
with blocks W_x gauges every edge to its scalar magnitude, leaving
the on-site energies scalar. Its Hamiltonian is therefore unitarily
equivalent to A_tree tensor I_fiber. At the root,

    <root|F(H)|root>=<root|F(A_tree)|root> I_fiber          (4)

for every spectral function F. Choosing any other vertex as root gives
the same conclusion for that return source. Cross-vertex propagators
can contain nontrivial unitaries; that does not change (4).

This proof permits arbitrary noncommuting data operations on successive
edges, arbitrary scalar energy ladders, and either sign of the scalar
potentials. It excludes only unitary edge decorations on a tree with
scalar potentials. Loops with nontrivial holonomy, nonunitary edge
weights, and logical-valued on-site energies are outside its hypotheses.
The Clifford cursor uses the last escape.

[numerical] ladder_checks.py also tests a nontrivial unitary decoration.
On the hop to site i+1 it uses exp(i*0.37*A tensor X_(i+1)), where
A exchanges X and Z syndromes. The hopping remains a two-site operator,
is parity even, and leaves the first pinned data bit alone. Diagonalizing
its entire one-cursor sector at n=1,2,3 gives the scalar weights (2).

## Exact finite-coupling control

[numerical] The undecorated model was diagonalized on two and three
cells at n=1,2, and two cells at n=3, with Delta=3,t=0.4. The data
registers are exact spectators. The syndrome space splits by the
bitwise XOR of its Bell labels; the code is in XOR zero. Every hopping
and Pauli-exchange term preserves this label. Restricting to that whole
sector is exact and retains all excitation numbers.

| n | Cells | Full syndrome dimension | XOR-zero dimension | Logical multiplicity |
| --: | --: | --: | --: | --: |
| 1 | 2 | 256 | 64 | 4 |
| 1 | 3 | 4096 | 1024 | 8 |
| 2 | 2 | 4096 | 1024 | 16 |
| 2 | 3 | 262144 | 65536 | 64 |
| 3 | 2 | 65536 | 16384 | 64 |

The physical Hilbert space also contains every data state, including
the pinned bits. Their omission from ED is an exact tensor-factor
reduction. The displayed multiplicity is the full low-band logical
multiplicity. Excited pin sectors remain separated by one unit.
Other XOR sectors start at least Delta-2t above the code. With
||V||<=3(M-1), the tested lambda<=0.08 cannot close this separation.

An independent linear solve for the second-order resolvent agrees with
(3) to at most 1.2e-15. At lambda=0.08,0.04,0.02, the error from
lambda^2(M-1)a0 has successive halving exponents between 3.0021 and
3.0038. Eigenpair residuals are below 1.3e-14.

These checks establish the scalar answer and its cubic remainder.
They cannot validate a nonzero logical hopping bond. Subtracting the
intended nonzero C19 target would instead leave an order-lambda^2
matrix error. No microscopic selector or d0 for the target follows
from the eight-state module count.

[provenance] C22 independently recorded the tree-gauge observation
while its progress file was initially absent here. Its later work
studies a different, charge-conserving bus. No claim of independent
verification of that bus is made in C21.
