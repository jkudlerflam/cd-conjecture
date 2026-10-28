# A charge conserving nearest neighbour circuit for the sparse encoding

## Statement

[proved] Given a simple bipartite graph explicitly as an edge list, the
isometry in [README](README.md) has a uniform, polynomial-size
nearest-neighbour **circuit** implementation. Every elementary gate
commutes with every total species charge. One convenient common local
alphabet has

\[
D=q+3+2(q+1)\lceil\log_2q\rceil=O(q\log q)
\tag{1}
\]

states, including charged data, neutral scratch, and a temporary head.
All scratch and temporary head occupations return to their prescribed
initial state. At the output, only the \(N+4q\) charged data sites of the
two-port encoding carry variable data.

The circuit is exact if algebraic one-qubit rotations on neutral
registers are allowed. With a fixed finite neutral universal gate set,
it approximates the isometry to any operator-norm accuracy \(\delta>0\)
with polynomial gate and workspace bounds, and preserves species
charges **exactly** even at finite approximation accuracy.

The later [EXACT_AMPLITUDES](EXACT_AMPLITUDES.md) result sharpens the
first statement: a constant number of algebraic opcode families, whose
matrices may depend on \(q\), suffices for the exact circuit with the
same \(O(q\log q)\) alphabet order. The
[ORACLE_SCHEDULER](ORACLE_SCHEDULER.md) and
[SWEEP_INTERFACE](SWEEP_INTERFACE.md) supply its reversible moving-head
implementation with only constant controller overhead.

Uniform means that a classical algorithm constructs the circuit from the
edge list, \(q\), and the desired precision in polynomial time. It does
not mean that one repeated Hamiltonian runs it without initialization.
This result supplies no initialization penalties, clock spectrum,
illegal-sector bound, cell selector or output-time port filter.
The circuit writes and may temporarily move the eventual port registers.

The [executable circuit generator](../controller_instantiation/circuit/README.md)
now instantiates both port stages, every Boolean computation, the local
gate lowering, and inverse cleanup. Its neutral auxiliary records start
and finish in the unique blank state: an explicit blank-to-bit0 involution
is applied before use and undone afterwards. This avoids the ambiguity
of verifying an ordinary zero bit using TYPE/species queries alone.
The generator is externally scheduled; its complete autonomous execution
and initialized program sector remain separate requirements.

## 1. The reversible map in the occupation basis

Fix an ordering of the \(N\) edges. For the left stage, the coherent branch
label \(b\) ranges over \(u+N\) values: an upper vertex \(\alpha\), or an
edge \(e=(\alpha,\beta)\). Prepare amplitudes

\[
a=\frac1{\sqrt{2u}}\quad\hbox{on upper vertices},\qquad
b_0=\frac1{\sqrt{2N}}\quad\hbox{on edges}.
\tag{2}
\]

The label lives in neutral binary scratch. It is independent of the
input logical occupation string \(\mathbf n\). Let \(P_e(\mathbf n)\)
be the parity of occupied modes preceding edge \(e\).

The map on each valid branch is:

| Branch | Output port | Output flag | Residual occupation | Sign |
| :--- | :--- | :--- | :--- | :--- |
| upper \(\alpha\) | \(1,\alpha\) | \(g_\alpha\) | \(\mathbf n\) | \(+1\) |
| edge \(e\), \(n_e=0\) | \(0,\beta\) | \(g_\alpha\) | \(\mathbf n\) with \(n_e=1\) | \((-1)^{P_e(\mathbf n)}\) |
| edge \(e\), \(n_e=1\) | \(0,\beta\) | \(h_e\) | \(\mathbf n\) with \(n_e=0\) | \((-1)^{P_e(\mathbf n)}\) |

This is exactly the Kraus formula in README, not an approximation to it.
The inverse is explicit:

* An upper port identifies the upper branch, and the payload is unchanged.
* A lower port \(0,\beta\) with flag \(g_\alpha\) identifies edge
  \((\alpha,\beta)\), which had \(n_e=0\).
* A lower port with flag \(h_e\) identifies edge \(e\), which had \(n_e=1\).

The indicated residual bit is toggled back in the last two cases. Valid
output words therefore recover both branch and input.

For the right stage, use the same basis map. Multiply the creation branch
sign by \((-1)^Q\), and the annihilation branch sign by \((-1)^{Q-1}\),
where \(Q\) is the input occupation number for that stage. These are
precisely the signs of \(f_e^\dagger P\) and \(P f_e\). The inverse
basis map is unchanged.

## 2. Replace formal charge weights by actual species words

The input charged string consists of the \(N\) logical occupation sites
and a fixed reference word containing two copies of each species. Its
length is \(L=N+2q\). The output consists of one port, one flag word of
length \(2q-1\), and \(N\) residual occupation sites. Use a fixed
canonical order for the species in every flag word.

The charge identities in README prove that input and output contain
exactly the same number of every species, separately on every valid
branch. The map is consequently a reversible rearrangement of existing
charged symbols, together with a phase. It creates no charged particle
in neutral scratch. For the two-stage construction reserve two reference
words at the start and apply the same procedure successively to the
residual logical data; total charged length is \(N+4q\).

All basis computations used below have polynomial complexity. Read a
species label in \(\lceil\log_2q\rceil\) bits. A payload occupation is
determined by which of its two permitted species occurs. The edge list
supplies the endpoint labels. The target flag has explicitly known
counts \(2S-e_{1,\alpha}\) or
\(2S+e_{1,\alpha}-2e_{0,\beta}\). Prefix parity and total parity take
linear time in \(N\). There is no enumeration of \(2^N\) occupations.

## 3. Read and route without breaking the symmetry

Introduce these mutually orthogonal local types:

* \(D_s\), one ordinary data type of charge \(e_s\), for each species.
* One blank and two neutral bit types.
* \(H_{j,c}\), a neutral head with bit index \(0\le j<\ell\), bit
  \(c\in\{0,1\}\), where \(\ell=\lceil\log_2q\rceil\).
* \(T_{j,c,s}\), a loaded head of charge \(e_s\).

Their count is (1). All gates below are permutations or neutral
one-qubit gates on subspaces of two adjacent sites, extended by the
identity on unused basis states.

**Species query.** The local reversible gate

\[
|D_s,H_{j,c}\rangle\longmapsto
|D_s,H_{j,c\oplus\operatorname{bit}_j(s)}\rangle
\tag{3}
\]

copies any desired bit of a species label into a neutral head. Its
control is diagonal on the charged site; it preserves every charge.
Copy the head bit into neutral working memory and apply (3) again to
clear the head. Neutral head-index changes are ordinary local basis
permutations.

**Controlled adjacent data transposition.** A three-register
controlled swap need not be assumed as an elementary gate. Place the
neutral control head between the two data registers. Use three
nearest-neighbour gates:

\[
\begin{aligned}
|D_s,H_{j,c},D_t\rangle
&\longleftrightarrow|\mathrm{blank},T_{j,c,s},D_t\rangle,\\
|\mathrm{blank},T_{j,1,s},D_t\rangle
&\longmapsto|\mathrm{blank},T_{j,1,t},D_s\rangle,\\
|\mathrm{blank},T_{j,c,s'},D_{t'}\rangle
&\longleftrightarrow|D_{s'},H_{j,c},D_{t'}\rangle.
\end{aligned}
\tag{4}
\]

The middle gate is the identity for \(c=0\). The first and last use the
same load/unload involution on the left pair. Equation (4) swaps the two
data symbols exactly when \(c=1\), and returns the head and blank
occupancies to their initial form. Each constituent gate conserves
every total charge. A fixed local gate (the middle line) is its own
inverse, including when \(s=t\).

Unconditional swaps of complete neighbouring site states also conserve
all species charges. They bring specified registers together and
restore their positions, so arbitrary finite circuit layouts can be
compiled to a line with a polynomial overhead. Two-site neutral
controlled gates, one-site neutral rotations, and exact neutral Toffoli
decompositions complete the ordinary reversible computation toolkit.
All such gates act on charge-zero states.

The loaded head is essential to this explicit two-site compilation.
Merely invoking a three-site controlled swap would leave a locality gap.

## 4. Clean reversible evaluation and cleanup

Here is a constructive reversible procedure.

1. Query the input species word \(X\) into neutral scratch. Retain the
   coherent branch label.
2. Reversibly compute the target word \(Y\) and its sign using the table
   in Section 1. Compute a sequence of adjacent transpositions that maps
   \(X\) to \(Y\): process target positions from left to right, select
   the first matching remaining source symbol, and move it left by
   adjacent swaps. Duplicated species are matched stably. This uses at
   most \(L(L-1)/2\) swaps.
3. Apply the swaps to the actual charged data using (4), with control
   bits in neutral scratch. Apply the computed sign by a neutral phase
   gate. The charged data are now \(Y\).
4. Reverse the neutral work that computed the swap plan, target word
   and sign, retaining just the original input-word copy and branch.
   This reversal uses the saved input and branch, so it does not undo
   the physical routing.
5. Query the new physical word \(Y\), compute its explicit inverse
   from Section 1, and XOR the recovered input and branch into the
   saved registers. They become zero. Reverse this inverse computation
   and its queries, using \(Y\), to clear the remaining work.

Each classical evaluation can be made reversible by retaining the
inputs and intermediate Boolean-operation outputs, copying out its
result, and reversing its gates. Thus this description is an algorithm
for a circuit, not an assumption that an irreversible classical map is
unitary. The decoder is used only on valid output words when proving
the isometry; it can be assigned any total Boolean extension elsewhere.
The complete circuit is a unitary on the entire enlarged space because
each constituent gate is unitary.

The resulting action on the valid input subspace is the desired
isometry, with branch and scratch reset. Since the proof applies to
every occupation basis state with its phase, it also applies coherently
to every superposition.

## 5. Amplitudes and resource bounds

For exact preparation of (2), use a binary tree over the \(u+N\) valid
labels. Each leaf has rational probability \(1/(2u)\) or \(1/(2N)\).
At a prefix, rotate the next neutral bit with amplitudes equal to the
square roots of its two child probability ratios. Empty subtrees need
no rotation. Controlled rotations and reversible prefix comparisons
give \(O((u+N)\log(u+N))\) elementary neutral gates with work bits; no
branch amplitude depends on the payload.

These rotations have algebraic matrix entries, but generally are not
exact words in a chosen finite universal gate set. This distinction
matters for a history-machine implementation. With a fixed inverse-closed
neutral universal set, the Solovay–Kitaev theorem supplies
polylogarithmic approximation overhead for each rotation
([Dawson and Nielsen, primary source](https://arxiv.org/abs/quant-ph/0505030)).
Allocate error \(\delta/G\) to each of \(G\) gates requiring approximation;
the telescoping product bound gives total circuit error at most
\(\delta\). All approximating gates remain exactly charge preserving,
because they act solely on neutral states. An approximation may create
a small component on invalid branch labels; applying the same exact
unitary transducer preserves the norm of that error.

Conservative explicit polynomial bounds suffice. Let
\(L=N+4q\) cover both stages. Producing and reversing a full Boolean
history of the graph lookup, target generation and stable routing plan
takes \(T=O(L^3\log L)\) gates and work bits. This overcounts simple
implementations but follows from at most three nested loops of length
\(L\), each using comparisons and updates on \(O(\log L)\) bits.
Reversible computation, plan execution and cleanup use \(O(T)\)
abstract-register gates and sites. Compiling every gate to nearest
neighbours by register swaps costs at most another factor \(O(T)\).
Thus one sufficient bound is

\[
\begin{array}{ll}
\text{neutral workspace:}&O(L^3\log L),\\
\text{nearest-neighbour gate count:}&O(L^6\log^2 L),
\end{array}
\tag{5}
\]

with an additional polylogarithmic factor in \(L/\delta\) if the
neutral rotations are compiled to a fixed finite universal set.
These are bounds on circuit length and cell length, not on local
alphabet dimension; (1) is independent of both.

If \(\|\widetilde W-W\|\le\delta\), then for any unit-norm port operator

\[
\|\widetilde W^\dagger F\widetilde W-W^\dagger F W\|\le2\delta.
\tag{6}
\]

For a unit-norm physical intercell link, the corresponding two-cell
compression error is at most \(4\delta\). Relative to the desired
logical hopping scale \(\eta\), choose \(\delta\ll\eta\).
No symmetry error is introduced by this choice.

## 6. The separate Hamiltonian contract

A circuit is scheduled externally and begins in the initialized domain.
A translation-invariant nearest-neighbour Hamiltonian must establish
those conditions energetically on every chain length. In particular,
it must enforce the neutral reference and branch preparation, the
charged reference words, valid graph/program data, a single legal
history, and an appropriate endpoint filter so that an intercell port
bond sees the desired code. It must also control competing illegal and
shorter history components and combine with the cell selector.

These obligations are addressed in the later
[initialized-core composition](../torus_bottleneck/EXISTENTIAL_CORE.md),
using the exact gate menu and interpreter linked above. The circuit
result itself removes general isometry synthesis and charge-conserving
gate universality as separate objections; its proof does not replace
the Hamiltonian construction or the finite-coupling infrared analysis.
