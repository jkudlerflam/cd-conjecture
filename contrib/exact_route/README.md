# Sparse fermionic encodings for two physical ports

## Claim and status

[proved] For every finite simple bipartite graph with \(u\) upper vertices,
\(v\) lower vertices and \(N>0\) edges, a physical port of dimension
\(q=u+v\) can compress its edge matrix units to \(N\) logical fermion
creation operators with visibility \(1/(2\sqrt{uN})\). Two distinct ports
can simultaneously compress to the left and parity-dressed right
operators needed for nearest-neighbour fermion hopping. The maps are
exact, preserve the full species torus, and have a sparse explicit
description; no Choi matrix square root or general isometry synthesis is
needed.

[proved] A [uniform charge-conserving nearest-neighbour circuit](CIRCUIT.md)
implements this sparse map with \(O(q\log q)\) data and head types and
polynomial gate and workspace bounds. It is exact with algebraic
neutral rotations, or arbitrarily accurate with a fixed finite neutral
universal gate set while conserving all charges exactly.

[proved] [Exact indexed amplitudes](EXACT_AMPLITUDES.md), a
[reversible oracle scheduler](ORACLE_SCHEDULER.md), and an explicit
[ordered-sweep interface](SWEEP_INTERFACE.md) implement the circuit
with a fixed number of opcode/control families and the same alphabet
order. The initialized history Hamiltonian is supplied separately by
[the companion core composition](../torus_bottleneck/EXISTENTIAL_CORE.md).

These are finite-cell construction results. The complete
[candidate counterexample](../research_continuation/PROPOSED_COUNTEREXAMPLE.md)
also uses the mobile wrapper and the stated infrared physics inputs.
The algebraic identities and finite checks here alone do not establish
a counterexample to \(c\le d-1\).

**Provenance and checking.** Derived on 2026-10-27 with AI assistance from
Codex GPT-6. A separate Codex agent independently checked the Kraus
normalization, the channel attenuation and the species-charge assignments.
`verify.py` performs independent dense finite-dimensional checks. Neither
these calculations nor the proof have received human review.

## 1. An arbitrary graph and one port

Let \(G=(U,V,E)\) be a finite simple bipartite graph, with
\(|U|=u\ge1\), \(|V|=v\ge1\), \(|E|=N>0\). Isolated vertices are allowed.
Assign one logical complex fermion \(f_e\) to each edge. Write
\(\mathcal F_N\) for their Fock space and \(P=(-1)^Q\) for total parity.
The physical port is

\[
B=\operatorname{span}\{|1,\alpha\rangle:\alpha\in U\}
\oplus\operatorname{span}\{|0,\beta\rangle:\beta\in V\}.
\]

Its edge matrix units are
\(F_e=|1,\alpha\rangle\langle0,\beta|\),
where \(e=(\alpha,\beta)\).

Take an auxiliary flag space \(A\) of dimension \(u+N\), with orthonormal
flags \(g_\alpha\) and \(h_e\). For positive \(a,b\) satisfying

\[
u a^2+N b^2=1,
\tag{1}
\]

define maps \(K_p:\mathcal F_N\longrightarrow A\otimes\mathcal F_N\) by

\[
\begin{aligned}
K_{1,\alpha}&=a|g_\alpha\rangle\otimes I,\\
K_{0,\beta}&=b\sum_{e=(\alpha,\beta)\in E}
\bigl(|g_\alpha\rangle\otimes f_e^\dagger
       +|h_e\rangle\otimes f_e\bigr).
\end{aligned}
\tag{2}
\]

Set

\[
V_L\psi=\sum_{p\in B}|p\rangle\otimes K_p\psi.
\tag{3}
\]

The flags in a fixed \(K_{0,\beta}\) are all distinct. The CAR therefore
give

\[
\sum_pK_p^\dagger K_p
=u a^2I+b^2\sum_{e\in E}(f_e f_e^\dagger+f_e^\dagger f_e)
=I.
\tag{4}
\]

Thus \(V_L\) is an exact isometry, with

\[
V_L^\dagger(F_e\otimes I)V_L=K_{1,\alpha}^\dagger K_{0,\beta}
=ab f_e^\dagger.
\tag{5}
\]

For a pair \((\alpha,\beta)\notin E\) this compression is zero. This
also shows why using a graph rather than all \(uv\) matrix units causes
no difficulty.

The coefficient \(ab\) is maximized within (1) at

\[
a^2=\frac1{2u},\qquad b^2=\frac1{2N},\qquad
\epsilon:=ab=\frac1{2\sqrt{uN}}.
\tag{6}
\]

This is an optimum for this explicit one-parameter construction, not a
universal optimum over all port encodings. Compared with the sufficient
generic-Choi choice \(1/(4qN)\) in the imported complete-bipartite
construction, (6) gives a larger rigorously achieved visibility.

## 2. The residual logical channel

The Heisenberg channel on the residual Fock factor is

\[
\mathcal D_L(O)=V_L^\dagger(I_{B\otimes A}\otimes O)V_L
=u a^2 O+b^2\sum_e(f_e O f_e^\dagger+f_e^\dagger O f_e).
\tag{7}
\]

For \(O=f_j^\dagger P\), every term with \(e\ne j\) equals \(O\),
since \(O\) commutes with \(f_e,f_e^\dagger\). The \(e=j\) term is zero
by nilpotence. Hence

\[
\mathcal D_L(f_j^\dagger P)
=[u a^2+(N-1)b^2]f_j^\dagger P
=(1-b^2)f_j^\dagger P.
\tag{8}
\]

For comparison, the undressed odd operator obeys

\[
\mathcal D_L(f_j^\dagger)
=[u a^2-(N-1)b^2]f_j^\dagger
=[1-(2N-1)b^2]f_j^\dagger.
\tag{9}
\]

Distinguishing (8) from (9) is essential to the two-port construction.

## 3. Two disjoint ports

Define another set of CAR generators on the same logical Fock space by

\[
g_e^\dagger=f_e^\dagger P,\qquad g_e=P f_e.
\tag{10}
\]

They have the same adjoints, anticommutators and number operators as
the original generators. Build \(V_R\) from (2) using \(g_e\) instead of
\(f_e\), with an independent port \(B_R\) and flag space \(A_R\). Apply
\(V_L\) first, and then apply \(V_R\) to its residual Fock factor:

\[
W=(I_{B_L\otimes A_L}\otimes V_R)V_L.
\tag{11}
\]

The order of the output factors is \(B_L,A_L,B_R,A_R,\mathcal F_N\).
Since \(V_R\) is an isometry, later encoding cannot change the already
encoded left port. Equation (8) determines the right port. Thus

\[
\boxed{
W^\dagger F_e^{(L)}W=\epsilon f_e^\dagger,
\qquad
W^\dagger F_e^{(R)}W=\epsilon(1-b^2) f_e^\dagger P.
}
\tag{12}
\]

The two physical port algebras commute. Their compressed images need
not commute: an isometric compression is not a homomorphism. Equation
(12) has been checked directly from the full joint isometry, independently
of the channel derivation.

Using different normalized choices for the two stages is possible.
With identical choices (6), the product of port visibilities is

\[
\eta=\epsilon^2(1-b^2)
=\frac1{4uN}\left(1-\frac1{2N}\right)>0.
\tag{13}
\]

On a chain of cells, the physical bond

\[
h_{\rm port}=\sum_{e\in E}
\bigl(t_e F_{e,j}^{(R)}(F_{e,j+1}^{(L)})^\dagger
 +\overline t_e(F_{e,j}^{(R)})^\dagger F_{e,j+1}^{(L)}\bigr)
\tag{14}
\]

compresses **exactly at first order** to

\[
\eta\sum_e\bigl(t_e f_{j,e}^\dagger P_j f_{j+1,e}
+\overline t_e f_{j+1,e}^\dagger P_j f_{j,e}\bigr).
\tag{15}
\]

With \(c_{j,e}=(\prod_{i<j}P_i)f_{j,e}\), this is ordinary
nearest-neighbour hopping
\(\eta\sum_e(t_e c_{j,e}^\dagger c_{j+1,e}+\mathrm{h.c.})\).
On a ring, the conventional parity-dependent fermion boundary condition
must be retained. Equation (15) says nothing about higher-order leakage
from a finite-gap microscopic cell.

There is no channel-count factor in the physical norm of (14).
For a simple graph, different edges use disjoint pairs of two-port
basis states:

\[
|0,\beta;1,\alpha\rangle
\ \longleftrightarrow\
|1,\alpha;0,\beta\rangle.
\tag{15a}
\]

Thus (14) is a direct sum of \(2\times2\) exchange blocks of eigenvalues
\(\pm|t_e|\), together with zero blocks, and
\(\|h_{\rm port}\|=\max_e|t_e|\). Choosing
\(t_e=-1/(2\cos k_e)\) with \(\pi/4<k_e<\pi/3\) gives
\(\|h_{\rm port}\|<1\), independently of \(N\).

## 4. Exact species covariance

Give each physical species its standard charge vector
\(e_{1,\alpha}\) or \(e_{0,\beta}\), and give logical creation \(f_e^\dagger\)
the charge

\[
w_e=e_{1,\alpha}-e_{0,\beta}.
\tag{16}
\]

Let \(S=\sum_{s=1}^q e_s\) and \(C_0=2S\). Assign the auxiliary flags
the charges

\[
\xi(g_\alpha)=C_0-e_{1,\alpha},\qquad
\xi(h_{\alpha\beta})=C_0+e_{1,\alpha}-2e_{0,\beta}.
\tag{17}
\]

Every nonzero term in (2) has total output charge equal to its input
logical charge plus \(C_0\). For example, the two lower-port branches
have charge increments

\[
e_{0,\beta}+\xi(g_\alpha)+w_e=C_0,
\qquad
e_{0,\beta}+\xi(h_e)-w_e=C_0.
\tag{18}
\]

The upper-port branch satisfies the same identity. Thus \(V_L\)
intertwines the complete \(U(1)^q\) species torus up to a fixed overall
charge. Since \(P\) is neutral, the same is true for \(V_R\), and \(W\)
has the fixed offset \(2C_0=4S\).

The upper-species count modulo two is logical fermion parity. Both
flag types in (17) have odd upper-species count, while \(C_0\) has even
upper-species count. The three branches therefore preserve parity
exactly, with even baseline offsets; no grading convention is suppressed.

All flag charge vectors have nonnegative integer entries summing to
\(2q-1\), so each flag is realized by one fixed-order word of that length
in the **same** \(q\)-species alphabet. Different flags have different
charge vectors and hence orthogonal words. Encode the logical occupation
of edge \((\alpha,\beta)\) in one further data site by

\[
|0\rangle_e\mapsto|0,\beta\rangle,
\qquad |1\rangle_e\mapsto|1,\alpha\rangle.
\tag{19}
\]

The original logical reference charge is
\(R=\sum_{\beta\in V}\deg(\beta)e_{0,\beta}\).
Consequently the complete two-port isometry \(W\) has a representation
inside \(N+4q\) ordinary data sites, with total charge

\[
R+4S+\sum_e n_e w_e.
\tag{20}
\]

This is a representation embedding and a word-length count, not a
nearest-neighbour Hamiltonian or an operational clock construction.

Total physical upper-species charge on these data sites restricts
exactly to \(Q+4u\). A uniform upper-species field therefore gives the
desired common logical chemical potential, plus a constant \(4u\) per
cell. A mobile-cell implementation may absorb that constant in its
marker chemical potential. Any additional charged controller reference
sites would contribute their own fixed offset, which must also be
included; the neutral workspace in [CIRCUIT](CIRCUIT.md) contributes none.

## 5. A sparse reversible description

For an occupation state \(|\mathbf n\rangle\), the one-port isometry
has only \(u+N\) nonzero branches:

* \(u\) upper branches leave the logical occupation unchanged, each
  with amplitude \(a\).
* Each edge contributes exactly one lower branch of amplitude magnitude
  \(b\): an empty mode is created and uses flag \(g_\alpha\), whereas an
  occupied mode is annihilated and uses flag \(h_e\). The sign is the
  usual Jordan–Wigner prefix parity.

The branch is recoverable from the output port, auxiliary flag, and
logical occupations. Thus this is a coherent reversible map after a
state-independent branch preparation. At (6), branch preparation is
an equal superposition of an upper-vertex register and an edge register,
each of total probability \(1/2\), with uniform amplitudes inside each
register. The right stage additionally uses the explicitly known total
parity sign from (10).

There are \(u+2N\) nonzero operator-valued Kraus entries at one port;
the joint encoding contains at most \((u+2N)^2\) products. Each product
involves at most two creation or annihilation operators and, for the
right stage, a parity operator. This replaces an unconstrained
exponentially large matrix description by a short, structured rule.
[CIRCUIT](CIRCUIT.md) implements coherent graph lookup, the parity scan,
charge-conserving rearrangement and scratch cleanup with
\(O(q\log q)\) local data and head types.
[EXACT_AMPLITUDES](EXACT_AMPLITUDES.md) gives an exact finite-cardinality
algebraic opcode menu. [ORACLE_SCHEDULER](ORACLE_SCHEDULER.md) proves an
existential reversible neutral interpreter with constant control
overhead, and [SWEEP_INTERFACE](SWEEP_INTERFACE.md) embeds each
interpreter step in one ordered sweep of a fixed two-site gate.
The [initialized-core composition](../torus_bottleneck/EXISTENTIAL_CORE.md)
addresses its legal history and endpoint readout as separate obligations
from these circuit/compiler lemmas.

The construction applies directly to the projective-plane incidence
graphs in [the graph selection-rule result](../torus_bottleneck/README.md).
For a plane of order \(p\), put \(s=p^2+p+1\). Then \(u=v=s\),
\(q=2s\), \(N=(p+1)s\), and

\[
\epsilon=\frac1{2s\sqrt{p+1}},\qquad
\eta=\frac{1-1/[2(p+1)s]}{4(p+1)s^2}.
\tag{21}
\]

The graph's absence of four-cycles suppresses flavour-changing quartic
operators independently of this encoding. The code and that selection
rule solve different obligations and can be combined.

## 6. Reproduction and scope

Run `python3 contrib/exact_route/verify.py`. The script requires numpy,
uses no network or external data, and writes `verification.json` next to
itself. It checks a single edge, a three-edge path, a four-cycle and a
graph with an isolated vertex. The largest Fock space has four modes.
It constructs the full joint isometry for every case and tests its
normalization, both port compressions, all missing-edge compressions,
the two channel identities, parity, and every species charge.

It also exhausts 376 occupation/branch inputs across the left and right
maps, checking the actual charged output words, the explicit inverse and
stable adjacent routing. At \(q=2,3,4,5\), every local load, exchange
and query gate is checked as a charge-preserving involutive permutation
on its complete two-site alphabet, and all 258 controlled-swap inputs
are checked through the three nearest-neighbour gates. The physical
link norm is verified separately. The recorded run took about
0.03 seconds; the largest floating-point error was
\(6.7\times10^{-16}\), and all discrete charge, inversion and routing
checks were exact. These tests do not claim an implemented generic
circuit compiler or an initialized history Hamiltonian.

Run sweep_checks.py for the fixed-gate moving-head interface: it tests
the full two-site unitary with complex number-conserving opcode blocks,
272 complete ordered sweeps, and loaded-token exchanges at \(q=2,3,5\).
Every reported discrepancy is zero. It explicitly includes alphabet
states with no head or multiple heads when checking the two-site gate.

Run amplitude_checks.py for the exact indexed-rotation construction:
256 uniform intervals and 1,024 upper/edge branch preparations pass,
with maximum error \(2.0\times10^{-15}\). At most one exceptional
rotation per prefix depth is used, as in the analytic proof.

Run helper_checks.py for the verifier's trusted boundary resource:
29,359 basis-query round trips pass exactly at \(q=2,3,4\), including
18,364 queries targeting malformed \(H\) or \(T\) records. Every
record is restored, the helper returns to its initialized state, and
every species count is preserved.

These are finite-dimensional consistency checks of the general proof;
they are not central-charge estimates. The sparse exact encoding is a
new constructive ingredient for the first-order route in
[the imported torus proposal](../chat_followup/torus_encoding.md).
Its use in the [candidate microscopic construction](../research_continuation/PROPOSED_COUNTEREXAMPLE.md)
must be assessed together with that construction's initialization,
mobile-wrapper and finite-coupling infrared arguments.
