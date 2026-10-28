# Exact branch preparation with a fixed number of opcode families

## Claim and scope

[proved] The sparse isometry has an exact circuit whose elementary gate
menu has constant cardinality, although the numerical matrices in that
menu depend on the graph parameter. Its local alphabet still has
\(O(q\log q)\) states. Thus neither a precision parameter nor a growing
number of machine control states is needed for the existential
history-machine application.

This does not claim exact synthesis from a fixed numerical gate set such
as Clifford+T. The allowed Hamiltonian family already depends on \(q\):
for example, the species-query matrix does. Including a block-diagonal
neutral rotation with algebraic entries uses the same form of dependence.

## 1. Uniform intervals have few exceptional rotations

Let \(1\le M\le 2^\ell\). We prepare

\[
|u_M\rangle=M^{-1/2}\sum_{x=0}^{M-1}|x\rangle
\tag{1}
\]

from \(\ell\) zero bits, from most to least significant bit. For a prefix
\(x\) of length \(r\), let its interval have width \(2^{\ell-r}\).
Its intersection with \([0,M)\) is either full, empty, or partial.
There is at most one partial prefix at each depth: the prefix containing
the upper boundary of the interval.

At a full prefix apply a Hadamard to the next zero bit. At a partial
prefix put

\[
t_r=M\bmod 2^{\ell-r},\qquad
A_r=\min(t_r,2^{\ell-r-1}),\qquad B_r=t_r-A_r.
\tag{2}
\]

Here the partial case means \(t_r\ne0\). Apply

\[
R_{M,r}=
\begin{pmatrix}
\sqrt{A_r/t_r}&-\sqrt{B_r/t_r}\\
\sqrt{B_r/t_r}&\sqrt{A_r/t_r}
\end{pmatrix}.
\tag{3}
\]

An empty prefix has zero amplitude and may receive the identity.
An induction on the prefix length proves that the squared amplitude of
every prefix is its number of valid completions divided by \(M\).
After the last step all \(M\) valid labels have amplitude \(M^{-1/2}\).
The construction uses at most \(\ell\) exceptional matrices (3).
In particular \(B_r=0\) gives the identity; no singular division occurs
because (3) is used only for a nonempty partial prefix.

The predicates identifying full and partial prefixes are Boolean
comparisons of \(r\) bits with known integers. Compute each predicate
reversibly, apply its controlled gate, then uncompute it. These predicates
depend on the earlier prefix bits, which that gate does not change.
The schedule is fixed even though its controls are coherent.

## 2. Upper and edge branches

Set \(\ell=\lceil\log_2\max(u,N)\rceil\), with the zero-bit convention when
\(\max(u,N)=1\). Prepare one branch-type bit in \(|+\rangle\).
Conditional on its value, prepare (1) with \(M=u\) or \(M=N\), using a
common \(\ell\)-bit label register. The resulting amplitudes are exactly
\(1/\sqrt{2u}\) and \(1/\sqrt{2N}\), as required in CIRCUIT.

There are at most \(2\ell\) exceptional rotation matrices across the two
cases. The same menu is reused for both ports. A graph with isolated
vertices causes no difficulty, since the upper branch is present for
every upper vertex.

## 3. A constant menu, including the index selector

Choose

\[
m=\max\{1,\lceil\log_2q\rceil,2\ell\}.
\tag{4}
\]

Extend CIRCUIT's neutral \(H_{j,c}\) and loaded \(T_{j,c,s}\) index range
to \(0\le j<m\). The alphabet count becomes

\[
D_{\rm exact}=q+3+2(q+1)m.
\tag{5}
\]

Since a simple bipartite graph has \(N\le q^2/4\) and \(u\le q\),
\(m=O(\log q)\). Assign the at most \(2\ell\) rotation matrices to fixed
index values and assign the identity to any unused value.

One controlled-rotation opcode acts on a neutral head and a neutral bit:

\[
|H_{j,c},d\rangle\longmapsto
\begin{cases}
|H_{j,c},0\rangle,&d=0,\\
\sum_{c'}(R_j)_{c'c}|H_{j,c'},1\rangle,&d=1.
\end{cases}
\tag{6}
\]

Its domains for different \(j\) are orthogonal. Extend it by the identity
on other two-site basis states. It is a unitary and conserves all species
charges on the entire local Hilbert space.

A single SHIFT opcode sends \(j\) to \(j+1\bmod m\), fixing \(c,s\).
A neutral type-conversion permutation exchanges ordinary bit \(c\)
with \(H_{0,c}\). Consequently a target bit can be converted, shifted
to any desired index, rotated using (6), shifted back and converted to
an ordinary bit. Every index selection uses a fixed sequence of at most
\(m\) shifts determined by the circuit description. A predicate on the
quantum data changes only the action of (6), never this sequence.

The species-query opcode is defined as before for
\(j<\lceil\log_2q\rceil\), and as the identity for the additional index
values. LOAD and EXCHANGE extend over all \(j<m\). These are still single
block-diagonal opcode families; there is no separate finite-control
state for every index or species. Inverses add at most a constant factor
to the opcode count.

Neutral controlled-H and Boolean gates, together with this indexed
rotation, implement Section 1 exactly. Thus the enlarged finite menu
implements all of CIRCUIT with algebraic matrix entries and no
approximation error. Exact amplitudes replace only its branch preparation;
the reversible transducer and species routing remain unchanged.

## 4. Verification and limits

The accompanying amplitude_checks.py verifies the prefix preparation for
every \(M\le256\), and verifies the two-branch amplitudes for every
\(1\le u\le16\), \(1\le N\le64\). It also counts exceptional rotations.
These are numerical checks of the induction, not a substitute for it.

The result supplies a constant-cardinality algebraic gate menu for the
[oracle scheduler](ORACLE_SCHEDULER.md). It does not provide initialized
history-Hamiltonian penalties or an infrared stability theorem.
