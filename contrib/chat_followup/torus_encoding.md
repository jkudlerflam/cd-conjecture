# A species-torus refinement of the first-order encoding proposal

Imported 2026-10-27 from the construction follow-up in the user's earlier
chat **“Prove c d conjecture.”** This records material absent from the
repository's existing joint-compression and infrared-control reports.
The retrieved source ends during its final RG discussion; this note does
not reconstruct the missing ending. The source assistant's model was not
identified in the retrieved transcript. Codex prepared this import and
rederived the elementary algebra below; there has been no independent
human review or numerical validation of this proposal.

**Status.** [proved] The finite-dimensional equivariant compression and
the integer-kernel selection rule below. [heuristic / open] Their
realization by the claimed small-alphabet nearest-neighbour Hamiltonian,
the full infrared conclusion, and the proposed violation of c–d. In
particular, the source's statement that only a constructive RG extension
remains is too strong: its interpreter is an outline, not a complete
microscopic construction.

The proposed family has a data register of dimension $b=2r$, $N=r^2$
logical complex fermions per cell, and an intended full physical onsite
dimension $d=O(r\log r)$. If the remaining steps worked, its gapless
logical theory would have total Virasoro central charge $c=N$. This
would challenge $c\le d-1$ for sufficiently large $r$. Because $d$
grows with $r$, this proposal **would not establish unbounded c at fixed
d**, which is the stronger objective of route A in [STATUS](../../STATUS.md).

## 1. Exact compression with an explicit abelian symmetry

Let the physical data register be

\[
B=\operatorname{span}\{|0,\beta\rangle,|1,\alpha\rangle:
1\le\alpha,\beta\le r\},\qquad
F_{\alpha\beta}=|1,\alpha\rangle\langle0,\beta|.
\]

Let $\mathcal K=\mathcal F_N$ be the Fock space of $N=r^2$ logical
fermions $f_a$, with $a=(\alpha,\beta)$, number $Q$, and parity

\[
P=(-1)^Q.
\]

Define Hermiticity-preserving maps $\Gamma_L,\Gamma_R:M_b\to\mathcal B(\mathcal K)$ on matrix units by

\[
\Gamma_L(F_a)=f_a^\dagger,\qquad
\Gamma_R(F_a)=f_a^\dagger P.
\]

Assign the adjoint matrix units the adjoint images, and assign all other
matrix units zero. In particular $\Gamma_L(I)=\Gamma_R(I)=0$. With

\[
\tau(A)=\operatorname{Tr}(A)/b,
\]

extend linearly from simple tensors the map

\[
\Phi(A\otimes B')=\tau(A)\tau(B')I+
\epsilon\{\tau(B')\Gamma_L(A)+\tau(A)\Gamma_R(B')\}.
\tag{1}
\]

This map is unital. In the unnormalized Choi convention, its first term
has Choi matrix $I/b^2$. Each charged matrix-unit/adjoint pair contributes
a Hermitian off-diagonal Choi block of norm one, so

\[
\|C_{\Gamma_L}\|,\|C_{\Gamma_R}\|\le N,
\qquad
\left\|C_\Phi-I/b^2\right\|\le 2N\epsilon/b.
\]

Consequently $0<\epsilon\le1/(4bN)$ ensures

\[
C_\Phi\ge I/(2b^2)>0.
\]

The Choi criterion and Stinespring dilation give an isometry

\[
V:\mathcal K\longrightarrow B_L\otimes B_R\otimes E,
\qquad \dim E\le b^2 2^N,
\]

with exact compressions

\[
V^\dagger F_a^{(L)}V=\epsilon f_a^\dagger,
\qquad
V^\dagger F_a^{(R)}V=\epsilon f_a^\dagger P.
\tag{2}
\]

The construction is real in occupation and species bases. Since
compression need not preserve products, the commuting physical port
algebras can have noncommuting logical images. This is a refinement of
the existing [C8 joint-compression result](../../codex/work_C8/RESULT.md),
not a new general mechanism or a local-parent-Hamiltonian theorem.

The refinement is its full species symmetry. Give the $2r$ physical
species independent phase rotations, with charge basis vectors

\[
e_{0,\beta},e_{1,\alpha},\qquad
w_{\alpha\beta}=e_{1,\alpha}-e_{0,\beta}.
\tag{3}
\]

Assign the logical creation operator $f_a^\dagger$ charge $w_a$.
Equation (1) is equivariant under $U(1)^{2r}$; the simultaneous overall
phase has trivial action on operators. A covariant Stinespring
environment can be taken inside

\[
\mathcal K\otimes\overline B\otimes\overline B.
\]

There is an explicit embedding of this representation into ordinary
data sites. Encode logical occupation $n_{\alpha\beta}=0,1$ by

\[
|0\rangle\mapsto|0,\beta\rangle,
\qquad |1\rangle\mapsto|1,\alpha\rangle.
\]

The $N$ sites contribute a reference charge

\[
R=r\sum_\beta e_{0,\beta}.
\]

Write $S=\sum_s e_s$. A dual species vector of weight $-e_s$ can be
represented, after the common shift $S$, by a fixed-order string of
$b-1$ sites containing every species except $s$. Two such strings
realize the two dual factors, up to the common shift $2S$. Thus the
environment fits into $N+2(b-1)$ data sites and the two ports bring the
total to $N+2b$. Every output has charge

\[
R+2S+\sum_a n_aw_a,
\]

matching an input of the $N$ logical occupation sites plus two reference
strings containing every species. The reference upper-species count is
$2r$, which is even, so physical upper-species parity equals logical
fermion parity. This verifies representation compatibility; it does not
implement the isometry by a repeated nearest-neighbour interaction.

On a chain of encoded blocks, define

\[
c_{j,a}=\left(\prod_{i<j}P_i\right)f_{j,a}.
\]

Then

\[
f_{j,a}^\dagger P_jf_{j+1,a}=c_{j,a}^\dagger c_{j+1,a}.
\]

The exchange $F_a^{(R)}\otimes (F_a^{(L)})^\dagger+\mathrm{h.c.}$
therefore compresses to $\epsilon^2$ times the corresponding logical
hopping. The physical exchange preserves every species count. On a ring
the usual parity-dependent fermionic boundary condition must still be
tracked.

## 2. The incidence-kernel selection rule

Define the integer charge map

\[
W:\mathbb Z^{r\times r}\longrightarrow\mathbb Z^{2r},
\qquad Wp=\sum_{\alpha,\beta}p_{\alpha\beta}w_{\alpha\beta}.
\]

It is the oriented incidence map of the complete bipartite graph:

\[
Wp=0\iff
\sum_\beta p_{\alpha\beta}=0\quad(\text{every row}),
\qquad
\sum_\alpha p_{\alpha\beta}=0\quad(\text{every column}).
\tag{4}
\]

For $r\ge2$, its rank is $2r-1$ and its kernel has rank $(r-1)^2$.
Every nonzero kernel vector has at least four nonzero integer entries:
a nonzero entry forces a second entry in its row and a second in its
column, and those entries cannot cancel all remaining row/column sums
without another entry. The rectangle

\[
\begin{pmatrix}1&-1\\-1&1\end{pmatrix}
\]

embedded in any two rows and columns achieves the minimum. Hence

\[
\min_{0\ne p\in\ker W}\|p\|_1=4,
\qquad
\min_{0\ne p\in\ker W}\|p\|_2^2=4.
\tag{5}
\]

For $r=1$ the kernel is zero. This symmetry is weaker than exact
$U(1)^N$ conservation of every logical flavour: rectangle processes
are allowed. It is nevertheless sufficient to forbid the dangerous
four-field Cooper transfer described below.

Take the free logical target

\[
H_f=\sum_{j,a}\{-t_a(c_{j,a}^\dagger c_{j+1,a}+\mathrm{h.c.})
+n_{j,a}\},
\qquad t_a=(2\cos k_a)^{-1},
\tag{6}
\]

with $k_a\in(\pi/4,\pi/3)$. Each band has simple Fermi points
$\pm k_a$ and positive speed $2t_a\sin k_a$. One common number term
can be obtained from total physical upper-species charge, whose code
restriction is $Q$ plus a constant.

For example, choose distinct primes $\mathfrak p_a$ and a sufficiently
small positive rational $\tau$, and set

\[
k_a=\pi\left(\frac14+\tau\sqrt{\mathfrak p_a}\right).
\]

Rational independence of $1,\sqrt{\mathfrak p_1},\ldots,\sqrt{\mathfrak p_N}$ gives $q\cdot k\notin2\pi\mathbb Z$ for every
nonzero integer $q$. Stability under perturbation only requires a
positive minimum mismatch for the **finite** set of momentum relations
that could be relevant or marginal near the free point; it does not
require all-order Diophantine stability.

For a local fermion monomial, let $u_a,v_a$ be the net numbers of right-
and left-moving creation operators, respectively. Its charge and
oscillating momentum are

\[
W(u+v),\qquad (u-v)\cdot k.
\]

For the potentially relevant/marginal monomials with at most four
fermion fields, nonresonance forces $u=v=p$, and neutrality forces
$Wp=0$. A nonzero $p$ then requires at least $2\|p\|_1\ge8$ fields,
contradicting the four-field restriction. In particular, for $a\ne b$,

\[
\psi^\dagger_{a,R}\psi^\dagger_{a,L}\psi_{b,L}\psi_{b,R}
\]

has nonzero charge $2(w_a-w_b)$ and is forbidden exactly. The remaining
quadratic and quartic operators, after normal ordering, are density
shifts, kinetic terms, and forward density/current interactions. This
selection rule removes a marginal pair-transfer sign problem present
when only total particle number is protected.

As a bosonic check, use the convention

\[
\Delta(m,n)=\tfrac14(m^T K m+n^T K^{-1}n),\qquad K>0.
\]

A nonoscillating channel-changing vertex has $m=0,n=2p,Wp=0$, so

\[
\Delta(0,2p)=p^T K^{-1}p\ge4\lambda_{\min}(K^{-1}).
\tag{7}
\]

Its minimum dimension is four at $K=I$; it stays strictly above two
whenever $\lambda_{\max}(K)<2$. Near the free theory, higher-dimension
operators stay irrelevant by continuity and a uniform positive lower
bound on the quadratic form. The potentially dangerous nonzero-momentum
vertices form a finite set there and must retain their mismatch. These
are statements about perturbations of a proposed Gaussian theory,
not a proof that a microscopic model flows to it.

Only $2r-1$ independent physical charges act on the $r^2$ flavours.
Thus the proposal does not meet the existing conditional bound's
hypothesis that microscopic commuting charges resolve *all* Gaussian
currents. Its additional low-energy flavour currents would be emergent.

## 3. The proposed interpreter and its unproved dimension count

The source takes $b=2^k$ (so $r$ runs through powers of two), appends
a fixed number of neutral work qubits per site, and proposes a reversible
interpreter whose charged primitives are species-preserving SWAPs and
the $k$ exact bit queries

\[
R_t=\sum_{s=0}^{b-1}|s\rangle\langle s|\otimes
X^{\operatorname{bit}_t(s)},\qquad 0\le t<k.
\tag{8}
\]

These primitives are real and conserve every species count. Neutral
controls can therefore preserve the torus exactly even when the target
isometry is only approximated. The proposed circuit strategy copies
species labels into neutral workspace, reversibly sorts the species
into a canonical order, performs a real transformation on the neutral
encoding within each count sector, and reverses the bookkeeping.
Neutral workspace can occupy more sites; it need not all be located on
the site holding a species register.

The claimed implementation uses the block length to encode the program,
initialization checks to leave exactly $2^N$ free logical input states,
and endpoint flags active only after both outputs have been prepared.
For an actually realized history isometry $W_h$, such a flag would give

\[
W_h^\dagger(\Pi_{\mathrm{flag}}O_{\mathrm{end}})W_h
=p_{\mathrm{end}}(V')^\dagger O_{\mathrm{end}}V',
\]

with $p_{\mathrm{end}}>0$ and an approximation $V'$ to $V$. The
first-order logical hopping scale would then be

\[
\eta=\lambda p_Lp_R\epsilon^2.
\]

The proposed alphabet estimate is

\[
d\le Cb(1+\log_2 b),\qquad C\text{ independent of }b
\text{ and program length}.
\tag{9}
\]

The gate list alone does **not** prove (9). The source supplies no complete
repeated two-site interaction, controller transition table, verified
alphabet count, initialization and illegal-sector analysis, or complete
two-endpoint readout schedule. Existing computation-history frameworks
were invoked in the source, but no checked specialization establishing
all of these properties was supplied. This remains the same category
of missing compiler/local-realization step identified in
[C8's progress report](../../codex/work_C8/PROGRESS.md) and
[STATUS](../../STATUS.md), despite the plausible logarithmic query index.

If (9) were eventually established and all $N$ modes survived, then

\[
\frac{c}{d}\ge\frac{b}{4C(1+\log_2 b)}\longrightarrow\infty.
\]

No value of $C$, threshold $r$, or explicit physical bond realizing
this ratio has been certified here.

## 4. What would still be required for a counterexample

1. **Microscopic realization.** Construct the same Hermitian nearest-
   neighbour bond on every site, with the full onsite space counted in
   (9), exact species symmetry, the intended code and flags, a positive
   cell gap, and no lower-energy unwanted sectors.
2. **Thermodynamic control.** Control errors from approximate encoding and
   excursions outside the code in a local norm, uniformly in the number
   of cells, at a fixed positive intercell coupling. The source proposes
   a massive-auxiliary-fermion enlargement for this step; the
   finite-dimensional compression itself supplies no such bound.
3. **Infrared identification.** Show that the resulting general weak
   interactions, including higher-degree or time-retarded effective
   kernels, remain in an $N$-boson phase, and connect that conclusion
   to the conjecture's unitary conformal spectrum. The selection rule
   helps with this step but is not the needed RG theorem.
4. **Ring lengths and translation.** The source proposes frozen cells
   and a chosen period-$p$ phase, or rings $L=pM$. That qualification
   does not resolve the repository's obstruction from frozen-marker
   degeneracy on incompatible rings. A mobile-cell implementation would
   require additional work preserving the species symmetry and the
   code during marker motion. No such implementation is provided here.

For a proposed family member, the intended order is to fix $r$, then
the finite code accuracy, program and cell length, then a sufficiently
small but positive coupling, and only afterwards take $L\to\infty$.
Allowing microscopic parameters to vary with $L$ would not establish
the claimed counterexample.

There are no numerics to reproduce in this note. Equations (1)–(8) give
the elementary algebraic checks in full. The new reusable point is the
species-charge assignment and its eight-field threshold; the global
status of the c–d conjecture is unchanged.
