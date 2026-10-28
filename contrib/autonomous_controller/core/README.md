# An explicit initialized autonomous qudit core

**Result.** [`core.py`](core.py) specifies a single two-site Hermitian matrix
\(h_c\) and two one-site boundary matrices \(a_c,b_c\), independent of chain
length. For every \(n\ge4\), their open-chain Hamiltonian

\[
H_c(n)=a_c^{(1)}+\sum_{j=1}^{n-1}(h_c)_{j,j+1}+b_c^{(n)}
\]

is positive, has exactly \(2^n\) zero-energy states, and has a strictly positive
finite-size gap. Its zero states are initialized computation histories of the
[integrated controller](../bridge/INTEGRATED_READOUT.md): it creates its binary
program from blank working memory, applies coherent quantum gates, exposes
selected records at the endpoints during safe READY intervals, restores them,
and uncomputes. At sufficiently long lengths a whole such cycle lies inside
the history.

This is a demonstration of **autonomous initialization and execution**, not
the high-channel sparse-fermion code. No central charge or c–d violation is
claimed. The \(2^n\)-dimensional ground space alone does not imply \(n\)
accessible gapless channels after coupling cells.

## 1. An actual local Hilbert space and matrices

One site has a 15-state clock and a work register. The work register consists
of five static tape marks, four dynamic classical bits, either no head or
one of 40 controls in a fresh/updated phase, and one neutral quantum qubit:

\[
\dim W=5\cdot2^4\cdot(1+2\cdot40)\cdot2=12,960,
\qquad \boxed{d_c=15\dim W=194,400.}
\]

The basis order is public in `encode`/`decode`. The functions
`bond_column(left,right)`, `left_column(index)` and `right_column(index)`
return every nonzero coefficient of the requested matrix column. They are
defined on every basis input, including malformed tapes and multiple heads.
All entries are complex integer multiples of \(1/4\), exactly representable
by the binary arithmetic used in the implementation. The repeated bond has
37,791,360,000 rows and the same number of columns; its sparse column rule is
the explicit specification of that finite matrix.

The work unitaries are the literal full-space columns in
[`local_gates.py`](../bridge/local_gates.py):

- \(G\) moves a properly oriented fresh head one site and applies its named
  `CNOT`, `CV`, `SWAP`, or identity gate. The reverse block carries the
  adjoint. Every other block is identity. Thus \(G=G^\dagger=G^{-1}\).
- \(F\) applies the completed classical collision permutation to an updated
  head and its inverse to a fresh head, toggling the phase. It is an on-site
  involution and acts identically on the quantum qubit.
- \(B\), the seed, swaps no head with a fresh `F:ADVANCE` head whenever the
  site's static mark is LEFT. It preserves the dynamic bits and quantum
  record, and is identity on the remaining states.

Every static tape mark is preserved by all three gates, on all inputs.

## 2. Local positive terms

The clock rules and allowed pairs are those of the proved
[15-state full-sweep clock](../../../codex/work_C8/sweeping_clock.md).
The source repeats their small combinatorial definition so matrix queries
need only the standard Python library. Each directed clock edge \(s\to t\)
is given the positive propagation block

\[
\frac12\begin{pmatrix}I&-U^\dagger\\-U&I\end{pmatrix},
\]

with the following work unitary:

| Clock edge | Work gate |
| --- | --- |
| Outward sweep | \(G\) on the traversed bond |
| Return sweep | \(F\) on the departed right site |
| Left sweep finish | \(F\) on the left endpoint |
| Initialization `I Z → L I` | \(B\) on the departed left site |
| Carry, turn, return, sweep entry, right bounce, activation | Identity |

The clock's on-site turn projector `C0 ↔ R1` is split equally between
adjacent bonds. Each endpoint term supplies the missing half. The same split
is used for the on-site initialization pin described below. This gives one
length-independent bond, rather than a list of position-dependent gates.

Clock-pair and clock-boundary defects have unit penalty. Static tape-pair
defects also have unit penalty, with allowed pairs

```
LEFT OUT_L; OUT_L BULK; OUT_L OUT_R;
BULK BULK; BULK OUT_R; OUT_R RIGHT.
```

The endpoint terms require static LEFT and RIGHT respectively. For each
\(n\ge4\), these constraints select the unique static word

```
LEFT OUT_L BULK^(n-4) OUT_R RIGHT.
```

They do not supply a program string: all dynamic program bits are initially
zero. Their subsequent values are produced by the controller.

Finally, at clock letter `I`, the on-site pin requires **no work head and all
four dynamic bits zero**. The quantum qubit remains unconstrained. The seed
acts only after the leftmost `I` pin has been encountered. During later
initialization moves it is identity, because their departed static mark is
not LEFT. This explicitly resolves the seed/initialization ordering issue.

## 3. Complete ground-space proof

Every summand is positive. Therefore a zero vector must satisfy all clock
and static constraints. The static word is unique. The complete legal clock
graph is one path with

\[
C_n=(2n+4)2^n+n-4
\]

vertices, making exactly \(2^n\) full spatial sweeps. This is the clock's
full legal sector, not a restriction to a preferred head configuration.

Write \(V_t\) for the product of work unitaries along the path from its
initial word \(I Z^{n-1}\) to time \(t\). All work matrices are defined
unitaries on the full work space, so propagation is exactly gauge-equivalent
to one half the path Laplacian tensored with that space. Its zero vectors
are consequently precisely

\[
|\Psi(\psi)\rangle=C_n^{-1/2}
\sum_{t=0}^{C_n-1}|t\rangle V_t|\psi\rangle.
\]

During the initialization prefix the `I` head visits every site. Its pins
require no work heads and zero dynamic bits at the initial time. The only
earlier nonidentity work gate is the seed at the left endpoint; it has
disjoint support from every later pin. Thus, after pulling the pins back to
\(t=0\), their joint kernel is exactly one fixed classical input tensored
with \((\mathbb C^2)^{\otimes n}\). The quantum inputs are arbitrary,
including entangled superpositions. The zero-space dimension is exactly
\(2^n\).

In particular, additional head configurations or periodic orbits of the
work machine cannot create extra zero histories. The clock has one legal
path, and every propagated input must satisfy all its initialization pins.

For an explicit gap, set

\[
\gamma_n=1-\cos(\pi/C_n).
\]

The summed initialization pin operator has norm at most one on the legal
clock space, since there is at most one `I` at each time. Its compression
to the constant-clock propagation kernel is at least \(1/C_n\) on every
rejected input. The positive-core pin lemma in
[COMPILER_INTERFACE.md](../../torus_bottleneck/COMPILER_INTERFACE.md)
therefore yields

\[
\boxed{\operatorname{gap}H_c(n)\ge
\frac{\gamma_n}{1+C_n+C_n\gamma_n}>0.}
\]

Illegal clock or static sectors have energy at least one: their defect
projectors are conserved by propagation, and all remaining terms are
positive. The displayed smaller bound therefore applies to the entire
Hamiltonian. It is a finite-size bound; no size-independent gap is asserted.

## 4. What the actual history executes

One outward/return sweep executes the ordered product of all \(G\) bonds,
then all \(F\) sites from right to left. On the initialized single-head
orbit this is exactly one step of the integrated controller. The explicit
G/F sequence also specifies every intermediate clock vertex, where readout
safety is checked separately.

The integrated controller's cycle takes

\[
T_{\rm cycle}(n)=2\bigl(n^2+8n-11-2\operatorname{popcount}(n-1)\bigr)+6n-8
\]

steps. The inequality \(2^n\ge T_{\rm cycle}(n)\) holds for all \(n\ge9\),
so those histories include a complete blank-to-program-to-readout-to-cleanup
cycle. Shorter histories remain well-defined initialized cores; they can
stop partway through a computation. In particular the \(n=4\) test stops
before program execution and says nothing about readout.

Independent complete-clock tests in the
[audit results](../audit/sweep_audit_results.json) find no READY intervals
at \(n=4,5,6\), positive endpoint intervals at \(n=7,8\), and complete
cycles at \(n\ge9\) in the tested range. At \(n=9\), the actual uniform
Hamiltonian history has weights \(658/11269\) and \(656/11269\) for the two
endpoint READY operators. These differ because the clock adds idle times;
they are not the weights of a uniformly sampled machine-gate list.

## 5. Matrix and history checks

Run `python3 contrib/autonomous_controller/core/checks.py`.
The saved [results](results.json) include:

- every column of both boundary matrices checked for Hermiticity;
- all 12,960 columns of \(F\) and \(B\) checked for involution and static
  mark conservation;
- sampled bulk Hermiticity covering every clock-pair label and every
  controller label, alongside arbitrary malformed work inputs;
- exhaustive clock-sector classification through \(n=4\), including all
  50,625 clock words at that length;
- all 16 accepted \(n=4\) basis histories checked against the **complete
  global Hamiltonian**, with exactly zero residual and identity Gram matrix
  after normalization;
- an \(n=6\) superposed input whose history includes coherent gate branching,
  again with exactly zero complete-Hamiltonian residual;
- an \(n=9\) input \(|111101111\rangle\) whose full history includes genuine
  coherent gate branching,
  a complete cleanup cycle, and both positive readout intervals, checked
  against the complete Hamiltonian with zero residual.

The full \(d_c^2\)-column bulk matrix is not enumerated. Its full-space
Hermiticity and positivity follow from the explicit involution blocks and
positive-term decomposition; the sampled checks test the implementation.
The ground-space completeness and gap above are proofs based on the complete
clock/initialization argument, rather than claims from a small diagonalization.

The remaining high-channel work is to replace this binary-digit demonstration
program by the actual sparse-circuit generator, compile its word-level
instructions to local finite transitions, and establish the required
fermionic boundary compressions and full resource count. This core closes a
concrete initialization/composition example, not that larger construction.
Its free input qubits and neutral CNOT/CV gates do not supply the indexed
\(N\)-flavor fermionic endpoint maps required by the candidate construction.
There is no fixed-\(d\) large-central-charge conclusion from this code.

**Provenance.** Implemented with GPT-6-family assistance on 2026-10-28 and
independently checked by another AI worker. No human review is claimed.
