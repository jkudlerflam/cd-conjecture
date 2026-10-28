# Explicit resource bounds for this circuit and preloaded dispatcher

**Scope.** These are conservative numerical bounds obtained from the loops in
[`sparse_circuit.py`](sparse_circuit.py), not fits to test instances. They count
the two-port circuit, its exact nearest-neighbor lowering, and the literal
cursor program consumed by the existing dispatcher. They **exclude** autonomous
program generation, verification, the history clock, transport to actual cell
endpoints, all other controller work, and the mobile wrapper. In particular,
the tape-length bound below is not a sufficient bound on the final core length.

Let the input be any nonempty finite simple bipartite graph, with `u,v>=1`,
`q=u+v`, and `1<=N<=uv`. Define

\[
L=N+4q,\qquad
w=\max\{1,\lceil\log_2q\rceil,
              \lceil\log_2\max(u,N)\rceil\},\qquad A=L^2w.
\]

Then `L>=9`, `1<=w<=L`, and the code's species and branch bit widths are at
most `w`. The local circuit alphabet obeys

\[
|\Sigma|\le q+3+4(q+1)w.
\]

The following bounds hold for the actual generated two-port circuit:

| Resource | Upper bound |
| --- | ---: |
| Neutral bit records `W` | `32 A` |
| Abstract gate instructions `G` | `128 A` |
| Charged plus neutral data sites `D=L+W` | `33 A` |
| Actual two-site gates after SWAP routing `J` | `70,000 A^2` |
| Literal cursor-program symbols `P` | `2,310,000 A^3` |
| Preloaded dispatcher tape sites `n=P+2+D` | `2,310,035 A^3` |
| Dispatcher forward/inverse cycle head steps | `30,000,000,000,000 A^6` |

The constants are intentionally loose. Circuit sites and program symbols are
counts of positions, not local Hilbert-space dimensions.

## Boolean building blocks

A new XOR output consumes one fresh bit and one CNOT per input. NOT-on-copy
consumes one bit and two gates. AND consumes one bit and one Toffoli (or a
CNOT in the identical-control case). OR consumes two bits and four gates.
Comparison of an `r>=1` bit word with a fixed integer consumes at most `2r`
new bits and `3r` gates. A zero-bit comparison consumes one bit and one gate.

These bounds count the code's actual new-output allocation policy. Reversing
a recorded Boolean block adds gates but allocates no further bits.

Each port processes `S=N+2q<=L` charged records. Its neutral branch register
uses at most `w+1` bits. The preparation loops contain at most

\[
(u+N)w\le Lw
\]

prefix predicates: at any depth the number of nonempty prefixes for an
interval of size `M` is at most `M`. Each predicate costs at most `2w+3`
new bits and, including its reversal and controlled rotation, `6w+9` gates.
There is one additional Hadamard per stage.

## The target, routing, and inverse loops

For one stage, use the following budgets. The target and inverse gate counts
are **one-way computation counts**; their explicit reversal is added below.

| Source routine | New bits | Gates |
| --- | ---: | ---: |
| `target_network` | `17 L w` | `7 L^2 w` |
| `routing_network` | `5 L^2 w` | `8 L^2 w` |
| `inverse_network` | `4 L^2 w` | `6 L^2 w` |

Here is the loop accounting behind these three rows.

For the target, the occupation comparisons cost at most `2Nw` bits and
`3Nw` gates. Upper/edge branch selectors cost at most `L(2w+2)` bits and
`L(3w+3)` gates. Creation/annihilation selectors cost `3N` bits and `4N`
gates. The port word uses at most `w` bits and `Lw` gates. The `2q-1` flag
positions use at most `Lw` bits and `2L^2w` gates, since their disjoint case
list has size `u+2N<=2L`. Residual occupations use at most `N(w+1)` bits and
`N(2w+2)` gates. Even for the more expensive right-port parity calculation,
the sign uses at most `3N+1` bits and `2L^2+4L` gates. Thus totals are at most

\[
6Lw+9L+w+1\le17Lw,
\]
\[
2L^2w+2L^2+9Lw+13L\le7L^2w.
\]

The routing loops have `S(S-1)/2<=L^2/2` iterations. An iteration performs a
word comparison, OR update, and a coherent word swap. These consume at most
`7w+2` bits and `12w+4` Boolean gates. At most `L` extra zero `found` bits
are allocated. Hence

\[
W_{\rm route}\le L+\tfrac12L^2(7w+2)\le5L^2w,
\qquad
G_{\rm route}\le\tfrac12L^2(12w+4)\le8L^2w.
\]

The inverse compares the flag word, of at most `Lw` bits, against `u+N<=L`
canonical flags. It also decodes the port, residual occupations, branch,
and fixed reference words. Counting those loops gives

\[
W_{\rm inv}\le2L^2w+8Lw+3L+w+1\le4L^2w,
\]
\[
G_{\rm inv}\le3L^2w+13Lw+6L\le6L^2w.
\]

These inequalities use only `L>=9` and `w>=1`; there is no enumeration of
`2^N` logical occupations in the generator.

## Combining both stages

Besides those three networks, each stage allocates the branch preparation,
one queried source word, and one queried output word. Counting all of them
without reusing any cleared bit gives

\[
W\le2[9A+22Lw+2Lw^2+w+1]\le32A.
\]

The target and routing computations are both reversed once, as is the
inverse decoder. The charged routing itself uses at most `L^2/2` abstract
controlled swaps per stage, followed by one phase gate. There are three
word-query passes and at most `Lw+w+1` cleanup CNOTs. Including preparation,
the full two-stage abstract instruction count obeys

\[
G\le84A+L^2+12Lw^2+26Lw+2w+6\le128A.
\]

## Nearest-neighbor lowering

A two-register operation is brought together by SWAPs and restored, costing
at most `2D` local instructions. The five-pair Toffoli decomposition costs
at most `10D`. A species query costs at most `2D+6w`; an indexed rotation
costs at most `2D+12w`. Gathering the controlled data-swap triple, applying
`LOAD,EXCHANGE,LOAD`, and restoring the route costs at most `6D+9`.
One-site operations cost at most three instructions, including the last-site
boundary case. Since `D>=L>=9` and `D>=w`, **every abstract instruction costs
at most `16D` two-site gates**.

There are two `INIT` passes over `W` neutral sites, costing at most `6W`
additional gates. Therefore

\[
J\le16DG+6W\le33A(2048A+6)<70{,}000A^2.
\]

This bound counts the actual routes in the source rather than assuming an
abstract nonlocal gate is already nearest-neighbor.

## The literal dispatcher program

The dispatcher lowerer starts its data cursor at address zero. Consecutive
local-gate addresses lie in `[0,D-2]`, so each gate emits at most `D-2` cursor
moves and one oracle symbol. Consequently

\[
P\le J(D-1)\le DJ\le2{,}310{,}000A^3.
\]

The program is executed forward and then backward by the same dispatcher;
a second copy of the program is unnecessary. The existing two endpoint
symbols give

\[
n=P+2+D\le2{,}310{,}035A^3.
\]

The dispatcher's already specified traversal bound is
`2(P+1)(2n+6)`, which is below `3*10^13 A^6` under these estimates.
No autonomous generator's control states or workspace are inferred here.

## A numerical projective-plane corollary

For any prime `p>=2`,

\[
L=(p+9)(p^2+p+1)<10p^3.
\]

Here `N>=q`, so `w=ceil(log2 N)`. Since `N<=3p^3<=2^(3p)` for integers
`p>=2`, we have `w<=3p`. The latter exponential inequality follows at
`p=2` and persists because multiplying `p` by `(p+1)/p` increases `p^3`
by at most `(3/2)^3<8`. It follows that

\[
A\le300p^7,\qquad
n\le2{,}310{,}035\,(300p^7)^3
 <6.3\times10^{13}p^{21}<2^{50}p^{21}.
\]

Thus a tape of `2^50 p^21` sites is a concrete sufficient capacity for
**this preloaded circuit/dispatcher layer**. This is not a selected core
length: a length-to-program generator, its workspace, verification,
endpoint transport, and the complete clock composition remain to be counted.
It supplies neither the final local dimension nor a violating prime.

## Reproduction

Run

```
python3 contrib/controller_instantiation/circuit/resource_checks.py
```

The script writes [`resource_results.json`](resource_results.json). It checks
the bounds against generated one-edge, isolated-vertex, complete-bipartite,
`p=2`, and `p=3` circuits. Exact local counts use the analytically counted
route lengths; only the tiny stream is expanded. The actual tiny dispatcher
program length is reconstructed independently from successive addresses.
Parameter-only checks through `p=1000` do not construct those larger graphs
or circuits. These checks verify the implementation against the stated
inequalities; the source-loop argument above establishes their general form.

**Provenance.** Derived with AI assistance by Codex on 2026-10-28. No human
review is claimed.
