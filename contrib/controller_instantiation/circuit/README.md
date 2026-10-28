# An executable sparse-port circuit, including cleanup

**Status, 2026-10-28.** This directory replaces the reversible-transducer
pseudocode in [the circuit construction](../../exact_route/CIRCUIT.md) by an
executable generator of an explicit finite gate sequence. Its output is still
an **externally scheduled circuit**. It does not supply the initialized,
length-reading autonomous controller or establish a counterexample to the
central-charge conjecture.

The generator uses only Python's standard library. Every schedule and gate
address depends on the input graph, never on an observed quantum value. The
classical-looking intermediate registers are coherent neutral quantum bits.

## What is now instantiated

[`sparse_circuit.py`](sparse_circuit.py) constructs all of the following:

1. The projective-plane incidence graph from a prime `p`, with explicit point,
   line, and edge orderings.
2. Exact upper/edge branch preparation, using a Hadamard, Boolean predicates,
   and indexed rotations whose squared entries are rational.
3. A reversible Boolean network computing the target charged word and the
   exact Fock signs. This includes the parity dressing in the right port.
4. A fixed bubble-routing schedule whose swap controls are coherently
   computed from the source and target words.
5. Physical, charge-preserving routing of the actual species records.
6. Reversal of every routing/target workspace computation.
7. A separately generated decoder of the output word, which erases the saved
   source word and branch label, followed by reversal of the decoder.
8. Both successive port encodings, with disjoint reference words and a shared
   residual Fock register.
9. A gate-by-gate lowering to a line using a fixed menu of **17 two-site
   opcode families**, with inverses included.

No unspecified Boolean function is called as a quantum oracle. The generated
Boolean portion is a list of NOT, CNOT, and Toffoli gates. Charged queries and
controlled transpositions are then lowered to the specific two-site gates
below. Temporary Boolean wires are retained and uncomputed explicitly.

## Initialized input and final output

For a graph with `q = u + v` species and `N` edges, the charged input is

```
N logical occupation records | two copies of each species | two copies of each species
```

A logical edge `(alpha,beta)` stores species `alpha` when occupied and species
`u+beta` when empty. Each stage maps its `N+2q` charged records to

```
port | canonical flag word of length 2q-1 | residual N occupation records
```

The two-stage output is

```
left port and flag | right port and flag | residual N occupation records
```

All neutral auxiliary sites begin in the **unique BLANK symbol**. Before the
circuit, the `INIT` permutation exchanges BLANK with the ordinary neutral bit
zero on each designated auxiliary site. Every such bit ends at zero; the final
`INIT` operations return all of them to BLANK. This input contract avoids
asking a TYPE-only initializer to distinguish neutral bit zero from bit one.
The trusted verifier scratch and its separately pinned helper belong to the
external controller, not this circuit's auxiliary records.

For each valid basis branch the target source multiset agrees exactly. The
routing algorithm processes target positions from left to right. For each
position it scans from the right, setting a coherent `found` bit when it sees
that target species and conditionally moving the matching symbol left. This
uses a predetermined sequence of at most `L(L-1)/2` adjacent charged swaps per
stage. It need not make classical decisions during execution.

The original branch and input are erased using the port and canonical flag:
upper ports identify an upper branch; lower ports with `g_alpha` identify a
creation branch; lower ports with `h_edge` identify an annihilation branch.
The residual bit recovers the original occupation. Inputs outside the promised
format need not decode correctly, but the circuit remains unitary there
because every elementary gate is a full-space unitary.

## Exact local alphabet and gates

Let

\[
\ell=\max\{1,\lceil\log_2\max(u,N)\rceil\},\qquad
m=\max\{1,\lceil\log_2 q\rceil,2\ell\}.
\]

The implementation uses

\[
\Sigma=\{D_s\}_{s<q}\sqcup\{\mathrm{blank},0,1\}
\sqcup\{H_{j,c}\}_{j<m,c<2}
\sqcup\{T_{j,c,s}\}_{j<m,c<2,s<q},
\]

so

\[
|\Sigma|=q+3+2(q+1)m=O(q\log q).
\]

Only `D_s` and `T_(j,c,s)` have species charge, namely `e_s`.

The 17 families are

```
SWAP X H Z CNOT CV CV_DAG TYPE SHIFT SHIFT_DAG
QUERY LOAD EXCHANGE ROT ROT_DAG CH INIT
```

`X,H,Z,CNOT,CH` have their usual actions within the ordinary neutral bit
subspaces and are the identity elsewhere. `CV` controls on an ordinary bit
one and applies

\[
V=\frac12\begin{pmatrix}1+i&1-i\\1-i&1+i\end{pmatrix},\qquad V^2=X,
\]

to a second ordinary neutral bit; `CV_DAG` is its adjoint. The exact five-gate
sequence `CV(b,t), CNOT(a,b), CV_DAG(b,t), CNOT(a,b), CV(a,t)` implements a
Toffoli with controls `a,b` and target `t`.

`SWAP` exchanges complete neighboring site states. `TYPE` exchanges ordinary
bit `c` with `H_(0,c)`. `SHIFT` increments the index on `H` and `T` modulo `m`.
`QUERY`, `LOAD`, and `EXCHANGE` act by

\[
\begin{aligned}
(D_s,H_{j,c})&\mapsto(D_s,H_{j,c\oplus\operatorname{bit}_j(s)}),\\
(D_s,H_{j,c})&\longleftrightarrow(\mathrm{blank},T_{j,c,s}),\\
(T_{j,1,s},D_t)&\longleftrightarrow(T_{j,1,t},D_s).
\end{aligned}
\]

The query is the identity for index values beyond the species bit width.
`EXCHANGE` is the identity when the loaded control bit is zero. All unspecified
blocks are the identity. A controlled data transposition is the three-gate
`LOAD, EXCHANGE, LOAD` sequence after gathering `(data,head,data)` together;
all gathering swaps are undone afterwards.

`ROT` applies the indexed real rotation `R_j` to the head's bit conditional on
an adjacent ordinary bit one. The exact rational pairs specifying
`R_j=[[sqrt(A),-sqrt(B)],[sqrt(B),sqrt(A)]]` are emitted in the manifest.
`ROT_DAG` is its inverse. `INIT` exchanges BLANK and ordinary bit zero.
One-site families act on the first factor of a two-site gate, with identity
on the second; the localizer handles the last site by a swap and restoration.

Each gate conserves every total species charge on the entire local Hilbert
space, not merely on the accepted data domain. Gate-family cardinality is
constant; the indexed matrix entries and the species alphabet depend on `q`.

## Concrete size of the first projective-plane instance

The generated [`p2_manifest.json`](p2_manifest.json) records the complete
abstract-program digest, exact rotation parameters, graph edge list, and gate
counts. For `p=2`:

| Quantity | Value |
| --- | ---: |
| Species `q` | 14 |
| Logical modes `N` | 21 |
| Local circuit alphabet dimension | 317 |
| Charged records | 77 |
| Neutral auxiliary records | 81,265 |
| Total circuit sites | 81,342 |
| Abstract gates | 287,798 |
| Two-site nearest-neighbor gates after routing | 1,213,761,392 |

These deliberately unoptimized counts are finite and reproducible. They
highlight the distinction between circuit workspace/length and local
alphabet dimension. They do **not** give the complete microscopic qudit
size, which must also count all controller and wrapper factors. `p=2` is a
construction check, not a violating example.

The manifest's local count is calculated analytically from the exact swap
routes rather than materializing a billion records.
[RESOURCE_BOUND.md](RESOURCE_BOUND.md) derives numerical polynomial upper
bounds from the source loops, including the preloaded dispatcher tape. A separate one-edge
fixture checks this counting method against its entire generated local stream.

## Verification

Run from the repository root:

```
python3 contrib/controller_instantiation/circuit/verify.py
```

[`verify.py`](verify.py) writes [`results.json`](results.json). The current
checks include:

- Exhaustive occupation/branch inputs for four small graphs, including an
  isolated upper vertex. Both compressed port matrices agree with
  `epsilon f_edge^dagger` and
  `epsilon (1-1/(2N)) f_edge^dagger P`, with maximum floating-point discrepancy
  below `2.3e-16`.
- Exact integer checks of all neutral scratch cleanup, species counts, and
  occupation signs in 2,592 `p=2` input/branch cases.
- Separate sparse-state execution of the actual branch-preparation gates,
  including all workspace cleanup and non-power-of-two amplitudes.
- Full-space bijectivity and charge preservation of the local permutation
  families on the one-edge alphabet.
- Nearest-neighbor routing with arbitrary untouched registers, checking both
  operand orders and the last-site boundary case.
- An independent coherent check of the controlled-square-root-X Toffoli
  decomposition.
- **End-to-end coherent execution of the actual 36,140-gate local stream**
  for both basis vectors of the one-edge two-port fixture. Its 288 sites
  begin with BLANK auxiliaries and end with every auxiliary BLANK. The
  amplitudes agree with the target isometry within `1.2e-16`.

The Boolean and charge checks use exact integer arithmetic. Rotation and
state-vector amplitude checks use floating point; their mathematical gate
specification is algebraic and exact. These finite checks support, but do not
replace, the circuit identity argument.

## Dispatcher interface and remaining work

The local stream is a Python iterator:

```python
from sparse_circuit import projective_plane, make_circuit, Localizer
circuit = make_circuit(projective_plane(2))
for opcode, left_site in Localizer(circuit).stream():
    # Apply this opcode to (left_site, left_site + 1).
    pass
```

`LOCAL_OPCODES` and `LOCAL_INVERSES` expose the fixed menu and adjoints. Site
addresses are zero-based. The companion dispatcher may replace address changes
by literal cursor moves; it need not introduce an internal control state for
each address or species.

To regenerate only the small manifest:

```
python3 contrib/controller_instantiation/circuit/sparse_circuit.py --p 2 --manifest /tmp/p2-manifest.json
```

Optional `--abstract-jsonl PATH` and `--local-jsonl PATH` write the actual gate
streams. The latter is deliberately not used for the billion-gate `p=2`
instance in the repository.

Still missing from this circuit directory are an autonomous, initialized
length-to-program generator; its full transition-table composition with the
history clock and verifier; the total microscopic dimension and a numerical
violating parameter choice; and the independent infrared stability argument.
A preloaded script alone does not discharge those obligations.

**Provenance.** Implemented with AI assistance by Codex on 2026-10-28 and checked
by the scripts above. No human review is claimed.
