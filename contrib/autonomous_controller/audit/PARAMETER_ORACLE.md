# Reading the chosen prime from the local matrices

The Hamiltonian family already permits its local matrices to depend on the
chosen prime `p`: the species alphabet, exact rotation blocks and indexed
query gates do so. Consequently **the controller need not recover `p` from
the physical cell length**. Two additional fixed-cardinality neutral gate
families can expose it through the existing indexed helper. This removes
one unnecessary compilation task; it does not generate the sparse circuit.

Status: the local gates and the extraction/cleanup circuit are explicit
whole-space permutations. A uniform compact-walk algorithm is supplied.
Its finite head/collision table is **not** instantiated here. This is a
separate interface simplification, not a claim that the complete autonomous
high-channel Hamiltonian has been compiled.

## 1. The two parameter-dependent gates

Use the already present neutral records `H_(j,c)`, `0 <= j < m`, `c in {0,1}`.
Let `p_j` be bit `j` of `p`, with zero padding. Define

\[
 Q_p\lvert H_{j,c}\rangle=\lvert H_{j,c\oplus p_j}\rangle,\qquad
 E_m\lvert H_{j,c}\rangle=\lvert H_{j,c\oplus[j=m-1]}\rangle.
\]

Both gates fix every non-H alphabet state. They are involutive permutation
matrices and commute with each species charge on the **whole** alphabet.
Each is one gate family, independently of `p` and `m`. Neither adds an
index register, a species register, or a control-state label for every `p`.
One-site gates can be represented as two-site gates tensored with identity.

For the current exact circuit,

\[
q=2(p^2+p+1),\quad N=(p+1)(p^2+p+1),\quad
m=\max(1,\lceil\log_2q\rceil,2\lceil\log_2N\rceil),
\]

so `m >= bit_length(p)` for all primes `p >= 2`. The complete work alphabet
remains `D(q)=q+3+2(q+1)m`. Two neutral digit/end tracks multiply this by
four if they are added afresh; existing trusted neutral scratch can supply
them. The uncompiled control overhead is not included in that factor.

## 2. Extract without trusting the other quantum records

The input pins already provide one `H_(0,0)` helper at RIGHT and zeroed
trusted neutral tracks. Assume the segment has at least `m` sites. Other
Sigma records are arbitrary, including malformed H and loaded T records.

At the helper's current site:

1. Apply `Q_p`, copy the helper bit into the local neutral digit track,
   and apply `Q_p` again. The helper bit is again zero.
2. Apply `E_m`, copy its bit into the local neutral end track, and apply
   `E_m` again. The helper bit is again zero.
3. If the trusted end bit is zero, increment the helper index with SHIFT,
   exchange its Sigma record with the left neighbor, and repeat.
4. When the end bit is one, transport the helper back to RIGHT by Sigma-only
   SWAPs, decreasing its index after each rightward move.

The end decision is a deterministic function of the known helper index;
it never depends on a logical occupation, species label, or coherent
encoding branch. The neutral output tracks stay at their physical sites.
No other H-like record is mistaken for the helper: its tensor factor is
followed by the prescribed adjacent swaps, as in the existing format audit.

The result is `p` written little-endian from the right endpoint, with one
end delimiter `m-1` sites to its left. The helper is back at `H_(0,0)` and
**every untrusted Sigma tensor factor is restored to its original site**.
Linearity makes the statement hold for arbitrary entangled quantum inputs.
The inverse word clears all parameter output and restores the input exactly.

The expanded reference word has `6m+4(m-1)=10m-4` primitive gates. The
compact algorithm needs neither a stored instruction per physical site nor
an index value in classical finite control; END supplies termination and
the immutable RIGHT mark supplies the return stopping condition. This is
an algorithmic resource statement, not a count of an expanded machine table.

The full word is a unitary on arbitrary inputs. The described adaptive
execution schedule relies on the stated helper and neutral-track pins.
It is not asserted to have the same trajectory on uninitialized tracks.

## 3. What this removes, and what remains

The circuit generator can take these extracted bits as its parameter
input. It no longer needs an integer-root routine for `n=C_0 p^K`, or a
length counter merely to learn `p`. Cell length still matters for fitting
the workspace, bounding runtime, selecting a cell size in the wrapper,
and checking every untrusted padding record. Compact endpoint scans and
the independent blank-to-length controller remain useful for those tasks.

The construction uses a **different** microscopic interaction for each
`p`, as the original sparse family already does. It does not claim a fixed
interaction or a fixed local dimension that can acquire unbounded `p`.
The added lookup matrices have finitely many explicitly specified 0/1
entries; no real-number digit extraction or precision assumption is used.

The remaining work is to compile extraction, graph/circuit generation,
format verification, execution and readout into one finite reversible
head table, include all its factors in the local dimension, and establish
the intended finite-coupling infrared behavior.

## 4. Reproduction and audit status

Run `python3 contrib/autonomous_controller/audit/parameter_checks.py`.
The checks enumerate every local gate column on the complete alphabet,
including loaded and malformed symbols, for five projective-plane sizes.
They exhaust the entire promised-helper basis of a three-site example,
including arbitrary values of both neutral target tracks. Larger random
malformed payloads test the END-based stopping algorithm and exact cleanup.

The source and proof were developed by a GPT-6-family audit agent in this
session. No independent human review is claimed.
