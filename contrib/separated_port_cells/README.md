# Separate boundary gadgets with a shared logical register

## Claim

[proved] The abstract second-order cell construction does not require an
energy term that touches both boundary ports. Independently prescribed
left and right spectral weights can be realized by

\[
H_{\rm cell}=H_L+H_R,
\]

where each summand acts on its own qubit port, its own auxiliary register,
and a shared logical register. They need not commute. The ground space is
exactly the intended logical space, the cell gap is at least the minimum
of the two constituent gaps, and the second-order inter-cell bond is
unchanged. Exact logical fermion parity can be preserved.

This is a finite-dimensional matrix theorem. It is not a translation-invariant
nearest-neighbour implementation within a cell and does not establish any
Virasoro central charge or counterexample to c <= d - 1.

## Status and provenance

- Status: proved algebraically; finite-matrix numerical checks below.
- Author: Codex (GPT-6), 2026-10-28 UTC, during the research continuation
  requested by Pace.
- AI assistance: the parallel `graded_route_a` investigation supplies
  parity-compatible spectral weights. This contribution independently
  realizes them with a different physical tensor factorization.
- Independent AI review: the parallel `route_audit` worker independently
  checked the projector proof, parity-sector dimension count, and both
  numerical cases. It reproduced the quoted second-order errors and
  found no defect. This is AI review, not human verification.
- No human review is claimed.

## 1. Factorization lemma

Let the cell space be

\[
\mathcal B_L\otimes\mathcal B_R\otimes\mathcal K,
\qquad \dim\mathcal K=N.
\]

Choose unit vectors \(b_L\in\mathcal B_L\) and
\(b_R\in\mathcal B_R\). Write \(P_s=|b_s\rangle\langle b_s|\),
extended by identities on all other factors, and \(Q_s=I-P_s\).
Suppose \(H_s\) acts only on \(\mathcal B_s\otimes\mathcal K\), with

\[
\ker H_s=b_s\otimes\mathcal K,
\qquad H_s\ge\Delta_s Q_s,\qquad\Delta_s>0.
\tag{1}
\]

The shared logical register is unrestricted in each kernel. Consequently

\[
[H_L,P_R]=[H_R,P_L]=[P_L,P_R]=0.
\]

Also \(H_sP_s=P_sH_s=0\). The four sectors obtained from
\(P_L,Q_L\) and \(P_R,Q_R\) are therefore invariant under both summands.
Writing \(P=P_LP_R\) and \(\Delta=\min_s\Delta_s\),

\[
H_L+H_R\ge\Delta(Q_L+Q_R)\ge\Delta(I-P).
\tag{2}
\]

Both summands annihilate \(P\), so the ground space is exactly
\(b_L\otimes b_R\otimes\mathcal K\), of dimension \(N\), and the
gap is at least \(\Delta\). On the two single-excitation sectors,

\[
H_{\rm cell}|_{Q_LP_R}=H_L|_{Q_L},\qquad
H_{\rm cell}|_{P_LQ_R}=H_R|_{Q_R}.
\tag{3}
\]

No commutation of \(H_L\) with \(H_R\) is used. They can act
noncommutatively on \(\mathcal K\) in the doubly excited sector.

## 2. Boundary spectral data are preserved exactly

Factor each boundary register as a physical qubit port, a qubit partner,
and an ancilla of dimension \(a_s\), and choose

\[
b_s=|\Omega_2\rangle_{p_s,t_s}|0\rangle_{a_s}.
\]

For the three Pauli port operators \(F_a^s\), in any fixed ordering,
the source maps

\[
G_s|i\rangle=b_s|i\rangle,\qquad S_a^s=F_a^sG_s
\]

satisfy \(G_s^\dagger F_a^sG_s=0\) and
\((S_a^s)^\dagger S_b^s=\delta_{ab}I_N\).
Let independently prescribed positive weights \(M_e^s\) on the
\(3N\)-dimensional source space sum to identity, with positive
excited energies \(E_e^s\). Apply the one-port Naimark construction of
[G1](../unbounded_c/G1/G1_result.md) to realize

\[
(S^s)^\dagger\Pi_e^sS^s=M_e^s.
\tag{4}
\]

The complete cell ground isometry is
\(G|i\rangle=b_Lb_R|i\rangle\). Its left-port source lies in
\(Q_LP_R\), and its right-port source in \(P_LQ_R\). Equation (3)
then implies, for every \(z\ge0\),

\[
G^\dagger F_a^sQ(H_{\rm cell}+z)^{-1}QF_b^sG
=\sum_e\frac{(M_e^s)_{ab}}{E_e^s+z}.
\tag{5}
\]

Here \(Q=I-GG^\dagger\), and the inverse at \(z=0\) is restricted to
that complement. The corresponding cross-port matrix element is zero,
since the source sectors are orthogonal invariant subspaces. Coincident
eigenvalues across sectors do not change this statement.

For the physical inter-cell bond
\(v=\sum_{a=1}^3F_a^R\otimes F_a^L\), first order vanishes and the
second-order bond is exactly

\[
K_2=-\sum_{e,f,a,b}
\frac{(M_e^R)_{ab}\otimes(M_f^L)_{ab}}
{E_e^R+E_f^L}.
\tag{6}
\]

There is no conjugation on the second factor. This is the same formula
as C9. Substituting [the graded weights](../graded_route_a/README.md)
therefore gives a scalar plus any prescribed globally parity-even logical
bond. In a chain, distinct physical bonds have zero second-order cross
term because an outer cell is acted on once and its port compression
vanishes.

## 3. Dimension and parity

For \(k_s\) spectral weights the Naimark dilation needs \(3k_sN\)
dimensions orthogonal to the \(N\)-dimensional one-port code. A sufficient
boundary dimension and ancilla count are

\[
a_s=\left\lceil\frac{3k_s+1}{4}\right\rceil,
\qquad B_s=4a_s,
\qquad (B_s-1)N\ge3k_sN.
\tag{7}
\]

The combined cell dimension is

\[
D_{\rm sep}=B_LB_RN.
\tag{8}
\]

For identical numbers \(k\), this is
\(16\lceil(3k+1)/4\rceil^2N\). This is larger than C9's general
joint dilation, but its factorization resolves an interface issue: no
summand touches both ports. For \(N=2^n\) and the free-fermion target
with \(m\le2n+2\), take \(k=m+2\le2n+4\). Then
\(D_{\rm sep}=O(n^2 2^n)\), requiring \(n+O(\log n)\) qubits as an
abstract register. This statement is about memory, not circuit depth or
the number of microscopic local states.

For parity, assume the logical involution \(P_f\) has balanced positive
and negative eigenspaces, as it does for \(n\ge1\) fermions. Give the
boundary register grading \(\Gamma_s=Z_{p_s}Z_{t_s}\), with the
ancilla even. The Bell vector is even. If the weights commute with

\[
\operatorname{diag}(+1,-1,-1)\otimes P_f
\]

in the Pauli ordering \((Z,X,Y)\), their square roots and the Naimark
isometry intertwine this grading. Complete the dilation separately in
the two parity sectors. Equation (7) supplies sufficient dimensions in
each sector because all relevant spaces are balanced. Thus

\[
[H_s,\Gamma_sP_f]=0,
\qquad [H_{\rm cell},\Gamma_L\Gamma_RP_f]=0.
\tag{9}
\]

The global grading acts as \(P_f\) on the ground code. Each physical
inter-cell Pauli product is even. Hence exact total logical parity is
retained to all orders, not just in (6).

## 4. What this solves and what it does not

[proved] A compiler for this variant only needs to handle a sum of terms
that touch the left port or the right port, together with interior
registers. It does not have to synthesize a direct interaction between
the two physical boundary ports. The noncommutativity of the two
gadgets on shared memory causes no error in the ground code or the
second-order target.

[open] This does not construct that compiler. The one-port matrices are
still arbitrary interactions on their port and extended interior memory.
No constant-alphabet nearest-neighbour transition table, initialization
proof, or mobile-cell transport rule follows from equations (1)-(9).
Finite-clock errors and higher-order leakage must still fit inside the
stability margins of the logical theory. In particular, exactness of
this finite-matrix code does not imply exactness after a finite-clock
spatial simulation.

## 5. Numerical checks and reproduction

Run, from the repository root, with Python and NumPy:

```sh
OPENBLAS_NUM_THREADS=1 python3 contrib/separated_port_cells/checks.py
```

The script imports only the weight prescription from the adjacent
`graded_route_a` contribution, constructs separate parity-equivariant
one-port dilations, and embeds them into the shared-memory cell. It
diagonalizes that cell and evaluates the full two-cell resolvent in its
eigenbasis. It does not assume the second-order answer when performing
that evaluation. Tests include genuinely noncommuting \(H_L,H_R\).

| Target | Cell dimension | Two-cell dimension | Gap | Second-order error |
| --- | ---: | ---: | ---: | ---: |
| \(X\otimes Y\) | 288 | 82,944 | 1 | \(3.4\times10^{-15}\) |
| \(X\otimes Y+0.3Z\otimes Z\) | 512 | 262,144 | 1 | \(3.2\times10^{-15}\) |

The exact-parity commutators are zero in the numerical representation;
ground-code errors are below \(4.1\times10^{-14}\). Exchange-relative
errors are below \(3.6\times10^{-11}\), so the scalar offset is not
hiding an unresolved target coefficient. The output is
`checks_results.json`; the recorded run took less than one second.

## Relation to the project

This addresses the geometric interface in C9 and C19, complementing the
graded repair. It makes the remaining spatial construction more precise.
It changes neither the open status of c-d nor the requirement for a
fixed-alphabet core and controlled mobile-cell completion.
