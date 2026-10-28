# Hidden free dualities do not enlarge the nearest-neighbor qubit class

## Claim and status

**[Proved algebraically below; existing exact free-chain input for the central-charge conclusion.]** Consider a uniform qubit chain on the full tensor product, with Hermitian nearest-neighbor Hamiltonian

\[
H=\sum_j\left[\sum_{a=X,Y,Z}B_a a_j+
\sum_{a,b=X,Y,Z}J_{ab}a_jb_{j+1}\right],\qquad B_a,J_{ab}\in\mathbb R.
\]

Suppose there is a faithful term-by-term fermionization in which each nonzero Pauli monomial in this expression becomes a nonzero scalar multiple of a **single Majorana bilinear**. The fermionization may be nonlocal, and the Majorana couplings need not be nearest-neighbor. Then the original Hamiltonian already admits a uniform on-site rotation to either:

1. an ordinary parity-preserving nearest-neighbor Jordan–Wigner quadratic chain; or
2. a commuting longitudinal Ising chain.

Consequently this construction cannot produce a finite-velocity compact unitary CFT with total Virasoro central charge \(c_T>1\). The free branch uses the exact quasi-free bound established in [C1](../../../codex/work_C1_astra/); the commuting branch has no propagating CFT modes.

In particular, a faithful Clifford duality followed by parity-preserving Jordan–Wigner cannot conceal a higher-central-charge quadratic chain inside a nearest-neighbor qubit Hamiltonian. This extends the repository's product-rotation closure to a precise class of nonlocal transformations. It is a closure result for an attempted counterexample mechanism, **not a proof of the full qubit conjecture**.

## Why a claw obstructs bilinears

Write a Majorana bilinear as \(i\gamma_r\gamma_s\), \(r\ne s\). Two distinct bilinears anticommute exactly when their unordered endpoint pairs have exactly one endpoint in common. Thus their anticommutation graph is a line graph.

Here a *claw* means four distinct operators \(C,A,B,D\), with \(C\) anticommuting with each of \(A,B,D\), while the three leaves commute pairwise. Such a configuration is impossible for distinct single bilinears: each leaf must use exactly one endpoint of the center; by the pigeonhole principle two leaves use the same endpoint. Since those leaves are distinct, they share exactly that endpoint and anticommute, a contradiction.

A faithful algebra map preserves this contradiction. Identifying two distinct local operators after imposing additional constraints is outside the theorem's hypothesis.

## Six local obstructions

The support of \(H\) is the subset of the twelve labels

\[
\{X,Y,Z,XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ\}
\]

whose coefficients are nonzero. A two-letter label \(ab\) at index \(j\) denotes \(a_jb_{j+1}\). Uniformity supplies all translates of every supported label.

The following six support templates each force a claw. Every row is checked directly from Pauli commutation, and is valid on the infinite chain or a sufficiently large periodic ring.

| Forbidden support | Center | Three pairwise commuting leaves |
|---|---|---|
| \(X,XY\) | \((XY)_0\) | \(X_1,(XY)_{-1},(XY)_1\) |
| \(XY,ZZ\) | \((XY)_0\) | \((XY)_{-1},(XY)_1,(ZZ)_{-1}\) |
| \(XY,YZ\) | \((XY)_0\) | \((XY)_{-1},(XY)_1,(YZ)_{-1}\) |
| \(X,Y,XX\) | \(Y_0\) | \(X_0,(XX)_{-1},(XX)_0\) |
| \(X,XX,YY\) | \((YY)_0\) | \(X_0,X_1,(XX)_{-1}\) |
| \(XX,YY,ZZ\) | \((XX)_0\) | \((YY)_{-1},(YY)_1,(ZZ)_{-1}\) |

Every permutation of the three axes and every reflection of these templates is forbidden as well. Reflection reverses each two-letter label. These transformations produce 37 distinct inclusion-minimal forbidden subsets: 24 pairs and 13 triples.

This argument does not use a numerical graph-recognition test. The saved exact enumeration independently checks its finite combinatorics.

## Exhaustion of the remaining supports

Suppose first that an off-diagonal exchange is present; rename axes so it is \(XY\).

* The first three templates exclude the fields \(X,Y\), the exchange \(ZZ\), and the directed-chain exchanges \(YZ,ZX\).
* The only remaining possible exchanges, besides \(XY\), are \(XX,YY,YX,XZ,ZY\).
* If \(XZ\) is also present, the same templates exclude the remaining field \(Z\), and exclude \(YY,YX,ZY\). The support is therefore contained in the single row \(\{XX,XY,XZ\}\).
* If \(ZY\) is present, reflection gives the analogous single column \(\{XY,YY,ZY\}\).
* Otherwise the support is contained in \(\{Z,XX,XY,YX,YY\}\), the standard quadratic plane with perpendicular field.

Now suppose no off-diagonal exchange is present.

* Three diagonal exchanges are forbidden by the last template.
* With exactly two diagonal exchanges, say \(XX,YY\), the fifth template excludes fields \(X,Y\). Only \(Z\) can remain. This is again the standard quadratic plane.
* With only \(XX\), the fourth template says that the field either lies along \(X\), or lies entirely in the transverse \(YZ\) plane. These are respectively the commuting Ising family and a transverse-field Ising family after a uniform rotation around \(X\).
* With no exchanges, arbitrary fields give independent on-site spins.

This exhausts all supports. The maximal allowed sets consist of the following sixteen sets:

| Type | Number of maximal supports | General form |
|---|---:|---|
| Arbitrary on-site field | 1 | \(\{X,Y,Z\}\) |
| Longitudinal commuting Ising | 3 | \(\{a,aa\}\) |
| One diagonal exchange and transverse field | 3 | \(\{b,c,aa\}\), \(\{a,b,c\}=\{X,Y,Z\}\) |
| One exchange row | 3 | \(\{aX,aY,aZ\}\) |
| One exchange column | 3 | \(\{Xa,Ya,Za\}\) |
| Quadratic plane and normal field | 3 | \(\{a,bb,bc,cb,cc\}\) |

A single exchange row has \(J=u v^T\) and \(B=0\), and a column has the transpose form. Choose a unit vector \(n\) perpendicular to both \(u\) and \(v\). Then

\[
Jn=J^Tn=0.
\]

Rotate \(n\) to the Jordan–Wigner quantization axis. Every exchange lies in its transverse plane, so the Hamiltonian is quadratic. This works also when \(u\) and \(v\) are parallel, by choosing any perpendicular \(n\).

For a single diagonal exchange and transverse field, choose \(n\) parallel to the field (or any perpendicular axis when the field vanishes). Again \(Jn=J^Tn=0\), and the field lies along \(n\). The remaining cases are explicitly quadratic or commuting. This proves the theorem.

## Why the Clifford corollary follows

A faithful Clifford map sends distinct Pauli monomials to distinct Pauli monomials, up to signs. A parity-preserving quadratic Majorana Hamiltonian has a Pauli expansion consisting of the distinct strings that represent its bilinears. Linear independence of Pauli monomials prevents cancellation of a nonquadratic image against an image of another supported term. Thus, if the transformed Hamiltonian is quadratic, its nonzero input terms satisfy the bilinear hypothesis above.

The same argument applies to faithful bulk Kramers–Wannier-type Pauli algebra maps wherever their images satisfy this hypothesis. Ordinary finite global parity/twist bookkeeping does not invalidate a local claw on four sites. Additional local projections, auxiliary gauge variables, or representations in which a Pauli term becomes a **sum** of bilinears require a separate argument and are not included.

Uniformity of the original coefficients matters: it ensures that all translates appearing in each claw occur in the Hamiltonian. A period-two compass decomposition, an enlarged unit cell, or a nonuniform chain must be analyzed separately.

## Reproduction and independent checks

Run from the repository root:

```sh
python3 contrib/direct_attack_20261028/alternatives/check_supports.py
```

The standard-library script checks all six displayed claws, independently constructs the Pauli anticommutation graph on a seven-site ring, and enumerates all 4,096 supports. It finds 37 minimal forbidden supports, 113 allowed supports, and exactly the sixteen maximal supports listed above. All arithmetic is Boolean/integer; there are no numerical tolerances. Output is [checks.json](checks.json). Runtime is below one second on the recorded machine.

The enumeration is a consistency check on the local analytic proof, not an inference about an unknown thermodynamic phase. No new central charge was fitted, and no independent human or second-agent verification is claimed in this file.

## Consequence for the direct search

This removes a concrete route to a simple \(d=2\) counterexample: taking a known free multichannel chain and attempting to compress its interaction range by a faithful Clifford or Pauli-monomial duality. A successful small-\(d\) alternative must evade at least one hypothesis, for example by being genuinely interacting or by using a non-Pauli encoding with a controlled low-energy reduction.

The repository's candidate at \(d=24\) deliberately uses such a low-energy encoding and is unaffected. Its finite-coupling matching problem remains a stronger route to an actual violation than the free-duality shortcut ruled out here. This is a research assessment, not a theorem that no other simpler construction exists.

Written 2026-10-28 by a GPT-6-based Codex research agent. The proof and exact support enumeration were developed in this run. The existing nearest-neighbor free-fermion bound and the product-rotation criterion are inherited from [C1](../../../codex/work_C1_astra/) and [the qubit reductions](../../d2_gap_reduction/README.md). No priority claim outside this repository is made.
