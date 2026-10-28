# A complete ferromagnetic rank-one subclass, and a larger gapped region

## Claim and status

Consider periodic qubit chains, with Pauli matrices of eigenvalues ±1,

\[
H_L=\sum_{j=1}^L\left[
\boldsymbol\sigma_j^T J\boldsymbol\sigma_{j+1}
+B\cdot\boldsymbol\sigma_j\right],\qquad J\in\mathbb R^{3\times3}.
\]

The central charge below is the total Virasoro central charge of the
propagating sectors of a unitary compact CFT with nonzero limiting velocities.

**[Proved, using the cited preprint gap theorem stated below.]** If
\(J=-\kappa nn^T\), \(\kappa\geq0\), \(\|n\|=1\), then arbitrary uniform
fields \(B\) obey \(c_T\leq1/2\) whenever a conformal scaling limit exists.
For \(\kappa>0\), the only conformal parameters are

\[
B\cdot n=0,\qquad \|B\|=\kappa,
\]

where the transverse-field Ising solution has \(c_T=1/2\).

This is a complete statement for **symmetric ferromagnetic rank-one coupling
matrices**. It does not classify general rank-one matrices \(J=uv^T\), nor
the antiferromagnetic sign \(J=+\kappa nn^T\).

## Primary-source input

**[Read.]** Bravyi, Gosset, Liu and Wong, *Efficient quantum algorithm for
Heisenberg spin systems*, [arXiv:2607.14401v1, Definition 1 and Theorem 1](https://arxiv.org/html/2607.14401v1),
prove the following finite-system result. Write

\[
H=-\sum_{(i,j)}\left[K_{ij}Z_iZ_j+
(X_i,Y_i)A_{ij}\binom{X_j}{Y_j}\right]
-\sum_i(\mu_i^xX_i+\mu_i^yY_i+\mu_i^zZ_i).
\]

If \(K_{ij}\geq\|A_{ij}\|_{\rm op}\) and
\(\mu_i^z\geq\mu>0\), then the ground state is unique and
\(E_1-E_0\geq2\mu\). The transverse fields are unrestricted. The result
holds on arbitrary finite graphs and therefore on periodic chains. Here
\(E_0\leq E_1\leq\cdots\), counting multiplicity.

The paper supplies the gap theorem; this note does not supply a new proof of
that theorem. The deductions below are rotations and exact Ising analysis.

## A rotation-invariant sufficient condition

**[Proved from the primary-source input.]** Suppose there is a real unit vector
\(n\) and real \(\lambda\) satisfying

\[
Jn=J^Tn=\lambda n,\qquad
\lambda\leq-\|PJP\|_{\rm op},\qquad P=I-nn^T.
\tag{1}
\]

For every \(B\) with \(B\cdot n\ne0\),

\[
\boxed{\quad E_1(H_L)-E_0(H_L)\geq2|B\cdot n|>0\quad}
\tag{2}
\]

uniformly in \(L\).

Proof. The two eigenvector identities give both vanishing cross blocks,
\(PJn=0\) and \(n^TJP=0\). Consequently
\(J=PJP+\lambda nn^T\). Choose an oriented orthonormal basis with last
axis either \(n\) or \(-n\), so that the field coefficient along its last
axis is \(-|B\cdot n|\). This basis change is implemented by the same
single-qubit unitary at every site. The transformed bond is a transverse
\(2\times2\) real block plus \(\lambda ZZ\). In the sign convention of the
cited theorem take \(K=-\lambda\), and take \(A\) to be minus that
transverse block. Its operator norm is \(\|PJP\|_{\rm op}\), so (1) supplies
the required dominance inequality. The transformed longitudinal field is
\(\mu^z=|B\cdot n|\); the other field components are allowed. The unitary
does not change the spectrum. This proves (2).

The bound excludes a finite-velocity conformal tower with positive central
charge. It holds even when the uniform field is misaligned with every
interaction axis and the spin Hamiltonian has no nontrivial uniform unitary
internal symmetry. Thus this criterion removes interacting points outside
the earlier free, uniform-U(1), and zero-field-XYZ loci.

For symmetric \(J\), condition (1) is particularly easy to check: a
negative eigenvalue must have magnitude at least that of each other
eigenvalue, and \(B\) must have nonzero projection onto a corresponding
eigenvector. If the dominant negative eigenspace is degenerate, any vector
in it with nonzero field projection suffices.

## Ferromagnetic rank-one proof and endpoints

Take \(J=-\kappa nn^T\) with \(\kappa>0\). A uniform rotation first makes
\(n\) the \(z\) axis. A rotation about \(z\), followed if needed by a
\(\pi\) rotation about \(x\), puts the Hamiltonian in the form

\[
H_L=-\kappa\sum_jZ_jZ_{j+1}-g\sum_jX_j-h\sum_jZ_j,
\quad
g=\|B-(B\cdot n)n\|\geq0,
\quad h=|B\cdot n|\geq0.
\tag{3}
\]

* If \(h>0\), (1) holds because \(PJP=0\). Equation (2) gives a unique
  ground state and a size-independent gap at least \(2h\), for every
  transverse field \(g\).
* If \(h=0\), (3) is exactly the transverse-field Ising chain. Its
  Jordan–Wigner quasiparticle dispersion is
  \[
  \varepsilon(k)=2\sqrt{\kappa^2+g^2-2\kappa g\cos k}.
  \tag{4}
  \]
  For \(g\ne\kappa\), it has a positive bulk excitation gap. In the
  ordered regime \(0\leq g<\kappa\), the two symmetry-breaking vacua
  require care: finite-ring splitting tends to zero, but the gap above
  this finite vacuum manifold stays positive. This is not a conformal
  spectrum. At \(g=\kappa>0\), (4) becomes
  \(4\kappa|\sin(k/2)|\). Its single massless real fermion gives the
  Ising CFT, \(c_T=1/2\), with velocity \(2\kappa\). The exact
  Ising solution and its continuum identification are the standard input
  already used in [the repository's qubit analysis](../../codex/work_C6_d2/c6.tex).

At \(\kappa>0\) and \(B=0\), the two all-aligned product states are the
ground states; on a ring of at least three sites, the first domain-wall
pair costs \(4\kappa\). There is no propagating conformal sector.

At \(\kappa=0\), spins are independent. For \(B\ne0\), the ground state
is a unique product and the gap is \(2\|B\|\). For \(B=0\), the
Hamiltonian is zero and every state is a ground state. This extensively
degenerate zero operator has no finite-velocity CFT scaling limit; assigning
a central charge from an arbitrarily chosen entangled ground state would
not meet the conjecture's spectral hypothesis.

## What this does not establish

For \(J=+\kappa nn^T\), rotating alternate sites reverses the Ising bond
sign but makes a nonzero uniform longitudinal field staggered. The
positive-field hypothesis of the gap theorem then fails. The mixed-field
antiferromagnetic Ising problem is therefore not closed by this argument.

For a nonsymmetric rank-one matrix \(J=\alpha uv^T\) with nonparallel
\(u,v\), a common eigenvector of \(J\) and \(J^T\) with nonzero eigenvalue would have to
lie in both \(\operatorname{span}(u)\) and \(\operatorname{span}(v)\),
which is impossible. A common null vector lies perpendicular to both;
then \(PJP=J\ne0\), violating (1). Thus the displayed gap criterion does
not accidentally classify this additional family. If the field is
perpendicular to both \(u\) and \(v\), that family instead lies in the
previously controlled quadratic Jordan–Wigner locus; general fields remain
outside the present proof.

## Provenance and checks

Written on 2026-10-28 by a GPT-6-based Codex subagent. The cited primary
paper was opened at its definitions, theorem, and proof. The algebraic
reduction was derived here; the coordinating Codex agent independently
opened the primary theorem. No human review or new numerical calculation
is claimed in this file. The finite-size gap and rotation checks are
maintained separately by the coordinator in this contribution folder.
