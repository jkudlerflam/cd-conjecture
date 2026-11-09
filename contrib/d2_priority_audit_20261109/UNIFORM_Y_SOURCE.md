# The transverse-DM symmetry permits an imaginary uniform source

AI-only analysis by OpenAI Codex (GPT-6 family), 2026-11-09. This corrects a
specific symmetry inference in `codex/work_C183/tune183.py` and gives a
finite-ring perturbative coefficient. It does not establish a bulk phase.

For the real-pairing, antisymmetric-facilitated family in C186's conventions,

\[
H=V\sum n_jn_{j+1}+\sum\{w b^\dagger_{j+1}b_j
+p b^\dagger_jb^\dagger_{j+1}+\mathrm{h.c.}\}
+f\sum(n_jX_{j+1}-X_jn_{j+1})-\mu\sum n_j,
\]

take real \(p,f\), \(w=(a+ib)/V\), and \(\mu=m/V\). The surviving
antiunitary is \(A=\mathcal P\mathcal I K\), with
\(\mathcal P=\prod_jZ_j\), inversion \(\mathcal I\), and conjugation \(K\).
Uniform \(X\) is odd, but uniform \(Y\) and \(Z\) are even. Therefore
\(\langle X\rangle=0\) does **not** test whether every uniform charge-one
source vanishes. For the vacuum-to-one-particle matrix element, the symmetry
says \(h_{1,0}=-h_{1,0}^*\): it may be purely imaginary. C183's check of the
real \(X\) component cannot establish its stated all-orders zero claim.

`check_transverse_symmetry.py` verifies this antiunitary directly for the
literal C183 Pauli bond. Its periodic small-ring results in
`transverse_symmetry.json` also exhibit nonzero uniform \(Y\), including a
large shell-dependent response at L=10. An expectation value is not itself
a source coefficient or a bulk phase diagnostic.

## Canonical effective-Hamiltonian calculation

Let \(\epsilon=1/V\) and \(H/V=H_0+\epsilon B+\epsilon^2 C\), with
\(H_0=\sum nn\), and let \(P\) project onto the blockade subspace. Write
the low-energy invariant subspace as the graph of \(X:P\to Q\), where

\[
H_{QP}+H_{QQ}X-XH_{PP}-XH_{PQ}X=0.
\]

Solving this equation recursively in \(\epsilon\), set
\(H_B=H_{PP}+H_{PQ}X\) and \(M=1+X^\dagger X\). The canonical Hermitian
effective Hamiltonian is

\[
H_{\rm can}=M^{1/2}H_BM^{-1/2}.
\]

This is the same polar-decomposition convention as C186's all-orders
`des_cloizeaux` function. `uniform_source_series.py` implements the recurrence
through order six. On periodic rings L=6,7,8 it gives

\[
\boxed{\langle j|H_{\rm can}^{\rm physical}|0\rangle
=i\frac{2p^3fb}{V^5}+O(V^{-6}).}
\]

Every coefficient of \(H_{\rm can}/V\) through order five is zero.
Equivalently, the leading term under this scaling is
\(h_{y,\rm eff}=2p^3 f\,\operatorname{Im}w/V^4\).
The displayed remainder is a fixed-finite-ring asymptotic statement; no
uniform-in-volume remainder bound has been proved here.

`uniform_source_exact.py` checks the coefficient using **exact Python integer
matrices and rational division**, not eigensolver subtraction. If
\(D=\operatorname{lcm}(1,\ldots,L)\), the graph coefficient can be stored as
\(X_n=A_n/D^n\), with Gaussian-integer \(A_n\). The final canonical coefficient
is evaluated with denominator \(8D^{n-1}\). All canonical Hermiticity residuals
are exactly zero. Saved results are in `uniform_source_exact.json`.

There is a short polynomial-completeness check for this finite-ring result.
The source vanishes when \(p=0\), is odd in \(f\) by particle parity, and
its imaginary part is odd in \(b\) by conjugation. Give \(p,f\) weight one
and \(a,b,m\) weight two. At weight six the only possible imaginary terms are

\[
i pfb(Ap^2+Bf^2+Ca+Dm).
\]

The exact evaluations \((p,f,a,b,m)=(1,1,0,1,0), (2,1,0,1,0),
(1,2,0,1,0),(1,1,1,1,0),(1,1,0,1,1)\) give
\(A=2,B=C=D=0\), separately for L=6,7,8. Lower possible imaginary weights
four and five are \(ipfb\) and \(ip^2fb\); their exact coefficients vanish.
The real part is forbidden by the antiunitary. These checks establish the
polynomial identity on these rings. Extending it to every L requires a
linked-cluster/support argument that is **not supplied here**.

The independent all-orders diagonalization `uniform_source_scan.py` checks
L=7,8 at V=8,16,32,64,128. For \(p=1,f=\sqrt{0.3}/2,b=1.035/2\), the leading
coefficient is \(2p^3fb=0.283446423505\ldots\). The L=8 imaginary matrix element
is approximately \(2.7898\times10^{-7}\) at V=16 and
\(8.3141\times10^{-12}\) at V=128; multiplying by V^5 approaches the predicted
coefficient. The projected-band singular value exceeds 0.94 throughout.
These are numerical cross-checks, not error-certified bounds.

## What extra tuning would mean

In Pauli conventions
\(p=J_x-J_y,\ f=2\Gamma_{xz},\ \operatorname{Im}w=2D_z,\ V=4J_z\),

\[
h_{y,\rm eff}=\frac{8(J_x-J_y)^3\Gamma_{xz}D_z}{V^4}+\cdots.
\]

For the C183 V=16 retune this is about \(+2.54\times10^{-7}\). Since
\(\langle1|Y|0\rangle=i\), adding a bare
\(h_y\sum_jY_j\) with \(h_y\simeq-2.54\times10^{-7}\) cancels the leading
vacuum-to-one matrix element in this effective-Hamiltonian convention.
One scalar real tuning suffices for this particular matrix element because
the real component remains symmetry-forbidden. This does not cancel all
density-dependent charge-one operators, prove emergent U(1), or prove that
such a cancellation is necessary for a two-mode liquid. Momentum selection
and the actual low-energy operator content still determine relevance.

The coefficient is tiny compared with \(|t_2|=1/(2V)=0.03125\). A large
\(\langle Y\rangle\) could reflect susceptibility, ordering, boundary response,
or a combination; it should not be converted directly into a bulk-gap claim.
In particular, on an open chain the pure Pauli DM bond and facilitated form
differ by edge fields:

\[
\Gamma\sum_{j=1}^{L-1}(X_jZ_{j+1}-Z_jX_{j+1})
=f\sum_{j=1}^{L-1}(n_jX_{j+1}-X_jn_{j+1})
+\frac f2(X_1-X_L),\qquad f=2\Gamma.
\]

The most useful immediate numerical work remains density/field tuning and
size/bond-dimension scaling of the actual microscopic chain. The tiny
counterterm is a controlled follow-up diagnostic if bulk evidence warrants
it, not a replacement for those checks.

Reproduce with the scientific Python environment used elsewhere in this
session, setting BLAS/OMP threads to one, and run the three `uniform_source_*`
scripts plus `check_transverse_symmetry.py`.
