# Exact gap on one whole wing of the staggered-charge family

## Statement and scope

Let \(L\geq4\) be even, impose periodic boundaries, and use Pauli matrices
with eigenvalues \(\pm1\). Consider

\[
H_L(a,c,b)=\sum_{j=1}^L
\left[a(X_jX_{j+1}-Y_jY_{j+1})+cZ_jZ_{j+1}+bZ_j\right].
\tag{1}
\]

The conserved additive charge is \(Q=\sum_j(-1)^jZ_j\): pair creation or
annihilation on a neighboring pair preserves this staggered charge.
Flipping alternate spins maps (1) to an XXZ chain in a staggered field.

**[Proved below.]** Set \(c_*=\sqrt{a^2+b^2/4}\). For \(c>c_*\), the
ground space consists of exactly the two Néel product states, its energy is
\(-cL\), and the spectral gap **above that two-dimensional ground space** is
exactly

\[
\boxed{\Delta_L=4(c-c_*)\qquad (L\text{ even},\ L\geq4).}
\tag{2}
\]

This exact finite-size gap holds for the entire indicated wing, at arbitrary
distance from its boundary. For \(a\ne0\), the boundary \(c=c_*\) has
exactly \(L+1\) ground states and an exact quadratic low-energy branch. It
does not have an ordinary finite-velocity CFT spectrum.

The boundary ground-space construction and its defect branch were already
proved in [C15, Section 2.4](../../codex/work_C15/RESULT.md). They are
rederived here to establish the gap throughout the adjoining wing, and are
not presented as new discoveries. The central charge of interest is the
total Virasoro charge of all propagating unitary compact CFT sectors.

## Positive decomposition and lower bound

Split the field evenly between bonds:

\[
h(a,c,b)=a(XX-YY)+cZZ+\frac b2(ZI+IZ).
\]

The two odd-parity states \(|01\rangle,|10\rangle\) have energy \(-c\).
On \(\operatorname{span}(|00\rangle,|11\rangle)\),

\[
h(a,c,b)+cI=
\begin{pmatrix}2c+b&2a\\2a&2c-b\end{pmatrix}.
\tag{3}
\]

Its eigenvalues are \(2c\pm2c_*\). Thus
\(K_*=H_L(a,c_*,b)+c_*L\) is positive semidefinite, and

\[
H_L(a,c,b)+cL
=K_*+(c-c_*)D_L,
\qquad D_L=\sum_j(I+Z_jZ_{j+1}).
\tag{4}
\]

Both summands annihilate the two Néel states. In a computational-basis
configuration, \(D_L\) equals twice the number of equal neighboring pairs.
On an even ring that number is even; it is zero exactly for the two Néel
configurations and otherwise at least two. If \(P_N\) projects onto their
span, then

\[
D_L\geq4(I-P_N),\qquad
H_L+cL\geq4(c-c_*)(I-P_N).
\tag{5}
\]

For \(c>c_*\), this proves both the exact two-dimensional kernel and the
gap lower bound in (2). The argument uses the even-ring parity constraint;
it is not an assertion about odd rings or open-chain boundary gaps.

## An invariant defect subspace saturates the bound

Fix a Néel state \(|N\rangle\) with alternating \(Z\)-eigenvalues
\(z_j=(-1)^j\), and define \(|j\rangle=X_j|N\rangle\). These \(L\)
orthonormal one-spin-defect states span an invariant subspace. Indeed, only
the two bonds neighboring the flipped spin have equal adjacent spins; the
pairing operator on either bond moves the defect by one site with amplitude
\(2a\). All other bonds annihilate that basis state. Relative to \(-cL\),

\[
(H_L+cL)|j\rangle
=(4c-2bz_j)|j\rangle+2a(|j-1\rangle+|j+1\rangle).
\tag{6}
\]

For two-site-cell momentum \(k=2\pi m/(L/2)\), this gives the Bloch matrix

\[
M(k)=\begin{pmatrix}
4c-2b&2a(1+e^{-ik})\\
2a(1+e^{ik})&4c+2b
\end{pmatrix},
\]

and exact excitation energies

\[
E_\pm(k)=4c\pm\sqrt{4b^2+16a^2\cos^2(k/2)}.
\tag{7}
\]

The allowed momentum \(k=0\) supplies an eigenstate orthogonal to both
Néel states with energy \(4c-4c_*\). Combined with (5), this proves equality
in (2). No dilute-gas approximation or claim about all other excitations
is needed.

## Boundary kernel and quadratic branch

Assume \(a\ne0\). Then \(c_*>0\) and the shifted local bond at the boundary
has rank one:

\[
h(a,c_*,b)+c_*I=4c_*|v\rangle\langle v|,
\qquad |v\rangle=A|00\rangle+B|11\rangle,
\tag{8}
\]

where

\[
A=\sqrt{\frac{2c_*+b}{4c_*}},\qquad
B=\operatorname{sgn}(a)\sqrt{\frac{2c_*-b}{4c_*}}.
\]

Both coefficients are nonzero and \(AB=a/(2c_*)\). Flip every odd site by
\(X\). The forbidden bras then alternate between
\(A\langle01|+B\langle10|\) and
\(A\langle10|+B\langle01|\). Apply the invertible diagonal filter
\(F=\operatorname{diag}(B,-A)\) on each even site of a trial state. On every
bond its forbidden bra becomes proportional to
\(\langle01|-\langle10|\). Hence the filtered constraints require symmetry
under every neighboring swap. Those swaps generate all permutations of
the sites, so their joint kernel is exactly the fully symmetric spin
subspace, of dimension \(L+1\). Inverting the filter and the odd-site
rotation proves the same ground-space dimension for (1). The filter is
used only to identify kernels; it is not a unitary spectral equivalence.

At \(c=c_*\), the lower branch in (7) is

\[
E_-(k)=4c_*-4\sqrt{c_*^2-a^2\sin^2(k/2)}
=\frac{a^2}{2c_*}k^2+O(k^4).
\tag{9}
\]

The \(k=0\) state joins the ground manifold, while the lowest nonzero
momenta give energies of order \(L^{-2}\). This branch and the growing
ground-space dimension rule out the ordinary finite-velocity conformal
description required by the conjecture. They cannot be counted as an
additional conformal sector.

To match C15's variables exactly, write \(n_j=(1-Z_j)/2\). Apart from a
constant, (1) has density interaction \(V=4c\), pair amplitude
\(\delta=2a\), and chemical potential \(\mu=4c+2b\). A global spin flip
can make \(b\leq0\), after which the boundary is
\(\mu=V-\sqrt{V^2-4\delta^2}\), as in C15. Its \(b=0\) endpoint is
included here; the same kernel proof remains invertible there.

## The other wing and the remaining two-parameter region

**[Proved using the primary-source input in the companion note.]** For
\(c\leq-|a|\) and \(b\ne0\), the
[rotation-invariant Suzuki–Fisher criterion](rank_one.md) applies with
\(n=\hat z\), \(\lambda=c\), and transverse block
\(\operatorname{diag}(a,-a)\). It gives a unique ground state and

\[
\Delta_L\geq2|b|.
\tag{10}
\]

Consequently, for \(a\ne0\) and \(b\ne0\), any as-yet-unexcluded
finite-velocity conformal point of this family must lie in

\[
-|a|<c<\sqrt{a^2+b^2/4}.
\tag{11}
\]

The slice \(c=0\) is quadratic Jordan–Wigner and is gapped for \(b\ne0\),
so it can also be removed. Its bulk quasiparticle energies are
\(\sqrt{4b^2+16a^2\sin^2p}\).

The signs of \(a\) and \(b\) can separately be reversed by uniform spin
rotations. Dividing by \(|a|\) leaves two parameters,
\(x=c/|a|\) and \(y=|b|/(2|a|)>0\). The remaining region is therefore
\(-1<x<\sqrt{1+y^2}\), with \(x=0\) removed. This is a reduction of the
open problem, not a proof that every critical point in that region is a
single Luttinger liquid.

## Exceptional parameters

* **\(b=0\), \(a\ne0\).** Alternate-site \(X\) rotations map (1) to
  \(a(XX+YY)-cZZ\), a zero-field XXZ chain. After making the XY coefficient
  positive, its anisotropy is \(-c/|a|\). The standard XXZ classification
  already used in the repository gives \(c_T=1\) for
  \(-|a|\leq c<|a|\); \(c=|a|\) is the quadratic ferromagnetic endpoint;
  the regions \(c>|a|\) and \(c<-|a|\) are gapped above their finite vacuum
  manifolds. In particular, (10) does not give a positive gap when \(b=0\),
  and the point \(c=-|a|\) is conformal.
* **\(a=0\).** The model is diagonal in the computational basis and has no
  propagating conformal modes. The proof of (2) remains valid for
  \(c>|b|/2\), but the boundary filter in (8) is singular. At
  \(c=|b|/2>0\), the boundary ground space instead consists of classical
  configurations avoiding one adjacent equal-spin pair and grows
  exponentially. At \(a=b=c=0\), the Hamiltonian is zero. These degenerate
  classical limits are not finite-velocity CFTs.

## Relation to the recent staggered-XXZ paper

**[Read; exact results and numerical evidence distinguished.]**
[Zhang and Song, arXiv:2609.29728v1](https://arxiv.org/html/2609.29728v1)
give the same degenerate boundary and polarized-wing gap. Their
spin-\(1/2\) normalization maps to ours with
\(J=4|a|\), \(\Delta=c/|a|\), and \(|h|=2|b|\), up to a constant.
Thus their \(\Delta=\sqrt{1+(h/J)^2}\) is precisely \(c=c_*\), and their
wing gap becomes (2). The surrounding Luttinger-liquid and gapped-region
phase diagram uses perturbation theory and DMRG; it is not imported here
as an exhaustive theorem. The proof above establishes the gap directly.

## Checks and provenance

Written on 2026-10-28 by a GPT-6-based Codex subagent. The local decomposition,
invariant subspace, and matching lower and upper bounds were independently
derived here. The boundary duplicates C15 and the wing gap agrees with the
primary paper cited above; no priority claim is made.

**[Numerical convention checks.]** An independent dense Kronecker-product
construction checked 24 full spectra: \(L=4,6\),
\((a,b)=(0.7,1.3),(-0.7,1.3),(0.7,-1.3),(0.7,0)\), and
\(c-c_*=0,0.013,0.4\). Every boundary had kernel dimension \(L+1\).
The largest error in (2) was \(3.56\times10^{-15}\); that in the complete
defect spectrum (7) was \(1.78\times10^{-15}\); invariant-subspace leakage
was zero. These are checks of factors and finite-size statements, not the
thermodynamic proof. Reproducible related tests are in the coordinator's
[checks.py](checks.py); the independent temporary Kronecker check is not a
second maintained script. No DMRG or human verification is claimed here.
