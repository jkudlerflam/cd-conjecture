# Exact closure of the on-site charge and rotation routes at d = 2

**Status: [proved] algebraic statements.** The statements about central charge
use the XXZ/XYZ infrared inputs already separated from algebraic proofs in
[C6](../../codex/work_C6_d2/c6.tex). They are not a proof for the general
interacting qubit chain.

The useful conclusion is that allowing arbitrary site-dependent charge axes
or spin rotations does not uncover a new uncontrolled zero-field class with
an additive on-site continuous symmetry. With a nonzero uniform field there
is exactly one extra class beyond the uniform U(1) reduction: the
pair-creation exchange with its staggered charge. Its interacting conformal
classification remains open in this argument.

## 1. Setup and local equations

Consider a periodic ring with L >= 4,

\[
H_L=\sum_{j=1}^L\left[\boldsymbol\sigma_j^T J
\boldsymbol\sigma_{j+1}+\mathbf B\cdot\boldsymbol\sigma_j\right],
\qquad J\in M_3(\mathbb R),\quad \mathbf B\in\mathbb R^3.
\]

Indices are cyclic. Constants and telescoping one-site bond terms have
already been removed. An additive on-site charge means

\[
Q_L=\sum_j\mathbf n_j\cdot\boldsymbol\sigma_j,
\qquad \mathbf n_j\in\mathbb R^3,
\]

with at least one nonzero coefficient. The coefficients need not have equal
norm, need not be periodic with a size-independent period, and may vanish at
some sites. This definition concerns an exact microscopic charge, not an
emergent continuum current or an operator with multispin densities.

Let \(\Omega(\mathbf n)\mathbf x=\mathbf n\times\mathbf x\). Independence of
Pauli strings on distinct sites and bonds gives the exact equivalence

\[
[Q_L,H_L]=0
\quad\Longleftrightarrow\quad
\boxed{\Omega(\mathbf n_j)J=J\Omega(\mathbf n_{j+1}),\qquad
\mathbf n_j\times\mathbf B=0\quad\text{for every }j.}
\tag{1}
\]

Thus no cancellation between a one-site commutator and a two-site
commutator, or between different bonds, has been assumed. The choice L >= 4
avoids small-ring conventions; it is enough for the thermodynamic question.
Open-chain end charges are not classified by the theorem below.

## 2. Exhaustive classification of the charges [proved]

The following list is exhaustive. All global changes of spin axes are
proper rotations.

### Rank zero

If J = 0, (1) only requires each \(\mathbf n_j\) to be parallel to B when
B is nonzero. If B = 0 all coefficients are arbitrary. These are product
Hamiltonians and supply no propagating conformal sector.

### Rank one

Write \(J=s\mathbf u\mathbf v^T\), with s != 0 and unit vectors u,v.
Equation (1) forces

\[
\mathbf n_j\parallel\mathbf u,\qquad
\mathbf n_{j+1}\parallel\mathbf v.
\]

Indeed, its two sides are
\(s(\mathbf n_j\times\mathbf u)\mathbf v^T\) and
\(-s\mathbf u(\mathbf n_{j+1}\times\mathbf v)^T\). Their left factors
are respectively orthogonal and parallel to u, so both vanish.
At each site the two neighboring bonds therefore require
\(\mathbf n_j\in\operatorname{span}(\mathbf u)\cap
\operatorname{span}(\mathbf v)\).

A nonzero charge exists precisely when u and v are collinear and, in
addition, B is parallel to their common axis or vanishes. In that case
\(J=\lambda\mathbf u\mathbf u^T\), and
\(\mathbf n_j=q_j\mathbf u\) with arbitrary real coefficients q_j.
These include coefficients that vanish. The Hamiltonian is a commuting
Ising exchange with a longitudinal field.

### Rank at least two, unequal singular values

A nonzero charge is possible only if J has a repeated singular value.
More precisely, its singular values must be \((s,s,t)\), with s > 0 and
t >= 0, t != s. Its left and right singular axes for the distinct value t
must coincide as unoriented lines. In axes with that line equal to z, J
then has the form

\[
J=\begin{pmatrix}sO&0\\0&\lambda\end{pmatrix},
\qquad O\in O(2),\quad |\lambda|=t.
\tag{2}
\]

Every charge is \(\mathbf n_j=q_j\mathbf e_z\), with

\[
q_{j+1}=\det(O)\,q_j.
\tag{3}
\]

In particular, all q_j are nonzero and have equal magnitude. The field must
be along z or vanish. The two orientations have different consequences:

* If det(O) = +1, the planar block is \(aI+b\Omega_2\). This is the
  usual XXZ exchange plus axial DM, with a uniform U(1) charge.
* If det(O) = -1, the planar block is a symmetric reflection times s.
  A common rotation around z diagonalizes it to diag(s,-s). This gives
  \(J=\operatorname{diag}(s,-s,\lambda)\), with staggered charge
  \(Q\propto\sum_j(-1)^jZ_j\). A nonzero charge on the ring requires
  even L.

If all three singular values are distinct, or if the distinct left and
right singular axes are different, no nonzero charge exists.

**Proof.** Transposing and multiplying the first equation of (1) gives

\[
[\Omega(\mathbf n_j),JJ^T]=0,
\qquad
[\Omega(\mathbf n_{j+1}),J^TJ]=0.
\tag{4}
\]

For rank J >= 2, a vanishing coefficient forces its neighbor to vanish:
\(J\Omega(\mathbf n)=0\) with nonzero n would put the two-dimensional
image of \(\Omega(\mathbf n)\) inside a kernel of dimension at most one.
The analogous statement applies in the other direction. Thus every
coefficient of a nonzero charge is nonzero.

A real symmetric matrix commuting with \(\Omega(\mathbf n)\), n != 0,
is scalar on the plane perpendicular to n and leaves its axis invariant.
When it is not scalar, n must be its unique distinct-eigenvalue axis.
Apply this to (4). The same site is the left endpoint of one bond and the
right endpoint of the preceding bond, so its nonzero coefficient must lie
along both distinct singular axes. They must coincide. In this common-axis
basis, the cross blocks of J vanish and its planar block has both singular
values equal to s, giving (2). Finally
\(\Omega_2O=\det(O)O\Omega_2\), which gives (3).

### Three equal nonzero singular values

Write

\[
J=\lambda R,\qquad R\in SO(3),\quad\lambda\ne0.
\tag{5}
\]

This includes improper orthogonal matrices: for J = sO, s > 0, take
\(\lambda=s\det O\) and \(R=(\det O)O\). The full solution of the bond
equation is

\[
\mathbf n_j=R^{-(j-1)}\mathbf n_1,
\qquad R^L\mathbf n_1=\mathbf n_1.
\tag{6}
\]

The field condition of (1) must also hold at every site. For B = 0 these
are the possible helical charges, including the uniform charge along a
rotation axis of R. In particular, the Hamiltonian already has a uniform
U(1), even when additional helical charges exist.

For B != 0, (6) is possible only if \(R\mathbf B=\mathbf B\) or
\(R\mathbf B=-\mathbf B\). The former gives a uniform charge and axial
XXZ plus DM. The latter requires a rotation by pi, gives an alternating
charge, and in axes with B along z again belongs to
\(J=\operatorname{diag}(s,-s,\lambda_z)\), after a planar rotation.
The alternating case requires even L. There is no further field-bearing
helical class.

**Proof.** For proper R,
\(R\Omega(\mathbf n)R^T=\Omega(R\mathbf n)\); substitute (5) into (1).
When B != 0 all coefficient vectors are parallel to B, and the rotation
preserves their nonzero norm. Its action on their common line is therefore
either +1 or -1. A proper three-dimensional rotation has eigenvalue -1
only at angle pi.

## 3. What this proves about the simplified conjecture

**Zero uniform field.** Every chain with a nonzero additive on-site charge
is either commuting, or reduces to XXZ by a product of on-site rotations,
with a possible boundary twist. The orientation-preserving case is the
standard spiral removal of axial DM. In the orientation-reversing case,
rotate alternate sites by pi about x:

\[
sXX-sYY+\lambda ZZ
\longmapsto sXX+sYY-\lambda ZZ.
\]

This is exact on even periodic rings. Consequently, using the standard
XXZ conformal classification already stated in C6, every unitary conformal
regime in this zero-field charge class has total Virasoro central charge
c_T <= 1. Product rotations preserve every bipartite entropy exactly.
Boundary twists require care with finite-size zero modes but do not add a
propagating bulk sector.

**Nonzero uniform field.** The uniform-charge class still reduces to XXZ
in a uniform longitudinal field. The remaining interacting class is

\[
H=\sum_j\left[a(X_jX_{j+1}-Y_jY_{j+1})+
cZ_jZ_{j+1}+BZ_j\right],\qquad aB\ne0.
\tag{7}
\]

The same alternating rotation produces

\[
\widetilde H=\sum_j\left[a(X_jX_{j+1}+Y_jY_{j+1})-
cZ_jZ_{j+1}+(-1)^jBZ_j\right].
\tag{8}
\]

This is XXZ with a **staggered** longitudinal field. Uniform-field Bethe
ansatz input cannot be applied to (8). If c = 0 the model is quadratic
after Jordan-Wigner and is already controlled. The general interacting
c != 0, B != 0 case is not settled by this classification. C6 already
exhibited (7); the new statement is its exhaustiveness among chains with
a nonuniform additive on-site charge and a uniform nonzero field.

The classification does not require an exact microscopic charge to exist
at every critical point. It gives no obstruction to emergent conserved
currents, multispin conserved densities, or neutral CFT sectors in models
without such a charge.

## 4. Arbitrary spin rotations do not enlarge the quadratic class [proved]

Suppose a product of arbitrary on-site rotations makes the Hamiltonian
quadratic under a standard parity-preserving Jordan-Wigner map. The
quadratic spin form has only transverse-transverse exchange and a field
along the local Jordan-Wigner z axis. Let n_j be that unit axis expressed
in the original physical spin coordinates. Then

\[
J\mathbf n_{j+1}=0,\qquad J^T\mathbf n_j=0,
\qquad\mathbf B\parallel\mathbf n_j\quad(B\ne0).
\tag{9}
\]

At every site the preceding bond also requires \(J\mathbf n_j=0\).
Therefore

\[
\mathbf n_j\in\ker J\cap\ker J^T.
\]

Pick any one such axis. It is a single uniform Jordan-Wigner axis for the
original Hamiltonian, including its field. Conversely, a common unit null
axis with B along it gives the quadratic form immediately. Hence arbitrary
periods, spirals, and nonperiodic choices of product rotations do not enlarge
the exact quadratic class

\[
\mathcal F=\{(J,B):\exists\,|\mathbf n|=1,
\ J\mathbf n=J^T\mathbf n=0,\ B\times\mathbf n=0\}.
\]

This statement concerns this direct local-rotation/Jordan-Wigner route.
It does not classify nonlocal dualities or finite-depth entangling circuits.

## 5. Product-rotation closure of homogeneous zero-field XYZ [proved]

Suppose arbitrary on-site rotations \(R_j\in SO(3)\) transform the
uniform zero-field bond J into the same symmetric bond K on every bulk
edge:

\[
R_jJR_{j+1}^T=K=K^T.
\tag{10}
\]

At a site with both neighboring bulk edges, multiplying the outgoing
bond by its transpose and the incoming bond's transpose by that bond gives

\[
R_jJJ^TR_j^T=K^2=R_jJ^TJR_j^T.
\]

Thus J must be **normal**, \(JJ^T=J^TJ\). No invertibility or
nondegeneracy assumption is needed.

Conversely, every real normal 3-by-3 J has the required bulk reduction.
Write J = S + A, with S symmetric and A antisymmetric. Normality is
equivalent to [S,A] = 0. If A = 0 a common spin rotation diagonalizes S,
giving XYZ. If A != 0, its kernel is its DM axis, and commutation forces
S to be scalar on the perpendicular plane. This is precisely axial XXZ
plus DM. Spiral rotations remove the DM in the bulk and produce a fixed
symmetric XXZ bond. A periodic ring can retain a boundary twist, so the
converse is a bulk statement or a statement on compatible rings.

Therefore the product-rotation closure of homogeneous zero-field XYZ is
exactly the normal-matrix locus. It is already the union of zero-field XYZ
and the zero-field axial U(1) class. A nonnormal exchange matrix cannot be
made into a homogeneous XYZ model by any choice of on-site spin axes.
Allowing an inhomogeneous target bond is outside this lemma.

## Provenance and checks

Derived on 2026-10-28 with Codex (a GPT-6-based assistant). Two parallel agents
independently derived the additive-charge classification and the quadratic
closure statement, and compared rank-deficient cases, vanishing charge
coefficients, improper rotations, and periodic compatibility. This is an
internal algebraic cross-check, not external peer review. The normal-matrix
closure proof was checked directly by multiplying neighboring bond
identities. No new lattice-to-CFT theorem is asserted. The parent
contribution supplies reproducible numerical checks of the linear charge
equations separately.
