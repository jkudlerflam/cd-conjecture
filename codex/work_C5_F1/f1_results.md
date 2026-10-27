# F1. Partons, gauge subtraction, and the missing implication

This is a closed-book analysis. No pre-existing files or external sources were read. The accompanying `f1.py` checks finite-dimensional representation data and rational central-charge formulas. It does not establish a continuum limit.

The result is conditional. An identified free-fermion coset gives a useful bond-rank bound, and the minimal one-parton representation gives `c <= d-1` within that class. Several larger representations obey the same bound. However, the assertion that every representation of a local `C^d` obeys the corresponding formal coset bound is false. A weighted U(1) representation of `C^2` has a formal coset with arbitrarily large central charge. An explicit current identity shows why that coset cannot describe its hard projection with the usual identification of operators.

This is a failure of the proposed representation-to-coset inference, not a counterexample to the c-d conjecture.

## 1. Scope, conventions, and assumptions

Write `c` for the chiral central charge of a nonchiral theory with `c_L=c_R`. For a product of gapless sectors, `c_T` is the sum in the question. A complex massless Dirac fermion contributes 1 and a nonchiral Majorana contributes 1/2. A gapped sector contributes zero.

Affine levels use the normalization in which long roots have squared length 2 and a complex SU(N) fundamental gives level 1. Abelian levels are stated explicitly through their current OPEs. Global quotients of gauge groups can affect sectors but do not alter the central-charge formulas used below.

A local fermionic parton representation consists of a finite Fock space `F_x`, a gauge action, and a specified charge or invariant sector

\[
\mathcal H_x=P_x\mathcal F_x\cong\mathbb C^d.
\]

The sector can transform by a fixed character before including a background gauge charge. Thus `N_x=1` is the invariant sector of the generator `N_x-1`. Constant background charges do not change current OPE levels. A representation is kinematic information. It specifies neither a Hamiltonian on the code nor its IR.

The following distinctions are necessary.

1. A free mean-field Hamiltonian is usually not invariant under independent gauge transformations at every site. Its hopping leaves the code. It is an ansatz in an enlarged space.
2. Hard projection means applying `P=prod_x P_x` to a state, or restricting operators to the code.
3. Lattice gauging introduces link variables and Gauss laws. Its Hilbert space generally includes electric flux states as well as matter.
4. A continuum current-algebra coset is a particular conformal construction. It is not the definition of either operation in items 2 and 3.

The common claim that these operations have the same IR is an additional dynamical assertion.

### Theorem 1. Conditional coset and bond-rank theorem

Consider a nearest-neighbor quadratic parton Hamiltonian, with any finite period, satisfying the supplied linear-crossing hypothesis. Let its actual gapless CFT be `F`, with central charge `c_F`. Suppose the following conditions hold.

(A1) The gapless fields possess a unitary affine gauge-current subalgebra

\[
\mathcal A_g=\left(\bigotimes_\alpha
\widehat{\mathfrak g}_{\alpha,k_\alpha}\right)
\otimes\widehat{\mathfrak u}(1)^u.
\]

Each simple level is positive. The Abelian level matrix on the `u` currents is positive definite. These are currents of the gapless theory, not just generators acting on massive orbitals. Vector gauging is free of local and global gauge anomalies. In particular, left and right levels match for each gauged factor.

(A2) The stress tensor splits as a unitary coset stress tensor plus the Sugawara gauge stress tensor,

\[
T_F=T_{F/\mathcal A_g}+T_g,
\qquad
c_g=\sum_\alpha\frac{k_\alpha\dim\mathfrak g_\alpha}
{k_\alpha+h^\vee_\alpha}+u.
\]

This is a conformal embedding statement about the actual massless channels. The notation `U(n_c)_1` is appropriate only when all `n_c` complex channels are present with the required conformal structure. General paired or partially gapped bands require their actual free CFT instead.

(A3) The interacting physical chain has an IR equivalent to the nonchiral coset `F/A_g`, possibly times gapped or finite topological sectors, and with no additional gapless sectors. This is an assumption about the chain, not a consequence of the on-site parton representation.

Then

\[
c_T=c_F-c_g\le R_{\rm mf}-c_g,
\qquad
R_{\rm mf}=\frac{\operatorname{rank}A_1}{2}.
\tag{1}
\]

For a number-conserving ansatz one can instead use

\[
R_{\rm mf}=\operatorname{rank}T_1\le n_c.
\tag{2}
\]

The rank is taken across one cell boundary, as in the supplied theorem. It need not be an integer in the Majorana convention.

Proof. The commuting Sugawara and coset stress tensors have additive central terms. Assumption A3 identifies the remaining stress tensor with the physical one, so `c_T=c_F-c_g`. Apply the supplied free bound `c_F<=R_mf`. Unitarity also gives `c_F>=c_g`. If the right side of (1) is negative, the assumptions are incompatible. One cannot repair an invalid embedding by replacing a negative bound by zero. This proves (1) and (2).

Here A2, once its current algebra is specified, is a mathematical coset construction. A3, and the claim that a chosen lattice gauging realizes A2, are the substantive dynamical assumptions.

There is a variant with `c_T<=c_F-c_g` if a further flow from the coset is established and the hypotheses of the unitary relativistic c-theorem hold. That flow is another assumption. A general collection of sectors with distinct velocities does not automatically admit one such application. For distinct velocities, apply the theorem to conformal blocks with compatible gauge currents and add their charges. Gauge currents mixing unequal-velocity modes require a separate argument for the proposed fixed point.

### What IR equivalence must mean

It must include the scaling operator algebra, the stress-tensor Ward identities and central terms, and the matching of physical gauge-invariant lattice operators to continuum operators. The finite-size low-energy spectrum, with each velocity identified, should realize the corresponding conformal sectors and multiplicities. Boundary conditions and possible gauge holonomy sectors must be treated consistently.

Equality of global symmetries, a shared anomaly, a similar mean-field dispersion, or good variational energy does not suffice. Matching a Casimir coefficient alone only matches `sum_a v_a c_a`, not `sum_a c_a`. Matching one wavefunction does not specify a dynamical exponent or a parent Hamiltonian. Matching its entanglement coefficient can identify the CFT central charge only after conformal ground-state scaling is established.

Finite-group identifications and global quotients can change sectors and operator content without changing the stress-tensor central term. They still matter for the operator matching in A3.

## 2. The minimal representation

### Corollary 2. Minimal, number-conserving partons

Take `d>=2` complex orbitals with `N_x=1`. Suppose their number-conserving nearest-neighbor mean field has an active diagonal U(1) current and satisfies A1-A3 with this charge sector removed. Then

\[
c_T=c_F-1\le\operatorname{rank}T_1-1\le d-1.
\tag{3}
\]

The first inequality follows from Theorem 1. The subtraction is 1 regardless of the positive U(1) level. If no charge current is present, this proof does not apply. In particular, a paired ansatz that Higgses the U(1) to a discrete group does not justify subtracting 1 merely because the original code was defined by `N_x=1`.

For `d=1`, the physical chain has one state and `c_T=0`. A uniform one-parton-per-site band is full and has no linear crossings. Subtracting a nonexistent charge current from that band would be incorrect.

### Saturation in the continuum construction

For `d>=2`, choose identical hopping `t>0` for all `d` flavors, each at filling `1/d`. There are Fermi points at `+/- pi/d` and nonzero velocity

\[
v_F=2t\sin(\pi/d).
\]

The hopping matrix has rank `d`; all `d` complex channels are gapless. The continuum theory has `c_F=d`. With the convention

\[
j_a(z)j_b(0)\sim\frac{\delta_{ab}}{z^2},
\qquad J=\sum_{a=1}^d j_a,
\]

the diagonal current has level `d` and

\[
T_{U(1)}=\frac{:JJ:}{2d},\qquad c_{U(1)}=1.
\]

The traceless bilinears generate `SU(d)_1` by Wick contraction. They commute with `J`. Their Sugawara central charge is

\[
c_{SU(d)_1}=\frac{d^2-1}{d+1}=d-1.
\]

The two central charges add to `d`, so the remaining stress tensor has zero norm in the unitary free theory. This proves the conformal stress-tensor decomposition

\[
U(d)_1\supset SU(d)_1\times U(1)_d,
\qquad U(d)_1/U(1)_d\simeq SU(d)_1
\tag{4}
\]

at the level of local stress tensors and central charge. Finite center identifications are understood in the full operator spectrum.

Thus the uniform sea saturates the conditional continuum bound. The supplied ULS result gives an actual nearest-neighbor physical realization with `c=d-1`.

There is a qualification to the request for equality for the projected sea. The on-site projection of that sea is a particular variational wavefunction. Equation (4) alone does not prove its infrared correlations, and it does not show that it is the exact ULS ground state. Identifying its critical limit with `SU(d)_1` is standard parton input, not proved here. The unconditional saturation statement for a nearest-neighbor chain comes from the ULS result supplied in the question. These two statements should not be conflated.

## 3. Non-Abelian families with a proved dimension inequality

The representation dimensions and conformal-embedding arithmetic below can be proved. Identifying a projected lattice state or a specified interacting Hamiltonian with each coset still invokes A3.

### 3.1 Rectangular SU(N) representations and color gauging

Use complex fermions `f_{a alpha}`, with `a=1,...,N`, `alpha=1,...,m`. There are `n_c=Nm` orbitals. Gauge `U(m)` on the color index and select the character `det^q`, with `1<=q<=N-1`. This fixes total occupation `qm` and imposes the color SU(m) singlet condition.

The exterior-algebra decomposition is

\[
\wedge^{qm}(\mathbb C^N\otimes\mathbb C^m)
=\bigoplus_{\lambda\vdash qm}
S_\lambda\mathbb C^N\otimes S_{\lambda^t}\mathbb C^m.
\]

The color character `det^q` has Young diagram `(q^m)`, so the physical representation has diagram `(m^q)`, a rectangle with `q` rows and `m` columns. This also explains the local representation rather than merely matching its central charge.

There are `m` copies of an SU(N) fundamental current and `N` copies of an SU(m) fundamental current. Their induced levels are `m` and `N`, respectively. The total-number current has level `Nm`. Hence

\[
U(Nm)_1\supset SU(N)_m\times SU(m)_N\times U(1)_{Nm},
\tag{5}
\]

\[
c_g=1+\frac{N(m^2-1)}{N+m},
\qquad
c_{\rm coset}=\frac{m(N^2-1)}{N+m}.
\tag{6}
\]

For `m=1`, the SU(m) factor is absent. The two commuting non-Abelian currents and the charge current have total central charge

\[
\frac{m(N^2-1)+N(m^2-1)}{N+m}+1=Nm.
\]

Unitarity again proves the stress-tensor equality underlying (5). This is a level and embedding calculation, independent of an assumed lattice IR.

The local dimension is the Weyl product

\[
D_q(N,m)=\prod_{i=1}^q\prod_{j=q+1}^N
\frac{m+j-i}{j-i}.
\tag{7}
\]

We now prove `c_coset<=D_q-1`, uniformly over this infinite family. Taking adjacent ratios in (7) gives

\[
\frac{D_{q+1}}{D_q}
=\frac{\binom{m+N-q-1}{N-q-1}}
{\binom{m+q}{q}}.
\]

For fixed positive `m`, `binom(m+s,s)` increases with `s`. The dimensions therefore increase toward the middle and decrease symmetrically afterward. Their minimum over `1<=q<=N-1` occurs at the endpoints:

\[
D_q\ge D_1=\binom{N+m-1}{m}.
\]

The latter counts nonnegative occupations of `m` objects in `N` boxes. It includes the configuration with all objects in box 1 and the `m(N-1)` distinct configurations obtained by placing `r=1,...,m` objects in one other box and leaving the rest in box 1. Consequently

\[
D_q-1\ge m(N-1),
\qquad
\frac{m(N^2-1)}{N+m}\le m(N-1),
\tag{8}
\]

where the last step uses `m>=1`. This proves the bound. For `m>1` the second inequality is strict. For `m=1`, saturation occurs at the fundamental or its conjugate, including coincidences at small `N`.

The cases `q=0,N` are excluded for a reason. Their local representations are one-dimensional, and the corresponding uniform occupation is empty or full. They do not have the `Nm` massless channels assumed in (5).

For `m=1`, the physical representation is `wedge^q C^N`, of dimension `binom(N,q)`, and the coset has `c=N-1<=binom(N,q)-1`. This includes the minimal representation at `q=1`.

### 3.2 Spin S and Takhtajan-Babujian

Set `N=2`, `m=2S`, `q=1` above. The local representation is `Sym^m C^2`, so

\[
d=m+1=2S+1,
\qquad n_c=2m=4S,
\qquad
c=\frac{3m}{m+2}=\frac{3S}{S+1}\le m=d-1.
\tag{9}
\]

The gauge algebra is `SU(m)_2` plus U(1). Its central charge is `1+2(m^2-1)/(m+2)`. The inequality follows from `m+2>=3` for integer `m>=1`. It is saturated only for `S=1/2`. At `S=1`, the TB value is `3/2<2`; at `S=3/2`, it is `9/5<3`. The large-S limit is `c -> 3`, while `d-1=2S` grows without bound.

The usual identification of the TB chain with `SU(2)_{2S}` is an IR input for this example. The representation, induced levels, and inequality in (9) have been derived here.

There is an equivalent central-charge construction with compact symplectic gauge group `Sp(2m)`, meaning the rank-m group with a `2m`-dimensional fundamental. Take `2m` complex fermions transforming in that fundamental. In their exterior algebra, the invariant states are

\[
|0\rangle,\quad \omega^\dagger|0\rangle,\quad\ldots,
\quad(\omega^\dagger)^m|0\rangle,
\]

where `omega^dagger=sum_j f^dagger_{2j-1}f^dagger_{2j}`. The maximal torus requires occupied indices to occur in symplectic pairs. Symplectic rotations mixing pairs then fix the invariant at each even degree up to a scalar. There are no odd-degree invariants. Thus the local invariant space has dimension `m+1`.

The pair operators and `(N-m)/2` act on these states as spin `m/2`. In the free CFT their currents have SU(2) level `m`, obtained by adding the level-one currents of the `m` pairs. The commuting symplectic bilinears have level 1. Consequently

\[
SO(4m)_1\supset Sp(2m)_1\times SU(2)_m,
\qquad
\frac{m(2m+1)}{m+2}+\frac{3m}{m+2}=2m.
\tag{10}
\]

This gives the same bound with a different code and gauge group. For `m=1`, four Majoranas with SU(2) gauge redundancy leave a physical two-state spin and a central charge of 1.

### 3.3 An SO(n) family with d=2

Take `n` complex fermions in the vector of SO(n), for even `n>=4`. The only invariant exterior forms are the vacuum and the volume form. To see this, use the SO(n) transformations that reverse the signs of two coordinate axes. For any nonempty proper occupied subset, choose one axis inside the subset and one outside. The corresponding exterior basis vector changes sign. Since these transformations act diagonally, its coefficient in an invariant vector must vanish. The local code therefore has two states, both of even fermion parity. Each complex vector is two real vectors, so the gauge algebra is `SO(n)_2` inside `SO(2n)_1`, with

\[
c_g=\frac{2[n(n-1)/2]}{2+n-2}=n-1,
\qquad c_{\rm coset}=n-(n-1)=1=d-1.
\tag{11}
\]

For `n=4`, interpret SO(4) as its two simple SU(2) factors; both have level 2 and their charges add to 3. Equation (11) is unchanged. This is an example in which many orbitals cause no difficulty because the gauge subtraction grows with them.

### 3.4 Product codes

If independent codes have dimensions `d_i` and coset charges `c_i<=d_i-1`, their product has `d=prod_i d_i` and

\[
\sum_i c_i\le\sum_i(d_i-1)\le\prod_i d_i-1.
\]

For two factors the last inequality is `(d_1-1)(d_2-1)>=0`; induction proves the general statement. A direct sum of local codes is different and does not justify adding their central charges.

## 4. Why an enumeration of all gauge representations cannot prove the claim

### Proposition 3. A two-state code with an unbounded formal one-U(1) coset

Let `n>=4` be even. Introduce `n` complex fermions and the integer U(1) charges

\[
q=(1,2,4,\ldots,2^{n-2},2^{n-1}-1),
\qquad Q=2^{n-1}-1.
\]

Define the local code by

\[
\sum_{a=1}^n q_a n_a=Q.
\tag{12}
\]

It has exactly two basis states,

\[
|A\rangle=|1,1,\ldots,1,0\rangle,
\qquad
|B\rangle=|0,0,\ldots,0,1\rangle.
\tag{13}
\]

Proof. If the last orbital is occupied, all others must be empty because all charges are positive. If the last orbital is empty, reaching `Q`, which is the sum of all remaining charges, requires every remaining orbital to be occupied. Both states have odd parity when `n` is even. Their relative transition operator is even, so the local code describes an ordinary two-state tensor factor without a varying physical fermion parity.

Take independent uniform nearest-neighbor hopping for all `n` species at half filling. Every species has two linear Fermi points and the same velocity. Thus

\[
R_{\rm mf}=n,\qquad c_F=n.
\]

The mean gauge density is `sum_a q_a/2=Q`, so there is no mismatch even at the level of average charge. The U(1) current has positive level

\[
\kappa=\sum_aq_a^2
=\frac{4^{n-1}-1}{3}+(2^{n-1}-1)^2.
\]

The free conformal embedding into this U(1) and its commutant exists. Its formal quotient has

\[
c\bigl(U(n)_1/U(1)_\kappa\bigr)=n-1>1=d-1.
\tag{14}
\]

For the smallest example, `q=(1,2,4,7)`, `Q=7`, `d=2`, `kappa=70`, and the formal charge is 3. This is a continuous, faithful, vector U(1) gauge representation. The obstruction is not restricted to finite-group gauging or to anomalous chiral constructions.

The projected half-filled sea is not identically zero. On an even ring, choose a filled set of consecutive single-particle momenta for each species, with a twist if needed to avoid a Fermi-level degeneracy. A physical configuration with `L/2` sites of type A and `L/2` sites of type B gives a product of nonzero Vandermonde determinants for the occupied positions of each flavor. Thus the failure cannot be dismissed as a zero-norm ansatz.

### An operator-level diagnosis of the failure

On the code (13), the exact local identity

\[
2n_1-3n_2+n_3=0
\tag{15}
\]

holds. In the proposed free coset, however, the current

\[
J_a=2j_1-3j_2+j_3
\]

has zero mixed OPE with the gauged current because `2*1-3*2+1*4=0`, and has positive self-OPE

\[
J_a(z)J_a(0)\sim\frac{14}{z^2}.
\tag{16}
\]

It therefore survives the one-U(1) coset. Equation (15) requires its proposed physical density to vanish, including its smooth scaling part. Equations (15) and (16) cannot both hold under the usual identification of parton densities with physical densities.

This proves failure of A3 with that identification. It does not prove that no unrelated physical chain could share the same abstract CFT. It also does not compute the actual IR of the projected sea. Those would be stronger claims.

The exact code has additional continuous redundancies. A full set of independent linear constraints is

\[
n_j-n_1=0\quad (j=2,\ldots,n-1),
\qquad n_1+n_n=1.
\tag{17}
\]

Their charge matrix has rank `n-1`. With all free species active, its current level matrix is the positive-definite Gram matrix `K=CC^t`. Gauging this entire torus removes `n-1` units of central charge and leaves 1. The original weighted charge lies in this constraint space.

Thus two choices of gauge generators can define the same finite-dimensional on-site code while giving different naive continuum quotients. The lattice hard constraint knows the binary occupation spectrum. Linearizing one charge-density constraint before projection loses some of its consequences. This is the precise noncommutation between taking a continuum limit and enforcing the hard constraint in this example.

It follows that the statement

\[
\dim(P\mathcal F)=d\quad\Longrightarrow\quad
n_c-c(G_g)\le d-1
\]

is false for arbitrary choices of `G_g`. Finite lists of familiar representations cannot fix this. Equation (14) supplies infinitely many counterexamples to that representation-only implication.

### A finite gauge group gives an even simpler warning

Two complex orbitals with the local even-parity constraint have code `span{|00>,|11>}`, so `d=2`. A finite Z2 orbifold of the free `c=2` CFT retains central charge 2, since it retains the stress tensor and adds twisted sectors rather than a positive-central-charge Sugawara gauge sector.

This is not a qubit nearest-neighbor counterexample. A lattice gauge theory has link states, while the strict on-site parity code has additional structure. In this example the spin SU(2) bilinears annihilate both code states. Completing to that continuous redundancy gives an SU(2)_1 subtraction and leaves 1. Finite-group orbifolding and hard parity projection are distinct operations.

## 5. A positive theorem for completed Abelian constraints

### Theorem 4. Occupation-code affine-dimension bound

Let `S={s^(1),...,s^(d)} subset {0,1}^p` be the allowed occupation configurations of a code. Remove any frozen orbitals and let `r` be the affine dimension of `S` in the remaining real occupation space. Suppose the following conditions hold.

(B1) The ansatz has one independent massless complex channel per remaining orbital, so `c_F=p`, with commuting number currents of positive level matrix.

(B2) All linear number constraints constant on `S` are gauged, including their backgrounds. Thus the gauge-charge row space is the orthogonal complement of `span{s^(i)-s^(1)}`.

(B3) The physical IR is the associated torus coset, as in A3.

Then

\[
c_T=r\le d-1.
\tag{18}
\]

Proof. The differences span a space of rank `r<=d-1`, since there are only `d-1` nontrivial differences. The independent linear constraints have rank `p-r`. Because the vertices are integral, this space admits an integer basis and defines a compact torus. Their current level matrix is positive definite by B1. Every independent U(1) current contributes 1 to the Sugawara charge, regardless of level or mixing. Hence

\[
c_T=p-(p-r)=r\le d-1.
\]

This theorem also gives an interpretation: the surviving bosons count independent occupation-density directions in the affine span of the local configurations. The minimal one-particle code consists of `d` unit vectors, whose affine dimension is `d-1`. The weighted two-state code has affine dimension 1.

The theorem has limits. Completing only linear constraints need not reproduce an arbitrary discrete allowed set, because additional nonlinear constraints may remain. Pairing or flavor mixing may destroy the assumed number currents. Non-Abelian singlets can be coherent superpositions of occupation configurations, so their full dimension is not counted by this argument. Most importantly, B3 remains a dynamical assumption. The theorem proves a useful completed-Abelian subclass, not a universal representation theorem.

## 6. Checks on the requested examples

| Physical chain or construction | Local d | Gapless theory and c | Parton or gauge interpretation | Status of comparison |
|---|---:|---|---|---|
| SU(d) ULS | d | SU(d)_1, d-1 | d complex partons; U(1)_d removes 1 | Saturation is supplied; projected-sea matching is separate |
| Spin-S TB | 2S+1 | SU(2)_{2S}, 3S/(S+1) | 4S complex partons; SU(2S)_2 plus U(1), or Sp(4S)_1 | Derived dimension inequality; standard TB IR input |
| Ising | 2 | Majorana, 1/2 | Jordan-Wigner free Majorana description | 1/2 < 1; no gauge subtraction is needed |
| Three-state Potts | 3 | Z3 parafermion, 4/5 | SU(2)_3/U(1) continuum coset | 4/5 < 2; not an on-site parton proof |
| Four-state Potts | 4 | c=1 | SU(2)_4/U(1) has the same charge | 1 < 3; matching full operator content requires more than c |
| Gapless t-J regime | 3 | spin and charge sectors, total 2 | A constrained electron description is not quadratic free hopping | Saturating example supplied; not every parameter choice has c=2 |
| Hubbard with two gapless modes | 4 | spin and charge sectors, total 2 | Two physical complex fermions | 2 < 3; interaction and filling determine which modes remain |
| SU(4)_1 on a spinful fermion site | 4 | SU(4)_1, 3 | Four minimal auxiliary partons, U(1)_4 removed | Saturates 3 and exceeds the two-orbital free bound |

For the Ising chain, a standard Jordan-Wigner representation has two Majoranas per site and a rank-one intersite Majorana block at criticality. The supplied free theorem then gives `c<=1/2`, with the single critical branch saturating it. The string in the spin-to-fermion map does not change this spectral statement, subject to the familiar parity-dependent boundary sectors.

The continuum parafermion calculation is

\[
\frac{SU(2)_N}{U(1)}:\qquad
c=\frac{3N}{N+2}-1=\frac{2(N-1)}{N+2}\le N-1.
\tag{19}
\]

For a physical N-state clock realization, `d=N`, so (19) is below the conjectured bound. At `N=2,3,4` it gives `1/2,4/5,1`. The critical Ising and three-state Potts identifications are standard example input, not derived here. The ordinary ferromagnetic N-state Potts transition is not generally a conformal example for `N>4`; one cannot label every N-state model by (19). Also, discrete Z_N gauging by itself does not produce the U(1) subtraction in (19). A particular continuum coset and a literal on-site code must be distinguished.

For t-J, the three states are empty, up, and down. Enforcing no double occupancy is a genuine interaction constraint, so applying the two-orbital free rank bound to that constrained Hamiltonian is invalid. The supplied gapless `c=2` example equals `d-1`. Its equality of central charge with `SU(3)_1` does not establish equality of CFTs: velocities, compactification data, and symmetries can differ.

For Hubbard, two gapless charge and spin sectors give `c=2`. At half filling with repulsive interaction the charge sector can be gapped, giving `c=1`; other regimes likewise require care. The stated two-mode value is an example, not a bound derived for all interacting Hubbard variants from the free theorem.

For the spinful-site SU(4) realization, label the four states `0, up, down, updown` by `a=1,...,4`, and write the tensor-product swap

\[
P_{x,x+1}=\sum_{a,b=1}^4 X_x^{ab}X_{x+1}^{ba},
\qquad X^{ab}=|a\rangle\langle b|.
\]

It is a nearest-neighbor operator on two four-dimensional sites. In fermionic coordinates the appropriate local parity signs must be included; the full swap conserves total fermion parity and charge. It contains interactions beyond a quadratic hopping Hamiltonian. Its physical two-orbital Fock description has free bound 2, whereas its four-auxiliary-parton description has coset charge `4-1=3`. There is no contradiction.

Finally, the induced-level criterion supplied in the background is another conditional subclass. If the physical on-site representation R embeds `G_ind(R)` in `SU(d)_1`, coset unitarity gives `c(G_ind(R))<=d-1`. Since

\[
\frac{d}{dk}\frac{k\dim G}{k+h^\vee}
=\frac{h^\vee\dim G}{(k+h^\vee)^2}>0,
\]

the same follows for `k<=ind(R)`. An arbitrary parton construction does not prove that inequality on the physical level. For spin 1, for example, the minimal SU(3)_1 construction induces SU(2)_4, with charge 2, whereas TB uses SU(2)_2, with charge 3/2. These are distinct consistent levels for the same local dimension.

## 7. What survives of the one-bond argument

Equation (1) is the surviving statement:

\[
\mathcal B(H_{\rm mf},\mathcal A_g)
=R_{\rm mf}-c_g,
\qquad c_T\le\mathcal B
\]

whenever A1-A3 hold. With no gauge sector it reduces to the supplied free bound. If the mean field saturates its rank bound, it equals the coset central charge. In the minimal charge-conserving construction it is `rank T_1-1`.

This is not an invariant of a physical bond. It depends on the auxiliary representation, the ansatz, which channels are massless, the embedding levels, and the dynamical identification of the quotient. A gauge generator acting only on massive fields has no massless Sugawara charge to subtract. A finite gauge group has no such positive-charge current sector. Spectator orbitals can change auxiliary hopping rank without changing the physical code or its IR.

Projection does not preserve a matrix rank relevant to the supplied determinant proof. In the minimal code,

\[
P f_{x,a}^\dagger f_{y,b}P=0\quad(x\ne y),
\tag{20}
\]

because the hop changes occupations `1,1` to `2,0` or annihilates the state. Thus `PH_mf P` has no intersite hopping even when the projected variational state is critical. Physical exchange can arise through virtual departures from the code, or can be specified directly, and is quartic or higher in partons. Equation (20) rules out identifying the projected hopping rank with the rank of the actual interacting bond.

The gauge subtraction is a stress-tensor subtraction across the entire low-energy theory. It is not a count of independent rows removed from one microscopic hopping matrix. The original Laurent-determinant argument has no general interacting determinant to which it can be applied after projection.

## 8. The missing lemmas for general nearest-neighbor chains

There are two separate gaps.

### Gap I. From a projected ansatz or lattice gauge theory to its proposed coset

One must show that hard local constraints and the continuum limit yield the asserted stress tensor and operator spectrum, with no missing constraints and no additional gapless sectors. The weighted family proves that merely specifying a gauge group and subtracting its affine charge is insufficient. The density-current contradiction (15)-(16) locates the failure without assuming the c-d conjecture.

Introducing dynamical gauge links does not solve this automatically. Gauss's law has the schematic form

\[
E_{x,x+1}-E_{x-1,x}
=\sum_a q_a n_{x,a}-Q.
\]

Matter outside the on-site code can be accompanied by electric flux. The physical Hilbert space is then not a tensor product of the specified `C^d` code sites. Conversely, imposing the strict code without links removes bare hopping as in (20). A strong-coupling limit that suppresses flux and generates exchange must be analyzed before borrowing a weak-coupling continuum coset. These limits need not agree.

### Gap II. Coverage and a dimension-controlled auxiliary bound

A sufficient missing theorem would be the following.

For every unitary conformal nearest-neighbor chain with local dimension d, there exists a nearest-neighbor free-parton ansatz and a complete active gauge-current algebra such that (i) the physical IR is its coset, or is bounded by a proved admissible flow from that coset, and (ii)

\[
R_{\rm mf}-c_g\le d-1.
\tag{21}
\]

The minimal version would use `d` number-conserving partons, an active U(1), and a controlled coset or c-decreasing flow. It would immediately prove the conjecture by (3). Neither existence nor control of that flow follows from the kinematic identity `C^d = {N_f=1}`.

For more general representations, one also needs a justified completion of all constraints and a dimension inequality for the resulting active coset. Theorem 4 proves such an inequality for its Abelian subclass, and (8)-(11) prove it for several non-Abelian families. Proposition 3 shows that it is false without qualification. No theorem covering all complete non-Abelian parton codes has been proved here.

A variational principle bounds the ground-state energy, not the coefficient of a logarithmic entropy or the number of low-energy conformal channels. A low variational energy density permits errors on scales relevant to a critical IR. It does not imply a central-charge ordering.

A gapped adiabatic-continuation theorem does not apply to a thermodynamic critical chain with gap of order `1/L`. A path through critical Hamiltonians can change its gapless content at multicritical points. Claiming that a path stays in one universality class would assume the needed conclusion. Likewise, Gutzwiller projection is noninvertible local filtering, not a controlled gapped adiabatic transformation or a known c-monotone operation on arbitrary lattice states.

For the very narrow class of exact charge-only minimal cosets, there are immediate examples outside the class. A uniform fully active d=2 construction gives `SU(2)_1` with c=1, whereas the Ising chain has c=1/2. Ising is reached by a different free construction, so this is not a counterexample to all parton methods. More generally, a generic-radius compact-boson theory need not be the fixed SU(2)_1 compactification obtained by the symmetric two-flavor charge-only quotient; changing the radius requires further dynamics or a different construction.

No class of CFTs has been proved here unreachable by every nearest-neighbor parton construction with arbitrary representations, pairing, and additional controlled interactions. The scope is too broad for that conclusion. The failure established here is more concrete: on-site representation data do not determine the correct continuum gauging, and there is no established universality theorem placing arbitrary nearest-neighbor chains in the class of Theorem 1 with (21).

## 9. Verification and conclusion

Run `python3 f1.py`. It reads no input files. It exhaustively enumerates the weighted-code sectors for `n_c=4,6,8,10,12`, verifies their common local parity and compatible half filling, checks the killed current and its positive coset level, and checks that single-parton hopping vanishes after projection. It also verifies the rectangular dimensions and induced central charges over 450 cases, together with the spin, symplectic, orthogonal, exterior-power, parafermion, and requested-example arithmetic. Independent checks construct the fermion-signed Fock-space generators and compute their common kernels by rational row reduction: the color and symplectic singlet spaces have dimension `m+1` for `1<=m<=5`, and the SO(n) singlet spaces have dimension 2 for `n=4,6`. All checks passed. The output is saved in `f1_checks.txt`. The proofs above extend the inequalities beyond these finite checks.

The precise outcome is (i) a proved conditional coset bond-rank bound, (ii) a minimal-representation `d-1` bound and several infinite families satisfying it, and (iii) a representation-only counterexample that identifies the extra premise needed for hard projection. No new violation of the nearest-neighbor c-d conjecture, and no proof for arbitrary interacting chains, is claimed.
