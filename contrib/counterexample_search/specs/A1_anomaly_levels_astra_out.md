| On-site representation | Invariant bonds, including \(1\) | Projective class | \(\operatorname{ind}(R)\) | Levels compatible with a center action of \(T\) in the diagonal WZW theory | Smallest pure WZW theory above \(d-1\) |
|---|---:|---|---:|---|---|
| (a) \(SU(2)\), spin \(1\) | \(3\) | Trivial | \(4\) | All \(k\), choosing \(T=1\) | \(SU(2)_5\), \(c=15/7\) |
| (b) \(SU(2)\), spin \(3/2\) | \(4\) | Order \(2\) | \(10\) | \(k\) odd | None |
| (c) \(Sp(4)\), \(\mathbf4\) | \(3\) | Order \(2\) | \(1\) | None | None |
| (d) \(SU(3)\), \(\mathbf3\) | \(2\) | Order \(3\) | \(1\) | \(3\nmid k\) | \(SU(3)_2\), \(c=16/5\) |
| (e) \(SU(2)^2\), \((\frac12,\frac12)\) | \(4\) | \((1,1)\), order \(2\) | \((2,2)\) | Both levels odd | Equal levels: \((3,3)\), \(c=18/5\). Unequal levels: \((1,5)\), \(c=22/7\) |
| (f) \(SO(N)\) vector, \(N\) odd | \(3\) | Trivial | \(2\) | All \(k\) | \(SO(N)_3\), \(c=\frac{3N(N-1)}{2(N+1)}\) |
| (f) \(SO(N)\) vector, \(N\equiv0\pmod4\) | \(3\) | Order \(2\) | \(2\) | \(k\) odd | Same |
| (f) \(SO(N)\) vector, \(N\equiv2\pmod4\) | \(3\) | Order \(2\) | \(2\) | \(4\nmid k\) | Same |
| (g) \(G_2\), \(\mathbf7\) | \(4\) | Trivial | \(2\) | All \(k\) | \(G_{2,4}\), \(c=7\) |

For the orthogonal rows, “center” means the axial center symmetry of the simply connected, diagonal affine theory. This distinction matters. The symmetry used at odd \(k\) for even \(N\) is the element \(-1\in Spin(N)\), which maps to the identity in \(SO(N)\). A literal quotient-group WZW model requires a separate choice of global form and modular invariant.

There are also two distinctions needed for the requested minima.

First, containing \(\widehat{\mathfrak g}_k\) does not imply that the full CFT is a tensor product with the diagonal \(G_k\) theory. Extensions can introduce new internal symmetries and change anomaly matching. This is essential for the parton construction.

Second, an absolute minimum over arbitrary additional unitary sectors is a broader CFT classification problem. Below I give exact minima for pure WZW theories and for products with one unitary Virasoro minimal model, together with explicit broader constructions. For \(Sp(4)\), a stronger result is possible: among unitary CFTs containing \(Sp(4)\) currents, the smallest total central charge strictly above \(3\) that can match the spinor anomaly is
\[
\boxed{c=\frac{16}{5}}.
\]
It is a local extension of \(Sp(4)_1\) times tricritical Ising, not their unextended tensor product.

All calculations below are closed-book. I ran short Python checks of the rational arithmetic, dimensions, indices, and finite-group congruences. No files or external sources were consulted.

---

The twisted-sector calculation can be organized uniformly.

Let \(\widetilde G\) be the simply connected cover, and write the center character of an irreducible representation as
\[
q(\lambda)\in P/Q.
\]
The on-site projective class for the adjoint group is
\[
r=q(R).
\]

For the diagonal invariant, the physical diagonal-group representation in a center-twisted sector is, with one orientation convention,
\[
\lambda\otimes(J_z\lambda)^*.
\]
Its center character is
\[
q(\lambda)-q(J_z\lambda).
\]
A simple current has an additive center character, so
\[
q(J_z\lambda)=q(J_z)+q(\lambda).
\]
Consequently every state in this twisted sector has the same projective class,
\[
-q(J_z).
\]
Reversing the twist replaces this by \(q(J_z)\). I will use the latter convention. Thus the matching equation is
\[
\boxed{q(J_z)=r.} \tag{1}
\]

For the WZW simple currents used here,
\[
J_z=k\omega_z,
\qquad
q(J_z)=k[\omega_z].
\]
The problem is therefore multiplication by \(k\) in a finite abelian group.

A trivial on-site class requires a trivial IR mixed anomaly. It does not require a nontrivial center action of translation. Choosing \(T=1\) is allowed.

---

For the bond count, if
\[
R\otimes R=\bigoplus_\alpha V_\alpha\otimes\mathbb C^{m_\alpha},
\]
then
\[
\operatorname{End}_G(R\otimes R)
=\bigoplus_\alpha \operatorname{Mat}_{m_\alpha}(\mathbb C),
\qquad
\dim_{\mathbb C}\operatorname{End}_G(R\otimes R)=\sum_\alpha m_\alpha^2.
\]
The real dimension of its Hermitian part is the same. All decompositions in this problem are multiplicity-free.

| Case | \(\operatorname{Sym}^2R\) | \(\Lambda^2R\) | Dimension check |
|---|---|---|---|
| (a) Spin \(1\) | Spin \(0\oplus2\) | Spin \(1\) | \(1+5=6,\quad3=3\) |
| (b) Spin \(3/2\) | Spin \(1\oplus3\) | Spin \(0\oplus2\) | \(3+7=10,\quad1+5=6\) |
| (c) \(Sp(4)\), \(\mathbf4\) | \(\mathbf{10}\) | \(\mathbf1\oplus\mathbf5\) | \(10,\quad1+5=6\) |
| (d) \(SU(3)\), \(\mathbf3\) | \(\mathbf6\) | \(\overline{\mathbf3}\) | \(6,\quad3\) |
| (e) \((\frac12,\frac12)\) | \((0,0)\oplus(1,1)\) | \((1,0)\oplus(0,1)\) | \(1+9=10,\quad3+3=6\) |
| (f) \(SO(N)\), vector | \(\mathbf1\oplus S^2_0V\) | \(\mathfrak{so}(N)\) | \(1+\left[\frac{N(N+1)}2-1\right],\quad\frac{N(N-1)}2\) |
| (g) \(G_2\), \(\mathbf7\) | \(\mathbf1\oplus\mathbf{27}\) | \(\mathbf7\oplus\mathbf{14}\) | \(1+27=28,\quad7+14=21\) |

Here are the mechanisms behind these decompositions.

For \(SU(2)\),
\[
j\otimes j=\bigoplus_{J=0}^{2j}J.
\]
The exchange parity of the total-spin-\(J\) summand is
\[
(-1)^{2j-J}.
\]
One way to see this is to construct that summand with \(2j-J\) antisymmetric spinor contractions. Each contraction changes sign under exchange. Substituting \(j=1\) and \(j=3/2\) gives the first two rows.

For \(Sp(4)\), the invariant symplectic form supplies the singlet in \(\Lambda^2\mathbf4\). Removing it leaves the \(\mathbf5\). The symmetric square is the adjoint \(\mathbf{10}\).

For \(SU(3)\), contraction with the invariant epsilon tensor identifies \(\Lambda^2\mathbf3\) with \(\overline{\mathbf3}\).

For the product group, each spin-\(\frac12\) singlet is antisymmetric and each triplet is symmetric. The exchange signs multiply. Thus singlet-singlet and triplet-triplet are symmetric.

For \(SO(N)\), symmetric tensors split into their trace and traceless part; antisymmetric tensors are the adjoint. In particular,
\[
\begin{aligned}
\mathbf5\otimes\mathbf5&=\mathbf1\oplus\mathbf{10}\oplus\mathbf{14},\\
\mathbf6\otimes\mathbf6&=\mathbf1\oplus\mathbf{15}\oplus\mathbf{20}',\\
\mathbf7\otimes\mathbf7&=\mathbf1\oplus\mathbf{21}\oplus\mathbf{27}.
\end{aligned}
\]
The general orthogonal statement here assumes \(N\ge5\). At \(N=4\), the antisymmetric square splits into self-dual and anti-self-dual parts, giving case (e).

For \(G_2\), the invariant cross product gives a surjection
\[
\Lambda^2\mathbf7\longrightarrow\mathbf7.
\]
Its complementary \(\mathbf{14}\) is the adjoint. The symmetric traceless tensors form the \(\mathbf{27}\).

As independent checks, the relevant Weyl dimension formulas are
\[
\dim_{Sp(4)}(a,b)
=\frac{(a+1)(b+1)(a+b+2)(a+2b+3)}6
\]
and
\[
\dim_{G_2}(a,b)
=\frac{(a+1)(b+1)(a+b+2)(a+2b+3)(a+3b+4)(2a+3b+5)}{120}.
\]
They give
\[
\dim_{Sp(4)}(1,0)=4,\quad
\dim_{Sp(4)}(0,1)=5,\quad
\dim_{Sp(4)}(2,0)=10,
\]
and
\[
\dim_{G_2}(1,0)=7,\quad
\dim_{G_2}(0,1)=14,\quad
\dim_{G_2}(2,0)=27.
\]
The script also verified every row against
\[
\dim\operatorname{Sym}^2R=\frac{d(d+1)}2,
\qquad
\dim\Lambda^2R=\frac{d(d-1)}2.
\]

---

The embedding indices can be computed directly from a long-root \(SU(2)\).

Normalize its Cartan generator \(H\) to have eigenvalues \(1,-1\) on a doublet. Since a long-root \(SU(2)\) has embedding index one,
\[
\operatorname{ind}(R)=\frac{\operatorname{tr}_R H^2}{2}. \tag{2}
\]
For spin \(j\), the eigenvalues are \(2m\), so
\[
\begin{aligned}
\operatorname{ind}(j)
&=\frac12\sum_{m=-j}^j(2m)^2\\
&=2\sum_{m=-j}^j m^2\\
&=\frac23j(j+1)(2j+1).
\end{aligned}
\]
Thus
\[
\operatorname{ind}(1)=\frac23(1)(2)(3)=4,
\]
and
\[
\operatorname{ind}\left(\frac32\right)
=\frac23\frac32\frac52(4)=10.
\]

For the remaining representations, restriction to a long-root \(SU(2)\) gives

\[
\begin{array}{c|c|c}
R & R|_{SU(2)} & \operatorname{ind}(R)\\ \hline
Sp(4):\mathbf4 & \mathbf2+2\mathbf1 &1\\
SU(3):\mathbf3 & \mathbf2+\mathbf1 &1\\
SU(2)^2:(\mathbf2,\mathbf2),\text{ each factor}&2\mathbf2&2\\
SO(N):V&2\mathbf2+(N-4)\mathbf1&2\\
G_2:\mathbf7&2\mathbf2+3\mathbf1&2
\end{array}
\]

Each doublet contributes \(1^2+(-1)^2=2\) to the trace in (2).

An independent calculation uses
\[
\operatorname{ind}(R)
=\frac{2\dim R\,C_2(R)}{\dim G}. \tag{3}
\]
The Casimirs are
\[
C_2(j)=j(j+1),\quad
C_2^{Sp(4)}(\mathbf4)=\frac54,\quad
C_2^{SU(3)}(\mathbf3)=\frac43,
\]
\[
C_2^{SO(N)}(V)=\frac{N-1}{2},
\qquad
C_2^{G_2}(\mathbf7)=2.
\]
For example,
\[
\frac{2(4)(5/4)}{10}=1,\qquad
\frac{2N[(N-1)/2]}{N(N-1)/2}=2,\qquad
\frac{2(7)(2)}{14}=2.
\]
The script checked agreement between the two index calculations.

---

The central-charge inputs are
\[
\begin{array}{c|c|c}
\mathfrak g&\dim\mathfrak g&h^\vee\\ \hline
\mathfrak{su}(2)&3&2\\
\mathfrak{sp}(4)&10&3\\
\mathfrak{su}(3)&8&3\\
\mathfrak{so}(N)&N(N-1)/2&N-2\\
\mathfrak g_2&14&4
\end{array}
\]
and
\[
c_k=\frac{k\dim\mathfrak g}{k+h^\vee}.
\]
For equal levels on \(SU(2)^2\), add the two charges:
\[
c_{k,k}=\frac{6k}{k+2}.
\]

A dagger marks a value strictly exceeding \(d-1\).

| Case | \(d-1\) | \(k=1\) | \(2\) | \(3\) | \(4\) | \(5\) | \(6\) | \(7\) | \(8\) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| (a) Spin \(1\) | \(2\) | \(1\) | \(3/2\) | \(9/5\) | \(2\) | \(15/7^\dagger\) | \(9/4^\dagger\) | \(7/3^\dagger\) | \(12/5^\dagger\) |
| (b) Spin \(3/2\) | \(3\) | \(1\) | \(3/2\) | \(9/5\) | \(2\) | \(15/7\) | \(9/4\) | \(7/3\) | \(12/5\) |
| (c) \(Sp(4)\), \(\mathbf4\) | \(3\) | \(5/2\) | \(4^\dagger\) | \(5^\dagger\) | \(40/7^\dagger\) | \(25/4^\dagger\) | \(20/3^\dagger\) | \(7^\dagger\) | \(80/11^\dagger\) |
| (d) \(SU(3)\), \(\mathbf3\) | \(2\) | \(2\) | \(16/5^\dagger\) | \(4^\dagger\) | \(32/7^\dagger\) | \(5^\dagger\) | \(16/3^\dagger\) | \(28/5^\dagger\) | \(64/11^\dagger\) |
| (e) Equal levels | \(3\) | \(2\) | \(3\) | \(18/5^\dagger\) | \(4^\dagger\) | \(30/7^\dagger\) | \(9/2^\dagger\) | \(14/3^\dagger\) | \(24/5^\dagger\) |
| (f) \(SO(5)\), vector | \(4\) | \(5/2\) | \(4\) | \(5^\dagger\) | \(40/7^\dagger\) | \(25/4^\dagger\) | \(20/3^\dagger\) | \(7^\dagger\) | \(80/11^\dagger\) |
| (f) \(SO(6)\), vector | \(5\) | \(3\) | \(5\) | \(45/7^\dagger\) | \(15/2^\dagger\) | \(25/3^\dagger\) | \(9^\dagger\) | \(105/11^\dagger\) | \(10^\dagger\) |
| (f) \(SO(7)\), vector | \(6\) | \(7/2\) | \(6\) | \(63/8^\dagger\) | \(28/3^\dagger\) | \(21/2^\dagger\) | \(126/11^\dagger\) | \(49/4^\dagger\) | \(168/13^\dagger\) |
| (g) \(G_2\), \(\mathbf7\) | \(6\) | \(14/5\) | \(14/3\) | \(6\) | \(7^\dagger\) | \(70/9^\dagger\) | \(42/5^\dagger\) | \(98/11^\dagger\) | \(28/3^\dagger\) |

For instance, the \(SO(7)\) row comes from
\[
\left(
\frac{21}{6},
\frac{42}{7},
\frac{63}{8},
\frac{84}{9},
\frac{105}{10},
\frac{126}{11},
\frac{147}{12},
\frac{168}{13}
\right).
\]
The table contains the reduced fractions.

The general orthogonal answer for \(k=1,\ldots,8\), with \(D=N(N-1)/2\), is
\[
\left(
\frac N2,\,
N-1,\,
\frac{3D}{N+1},\,
\frac{4D}{N+2},\,
\frac{5D}{N+3},\,
\frac{6D}{N+4},\,
\frac{7D}{N+5},\,
\frac{8D}{N+6}
\right).
\]

The threshold calculations are
\[
\frac{3k}{k+2}>2
\iff 3k>2k+4
\iff k>4,
\]
\[
\frac{10k}{k+3}>3
\iff 10k>3k+9
\iff k>\frac97,
\]
\[
\frac{8k}{k+3}>2
\iff 8k>2k+6
\iff k>1,
\]
\[
\frac{6k}{k+2}>3
\iff 6k>3k+6
\iff k>2,
\]
and
\[
\frac{14k}{k+4}>6
\iff 14k>6k+24
\iff k>3.
\]
For \(SO(N)\),
\[
\frac{kN(N-1)}{2(k+N-2)}>N-1
\]
is equivalent, for \(N\ge5\), to
\[
kN>2k+2N-4
\iff (N-2)k>2(N-2)
\iff k>2.
\]
For spin \(3/2\),
\[
c(SU(2)_k)=3-\frac6{k+2}<3
\]
at every finite level.

As an independent arithmetic check, every table entry was recomputed using
\[
\frac{kD}{k+h}=D-\frac{Dh}{k+h}.
\]

---

(a) For spin \(1\), the central element of \(SU(2)\) acts as
\[
(-1)^{2j}=(-1)^2=1.
\]
The representation is therefore linear under \(SO(3)\), and there is no LSM anomaly.

The \(SU(2)_k\) simple current sends
\[
j\longmapsto \frac k2-j.
\]
In the twisted sector, the diagonal \(SU(2)\) center acts by
\[
(-1)^{2j}\,(-1)^{2(k/2-j)}
=(-1)^k.
\]
Thus a nontrivial center action of translation has trivial mixed anomaly only for even \(k\).

But \(T=1\) has trivial mixed anomaly at every \(k\). Therefore no level is excluded for the spin-\(1\) chain unless one separately fixes a nontrivial center action of translation.

The smallest pure WZW violation is
\[
k=5,\qquad c=\frac{15}{7}=2+\frac17.
\]
If \(T\) is required to be the nontrivial center element, the first one instead is
\[
k=6,\qquad c=\frac94.
\]

For an independent subgroup check, restrict \(SO(3)\) to the commuting \(\pi\)-rotations about two perpendicular axes. Their \(SU(2)\) lifts anticommute, with commutator \(-1\). On spin \(j\), the resulting projective commutator is \((-1)^{2j}\). It is \(+1\) on the site and \((-1)^k\) in the center-twisted WZW sector, reproducing the result.

---

(b) For spin \(3/2\),
\[
(-1)^{2j}=(-1)^3=-1.
\]
The projective class under \(SO(3)\) has order two.

The same twisted-sector calculation gives matching precisely when
\[
(-1)^k=-1,
\qquad\text{or}\qquad
k\ \text{odd}.
\]
Among \(1,\ldots,8\), the admissible levels are
\[
1,3,5,7.
\]

The perpendicular-\(\pi\)-rotation subgroup gives the same sign test independently.

There is no pure WZW central-charge violation because \(c(SU(2)_k)<3\). A small explicit product candidate is
\[
SU(2)_7\otimes M(4,5),
\]
with translation acting through the \(SU(2)\) center and trivially on tricritical Ising. Its charge is
\[
c=\frac73+\frac7{10}
=\frac{70+21}{30}
=\frac{91}{30}
=3+\frac1{30}.
\]
The extra sector is symmetry-neutral, so it leaves the already correct anomaly unchanged.

This is the smallest charge among products of an admissible diagonal \(SU(2)_k\) and one unitary Virasoro minimal model. The proof is given below.

---

(c) For \(Sp(4)\), use Dynkin labels \((a,b)\), with
\[
\mathbf4=(1,0),\qquad \mathbf5=(0,1).
\]
The center character is
\[
q(a,b)=a\pmod2.
\]
In particular, the on-site \(\mathbf4\) has the nontrivial order-two class under \(PSp(4)\cong SO(5)\).

At level \(k\), the integrable weights satisfy
\[
a+b\le k.
\]
The center simple current acts as
\[
J(a,b)=(a,k-a-b).
\]
Therefore
\[
q(J(a,b))=a=q(a,b).
\]
Every representation in the center-twisted sector has trivial projective class:
\[
q(J(a,b))-q(a,b)=0.
\]
Equivalently,
\[
J(0,0)=(0,k)=k\omega_2
\]
is center-neutral for every \(k\).

It follows that
\[
\boxed{\text{No diagonal }Sp(4)_k\text{ matches the spinor-site anomaly by a center action of }T.}
\]

The principal \(SU(2)\) restriction gives a strong independent check. The \(\mathbf4\) restricts to spin \(3/2\), so the site remains projective under the corresponding \(SO(3)\). The embedding index is \(10\), hence the induced affine level is
\[
K=10k.
\]
The \(Sp(4)\) center is the image of the \(SU(2)\) center. Its twisted-sector sign is consequently
\[
(-1)^K=(-1)^{10k}=1,
\]
whereas the site requires \(-1\).

At level one, the conformal-charge check is
\[
c(SU(2)_{10})=\frac{30}{12}=\frac52
=\frac{10}{4}=c(Sp(4)_1).
\]

A symmetry-neutral extra tensor factor cannot fix this. In an ordinary product, the projective charge in a product-twisted sector is the sum of the charges contributed by the factors. A sector on which \(PSp(4)\) acts trivially contributes zero.

One sufficient literal product completion is an additional \(SU(4)_1\) sector, with \(Sp(4)\subset SU(4)\). Its fundamental simple current restricts to the \(\mathbf4\), so it supplies the missing projective class. For example,
\[
Sp(4)_1\otimes SU(4)_1
\]
has
\[
c=\frac52+3=\frac{11}{2}.
\]
Here the physical \(Sp(4)\) acts diagonally, and its diagonal affine level is \(1+1=2\). A literal \(Sp(4)_2\) factor tensored with this anomaly-carrying sector has \(c=4+3=7\).

There is a much smaller completion if a local extension is allowed.

At level one, denote the \(Sp(4)\) primaries by
\[
1,\quad v=\mathbf5,\quad s=\mathbf4.
\]
Their weights and fusion rules are
\[
h_v=\frac12,\qquad h_s=\frac5{16},
\]
\[
v^2=1,\qquad vs=s,\qquad s^2=1+v.
\]
Thus
\[
d_v=1,\qquad d_s=\sqrt2.
\]
The center charges are
\[
q(v)=0,\qquad q(s)=1.
\]

For tricritical Ising \(M(4,5)\),
\[
c=1-\frac6{4\cdot5}=\frac7{10}.
\]
The Kac formula is
\[
h_{r,s}=\frac{(5r-4s)^2-1}{80}.
\]
Take
\[
f=(3,1),\qquad \sigma'=(2,1).
\]
Their weights are
\[
h_f=\frac{11^2-1}{80}=\frac{120}{80}=\frac32,
\]
\[
h_{\sigma'}=\frac{6^2-1}{80}=\frac{35}{80}=\frac7{16}.
\]
Their fusion rules include
\[
f^2=1,\qquad f\sigma'=\sigma',
\qquad (\sigma')^2=1+f.
\]
Hence \(d_f=1\) and \(d_{\sigma'}=\sqrt2\). Independently,
\[
d_{\sigma'}=\frac{\sin(2\pi/4)}{\sin(\pi/4)}=\sqrt2.
\]

Now extend
\[
Sp(4)_1\otimes M(4,5)
\]
by the order-two simple current
\[
B=(v,f).
\]
Its conformal weight is integral:
\[
h_B=\frac12+\frac32=2.
\]
It therefore gives a local bosonic extension.

The object
\[
X=(s,\sigma')
\]
is fixed by \(B\):
\[
B X=(vs,f\sigma')=(s,\sigma')=X.
\]
It is local with respect to the extension because
\[
Q_B(X)=h_B+h_X-h_{BX}=2=0\pmod1.
\]
Its weight and quantum dimension are
\[
h_X=\frac5{16}+\frac7{16}=\frac34,
\qquad
d_X=\sqrt2\,\sqrt2=2.
\]
Fixed-point resolution splits it into two extended primaries \(X_+\) and \(X_-\), each of quantum dimension
\[
d_{X_\pm}=\frac{d_X}{2}=1.
\]
They are therefore simple currents.

Both carry the odd \(Sp(4)\) center character inherited from \(s\). Translation acting through either associated simple-current symmetry gives the required projective class in its twisted sector, by equation (1).

The total central charge is
\[
\boxed{\frac52+\frac7{10}=\frac{16}{5}=3+\frac15.}
\]

This is also a minimum above \(3\), under the stated current-algebra assumption. If
\[
3<c<\frac{16}{5},
\]
the \(Sp(4)\) level must be one, since
\[
c(Sp(4)_2)=4.
\]
The commuting coset Virasoro algebra would then have
\[
\frac12<c_{\rm coset}=c-\frac52<\frac7{10}.
\]
Unitary Virasoro representation theory permits no central charge in this interval. The consecutive nontrivial values below one begin at \(1/2\) and \(7/10\).

For comparison, replacing tricritical Ising by Ising gives the familiar extension
\[
Sp(4)_1\otimes\text{Ising}\ \longrightarrow\ SU(4)_1,
\qquad c=\frac52+\frac12=3.
\]
The new extending fields have weight \(1/2+1/2=1\); their five currents augment the ten \(Sp(4)\) currents to the fifteen \(SU(4)\) currents. This independently checks the central charge and explains how the parton theory supplies the missing symmetry.

---

(d) For \(SU(3)\), a representation with labels \((a,b)\) has triality
\[
q(a,b)=a+2b\pmod3.
\]
The on-site fundamental has \(q=1\), so the LSM class has order three.

A center simple current acts by
\[
J(a,b)=(k-a-b,a).
\]
Its triality is
\[
\begin{aligned}
q(J(a,b))
&=(k-a-b)+2a\\
&=k+a-b\\
&=k+(a+2b)-3b\\
&=q(a,b)+k\pmod3.
\end{aligned}
\]
For translation acting through the \(m\)-th center power, matching requires
\[
mk=1\pmod3. \tag{4}
\]
This has a solution precisely when \(3\nmid k\).

Thus the allowed levels through eight are
\[
1,2,4,5,7,8.
\]
For \(k=2\), choose \(m=2\), since
\[
2\cdot2=4=1\pmod3.
\]
If a particular generator of the center were fixed in advance, only one residue class would match. Allowing the center element to be chosen is why both nonzero residue classes survive.

The principal \(SU(2)\) restriction sends \(\mathbf3\) to spin \(1\), has index \(4\), and induces level \(4k\). Its center maps to the identity in \(SU(3)\). The restricted site anomaly and the induced \(SU(2)\) parity are both trivial. At \(k=1\),
\[
c(SU(2)_4)=\frac{12}{6}=2=c(SU(3)_1).
\]
This check does not detect the order-three anomaly.

A subgroup check that does detect it uses the clock and shift matrices \(C,S\in SU(3)\), satisfying
\[
CS=\omega SC,\qquad \omega^3=1.
\]
Their images commute in \(PSU(3)\). On a triality-\(q\) representation, their projective commutator is \(\omega^q\). It is \(\omega\) on the site and \(\omega^{mk}\) in the twisted sector, reproducing (4).

The smallest pure WZW target is
\[
SU(3)_2,\qquad c=\frac{16}{5}.
\]
A smaller product target is
\[
SU(3)_1\otimes\text{Ising},\qquad c=2+\frac12=\frac52,
\]
with the anomaly entirely supplied by the \(SU(3)_1\) factor.

---

(e) The center character of the on-site representation is
\[
q\left(\frac12,\frac12\right)=(1,1)\in\mathbb Z_2^2.
\]
It has order two.

Although \((-1,-1)\) acts trivially on the site, the two individual central elements act as \(-1\). The faithful on-site group is \(SO(4)\), while the adjoint group in the question is
\[
SO(3)\times SO(3).
\]
Relative to the latter, the projective class is nontrivial.

A product affine algebra has two independent levels. Write them as \(k_1,k_2\). If translation uses center powers \(m_1,m_2\in\{0,1\}\), its twisted-sector class is
\[
(m_1k_1,m_2k_2)\pmod2.
\]
Matching \((1,1)\) requires
\[
k_1,k_2\ \text{odd},
\qquad m_1=m_2=1.
\]
For equal levels, this means \(k=1,3,5,7\) through eight.

Restriction to either factor gives the spin-\(\frac12\) anomaly test independently. Restriction to the diagonal \(SU(2)\) gives
\[
\left(\frac12,\frac12\right)\big|_{\rm diag}=0\oplus1,
\]
which is center-neutral. The diagonal affine level is \(k_1+k_2\), even when both levels are odd, as required.

For equal levels,
\[
c_{1,1}=2,\qquad c_{3,3}=\frac{18}{5}>3.
\]
Without imposing equality,
\[
c_{1,3}=1+\frac95=\frac{14}{5}<3,
\]
whereas
\[
c_{1,5}=1+\frac{15}{7}=\frac{22}{7}=3+\frac17.
\]
If both levels are at least three, the charge is at least \(18/5\). Therefore \((1,5)\), or its interchange, is the smallest pure WZW target.

The parton embedding has levels \((2,2)\), and
\[
c(SU(2)_2)+c(SU(2)_2)=\frac32+\frac32=3=c(SU(4)_1).
\]
Its diagonal \(SU(2)_2^2\) tensor product fails the center test. The extended \(SU(4)_1\) theory succeeds because its fundamental simple current restricts to \((\frac12,\frac12)\).

---

(f) For odd \(N\), \(SO(N)\) has trivial center, and the vector is a linear representation of the adjoint group. There is no LSM anomaly.

For even \(N\),
\[
Z(SO(N))=\{1,-I\},
\]
and \(-I\) acts as \(-1\) on the vector. Thus the vector is projective under
\[
PSO(N)=SO(N)/\{1,-I\},
\]
with an order-two class.

It is useful to derive the complete center-twist test on the simply connected cover.

For \(B_n=\mathfrak{so}(2n+1)\), choose roots in the usual \(e_i\) coordinates. The root lattice is
\[
Q=\mathbb Z^n.
\]
The vector weight \(\omega_1=e_1\) lies in \(Q\), so it has trivial center character. The nontrivial center simple current has highest weight \(k\omega_1\), also neutral. This agrees with the absence of a vector-site anomaly for odd \(N\).

For \(D_n=\mathfrak{so}(2n)\),
\[
Q=\left\{(x_1,\ldots,x_n)\in\mathbb Z^n:\sum_i x_i\text{ even}\right\}.
\]
Introduce the classes
\[
v=[e_1],
\qquad
s=\left[\frac12(e_1+\cdots+e_n)\right],
\]
\[
c=\left[\frac12(e_1+\cdots+e_{n-1}-e_n)\right].
\]
They obey
\[
2v=0,\qquad s-c=v.
\]

If \(n\) is even, then
\[
2s=[e_1+\cdots+e_n]=0,
\]
and similarly \(2c=0\). Therefore
\[
P/Q=\mathbb Z_2^2,\qquad v=s+c.
\]

If \(n\) is odd, then
\[
2s=[e_1+\cdots+e_n]=v.
\]
Consequently
\[
P/Q=\mathbb Z_4,\qquad v=2s,\qquad c=-s.
\]

The possible center-twist charges at level \(k\) are
\[
0,\quad kv,\quad ks,\quad kc.
\]
The site requires \(v\).

For \(N=4r\), multiplication by an even \(k\) annihilates \(\mathbb Z_2^2\), while odd \(k\) preserves it. Hence
\[
\boxed{N\equiv0\pmod4:\quad k\text{ odd}.}
\]

For \(N=4r+2\), solve
\[
km=2\pmod4.
\]
There is a solution for odd \(k\), using \(m=2\), and for \(k\equiv2\pmod4\), using \(m=1\) or \(3\). There is none for \(k\equiv0\pmod4\). Hence
\[
\boxed{N\equiv2\pmod4:\quad 4\nmid k.}
\]

In particular,
\[
\begin{array}{c|c}
N&\text{admissible levels through eight}\\ \hline
5&1,2,3,4,5,6,7,8\\
6&1,2,3,5,6,7\\
7&1,2,3,4,5,6,7,8
\end{array}
\]

The global-form qualification is substantive. The simple current of weight \(k\omega_1\) corresponds to \(-1\in Spin(2n)\), the kernel of \(Spin(2n)\to SO(2n)\). It is nevertheless a nontrivial axial symmetry of the diagonal affine theory. The spinor simple currents correspond to lifts of \(-I\in SO(2n)\).

If one only permits lifts of the nontrivial \(SO(2n)\) center element, the available charges are \(ks,kc\). These never equal \(v\) for \(N=4r\), and equal \(v\) only for \(k\equiv2\pmod4\) when \(N=4r+2\). But this restriction discards an actual internal symmetry of the simply connected diagonal CFT. Passing instead to a quotient WZW model changes the modular invariant, so its twisted sectors must be analyzed afresh.

Several independent checks are available.

For \(N=6\),
\[
Spin(6)\cong SU(4),\qquad V=\Lambda^2\mathbf4.
\]
Its \(SU(4)\) four-ality is \(2\). The matching equation is therefore
\[
km=2\pmod4,
\]
as above.

For \(N=4\), the cover is \(SU(2)^2\), and its vector is \((\mathbf2,\mathbf2)\). The equal-level oddness condition agrees with case (e).

For general \(N\), a three-coordinate \(SO(3)\) subgroup sees
\[
V=\mathbf3+(N-3)\mathbf1,
\]
so the restricted site anomaly is trivial. Its \(SU(2)\) embedding index is \(2\), obtained by dividing the spin-\(1\) index \(4\) into \(SU(N)\) by the orthogonal index \(2\). The induced level \(2k\) is even. This checks the restriction, although this subgroup does not detect the full even-\(N\) projective class.

Finally,
\[
SO(N)_2\subset SU(N)_1
\]
is conformal at the level of central charges:
\[
c(SO(N)_2)
=\frac{2[N(N-1)/2]}{N}
=N-1
=c(SU(N)_1).
\]
The fundamental \(SU(N)_1\) simple current restricts to the vector. Thus the extended parton theory matches the vector-site anomaly even in the \(N=4r\) cases where the diagonal \(SO(N)_2\) center test fails.

The smallest pure WZW violation is always level three:
\[
c_3=\frac{3N(N-1)}{2(N+1)}.
\]
Its excess is
\[
\begin{aligned}
c_3-(N-1)
&=(N-1)\left[\frac{3N}{2(N+1)}-1\right]\\
&=\frac{(N-1)(N-2)}{2(N+1)}>0.
\end{aligned}
\]
For \(N=5,6,7\), these charges are respectively
\[
5,\qquad \frac{45}{7},\qquad \frac{63}{8}.
\]

A smaller explicit completed target, valid for every \(N\ge5\), is
\[
SU(N)_1\otimes\text{Ising},
\qquad
c=N-1+\frac12=N-\frac12.
\]
It contains \(SO(N)_2\), but generally uses its extended modular invariant and an \(SU(N)\) center action of translation.

---

(g) The center of \(G_2\) is trivial. Its \(\mathbf7\) is therefore linear under the adjoint group, and there is no LSM anomaly.

There is no nontrivial center twist to choose. Taking the internal action of translation to be the identity gives the required trivial anomaly at every level.

A useful subgroup check is the principal \(SU(2)\). The \(\mathbf7\) restricts to spin \(3\). Its index into \(SU(7)\) is
\[
\frac23(3)(4)(7)=56.
\]
Dividing by \(\operatorname{ind}(G_2\subset SU(7))=2\) gives the principal embedding index
\[
28.
\]
The induced affine level is \(28k\), always even, and the site is integer-spin. Both restricted anomalies are trivial.

At level one, the charge equality is
\[
c(SU(2)_{28})
=\frac{84}{30}
=\frac{14}{5}
=c(G_{2,1}),
\]
providing a conformal-embedding cross-check.

The smallest pure WZW violation is
\[
G_{2,4},\qquad c=\frac{56}{8}=7.
\]
A product with the smallest nontrivial unitary neutral sector gives
\[
G_{2,3}\otimes\text{Ising},
\qquad c=6+\frac12=\frac{13}{2}.
\]

This last number is not a global minimum over arbitrary neutral sectors. For example,
\[
G_{2,2}\otimes M(3,4)\otimes M(6,7)
\]
has
\[
c=\frac{14}{3}+\frac12+\frac67
=\frac{196+21+36}{42}
=\frac{253}{42}
=6+\frac1{42}.
\]
An even closer explicit product is
\[
G_{2,1}\otimes M(3,4)\otimes M(7,8)^{\otimes2}\otimes M(8,9),
\]
with
\[
\begin{aligned}
c
&=\frac{14}{5}+\frac12+2\frac{25}{28}+\frac{11}{12}\\
&=\frac{1176+210+750+385}{420}\\
&=\frac{2521}{420}
=6+\frac1{420}.
\end{aligned}
\]
These are anomaly-compatible CFTs. Their realization by the specified nearest-neighbor chain is not established.

---

An outer automorphism does not provide a general escape from these conclusions.

Microscopically,
\[
T\,U(g)\,T^{-1}=U(g)
\]
for every \(g\in G\). The IR action must preserve this relation. If an internal action of \(T\) implements an outer automorphism \(\phi\) of the physical continuous symmetry, it instead gives
\[
T\,U(g)\,T^{-1}=U(\phi(g)).
\]
Multiplying it by a center element does not change this conjugation action.

Thus a genuine outer automorphism of the physical \(G\) is incompatible with the translation specified in the question. Examples include exchanging the two independently acting \(SU(2)\) factors, charge conjugation of \(SU(3)\), and the nontrivial diagram automorphisms of the orthogonal algebras.

If the microscopic symmetry were changed to a semidirect product, outer-automorphism twisted sectors would be relevant, but that would be a different anomaly problem.

Extra sectors allow more possibilities. An automorphism of the full CFT may exchange sectors or act on an enlarged current algebra while fixing the physical \(G\) pointwise. Such an operation is allowed. Its anomaly must be calculated in its own twisted sector; it is not determined by the center test for the isolated diagonal \(G_k\) factor. The \(Sp(4)\) extension above is an explicit realization of this broader mechanism.

For an ordinary tensor product with a factorized symmetry action,
\[
r_{\rm total}=r_{\rm WZW}+r_{\rm extra}.
\]
A symmetry-neutral factor has \(r_{\rm extra}=0\). A sector carrying the continuous symmetry can contribute a nonzero class, but if the physical group acts diagonally, its affine levels add as well.

---

There is a useful precise version of the “minimal extra sector” calculation.

For an ordinary nontrivial unitary sector with \(c<1\), Virasoro unitarity gives
\[
c_m=1-\frac6{m(m+1)},\qquad m=3,4,\ldots.
\]
The first values are
\[
\frac12,\quad \frac7{10},\quad \frac45,\quad \frac67,\quad
\frac{25}{28},\quad \frac{11}{12}.
\]

Restricting to a pure diagonal WZW theory or its product with one symmetry-neutral Virasoro minimal model gives the following minima above the proposed bound.

| Case | Smallest theory in this specified class | Total \(c\) |
|---|---|---:|
| (a) Spin \(1\) | \(SU(2)_5\), \(T=1\) | \(15/7\) |
| (b) Spin \(3/2\) | \(SU(2)_7\otimes M(4,5)\) | \(91/30\) |
| (c) \(Sp(4)\), \(\mathbf4\) | None | Neutral products cannot repair the anomaly |
| (d) \(SU(3)\), \(\mathbf3\) | \(SU(3)_1\otimes M(3,4)\) | \(5/2\) |
| (e), independent levels | \(SU(2)_1\otimes SU(2)_5\) | \(22/7\) |
| (e), equal levels | \(SU(2)_3\otimes SU(2)_3\) | \(18/5\) |
| (f), odd \(N\) or \(N\equiv2\pmod4\) | \(SO(N)_2\otimes M(3,4)\) | \(N-\frac12\) |
| (f), \(N\equiv0\pmod4\), \(N\ge8\) | \(SO(N)_3\) | \(\frac{3N(N-1)}{2(N+1)}\) |
| (g) \(G_2\) | \(G_{2,3}\otimes M(3,4)\) | \(13/2\) |

For spin \(3/2\), the only levels needing separate examination are \(5,7,9,11\). Levels \(1,3\) cannot exceed \(3\) with an extra charge below one. The first qualifying minimal models give
\[
\begin{array}{c|c|c}
k&m&c_{\rm total}\\ \hline
5&7&\frac{15}{7}+\frac{25}{28}=3+\frac1{28}\\
7&4&\frac73+\frac7{10}=3+\frac1{30}\\
9&4&\frac{27}{11}+\frac7{10}=3+\frac{17}{110}\\
11&3&\frac{33}{13}+\frac12=3+\frac1{26}
\end{array}
\]
For larger odd \(k\), the Ising addition already qualifies and the total increases monotonically. The minimum is therefore \(3+1/30\).

For spin \(1\), levels \(2,3,4\) give smallest qualifying totals \(11/5,23/10,5/2\), respectively; each exceeds \(15/7\). Level one plus a sector below one does not exceed two.

For the orthogonal theories, level one plus a sector below one cannot exceed \(N-1\) for \(N\ge5\). Whenever level two matches the anomaly, adding Ising gives the minimum in this class. When level two fails, a neutral factor cannot repair it.

The extended \(Sp(4)\) theory at \(16/5\), and the extended \(SU(N)_1\otimes\)Ising theories at \(N-\frac12\), belong to the broader current-algebra setting rather than the literal diagonal-product class in this table.

Without specifying the allowed extra sectors and modular invariants, there is no justified universal list of global minima for the other cases. The two \(G_2\) examples demonstrate why identifying “add Ising” with “smallest total central charge” would be incorrect.

---

For S1, the conclusions about the proposed pure diagonal WZW targets are as follows.

| Proposed target | Excluded by the mixed anomaly? | Reason |
|---|---|---|
| \(Sp(4)_2\) on \(\mathbb C^4\) | Yes, as a standalone diagonal WZW theory | Every center twist is spinor-neutral |
| \(SU(2)_k\), \(k\ge5\), on spin \(1\) | No | Choose \(T=1\); odd \(k\) only fails if a nontrivial center action is separately imposed |
| \(SU(3)_2\) on \(\mathbb C^3\) | No | The squared center generator gives \(2\cdot2=1\pmod3\) |
| \(SO(N)_3\) on the vector | No, with the covering-group axial symmetry described above | Odd \(N\) has no site anomaly; for even \(N\), the vector simple-current charge is \(3v=v\) |
| \(G_{2,k}\), \(k\ge4\), on the \(\mathbf7\) | No | The center and site projective class are trivial |

The smallest pure WZW violations are the ones in the first table. Useful smaller completed targets are
\[
Sp(4):\quad c=\frac{16}{5}
\]
through the tricritical-Ising extension,
\[
SU(3):\quad SU(3)_1\otimes\text{Ising},\quad c=\frac52,
\]
and
\[
SO(N):\quad SU(N)_1\otimes\text{Ising},\quad c=N-\frac12.
\]
For \(G_2\), the examples above reach \(c=6+1/42\) and \(c=6+1/420\), without a claim of global optimality.

The \(Sp(4)\) construction at \(16/5\) has affine level one. It does not rehabilitate the standalone \(Sp(4)_2\) target.

---

For S2, anomaly matching does not imply
\[
k\le\operatorname{ind}(R).
\]
The matching conditions are congruences, and their nonempty solution sets here contain arbitrarily large levels. In other cases, such as the \(Sp(4)\) spinor, the relevant center action fails at every level.

Nor is the diagonal theory \(G_{\operatorname{ind}(R)}\) always anomaly-compatible. Explicit counterexamples are

\[
\begin{array}{c|c|c}
R&\text{Parton-induced affine level}&\text{Failure of diagonal center matching}\\ \hline
SU(2)\text{ spin }3/2&10&\text{Even level, but the site class is odd}\\
Sp(4)\ \mathbf4&1&\text{Center twists are neutral at every level}\\
SU(2)^2\ (\frac12,\frac12)&(2,2)&\text{Both required projective charges vanish}\\
SO(4r)\text{ vector}&2&\text{Multiplication by }2\text{ annihilates the relevant center classes}
\end{array}
\]

What is anomaly-compatible is the full parton theory \(SU(d)_1\), with its appropriate translation action.

Indeed, its fundamental simple current transforms as \(\mathbf d\). Restricting \(\mathbf d\) to \(G\) gives \(R\), so its twisted-sector projective character restricts to
\[
q(R),
\]
the microscopic class. This argument concerns the full \(SU(d)_1\) spectrum and symmetry action. It does not imply that the embedded \(G_{\operatorname{ind}(R)}\) factor, with a diagonal invariant and translation restricted to its own center, matches by itself.

All anomaly exclusions above follow from the twisted-sector charge calculation. The central charges, dimensions, and indices were checked independently as described. None of the admissible CFTs has been shown here to arise from the specified nearest-neighbor Hamiltonian. Their use as counterexample targets is a search proposal, not a proof of a violation of the central-charge conjecture.

[gpt_solve: backend=codex model=gpt-6-astra effort=high transcript=20261027T160801-70730]
