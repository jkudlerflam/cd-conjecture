| Question | Result |
|---|---|
| Extended chiral algebra | Eight primaries. Its modular category is \(\overline{\mathrm{Fib}}\boxtimes\overline{\mathrm{Sem}}\boxtimes\overline{\mathrm{Sem}}\). This is a statement about modular data, not a factorization of the CFT Hamiltonian. |
| Fixed-point resolution | Each of \((s,\sigma)\) and \((s,\sigma')\) splits into two extended primaries. |
| TCI chiral-parity anomalies | The parity reversing the left supercurrent has \(\nu=7=-1\pmod8\). The parity reversing the right supercurrent has \(\nu=1\). |
| Translation inherited from the \(c=3\) phase | Five-Majorana left parity times TCI right parity. Its internal action has order four. It exchanges the resolved Ramond modules. |
| Relevant symmetric perturbations | Two, both with \(\Delta=6/5\): \(\varepsilon'\) and \(\varepsilon\sum_a i\psi_L^a\psi_R^a\). This count is the same with Sp(4) or only SU(2). |
| Symmetric scalar marginal perturbations | One with Sp(4), two with SU(2). Velocity changes are additional, nonscalar marginal deformations. |
| Local codimension | Two, within a basin where the marginal couplings flow to zero. A generic one-parameter Sp(4) family cannot reach the point. A two-parameter SU(2) family has enough parameters, but existence is not guaranteed. |
| Leading symmetric RG eigenvalues | \(y_1=y_2=4/5\). A trajectory controlled by these directions has \(\xi\sim|\delta|^{-5/4}\), up to marginal corrections. |
| Useful spectroscopic signals | Ramond scalar towers at \(\Delta=7/10,3/2\), at momenta \(\pm\pi/2\). Pure-TCI singlets at \(\Delta=1/5,6/5\), at \(\pi,0\), with energy ratio \(6\) even when the velocities differ. |
| Limit of what can be concluded | Symmetry, spectrum, and perturbative OPEs do not determine the full strong-coupling phase diagram or establish a first-order boundary. Those conclusions require additional dynamics and a lattice-to-continuum map. |

I use the diagonal bosonic theory specified by the diagonal GSO projection. “Order four translation” below means the internal part of translation. The actual lattice translation also moves the continuum coordinate.

All calculations below are closed-book. I ran in-memory Python calculations using rational arithmetic for weights and spins, Kac fusion rules, modular transformations, and free-fermion character counting. No files or external sources were consulted.

1. Construction of the extension

Write
\[
e=\varepsilon,\qquad f=\varepsilon',\qquad g=\varepsilon'',
\qquad \varphi=\frac{1+\sqrt5}{2}.
\]

First check the input weights independently.

For \(M(4,5)\),
\[
h_{r,s}=\frac{(5r-4s)^2-1}{80},
\qquad
(r,s)\sim(4-r,5-s).
\]
The assignments
\[
1=(1,1),\quad e=(1,2),\quad f=(1,3),\quad
g=(1,4),\quad \sigma=(2,2),\quad \sigma'=(2,1)
\]
give
\[
0,\quad \frac1{10},\quad \frac35,\quad
\frac32,\quad \frac3{80},\quad \frac7{16}.
\]

For \(\mathfrak{so}(5)_1\), take long roots to have squared length two,
\[
\rho=\frac{3e_1+e_2}{2},\qquad
\lambda_v=e_1,\qquad \lambda_s=\frac{e_1+e_2}{2}.
\]
Then
\[
C_2(\lambda)=\frac{(\lambda,\lambda+2\rho)}2
\]
gives \(C_2(v)=2\), \(C_2(s)=5/4\). Since \(h^\vee=3\),
\[
h_v=\frac2{1+3}=\frac12,\qquad
h_s=\frac{5/4}{1+3}=\frac5{16}.
\]
These also follow from a single NS fermion and the Ramond ground-state energy of five Majoranas.

The central charge is
\[
c=\frac{10}{4}+\left(1-\frac6{4\cdot5}\right)
=\frac52+\frac7{10}=\frac{16}{5}.
\]

The Sp(4) fusion rules are
\[
v^2=1,\qquad vs=s,\qquad s^2=1+v.
\]
The TCI fusion rules needed below, including all nontrivial products, are
\[
\begin{array}{c|ccccc}
 &e&f&g&\sigma&\sigma'\\ \hline
e&1+f&e+g&f&\sigma+\sigma'&\sigma\\
f&&1+f&e&\sigma+\sigma'&\sigma\\
g&&&1&\sigma&\sigma'\\
\sigma&&&&1+e+f+g&e+f\\
\sigma'&&&&&1+g
\end{array}
\]
They follow from the two truncated SU(2) fusion ranges in the Kac formula. I separately checked them with the minimal-model modular \(S\) matrix and Verlinde’s formula.

The extension current is
\[
J=(v,g),\qquad h_J=2.
\]
Its monodromy charge is
\[
Q_J(a,x)=h_J+h_a+h_x-h_{va}-h_{gx}\pmod1.
\]
For the two factors,
\[
Q_v(1)=Q_v(v)=0,\qquad Q_v(s)=\frac12,
\]
and
\[
Q_g(1)=Q_g(e)=Q_g(f)=Q_g(g)=0,\qquad
Q_g(\sigma)=Q_g(\sigma')=\frac12.
\]

Consequently, the local sectors are the eight NS combinations and the two Ramond combinations. The NS combinations form four length-two orbits. Both Ramond combinations are fixed points.

Here is the resulting extended spectrum. The notation in the last column will make the fusion rules compact.

| Extended primary | Restriction to Sp(4)\(_1\times\)TCI | Lowest \(h\) | Quantum dimension | Fusion label |
|---|---|---:|---:|---|
| \(A\) | \((1,1)\oplus(v,g)\) | \(0\) | \(1\) | \(1\) |
| \(B\) | \((v,1)\oplus(1,g)\) | \(1/2\) | \(1\) | \(b\) |
| \(C\) | \((1,e)\oplus(v,f)\) | \(1/10\) | \(\varphi\) | \(b\tau\) |
| \(D\) | \((v,e)\oplus(1,f)\) | \(3/5\) | \(\varphi\) | \(\tau\) |
| \(R_1,R_2\) | Each restricts to \((s,\sigma')\) | \(3/4\) | \(1\) each | \(r_1,r_2\) |
| \(X_1,X_2\) | Each restricts to \((s,\sigma)\) | \(7/20\) | \(\varphi\) each | \(r_1\tau,r_2\tau\) |

For example,
\[
h_{s\sigma}=\frac5{16}+\frac3{80}=\frac7{20},
\qquad
h_{s\sigma'}=\frac5{16}+\frac7{16}=\frac34.
\]

The two summands of \(D\) both start at \(h=3/5\). Its chiral ground space is therefore \(5\oplus1\), of dimension six. This degeneracy will matter when counting scalar perturbations.

At a fixed point, the stabilizer is \(\mathbb Z_2\). Its two one-dimensional representations give the two resolutions. Each resolved module restricts to one copy of the original product module. Thus
\[
\chi_{R_1}=\chi_{R_2}=\chi_s\chi_{\sigma'},
\qquad
\chi_{X_1}=\chi_{X_2}=\chi_s\chi_\sigma.
\]
Equal characters do not mean equal modules. Resolving their modular transformation requires the stabilizer data.

The complete fusion algebra is
\[
b^2=r_1^2=r_2^2=1,\qquad r_1r_2=b,
\qquad
\tau^2=1+\tau,
\]
with all these labels commuting. Thus, for example,
\[
D^2=A+D,\qquad
C^2=A+D,\qquad
CD=B+C,
\]
\[
X_1^2=A+D,\qquad
X_1X_2=B+C,\qquad
R_iD=X_i.
\]

There is a useful check on the fixed-point fusion. If \(r_i^2=b\), the pointed group would be \(\mathbb Z_4\). But a quadratic twist would then require
\[
\theta_b=\theta_{r_i}^{\,4}=(-i)^4=1,
\]
whereas \(\theta_b=-1\). Hence \(r_i^2=1\), and the pointed group is \(\mathbb Z_2\times\mathbb Z_2\).

The modular data can be written explicitly. In pointed order \(1,b,r_1,r_2\),
\[
S_{\rm pt}=\frac12
\begin{pmatrix}
1&1&1&1\\
1&1&-1&-1\\
1&-1&-1&1\\
1&-1&1&-1
\end{pmatrix},
\]
while
\[
S_{\rm Fib}=\frac1{\sqrt{1+\varphi^2}}
\begin{pmatrix}
1&\varphi\\
\varphi&-1
\end{pmatrix}.
\]
Then
\[
S=S_{\rm pt}\otimes S_{\rm Fib},
\qquad
T_{\alpha\alpha}
=\exp\!\left[2\pi i\left(h_\alpha-\frac2{15}\right)\right].
\]
The Fibonacci twist here is
\[
\theta_\tau=e^{6\pi i/5},
\]
which is the conjugate Fibonacci choice.

The total squared quantum dimension is
\[
\mathcal D_{\rm ext}^2=4+4\varphi^2
=14.47213595499958.
\]
Independently, before extension,
\[
\mathcal D_{\rm prod}^2
=4\,[4(1+\varphi^2)].
\]
Extending by an algebra of dimension two divides this by four, giving the same answer.

The numerical checks gave
\[
S^2=1,\qquad (ST)^3=1,
\]
with maximum residuals \(0\) and \(5.6\times10^{-16}\). All extended Verlinde coefficients agreed with the fusion algebra above to \(2.3\times10^{-16}\). The Gauss sum agreed with
\[
e^{2\pi i c/8}=e^{4\pi i/5}
\]
to \(1.2\times10^{-16}\).

2. Translation and the TCI anomaly

The chirality in the phrase “chiral fermion parity” needs care in TCI. It refers to the supercurrent whose sign is reversed. It does not identify the anomaly by counting an elementary free fermion.

Let \(g_L\) reverse the left supercurrent and preserve the right one. The fermionic TCI NS partition function is
\[
Z_{\rm NS}
=|\chi_1+\chi_g|^2+|\chi_e+\chi_f|^2.
\]
A candidate \(g_L\) insertion has the form
\[
Z_{\rm NS}^{g_L}
=(\chi_1-\chi_g)(\bar\chi_1+\bar\chi_g)
+\eta(\chi_e-\chi_f)(\bar\chi_e+\bar\chi_f).
\]
The vacuum fixes the first sign. Modular transformation fixes the second:
\[
\eta=-1.
\]
Indeed, the other sign produces negative coefficients in the spatially twisted expression.

For \(\eta=-1\), modular \(S\) gives
\[
Z_{g_L}
=\sqrt2\left[
\chi_{\sigma'}(\bar\chi_1+\bar\chi_g)
+\chi_\sigma(\bar\chi_e+\bar\chi_f)
\right].
\]
The common \(\sqrt2\) is the familiar normalization associated with the odd-anomaly fermionic defect.

The four ground-state spin differences in this expression are
\[
\frac7{16},\quad
\frac7{16}-\frac32=-\frac{17}{16},\quad
\frac3{80}-\frac1{10}=-\frac1{16},\quad
\frac3{80}-\frac35=-\frac9{16}.
\]
All obey
\[
s=-\frac1{16}\pmod{\frac12}.
\]
With the convention in which a free left Majorana has \(\nu=1\),
\[
s=\frac{\nu}{16}\pmod{\frac12},
\]
so
\[
\nu(g_L)=7,\qquad \nu(g_R)=1.
\]
Here \(g_R=(-1)^F g_L\), and its twisted spins are the negatives of the ones above.

The five free left Majoranas contribute \(\nu=5\). Therefore the inherited translation must be
\[
t=(-1)^{F_{5,L}}\,g_R,
\qquad
\nu(t)=5+1=6\pmod8.
\]
Choosing TCI \(g_L\) instead would give \(5+7=4\pmod8\), which does not match the stipulated continuation from the six-Majorana phase.

For later use, denote the TCI superfield components by
\[
\phi=(e,e),\qquad F=(f,f),
\]
\[
\lambda_L=(f,e),\qquad \lambda_R=(e,f).
\]
Their weights are
\[
\phi:\left(\frac1{10},\frac1{10}\right),\quad
F:\left(\frac35,\frac35\right),
\]
\[
\lambda_L:\left(\frac35,\frac1{10}\right),\quad
\lambda_R:\left(\frac1{10},\frac35\right).
\]
The matched translation acts as
\[
t:\quad
\phi\mapsto-\phi,\quad F\mapsto F,\quad
\lambda_L\mapsto-\lambda_L,\quad
\lambda_R\mapsto\lambda_R,
\]
and
\[
\psi_L^a\mapsto-\psi_L^a,\qquad
\psi_R^a\mapsto\psi_R^a.
\]

The extension generators transform differently from the left-left assignment:
\[
\psi_L^aG_L\mapsto-\psi_L^aG_L,\qquad
\psi_R^aG_R\mapsto-\psi_R^aG_R.
\]
Thus \(t\) is an automorphism of the extended algebra that changes the sign of the added generators. It exchanges the two fixed-point resolutions.

For anomaly \(\nu=6\), the lift to the bosonized local theory obeys
\[
t^2=+1\quad\hbox{on NS local fields},\qquad
t^2=-1\quad\hbox{on Ramond local fields}.
\]
Equivalently, \(t^2\) is the quantum \(\mathbb Z_2\) of the bosonization. Hence the internal symmetry has order four.

The Ramond momentum eigenfields are linear combinations of the two resolved-module sectors, with
\[
t=\pm i,\qquad K=\pm\frac{\pi}{2}.
\]
They should not be identified with \(R_1,R_2\) or \(X_1,X_2\) individually.

For NS fields, the complete branchwise translation rule is as follows. The entries refer to left and right product-algebra branches; descendants retain these signs.

| Extended sector | Branches | Same-branch pairing | Different-branch pairing |
|---|---|---:|---:|
| \(A\) | \((1,1),(v,g)\) | \(+\) | \(-\) |
| \(B\) | \((v,1),(1,g)\) | \(-\) | \(+\) |
| \(C\) | \((1,e),(v,f)\) | \(-\) | \(+\) |
| \(D\) | \((v,e),(1,f)\) | \(+\) | \(-\) |

The signs \(+\) and \(-\) mean internal momenta \(0\) and \(\pi\).

For comparison, monodromy with the extended simple current \(r_1\) gives
\[
t_{\rm mono}(\alpha)
=\frac{\theta_{r_1\alpha}}{\theta_{r_1}\theta_\alpha}.
\]
It acts by
\[
(A,B,C,D)\mapsto(+,-,-,+),
\]
and gives opposite signs to the two resolved Ramond sectors. Its order is two. It is constant across each extended module, unlike the matched translation above.

This alternative is the assignment associated with five-Majorana left parity times TCI \(g_L\), up to interchanging the resolutions. It has anomaly \(4\), not \(6\). Two direct distinctions are:

- Its Ramond towers occur at \(0,\pi\), rather than \(\pm\pi/2\).
- It leaves the extension generators invariant, whereas the inherited translation makes them odd.

3. All scalar local operators through dimension two

I give time-reversal signs for Hermitian components in a real tensor basis. At nonzero momenta \(\pm\pi/2\), inversion and time reversal exchange the two momentum sectors, so a simultaneous assignment of \(t\), \(I_s\), and \(I_b\) eigenvalues is impossible.

Define
\[
M_{ab}=i\psi_L^a\psi_R^b,\qquad
M=\sum_a M_{aa}.
\]
Its symmetric and antisymmetric parts transform as
\[
M_{(ab)}:\ 1\oplus14,\qquad
M_{[ab]}:\ 10.
\]
The principal-SU(2) restrictions are
\[
4\downarrow=[3/2],\qquad
5\downarrow=[2],\qquad
10\downarrow=[1]\oplus[3],
\]
\[
14\downarrow=[2]\oplus[4].
\]

The two vector operators at dimension \(6/5\) can be written
\[
K_\pm^a=i\psi_L^a\lambda_R
\ \pm\ i\lambda_L\psi_R^a.
\]

The reflection and time-reversal matching used here is
\[
I_s:\ \psi_L\leftrightarrow\psi_R,\qquad
\lambda_L\leftrightarrow\lambda_R,
\]
\[
\Theta:\ \psi_L\leftrightarrow\psi_R,\qquad
\lambda_L\leftrightarrow-\lambda_R,
\]
together with
\[
I_s\phi=-\phi,\quad \Theta\phi=\phi,\qquad
I_sF=F,\quad \Theta F=F.
\]
The relative minus sign for the singlet fermion under \(\Theta\) is fixed by the SU(4) matching: the currents \(i\psi^a\psi^6\) form the physical spin-two multipoles and are time-reversal even, whereas the Sp(4) currents have spins one and three and are time-reversal odd.

For example,
\[
I_sM_{ab}=-M_{ba},\qquad
\Theta M_{ab}=M_{ba}.
\]
The singlet \(M\) therefore has the same parities as a staggered invariant bond operator. The same matching assigns these parities to \(\phi\). Since \(e\times e\) contains \(f\), \(F\) is even under all three discrete lattice symmetries.

The complete scalar list is:

| Operator family | \(\Delta\) | Sp(4) representation | Principal SU(2) spins | \(t\) | \(I_s\) | \(I_b\) | \(\Theta\) |
|---|---:|---|---|---:|---:|---:|---:|
| Identity | \(0\) | \(1\) | \(0\) | \(+\) | \(+\) | \(+\) | \(+\) |
| \(\phi\) | \(1/5\) | \(1\) | \(0\) | \(-\) | \(-\) | \(+\) | \(+\) |
| \(X\) Ramond pair | \(7/10\) | Two copies of \(1\oplus5\oplus10\) | Two copies of \(0,1,2,3\) | \(\pm i\) | Exchanges momenta | Exchanges momenta | See below |
| \(M_{(ab)}\) | \(1\) | \(1\oplus14\) | \(0,2,4\) | \(-\) | \(-\) | \(+\) | \(+\) |
| \(M_{[ab]}\) | \(1\) | \(10\) | \(1,3\) | \(-\) | \(+\) | \(-\) | \(-\) |
| \(F\) | \(6/5\) | \(1\) | \(0\) | \(+\) | \(+\) | \(+\) | \(+\) |
| \(\phi M_{(ab)}\) | \(6/5\) | \(1\oplus14\) | \(0,2,4\) | \(+\) | \(+\) | \(+\) | \(+\) |
| \(\phi M_{[ab]}\) | \(6/5\) | \(10\) | \(1,3\) | \(+\) | \(-\) | \(-\) | \(-\) |
| \(K_+^a\) | \(6/5\) | \(5\) | \(2\) | \(-\) | \(-\) | \(+\) | \(-\) |
| \(K_-^a\) | \(6/5\) | \(5\) | \(2\) | \(-\) | \(+\) | \(-\) | \(+\) |
| \(R\) Ramond pair | \(3/2\) | Two copies of \(1\oplus5\oplus10\) | Two copies of \(0,1,2,3\) | \(\pm i\) | Exchanges momenta | Exchanges momenta | See below |
| Symmetric part of \(J_L^A J_R^B\) | \(2\) | \(1\oplus5\oplus14\oplus35_S\) | Given below | \(+\) | \(+\) | \(+\) | \(+\) |
| Antisymmetric part of \(J_L^A J_R^B\) | \(2\) | \(10\oplus35_A\) | Given below | \(+\) | \(-\) | \(-\) | \(-\) |

Here
\[
35_S\downarrow=[0]\oplus[2]\oplus[3]\oplus[4]\oplus[6],
\]
\[
35_A\downarrow=[1]\oplus[2]\oplus[3]\oplus[4]\oplus[5].
\]
In Sp(4) Dynkin labels with \(4=(1,0)\), these are
\[
35_S=(4,0),\qquad 35_A=(2,1).
\]

The full current-current spin content is
\[
10\otimes10\downarrow
=2[0]\oplus2[1]\oplus4[2]\oplus3[3]
\oplus3[4]\oplus[5]\oplus[6].
\]
Its dimension is
\[
2+6+20+21+27+11+13=100.
\]

For the Ramond rows, choose momentum fields \(O_{+,J}\), \(O_{-,J}=O_{+,J}^\dagger\), and a phase convention such that
\[
tO_{\pm,J}=\pm i\,O_{\pm,J},
\]
\[
I_bO_{+,J}=O_{-,J},\qquad
I_sO_{+,J}=iO_{-,J}.
\]
The transformations of \(O_-\) follow by conjugation. These equations give \(t=I_bI_s\).

The standard continuation of the lattice multipoles gives
\[
\Theta O_{+,J}=(-1)^J O_{-,J}.
\]
Thus the Hermitian real and imaginary quadratures are time-reversal even for \(J=0,2\), and odd for \(J=1,3\). They have opposite bond-inversion parities; site inversion mixes them. The \(f\times\sigma\) fusion channel containing \(\sigma'\), together with the evenness of \(F\), gives the same assignment for both Ramond towers.

The Ramond multipole matching assumes the usual nonzero lattice overlaps continued from SU(4). Modular data alone would not establish those overlaps. The absence of a translation-invariant Ramond perturbation does not depend on that assumption.

To see that the table is exhaustive, inspect the chiral weights:

- \(A\) contributes the identity and \(h=\bar h=1\) current products.
- \(B\) contributes \((v,1)_L(v,1)_R\), at \(\Delta=1\).
- \(C\) contributes only \((1,e)_L(1,e)_R\) below scalar dimension \(2\).
- \(D\) has the full \((5\oplus1)_L\otimes(5\oplus1)_R\) ground space at \(\Delta=6/5\), giving \(25+1+5+5=36\) fields.
- \(X_i,R_i\) contribute their \(4_L\otimes4_R\) ground spaces.
- The first scalar descendant of \(\phi\) has dimension \(1/5+2=11/5\), already above the cutoff.

The relevant invariant perturbation is therefore
\[
\delta\mathcal L_{\rm rel}
=uF+w\,M\phi.
\]
Both operators have
\[
\Delta=\frac65.
\]
There is no additional SU(2) singlet inside \(14\), \(10\), or \(5\). Reducing Sp(4) to SU(2) consequently does not increase the relevant count.

With Sp(4), the scalar marginal perturbation is
\[
g\sum_{A=1}^{10}J_L^A J_R^A.
\]
With only SU(2), it splits into
\[
g_1\sum_{m=-1}^{1}(-1)^mJ^{(1)}_{L,m}J^{(1)}_{R,-m}
+
g_3\sum_{m=-3}^{3}(-1)^mJ^{(3)}_{L,m}J^{(3)}_{R,-m}.
\]
There are two independent singlets, as the branching calculation also shows.

Thus the relevant codimension is two. This statement presumes entry into a marginally irrelevant basin. A marginally relevant coupling can obstruct reaching the fixed point without appearing in the relevant-operator count.

4. What the local RG analysis determines

For a scalar perturbation of dimension \(\Delta\),
\[
y=2-\Delta.
\]
The full scalar spectrum therefore has the following linear eigenvalues:
\[
\begin{array}{c|cccccc}
\Delta&1/5&7/10&1&6/5&3/2&2\\ \hline
y&9/5&13/10&1&4/5&1/2&0
\end{array}
\]
Only the two \(y=4/5\) directions are invariant under all the lattice symmetries.

Their leading OPE structure follows directly from
\[
f\times f=1+f,\qquad
f\times e=e+g,\qquad
e\times e=1+f,
\]
and \(v\times v=1\).

Let
\[
\mathcal O_u=F,\qquad \mathcal O_w=M\phi,\qquad
\mathcal O_g=\sum_{a<b}J_L^{ab}J_R^{ab}.
\]
Then:

| Product | Generated couplings among those retained |
|---|---|
| \(\mathcal O_u\mathcal O_u\) | Identity and \(u\) |
| \(\mathcal O_u\mathcal O_w\) | \(w\) |
| \(\mathcal O_w\mathcal O_w\) | Identity, \(u\), \(g\), and stress-tensor terms |
| \(\mathcal O_g\mathcal O_w\) | \(w\) |
| \(\mathcal O_g\mathcal O_g\) | \(g\) |

In particular, \(w^2\) does not generate \(w\): the Sp(4) vector fusion does not return \(v\).

The current-current term in \(M^2\) can be seen without fusion notation. With
\[
J_L^{ab}=i\psi_L^a\psi_L^b,
\]
normal ordering gives
\[
:M^2:
=\sum_{a,b}\psi_L^a\psi_L^b\psi_R^a\psi_R^b
=-2\sum_{a<b}J_L^{ab}J_R^{ab}.
\]
Meanwhile the double contraction gives \(5/|z|^2\), in unit-normalized fermion conventions. Combining these terms with \(e\times e=1+f\) gives the \(u\) and \(g\) channels above.

A schematic beta function is consequently
\[
\dot u=\frac45u+a\,u^2+b\,w^2+\cdots,
\]
\[
\dot w=\frac45w+c\,uw+d\,gw+\cdots,
\]
\[
\dot g=e\,g^2+f\,w^2+\cdots.
\]
The coefficients depend on normalization, and quadratic terms involving relevant couplings also depend on the choice of RG coordinates. They should not be used to infer a global phase portrait without further analysis.

The familiar dynamical interpretation of the TCI \(f=\phi_{1,3}\) perturbation supplies a useful reference trajectory. For \(w=0\), one sign of \(u\) gives the massless flow
\[
\mathrm{TCI}\longrightarrow\mathrm{Ising},
\qquad
c:\frac7{10}\longrightarrow\frac12.
\]
The spectator Sp(4) sector then gives \(c=3\). The other sign gives a phase with broken chiral \(\mathbb Z_2\), interpreted here as translation breaking.

A supersymmetric Landau-Ginzburg picture explains this distinction through
\[
W(\phi)=a\phi+b\phi^3,\qquad
V(\phi)\propto(a+3b\phi^2)^2.
\]
One sign has two minima exchanged by \(\phi\mapsto-\phi\); the other has the massless Majorana associated with the supersymmetry-breaking branch. This is a physical explanation of the standard deformation, not a derivation of its full RG flow from fusion rules.

There are two qualifications in the glued theory.

First, on the translation-broken branch, \(w\ne0\) produces
\[
m_5\propto w\langle\phi\rangle
\]
for the five Majoranas. Generically this gaps the Sp(4) sector as well. At \(w=0\), the ordered singlet sector can coexist with gapless Sp(4), provided its marginal coupling is harmless.

Second, on the massless branch,
\[
\phi\longrightarrow \text{Ising energy},
\]
so \(M\phi\) becomes the marginal mixed interaction between the five Majoranas and the singlet Majorana. It cannot be discarded merely because the TCI sector has flowed to Ising.

In the \(c=3\) theory, choose a normalization and sign for the two Sp(4)-symmetric Gross-Neveu couplings such that positive isotropic coupling is marginally relevant. Their one-loop equations are
\[
\dot g_5=3g_5^2+g_\times^2,\qquad
\dot g_\times=4g_5g_\times.
\]
The factors have a direct flavor-counting origin: three other flavors in the five-component block, one singlet, and four flavors contributing to the mixed vertex.

These equations have a massless basin
\[
g_5<0,\qquad |g_\times|\le -g_5
\]
at weak coupling. A check is
\[
\frac{d}{d\ell}(g_5^2-g_\times^2)
=6g_5(g_5^2-g_\times^2),
\]
so the rays \(g_\times=\pm g_5\) are invariant. Outside the basin the perturbative flow runs toward strong coupling.

This gives a controlled reason for a marginal, BKT-like exit from the \(c=3\) phase when the phase boundary is reached at weak coupling. It does not prove that every lattice boundary is of that type. A first-order transition can preempt it at strong coupling.

Accordingly, the strongest justified conclusions are:

- The glued point has two symmetric relevant directions with the same eigenvalue \(4/5\).
- A decoupled reference direction connects to the \(c=3\) phase.
- The translation-broken branch generically masses the five Majoranas when \(w\ne0\).
- Away from the multicritical point, a weak-coupling exit from \(c=3\) is controlled by marginal interactions.
- Which microscopic trajectories instead encounter a first-order boundary is not fixed by the stated CFT data.

5. Periodic-ring spectroscopy, including unequal velocities

Take \(L\in4\mathbb Z\) and the vacuum sector. Other length classes can select different finite-size sectors.

For a product-algebra state,
\[
E-E_0=\frac{2\pi}{L}
\left(v_5\Delta_5+v_T\Delta_T\right)+o(L^{-1}),
\]
\[
P-P_0=K+\frac{2\pi}{L}s,
\]
where
\[
\Delta_5=h_5+\bar h_5+n_5+\bar n_5,
\]
\[
\Delta_T=h_T+\bar h_T+n_T+\bar n_T,
\]
\[
s=h_5+h_T+n_5+n_T
-\bar h_5-\bar h_T-\bar n_5-\bar n_T.
\]
The internal momentum \(K\) is fixed by translation as above.

The ground-state Casimir term is
\[
E_0(L)=e_\infty L-\frac{\pi}{6L}
\left(\frac52v_5+\frac7{10}v_T\right)+o(L^{-1}).
\]
Unequal velocities preserve this product-sector bookkeeping but remove the common relativistic velocity. In particular, an extended-algebra multiplet need not remain energy-degenerate across its product branches.

The scalar ground levels are:

| Family | \((\Delta_5,\Delta_T)\) | Total \(\Delta\) | \(K\) | SU(2) spins |
|---|---:|---:|---:|---|
| \(\phi\) | \((0,1/5)\) | \(1/5\) | \(\pi\) | \(0\) |
| \(X\), each momentum | \((5/8,3/40)\) | \(7/10\) | \(\pm\pi/2\) | \(0,1,2,3\) |
| \(M_{ab}\) | \((1,0)\) | \(1\) | \(\pi\) | \(0,1,2,3,4\) |
| \(F\) | \((0,6/5)\) | \(6/5\) | \(0\) | \(0\) |
| \(M_{ab}\phi\) | \((1,1/5)\) | \(6/5\) | \(0\) | \(0,1,2,3,4\) |
| \(K_+,K_-\) | \((1/2,7/10)\) | \(6/5\) | \(\pi\) | Two spin-two multiplets |
| \(R\), each momentum | \((5/8,7/8)\) | \(3/2\) | \(\pm\pi/2\) | \(0,1,2,3\) |
| \(J_LJ_R\) | \((2,0)\) | \(2\) | \(0\) | \(2[0]+2[1]+4[2]+3[3]+3[4]+[5]+[6]\) |

For instance, the three nominally \(\Delta=6/5\) energies are proportional to
\[
\frac65v_T,\qquad
v_5+\frac15v_T,\qquad
\frac12v_5+\frac7{10}v_T.
\]
They coincide only when \(v_5=v_T\).

The following construction specifies every level through \(\Delta=11/5=2.2\), including descendants, momenta, multiplicities, and unequal-velocity energies.

Use these chiral branches:

| Sector | Branch notation | \((h_5,h_T)\) | Ground spin |
|---|---|---:|---|
| \(A\) | \(a=(1,1)\) | \((0,0)\) | \(0\) |
|  | \(j=(v,g)\) | \((1/2,3/2)\) | \(2\) |
| \(B\) | \(p=(v,1)\) | \((1/2,0)\) | \(2\) |
|  | \(q=(1,g)\) | \((0,3/2)\) | \(0\) |
| \(C\) | \(e_0=(1,e)\) | \((0,1/10)\) | \(0\) |
|  | \(e_1=(v,f)\) | \((1/2,3/5)\) | \(2\) |
| \(D\) | \(d_5=(v,e)\) | \((1/2,1/10)\) | \(2\) |
|  | \(d_1=(1,f)\) | \((0,3/5)\) | \(0\) |
| \(X_i\) | \(x=(s,\sigma)\) | \((5/16,3/80)\) | \(3/2\) |
| \(R_i\) | \(r=(s,\sigma')\) | \((5/16,7/16)\) | \(3/2\) |

Here “grade” means an integer descendant level above the indicated product branch. The allowed left-right pairs are:

| Sector | Pair | Allowed total grade \(n_L+n_R\) | Total dimensions | Internal \(K\) |
|---|---|---:|---|---|
| \(A\) | \(a\bar a\) | \(0,1,2\) | \(0,1,2\) | \(0\) |
| \(A\) | \(j\bar a,\ a\bar j\) | \(0\) | \(2\) | \(\pi\) |
| \(B\) | \(p\bar p\) | \(0,1\) | \(1,2\) | \(\pi\) |
| \(B\) | \(p\bar q,\ q\bar p\) | \(0\) | \(2\) | \(0\) |
| \(C\) | \(e_0\bar e_0\) | \(0,1,2\) | \(1/5,6/5,11/5\) | \(\pi\) |
| \(C\) | \(e_0\bar e_1,\ e_1\bar e_0\) | \(0,1\) | \(6/5,11/5\) | \(0\) |
| \(C\) | \(e_1\bar e_1\) | \(0\) | \(11/5\) | \(\pi\) |
| \(D\) | \(d_5\bar d_5,\ d_1\bar d_1\) | \(0,1\) | \(6/5,11/5\) | \(0\) |
| \(D\) | \(d_5\bar d_1,\ d_1\bar d_5\) | \(0,1\) | \(6/5,11/5\) | \(\pi\) |
| \(X\) | \(x\bar x\) | \(0,1\) | \(7/10,17/10\) | \(\pm\pi/2\) |
| \(R\) | \(r\bar r\) | \(0\) | \(3/2\) | \(\pm\pi/2\) |

For each total grade, include every distribution between left and right, and between the Sp(4) and TCI factors, using the following chiral characters.

| Sp(4) module | Grade \(0\) | Grade \(1\) | Grade \(2\), when needed |
|---|---|---|---|
| \(1\) | \([0]\) | \([1]\oplus[3]\) | \([0]\oplus[1]\oplus2[2]\oplus[3]\oplus[4]\) |
| \(v\) | \([2]\) | \([1]\oplus[2]\oplus[3]\) | Not needed |
| \(s\) | \([3/2]\) | \([1/2]\oplus[3/2]\oplus[5/2]\oplus[7/2]\) | Not needed |

The corresponding dimensions are
\[
1,\ 10,\ 30;\qquad 5,\ 15;\qquad 4,\ 20.
\]
I checked these by explicitly enumerating free-Majorana Fock states. For example, at vacuum grade two,
\[
5\otimes5\ \oplus\ \wedge^4 5
=(1\oplus10\oplus14)\oplus5,
\]
giving dimension thirty.

The TCI descendant multiplicities required here are:

| TCI module | Grade \(0\) | Grade \(1\) | Grade \(2\), when needed |
|---|---:|---:|---:|
| \(1\) | \(1\) | \(0\) | \(1\) |
| \(e\) | \(1\) | \(1\) | \(1\) |
| \(f\) | \(1\) | \(1\) | Not needed |
| \(g\) | \(1\) | Not needed | Not needed |
| \(\sigma\) | \(1\) | \(1\) | Not needed |
| \(\sigma'\) | \(1\) | Not needed | Not needed |

The \(e\) module has one state, rather than two, at grade two because of its level-two null vector.

Finally, for every left-right pair,
\[
[j_L]\otimes[j_R]
=\bigoplus_{J=|j_L-j_R|}^{j_L+j_R}[J].
\]
Together with the energy and momentum formulas, these tables are an exhaustive level specification through \(2.2\), without suppressing multiplets that split when \(v_5\ne v_T\).

Some prominent nonscalar levels illustrate the rules:

- Sp(4) currents have \(\Delta=1\), spins \(1,3\), and momenta \(\pm2\pi/L\).
- The two TCI descendants of \(\phi\) have \(\Delta=6/5\), spin zero, and momenta \(\pi\pm2\pi/L\).
- Acting on \(\phi\) with a Sp(4) current gives spins \(1,3\), the same total dimension \(6/5\), and momenta \(\pi\pm2\pi/L\), but a different velocity dependence.
- The \(C\)-sector cross branches give spin two at \(\Delta=6/5\), momenta \(\pm2\pi/L\), and energy proportional to \(v_5/2+7v_T/10\).
- The extension generators have \(\Delta=2\), spin two, and momenta \(\pi\pm4\pi/L\).
- The \(B\)-sector cross branches have \(\Delta=2\), spin two, and momenta \(\pm2\pi/L\).
- The \(X\) descendants have \(\Delta=17/10\), momenta \(\pm\pi/2\pm2\pi/L\). A Sp(4) descendant has spin content
  \[
  [0]\oplus3[1]\oplus4[2]\oplus3[3]\oplus2[4]\oplus[5],
  \]
  while a TCI descendant has \([0]\oplus[1]\oplus[2]\oplus[3]\). Their energy increments are respectively \(2\pi v_5/L\) and \(2\pi v_T/L\).

There are several velocity-independent TCI tests. Writing energies relative to the ground state,
\[
\frac{E_F}{E_\phi}=6,
\qquad
\frac{E_{T_{\rm TCI}}}{E_\phi}=10,
\]
\[
\frac{E_{\text{pure-TCI }e\text{ level at }11/5}}{E_\phi}=11,
\qquad
\frac{E_{\partial F}}{E_F}=\frac{11}{6}.
\]
The pure-TCI first descendants of \(\phi\) also have the same leading energy as \(F\), but occur at \(\pi\pm2\pi/L\) instead of zero momentum.

In the \(c=3\) phase, replacing TCI by Ising gives instead
\[
\Delta_{\text{singlet mass}}=1,\qquad
\Delta_{\text{Ramond}}=
2\left(\frac5{16}+\frac1{16}\right)=\frac34.
\]
The mass remains at \(\pi\), and the Ramond towers remain at \(\pm\pi/2\). The pure-Ising stress-to-mass energy ratio is \(2\), rather than the TCI ratio \(10\).

At SU(4)\(_1\), the velocities coincide. The leading scalar towers are
\[
\Delta=\frac34,\quad K=\pm\frac{\pi}{2},\quad
J=0,1,2,3,
\]
and
\[
\Delta=1,\quad K=\pi,\quad
6\otimes6\downarrow
=2[0]\oplus[1]\oplus3[2]\oplus[3]\oplus[4].
\]
The chiral currents comprise spins \(1,2,3\) at \(\Delta=1\). At the glued TCI point, only spins \(1,3\) remain in the dimension-one current algebra.

The calculation therefore supports a definite operator-theoretic conclusion: the proposed glued point is a two-relevant-direction multicritical point, with inherited order-four translation and a distinctive two-velocity spectrum. It does not establish that the one-parameter Sp(4) bond family realizes it, or determine a first-order phase boundary. Those are the remaining lattice and strong-coupling questions.

[gpt_solve: backend=codex model=gpt-6-astra effort=high transcript=20261027T181712-15382]
