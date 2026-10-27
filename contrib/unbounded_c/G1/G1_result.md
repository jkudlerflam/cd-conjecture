# Exchanging large logical spins through one small bond

## Result and scope

For every finite \(N>d\geq 2\), there is a construction with

\[
H^{(1)}=0,\qquad H^{(2)}=a_2 I+b_2P_N,\qquad b_2>0.
\]

The sign can also be reversed. Thus the smallest possible order is \(k=2\), and there is no finite upper bound on \(N\) at second order for any fixed \(d\geq2\). In particular, a universal lower bound growing like \(\log_d N\) is false in the stated setting. The construction uses the freedom to engineer the excited spectrum and its matrix elements, not a large port.

The explicit examples below have \((d,N)=(2,3),(2,4),(3,4)\). Their local Hilbert dimensions are 22, 40, and 69. Their gaps are all 1. The first example has the rational coefficient

\[
H^{(2)}=-\frac{7177}{8192}I+\frac{27}{8192}P_3.
\]

These statements are proved below. Full-Hamiltonian numerical diagonalization verifies the leading swap and its eigenspaces. At finite nonzero coupling the multiplets have higher-order internal splittings. Exact degeneracy at finite coupling is not claimed or needed for the specified perturbative property.

No exact nearest-neighbor chain realization is supplied. Section 8 identifies what such a realization would still need to establish. That is the unresolved part of the locality question.

This work is closed book. No external sources or preexisting files were consulted. We use the convention

\[
H(\lambda)=H_0+\lambda h.
\]

The question defines \(V=\lambda h\) and later writes \(H_0+\lambda V\). We interpret the requested numerical test as \(H_0+\lambda h\), with a single perturbation parameter. All reported orders refer to this convention.

## 1. A spectral-weight realization lemma

The following elementary dilation is the local ingredient.

Lemma. Suppose a physical block contains mutually orthogonal isometries

\[
G:\mathbb C^N\longrightarrow\mathbb C^D,\qquad
S_a:\mathbb C^N\longrightarrow\mathbb C^D,\quad a=1,\ldots,r.
\]

Write \(S=[S_1\ \cdots\ S_r]\). Given positive matrices

\[
M_e>0\quad\hbox{on }\mathbb C^r\otimes\mathbb C^N,
\qquad \sum_{e=1}^{L}M_e=I_{rN},
\]

and positive energies \(E_e\), assume \(D-N\geq LrN\). There is a positive local Hamiltonian with kernel exactly \(\operatorname{ran}G\), and with excited spectral projectors \(\mathcal P_e\) satisfying

\[
S^\dagger\mathcal P_e S=M_e.
\]

Its nonzero energies are \(E_e\), together with an arbitrary positive padding energy if necessary.

Proof and explicit recipe. Set \(q_0=D-N\), and form the \(q_0\times rN\) matrix

\[
W=\begin{pmatrix}\sqrt{M_1}\\ \vdots\\ \sqrt{M_L}\\0\end{pmatrix}.
\]

Then \(W^\dagger W=\sum_eM_e=I\). Complete its columns to a unitary \(U\), with \(U[:,1:rN]=W\). Choose a physical orthonormal frame \(F_Q\) for \(G^\perp\) whose first \(rN\) columns are \(S\). Define

\[
D_E=\operatorname{diag}(E_1I_{rN},\ldots,E_LI_{rN},E_{\rm pad}I),
\qquad
H_{\rm loc}=F_QU^\dagger D_EUF_Q^\dagger.
\tag{1}
\]

It annihilates \(G\), and is positive definite on \(G^\perp\). The excited eigenvectors are the columns of \(F_QU^\dagger\). Let \(R_e\) select the corresponding \(rN\) spectral coordinates, including padding coordinates if their energy coincides with \(E_e\). The padding rows of \(W\) vanish, so

\[
S^\dagger\mathcal P_eS=W^\dagger R_eW=M_e.
\]

This proves all claims. In the code, positive square roots and full QR decompositions implement (1). The supplied isometry columns replace the first QR columns, preserving the prescribed frame. Different completions can change higher-order terms, but not the second-order result.

## 2. Construction for every \(d\geq2\) and every finite \(N\geq2\)

### 2.1 Port and ground code

Set \(r=d^2-1\). Let \(F_a\) be a traceless Hermitian basis on the port with

\[
\operatorname{tr}(F_aF_b)=d\delta_{ab}.
\]

Use the same port interaction for all logical dimensions,

\[
h=\sum_{a=1}^{r}F_a\otimes F_a=dP_d-I_{d^2}.
\tag{2}
\]

For \(d=2\), this is \(X\otimes X+Y\otimes Y+Z\otimes Z\). Its norm is \(d+1\), independent of \(N\).

In the auxiliary space reserve orthonormal vectors \(|a,i\rangle_K\), with \(a=1,\ldots,d\) and \(i=1,\ldots,N\). Define the ground isometry by

\[
G|i\rangle=\frac1{\sqrt d}\sum_{a=1}^{d}|a\rangle_{\rm port}|a,i\rangle_K.
\tag{3}
\]

For any port operator \(O\),

\[
G^\dagger(O\otimes I)G=\frac{\operatorname{tr}O}{d}I_N.
\tag{4}
\]

In particular, \(G^\dagger F_aG=0\). Define \(S_a=(F_a\otimes I)G\). Equation (4) gives

\[
G^\dagger S_a=0,\qquad S_a^\dagger S_b=\delta_{ab}I_N.
\tag{5}
\]

Thus a port operator reaches \(r\) orthogonal excited syndrome spaces, each carrying the logical space. Both blocks use (3), but their excited Hamiltonians will differ. Equality \(H_A=H_B\) is not assumed in the question and is not imposed here.

### 2.2 Energy coordinates

Let \(m=N^2-1\) and \(c=r^2-1\). Choose

\[
q=\left\lceil\frac{m}{c}\right\rceil,\qquad L=q+1,
\qquad E_e=3^{e-1},\qquad w_e=\frac1L.
\tag{6}
\]

The Cauchy matrix

\[
C_{ef}=\frac1{E_e+E_f}
\tag{7}
\]

is strictly positive definite. Indeed, for real \(v\ne0\),

\[
v^TCv=\int_0^\infty\left(\sum_ev_e e^{-E_et}\right)^2dt>0.
\tag{8}
\]

For the strict inequality, a vanishing integral would make the continuous integrand vanish identically. Ordering the distinct energies, multiplying by the slowest exponential's inverse, and taking \(t\to\infty\) forces its coefficient to vanish. Repeating forces every coefficient to vanish.

Let \(Z\in\mathbb R^{L\times q}\) be an orthonormal basis of the sum-zero subspace, and set

\[
X=Z(Z^TCZ)^{-1/2}.
\tag{9}
\]

Then

\[
\sum_eX_{ej}=0,\qquad X^TCX=I_q.
\tag{10}
\]

For a fully specified choice, the \(j\)-th column of \(Z\), \(j=1,\ldots,L-1\), has its first \(j\) entries equal to \(1/\sqrt{j(j+1)}\), entry \(j+1\) equal to \(-j/\sqrt{j(j+1)}\), and all remaining entries zero.

### 2.3 Positive spectral weights

Let \(A_p\), \(p=1,\ldots,c\), be a traceless Hermitian orthonormal basis of \(r\times r\) matrices, and \(T_\mu\), \(\mu=1,\ldots,m\), such a basis of \(N\times N\) matrices. Their normalizations are

\[
\operatorname{tr}(A_pA_{p'})=\delta_{pp'},\qquad
\operatorname{tr}(T_\mu T_\nu)=\delta_{\mu\nu}.
\]

Assign to each logical index \(\mu\) a distinct pair \((p_\mu,j_\mu)\in\{1,\ldots,c\}\times\{1,\ldots,q\}\). There are enough pairs by (6). The code assigns them in lexicographic order with \(p\) running fastest. It uses the normalized real and imaginary off-diagonal matrix units, followed by the usual traceless diagonal matrices, for both Hermitian bases.

Define

\[
D_{A,e}=\sum_{\mu=1}^{m}X_{e j_\mu}A_{p_\mu}\otimes T_\mu,
\qquad
D_{B,e}=\sum_{\mu=1}^{m}X_{e j_\mu}A_{p_\mu}^{T}\otimes T_\mu.
\tag{11}
\]

The transpose acts on the syndrome matrix only. Both \(D_{A,e}\) and \(D_{B,e}\) are Hermitian. Put

\[
M_* =\max_e\{\|D_{A,e}\|,\|D_{B,e}\|\},\qquad
\tau=\frac{3}{4LM_*}>0,
\tag{12}
\]

and, for \(s\in\{-1,+1\}\), take

\[
M_{A,e}=w_eI+\tau D_{A,e},\qquad
M_{B,e}=w_eI+s\tau D_{B,e}.
\tag{13}
\]

These are legitimate spectral weights, not assumed ones. Equations (10) and (11) imply that they each sum to \(I\). Equation (12) gives

\[
M_{A,e},M_{B,e}\ \geq\frac1{4L}I>0.
\tag{14}
\]

Choose the actual physical dimension and auxiliary dimension as

\[
D=d\left\lceil\frac{N+LrN}{d}\right\rceil,\qquad K=D/d.
\tag{15}
\]

Since \(L\geq2\), \(D\geq d^2N\), so the code (3) and the syndrome spaces (5) fit. Also \(D-N\geq LrN\), as required by the realization lemma. Apply (1) separately to (13), with \(E_{\rm pad}=E_L\). This completes the explicit definition of \(H_A,H_B,h\).

The ground degeneracy is exactly \(N\), not at least \(N\). All other eigenvalues are at least 1.

## 3. Calculation of the effective Hamiltonian

Equation (4) immediately gives

\[
H^{(1)}=\sum_a(G^\dagger F_aG)\otimes(G^\dagger F_aG)=0.
\tag{16}
\]

Moreover \(h\) maps a ground-ground state into a state with both blocks excited. If \(M_{A,e}^{ab}\) denotes the logical \(N\times N\) block of \(M_{A,e}\) in syndrome row \(a\), column \(b\), the second-order resolvent therefore gives

\[
H^{(2)}=-\sum_{e,f}C_{ef}\sum_{a,b=1}^{r}
M_{A,e}^{ab}\otimes M_{B,f}^{ab}.
\tag{17}
\]

To see the indices directly, the two insertions of the bond contribute \(F_a\) on the bra side and \(F_b\) on the ket side in both blocks. Their local matrix element is

\[
G^\dagger F_a\mathcal P_e F_bG
=S_a^\dagger\mathcal P_eS_b=M_e^{ab}.
\]

There is no complex conjugate on the second block in (17). This is why the syndrome transpose in (11) matters.

Expand (17) using (13). The constant contribution is

\[
-r\,w^TCw\,I.
\]

Every term linear in \(\tau\) vanishes because it contains \(\sum_a(A_p)_{aa}=0\). For the quadratic term,

\[
\sum_{a,b}(A_p)_{ab}(A_{p'}^T)_{ab}
=\operatorname{tr}(A_pA_{p'})=\delta_{pp'},
\]

and

\[
\sum_{e,f}C_{ef}X_{ej}X_{fj}=1,\qquad
\sum_{e,f}C_{ef}X_{ej}X_{fj'}=0\quad(j\ne j').
\]

Since the pairs \((p_\mu,j_\mu)\) are distinct, the quadratic contribution reduces to

\[
-s\tau^2\sum_{\mu=1}^{N^2-1}T_\mu\otimes T_\mu.
\]

Completeness of a Hermitian Hilbert-Schmidt basis gives

\[
\sum_{\mu=1}^{N^2-1}T_\mu\otimes T_\mu=P_N-\frac1NI.
\tag{18}
\]

Consequently,

\[
\boxed{\quad
H^{(2)}=\left(-r\,w^TCw+\frac{s\tau^2}{N}\right)I
-s\tau^2P_N.
\quad}
\tag{19}
\]

Taking \(s=-1\) gives the antiferromagnetic sign. Taking \(s=+1\) gives the ferromagnetic sign. No traceless one-block term remains, and the connected coefficient matrix is \(-s\tau^2I_{N^2-1}\), of rank \(N^2-1\).

This is an algebraic proof. No approximation to the resolvent, variational assumption, or numerical tuning is involved.

## 4. Small examples

### 4.1 Two-dimensional port and a logical qutrit

Here \(r=N=3\), \(L=2\), \(E_1=1,E_2=3\), and

\[
C=\begin{pmatrix}1/2&1/4\\1/4&1/6\end{pmatrix},\qquad
X=\sqrt6\begin{pmatrix}1\\-1\end{pmatrix},\qquad
w^TCw=\frac7{24}.
\]

Take \(A_\mu=T_\mu\). Write \(\Phi_3=|\Omega_3\rangle\langle\Omega_3|\), where \(|\Omega_3\rangle=3^{-1/2}\sum_i|ii\rangle\). Let

\[
F=P_3-I/3,\qquad J=3\Phi_3-I/3,\qquad \eta=9/64.
\]

For the antiferromagnetic example the four weights are

\[
M_{A,1}=I/2+\eta F,\quad M_{A,2}=I/2-\eta F,
\]
\[
M_{B,1}=I/2-\eta J,\quad M_{B,2}=I/2+\eta J.
\tag{20}
\]

The eigenvalues of \(F\) are \(2/3,-4/3\), and those of \(J\) are \(8/3,-1/3\). Thus all four weights are positive, and the smallest eigenvalue among them is \(1/8\). This checks the positivity directly, without relying on a numerical norm estimate.

Here \(\tau=\eta/\sqrt6=9/(64\sqrt6)\). Equation (19) yields

\[
a_2=-\frac78-\frac{\tau^2}{3}=-\frac{7177}{8192},\qquad
b_2=\tau^2=\frac{27}{8192}>0.
\]

Use \(K=11\), hence \(D=22\), in (3) and (1). The local spectrum is 0 with multiplicity 3, 1 with multiplicity 9, and 3 with multiplicity 10. Equations (1), (2), and (20) specify the complete Hamiltonians.

At second order the antisymmetric triplet has energy \(-1801\lambda^2/2048\), while the symmetric sextet has energy \(-3575\lambda^2/4096\). Their difference is \(27\lambda^2/4096\). Both shifts are negative, consistent with the general negative semidefiniteness of \(H^{(2)}\).

### 4.2 Other explicit instances

The basis order and formulas (6)-(15), implemented in `gadget.py`, fix the remaining examples. Values are rounded here, not in the construction.

| \(d,N\) | \(K\) | \(D=dK\) | Excited energies | \(\tau\) | \(a_2\) | \(b_2>0\) |
| :-- | --: | --: | :-- | --: | --: | --: |
| 2,3 | 11 | 22 | 1,3 | 0.05740991585 | -0.8760986328125 | 0.0032958984375 |
| 2,4 | 20 | 40 | 1,3,9 | 0.02196731124 | -0.529750270320 | 0.000482562763192 |
| 3,4 | 23 | 69 | 1,3 | 0.08353830914 | -2.33507799561 | 0.00697864909480 |

With two excited energies this recipe supports every \(2\leq N\leq r=d^2-1\). That is the capacity of this particular two-energy recipe, not a universal upper bound. For \(N=r\), use \(F=P_r-I/r\), \(J=r\Phi_r-I/r\), and \(\eta=3r/[8(r^2-1)]\) in (20). Then \(\tau=\eta/\sqrt6\). In particular, \(d=3,N=r=8\) has \(\eta=1/21\), \(D=138\), and \(b_2=1/2646\). This last instance is an algebraic consequence of the formulas and was not numerically diagonalized.

## 5. An even narrower bond also suffices

For completeness, unbounded \(N\) does not require the several Pauli channels in (2). A single product bond works. This alternative is less numerically convenient, and is not used for the numerical tables.

For any \(d\geq2\), use

\[
G|i\rangle=|0\rangle|i\rangle_K,\quad
S|i\rangle=|1\rangle|i\rangle_K,\quad
F=|0\rangle\langle1|+|1\rangle\langle0|,\quad h=F\otimes F.
\]

For \(d>2\), let \(F\) vanish on the other port levels. Set \(L=N^2+1=m+2\), choose distinct positive energies, and let \(C\) be (7), \(w_e=1/L\). The vectors \(\mathbf 1\) and \(Cw\) are linearly independent: the entries of \(Cw\) strictly decrease when the corresponding \(E_e\) increases.

Choose a real \(L\times m\) matrix \(Z\) spanning the Euclidean orthogonal complement of both vectors and set \(X=Z(Z^TCZ)^{-1/2}\). Then

\[
\mathbf 1^TX=0,\qquad w^TCX=0,\qquad X^TCX=I_m.
\]

Let

\[
M_{A,e}=w_eI+\epsilon\sum_\mu X_{e\mu}T_\mu,\qquad
M_{B,e}=w_eI+s\epsilon\sum_\mu X_{e\mu}T_\mu.
\]

Choosing

\[
0<\epsilon<\frac{\min_e w_e}{\max_e\|\sum_\mu X_{e\mu}T_\mu\|}
\]

makes every weight positive. They sum to identity. The realization lemma with \(r=1\), and dimension \(D=d\lceil(L+1)N/d\rceil\), gives actual physical blocks. The same resolvent calculation now gives

\[
H^{(1)}=0,\qquad
H^{(2)}=\left(-w^TCw+\frac{s\epsilon^2}{N}\right)I-s\epsilon^2P_N.
\tag{21}
\]

Here the linear terms vanish by \(w^TCX=0\). Thus even an operator-Schmidt-rank-one bond permits both signs for arbitrarily large \(N\). The scalar offset in (21) is necessary; it is consistent with the negative semidefiniteness of the second-order operator.

## 6. Numerical verification

### Procedure

Run

```sh
python3 gadget.py
```

Only NumPy, SciPy, and the Python standard library are used. The program constructs the physical local matrices, checks \(H^{(1)}\), and computes \(H^{(2)}\) directly from the physical ground-to-excited transition amplitudes and denominators \(E_e+E_f\). It also checks the opposite sign at second order.

For the antiferromagnetic examples it diagonalizes the full operator \(H_A\otimes I+I\otimes H_B+\lambda h\), using a matrix-free Hermitian eigensolver. This is a numerical diagonalization of the full finite Hilbert space, not a perturbative truncation. Two extra eigenpairs check isolation of the \(N^2\)-dimensional low-energy band. A QR and Rayleigh-Ritz step restores orthonormality inside the tightly clustered eigenspaces before extracting the effective Hamiltonian. This matters because the complex iterative solver can return slightly nonorthogonal vectors within nearly degenerate clusters.

Let \(\Psi\) contain the \(N^2\) orthonormal low eigenvectors and \(\mathcal E\) their energies. Form the ground overlaps

\[
B=(G_A\otimes G_B)^\dagger\Psi,
\qquad W=(BB^\dagger)^{-1/2}B.
\]

The code evaluates \(W\) using the SVD polar factor. The canonical, or des Cloizeaux, effective Hamiltonian is

\[
H_{\rm dc}=W\mathcal EW^\dagger.
\tag{22}
\]

The least-squares scalar and swap coefficients are

\[
b_{\rm fit}=\frac{\operatorname{tr}[(P_N-I/N)H_{\rm dc}]}{N^2-1},\qquad
a_{\rm fit}=\frac{\operatorname{tr}H_{\rm dc}}{N^2}-\frac{b_{\rm fit}}N.
\tag{23}
\]

For positive \(b_2\), the lowest \(N(N-1)/2\) states are compared with the antisymmetric projector \((I-P_N)/2\), using their polar-projected ground components. This checks the eigenspaces in the fixed logical identification, in addition to checking their energies.

### Results

The direct second-order Frobenius errors \(\|H^{(2)}-a_2I-b_2P_N\|_F\) were respectively \(2.1\times10^{-15}\), \(2.9\times10^{-15}\), and \(5.8\times10^{-15}\). First-order norms were below \(5\times10^{-33}\). The ferromagnetic second-order checks were also below \(7\times10^{-15}\).

The exact-diagonalization sequence is \(\lambda=0.02,0.01,0.005,0.0025,0.00125\). Representative fitted coefficients are

| \(\lambda\) | \(b_{\rm fit}/\lambda^2\), \(d=2,N=3\) | \(d=2,N=4\) | \(d=3,N=4\) |
| --: | --: | --: | --: |
| 0.02 | 0.003402730 | 0.000482304 | 0.007191320 |
| 0.01 | 0.003348940 | 0.000482458 | 0.007085366 |
| 0.005 | 0.003322326 | 0.000482517 | 0.007032093 |
| 0.0025 | 0.003309089 | 0.000482541 | 0.007005391 |
| 0.00125 | 0.003302488 | 0.000482552 | 0.006992025 |
| Predicted limit | 0.003295898 | 0.000482563 | 0.006978649 |

At the smallest coupling, the lower and upper cluster sizes are \(3,6\), \(6,10\), and \(6,10\), as required. The cluster diagnostics are

| \(d,N\) | Mean symmetric-minus-antisymmetric splitting | Antisymmetric width | Symmetric width | Antisymmetric projector Frobenius error |
| :-- | --: | --: | --: | --: |
| 2,3 | \(1.0320\times10^{-8}\) | \(1.918\times10^{-11}\) | \(5.082\times10^{-11}\) | \(1.167\times10^{-3}\) |
| 2,4 | \(1.5080\times10^{-9}\) | \(9.542\times10^{-12}\) | \(1.804\times10^{-11}\) | \(5.650\times10^{-3}\) |
| 3,4 | \(2.1850\times10^{-8}\) | \(1.233\times10^{-10}\) | \(1.856\times10^{-10}\) | \(1.228\times10^{-3}\) |

Log-log fits over all five couplings give

| \(d,N\) | Mean splitting power | \(\|H_{\rm dc}-a_{\rm fit}I-b_{\rm fit}P_N\|_F\) power | \(\|H_{\rm dc}-\lambda^2H^{(2)}\|_F\) power | Projector-error power |
| :-- | --: | --: | --: | --: |
| 2,3 | 2.010 | 3.012 | 3.001 | 1.000 |
| 2,4 | 2.000 | 3.008 | 3.000 | 1.009 |
| 3,4 | 2.010 | 3.011 | 2.996 | 1.002 |

The fitted residuals at the smallest coupling are approximately \(4.78\times10^{-11}\), \(2.01\times10^{-11}\), and \(2.03\times10^{-10}\). Thus the effective matrix agrees with a scalar plus swap through order two; the remaining nonswap component starts at order three for these completions. Both internal widths scale as \(\lambda^3\). The two cluster projectors converge linearly in \(\lambda\), as expected when an \(O(\lambda^3)\) perturbation mixes sectors separated by \(O(\lambda^2)\).

The full numerical records, including all low levels, swap expectations, isolation gaps, and eigenpair residuals, are in `gadget_results.json`. A dense full-matrix cross-check is reproducible with

```sh
python3 gadget.py --cases '2,3' --lambdas 0.005 --method dense --output gadget_dense_check.json
```

The dense and iterative low spectra agreed to better than \(10^{-14}\). No full-spectrum diagonalization for the larger examples is claimed; their complete operators were diagonalized for the required low eigenpairs. Script runtime was well below one minute in total for the reported runs and reruns, within the 40-minute limit.

### Checks and limitations

The local gap is 1 by construction. The full low band is isolated for \(|\lambda|<1/[2(d+1)]\), since the unperturbed gap is 1 and \(\|h\|=d+1\). All tested couplings satisfy this sufficient bound.

Multiplying both local Hamiltonians by an energy scale \(\Delta\) divides \(a_2,b_2\) by \(\Delta\), as required by second-order perturbation theory. Taking \(\tau\to0\) removes the connected term and leaves only a scalar shift. These are algebraic consistency checks.

The numerical results support the asymptotic orders; the proof of the leading term is (17)-(19). The finite-coupling widths show that these particular completions do not preserve exact symmetric and antisymmetric degeneracies to all orders. No higher-order cancellation is assumed.

## 7. Second-order size bound and lower bound on perturbative order

Theorem. Under the assumptions in the question, the minimum achievable order for every \(N>d\geq2\) is 2. The achievable logical dimension at this order is unbounded for each fixed \(d\geq2\).

Proof. The given first-order rank bound excludes \(H^{(1)}=aI+bP_N\) with \(b\ne0\) when \(N>d\). Hence \(k\geq2\). Equations (1)-(19) provide \(k=2\) for every finite \(N\). Equation (21) is an alternative proof even for a product bond. There is therefore no finite upper bound on \(N\), and no universal function \(g(N,d)\) growing with \(N\) that lower-bounds \(k\). For \(d=1\), the bond is a scalar on the entire system and no connected term occurs at any order. This is a separate trivial obstruction.

The failure of port counting can be located in the resolvent. In the present construction,

\[
H^{(2)}=-\int_0^\infty dt\sum_{a,b}
R_{A,ab}(t)\otimes R_{B,ab}(t),
\]
\[
R_{A,ab}(t)=S_a^\dagger e^{-H_At}S_b
=\sum_e e^{-E_et}M_{A,e}^{ab},
\tag{24}
\]

and likewise for \(B\). At a fixed time there are only \(r^2\) product terms. Integration does not preserve that rank bound. The \(L\) distinct exponential functions provide additional independent logical directions, and the full-rank matrix \(1/(E_e+E_f)\) pairs them. This is already visible for the single-product example, which has just one product at each fixed time.

Equivalently, between the two port insertions, functions of the unrestricted \(H_A,H_B\) can act on arbitrarily large excited spaces. Such internal evolution costs no extra power of the inter-block parameter \(\lambda\). Treating it as if it were a fixed scalar denominator discards the resource used by the construction.

If additional small parameters were attached to internal block couplings, the logical exchange might be high order in those parameters while remaining second order in \(\lambda\). Those parameters are absent from the stated assumptions.

## 8. Resources and the locality question

For fixed \(d\), the explicit many-channel construction has

\[
L=1+\left\lceil\frac{N^2-1}{(d^2-1)^2-1}\right\rceil=O(N^2),
\qquad D=O(N^3),\qquad K=O(N^3).
\]

The chosen energies range from 1 to \(3^{L-1}\). This growing ratio is a choice, not a mathematical necessity: any distinct energies, even within a fixed interval such as \([1,3]\), make (8) strictly positive. Closely spaced energies make the whitening increasingly ill conditioned and can force a very small positive \(\tau\), hence a very weak exchange. No \(N\)-independent lower bound on \(b_2\) is proved or asserted. The theorem is an existence statement at finite \(N\), not an efficient, robust scaling claim.

Encoding a Hilbert space of dimension \(D=O(N^3)\) into qudits of fixed dimension needs only \(O(\log N)\) sites as a matter of storage capacity. That observation does not make (1) a nearest-neighbor Hamiltonian. Its spectral projectors and the unitary completion generally act on the whole block. The arbitrary positive-matrix dilation is the precise step where locality has been used up.

An approximate realization by a longer chain with ancillas and large penalty energies appears plausible through Hamiltonian simulation or a computation-clock construction. This is a proposal, not a proved reduction here. Such a reduction would need to establish all of the following together: an exactly \(N\)-dimensional ground space, an appropriate gap, the intended action of the physical last-site port, and the matrix-valued boundary resolvent (24). Matching a low-energy spectrum alone would not suffice.

For an exact chain construction, approximate simulation errors are particularly relevant. A nonscalar \(H^{(1)}\) error, however small, dominates the desired exchange as \(\lambda\to0\). If (4) is enforced exactly, this first-order problem is removed, but nonswap errors at second order still violate the exact claim unless canceled. No such cancellation scheme is assumed.

Thus the established resource statement is that the unrestricted-block construction keeps \(k=2\) while growing auxiliary dimension and spectral complexity. A local implementation could require substantially more than the storage-length lower bound, additional energy ratios, and precision that increases with \(N\). No quantitative local-chain upper bound, or exact local-chain realization, has been proved here. Deriving one would be a separate problem about boundary spectral engineering.

The main conclusion is therefore sharp for the stated setting: second-order antiferromagnetic swap is possible for arbitrarily large logical dimension through a fixed small port. The remaining directions are (i) an exact nearest-neighbor realization, and (ii) bounds on exchange strength or robustness when internal energy scales and chain interactions are constrained.
