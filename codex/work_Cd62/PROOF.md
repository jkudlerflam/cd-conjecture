# A length comparison in the ferromagnetic strip

## [proved] Statement

Put \(S^\mu=\sigma^\mu/2\). Let \(F_N(m)\) be the ground energy in the sector with \(m\) down spins of the periodic Hamiltonian

\[
 H_N^0=-J\sum_{x=1}^N
 (S_x^x S_{x+1}^x+S_x^y S_{x+1}^y+\cos\eta\,S_x^z S_{x+1}^z),
 \qquad J>0 .
\]

Use the sum of both oriented bonds at \(N=2\), as in the question. For every \(L\ge2\), \(0\le n\le L\), and \(0<\eta\le\pi/4\),

\[
 \boxed{\quad L\,|F_{2L}(2n)-2F_L(n)|\le\pi J\sin\eta .\quad}       \tag{1}
\]

This is the range \(-1<J_z/J\le-1/\sqrt{2}\). The coefficient in (1) tends to zero at the ferromagnetic endpoint. The proof compares the two finite rings without taking a thermodynamic limit.

## [proved] Finite Bethe equations and their ground state

Write \(c=\cos\eta\). Until the kernel estimate we allow \(0<\eta\le\pi/2\). Use the odd functions

\[
 p(x)=2\arctan(\tan(\eta/2)\tanh x),\qquad
 \theta(x)=2\arctan(\cot\eta\,\tanh x).
\]

Their derivatives divided by \(2\pi\) are

\[
 a(x)=\frac{\sin\eta}{\pi(\cosh 2x+c)},\qquad
 T(x)=\frac{\sin 2\eta}{\pi(\cosh 2x-\cos 2\eta)}.
\]

Here \(a>0\) and \(T\ge0\). For \(m\le N/2\) define the centered levels \(t_i^*=(i-(m+1)/2)/N\), \(1\le i\le m\). The counting map and one-particle energy are

\[
 f_i(y)=\frac{p(y_i)}{2\pi}+\frac{1}{2\pi N}\sum_{j=1}^m\theta(y_i-y_j),
 \qquad
 \epsilon(x)=J(c-\cos p(x))=-\pi J\sin\eta\,a(x).
                                                               \tag{2}
\]

We first establish existence also for the levels used between the two rings. Suppose the levels \(t_i\) are ordered, symmetric about zero, and the sum of their largest \(r\) elements is at most \(r(m-r)/(2N)\) for every \(1\le r\le m\). The centered levels meet these bounds with equality. The system \(f(y)=t\) has a unique real solution. Its Jacobian \(B\) acts by

\[
 (B u)_i=a(y_i)u_i+\frac1N\sum_j T(y_i-y_j)(u_i-u_j).
                                                               \tag{3}
\]

This symmetric matrix is positive definite by its quadratic form. A potential with gradient \(f-t\) is the sum of the primitives of \(p/(2\pi)\) at \(y_i\), plus \(1/N\) times the sum of the primitives of \(\theta/(2\pi)\) at \(y_i-y_j\) over \(i<j\), minus \(\sum_i t_i y_i\). All primitives start at zero. It is coercive, as follows. Up to a bounded term for fixed \(N,m,\eta\), its first two terms are

\[
 \alpha\sum_i|y_i|+\frac{b}{N}\sum_{i<j}|y_i-y_j|,
 \qquad \alpha=\frac{\eta}{2\pi},\qquad b=\frac12-\frac{\eta}{\pi}.
\]

Indeed the two odd functions approach their limiting values exponentially. Sorting the \(y_i\) in the same order as the \(t_i\) maximizes the subtracted linear term. At a positive height with \(r\ge1\) coordinates above it, the contribution per unit height to the resulting piecewise linear potential is at least

\[
 \alpha r+\frac bN r(m-r)-\frac{r(m-r)}{2N}
 =\alpha r\left(1-\frac{2(m-r)}N\right)>0 .
\]

The same bound holds at negative heights by symmetry of the levels. Integrating these inequalities gives a positive multiple of \(\max_i|y_i|\), since there are only \(m\) possible nonzero values of \(r\). Strict convexity now gives existence and uniqueness. Equation (2) orders the roots in the order of the levels, including ties, since \(p\) is strictly increasing and \(\theta\) is nondecreasing. The implicit function theorem gives smooth dependence on the levels along any of the paths below and real analytic dependence on \(\eta\) in the open interval. The Jacobian remains invertible at \(\pi/2\).

For the centered levels, the sector ground energy is

\[
 F_N(m)=-NJc/4+\sum_{i=1}^m\epsilon(y_i).                      \tag{4}
\]

Here is a verification of the ground-state identification. In the spin-position basis write a plane wave for ordered positions \(z_1<\dots<z_m\) as a sum over permutations of the momenta \(p(y_i)\). The hopping matrix element is \(-J/2\). Relative to \(-NJc/4\) the diagonal is \(Jcm\) minus \(Jc\) per occupied nearest-neighbour pair. In a free plane-wave equation the two forbidden hops at adjacent positions \(z,z+1\) give the matching condition

\[
 \psi(\ldots,z,z,\ldots)+\psi(\ldots,z+1,z+1,\ldots)
     =2c\,\psi(\ldots,z,z+1,\ldots).
\]

This accounts additively for each occupied bond, including consecutive occupied bonds. With \(A=e^{ip(u)}\), \(D=e^{ip(v)}\), exchanging two adjacent momenta multiplies the coefficient by

\[
 -\frac{1+AD-2cD}{1+AD-2cA}=-e^{i\theta(u-v)} .
\]

For the equality, substitute \(A=(1+i\tan(\eta/2)\tanh u)/(1-i\tan(\eta/2)\tanh u)\) and use the subtraction formula for \(\tanh\). The ratios are scalars and reversing an exchange inverts its ratio, so they consistently define all coefficients starting with coefficient one. Taking \(u=y_i\) once around the ring imposes

\[
 e^{iNp(y_i)}=(-1)^{m-1}\exp\!\left(-i\sum_j\theta(y_i-y_j)\right),
\]

which holds at the centered levels. Thus wherever this vector is nonzero it is a periodic eigenvector of energy (4).

At \(\eta=\pi/2\) the wavefunction is a determinant at consecutive momenta \(2\pi t_i^*\). The Vandermonde formula makes it a fixed nonzero constant times \(\prod_{k<l}\sin(\pi(z_l-z_k)/N)\), positive after one choice of phase. The spin configuration graph is connected for \(0<m<N\), since adjacent particle-hole swaps connect all configurations. All its hopping matrix elements are negative. The Perron-Frobenius theorem applied to a scalar shift of \(-H_N^0\) makes this eigenvalue the simple ground energy. The same simplicity holds throughout the interval. Near the free endpoint the Bethe vector stays nonzero by continuity. The characteristic polynomial of the sector matrix evaluated at the expression (4) is real analytic in \(\eta\), so it vanishes on the entire connected interval. Consequently (4) is an eigenvalue everywhere, even if a particular plane-wave formula were to vanish at an isolated parameter. A continuous eigenvalue coinciding with a simple ground eigenvalue continues to coincide locally, and the coincidence set is closed. This proves (4). The empty sector is immediate. Spin reversal gives \(F_N(m)=F_N(N-m)\).

## [proved] A total variation inequality for the Gaudin equation

The following argument works at arbitrary real locations \(y_j\), including repeated locations. Write \(\mu=N^{-1}\sum_j\delta_{y_j}\), \(k(x,y)=T(x-y)/a(x)\), and \(\bar k(x)=\int k(x,y)\,d\mu(y)\). For the ordered energy derivative in (2), set

\[
 w(x)=\epsilon'(x)/a(x)=2\pi J\sin p(x),\qquad
 v_i=(B^{-1}\epsilon')_i,
\]

where the vector on the right has entries \(\epsilon'(y_i)\) before multiplication by \(B^{-1}\). Extend the samples by

\[
 v(x)=\frac{w(x)+\int k(x,y)v(y)\,d\mu(y)}{1+\bar k(x)}.      \tag{5}
\]

On the right use the solved samples \(v_j\). Equation (3) ensures \(v(y_j)=v_j\), also at ties. For fixed finite data the functions in (5) are smooth with finite limits at both ends. They are rational functions of \(e^{2x}\) there, with removable endpoint limits, so their derivatives are integrable. All the \(k(\cdot,y_j)\) are bounded for these data.

A sufficient condition for \({\rm Var}(v)\le{\rm Var}(w)\) is the pair of half-line inequalities

\[
 \begin{split}
 \int_s^\infty |\partial_x k(x,y)|\,dx&\le k(s,y) \qquad (s\ge y),\\
 \int_{-\infty}^s |\partial_x k(x,y)|\,dx&\le k(s,y) \qquad (s\le y).
 \end{split}                                                    \tag{6}
\]

Differentiate (5) after multiplying by \(1+\bar k\). The result is

\[
 (1+\bar k(x))v'(x)=w'(x)+\int\partial_x k(x,y)(v(y)-v(x))\,d\mu(y).
\]

Take absolute values, bound \(|v(y)-v(x)|\) by the integral of \(|v'(s)|\) between its arguments, and integrate over \(x\). For fixed \(s>y\) the \(x\) integral runs from \(s\) to infinity. For fixed \(s<y\) it runs from minus infinity to \(s\). The single point \(s=y\) does not affect the integral in \(s\). Thus (6) bounds the double-integral contribution by \(\int\bar k(s)|v'(s)|\,ds\). Subtracting it from both sides proves the assertion. In particular, a large density of rapidities costs nothing in this estimate.

## [proved] The half-line inequalities in the strip

Take \(0<\eta\le\pi/4\), so \(C=\cos 2\eta\) belongs to \([0,1)\). For the right half-line put \(u=2(x-y)\ge0\), \(P=\cosh 2y\) and \(Q=\sinh 2y\). Apart from the positive factor \(2c\), \(k(x,y)\) is

\[
 g(u)=\frac{P\cosh u+Q\sinh u+c}{\cosh u-C},\qquad
 g(\infty)=P+Q .
\]

The derivative has the sign of \(Q(1-C\cosh u)-(CP+c)\sinh u\). For \(Q<0\) this is at most \(-|Q|+C|Q|e^{-u}-c\sinh u<0\), since \(P\ge |Q|\). For \(Q=0\) the function is nonincreasing as well. In these cases its right-tail variation is at most its starting value. For \(Q>0\) the derivative numerator is strictly decreasing, starts positive and ends negative. There is a single maximum.

It suffices to show \(2(g(u)-g(0))\le P+Q\) at every \(u\). Indeed, for a starting point before the maximum the remaining variation is twice the maximum minus the starting and limiting values. After the maximum it is simply the starting value minus the limit. To prove the inequality set \(z=\tanh(u/2)\). Direct subtraction gives

\[
 g(u)-g(0)
 =\frac{2Qz-2(CP+c)z^2/(1-C)}{(1-C)+(1+C)z^2}.
\]

After clearing the denominator, the required nonnegative expression is

\[
 (P+Q)((1-C)+(1+C)z^2)-4Qz+4(CP+c)z^2/(1-C).
\]

Using \(P\ge Q>0\), it is at least

\[
 2Q\big[(1-C)-2z+(1+C+2C/(1-C))z^2\big]\ge0.
\]

The last inequality follows from \((1-C)(1+C+2C/(1-C))=1+2C-C^2\ge1\). This proves the right half of (6). The identity \(k(-x,-y)=k(x,y)\) proves the left half. Finally \(p\) increases from \(-\eta\) to \(\eta\), so \(w\) is increasing with variation \(4\pi J\sin\eta\). We have proved

\[
 {\rm Var}_{\mathbb R}(v)\le 4\pi J\sin\eta .                   \tag{7}
\]

The range endpoint has a concrete explanation for this proof. For \(\pi/4<\eta<\pi/2\), the right half of (6) at \(s=y\) would require \(k(+\infty,y)\le2k(y,y)\). The ratio of the left side of this necessary condition to the right side tends to \(1-C>1\) as \(y\to+\infty\), using the same formula for \(g\). Thus the kernel condition itself fails there.

## [proved] Splitting each quantum number

By spin reversal take \(n\le L/2\). On the \(2L\) system take \(m=2n\) and use the levels

\[
 t_{2j-1}(s)=t_j^{\rm old}-s/(4L),\qquad
 t_{2j}(s)=t_j^{\rm old}+s/(4L),\qquad
 t_j^{\rm old}=\frac{j-(n+1)/2}{L},\qquad 0\le s\le1 .
\]

At \(s=1\) these are the centered levels for \(2L,2n\). At any smaller \(s\) they are obtained by moving each pair toward its mean. Sums of the largest \(r\) elements either contain whole pairs or end with a single upper member, and therefore obey the coercivity bound above. At \(s=0\) the unique solution consists of two copies of each old root, since doubling the number of terms in (2) cancels doubling \(N\). It need only solve the counting equations at this endpoint.

Evaluate the expression on the right of (4) along the path and call it \(\mathcal E(s)\). Its endpoints are \(2F_L(n)\) and \(F_{2L}(2n)\). Differentiating \(f(y(s))=t(s)\), using the symmetry of \(B\), gives

\[
 \mathcal E'(s)
 =(\epsilon')^{\mathsf T}B^{-1}t'(s)
 =\frac{1}{4L}\sum_{j=1}^n\big(v(y_{2j})-v(y_{2j-1})\big).
\]

The paired intervals have disjoint interiors. The absolute value of the sum is at most the variation in (7). Integrating over \(0\le s\le1\) gives (1), including odd lengths and arbitrarily dilute sectors. For \(n=0\) both sides of the energy comparison are equal.

## [proved] Field, hopping phase and the optimized bond scale

Choose a basis diagonalizing the nontrivial on-site qubit charge. Up to gauge the invariant bonds are exactly the family

\[
 h=J_\perp(S^x\otimes S^x+S^y\otimes S^y)
   +D(S^x\otimes S^y-S^y\otimes S^x)+J_zS^z\otimes S^z
   -\tfrac{B_z}{2}(S^z\otimes1+1\otimes S^z).
\]

Indeed the charge sectors have dimensions one, two and one. The four diagonal entries give a scalar, two individual \(z\) fields and \(J_z\), and the difference of the fields is gauge. The sole complex off-diagonal element gives \(J_\perp,D\). Put \(J=\sqrt{J_\perp^2+D^2}>0\) and \(J_z=-J\cos\eta\). The forward hopping matrix element has the form \(-Je^{i\alpha}/2\) for a fixed \(\alpha\), the backward element its conjugate. Site-dependent rotations about \(z\) put the hopping phase on an \(N\) ring into uniform increments \(\Phi_N/N\), with \(e^{i\Phi_N}=e^{iN\alpha}\) and \(|\Phi_N|\le\pi\). This follows by assigning the successive differences of the site phases, whose sum vanishes modulo \(2\pi\).

In sector \(m\) remove the field energy \(B_z(m-N/2)\), and denote the remaining energy by \(F_N^{\Phi_N}(m)\). Replacing wavefunction entries by their absolute values gives the lower inequality below. For the upper one use the real ground state at zero phase as a trial state. Its imaginary hopping expectation vanishes, and the real hopping operator on each bond has norm \(J/2\). Hence

\[
 0\le F_N^{\Phi_N}(m)-F_N(m)
 \le \frac{NJ}{2}(1-\cos(\Phi_N/N))\le\frac{J\pi^2}{4N}.
\]

Let \(n\) minimize the actual field and phase dependent energy at \(L\). Use sector \(2n\) at \(2L\). The field energies cancel exactly. For the short ring use the lower phase inequality and for the long ring the upper one. Then

\[
 L(E_{2L}-2E_L)\le \pi J\sin\eta+J\pi^2/8.
\]

Compress any gauge representative of the two-site bond to the span of the two opposite-spin basis vectors. The gauge additions leave its off-diagonal entry unchanged, since those vectors differ at both sites. Its magnitude is \(J/2\). A Hermitian two-by-two matrix with this entry has eigenvalue spread at least \(J\), and the full spread dominates that of the compression by the variational principle. Thus \(M\ge J/2\), proving

\[
 \boxed{\quad L(E_{2L}-2E_L)/M\le 2\pi\sin\eta+\pi^2/4 .\quad}   \tag{8}
\]

At prescribed vacuum down-spin density \(\rho\), one precise exact-compatibility convention for this class is \(L\rho\in\mathbb Z\) and \(e^{iL\alpha}=1\). This allows the vacuum charge and introduces no hopping twist on either ring. For example it selects even rings at zero phase and half filling. Equation (8) holds on the larger test set of all \(L\ge2\), including when exact compatibility selects no length for a given parameter. No conformal expansion is required.

For this same strip (8) also gives the one-sided Casimir variant at all lengths with \(K=(2/3)(2\pi\sin\eta+\pi^2/4)\). To see this, put \(C=2\pi\sin\eta+\pi^2/4\) and telescope \(E_{2^kL}/2^k\). Its \(k\)-th increment is at most \(CM/(2L\,4^k)\). The limit equals \(Le_0\). Indeed, the open-chain energies are subadditive up to the norm of one bond by product gluing, with a linear lower bound. Adding that norm to each open energy gives a subadditive sequence, so the energy densities converge. Closing the ring changes its energy by at most one bond norm. Summing the geometric series gives the stated \(K\).

## [assumed] Completing the uniform U(1) qubit class with supplied results

Take \(\delta=1-1/\sqrt 2\) in the away-from-endpoint bound of the supplied report, section 9.94. Combining (8) with the bounds in sections 9.91, 9.94, 9.104 and 9.106 gives the finite class constant

\[
 \max\{2\pi+\pi^2/4,\,6\cdot10^{15}(1-1/\sqrt 2)^{-25/3},
       \,4078,\,10781\}.
\]

The quoted ranges are respectively supplemented by zero defect at \(J_z/J\le-1\), the away bound from \(-1+\delta\) to zero, \(4078\) from zero to one, and \(10781\) above one. At \(J=0\) the bond is diagonal after removing gauge terms. Repeating a minimizing classical ring configuration shows \(E_{2L}\le2E_L\). The use of the supplied theorems on the complementary ranges is explicit here. Only the strip proof was derived above.

# A velocity budget for quadratic exchange through isolated cells

## [proved] Scope of the cell theorem

Consider identical cells of \(n\ge2\) consecutive sites of local dimension \(d\). The unperturbed Hamiltonian is a sum of independent open-cell Hamiltonians. Each has one nearest-neighbour term per internal bond, of norm at most \(J_*\) up to scalar energies. Subtract the cell ground energy, denote the resulting nonnegative operator by \(H_c\), its ground projection by \(P_c\), and suppose its excited energies are at least \(\Delta>0\). A bond \(V_i\) of norm \(g\le\min(\Delta,J_*)\) couples the two sites at interface \(i\). The bonds are translates by one cell. Write \(r\le d^2\) for the minimum number of product operators in \(V_i\) across that interface. The bound follows by expanding in any basis of single-site operators on the left.

On two cells put \(H=H_c\otimes1+1\otimes H_c\), \(P=P_c\otimes P_c\), and let \(H^+\) be the inverse on \(1-P\), zero on \(P\). Either value of \(\sigma\) is allowed in the following effective Hamiltonian.

\[
 K=\sum_i W_i,\qquad W_i=\sigma P V_iP-PV_i H^+V_i P,\qquad
 \sigma\in\{0,1\}.                                               \tag{9}
\]

Each \(W_i\) acts on the two indicated cell ground spaces. Assume that one fixed identification of each ground space with the Fock space of \(m\) complex modes makes each \(W_i\) a number-conserving quadratic operator plus a scalar. Nearest-cell hopping includes the usual parity on the left cell in the ordinary tensor product. The Fermi energy can be chosen arbitrarily. Equation (9) covers the second-order exchange through a virtual pair when only a second action on the same interface can return the system to the ground space. For an application of the theorem, (9) and the quadratic identification both need verification.

At each linear Fermi crossing let \(u_j\) be the absolute group velocity in original lattice sites per unit time, counting left and right movers separately. There is a universal explicit \(C_*\) such that

\[
 \sum_j u_j\le C_* d^4 J_* .                                    \tag{10}
\]

The estimate is independent of \(n,m,\Delta\), and of the largest excited energy of a cell. There is also a velocity-resolved bound growing at most as \(d^4[1+\log_+^2(J_*/u)]\) below a universal multiple of \(J_*\), where \(\log_+ t=\max(0,\log t)\). For linear complex-fermion sectors the free-chain formula in section 9.71 of the supplied report identifies the effective Casimir coefficient with half the sum in (10). If the cell construction has \(J_*\) bounded by a specified multiple of the microscopic \(M\), (10) gives the corresponding budget in units of \(M\).

The two estimates needed for (10) are a bound on the number of fast crossings and a bound on the largest speed when the interfaces are far apart on the scale \(J_*/\Delta\).

## [proved] A resolvent approximation and a count of crossings

Set \(\delta=\min(\Delta,J_*)\), \(x=n\delta/J_*\), and \(\lambda=\sigma g+g^2/\delta\). We can take \(g>0\). For each integer \(h\ge1\), the second-order term in (9) has an approximation of operator Schmidt rank at most \(r^2(4h^2+1)\) and operator norm error at most \(3(g^2/\delta)2^{-h}\).

To see this, approximate \(1/z=\int_0^\infty e^{-tz}\,dt\) for \(z\ge1\). Discard \(t<2^{-h}\) and \(t>b_h=2^{\lceil\log_2(h+1)\rceil}\), losing at most \(2\cdot 2^{-h}\). The retained range consists of at most \(2h\) dyadic intervals. On each \([b,2b]\) approximate \(e^{-tz}\) in \(t\) by its Taylor polynomial of degree \(2h-1\) at \(3b/2\). The pointwise remainder is at most \(2^{-2h}\), since \(\sup_{z\ge0}z^k e^{-bz}=(k/(eb))^k\) and \(k!\ge(k/e)^k\). The integrated error is at most \(b_h2^{-2h}\le2^{-h}\). Write \(f_h(z)\) for the resulting integral approximation.

For \(z=H/\delta\), each integrated polynomial is \(e^{-3bz/2}\) times a polynomial of degree at most \(2h-1\) in \(z\). In its expansion in \(H_c/\delta\) on the two cells, collect terms by the power on the left. There are at most \(2h\) products. Thus \(f_h(H/\delta)-f_h(0)P\) has rank at most \(4h^2+1\) across the interface. On division by \(\delta\) its error from \(H^+\) has the asserted scalar bound by the spectral theorem. Multiplication on each side by \(V_i\) multiplies the rank bound by \(r^2\), and compression by the product \(P\) costs no rank. Adding the first-order term costs at most \(r\) more.

Let \(T_0\) be the \(m\)-by-\(m\) nearest-cell hopping matrix in \(K\). The quadratic hopping coefficients can be extracted in normalized Hilbert-Schmidt inner product on the two-cell Fock space. On the left use the orthonormal operators \(\sqrt2\,a_j^\dagger(-1)^{\hat n}\), on the right \(\sqrt2\,a_j\), with \(a_j\) the cell annihilation operators. Orthonormality follows from \(2^{-m}{\rm tr}(a_j^\dagger a_k)=\delta_{jk}/2\) and unitarity of parity. The coefficient submatrix is \(T_0/2\). Terms on one cell and the conjugate hopping term are orthogonal to that block. A sum of \(R\) products has a coefficient submatrix of rank at most \(R\). The Frobenius norm of the extracted error is at most the normalized Hilbert-Schmidt norm of the operator error, hence at most its operator norm. In decreasing order of singular values this proves

\[
 s_{R+1}(T_0)\le6\lambda\,2^{-h},\qquad R=r+r^2(4h^2+1),       \tag{11}
\]

with missing singular values taken as zero. Full exact rank at second order is consistent with (11).

The one-cell Bloch matrix has the form \(H(k)=H_b+e^{ik}T_0+e^{-ik}T_0^\dagger\). Subtract the Fermi energy from \(H_b\). For a normalized zero eigenvector \(q_j\) at \(k_j\) write \(z_j=e^{ik_j}\) and set

\[
 X_j=(q_j,z_jq_j),\qquad
 \mathcal J=\begin{pmatrix}0&iT_0\\-iT_0^\dagger&0\end{pmatrix}.
\]

Then \(X_j^\dagger\mathcal J X_j=\epsilon'_j(k_j)\) and \(\|X_j\|^2=2\). Distinct crossings are \(\mathcal J\)-orthogonal. For different phases use \(q_i^\dagger(H(k_j)-H(k_i))q_j=0\) and divide by \(z_j-z_i\). At equal phase diagonalize the compression of \(H'(k)\) to the zero eigenspace. Its eigenvalues are the first-order slopes, as follows by expanding the eigenvalue equation to first order and projecting to that eigenspace.

Let \(p\) columns have positive slope at least \(s>0\), and form their column matrix \(X\). Then \(X^\dagger\mathcal J X\ge sI\). Project onto the positive eigenspaces of \(\mathcal J\) with eigenvalues at least \(s/4\), using a projection \(D\) of rank \(a\). These positive eigenvalues are the singular values of \(T_0\). With \(Q\) the orthogonal projection onto \(\ker(DX)\) in column space,

\[
 s\,{\rm rank}(Q)\le {\rm tr}(QX^\dagger\mathcal J XQ)
 \le (s/4){\rm tr}(X^\dagger X) = sp/2.
\]

Since \({\rm rank}(Q)\ge p-a\), there are at most \(2a\) such positive crossings. Applying the argument to \(-\mathcal J\) gives at most \(4a\) crossings in both directions. Physical speeds are \(n|\epsilon'_j|\). Therefore, with \(N_F(u)=\#\{j:u_j\ge u\}\),

\[
 N_F(u)\le (16h^2+8)r^2,\qquad
 h=\max\{1,\lceil\log_2(48n\lambda/u)\rceil\}.                 \tag{12}
\]

Indeed (11) then bounds \(s_{R+1}\) strictly below \(u/(4n)\).

## [proved] Commutation of distant projected interfaces

On three cells let \(P_{\rm all}\) be their product ground projection and \(\tau_t\) the dynamics of their decoupled Hamiltonian. We construct a representative of each \(W_i\) that commutes with \(P_{\rm all}\) on the full Hilbert space. This step matters because compression alone need not preserve commutators of separated operators.

Choose a smooth even cutoff \(\chi(\omega)\) equal to one for \(|\omega|\le1/2\) and zero for \(|\omega|\ge1\). An explicit choice is \(b(1-|\omega|)/(b(1-|\omega|)+b(|\omega|-1/2))\) with \(b(t)=e^{-1/t}\) for \(t>0\) and zero otherwise. Let \(\ell,f\) have Fourier transforms \(\chi\) and \((1-\chi(\omega))/|\omega|\), respectively, defining the latter as zero near zero. Use \(\widehat f(\omega)=\int f(t)e^{it\omega}\,dt\). Both inverse transforms are real and even and

\[
 B_f=1+\int_{\mathbb R}(1+|t|)^8(|\ell(t)|+|f(t)|)\,dt<\infty .
\]

Indeed both frequency functions are square integrable, so their inverse transforms are integrable on bounded intervals by Plancherel and Cauchy-Schwarz. Their tenth derivatives are integrable, giving \(O(|t|^{-10})\) decay away from zero by integration by parts.

For \(V=V_i\) write \(R_V=\delta^{-1}\int f(t)\tau_{t/\delta}(V)\,dt\). Between a ground state and an excited state of energy \(E\) this multiplies the matrix element by \(1/E\), and its ground-to-ground matrix elements vanish. Consequently

\[
 D_i=\int\ell(t)\tau_{t/\delta}\big(\sigma V-(V R_V+R_V V)/2\big)\,dt
\]

is self-adjoint, commutes with \(P_{\rm all}\), and restricts there to \(W_i\). The last filtering multiplies ground-to-excited elements by zero and ground-to-ground elements by one. Because the cells are independent under \(\tau\), the inverse energy seen after applying \(V_i\) to a ground state is exactly the two-cell one in (9).

For clarity, the locality estimate used for two neighbouring interfaces, whose bond indices are \(n\) apart, is

\[
 \|[\tau_t(V_i),V_{i+1}]\|
 \le 2g^2\min\{1,e^{\,1-n+6eJ_*|t|}\}.                        \tag{13}
\]

Here is a path derivation. For fixed \(B\) let \(C_B(X,t)\) be the supremum of \(\|[\tau_t(A),B]\|/\|A\|\) over operators supported on \(X\). Its bound is \(C_B(X,0)\) plus \(2\sum_{Z:Z\cap X\ne\varnothing}\|h_Z\|\int_0^{|t|}C_B(Z,s)\,ds\), with the sign of time reversed if needed. Indeed only those \(h_Z\) enter \([H,A]\). Differentiating the commutator and applying the Jacobi identity gives a homogeneous commutator with \(\tau_t(\sum_Z h_Z)\), removable by unitary conjugation, and source terms bounded by the stated integrand. Iteration from a bond makes at most three choices per step, including the current bond. To end on a bond touching the other interface needs at least \(n-1\) steps. The iterated integrals give the exponential series with terms \((6J_*|t|)^j/j!\) for \(j\ge n-1\). Multiplying each term by \(e^{j-(n-1)}\ge1\) proves (13). The iterated remainder tends to zero by the same factorial bound. Scalars in the cell energies have no effect.

Each term in \(D_i\) is an integral of either one translated \(V_i\) or two translated \(V_i\)'s. In dimensionless time their arguments are \(t,t+s\), with weights \(\ell(t),f(s)\) as applicable. In a commutator of two such products there are at most four single commutators by the product rule. When each of the at most four integration variables has absolute value at most \(x/(48e)\), every physical time difference is at most \(n/(12eJ_*)\). Equation (13) contributes at most the factor \(e^{1-n/2}\). On the remaining region one of the integrals has tail mass at most \(B_f(48e)^8/(1+x)^8\). All other integral norms are at most \(B_f\). Using \(x\le n\) and \((1+n)^8 e^{-n/2}\le16^8\), we obtain

\[
 \|[W_i,W_{i+1}]\|\le\|[D_i,D_{i+1}]\|\le
 \gamma=\frac{C_8\lambda^2}{(1+x)^8},\qquad
 C_8=8 B_f^4\big(e\,16^8+4(48e)^8\big).                       \tag{14}
\]

The first inequality uses the exact invariant ground space for both representatives. In the second, a product of \(p\) factors from one interface and \(q\) from the other costs at most \(2pq\,g^{p+q}\) before denominators and filter norms. Here \(p,q\le2\), the two orders of the anticommutator have total weight one, and the coefficients sum to \(\lambda^2\). This verifies the constant also when the first-order term is present.

## [proved] Largest speed from the interface commutator

Subtracting its vacuum scalar, restrict each quadratic \(W_i\) to the one-particle space to obtain \(w_i\), supported on cells \(i,i+1\). The commutator norm for neighbours is at most \(\gamma\) by (14). Group two cells for Bloch analysis. The even-bond sum has a constant Bloch matrix \(a\). The odd-bond sum has \(b(k)=D(k)b(0)D(k)^\dagger\) where \(D(k)\) puts the phase \(e^{ik}\) on one subcell and one on the other, up to the sign convention for \(k\). This holds since the odd bonds also pair disjoint cells. Thus \(\|\partial_k e^{-itb(k)}\|\le2\) for all real \(t\).

The commutator \(C(k)=[a,b(k)]\) has norm at most \(4\gamma\). In position space it is the sum of \([w_{2j},w_{2j\pm1}]\). At fixed sign and fixed parity of \(j\), the supporting triples are disjoint and form a direct sum in the one-particle space. This gives four sums of norm at most \(\gamma\). Since \(C(k)\) has Fourier degree one, its two nonconstant Fourier coefficients each have norm at most \(4\gamma\) by integration. In particular \(\|C'(k)\|\le8\gamma\).

Write \(U=e^{-it(a+b(k))}=U_aU_b V\), with \(U_a=e^{-ita}\), \(U_b=e^{-itb(k)}\). The generator of \(V\) is

\[
 G(t,k)=U_b^\dagger(U_a^\dagger b(k) U_a-b(k))U_b
       =U_b^\dagger\left(\int_0^t i e^{isa} C(k)e^{-isa}\,ds\right)U_b.
\]

For \(t\ge0\), differentiation in \(k\) and the preceding norm bounds give \(\|\partial_k G(t,k)\|\le24\gamma t\). Integrating the differentiated unitary equation gives \(\|\partial_k V(t,k)\|\le12\gamma t^2\), so \(\|\partial_k U\|\le2+12\gamma t^2\). On an eigenvector of \(a+b(k)\) the diagonal matrix element of \(U^\dagger\partial_k U\) is \(-it\) times that of \(a'+b'\), by differentiating the exponential under its time integral. At degeneracy one uses the slope eigenvectors as before. Taking \(t=\gamma^{-1/2}\) bounds each two-cell slope by \(14\sqrt\gamma\). If the commutator vanishes one may take \(t\) arbitrarily large. The two-cell wave number is twice the one-cell wave number modulo \(2\pi\). Thus in physical units every speed obeys

\[
 u_j\le U_*=\frac{B n\lambda}{(1+x)^4},\qquad B=28\sqrt{C_8}.   \tag{15}
\]

## [proved] Summation and boundary of the argument

For \(z=n\lambda/u\) the integer \(h\) of (12) obeys \(h^2\le256(1+\sqrt z)\). For \(z\le1\) use \(h\le 1+\log_2 48\) or \(h=1\). For \(z\ge1\) use \(\log z\le4z^{1/4}\). Equations (12) and (15) now yield

\[
 \sum_j u_j=\int_0^{U_*}N_F(u)\,du
 \le5000 r^2\big(U_*+2\sqrt{n\lambda U_*}\big)
 \le 10000 r^2 J_*(B+2\sqrt B),
\]

since \(n\lambda\le2J_*x\). This gives an explicit \(C_*\) for (10). To see the sharper resolved count, a nonempty count requires \(u\le 2BJ_*/(1+x)^3\) by (15). Substituting \(n\lambda/u\le(2J_*/u)(2BJ_*/u)^{1/3}\) in (12) proves the asserted squared-logarithm bound, including a velocity cutoff \(2BJ_*\).

This rules out an unbounded harmonic tower in the effective Hamiltonian (9) at a common bounded \(J_*/M\). The spectral approximation controls how many fast crossings one interface can carry. The commutator estimate controls transmission from that interface to the next one. An additional intracell term mixing the two ends would enter the even or odd Bloch matrix in the commutator argument and would require its own locality estimate. Resonances outside the retained ground space, general interacting payloads, additional multi-interface terms and higher orders are outside (9). In particular, a residual norm estimate of \(\varepsilon\) per cell gives by direct energy comparison only \(4\varepsilon L^2/n\) for the error in the scaled doubling defect on a physical ring of \(L\) sites divisible by \(n\). This grows even at fixed \(n\). A microscopic result needs a length-uniform estimate for that remainder or a direct ring argument.
