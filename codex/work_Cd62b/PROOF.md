# The velocity budget for cells beyond second order

**Neither (P) nor (D).** There is an all-order budget for the Hamiltonian obtained by solving each interface exactly in isolation, if the resulting chain is quadratic. There is a quadratic budget for the full expansion through *any fixed order*, including terms spanning several interfaces. For interacting truncations the argument bounds the largest speed. The hypothesis needed to sum the interacting speeds does not follow from the cell gap.

## 1. Conventions and the effective Hamiltonian [proved]

Write \(N=L/n\) for the number of cells and \(P=\bigotimes_c P_c\). Combine the terms on each physical internal bond before measuring \(J_*\), as in the supplied quadratic theorem. Subtract each cell ground energy. Put \(H_0=\sum_c H_c\), \(B=H_0/\Delta\), \(q=g/\Delta\), \(W_i=V_i/g\). The case \(g=0\) has no transmission. Write \(r\le d^2\) for a bound on the ordinary operator Schmidt rank of a port across its bond. Expanding in a matrix-unit basis on its left physical site gives this bound. One may use any positive lower bound on the true cell gap for \(\Delta\). The nontrivial cells have \(n\ge2\) and a nonempty excited space. Completely flat cells of length at least two give disjoint physical interface bonds when the full cell Hilbert space is retained. If one-site cells with on-site terms bounded by \(J_*\) are allowed, a nonempty cell excitation spectrum gives \(\Delta\le2J_*\). The short-cell estimates below then apply directly. For the velocity results use a periodic array of identical cells and translate the couplings by a cell. The local operator estimates allow nonidentical cells with the same bounds.

Use the sector decomposition when it gives independent cell spaces at zero port. The isolated cell dynamics and each retained port must act within this sector, so turning ports on separately has the stated meaning there. The physical Lieb–Robinson estimates below can be made on a local extension and then restricted to the invariant sector. A later restriction to a conserved total charge can be made after the tensor product calculation. An arbitrary nonlocal subspace with no such cell decomposition supplies no product resolvent. Periodic vacuum comparisons refer to complete periods with the infrared vacuum allowed, with physical speed in sites per unit time. The operator estimates do not require a choice of vacuum subsequence.

For a set \(X\) of \(b\) interfaces let \(W_X=\sum_{i\in X}W_i\). Unused cells can be omitted. In dimensionless energy units define
\[
 P_X(z)=\frac{1}{2\pi i}\oint_{|\zeta|=1/2}
       (\zeta-B-zW_X)^{-1}\,d\zeta,\qquad
 U'_X=[P'_X,P_X]U_X,\quad U_X(0)=1 .
\]
Define \(K_X(z)=P U_X^{-1}(B+zW_X)U_X P\), with \(K_\varnothing=0\). The dimensional energy is \(\Delta K_X(q)\). This fixes the ground-space basis to all orders. From \(P_X^2=P_X\) one gets \([[P'_X,P_X],P_X]=P'_X\), so \(U_X P U_X^{-1}=P_X\). The generator has zero diagonal blocks relative to \(P_X\). For real \(z\), \(U_X\) is unitary. Differentiation, using \([P_X,B+zW_X]=0\), gives
\[
 K_X'(z)=P U_X^{-1}W_X U_X P .                                      \tag{1}
\]
These facts only need an analytic germ at zero. For a finite cluster the resolvent series converges for \(2b|z|<1\). On \(|z|\le1/(8b)\), the contour expansion gives \(\|P_X\|\le4/3,\ \|P'_X\|\le 2b/(1-1/4)^2\). Integration of the transport equation on a radius bounds both transport norms by \(e^2\). Thus
\[
 \|[z^k]K_X\|\le e^4 b(8b)^{k-1},\qquad k\ge1.                     \tag{2}
\]
Indeed \(K_X\) on the circle is bounded by \(e^4 b|z|\) by (1). For real \(|z|<1/(2b)\), (1) also gives \(\|K_X(z)\|\le b|z|\).

Define the contribution with precisely the interface set \(X\) by
\[
 C_X(z)=\sum_{Y\subseteq X}(-1)^{b-|Y|} K_Y(z).
\]
On disconnected components touching disjoint cells the low projection factorizes for sufficiently small \(z\). The product of the two horizontal transports on the initial low space is horizontal, since each factor has zero low-to-low derivative. Uniqueness of the horizontal lift gives additivity of \(K\). Equivalently, on the low space the product solves the displayed transport equation. Consequently \(C_X=0\) for disconnected \(X\). Each monomial of degree \(k\) in the ports uses at most \(k\) labels. Inclusion and exclusion cancels a monomial missing any label of \(X\). Summing \(C_X\) recovers \(K\) by binomial inversion. This proves the linked expansion on every finite ring as an identity of Taylor coefficients. At order \(k\) the support is an interval of at most \(k\) interfaces, or the ring itself if it winds. Each contribution has norm at most
\[
 2^b e^4\,\Delta\, b(8b)^{k-1}q^k,\qquad b\le k .                 \tag{3}
\]
At most \(\sum_{b\le k}(b+1)\) such intervals touch a fixed cell, up to duplications which only reduce the count. A one-cell ring has a self-interface and can instead be treated directly by (2). Ranks across a two-sided cut below concern distinct cells. This construction takes the Taylor coefficient before increasing the ring length.

## 2. Resolvent words and ranks [proved]

Let \(R=B^+\), the inverse off the ground space, zero on it. In the residue for \([z^j]P_X(z)\) there are \(j\) ports and \(j+1\) unperturbed resolvents. Use
\[
 (\zeta-B)^{-1}=P/\zeta-\sum_{m\ge1}\zeta^{m-1}R^m.
\]
If \(m=0\) denotes choosing the pole, the residue condition is exactly \(\sum m=j\). There are \(\binom{2j}{j}\le4^j\) choices. Every resulting word has \(j\) factors \(W_X\), exactly \(j\) copies of \(R\), and some \(P\)'s, with coefficient of modulus one before aggregation.

For clarity, count words through the transport, including the normalization of the low states. Write \(U_c=[z^c]U_X\). Its recursion is
\[
 U_c=\sum_{a+l+t=c,\ a\ge1}(a/c)[P_a,P_l]U_t,\qquad P_a=[z^a]P_X.
\]
The inverse transport has the commutator on the other side and the opposite sign. In both recursions the scalar weight along each branch has modulus at most one. If \(p(x)=1/(1-4x)\), the two word counts have the coefficientwise majorant \(u(x)=1/[1-2(p(x)-1)p(x)]\), by this recursion. At \(x=1/32\), \(u(x)<2\). By (1), \([z^k]K_X\) therefore has at most \(4\cdot32^{k-1}\) uncollected words, each with coefficient at most one, \(k\) ports and exactly \(k-1\) reduced inverses. This avoids an independent sum over energy exponents at each inversion.

Here are two different rank estimates. The first is in operator norm. For any cut into two collections of cells, \(B=B_l+B_r\) is a tensor sum and \(P=P_l\otimes P_r\). There are operators \(R_h\) of Schmidt rank at most \(4h^2+1\) satisfying \(\|R-R_h\|\le3\,2^{-h}\), \(h\ge1\). The construction from the supplied quadratic proof also works on this larger cut. For \(E\ge1\) integrate \(e^{-tE}\) from \(2^{-h}\) to \(2^{\lceil\log_2(h+1)\rceil}\), losing at most \(2\,2^{-h}\). Use at most \(2h\) dyadic intervals \([a,2a]\). Taylor expansion of degree \(2h-1\) about \(3a/2\) has pointwise error at most \(2^{-2h}\) by the bound on \(E^{2h}e^{-aE}\) and \((2h)!\ge(2h/e)^{2h}\). The total integrated error is at most \(2^{-h}\). Each integrated polynomial times its exponential has at most \(2h\) products across the cut, by collecting powers of \(B_l\). Subtract the value of the approximation at zero times \(P\).

For example, a rank bound at order \(k\) for one interface is
\[
 4\cdot32^{k-1}r^k(4h^2+1)^{k-1}.
\]
The dimensional error is at most \(4\cdot32^{k-1}\Delta q^k (k-1)3\,2^{-h}(1+3\,2^{-h})^{k-2}\) for \(k\ge2\), by telescoping the word. The first order is exact. On a cluster replace \(r\) by \(r+2\) on an interval cut crossing one port, since all other ports sum to one term on each side. Replace \(q^k\) in the error by \((bq)^k\), and sum over subsets for \(C_X\). On a cut crossing two ports use \(2r+2\). Single-cell terms are obtained, for example, by taking normalized partial traces on the other ground spaces. These maps are contractions. On two cells, subtraction of the two one-cell parts and addition of the scalar costs at most four times the original norm. The unperturbed compressed cell energy was scalar.

The second estimate controls an entire series. Use normalized trace inner products on both operator spaces, so the coefficient matrix of an operator of norm \(w\) has Frobenius norm at most \(w\). Set \(R_0=0\), \(E_h=R_h-R_{h-1}\). Then
\[
 R=\sum_{h\ge1}E_h,\quad \|E_h\|\le9\,2^{-h},\quad
 {\rm Schmidt\ rank}(E_h)\le d_h=10(h+1)^2 .
\]
Replace the \(k-1\) inverses in every word by this telescopic series. The powers of \(h\) now enter a product of convergent one-variable sums.

For \(0<p\le1\) a coefficient submatrix of rank \(s\) and Frobenius norm \(w\) has \(\sum_j s_j^p\le s^{1-p/2}w^p\), where \(s_j\) are its singular values. These upper bounds add over different summands. Here is a proof of the addition step. Collect the individual singular value decompositions as \(T=\sum_j t_j a_j b_j^\dagger\) for unit vectors. The Hermitian dilation with off-diagonal block \(T\) is bounded above by \(\sum_j t_j u_j u_j^\dagger\), \(u_j=(a_j,b_j)/\sqrt2\), since the omitted part is negative semidefinite. Its positive eigenvalues are the singular values of \(T\). By min-max they are bounded by the corresponding eigenvalues of the upper bound, whose nonzero spectrum equals that of the Gram matrix of columns \(\sqrt{t_j}u_j\). Scalar concavity applied to the spectral decomposition at each diagonal entry of that Gram matrix bounds the trace of its \(p\)-th power by \(\sum_j t_j^p\). This proves the claim for finite sums. Limits follow by continuity in fixed finite dimension, with bounds independent of dimension.

Define
\[
 S_p=\sum_{h\ge1}d_h^{\,1-p/2}(9\,2^{-h})^p,\qquad
 b_p(q,r)=32 r^{\,1-p/2} S_p q^p .                              \tag{4}
\]
Across one isolated interface the coefficient matrix of its exact dimensional energy, or any submatrix in local orthonormal operator bases, satisfies
\[
 \sum_j s_j^p\ \le\ \frac{4r^{\,1-p/2}}{1-b_p(q,r)}\,g^p
 \quad\text{if }b_p(q,r)<1\text{ and }q<1/8.                    \tag{5}
\]
Indeed, at order \(k\) the sum over telescopic blocks is at most \(4r^{1-p/2}g^p b_p^{k-1}\). Apply the bound first to finite telescopes and finite Taylor sums, then take limits using (2). The same order-by-order argument on a \(b\)-port cluster uses \(r+2\) and \(b g\) on an open cut. It gives a finite bound for each fixed \(k,p\), irrespective of convergence of the majorant over \(k\). At \(p=1\), \(S_1=27\sqrt{10}\). Thus
\[
 q\le q_0(d)=\frac{1}{64\,(27\sqrt{10})\,d}                      \tag{6}
\]
gives \(b_1\le1/2\). By continuity some \(p(d)<1\) works with \(b_p\le3/4\) throughout this range. One explicit choice uses \(r=d^2,\ q=q_0(d)\), \(M_0=\max(2,2b_{1/2})\) and \(p=1-\log(3/2)/(2\log M_0)\). Log convexity of \(b_p\) follows from Hölder applied to (4). This improves the power of \(d\) in the suggested threshold for *this coefficient norm*. The factor \(r^{1-p/2}\) uses the matrix Frobenius bound. Equation (5) alone does not convert an interacting coefficient matrix into an infrared spectrum. For general operators one can also get an all-order approximation in operator norm at a smaller threshold. Use \(d_h\) in place of \(d_h^{1-p/2}\) and \(r\) in place of \(r^{1-p/2}\) in the counting. If the resulting sum of rank times dimensionless weight to the \(p\)-th power is \(B_0<\infty\), discard weights below \(t\). The error is at most \(g B_0 t^{1-p}\) and the retained rank at most \(B_0t^{-p}\), for \(p<1\).

### Exact second-order rank can grow without bound [proved]

Take \(d=4\), a static gate bit and a hard-core token per site. In a cell charge energy \(m>2\) per token and hop a token with amplitude \(-1\) only across two open gates. Give each interface the port \(g(\sigma^+\otimes\sigma^++\sigma^-\otimes\sigma^-)\) on its two tokens. Its norm is \(g\) and its Schmidt rank is two. The cell vacuum is exactly degenerate over all gate patterns. In the sector of \(t\) tokens the hopping row sum is at most \(2t\), so the gap is at least \(m-2\). All terms are bounded local terms.

Choose patterns for which the run of open gates at the port has length \(a\) on the left and \(b\) on the right, from 1 to \(n\). After pair creation the two walkers live on paths of those lengths. The diagonal code matrix elements of the second-order term, divided by \(-g^2\), contain the matrix
\[
 G_{ab}=\int_0^\infty f_a(t) f_b(t)\,dt,\qquad
 f_a(t)=\langle1|e^{-t(mI-\mathcal A_a)}|1\rangle
\]
with \(\mathcal A_a\) path adjacency. The slowest exponential in \(f_a\) has rate \(m-2\cos(\pi/(a+1))\), strictly decreasing with \(a\), and positive coefficient from the first sine eigenvector. Taking the slowest surviving exponential in any linear relation proves independence. The Gram matrix has rank \(n\). Evaluating a product operator on product diagonal states gives a rank-one matrix, hence the effective Schmidt rank is at least \(n\) at order two. This is a classical interaction on the code. It gives no propagating tower. It shows why the approximation or the weighted singular value estimate is essential even at fixed order.

### A connected two-interface term at fourth order [proved]

A short cell need not suppress a process using both ends. In each two-qubit cell fix the conserved sector \(X_LX_R=+1\). The two boundary operators \(X_L,X_R\) then act as the same flag Pauli \(X\), and \(Z_LZ_R\) acts as its anticommuting Pauli \(Z\). Take \(H_c=(1-Z_LZ_R)/2\) and the physical port \(z X_R X_L\) between neighbouring cells. A free spectator qubit on a site of each cell gives an exactly degenerate code if desired. On a three-cell interval the model reduces to three flags with energy one per flip and perturbation \(z(X_1X_2+X_2X_3)\).

For a single active port the energy is \(1-\sqrt{1+z^2}\). For both active ports the vacuum couples with amplitude \(\sqrt2z\) to \((|110\rangle+|011\rangle)/\sqrt2\), which couples with the same amplitude to \(|101\rangle\). The latter two states have unperturbed energy 2. Eliminating these two rows gives \(E=2z^2/[E-2-2z^2/(E-2)]\). Expansion at zero gives \(E=-z^2+z^6/4+O(z^8)\). Subtracting the two isolated energies leaves \(-z^4/4+3z^6/8+O(z^8)\). Thus even a scalar linked energy can be missed by summing the isolated answers. The script *check_resolvent.py* checks the coefficients by both contour residues and the eigenvector recursion in rational arithmetic. It also computes the first five gate-path Gram pivots above by solving the two-particle linear systems at mass three.

## 3. Transmission between the ends of a cell [proved]

Set \(x=n\min(1,\Delta/J_*)\). For any positive integer \(m\) let \(f_m=(1+x)^{-m}\). Constants written \(c_m\) can depend on \(m\) and fixed scalar filter functions. They are independent of \(n,D,N,d,g,\Delta,J_*\). Four consequences in the horizontal gauge are as follows. Here \(E_X=\Delta K_X(q)\) and \(E_i=E_{\{i\}}\).

1. If \(b g\le\Delta/8\), then \(\|E_X-\sum_{i\in X}E_i\|\le c_m b^2(g^2/\Delta) f_m\). This compares exact finite-cluster Hamiltonians.
2. If \(g\le\Delta/8\), then \(\|[E_i,E_{i+1}]\|\le c_m g^2 f_m\), to all orders on the isolated interfaces.
3. If \(b\ge2\), the right side of (3) for a fixed Taylor order can be replaced by \(c_{k,m}2^b \Delta (bq)^k f_m\), even if \(b g>\Delta/8\). For adjacent single-interface coefficients of orders \(j,k\), the commutator is at most \(c_{j,k,m}\Delta^2 q^{j+k} f_m\).
4. The same commutator bound applies to the two induced one-cell terms from opposite ends. For the normalized traces over the two outside cells, the commutator of the traced terms equals the trace over both outside cells of the original commutator. The middle factors retain their order. On a two-cell ring take independent copies of the outside cell when computing the two one-cell traces, which gives the same two operators. One-cell parts of a connected multi-interface term are bounded by its norm.

Here is the proof, with the real-time and analytic estimates separated. Take a real perturbation path \(H_s=H_0+s\sum_{i\in X}V_i,\ 0\le s\le1\). In the sector under consideration the low energies have width at most \(2b g\), and the separation from the rest is at least \(\Delta-2b g\), by min-max. Choose an even smooth cutoff \(\chi\) equal to one for \(|\omega|\le1/4\) and zero for \(|\omega|\ge1/2\). For an explicit choice interpolate with \(a(t)=e^{-1/t}\) for \(t>0\), zero otherwise, using \(a(1/2-|\omega|)/(a(1/2-|\omega|)+a(|\omega|-1/4))\). Take scalar kernels with
\[
 \widehat F(\omega)=-(1-\chi(\omega))/\omega,\qquad
 \widehat l(\omega)=\chi(\omega),\qquad
 \widehat F(\omega)=\int F(t)e^{it\omega}\,dt .
\]
The first transform is zero near zero. All absolute moments of both kernels are finite. Their transforms are square integrable, and arbitrarily high positive derivatives are integrable. Plancherel on a bounded time interval and integration by parts outside it give the assertion.

With \(\tau_t^H(A)=e^{itH}Ae^{-itH}\), put
\[
 D_X(s)=\Delta^{-1}\int F(t)\tau_{t/\Delta}^{H_s}\Big(\sum_{i\in X}V_i\Big)dt,\qquad
 T'_X=D_XT_X,\quad T_X(0)=1 .
\]
The odd real frequency multiplier makes \(D_X\) anti-Hermitian. For a high state \(b'\) and low state \(a'\), its matrix element is \(V_{b'a'}/(E_{a'}-E_{b'})\). This is the high-low block of the horizontal generator, by differentiating the spectral projection or its contour. Its low-low block vanishes. Thus \(T_X(s)P=U_X(sq)P\) by the transport equation on these columns, even though the two generators can differ on the excited block. In particular
\[
 E_X=P O_X P,\qquad
 O_X=\int_0^1 T_X(s)^\dagger\Big(\sum_{i\in X}V_i\Big)T_X(s)\,ds.       \tag{7}
\]

For the localization estimate, form a box of radius \(\lfloor n/16\rfloor\) about each interface, including both port sites. Delete the intracell bonds at its two outside edges. For bounded small \(n\), use \(\|T_X(s)-I\|\le bq\|F\|_1\) directly in (7), and \(\|T_i(s)-I\|\le q\|F\|_1\) on a single port. Comparing both expressions with \(P\sum V_iP\) costs \(O(b^2g^2/\Delta)\). The commutator estimate in this case follows from the single-port norm. Enlarging \(c_m\) absorbs \(f_m\) since \(x\le n\). Otherwise the boxes are disjoint. A time-evolved port in the cut Hamiltonian sees only its own box. A direct Lieb–Robinson path estimate compares this observable with the full evolution, uniformly in \(b\). To cover strong boundary bonds, group both endpoints of each port into a single site and go to the interaction picture of the on-site terms on these grouped sites. The remaining interactions are time-dependent nearest-neighbour bonds of norm at most a fixed multiple of \(J_*\), with bounded multiplicity. Support distances between the port and the cut are bounded below by \(c n\) for a numerical \(c>0\). For a norm-one observable on the port and a cut bond \(h\), the commutator with the cut evolution at time \(t'\) is bounded by
\[
 c'\|h\|\min(1,e^{-c n+c'J_*|t'|}).
\]
Indeed iteration of the commutator integral equation gives at most three choices per step on an open bond chain, a factor \(2\|h_Z\|\) per step and ordered time volume \(|t'|^j/j!\). A path must take at least the support distance in steps, up to one endpoint. Weighted summation with \(e^j\) gives the displayed bound. Bounded multiplicities only change \(c'\). On-site interaction-picture unitaries preserve grouped support. Duhamel comparison with the full evolution integrates this bound for the two missing box bonds. Other deleted bonds commute with the cut evolution of this observable. For \(|t'|\le c'' n/J_*\) it gives a relative error at most \(C e^{-c''' n}\), absorbing the integration factor \(n\). For later times use the relative bound two. This proof also works for an operator on an inner box and a cut at a larger radius, with at worst an extra factor \(C n\) before the exponential by summing its starting sites.

Integration with \(F(t)\), using \(t'=t/\Delta\), costs at most \(c_m f_m\) in relative error. Indeed the late tail starts at a fixed multiple of \(n\Delta/J_*\), and \(x\le n,\ x\le n\Delta/J_*\). Denote the isolated box transport generators by \(D_i^{\rm box}\). We obtain
\[
 \big\|D_X-\sum_{i\in X}D_i^{\rm box}\big\|
 \le c_m b(g/\Delta) f_m .
\]
The box transports commute on disjoint supports. By Duhamel for unitaries over \(0\le s\le1\), their product approximates \(T_X\) to the same accuracy. Inserting this in (7) changes \(O_X\) by at most \(2c_m b^2(g^2/\Delta)f_m\), leaving a sum of single-box operators \(O_i^{\rm box}\). The identical boxes are used in the single-port problem. Subtracting those single-port estimates proves assertion 1. In particular the error starts at the scale \(g^2/\Delta\). For a proper connected increment with \(b\ge2\), inclusion and exclusion costs at most an additional factor \(2^b\).

Compression would by itself spoil commutativity of disjoint operators. For a single port take the representative
\[
 \bar O_i=\int l(t)\tau_{t/\Delta}^{H_0}(O_i)\,dt .
\]
It commutes with \(P\) in the sector and has compression \(E_i\). Replace \(O_i\) by \(O_i^{\rm box}\) as above, then cut the last evolution at a box of radius \(\lfloor n/4\rfloor\). The same commutator estimate and the moments of \(l\) approximate \(\bar O_i\) to \(c_m g f_m\) by an operator in the larger box, of norm at most \(c g\). The larger boxes for adjacent interfaces are disjoint. The commutator of the exact representatives is consequently at most \(c_m g^2 f_m\). Restriction to an invariant subspace cannot increase this norm. This proves assertion 2. The cut-box operators only serve for norm comparison and need not preserve the sector.

For assertion 3, use \(y=bz\). The analytic function \(F_X(y)=C_X(y/b)\) is bounded by \(c2^b\) for \(|y|\le1/8\) by (1) and (2). On the real interval \(0\le y\le1/8\), assertion 1 and inclusion and exclusion bound it by \(c_s 2^b y^2 f_s\) for every \(s\), by taking dimensional ports of norm \(\Delta y/b\). To extract coefficient \(k\), use the \(k\)-th forward difference with step \(h=(1+x)^{-m}/(16k)\). It annihilates powers of degree less than \(k\) and gives \(k!h^k\) on \(y^k\). The higher Taylor terms, after division by \(k!h^k\), cost at most \(c_k 2^b h\). For instance bound coefficient \(j\) by \(c2^b8^j\), each difference on \(y^j\) by \(2^k(kh)^j\), and use \(8kh\le1/2\). The difference of the real values costs at most \(c_{k,s}2^b h^{-k} f_s\). Choosing \(s=m(k+1)\) gives the required decay. Restoring units gives the stated coefficient estimate. For the commutator apply the two forward differences in independent coupling variables to the bounded analytic commutator of the two dimensionless single-port energies. On a fixed real square its norm is at most \(c_s f_s\) by the proof of assertion 2 with separate port strengths. The differences kill all powers below \(j\) in the first variable or \(k\) in the second. At the first surviving total degree \(j+k\), take \(s=m(j+k+1)\) and the common step \((1+x)^{-m}/(16(j+k))\). The same geometric remainder proof applies. This establishes arbitrary inverse powers of the cell length at each fixed order. The constants from the moments and differences depend on the order.

The exact finite-cluster estimate applies up to \(bq\le1/8\). At arbitrary \(b\), the coefficient estimates hold at each fixed order.

### An all-order change of basis on the full ring [proved]

There is also a finite-coupling statement for *arbitrarily many* ports. For \(q\le1/8\) there is a unitary \(\mathcal U\) consisting of at most three layers of disjoint two-cell unitaries, such that in dimensional units
\[
 \mathcal U^\dagger H\mathcal U=H_{\rm ref}+\mathcal R,\qquad
 H_{\rm ref}=H_0+\sum_i A_i,\qquad
 [H_{\rm ref},P]=0,\quad P H_{\rm ref}P=\sum_i E_i.
\]
All \(A_i\) act on the adjacent pair, \(\|A_i\|\le C_l g\) with \(C_l=\int|l(t)|dt\). The remainder has terms supported on intervals of a fixed number of cells, with
\[
 \sup_j\sum_{Z\ni j}\|\mathcal R_Z\|
      \le c_m (g^2/\Delta)(1+x)^{-m} .                           \tag{10}
\]
The same unitary and remainder work for every \(m\). The construction uses the gap on a pair of cells. The remainder includes matrix elements between the code and its complement.

Use \(T_i=T_{\{i\}}(1)\) from (7), on its pair. Set \(\Phi_i(O)=T_i^\dagger O T_i\),
\(A_i=\Phi_i(H_0+V_i)-H_0\) and \(B_i=A_i-\Phi_i(V_i)\). The pair transport implies \([A_i,P_iP_{i+1}]=0\) and the asserted compression. Moreover differentiating the conjugated Hamiltonian along the single-port path gives
\[
 A_i=\int_0^1 T_i(s)^\dagger
        \left[\int l(t)\tau_{t/\Delta}^{H_0+sV_i}(V_i)\,dt\right]T_i(s)\,ds .
\]
The middle operator from differentiation is \(V_i+[H_0+sV_i,D_i(s)]\). At energy difference \(E\) its multiplier on \(V_i\) is \(1+(E/\Delta)\widehat F(E/\Delta)=\chi(E/\Delta)\). This identity is valid on the entire pair. If there are additional invariant sectors on the microscopic extension, the norm identity holds there as well, without requiring a gap there.

We have \(\|T_i-I\|\le\|F\|_1 q\), and the single-box approximation to \(T_i\) has error \(c_m q f_m\). The last display also gives a box approximation to \(A_i\) with error \(c_m g f_m\), by cutting the \(l\) evolution starting from the port and using the box transport. The same bounds on norm and box error, up to constants, hold for \(B_i,\Phi_i(V_i),V_i\). Thus for \(j\ne i\) and any such operator \(O_i\),
\[
 \|\Phi_j(O_i)-O_i\|=\|[O_i,T_j-I]\|\le c_m g q f_m
\]
if the pairs meet, and the difference is zero if they do not. Use disjoint box approximants for \(O_i\) and \(T_j-I\), keeping the \(O(q)\) norm of both the second operator and its approximant. At bounded small \(n\) the same estimate uses just the norm bounds and an enlarged constant.

Compose the \(\Phi_i\)'s in the three disjoint layers. A single layer sends \(H_0\) exactly to \(H_0+\sum_{i\ {\rm in\ layer}}B_i\), since \(H_0\) is a sum over cells. Later layers change each resulting \(B_i\) by at most a fixed number of the displayed errors. Telescope around \(B_i\) itself and use norm preservation under later conjugations. Similarly, for an original \(V_i\) compare the action before its own layer with \(V_i\), and the action after that layer with \(\Phi_i(V_i)\). In the own layer other gates are disjoint. An error is initially supported near the two meeting pairs. Each remaining layer enlarges that support by at most one cell at each end. There are a bounded number of errors per original port and a bounded number touching any cell. Adding \(B_i+\Phi_i(V_i)=A_i\) proves (10). The bounded multiplicity convention covers short rings with repeated neighbours. A single self-port can be diagonalized as one cluster.

The reference ground energy also has a uniform interpretation. Let \(e_{\rm iso}(N)\) be the bottom of \(\sum E_i\) on the code. On each cell use the trace-preserving map which keeps the \(P_i\) block, drops its off-diagonal blocks with \(Q_i=1-P_i\), and replaces the \(Q_i\) block by its trace times a fixed code density matrix. Let \(\mathcal E\) be the product of these maps. Its adjoint is unital and positive, also when tensored with the other cells (the reset has Kraus operators). For each pair, \(\mathcal E^*(E_i)\) and \(A_i\) agree on \(P_iP_{i+1}\), preserve that subspace, and have norms at most \(C_l g\). Their difference is bounded below by \(-2C_lg(1-P_iP_{i+1})\). Since \(1-P_iP_{i+1}\le Q_i+Q_{i+1}\),
\[
 H_{\rm ref}\ \ge\ \mathcal E^*(\textstyle\sum_i E_i)
            +(\Delta-4C_l g)\sum_i Q_i
 \ \ge\ e_{\rm iso}(N) I+(\Delta-4C_l g)\sum_i Q_i.             \tag{11}
\]
For \(4C_l g<\Delta\), the ground energy of the reference on the product sector is exactly \(e_{\rm iso}(N)\), and the bottom on \(1-P\) lies at least \(\Delta-4C_l g\) higher. High-energy code states can overlap complementary states in energy. At prescribed total charge (11) still gives the unrestricted lower bound. A charge-compatible reset or a separate comparison would be needed to replace it by the prescribed-charge ground energy. Equation (10) and its norm estimate hold before or after restricting a conserved charge respected by the transports.

## 4. Propagation speed with interactions [proved]

The following local estimate applies on an open chain, a ring of either parity, or the infinite cell chain. Suppose \(H=\sum_i W_i+E\), the \(W_i\) are Hermitian two-cell terms and \(\|[W_i,W_{i+1}]\|\le\gamma\). The terms of \(E\) have range at most \(s\) cells and \(\sup_j\sum_{Z\ni j}\|E_Z\|\le\epsilon\). Then there is a Lieb–Robinson bound with velocity, in cells per unit time,
\[
 v_{\rm LR}\le c_s(\sqrt\gamma+\epsilon).                           \tag{8}
\]
There is no local dimension factor or dependence on \(\|W_i\|\) in this velocity.

Color the bonds with at most three colors, each a disjoint layer \(A_a\). On an interval of time \(t\) take \(U_0(t)=\prod_{a=1}^3 e^{-it A_a}\) in this fixed order. Its generator \(iU_0'U_0^\dagger\) is the sum of the \(A_a\) conjugated by the preceding layers. The difference from \(\sum_a A_a\) is a sum of conjugates of integrals of \(i[A_b,A_a]\) over an interval of length \(|t|\). This follows by telescoping the conjugations and differentiating a single one. Each elementary commutator comes from adjacent bonds, has norm at most \(\gamma\), and there are at most a fixed number per cell. Conjugating by one disjoint layer adds at most one cell on each end of a support. Hence the exact remaining generator in \(U(t)=U_0(t)V(t)\), including \(U_0^\dagger E U_0\), has bounded range and density at most \(c_s(\gamma |t|+\epsilon)\).

Put \(a=\sqrt\gamma+\epsilon>0\) and reset the factorization at intervals \(t_0=1/a\). For any step of length at most \(t_0\), implement the residual evolution first, then the three final disjoint layers in reverse product order. Each bond unitary in a final layer has a Hermitian logarithm of norm at most \(\pi\), so allow time \(t_0\) for that layer with generator density at most \(\pi/t_0\). This exactly implements each step with bounded range, density \(c_s a\), and total simulated time at most \(4(t+t_0)\) for a physical interval of length \(t\ge0\). For interactions of range \(R\) and density \(h\), the commutator path iteration after the first step has total next-bond norm at most \(Rh\) and needs at least \({\rm dist}/R-O(1)\) steps. The exponential series gives a velocity at most \(2e R^2 h\), with a support-dependent prefactor. The extra \(t_0\) only changes that prefactor. This proves (8). For \(a=0\) the layers commute and have bounded propagation range for all times.

Fix a truncation order \(K\). Gather all single-interface contributions through this order into \(W_i\) and put the multi-interface contributions in \(E\). Section 3 and the interval count give \(\gamma\le c_{K,m}g^2 f_{2m}\) and \(\epsilon\le c_{K,m}g f_m\) for \(q\le1\). For the sum of the *exact* isolated-interface energies take \(E=0\) and use assertion 2 at \(q\le1/8\). Thus in physical sites per unit time both cases have
\[
 U_*\le C_m n g(1+x)^{-m}.                                         \tag{9}
\]
For the truncations \(C_m\) also depends on \(K\). This bound includes the one-cell contributions at opposite ends. A separate, arbitrarily prescribed non-scalar Hamiltonian on the code need not obey it. In a local interacting infrared theory (8) limits a mode visible in the propagation of local operators. It does not count such modes.

For the quadratic applications one can use (9) rigorously for every band speed. To see this without a locality assumption about a Jordan–Wigner string, use the product empty Fock vacuum, an eigenstate of the full number-conserving quadratic Hamiltonian. A normalized one-particle packet of finite cell support can be made from it by a norm-one spin operator \(D\) on that interval. The projector \(B_y\) onto nonempty cell \(y\) kills the vacuum at every time in the Heisenberg picture. Thus the weight of the packet at \(y\) is at most \(\|[\tau_t(B_y),D]\|\). Summing the exponential tail bounds the weight outside every cone slightly larger than (8) by a quantity tending to zero. Finite packets are dense, so the same asymptotic conclusion holds for general packets. Conversely, for a smooth packet on a differentiable Bloch eigenbranch, the limiting distribution of cell position divided by time is the distribution of its group velocity under the packet's squared Fourier amplitude. Indeed the characteristic function is the overlap of the evolved Fourier packet with its translate by \(\theta/t\). The phases converge to \(\exp(i\theta e'(k))\) up to the Fourier orientation, by dominated convergence. Taking a packet in a small regular interval rules out a slope outside the cone. One-sided limits give the linear slopes at crossings, including branches with identical energies. Analytic local eigenbranches can be obtained at a touching by diagonalizing the first non-scalar Taylor coefficient on its eigenspace, splitting by its separated eigenprojections and iterating if necessary. Identically degenerate branches have the same velocity.

## 5. An all-order quadratic budget, and a full fixed-order budget [proved]

Assume a fixed on-cell identification with a Fock space in which the Hamiltonian under consideration is number-conserving and quadratic up to a scalar. Consider either of the following Hamiltonians.

- The sum of the exact isolated-interface energies \(E_i\), at \(q\le q_0(d)\) of (6).
- The full linked truncation through any specified order \(K\), at \(q<1\). Include its one-cell terms and every linked multi-interface term. Constants can depend on \(K\).

For either Hamiltonian there are \(a<1\) and constants depending only on \(d\) in the first case, and on \(K,d\) in the second, such that
\[
 N_F(u)\le C(J_*/u)^a,\qquad \sum_{\rm chiral}u_j\le C' J_* .       \tag{12}
\]
Here \(N_F\) counts the linear Fermi crossings of a given physical speed or higher. The Fermi level is arbitrary. The free Casimir coefficient is half the chiral sum. Only the total Hamiltonian needs to be quadratic. The rank estimates can be applied before summing its terms.

Use a supercell of \(s=K+1\) cells in the second case and a single cell (\(s=1\)) in the first. Every support fits within at most two successive supercells. A supercell boundary has only a \(K\)-dependent number of intervals crossing it. The Bloch analysis uses an infinite periodic array. The finite-ring local estimates hold at every length in cells. Let \(T\) be the hopping between two supercells. Use normalized traces and extract its entries by the left operators \(\sqrt2\,a_j^\dagger(-1)^{\hat n}\) and the right operators \(\sqrt2\,a_k\). They are orthonormal on their respective Fock spaces, by the canonical anticommutation relations and the normalized trace. With the tensor-product representation of fermions the extracted block of the quadratic Hamiltonian is \(T/2\). Terms supported solely inside either supercell give zero there. Embedded local terms and their normalized coefficient submatrices have the norm and rank bounds of Section 2, since tensoring identities changes neither the normalized Frobenius norm nor the rank. Thus (5) gives \(\sum_j s_j(T)^p\le C(d)g^p\) for \(p=p(d)<1\) in the isolated-interface problem. At finite \(K\) the same estimate holds for *every* fixed \(0<p\le1\) with \(C(K,d,p)\) by summing the finite word bounds up to \(K\), using \(b\le K,\ q\le1\), and applying the dilation argument. This step allows higher-degree terms that cancel in the total Hamiltonian.

For completeness, the crossing count in the supplied quadratic proof is a flux calculation. In supercell wave number \(k\), write \(h(k)=h_0+T e^{ik}+T^\dagger e^{-ik}\) relative to the Fermi level. For a normalized crossing eigenvector \(w_j\), set \(X_j=(w_j,e^{ik_j}w_j)\) and
\[
 \mathcal J=\begin{pmatrix}0&iT\\-iT^\dagger&0\end{pmatrix}.
\]
Its quadratic form on \(X_j\) is the signed slope and \(\|X_j\|^2=2\). At distinct phases \(X_i^\dagger\mathcal J X_j=0\) by subtracting \(w_i^\dagger h(k_i)w_j\) from \(w_i^\dagger h(k_j)w_j\) and dividing by \(e^{ik_j}-e^{ik_i}\). At a repeated phase diagonalize the compression of \(h'\) to the zero space. Suppose \(b'\) columns have positive slope at least \(t\). Let \(D_+\) be the spectral projection of \(\mathcal J\) for \([t/4,\infty)\). Its rank is \(\#\{s_j(T)\ge t/4\}\). Let \(Q\) project onto \(\ker(D_+X)\) in column space. Then
\[
 t\,{\rm rank}(Q)\le{\rm tr}(QX^\dagger\mathcal J XQ)
       \le (t/4){\rm tr}(X^\dagger X)=tb'/2 .
\]
Since \({\rm rank}(Q)\ge b'-{\rm rank}(D_+)\), this bounds \(b'\) by twice that rank. Apply also to \(-\mathcal J\). The physical speed is \(sn\) times the absolute slope, giving
\[
 N_F(u)\le 4\,\#\{s_j(T)\ge u/(4sn)\}\le C(ng/u)^p.                \tag{13}
\]

Combine (13) with (9). For every \(\alpha>p\),
\[
 \sum_j u_j^\alpha
 \le \frac{C\alpha}{\alpha-p}(ng)^p U_*^{\alpha-p}
 \le C_{\alpha,m}(ng)^\alpha (1+x)^{-m(\alpha-p)}.                \tag{14}
\]
For \(\Delta\le J_*\), \(ng\le n\Delta=xJ_*\). Otherwise \(x=n\) and the spectral width of a cell is at most \(2(n-1)J_*\), so \(ng\le2J_*x^2\). The same width estimate holds for a restricted local cell with the stated internal norm bound. For a one-site cell with an on-site term bounded by \(J_*\), use \(\Delta\le2J_*\) directly. In both cases \(ng\le2J_*(1+x)^2\). Choose an integer \(m\) with \(m(\alpha-p)\ge2\alpha\) to bound (14) by \(C J_*^\alpha\). Take \(\alpha=1\) for the sum and \(\alpha=(1+p)/2\) for the count via \(N_F(u)u^\alpha\le\sum u_j^\alpha\). This proves (12).

A sublinear count suffices for the proposed Casimir bound by integration up to a speed cutoff. The converse fails for the example of \(m\) sectors, each with \(c=1\) and speed \(J_*/m\). The sum is fixed while every uniform power count with exponent below one fails. This observation is only about the two inequalities.

### A finite-coupling quadratic perturbation of the reference [proved]

The comparison (10) also controls the *sum* of speeds of the full compression \(P\mathcal U^\dagger H\mathcal U P\) under two algebraic hypotheses. Suppose the sum of the isolated-interface energies and the full compression are both translation-periodic number-conserving quadratic operators plus scalars in the same cell Fock spaces. Use a bounded-period coloring for the dressing on the infinite array, for example even and odd bonds. Then for \(q\le q_0(d)\) the full compression satisfies \(A_{\rm free}\le C(d)J_*\), including its finite-coupling terms spanning several interfaces. The comparison can allow changes larger than the narrowest bandwidth.

To prove this, write the compressed residual as \(\sum_Z P\mathcal R_Z P\), supported on bounded cell intervals with density \(\epsilon_{\rm cell}\) from (10). Its *sum* is quadratic by subtraction. Group a fixed number \(s\) of cells, a multiple of the dressing period large enough that each interval meets at most two groups. Denote the residual hopping across a cut by \(T_1\). On the sum of the residual terms across this cut take the Fourier component of charge \(+1\) on the left group and \(-1\) on the right. This map averages conjugations by the two independent charge phases with a scalar phase weight, so it contracts the norm. All one-sided terms vanish. Since the full residual is quadratic, the result is exactly the directed hopping \(\sum_{ab}(T_1)_{ab}a^\dagger_{L,a}a_{R,b}\), of norm at most \(c_s\epsilon_{\rm cell}\). Rotate the two sets of orbitals to singular vector bases. In the Slater state occupying \((|L,j\rangle+|R,j\rangle)/\sqrt2\) for every singular pair, the expectation of this hopping is \(\sum_j s_j(T_1)/2\). Consequently \(\|T_1\|_1\le2c_s\epsilon_{\rm cell}\).

Write \(T_0\) for the hopping of the reference across the same cut. The number of singular values above threshold \(t\) for \(T_0+T_1\) is at most the sum of the two counts at \(t/2\). Indeed on the intersection of the complements of the counted right singular spaces the norm of the sum is less than \(t\). Use the crossing calculation of (13) and integrate. The \(T_0\) contribution is bounded by \(C(ng)^p U_*^{1-p}\). Here (8) applied to the complete compression, with its residual treated as bounded-range terms, still gives \(U_*\le C_m ng f_m\). The \(T_1\) contribution integrated over all positive speeds is at most \(C_s n\|T_1\|_1\), by the elementary integral of the singular value count. Equation (10) and \(ng\le2J_*(1+x)^2\) bound this contribution by \(C_s'J_*\). This proves the assertion. It establishes an infrared budget for the compression despite the loss in the energy comparison (17). To identify this compression with the microscopic infrared still requires control of the off-block part of \(\mathcal R\).

## 6. What the interaction matrices would have to satisfy

**[proved] A conditional stiffness estimate.** In a reciprocal multicomponent Gaussian liquid without mixed current-density terms, normalize the positive stiffness matrices \(V_J,V_N\) so the positive speeds are the square roots of the eigenvalues of \(V_J^{1/2}V_N V_J^{1/2}\). Use physical sites as length units, and let \(V={\rm diag}(u_j)\) hold positive bare band speeds, one per boson. Suppose
\[
 0\le V_J\le b V,\qquad 0\le V_N\le b(V+n\epsilon I).              \tag{15}
\]
Here \(b\) is a uniform dimensionless constant and \(\epsilon\) is an energy in cell units. Then
\[
 \sum_j v_j\le b\sum_j\sqrt{u_j(u_j+n\epsilon)}
      \le b\left(\sum_j u_j+\sqrt{n\epsilon}\sum_j\sqrt{u_j}\right). \tag{16}
\]
To check the noncommuting case, first increase \(V_N\) in the congruence using min-max. Interchange the two square-root factors by equality of the nonzero spectra of \(DD^\dagger\) and \(D^\dagger D\), then increase \(V_J\). The two upper comparison matrices commute. The final inequality is scalar. The wave equation for the canonical density and phase fields, \(\partial_t^2\phi=V_J V_N\partial_x^2\phi\) in the corresponding convention, gives the speed formula. For one mode it allows growth of order \(\sqrt{u_j n\epsilon}\).

For bare spectra covered by (13) with \(p<1/2\), even the square-root term is uniformly summable under \(\epsilon\le c g\). Use (14) with \(\alpha=1/2\) *before* eliminating \(x\). Multiplication by \(\sqrt{n\epsilon}\) gives \(C ng(1+x)^{-m(1/2-p)}\), bounded by \(C'J_*\) on taking \(m(1/2-p)\ge2\). The finite-order result provides every such \(p\). For the all-order isolated model a sufficient explicit smaller threshold is
\[
 q\le\min\{q_0(d),[64(d^2)^{3/4}S_{1/2}]^{-2}\}.
\]
Then \(b_{1/2}\le1/2\) and continuity gives a uniform \(p<1/2\) at that threshold. A geometric bare tower also has finite half moment at fixed first speed, by summing a geometric series. Applying this result requires establishing (15) and \(\epsilon\le cg\) for the liquid.

The half-moment condition has a sharp obstruction at the Gaussian level. In physical units take \(u_j=J_*/j^2\), \(V_J=V\), and \(V_N=V+\kappa I\) for any fixed \(\kappa>0\). The bare sum is bounded and the bare count has exponent \(1/2\). The exact positive speeds of this stable Gaussian theory obey \(v_j\ge\sqrt{\kappa J_*}/j\), so the sum through \(m\) components is unbounded as \(m\) grows. This identifies a specific stiffness pattern which would break the budget if it could be maintained microscopically. Realizing the common density correction as a bounded interaction through one physical port is the missing construction.

**[proved] A second conditional estimate allows unequal chiral speeds.** Normalize the chiral commutator matrix to \(\Sigma={\rm diag}(I,-I)\). For a stable Gaussian energy with positive matrix \(G\) in those coordinates, the signed velocities are the eigenvalues of \(\Sigma G\). They coincide with those of \(G^{1/2}\Sigma G^{1/2}\) (use a limit if necessary). Half the sum of absolute speeds is at most \({\rm tr}(G)/2\). Indeed \(-G\le G^{1/2}\Sigma G^{1/2}\le G\), so tracing the positive and negative spectral subspaces separately gives the inequality. A bound on the trace of the diagonal chiral stiffnesses would control the whole velocity sum, including arbitrarily large allowed mixing.

**[proved] Why the cell hypotheses do not give the needed matrix identification.** Even the first-order compression can be any two-site interaction of a fixed smaller local dimension. Take cells of two sites of dimension \(d=s+1\), with a distinguished blank \(|0\rangle\) and \(s\) data states. Encode \(|a\rangle\) as \((|a,0\rangle+|0,a\rangle)/\sqrt2\), and let the one internal bond be \(\Delta(1-P_c)\). On either physical end the compression of a data matrix unit is half the corresponding logical matrix unit. Put any Hermitian two-\(s\)-state operator \(h\) on the boundary data-data subspace and extend it by zero. It has norm \(\|h\|\) and compression \(h/4\). Rescaling the port makes \(q\) arbitrarily small. This realizes arbitrary interacting first-order dynamics at its own energy scale. Stability of any particular infrared point under the higher orders is a separate hypothesis. The interactions can already have the full energy scale of the leading encoded Hamiltonian.

If one-site cells with on-site penalties are included, there is an exact version of the reduction. Take a site \(\mathbb C^s\oplus{\rm span}\{|e\rangle\}\), \(H_c=\Delta|e\rangle\langle e|\), and extend \(\eta h\) on two data sites by zero on the other blocks. In the conserved all-data sector the full microscopic chain is exactly \(\eta\sum h\). Fixing any allowed positive \(q=\eta\|h\|/\Delta\), a universal (P) for this convention with \(J_*=\Delta\) would imply \(A(h)/\|h\|\le C(s+1)/q\) for arbitrary Hermitian bonds \(h\) with the stated infrared. The two-site-cell construction above gives first-order universality under the stricter open-cell convention.

Nor does a norm \(\epsilon\) per code cell automatically control a sum over species. For example the on-cell operator \(\epsilon|\mathrm{vac}\rangle\langle\mathrm{vac}|\) on \(m\) orbitals has norm \(\epsilon\). Its mixed second occupation differences at the empty configuration equal \(\epsilon\) for every pair of distinct orbitals. The off-diagonal difference matrix \(\epsilon({\bf 1}{\bf 1}^T-I)\) has eigenvalue \((m-1)\epsilon\). Consequently, bounding each coefficient this way can lose a factor growing with the number of orbitals. The example leaves the stiffness of a propagating liquid undetermined.

**[physics-level]** Linearizing a weak, nonresonant fermion model and retaining positive forward-scattering matrices gives a Gaussian liquid. Equations (15) or a chiral trace budget would then give the corresponding velocity bound. General cell data specify no bare fermions, lower bound on all relative bandwidths, exclusion of relevant vertices, or renormalized stiffness trace. Counting operator Schmidt coefficients of a vertex does not count bosons produced by interacting renormalization. One cannot apply (13) to those bosons. Gapping a sector removes its contribution in a specified liquid. A bound along a general flow would also have to control the surviving or newly emerging speeds.

## 7. From a reference Hamiltonian to the chain [proved]

Two general speed bounds do apply directly to the microscopic dynamics. In the interaction picture of the independent cells the time-dependent ports have norm \(g\) and occupy pairs of cells, giving a Lieb–Robinson speed at most a universal multiple of \(ng\) in sites per unit time. Grouping the port endpoints as in Section 3 also gives a speed of order \(J_*\) for \(n\ge2\), even for strong disjoint port bonds. If \(n=1\) with an on-site gap, \(g<\Delta\le2J_*\) suffices for the same conclusion. These bound individual locally propagating speeds. They contain no multiplicity estimate.

For the all-order comparison (10), if \(q\le1/8\) and \(4C_lg<\Delta\), (11) and the variational principle give
\[
 |E_0(H_L)-e_{\rm iso}(N)|\le N\epsilon_{\rm cell},\qquad
 \epsilon_{\rm cell}=c_m(g^2/\Delta)(1+x)^{-m}
\]
on the product cell sector. Charging conventions are as in (11). In fact the dressed microscopic ground state satisfies
\[
 \frac1N\sum_i\langle Q_i\rangle_{\rm dressed}
       \le \frac{2\epsilon_{\rm cell}}{\Delta-4C_l g}.
\]
Use (11) as a lower bound in that state, subtract \(N\epsilon_{\rm cell}\) for the residual, and use a code ground state as the variational upper bound. Thus even the density of excited cells is controlled at finite coupling. In the physical scaled doubling defect the energy comparison gives
\[
 L\left|(E_{2L}-2E_L)-(e_{\rm iso}(2N)-2e_{\rm iso}(N))\right|
        \le 4\epsilon_{\rm cell}L^2/n .                          \tag{17}
\]
At fixed parameters the infrared limit sends \(L\) to infinity first. Even arbitrary inverse powers of \(x\) do not make (17) a bound on the Casimir coefficient. The velocity route in Sections 5 and 6 is direct, but its quadratic or stiffness hypotheses concern the actual infrared Hamiltonian to which it is applied.

One cannot simply take a contour around the entire perturbed code band uniformly in \(L\). For an explicit example, use two-site cells with independent logical spins and excitation flags. Give each flipped flag energy \(\Delta\) on the internal bond and take \(V_i=g Z_{\rm right}Z_{\rm left}\) on the two logical spins at each boundary. The \(N\) pairs are disjoint. The code energies include all \(-Ng+2jg,\ 0\le j\le N\), and the one-flag sector has the same list shifted by \(\Delta\). If \(\Delta=2kg\) and \(N\ge k\), the lists intersect. This happens at arbitrarily small fixed \(q\) along choices of \(k\). For a single interface, levels \(+s g\) in the code and \(\Delta-s g\) with a flag already meet on \(0\le s\le1\) when \(g\ge\Delta/2\). The example is exactly diagonal. Its invariant blocks survive the spectral overlap. The bound (3) grows with order \(k\), leaving convergence of the many-cell series open. Equation (10) bypasses this series at the price of a residual. The intermediate coupling window remains open for unbounded \(A\).

There is also a convergence obstruction without a real crossing. Put a flag of energy \(\Delta\) at a cell end and use a port \(gX_{\rm flag}\otimes Z\) diagonal in a conserved neighbour spin. Additional free logical spins give any spectator degeneracy allowed by the cell size. In a fixed neighbour-spin sector the isolated low eigenvalue is
\[
 \frac{\Delta}{2}\big(1-\sqrt{1+4q^2}\big).
\]
It is smooth for real \(q\). Its Taylor series has radius \(1/2\) because of the square-root singularities at \(q=\pm i/2\). Thus even an isolated term can require a different expansion in the upper part of \(g<\Delta\). This example has no velocity tower.

## 8. The frozen nine-letter chain

**[proved] Cell dictionary.** Use \(b=m+1\) for the physical length, where \(m\) is called \(n\) in the supplied construction notes. In the task's notation \(n=b\), \(D_{\rm code}=2^m\). The local labels comprise four oriented bit states and five head states. In the sector with no bare \(LR\) bond, fix the antikink indicators \(B_i=[x_i\notin L][x_{i+1}\notin R]\) at spacing \(b\). The moving head is allowed to change position. The rule \(L_a H\leftrightarrow H R_a\) keeps the indicator on its bond zero and leaves both neighbouring indicators unchanged. The pairs replaced on the three affected bonds never include an \(LR\) in either direction. The other transitions keep the head/bit orientations. Between two fixed indicators a legal word is a run of \(L\)'s, one head, and a run of \(R\)'s. The tensor product over these legal cell spaces is consequently exact.

The conditional rail switches on a boundary reduce to on-site end switches because the neighbouring cell endpoint already satisfies the condition \(y\notin L\) or \(y\notin R\). The charge-conversion term lies inside one cell. Put it, the switches, the head motion and all zero-field on-site energies in \(H_c\). Only the touching pair-port remains a coupling. It joins \(|WW\rangle\) to two orthogonal unit vectors with amplitudes of magnitude \(|\lambda|\). Thus \(g_{\rm port}=\sqrt2|\lambda|\). The parameter called \(g\) in the code is the intracell conversion amplitude \(g_{\rm conv}\). The forbidden-pair penalty vanishes. On the ring the field in the nine-letter bond is \(fQ/2\) for the doubled charge in the Python file. Work first at \(f=0\) and add this scalar in the desired fixed-total-charge sector, or use the field to select the density and linearize the grand-canonical dispersion there. At nonzero \(f\) the unrestricted isolated cell would have a split code if this field were included in \(H_c\). The sector or chemical-potential convention matters for exact degeneracy.

The neutral head sees a path of \(b\) sites for every bit string. Its ground energy is \(\nu-2\cos(\pi/(b+1))\). Its edge amplitude is \(\alpha=\sqrt{2/(b+1)}\sin(\pi/(b+1))\), and its first clock gap is
\[
 \delta_{\rm clk}=2\big[\cos(\pi/(b+1))-\cos(2\pi/(b+1))\big].
\]
For the tabulated \(m=9\) matched schedule take \(b=10,\ D_{\rm code}=512,\ g_{\rm conv}=2^{-163},\ \lambda=2^{-359}\). The token diagonal relative to the neutral ground is 4 on the bright rail and 1004 on the dark rail. The on-site cost parameter called Delta in the Python bond is \(\Delta_0=4-2\cos(\pi/11)\), with a further dark-rail surcharge of 1000. With the closing link of magnitude at most one, the hopping row sum on either rail including its end switch is at most two. The conversion adds at most \(|g_{\rm conv}|\) in a row. Gershgorin puts every charged-head state above \(2-|g_{\rm conv}|\) relative to the neutral ground. Thus the class gap is exactly \(\Delta=\delta_{\rm clk}\). Each internal hopping block on a physical bond is a weighted matching of norm one. Assign each on-site energy to one adjacent internal bond, including the end switches and assigning at most two sites to a bond. The bound
\[
 J_*=2(|\nu|+|\Delta_0|+1000)+4+|g_{\rm conv}|+|t_d|
\]
suffices. One may simply use \(J_*=2020\) at the selector midpoint \(|\nu|<2\) in the supplied schedule. This member lies far below (6). The \(g_{\rm conv}\)-independent return at order \(\lambda^2\) is scalar because the bit signs along the head paths are removable by the same parity gauge on both rails, including their closing link. The hopping below is the supplied \(\lambda^2g_{\rm conv}^2\) coefficient. Exact isolated-interface energies at finite \(g_{\rm conv},\lambda\) are allowed to contain nonquadratic terms.

**[proved] Other parameter choices.** At general zero-field parameters the gap is the smaller of \(\delta_{\rm clk}\) and the lowest charged-head energy above the neutral ground, provided the latter is positive. In particular the moderate choice with bright diagonal 2 and dark surcharge \(D_0=100\) in the README has a different gap test. Write \(\eta=2-2\cos(\pi/(b+1))\). At zero conversion the two isolated rails have bottoms \(\eta\) and \(\eta+D_0\). The matrix of end switches has norm one when \(|t_d|\le1\). Applying the two-by-two variational bound to the norms of the bright and dark components, then restoring the conversion of norm at most \(|g_{\rm conv}|\), puts the charged bottom above or equal to \(\eta+(D_0-\sqrt{D_0^2+4})/2-|g_{\rm conv}|\). A trial state on one bright rail gives the upper bound \(\eta\). For \(b=10\), \(\eta<\delta_{\rm clk}\). When the lower bound is positive, it is the charged level that sets the true cell gap. These members are in the same fixed-sector class whenever the port norm is below that gap.

**[physics-level] An interacting bound for the clock family given its mode count.** The microscopic cell interaction-picture estimate in Section 7 has velocity at most \(6e\,b g_{\rm port}\) because the commutator iteration has at most three next bonds, with factor \(2g_{\rm port}\) and cell distance one per step. With a walker hopping amplitude \(t_0\) the clock gap is at most \(3\pi^2 t_0/(b+1)^2\), by \(2(\cos x-\cos 2x)=4\sin(3x/2)\sin(x/2)\). If the frozen liquid has total locally propagating central charge at most \(b-1\), as predicted for these members, then at *any* \(g_{\rm port}<\Delta\le\delta_{\rm clk}\)
\[
 A\le 6e\,b(b-1)g_{\rm port}\le18e\pi^2 t_0 .
\]
This conclusion survives arbitrary renormalization of those velocities. More generally a microscopic budget follows if \(n\Delta\,c_{\rm tot}\le C J_*\) is known separately. The cell hypotheses alone impose neither this gap scaling nor this bound on emergent central charge. In the fixed-pattern sector competing lattice periods need not be compared. To assert that sector is the unrestricted ground state at finite port still requires the selector stability input in the notes.

**[numerical] Band widths and actual velocities at the matched nine-orbital parameters.** The following computation uses only the scalar matrices specified in the supplied notes. The calculated speeds are those of the leading logical Bloch Hamiltonian. In units
\[
 E_u=\lambda^2\alpha^4\beta^2 g_{\rm conv}^2,\qquad
 \beta=1/d_{10},\quad d_0=1,\ d_1=1004,\ d_{j+1}=1004d_j-d_{j-1},
\]
set \(t_d=\beta 2^{-35}, f/E_u=2^{-83}\). Use the fluxes \(\phi_A=\pi/(2m),\phi_B=0,\psi_A=\pi,\psi_B=\pi/2\). The transformation \(H_s=D_s H_b D_s^\dagger\) with \(D_s\) equal to one on bright sites and \(e^{ij\phi_s}\) on dark site \(j\) gives a common scalar graph \(H_b\). Both rails and the left switch have real hopping \(-1\), and the entry \(H_b[X_m,B_m]\) on the right is \(-i t_d\). Here \(B_j\) denotes a bright graph vertex. This equality follows from \(\psi_s-m\phi_s=\pi/2\).

Let \(R=(H_b\otimes I+I\otimes H_b)^{-1}\), \(\Omega=|X_m,B_0\rangle\), \(l_j=\Pi^L_j R\Omega,\ r_j=\Pi^R_j R\Omega\) for the bright projectors. Gauge factors cancel at each conversion because it changes species at a bright site. At \(\Omega\) the left \(A\) species has gauge \(i\) and the left \(B\) species gauge one. Multiplying the prefactor \(i\) in the supplied \(T\) formula gives \(-1\). Its two orderings become conjugates. The on-cell formulas give the difference of conjugates. Thus the dimensionless matrices actually used are
\[
 T_{jk}=-2\,{\rm Re}\langle l_j,Rr_k\rangle/\beta^2,\quad
 \mu_{jk}=-2i\,{\rm Im}\big(\langle l_j,Rl_k\rangle+
                                      \langle r_j,Rr_k\rangle\big)/\beta^2.
\]
The physical eigenvalues are \(E_u\) times the eigenvalues of
\[
 h(k)=\mu+2^{-83}I+(T+T^T)\cos k+i(T-T^T)\sin k .
\]
At a simple zero write \(s=\partial_k e(k)\) for the dimensionless slope. For an unnormalized null vector \(w\), both calculations below give \(s\):
\[
 \frac{w^\dagger h' w}{w^\dagger w},\qquad
 \frac{\partial_k\det h}{\sum_j \det h_{\widehat j,\widehat j}}.
\]
The second equality follows from differentiating \(\det(h(k)-e(k)I)=0\). Thus no differentiation of nearly equal energy samples is required. A congruence used to whiten a Bloch matrix preserves its zeros but changes the normalization for its velocities.

The standard-library script *check_d9.py* realifies the positive pair operator, factors it by banded Cholesky, constructs the two matrices above and brackets the roots on two rational half-angle charts. Each of the 18 sign changes gives a distinct root numerically. Multiplication by \(z^9\), \(z=e^{ik}\), gives a polynomial of degree at most 18, which explains why this is the full count when the brackets hold. The values below are decimal floating-point results checked at two precisions and by a null-vector calculation.

| Quantity | Recalculated value |
| --- | ---: |
| Clock gap | \(0.236478881567\) |
| \(\log_{10} E_u\) | \(-377.99135921\) |
| Fermi crossings | 9 positive, 9 negative |
| Smallest angular root and pair-sum spacings | \(0.01680848,\ 0.0000883304\) |
| Linear fit of \(\log_{10}\) eigenvalues of \(-(T+T^T)/2\) against descending orbital index | \(-2.73613\) |
| \(\log_{10}(\text{smallest}/\text{largest})\) for these eigenvalues | \(-21.93059\) |
| Smallest and largest \(|s|\) at the actual roots | \(1.54637\,10^{-24},\ 1.16662\,10^{-22}\) |
| \(\log_{10}(\min|s|/\max|s|)\) | \(-1.87762\) |
| \(A_{\rm free}/(b E_u)=\tfrac12\sum_{\rm chiral}|s|\) | \(3.91992\,10^{-22}\) |

The width fit reproduces the quoted \(r\simeq10^{-2.74}\). A hypothetical geometric speed ladder of nine bands with precisely that ratio would sum to \(1.001823\,v_1\). The actual crossing slopes in the specified Bloch matrix span less than two decades. In particular no pairing and averaging of those absolute slopes can give a spread of 22 decades. Different eigenvectors at the Fermi points account for the distinction between slopes and eigenvalues of the symmetric hopping. The physical conversion is \(b E_u|s|\), since \(k\) is conjugate to the cell index. Both the hypothetical ladder and the directly calculated free Casimir coefficient are compatible with a velocity budget. The computed value is for the leading quadratic Hamiltonian. The microscopic many-cell stability and infrared mode count at finite port remain inputs, as acknowledged in both construction files.

## 9. Remaining statements

**[proved]** The local word, rank, long-cell and finite-coupling comparison bounds above do not grow with the ring length. The exact comparison includes all perturbative orders in its remainder. It identifies a specific interacting reference Hamiltonian, plus an error small per cell. A bound on its Schmidt spectrum and a speed cutoff are separate pieces of information.

**[assumed]** To deduce (P) through an interacting liquid one still needs a velocity-resolved count, a stiffness comparison such as (15) with a sufficiently light bare tail, or a bound on the chiral stiffness trace for the actual infrared theory. One also needs a uniform way to include leakage and the residual (10) in that infrared bound. A direct ring proof would instead need cancellation on the \(1/L\) energy scale in (17). None of these inputs is supplied by local block diagonalization. The frozen test case has a separate conditional bound from its clock gap and predicted number of modes. No family with unbounded microscopic \(A/J_*\) is exhibited here.
