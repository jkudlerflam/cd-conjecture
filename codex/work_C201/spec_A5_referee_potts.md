# Task: adversarial referee of a claimed proof (periodic Potts chain inside twisted XXZ rings), with an independent second derivation of its key steps

You are an independent referee. Below are (I) definitions and claims, (II) a claimed proof produced by another model, verbatim, and (III) your tasks. Assume the proof may contain an error and try to find it.

## (I) Definitions and claims

Fix an integer Q in {2, 3, 4} and L >= 2. Potts ring: (C^Q)^{⊗L}, sites j in Z_L, basis |s>, s in Z_Q. E_j = Q^{-1/2} sum_{a,b} |a><b| on site j, F_{j,j+1} = sqrt(Q) sum_s |s s><s s| on sites j, j+1. N = 2L, e_{2j-1} = E_j, e_{2j} = F_{j,j+1}, indices mod N. H_P^{(s)}(L) = -s sum_{i=1}^{N} e_i, s = ±1. sqrt(Q) = 2 cos(gamma), a = cos(gamma).
Twisted XXZ ring on N spins 1/2: H_X(N; Delta, phi) = (1/2) sum_{i=1}^{N-1} [sigma^x_i sigma^x_{i+1} + sigma^y_i sigma^y_{i+1} + Delta (sigma^z_i sigma^z_{i+1} - 1)] + (1/2)[2 e^{i phi} sigma^+_N sigma^-_1 + 2 e^{-i phi} sigma^-_N sigma^+_1 + Delta (sigma^z_N sigma^z_1 - 1)], and Spec(N; Delta, phi, m) is its spectrum in the sector S^z = m.
Claim 1. Every eigenvalue of H_P^{(s)}(L) lies in the union of Spec(2L; s a, phi, m) over the finite list: m = 0 with phi in {2 gamma, -2 gamma, pi}; m = j in {1,...,L} with phi = 2 pi k / j, k integer.
Claim 2. E_P^{(+)}(L) := min spec H_P^{(+)}(L) = min Spec(2L; a, 2 gamma, 0).
Numerical status (exact diagonalization by the manager): Claim 1 holds for Q = 2, 3 at L <= 5 and Q = 4 at L <= 4, both signs, to 5e-10. Claim 2 holds for Q = 3 at L <= 7, Q = 4 at L <= 6, Q = 2 at L <= 7, to 1e-13. So the statements are very likely true and the question is whether the proof is complete.

## (II) The claimed proof, verbatim (Sections 1 to 4 of the other model's output)

## 1. The Potts generators

[rigorous] Write \(\delta=\sqrt Q\), \(u=Q^{-1/2}\sum_a|a\rangle\), \(P=|u\rangle\langle u|\), and \(D=\sum_a|aa\rangle\langle aa|\). Then \(E=\delta P\) and \(F=\delta D\). Both are Hermitian, and each squares to \(\delta\) times itself. For a site at the left end of a bond,
\[
(P\otimes1)D(P\otimes1)=Q^{-1}(P\otimes1),
\qquad
D(P\otimes1)D=Q^{-1}D .
\]
The first identity follows from \(P|a\rangle\langle a|P=P/Q\) and summation over \(a\). For the second, the matrix element between \(|aa\rangle\) and \(|bb\rangle\) is \(\langle a|P|b\rangle\delta_{ab}=\delta_{ab}/Q\), and all other matrix elements vanish. Exchanging the sites proves the right-end identities. Multiplication by \(\delta^3\) gives \(e_i e_{i\pm1}e_i=e_i\).

[rigorous] Two site generators on distinct sites commute. All bond generators are diagonal in the color basis and commute, including the two coincident bond operators when \(L=2\). A site generator and a bond generator with nonincident supports commute. These cases exhaust the pairs at cyclic distance at least two. Finally, summing the bond \(h\) in the question counts each site term twice with coefficient \(1/2\), giving \(H_P^{(s)}=-s\sum_i e_i\).

## 2. The annular representations that occur in the Potts chain

[rigorous] The issue with the twist is different in the sectors with and without through lines. With through lines, the relevant operator is a cyclic permutation. Without through lines, it is the adjacency matrix of a graph. The following construction identifies both operators before making any use of the XXZ correspondence.

### 2.1. A representation on closed walks

[rigorous] Let \(G\) be the star with center \(o\) and leaves \(1,\ldots,Q\). Write \(A\) for its adjacency matrix and set
\[
\mu_o=\delta,\qquad \mu_a=1\quad(1\le a\le Q).
\]
Then \(A\mu=\delta\mu\), since \(Q=\delta^2\). For even \(n>0\), let \(V_n\) have orthonormal basis the closed walks of length \(n\) with a marked starting position. Thus its basis labels are cyclic strings \(a_1,\ldots,a_n\), with \(a_i\) adjacent to \(a_{i+1}\), and rotations of a string are retained as distinct labels whenever the strings differ. Let \(V_0=\mathbb C^{Q+1}\), with basis the graph vertices.

[rigorous] The following explicit rule assigns a map \(T_D:V_r\to V_n\) to a diagram \(D\) of nonintersecting strings in an annulus, with \(r\) marked string endpoints on the inner boundary and \(n\) on the outer boundary. Endpoints are paired by strings, and closed strings are allowed. A boundary walk labels the gaps between consecutive endpoints. An unmarked boundary circle carries one vertex label. Color every complementary region by a graph vertex, requiring adjacent colors across each string and agreement with both boundary labels. For a region \(R\), let \(b(R)\) count its incident marked-boundary gap intervals, with multiplicity. An entire unmarked boundary circle contributes zero to this count. The matrix element of \(T_D\) is the sum over compatible colorings of
\[
\prod_R \mu_{\operatorname{col}(R)}^{\,\chi(R)-b(R)/2}.
\tag{1}
\]
Here regions are taken after cutting along the strings, so the two sides of a cut are counted separately. In particular, all matrix entries are nonnegative.

[rigorous] This rule respects both isotopy and composition. Isotopy preserves the regions, their Euler characteristics and their boundary incidences. On gluing two annuli, compatible colorings are in bijection with pairs of colorings agreeing on the common boundary walk. Gluing along one interval subtracts one from the sum of the Euler characteristics and removes two gap incidences. The exponent in (1) is therefore additive. Gluing along an unmarked circle changes neither quantity, since its Euler characteristic is zero. Applying this argument to all gaps gives \(T_{D_2D_1}=T_{D_2}T_{D_1}\). Straight strings give the identity, since each intervening strip has Euler characteristic one and two boundary gaps.

[rigorous] Three evaluations will be used below.

1. Inserting an innermost contractible circle in a region of color \(a\) removes a disk from that region and creates a disk of some color \(b\sim a\). Its total multiplier is \(\sum_{b\sim a}\mu_b/\mu_a=\delta\). Nested contractible circles can therefore be removed from the inside out.
2. A single essential circle from \(V_0\) to \(V_0\) has matrix \(A\). Its two regions are annuli, both exponents in (1) vanish, and the two labels must be adjacent.
3. A cup inserting the backtrack \(a,b,a\) has coefficient \(\sqrt{\mu_b/\mu_a}\). Its disk of color \(b\) contributes exponent \(1/2\), and the added boundary incidence in the surrounding region gives the factor \(\mu_a^{-1/2}\). The reversed cap has the transposed matrix. Consequently the cap followed by a cup at position \(i\) acts by
\[
e_i|\ldots,a,b,c,\ldots\rangle
 =\mathbf1_{a=c}\sum_{b'\sim a}
 \frac{\sqrt{\mu_b\mu_{b'}}}{\mu_a}
 |\ldots,a,b',a,\ldots\rangle .
\tag{2}
\]

[rigorous] In \(V_N\), take the subspace where the odd positions are leaves and the even positions are centers. Identify its \(L\) leaf labels with the Potts colors. At an odd position the matrix in (2) has every entry \(1/\delta\), hence is \(E_j\). At an even position it is zero unless the neighboring leaves coincide, in which case it multiplies by \(\delta\), hence is \(F_{j,j+1}\). The subspace and its complementary parity subspace are invariant under every \(e_i\). This realizes the physical Potts chain without assuming that one-step translation acts within its Hilbert space.

### 2.2. Filtration by through lines

[rigorous] Fix \(N\). For \(r=0,2,\ldots,N\), let \(\mathcal F_r\subset V_N\) be the span of \(T_D V_k\) for all diagrams \(D:k\to N\) and even \(k\le r\). Put \(\mathcal F_{-2}=0\). Composition makes each \(\mathcal F_r\) invariant under the \(e_i\), and \(\mathcal F_N=V_N\) by the identity diagram.

[rigorous] Consider \(r>0\) modulo \(\mathcal F_{r-2}\). A diagram \(r\to N\) pairing two inner endpoints represents zero in this quotient. Indeed, an inner-to-inner arc cuts off a disk on the side away from the outer boundary. After removing closed contractible loops, choose an innermost such arc. Its endpoints are consecutive on the inner boundary, and it can be placed in an inner collar. Cutting off this cap factors the map through \(V_{r-2}\). Thus the surviving diagrams have \(r\) through lines. An essential circle cannot coexist with a through line because it separates the two boundaries.

[rigorous] For \(r>0\), define \(W_{r/2,z}\) by annular diagrams with \(r\) through lines, imposing the following rules. Concatenate at the \(N\)-point boundary, remove each contractible circle with weight \(\delta\), set a diagram that joins two through lines to zero, and give one cyclic shift of the inner attachments the scalar \(z\). Choose the shift direction corresponding to leftward seam passage. These rules specify the winding convention for \(z\). To see why the shift is the only extra datum, fix the outer arcs and lift the annulus to a strip. The ordered sequence of inner through endpoints must connect to the ordered sequence of outer through endpoints by an increasing bijection commuting with period \(r\). Such a bijection is \(k\mapsto k+t\) for one \(t\in\mathbb Z\). With the attachments fixed, the disjoint arcs in the cut-open strip are unique up to isotopy, as follows by successively cutting along an outermost arc. Thus shifting the inner attachments accounts also for full turns.

[rigorous] On \(V_r\), that inner shift acts as the cyclic permutation \(R_r\) of the walk labels. In (1) its intervening regions are strips and have weight one. In particular \(R_r^r=1\). For an eigenvector \(R_r\eta=z\eta\), the assignment
\[
D\longmapsto T_D\eta \pmod{\mathcal F_{r-2}}
\tag{3}
\]
respects every defining rule of \(W_{r/2,z}\). It intertwines all \(e_i\) by composition. Since \(R_r\) has an eigenbasis, the direct sum of the maps (3), over this basis, surjects onto \(\mathcal F_r/\mathcal F_{r-2}\). Every parameter used here satisfies \(z^r=1\).

[rigorous] For \(r=0\), define \(W_{0,z}\) using outer pairings in the annulus, with contractible loop weight \(\delta\) and essential loop weight \(z+z^{-1}\). Pairings distinguish the side on which the hole lies. Every essential circle can be moved into the inner collar after contractible circles are removed, since the remaining arcs have both endpoints on the outer boundary. Hence for \(A\eta=\alpha\eta\) the map \(D\mapsto T_D\eta\) intertwines \(W_{0,z}\) with its image whenever \(z+z^{-1}=\alpha\). The maps for an eigenbasis of \(A\) span \(\mathcal F_0\).

[rigorous] The star adjacency matrix has eigenvalues \(\delta,-\delta,0\). On vectors constant on the leaves its matrix in center-value and leaf-value coordinates is \(\left(\begin{smallmatrix}0&Q\\1&0\end{smallmatrix}\right)\). It kills the \(Q-1\) dimensional space of leaf vectors with zero sum and zero center coordinate. For \(Q\le4\) all these eigenvalues lie in \([-2,2]\), so each admits \(z=e^{i\arccos(\alpha/2)}\). In particular the Perron weight \(\alpha=\delta\) admits \(z=e^{\pm i\gamma}\). This is the required reality argument in the zero-through-line sector. Finite-order translation alone would give no bound on \(\alpha\).

[rigorous] Write \(B=\sum_{i=1}^N e_i\). In a basis adapted to the filtration, the characteristic polynomial of \(B\) is the product of those of its successive quotients. Under an intertwining surjection, the characteristic polynomial on the image divides that on the domain, by taking a basis adapted to the kernel. Therefore every eigenvalue of the physical \(B\) occurs on at least one of the \(W_{r/2,z}\) just constructed. This argument applies to reducible link modules as well.

## 3. The XXZ dictionary and the exceptional parameters

[rigorous given cited theorem M] The external representation theorem used here is the following specialization of Morin-Duchesne and Saint-Aubin, Theorems 4.5 and 4.6, with their Definition 2.6 and equations (8), (9). For even \(N\), even \(0\le r\le N\), and nonzero parameters \(u_*,v\), take contractible weight \(u_*^2+u_*^{-2}\), essential weight \(v^N+v^{-N}\), and weight \(v^d\) for total leftward defect displacement \(d\) measured in lattice spacings. On \(S^z=r/2\) let the spin generator vanish on parallel spins and have antiparallel block, in order \(\uparrow\downarrow,\downarrow\uparrow\),
\[
\bar e_i=
\begin{pmatrix}
u_*^2&v^2\\
v^{-2}&u_*^{-2}
\end{pmatrix}
\tag{4}
\]
with cyclically ordered sites also at the seam. The link representation and this spin representation are intertwined by a square matrix which is invertible whenever
\[
(-u_*^2)^{k+r/2}v^N-
(-u_*^2)^{-(k+r/2)}v^{-N}\ne0
\quad\text{for every }1\le k\le (N-r)/2 .
\tag{5}
\]
Only this sufficient condition for invertibility is needed. [A. Morin-Duchesne and Y. Saint-Aubin, J. Phys. A 46, 285207 (2013), full text](https://arxiv.org/pdf/1203.4996)

[rigorous] The seam and displacement conventions agree with \(z=v^N\). Indeed, let \(p(w)\) be the sum of the defect coordinates of a link state \(w\), each taken in \(1,\ldots,N\). In a transition \(w\to w'\), the total lifted leftward displacement is
\[
p(w)-p(w')+Nt,
\]
where \(t\) counts signed leftward seam passages. Rescaling the basis vector \(w\) by \(v^{-p(w)}\) removes the endpoint factor and leaves \(v^{Nt}=z^t\). For \(r=0\) the same identification follows directly from the essential weight. Thus in the convention of Section 2,
\[
2j=r,\qquad m=j,\qquad z=v^N,\qquad e^{i\phi}=z^2 .
\tag{6}
\]

[rigorous given cited theorem M] Fix \(u_*^2=e^{i\gamma}\), including \(\gamma=0\). For this fixed \(u_*\), condition (5) excludes only finitely many \(v\ne0\). Outside this finite set the two matrices for \(B\) are similar. Their characteristic polynomial coefficients are Laurent polynomials in \(v\), since the link matrices use only finite displacements and loop factors, and the spin matrices are given by (4). The coefficients agree for all other \(v\), hence agree identically. This proves equality of characteristic polynomials at the exceptional twists themselves, at fixed \(\gamma\). In particular, no invertible intertwiner at a root of unity is being assumed.

[rigorous] Take a unit-modulus \(N\)-th root \(v\) of the \(z\) found in Section 2, and put \(a=\cos\gamma\). Conjugate the spin matrices by \(O=v^{\sum_{\ell=1}^N\ell\sigma_\ell^z}\). On an interior bond the coefficient of \(\sigma_i^+\sigma_{i+1}^-\) becomes \(v^2v^{2i-2(i+1)}=1\). On the seam it becomes \(v^2v^{2N-2}=v^{2N}=e^{i\phi}\). The reverse hoppings have conjugate coefficients. The diagonal part of (4) on each bond is
\[
\frac a2(1-\sigma_i^z\sigma_{i+1}^z)
+\frac{u_*^2-u_*^{-2}}4(\sigma_i^z-\sigma_{i+1}^z).
\]
The second term sums to zero. It follows that \(O(\sum_i\bar e_i)O^{-1}=H_X(N;-a,\phi)\), with no additive constant.

[rigorous] Let \(U=\prod_{i\ {\rm odd}}\sigma_i^z\), which implements the odd-site rotation up to a scalar phase. Every bond, including \((N,1)\), contains one odd site. Thus \(U\) reverses the sign of every hopping, preserves every diagonal term, and preserves \(S^z\). It follows that
\[
U[-H_X(N;-a,\phi)]U^{-1}=H_X(N;a,\phi).
\tag{7}
\]
The twist is unchanged.

[rigorous given cited theorem M] Combining (3) through (7) proves Claim 1. More explicitly, a sufficient finite set of sectors and twists is
\[
\begin{array}{c|c|c}
\text{through lines}&m&\text{allowed seam phases}\\ \hline
r=0&0&e^{i\phi}=z^2,\quad z+z^{-1}\in\{-\sqrt Q,0,\sqrt Q\}\\
r=2,4,\ldots,2L&r/2&e^{i\phi}=z^2,\quad z^r=1 .
\end{array}
\tag{8}
\]
Some of these modules can have zero image. For either sign, inclusion immediately gives
\[
E_P^{(s)}(L)\ \ge\
\min_{m,\ \phi\in\mathbb R}\min\operatorname{Spec}(2L;s a,\phi,m).
\tag{9}
\]

## 4. The ferromagnetic ground state

[rigorous] Use \(r=0\), \(\eta=\mu\), and \(z=e^{i\gamma}\), so every loop has weight \(\delta\). Let \(C:W_{0,z}\to\mathcal H_P\) be the diagram map followed by projection to the physical parity subspace. It intertwines \(B\). Every column of \(C\) is nonnegative and nonzero. To verify the latter assertion, take any loopless basis pairing. Its arcs cut the annulus successively into regions, each cut separating off a disk. Color neighboring regions alternately by the center and one fixed leaf. The boundary gaps then alternate between center and leaf. Interchanging these two assignments if necessary gives the physical parity. The inner region has a strictly positive coefficient in \(\mu\), and every factor in (1) is positive.

[rigorous] Denote the physical and link matrices of \(B\) by \(B_P,B_W\). The matrix \(B_P\) is real symmetric and entrywise nonnegative. A unit vector maximizing its Rayleigh quotient can be replaced by its componentwise absolute value, so it has a nonnegative top eigenvector \(x\). In fact \(x>0\). Each single-site color change is a positive matrix element of \(B_P\), and these changes connect all configurations. A zero coordinate of \(x\) would, by the eigenvector equation, force every neighboring coordinate to vanish and then force \(x=0\).

[rigorous] Let \(\lambda=\lambda_{\max}(B_P)\) and \(y=C^{\mathsf T}x\). The column property gives \(y>0\), and intertwining gives \(y^{\mathsf T}B_W=\lambda y^{\mathsf T}\). The link matrix \(B_W\) is nonnegative because both loop weights are positive. For every complex vector \(f\),
\[
\sum_i y_i|(B_W f)_i|
\le \sum_{i,j}y_i(B_W)_{ij}|f_j|
=\lambda\sum_j y_j|f_j|.
\]
Hence every eigenvalue of \(B_W\) has modulus at most \(\lambda\). Also \(\lambda\) itself is an eigenvalue, by the left eigenvector equation. This proves that the largest link eigenvalue is the largest physical eigenvalue, even if \(C\) has a kernel.

[rigorous given cited theorem M] The characteristic polynomial identity and (7) now give the finite-volume equality
\[
\boxed{\displaystyle
E_P^{(+)}(L)=\min\operatorname{Spec}(2L;a,2\gamma,0)
            =\min\operatorname{Spec}(2L;a,-2\gamma,0).}
\tag{10}
\]
The second equality also follows by complex conjugation in the spin basis. Equation (10) supplies the requested explicit upper bound for \(s=+1\).


## (III) Tasks

1. Referee Section 2.1. The rule (1) assigns to an annular diagram a matrix between closed-walk spaces of the star graph, with region weights mu^{chi(R) - b(R)/2}. Check, with full detail, that it is well defined, invariant under isotopy, and multiplicative under gluing of annuli (including gluing along an unmarked circle, the treatment of regions that are annuli, and the normalization that makes straight strings the identity). Verify evaluations 1 to 3 and formula (2), and that on walks with leaves at odd positions formula (2) reproduces E_j and F_{j,j+1}. If the rule needs a correction (for example in the exponent or in the treatment of the inner boundary), give it and check that the conclusions survive. If you know a published source for this construction (for example Jones's planar algebra of a bipartite graph, or the annular Temperley-Lieb modules of Graham and Lehrer and of Jones), state precisely which theorem covers which step.
2. Referee Section 2.2: the filtration by through lines, the claim that diagrams with an inner-to-inner arc factor through fewer through lines, the definition of the modules W_{r/2,z} and the winding convention, the claim that the inner shift acts on V_r as the cyclic rotation of walk labels with unit weight, the surjection (3), the r = 0 case with essential loop weight an eigenvalue of the adjacency matrix (eigenvalues sqrt(Q), -sqrt(Q), 0), and the characteristic-polynomial argument. In particular decide whether every (j, z) in the list actually needs to be allowed and whether any eigenvalue could come from a non-unimodular z.
3. Referee Section 3: the precise statement of the Morin-Duchesne and Saint-Aubin theorem that is used (J. Phys. A 46, 285207 (2013), arXiv:1203.4996), whether the parameters and the seam convention are matched correctly (z = v^N, e^{i phi} = z^2, m = r/2), the Laurent-polynomial continuity argument at the exceptional parameters, the gauge transformation to H_X(N; -a, phi), and the odd-site rotation. Check the dictionary against the smallest case you can do by hand (N = 2 or N = 4, Q = 2), including the sign of Delta.
4. Referee Section 4: nonnegativity and nonvanishing of the columns of the map C, the Perron-Frobenius argument with a possibly non-injective intertwiner, and the conclusion (10).
5. Independent derivation. Give a second, self-contained proof of Claim 2 by a route that does not use the star-graph walk representation, for example (a) the Fortuin-Kasteleyn or loop representation of the Potts transfer operator and the standard map to the six-vertex model with a seam, or (b) a direct construction of an intertwiner from the zero-through-line link module with all loop weights sqrt(Q) into the Potts space (for instance sending a link pattern to a product of delta functions and uniform sums) together with the same Perron-Frobenius step, or (c) the theorem of Aufgebauer and Klumper, J. Stat. Mech. (2010) P05018, stated precisely. Whatever route you choose, show every step or cite precisely.
6. Verdict: correct as stated, correct after repairs (list them), or wrong (give the broken step). State which parts are fully rigorous and which rest on which citation.

Tag steps [rigorous], [rigorous given cited theorem X] with X stated precisely, or [not proved]. Plain declarative prose, no em-dashes, no probability estimates.
