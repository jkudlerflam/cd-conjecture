# Finite construction and conditional infrared argument

All couplings below are fixed for a chosen n before any thermodynamic limit.
Statements about the central charge are conditional as specified in §8.
The finite algebraic derivations have not had independent review.

## 1 Alphabet and complete bond

There are twenty material letters

\[
 L_a,R_a,V,D_{r,s},X_{r,s},E_a,M_v,M_{r,s},
 \quad a=0,1,\quad r=1,2,\quad s=-,+.
\]

The incoming and outgoing arrows take values 0 and 1:

| Letters | Arrow |
| --- | --- |
| L_a, M_v, M_(r,s) | 0→0 |
| R_a | 1→1 |
| V, D_(r,s), X_(r,s) | 0→1 |
| E_a | 1→0 |

A neighboring pair is allowed exactly when its middle arrows match.
Every other pair has a positive penalty K. The material markers are
E_a,M_v,M_(r,s). Their positions are conserved in the frozen model.

The additive integer charge is 2a on L_a,R_a,E_a, zero on V,M_v, and
s=±1 on D_(r,s),X_(r,s),M_(r,s). Let z_-(a)=(-1)^a and z_+(a)=1,
w_1=w, w_2=wv, φ_-=π/(2n), φ_+=0.

For each a, the following vacuum transitions are positive rank-one terms:

\[
 (q|V R_a\rangle-|L_a V\rangle)(q\langle V R_a|-\langle L_a V|),
\]
\[
 (q|V E_a\rangle-|L_a M_v\rangle)(q\langle V E_a|-\langle L_a M_v|).
\]

Add q²|M_v⟩⟨M_v| onsite. The remaining off-diagonal rules, with their
Hermitian conjugates, are:

| Input | Output | Amplitude |
| --- | --- | --- |
| D_(r,s) R_a | L_a D_(r,s) | −w_r z_s(a) |
| X_(r,s) R_a | L_a X_(r,s) | −exp(iφ_s) z_s(a) |
| X_(r,s) E_a | L_a M_(r,s) | −exp(iφ_s) z_s(a) |
| m D_(r,s), any marker m | m X_(r,s) | −1 |
| D_(r,−) R_1 | D_(r,+) R_0 | g |
| D_(r,−) E_1 | D_(r,+) E_0 | g |
| M_v V | M_(1,−) D_(1,+) | −iλ |
| M_v V | M_(1,+) D_(1,−) | λ |
| M_v V | M_(2,−) D_(2,+) | −λ√η |
| M_v V | M_(2,+) D_(2,−) | λ√η |

There is no closing link from D to the terminal marker. Each excited bus
is a real tree after gauging away the displayed phases.

The onsite operator adds Δ to each excited letter, an additional D to
X_(r,s),M_(r,s), μ_sel to every marker, and f to L_1. Split this onsite
operator equally over adjacent bonds. Together with the vacuum terms,
transition table and K times the forbidden-pair projector, this specifies
the entire repeated bond. No transition is discarded at finite K.

All terms conserve the displayed charge. Every off-diagonal move has the
same incoming left arrow and outgoing right arrow before and after the
move, and both its input and output middle arrows match. Therefore every
forbidden-bond projector commutes with the complete chain Hamiltonian.

## 2 Complete legal cell space

For a payload a_0,...,a_(n-1) and a head at j<n, the word is

\[
 L_{a_0}\cdots L_{a_{j-1}}\,T\,
 R_{a_j}\cdots R_{a_{n-2}}\,E_{a_{n-1}},
\]

where T is V or one of the eight D,X letters. At the terminal position,
the word is L_(a_0)...L_(a_(n-1)) M, with M=M_v or one of four excited
terminal markers. Thus the dimension is (9n+5)2^n.

This list is exhaustive. Starting after a marker, the arrow is 0. L's
preserve it; a V,D,X changes it to 1; thereafter only R's can occur until
an E marker restores 0. Alternatively, an all-L word ends at an M marker.
There cannot be a second head in the same material cell. The last bit is
stored on E while the head occupies a physical data site and moves to the
last L when the head reaches the marker.

At λ=f=0, the vacuum block is Q_n⊗I_(2^n), with Q_n from the README.
Every excited scalar tree has degree at most two, off-diagonal magnitudes
at most one, and diagonal Δ or Δ+D. Conversions have norm at most |g|
per cell. Choose Δ=4+E_n, 0<w≤1, D≥0 and |g|≤1/2. Excited states of the
selected cell then cost at least 2−|g| above its vacuum. The vacuum
excitation gap is at least 1−2q. Thus the complete selected-cell code is
exactly 2^n-dimensional and has an auxiliary gap at least 1/2 for q≤1/4.

The same positive excited-state estimate applies to cells of other
lengths, because their vacuum energy is at most q². Their lowest state
is therefore also their vacuum path ground state.

## 3 Uniform length selection

Write Q_k=B_k†B_k, where B_k has q on the diagonal and −1 on the
superdiagonal. The inverse has (B_k^-1)_ij=q^{-(j-i+1)} for j≥i. Hence

\[
 \|B_k^{-1}\|^2\le\|B_k^{-1}\|_F^2
 =q^{-2k-2}\sum_{s=0}^k(s+1)q^{2s}
 \le q^{-2k-2}(1-q^2)^{-2}.
\]

This proves E_k≥(1−q²)²q^{2k+2}. The trial vector with amplitudes q^j
annihilates every edge term and gives E_k≤q^{2k+2}/Σ_(j=0)^k q^{2j}
≤q^{2k+2}. Consequently E_(k+1)≤rE_k with r≤16/225.

Choose μ_sel=−(n+2)E_n. Subtracting e_0=−E_n per material site leaves
δ_k=E_k+(k−n−1)E_n for a cell of length k+1. It vanishes at k=n.
For k=n+1 it equals E_(n+1); for k≥n+2 it exceeds E_n. For k<n put
m=n−k≥1. Then

\[
 \delta_k\ge [r^{-m}-(m+1)]E_n>E_n>E_{n+1}.
\]

Thus n is selected against every length, with a gap E_(n+1) per wrong
cell relative to the common bulk energy. Marker-free legal rings consist
of all-L or all-R backgrounds and have zero energy, above e_0L.

The gap of Q_k above its lowest eigenvalue follows by deleting its first
row and column. The resulting Toeplitz matrix has smallest eigenvalue
at least (1−q)². Interlacing therefore gives a vacuum gap of at least
(1−q)²−E_k≥1−2q.

For illegal sectors the bad links are static cuts: no off-diagonal term
acts across them. Complete cells inside the resulting legal segments
obey the preceding estimate. Boundary fragments have nonnegative vacuum
terms and nonnegative excited contribution before marker chemical
potentials. There are at most two boundary fragments per cut; each has
at most one marker after complete cells are removed. A conservative
lower bound is therefore

\[
 H_{\lambda=f=0}\ge e_0L+(K-2|\mu_{\rm sel}|)N_{\rm bad}.
\]

Taking K>2|μ_sel| excludes illegal sectors on compatible rings. This is
a zero-port statement. Extending the selected-pattern comparison to the
interacting ground state at λ≠0 is retained as an infrared hypothesis.

## 4 Two charged paths and the quadratic coefficient

Gauge the minus-charge carrier by the parity of the bits it has crossed,
and remove the express-rail phases. The conversion then contains the
logical Jordan-Wigner operator c_i. The terminal minus-charge source has
factor −iP, where P is the left cell's payload parity. Thus the two bus
sources in the pair code have coefficients

\[
 S_r=a_r|+, -\rangle+b_rP|-,+\rangle,
 \qquad (a_1,b_1)=(1,-1),\quad(a_2,b_2)=(1,i),
\]

with an overall √η on bus 2 and endpoint factor v_nv_0 on both buses.

The scalar tree for bus r has vertices D_0,...,D_(n−1), X_0,...,X_(n−1),M.
The D diagonal is 4, the X,M diagonal is 4+D, D hops have magnitude w_r,
and all other path hops have magnitude one. Let

\[
 R_r=(H_r\otimes I+I\otimes H_r)^{-1},\quad
 x_r=R_r|M,D_0\rangle,
\]
\[
 U_i=(|D_i\rangle\langle D_i|\otimes I)x_r,\quad
 V_j=(I\otimes|D_j\rangle\langle D_j|)x_r,\quad
 C^{(r)}_{ij}=U_i^T R_r V_j.
\]

The two one-conversion sectors have operator parts

\[
 Y_{--}=a\sum_iU_i f_{L,i}^\dagger+b\sum_jV_jf_{R,j}^\dagger,
\]
\[
 Y_{++}=P\left[-b\sum_iU_i f_{L,i}+a\sum_jV_jf_{R,j}\right].
\]

The λ²g² term is minus the sum of Y†RY over these sectors, times
|v_0v_n|². The CAR convert it to a scalar and quadratic hopping. As the
tree resolvent is real and |a|=|b|, nonscalar onsite terms cancel. The
hopping coefficient is 2b* a C. Consequently

\[
 T_{\lambda^2g^2}=-2|v_0v_n|^2[C^{(1)}+i\eta C^{(2)}].
\]

The physical coefficient includes λ²g². Different buses have no cross
term because their excited states are orthogonal and propagation preserves
the bus. Distinct spatial bonds have no second-order cross term because
an outer cell remains excited. This is a full-code identity, not merely
a one-particle compression. Higher conversion orders are not Gaussian
in general and are not discarded in the microscopic model.

## 5 The all-n Gram limit and geometric readout

For fixed n and w, eliminate the n+1 high vertices M,X_(n−1),...,X_0.
The M-to-D_0 amplitude has leading factor D^{-(n+1)}. Let
R_0=[A(w)⊗I+I⊗A(w)]^-1 on the two low rails. In the leading term both
high-rail traversals have finished before the intermediate low propagation,
so D^{2n+2}C_ij tends to

\[
 \langle00|R_0(P_i\otimes I)R_0(I\otimes P_j)R_0|00\rangle.
\]

Its Laplace representation is the integral of the product of the two
return-filter integrands over 0<s<r<t. The full product defining Γ_ij
integrates over the square 0<s,r<t. Time reversal (s,r)→(t−s,t−r)
interchanges the two triangles without changing the integrand, because
A(w) is real symmetric. Consequently the two orderings are equal and

\[
 2D^{2n+2}C^{(r)}\longrightarrow\Gamma(w_r),
\]

with Γ defined in the README. The high-vertex expansion is a finite
resolvent expansion at fixed n; extra visits to high vertices or low-to-high
self-energies have higher powers of 1/D. Finite-tree numerical checks are
retained separately.

The lowest nonzero Taylor power of f_i(t;w) is
w^{2i}t^{2i+1}/(2i+1)!, so the f_i are linearly independent and Γ>0.
After dividing entry (i,j) by w^{2i+2j}, C(w) extends positively to w=0,
where

\[
 C(0)_{ij}=\frac{(2i+2j+2)!}
 {(2i+1)!(2j+1)!\,8^{2i+2j+3}}.
\]

For fixed n, the vacuum recurrence or nondegenerate perturbation theory
at q=0 gives v_j=q^j[1+O(q²)]. The onsite L_1 readout has coefficients
γ_i=Σ_(j=i+1)^n v_j²=q^{2i+2}[1+O(q²)]. Choose q=(wu)², and set its
physical coefficient f=Wν/q², where

\[
 W=\lambda^2g^2|v_0v_n|^2D^{-2n-2}.
\]

After the D_w congruence its diagonal is
E(q)_ii=γ_i/(q²w^{4i})=u^{4i}[1+O(q²)]. This is a retained perturbation;
its leakage outside the vacuum code is a higher-order matching correction.

## 6 Why the band argument holds for every n

The argument in C39 §§4–5 depends on the leading geometric behavior of
the readout, not on its being exactly geometric. Here are the needed steps.

1. Fix n and rational 0<w<1/2. Put A=C(w)>0 and
   B(v)=D_v C(wv)D_v. Expanding det(B−tA) in minors shows that the
   generalized eigenvalues have distinct valuations
   r_j=κ_jv^{4j}[1+o(1)], κ_j>0, j=0,...,n−1. The coefficient with k
   rows from B has unique least exponent 2k(k−1), supplied by the
   first k rows and columns.
2. Normalize the j-th eigenvector x_j by (x_j)_j=1. At v=0 its limit
   has zero entries i<j and satisfies (Ax_j)_i=0 for i>j. In particular
   x_j has a nonzero j-th entry. Choose η>0 small enough that the
   centers θ_j=−atan(ηr_j) are all in a small interval about zero.
3. The slopes of the two roots under νE(q) are opposite and have magnitude
   s_j=x_j^T E(q)x_j/[2 sqrt(1+η²r_j²) x_j^T A x_j]. At v=0 and u→0,
   s_j is a strictly positive constant times u^{4j}[1+o(1)]. Choose u
   small so successive limiting slopes have ratio below 1/64. This
   fixes q and the length selector. Next choose v small so the actual
   slopes and absolute centers both have successive ratios below 1/32.
4. At ν=0 the roots are θ_j±π/2. Any difference of two unordered pair
   sums either has a nonzero center combination, or uses the same
   flavor multiset. Nonzero center combinations with coefficient
   ℓ1-norm at most four cannot vanish because of the geometric
   separation. For the same flavor multiset, a nonzero slope combination
   remains, again with ℓ1-norm at most four. It cannot vanish either.
5. There are finitely many comparisons at this fixed n. Analytic
   continuation to sufficiently small nonzero ν therefore preserves
   all previously nonzero differences and splits every remaining one.
   The 2n roots remain simple. Thus both δ_2 and δ_4 are strictly positive.

Repeated labels are included in the pair sums, so this excludes pair
transfer, quartic exchange and umklapp of the leading quadratic reference.
It does not assert rational independence of arbitrarily long combinations.
High-order operators require the separate interacting scaling-dimension
argument, as in the parent construction.

The order of choosing parameters is n,w,η; u and q; v; ν; sufficiently
large finite D; sufficiently small finite g; and sufficiently small finite
λ. The readout f and Δ are then the displayed functions of these values.
K may be chosen at any later stage because the legal dynamics is exactly
independent of K. At each stage the preceding positive finite margins
permit strict inequalities. Root stability here controls finite-order
quadratic matching, not an extensive many-body spectrum.

## 7 The two-color mobile lift

Add holes V_0,V_1 with arrows 0→0 and 1→1. For each material letter x
with arrow a→b, add

\[
 -t(|xV_b\rangle\langle V_a x|+\mathrm{h.c.}).
\]

On adjacent material sites use J times the full material bond. On any
pair involving a vacancy put JK on a forbidden arrow match. Add
−μ_occ(N_occ,left+N_occ,right)/2. This gives a 22-state nearest-neighbour
Hermitian bond, with conserved occupied number, marker number and Q.

At positive occupancy a legal material word has exactly one vacancy
coloring for any vacancy placement. Between successive material letters,
all vacancies must have the common outgoing/incoming arrow. Hence colors
supply no independent extensive state count or low-energy mode. The
legal hopping amplitude is exactly −t, with no multiplicity. Local
arrow preservation still holds on the entire physical Hilbert space.

In squeezed coordinates the charge motion is the usual hard-core vacancy
fluid, with a cyclic material-word twist on a ring. The leading charge-sea
compression of J h_mat is α(ρ)J times the squeezed material Hamiltonian,
where α(ρ)=ρ−sin²(πρ)/(π²ρ)>0. As in the existing vacancy construction,
this is a first-order compression, not a uniform finite-J theorem.

Because the legal space is exactly invariant, K has no effect on its
Hamiltonian. It may be chosen after J/t to penalize illegal cuts. A large
JK therefore introduces no t²/(JK) correction on the legal subspace, in
contrast to a penalty that is crossed by the hopping operator. This does
not eliminate the material wrong-length defect problem inside that space.

The finite-temperature, defect and continuum arguments cannot be replaced
by the literal coloring check. In particular the charge sector is gapless
and must be retained while integrating out material auxiliary excitations.

## 8 Exact remaining assumptions for c=n+1

For each fixed n the following statements would complete the candidate:

- The full material ground state at sufficiently small nonzero λ is in
  the selected length-n sector, with positive bulk-subtracted costs for
  wrong lengths and grammar defects after dressing.
- For a density in a suitable incommensurate window, sufficiently small
  fixed J/t gives a mobile matching with corrections small relative to
  every retained payload scale, while preserving positive defect costs.
  The n-dependent microscopic couplings are already fixed when J/t
  is chosen; no uniform estimate in n is required for unboundedness.
- The nonresonant quadratic reference and the controlled symmetric
  corrections flow to n compact payload bosons, jointly with the vacancy
  boson, with positive velocities and the expected conformal entropy
  coefficient. This includes control of high-order vertices and of
  conformal observables, not just a count of Fermi crossings.

These are substantive thermodynamic inputs. The present calculation
establishes a smaller alphabet and an all-n finite-band mechanism to which
those inputs could be applied. It neither proves those inputs nor measures
the microscopic c at any n.
