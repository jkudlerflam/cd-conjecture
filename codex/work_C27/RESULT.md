# C27. Gluing a common-parity payload to the mobile wrapper

## 1. Result and scope

[proved] For C8's projected wrapper with identity payload transport, the finite-ring gluing is fixed by ordinary cyclic permutation of cell tensor factors. With one common Jordan-Wigner grading, every payload flavour has the same fermionic seam sign

\[
\tau=-\sigma=(-1)^{N_{\rm tot}-1},\qquad
\sigma=(-1)^{N_{\rm tot}}.
\tag{1}
\]

[proved] There is no additional graded permutation in the wrapper quotient. In the free reference, the complete spectrum and the common-parity winding lattice are given in Sections 3-5. For a parity-only interacting payload the exact quotient survives, while separate flavour occupations cease to be exact eigenstate labels.

[heuristic, standard bosonization input] In physical coordinates the parity-even operator with payload labels \((A_\alpha,B_\alpha)\) and integer physical harmonic \(D\) has translation character and length dressing

\[
\begin{split}
q_{\rm pay}&=\sum_\alpha(B_\alpha\theta_\alpha-A_\alpha k_\alpha),\\
Q_T&=\rho(q_{\rm pay}+2\pi D)\pmod{2\pi},\\
\lambda_{\rm phys}
&=\left(2\rho\left[D+\frac{q_{\rm pay}}{2\pi}\right],
A_1,\ldots,A_n;0,B_1,\ldots,B_n\right),\\
\rho&=(b+\nu)^{-1}.
\end{split}
\tag{2}
\]

[proved, conditional on the continuum identification in (2)] If the vertex is invariant under physical translation, write its unreduced momentum as \(2\pi t\), with \(t\in\mathbb Z\). Then

\[
q_{\rm pay}-2\pi t\nu\in2\pi\mathbb Z,
\qquad
\lambda_{\rm phys}=(2t,A;0,B).
\tag{3}
\]

[proved, conditional] Equation (3) is C20's cell resonance test and vertex vector, with \(m=t\). In particular, at the free reference its dimension is

\[
\Delta=t^2+\frac14\sum_\alpha(A_\alpha^2+B_\alpha^2).
\tag{4}
\]

[proved, scope] G7 is therefore resolved for the specified projected identity-transport implementation at the requested physics standard. Its finite-dimensional part is proved below. The continuum operator identification is stated as a bosonization input, not as a theorem about the unfinished microscopic chain. No obstruction here kills route A or this projected implementation. G1/G2, second-order wrapper composition, and microscopic infrared control remain separate obligations. A different transport rule with payload-dependent holonomy would require a different calculation.

[numerical] Independent physical-graph matrices of dimensions 144, 448 and 512 agree with the quotient spectrum, including payloads with pairing and flavour conversion. The largest full-spectrum error is \(2.67\times10^{-14}\). Translation-resolved free payload checks reach dimension 1024. These are finite-spectrum checks, not central-charge estimators.

## 2. Inputs, notation and conventions

[assumed] The implementation is the projected mobile wrapper of C8, with cells of physical lengths \(b\) and \(b+1\), identical internal code spaces \(\mathcal F_n\cong\mathbb C^{2^n}\), and marker hops that transport each code state by the identity. A physical ring has

\[
L=bM+r,\quad 0<r<M,\quad \nu=r/M,\quad a=b+\nu=L/M.
\tag{5}
\]

[proved, inherited input] C8's `joint_fluid.md` supplies the reference cover and its free quotient, while `tagged_vacancy.md` and the wrapper construction supply identity transport in the projected code. C19 Section 2 supplies the local grading and the \(\eta\otimes\gamma\) implementation of a common fermion string. C13 Section 6.2 already records the free common-spin-structure formulas. They are used here, with the translation sign proved explicitly and with their relation to C20's physical local-operator test supplied. The source spectral filters, scale hierarchy and wrapper matrix elements are not recomputed.

[proved, convention] Write \(S(z_0,\ldots,z_{M-1})=(z_1,\ldots,z_{M-1},z_0)\). Translation characters below mean coefficient functions satisfying \(\chi(Sz)=e^{2\pi i m/M}\chi(z)\). Physical coefficient functions have \(F(x,\eta,z)=e^{iPx}\phi(\eta,z)\), with \(P=2\pi\ell/L\). Reversing the active-versus-passive translation convention reverses all momentum signs and leaves every resonance condition unchanged.

[proved, convention] Payload orbital coefficients are \(e^{-ipj}\). We call \(\Theta_\alpha\) the center of the dispersion in this convention:

\[
\epsilon_\alpha(p)=-2t_\alpha\cos(p-\Theta_\alpha)+\mu_\alpha.
\tag{6}
\]

[proved] The local bond \(-t_\alpha e^{i\Theta_\alpha}c^\dagger_{j\alpha}c_{j+1,\alpha}+\mathrm{h.c.}\) has (6) with these coefficient functions. C19 uses the opposite orbital-momentum convention in its displayed dispersion. In the bosonization convention of Sections 1 and 6, annihilation fields carry the factor \(e^{+i(\theta_\alpha+s k_\alpha)y}\), so \(\theta_\alpha=-\Theta_\alpha\), with the branch labels reversed as needed. We distinguish these two centers explicitly. Reversing the whole momentum convention does not change any nonresonance condition. The operator character in (2) is specified by its oscillatory factor; it need not have the sign of a state momentum difference defined with coefficient functions.

## 3. Common Jordan-Wigner grading and ordinary translation

### 3.1 The seam sign

[proved] Let \(a_{j\alpha}\) be the local fermion operators within cell \(j\), including the within-cell ordering of the \(n\) modes. Operators on different cell tensor factors commute. Set

\[
P_j=(-1)^{\sum_\alpha a^\dagger_{j\alpha}a_{j\alpha}},\qquad
c_{j\alpha}=\left(\prod_{h<j}P_h\right)a_{j\alpha},\qquad
\Sigma=\prod_jP_j.
\tag{7}
\]

[proved] The operators \(c\) satisfy the canonical anticommutation relations. Within a cell this is the local CAR algebra. For \(j<h\), the string of \(c_h\) contains \(P_j\), which anticommutes with the odd operator in \(c_j\); moving the two operators past one another therefore gives one minus sign. All remaining factors commute. On an interior bond,

\[
c^\dagger_{j\alpha}c_{j+1,\beta}
=a^\dagger_{j\alpha}P_j a_{j+1,\beta}.
\tag{8}
\]

[proved] On the tensor-product seam the same repeated bond contains \(a^\dagger_{M-1,\alpha}P_{M-1}a_{0\beta}\). Move the full parity through the odd creation operator to obtain

\[
a^\dagger_{M-1,\alpha}P_{M-1}a_{0\beta}
=-\Sigma\,c^\dagger_{M-1,\alpha}c_{0\beta}.
\tag{9}
\]

[proved] In detail, \(c^\dagger_{M-1}c_0=a^\dagger_{M-1}\Sigma P_{M-1}a_0=-\Sigma a^\dagger_{M-1}P_{M-1}a_0\), which gives (9) after multiplying by \(-\Sigma\). Consequently in the sector \(\Sigma=\sigma\), periodic tensor-product hopping is fermionic hopping with \(c_{M\alpha}=\tau c_{0\alpha}\), where \(\tau=-\sigma\). The allowed orbital momenta are

\[
p_{\alpha I}=\frac{2\pi}{M}(I+\beta_\sigma),\qquad
\beta_\sigma=\begin{cases}0,&\sigma=-1,\\1/2,&\sigma=+1,\end{cases}
\quad I\in\mathbb Z_M.
\tag{10}
\]

[proved] The shift is common to all flavours because (7) contains \(P_j\), the total cell parity. Independent flavour-dependent seam signs describe a different spin-chain realization. They cannot replace (10) for this task.

### 3.2 Why the ordinary cyclic shift has no extra sign

[proved] Consider an occupation configuration with \(N_0\) occupied modes in the cell that crosses the reference cut and total occupation \(N\). In a fermionic representation, moving that cell's creation operators past all the other cells gives the reordering sign \((-1)^{N_0(N-N_0)}\). The twisted seam gives \(\tau^{N_0}=(-1)^{N_0(N-1)}\). Their product is

\[
(-1)^{N_0(N-N_0)+N_0(N-1)}
=(-1)^{N_0(2N-N_0-1)}=1.
\tag{11}
\]

[proved] The last equality holds because \(N_0(N_0+1)\) is even. Thus twisted fermionic translation acts as ordinary cell permutation on the occupation basis. This is a basis-by-basis identity, so it also holds on arbitrary superpositions of occupations of the same total parity. It does not require total number conservation.

[proved] For a free Slater state, ordinary cell translation therefore has character

\[
e^{2\pi i m_{\rm pay}/M},\qquad
m_{\rm pay}=\sum_{\alpha,I\ {
m occupied}}I+\beta_\sigma N
\pmod M.
\tag{12}
\]

[proved] This is an integer. If \(\beta=1/2\), the allowed sector has even \(N\); if \(\beta=0\), integrality is immediate. Equivalently one may use \(\beta=(N-1)/2\) before reducing modulo one, since \(N(N-1)/2\) is integral.

[proved] Equations (9)-(11) also determine the seam convention for every parity-even local fermion monomial. One can translate a monomial from an interior support to the seam using the just-proved equality of translation operators. Adding a second independent seam sign for pairing terms would violate this equality.

### 3.3 The finite local operator algebra

[proved] An even monomial in the \(c,c^\dagger\) on a finite interval is a tensor-product operator on that interval: all Jordan-Wigner factors outside it occur an even number of times and cancel. Conversely, expand any parity-even tensor-product operator on a finite cell interval in products of local matrix units. Substituting (7) expresses each even product as a finite-interval even CAR polynomial. These two algebras coincide.

[proved] A parity-odd tensor-product operator instead carries a semi-infinite CAR string in this representation. Such operators exist as observables, but the exact parity symmetry forbids them as Hamiltonian perturbations. It would be incorrect to add their twist fields to the parity-preserving perturbation list, or to exclude odd-parity states from the full Hilbert space. The full spectrum includes both \(\sigma\) sectors, with (10) in each.

## 4. The physical quotient and exact finite spectrum

### 4.1 Derivation of the character constraint

[proved] The reference cover has coordinates \((x,\eta,z)\), with \(x\in\mathbb Z_L\), \(\eta_j\in\{0,1\}\), \(\sum\eta_j=r\), and \(z_j\in\mathcal F_n\). Forgetting which marker is the reference imposes

\[
(x,\eta,z)\sim(x+b+\eta_0,S\eta,Sz).
\tag{13}
\]

[proved] The \(M\) choices of reference act freely: a proper partial circuit of the ring advances the marked position by a strictly positive distance smaller than \(L\). The quotient dimension is therefore

\[
\dim\mathcal H_{L,M}=\frac LM\binom Mr\,2^{nM}.
\tag{14}
\]

[proved] At \(P=2\pi\ell/L\), (13) gives

\[
\phi(S\eta,Sz)=e^{-iP(b+\eta_0)}\phi(\eta,z).
\tag{15}
\]

[proved] Define \(Q_\eta=\sum_jj\eta_j\), so that \(Q_{S\eta}=Q_\eta-r+M\eta_0\), and set \(\phi=e^{-iP Q_\eta/M}\psi\). Substitution into (15) gives, without any approximation,

\[
\psi(S\eta,Sz)=e^{-iPL/M}\psi(\eta,z)
=e^{-2\pi i\ell/M}\psi(\eta,z).
\tag{16}
\]

[proved] In a separated length/payload eigenstate, this is

\[
m_{\rm gap}+m_{\rm pay}+\ell=0\pmod M.
\tag{17}
\]

[proved] Neither ordinary cell permutation nor physical translation changes total parity. Thus these translations preserve the common spin-structure sector; they do not permute it with the opposite-parity sector. For an interacting squeezed Hamiltonian only the simultaneous shift need be conserved. The exact statement is restriction to its character \(m_{\rm total}=-\ell\pmod M\), as in C13 Section 6.1. The additive energy formula below requires the separated reference. No exact independent flavour charge is required for (16).

### 4.2 The complete separated reference

[proved] The gap XX chain contains \(r\) hard-core particles. With the gap gauge in (16), its occupied wave numbers and energy are

\[
\begin{split}
q_j&=\frac{2\pi}{M}\left(I_j+\frac{r-1}{2}\right),\\
m_{\rm gap}&=\sum_jI_j+\frac{r(r-1)}2\pmod M,\\
E_{\rm gap}&=-2J\sum_{j=1}^r\cos\left(q_j+\frac{2\pi\ell}{LM}\right).
\end{split}
\tag{18}
\]

[proved, inherited diagonalization] These are C8's gap-chain formulas, with the same coefficient convention. For any ordinary translation-invariant payload Hamiltonian, diagonalize it simultaneously with the ordinary cyclic shift and write its eigenpairs as \((E_\gamma,m_\gamma)\). Then the full separated projected spectrum, with multiplicity, is

\[
E=\kappa M+E_\gamma+E_{\rm gap},\qquad
m_{\rm gap}+m_\gamma+\ell=0\pmod M,
\quad 0\le\ell<L.
\tag{19}
\]

[proved] Completeness follows from Fourier transformation in \(x\), the invertible diagonal gauge, the quotient restriction (16), and complete diagonalization of the two commuting factors. These operations are unitary restrictions or decompositions of the entire quotient Hilbert space. There is no selection of only a convenient family of candidate states.

[proved] For the number-conserving free payload, (10) and (12) make (19) explicit. Choose occupied orbital subsets for every flavour, retain only those with parity \(\sigma\), and use

\[
E_\gamma=\sum_{\alpha,I\ {\rm occupied}}\epsilon_\alpha(p_{\alpha I})
-\frac M2\sum_\alpha\mu_\alpha.
\tag{20}
\]

[proved] Take the union over both \(\sigma\). For each occupation count vector exactly one seam grid is allowed. Summing over those counts gives \(\sum_{N_1,\ldots,N_n}\prod_\alpha\binom M{N_\alpha}=2^{nM}\), so no extra factor of two is introduced by mentioning the two spin structures.

## 5. Winding lattice and parity sectors

### 5.1 Exact finite sea labels

[proved] In a consecutively occupied sea of \(N_\alpha\) orbitals, the center in orbital units is

\[
c_\alpha\in\mathbb Z+\beta_\sigma+\frac{N_\alpha-1}{2}.
\tag{21}
\]

[proved] This follows by averaging \(I_0+\beta,\ldots,I_0+\beta+N_\alpha-1\). Summing the occupied orbital labels gives \(N_\alpha c_\alpha\). For the length sea the identical formula has \(\beta_0=(r-1)/2\), so its center \(j_0\) is integral. The two half-integer offsets cancel modulo one.

[proved] With \(N=\sum N_\alpha\), (21) is equivalently

\[
2c_\alpha=N+N_\alpha\pmod2.
\tag{22}
\]

[proved] Let \(K_{\rm osc}\) be the signed total change in orbital sum from particle-hole partitions, including length and payload. The ring constraint is

\[
rj_0+\sum_\alpha N_\alpha c_\alpha+K_{\rm osc}+\ell=0\pmod M.
\tag{23}
\]

[proved] Its non-\(\ell\) part is integral. Directly,
\(\sum N_\alpha c_\alpha\) differs from an integer by
\([N^2+\sum N_\alpha^2]/2\), which is integral because
\(N^2\equiv N\equiv\sum N_\alpha^2\pmod2\).

[proved] Equations (5), (21), (23), integer occupations \(0\le N_\alpha\le M\), and distinct occupied orbitals are a finite-volume specification. For low-energy towers one uses interior fillings and bounded particle-hole partitions. That restriction is not part of the full spectrum (19).

### 5.2 Relative winding congruences

[proved] Between two reference seas, let \(Q_\alpha=\delta N_\alpha\) and \(j_\alpha=\delta c_\alpha\). Subtracting (22) gives

\[
2j_\alpha=Q_{\rm tot}+Q_\alpha\pmod2,
\qquad Q_{\rm tot}=\sum_\alpha Q_\alpha.
\tag{24}
\]

[proved] For a parity-preserving perturbation \(Q_{\rm tot}\) is even, and (24) reduces to

\[
Q\in\mathbb Z^n,\quad \sum_\alpha Q_\alpha\in2\mathbb Z,
\qquad 2j_\alpha=Q_\alpha\pmod2.
\tag{25}
\]

[proved] If a state in the opposite total-parity sector is compared instead, every flavour grid changes its spin structure. The additional \(Q_{\rm tot}\) in (24) must then be retained. This is the common-grading information absent from a product of independently graded XX chains. It affects ring sectors but introduces no extra restriction on the number \(n\) of flavours.

[heuristic, continuum identification] The even operator lattice below is the operator version of (25). The odd-sector states can involve the familiar fermionic twist description in a continuum representation. They are included through (24), not through new parity-even Hamiltonian vertices.

### 5.3 Explicit current gluing

[proved] Define

\[
D_\alpha=\frac{M\Theta_\alpha}{2\pi}-c_\alpha,\qquad
f_\alpha=N_\alpha/M,\qquad
\chi_0=\sum_\alpha\frac{N_\alpha\Theta_\alpha}{2\pi}.
\tag{26}
\]

[proved] Lift (23) to
\(\ell=Mt-rj_0-\sum N_\alpha c_\alpha-K_{\rm osc}\), with \(t\in\mathbb Z\), and put \(D_0=bj_0+t\). Then

\[
W_0=j_0+\frac{\ell}{L}
=\frac{D_0-\chi_0+\sum_\alpha f_\alpha D_\alpha-K_{\rm osc}/M}{a},
\qquad P=2\pi W_0\pmod{2\pi}.
\tag{27}
\]

[proved] To verify (27), multiply its middle expression by \(L\):
\(Lj_0+\ell=M(bj_0+t)-\sum N_\alpha c_\alpha-K_{\rm osc}\), then substitute (26) and divide by \(L=aM\). The allowed payload winding cosets are

\[
D_\alpha\in\mathbb Z+\frac{M\Theta_\alpha}{2\pi}
-\beta_\sigma-\frac{N_\alpha-1}{2}.
\tag{28}
\]

[proved] For exact finite counting use (21)-(23), with \(\ell\) in its stated range. Treating all unrestricted integer representatives in (27) as distinct finite states would overcount orbital boosts around the Brillouin zone. In the tower limit, (27)-(28) display the current lattice and its affine offsets.

[proved, asymptotic expansion of finite sums] The filled-sea energy for flavour \(\alpha\) is
\(-2t_\alpha\sin(\pi N_\alpha/M)\cos(2\pi D_\alpha/M)/\sin(\pi/M)\), apart from the field term. The analogous length cosine has center displacement \(2\pi W_0/M\). Expanding \(1-\cos z=z^2/2+O(z^4)\) at bounded currents gives

\[
E_{\rm current}=\frac{2\pi}{L}
\left[v_0W_0^2+\sum_\alpha v_\alpha D_\alpha^2\right]+O(L^{-3}),
\quad
v_0=2Ja\sin(\pi\nu),\quad
v_\alpha=2t_\alpha a\sin(\pi f_\alpha).
\tag{29}
\]

[proved] This recovers C13 Section 6.2. It shows why dropping \(\chi_0\), the common seam coset, or the integer congruence changes the current spectrum. The integer oscillator term in (27) is necessary for exact momenta even when it only affects the current energy beyond leading order.

[assumed] A parity-only interacting infrared liquid retains this compactification if the nonconserving perturbations flow away as in C20. The microscopic exact statements are only total parity, the simultaneous translation character, and (16). Separate \(N_\alpha\), free sea labels and additive energies are not exact interacting quantum numbers.

## 6. Physical translation on the local vertex lattice

### 6.1 Which labels are local

[heuristic, standard bosonization input] Adopt C20's payload fields and chiral convention. If \(h_{\alpha R},h_{\alpha L}\) count annihilations minus creations in a monomial, its vertex and oscillatory momentum are

\[
\begin{split}
A_\alpha&=-h_{\alpha R}+h_{\alpha L},\qquad
B_\alpha=h_{\alpha R}+h_{\alpha L},\\
V_{A,B}&=\exp i\sum_\alpha(A_\alpha\phi_\alpha+B_\alpha\vartheta_\alpha),\\
q_{\rm pay}&=\sum_\alpha(B_\alpha\theta_\alpha-A_\alpha k_\alpha).
\end{split}
\tag{30}
\]

[proved, integer algebra] The inverse formulas are \(h_R=(B-A)/2\), \(h_L=(B+A)/2\). Thus the parity-even local lattice is

\[
A_\alpha,B_\alpha\in\mathbb Z,\qquad
A_\alpha=B_\alpha\pmod2,\qquad
\sum_\alpha B_\alpha=0\pmod2.
\tag{31}
\]

[proved] Section 3.3 shows that the common Jordan-Wigner grading introduces no additional locality condition on these even monomials. There is no condition \(\sum B_\alpha=0\) as an integer. Pairing and flavour conversion must be retained. Conversely, every vector in (31) has integer chiral counts and an even finite CAR representative, using point splitting when identical fermions coincide.

[assumed] Exact marker-number conservation is retained, as in C8/C13. At fixed \(L,M\), it also fixes \(r=L-bM\). Perturbations therefore contain no vertex of the length dual field. Breaking marker-number conservation would enlarge the perturbation lattice and is outside this calculation.

### 6.2 Deriving the physical harmonic and its dressing

[heuristic, standard hydrodynamic input] Let \(y\) be the continuous cell coordinate, with integer values at markers. The smooth gap density is \(\eta(y)=\nu-\partial_y\phi_\ell/\pi\). Since the marker spacing is \(b+\eta\), integrate to obtain

\[
x(y)=a y-\frac{\phi_\ell(y)}{\pi},\qquad
y(x)=\rho x+\frac{\rho}{\pi}\phi_\ell(y(x)).
\tag{32}
\]

[heuristic] The second equation is the inverse displacement relation, to the hydrodynamic order at which the first is used. Local Jacobians contribute smooth density factors and descendants. The leading primary and its oscillatory phase follow from (32). The length field in the following formulas is normalized as the squeezed free XX field, with free scaling parameter one. Rescaling it to the physical marker-density field instead gives C8's \(K=\rho^2\).

[heuristic] A local cell harmonic has the form

\[
e^{i(q_{\rm pay}-2\pi m\nu)y}
 e^{2im\phi_\ell(y)}V_{A,B}(y),\qquad m\in\mathbb Z.
\tag{33}
\]

[proved, distributional identity used in the continuum interpolation] Anchoring an operator to the discrete markers gives an integer comb in \(y\). Its Fourier expansion is \(\sum_{j\in\mathbb Z}\delta(y-j)=\sum_{s\in\mathbb Z}e^{2\pi i s y}\). Hence the same anchored operator has reciprocal branches \(s\in\mathbb Z\). Treating \(q\) and \(q+2\pi\) as the same physical momentum before doing this step loses their different displacement dressing.

[heuristic, substitution into (32)] Put
\(q=q_{\rm pay}-2\pi m\nu+2\pi s\). The branch of (33) becomes

\[
e^{i\rho q x}
\exp i\left[\left(2m+\frac{\rho q}{\pi}\right)\phi_\ell
+\sum_\alpha(A_\alpha\phi_\alpha+B_\alpha\vartheta_\alpha)\right].
\tag{34}
\]

[proved, algebra] Because \(b\) is integral, \(D=s+bm\) is an integer and

\[
\begin{split}
\rho q&=\rho(q_{\rm pay}+2\pi D)-2\pi m,\\
2m+\frac{\rho q}{\pi}
&=2\rho\left(D+\frac{q_{\rm pay}}{2\pi}\right).
\end{split}
\tag{35}
\]

[proved, conditional] On physical lattice sites the factor \(e^{-2\pi imx}\) is one. Equations (34)-(35) therefore give (2). The labels \(m\) and \(s\) do not supply two independent physical harmonic families. Their primary vertex and character depend on the single integer \(D\). Different microscopic representatives can change amplitudes or descendants of the same primary.

[heuristic, locality under the projection] A physical operator of fixed finite range intersects finitely many cells because \(b\) is fixed for the chosen Hamiltonian. Its projection is a finite sum of such anchored operators, with local gap-dependent coefficients. Those coefficients have the integer harmonics in (33). This is the usual continuum expansion of the complete even local algebra, not a claim that every allowed harmonic has a nonzero coefficient in a particular bond.

[proved, finite selection-rule consistency] At the discrete level, a homogeneous projected perturbation is a sum of translated finite-cell terms with local gap indicators. It has no dependence on the arbitrary reference-cell label. In the fixed-\(P\) cover it therefore preserves the simultaneous cyclic character (16). This independent finite statement is consistent with the continuum reciprocal-branch accounting above; physical translation is not ordinary one-cell translation on an isolated payload factor.

### 6.3 Resonances and the C20 dimension test

[proved, conditional on (2)] Choose the real representative
\(\widetilde Q=\rho(q_{\rm pay}+2\pi D)\). Physical translation invariance means \(\widetilde Q=2\pi t\), \(t\in\mathbb Z\). Multiplying by \(a\) gives

\[
q_{\rm pay}-2\pi t\nu=2\pi(bt-D)\in2\pi\mathbb Z.
\tag{36}
\]

[proved, converse] If the cell condition holds, write
\(q_{\rm pay}-2\pi t\nu=2\pi v\), \(v\in\mathbb Z\), and choose \(D=bt-v\). Then \(\widetilde Q=2\pi t\). This proves equivalence in both directions, including umklapp rather than only literal zero momentum.

[proved, conditional] At that resonance, the coefficient of \(\phi_\ell\) in (2) is \(2t\). Thus the resonant lattice, in the same normalized fields as C20, is

\[
\Lambda_{\rm res}=\{(2t,A;0,B):t\in\mathbb Z,
\ A_\alpha=B_\alpha\pmod2,
\ \sum B_\alpha\in2\mathbb Z,
\ q_{\rm pay}-2\pi t\nu\in2\pi\mathbb Z\}.
\tag{37}
\]

[heuristic, Gaussian input] At the decoupled reference the scaling metric for these squeezed normalized fields is the identity, giving (4). For weak interactions use the full scaling metric in these same field coordinates:

\[
\Delta=\frac14(2t,A;0,B)^T\mathcal G(2t,A;0,B).
\tag{38}
\]

[proved, conditional] A common rescaling of the spatial coordinate changes velocities but not this scaling dimension. Changing field coordinates requires transforming both \(\lambda\) and \(\mathcal G\). In particular, replacing the physical character while leaving out the length dressing in (2) would give an incorrect dimension test.

[proved, conditional conclusion] C20's finite list near \(\mathcal G=I\), its parity restrictions, and its zero-momentum test can therefore be used for this implementation. There is no new obligatory resonant pairing or mixed vertex from the common seam. Nonresonant physical dimensions need not equal the dimensions assigned to individual cell harmonics. No claim about that stronger equality is needed for the resonance classification.

[heuristic, quantitative limitation] The scale that stops a detuned perturbation is its physical mismatch \(\operatorname{dist}(\widetilde Q,2\pi\mathbb Z)\), in physical lattice units. Near a specified resonance branch, a cell mismatch \(\delta q\) gives physical mismatch \(\rho\delta q\). At fixed finite \(b\) this is a finite conversion factor, but it must be included when applying the RG amplitude bounds. No new microscopic coupling window is inferred from the finite gluing calculation alone.

## 7. Checks and reproducibility

### 7.1 Analytic checks

[proved, conditional] For \(n=0\), (2) gives \(Q_T=2\pi\rho D\) and free dimension \(\rho^2D^2\). This agrees with C8's physical marker current spectrum and \(K=\rho^2\). It would fail if the reciprocal-branch dressing were omitted.

[proved, conditional] For a rational gap density \(\nu=p/q\) in lowest terms, the first pure-length resonance is \(t=q\), \(D=bq+p\), and its free dimension is \(q^2\). Indeed \(D=t(b+\nu)\) is integral first at \(t=q\). For example \(b=2,\nu=1/2\) has \(\rho=2/5\); the first marker harmonic has dimension \(4/25\) but nonzero momentum \(4\pi/5\), while the first translation-invariant pure-length harmonic has \(D=5,t=2\), dimension four. This reproduces the gap-chain umklapp criterion.

[proved, conditional] For opposite-branch pairing of one flavour, \(A=0,B=2\), so \(q_{\rm pay}=2\theta\). At \(\theta=0\), the \(t=0\) pairing is resonant with free dimension one. Displacing the sea removes this resonance unless a mixed relation \(2\theta-2\pi t\nu\in2\pi\mathbb Z\) holds. The common Jordan-Wigner grading does not protect against pairing by itself.

[proved] At \(n=1\), (22) gives integral \(c_1\) in both total-parity sectors, as for an ordinary periodic XX spin chain. At \(n=2\) the two centers are coupled to the other flavour's number parity. For example \((N_1,N_2)=(1,0)\) uses the periodic grid for both flavours, while \((1,1)\) uses the antiperiodic grid for both. These checks distinguish the common grading from two independent XX chains.

[proved] There is no \(L\bmod b\) obstruction in (5). Every integer \(M\) with \(L/(b+1)<M<L/b\) gives an interior length sector, and this interval contains integers for every sufficiently large \(L\). This says nothing about uniqueness of the interacting global ground state or convergence to one fixed vacuum sector on all ring lengths.

### 7.2 Finite matrices

[numerical, estimator] `gluing_checks.py` uses dense Hermitian diagonalization and maximum absolute eigenvalue error, with full multiplicities. It constructs local cell annihilators with within-cell fermion signs, then builds the periodic payload as an ordinary tensor-product matrix. Translation sectors are obtained from orbit Fourier bases, including short cyclic orbits. A separate free prediction uses the common seam grid (10). The code asserts block invariance before comparing spectra.

[numerical] Free payload results are

| M | n | Full dimension | Maximum translation-resolved energy error |
| :-- | :-- | :-- | :-- |
| 3 | 1 | 8 | 8.89e-16 |
| 3 | 2 | 64 | 2.23e-15 |
| 4 | 2 | 256 | 1.16e-14 |
| 5 | 2 | 1024 | 2.31e-14 |

[numerical, implementation] The independent physical graph enumerates actual ordered marker positions on the ring and all payload words. A marker crossing the physical origin rotates its attached payload labels when the marker list is sorted again. The payload bond is applied directly to cyclically adjacent markers. No quotient formula is used to build this matrix. The independent prediction uses the payload translation blocks and (18)-(19).

[numerical] The physical graph results, with \(b=2\), are

| L | M | r | n | Dimension | Free-payload error | Parity-only-payload error |
| :-- | :-- | :-- | :-- | :-- | :-- | :-- |
| 7 | 3 | 1 | 2 | 448 | 1.96e-14 | 2.05e-14 |
| 8 | 3 | 2 | 2 | 512 | 2.40e-14 | 2.67e-14 |
| 9 | 4 | 1 | 1 | 144 | 9.33e-15 | 6.89e-15 |

[numerical, parameters] The marker hopping is \(J=0.67\); payload hoppings are \((0.8,1.1)\), spectral centers \(\Theta=(0.31,-0.57)\), and fields \((0.13,-0.21)\), truncated to one flavour when needed. The parity-only deformation adds a bond pair term of magnitude 0.073 and phase 0.4. For two flavours it also adds on-cell conversion of magnitude 0.091 and phase 0.2, and on-cell pairing of magnitude 0.061, with Hermitian conjugates. These terms break payload number conservation and the separate flavour parities while preserving total parity. They are a check of the quotient for a parity-only matrix, not a selected gapless parameter set.

[numerical] An additional finite enumeration checks integrality of \(\sum N_\alpha c_\alpha\) and (22) for all \(0\le N_\alpha\le M\) and five consecutive center representatives per flavour. It covers 20, 400 and 625 assignments for \((M,n)=(3,1),(3,2),(4,2)\). The identities themselves are proved in Section 5.

[proved, reproducibility] Run from this folder with Python, numpy and scipy:

    python3 -B gluing_checks.py

[proved, file behavior] The script sets single-thread BLAS, writes `checks.json` in this folder, and imports no neighboring implementation. It uses no network. Every matrix calculation completed within minutes. `PROGRESS.md` records the continuation state.

## 8. Remaining obligations and task accounting

[proved, scope] The finite gluing calculation does not require flavour U(1), a particular filter spectrum, or a particular value of the logical bandwidth. It needs the stated ordinary cell tensor product, one common parity, and identity transport. At the separated reference it works for an arbitrary translation-invariant payload Hamiltonian, including the parity-only matrices tested here.

[assumed] To use it for the final route-A chain, the spatial core must actually realize this grading and transport, the second-order bond must be equalized as required by G3, and local elimination of the gapped sectors must leave the same symmetry and locality structure. The existence of that microscopic construction is not supplied by this report.

[heuristic] Given those inputs and the usual Gaussian infrared identification, the gluing does not change C20's parity-even resonance classification. It couples zero modes and physical momenta, while retaining the proposed \(n+1\) oscillator families. A proof of a thermodynamic CFT, a determination of interacting all-ring uniqueness, and a central-charge measurement are not part of this calculation.

[proved, source coordination] C21/C22 progress and C8 run-8 progress were read. Their compiler, bus and readout work was not duplicated. C25/C26 progress appeared during this calculation and was read before finalization. C25 is computing the gated wrapper resolvents, including the finite-vacancy change of the core spectrum; C26 is computing the Sidon and fourth-order corrections. Neither changes the kinematic assumptions here, and neither task was duplicated. The three named Claude contribution progress files remained absent. The assembly, C13 and C20 are used with their stated qualifications rather than promoted to microscopic theorems.

[proved, workspace accounting] Only this folder was modified. No subagents, network access, cluster computations, repository synchronization, or writes to project STATE.md were used. The user's explicit folder restriction prevents updating that project-level file; the substantive memory is kept in PROGRESS.md here.
