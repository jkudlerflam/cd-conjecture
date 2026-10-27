# Gapped cells above a gapless logical chain

[proved] Several requested statements are false. The relevant perturbative parameter is \(t=|\lambda|/\Delta\), not \(W/\Delta\): \(h_1\) can vanish while \(V\) has norm one. Below, “code block” means the all-code projection \(P\). It does not mean a fixed excitation-number sector. No quadratic assumption on \(h_1\) is needed for A–C.

## A. A corrected local Schrieffer–Wolff theorem

[proved] Excitation-number diagonalization can fail already at first order. Give the auxiliary factor three even states of energies \(0,\Delta,2\Delta\), tensor an arbitrary logical Fock space, and include
\[
 V=|1,1\rangle\langle2,0|+|2,0\rangle\langle1,1|.
\]
It has norm one and zero code compression. The displayed states have equal energy and excitation numbers two and one. Between them, \([S_1,H_0]=0\), so no \(S=O(\lambda/\Delta)\) can cancel this first-order matrix element. Here \(H_0=\sum_jH_c(j)\).

[proved] Two other corrections are necessary. Order \(m\) connected terms can occupy \(m+1\) cells, not \(m\). For an explicit second-order example, use even auxiliary qubits with \(H_c=\Delta q\), logical \(A_j=n_j\), and
\[
 V_j=g(A_j\sigma^x_{j+1}+\sigma^x_j A_{j+1}),\qquad 0<g\le1/2.
\]
The code compression vanishes. Creating and removing the excitation at \(j\) gives
\(-g^2\lambda^2(A_{j-1}+A_{j+1})^2/\Delta\), including the three-cell-span term \(-2g^2\lambda^2 A_{j-1}A_{j+1}/\Delta\). Also, the exact BCH series contains connected strings of arbitrarily many overlapping generators. Finite-order support counting therefore does not give the requested finite-range exact remainder.

[proved] For an explicit tail mechanism, take free auxiliary fermions and a hybridization with Bloch amplitude \(q(p)=(1+e^{ip})/2\). The first SW generator rotates the logical/auxiliary two-level matrix by angle \(\theta=\lambda|q|/\Delta\). Its remaining off-diagonal amplitude, up to the phase of \(q\), is
\(\lambda|q|\cos(2\theta)-(\Delta/2)\sin(2\theta)\).
Its Taylor series contains nonzero powers \(|q|^{2r+1}\) for arbitrarily large \(r\), hence arbitrarily long spatial ranges. This illustrates the issue with exact exponentiation, independently of the excitation-number obstruction.

[proved modulo a cited standard lemma] Here is a uniform replacement for every integer \(k\ge1\), using the Cauchy coefficient bound stated below. Write a translation-invariant interaction as a sum of translated interval densities, and set
\[
 \|\Phi\|_\mu=\sum_I e^{\mu(|I|-1)}\|\phi_I\|,
\]
where \(I\) runs over density types with one fixed translation anchor. Multiple densities on the same interval are allowed. Take
\[
 C=256e^2,\qquad C_*=1280eC,\qquad
 t\le[32C(k+1)]^{-1},\qquad \mu=(k+1)^{-1}.
\]
There is an even, translation-invariant, anti-Hermitian generator with density supported on at most \(k+1\) cells and norm at most \(2t\), such that
\[
 e^{S_k}He^{-S_k}=B_k+R_k,\quad
 B_k=H_0+\sum_{m=1}^k\lambda^m\Delta^{1-m}D_m,
\]
\[
 [B_k,P]=0,\qquad
 \|R_k\|_\mu\le C_*^k k!\,|\lambda|^{k+1}/\Delta^k.
\]
The densities of \(D_m\) have range at most \(m+1\). Their code compressions define the effective Hamiltonian. The remainder is exponentially local, not finite range. The constants are deliberately loose. On short rings, use connected cyclic supports, capped at \(M\); the same upper bounds hold.

[proved] Construction and estimates follow. For an interval \(I\), put
\[
 P_I=\prod_{j\in I}GG^\dagger_j,\quad Q_I=1-P_I,\quad
 L_I=Q_I(H_{0,I}|_{Q_I})^{-1}Q_I.
\]
For Hermitian \(A\) supported in \(I\), define
\[
 \mathcal O_I A=P_IAQ_I+Q_IAP_I,
 \quad\mathcal D_I A=A-\mathcal O_I A,
 \quad\mathcal T_I A=L_IAP_I-P_IAL_I.
\]
Block matrices give
\[
 [\mathcal T_I A,H_0]=-\mathcal O_I A,
 \quad\|\mathcal O_I A\|,\|\mathcal D_I A\|\le\|A\|,
 \quad\|\mathcal T_I A\|\le\|A\|/\Delta.
\]
Thus no denominators between excited energies occur.

[proved] Temporarily divide energies by \(\Delta\), write \(V=\sum_jV_j\), and use the signed variable \(z=\lambda/\Delta\). In this construction \(S_i\) denotes the degree-\(i\) coefficient; the full generator is their sum \(\sum_{i=1}^kz^iS_i\). Suppose the coefficients below degree \(m\) are constructed. Let \(A_m\) be the coefficient of \(z^m\) in
\(e^{\operatorname{ad}(\sum_{i<m}z^iS_i)}(H_0/\Delta+zV)\).
Apply the preceding maps density by density, using the dimensionless cell Hamiltonian, and set \(S_m=\mathcal T A_m\), \(D_m=\mathcal D A_m\). This cancels precisely the matrix elements between the local all-code block and its complement. Each \(D_m\) commutes with global \(P\); transitions within its local \(Q_I\) block, including excitation-number-changing transitions, remain in \(B_k\). Global \(P\)-to-\(1-P\) leakage begins in \(R_k\).

[proved] Disjoint even operators commute. Consequently every nonzero nested commutator has connected support, and its range is at most its perturbative degree plus one. Two anchored densities of degrees \(a,b\) have at most \(a+b+1\) overlapping translations, so their commutator density norm is bounded by \(2(a+b+1)\) times the product of norms. All projections and resolvents above are even and translation covariant. This proves the support, parity, and translation assertions inductively.

[proved] Here are explicit coefficient bounds, including control independent of the largest cell energy. Let \(a_1=1\). A majorant for \(\|A_m\|\), hence also \(\|S_m\|\), is
\[
 a_m\le\sum_{q=1}^{m-1}\frac{[2(m+1)]^q}{q!}
 \sum_{i_1+\cdots+i_q=m-1}\prod a_{i_l}
 +\sum_{q=2}^{m}\frac{[2(m+1)]^{q-1}}{q!}
 \sum_{i_1+\cdots+i_q=m}\prod a_{i_l}.
\]
All indices are positive and below \(m\). In the second sum the innermost commutator with \(H_0/\Delta\) is replaced by \(-\mathcal O A_i\). This avoids its potentially large norm. The elementary inequalities
\[
 \sum_{i_1+\cdots+i_q=N}\prod i_l!
 \le4^{q-1}(N-q+1)!,\qquad
 (m+1)^p(m-p)!/m!\le e^{2p}
\]
give \(a_m\le C^{m-1}m!\) by induction. For the first inequality, convolve repeatedly using
\(\sum_{i=1}^{N-1}i!(N-i)!\le4(N-1)!\): the two endpoint terms contribute \(2(N-1)!\), and every interior term is at most \(2(N-2)!\). For the second, the geometric mean of the largest \(p\) integers up to \(m\) is at least \((m!)^{1/m}\ge m/e\). The induction's two contributions, divided by \(C^{m-1}m!\), sum to at most
\[
 (e^u-1)/4+(e^u-1-u)/u<1,\qquad u=8e^2/C=1/32.
\]
Successive generator terms have ratio at most \(Ct(k+1)\le1/32\), proving \(\|s\|\le2t\).

[proved modulo a cited standard lemma] To bound the exact BCH tail, use the complex disk
\( |z|\le\rho=[16C(k+1)]^{-1}\).
There \(\|S_k(z)\|_\mu\le2e|z|\), and
\(\|[S_k(z),H_0/\Delta]\|_\mu\le2e|z|\).
For densities \(A,B\) of range at most \(k+1\), repeated overlap counting gives
\[
 \|\operatorname{ad}_A^q B\|_\mu
 \le(2k\|A\|_\mu)^q (q+2)!\,\|B\|_\mu/2.
\]
Summing with \(1/q!\), using \(\sum_{q\ge0}\binom{q+2}{2}x^q=(1-x)^{-3}\), bounds the transformed Hamiltonian minus \(H_0/\Delta\) by \(20|z|\). The standard Banach-valued Cauchy coefficient bound says that an analytic \(F(z)=\sum F_mz^m\) with \(\sup_{|z|=\rho}\|F(z)\|\le K\) satisfies \(\|F_m\|\le K\rho^{-m}\). Applying it and summing the geometric tail for \(t\le\rho/2\) gives \(40\rho^{-k}t^{k+1}\). Finally \((k+1)^k\le(2e)^k k!\) gives the stated \(C_*\). This also proves convergence in the weighted interaction norm.

## B. What the remainder controls

[proved] Set \(\epsilon=\|R_k\|_0\). Then \(\|R_k\|\le M\epsilon\). The Rayleigh principle gives
\[
 |E_0(H)/M-E_0(B_k)/M|\le\epsilon.
\]
For \(H_u=B_k+uR_k\), differentiation gives
\(\partial_u[-\beta^{-1}\log\operatorname{Tr}e^{-\beta H_u}]=\langle R_k\rangle_u\).
Integration proves the same bound for free energy per cell, uniformly in \(M,\beta>0\). Unitary invariance identifies \(H_1\) with \(H\). All excited blocks of \(B_k\) must be retained in its partition function.

[proved modulo a cited standard lemma] For observables the required comparison is dressed. Put \(\widetilde O=e^{S_k}Oe^{-S_k}\). Duhamel differentiation yields
\[
 \langle O\rangle_H-\langle\widetilde O\rangle_{B_k}
 =-\beta\int_0^1du\int_0^1dv\,
 [\operatorname{Tr}(\rho_u^vR_k\rho_u^{1-v}\widetilde O)
 -\langle R_k\rangle_u\langle\widetilde O\rangle_u].
\]
Use the standard Schatten Hölder inequality \(\|ABC\|_1\le\|A\|_p\|B\|_\infty\|C\|_q\), \(p^{-1}+q^{-1}=1\), with \(A=\rho^v,C=\rho^{1-v}\). Their indicated norms are one. Thus each trace term is bounded by the product of the observable norms, and
\[
 |\langle O\rangle_H-\langle\widetilde O\rangle_{B_k}|
 \le2\|O\|\min\{1,\beta M\epsilon\}.
\]
For even \(O\) on \(L\) consecutive cells,
\(\|\widetilde O-O\|\le\|[S_k,O]\|\le4(L+k)t\|O\|\).
This is an additional dressing error if one compares undressed observables.

[proved] A volume-independent estimate requires a response bound. For example, if the absolute Duhamel covariance of \(O\) and a translated range-\(\ell\) density is at most
\(A\|O\|\epsilon e^{-\operatorname{dist}(\operatorname{supp}O,I)/\xi}\)
uniformly along the interpolation, summing translations gives
\[
 |\delta\langle O\rangle|\le
 \beta A\|O\|\epsilon
 [L+\ell-1+2/(e^{1/\xi}-1)].
\]
For a tail, sum this expression over its densities. This statement is conditional on the specified covariance bound, not a consequence of a per-cell norm alone.

[proved] No bound tending to zero with \(\epsilon\), depending only on \(\beta,L,\epsilon\), holds across arbitrary local Hamiltonians. The classical Ising chain with coupling \(J\) has zero magnetization at zero field, whereas a field \(\epsilon\) gives, in infinite volume,
\[
 m=\sinh(\beta\epsilon)/
 \sqrt{\sinh^2(\beta\epsilon)+e^{-4\beta J}}.
\]
This follows by differentiating the logarithm of the larger transfer-matrix eigenvalue
\(e^{\beta J}[\cosh(\beta\epsilon)+\sqrt{\sinh^2(\beta\epsilon)+e^{-4\beta J}}]\).
Taking \(e^{-2\beta J}\ll\beta\epsilon\to0\) gives \(m\to1\). Additional Hamiltonian-dependent clustering constants are indispensable.

[proved] Likewise, an allowed remainder \(\epsilon\sum_j b_j^\dagger b_j\) shifts each one-particle energy by \(\epsilon\). In a hopping chain adjacent levels near a simple Fermi point differ by \(2\pi v_F/M+O(M^{-2})\). For fixed \(\epsilon>0\), the shift exceeds this spacing for large \(M\), moves the Fermi point, and changes the finite-volume levels crossing zero. Energy-density accuracy is not infrared spectral accuracy.

## C. Exact removal of leakage

[proved] Gaplessness alone neither proves nor disproves exact block diagonalization. If the interaction already preserves auxiliary excitation number, \(U=1\) works even with a gapless logical chain. Conversely, an exact exponentially local \(U=e^S\) would supply the exact conserved quantity \(UNU^\dagger\), with \(N=\sum_jQ_j\), and hence additional conserved sectors. The assumptions supply no such conservation law. The first-order resonance above already excludes a universally valid near-identity construction.

[heuristic] With logical bandwidth \(W\), a cell excitation can distribute energy \(\Delta\) among order \(\Delta/W\) logical excitations. A generic nonzero on-shell matrix element into this continuum gives decay, even when it first appears at very high perturbative order. This motivates approximate, rather than exact, conservation. A golden-rule calculation alone is not a proof of absence of every exponentially local conserved quantity: absence of a quasiparticle pole is insufficient for that conclusion. No generic impossibility theorem is asserted here.

[proved modulo a cited standard lemma] There is a finite-volume low-energy exception. If \(2M|\lambda|<\Delta\), the lowest \(d^M\) eigenvalues, \(d=2^n\), lie in \([-M|\lambda|,M|\lambda|]\), while the others lie above \(\Delta-M|\lambda|\). Use the finite-dimensional min–max principle: the \(r\)-th increasing eigenvalue of \(A\) is \(\min_{\dim L=r}\max_{\psi\in L,\|\psi\|=1}\langle\psi,A\psi\rangle\). Since \(\|\lambda V\|\le M|\lambda|\), every such maximum, and therefore every eigenvalue, changes by at most \(M|\lambda|\). The isolated subspace has the code's dimension, so a unitary maps it to the code and removes leakage from that subspace exactly. This proves neither locality uniform in \(M\) nor conservation of every excitation-number sector. The condition disappears in the thermodynamic limit at fixed coupling.

## D. Imaginary time: obstructions and a sufficient construction

[proved] For \(n\ge1\), the requested graded identification exists if and only if
\[
 D=Ld\text{ for an integer }L,\qquad D_+=D_-=Ld/2.
\]
Indeed any graded tensor product with \(\mathcal F\) has these parity multiplicities. Conversely take an entirely even \(L\)-dimensional auxiliary space, map the graded code to \(|0\rangle\otimes\mathcal F\), and complete orthonormal bases separately in the two parity sectors. This extends to an even unitary. If necessary, embed into a larger space with
\(L=\lceil2\max(D_+,D_-)/d\rceil\).
Finite-energy added states change the thermal trace. An exact embedding of the original thermal problem requires restricting the trace or assigning infinite energy to the added states.

[proved] When factorization exists, let \(\widehat h_1\) act on logical factors in every auxiliary sector and define
\[
 K_j=V_j-\widehat h_{1,j},\qquad
 H_{\rm int}=\lambda\sum_jK_j.
\]
Then \(P_{j,j+1}K_jP_{j,j+1}=0\) and
\(\|K_j\|\le1+\|h_1\|\le2\).
Its block decomposition consists of vacuum-to-excited, excited-to-vacuum, and excited-to-excited terms. This proves the requested algebraic property. However \(H_c\) need not act only on the auxiliary factor, and generally \([H_c,H_{\rm Gauss}]\ne0\).

[proved] The latter obstruction invalidates the proposed Gaussian reference expansion. Duhamel's formula is still exact:
\[
 Z=\sum_{p\ge0}(-1)^p
 \int_{0<t_1<\cdots<t_p<\beta}\!
 \operatorname{Tr}[e^{-(\beta-t_p)A}H_{\rm int}
 e^{-(t_p-t_{p-1})A}\cdots H_{\rm int}e^{-t_1A}]\,d\mathbf t,
\quad A=H_{\rm Gauss}+H_c.
\]
The series converges in every finite volume; insert observables at their specified times for Schwinger functions. Its reference trace is generally not a Gaussian logical trace times an auxiliary trace.

[proved] More decisively, the requested coupling-small kernel theorem is false at finite temperature. Take two logical modes per cell and one auxiliary occupation \(q\), with
\[
 H_c=q(\Delta+u n_1n_2),\qquad u>0.
\]
At \(\lambda=0\), tracing the auxiliary gives logical weights \(A=1+e^{-\beta\Delta}\) on \(00,10,01\), and \(B=1+e^{-\beta(\Delta+u)}\) on \(11\). Consequently
\[
 \langle n_1n_2\rangle-\langle n_1\rangle\langle n_2\rangle
 =A(B-A)/(3A+B)^2\ne0.
\]
All cross-flavor and anomalous two-point functions vanish, so Wick factorization would make this covariance zero. The logical state is therefore non-Gaussian although every positive power of \(\lambda/\Delta\) vanishes. No normalization constant fixes this. A corrected general statement must at least include thermal errors depending on \(\beta\Delta\), or impose auxiliary/logical factorization of \(H_c\), or specify a zero-temperature limit.

[proved] At finite time discretization, a graded factorization does permit an exact Grassmann representation: insert fermionic coherent-state resolutions in the logical factors and trace the remaining finite auxiliary transfer matrices. Extract the chosen Gaussian quadratic form. The residual polynomial has a nonzero constant term, since at zero Grassmann variables it is the trace of a power of the positive auxiliary matrix \(\langle0_F|e^{-\delta\tau H}|0_F\rangle\). Thus its logarithm exists as a finite nilpotent expansion. This proves an algebraic representation, not a volume-uniform kernel estimate.

[heuristic] For general excited-cell dynamics, a promising corrected proof expands auxiliary excursions, keeping logical operators inside them. Excited segments obey 
\(\|Qe^{-sH_c}Q\|\le e^{-\Delta s}\).
Connected components joined by insertions form polymers; winding components describe thermal excitations. The logarithm involves connected clusters of polymers, not only individual polymers. Establishing a summable cluster expansion with uniform bounds on all external-field coefficients is the missing step here. Winding terms cannot be discarded at finite temperature.

[proved] Fermionic signs alone do not perform that step. If a contraction matrix has Gram form \(C_{ij}=\langle f_i,g_j\rangle\), then
\[
 |\det C|\le\prod_i\|f_i\|\prod_j\|g_j\|.
\]
Express the determinant as the inner product of exterior products and apply Cauchy–Schwarz. Gram–Schmidt expresses each exterior-product norm as a product of orthogonal-component lengths, each bounded by the original vector's length, proving the inequality. It avoids factorial growth from summing Wick pairings. It does not bound the auxiliary history sum or coefficients of uncontracted external fields. Those require separate estimates. No general polymer theorem with the requested constants has been proved above.

## E. An exact interacting Grassmann theorem under explicit extra hypotheses

[proved modulo a cited standard lemma] A useful sufficient case allows genuine higher-degree logical interactions. Take \(\mathcal K=\mathcal F(a)\mathbin{\widehat\otimes}\mathcal F\), \(G\psi=|0_a\rangle\otimes\psi\), one auxiliary fermion \(a_j\) per cell, real \(\lambda\), and an even, translation-invariant quadratic \(h_1\). Assume
\[
 H_c(j)=\Delta a_j^\dagger a_j,\qquad
 V_j=h_{1,j}+g(a_j^\dagger B_j+B_j^\dagger a_j),
\]
where \(B_j\) is an odd logical polynomial on cells \(j,j+1\), \(\|B_j\|\le1\), translated identically, and \(g\) is real. Require \(\|h_1\|+2|g|\le1\). Let \(\mathcal B_j\) be its normally ordered Grassmann symbol. Then, for every \(M,\beta>0\),
\[
 Z=Z_{\rm Gauss}(1+e^{-\beta\Delta})^M
 \int P_{\rm Gauss}(d\psi)\,e^{\mathcal V_{\rm eff}(\psi)},
\]
\[
 \mathcal V_{\rm eff}
 =\lambda^2g^2\sum_j\int_0^\beta d\tau d\tau'\,
 \overline{\mathcal B}_j(\tau)C_\Delta(\tau-\tau')\mathcal B_j(\tau'),
 \quad C_\Delta=(\partial_\tau+\Delta)^{-1}_{\rm antiperiodic}.
\]
All Schwinger functions of the tensor-factor logical fermions equal normalized integrals with the corresponding logical insertions and the coherent-state equal-time convention: these insertions contain no auxiliary variables and pass through their integration. The action is even and space-time translation invariant, with the usual fermionic signs across the thermal seam.

[proved modulo a cited standard lemma] The standard lemma used here is the finite-dimensional fermionic coherent-state time-slicing formula: the integral with antiperiodic boundary conditions, kinetic factor \(\sum_l\bar\psi_l(\psi_{l-1}-\psi_l)\), and interaction factor \(-(\beta/N)\sum_l H_{\rm normal}(\bar\psi_l,\psi_{l-1})\) tends to \(\operatorname{Tr}e^{-\beta H}\). Insert coherent-state resolutions and use that the normal-symbol step differs from the exact step by \(O((\beta/N)^2)\) in operator norm, at fixed finite volume. The accumulated error vanishes. At every such approximating discretization the auxiliary integral is Gaussian. The Berezin identity
\[
 \int d\bar a\,da\,e^{-\bar a A a+\bar\eta a+\bar a\eta}
 =\det A\,e^{\bar\eta A^{-1}\eta}
\]
follows by translation of Grassmann variables and expansion of the exponential. Applying it with the odd polynomial sources proportional to \(\mathcal B\) proves the theorem. It bypasses the auxiliary polymer combinatorics entirely.

[proved] Here are explicit kernel estimates. Use ordered-monomial coefficients without factorial denominators. Define their weighted kernel norm as the sum of coefficient total variations, integrating and summing every field coordinate, divided by \(\beta M\), with weight
\(e^{\kappa_t\operatorname{diam}_\beta\{\tau_i\}+\kappa_x\operatorname{diam}\{x_i\}}\).
Equal-time delta functions are measures in this definition. For \(L=2n\) logical modes, expansion in matrix units gives symbol coefficient sum at most
\(b_n=8^{2n}\): there are \(4^L\) matrix units, each expands into at most \(2^L\) normally ordered monomials, and each matrix entry is bounded by the operator norm.

[proved] For \(0<\tau<\beta\),
\[
 C_\Delta(\tau)=e^{-\Delta\tau}/(1+e^{-\beta\Delta}),\qquad
 \int_0^\beta e^{\kappa_t d_\beta(\tau,0)}|C_\Delta(\tau)|d\tau
 \le(\Delta-\kappa_t)^{-1}.
\]
Every effective monomial uses two times and cells within one bond. Hence
\[
 \sum_m\|\mathcal V_{2m}\|_{\kappa_t,\kappa_x}
 \le\frac{\lambda^2g^2 b_n^2e^{\kappa_x}}{\Delta-\kappa_t}.
\]
There are no terms with \(2m>8n-2\). In particular, for \(0<t\le1\), choosing
\(\kappa_t=\Delta/2\), \(\kappa_x=\tfrac12\log(1/t)\) gives the uniform bound
\[
 \sum_m\|\mathcal V_{2m}\|_{\kappa_t,\kappa_x}
 \le2g^2b_n^2\Delta t^{3/2}.
\]
This is also a bound for each field number separately. The unweighted bound is \(\lambda^2g^2b_n^2/\Delta\). No smallness is needed for the identity. A concrete small-coupling formulation takes \(t\le c_0=1/16\), which also ensures \(W/\Delta\le1/16\); at \(t=0\) the logical action is Gaussian. Smallness relative to the logical velocities is a separate condition stated below.

[proved] The second-order vacuum SW operator, in the original general setting, is
\[
 -\lambda^2PVQ(H_0|_Q)^{-1}QVP.
\]
It follows by inserting \(S_1=\mathcal T V\) into the second-order BCH terms: their code compression is \(\tfrac12P[S_1,\mathcal O V]P\). The same expression follows from an ordered virtual excursion because
\(\int_0^\infty Qe^{-sH_0}Q\,ds=QH_0^{-1}Q\).
In the sufficient model, the ordered zero-duration operator limit gives \(-\lambda^2g^2\sum_jB_j^\dagger B_j/\Delta\). For nonlinear \(B\), merely replacing \(C_\Delta\) by \(1/\Delta\) in bare external-field kernels need not reproduce this operator: contact contractions must also be matched. For example, \(B=b_1n_2\) has \(B^\dagger B=n_1n_2\), whereas the coincident Grassmann product \(\overline{\mathcal B}\mathcal B\) vanishes because it repeats the mode-2 variables. The ordered, point-split expression retains the contractions. Thus agreement holds for the matched vacuum effective operator, not a naive multiplication of symbols. The sign in \(\mathcal V_{\rm eff}\) is the sign of minus the Hamiltonian. Finite-temperature winding histories are additional to the vacuum argument.

[proved] Retardation is quantitatively small at low frequency: for a positive excited Hamiltonian \(A\ge\Delta\) and \(|z|<\Delta\), the resolvent identity gives
\[
 \|(A-z)^{-1}-A^{-1}\|
 \le |z|/[\Delta(\Delta-|z|)].
\]
Thus derivative terms carry extra external frequencies divided by \(\Delta\).

[heuristic] Under \(1+1\)-dimensional scaling near simple Fermi points, each extra time derivative lowers relevance by one power of the infrared scale. This explains why the retarded corrections are irrelevant when the instantaneous vertex is marginal. It does not imply that every induced instantaneous vertex is harmless.

[proved] To apply a supplied RG black box, one must additionally verify its \(2n\) non-resonant simple Fermi points and make the even, translation-invariant, exponentially decaying interaction small relative to the logical velocities, not merely relative to \(\Delta\). Any relevant quadratic counterterms must be included in that verification. For example, with minimum velocity \(v_*\), the sufficient model meets a weighted smallness threshold \(\varepsilon_{\rm RG}\) if
\(2g^2b_n^2\Delta t^{3/2}/v_*\le\varepsilon_{\rm RG}\).
Evenness and translation invariance alone do not guarantee a Luttinger liquid: a small translation-invariant pairing term gives spectrum
\(\sqrt{\varepsilon(k)^2+4\delta^2\sin^2k}\), gapping ordinary Fermi points. The black box therefore needs symmetry or explicit exclusions of relevant pairing/mass channels, as well as its non-resonance hypotheses. The general Spec B does not prove those inputs. The theorem in E proves the integration step under its stated additional hypotheses, while the arbitrary-cell, uniformly small polymer theorem remains unproved and its original finite-temperature formulation is false.
