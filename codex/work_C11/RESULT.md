# C11. Infrared fate of route B

The first-order free-Majorana Hamiltonian does not establish infrared stability. A positive coefficient obtained by projecting onto the isotropic Gross–Neveu operator is also insufficient when the actual symmetry is smaller than O(r).

There are three distinct results below.

1. The supplied one-port Hamiltonian admits an explicit second-order Schrieffer–Wolff calculation. The four-Majorana interaction starts at order lambda^2 g^4/Delta^5. Its complete leading flavor tensor is given below, not just its O(r)-invariant projection.
2. For a specified, natural two-port completion, this is also the leading small-g bond tensor. In the family q_+=1, q_-=m, the non-Abelian flavor coupling changes to the marginally relevant sign at m=7. A common chemical potential does not cure that instability. Thus this flavor-symmetric, weak-dressing family does not establish unbounded stable c.
3. A conditional way around this obstruction is to retain exact U(1)^m flavor charges and introduce distinct, small chemical potentials. The dangerous exchange and umklapp processes then oscillate. Weak-coupling field theory gives an open multicomponent Luttinger-liquid region with c=m. This requires additional microscopic symmetry and implementation conditions. It is not a proof that the history-cell construction as supplied already satisfies them.

The joint two-port excited-state Hamiltonian is not specified by the two port compressions in the question. Consequently there is no unique, unconditional finite-g answer for the stated route. I distinguish the given one-port model, a specific two-port completion, and the additional assumptions used for infrared conclusions throughout.

## 1. Conventions and the missing datum

Write tau(X)=tr(X)/q and normalize the Hermitian odd port basis by

\[
 \tau(F_aF_b)=\delta_{ab},\qquad \{\gamma_a,\gamma_b\}=2\delta_{ab}.
\]

This normalization is implicit in the supplied r-dimensional bright state and compression formula. Set Delta>0. Signs of g can be absorbed in the dressing convention. Initially take g>=0 and lambda>0. Reversing lambda reverses the hopping and exchanges the assignment of right and left movers; it does not reverse the second-order interaction.

For a cell with ground isometry G, ground energy e_0, and projector Pi=GG^dagger, the second-order answer requires

\[
 Q(H_{\rm cell}-e_0)^{-1}Q,\qquad Q=1-\Pi,
\]

and both left-left, right-right, and left-right matrix elements of this resolvent. Knowing G^dagger F_L G and G^dagger F_R G does not give these matrix elements. Adding a cell-local term Q X Q leaves the ground states and both first-order compressions unchanged, provided it does not close the cell gap. It changes the second-order answer.

The Hamiltonian explicitly supplied in the question specifies the one-port resolvent. To give a concrete chain calculation, I also examine the following completion, called completion C below.

\[
 H_C=\Delta(Q_L+Q_R)+g(B_L+B_L^\dagger+B_R+B_R^\dagger),
 \tag{1}
\]
\[
 B_L=\sum_a |a\rangle_L\langle0|\otimes1_R\otimes\gamma_a,
 \qquad
 B_R=\sum_a1_L\otimes|a\rangle_R\langle0|\otimes\eta_a.
\]

Here each auxiliary space is a separate port-partner pair, |0>=|Omega>, |a>=|Omega_a>, Q_L and Q_R penalize all states orthogonal to their respective |0>, and eta_a=i gamma_a P. This is an additional modeling assumption. It is not determined by the quoted compression identities.

There is a useful check on this distinction. In completion C the low-energy invariant block connected to |00> is

\[
 \begin{pmatrix}
 0&g\sqrt{2r}&0\\
 g\sqrt{2r}&\Delta&g\sqrt{2r-2}\\
 0&g\sqrt{2r-2}&2\Delta
 \end{pmatrix}\otimes1_{2^{r/2}}.                         \tag{2}
\]

To obtain it, define the normalized singly excited vectors L and R by attaching gamma_a or eta_a to |a0> or |0a>. Define double-excitation isometries

\[
 U={1\over r}\sum_{ab}|ab\rangle\gamma_a\eta_b,
 \quad V={1\over r}\sum_{ab}|ab\rangle\eta_b\gamma_a.
\]

Their overlap is U^dagger V=(r-2)/r. The symmetric combinations (L+R)/sqrt(2) and (U+V)/sqrt(2+2(r-2)/r), together with |00>, give (2). The remaining blocks have positive excitation gaps at sufficiently small g, by continuity from g=0. Thus this block gives an exactly degenerate low-energy multiplet in that regime.

If its real, normalized ground eigenvector is (w_0,w_1,w_2), the port compression is

\[
 \alpha_C=\sqrt{2/r}\,w_0w_1+{\sqrt{2r-2}\over r}w_1w_2
 =-{2g\over\Delta}+{4(r+1)g^3\over\Delta^3}+O(g^5).
 \tag{3}
\]

The stated one-port alpha instead expands as -2g/Delta+4rg^3/Delta^3+O(g^5). The two models agree to leading order but not at finite dressing. This does not disprove the existence of some other two-port realization with the stated alpha. It shows why its excited-state Hamiltonian must be supplied before importing a one-port resolvent.

## 2. Exact one-port spectral tensors

Put

\[
 D=\sqrt{\Delta^2+4rg^2},\quad d={\Delta+D\over2},\quad
 u={\Delta\over D}=\cos2\theta,\quad x=\sin^2\theta={1-u\over2}.
\]

Let V=(1/sqrt(r)) sum_a |a> gamma_a. The ground and upper bright isometries are

\[
 G=\cos\theta\,|0\rangle-\sin\theta\,V,
 \qquad U=\sin\theta\,|0\rangle+\cos\theta\,V.
\]

Measured from the ground energy, the upper bright multiplet has gap D. All states orthogonal to |0> and V have gap d. Denote the corresponding projectors by Pi_D and Pi_d, and the ground projector by Pi_0.

Define the port spectral tensors

\[
 E^s_{ab}=G^\dagger F_a\Pi_sF_bG,\qquad s=0,D,d.
\]

Oddness of F_a implies that every trace of three F's vanishes. Direct multiplication gives

\[
 E^0_{ab}={1-u^2\over r}\gamma_a\gamma_b,
 \qquad
 E^D_{ab}={u^2\over r}\gamma_a\gamma_b,                 \tag{4}
\]
\[
 E^d_{ab}=(1-x)\delta_{ab}1+{x\over r}K_{ab}
             -{1\over r}\gamma_a\gamma_b,
 \quad
 K_{ab}=\sum_{cd}\tau(F_cF_aF_bF_d)\gamma_c\gamma_d.   \tag{5}
\]

For example, G^dagger F_a U=(u/sqrt(r)) gamma_a, which proves the second expression in (4). Equation (5) follows by subtracting E^0+E^D from G^dagger F_aF_bG. These formulas include the dark states. Dropping those states gives the wrong sign and even the wrong small-g limit.

For a single bond between two such endpoint resolvents, with its right Clifford basis replaced by eta when appropriate, the complete second-order operator is

\[
 H^{(2)}_{\rm bond}=-\lambda^2
 \sum_{ab}\sum_{(s,t)\ne(0,0)}
 {E^s_{ab}\otimes E^t_{ab}\over e_s+e_t},
 \quad (e_0,e_D,e_d)=(0,D,d).                         \tag{6}
\]

There is no complex conjugation between the two factors in this expression. The sum over a,b makes the answer Hermitian. Equation (6) is exact to order lambda^2 at arbitrary g/Delta for the specified one-port endpoints.

### Extracting all four-Majorana terms

Use pair labels A=(c,d), c<d, and set Q_A=gamma_c gamma_d. Write

\[
 J_{ab,A}=\delta_{ac}\delta_{bd}-\delta_{ad}\delta_{bc},
\]
\[
 L_{ab,A}=\tau(F_cF_aF_bF_d-F_dF_aF_bF_c).
\]

The coefficients of Q_A in the three E tensors are

\[
 A_0J_{ab,A},\quad B_0J_{ab,A},\quad C_0J_{ab,A}+zL_{ab,A},
\]
\[
 A_0={4x(1-x)\over r},\quad B_0={(1-2x)^2\over r},
 \quad C_0=-{1\over r},\quad z={x\over r}.
\]

Therefore the coefficient matrix in

\[
 H^{(2)}_4=\sum_{A,B}\kappa_{AB}Q_{j,A}Q_{j+1,B}
\]

is

\[
 \kappa=-\lambda^2\{f_{JJ}J^TJ+f_{JL}(J^TL+L^TJ)+f_{LL}L^TL\}, \tag{7}
\]
\[
 f_{JJ}={2A_0B_0\over D}+{2A_0C_0\over d}
       +{B_0^2\over2D}+{2B_0C_0\over D+d}+{C_0^2\over2d},
\]
\[
 f_{JL}={A_0z\over d}+{B_0z\over D+d}+{C_0z\over2d},
 \qquad f_{LL}={z^2\over2d}.
\]

The transpose contracts the port indices a,b without conjugation. Expanding at small g,

\[
 x={rg^2\over\Delta^2}+O(g^4),
\]
\[
 f_{JJ}=-{119x^2\over4\Delta r^2}+O(x^3),\quad
 f_{JL}={7x^2\over4\Delta r^2}+O(x^3),\quad
 f_{LL}={x^2\over2\Delta r^2}+O(x^3).
\]

Consequently

\[
 \boxed{\kappa_{AB}={\lambda^2g^4\over\Delta^5}
 \left[{119\over4}J^TJ-{7\over4}(J^TL+L^TJ)-{1\over2}L^TL\right]_{AB}
 +O(\lambda^2g^6/\Delta^7).}                         \tag{8}
\]

The constants in the remainder can depend on r and q. The small-dressing expansion requires rg^2/Delta^2 small, not just g/Delta small uniformly in r. Independently, SW requires the local bond norm |lambda| ||sum_a F_a tensor F_a|| to be small compared with the cell excitation gap. The continuum conclusions require the induced interactions to be weak compared with the hopping velocity. All parameter statements below are for each fixed finite r; no uniform large-r perturbation bound is assumed.

### Simplification for a complete graded port basis

For q_+=p, q_-=m, q=p+m, take

\[
 F_{i\alpha,x}=\sqrt{q/2}(E_{i\alpha}+E_{\alpha i}),
 \quad
 F_{i\alpha,y}=\sqrt{q/2}(-iE_{i\alpha}+iE_{\alpha i}).
\]

Let C_A=[F_d,F_c], t_A=tr_+(C_A), and S_AB=-tau(C_A C_B). Cyclicity of the trace gives

\[
 L_{ab,A}=\tau(C_AF_aF_b),\qquad J^TL=S.
\]

Completeness of the odd basis gives, for a block-diagonal X,

\[
 \sum_aF_aXF_a=q\,\operatorname{diag}
 (\operatorname{tr}_-X\,1_p,\operatorname{tr}_+X\,1_m).
\]

Since tr(C_A)=0, this proves

\[
 L^TL=-2tt^T,\qquad J^TJ=2\,1.
\]

Thus (8) becomes the particularly useful formula

\[
 \boxed{\kappa={\lambda^2g^4\over\Delta^5}
       \left({119\over2}1-{7\over2}S+tt^T\right)+O(g^6).} \tag{9}
\]

Here t is imaginary, so tt^T is a real, negative semidefinite matrix. Replacing it by tt^dagger would reverse an important contribution.

There are earlier, bilinear terms. The coefficient of i gamma_c gamma_d on one endpoint of one bond is

\[
 h_A=-{\lambda^2g^2\over2i\Delta^3}\sum_a L_{aa,A}+O(g^4).
 \tag{10}
\]

Indeed, only the scalar part of the other endpoint tensor is needed at this order. Since

\[
 \sum_aF_a^2=q\,\operatorname{diag}(m1_p,p1_m),
 \qquad \sum_aL_{aa,A}=(m-p)t_A,
\]

each pair (i alpha,x; i alpha,y) has

\[
 h_A={q(m-p)\over2}{\lambda^2g^2\over\Delta^3}
\]

per incident bond. In the bulk this contribution is doubled. Balanced grading removes this leading charge term; it does not remove the marginal tensor (9).

## 3. Adjacent bonds and completion C

For a fully specified cell, equations (6) and the following adjacent-bond expression are the complete second-order linked-cluster answer. Nonadjacent bonds cannot remove each other's cell excitations and give no second-order linked term.

Let R_cell=Q(H_cell-e_0)^(-1)Q and define

\[
 M_{ab}=G^\dagger(F_{L,a}R_{\rm cell}F_{R,b}
                    +F_{R,b}R_{\rm cell}F_{L,a})G.
\]

Then

\[
 H^{(2)}_{\rm adjacent}=-\lambda^2\sum_{j,ab}
 (G^\dagger F_{R,a}G)_{j-1}\,M_{j,ab}\,
 (G^\dagger F_{L,b}G)_{j+1}.                         \tag{11}
\]

In particular, one cannot infer this term from separate port compressions.

In completion C, inactive-port dressing does not change (8) or (10) at their displayed orders. A nontrivial bilinear response at an active port needs two insertions of its dressing g. Two insertions at the inactive port alone contract to sum_a eta_a^2=r and are scalar. Dressing both ports adds two more powers of g. The leading four-Majorana bond term therefore has two active-port insertions at each endpoint, giving precisely the one-port tensor at order g^4.

The possible order-g^2 contribution to M_ab cancels. For a!=b the two logical endpoint couplings gamma_a and eta_b commute, and at this order the corresponding auxiliary responses factorize. For a=b they anticommute, and the Hermitian cross response contains their vanishing anticommutator. The first nonzero M_ab is order g^4. With the two outside port compressions, (11) therefore starts at order lambda^2 g^6/Delta^7. This power counting applies to completion C, not an arbitrary cell with the same compressed ports.

I checked these cancellations and the full second-order formula numerically for r=2 and r=4. The code uses the entire auxiliary Hilbert space, including the dark states. For r=2, g/Delta=0.1, Delta=1, completion C gives

\[
 \alpha_C=-0.189030037122838,
 \quad \kappa_{12,12}/\lambda^2=0.003163787106,
\]

and the adjacent-bond contribution to the continuum coefficient is

\[
 G_{\rm adjacent}/(a\lambda^2)=0.00442661550526.
\]

The full coefficient is G/(a lambda^2)=0.0550472092013. It differs from using the isolated one-port answer at this finite g.

A stronger numerical check compares exact low-energy spectra with the full SW Hamiltonian at g/Delta=0.3. The maximum spectral errors divided by lambda^3 are

| Number of cells | lambda=0.004 | lambda=0.002 | lambda=0.001 |
| :-- | --: | --: | --: |
| 2 | 0.171706 | 0.171612 | 0.171566 |
| 3 | 0.318275 | 0.318021 | 0.317894 |

The three-cell check includes (11). The expected O(lambda^3) remainder is observed. This checks the perturbative operator calculation, not the thermodynamic infrared phase.

## 4. Continuum normalization and the sign convention

Use global Jordan–Wigner Majoranas Gamma_{j,a}, with {Gamma,Gamma}=2 delta, and cell spacing a. The first-order hopping is

\[
 H_1=t\sum_{j,a}i\Gamma_{j,a}\Gamma_{j+1,a},
 \quad t=\lambda\alpha^2.
\]

Normalize the continuum fields by {R_a(x),R_b(y)}=delta_ab delta(x-y), and likewise for L. Then

\[
 \Gamma_{j,a}=\sqrt{2a}\,[R_a(x)+(-1)^jL_a(x)]+\cdots,
 \quad v=4|t|a.
\]

For t>0 the continuum kinetic Hamiltonian is

\[
 H_0={iv\over2}\int dx\sum_a(R_a\partial_xR_a-L_a\partial_xL_a).
\]

The nonderivative staggered part of T_{j,a}=i Gamma_{j,a} Gamma_{j+1,a} is -4ia(-1)^j R_aL_a. Hence

\[
 \kappa\sum_{j,a<b}T_{j,a}T_{j,b}
 \longrightarrow 16a\kappa\int dx\sum_{a<b}R_aR_bL_aL_b.
\]

In the convention requested in the question,

\[
 G=16a\kappa.                                        \tag{12}
\]

To fix the RG sign without terminology ambiguities, set epsilon_a=iR_aL_a. Then R_aR_bL_aL_b=epsilon_a epsilon_b. With propagators normalized as 1/(2 pi z),

\[
 \epsilon_a(z)\epsilon_a(0)\sim{1\over4\pi^2|z|^2}.
\]

For O=sum_{a<b} epsilon_a epsilon_b,

\[
 O(z)O(0)\supset{2(r-2)\over4\pi^2|z|^2}O(0).
\]

The minus sign in the second cumulant of the Euclidean action gives

\[
 {dG\over d\ell}=-{r-2\over2\pi v}G^2+O(G^3).        \tag{13}
\]

Thus G>0 is marginally irrelevant for r>2, and G<0 is marginally relevant. In an isotropic theory the latter is the attractive interaction -|G| sum epsilon_a epsilon_b in the mass channel. The associated strong-coupling scale is exponentially small in v/|G|. Dynamical dimerization is compatible with exact translation symmetry because translation can break spontaneously.

For a general tensor kappa_AB, define

\[
 A_{ab}=R_aR_b,\quad B_{ab}=L_aL_b,\quad
 C_{ab}=R_aL_b-R_bL_a.
\]

The complete tree-level, normal-ordered quartic continuum limit of the neighboring-cell term is

\[
 4a\sum_{A,B}\kappa_{AB}
 [(A_A+B_A)(A_B+B_B)-C_AC_B].                        \tag{14}
\]

This includes mixed current interactions and pure chiral quartics. Short-distance contractions also renormalize the ground energy and quadratic kinetic coefficients. Those contractions do not change the displayed tree-level four-fermion amplitudes. Equations (9), (10), and (14) specify all leading interactions for any grading, without assuming O(r).

## 5. q=2, r=2

Take F_1=X and F_2=Y. Then K_ab=2 gamma_a gamma_b. For a!=b the coefficients in (4)-(5) reduce to

\[
 A={1-u^2\over r},\quad B={u^2\over r},\quad C=-{u\over r}.
\]

The two ordered choices a,b give

\[
 \kappa=-2\lambda^2\left[
 {2AB\over D}+{2AC\over d}+{B^2\over2D}
 +{2BC\over D+d}+{C^2\over2d}\right].
\]

Using d=D(1+u)/2 yields

\[
 \boxed{\kappa={\lambda^2u(1-u)^2(3u^3+18u^2+38u+24)
                  \over Dr^2(1+u)(3+u)}>0}
 \quad (0<u<1),\quad r=2.                           \tag{15}
\]

In particular,

\[
 \kappa={83\over2}{\lambda^2g^4\over\Delta^5}+O(g^6),
 \qquad G=664a{\lambda^2g^4\over\Delta^5}+O(g^6).
\]

The coefficient is even in g and lambda. It cannot be reversed by changing their signs. Taking g=0 or lambda=0 also removes the hopping, so this does not produce a nontrivial critical zero-interaction chain.

Direct diagonalization of the one-port two-cell resolvent gives kappa/lambda^2=0.0000329878154448, 0.00339014278434, and 0.0800818811814 at g/Delta=0.03, 0.1, and 0.3. Equation (15) agrees within 1e-16 in these tests. Completion C has the same leading coefficient but different higher powers of g.

For r=2, (13) vanishes. Combine the two Majoranas into one Dirac fermion. If rho_R=iR_1R_2 and rho_L=iL_1L_2, the interaction is -G rho_R rho_L. With exact U(1), this is the Abelian current interaction and changes the Luttinger parameter rather than generating a mass. To first order, after absorbing same-chirality velocity corrections,

\[
 K=1+{G\over2\pi v}+O(G^2/v^2)>1.
\]

A lattice sign check is useful. With f_j=(Gamma_{j,1}+i Gamma_{j,2})/2,

\[
 i\Gamma_{j,1}\Gamma_{j,2}=2n_j-1,
 \quad T_{j,1}T_{j,2}=-(2n_j-1)(2n_{j+1}-1).
\]

Positive kappa is attractive nearest-neighbor density interaction, consistent with K>1.

At commensurate filling the two-particle umklapp is cos(4 phi), not a dimension-two product of identical chiral fermions without derivatives. Its scaling dimension is 4K. The local product (R^dagger L)^2 vanishes without point splitting because there is only one Dirac flavor. At K=1 the nonzero lattice umklapp has dimension four. It is irrelevant for K>1/2 and is therefore irrelevant in the weak-coupling region found here.

Exact U(1) is a sufficient protection against pairing vertices. If it is relaxed, one-cell translation still forbids the dimension-one R-L pairing mass in the original Majorana momentum convention. A possible four-particle dual vertex has dimension 4/K and is irrelevant near K=1, but the full reduced-symmetry perturbation set must then be reconsidered. The clean conclusion here is an open c=1 Luttinger-liquid region with exact U(1), weak interactions, and no band-edge crossing. No critical-point tuning is required in that region.

## 6. q=3, r=4: anisotropy matters

Take the normalized odd basis X_1,Y_1,X_2,Y_2 coupling the one-dimensional plus block to each of the two minus states. In pair order (12,13,14,23,24,34), (9) gives

\[
 {\kappa\over\lambda^2g^4/\Delta^5}=
 \begin{pmatrix}
 29.5&0&0&0&0&-19.5\\
 0&54.25&0&0&-5.25&0\\
 0&0&54.25&5.25&0&0\\
 0&0&5.25&54.25&0&0\\
 0&-5.25&0&0&54.25&0\\
 -19.5&0&0&0&0&29.5
 \end{pmatrix}.                                    \tag{16}
\]

The bulk bilinear is

\[
 H_{\rm charge}=h\sum_j i(\gamma_{j,1}\gamma_{j,2}
                              +\gamma_{j,3}\gamma_{j,4}),
 \quad h=3{\lambda^2g^2\over\Delta^3}+O(g^4).        \tag{17}
\]

This term is allowed by U(2). It is not a forbidden momentum-pi mass.

Write the mixed continuum interaction as sum_AB M_AB R_A L_B and set U=a lambda^2 g^4/Delta^5. Equation (14) gives

\[
 {M\over U}=
 \begin{pmatrix}
 472&0&0&0&0&-240\\
 0&868&0&0&-156&0\\
 0&0&868&156&0&0\\
 0&0&156&868&0&0\\
 0&-156&0&0&868&0\\
 -240&0&0&0&0&472
 \end{pmatrix},                                    \tag{18}
\]

and the pure chiral term is

\[
 -72U(R_1R_2R_3R_4+L_1L_2L_3L_4).                  \tag{19}
\]

The O(4)-invariant projection of (18) is 736U. Reporting only that positive number would miss the instability.

Use the orthonormal self-dual combinations (12+34,13-24,14+23)/sqrt(2), and the anti-self-dual combinations with opposite signs. In this basis (18) is diagonal:

\[
 (g_z,g_\perp,g_\perp)=(232,1024,1024)U,
 \qquad (g_s,g_s,g_s)=(712,712,712)U.                \tag{20}
\]

The second sector is marginally irrelevant. The first obeys the axial SU(2) equations, to leading order at the common free velocity,

\[
 {dg_z\over d\ell}=-{g_\perp^2\over\pi v},\qquad
 {dg_\perp\over d\ell}=-{g_zg_\perp\over\pi v}.
\]

The invariant g_z^2-g_perp^2 is negative. Therefore at commensurate filling this sector flows to strong coupling. Positivity of all eigenvalues of the bare current matrix is not a stability criterion for an anisotropic non-Abelian interaction. The chiral terms change velocities at weak coupling and do not remove this leading BKT conclusion.

There is a necessary qualification. The transverse self-dual currents are chiral pairs of the two Dirac flavors. Their interaction transfers a pair between right and left Fermi points, so it is a two-particle umklapp. Equation (17) shifts the filling. After the gauge transformation f_j=i^j c_j, the quadratic dispersion for each flavor is

\[
 \epsilon(k)=-4t\cos(ka)+2h,
 \quad t={4\lambda g^2\over\Delta^2}+O(g^4).
\]

Hence |delta k_F|=|h|/(2|t|a)=3|lambda|/(8a Delta)+higher orders. The umklapp wavevector is 4 delta k_F. Its detuning is polynomial in lambda, whereas the commensurate strong-coupling scale is exponentially small in Delta^3/(|lambda|g^2). At fixed small g and sufficiently small lambda, detuning stops the umklapp flow first.

Thus completion C is consistent with a weak-coupling c=2 phase at its freely selected, shifted filling: one charge boson and the marginally irrelevant SU(2)_1 flavor sector. If the total charge is constrained to the commensurate value, or a counterterm is adjusted to restore that value, the detuning argument does not apply. The leading commensurate theory instead has a charge-sector instability and an expected remaining c=1 flavor sector. These are perturbative continuum conclusions, not a nonperturbative lattice phase proof.

## 7. The family q_+=1, q_-=m

This family has m complex fermions and an exact U(m) symmetry in completion C. The port basis need not, and generally does not, realize O(2m). Put q=m+1.

For clarity, first replace (9) by a general matrix

\[
 \kappa=A1+B S+Ctt^T.
\]

The leading coefficients are A=(119/2)lambda^2g^4/Delta^5, B=-(7/2)lambda^2g^4/Delta^5, C=lambda^2g^4/Delta^5. For Cartan pair labels H_i=(X_i,Y_i), direct commutators give

\[
 \kappa_{H_iH_i}=A+2qB-q^2C,
 \quad \kappa_{H_iH_j}=qB-q^2C\quad(i\ne j).
\]

For cross-flavor pair labels the diagonal is A+qB/2 and the paired off-diagonal entries are plus or minus qB/2. Substitution into (14) gives

\[
 M_{H_iH_i}=16aA+32aqB-16aq^2C,
 \quad M_{H_iH_j}=16aqB-8aq^2C.
\]

The uniform Cartan vector therefore has eigenvalue 16aA+16aq^2B-8aq^3C. A Cartan difference has eigenvalue 16aA+16aqB-8aq^2C. U(m) symmetry extends the latter to all SU(m) generators. The remaining two-form currents are the real and imaginary pair currents, with eigenvalue 16aA+8aq^2C. This derives, rather than fits, the following complete mixed-current spectrum:

\[
 \boxed{\begin{aligned}
 G_c/U&=952-8q^2(q+7), &&\text{multiplicity }1,\\
 G_s/U&=8(119-7q-q^2), &&\text{multiplicity }m^2-1,\\
 G_p/U&=8(119+q^2), &&\text{multiplicity }m(m-1).
 \end{aligned}}                                    \tag{21}
\]

Here U=a lambda^2g^4/Delta^5. For m=1 only the charge eigenvalue exists. The pure chiral coefficient on each four-field set belonging to two distinct complex flavors is -8q^2 U. Such terms are U(m)-invariant combinations of charge and flavor chiral stress tensors, and modify their velocities.

The generated common bulk chemical potential is

\[
 h=(m^2-1){\lambda^2g^2\over\Delta^3}+O(g^4).        \tag{22}
\]

At an unconstrained shifted filling, G_p becomes oscillatory. G_c is an Abelian forward interaction and can have either sign without opening a gap at sufficiently weak coupling. The remaining question is the SU(m) flavor coupling G_s. In this convention its negative sign is marginally relevant.

| m | q | G_c/U | G_s/U | G_p/U |
| --: | --: | --: | --: | --: |
| 1 | 2 | 664 | absent | absent |
| 2 | 3 | 232 | 712 | 1024 |
| 3 | 4 | -456 | 600 | 1080 |
| 4 | 5 | -1448 | 472 | 1152 |
| 5 | 6 | -2792 | 328 | 1240 |
| 6 | 7 | -4536 | 168 | 1344 |
| 7 | 8 | -6728 | -8 | 1464 |

For m>=7 the leading flavor coupling is negative. The common detuning in (22) does not make SU(m) flavor-exchange backscattering oscillatory, because all flavors still have the same Fermi momentum. The flavor sector therefore has a perturbative instability. A flavor gap with a remaining charge Luttinger liquid is the expected weak-coupling outcome away from commensurability, rather than c=m. The strong-coupling phase identification is an expectation from the continuum model, not established by the finite-cell numerics.

This is a concrete obstruction to taking the small-g, flavor-symmetric family to arbitrarily large m. At m=7 the negative coefficient is small, so g and lambda must be sufficiently small for the displayed leading sign to dominate corrections. For every fixed m>=7 such a regime exists within the perturbative expansion.

There is a further conditional finite-g result for the isolated one-port endpoints. Equations (7) and (9) retain the form A1+B S+Ctt^T, with A=-2lambda^2 f_JJ, B=-2lambda^2 f_JL, C=2lambda^2 f_LL. Their SU(m) coefficient is

\[
 {G_s\over a\lambda^2}=
 {-8x^2(2x-1)\over\Delta r^2(x-2)(x-1)}
 [q^2x-2q^2+16qx-14q+96x^4-480x^3+880x^2-736x+238].
 \tag{23}
\]

The prefactor is positive for 0<x<1/2. At q=8 the bracket is twice

\[
 48x^4-240x^3+440x^2-272x-1.
\]

On x in [0,1/2], its Bernstein coefficients are -1, -35, -152/3, -111/2, -54, all negative. The bracket decreases further with q>=8. Thus varying g/Delta cannot repair the isolated-endpoint SU(m) sign for q>=8. Equation (23) must not be promoted to an all-g statement about completion C or another unspecified mirror construction.

## 8. Relevant operators and exact symmetries

One-cell translation acts on the continuum fields as R_a -> R_a and L_a -> -L_a, up to the spatial shift. It forbids every iR_aL_b quadratic mass. Fermion parity forbids odd fermion polynomials. Neither statement forbids the marginal interactions above.

The full list of polynomial perturbations of dimension at most two near the free theory includes the following.

- Same-chirality bilinears iR_aR_b and iL_aL_b have dimension one. They are uniform flavor currents, and are allowed when the exact flavor group has an invariant generator. Equation (10) is an explicit example. They can shift Fermi points without gapping the system. Translation alone does not eliminate them.
- Quadratic derivative terms have dimension two and change velocities or mix chiral flavors.
- Quartics with two R and two L fields are marginal. Their full tensor, not its trace, must be followed under RG.
- Quartics with four R or four L fields are also dimension two when four distinct flavors are available. In the U(m) examples above they are chiral velocity terms.
- Quartics with an odd number of L fields are forbidden by one-cell translation at the original commensurate point.

The restriction to polynomial local fermion operators uses exact total fermion parity and the local even operator algebra after Jordan–Wigner. An isolated Ising spin or disorder field is not an additional allowed local, parity-even perturbation of this algebra. If the microscopic model does not realize this parity/locality structure, that exclusion must be rechecked.

Exact O(r), if it is actually implemented, forbids invariant flavor currents and restricts the nonchiral quartic to the isotropic Gross–Neveu interaction. For r=4, SO(4) alone also permits the epsilon-tensor distinction between its two SU(2) sectors. An O(4) reflection removes that distinction. For r=2, SO(2)=U(1) permits a chemical potential, whereas charge conjugation can forbid it.

The graded port space generally has U(q_+) x U(q_-) symmetry, with possible additional discrete symmetries. Orthogonality of its r odd operators does not imply a physical O(r) symmetry. Its fourth-order trace tensor is the obstruction, explicitly visible in (16).

No action of time reversal on the ports, logical Clifford algebra, or clock is supplied in the question. It therefore cannot be used to discard operators without an additional assumption. Ordinary number-conserving time reversal after the hopping gauge transformation permits chemical potentials and the interactions discussed here. Time reversal does not fix the Gross–Neveu sign. The real, reflection-symmetric continuum tensors above do not require time-reversal breaking.

## 9. What can be adjusted, and what counts as a tuning

Changing the sign of lambda or g does not reverse (8). A basis rotation within a fixed grading changes coordinates on the tensor, not its invariant spectrum. Changing the grading changes the symmetry and the charge term, but does not provide a general sign theorem. In particular, the positive isotropic projection in the q=3 example is not enough.

Ground-preserving modifications of excited cell energies can change the answer. For example, keep the one-port G and u fixed in the r=2 calculation but increase only the dark gap d. The expression before (15) still applies with independent positive D,d. As d tends to infinity,

\[
 \kappa\longrightarrow-2\lambda^2
 \left({2AB\over D}+{B^2\over2D}\right)<0,
\]

whereas the undeformed value is positive. By continuity a zero exists between them. Such a deformation is implemented abstractly by adding a multiple of Pi_d. It preserves the ground compression. This example proves that the compression does not determine the sign. It does not prove that this projector deformation is an allowed physical single-site term in the encoded chain, nor that its zero cancels adjacent-bond and higher-order terms in a two-port realization.

If exact O(r) is imposed and a fully calculated microscopic model lies strictly on the G>0 side, sufficiently small symmetry-preserving corrections stay in that basin. This is an inequality, not a critical tuning. Setting G=0 instead is generally a tuning of the renormalized coupling; no symmetry listed in the question forbids the O(r)-invariant quartic. Cancelling only its leading SW coefficient does not establish an exactly free infrared theory.

### A sufficient detuned construction with smaller symmetry

For q_+=1, q_-=m, completion C has U(m), hence all m Cartan U(1) charges. Assume the microscopic history implementation preserves these charges exactly and permits adding cell-local logical terms whose leading compression is

\[
 H_h=\sum_{j,i}h_i\,i\Gamma_{j,2i-1}\Gamma_{j,2i}
     =\sum_{j,i}h_i(2n_{j,i}-1).                     \tag{24}
\]

These assumptions concern the entire Hamiltonian, including its virtual states and clock implementation. First-order charge conservation alone is insufficient. A logical occupation term on the workspace commutes with the combined Cartan charge of completion C; realizing its cell-local physical implementation remains a construction requirement.

Choose distinct positive h_i, with all h_i and their spread much smaller than |t|, so that every band remains partially filled. The free Fermi momenta satisfy

\[
 \cos(k_{F,i}a)={h_i\over2t}.
\]

The dangerous number-conserving four-fermion processes then have phases

\[
 R_i^\dagger L_iL_j^\dagger R_j:
       \quad 2(k_{F,i}-k_{F,j})x,
\]
\[
 R_i^\dagger L_iR_j^\dagger L_j:
       \quad 2(k_{F,i}+k_{F,j})x
       \quad\text{modulo }2\pi/a.
\]

Distinct shifts eliminate the first resonance. A common shift to the same side of half filling eliminates the second. Choosing detunings polynomially small in the microscopic weak couplings cuts off their RG flow well before any exponentially small marginally generated gap.

A Cooper-pair transfer between different flavors can retain zero total momentum even when their fillings differ. Exact U(1)^m forbids it. Total U(1) or fermion parity alone is not the sufficient symmetry being invoked here.

After these exclusions, the nonoscillating weak interactions are forward density interactions and chiral kinetic terms. Bosonization gives a positive quadratic form in the m density and phase fields. For sufficiently weak couplings this form remains positive, giving m massless bosons and c=m. Higher local vertices have dimension greater than two near the free point, apart from the resonances already excluded. Accidental higher-order commensurabilities do not create a relevant operator arbitrarily close to the free point.

This is a controlled weak-coupling field-theory argument, not a rigorous many-body spectral proof. For each fixed m it describes an open region. One can choose a fixed distinct positive list of h_i ratios and vary one overall small amplitude within that region. There are m chemical-potential coefficients to implement, but no m equations that must be finely tuned to criticality. The first-order target is now m Dirac bands with distinct fillings, rather than coincident Majorana cones at 0 and pi.

## 10. Verdict and limitations

Under the stated translation and parity conditions alone, infrared stability of route B is not established. The exact flavor symmetry and virtual-state Hamiltonian matter.

For r=2, a weakly interacting c=1 phase is supported by the explicit calculation. With exact U(1), umklapp is irrelevant near the free point and the interaction is an Abelian Luttinger coupling. No critical tuning is needed within that weak-coupling phase.

For r=4 in completion C, the commensurate theory has a marginally relevant anisotropic sector. The generated common chemical potential can remove its umklapp instability when the ground state is free to select its filling. This distinction disappears if the charge is constrained to its commensurate value.

For the flavor-symmetric q_+=1, q_-=m family at small dressing, a non-Abelian marginally relevant coupling survives common detuning for m>=7. The first-order route is therefore not automatically infrared-stable at unbounded r. Neither translation nor the positive isotropic projection resolves this obstruction.

A sufficient proposed repair is exact U(1)^m together with distinct noncommensurate fillings and weak interactions. It uses an open parameter region rather than cancellation of a Gross–Neveu coefficient. Its symmetry and cell-local implementation must be supplied before it becomes a completed construction.

The following work remains necessary for an unconditional claim about the route in the question: (i) specify the joint two-port Hamiltonian and the microscopic action of its symmetries, (ii) implement any stabilizing chemical-potential terms and exact charges in the history cells, and (iii) control the remaining finite-clock and higher-order corrections inside the chosen stable region. Small clock errors are sufficient only after their symmetry content and size relative to the leading stability margins are controlled. Approximate symmetry alone does not exclude a marginally relevant perturbation at arbitrarily long distances.

I have not proved infrared stability for an unspecified mirror completion, at arbitrary g, for all gradings, or at finite clock coupling. No large-system numerical phase determination was attempted. All computations were local dense linear algebra or small-chain sparse diagonalization, and no pre-existing files or external sources were read.

## Reproducibility

All files listed here were created for this calculation.

- `c11_calculation.py` constructs the one-port and completion-C cells, spectral tensors, and bond SW operator.
- `leading_tensor.py` checks the leading q=2 and q=3 flavor tensors and charge terms.
- `continuum_leading.py` expands the q=3 tensor into chiral fields and its two SU(2) sectors.
- `general_um.py` checks the full continuum spectrum for the U(m) family.
- `chain_sw.py` includes adjacent-bond processes and reports finite-g continuum coefficients for r=2,4.
- `validate_sw.py` compares second-order SW spectra with exact two-cell and three-cell spectra.
- `calculation_output.txt`, `chain_output.txt`, `general_um_output.txt`, and `validation_output.txt` contain the numerical outputs.
- `PROGRESS.md` records the calculation and the distinctions between proved algebra, numerical checks, and continuum assumptions.

Run the Python scripts with Python 3, NumPy, and SciPy. The exact rational coefficients in the leading tensors are also derived in the text; the numerical checks are not their sole justification.
