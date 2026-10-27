[proved] Port-idle history-state cells admit an exact code and exact port compression, together with quantitative convergence of the entire boundary resolvent and the second-order bond. Two qualifications are necessary. Signed energy terms require a sufficiently large clock scale, and universality for an arbitrary target may require extra initialized workspace. The syndrome-sector targets needed in the stated finite-dimensional gadget admit a simpler, structurally protected construction.

[open] These results do not establish the proposed fixed-dimension, translation-invariant nearest-neighbor embedding. In particular, the legal-sector Hamiltonian, its full gap, the length-selection mechanism, and thermodynamic control still need proofs. An exact bound on boundary-response complexity depending only on local dimension is false, even at a fixed positive gap. A finite-precision bound on a sufficiently strong logical swap does follow from locality, as derived below.

[proved] This analysis is closed book. No network or files outside this directory were consulted. The scripts set BLAS thread counts to one before importing numpy or scipy. All norms below are operator norms unless specified otherwise. Each status label governs the following derivation or table until the next label.

[proved] 1. The history isometry and the clock gap.

Write L=T+1. If E embeds the initialized input space into the full input space, define

\[
 \eta\psi=L^{-1/2}\sum_{t=0}^{L-1}|t\rangle W_tE\psi,
 \qquad P=\eta\eta^\dagger,\quad Q=1-P.
\]

E acts only on the workspace, apart from identity factors on p,s. Below G is the code isometry on this initialized input space. Thus the termwise annihilation hypothesis means h_t W_t E G=0. Without extra ancillas, E can be the identity.

[proved] For every operator F on p, [F,U_t]=[F,Pi_bad]=0. Each clock matrix commutes with F as well. Term by term in the Hamiltonian,

\[
 [F,H_{\rm prop}]=[F,H_{\rm init}]=0,
 \qquad F\eta=\eta F.
\]

The same identity holds on the history space before imposing initialization. It also holds on the initialized history space because F preserves initialized workspace. In particular, F has no matrix elements from a history state to the orthogonal complement of the history space. This statement does not require F to be Hermitian.

[proved] Put H_0=H_prop+H_init and

\[
 V=H_{\rm pen},\quad
 A=\eta^\dagger V\eta=H_A,\quad
 M=\|V\|=L\max_t\|h_t\|,\quad
 B=\|A\|\leq\sum_t\|h_t\|\leq M.
\]

The clock-controlled unitary \(\mathcal W=\sum_t|t\rangle\langle t|\otimes W_t\) transforms H_0 into

\[
 \mathcal W^\dagger H_0\mathcal W
 =\frac{J_c}{2}\sum_{t=1}^{L-1}
 (|t\rangle-|t-1\rangle)(\langle t|-\langle t-1|)\otimes1
 +J_c|0\rangle\langle0|\otimes\Pi_{\rm bad}.
\]

Its kernel is exactly P. A uniform, deliberately conservative bound is

\[
 QH_0Q\geq\gamma Q,\qquad \gamma=\frac{J_c}{2L^2}.
 \tag{1}
\]

[proved] To check (1), first restrict to initialized workspace. The nonzero path eigenvalues are J_c[1-cos(k pi/L)], k=1,...,L-1. This follows by substituting cos[k pi(t+1/2)/L] in the difference equation and the two endpoint equations. The smallest is at least 2J_c/L^2, using 1-cos x >= 2x^2/pi^2 for 0<=x<=pi. In a bad-initialization sector, write a vector as x_0,...,x_{L-1}. Its energy divided by J_c is

\[
 \|x_0\|^2+\frac12\sum_{t=1}^{L-1}\|x_t-x_{t-1}\|^2.
\]

Since x_t=x_0+sum_{j=1}^t(x_j-x_{j-1}), Cauchy-Schwarz gives

\[
 \sum_t\|x_t\|^2
 \leq 2L\|x_0\|^2+L(L-1)\sum_{j=1}^{L-1}\|x_j-x_{j-1}\|^2
 \leq2L^2\left(\|x_0\|^2+\frac12\sum_j\|x_j-x_{j-1}\|^2\right).
\]

Combining the two sectors proves (1), including L=1 when the propagation term is absent. No assumption about how Pi_bad transforms under later gates is needed, since the initialization term is at t=0.

[proved] 2. Exact ground space, and the needed positivity qualification.

If h_t W_t E G=0 for every t, then

\[
 V\eta G=0,\qquad H\eta G=0.
 \tag{2}
\]

There is no ground-space rotation or normalization correction. However, (2) alone does not prove that these zero modes are the ground states.

[proved] If every h_t is positive semidefinite, then for every J_c>0,

\[
 \ker H=\eta\ker A=\operatorname{ran}(\eta G).
 \tag{3}
\]

Indeed, a vector has zero energy for H_0+V if and only if it belongs to both kernels. The first condition makes it eta psi. The second is equivalent to <psi,A psi>=0, and positivity of A identifies this with ker A. This also proves that for positive terms, A G=0 already implies termwise annihilation.

[proved] For signed h_t, introduce the blocks on P and Q,

\[
 C=QV\eta,\quad D=QHQ,\quad
 \xi=\|C\|\leq M,\quad d=\gamma-M.
\]

Then D>=dQ. Let A>=0 have kernel G and gap Delta. Equation (2) implies C G=0. If

\[
 d>0,\qquad \delta_0:=\frac{\xi^2}{d}<\Delta,
 \tag{4}
\]

then H>=0 and (3) holds. On the orthogonal complement of eta G, its gap is at least

\[
 g_*=
 \frac{\Delta+d-\sqrt{(\Delta-d)^2+4\xi^2}}2>0.
 \tag{5}
\]

To prove these statements, decompose a vector orthogonal to eta G as eta x+y with x perpendicular to G and y in Q. Its energy is bounded below by

\[
 \Delta\|x\|^2-2\xi\|x\|\|y\|+d\|y\|^2.
\]

The two eigenvalues of this scalar quadratic form give (4)-(5). Alternatively, completing the square in y gives the Schur complement A-C^dagger D^{-1}C, which is positive on G-perpendicular under (4). A convenient sufficient condition for g_*>=Delta/2 is

\[
 \gamma\geq M+\frac{\Delta}{2}+\frac{2\xi^2}{\Delta}.
 \tag{6}
\]

The verification uses 2xi ||x|| ||y|| <= (Delta/2)||x||^2+(2xi^2/Delta)||y||^2.

[proved] The unconditional ground-space assertion in the question is false for signed terms. Take T=1, no gates or initialization, a syndrome projector Q_Omega=1-|Omega><Omega|, and

\[
 h_0=2\Delta Q_\Omega,\qquad h_1=-\Delta Q_\Omega.
\]

Both terms annihilate the code and A=Delta Q_Omega. On any syndrome vector the full clock matrix is

\[
 \begin{pmatrix}J_c/2+4\Delta&-J_c/2\\-J_c/2&J_c/2-2\Delta\end{pmatrix},
 \qquad \det=J_c\Delta-8\Delta^2.
\]

It has a negative eigenvalue for J_c<8Delta, and an additional zero mode at equality. Thus either positivity of the individual penalties or a scale condition such as (4) is necessary.

[proved] There is a useful schedule-dependent improvement over xi<=M:

\[
 \xi^2=\left\|
 L E^\dagger\sum_t W_t^\dagger h_t^2 W_t E-A^2
 \right\|.
 \tag{7}
\]

This follows from C^dagger C=eta^dagger V^2 eta-A^2. The matrix inside the norm is positive semidefinite. In particular, if V preserves P, then xi=0 and the history-space resolvent is already exact at finite J_c whenever the indicated ground space is the ground space.

[proved] 3. Exact port compression and the full boundary resolvent.

Under either (3)'s positive-term hypothesis or (4), the ground isometry is exactly G_sim=eta G. Consequently,

\[
 G_{\rm sim}^\dagger F G_{\rm sim}=G^\dagger F G.
 \tag{8}
\]

This is exact for every J_c for which the stated ground space is valid. There is no factor 1/L. In particular, all first-order port matrix elements vanish when G^dagger F_a G=0.

[proved] Let m=q^2-1 and define S from a direct sum of m logical spaces to the target input space by

\[
 S(x_1,\ldots,x_m)=\sum_a F_aGx_a.
\]

The cancellation condition for the full traceless port basis implies

\[
 G^\dagger X_pG=\frac{\operatorname{tr}X_p}{q}1_N
 \quad\hbox{for every port matrix }X_p,
\]

and hence

\[
 G^\dagger F_aF_bG=\delta_{ab}1_N,
 \qquad S^\dagger S=1_{mN},\qquad G^\dagger S=0.
 \tag{9}
\]

This uses linearity for non-Hermitian X_p as well. Thus the block boundary function has the particularly simple form

\[
 R_A(z)=S^\dagger(A+z)^{-1}S,\qquad
 R_{\rm sim}(z)=S^\dagger\eta^\dagger(H+z)^{-1}\eta S.
 \tag{10}
\]

The second identity uses F_a eta G=eta F_a G, not merely spectral matching.

[proved] Take z with x=Re z>=0. At z=0 use the inverse on the excited subspace. Let

\[
 \delta(z)=\frac{\xi^2}{d+x},\qquad d>0,\qquad\delta(0)<\Delta.
\]

Then, in the norm of the entire mN by mN block matrix,

\[
 \boxed{\quad
 \|R_{\rm sim}(z)-R_A(z)\|
 \leq
 \frac{\delta(z)}{(\Delta+x)(\Delta+x-\delta(z))}.
 \quad}
 \tag{11}
\]

This is uniform on every fixed interval 0<=z<=Z, and in fact on the full closed right half-plane, with the stated reduced-inverse convention. A uniform bound on an interval containing an excited-state pole cannot be asserted. The same proof works for real negative z as long as d+z>0 and Delta+z>xi^2/(d+z).

[proved] To prove (11), restrict to the orthogonal complement of the exact ground space. Block inversion gives

\[
 \eta^\dagger(H+z)^{-1}\eta
 =\left[A+z-C^\dagger(D+z)^{-1}C\right]^{-1}.
\]

Set E_z=C^dagger(D+z)^{-1}C. Since ||(D+z)^{-1}||<=1/(d+x), we have ||E_z||<=delta(z). Also ||(A+z)^{-1}||<=1/(Delta+x). A Neumann-series bound, followed by the resolvent identity, gives

\[
 \|[A+z-E_z]^{-1}-(A+z)^{-1}\|
 \leq\frac{\delta(z)}{(\Delta+x)(\Delta+x-\delta(z))}.
\]

Compression with the isometry S proves (11). Thus all excursions into clock-excited and initialization-bad states have been included through the full block D. They were not discarded. Port idleness removes direct Q-components of the boundary source, but it does not remove the self-energy C^dagger(D+z)^{-1}C.

[proved] For other collections of boundary operators, the right side of (11) is multiplied by ||S_perp||^2, where S_perp is their concatenated excited source. Ground components give identical 1/z contributions in the two models and cancel in the difference. An individual a,b block can also be bounded using ||F_a|| ||F_b|| instead of (9).

[proved] 4. The two-cell second-order exchange.

For two possibly different cells with zero port compression, define

\[
 X=\sum_{a=1}^m(F_aG_A)\otimes(F_aG_B).
\]

Here X is a map from the tensor product of the two logical spaces. Equation (9) gives X^dagger X=m1. The source lies in the excited subspace of each individual cell. Thus

\[
 K_{2,A}=-X^\dagger(A_A\otimes1+1\otimes A_B)^{-1}X,
\]

and the full simulator expression replaces A_A,A_B by H_A,H_B and X by (eta_A tensor eta_B)X. These are the same second-order operators as the spectral formula in the question.

[proved] Set

\[
 \Delta_\Sigma=\Delta_A+\Delta_B,\quad
 d_2=\min(\gamma_A,\gamma_B)-M_A-M_B,\quad
 \delta_2=\frac{\xi_A^2+\xi_B^2}{d_2}.
\]

If d_2>0 and delta_2<Delta_Sigma, then

\[
 \boxed{\quad
 \|K_{2,\rm sim}-K_{2,A}\|
 \leq\frac{m\,\delta_2}{\Delta_\Sigma(\Delta_\Sigma-\delta_2)}.
 \quad}
 \tag{12}
\]

For identical cells, d_2=gamma-2M and delta_2=2xi^2/d_2.

[proved] For completeness, restrict the pair Hilbert space to the tensor product of the two exact ground complements. This subspace is invariant because each individual ground space is exact. The low block is P'_A tensor P'_B, where P'_j=P_j-G_sim,j G_sim,j^dagger. Its Hamiltonian is A_A tensor1+1tensor A_B and its gap is Delta_Sigma. The complementary clock block is bounded below by d_2. The low-to-high coupling has two orthogonal ranges, Q_A tensor P'_B and P'_A tensor Q_B, so its squared norm is xi_A^2+xi_B^2. Apply the same Schur-complement estimate at z=0 and multiply by ||X||^2=m. This proof includes states in which one or both clocks are excited.

[proved] The physical second-order energy error is lambda^2 times (12). This is a statement about K2, not an assertion that the full finite-lambda two-cell Hamiltonian equals its second-order truncation. For a finite pair, perturbation theory additionally requires |lambda| ||v|| small compared with the simulator cell gap. Here the port completeness identity gives v=q Swap_q-1 and ||v||=q+1 for q>=2.

[proved] 5. An explicit hierarchy and its limitations.

For identical cells, choose desired absolute errors epsilon_R and epsilon_K. The following sufficient prescription uses only the specified schedule norms:

\[
 D_{\rm req}=\max\left\{
 \Delta+\frac{2M^2}{\Delta},
 \frac{2M^2}{\epsilon_R\Delta^2},
 \frac{mM^2}{\epsilon_K\Delta^2}
 \right\},
 \qquad
 J_c\geq2L^2(2M+D_{\rm req}).
 \tag{13}
\]

It guarantees an exact ground code, a simulator gap at least Delta/2, and errors at most epsilon_R and epsilon_K in (11)-(12). Indeed, d_2>=D_req, d>=M+D_req, delta(0)<=Delta/2 and delta_2<=Delta. Therefore (11) is at most 2M^2/(d Delta^2), and (12) at most mM^2/(d_2 Delta^2). Replacing M by the actual xi in the quadratic terms improves (13).

[proved] With h_max=max_t||h_t||, the scale appearing in these estimates is L^2 M^2=L^4 h_max^2. Thus, at fixed schedule, both errors are O(L^4 h_max^2/(J_c Delta^2)), with the displayed constants specifying that notation. If a construction supplies h_max<=c_q B, this becomes O(c_q^2 L^4 B^2/(J_c Delta^2)). One cannot replace schedule norms by B alone without such a hypothesis. The signed two-time counterexample keeps B=Delta while making the individual canceling terms arbitrarily large.

[proved] The dimensions of (11)-(13) are consistent. R and K2 have units of inverse energy for dimensionless F_a and v. Delta,d,M,xi,J_c have units of energy. Consequently the required scale M^2/(epsilon Delta^2) has units of energy. When xi=0 the resolvent and K2 errors vanish, and when J_c tends to infinity at fixed L and h_t they vanish as 1/J_c.

[proved] 6. Universality on the syndrome sector.

Choose an orthonormal ps basis |0>=|Omega>, |a>=(F_a tensor1)|Omega>, a=1,...,m. Let Q_Omega=sum_a |a><a|. For a workspace of n qubits and dimension D=2^n, any Hermitian operator supported on syndrome tensor workspace has an expansion

\[
 H_{\rm syn}=\sum_{\mu=1}^{m^2}\sum_{\nu=1}^{D^2}
 c_{\mu\nu}K_\mu\otimes P_\nu.
 \tag{14}
\]

The K_mu are a Hermitian basis for matrices on the syndrome space, extended by zero on |Omega>. The P_nu are the n-qubit Pauli strings. All c_mu,nu are real, and K_mu|Omega>=0=<Omega|K_mu. The expansion follows by taking tensor products of two real bases of Hermitian matrices. No syndrome label is measured by the computation.

[proved] Every nonidentity Pauli string is V_nu^dagger Z_c V_nu for a workspace-only Clifford circuit V_nu. To construct it, rotate each nonidentity single-qubit Pauli into Z, and compute the parity of its support onto a chosen support qubit using CNOTs. Conjugation of that qubit's Z by the parity circuit is the product of the supported Z operators. A SWAP can move the chosen qubit to the designated c. The identity string uses Y_c=1. Hence each term in (14) is

\[
 V_\nu^\dagger h_{\mu\nu}V_\nu,
 \qquad h_{\mu\nu}=c_{\mu\nu}K_\mu\otimes Z_c
\]

or the corresponding identity term. Every such local energy term annihilates |Omega> tensor every workspace state. The terms may be signed even when their sum is positive.

[proved] There are at most m^2D^2 terms. Each conjugating circuit has O(n) elementary gates with unrestricted workspace connectivity, or O(n^2) gates with nearest-neighbor routing along the workspace. Run a circuit, put its energy term at the endpoint, and run its inverse before the next term. All other history times have h_t=0. If the chosen circuits have lengths ell_j, one may use

\[
 T\leq\sum_{j=1}^{M_{\rm terms}}(2\ell_j+1),
 \qquad
 T=O_q(D^2\max\{1,n^2\}).
 \tag{15}
\]

The added idle slot permits separate energy terms even when a circuit has zero length. Several terms with the same endpoint can instead be combined. Padding non-power-of-two workspace dimension into qubits changes these bounds only by constant factors, provided unused states are penalized appropriately.

[proved] The gadget auxiliaries and the FK initialization ancillas must be distinguished. A gadget auxiliary that participates in the target syndrome resolvent remains part of the target workspace. It must not be projected out by E. For the usual code G_w|i>=|i>_logical tensor|0...0>_aux, choose, for example,

\[
 H_A=Q_\Omega H_{\rm syn}Q_\Omega
       +\Pi_\Omega\otimes a\sum_j |1\rangle\langle1|_{{\rm aux},j},
 \tag{16}
\]

where the syndrome block is positive definite on all retained workspace, and a>0. Each auxiliary term can be routed to c using workspace SWAPs, with an exact return to the original frame. It annihilates the code term by term. The kernel of (16) is precisely |Omega> tensor logical tensor initialized gadget auxiliaries. An additional multiple of Q_Omega tensor1 may be used to set a syndrome gap. Thus the arbitrary syndrome matrices used in the finite-dimensional gadget supplied in the question are included without asserting unrestricted universality on an unchanged workspace.

[proved] For fixed q, the K_mu can be normalized so ||K_mu|| is bounded and each coefficient in (14) is at most c_q||H_syn||. For example, choose a Hilbert-Schmidt orthogonal basis on the fixed m-dimensional syndrome space. The coefficient is a normalized trace against K_mu tensor P_nu; the bound follows from ||K_mu||_1=O_q(1) and ||P_nu||_1=D. This is one setting in which h_max<=c_q B, after accounting for the separately chosen auxiliary penalties.

[proved] The ps pair cannot be replaced by a readable syndrome stored in s alone. For a qubit port, the four Bell states, including |Omega> and its three Pauli syndromes, all have reduced state 1_s/2. Every observable on s therefore has the same expectation in these four states, so no measurement on s can distinguish them. The construction above uses joint ps energy operators K_mu. It does not require a workspace gate to extract a syndrome label from s.

[proved] 7. Literal fixed-workspace universality fails, but initialized dilation repairs it.

Here is a counterexample when c must be one workspace qubit and no extra initialized workspace is allowed. Let the workspace have dimension four, G_w span |00>,|01>, and

\[
 H_A=Q_\Omega\otimes1_4
 +\Pi_\Omega\otimes\operatorname{diag}(0,0,1,2).
 \tag{17}
\]

This is positive with the required two-dimensional ground code. Suppose it were a sum W_t^dagger h_t W_t with every h_t W_t(|Omega> tensor G_w)=0. Define the Hermitian one-qubit matrix B_t=<Omega|h_t|Omega>. Then (B_t tensor1) W_t G_w=0. If B_t is nonzero it must have rank one, since its kernel must contain a two-dimensional workspace subspace. Its full workspace kernel therefore equals W_t ran G_w. Consequently

\[
 W_t^\dagger(B_t\otimes1)W_t=b_t\Pi_{G_w^\perp}.
\]

The Omega-Omega block of the entire sum would be a scalar multiple of Pi_{G_w-perp}, contradicting the two different nonzero eigenvalues in (17). Allowing either workspace qubit as c does not change this rank argument. For a fixed c dimension r, an analogous obstruction occurs at code rank (r-1)D/r when D/r>1. Thus a universal statement without a workspace-enlargement qualification is too strong.

[proved] With one extra initialized workspace qudit of dimension q^2+1, every positive target H_A on the original ps tensor workspace does have a representation with positive local energy terms, each annihilating the code. This assertion is about the initialized effective input, exactly as in the FK model. The full enlarged input contains extra states penalized by H_init.

[proved] To construct it, let the original workspace dimension be D and write

\[
 H_A=\sum_{j=1}^{R}\lambda_j|v_j\rangle\langle v_j|,
 \qquad \lambda_j>0,\quad R\leq q^2D.
\]

Every v_j is orthogonal to ran G. Introduce c of dimension q^2+1, initialize it to |0>, and retain a D-dimensional register r for the original workspace. For each j take a Schmidt decomposition

\[
 |v_j\rangle=\sum_{\ell=1}^{r_j}\sigma_{j\ell}
 |a_{j\ell}\rangle_{ps}|s_{j\ell}\rangle_w,
 \qquad r_j\leq q^2.
\]

Choose a workspace unitary V_j that sends the initialized states |0>_c|s_jell>_r to |ell>_c|0>_r. Complete the s_jell to a workspace basis and send the remaining D-r_j initialized basis vectors to |0>_c|k>_r, k=1,...,D-r_j. These destinations are orthonormal; an orthonormal completion extends the map to a unitary on c,r.

[proved] Define

\[
 |v'_j\rangle=\sum_{\ell=1}^{r_j}\sigma_{j\ell}
 |a_{j\ell}\rangle_{ps}|\ell\rangle_c,
 \qquad h_j=\lambda_j|v'_j\rangle\langle v'_j|_{ps,c}.
\]

Writing E for initialization of c, direct application to an arbitrary original input gives

\[
 (h_j\otimes1_r)V_jE
 =\lambda_j(|v'_j\rangle\otimes|0\rangle_r)\langle v_j|.
 \tag{18}
\]

The components outside the selected Schmidt subspace have c=0 and are killed by h_j. Equation (18) proves both

\[
 E^\dagger V_j^\dagger h_jV_jE=\lambda_j|v_j\rangle\langle v_j|,
 \qquad h_jV_jEG=0.
\]

Summing establishes the general representation. It also proves exact ground-space equality for every J_c>0, since these terms are positive. The extra local qudit has dimension five for q=2. This statement does not assert that a single two-level c suffices without further construction.

[proved] The schedule length is finite and polynomial in the original matrix dimension when continuous elementary gates are allowed. A unitary on (q^2+1)D workspace dimensions can be eliminated by O_q(D^2) two-level rotations and phases. Each two-level rotation between computational basis vectors can be implemented by a Gray-code sequence of controlled basis changes, a controlled one-qubit rotation, and the inverse sequence. With O(log D) clean scratch bits, decomposing the controls and routing on a workspace line costs O_q(log^3 D) elementary gates per rotation. The scratch bits are returned to their initialized states. Thus one explicit loose count is

\[
 \ell_j=O_q(D^2\max\{1,\log^3D\}),\qquad
 T\leq\sum_j(2\ell_j+1)
   =O_q(D^3\max\{1,\log^3D\}).
 \tag{19}
\]

The controlled rotations here have arbitrary angles. For a fixed discrete gate set, approximation cost and computability of the angles must be added. This abstract count does not prove that a single fixed microscopic translation-invariant interaction can program arbitrary real coefficients exactly.

[proved] 8. What approximate compilation preserves.

For terms K_ps tensor Y_c with K|Omega>=0=<Omega|K, every workspace-only compilation, however inaccurate as an approximation to the intended workspace unitary, still annihilates the entire Omega sector. Thus the code stays exact. For Pi_Omega tensor Y_c terms, exactness instead requires Y_c W_t G_w=0. An exactly reversible computation of a code-validity flag can enforce this condition. Its uncomputation must also be exact in the circuit model. Running a compiled circuit and its literal adjoint cancels exactly as a matrix identity.

[proved] A generic approximation to the V_j in the positive spectral dilation need not preserve (18)'s cancellation on G. That construction proves exact representability, not automatic protection against arbitrary compilation. Likewise, an arbitrary code G_w need not have a flag computable by exact classical reversible logic. The protected syndrome construction and simple initialization flags avoid this problem for the stated gadget class. Physical errors that perturb the protected ps matrices or the flag logic are outside this exact-compilation assertion.

[proved] If the compiled unitaries satisfy ||Wtilde_j-W_j||<=e_j, then with unchanged h_j,

\[
 \|\widetilde A-A\|\leq\epsilon_A:=2\sum_j\|h_j\|e_j.
 \tag{20}
\]

This follows by adding and subtracting Wtilde_j^dagger h_j W_j in each conjugation. If the structural annihilation condition remains exact and epsilon_A<Delta, the target code stays exact and the compiled target gap is at least Delta-epsilon_A. In addition to the FK error, the boundary-resolvent error from compilation is at most

\[
 \frac{\epsilon_A}{(\Delta+x)(\Delta-\epsilon_A+x)}.
\]

For two identical compiled targets the corresponding K2 error is at most m epsilon_A/[2Delta(Delta-epsilon_A)]. These bounds use the same sources and follow directly from the resolvent identity. If the annihilation condition is broken, neither exact ground-space equality nor exact first-order cancellation follows from merely increasing J_c.

[proved] 9. An unavoidable clock tower in this FK model.

Termwise annihilation gives more than an exact ground state. For each logical vector psi, the L-dimensional subspace spanned by |t> W_t E G psi has H_pen=0 identically and H_init=0 at its initialized endpoint. Propagation restricts to the ordinary path Laplacian. Thus the simulator contains the exact energies

\[
 J_c[1-\cos(k\pi/L)],\qquad k=0,\ldots,L-1,
 \tag{21}
\]

each with at least the logical multiplicity. For L>=2,

\[
 \operatorname{gap}(H)\leq J_c[1-\cos(\pi/L)]
 \leq\frac{\pi^2J_c}{2L^2}.
\]

When H has the claimed exact ground space, the variational principle also gives gap(H)<=Delta. The clock states in (21) are invisible to a port source from the ground when G^dagger F_aG=0, since their input component lies in G rather than F_aG. They nevertheless limit the full cell gap. In a faithful local realization where J_c is comparable to the local propagation-term norm, keeping a fixed gap by increasing J_c increases the microscopic norm. After normalization, the gap-to-norm ratio is at most a constant times L^{-2}. The clock bound itself is exact; relating J_c to local term norms is an additional requirement on the embedding.

[proved] 10. Explicit failures when a gate acts on p.

For a compression failure, use a single logical qubit w, T=1, and U_1=Swap_{p,w}. Take the input code |Omega>_ps tensor|psi>_w. Choose h_0=(Delta/2)Q_Omega and h_1=U_1 h_0 U_1^dagger, both allowed to act on p,s,w. Then A=Delta Q_Omega, both penalties are positive and annihilate their code frames, and the history code is exact for every J_c>0. For any traceless port Pauli F,

\[
 G_{\rm sim}^\dagger F_pG_{\rm sim}
 =\frac12G^\dagger(F_p+F_w)G=\frac12F_{\rm logical}.
 \tag{22}
\]

The error is independent of J_c. For two such cells the first-order exchange is already nonzero, namely (lambda/4) sum_a sigma_a tensor sigma_a on the logical qubits.

[proved] A still smaller spectral counterexample uses T=1 and U_1=Z_p, with an arbitrary spectator logical space. Put h_0=(Delta/2)Q_Omega and h_1=Z_p h_0 Z_p. In the propagation frame,

\[
 H=\frac{J_c}{2}(|1\rangle-|0\rangle)(\langle1|-\langle0|)
       \otimes1+1_{\rm clock}\otimes\Delta Q_\Omega.
\]

The target is A=Delta Q_Omega and the exact port compression is zero. However, the X and Y physical port sources are in the antisymmetric clock mode, while the Z source is in the uniform clock mode. Consequently

\[
 R_{\rm sim}(z)=
 \operatorname{diag}\left(\frac1{\Delta+J_c+z},
                         \frac1{\Delta+J_c+z},
                         \frac1{\Delta+z}\right)\otimes1_N,
\]

whereas R_A(z)=1_3 tensor1_N/(Delta+z). The error tends to 1/(Delta+z), not zero. The two-cell operators are

\[
 K_{2,\rm sim}=-\left[\frac1{\Delta+J_c}+\frac1{2\Delta}\right]1_{N^2},
 \qquad K_{2,A}=-\frac3{2\Delta}1_{N^2}.
 \tag{23}
\]

Their difference tends to 1/Delta in norm. This example is scalar on the logical space, but it already disproves resolvent preservation from spectra alone. The random example below has a nonzero connected logical interaction.

[proved] The dilution comparison requires counting affected history frames, not just gates that touch p. For a general circuit,

\[
 \eta^\dagger F\eta=\frac1L\sum_t W_t^\dagger F W_t.
\]

If only r of the cumulative W_t fail to commute with F, its difference from F is at most 2||F|| r/L. On the input code one also has

\[
 \|F\eta G-\eta FG\|\leq2\|F\|\sqrt{r/L},
\]

because the squared norm is a sum over orthogonal clock frames. A single port gate can affect every subsequent W_t, so r need not be one. Padding with genuinely port-restored idle frames can reduce these errors. Raising J_c at fixed schedule cannot reduce their limiting history-space part. For port-idle circuits r=0 identically, which removes this dilution mechanism.

[numerical] 11. Exact diagonalization.

The reproducible main experiment is experiment.py. It uses q=2, two workspace qubits, T=5, and random Haar workspace gates with seed 20261027. One workspace qubit is logical and one is an FK initialization ancilla. At each time h_t=K_t tensor Y_t tensor1_anc, where K_t is a positive random three-by-three syndrome matrix extended by zero on Omega, and Y_t is a positive one-qubit matrix. The full Hilbert-space dimension is 96, the initialized target dimension is eight, and the code dimension is two. The observed target gap and norm are

\[
 \Delta=1.10482463097,\qquad B=3.38030172496,\qquad M=6.
\]

The target K2 eigenvalues are approximately -1.08482110, -0.71471869, -0.71250776, -0.52917274. After subtracting the scalar and both one-cell terms by partial traces, the remaining interaction has norm 0.04784196051.

[numerical] R errors below are norms of the complete six-by-six block boundary matrix, maximized over z in {0.05,0.2,1,3}. K2 is evaluated using every excited eigenstate of each full simulator. No excited-state truncation is used. The source amplitudes x[e,a,i]=<e|F_a|G_i> are contracted as v[e,f,i,j]=sum_a x[e,a,i]x[f,a,j], and K2=-sum_ef v_ef^dagger v_ef/(E_e+E_f). Thus the calculation includes clock-excited intermediate states directly.

[numerical] The idle-port results are:

| J_c | Ground dimension | Ground splitting | Port compression error | R error | K2 error | Full gap |
| --- | --- | --- | --- | --- | --- | --- |
| 10 | 2 | 2.86e-15 | 1.91e-17 | 2.23910e-1 | 2.02542e-1 | 0.340742 |
| 100 | 2 | 5.03e-14 | 1.91e-17 | 4.79501e-2 | 5.43245e-2 | 1.047129 |
| 1000 | 2 | 3.52e-12 | 1.91e-17 | 6.22033e-3 | 7.78059e-3 | 1.098217 |
| 10000 | 2 | 3.53e-11 | 1.91e-17 | 6.48094e-4 | 8.22875e-4 | 1.104153 |
| 100000 | 2 | 3.75e-10 | 1.91e-17 | 6.50921e-5 | 8.27810e-5 | 1.104757 |
| 1000000 | 2 | 3.11e-10 | 1.91e-17 | 6.51209e-6 | 8.28309e-6 | 1.104818 |

[numerical] Fits to the last three points give powers -0.998959 for R and -0.998571 for K2. The history-source leakage is 1.26e-16 and the largest termwise code residual is 2.78e-17. Ground multiplicity is counted with an absolute energy threshold of 1e-7, far below the smallest observed excited energy. The analytic splitting is zero. The displayed nonzero splittings are diagonalization roundoff, which grows with the large clock matrix entries. At lambda=10^{-3}, the second-order energy error at J_c=10^6 is 8.28309e-12. No claim about higher orders in lambda is inferred from this coefficient calculation.

[numerical] An independent diagonalization of the 64-dimensional two-target-cell Hamiltonian reproduces the spectral-contraction K2 to 2.40e-15. The measured coupling norm is xi=2.38223851520. The script also evaluates (11)-(12) and verifies the bounds at each sampled J_c where their hypotheses hold; those bound values are saved in numerics.json.

[proved] The active-port comparison inserts one Z_p gate at step three and conjugates h_t by Z_p for t>=3. These energy terms still have the allowed support and still annihilate the history code. The effective target A is unchanged. In fact the entire active Hamiltonian is unitarily conjugate to the idle Hamiltonian at every J_c, using a clock-controlled Z_p on its last three clock frames. Thus the spectra, ground multiplicities, and gaps agree. The physical port observables have not been conjugated along with the Hamiltonian, so their responses can differ.

[numerical] The active-port ground splitting, dimension, gap, and port-compression error equal the corresponding entries above to displayed precision. Its boundary source has history-space leakage norm one. Its errors are:

| J_c | R error relative to target | K2 error relative to target | R error relative to active limit | K2 error relative to active limit |
| --- | --- | --- | --- | --- |
| 10 | 0.426597 | 0.328744 | 3.85856e-1 | 3.66789e-1 |
| 100 | 0.697203 | 0.626044 | 6.41801e-2 | 6.87684e-2 |
| 1000 | 0.747058 | 0.686128 | 7.15246e-3 | 8.04607e-3 |
| 10000 | 0.752510 | 0.693217 | 7.25648e-4 | 8.22837e-4 |
| 100000 | 0.753061 | 0.693942 | 7.26739e-5 | 8.24793e-5 |
| 1000000 | 0.753116 | 0.694015 | 7.26845e-6 | 8.24983e-6 |

[proved] There are three unrotated and three Z-rotated frames. Therefore the history-projected X and Y sources cancel, while the Z source is unchanged. The limiting boundary function has only its ZZ source block, and the limiting K2 uses only the Z source. The complementary clock sources have energies diverging with J_c, so they give no limiting contribution.

[numerical] Evaluating those limiting operators independently gives errors 0.753121915807 for R and 0.694023134120 for K2. The residual convergence to those limits has fitted powers -0.999642 and -0.999434. The fit of the uncorrected errors instead has powers 0.000175 and 0.000250, consistent with nonzero limits.

[numerical] checks.py separately verifies the signed two-time counterexample, constructs the positive ancillary dilation for a random target with Omega/syndrome mixing, and checks the pole formula used below. The dilation reconstruction error is 1.67e-14, the largest termwise code residual is 2.01e-15, and the largest unitary residual is 1.46e-15. These checks support the formulas but do not replace their proofs. Raw outputs are in numerics.txt, numerics.json, and checks.txt. Both scripts run with python3 and only numpy/scipy plus the standard library.

[open] 12. Audit of the proposed translation-invariant embedding.

The list below separates consequences of a faithful embedding from steps that remain unproved here. Calling a step routine describes its expected engineering difficulty, not a proof that the proposed Hamiltonian already has it. "Open" here means not established from the supplied assumptions or this analysis; it does not assert that every listed proof obligation is a new unsolved problem in the literature.

[proved] Static-register identity is sufficient for exact port idleness. If every propagation, clock, initialization, and legal-configuration term acts as identity on the static p and s registers, then all of them commute with every physical port F. Transporting a workspace qudit to the boundary is compatible with this condition. Transporting the port or its partner is not. A legal-sector encoding must also intertwine the physical F with the abstract F, with no rescaling.

[open] The proposed computation construction has not been specified by a local alphabet and a complete list of nearest-neighbor terms. One needs an invariant legal sector unitarily equivalent to a path clock tensor input, with its input-independent history length and correct endpoint conditions. Multiple heads, no head, multiple counters, cycles, malformed programs, overlapping clocks, premature halting, and histories cut off by a marker all require analysis. A proof for the accepting computation sector alone does not give the full gap. An abstract reference to a Gottesman-Irani-type construction does not supply these missing properties.

[heuristic] Local constraints for a single head and a well-formed sweep, exact classical reversible counters, and a finite library of workspace gates are standard tasks in a history-state design. They are plausible ingredients. Their compatibility with this particular static-port geometry, with frustration freeness, and with a uniform lower bound on all illegal-sector energies still needs to be checked in one explicit construction.

[open] Legal-state constraints must permit every port syndrome, not just |Omega>. If an allegedly harmless legality projector pins the static pair to |Omega> inside H_0, it fails to commute with F and sends a boundary source outside the legal history sector. That would invalidate the port-idle resolvent proof. Pair initialization belongs in the energy construction, with the syndrome degrees of freedom retained in the simulated input, rather than being hidden in a port-dependent clock constraint.

[heuristic] The proposed two-site boundary geometry can support a local h on p,s,c if p is on the marker site, s is on its neighbor, and c is a computation-track register on one of these sites. The term is then a two-site operator even though it acts on three registers. A mirror arrangement can do the same at the other end. The complete selector rules must distinguish left and right boundaries without creating forbidden orientations or double assignments. Unused static p and s registers also need local pinning terms, otherwise they introduce extra ground degeneracy.

[proved] There is no abstract two-port obstruction. Use a code Omega_L tensor Omega_R tensor G_w and a sum of a left syndrome Hamiltonian and a right syndrome Hamiltonian, each supported on its own ps pair and the shared workspace, together with workspace-initialization penalties. Arrange each syndrome Hamiltonian to vanish on its own Omega sector and preserve that decomposition. On a left-only syndrome source the right Hamiltonian vanishes, and conversely. The two desired boundary resolvents can therefore be programmed independently in those sectors. Positivity on the sector where both ports are excited can be imposed by taking the syndrome Hamiltonians positive. This abstract observation does not furnish the local boundary schedule in a chain.

[open] A fixed local interaction has only a finite menu of matrix elements. The schedule in (14) or (18) contains target-dependent coefficients and, in general, arbitrary angles. A uniform algorithm can print finite descriptions of computable targets and compile them to a fixed gate library, but this is a separate programmable-Hamiltonian construction. It cannot print arbitrary uncomputable real matrices. Repetition, rational weights, controlled flags, or further simulation gadgets must account quantitatively for coefficients, their precision, runtime, and energy scales. Exactness of code annihilation must survive this implementation, as in the protected forms above.

[proved] Coding N through b only supplies classical information if an explicit map b -> N has been fixed. A b-site cell of local dimension d cannot have code dimension greater than d^b, so b>=log N/log d. Likewise, a legal path with T+1 mutually orthogonal configurations requires enough track Hilbert space to hold them. The abstract polynomial-in-N circuit count can therefore correspond to a runtime exponential in b, and a very small clock gap. No claim of polynomial runtime in b follows from (15) or (19).

[open] Frozen, conserved markers define invariant sectors; they do not energetically select one of them. A length-selection proof must compare the desired b(N) sector against all other lengths, mixed-length tilings, marker-free configurations, adjacent markers, and any illegal selector words. It must specify where the N-dependent parameter lives. A frustration-free accepting computation with zero energy at many lengths does not by itself favor one length. Boundary terms and a marker chemical potential can change energy per cell and per site differently, so both comparisons matter. No such calculation has been provided here.

[open] Even after a preferred isolated-cell length is found, coupling cells changes their energy by order lambda^2 times K2 per boundary. This can reorder near-degenerate length sectors. The length-selection margin must dominate that shift and higher-order corrections, or the selection calculation must include them. Interfaces between different selected lengths and marker placements require a domain-wall estimate. Superselection of marker positions does not replace this energetic argument for the unrestricted ground state.

[proved] The factor L in H_pen is an energy scale. In a uniform history, a physical term of norm j applied at a single clock time contributes only j/L to the input effective Hamiltonian. Port idleness eliminates dilution of the boundary operator, not this time-averaging of energy. Implementing a fixed nonzero h_t with L h_t as a bounded-strength local microscopic term is impossible for unbounded L unless the target scale is reduced, the schedule is redistributed, or the microscopic coupling is allowed to grow.

[proved] Equations (13) and (21) expose the normalization cost. The large J_c needed for accurate simulation can be chosen in a family of Hamiltonians with N-dependent strengths. After normalizing the largest microscopic bond norm, the cell gap and the target bond scale generally shrink with T and N. In a single fixed Hamiltonian, one cannot choose a new bare J_c independently for each of infinitely many cell lengths. A construction claiming all lengths at once must generate the scale separation internally or reduce the effective target scale with length.

[open] The path-clock lower bound is not a bound on every sector of a translation-invariant computation Hamiltonian. Illegal configurations can have their own long paths or zero modes. The abstract proof requires a lower bound on the entire orthogonal complement that a physical boundary operator or a local interaction can reach. If the encoding is not an exact invariant sector, its leakage and its effect on R(z) need an additional Schur-complement analysis. Spectral matching within a selected low-energy band is insufficient.

[proved] Exact first-order cancellation survives an exact embedding satisfying the static-port intertwining and termwise annihilation conditions. Approximate target spectral matching alone does not imply it. If compilation breaks annihilation, a small residual first-order term scales as lambda, whereas the desired exchange scales as lambda^2. It must therefore be controlled relative to the latter scale, not merely relative to the cell gap.

[open] A global finite-volume perturbation estimate using ||sum bonds lambda v|| grows with the number of cells and does not establish a thermodynamic effective Hamiltonian. One needs a local block-diagonalization or linked-cluster argument with constants independent of chain length. Its output should specify a local interaction norm or energy-density error. In such an argument one expects higher-order local terms of size O(lambda^3/g^2), with fixed coordination and bounded port interactions, but that estimate has not been proved for the proposed embedding here. An extensive O(number of cells) total error is compatible with a controlled energy density, not with a uniform global operator-norm error.

[open] A small local Hamiltonian error is not by itself a proof of the intended long-distance theory. A gapless logical target can be destabilized by a symmetry-breaking or relevant perturbation of arbitrarily small fixed strength. If the eventual goal is a particular SU(N) chain or continuum theory, one needs exact protecting symmetries, a stability result, or additional tuning. Taking J_c to infinity before the thermodynamic limit is a family of limiting models, not automatically one fixed microscopic chain with the desired limit.

[open] The substantive remaining embedding problems are therefore an explicit legal/illegal-sector construction with static ports, a programmable local energy interface at controlled precision, global energetic selection of b(N), and a uniform thermodynamic argument. Counter arithmetic, workspace routing, and finite gate compilation appear routine once these conditions are fixed. None of the abstract results above establishes the complete proposed translation-invariant construction.

[proved] 13. Why an exact bound depending only on d is false.

Even fixed gap and bounded bond norm do not bound the exact number of visible boundary poles. Here is a translation-invariant, nearest-neighbor bulk example with local dimension four. At each site regard the physical space as a qubit pair p,s, and use the basis |0>=|Omega>, |a>=sigma_a|Omega>, a=1,2,3. Let n=1-|0><0|. On a length-b open cell put

\[
 H_b=\mu\sum_{j=1}^b n_j
 -t\sum_{j=1}^{b-1}\sum_{a=1}^3
 \bigl(|a,0\rangle\langle0,a|+|0,a\rangle\langle a,0|\bigr)_{j,j+1},
 \qquad \mu>2|t|>0.
 \tag{24}
\]

On-site terms are included in the usual range-at-most-one convention; they can also be assigned to neighboring bonds with the corresponding open-end terms. The same local terms are used for every b. The vacuum |0>^{tensor b} is the unique zero-energy state and has zero compression for all three Pauli operators on the first p register.

[proved] On each bond the hopping operator has norm |t| on the one-excitation subspaces and vanishes elsewhere. Its negative part is bounded by |t|(n_j+n_{j+1}). Therefore

\[
 H_b\geq(\mu-2|t|)\sum_j n_j,
\]

which proves uniqueness of the vacuum and a gap at least mu-2|t| uniformly in b. The microscopic term norms are also bounded uniformly in b. A port Pauli at the first site creates one excitation of its corresponding species. In that sector the Hamiltonian is the tridiagonal hopping matrix, with energies and boundary weights

\[
 E_k=\mu-2t\cos\frac{k\pi}{b+1},\qquad
 w_k=\frac{2}{b+1}\sin^2\frac{k\pi}{b+1},\qquad k=1,\ldots,b.
\]

Substitution of the normalized sine eigenvectors proves these formulas. All w_k are positive. Hence

\[
 R_{ab}(z)=\delta_{ab}\sum_{k=1}^b\frac{w_k}{E_k+z}.
 \tag{25}
\]

For t nonzero the b energies are distinct, so there are b visible poles. There is no bound f(d) on this exact analytic complexity, even with a unique ground state, exact zero port compression, and a gap-to-local-norm ratio bounded below.

[proved] Code dimension is an even weaker candidate for such a bound. Tensor a free r-dimensional register into each site of (24), with the Hamiltonian acting as identity on those registers. The local dimension 4r stays fixed, the cell code dimension is N=r^b, the gap stays positive, the port compression is still zero, and R_ab is (25) times 1_N. These logical degrees of freedom are spectators, so this example does not claim a programmable swap on them. It shows why code dimension cannot be identified with the number of independently accessible logical observables. Conserved marker sectors can delimit these open cells using a fixed selector alphabet, but this observation does not establish selection of one preferred b.

[proved] 14. A port-strength tradeoff that does not require locality.

There is a useful necessary condition on a swap gadget even for an arbitrary cell matrix. Define the operator Schmidt nuclear norm of a two-logical-cell operator as the sum of its operator Schmidt coefficients, using Hilbert-Schmidt orthonormal operator bases on each side. For a tensor product X tensor Y this norm is ||X||_HS ||Y||_HS, and it is subadditive.

[proved] Under zero port compression, (9) implies for the imaginary-time boundary block matrix

\[
 R(t)=S^\dagger e^{-tH_A}S,
 \qquad 0\leq R(t)\leq e^{-\Delta t}1_{mN}.
\]

Using the representation in the question and subadditivity,

\[
 \begin{aligned}
 \|K_2\|_{{\rm Schmidt},1}
 &\leq\int_0^\infty\sum_{a,b}\|R_{ab}(t)\|_{\rm HS}^2\,dt\\
 &=\int_0^\infty\operatorname{tr}R(t)^2\,dt
 \leq\frac{mN}{2\Delta}.
 \end{aligned}
 \tag{26}
\]

The middle equality uses R_ba=R_ab^dagger. For different cells, Cauchy-Schwarz gives mN/(Delta_A+Delta_B) instead. This proof includes all target excitations, regardless of auxiliary dimension.

[proved] In a Hermitian orthonormal logical operator basis with T_0=1/sqrt(N), the swap is P_N=sum_mu T_mu tensor T_mu. Thus alpha 1+beta P_N has operator Schmidt coefficients |alpha N+beta| and N^2-1 copies of |beta|. Equation (26) yields, for N>=2,

\[
 |\alpha N+\beta|+(N^2-1)|\beta|
 \leq\frac{mN}{2\Delta},
 \qquad
 \boxed{\quad |\beta|\Delta\leq\frac{mN}{2(N^2-1)}.\quad}
 \tag{27}
\]

For a qubit port m=3. This is compatible with the arbitrary-matrix universality supplied in the question, because that statement does not fix both the target gap and the swap coefficient uniformly in N. It rules out a nonvanishing beta Delta at unbounded N with a fixed port. Adding a large scalar does not evade the inequality. The same conclusion applies to the simulator with its appropriate excitation gap, or a lower bound on the gap in the boundary-source spectral support.

[proved] Negativity of K2 gives a slight improvement. Since alpha<=-|beta|, (27) implies beta Delta<=mN/[2(N-1)(N+2)] for beta>0, and |beta|Delta<=m/[2(N+1)] for beta<0. The main large-N conclusion remains a coefficient at most of order 1/(N Delta).

[proved] 15. A finite-precision locality statement for a strong logical swap.

A meaningful further bound must specify precision and the strength of what is visible. The following conditional form makes those quantities explicit. Let two disjoint local cells have exact product ground code G_A tensor G_B, a gap at least g, and a local boundary interaction v of norm nu with zero ground compression. Assume their uncoupled real-time evolution obeys the locality estimate

\[
 \|\tau_t(v)-\tau_t^{(\ell)}(v)\|
 \leq C\nu e^{-\mu(\ell-u|t|)}.
 \tag{28}
\]

Here tau^(ell) evolves with the Hamiltonian truncated to ell physical sites on each side of the boundary, C,mu,u are independent of cell length, and u is a propagation speed. Finite-range bounded interactions have an estimate of this form, with u proportional to the local interaction norm and constants depending on range and coordination. One obtains it by iterating the commutator evolution equation: only chains of overlapping interaction supports can connect separated supports, each commutator contributes at most twice the term norm, and the time-ordered integrals supply 1/n!. Summing the path bound with an exponential distance weight gives (28). The result below needs only the displayed inequality and does not assume an area law or ground-state uniqueness within the code.

[proved] There is an operator O_ell supported on those 2ell sites such that its compression approximates the exact K2 with error at most

\[
 \frac{\nu^2\omega}{g^2}
 +\frac{\nu^2}{\omega}\left[
 C e^{-\mu\ell/2}+e^{-\omega\ell/(2u)}\right]
 \qquad (\omega>0).
 \tag{29}
\]

To prove this, put H=H_A tensor1+1tensor H_B and write G_pair for its ground isometry. The source v G_pair is orthogonal to the ground space. On its excited support,

\[
 \|H^{-1}-(H+i\omega)^{-1}\|\leq\omega/g^2.
\]

The full regulated inverse has the convergent real-time representation

\[
 (H+i\omega)^{-1}=-i\int_0^\infty e^{-\omega t}e^{iHt}\,dt.
\]

Since H G_pair=0, the compression of v e^{iHt}v equals that of v tau_t(v). Truncate the integral at tau=ell/(2u), and replace tau_t(v) by tau_t^(ell)(v). The discarded integral has norm at most nu^2 e^{-omega tau}/omega. On the retained interval (28) bounds the replacement error by C nu^2 e^{-mu ell/2}/omega. The resolvent regularization supplies the first term in (29). Multiplying the truncated integral by i gives O_ell for the sign convention K2=-G_pair^dagger v H^{-1}v G_pair. Hermitizing O_ell cannot increase its error from the Hermitian K2 and does not change its support.

[proved] Suppose now K2=alpha1+beta P_N, beta nonzero, and fix 0<epsilon<1. Define eta_vis=|beta|g/nu^2. Set

\[
 \omega=\frac{\epsilon|\beta|g^2}{3\nu^2}
       =\frac{\epsilon\eta_{\rm vis}g}{3}.
\]

Choosing an integer ell at least

\[
 \max\left\{
 \frac{2}{\mu}\log\frac{9C}{\epsilon^2\eta_{\rm vis}^2},
 \frac{6u}{\epsilon\eta_{\rm vis}g}
       \log\frac{9}{\epsilon^2\eta_{\rm vis}^2}
 \right\}
 \tag{30}
\]

with negative logarithms replaced by zero makes (29) at most epsilon|beta|. A trivial larger choice can always include the physical port support. If a cell is shorter than this ell, its whole Hilbert space already has dimension at most d^ell.

[proved] Every operator supported on ell sites in each cell has operator Schmidt rank at most d^(2ell). Compression by separate cell isometries cannot increase this rank. Subtracting alpha1 adds rank at most one. On the other hand, beta P_N has N^2 equal operator Schmidt coefficients |beta|. A rank-r approximation to it has squared Hilbert-Schmidt error at least |beta|^2(N^2-r), by keeping at most r of these singular values. An operator-norm error epsilon|beta| on an N^2-dimensional logical Hilbert space has Hilbert-Schmidt error at most N epsilon|beta|. Therefore

\[
 \boxed{\quad
 (1-\epsilon^2)N^2\leq d^{2\ell}+1,
 \quad}
 \tag{31}
\]

with ell given by (30). If the actual K2 is already within epsilon_0|beta| of a swap plus scalar, the same proof uses epsilon_0 plus the locality-approximation error in place of epsilon, provided their sum is below one.

[proved] Equations (29)-(31) give a coarse but rigorous finite-precision tradeoff under (28). At fixed d, local-norm-to-gap ratio, and visibility strength eta_vis bounded away from zero, arbitrarily large robust swaps are impossible. Both inverse gap and inverse visibility enter. Equation (27) already forces eta_vis to decrease with N for a fixed port, so one must not silently treat a large-N gadget's swap coefficient as order 1/g. The bound does not prohibit many tiny residues, exact additional poles, or arbitrarily many spectator logical states. Those exclusions are necessary in view of (24)-(27).

[proved] The same locality proof applies directly to a one-cell boundary block R_ab(z) for real z>=0. Replace the two copies of v by F_a and F_b, include the phase exp(i z t) in the regulated integral, and replace nu^2 by ||F_a|| ||F_b|| in (29). The regularization error is at most ||F_a|| ||F_b|| omega/(g+z)^2. Thus, at fixed precision, every block throughout a fixed positive z range is approximated by compression of an operator in the same ell-site boundary operator space, whose dimension is at most d^(2ell). For the full block-matrix norm one may allocate an error of epsilon/m to each block. This is a finite-precision bound on the space of linear boundary responses, not on the exact number of poles or the algebra generated by arbitrarily many products of responses.

[proved] A related area-law inference is conditional on different hypotheses. If an actual one-dimensional ground state is unique, uniformly gapped, and subject to an entropy bound S<=C(d,j/g), then a cut carrying an exact maximally entangled logical N-level pair must satisfy log N<=C(d,j/g). The latter step follows because local encoding isometries preserve the Schmidt coefficients of that pair. This is not a bound on arbitrary degenerate cell codes, nor on a gapless logical chain. No area-law theorem or stability of a specific gapless target is proved here, so such an inference cannot finish the embedding argument.

[proved] The abstract port-idle interface therefore removes the boundary dilution obstruction, with exact first-order cancellation and controlled full resolvents. It supplies sufficient syndrome-sector universality and a general initialized dilation, while identifying the false stronger claims. The remaining directions are (i) an explicit static-port local computation with a full-sector gap proof, (ii) a length-selection calculation including inter-cell energy shifts, and (iii) a uniform effective-interaction argument at the precision and symmetry required by the intended thermodynamic theory.
