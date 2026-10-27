# C8 result

2026-10-27. Closed-book work by Codex, GPT-6. No network access or subagents. All project changes are confined to this folder.

## Verdict

Resumed-run update. The whole-class verdict is still nonconclusive. Sections 10-14 contain the new work and supersede the earlier next-step ranking.

The exact product-code route is now excluded for the parity-dressed hopping target when its endpoints are disjoint. A quantitative leakage bound and a gap/length tradeoff explain why the perturbative route remains different. The full flat-parent second-order calculation has an attractive singlet mass-channel coefficient, including a contribution from adjacent bonds. An explicit parity-preserving two-energy parent reverses that coefficient. Its complete 324-dimensional example and checks are provided, but no finite-coupling CFT is claimed.

The all-ring issue is also sharper. Conserved markers force at least L ground states on large prime rings under stated energy-density conditions. An explicit nearest-neighbor marker hop at the crossing of two preferred cell lengths instead produces an XX model of cell lengths at first order. Projected graph calculations distinguish a finite-density length fluid from one quadratic forced defect. A payload code compatible with both lengths has not been constructed.


Not conclusive under the definition in the task. I have neither an explicit microscopic counterexample, a bound for the full class, nor a construction with unbounded central charge at fixed microscopic dimension.

The main new results are finite-dimensional theorems about the proposed block route.

1. Both boundary sites can simultaneously compress to the required scaled Clifford operators, including the parity dressing that makes the inter-block interaction genuinely quadratic in fermions. A constructive dilation proves this. Thus simultaneous boundary compatibility does not close the first-order window.
2. A single one-way gate sweep on an otherwise unrestricted register cannot implement that code. One end acts on only two input sites, regardless of block length.
3. G1's second-order exchange construction extends to a translation-invariant chain of identical blocks with distinct left and right ports. There are no three-block terms at this order. An explicit 48-dimensional block with two qubit ports produces an antiferromagnetic logical qutrit swap. This is an extension of G1, not a new claim to its unrestricted-block result.
4. Exact large microscopic symmetry cannot protect the proposed exchange through a small port. If a simple symmetry algebra has dimension greater than q^2-1, where q is the port dimension, it acts trivially on the interaction algebra. For a single irreducible ground code this forbids every effective interaction, at all orders. In particular this applies to an SU(N) fundamental code with N>q.
5. Generic incommensurate free-fermion targets provide a possible alternative to exact large symmetry. Their low-dimension perturbations can be classified. The operator classification is conditional on the stated bosonization input; the inference of a stable microscopic phase remains an RG argument, not a theorem proved here.

The outstanding problem is now particularly concrete: realize the joint boundary code, or the joint boundary spectral weights, with a repeated microscopic nearest-neighbor bond. Neither an arbitrary parent projector nor an arbitrary matrix on a block meets that requirement.

The numerical checks below concern the abstract codes, the unrestricted-block gadget and the exact logical free-fermion target. They are not counterexample numerics.

## 1. Simultaneous boundary Clifford compression

Status: proved by explicit finite-dimensional construction.

Let q>=2 and let m=2n<=q^2-1. On D=2^n dimensions choose Hermitian Clifford matrices

\[
\{\Gamma_a,\Gamma_b\}=2\delta_{ab}I,\qquad
P=Z^{\otimes n},\qquad \{P,\Gamma_a\}=0.
\]

One concrete convention is

\[
\Gamma_{2j-1}=Z^{\otimes(j-1)}\otimes X\otimes I^{\otimes(n-j)},\qquad
\Gamma_{2j}=Z^{\otimes(j-1)}\otimes Y\otimes I^{\otimes(n-j)}.
\]

Put R_a=i Gamma_a P. These matrices are Hermitian and obey the same Clifford relations. For example,

\[
R_a^\dagger=-iP\Gamma_a=i\Gamma_aP,
\qquad R_aR_b=\Gamma_a\Gamma_b.
\]

Choose m traceless Hermitian q by q matrices F_a with

\[
\operatorname{tr}(F_aF_b)=q\delta_{ab},\qquad \|F_a\|\leq\sqrt q.
\]

The normalized real and imaginary off-diagonal matrix units and the usual traceless diagonal matrices provide such a basis. The norm bound also follows directly from tr(F_a^2)=q.

Theorem 1. There is an isometry

\[
V:\mathbb C^D\longrightarrow
\mathbb C^q_L\otimes\mathbb C^q_R\otimes\mathbb C^{q^2D}_E
\]

such that, simultaneously for every a,

\[
V^\dagger(F_a\otimes I\otimes I)V=\alpha\Gamma_a,
\qquad
V^\dagger(I\otimes F_a\otimes I)V=\alpha R_a,
\qquad
\alpha=\frac1{4m\sqrt q}.
\tag{1}
\]

The environment dimension is an upper bound; unused dimensions may be padded. Its tensor factors can be reordered to place E between the two boundary sites.

Proof. Define T:M_(q^2) -> M_D by

\[
T(A)=\frac{\operatorname{tr}A}{q^2}I+
\frac\alpha{q^2}\sum_a
\left[
\operatorname{tr}((F_a\otimes I)A)\Gamma_a+
\operatorname{tr}((I\otimes F_a)A)R_a
\right].
\tag{2}
\]

Tracelessness gives T(I)=I. With Choi convention

\[
J_T=\sum_{i,j}|i\rangle\langle j|\otimes T(|i\rangle\langle j|),
\]

we obtain

\[
J_T=\frac1{q^2}\left[I+
\alpha\sum_a\left((F_a^T\otimes I)\otimes\Gamma_a+
(I\otimes F_a^T)\otimes R_a\right)\right].
\tag{3}
\]

The norm of the perturbation inside the brackets is at most

\[
2\alpha\sum_a\|F_a\|\leq2\alpha m\sqrt q=\frac12.
\]

Consequently J_T >= I/(2q^2)>0, so T is completely positive. This invokes only the finite-dimensional Choi criterion, which can also be verified directly by the factorization below. The trace orthogonality of the F_a gives the two images in (1).

To make the dilation explicit, diagonalize J_T=S S^dagger. Regard its row index as (i,u), with i an input index of size q^2 and u a logical index of size D. Define

\[
V_{(i,e),u}=\overline{S_{(i,u),e}}.
\tag{4}
\]

Then

\[
(V^\dagger(A\otimes I)V)_{uv}
=\sum_{i,j,e}S_{(i,u),e}A_{ij}\overline{S_{(j,v),e}}
=T(A)_{uv}.
\]

Putting A=I proves V^dagger V=I. The number of columns e of S is at most q^2D. This proves the dimension bound and the theorem.

The result allows both ends to refer to the same logical Clifford algebra. It does not rely on assigning independent logical degrees of freedom to the two ends. In particular, the commutativity of the two physical boundary algebras does not imply that their compressions commute.

There is a useful elementary bound on how sharp simultaneous compressions can be. Suppose physical Hermitian contractions A and B commute and their code images are alpha Gamma and beta R, with Gamma and R anticommuting Hermitian unitaries. Let Q=I-VV^dagger. Then

\[
[V^\dagger AV,V^\dagger BV]
=V^\dagger BQAV-V^\dagger AQBV.
\]

Moreover

\[
\|QAV\|^2\leq1-\alpha^2,\qquad
\|QBV\|^2\leq1-\beta^2.
\]

Taking norms gives

\[
2|\alpha\beta|\leq2\sqrt{(1-\alpha^2)(1-\beta^2)},
\qquad \alpha^2+\beta^2\leq1.
\tag{5}
\]

Thus unit-strength incompatible copies are impossible, but arbitrarily small nonzero scales are allowed. Equation (1) supplies them constructively. This tradeoff concerns amplitudes, not a bound on the number of logical channels.

### A circuit-friendly version for power-of-two ports

Status: proved. This further construction avoids an arbitrary Choi square root when describing the code, although it uses a larger environment.

If q=2^s, choose the F_a from the nonidentity s-qubit Pauli matrices, so that F_a^2=I. For each boundary label e=L,R and generator a, let O_(L,a)=Gamma_a and O_(R,a)=R_a, and let

\[
P_{e,a,\sigma}=(I+\sigma O_{e,a})/2,\qquad \sigma=\pm1.
\]

Prepare on the selected boundary a purification of

\[
\rho_{a,\sigma}=(I+\sigma F_a)/q,
\]

and on the other boundary a maximally entangled q-dimensional pair, whose boundary marginal is I/q. Denote their combined purification by |phi_(e,a,sigma)>. A coherent isometry is

\[
V_{\rm mix}|\psi\rangle=
\frac1{\sqrt{2m}}\sum_{e,a,\sigma}
|e,a,\sigma\rangle_{\rm flag}
|\phi_{e,a,\sigma}\rangle\otimes P_{e,a,\sigma}|\psi\rangle.
\tag{5a}
\]

Orthogonality of the flags and sum_sigma P_(e,a,sigma)=I prove V_mix^dagger V_mix=I. On the selected boundary,

\[
\operatorname{tr}(F_b\rho_{a,\sigma})=\sigma\delta_{ab},
\]

whereas all nonidentity Pauli expectations on the other boundary vanish. Hence

\[
V_{\rm mix}^\dagger F_b^{(L)}V_{\rm mix}=\Gamma_b/(2m),\qquad
V_{\rm mix}^\dagger F_b^{(R)}V_{\rm mix}=R_b/(2m).
\tag{5b}
\]

This is the same joint map (2), at the larger scale alpha=1/(2m). Complete positivity follows directly from (5a), without a triangle-inequality estimate. Coherent recording of a Pauli measurement, preparation of a Pauli eigenspace purification and their controlled versions are finite quantum circuits. When 2m is a power of two, the uniform branch flag is prepared with Hadamard gates. The branch controls can be compiled into a finite circuit; no efficient or uniform microscopic Hamiltonian implementation is asserted.

With two q-dimensional purifiers, the flag, and the retained logical register, this version uses total dimension at most 4m q^4 D. Thus q=16, m=128, n=64 fits in 23 q-state spectator sites and has alpha=1/256. Adding the four-letter selector would give d=64 and the intended c=64>d-1=63. The smaller Choi dilation needs only 20 spectator sites for these same parameters. Both are abstract codes. The explicit circuit is useful input to a future history-state construction, but does not enforce the code with a repeated nearest-neighbor bond.

### Resource accounting and its limitation

The total block Hilbert space used in Theorem 1 has dimension at most q^4 2^n. It therefore embeds in b q-state sites whenever

\[
b\geq\left\lceil4+\frac{n\log2}{\log q}\right\rceil.
\tag{6}
\]

This is a storage bound, not a local-Hamiltonian construction.

For example q=9, m=80, n=40 fits in b=17 spectator sites. Tensoring with C7's four-state selector would give microscopic d=36, whereas the intended logical target below has c=40>35. Here alpha=1/960 and the compressed bond coefficient is reduced by alpha^2=1/921600. All these are finite numbers. No 36-state bond realizing this example has been constructed.

In particular, C7's unperturbed spectator ground space is the entire q^b-dimensional register. It is not the subspace ran(V). Simply placing the isometry into that register does not enforce the code. Choosing a nonlocal parent projector I-VV^dagger would enforce it, but would violate the microscopic nearest-neighbor requirement. This is the unresolved step in this resource example.

## 2. The parity-dressed target is a genuine free theory

Status: exact mapping, using the standard continuum identification of a linear free-fermion dispersion with a unitary free-fermion CFT. Numerical checks are separate below.

On a chain of D-dimensional logical sites define

\[
H_{\rm log}=J\sum_j\sum_{a=1}^{2n} R_{j,a}\Gamma_{j+1,a}.
\tag{7}
\]

This is Hermitian because operators on different sites commute. Define global Majoranas

\[
\gamma_{j,a}=\left(\prod_{k<j}P_k\right)\Gamma_{j,a}.
\]

Their anticommutators are canonical. In the bulk,

\[
i\gamma_{j,a}\gamma_{j+1,a}
=i\Gamma_{j,a}P_j\Gamma_{j+1,a}=R_{j,a}\Gamma_{j+1,a}.
\tag{8}
\]

Thus (7) is a sum of 2n independent Majorana hopping chains. Pairing adjacent species into complex fermions c_(j,l) gives

\[
i\gamma_{j,2l-1}\gamma_{j+1,2l-1}
+i\gamma_{j,2l}\gamma_{j+1,2l}
=2i(c_{j,l}^\dagger c_{j+1,l}-c_{j+1,l}^\dagger c_{j,l}).
\]

Up to reversing momentum, each dispersion is 4J sin k, with two simple Fermi points and velocity 4|J|. There are n complex massless fermions, hence

\[
c=n=m/2.
\tag{9}
\]

The periodic spin chain imposes a common parity-dependent fermion boundary condition. It does not eliminate a bulk oscillator species. This finite parity projection and the associated spin structures change sectors, not the local stress-tensor coefficient.

For J>0 and M logical sites divisible by four, the even-parity antiperiodic sector contains the filled negative-energy sea. Its total occupation nM/2 is even, so the sector is self-consistent. Summing the sine dispersion gives

\[
E_0(M)=-4Jn\csc(\pi/M)
=-\frac{4Jn}{\pi}M-\frac{2\pi Jn}{3M}+O(M^{-3}).
\tag{10}
\]

The periodic-sector minimum is -4Jn cot(pi/M), which is higher. There are zero modes in that sector allowing the required odd parity. Their energy difference is

\[
4Jn\tan(\pi/(2M)),
\]

and multiplication by M/(2 pi v) tends to n/4. This agrees with the total parity-twist scaling dimension. Ordinary particle-hole excitations near each Fermi point have gaps 2 pi v/M times integer oscillator levels, separately for every flavor. Thus this target has finite velocities, conformal towers and z=1, rather than a quadratic clock mode.

Compressing the boundary interaction

\[
h_{\rm port}=\sum_a F_a^{(R)}\otimes F_a^{(L)}
\]

using Theorem 1 gives alpha^2 times the bond in (7), with J=alpha^2 for unit microscopic coupling. Its connected operator Schmidt rank is m. The undressed expression sum Gamma_a tensor Gamma_a would in general be a different, interacting spinor chain, so the right-boundary parity factor cannot be omitted.

This target itself has on-site dimension D=2^n and does not violate c<=D-1. A violation would require its realization inside microscopic sites of the smaller dimension q, including the selector overhead and a local code Hamiltonian.

For other logical ring lengths, the free problem still has flux/parity sectors and linear Fermi points. This does not solve the separate issue of microscopic rings whose length is not divisible by the selector period. Such rings force selector defects, as C7 already explains.

## 3. Obstruction for a one-pass history head

Status: proved for the specified circuit family, not for general computation clocks.

Consider ell q-state data sites with an unrestricted input Hilbert space. Let U_j act on data sites j and j+1. Define W_0=I and W_t=U_t...U_1 for 1<=t<=ell-1. A legal propagation Hamiltonian whose ground space is the history code

\[
V_{\rm hist}|\psi\rangle
=\sum_{t=0}^{\ell-1}\sqrt{p_t}|t\rangle W_t|\psi\rangle,
\qquad \sum_t p_t=1,
\tag{11}
\]

has the following limitation for a data-only operator O_1 on its left endpoint.

Theorem 2. Its compressed image is supported on only the first two input data sites,

\[
V_{\rm hist}^\dagger O_1 V_{\rm hist}
=p_0 O_1+(1-p_0)U_1^\dagger O_1 U_1.
\tag{12}
\]

Proof. Orthogonality of the clock states removes unequal t terms. For t>=1,

\[
W_t^\dagger O_1W_t
=U_1^\dagger\cdots U_t^\dagger O_1U_t\cdots U_1.
\]

Every U_j with j>=2 commutes with O_1. Remove them starting from the center of this expression. The result is U_1^dagger O_1 U_1, independent of t. Summation proves (12).

Suppose m=2n nonzero scaled Clifford generators were among these images. Dividing out their scales gives a representation of the complex Clifford algebra Cl_(2n)=M_(2^n) inside M_(q^2). Every representation of this full matrix algebra is a direct sum of its 2^n-dimensional defining representation. Therefore

\[
2^n\mid q^2,\qquad n\leq2\nu_2(q)\leq2\log_2q.
\tag{13}
\]

Here nu_2(q) is the exponent of two in q. The factors on the remaining input sites are identities, so they do not enlarge the representation on the first two sites. This proves the bound. An odd number of Clifford generators would allow one additional generator; the target of Section 2 uses an even number.

A long cell and gates entangling much of the register do not change this result. It is the ordering of the one-pass gates that limits the endpoint. Reversing the sweep moves the limitation to the other end.

The theorem assumes the full input register is a ground code and the boundary observable is data-only, diagonal in the history clock. It does not cover arbitrary clock-changing boundary operators, additional input-code constraints, multiple sweeps, or a different computation graph. In C7's frozen-marker setting the spectator boundary observables have the required clock-diagonal form. This rules out a straightforward one-pass spectator-gate extension; it is not a general no-go for Theorem 1's code.

## 4. Identical blocks can realize G1's second-order exchange

Status: proved. The unrestricted-block swap construction itself is due to the existing G1 contribution. The new points here are simultaneous left/right realization, identical blocks, and the absence of unwanted second-order multi-block terms.

G1 is at `../../contrib/unbounded_c/G1/G1_result.md`. It already disproves an unrestricted claim that diagonal terms must precede noncommuting exchange. The current parent README incorporates that result. The earlier universal hierarchy objection cannot be used for unrestricted blocks.

### 4.1 Joint spectral-weight realization

Let G:C^N -> H be a ground isometry. Suppose there are two families of orthogonal excitation isometries

\[
S_{R,a}=F_{R,a}G,\qquad S_{L,a}=F_{L,a}G,
\]

orthogonal to G and to one another, including between the left and right families. Write S for their concatenation, with s columns. Given positive matrices M_e on C^s with sum M_e=I, and positive energies E_e, the following realizes all their boundary spectral weights in one block.

Let Q dimension be dim(H)-N, at least s times the number of energies. Form W by vertically stacking sqrt(M_e) and padding zeros. Then W^dagger W=I. Complete W to a unitary U on Q. Choose a physical orthonormal frame F_Q for G^perp whose first s columns are S, and put

\[
H_0=F_QU^\dagger\operatorname{diag}(E_1I_s,E_2I_s,\ldots,E_{\rm pad}I)UF_Q^\dagger.
\tag{14}
\]

It has kernel ran(G), a positive gap, and

\[
S^\dagger\Pi_eS=M_e.
\tag{15}
\]

To verify the last identity, the excited eigenvector matrix is F_QU^dagger, and its overlap with S is W^dagger. Selecting the rows for energy E_e produces sqrt(M_e)^dagger sqrt(M_e)=M_e. This is the same spectral-weight dilation as G1, applied jointly to two orthogonal port syndromes.

In particular, independently prescribed right and left weights can be combined as

\[
M_e=M_{R,e}\oplus M_{L,e}.
\tag{16}
\]

The block then need not be reflection invariant, but it is the same block at every site of the effective chain.

### 4.2 An explicit N=3 example

Use a qubit right port, a qubit left port, one qubit Bell partner for each port, and a three-dimensional logical register. Thus dim(H)=48. In physical order R,L,Raux,Laux,logical, set

\[
G|i\rangle=\frac12\sum_{r,l=0}^1|r,l,r,l,i\rangle.
\tag{17}
\]

Let F_(R,a) and F_(L,a) be the three Pauli matrices on the respective ports. Bell-state expectations give

\[
G^\dagger F_{R,a}G=G^\dagger F_{L,a}G=0,
\]

\[
G^\dagger F_{R,a}F_{R,b}G=\delta_{ab}I,
\qquad G^\dagger F_{R,a}F_{L,b}G=0,
\tag{18}
\]

and similarly on the left. Thus the six syndrome isometries used in (14) are orthogonal.

Let P_3 be the qutrit swap, Phi_3 the projector onto 3^(-1/2) sum_i |ii>, and define

\[
F=P_3-I/3,\qquad K=3\Phi_3-I/3,\qquad \eta=9/64.
\]

At energies E_1=1 and E_2=3 choose

\[
M_{R,1}=I/2+\eta F,\quad M_{R,2}=I/2-\eta F,
\]

\[
M_{L,1}=I/2-\eta K,\quad M_{L,2}=I/2+\eta K.
\tag{19}
\]

These matrices act on syndrome C^3 tensor logical C^3. Their smallest eigenvalue is 1/8. Combine them by (16), which gives two weights of size 18. The excited physical space has dimension 45, so (14) fits with nine padding states at energy 3. The complete local spectrum is 0 with multiplicity 3, 1 with multiplicity 18, and 3 with multiplicity 27.

For a chain of these identical blocks, take

\[
H(\lambda)=\sum_j H_{0,j}
+\lambda\sum_j\sum_{a=1}^3 F_{R,a,j}F_{L,a,j+1}.
\tag{20}
\]

The first-order effective Hamiltonian vanishes by (18). The second-order contribution from one bond is

\[
K^{(2)}=-\sum_{e,f}\frac1{E_e+E_f}
\sum_{a,b=1}^3 M_{R,e}^{ab}\otimes M_{L,f}^{ab}.
\tag{21}
\]

The superscripts are syndrome indices; each entry is a logical matrix. For clarity, the two port insertions on one block give

\[
G^\dagger F_a\Pi_eF_bG=S_a^\dagger\Pi_eS_b=M_e^{ab}.
\]

There is no complex conjugation on the second factor of (21).

Choose an orthonormal traceless Hermitian qutrit basis T_mu. Then

\[
F=\sum_\mu T_\mu\otimes T_\mu,
\qquad K=\sum_\mu T_\mu^T\otimes T_\mu.
\]

The constant part of (21) is

\[
-\frac34\left(\frac12+\frac14+\frac14+\frac16\right)I=-\frac78I.
\]

The terms linear in eta vanish by tracelessness in syndrome space. For the quadratic term, contraction of a T_mu with a transposed T_nu gives delta_(mu,nu), and the energy coefficient is

\[
\eta^2\left(\frac12-\frac14-\frac14+\frac16\right)=\frac{\eta^2}{6}=\frac{27}{8192}.
\]

Therefore

\[
K^{(2)}=-\frac78I+\frac{27}{8192}(P_3-I/3)
=-\frac{7177}{8192}I+\frac{27}{8192}P_3.
\tag{22}
\]

The sign is antiferromagnetic. There is no nonscalar one-block term.

Two different bonds cannot contribute at this order. After the first bond acts, both of its blocks are excited, since each P F_a P vanishes. The resolvent of sum H_0 preserves the ground/excited decomposition of each block. A distinct second bond leaves at least one of the first bond's outer blocks untouched and excited, so final projection onto all ground spaces annihilates the result. This includes neighboring bonds. For ordinary rings with at least three blocks, the second-order chain Hamiltonian is consequently the sum of (22), not a Hamiltonian with additional three-block interactions.

The formula concerns a perturbative coefficient at each finite chain length. It is not an assertion of a global, volume-independent isolated band containing all 3^M logical states at finite lambda. Such a band need not stay isolated as M grows. A local effective expansion or a stability theorem would be needed to use it in the thermodynamic limit.

### 4.3 Arbitrary rank at second order through a product bond

The identical-block extension also applies to G1's rank-one construction, and more generally to an arbitrary connected logical bond K on C^N tensor C^N.

Here is a self-contained algebraic version. Take two qubit ports with

\[
G|i\rangle=|00\rangle|i\rangle,
\qquad S_R=X_RG,\qquad S_L=X_LG,
\qquad V_{j,j+1}=X_{R,j}X_{L,j+1}.
\]

Padding the auxiliary space as necessary is allowed in this subsection. Let m=N^2-1, and expand

\[
K=\sum_{\mu,\nu=1}^m J_{\mu\nu}T_\mu\otimes T_\nu
\]

in an orthonormal traceless Hermitian basis. Hermiticity makes J real. Choose L=m+2 distinct positive energies and define

\[
C_{ef}=(E_e+E_f)^{-1},\qquad w_e=1/L.
\]

The Cauchy matrix C is positive definite since

\[
u^TCu=\int_0^\infty\left(\sum_eu_e e^{-E_et}\right)^2dt>0
\]

for any nonzero real u. Linear independence of distinct exponentials proves strictness. The vectors 1 and Cw are independent because (Cw)_e strictly decreases with E_e. Their joint orthogonal complement has dimension m. Whitening C on this complement gives an L by m real matrix X satisfying

\[
1^TX=0,\qquad w^TCX=0,\qquad X^TCX=I.
\]

Set

\[
A_e=\sum_\mu X_{e\mu}T_\mu,
\qquad B_e=-\sum_{\mu,\nu}X_{e\mu}J_{\mu\nu}T_\nu,
\]

\[
M_{R,e}=w_eI+\epsilon A_e,\qquad
M_{L,e}=w_eI+\epsilon B_e.
\]

For sufficiently small positive epsilon every matrix is positive, and each family sums to I. Use the joint dilation (14)-(16). Applying (21) with just one syndrome channel, the linear terms vanish by w^TCX=0, and the quadratic term is epsilon^2 K. Hence

\[
H^{(1)}=0,\qquad H^{(2)}=-w^TCw\,I+\epsilon^2 K.
\tag{23}
\]

This includes a positive swap modulo its scalar part, with connected rank N^2-1, through a single product bond. Choosing small epsilon is consistent with the necessary negative semidefiniteness of the full second-order operator because of its scalar offset.

The port dimensions and the interaction rank stay fixed; the block dimension, spectral data and precision grow with N. The construction is not a fixed-microscopic-d nearest-neighbor chain because (14) is an arbitrary whole-block Hamiltonian. It also does not establish the IR of the finite-lambda chain. In particular, a large second-order rank is not by itself a counterexample to c<=r(h).

### 4.4 The actual identical-port sign obstruction

For comparison, suppose the two ends use the same single Hermitian port operator F and the same block Hamiltonian H_0. Assume G^dagger F G=0. Define

\[
R(t)=G^\dagger F e^{-tH_0}F G.
\]

It is Hermitian and positive. The bond F tensor F produces

\[
K^{(2)}=-\int_0^\infty R(t)\otimes R(t)\,dt.
\tag{24}
\]

Expand R(t)=r_0(t)I+sum_mu r_mu(t)T_mu. The real connected coefficient matrix is

\[
J_{\mu\nu}=-\int_0^\infty r_\mu(t)r_\nu(t)\,dt,
\]

which is negative semidefinite. It cannot equal the positive definite coefficient matrix of an antiferromagnetic swap. This obstruction assumes the same port channel on both sides. The independently engineered left and right ports in (23) replace the Gram matrix by a cross-Gram matrix and remove the sign restriction.

## 5. An all-orders obstruction to exact symmetry protection

Status: proved. No assumptions about nearest-neighbor structure inside the blocks are used.

Theorem 3. Consider identical finite-dimensional blocks with ground projector P, and a Hermitian boundary interaction v acting through a q-dimensional right port and a q-dimensional left port. Suppose a compact connected group with simple Lie algebra g acts by an on-block unitary U(g) such that

(i) H_0 is invariant;

(ii) the connected bond v is invariant under U(g) tensor U(g);

(iii) dim(g)>q^2-1;

(iv) the ground space carries a single irreducible representation of the group, with no multiplicity.

Then each Schmidt operator in the boundary interaction commutes with the on-block group action. The full chain has an independent copy of that symmetry on each block. Every effective Hamiltonian on the product ground code is scalar, to all perturbative orders and also for symmetry-preserving exact low-energy reductions whenever they are defined.

Here connected means that one-block and scalar terms have been removed. Those can be absorbed into H_0; the assumption requires H_0 so obtained to be invariant.

Proof. Use a minimal operator Schmidt decomposition of v. The right Schmidt span is the range of the linear map obtained by contracting v against arbitrary operators on the second block. Invariance of v under U tensor U implies that conjugation by U preserves this span. The same holds for the left span. This statement is independent of degeneracies among Schmidt coefficients.

Let A_R be the unital star algebra generated by the right Schmidt span. Conjugation by U therefore acts on A_R by star automorphisms. Since the interaction acts only on the right port, A_R is represented inside M_q tensor I_aux. As an abstract finite-dimensional C-star algebra it has a decomposition

\[
\mathcal A_R\simeq\bigoplus_\alpha M_{s_\alpha},
\qquad \sum_\alpha s_\alpha\leq q.
\]

The inequality follows from the faithful representation on the q-dimensional port, allowing multiplicities of each simple summand. A connected automorphism group cannot permute its finitely many central projections. On each full matrix summand every star derivation is inner. Thus the real Lie algebra of star derivations is

\[
\operatorname{Der}_*(\mathcal A_R)\simeq
\bigoplus_\alpha\mathfrak{pu}(s_\alpha),
\]

and

\[
\dim\operatorname{Der}_*(\mathcal A_R)
=\sum_\alpha(s_\alpha^2-1)\leq q^2-1.
\tag{25}
\]

Differentiating the group action gives a Lie algebra homomorphism

\[
\mathfrak g\longrightarrow\operatorname{Der}_*(\mathcal A_R).
\]

Its kernel is an ideal. Since g is simple, a nonzero homomorphism would be injective. The dimension inequality forbids injectivity, so the homomorphism is zero. Connectedness then makes the group action on A_R trivial. Repeat on the left.

Consequently each boundary Schmidt operator commutes with U(g), not merely the diagonal bond. H_0 also commutes with it. One can therefore rotate any single block independently without changing the full Hamiltonian. Ground projectors, excited resolvents and every term of an effective perturbation expansion commute with these independent rotations. The tensor product of the single irreducible ground codes is irreducible as a representation of the product group. Schur's lemma makes every such effective term scalar. The same argument applies to an equivariant exact reduction. This proves the theorem.

Corollary 3a. An exact SU(N)-protected exchange of fundamental logical spins through a q-dimensional port is impossible for N>q, at every order. Here dim su(N)=N^2-1. This extends the first-order rank obstruction to all orders under the exact symmetry hypothesis.

Corollary 3b. If a spinor code for 2n Majoranas is to be protected by an exact simple Spin(2n) symmetry of the block and bond, a necessary condition for nontrivial dynamics is

\[
n(2n-1)\leq q^2-1.
\tag{26}
\]

For the large n of interest, so(2n) is simple. This excludes exact microscopic Spin(80) protection of the q=9, n=40 resource example by a large margin. The Clifford representation on 2^n dimensions splits into parity spinors under Spin(2n), so one must apply the scalar conclusion separately to a single irreducible ground code. Even without hypothesis (iv), the preceding proof already forces all port operators to be symmetry invariant, preventing them from transforming as the desired vector Gamma_a. The dimension obstruction to the proposed exact protection therefore still applies to that full Clifford code.

What the theorem does not say is as important as its hypotheses. It does not forbid emergent SU(N) or Spin(2n), approximate symmetry, reducible codes with neutral multiplicity spaces, or a critical theory without a large internal symmetry. G1's gadgets are not invariant under the needed exact diagonal SU(N), so there is no contradiction. The theorem is not a bound on the central charge of the full class.

## 6. A possible stable target without exact large symmetry

Status: the free lattice target is explicit and exactly solvable. The perturbation classification uses standard Abelian bosonization as an assumed continuum input. Stability against all sufficiently small microscopic perturbations is an RG inference, not a lattice theorem proved here.

The highly symmetric hopping point in (7) should not be assumed stable. One can instead use n complex flavors with distinct incommensurate Fermi points. In the fermion description take

\[
H_*=\sum_{j,a}\left[
 e^{-i\theta_a}c_{j,a}^\dagger c_{j+1,a}
+e^{i\theta_a}c_{j+1,a}^\dagger c_{j,a}
-2\cos k_a\,c_{j,a}^\dagger c_{j,a}\right].
\tag{27}
\]

The dispersion is 2 cos(k-theta_a)-2 cos k_a, with Fermi points theta_a +/- k_a and velocities of magnitude 2 sin k_a. For the particular n=40 resource example, choose distinct primes p_a and q_a among the first 80 primes and set

\[
\theta_a=\sqrt{p_a}/100,
\qquad k_a=\pi/2+\sqrt{q_a}/100.
\tag{28}
\]

These k_a lie strictly between zero and pi. The target has c=n and nonzero velocities. In logical spin variables it is a nearest-neighbor Hermitian chain on dimension 2^n, with the total-parity Jordan-Wigner convention of Section 2. The intersite terms are real linear combinations of R_a tensor Gamma_b, and the chemical potentials are on-site Clifford bilinears. Thus it uses at most 2n boundary Schmidt channels and is algebraically compatible with Theorem 1. Enforcing the required internal logical chemical potentials by a repeated microscopic bond is part of the unresolved local realization problem.

For these Fermi points there is no nontrivial integer momentum resonance,

\[
\sum_a\left[u_a(\theta_a+k_a)+v_a(\theta_a-k_a)\right]\in2\pi\mathbb Z
\quad\Longrightarrow\quad u_a=v_a=0\ \text{for every }a.
\tag{29}
\]

Indeed, the left side is an algebraic combination of distinct square roots of primes plus (sum_a(u_a-v_a))pi/2. Transcendence of pi first forces the coefficient of pi after subtracting 2pi times the integer to vanish. Rational linear independence of distinct prime square roots then gives u_a+v_a=0 and u_a-v_a=0 for each a. Equation (29) follows.

The standard free-boson/fermion operator input is the following.

* A local parity-even lattice operator becomes a sum of integer-charge fermion vertex operators and their derivatives. At the free point a vertex with chiral integer charges Q_(a,R),Q_(a,L) has dimension one half of the sum of their squares, plus a nonnegative descendant level. Its lattice momentum is the corresponding charge-weighted sum of Fermi points.
* A local parity-odd spin operator changes the common Jordan-Wigner spin structure. Its lightest continuum fields twist all n complex flavors, with total dimension n/4. Equivalently, the chiral Ramond weights are n/8 on each side. Derivatives or higher charges do not lower that dimension.
* Momentum-zero current products are the exactly marginal forward-scattering interactions of a multicomponent Luttinger liquid, while current terms shift the Fermi data.

Using these inputs, (29) eliminates every nonzero-charge parity-even vertex at a translation-invariant coupling. The remaining dimension-one operators are currents. At dimension two one has current products, stress-tensor terms and derivatives. For n>=9 all parity-odd twist fields have dimension greater than two at the free point.

At the perturbative level this removes relevant operators without imposing an exact internal symmetry. Current terms move the Fermi points and current products move the Gaussian parameters. Only finitely many charge vectors have dimension at most two in a neighborhood of the free point, so excluding their resonance hyperplanes is an open condition; one need not demand an open set of full Diophantine independence at every charge. Twist dimensions also retain a positive margin above two for sufficiently small changes when n>=9.

This is a more specific stability proposal than assuming that arbitrary errors around an SU(N) point are harmless. It does not prove that a microscopic encoded chain realizes this phase. It still needs control of the effective local errors and a stability theorem applicable to those errors. No such theorem is supplied here, and no finite-lambda encoded chain has been measured.

### A conditional challenge to c<=r(h)

There is a useful distinction between the c-d conjecture and the stronger bond-rank proposal. For the latter, the whole unrestricted block may itself be declared a physical site. This produces a large but finite microscopic d, and removes the need to implement its H_0 by smaller sites. It does not remove the infrared stability problem.

The rank-one construction (23) can be made more economical for a target bond of connected Schmidt rank r. Write that bond as sum_(a=1)^r A_a tensor B_a, with traceless Hermitian factors, absorbing real singular values and signs into the factors. Use r whitened columns, L=r+2 energies, and right/left spectral weights built from A_a and -B_a. The same proof gives the target at second order. A total physical dimension

\[
d=4\left\lceil\frac{(2L+1)N}{4}\right\rceil
\tag{30}
\]

fits the two qubit ports, the N-dimensional ground code and all prescribed excited weights. The argument uses no basis directions outside the target's Schmidt span.

For the n=9 target (27), N=512 and r=18, so L=20 and d=20992 suffice. Choose, for example, E_e=e, w_e=1/20, form the whitening matrix on the complement of 1 and Cw, and take any epsilon less than the explicit positivity threshold

\[
\epsilon<\frac{1/20}{\max_e\{\|\sum_aX_{ea}A_a\|,\|\sum_aX_{ea}B_a\|\}}.
\]

This gives a finite, algorithmically specified on-site H_0 of the form (14). Add the target chemical-potential operator b on the code by the microscopic on-site term lambda^2 epsilon^2 G b G^dagger. The microscopic bond is

\[
h_\lambda=\frac12(H_0+\lambda^2\epsilon^2GbG^\dagger)\otimes I
+\frac12 I\otimes(H_0+\lambda^2\epsilon^2GbG^\dagger)
+\lambda X_R\otimes X_L.
\tag{31}
\]

Its connected operator Schmidt rank is exactly one for lambda nonzero. The repeated site and bond are identical throughout the chain. At second order its logical Hamiltonian is a scalar plus lambda^2 epsilon^2 times the c=9 target, including the chemical potentials.

The following conditional implication is rigorous as an implication, but its two premises are not proved here.

(A) For this product-block construction there is a thermodynamic low-energy reduction with local errors o(lambda^2 epsilon^2), at fixed epsilon and fixed block dimension, preserving the gap to auxiliary modes. This must be a local statement, not a global operator-norm approximation whose error grows uncontrolled with length.

(B) The explicit n=9 incommensurate logical target lies in an open c=9 phase under the quasi-local, translation-invariant perturbations delivered by (A), including the parity-odd spin perturbations.

If (A) and (B) hold, sufficiently small fixed nonzero lambda in (31) is a counterexample to c<=r(h), with c=9 and r(h)=1. To prove the implication, divide the effective Hamiltonian by its positive overall scale lambda^2 epsilon^2. Premise (A) places it arbitrarily close in the required local norm to H_*, and (B) keeps c=9. The auxiliary sector remains gapped by (A). The rank computation of the actual microscopic bond is independent of either premise.

This supplies a concrete conditional route against the rank conjecture, beyond observing that second-order rank can grow. Neither premise is replaced by the finite-volume checks of Section 7. In particular, no value of lambda has been certified to have a c=9 infrared spectrum in the actual d=20992 chain. It is not a verified counterexample to c<=r(h), and even a successful verification would leave c<=d-1 intact for this example.

## 7. Numerical checks and reproduction

All code uses Python, numpy and scipy. Run from this folder,

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 c8_checks.py --ed
```

It writes `checks.json` and `two_port_N3.npz`. The latter contains the complete 48 by 48 H0, the 48 by 3 ground isometry, the six physical port operators and the physical second-order matrix. The construction is specified independently by (14), (17)-(19), with the deterministic QR completion implemented in the script. Different admissible completions have the same second-order matrix but can differ at higher order.

The run takes only a few seconds on one machine. No cluster calculation was needed. The numerical checks are consistency checks on the proofs, not their justification.

### Joint compression

| q | m | Logical dimension | Minimum Choi eigenvalue | Largest isometry/compression error |
| :-- | :-- | :-- | :-- | :-- |
| 2 | 2 | 2 | 0.1875 | 2.3e-16 |
| 3 | 4 | 4 | 0.0870548499 | 3.2e-16 |
| 3 | 8 | 16 | 0.0932887158 | 3.7e-15 |

Both boundary maps were checked for every generator. The Clifford identities have zero residual in these explicit Pauli representations. The Choi lower bounds from (3), respectively 1/8 and 1/18, are respected.

The stronger Pauli-mixture version was also checked at (q,m)=(2,2) and (4,8). Its Choi minima were 0.0732233047 and 0.0263202057, with compression errors below 5.0e-16 and 6.1e-15. These use alpha=1/4 and 1/16, respectively.

### One-pass obstruction

For q=2 and five data sites, independently sampled unitary two-site gates were applied in monotone order. Every endpoint conjugate after the first gate agreed with U_1^dagger O_1 U_1 to less than 4.5e-15. The script uses a fixed random seed.

### Identical-block second order

The complete physical 48-dimensional construction gives

\[
\|H^{(2)}+7177I/8192-27P_3/8192\|_F<8.4\times10^{-16}.
\]

The local first-order compressions are zero, the ground-kernel residual is below 4.6e-16 and the orthonormal-frame residual is below 1.4e-15. The smallest spectral-weight eigenvalue is 1/8. This is a direct physical transition-amplitude calculation using every excited eigenvector and its denominator, independently of the contracted weight formula (21).

No finite-coupling CFT estimate was made for the 48-state block chain. It would not be a counterexample to c<=d-1 in any case, since its physical site dimension, if the whole block were a site, would be 48.

### Logical free target

Two independent estimators were calibrated on the exact n=1 chain, with n=2 checked directly in the spin representation.

For logical sizes M=40,48,64,80,96,128, the code computes the occupied-momentum correlation matrix and the half-chain entropy from its eigenvalues. Fitting

\[
S(M/2)=s+(c/3)\log M+b/M^2
\]

gives c=1.000000235 for one flavor. Separately, fitting E_0(M)/M to a constant plus M^(-2) and M^(-4), with velocity 4, gives c=0.999999661. Particle-hole oscillator velocities rise from 3.94669 at M=40 to 3.99478 at M=128, consistent with 4. For n flavors the exact energies and interval entropies multiply by n.

Independent sparse diagonalization of the periodic logical spin Hamiltonian gave

| n | M | E0 | Half-chain entropy |
| :-- | :-- | :-- | :-- |
| 1 | 4 | -5.65685424949 | 0.83299106140 |
| 1 | 8 | -10.45250371901 | 1.04373803047 |
| 2 | 4 | -11.31370849898 | 1.66598212280 |
| 2 | 8 | -20.90500743802 | 2.08747606094 |

The energies agree with (10), and the n=2 entropies are twice the n=1 values. Eigensolver residuals are below 1.1e-13 in the reported runs. Small changes in these last digits are possible because the eigensolver start vector is not fixed.

The target calculations validate the parity dressing, normalization, velocity and central-charge counting. They do not validate a microscopic q-state or 4q-state construction.

## 8. Approaches that remain incomplete or fail

The simultaneous-boundary obstruction fails as a general algebraic argument, by Theorem 1. Imposing complete positivity on both ends does not force a logarithmic number of channels when arbitrarily small compression scales are permitted.

The one-pass history-head implementation fails by Theorem 2. A long register by itself does not give both boundaries access to many Clifford directions.

The proposed universal second-order diagonal hierarchy is false for unrestricted blocks, by G1 and Section 4. It may still express a difficulty for a specific local construction, but it needs a proof that uses that construction's locality. A counting argument based only on the number of ports or insertions is insufficient.

Exact large-symmetry protection fails by Theorem 3. Protecting a large SU(N) or Spin(2n) code through a small port is not achieved merely by choosing a complicated on-block representation.

Generic free targets at a symmetric point cannot be declared stable without studying allowed perturbations. Section 6 supplies a possible incommensurate alternative, but its use in a microscopic construction remains conditional.

A computation-clock or universal-Hamiltonian theorem that simulates a finite target spectrum is not yet enough. It must also control the physical endpoint operators or their matrix-valued resolvents, the ground-code splitting, and errors locally in the thermodynamic limit. The present work has not established those properties for a translation-invariant fixed-d implementation. Invoking universality without them would leave the principal step unproved.

Finally, the q=9, b=17 storage count says nothing about incompatible microscopic ring lengths. The frozen C7 selector leaves static forced defects there. Neither this run nor the abstract code construction makes their full low-energy spectrum conformal.

## 9. Run 1 next steps, superseded by Section 14

1. Construct a local block Hamiltonian realizing a specified joint boundary map. The first concrete target is (1), with the generic quadratic target (27) and its on-site chemical potentials. A two-sweep or programmable history head is a candidate because the one-pass theorem does not cover it. The construction must state the same microscopic h on every internal bond and show the code is its low-energy space. A generic Stinespring isometry or whole-block projector does not suffice.
2. For any such code, prove a local effective-Hamiltonian expansion with controlled errors and compare those errors with the scale alpha^2 lambda. Establish its IR using an applicable multicomponent Luttinger stability result, or tune every dangerous coupling and check the CFT numerically with independent estimators. The free target's conformality alone does not prove this step.
3. In parallel to a first-order code construction, attempt local realization of the joint matrix-valued boundary resolvents in Section 4. The algebraic second-order hierarchy has been removed, so locality and stability are the issues. Exact large SU(N) protection is excluded by Theorem 3 and should not be used as a shortcut.
4. Handle every microscopic ring length, or establish the intended admissible-sector definition explicitly. A successful commensurate subsequence with immobile forced defects elsewhere does not settle the stronger all-ring formulation in the task.
5. For an upper bound, seek a locality restriction on these joint maps or boundary resolvents that goes beyond their complete positivity and port rank. Theorem 1 and (23) show why those two conditions alone cannot provide the missing restriction. No new unconditional c<=f(d) bound follows from this run.

For the separate rank conjecture, a useful intermediate target is to prove or disprove premises (A) and (B) for (31). It avoids the fixed-small-d block-realization problem altogether. A resolution could settle c<=r(h) while leaving the original c-d question open.

The c-d conjecture and the existence of any fixed-d bound remain open after these calculations. The run narrows the missing step to a local realization and its infrared control, with explicit algebraic targets and two restricted no-go theorems available for the next attempt.


## 10. Exact product encodings and a locality tradeoff

2026-10-27, resumed C8 run. All results in Sections 1 and 2 are proved finite-dimensional statements. Section 3 is a conditional quantitative bound using the stated Lieb-Robinson estimate. None is a bound on the full class of critical chains.

### 1. The active boundary algebras of an exact product encoding commute

Partition a periodic chain into M>=4 blocks. Each block has two distinct boundary sites of dimension at most q, and all interblock interactions are supported on these sites. Write

    H = sum_j B_j + sum_j v_(j,j+1),
    V = tensor_j V_j,    P_j = V_j V_j^dagger,    P = tensor_j P_j.

The internal block terms B_j may be arbitrary. The codes need not be ground spaces of B_j. Suppose [H,P]=0. Let K_(j,j+1) be the connected part of the compressed bond, obtained by subtracting its normalized partial traces and scalar. Its left and right operator Schmidt spans are unambiguously defined as traceless real vector spaces of Hermitian logical operators.

Theorem. Every element in either Schmidt span has a Hermitian physical lift on the corresponding boundary site that preserves its block code. The two active logical boundary algebras within a block consequently commute. In particular, the connected effective bonds on successive block links commute.

Proof. Fix block j and apply Q_j=1-P_j on its output and P_i on every other output to H V. Invariance says the result is zero. Every term away from block j vanishes because Q_j V_j=0. What remains is the internal term and the two incident bonds. In logical input/output coordinates for the other blocks it has the form

    E_j = Q_j B_j V_j tensor I
          + sum_alpha Q_j A^R_alpha V_j tensor C_(j+1,alpha)
          + sum_beta D_(j-1,beta) tensor Q_j A^L_beta V_j = 0.       (1)

Here the A's are physical boundary operators and C,D are their neighbors' compressed operators. Let phi be any Hermitian linear functional on logical block j+1 satisfying phi(I)=0. Apply phi there and normalized traces on all other logical blocks except j. Equation (1) gives

    Q_j A_phi V_j=0,
    A_phi = (id tensor phi)((I tensor V_(j+1)^dagger)
                           v_(j,j+1)(I tensor V_(j+1))).          (2)

A_phi is a Hermitian operator on the right boundary site of block j. Its Hermiticity also gives V_j^dagger A_phi Q_j=0, so it commutes with P_j. Its compression is the phi-contraction of the compressed bond. Subtracting its scalar part removes any contribution from a one-block term on j+1. All elements of the connected Schmidt span are obtained by these contractions. This proves the right lift statement. Applying the same argument to block j-1 proves the left one.

The physical left and right lifts have disjoint supports and commute. Since both preserve P_j, restriction is multiplicative on them, hence their logical restrictions commute. The algebras generated by these restrictions commute as well. Expanding two successive connected bonds in their Schmidt spans now proves that the bonds commute. QED.

The assumption M>=4 avoids repeated edges and identifications of the two neighboring blocks. The result needs no translation invariance, gap, frustration freeness, or perturbation expansion. It does need exact invariance of the entire tensor product code, rather than only a few low-energy eigenvectors or a code dressed across block boundaries.

Corollary 1. The hopping target in RESULT.md Section 2 cannot be the exact product-code restriction of such a chain with separated endpoint sites, for any nonzero J and any port dimension.

Indeed, that target's bond is

    K = J sum_(a=1)^(2n) R_a tensor Gamma_a,
    R_a=i Gamma_a P_log.

The two active boundary algebras contain R_a and Gamma_a, whereas

    [R_a,Gamma_b]=-2i delta_ab P_log.                            (3)

They do not commute. This is stronger than a rank restriction, but applies only to an exact product encoding. It does not invalidate the abstract joint compression or its first-order perturbative use.

Corollary 2. Suppose an exactly encoded ladder uses m=2n linearly independent Clifford generators at each boundary, with nonzero coefficients in its connected interaction. Then q>=2^n. If its critical point consists of m independent Majorana CFTs, its central charge satisfies c=m/2<=log_2(q).

For the proof, let A be the unital star algebra generated by the Hermitian reducing lifts at a boundary. It is a subalgebra of M_q, so its simple matrix summands have dimensions at most q. Restriction to the invariant code is a star homomorphism, and its image is a quotient of A. The image contains the unital Clifford algebra Cl_(2n)=M_(2^n). Every nonzero finite-dimensional representation of this matrix algebra has dimension a positive multiple of 2^n. At least one simple summand of the image, hence of A, must therefore have size at least 2^n. This implies q>=2^n. The central-charge statement additionally uses the specified free ladder identification.

For this corollary the left and right ladder Majoranas can be different logical generators, so their parity-dressed interaction algebras may commute. Unlike Corollary 1, this is a dimension bound rather than a universal exclusion. Odd numbers of generators require separate Clifford representation bookkeeping and are not asserted in this corollary.

### 2. A quantitative leakage lower bound

Drop invariance and retain the hopping target of Corollary 1 as the connected compressed bond. Define E_j as the one-block leakage operator in (1), with no assumption that it vanishes. Take the dual functionals

    phi_a(X)=tr(Gamma_a X)/D,
    psi_a(X)=tr(R_a X)/D,

on the appropriate neighbor. They have norm one as functionals on the operator norm, and annihilate I. The lifts extracted in (2), after harmless scalar subtraction, satisfy

    V_j^dagger A_a V_j=J R_a,
    V_j^dagger B_a V_j=J Gamma_a,
    ||Q_j A_a V_j||<=||E_j||,
    ||Q_j B_a V_j||<=||E_j||.                                  (4)

The normalized partial traces used to obtain these inequalities are contractions. Complexifying the maps does not alter their operator norms.

For commuting Hermitian physical A,B and an arbitrary code isometry V,

    [V^dagger A V,V^dagger B V]
       = V^dagger B Q A V - V^dagger A Q B V.

Consequently

    ||[V^dagger A V,V^dagger B V]||
       <=2 ||Q A V|| ||Q B V||.                                (5)

Equations (3)-(5) yield

    ||Q_j A_a V_j|| ||Q_j B_a V_j|| >= |J|^2,
    ||E_j|| >= |J|.                                            (6)

For unequal neighboring target coefficients the last bound is instead the square root of the product of their magnitudes. Internal block terms cannot cancel this extracted leakage because their contribution is proportional to the identity on the neighbor.

This is an obstruction to making the undressed product code arbitrarily close to an invariant subspace relative to its hopping coefficient. It is not an obstruction to perturbation theory: a large excitation penalty can make the resulting energy corrections small despite leakage amplitudes of order J or larger.

### 3. Gap, endpoint separation, and compression scale

Assume a block Hamiltonian has ground projector P, ground energy zero, and a gap Delta>0. Suppose its Heisenberg evolution obeys, for normalized endpoint operators A and B separated by r sites,

    ||[A(t),B]|| <= C exp[-mu(r-v|t|)],                         (7)

as well as the trivial bound 2. Here v>0 has units of sites/time, and C,mu are positive constants. This is the only dynamical assumption in the following bound. No unique ground state is required.

For sigma>0 define the Gaussian energy filter

    A_sigma = integral dt exp[-t^2/(2 sigma^2)] A(t)
                         /sqrt(2 pi sigma^2).

Since every ground energy is zero,

    P A_sigma P=PAP,
    ||P A_sigma Q||, ||Q A_sigma P||
       <= exp[-sigma^2 Delta^2/2].                              (8)

The inequality follows by applying the scalar spectral filter exp[-sigma^2 H_Q^2/2] on the excited side of PAQ. Therefore

    ||[PAP,PBP]||
      <= integral dt g_sigma(t)||[A(t),B]||
           +2 exp[-sigma^2 Delta^2/2].                          (9)

Split the integral at T=r/(2v). Inside, (7) bounds it by C exp[-mu r/2]. Outside, the Gaussian probability is erfc[T/(sqrt(2)sigma)]<=exp[-T^2/(2 sigma^2)]. Using the bound 2 there and choosing sigma^2=T/Delta proves

    ||[PAP,PBP]|| <= C exp[-mu r/2]
                       +4 exp[-Delta r/(4v)].                  (10)

If the two compressions are alpha Gamma and beta R with anticommuting Hermitian unitaries Gamma,R, then

    2|alpha beta| <= C exp[-mu r/2]
                        +4 exp[-Delta r/(4v)].                 (11)

For example, whenever r >= (2/mu) log[C/|alpha beta|], it follows that

    Delta <= (4v/r) log[4/|alpha beta|].                        (12)

The derivation controls the constants in terms of the specified Lieb-Robinson estimate. For a microscopic nearest-neighbor Hamiltonian that estimate can be supplied by a standard Lieb-Robinson theorem. No numerical value of v for a proposed microscopic code is claimed here.

A fixed-gap, fixed-velocity local parent therefore cannot maintain nonvanishing incompatible endpoint compressions as its length tends to infinity. This does not close the construction route, because its block gap and compression amplitudes may shrink with block length. In particular the C7 selector gap already shrinks much faster than 1/r, so (12) creates no contradiction with it.


## 11. Second-order stability and an explicit spectral repair

2026-10-27, resumed C8. The finite-dimensional construction and all effective-Hamiltonian coefficients below are proved. Their continuum interpretation uses the stated free-fermion operator expansion and perturbative RG. No finite-coupling microscopic CFT or gap is proved.

### 1. A graded joint code with controlled second moments

Split a q-dimensional port into grading eigenspaces of dimensions q_+,q_-. Choose m=2n<=2 q_+ q_- traceless Hermitian matrices F_a odd under the grading, normalized by tr(F_a F_b)=q delta_ab. Normalized real and imaginary matrix units between the two grading sectors suffice. Let Gamma_a,R_a=i Gamma_a P_log be the logical operators from RESULT.md. Use the joint completely positive map of that report, with only these odd F_a and a sufficiently small alpha, for example alpha=1/(4m sqrt(q)). Its Stinespring isometry V satisfies

    V^dagger F_La V=alpha Gamma_a,
    V^dagger F_Ra V=alpha R_a,
    V^dagger F_La F_Lb V=V^dagger F_Ra F_Rb V=delta_ab I,
    V^dagger F_La F_Rb V=0.                                    (1)

The second line follows because F_a F_b is grading even. Its only surviving term under the channel is its normalized trace. The last line follows directly from the joint channel, which has no cross-port terms.

This code can preserve an exact total parity. In the positive-square-root Choi dilation, write z for the diagonal port grading. The Choi matrix commutes with z_L z_R tensor P_log. The physical block grading

    Omega = (z_L tensor z_R) tensor (z_L tensor z_R tensor P_log)

then obeys Omega V=V P_log, and both physical port families are odd under Omega. Here the second parenthesis acts on the Choi environment. Thus an Omega-invariant block parent and the port bond sum_a F_Ra tensor F_La preserve product_j Omega_j. This construction does not require a large continuous microscopic symmetry.

For example, q=16 admits m=128 odd channels, with intended c=64. It exceeds d-1 for a d=4q selector by one unit, and exceeds d-1 for the d=3q+1 selector by sixteen units. These are resource comparisons, not microscopic Hamiltonians.

### 2. The flat parent, including adjacent bonds

Let P=VV^dagger, Q=1-P, and take an unrestricted block parent H_0=Delta Q. Couple identical blocks by lambda v, where

    v_(j,j+1)=sum_a F_Ra,j F_La,j+1,
    T_j=sum_a R_a,j Gamma_a,j+1.

The first-order target is lambda alpha^2 sum_j T_j. A finite-volume degenerate perturbation expansion on the product code gives

    H_eff^(2)=(lambda^2/Delta) {
       -M m(1/2+alpha^2) I
       +alpha^4 [(3/2)sum_j T_j^2
                          +sum_j {T_j,T_(j+1)}] }.              (2)

The displayed scalar is for a ring of M>=4 blocks. For an open chain it is proportional to the number of bonds, with only existing adjacent pairs included.

Proof for one bond. On two blocks the inverse excitation energy, after removing Delta, is

    R=P tensor Q + Q tensor P +(1/2) Q tensor Q.

Equation (1) gives Pv^2P=mI, while the intermediate projection onto both codes gives alpha^4 T^2. Intermediate states with only one block excited contribute, in total,

    2m alpha^2 I - 2alpha^4 T^2.

Combining the one- and two-excitation denominators, with the perturbative minus sign, yields

    -PvRvP=-m(1/2+alpha^2)I+(3/2)alpha^4 T^2.                   (3)

For adjacent bonds, the two outer blocks must remain in their codes, so only the common block can be excited. Its relevant matrix element is

    V^dagger F_La Q F_Rb V=-alpha^2 Gamma_a R_b.                 (4)

The two outer matrix elements supply alpha^2, the denominator is Delta, and the minus sign in the second-order formula cancels the minus sign in (4). The two orderings therefore give alpha^4 {T_j,T_(j+1)}/Delta. Disjoint bonds cannot return to the product code after creating an excitation, so they give zero. This proves (2).

Keeping only (3) misses the adjacent-bond term and gives the wrong sign in the staggered mass channel.

### 3. Continuum interpretation and its limitation

At the symmetric hopping point, each Majorana lattice field has linear modes at momenta zero and pi. Use the standard free-field expansion

    gamma_(j,a) ~ psi_Ra(x)+(-1)^j psi_La(x).

The alternating part of the bond T_j is a nonzero real multiple of

    (-1)^j Epsilon(x),    Epsilon=sum_a i psi_Ra psi_La.

Therefore a lattice correction

    a sum_j T_j^2 + b sum_j {T_j,T_(j+1)}

has a staggered mass-channel coefficient proportional to a-2b. Contact contractions alter scalar and lower-dimensional terms but do not change the coefficient of the cross-species quartic operator in this expansion.

For (2), a-2b=-alpha^4/2. Thus the coefficient of Epsilon^2 is negative. This is the attractive sign of the SO(m) Gross-Neveu interaction. One can identify the one-loop tendency without assuming a lattice mass theorem. Put epsilon_a=i psi_Ra psi_La and O=sum_(a<b) epsilon_a epsilon_b. The contraction of a common species in the product of two summands gives a positive O/|z|^2 term. For a given surviving pair, there are m-2 choices of contracted species. Hence the self-OPE coefficient is positive for m>2. With Euclidean action S_0+g integral O, integrating out a short-distance shell gives dg/dell=-C(m-2)g^2+O(g^3), C>0. Negative g runs to larger negative magnitude, while positive g runs toward zero at this order.

For m=2 the factor m-2 vanishes and the interaction is the usual Gaussian marginal deformation. For m>2 the flat parent is therefore not justified as a stable realization of the symmetric free target. Standard Gross-Neveu field theory associates the attractive flow with mass generation and spontaneous dimerization. In this report that mass-generation identification is assumed continuum input, and the conclusion about the actual microscopic chain remains heuristic. Higher orders, symmetry breaking, and shifted Fermi points have not been controlled.

### 4. Two-energy spectral weights reverse the coefficient

The negative sign in Section 3 is not forced by the joint code. It can be changed explicitly while preserving (1), the code gap, and total parity.

Let S_s=Q F_s V, with s running through the m left and m right channels. Write

    U_L=[Gamma_1 ... Gamma_m],    U_R=[R_1 ... R_m],
    U=[U_L U_R].

The syndrome Gram matrix is

    G=S^dagger S=I-alpha^2 U^dagger U.                          (5)

Since U U^dagger=2m I, this is positive when 2m alpha^2<1. Prescribe excited spectral weights on the syndrome source as

    M_e = p_e I + [ s_e U_L^dagger U_L   t_e U_L^dagger U_R
                    t_e U_R^dagger U_L   s_e U_R^dagger U_R ].   (6)

Their sum equals G provided

    sum_e p_e=1,    sum_e s_e=sum_e t_e=-alpha^2.                (7)

Positivity is transparent. On the dark subspace the eigenvalue is p_e. On the two bright copies of the D-dimensional logical space the eigenvalues are

    p_e+m(s_e+t_e),    p_e+m(s_e-t_e).                          (8)

Thus it suffices that p_e>=0 and p_e+m s_e>=m|t_e|.

Use two excitation energies and weights

    E_1=Delta,           E_2=3Delta,
    p_1=p_2=1/2,
    s_1=s_2=-alpha^2/2,
    t_1=0,              t_2=-alpha^2.                          (9)

The smallest weight eigenvalue is 1/2-(3/2)m alpha^2, which is strictly positive for the alpha chosen above. Their sum is (5).

For completeness these are spectral weights of an actual finite matrix parent, not merely formal coefficients. Form the positive normalized measure

    Mtilde_e=G^(-1/2) M_e G^(-1/2),    sum_e Mtilde_e=I.

Its square-root dilation W=stack_e sqrt(M_e) G^(-1/2) is an isometry. The physical frame S G^(-1/2) is also an isometry. Extend both to orthonormal frames, map one frame to the other, and assign energies E_e to the corresponding coordinate blocks. Padding states can have energy 3Delta. This is the same finite-dimensional spectral dilation as RESULT.md Section 4, now with a nonidentity syndrome Gram matrix. It leaves ker H_0=ran V and gap Delta. The dilation needs at most 4mD excited dimensions, and the original q^4 D Choi block has enough room for the stated odd-channel family.

Every M_e commutes with I_(2m) tensor P_log because its nonscalar entries are products of two odd logical operators. The dilation and its orthonormal completions can consequently be chosen separately in the two parity sectors. The resulting H_0 commutes with Omega. The supplied Python implementation performs these separate completions.

To calculate the effective coefficients, remove Delta from the energies and define

    C_ef=1/(E_e+E_f),
    p_inv=sum_e p_e/E_e,    s_inv=sum_e s_e/E_e,
    t_inv=sum_e t_e/E_e.

For general weights (6), the same-bond correction is a scalar plus A T_j^2 and the adjacent correction is B {T_j,T_(j+1)}, with

    A=-2alpha^2 s_inv-s^T C s,
    B=-alpha^2 t_inv,                                         (10)

and the scalar per bond is

    -m[p^T C p+2p^T C s+2alpha^2 p_inv].                       (11)

The first term in A comes from single-block virtual excitations. The second comes from two simultaneously excited blocks. Formula B is the cross-port resolvent on their common block. These contributions are obtained by the same insertion of spectral projectors used in (3)-(4).

For (9), in units Delta=1,

    s_inv=-(2/3)alpha^2,    t_inv=-(1/3)alpha^2,
    s^T C s=(7/24)alpha^4.

Hence

    A=(25/24)alpha^4,    B=(1/3)alpha^4,
    A-2B=(3/8)alpha^4>0.                                     (12)

The scalar is -m[7/24+(3/4)alpha^2]. The leading hopping remains lambda alpha^2 sum T_j. The sign of the singlet Gross-Neveu coupling has been reversed. No first-order splitting was added.

A useful more general one-parameter choice is t_1=(k-1/2)alpha^2 and t_2=-(k+1/2)alpha^2, with the same p,s. While (8) remains positive, the mass-channel coefficient is

    (A-2B)/alpha^4 = 4k/3-7/24.

It crosses zero at k=7/32. This illustrates that the adjacent-bond resolvent, not only the single-port spectral weights, controls the perturbative mass-channel sign.

The positive sign in (12) addresses one dangerous continuum coupling. It does not prove a stable lattice CFT. Higher orders can break the accidental continuous symmetry, generate other marginal couplings, and shift Fermi points. Exact total parity is retained, but it alone does not establish infrared stability.

### 5. Reproduction and checked values

Run

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 c8_resume_checks.py

This writes resume_checks.json, spectral_parent_checks.json, and spectral_parent_q3_m4.npz. The NPZ contains the complete 324-dimensional H_0, code isometry, port matrices, grading, spectral weights, and second-order operators. The physical order is L,R,environment; a port matrix acts on its named q-dimensional factor and as identity elsewhere. Reordering the environment between the two ports makes them geometric endpoints, but does not make H_0 nearest neighbor.

For the flat-parent three-block checks, (q,m)=(2,2) and (3,4), the maximum isometry/compression/second-moment residual is 1.2e-15. The full second-order residual is at most 4.9e-14, and the adjacent-bond residual is at most 2.1e-17. These are direct transition-amplitude contractions in the subspace reached from the code by one port insertion. The flat parent preserves this subspace, so the reduction is exact for the coefficient being checked.

For the two-energy q=3,m=4 parent, the dimension is 324, the gap is 1 to numerical precision, and the smallest spectral-weight eigenvalue is 0.4921875. Kernel, parity, spectral-weight, same-bond and adjacent-bond residuals are all below 1.0e-14. The measured coefficient ratios are A/alpha^4=25/24, B/alpha^4=1/3, and (A-2B)/alpha^4=3/8.

One implementation check initially failed because a spectral factor of the Choi matrix was used while the parity formula assumed its positive square root. Replacing that factor by the positive square root makes the environment grading match the written construction. The corrected full run passes. This was a dilation-basis mismatch, not an adjustment of the desired spectral weights.

No central charge was estimated for the finite-lambda block chain. These checks validate the algebra and the sign-changing construction, not a counterexample to c-d or a proof of microscopic criticality.


## 12. Frozen selectors and incompatible rings

2026-10-27, resumed C8. The frozen-marker theorem is proved under the explicit energy separation hypotheses below. The mobile-defect statement is proved only for the stated finite-band model. Neither is a theorem excluding every possible spontaneously blocked critical chain.

### 1. Frozen marker patterns force unbounded ground-state multiplicity

Let H_L be one-site translation invariant. Suppose every site has a binary marker projector M_i and

    [H_L,M_i]=0    for every i.                                (1)

Thus marker patterns label orthogonal invariant sectors. The marker can coexist with arbitrary quantum spectator degrees of freedom. Assume there are constants e in R, eta>0, C_0,C_1, independent of L, with the following properties for all sufficiently large L.

(a) Each constant marker pattern, all zeros or all ones, has energy at least (e+eta)L-C_0 in every state.

(b) Some state with a nonconstant marker pattern has energy at most eL+C_1.

Here e need not be positive. The notation eta>0 is the energy-density separation, not a conformal gap. These hypotheses formalize a stable marker crystal and a bounded cost for fitting it on an incompatible ring.

Theorem. For every sufficiently large prime ring length L, the ground-state multiplicity is at least L.

Proof. For L>(C_0+C_1)/eta, (b) lies below both constant-pattern sectors in (a). A ground state can be chosen in a joint eigenspace of all M_i because of (1). Its marker pattern is nonconstant. On a prime ring any nonconstant binary pattern has translation orbit of length L: its orbit length divides L and cannot be one. Translating the chosen ground state through this orbit produces L states of the same energy in mutually orthogonal marker sectors. They are linearly independent. QED.

This ground-state multiplicity is incompatible with convergence of the complete low-energy spectrum to a fixed finite collection of compact unitary CFT sectors with finite vacuum multiplicity. It does not become a larger central charge. Restricting by hand to one marker sector or to compatible lengths changes the spectral question.

The theorem is not a statement about CFTs with an explicitly allowed infinite vacuum multiplicity or a divergent torus trace. Such a convention would need to be specified separately and would also change the modular counting assumptions in C7.

#### Application to the C7 selector

For a fixed selected period b, let e=e_b<0 be the optimal cell energy density from C7. A large ring of any length can be tiled by cells of size b except for one cell whose length lies between b and 2b-1. The excess energy over e_b L is bounded by

    C_b = max_(0<=r<b) [E_(b+r)-e_b(b+r)],                     (2)

which is finite for this fixed b. This proves (b) for the bare selector. The no-marker legal sectors have energy zero. C7's forbidden-bond estimate excludes lower-energy illegal sectors; for the bare selector the all-marker pattern is also uniformly costly. Thus (a) holds with a positive eta for the two constant patterns, and the theorem applies.

An encoded spectator construction inherits the same conclusion if it retains exact marker conservation and the two energy-density properties after adding its couplings. It is not necessary that the selected nonconstant ground pattern consist of only one defective cell. Primality alone gives the full orbit as soon as the constant patterns lose.

Consequently a proof that a frozen-marker encoding has the desired CFT on lengths divisible by b cannot meet the stronger all-ring requirement by appealing to a bounded defect energy. The defect position supplies an unbounded exact degeneracy. This is independent of whether the first-order code or second-order spectral gadget has been realized successfully.

### 2. Making an isolated defect mobile does not automatically fix this

Consider the following narrower effective model. A forced defect has a position x on a ring and a finite internal space of dimension s independent of L. Its translation-invariant Hamiltonian has finite-range hopping matrices,

    H_def=sum_x sum_(r=-R)^R |x+r><x| tensor A_r,
    A_(-r)=A_r^dagger.                                        (3)

Exponential hopping decay with an analytic Bloch matrix is equally sufficient. The Bloch matrix is

    h(k)=sum_r exp(ikr) A_r,

a Hermitian analytic matrix for real k. Let E_* be the minimum of its lowest band over the Brillouin zone.

Proposition. Such a defect has finite-size excitation energies above its ring ground state bounded by O(L^(-2)), or an increasing exact degeneracy. A nonzero linear cone cannot be the bottom of its finite-band spectrum.

Proof. Choose a normalized eigenvector u of h(k_*) at a minimizing momentum, with eigenvalue E_*. The scalar analytic function f(k)=u^dagger h(k)u obeys f(k)>=E_* for every real k and f(k_*)=E_*. Thus f'(k_*)=0 and f(k)<=E_*+C(k-k_*)^2 near k_*. For any fixed integer n, choose n distinct allowed ring momenta within O(n/L) of k_*. Their plane waves with internal vector u are orthogonal trial states, all with energy at most E_*+C_n/L^2. The Hamiltonian is diagonal in momentum, so the same bound holds throughout their span. The min-max principle gives at least n eigenvalues below that upper bound. The ring ground energy is at least E_*, proving the assertion. Since n is arbitrary, a fixed finite vacuum multiplicity cannot absorb all these levels after conformal rescaling by L. No assumption about analytic eigenvalue branches is needed. QED.

If this defect is decoupled from a critical spectator theory, the product spectrum contains extra quadratic or flatter levels alongside the conformal towers and is not purely a fixed finite collection of z=1 CFT spectra. This proves why an ordinary massive defect hopping term does not repair the frozen-selector construction.

The proposition does not apply without further work to a defect coupled to a gapless many-body bath. Its internal low-energy space may not have fixed finite dimension, and its effective dispersion need not follow from the analytic finite-band model (3). No universal theorem excluding a successful interacting repair is claimed. A complete counterexample must either control this coupling or avoid the forced-defect mechanism.


## 13. A mobile two-length selector

2026-10-27, resumed C8. The microscopic move, its conservation laws, and its first-order projection are proved. The projected-model spectra below are numerical. No conformal finite-u microscopic simulator, logical payload transport, or c-d counterexample is claimed.

### 1. An explicit local move

Use C7's four letters #,L,X,R and its allowed arrow pairs, with

    #=(R,L),   L=(L,L),   X=(L,R),   R=(R,R).

Add a marker hopping term to the same nearest-neighbor bond,

    h=K Q + (mu/2)(M tensor I+I tensor M)
      -t (|LX><XR|+|XR><LX|)
      -u (|#L><R#|+|R#><#L|),                                 (1)

where M=|#><#| and Q is the forbidden-pair projector. All coefficients are independent of ring length. The new move has identical external arrow endpoints on its two sides. It therefore preserves the allowed/forbidden status of its own bond and both adjacent bonds, just as C7's head move does. Every Q_i still commutes with H. Individual marker positions no longer commute with H, but total marker number sum_i M_i does.

The move #L -> R# shifts a marker one site to the right. The head in the cell to its right must not occupy the immediately adjacent site, since that site must be L. Its Hermitian conjugate moves the marker left.

### 2. Projection at a crossing

Fix b>=5 and set the unperturbed chemical potential to the C7 crossing between periods b and b+1,

    mu_c=2t-delta_b,
    delta_b=(b+1)a_b-b a_(b+1),
    a_l=2t[1-cos(pi/l)].                                      (2)

Then E_b/b=E_(b+1)/(b+1)=e_c<0. C7's convexity argument shows these are the only minimizing lengths. Consequently every tiling by b- and (b+1)-cells, with each head in its clock ground state, has energy e_c L. For fixed b the gap above this entire manifold is positive on ring lengths that admit such a tiling. It is the minimum of the internal clock gaps, the excess costs of other cell lengths, and the forbidden-sector gap. The manifold can be exponentially large; that does not invalidate the gap statement above the whole manifold.

For all sufficiently large L a tiling exists: choose an integer M between L/(b+1) and L/b. The number of long cells is then

    r=L-bM,    0<=r<=M.                                       (3)

In particular, every L>=b(b+1) admits at least one such M. Primality of L is no longer an obstruction to the tiling manifold.

The normalized head wavefunction in a cell of length l is

    phi_l(i)=sqrt(2/l) sin(pi i/l),    i=1,...,l-1.

A marker hop stays within the minimizing manifold only when it exchanges neighboring lengths (b,b+1) and (b+1,b). The left head overlap is

    A_b=sum_(i=1)^(b-1) phi_b(i) phi_(b+1)(i).                  (4)

On the right, shifting the marker removes its first site, so the overlap is sum_(j=2)^b phi_(b+1)(j) phi_b(j-1). Reflection of the sine wavefunctions makes this also A_b. Thus the projected off-diagonal matrix element is

    -u A_b^2.                                                 (5)

There is no diagonal first-order term. For b=5,8,16, A_b is respectively 0.9412156469, 0.9754576110, and 0.9934958020.

Within a fixed marker-number sector, label a cell by eta_j=0 or 1 according to its length b+eta_j. In the interior, (5) is the XX hopping of eta particles,

    H_gap=-J_b sum_j (tau_j^+ tau_(j+1)^-+h.c.),
    J_b=u A_b^2,    sum_j eta_j=r.                             (6)

On a ring one must also keep the reference marker position. A hop across the reference boundary shifts that position, so a fixed physical-translation sector gives a twist of the gap-chain boundary condition and a consistency constraint from relabeling the reference marker. Dropping this bookkeeping would lose finite-size levels. The script below instead builds the physical marker-pattern graph directly, with no reference-marker redundancy.

At physical momentum zero, the periodic XX ground wavefunction on gap sequences lifts to a positive wavefunction on physical marker patterns by cyclically identifying the reference. The usual positive fixed-r XX ground state is cyclically invariant, so the lift is well defined. It has energy

    E_XX(M,r)=-2J_b sin(pi r/M)/sin(pi/M),    0<r<M.             (7)

Every physical-pattern graph component admits this positive eigenvector. On a component where the hopping graph is connected, Perron-Frobenius makes it the unique component ground state. Connectivity is checked in the examples below; a general connectivity classification is not assumed in this note.

Equation (6) is a first-order coefficient, not a volume-uniform theorem for (1) at fixed u. Excited clock states and excluded lengths give higher-order terms. Their local control and the infrared phase need a separate argument.

### 3. Why a finite-density length fluid is different from one forced defect

In the thermodynamic approximation r/M=nu, the kinetic energy density obtained from (7) is

    e(nu)=e_c-(2J_b/pi) sin(pi nu)/(b+nu),    0<=nu<=1.         (8)

It is strictly lower in the interior than at either pure cell length. Its minimizing nu_* satisfies

    tan(pi nu_*)=pi(b+nu_*),    0<nu_*<1/2.                    (9)

To see uniqueness, differentiate sin(pi nu)/(b+nu). Its numerator is pi(b+nu)cos(pi nu)-sin(pi nu), whose derivative is -pi^2(b+nu)sin(pi nu)<0 on (0,1). It starts positive and ends negative, with its zero before 1/2. The minimum of (8) is therefore nondegenerate. Large b gives nu_*=1/2-1/[pi^2(b+1/2)]+O(b^-3).

This leading description has an extensive number of each cell length. It does not consist of an otherwise rigid crystal with one massive forced defect. That distinction is how it can avoid the finite-band defect obstruction in ring_obstruction.md. Changing L can be absorbed by small changes in M and r around their thermodynamic densities.

A plausible continuation is a c=1 Luttinger liquid of cell lengths. The XX description has exact marker-number conservation, and weak local corrections at an interior density have the usual spinless-fermion stability argument. Here that is only a continuum/RG proposal. Neither a thermodynamic local elimination theorem nor finite-u microscopic conformal estimators have been supplied. The first-order graph spectrum is evidence for the mechanism, not a proof of its completion.

### 4. Numerical graph checks

Run

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 mobile_selector.py

The code generates every binary marker pattern whose consecutive spacings are b or b+1, applies the legal marker hops, and diagonalizes the resulting sparse graph Hamiltonian. Energies are in units J_b. It also exhaustively checks the legality indicators of both local moves against every choice of their neighboring letters. All local checks pass.

For b=5 and nu=1/2 the results are

| L | M | r | Graph dimension | Components | E0 | L gap |
| :-- | :-- | :-- | :-- | :-- | :-- | :-- |
| 33 | 6 | 3 | 110 | 1 | -4.0000000000 | 2.385412 |
| 44 | 8 | 4 | 385 | 1 | -5.2262518595 | 2.340611 |
| 55 | 10 | 5 | 1386 | 1 | -6.4721359550 | 2.320290 |
| 66 | 12 | 6 | 5082 | 1 | -7.7274066103 | 2.309359 |

The ground energies agree with (7) within 3e-14. The displayed first gaps are consistent with a nonzero 1/L scale. They are not two independent central-charge estimators, and no central charge is inferred numerically from them alone. The density 1/2 is a convenient tested sector, not the optimal density (9) at mu=mu_c.

For comparison, keep only one long cell, r=1, with M=6,8,10,12. The physical graph is a cycle of length L=31,41,51,61. The ground energy is -2 and the gap is 4 sin^2(pi/L). Accordingly L^2 gap increases from 39.34345 to 39.44353, approaching 4pi^2. This is the quadratic isolated-defect behavior proved in ring_obstruction.md.

### 5. What this changes for the large-c route

The explicit move (1) offers a way to address the all-ring problem without keeping defects immobile or assigning an ordinary massive hopping band to only O(1) defects. It replaces a unique cell length by a finite-density quantum mixture of two lengths.

A logical payload is an additional problem. One needs compatible code isometries V_b and V_(b+1), both with the same logical space and the required boundary compressions. In addition, the marker hop must act between neighboring codes as a scalar multiple of the identity on their payloads, or as a known harmless transport unitary. A generic pair of abstract Stinespring codes does not satisfy that condition. Adding an idle physical site at one end of one code is not enough unless both possible directions of cell growth and shrinkage are handled.

Even with this transport property, the coupled length and payload modes require an infrared analysis. Their central charges cannot be added merely because both isolated leading targets are critical. This note resolves neither the repeated-bond code construction nor the coupled conformal spectrum. It supplies a concrete alternative to the frozen selector and a testable extra code-compatibility condition.


## 14. Final resumed-run status and next steps

No conclusive answer under the task's definition was reached. There is no verified c>d-1 chain, no bound for the full class, and no construction of unbounded c at fixed d. The constructions in Sections 10-13 are finite-dimensional or first-order effective statements with their hypotheses stated. The RG interpretations are separate assumptions or heuristics.

The additional local-parent search was also negative. `local_parent_scan.py` tests the smallest q=2,m=2 positive-square-root Choi code on five qubits, for all six orders of its three environment qubits between its two ports. For each order it solves the real linear constraints (H-E)V=0 for every uniform Hermitian two-site bond, including all oriented one-site terms and excluding only the scalar. Each 15-column constraint matrix has full rank, with smallest singular value about 1.408. Thus none of these specific embeddings admits a nonconstant repeated-bond parent, numerically. Allowing 51 independent nearest-neighbor Pauli coefficients leaves no solution in five orders and two null directions in the remaining order. Those directions have zero code energy and are traceless nonzero Hamiltonians, so their nonzero combinations cannot have the code as a ground space. This is a numerical exclusion of the tested embeddings. An environment unitary, a different dilation, extra ancillas, or a different code was not searched.

Reproduce that check with `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 local_parent_scan.py`. It writes `local_parent_scan.json`. The script `c8_resume_checks.py` reproduces the new leakage, flat-parent, and spectral-parent checks. The script `mobile_selector.py` reproduces the length-fluid graph checks. All three use only numpy and scipy and complete in seconds on this machine. Python compilation checks passed. The pre-existing run 1 checks were not repeated because their code and claims were not changed.

The next steps are now ranked as follows.

1. Construct compatible local payload codes at lengths b and b+1, with the same logical space, both required boundary maps, and a controlled action of the marker hop between the two lengths. Even a small example is useful. This combines the remaining locality problem with the all-ring requirement; solving only a frozen single-length code is insufficient for the stronger spectrum definition.
2. Prove a thermodynamic local effective-Hamiltonian reduction for a candidate and analyze the coupled length and payload modes. The positive second-order singlet sign in Section 11 removes one obstruction, but higher-order anisotropies and relevant couplings still need control. An incommensurate target with exact parity is another option.
3. For the separate c<=r(h) conjecture, continue the explicit rank-one construction in Section 6. Its two missing premises, local low-energy reduction and infrared stability, remain unproved. This route bypasses small-site programmability and does not require a selector, but even success would not settle c<=d-1 for that large physical d.
4. For an upper bound, strengthen the locality tradeoff in Section 10 into a restriction that cannot be evaded by shrinking the block gap and exchange scale. No such uniform restriction was obtained. Port rank, exact-encoding arguments, and finite-order signs alone do not bound the general conformal class.

No network search, cluster run, full finite-coupling microscopic CFT calculation, or independent referee check was performed in the resumed run. No project files outside this folder were modified. An incidental temporary stdout capture in /private/tmp was created and removed during the local-parent check. The project-level STATE.md remains untouched as required by the task's read-only restriction.
