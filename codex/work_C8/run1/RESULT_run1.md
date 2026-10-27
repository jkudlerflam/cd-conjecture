# C8 result

2026-10-27. Closed-book work by Codex, GPT-6. No network access, no subagents, and no files outside this folder were modified.

## Verdict

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

G1 is at `../../contrib/unbounded_c/G1/G1_result.md`. It already disproves an unrestricted claim that diagonal terms must precede noncommuting exchange. Its parent README still describes that hierarchy as an obstruction without incorporating the result.

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

## 9. Ranked next steps

1. Construct a local block Hamiltonian realizing a specified joint boundary map. The first concrete target is (1), with the generic quadratic target (27) and its on-site chemical potentials. A two-sweep or programmable history head is a candidate because the one-pass theorem does not cover it. The construction must state the same microscopic h on every internal bond and show the code is its low-energy space. A generic Stinespring isometry or whole-block projector does not suffice.
2. For any such code, prove a local effective-Hamiltonian expansion with controlled errors and compare those errors with the scale alpha^2 lambda. Establish its IR using an applicable multicomponent Luttinger stability result, or tune every dangerous coupling and check the CFT numerically with independent estimators. The free target's conformality alone does not prove this step.
3. In parallel to a first-order code construction, attempt local realization of the joint matrix-valued boundary resolvents in Section 4. The algebraic second-order hierarchy has been removed, so locality and stability are the issues. Exact large SU(N) protection is excluded by Theorem 3 and should not be used as a shortcut.
4. Handle every microscopic ring length, or establish the intended admissible-sector definition explicitly. A successful commensurate subsequence with immobile forced defects elsewhere does not settle the stronger all-ring formulation in the task.
5. For an upper bound, seek a locality restriction on these joint maps or boundary resolvents that goes beyond their complete positivity and port rank. Theorem 1 and (23) show why those two conditions alone cannot provide the missing restriction. No new unconditional c<=f(d) bound follows from this run.

For the separate rank conjecture, a useful intermediate target is to prove or disprove premises (A) and (B) for (31). It avoids the fixed-small-d block-realization problem altogether. A resolution could settle c<=r(h) while leaving the original c-d question open.

The c-d conjecture and the existence of any fixed-d bound remain open after these calculations. The run narrows the missing step to a local realization and its infrared control, with explicit algebraic targets and two restricted no-go theorems available for the next attempt.
