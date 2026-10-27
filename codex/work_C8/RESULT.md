# C8 result

2026-10-27. Closed-book work by Codex, GPT-6. No network access or subagents. All project changes are confined to this folder.

## Verdict

Fourth-run update. The verdict remains nonconclusive under the definition
in the task. There is no verified microscopic counterexample, no bound
for the full class, and no fixed-d construction with unbounded c.
Sections 20-24 contain this run's new results and supersede the earlier
next-step rankings. All earlier proofs and checks remain below.

The strongest new obstruction is a boundary-Schmidt-rank bound. If a
code contains one state of Schmidt rank s across an active q-dimensional
port, an exact scaled Clifford family has at most 2qs-s^2 generators.
One product vector therefore forces c_target<=q-1 for the proposed even-
channel Majorana target. This excludes complete ground codes of positive
bonds with matrix rank below q, and gives a robust obstruction to their
approximate joint ports. Commuting-projector parents and the specified
injective two-site MPS parent class are also excluded as implementations
of the proposed high-channel code. These are restricted construction
no-gos, not a general c bound.

The open-chain frustration-free witness length is now explicit and sharp
for qubits, K_FF(2)=4. A positive bond with a three-site zero state and no
four-site zero state proves sharpness. A q-state extension gives
K_FF(q)>=q+2 for q>=2. Neither statement covers frustrated conformal chains
or marker-dependent boundary penalties.

A constructive test supplies a repeated positive graded-swap parent on
q=n+1 states, with a complete 2^n-dimensional ground code and controlled
maps at both ends. Its normalized error from 2n Clifford generators is
1-sqrt(1-(n-1)/N) on a cell of length N. It only approaches c_target=q-1
before paying selector costs, so it does not produce a violation.

A separate conditional bound uses the connected bond itself. If the
Schmidt singular values have multiplicities n_j, the dynamical image of
any exact on-site symmetry algebra embeds in the direct sum of so(n_j).
When the corresponding microscopic currents exhaust the IR stress tensor,
this gives c<=min(d^2-1,sum_j n_j(n_j-1)/2). Gaussian modes resolved by
commuting on-site charges obey c<=sum_j floor(n_j/2)<=floor(r(h)/2).
Neutral and emergent sectors remain outside these hypotheses.

Earlier work supplies joint two-port dilations, second-order synthesis
of arbitrary logical bonds, an exact d=4 mobile selector with c=1,
and conditional vacancy transport of frustration-free codes. It also
identifies the finite-coupling stability and incompatible-ring issues.
The unresolved construction still needs a sufficiently entangled
high-channel code with a repeated local parent, number-sector selection,
and thermodynamic infrared control. The new tests concern only finite
codes and symmetry algebras. No new microscopic CFT was measured.

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


## 15. Exact finite-hopping selector and all-ring conformality


Status: the Hamiltonian mappings, sector bound, entropy identity, and finite-size
energies below are proved. The identification of a free Fermi sea with its
compact c=1 CFT is the stated standard continuum input. This supplies a
microscopic c=1 model at d=4, not a counterexample to c <= d-1.

## 1. Hamiltonian and legal arrows

Use the four letters #,L,X,R with arrows

    #=(1,0), L=(0,0), X=(0,1), R=(1,1).

Write Q for the projector onto pairs whose touching arrows disagree and
M=|#><#|. Set a>0 and

    h = K Q + (mu/2)(M tensor I + I tensor M)
        -a(|LX><XR| + |XR><LX|)
        -a(|#L><R#| + |R#><#L|).                         (1)

This is the preceding mobile-selector Hamiltonian at t=u=a, with no small-u
projection. It is one-site translation invariant and nearest neighbor.
Every Q_i commutes with H because each allowed move preserves its three
affected legality indicators. Thus the full Hilbert space splits by the
binary violation pattern epsilon_i=Q_i.

In the legal sector let b_i be the right arrow on site i. Its left arrow is
b_(i-1). Every binary sequence b gives one legal basis vector, including
on odd rings. The legal Hilbert space has dimension 2^L.

A hop at bond i flips b_i precisely when b_(i-1) != b_(i+1). Both orientations
have the same amplitude -a. Let Z_i|b_i> = (-1)^b_i |b_i>. The restricted
Hamiltonian is therefore

    H_legal = -(a/2) sum_i X_i(1-Z_(i-1) Z_(i+1))
              +(mu/4) sum_i (1-Z_(i-1) Z_i).             (2)

The chemical term follows because # counts the 1-to-0 transitions. On a
ring they are half of all transitions. This step would fail on an open
interval without a boundary term.

## 2. All forbidden sectors stay uniformly above the ground state

In a fixed violation sector, the left arrow of site i is

    a_i = b_(i-1) xor epsilon_(i-1).

Consequently every violation sector also has dimension 2^L and the same
binary-coordinate Hilbert space. The comparison is not a variational
identification of a few basis states.

For a violation on bond i, that bond's hopping is absent. This changes the
operator by norm at most a. The hopping coefficient on bond i+1 can also
change, because its left outside arrow a_(i+1) has flipped. The two
conditional flips differ by norm at most a. All other kinetic terms agree,
apart from the analogous effects of other violations. The marker projector
on site i+1 changes by norm at most one. It costs at most |mu| in operator
norm. The triangle inequality gives, for k violations,

    ||H_epsilon - K k I - H_legal|| <= k(2a+|mu|).        (3)

This remains valid when affected neighborhoods overlap. Thus

    E_min(epsilon) >= E_min(legal) + k(K-2a-|mu|).        (4)

Choose K>2a+|mu|. The legal ground state is then a ground state of the full
d=4 chain, and every illegal sector is separated from it by a positive,
L-independent gap at least K-2a-|mu|. This includes all ring lengths L>=3.
The small L=2 convention, in which the same unordered edge occurs twice,
is irrelevant to the limit and is not included.

## 3. Free domain walls, including both ring twists

Define wall occupations w_i=b_(i-1) xor b_i. Their total W is even on a ring.
Every even wall pattern has two preimages, related by the global flip
G=product_i X_i. The move in (2) exchanges w_i,w_(i+1)=(1,0) and (0,1).
It is nearest-neighbor hard-core hopping with amplitude -a, and the chemical
potential is mu W/2.

Fix the orientation b_0 when reconstructing b from w. A hop across that
reference changes b_0; every other wall hop leaves it fixed. In a sector
G=g, g=+1 or -1, that crossing has amplitude -ag. This is the hard-core
particle boundary twist. After Jordan-Wigner, even W gives fermion
boundary phase -g. The complete spectrum is the union over

    W even, g=+1,-1,
    k_j = 2pi(j+phi_g)/L,
    phi_+=1/2, phi_-=0,
    E = sum_(occupied k) [-2a cos(k) + mu/2].             (5)

Both twists are necessary. Counting all these states gives
2 times 2^(L-1)=2^L, the whole legal Hilbert space. This also proves there
are no omitted extensive zero-energy sectors.

For |mu|<4a the dispersion has two simple Fermi points

    cos(k_F)=mu/(4a),       v=2a sin(k_F)>0.               (6)

Linearization at them gives a massless complex fermion. The even-W
projection and the two twists impose finite global sector restrictions;
they change the compactification and allowed towers but not the local
stress tensor or its c=1 central charge. This last free-fermion continuum
identification is the standard input. At |mu|=4a the slope vanishes and the
band edge is quadratic, so those endpoints are excluded.

On odd or otherwise incommensurate rings, W and the filled momenta are
rounded within O(1) of their thermodynamic values. Their energies differ
at order 1/L. Each Fermi level has at most two orbitals with a given energy,
so ground-state multiplicities stay bounded, independently of L. No
immobile defect with L possible positions occurs. Irrational k_F gives
the usual size-dependent charge-sector offsets of a lattice Fermi sea;
there is no claim that its vacuum offset is identical for all L.

## 4. Energies, towers, and entropy at mu=0

Take L divisible by four. The positive ground state is in g=+1 with
W=L/2. Summing the occupied cosines gives

    E0(L) = -2a/sin(pi/L)
          = -(2a/pi)L - (pi a)/(3L) + O(L^-3).           (7)

With v=2a, the Casimir coefficient is -pi c v/(6L), giving c=1.
The lowest state in g=-1 has energy -2a cot(pi/L). Its gap is

    2a tan(pi/(2L)) = pi a/L + O(L^-3),                  (8)

so it has scaling dimension x=1/4. Particle-hole excitations near either
Fermi edge have energy 2pi v/L times an integer, with the ordinary
fermionic partitions producing the oscillator towers. Changing even W
by 2 has leading energy 4pi a/L, or x=1. Thus both the first intersector
gap and the oscillator scale are linear in 1/L.

There is also an exact entanglement identity. Let |psi_w> be an even-wall
state, and take its g=+1 legal-arrow lift. For a nonempty proper interval
A of physical sites, write n_A for its wall parity and s for its initial
arrow. The basis change within A identifies its allowed letter strings
with (w_A,s), and the complement with (w_B,s xor n_A). The lifted state is

    (1/sqrt(2)) sum_(w_A,w_B,s) psi(w_A,w_B)
       |w_A,s>_A |w_B,s xor n_A>_B.                     (9)

Tracing B forces equal w_B. Global even parity then forces equal n_A
on the bra and ket, and hence equal s. Therefore

    rho_A^physical is unitarily equivalent to
       (rho_A^walls/2) direct-sum (rho_A^walls/2),
    S_A^physical = S_A^walls + log(2).                  (10)

Zero eigenvalues from illegal interval strings do not affect the entropy.
This proves equality of the logarithmic entropy coefficients, rather
than estimating the effect of the encoding by an unspecified constant.

## 5. Independent checks and limitations

Run `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 equal_hopping_selector.py`.
The output is `equal_hopping_selector_checks.json`.

The script compares every eigenvalue of (2) at L=5,6,7,8,9 and mu=0.37 with
the independent occupation sum (5), with maximum error 1.96e-14. It tests
(3) for all 15,31,63 nonempty violation patterns at L=4,5,6 respectively. It constructs the many-body Slater wavefunction and
computes physical-site reduced density matrices, without storing a 4^L
state vector. The XX wall state is the known c=1 calibration for entropy.

At L=20, fitting the interval entropy to log[(L/pi)sin(pi ell/L)] for
ell=5,...,10 gives c=1.003981. The half-chain entropy slope from L=16 to 20
gives c=0.992700. Each equals its XX control to numerical precision, as
(10) predicts. The Casimir fit at L=40,48,64,80,96,128,160, with independently
known v=2 and a 1/L^4 correction in E0/L, gives c=0.999999704.

This resolves the finite-u CFT question only at the equal-hopping point.
It does not prove that the long-cell, small-u regime is in the same phase,
and it contains no high-dimensional logical code. Moving from this
soluble selector to a large-c payload still requires a microscopic
Hamiltonian and a separate stability analysis.


## 16. Two-length code compatibility by padding


Status: proved finite-dimensional code construction. No local parent or
finite-coupling CFT is asserted.

## 1. Padding an arbitrary code

Suppose V:C^D -> (C^q)^(tensor k) is an isometry, k>=2, whose first and last
site compressions are

    V^dagger F_a^(first) V = A_a,
    V^dagger F_a^(last) V = B_a.                          (1)

No commutation relation between A and B is required. Add one new local
vacancy symbol |0> orthogonal to C^q, so the new site dimension is q+1.
Extend F_a by zero on that vacancy. Regard V as living on k occupied sites.
For p=1,2 define

    W_p = (V tensor |0>^p + |0>^p tensor V)/sqrt(2).       (2)

The two summands have different, orthogonal occupation patterns. Hence
W_p^dagger W_p=I. In the first summand only the left physical endpoint is
occupied; in the second only the right endpoint is occupied. A one-site
operator cannot connect the two occupation patterns. It follows that

    W_p^dagger F_a^(first) W_p = A_a/2,
    W_p^dagger F_a^(last) W_p = B_a/2,    p=1,2.          (3)

Both lengths k+1 and k+2 therefore have the SAME boundary maps. This is
stronger than padding a single code on just one side, which would erase
one of its ports.

## 2. Moving the boundary preserves both logical states after projection

Place a short cell followed by a long cell on 2k+3 physical sites. Move the
boundary one site to the right and leave the spectator variables untouched.
The old and new isometries are W_1 tensor W_2 and W_2 tensor W_1. Then

    (W_2 tensor W_1)^dagger (W_1 tensor W_2)
         = (1/4) I_(D tensor D).                         (4)

To prove this, expand each side into its four occupation patterns. Label
the first site zero. The two cores on the input start at

    (0 or 1, k+1 or k+3),

and the two output cores start at

    (0 or 2, k+2 or k+3).

The only common occupied-site pattern has core starts (0,k+3). In it the
site transferred across the boundary is vacant. The overlap of each
unchanged core is V^dagger V=I, and the four factors 1/sqrt(2) give 1/4.
All other terms vanish by orthogonality of vacancy and occupied symbols.
The same statement holds for a hop in the opposite direction by taking
the adjoint. Adding a projector onto a vacancy at the transferred site
does not change (4).

This is a projected transition, not an exact invariant encoding. Its
amplitude is less than one; a physical hop also has components outside
the desired code spaces. It does not evade the earlier exact-product-code
theorem, whose hypotheses require that leakage to vanish.

## 3. Leading interaction is independent of cell length

For a boundary interaction sum_a F_a tensor F_a, every pair p,p' in {1,2}
has the same compression

    (W_p tensor W_p')^dagger h_boundary (W_p tensor W_p')
         = (1/4) sum_a B_a tensor A_a.                   (5)

If A_a=alpha Gamma_a and B_a=alpha i Gamma_a P, this is the previous
quadratic Majorana target with hopping coefficient lambda alpha^2/4.
If selector clock overlaps contribute a scalar A_b^2 as in the mobile
selector, the marker-hopping amplitude acquires the additional factor
1/4 from (4). The two length variables carry XX hopping tensor the
identity on the payload. The payload interaction (5) is independent
of both length variables.

Thus in each fixed-cell-number sector, and before higher-order corrections,
the bulk effective Hamiltonian is a sum of an XX length chain and the
logical target on the same cyclic cell index. On a physical ring their
translation and twist sectors are correlated by cyclic relabeling of the
cells. They are not unrestricted tensor products of finite-size spectra.
The local stress tensors would still add if this leading description
were the actual stable infrared theory. That last premise is not proved.

There is no local-parent theorem hiding in (2). The projectors
I-W_p W_p^dagger are arbitrary (k+p)-site operators. Even if V had a
nearest-neighbor parent, making the coherent two-placement ground spaces
in (2) is another problem. Hopping vacancies through the core introduces
additional positions and modifies its boundary moments. The construction
only removes a linear-algebra incompatibility between two lengths.

## 4. Dimension bookkeeping and a check

One may tensor this new spectator with the four-letter selector, at
d=4(q+1). For example the graded construction in the earlier C8 report
allows q=17 and m=2*8*9=144 odd channels. Its logical dimension is 2^72,
its formal target charge is m/2=72, and its Choi dilation fits in

    k = ceil(log_q(q^4 2^72)) = 22

spectator sites. The padding gives lengths 23 and 24 and d=72. The formal
payload alone would exceed d-1=71; a length mode would add one more charge
if it survived. These numbers establish no microscopic example. In
particular a uniform two-site parent for either padded code has not been
given, and the literal selector length counts must be aligned with the
chosen placement of its markers before constructing such a parent.

Run `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 padding_transport.py`.
It uses the earlier q=2,m=2 simultaneous Choi code, with its environment
placed between its endpoints, as a k=5 core. It checks both padded
isometries, both sets of endpoint compressions, the full two-code overlap
(4), and all four interactions (5). These are matrix identities on the
physical spectator spaces, not tests of the formulas substituted into
themselves. The output is `padding_transport_checks.json`.


## 17. Local vacancy parent and the exact-cat obstruction


Status: finite-dimensional theorems, proved below. The explicit local
example uses a ferromagnetic code and has no large-c claim. The final
spin-charge discussion is conditional.

## 1. The exact cat padding cannot be an isolated local ground code

There is a necessary qualification to `padding_transport.md`. Its exact
two-placement code has no strictly local parent with only that code as
ground space, once the core is longer than the interaction range.

Let A_p=V tensor |0>^p and B_p=|0>^p tensor V. These isometries have
orthogonal occupation supports. Suppose every Hamiltonian term omits at
least one site where the two supports have different occupation. Then

    A_p^dagger H B_p = 0.                                  (1)

For the proposed p=1,2 padding and nearest-neighbor terms, k>=2 suffices:
already the first and last sites disagree, and no bond contains both.
Matrix elements vanish because at an untouched disagreeing site one
branch is vacant and the other occupied.

Suppose every vector W_p psi=(A_p psi+B_p psi)/sqrt(2) is a ground vector
of H with ground energy E. Positivity of H-E gives

    0 = <W_p psi|(H-E)|W_p psi>
      = ( <A_p psi|(H-E)|A_p psi>
          +<B_p psi|(H-E)|B_p psi> )/2.                    (2)

Each summand is nonnegative and must vanish. A positive matrix with zero
quadratic form annihilates the vector. Hence every A_p psi and every
B_p psi separately belongs to the ground space. The degeneracy is at
least 2D, not D. The antisymmetric combinations are also ground vectors.

This argument does not assume translation invariance, frustration
freeness, or commuting terms. The pure two-placement construction gives
coherent projected transport algebraically, but a local realization must
include intermediate vacancy positions or some other dressing. The next
construction does so.

## 2. A local extension theorem for any frustration-free core

Assume the k-site core has a uniform open-chain parent

    H_core = sum_(j=1)^(k-1) h0_(j,j+1),  h0>=0,
    ker(H_core)=range(V).                                 (3)

The positivity of each summand implies h0_j V=0. This is the extra
hypothesis. A general frustrated parent does not satisfy it.

Enlarge each site from C^q to C|0> direct-sum C^q. Write n=|0><0| for
vacancy number and extend h0 by zero outside the two-occupied-site sector.
On N=k+p sites, for p=1 or 2, use the SAME repeated bond

    h_ext = U(n tensor I + I tensor n - 2n tensor n)
            -w sum_(a=1)^q (|0,a><a,0|+|a,0><0,a|)
            +g h0,          U,w,g>0.                     (4)

Restrict to exactly k occupied sites. The core's q-state labels cannot
pass each other under vacancy hopping on an open chain. Therefore this
sector has the exact ordered-coordinate decomposition

    C^(binomial(N,p)) tensor (C^q)^(tensor k),             (5)

where the first factor specifies hole positions and the second lists the
occupied labels in their physical order. This is a unitary change of
basis, not an approximation or a low-energy restriction.

The first two terms of (4) become H_charge tensor I. The charge matrix
has diagonal energy U times the number of occupied-vacant interfaces
and off-diagonal -w for every allowed nearest-neighbor hole hop. Its
configuration graph is connected. Perron-Frobenius gives a unique,
strictly positive normalized ground vector chi_p and energy e_p. Reflection
invariance makes chi_p reflection symmetric.

Let D_j be the diagonal charge projector onto configurations in which
the j-th and (j+1)-th occupied labels are physically adjacent. The full
Hamiltonian in (5) is

    H_ext = H_charge tensor I + g sum_j D_j tensor h0_j.   (6)

Since all terms in the second sum are positive, H_ext>=e_p. The states

    V_p psi = chi_p tensor V psi                          (7)

attain e_p, because every h0_j annihilates V.

Conversely, a ground state must minimize the first term and annihilate
the second, so it has the form chi_p tensor phi. Since chi_p is strictly
positive, <chi_p|D_j|chi_p> is strictly positive for every j. Its energy
above e_p is zero only if every <phi|h0_j|phi>=0, hence phi is in the
core kernel. Thus (7) is the COMPLETE ground space in the fixed-k sector.

For every fixed k,p,U,w,g>0 the gap above it is positive, by finite
dimension and the exact kernel characterization. This is not a bound
uniform in k. Its small-w splitting can be very small, which matters
when coupling the cells. The theorem makes no claim about the ordering
of different particle-number sectors; a full selector must control them.

The construction preserves an arbitrary frustration-free core code, not
just a classical or product code. It uses a single bulk bond at both
lengths, with no length-dependent endpoint term.

## 3. Endpoint maps and scalar boundary transport

Extend a core one-site observable F by zero on the vacancy. If physical
site 1 is occupied, it contains the first ordered label, regardless of
all hole positions. A number-preserving endpoint observable leaves those
positions unchanged. Consequently

    V_p^dagger F_first V_p = r_p A,
    V_p^dagger F_last V_p = r_p B,
    r_p = 1-Prob_(chi_p)(first site vacant),               (8)

where A,B are the original core endpoint compressions. Equality of the
two coefficients uses reflection of chi_p, not reflection of the core.

The leading boundary-shift overlap is also scalar for arbitrary V.
On its input, the left cell has p=1 and the right p=2; after shifting
their cut right, they have p=2 and p=1. To keep k occupied labels in
each cell, the transferred site must be vacant. No ordered logical label
crosses the cut. Every surviving core overlap is V^dagger V=I. Thus

    (V_2 tensor V_1)^dagger (V_1 tensor V_2)
        = beta I_(D tensor D),                           (9)

with beta=a_left a_right>0. For N1=k+1, N2=k+2,

    a_left = sum_(s=0)^(N1-1) chi_1({s}) chi_2({s,N1}),
    a_right = sum_(s=0)^(N1-1) chi_1({s}) chi_2({0,s+1}). (10)

All amplitudes are real and positive. Reflection gives a_left=a_right.
These formulas are in the physical spectator coordinates after the cut
has moved. They do not assert that the microscopic hop has no leakage.

For fixed k, as w/U tends to zero the charge states approach the equally
weighted superposition of the two endpoint hole droplets. To check this,
at w=0 these are the only configurations with one interface, with energy
U; all others have at least two interfaces. In the reflection-even sector
the ground vector at w=0 is unique, and the gap there is at least U.
Finite-dimensional perturbation theory gives an excited-configuration
amplitude O(w/U). Reflection fixes the two endpoint weights to be equal.
It follows that

    r_1,r_2 = 1/2 + O((w/U)^2),
    beta -> 1/4.                                         (11)

Constants can depend on k. The odd reflection partner is separated only
by a high-order tunneling splitting; it is not controlled by the U gap
in the even sector. The coupling between cells must be below this actual
gap, not merely below U.

Unlike exact cat padding, the local model generally has r_1 != r_2 at
finite w. Hence its payload bond is r_p r_p' times the core bond. The
length/payload coupling can be made small at fixed k, but is not zero.
Its infrared effect requires analysis. No tuning that makes r_1=r_2
with this same bond has been proved.

## 4. Explicit check with a quantum core

`vacancy_parent.py` takes q=2, k=3, and h0=(I-swap)/2. The core kernel is
the four-dimensional symmetric spin-3/2 representation. It is quantum:
the compressed X,Y,Z endpoint operators do not commute. This example
does not have an extensive Clifford family or a large target charge.

For w=0.15,0.3,0.5, U=1, g=0.7 the script builds the repeated physical
bond (4) on both lengths, diagonalizes each fixed-k sector, and checks
the entire ground code (7), its degeneracy, the endpoint maps (8), and
the two-cell transition (9). The data are in `vacancy_parent_checks.json`.

Run `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 vacancy_parent.py`.

## 5. Why this still does not finish the conjecture

The theorem reduces two-length local compatibility to the one-length
problem in the frustration-free subclass. No uniform nearest-neighbor
frustration-free core with the required many-channel endpoint maps has
been constructed here. The earlier Choi code has an arbitrary parent
projector and cannot be inserted into (3) on that basis alone. Its six
tested local embeddings already failed the parent test in the prior run.

The contribution on programmable cells proves a useful finite-history
interface, but does not provide the microscopic repeated-bond core
needed by (3), with its required endpoint compressions. Neither that
interface nor this extension establishes a stable large-c chain.

Vacancy doping of an entire hypothetical counterexample is another
conditional idea. On C|0> direct-sum C^d, nearest-neighbor particle-vacancy
hopping preserves the order of all occupied labels. At large hopping
compared with occupied-label interactions, first-order charge projection
gives a squeezed label chain. If spin-charge separation and infrared
stability held, it would add one charge mode, taking (d,c) to (d+1,c+1)
and preserving c-(d-1). Its particle number could adjust to a multiple
of a preferred cell period on every large physical ring.

This is not an all-ring completion theorem. The charge sector is gapless,
so finite-volume perturbation theory has denominators O(1/L) and cannot
be used uniformly as L grows. Exact factorization (5) at zero label
interaction does not make the interacting squeezed Hamiltonian exact.
The earlier direct-sum contribution already uses this mechanism at the
RG level. Without an additional thermodynamic stability input, it does
not close the present construction, and no further numerical claim is
made for it.


## 18. Finite witness for uniform frustration freeness


Status: proved, using the Hilbert basis theorem for a polynomial ring in
finitely many variables. This is not a bound on the central charge of a
general frustrated chain.

## 1. Statement

Fix the site dimension q. There exists a finite integer K_FF(q) such that
for every positive semidefinite Hermitian bond h on C^q tensor C^q,

    H_N(h) = sum_(i=1)^(N-1) h_(i,i+1),

having a zero-energy state on the open chain of length K_FF(q) implies
having a zero-energy state on every finite open-chain length. There are
no additional boundary terms. The integer is uniform over all bond
coefficients, gaps, and energy scales.

Equivalently, a fixed bond with a finite maximum frustration-free open
length has maximum length less than K_FF(q). No numerical value or
effective algorithm for K_FF(q) is supplied here.

## 2. Proof

Parameterize a Hermitian q^2 by q^2 matrix h by its q^4 real coefficients
x=(x_1,...,x_(q^4)), using real diagonals and the real and imaginary parts
of upper-triangular entries. For each N>=2 define

    f_N(x)=det(H_N(h(x))).                                (1)

This is a polynomial with real coefficients. Complex off-diagonal entries
cause no problem: the symbolic determinant is real on all real x by
Hermiticity, so its imaginary polynomial vanishes identically.

Let I be the ideal generated by all f_N, N>=2, in the polynomial ring
R[x_1,...,x_(q^4)]. The Hilbert basis theorem states that this ring is
Noetherian, meaning every ideal has a finite generating set. Each member
of such a finite generating set for I is itself a finite polynomial
combination of the original f_N. Taking the largest index used in all
these finitely many expressions gives a finite K with

    I = <f_2,...,f_K>.                                    (2)

In particular every later determinant has an exact polynomial identity

    f_N = sum_(j=2)^K p_(N,j)(x) f_j(x).                  (3)

Now specialize to any h>=0 with ker H_K nonzero. Take a normalized zero
vector on K sites. Every term h_(i,i+1) is positive and has zero
expectation, hence annihilates this vector. Its reduced density matrix
on the first j sites has zero expectation of H_j for every j<=K.
Since H_j>=0, the support of that density matrix lies in ker H_j and
is nonempty. Therefore f_j(h)=0 for all 2<=j<=K.

Equation (3) implies f_N(h)=0 for every N. Each H_N is positive, so a
zero determinant means it has a zero-energy ground state. The definition
of the polynomial sequence depends only on q, making K uniform. QED.

Positivity is used in passing from one length K to all shorter lengths.
The finite-generation identity alone does not make a single vanishing
determinant imply all later determinants for an arbitrary Hermitian h.

## 3. Implication for one proposed length-selection mechanism

The local vacancy extension in `vacancy_parent.md` has the positive
uniform core bond h0. A proposal to choose the number of occupied labels
by taking the LONGEST frustration-free core therefore cannot make that
number unbounded at fixed q. Its cutoff is below K_FF(q).

For a first-order Clifford target living on that core, its code dimension
D obeys D<=q^k. An even family of m Clifford generators needs
D>=2^(m/2). Hence in this particular longest-core construction,

    c_target=m/2 <= log_2 D <= k log_2 q
                  < K_FF(q) log_2 q.                    (4)

This is a bound on the stated free target and selection mechanism. It is
not a general c bound at fixed q. It does not bound the charge of an
arbitrary interacting target encoded on D dimensions.

The simpler attempt to select a partially occupied cell with just the
charge attraction and a vacancy chemical potential also fails whenever
the fully occupied core at that length is frustration free. For U>=w,
the mixed occupancy block of the charge bond is

    [[U,-w],[-w,U]] >= 0,

and the two equal-occupancy blocks have zero charge energy. With nu sum n
added, nu>0 favors the zero-hole, zero-core-energy state; nu<0 favors the
all-hole state, which saturates the lower bound nu N; and at nu=0 both
extremes have zero energy while a partial droplet at w<U costs energy.
This explains why the fixed-p hypothesis in the vacancy extension cannot
be silently removed. A distinct number-selection mechanism is needed.

## 4. Why it does not settle the main question

A boundary computation may impose different conditions at the two ends.
Then frustration freeness at length N need not imply it at smaller
lengths. The nested-kernel step of the proof fails, even though each
individual determinant is still a polynomial. Such boundary conditions
can be supplied by marker-dependent bonds in a larger translation-
invariant chain. They are outside this theorem's uniform open-bulk
hypothesis.

Likewise, the C7 selector minimizes negative, length-dependent kinetic
energies and is not a test for zero energy of a positive uniform core
bond. It is not constrained by K_FF(q). Nor are continuum finite-size
asymptotics polynomial equalities of the form (1). There is no justified
step from this Noetherian argument to a universal conformal onset length
K(d), and none is claimed.

The theorem isolates one way in which a proposed unbounded-length
frustration-free construction can fail. Boundary conditions, frustrated
energetics, and approximate codes remain possible ways around it.


## 19. Third-run verdict, checks, and next steps

No conclusive answer under the task's definition was reached. In particular,
there is no verified h with c>d-1, no full-class bound c<=f(d), and no proof
of unbounded c at fixed d. The exact c=1 selector and the local vacancy
extension are proved restricted constructions. The q=17, d=72, c=72
numbers remain resource accounting, not a candidate Hamiltonian with a
verified continuum limit.

The new computations are small and reproducible with Python, numpy, and
scipy. No cluster calculation, network access, or subagents were used.
Only files in this folder were changed.

- `equal_hopping_selector.py` checks complete spectra at L=5 through 9,
  every nonempty illegal pattern at L=4,5,6, entropies at L=8,12,16,20,
  and analytic finite-size energies at L=40 through 160. The entropy
  calibration is the XX wall chain, with an exact additive log(2)
  relation to the physical entropy. The independent Casimir estimator
  uses the velocity from the dispersion.
- `padding_transport.py` checks the q=2,m=2 Choi code at padded lengths
  6 and 7, including the full two-code overlap and all four choices of
  adjacent lengths.
- `vacancy_parent.py` constructs a positive repeated physical bond for
  a three-qubit symmetric core, at lengths 4 and 5. It checks the ground
  code, spectrum, endpoint maps, and coherent projected transport for
  three hopping strengths. Its first run failed a 1e-10 endpoint
  tolerance because a nearly degenerate unconstrained eigenvector had
  a reflection asymmetry of numerical origin. Enforcing the proved
  even-reflection sector fixed the issue, with the physical eigenvector
  residual retained as a check. This failure and repair are recorded
  in PROGRESS.md.

Run each with `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 <script>`.
The corresponding `*_checks.json` files contain the data. The older
code/spectral-parent checks were not repeated, since their implementations
were unchanged. The new Python files also passed compilation.

Next steps, ranked:

1. Construct a high-channel core with an explicit repeated microscopic
   bond and controlled endpoint maps. The new vacancy theorem can then
   handle two lengths if the core is frustration free. A canonical
   Choi projector alone does not satisfy this hypothesis, and the
   exact cat code cannot be an isolated nearest-neighbor ground space.
2. Select the occupied-number sectors, not just the total physical
   cell length, in the combined selector/core model. Uniform bulk
   frustration-free cutoffs are bounded by Section 18; boundary
   computation or frustrated selection is needed to evade that
   restriction. The finite-hopping vacancy gap, which can be a very
   small tunneling splitting, must control the inter-cell scale.
3. Establish volume-uniform effective interactions and the infrared
   stability of the coupled length/payload chain. Length-dependent
   endpoint strengths, other higher-order terms, and the correlated
   ring twists must all be included. The exact t=u selector does not
   prove the long-cell small-u phase is conformal.
4. Independently pursue the conditional rank-one physical-site model
   in Section 6, which bypasses microscopic block and ring-length
   issues but still needs thermodynamic reduction and infrared
   stability. Establishing it would challenge c<=r(h), not by itself
   settle c<=d-1 or any of the three conclusive outcomes.
5. For a full-class upper bound, obtain a restriction on neutral
   stress-tensor sectors or on the conformal onset scale that survives
   arbitrarily small gaps. The finite-frustration-free-witness theorem
   supplies neither, and should not be used as if it did.




## 20. Boundary Schmidt rank and local-parent obstructions

Status: proved finite-dimensional statements. These exclude a class of
first-order block constructions, not the full c-d conjecture.

### 1. A single low-Schmidt-rank code vector limits all channels

Let V:C^D -> C^q tensor E be an isometry, P=VV^dagger, and let F_a be
Hermitian operators on C^q. Suppose

    V^dagger (F_a tensor I) V = alpha_a Gamma_a,
    alpha_a real and != 0,  {Gamma_a,Gamma_b}=2 delta_ab I.          (1)

If the code contains a normalized vector psi of Schmidt rank s across
this boundary cut, then

    m <= 2 q s - s^2.                                    (2)

The scales alpha_a may differ and may be arbitrarily small.

Proof. Regard the physical Hilbert space as a real inner-product space
with (u,v)_R=Re<u,v>. Write phi=V^dagger psi. The m vectors

    v_a=P(F_a tensor I)psi/alpha_a=V Gamma_a phi

are orthonormal in this real inner product, because

    Re<Gamma_a phi,Gamma_b phi>
      = <phi|{Gamma_a,Gamma_b}/2|phi>=delta_ab.

They belong to the real image of the linear map

    Herm(q) -> range(P),  F -> P(F tensor I)psi.

Projection cannot increase its real dimension. Before projection, its
kernel consists of Hermitian F annihilating the s-dimensional Schmidt
support S in C^q. In an orthogonal decomposition S direct-sum S^perp,
only the (q-s) by (q-s) bottom block of F may be nonzero. The kernel
therefore has real dimension (q-s)^2. The domain has dimension q^2,
so the image dimension is q^2-(q-s)^2=2qs-s^2. This proves (2).

In particular, a code containing even ONE vector that factors across
the active port satisfies m<=2q-1. For an even number m of Majorana
channels, m<=2q-2, and its first-order target has

    c_target=m/2 <= q-1.                                 (3)

This conclusion does not assume a local parent, frustration freeness,
a gap, or exact invariance of the code under the boundary bond. It
applies to scaled compressions, which are the setting of the proposed
loophole. It is only a bound for the specified Clifford target.

For a proposed m-channel code every nonzero code vector must have

    Schmidt rank >= ceil(q-sqrt(q^2-m))                   (4)

at EACH port that realizes the channels. The integer ceiling is taken
on the exact expression. In particular m>2q-1 requires a completely
entangled code subspace across each boundary cut.

If there is also a Hermitian physical lift F_0 of the logical parity,
anticommuting with every Gamma_a and squaring to I after normalization,
apply (2) to the m+1 generators. For s=1 this still gives m<=2q-2.

### 2. Uniform positive bonds of rank less than q

Let h>=0 act on C^q tensor C^q and have matrix rank r<q. Let G_N be the
COMPLETE kernel of the open-chain sum of h, without boundary penalties.
For every unit vector a at the first site, G_N contains a product state
with that first-site vector. The same statement holds at the last site.

Proof. Choose r vectors w_mu spanning range(h). For a fixed nonzero a,
the conditions

    <w_mu|a tensor b>=0,  mu=1,...,r,

are r homogeneous complex linear equations for q entries of b. There
is a nonzero solution b. Repeat this step to choose the next site. Every
bond of the resulting normalized product state is in ker(h). For a
specified last site, repeat the argument in reverse. QED.

Consequently (3) excludes a first-order Majorana violation from every
such complete ground code. There is a stronger statement.

Suppose a Hermitian F has compressed action A on G_N and A^2=alpha^2 I
with alpha>0. Every eigenvector of F can be extended to a ground product
state, so the maximum and minimum eigenvalues of F occur as Rayleigh
quotients of A. Hence ||F||<=alpha. Conversely ||A||<=||F||, so
||F||=alpha and F^2<=alpha^2 I. For psi in G_N,

    alpha^2 ||psi||^2 = ||P F psi||^2
       <= ||F psi||^2 <= alpha^2 ||psi||^2.

Equality forces (I-P)F psi=0. Thus every such F preserves the code.
Because every boundary vector occurs in a code product state, a boundary
operator vanishing on the whole code is the zero physical operator.
The Clifford identities of A_a therefore imply the same identities for
F_a/alpha_a themselves on C^q. For even m this implies

    2^(m/2) divides q,   m/2 <= nu_2(q) <= log_2 q.        (5)

Here the representation fact used is that the complex even Clifford
algebra is M_(2^(m/2)), whose finite-dimensional unital representations
are direct sums of its defining representation.

Furthermore both endpoint lifts preserve G_N. Since the physical
endpoint operators commute when N>=2, their logical restrictions
commute. Thus this subclass cannot realize even one pair of the desired
joint maps Gamma_a and i Gamma_a P_logical, whose commutator is nonzero.
This conclusion does not assume that the eventual inter-cell Hamiltonian
preserves the code exactly. It is forced already by the low-rank parent
and the exact involution compressions.

These claims concern the entire ground kernel. Selecting a smaller
subcode by new terms can destroy the product-extension property. Extra
boundary penalties, higher-rank bulk constraints, frustrated parents,
and approximate rather than exact Clifford maps are not excluded.

### 3. A quantitative version for approximate Clifford relations

Let A_a=V^dagger F_a V/alpha_a be Hermitian. For any unit code vector phi
of boundary Schmidt rank s, form its real Gram matrix

    M_ab=Re<A_a phi,A_b phi>.

It is positive semidefinite and rank(M)<=R=2qs-s^2. If
||{A_a,A_b}/2-delta_ab I||<=epsilon entrywise, then
|M_ab-delta_ab|<=epsilon and ||M-I||_op<=m epsilon.
For m>R, M has a zero eigenvalue, so

    epsilon >= 1/m.                                      (6)

A sharper, readily tested formulation is ||M-I||_op>=1.
Thus arbitrarily accurate normalized Clifford relations cannot bypass
the rank obstruction at fixed m. This does not give a lower bound on
unnormalized errors independent of the small alpha_a.

For the low-rank-parent setting, if ||A^2-alpha^2 I||<=eta for the
UNNORMALIZED compression A=V^dagger F V, extension of extremal physical
eigenvectors gives ||F||^2<=alpha^2+eta. Since A^2>=alpha^2-eta,

    ||(I-P) F V||^2 = ||V^dagger F^2 V - A^2|| <= 2 eta.   (7)

For disjoint ports F,B, their normalized compressions obey
||[A_F,A_B]||<=2 ||Q F V|| ||Q B V|| after using [F,B]=0.
If their square defects are eta_F, eta_B, this gives
||[A_F,A_B]||<=4 sqrt(eta_F eta_B). This is an unnormalized
inequality; the target commutator has size 2|alpha_F alpha_B|.
Therefore a parity-dressed pair requires

    sqrt(eta_F eta_B) >= |alpha_F alpha_B|/2.              (8)

The required relative error is of order one. The low-rank FF-parent
obstruction is robust against improving the normalized accuracy.

### 4. Two additional local-parent classes

A commuting-projector parent cannot supply the joint noncommuting ports.
Let P be the product of its commuting ground-space projectors. Assume
no individual projector touches both active endpoints. For endpoint
operators F and B, partition the projectors into those touching only
the first endpoint, only the last endpoint, and neither. The last
family commutes with both operators; the first family commutes with B;
the second commutes with F. These commutations and idempotence give

    P F P B P = P F B P = P B P F P.                      (9)

Hence the compressed endpoint operators commute, ruling out the
parity-dressed Clifford pair. The argument does not require F or B to
preserve P, or exact invariance under the later inter-cell perturbation.
For a nearest-neighbor commuting-projector parent, three or more sites
ensure that no parent term touches both endpoints. This statement also
covers commuting positive bonds after replacing them by kernel
projectors. It does not cover a noncommuting frustration-free parent.

A separate dimension obstruction applies to a standard injective MPS
code with a nonzero nearest-neighbor frustration-free parent. State
its assumptions explicitly. The virtual bond dimension is D_v, the
two-site MPS map from M_(D_v) to C^q tensor C^q is injective, its image
is annihilated by a nonzero positive two-site parent term, and the
complete open-chain code for the cell is the D_v^2-dimensional MPS
edge space. Then

    D_v^2 < q^2,  hence D_v<=q-1.                         (10)

The strict inequality follows because otherwise the injective image
would fill the two-site physical space and no nonzero positive term
could annihilate it. An even m=2n Clifford family on the complete edge
space requires 2^n to divide D_v^2, so

    n <= 2 nu_2(D_v) <= q-1.                             (11)

For the last inequality put k=nu_2(D_v). If k=0 it is immediate. If
k>=1, 2k<=2^k<=D_v<=q-1. Thus this MPS-parent subclass cannot provide
c_target=n>q-1, even without inspecting its endpoint maps. The result
uses the complete edge-space dimension; a differently selected subcode
or a parent with larger interaction range is not covered by (10).

These exclusions identify features needed by the unresolved
frustration-free route. A useful complete code must have no product
vector at either active boundary, use noncommuting constraints, avoid
the low-rank product-extension setting, and go beyond the specified
injective two-site MPS parent class. None of these requirements proves
that such a code cannot exist.


### 5. Consistency with the earlier joint Choi code

The earlier unrestricted dilation is not excluded by the new theorem.
For a unit logical state phi, its left-port reduced state is

    rho_L=(I+alpha sum_a <Gamma_a>_phi F_a)/q.

The Clifford expectation vector obeys sum_a <Gamma_a>^2<=1. Indeed,
for real u, (sum_a u_a Gamma_a)^2=||u||^2 I; taking u equal to the
expectation vector proves the claim. Since tr(F_a F_b)=q delta_ab,

    ||sum_a <Gamma_a> F_a|| <= ||sum_a <Gamma_a> F_a||_HS
                             <=sqrt(q).

Therefore every eigenvalue of rho_L lies between
(1-alpha sqrt(q))/q and (1+alpha sqrt(q))/q. The same proof applies
to the right-port Clifford family. At the earlier
alpha=1/(4m sqrt(q)), every nonzero code vector has full Schmidt rank
q at BOTH ports, with entropy at least log(q)-log(1+1/(4m)).
This is an analytic check on the scope of (2), not a new numerical
claim or a construction of a local parent. For the q=9,m=80 resource
example, (4) requires rank at least 8 at each port for every code
vector; the Choi code has rank 9.


## 21. The sharp qubit frustration-free witness

Status: proved. This concerns a positive uniform bond on an OPEN chain
with no boundary terms. It does not classify conformal frustrated chains.

The finite-witness integer from `finite_ff_witness.md` has the sharp value

    K_FF(2)=4.                                           (1)

Thus a positive two-qubit bond that is frustration free on four open
sites is frustration free at every length. A finite longest FF segment
can have two or three sites, but cannot be longer than three.

### 1. Reduce to a two-dimensional allowed space

Put S=ker(h) subset C^2 tensor C^2. The zero space of H_N is

    G_N = intersection_j C^2 tensor ... tensor S_(j,j+1)
                              tensor ... tensor C^2.     (2)

Only S matters, not the nonzero eigenvalues of h.

If dim(S)=0, there is no zero vector at N=2. If dim(S)>=3,
rank(h)<=1<2, so the product-extension proof in
`boundary_rank_obstruction.md` supplies a zero vector at every length.
Alternatively a single homogeneous quadratic equation <w|a,a>=0 has a
nonzero solution over C, giving the repeated product a^tensor N.

If dim(S)=1, write S=span(v). An extension to three sites requires
v_12 tensor x_3=y_1 tensor v_23 for nonzero vectors x,y. Taking Schmidt
ranks across site 1 versus sites 2,3 forces v to have rank one. Thus
v=a tensor b. The same equation then forces a proportional to b. This
case has the repeated product a^tensor N at all lengths. Otherwise
G_3=0.

It remains to examine dim(S)=2.

### 2. The product vectors in a two-dimensional two-qubit space

Identify a vector with its 2 by 2 coefficient matrix. Product vectors
are the zeros of the determinant. On a projective line P(S), determinant
is a homogeneous quadratic polynomial. There are three cases over C.

(a) It vanishes identically. Then S is a line entirely of product
vectors, necessarily a tensor C^2 or C^2 tensor b. To verify the last
statement, take two independent products a tensor b and c tensor d.
If both a,c and b,d were independent, their sum would have matrix rank
two. Hence one pair is dependent. The subspace contains a tensor a or
b tensor b and is FF at all lengths.

(b) There are two distinct product directions. Write

    S=span(a_0 tensor b_0, a_1 tensor b_1).               (3)

Both a_0,a_1 and b_0,b_1 are bases, since otherwise case (a) applies.
For a three-site vector in S_12 tensor C^2, expand the third site in
the b basis:

    psi=sum_(i,j) c_ij a_i tensor b_i tensor b_j.

Membership in C^2 tensor S_23 says that the slice with last vector b_j
has middle factor a_j. Independence of the first-site a_i then gives

    c_ij !=0 only if b_i is proportional to a_j.          (4)

Define a directed graph on the two bond labels, with edge i->j if
b_i is proportional to a_j. Each vertex has at most one incoming and
one outgoing edge. Equation (4), iterated along the chain, says that
G_N is spanned by the product vectors associated to graph paths of
N-1 bond labels. The iteration is valid because distinct paths have
distinct first labels and hence independent first-site factors.

If a directed cycle exists, following it gives a product ground state
at every open length. If no cycle exists, a path visits at most the two
vertices. Therefore G_N=0 for N>=4. In particular G_4 nonzero forces
FF at every length.

(c) There is one double product direction. The line is tangent to the
rank-one quadric. It can be written

    S=span(a tensor b, a tensor d+c tensor b),            (5)

where a,c and b,d are bases. To see the normal form, transform the
unique product vector to |00>. For a second matrix M, the determinant
of x|00>+y M has a double root at y=0 only if M_11=0. Subtract its
00 entry, leaving the nonzero 01 and 10 entries in (5). Neither can
vanish, because the line is not all-product.

If b is proportional to a, the repeated product a^tensor N works at
all lengths. Otherwise any element of S_12 tensor C^2 is

    psi=a tensor b tensor x
         +(a tensor d+c tensor b) tensor y.

Its slice with first vector c is b tensor y. If psi also satisfies
the 23 constraint, this product belongs to S. The only nonzero product
direction in S is a tensor b, and b is not proportional to a, so y=0.
The first-site a slice now forces b tensor x in S, giving x=0. Thus
G_3=0.

Every case proves the upper bound K_FF(2)<=4.

### 3. Sharpness and a general lower bound

Let |+>=(|0>+|1>)/sqrt(2), and use the explicit rank-two projector

    h=I-|0,1><0,1|-|1,+><1,+|.                           (6)

The two subtracted vectors are orthonormal. Its allowed-space graph
is the single path from |0,1> to |1,+>. The three-site ground space is
span(|0,1,+>), but G_4=0. Hence K_FF(2)>3, proving (1).

The same idea works at any q>=2. Choose v=q^(-1/2) sum_j |j>, and let
S be the span of the q orthonormal products

    |0,1>, |1,2>, ..., |q-2,q-1>, |q-1,v>.

Put h=I-P_S. The left factors are a basis, and the right factors
|1>,...,|q-1>,v are another basis. The graph proof of case (b) applies
unchanged. Its only maximal path has q bond labels, so

    G_(q+1)=span(|0,1,...,q-1,v>),  G_(q+2)=0.

Consequently the general finite-witness length obeys

    K_FF(q)>=q+2, q>=2.                                  (7)

The nonorthogonal final local letter matters. A directed graph built
only from an orthonormal physical alphabet has a smaller acyclic path
bound. There is no contradiction, since S here is a quantum subspace.
No matching upper bound for q>=3 is proved.

### 4. Scope for the block construction

This makes the earlier noetherian obstruction explicit for qubits. A
longest-FF-segment selector using only a uniform positive qubit bond
can select at most three occupied sites. Boundary penalties and
frustrated kinetic selectors remain outside the theorem. Even within
this subclass, (1) is a statement about exact zero modes, not about
central charge or conformal onset lengths of arbitrary chains.


## 22. A uniform graded-swap core with controlled ports

Status: exact finite-dimensional construction and port formulas, proved
below. Its large-cell Clifford approximation is controlled. No finite-
coupling CFT or counterexample is asserted.

### 1. Repeated bond and complete ground code

Take q=n+1 local states, one even vacuum |0> and n odd colors |a>,
a=1,...,n. Define the graded swap and the positive bond

    P_g|a,b>=(-1)^(p(a)p(b))|b,a>,   h0=(I-P_g)/2,
    p(0)=0, p(a)=1 for a>0.                               (1)

On an open chain of N>=n sites, the complete zero space of the repeated
h0 has dimension 2^n. Label its orthonormal basis by subsets S of colors.
For r=|S|, its vector |S;N> is the signed equal superposition of every
configuration with precisely these r colors and N-r zeros. The sign is
the parity of the spatial permutation of the occupied colors relative
to their increasing order. Every amplitude has magnitude

    [N!/(N-r)!]^(-1/2).                                   (2)

Proof of completeness. Each adjacent graded swap must have eigenvalue
+1 in a zero vector. These swaps move any configuration with distinct
colors into any other configuration with the same color set. They fix
all relative coefficients to (2). A repeated odd color can be brought
to adjacent equal letters, where the graded swap is -1. Its coefficient
must vanish. There are no further constraints on the distinct-color
sets, and (2) explicitly satisfies every swap. Thus there is one ground
vector per subset. QED.

This is a genuinely uniform, nearest-neighbor, positive parent. Its
rank count includes n forbidden identical-color pairs, n vacuum-color
antisymmetric pairs, and n(n-1)/2 color-color symmetric pairs. Hence

    rank(h0)=n(n+3)/2.                                   (3)

### 2. Both endpoint maps

Identify the logical subset basis with the Fock space of n fermions.
Let f_a remove color a with the sign (-1)^(number of smaller colors),
let R=sum_a f_a^dagger f_a, and P=(-1)^R. Define

    g_N(R)=sqrt(N-R)/N,
    a_a=g_N(R) f_a,
    A_(a,x)=a_a+a_a^dagger,
    A_(a,y)=-i(a_a-a_a^dagger).                            (4)

The physical port operators are

    X_a=|0><a|+|a><0|,
    Y_a=-i|0><a|+i|a><0|.                                (5)

Their compressions at the first site are A_(a,x), A_(a,y).
To check the coefficient, remove a from a state with r occupied colors.
The surviving configurations with a fixed at site 1 number
(N-1)!/(N-r)!. Dividing by both normalizations from (2) gives

    [(N-1)!/(N-r)!] /
       sqrt([N!/(N-r)!][N!/(N-r+1)!])
      = sqrt(N-r+1)/N.                                  (6)

The fermionic sign is the number of colors smaller than a. This is
just the action of g_N(R)f_a, since R is evaluated after removal.
Creation follows by Hermitian conjugation.

At the last site the removal sign counts the colors greater than a,
so it differs by (-1)^(r-1). The two last-site compressions are

    B_(a,x)=-i A_(a,y) P,
    B_(a,y)= i A_(a,x) P.                                (7)

Consequently the physical inter-cell bond

    v=sum_a [Y_(a,right) tensor X_(a,left)
             -X_(a,right) tensor Y_(a,left)]              (8)

projects to sum_mu i A_mu P tensor A_mu, with the intended common
fermion-parity dressing. These are exact endpoint formulas, but for
n>1 the A_mu are not exact scaled Clifford matrices.

### 3. Controlled approximation and why it cannot violate c-d

Let Gamma_(a,x)=f_a+f_a^dagger and
Gamma_(a,y)=-i(f_a-f_a^dagger). For N>=n,

    ||sqrt(N) A_mu-Gamma_mu||
       = 1-sqrt(1-(n-1)/N).                              (9)

Each A_mu acts in independent 2 by 2 blocks joining a subset not
containing a to that subset with a added. The normalized amplitude in
such a block is sqrt(1-r_other/N), where 0<=r_other<=n-1. This proves
the equality in (9). For n=1 the compressions are exact at every N.

Thus this core realizes 2n approximate Clifford channels, with error
O(n/N), and both ports have the desired parity relation. But its free
first-order target has only c_target=n=q-1, before paying any selector
cost. The code dimension 2^n itself rules out more than 2n even Clifford
generators. It cannot access the proposed window above q-1.

Allowing p>=1 even colors and n odd colors in the same graded-swap
construction makes the code dimension grow with N. It does not repair
the exact-channel budget. The complete ground space contains the
constant product of any one even color. The boundary-Schmidt theorem
then bounds every even exact Clifford family by 2q-2 channels, with
q=p+n. Its normalized Gram obstruction also prevents arbitrarily
accurate larger families at fixed channel count on that entire code.

Removing all even colors instead gives the exterior power of C^q.
It has no zero states beyond N=q and so does not make arbitrarily long
cores. Selecting a proper subspace of any of these codes by new terms
is a different construction; neither this exclusion nor the exact
endpoint formulas establish what those terms can achieve.

### 4. Checks and the remaining issue

`graded_core_checks.py` constructs the physical code, repeated bond,
and endpoint matrices at (n,N)=(1,3),(2,3),(2,5),(3,4),(3,5). It verifies
the code residual, both ports, the error formula (9), and the physical
bond rank. At sizes with Hilbert dimension at most 256 it also checks
the complete ground multiplicity by dense diagonalization.

This is a useful local replacement for an abstract code in the vacancy
construction, but it only saturates the unblocked channel budget. It
therefore stalls as a counterexample route. A useful high-channel core
must avoid the product vectors of this family and have a larger logical
space, while retaining a repeated local parent and controlled ports.


## 23. Bond-rank bounds for microscopic current sectors

Status: proved finite-dimensional symmetry statement. The central-charge
corollary assumes the stated CFT current correspondence and Sugawara
completeness. Neutral or emergent sectors remain uncontrolled.

### 1. Canonical connected bond

Remove the scalar and one-site terms from a Hermitian bond and write

    H=sum_i b_i + sum_i v_(i,i+1),
    tr_left(v)=tr_right(v)=0,
    v=sum_(a=1)^r s_a A_a tensor B_a,  s_a>0.              (1)

The A_a and B_a are orthonormal traceless Hermitian operators for the
Hilbert-Schmidt inner product. This is a REAL singular value decomposition
of the Hermitian coefficient matrix. The positive s_a are the connected
operator Schmidt singular values; r is the connected Schmidt rank r(h).
Let n_1,...,n_t be the multiplicities of its distinct singular values.

Suppose a compact connected on-site group G, acting by U(g) on C^d,
preserves b and v under simultaneous conjugation. This holds
automatically if the uniform on-site action preserves the full ring
Hamiltonian for every length at least three. After centering as in (1),
one-site and distinct two-site traceless supports are orthogonal in
Hilbert-Schmidt inner product, so their changes under conjugation
cannot cancel. Let g be its Lie algebra after removing scalar generators.

There is an ideal k of g such that every generator in k commutes with
b and with every A_a and B_a individually. Moreover

    g/k embeds in direct_sum_j so(n_j).                  (2)

In particular

    dim(g/k) <= sum_j n_j(n_j-1)/2 <= r(r-1)/2.            (3)

These are statements about the interaction channels, not about the full
on-site symmetry dimension. The kernel can be large.

### 2. Proof

The real left Schmidt space span_R{A_a} is intrinsic. It is the image
of partial contractions of v on the right. Since simultaneous on-site
conjugation preserves v, conjugation preserves this space, and likewise
the right space. Unitary conjugation preserves the Hilbert-Schmidt
inner product, so its matrices O_L(g),O_R(g) are orthogonal. They lie
in SO(r) for connected G. Invariance of v is

    O_L S O_R^T=S,   S=diag(s_1,...,s_r).                 (4)

Multiplying (4) by its transpose gives

    O_L S^2 O_L^T=S^2.

Thus O_L preserves each eigenspace of S^2, and its infinitesimal action
lies in the direct sum of so(n_j). Since S is positive, O_L commutes
with S itself. Equation (4) then gives O_R=O_L in the paired Schmidt
bases. In particular their infinitesimal kernels coincide.

Let k be this common kernel. A generator in k commutes with every
Schmidt operator on each side. It also commutes with b by hypothesis.
At any site i its local generator Q_i therefore satisfies

    [H,Q_i]=0                                            (5)

on every ring. The kernel is an ideal, and the quotient representation
is faithful by its definition. This proves (2) and (3).

The use of singular-value multiplicities is material. A connected bond
with pairwise distinct singular values has no nontrivial continuous
action on its interaction space, even if there are local conserved
on-site observables. Degeneracies permit the continuous rotations.

### 3. Conditional central-charge and torus bounds

Assume the following infrared input.

(i) The relevant microscopic charge densities have scaling limits
whose nonzero dynamical components supply the independent affine
currents of a unitary CFT, with nonzero finite velocities.

(ii) These currents exhaust the stress tensor, possibly as independent
factors with different velocities. There is no positive-charge coset.

A generator in k cannot supply a nonzero dynamical affine current under
(i). Equation (5) says its density is stationary at each individual
lattice site, so its dynamical spectral measure is supported at omega=0
at every spatial momentum. A finite-velocity affine current has
nonzero-momentum spectral weight at omega=v|k|. This weight cannot be
the scaling limit of that identically stationary density. Removing
these inactive generators, the joint current algebra is an image of
g/k and also an image of the microscopic algebra in pu(d).

For a compact reductive current algebra l, the supplied Sugawara input is

    c_Sug=rank(abelian part)
           +sum_simple k_a dim(l_a)/(k_a+h_a^vee)
          <=dim(l).                                     (6)

The inequality is strict on any nonzero simple factor. Assumption (ii)
sets c=c_Sug. Therefore

    c <= min[d^2-1, sum_j n_j(n_j-1)/2].                  (7)

No affine-level bound is used. For r<=3 this proves c<=r within this
current-complete class. For r=1 it gives c=0, so any critical rank-one
chain with positive c must have neutral or emergent sectors outside
these hypotheses. This is consistent with the conditional rank-one
construction in RESULT.md Section 6, whose propagating U(1)s are
emergent rather than exact microscopic charge-transfer symmetries.

There is a stronger conclusion if the infrared is Gaussian and all its
bosons are resolved by commuting microscopic on-site charges. The
abelian image in (2) has dimension at most the rank of its compact
ambient Lie algebra. Commuting real skew-symmetric matrices can be
simultaneously put into independent 2-plane rotation blocks, so

    c <= min[d-1, sum_j floor(n_j/2)] <= floor(r/2).       (8)

The identification c=number of independent resolved currents is the
same Gaussian-current completeness hypothesis used in C5/C7. Equation
(8) controls these modes by the actual bond rank, and holds regardless
of perturbative order or cell size when applied at the microscopic
site level. A global on-site symmetry alone is insufficient, because
it need not resolve every infrared sector.

### 4. Examples and limits

For an XX qubit bond, its two equal singular values give n_1=2 and
c<=1 under current completeness. For XXZ with three nonzero channels
and unequal transverse and longitudinal strengths in magnitude, the
multiplicities are 2 and 1, again giving c<=1. This includes arbitrary
levels in principle; the active algebra is abelian.

For an isotropic Heisenberg bond, all three singular values coincide,
so (7) gives c<3 for an SU(2) current-complete theory. It does not derive
the sharper spin-1/2 bound c<=1 or an affine-level ceiling.

For two decoupled XX chains in one d=4 site with unequal strengths,
r=4 and the multiplicities are 2,2. Both (7) and (8) give c<=2, matching
the two charge modes. At equal strengths the dimension bound becomes
6, while the torus bound remains 2.

An Ising critical chain is not a contradiction at rank one. Its c=1/2
sector has no microscopic continuous affine currents. A neutral coset
or emergent currents remain the missing part of any general rank bound.
No conclusion about those sectors is drawn here.

`schmidt_symmetry_checks.py` solves the full on-site infinitesimal
symmetry equations for explicit bonds, computes their induced action
on the actual Schmidt space, and checks dimensions and singular-value
multiplicities. These are algebraic checks, not measurements of c.


## 24. Fourth-run verification, limitations, and next steps

This run remains nonconclusive. It followed three distinct directions
beyond their first partial results: local-code obstructions, an explicit
higher-rank frustration-free construction, and microscopic current
bounds. The constructive graded-swap route stalls at c_target=q-1;
the other results exclude restricted mechanisms or require current
completeness. None controls every neutral infrared sector.

Reproduce the new checks from this folder with

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 boundary_ff_checks.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 graded_core_checks.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 schmidt_symmetry_checks.py
python3 -m py_compile boundary_ff_checks.py graded_core_checks.py schmidt_symmetry_checks.py
```

Each numerical script writes its correspondingly named JSON file.
No previous C8 numerical suite was rerun merely for orientation.

The boundary checks use q=2,3,4, N=4 low-rank parents. The largest
compressed-extremum error is 2.67e-15 and the product-state residual is
below 2.77e-16. The real-image rank formula is checked for every s at
q=2,3,4,5. Six types of two-dimensional qubit allowed spaces, each in
three common invertible physical bases, have the predicted ground
multiplicities at N=3,4,5. The sharp-witness family is checked for
q=2,3,4 up through dimension 4096. At the first frustrated lengths
N=4,5,6 its lowest energies are respectively 0.0761204675,
0.2218619272, and 0.3169872981. These positive finite-size values
support the explicit kernel proof, not a continuum conclusion.

The graded-core checks cover (n,N)=(1,3),(2,3),(2,5),(3,4),(3,5), up
to physical dimension 1024. The endpoint formula errors are below
2.36e-16, code isometry errors below 1.40e-15, and every reported
normalized Clifford error agrees with the exact formula. Complete
spectra through dimension 256 confirm the 2^n ground multiplicity.
The symmetry checks cover seven models, including the critical Ising
algebra as a neutral-sector control. Their largest Schmidt-space
invariance residual is below 3.93e-16. All three scripts compiled.

The proofs were audited for the hypotheses that matter: complete
versus selected ground spaces, open versus periodic chains, real
versus complex dimensions, exact versus approximate Clifford
relations, and microscopic versus emergent currents. The earlier
Choi code has full boundary Schmidt rank for every code vector, so
the new obstruction does not contradict its existence. No independent
referee or subagent check was performed. No network, cluster job,
commit, push, or change outside this folder was made.

The next steps, ranked, are:

1. Construct a high-channel repeated-bond core that avoids the new
   exclusions. For an exact first-order family with m channels, every
   code vector needs boundary Schmidt rank at least
   ceil(q-sqrt(q^2-m)) at both ports. A noncommuting, higher-rank FF
   history code with boundary initialization remains possible, but
   no explicit uniform implementation with its illegal-sector gap
   has been supplied. Approximate codes must beat the normalized
   Gram obstruction, not only an absolute error tolerance.
2. If such a core is found, select the occupied-number sectors in
   both mobile cell lengths and control the resulting coupled
   length/payload theory at finite coupling. The graded-swap code
   and chemical-potential-only vacancy model do not solve this.
3. Separately, establish the two thermodynamic and infrared premises
   for the explicit rank-one, d=20992 construction of Section 6.
   The new symmetry bound explains why an emergent-current argument
   is necessary there; it supplies no stability theorem.
4. For an upper bound in the full class, control neutral or emergent
   stress-tensor sectors, or establish uniform conformal visibility.
   The Schmidt-multiplicity bound sharpens the existing microscopic
   current budget without addressing this remaining sector.
5. Determine K_FF(q) for q>=3 if pursuing longest-FF-segment selectors.
   The q+2 lower bound is explicit, but neither its sharpness nor a
   useful upper bound at higher q has been established. This is a
   restricted algebraic question, not a substitute for a CFT bound.
