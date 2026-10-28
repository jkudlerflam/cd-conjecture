# C24. Independent review of the species-graph candidate

[heuristic, referee verdict] I would not yet accept an unqualified claim that this is an established refutation of c <= d - 1. I would accept the graph argument, the sparse port construction, the charge-preserving circuit, and the conditional Gaussian stability argument. I found no fatal algebraic error, unavoidable relevant resonance, or contradiction in the fixed-coupling hierarchy. The proposal is substantially stronger than a large-code-space counting argument.

[assumed, outstanding inputs] Two steps still carry the conclusion. The complete autonomous sparse-encoding controller must satisfy the stated initialization and readout contract, and the mobile microscopic model must admit the claimed local low-frequency elimination with corrections o(s) relative to its first-order Hamiltonian. The notes give a plausible construction argument for the former and a physical expansion argument for the latter. The supplied executable checks establish neither complete statement. These are serious verification gaps, rather than demonstrated failures. I explain below what would close them without demanding a rigorous lattice-to-CFT theorem.

[proved, scope] All new files and all reproduced outputs are inside work_C24. No network was used and no source outside this folder was modified. The project STATE.md was read but not updated because the task expressly prohibits writes outside this folder. PROGRESS.md is the local working memory. The reviewed source tree changed during the review. The requested contribution files and C8/C13/C20 interfaces are the basis of this report; the later MICROSCOPIC_HAMILTONIAN.md was also inspected. Later controller prototypes mentioned by the changing STATUS.md have not been audited here. source_manifest.json records the source hashes at the end of the reading phase, not an immutable snapshot of their earlier contents.

[proved, terminology] A claim tagged “proved” below is either finite algebra or a deduction under explicitly named assumptions. “Numerical” refers to the stated finite checks, not a proof for arbitrary size. “Heuristic” denotes a physical inference or referee assessment. “Assumed” identifies an input not established by this review. In particular, bosonization and the standard identification of stable Gaussian modes with compact bosons are assumed continuum inputs.

## 1. Selection rule and the six-field threshold

### 1.1 The graph and charge lattice

[proved] Put s = p^2 + p + 1. The one-dimensional subspaces of F_p^3 number (p^3 - 1)/(p - 1) = s. The same count applies to lines, represented by nonzero covectors modulo rescaling. A line contains (p^2 - 1)/(p - 1) = p + 1 points. Two distinct points span one two-dimensional subspace, hence determine one line. The bipartite incidence graph therefore has

\[
q=2s,\qquad N=(p+1)s,
\]

and no four-cycle. Three noncollinear points and their joining lines give a six-cycle. The graph is connected and has girth six. Its oriented incidence matrix B has rank q - 1 and integer cycle rank N - q + 1. I use B for this matrix to distinguish it from the encoding isometry W.

[proved] Give creation on edge e = (alpha,beta) charge w_e = e_alpha - e_beta, with alpha in the upper part. For a normal-ordered monomial let z_e be creations minus annihilations of that flavor. Its species charge is Bz. If Bz = 0 and z is nonzero, reverse edges with negative coefficients and regard their absolute coefficients as an integer flow. Conservation prevents a directed path from ending. A finite nonempty support therefore contains a directed cycle. Removing loops from the path gives a simple cycle containing at least six distinct edges. Consequently

\[
z\ne0,\quad Bz=0\quad\Longrightarrow\quad \|z\|_1\ge6.
\]

[proved] A monomial containing r fermion fields has ||z||_1 <= r. Thus every neutral monomial with at most four fields preserves every flavor. In fact the implication holds through five fields, and total upper-species U(1) already forbids an odd number of fields in a neutral monomial. This is a selection rule for operators. It does not promote the microscopic symmetry to U(1)^N. A six-cycle conversion changes three occupations into three different occupations while preserving all q species charges.

[numerical, p = 2,3,5,7] The independent generator in independent_checks.py verifies the incidence identity A A^T = p I + 11^T and checks uniqueness of every one-edge and two-edge charge multiset. The two-edge tests cover 231, 1,378, 17,391 and 104,196 multisets, respectively. A neutral four-field conversion would require two different multisets with the same charge. None occurs. The copied torus_bottleneck script additionally enumerates all nonzero integer vectors of l1 norm at most four for p = 2.

### 1.2 Scaling dimensions, including actual resonant six-field operators

[assumed] Use integer chiral creation counts u,v and the standard free-fermion bosonization normalization. Define z = u + v and m = u - v. Then m = z modulo two coordinatewise, the momentum is m dot k for centered bands, and

\[
\Delta_0=\frac{\|m\|_2^2+\|z\|_2^2}{4}.
\]

[proved under this normalization] For an integer pair m_e,z_e of the same parity,

\[
m_e^2+z_e^2\ge2|z_e|.
\]

For even nonzero z_e, z_e^2 >= 2|z_e|. For odd z_e, m_e^2 >= 1 and z_e^2 + 1 >= 2|z_e|. Summing gives Delta_0 >= ||z||_1/2 >= 3 whenever Bz = 0 and z is nonzero. Equality is attained by a primitive six-cycle with six z_e = +/-1 and six m_e = +/-1. Each distinct fermion contributes 1/2.

[proved under this normalization] If z = 0, write m = 2u. Then Delta_0 = ||u||^2. The only nonconstant vertices of dimension at most two are u = +/-e_a and u = +/-e_a +/-e_b with a != b. These are single backscattering and products of two backscatterings. Same-flavor double backscattering has u = 2e_a and dimension four, not two. Its point-split microscopic representative includes the Pauli suppression already encoded in the vertex dimension. One must not count those derivatives twice.

[proved] A literal six-field expression can contain contractions or current factors multiplying a lower-field vertex. Normal ordering separates those cases. They either reduce to the previously classified primaries or add descendant dimension. A zero-vertex six-field current product has dimension three. Thus “all other resonances start at dimension three” is a claim about normal-ordered primaries and their descendants, not a count of symbols in an unnormal-ordered lattice expression.

[numerical, p = 2,3] I explicitly searched every six-cycle and all 64 choices of m_e = +/-1 on each cycle, with one fixed orientation of z per cycle. At p = 2 there are 28 cycles and 18 non-oscillating chiral assignments. At p = 3 there are 234 cycles and 220 assignments. These are counts of assignments in the stated enumeration, not counts of independent cosine interactions. The arithmetic progression really does leave six-field resonances.

[proved, explicit example] In the independent p = 2 edge ordering, one cycle has edge labels

\[
(1,6,5,8,9,3),\quad
z=(1,-1,1,-1,1,-1),\quad
m=(-1,-1,1,1,-1,1).
\]

The representative order and incidence list are fixed by independent_checks.py. Here Bz = 0, sum m = 0 and

\[
-1-6+5+8-9+3=0.
\]

It follows directly that m dot k = 0 for k_a = pi/4 + a pi/[12(N+1)]. This is an allowed six-field operator of free dimension three. Momentum detuning alone does not remove it.

[assumed] The interacting point is not supplied as a numerical stiffness matrix. The candidate supplies an inequality for its full scaling-dimension form. It is therefore possible to give bounds on interacting dimensions, not their individual numerical values. For a density/current-decoupled payload theory,

\[
\Delta(m,z)=\frac14(m^T K m+z^T K^{-1}z).
\]

With length mixing and density-current couplings, the corresponding full metric must be used.

[proved under the stated metric inequality] If Delta >= kappa Delta_ref for the complete joint vertex lattice and kappa > 2/3, every previously unexcluded resonant vertex has Delta >= 3 kappa > 2. The six-cycle example above is therefore irrelevant, including when its momentum is zero. At kappa = 2/3 the lower bound permits marginality and is insufficient. Strictness matters.

[proved, a useful countercheck] Positivity of stiffness alone would not suffice. For the displayed resonant cycle, m dot z = -2 and ||m||^2 = ||z||^2 = 6. Let e_1 = m/sqrt(6), and let e_2 be the normalized component of z perpendicular to e_1. Choose K to be the identity on the orthogonal complement and to have the block

\[
K\big|_{\{e_1,e_2\}}=
\begin{pmatrix}
1/3&-2\sqrt2/3\\
-2\sqrt2/3&3
\end{pmatrix}.
\]

Its determinant is 1/9 and its trace is 10/3, so it is positive. It obeys Km = -z. Hence Delta = (2 + 2)/4 = 1. This gives a relevant resonant six-field operator at a positive Gaussian point far from the proposed neighborhood. It is not a counterexample to the candidate's stronger kappa condition. It shows why that condition cannot be replaced by positivity or by separate one-mode checks.

### 1.3 Payload and mixed resonances

[proved] For the proposed momenta,

\[
\pi/4<k_a<\pi/3,\qquad
\delta_{\rm pay}=\frac{\pi}{6(N+1)}.
\]

Single backscattering has momentum 2k_a in (pi/2,2pi/3). For a != b the opposite-sign pair has |2k_a - 2k_b| >= delta_pay and lies below pi/6 in absolute value. The same-sign pair has momentum 2k_a + 2k_b in (pi,4pi/3), with distance greater than 2pi/3 from 2pi Z. These exhaust the non-Gaussian payload primaries through dimension two after applying Bz = 0. Total upper-species charge forbids pairing exactly, so C13's parity-only pairing obstruction does not apply.

[assumed, physical-coordinate dictionary] Let a_0 = b + nu be mean cell length and y(x) = x/a_0 - phi_0(x)/pi the coarse-grained cell count. For a payload vertex put S = m dot k/(2pi). Expanding the marker comb gives harmonics exp[2pi i(D+S)y], D integer. Their physical momenta and reference dimensions are

\[
P_{\rm phys}=\frac{2\pi(D+S)}{a_0}\pmod{2\pi},\qquad
\Delta_0(D,m,z)=\frac{(D+S)^2}{a_0^2}+\Delta_{\rm pay}.
\]

[proved within this dictionary] A uniform physical resonance requires D + S = h(b + nu), h integer. Eliminating the integer D gives

\[
m\cdot k-2\pi h\nu=0\pmod{2\pi},\qquad
\Delta_{0,\rm res}=h^2+\Delta_{\rm pay}.
\]

This reproduces C13's formula and identifies why an independent length fermion with an arbitrarily chosen compactification is not a valid substitute. In particular the small physical marker parameter K_marker = 1/(b+nu)^2 does not make a resonant high-harmonic umklapp relevant.

[proved within this dictionary] With nu in [5/12,7/12], 2pi nu lies in [5pi/6,7pi/6]. Thus all mixed dimension-two candidates obey

\[
\operatorname{dist}(\pm2k_a\pm2\pi\nu,2\pi\mathbb Z)\ge\pi/6.
\]

For differences, the nearest endpoint is 5pi/6 - 2pi/3 = pi/6. For sums, the range is contained in (4pi/3,11pi/6), whose nearest possible approach to 2pi is pi/6. Conjugation covers the remaining signs. Pure h = +/-1 length resonance requires nu integer, excluded here. A nonzero h with |h| >= 2 contributes at least four to the dimension.

[proved within this dictionary] Mixed non-oscillating operators of higher dimension are allowed. For example h = 1 and two same-sign payload backscatterings resonate when nu = (k_a+k_b)/pi. Some such values lie inside the chosen density interval. Their free dimension is three. At nu = 1/2 the pure h = 2 length umklapp is resonant and has dimension four. A six-cycle dressed by h != 0 has dimension at least 3 + h^2. The kappa bound therefore excludes relevance or marginality in all these cases.

[proved, scale qualification] The quoted pi/[6(N+1)] and pi/6 margins are cell-coordinate mismatch numerators. Near a specified physical reciprocal vector the actual mismatch is

\[
P_{\rm phys}-2\pi h
=\frac{m\cdot k-2\pi h\nu+2\pi(D-hb)}{b+\nu}.
\]

The physical momentum margin is smaller by b + nu. At fixed b this is positive. An effective cutoff of order one cell compensates this conversion in dimensionless RG units. Either convention is legitimate, but the cutoff and momentum units must be stated together. At fixed p,b only finitely many D,m,z can have dimension near two, so small changes of fillings and the metric can be controlled on a finite list. No all-order Diophantine condition is needed.

[proved, attribution of this qualification] C13 already states the physical-momentum conversion and its cancellation against the physical velocity in the mismatch energy. This is a clarification needed when using the shorter continuation summary, not a newly discovered failure of C13 or an unavoidable extra obstruction at large b.

## 2. Finite encoding and charge-preserving preparation

### 2.1 Direct derivation of both ports

[proved] With u upper vertices, define Kraus maps

\[
K_{1,\alpha}=a|g_\alpha\rangle I,\qquad
K_{0,\beta}=b_0\sum_{e=(\alpha,\beta)}
\left(|g_\alpha\rangle f_e^\dagger+|h_e\rangle f_e\right),
\]

where ua^2 + N b_0^2 = 1. Within a fixed lower-port Kraus map all flags are distinct because the graph is simple. Orthogonality of the flags and the CAR give

\[
\sum_p K_p^\dagger K_p
=ua^2I+b_0^2\sum_e(f_ef_e^\dagger+f_e^\dagger f_e)=I.
\]

For F_e = |1,alpha><0,beta|,

\[
V_L^\dagger F_e V_L=K_{1,\alpha}^\dagger K_{0,\beta}
=ab_0 f_e^\dagger.
\]

A missing edge has zero compression. Choosing a^2 = 1/(2u) and b_0^2 = 1/(2N) gives epsilon = 1/(2 sqrt(uN)).

[proved] The residual Heisenberg channel is

\[
\mathcal D_L(O)=ua^2O+b_0^2\sum_e
(f_eOf_e^\dagger+f_e^\dagger Of_e).
\]

For O = f_j^dagger P and e != j, O commutes with both f_e and f_e^dagger, since each factor f_j^dagger and P separately anticommutes with them. The e contribution is O. The j contribution vanishes by nilpotence. Thus D_L(O) = (1-b_0^2)O. For the undressed f_j^dagger, the e != j contributions instead have a minus sign, giving [1-(2N-1)b_0^2]f_j^dagger. Confusing these channels would invalidate the right port.

[proved] Set g_e^dagger = f_e^dagger P and g_e = P f_e. Direct use of P f_e = -f_e P gives the CAR and the same number operators. Build V_R from these generators and compose W = (I tensor V_R)V_L. The isometry of V_R leaves the left compression unchanged, while the channel identity gives

\[
W^\dagger F_e^L W=\epsilon f_e^\dagger,\qquad
W^\dagger F_e^R W=\epsilon(1-1/(2N))f_e^\dagger P.
\]

For cell Jordan-Wigner operators c_{j,e} = (product_{i<j}P_i) f_{j,e}, c_{j,e}^dagger c_{j+1,e} = f_{j,e}^dagger P_j f_{j+1,e}. This checks the required bond sign and dressing. The product visibility is eta = [1-1/(2N)]/(4uN).

[proved] The physical hopping matrix is a direct sum of two-dimensional exchanges |lower beta;upper alpha> <-> |upper alpha;lower beta>, one block per edge. Distinct edges give distinct basis pairs, even when they share one endpoint. Its norm is max_e |t_e|, not sum_e |t_e|. Compensator norms and inverse history weights still enter the full wrapped norm.

[numerical, N = 1,3,4,6] The independent sparse-branch implementation constructs W without using the candidate's Kraus code. It checks normalization, both compressions, and every charge assignment. The largest test has 64 input states and 5,184 nonzero output rows. The maximum entrywise residual is 6.7e-16. These tests extend the supplied dense two-port tests, whose largest Fock space has four modes. Neither test constructs a full projective-plane code, whose smallest case has N = 21.

### 2.2 Charges, grading, and circuit resources

[proved] Put S_species = sum_s e_s and C_0 = 2 S_species. The flag charges

\[
\xi(g_\alpha)=C_0-e_\alpha,\qquad
\xi(h_{\alpha\beta})=C_0+e_\alpha-2e_\beta
\]

are nonnegative integer vectors of total weight 2q - 1. The lower creation branch has total incremental charge e_beta + xi(g_alpha) + (e_alpha-e_beta) = C_0. The annihilation branch gives e_beta + xi(h_e) - (e_alpha-e_beta) = C_0. The upper branch also gives C_0. Each stage therefore preserves charge after including a reference word with two copies of each species.

[proved] Encoding occupation e by lower beta for zero and upper alpha for one uses N sites. Two port stages use two additional blocks of 2q sites, totaling N + 4q. The initial lower-species baseline is R = sum_e e_beta. The complete physical charge is R + 4 S_species + sum_e n_e w_e. Its upper-species component is Q_logical + 4u, so physical upper-species parity equals logical parity. There is no hidden odd offset per cell.

[proved] On each occupation basis input a stage has u + N branches. A port/flag/residual triple identifies the branch and recovers the input occupation, including whether that edge was initially empty or full. For the right stage the additional signs are (-1)^Q for creation and (-1)^(Q-1) for annihilation. Thus the branch map is reversible on its valid domain with the required fermionic signs.

[proved, circuit construction] Input and output species words have the same multiplicities. A reversible Boolean computation can calculate the target word and a stable adjacent-swap plan. Apply this plan to the physical charged word, reverse its neutral computation, and erase the saved input and branch by calculating the inverse from the physical output. This is coherent compute, use, and uncompute, not an irreversible sorting operation on a superposition.

[proved] The local LOAD, controlled EXCHANGE, and UNLOAD sequence acts on D_s,H_{j,c},D_t as a controlled swap of the two charged records. LOAD transfers the existing charge into T_{j,c,s}; EXCHANGE trades s and t between two sites; UNLOAD returns the helper. Each gate is an involution or unitary block on the full local space and preserves every species count. Two independent species occur on two sites, not in one on-site register. Neutral scratch and Boolean gates introduce only a constant number of tracks.

[proved, exact amplitude construction] For a uniform interval [0,M), at each binary prefix depth there is at most one partial subtree. Its nonempty children contain A and B leaves. The rotation with first column (sqrt(A/(A+B)),sqrt(B/(A+B))) assigns each prefix probability equal to its number of completions divided by M. Induction gives a uniform state. At most ceil(log_2 M) exceptional rotations suffice. The two upper/edge branches require O(log q) indexed blocks because N <= q^2/4. A single block-diagonal rotation opcode and SHIFT can select them. The matrices depend on p, which is allowed; their number as opcode families is constant.

[proved, resource accounting conditional on a fixed interpreter] With m_idx = O(log q), the circuit record count is

\[
D(q)=q+3+2(q+1)m_{\rm idx}=O(q\log q).
\]

If a fixed neutral controller has B_neutral tape states and K_head control/opcode states, the moving-head overlay costs (1+4K_head)B_neutral D(q), multiplied by fifteen for the outer clock. Static marks, Q scratch tracks and READY/FAIL multiply by constants. No q^2 on-site factor is forced by the stated operations. The missing numerical value of C is not itself an asymptotic gap.

[assumed, compiler input] Turning the uniform circuit-generating and navigation algorithms into one fixed reversible one-head transition function uses a standard reversible-machine compilation result. I checked the resource use described in the notes: the preparation workspace is polynomial in p, full-length scans can use counters and frontier marks instead of an n log n instruction list, and n = C_0 p^K can dominate the p-dependent workspace. Polynomial runtime then fits within 2^n sweeps. I have not independently constructed the complete resulting transition function. A standard universality citation alone is insufficient, but the supplied algorithms do give substantially more than such a citation.

## 3. Initialized history core and the precise remaining compiler check

[proved, local sweep interface] The fresh/updated and left/right tags make the domains of the two-site head moves orthogonal. A move and its adjoint form the block [[0,U^dagger],[U,0]], which is unitary. During a left-to-right sweep a right head moves once and cannot match another forward rule; a left head moves behind the sweep. An on-site collision on the unique head then updates the controller. The overlay carries only the charge of its current tape record. This verifies the claimed constant-factor interface.

[proved, clock graph] The fifteen-letter legal words split into the counter words, the n initialization words, and outward/return sweep words. Their counts are 4(2^n-1), n and 2n 2^n, respectively. Inserting one full sweep before each counter increment joins them into one path of

\[
C_n=(2n+4)2^n+n-4
\]

vertices and 2^n sweeps. Every rewrite has the same permitted external neighbors on its two sides, so forbidden-pair indicators are conserved. Positive unitary-decorated propagation therefore cannot lower an illegal clock sector below its legality penalty. The nine-letter binary clock alone would not provide this sweep budget; the supplied fifteen-letter remedy is essential.

[numerical, n = 1,2,3,4] The copied sweeping_interface_checks.py exhausts all clock words, including 50,625 at n = 4, and verifies the complete path and preservation of each legality indicator. Its overlapping noncommuting three-qubit sweep test passes. This is a test of the outer clock and arbitrary work-gate decoration, not a test of the full sparse-encoding interpreter.

[proved conditional on the specified initialization] During the initial I sweep, all data gates are identity. Pins can therefore fix classical tracks, neutral scratch, READY/FAIL, absence of interpreter overlays, and one neutral boundary query helper in a common initial gauge. The unique boundary seed creates one head. Immutable endpoint marks are on separate tracks and can be fixed by positive local conditions. These arrangements avoid borrowing an untrusted quantum slot as the verifier's own initialized scratch.

[proved conditional on a correct verifier] The desired predicate admits exactly two species states at each of N designated occupation slots and one prescribed state at every other quantum-data slot. Thus its accepted dimension is 2^N. TYPE must reject malformed H/T or blank records before interpreting species bits. A scan can count violations into ceil(log_2(n+1)) trusted bits, copy the nonzero test to FAIL, and reverse the count. This checks padding as well as the active data. FAIL must be touched only at the final predicate copy and its inverse, so valid histories never incur a transient FAIL penalty from a Boolean decomposition.

[proved, acceptance lemma] Let H_init >= gamma(I-P), with P its initialized history kernel. Suppose A >= 0, ||A|| <= B_A, A P_good = 0 and PAP >= a(P-P_good). For psi perpendicular to P_good, write x = P psi and y = (I-P)psi. Then

\[
\langle H_{\rm init}+A\rangle\ge
\gamma\|y\|^2+
\left(\sqrt a\|x\|-\sqrt{B_A}\|y\|\right)_+^2.
\]

Where the bracket is positive the two-dimensional quadratic form has determinant a gamma and trace a + B_A + gamma. Its smaller eigenvalue is at least determinant/trace. Where it is negative, ||y||^2 >= a/(a+B_A), which gives a stronger bound. Therefore the complete kernel is P_good and

\[
\operatorname{gap}(H_{\rm init}+A)
\ge\frac{a\gamma}{a+B_A+\gamma}.
\]

[proved conditional on the verifier] An all-time FAIL penalty gives a >= eta_FAIL/C_n because every rejected initialized basis input is flagged at least once. The compressed expectation is diagonal before preparation, and preparation preserves FAIL. Thus the bound applies to superpositions too. B_A can be taken as eta_FAIL n. No commutation of all pulled-back pins on malformed controller states is needed. A terminal-only pin would fail when the inverse computation clears FAIL.

[numerical] Independent random PSD tests of this lemma pass in 60 instances of dimensions 8,15,25. An independent 144-dimensional, nine-vertex history toy has a four-dimensional good kernel with all-time FAIL, versus an eight-dimensional kernel with terminal-only FAIL. This tests the mechanism and its failure mode. It does not verify the claimed 2^N kernel of the full compiler.

[proved conditional on the readout schedule] Let U_p prepare and route both ports with READY zero. Set endpoint READY only while the interpreter is present, exclude the head-present subspace from the port observable, perform neutral identity-data travel during readout, clear READY before inverse preparation, and execute the literal inverse word. At every contributing elementary history vertex the same prepared port is exposed. Hence

\[
G_p^\dagger\widehat F_LG_p=p_L W^\dagger F_LW,
\qquad
G_p^\dagger\widehat F_RG_p=p_R W^\dagger F_RW,
\]

where p_L,p_R > 0 if one full cycle fits. They need not be equal or active simultaneously. The actual inverse eliminates accumulated preparation error between cycles. Truncation of a final cycle is harmless under the same per-vertex readout property. An absorbing halt would not be a unitary alternative.

[heuristic, serious verification gap] These are sound conditional lemmas. Their end-to-end instantiation is the remaining finite check. The requested files do not provide an executable complete controller for length/parameter acquisition, every format query, sparse preparation, output routing, READY timing, inverse preparation, and reset, all using one common fixed neutral alphabet. The supplied tests exercise separate components. This does not demonstrate that no such controller exists, and I do not find a resource obstruction to it. It does mean the complete initialized-core claim has not been independently verified at the same level as the port identity.

[heuristic, sufficient repair] A complete symbolic finite controller with a proof of its register invariants and uniform resource bounds would suffice. It need not print a huge dense h matrix. An executable compiler or interpreter would make that proof auditable: start from the actual pinned configuration, check every elementary readout vertex, verify scratch reset, and establish that all other initialized basis formats trigger FAIL. Small representative cases should test superpositions and malformed types, while the proof covers arbitrary p. Merely running a preloaded gate list does not test autonomous initialization or the full local alphabet bound.

## 4. Matching the C8, C13 and C20 interfaces

[proved conditional on the complete core] The core has the required positive uniform form sum h_c + a_c + b_c. Positive on-site penalties can be split between bulk bonds and endpoints without negative terms. Positive terms annihilate its zero-energy code individually. C8 does not require boundary-free cores, reflection symmetry, fixed ports throughout the preparation, or the same code at two different occupied lengths. It requires one complete core at one occupied length k.

[proved] C8's occupied-site helix filter selects that same occupied length k inside physical lengths b = k+1 and b+1 = k+2. The two physical lengths differ by one neutral vacancy, not by rerunning the computation on k+1 occupied records. This distinction resolves a potentially serious compiler/wrapper mismatch. The history length, charge offset and logical code are identical in both cells.

[proved from the displayed wrapper construction] The alphabet count is

\[
d_{\rm occupied}=2d_{\rm core},\qquad
d_{\rm arrow}=3d_{\rm occupied}+3=6d_{\rm core}+3,
\]

followed by four selector states and two further helix qubits. Therefore d_mobile = 16(6d_core+3) = 96d_core+48. The first filter qubit is already included in d_arrow and must not be omitted or counted again.

[proved] The arrow-wrapper vacancy exchanges transport an occupied internal record unchanged. The lifted core terms preserve total charge. Diagonal tag/legality projectors and every helix label are neutral. Marker conversion changes H_0 into H_2 on a neutral vacancy. The endpoint compensators are neutral matrices, so the lifted port retains w_e and a right-left hopping link is neutral. Every term preserves upper-species parity as well as the full species torus. The extra graded route-A gadget is not needed to supply this symmetry.

[proved] The charged reference offset per cell remains 4u. Thus the physical uniform upper-species field restricts to sum_j Q_j + 4u M. Subtracting 4u times the same coefficient in the marker field supplies the common logical density term. The extra vacancy is neutral, so no length-dependent offset appears. The logical flavor charges need not each have their own microscopic chemical potential: independent hopping amplitudes with a common density field already give the desired k_a.

[proved for the first-order model] The half-density marker tuning in MICROSCOPIC_HAMILTONIAN.md also checks directly. A payload dispersion 1-2t_a cos k, with t_a = 1/(2 cos k_a), has energy per cell e_a = (k_a-tan k_a)/pi. Set J_length = 1/2 and the remaining marker field to zeta = 1/pi - sum_a e_a. The first-order correction per physical site is then

\[
\frac{\sum_a e_a+\zeta-\sin(\pi\nu)/\pi}{b+\nu}
=\frac{1-\sin(\pi\nu)}{\pi(b+\nu)}.
\]

It is nonnegative and vanishes only at nu = 1/2 in [0,1]. The minimum is interior and has positive curvature. Thus the desired reference density can be selected by an available microscopic marker field; it is not an externally imposed thermodynamic sector. Continuity of this choice at finite s still uses the local effective-action input.

[proved, equalizer algebra] If the two endpoint helix states u_1,u_2 have overlap deficit chi = 1-|<u_1,u_2>|^2, the matrix

\[
E_L=\frac{r_1^{-1}(I-|u_2\rangle\langle u_2|)
+r_2^{-1}(I-|u_1\rangle\langle u_1|)}{\chi}
\]

has expectation 1/r_p in u_p. Multiplication by it cancels the endpoint occupation factor for both lengths. Both r_p are positive by positivity of the finite vacancy ground state. Its norm is finite at fixed b and must be included in B_p.

[proved, transport algebra] A cut shift can preserve k occupied core records in both cells only by transferring a vacancy. H_0 to H_2 is therefore the required conversion; the identity gives zero overlap. The unchanged ordered occupied sequence leaves the logical state untouched. For the two additional helix tracks the unchanged-site overlap is exp(i b phi) cos(phi/2)^(2b), canceled in phase by the transferred-site operator. This gives a strictly positive scalar transport coefficient for 0 < phi < pi. At phi = 1/b it remains nonzero and the compensator grows polynomially in b when r_p stay bounded away from zero. No global copying of the encoded state is used.

[proved from the C8 bounds] For a legal fixed marker pattern the cells are independent at s = 0. Selected cells have finite gaps; wrong lengths have energy excess at least gamma(1+|ell-b|). Illegal selector sectors cost at least U per violation. Long marker-free fragments have positive volume cost relative to the negative selected energy density. Thus the unperturbed result is stronger than a single global spectral gap. The displayed defect-energy bound in IR_AUDIT.md follows in each legal marker sector. This is useful input for a local elimination argument.

[numerical] Local copies of tagged_vacancy_checks.py, length_filter_checks.py, arrow_wrapper_checks.py, selector_wrapper_checks.py, helix_equalizer_checks.py and joint_fluid_checks.py all pass. The two-cell helix overlap test uses dimension 262,144. The independent physical marker/payload graph spectra in the last script have dimensions 162,340,1800. These are finite wrapper and first-order checks, not simulations of the autonomous high-channel core or of its finite-s infrared limit.

[proved conditional on the projected Hamiltonian] The joint ring quotient rotates both gap and payload coordinates. Its momentum constraint does not eliminate either family of particle-hole oscillators. C8's spectral reduction applies to any translation-invariant logical payload, including the present even fermionic Hamiltonian in its spin representation. The common Jordan-Wigner parity and cyclic cell permutation affect global twists, which must be retained. The one-XX-payload test does not independently determine all many-flavor zero-mode sectors.

[assumed, continuum interface] C20 treats product compactification as an assumption and explicitly leaves physical gluing unspecified. It cannot independently certify the physical mixed-resonance formula. C13 supplies the marker-comb and squeezed-coordinate dictionary used above. Finite-range, marker-conserving, parity-even payload operators give finite fermion polynomials in squeezed order; no local odd Jordan-Wigner twist is forced. Extending this dictionary to the integrated quasi-local action is still tied to the local elimination input. C20's full Gaussian metric, rather than a product formula for separate K_a, is the appropriate test once mixed quadratic terms appear.

[proved, scope of the graded-route-A checks] The equation P M_e[a,b] P = tau_a tau_b M_e[a,b] makes diagonal syndrome weights logically even, so the stated repair of C9's odd hopping prescription is necessary. Real off-diagonal even/odd syndrome blocks have the stated parity and entrywise contraction. The copied route-A resolvent and separated-port scripts pass. These results strengthen a separate second-order route; they neither prove nor obstruct the first-order sparse controller. Importing their fixed-alphabet claim would be an interface error, but the present candidate does not need to do so.

## 5. Infrared control and order of limits

### 5.1 Positive reference and full scaling metric

[assumed] Adopt the physical Gaussian fields and symplectic normalization of C13/C20. With f_a = k_a/pi, a_0 = b+nu and positive physical velocities v_i, put

\[
T=\begin{pmatrix}a_0&0\\-f&I\end{pmatrix},\qquad
A=T^T\operatorname{diag}(v_i)T,\qquad
B=T^{-1}\operatorname{diag}(v_i)T^{-T}.
\]

[proved] T is invertible, so both A and B are positive. AB is similar to diag(v_i^2). The reference Luttinger matrix K_ref = T^{-1}T^{-T} satisfies K_ref A K_ref = B. For N = 0 this gives K_marker = 1/a_0^2, recovering the length-fluid limit. These checks establish positivity and the geometric normalization of the reference, not of the unknown dressed microscopic coefficients.

[proved] With density-current coupling C, energetic stability requires H_G = [[A,C],[C^T,B]] > 0, equivalently A > 0 and B-C^T A^{-1} C > 0. The norm bound ||H_ref^{-1/2} delta H H_ref^{-1/2}|| < 1 is sufficient. It is a relative norm, not a comparison to the fastest velocity. In particular large b makes the raw physical matrices ill-conditioned.

[assumed, standard Gaussian formula] With J = [[0,I],[I,0]], the vertex scaling metric is

\[
\mathcal M=H_G^{-1/2}
\left|H_G^{1/2}JH_G^{1/2}\right|H_G^{-1/2},\qquad
\Delta(\lambda)=\tfrac14\lambda^T\mathcal M\lambda.
\]

[proved under this formula] At a fixed positive reference, the metric is continuous because the finite-dimensional chiral velocity matrix has no zero eigenvalue. There is consequently a smaller open neighborhood in which M >= kappa M_ref for any chosen 2/3 < kappa < 1. Positivity alone is a larger condition, as the explicit six-cycle countercheck demonstrates. At a stable point each oscillator family remains, irrespective of unequal velocities.

[numerical, illustrative Gaussian matrices] independent_checks.py verifies K_ref A K_ref = B for (N,b) = (1,4),(6,20),(21,100). The smallest raw reference eigenvalues are approximately 0.0492, 0.00238 and 0.0000990. Small random relative perturbations including C retain normalized scaling-metric lower bounds approximately 0.994,0.994,0.993. These are continuity checks on constructed Gaussian matrices, not calculations of the candidate's dressed stiffness. The copied C20 metric and vertex-shell counting checks also pass for up to eight payload flavors.

### 5.2 Fixed nonzero coupling is compatible with all finite margins

[proved conditional on local effective bounds] Fix p, then a finite exact core and occupied length, then the wrapper. All port time fractions, vacancy overlaps, compensator norms and finite gaps are fixed positive or finite numbers. Normalize V_p so that P V_p P gives the intended length/payload Hamiltonian in common units and let its local norm be bounded by B_p. Write H(s) = H_base + s V_p. At any fixed connected perturbative order, corrections begin at s^2, while all retained velocities are O(s).

[proved conditional on a uniform local bound] If the eliminated action obeys ||delta H_eff||_local <= C_p s^2 B_p^2/Delta_p, then its relative perturbation is at most C_p s B_p^2/Delta_p. Choose a positive dimensionless tolerance epsilon_p,b inside all normalized stiffness, scaling-dimension, momentum-shift and oscillatory-RG margins. A sufficient choice has the form

\[
0<s<\min\left\{
\frac{\Delta_p}{4C_pB_p},
\frac{\epsilon_{p,b}\Delta_p}{4C_pB_p^2}
\right\}.
\]

Changing harmless numerical constants does not affect the argument. For every fixed finite p,b this interval is nonempty, even if its upper endpoint is extremely small. No numerical value of s is established until the input bounds are supplied, but there is no algebraic requirement that s vanish with ring length.

[heuristic, oscillatory RG qualification] A nonzero mismatch does not stop a relevant perturbation from growing before its oscillation length. For a dimensionless coupling g and detuning Q in cutoff units, the usual requirement is |g| |Q|^{-(2-Delta)} much less than one; for a marginal channel one also controls |g| log(1/|Q|). At fixed p,b all relevant detunings have positive lower bounds. These conditions can therefore be included in epsilon_p,b. Quoting a mismatch alone without a coupling comparison would be insufficient. The notes' small-error requirement can be strengthened this way without changing the asymptotic alphabet count.

[proved, exact versus approximate core] The algebraic-rotation construction has no preparation error. In the optional approximate circuit, each actual gate remains unitary and the history parent still has exact zero-energy histories. Approximation changes port matrices, not automatically the code energies. For ||W'-W|| <= delta and ||F|| <= 1, the compression error is at most 2delta; the two-cell bond error is at most 4delta. Common history readout factors multiply signal and error and cancel in their ratio. Thus delta must be chosen small relative to eta times the desired first-order tolerance before selecting s. A genuine fixed core splitting would instead require splitting much less than s eta and would become worse relative to the target as s decreases. SCALE_HIERARCHY.md correctly distinguishes these cases.

[proved, order of limits] The viable order is fixed p, fixed finite circuit/core length and couplings, fixed wrapper, fixed small s > 0, then L to infinity. Finally one compares the separately defined models as p grows. The claim does not exchange the p and L limits or establish a single fixed-d family with unbounded c. Since N grows as p^3 while the conditional d bound grows as p^2 log p, the claimed c/d divergence follows once c = N+1 and the uniform compiler constant are established.

### 5.3 What the local elimination argument does and does not establish

[proved] The complete retained manifold has an extensive bandwidth after adding s P V_p P. It is not a globally isolated spectral band at fixed s and arbitrarily large L. A global Schrieffer-Wolff band argument or an O(L s^2) spectral error cannot control conformal gaps of order 1/L. The candidate correctly rejects that argument.

[proved] In a fixed legal marker sector, the static unperturbed bound has the form

\[
H_{\rm base}-eL\ge
\Delta_*\sum_{j\ {\rm selected}}Q_j
+\gamma\sum_{j\ {\rm wrong}}(1+|\ell_j-b|).
\]

Selected Q_j act on disjoint finite cells. Wrong-cell costs grow with their length. This controls unperturbed propagators with specified defects and excludes a zero-cost arbitrarily long cell. Exact species and marker symmetries are preserved by virtual corrections. These are meaningful, model-specific inputs for a physical integration argument.

[heuristic, serious verification gap] IR_AUDIT.md equations (9)-(10) claim geometric local bounds on connected excursions and exponential fast-time kernels. Those bounds do not follow from the displayed defect inequality by the short counting argument alone. One must define the connected reorganization when markers move, separate retained propagation without expanding its arbitrarily small denominators, and control insertions on regions visited earlier in an excursion. Such insertions need not touch a currently active defect and need not commute with earlier vertices. Saying that all distant retained insertions cancel does not by itself give the claimed order-independent constant C_p. An order-r coefficient being finite for every fixed r also does not establish one nonzero convergence interval for all r.

[heuristic] None of this identifies an instability. Integrating out finite-gap local machinery is a conventional and plausible physical procedure, and the defect-energy bound makes it more credible here than in an arbitrary extensively degenerate model. Nevertheless the needed quasi-local o(s) action for this moving-core system is the substantive microscopic-to-infrared input. Calling the proposed geometric estimate “standard” is not an independent verification of its applicability. The notes themselves acknowledge it is not proved.

[heuristic, sufficient repair at the physics standard] I would accept a controlled low-frequency derivation that constructs the retained/defect variables for this wrapper, states the local norm and cutoff, demonstrates why connected irreducible fast kernels have a positive dressed threshold and decaying spatial tails, and bounds their leading correction relative to s uniformly in L. A rigorous exact block diagonalization of the entire retained band is neither required nor appropriate. A conventional physical resummation with explicit assumptions and a checked small parameter would suffice. Two- and three-cell second-order calculations using actual compiled data would provide useful checks of the normalization and defect gaps, but cannot alone replace the uniform low-frequency argument.

[heuristic, conditional infrared conclusion] If that local-action input holds, the remaining RG step is conventional. There is no symmetry-allowed non-Gaussian relevant or marginal resonance inside the stated metric/momentum neighborhood. Currents shift equilibrium densities and currents, and quadratic current products change the Gaussian matrices. All other resonant vertices have a strict irrelevance margin. Hence a sufficiently small symmetry-preserving deformation is predicted to have N+1 compact-boson modes. No specific sign of an unknown quartic coefficient is needed for this implication.

[assumed, convention for c] The repository permits a collection of conformal sectors with different velocities and defines c as their sum. Under this convention each stable boson contributes one. A Casimir coefficient alone measures sum_i c_i v_i and is not an independent measurement of the unweighted c. Finite-size charge rounding and ring gluing can move zero-mode offsets with L. They do not remove oscillator families, but they also do not prove convergence to one fixed untwisted vacuum tower along every integer L. If that stronger convention were imposed, a separate construction would be needed. I do not impose it in this review.

## 6. Severity and acceptance conditions

[heuristic, classification] No fatal flaw was found in the reviewed construction. The following table separates missing verification from errors in statements or bookkeeping.

| Status | Severity | Finding | What would settle it |
| :-- | :-- | :-- | :-- |
| Assumed end-to-end finite implementation | Serious | The full autonomous sparse controller and its all-input, all-history-vertex contract have not been independently instantiated. The component lemmas pass. | A complete symbolic controller/compiler with invariant and resource proof, preferably executable component-to-core tests. A dense printed Hamiltonian is unnecessary. |
| Heuristic local elimination | Serious | The connected-excursion bounds and uniform o(s) low-frequency action are asserted, not established by the static defect count or finite-order checks. | A controlled physical derivation for the mobile projectors and actual microscopic scales. A new rigorous CFT theorem is unnecessary. |
| Proved unit conversion | Fixable | Momentum margins are sometimes described without specifying cell versus physical units. | State Q_phys = Q_cell/(b+nu) and the effective cutoff when imposing RG bounds. |
| Proved logical distinction | Fixable | Positive stiffness alone does not guarantee the six-field irrelevance bound. | Retain and verify the full relative scaling-metric condition. The explicit relevant six-cycle countercheck shows why. |
| Numerical scope | Fixable | Finite script passes can be read as stronger evidence for the complete core or thermodynamic phase than they supply. | Keep circuit, clock, initialized-core and finite-coupling checks distinct. None of the existing tests measures the candidate's c. |
| Assumed ring convention | Fixable if left implicit | C20 alone does not supply physical compactification; C13's dictionary and zero-mode conventions are needed. | State the physical gluing and unequal-velocity, moving-sector convention in the final claim. |
| Unexpanded numerical C | Cosmetic for existence | No optimized compiler constant, first violating prime or manageable matrix has been supplied. | A verified uniform finite controller makes finite C sufficient for the asymptotic argument. Numerical optimization is optional. |
| Not a defect | None | Six-field resonances exist, visibility/gaps can be tiny, and the graph does not give exact U(1)^N. | The full metric bound, fixed-p hierarchy and exact U(1)^q already address these issues. |

[heuristic, final recommendation] The appropriate current statement is a conditional physics-level counterexample construction with independently checked finite algebra and a viable protected infrared target. The review does not produce a competing instability or disprove the family. It also does not certify a complete microscopic refutation. To make the latter acceptable, close the two serious verification gaps above and present the combined Hamiltonian and inference in one consistent set of conventions. An optimized dimension, a simulation at a violating prime, and a mathematical proof of the interacting CFT limit are not prerequisites at the stated physics standard.

## 7. Reproducibility and limits

[numerical] All ten requested contribution checks pass from local copies: exact_route/verify.py, torus_bottleneck/checks.py, graded_route_a/checks.py, separated_port_cells/checks.py, graded_route_a/sweeping_interface_checks.py, exact_route/sweep_checks.py, exact_route/amplitude_checks.py, exact_route/helper_checks.py, graded_route_a/clock_schedule_checks.py, and research_continuation/ir_audit_checks.py. The six C8 checks named above also pass, as does C20/check.py. Each run took seconds on this machine. Captured outputs for the contribution and C8 scripts are beside their local copies under reproduction/.

[numerical] New independent checks are in independent_checks.py, with output in independent_results.json and independent_checks.stdout.txt. They cover sparse two-port matrices and charges, projective incidence and all two-edge charge multisets, actual resonant six-cycles, the acceptance-gap lemma and terminal-pin failure, the physical Gaussian reference, and an explicit positive-stiffness relevant-cycle countercheck. Run `python3 independent_checks.py` from this folder. Only NumPy, SciPy and the Python standard library are used.

[proved, limits of the review] I did not execute a complete autonomous projective-plane core, diagonalize its full initialized Hamiltonian, compute its dressed finite-s stiffness, or measure its central charge. I did not infer any of those results from the passing small-matrix tests. Those limits are part of the verdict.
