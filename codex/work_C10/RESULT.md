# C10. An explicit finite-program track, and the obstruction to the stronger claim

This is a partial construction with a concrete obstruction, not a fixed-dimension counterexample to the c–d conjecture.

I constructed a translation-invariant nearest-neighbor track for any hard-coded finite nearest-neighbor gate program. I prove its path-clock spectrum, initialization bound, illegal-sector bound, and preservation of the supplied selector ordering. A three-time program implements the left Route B Hamiltonian under history compression using only a workspace SWAP gate, with Bell preparation appearing as an energy term. The compact version has k_c=12 and d0=96. Section 11 gives a stronger version with k_c=24 and d0=192 whose propagation leaves every physical port register idle even in illegal sectors. Both actual two-site bonds are implemented, in `build_bond_c10.py` and `build_strict_bond_c10.py`, respectively.

Two requested conclusions are not established. First, the hard-coded alphabet grows with program length and workspace size. I have not supplied a fixed-alphabet decoder of N from b or a coupling phase. Second, selected-time energies generally do not have uniform propagation histories as their exact ground states at finite clock coupling. This is an algebraic obstruction, demonstrated below in the explicit Route B program. Exact port compression holds on the ideal history code, but not on that program's actual finite-J ground space.

Everything below was derived from the supplied context. No existing files, outside references, or external data were read. Only files created for this task were accessed. The supplied selector and static Route A/B conclusions are inputs, not conclusions re-derived here.

## 1. Conventions and scope

The distance between consecutive selector markers is b. There are m=b-1 interior sites. The chain length is denoted by L_chain, and the number of clock times by n=T+1. This avoids confusing the physical length with the clock length in the tower formula.

The local space is

    C^4_selector tensor C^(k_c)_computation tensor C^q_data.

The selector alphabet remains (#,L,X,R). The computation alphabet uses A and R_j as new symbols, unrelated to the selector letters. The data register at an actual marker is pinned to a fixed state |0> by a positive onsite term. This removes an otherwise free q-fold degree of freedom at every marker.

The prescribed circuit is first written as an adjacent head walk

    j_0,j_1,...,j_T in {1,...,m},     |j_(t+1)-j_t|=1.

The gate W_t acts on the bond traversed in step t. Its tensor order is always physical left-to-right order, including on a leftward step. A desired nearest-neighbor gate list can be put in this form by moving the head to each gate bond with identity gates, traversing that bond with the desired gate, and continuing to the next bond. Energy insertions can likewise be scheduled when the head is on an endpoint of their support bond. Additional identity traversals supply additional clock times if needed.

In the valid cell, sites 1 and m are the ports. Every W_t on a bond incident to either of these sites is identity on its port factor. In particular one can use identity for every traversal of a port bond. All nontrivial circuit gates act on sites 2 through m-1.

This is a compilation of prescribed nearest-neighbor circuits. Arbitrary real unitaries are not asserted to have exact finite decompositions into an arbitrarily preselected finite universal gate set. The Route B Pauli implementation in Section 7 does have an exact finite Clifford gate set.

## 2. Alphabet and allowed pairs

For the hard-coded program, use

    C = {M} union {A_j,R_j : 1 <= j <= m} union {H_t : 0 <= t <= T}.

Thus

    k_c = 1 + 2m + (T+1) = 2m+T+2.

Associate position j to A_j and R_j, and position j_t to H_t. The allowed ordered computation pairs are exactly

    A_j A_(j+1), R_j R_(j+1)                  1 <= j < m,
    M A_1, R_m M,
    A_(j_t-1) H_t                            j_t > 1,
    M H_t                                    j_t = 1,
    H_t R_(j_t+1)                            j_t < m,
    H_t M                                    j_t = m.

Every unlisted ordered pair is forbidden, including MM, A_j R_k, R_j A_k, and H_t H_u. Unused letters are retained so the size formula is literal.

Add an onsite synchronization penalty for

    selector # with computation != M,
    selector != # with computation M.

No condition distinguishes selector L, X, and R.

### Legal-string lemma

Fix two correctly synchronized markers with ell interior sites, and suppose all computation pairs between them are allowed.

Every allowed pair not involving M increases the position index by one. The left boundary forces the first position index to be 1, and the right boundary forces the last to be m. Therefore ell=m.

The left boundary permits an A or a head, never an R. After an A, only an A or a head is allowed. After a head, only an R or the right marker is allowed. After an R, only an R or the right marker is allowed. A string with only A's cannot terminate at the right marker. A string with only R's cannot start at the left marker. It follows that there is exactly one head.

Consequently the complete list of legal cell strings is

    C_t = A_1 ... A_(j_t-1) H_t R_(j_t+1) ... R_m,    0 <= t <= T.

In particular, the argument excludes no-head strings, multiple-head strings, and all incorrect cell lengths. It is not a counting assumption about the low-energy sector.

## 3. Transition rules and the full bond

If the head moves right, use

    H_t R_(j_t+1)  <->  A_(j_t) H_(t+1),

with W_t on the forward transition. If it moves left, use

    A_(j_t-1) H_t  <->  H_(t+1) R_(j_t).

Again W_t acts on the forward transition. Write a_t and b_t for the left and right two-letter patterns in these rules, respectively. The propagation bond contribution is

    P_t = (J/2) [ |a_t><a_t| tensor I + |b_t><b_t| tensor I
                 - |b_t><a_t| tensor W_t
                 - |a_t><b_t| tensor W_t^dagger ].                 (3.1)

For data vectors x,y, its quadratic form is (J/2)||y-W_t x||^2. Thus P_t is positive semidefinite for either direction of head motion.

Let p=|#><#| on the selector, chi=1-p, and M=|M><M| on the computation track. Define the following onsite operators, with omitted identities understood:

    S_sync = p(1-M) + chi M,
    S_marker = p tensor I_comp tensor (I_data-|0><0|),
    S_head = chi tensor sum_t |H_t><H_t| tensor I_data,
    A_site = kappa S_sync + kappa_M S_marker + C S_head - rho p.    (3.2)

Let Q_comp be the projector onto the forbidden computation pairs in Section 2. It acts independently of the selector and data registers.

Each scheduled local energy e_(t,a) is attached to the legal adjacent computation pattern u_(t,a) containing H_t on its support bond. For example, an energy on the head's left bond has pattern A_(j_t-1) H_t; an energy on its right bond has pattern H_t R_(j_t+1). Several terms may be attached to the same clock time. Define

    E_bond = sum_(t,a) |u_(t,a)><u_(t,a)| tensor e_(t,a).           (3.3)

Initialization checks are attached by the same prescription, with positive workspace projectors and their own positive coefficients. Checks on the head's own data register can instead be onsite and distributed equally over its two bonds. No Bell-state condition is included in Q_comp or S_sync.

The complete translation-invariant bond is

    h = h_selector
        + kappa Q_comp
        + chi_i chi_(i+1) [ sum_t P_t + E_bond + I_bond ]
        + (A_site tensor I + I tensor A_site)/2.                  (3.4)

Here I_bond denotes the optional positive initialization terms. The selector bond is the supplied one, written with its original chemical potential mu:

    h_selector = K Q_selector
                 + (mu/2)(p_i+p_(i+1))
                 - t_sel (|LX><XR|+|XR><LX|),                    (3.5)

tensored with identity on the other tracks. The coefficient J is the clock coupling and is unrelated to t_sel. Throughout, J,kappa,kappa_M,nu and Delta are positive, t_sel>0, K>4 t_sel, and g and mu are real. The symbol nu is used for the extra pins in Section 11.

The factors chi_i chi_(i+1) prohibit every propagation and time-energy term across an actual selector marker. Marker-data pinning also never involves a neighboring register. Hence the computation decomposes into finite segments for each fixed marker pattern, even in the presence of malformed computation strings.

The values of C and rho are specified in Section 8. All coefficients, patterns, and matrices in (3.4) are position independent. On a ring, the half-onsite convention gives each onsite term exactly once. For an open chain one adds the missing endpoint halves, or fixes exterior markers.

### Commutation with the selector

Every computation term is either independent of the selector or depends only on p and chi. The selector hopping changes LX into XR and its reverse, so it never changes whether either site is a marker. Both p and chi commute with each one-site matrix factor of that hopping. The selector penalties and chemical potential are diagonal and also commute with the computation terms.

This proves [H_selector,H_comp]=0, including terms on overlapping bonds. Same-bond commutation alone would not have sufficed.

The circuit on the legal cell is port-idle. There is a distinction in illegal configurations: a workspace position label placed next to a real marker can make a propagation rule act on that site's data before its penalty is paid. The construction does not assert that every propagation term commutes with every marker-defined port operator on the entire illegal Hilbert space. The exact circuit and port statements for this compact version apply to the invariant legal subspace. Section 11 repairs this issue and establishes port-idle propagation on the entire Hilbert space.

## 4. Exact propagation histories and initialization

Set

    U_0=I,       U_t=W_(t-1)...W_0,
    V psi = n^(-1/2) sum_(t=0)^T |C_t> tensor U_t psi.             (4.1)

Each legal configuration has only the propagation edges to its predecessor and successor. The head labels H_t are distinct even when the head revisits the same physical position. Therefore the legal graph is a path with n vertices, not a cycle.

Conjugating by sum_t |C_t><C_t| tensor U_t gives

    H_prop = J L_n tensor I,

    L_n = (1/2) sum_(t=0)^(n-2)
                   (|t>-|t+1>)(<t|-<t+1|).                     (4.2)

Its quadratic form vanishes exactly on constant vectors, proving ker(H_prop)=ran(V).

For completeness, the eigenvectors of L_n can be chosen proportional to

    f_t^(k) = cos[(k pi/n)(t+1/2)],    0 <= k <= n-1.

The interior second difference gives (1-cos(k pi/n))f_t^(k). At the left endpoint, (f_0-f_1)/2 gives the same eigenvalue by the cosine subtraction identity. The right endpoint follows by reflection. These n independent vectors exhaust the spectrum. Thus the clock tower and its gap are

    E_k = J[1-cos(k pi/n)],
    Gamma_clock = J[1-cos(pi/n)].                                (4.3)

### Correct initialization

An initial identity-gate sweep can check each intended ancilla when the head reaches it. Let Pi_a be the initial workspace projectors being penalized, let the checks occur at times tau_a, and assume the Pi_a commute. Product-state ancilla checks satisfy this assumption. The ports and Bell partner are not included in these checks.

Because the sweep precedes all nontrivial gates, the rotated initialization Hamiltonian is

    H_in = sum_a kappa_a |tau_a><tau_a| tensor Pi_a.               (4.4)

Its zero space within the propagation kernel is exactly V applied to the correctly initialized input space. This proves the requested exact history-ground statement for propagation plus these initialization terms.

Here is an explicit gap bound. Set kappa_0=min_a kappa_a. Decompose the input space into simultaneous eigenspaces of the Pi_a. A wrong-input eigenspace has at least one check equal to 1, at some anchor time a. For any normalized clock vector f,

    f_t = f_a + sum_(edges between a and t) signed differences,

so Cauchy-Schwarz gives

    |f_t|^2 <= 2|f_a|^2
               + 2|t-a| sum_(edges between a and t) |f_(s+1)-f_s|^2.

Summing and using sum_t |t-a| <= n(n-1)/2 yields

    1 <= 2n |f_a|^2 + n(n-1) sum_s |f_(s+1)-f_s|^2.

Its energy is at least

    kappa_0 |f_a|^2 + (J/2) sum_s |f_(s+1)-f_s|^2.

Consequently, for n>=2,

    Gamma_input >= kappa_0 J / [2n(J+(n-1)kappa_0)].              (4.5)

The full gap of H_prop+H_in is at least min(Gamma_clock,Gamma_input). For n=1 the wrong-input gap is kappa_0.

The commuting-check assumption is used in this proof and is not dropped. General additional time energies can mix initialized and uninitialized input spaces. Calling them separate sectors then requires a further commutation proof. The short example below has no initialization ancillas, so this issue is absent there.

## 5. What arbitrary time energies actually do

Within the legal sector, write

    D = sum_t |C_t><C_t| tensor h_t,
    A_t = U_t^dagger h_t U_t,
    H_bar = (1/n) sum_t A_t.                                    (5.1)

Then V^dagger D V=H_bar. This equality is exact. It does not imply that ran(V) is invariant under D.

For an individual history V psi to be an eigenvector of H_prop+D with eigenvalue E, the independent clock components require

    A_t psi = E psi                  for every t.                (5.2)

Conversely, (5.2) is sufficient. Thus it is a necessary and sufficient condition, not a perturbative estimate.

For the entire history space to be invariant, all A_t must be one and the same operator. To see this, invariance implies D V=V A for some A. Comparing each clock component gives A_t=A. The converse is immediate.

A useful measure of the failure is

    ||(1-VV^dagger) D V psi||^2
        = (1/n) sum_t ||(A_t-H_bar)psi||^2.                      (5.3)

To derive it, D V psi has clock components U_t A_t psi/sqrt(n), while V H_bar psi has components U_t H_bar psi/sqrt(n). Their difference is the orthogonal complement of the history projection.

There are two sufficient special cases for exactness.

1. Each h_t is positive and annihilates the desired propagated code U_t G. Then the sum is frustration-free, and its ground space is the intersection of the propagation kernel with those energy kernels.
2. Every A_t equals the same target Hamiltonian. Then the legal Hamiltonian is unitarily equivalent to J L_n tensor I plus I_clock tensor H_target.

Neither condition is automatic for selected-time energies. In the Route B example, one time has h_t=0 and the target ground energy is negative. Equation (5.2) would force E=0 at the idle time, whereas the history variational state has negative energy. Therefore its uniform target histories cannot be exact ground eigenvectors for any finite J. Adding C at every clock time changes both energies by C and does not repair the contradiction. This is compatible with the supplied 1/J low-energy convergence; it is the stronger exact finite-J ground-history identification that fails for this program.

### Port operators

For a port operator F, port-idle gates give [F,U_t]=0 and hence

    (I_clock tensor F)V=VF,
    V^dagger(I_clock tensor F)V=F.                               (5.4)

If G is the ideal static code, then

    (VG)^dagger(I_clock tensor F)(VG)=G^dagger F G.                (5.5)

Both statements are exact. They concern the propagation histories and the ideal history code. They do not identify VG with the actual finite-J ground space in the presence of arbitrary D. Port operators also need not preserve the smaller energetic code G itself.

### A controlled finite-J replacement

For the uninitialized full history space, let P=VV^dagger, Q=1-P, v=||D||=max_t||h_t||, and Gamma=Gamma_clock. In an eigenvalue window |E|<=v, if Gamma>2v, then

    Q(H_prop+D-E)Q >= Gamma-2v.

Eliminating the Q component gives the exact Schur complement

    [PDP - PDQ (QHQ-E)^(-1) QDP] x = E x,

with

    ||PDQ (QHQ-E)^(-1) QDP|| <= v^2/(Gamma-2v),
    ||Q psi|| <= [v/(Gamma-2v)] ||P psi||.                       (5.6)

These are bounds on the correction and history leakage, not an exact-history theorem. With initialization one must use the appropriate separated space and its proved gap, and control whether D mixes the input sectors. No thermodynamic critical-stability theorem follows merely from (5.6).

## 6. A literal short Route B program

Take q=2, m=4, b=5. The data registers are

    site 1 = left port p,
    site 2 = partner s,
    site 3 = one logical qubit l,
    site 4 = idle right port r.

Only the left Route B gadget is implemented here. The right port is a free spectator. The two-ended mirror construction is not tested by this example.

Use F_1=X, F_2=Y and gamma_1=X_l, gamma_2=Y_l. With

    P_Omega = (I + X_p X_s - Y_p Y_s + Z_p Z_s)/4,

Pauli anticommutation gives

    {X_p,P_Omega} = (X_p+X_s)/2,
    {Y_p,P_Omega} = (Y_p-Y_s)/2.

For example, in the first identity the I and X_p X_s terms each contribute twice their product with X_p, while the Y_p Y_s and Z_p Z_s anticommutators vanish. The second identity follows in the same way, with the minus sign of the YY Bell term. Thus

    B+B^dagger = [X_p X_l + Y_p Y_l + X_s X_l - Y_s Y_l]/2,

and the target is

    H_L = Delta(1-P_Omega)
          + (g/2)(X_p X_l+Y_p Y_l+X_s X_l-Y_s Y_l).              (6.1)

Use three clock times at positions (2,3,2). The alphabet is

    M, A1,A2,A3,A4, R1,R2,R3,R4, H0,H1,H2.

The allowed pairs, explicitly, are

    M A1, R4 M,
    A1 A2, A2 A3, A3 A4,
    R1 R2, R2 R3, R3 R4,
    A1 H0, H0 R3,
    A2 H1, H1 R4,
    A1 H2, H2 R3.

The propagation rules are

    H0 R3 <-> A2 H1        with SWAP on the two data registers,
    A2 H1 <-> H2 R3        with identity on the data registers.

In the legal cell both rules act on bond (2,3). Thus U_0=I and U_1=U_2=SWAP_(s,l).

There are three energy-pattern terms in the uniform bond:

    |A1 H0><A1 H0| tensor 3 Delta(1-P_Omega),
    |H0 R3><H0 R3| tensor (3g/2)(XX-YY),
    |A1 H2><A1 H2| tensor (3g/2)(XX+YY).                         (6.2)

As in (3.4), each is multiplied by chi_i chi_(i+1). Consequently the legal-sector h_t are

    h_0 = 3 Delta(1-P_Omega,ps) + (3g/2)(X_s X_l-Y_s Y_l),
    h_1 = 0,
    h_2 = (3g/2)(X_p X_s+Y_p Y_s).

Their compression is

    (h_0 + SWAP_(s,l) h_2 SWAP_(s,l))/3 = H_L.                   (6.3)

This proves the desired algebraic energy synthesis using the finite gate set {I,SWAP}. The Bell term is an energy of strength 3 Delta, with no Bell restriction in the allowed pairs, input checks, or synchronization constraints. No circuit gate acts on either port.

For the supplied Route B isometry G, choose the sign convention

    theta = (1/2) arctan(2g sqrt(2)/Delta),
    G = cos(theta) A - [sin(theta)/sqrt(2)] B A,

where A maps a logical state to |Omega> tensor that state, leaving the right spectator unchanged. The supplied target compression is

    G^dagger X_p G = alpha X_l,
    G^dagger Y_p G = alpha Y_l,
    alpha = -sin(2 theta)/sqrt(2).

Equation (5.5) gives these same identities on VG without approximation. Section 10 shows that the actual finite-J ground code has a different compression.

## 7. Exact finite gate-set energy synthesis beyond the small example

There is a useful finite-program extension, still with a growing computation alphabet.

Restrict this extension to q=2^u and choose the odd F_a from an orthonormal Hermitian Pauli basis compatible with the chosen grading. Take the logical Clifford generators in their usual Pauli-string representation. These are explicit choices for this extension, not assumptions about every possible F_a in the context.

The Bell projector has the expansion

    P_Omega = q^(-2) sum_(A in Pauli_u) A_p tensor A_s^T.

Substitution into Delta(1-P_Omega)+g sum_a {F_a,P_Omega} tensor gamma_a gives a finite real linear combination of Hermitian Pauli strings. Anticommuting F_a,A pairs contribute zero to the anticommutator. Commuting pairs produce Hermitian Pauli products with their computable signs. Each surviving term has the form

    c P_port tensor Q_workspace.

No term requires a gate on the port. If Q_workspace is nonidentity, a workspace Clifford circuit C_Q can be chosen so that

    C_Q^dagger Z_partner C_Q = Q_workspace.                     (7.1)

Here is a constructive reason. Single-qubit Clifford gates map every nonidentity factor of Q_workspace to Z. A sequence of CNOT gates combines the resulting product of Z's into a single Z, using the conjugation rule Z_target -> Z_control Z_target. Workspace SWAP gates bring its location to a selected partner qubit. Reversing that sequence gives the required direction in (7.1). Gates on nonadjacent workspace qubits are routed by nearest-neighbor SWAPs, without using either port.

Run C_Q, return the head to the partner, and insert the local two-site energy

    n c P_port tensor Z_partner.

Run C_Q^dagger before synthesizing the next monomial. Its history compression is c P_port tensor Q_workspace. If Q_workspace=I, use the port onsite energy P_port instead. Identity port factors are allowed. A direct Bell energy can be scheduled separately at a time with U_t=I. All intermediate times without an insertion have zero energy.

For fixed q, this uses a finite collection of port-Pauli energy primitives and the finite gate set H,S,S^dagger,CNOT,SWAP on workspace qubits. Gates within a q-dimensional site can be carried by an adjacent head transition with identity on the other data register. This proves exact synthesis of the compressed finite Route B Hamiltonian for these Pauli choices. It does not prove that the synthesized energy is frustration-free, or that its finite-J ground space equals the uniform history code.

The factor n in each energy is chosen after the finite schedule is fixed. It compensates the uniform history weight 1/n. The resulting couplings can be large and depend on the finite program, but not on L_chain.

## 8. Illegal sectors and preservation of the selector margin

The bounds in this section apply to the explicit hard-coded alphabet. They do not supply the missing bounds for a universal fixed-alphabet machine.

### Making the time-energy terms nonnegative in every malformed string

For each scheduled local term choose c_(t,a)>=0 such that

    e_(t,a) >= -c_(t,a) I.

Choose

    C >= max_t sum_a c_(t,a).                                   (8.1)

In an arbitrary computation basis string, each occurrence of H_t activates at most the local energy terms assigned to that head, on its two adjacent bonds. Some patterns may be absent, but their absence cannot increase the total negative lower bound. Therefore all activated energies associated with that head, plus C times that head's onsite projector, sum to a positive operator. This remains true when different heads' data supports overlap: operator lower bounds add without requiring commutativity.

Thus the entire time-energy Hamiltonian plus C S_head is positive on all strings, not just legal strings. Propagation, input checks, synchronization penalties, and marker-data pinning are also positive before the marker counterterm.

For the short program, the Bell term is positive and

    ||(3g/2)(XX-YY)|| = ||(3g/2)(XX+YY)|| = 3|g|.

Each head activates at most one of these potentially negative terms. Hence C=3|g| suffices.

### Illegal-cell lower bound

For a fixed actual selector interval, assign all its internal and endpoint-pair penalties to that interval and half of each endpoint's onsite synchronization penalty to it. Interior synchronization penalties are assigned wholly to that interval. Across the whole ring this counts every term once.

A cell is good when its endpoint marker labels, all its interior synchronization conditions, and every pair condition are satisfied. By Section 2 it then has length m and one of the strings C_t.

If a cell is bad, at least one diagonal condition is violated. A pair or interior synchronization violation costs kappa; an endpoint synchronization violation costs kappa/2 under this accounting. The diagonal penalty is therefore at least kappa/2 on the span of bad basis strings.

The good span is invariant. A propagation rule in a good cell produces the adjacent good history string and never crosses an actual marker. The energies and input checks do not change computation labels. By Hermiticity the bad complement is invariant as well. Therefore the penalty bound is a spectral bound, not just a classical expectation estimate.

Since all the remaining terms before the counterterm are positive,

    H_cell,bad >= kappa/2.                                      (8.2)

With correct endpoint marker labels fixed, the stronger bound is kappa. In particular (8.2) covers no head, several heads, wrong lengths, and a propagation history cut by an actual marker. Cycles in an illegal transition component do not evade it, because their possibly zero propagation energy cannot cancel a positive legality penalty.

A ring with no actual markers cannot have a zero-penalty computation string. Either a computation M creates a synchronization violation, or every pair would have to increase a position index around a closed ring, which is impossible. Such a ring has at least kappa of penalty before adding the other positive computation terms.

Marker data are never touched by propagation or time energies. Exciting them costs kappa_M independently. This is a separate positive gap, not a Bell or workspace legality condition.

### Counterterm that preserves the selector ordering

Let e0 be the ground energy of a good cell before the marker counterterm, including the C shift and any initialization penalties. It is nonnegative. Choose an explicit number rho with

    rho >= e0,             kappa/2 > rho.                       (8.3)

For a general finite program, a valid-input uniform history trial state proves

    e0 <= C+v,       v=max_t ||h_t||,

so rho=C+v is sufficient. For the short Route B program its ideal ground history gives

    e0 <= C+E_minus <= C,
    E_minus = [Delta-sqrt(Delta^2+8g^2)]/2 <= 0,

and rho=C suffices. No numerically determined eigenvalue is needed as a coupling.

The term -rho p subtracts rho per actual marker, hence once per cell. Let

    a = (e0-rho)/b <= 0.

A good cell has marker spacing b and ground energy a b after this subtraction. A bad cell of any marker spacing B has energy at least

    kappa/2-rho >= 0 >= a B.

Summing over all actual cells gives

    H_comp >= a L_chain,                                       (8.4)

and an all-good b-periodic arrangement attains the lower bound. More quantitatively, a bad cell contributes at least

    kappa/2-rho-a B >= kappa/2-rho

above the comparison density a. The marker-free sector also obeys (8.4), since its computation energy is nonnegative and a<=0.

For chain lengths divisible by b, the supplied selector ground arrangement and the computation ground arrangement can therefore be attained simultaneously. Relative to their common ground energy, computation adds a nonnegative amount to every selector sector. Its addition cannot reduce the supplied selector selection margin, including the order-t_sel/b^4 wrong-length margin. It can increase that margin. The supplied walker excitation energies are unchanged because the tracks commute and the computation is blind to L,X,R.

This is a comparison with the supplied selector Hamiltonian at its original mu. The coefficient of the bare marker term in the total physical bond is mu-rho. Merely keeping mu unchanged while adding an arbitrary cell ground energy would not justify preservation of the selector margin.

If desired, choosing rho=e0 makes every good cell have zero energy. That option is implemented, but it is unnecessary for the final construction and requires specifying e0. For incommensurate finite rings an all-good tiling is absent; (8.4) remains a bound, but I do not claim an all-good ground state there.

### Gap summary for this track

For propagation and commuting initialization alone, the gap above the desired initialized histories is bounded below by

    min{Gamma_clock, Gamma_input, kappa/2, kappa_M},

in addition to the supplied selector gaps, with absent initialization checks omitted. Wrong marker patterns retain at least their supplied selector margin.

For the program with time energies and counterterm, every bad computation cell remains separated from the all-good ground energy by at least

    kappa/2-rho > 0,                                            (8.5)

and marker-data excitations still cost kappa_M. A good cell's actual internal gap must be proved or calculated for its actual Hamiltonian. The bare tower (4.3) is not the full spectrum after arbitrary time energies are added. Claiming otherwise would repeat the exact-history error in Section 5.

## 9. What is missing for a fixed alphabet

The present allowed pairs hard-code the position labels j and the time labels t. Changing b or the program changes the alphabet and the bond table. In particular,

    d0 = 4 q (2m+T+2),                                         (9.1)

and this is not independent of a growing target N.

The finite gate set in Section 7 does not solve this problem. A fixed quantum gate set and a fixed computation alphabet are different requirements.

A potential replacement would use a fixed head alphabet and a counter/program tape inside the cell, decode N from b, and simulate the chosen program. That replacement has not been constructed here. Its missing proof can be located fairly sharply.

Local domain-wall rules can enforce one head without storing t in its label. They do not, by themselves, prove that every locally allowed tape belongs to the intended computation orbit. An unintended component may never visit an initialization-check configuration. If it is a path or cycle of identity-gate transitions with no violated check, its uniform history has zero energy. For nonidentity gates a cycle's holonomy may alter this conclusion, but there is no general positive bound merely from the presence of a cycle.

A completed construction therefore needs an explicit reversible decoder and transition table, plus a proof for every locally admissible malformed tape that its entire orbit either meets a penalty or has a quantitatively positive propagation/check energy. The proof must include reset behavior, counter overflow, insufficient tape length, marker truncation, and wrong lengths. It must provide an input-to-history gap before invoking a large J. I have not supplied these pieces, nor a coupling-phase decoder.

There is also a separate choice about requirement (b). If exact uniform ground histories are retained, the energy program needs an additional exact mechanism such as a frustration-free implementation annihilating the desired code at every occupied time. The Pauli energy compiler above does not provide that mechanism. If the intended requirement instead permits an approximate low-energy history band, estimates such as (5.6) apply, and finite-J changes of the port operators must be carried through the rest of the construction.

## 10. Numerical verification

All computations used numpy/scipy locally and completed in seconds. No large computation was attempted. Reproduce them with

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 verify_c10.py
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 build_bond_c10.py
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 verify_full_cell_c10.py
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 build_strict_bond_c10.py

The executable named `python` was unavailable; the successful runs used `python3`.

### Legality and the actual bond

All 12^4=20,736 four-site computation strings were checked, including interior M labels. Exactly three have zero penalty:

    A1 H0 R3 R4,
    A1 A2 H1 R4,
    A1 H2 R3 R4.

Their edges are exactly 0--1 and 1--2. A transfer-matrix count for lengths 1 through 8 finds no legal length other than 4, consistent with the all-length proof.

The implemented local dimension is 96, so the two-site matrix has shape 9216 by 9216. It has 10,470 nonzero entries. The computed Hermiticity error and same-bond selector commutator are both zero. The overlap-bond commutation is supplied by the proof in Section 3.

The final default parameters are

    J=100, kappa=20, kappa_M=20,
    K=5, t_sel=1, mu=0.9,
    Delta=1, g=0.2, C=rho=0.6.

Thus K>4 t_sel, and the bare marker coefficient in the total bond is 0.3. The numerical value mu=0.9 was used to test the bond, not to perform a new selector optimization. The period-selection conclusion is conditional on the supplied choice of selector parameters selecting b=5. The preservation proof applies to any such choice, and no selector optimization was rerun.

An independent contraction of the physical bond over the two boundary markers and the legal clock basis agrees with the 48-dimensional legal-sector matrix. With the final explicit counterterm the measured difference is zero at the tested floating-point precision. The earlier exact-counterterm version differed by 1.42e-14.

### Propagation and exact energetic Bell preparation

For J=7 and n=3 the bare clock eigenvalues are

    0, 3.5, 10.5,

matching J[1-cos(k pi/3)].

With only the positive energy 3 Delta(1-P_Omega) at time zero, the ground space is the uniform propagation history of a Bell pair and arbitrary logical and right-port qubits. Its dimension is four. The measured eigenvector residual is 4.4e-16, and the ground-projector error is 2.1e-14. Its gap is 0.635209481, above the bound 0.269230769 from (4.5), with kappa_0=3 Delta.

This is an exact history-ground example with Bell preparation by energy. It is also the g=0 limit of the Route B program. It does not include g(B+B^dagger).

### Route B compression and its finite-J failure

For Delta=1 and g=0.2,

    E_minus = -0.0744562646538,
    alpha   = -0.348155311911.

The Pauli expansion (6.1) has operator-norm error 7.9e-17. The history-compressed energy equals H_L to 3.8e-16, and both ideal port compressions agree to 5.6e-17. Equation (5.4) is satisfied with zero measured error.

However, the proposed uniform ground histories have eigenvector residual

    ||(H_prop+D)VG-E_minus VG|| = 0.328264163522.

The residual lies outside the propagation histories and does not decrease with J, because H_prop annihilates VG. The distance of the true ground space from VG does decrease with J, as perturbation theory predicts.

| J | Actual ground energy | Leakage out of full history space | Port-compression error |
|---:|---:|---:|---:|
| 10 | -0.09669107634 | 6.82933e-2 | 5.65939e-2 |
| 100 | -0.07658244945 | 6.51381e-3 | 6.96736e-3 |
| 1000 | -0.07466816055 | 6.49206e-4 | 7.08510e-4 |
| 10000 | -0.07447744731 | 6.49002e-5 | 7.09677e-5 |

Here the port error is the operator norm after aligning the actual four-dimensional ground basis with the ideal basis by polar decomposition. It is not just a basis artifact. At J=100 the actual compressed X_p and Y_p eigenvalues have magnitudes 0.355122674, whereas the ideal magnitude is 0.348155312. The invariant square test gives

    ||F_eff^2-alpha^2 I|| = 0.00489999224.

The four low eigenvalues are equal within numerical precision. Two of those degrees of freedom come from the idle right spectator. An analytic degeneracy proof for the finite-J synthesized Hamiltonian is not claimed. The gap above these four states at J=100 is 1.038144972.

The second-order correction computed directly with the clock pseudoinverse is

    -G_hist^dagger D U_diag (L_3^+ tensor I) U_diag^dagger D G_hist
       = -0.211818898541 I

to numerical precision. Thus its ground-energy prediction is E_minus-0.211818898541/J+O(J^(-2)), consistent with the scan. This is a numerical consistency check, not an all-orders equality.

### Every quantum sector of the four-site cell

`verify_full_cell_c10.py` decomposes all 20,736 computation strings into invariant transition components, includes all 16 data states for each, and diagonalizes every resulting block. This covers a quantum space of dimension 331,776 with correct fixed endpoint marker labels and marker data.

There are 19,449 singleton computation components, 426 components with three strings, and one with nine strings. The largest quantum matrix diagonalized has dimension 144. The nine-string component includes cycles arising from two independent heads, so the illegal-cycle bound is tested rather than bypassed by inspecting only single-head sectors.

The minimum illegal energy is 20 before the counterterm. The all-A and all-R no-head strings attain it. After the final rho=0.6 counterterm the minimum is 19.4, equal to the stronger fixed-endpoint bound kappa-rho. The general shared-boundary bound is kappa/2-rho=9.4.

The good-cell ground energy before the counterterm is 0.523417550545. After it, the ground energy is -0.076582449455. The illegal-to-good separation in this fixed cell is therefore about 19.47658245. These distinct energy conventions should not be conflated.

Endpoint-marker mismatches, arbitrary other lengths, and the full selector Hilbert space are not included in this finite ED. They are covered by the analytic synchronization, length, and commuting-track arguments. No two-cell intercell interaction or many-cell critical spectrum was diagonalized.

## 11. A version with port-idle propagation in every sector

The compact version economizes on data registers by using a site's q-dimensional register either as a port or as workspace, according to its legal position. There is a direct way to remove the illegal-sector qualification without introducing nonlocal terms.

Give every site a q-dimensional workspace register z inside its computation alphabet and retain its q-dimensional register p as the physical port register. The local factors are now

    selector s, clock letter a, workspace z, physical register p.

The full computation alphabet is

    C_strict = {(a,z) : a in C, 0 <= z < q},
    k_c,strict = q(2m+T+2),
    d0,strict = 4 q^2(2m+T+2).                                (11.1)

A two-site pair ((a,z),(b,z')) is allowed if and only if (a,b) is allowed in Section 2, for every z,z'. Synchronization tests a=M, independently of z. In particular, neither workspace data nor Bell entanglement is a legality condition. The legal clock-pattern proof applies unchanged to the a factors, while the z factors are quantum workspace.

Every propagation gate acts on the z registers, tensor identity on both physical p registers. The clock-pattern transitions and the chi_i chi_(i+1) marker guards are unchanged. Thus, as an operator identity on the whole Hilbert space,

    [H_prop, F_(p,i)] = 0

for every site i and every physical-register operator F. Since marker projectors also commute with propagation, the same is true of a port operator conditioned on an adjacent actual selector marker. Illegal workspace labels can no longer violate port idleness.

### Explicit bond prescriptions for the short program

Take q=2. Keep the twelve clock letters, fourteen allowed clock pairs, and two transition patterns from Section 6. The actual alphabet has 24 letters (a,z). The two gates are SWAP of z_left,z_right and identity. Their propagation terms are (3.1), with those gates acting on z and identity on p.

On a bond in local register order (z_left,p_left,z_right,p_right), the three energy matrices are

    A1 H0: 3 Delta(1-P_Omega on p_left,z_right),
    H0 R3: (3g/2)(X_(z,left) X_(z,right)-Y_(z,left) Y_(z,right)),
    A1 H2: (3g/2)(X_(p,left) X_(z,right)+Y_(p,left) Y_(z,right)).  (11.2)

They are controlled only by the displayed clock-letter pair and by chi_left chi_right. The unused registers carry identity in these matrices. The Bell partner is z at site 2; the logical qubit is z at site 3. The left and right ports are p at sites 1 and 4.

To remove spectator degeneracies, set Q_z=I-|0><0| on z and Q_p=I-|0><0| on p. Add the positive onsite term

    H_pin,site = nu p_selector# (Q_z+Q_p)
        + nu chi [ sum_(a: position(a) in {1,m}) |a><a| Q_z
                 + sum_(a: 2 <= position(a) <= m-1) |a><a| Q_p ]. (11.3)

The position of H_t is j_t. M has no position and is omitted from the two sums. The first term replaces the compact marker-data pinning term. All these are finite energy penalties.

For the short example this pins z at sites 1 and 4, p at sites 2 and 3, and both registers on actual markers. The active four qubits are (p_1,z_2,z_3,p_4). The pinning is independent of Bell preparation.

The full bond is now (3.4), with the lifted clock and energy operators just specified, marker pinning replaced by (11.3), and the same C S_head-rho p_selector# counterterm. Distribute all onsite terms equally between adjacent bonds. This is an explicit nearest-neighbor construction on the 24-letter computation alphabet.

### Spectral relation to the compact version

Each propagation transition preserves the position index attached to each physical site: H_t and the A or R letter replacing it have the same index. Thus the full string of position indices, including the locations of M, is conserved. This statement holds in illegal sectors as well.

Within a fixed index string, identify the compact data register with p on indices 1,m and with z on indices 2 through m-1. At a computation M with no active interactions, either factor can be chosen. Every term in the compact Hamiltonian acts on these active registers in the strict version. The other register is unused. Every pin in (11.3) acts on an unused register, except that the actual-marker pin on the active register reproduces the compact marker-data pin.

Consequently each such block of the strict Hamiltonian is unitarily equivalent to the corresponding compact block plus nonnegative, decoupled unused-register pins. Setting all unused registers to |0> attains the compact spectrum. This proves that the entire legality and bad-cell energy comparison transfers, not just the legal ground energy. Excitations of a pinned unused register cost nu. Unpinned unused factors at a misplaced computation M can add degeneracy only inside an already penalized illegal sector.

In a good cell, the initialized active Hamiltonian is exactly the compact one, tensor the fixed unused-register state. Its ideal code compression, actual finite-J compression, and actual ground energy are therefore the same as those in Section 10. The good-cell gap gains an additional possible scale nu, with no other change. The bad-cell lower bound remains kappa/2-rho, and the selector-preservation proof is unchanged.

### Numerical check of the strict physical bond

`build_strict_bond_c10.py` constructs this bond using

    k_c=24, q=2, d0=192, nu=20,

with the other parameters from Section 10. Its two-site matrix is 36,864 by 36,864 and has 42,168 nonzero entries. Hermiticity and same-bond selector commutation have zero measured error. The propagation commutators with both X_p and Y_p on either endpoint are zero on the full two-site matrix, with no legality projection.

Independently projecting that actual bond onto the three legal clock states and the sixteen active-register states gives the compact 48-dimensional matrix to 1.12e-16. The actual ground energy is -0.076582449455, as before. These results are in `strict_bond_verification.json`.

The enlarged illegal quantum Hilbert space was not directly diagonalized. The block equivalence above transfers the exhaustive compact-sector calculation: its minimum illegal energy with correct fixed endpoint markers is still 20 before the counterterm and 19.4 afterward. The all-sector shared-boundary bound remains 9.4. This transfer is an analytic argument, not an additional ED claim.

The strict version repairs all-sectors port idleness. It does not repair the growing alphabet or the finite-J exact-history obstruction.

The equal-dimension register duplication used here also has a dimension cost for Route B. The supplied flavor count gives c=q_+q_-<=q^2/4, whereas this implementation has d0=4q^2(2m+T+2). Thus this particular duplicated family cannot produce a Route B violation by increasing q. It is a literal port-idle track and a finite-program verification, not a dimension-efficient Route B construction. This cost is not a no-go theorem for other encodings; fixed qubit ports in Route A also do not incur a growing-q cost.

## 12. Dimension and unresolved claims

For the implemented short program with port-idle propagation on every sector,

    k_c=24, q=2, d0=4*24*2=192.

The compact legal-sector-only version has k_c=12 and d0=96.

Its left Route B example has r=2. It is a verification gadget, not a large-central-charge construction. For the general hard-coded track with the stronger port condition, the proved formula is d0=4q^2(2m+T+2), which still grows with the program.

The completed parts are an explicit finite bond, exact path histories before incompatible energies are added, exact ideal port compression, exact energetic Bell preparation in the frustration-free example, an all-illegal-sector penalty bound for this alphabet, and a selector-preservation comparison using explicit couplings.

The missing parts are (i) a fixed-alphabet program/counter decoder with the required all-sector proof, (ii) an exact finite-J realization of the desired Route B ground history code, or an explicitly approximate reformulation with its downstream error control, and (iii) a tested two-ended or Route A compiler giving an N-independent d0. No conjecture violation follows from the present partial construction.

The two next technical tasks are separate. A uniform machine must replace the position and time labels and prove that no unchecked malformed computation component survives. An exact-history version must also supply time energies satisfying (5.2), usually by an additional frustration-free mechanism. Increasing J alone does not prove either statement.

## Files produced

- `PROGRESS.md` records the work and the change from the exact-energy to explicit-upper-bound counterterm.
- `RESULT.md` is this derivation and scope statement.
- `build_bond_c10.py` constructs and checks the compact sparse physical bond.
- `build_strict_bond_c10.py` constructs and checks the strict all-sectors port-idle bond.
- `verify_c10.py` checks legality, clock energies, Bell preparation, and Route B compression.
- `verify_full_cell_c10.py` diagonalizes every quantum computation component of the small cell and independently checks the physical-bond projection.
- `verification.json`, `bond_verification.json`, `full_cell_verification.json`, and `strict_bond_verification.json` contain the numerical results.
- `verification_run.txt` is the saved output of the updated small-program verification.
