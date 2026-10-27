# C9 result

## Verdict

C9 remains unresolved. I have not constructed the required translation-invariant nearest-neighbor Hamiltonian with an explicit fixed local dimension d0, and I have not proved that such a Hamiltonian is impossible.

There is a precise obstruction to one proposed completion. If the microscopic symmetry is a product of on-cell SU(N) actions, the code in each cell is fundamental, and neighboring cells interact only through q-dimensional boundary ports, then N > q forces every logical interaction to be scalar. This applies to arbitrary port interactions, not just the port swap, under the symmetry assumptions stated below. Thus imposing this exact microscopic SU(N) symmetry cannot protect the large-N exchange gadget.

There is also a positive finite-dimensional result. Two independent qubit ports suffice to synthesize any Hermitian logical nearest-neighbor bond at second order, with an exactly degenerate N-dimensional code, vanishing first order, and a finite cell gap. The construction includes antiferromagnetic swap with a positive coefficient. It is an arbitrary matrix Hamiltonian on a cell, not a local implementation within that cell. A conservative dimension bound is

    D = 4 ceil((6 N^3 + 13 N)/4).

The remaining microscopic problem can therefore be stated more broadly and more sharply. Realize this two-port spectral gadget by a finite programmable cell, preserving the boundary matrix elements throughout the relevant excited band. Matching its spectrum alone is insufficient. The remaining infrared problem is independent: either protect the target beyond perturbative order, or simulate a lattice model in an open stable critical region.

The first-order Majorana alternative admits explicit joint two-port encodings. I construct and check small examples below. It requires different left and right Clifford operators, an even number of channels for the direct on-cell Fock-space construction, and protection against higher-order interactions. A maximal ungraded three-channel example violates logical fermion parity already at second order.

All statements below are derived closed book. Only files created during this task were read. The two ingredients supplied in the question are accepted as premises. No external universality theorem or CFT stability theorem is silently assumed.

## 1. What a cell simulator actually has to preserve

Let G map the logical code into an abstract cell, let P = GG†, and let H_A have kernel P and positive excited energies E_e. For a physical boundary operator F_a, define

    S_a = F_a G,
    M_e[a,b] = S_a† P_e S_b.

Here P_e is an excited spectral projector. In the zero-first-order case, the second-order bond is

    H_eff^(2)
      = -lambda^2 sum_(e,f,a,b)
          M_e^R[a,b] tensor M_f^L[a,b] / (E_e + E_f).             (1)

This follows directly by inserting the spectral projectors in

    -P_AB V (H_A + H_B)^(-1)|_(P_AB perpendicular) V P_AB.

The a,b ordering in both factors is the same because the physical coupling is sum_a F_a tensor F_a. Off-diagonal blocks can be non-Hermitian individually; the full sum is Hermitian. The constructions below use block-diagonal M_e, avoiding this bookkeeping issue.

A simulator must preserve the spectral weights M_e, or the equivalent port resolvent functions

    G† F_a Q (H_A + z)^(-1) Q F_b G,

at the energies relevant to (1). It is not enough to approximate eigenvalues, ground-state energy, or ground-space dimension.

For example, consider a two-time history isometry

    W|psi> = (|0> |psi> + |1> X|psi>)/sqrt(2).

It satisfies

    W† (1 tensor Z) W = (Z + XZX)/2 = 0.

One can assign any desired spectrum on ran W and put its orthogonal complement at arbitrarily high energy. Nevertheless, the physical Z observable has lost the target matrix elements. This is an exact finite-dimensional example of why a spectral simulation statement does not imply port simulation.

There is also a geometry issue. A one-port cell cannot use that same physical boundary site at both ends of a long one-dimensional cell. Two simultaneous geometric ports are needed. The next construction repairs this at the abstract-cell level.

## 2. Any logical bond through two qubit ports, at second order

### 2.1 Code and two orthogonal syndrome families

The construction works for any q >= 2. Later set q = 2 and r = q^2 - 1 = 3.

Take two q-dimensional ports, L and R, and an auxiliary space containing two q-dimensional registers and an N-dimensional logical register. Define

    G|i> = (1/q) sum_(a,b=1)^q
             |a>_L |b>_R |a,b,i>_aux.                           (2)

The initial ambient dimension is q^4 N. More unused auxiliary states can be added.

For traceless Hermitian F_a with tr(F_a F_b) = q delta_ab,

    G† F_a^L G = G† F_a^R G = 0,
    G† F_a^L F_b^L G = delta_ab 1_N,
    G† F_a^R F_b^R G = delta_ab 1_N,
    G† F_a^L F_b^R G = 0.                                      (3)

These identities follow by taking the normalized trace in each independent maximally entangled pair in (2).

Thus the columns of

    S = [F_1^L G, ..., F_r^L G, F_1^R G, ..., F_r^R G]

are an orthonormal frame of dimension m_syn = 2 r N in the orthogonal complement of G. In particular, the two ports can have independently specified positive spectral weights on their syndrome spaces.

### 2.2 Positive spectral-weight construction

Write the desired dimensionless logical bond as

    K_target = sum_(i,j=1)^m K_ij A_i tensor B_j,                 (4)

where A_i and B_j are Hermitian and K is real. One may use complete Hermitian bases, including the identity, with m = N^2. A smaller operator-Schmidt decomposition often reduces m.

Choose k = m + 2 distinct positive dimensionless numbers e_e and probabilities p_e > 0 with sum_e p_e = 1. The actual excited energies will be Delta e_e. Define

    C_ef = 1/(e_e + e_f).

C is strictly positive definite. For a nonzero real vector x,

    x^T C x = integral_0^infinity (sum_e x_e exp(-e_e t))^2 dt > 0.

The strict inequality follows from linear independence of exponentials with distinct exponents.

Consider the real subspace

    V = {z : 1^T z = 0 and (C p)^T z = 0}.

It has dimension at least k - 2 = m. Choose a k by m matrix Z whose columns lie in V and are orthonormal in the C metric. Explicitly, if B spans any m-dimensional subspace of V, set

    Z = B (B^T C B)^(-1/2).

Then

    1^T Z = 0,    p^T C Z = 0,    Z^T C Z = 1_m.                (5)

Only one syndrome label, a = 1, needs a nontrivial logical weight. On the rN-dimensional right and left syndrome spaces, respectively, choose block-diagonal weights

    M_e^R[1,1] = p_e 1_N + rho sum_i Z_ei A_i,
    M_e^L[1,1] = p_e 1_N - rho sum_j (Z K)_ej B_j,

    M_e^R[a,a] = M_e^L[a,a] = p_e 1_N       for a > 1,
    M_e^R[a,b] = M_e^L[a,b] = 0             for a != b.          (6)

Their sums over e are the identity by the first relation in (5).

Positivity is explicit. Let

    a_e = ||sum_i Z_ei A_i||,
    b_e = ||sum_j (Z K)_ej B_j||.

Choose any rho > 0 satisfying

    rho max(a_e,b_e) <= p_e/2     for every e.                  (7)

A zero denominator imposes no restriction. Every weight is then bounded below by p_e/2 on its modified block and by p_e on all other blocks. All choices are finite for each N, even if rho is very small.

### 2.3 Realizing the weights by an actual finite matrix H_A

Combine the left and right weights into

    M_e = M_e^L direct_sum M_e^R.

They form a positive operator-valued measure on the m_syn-dimensional syndrome source. Define the isometry

    W x = direct_sum_e sqrt(M_e) x

into a space of dimension k m_syn. The equality W†W = 1 follows from sum_e M_e = 1.

Suppose the physical cell dimension obeys

    D - N >= k m_syn.

There exists an isometric embedding J of the dilation space into G perpendicular such that J W = S. To construct it, extend W to an orthonormal basis of its dilation space, extend S to an orthonormal basis of a k m_syn-dimensional subspace of G perpendicular, and map the first basis to the second.

Let Pi_e project onto the e-th coordinate block of the dilation. Set

    P_e = J Pi_e J†,
    H_A = Delta sum_e e_e P_e
          + Delta e_* (1 - GG† - JJ†),                         (8)

where e_* > 0. Then

    ker H_A = ran G,
    S† P_e S = W† Pi_e W = M_e,
    gap(H_A) = Delta min({e_e}, e_*)

when the corresponding eigenspaces are present. This proves the required finite-dimensional dilation directly.

For q = 2, m_syn = 6N. Taking m = N^2, k = N^2 + 2, and keeping D divisible by four so it factors into the two ports and an auxiliary space gives

    D = 4 ceil(max(16N, (6k+1)N)/4)
      = 4 ceil((6N^3 + 13N)/4),       N >= 1.                   (9)

This is a sufficient bound, not an optimized dimension.

### 2.4 Effective bond and sign

Couple the right port of one cell to the left port of the next by

    lambda v,    v = sum_(a=1)^r F_a^R tensor F_a^L.

The first order vanishes by (3). Inserting (6) into (1), the scalar contribution is

    -(lambda^2/Delta) r p^T C p.

Each term linear in rho vanishes by p^T C Z = 0 and symmetry of C. The remaining quadratic contribution is

    +(lambda^2/Delta) rho^2
        sum_(i,j) (Z^T C Z K)_ij A_i tensor B_j
      = (lambda^2/Delta) rho^2 K_target.                       (10)

Therefore

    H_eff^(2) = (lambda^2/Delta)
                  [-r p^T C p 1 + rho^2 K_target].             (11)

The result has the correct energy units, lambda^2/Delta. It also remains negative semidefinite as required by second-order perturbation theory, because its original expression is minus a positive resolvent quadratic form. Positivity of (6) ensures that the scalar term dominates any positive eigenvalue of rho^2 K_target. For the numerical example below, its eigenvalues are -0.42463235 plus or minus 0.00063633.

The sign in front of K_target is positive. Choosing K_target = P_N gives antiferromagnetic exchange. Choosing its negative gives the opposite sign. Arbitrary on-cell terms can also be included in (4) using identity factors; no order-lambda^2 splitting of the isolated code is required.

For a chain, two different bonds give zero cross terms at second order. If the bonds touch, each outer cell is acted on only once, and its compression of F_a is zero. If they do not touch, the same argument applies to all four outer cells. Thus (11) is the entire nonconstant second-order chain Hamiltonian, summed over neighboring cells.

This result is compatible with the symmetry obstruction below. The microscopic cell weights are not SU(N) covariant. SU(N) symmetry for K_target = P_N is an accidental property of the second-order coefficient.

### 2.5 Scales and limitations

With q = 2 and Pauli F_a, ||v|| = 3. The isolated-cell code is exactly degenerate, and its gap is fixed by (8). Define

    J = rho^2 lambda^2/Delta.

For a fixed finite cluster, ordinary degenerate perturbation theory gives a third-order bound of the form

    ||remainder|| <= C |lambda|^3 ||v||^3 / Delta^2,

provided |lambda| ||v|| is sufficiently small compared with the cell gap. The constant depends on the gap convention and cluster size. To make this small relative to J, a sufficient hierarchy is

    |lambda|/Delta << rho^2/(C ||v||^3).                        (12)

The small rho produced by (7) matters in this inequality. It cannot be omitted when comparing higher orders to the desired exchange.

In an infinite chain these estimates must be organized per connected cluster or per local interaction term. A bound requiring |lambda| L << Delta is not acceptable. Nor does a finite-cluster calculation prove an exact global low-energy band separated from all excited bands: the logical band has an extensive bandwidth. A completed microscopic construction would need a local effective-Hamiltonian argument and a bound on the cost of creating machinery excitations.

Equation (11) is proved only to second order. No exact higher-order symmetry or criticality follows from it.

## 3. An exact SU(N) protection obstruction

### 3.1 Statement and assumptions

Suppose a cell has a continuous unitary representation U(g) of SU(N), and assume all of the following.

1. H_A commutes with U(g).
2. The code transforms as the irreducible fundamental representation of SU(N).
3. An inter-cell bond acts only on a q_R-dimensional right port and a q_L-dimensional left port.
4. The bond is invariant under U(g) tensor U(g).
5. N > max(q_R,q_L).

Then the interacting port operators commute with U(g) on each cell separately. The full chain has an independent SU(N) action in each cell. Every effective interaction on the tensor product of fundamental code spaces is scalar, at every perturbative order.

For a microscopic nearest-neighbor realization, q must count the full interacting boundary Hilbert space, including any control states used by the bond. It is at most d0 for a single-site boundary. Thus the obstruction still applies in the desired regime N > d0.

The assumptions are essential. This is not a no-go theorem for emergent SU(N), for a symmetry whose microscopic action is not a product over cells, for a different protecting group, or for a stable critical theory requiring no such symmetry.

### 3.2 Proof

Take a minimal operator-Schmidt expansion of the bond,

    v = sum_(mu=1)^s A_mu tensor B_mu,

with both families linearly independent. The left and right operator-Schmidt spaces are intrinsic to v. Invariance under U(g) tensor U(g) therefore implies

    Ad_U(g) span{A_mu} = span{A_mu},
    Ad_U(g) span{B_mu} = span{B_mu}.                             (13)

Let A_R be the unital star algebra generated by the A_mu. It is a subalgebra of M_(q_R), tensored with the identity on the rest of the cell. Equation (13) implies that conjugation by U(g) induces a continuous automorphism of A_R. The same holds for A_L.

Every finite-dimensional unital star subalgebra of M_q has the form

    A = direct_sum_alpha (M_(n_alpha) tensor 1_(m_alpha)),
    sum_alpha n_alpha m_alpha = q.

A connected automorphism fixes its central blocks and acts by inner automorphisms on each full matrix block. Its Lie algebra is therefore a sum of su(n_alpha), of dimension

    sum_alpha (n_alpha^2 - 1) <= q^2 - 1.                       (14)

The derivative of the SU(N) action is a Lie algebra homomorphism from su(N) into this algebra. Since su(N) is simple, the kernel is either zero or all of su(N). A nonzero map would be injective, which is impossible when

    N^2 - 1 > q^2 - 1.

Thus the derivative vanishes. Connectedness of SU(N) implies that its entire action on A is trivial. Hence every A_mu and B_mu commutes with the on-cell symmetry.

Each bond consequently commutes with independent transformations of its two cells, and H_A does too. The full Hamiltonian commutes with SU(N) independently at every cell. The tensor product of fundamental code spaces is irreducible under this product group. By Schur's lemma, an operator on it commuting with the product group is scalar.

For a fixed finite number of cells and sufficiently small coupling, the spectral subspace continuing the unperturbed code has the same representation content. Its Hamiltonian is scalar. Equivalently, every coefficient of symmetry-respecting degenerate perturbation theory is scalar. This proves the obstruction without assuming convergence of an infinite-chain expansion.

For the full port swap, the argument is especially direct. Its Schmidt space is all traceless q by q matrices. U(g) would induce a projective q-dimensional representation of SU(N). The same dimension argument makes it trivial for N > q.

### 3.3 A smaller subgroup is not an automatic repair

Here is a related restricted obstruction. Suppose the port is a qubit transforming as spin 1/2, the code is a spin-j irreducible SU(2) representation, and H_A is exactly SU(2) invariant. Each Pauli port operator transforms as spin 1. Therefore every operator

    G† F_a P_e F_b G

transforms in spin 1 tensor spin 1, containing tensor ranks only 0, 1, and 2. The second-order interaction can use only those ranks on each logical cell.

On a spin-j code, the swap has nonzero operator-Schmidt components in every tensor rank from 0 through 2j. This follows by choosing a complete orthonormal operator basis organized by SU(2) tensor rank in the identity

    P_N = sum_mu T_mu tensor T_mu†,

where tr(T_mu† T_nu) = delta_mu,nu and N = 2j + 1. All ranks must appear to span End(C^N).

For N >= 4, rank 3 is already needed. Exact SU(2) symmetry with a spin-1/2 port therefore cannot yield P_N at second order on such a code. Other subgroup representations require their own analysis; this example is not a general subgroup no-go theorem.

### 3.4 A separate first-order rank bound

At first order, any encoded interaction through q-dimensional ports has operator-Schmidt rank at most q^2. Local compression cannot increase operator-Schmidt rank.

For b != 0, the bond a 1 + b P_N has rank N^2, except at the scalar cancellation a = -b/N, where it has rank N^2 - 1. Both exceed q^2 when N > q. Thus direct first-order SU(N) exchange through a smaller port is impossible, regardless of symmetry.

Second order escapes this bound because the arbitrarily many excited spectral weights supply additional logical operators. Section 2 exhibits the mechanism explicitly.

## 4. First-order Majorana route

### 4.1 Correct left and right boundary operators

Let r = 2m be even and use a code of dimension N = 2^m. Let gamma_a be Hermitian Clifford generators and P be their fermion parity,

    {gamma_a,gamma_b} = 2 delta_ab,
    P^2 = 1,    P gamma_a P = -gamma_a.

Define

    eta_a = i gamma_a P.

Each eta_a is Hermitian and squares to one. On a chain of logical cells, define global Majoranas by the Jordan-Wigner transformation

    Gamma_(j,a) = (product_(k<j) P_k) gamma_(j,a).

Direct multiplication gives

    i Gamma_(j,a) Gamma_(j+1,a)
      = (i gamma_(j,a) P_j) gamma_(j+1,a)
      = eta_(j,a) gamma_(j+1,a).                               (15)

Consequently the required two boundary compressions are

    G† F_a^L G = s gamma_a,
    G† F_a^R G = s eta_a.                                      (16)

Using gamma_a at both ends does not give (15). The missing parity factor changes the interacting spin Hamiltonian.

With (16), first order is

    lambda s^2 sum_(j,a) i Gamma_(j,a) Gamma_(j+1,a).             (17)

For an even number of cells on a ring one also fixes the usual global parity-dependent fermion boundary condition. Each independent uniform Majorana chain has gapless points at momenta 0 and pi, which form one nonchiral massless Majorana field. Its central charge is 1/2. Thus the exact free Hamiltonian (17) has c = r/2. This statement concerns (17), not the full gadget Hamiltonian.

### 4.2 A joint completely positive encoding exists

It is necessary to realize (16) with two commuting physical port algebras in the same cell. Constructing the two compressions separately does not establish that.

Choose r <= q^2 - 1 orthonormal traceless Hermitian port matrices F_a. Let Q = q^2 and define matrices on the two-port space and on the code by

    B_a = F_a tensor 1,       T_a = gamma_a,
    B_(r+a) = 1 tensor F_a,   T_(r+a) = eta_a.

Define the map from physical port observables to logical observables,

    Phi(A) = tr(A)/Q 1_N + (s/Q) sum_(nu=1)^(2r) tr(B_nu A) T_nu.
                                                                    (18)

It is unital and has the desired compressions. Its Choi matrix is

    J_Phi = (1/Q) [1_(QN) + s sum_nu B_nu^T tensor T_nu].        (19)

All T_nu have norm one. Hence J_Phi is positive whenever

    2 s sum_a ||F_a|| <= 1.                                    (20)

Strict inequality gives a positive definite Choi matrix. This proves complete positivity at a nonzero scale for every finite r and q.

For completeness, diagonalize J_Phi = sum_mu j_mu |v_mu><v_mu| and define

    G_(i,mu),alpha = sqrt(j_mu) conjugate(v_mu(i,alpha)).

Then

    G† (A tensor 1_env) G = Phi(A),    G†G = Phi(1) = 1.

Thus this is an actual isometric encoding, with cell dimension

    D = Q rank(J_Phi) <= q^4 N.                                (21)

It also satisfies

    Phi(F_a tensor F_b) = 0.                                   (22)

Take H_A = Delta(1-GG†). Its code is exactly degenerate and its gap is exactly Delta. This solves the abstract joint-encoding problem.

### 4.3 Odd r leaves unused Majoranas in this construction

A graded representation of an odd number r of Clifford generators requires an additional generator. A convenient representation uses r+1 generators on a code of dimension 2^((r+1)/2), with P the parity of that even Clifford algebra.

After the Jordan-Wigner transformation, the additional global Majorana in each cell commutes with (17), because it anticommutes with both factors of every hopping bilinear. For an even number M of cells, these M unused Majoranas generate a Clifford algebra and force at least 2^(M/2)-fold degeneracy of every level.

This is not a pure r-channel critical ladder with gapped machinery. The direct construction must use even r, or supply and analyze an additional mechanism that removes the unused modes. It does not rule out every possible odd-channel construction with constraints or a different unit cell.

For the direct even-r construction,

    r_max = 2 floor((q^2-1)/2),
    c_first_order = floor((q^2-1)/2).                           (23)

This is finite at fixed q. Merely increasing the cell length cannot make the central charge of this particular first-order free ladder unbounded.

### 4.4 Exact physical parity reduces the available operator count

Logical parity at first order need not be a microscopic symmetry. To implement it exactly through a physical port grading, let the port grading have q_+ positive and q_- negative eigenvalues. Its odd Hermitian operators are off-diagonal between these eigenspaces, a real vector space of dimension

    r_odd = 2 q_+ q_-.

Using only an orthonormal basis of these odd matrices in (18) gives

    Phi((u tensor u) A (u tensor u)) = P Phi(A) P,

where u is the port grading. A covariant Stinespring encoding can be chosen, so the projector Hamiltonian H_A and the inter-cell bond have exact global parity.

This covariance can be seen directly from (19). Its Choi matrix commutes with (u tensor u)^T tensor P. Choose the Choi eigenvectors to diagonalize that involution as well. Their eigenvalues define an auxiliary involution U_env. The isometry then satisfies

    [(u tensor u) tensor U_env] G = G P.

Thus GG†, and hence H_A, commute with that physical cell involution. The product of two odd port operators on a bond is invariant.

The maximum first-order central charge with this grading is therefore

    c = q_+ q_- <= floor(q^2/4).                                (24)

One cannot make every traceless port matrix odd under an automorphism of M_q. If an automorphism sent F and G to their negatives, it would fix [F,G] by multiplicativity, whereas inverting every traceless matrix would negate [F,G]. A nonzero commutator gives a contradiction. This explains why the ungraded r = q^2-1 construction generally fails to have exact physical parity.

Exact parity still allows many quartic interactions. It does not by itself prove stability of a large free Majorana ladder.

### 4.5 An explicit non-free second-order correction

The numerical r = 2 cell below uses q = 2, F_1 = X, F_2 = Y, and maps the unused Z direction to zero. Let

    T = sum_(a=1)^2 eta_a tensor gamma_a.

For either port, (18) gives

    C_ab = G† F_a F_b G = delta_ab 1,
    K_ab^R = C_ab - s^2 eta_a eta_b,
    K_ab^L = C_ab - s^2 gamma_a gamma_b.

Intermediate states with one excited cell have denominator Delta, while those with both cells excited have denominator 2 Delta. Therefore the coefficient of lambda^2 is

    H_2 = -(1/Delta) sum_ab [K_ab^R tensor s^2 gamma_a gamma_b
                            + s^2 eta_a eta_b tensor K_ab^L]
          -(1/(2 Delta)) sum_ab K_ab^R tensor K_ab^L.

Expanding gives

    H_2 = -(r(1+2s^2)/(2 Delta)) 1
          + (3s^4/(2 Delta)) T^2,       r=2.                   (25)

For gamma_1 = X, gamma_2 = Y and P = Z,

    eta_1 = Y,    eta_2 = -X,
    T = Y tensor X - X tensor Y,
    T^2 = 2(1 - Z tensor Z).

There is a nonconstant density interaction, not merely a change in the hopping velocity.

For the ungraded r = 3 example, the Pauli multiplication law gives

    C_ab^R = delta_ab 1 + i s epsilon_abc eta_c,
    C_ab^L = delta_ab 1 + i s epsilon_abc gamma_c.

The analogous expansion contains the parity-odd contribution

    H_(2,odd) = -(i s^3/(2 Delta)) sum_(a,b,c) epsilon_abc
                    [eta_c tensor gamma_a gamma_b
                     + eta_a eta_b tensor gamma_c].            (26)

The other terms are parity even. The operator in (26) is nonzero in the explicit representation used below. It breaks P tensor P. A logical parity grading stated only for the first-order code is therefore insufficient.

On a chain, the first-order construction also has second-order processes involving adjacent bonds. Unlike Section 2, each outer port now has a nonzero ground-space compression. The two-cell calculation does not enumerate those longer-range terms and makes no claim that they vanish.

## 5. Where a programmable-cell construction is still missing

### 5.1 Capacity is not a Hamiltonian construction

The spontaneous-cell premise supplies large registers at fixed physical dimension. For example, a spectator qubit supplies exponentially many states as the cell length grows. This proves a capacity statement. It does not show that the restriction of a fixed nearest-neighbor interaction to that register can equal the matrices (8) or (21).

In particular, D in (9) is the dimension of an abstract cell, not the dimension of one physical site. Setting d0 = 8 on the basis of spectator capacity would not implement H_A. No such identification is made here.

The same distinction applies to a program written along a cell. A fixed alphabet can store an arbitrarily long program, but a translation-invariant Hamiltonian must select the intended program and exclude other program sectors. Supplying a favorable sector by hand does not show that it is a ground sector.

A useful elementary restriction is that a purely classical diagonal translation-invariant nearest-neighbor Hamiltonian on a fixed alphabet cannot uniquely select arbitrarily large periods by changing its pair energies. A periodic configuration is a closed walk in a finite weighted graph. Decomposing a closed walk into simple cycles shows that its mean energy is at least the minimum simple-cycle mean. Some minimizing simple cycle has length at most the alphabet size. Longer minimizing words can occur through degeneracy, but they cannot exclude all shorter minimizing cycles. The quantum walker in the supplied premise evades this classical argument. It does not automatically select a computational program on an additional track.

### 5.2 What a finite Feynman-Kitaev implementation would buy

The following is a conditional calculation on an already legal, initialized history sector. It is not a construction of that sector by a fixed translation-invariant nearest-neighbor rule.

Let U_t denote the accumulated circuit at clock time t, with t = 0,...,T. The history isometry is

    W|psi> = (1/sqrt(T+1)) sum_t |t> U_t|psi>.

A path propagation Hamiltonian with edge terms

    (J_c/2) (|t> tensor 1 - |t+1> tensor U_(t+1,t))
              (Hermitian conjugate)

is unitarily equivalent to a path Laplacian. Its first nonzero energy is

    delta_clock = J_c [1 - cos(pi/(T+1))]
                ~ J_c pi^2/[2(T+1)^2].                        (27)

Suppose a circuit, using initialized ancillas, implements a positive-output test whose compressed output projector is H_A/B, with B >= ||H_A||. The abstract unitary dilation is explicit. With A = sqrt(1-H_A/B) and C = sqrt(H_A/B), the block matrix U = [[A,-C],[C,A]] is unitary because A and C commute and A^2+C^2=1. Acting on an initialized zero ancilla, its output-one probability operator is H_A/B. Compiling this unitary with local gates to a prescribed error is a separate step.

Applying (T+1) B times the output projector at the designated clock time gives H_A after compression into the uniform history. Its uncompressed norm is at most (T+1)B. If it is small compared with delta_clock, the error from mixing with clock excitations is bounded in a fixed finite cell by a constant times

    (T+1)^2 B^2 / delta_clock
      = O(B^2 (T+1)^4/J_c).                                   (28)

Thus one can make that error at most epsilon_H by taking, for example,

    J_c >> B^2 (T+1)^4/epsilon_H,
    J_c >> B (T+1)^3.                                         (29)

All these scales may be finite functions of N if T depends only on the cell, not on L. In a positive propagation-plus-test construction, zero-energy histories of target kernel states remain exact ground states in the legal history sector. Initialization and illegal-sector penalties still need their own proof.

Boundary observables need a separate argument. Suppose the physical port coincides with the target port during T_idle+1 idle clock times, and there are T_act active times. Then

    ||W† O_physical W - O_target||
      <= 2 ||O|| T_act/(T+1).                                 (30)

This is just the triangle inequality applied to the average of U_t† O U_t. Long idling can therefore make port compression accurate. It generally does not make it exact. In particular, it need not preserve zero first order exactly.

Equations (27)-(30) explain why arbitrarily large but finite energy scales could solve the finite-cell accuracy problem. They do not provide the missing translation-invariant clock and program.

### 5.3 The unresolved local construction requirements

An implementation of this approach must give, with a fixed alphabet and fixed local dimension:

- A unique intended computational sector inside each selected cell, including initialized workspace and the intended program or an algorithm that generates it.
- Local transition rules whose legal histories implement the desired circuit, with a bound on every clock and workspace gap.
- Penalties excluding malformed configurations, extra heads, disconnected computation orbits, and incorrect program histories. In particular, an orbit that never reaches an initialization test must not supply an unintended low-energy sector.
- Two accessible geometric boundary ports, with the observable estimate valid throughout the entire target band used in (1), not only on the ground code.
- Preservation of marker sectors, or a new proof of cell selection if marker positions are allowed to move.
- Bounds uniform in the number of cells, with couplings depending on N but not on L.

I do not provide such local rules here. Invoking a general universal-simulation construction without these properties would assume the main missing step.

A global history clock whose duration grows with L does not solve this finite-cell problem. At fixed propagation coupling its gap falls as its inverse duration squared. Increasing its coupling with L is prohibited by the question, and allowing those low-energy clock modes requires a new infrared analysis. The computation must be confined and reset within cells of fixed length b(N).

## 6. Error scales and the infrared obstruction

### 6.1 Sufficient finite-cell accuracy conditions

Let Delta be a lower bound on the intended gadget gap, B its upper spectral scale, Lambda the separation to additional simulator states above the entire relevant band, and lambda the port coupling. Assume, conditionally, that a simulator provides

    Hamiltonian error <= epsilon_H,
    boundary-observable error <= epsilon_F,

on that full band, with correctly identified code and excited states. Here epsilon_H has energy units and epsilon_F is dimensionless. Code-projector rotation can be absorbed into epsilon_F; generically it also contributes O(epsilon_H/Delta).

For fixed N, a schematic local perturbative bound takes the form

    ||delta h_eff|| <= C_N [ epsilon_H
                           + |lambda| epsilon_F
                           + lambda^2 epsilon_H/Delta^2
                           + lambda^2 epsilon_F/Delta
                           + lambda^2/Lambda
                           + |lambda|^3/Delta^2 ].             (31)

Constants include port norms, finite coordination, gap ratios, and the chosen target normalization. The term lambda^2/Lambda accounts for port transitions to simulator states outside the target band, even if their matrix elements are order one. This term is why simply putting the clock above the logical exchange J can be inadequate; it must be high enough relative to the full gadget process.

For the construction in Section 2, J = rho^2 lambda^2/Delta. A sufficient hierarchy for a relative local error at most a small eta is

    epsilon_H << eta J/C_N,
    epsilon_F << eta J/(C_N |lambda|),
    |lambda|/Delta << eta rho^2/C_N,
    Lambda >> C_N Delta/(eta rho^2).                           (32)

Also require |lambda| times the physical bond norm much smaller than the marker-defect, walker, clock, and workspace excitation costs. These are finite-cell or local-cluster requirements, not a proof of criticality.

If each simulator ground state remains exactly at zero by a frustration-free construction, its code splitting can vanish exactly. If not, the first inequality in (32) is needed to keep the splitting below the exchange. Approximate vanishing of port compression requires the second inequality.

### 6.2 Finite smallness does not imply asymptotic criticality

Let a permitted continuum perturbation have scaling dimension x < 2 and bare dimensionless coefficient g. Linearized renormalization gives

    dg/dell = (2-x) g + O(g^2),
    xi/a ~ |g|^(-1/(2-x)).                                    (33)

For a marginally relevant perturbation, with dg/dell = A g^2 + higher orders and A g > 0,

    xi/a ~ exp[1/(A g)].                                       (34)

These are field-theory estimates, not lattice gap theorems. They suffice to show the problem with an accuracy argument: for any fixed nonzero g, making it very small can give an enormous crossover length without giving an infinite one.

The order of limits in the question is fixed N and fixed microscopic couplings first, followed by L -> infinity. One cannot use epsilon(N) > 0 or lambda(N) > 0 as a substitute for a stability argument, however small they are. Nor can a finite-size numerical fit distinguish an arbitrarily long crossover from the asymptotic critical theory.

### 6.3 What translation does and does not protect at SU(N)_1

Under the usual fundamental-chain identification, one-cell translation acts by the center on the SU(N)_1 primaries. The primary labeled by the rank-p antisymmetric representation has total scaling dimension

    x_p = p(N-p)/N

and translation phase exp(2 pi i p/N). This is the standard WZW/translation identification, used here as a field-theory input, not derived from a proposed microscopic simulator.

Thus it would be too strong to say that every translation-invariant SU(N)-breaking correction necessarily turns on the most relevant WZW primary. One-cell translation can forbid those nonzero-momentum terms. Uniform currents and anisotropic current-current interactions still require analysis. In particular, translation does not determine the signs or flows of all marginal couplings.

The exact microscopic SU(N) protection strategy of Section 3 is unavailable for N > q. Other strategies remain possible: a suitable subgroup, an exactly constrained effective Hamiltonian, a carefully tuned critical manifold, or an emergent symmetry with a demonstrated basin of attraction. None is established by the second-order calculation.

### 6.4 What a Narain alternative would require

There is an important terminology issue in the supplied premise. A Narain theory always has dimension-one Abelian currents and dimension-two Gaussian deformations. Thus 'no primary of dimension <= 2 except the vacuum' cannot mean all Virasoro primaries. It can mean non-vacuum primaries of the extended U(1)^c current algebra, or equivalently the nonzero-charge vertex operators relevant here. The currents are in its vacuum module.

Accept the stated existence of suitable Narain charge lattices for every c >= 31. To use one, choose a point with a strict vertex-operator gap

    x_vertex >= 2 + kappa,    kappa > 0.                        (35)

The remaining low-dimension vacuum-module terms are currents and Gaussian quadratic terms. Linear currents shift backgrounds, and Abelian left-right current bilinears deform the Gaussian moduli. Small positive-definite quadratic deformations can preserve (35). This is the perturbative field-theory reason a stable open critical region is plausible.

A completed route still needs the following inputs, none proved here:

1. For each desired c, a finite-dimensional local lattice Hamiltonian whose infrared theory is such a Narain point, together with an open stability region for the microscopic perturbations under consideration.
2. Quantitative control of the distance to the boundary of that region, including the allowed current and velocity deformations. The existence of the continuum charge lattice alone does not provide a lattice Hamiltonian or this bound.
3. The fixed-dimension programmable-cell implementation of Section 5, with the local errors in (31) inside that stability region uniformly in L.

Section 2 removes an additional possible bottleneck: it can synthesize any finite-dimensional logical bond at second order, rather than only P_N. If a stable target requires a finite interaction range, finitely many target sites can first be grouped into a logical site so that the range becomes nearest neighbor. This grouping increases N and does not solve the microscopic implementation.

The ordinary SU(N) permutation chain is not by itself a construction of a high-minimum-dimension Narain lattice point. Switching the name of the target CFT without changing and analyzing its lattice Hamiltonian would not close the gap.

## 7. Spontaneous cells, gaps, and commensurability

### 7.1 What follows from the supplied selector

Accept the selector in the question, with K > 4t and mu chosen in the interval selecting period b. Its marker positions are conserved. In a fixed cell the walker is a path particle. From the supplied eigenvalue -2t cos(pi/b), the next path eigenvalue gives

    delta_walker = 2t [cos(pi/b) - cos(2pi/b)]
                 ~ 3 pi^2 t/b^2.                             (36)

For each fixed b this is positive, and t may be increased as a function of N to keep it above the inter-cell scale. Defect-selection energies can be smaller than the walker gap and must be compared separately with the coupling. The selector premise supplies a positive gap for each chosen b; it supplies no uniform-in-b lower bound, and none is assumed here.

Giving every letter a q-dimensional spectator costs 4q states per physical site. In a fixed legal marker sector, a cell can include its starting marker and the b-1 sites after it. Its left port can sit on that marker, and its right port on the site immediately before the next marker. A uniform bond term conditioned on an R# or X# pair then acts only across cell boundaries. Since every letter has its spectator, the port exists even when X reaches that boundary site. This implements the geometry of a boundary bond, provided the required code has already been realized inside the cell.

The 3q+1 variant needs additional care at a boundary that can be occupied by the spectator-free X. Its nominal dimension count cannot simply be substituted into this two-port construction without specifying how the port survives that occupation.

### 7.2 Adding computation must preserve selection

The original selector's spectrum cannot simply be carried over after arbitrary computational terms are added. If the same walker controls the computation, its hopping amplitudes and its cell ground state may change. If a separate program track is used, its ground energy may depend on the cell length and change the selected period.

A sufficient starting point would be nonnegative program terms with a common zero-energy state with the selector in the intended b-cell, together with proof that other sectors cannot lower their energy below the selector comparison. Once inter-cell couplings are added, one also needs a local defect-energy estimate. A crude extensive norm bound O(|lambda| L) cannot prove that one isolated defect remains unfavorable.

These are additional obligations on the missing local construction. The finite matrices in Sections 2 and 4 do not meet them by themselves.

### 7.3 Rings with L divisible by b

Conditionally on a successful construction, take L = b M and choose one of the b translated marker sectors. The effective logical chain has M sites. Broken translation gives finitely many marker sectors, not an additional propagating gapless mode.

If the logical chain has velocity v_cell measured in cells per unit time, its physical velocity is b v_cell. Thus

    -pi v_cell c/(6M) = -pi (b v_cell) c/(6L).

The central charge is unchanged by blocking. For a fundamental SU(N) chain, the finite-size representation sector can also depend on M modulo N; this does not change the proposed bulk central charge. Taking M divisible by N is an available further subsequence.

### 7.4 Rings with L not divisible by b

A uniform b-cell tiling is impossible. For fixed b, the remainder can be accommodated by a finite number of anomalous cell lengths, with a number bounded independently of L. Which pattern minimizes the energy requires the actual cell and defect energies. It need not be literally one anomalous cell in every remainder sector.

If the new interaction preserves the exact marker conservation, these defects are immobile. Their positions can produce a large static degeneracy. They can also act as impurity bonds, cut an effective chain, or bind local degrees of freedom. Their finite-size spectrum need not equal that of a clean periodic CFT.

At zero defect density, a stable critical bulk away from the defects would still have the same bulk central charge. This is a conditional bulk statement, not a claim that the incommensurate rings have the clean-ring Casimir term. A valid counterexample could state its clean critical limit along L in b Z while specifying these other sectors.

## 8. Dimension counts and comparison of the routes

No explicit physical d0 is established by this work. The explicit dimensions D = 32, 40, and 64 below are whole-cell matrix dimensions. The qubit port dimension q = 2 is not a physical-site dimension for a realization of those whole-cell matrices.

For the second-order route, the abstract result is particularly favorable: q can remain 2 while N is arbitrary, and (9) is polynomial in N. If a fixed-d0 programmable cell could realize it with controlled errors and a stable target, the port itself would place no bound on N or on the target central charge. For the SU(N)_1 target, an actual construction would violate c <= d0-1 when N > d0. An arbitrary-N construction with a valid infrared argument would give unbounded c at fixed d0. Neither implication has been established microscopically here.

For the direct first-order free ladder, the number of independent odd logical operators compressed from one port is at most q^2-1, since the identity must compress to the identity. Thus the route in Section 4 has a finite c bound at fixed q. It can potentially beat a linear c-d bound, but cannot by itself produce unbounded c by enlarging the cell at fixed q.

The following are only algebraic dimension comparisons. They do not include a verified programmable-cell implementation.

| Encoding | Port q | Channels r | First-order c | Nominal selector count | Comparison |
| --- | ---: | ---: | ---: | ---: | --- |
| Ungraded, even r | 7 | 48 | 24 | 3q+1 = 22 | 24 > 21, but the spectator-free boundary and parity still need repair |
| Ungraded, even r | 8 | 62 | 31 | 4q = 32 | Equality with d-1, no violation |
| Ungraded, even r | 9 | 80 | 40 | 4q = 36 | 40 > 35, before machinery and stability |
| Exact physical parity | 16 | 128 | 64 | 4q = 64 | 64 > 63, before machinery and stability |

The graded construction avoids the odd-channel degeneracy and makes fermion parity exact at the abstract-cell level. It sacrifices roughly half the ungraded channel count. It still allows interactions that need an infrared analysis. Any extra local program or clock alphabet changes the physical d0 and must be included in the inequality.

The practical comparison is therefore as follows. The second-order route has unbounded abstract logical capacity through fixed qubit ports and can synthesize arbitrary target bonds, but it needs detailed excited-state simulation and cannot use the straightforward exact SU(N) protection. The first-order route needs only ground-code compression to obtain its leading bond and has an explicit joint encoding, but its channel count is bounded by the port algebra and its higher-order interactions are already visible in the smallest checks.

## 9. Numerical checks

### 9.1 Scope and reproducibility

The code uses only NumPy and SciPy. All diagonalizations are finite-dimensional calculations on the abstract cells defined here. They are not calculations on a translation-invariant physical nearest-neighbor realization inside a cell.

Run from this directory with

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 checks.py
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 spectral_cell_check.py

The first script writes checks_results.json. The second writes spectral_cell_results.json. The saved stdout files also contain the reports. The shell did not have a `python` executable; the initial invocation failed with `command not found`, and all completed computations used `python3`.

The large two-cell Hilbert spaces are treated by sparse low-eigenvalue diagonalization with exact matrix-vector products, not by truncating the physical cell Hilbert space. Residual norms of the returned eigenpairs were also checked. The largest reported residual is about 1.1e-12. The tested perturbative errors are larger than this numerical floor.

No numerical calculation establishes a thermodynamic central charge or a stability basin.

### 9.2 Joint CP cells

For q = 2, take s = 0.1 and Delta = 1 in (18)-(21).

| Quantity | r = 2 | r = 3 |
| --- | ---: | ---: |
| Logical dimension N | 2 | 4 |
| Abstract cell dimension D | 32 | 64 |
| Two-cell Hilbert dimension | 1024 | 4096 |
| Smallest Choi eigenvalue | 0.1792893219 | 0.1439339828 |
| Isometry error | 2.3e-16 | 7.5e-16 |
| Largest single-port compression error | 9.9e-17 | 7.9e-16 |
| Largest cross-port product compression | 3.2e-17 | 1.2e-15 |
| Exact abstract cell gap | 1 | 1 |
| Norm of [H_2, P tensor P] | 2.0e-17 | 0.012 |

The gap is exact analytically because H_A = 1-GG†. The numerical isometry test checks that the constructed G realizes that projector to roundoff. The r=3 parity commutator agrees with the nonzero term (26).

For each cell, the first-order matrix agrees with s^2 T. I compared the lowest N^2 exact two-cell eigenvalues against the eigenvalues of lambda H_1 + lambda^2 H_2. The maximum absolute spectral errors are

| lambda | r = 2 error | r = 3 error |
| --- | ---: | ---: |
| 0.04 | 1.05e-7 | 1.015e-4 |
| 0.02 | 1.20e-8 | 1.265e-5 |
| 0.01 | 1.43e-9 | 1.578e-6 |
| 0.005 | 1.75e-10 | 1.971e-7 |

For r=3, error/lambda^3 approaches 1.5765. For r=2 it approaches about 0.0014 over the tested range. The explicit r=2 expression (25) agrees with the direct second-order matrix to 4.5e-16. Its non-scalar part has norm 0.0003, so the non-free correction is resolved numerically.

On four open logical cells, direct construction of the Jordan-Wigner operators agrees with (15) exactly in floating-point arithmetic. The r=2 ladder has a nondegenerate ground state with energy -4.47213595499958 in unit hopping. The r=3 ladder has a fourfold ground state with energy -6.70820393249938. Every unused global Majorana commutes with its Hamiltonian, with computed commutator norm zero. This confirms the unused-mode mechanism in Section 4.3.

### 9.3 Two-port spectral gadget

The script spectral_cell_check.py implements Section 2 for

    q = 2,    N = 2,    m = 1,    k = 3,    D = 40,
    (e_1,e_2,e_3) = (1,4,16),    p_e = 1/3,
    A_1 = X,    B_1 = Z,    K_11 = 1.

The unused complement is assigned energy 32. The construction gives

    rho = 0.02522559045216786,
    rho^2 = 0.0006363304136605024,
    -r p^T C p = -0.42463235294117646.

The exact second-order prediction in units Delta=1 is

    H_2 = -0.42463235294117646 1
          + 0.0006363304136605024 X tensor Z.                  (37)

Checks gave

    G†G - 1                         norm 0,
    S†S - 1                         norm 0,
    sum_e M_e - 1                   norm 1.4e-17,
    max_e ||S†P_eS - M_e||                3.7e-16,
    ||H_A G||                            1.8e-15,
    first excited energy                 0.999999999999987,
    first-order norm                     0,
    error in (37)                        1.2e-15.

The smallest weight eigenvalue is 1/6, so this is not a numerically marginal positivity construction.

The two-cell Hilbert dimension is 1600. Exact low-spectrum comparisons give

| lambda | Maximum second-order spectral error | Error/lambda^3 | Exact logical bandwidth | Predicted bandwidth |
| --- | ---: | ---: | ---: | ---: |
| 0.04 | 7.09e-6 | 0.11080 | 2.1771e-6 | 2.0363e-6 |
| 0.02 | 8.85e-7 | 0.11068 | 5.2674e-7 | 5.0906e-7 |
| 0.01 | 1.11e-7 | 0.11062 | 1.2948e-7 | 1.2727e-7 |
| 0.005 | 1.38e-8 | 0.11059 | 3.2093e-8 | 3.1817e-8 |

The next eigenvalue remains above 0.983 throughout this range. The logical bandwidth tends to 2 rho^2 lambda^2, and the absolute spectral residual is cubic. The small splitting of the predicted degenerate pairs at nonzero lambda is a higher-order effect, not evidence of exact all-order target dynamics.

Both scripts ran in seconds. Repeated runs used only a small fraction of the requested runtime budget.

## 10. Requirement-by-requirement status

| Requirement | Established here | Remaining gap |
| --- | --- | --- |
| (a) Programmable local cell | Exact abstract two-port spectral gadget for every N, including any target bond at second order; explicit finite-clock accuracy requirements | Fixed-d0 translation-invariant nearest-neighbor rules selecting and running the correct program, excluding every unwanted sector, and preserving port observables |
| (b) Degeneracy, first order, machinery gaps | Exact code degeneracy, zero first order, and finite gap for the abstract second-order cells; joint scaled Clifford encodings for the first-order route | Gaps and degeneracy for actual local machinery; compatibility with the selector; uniform local effective-Hamiltonian control |
| (c) Error stability | Exact on-cell SU(N) protection is obstructed for N>q; explicit higher-order Majorana corrections; quantitative scale hierarchy | A protecting mechanism outside the no-go assumptions, a tuned critical manifold, or a stable lattice Narain target with controlled error basin |
| (d) Commensurability | Clean conditional statement on L in b Z; static finite-density-zero defects otherwise | Actual impurity spectra and selection bounds after the missing local implementation is added |
| (e) Explicit d0 and c-d consequence | Explicit whole-cell dimensions and port dimensions, with conditional central-charge counts | No explicit microscopic d0 and no counterexample |

The finite-dimensional port problem is more flexible than the original swap gadget alone suggests. The obstacle is now concentrated in two places: (i) a genuinely local programmable cell with the required boundary-resolvent guarantees, and (ii) an infrared argument that survives every finite implementation error. Neither can be replaced by register capacity or finite-size agreement.
