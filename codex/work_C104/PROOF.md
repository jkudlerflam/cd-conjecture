# C104. Conserved charges and the velocity-weighted bound

## 1. Status and hypotheses

[Assumed] The chain, gauge-optimized half-spread M, infrared sectors and A = sum_a c_a v_a have the meaning in ../../specs/_bound_background.md, including its updates. Stiffness limits exist along vacuum ring sequences on which the usual twist response is defined. The lattice-to-CFT current-response dictionary is an explicit input. A finite-size CFT spectrum alone does not prove convergence of the second derivative with respect to flux.

[Proved, conditional on that dictionary] This note proves a crossover-free bound for full simple WZW factors seen by a conserved on-site semisimple symmetry,

    A_WZW <= (4 pi/3) d(d^2 - 1) M.                         (1)

An upper bound on affine levels is not needed. The sharper version is (21) below. Multiple WZW factors with different levels and velocities can couple diagonally to the same microscopic simple Lie algebra.

[Proved, limitation of the argument] Quantization and locality of compact-boson vertex operators do not yield a lower bound on their current weights. Sections 3 and 4 demonstrate this within the continuum theory, including the proposed scaling-dimension argument. This is not a fixed-d nearest-neighbor lattice counterexample to a lower bound on K. Such a lattice lower bound, or a different sufficient aggregate inequality, remains open here. Consequently the requested theorem for every charged Gaussian infrared is not proved.

[Assumed] When left and right velocities differ, the stiffness dictionary below uses their arithmetic mean, consistent with A. The nonchiral sectors have matched left and right current levels and the on-site symmetry is vector-like.

## 2. The lattice budget and normalizations

[Proved, imported lattice lemma] For commuting traceless one-site generators q^alpha with

    tr_d(q^alpha q^beta) = delta_alpha,beta,

and conserved total charges, the pairing note proves

    D_alpha <= -<[q^alpha,[q^alpha,h_tilde]]>,
    sum_alpha D_alpha <= 4 M.                              (2)

The generators on the right act on the first site. Use an invariant centered representative with ||h_tilde|| = M, or a sequence approaching M. Group averaging makes this compatible with the symmetry. The sign in the expectation in (2) is minus. The displayed torus formula in the Claim of the pairing README has the opposite sign, whereas its derivation has the correct sign.

[Proved] If Q is an integer-valued one-site charge, set

    Q_c = Q - tr(Q)/d,
    S_Q = tr(Q_c^2),
    q = Q_c/sqrt(S_Q).

For a nontrivial Q, twisting by q is twisting by Q with flux divided by sqrt(S_Q). Hence

    D_q = D_Q/S_Q.                                         (3)

The normalized q generally does not have integer eigenvalues. Integer charge labels and Hilbert-Schmidt normalization must not be conflated.

[Proved] Having at most d integer eigenvalues does not imply a window of width d. For example, diag(0,1,N) is primitive and has three values for every N > 1. We also consider the stronger, explicitly requested window assumption W = max Q - min Q <= d. It implies

    S_Q <= d W^2/4.                                        (4)

To verify (4), put a <= lambda_i <= b. Averaging (lambda_i-a)(b-lambda_i) >= 0 gives Var(lambda) <= (mean-a)(b-mean) <= (b-a)^2/4, then multiply by d. This controls normalization, not an infrared stiffness from below.

[Assumed, Gaussian response dictionary] For canonically normalized Gaussian modes a and the orthonormal microscopic torus, let Z_a be the row of charge couplings. Then

    D_alpha,beta = (1/pi) sum_a v_a Z_a,alpha Z_a,beta,
    pi tr D = sum_a v_a w_a,       w_a = |Z_a|^2.           (5)

For a Gaussian infrared A_G = sum_a v_a. Therefore (2) gives only

    sum_a v_a w_a <= 4 pi M.                               (6)

[Proved, conditional statement] If w_a >= w_*(d) > 0 for every contributing mode, then A_G <= 4 pi M/w_*(d). Merely assuming w_a > 0 for each mode does not supply a uniform w_*(d).

## 3. What charge quantization says about a single boson

[Assumed, continuum normalization] Take fields theta, phi with theta periodic by 2 pi, phi periodic by pi, and Hamiltonian

    H = (v/(2 pi)) integral dx [ K (partial_x theta)^2
                               + K^(-1) (partial_x phi)^2 ],
    rho = -partial_x phi/pi,
    [phi(x),partial_y theta(y)] = i pi delta(x-y).          (7)

The charge Q shifts theta. With a consistent choice of the sign of Q, exp(i m theta) has charge m. Only the sign convention changes if rho is reversed. The local bosonic vertex lattice can be chosen as

    V_m,n = exp(i m theta + 2 i n phi),    m,n integers.

[Proved within (7)] Integrating out phi gives the theta action K/(2 pi) integral (partial theta)^2 after rescaling Euclidean time by v. The logarithmic propagator is

    <theta(z)theta(0)> = -(1/(2K)) log |z| + constant.

Thus exp(i m theta) has dimension m^2/(4K). The dual action gives dimension n^2 K for exp(2 i n phi). More explicitly,

    h_m,n = (m/(2 sqrt(K)) + n sqrt(K))^2/2,
    hbar_m,n = (m/(2 sqrt(K)) - n sqrt(K))^2/2,
    x_m,n = m^2/(4K) + n^2 K,
    h_m,n - hbar_m,n = m n.                                (8)

The mutual monodromy exponent of two vertices is m n' + n m', an integer independent of K. Appropriate cocycles give the bosonic local vertex algebra. Every K > 0 defines a unitary local compact-boson theory with the same integer labels. The radius deformation changes dimensions continuously. Standard compact-boson quantization and its radius modulus are reviewed in [Ginsparg, Sections 8.1 and 8.6](https://arxiv.org/pdf/hep-th/9108028).

[Proved within (7)] A flux Phi makes theta(L)-theta(0) = Phi on the vacuum branch. Minimizing (7) gives

    E_0(Phi)-E_0(0) = v K Phi^2/(2 pi L),
    D_Q = L partial_Phi^2 E_0 = v K/pi.                   (9)

Changing the total charge by Delta Q instead gives

    E(Delta Q)-E(0) = pi v (Delta Q)^2/(2 K L),
    kappa_Q = K/(pi v).                                   (10)

This independently checks the factors of pi and K in (8) and (9).

[Proved] The minimum charge-one dimension in this vertex lattice is

    x_1,min = 1/(4K).                                     (11)

There is no upper bound on (11) from charge integrality or mutual locality. It diverges when K tends to zero. Even at a fixed K, the dimensions of all charge-one vertices are unbounded as |n| grows. The existence of a charge-one operator is not the existence of a charge-one operator of uniformly bounded dimension.

[Proved] A bounded lattice operator does not close this gap. A possible infrared expansion O_lat = b V_1,0 + higher fields gives a correlator proportional to |b|^2 |x|^(-1/(2K)) beyond the matching scale. The bound ||O_lat|| <= 1 puts no upper bound on the decay exponent. Neither b nor the matching scale is uniformly fixed by the one-site charge range. A charge-one lattice operator can also couple only to massive states, so its existence does not even guarantee a charge-one gapless primary.

[Proved, sufficient extra hypothesis for one sector] If a charge-one gapless primary of dimension at most X(d) is supplied independently, (8) implies K >= 1/(4X(d)). Combining (3), (4), and (9),

    w = K/S_Q >= 1/[d W^2 X(d)],
    v <= 4 pi d W^2 X(d) M.                              (12)

The coefficient in (12) uses the general torus budget and is not optimized. The required upper bound on X(d) is new input, not a consequence of charge quantization. For several modes, one would need a corresponding mode-resolved statement or a frame bound on their charge couplings.

[Proved, stability does not replace the missing hypothesis] If translation and symmetry allow an umklapp cos(2p phi), its dimension is p^2 K. Stability against a nonzero coefficient requires p^2 K >= 2, away from the relevant side. This is a constraint on a specified stable phase with a specified allowed p. It is not a theorem for all critical Hamiltonians: an allowed coefficient can be tuned to zero, and incommensurate translation quantum numbers can forbid the operator. The charge window does not bound the commensurability denominator.

[Proved, scope] The K -> 0 continuum family establishes failure of the proposed continuum inference. It does not establish that this family has a one-site invariant nearest-neighbor realization at one fixed finite d. No such realization is asserted in this note.

## 4. Several charged modes and the actual missing estimate

[Proved within Gaussian field theory] Integer charge labels also do not give a mode-by-mode coupling bound. For example, take two compact pairs (phi_i,theta_i), each with the periods of Section 3, and

    H = (1/(2 pi)) integral dx [ (partial theta)^T B (partial theta)
                                + (partial phi)^T B (partial phi) ],
    B = O diag(v_1,v_2) O^T,
    Q = -(1/pi) integral partial_x phi_1.

Here O is a rotation by epsilon and v_1 != v_2 are positive. The normal modes have K_a = 1 and couplings to Q equal, up to signs, to cos(epsilon) and sin(epsilon). Both are nonzero for epsilon != 0, but the second squared weight tends to zero. All charges of local vertices in the original basis remain integer. The transformed compactification lattice is generally oblique, so this example concerns coupled toroidal modes, not a claimed product of two independently compactified circles in the normal-mode basis.

[Proved, qualification] This example has a second conserved continuum U(1). If both currents were supplied as suitably normalized microscopic charges, the sum of their weights for each mode would equal one. It demonstrates the insufficiency of coupling nontrivially to a specified charge, and is not a counterexample for the full microscopic conserved torus. The single-boson calculation in Section 3 already shows that integral charge labels do not constrain K.

[Proved] A rank count cannot fix the general problem. One charge can have nonzero projection onto many modes. The number of nonzero rows of a rectangular Z is not bounded by its column rank. Replacing “every mode has some charge coupling” by “there is an independent microscopic conserved charge for every mode” adds a new hypothesis.

[Proved, sufficient aggregate replacement] The precise estimate needed by this route can be weaker than a pointwise lower bound. It suffices to prove

    sum_a v_a <= G(d) sum_a v_a w_a                       (13)

for the microscopic realizations under consideration. Then (6) gives A_G <= 4 pi G(d) M. No such estimate is established here. Small w_a could be harmless if the corresponding velocities are sufficiently small. An unweighted count of sectors discards this possibility.

[Proved, normalization audit] The matrix compressibility formula stated in Section 2 of the pairing README, Z^(-1) V^(-1) Z^(-T)/pi, fails its own one-boson check. In the convention of (5), where the normalized mode charge densities have susceptibility 1/(pi v_a), the susceptibility of the microscopic charges is

    kappa = Z^T V^(-1) Z/pi.                              (14)

For one charge Z = sqrt(K), this is K/(pi v), as in (10). None of the results here use an inverse-Z compressibility formula.

## 5. WZW response with the on-site trace kept explicit

[Assumed] A compact semisimple algebra g = direct_sum_j g_j acts on C^d through R, and the bond is invariant under its diagonal action. Discard ideals acting trivially. The counted infrared sectors are full unitary g_j,k_a WZW factors, with positive integer k_a and velocities v_a. The microscopic g_j current couples to the diagonal sum of the standard currents of these factors. Its action on a counted factor may differ by an automorphism, which preserves the invariant form. A larger current algebra seen only through a proper subalgebra requires a separate coset accounting and is not silently included.

[Definition] For g_j choose a canonical invariant metric and Hermitian generators T_j^A such that a long-root su(2) subalgebra has its usual spin generators, [S^x,S^y] = i S^z, with S^z having eigenvalues +/-1/2 in its doublet. Each basis generator has the same canonical norm as this S^z. Define I_j by

    tr_R(T_j^A T_j^B) = I_j delta_AB.                     (15)

For the fundamental of SU(N), I_j = 1/2. For an SU(2) spin-S site, I_j = S(S+1)(2S+1)/3. This definition fixes any convention-dependent factor in the name “Dynkin index.” The trace is the ordinary trace, not tr/d.

[Assumed, WZW stiffness dictionary] In this canonical basis,

    D_T_j^A = (1/(2 pi)) sum_{a in j} k_a v_a.             (16)

The chiral current two-point coefficient is k/2 in the same basis; the left and right vacuum spectral flows under Phi each contribute (k/2)(Phi/(2 pi))^2/2 to their conformal energy. Multiplying their sum by 2 pi v/L gives k v Phi^2/(4 pi L), hence (16). This checks the continuum flux coefficient. Identification with the lattice stiffness remains the explicit matching assumption.

[Proved from (15) and (16)] The Hilbert-Schmidt normalized generators t_j^A = T_j^A/sqrt(I_j) obey

    D_t_j^A = (1/(2 pi I_j)) sum_{a in j} k_a v_a.         (17)

The formula D_t = k v/pi for tr_R(t^2)=1 is valid when I_R=1/2, not in a general on-site representation. This is the representation-normalization correction to the pairing README.

[Assumed, standard WZW input] The Sugawara central charge is

    c(g_j,k) = k m_j/(k+h_j),
    m_j = dim g_j,       h_j = h_j^vee.

The current algebra, Sugawara construction, and Lie-algebra constants can be found in [Eberhardt, Sections 3.2 to 3.4](https://www.conferences.itp.phys.ethz.ch/esi-school/Lecture_notes/WZW%20models.pdf). We use h_j >= 2 and k >= 1.

[Proved, one simple factor via its Casimir superoperator] Let C_t X = sum_A [t^A,[t^A,X]]. Adding the lattice variational inequalities gives

    sum_A D_t^A <= ||C_t h_tilde||.

Consequently, even for several levels and velocities coupled to this same factor,

    sum_a c(g,k_a)v_a
      <= [m/(1+h)] sum_a k_a v_a
      = [2 pi I/(1+h)] sum_A D_t^A
      <= [2 pi I/(1+h)] ||C_t h_tilde||.                  (18)

For a single level k, replace 1+h by k+h. Since C_t = C_T/I, the last expression then equals 2 pi ||C_T h_tilde||/(k+h). The index cancels only if the superoperator uses the canonical, rather than the on-site orthonormal, generators consistently. The operator-norm estimate need not assume that h_tilde lies in a single irreducible component of End(R).

## 6. A uniform polynomial WZW bound

[Proved] Use only Cartan generators in (17). Cartans for different simple factors commute and are trace orthogonal. Trace orthogonality follows by decomposing R into sums of tensor products of irreducibles: every generator of a simple Lie algebra has zero trace in each irreducible, so a cross-factor trace vanishes. Thus all normalized Cartans together are an admissible torus for (2).

Let r_j = rank g_j, and define the nonnegative Cartan stiffness trace

    S_j = sum_{i=1}^{r_j} D_t_j^i
        = r_j/(2 pi I_j) sum_{a in j} k_a v_a.

Then

    A_j = sum_{a in j} [k_a m_j/(k_a+h_j)] v_a
        <= m_j/(1+h_j) sum_{a in j} k_a v_a
        = B_j S_j,
    B_j = 2 pi I_j m_j/[r_j(1+h_j)].                      (19)

All S_j are nonnegative by the assumed unitary WZW dictionary. Equation (2) therefore gives

    A_WZW <= 4 M max_j B_j.                              (20)

Equivalently the representation-sensitive result is

    A_WZW <= 8 pi M max_j [I_j dim(g_j)/
                           (rank(g_j)(1+h_j^vee))].      (21)

This is already a crossover-free bound for any specified on-site representation, independent of the number of infrared copies and of their levels.

[Proved, index bound] Restrict R to a long-root SU(2) of g_j. Write this restriction as a sum of spin-s multiplets, dimensions n_l = 2s_l+1, with sum_l n_l = d. The chosen S^z is one of the canonically normalized generators, so

    I_j = tr_R((S^z)^2)
        = sum_l s_l(s_l+1)(2s_l+1)/3
        = (1/12) sum_l (n_l^3-n_l)
        <= d(d^2-1)/12.                                  (22)

The last inequality follows from sum n_l^3 <= (sum n_l)^3, including multiplicities and singlets. This proof does not identify the on-site representation with an affine primary.

[Proved, finite Lie-algebra inequality] For every compact simple algebra,

    dim g/[rank(g)(1+h^vee)] <= 2.                        (23)

For simply laced A_r,D_r,E_6,E_7,E_8 the ratio is one. For B_r it is (2r+1)/(2r), for C_r it is (2r+1)/(r+2), for F_4 it is 13/10, and for G_2 it is 7/5. These expressions follow directly from the root-system dimensions and dual Coxeter numbers in the cited table. The classical ranges are B_r with r >= 2, C_r with r >= 1, D_r with r >= 4, with the usual low-rank identifications handled by A. In particular (23) is uniform in the algebra.

[Proved, conditional theorem] Combining (20), (22), and (23) yields

    max_j B_j <= 4 pi max_j I_j <= pi d(d^2-1)/3,
    A_WZW <= (4 pi/3) d(d^2-1) M,                        (24)

which proves (1). For only simply laced factors, the same argument gives the improved coefficient 2 pi/3 in place of 4 pi/3. For d=1 there is no nontrivial symmetry factor and A_WZW=0.

[Proved, why no level upper bound occurs] The important ratio is c/k = dim g/(k+h^vee), which decreases as k increases. The inequality k >= 1, not an upper bound on k, produces (19). Slow sectors at large k are therefore compatible with the result. Equation (24) bounds a velocity-weighted sum and does not bound the total central charge or the number of copies.

[Proved, limitation of a representation-to-level argument] A microscopic representation fixes the zero-mode Lie algebra and I_R. The affine level multiplies the derivative term in the current commutator. The zero modes do not see this term. Integrability of a specified affine primary of spin s imposes 2s <= k in SU(2), which is in the opposite direction from a level upper bound, and an on-site state need not be that primary. Constructions using free-fermion conformal embeddings add a specific chiral realization whose index fixes a level. Such a realization is not part of the present hypotheses.

[Assumed, literature input] Under the symmetry and translation assumptions in Furuya and Oshikawa's analysis, the SU(2) anomaly preserves the parity class of k under RG. This is a congruence condition, not an upper bound on k. See [their paper](https://arxiv.org/abs/1503.07292). We neither assert nor need a general theorem that fixed-d nearest-neighbor chains realize arbitrary levels.

## 7. Neutral stress tensors and mixed infrared theories

[Proved within the infrared decomposition] A U(1)-neutral sector contributes zero to D for that symmetry, but contributes c v to A. Being neutral means that the symmetry has no current component in that factor, not that every local observable in it has zero scaling dimension.

[Proved example, using the solvable XX and Ising infrared dictionaries] A d=4 site can be C^2_A tensor C^2_B. Take

    H_epsilon = epsilon H_XX,A + H_Ising,B,
    Q_x = n_A,x tensor 1_B,         eigenvalues 0,1.

Both summands are translation-invariant and nearest-neighbor, with the transverse Ising field split equally between adjacent bonds. At their critical points,

    D_Q = epsilon v_XX/pi,
    A = epsilon v_XX + v_Ising/2.

As epsilon tends to zero through positive values, A/D_Q diverges. All on-site charges remain in a window of width one and d stays four. This refutes control of the entire A by this charge stiffness alone when an independent neutral factor is allowed. It does not violate any proposed A/M bound.

[Proved, distinction for parafermions] The standard decomposition of SU(2)_k has a Cartan boson and a Z_k parafermion coset,

    c_SU(2)_k = 1 + 2(k-1)/(k+2) = 3k/(k+2).

When the full SU(2)_k current algebra is realized by the on-site symmetry and these pieces share its WZW velocity, (24) bounds their combined contribution. At k=2 this includes the Ising coset. It would be incorrect to discard that contribution merely because it is neutral under the Cartan U(1). The non-Cartan currents and the Sugawara stress tensor supply the relation.

[Proved, scope of the remaining gap] An independent Ising factor, a standalone parafermion theory with only discrete symmetry, or a coset stress tensor not fixed by the counted symmetry currents is not controlled by (24). If only the Cartan U(1) survives microscopically, the WZW relation between its stiffness and the parafermion velocity cannot be assumed. Different neutral and charged velocities require separate accounting.

[Proved, conditional combined estimate] Suppose the commuting set consists of the normalized central U(1) generators and all the Cartans above, and the Gaussian part couples only to the central U(1)s while the WZW parts have the stated semisimple couplings. If every Gaussian row weight obeys w_a >= w_* > 0, (5) and (19) give the shared-budget inequality

    A_G + A_WZW <= 4 M max{ pi/w_*, max_j B_j }.           (25)

This is sharper than adding two independently exhausted budgets. Any residual neutral contribution A_neutral must be added to the right and requires a new bound. Neither w_*(d) nor such a neutral bound has been obtained here.

## 8. Consistency checks and conclusion

[Proved] Scaling h by a positive lambda multiplies M, all velocities, all stiffnesses and both sides of (24) by lambda. No temperature, length, or inverse velocity enters the constant. In lattice units the dimensions therefore agree.

[Proved, dictionary check] For the spin-1/2 Heisenberg chain, I=1/2, k=1, v=pi/2. Equations (16) and (17) give D_Sz=1/4 and D_sqrt(2)Sz=1/2. Since c=1, A=pi/2=pi D_sqrt(2)Sz, agreeing with the exact sector conversion (19). The norm estimate on stiffness is looser than its actual value.

[Proved, dictionary check] For a spin-1 SU(2)_2 sector, I=2. The same equations give D_Sz=v/pi and D_Sz/sqrt(2)=v/(2 pi). The latter is a factor four below k v/pi. This exposes the fundamental-versus-on-site trace issue directly. For a fundamental SU(N) site, I=1/2 and the torus trace is (N-1)k v/pi, as in the pairing note.

[Numerical] check.py performs only small deterministic matrix and arithmetic checks of (3), the spin index in (22), the Casimir and torus estimates, the Lie-algebra ratios in (23), and the Heisenberg normalization. These checks do not establish a lattice-to-CFT limit or a lower bound on Gaussian K. See CHECKS.md for their output.

[Proved, final status] The WZW subclass has the explicit polynomial bound (24), conditional on its stated infrared response identification. The general charged Gaussian theorem is still missing (13) or a stronger lower bound on weights. Charge quantization and vertex locality alone do not provide it. Independent neutral stress tensors require additional input even after that Gaussian gap is filled.
