# C121. Coarse seam measurements and the missing reconnection estimate

2026-10-29. Natural logarithms. Claim tags specify the status of each result. All measurements below are nonselective unless a history is explicitly fixed.

## 1. Setting and what is obtained

[Assumed] Use the translation-invariant nearest-neighbor Hamiltonians H_L and optimized bond scale M of the background specification. Write Z_L = Tr exp(-beta H_L), R_L = 2 log Z_L - log Z_(2L), and Delta_L = E_(2L)-2E_L. The physical on-site dimension is d. Vacuum-compatible lengths are required when extracting A. Scalar bond shifts and coboundaries leave every ratio in this note invariant.

[Proved, inherited] C116 and C118 establish that an all-temperature estimate R_L <= K(d)+k M f(d) beta/L implies A <= (4k/pi) M f(d) under the vacuum Casimir expansion. They also establish the general bound |R_L| <= 4 beta M and the stronger bounds for classical and XX chains. These are inputs, not rederived here.

[Proved here] A seam pinching can be inserted at arbitrary times at a cost at most log d per insertion in the log partition trace, in the positive channel construction below. No gap, Markov property, stoquasticity or Lieb-Robinson bound is needed for this assertion. The recorded histories obey Cauchy-Schwarz. Their squared weights do not obey the proposed unconditional comparison with the reconnected ring. Coarse pinching leaves quantum propagation through the seam between insertions.

[Open] A uniform reconnection estimate with a controlled exponential cost per coarse slice is stated in Section 8. Neither this lemma nor a general bound on A is proved. The failure of the zero-cost reconnection claim is not a counterexample to the target bound.

## 2. Coherent Euclidean histories and exact Cauchy-Schwarz

[Proved definitions and identities] Let P_s = |s><s| tensor I_rest be the d seam projectors, let N >= 1, tau = beta/N, and define positive operators

    B_s = exp(-tau H_L/2) P_s exp(-tau H_L/2),
    a_sbold = Tr(B_(s_N) ... B_(s_1)),
    S_coh = sum_sbold |a_sbold|^2.                              (1)

The history sbold ranges over d^N labels. Cyclically moving the half propagators gives the equivalent expression Tr[P_(s_N) exp(-tau H_L) ... P_(s_1) exp(-tau H_L)]. Summing independently over every label replaces each projector sum by I, so

    sum_sbold a_sbold = Z_L.                                   (2)

In particular, merely inserting resolutions of identity has no error and does not implement a measurement channel. Products of positive operators in (1) need not have positive or real traces for a general complex Hamiltonian and N >= 3.

[Proved operator Cauchy-Schwarz] Introduce the diagonal operator C = diag(a_sbold) on the classical history register of dimension D = d^N. Hilbert-Schmidt Cauchy-Schwarz gives

    |Tr C|^2 <= Tr I_D * Tr(C^dagger C),
    Z_L^2 <= d^N S_coh.                                        (3)

This is the exact quantum counterpart of the first classical step after the quantum output has been traced. It remains true with complex histories. Applying Cauchy-Schwarz before tracing the unmeasured quantum output does not remove that output's dimension.

[Open candidate, false without further restrictions] The classical second step would require

    S_coh <= Z_(2L).                                          (4)

There is no identity identifying S_coh with a positive sub-sum of the doubled-ring trace. For N = 1, onsite symmetry in XX, Heisenberg and SU(3) gives a_s = Z_L/d, and hence

    log[S_coh/Z_(2L)] = R_L-log d.                             (5)

C118 already shows that R_L grows without bound at fixed L as beta increases whenever Delta_L > 0. Thus (4) cannot be an unconditional positivity theorem. Section 6 gives a finite, analytic XX witness with tau a fixed multiple of L/M, as well as tests at tau = L/M.

## 3. A positive pinching channel with an exact error bound

[Proved construction] Define

    D(X) = sum_s P_s X P_s,
    A_t(X) = exp(-tH/2) X exp(-tH/2).

Here H can be any finite-dimensional Hermitian matrix. Apply the imaginary-time maps for total time beta, inserting N pinchings at the midpoints of N equal intervals of length tau = beta/N. Explicitly, put

    F = exp(-tau H/4),       E = exp(-tau H/2),
    W_sbold = F P_(s_N) E P_(s_(N-1)) E ... E P_(s_1) F,
    q_sbold = Tr(W_sbold W_sbold^dagger) >= 0,
    Q_(H,N) = sum_sbold q_sbold.                               (6)

For N = 1, W_s = F P_s F. In channel notation, start from I, apply A_(tau/2), then D, then N-1 copies of A_tau followed by D, and finally A_(tau/2), and take the trace. This gives (6).

[Proved distinction] The trace in (6) contains two reflection-paired copies of each projector, one on each side of the evolving operator. It is a positive measurement-channel partition functional, not the coherent sum (2). Its N record labels correspond to 2N projector appearances on the unfolded Euclidean contour. This fixed factor does not change the desired beta M/L scaling, but the distinction matters for spatial gluing. This functional is also not a superoperator trace, nor is it Tr[(D(exp(-tau H)))^N].

### 3.1 Lower bound

[Proved] For every X >= 0, pinching satisfies

    D(X) >= X/d.                                               (7)

For any vector v, Cauchy-Schwarz on the d vectors X^(1/2)P_s v gives

    <v,Xv> = ||sum_s X^(1/2)P_s v||^2
           <= d sum_s ||X^(1/2)P_s v||^2 = d <v,D(X)v>.

Each A_t and D preserves operator order. Replace the first D by d^(-1) times the identity map using (7), propagate this order through the remaining maps, and repeat for every insertion. The maps without measurements compose to A_beta. Taking traces yields

    Q_(H,N) >= d^(-N) Tr exp(-beta H).                          (8)

This is a lower bound per insertion even if other measurements are already present. No statement that adding an insertion always decreases Q is needed.

### 3.2 Upper bound

[Proved] Let omega = exp(2 pi i/d), and U = sum_s omega^s P_s. Summing the geometric series gives

    D(X) = (1/d) sum_(r=0)^(d-1) U^r X U^(-r).                 (9)

Expand every D using (9). Each resulting term in Q is ||V||_2^2, where V is an ordered product of factors exp(-t_j H/2), interspersed with unitaries, and t_j >= 0 with sum_j t_j = beta. For t_j > 0 choose p_j = 2 beta/t_j. Then sum_j 1/p_j = 1/2, and Schatten Holder gives

    ||V||_2 <= product_j ||exp(-t_j H/2)||_(p_j)
             = product_j Z_H^(t_j/(2 beta)) = Z_H^(1/2).        (10)

The unitary factors have operator norm one. Zero-time factors can be omitted. Averaging (10) over the d^N terms of (9) gives Q <= Z_H. At beta = 0 the equality follows directly. Combining with (8),

    d^(-N) Z_H <= Q_(H,N) <= Z_H,
    0 <= epsilon_(H,N) := log Z_H-log Q_(H,N) <= N log d.        (11)

The same proof works at arbitrary unequal insertion times and with different projective bases at different times. If insertion j has r_j projectors, replace d^N by product_j r_j. Local rank-one projectors are not rank one on the full ring; their number, rather than their ranks, enters (7).

[Proved sharpness] Suppose H has a unique ground vector psi with <psi|P_s|psi> = 1/d. At fixed N, let tau tend to infinity. Factoring the ground energy out of each propagator gives

    W_sbold = exp(-beta E_0/2) d^(-N) |psi><psi| + o(exp(-beta E_0/2)),
    Q_(H,N)/Z_H -> d^(-N).                                    (12)

Thus the factor d per insertion cannot generally be improved. This is a fixed finite-system, fixed-N limit. It does not assert uniform ground dominance at tau proportional to L.

### 3.3 The proposed time scale

[Proved] Fix a dimensionless a > 0. For M > 0 and beta > 0 choose

    N = max(1, ceil(beta M/(a L))),       tau = beta/N.          (13)

Then N <= 1+beta M/(a L). If beta M/(a L) >= 1, one also has a L/(2M) <= tau <= a L/M. Therefore (11) gives the desired error budget

    epsilon_(L,N) <= log d + (log d/a) beta M/L.                (14)

This conclusion holds for every a, without choosing a seam level spacing. For the reconnected ring use the same beta, N and tau, with both seam sites 0 and L pinched. Each insertion then has d^2 projectors, giving

    0 <= epsilon_(2L,N) <= 2N log d.                           (15)

For M = 0, the ring Hamiltonians are scalar up to a telescoping gauge term and R_L = 0. The beta = 0 case also has R_L = 0.

## 4. Positive-history Cauchy-Schwarz and the missing second step

[Proved] Write Q_L for (6) with one seam on the L-ring, and Q_D for the same channel with both seam sites pinched on the reconnected 2L-ring. Define

    S_pin = sum_sbold q_(L,sbold)^2,
    G_(L,beta,N) = log[S_pin/Q_D].                             (16)

The positive diagonal history-register operator diag(q_sbold) obeys the same Hilbert-Schmidt inequality as in (3), so

    Q_L^2 <= d^N S_pin.                                       (17)

The exact partition-function identity and (11), (15), (17) give

    R_L = 2 epsilon_L-epsilon_D + log(Q_L^2/Q_D),
    R_L <= 3N log d + G_(L,beta,N).                            (18)

The sign of epsilon_D is favorable, so it is discarded rather than paid again. This is a valid, explicit reduction using only positive histories.

[False candidate] The intended positivity step is S_pin <= Q_D, or G <= 0. It is false on all three calibration models at tau = L/M with N = 1, as recorded in CHECKS.md. It also fails on XX if Q_D is replaced by the larger unmeasured Z_(2L) once the proportionality constant a is large enough. Entrywise nonnegative imaginary-time evolution does not establish this comparison.

[Proved distinction] The squared q weights describe two independently evolving L-rings with the same recorded labels. They do not describe the reconnected ring with every intervening microscopic seam variable fixed equal. The latter would require access to the unrecorded seam evolution between the coarse slices. A local measurement record is not the full spatial cut of a Euclidean tensor network.

## 5. Why stoquasticity does not supply the missing separator

[Proved local witness] Consider the three-site open Hamiltonian

    H = -J (X_A X_c + X_c X_B),       J > 0,

where X is a Pauli matrix and the seam c is measured in the Z basis. H is stoquastic in that basis, and its two bond terms commute. Nevertheless, conditioning on seam state 0 at both ends of an interval t gives

    <0|exp(-tH)|0>_c
      = cosh[tJ(X_A+X_B)]
      = cosh(tJ)^2 I_A tensor I_B + sinh(tJ)^2 X_A tensor X_B.  (19)

The last operator has operator Schmidt rank two for t > 0, since I and X are linearly independent on each side and both coefficients are nonzero. It does not factor into separate left and right evolutions. Positivity of its matrix entries does not remove this conditional coupling. This model does have a classical X basis; the witness refutes a basis-independent inference from a recorded label to spatial factorization. For the XX chain, seam couplings contain both X_c and Y_c. Any seam operator commuting with the Hamiltonian must commute with both, so no nontrivial one-site projective basis makes the evolution diagonal at all times.

[Proved history accounting] In a finite stoquastic Trotter network let x denote a complete microscopic seam history, let B be the nonnegative spatial propagation matrix across width L, and group histories by their sampled coarse record c(x). Then

    a_s = sum_(x:c(x)=s) B_xx,
    sum_s a_s^2 = sum_(x,y:c(x)=c(y)) B_xx B_yy,
    Z_(2L) = sum_(x,y) B_xy B_yx.                              (20)

The microscopic classical positivity proof compares sum_x B_xx^2 to the last sum. The grouped expression in (20) contains additional terms B_xx B_yy with x != y. Nonnegativity supplies no comparison of these terms with B_xy B_yx. For instance, B = I_2 and grouping its two labels together gives a squared grouped trace four while Tr B^2 = two. The local chain content needed to exclude this behavior at the selected coarse scale is a further estimate.

[Proved and numerical] On even rings, the alternating Z rotation changes antiferromagnetic XX and Heisenberg to stoquastic Hamiltonians. It commutes with the measured Z projectors, so all the coherent and pinched trace functionals here are unchanged. Consequently the XX and Heisenberg failures of the zero-cost gluing step are stoquastic counterchecks. SU(3) supplies a model with no stoquastic assumption. Stoquasticity removes the sign issue in (1), but not the missing between-slice histories.

## 6. An analytic failure of zero-cost gluing and the Zeno limit

### 6.1 A finite XX witness

[Proved] For even XX rings, E_L = -csc(pi/L), as recorded in C118/C120. At L = 4,

    Delta_4 = 2 sqrt(2)-csc(pi/8) = 0.2153011949... > 0.

The elementary ground-state and dimension bounds give

    Z_4 >= exp(-beta E_4),
    Z_8 <= 2^8 exp(-beta E_8),
    R_4 >= beta Delta_4-8 log 2.                               (21)

Take beta = 64, N = 1, M = 1/2, so tau = 8 L/M. The right side of (21) is greater than 8.23, and hence exceeds 3 log 2. All numerical inequalities in this paragraph follow directly from the displayed elementary expressions, with margins greater than six.

[Proved] For the pinched histories, q_0 = q_1 by spin-flip symmetry, so S_pin = Q_L^2/2. By (11), Q_L >= Z_L/2, and Q_D <= Z_8. Therefore

    S_pin/Q_D >= Z_4^2/(8 Z_8) = exp(R_4)/8 > 1.               (22)

Even S_pin/Z_8 > 1. The coherent version follows from (5). Thus no-cost gluing fails at a spacing proportional to L/M, already in a stoquastic chain. This does not exclude a theorem with a sufficiently small universal a or a finite cost per slice.

### 6.2 Continuous measurement changes the Hamiltonian

[Proved, fixed finite system] Fix beta and send N to infinity in (6). Between pinchings,

    A_tau(X) = X - (tau/2){H,X} + O(tau^2).

On the range of D, D(X) = X and the bimodule property of pinching gives

    D({H,X}) = {D(H),X}.

The product formula on this finite-dimensional operator space therefore converges to the evolution generated by -{D(H), .}/2. Endpoint half intervals vanish in the limit. Starting from I and taking the trace yields

    lim_(N -> infinity) Q_(H,N) = Tr exp[-beta D(H)].           (23)

Thus continuously classicalizing the seam is a different operation from keeping a finite number of coarse records.

[Proved, XX] Z-basis pinching on one site deletes its two adjacent XX hopping bonds. The site becomes isolated and the other L-1 sites form an open XX chain. Its single-particle energies are cos(pi j/L), j = 1,...,L-1. For even L, summing the occupied levels gives

    E_Zeno = 1/2 - (1/2) cot(pi/(2L)),
    E_Zeno-E_L = 1/2 + (1/2) tan(pi/(2L)) -> 1/2 = M.          (24)

The last equality follows from csc(2x) = [cot x + tan x]/2. The isolated spin and the open-chain zero mode give total ground multiplicity four, while the even periodic XX ground state is unique. At fixed L,

    log Z_L-log Z_Zeno = beta(E_Zeno-E_L)-log 4+o(1)            (25)

as beta tends to infinity. This costs order beta M. At beta proportional to L the same conclusion follows from the explicit full traces, since the shifted excitation traces are finite functions of beta/L. Sparse pinching has bounded cost but leaves seam propagation. Continuous pinching removes that propagation but has the extensive-in-beta seam-energy cost.

## 7. Exact large-L free-fermion check of one coarse pinch

[Proved] For one midpoint pinching in the XX Z basis, write A = exp(-beta H_L/2), and Z_c = I-2 n_c. Since D(X) = (X+Z_c X Z_c)/2,

    Q_(L,1) = Tr[D(A)^2]
            = (1/2)[Z_L + Tr(A Z_c A Z_c)].                   (26)

Define C_L = Tr(A Z_c A Z_c)/Z_L. It is nonnegative because Tr(A Z_c A Z_c) = Tr[(A^(1/2) Z_c A^(1/2))^2] >= 0, and is at most one by Hilbert-Schmidt Cauchy-Schwarz. Hence Q/Z = (1+C_L)/2.

[Proved finite-grid formula] Keep both physical Jordan-Wigner parity sectors. For s = 0,1/2, let eps_(s,j) = cos[2 pi(j+s)/L], and put

    P_s = product_j (1+exp(-beta eps_(s,j))),
    g_s = sum_j 1/[2 cosh(beta eps_(s,j)/2)].

Then the full spin trace Z_L is the C118 parity-projected expression, and

    C_L = [2/(L^2 Z_L)] sum_(s=0,1/2) P_s g_s^2.              (27)

The factors P_s/Z_L are evaluated with their filled-sea energies removed in the code. This avoids subtracting bulk quantities of order beta L.

[Proved derivation] On either unprojected grid the half-filled Gaussian trace has <n_c> = 1/2. Wick's identity applied to the two density insertions at imaginary-time separation beta/2 gives

    Tr(A Z_c A Z_c)_grand / P_s = 4 g_s^2/L^2.                 (28)

One can verify Wick's step directly. Write n_c = L^(-1) sum_(k,q) c_k^dagger c_q. The diagonal part contributes <N_f^2>/L^2. The off-diagonal terms give L^(-2) sum_(k != q) exp[(beta/2)(eps_k-eps_q)] <n_k(1-n_q)>. Combining the diagonal variance with this sum yields L^(-2)[sum_k exp(beta eps_k/2) f_k][sum_q exp(-beta eps_q/2)(1-f_q)] = g_s^2/L^2; the mean-density part cancels against the constant and linear terms in Z_c Z_c.

[Proved parity contribution] The same calculation with a parity insertion uses the unnormalized factor product_j(1-exp(-beta eps_j)) and formal occupations f_j = -x_j/(1-x_j), x_j = exp(-beta eps_j). When there are no zero modes, energies occur in opposite pairs, so sum_j f_j = L/2 and sum_j sqrt(x_j)/(1-x_j) = 0. Both terms in the Z_c correlator vanish. If the grid contains zero modes, there are two. A linear combination of their plane waves has zero amplitude at site c. It is fixed by both exp(-beta h/2) and the one-particle reflection implementing Z_c. Therefore the parity-inserted Gaussian trace det(I-exp(-beta h/2) R_c exp(-beta h/2) R_c) has a zero eigenvalue and vanishes. Thus each physical parity projection supplies one-half of (28). Summing the two physical grids proves (27), including by continuity at beta = 0.

[Numerical] Equation (27), including both parity sectors, agrees with direct ED for L = 4,6,8 and multiple temperatures to below 4e-15. Large-ring computations extend to L = 8192 and paired length 16384 for R_L. The one-pinch loss approaches log 2 at fixed nonzero beta/L as L grows. At L = 8192 and beta = 4L, R_L = 3.066133108 and log[S_pin/Z_(2L)] = R_L + 2 log(Q_L/Z_L)-log 2 = 0.986691567 > 0. Thus failure of zero-cost gluing is not confined to small rings. No large-L interacting channel computation or multi-pinch free-fermion determinant is claimed.

## 8. The precise remaining lemma and conditional theorem

[Open lemma P] There exist explicit a > 0 and finite nonnegative b_0(d), b_1(d), independent of h, L, beta and every crossover scale, such that the positive-history functional (16), with N chosen by (13), satisfies

    sum_sbold q_(L,sbold)^2
      <= exp[b_0(d)+b_1(d) N] Q_D.                            (29)

One may specify a basis-selection rule depending on h as part of the lemma. If that rule is used, it must be the same rule on the L and 2L systems. The stronger b_0 = b_1 = 0 claim is false at a = 1 in the reported examples. No conclusion excluding all sufficiently small a is asserted.

[Proved conditional theorem] If lemma P holds, (18) and N <= 1+beta M/(aL) imply

    R_L <= b_0(d)+3 log d+b_1(d)
           + [(3 log d+b_1(d))/a] M beta/L.                   (30)

Consequently, under the inherited vacuum Casimir expansion,

    A <= [4/(pi a)] [3 log d+b_1(d)] M.                       (31)

This identifies all constants. The pinching error has already been proved and is not part of the open lemma.

[Open alternative lemma C] A smaller sufficient constant would follow directly from the coherent histories if, at (13),

    sum_sbold |a_(L,sbold)|^2
      <= exp[c_0(d)+c_1(d) N] Z_(2L).                         (32)

Then (3) gives K = c_0+log d+c_1 and k f(d) = (log d+c_1)/a. For stoquastic models the absolute values can be dropped. Equation (20) states the corresponding conditional microscopic-history comparison. The absence of a sign problem alone does not prove (32).

[Proved equivalence up to the known insertion cost] Since q_s >= 0, one also has S_pin <= Q_L^2. Apply Q_L <= Z_L and Q_D >= d^(-2N) Z_(2L) to obtain G <= R_L+2N log d. Together with (18),

    R_L-3N log d <= G <= R_L+2N log d.                         (33a)

Thus lemma P is comparable in strength to the original target once N is chosen as in (13). It is a concrete unnormalized reconnection problem, not an independent bound on that problem. The pinching theorem alone does not control its remaining slope.

[Proved diagnostic] At fixed finite L with a unique maximally mixed one-site ground marginal, the cold-slab limit of (12) gives S_pin = exp(-2 beta E_L) d^(-3N)(1+o(1)). If r_st = <psi_(2L)|P_s^(0) P_t^(L)|psi_(2L)> and u = sum_(s,t) r_st^2, then Q_D = exp(-beta E_(2L)) u^N(1+o(1)). Hence

    G = beta Delta_L-3N log d-N log u+o(1).                   (33)

The same argument on coherent histories gives log(S_coh/Z_(2L)) = beta Delta_L-N log d+o(1). These formulas are fixed-N, cold-slab limits. They expose the retained doubling-energy term; they are not uniform estimates at fixed tau/L or proofs of (29).

## 9. Markov and locality inputs

[Assumed published context] Kato and Brandao prove approximate Markov properties for spatial tripartitions of one-dimensional Gibbs states, with conditional mutual information decaying with the square root of the buffer length. This is a statement about physical spatial regions, not an independent tensor factor for each Euclidean-time sample. See [Quantum Approximate Markov Chains are Thermal](https://arxiv.org/abs/1609.06636). Applying such a theorem to (29) would require a construction of the time-history state and a multiplicative comparison of its unnormalized contractions. Neither follows from the existence of a spatial recovery channel.

[Assumed published context and inference] Perez-Garcia and Perez-Hernandez distinguish real- and complex-time locality estimates, and describe Araki's finite-range complex-time localization region as growing exponentially with the modulus of time. See the introduction and Section 2 of [Locality Estimates for Complex Time Evolution in 1D](https://pmc.ncbi.nlm.nih.gov/articles/PMC10113381/). Thus a real-time light cone cannot be analytically continued into a claim that coarse Euclidean propagation over tau proportional to L/M leaves only a single classical seam label. This observation does not preclude a stronger trace-specific estimate.

[Proved limitation of the level-spacing premise] The velocity bound v_a <= 4 e M is an upper bound. It gives no lower bound of order M/L on all finite-size gaps. The slow-sector examples recorded in the background have arbitrarily small velocities, and the number of infrared sectors is not bounded by d. Thus a universal seam level spacing of order M/L is not an available hypothesis. Fortunately it is unnecessary for (11). The missing use of locality is instead (29) or (32), where the full finite-size reconnection information remains.

[Scope] The general velocity-weighted bound remains open. The new unconditional theorem is the arbitrary-spacing pinching estimate (11), with an explicit positive-history reduction (18). The zero-cost spatial gluing step is disproved, including for stoquastic XX. All numerical sizes, failed candidates and reproduction commands are listed in CHECKS.md. No heavy local computation, cluster submission, external-file modification or independent review was performed.
