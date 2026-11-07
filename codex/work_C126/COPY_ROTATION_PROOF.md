# C126 copy-rotation proof and obstruction

2026-11-06. This file records an independent derivation in a directory with concurrent edits. Claims marked proved refer to the stated hypotheses. No independent referee review has been performed.

## 1. Seam replacement and the two limits

[Proved] Label sites by (x,r), with 0<=x<L and 0<=r<N. The untwisted Hamiltonian is

    H0 = sum_r sum_(x=0)^(L-2) h_((x,r),(x+1,r))
         + sum_r h_((L-1,r),(0,r)).

Replace the final sum by sum_r h_((L-1,r),(0,r+1)), with copy indices modulo N. Traversing its graph visits

    (0,0),...,(L-1,0),(0,1),...,(L-1,N-1),

then returns to (0,0). Tensor-factor relabeling therefore identifies Htw with H_(NL), including the orientation of a nonsymmetric Hermitian bond. In particular

    min spec Htw - min spec H0 = E_(NL)-N E_L.                 (1)

The seam term is conjugated by the on-site copy cycle at x=0. The whole Hamiltonian is not conjugated. Whole-Hamiltonian conjugation would change both incident bond sums and preserve the spectrum.

[Proved] Existence of the thermodynamic energy density gives, at fixed L,

    lim_(N->infinity) [E_(NL)-N E_L]/N = L e0-E_L.            (2)

There is no uniform-CFT approximation in this limit. The statement follows directly from E_n/n -> e0.

[Assumed for the CFT check] On the selected vacuum sequence and its fixed integer multiples,

    E_L = e0 L - pi A/(6L) + o(1/L),
    A = sum_a c_a v_a.

The sequence is selected before estimating the defect. For the numerical models it is even L for spin-1/2 antiferromagnetic XXZ and L divisible by 3 for the SU(3) fundamental swap. Every L>=3 is compatible in the ferromagnetic Ising example below. No lower bound on the onset of CFT scaling is assumed.

[Proved conditional on this expansion] At fixed N,

    E_(NL)-N E_L = pi A/(6L) (N-1/N) + o(1/L).               (3)

It is not an exact lattice identity without the remainder. For N=2 it gives pi A/(4L). The large-N limit of the continuum expression divided by N is pi A/(6L), agreeing with the Casimir deficit in (2) after subsequently taking large L.

[Proved within CFT] Uniformization w=z^(1/N) has Schwarzian

    {w,z} = (1-1/N^2)/(2z^2).

The stress tensor on each sheet has expectation (c/12){w,z}. Summing N sheets gives h_N/z^2 with h_N=c(N-1/N)/24. The antiholomorphic computation is identical. The cylinder excitation energy above the N-copy vacuum is

    (2pi v/L)(h_N+hbar_N) = pi c v/(6L)(N-1/N).

Sector energies add when their velocities differ. This is the branch-point twist construction of [Calabrese and Cardy](https://arxiv.org/abs/hep-th/0405152) and [Cardy, Castro-Alvaredo and Doyon](https://arxiv.org/abs/0706.3384).

## 2. A finite-size quasi-free bound, including pairing

[Proved] The bound below applies to nearest-neighbor even quadratic fermion Hamiltonians with finitely many modes per site, and to their tensor-site spin realizations with the usual Jordan-Wigner parity boundary condition. On a fixed fermionic boundary condition independent of parity, the unsigned odd-N version alone suffices for the Casimir conclusion. A tensor-spin copy permutation must not be identified with a graded fermionic permutation as an operator.

Write the Majorana form with an even number of Majoranas per physical site as

    h_hop(x,x+1) = (i/2) gamma_x^T B gamma_(x+1),
    onsite term = (i/4) gamma_x^T A0 gamma_x,

with real B and antisymmetric real A0. All matrices in the following copy rotations act on the copy index, not the internal Majorana index. The sum of N identical quadratic Hamiltonians is invariant under O(N) rotations of the copies. With conserved particle number it has the larger U(N) symmetry. General pairing does not preserve U(N).

[Proved, using the bond-scale lemma of REPORT 9.71] Put t=||h_hop||. Then

    t = ||B||_1/2 <= M.                                      (4)

In number-conserving notation t=||B_complex||_1. The gauge-independent lower bound follows from two extremal hopping dimer states with identical maximally mixed one-site marginals. Their energy difference is 2t. On-site terms and telescoping gauge terms have zero difference between these states, so every representative has half-spread at least t. This is the inherited lemma, not a new estimate of M.

[Proved] Choose a ground vector psi_L with definite fermion parity p=+1 or -1. In the spin realization the equivalent quadratic ring has boundary sign eta=-p. Let Omega be its N-copy graded product. Its energy is N E_L and its total parity is p^N. Cross-copy two-point functions vanish because each copy has definite parity, while the diagonal copy covariances are identical.

Define R_N to be a real signed cyclic shift whose last link has sign (-1)^(N-1). Thus

    R_N^N = (-1)^(N-1) I,   det R_N=1,
    eigenangles theta_k = (2pi/N)[k-(N-1)/2],  0<=k<N.       (5)

All angles lie strictly between -pi and pi and occur in opposite pairs, with a zero angle for odd N. Therefore R_N lies in SO(N) and has a real antisymmetric logarithm K. The special role of the sign is essential for even N.

Insert R_N into the copy seam of the quadratic Hamiltonian, retaining eta at every original closing bond. Unfolding gives a single NL ring with total boundary sign

    eta^N (-1)^(N-1) = (-p)^N (-1)^(N-1) = -p^N.           (6)

This is the physical spin boundary condition for the conserved total parity p^N. Intermediate signs along the unfolded ring can be removed by c_j -> -c_j on appropriate intervals, which also works for pairing. Consequently a trial vector in this parity sector gives an upper bound on the unrestricted physical E_(NL). For odd N and eta=+1 or -1 fixed independently of parity, the unsigned cyclic shift instead has eta^N=eta and gives the same fixed-boundary long ring.

[Proved] Spread the rotation with O_x=exp(xK/L), x=0,...,L-1. Its even Gaussian implementers preserve total parity. Conjugating the seam Hamiltonian by the product of these implementers puts exp(K/L) on every bond, because

    O_x^T O_(x+1)=exp(K/L),
    O_(L-1)^T R_N O_0=exp(K/L).

On-site quadratic terms are invariant. Denote the one-copy average hopping expectation, including its boundary sign, by e_hop. Factorization of copy covariances makes the hopping expectation on the N-copy ring proportional to tr exp(K/L). Hence the trial excess is exactly

    T_(N,L) = L e_hop [sum_k cos(theta_k/L)-N].              (7)

This formula requires no assumed stiffness limit and no nondegeneracy. If a translation choice is needed, choose a simultaneous ground and translation vector; alternatively e_hop is the average over the ring. The current terms cancel between opposite angles. Only this cancellation and identical copy covariances are used.

[Proved] Since |e_hop|<=t<=M, 1-cos u<=u^2/2, and

    sum_k theta_k^2 = pi^2 (N-1/N)/3,

the variational principle gives

    E_(NL)-N E_L <= T_(N,L)
                  <= pi^2 M/(6L) (N-1/N).                  (8)

The norm estimate holds even if e_hop is positive. In that case the trial excess in (7) is nonpositive. This finite-size statement does not assume a conformal regime. For the standard spin parity realization it holds for all integer N>=2. At fixed fermionic periodic or antiperiodic boundary condition, use odd N if the even-N boundary convention is not part of the model.

[Proved] In the spin realization, N=2 gives

    L [E_(2L)-2 E_L] <= (pi^2/4) M.                         (9)

For either realization, divide (8) by N and let N tend to infinity through the allowed integers. Equation (2) gives

    E_L >= e0 L - pi^2 M/(6L).                              (10)

Equations (3) and (8), or the Casimir limit of (10), imply

    A <= pi M.                                              (11)

[Proved comparison with inherited results] The sharp asymptotic bound A<=2M in REPORT 9.71 is stronger than (11). Its lifted Fermi-vector trace lemma uses more spectral information. The copy-spreading estimate instead supplies a uniform finite-length bound with C=pi^2/4 in (9). C119's coefficient 10pi/3 belongs to the all-temperature estimate R_L(beta)<=3 log d+(10pi/3)M beta/L for its specified free families. Equation (9) improves that coefficient only at zero temperature and does not prove an all-temperature bound. The sharp free asymptotic doubling coefficient is pi/2, inherited from A<=2M, not pi^2/4.

[Proved normalization caution] The stiffness budget sum_alpha D_alpha<=4M in the pairing contribution uses Hilbert-Schmidt normalized conserved charges in the physical on-site representation. A replica rung has dimension d^N and its bond scale is NM. The unnormalized copy-rotation generators cannot be inserted into that theorem using the original M. Equation (7) evaluates their total cost directly, avoiding representation-index and normalization ambiguities. In the interacting case the needed generators are not conserved in the first place.

[Proved check] For the half-filled XX spin chain, E_L=-csc(pi/L) on even rings and e_hop=E_L/L. Equation (7) at N=2 becomes

    T_(2,L) = -2 E_L [1-cos(pi/(2L))]
            = 2 csc(pi/L)-csc(pi/(2L))
            = E_(2L)-2E_L.

Thus the trial is exact in this example. All factors of two match the M=1/2 convention.

## 3. A ground-state stabilizer obstruction and an orbit obstruction

[Proved] There is already an obstruction to an exact continuous ground-state symmetry on the even Heisenberg rings L>=4. Their singlet ground vector has two-site marginal

    rho_xy = [I tensor I + r sum_(a=1)^3 sigma_a tensor sigma_a]/4,
    r = 4 E_L/(3L),   -1<r<0.

For two copies, use the on-rung operator basis sigma_a tensor sigma_b, with sigma_0=I. The two-rung correlation matrix is

    T_copy = diag(1,r,r,r) tensor diag(1,r,r,r).

Its eigenvalues are 1, r of multiplicity 6, and r^2 of multiplicity 9. The r eigenspace is the Lie algebra

    span{sigma_a tensor I, I tensor sigma_a} = su(2) + su(2),

with an irrelevant factor of i in the anti-Hermitian convention. If U^tensor L preserves psi_L^tensor 2 up to phase, its adjoint action on the rung operator space must preserve T_copy and hence this eigenspace. It therefore induces an automorphism of the direct sum of two simple Lie ideals. A continuous automorphism path starting at the identity preserves each ideal separately. The copy swap exchanges them. Thus no continuous path of uniform on-site unitaries that exactly fixes this ground vector can connect identity to copy swap. This includes arbitrary U(4), not just isotypic phases. At L=2, r=-1 and the spectral argument needs separate treatment; it is not used.

[Proved scope] This disproves an exact ground-sector symmetry path even for a unique critical singlet vacuum. Approximate low-energy paths would require quantitative control of their leakage and of the energy of spatially varying rotations. The stabilizer argument alone is not a lower bound on every such approximate construction.

[Proved] The same argument works at every finite N. The correlation matrix is T^tensor N, with eigenvalues r^k. Since 0<|r|<1, the r eigenspace consists only of the 3N one-copy Pauli operators. It is su(2)^direct-sum N. A continuous normalizer path cannot permute its simple ideals, whereas an N-cycle does. Increasing N does not create an exact connected ground-state stabilizer containing the copy cycle.

[Proved theorem] Let H_L=-J sum_x Z_x Z_(x+1), J>0, with L>=3, and let

    psi_L = (|0>^L+|1>^L)/sqrt(2),
    Omega = psi_L tensor ... tensor psi_L   (N copies).

The on-site dimension of the original chain is d=2 and M=J. For every N>=2,

    inf_(U_0,...,U_(L-1))
      <Omega| U^dag Htw U |Omega> - N E_L = JN,
    U = tensor_x U_x,   U_x in U(2^N).                     (12)

The true twist cost is zero, E_(NL)-N E_L=0. All lengths in this example are vacuum compatible. No restriction to a particular logarithm, isotypic decomposition, or smooth profile appears in (12).

[Proved, step 1] Regrouping by rungs, put D=2^N and write

    Omega = D^(-1/2) sum_(a in {0,1}^N) |a>^L.

After U, every two-rung marginal is

    D^(-1) sum_a |U_x a><U_x a| tensor |U_y a><U_y a|.

Tracing at least one other rung kills all unequal a terms. This explains the hypothesis L>=3.

[Proved, step 2] The Ising energy depends on U_x only through

    W_x(b,a)=|<b|U_x|a>|^2.

Each W_x is doubly stochastic. The energy is affine in each W_x separately when all other matrices are fixed, since no bond has both ends at the same rung. Enlarge the domain to all doubly stochastic matrices. A minimum of this multilinear function can be chosen at an extreme point in every variable, by replacing one variable at a time at a global minimum. The extreme points are permutation matrices. They are allowed unitary matrices, so this relaxation is tight.

[Proved, step 3] It remains to minimize over bijections f_x of the N-bit strings. For a label a, the energy excess is twice J times

    sum_(x=0)^(L-2) dist(f_x(a),f_(x+1)(a))
          + dist(f_(L-1)(a),P f_0(a)),

where dist is Hamming distance and P cycles the N coordinates. The triangle inequality bounds this below by dist(f_0(a),P f_0(a)). Averaging uniformly over a removes the bijection f_0. For N>=2, each neighboring coordinate pair disagrees in half the strings, so

    D^(-1) sum_b dist(b,Pb) = N/2.

The lower bound is therefore JN. Taking every f_x to be the identity attains it. This proves (12).

[Proved consequence for route (a)] A continuous on-site path from identity to the copy cycle cannot in general be made into an O(M/L) twist trial merely by requiring symmetry on the relevant ground sector. Equation (12) includes all products of on-site unitaries and starts from an exact permutation-invariant ground vector of H0. Any proposed low-energy path theorem valid for every such vector is false. For S_2, the central isotypic path is, up to an overall phase, exp(i theta S_x), and is a special case.

[Proved scope qualification] The theorem does not exclude choosing a different vector in a degenerate ground space. Starting with |0>^L in every copy makes the twist free. Nor does it exclude nonunitary filters, coherent sums of different rung-unitary products, or operations coupling distinct rungs. Thus it disproves a universal on-site-unitary construction for the specified starting vector, not the desired energy bound or every possible choice of initial vacuum.

[Proved consequence for route (b)] Dividing (12) by N leaves J, independent of N. Large N does not repair this obstruction. The copies are already an exact product as whole chains, and the state is exactly exchangeable in the copy direction. Thus exchangeability or a de Finetti description of fixed-copy marginals cannot by itself imply the required energy estimate. Such information is present in a counterexample to that inference. A cycle uses all N copies.

[Proved distinction] An on-site reduced state rho_x^tensor N commutes with every on-site copy permutation and with its logarithm. This does not imply invariance of the two-rung state under the sum of the two on-site logarithms. In general

    [rho_(xy)^tensor N, sigma_x sigma_y]=0

does not imply

    [rho_(xy)^tensor N, log(sigma_x)+log(sigma_y)]=0.

The logarithm is not additive under the tensor representation. The Ising theorem and the numerical contractions provide explicit violations of the proposed energy consequence. For fermionic Gaussian copies, quadratic copy covariance supplies the additional continuous invariance used in section 2.

## 4. Finite-group LSM and ordered vacua

[Proved] The on-site S_N action on (C^d)^tensor N is an ordinary representation. It does not by itself carry the projective representation per unit cell used in discrete-symmetry LSM obstructions. A one-site vector |v>^tensor N is invariant, and a tensor product of such vectors is an allowed trivial symmetric ground state. There is no replica-only LSM anomaly that forces a small excited-state energy.

[Literature input, assumed only as context] The discrete-symmetry twisted-boundary argument of [Yao and Oshikawa](https://arxiv.org/abs/2010.09244) establishes degeneracy under a projective on-site representation and uses robustness of gapped phases to connect boundary conditions. It does not supply a universal O(M/L) upper bound on a permutation-twist energy. This route would need an additional quantitative estimate.

[Proved] For a classical ordered vacuum with period p dividing L, put the same ordered pattern in every copy. Its tuple of order labels is fixed by the copy cycle, so the seam does not force a domain wall. In particular the Ising twisted ground state has zero cost on all lengths. More generally any identical extremal phase omega^tensor N is fixed by copy permutation, so permutation holonomy alone does not force an interface between distinct phases.

[Proved distinction] Vacuum compatibility does not prevent a domain-wall cost in the specified trial orbit. Equation (12) has cost JN on compatible lengths. A symmetric cat in each copy contains unequal order labels among copies, and on-site unitaries cannot remove their cost. The antiferromagnetic Ising version follows by a staggered Z-basis flip on even L. It has the same obstruction on precisely the period-compatible sequence.

[Assumed, not a uniform theorem] In a fixed gapped quantum ordered phase one expects unrestricted finite-size defects to become small once L exceeds the correlation length. No bound uniform over that correlation length follows here. The absence of an order-parameter interface does not prove the C(d) estimate in the crossover regime. No physical-chain counterexample to the target is supplied.

## 5. Fractional swaps on the interacting ED examples

[Proved] For two copies define S_x as the ordinary tensor swap and U_x(theta)=cos(theta)I+i sin(theta)S_x. Let rho be the two-site marginal of a real translation ground vector, R the corresponding two-copy two-rung marginal, and B=h_A+h_B. Both R and B commute with T=S_x S_y. On its T=+1 and -1 sectors, U_x(a)U_y(b) reduces to exp(i(a+b)S_x) and exp(i(a-b)S_x), respectively. Reality eliminates the sine-cosine cross expectation. Therefore

    F(a,b) = tr[R U^dag B U]
           = 2e_L + k_plus sin^2(a+b) + k_minus sin^2(a-b).  (13)

The coefficients can be obtained without fitting,

    k_plus  = F(pi/4, pi/4)-2e_L,
    k_minus = F(pi/4,-pi/4)-2e_L,
    k_plus+k_minus = 2[b_L-e_L],
    b_L=tr[h (rho_x tensor rho_x)].

This is an exact two-rung calculation even though the full state is entangled. Other rung unitaries cancel from each local expectation.

[Proved] For the permutation-twisted ring set theta_L=theta_0+pi/2. The full trial excess is

    T(theta)=sum_(x=0)^(L-1)
       [k_plus sin^2(theta_x+theta_(x+1))
        +k_minus sin^2(theta_x-theta_(x+1))].               (14)

The uniform ramp theta_x=pi x/(2L) gives

    T_ramp = (L/2) k_plus + L k_minus sin^2(pi/(2L)).        (15)

Thus the missing symmetry appears as a bulk potential, not just a gradient stiffness. A uniform intermediate angle costs k_plus sin^2(2theta) per untwisted bond.

[Proved] When both coefficients are positive, every profile, including the global optimizer, obeys

    T(theta) >= 2 min(k_plus,k_minus).                      (16)

Indeed each bond is at least m[1-cos(2theta_x)cos(2theta_(x+1))], where m is the smaller coefficient. Put z_x=cos(2theta_x), with z_L=-z_0. The sum of 1-z_x z_(x+1) is multilinear on the cube [-1,1]^L. A minimum occurs at vertices. Antiperiodicity forces at least one disagreeing edge, giving a minimum of 2. This proves (16). If the coefficients stay bounded away from zero, no optimization of this family can have a cost tending to zero.

[Numerical] All nontrivial ED rows L>=3 have positive coefficients. The uniform-ramp ratios grow rapidly. The optimized ratios also grow, despite allowing arbitrary real angles and profiles outside the interval [0,pi/2]. Optimization values are upper bounds on the best family cost, not global certificates; (16) is the analytic certificate against vanishing cost when the limiting coefficients remain positive.

[Proved, uniform Heisenberg obstruction] For the singlet Heisenberg ground state put a=-e_L. Spin contractions in (13) give

    k_plus=a+4a^2/3,  k_minus=a-4a^2/3.

A dimer-product variational state gives e_L<=-3/8. The three-site operator S_1.S_2+S_2.S_3 has minimum -1, and summing its translates gives 2H_L, so e_L>=-1/2 for even L>=4. Therefore 3/8<=a<=1/2 and k_minus>=1/6. Equation (16) gives

    inf_theta T(theta) >= 1/3,  even L>=4.                  (18)

This excludes an O(M/L) upper bound within the entire fractional-swap family without relying on the optimizer or on a limiting coefficient. Conditional only on the inherited Heisenberg Casimir asymptotics D(L)~pi^2/(8L), its ratio to D(L) diverges at least as 8L/(3pi^2).

[Proved, inherited comparator] C114 gives the coherent seam state sum_j phi_j, where phi_j=(product_(x<j) S_x)Omega. Its trial excess is

    T_coh = q_L/S_L,
    q_L=2(b_L-e_L),  S_L=sum_(r=0)^(L-1) tr(rho_[0,r)^2).

For real bonds and ground vectors this is the optimum over the entire seam-position span. This established result is used for comparison, not offered as a new route. Its failure to achieve the target on XX was already proved in C114. In the Ising cat example S_L=(L+1)/2 and q_L=2J, so T_coh=4J/(L+1). Coherent recombination escapes the obstruction (12).

[Heuristic, CFT scaling] If interval second-Renyi entropies have their usual fixed-point form, tr rho_r^2 scales as [L sin(pi r/L)]^(-c/4), apart from endpoint and logarithmic corrections. For c<4, S_L then scales as L^(1-c/4), so T_coh/D grows as L^(c/4) in a critical chain with nonzero limiting q_L. This explains the slow numerical growth for c=1 and c=2. It is not a uniform lattice proof. In a gapped ordered chain an increasing ratio to an exponentially small D does not by itself rule out an O(M/L) upper bound. The ratio diagnostic must be interpreted with this distinction at Delta=2.

## 6. What remains open

[Proved] The tensor-copy identity and its CFT normalization are valid. The continuous copy-rotation proof yields A<=pi M, a finite-size C<=pi^2/4 in the standard quasi-free spin realization, and the Casimir variant with K=pi^2/6. The sharper A<=2M is inherited. A universal on-site rung-unitary construction starting from every exact copy-product vacuum is false, for every N, including N tending to infinity.

[Assumed sufficient missing lemma] Use the same integrated seam estimate as PROOF.md section 9. Set H(s)=H0+sV, where V is the seam replacement, epsilon(s)=min spec H(s), q_*=epsilon'_+(0), and mu=-epsilon'' on (0,1). The single required estimate is

    integral_(0,1) (1-s) dmu(s) >= q_*-B(d) M/L.             (17)

[Proved conditional closure] Concavity and integration of the right derivative give

    E_(2L)-2E_L = q_*-integral_(0,1)(1-s)dmu(s).

Thus (17) gives the target and A<=(4/pi)B(d)M. At a smooth nondegenerate point the measure has density 2 sum_(n>0) |<n|V|0>|^2/(E_n-E_0). At a degeneracy it includes derivative jumps. This formulation permits phase selection and correlated trial states. For example q_*=0 for the Ising ordered ground space even though the prescribed cat-product seam expectation is 2J. The estimate is equivalent to the target after the displayed identity, not a weaker lemma already established by the replica construction. Replacing it by an on-site low-energy-symmetry lemma is invalid by (12).

[Proved classification of the obstruction] The cat-state obstruction is a new microscopic obstruction to this route, rather than one of lessons 1, 3, or 5. Treating the fractional-swap bulk potential as a harmless local response would repeat lesson 2. Calling (17) a solved reduction rather than the missing energy estimate would repeat lesson 4. None of these results establishes or disproves finite C(d) for general interacting chains.
