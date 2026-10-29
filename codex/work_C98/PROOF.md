# C98. Casimir energy, gluing, and two-seam response

2026-10-29. Every new file is in work_C98. Natural logarithms are written as log; log_2 is explicit.

## 1. Scope and normalization

[Assumed] The microscopic Hamiltonian is H_L=sum_j h_(j,j+1) on a periodic chain of local dimension d. Let

    M = inf_o inf_alpha ||h+o tensor 1-1 tensor o+alpha 1||.

The thermodynamic ground energy density is e_0. On a sequence of vacuum-compatible lengths, closed under doubling, assume

    E_L = e_0 L-K/L+o(1/L),       K=pi A/6,
    A=sum_a c_a v_a.

The requested sharp bound is K<=pi^2 M log_2(d)/6. Its natural-log coefficient is pi^2/(6 log 2). The corresponding coefficient for A is pi/log 2.

[Proved caveat] An actual vacuum sequence matters. If the lowest allowed CFT sector has scaling dimension x_a(L), the energy instead includes +2pi sum_a v_a x_a(L)/L. Nonnegative sector offsets can hide A in a lower bound on the raw ground energy. Incommensurate charge rounding can prevent a single coefficient from existing. All coefficient extractions below use the stated vacuum expansion. The finite-dimensional identities do not require it.

[Assumed inherited results] C14 supplies the symmetric SU(N) chain normalization and Bethe/WZW input. C35 proves the logarithmic inequality throughout that family and the zero-field critical XXZ line. These results are used as calibration, not rederived. The dimension-independent alternative A<=C M is already false: at symmetric rank k=N, A/M=pi(N^2-1)/(N H_N) is unbounded. This does not falsify the logarithmic statement.

## 2. What ring doubling would have to prove

[Proved] Define

    D_L=E_(2L)-2E_L.

The vacuum expansion gives

    D_L=3K/(2L)+o(1/L)=pi A/(4L)+o(1/L).                 (1)

Thus a useful variational comparison is an upper bound on D_L. A lower bound or positivity statement has the wrong direction for bounding A from above.

[Proved without a CFT assumption] Put e_L=E_L/L. For a finite integer n,

    e_(2^n L)-e_L = sum_(j=0)^(n-1) D_(2^j L)/(2^(j+1)L).

The thermodynamic limit therefore gives the exact identity

    e_0 L-E_L = sum_(j>=0) 2^(-j-1) D_(2^j L).          (2)

The series is absolutely convergent because cutting and reconnecting a fixed number of bounded bonds bounds D_L uniformly, as below. If D_n<=B/n for all dyadic n>=L, then

    e_0 L-E_L <= (B/(2L)) sum_(j>=0)4^(-j)=2B/(3L).    (3)

For the target, B=pi^2 M log_2(d)/4. If only the asymptotic statement is needed, (1) suffices; no uniform control of the little-o term across models is asserted.

## 3. Unconditional locality estimates

[Proved] Choose any gauge representative with norm at most M+epsilon, and abbreviate its norm by m. A scalar shift changes E_L and e_0 L equally. A coboundary cancels on every ring. At the end take epsilon to zero.

Choose a translation-invariant ground density matrix rho_L, obtained by translation averaging if necessary. Its one-site marginal is sigma_L and every bond expectation is e_L. Set

    b_L=tr[h (sigma_L tensor sigma_L)].

Tiling an arbitrarily long chain by independent copies of rho_L retains L-1 internal bonds per block and replaces the remaining one by a product bond. Hence

    e_0 <= [E_L-e_L+b_L]/L,
    e_0 L-E_L <= b_L-e_L <=2m.                         (4)

Conversely take a translation-invariant infinite-chain state minimizing the energy density, restrict it to L sites, and close the chain. The first L-1 bonds have energy e_0; the closing expectation is between -m and m. Since |e_0|<=m,

    E_L-e_0 L <=2m.

Consequently

    |E_L-e_0 L|<=2M.                                  (5)

Existence of the minimizing state follows from compactness of local density matrices and translation averaging. A finite-range bounded interaction has the same bulk energy density with open or periodic boundaries, since their energies differ by at most one bond norm in this chosen gauge.

[Proved] A product of two L-ring ground states is a trial state on a 2L ring after replacing two closing bonds. Its energy is

    2E_L+2(b_L-e_L).

Therefore

    D_L<=2(b_L-e_L)<=4M.                               (6)

The reverse operator comparison gives D_L>=-4M: the two topologies differ by four signed local terms, each of norm at most m. These inequalities are gauge optimized by the preceding limiting argument. They do not have a factor 1/L.

[Proved, XX calibration] For h=SxSx+SySy and even L,

    M=1/2,    e_0=-1/pi,    E_L=-csc(pi/L).

The one-site marginal is 1/2, so b_L=0. The trial excess in (6) is

    -2e_L=2/[L sin(pi/L)] -> 2/pi.

The actual excess from (1) is pi/(4L)+O(L^-3), since A=1. Thus the simplest ring gluing misses the desired order even on XX. Averaging over the position of the seam does not improve the total energy, since the state is translation invariant and energy is linear in the density matrix.

## 4. Open-chain closure and finite-size counterchecks

[Proved, free-fermion sum] In the same XX normalization, the even open-chain energy is

    E_open(L)=1/2-[2 sin(pi/(2(L+1)))]^-1.              (7)

The one-particle energies are cos(j pi/(L+1)), j=1,...,L. Filling the negative half and summing the arithmetic progression of cosines gives (7). On the spin ring, the Jordan-Wigner parity fixes the periodic or antiperiodic momentum grid; filling L/2 consecutive momenta gives E_L=-csc(pi/L) for either even-length congruence class.

Expanding these exact expressions yields

    E_ring(L)=-L/pi-pi/(6L)+O(L^-3),
    E_open(L)=-L/pi+(1/2-1/pi)-pi/[24(L+1)]+O(L^-3).

It follows that

    E_ring(L)-E_open(L)
      =-(1/2-1/pi)-pi/(8L)+O(L^-2).                  (8)

The raw cost of closure has a nonzero limiting boundary contribution. An O(M log d/L) bound on that raw cost is false, already at d=2.

[Proved conditional on boundary CFT] If instead one knows

    E_open(L)=e_0 L+e_s-pi A/(24L)+o(1/L)

with the appropriate boundary vacuum, then

    E_ring(L)-E_open(L)+e_s=-pi A/(8L)+o(1/L).         (9)

A lower bound on this subtracted closure cost with coefficient pi^2 M log_2(d)/8 would prove the target. The value e_s and the boundary sector are additional data. A coboundary changes open-chain endpoint terms, so one must use the same bond representative to define both E_open and e_s. Neither a norm bound on one closing bond nor its endpoint dimension estimates the subtracted quantity in (9).

[Proved finite-size counterexample] The sharp inequality with no remainder, for every L, fails for the Heisenberg chain h=S.S. Here M=1/2, e_0=1/4-log 2, and E_4=-2. Therefore

    4e_0-E_4=3-4 log 2=0.227411277760...,
    pi^2 M log_2(2)/(6*4)=pi^2/48=0.205616758356....

The first number is larger. This refutes the finite-L strengthening, not the requested asymptotic statement. ED also gives

    E_8-2E_4=0.348906591063... > pi^2/32=0.308425137534...,

so the sharp no-remainder doubling bound fails too. No conclusion about eventual inequalities at arbitrarily large L is drawn from these small sizes.

## 5. Entropy, Schmidt rank, and MPS gluing

[Proved] A local interaction crossing a cut acts on d by d endpoint spaces. This does not bound the Schmidt rank of the many-body state across that cut by d. A bipartition of an L-site pure state into n and L-n sites has rank at most d^min(n,L-n), which grows with the region size.

[Proved on XX] At half filling, restrict the occupied Fourier orbitals to a contiguous block of L/2 sites. The resulting square Fourier submatrix is an invertible Vandermonde matrix. The analogous restriction of the unoccupied orbitals is also invertible. Hence the block correlation matrix C and 1-C are strictly positive. Every one of its L/2 eigenvalues lies strictly between zero and one. The Gaussian reduced state consequently has all 2^(L/2) eigenvalues nonzero, so

    Schmidt rank =2^(L/2).

The contiguous spin bipartition has the same reduced-state spectrum as the fermionic calculation. This proves the rank claim without interpreting small numerical singular values as exact zeros.

[Proved distinction] In an open-boundary MPS, a cut crosses one virtual bond and the rank is at most chi. In a periodic MPS, a connected bipartition crosses two virtual bonds and the rank is at most chi^2. Neither chi<=d nor chi bounded independently of L follows from locality. Entropy can be much smaller than log rank. An entropy bound by itself also does not supply a particular truncation error at the accuracy required below.

[Proved sufficient condition, not established universally] Let rho_L^(2) be the averaged nearest-neighbor marginal of an L-ring ground state. Suppose one can construct a state sigma_(2L) on a 2L ring whose averaged nearest-neighbor marginal obeys

    ||sigma_(2L)^(2)-rho_L^(2)||_1
         <=b log_2(d)/L^2+o(L^-2).                   (10)

Averaging over translations can make sigma_(2L) translation invariant. Both pair marginals have identical left and right one-site marginals internally, so any coboundary has zero expectation on their difference. Hölder's inequality therefore uses the optimized M and gives

    E_(2L)<=tr[H_(2L) sigma_(2L)]
           <=2E_L+2b M log_2(d)/L+o(1/L).

Combining with (1) proves

    A <=(8b/pi) M log_2(d).                           (11)

The sharp target requires b<=pi^2/8. An energy-weighted version of (10) is weaker, but without a constructive state it merely restates the needed variational inequality. An MPS or purification argument becomes useful only after it supplies the local error O(L^-2) with a coefficient controlled by log d, including its dependence on all preparation or correlation lengths. No such estimate is proved here.

## 6. Why a smooth partial swap is not a free twist

[Proved finite-dimensional obstruction] Put two independent replicas A and B on each site and let S_j swap their two d-dimensional local spaces. The doubled bond h^A+h^B commutes with S_j S_(j+1). It need not commute with S_j+S_(j+1). Discrete copy exchange therefore does not justify distributing a reconnection as a continuously varying symmetry twist.

For a direct calculation, act on both replicas with

    U(theta)=product_j exp(i theta S_j).

Suppose the single-copy nearest-neighbor marginal has zero magnetization and diagonal correlations r_a=<S_j^a S_(j+1)^a>. For a qubit bond h=sum_a J_a S^a tensor S^a, put e=sum_a J_a r_a. Expanding exp(i theta S)=cos(theta)+i sin(theta)S gives

    Delta e_doubled(theta)
      =sin^2(2theta) [-e+4 sum_a J_a r_b r_c],         (12)

where b,c are the two axes different from a.

[Proved derivation of (12)] Under the local rotation, a spin in one copy becomes a linear combination with coefficients cos^2(theta), sin^2(theta) of the two spins, plus or minus 2 sin(theta)cos(theta) times their cross product. The copies have zero one-point functions, so the mixed linear terms and the terms with one cross product vanish. For each component a, contracting two cross products in the product of replica marginals gives 2 r_b r_c. Summing the two layers leaves

    2e[cos^4(theta)+sin^4(theta)]
      +16 sin^2(theta)cos^2(theta) sum_a J_a r_b r_c.

Subtracting 2e and using 4 sin^2(theta)cos^2(theta)=sin^2(2theta) gives (12).

[Proved conditional on the known bulk correlators] For Heisenberg, rotational invariance gives r_a=e/3 and e=1/4-log 2. Thus

    Delta e_doubled(theta)=q_H sin^2(2theta),
    q_H=-e+4e^2/3=0.704986412078....                   (13)

For XX, r_x=r_y=-1/(2pi) and r_z=-1/pi^2. The latter follows by Wick contraction at half filling. Consequently

    q_XX=1/pi+4/pi^3=0.447316023917....                (14)

These are positive costs per doubled bond at intermediate angles, not gradient energies. For an angle ramp theta_j=pi j/(2L), the sum of the uniform-angle terms is q L/2+O(1). Replacing theta_(j+1) by theta_j in each bond expectation changes it by at most a constant times M |theta_(j+1)-theta_j|; the total error is O(M). Reconnecting a fixed number of seams also changes only O(M) energy. Thus this particular smooth-swap trial state costs order L, not order 1/L. This calculation does not exclude more general entangled gluing states.

## 7. Reflection positivity fixes a sign, not the coefficient

[Assumed published theorem in its stated class] The Lieb-Schupp comparison gives E_LL+E_RR<=2E_LR under the required antiferromagnetic crossing-bond structure. For the open Heisenberg chain it gives 2E_(m+n)>=E_(2m)+E_(2n). This is a model-dependent reflection inequality, not an implication of arbitrary Hermitian locality. A primary account with the hypotheses, the formula, and its parity limitation is [Wojtkiewicz and Skolasinski, arXiv:1410.5322, Sections 2 and 3.1](https://arxiv.org/html/1410.5322). The related reflection structure is in [Lieb and Schupp, math-ph/9910037](https://arxiv.org/abs/math-ph/9910037).

[Proved logical obstruction] Even granting the favorable concavity properties for a ring-energy sequence, they cannot bound K from the microscopic norm ceiling. Consider the abstract sequence

    E_L=e_0 L-B xi/(L+xi),        B>0, xi>0.           (15)

This is an energy sequence used to test an inference, not a constructed Hamiltonian. Its deficit obeys 0<e_0 L-E_L<B. Its energy is concave as a function of continuous L, its energy density increases to e_0, and D_L>0. The deficit is even completely monotone, since all alternating derivatives of B xi/(L+xi) are nonnegative. Nevertheless

    E_L=e_0 L-(B xi)/L+O(L^-2),

and the coefficient B xi is arbitrarily large at fixed B. Setting B=M or 2M respects the elementary norm bounds. Reflection midpoint inequalities, monotonicity, and a uniform O(M) boundary ceiling therefore do not by themselves yield an O(M log d) Casimir coefficient. This does not prove that a reflection-positive Hamiltonian violates the conjecture.

[Proved conditional scale bound] A quantitative hypothesis would repair this inference. For example, if for some L_0

    e_0 L_0-E_(L_0) >= alpha K/L_0,       alpha>0,

then (5) implies K<=2M L_0/alpha and A<=12M L_0/(pi alpha). This proves a logarithmic bound if L_0/alpha<=C log d. It leaves the same uncontrolled onset length encountered in other campaign routes. Fundamental SU(N) vacuum rings have N dividing L, so one must not assume a vacuum-compatible L_0=O(log d) for that family.

## 8. An exact two-seam susceptibility representation

[Proved, finite-dimensional] Start with two L-site rings, with doubled Hamiltonian H_00. Choose two separated cuts x and y of the rung circle. At either cut define V_x by replacing the two straight bonds by the two crossed bonds,

    V_x = h_(xA,(x+1)B)+h_(xB,(x+1)A)
          -h_(xA,(x+1)A)-h_(xB,(x+1)B),               (16)

and similarly for V_y. The orientation of h is kept in every term. After centering and gauge optimization, ||V_x||,||V_y||<=4M. In fact scalar shifts and coboundaries cancel in V_x itself.

Let H(s,t)=H_00+s V_x+t V_y for 0<=s,t<=1. Switching either cut alone produces a single ring of length 2L. Switching both gives two L-site rings, related to H_00 by swapping the replicas on the segment between the cuts. Hence the corner ground energies are

    E(0,0)=E(1,1)=2E_L,
    E(1,0)=E(0,1)=E_(2L).                             (17)

The same statements hold for free energies, with F(0,0)=F(1,1)=2F_L(beta) and F(1,0)=F(0,1)=F_(2L)(beta).

[Proved] At finite beta, define

    F_beta(s,t)=-beta^-1 log tr exp[-beta H(s,t)],
    chi_xy^(beta)(s,t)
      =Re integral_0^beta d tau <delta V_x(tau) delta V_y>_(s,t,beta).

The delta subtracts the Gibbs expectation, and V_x(tau)=exp(tau H)V_x exp(-tau H). Differentiation of the finite-dimensional Gibbs trace gives

    partial_s partial_t F_beta=-chi_xy^(beta).

Integrate over the parameter square and use (17):

    F_(2L)(beta)-2F_L(beta)
       = (1/2) integral_0^1 ds integral_0^1 dt chi_xy^(beta)(s,t).    (18)

For each fixed L, F_beta converges uniformly to the ground energy as beta tends to infinity. The error is at most 2L log d/beta on the 2L-site space. Thus

    D_L=(1/2) lim_(beta->infinity) integral_[0,1]^2 chi_xy^(beta).   (19)

The order of limits is fixed: beta tends to infinity at fixed L before taking L to infinity. Equation (19) remains valid at level crossings. It does not assume that a ground vector can be differentiated.

[Proved where the ground state is nondegenerate] Spectral decomposition identifies the zero-temperature integrand as

    chi_xy=2 Re sum_(n>0)
       <0|V_x|n><n|V_y|0>/(E_n-E_0)
      =2 Re integral_0^infinity d tau <delta V_x(tau) delta V_y>_0.  (20)

This is an imaginary-time connected response between spatially separated reconnections. It is not an equal-time entanglement entropy or a real-time entangling rate.

[Proved sufficient extra hypothesis] Define the positive-part response budget

    R_L=limsup_(beta->infinity) integral_[0,1]^2 [chi_xy^(beta)]_+.

If, on the vacuum-compatible sequence,

    limsup_(L->infinity) L R_L <= B M log_2(d),                     (21)

then (19) and (1) give

    A <=(2B/pi) M log_2(d).                                      (22)

In particular B=pi^2/2 proves the requested sharp bound. Using the integral of the absolute value instead of the positive part is a stronger sufficient condition. No pointwise 1/L decay in s,t is required; such a pointwise assumption would unnecessarily constrain the possibly narrow crossover regions of the interpolation.

[Open] Condition (21) is the most targeted additional hypothesis obtained by this response route. It is not proved from locality, M, and d. The signed version of (21) is just a reformulation of the Casimir bound through (19); that would be circular as a proof. The positive-part version is stronger and identifies an independently estimable correlation function, but it is not claimed to be a weakest physical hypothesis among all possibilities.

[Proved obstruction to the elementary estimate] Kubo covariance Cauchy-Schwarz and the spectral logarithmic-mean bound give

    |chi_xy^(beta)|
       <=sqrt(chi_xx^(beta) chi_yy^(beta))
       <=beta ||V_x|| ||V_y|| <=16 beta M^2.                      (23)

For example chi_xx<=beta Var(V_x)<=beta ||V_x||^2 follows by bounding the logarithmic mean of the Gibbs weights by their arithmetic mean. This bound diverges in the ground-state limit. Finite d bounds the thermal entropy by 2L log d, but does not bound the inverse excitation energies in (20). Lieb-Robinson locality alone is a real-time statement and does not remove these inverse energies. A uniform gap could control them, at the cost of introducing a new scale and excluding a nonzero conformal Casimir coefficient in the usual gapped setting.

## 9. Solvable-chain checks and status

[Assumed continuum inputs, algebraic calibration] The critical XXZ family h=SxSx+SySy+Delta SzSz has M=1/2 for -1<=Delta<=1. For -1<Delta<1, writing Delta=cos gamma gives c=1 and v=pi sin gamma/(2gamma), so

    A/(pi M log_2 d)=sin gamma/gamma<=1.

The Heisenberg limit saturates. The endpoint Delta=-1 is not included as a relativistic CFT. Fundamental SU(N), h=P, has M=1 and ratio 2(N-1)/(N log_2 N). Symmetric SU(N) rank k has the C14/C35 ratio

    A/(pi M log_2 d)
       =2k(N^2-1)/[N(N+k)H_k log_2 binom(N+k-1,k)].                (24)

C35 already proves (24)<=1 for all N,k. In particular the SU(8), k=2 counterexample to the dimension-independent conjecture has ratio 0.406195448 here.

[Numerical] checks.py performs XXZ ED at Delta=-1/2,0,1/2,1 and L=4,6,8,10,12, open XX checks at the same sizes, and fundamental SU(N) ED at (N,L)=(2,4/6/8),(3,3/6/9),(4,4/8). It evaluates the inherited symmetric-family formulas, the exact partial-swap identity, and the telescoping obstruction sequence. The largest space is 4^8=65536; no large computation was run. All eigenpair residuals are below 6e-12. CHECKS.md records the numbers and the two expected failed finite-size inequalities.

[Numerical] response_check.py checks (18) for two three-site XX rings at beta=2 and 8 by a 24 by 24 Gauss-Legendre integration of the spectral Kubo formula. Odd lengths here test the identity, not the conformal vacuum coefficient. RESPONSE_CHECK.md records the residuals and signs.

[Conclusion, open universal question] No universal logarithmic Casimir bound is proved or falsified. The new rigorous statements are the explicit variational reductions, the failure mechanisms of the proposed elementary gluing constructions, and the exact two-seam representation (19) with the sufficient response condition (21). The main unresolved step is to bound that integrated response, or to construct an extension satisfying (10), without introducing an uncontrolled crossover length. Reflection positivity, endpoint dimension, and Schmidt rank alone do not complete this step.
