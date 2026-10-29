# C108. Charge curvature, stiffness, and propagation

## 1. Scope and conclusions

[Assumed] The microscopic system and gauge-optimized bond scale M are those of ../../specs/_bound_background.md. Choose a charge-invariant centered representative h with norm at most M+eta and take eta to zero at the end. The pairing note proves that this can be done. On a ring the gauge changes only a scalar energy, and on the infinite chain the telescoping term does not change the dynamics of a local observable. The lattice spacing and hbar equal one.

[Proved, lattice statements] For a nonconstant on-site charge q with endpoint spectral gaps delta_- and delta_+, the convex canonical energy has interior chemical potentials in

    -4M/delta_- <= mu <= 4M/delta_+ .                    (1)

Its distributional second derivative is a positive measure with total mass at most 4M(1/delta_- + 1/delta_+). In an integer normalization this is at most 8M. This includes any Mott cusp. It does not bound the ordinary second derivative pointwise.

[Proved, conditional theorem] Suppose a charged Gaussian mode appears with nonzero coefficient in the real-time scaling limit of a local conserved density. Then

    v_a <= v_* := 6 exp(1) M.                             (2)

For one non-drifting Luttinger liquid, with the stiffness and susceptibility matching assumed in the question, this gives

    D e'' = v^2 <= 36 exp(2) M^2.                        (3)

The coefficient is deliberately loose. The locality estimate is unconditional; its identification with an infrared velocity requires the stated local real-time matching. No uniform crossover length or lower bound on K is used.

[Proved, conditional subclass] If the m Gaussian modes have linearly independent microscopic charge-coupling rows, then m <= r <= d-1, where r is the dimension of the conserved on-site torus modulo the identity. Consequently

    A_G = sum_a v_a <= 6 exp(1) M (d-1).                 (4)

Nonzero coupling of every mode is weaker than row independence. Equation (4) has not been proved when many modes share the same microscopic charges. Independent neutral sectors are not covered by this counting argument.

## 2. Convex canonical energy

[Definition] Let q_- and q_+ be the smallest and largest distinct eigenvalues of q. The endpoint gaps are

    delta_+ = q_+ - max(spec(q) without q_+),
    delta_- = min(spec(q) without q_-) - q_- .

Degeneracies are allowed. The canonical energy e(n) is the thermodynamic energy per site at density n, including phase-separated states. Equivalently it is the closed convex energy obtained from the grand-canonical pressure. Exact finite-ring energies need not be convex in their integer charge label. We use their thermodynamic convex envelope, which is the physical canonical energy for a short-range chain.

[Proved] Joining long open chains at densities n_1 and n_2 changes energy by O(M) per join. A chain with fractions lambda and 1-lambda of the two phases therefore has energy density no greater than lambda e(n_1)+(1-lambda)e(n_2) at density lambda n_1+(1-lambda)n_2. Taking the thermodynamic limit proves convexity. The centered bond norm gives -M <= e(n) <= M. Every density in [q_-,q_+] can be approached by mixtures or macroscopic blocks of endpoint-charge states. The finite number of joins has vanishing energy per site.

[Proved, warning] Boundedness and convexity alone do not bound endpoint slopes. For example, -sqrt(n) is convex and bounded on [0,1] but has an unbounded negative slope at zero. Locality and charge conservation enter the next argument.

## 3. A local replacement channel and the saturation fields

[Proved] Let P project onto the q_+ eigenspace, let R=1-P, and choose a fixed unit vector |p> in P. For an orthonormal basis |b> of R, define the trace-preserving completely positive map

    T(rho) = P rho P + |p><p| tr(R rho),
    Kraus operators P and |p><b|.

The adjoint is unital. This map removes all nonmaximal charge states, retains the entire maximal-charge subspace, and does not divide by the norm of a charge-raising vector. That last feature avoids the small-denominator problem in a single-insertion variational argument.

[Proved] Charge conservation implies that h commutes with P tensor P. Indeed, this is the eigenspace of the two-site charge with maximal eigenvalue 2q_+. Let

    P_2 = P tensor P,
    B = (T^* tensor T^*)(h) - h.

The adjoint channel is block diagonal in the P,R labels at each site. Its restriction to P_2 agrees with h. Thus B=P_2^perp B P_2^perp. Unital positivity gives norm((T^* tensor T^*)(h)) <= M, hence norm(B) <= 2M. As an operator inequality,

    B <= 2M(1-P tensor P)
      <= 2M(R tensor 1 + 1 tensor R).                    (5)

The same argument works with M+eta before taking the infimum.

[Proved] Apply T on every site of a finite ring to any state rho. The output has maximal charge at every site. Let E_+(L) be the ground energy restricted to that subspace. Each site belongs to two bonds, so (5) yields

    E_+(L) <= tr(H T^tensor L(rho))
            <= tr(H rho) + 4M sum_i tr(rho R_i).

Because q_+-q >= delta_+ R,

    tr(H rho) >= E_+(L)
                - (4M/delta_+) [Lq_+ - tr(rho Q)].       (6)

No assumption about the state, its entanglement, or a nondegenerate endpoint eigenspace enters this inequality. Replacing toward the minimal eigenspace gives its lower-end counterpart.

[Proved] In the thermodynamic limit,

    e(n) >= e(q_+) - (4M/delta_+)(q_+-n),
    e(n) >= e(q_-) - (4M/delta_-)(n-q_-).                (7)

For a convex function and an interior subgradient mu at n,

    [e(n)-e(q_-)]/(n-q_-) <= mu
      <= [e(q_+)-e(n)]/(q_+-n).

Insert (7) to get (1). Equivalently, H-mu Q is saturated at the upper endpoint when mu>4M/delta_+, and at the lower endpoint when mu<-4M/delta_-. The external normal rays of the subdifferential at the endpoints are not part of the claimed bounded interval.

[Proved, normalization] For an integer-valued q, delta_+ and delta_- are at least one, even when some integer eigenvalues are missing. Thus the range length is at most 8M. For q -> a q+b with a>0, n and the endpoint gaps scale by a, while mu scales by 1/a. A universal formula proportional to the charge spread or the transfer bandwidth sigma_q does not have this covariance. For a diagonal bond sigma_q can vanish despite a nonzero chemical-potential range. The endpoint-gap formulation avoids both problems.

## 4. What can be averaged

[Proved] Write nu=d mu for the distributional second derivative of e on the open density interval. Convexity gives nu>=0, and (1) gives

    nu((q_-,q_+)) <= B_q,
    B_q := 4M(1/delta_- + 1/delta_+).                    (8)

The absolutely continuous part has density e'' almost everywhere. There may also be atoms at cusps and a singular continuous part. In particular,

    integral_I e''(n) dn <= B_q                          (9)

for any density interval I. If its length is ell, the mean curvature is at most B_q/ell. The integral, rather than an average over an arbitrarily short window, is bounded uniformly.

[Proved, conditional on the one-mode response dictionary at almost every density in a measurable set J] If J consists of non-drifting single-mode Luttinger states, D(n)<=sigma_q^2 M and v(n)^2=D(n)e''(n) give

    integral_J v(n)^2 dn <= sigma_q^2 M B_q,
    integral_J v(n) dn <= sqrt(|J| sigma_q^2 M B_q).      (10)

Here h and q are fixed and the canonical density varies. In a grand-canonical description the added field does not alter the double commutator controlling D. The bound uses M for the original bond h, not the larger norm obtained by incorporating an arbitrarily large field into that bond.

[Proved, obstruction to a purely convexity-based pointwise argument] On [0,1], set x=n-1/2 and

    f_epsilon(n) = M sqrt(x^2+epsilon^2) - c_epsilon,

where c_epsilon centers its range. Then |f_epsilon|<=M/4, |f_epsilon'|<=M, and

    f_epsilon'' = M epsilon^2/(x^2+epsilon^2)^(3/2),
    f_epsilon''(1/2)=M/epsilon.

Its integrated curvature is less than 2M for every epsilon>0. Thus both energy and slope bounds allow unbounded pointwise curvature. Taking a constant D_0=M also respects a scalar stiffness budget with sigma_q=1, but D_0 f_epsilon''(1/2)=M^2/epsilon. The joint local response ansatz f_epsilon(n)+D_0(1-cos a) has a positive diagonal quadratic response at a=0 and does not remove this freedom. These are response functions, not a claimed lattice realization. They show why the separate static bounds or Hessian positivity cannot prove (3).

## 5. Mott edges and actual singularities

[Assumed, non-drifting Luttinger normalization] Use

    H_LL = (v/(2pi)) integral [K theta_x^2 + K^(-1) phi_x^2] dx,
    rho = -phi_x/pi.

A uniform density shift costs pi v (delta n)^2/(2K), and a uniform twist costs vK a^2/(2pi). Hence

    e''=pi v/K,  D=vK/pi,  kappa=K/(pi v),  D e''=v^2.   (11)

These identities concern a smooth compressible branch and the corresponding thermodynamic responses. They do not multiply an atom in nu by a vanishing stiffness.

[Proved from the stated Mott-edge behavior] If v tends to zero and K tends to a finite positive K_c, then e'' tends to zero, kappa diverges, and D tends to zero. A generic dilute soliton gas with charge s, effective mass m_s, and soliton density rho_s=|delta n|/|s| illustrates the mechanism:

    e(n)-e(n_c)-mu_c delta n
       = pi^2 |delta n|^3/(6m_s |s|^3) + o(|delta n|^3),
    e'' = pi^2 |delta n|/(m_s |s|^3) + o(|delta n|),
    v = pi |delta n|/(m_s |s|) + o(|delta n|),
    K -> s^2,
    D = |s| |delta n|/m_s + o(|delta n|).

The cubic energy follows by integrating k^2/(2m_s) over the filled interval [-pi rho_s,pi rho_s]. Twice differentiating and inserting in (11) checks every coefficient. The dilute-soliton description is an infrared input, not a theorem for every plateau edge.

[Assumed, literature input] The conventional commensurate-incommensurate Mott transition has this vanishing velocity, divergent compressibility, and finite limiting Luttinger parameter. See Cazalilla et al., Section VI.A.3, especially Eq. (160) and the following paragraph in the linked version, [One dimensional Bosons](https://arxiv.org/html/1101.5337v1#S6.SS1.SSS3). A plateau corresponds to a jump of mu as a function of density, hence an atom of nu at its commensurate density.

[Proved, exact lattice cusp] Take d=3, q=diag(0,1,2), and

    h = (U/2)[(q-1)^2 tensor 1 + 1 tensor (q-1)^2], U>0.

The uncentered bond has eigenvalues 0,U/2,U and M=U/2. No gauge lowers this: its expectation on uniform |0> and |1> product bonds differs by U, and a telescoping gauge term has zero expectation on either. After centering,

    e(n)=U |n-1|-U/2,  0<=n<=2,
    nu=2U delta_1,  D=0.

The cusp has mass 2U=4M and obeys (8). There is no charged Luttinger sector here and no ordinary finite e'' at n=1. This example prevents interpreting zero compressibility inside a plateau as a divergent Luttinger velocity.

[Proved, status] No fixed-d nearest-neighbor family with unbounded smooth e''/M for a fixed integer charge is constructed here. No pointwise upper bound on that smooth curvature alone is proved. The abstract convex example does not settle this microscopic question. A pointwise bound on the product is available independently, as follows.

## 6. A dimension-independent lattice speed

[Proved] For one-site observables A_x,B_y at separation r>=1, the usual commutator recursion gives

    C_B(X,t) <= C_B(X,0)
       + 2 sum_{b: b intersects X} ||h_b||
             integral_0^|t| C_B(b,s) ds,
    C_B(X,t)=sup_{A supported in X} ||[tau_t(A),B]||/||A||.

One obtains this by differentiating the commutator, separating the homogeneous commutator evolution, which preserves its norm, and integrating the remaining source term. This recursion and its iterated interaction paths are Eqs. (2.9)-(2.12) in [Nachtergaele, Ogata, Sims](https://arxiv.org/pdf/math-ph/0603064). The following explicit constant is the elementary one-dimensional path count.

[Proved] The first bond in a path has at most two choices. Each later bond intersects the preceding bond and has at most three choices, including repetition. Reaching site y from site x requires at least r bonds. With all bond norms at most M, iteration therefore yields

    ||[tau_t(A_x),B_y]||
       <= (4/3)||A_x||||B_y|| sum_{k>=r}(6M|t|)^k/k!.

For any a>0, 1_{k>=r}<=exp(a(k-r)), so

    ||[tau_t(A_x),B_y]||
       <= (4/3)||A_x||||B_y||
          exp[-a r + 6M exp(a)|t|].                      (12)

The velocity in (12) is 6M exp(a)/a. Its minimum is 6 exp(1)M at a=1. Take the infinite-chain limit before any large-time scaling. The estimate is state independent. One can instead work in rings whose circumference grows faster than the distances in question.

[Proved] The use of the optimized representative is legitimate. For a local A the sum of commutators with o_i-o_{i+1} telescopes to zero, and scalar shifts do not affect dynamics. Thus the bound depends on M, not a gauge-dependent original bond norm. If a density is selected by H-mu Q, then q_x commutes with Q and its Heisenberg evolution is the same under H and H-mu Q. Large chemical potentials do not enlarge the density propagation cone.

## 7. Passing the cone to the Gaussian limit

[Assumed, explicit real-time matching] For each linear combination q^u of the conserved local densities, the rescaled expectation of its commutator has a distributional scaling limit. Its smooth long-wavelength part is the Gaussian current commutator with the same charge-coupling matrix Z that enters stiffness and susceptibility. This assumption is stronger than fitting a finite-size Casimir term or postulating static response coefficients alone. Oscillating lattice components are removed by smooth spatial smearing. Only convergence for test functions supported away from (x,t)=(0,0) is needed.

[Proved from (12)] Let a test function have compact support in the strict region |x|>v_*|t|, separated from its boundary by a positive margin. At lattice coordinates (lambda x,lambda t), (12) bounds every sampled commutator by a constant times exp(-c lambda), uniformly on that support. Smooth integration, the spatial Riemann sum, and any fixed power of lambda from field normalization add only polynomial factors. Their limit is zero. Therefore any distributional scaling limit of these local commutators is supported in |x|<=v_*|t|. No estimate on how soon the scaling limit becomes accurate is required.

[Proved within the Gaussian theory] A canonically normalized boson density j=-phi_x/pi has chiral currents with opposite current commutators. With the convention [phi(x),theta_y(y)]=i pi delta(x-y), this gives

    [j(x,t),j(0,0)]
       = (i/(2pi))[delta'(x+vt)-delta'(x-vt)] .           (13)

For q^u=sum_a (Z_a u) j_a, the right side is summed over a with nonnegative weights (Z_a u)^2. Distinct positive velocities give distinct rays for t>0. On a fixed ray equal velocities add with the same sign. If Z_a is nonzero, choose u with Z_a u nonzero. A test function near x=v_a t, t>0, detects the delta derivative. Hence v_a>v_* would contradict the support conclusion above. This proves (2), and (11) then proves (3).

[Proved, consequence] At any smooth one-mode point with e''>0,

    D <= 36 exp(2) M^2/e''.                              (14)

Thus a family with e''/M tending to infinity would necessarily have D/M tending to zero, provided its local real-time Luttinger identification holds at each member. The statement is pointwise in the Hamiltonian and does not interchange the family limit with the infrared limit.

[Proved, scope] This is a propagation proof of the product inequality, not a new joint twist-density variational inequality. The static hypotheses in Section 4 do not imply it. It is not an unconditional theorem that a Casimir fit alone implies the real-time matching in this section.

## 8. Several charges and the correct matrix product

[Assumed, Gaussian response matching without drift] Let there be m modes, r independent microscopic charges, V=diag(v_a)>0, and an m by r real charge matrix Z. Then

    D = Z^T V Z/pi,
    kappa = Z^T V^(-1) Z/pi.                             (15)

The inverse-Z expression printed in the original pairing README is inconsistent with its scalar normalization. Equation (15) is the correction already recorded in C104. For one mode Z=sqrt(K), it reproduces (11).

[Proved] If Z is square and invertible, E=Hess_n e=kappa^(-1) and

    D E = Z^T V^2 Z^(-T).

Although D E need not be symmetric, it is similar to the positive symmetric matrix E^(1/2) D E^(1/2). Its eigenvalues are v_a^2, and

    sum_a v_a = tr sqrt(E^(1/2) D E^(1/2)).              (16)

For full row rank m<=r, restrict the charges to the m-dimensional active response subspace. The same calculation applies there. Any gapped charge combinations with singular static response must first be excluded. The inverse Hessian here refers to the active compressible variables, with other densities allowed to relax consistently with (15).

[Proved, conditional on Section 7] Each nonzero row has v_a<=v_*. Consequently, whether or not Z is square,

    D <= v_*^2 kappa                                    (17)

as a positive-matrix inequality. This follows term by term from v_a <= v_*^2/v_a. It gives an upper bound v_*^2 on the generalized response eigenvalues. In the full-row-rank case these are all m physical squared velocities, so (16) proves A_G<=m v_*.

[Proved] The number r of independent commuting Hermitian on-site charges modulo the identity is at most d-1, since they can be diagonalized simultaneously. If the rows of Z are independent, m<=r and (4) follows. This is a crossover-free linear bound in d for that subclass. It requires no lower singular-value bound on Z. The torus stiffness budget tr D<=4M remains valid, but is not needed once locality and the full-row-rank hypothesis are supplied.

[Proved, obstruction to the proposed general sector count] If m>r, every row of Z can be nonzero. For any charge vector u with nonzero susceptibility,

    (u^T D u)/(u^T kappa u)
       = [sum_a v_a (Z_a u)^2]
         /[sum_a (Z_a u)^2/v_a]

is a weighted mean of v_a^2, with weights (Z_a u)^2/v_a. It does not enumerate the modes. For example, with one charge, m modes all at speed v_0, and Z_a=1/sqrt(m),

    D=v_0/pi,  kappa=1/(pi v_0),  D/kappa=v_0^2,
    A_G=m v_0.

The response matrices are independent of m. Every mode is charged. A compact-boson realization of these continuum weights takes the diagonal integer charge of m bosons and K_a=1/m; integer vertex charges and locality remain intact. This is not a fixed-d lattice construction. It proves that response ranks, compactification, and the two budgets alone cannot justify m<=r. Locality bounds the speeds, but not the number of nonzero rows.

[Proved, remaining target] To extend (4) to all charged Gaussian infrared theories, one needs an additional microscopic aggregate constraint on the multiplicity and weights, or a proof that a specified subclass has full row rank. No such general constraint is proved here. The supplied fixed-d slow-mode constructions make an unconditional count especially inappropriate. This run does not claim a counterexample to A<=C M f(d).

## 9. Unequal left and right speeds

[Proved, qualification] The scalar identity in the question assumes the non-drifting Luttinger Hamiltonian in Section 5. In a more general Gaussian Hamiltonian add 2u theta_x phi_x inside that integrand with v>|u|. The two propagation speeds are v_R=v+u and v_L=v-u. At fixed charge the flux curvature remains D=vK/pi. At a nonzero density shift, thermodynamic minimization over the uniform theta gradient instead gives

    e'' = pi(v^2-u^2)/(vK),
    D e'' = v^2-u^2 = v_R v_L.

To see this, complete the square vK(theta_x+u phi_x/(vK))^2+(v^2-u^2)phi_x^2/(vK). If winding or current is held fixed instead, the response being measured is different. Thus the arithmetic-mean velocity squared is not generally D e'' after current relaxation. The lattice cone bounds both v_R and v_L, so their mean, their product, and the conditional bound (4) are still bounded by the corresponding powers of v_*. The simple matrix identities (15)-(16) are restricted to the non-drifting case.

## 10. Checks and open statements

[Proved, exact free check] For a spinless nearest-neighbor hopping chain with dispersion -2t cos k and charge n in [0,1], M=t and

    e(n)=-(2t/pi)sin(pi n),
    mu=-2t cos(pi n),  e''=2pi t sin(pi n),
    v=2t sin(pi n),  D=2t sin(pi n)/pi.

Thus D e''=v^2, the chemical-potential range is 4t<=8M, and the integrated curvature is 4t. At the dilute endpoints e'' and D both vanish. The velocity bound is loose but valid. M=t follows from the charge-zero hopping block with eigenvalues +/-t and its equal one-site marginals, which make its extremal expectation difference invariant under telescoping gauges.

[Proved, exact XXZ check using the inherited formulas] At zero magnetization, Delta=cos gamma, J=1, and 0<gamma<pi,

    M=1/2,
    v=pi sin(gamma)/(2gamma),
    K=pi/[2(pi-gamma)],
    e''=pi(pi-gamma)sin(gamma)/gamma,
    D=pi sin(gamma)/[4gamma(pi-gamma)].

As gamma tends to zero these give v=pi/2, K=1/2, e''=pi^2, D=1/4. At gamma=pi/2 they give the free values v=1, K=1, e''=pi, D=1/pi. Writing epsilon=pi-gamma near the ferromagnetic point, v is asymptotic to epsilon/2, e'' to epsilon^2, and D tends to 1/4. The curvature vanishes there even though the bond derivative of v diverges, as found in the earlier campaign. None of these checks supplies a universal curvature bound.

[Numerical, validation scope] check.py checks the local replacement-channel operator inequality on small charge-invariant matrices, the Gaussian matrix eigenvalue identity, the multiple-mode response example, and the elementary scalar formulas. These are small finite-dimensional checks of the derivation. No many-body scaling fit, heavy numerical calculation, or cluster job is needed.

[Proved, final status] The chemical-potential and integrated-curvature bounds are lattice results. The pointwise product bound follows from locality under an explicit real-time Gaussian matching assumption. The linear-in-d bound follows when the Gaussian modes are independently resolved by microscopic charges. A general smooth-curvature bound and a bound on the sum over arbitrarily many modes sharing those charges remain open.
