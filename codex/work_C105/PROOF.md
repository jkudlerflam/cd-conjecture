# C105 derivations

## 1. Scope and conventions

[Assumed] The uniform chain has the infrared sectors and Casimir expansion stated in specs/_bound_background.md. Write A = sum_a c_a v_a. The symbol C in the requested inequality denotes a universal constant, not the sum of central charges. The interaction scale M is the gauge-optimized half-spread.

[Assumed] For the smooth energy density in the appropriate gauge, the infrared identification is

    delta_smooth(x) = A/(24 pi x^2) + o(x^(-2)),
    Delta_N = e0 W_N - E0(K_w) = A/(12 pi) log N + O(1),
    w_x = x(N-x)/N,  1 <= x <= N-1,
    W_N = (N^2-1)/6.

The second identity is the integrated version used below. The raw lattice bond observable requires a separate check. Section 2 supplies that check for XX and shows why it cannot be identified with the smooth stress-energy density without a prescription for oscillating terms.

[Proved] Shifting the bond by a scalar changes both e0 W and E0(K_w) by the same amount. Delta_N is invariant under this shift. A discrete divergence is different. For h_i -> h_i + o_i - o_(i+1),

    K_w -> K_w + w_0 o_0
                     + sum_(i=1)^(N-2) (w_i-w_(i-1)) o_i
                     - w_(N-2) o_(N-1).

It leaves the uniform periodic Hamiltonian unchanged but generally changes the weighted ground state. For a linear ramp it supplies a uniform interior field. A bound in terms of the gauge-invariant M therefore needs a specified representative of h and an infrared matching statement in that representative.

## 2. The exact half-line XX spectral formula

[Proved] Use sites n = 0,1,... and positive bond coordinate x = n+1. The bond operator is

    h_x = (c_(x-1)^dagger c_x + c_x^dagger c_(x-1))/2,
    e0 = -1/pi,  M = 1/2,  A = 1.

The zero-weight site in the convention sum_(x>=0) x h_(x,x+1) decouples. Removing it gives these coordinates without changing the large-distance question. The one-particle Rindler operator J has J_(n,n+1) = (n+1)/2 and zero diagonal. Carleman's criterion applies because sum_n 2/(n+1) diverges, so the Jacobi operator has a unique self-adjoint closure.

[Proved, using the cited orthogonal-polynomial identities] Define

    P_(-1)(lambda) = 0,  P_0(lambda) = 1,
    (n+1) P_(n+1)(lambda) = 2 lambda P_n(lambda) - n P_(n-1)(lambda).

These are the Meixner-Pollaczek polynomials with parameters 1/2 and pi/2. Their orthogonality in this normalization is

    integral_R P_m(lambda) P_n(lambda) sech(pi lambda) dlambda = delta_mn.

Consequently the spectral transform sends J to multiplication by lambda. The Rindler ground-state covariance is the negative spectral projector 1_(J<0). The uniform infinite-chain covariance restricted to the half-line is

    C_U = (1 + exp(2 pi J))^(-1).

The latter exact half-line relation is proved in Eisler, J. Stat. Mech. (2025) 013101, Section 3 and Appendix B, with the sign of hopping changed by a staggered unitary. The recurrence, orthogonality and generating function are also in DLMF 18.19, 18.22 and 18.23. Sources are linked in Section 9.

[Proved] Since the bond energy is the real adjacent covariance, subtracting the two states gives the convergent formula

    delta_x := e0 - <h_x>_B
       = integral_R g(lambda) P_(x-1)(lambda) P_x(lambda)
                        sech(pi lambda) dlambda,
    g(lambda) = 1/(1+exp(2 pi lambda)) - 1_(lambda<0)
              = sign(lambda)/(1+exp(2 pi |lambda|)).

P_n has parity n. Therefore

    delta_x = 2 integral_0^infinity
                P_(x-1)(lambda) P_x(lambda)
                sech(pi lambda)/(1+exp(2 pi lambda)) dlambda.       (2.1)

This formula involves no finite-chain cutoff. It is the quantity evaluated in xx_halfline.py.

[Proved] Truncating the full-line integral to |lambda| <= T has an absolute error at most exp(-2 pi T). Indeed, |g| <= exp(-2 pi T) on the omitted set and Cauchy-Schwarz with the orthogonality measure gives integral |P_n P_(n+1)| sech(pi lambda) <= 1. At T = 12 the bound is 1.80e-33. This controls the spectral tail, not the quadrature error inside the retained interval.

[Proved] The first bond has a closed form. With G the Catalan constant,

    <h_1>_B = -2 integral_0^infinity lambda sech(pi lambda) dlambda
             = -4G/pi^2,
    delta_1 = 4G/pi^2 - 1/pi = 0.05291698652698... .

The second equality follows by expanding sech(pi lambda) = 2 sum_(k>=0) (-1)^k exp(-(2k+1)pi lambda) and integrating lambda. In particular, the proposed sharp all-bond bound delta_1 <= M/24 = 1/48 fails. This finite-distance observation does not by itself refute an asymptotic bound with an arbitrary constant.

### 2.1 What the large-distance expansion does and does not establish

[Proved, fixed spectral parameter] The generating function is

    sum_(n>=0) P_n(lambda) z^n
      = (1-iz)^(-1/2+i lambda) (1+iz)^(-1/2-i lambda).

Expanding its two algebraic singularities at z = i and z = -i gives, uniformly for lambda in a fixed compact real interval,

    P_n(lambda) = (2n)^(-1/2) [ i^(-n) (2n)^(i lambda)/Gamma(1/2+i lambda)
                              + i^n (2n)^(-i lambda)/Gamma(1/2-i lambda) ]
                  + O_lambda(n^(-3/2)).                          (2.2)

For example, the singularity at z = i contributes the coefficient of

    2^(-1/2+i lambda) (1-z/i)^(-1/2-i lambda).

Using [z^n](1-z/i)^(-alpha) = i^(-n) Gamma(n+alpha)/(Gamma(alpha)Gamma(n+1)) and the gamma-ratio expansion gives its term in (2.2). The conjugate singularity supplies the other term. Taylor expansion of the analytic factor at each singularity gives the indicated next power of n.

[Heuristic, asymptotic integration not completed here] Let L_n = log(2n). Multiplying the leading terms for n and n+1, the oscillatory contribution to (2.1) is

    (-1)^n/n integral_R g(lambda) sech(pi lambda)
        Im[ exp(2i lambda L_n)/Gamma(1/2+i lambda)^2 ] dlambda.

For lambda > 0 set

    F(lambda) = sech(pi lambda)/[(1+exp(2 pi lambda))Gamma(1/2+i lambda)^2].

F(0) = 1/(2 pi). The full integral is twice the imaginary part of integral_0^infinity F(lambda) exp(2i lambda L_n) dlambda. Integration by parts supplies the leading value 1/(2 pi L_n). It predicts

    delta_x = (-1)^(x-1)/(2 pi x log x) [1 + O(1/log x)]
                 + a smooth contribution of order x^(-2).       (2.3)

The derivative of the gamma phase replaces log(2n) by log(2n)-psi(1/2) = log(8n)+EulerGamma at the next level of approximation.

[Heuristic, same qualification] The nonoscillatory product in (2.2) comes from the phase difference pi/2 - lambda/n. Its contribution is lambda/[n^2 |Gamma(1/2+i lambda)|^2]. Since |Gamma(1/2+i lambda)|^2 = pi sech(pi lambda), the resulting smooth coefficient is

    1/pi integral_R lambda g(lambda) dlambda
      = 2/pi integral_0^infinity lambda/(1+exp(2 pi lambda)) dlambda
      = 1/(24 pi).

[Assumed only for upgrading (2.3) to a theorem] An integrable remainder estimate for (2.2) uniform on spectral windows that grow with n is still needed. A fixed compact-window Darboux expansion cannot on its own justify inserting (2.2) into the full spectral integral. No theorem asserting divergence of the raw scaled profile is claimed in this report. The exact formula (2.1), the first-bond counterexample and the numerical results do not depend on that unresolved remainder estimate.

[Numerical] The exact half-line quadrature gives x^2 delta_x = 0.84788490, 4.95491062, 16.93641809, 111.42982189 at odd x = 33,257,1025,8191. Define the centered filter

    delta_bar_x = (delta_(x-1) + 2 delta_x + delta_(x+1))/4.

At the same coordinates x^2 delta_bar_x is 0.01280006, 0.01321859, 0.01325355, 0.01326197. The smooth prediction is 0.0132629119243. Doubling the quadrature order from 400 to 800 changes any computed delta_x through x = 8193 by at most 2.08e-15. The computation took 0.12 seconds and allocated only a few one-dimensional arrays.

[Proved] A fixed kernel can annihilate a constant-amplitude alternating component without annihilating the smooth component. For a smooth envelope b(x), the filter applied to (-1)^x b(x) is (-1)^x[2b(x)-b(x-1)-b(x+1)]/4. If b(x) is of order 1/(x log x) with the corresponding derivative bounds, its filtered contribution is of order 1/(x^3 log x). This explains why the raw and filtered profiles test different assertions. The analogous SU(3) filter used numerically is (1,2,3,2,1)/9, with a double zero at momenta plus or minus 2pi/3.

## 3. A logarithmic lower-energy certificate for XX

[Proved] There is a nonlocal positive-operator certificate that avoids the layer boundary loss for XX. It works for arbitrary finitely supported nonnegative site amplitudes a_i, with bond weights w_i = a_i a_(i+1). Extend a by zero outside its support. Then

    E0(K_w) >= -(1/pi) sum_i a_i^2,
    Delta(w) <= (1/(2 pi)) sum_i (a_(i+1)-a_i)^2.                 (3.1)

The sums in the second expression include both exterior zero-to-nonzero steps.

[Proved] To see this, let T be the uniform one-particle hopping operator on the full line, with T_(i,i+1) = 1/2. Its Fourier symbol is cos k, so

    -|T| <= T <= |T|,
    <i||T||i> = integral_(-pi)^pi |cos k| dk/(2 pi) = 2/pi.

Let D_a be diagonal multiplication by a, and restrict to the finite support. The weighted hopping matrix is J_w = D_a T D_a. Congruence preserves positivity, hence

    -Q <= J_w <= Q,  Q = D_a |T| D_a >= 0,
    Tr Q = (2/pi) sum_i a_i^2.

If u_j is an eigenbasis of J_w with eigenvalues lambda_j, its two operator bounds give |lambda_j| <= <u_j,Q u_j>. Summing gives Tr|J_w| <= Tr Q. Since Tr J_w = 0, filling all negative one-particle energies gives

    E0(K_w) = -(1/2) Tr|J_w| >= -(1/pi) sum_i a_i^2.

Subtract this from e0 W = -(1/pi) sum_i a_i a_(i+1). Finally,

    sum_i a_i^2 - sum_i a_i a_(i+1)
       = (1/2) sum_i (a_(i+1)-a_i)^2,

which proves (3.1). This is a lower bound on the ground energy, with the required variational direction.

[Proved] Quadratic interval weights admit an explicit, slowly varying factorization. For n >= 0 set

    r_n = 2 [Gamma((n+2)/2)/Gamma((n+1)/2)]^2,
    a_i^2 = r_i r_(N-1-i)/N,  0 <= i <= N-1.

The gamma recurrence gives r_n r_(n+1) = (n+1)^2. Therefore

    a_i a_(i+1) = (i+1)(N-1-i)/N,

which is the desired weight. Stirling expansion gives

    r_n = n + 1/2 + 1/[8(n+1/2)] + O((n+1)^(-3)).

Writing b_i = i+1/2, the base term obeys

    (1/N) sum_i b_i b_(N-1-i) - W_N = 1/4.

The two first corrections give

    (1/(8N)) sum_i [b_(N-1-i)/b_i + b_i/b_(N-1-i)]
      = (1/4) sum_i 1/(i+1/2) - 1/4
      = (1/4) log N + O(1).

The summed remainders are O(1), including the finitely many end sites. Thus

    sum_i a_i^2 - W_N = (1/4) log N + O(1),
    Delta_N <= (1/(4 pi)) log N + O(1).                         (3.2)

[Proved, conditional on the infrared matching] Equation (3.2) gives A <= 3 = 6M for unit XX. This is weaker than the already proved quasi-free theorem A <= 2M in the background. Its use here is methodological. A lower-energy certificate with logarithmic loss exists, and its positive comparison Q is nonlocal. It does not arise by adding independent uniform open-chain lower bounds.

[Numerical] Equation (3.1) passes direct Gaussian diagonalization for N = 8,16,32,64,128,256. At N = 256 the true defect is 0.2077077035 and the certificate allowance is 0.5131771861. The factorization residual is below 5.9e-12. These are small checks, not the unrun large finite-chain job.

[Proved, candidate rejected] A naive interacting extension would replace -1/pi in the first inequality of (3.1) by the centered uniform energy density e0. Equivalently it would assert Delta <= -e0 sum_i (a_(i+1)-a_i)^2/2. It fails already for four Heisenberg sites. Take h = S.S + 1/4 = P/2, with e0 = 1/2-log(2), a = (3/4,1,1,3/4), and w = (3/4,1,3/4). Here sum_i a_i^2 = 25/8. On the two-dimensional singlet subspace, in the basis formed by the product of the two outer singlets and the two outer triplets coupled to total spin zero, the centered weighted Hamiltonian is

    K_singlet = [ -1/2     sqrt(3)/4 ]
                [ sqrt(3)/4     1/2 ].

Its lower eigenvalue is -sqrt(7)/4. Thus a singlet trial vector already gives

    E0(K_w) <= -sqrt(7)/4 = -0.6614378278...
             < (25/8)(1/2-log(2)) = -0.6035849392... .

This contradicts the proposed lower-energy bound without needing to prove that the singlet is the full ground state. Equivalently the trial-state deficit is 0.1785698764, while the proposed allowance is only (5/8)(log(2)-1/2) = 0.1207169878.

[Numerical] The failure persists at N = 12, where the Heisenberg defect is 0.2213542002 and the allowance is 0.1669331955. The same small test passes for XX, XXZ at anisotropy 1/2, and SU(3), but those passes do not establish the extension. The scalar centering matters for the proposed right side, although Delta itself is invariant.

## 4. Positive uniform-interval decompositions cannot telescope

[Proved] Suppose a nonnegative bond weight is represented as

    w_i = sum_I lambda_I 1_(i in I),  lambda_I >= 0,

where I ranges over arbitrary intervals of bonds, with no nesting restriction. Put g(n) = n e0 - E0(H_open on n bonds). Superadditivity of the minimum eigenvalue gives

    Delta(w) <= B(lambda) := sum_I lambda_I g(|I|).              (4.1)

At a bond i_* where w is largest,

    max_i w_i = sum_(I containing i_*) lambda_I <= sum_I lambda_I.

Thus if g(n) >= g_* > 0 for every allowed n,

    B(lambda) >= g_* max_i w_i                                (4.2)

for every positive interval decomposition, including the optimizer. Allowing overlaps can change constants, but cannot remove this linear loss.

[Proved] For XX this has an explicit uniform constant. Let q = n+2 and z = pi/(2q). Summing the filled open-chain single-particle energies cos(k pi/q) gives

    E0_open(n) = 1/2 - (1/2)csc z,  n odd,
    E0_open(n) = 1/2 - (1/2)cot z,  n even.

Consequently, with b = 2/pi - 1/2,

    g(n) = b + [csc z - 1/z]/2,  n odd,
    g(n) = b + [cot z - 1/z]/2,  n even.

For odd n, sin z <= z gives g(n) >= b. For even n >= 2, z <= pi/8. The derivative of cot z - 1/z is -csc(z)^2 + 1/z^2 < 0, so its minimum on this range is at z = pi/8. Therefore

    g(n) >= g_* = 1/sqrt(2) - 2/pi = 0.07048700881897... > 0.

Equality is attained at n = 2. For the quadratic weight max w is N/4 for even N, so even the optimal positive interval allowance is at least g_* N/4. This cannot prove an O(log N) allowance. Signed decompositions, boundary operators and nonlocal positive certificates are outside this obstruction.

[Numerical] Linear programming over every interval, with exact-diagonalization block energies, satisfies Delta <= B_opt <= B_layer for every tested chain and weight. The N = 12 quadratic XX values are 0.1312278503 <= 0.3371810216 <= 0.4547805842. The improvement over nested layers is real, but it does not change the obstruction (4.2).

## 5. Fixed-range local-marginal SDP certificates also have a bulk loss

[Proved] Consider the usual r-site marginal relaxation with positivity, unit trace, and equality of the left and right (r-1)-site marginals. Assume the uniform ground state has a strictly positive r-site marginal rho_r and that its energy e0 is below the maximally mixed bond energy e_mix. The maximally mixed r-site matrix u_r = I/d^r satisfies the same affine constraints. Let m_r be the least eigenvalue of rho_r. For sufficiently small t > 0,

    sigma_r = (1+t)rho_r - t u_r

is positive and satisfies all the constraints. A concrete choice is t = d^r m_r/2. Indeed its smallest eigenvalue is

    (1+t)m_r - t/d^r = m_r/2 + t m_r > 0.

Its common bond energy is

    e_sigma = (1+t)e0 - t e_mix = e0 - t(e_mix-e0) < e0.         (5.1)

Thus this finite-range relaxation is strictly below the true uniform energy for every such r.

[Proved] Use sigma_r on every r-site window of a weighted interval, and use its reductions on shorter windows near the ends. All overlaps match. This feasible point of the local-marginal relaxation has weighted energy e_sigma W. By weak duality, no dual certificate from this relaxation can have a lower-energy bound exceeding e_sigma W. An attempted bound E0(K_w) >= e0 W - B therefore requires

    B >= (e0-e_sigma) W = t(e_mix-e0) W.                        (5.2)

For fixed r and the quadratic weight, W grows as N^2. Hence no fixed-r certificate using only these local consistency and positivity constraints can prove B = O(log N). This includes position-dependent dual multipliers. It does not exclude r increasing with N, nonlocal constraints, stationarity constraints for the weighted ground state, or special fermionic covariance certificates such as Section 3.

[Proved] The hypothesis holds for every finite r in XX. The restricted Fermi covariance has quadratic form

    z^dagger C_r z = integral_F |sum_(j=0)^(r-1) z_j exp(ikj)|^2 dk/(2 pi).

A nonzero finite trigonometric polynomial cannot vanish on a momentum interval. Both the Fermi interval and its complement have positive length. Therefore 0 < C_r < I. All eigenvalues nu_j lie strictly between zero and one, and the Gaussian density matrix has least eigenvalue product_j min(nu_j,1-nu_j) > 0. Here e_mix = 0 and e0 = -1/pi, so (5.1) is a strict improvement for every finite r.

[Numerical] The explicit feasible improvement e0-e_sigma is 0.02101564467 for r = 2, 0.001581496699 for r = 3, and 1.635184903e-5 for r = 4. At r = 8 it is about 7.45e-19 and the energy difference is lost in double-precision subtraction. Strict positivity for every finite r is established analytically above, not inferred from that last floating-point energy.

## 6. What the local deformation route still needs

[Proved] A trial vector formed from the uniform ground state gives E0(K) <= <K>_trial. It supplies a lower bound on the deficit, not the upper bound requested here. A useful variational argument must provide a lower bound valid for all vectors, or a controlled inverse transformation taking the unknown weighted ground state back to a uniform-chain comparison state.

[Proved, finite system with a simple ground eigenvalue along the path] For H(s) = H0 + sV, let E(s) be the ground energy, |0_s> the ground vector and |m_s> excited eigenvectors. Differentiation of the eigenvalue equation, with <0_s|partial_s 0_s> = 0, gives

    partial_s |0_s> = -sum_(m>0) |m_s> <m_s|V|0_s>/(E_m(s)-E(s)),
    E''(s) = -2 sum_(m>0) |<m_s|V|0_s>|^2/(E_m(s)-E(s)).

Integrating twice gives

    E(0)+E'(0)-E(1)
      = 2 integral_0^1 (1-s) sum_(m>0)
          |<m_s|V|0_s>|^2/(E_m(s)-E(s)) ds.                    (6.1)

Reflection symmetry can force E'(0) = <V>_0 = 0 when V is odd. It does not bound the inverse-energy denominators in (6.1). With a uniform finite-volume gap gamma, the expression is at most ||V||^2/gamma. No such gap is among the critical-chain assumptions.

[Heuristic, unresolved route] A Lieb-Robinson estimate controls real-time spatial spreading. A localized adiabatic continuation also needs control of its low-frequency filter, which is the missing inverse-energy response in (6.1). I have not obtained a bound on that response by M f(d), nor an operator identity that cancels its dangerous part. This is an obstruction to the attempted estimate, not a theorem that every localized-deformation proof fails.

## 7. Numerical checks and their limits

[Numerical] checks.py uses deterministic sparse diagonalization in balanced charge sectors, N <= 12 for XX, Heisenberg and XXZ(1/2), and N <= 9 for SU(3). These are small checks with at most 1680 basis vectors. The full suite, including interval linear programs, took 0.18 seconds. Its largest eigenvector residual was 6.67e-12 and the local-energy sum agrees with the eigenvalue to below 3.6e-14. Sector restriction is stated explicitly. These numerical eigenvalues are not certified lower bounds.

[Numerical] The last two quadratic sizes give defect slopes 0.02396377 for XX, 0.03964901 for Heisenberg, 0.03189360 for XXZ(1/2), and 0.10340385 for SU(3). Their infrared predictions are 0.02652582, 0.04166667, 0.03445806, and 0.11111111. This is consistent with the integrated matching but does not prove the asymptotic statement. At these sizes the raw profiles oscillate strongly, and even filtered finite-interval values can exceed the asymptotic coefficient.

[Proved, infrastructure observation] Both short SSH attempts to Typhon failed before authentication with "Connection timed out during banner exchange". No remote job was submitted. The large finite-chain Gaussian script and a SLURM job file are prepared, but not run. The script rejects --cluster unless SLURM_JOB_ID is set. No heavy calculation was moved to this machine or a login node. The exact half-line calculation reaches bond 8193 using a 0.12-second scalar recurrence rather than a large finite matrix.

## 8. What is established

[Proved] The new lower-energy certificate (3.1) has logarithmic loss for XX quadratic weights. Positive sums of arbitrary uniform intervals cannot have logarithmic allowance even for XX. Fixed-range local-marginal SDP duals cannot have logarithmic allowance for XX. A naive interacting extension of the XX certificate is false by an explicit four-site Heisenberg witness.

[Heuristic] The raw Rindler bond deficit has a staggered 1/(x log x) contribution in XX. The exact spectral formula and its numerical evaluation support this mechanism, but the integrated uniform remainder estimate needed to prove the asymptotic has not been completed here.

[Assumed, unresolved target] After specifying the energy-density gauge and a smooth observable or integrated formulation, the general bound A <= C M f(d) remains unproved. None of the obstructions above refutes it.

## 9. Sources

[Assumed] Local research context was supplied by specs/_bound_background.md with all updates, contrib/fable_bound_modular/README.md, and REPORT.md Sections 9.68 through 9.75. The quasi-free bound A <= 2M is cited from that context and not rederived.

[Proved in the cited literature] The exact half-line free-fermion spectral relation used in Section 2 is from [Viktor Eisler, On the Bisognano-Wichmann entanglement Hamiltonian of nonrelativistic fermions](https://arxiv.org/abs/2410.16433), J. Stat. Mech. (2025) 013101, Section 3 and Appendix B. The orthogonal-polynomial normalization, recurrence and generating function are recorded in [DLMF 18.19](https://dlmf.nist.gov/18.19), [DLMF 18.22](https://dlmf.nist.gov/18.22), and [DLMF 18.23](https://dlmf.nist.gov/18.23). The Boulware asymptotic expansion in Section 2.1 is this calculation's inference, not a theorem attributed to those sources.
