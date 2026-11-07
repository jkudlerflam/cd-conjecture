# C129 proof and literature audit

## 1. Conventions

[Proved] Let S^a=sigma^a/2 and use the oriented ring sum, including both oriented bonds when L=2. Write J=sqrt(J_perp^2+D^2). For J>0 set Delta=J_z/J. Let e_L(n,Phi) be the field-free XXZ sector ground energy with flux Phi, and set F_L(n)=min_Phi e_L(n,Phi). A convenient convention for the uniformly distributed flux is

    H_L(a) = sum_x { (J/2)[exp(ia) S_x^+ S_(x+1)^-
                              + exp(-ia) S_x^- S_(x+1)^+]
                         + J_z S_x^z S_(x+1)^z },   Phi=La.

Fluxes differing by 2pi are unitarily equivalent in each sector. This remains true for spin 1/2 because the possible minus sign of a 2pi on-site rotation is a scalar in its adjoint action. The directed-bond convention matters at L=2. The two DM terms then cancel in the physical Hamiltonian, as they must.

## 2. Step 1, removal of the field

[Proved] Every bond commutes with S_x^z+S_(x+1)^z. In the sector with n down spins, S_tot^z=L/2-n, so

    E_L(B)=min_(0<=n<=L) {e_L(n,Ltheta)-B(L/2-n)}.

Choose any minimizing n_*. The variational principle on the 2L ring, restricted to sector 2n_*, gives

    E_(2L)(B) <= e_(2L)(2n_*,2Ltheta)-B(L-2n_*).

Subtracting twice the preceding equality gives

    E_(2L)(B)-2E_L(B)
      <= e_(2L)(2n_*,2Ltheta)-2e_L(n_*,Ltheta).                 (2.1)

No differentiability in B or uniqueness of n_* is required. Spin flip sends (n,Phi) to (L-n,-Phi), while complex conjugation preserves n and reverses Phi. Consequently e_L(n,Phi)=e_L(L-n,Phi). It suffices to bound sectors n<=L/2 even when the physical minimizing sector is on the other side of half filling.

The field drops out of the numerator estimate. It does not drop out of the definition of M. Section 4 supplies the comparison needed for the denominator.

## 3. Step 2, the DM term and flux

[Proved] The transverse bond is

    (J_perp+iD) S_x^+ S_(x+1)^-/2
      +(J_perp-iD) S_x^- S_(x+1)^+/2.

For J>0 let theta=atan2(D,J_perp), so J_perp+iD=J exp(i theta). Number sites x=0,...,L-1 and define U=exp(-i theta sum_x x S_x^z). In U^dagger H U, the operator S_x^+ S_y^- acquires exp(i theta(x-y)). Every bulk hopping therefore becomes J/2, while the closing hopping becomes (J/2)exp(iLtheta). The diagonal interaction and field are unchanged. This proves the claimed mapping with flux Phi=Ltheta, and flux 2Phi on the doubled ring. Using UHU^dagger instead would reverse the chosen rotation convention.

This argument works at every L>=2. In particular it does not require an even ring or a sublattice rotation.

## 4. The optimized scale, including the field

[Proved] The gauge-optimized half-spread of the physical bond is already attained without a telescoping correction. Its four eigenvalues are

    J_z/4-B/2, J_z/4+B/2, -J_z/4-J/2, -J_z/4+J/2.            (4.1)

To check optimality, average any correction o tensor 1-1 tensor o over simultaneous rotations about z. Spectral diameter is convex and unitarily invariant, so this cannot increase it. Only o=t S^z need be considered, up to an irrelevant scalar. The parallel-spin eigenvalues in (4.1) do not change, and the two antiparallel eigenvalues become

    -J_z/4 +/- sqrt(J^2/4+t^2).

Increasing |t| can only enlarge the interval containing all four eigenvalues. The minimum is t=0. Thus M(B) is one half the range of (4.1). This proof also applies when J=0.

In particular, with M_0=M(0),

    M_0 = J/2                         if |Delta|<=1,
          J(1+|Delta|)/4              if |Delta|>=1,
    M(B)>=M_0>=J/2.                                             (4.2)

The field only expands the interval between the two parallel-spin eigenvalues, which proves the comparison. The quoted values 0.5, sqrt(1.09)/2, and 0.625 follow at B=0. The assertion M=J/2 throughout |Delta|<=1 would be false for sufficiently large B.

## 5. Step 3, exact twist cost

[Proved] Work in the uniform gauge H_L(a). Let Phi_m minimize the sector ground energy. Choose a lift a_m=Phi_m/L and any normalized ground vector psi of H_L(a_m) in sector n. Write its transverse part as T(a_m), and its derivative with respect to a as I(a_m). Direct expansion of the two hopping phases gives the operator identity

    H_L(a_m+alpha)-H_L(a_m)
      = (cos alpha-1)T(a_m)+sin alpha I(a_m).                    (5.1)

The scalar function g(alpha)=<psi,H_L(a_m+alpha)psi> has a minimum at zero. Indeed it is at least the ground energy at that flux, which is at least the globally minimal ground energy g(0). Hence g'(0)=<I(a_m)>=0. This argument also works when the minimizing eigenvalue is degenerate. No derivative of an eigenvalue branch is assumed.

Each transverse bond has norm J/2, so ||T||<=JL/2. The variational principle and (5.1) now give, for a representative delta of Phi-Phi_m with |delta|<=pi and alpha=delta/L,

    0 <= e_L(n,Phi)-F_L(n)
      <= (1-cos alpha)(-<T>)
      <= (JL/2)(1-cos(delta/L))
      <= (JL/2)(1-cos(pi/L)) <= J pi^2/(4L).                    (5.2)

The trial state is psi in the uniform gauge, with the corresponding on-site gauge transformation if a closing-bond gauge is used. Minimality also implies <T><=0 by considering small nonzero alpha in (5.1).

There is an explicit remainder if a Taylor statement is desired. For |delta|<=pi and L>=2,

    (JL/2)(1-cos(delta/L))
      = J delta^2/(4L)-J delta^4/(48L^3)+R_6,
    0 <= R_6 <= J |delta|^6/(1440L^5).                          (5.3)

For the trial expectation itself, its difference from -alpha^2<T>/2 has absolute value at most J|delta|^4/(48L^3). The quadratic upper bound in (5.2) needs no remainder.

Combining (2.1) and (5.2) at size 2L, and using e_L>=F_L, yields

    E_(2L)-2E_L <= F_(2L)(2n_*)-2F_L(n_*)
                      + JL[1-cos(pi/(2L))]
                <= F_(2L)(2n_*)-2F_L(n_*)+pi^2 M(B)/(4L).       (5.4)

The twist cost is paid only on the doubled ring.

## 6. A minimizing flux and positivity

[Proved] F_L(n) is the sector ground energy of

    K_L = sum_x {-J(S_x^x S_(x+1)^x+S_x^y S_(x+1)^y)
                       +J_z S_x^z S_(x+1)^z}.                 (6.1)

Thus Phi_m=Lpi modulo 2pi is a minimizer for all n and Delta. To prove this, write a vector in the spin configuration basis and replace every component by its modulus. Each directed hopping contribution to its energy is bounded below by minus its hopping magnitude times the product of the two component moduli. Summing gives <psi,H_L(a)psi> >= <|psi|,K_L|psi|>. At a=pi this lower bound is attained as a Hamiltonian. At L=2 the argument applies to the two hopping contributions separately, so it still holds.

For 0<n<L and J>0 the nearest-neighbor exchange graph is connected. Adding a scalar to -K_L gives an irreducible nonnegative matrix. Perron-Frobenius therefore gives a unique, strictly positive ground vector in every such sector for all real Delta. The fully polarized sectors are one-dimensional.

This proves ground-state uniqueness at the minimizing flux without invoking a twisted version of Yang-Yang. It does not prove that the published even-length asymptotic results apply at odd length. At even L the minimizing flux is zero, and the usual staggered unitary relates (6.1) to the periodic antiferromagnetic transverse convention. At odd L the minimizing flux is pi.

## 7. Free fermions, all lengths and all sectors

[Proved] Set Delta=0. Jordan-Wigner applied to (6.1) gives single-particle dispersion -J cos k, with boundary condition exp(ikL)=(-1)^(n-1). For odd n the momenta are integers times 2pi/L, and the occupied block is centered on zero. For even n they are half-odd-integer multiples of 2pi/L, and the occupied block is again centered on zero. Summing the geometric progression gives

    F_L(n)=-J sin(pi n/L)/sin(pi/L).                             (7.1)

This includes n=0 and n=L by continuity of the numerator, and L=2 with the directed-bond convention. For a general flux the momentum grid is displaced. The n lowest levels form a consecutive block, whose center can be chosen within pi/L of zero. If its displacement is beta, the energy is the right side of (7.1) times cos beta. At the edge of that interval two blocks tie, which is the open-shell degeneracy. Minimizing the flux sets beta=0 for either parity. Thus no open-shell term has been discarded.

Put rho=n/L and x=pi/L. Then

    L[F_(2L)(2n)-2F_L(n)]
      = J sin(pi rho) L[2 csc x-csc(x/2)]
      = J sin(pi rho) L tan(x/4) sec(x/2)
      = (pi J/4) sin(pi rho)[1+7pi^2/(48L^2)+O(L^-4)].          (7.2)

The Taylor remainder is uniform in rho because 0<=sin(pi rho)<=1. Both tan(x/4)/(x/4) and sec(x/2) increase for 0<x<=pi/2. It follows that the maximum over L>=2 is at L=2. Therefore

    F_(2L)(2n)-2F_L(n) <= (8-4sqrt(2)) M_0/L.                  (7.3)

The maximum in units of J is 4-2sqrt(2)=1.171572875... at L=2,n=1. At L=4 it is 0.861204779... at n=2, and the limit at half filling is pi J/4. Combining (7.3) with (5.4) proves the explicit, nonoptimal original-class bound

    L(E_(2L)-2E_L) <= [8-4sqrt(2)+pi^2/4] M(B)                 (7.4)

for J_z=0, all DM angles, fields, and L>=2.

## 8. Ferromagnetic and classical regions

[Proved] If J_z<=-J, the minimum eigenvalue of each field-free physical bond is J_z/4, attained on both parallel-spin vectors. Also -B S_tot^z>=-|B|L/2. A fully polarized vector in the favored direction attains both lower bounds simultaneously, independently of the DM angle. Hence

    E_L(B)=L[J_z/4-|B|/2],   E_(2L)-2E_L=0.                    (8.1)

This proves the original target with constant zero on Delta<=-1. At Delta=-1, (6.1) is the isotropic ferromagnet. Its symmetric spin-L/2 multiplet has one vector in every sector and saturates the local lower bound, giving F_L(n)=-JL/4 and sectorwise doubling defect zero.

[Proved] If J=0, the Hamiltonian is diagonal in the spin basis. Repeating a minimizing length-L cyclic spin word twice produces a length-2L word with twice its energy. Thus E_(2L)<=2E_L, including antiferromagnetic Ising interactions on odd rings and arbitrary B.

[Assumed, unresolved] Equation (8.1) does not prove the proposed sectorwise droplet inequality for Delta<-1. That stronger inequality is unnecessary for the original target in this region.

## 9. A dilute-sector bound uniform in Delta>=-1

[Proved] Write eta_x=1/2-S_x^z and A=sum_x eta_x eta_(x+1). In sector n, (6.1) becomes

    K_L=J Delta L/4-J Delta n+T_L+J Delta A,
    T_L=-(J/2)sum_x(S_x^+ S_(x+1)^-+S_x^- S_(x+1)^+).

On a single bond,

    -(J/2)(S^+ S^-+S^- S^+)+J Delta eta_x eta_y
                    +(J/2)(eta_x+eta_y) >= 0                  (9.1)

for Delta>=-1. On 00 the eigenvalue is zero, on the one-particle block the eigenvalues are 0 and J, and on 11 it is J(Delta+1). Summing gives

    F_L(n)>=J Delta L/4-J(Delta+1)n.                            (9.2)

Use the free ground state with k=2n particles on m=2L sites as a trial state for F_(2L)(2n). Its occupied momenta are centered on zero as in Section 7. Set

    A_k=sin(kpi/m)/sin(pi/m),
    d_k=k-A_k=sum_(occupied p)(1-cos p).

Its hopping energy is -J A_k. Its nearest-neighbor density expectation, by the two-point determinant of a Slater state, is

    <sum_x eta_x eta_(x+1)>=(k^2-A_k^2)/m.                     (9.3)

For completeness, <eta_x>=k/m and <c_x^dagger c_(x+1)>=A_k/m, so each distinct neighboring pair has expectation (k/m)^2-(A_k/m)^2. Here m=2L>=4, so all these pairs are distinct sites. Summing proves (9.3).

Subtract twice (9.2) from this trial energy. For Delta_+=max(Delta,0),

    F_(2L)(2n)-2F_L(n)
       <= J d_k + J Delta_+ (k^2-A_k^2)/m.                     (9.4)

Since k<=m/2, 0<=A_k<=k. The elementary inequality 1-cos p<=p^2/2 and the sum of the centered squared momenta give

    d_k <= pi^2 k(k^2-1)/(6m^2),
    k^2-A_k^2=(k-A_k)(k+A_k)<=2k d_k.

Using J/M_0<=2 and J Delta_+/M_0<=4, substituting m=2L,k=2n gives

    L[F_(2L)(2n)-2F_L(n)]/M_0
      <= pi^2 n(4n^2-1)/(6L)
                         +2pi^2 n^2(4n^2-1)/(3L^2)
      <= 2pi^2 n^3/L.                                        (9.5)

The last step uses n/L<=1/2. Thus the sectorwise constant is 2pi^2 on n^3<=L, uniformly for every Delta>=-1. With the physical flux and field, add pi^2/4 to this constant using (5.4), provided the minimizing sector satisfies min(n_*,L-n_*)^3<=L.

For every fixed n, the stronger first line of (9.5) tends to zero as L increases. At n=1 it gives pi^2/(2L)+2pi^2/L^2. A constant valid for all L>=2n is pi^2(4n^2-1)/4. These are bounds on the doubling defect, not a derivation of the dilute-magnon asymptotic coefficient.

[Assumed, unresolved] This does not cover n/L small with n^3/L unbounded. In particular it does not establish suggested piece (b) in its full density-uniform form.

[Proved] There is also an exact saturation region. For B>=J(Delta+1) and Delta>=-1, (9.2), the diamagnetic inequality, and the sector field term give

    e_L(n,Ltheta)-B(L/2-n)
       >= J Delta L/4-BL/2+n[B-J(Delta+1)].

The fully polarized sector attains this lower bound at n=0. Spin flip gives the corresponding statement for B<0. Hence the original defect is zero whenever |B|>=J+J_z with J_z>=-J, for every L and DM angle.

## 10. An explicit ferromagnetic crossover band

[Proved] Put Delta=-1+epsilon, epsilon>=0. In sector n the lower bound (9.2) is

    F_L(n)>=J Delta L/4-J epsilon n.

For the doubled ring use the normalized equal-amplitude sum over all configurations with k=2n down spins on m=2L sites. This is a vector in the ferromagnetic symmetric multiplet, so K_m(-1) has energy -Jm/4. Counting configurations gives

    <eta_x eta_y>=k(k-1)/[m(m-1)]  for x!=y,
    <sum_x S_x^z S_(x+1)^z>=m/4-k+k(k-1)/(m-1).

The trial energy is therefore J Delta m/4-J epsilon k+J epsilon k(k-1)/(m-1). Subtracting twice the lower bound gives

    F_(2L)(2n)-2F_L(n) <= J epsilon (2n)(2n-1)/(2L-1),
    L[F_(2L)(2n)-2F_L(n)]/M_0
       <= 2epsilon L(2n)(2n-1)/(2L-1) <= epsilon L^2.          (10.1)

The last inequality uses 2n<=L and 2(L-1)<=2L-1. It follows that the original physical defect is bounded by (A+pi^2/4)M(B)/L throughout epsilon L^2<=A. This includes odd L and all sectors. At epsilon=0 the stronger exact result in Section 8 applies.

[Assumed, unresolved] For epsilon L^2 growing without bound this estimate is insufficient. Combining it with pointwise Bethe asymptotics does not justify a uniform bound on the intervening scales.

## 11. Literature audit and normalization

[Proved, literature statement] The primary source for the expansion is [Kozlowski, arXiv:1508.05741v2](https://arxiv.org/html/1508.05741). The paper assumes even L and uses Pauli matrices, so its exchange parameter is J/4 in the present convention. Its Proposition 3.4 concerns finite Fermi boundaries. Theorem 3.5 supplies half-filled convergence without a rate for -1<Delta<=1. Proposition 4.2 and Corollary 4.3 already give the energy correction below half filling. Theorem 3.2 explicitly makes its leading error and size threshold uniform on compact density ranges below 1/2 at fixed Delta. It does not state uniformity in Delta. The discussion after Lemma 4.1 explicitly leaves energy control near half filling unresolved.

These corrections to TOOLBOX Sections 7 and 8.3 are consequential. The density interval proposed in piece (a) contains a boundary not covered by the finite-size expansion. Even at Delta=0, where Section 7 proves the answer directly, one cannot obtain the half-filled Casimir term from Theorem 3.5 alone.

[Proved, literature statement with access limitation] The [Yang-Yang paper](https://journals.aps.org/pr/abstract/10.1103/PhysRev.150.321) establishes the fixed-sector ground-state Bethe hypothesis on a cyclic chain. Its original PDF was inaccessible through the browsing tool. The Delta range checked here is therefore the range restated in Kozlowski's introduction, Delta>=-1 in the antiferromagnetic transverse convention, with the endpoint understood by a limit. The real-root uniqueness argument quoted there covers -1<Delta<=0, with continuation identifying the ground branch beyond zero. This is not a theorem for arbitrary flux, and this run does not claim to have checked such an extension in the original 1966 text. Section 6 supplies the needed uniqueness and flux minimization directly.

## 12. Constants and the uniformity issue

The following is an audit of what would have to be controlled to make the finite-boundary proof uniform. Equation numbers below refer to the HTML rendering, whose numbering differs by one from the PDF in parts of Sections 3 and 4. The named results and displayed quantities identify the statements without relying on that discrepancy.

[Proved, dependency audit] The inputs entering Proposition 3.4 are the anisotropy, limiting density, fixed number of particle-hole excitations, bounded charge shift s, expansion order r, distance of particle quantum numbers from their admissible endpoints, and the smaller complex strip on which the remainder is measured. For the sector ground state, the excitation count and s are zero. The relevant constants are not a single C(Delta,rho). They arise in the following steps.

1. Root localization uses Lambda from Proposition 2.2. The available interval for the massless particle quantum numbers has margin proportional to (pi-zeta)(1/2-rho)/pi. Root bounds and this margin enter before the contour argument.
2. The complex strip uses the printed kappa_Delta in (3.11). It is zeta/4 for Delta=cos zeta or cosh zeta, and 1/4 at Delta=1. Bounds B on the counting function appear in (3.19). Their inputs include Lambda and the distances to poles of the bare momentum and scattering kernel.
3. Inverse-function neighborhoods use eta, the horizontal enlargement, and the contour height alpha. The contour separation c_0 in (3.57) depends on the limiting counting function. Its inverse derivative is controlled by a positive lower bound on the root density over the relevant compact rapidity interval.
4. Resolvent sup norms and derivative bounds enter (3.45)-(3.46), the inverse-function comparison (3.58)-(3.59), and the contraction estimate for the nonlinear integral equation. Their numerical values depend on Delta and the integration interval.
5. Endpoint expansions divide by rho_Bethe(q|q), explicitly in (3.105)-(3.108). Repeated differentiation and inversion amplify this dependence.
6. The higher-order remainder involves derivatives of R(lambda,xi_hat^(-1)(z))/xi_hat'(xi_hat^(-1)(z)) in (3.122)-(3.134), contour separation, and integrals of powers of t against log(1+exp(-t)). The contour part away from the real axis is exponentially small with an exponent set by alpha. This is where fixed-parameter analyticity becomes a quantitative remainder estimate.
7. Theorem 3.5 instead introduces Lambda_epsilon from tail tightness, further constants C,C',C'', and a size threshold L_0 depending on the requested accuracy. Its proof takes L to infinity before removing rapidity cutoffs. It supplies no bound relating L_0 to epsilon that would imply an O(L^-1) energy correction.

[Proved] Arbitrary constants chosen to make an inequality true have no canonical continuity property. The paper does not define continuous functions C(Delta,rho) and L_0(Delta,rho). In particular L_0 is integer-valued and is not expected to be continuous. The useful question is whether there are locally bounded majorants, not whether the choices of constants are continuous.

[Proved, elementary analytic input] On a compact set of finite Fermi boundaries in one rapidity chart, the integral-equation data vary continuously in operator norm. Rescale a finite interval to [-1,1]. The continuous kernel and its derivatives then vary uniformly. If A_0=1+K_0 is invertible and ||A_0^(-1)(A-A_0)||<1, the Neumann series for [1+A_0^(-1)(A-A_0)]^(-1) proves continuity of A^(-1). Thus the dressed functions and their derivatives vary continuously. The strict derivative of the density-to-boundary map gives local continuity of q by the inverse function theorem. Positive root density on a common compact interval then has a positive lower bound. These observations control the linear resolvent and inverse-function parts locally, wherever their stated hypotheses hold.

[Assumed, unresolved parameter extension] This does not by itself give a common L_0 for the actual discrete Bethe roots. The initial normal-family convergence argument in Theorem 3.2 must be repeated for sequences with Delta varying, using a common root bound, strip, and limit identification. Only after that step can the contraction and Taylor estimates be made uniform. Density uniformity at fixed Delta is stated in the theorem. Uniformity in Delta is not. It appears accessible away from the singular boundaries, but this run does not promote that extension to a proved theorem or an explicit constant.

At Delta=1, the rapidity variables used on the two sides need rescaling before the kernels have a regular limit. The printed strip widths shrink to zero in those coordinates. A covering argument that treats those widths as already uniform across Delta=1 is invalid.

[Proved, additional analytic-domain check] There is also a literal strip-width issue for negative Delta in the printed formulation. For Delta=cos zeta, the scattering derivative is

    theta'(lambda)=sin(2zeta)/[sinh(lambda+izeta)sinh(lambda-izeta)].

For pi/2<zeta<pi its nearest nonremovable poles are at lambda=+/-i(pi-zeta). The finite counting derivative contains shifted copies at lambda_a+/-i(pi-zeta). Thus a general counting function cannot be holomorphic on the whole printed strip |Im lambda|<2kappa_Delta=zeta/2 when zeta>2pi/3. The bare momentum has poles at +/-izeta/2. A safe common strip must be strictly narrower than min(zeta/2,pi-zeta). This is a directly checked domain restriction, not a refutation of the real-axis asymptotic result. Shrinking the strip repairs this particular issue on compact subsets away from Delta=-1. It makes the loss of uniform complex width as Delta approaches -1 explicit.

[Proved] At half filling in the massless regime, q tends to infinity and rho_Bethe(q|q) tends to zero. The endpoint denominators in item 5 and the common compact rapidity domain in item 3 are lost. Compactness of the parameter rectangle in (Delta,rho) does not restore these hypotheses. At Delta=-1, the scattering kernel and the corresponding inverse-operator bounds also become singular. These are separate obstructions to using the stated asymptotics on the requested closed parameter region.

## 13. The energy term and what it implies

[Proved from the cited finite-boundary expansion] Set the charge shift and the particle-hole count to zero. Then xi_1=1/2, so both first endpoint shifts vanish. Corollary 4.3 reduces to

    F_L(n)=L e(Delta,n/L)-a(Delta,n/L)/L+O(L^-2),
    a=epsilon_eff'(q)/(12 rho_Bethe(q|q)).                      (13.1)

This statement is used here on the paper's even-length domain and below half filling. One may choose the auxiliary field whose dressed energy vanishes at q, because changing the field in a fixed sector adds only an extensive term. At this field epsilon_eff=epsilon. A particle-hole displacement by one quantum number changes rapidity by 1/[L rho_Bethe(q|q)] and physical momentum by 2pi/L. Consequently the physical velocity is

    v=epsilon'(q)/[2pi rho_Bethe(q|q)],
    a=pi v/6.                                                  (13.2)

The factor 2pi is needed because the paper's dressed momentum is normalized by the quantum number density. The Pauli-to-spin factor J_paper=J/4 must also be included. No twisted charge formula is required for the target after Section 6.

As an independent normalization check, for Delta=0 below half filling the expansion of (7.1) gives a=pi J sin(pi rho)/6. The group velocity of -J cos k at k=pi rho is J sin(pi rho). This agrees with (13.2).

[Proved, conditional finite-size implication] Suppose a region has an actual estimate

    |F_L(n)-L e(Delta,n/L)+pi v(Delta,n/L)/(6L)|
                      <= k M_0/L^2  for every L>=L_0.           (13.3)

Then the densities of (L,n) and (2L,2n) agree, the extensive terms cancel, and

    L[F_(2L)(2n)-2F_L(n)]/M_0 <= pi v/(4M_0)+9k/(4L).

Using the recorded C109 bound v<=4e M_0 gives pi e+9k/(4L). For smaller L there is a safe sector-preserving gluing bound 4M_0. Indeed take a tensor product of two sector ground vectors, change the two seam bonds, and bound each old-to-new bond expectation difference by its spectral diameter 2M_0. Therefore a valid sector constant under (13.3) is

    max(4L_0, pi e+9k/(4L_0)),                                 (13.4)

with the obvious rounding if L_0 is not an integer. An additional pi^2/4 covers physical fluxes and fields. This formula displays why a statement with unspecified k and L_0 is not the explicit constant requested in C129.

[Proved, fixed-parameter consequence] For each fixed Delta>-1 and each fixed density strictly below 1/2 on a sequence of admissible even lengths, the expansion implies bounded L times the sector doubling defect. On compact density intervals within the finite-boundary regime at fixed Delta, the uniform leading estimates and the bounded contour ingredients give the same qualitative conclusion. This uses a nonexplicit threshold and nonexplicit remainder constants. It gives no global supremum over Delta and does not include half filling.

[Assumed, unresolved] Extending the asymptotics to the odd-length minimizing flux requires checking the Bethe-vector identification and the expansion with that boundary condition. The centered logarithmic equations used in the numerical code pass the odd-length ED tests, but a numerical match is not a substitute for the extension. Even if this extension and local parameter uniformity were supplied, piece (a) as written would still require new control at half filling for -1<Delta<=1.

## 14. The remaining literature and the proposed decomposition

[Proved, literature audit] [DKKMT, arXiv:2012.11675](https://arxiv.org/html/2012.11675), Theorem 2, has a per-site O(1/N) error and excludes a range near n=N/2. Theorem 22 improves the error to O(1/[zeta(N-2n)N]) under additional hypotheses, again excluding the half-filled endpoint. It does not give the uniform half-filled Casimir bound. Taking a quantum Hamiltonian limit also requires control of the transfer-parameter derivative of the error. A bound on the transfer eigenvalue alone does not justify that derivative. Its condensation estimates are useful inputs for a parameter-uniform argument, but no such completed argument is supplied here.

[Proved, literature audit] [Aizenman-Duminil-Copin-Warzel, arXiv:2002.02543](https://arxiv.org/html/2002.02543), Theorem 1.4, proves two infinite-volume XXZ ground states with Neel order for Delta>1. This statement is not a uniform finite-ring energy estimate, and the theorem for the spin-S projector model elsewhere in the paper should not be silently substituted for an XXZ Casimir estimate. Neel order by itself does not give a bound uniform as Delta decreases to 1.

[Proved] Sections 2-6 give the twist estimate directly, so no unverified numerical coefficient from Tasaki or the Kohn budget is needed. Those results motivate the same variational mechanism and are not used as additional assumptions.

The proposed pieces have the following status.

(a) [Assumed, unresolved in the requested form] The energy term below half filling is already a corollary in the literature. Parameter-uniform constants are not extracted here. The half-filled part of the proposed compact set needs a different estimate.

(b) [Proved in a restricted range] Section 9 covers n^3<=L, and all fixed-n sequences, uniformly in Delta>=-1. The whole small-density region is unresolved.

(c) [Proved in a restricted range] Section 10 covers (Delta+1)L^2<=A with sector constant A. The full crossover remains unresolved.

(d) [Assumed, unresolved uniformly] The finite-boundary energy expansion applies away from half filling at fixed Delta>1. Neither the Heisenberg half-filled estimate nor a bound uniform down to Delta=1 follows from the cited ordering theorem. Large-Delta and half-filled finite-size constants have not been extracted.

(e) [Proved for the original target] The original unrestricted defect is zero for Delta<=-1 by Section 8. The stronger droplet-sector lemma was not proved.

[Proved] None of the elementary reductions needs a vacuum-compatible subsequence. What remains unproved is the core bound itself, including the extension of the cited even-length asymptotics where needed. Thus this run does not establish C_U1<infinity or the stated partial success (a)+(b)+(e).
