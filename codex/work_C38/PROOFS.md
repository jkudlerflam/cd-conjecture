# C38. What a scale hierarchy does and does not prove

[Proved, scope] Bandwidth separation alone does not imply separation of Fermi momenta or suppression of cross-flavor effects relative to the slow band. This note gives a quantitative obstruction in the natural resolvent limit of the folded bus, an explicit all-n lacunary target, and sufficient matrix bounds for preserving that target. The realization of those target bounds by the 23-letter frozen core or C32's 28-letter mobile core remains unproved. This is a failure of the proposed inference, not a no-go theorem for either family.

[Assumed, infrared scope] Operator dimensions below use the free-fermion dictionary and its weakly interacting Gaussian continuation, as in C20 and ir_control. All matrix and arithmetic statements are independent of that input. Exact total U(1) is used throughout. We do not attempt the parity-only operator classification. Density interactions remain allowed and are not called unwanted resonances.

## 1. Energy separation is not a gap for eliminating a flavor

[Proved] Let S be any invertible positive diagonal matrix, independent of momentum, and let F(k) be Hermitian. Then

    h(k)=S F(k) S,
    det h(k)=(det S)^2 det F(k).

Consequently h and F have the same Fermi zeros, with the same determinant multiplicities and inertia at every k. Congruence is not unitary diagonalization and need not preserve velocities, but it preserves the zero set. Taking S_aa arbitrarily far apart cannot change any resonance of that zero set. If F is diagonal with nonzero entries of order one at some k, the energies of h at that k have arbitrary prescribed scale separation.

[Proved] For example, h_a(k)=t_a[mu_a-2cos(k-theta_a)], with any positive geometric sequence t_a, has roots determined by mu_a and theta_a, independently of t_a. Choosing every theta_a=0 retains all inter-flavor Cooper-transfer resonances at every ratio t_a/t_b. Choosing the target in Section 4 removes them. Thus the same energy hierarchy is compatible with both outcomes.

[Proved] The cross-scale Schur complement is already instructive for two flavors:

    h(k) = [[t_1 f_1(k), sqrt(t_1 t_2) e(k)],
            [sqrt(t_1 t_2) conj(e(k)), t_2 f_2(k)]].

Where f_1 is nonzero, eliminating its row gives the exact scalar equation

    t_2 [f_2(k)-|e(k)|^2/f_1(k)] = 0.                  (1)

The relative correction is |e|^2/f_1, with no factor t_2/t_1. For t_1=1, t_2=Q^(-4), e=1/10, and f_1=1 at the slow reference root, the correction divided by t_2 is 1/100 for every Q. The off-diagonal entry divided by the fast bandwidth is only 1/(10Q^2). It still does not become small relative to the slow band as Q grows.

[Proved] Each fast flavor itself crosses the Fermi level. Its resolvent therefore has a pole at its own momentum even when its bandwidth is large. One may use (1) in a neighborhood of a different, already separated root, but not as a globally gapped elimination. A lower bound on |f_1| in that neighborhood is a momentum-separation hypothesis. Inferring it from the desired conclusion would be circular.

[Assumed, physical interpretation] Interacting fast flavors likewise remain low-energy fields. Different velocities do not remove momentum-conserving Cooper transfer from the marginal free-fermion operator list. This is why the matrix calculation cannot be repaired by integrating out an entire wide band as a massive mode.

## 2. The natural Gram hierarchy has unsuppressed cross-scale terms

[Proved, inherited setup] In the large-express-energy tree limit of both folded-bus cores the leading hopping is a nonzero scalar times the Gram matrix

    Gamma_ij = integral_0^infinity f_i(t) f_j(t) dt,
    f_j(t) = integral_0^t K_0j(u) K_j0(t-u) du.          (2)

For the frozen stub, K(t)=exp[-t(a I-A_path)]. For the mobile positive core, K(t)=exp[-t(m I+B)], with the bounded tridiagonal B specified in C32 PROOFS (19). Its nearest-neighbor off-diagonal entries are -1. Both ports reach the same stub endpoint. We use this existing Gram limit, not a new claim about exact equality at finite express energy.

[Proved] Number the stub sites j=0,...,n-1. In the frozen case rescale t=u/a. The unique shortest path gives

    [exp(u A_path/a)]_0j = (u/a)^j/j! + O_n(a^(-j-2))

on bounded u intervals. Substituting in the convolution gives

    f_j(u/a)
      = a^(-2j-1) exp(-u) u^(2j+1)/(2j+1)!
        + O_n(a^(-2j-3)).                              (3)

The missing odd correction follows from bipartiteness of the path. For the mobile core replace a by m. Its diagonal entries in B allow one extra step, so the relative remainder is O_n(1/m), rather than O_n(1/a^2), with the same leading term.

[Proved, control of integration] These are also asymptotic expansions of (2), not just short-time formal series. For any fixed n, the matrix-exponential series after its first nonzero power is bounded by a fixed polynomial in u times exp(2u/a) in the frozen case. For the mobile core replace 2 by ||B||. Choose a or m larger than four times that fixed norm. The product in the Gram integral is then bounded by a fixed polynomial times exp(-u). Taylor remainders have the same integrable bound after extracting their displayed powers. Dominated convergence, or the integrated remainder estimate, proves the expansions entry by entry. All constants can depend on n.

[Proved] Put S_a=diag(1,a^(-2),...,a^(-2n+2)). Equations (2)-(3) give

    Gamma = a^(-3) S_a [C+O_n(a^(-2))] S_a,             (4)
    C_ij = (2i+2j+2)! /
           [(2i+1)!(2j+1)! 2^(2i+2j+3)].               (5)

For the mobile core the corresponding expression has a=m and O_n(m^(-1)). C is positive definite because it is the Gram matrix of the distinct monomials exp(-u)u^(2j+1)/(2j+1)!. Thus there are positive finite constants c_n,C_n with

    c_n a^(-3) S_a^2 <= Gamma <= C_n a^(-3) S_a^2

for sufficiently large a. The min-max principle places its ordered eigenvalues between constants times a^(-3-4j). Their adjacent ratios can therefore be made arbitrarily small at each fixed n. This is an actual hierarchy, not a numerical exponential fit.

[Proved] However the normalized off-diagonal entries do not vanish:

    Gamma_ij / sqrt(Gamma_ii Gamma_jj)
         -> C_ij/sqrt(C_ii C_jj).                       (6)

For the first two indices the exact rational matrix is

    C = [[1/4, 1/8], [1/8, 5/32]],
    C_01^2/(C_00 C_11)=2/5,
    C_11-C_01^2/C_00=3/32.                              (7)

At n=2, the slow eigenvalue is (3/32)a^(-7)(1+o(1)), whereas ignoring the fast-site cross term would give (5/32)a^(-7). The relative error of that approximation is not small as the scale ratio tends to zero. This is also the leading Schur pivot in the explicit natural hierarchy. For n>2, (7) describes the first elimination step, without claiming it is the final slowest eigenvalue.

[Proved] The adjacent normalized correlations become stronger at large index:

    C_(j,j+1)^2/[C_jj C_(j+1,j+1)]
      = (4j+4)(4j+3)/[(4j+6)(4j+5)] -> 1.              (8)

This follows by cancellation of factorials in (5). Strong scale separation and strong normalized cross correlations coexist in this resolvent mechanism.

[Proved limitation] One can diagonalize the complete Gram matrix and call its eigenvectors flavors. That correctly incorporates (7). It does not diagonalize the loop, onsite, and finite-express corrections simultaneously. The resulting matrices require their own bounds and resonance tests. Equation (7) disproves an argument that identifies flavors with successive stub distances and treats their mutual contributions as relatively vanishing. The congruence argument of Section 1 explains why changing to the eigenbasis does not, by itself, prove nonresonance either.

[Numerical, n=2, Lyapunov dimension 8] checks.py evaluates (2) without a small-time truncation at a or m=8,16,32,64,128. For the frozen path the squared correlation is 0.4107004, 0.4026610, 0.4006644, 0.4001660, 0.4000415. For the positive mobile stub it is 0.3290435, 0.3582882, 0.3772264, 0.3880779, 0.3938973. These approach the proved value 2/5. They are double-precision consistency checks, not error certificates for the full microscopic bond.

## 3. The loop displacement depends on cross-scale matrix elements

[Proved, inherited structure] In the real-hopping subfamilies used in C33 and the small-loop studies, the leading hopping T is real, and the onsite quadratic matrix is imaginary antisymmetric. A common real field is added separately. Write

    T=-S+K,  S=S^T>0,  K=-K^T real,
    M=iA,    A=-A^T real.

Then

    h(k)=epsilon I-2S cos k+i[A+2K sin k].               (9)

The real source formulas in C33 fast_loop.py explicitly have this form. More general allowed phases need not preserve it. Nothing below asserts (9) for every parameter choice of the bond.

[Proved] The strict Hermitian Gram approximation sets K=0. At any A the resulting h(k) is even in k. Every pair of reflected roots then has center zero modulo pi. If instead A=0, complex conjugation relates h(k) to h(-k), with the same result for the zero set. In this class, a nonzero loop onsite matrix together with a finite-express antisymmetric hopping can displace centers. Keeping A while discarding K because it is a small correction destroys the mechanism one is trying to prove.

[Proved] This can be quantified in a two-flavor eigenbasis S=diag(s_1,s_2). With A_12=a and K_12=kappa,

    det h(k)=(epsilon-2s_1 cos k)(epsilon-2s_2 cos k)
                -(a+2kappa sin k)^2.                   (10)

The term odd in k is -4a kappa sin k. Arbitrarily strong separation between s_1 and s_2 does not create it. It vanishes whenever a or kappa vanishes.

[Proved, local perturbation formula] For distinct s_a>0 and 0<epsilon<2 min s_a, let

    k_a=arccos[epsilon/(2s_a)].

Take A=t A_1 and K=u K_1 as two independently small matrix perturbations. At each unperturbed root +/-k_a the other diagonal entries are nonzero. The zero is simple, so its root and center have Taylor expansions in t,u. The first-order shift vanishes since real antisymmetric matrices have zero diagonal in this real eigenbasis. The mixed term in the center is

    theta_a(t,u)
      = (2tu/epsilon) sum_(b!=a)
            (A_1)_ab (K_1)_ab/(s_a-s_b)
          + terms of total degree at least three.       (11)

Terms proportional to t^2 alone or u^2 alone do not shift the center, by the two reflection symmetries above. The remainder constant in (11) depends on all unperturbed root separations and velocities; it is not uniform in n or epsilon.

[Proved, derivation] At k=+/-k_a the other diagonal entry is

    d_b=epsilon(1-s_b/s_a).

The second-order scalar Schur term is

    -sum_b |i[t(A_1)_ab+2u(K_1)_ab sin k]|^2/d_b.

Its mixed numerator is 4tu(A_1)_ab(K_1)_ab sin k. Divide by d_a'=2s_a sin k when solving for the shifted root. The sine cancels and the same mixed shift occurs at both signs. Averaging the two roots gives (11). For n=2 it agrees directly with (10).

[Proved] Formula (11) is a signed sum of matrix-element products and spectral denominators. It is not a function of the ordered s_a alone. Contributions from s_b>s_a and s_b<s_a have opposite denominator signs. They can cancel. Summing (11) over a gives zero, since every pair cancels. None of this proves that individual centers vanish or coincide in the actual core; it identifies the missing all-n inequality that an argument from the hierarchy would have to supply.

[Proved] In particular, for normalized elements A_ab=alpha sqrt(s_a s_b), K_ab=kappa_0 sqrt(s_a s_b), their pair contribution is

    2 alpha kappa_0 s_a s_b/[epsilon(s_a-s_b)].

When s_a/s_b grows and epsilon is proportional to s_b, this tends to a nonzero multiple of alpha kappa_0. The small energy ratio does not suppress the relative root displacement. Bounds on absolute loop amplitudes sufficient to preserve crossings do not establish its sign, differences between flavors, or its separation from higher corrections.

## 4. A lacunary target that does solve the arithmetic problem

[Definition] Fix r=1/16 and b=1/64. For each n>=1 prescribe normalized momenta x=p/(2pi) by

    x_(a,L)=-1/4+b r^(2a-2),
    x_(a,R)=+1/4+b r^(2a-1),       a=1,...,n.           (12)

All momenta are measured in the cell coordinate. The offsets form one geometric sequence indexed by the full list of 2n endpoints, not merely by centers or bandwidths. Let

    d_n=2pi(1-2r)b r^(2n-1)
       =(7pi/256) 16^(-(2n-1)).                         (13)

[Proved] Every two distinct endpoints are separated modulo 2pi by at least d_n. Every two distinct unordered endpoint pairs, allowing repeated endpoints, have different sums modulo 2pi, separated by at least d_n. Thus this stronger finite test includes pair transfer, exchange, and umklapp, as well as quadratic flavor conversion.

[Proved, pair-sum proof] Cancel common endpoints from two pair multisets. Let b r^j be the largest remaining offset. If it appears once, the opposite side contains at most two offsets, each at most r times it. Its signed contribution therefore has magnitude at least (1-2r)b r^j. If it appears twice the lower bound is stronger. Since j<=2n-1, the offset difference is at least (1-2r)b r^(2n-1). Its magnitude is at most 2b. The difference of the base sums in (12) is a multiple of 1/2. For an integer base difference the circular distance is the offset difference. For a half-integer base difference the distance is at least 1/2-2b. Both exceed the bound in (13) divided by 2pi. The endpoint proof is the same argument with one term on either side.

[Proved, general lacunarity statement] The dominance part of this proof works for any decreasing positive sequence with successive ratio at most r<1/2, for differences of two sums with two terms each. Modular wraparound must still be controlled separately. Lacunarity of hopping eigenvalues is not lacunarity of these endpoints. Lacunarity of the centers alone also leaves the width-dependent comparisons uncontrolled.

[Proved] The centers and half-widths of (12) are

    theta_a=pi b(1+r)r^(2a-2),
    k_a=pi/2-pi b(1-r)r^(2a-2).                         (14)

The centers are distinct with ratio r^2. The widths lie between pi/2-pi/64 and pi/2. In particular sin k_a>1/2. A finite explicit reference Hamiltonian is

    h_a(k)=2t_a[cos k_a-cos(k-theta_a)],
    t_a>0 arbitrary.                                   (15)

Its roots are (12) and its velocities have magnitude 2t_a sin k_a>t_a. Taking t_a=Q^(-(a-1)) with any Q>2 adds a geometric energy hierarchy without changing the resonance proof.

[Proved, arithmetic certification] checks.py enumerates the normalized momenta and all pair sums with fractions.Fraction at n=1,2,3,4,8,16,23,32,64,128. It sorts the pair sums on the rational circle and checks all adjacent gaps, including the wrap gap. The largest test has 256 endpoints and 32896 pair sums. All tests have minimum pair-sum gap equal to 240/7 times our conservative lower bound. This is exact finite arithmetic supporting the all-n proof, not numerical evidence for realization by a folded bus.

[Assumed, missing realization] Implementing (15), or a sufficiently close approximation in the norm below, requires control of both the flavor-dependent fields and hopping phases in the band basis. C26 supplies such control at the abstract reflection-cell level. No fixed-alphabet implementation of (12) by the 23-letter or C32 core is proved here. Increasing the allowed fixed local dimension does not by itself specify that implementation. C39 is investigating local position-dependent readout, which is a separate task.

## 5. Explicit control of all quadratic cross terms

[Definition] Let t=diag(t_1,...,t_n), and let h_0 be (15). For a translation-invariant C^1 Hermitian correction R(k), including every inter-flavor entry, put

    F(k)=t^(-1/2) R(k) t^(-1/2),
    eta=max(sup_k ||F(k)||, sup_k ||F'(k)||).

[Proved] The sufficient bound

    eta <= d_n^2/100000                                  (16)

preserves exactly 2n simple Fermi zeros. Each moves by at most 2eta radians, and all endpoint and pair-sum margins remain at least d_n/2. Every physical Fermi velocity has magnitude at least min_a(t_a)/2. The result allows arbitrary C^1 range of R, not just nearest-neighbor perturbations.

[Proved, detailed bounds] Write the congruent normalized symbol as D(k)+F(k), with D_aa=2[cos k_a-cos(k-theta_a)]. Around each reference root take radius rho=d_n/16. These intervals are disjoint. On its own interval, |D_aa'|>=1 because all endpoints lie within pi/32 of +/-pi/2 and theta_a is small. At any k, factorization gives

    |D_bb(k)|=4 |sin((k-p_bL)/2) sin((k-p_bR)/2)|.

The two roots of a flavor have circular separation at least pi-pi/32. If z is the distance to the nearer root, the distance to the farther is at least half that separation. Using sin(z/2)>=z/pi for 0<=z<=pi, the displayed product is at least z/2. On an a-root interval the other roots are at least 15d_n/16 away, hence |D_bb|>=d_n/4 for b!=a. Outside all root intervals every |D_bb|>=rho/2.

[Proved] Since eta<d_n/8, the complementary block D_B+F_BB has inverse norm at most 8/d_n. Its scalar Schur correction is

    sigma_a=F_aa-F_aB (D_B+F_BB)^(-1) F_Ba.

Thus

    |sigma_a| <= eta+8eta^2/d_n <= 2eta,
    |sigma_a'| <= eta+16eta^2/d_n
                    +64eta^2(2+eta)/d_n^2 < 1/4.         (17)

Here ||D_B'||<=2, and differentiation of an inverse gives the last term. The estimates use only (16) and 0<d_n<1.

[Proved] The scalar D_aa+sigma_a is strictly monotone on its interval, its endpoint values have opposite signs since 2eta<rho, and it therefore has one simple zero. Its displacement is at most 2eta by |D_aa'|>=1. Outside the intervals the diagonal inverse bound and eta<rho/2 exclude any zero. Moving four endpoints changes a quartic mismatch by at most 8eta, which is below d_n/2. This also preserves distinct centers because their doubled differences are among the tested pair sums.

[Proved, velocity bound] At a perturbed root take a normalized-symbol kernel vector xi with xi_a=1. Equation (17) gives ||xi_B||<=8eta/d_n. The physical kernel vector is proportional to t^(-1/2)xi. Its eigenvalue derivative is

    e'(p)= [D_aa'(p)+sigma_a'(p)]/[xi^dag t^(-1) xi].

The numerator has magnitude at least 3/4. The denominator is at most [1+(8eta/d_n)^2]/min t_a. This proves the stated min(t_a)/2 bound. The formula follows by differentiating the Schur identity; the complementary residual is zero, so derivatives of xi do not contribute.

[Proved] A more conservative sufficient condition in an ordinary norm is

    max(sup ||R||,sup ||R'||)
          <= (min t_a) d_n^2/100000.                    (18)

This follows by multiplying R on both sides by t^(-1/2). It explicitly pays for the slowest scale. No factor representing fast/slow energy separation makes a normalized cross entry small automatically.

[Proved] Equation (17) is the required treatment of cross-scale quadratic terms. If instead one already has only entrywise bounds |F_ab|,|F_ab'|<=e, their row-sum norm can be as large as n e, so n e<=d_n^2/100000 is a sufficient replacement. Keeping only each individual matrix element small and ignoring their sum is not sufficient.

[Proved, relative accuracy] Since the smallest endpoint offset is 2pi b r^(2n-1), the displacement-to-offset ratio under (16) is bounded by

    2eta/[2pi b r^(2n-1)]
       <= 2(1-2r)d_n/100000.

It is exponentially small in n for this designed target. This accuracy is achieved by imposing the correction budget. It is not a consequence of the arbitrary energy ratios t_a/t_b in (15).

## 6. Frozen and mobile conclusions

[Proved, frozen obstruction] For the frozen core the tree limit has coincident centers. The leading Gram hierarchy does not repair them. The loop and finite-express matrices determine their displacements through information not present in the Gram eigenvalues. C33's finite-n certificates and C37's algebraic certificates establish additional information at their stated sizes. Neither is an all-n consequence of the hierarchy, and C38 does not redo those certificates.

[Proved, frozen conditional result] Given an all-n construction of the target (12) satisfying (16), or any other all-n certified target with analogous root neighborhoods and margins, the quadratic part has the required nonresonance. One then has a nonzero absolute error budget at each fixed n. Finite-cluster higher-conversion and port remainders can be reduced below it by the ordered small-coupling choices in C32/C33. This does not establish their all-order thermodynamic continuation.

[Proved, audit of the crossing-only schedule] C32's sufficient schedule chooses eta_n=g^2 a_n and t_d<=eta_n m^2/(100000 n). In the quadratic symbol divided by g^2, the field is fixed at a_n/4, the finite-conversion remainder is bounded by a constant times g^2, and the loop is constrained to be O_n(g^2). If simple tree roots are available, their analytic center differences vanish at t_d=0 and are at most C_n |t_d| near zero. Thus a conservative bound on the finite-conversion correction and an upper bound on the pair-transfer margin can both be O_n(g^2). Decreasing g alone does not improve their ratio under that tied schedule. This is a limitation of those norm estimates, not proof that the actual finite-g corrections destroy the margin. Symmetry-correlated corrections could be much smaller in the relevant pair sum. A successful completion must prove that correlation or choose the loop at the leading-coefficient stage and then choose g using its certified margin. The crossing inequalities themselves have no lower bound on that margin.

[Assumed, interacting condition] Conversion of those local remainders to dressed Fermi-point and Gaussian-coupling bounds requires the local elimination and matching assumptions of C26/ir_control. For example a dressed shift less than d_n/16 preserves at least 3d_n/4 of the margin. Oscillatory quartics must stay weak until that momentum cutoff, and the stiffness metric must keep the higher-dimension operators irrelevant. These are additional requirements, not consequences of the finite root calculation. C40 is treating the velocity-dependent RG budget; no new RG theorem is claimed here.

[Proved, mixed arithmetic conditional on (12)] For the mobile chain take any length density

    1/5 <= nu <= 3/10,

for example the irrational value 1/4+sqrt(2)/100. At the free reference with exact length conservation and total payload U(1), every mixed non-Gaussian operator of dimension at most two is a length harmonic with |m|=1 times a payload bilinear. Its momentum divided by 2pi is +/-nu+x_i-x_j. The endpoint differences in (12) are within b of an integer or a half-integer. Therefore every mixed mismatch is at least

    2pi(1/5-b)=59pi/160.                                (19)

Pure length |m|=1 has mismatch at least 2pi/5. Terms with |m|>=2 or a length harmonic times a payload quartic have free dimension at least three (the pure |m|=2 dimension is four). They remain irrelevant if the full dimension metric is bounded below by (1-epsilon)I with epsilon<1/3. This uses the C20 compactification and dimension convention. Density-current terms remain in the Gaussian sector.

[Proved, mobile error accounting] C32 controls all four ordered length-pair bonds by

    s K2_(p,q)=K2_bare+E_(p,q),
    ||E_(p,q)|| <= E_fast(s/delta)
                   + x^2 zeta^2 L_bond/a_0^2,           (20)
    zeta=max_p w(k+p-1)/[U-w(k+p-1)].

Here a_0 is C32's positive inverse bound, not the mass a of Section 2. We retain its notation x,L_bond,delta by reference to C32 PROOFS (12)-(17). Both terms, and the scalar parts, are included. If p,q are encoded by length occupations N_j,N_(j+1), the residual is identically

    E_11+(E_21-E_11)N_j+(E_12-E_11)N_(j+1)
           +(E_22-E_21-E_12+E_11)N_j N_(j+1).            (21)

For max||E_pq||<=e the four coefficient norms are bounded by e,2e,2e,4e. Thus a sufficient triangle bound for their sum is 9e. Their scalar payload components are genuine length interactions and cannot all be removed by a constant shift. Formula (21) explicitly retains the length-payload cross terms.

[Proved, conditional parameter order] At fixed n, once a reference nonresonance and continuum tolerance is known, choose the right side of (20) below that tolerance divided by the finite extraction/matching constants and by 9. First choose w/U small but nonzero; next choose s/delta small but nonzero. Delta includes the actual charge tunneling gap of order U(w/U)^(2n+4), not U. Then choose marker hopping u proportional to lambda^2/s to set the length bandwidth, and finally choose lambda small. This is the existing C32 order. None of these steps provides a payload margin if it was zero before the wrapper was added.

[Assumed, mobile infrared scope] At finite wrapper ratios the residual operators in (21) couple gapless length and payload fields, so they are not one numerical Bloch correction to which (16) can be applied literally. Their continuum matching and dressed-root control are additional inputs. Conditional on that control, (19) and the payload margins protect the weakly coupled (n+1)-mode theory. The physical translation and winding-sector identification must also be the inherited C32/C8 gluing. No claim of a unique ground state on every ring or a thermodynamic central charge is proved by this note.

[Proved, final obstruction] The missing finite-dimensional statement for both actual cores is a lower bound on the nontrivial pair-sum functions of their full loop-dependent quadratic pencil, together with simple roots. The energy hierarchy supplies positive eigenvalues and permits small corrections after such a bound is known. It supplies neither the bound nor a sign/noncancellation theorem for (11). Therefore it does not complete the requested all-n frozen or mobile proof.
