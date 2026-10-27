P2 report. WZW levels and integrable nearest-neighbor chains

Closed-book investigation, October 21, 2026

No preexisting files were read or listed. No external references were consulted. The supplied background is treated as input. Algebraic statements proved below are distinguished from assumptions about infrared limits. No counterexample satisfying all the hypotheses was found, and neither proposed universal bound was proved.

The strongest results are (i) a complete answer when the invariant two-site algebra is spanned by the identity and the swap, (ii) a representation-theoretic bound covering level-matched WZW fusion families, (iii) a conditional theorem for ordinary Gaussian Bethe seas, (iv) an analytic integrability obstruction reducing invariant orthogonal and symplectic vector bonds to three branches, and (v) an exact exclusion of the mixed graded-permutation Hamiltonian at zero chemical potential. The unresolved steps are stated explicitly.

1. What the level bound would and would not establish

Write D = dim g, H = h^vee, and n = d-1 for a simple Lie algebra. The supplied WZW formula gives

    c(g_k) = kD/(k+H).

Consequently, if D > n,

    c(g_k) <= n  <=>  k(D-n) <= nH,
    k_max = floor[nH/(D-n)].                                      (1)

If D <= n, every positive level satisfies the central-charge inequality. In that case k_max should be interpreted as infinity; the set in the question has no finite maximum.

The sharp examples give the following values. SO(N) is taken with N >= 5 here to avoid low-rank normalization and semisimplicity exceptions.

| Lie algebra and on-site representation | d | k_max | embedding index supplied in the question |
| --- | ---: | ---: | ---: |
| su(2), spin 1/2 | 2 | 1 | 1 |
| su(2), spin 1 | 3 | 4 | 4 |
| so(N), vector | N | 2 | 2 |
| sp(4), fundamental | 4 | 1 | 1 |
| su(N), fundamental | N | 1 | 1 |
| su(N), antisymmetric square | N(N-1)/2 | N-2 | N-2 |
| e6, 27 | 27 | 6 | 6 |

For su(2) spin S >= 3/2, c(su(2)_k) < 3 <= 2S = d-1 for every finite k. A level bound there would be stronger than needed for a pure WZW theory.

There is an additional logical issue with spectator sectors. The implication

    k <= ind(R)  =>  c(g_k) <= d-1

bounds the charged WZW sector. It does not imply

    c_total = c(g_k) + c_spectators <= d-1

without another bound on the spectators. The GKO statement compares g_ind(R) with su(d)_1. It does not identify an arbitrary infrared spectator theory with the coset in that embedding.

A concrete consistency check is the d=4 ULS chain, restricted to its on-site Sp(4) subgroup. The given embedding construction fixes level 1, with

    c(sp(4)_1) = 10/4 = 5/2,
    c(su(4)_1) - c(sp(4)_1) = 3 - 5/2 = 1/2.

The full chain saturates c_total=3. Its charged WZW sector has c=5/2. This example uses only the supplied ULS identification and the supplied current embedding; no independent identification of the c=1/2 coset is needed. It illustrates why the charged-sector and total-central-charge questions must remain separate.

Global-form convention. Levels below use the simple Lie algebra with long roots of squared length 2. A spin-1 representation is faithful for SO(3), not SU(2). Referring to its su(2) current algebra avoids confusing a quotient of the symmetry group with a change in affine normalization.

2. A complete nearest-neighbor result

Theorem 1. Suppose the unitary on-site representation satisfies

    End_G(V tensor V) = span{I, P},

where P is the ordinary swap. Every nontrivial G-invariant nearest-neighbor Hamiltonian with a unitary compact CFT limit has c=d-1. Its microscopic G currents, when restricted from the ULS fixed point, have the embedding level ind(R).

Proof. First one can take the bond h itself to be G invariant. If only the full periodic Hamiltonian is stated to be invariant, replace h by its Haar average. The sum of the averaged bonds equals the original Hamiltonian, so possible telescoping choices of bond representative do not affect the argument.

The commutant hypothesis then gives

    h = a I + b P,     a,b real.

The additive constant a is irrelevant. If b>0, this is the ULS chain with a positive rescaling of energy. Its infrared theory su(d)_1 and c=d-1 are supplied facts. Restriction of its currents to G has the supplied embedding index.

If b<0, write b=-J with J>0. Since each P has eigenvalues +/-1, H >= -JL. A constant-color product state attains the bound. In the subspace with one site of a second color, a plane wave of momentum q has energy above that state

    epsilon(q) = 2J(1-cos q) = J q^2 + O(q^4).

Thus the low-energy spectrum contains a quadratic branch and does not have the assumed compact CFT scaling with positive velocity. If b=0, the Hamiltonian is flat and also lacks the stated CFT limit. This exhausts the cases. QED.

For a spin-1/2 on-site doublet,

    2 tensor 2 = 1 plus 3,

with each summand occurring once. The two projectors span the same algebra as I and P. Therefore this proves the requested k=1 result for all SU(2)-invariant nearest-neighbor spin-1/2 chains, using the Heisenberg/ULS infrared identification supplied in the question. It is not an anomaly argument.

The same proof applies to the SU(d) fundamental representation. It also applies to any other on-site group representation with the stated two-dimensional commutant.

The harder sharp cases have three-dimensional invariant bond algebras.

* Spin 1: h = a I + b X + c X^2, X = S_i dot S_{i+1}. The eigenvalues of X in total-spin 0,1,2 are -2,-1,1. Hence these three polynomials span all three projectors.
* SO(N) vector, N>=5: V tensor V splits into the singlet, antisymmetric tensor, and symmetric traceless tensor. Thus h = a I + b P + c K, with K = N |Omega><Omega| and |Omega> the normalized invariant singlet.
* Sp(2n) fundamental, n>=2: the decomposition is the singlet, antisymmetric traceless tensor, and symmetric tensor. Again h = a I + b P + c K, now with the singlet in the antisymmetric subspace.

Removing the additive constant and an overall positive scale leaves a one-parameter family in each case. This is a reduction of the problem, not a classification of its critical points. Counting relevant continuum couplings would only give a genericity argument. It cannot exclude an accidental higher-level point in such a family. No complete infrared classification of these families is proved here.

3. A dimension theorem for level-matched WZW sectors

The following result does not assume a Bethe ansatz. It is useful once a fusion/TBA analysis has identified the levels.

Let g be simple, and let an irreducible on-site representation have highest weight

    lambda = sum_i a_i omega_i,      s(lambda) = sum_i a_i.

For a reducible representation R = direct_sum_mu V_(lambda_mu), count multiplicities and set

    s(R) = sum_mu s(lambda_mu).

Theorem 2. Suppose a collection of WZW sectors based on g has positive integral levels k_1,...,k_p and

    K = sum_a k_a <= s(R).

Then

    sum_a c(g_(k_a)) <= dim R - 1.                               (2)

The theorem bounds the total only if these sectors account for the entire central charge. Unspecified neutral sectors are not covered. The level-budget hypothesis K<=s(R) is an assumption about the infrared theory, not a conclusion from nearest-neighbor locality or symmetry.

Proof, step 1: a root-string lower bound on dimension.

For every positive root alpha, the highest-weight vector generates an sl(2)_alpha string containing the weights

    lambda - j alpha,     1 <= j <= <lambda, alpha^vee>.

These are genuine nonzero weight spaces: the vector is a highest-weight vector for this sl(2), of nonnegative integral highest weight <lambda,alpha^vee>.

The weights listed for distinct positive roots are distinct. Equality j alpha = j' beta would make alpha and beta proportional. A reduced root system has no distinct positive proportional roots. Together with the highest weight itself, this proves

    dim V_lambda - 1
       >= sum_(alpha>0) <lambda,alpha^vee>
        = 2 <lambda,rho^vee>
        = sum_i C_i a_i,                                        (3)

where C_i = 2<omega_i,rho^vee>.

This counts weights, not multiplicities, so multiplicities can only strengthen the bound.

Proof, step 2: compare C_i with c(g_1).

For the Cartan convention A_ij=<alpha_i,alpha_j^vee>, write rho^vee=sum_j x_j alpha_j^vee. The defining equations <alpha_i,rho^vee>=1 give

    A x = 1,       C_i = 2 x_i.                                 (4)

Thus the required numbers follow from the finite Cartan matrices. The result is

| Type | minimum C_i | c(g_1) |
| --- | ---: | ---: |
| A_r | r | r |
| B_r, r>=2 | min(2r,r(r+1)/2) | r+1/2 |
| C_r, r>=2 | 2r-1 | r(2r+1)/(r+2) |
| D_r, r>=4 | 2r-2 | r |
| E6 | 16 | 6 |
| E7 | 27 | 7 |
| E8 | 58 | 8 |
| F4 | 16 | 26/5 |
| G2 | 6 | 14/5 |

For completeness, the classical solutions of (4) are

    A_r: C_i = i(r+1-i).
    B_r: C_i = i(2r+1-i), i<r;  C_r = r(r+1)/2.
    C_r: C_i = i(2r-i).
    D_r: C_i = i(2r-1-i), i<=r-2;
         C_(r-1) = C_r = r(r-1)/2.

Substitution into A C=2 verifies each expression. For the simply-laced exceptional diagrams, take a central node with arms of lengths (2,2,1), (2,3,1), and (2,4,1). The corresponding vectors C, ordered first by the central node and then along the arms, are

    E6: (42,30,16,30,16,22),
    E7: (96,66,34,75,52,27,49),
    E8: (270,182,92,220,168,114,58,136).

For F4 take the path Cartan matrix with A_23=-2 and A_32=-1, giving C=(22,42,30,16). For G2 take A=((2,-1),(-3,2)), giving C=(6,10). These are finite algebraic checks, not assumptions about a continuum limit.

Every row obeys min C_i >= c(g_1). For C_r the difference is

    (2r-1) - r(2r+1)/(r+2) = 2(r-1)/(r+2) >= 0.

The other inequalities follow directly from the table. Combining with (3) gives

    dim V_lambda - 1 >= s(lambda) c(g_1).                        (5)

For a sum of M irreducible summands,

    dim R - M >= s(R) c(g_1).                                   (6)

Proof, step 3: compare levels.

For each positive integer k,

    c(g_k) = k dim g/(k+h^vee)
           <= k dim g/(1+h^vee)
            = k c(g_1).

Therefore

    sum_a c(g_(k_a))
       <= K c(g_1)
       <= s(R) c(g_1)
       <= dim R - M
       <= dim R - 1.

This proves (2). QED.

A semisimple version follows by applying the root-string count to the union of the positive roots of all simple factors. Distinct factors change different components of the highest weight, so their counted weights remain distinct. One obtains the same conclusion if the level budget is satisfied separately for every factor and all central charge is included in those factors.

Theorem 2b. For every representation R of a simple Lie algebra in this normalization,

    ind(R) >= s(R).                                              (9a)

Consequently, the level-budget assumption in Theorem 2 also implies the proposed index bound for each of its WZW sectors, and for their diagonal current algebra.

Proof. Let theta be the highest long root and m=<lambda,theta^vee> for an irreducible summand. The coroot theta^vee has strictly positive integral coefficients in the simple-coroot basis, so

    m >= sum_i a_i = s(lambda).

Restrict the representation to the long-root su(2). This subalgebra has embedding index 1 in g with the long-root normalization. The highest-weight vector generates a spin m/2 multiplet. The index of a spin m/2 multiplet, normalized to 1 for a doublet, is obtained by summing the squared t^z eigenvalues,

    I_(m/2) = 2 sum_(j=-m/2)^(m/2) j^2
            = m(m+1)(m+2)/6.

The remaining su(2) multiplets contribute nonnegative indices. Restriction through an index-1 subalgebra preserves the index, so

    ind(V_lambda) >= m(m+1)(m+2)/6 >= m >= s(lambda).

Indices add under direct sums, proving (9a). The normalization of the long-root subgroup is essential; a short-root subgroup would carry an extra embedding factor. QED.

Thus, for the standard level-matched rectangular fusion identification k_IR=k and lambda=k omega_a, the actual stronger statement

    k_IR <= ind(R)

is established. The unproved part in extending this beyond that family remains the identification or upper bound of k_IR from the microscopic Hamiltonian. For the spin-S fusion family the index can be written explicitly as

    ind(R) = (2/3) S(S+1)(2S+1)
           = k(k+1)(k+2)/6,       k=2S,

which is at least k. At S=1 this gives k_IR=2 and ind(R)=4, consistent with the distinct spin-1 ULS point at level 4.

Application to rectangular SU(N) fusion representations. Let lambda=k omega_a. Its Young diagram is a rectangle of width k and height a. Equation (3) becomes

    d - 1 >= k a(N-a) >= k(N-1).

Consequently,

    c(su(N)_k)
      = k(N^2-1)/(k+N)
      <= k(N-1)
      <= d-1.                                                   (7)

Assumption used when applying this to an actual integrable chain: its full infrared theory is su(N)_k, or another collection obeying the level budget, with no additional central charge. The statement that the usual antiferromagnetic rectangular fusion chains have this infrared identification is a Bethe-ansatz/TBA input, not derived here from the Yang-Baxter equation alone. Once that input is supplied, (7) is a proof of the desired inequality for the whole family, including nonsymmetric rectangles.

For SU(2), k=2S and d=2S+1=k+1. The Takhtajan-Babujian identification gives

    c = 3k/(k+2),
    (d-1)-c = k(k-1)/(k+2) >= 0.                                (8)

The inequality is proved. The infrared identification is the additional input just stated.

This level budget is not a proposed universal microscopic upper bound. The spin-1 ULS chain has s(R)=2 but the principal su(2) embedding has level 4. The sufficient condition fails there while c=2 still saturates d-1. This check prevents confusing Theorem 2 with a solution of the level-bound problem.

4. Why the proposed WZW methods do not close the gap

4.1 The finite lattice charge algebra has no level

Let q_j^a be the on-site generators. They obey

    [q_j^a,q_l^b] = i f^(ab)_c delta_(jl) q_j^c.

For Fourier modes Q_n^a=sum_j exp(2 pi i n j/L) q_j^a, the algebra is

    [Q_n^a,Q_m^b] = i f^(ab)_c Q_(n+m)^c,

with momenta understood modulo L and no central term. This identity is independent of the Hamiltonian and contains no affine level.

At a nonchiral WZW fixed point the physical smooth charge density is J_L+J_R. The equal-time derivative-of-delta terms in the left and right current commutators have opposite signs and cancel in this sum. Hence the exact lattice algebra is compatible with every equal left/right level.

Projection onto a low-energy subspace does not preserve the algebra. If Pi is the projection, then

    [Pi A Pi, Pi B Pi]
      = Pi[A,B]Pi
        - Pi A(1-Pi)B Pi + Pi B(1-Pi)A Pi.                      (9)

The omitted states enter precisely through the extra terms. Chiral currents also require a dynamical separation into density and current. Neither operation is determined by the finite-dimensional trace ind(R). Thus no Schwinger-term proof of the upper bound follows from the on-site commutator.

4.2 Susceptibility and stiffness isolate the level, but do not bound it

Choose an SU(2) Cartan normalization with tr_fund(t^z t^z)=1/2. Define the zero-temperature stiffness D_K by

    E_0(Phi)-E_0(0) = D_K Phi^2/(2L) + o(Phi^2/L).

For one WZW sector with velocity v and no other gapless charged sector, the continuum response is

    chi = k/(2 pi v),     D_K = k v/(2 pi),
    2 pi sqrt(chi D_K) = k.                                    (10)

These formulas are continuum input with the normalization just specified. They are a possible way to measure k. They are not a microscopic upper bound.

For comparison, spreading a twist over the bonds gives

    H'(0)  = (i/L) sum_j [q_j,h_(j,j+1)],
    H''(0) = -(1/L^2) sum_j [q_j,[q_j,h_(j,j+1)]].

For a nondegenerate finite-size ground state, second-order perturbation theory gives

    E_0''(0) = <H''(0)>
               - 2 sum_(n>0) |<n|H'(0)|0>|^2/(E_n-E_0).

Therefore

    D_K <= -(1/L) sum_j <[q_j,[q_j,h_(j,j+1)]]>
         <= 4 ||q||^2 ||h||.                                  (11)

One can shift q by a scalar to optimize the last bound. Degeneracies and the order of thermodynamic limits require the usual care in defining the stiffness; they do not supply the missing susceptibility bound. No bound on chi D_K in terms of ind(R) follows from (11). In particular, a uniform finite-size susceptibility of a conserved total charge must not be computed by a naive nondegenerate Lehmann formula: its thermodynamic response involves level crossings and an order of limits.

If several charged sectors with velocities v_a contribute, the response instead gives

    chi D_K = (1/(4 pi^2)) (sum_a k_a/v_a)(sum_a k_a v_a)
            >= (sum_a k_a)^2/(4 pi^2),                          (12)

by Cauchy-Schwarz. The product is velocity independent only when there is a single charged velocity, or all charged velocities coincide.

A static sum rule has the same obstruction. With the normalization above, a WZW Cartan structure factor has

    S(q) = k |q|/(4 pi) + o(|q|).

But the on-site variance fixes the integral, not the slope,

    integral_(-pi)^pi dq/(2 pi) S(q) = Var(q_j).

If the linear law holds approximately up to q_*, its contribution is approximately k q_*^2/(8 pi^2). A bound on the variance therefore bounds k times the square of the scaling window. There is no supplied lower bound on q_*. This is the momentum-space version of the conformal-window obstruction in the background.

4.3 Anomalies give no upper bound on a nonchiral level

For the microscopic vector G symmetry, the perturbative anomaly is proportional to k_L-k_R. It vanishes at every nonchiral level. Mixed translation/internal-symmetry anomalies can constrain a residue class of the level after specifying the continuum action of translation. Such a congruence leaves arbitrarily large admissible levels. No upper bound follows from anomaly matching alone.

4.4 Partons supply a construction, not a universal classification

The supplied d-parton construction produces the embedding g_ind(R) inside su(d)_1. To turn it into a universal level bound would require proving that every allowed critical chain is obtained from that construction in a way that cannot increase the surviving G level. No such statement follows from an on-site Hilbert-space identification. Gutzwiller projection is a constrained, normalized projection of a state, not an argument that every possible critical state is a conformal subtheory of that particular free-parton theory.

4.5 An exact test showing why nearest-neighbor locality must enter

For a spin-1/2 chain, consider

    H_r = J sum_j S_j dot S_(j+r),      J>0,

with L divisible by r. This is a range-r chain, so r>1 lies outside the question. It decomposes into r independent antiferromagnetic Heisenberg chains. Using the supplied Heisenberg infrared identification, the total central charge is r. The microscopic uniform SU(2) current is the sum of the r level-1 currents, so its affine level is

    k_diag = r,        ind(on-site doublet) = 1.

All on-site charge-algebra and representation data are unchanged as r grows. A proof using only those data would therefore give a false result. The precise nearest-neighbor restriction has to enter somewhere substantive.

The diagonal current algebra and its coset need not form an unrestricted tensor-product modular invariant; the affine-level statement does not require that stronger claim. For r=2, the diagonal su(2)_2 stress tensor has c=3/2 and the commuting coset has c=1/2.

Blocking r sites gives d_block=2^r and embedding index r 2^(r-1). The same example then obeys both the index inequality and r<=2^r-1. Thus it is consistent with the blocked nearest-neighbor conjecture.

A separate false positive must also be avoided. The d=3 ULS chain admits both the regular su(2) embedding, with on-site representation 2 plus 1 and index 1, and the principal embedding, with on-site spin 1 and index 4. The level-4 current algebra belongs to the second microscopic embedding. It is not a violation for the first representation. The G symmetry used to define the infrared level must be the same one used to define R.

5. Integrable chains and what can be proved from Bethe-sea counting

Theorem 3, conditional Gaussian-sea version. Suppose a critical chain has r independent ordinary Gaussian modes, each arising from a single occupied Bethe interval with two simple Fermi edges, with no further gapless auxiliary/string/neutral sector. Suppose also that these modes are indexed by independent on-site Cartan charges of a faithful compact semisimple representation on C^d. Then

    c_total = r <= d-1.                                        (13)

Proof of the thermodynamic part. Under the stated Gaussian and regular-dressing assumptions, a Fermi edge at lambda=B has dressed excitation energy epsilon'(B)(lambda-B) and dressed momentum derivative p'(B). The density of states per unit length and rapidity is p'(B)/(2 pi), and

    v = |epsilon'(B)|/p'(B).

The thermal edge contribution relative to zero temperature is

    -T integral d lambda [p'(B)/(2 pi)]
          log(1 + exp(-|epsilon'(B)(lambda-B)|/T))
      = -pi T^2/(12 v),

using integral_0^infinity dx log(1+exp(-x))=pi^2/12. Two edges contribute -pi T^2/(6v). In a Gaussian finite-size spectrum these are the left and right parts of one compact boson, with c=1. Summing the independent modes gives c_total=r. If the two edge velocities differ, they appear separately in the thermal coefficients; the Gaussian mode still has one left and one right chiral central charge. A common relativistic velocity is not needed for the count across decoupled Gaussian sectors.

This step assumes that the quasiparticle dressing and auxiliary sectors introduce no additional infrared entropy. Without that assumption, the edge integral alone is not a complete TBA calculation. Also, a thermal coefficient alone does not identify c unless the velocities and the mode interpretation have been fixed.

Proof of the representation-theoretic part. The Cartan generators are commuting Hermitian matrices and can be simultaneously diagonalized. They are traceless because a representation of a semisimple Lie algebra has tr R([X,Y])=0 and the algebra is equal to its commutator algebra. Faithfulness makes the Cartan map injective. The real vector space of traceless diagonal d by d Hermitian matrices has dimension d-1. Hence

    rank g <= d-1.

With one Gaussian mode per independent Cartan direction, r<=rank g<=d-1. QED.

This proves the desired bound for the specified Gaussian-sea class. It does not prove that every regular Hermitian R-matrix chain belongs to that class.

Why nesting rank is insufficient. For a level-k SU(2) WZW theory,

    c = 3k/(k+2) = 1 + 2(k-1)/(k+2).                           (14)

The rank-one Cartan boson accounts for the first term. The remaining parafermion central charge is nonzero for k>1. Taking the usual spin-S fusion identification k=2S as input, the spin-1 Takhtajan-Babujian chain already has c=3/2 although su(2) has rank one. The spin-1 Hamiltonian can be written, up to scale and a constant, as X-X^2, so it is nearest neighbor on d=3. This example does not violate c<=d-1. It invalidates the premise that a nesting level contributes at most one unit of c.

Dilogarithm sums make the same issue visible. A TBA central charge is a sum or difference of Rogers dilogarithm terms evaluated at constant Y-system solutions. Each normalized term lies between zero and one, but the number of string and auxiliary nodes is not bounded by d-1 merely because the local R matrix acts on C^d tensor C^d. Even rank-one models can have arbitrarily long strings. Cancellation and truncation identities depend on the particular TBA system. The Yang-Baxter equation as stated supplies neither a finite-node bound nor the identity needed to compare the complete sum with d-1.

Theorem 2 replaces the false nesting count by a valid dimension comparison for the usual level-matched WZW fusion families. Its remaining model-dependent input is the infrared level identification and the absence of uncounted sectors.

A closure property gives another rigorous subclass. Suppose independent regular R-matrix chains act on local spaces of dimensions d_1,...,d_p and each satisfies c_a<=d_a-1. On the local tensor product, the product of their commuting R matrices is again regular and solves the Yang-Baxter equation, and its logarithmic derivative is the sum of the factor Hamiltonians. The full central charge is sum_a c_a, while

    sum_a (d_a-1) <= product_a d_a - 1.

For two factors the difference is (d_1-1)(d_2-1)>=0, and induction proves the general formula. Therefore this operation cannot generate a violation from factors already known to satisfy the bound. Different positive energy rescalings, and thus different velocities, do not change this conclusion.

6. An analytic integrability obstruction for invariant vector bonds

There is a further restriction in the intersection of the WZW test classes and the integrable class. This restriction is derived from the Yang-Baxter equation itself.

Theorem 4. Suppose check R(u) is analytic near u=0, check R(0)=I, and obeys the additive braid Yang-Baxter equation

    R_12(u) R_23(u+v) R_12(v)
      = R_23(v) R_12(u+v) R_23(u).                              (15)

Suppose its logarithmic derivative, after a scalar shift, is the invariant bond

    h = b P + c K,

with K=d|Omega><Omega| the positive singlet contraction. Then the following necessary conditions hold:

    orthogonal vector, d>=3:    bc[2b+(d-2)c] = 0;
    symplectic vector, even d>=4: bc[2b-(d+2)c] = 0.             (16)

Here the orthogonal statement is for bonds of the displayed form; in the SO(4) case this need not exhaust all SO(4)-invariant bonds. The symplectic notation is Sp(d) with fundamental dimension d.

Assumptions and scope. Analyticity and the one-parameter additive equation are used. General non-difference two-parameter R matrices are not covered. The bond itself is assumed invariant. Haar averaging a merely globally invariant Hamiltonian does not by itself prove that the averaged bond is still the logarithmic derivative of such an R matrix.

Proof, step 1: obtain the third-order necessary identity.

Put v=-u in (15). The resulting equality is

    [R(u)R(-u)]_12 = [R(-u)R(u)]_23.

An operator supported on sites 1,2 that is also supported on sites 2,3 must be supported on site 2. Applying the same equation with u replaced by -u shows that R(u)R(-u) is both I tensor M and N tensor I. Thus it is scalar. Near u=0 one can divide R(u) by an even analytic scalar so that

    R(u)R(-u)=I.

This scalar normalization does not change the nonconstant part of h. The Taylor expansion then has the form

    R(u)=I+u h+(u^2/2)h^2+(u^3/6)r+O(u^4).

Set A=h_12 and B=h_23. Twice the coefficient of u^2 v on the left side of (15) is

    r_23 + A^2 B + A^3 + 2 A B^2 + B^2 A + 2 A B A.

Twice the corresponding coefficient on the right is

    r_12 + B A^2 + B^3 + 2 A^2 B + A B^2 + 2 B A B.

Subtracting and rearranging gives the necessary identity

    [A+B,[A,B]] = Q_23-Q_12,      Q=r-h^3.                      (17)

No assumption about the symmetry of R or Q is needed for this identity.

Proof, step 2: a matrix element that annihilates every right-hand side of (17).

Write

    s=P_12, t=P_23, e=K_12, f=K_23.

Expanding the left side as a polynomial in b,c gives

    [A+B,[A,B]] = b^3 C0 + b^2 c C1 + bc^2 C2 + c^3 C3,

where

    C0=[s+t,[s,t]]=2(t-s),
    C1=[e+f,[s,t]]+[s+t,[e,t]+[s,f]],
    C2=[s+t,[e,f]]+[e+f,[e,t]+[s,f]],
    C3=[e+f,[e,f]]=2(f-e).

The first equality uses the permutation braid relations. For the last, use e^2=de, f^2=df, efe=e, and fef=f. These singlet-contraction relations follow directly from the normalized maximally entangled vector. For example, expanding the double commutator cancels the terms d(ef+fe), leaving 2(f-e).

For the orthogonal case use basis labels 0,1,2,...,d-1 and

    K|a,b> = delta_ab sum_l |l,l>.

Take the matrix element from |2,1,2> to <0,0,1|. Since both the first and third labels change, the matrix element of Q_23-Q_12 vanishes for every two-site Q. The matrix elements of C0 and C3 also vanish. Direct application of P and K gives

    <[e+f,[s,t]]> = 1,
    <[s+t,[e,t]+[s,f]]> = 1,
    <[s+t,[e,f]]> = -1,
    <[e+f,[e,t]+[s,f]]> = d-1.

The brackets in this display denote the selected matrix element. To make the contractions explicit, the only nonzero word contributions in these four expanded commutators are respectively

    est: +1;
    set: +1;
    -efs: -1;
    eet-efs: d-1.

For example, t|2,1,2>=|2,2,1>, e sends this to sum_l|l,l,1>, and a second e multiplies by d. These facts give <C1>=2 and <C2>=d-2. Equation (17) therefore implies the first polynomial in (16).

For the symplectic case write d=2n and use the real antisymmetric form

    J_(a,a+n)=1, J_(a+n,a)=-1,      0<=a<n,
    K|a,b> = J_ab sum_(i,j) J_ij |i,j>.

Take the matrix element from |1,0,n+1> to <0,0,n|. Again both outer labels change, so the right side of (17), C0, and C3 have zero matrix element. In the same four groups as above, the nonzero word contributions are

    -fts: +1;
    -2sfs+2tet-tfs: -2+2+1 = 1;
    fet: -1;
    fet-ffs: -1-d.

These follow from the singlet vector sum_l(|l,l+n>-|l+n,l>). For example, ts sends the ket to |0,n+1,1>, so fts has matrix element -1. Also f^2= d f gives <ffs>=d. Thus <C1>=2 and <C2>=-d-2, proving the second polynomial in (16). QED.

The possibilities are therefore

    orthogonal:  c=0, or b=0, or c=-2b/(d-2);
    symplectic:  c=0, or b=0, or c= 2b/(d+2).

In particular, the spin-1 invariant bond has only c=0, b=0, or c=-2b under these analytic integrability assumptions. The last is the Takhtajan-Babujian bond up to a constant and scale. The Sp(4) fundamental has only c=0, b=0, or c=b/3. The b=0 branch is the singlet Temperley-Lieb bond. These are necessary restrictions, not a proof that every sign and scale on every branch has a compact CFT limit.

Some signs can be discarded without a Bethe ansatz. A positive multiple of K is frustration free: for the orthogonal case take a constant complex null vector v with sum_a v_a^2=0, and for the symplectic case take any constant vector v. Then K|v,v>=0. The negative-b rational branches are also frustration free after adding a constant. Their bond eigenvalue minimum is attained by symmetric traceless products in the orthogonal case, or symmetric products in the symplectic case. By the supplied frustration-free fact, if these branches are gapless they do not have the stated compact-CFT scaling. This leaves the antiferromagnetic permutation branch, the positive-b rational branch, and the negative-K branch as the remaining candidates in these invariant analytic families. The negative-K branch's gap and the general rational branch's full infrared identification are not derived here.

The accompanying integrability_check.py verifies the full cubic coefficient matrices with exact integer arithmetic for orthogonal dimensions 3 through 6 and symplectic dimensions 4 and 6. Those checks support, but are not needed in place of, the general selected-matrix-element proof above.

7. A graded integrable family that can be excluded exactly

The following is a full lattice calculation, with no infrared assumption.

Give the d=m+n basis states parities p(a) in {0,1}, with m even and n odd states, and define

    P_g |a,b> = (-1)^(p(a)p(b)) |b,a>.

Then P_g is a Hermitian involution. On neighboring factors it obeys the ordinary braid relation: both sides reverse three labels and accumulate the same sum of pairwise parity products. It follows by expanding the additive braid Yang-Baxter equation that

    check R(u) = I + u P_g

is a regular solution. For clarity, with s=P_(g,12), t=P_(g,23), the two sides of

    (1+u s)(1+(u+v)t)(1+v s)
      = (1+v t)(1+(u+v)s)(1+u t)

both expand to

    (1+uv) I + (u+v)(s+t)
      + u(u+v) st + v(u+v) ts
      + uv(u+v) sts,

after using s^2=t^2=I and tst=sts. Its logarithmic derivative at u=0 is h=P_g.

Theorem 5. If m,n are both nonzero, neither sign of the zero-field graded-permutation Hamiltonian has the full low-energy limit of a unitary compact CFT with positive velocities.

Proof for H=J sum P_g with J>0. Each bond is bounded below by -J. Choose a fixed odd basis state f. The state |f> tensor ... tensor |f> has energy -JL and saturates the lower bound.

Now replace f at one site j by a fixed even state b, and denote this state by |j>. For L>=3 the two bonds touching the defect swap it to its neighbors with positive sign; all other bonds contribute -1. Therefore

    (H+JL)|j> = J(2|j> + |j-1> + |j+1>).

A plane wave has

    epsilon(q) = 2J(1+cos q).

For even L the minimum is at q=pi, and near it

    epsilon(pi+p) = J p^2 + O(p^4).

In particular the first nonzero branch energy is O(L^-2), incompatible with a compact CFT spectrum with fixed positive velocities. This is an actual branch of the full Hamiltonian, not a variational estimate.

For H=-J sum P_g, use a constant even background, which again has energy -JL and saturates the lower bound. A defect of a different parity hops with amplitude -J and has

    epsilon(q) = 2J(1-cos q),

with the same quadratic conclusion. QED.

This result applies to the ordinary tensor-product Hamiltonian whose local matrix is the signed swap just defined. It does not assume that every convention for a graded periodic transfer matrix has the same boundary implementation.

Scope. Chemical potentials, fixed-density restrictions, trigonometric deformations, and other superalgebra R matrices are not covered. A chemical potential that favors finite particle density can remove the polarized vacuum from the ground-state sector. The supplied supersymmetric t-J example belongs to the broader discussion of graded chains and cannot be excluded merely by the words "graded" or "superalgebra." A positive-definite Hermitian lattice model need not have a nonunitary continuum merely because its algebraic organization uses a superalgebra.

8. Candidate families and the remaining work

| Family | What is established here | What remains unproved |
| --- | --- | --- |
| Fundamental SU(d) permutation chain | Complete nearest-neighbor answer by Theorem 1; c=d-1 in the AF case | Nothing for this family beyond the supplied ULS infrared fact |
| Spin-1/2 SU(2), arbitrary symmetric NN bond | All such bonds reduce to the preceding case; k=1 | Nothing beyond the supplied Heisenberg infrared fact |
| Spin-1 SU(2), SO(N) vector, Sp(4) fundamental | Invariant bond space reduces to three projectors; analytic additive R-matrix derivatives reduce further to three branches; thresholds are k=4,2,1 | Exclusion of every higher-level critical point across the general one-parameter bond family, and full IR identification of the remaining integrable branches |
| Rectangular SU(N) fusion chains | Algebraic proof c<=d-1 once the standard full su(N)_k infrared identification is supplied | A derivation of that identification is not included here; arbitrary regular R matrices are not classified |
| Other level-matched simple-group fusion families | Theorem 2 applies whenever the full level budget is proved | Identification of levels, additional sectors, and all massless modes |
| Ordinary Gaussian nested Bethe seas | Theorem 3 proves the bound with one charge mode per Cartan direction | Deriving those hypotheses from an arbitrary regular R matrix |
| Rational mixed gl(m|n) signed-swap model, zero field | Theorem 5 exhibits an exact quadratic branch for either sign | Finite-density, deformed, and more general graded models |
| osp(m|2n) and other superalgebra R matrices | No universal rank or grading shortcut is valid | Hermiticity, compactness, full TBA entropy, and central charge must be checked case by case |
| Twisted a_2^(2), d=3 | Rank-one counting alone is insufficient | No classification or full infrared derivation is given here |

For twisted or graded models, compactness must be checked independently of Hermiticity. A useful conditional diagnostic is the noncompact cigar formula

    c = 2 + 6/(kappa-2),
    h_min = 1/[4(kappa-2)],
    c_eff = c-24 h_min = 2,        kappa>2.

If a proposed d=3 continuum limit has this spectrum, its formal c>2 does not give a counterexample under the compactness assumption. The spectrum is noncompact and the normalizable ground-state finite-size coefficient involves c_eff. The association of any particular a_2^(2) regime with this continuum, and Hermiticity of its Hamiltonian in that regime, are not established in this report. These formulas are a conditional diagnostic, not an asserted identification of an explicit lattice counterexample.

The decisive missing statements are now localized.

(A) One needs a nearest-neighbor-sensitive inequality on the infrared affine level, or a complete classification of the small invariant bond families. On-site traces, vector anomalies, current commutators, and response sum rules do not provide it. If neutral sectors are allowed, a separate bound on their central charge is also required.

(B) One needs a theorem connecting an arbitrary finite-dimensional regular Hermitian R matrix to its complete massless TBA/Y-system data and bounding the resulting dilogarithm sum. Counting Cartan generators or nesting rank proves only the Gaussian subclass. Fusion already shows why the general inference fails.

(C) No explicit nearest-neighbor violation of c_total<=d-1 or of k<=ind(R) was obtained. The range-r spin chain is a counterexample only to an on-site-only level argument, and it was not treated as a counterexample to the stated problem.

9. Checks and limits of this investigation

The accompanying checks.py was created in this session and run locally. It uses exact rational arithmetic for the Cartan-matrix calculations. It verified all classical min-C formulas for ranks 2 through 12, the exceptional vectors displayed above, the elementary SU(N) level/dimension inequality for 2<=N<=29 and 1<=k<=29, and that sp(2n)_2 exceeds d-1 for 2<=n<=29. The analytic arguments establish these inequalities beyond the tested ranges.

It also checked the exact mixed graded-permutation one-defect dispersion for periodic lengths 4 through 9 and every allowed momentum. All checks passed. Runtime was below one second. No many-body ground-state extrapolation or central-charge fit was performed. Such a fit would not certify the absence of a higher-level point, and no specific candidate violation was found that would justify using the requested ED budget.

The pure spin-1/2 class is settled under the supplied Heisenberg input. The level-matched fusion and Gaussian-sea classes obey the bound under their stated infrared assumptions. The universal WZW level bound and the universal regular-R-matrix bound remain open in this investigation.
