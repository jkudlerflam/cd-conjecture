# C29 derivation

[proved, scope] The main obstruction below concerns balanced charge-conserving rails, including arbitrary position-dependent Peierls phases. The partial repair uses two such buses and works for the leading conversion coefficient at n=2. It is not an arbitrary-n construction and it is not an exactly Gaussian order-lambda-squared bond at fixed conversion strength.

## 1. Nonuniform flux in one C22 bus

[assumed, inherited definitions] Use C22 RESULT sections 2 to 4 for its arrow grammar, charge convention, neutral Z-string gauge, and vacuum walker. The hypotheses retained here are an open path with vertices 0,...,n, data on vertices 0,...,n-1, three excited rails with identical scalar diagonal d and hopping magnitude t, balanced local conversions g c_i and g c_i^dagger, and the same three port channels. Only their link phases are changed. No uniformity of the phases is assumed in this section.

[proved] Let phi_(a,l) be the phase for a hop l to l+1 on rail a=0,+,-. Define

    K_i = sum_(l<i) [phi_(0,l)-(phi_(+,l)+phi_(-,l))/2],
    X_i = sum_(l<i) [phi_(-,l)-phi_(+,l)]/2,
    K_0=X_0=0.

Each scalar rail is gauge equivalent to A=d I-t A_path. If p_(a,i)=sum_(l<i) phi_(a,l), its heat kernel is

    [exp(-tau A_a)]_(s,i)
      = exp[i(p_(a,s)-p_(a,i))] [exp(-tau A)]_(s,i).

This follows by conjugating A by diag(exp(i p_(a,i))). It uses the path, not translation invariance.

[proved] Write k_si(tau)=[exp(-tau A)]_(s,i) and

    f_s(i,tau)=integral_0^tau k_si(tau-u) k_is(u) du,
    F_ij(d)=integral_0^infinity f_n(i,tau) f_0(j,tau) d tau.

For d>2t cos(pi/(n+2)), the integrals converge. The entrywise expansion of exp(t tau A_path) makes f_s(i,tau)>0 for tau>0 when t>0. Thus every F_ij is positive.

[proved] The first Duhamel conversion in the 0-to-+ source return has coefficient -g f_s(i,tau) times the relative rail gauge and c_i or c_i^dagger as set by the source orientation. The negative-charge channel gives its charge-conjugate operator. In the two-cell bond, one conversion in each cell gives the two phases

    exp[i(X_n-X_i+X_j)] exp[+i(K_n-K_i-K_j)],
    exp[i(X_n-X_i+X_j)] exp[-i(K_n-K_i-K_j)].

Their relative sign is minus. Indeed the right neutral source is |n,0> P, the charged right sources are |n,+> and |n,->, and P c_i=-c_i P. Subtracting the two phases therefore gives

    T_ij = 2i exp(i phi) F_ij(d)
           exp[i(X_n-X_i+X_j)] sin(K_n-K_i-K_j).       (1)

Here phi is the inter-cell charged-pair phase. The logical term is

    lambda^2 alpha_n^4 g^2
       sum_ij [T_ij c_i^dagger P tensor c_j + h.c.].   (2)

The convention in the scripts maps this directly to global Jordan-Wigner fermion hopping. At uniform phases, (1) reduces to C22 equation (16).

[proved] Terms with both conversions in one cell are scalar at this order. For each ordered pair i,j, exchange the two conversion times and the initial neutral/charged source species. Equality of the scalar rail kernels and their diagonal gauges makes the two scalar coefficients equal after reversing the order. The two logical factors are c_i^dagger c_j and c_j c_i^dagger. Their sum is delta_ij I. The negative-charge channel obeys the same identity. Consequently there is no induced quadratic chemical potential in the balanced equal-rail family. This argument allows arbitrary link phases but not unequal rail scalar Hamiltonians or unbalanced conversions.

[proved, independent evaluation] With H_pair=A_pair+g C_pair, and B the three pair-source maps, the same coefficient is

    [g^2](-B^dagger H_pair^-1 B)
       = -B^dagger A_pair^-1 C_pair A_pair^-1
                         C_pair A_pair^-1 B.           (3)

Expanding the inverse by its convergent Neumann series proves (3); analytic continuation within the positive interval gives the same Taylor coefficient. The scripts evaluate (3) and (1) independently.

## 2. A flux-only normality obstruction

[proved] At n=2 write T=[[a,b],[c,d']]. The two off-diagonal sines in (1) are both sin(K_2-K_1)=sin(kappa_1). The factors involving X have unit modulus. Hence

    |b|=2 F_01 |sin(kappa_1)|,
    |c|=2 F_10 |sin(kappa_1)|.                          (4)

For the three-vertex scalar path, reflection gives

    F_10-F_01
      = integral f_0(1,tau)[f_0(0,tau)-f_0(2,tau)] d tau >0. (5)

To prove strictness directly, put a0=sqrt(2)t. The required heat entries are

    k_00(tau)=exp(-d tau)[1+cosh(a0 tau)]/2,
    k_02(tau)=exp(-d tau)[cosh(a0 tau)-1]/2.

The first is strictly larger than the second. Both are nonnegative, so their self-convolutions obey f_0(0,tau)>f_0(2,tau). The other factor in (5) is positive. This proof works throughout the positive path spectrum, d>sqrt(2)t.

[proved] The (0,0) entry of [T,T^dagger] is |b|^2-|c|^2. Normality therefore forces sin(kappa_1)=0 and b=c=0. If kappa_1=m pi, then

    sin(K_2)=(-1)^m sin(K_1),
    sin(K_2-2K_1)=-(-1)^m sin(K_1).

The diagonal X factors in (1) both equal exp(i X_2). Thus both nonzero diagonal hoppings have the same phase modulo pi. They cannot furnish two independently displaced flavor centers. For zero sine there is no hopping at all.

[proved, scope] This excludes any mechanism whose only effect is to prescribe arbitrary phases on these same balanced rails. A frozen classical density track that only changes these phases also falls within the theorem, even if it supplies a different phase on every link. It does not exclude a fluctuating second species, unequal scalar rail dynamics, more bus types, different graph geometry, or an entangled code. It is a restricted spatial-resolvent obstruction, not an obstruction to computation-free cells in general, and not a no-go for all fixed-alphabet Hamiltonians.

## 3. Two virtual bus types on 26 letters

[proved, definition of the microscopic model] Use the C22 symbols but let the token species be

    s in {v,A0,A+,A-,B0,B+,B-}.

There are L_a,R_a (4 letters), X_(s,a) (14 letters), M_empty and M_s (8 letters). Thus d0=26. The arrows are unchanged. The allowed words in a cell containing n data sites are L^j X_s R^(n-j-1) M_empty and L^n M_s. There is one token, not an optional truncation to one token. The complete legal dimension is D_n=7(n+1)2^n.

[proved, rule table] All forward moves have amplitude -t times a phase,

    X_(s,a) R_b -> L_a X_(s,b),
    X_(s,a) M_empty -> L_a M_s,

plus their Hermitian conjugates. Multiply neutral excited moves (A0,B0) by (-1)^a. All vacuum moves and all charged excited moves have phase zero. At n=2 the neutral bulk-to-bulk move has phase pi/2, and the neutral bulk-to-marker move has phase zero. This is one repeated bond that distinguishes a marker locally; it is not a position-indexed Hamiltonian. For the separate n=1 control, the terminal neutral phase is pi/3. Couplings may depend on n, as in route A.

[proved, remaining terms] For each r=A,B include

    g |X_(r+,0)><X_(r0,1)|
      + g |X_(r-,1)><X_(r0,0)| + h.c.

Give each excited type r the onsite energy Delta_r=d_r-2t cos(pi/(n+2)), with d_A=3, d_B=5, t=1. Give every marker the C22 selector chemical potential and put penalty K=8 on forbidden arrow bonds. Onsite terms are split equally across adjacent bonds. All moves preserve external arrows, so the forbidden-bond projectors commute with the Hamiltonian. All terms preserve data occupation plus bus charge, with charges 0,+1,-1 for r0,r+,r-. Exact on-site U(1) and parity follow.

[proved] The vacuum walker is unchanged. Its wavefunction is

    v_j=sqrt(2/(n+2)) sin(pi(j+1)/(n+2)),
    G psi=sum_j v_j |j,v> psi,
    alpha_n=v_0=v_n.

It gives an exact 2^n code. The excited types are orthogonal blocks. Each has a gap at least Delta_r-|g| above the vacuum ground, because its hopping is a direct sum of unitary-gauged paths and the local conversion norm is |g|. The gap of the whole cell is at least the minimum of these bounds and the vacuum-walker gap. Adding the second type does not change the selector energy or its zero-coupling length choice. The C22 selector and defect argument applies with min_r Delta_r>|g|. This statement inherits that frozen-selector argument; it is not a new mobile-selector proof.

[proved, ports] For each r define

    F_(r,0)=|r0><v|+|v><r0|,
    F_(r,+)=|r+><v|+|v><r-|,
    F_(r,-)=F_(r,+)^dagger.

Use these on both marker and bulk token letters, acting trivially on the data bit. On an oriented marker-to-bulk bond put

    V=sum_r w_r [F_(r,0)^R F_(r,0)^L
           +exp(i phi_r) F_(r,+)^R F_(r,-)^L
           +exp(-i phi_r) F_(r,-)^R F_(r,+)^L],         (6)

with (w_A,w_B)=(1,8) and (phi_A,phi_B)=(pi/2,0). Every port takes v to an excited species, so G^dagger F G=0.

[proved, optional mass probe used in ED] Include lambda eta E on every site, where E is the projector onto an excited token. This annihilates the code and is local. The final ED uses eta=10. It changes neither PVP=0 nor the second-order coefficient. It exposes a nonzero cubic remainder; without it the open-chain model is even in lambda. This term is explicitly part of the tested Hamiltonian.

## 4. Exact second order and its controlled conversion expansion

[proved] In the neutral Z-string gauge, H_r has scalar path hopping on each of its three rails and conversions g c_i, g c_i^dagger. Sources divided by alpha_n are

    S_(r,L,a)=|0,ra> I,
    S_(r,R,0)=|n,r0> P,
    S_(r,R,+/-)=|n,r+/- > I.

Let B_r=S_(r,R,0) tensor S_(r,L,0)
       +exp(i phi_r) S_(r,R,+) tensor S_(r,L,-)
       +exp(-i phi_r) S_(r,R,-) tensor S_(r,L,+).

The complete second-order bond, with g fixed, is

    K2(g)=-alpha_n^4 sum_r w_r^2 B_r^dagger
                   (H_r tensor I+I tensor H_r)^-1 B_r. (7)

The physical contribution is lambda^2 K2(g). The inverse is positive. Each B_r has norm sqrt(3). Cross terms between A and B vanish because H0 preserves bus type and their excited sectors are orthogonal. Cross terms between distinct spatial bonds vanish because one outer cell is left excited after two different bond actions. The mass probe annihilates P and has no P-to-Q matrix elements, so it contributes nothing to (7).

[proved] Expanding (7) in g gives a scalar at g=0, followed by the quadratic bond with

    T=sum_r w_r^2 T^(r),                               (8)

where each T^(r) is (1). Odd powers of g vanish. The operator assigning a minus sign to every neutral excited rail changes C_pair to -C_pair and leaves the pair-source space invariant with positive sign. Thus all odd matrix elements of the inverse series between B_r sources vanish.

[proved, explicit error bound] Put delta_r=d_r-2t cos(pi/(n+2))>0. The pair A_r has smallest eigenvalue 2 delta_r and ||C_pair||<=2. For |g|<delta_r, the sum after g^2 obeys

    ||K2(g)-K2(0)-alpha_n^4 g^2 K_(2,2)||
      <= alpha_n^4 sum_r [3 w_r^2/(2 delta_r)]
                 (|g|/delta_r)^4/[1-(g/delta_r)^2].    (9)

Here K_(2,2) is the full g^2 Taylor coefficient divided by alpha_n^4, including its scalar same-cell contribution. This is a conservative absolute bound, not the measured size of the nonquadratic part. It follows by summing the even Neumann terms. For a fixed two-flavor reference with positive hopping and momentum margins, (9) divided by its g^2 logical scale goes to zero as g goes to zero. There is no uniform-in-n assertion.

[proved, finite-cluster cubic coefficient] On an open M-cell chain let H=H0+lambda(V+eta sum E_x)+lambda^2 W, where W preserves the code exactly as in the next section. Let R=Q H0^-1 Q. Since P(V+eta E)P=0, the canonical effective coefficients are

    H_eff=lambda^2(PWP-PVRVP)
                 +lambda^3 PV R (V+eta E) R VP+O(lambda^4).

The mixed terms with W vanish because QWP=0. The pure three-V term vanishes: the product of excited-token signs on alternating cells sends V to -V and is identity on the code. Between two bond creations E counts two excited cells. Therefore the cubic coefficient is the sum of bond terms

    K3=2 eta alpha_n^4 sum_r w_r^2 B_r^dagger
                        (H_r tensor I+I tensor H_r)^-2 B_r. (10)

The same outer-cell argument excludes mixed spatial bonds. It is positive. The scripts compare with both (7) and (10).

[proved, finite remainder] For each fixed cluster the code is an isolated finite-dimensional spectral subspace. On the circle |z|=gap(H0)/2, the resolvent Neumann series converges if |lambda| ||V+eta E||+lambda^2 ||W||<gap(H0)/2. Its contour integral yields an analytic spectral projector. The overlap with P is invertible near zero, so its polar isometry and the canonical effective matrix are analytic. The Taylor expansion above then has a bounded O(lambda^4) remainder in that neighborhood. This justifies the O(lambda^3) error after second order. The radius and constants depend on the cluster; no thermodynamic band argument is claimed.

## 5. Normal exchange with two independent phases at n=2

[proved] For both buses set kappa_0=pi/2 and kappa_1=0, with chi=0. Then K_2=K_1=pi/2. Equation (1) becomes

    T^(r)=2i exp(i phi_r) diag(F_00(d_r),-F_11(d_r)).    (11)

Every linear combination (8) is diagonal and hence normal. This is exact for the g^2 coefficient. It is stronger than imposing normality at one numerically tuned parameter point.

[proved, spatial weights in rational form] To evaluate the two diagonal weights without a numerical rank argument, set t=1, a=sqrt(2), c=cosh(a tau), s=sinh(a tau). The scalar convolutions are

    f_0(0,tau)=exp(-d tau)(2 tau+tau c+5s/a)/8,
    f_2(0,tau)=exp(-d tau)(2 tau+tau c-3s/a)/8,
    f_0(1,tau)=f_2(1,tau)=exp(-d tau)(tau c-s/a)/4.

These follow by convolving (1+c)/2, (c-1)/2, and s/sqrt(2), respectively. Put z=2d and define

    I20=2/z^3,
    I2c(b^2)=2z(z^2+3b^2)/(z^2-b^2)^3,
    I2cc=[I20+I2c(4a^2)]/2,
    I1s/a=2z/(z^2-a^2)^2,
    I1cs/a=2z/(z^2-4a^2)^2,
    I0ss/a^2=2/[z(z^2-4a^2)].

For example these are the integrals of tau^2, tau^2 cosh(b tau), tau^2 c^2, tau s, tau c s and s^2, including exp(-z tau). Elementary integration of exponentials proves every formula. Multiplying the convolutions gives

    F00=[4I20+4I2c(a^2)+I2cc+4I1s/a+2I1cs/a-15I0ss/a^2]/64,
    F11=[I2cc-2I1cs/a+I0ss/a^2]/16.                    (12)

Substitution gives

    F00(3)=29543/647101056,    F11(3)=43/296352,
    F00(5)=17656143/22902966128000, F11(5)=123/48668000.

The determinant is strictly positive,

    F00(3) F11(5)-F11(3) F00(5)
       =76553167/21709765967228073000 >0.             (13)

Thus the two spatial readout vectors are not proportional.

[proved] At the stated phases and weights, the hopping eigenvalues are

    t_0=-2 F00(3)+128 i F00(5),
    t_1= 2 F11(3)-128 i F11(5).                        (14)

Equation (13) implies Im(t_0 conjugate(t_1)) !=0, so their phases differ modulo pi. This proves displaced relative centers without rounding a small phase difference.

[proved, local independent tunability] Write the diagonal entries, ignoring their fixed signs and common i, as a_i exp(i phi_A)+b_i exp(i phi_B), with a_i=2 F_ii(3), b_i=128 F_ii(5). At relative phase pi/2 the Jacobian of (arg t_0,arg t_1) with respect to (phi_A,phi_B) has rows

    [a_i^2/(a_i^2+b_i^2), b_i^2/(a_i^2+b_i^2)].

Its determinant is nonzero because a_0/b_0 != a_1/b_1 by (13). The inverse function theorem therefore gives independent phase control in a neighborhood. Normality is preserved throughout that neighborhood by (11).

## 6. Local fields and the n=2 chain

[proved] A marker-to-data diagonal bond can count the occupation of the first data site, and a data-to-marker diagonal bond can count the last. Markers include both M_empty and every occupied M_s, so these operators do not depend on the token position. At n=2 these are the two data occupations. Give them physical coefficients lambda^2 alpha_2^4 g^2 mu_0 and lambda^2 alpha_2^4 g^2 mu_1. They are separately adjustable nearest-neighbor terms on the same 26-state alphabet. Since the vacuum walker factor is unchanged, W G=G W_code exactly. They do not split H0 or produce an order-lambda compression.

[proved] In global fermions f_(x,a), the leading logical chain after subtracting scalar energies is

    H_target=J sum_x {sum_(a=0,1)
       [t_a f_(x,a)^dagger f_(x+1,a)+conjugate(t_a) f_(x+1,a)^dagger f_(x,a)
                    +mu_a f_(x,a)^dagger f_(x,a)]},
    J=lambda^2 alpha_2^4 g^2,
    alpha_2=1/2,
    mu_0=|t_0|/2,   mu_1=-7|t_1|/10.                  (15)

A different positive energy normalization may absorb a common scale in t and mu. For periodic chains the usual fermion-parity boundary twist must be retained when converting from local tensor cells. Open-chain numerics need no such twist.

[proved] The two Bloch energies in units of J are

    epsilon_a(k)=mu_a+2|t_a| cos(k+theta_a),
    theta_a=arg t_a.

Because |mu_a|<2|t_a|, each has two simple zeros

    k_(a,+/-)=-theta_a +/- arccos[-mu_a/(2|t_a|)] mod 2pi.

Their velocities are nonzero, 2|t_a| sin(kF_a). The exact quadratic reference therefore has two gapless bands. Equation (13) proves distinct centers modulo pi. The script gives numerical Fermi-point locations and checks the remaining pair-sum mismatches.

[proved, explicit finite margins] interval_certificate.py encloses arctan(64 F00(5)/F00(3)) in (0.8241,0.8242), arctan(64 F11(5)/F11(3)) in (0.8396,0.8397), and the two kF values in (1.8234,1.8235) and (1.2132,1.2133). Every comparison is rational: sine, cosine and arctangent are bounded by consecutive alternating-series partial sums, and pi is enclosed using 16 arctan(1/5)-4 arctan(1/239). Interval addition then gives disjoint intervals for all four Fermi points. Checking all six bilinear differences gives a circular separation greater than 0.088. Checking all fifteen differences between the six unordered distinct-branch pair sums gives a circular separation greater than 0.030. These are certified lower bounds, not fits to a momentum grid. Pair sums that coincide as the same pair are the allowed forward interactions and are excluded from that comparison.

[heuristic, scope] At fixed n=2, sufficiently small g and then sufficiently small lambda make the analytic higher-order terms small relative to any fixed positive quadratic and quartic mismatch margins of this reference. Exact total U(1) forbids pairing. Standard weak-coupling bosonization would then give a two-boson phase if the Gaussian stiffness stays positive. No interacting thermodynamic calculation, central-charge estimator, or proof of the microscopic CFT limit is supplied here.

## 7. What this does not repair

[numerical, n=2] At finite g the full bond (7) is not exactly quadratic, and its fitted hopping is not exactly normal. checks.py measures the nonquadratic remainder and normality defect at g=.2,.1,.05. The nonquadratic remainder decreases approximately as g^4. The proved bound (9), not an assumption that these terms vanish, is the available control.

[proved, n=3 failure of these rules] At n=3, the same flux pattern gives T_02=T_20=0. With the two inter-cell phases in quadrature,

    |T_10|^2-|T_01|^2
      =4[F_10(3)^2-F_01(3)^2]+128^2[F_10(5)^2-F_01(5)^2]>0.

The strict inequalities have an exact rational verification,

    F_10(3)-F_01(3)=429761/8546653500,
    F_10(5)-F_01(5)=1180141/1487796700500.

Here is a finite algebraic certificate for the calculation. On the scalar two-path space set R=(A tensor I+I tensor A)^-1, v=R|n,0>, P_i=|i><i| tensor I and Q_j=I tensor |j><j|. Differentiating this inverse with respect to independent potentials at i and j gives

    F_ij=2 v^T P_i R Q_j v.

Equivalently this follows by integrating the two Duhamel time insertions defining F. At integer d the 16-by-16 matrix is integer and positive definite. exact_weights.py computes its inverse by rational Gaussian elimination, and reports every F entry. Thus [T,T^dagger]_00<0, proving nonnormality. The relative commutator norm 1.21 is numerical. No optimization over all two-bus geometries has been performed.

[proved, scope] First/last density fields give two independent controls, not a spatially addressed field for every site of an arbitrarily long cell. The construction has supplied neither a general-n normal exchange nor a fixed-alphabet mechanism prescribing all the necessary general-n Fermi points. The one-bus theorem in section 2 cannot exclude these more general mechanisms, since section 5 already escapes it at two flavors.

[assumed, mobile accounting only] If a 26-letter core satisfying all mobile-wrapper hypotheses were available, the cited wrapper formula would give 96*26+48=2544 letters. No such composition is proved here. The present vacuum adjacency walker is not supplied as a positive frustration-free local core for that wrapper, adjacent lengths change the number of stored data bits, and a length-independent second-order readout has not been constructed. Therefore 2544 is conditional alphabet arithmetic, not a realized mobile d0. The frozen demonstration has actual d0=26; the successful arbitrary-n frozen and mobile d0 remain unknown.
