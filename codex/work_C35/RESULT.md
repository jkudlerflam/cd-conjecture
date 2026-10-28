# C35. The logarithmic Casimir bound

2026-10-28. Closed-book work. No network, subagents, or writes outside work_C35.

## 1. Answer and scope

[Proved, conditional on the stated Bethe/CFT inputs] The proposed inequality

    A <= pi M log_2 d,        M = gauge-minimized bond half-spread,

holds for every symmetric SU(N) rank-k rational antiferromagnet, N>=2, k>=1. Within that family equality occurs only for the spin-1/2 Heisenberg chain. The next largest normalized ratio is the fundamental SU(3) chain, with 0.8412396714. Every k>=2 member has ratio strictly below 3/4. Thus the family which killed the dimension-independent bound leaves substantial slack in the logarithmic one.

[Proved, conditional on additional fusion inputs stated in Section 4] The same argument extends to rectangular SU(N) representations k omega_r, including all antisymmetric fundamental representations. No complete classification of arbitrary higher representations is claimed.

[Proved, conditional on supplied or explicitly assumed continuum data] The bound holds for the SO(N) vector rational chain, the full critical zero-field spin-1/2 XXZ line, finite-density gl(1|1), and the quadratic Majorana class established in C14. XXZ approaches equality continuously at its Heisenberg endpoint. Uniform spin twists give equivalent equality examples, so “only Heisenberg” should be understood modulo such transformations and rescaling, and only within the checked classes.

[Assumed/open] There is no universal proof here. The general Sp(2n) rational branch and finite-density higher-rank graded chains remain unresolved. Their missing inputs are identified below. Small-ring Sp computations expose a size-sector obstruction rather than resolve the continuum limit. They do not furnish a counterexample.

[Heuristic] The logarithmic bound is a viable replacement with an optimal coefficient if it is universal. Thermal entropy and elementary energy-current sum rules do not prove it. They leave an explicit low-temperature susceptibility estimate missing. Route A respects it in the recorded weak-coupling regimes, conditional on the existence and infrared identification of the proposed microscopic chains. The new bound would require an unbounded-c construction at fixed d to have a conformal window closing at least as c^{-2}.

Labels apply to each paragraph or displayed derivation under them. “Proved, conditional” means that the algebraic implication is proved, while the separately labeled continuum or fusion input is assumed. Numerical results always include their sizes. An obstruction to an argument is not a counterexample to the conjecture.

## 2. Conventions and finite-dimensional facts

[Definition] Log denotes the natural logarithm. Define

    M(X) = (lambda_max(X)-lambda_min(X))/2
         = inf_t ||X-t I||_op,
    M = inf_o M(h + o tensor I - I tensor o),
    Q = A/(pi M log_2 d).

The infimum is over Hermitian o. All velocities use microscopic site spacing one, unless labeled logical. An actual vacuum sector, or subtraction of the CFT zero-mode energy, is needed when extracting A from a ring. At incommensurate filling the raw 1/L ground-energy coefficient can include a charge-rounding term.

[Proved, finite-dimensional] M(X) is convex and unitarily invariant. If the ordinary site swap S obeys ShS=h, then

    (h_o + S h_o S)/2 = h,
    M(h) <= M(h_o).

Thus reflection-symmetric bonds already minimize the half-spread. Alternatively, simultaneous group twirling removes a coboundary for an invariant bond in an irreducible on-site representation. These statements do not require an assumption about the ground state.

[Proved, finite-dimensional] A useful lower bound which also applies to oriented cell bonds is the following. For orthogonal local basis labels a!=c and b!=e, a coboundary has zero matrix element between |a,b> and |c,e>. Centering by a scalar does too. Therefore

    M >= |<a,b|h|c,e>|.

In particular, an isolated selector move changing both letters, of amplitude t, gives M>=|t| even when the bond is not reflection symmetric. This avoids replacing the optimized M by an unproved estimate from the largest coefficient.

## 3. Symmetric SU(N)_k, all N and k

[Assumed, inherited from C14 Sections 2 to 3] In C14's normalization the antiferromagnetic symmetric-rank-k rational chain has

    d = binom(N+k-1,k),    M=H_k=sum_(j=1)^k 1/j,
    c = k(N^2-1)/(N+k),    v=2pi/N.

The local spectrum is h=sum_(j=0)^k 2(H_(k-j)-H_k) Pi_j, with both endpoints present, so its half-spread is H_k. Reflection gives gauge optimality as in Section 2. The channel decomposition and rational fusion identification are the supplied C14 construction. The Bethe/string thermodynamics and unitary SU(N)_k scaling limit are continuum inputs. C14 and the supplied survey provide separate numerical checks; this work does not rerun their TBA or finite-size Bethe solvers.

[Proved from these inputs] Set

    F(N,k) = A/(pi M) = 2k(N^2-1)/[N(N+k)H_k].              (1)

It remains to prove F<=log_2 binom(N+k-1,k).

### 3.1 The binomial inequality

[Proved, elementary] For positive integers a,b,

    log_2 binom(a+b,a) >= 2ab/(a+b).                       (2)

By symmetry suppose a<=b. The product formula gives

    log binom(a+b,a) = sum_(j=1)^a log(1+b/j)
                    >= a log(1+b/a).

For x>=1, let g(x)=log(1+x)-2 log(2)x/(1+x). Then g(1)=0 and

    g'(x) = [1+x-2 log(2)]/(1+x)^2 > 0.

Substitute x=b/a to prove (2). Equality requires a=b=1. This proves the finite-dimensional counting inequality without asymptotics.

### 3.2 k>=2

[Proved] Apply (2) with a=k,b=N-1. Since H_k>=3/2 and (N+1)/N<=3/2,

    Q <= [(N+1)/(N H_k)] [(N+k-1)/(N+k)] < 1.             (3)

Both factors are displayed to make the strictness explicit. In fact Q<3/4 for every k>=2. To see the uniform margin, write the right side as

    [(N+1)/N] f_k,    f_k=(1-1/(N+k))/H_k.

The inequality f_(k+1)<=f_k follows from

    H_k (k+1) <= (N+k-1)(N+k+1).

Its left side is at most k(k+1), which is smaller than the right side for N>=2. Hence (3) is largest at k=2, where it is

    (2/3) (N+1)^2/[N(N+2)]
      = (2/3)[1+1/(N(N+2))] <= 3/4.

The strictness of (2) for a=k>=2,b=N-1 rules out equality even in the last case.

### 3.3 k=1 and the closest member

[Proved] At k=1, (1) is 2(N-1)/N. Let

    f(x)=log x-2 log(2)(1-1/x).

Then f(2)=0 and f'(x)=[x-2 log(2)]/x^2>0 for x>=2. Thus equality holds only at N=2. Moreover (1-1/x)/log x decreases for x>1, since its derivative has the sign of log x-x+1. The k=1 runner-up is N=3, with

    Q(3,1)=4/[3 log_2 3]=0.8412396714.

Together with Q<3/4 for k>=2, this proves the runner-up statement throughout the symmetric family, not just on a scan.

[Numerical, formula check only] `checks.py` scans 3,998,000 parameter pairs, 2<=N<=2000 and 1<=k<=2000. It finds maximum Q=1 at (2,1), and the next largest value above at (3,1). It checks formulas, not the Bethe ansatz. Selected values are

| N,k | log d | A/(2pi M) | Q |
|:--|--:|--:|--:|
| 2,1 | 0.693147 | 0.5 | 1 |
| 3,1 | 1.098612 | 0.666667 | 0.841240 |
| 2,2 | 1.098612 | 0.5 | 0.630930 |
| 7,2 | 3.332205 | 1.015873 | 0.422633 |
| 8,2 | 3.583519 | 1.05 | 0.406195 |
| 10,10 | 11.433644 | 1.690015 | 0.204909 |
| 100,100 | 135.060089 | 9.637818 | 0.098925 |
| 1000,1000 | 1381.574846 | 66.795998 | 0.067024 |

### 3.4 What the large-d envelope does and does not require

[Proved asymptotics from (1)] Put k=alpha N+O(1) with fixed alpha>0 and let N grow. Stirling and H_k=log k+gamma+o(1) give

    log d = N f(alpha)+O(log N),
    f(alpha)=(1+alpha)log(1+alpha)-alpha log alpha,
    A/(2pi M) ~ alpha N/[(1+alpha)log N].

Consequently

    [A/(2pi M)] [log log d/log d]
        -> alpha/[(1+alpha) f(alpha)].                   (4)

[Proved] The coefficient in (4) is maximized at alpha=1, with value 1/(4 log 2). Put p=alpha/(1+alpha). The coefficient is p(1-p)/h(p), where h is binary entropy. The inequality h(p)>=4 log(2)p(1-p) follows, for x=2p-1, from

    h(p)=log 2-sum_(m>=1) x^(2m)/[2m(2m-1)],
    h(p)-log(2)(1-x^2)
      =sum_(m>=2) [x^2-x^(2m)]/[2m(2m-1)] >=0.

For 0<p<1 equality occurs at p=1/2. Thus the balanced rank/level ray is best among fixed-ratio rays. In particular k=N gives

    A/M ~ [pi/(2 log 2)] log d/log log d,
    Q ~ 1/[2 log log d].                                 (5)

[Proved, order of the full symmetric envelope] The envelope over all d(N,k)<=D is Theta(log D/log log D). For the upper bound let m=min(k,N-1) and L=log D. The product proof above gives log d>=m log 2, while

    A/(2pi M) <= min(k,N)/H_k <= (m+1)/H_k,
    H_k >= log(k+1) >= log(m+1).

If m<=sqrt L the ratio is at most sqrt L+1. Otherwise it is at most 2(L/log 2+1)/log L. The balanced sequence supplies the lower bound using (5). This proves the order uniformly. It does not prove the optimal leading constant of the unrestricted envelope, although (4) determines it on every fixed-ratio ray.

[Conclusion, proved for this family] The symmetric counterexamples force growth at least log d/log log d, but do not force a full log d. The logarithmic conjecture is therefore sharp at small d and not asymptotically saturated by this family.

## 4. Higher SU(N) representations

### 4.1 Rectangular on-site representations

[Assumed fusion inputs, beyond C14's symmetric proof] For the standard rational antiferromagnetic chain in R=k omega_r, 1<=r<=N-1, use the following two inputs.

1. Its complete infrared theory is SU(N)_k with velocity 2pi/N in the same spectral-parameter convention. The source is at nesting node r; the rank kernel has denominator sinh(N|omega|/2), giving this endpoint velocity. The identification of all auxiliary sectors remains the standard fusion/TBA input, as in C4-P2 Section 5.
2. The fused R matrix restricts on the highest-weight alpha_r string to the spin-k/2 rational R matrix. Equivalently, on pair channels with highest weights 2k omega_r-j alpha_r, 0<=j<=k, the braided eigenvalue ratios are products of (ell-u)/(ell+u), ell=k-j+1,...,k.

The second statement is the usual rank-one fusion restriction. It is made explicit because the supplied C14 calculation proves it only for r=1. This report does not supply a construction of the full rectangular Yangian R matrix. It must not be read as an independently proved integrability theorem for arbitrary rectangles.

[Proved conditional local implication] The logarithmic derivative on these channels has eigenvalues differing by 2H_k from j=0 to j=k. Hence the full bond half-spread is at least H_k. Group twirling shows that a coboundary cannot reduce it. The rank-one string is a legitimate finite subspace: the weight space with total lowering only along alpha_r cannot contain a lowering root involving any other simple root. The extra fusion input fixes the action on that subspace, not just its dimension.

[Proved, finite-dimensional representation theory] The Weyl dimension formula is

    d_r(k)=product_(i=1)^r product_(j=r+1)^N (1+k/(j-i)).  (6)

It obeys d_r=d_(N-r), and cancellation of the factors common to adjacent heights gives

    d_(r+1)(k)/d_r(k)
      = product_(t=1)^(N-r-1)(1+k/t)
        /product_(t=1)^r(1+k/t).

This is at least one up to the middle height. Therefore

    d_r(k)>=d_1(k)=binom(N+k-1,k).

The central charge and velocity in the assumed infrared input are independent of r. Combining d_r>=d_1 and M_r>=H_k with Section 3 proves Q_r<=Q_1<=1. Equality requires k=1,N=2. Conjugate rectangles give the same comparisons.

[Numerical, finite formula check] `checks.py` verifies the dimension inequality for 2<=N<=60, 1<=k<=30 and every allowed height. It is not a check of the fusion input or the continuum identification.

### 4.2 Antisymmetric representations and a normalization check

[Assumed, standard antisymmetric fusion spectrum] For R=wedge^r C^N put s=min(r,N-r). The pair channels, after replacing r by s through conjugation, have diagrams [2^(s-j),1^(2j)], j=0,...,s. Their braided factors are product_(ell=1)^j (ell-u)/(ell+u), so the derivative eigenvalues are -2H_j. This gives M=H_s, stronger than the M>=1 needed in Section 4.1. With c=N-1,v=2pi/N,

    A/(pi M)=2(N-1)/(N H_s),
    d=binom(N,r)>=N.

The dimension inequality is proved by the adjacent-rank binomial ratios. The Casimir inequality then follows directly from the fundamental case.

[Proved consistency check conditional on those spectra] SU(4) in its six-dimensional antisymmetric representation gives local energies 0,-2,-3 and M=3/2. SO(6)'s vector bond in Section 5, shifted by -I, has the same energies, with c=3 and v=pi/2. SU(2) at arbitrary k recovers the spin-S=k/2 Takhtajan-Babujian family already included in Section 3. These checks fix the different harmonic orderings for symmetric and antisymmetric fusion.

[Assumed/open] Nonrectangular evaluation representations can have multiplicities and multiple seas. Neither their full local spectra nor a complete velocity/central-charge dictionary is supplied here. The rectangular result is not a theorem covering them.

## 5. Other solvable families

### 5.1 SO(N) vectors

[Proved, finite-dimensional] For N>=3 take

    h=P-2K/(N-2),    K=N|Omega><Omega|,
    |Omega>=N^(-1/2) sum_a |a,a>.

The symmetric-traceless, antisymmetric and singlet energies are respectively

    1, -1, -(N+2)/(N-2).

Their half-spread is M=N/(N-2), and reflection proves gauge optimality.

[Assumed infrared input, as in C14] The antiferromagnetic rational point has c=N/2 and v=2pi/(N-2). Therefore

    A/(pi M)=1,    Q=1/log_2 N<1.                        (7)

[Proved checks from the inputs] At N=3 this is TB with twice C14's rank-two SU(2) normalization, M=3,v=2pi,c=3/2. At N=4 the SO(4) vector is (2,2) under SU(2) times SU(2). Since K=(I-P_1)(I-P_2), P=P_1 P_2, the bond is P_1+P_2-I. Thus M=2 and A=2pi, in agreement with two Heisenberg copies. SO(2) is excluded because the rational formula is singular.

### 5.2 XXZ and twists

[Proved, finite-dimensional] Normalize the spin-1/2 bond as

    h=(sigma_x tensor sigma_x+sigma_y tensor sigma_y
                      +Delta sigma_z tensor sigma_z)/2.

For -1<=Delta<=1 its eigenvalues are Delta/2 twice and -Delta/2 +/-1. Thus M=1 throughout this interval, including after coboundary minimization by reflection.

[Assumed Bethe infrared input] For -1<Delta<=1, put Delta=cos gamma, 0<=gamma<pi. At zero field c=1 and v=pi sin(gamma)/gamma, with the continuous endpoint value v=pi. Therefore

    Q=sin(gamma)/gamma<=1,

with equality only at gamma=0. The Delta=-1 endpoint has quadratic modes and is not included as a CFT.

[Proved transformation statement, conditional on the same thermodynamic limit] A uniform phase on spin-flip hopping is removed by a site-dependent z rotation, leaving a boundary twist on a ring. The two-site spectrum is unchanged. To see gauge optimality even though ordinary reflection is lost, conjugate the bond swap by the opposite z rotations that generate the phase. This twisted swap fixes h and sends an arbitrary coboundary into another coboundary with its two terms reversed and rotated. Alternatively twirl the coboundary over the preserved on-site U(1), leaving only o proportional to sigma_z, and then apply the twisted swap. Convexity removes it. The thermodynamic velocities remain those above, although finite-ring zero modes can be twisted. Thus the equality class includes these twisted Heisenberg representatives when A means the vacuum anomaly coefficient.

### 5.3 Free fermions and a finite-density graded chain

[Proved conditional on the quadratic continuum identification, inherited from C14 Section 5] Every even nearest-neighbor quadratic Majorana/BdG chain satisfies A<=pi M. C14 proves the optimized-norm statement through parity projections, a trace-norm bound on intersite hopping, and number-conserving doubling. No assumption that the optimal gauge is quadratic is needed. Since d>=2 this implies the proposed logarithmic inequality. A grouped microscopic period must use grouped-site velocities and its own bond norm.

[Proved, finite-dimensional] The gl(1|1) graded swap with a chemical potential is, after a scalar and field reparametrization, the hopping bond

    h=c_1^dagger c_2+c_2^dagger c_1
                          -(mu/2)(n_1+n_2-1).

Its two-site eigenvalues are +/-1 and +/-mu/2. The ordinary qubit matrix is reflection symmetric, so M=max(1,|mu|/2). The on-site chemical potential has been included in the norm.

[Proved within free dispersion theory] The symbol is 2 cos q-mu. For |mu|<2 it has two simple Fermi points, c=1 and

    v=sqrt(4-mu^2),    Q=sqrt(4-mu^2)/pi<=2/pi<1.

At |mu|=2 the branch is quadratic. This gives an analytic finite-density graded example rather than treating all graded chains as excluded.

[Proved exclusion, inherited and directly visible] For a mixed graded swap P_g at zero chemical potential, a uniform odd-color product attains the minimum of +P_g. Replacing one color by an even color gives energy 2+2cos q relative to that product. For -P_g the even-color product has the analogous 2-2cos q branch. These quadratic low-energy branches violate the pure z=1 CFT premise. They are not counterexamples to the proposed inequality.

[Assumed/open] Higher-rank finite-density t-J and gl(m|n) chains need a dressed-velocity calculation at the specified chemical potentials, with those potentials included as half an on-site term on each bond. C14 supplies no such analytic solution. The casimir_survey entry at d=3 and density 2/3 uses only L=6,12 and reports A/(2pi M) about 0.42, hence Q about 0.53. This is inherited numerical evidence at those sizes, not an all-density check.

[Proved exclusion of the displayed realization, inherited from C4-D1 Section 6] The conventional mixed orthosymplectic rational density in those notes is not Hermitian in the positive tensor-product inner product. Its anti-Hermitian entries change both site colors, so a coboundary cannot cancel them. This excludes that displayed realization only. It does not exclude another Hermitian realization of a superalgebra model.

### 5.4 Sp(2n), the unresolved analytic branch

[Proved, finite-dimensional] In the fundamental representation let

    h=P+K/(n+1),    K=|J><J|,
    |J>=sum_(a=1)^n (|a,a+n>-|a+n,a>).

Then K has eigenvalue 2n on the symplectic singlet. The symmetric, antisymmetric-traceless and singlet energies are

    1, -1, (n-1)/(n+1).

For n>=2, M=1 by reflection. For n=1 the traceless antisymmetric sector is absent, M=1/2 and the bond is P/2+I/2. It is the rescaled Heisenberg equality case.

[Assumed/open, precise obstruction] For n>=2, the local norm is settled, but the supplied sources do not establish the entire infrared sector content and its velocities. C4-D1 explicitly leaves this open. C4-D3 and the later survey find split-velocity behavior in Sp(4)-invariant chains. Inserting c=n(2n+1)/(n+2) and an assumed common velocity would miss any additional sector. The analytic task for this branch reduces to proving

    sum_a c_a v_a <= pi log_2(2n)

from its actual nested Bethe solution, not from its global symmetry.

[Proved conditional sufficient criterion, not an identification of this chain] If an Sp(2n) fundamental chain were Sugawara complete at level k>=1, C14's current sum rule with C_2(fund)=(2n+1)/4 would give

    A/(pi M)<=2(2n+1)/(k+n+1)
             <=4-6/(n+2).

For n>=4 this is at most log_2(2n). Equality of the last comparison occurs at n=4; its difference increases thereafter since 1/[n log 2]-6/(n+2)^2>0. This conditional shortcut cannot control a neutral sector and is not used to certify the rational Sp branch.

[Numerical, new local checks] `checks.py` diagonalizes SO bonds at d=3,...,8 and Sp bonds at d=2,4,6,8. It verifies the spectra above and verifies that swap averaging removes random unrestricted Hermitian coboundaries. Largest local matrix size is 64. These checks supplement the algebraic proof.

[Numerical, new ring calculation with a failed Casimir estimator] `sp_ed.py` diagonalizes the zero-Cartan sectors for Sp(4) at L=4,6,8,10, Sp(6) at L=4,6,8, and Sp(8) at L=4,6. The largest sector has 63,504 states. The maximum eigenvector residual over all runs is 8.04e-11, and the combined eigensystem runs take about 1.3 seconds on this machine. Full-Hilbert-space checks at (d,L)=(4,4),(4,6),(6,4),(6,6),(8,4) verify that the same ground energies occur there. Larger rings still assume that the zero-Cartan sector contains the ground state. The Hamiltonian is built directly from P and K, independently of the survey implementation.

[Numerical, failed check reported explicitly] Fitting E0/L=e+a/L^2 and setting A=-6a/pi is not stable on these rings. For Sp(4), L=6,8,10 gives A=-3.21762, while L=8,10 gives A=15.92244, an apparent Q=2.53414. Adding a 1/L^4 coefficient produces A=53.49562. Sp(6) at L=6,8 similarly gives an apparent Q=2.09057. These fits are rejected. For Sp(4), restricting to L=4,8 gives Q about 0.780, but two small sizes do not determine an asymptotic coefficient. Strong oscillation between length classes is directly visible in the energies. Its assignment to particular CFT zero-mode sectors remains unproved here.

[Conclusion of this check] No Sp violation or general Sp confirmation follows. A finite-size Bethe treatment identifying the vacuum length classes, or sufficiently large commensurate sizes, is the concrete missing calculation. Repeating unconstrained fits to these small rings would not resolve it.

## 6. Stacking

[Proved, finite-dimensional] Suppose the site is a tensor product of factors of dimensions d_i, with a decoupled bond sum h=sum_i h_i. Let E_i be normalized partial trace over every other factor at both sites, followed by adjoining their identities. It is a unital positive map. For every Hermitian X, its spectral interval is contained in that of X, hence M(E_i X)<=M(X). Applied to a general, possibly entangling on-site coboundary O,

    E_i(h+O tensor I-I tensor O)
       = h_i+o_i tensor I-I tensor o_i+scalar.

Therefore M_stack>=max_i M_i. If each factor satisfies the logarithmic bound, then

    A_stack=sum_i A_i
      <=pi sum_i M_i log_2 d_i
      <=pi M_stack sum_i log_2 d_i
       =pi M_stack log_2(product_i d_i).

The continuum additivity of A is the decoupled-product input. The norm inequality itself is finite-dimensional and covers arbitrary coboundary gauges. For reflection-symmetric factors the stronger equality M_stack=sum_i M_i follows from spectral addition and reflection. For example, m identical Heisenberg copies have Q=1/m, not equality at d=2^m.

## 7. Thermal and energy-current routes

### 7.1 What entropy proves

[Proved, finite-dimensional] For a chain of L d-state sites, positivity of relative entropy against I/d^L gives S<=L log d. Thus the thermodynamic entropy density satisfies s(T)<=log d. Center a representative bond at its spectral midpoint. Its periodic sum has spectrum in [-LM,LM], since coboundaries telescope, and its total thermal energy variation per site is at most 2M. These arguments also hold by taking a sequence of representatives approaching the optimized infimum.

[Assumed conformal window] Let B=sum_a c_a/v_a and suppose

    s(T)>=(1-epsilon) pi B T/3,    0<T<=T_*.

The finite-dimensional entropy bound then gives

    B <= 3 log d/[pi(1-epsilon)T_*].                     (8)

It controls B, not A. A velocity bound v_a<=v_LR gives A<=v_LR^2 B, which still depends on T_*. Taking T to zero removes rather than strengthens this estimate.

[Proved obstruction to inference from these thermodynamic inequalities alone] Fix d and a nonzero velocity v_0 below any proposed velocity cap. For arbitrarily large B define the thermodynamic test function

    T_0=3 log d/(pi B),
    s(T)=log d tanh(T/T_0).

It has the required low-T slope pi B/3, positive heat capacity, and s<=log d. Its total thermal energy rise is

    integral_0^infty T s'(T)dT
       =T_0 log d integral_0^infty x sech^2 x dx
       =T_0 log d log 2,

which becomes arbitrarily small as B increases. Assigning c=Bv_0 would give A=Bv_0^2 arbitrarily large while respecting these inequalities. This is a test of the inference, not a proposed lattice Hamiltonian or a claim that its entire equation of state is realizable in finite dimension. It shows why entropy capacity, energy capacity and a velocity cap alone cannot prove the desired result. Even the chiral current pi c T^2/12 becomes small at T of order T_0.

### 7.2 The quantity a thermal-current proof would need

[Proved CFT calculation, conditional on the conformal description] In one sector, J_E=v(H_R-H_L). Each chirality has energy density pi c T^2/(12v), so thermal differentiation gives

    Var(H_R)/L = pi c T^3/(6v),
    Var(J_E)/L = pi c v T^3/3.

Independent sectors add. Thus

    Var(J_E)/L = pi A T^3/3.                             (9)

A sufficient new estimate, if proved for the appropriate low-energy lattice current weight, would be

    limsup_(T->0) Var(J_E)_IR/(L T^3)
       <= [pi^2/(3 log 2)] M log d.                     (10)

Then (9) would yield the proposed coefficient. Equation (10) is a missing susceptibility theorem. Entropy bounded by log d is not a proof of it. On the lattice the total current can have high-energy matrix elements or nonzero ground-state fluctuations, so replacing the IR weight in (10) by the unprojected equal-time variance is an additional and generally invalid step.

[Proved local estimate; assumed use of a lattice thermal sum rule] The centered bond current j=i[h_12,h_23] obeys ||j||<=2M^2 by the commutator norm inequality. A thermal contact sum rule built from nested energy commutators has a local coefficient of order M^3. In the convention D_th=beta^2 Var(J_E)/(2L), the CFT ballistic contribution is pi A T/6. Bounding the entire contact term in norm gives a weight of order M^3/T, hence only

    A <= constant M^3/T^2.

The detailed thermal sum-rule form is a framework input, not derived here. Its dimensional mismatch already locates the failure: one needs a positive low-temperature estimate of order M T^2 log d for the relevant contact term, after handling regular UV weight. Subtracting its zero-temperature value does not preserve spectral positivity. At a nonintegrable lattice fixed T, scattering can broaden the current peak, so a literal Drude delta function is not automatically the CFT weight either.

[Proved CFT diagnostic, inherited normalization from C14] A pure stress-tensor lattice density has inverse-frequency moment

    m_(-1)^CFT(q)=A q^2/(24pi).

A bound on its total positive inverse moment would be sufficient, but inverse excitation energies are not bounded by a local operator norm as the gap closes. Density improvements h_j -> h_j+o_j-o_(j+1) can change this q^2 coefficient through dimension-one currents without changing A. An argument using this moment must first choose and justify a stress-tensor representative. Minimizing M does not perform that projection.

### 7.3 A possible mechanism, and the missing assumption

[Heuristic] A dimension-logarithmic coefficient could arise from a bound on how much low-energy stress response a d-state bond can transmit. Entropy production bounds across a cut have the right M log d scale. Relating that rate to A, however, requires a state whose entanglement-growth rate samples all conformal modes at a fixed microscopic smoothing length. A continuum quench gives a rate proportional to A divided by that length. The CFT is only valid once the length exceeds its onset length. Using it at one lattice spacing would add the very uniform-window assumption the conjecture is meant to avoid.

[Conclusion] Thermal/current and entanglement-rate arguments suggest a useful target estimate, (10), but currently provide no physics-level derivation of a universal C independent of the conformal window. The log form remains a conjecture rather than a consequence of s<=log d.

## 8. Consequences if A<=C M log d

[Proved conditional implication] Assume C is a dimension-independent constant in natural-log units and (8) holds. Cauchy-Schwarz gives

    c^2=(sum_a sqrt(c_a v_a) sqrt(c_a/v_a))^2 <= A B.

With R=M/T_*,

    c^2 <= (3C/pi) R (log d)^2/(1-epsilon).              (11)

For the proposed sharp value C=pi/log 2 this becomes

    c^2 <= 3 R (log d)^2/[(1-epsilon)log 2],
    R >= (1-epsilon)log(2)c^2/[3(log d)^2],
    T_*/M <= 3(log d)^2/[(1-epsilon)log(2)c^2].           (12)

The last line constrains an existing conformal chain of charge c. Conversely, a sufficient condition for c<=d-1 is

    T_*/M >= 3(log d)^2/
                    [(1-epsilon)log(2)(d-1)^2].         (13)

At epsilon=0.1 substitute 1-epsilon=0.9. The earlier failed 2pi conjecture gave one power of log d in (11); the replacement gives two. Keeping the old coefficient-6 conclusion would be an algebraic error.

[Proved further consequences] The central-charge-weighted mean velocity satisfies A/c<=C M log d/c. Thus at fixed d, unbounded c requires this average to vanish at least as 1/c. If all velocities are bounded below by v_min>0, then c<=C M log d/v_min. Without a lower velocity bound or a lower bound on T_*/M, the conjecture supplies no universal upper bound on c at fixed d. It is compatible with the intended fixed-d Route A program.

## 9. Route-A cells

### 9.1 What is established without continuum assumptions

[Proved, finite-dimensional second order] Let P project onto two decoupled cell codes, Q=I-P, and Q(H_0-E_0)Q>=Delta Q. Write the interface as lambda V with ||V||<=g and PVP=0. Then

    K_2=-lambda^2 P V Q[Q(H_0-E_0)Q]^(-1)Q V P,
    -lambda^2 g^2/Delta P <= K_2 <=0,
    M(K_2)<=lambda^2 g^2/(2Delta).                       (14)

This follows from the inverse bound 0<=R<=Q/Delta and (QVP)^dagger(QVP)<=g^2 P. It is a statement on a finite pair, not a global isolated low-energy-band theorem. The recorded cell constructions make cross-interface three-cell terms vanish at this order. If a construction does not, (14) alone does not certify a nearest-neighbor logical target. Retained on-cell fields must be included separately; they often are not negative semidefinite.

[Proved rescaling conditional on a logical CFT] A frozen cell of length b changes a logical ring length L_c into physical length L=b L_c. Its oscillator Casimir term obeys

    -pi A_log/(6L_c)=-pi (b A_log)/(6L),
    A_phys=b A_log.                                     (15)

Thus the required test is b A_log<=pi M_phys log_2 d_0. The microscopic dimension is d_0, not the code dimension 2^n or N. Applying a logical bound with log dim(code) and stopping there would not establish the microscopic conjecture.

### 9.2 Explicit small-coupling conditions

[Proved for the specified projected free target] For n diagonal hopping bands t_a e^{i theta_a}, with fields placing them inside their bands,

    v_a=2|t_a| sin k_F,a <=2|t_a|,
    A_payload,phys <=2b sum_a |t_a|.                    (16)

For the pure cross-hopping bond the many-body half-spread is sum_a |t_a|, since a two-site singular-value decomposition gives independent modes with energies +/-|t_a| and zero. The projection bound of C14 also gives this as a lower bound on the full quadratic bond norm including fields and optimized coboundaries. For more general quadratic targets one can use its universal A_log<=pi M_log instead.

[Proved sufficient condition for a projected model; heuristic extension to weak interactions] If |t_a|<=W, a sufficient payload condition is

    2bnW <= pi M_phys log_2 d_0.                        (17)

Taking lambda sufficiently small at each fixed n,b enforces this whenever W is of order lambda^2/Delta times fixed port factors and M_phys has a nonzero lower bound independent of lambda. The selector matrix element in Section 2 supplies such a lower bound. An interacting infrared theory close to this free target has nearby velocities only under the perturbative RG assumptions of the construction. This last step is not finite-dimensional algebra.

[Proved check for the old swap target, conditional on its CFT] For effective SU(N) fundamental exchange J, (15) gives

    A_phys=2pi bJ(1-1/N).

The C7 selector condition J much smaller than t/b^4, together with M_phys>=t from its two-letter move, makes A_phys/(pi M_phys) much smaller than 2/b^3. This is the recorded old construction's scale hierarchy, not a claim that every possible cell construction must obey it. The SU(N) target's stability problem remains as recorded in STATUS.

### 9.3 The mobile length mode must be counted

[Assumed exact projected-fluid input, from C8 length_fluid Sections 2 to 4] If cell lengths are b or b+1, with fraction nu of long cells and length-exchange hopping J_len, the physical mean spacing is b_bar=b+nu. Its c=1 velocity is

    v_len=2J_len b_bar sin(pi nu).

The equal-payload construction has all payload velocities rescaled by b_bar too. Consequently a sufficient leading free condition is

    2 b_bar [J_len+sum_a |t_a|]
        <=pi M_phys log_2 d_0.                         (18)

The gluing of ring momenta changes allowed zero modes, not this sum of oscillator anomaly coefficients. An apparent constant from raw energies on incompatible rings need not equal it.

[Heuristic, using the recorded hierarchy] In routeA_assembly Section 3, J_len is chosen comparable to W, and all payload and marker scales lie below U=O(t/b^4). In addition, storing n independent logical complex fermions in a b-site cell requires 2^n<=d_0^b, so n<=b log_2 d_0. Even the loose estimate W<=U then makes (18)'s ratio O(b^{-2}) at fixed d_0. The actual second-order and visibility factors reduce it much further. This is how the length mode can be added without threatening the bound. It does not follow merely from the fact that the cells are mobile.

[Proved distinction] Small lambda/Delta alone, with b allowed to grow independently, does not imply (17) or (18). One must compare the physical rescaled sum of velocities with M_phys. The recorded constructions do this through their selector hierarchy and n-dependent choices of small couplings. A hypothetical future construction with a different hierarchy must be checked anew.

### 9.4 Check against the actual recorded versions

[Heuristic scale arithmetic on the assembly assumptions] routeA_assembly gives W=7.2e-22 at n=2,b=10 and W=2.0e-27 at n=3,b=12 in units t=1 for its less restrictive worked hierarchy. Treating W as the largest hopping and taking J_len<=W, (18) gives the conservative bounds

    A_phys <=2(b+1)(n+1)W
        =4.752e-20       for n=2,
        =2.08e-25        for n=3.

These are not parameters of a completed spatial compiler. They show that the arithmetic example has ample room even with M_phys>=1 and the weakest nontrivial log_2 d_0>=1. Its more restrictive hierarchy only decreases these numbers.

[Numerical input from the spatial-cell README, n=1,2,3; heuristic infrared use] The 23-letter folded-bus construction has b=n+1 and a reported hopping norm at lambda=1 of 3.6e-7,1.1e-8,1.3e-10. For an n by n hopping matrix T, ||T||_*<=n||T||_op. C14's hopping trace-norm estimate gives A_log<=pi||T||_* for its quadratic projection. Using M_phys>=t=1 therefore bounds the quadratic payload ratio by

    Q_quad <= b n ||T|| lambda^2/log_2 23,

which gives 1.59e-7 lambda^2, 1.46e-8 lambda^2 and 3.45e-10 lambda^2, respectively. If the table's norm is a Frobenius norm, it still upper-bounds the operator norm used here. The table also records non-Gaussian terms, so these are estimates for the quadratic projection, not proved bounds on the full interacting chain's A. The chosen small couplings and the claimed perturbative infrared regime are consistent with the conjecture. No new thermodynamic simulation was performed here.

[Assumed/open status audit] The 48-letter central-switch compile variant shows n conducting quadratic bands through n=4 and displaced centers at n=2, but its equal-arm selector, general phases and infrared control remain open. C28 now gives a positive frustration-free 17-letter initialized core and a two-mode doped candidate, while proving a restricted field-rank obstruction for n>=5. It expressly does not complete the all-n fixed-alphabet construction. These later facts supersede a literal reading of older “assembled end to end” language. None establishes a finite-coupling all-n microscopic CFT for which one could claim a universal, rigorous A measurement.

[Conclusion, conditional physics statement] All recorded controlled Route-A regimes are compatible with the logarithmic bound. The mechanism is the suppression of the physical sum of velocities, including the length fluid, by the second-order scale and the selector hierarchy. Frozen-cell degeneracy still puts those versions outside a strict unique-vacuum interpretation of the original chain class. On compatible rings the payload comparison (15) remains meaningful.

## 10. Reproduction, limits, and next steps

[Numerical] Run from this folder with numpy and scipy only.

    OPENBLAS_NUM_THREADS=1 python3 checks.py
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 sp_ed.py

`checks.json` contains the formula scan and local-algebra checks. `sp_ed.json` contains every ring energy, eigenvector residual, full-space control, and the rejected fit diagnostics. No result depends on a large local computation, cluster job, downloaded source, or imported project code. The scripts write only in this folder.

[Assumed/open] The remaining substantive tasks are (i) derive the complete Sp(2n) rational thermodynamics and its vacuum sectors, (ii) solve the finite-density higher-rank graded comparison with chemical potentials included, (iii) independently establish the rectangular rank-one fusion input in the normalization used here, and (iv) prove a low-energy stress/current susceptibility estimate such as (10), or find a chain that violates it in the required gauge-invariant form. Arbitrary nonrectangular representations are also outside the completed survey.

[Summary of status] The all-N,k symmetric-family inequality and its closest member are settled analytically given C14's continuum input. The thermal argument and the remaining integrable families have located gaps. The cells do not currently pressure the conjecture, because their large central charge is carried by very small velocities. Project STATE.md was not edited, as the user prohibited writes outside this folder. PROGRESS.md is the local resumable record.
