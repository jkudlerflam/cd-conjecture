# C110. Vacuum prescriptions, gauge covariance, and a gapped cutoff obstruction

2026-10-29. Claims are tagged [Proved], [Numerical], [Assumed], [Heuristic], or [Open]. Tags govern the argument that follows them. This note has not had an independent referee review.

[Open] This investigation does not prove a finite universal F(d), construct a fixed-d family with divergent physical A, or define a finite upper semicontinuous extension agreeing with physical A on every conformal bond. It establishes obstructions to specific proposed prescriptions and a finite, gauge-invariant family of upper semicontinuous cutoff observables. Identifying an infimum of those observables with physical A remains open.

## 1. Inherited conventions and the logical target

[Proved, inherited from C107] Write G for the scalar shifts and on-site coboundaries,

    G = {alpha I + o tensor I - I tensor o}.

The bond scale is the quotient norm M([h]) = inf_(g in G) ||h+g||. The quotient unit ball and unit sphere are compact. For periodic ground energies E_L and the thermodynamic energy density e,

    Q_L(h) = (6/pi) L [L e(h)-E_L(h)]

is gauge invariant and obeys

    |Q_L(h)-Q_L(k)| <= (12/pi) L^2 M([h]-[k]).

These are C107 PROOF Sections 1 and 2. They are not re-derived here. The physical target on a conformal bond is A=sum_a c_a v_a, with vacuum energy E_L=Le-pi A/(6L)+o(1/L). Gapped bonds must be assigned zero.

[Proved] A finite real-valued upper semicontinuous A on the quotient unit sphere would be bounded there. Positive homogeneity would then give A(h)<=F(d)M(h) for M(h)>0. At M=0 the uniform Hamiltonian is a scalar and the desired value is zero. If the physical A is unbounded on the sphere, compactness supplies a convergent sequence with A diverging, as established in C107 Section 7.

[Proved scope] Failure of upper semicontinuity alone is not the divergent-A alternative. A bounded function can fail upper semicontinuity. Likewise, failure of one proposed vacuum prescription says nothing by itself about boundedness of physical A. The two requested outcomes are sufficient resolutions of the bound, not an exhaustive list of possible outcomes of a test of prescriptions.

## 2. The direction of twist optimization

[Proved] Fix L and a collection of sectors or boundary conditions indexed by theta. If their energies are E_L(theta) and the common bulk density is e, then

    inf_theta Q_L(theta) = (6/pi)L[Le-sup_theta E_L(theta)],
    sup_theta Q_L(theta) = (6/pi)L[Le-inf_theta E_L(theta)].

Thus an infimum of Q_L selects the highest ground energy among the boundary conditions. Minimizing energy requires a supremum of Q_L. This distinction applies to the sign convention in the task and in C107.

[Proved] The error is physical already for the ferromagnetic XX chain at d=2,

    H = -t sum_j (sigma_j^+ sigma_(j+1)^- + h.c.),  t>0,

with a U(1) flux phi through the seam. Take -pi<=phi<=pi and even L. In a sector of N spins down, the Jordan-Wigner boundary condition at zero flux is periodic for odd N and antiperiodic for even N. The N filled momenta can therefore be centered at phi/L, modulo a displacement by 2pi/L. Choosing the displacement with smallest absolute center gives phi/L in this flux interval. Summing the occupied cosines gives

    E_L(N,phi) = -2t cos(phi/L) sin(pi N/L)/sin(pi/L).

At phi=+-pi the two centers are degenerate. Since cos(phi/L)>0 and L is even, minimizing over N selects N=L/2. Hence

    E_L(phi) = -2t cos(phi/L)csc(pi/L),
    e = -2t/pi.

The expansion of the product is

    cos(phi/L)csc(pi/L)
       = L/pi + [pi/6-phi^2/(2pi)]/L + O(L^-3),

uniformly for phi in the compact interval [-pi,pi]. Therefore

    Q_L(phi) = 2t - 6t phi^2/pi^2 + O(L^-2).

The physical value is A=2t, inherited from C107's two Fermi velocities and c=1. Nevertheless,

    lim_(L even) inf_phi Q_L(phi) = -4t,
    lim_(L even) sup_phi Q_L(phi) = 2t.

Replacing the infimum by its positive part gives zero, not 2t. The twist is an actual conserved U(1) twist, so this example does not rely on an arbitrary impurity at the seam. It already uses the vacuum-compatible even residue class.

[Numerical] Nine spin matrices at L=4,6,8 and phi=0,pi/3,pi agree with the formula. These checks appear in CHECKS.md.

## 3. Infima over lengths do not repair finite-size errors

[Proved, using C107 Section 5] For the critical XY family in C107, at fixed gamma>0 on even lengths,

    Q_L = gamma + a_gamma/L^2 + O(L^-4),
    a_gamma = 7pi^2(4gamma^2-3)/(240gamma).

At gamma=1/2, a_gamma<0 and physical A=gamma. For every finite cutoff N, some even L>=N satisfies Q_L<A. Consequently,

    inf_(L>=N, L even) Q_L < A

for every N. This remains true on any infinite residue class contained in the even lengths. The vacuum is already correctly selected in this example. The problem is the finite-size correction.

[Proved] Taking the supremum over N recovers the asymptotic lower limit, but loses the desired elementary upper semicontinuity argument:

    liminf_(L even) Q_L = sup_N inf_(L>=N, L even) Q_L.

Each inner infimum is upper semicontinuous. An increasing supremum of upper semicontinuous functions need not be upper semicontinuous. This is not remedied by removing a finite number of small rings. C107 also rules out a universal C M/L^2 correction restoring decreasing monotonicity for these estimators.

[Proved scope] An infimum over all finite scales could yield the desired topology if every finite-scale value were an upper bound on A and the bound were asymptotically sharp. The raw vacuum Casimir quantities fail that necessary upper-bound property. No conclusion about a different, rigorously corrected observable follows from this example.

## 4. A twisted seam is not automatically a map on the quotient

[Proved] Orient the seam as (L,1) and define a literal on-site seam twist by

    H_L^U(h) = sum_(j=1)^(L-1) h_(j,j+1)
                  + (I tensor U) h_(L,1) (I tensor U^dagger).

For h'=h+o tensor I-I tensor o, the interior coboundaries sum to o_1-o_L. The twisted seam coboundary is o_L-Uo_1U^dagger. Thus

    H_L^U(h')-H_L^U(h) = o_1-Uo_1U^dagger.                 (1)

The scalar shift is harmless, but this boundary operator is generally nonzero. A finite O(1) seam energy becomes O(L) in Q_L.

[Proved] An explicit gauge orbit consists of

    h_r = -(Z tensor I+I tensor Z)/2
              + r(Z tensor I-I tensor Z),   r real.

Every member gives the same periodic Hamiltonian H=-sum_j Z_j, with unique product ground state, gap 2, e=-1, and physical A=0. Its quotient norm is M=1. Indeed, the r=0 representative has norm 1, while the expectations in |00> and |11> differ by 2 for every coboundary and scalar shift, proving the matching lower bound.

Put a=-1/2+r and b=-1/2-r. Twisting the seam leaves sites 2 through L with field -Z and gives site 1 the local operator aZ+bUZU^dagger. Thus

    E_L^U(h_r) = -(L-1)-||aZ+bUZU^dagger||.

Optimizing over all U in U(2),

    max_U ||aZ+bUZU^dagger|| = |a|+|b| = max(1,2|r|).

It follows that

    sup_U Q_L^U(h_r)
       = (6/pi)L[max(1,2|r|)-1].                          (2)

For r=0 this is zero. For r=1 it diverges as 6L/pi. These are the same point of the gauge quotient, and the physical A remains zero. This is not a divergent-A family.

For r=0 and U=X, the local field at site 1 cancels completely. Then E_L=-(L-1) and Q_L=-6L/pi. Since the norm in the preceding formula is nonnegative, this is also the infimum over U. Thus unrestricted twist infima fail even the stipulated gapped value.

[Proved scope] Equation (1) identifies a required repair, not a no-go theorem for symmetry twists. One may choose a fixed canonical representative, restrict to an appropriate symmetry-compatible gauge, or include an explicitly specified boundary counterterm. A genuine symmetry twist requires the symmetry action and its allowed representatives to be specified. Minimizing over arbitrary U without such a convention is not a quotient-invariant prescription. Canonical gauge fixing alone does not prove that the resulting seam has no nonuniversal defect energy.

## 5. Smooth weights and an exactly solvable gauge artifact

[Proved] For periodic weights w_j, set

    K_w(h) = sum_j w_j h_(j,j+1),
    delta_j = w_j-w_(j-1).

A coboundary changes this Hamiltonian by

    K_w(h+o tensor I-I tensor o)-K_w(h) = sum_j delta_j o_j. (3)

So nonuniform weighting is not itself gauge invariant.

Take a second gauge orbit of the same gapped product chain,

    g_r = -(Z tensor I+I tensor Z)/2
               + r(X tensor I-I tensor X).

Let b_j=(w_j+w_(j-1))/2>0 and W=sum_j w_j=sum_j b_j. Equation (3) gives independent one-site operators,

    K_w(g_r) = sum_j [-b_j Z_j+r delta_j X_j],
    E_w(g_r) = -sum_j sqrt(b_j^2+r^2 delta_j^2).

Since e=-1, the deficit is

    D_w(g_r) := eW-E_w(g_r)
             = sum_j [sqrt(b_j^2+r^2 delta_j^2)-b_j].       (4)

Let w_j=f(j/L), where f is a fixed smooth positive nonconstant function of period one. Taylor expansion, with f bounded away from zero, gives

    sum_j delta_j^2/b_j
       = L^-1 integral_0^1 f'(x)^2/f(x) dx + O(L^-2).

Expanding each square root in (4), with fixed r, therefore yields

    D_w(g_r) = r^2/(2L) integral_0^1 f'^2/f dx + O(L^-2).  (5)

The physical A is zero throughout this gauge orbit, while the coefficient in (5) is arbitrary. Thus the problem can occur at the same order as a proposed smooth conformal coefficient, and is not only a leading boundary-energy subtraction. Fixing r before taking L to infinity suffices for this conclusion.

## 6. A gapped counterexample in the canonical gauge

[Proved] Gauge fixing to the Hilbert-Schmidt orthogonal complement of G does not by itself remove the smooth response. Consider the spin-1/2 bond

    h_s = -(Z tensor I+I tensor Z)/2
               + s(X tensor Z-Z tensor X),
    0<|s|<1/2.                                           (6)

It has trace zero and equal partial traces, both equal to -Z. For any o,

    Tr[h_s(o tensor I-I tensor o)]
       = Tr[o(Tr_2 h_s-Tr_1 h_s)] = 0.

Hence (6) is the canonical Hilbert-Schmidt representative. Also

    1 <= M(h_s) <= ||h_s|| <= 1+2|s|.

The lower bound follows from the difference of expectations in |00> and |11>, as in Section 4. Thus its physical bond scale is finite and nonzero.

### 6.1 Unique ground state and a uniform gap

[Proved] Let n_j=(I-Z_j)/2 and N=sum_j n_j. The uniform ring Hamiltonian satisfies

    H_s+L = 2N+V,
    V = -2s sum_j [X_j n_(j+1)-n_j X_(j+1)].              (7)

The linear X terms cancel in the periodic sum. Since X_j commutes with n_k for j different from k, and -n_k<=X_j n_k<=n_k,

    V >= -4|s| N,
    H_s+L >= (2-4|s|)N.                                  (8)

The all-up product vector Omega has N Omega=0 and V Omega=0. It is therefore a ground vector of energy -L. On its orthogonal complement N>=1, so (8) proves uniqueness and a gap at least

    Delta >= 2-4|s| > 0                                   (9)

for every L>=3. The same positive-form estimate applies to finite excitations of the infinite product state. This is a gapped chain with e=-1 and physical A=0, not a nonconformal gapless endpoint.

### 6.2 An upper bound on the weighted deficit

[Proved] Let the positive periodic weights, b_j, delta_j, and W be as in Section 5. Now

    K_w+W = 2 sum_j b_j n_j + s sum_j delta_j X_j + V_w,
    V_w = -2s sum_j w_j[X_j n_(j+1)-n_j X_(j+1)].         (10)

Using the same operator bound on each term, and w_j>0,

    V_w >= -4|s| sum_j b_j n_j.

Set a_s=2-4|s|>0. The remaining one-site operators have lower eigenvalues

    lambda_j = a_s b_j/2
                   -sqrt((a_s b_j/2)^2+s^2 delta_j^2)
             >= -s^2 delta_j^2/(a_s b_j).

Here the inequality is sqrt(u^2+v^2)<=u+v^2/(2u) for u>0. Consequently

    E_w+W >= -s^2/a_s sum_j delta_j^2/b_j,
    D_w := -W-E_w <= s^2/a_s sum_j delta_j^2/b_j.         (11)

The unrotated product vector has expectation -W, proving D_w>=0.

### 6.3 A strictly positive lower bound of order 1/L

[Proved] Take w_j=f(j/L) as in Section 5. For all sufficiently large L, choose a product trial state with Bloch components

    x_j = -s delta_j/b_j,
    z_j = sqrt(1-x_j^2),
    y_j = 0.

Its energy is exactly

    E_trial = -sum_j b_j z_j
                +s sum_j w_j[x_j z_(j+1)-z_j x_(j+1)].    (12)

Uniformly in j, x_j=O(L^-1) and z_j=1-x_j^2/2+O(L^-4). The onsite contribution relative to -W is

    sum_j b_j(1-z_j)
       = (1/2)sum_j b_j x_j^2+O(L^-3).

The linear interaction contribution telescopes to

    s sum_j w_j(x_j-x_(j+1))
       = s sum_j delta_j x_j.

The remaining interaction terms contain at least three x factors and have total size O(L^-2). Substituting x_j yields

    E_trial+W = -s^2/2 sum_j delta_j^2/b_j+O(L^-2).

The variational principle gives E_w<=E_trial. Define

    I_f = integral_0^1 f'(x)^2/f(x) dx > 0.

Combining this trial bound with (11) proves

    (s^2/2) I_f <= liminf_(L to infinity) L D_w
               <= limsup_(L to infinity) L D_w
               <= [s^2/(2-4|s|)] I_f.                    (13)

In particular the smooth energy deficit is of order 1/L and has a strictly positive asymptotic lower coefficient, despite the uniform gap (9). No existence of the coefficient as an ordinary limit is needed to refute its identification with zero.

[Proved check] For f(x)=1+b cos(2pi x), 0<|b|<1,

    I_f = 4pi^2[1-sqrt(1-b^2)].

This follows by writing sin^2(theta)/(1+b cos(theta)) as a linear polynomial in cos(theta) plus a multiple of (1+b cos(theta))^-1 and integrating. The integral of (1+b cos(theta))^-1 over a period is 2pi/sqrt(1-b^2). For s=0.2 and b=0.4, the lower and upper coefficients in (13) are 0.065917004131 and 0.109861673552.

[Numerical] Exact matrices at L=4,6,8 satisfy the gap and both weighted variational bounds. The dimension never exceeds 256. Scalar evaluations of the product trial state at L=64,128,256,512 approach the lower coefficient. These computations verify signs and normalizations, not an ordinary-limit claim or a universal bound.

### 6.4 What this example does and does not exclude

[Proved] The observable eW-E_w in a fixed canonical gauge cannot generally be identified with the conformal vacuum coefficient just from its 1/L scaling. Equation (13) is a gapped lattice contribution at the same order. Subtracting the constant-profile ring Casimir energy does not help here, since that energy is exactly -L.

[Proved] This example does not exclude optimization over gauges. In fact it has an explicit repair. Set

    h'_s = h_s-s(X tensor I-I tensor X)
         = -(Z tensor I+I tensor Z)/2
               -2s[X tensor n-n tensor X].               (14)

For one bond,

    h'_s+I >= (1-2|s|)(n tensor I+I tensor n) >= 0,

and h'_s Omega=-Omega. Hence every positive weighted Hamiltonian built from (14), with either open or periodic support, has ground energy -W and D_w=0. The physical Hamiltonian is unchanged. The distinction is between the Hilbert-Schmidt canonical gauge and a gauge adapted to its product vacuum.

[Open] It remains possible that an appropriately optimized and normalized smooth-cutoff prescription has the desired physical limit. Equations (5) and (13) show why gauge covariance and the removal of gapped gradient response must be proved as part of such a prescription. The assertion that any smooth ground-energy coefficient is already A is false.

## 7. A gauge-invariant upper semicontinuous cutoff observable

[Proved] There is a finite-scale construction that handles gauge covariance without selecting a canonical representative. Let w_1,...,w_n be fixed nonnegative weights, with W=sum_j w_j>0, on an open segment of n+1 sites. Define

    D_w([h]) = inf_(o=o^dagger)
        { W e(h)-lambda_min(sum_(j=1)^n
                  w_j [h+o tensor I-I tensor o]_(j,j+1)) }. (15)

Scalar shifts cancel between the two terms. Replacing h by a coboundary merely translates the optimization variable o. Thus (15) is defined on the quotient.

A translation-invariant infinite-volume energy-minimizing state has bond expectation e(h), and its expectation of every coboundary is zero. Such a state can be obtained as a local weak limit of translation-averaged periodic ground states. Its restriction to the open segment is a trial state with weighted energy W e(h) for every o. Therefore every expression in braces is nonnegative.

Choose a norm-minimizing representative, including a scalar shift, with ||h||=M. Then |e(h)|<=M and the weighted Hamiltonian has norm at most WM. Evaluating at this representative gives

    0 <= D_w([h]) <= 2W M([h]).                            (16)

For fixed o the expression in braces is continuous in h, since e is Lipschitz and the finite matrix minimum eigenvalue is continuous. The infimum over o is upper semicontinuous, even though the optimization space is not compact. It is finite by (16). This upper semicontinuity also holds on the quotient, using any fixed continuous linear section of the finite-dimensional quotient. Positive scaling of h is accompanied by scaling o, so D_w is positively homogeneous.

For the gapped family (6), the representative (14) proves D_w=0 for every positive open weight profile. This is a useful exact check on the optimization.

[Proved conditional template] Let {w^(j),eta_j} be any nonempty fixed family of finite weight profiles and positive normalizations. Then

    B([h]) = inf_j eta_j D_(w^(j))([h])                   (17)

is finite, nonnegative, positively homogeneous, gauge invariant, and upper semicontinuous. For any chosen index j_0,

    B([h]) <= 2 eta_(j_0) W_(j_0) M([h]).                 (18)

Thus a proof that physical A<=B would yield a universal bound. A proof of equality on conformal bonds, together with the stipulated values elsewhere, would give more.

[Open] Neither inequality A<=B nor equality A=B has been established for a suitable fixed family. This is the unresolved substantive step. Selecting an arbitrary family in (17), including trivial choices, constructs many bounded semicontinuous functions but does not make any of them the physical Casimir coefficient. It has also not been proved that (15), suitably scaled, vanishes on every gapped bond after optimization over on-site coboundaries alone.

[Proved scope] Suppose a family of normalized observables R_L tends to A pointwise. Writing A as inf_L R_L additionally requires R_L>=A at every retained scale, or another justified upper-envelope construction. Otherwise one must take a tail limit, and the elementary infimum argument no longer proves upper semicontinuity. Section 3 demonstrates this issue for an actual vacuum Casimir family. Section 7 supplies a possible class of upper observables, not the missing comparison theorem.

## 8. What would finish the problem

[Proved, inherited from C107 Section 7] On the normalized conformal subset C, define the local upper envelope

    A^*(h) = inf_(r>0) sup {A(g): g in C, M(g-h)<r},

with empty supremum zero. This is an extended-real upper semicontinuous majorant. Its finiteness everywhere is equivalent to a finite global bound on C. Exact agreement with physical A on C requires relative upper semicontinuity in addition. At a gapped point, the requirement A^*=0 is stronger still unless a stable zero-A neighborhood has been proved.

[Proved] A finite assignment at nonconformal points does not establish these properties. Conversely, assigning the local upper envelope makes the topology automatic but can assign infinity. This is why a definition alone does not settle the bound.

[Proved conditional clarification, inherited from C107 Section 8] A Lieb-Robinson velocity ceiling for the physical scaling modes implies that divergent A requires divergent total central charge. It need not require any velocities to stay away from zero. A formal list of k sectors with c_a=1 and v_a=k^-1/2 has A=sqrt(k) while all velocities vanish. This list is not a lattice construction. A purported divergent lattice example must be checked through the velocity-weighted sum itself.

[Open] The remaining targets are (i) an actual comparison of physical A with finite gauge-invariant upper observables such as (17), or an independent proof of local boundedness at every quotient point, and (ii) an actual fixed-d bond sequence violating that local boundedness. Neither has been obtained here. There is no new proof of finite F(d), no physical counterexample to it, and no numerical evidence of divergence in this investigation.

## 9. Files and verification

[Work record] All changes are inside work_C110. PROGRESS.md records the work, RESULT.md gives the concise outcome, STATE.md is the local memory handoff, and checks.py generates CHECKS.md. No parent memory or shared report was modified. No external literature claim is needed for the new arguments, which are derived directly above. No network shell commands, heavy numerics, or cluster jobs were used.
