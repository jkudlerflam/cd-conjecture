# C60. A source-majorant obstruction and the scale-zero test

[proved, verdict with convention specified] Lemma N cannot hold with the interaction prefactor in C55 (30a), when the reference is the linear continuum theory, the crossover scale is fixed, and the finite-crossover term means the contribution of the finitely many smooth integrations above that scale. A free, off-half-filled nearest-neighbor chain already contradicts this version. This is a counterexample to a proposed correlation majorant, not to a Luttinger liquid or to the source-free endpoint lemma.

[proved, second result] The matched action satisfies the degree powers displayed at the start of the BM infrared construction after any fixed finite initialization. Time nonlocality at scale zero is allowed in that construction. The degree powers alone are not the hypothesis of an abstract theorem establishing Lemma N for replacement kernels.

[scope] The counterexample below uses the ordinary interpretation of a finite smooth crossover. If that phrase instead permits an arbitrary algebraically decaying function, or if the reference retains the full free lattice dispersion and point-split observable, the counterexample does not apply to that different assertion. Section 7 gives the distinction and repairs. No counterexample to the physical demonstrator is asserted.

## 1. Inputs and status

[imported] The microscopic action is C55's tuned nearest-neighbor demonstrator in hopping units, with the hopping ratio chosen before the physical port coupling. The trial filling belongs to a compact interval I inside (0,pi), away from pi/2. Massive matching and the exponential terminal-tree norm are imported from the finite-cell construction and TREE_NORM_EXTENSION. The physical order/degree hypothesis and its smallness reparameterization are imported from DEGREE_REPARAMETERIZATION. C55 Sections 2 through 4 supply the off-half-filled ultraviolet, compact-patch, and first-shell initialization. No microscopic matching result is rederived.

[read] C55 proves the complete Gaussian response of one initial endpoint, with a shell increment O(epsilon E_h^2). C56 leaves the same nonlinear comparison open. MARKED_COVARIANCE and REAL_CUTOFF_MARKS distinguish fixed analytic cutoffs from real velocity-dependent cutoffs. NORMALIZATION_RESOLVENT solves dependent ratios conditional on raw marked-tree bounds. COUNTERTERM_INDUCTION solves the sequence problem conditional on beta bounds. Those conditional implications are not promoted to unconditional RG results here.

[read, audit boundary] RETARDED_ENDPOINT_MAJORANT now claims a source-free endpoint theorem, including bounded coordinate loads, derivative carriers, and an anchor-incidence integration argument. Its stated marks are history/normalization marks and endpoint replacements satisfying its kernel hypothesis. Density and pair generating-function insertions are not among its proved conclusions. In particular, an unsuppressed bare observable is not an interaction endpoint with its epsilon-squared budget. This distinction is decisive for the counterexample.

[convention] All coefficients use one antisymmetric coefficient convention, with the degree factorial included in the pinned coefficient norm. Scalar source-free kernels are tracked separately. Claims marked proved below are the indicated statements or implications, not a proof of an interacting infrared limit. Files are confined to work_C60.

## 2. One fixed cutoff convention

[definition] Fix the trial Fermi point p and a real positive reference speed c. Choose a smooth even cutoff chi, equal to one on a neighborhood of zero, with support strictly inside the two disjoint branch patches. Choose its plateau and support so that

    chi_j(q)=chi(2^(-j)q),
    chi_j=1 on supp chi_(j-1),
    f_j=chi_j-chi_(j-1),                 j<=0.             (1)

The fixed speed c may enter the elliptic norm defining chi. It does not run with j or with a history coordinate. A trial-speed selection is a separate finite-dimensional real parameter problem. Write

    D_omega(q)=-i q_0+E_omega(q_1),
    E_omega(q_1)=omega c q_1+r_omega(q_1),
    r_omega(q_1)=O(q_1^2).                              (2)

This is the fixed-dispersion convention used for the BM comparison. A residual spatial derivative remains an interaction coordinate. It is not silently incorporated into a moving cutoff.

[proved, normalization identity] Suppose a quadratic absorption at step j is z_j D_omega. Starting from chi_j/(Z_j D_omega), the absorbed low covariance is

    chi_j/[Z_j D_omega(1+z_j chi_j)].

Set Z_(j-1)=Z_j(1+z_j), and subtract the retained covariance chi_(j-1)/(Z_(j-1)D_omega). Its shell, normalized by Z_(j-1), is exactly

    G_j=(1+z_j) f_j/[D_omega(1+z_j chi_j)].             (3)

Indeed, on supp chi_(j-1), chi_j=1, so both denominators agree there after multiplying by 1+z_j. Outside that support there is nothing to subtract. Absorbing z_j D rather than just its linear Taylor polynomial is an exact add-and-subtract operation. The compensating quadratic curvature term belongs to the irrelevant interaction. It cannot be discarded.

[proved, covariance consequences] On a fixed small complex disk for z_j, |1+z_j chi_j| is bounded below. D has magnitude comparable to |q| on each sufficiently small patch annulus. Momentum derivatives of (3) therefore have the usual shell bounds. A derivative in z_j differentiates only bounded rational factors, and a higher-scale ratio derivative differentiates the explicit product of 1+z's. No derivative of a moving smooth cutoff occurs. The determinant cofactor and triangular-resolvent arguments in the supplied notes apply to these rational factors. They still require a raw tree bound to construct the actual z_j.

[convention] For the counterexample z_j=0, all running interactions vanish, and (3) reduces to f_j/D. Thus it is independent of any unproved normalization equation. Work on the infinite spatial lattice at zero temperature. Equivalently, at every fixed infrared depth first remove L and beta, and only then send the depth to infinity. No last-shell estimate at arbitrary aspect ratio is used.

## 3. Degree powers survive the finite initialization

[proved conditional on the imported finite-integration theorem] Let the initial degree-2m action coefficients obey

    ||K_(2m)||_(H,q) <= A^(2m) epsilon^d_m,
    d_m >= max(1,m-1).                                  (4)

The norm is pinned, includes coefficient/label sums, and can have any fixed finite moment order supplied by the matching theorem. Local quadratic counterterms of order epsilon are allowed. The residual hierarchy in DEGREE_REPARAMETERIZATION is stronger, d_m=max(2,m-1). After a fixed finite sequence of the integrations in C55, the output obeys

    ||K'_(2m)||_(rho,q') <= C^(2m) epsilon^max(1,m-1),    (5)

with rho>0 fixed, q' fixed, and C independent of m and infrared depth. This is a finite-initialization assertion, not an iteration of a radius loss over infinitely many scales.

Proof. Freeze the physical coupling and all time rescalings. Introduce the artificial analytic family K_(2m)(z)=(z/epsilon)^d_m K_(2m). At a sufficiently small fixed radius |z|<=r, its coefficient norm at a suitably reduced field radius is uniformly small. C50's finite-integration bound gives an analytic output in a fixed output norm, uniformly on this disk. The finite number of integrations spends a finite number of norm reserves.

Consider a connected term with n initial vertices of degrees 2m_i and c contracted pairs. If its external degree is 2m, then

    m=sum_i m_i-c,
    c>=n-1.                                            (6)

Connectedness gives the second inequality. The power of z is Q=sum_i d_(m_i). Since d_(m_i)>=m_i-1 and d_(m_i)>=1,

    Q>=sum_i m_i-n=m+c-n>=m-1,
    Q>=n>=1.                                           (7)

Thus the degree-2m output has no Taylor coefficient below r_m=max(1,m-1). If its uniform norm on |z|=r is at most B rho^(-2m), Cauchy and the vanishing orders give, for epsilon<=r/2,

    ||K'_(2m)(epsilon)||
       <= B rho^(-2m) sum_(Q>=r_m)(epsilon/r)^Q
       <= 2B rho^(-2m) r^(-r_m) epsilon^r_m.            (8)

Since r_m<=m for m>=1, the fixed constants in (8) are bounded by C^(2m). This proves (5). The argument works for finite coefficient measures and does not assert analyticity in the physical port coupling through zero after the bandwidth-dependent time rescaling.

[proved, moment qualification] A compact smooth Fourier crossover does not preserve literal exponential decay in position. Use the fixed polynomial moments supplied by C55, or retain the finite ultraviolet integration as part of the endpoint's bounding representation. The degree proof is valid in either available finite-integration norm. It does not upgrade a polynomial moment bound to an exponential one.

[proved, Fourier consequence] The absolute value of a pinned kernel's Fourier transform is at most its total variation. Equation (5) therefore gives the degree-graded Fourier bound with exponent max(1,m-1). Finitely many fixed derivatives in external momenta follow from the corresponding fixed moments. Symmetries and momentum support leave the localized quadratic and single forward-quartic terms described in C55. These conclusions do not use equal-time support of K.

## 4. What the BM shortcut establishes

[read, primary-source check] BM [0907.2837, Section 3, (63)-(65)](https://arxiv.org/html/0907.2837) starts its infrared construction with all-degree kernels in independent spacetime arguments. Time nonlocality is permitted there. The displayed Fourier degree estimate agrees with (5). The subsequent analyticity assertion invokes an earlier tree theorem, rather than stating a replacement-kernel theorem with that Fourier bound as its sole hypothesis.

[read, primary-source check] In [cond-mat/0004226, Theorem 3.12 and (3.83), (3.101)-(3.105)](https://arxiv.org/pdf/cond-mat/0004226), the proof keeps coordinate zeros, derivative carriers, weighted position integrals, and scale gains together. A bound on Fourier amplitudes alone does not supply these estimates. The endpoint types in that theorem are those of its constructed expansion.

[proved, why the distinction is substantive] For example, a translation-invariant retarded quartic joining two different branch densities at times separated by T can have coefficient a and Fourier dependence a cos(T q_0). Its Fourier supremum is at most |a| independently of T. Its second derivative at q_0=0 has magnitude |a|T^2. Its first absolute position moment is of order |a|T. Adding its time-reflected copy preserves the reality symmetry. It conserves charge, and uses one right and one left density, so Grassmann antisymmetry does not make it vanish.

Consequently a Fourier degree bound cannot provide uniform Taylor remainder constants. This example is only a counterexample to that norm implication. The actual matched kernels have a stronger locality norm that excludes an arbitrary T at fixed amplitude. Section 3 verifies the grade hierarchy for them; it does not establish the all-scale localization estimate by invoking (65) alone.

[audit boundary] The newest endpoint note addresses this stronger norm. Its edge estimate product_e r_e^b_e <= C^p exp(T_tree) pays loads once they are bounded. Its incidence argument controls anchor integrations once interpolated line arguments have been assigned. The additional issue is simultaneous allocation of scale gains to all nested remainders, including the product-rule terms when a derivative acts on an already present coordinate zero. I did not establish or refute its complete source-free assertion in this run. No failure of its actual strong endpoint norm is claimed.

## 5. An explicit free counterexample to the interaction-prefactored source bound

[definition] Consider the normalized free chain

    e_p(k)=2(cos p-cos k),             0<p<pi,

with the sea |k|<p filled. Take p=pi/3, which is away from half filling and the band edges. At distinct equal-time sites,

    F(x)=<f_x^dagger f_0>=sin(px)/(pi x).                (9)

The same real even function gives <f_x^dagger f_y>=F(x-y). For x!=y the complementary contraction is <f_x f_y^dagger>=-F(x-y).

Use the charged nearest-neighbor pair observable appearing in C55's source-overlap calculation,

    P_x=f_x f_(x+1),
    P_0^dagger=f_1^dagger f_0^dagger.

[proved, exact correlation] For every integer x>=2, Wick's rule gives

    C_pair(x)=<P_x P_0^dagger>
             =F(x)^2-F(x-1)F(x+1).                     (10)

The sign follows directly from <ABCD>=-<AC><BD>+<AD><BC> when <AB>=<CD>=0. Both complementary contractions carry a minus sign, whose product is positive.

Use

    sin[p(x-1)] sin[p(x+1)]=sin^2(px)-sin^2 p.

Substitution in (10) yields

    C_pair(x)
      =sin^2 p/[pi^2(x^2-1)]
           -sin^2(px)/[pi^2 x^2(x^2-1)]
      =sin^2 p/(pi^2 x^2)+D_p(x),

    D_p(x)=[sin^2 p-sin^2(px)]/[pi^2 x^2(x^2-1)].       (11)

This is an exact identity, not a large-distance fit.

[proved, amplitude matching] The free continuum reference has X_pair=1. Expanding the lattice pair in its two chiral fields gives the local coefficient -2i sin p, as in C55 (41). Each chiral propagator has coefficient 1/(2pi), so the matched leading pair correlation is

    C_pair^ref(x)=sin^2 p/(pi^2 x^2).                   (12)

Alternatively (12) follows without any chiral normalization convention from lim_(x->infinity) x^2 C_pair(x)=sin^2 p/pi^2. A constant change of source amplitude cannot remove D_p while retaining this limit.

[proved, nonzero algebraic remainder] At p=pi/3 and x=3m,

    D_(pi/3)(x)=3/[4pi^2 x^2(x^2-1)],
    lim_(m->infinity) x^4 D_(pi/3)(x)=3/(4pi^2)>0.      (13)

In particular D is not rapidly decreasing. For any M>4, |x|^M |D(x)| is unbounded on this subsequence.

[proved, even the zero-momentum pair channel has this issue] Writing sin^2(px)=(1-cos(2px))/2 separates

    D_p(x)=[sin^2 p-1/2]/[pi^2 x^2(x^2-1)]
             +cos(2px)/[2pi^2 x^2(x^2-1)].             (14)

At p=pi/3 the nonoscillating numerator is 1/4. Thus removal of the subleading oscillatory channel alone does not eliminate the free x^(-4) correction. The obstruction is also a descendant contribution within the leading pair channel.

## 6. Why a finite smooth crossover cannot absorb (13)

[proved] Fix any smooth momentum partition with patches equal to one near the two Fermi points. For equal-time correlations, the high part of the filled-sea multiplier is smooth and periodic in k. Its factor vanishes on neighborhoods of the only jumps, at +/-p. Repeated integration by parts therefore gives, for every integer M,

    |F_high(x)|<=C_M(1+|x|)^(-M).                        (15)

Write F=F_low+F_high in (10). Every term in C_pair-C_pair,low contains at least one F_high at x or x+/-1. The other factor is bounded. Thus that difference satisfies (15) for every M. The corresponding statement holds for any finite number of smooth shells above a fixed positive scale. Smooth continuum cutoffs away from the two chiral singularities have the same property.

It follows that changing the finite crossover changes (11) only by a rapidly decreasing sequence. The low-energy remainder still has the nonzero subsequential limit (13). This also follows directly from the low-frequency nonanalyticity associated with that algebraic tail; no ultraviolet contact term can cancel it.

[proved, contradiction to the proposed bound in its finite-crossover meaning] At zero interaction the small action parameter in C55 is epsilon=0. Curvature of the fixed nearest-neighbor free dispersion and point splitting of P are still present. Suppose C55 (30a) held for this continuum comparison:

    D_p(x)=F_cross(x)+sum_(h<=0) R_h(x),

    |R_h(x)|<=C epsilon E_h^(2+alpha)
                       /[1+(E_h |x|)^M].              (16)

The actual finite-crossover contribution F_cross satisfies (15). With epsilon=0, (16) forces every R_h to vanish. It would make D_p rapidly decreasing, contradicting (13). The contradiction occurs in one channel and at an allowed filling. No determinant estimate, counterterm choice, normalization convention, or history mark is involved.

[qualification on the wording of C55] If the only condition imposed on F_cross were an arbitrary decay bound with one chosen M<4, then D_p itself could be renamed F_cross. The displayed inequality in that purely existential interpretation would not be contradicted by (13). Such a function is not generated by a finite smooth crossover, as (15) proves. For any version requiring M>4, even that weaker interpretation is contradicted. This qualification matters: the counterexample refutes the RG finite-crossover majorant, and does not manufacture a contradiction by imposing a different meaning on an unspecified remainder.

[scope] This is a uniform perturbative-neighborhood obstruction. The zero-input action belongs to the Gaussian boundary of the proposed small-kernel class. It is not asserted that the exact port model at a particular nonzero coupling is free. This boundary is compatible with the tuned schedule: keep the cell fixed, take Jz=1, Jp=s, and the physical port coupling lambda=s^2. Then |t| is proportional to s, the normalized static interaction is O(s), and the matching error bound lambda/|t| is O(s). The normalized free dispersion remains the one used above. Allowing constants to diverge as the normalized interaction tends to zero would avoid the test but would not be the uniform small-coupling lemma used in the stated induction.

[definition of the crossover term used in the test] A scalar two-source contribution produced above the fixed starting scale, or a finite sum of correlation shells at fixed positive energies, is a finite-crossover contribution here. An irrelevant one-source vertex generated there and subsequently propagated through arbitrarily low scales belongs to the infrared remainder. Including the latter in an unrestricted “crossover” would conceal the very bound being tested.

## 7. Repairs and the remaining nonlinear task

[proved, necessary distinction] Bare source Taylor remainders and fixed lattice curvature need not carry an interaction factor. At the free point they are O(1) in a fixed crossover convention. This is consistent with the supplied IR_SOURCE_REVIEW, which separates bare source remainders from interaction-induced source corrections. The strong residual endpoint budget in RETARDED_ENDPOINT_MAJORANT cannot be assigned to them.

[possible corrected statements, not proved here] There are three different ways to avoid the contradiction.

1. Replace the coefficient C epsilon in the full correlation remainder by C(A_free+epsilon), where A_free controls the fixed curvature and bare-source irrelevant terms and need not vanish at the free point.
2. Subtract the full free lattice-to-reference correlation difference, including its algebraic descendants, before requiring an interaction-proportional bound. Any extension of this subtraction to interacting anomalous powers must itself be proved. The free calculation alone does not prove that extension.
3. Compare first with an instantaneous lattice reference retaining the same free dispersion and the same point-split source. At zero interaction the two correlations then agree. Extracting that reference's continuum asymptotics is a separate step with its own free subleading terms.

One could instead make the crossover scale depend on the interaction and budget curvature as part of the small parameter. That changes the uniformity and cutoff convention. It is not the fixed-crossover statement tested here.

[unresolved] None of these changes proves the corrected nonlinear endpoint and source comparison. In particular, the source-free bounded dimension d_R(p)>=1 does not apply unchanged to a source cluster: with one bilinear source, the two-fermion source vertex is marginal and must be normalized; with two sources, their separation must be retained. Treating source differentiation as selection of an interaction endpoint would erase this distinction.

[precise source-free routing target] To make the remaining localization check concrete, use a reduced forest whose edges have scale lengths ell_v=h_v-h_(v')>0, with v' the parent toward the root. Let z_v=2 for a quadratic remainder, 1 for a quartic remainder, and 0 otherwise. A legal routing term A has derivative order a_f on field f, evaluated at its first contraction scale s(f), and coordinate order b_l on a covariance line of scale s(l). A derivative on a root external field is evaluated at the retained root scale. Endpoint-edge moments are paid at the fixed initial scale, set to zero in these units. The needed paired scale allocation is

    sum_(legal A) |c_A|
      2^[sum_f a_f s(f)-sum_l b_l s(l)]
       <= C^(sum_i p_i+n)
            product_v 2^(-z_v ell_v).                 (17)

The sum includes the product-rule terms in which derivatives consume coordinate zeros. It must be realized by the same routing for which the endpoint loads and shifted anchor integrals are bounded. A fixed loss C per endpoint field is allowed; C per unit scale depth is not. On uncompressed forests ell_v=1. Existing derivative factors in running spatial-quadratic endpoints are counted in their engineering dimensions, not counted again as localization derivatives in (17).

[conditional implication] Once (17) is established for the proposed arbitrary endpoints, ordinary source-free power counting gives the edge exponent

    d_R(p)=p/2-2+z(p),
    d_R(2)=d_R(4)=1,          d_R(p)=p/2-2 for p>=6.

The degree-uniform ratio bound in MARKED_COVARIANCE spends only a fixed fraction of d_R. Positive reserves then sum fields and scales and supply short memory to an initial endpoint or covariance mark. NORMALIZATION_RESOLVENT preserves that decay after eliminating z. C55's linear endpoint theorem checks the single-initial-endpoint sector and does not prove (17) for multiple vertices. The pure-forward reference cancellation concerns the summed reference beta, not an absolute bound on each reference graph.

[scope of (17)] Equation (17) is a focused diagnostic of the source-free assembly, not a claim that it alone proves N3. Full Lemma N also needs its source version and a corrected correlation remainder such as one of the three alternatives above. Since the requested counterexample alternative is met by (11)-(16), no claim of a single sufficient unresolved inequality for an otherwise proved Lemma N is made.

## 8. Checks and outcome

[proved, consistency checks] Charge conservation makes <P_x>=0, so no disconnected term was omitted. Reflection makes F even and fixes the signs in (10). The leading coefficient in (12) agrees with C55's nonzero pair overlap. The choice p=pi/3 satisfies all stated filling exclusions. The x^(-4) term is subleading to the free x^(-2) power, so it leaves X_pair=1 unchanged. It is a correction-to-scaling obstruction to an epsilon-prefactored remainder, not a new marginal operator or instability. The argument does not use half filling, a band edge, or a finite-volume zero mode.

[numerical, supplementary] checks.py compares the Wick determinant, the closed expression (11), and the limit (13), and checks the connected degree exponent inequality (7) on finite examples. Its output is in checks.json. These checks are controls of exact formulas, not simulations of the RG or evidence for the unresolved forest inequality.

[final status] A complete proof of Lemma N as written is not supplied, because its continuum source-majorant version has the explicit free obstruction above. The degree-graded scale-zero route passes its algebraic hierarchy and time-nonlocality tests but does not by itself establish the localization and source comparison. The nonlinear strong endpoint majorant is neither proved anew nor refuted here. No Luttinger theorem, thermal coefficient, or central-charge theorem for the port model is claimed.
