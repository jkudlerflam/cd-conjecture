# C51. Central charge after a multichannel Luttinger limit

## 1. Result, assumptions and provenance

[proved conditional on G below] Put N=n+1. A complete compact Gaussian infrared theory of N canonical bosons, with any strictly positive quadratic Hamiltonian including density-current mixing, has

    c_R = c_L = N,                 c_cd = (c_R+c_L)/2 = N.       (1)

Here c_cd is the convention in which one nonchiral compact boson has c=1. It is the total central charge of the collection of modes, as used in the local statement of the c-d conjecture. It is not c_R+c_L=2N. The central terms do not depend on the velocities, compactification radii or finite-size background charges.

[proved conditional on G] If the positive speeds are u_(a,R),u_(a,L), the Gaussian thermodynamic free energy in physical coordinates is

    f_G(T) = e_G - A T^2,
    A = (pi/12) sum_(a=1)^N [1/u_(a,R) + 1/u_(a,L)].           (2)

If u_(a,R)=u_(a,L)=v_a, this is the requested A=(pi/6) sum_a 1/v_a. There is no imposed symmetry in the two-bus construction that establishes these pairwise equalities. A thermal coefficient is velocity weighted. Equation (2) alone is not a measurement of the unweighted number (1).

[assumed, G, complete Gaussian scaling limit] The following is the exact meaning of a complete Luttinger limit used for (1).

1. There are N compact canonical bosons. Their 2N equal-time current modes I_m obey

       [I_m^A,I_l^B] = m S_AB delta_(m+l,0),
       (I_m^A)^dagger=I_(-m)^A,        S=diag(1_N,-1_N).

   A positive level matrix can first be normalized to this form. The minus signs describe left movers with the same spatial Fourier convention as the right movers, not negative-norm states.
2. The limiting Hamiltonian is the normal-ordered quadratic form with a real symmetric positive definite 2N-by-2N matrix V. It acts on the positive-energy compact-boson representation. Momentum has the usual signed quadratic form S. The matrix V and all velocities are in physical length units.
3. The limiting Hilbert space and local field algebra are the compact-boson ones, including their charged vertex sectors. There is no further massless tensor factor or stress tensor outside them. The microscopic identification covers this full infrared theory, rather than only a selected collection of source correlations.

[assumed, finite-ring version G_ring] When ring spectra are discussed, the compact zero modes form a full-rank lattice, or a fixed finite union of lattice cosets with fixed finite multiplicities. After canonical changes of basis their chiral charges are

    x = B(m+alpha_L),       m in Z^(2N),                       (3)

with B invertible, and alpha_L taken modulo Z^(2N). Finite-index parity restrictions can be absorbed in B and the finite coset list. Each allowed charge sector carries its stated full oscillator Fock space. More general oscillator projections need their own character calculation and are not silently included in (3). Offsets can depend on L and need not converge. The finite-volume background is measured relative to fixed equilibrium densities and any chemical potentials already in H.

[proved, why this is the compact-boson lattice] For a field valued in R^N/(2pi Lambda), winding around a spatial circle lies in Lambda. Wave functions of its constant coordinate are periodic on that torus, so their Fourier momenta lie in the dual lattice Lambda*. Thus momentum and winding supply rank 2N discrete data. An invertible real canonical transformation sends these data to a full-rank lattice of chiral charges. Constant background fields shift it. This explains (3) for compact bosons; it does not determine the microscopic wrapper's particular B or selection rules.

[assumed, distinct microscopic transfers] There are two additional, alternative contracts in this note. F in Section 5 transfers finite-size spectra and multiplicities with errors o(1/L). P in Section 7 transfers the thermal coefficient by controlling scalar RG terms. Neither is built into a statement only about decay exponents and beta functions. G identifies the central charge of the complete scaling theory without needing either lattice entropy convergence or a thermal-pressure theorem.

[proved, status of C50 at the final cross-read] C50 PROGRESS.md and PROOF.md appeared during this work and have been read. They prove finite-scale integration, localization and conditional recurrence/source estimates, and explicitly leave the all-scale Luttinger theorem open. Their finite-scale coefficient norm retains a scalar term, but no temperature-sensitive scalar bound of the form P is established. They do not establish G or F either. This note assumes the requested eventual Gaussian conclusion where stated. The completeness clause of G and the scalar estimate P remain substantive requirements, not consequences of calling the current partial output "Luttinger".

[proved, relation to prior work] C46 PROOF Sections 7 and 8 already give (2) for a diagonal Gaussian reference, its simple-root free-fermion counterpart, a common-bond quotient-pressure theorem, and the one-current Sugawara calculation. Those conclusions are inherited, not claimed as new. The additions here are the general canonical diagonalization, a uniform treatment of shifted compact sectors, a finite-size multiplicity and transfer statement, and a sufficient scale-by-scale scalar estimate with its full summation proof. All proofs below are algebraic or analytic and use no external theorem beyond the explicitly assumed infrared inputs.

## 2. Positive Gaussian forms and arbitrary velocity ratios

[proved, diagonalization lemma] For V>0 let

    W = V^(1/2) S V^(1/2).

W is real symmetric and nonsingular. Congruence with S shows that it has N positive and N negative eigenvalues. Choose an orthogonal O such that

    O^T W O = D = diag(u_(1,R),...,u_(N,R),
                      -u_(1,L),...,-u_(N,L)),    all u>0.

Put U=|D| and

    R = V^(-1/2) O U^(1/2).                                  (4)

Then

    R^T V R = U,              R^T S R = S.                    (5)

The first identity follows by direct multiplication. For the second, W^(-1)=V^(-1/2) S V^(-1/2), so

    R^T S R = U^(1/2) D^(-1) U^(1/2) = S.

Also S V R = R D. Thus the eigenvalues are the signed propagation speeds, not eigenvalues of V read in a noncanonical basis. This proof uses positivity only. Distinct speeds, comparable speeds and a lower bound uniform in N are unnecessary.

[proved, currents and positivity] Set K_m=R^(-1)I_m. Equations (5) imply R^(-1) S R^(-T)=S, so the K currents have the same commutators. Define

    J_m^a = K_m^a,             bar J_m^a = K_(-m)^(N+a).

Both now obey [J_m^a,J_l^b]=m delta_ab delta_(m+l,0), and all barred and unbarred algebras commute. In a sector of fixed real charges x_i, positive modes annihilate its cyclic vector |x>. The normalized creation operators are J_(-m)/sqrt(m). The quadratic Hamiltonian is a sum of harmonic oscillators with frequencies 2pi u_i m/L, together with its zero-mode quadratic form. Negative left-moving Fourier signature has become positive left-moving energy.

[proved, effect on compactness] The same R^(-1) sends the original charge lattice to the chiral charge lattice. It need not send it to a product of N two-dimensional lattices. The zero-mode sectors may therefore couple different channels even though the oscillator Hamiltonian is diagonal. Factorization of the full compact theory into N independently compactified scalar Hilbert spaces is not required for the sum of central terms.

[proved, unequal right and left speeds] For any chosen pair write

    v_a=(u_(a,R)+u_(a,L))/2,
    w_a=(u_(a,R)-u_(a,L))/2.

Its Hamiltonian is proportional to v_a(L_0+bar L_0-1/12)+w_a(L_0-bar L_0), with v_a>|w_a|. The second term is a drift. One common change of reference frame cannot remove arbitrary different w_a. Even if all w_a vanish, distinct v_a prevent a single spacetime conformal coordinate for the full dynamics. The collection still has commuting Virasoro algebras and a well-defined sum (1). This is the multi-velocity definition of c used here, not an assertion of one common Lorentz symmetry.

## 3. Central terms, including nonzero charge sectors

[proved, algebra on a dense domain] For a single positive current algebra define, on vectors with finitely many oscillator excitations,

    L_m = (1/2) sum_(r in Z) :J_(m-r) J_r:,
    L_0 = x^2/2 + sum_(r>0) J_(-r) J_r.                      (6)

Normal ordering puts negative modes to the left. The sums have only finitely many nonzero actions on each such vector. For m>0 this reads

    L_m = x J_m + (1/2) sum_(r=1)^(m-1) J_(m-r)J_r
                      + sum_(r>=1) J_(-r)J_(m+r),
    L_(-m)=L_m^dagger.                                       (7)

Commuting either current in (6) with J_l gives

    [L_m,J_l]=-l J_(m+l).                                    (8)

It follows by the Jacobi identity that C_(m,l)=[L_m,L_l]-(m-l)L_(m+l) commutes with every current. On a fixed-charge Fock space any such operator is scalar: its value on the cyclic vector is proportional to that vector because positive currents annihilate it, and commutation with all creation operators fixes its value everywhere. Its L_0 degree is -(m+l), so this scalar vanishes for m+l!=0.

[proved, the remaining scalar without assuming a zero-charge ring vacuum] For m>0, (7) gives

    L_(-m)|x> = x J_(-m)|x>
                   + (1/2) sum_(r=1)^(m-1) J_(-r)J_(-(m-r))|x>.

The one- and two-oscillator terms are orthogonal. Contraction of the latter, including the coincident-mode case, gives

    ||L_(-m)|x>||^2 = m x^2 + (1/2) sum_(r=1)^(m-1) r(m-r)
                    = m x^2 + (m^3-m)/12.                  (9)

Here sum r(m-r)=m sum r-sum r^2=(m^3-m)/6. Since L_m|x>=0 and 2m<L_0>=m x^2, subtraction of 2m L_0 removes the charge contribution. Thus on every charge sector

    [L_m,L_l]=(m-l)L_(m+l)
                   +(m^3-m) delta_(m+l,0)/12.              (10)

This proves c=1 even when the finite ring has no zero-charge ground sector. It also proves unitarity on the Fock representation.

[proved, total central charge] Apply (6)-(10) independently to all N right and N left currents. Cross commutators vanish. Therefore L_m^tot=sum_a L_m^a and bar L_m^tot=sum_a bar L_m^a have c_R=c_L=N. Their central terms are insensitive to the zero-mode lattice and its shifts. Under the completeness clause of G these are the full infrared central charges, proving (1). Without that clause this calculation identifies only the Gaussian current subtheory.

[proved, cylinder constant in the Gaussian convention] The oscillator zero-point energy of one chirality, regulated by exp(-epsilon m), is

    (pi u/L) sum_(m>=1) m exp(-epsilon m)
       = (pi u/L)[epsilon^(-2)-1/12+O(epsilon^2)].

To verify the series, sum the geometric series and differentiate it to obtain exp(-epsilon)/(1-exp(-epsilon))^2, then expand at zero. With a fixed physical momentum cutoff, epsilon=2pi a/L, the divergent term is proportional to L and is absorbed into e_G L. The remaining term is -pi u/(12L), or (2pi/L)(-u/24). This fixes the Gaussian cylinder convention. Showing that the microscopic finite-size scalar agrees with it is part of F, not a consequence of normal ordering a lattice Hamiltonian.

[proved, what is not an irreducible Virasoro multiplicity claim] The oscillator character of a charge sector is the Heisenberg Fock character. It is not always the character of one irreducible Virasoro module. For example, at x=0, L_(-1)|0>=0 while J_(-1)|0> has positive norm, so the level-one current state is not a Virasoro descendant of the vacuum. The complete Fock multiplicities below are proved directly. An additional decomposition into irreducible Virasoro highest-weight modules is not needed for (1) and is not asserted without a separate representation-theoretic proof.

## 4. Shifted sectors, oscillator multiplicities and irrational density

[proved under G_ring, exact Gaussian spectrum] Enumerate the 2N positive speeds by i and put

    Q_alpha(m) = (1/2) [B(m+alpha)]^T U [B(m+alpha)],
    C_U = (1/24) sum_i u_i.

For each allowed charge sector and each finite set of occupations n_(i,r)>=0, r>=1,

    E_G = e_G L + (2pi/L)
                  [Q_(alpha_L)(m) + sum_(i,r) u_i r n_(i,r) - C_U].    (11)

Let k_i=sum_r r n_(i,r). The multiplicity at fixed m and fixed vector (k_1,...,k_(2N)) is

    product_i p(k_i),
    sum_(k>=0) p(k) z^k = product_(r>=1) (1-z^r)^(-1).       (12)

Proof: for each oscillator r its occupation contributes 1+z^r+z^(2r)+..., and multiplying counts partitions of k. Different currents have independent creation operators. If weighted energies coincide, sum these multiplicities over all coincident labels and over any repeated coset sectors. Unequal or rationally independent speeds generally destroy degeneracies between different k vectors, but never remove the states counted in (12).

[proved, bounded zero-mode rounding] Let A=B^T U B/2>0 and alpha range over [0,1]^(2N). Choosing m=0 gives Q_alpha(0)<=C for a fixed C. Therefore min_m Q_alpha(m)<=C uniformly. If m minimizes, lambda_min(A)|m+alpha|^2<=C, so all minimizers lie in one bounded set after reducing alpha to this cube. Rounding costs O(1/L), not a volume-independent gap. The ground energy is

    E_(G,0)(L) = e_G L + (2pi/L)
                      [min_m Q_(alpha_L)(m) - C_U].          (13)

There need not be one minimizer, either at isolated sizes or at limiting values of the offset.

[proved, source of the moving offset] Suppose an equilibrium charge density is rho and the finite-volume charge is an integer q. The difference q-L rho equals m+alpha_L after subtracting an integer part of L rho. Several charges, background fluxes and parity conditions give the vector version of this construction. Irrational rho does not require a sequence of exactly matching rational densities. Alpha_L stays on a compact torus, and the bounded rounding estimate applies to every L. For the wrapper, M is an integer and L=bM+r. At an interior equilibrium density nu_0, choosing M nearest L/(b+nu_0) gives

    M-L/(b+nu_0)=O(1),
    r/M-nu_0=O(1/L).

For all sufficiently large L, this choice satisfies 0<=r<=M because the interior density leaves an order-L reserve at either endpoint. Thus there is no arithmetic obstruction to the length sector. Its curvature and its couplings to the other charges determine the entries of Q, rather than introducing an additional oscillator.

[proved, why a raw Casimir fit can fail] Equation (13) contains an L-dependent term of the same order as the Casimir term. It cannot be included in an o(1/L) error. A one-coordinate slice Q_alpha(m)=a(m+alpha)^2 already gives min Q=a dist(alpha,Z)^2. At alpha_L equal to the fractional part of L rho for irrational rho, this is bounded but nonconvergent. The nonconvergence follows without an equidistribution assumption: irrational rotations are dense. One elementary proof chooses, by the pigeonhole principle, a nonzero multiple of rho modulo one arbitrarily close to zero; its successive additions approximate any target point, and arbitrarily late orbit tails have the same property. Hence subsequences approach offsets 0 and 1/2, giving different limits. A compact model with additional unshifted coordinates has the same example. The appropriately corrected quantity is

    L[E_(G,0)(L)-e_G L] - 2pi min_m Q_(alpha_L)(m)
            = -(pi/12) sum_i u_i.                           (14)

A naive effective c extracted before the zero-mode subtraction need not converge.

[proved, precise subsequential spectrum statement] Every sequence of L has a subsequence on which alpha_L converges on the torus. Use consistent representatives along that subsequence; changing representatives only relabels m. Below any fixed bound on L(E-e_G L), positivity bounds |m+alpha_L| and all k_i uniformly. There are finitely many such labels, since u_min>0 and p(k)<infinity. Each energy in (11) is continuous in alpha. Thus the ordered levels, including multiplicities, converge along the subsequence, away from any chosen cutoff level that belongs to the limiting spectrum. Levels can cross, so no fixed label for the unique ground sector or uniform gap between competing zero sectors is claimed. Formula (11), with alpha_L retained, is the all-ring assertion.

[proved, interpretation as towers] In each moving charge sector the commuting Virasoro actions in Section 3 and the full oscillator counts (12) give the Gaussian conformal towers. They are current-algebra towers which carry Virasoro representations. With several speeds the physical energies are weighted sums of the separate conformal levels, not multiples of one common velocity. For a local field with weight h_i in each chiral coordinate, those coordinates are x-u_(a,R)t and x+u_(a,L)t. There is no reason to infer a single spacetime Virasoro symmetry of the common Hamiltonian from the equal-time sum of the algebras.

## 5. What transfers the tower statement to a lattice

[assumed, F, finite-size transfer] Fix all microscopic couplings and N. For every finite rescaled-energy cutoff K, suppose the microscopic low-energy space can be identified with the Gaussian space (11) below a slightly larger cutoff, including multiplicities and the actual physical sewing rules. Require an error bounded by epsilon_K(L)/L with epsilon_K(L)->0 on this space, including the scalar energy term after subtracting the actual bulk energy e_0 L. Require also that there are no additional microscopic states below the target cutoff. The latter condition is a spectral exhaustion requirement. Boundary cutoffs can be chosen away from a limiting eigenvalue, or one can state the identification for ordered finite collections of levels with a buffer.

[proved, transfer theorem] Under F, every fixed collection of low-lying eigenvalues is given by (11) with e_G replaced by e_0 and error o(1/L), and has the multiplicities in (12) in the limit. To prove this, compare the two finite matrices on the identified buffered space. For self-adjoint A and E with ||E||<=epsilon, every Rayleigh quotient changes by at most epsilon. Applying the minimum-over-subspaces characterization of each ordered eigenvalue gives |lambda_j(A+E)-lambda_j(A)|<=epsilon. Multiplying the F bound by L gives convergence. Spectral exhaustion prevents missing levels in the comparison. It is not legitimate to infer this exhaustion from a few trial states. At an exact Gaussian degeneracy the microscopic levels may split by o(1/L); their combined limiting multiplicity is fixed.

[proved, emergent rather than exact charge labels] The individual Gaussian charges in (3) need not be exactly conserved by the microscopic Hamiltonian. Exact total payload U(1) and marker number do not imply separate U(1) symmetries for every normal mode. Under F, irrelevant terms may mix the asymptotic charge sectors and split crossings by o(1/L). The spectral comparison still gives the limiting multiplicities. No claim of an exact finite-lattice decomposition into all Gaussian charge sectors is being made.

[proved, a weaker bound is insufficient] If the only microscopic estimate is an interaction norm ||R||_loc<=epsilon at fixed couplings, its total operator-norm bound is generally of order epsilon L. Multiplication by L makes the bound worse. C45's fixed local remainder estimate therefore does not prove F. Likewise, an O(1) seam bound proves equality of pressures in the thermodynamic limit but not an o(1/L) spectral comparison.

[assumed, application boundary for the wrapper] C44 gives identity transport and an explicit cyclic gluing spectrum for its common-bond reference. C46 proves the reference quotient changes log Z by at most O(beta)+O(1). Those are sufficient for that reference pressure. Extending (3), (11) and F to the full interacting finite wrapper requires a finite-volume identification of the dressed charge lattice and the sewing of oscillator states. We do not replace that work by a statement that any projection preserves every tower. Fixing one total physical momentum, for example, can discard oscillator levels even when the full physical Hilbert space retains them. The full thermal trace here sums all momenta.

## 6. Thermodynamic limit and zero modes

[proved, uniform theta bound] Put d=2N and t=2pi beta/L. For a single lattice coset define

    Theta_alpha(t)=sum_(m in Z^d) exp[-t Q_alpha(m)].

There are constants C_0,C_1,c>0, independent of alpha modulo integers and of t, such that

    exp(-C_0 t) <= Theta_alpha(t)
                  <= C_1 (1+t^(-1/2))^d.                  (15)

For the lower bound retain m=0 with alpha in a bounded fundamental cell. For the upper bound use Q_alpha(m)>=c|m+alpha|^2 and factor the sum into d one-dimensional sums. For any shifted lattice, a nonnegative unimodal function has sum bounded by twice its maximum plus its integral, by comparing the terms on each side of its maximum to adjacent integrals. Hence sum_(k in Z) exp[-ct(k+a)^2]<=2+sqrt(pi/(ct)). This proves (15). A fixed finite union of cosets and fixed finite multiplicities only changes the constants. The estimate is uniform through zero-mode crossings and irrational sequences.

[proved, disappearance from the pressure] At fixed beta>0, (15) gives |log Theta_(alpha_L)(2pi beta/L)|=O(1+log L). Its contribution -log Theta/(beta L) therefore vanishes as L->infinity, uniformly over all offsets. The finite-size Casimir contribution also vanishes per length. Fixing a bounded number of charges to their rounded equilibrium values removes some sums and leaves the same conclusion, provided the corresponding oscillator modes are still allowed. Fixing an extensive different density changes the reference V and e_0 and is a different problem.

[proved, oscillator calculation] For one chirality,

    log Z_osc = -sum_(r>=1) log(1-exp[-2pi beta u r/L]).

Let g(x)=-log(1-exp(-x)). It is positive, decreasing and integrable, with g(x)=-log x+O(x) near zero and exponential decay at infinity. For mesh a, monotonicity bounds a sum_(r>=1)g(ar) between integral_a^infinity g and integral_0^infinity g. The difference is at most integral_0^a g, which tends to zero. Consequently

    lim_(L->infinity) log Z_osc/L
          = (1/(2pi beta u)) integral_0^infinity g(x) dx
          = pi/(12 beta u).                               (16)

For the last equality, expand g(x)=sum_(k>=1)e^(-kx)/k and integrate its nonnegative terms to get sum k^(-2)=pi^2/6. C46 gives a direct Fourier-series evaluation of this sum. Combining (15)-(16) for all chiralities proves (2), with L->infinity before beta->infinity. A fixed positive momentum cutoff changes the result by an exponentially small low-T tail, since all u_i>0 are fixed.

[proved, fixed-aspect limit has a different use] If beta/L is held fixed, Theta_alpha stays of order one rather than tending to a negligible universal constant in log Z. It depends on the moving sector. Thus a torus partition function at fixed aspect ratio need not have a single limit along all L. The subsequential statement of Section 4 is the appropriate one. The pressure proof instead takes the thermodynamic limit at each fixed beta and then T->0. These limits must not be exchanged merely because oscillator gaps scale as 1/L.

[proved conditional on ensemble control, physical-coordinate conversion] In cell coordinates let the mean physical cell length be ell=b+nu. Then L=ell M, k_phys=k_cell/ell and u_phys=ell u_cell. Dividing the free energy per cell by ell sends its coefficient to

    A_phys=A_cell/ell=(pi/12) sum_i 1/u_(i,phys).             (17)

This is covariance under the length conversion, not an additional boson or a factor of ell in c. For the common-bond reference the pressure transfer is already C46's theorem Q. For the dressed microscopic theory it requires that the fixed physical ensemble is covered by G/P or by a separate comparison.

[proved, density optimization does not add a leading coefficient under uniform control] Suppose on a compact set of densities

    g_T(rho)=g_0(rho)-A(rho)T^2+r_T(rho),
    sup_rho |r_T(rho)|=o(T^2),

with continuous g_0,A and a unique minimizing rho_0 for g_0, and all relevant thermal minimizers in that set. Then min g_T=g_0(rho_0)-A(rho_0)T^2+o(T^2). Indeed compactness and uniform convergence force rho_T->rho_0. Minimality gives

    0 <= g_0(rho_T)-g_0(rho_0)
       <= T^2[A(rho_T)-A(rho_0)]+o(T^2)=o(T^2).

Substitution finishes the proof. This is C46's minimization lemma, equally valid for several densities. If there are several ground-energy minimizers, the leading coefficient instead selects the largest A on that minimizing set, under the same compactness and uniformity assumptions. No all-ring ground-state uniqueness is needed for the unique-density version.

[proved, summing charge sectors] At fixed L, a fixed number of bounded integer charges has polynomially many extensive sectors. For positive sector partition functions, max Z_sector <= sum Z_sector <= (# sectors) max Z_sector. The corresponding free-energy error is at most T log(# sectors)/L, which vanishes at fixed T. Together with uniform sector thermodynamics and a dense density mesh this gives the minimization just used. It does not establish the required uniform sector thermodynamics for the interacting wrapper; that remains an application hypothesis.

## 7. A quantitative scalar RG criterion for the lattice thermal coefficient

[proved, purpose of this section] A correlation theorem often normalizes the generating functional by its value at zero source. That removes the scalar information needed for f(T). The following theorem states a sufficient scalar estimate and proves the thermal conclusion by summing it. It does not assert that correlation estimates automatically give this scalar estimate.

[assumed, P, scalar RG contract] Work in the physical thermodynamic ensemble at fixed microscopic parameters. Let q_h=E_* 2^(-h), h=0,1,..., where E_*>0 is below the last linearization and detuning crossover. Let f_G be the Gaussian pressure with the final, dressed V. Suppose the exact, ground-subtracted difference has an absolutely convergent representation

    [f(T)-e_0]-[f_G(T)-e_G]
       = r_UV(T) + sum_(h>=0) delta_h(T),                   (18)

where delta_h is a difference between the finite-T and zero-T scalar contributions on scale h after the dressed Gaussian comparison has been subtracted. For some p>2 and a bounded nonnegative sequence eta_h->0, suppose

    |r_UV(T)| <= C E_*^2 (T/E_*)^p,
    |delta_h(T)| <= C q_h^2 eta_h min{1,(T/q_h)^p}.          (19)

Constants may depend on N, every fixed velocity ratio and the chosen couplings. They need not be uniform as N->infinity or u_min->0. An additional exponentially small massive thermal term is harmless. In the power-law case require eta_h<=C (q_h/E_*)^omega, omega>0. If densities will be optimized, (18)-(19) must hold uniformly in their compact neighborhood.

[proved, scalar transfer theorem] Under P,

    f(T)=e_0-(pi/12) [sum_i 1/u_i] T^2+o(T^2).                (20)

In the power-law case, for T<=E_*/2, the error is bounded by a constant times

    T^2 (T/E_*)^min(omega,p-2),                  p != 2+omega,
    T^2 (T/E_*)^omega [1+log(E_*/T)],             p = 2+omega. (21)

In these expressions sum_i 1/u_i means the sum over all 2N chiral speeds; the multiplication by T^2 is outside the sum. With paired speeds (20) is the free-energy formula in the question.

[proved, summation with a power] Let H be the integer with q_H<=T<2q_H. For h>=H, (19) gives

    sum_(h>=H)|delta_h|
       <= C E_*^(-omega) sum_(h>=H) q_h^(2+omega)
       <= C' E_*^(-omega) T^(2+omega).                      (22)

For h<H it gives

    sum_(h<H)|delta_h|
       <= C T^p E_*^(-omega) sum_(h<H) q_h^(2+omega-p).       (23)

If 2+omega-p>0, the geometric sum is controlled by E_*^(2+omega-p), giving O(E_*^(2-p)T^p). If 2+omega-p=0, there are H terms and the bound is O(E_*^(-omega)T^(2+omega)log(E_*/T)). If 2+omega-p<0, the last term controls the sum and gives O(E_*^(-omega)T^(2+omega)). Together with r_UV these are (21), hence o(T^2). Adding (2) proves (20).

[proved, summation with no prescribed rate] Divide (19) by T^2 and write h=H+j. Apart from constants independent of H, the lower-energy terms j>=0 are bounded by 2^(-2j) eta_(H+j), and the higher-energy terms j<0 by 2^((p-2)j) eta_(H+j). These are summable weights because p>2. For any fixed finite range of j, eta_(H+j)->0. The remaining tails are uniformly small since eta is bounded and the geometric weights are summable. This proves the sum divided by T^2 tends to zero. The UV term also tends to zero after division by T^2. Thus even logarithmic approach to the Gaussian point is enough if the scalar thermal comparison obeys (19).

[proved, why the hypotheses cannot be shortened to running couplings] The factor q_h^2 eta_h alone gives an O(E_*^2) bound on the sum, not o(T^2). If p=2, one fixed high-energy scale can contribute delta_0(T)=kappa T^2 for T<E_*, consistently with a bounded eta_0 and eta_h=0 for all h>0. The running correction vanishes identically in the far infrared but the proposed coefficient has shifted. This is an obstruction to that abstract estimate implying (20), not a claim that arbitrary scalar terms arise in a local Hamiltonian. The power p>2, or a proven cancellation replacing it, is load-bearing. Two time moments in an effective interaction norm by themselves do not state (19) with p>2.

[assumed, what an extension of C50 must prove to use P] Keep the zero-source normalization through every RG step. Compare the scalar part to the Gaussian model with the limiting velocities, subtracting the common ground-energy renormalization. Then prove a temperature-sensitive bound of the form (19), or another summable bound whose sum is o(T^2). It must include the final thermal scale, any temperature dependence of matched couplings, and the physical cyclic/length ensemble. Estimates of amputated vertices or beta functions with the scalar term dropped do not provide (18). This note proves the summation and coefficient identification, not the missing determinant, Ward-identity or thermal-winding estimates that would establish P for the microscopic action.

[proved, independent logical roles] G proves the Virasoro charge of the complete limit. P proves a lattice thermal coefficient using the matched Gaussian reference. P alone need not identify a local stress tensor; G alone, phrased as an infinite-line scaling limit, need not imply uniform finite-temperature or finite-size convergence. F is a separate way to identify finite-size towers. Any claim that all three follow from a particular constructive theorem must check what that theorem actually controls.

## 8. What weaker correlation information does and does not identify

[proved, a correlation formulation sufficient for the current subtheory] A more direct way for C50 to supply the current part of G is to prove convergence of all smeared current vacuum correlation functions, including their normalizations and real-time commutators, to the Gaussian ones after the canonical transformation (4). Equality of these moments identifies the current Hilbert space: map each polynomial of smeared currents applied to the limiting cyclic vector to the same polynomial in Gaussian currents. The moment identities preserve every inner product, hence all null vectors, so this map extends to a unitary between the closures of these polynomial spaces. This is a direct Gram-matrix argument. Two-point exponents alone do not provide the required moment identities. To identify the full compact theory, the analogous identification must cover the charged vertex sectors and exhaust the infrared space.

[proved, local stress in the identified Gaussian subtheory] In each normalized chiral coordinate, <J(z)J(w)>=(z-w)^(-2). Define T(z)=(1/2):J(z)^2: in the limiting Gaussian theory. Wick contraction gives

    T(z)T(w) = 1/[2(z-w)^4] + 2T(w)/(z-w)^2
                            + partial T(w)/(z-w) + regular.

There are two double contractions, each with factor 1/4, and four single contractions, each with factor 1/4. Expanding the remaining normal-ordered J(z)J(w) gives 2T(w)+(z-w) partial T(w). Thus the local stress has c=1, consistent with the mode proof. This constructs a stress tensor within the limit. It does not interchange the lattice scaling limit with a coincident-point composite limit or assert convergence of a particular microscopic stress operator. Completeness is still needed to identify this stress as accounting for all modes.

[proved, normalized sources lose the scalar] If Z_beta[J] is replaced by exp(L a(beta)) Z_beta[J], the normalized generating functional Z_beta[J]/Z_beta[0] is unchanged. The free energy shifts by -a(beta)/beta. Taking a(beta)=kappa/beta shifts the T^2 coefficient. This repeats C46's information-loss observation; it is not a temperature-dependent modification proposed for the chain. It explains why a proof based only on normalized source correlations cannot reconstruct the zero-source scalar without an additional argument.

[proved, a temperature-independent spectator example] Tensor a given Gaussian N-boson model with an independent compact boson, and couple all selected sources only to the first factor. Then Z[J]=Z_N[J] Z_1[0], so every normalized correlation of those sources, at every temperature, agrees with the N-boson model. Its inferred running interactions also agree in that sector. Sections 3 and 6 give central charge N+1 and an extra contribution pi T^2/(6v_s) to the magnitude of the thermal free-energy term for a symmetric spectator of speed v_s. This is an actual Hamiltonian-level example showing why selected correlations do not exclude a neutral sector. It does not satisfy the completeness clause of G and is not proposed as a competing microscopic realization of C44.

[proved, weaker current bound under a CFT premise] Suppose a unitary chiral CFT has N independent currents of positive level, normalized as in Section 3, and its full stress generators obey [L_m^full,J_l^a]=-l J_(m+l)^a. Let L_m^J be their summed quadratic stress. Then K_m=L_m^full-L_m^J commutes with every current, hence with L_l^J. Subtracting the algebras shows that K is Virasoro with central charge c_full-N. On the vacuum,

    ||K_(-2)|0>||^2 = (c_full-N)/2 >= 0.

The equality follows from [K_2,K_(-2)]=4K_0+(c_full-N)/2 and annihilation of the vacuum by the appropriate modes. Thus c_full>=N. This is C46's conditional Sugawara bound, not a replacement for identifying the complete infrared theory. Applied to both chiralities it gives c_cd>=N. It still requires dimension-one chiral currents and their stress transformation; the microscopic single payload U(1) is not itself N such currents.

[proved, thermal degeneracy of central-charge interpretations] One symmetric boson of speed 1 has A=pi/6. Two symmetric bosons of speed 2 have the same A but c_cd=2 rather than 1. Therefore no unweighted c can be read from a single thermal coefficient without velocities and a mode/stress identification. With several velocities there is no justified replacement of them all by the smallest or largest speed in an equality for c.

[proved, no interval-entropy claim] None of G stated only as a local scaling limit, F stated only as a spectral approximation, or P stated only as a pressure estimate includes trace-norm/entropy control of interval density matrices. The interval coefficient S(ell)=(N/3)log ell+O(1) for the microscopic chain is not proved in this note. Entropy convergence would require a further argument. The definition of c delivered by G is the Virasoro central term of the full Gaussian scaling theory, with the optional operational thermal signature (20) under P.

## 9. Application, limit order and remaining work

[proved conditional on the stated inputs] Fix n and freeze the nonzero microscopic couplings selected in C45. Assume C50 establishes G for the complete physical infrared theory with N=n+1 modes, including the length mode. Then Sections 2 and 3 prove c_cd=n+1. No common speed, exactly rational density, convergence of the same ground sector along all rings, or uniqueness of the finite-ring ground state is needed. If C50 also supplies P, Sections 6 and 7 prove the microscopic formula (20). If it supplies F and the physical charge-lattice identification, Sections 4 and 5 prove the finite-size formula (11) with its moving sectors and multiplicities.

[proved, order of limits] In the thermal route, first L->infinity at fixed T for this fixed n and fixed bond, then T->0 below every fixed crossover, including the smallest speed/detuning scales. In the spectral route, first fix the number of rescaled levels and let L->infinity for that fixed bond, allowing the offsets alpha_L to move or taking subsequences. Only after each individual infrared limit may n increase. Constants may deteriorate without bound with n and speed ratios. No inference here requires uniform access to all channels at one common laboratory temperature or system size.

[assumed, remaining work for the actual construction] The completed microscopic statement is still conditional. C50's inspected partial proof explicitly leaves its all-scale comparison and source bounds unproved. Its requested eventual correlation and beta-function conclusions must be upgraded to either complete local Gaussian current/field identification for (1), or temperature-sensitive scalar control for (20), or full finite-size matching for (11). The actual compactification lattice and physical gluing need identification for a claim of exact lattice tower multiplicities. C48/C49's decoupling and excluded-sector control are also assumed in the intended route. C46's pressure comparison for the block-diagonal reference does not supply the corresponding full finite-coupling comparison.

[proved, concrete deliverable and boundary] This closes the continuum identification once "a rigorously controlled multichannel Luttinger liquid" means the complete Gaussian theory G. It gives an explicit sufficient analytic criterion P for the additional lattice thermal statement, with a rate and a proof. It does not prove P, F or microscopic completeness from power-law correlations and a vanishing beta function alone. Accordingly it is not an unconditional c=n+1 theorem for the 4464-state bond and does not change the project's conditional counterexample status.

## 10. Reproducible checks and write scope

[numerical] Run `python3 -B checks.py`. It uses only numpy and the Python standard library and writes checks.json in this directory. The largest oscillator matrix has dimension 139. It checks (4)-(5) for mixed positive matrices at N=1,2,3, with maximum entrywise residual 3.56e-14. A separate diagonal fixture allows coincident speeds and a speed ratio 9e7. This is an algebra fixture, not a conditioning certificate for C42's microscopic parameters.

[numerical] Truncated Fock operators at charges 0, 0.37 and sqrt(2) satisfy all Virasoro commutators with |m|,|l|<=3 on levels through three, with an adequate cutoff buffer and maximum error 6.22e-15. The directly enumerated partitions agree with the integer recurrence. These controls test the normalization and nonzero-charge subtraction in (9)-(10), not a numerical inference of c from the lattice.

[numerical] An abstract two-coordinate moving-charge fixture over L=128,...,512 has uncorrected L(E_0-e_0 L) ranging from -0.6210973883 to 2.3296086372. Subtracting its quadratic zero-mode term leaves the constant -0.6283185307 with zero reported floating-point residual. At beta=2.3 and L=256,1024,4096 its zero-mode pressure magnitudes decrease as 0.00690604, 0.00231512, 0.000725933. The bounds (15) pass. The fixture is not the wrapper's unidentified compactification lattice.

[numerical] Synthetic shell corrections saturating (19) pass the rates (21) for omega=0.4 and p=2.25,2.4,4. The equality case has the predicted logarithm. For eta_h=1/(h+1), p=4 and T/E_*=2^(-k), k=8,32,128,512, the error divided by T^2 decreases through 0.187422,0.0505465,0.0129206,0.00324887. The p=2 single-shell control stays constant, as the analytic obstruction predicts. No microscopic scalar RG estimate is tested by these fixtures.

[proved, limitations and scope] No interval entropy, microscopic finite-size spectrum, or microscopic thermal coefficient was computed. The numerical controls are supplementary to the proofs of the conditional implications. The first check run emitted an overflow warning from an unused numpy branch, which was fixed by selecting the exponent before exponentiation; the final run has no warnings and all assertions pass. All writes are inside work_C51. No network, subagents, cluster computation or project-level edits were used. STATE.md here is the handoff because the task forbids writing the parent project STATE.md.
