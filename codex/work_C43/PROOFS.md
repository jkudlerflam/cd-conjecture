# C43 derivations

[proved, scope] The finite statements below concern the literal C39 core bond and its vacancy lift, or are explicitly conditional statements for an already selected squeezed period. They do not identify the initialized core with an autonomous selector. The continuum dictionary and thermodynamic matching are labeled separately. All calculations use local files and numpy, scipy, or the Python standard library.

## 1. What is being completed

[proved, inherited finite input] C39 PROOFS Sections 1 and 7 specify a 46-letter initialized core with an open word A, n data sites, B. The positive boundary operators enforce A on the left and B on the right. Its code has dimension 2^n. The core alphabet alone does not select n or initialize such cells on a periodic chain. C39 explicitly counts the complete wrapper separately as 96*46+48=4464. This distinction is the obstruction to the proposed direct substitution.

[proved, definition] Write h_core for the repeated 2116 by 2116 bond in C39/literal.py, without its two external boundary operators. Let V be the sum of its two bus port operators from C39 PROOFS (5), supported on the ordered pair B,A. Add the nonnegative conditional data readout R with coefficient r lambda^2, r>=0. The concrete repeated occupied bond audited here is

    h = h_core + lambda V + (r lambda^2/2)(R tensor I + I tensor R).

[proved, definition] On C|0> direct-sum C^46, with occupied projector n_occ, define

    h_47 = -t sum_a (|0,a><a,0| + |a,0><0,a|)
           + J iota_2 h iota_2^dagger
           - (mu/2)(n_occ tensor I + I tensor n_occ).

[proved] Summing this bond on a ring gives a Hermitian, one-site translation invariant, nearest-neighbor Hamiltonian of local dimension 47. The complete occupied bond, including its half-site terms, is gated on occupied pairs. A physical single-site version of those terms would be a different lift with different relative first-order coefficients. There can be an extensive number of vacancies. The phrase one vacancy refers to one added local letter.

[proved] Assign charge zero to the vacancy and use C39's integer charge on the 46 occupied letters. A swap moves the charge with its letter; every matrix element of h conserves the sum of its two charges. Therefore occupied number N and the original charge Q_pay commute with H exactly. The parity exp(i pi Q_pay) is exact as well. It is unrelated to the hard-core seam sign (-1)^(N-1). The numbers of A and B markers are also conserved by these rules. Separate logical flavor numbers need not be conserved.

## 2. Exact squeezing and its limited factorization

[proved] Fix 0<N<L. Every physical basis vector is specified uniquely by ordered occupied positions x_1<...<x_N and an ordered internal word a=(a_1,...,a_N). Hence

    H_(L,N) = C^{binomial(L,N)} tensor (C^46)^{tensor N}

as a vector space with its orthonormal basis. There is no extra factor N for choosing a reference particle. A vacancy hop away from the cut changes x only. A hop from L-1 to 0 sends the internal ket to

    T|a_1,...,a_N> = |a_N,a_1,...,a_(N-1)>.

[proved] In this basis the exact Hamiltonian is

    H = H_hop[T] - mu N + J sum_(j=1)^N A_j(x) h_(j,j+1),
    A_j(x) = 1 if x_(j+1)=x_j+1,
    x_(N+1)=x_1+L.

[proved] The last h acts on the ordered internal seam N,1. The hopping seam carries T. This follows by applying each microscopic bond to a basis word. It proves the identity for arbitrary h, including illegal words and interacting h. The indicators A_j fluctuate. Consequently this identity does not say that the finite-J Hamiltonian is a sum of two independent chains.

[proved] At J=0, diagonalize T with eigenvalue e^(iQ), Q=2 pi s/N. The charge problem is hard-core hopping with seam phase e^(iQ). On mapping the ordered hard-core amplitudes to antisymmetric fermion amplitudes, moving one particle past the other N-1 contributes (-1)^(N-1). With orbital coefficient e^(-ikx), the allowed momenta and energies are

    k_I = [2 pi I + Q + pi(N-1)]/L,
    E_charge = -2t sum_(I occupied) cos k_I - mu N.

[proved] All distinct orbital subsets of size N must be included, and each Q sector is repeated with its actual multiplicity in the internal cyclic representation. This is a unitary diagonalization of the whole fixed-N space, not an extra choice of spin structure. Reversing the orientation of T reverses Q and all orbital momentum conventions together.

[proved] An internal cyclic projector P_Q obeys

    P_Q h_j P_Q = (1/N) P_Q H_s P_Q,
    H_s = sum_j h_j.

[proved] To derive this, use h_j=T^(j-j0) h_j0 T^(-(j-j0)), with an orientation-dependent sign in the exponent that cancels between equal Q sectors. All diagonal-sector compressions are equal, and their sum is P_Q H_s P_Q. For a chosen normalized charge state f_Q the diagonal block of the first derivative at J=0 is therefore

    P_Q <f_Q|sum_j A_j h_j|f_Q> P_Q
        = alpha_(L,N,Q) P_Q H_s P_Q,
    alpha_(L,N,Q) = <f_Q|sum_j A_j|f_Q>/N.

[proved] Individual rooted contact expectations need not be equal. This argument uses only their sum and the internal cyclic projector. Cross-Q blocks contain <f_Q|A_j|f_Q'> P_Q h_j P_Q' and are not generally zero. Degenerate charge seas require degenerate perturbation theory. Neither these blocks nor low-energy charge excitations may be discarded in a uniform thermodynamic construction.

[proved] For any consecutively filled Slater sea on the finite twisted charge grid,

    <n_x> = N/L,
    |<c_x^dagger c_(x+1)>| = sin(pi N/L)/[L sin(pi/L)],
    alpha_(L,N) = N/L - sin^2(pi N/L)/[N L sin^2(pi/L)].

[proved] The second identity is the finite geometric sum over occupied orbitals. Wick's two-by-two determinant gives <n_x n_(x+1)>=(N/L)^2-|<c_x^dagger c_(x+1)>|^2. Sum over x and divide by N to obtain the third. At a sea degeneracy this formula applies to each consecutive Slater representative, not automatically to arbitrary superpositions.

[proved, thermodynamic free-charge limit] At N/L -> rho in (0,1),

    alpha(rho) = rho - sin^2(pi rho)/(pi^2 rho) > 0,
    v_c = 2t sin(pi rho).

[proved] Positivity follows from sin(pi rho)<pi rho. At rho=1, alpha=1; as rho tends to zero, alpha=pi^2 rho^3/3+O(rho^5). At half occupancy alpha=1/2-2/pi^2. With exactly one vacancy, the finite formula gives alpha=(N-1)/N, consistent with the deterministic number N-1 of occupied contacts.

[assumed, infrared identification] The interior free charge sea supplies a c=1 boson with K_c=1. Applying J alpha H_s to a critical payload at fixed J>0 requires local low-frequency matching with the soft charge mode retained. It does not follow from finite-ring nondegenerate perturbation theory, whose charge denominators vanish with L.

## 3. A finite obstruction for the actual 47-state bond

[proved] Take a=4, |g|<1, 0<q<1, w,wv<=1, D>=0 and positive grammar penalty K, as allowed in C39. Its positive local decomposition proves h_core>=0. The pair block B,A is invariant under h_core. It has grammar energy K, nonnegative excited diagonal terms, and no vacuum kinetic edges or data conversions. Each F in a port has norm one, so

    ||V|| <= 3(1+sqrt(eta)).

[proved] If 3|lambda|(1+sqrt(eta))<K, the B,A block of h has strictly positive lower bound K-3|lambda|(1+sqrt(eta)). Every other pair block keeps the original positive decomposition, and the readout is positive. Thus h>=0. The condition is compatible with C39's final small-lambda choice for every n. It is not restricted to the diagnostic ED parameters.

[proved] The periodic squeezed Hamiltonian H_s=sum_j h_j, for N>=2, has kernel

    ker H_s = span{|L_(a1)...L_(aN)> : a_j=0,1}
              direct-sum
              span{|R_(a1)...R_(aN)> : a_j=0,1}.

[proved, kernel completeness] A zero vector is annihilated by every positive local bond. Every A has an incoming penalized pair, and every B has an outgoing penalized pair. The B,A exception remains strictly positive by the preceding bound, so no A or B can occur. On the remaining letters the positive onsite excited mass is at least 1-|g| per excited token in the ring, so excited tokens are absent. The remaining symbols are L_a, R_a, d_(v,a). Their grammar arrows are L to L, R to R, and L to R respectively. A cyclic legal word cannot contain a one-way L-to-R transition without a reverse transition, which is unavailable. It is therefore entirely L or entirely R. All such words are annihilated by every bond, port and readout. The diagonal grammar penalties preclude canceling this obstruction by superposition.

[proved] The charge hopping ground energy in any Q sector is bounded below by its Q=0 value. One direct derivation fills the N lowest cosine orbitals. Their consecutive-sea sum gives

    E_min(Q) = -A_(L,N) cos(dist(Q,2 pi Z)/L) - mu N,
    A_(L,N) = 2t sin(pi N/L)/sin(pi/L) > 0.

[proved] The hard-core parity shift cancels the parity of the consecutive-sea center, leaving center j+Q/(2pi), j integer. Minimizing its distance from zero yields this formula. The minimum is attained only at Q=0 modulo 2pi. This also follows from the strict flux version of the finite hopping diamagnetic bound. The Q=0 charge vector is unique and has strictly positive amplitude on every ordered position configuration by connectivity and Perron-Frobenius.

[proved] Since J sum A_j h_j>=0, the full fixed-N ground energy is at least -A_(L,N)-mu N. It attains this bound by tensoring the Q=0 charge ground vector with any cyclic-invariant all-L or all-R word superposition. Conversely, saturation requires Q=0 and the unique charge vector. Because that vector is nonzero on every configuration and each squeezed neighbor pair can be a physical contact, positivity forces h_j chi=0 for every j. The preceding kernel computation is therefore exhaustive.

[proved] For 2<=N<L the exact ground multiplicity is

    g_N = (2/N) sum_(r=0)^(N-1) 2^gcd(N,r).

[proved] In each binary L or R word space, the invariant dimension is Tr[(1/N)sum_r T^r]. A word fixed by rotation r is constant on its gcd(N,r) cycles, giving trace 2^gcd(N,r). This proves the formula, including the factor two for L and R. It is at least 2^(N+1)/N. Neither doping nor the readout lifts this extensive degeneracy.

[proved] At N=1 no occupied pair acts and the species multiplicity is 46. At N=0 there is a single vacuum. At N=L hopping vanishes and the all-L/all-R kernel has dimension 2^(L+1). For chemical potentials giving interior free-gas density, the grand-canonical minimizing N is interior on all sufficiently large rings, so the exponential degeneracy persists even if the winning N is unique. The weak ports cannot remedy it because their support contains endpoints absent from this ground space.

[proved, consequence] The literal dimension-47 completion has no selected n-mode payload, no specified spontaneous period b, and no unique large-ring ground state. It cannot establish c=n+1. This is an obstruction to this bond, not a no-go theorem for an enlarged autonomous selector. Importing C39's initialization by fixing endpoints or by restricting to legal cells would change the task's microscopic periodic chain.

## 4. Gluing a selected period, conditional on supplying a selector

[assumed] For this section only, replace h by a complete autonomous occupied bond whose low-energy squeezed configurations have period b, N=bM, a positive cost for defects of that period, and n critical logical bands per cell. Such a bond is not supplied by C39's 46 letters. Its one-vacancy lift would have dimension d_occ+1. The equations below describe what a successful additional selector would need to realize.

[proved, finite kinematics under that hypothesis] On the union of b offset sectors, cyclic translation T by one occupied site satisfies T^b=T_cell. If T_cell has eigenvalue e^(iP), the eigenvalues of T are

    Q_r = (P+2pi r)/b,    r=0,...,b-1.

[proved] This is diagonalization of a cyclic shift through b offsets whose last step applies T_cell. Its characteristic polynomial is z^b-e^(iP). Each root enters the charge seam from Section 2. The b offsets are physical states reorganized into twist sectors, not a factor b of identical vacua to be deleted by convention. Generally their charge energies differ by order 1/L. Their finite-J mixing is part of the full glued theory.

[proved, free reference] The charge sea current label is

    W_c=j+Q_r/(2pi)=(k+P/(2pi))/b,    k=bj+r in Z.

[proved] Expanding the exact sea sum gives its current energy (2pi v_c/L) W_c^2+O(L^-3) for bounded W_c. The number sector is N=bM. Thus low-energy number changes are multiples of b while currents occur in steps 1/b, with an additional payload momentum offset. These two restrictions are the compactification data. Treating the charge boson as an independent ordinary spinless fermion with only integer currents loses the period branches.

[proved, common logical grading] The logical cell tensor product has one total parity per cell. With F total logical fermion number, its ordinary cyclic translation is represented by fermions with seam sign (-1)^(F-1), common to all flavors. To see this, moving F_0 fermions in a translated cell past the other F-F_0 gives (-1)^(F_0(F-F_0)); the seam gives (-1)^(F_0(F-1)). Their product is one since F_0(F_0+1) is even. This is C27's finite algebra. There is no additional graded cell permutation and no constraint on whether n itself is even or odd.

[proved, free logical reference under the selected-cell hypothesis] More explicitly, in M cells choose distinct orbital subsets for each flavor on the common grid p_I=2pi(I+beta)/M, beta=(F-1)/2 modulo one, with F equal to the total number of selected orbitals. The ordinary cell translation momentum is P=sum_(alpha,I occupied) p_I modulo 2pi. It is an allowed M-th-root phase because beta F is an integer modulo one. For displaced bands the energy is the sum of the corresponding dispersions evaluated on this grid; the displacement does not alter the seam sign. Combine this P with Q_r above and the charge grid of Section 2, with N=bM. This specifies the separated-reference ring gluing without choosing independent flavor spin structures. At finite J it is a reference tower construction, not an exact additive spectrum of the microscopic adjacency-gated chain.

[assumed, standard hydrodynamic dictionary] Let y count occupied sites and u count logical cells. Normalize the charge field by

    y(x)=rho x-phi_c(x)/pi,    u(x)=y(x)/b.

[proved conditional on that dictionary] The first marker order harmonic is exp(2pi i u)=exp(2pi i rho x/b) exp(-2i phi_c/b). At K_c=1 its dimension is 1/b^2 and its two-point function decays with exponent 2/b^2. Physical positional long-range order is then lost, while squeezed order remains encoded in the twists. Its wavevector is 2pi rho/b; the average physical spacing is b/rho, not an exact broken-translation period at generic rho.

[proved conditional on a stable defect gap and smooth bulk energy] If the preferred occupied number is rho L, choosing the nearest multiple of b changes N by at most b/2. If e(rho) is the energy per physical site and e'' is bounded by C in the local interval, Taylor expansion at the chemical-potential minimum costs at most C b^2/(8L)+O(L^-2). A competing wrong-period defect has an L-independent positive cost by hypothesis. For large L, defect-free sectors therefore exist at lower cost, regardless of L modulo b. This argument uses grand-canonical minimization. A fixed incompatible N is a different problem. Interacting Gaussian shell energies can add bounded O(1/L) terms without changing this comparison.

[proved, limitation] The gluing and rounding comparison do not prove uniqueness within the winning sector or absence of crossings between sectors. Equal twist minima and conserved-charge crossings are possible. C39 has complex port phases, so a positive-amplitude Perron-Frobenius argument for its complete payload is unavailable. Genericity of parameters alone is not a proof of a unique ground state on every sufficiently large ring. Even unique finite-ring ground states can have size-dependent charge and current offsets, rather than one fixed untwisted vacuum tower along every integer L.

## 5. Complete neutral vertices and C39's margins

[assumed, continuum operator dictionary] In the selected-period theory, label a payload vertex by signed integer chiral counts z_(alpha,s), with annihilations positive. Exact original U(1) requires sum z=0, stronger than the even-sum parity constraint. Its cell momentum and free dimension are

    q_pay = sum_(alpha,s) p_(alpha,s) z_(alpha,s),
    Delta_pay = (1/2) sum_(alpha,s) z_(alpha,s)^2.

[proved within the dictionary] An even finite-interval tensor-product operator is an even finite CAR polynomial: Jordan-Wigner strings outside the interval cancel. Conversely every such polynomial is supported on that interval. Repeated chiral counts are represented by point-split fermions; their primary dimensions already include the Pauli cost. Exact occupied-number conservation forbids charge-dual vertices. Derivatives add nonnegative integer dimensions. Odd payload twist operators cannot appear in this U(1)-preserving Hamiltonian.

[proved conditional on the dictionary] A cell-local operator has all integer cell-comb harmonics m. Substituting u(x) gives

    exp[i(q_pay+2pi m)u(x)] V_z(x)
      = exp[2pi i rho w x] exp[-2i w phi_c(x)] V_z(x),
    w=(m+q_pay/(2pi))/b,
    k_phys=2pi rho w mod 2pi,
    Delta_ref=w^2+(1/2) sum z^2.

[proved] This follows from sum_j delta(u-j)=sum_m exp(2pi i m u). Replacing q_pay by q_pay+2pi changes the harmonic label, not the set of physical operators. The fractional w values are essential. Formulae for the two-length wrapper's independent length harmonics cannot be inserted here without this change of normalization.

[proved conditional on the dictionary] Physical translation invariance requires rho w=ell, ell integer. At resonance,

    Delta_ref,res = ell^2/rho^2 + (1/2) sum z^2.

[proved conditional on the dictionary] Choose rho=1/2. Every nonzero ell then contributes at least four. At ell=0, q_pay=0 modulo 2pi. Neutral nonzero integer z has even squared norm. The only shells below six have squared norm two or four. Norm two is a creation-annihilation bilinear, excluded by distinct Fermi points. Norm four with sum z=0 consists of two creations and two annihilations on distinct branches; its resonances are equal unordered pair sums. Coincident-index representatives either have larger squared norm or are neutral current descendants. C39 certifies distinct unordered pair sums including repeated labels, so these comparisons are all excluded except forward scattering, whose vertex charge is zero.

[proved conditional on the dictionary] All other nonzero neutral payload vertices have Delta_pay>=3. Derivatives cannot lower it. In the zero-vertex sector, the identity, linear currents, and quadratic currents renormalize energy, density, current and Gaussian stiffness. They do not supply a non-Gaussian cosine. Thus every exactly resonant non-Gaussian vertex at half occupancy has reference dimension at least three.

[proved conditional on the dictionary] Let delta_int=min(delta2,delta4)>0 be C39's fixed-n margin, reduced if necessary so delta_int<=2pi. For every non-Gaussian primary with Delta_ref<3, |w|<sqrt(3), hence |k_phys|<pi sqrt(3). Its distance from any nonzero physical reciprocal vector is at least pi(2-sqrt(3)). Its distance from zero is at least delta_int/(2b) for a nonforward bilinear or quartic, or pi/b for a pure nonconstant marker harmonic. Therefore

    delta_phys >= min(delta_int/(2b), pi(2-sqrt(3))) > 0.

[proved] The test is finite despite the arbitrary integer m. The dimension cutoff bounds both sum z^2 and |m+q_pay/(2pi)|. Higher shells retain an irrelevance bound if the full scaling metric G in these charge/payload field coordinates obeys G>=gamma I with gamma>2/3. For example G>=4I/5 gives Delta>=12/5 for every reference-dimension-three vertex. This metric condition includes density-current mixing and is separate from positivity of the Gaussian energy matrix. It is not a bound on a bare hopping coefficient.

[proved conditional on the same dictionary] A small density interval can replace half occupancy. For the dimension-three separation choose rho_max<1/sqrt(3), giving the nonzero reciprocal margin 2pi(1-rho_max sqrt(3)); the zero reciprocal margin is at least rho_min delta_int/b. This is a finite open exclusion, not a requirement of irrationality or of infinitely many Diophantine bounds. Half occupancy is especially useful for the charge momentum threshold below.

[proved, free-charge threshold] At half occupancy write each orbital relative to its closest Fermi point as k=s*pi/2+d, |d|<=pi/2. Since sin|d|>=2|d|/pi, its particle or hole cost is at least 4t|d|/pi. For any number-conserving collection of particle-hole changes, the Fermi-point part of the momentum is an integer multiple of pi. The triangle inequality then gives

    E_exc >= (4t/pi) dist(P_exc,pi Z).

[proved, thermodynamic application at fixed b] A change of internal twist by (q_pay+2pi m)/b changes the continuously followed charge sea momentum by (q_pay+2pi m)/(2b) modulo pi, although its sea energy changes only by O(1/L). If the full perturbation preserves physical momentum, a nonforward process must restore this mismatch with charge excitations. Their free threshold is at least 2t delta_int/(pi b). This includes arbitrary finite-energy collections of particle-hole pairs. Rounding N-L/2 by O(b) changes the limiting threshold by zero; bounded-energy current boosts grow at most as sqrt(L), so their rounding correction also vanishes. The threshold concerns momentum restoration, not an artificial gap between charge twist vacua. The collaborator's CHARGE_THRESHOLD.md supplies the same independent derivation.

## 6. Lacunary schedule and the missing matching

[proved, inherited C39 input] For each fixed n, C39 selects w,u,v,eta,mu_field, then D, then g, then lambda, with all parameters independent of L. Its strictly positive delta_int and minimum velocity can be very small. Its certified limiting examples have delta4>5.11e-16, 1.23e-28, 2.36e-41 at n=1,2,3. These numbers are inherited, not recomputed or promoted to finite-D ED margins here.

[proved, first-order scaling] If a complete occupied bond h has those properties, multiplication by J alpha scales its energy gaps and velocities uniformly and leaves its Fermi momenta unchanged. In particular, rescaling a whole core containing port coefficient lambda gives a second-order scale (J alpha lambda)^2/(J alpha)=J alpha lambda^2. Keeping the field inside the same gated bond preserves its relative coefficient. The scale hierarchy is therefore compatible with first-order squeezing.

[assumed, sufficient finite-J matching] After fixing the complete occupied h, assume local low-frequency integration of massive machinery and high-energy charge channels, with the charge boson and all payload bands retained. The normalized correction to this effective action must tend to zero as J/t tends to zero at fixed h. It must preserve the exact charges, the compactification above, and a positive period-defect cost. Corrections must be small relative to the slowest payload stiffness and physical detunings, not merely relative to J or the fastest band. A finite-size extensive Schrieffer-Wolff remainder does not imply this assumption.

[heuristic, ordered hierarchy] First use C39's small g to suppress its g^4 interaction relative to its g^2 hopping, then lambda to suppress higher port orders, including the field leakage. Complete and stabilize the selector at this fixed n. Finally take J/t sufficiently small at these fixed choices, and adjust the occupied chemical potential to half occupancy. The charge cutoff from the previous section is positive at fixed n even when tiny. Smooth charge-density times payload-current couplings are Gaussian marginal terms; charge density times payload energy has dimension three. Residual nonforward terms are cut off by delta_phys if their amplitudes stay weak until that scale. This is the usual weak-coupling mechanism, conditional on the preceding matching input. No uniform-in-n coupling window is asserted.

[heuristic, conditional conclusion] With a supplied autonomous selector, this matching, positive Gaussian energy, and the stated scaling-metric and momentum margins, the infrared theory has n+1 compact bosons. Unequal velocities do not change their total central charge. The literal 47-state bond fails the selector premise, so this implication cannot establish c=n+1 for it. The proof of irrelevance is useful but cannot repair a wrong microscopic ground space.

[proved, reference circumference check] The logical cell circumference is M=rho L/b. A logical velocity v_alpha measured per cell in the occupied reference becomes (b/rho) J alpha(rho) v_alpha in physical lattice units. The charge velocity is 2t sin(pi rho). These are n+1 finite positive speeds at fixed n and fixed nonzero parameters if all the hypotheses hold. Comparing a Casimir coefficient using only one of these speeds would not measure the total central charge. No such fit is attempted here.

## 7. Microscopic finite checks

[numerical, implementation] microscopic.py reproduces C39's 46 physical letters and local rules, including both buses, neutral phases, conversions, grammar terms, weak B,A ports, and conditional readout. It uses no neighboring imports. It then constructs the literal 47-state lift. checks.py constructs physical ring matrices by applying the repeated bond, independently constructs the ordered-coordinate matrix, and rejects any leakage when using an invariant sector.

[numerical, parameters] The diagnostic values are a=4, D=.2, w=.8, v=.65, g=.4, q=.8, eta=.6, K=10, lambda=.03, t=1. Full three-site rings use J=.07, mu=.17 and readout coefficient .03 lambda^2. These values resolve matrix effects and satisfy the positivity bound. They are not the lacunary witness parameters and do not test the extremely small certified margins.

[numerical, full microscopic ED] For each n=1,2 the complete L=3 Hilbert space has dimension 47^3=103823. Connected components of the nonzero off-diagonal matrix entries are exact invariant subspaces, and every component is diagonalized, preserving every eigenvalue and multiplicity. There are 51860 components with largest block 48. No occupied-number, charge, endpoint, token or legality sector is removed. Both spectra have ground energy -2.34, multiplicity six, and next energy -2.17. The occupied bond has numerical minimum eigenvalue zero and its charge commutator vanishes entrywise.

[proved, check of the ED ground] The fixed-N lower energies at L=3 are 0 for N=0, -2-.17 for N=1, -2-.34 for N=2, and -3*.17 for N=3. Thus N=2 wins and the exact multiplicity is 2*(2^2+2)/2=6. This proves the reported finite degeneracy independently of floating-point eigensolvers.

[numerical, larger microscopic sectors] At n=1, (L,N)=(4,3), dimension 648, and at n=2, (L,N)=(5,4), dimension 3200, the retained sector consists of one complete initialized cell, all cyclic offsets, and one mobile vacancy. It includes all 54 or 160 initialized-cell states, respectively, rather than just the 2^n code states. These are exact invariant sectors of the actual microscopic bond. No boundary field is added. They are not the full 47^L spaces or their global ground sectors. For this check only, the readout is set to zero to retain the data degeneracy, and mu=0.

[proved, sector interpretation] In that one-cell squeezed ring the B,A seam contributes K and the missing endpoint diagonals; its port is zero since a one-token cell never has active tokens at both endpoints. The squeezed restriction is consequently the initialized H_core+K. This makes a clean first-order test, but it does not test intercell port transport in a selected many-cell liquid. The full L=3 tests do include the port-active pair blocks, including illegal states.

[numerical] Independent physical and squeezed matrices agree entrywise at J=.02 in both sectors. All J=0 eigenvalues, with multiplicity and all cyclic twists, agree with the free-charge formulas within 8.9e-16 and 1.7e-15. The Q=0 sea projection agrees with alpha_(L,N) H_s within 1.8e-15 and 3.6e-15. The interaction does not commute with internal translation alone, its Frobenius commutator norms are recorded in checks.json. This is a direct finite check against claiming exact factorization at J>0.

[numerical, finite-J low spectra] Sparse ED at J=.02 and .01 gives maximum eigenvector residual below 3.2e-11. With E_1 the lowest eigenvalue of alpha H_s, the remainders (E_0(J)+2-J E_1)/J^2 are -7.56795, -7.48875 at n=1 and -8.12588, -7.96908 at n=2. The nonzero remainder verifies that the first-order product is not the exact finite-J spectrum. It is consistent with finite-dimensional perturbation theory and is not a thermodynamic extrapolation.

[proved, reproducibility and limits] Run `python3 -B checks.py` in this folder. It writes checks.json and sets single-thread BLAS before importing numpy. No computation requires more than minutes. The checks do not measure a central charge, establish a selected many-cell ground sector, diagonalize the 4464-state wrapper, or prove interacting every-ring uniqueness. Those claims are not inferred from these finite spectra.

## 8. Fallback and final dimension

[proved, accounting] No autonomous selector with fewer letters is constructed here. Merely adding unused letters to the direct sum leaves the obstruction. Changing the grammar, adding a length-energy mechanism, and excluding all other sectors would be new finite construction work. Lifting an already complete 4464-state frozen wrapper would instead cost 4465 states and supply no advantage over its existing mobile version.

[assumed, inherited finite wrapper; proved arithmetic] The requested fallback is C39's mobile wrapper, with onsite space

    C^4 tensor [(C^3 tensor C^2 tensor C^46) direct-sum C^3]
         tensor C^2 tensor C^2,
    d0=4*(3*2*46+3)*2*2=4464.

[assumed, inherited scope] This uses C39/C32/C8's initialized wrapper, occupation filters, global exclusion, compensated endpoints, and identity transport. The two selected physical lengths n+3,n+4 carry the same n-mode code. Their noncommuting source-resolvent approximation and all-n core margins are inherited as stated in C39, not independently re-proved by the present vacancy ED.

[heuristic, final infrared status] Under C39's local massive-sector elimination and weak-coupling Gaussian stability inputs, that fixed-d0 family remains a candidate with c=n+1 for every finite n. The couplings depend on n and are fixed before L tends to infinity. This task does not turn the conditional candidate into an unconditional counterexample, and does not prove unique ground states on every large ring for the fallback. The proposed 47-state reduction fails by the finite ground-space obstruction above. The general species-independent vacancy mechanism remains viable for a separately completed selected-period bond.
