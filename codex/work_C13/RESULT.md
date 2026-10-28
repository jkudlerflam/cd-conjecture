# C13. Mobile cells at finite coupling

Status: conditional Gaussian phase, with unresolved microscopic input. This is not a verified counterexample to the c-d conjecture.

[heuristic] The mobile construction can support n+1 gapless bosons at fixed nonzero coupling. The finite-coupling operator analysis does not reveal an unavoidable mixed instability. It does require more than distinct payload fillings. The length density must avoid additional resonances, pairing must be excluded by a conserved charge or by displaced Fermi-sea centers, the complete stiffness matrix must stay positive, and the corrections must be small compared with the smallest retained energy scale. These are substantive conditions, not consequences of a small coefficient in the microscopic bond alone.

[proved] Three statements in the question need qualification. (i) A unique ground state on each finite ring does not imply convergence, along every integer L, to the same CFT vacuum sector. Irrational charge rounding produces moving zero-mode offsets. (ii) Generic uniqueness for all integer L requires avoiding a countable set of sector crossings. It need not hold on an open parameter neighborhood. (iii) Compactness requires bounded or normalized bonds. The c-theorem alone does not establish local boundedness of lattice IR central charge without control of the crossover scale.

[assumed] This report uses standard bosonization, locality of the low-frequency effective action after integrating out gapped cell excitations, and the usual identification of stable Gaussian modes with compact-boson towers. These continuum inputs are stated explicitly. They are not proved here. The initialized high-channel microscopic core is still missing from the supplied construction, so its particular finite-coupling coefficients cannot be computed.

[numerical] Python scripts and JSON data in this folder provide two accurate estimators for the soluble reference and a locally unitary, interacting deformation. Small physical-graph ED with independent density interactions and conditional hopping is consistent with continuity, but does not provide two accurate central-charge estimates for the generic interacting model. The free physical-site entropy calibration itself fails at those small sizes. This limitation is part of the verdict.

## 1. Coordinates and exact symmetries

[proved] Write M for the number of cells, r=L-bM for the number of long cells, nu=r/M, a=b+nu, and rho=1/a. Payload flavour a will instead be indexed by Greek alpha when confusion is possible. Its filling per cell is f_alpha, with Fermi points

    p_(alpha,R)=theta_alpha+k_alpha,
    p_(alpha,L)=theta_alpha-k_alpha,
    f_alpha=k_alpha/pi.

Filling is defined by the occupied interval; reversing the sign convention for the dispersion reverses the corresponding chiral labels. All momenta in this paragraph are in cell-index units.

[proved, from the displayed microscopic moves] In selector_wrapper.md and helix_equalizer.md, every marker move exchanges one marker with one nonmarker, while the port term leaves selector letters unchanged. Thus M is conserved exactly, including before projection. One-site physical translation is exact. The selector-violation projectors are also conserved for these moves. They label gapped unwanted sectors, not additional critical charges.

[assumed] Exact payload parity requires the grading of the initialized core and every microscopic term to agree with logical fermion parity. The supplied interface preserves such a grading if the core has it. There is no proof of exact U(1)^n flavour conservation for an arbitrary synthesized core. It must not be inferred from the leading quadratic target. The analysis below first permits only marker U(1) and payload Z2. Where U(1)^n is imposed, its additional exclusions are stated.

[proved] The length fermion number r is conserved at fixed M and L, but it is not another independent global charge. Changing M at fixed physical L changes r by -b delta M. Counting both as independent charges would introduce a spurious boson.

## 2. What finite coupling generates

[proved, finite-order perturbation theory] Let P be the complete two-length code manifold, H0 its parent, and Delta the gap to an elementary excluded cell, clock, wrapper, or legality excitation. Include the full physical norms of compensators in V, and denote its local scale by w. With the base energy subtracted,

    H_eff = P V P - P V Q (Q H0 Q)^(-1) Q V P + ...,
    Q=1-P.

The second expression is the second-order degenerate perturbation formula. At fixed order, connected products of local terms give connected effective terms. Two overlapping interface bonds can already act on three cells. Bounds at order k have the form C_k w^k/Delta^(k-1), with constants depending on fixed b, the local alphabet, and the chosen local norm. This statement does not bound an extensive operator norm uniformly in L.

[proved, operator parametrization] A complete way to list the allowed corrections, without pretending to know their coefficients, is to choose a basis O_A of parity-even logical operators and partial length-transfer operators on each finite connected cell interval. Expand every linked-cluster term in that basis. The diagonal subset is

    sum_X F_X(eta_X) O_X,

and the off-diagonal subset is

    sum_(X,eta,eta') |eta'_X><eta_X| tensor O_(X;eta',eta),

with sum eta'=sum eta at fixed M, the physical reference-marker transport included, and Hermitian conjugates paired. Translation makes the coefficients independent of the initial cell label. This parametrization includes every finite-order term; it is not a finite list of numerical coefficients.

[proved, consequences of that parametrization] In particular it includes the following independent mechanisms.

| Microscopic process | Effective correction |
| :-- | :-- |
| Two marker hops, including a virtual wrong length | Length chemical shifts, length density interactions, longer hops, multi-length exchange |
| Two port insertions on one bond | Renormalized quadratic payload terms, density interactions, pairing, general even four- and higher-fermion terms |
| Port insertions on adjacent bonds | Three-cell payload terms and next-neighbor logical hopping |
| A marker hop and a port insertion | Payload-dependent marker amplitude or phase, and length-dependent logical hopping |
| Nonidentical virtual resolvents for b and b+1 | Residual eta dependence despite exact first-order helix equalization |
| Clock/core splitting or approximate initialization | Logical on-cell fields, including flavour mixing and allowed anomalous bilinears |
| Longer connected sequences | Arbitrary symmetry-allowed quasi-local multi-cell terms |
| Time-resolved excursions outside P | Retarded kernels and frequency-dependent vertices |

[heuristic] For a stable gapped eliminated sector, the last kernels decay on time scale Delta^(-1). Expanding at frequency much smaller than Delta gives local operators plus time derivatives suppressed by omega/Delta. This is the usual low-energy integration argument. A general uniform excursion expansion for this moving, highly degenerate core has not been supplied. A finite-order SW remainder of size L epsilon does not control levels of order 1/L.

[assumed, sufficient hierarchy] At each fixed n and b, require w/Delta much smaller than one, a positive dressed cost for every excluded local sector, and a quasi-local correction epsilon_eff small compared with the smallest retained stiffness and the resonance scales below. If W is the payload scale and J the marker hopping, a useful preliminary requirement is

    w^2/Delta + epsilon_core << min(J,W).

It is stronger than w<<Delta. For route A, W can itself be of order lambda^2/Delta. The same-order SW terms intended to synthesize W must be included in the target, and mixed u lambda/Delta terms must be bounded separately. One cannot label every second-order term an error smaller than a second-order target. Small gap denominators in the wrapper and the norms of the helix compensators also belong in this accounting.

## 3. The Gaussian theory and its positivity condition

[assumed, bosonization conventions] For an ordinary squeezed flavour use

    psi_(alpha,s)(y) ~ exp(i p_(alpha,s) y)
                      exp[-i(s phi_alpha-theta_alpha^B)],
    s=+1,-1,
    density_alpha=f_alpha-(partial_y phi_alpha)/pi+oscillations.

The superscript B distinguishes the bosonic dual field from the band-center angle. In physical coordinates let phi_0 count markers,

    y(x)=rho x-phi_0(x)/pi,
    density_marker=rho-(partial_x phi_0)/pi.

All fields are compact. Marker conservation prohibits vertices that change total marker charge.

[proved, from the first-order spectrum] Let v0=2J a sin(pi nu) and v_alpha=2t_alpha a sin(pi f_alpha). At zero band-center angles and for independent periodic XX payloads, define physical charge offsets

    Q0=M-rho L,
    Q_alpha=N_alpha-f_alpha rho L.

The squeezed payload charge offset is Q_alpha-f_alpha Q0. The length offset is -a Q0. The charge energy is consequently

    E_charge = pi/(2L) [v0 a^2 Q0^2
                  + sum_alpha v_alpha (Q_alpha-f_alpha Q0)^2].       (1)

For example, differentiating the exact marker bulk energy gives
f_marker''(rho)=2 pi J a^3 sin(pi nu)=pi v0 a^2, which checks the first term and its normalization.

[proved, canonical transformation] Define

    T = [[a, 0], [-f, I_n]],
    A = T^T diag(v0,v1,...,vn) T,
    B = T^(-1) diag(v0,v1,...,vn) T^(-T).

Then E_charge=pi Q^T A Q/(2L), and the dual winding energy is

    E_current = 2 pi/L [v0/a^2 (D0+sum_alpha f_alpha D_alpha)^2
                         +sum_alpha v_alpha D_alpha^2].             (2)

The signs of individual D_alpha depend on the Slater convention. Equations (1),(2) use the same convention. The eigenvalues of AB are v_alpha^2, including v0^2. Thus the geometric mixing has not created or removed a mode. At n=0 this gives K_marker=1/a^2, in agreement with length_fluid.md. Assigning K_marker=1 would give the wrong charge and current gaps.

[heuristic] General finite-coupling corrections give

    H_G = (1/2 pi) integral dx [phi'^T A phi'
                                  + theta'^T B theta'
                                  +2 phi'^T C theta'].

Linear currents are allowed and shift the equilibrium fillings, band centers, and persistent currents. Quadratic density-density, current-current, and density-current terms are allowed, including every length-payload cross term. They are marginal Gaussian deformations. Spatial total derivatives integrate to zero. Time-reversal or inversion would restrict C, but neither is assumed for the parity-only displaced target.

[proved, quadratic stability criterion] The necessary energetic positivity condition for this Gaussian Hamiltonian is

    G=[[A,C],[C^T,B]] > 0.

Equivalently A>0 and B-C^T A^(-1) C>0. For C=0 the speeds are the positive square roots of the eigenvalues of AB. In general the eigenvalues of [[C^T,B],[A,C]] give n+1 right and n+1 left velocities, with opposite signs. Positivity gives real nonzero velocities. A convenient sufficient smallness test is

    norm(G_ref^(-1/2) delta G G_ref^(-1/2)) < 1.

A bare coupling small compared with the fastest velocity need not satisfy this when a logical mode is extremely slow.

[heuristic] A positive G and absence of surviving non-Gaussian relevant or marginal vertices give n+1 compact bosons. With different velocities this means a product or Gaussian collection of conformal sectors, as allowed by STATUS.md. It is not one Lorentz-invariant theory with a common light cone unless the velocities are also matched.

## 4. Complete potentially dangerous operator classification

### 4.1 Why the small physical length K does not invalidate power counting

[assumed, bosonization locality dictionary] An even payload vertex is labeled by integer chiral charges q_i, i=1,...,2n, with sum_i q_i even. Its cell momentum is p.q. At the free point its payload dimension is

    Delta_pay(q)=sum_i q_i^2/2.

The parity condition excludes odd Jordan-Wigner twist fields from the Hamiltonian. With exact U(1)^n, also impose q_(alpha,R)+q_(alpha,L)=0 for each alpha. Marker conservation sets the marker charge to zero, but permits its density harmonics.

[proved within this dictionary] A physical harmonic dressing has integer label D. Put S=p.q/(2 pi). The first-order momentum gluing gives, up to an overall Fourier sign,

    k_phys = 2 pi (D+S)/a mod 2 pi,
    Delta(D,q) = (D+S)^2/a^2 + sum_i q_i^2/2.                         (3)

To see the field dependence explicitly, use the marker comb

    sum_j delta(y-j)=sum_(D in Z) exp(2 pi i D y).

Multiplication by the payload phase exp(i p.q y), followed by y=rho x-phi_0/pi, gives

    exp[i(2 pi D+p.q)rho x] exp[-2i(D+S)phi_0].

The second factor has dimension K_marker(D+S)^2, with K_marker=1/a^2. The independent squeezed payload factor contributes sum q_i^2/2. This proves (3) under the stated bosonization dictionary. Inserting the payload translation momentum into the current twist in Section 6 gives the same quadratic contribution. The sign change D->-D and q->-q changes neither selection nor dimension.

[proved] Many pure marker harmonics have dimension below two when a is large. A uniform Hamiltonian term, however, must have k_phys=2 pi h for an integer h. On that resonance,

    D+S=h(b+nu),
    p.q-2 pi h nu=0 mod 2 pi,
    Delta_res=h^2+sum_i q_i^2/2.                                     (4)

The integer hb is absorbed by D. Thus the resonant length contribution is h^2, not h^2/a^2. This is the reason that adding an unconstrained free length fermion without the gluing is misleading, while the final low-order resonance list remains finite.

[heuristic] Near a resonance the dimension and mismatch change continuously. At fixed b and n, restrict to a sufficiently small Gaussian neighborhood and keep a strict mismatch margin for every candidate that can approach dimension two. Very high q remain irrelevant because their squared norm is large. Very high physical harmonics remain irrelevant because |D+S|/a is large. No Diophantine condition on every arbitrarily high-order vertex is necessary in this weak-coupling neighborhood.

### 4.2 The list

[proved, free-point charge enumeration] An even integer q has even squared norm. Norms below six are 0, two, and four. Norm two consists of +/-e_i +/-e_j with i different from j. Norm four consists of +/-2e_i or four distinct signed unit vectors. Repeated-index products reduce to these classes or to neutral descendants. The doubled charge is not killed by Pauli exclusion: point splitting produces psi_i partial psi_i, of dimension two.

[heuristic, dimensions use bosonization] Combining this enumeration with (4) gives the complete list at a possible uniform resonance through dimension two.

| Class | Resonance condition in cell momenta | Free dimension | Role |
| :-- | :-- | :-- | :-- |
| Identity | None | 0 | Bulk energy |
| Neutral currents | None | 1 | Shift densities and band centers |
| Neutral quadratic currents | None | 2 | A, B, C and velocities |
| Payload bilinear q=+/-e_i+/-e_j | p.q=0 mod 2 pi | 1 | Mass, mixing, or pairing |
| Bilinear descendant | Same | 2 | Derivative bilinear |
| Doubled payload branch | 2p_i=0 mod 2 pi | 2 | Same-branch derivative pair |
| Four distinct payload branches | sum +/-p_i=0 mod 2 pi | 2 | Non-Gaussian quartic |
| Length backscatter times payload bilinear | p.q +/-2 pi nu=0 mod 2 pi | 2 | New mixed non-Gaussian channel |
| Pure length h=+/-1 | nu=0 or 1 mod 1 | 1 | Endpoint, outside the interior fluid |
| Pure length | h nu integer | h^2 | h>=2 is irrelevant near the free point |
| All remaining even sectors | As in (4) | At least 3 | Irrelevant in a small neighborhood |

[proved] For number-conserving flavour backscattering the new familiar conditions include

    2 pi nu +/-2 k_alpha = 0 mod 2 pi.

Without flavour U(1), also include

    +/-2 pi nu +/-p_i +/-p_j = 0 mod 2 pi,

with both independent signs and distinct branch indices. These include mixed pairing processes, not just differences of 2k_F. Exact marker U(1) excludes creation of one length particle paired with one payload particle. The permitted length bilinear has one creation and one annihilation branch, hence momentum +/-2 pi nu.

[proved, obstruction to a weaker hypothesis] Distinct incommensurate widths k_alpha do not alone forbid uniform pairing. For an unshifted band theta_alpha=0, psi_R psi_L carries p_R+p_L=0 and bosonizes to exp(2i theta_alpha^B). It has dimension one at the free point. A parity-even nearest-neighbor p-wave pairing term is allowed and gaps that payload. The supplied ir_control analysis already avoids this by displacing the band centers. That qualification must survive in the mobile construction. It is not an objection to the fully displaced target.

[numerical] additional_checks.py tests this elementary obstruction for k_F=pi sqrt(2)/4, hopping one, and pair amplitude Delta=0.02. It samples 200001 momenta in

    E(k)^2=(-2 cos k+2 cos k_F)^2+4 Delta^2 sin^2 k.

The sampled minimum is 0.0358390 (in hopping units), recorded in additional_results.json. This is a free BdG check, not an interacting gap estimate.

### 4.3 Extended Sidon margins and an open nonresonant set

[proved, using the supplied Sidon construction] Retain the payload margins from ir_control,

    delta2=min distance(2p_i, +/-p_i+/-p_j; 2 pi Z),
    delta4=min distance(sum of four distinct signed p_i; 2 pi Z).

Its explicit scaffold gives delta2>=pi/(6n), delta4>=pi/(576n^2) for the stated n range, before the arbitrarily small irrational displacement. Do not recompute or replace those already established constants. Add

    delta_mix=min_(i<j,s,t,u=+/-1)
                    distance(s p_i+t p_j+u 2 pi nu; 2 pi Z).         (5)

Pure length endpoints are separately excluded by nu in [1/4,3/4].

[proved, elementary extension] Put N=2n. At fixed payload momenta, the forbidden nu values in a unit interval lie in a set of at most K=4 binomial(N,2) points. Split [1/4,3/4] at these points. One interval has length at least 1/[2(K+1)]. Its midpoint has angular distance at least pi/[2(K+1)] from every forbidden value. A sufficiently small irrational displacement preserves the conservative margin

    delta_mix >= pi/[4(K+1)].                                      (6)

The endpoints of the chosen interval protect it even from forbidden points outside [1/4,3/4]. Therefore adding the length mode costs another order n^(-2) margin, without reducing the payload order n^(-1) quadratic spacing. This construction respects length number conservation; applying an unrestricted Sidon test to both +/-pi nu would incorrectly reject their symmetry-forbidden zero-momentum pair.

[numerical] additional_results.json explicitly enumerates the signed constraints for n=1,2,4,8. The resulting (delta2,delta4,delta_mix) are approximately

| n | delta2 | delta4 | delta_mix |
| :-- | :-- | :-- | :-- |
| 1 | 1.374447 | Not applicable | 0.719948 |
| 2 | 0.698132 | 0.0232711 | 0.349066 |
| 4 | 0.317689 | 0.00588313 | 0.1000132 |
| 8 | 0.122856 | 0.00195009 | 0.0438770 |

These are floating-point checks of the explicit algorithm, not interval-certified minima. The analytic interval argument is the existence proof. Small irrational shifts remove exact rational scaffold relations while retaining the finite margins.

[heuristic, RG smallness] A vertex of dimension Delta and initial dimensionless strength y grows as

    y(ell)=y(0) exp[(2-Delta) ell]

until its oscillation wavelength is reached. In physical units the mismatch is the cell mismatch divided by a; the physical velocity is the cell velocity times a. Thus the mismatch energy is still of order w_min delta, with w_min the smallest cell-coordinate velocity. A sufficient criterion is

    |y(0)| (E_UV/E_mismatch)^max(0,2-Delta) << 1.

For quadratic masses this reduces to epsilon_eff << w_min delta2. For near-marginal vertices also control their running before cutoff, for example g log(E_UV/E_mismatch)<<1 together with the small scaling-dimension shift. Dressed momenta must retain (5), not only the bare ones. Their changes must be smaller than the chosen margins, or the counterterms must retune them.

[heuristic] These strict inequalities and G>0 define an open sufficient Gaussian region at every fixed finite n, assuming the dressed parameters depend continuously on microscopic couplings. The region can shrink with n, b, the core gap, and the slowest logical velocity. No uniform neighborhood at n tending to infinity is claimed.

## 5. Ways to gap or destabilize the length sector

[heuristic, commensurate locking] At nu=p/q in lowest terms, the first pure length umklapp is cos(2q phi_gap), with dimension q^2 K_gap near the uncoupled reference. At the free reference K_gap=1. Interior q>=2 umklapps are irrelevant; at half filling the dimension is four, becoming relevant only when the appropriate dressed dimension crosses two. The physical marker density is rho=q/(bq+p). Its smallest lattice-invariant harmonic has index bq+p and dimension

    (bq+p)^2 K_marker = q^2,

which is a consistency check on Section 4. A large b does not make this allowed umklapp relevant at weak coupling. With payload mixing, use the full vertex metric rather than one scalar K.

[heuristic, mixed locking] If nu matches a payload filling, the mixed operator can be marginal already at the free point. Its sine-Gordon equations have the form dy/dell=(2-x)y and dx/dell=-C y^2, C>0 after normalization. There are arbitrarily weak couplings on the relevant side that pin a relative field and remove one boson. Distinct payload fillings do not exclude a coincidence with nu. Avoid (5), or analyze the particular resonant interaction and its sign. An attractive spin-charge interpretation is not a general protection.

[proved, Gaussian thermodynamic criterion] Phase separation corresponds to loss of convexity of the bulk energy as a function of the conserved physical densities. At the free point the length curvature is positive,

    f''(rho)=2 pi J a^3 sin(pi nu)>0.

With payloads, positivity must hold for the full density Hessian, including its Schur complement after the payload densities relax. A negative Schur complement cannot be repaired by checking the isolated marker curvature. It is the thermodynamic part of G>0. Strict positivity at a fixed interior reference persists under sufficiently small continuous corrections.

[heuristic, endpoints and feedback] Keep nu and every f_alpha bounded away from zero and one, maintain the positive Hessian, and keep the velocities nonzero. At a band edge the linearization fails and a dilute quadratic mode replaces a boson. Payload ordering can feed back through F(eta) O_pay or density-assisted hopping. If a resonant payload vertex orders, its expectation can supply the missing momentum to a length locking term. The sufficient conditions here exclude that first ordering instability by the complete parity-even operator test; the absence of a bare length potential alone would not exclude feedback.

## 6. Towers, spin structures, and the ring constraint

### 6.1 Exact kinematic statement, including interactions

[proved] The reference cover is

    (x,eta,z) ~ (x+b+eta_0,S eta,S z),
    x in Z_L, sum eta=r.

After Fourier transforming x at P=2 pi ell/L and applying the gap gauge of joint_fluid.md, the simultaneous cell shift must have eigenvalue exp(-2 pi i ell/M). Thus the exact interacting replacement of equations (3),(4) there is

    Spec H_L = union_(M,ell) Spec H_sq(M,r;P)
                    restricted to m_total+ell=0 mod M.              (7)

Here H_sq includes all retained length-payload interactions and the flux P induced by physical marker transport. Only m_total is conserved in general. Separate m_gap and m_payload cease to be quantum numbers. Equation (7) is a Hilbert-space quotient, so adding correlated hops does not delete it. The additive cosine energy formula, unlike (7), does not survive generic interactions.

### 6.2 An explicit free zero-mode lattice

[proved, for the number-conserving reference] The following gives the lattice including displaced bands and the usual fermion spin structures. Let beta_alpha be the boundary shift of flavour alpha in units of 2 pi/M. For independent periodic XX payloads beta_alpha=(N_alpha-1)/2 mod 1. For n local fermion flavours obtained with one common Jordan-Wigner grading,

    beta_alpha=(N_total-1)/2 mod 1,
    N_total=sum_alpha N_alpha.

An externally imposed spin structure would change this shift and must be specified. The microscopic tensor-product ring and grading fix it; it is not a free extra choice.

[proved] A consecutive sea of N_alpha orbitals has center

    c_alpha in Z+beta_alpha+(N_alpha-1)/2,
    d_alpha=c_alpha-M theta_alpha/(2 pi).

The length sea has center j0 in Z. Its apparent fermion boundary shift (r-1)/2 cancels the half-integer sea-center shift, explaining why there is no additional restriction on even or odd r in the XX length sector. Let K_osc be the signed total cell momentum of the particle-hole partitions, in orbital units. Then the reference constraint is

    r j0+sum_alpha N_alpha c_alpha+K_osc+ell=0 mod M.                 (8)

The total on the left excluding ell is integral, including in the common spin structure. Equations (8), r=L-bM, and integer N_alpha are a direct zero-mode and oscillator specification without suppressing any ring congruence.

[proved] To display the quadratic current energy, choose D_alpha=-d_alpha and D0=b j0+t, where t is the integer lifting the congruence. Define

    chi0=sum_alpha N_alpha theta_alpha/(2 pi),
    D_alpha in Z+M theta_alpha/(2 pi)
                         -beta_alpha-(N_alpha-1)/2,
    D0 in Z,
    W0=[D0-chi0+sum_alpha (N_alpha/M)D_alpha-K_osc/M]/a.              (9)

The length energy involves W0^2; the payload energies involve D_alpha^2. The physical momentum is P=2 pi W0 mod 2 pi. Expanding the exact filled-sea cosines gives

    E_current = (2 pi/L)[v0 W0^2+sum_alpha v_alpha D_alpha^2].        (10)

At zero band centers and independent XX spin structures, this reduces to (2) at leading order. The oscillator contribution in W0 shifts current energies only beyond leading order for bounded winding labels, but is needed to assign momenta. Terms involving O(1) density offsets in N_alpha/M give the standard charge-current corrections to P.

[numerical, tower_results.json] A separate filled-sea test uses n=1,2, M=40,80,160,320, b=4, band-center angles (0.31,-0.57), and both independent and common Jordan-Wigner spin structures. All integer congruences in (8) hold, and the two expressions for W0 agree within 8.4e-17. The exact cosine energy minus its Casimir-plus-current expansion is of order L^(-3) with bounded, shell-dependent coefficients. At M=320 its magnitude is at most 5.64e-6 in the tested states. This checks the finite-size expansion and the algebraic gluing, not a full microscopic diagonalization.

[proved] No universal restriction on the parity of n follows from this construction. The relevant parity is N_total for a common fermionic grading, or N_alpha for independent spin-chain gradings. Nor is L mod b a surviving crystallographic obstruction. It only determines which integers M satisfy 0<=L-bM<=M. For an interior density, such sectors occur for every sufficiently large L.

[heuristic, parity-only continuation] If flavour U(1)^n is broken microscopically but all its number-changing vertices oscillate or are irrelevant, N_alpha become emergent low-energy labels rather than exact finite-volume conserved charges. The fixed microscopic parity projection and the compactification lattice remain. Equation (7) remains kinematic; assigning exact separate N_alpha to every interacting eigenstate would be incorrect. Their tower description is the RG prediction. The core grading and possible extra exact symmetries must be supplied to specify the microscopic modular invariant completely.

### 6.3 Energies and central charge

[heuristic, Gaussian spectrum] With the equilibrium charge and current offsets retained, the general finite-size result is

    E = L e_infinity
        -pi/(12L) sum_j(u_(j,R)+u_(j,L))
        + E_zero(Q,D;offsets,G)
        +(2 pi/L) sum_j[u_(j,R) N_(j,R)+u_(j,L) N_(j,L)]
        +o(1/L).                                                   (11)

There are n+1 oscillator families on each side; N_(j,s) are nonnegative partition levels. The zero-mode quadratic form is obtained by integrating H_G with the lattice (7)-(9). At C=0 it is pi Q^T A Q/(2L)+2 pi D^T B D/L, with the twists included. Irrelevant interactions add subleading corrections; residual marginal running, if present at a resonance, must not be hidden in the o(1/L) notation.

[proved, reference consistency check] At n=0, theta=0, this is the length-fluid formula with K=rho^2. At independent XX payloads, zero band centers and compatible equilibrium counts, it gives

    E0 = kappa M - sum_(j=0)^n 2t_j sin(pi N_j/M)/sin(pi/M),
    N0=r,

and expanding 1/sin(pi/M)=M/pi+pi/(6M)+O(M^(-3)) gives -pi sum_j v_j/(6L). These checks fix the factors of two in (1),(2),(10),(11).

[heuristic, how to measure c] A Casimir coefficient with unequal velocities measures sum_j c_j v_j, not c_total times an arbitrary average velocity. Determine the speeds independently, resolve their oscillator branches, and compare each contribution. Alternatively interval entropy has coefficient c_total/3, and the low-temperature free energy has coefficient -pi T^2 sum_j(1/v_j)/6 when left and right speeds match. For unequal chiral speeds replace 1/v by (1/u_R+1/u_L)/2. These diagnostics require removal of charge-rounding and twist energies before a vacuum Casimir fit. Here each stable compact boson has c_j=1, so c_total=n+1.

## 7. What uniqueness does and does not mean

[proved, a sufficient subclass] The reference physical marker graph is connected whenever 0<r<M. Open-cell exchanges connect all gap strings, and a reference-marker move changes its origin by one site. Independent payload exchanges connect each fixed payload-charge sector. Their combined graph is connected. If every allowed off-diagonal amplitude is strictly negative and all additional terms are real diagonal potentials, Perron-Frobenius therefore gives a unique ground vector in every fixed-charge sector. Positive density-dependent hopping coefficients preserve this conclusion. The graph model in Section 8 satisfies these conditions.

[proved] Add a chemical potential vector multiplying the conserved charges. Within one fixed-charge sector it only shifts the spectrum by a scalar. Equality of two distinct sector ground energies is an affine hyperplane in that vector. Avoid the countable union of such hyperplanes over integer L and charges. Then the global finite-ring ground state is unique for every L for which its winning sector has the connected interior graph. This is an existence and genericity statement, not a proof that a particular decimal parameter avoids every crossing. Thermodynamic convexity and an interior minimizer place the winner in the interior for sufficiently large L in the controlled phase.

[heuristic] A displaced parity-only target generally has complex amplitudes, so this positivity proof is unavailable. One-site translation and Z2 have only one-dimensional irreducible representations and do not require a degeneracy, but that fact is not a uniqueness proof. Exact crossings between different momentum or parity sectors can remain. Generic parameters are expected to avoid their countable crossing loci when no additional antiunitary or projective symmetry enforces degeneracy. The actual microscopic core must be checked for those symmetries. Fine-tuned chemical potentials can fail even in the soluble first-order model.

[proved, strict all-L vacuum obstruction] Even a single length boson illustrates a stronger issue. Write x_L=rho L, and let m be its nearest integer. At leading order the competing charge energies are

    E(m)-E_bulk = A(m-x_L)^2/L + common Casimir,
    A=pi v/(2K)>0.

If frac(x_L) lies just below 1/2, the next charge sector has gap

    Delta E = A[1-2 frac(x_L)]/L + o(1/L).

For irrational rho, the fractional parts of rho L are dense. There are subsequences with L Delta E tending to zero, although every finite ring can have a nondegenerate vacuum. The o(1/L) correction does not remove this subsequence conclusion. The rescaled zero-mode spectrum also moves on other subsequences.

[proved, interpretation] Thus generic irrational density does not give convergence of the entire rescaled spectrum, along every integer L, to one fixed untwisted compact CFT with a uniformly isolated identity state. It gives CFT towers in size-dependent background charge and flux sectors. Fixing a subsequence or recording the twists repairs the comparison. If “strict sense” in the question excludes these moving sectors, the proposed generic incommensurate construction fails that stronger definition already at first order. If it means finite-L uniqueness plus the usual Luttinger-liquid spectrum with its specified offsets, the conditional positive result above is the appropriate statement. These two meanings should not be conflated.

## 8. Numerical tests and reproducibility

### 8.1 Hamiltonians and scope

[proved, explicit toy definition] The numerical model is the physical projected graph, with cells of lengths b and b+1 and n independent hard-core payload species attached to the ordered cells. It retains the joint reference-marker quotient exactly. Its unperturbed squeezed Hamiltonian is

    H0=-sum_(j,s=0,...,n) t_s(c^dag_(j,s)c_(j+1,s)+h.c.),
    (t0,t1,t2)=(1,0.73,0.51),

where species zero is the length bit eta. Payload hopping is the periodic XX hopping of each species. These are independently Jordan-Wigner fermionized species with exact U(1) charges. Their global spin structures differ from a parity-only n-flavour core with one common Jordan-Wigner string. This is an explicit test of a symmetry-protected subclass, not a numerical test of all parity-allowed operators.

[proved] The generic graph deformation adds

    0.05 sum_j eta_j eta_(j+1)
    +0.04 sum_(j,alpha) n_(j,alpha)n_(j+1,alpha)
    +sum_(j,alpha) [0.07/alpha](eta_j-nu_ref)(n_(j,alpha)-f_(alpha,ref))
    +0.035 sum_j n_(j,1)n_(j,2)                 [only n=2].

Alpha starts at one in this expression. A marker hop across cells j-1,j is multiplied by

    1+0.025 sum_alpha [n_(j-1,alpha)+n_(j,alpha)-2 f_(alpha,ref)]/alpha,

and a payload hop across j,j+1 by

    1+0.03(eta_j+eta_(j+1)-2 nu_ref).

These amplitudes stay positive. The centering constants are fixed at nu_ref=0.37, f_(1,ref)=0.27, f_(2,ref)=0.16 on every size. Thus the generic ED uses one fixed projected Hamiltonian per n, in the stated charge sectors. No microscopic wrapper assertion is inferred from that sequence.

[proved, soluble finite-interaction control] Also define

    U_gamma=exp[i 0.37 sum_(j,alpha) eta_j n_(j,alpha)/alpha],
    H_gamma=U_gamma H0 U_gamma^dagger.

This produces length-dependent payload hopping phases and payload-dependent marker phases at finite strength. It is an interacting Hamiltonian in the original variables, exactly isospectral to H0. In cell coordinates U_gamma is on-site and preserves interval entropy. In physical coordinates it is a finite-depth diagonal circuit of bounded range for fixed b, since eta_j tests whether the next marker is b or b+1 sites away. Its physical interval entropy differs from H0 by a bounded boundary contribution, so it preserves an existing logarithmic coefficient. This is a special, soluble direction in interaction space. It does not substitute for independent generic interactions or leakage.

### 8.2 Two accurate calibrated estimators

[numerical, controls.json] For n=0,1,2, use the fixed irrational target fillings

    (nu,f1,f2)=(sqrt(2)/4,sqrt(3)/7,sqrt(5)/11),

rounded to integer counts at M=40,64,96,128,192, with b=4. The first estimator fits the sum of free correlation-matrix interval entropies over M/4<=ell<=M/2 to the periodic chord logarithm. It is a cell-coordinate entropy estimator. The second uses the exact physical energy and the independently known velocities, subtracting the bulk energy at the actual finite-size densities. This subtraction removes the charge-shell energy rather than mistaking it for a Casimir term. The reported weighted estimator is

    c_Cas=(n+1)[-6L(E0-E_bulk)/(pi sum_j v_j)].

It tests whether the weighted coefficient is the predicted sum of n+1 unit-charge sectors. It would not by itself determine the total c of an unknown unequal-velocity model.

| Payload n | M | Physical L | c from cell entropy | Weighted c from Casimir |
| :-- | :-- | :-- | :-- | :-- |
| 0 | 192 | 836 | 1.0000812 | 1.0000312 |
| 1 | 192 | 836 | 2.0003014 | 2.0000625 |
| 2 | 192 | 836 | 3.0006723 | 3.0000937 |

[proved plus numerical] The same two values apply to H_gamma: its spectrum is unchanged, and its cell-interval entropy is unchanged. The construction of the physical graph and the physical entropy tests below check the distinction between this statement and claiming that squeezed and physical entropies are identical. They are not identical.

[numerical, controls.json] At fixed chemical potentials, the exact free grand-canonical ground energy was minimized over every allowed M and every payload count on every ring L=40,...,400, for n=1 and n=2. No sector degeneracy occurred, all winners had 0<r<M, and the smallest difference between the two best charge sectors in this range was 8.80269e-5. This includes all primes in the interval. Perron-Frobenius supplies uniqueness inside the winning sector. These conclusions also hold for H_gamma by unitary conjugation. This is a finite range check, not a test of all integers L.

### 8.3 Physical-graph ED with independent interactions

[numerical, graph_results.json] The graph calculation uses b=2 and (M,r)=(6,2),(7,3),(8,3), giving physical L=14,17,19. Payload counts are N1=2 and, for n=2, N2=1. The largest Hilbert space has dimension 29792. The free Slater-product wavefunction, lifted with the factor sqrt(M/L), has Hamiltonian residual below 2.1e-15 in the final data. Reduced density matrices are constructed in the physical site basis, with empty sites and payload-bearing markers, and are blocked by charges in the physical interval. No gap-chain entropy is substituted here.

[numerical, matched controls] The physical-site chord fits give the following. The substantial offset already present in the free calibration makes an accurate interacting c determination impossible at these sizes.

| n | L | Free physical c fit | Interacting physical c fit |
| :-- | :-- | :-- | :-- |
| 1 | 14 | 2.60478 | 2.60680 |
| 1 | 17 | 2.50834 | 2.50765 |
| 1 | 19 | 2.44074 | 2.44138 |
| 2 | 14 | 3.75100 | 3.75400 |
| 2 | 17 | 3.66349 | 3.66345 |
| 2 | 19 | 3.56429 | 3.56558 |

[numerical] Renyi-2 fits and the physical entropy after U_gamma are retained as additional diagnostics in the JSON. They are not advertised as independent successful determinations of the generic interacting central charge. There are too few particles and too few cells for that claim. The matched entropy changes are small, which is consistent with weak Gaussian deformation but cannot exclude a long crossover to a small gap.

[numerical, contiguous_rings.json] Separate ED covers every physical ring L=14,...,19, including primes 17 and 19, with M=round(L/2.4), N1=2, and N2=1 when present. It searches the full fixed-charge physical graph, including all momenta. Ground states are unique, as independently guaranteed by graph connectivity and the positive hopping amplitudes. The products L times the first gap lie between 4.00 and 5.21. This short range does not determine a dynamical exponent. These are fixed-charge tests; no claim of global interacting-sector minimization is made.

[numerical-method correction] An initial Krylov vector equal to the translation-invariant free ground state missed a nonzero-momentum excited sector in one repeat. Both scripts were corrected to use a fixed-seed random vector before regenerating the data. A translation-invariant starting vector is unsuitable for a complete low-energy spectrum search. This correction affects the gap diagnostic, not the exact free reference formulas.

### 8.4 Reproduction and limits

[proved, file inventory] Run from this folder with Python 3, numpy, and scipy installed:

    python3 test_fluids.py
    python3 additional_checks.py

Both set single-thread BLAS before importing numpy. The first writes controls.json and graph_results.json. The second writes additional_results.json and contiguous_rings.json. An additional direct check of the displaced-band current lattice is reproduced with `python3 tower_checks.py`, which writes tower_results.json. To rerun only the large analytic/correlation-matrix controls:

    python3 test_fluids.py --controls-only

[numerical] Individual graph calculations take seconds to tens of seconds at the stated sizes on this machine. Each complete suite takes minutes. No network, DMRG package, cluster job, or computation outside this folder was used. No original wrapper, compiler, or finite-gap cell excitation was numerically simulated. A two-estimator test of that microscopic finite-coupling chain remains undone because that chain has not been specified and the small projected-graph tests are not accurate enough to stand in for it.

## 9. The large-n accumulation point

[proved] Hermitian bonds on a fixed finite-dimensional two-site Hilbert space form a finite-dimensional vector space. Only a bounded sequence has a guaranteed convergent subsequence. Normalize the bonds, for example by norm after removing an additive scalar, before invoking compactness. The supplied wrapper parameters depend strongly on b, so this normalization changes the comparison of all physical energy scales.

[proved, bare selector asymptotics] The source gives the crossing

    mu_c=2t-delta_b,
    delta_b=(b+1)a_b-b a_(b+1),
    a_b=2t[1-cos(pi/b)].

Expanding the displayed formula, rather than changing its convention, gives

    a_b=pi^2 t/b^2+O(t/b^4),
    delta_b=3 pi^2 t/b^2+O(t/b^3),
    e_cross=a_(b+1)-a_b=-2 pi^2 t/b^3+O(t/b^4).

Thus b tends to infinity as mu tends to 2t from below. At the threshold, a long interval costs

    E_l=2t-2t cos(pi/l)=pi^2 t/l^2+O(t/l^4).

The head dispersion near its minimum is t k^2. This is the candidate z=2 mechanism in the bare selector. The selected-cell density tends to zero, and the scale separating nearby cell lengths also closes, of order t/b^4 in the frozen crossing analysis.

[heuristic] A controlled encoding must then send its logical couplings below the closing cell scales, often much further below them. Its physical logical velocities, roughly cell length times logical hopping, can vanish. If the normalized family retains the selector kinetics while the interface and code scales vanish, the limit is consistent with a dilute z=2 threshold plus collapsed logical scales. It is not an n=infinity compact CFT with finite nonzero velocities.

[proved, distinction between constructions] The exact equal-hopping d=4 selector is a different solvable family. Its free wall dispersion is -2a cos k+mu/2, so its dilute threshold is |mu|=4a, not mu=2t. The latter belongs to the long-cell small-marker-hopping construction. These thresholds cannot be identified by transferring formulas from one family to the other.

[assumed or unresolved] The complete fixed-alphabet high-channel bond h_n, its normalization, and its convergent subsequence have not been supplied. Consequently the full accumulation Hamiltonian cannot be identified from the present notes. Depending on the relative normalization, its dominant terms could instead leave an extensive zero-energy code space, static marker sectors, or more than one vanishing dynamical scale. The z=2 selector calculation is evidence for a natural candidate, not a proof that every accumulation point has only z=2 dispersion.

[proved, logical limitation of the c-theorem argument] The continuum c-theorem compares endpoints of a suitable unitary relativistic RG flow. Closeness of two microscopic bonds does not itself show that the second chain flows out of the first chain's IR CFT. The perturbation can act above an onset length that diverges along the sequence. Therefore “c is locally bounded near any lattice bond with a finite-c IR CFT” needs an additional uniform crossover or basin hypothesis. Under that hypothesis a finite-c conformal accumulation point would indeed be excluded. Without it, neither nonconformality nor its particular exponent follows from the c-theorem alone.

## 10. Verdict and ranked next steps

[heuristic, positive conditional result] At fixed finite n and b, exact marker conservation, exact payload parity, displaced and jointly nonresonant Fermi points, a locally controlled gapped-sector elimination, and positive Gaussian stiffness give a consistent n+1-boson phase. The additional length mode does not force a relevant mixed resonance. The explicit extension (5),(6) provides an open weak-coupling parameter region at the level of standard bosonization.

[proved and unresolved] This does not establish the full question for the supplied microscopic candidate. A complete high-channel bond is absent, its dressed coefficients and symmetries have not been calculated, and generic physical-graph ED at the available sizes does not meet the two-accurate-estimator standard. Finite-ring uniqueness is proved in a positive-hopping subclass and checked over the stated ranges, but it is not proved for the parity-only microscopic core. A fixed untwisted vacuum-sector limit on every integer L is already obstructed by generic irrational rounding.

[heuristic, ranked continuation]

1. Specify one complete initialized mobile core and its grading, including all finite coefficients. Compute its two- and three-cell second-order effective terms, with the actual smallest gap and compensator norms. This is the input needed to decide whether the conditional region contains the intended microscopic family.
2. Retune the dressed payload momenta and nu using the extended margins. Check the full density/current stiffness matrix, including its slowest eigenvalue, instead of only checking individual velocities or isolated compressibilities.
3. Test the explicit interacting family at enough cells for the free physical-entropy calibration to be accurate. Use entropy together with independently measured velocities and shell-corrected Casimir or low-temperature free energy. The present small graphs do not justify a central-charge claim for generic corrections.
4. State the all-ring convention. If size-dependent charge and spin-structure sectors are allowed, retain (7)-(11) in the definition. If one fixed identity sector with a uniformly isolated vacuum is required, a different filling/sector mechanism is needed.
5. Only after a fixed-d bond family is explicit, normalize it and determine its accumulation Hamiltonian. The bare-selector z=2 limit is a useful guide, not an identified limit of the unfinished compiler.

[proved, scope of this run] PROGRESS.md records the resumable state. All new code, data, and this report are in work_C13. Project STATE.md and every source file outside this folder were left unchanged under the task's restriction. No source-draft claim was promoted to an independently checked microscopic theorem.
