# C41. One vacancy, a period-b payload, and the all-ring question

[heuristic, verdict] The one-vacancy lift has a consistent weak-coupling route to n+1 gapless Gaussian modes at each certified fixed n, including n=23 and d=24. The period-b order does not supply an unavoidable mixed instability. Its b offsets become charge flux sectors, and adjusting the occupied number removes the incompatible-ring defect at a cost of order 1/L. An especially useful sufficient density window is 2/5 < rho < 3/5, where rho is the occupied density. In that window the fractional marker harmonics are harmless after physical translation is imposed, under the quantitative Gaussian conditions below.

[assumed, qualification] This conclusion uses the payload infrared transfer and a weak-coupling local matching assumption for the vacancy lift. Neither follows from the finite-ring projection alone. The supplied certificates, now independently rerun by the campaign route_a worker at n=23, prove the quadratic reference and finite corrections through n=23, not the full interacting central charge. This report does not turn those certificates into an unconditional microscopic theorem or a central-charge measurement.

[proved, qualification] The stronger assertion in the question, one fixed CFT vacuum and a proved unique ground state on every large ring, does not follow. The separated reference already has charge-shell offsets that vary with L and can have flux degeneracies. The full complex bond has no supplied positivity proof of uniqueness. Unequal velocities also preclude a single common Lorentz velocity without further tuning. The positive conclusion concerns one glued Gaussian theory with n+1 oscillator families and size-dependent charge/flux sectors, in the convention used in C13, C20, C27 and C40.

## 1. Inputs and what is inherited

[proved, source record] I read `../../STATUS.md`, the campaign README and its doping, doping_check, route_a and sparse_core reports, the original 23-letter specification, the mobile-cell toy report, C13, C20, C27, C33, C37 and C40. The requested local `doping/README.md` is absent. Its intended content is in `../../contrib/counterexample_campaign/doping/README.md`. The campaign infrared report, its charge-threshold proof, the doping_check locality and threshold reviews, and the route_a C33 audit appeared during the final check and were then read. Sections 3-5 and the local numerical checks were derived before those reports appeared. Their continuum dictionary agrees with (16). The charge-threshold result added to Section 6 is attributed to their new proof and checked algebraically here.

[proved, inherited finite results] C33 certifies an n=23 quadratic reference with 46 simple roots, delta_2 > 0.0140193 and delta_4 > 3.1558e-7, including finite-conversion corrections. C37 certifies the corresponding existence statement for every n=1,...,23, with a different parameter schedule. These schedules must not be mixed. The campaign route_a proof supplies all-length internal-cell bounds and excludes wrong-length and illegal sectors on compatible frozen rings. For its n=23 parameters, 0 < |lambda| <= 2^-20 and K=10 suffice for its global-sector bound. Infrared control can require a smaller lambda. These source certificates were read, not rerun.

[assumed, payload input P] Fix one such original bond h, including all its small but nonzero couplings, inside its interacting Gaussian basin. Let H_s(M)=sum_j h_(j,j+1). On M=bN occupied sites its selected sector is a ring of N cells with n gapless fermion flavors, b=n+1, exact total payload U(1), simple roots p_i in cell units, and strict quadratic/quartic margins. Let Delta_s>0 denote a reserved defect/core cost in h units. Let v_*^cell>0 denote a reserved minimum payload speed. The claim that these are the dressed infrared data, with uniform local control, is C40's thermodynamic input, not an additional finite certificate.

[proved, definition] Write n_x=1 on the 23 occupied letters and zero on the added vacancy. The proposed Hamiltonian is

    H_24 = -t sum_(x,a) (|0,a><a,0| + h.c.)
           + J sum_x h_occ(x,x+1) - mu sum_x n_x,
    t>0, J>0.

[proved] This has d=24, one repeated nearest-neighbor Hermitian bond, occupied-number U(1), the inherited payload U(1), and marker-number U(1). Each assertion follows term by term: swaps transport the whole occupied letter, h acts only inside its original block, and the chemical term is diagonal. The entire h, including its half-site terms, must be lifted into the occupied-pair block. Lifting just its off-diagonal terms would give a different first-order internal chain. The vacancy density is 1-rho. It is finite here, not a single hole on a growing ring.

## 2. Exact squeezing and the first-order operator

[proved] Fix 2 <= M < L. Every physical basis state is uniquely specified by occupied positions x_1<...<x_M and their ordered internal word sigma. This is a unitary identification with

    C^(choose(L,M)) tensor (C^23)^tensor M.

[proved] Let R rotate the internal word to the right. A hop inside the chosen physical cut leaves that word unchanged. A hop L-1 to 0 applies R, and its inverse applies R^-1. If B includes the coefficient -t of that cut hop, the exact kinetic operator is

    T_bulk tensor I + B tensor R + B^dagger tensor R^-1.       (1)

[proved] An occupied physical bond joins consecutive squeezed labels. Set A_j(x)=1 when x_(j+1)=x_j+1, j<M, and A_M(x)=1 when x_1=0 and x_M=L-1. With oriented cyclic internal bonds,

    V = sum_j A_j tensor h_(j,j+1).                           (2)

[proved] Equations (1)-(2) follow on every basis vector, so are operator identities. No charge averaging, large-L limit, or assumption about h enters them. The original marker pattern and arrow violations, viewed cyclically in the occupied word, are conserved by these moves. A boundary hop changes their reference origin, not their cyclic arrangement. Individual physical bonds need not remain occupied or legal.

[proved] In the R=e^(i theta) sector the charge problem is a hard-core hopping ring with flux theta. Fermionizing it adds the sign (-1)^(M-1), because c_0^dagger c_(L-1) must pass M-1 occupied orbitals. A consistent orbital convention is

    k_a = [2pi I_a - theta - pi(M-1)]/L,
    E_c = -2t sum_a cos(k_a).                                 (3)

[proved] For consecutive I_a=ell,...,ell+M-1, summing the geometric series gives

    E_c(theta,ell) = -2t sin(pi M/L)/sin(pi/L)
                         * cos[(2pi ell-theta)/L].             (4)

[proved] The ground sea chooses the representative 2pi ell-theta of minimum absolute value. Both representatives are retained at theta=pi. At J=0 the whole internal Hilbert space is degenerate apart from these charge twists. In particular J=0 itself is not a c=n+1 theory. The order of limits is fixed J>0, then large L and low energy.

[proved] Let Pi_r project onto R=e^(2pi i r/M), and retain all minimizing seas in each cyclic sector. Their physical translation phases, with the convention of (3), are exp[2pi i(r-M ell)/L]. The integers r-M ell are distinct and lie between -M/2 and M/2, including the two half-flux endpoints. Since M<L, their phases are distinct. Translation invariance therefore makes V block diagonal between these retained sea branches.

[proved] Inside one cyclic block,

    Pi_r h_j Pi_r = Pi_r H_s Pi_r/M.

To prove it, conjugate h_j by powers of R and use the same scalar R eigenvalue on the two sides of the matrix element, then sum over j. Thus individual cut-dependent expectations of A_j need not coincide. Only their sum is used.

[proved] In a consecutive Slater sea the nearest-neighbor contraction has modulus sin(pi M/L)/[L sin(pi/L)]. Wick contraction of n_x n_(x+1), followed by summing the L bonds, gives

    alpha_(L,M) = M/L - sin^2(pi M/L)/[M L sin^2(pi/L)],
    P_sea V P_sea = direct_sum_(r,ell) alpha_(L,M) Pi_r H_s Pi_r.
                                                                  (5)

[proved] At fixed 0<rho=M/L<1,

    alpha(rho) = rho - sin^2(pi rho)/(pi^2 rho) > 0,
    v_c = 2t sin(pi rho),
    e_c(rho) = -2t sin(pi rho)/pi.                              (6)

The strict positivity follows from sin(pi rho)<pi rho. At half filling alpha=1/2-2/pi^2. Near rho=0, alpha=pi^2 rho^3/3+O(rho^5), so the dilute limit loses the desired uniform internal scale. At rho=1 the charge velocity vanishes. Both endpoints are excluded.

[proved, limitation] Equation (5) is an exact compression. It is not a uniform-in-L Schrieffer-Wolff reduction, since charge particle-hole gaps vanish as t/L. Also its different cyclic branches are not one exactly degenerate eigenspace. Calling it the first-order reference is legitimate; estimating the omitted theory with a denominator t/L at fixed J is not.

## 3. What becomes of the b offsets

[proved] For M=bN, the selected marker patterns form b mutually orthogonal offset subspaces of the squeezed internal Hilbert space. R permutes these subspaces cyclically. On returning after b applications, it translates the N-cell payload ring by one cell. Write that cell translation as S and choose a payload eigenstate S=e^(i Q), with Q an allowed N-cell momentum.

[proved] On its b-dimensional induced offset space, R is a cyclic shift with closing phase e^(i Q). Its characteristic equation is z^b=e^(i Q). Therefore

    theta_r = (Q+2pi r)/b,       r=0,...,b-1.                  (7)

This can also be checked by forming the Fourier sum of the b offset states with coefficients exp(-i theta_r a). This is an ordinary tensor permutation. The common payload Jordan-Wigner seam is already part of the definition of S and Q, as in C27. It must not be inserted a second time into (7).

[proved] Combining r with the integer charge boost gives the charge-current lattice

    D_c = [k+Q/(2pi)]/b,       k in Z,                        (8)

up to an overall momentum-sign convention. Expand (4), keeping D_c bounded. With E_c(0) denoting the zero-flux sea energy,

    E_c - E_c(0) = (2pi v_c/L) D_c^2 + O(L^-3).              (9)

[proved, separated reference] Equations (7)-(9) replace b degenerate copies by one finer winding lattice. For Q=0, the nearest nonzero winding costs 2pi v_c/(b^2 L)+O(L^-3). It is a conformal-scale excitation, not an exponentially split crystal vacuum. At Q=pi the two nearest windings are equal in this separated reference, illustrating why uniqueness cannot be concluded just from (7).

[proved, finite-J distinction] R alone need not commute with the full interaction (2), because its A_j coefficients depend on the positions. Q and the separate twist branches are therefore labels of the separated reference, not exact conserved quantum numbers of the coupled dynamics. What survives exactly is the cyclic boundary holonomy in (1) and the physical Hilbert-space identification. Charge-payload mixing must respect that identification. The Gaussian continuation of its winding lattice is part of the matching assumption.

[proved, finite toy] For a fixed repeating internal color word and no other internal dynamics, the hopping graph on positions and offsets is connected when 0<M<L. Interior hops connect all position configurations at a fixed offset. A configuration containing L-1 but not 0 admits a cut hop, changing the offset by one, after which interior hops can return to the original position configuration. Repetition connects every offset. All edges have negative amplitude. Perron-Frobenius then proves a unique ground state for this finite toy. It does not prove uniqueness after adding the original complex h.

[heuristic, local order] Let phi_c count occupied particles in the continuum,

    y(x)=rho x-phi_c(x)/pi,
    H_c=(v_c/2pi) integral dx [(partial phi_c)^2+(partial theta_c)^2]

at the free reference. The marker spacing is b occupied labels. Its first physical harmonic is proportional to

    exp[2pi i y(x)/b]
       = exp[2pi i rho x/b] exp[-2i phi_c(x)/b].               (10)

[proved conditional on this Gaussian dictionary] This harmonic has dimension 1/b^2 and a two-point decay |x|^(-2/b^2). Its expectation vanishes in the translation-invariant liquid. For b=24 the exponent is 1/288, so finite systems can look almost ordered. The squeezed string order survives, while ordinary physical positional order has algebraic correlations. The offset is not an additional independent critical field or a residual local b-fold choice of vacuum.

[proved conditional] Below the defect threshold the allowed occupied-number changes are multiples of b. Rescale phi_marker=phi_c/b and theta_marker=b theta_c to preserve the canonical commutator. The marker liquid has K_marker=K_c/b^2, with K_c=1 at the reference. Its number quantum is one marker, or b occupied sites. Its integer windings reproduce (8), including the payload-dependent shift. Charge conservation forbids theta_c vertices in the Hamiltonian, even though charge-changing operators occur in the spectrum. Exact marker conservation introduces no second charge boson because M=bN in this low-energy sector.

## 4. Occupied-number adjustment and the chemical potential

[proved, reference thermodynamics] Let e_s be the bulk ground energy per occupied site of the selected H_s, with all payload densities optimized in that h. Multiplying the internal energy by J alpha gives the first-order bulk energy density

    e(rho) = -2t sin(pi rho)/pi + J e_s F(rho),
    F(rho) = rho^2-sin^2(pi rho)/pi^2.                        (11)

[proved, differentiation of (11)] The leading chemical potential and relaxed occupied-density curvature are

    mu = -2t cos(pi rho)
             + J e_s [2rho-sin(2pi rho)/pi],
    e''(rho) = 2pi t sin(pi rho)+4J e_s sin^2(pi rho).         (12)

[proved] If e_s<0, J|e_s|<pi t/2 makes this scalar reference curvature positive throughout 0<rho<1. This is only a sufficient reference condition. The full dressed density/current matrix, including payload susceptibilities, must also be positive. At rho=1/2, mu=J e_s and e''=2pi t+4J e_s to this order. Taking mu=0 and small J/t instead gives rho=1/2-J e_s/(2pi t)+o(J/t), still in the interior window.

[assumed, finite-J thermodynamics] In the homogeneous phase the true relaxed energy is C^2 in the density window, with (11)-(12) as its first-order asymptotics and small derivative corrections. Let rho_* solve e'(rho_*)=mu for this true energy, and let C_rho bound its positive second derivative nearby. This assumption excludes phase separation. The exact condition on mu is e'(2/5)<mu<e'(3/5), with a strict interior reserve. Equation (12) supplies the leading endpoints, not an exact finite-J equation of state.

[proved conditional on this smooth thermodynamics] Choose M_b=b round(rho_* L/b). Then |M_b-rho_* L|<=b/2. Taylor's formula for e(rho)-mu rho gives

    L[e(M_b/L)-mu M_b/L-e(rho_*)+mu rho_*]
        <= C_rho b^2/(8L).                                 (13)

If a third derivative is bounded, the leading term is e''(rho_*)(M_b-rho_*L)^2/(2L), with O(b^3/L^2) error. The bound follows also from the integral remainder with only a bounded second derivative. For fixed b and interior rho_*, M_b is admissible on every sufficiently large physical ring, including prime L.

[proved, inherited comparison at first order] In the unscaled squeezed chain a wrong cell costs a positive amount above the bulk baseline. The route_a proof gives gamma_b-3B per wrong cell before the O(B) finite-ring comparison, with B=sqrt(3)|lambda|, and gamma_b-5B>0 on compatible rings. Illegal cells have a similar positive margin. Applying (5) multiplies these margins by J alpha_(L,M). The first-order defect cost is therefore of order J, not J/L.

[assumed, needed comparison for the lift] At fixed sufficiently small J/t, the dressed defect sectors retain a positive cost Delta_def >= J alpha_min Delta_s/2 relative to the optimized homogeneous bulk. All order-one defect boundary self-energies must be included in this positive Delta_def; they are not discarded as finite-size errors. This is a local sector-stability requirement. The frozen-sector proof by itself is not a proof of this statement for dynamically occupied bonds. A vacancy can reduce the exposure of an unfavorable internal bond, and its energetic response must be included in matching.

[proved conditional] If homogeneous finite-size charge, payload and flux corrections are bounded by C_shell/L, comparing a compatible trial state with an incompatible sector gives compatibility for all

    L > [C_rho b^2/8 + 2 C_shell]/Delta_def,                 (14)

with the constant enlarged to cover the remainder bounds. The shell term includes the finite number of payload modes, common parity, and (8). A uniform O(1/L) shell bound is part of the assumed Gaussian spectrum, not a consequence of a bulk thermodynamic limit alone.

[heuristic] Thus the crossover length scales at least as t b^2/(J alpha Delta_s), potentially much larger after the payload hierarchy is included. Below it a wrong squeezed length can win. The grand-canonical Hamiltonian minimizes over M, but it does not dynamically change its exactly conserved M. If one externally fixes an incompatible canonical M, its defect does not disappear. At a fixed mu there can also be crossings between neighboring compatible M sectors, and arbitrarily small shell gaps along subsequences. Removing the geometric obstruction does not remove these usual Luttinger-liquid shell effects.

## 5. The full neutral vertex lattice

[assumed, continuum locality input] After eliminating the massive internal excitations, local neutral operators are built from the payload's even CAR algebra on the N-cell ring, the charge displacement, and the cell comb. Defect-creating operators do not become additional gapless primaries. This is the local projection and bosonization assumption; exact conservation of occupied and marker numbers and the selected cyclic marker sequence support it. The argument below includes every harmonic in this algebra, rather than retaining only the ordinary integer harmonics of the uncolored charge gas.

[proved, integer algebra] Label payload annihilations minus creations by z_i in Z, i=1,...,2n, with

    sum_i z_i=0,       q_z=sum_i z_i p_i,
    d_z=sum_i z_i^2/2.                                      (15)

The equality is exact total payload U(1), stronger than C20's parity condition. The numbers d_z are integers because z_i^2=z_i modulo two. At the free reference d_z is the primary dimension. Repeated identical branches require point splitting; the squared-charge formula already counts its Pauli suppression.

[heuristic, derivation of the physical dressing] The continuous cell coordinate is u=y/b. An operator with payload cell phase exp(i q_z u), sampled at integer cells, has reciprocal branches exp[i(q_z+2pi m)u], m in Z, by the integer-comb Fourier identity. Substitute y=rho x-phi_c/pi. Its primary, momentum and free dimension are

    w_(m,z) = [m+q_z/(2pi)]/b,
    V_(m,z) = exp[-2i w_(m,z) phi_c] V_z,
    Q_(m,z) = 2pi rho w_(m,z) modulo 2pi,
    D_0(m,z) = w_(m,z)^2+d_z.                              (16)

[proved conditional] Replacing a representative q_z by q_z+2pi shifts m and leaves the family unchanged. The b choices of m modulo b are the fractional harmonics associated with the b squeezed offsets. A local operator attached to another position inside the cell changes its form factor, not this lattice. At z=0,m=1, (16) reproduces (10). At b=1 it gives the familiar momentum dressing of a squeezed spin operator by the particle displacement.

[proved conditional] A uniform Hamiltonian vertex must satisfy

    rho w_(m,z)=ell,       ell in Z.                         (17)

At resonance its free dimension becomes

    D_res=ell^2/rho^2+d_z.                                  (18)

This is the crucial distinction between an observable with a small dimension and a uniform relevant perturbation. The dimension 1/b^2 of the first marker harmonic cannot be used as the dimension of an allowed lattice umklapp without imposing (17).

[proved conditional, classification] The complete classes are as follows. All momenta and dimensions in the table are those of (16), so every entry includes every m, not just m=0.

| Payload part | d_z | Physical momentum | Free dimension |
| :-- | --: | :-- | :-- |
| Identity, z=0 | 0 | 2pi rho m/b | (m/b)^2 |
| Neutral bilinear, z=e_i-e_j, i!=j | 1 | rho(q_z+2pi m)/b | w^2+1 |
| Nontrivial two-in/two-out, four distinct branches | 2 | rho(q_z+2pi m)/b | w^2+2 |
| All remaining nonzero neutral primaries | >=3 | rho(q_z+2pi m)/b | w^2+d_z |
| Current times a listed primary | add 1 | same | add 1 |
| Zero-primary currents and current products | 1 or 2 | 0 | 1 or 2 |

[proved] Exhaustiveness through d_z=2 follows by enumerating integer squared norms. Norm two requires +1 and -1 on distinct branches. Norm four with zero sum requires two +1 and two -1 on four distinct branches. Coincident incoming/outgoing branches cancel to a bilinear or to the zero-primary sector. A repeated incoming branch with two distinct outgoing ones has norm six, hence dimension three. Exact U(1) excludes anomalous pairing at every order. Current products in the zero-primary sector are Gaussian, including density-current and charge-payload products. Linear currents change backgrounds. Derivative descendants of the nonzero vertices have the same mismatch and higher dimension.

[proved conditional] With ell=0, (17) requires q_z=0 modulo 2pi and w=0. The inherited delta_2 and delta_4 exclude the nonzero d_z=1 and d_z=2 cases. Equal incoming/outgoing multisets are the allowed density interactions, not an omitted non-Gaussian resonance. With ell!=0, (18) is bounded below by 1/rho^2, regardless of b, flavor number or payload momenta.

### A sufficient open window, including all harmonics

[assumed, quantitative Gaussian condition] Normalize the retained fields so that the separated charge and free payload have identity scaling metric. Require the dressed scaling metric G to obey G >= (4/5)I, and the full energy stiffness H to be strictly positive. These are different requirements. A sufficient relative stiffness bound inherited from C40 is

    eta=||H_ref^(-1/2)(H-H_ref)H_ref^(-1/2)|| < 9/41,

because sqrt[(1-eta)/(1+eta)]>4/5. The reference must include all actual speed ratios and the correct compactification. Continuum dimensions are one quarter of the vertex vector contracted with G, not functions of an isolated charge K alone.

[proved conditional] On 2/5<rho<3/5, every ell!=0 resonant primary has

    D >= (4/5)/rho^2 > 20/9 > 2.                            (19)

Every primary with d_z>=3 has D >=12/5>2. The remaining zero-ell nonzero primaries are excluded by the payload margins. Thus all uniform non-Gaussian vertices are irrelevant. The surviving marginal terms are quadratic currents and are Gaussian. No extra irrationality or mixed-density exclusion is needed in this particular window. It includes half filling.

[proved conditional, finite detuning list] A candidate that could have D<=2 must have D_0<=5/2. Hence d_z<=2 and |w|<=sqrt(5/2). There are finitely many z and finitely many m for each z, even for b=24. Their unreduced physical momenta obey |Q|<2pi(3/5)sqrt(5/2)<2pi. Excluding the zero vertex, a uniform lower bound on their reciprocal mismatch is

    delta_phys >= min{ (2/5) min(delta_2,delta_4,2pi)/b,
                       2pi[1-(3/5)sqrt(5/2)] }.             (20)

For the C33 n=23 margins and b=24, the first entry exceeds 5.2596e-9 in physical inverse-site units. This is a deduction from its reported lower bound, not a new root computation. Dressed momenta must retain a reserved fraction of that bound. Detuned vertices can renormalize couplings before their oscillation cuts off; the amplitude condition in Section 6 still applies.

[proved conditional, checks away from the safe window] For a pure marker harmonic, rho m/b=ell. At rho=1/2, the first solution is m=2b, ell=1, with free dimension four. At rho=b/(b+1), m=b+1 is uniform and has free dimension (1+1/b)^2, less than two for b>=3. Such a term is permitted and can lock the fluid if present. Its coefficient and sign are not determined here. For arbitrary rational rho=p/q in lowest terms, the first pure harmonic is m=bq/gcd(p,b), ell=p/gcd(p,b), with dimension [q/gcd(p,b)]^2. These checks show why arbitrary vacancy density is not protected, and why large b alone is neither a proof of instability nor a proof of safety.

## 6. Finite J/t and conditions on the actual Hamiltonian

[proved, free charge, late source cross-check] The campaign's new `infrared/CHARGE_THRESHOLD.md` and `doping_check/CHARGE_THRESHOLD.md` strengthen the matching argument at half filling. Here is the complete short proof. Write every orbital momentum as k=s*pi/2+d relative to a closest Fermi edge, with |d|<=pi/2. Then |epsilon(k)|=2t*sin|d| >=(4t/pi)|d|. A number-preserving particle-hole excitation has equal numbers of particles and holes, so its signed Fermi-edge contributions sum to an integer multiple of pi. If P is its momentum, triangle inequality gives

    E_exc >= (4t/pi) sum |d| >= (4t/pi) dist(P,pi Z).        (21a)

This applies to any occupation set, with any number of particle-hole pairs, on an exactly half-filled free ring. Endpoint zero modes cost zero and do not affect the inequality.

[proved, free charge] Changing the internal twist from theta to theta' shifts the charge sea's physical momentum by (theta'-theta)/2 modulo pi at M=L/2. To conserve the full physical momentum, the target charge state must compensate that shift. The range of the sea energies in (4) is 2t*tan(pi/(2L)). Combining this fact with (21a) yields

    E_target-E_initial >= (2t/pi) dist(theta'-theta,2pi Z)
                           -2t tan(pi/(2L)).               (21b)

For theta'-theta=(q_z+2pi m)/b, minimization over the offset branches gives

    E_target-E_initial >= [2t/(pi b)] dist(q_z,2pi Z)
                           -2t tan(pi/(2L)).               (21c)

Thus the low-order nonforward payload channels have a limiting charge threshold greater than 8.371e-9*t at the C33 n=23 point. The threshold comes from momentum restoration, not from the energy difference of unconstrained twisted seas. It does not gap the forward charge fluid.

[proved, inherited rounded-ring bound, algebra checked] The doping_check proof also covers M=L/2+a with bounded a. Put A=ceil(|a|)+1. Relative to the sea filling negative-energy orbitals, the desired M-particle ground sea differs on at most A edge orbitals whose total Fermi-edge distance D_S is at most 2pi*A^2/L. Applying the same triangle inequality to the target and reference occupation sets gives an additional energy error at most (4/pi+2)t*D_S. The sea-momentum change differs from its half-filled value by at most 2pi*|a|/L. Consequently (21b) is replaced by

    E_target-E_initial >= (2t/pi) dist(theta'-theta,2pi Z)
       -t[8|a|+(8+4pi)A^2]/L -2t tan(pi/(2L)).              (21d)

Every error vanishes uniformly for bounded a, including odd physical L and nearest-multiple-of-b rounding at rho_*=1/2. This is a finite free-fermion inequality, not an interacting matching theorem.

[proved, scope of this added bound] A global momentum threshold of this form is special to the half-filled reference. At irrational density, arbitrarily high current harmonics have dense momenta modulo 2pi. The open-window argument in Section 5 instead bounds their scaling dimensions and requires controlled local matching. It does not assert (21a) throughout that window. The collaborator's finite-momentum argument supplies extra support for choosing the reference density exactly 1/2, while leaving the strict Gaussian phase conditions open under small dressed changes.

[heuristic, mechanism] Fast charge motion averages the exposure of every internal bond, giving J alpha H_s. Its deviations couple charge strain and density harmonics to internal energy operators. At the free Gaussian reference the smooth density response is

    chi_c(q,i omega)=K_c v_c q^2/[pi(omega^2+v_c^2 q^2)]
                  <= K_c/(pi v_c).                          (21)

[proved conditional on Gaussian charge response] Completing the Gaussian integral shows that a coupling J times charge density times a retained internal operator induces a kernel bounded by a constant times J^2/v_c. Equation (21) keeps its frequency dependence. It does not replace the charge mode by a massive denominator or justify dropping it. A source proportional to the total internal energy has only scalar and stress-tensor components at the fixed point; smooth strain times a payload stress tensor has dimension three. Other allowed current couplings belong in H.

[heuristic] Nonzero momentum channels are cut off by (20), or are uniformly irrelevant by (19). At a fixed h and fixed b these give finite matching constants, possibly extremely large. Charge-induced corrections are expected to be o(J) in a local low-frequency matching norm as J/t tends to zero, often of order J^2/t with logarithms. They can then be made small relative to J times the slowest internal scale. This argument is stronger than a charge-sea projection but is still a Wilsonian argument. It does not prove an all-order microscopic remainder estimate.

[assumed, explicit matching input M] Keep the low-energy charge field, every payload mode, and the leading mixed couplings explicitly. Subtract those retained couplings, the scalar charge corrections and J alpha H_s before defining R. Assume symmetry-preserving local matching supplies this remaining induced correction and controlled density/filling derivatives, with

    ||R||_LE / J -> 0 as J/t -> 0,

at the chosen fixed h, rho and b. Require the analogous control of the quadratic symbol and of the defect excess energy. A concrete sufficient majorant is

    ||R||_LE <= C_h (J^2/t)[1+|log(J/t)|^p],                 (22)

with finite C_h and p at fixed h. Equation (22) is a sufficient assumption, not a calculated bound for this 24-letter model. Irreducible kernels must be used, with retained charge/payload propagation subtracted. A finite-order cluster bound or an extensive O(LJ^2) estimate does not prove it. This is the new transfer obligation in addition to P.

[proved conditional, normalization of retained mixed terms] The leading mixed couplings need not be o(J) in an unnormalized local coefficient norm. For a Gaussian cross coefficient J*c_h, with c_h fixed by h, the relative energy-matrix entry is of order

    J*c_h/sqrt(v_c * J alpha v_s)
        = O(sqrt(J/t)) * c_h/sqrt(alpha v_s),

where v_s is a fixed internal speed in consistently converted coordinates. It tends to zero at fixed h even if c_h/v_s is very large. Payload-only induced errors bounded by (22), divided by the slow payload stiffness J alpha v_s, also tend to zero. Pure charge corrections have relative size O(J/t). This is why retaining the mixed fields is consistent with a small normalized Gaussian perturbation, and why demanding that every raw mixed coefficient be o(J) would be the wrong matching criterion. Oscillatory mixed vertices are controlled at their finite mismatch scale, not by this Gaussian calculation alone.

[proved conditional, sufficient order of choices] First fix n and a payload h satisfying P, its finite certificates, its reserved interaction and momentum margins, and the global selector bounds. Then choose any t>0 as the external energy unit. Choose a strictly interior density in (2/5,3/5), and take one fixed J>0 sufficiently small that all the following hold.

1. J||h||/t is small and the defect bound in Section 4 holds. Small J/t alone is not the norm condition when h contains large penalties.
2. The full stiffness is positive, G >=(4/5)I, and the roots remain within their reserved quadratic and quartic margins. The correction ||R||_LE/J must be small compared with alpha times the smallest payload scale, not merely compared with an order-one coefficient of h.
3. Matching through every payload bandwidth remains perturbative. With E_det no greater than the smallest linearization energy and v_min delta_phys/32, put ell_run=log(E_top/E_det). The collective marginal norm x must obey A x ell_run<1/2, with its induced metric and filling shifts inside the reserves. Each potentially relevant detuned vertex of dimension D and bare dimensionless amplitude y must obey |y|(E_top/E_det)^max(0,2-D)<<1. These are sufficient weak-RG inequalities, with matching constants evaluated for the actual dispersion. Marginal logarithms and flavor sums cannot be dropped.
4. Set mu in the true interval e'(2/5)<mu<e'(3/5), away from its ends. The leading formula is (12). Taking mu=0 works for sufficiently small J/t if the stated smooth matching holds, since rho then tends to 1/2.

[proved conditional] Under P, M and conventional continuous weak-coupling RG, these are simultaneously satisfiable for each fixed certified n. The payload margins and speeds are positive constants after h is chosen. The extra correction divided by J tends to zero by M. Even a logarithmic increase of the charge-to-payload scale ratio can be dominated by a sufficiently small power of J/t in (22). All final inequalities can be strict and persist on an open neighborhood, with density still in the safe window. No coefficient is sent to zero with L.

[proved, limitation on requested exact conditions] There is no supplied numerical value of C_h, no exact equation of state at finite J, and no proved defect/matching remainder for h_24. Consequently no justified explicit numerical interval 0<J/t<j_* for the full n=23 bond can be reported. Equations (12), (17)-(22), the strict stiffness and root margins, and P/M state the reference formulas and sufficient conditions without disguising an unknown microscopic threshold as an exact result. An assertion of unconditional acceptance would need to discharge this matching obligation.

## 7. Towers, uniqueness, and acceptance

[heuristic, Gaussian spectrum] In the stable phase there are n+1 right and n+1 left oscillator families. At leading separation the payload physical speeds are

    u_(a,s) = (b/rho) J alpha(rho) v_(a,s)^cell,

and the charge speed is v_c. The factor b/rho is the physical length per cell; it follows by converting a cell wave number into a physical wave number. Finite J dresses these values and mixes the Gaussian fields. The resulting spectrum has the form

    E-mu M = L omega_infinity
       - pi/(12L) sum_(a=0)^n [u_(a,R)+u_(a,L)]
       + E_zero/L
       + (2pi/L) sum_(a,s) u_(a,s) N_(a,s) + o(L^-1).        (23)

[assumed, content of E_zero] The positive quadratic zero-mode form includes M=bN, the offset M-rho_*L, the charge winding (8), the payload cell momentum Q, and its common Jordan-Wigner parity cosets. Exact only total payload U(1) is required; the individual flavor charges are emergent low-energy labels when conversion vertices flow away. Equation (23) is the usual Gaussian infrared prediction, not an exact finite-L formula for h_24.

[proved conditional] The b flux branches add allowed zero modes, not another factor of b in the oscillator central charge. Each compact boson contributes c=1, so the mode central-charge sum is n+1. The Casimir coefficient weights each contribution by its velocity. It cannot by itself measure the unweighted total c. At fixed n=23 the proposed value is 24, exceeding d-1=23 in the convention that allows the sum of conformal sectors with unequal velocities.

[proved] Exact finite-ring uniqueness needs more. Chemical-potential crossings between different exact M sectors occur at mu=[E_L(M')-E_L(M)]/(M'-M). Their union over finite L and charges is countable, so avoiding them removes those ties, provided each winning fixed-M sector is already nondegenerate. This does not prove that last premise for the complex internal Hamiltonian. Momentum-sector ties are another issue, as the Q=pi example shows. No symmetry or positivity argument in the supplied data excludes every such tie for all L. Generic uniqueness is a reasonable heuristic, not a finite-dimensional theorem established here, and its all-L parameter set need not be open.

[proved conditional, shell obstruction] Even when uniqueness holds at each L, a density with irrational rho_*/b gives a nearest-integer offset N-rho_*L/b whose fractional part does not converge. In a quadratic charge spectrum, approaching a half-integer makes the gap between neighboring number sectors o(1/L) along a subsequence. Thus a fixed untwisted identity-sector spectrum on all integer L is already too strong for the reference liquid. Specifying the affine sectors or a subsequence is necessary. This is separate from the unwanted finite-energy wrong-cell defect, which (14) removes in the grand-canonical ground state.

[heuristic, acceptance at fixed n] A careful physicist can accept the n+1-mode conclusion as the standard weak-coupling inference under P and M. The period-b objection is resolved at that level, including its fractional local harmonics and its flux sectors. The d=24 proposal is substantially more concrete than the earlier d=2736 mobile wrapper because it retains the certified original h. It is not yet an established unconditional c=24 microscopic counterexample with proved all-ring uniqueness. The remaining dynamical step is finite-J transfer with the actual hierarchy. The new free-charge momentum-restoration bound supports that transfer, but does not prove the required local remainder or dressed defect stability.

[proved, all-n scope] The same kinematics and density-window argument work for every fixed finite b. What remains for all n is a simple-root, nonresonant reference for the actual original h at every n, then compatible payload and vacancy matching hierarchies for each member. C33's all-n positivity/crossing theorem does not prove the needed pair-sum nonresonance. C37 stops at n=23. No uniform-in-n neighborhood, velocity bound, onset length or finite-c accumulation theory follows. If these fixed-n obligations were solved for all n, this same d=24 alphabet would yield unbounded mode central charge at fixed d under the stated infrared convention. That is a conditional implication, not a result obtained here.

## 8. Independent checks and task accounting

[numerical, finite matrices] `checks.py` independently constructs hopping matrices on positions and b offsets, then compares their complete spectra with the union of twisted charge spectra. Cases (L,M,b)=(5,2,2),(7,4,2),(7,3,3),(8,4,4),(9,6,3) have dimensions 20,70,105,280,252. The largest error is 1.51e-14. All have one ground state, in agreement with the finite-toy proof above. These tests contain no n=23 payload.

[numerical, finite projection] Independent random complex Hermitian two-label bonds were used at (L,M)=(5,3),(6,4),(6,3), full dimensions 80,240,160. The complete sea projection, including both half-flux seas at M=4, agrees with (5) to at most 3.56e-15. No neighboring implementation was imported. These are checks of the proved basis identities and projection, not central-charge estimators.

[numerical, asymptotics] At b=24 and payload cell momentum Q=0.71, the exact free charge flux energy agrees with (9) at L=96,192,384,768, with relative errors 1.78e-4,4.44e-5,1.11e-5,2.78e-6. The factor-of-four decrease is the expected relative O(L^-2) error. A separate number-rounding toy uses a fixed artificial incompatible-sector penalty 0.01. It selects incompatible M at L=4000 and 32000, but compatible M at L=8000,16000,64000. This deliberately illustrates that compatibility need not improve monotonically before the uniform threshold. It is not a simulation of the microscopic defect gap.

[proved, exact arithmetic check] Standard-library Fraction arithmetic verifies (4/5)*3=12/5>2 and (4/5)/(3/5)^2=20/9>2. The general proofs are in Section 5. No floating-point sign test is used to prove a finite-dimensional spectral or positivity assertion.

[proved, reproducibility and scope] Run `python3 -B checks.py` in this folder. It uses only Python, NumPy and SciPy and writes `checks.json` here. All calculations completed within seconds. `PROGRESS.md` and `STATE.md` contain the handoff. No network, subagents, cluster calculations, source modifications outside this folder, repository operations, or microscopic central-charge measurement were used. The project-level STATE file was left unchanged under the explicit folder restriction.
