# C15. Unitary on-site Z2 qubit chains

## Verdict and status

[Proved / assumed standard inputs / heuristic, separated below.] This run does not establish c <= 1 on the whole four-parameter stratum. It establishes the bound on the quadratic locus, on the supplied XXZ and zero-field XYZ loci, and in regular perturbative scaling regimes descending from their c <= 1 continuum theories. It also controls ordinary large-field and large-|V| corners, with qualifications at their degeneracy surfaces. No c > 1 candidate survives the tests here. This is not an exclusion of an isolated interacting multicritical point at intermediate coupling.

Two tempting general exclusion arguments fail.

1. [Proved counterexample to the premise.] K >= 1/2 is not a bound on all nearest-neighbor spinless-fermion liquids. Strong repulsion away from half filling gives K=(1-rho)^2 and hence K down to 1/4. Atom-pair conversion can therefore be irrelevant on scaling-dimension grounds. The obstruction is realizing the additional stable pair fluid in the same microscopic chain.
2. [Proved conditional operator analysis.] Symmetry does not force an Ising times Luttinger point to have codimension greater than four. Its operator count depends on how parity and translation act on the two sectors, and on the local operator algebra. An Ising boundary touching a liquid is also not evidence that their central charges add.

[Proved.] A singular large-repulsion corner deserves particular attention. With hopping zero, virtual pairing generates hopping by two sites and two minima. Its leading Hamiltonian is a ferromagnetic degeneracy problem, not two independent liquids. An exact line through this corner has L+1 ground states on even rings and a quadratic excitation branch. This disposes of the exact line as a conformal counterexample, but does not classify every perturbation around it.

All work is local to this folder. No network, external-source download, cluster calculation, DMRG, or edit to the parent project was performed. Literature statements below are inherited from the supplied background, C6, and `contrib/d2_classification/README.md`, not independently verified against the papers. This is an unrefereed calculation.

## 1. Conventions and the free scaling regimes

[Proved by the on-site rotation and Jordan-Wigner map.] Use

    H = sum_j [Jx X_j X_{j+1} + Jy Y_j Y_{j+1} + Jz Z_j Z_{j+1}
               + D (X_j Y_{j+1} - Y_j X_{j+1}) - h Z_j].

The symmetric XY cross term has been removed by a uniform rotation about z. Define

    a = Jx+Jy, b = Jx-Jy, V=4Jz,
    tau = -a+2iD, delta=b, mu=V-2h, mu_bar=-2h.

With n=(1-Z)/2, the fermion Hamiltonian, up to a constant, is

    H = sum_j [-tau c_j^dagger c_{j+1} - tau* c_{j+1}^dagger c_j
          + delta(c_j^dagger c_{j+1}^dagger + c_{j+1} c_j)
          + V n_j n_{j+1} - mu n_j].                         (1)

Equivalently the density terms are V(n_j-1/2)(n_{j+1}-1/2)-mu_bar(n_j-1/2). The sign of real delta is removable by a uniform rotation. There are five real coefficients and four after scale. The spin form, including periodic spin boundary conditions, is what all full-chain ED scripts implement. Jordan-Wigner boundary conditions depend on parity.

### 1.1 Exact free regions

[Proved.] At V=0 the BdG eigenvalues, up to k reversal, are

    lambda_+/- = 4D sin k +/- 2 sqrt[(a cos k+h)^2+b^2 sin^2 k].    (2)

The supplied free-fermion theorem gives c <= 1 whenever the infrared is conformal. This holds even at coincident zeros if their conformal description exists. It does not assign a central charge to a nonlinear or flat branch.

[Proved by solving (2).] Put s^2=4D^2-b^2. For s^2>0, zeros obey

    (a x+h)^2 = s^2(1-x^2), x=cos k,
    x = [-a h +/- s sqrt(a^2+s^2-h^2)]/(a^2+s^2).               (3)

The regular chiral liquid occupies |h|<sqrt(a^2+s^2), and has c=1. It has nonzero linear slopes in the interior, including h=+/-a, where the representation of some zeros changes. At |h|=sqrt(a^2+s^2), a pocket closes with a quadratic edge. Outside, the ground state is gapped.

[Proved.] For s^2<0, the only zeros are at sin k=0 and h=+/-a. If a is nonzero, each surface is a single Ising theory, c=1/2. If a=h=0, both k=0 and k=pi are massless and their total c is one. The velocities are nonzero because |b|>2|D|. The remaining points are gapped.

[Proved.] The special case b=D=0 is the usual number-conserving cosine band. For a nonzero it has c=1 at |h|<|a| and quadratic saturation at |h|=|a|. If s^2=0 but b is nonzero, |h|<|a| has a quadratic touching branch, and h=+/-a generically has a cubic branch. The case a=h=0 on this threshold has a flat zero branch. These are not regular relativistic starting points. If |h|>|a| the system is gapped. The identically zero Hamiltonian is excluded from any CFT discussion.

[Check.] At the Ising/chiral collision a=-1, b=-1/2, D=1/4, h=1, the small positive momentum branch is

    e_slow(q) = sqrt[4(1-cos q)^2+sin^2 q] - sin q
              = q^3/2 + O(q^5).                              (4)

A second linear branch remains, but this is not Ising plus an independent finite-velocity liquid. The free rank bound and the cubic expansion agree on that point.

### 1.2 What the weak-coupling argument actually covers

[Assumed standard RG, with an explicit domain of applicability.] Fix a regular free point and choose a continuum cutoff below every discarded band or momentum scale. Retain all its linear modes. Let v_min>0 be their smallest velocity magnitude and let E_sep>0 be the energy below which no additional low-energy region enters. A sufficient perturbative regime has dimensionless interactions |V|/(2 pi v_min) small and interaction-induced shifts small compared with E_sep. At critical surfaces one retunes the allowed mass. For each compact regular piece, this gives a local scaling regime. There is no uniform neighborhood extending to v_min=0 or E_sep=0.

[Assumed standard RG.] The regular liquid has one compact boson. Forward scattering changes its stiffness and velocities. Momentum-allowed pairing and umklapp are sine-Gordon perturbations; they can remain irrelevant, become marginal, or gap degrees of freedom. The regular Ising point has one Majorana; after its mass is tuned, the weak residual interactions cannot supply a new independent low-energy band. At the two-Majorana free point, the marginal four-fermion interaction is the usual Ashkin-Teller/Gaussian interaction, with velocity anisotropies retained if present. This is a continuum theory with at most the original c=1 degrees of freedom, within its weak scaling regime.

[Proved conditional on a relativistic RG trajectory.] If this continuum theory provides a unitary Lorentz-invariant ultraviolet fixed point for the subsequent flow, the c-theorem gives

    c_IR <= c_continuum_UV <= 1.                              (5)

A common drift of left/right velocities can be removed by a moving frame. The c-theorem alone does not apply to arbitrary unequal species velocities without further input. In those cases the weak one-boson or two-Majorana analysis just given is the standard physics input; one should not claim a general Lorentz-invariant theorem for the lattice.

[Scope.] Thus "perturbatively connected" must mean connected by this specified continuum RG construction, including its relevant deformations. It cannot mean reachable along an arbitrary path of lattice Hamiltonians. Gapped free points have stable gapped neighborhoods. Quadratic saturation, cubic collisions, flat bands, ferromagnetic degeneracies, and flows requiring additional soft modes fall outside the weak-coupling extension. The same regular-scaling qualification applies to perturbing known interacting c=1 or c=1/2 lines.

## 2. Strong coupling

### 2.1 Classical large-density-interaction limits

[Proved at zero kinetic and pairing couplings.] For V>0, the uncentered chemical potential in (1) gives empty sites for mu<0, the two alternating states for 0<mu<2V, and a filled chain for mu>2V. At mu=0 every configuration without adjacent occupied sites is degenerate. At mu=2V there is the particle-hole-related constrained-hole degeneracy.

[Assumed standard gapped perturbation theory.] Compact interiors of these classical phases remain gapped when |tau|,|delta| are small compared with their local excitation gaps. The two CDW ground states are a broken discrete symmetry, not an extra CFT. A large field with a finite gap to all spin flips is likewise a stable polarized phase. These claims do not include field values tuned to the classical crossing surfaces.

[Proved at the classical level; controlled perturbative conclusion.] For V<0 the empty and full states cross at mu=V. A single finite occupied domain of N sites in the empty background costs

    E_domain = V(N-1)-mu N = |V| at mu=V.                     (6)

There are two domain walls and no energy gain from forming a gas of separated adjacent dimers. For |V| much larger than both kinetic couplings, the wall gap survives small perturbations. The crossing is first order, with gapped bulk excitations within either phase. A cat state of the two vacua contributes at most an O(1) entropy and does not add a conformal sector.

### 2.2 Regular repulsive boundary: one constrained liquid

[Proved leading projection.] Near mu=0 with V large and positive, let P forbid adjacent occupied sites. Then

    P H P = -sum_j P(tau c_j^dagger c_{j+1}+h.c.)P - mu N.     (7)

Nearest-neighbor pair creation has zero matrix element within P. Corrections arise at orders |tau|^2/V, |tau delta|/V and |delta|^2/V, with anomalous processes possible at higher orders as well. In this subsection ordinary hopping is nonzero and dominates those corrections. Its phase is gauged to a boundary twist.

[Proved in the projected model.] Ordered particle coordinates x_i with x_{i+1}>=x_i+2 become free hard-core coordinates y_i=x_i-(i-1) on a squeezed length L-N, up to the ring's global boundary condition. Put f=1-rho and k_F'=pi rho/f. The thermodynamic energy and sound velocity are

    e(rho) = -(2|tau| f/pi) sin k_F',
    v(rho) = 2|tau| sin k_F'/f.                              (8)

Differentiating, using dk_F'/d rho=pi/f^2, gives

    e''(rho) = 2 pi |tau| sin k_F'/f^3.

The Gaussian compressibility relation e''=pi v/K then gives

    K=(1-rho)^2, 0<rho<1/2.                                 (9)

There is one squeezed Fermi sea and one boson, c=1. The endpoints have zero velocity and are not covered by the interior CFT.

[Proved consequence.] For rho>1-1/sqrt(2), K<1/2, and K approaches 1/4 at the CDW edge. The commonly used K>=1/2 statement refers to the gapless half-filled XXZ line with anisotropy at most one, not to doping an Ising-like CDW. Equation (9) alone supplies a nearest-neighbor counterexample to that broader statement.

[Assumed bosonization and controlled perturbation theory.] At an incommensurate density, uniform pairing gives cos(2 theta) with dimension 1/K. It is irrelevant when K<1/2. A p-particle commensurability permits cos(2p phi) with dimension p^2 K; it must also be irrelevant or absent for the liquid to survive. For example at rho=1/3, K=4/9, pairing has dimension 9/4 and the first density umklapp has dimension 4. Both are irrelevant. Small finite-V corrections consequently give a controlled one-boson floating regime, even with only microscopic Z2. If pairing is relevant, the one-boson theory loses a mode rather than acquiring a second one.

[Scope.] Keep rho in a compact interior interval and the generated couplings small at the continuum cutoff. This establishes c<=1 there. It does not cover |tau| of order |delta|^2/V or smaller, nor the zero-velocity endpoints. Adjacent pairs in this repulsive regime cost order V and do not constitute a second low-energy fluid.

### 2.3 Singular repulsive boundary: vanishing ordinary hopping

[Proved by second-order degenerate perturbation theory.] Set tau=0, delta real, mu=O(delta^2/V), and V much larger than |delta|. In units delta^2/V, the effective Hamiltonian in P is

    diagonal(n) = -sum_{j: n_j=n_{j+1}=0} 1/(1+n_{j-1}+n_{j+2}) - m N,
    hopping by two sites = -1/2 when the final configuration is in P,
    m = mu V/delta^2.                                       (10)

To obtain (10), a pair can be created only on a 00 bond. Its interaction cost is V(1+n_{j-1}+n_{j+2}). Removing the same pair gives the diagonal term. If exactly one exterior neighbor was occupied, removing the other adjacent pair moves that particle by two sites, through a state of energy 2V. If both exterior neighbors were occupied, the shifted configuration violates the constraint. This accounts for every second-order process.

[Proved.] Let d_i>=2 be consecutive cyclic particle separations and N_r their multiplicities at d_i=r. For N>0 the diagonal term at m=0 is

    D(n)=-L+3N-N_2-N_3/3.

A gap with d>=4 allows two directed hops, whereas gaps 2 and 3 allow none. Thus the off-diagonal graph has weighted degree N-N_2-N_3, each edge having weight 1/2. Therefore

    H_eff(m=0) = -L+2N + graph_Laplacian + (2/3)N_3.          (11)

Both last terms are positive. If every particle lies on one sublattice, every gap is even and the uniform vector over that connected configuration graph has zero Laplacian energy. Hence E_min(N)=-L+2N for every 0<=N<=L/2 on even rings. At m=2 all these sectors are degenerate. The two one-particle minima are not evidence for two stable liquids; this leading theory is at a ferromagnetic degeneracy.

[Check.] On one occupied sublattice the same Hamiltonian is

    H_eff+L = 3N - (1/2)sum(b_l^dagger b_{l+1}+h.c.) - sum n_l n_{l+1}.

Its effective attraction is -2 times its hopping, the ferromagnetic t-V threshold. The full second-order operator was also checked against -P A Q H_0^{-1} Q A P at L=6, with zero matrix-element discrepancy.

[Proved nonconservation warning.] Total N is accidental at this order. The exact symmetry at tau=0 is staggered charge sum_j(-1)^j n_j. At third order, 0000 -> 1100 -> 1111 -> 1001 and its reversed creation order generate a source for a separated opposite-sublattice pair, with matrix element 2 delta^3/(3V^2) for an isolated four-site interval. The two virtual excitation energies are V and 3V. `operator_checks.py` checks this coefficient at L=8. Thus the second-order two-species particle-number description is not an exact two-U(1) model.

### 2.4 An exact nonconformal line through the singular corner

[Proved.] For tau=0 and 0<2|delta|<V, define

    mu_f = V-sqrt(V^2-4 delta^2).                            (12)

Take delta>0 by a rotation. A local bond, with the chemical potential equally divided between its sites, has the even-parity block

    [[0, delta], [delta, V-mu_f]],

and the odd states |01>,|10> both have energy -mu_f/2. The lower even eigenvalue is also -mu_f/2 because

    mu_f(2V-mu_f)=4 delta^2.

Consequently

    h_bond + mu_f/2 = V |v><v|,
    |v> = A |00> + B |11>,
    A=sqrt(mu_f/(2V)), B=sqrt(1-mu_f/(2V)).                 (13)

This is positive and rank one. On even rings, flip odd spins. The forbidden bras alternate between A<01|+B<10| and A<10|+B<01|. Apply the invertible filter diag(B,-A) on the even sites of a trial state. Each forbidden bra becomes proportional to <01|-<10|. The resulting constraints require symmetry under every adjacent swap. Their joint kernel is the fully symmetric spin subspace, dimension L+1. Inverting the filter proves that the original ground space also has dimension L+1. This proves the ground-space count, not a spectral equivalence under a nonunitary filter. At delta=0 the filter is singular: the forbidden vector becomes |11> and the ground space is the exponentially large hard-blockade space, not the L+1 space. That endpoint is excluded from this count.

[Proved exact spectral check.] Start from a Neel ground state. A state with one particle removed from its occupied sublattice costs mu_f. A state with one particle added to its empty sublattice costs 2V-mu_f. Pairing connects the neighboring defect states with Bloch matrix element delta(1+e^{ik}), where k is momentum per two-site cell. This sector is closed, giving

    E_+/- (k)=V +/- sqrt[(V-mu_f)^2+4 delta^2 cos^2(k/2)].    (14)

At (12),

    E_-(k)=V-sqrt[V^2-4 delta^2 sin^2(k/2)]
          = delta^2 k^2/(2V)+O(k^4).                       (15)

The extensive ground-space dimension and this quadratic branch exclude an ordinary finite-velocity CFT description. They are not a c=2 liquid. Expansion of (12) gives mu_f=2 delta^2/V+O(delta^4/V^3), in agreement with (11).

[Open.] The exact line is resolved, but its punctured neighborhood is not covered by a relativistic c-theorem. Hopping, chirality, higher-order virtual processes and field detuning compete there. The leading result does not justify an exhaustive classification of all nearby interacting critical surfaces.

### 2.5 Ferromagnetic, saturation, and XYZ-like limits

[Assumed standard XXZ solution, supplied by C6.] At delta=0, remove the phase of tau by a spiral rotation. The full t-V chain in a uniform chemical potential has a single massless Bethe sea whenever its ground state is conformal, hence c=1. This statement includes all longitudinal fields and both signs of V. A boundary twist changes zero-mode energies, not bulk c. Ferromagnetic V=-2|tau| at the appropriate coexistence field has a quadratic magnon and a macroscopic ground multiplet. Saturation has a quadratic dilute-particle threshold. Neither is itself a conformal exception to c<=1.

[Assumed standard dilute-gas analysis.] With one isolated nondegenerate quadratic minimum, a stable dilute gas gives one LL at nonzero density, and the vacuum threshold has z=2. Anisotropy can reduce this to an Ising transition or open a gap. This is a useful local classification when the isolated-species assumption is checked. It is not an exclusion theorem for multiple bound-state thresholds or singular ferromagnetic points. Equations (10)-(15) show why that qualification matters.

[Assumed standard XYZ solution, supplied by C6.] At D=h=0 the zero-field XYZ model is controlled. Its regular gapless lines reduce, with allowed axis or staggered rotations, to one Gaussian/XXZ theory, including the antiferromagnetic SU(2)_1 point. The ferromagnetic endpoints are nonconformal. Dominant single-axis exchange has gapped ordered phases. Small perturbations of a regular XYZ critical line are perturbations of one boson; small perturbations of a gapped region remain gapped. A generic finite field and D do not retain XYZ integrability.

[Inherited literature input, not independently reread.] The Chepiga-Mila phase diagram cited in the background contains floating c=1 phases, Ising boundaries, and an eight-vertex c=1 locus. The repository's account does not establish exhaustiveness for arbitrary chiral hopping phases. An Ising line ending on a Gaussian or saturation boundary does not imply a product CFT: it may describe the same degrees of freedom, or terminate at zero velocity.

## 3. Atom-pair conversion and the best two-fluid attempts

### 3.1 Scaling dimension and momentum selection

[Proved within the specified Gaussian theory.] For independent fluids use

    H_G = sum_a v_a/(2 pi) int dx [K_a (partial theta_a)^2
                                      + K_a^{-1} (partial phi_a)^2],
    [phi_a(x),partial_y theta_b(y)]=i pi delta_ab delta(x-y).

A vertex exp(i m theta_a) has dimension m^2/(4K_a), and exp(i n phi_a) has dimension n^2 K_a/4. The spinless atomic right-left pair is proportional to exp(2i theta_1); at finite density its adjacent-site p-wave amplitude is nonzero. A pair annihilation field has leading vertex exp(i theta_2). Conversion therefore contains

    g cos(2 theta_1-theta_2 + delta_q x + alpha),
    d_conv = 1/K_1 + 1/(4K_2).                              (16)

This conserves total charge N_1+2N_2 and is allowed even with exact U(1). It is relevant when d_conv<2 and nonoscillatory. Relevance pins one relative phase and leaves at most one charge mode, provided no other gapless sector is present. With only Z2, additional charge-nonconserving sources must also be checked.

[Proved qualification.] Interfluid forward scattering generally invalidates the sum of two independent dimensions. Define the phase covariance matrix M by

    <theta_a(x) theta_b(0)> = -(M_ab/2) log|x| + O(1).

Then d_conv=(1/4) l^T M l with l=(2,-1), hence M_11-M_12+M_22/4. In the independent case M=diag(1/K_1,1/K_2). Cross density/current couplings can change the relevance criterion; separate K values alone are then insufficient.

[Proved consequence.] Even in the independent case there is no universal relevance theorem. If K_2>1/8, conversion is irrelevant for

    K_1 < 4K_2/(8K_2-1).                                   (17)

For K_2=1 this threshold is 4/7. In particular K_1=1/2 already makes (16) greater than two for every finite positive K_2. If K_2<=1/8, (16) is greater than two for every finite positive K_1. The NN repulsive example (9) invalidates any attempted exclusion based on K_1 never falling below 1/2. Its adjacent pairs are, however, high-energy objects.

[Proved momentum rule.] If the atomic Fermi points are k_R,k_L and the pair-fluid minimum is P_2, then delta_q=k_R+k_L-P_2 modulo 2 pi. A nonzero mismatch makes the displayed leading vertex oscillatory. Commensurate harmonics or same-branch derivative vertices must still be considered. A pair band with two inequivalent minima can create further conversion and interminimum-scattering channels. One cannot conclude "conversion is allowed by parity, therefore it is nonoscillatory."

### 3.2 Exact two-body test with nearest-neighbor attraction

[Proved.] Set delta=0, write tau=t exp(i alpha) with t>0, and use an infinite chain. For two particles at separation r>=1 and total momentum P, the relative hopping is

    T_P=2t cos(P/2+alpha).

The equations are

    (E+2mu) f_r = -T_P(f_{r-1}+f_{r+1}), r>=2,
    (E+2mu-V) f_1 = -T_P f_2.

Putting f_r=z^{r-1} gives z=-T_P/V and

    E_pair(P)=V + (4t^2/V)cos^2(P/2+alpha)-2mu,
    |2t cos(P/2+alpha)|<|V|.                                (18)

For attractive V, the minimum is at P=-2alpha whenever that minimum is a normalizable bound state, which requires |V|>2t. There is a single minimum modulo 2 pi. Expanding the square gives a single cosine of P+2alpha, not independently tunable minima at +/-Q. The atomic band minimum is -alpha, and its two Fermi points in the U(1) liquid sum to -2alpha. Thus the hopping phase shifts both fluids together and does not make their leading conversion oscillate.

[Proved caveat.] For -2t<V<0 there are still bound states at momenta sufficiently far from -2alpha. These are not low-energy pair-fluid minima below the global two-atom bottom. On approaching their unbinding edge their energy is 2V>-4t, whereas two atoms can have energy arbitrarily close to -4t. Ignoring this momentum restriction would falsely infer a dilute paired ground phase for every attraction.

[Proved strong-attraction obstruction; assumed exact XXZ completeness at all densities.] For V<-2t, the minimum pair energy per particle is

    (V+4t^2/V)/2,

but a long occupied domain has energy per particle tending to V, which is smaller because (V^2-4t^2)/(2V)<0. The vacuum-to-dense transition occurs at mu=V before the single-atom threshold mu=-2t. Larger clusters preempt a dilute dimer gas. The full XXZ result rules out a separate atom-plus-pair ground-state liquid everywhere on delta=0, not just in the dilute limit.

[Scope.] Equation (18) is not an exact two-body solution at nonzero delta, where particle number is not conserved. It remains the controlled starting point for small delta compared with finite binding or wall gaps. At intermediate coupling, renormalized composite dispersions can acquire more harmonics, and the single-cosine result does not exclude two minima there.

### 3.3 What candidate remains

[Controlled exclusion.] The most direct candidate would combine an attractive adjacent-pair band with a stiffly repulsive atomic liquid so that (16) is irrelevant. The microscopic NN density interaction has one sign: its controlled attractive regime clusters, while its controlled repulsive regime penalizes adjacent pairs. There is no independently adjustable pair hopping, pair repulsion or atom-pair conversion amplitude. The exact U(1) solution and the strong-coupling analysis eliminate this construction in their respective domains.

[Heuristic candidate search, not an existence claim.] The best unresolved location is where these controlled descriptions break down: a finite-coupling composite threshold, or the neighborhood of (12), where virtual hopping creates two minima and degeneracies make higher-order terms compete. The leading theory and exact line there are ferromagnetic, so they are not a constructed two-fluid counterexample. An interacting two-fluid or Ising-plus-fluid point away from them would require an independent microscopic demonstration. The ED probes below provide no such demonstration.

## 4. Could a Majorana coexist with a liquid?

[Proved conceptual distinction.] Ising x Gaussian has c=3/2 if both sectors are independent in the infrared and both velocities are nonzero. The word "meets" in a phase diagram does not establish either condition. A one-boson theory can have Ising and Gaussian critical boundaries without ever having three massless Majoranas. The free example (4) instead loses a velocity at its apparent meeting point.

### 4.1 The relevant-operator count is conditional

[Proved within a chosen product operator algebra.] Write the Ising energy as epsilon=i xi_R xi_L, dimension one, and its order field as sigma, dimension 1/8. With the boson convention above, the following are possible operators. Every row requires locality and the stated momentum/parity permission.

| Operator | Dimension | Permission / consequence |
|---|---:|---|
| epsilon | 1 | Ising mass is generically allowed, requiring one tuning. |
| cos(2 theta), sin(2 theta) | 1/K | Atomic pair source, allowed by fermion parity when its momentum is zero. Relevant for K>1/2. |
| cos(2p phi), sin(2p phi) | p^2 K | Density locking only at the corresponding commensurability. Relevant for p^2 K<2. |
| epsilon times exp(2i theta) and conjugate | 1+1/K | Relevant for K>1 when local and nonoscillatory. |
| epsilon times exp(2ip phi) and conjugate | 1+p^2 K | Relevant for p^2 K<1 when momentum permits. |
| epsilon partial_x phi or epsilon partial_x theta | 2 | Marginal couplings; their RG flow and the velocities must be checked. |
| sigma | 1/8 | Allowed only if no protecting microscopic symmetry forbids the longitudinal field. |
| sigma exp(i theta) and conjugate | 1/8+1/(4K) | A possible mixed order operator if both factors are parity odd and their momenta cancel. Relevant for K>2/15. |
| exp(i theta_pair) and conjugate | 1/(4K_pair) | A source for a parity-even pair fluid, unless translation or a stronger symmetry forbids it. |

[Qualification.] This is a diagnostic list, not a universal local spectrum. In a fermionic embedding, a bare exp(i theta) need not be a local atomic operator. An Ising disorder field is not automatically local either. Projection and spin structures can remove or combine rows. Conversely, a local mixed fermion bilinear may be allowed. For a Gaussian fermion of dimension (K+1/K)/4 coupled to a Majorana, its dimension is 1/2+(K+1/K)/4, often relevant, unless their momenta differ. The physical operator embedding has to be specified before counting.

[Proved conditional codimension examples.] Suppose translation makes every nonderivative mixed vertex oscillatory, parity forbids sigma, and all boson locking vertices are irrelevant. Then the only relevant scalar is epsilon. The product point has codimension one in that continuum theory, subject to marginal stability. If one nonoscillatory real boson-locking amplitude is relevant as well, there are at least two independent tunings. If a complex mixed order amplitude is allowed, its real and imaginary components may supply two more relevant directions. Additional reflection or antiunitary symmetries can fix phases and reduce that count. Treating every complex amplitude as two independent microscopic directions without checking these relations would also be wrong.

[Conclusion.] Four available microscopic parameters do not exclude these counts. Conversely, a continuum count <=4 does not prove that the microscopic family intersects the corresponding stable manifold. One must exhibit a point and compute the coupling map, its rank, the momentum assignments and the marginal flow. None is supplied by Z2 or by having one orbital per site. There is therefore no justified universal codimension for "Ising x LL somewhere in this family," and no established reachability result in this run. This is a located obstruction to task 3, not a proof of absence.

[Controlled exclusion.] A product with c=3/2 cannot be the infrared of the regular weak continuum flows covered by (5). It does not occur on the exact XXZ/XYZ loci under the supplied exhaustive inputs, on the free locus, or on the nonconformal lines (4), (12). What remains is an interacting finite-velocity point outside those controlled basins, including possibly a perturbation of a singular point whose scaling theory changes.

## 5. Numerical tests and their limits

[Numerical.] `checks.py` uses sparse Lanczos in both exact spin-parity sectors. These sectors contain all momenta, so the global ground-state search does not assume a momentum. No translation-sector tower assignment is made. Entropies use the normalized ground eigenvector in its selected parity sector. Parameters below use the spin Hamiltonian in section 1, with field -hZ.

The fixed-size estimator fits

    S_L(l)=(c_chord/3) log[(L/pi) sin(pi l/L)]+s0

over max(2,floor(L/4))<=l<=L/2. The size estimator fits

    S_L(L/2)=s0+(c_half/3)log L+b/L^2.

The independent energy calibration fits

    E0(L)/L=e_infinity+A/L^2+B/L^4,
    c_E=-6A/(pi v).

The controls have known v=2. For a candidate with unknown or multiple velocities, -6A/pi is only a velocity-weighted coefficient and is not reported as c. No Casimir interpretation is assigned to a nonlinear branch.

### 5.1 Calibration and direct probes

[Numerical, L=10,12,14,16.] The dense Kronecker Hamiltonian and bit implementation agree at L=4 to 1.7e-16. Largest ground-state residual in `checks.py` is 1.9e-13; including the singular-corner probe it is 1.5e-12. Calibration results are

| Control | c_chord at 16 | c_half | c_E |
|---|---:|---:|---:|
| Ising, (-1,0,0,0,1) | 0.500721 | 0.500006 | 0.499995 |
| XX, (-.5,-.5,0,0,0) | 1.007044 | 1.000219 | 0.999839 |

These are two physically different estimators, entanglement and Casimir energy. Four-size fits are checks, not statistical error estimates.

[Numerical.] Three probes at the same sizes give

| Probe, (Jx,Jy,Jz,D,h) | c_chord(10) | c_chord(16) | c_half | last size-step c |
|---|---:|---:|---:|---:|
| Free collision, (-.75,-.25,0,.25,1) | .509028 | .502830 | .498265 | .510923 |
| Untuned interacting collision, (-.75,-.25,.1,.25,1.2) | .496929 | .491005 | .439264 | .466576 |
| Attractive pair probe, (-.55,-.45,-.55,.1,0) | 1.122748 | 1.027156 | 1.001624 | 1.009882 |

The last column is 3[S_16(8)-S_14(7)]/log(16/14). The attractive probe is near the U(1) attraction threshold with pairing and chirality added. Its finite-size entropy excess decreases toward one. This is evidence against that particular c>1 candidate, not an exact critical assignment. The interacting collision was not tuned to a thermodynamic multicritical point, so its negative result does not exclude one nearby.

[Exact dispersion evaluated numerically.] At the free collision, evaluate (4) at q=2pi/L:

| L | L e_slow | L^3 e_slow |
|---:|---:|---:|
| 32 | .119966 | 122.845548 |
| 64 | .030207 | 123.727263 |
| 128 | .007565 | 123.950458 |
| 256 | .001892 | 124.006433 |

The last column tends to 4pi^3. This spectral check shows why an entropy estimate near 1/2 is not sufficient to classify the collision as a relativistic Ising point, much less as c=3/2.

### 5.2 Pair binding and the low-K liquid

[Numerical, periodic spin chain L=24, fixed N=1,2,3,4.] At delta=D=0 and t=1, the energy per particle is

| V | E1 | E2/2 | E3/3 | E4/4 |
|---:|---:|---:|---:|---:|
| -1.8 | -2 | -1.993597 | -1.986522 | -1.978980 |
| -2.0 | -2 | -2 | -2 | -2 |
| -2.2 | -2 | -2.012403 | -2.025396 | -2.039249 |
| -4.0 | -2 | -2.500000 | -2.888889 | -3.142857 |

The fixed-N spin ring includes the correct parity-dependent Jordan-Wigner boundary twist. The few-particle results are consistent with the threshold and clustering argument, not with pair stabilization against larger clusters. They are not a finite-density central-charge measurement.

[Numerical, exact projected repulsive gas at rho=1/3.] The structure-factor estimate is K_L=2pi S(2pi/L)/(2pi/L). The independent Casimir estimate uses the known thermodynamic energy and velocity (8), with E0(L)=Le_infinity-pi v c/(6L)+... .

| L,N | K_L | c_E | c_chord |
|---|---:|---:|---:|
| 12,4 | .453251 | 1.018289 | 1.176465 |
| 15,5 | .449292 | 1.011636 | 1.132326 |
| 18,6 | .447377 | 1.008054 | 1.067122 |

All approach K=4/9 and c=1. The chord estimates above one are ordinary short-chain corrections in a model whose one-sea solution is known. This checks the low-K loophole without confusing it with a second fluid.

### 5.3 The two-minimum corner

[Numerical.] `singular_corner.py` diagonalizes (10) at L=24, N<=6, resolved by staggered charge Q. For Q=N, (E+L)/N is 2 to within 3e-14 for every tested N. The mixed sectors lie higher. Examples are 2.034263 for (N,Q)=(2,0), 2.058150 for (4,0), and 2.074802 for (6,0). This tests the positive decomposition (11).

[Numerical.] `factorization.py` tests the full original spin chain at V=4, delta=1, tau=0 and h=sqrt(3). At L=4,6,8, its ground degeneracies are 5,7,9 and its ground energies are -L, agreeing with (13). The defect branch (14) occurs in the dense spectrum to within 3.2e-15. It need not be the first excited level. This tests both the ground-space argument and the nonconformal spectral mechanism.

[Numerical, nearby finite-hopping candidate.] `corner_probe.py` uses (Jx,Jy,Jz,D,h)=(.49,-.51,2,.01,sqrt(15)), giving V=8, delta=1 and tau=.02+.02i at the field of the tau=0 factorization line. No critical tuning is assumed.

| L | c_chord | S_half |
|---:|---:|---:|
| 10 | .867843 | .923248 |
| 12 | .522578 | .907733 |
| 14 | .382496 | .924079 |
| 16 | .262705 | .920164 |

The last half-entropy step gives -0.087956. The ground parity changes between L=12 and 14. These sizes fail to support logarithmic entropy growth, much less c>1. A blind three-parameter half-entropy fit would return about .495 despite that failure, illustrating why such a fit is not a criticality test. This point does not eliminate a finely tuned critical surface elsewhere nearby.

None of these calculations is a DMRG extrapolation or a full four-parameter search.

## 6. Established set, residual set, and next steps

[Scope statement.] Denote the axial parameter space by A. A controlled set C consists of the following, with their displayed hypotheses retained:

1. The whole quadratic locus V=0, restricted to points that possess the conformal infrared assumed by the question.
2. The uniform U(1) locus delta=0 and the zero-field XYZ locus D=h=0, using the supplied complete XXZ/XYZ input. Include their rotated equivalents.
3. Regular weak continuum basins described in section 1.2, and regular continuum perturbations of the known c<=1 interacting loci, with no new soft sector and with the stated RG assumptions.
4. Stable gapped neighborhoods of the large-field, strong-attraction and strong-repulsion classical phases.
5. The regular strongly repulsive constrained-liquid regimes of section 2.2, including the small perturbations whose allowed locking operators were checked.

On conformal members of C, c<=1 is established at the stated physics standard. Items 3 and 5 are operational scaling domains, not explicit global algebraic inequalities in four bare parameters. A numerical radius for them has not been computed. The exact degeneracies (4), (12), the ordinary ferromagnetic points and saturation thresholds are separately identified as nonconformal, rather than counted as additional proved CFTs.

[Open residual.] The remaining conformal candidates lie in A outside C and outside these explicitly nonconformal loci. More concretely this includes finite-coupling nonintegrable points not reached by a controlled one-boson/Majorana scaling regime; possible composite thresholds where the two-body or small-interaction descriptions fail; and neighborhoods of zero-velocity or macroscopic-degeneracy surfaces, particularly |tau| comparable to or below |delta|^2/V near (12). The lack of a uniform weak-coupling radius near these surfaces is part of the residual problem. No exhaustive phase classification, no general lower bound on an unpaired-fluid K, and no symmetry-only codimension obstruction closes this set.

[Ranked next steps.]

1. Resolve perturbations of (12). Derive the competing higher-order terms and the effect of complex tau within its low-energy manifold. The exact projector form and defect spectrum provide calibration points. Check whether a second finite velocity can appear before committing to CFT fits.
2. For a candidate Ising/liquid collision at intermediate coupling, tune the Ising mass and track the separate putative liquid mode. Determine their translation and parity actions, enumerate the genuinely local relevant vertices, and compute the rank of the microscopic coupling map. Entropy alone cannot do this.
3. Use momentum-resolved ED and then converged DMRG for any candidate that passes the preceding checks. Require entropy scaling plus identifiable conformal towers or a Casimir fit with independently determined velocities. Avoid assigning a single velocity to a possible product theory.
4. Independently read and map Chepiga-Mila's conventions to the chiral family before claiming that their phase diagram covers the residual set. No network was used here, so that source check remains undone.

## 7. Reproduction and files

[Numerical protocol.] Run in this directory with Python, NumPy and SciPy:

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 checks.py
    OPENBLAS_NUM_THREADS=1 python3 summarize.py
    OPENBLAS_NUM_THREADS=1 python3 singular_corner.py
    OPENBLAS_NUM_THREADS=1 python3 factorization.py
    OPENBLAS_NUM_THREADS=1 python3 operator_checks.py
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 corner_probe.py

`checks.json`, `summary.json`, `singular_corner.json`, `factorization.json` and `corner_probe.json` contain the records. `PROGRESS.md` is the resumable state. Only small local ED was used. The largest Hilbert spaces are the full spin parity sectors at L=16; the fixed-N routines temporarily enumerate L=24 bit strings. Each computation completed on one machine in minutes or less.

Local sources used for orientation and inherited inputs were `../../STATUS.md`, `../../specs/_d2_background.md`, `../work_C6_d2/c6.tex`, and `../../contrib/d2_classification/README.md`. Sibling progress files were read again at completion. C16 independently treats constrained hopping and the general frustration-free kernel classification; C17 independently audits conversion, anomalies and microscopic sublattice decoupling. Those results were not imported as new numerical checks here. This report's detailed pairing-only effective Hamiltonian, exact line and candidate tests complement their broader analyses. The derivations and numerical claims are unrefereed.
