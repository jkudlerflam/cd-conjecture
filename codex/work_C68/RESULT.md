# C68. Antiunitary and symmetry-free qubit residuals

2026-10-28

[heuristic, verdict] The requested global statement needs correction. There are gapped regions in all three strata, and a protected ordinary Ising transition is locally codimension one with c=1/2. But incommensurate one-boson liquids can occupy open regions even without an internal symmetry. An unprotected Ising point instead requires two relevant coefficients to vanish. Consequently neither "generic points are gapped" nor "the gapless locus is codimension one" follows from the symmetry classification. At the requested theoretical-physics standard, the weak-coupling construction below argues against both statements as global claims.

[proved, symmetry-level conclusion] The surviving symmetries do not prohibit a tuned c>1 multicritical theory. They impose selection rules and, in the time-reversal branch, anomaly matching. They do not impose an upper bound on c. There is no supported c>1 nearest-neighbor qubit point in the supplied results. The most economical mechanism still worth testing is an independent neutral Ising sector becoming critical inside a calibrated floating liquid, with its relevant mixed couplings also controlled. This is a search target, not an identified multicritical point.

[assumed, scope] The chain is H=sum_j[sigma_j^T J sigma_(j+1)+B.sigma_j], with one-site translation, the full qubit tensor product, and Pauli rather than spin matrices. Central charge is the ordinary nonchiral Virasoro value. All claimed CFT modes must have nonzero limiting velocities. "No symmetry" below means no uniform internal unitary or antiunitary symmetry; translation remains exact. Spatial and site-dependent symmetries must be checked separately in a specific example. "Gapped" allows a finite broken-symmetry ground-state multiplet and means a positive bulk gap above it.

[assumed, provenance] The exact XXZ/XYZ continuum identifications, the free range-one theorem, bosonization, and Kramers LSM are inherited inputs from C6 and C63. No external sources were fetched. C63 RESULT.md was read in full, followed by the C6 analysis and the gap-reduction README. The supplementary local sources are listed at the end. Tags apply to the paragraph or displayed derivation they introduce. New algebra is distinguished from conditional continuum reasoning and inherited numerical evidence.

## 1. The three strata and what LSM actually says

[proved, inherited algebra] Write J=S+A(D), S^T=S, A(D)_ab=epsilon_abc D_c. Uniform spin symmetries satisfy

    G_u={R in SO(3): R S R^T=S, R D=D, R B=B},
    G_a={O in O(3): det O=-1, O S O^T=S, O D=-D, O B=B}.

The sign in the antiunitary condition on D is essential. If G_u is trivial, two distinct antiunitaries cannot survive, since their product would be a nontrivial unitary. The single improper Pauli action must square to I. It is either -I or a reflection across a plane. Thus the requested strata are

    T: B=0, G_u={I}, physical time reversal;
    K: B!=0, G_u={I}, S n parallel n, D parallel n, B.n=0;
    N: G_u={I}, G_a empty.

T squares to -1 on a qubit. K squares to +1. After uniform rotations and positive scale their maximal dimensions are 5, 5, and 8. These are dimensions of bare parameter strata, not dimensions of their critical subsets. With simple-spectrum S, its orthogonal stabilizer consists of the eight principal-axis sign choices, giving an exhaustive finite symmetry test. Repeated eigenvalues require the invariant equations above.

[assumed, LSM input and consequence] One Kramers doublet per primitive translation cell forbids a unique, symmetry-preserving gapped state in T. A gapped T phase must break time reversal or translation, or otherwise fail the premises of this statement. The ordered easy-axis example below satisfies LSM by breaking time reversal. K has no Kramers-per-cell obstruction, and N allows a unique product phase. Translation alone does not impose the missing spin anomaly.

[proved, limitation] LSM does not force gaplessness throughout T, does not say that a symmetric critical state is Ising, and does not forbid extra anomaly-free degrees of freedom. Tensoring an anomaly-matching theory with an anomaly-free neutral CFT leaves the anomaly unchanged and increases c. In particular a putative pure-Ising description on T must pass anomaly matching, or be understood inside a background that already breaks a protecting symmetry. The usual Ising mass argument by itself does not settle that issue.

[proved] In an antiunitary-invariant state, the antiunitary pairs translation eigenvalues exp(ik) and exp(-ik), so excitation energies at k and -k coincide. This relates momenta, not independent species. It cannot be used to double or halve c without identifying the fields. In a broken phase it can instead relate spectra above different vacua.

## 2. Gapped regions exist, but their existence is not genericity

[proved, inherited C34 certificate] An explicit anchor in each stratum is

    H=sum_j[a(XX+YY)+w ZZ+d_z(XY-YX)+d(YZ-ZY)+b Z],
    w<0, r=sqrt(a^2+d_z^2), kappa=|w|-r-2|d|.

For b=0 and kappa>0 the two polarized states span the ground space and the gap above them is at least 4 kappa. For b!=0 and kappa>=0 the field selects one polarized state and the gap is at least 2|b|. The proof is the C34 operator inequality

    H-wL >= 2 kappa W+b sum_j Z_j,
    W=sum_j(1-Z_j Z_(j+1))/2.

On a periodic ring W>=2 away from the two polarized states. After subtracting wL-|b|L, the field is 2|b|N_bad>=2|b|(I-P_favored). No commutation of W with H is needed. These are bulk gap bounds, not finite-size extrapolations.

[proved, examples] Set a=0.15, w=-1, d=0.10. With d_z=0.20 and b=0 the point is in T and kappa=0.55, giving gap>=2.20. With d_z=0 and b=0.10 it is in K and kappa=0.65, giving gap>=0.20. With d_z=0.20 and b=0.10 it is in N and kappa=0.55, again giving gap>=0.20. In the first and third cases D has components both along and perpendicular to the unique axis of S, so no proper stabilizer survives. In the second case the reflection normal is x, while D is along x and B along z.

[assumed, local stability] Standard stability of a one-dimensional ordered easy-axis phase, or of a unique gapped product phase, extends strict examples into neighborhoods in the corresponding strata. This input supplies open gapped regions but no explicit radius here. The separate common-axis gap theorem quoted in the gap-reduction README supplies further certified regions. Neither certificate covers every point.

[heuristic, correct local genericity statement] Near an isolated commensurate fixed point with an allowed relevant scalar, an untuned perturbation usually produces a mass. That is a statement about its basin and relevant-coupling coordinates. It is not a statement that the gapped set is open dense in an entire stratum. Section 4 gives a competing local mechanism with an open gapless basin.

## 3. Ordinary critical surfaces and Majorana counting

[assumed, Ising operator algebra] A local spin embedding of the Ising CFT has the relevant scalar primaries sigma, dimension 1/8, and epsilon, dimension 1. Its deformation is

    delta H=int dx [h_sigma sigma+m epsilon+...].

The disorder field is not an additional independent local scalar to tune in the same ordinary spin embedding. A dual embedding can exchange order and disorder. Locality must be specified before counting them. Descendants do not supply additional independent uniform scalar relevant perturbations after total derivatives and redundant terms are removed.

[proved, conditional codimension] If sigma is odd under an exact antiunitary, its Hermitian uniform source is forbidden. Multiplying sigma by i does not evade this, because i sigma is anti-Hermitian. A translation-odd sigma is also forbidden in a one-site-uniform Hamiltonian, even in N. If sigma is uniform and even, its source is allowed. Thus a regular protected Ising point needs m=0, while an unprotected one needs both m=0 and h_sigma=0. In T this count is additional to, not a replacement for, anomaly matching. In K, symmetry-breaking Ising transitions can use an antiunitary-odd order parameter. In N, a transition that spontaneously doubles the translation cell can use a translation-odd order parameter.

[proved, conditional geometry] For a smooth matching map p -> (g_1,...,g_r), the simultaneous zero set has codimension r if it exists and the differentials dg_i are independent there. One allowed Ising mass then gives a hypersurface. Two independent Ising sources give codimension two. Symmetry does not prove that the zero set exists, that the derivatives are independent, or that there are no other soft modes. A first-order level crossing is not a propagating critical sector even though a finite-ring gap can vanish.

[proved, inherited free calculation] At a regular quadratic point use

    H=sum_j[Jx XX+Jy YY+D(XY-YX)+g Z],
    t=Jx+Jy, delta=Jx-Jy,
    lambda_+/- (k)=4D sin k
                 +/-2 sqrt[(t cos k-g)^2+delta^2 sin^2 k].

For 4D^2<delta^2, a simple closing at g=+/-t has one Majorana sector, c=1/2. Simultaneous regular zeros at 0 and pi can give two Majoranas, total c=1. In the regular pocket regime 4D^2>delta^2 there is one quadratic channel with total c=1. Particle-hole-related BdG roots are not additional independent sectors. The supplied range-one rank bound gives c<=1 on the full quadratic locus, including exceptional regular coincidences. Cubic or quadratic zeros at velocity-cancellation and saturation boundaries are not additional finite-velocity CFTs.

[assumed, controlled weak continuum matching] If interactions are weak relative to the retained velocities and to the scale separating all discarded thresholds, a complete one-Majorana theory acquires its allowed magnetic and thermal perturbations, not a new propagating species. A complete two-Majorana theory acquires masses and the dimension-two epsilon_1 epsilon_2 interaction, giving the familiar Gaussian/Ashkin-Teller possibilities or massive flows. This is the appropriate local one-sector argument. "One sector" must mean one complete low-energy theory or fluid, not literally one Majorana, since two Majoranas already occur at c=1.

[proved, scope] These statements do not extend the free polynomial-root bound to arbitrary interacting points. Products of neighboring bonds generate three-site and effective next-neighbor operators; transverse spin sources need not be local even fermion polynomials under Jordan-Wigner. A new threshold or competing vacuum can invalidate the original continuum matching. A c-theorem bound starting at c<=1 is available only for an actual complete unitary relativistic continuum flow from such a UV theory, not for an arbitrary path through lattice couplings.

## 4. An explicit floating-region argument in T, K, and N

[proved, exact lattice anchor] Start with the nearest-neighbor bond

    H0=sum_j[XX+YY+0.4 ZZ+0.5(XY-YX)].

A spiral spin rotation maps its bulk to zero-field XXZ with transverse coupling r=sqrt(1.25), anisotropy Delta=0.4/r, and pitch alpha=arctan(0.5), up to reversing alpha by convention. Boundary twists affect finite-ring zero modes but not the bulk operator dimensions. Product rotations preserve bipartite entropies.

[assumed, XXZ and bosonization input] Put gamma=arccos(Delta). The complete critical theory at this anchor is one compact boson with

    K=pi/[2(pi-gamma)]=0.8110785334,
    v=2 pi r sin(gamma)/gamma=5.4442255453,
    alpha=0.4636476090.

Use H_LL=(v/2pi) int dx[K (partial_x theta)^2+K^(-1)(partial_x phi)^2], with the standard canonical commutator. Local bosonic vertices of the spin-chain continuum algebra have the form

    V_(m,n)=exp(i m theta+2i n phi),
    x_(m,n)=m^2/(4K)+n^2 K,       m,n integers.

This is the ordinary XXZ spin operator lattice. Half-integer choices would introduce nonlocal spinon fields and cannot be inserted as arbitrary local Hamiltonian perturbations. A subset of the integer vertices in a given embedding only weakens the instability test.

[assumed, translation action derived from the spin expansion] In the untwisted zero-magnetization chain, the leading transverse spin field is (-1)^j exp(i theta), the leading longitudinal staggered field is (-1)^j cos(2phi), and the uniform subleading transverse field includes exp(i theta) cos(2phi). Undoing the spiral rotation multiplies a charge-m term by exp(i m alpha j). Physical one-site translation can therefore be represented by

    theta -> theta+pi+alpha,
    phi -> phi+pi/2,
    Q_(m,n)=m(pi+alpha)+n pi mod 2pi.

The choices of cosine versus sine and alpha versus -alpha are conventions and do not affect the following nonresonance test. A uniform lattice perturbation can couple to these fields with spatial phases, but a nonzero Q cuts off its long-distance RG growth.

[proved, conditional vertex enumeration] Positivity of x bounds the relevant list by |m|<=2 and |n|<=1 at the displayed K. Up to conjugation, every vertex with x<=2 is

    (m,n)      dimension             momentum modulo 2pi
    (0,1)      K=0.81107853          pi
    (1,0)      1/(4K)=0.30823156     pi+alpha
    (1,1)      K+1/(4K)=1.11931009   alpha
    (1,-1)     K+1/(4K)=1.11931009   alpha
    (2,0)      1/K=1.23292623        2alpha.

None has zero momentum. The nearest resonance is separated by alpha. The first neutral translation-invariant vertex is cos(4phi), with dimension 4K=3.24431413>2. Derivatives of a charged vertex carry the same nonzero momentum. Dimension-one uniform currents change the background density or pitch; quadratic currents renormalize the Gaussian theory. They are not pinning potentials. No extra local relevant nonderivative term has been omitted within this continuum operator algebra.

[heuristic, weak RG conclusion] For a perturbation g V with dimension x<2 and fixed momentum mismatch Q, its growth g(l) approximately equals g exp[(2-x)l] only until the running momentum reaches the cutoff. In lattice cutoff units this gives a finite condition of the form |g| << const |Q|^(2-x) under which it stays weak before oscillations stop the growth. Every potentially relevant vertex in the list has a finite mismatch, and the first nonoscillating neutral vertex is separated from marginality. Small changes of K and of the momenta preserve these strict inequalities. Products of oscillating operators can renormalize currents and generate neutral harmonics, but the first nonconstant allowed neutral harmonic at the anchor is the irrelevant one already counted. Higher integer charges have increasing scaling dimensions. This gives an open floating basin under standard weak-coupling RG, not merely a line with an exactly conserved microscopic charge.

[proved, membership of three perturbation families] For any sufficiently small nonzero eta use

    S_eta=diag(1+eta,1-eta,0.4),

    T family: D=(eta,2eta,0.5),   B=0;
    K family: D=(0,0,0.5),       B=(eta,2eta,0);
    N family: D=(eta,2eta,0.5),   B=(3eta,4eta,5eta).

Each tends to H0. The three eigenvalues of S_eta are distinct. In T all components of D are nonzero, so only the identity proper rotation and the improper action -I survive. In K the only improper action is diag(1,1,-1), while no nonidentity proper rotation fixes both vectors. In N all components of D and B are nonzero, eliminating every nontrivial proper and improper sign action. These are genuine points of the requested residual strata, not disguised quadratic or uniform-charge points.

[heuristic, implication] The same finite nonresonance margins give relatively open c=1 regions around sufficiently weak members of these families. In N the field can shift the density and pitch; it cannot make a momentum initially separated from zero resonant under an arbitrarily small change. In T the liquid inherits the Kramers anomaly from its XXZ anchor, so it also passes the LSM check. No new anomaly is introduced on K or N. There is no explicit rigorous perturbation radius, and eta=0.01 below is only an algebraic test point, not a simulated or certified thermodynamic liquid. The conclusion is at the requested physics standard, conditional on complete weak bosonization and its RG stability.

[proved, consistency checks within the continuum input] At Delta=0 the formulas give K=1 and v=4r, matching the Pauli-normalized free XX band. As alpha tends to zero, the pair vertex (2,0) loses its momentum mismatch and has dimension 1/K<2 at this K. The floating-region argument then fails in the expected way, and exchange anisotropy can gap the ordinary XXZ liquid. At K=1/2 the neutral umklapp becomes marginal, so the argument also correctly loses its strict irrelevant-operator margin at the antiferromagnetic endpoint. Neither limiting case is silently included in the open neighborhood claimed above.

[numerical, inherited support] The independent local follow-up reports c=1.002 from finite-entanglement scaling for point 220 through correlation length 730, and c=1.008 for point 1282 through correlation length 155. The latter also has low-energy scaling consistent with z=1 and incommensurate momentum about 0.300 times 2pi. These are supplied numerical conclusions, not new calculations here. They support the floating mechanism beyond the small-perturbation construction, but do not themselves prove open neighborhoods. The symmetry test in this run confirms that 1282 has no uniform internal symmetry.

[heuristic, corrected global picture] The gapless set can have interiors with c=1, codimension-one protected Ising surfaces, codimension-two unprotected Ising points or manifolds, and higher multicritical subsets. The boundaries of floating phases can be BKT locking transitions, commensurate-incommensurate thresholds, or first-order boundaries. A BKT separatrix is not the whole gapless set. A quadratic onset has z=2 and cannot be counted as an additional ordinary CFT sector. "Generic gapped away from specified liquid basins and critical sets" is a useful conditional description, but identifying all such basins remains a dynamical classification problem.

## 5. Relevant operators at the remaining candidate fixed points

[assumed, candidate catalog] Ising and the single boson are the established regular theories near the controlled slices. More strongly interacting minimal models, products, and irreducible c>1 theories are possible candidate operator algebras, not established nearest-neighbor qubit realizations. There is no exhaustive classification of all c>1 CFT embeddings here. The following checks therefore cannot be promoted to a complete fixed-point exhaustion.

[proved, conditional Gaussian analysis] For a commensurate boson, cos(2phi) has dimension K, cos(2theta) dimension 1/K, and cos(4phi) dimension 4K. A relevant zero-momentum vertex pins a field and normally gaps that mode. Antiunitarity can forbid one Hermitian combination without forbidding its independent allowed partner, depending on its action on phi and theta. Translation can forbid both. At an incommensurate point one must first test momenta, as in Section 4; a low scaling dimension alone does not imply a gap. Competing dual pinning terms still act on the same canonical pair. They do not by themselves give a boson plus an independent Ising sector.

[assumed, minimal-model data] Tricritical Ising has relevant even fields of dimensions 1/5 and 6/5 and odd fields of dimensions 3/40 and 7/8 in its ordinary spin embedding. A symmetry reversing its order field removes the odd fields but leaves two thermal tunings. Without that symmetry and with uniform local odd fields, all four are allowed. Critical three-state Potts has a thermal field of dimension 4/5 and a conjugate magnetic pair of dimension 2/15. With no microscopic Z3, the magnetic sources need not vanish. A conjugation antiunitary can constrain the two real source components but does not generically remove both. Translation can supply a cyclic action in a period-three ordering problem. Additional local Potts perturbations and their momenta must be tested in that embedding, so the displayed fields give necessary tunings, not a universal complete count. These models have c<1 and illustrate why not every exceptional critical point has one relevant direction.

[proved, conditional product analysis] For Ising times a compact boson,

    c=1/2+1=3/2,
    x(epsilon)=1, x(sigma)=1/8,
    x(epsilon V)=1+x(V), x(sigma V)=1/8+x(V).

Thus epsilon V is relevant when x(V)<1, and sigma V when x(V)<15/8. Momentum and symmetry must be tested on the product, since two individually forbidden fields can form an allowed product. If sigma has momentum pi, a boson harmonic at pi can neutralize it. At the zero-magnetization anchor, that harmonic has dimension K, so such a product has dimension K+1/8=0.93607853 when its internal transformation permits it. The associated coupling cannot be ignored.

[proved, further mixed terms] Sigma times a dimension-one current has dimension 9/8 and is relevant when local and invariant. An antiunitary-odd sigma can multiply an antiunitary-odd current to form an even perturbation. Epsilon times a current has dimension two and requires a marginal RG analysis, including the actual velocities. Bilinear mixing between fermionic fields has to pass the spin-locality and momentum tests. A uniform sigma requires its own magnetic tuning if unprotected. These observations invalidate the assertion that a single neutral mass is always the only condition for a c=3/2 point inside a liquid.

[heuristic, economical multicritical mechanism] If sigma is translation-odd, all compensating low-dimension boson harmonics are nonresonant, and the marginal mixed terms are harmless, an independent Ising mass can cross zero inside a floating phase with only one relevant tuning. The total theory then has two stress tensors and c=3/2, despite codimension one. If the Ising order field is uniform and unprotected, at least its magnetic and thermal couplings must be tuned, often with additional mixed terms. Neither case is excluded just by having five or eight microscopic parameters. A concrete matching map is needed to determine reachability and rank.

[assumed, three-Majorana candidate] At a local three-Majorana realization the operators i xi_R^a xi_L^b have dimension one. Translation and internal actions select which entries of the real mass matrix M_ab are local and invariant. Without a protecting relation they can split the massless modes, so simultaneous criticality is unstable to generic perturbations. To retain all three modes all allowed relevant mass components must vanish, not just det M. Quartic terms are marginal and can further destabilize the point. The exact list of independent tunings depends on momentum and locality. The free range-one qubit theorem forbids constructing this spectrum as three microscopic quadratic channels, but does not forbid an emergent interacting realization.

[assumed, other c>1 examples] At an SU(2)_k WZW point the diagonal singlet in the spin-one primary has dimension 4/(k+2) for k>=2. Where translation permits it, it is a relevant thermal perturbation even if full SU(2) is retained. With less symmetry more components and primaries can be allowed. For a Z_k parafermion candidate the first thermal dimension is also 4/(k+2), while order fields can have dimensions l(k-l)/[k(k+2)]. At k=5, for example, c=8/7 and the first order dimension is 4/35. Emergent symmetries require tuning their permitted breaking fields. A relevant perturbation makes a fixed point unstable, not impossible. Orbifolds and extensions can change the local operator list. No inference excluding all irreducible c>1 theories follows from the familiar Ising and Gaussian list.

## 6. Is c>1 symmetry-permitted?

[proved, distinction] Symmetry-permitted means that a continuum theory can obey the specified symmetries, with every allowed relevant coupling set to its critical value, and match the required anomaly. It does not mean that some nearest-neighbor qubit bond is known to reach it. Absence of protection and impossibility are different claims.

[assumed, explicit continuum compatibility in T] Take the spin-half XXZ boson with its translation and time-reversal action, and tensor it with an anomaly-free neutral Ising theory on which time reversal acts by complex conjugation. The boson carries the Kramers-per-cell obstruction. The neutral factor changes c to 3/2 without changing that obstruction. Its thermal and magnetic fields and any allowed mixed fields must be tuned or forbidden. This is a consistent continuum symmetry assignment. It is not a construction on one physical qubit per nearest-neighbor cell.

[proved, implication for K and N] K has no such projective-cell obstruction, so an anomaly-free Ising-boson product can carry a square-positive antiunitary. N imposes no internal selection rule at all. Neither symmetry excludes the same abstract product, nor an irreducible c>1 theory. An emergent higher symmetry at an isolated infrared point does not require the microscopic Hamiltonian to have that symmetry. Conversely exact microscopic SU(2) on a two-qubit bond is only identity plus swap and reduces to Heisenberg, as C6 proves. That microscopic closure does not apply to emergent SU(2).

[proved, what is missing] The absence of a lattice realization is not settled by counting relevant couplings. A proof would need a rank obstruction for every matching map, a complete dynamical phase classification, or another microscopic theorem controlling all neutral sectors. C63 already identifies this missing step. In particular even excluding Ising times boson would leave irreducible theories with 1<c<3/2 when there is no exact current algebra requiring a charged c=1 factor.

## 7. Numerical target and a test that could decide it

[heuristic, priority] No existing point can honestly be nominated as a supported c>1 candidate. The strongest mechanism to test is a neutral transition in the interior of a verified floating phase, away from a zero-velocity saturation or factorization boundary. Point 1282 is a concrete symmetry-free starting point, chosen for its independently calibrated one-boson scaling and cleaner reported z estimate. It is a c=1 control, not itself an excess-central-charge signal.

[numerical, parameter extraction and algebra only] In the conventions of Section 1, the supplied parameter record is

    J0 = [[ 0.548197257891, -0.260919803070, -0.465109070557],
          [ 0.260919803070, -0.038083786052, -0.026441594586],
          [ 0.465109070557,  0.026441594586, -0.327912890974]],
    B0 = [-0.072915926520, 0.039770388744, -0.115735099721].

The unrounded coefficients are saved in checks.json. Its S is diagonal with three distinct eigenvalues, and D and B have three nonzero components. The finite sign test verifies membership in N. This run extracted this one parameter record; it did not read the simulation data or rerun VUMPS.

[heuristic, concrete search family] A reproducible first continuation is

    J(s)=J0+s e_z e_z^T,
    B(t)=B0+t e_z.

Use s to change an exchange channel and t to control a second relevant response, with (s,t)=(0,0) the calibrated liquid. This stays in N for generic nearby values. There is no derived location or evidence that this particular plane contains a neutral transition; it is a bounded starting search rather than a prediction. Track all local correlation branches while staying in the c=1 region. If an independent neutral branch softens, use additional tangent directions of S, D, and B to tune the relevant response matrix. Two coordinates need not suffice once mixed terms are found.

[heuristic, antiunitary counterpart] On K the underresolved tilted-field easy-axis XXZ boundary in the supplied multicritical study remains a useful targeted continuation. Take S=diag(1+eta,1-eta,2), D=(0,d,0), B=(h_x,0,h_z), with nonzero transverse field and generic small eta,d. The sole antiunitary reflects y. At eta=d=0 this connects to the studied real tilted-field slice. Determine its ordered boundary with larger systems before asking whether a neutral critical sector survives inside a liquid. The old L<=16 map did not locate the junction, and its anticipated commensurate-incommensurate onset is not evidence for c>1. On T there is no comparably calibrated multicritical seed in the supplied work; Section 4 provides a liquid anchor for a future boundary continuation.

[assumed, measurements that discriminate] A successful Ising-boson candidate must show (i) a converged bulk entropy coefficient approaching 3/2 rather than a short-ring chord overshoot, (ii) an independently identified Ising tower including its local-field assignments and descendants, (iii) a boson oscillator and zero-mode tower, with nonzero limiting velocities for both, and (iv) control of every relevant invariant product and the marginal mixed flow. A second correlation peak can be a harmonic of one boson. A small gap can be a moving zero mode. Neither is an extra sector.

[proved, response and energy diagnostics] For an operator of scaling dimension x, the matrix elements of its spatial integral scale as L^(1-x) between fixed conformal levels. Measuring independent bond perturbations in the candidate low-energy subspace can therefore identify relevant response directions before retuning them. For independent sectors with possibly unequal chiral velocities,

    E0(L)=e_inf L-pi/(12L) sum_a c_a(v_(a,R)+v_(a,L))
          +zero-mode terms+o(1/L).

A single assumed velocity does not extract the total c. The supplied C63 factorization boundary gives a direct warning: a half-system entropy-size estimator tends to 3/2 while an exact excitation branch is quadratic and the ground-space dimension grows with L. The C6 symmetry-free boundary also has finite-size evidence for quadratic softening. These are poorer c>1 targets than an independent neutral mode inside a finite-velocity liquid.

## 8. Checks, limitations, and handoff

[numerical, new checks] Running python3 checks.py in this directory verifies the three explicit residual symmetry assignments at eta=0.01, enumerates every Gaussian vertex with dimension at most two at the anchor, checks their nonzero momentum mismatches, and verifies the absence of a uniform internal symmetry at point 1282. It also checks the spiral rotation directly on the two-spin bond, with maximum entrywise error 1.7e-16, and the free-XX normalization. The five representative vertex pairs are (0,1), (1,-1), (1,0), (1,1), and (2,0). The smallest momentum mismatch is 0.4636476090 and the first neutral umklapp dimension is 3.2443141338. All assertions passed. These are finite algebra and formula checks, not new thermodynamic evidence.

[proved, completion scope] New files are RESULT.md, PROGRESS.md, STATE.md, checks.py, and checks.json. All writes and Python execution stayed in this folder. No network command, cluster job, large local simulation, or parent-file modification was performed. The project-wide STATE.md is outside the authorized write scope, so the handoff is recorded in the local STATE.md instead. No independent referee has checked this analysis.

[assumed, local sources] The sources used are ../work_C63/RESULT.md, ../work_C6_d2/c6.tex, ../../contrib/d2_gap_reduction/README.md, the relevant gap and residual arguments in ../work_C34/RESULT.md, and the README summaries in ../../contrib/d2_followup_220_1282/, ../../contrib/d2_classification/, ../../contrib/d2_multicritical/, and ../../contrib/d2_two_fluid/. C16 and C17 were consulted for their anomaly and operator-count qualifications. The numerical statements are attributed to those local summaries, not promoted to new proofs or a fresh literature audit.

[heuristic, final status] The strongest conclusion is a local one. Familiar protected critical surfaces have one Ising sector, and controlled liquid neighborhoods have one boson. Open gapped regions also exist. The full gapless locus is not restricted to hypersurfaces, and fine-tuned c>1 remains symmetry-permitted but microscopically unestablished. The next useful computation is a neutral-sector search inside a calibrated floating phase with a measured relevant-coupling response, not another unconstrained maximization of a short-size entropy fit.
