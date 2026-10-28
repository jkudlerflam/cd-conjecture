# C63. A stratified reduction and a test of coincident sectors

2026-10-28

This calculation does not prove the full qubit bound and does not produce a supported c>1 candidate. It gives a covering with the unresolved strata explicitly retained, audits the proposed sector-counting lemma, and tests the reduced staggered-field problem with small exact diagonalization. The most concrete additional result is an exact entropy false positive on its factorization boundary. A size-entropy estimator tends to 3/2 there, while an exact excitation branch has z=2.

[proved] The proposed statement that a nearest-neighbor qubit chain cannot have two simultaneous critical sectors is false as written. Two massless Majoranas already occur in the nearest-neighbor XY chain and have total c=1. The useful question is whether the total conformal content can exceed one. Excluding Ising times Luttinger would address an important mechanism, but would not exclude irreducible interacting CFTs with 1<c<3/2.

[assumed, scope] Central charge means the ordinary nonchiral Virasoro value, so a critical Ising chain has c=1/2, and independent sectors add. All modes being counted must have nonzero limiting velocities. Unitary compact continuum theories and the standard XXZ/XYZ identifications are the same inputs as in C6 and C15 to C17. Existing exact gap theorems are inherited from the supplied notes, not independently checked against external papers. All work was offline, using Python, NumPy, and SciPy. No file outside this folder was modified.

The tags [proved], [numerical], [heuristic], and [assumed] apply to the paragraphs or tables they introduce. An unresolved row with a proposed mechanism is tagged [heuristic], not promoted to a proof. These calculations have not had an independent referee.

## 1. A covering of the full parameter family

[proved, inherited algebra] Remove constants and telescoping bond terms and write

    H = sum_j [sigma_j^T J sigma_(j+1) + B.sigma_j],
    J = S + A(D),  S^T=S,  A(D)_ab=epsilon_abc D_c.

There are twelve real coefficients. Uniform proper rotations remove three generic coordinates and positive scale removes one. A convenient normalization is

    tr S^2 + 2 |D|^2 + |B|^2 = 1.

The quotient P has generic dimension eight. Negative scale is not an equivalence of ground-state problems. Repeated exchange eigenvalues give lower-dimensional strata and larger stabilizers. The invariant tests below remain valid there.

[proved, inherited algebra] Define

    F = {exists unit n: S n=0, D parallel n, B parallel n},
    U = {exists unit n: S=a(I-n n^T)+u n n^T,
                         D parallel n, B parallel n},
    X = {D=B=0}.

F is the complete ordinary quadratic Jordan-Wigner locus, including arbitrary product rotations before that map. U is the uniform additive-charge locus. X is zero-field XYZ. Their maximal dimensions are 3, 3, and 2.

[proved, inherited gap certificates] Let G include both of the following rotation-invariant families.

1. There is a common left and right axis n with

       J n=J^T n=lambda n,
       lambda <= -||(I-n n^T) J (I-n n^T)||op,
       B.n != 0.

   The imported Bravyi-Gosset-Liu-Wong theorem gives a unique ground state and gap at least 2|B.n|, at every size. This includes symmetric ferromagnetic rank-one exchange in arbitrary longitudinal and transverse fields.

2. After rotation, the Hamiltonian has C34's form

       a(XX+YY)+w ZZ+d_z(XY-YX)+d(YZ-ZY)+b Z,
       w<0,  kappa=|w|-sqrt(a^2+d_z^2)-2|d|.

   For b=0 and kappa>0 the gap above the two polarized vacua is at least 4 kappa. For b!=0 and kappa>=0 the ground state is unique and the gap is at least 2|b|. This elementary certificate reaches parts of all three residual branches with trivial unitary symmetry.

These are sufficient tests. Failing a test does not imply criticality. Known stable neighborhoods of uniformly gapped Hamiltonians can be added when their stability hypotheses are established.

[proved, inherited charge classification] Outside U and the zero-field controlled classes, an exact additive charge sum_j n_j.sigma_j leaves only

    A_st = {H(a,u,b)=sum_j [a(XX-YY)+u ZZ+b Z]}.

For a!=0 and b!=0, rotations and positive scale set a=1 and b>0. Let r=sqrt(1+b^2/4). The wings u<=-1 and u>r are gapped, and u=r has L+1 vacua and a quadratic branch on even rings. The line u=0 is quadratic and gapped for b>0. The unresolved part is

    A_st^o = {-1<u<r, b>0, u!=0}.

This is a two-dimensional subfamily of the four-dimensional unitary-Z2 class. It is not the generic residual of the eight-parameter problem.

[proved, covering] Read the following table in order and remove all preceding rows from each row. For the staggered row, split off its stated wings and boundary before taking A_st^o. The resulting pieces are disjoint and cover P. Dimensions are maximal dimensions before these removals, not codimensions of unknown critical loci.

| Piece | Dimension | Mechanism controlling the settled portion | Status on the whole piece |
| --- | ---: | --- | --- |
| F | 3 | Range-one quadratic Majorana count, c<=1; nonlinear zeros are outside the finite-velocity claim | [proved], using the supplied free theorem |
| U | 3 | Spiral gauge transformation to uniform-field XXZ; one Bethe sea when conformal | [assumed], standard complete XXZ input |
| X | 2 | XYZ gap classification; regular critical loci reduce to one Gaussian theory | [assumed], standard exact-model input |
| G | up to 6 within its explicit families | Uniform gap or gap above finitely many vacua | [proved], with the stated imported theorem for certificate 1 |
| A_st, its settled wings and boundary | 2 | Gap theorem, positive-bond decomposition, quadratic defect branch | [proved], inherited |
| A_st^o | 2 | One-boson sine-Gordon theory at weak staggered field; small ED in Section 5 | [assumed] in regular weak-field regimes; [numerical] at the tested points; globally unresolved |
| Exactly one nontrivial uniform proper spin rotation, a pi rotation | 4 | Interacting spinless pairing chain; weak Majorana/boson theories and specific gapped corners | [heuristic] global one-sector exhaustion is missing |
| Trivial unitary symmetry, time reversal only | 5 | Kramers LSM constraint; C34 controls an ordered wing, regular continuum neighborhoods remain useful | [heuristic] neutral-sector exclusion is missing |
| Trivial unitary symmetry, plane-reflection antiunitary only | 5 | Real spin-plane symmetry; gap certificates and some sampled tilted-field slices | [heuristic] neutral-sector exclusion is missing |
| No nontrivial uniform internal symmetry | 8 | Regular one-mode scaling regimes, gap certificates, and tested floating phases | [heuristic] no complete infrared exhaustion |

[proved, completeness of the last four rows] Uniform unitary and antiunitary actions satisfy

    G_u = {R in SO(3): R S R^T=S, R D=D, R B=B},
    G_a = {O in O(3), det O=-1:
                          O S O^T=S, O D=-D, O B=B}.

A proper rotation other than a pi rotation forces the continuous axial symmetry already in U. Two different surviving pi axes force D=B=0, already in X. If G_u is trivial, two distinct antiunitaries are impossible because their product would be a nontrivial unitary. The surviving improper involution is either -I, requiring B=0, or a plane reflection with normal n, requiring

    S n parallel n, D parallel n, B perpendicular n.

If neither exists, the last row applies. This proves an algebraic covering even at repeated eigenvalues. It does not prove a dynamical mechanism for every point in an unresolved row. Spatial symmetries and staggered charges are additional restrictions, not omitted uniform internal symmetries.

[proved, logical status] A successful global covering argument would need a controlled neighborhood at every remaining point, or a complete phase classification on each remaining row. Neither is supplied by the existing work. Compactness of P can extract a finite subcover only after those neighborhoods have been established. It does not create them or bound the crossover length. The table is the covering that can currently be justified, with its missing assumptions visible.

## 2. Majorana counting and its perturbative domain

### 2.1 The whole quadratic locus

[proved, C6/C15 conventions] On F use

    H=sum_j [Jx XX+Jy YY+D(XY-YX)+g Z].

With n_j=c_j^dagger c_j=(1-Z_j)/2, its quadratic terms are hopping with coefficient Jx+Jy-2iD and real pairing Jx-Jy. Put

    t=Jx+Jy, delta=Jx-Jy.

The BdG eigenvalues, up to reversing k, are

    lambda_+/- (k)=4D sin k
          +/- 2 sqrt[(t cos k-g)^2+delta^2 sin^2 k].       (1)

There are two real Majoranas per site and the intersite Majorana matrix has rank at most two. The inherited range-one theorem gives c<=rank A(1)/2<=1 whenever the limit is conformal.

For an explicit mode count, put s^2=4D^2-delta^2.

* If s^2<0, a zero needs sin k=0 and g=t cos k. For t!=0 this is a single Ising point on g=+/-t. At t=g=0, both 0 and pi are massless, and the two Ising contributions sum to one.
* If s^2>0, zeros solve

      (t x-g)^2=s^2(1-x^2), x=cos k,
      x=[t g +/- s sqrt(t^2+s^2-g^2)]/(t^2+s^2).

  In the regular liquid they describe the available chiral Fermi edges of one quadratic channel. Particle-hole-related BdG zeros are not additional physical species. The total c is one. At |g|=sqrt(t^2+s^2) a pocket closes with a quadratic edge.
* If s^2=0 and delta!=0, one velocity can vanish. The generic |g|<|t| touching is quadratic, and g=+/-t can be cubic. Flat limits also occur. They cannot serve as regular relativistic starting points.
* For delta=D=0, the interior of the single cosine band has c=1 and its saturation edge has z=2. Completely flat exceptional Hamiltonians have no ordinary propagating CFT interpretation.

These cases include every quadratic point, but only the regular cases admit the following perturbative argument.

### 2.2 A concrete counterexample to the overly strong lemma

[proved mapping; assumed Ising infrared identification] Take the allowed bond

    H=-J sum_j X_j X_(j+1)-g sum_j Z_j Z_(j+1).

With dual variables tau_j^z=X_j X_(j+1) and Z_j=tau_(j-1)^x tau_j^x,

    H_dual=-J sum_j tau_j^z-g sum_j tau_(j-1)^x tau_(j+1)^x.

At J=g>0 the even and odd dual lattices are both critical transverse-field Ising chains. Global duality constraints and twisted sectors affect the finite-ring state organization, but not the sum c=1/2+1/2. A uniform rotation also identifies the original Hamiltonian as XY. Independence here means two propagating Majorana sectors in the bulk, not two tensor factors on alternating physical sites. Literal physical-sublattice decoupling with range one would force every connected bond coefficient to vanish, as C17 proves.

### 2.3 Which perturbations can introduce a new sector?

[assumed, controlled RG] Fix a regular free point with smallest retained speed v_min>0 and a nonzero matching scale E_sep separating the retained low-energy momentum regions from the rest. Weak local interactions relative to v_min, with renormalizations small relative to E_sep, admit a continuum expansion in the retained fields.

At a single Ising point, the thermal mass has dimension one. A symmetry-allowed magnetic perturbation has dimension 1/8 when it is local in the spin embedding. Retuning a mass or adding a magnetic perturbation uses the existing Ising sector. The leading local scalar four-fermion contact built from one Majorana species vanishes without derivatives; the remaining irrelevant terms do not create a third fermion branch.

At two-Majorana points, the masses m_a i xi_R^a xi_L^a, a=1,2, are relevant, while epsilon_1 epsilon_2 has dimension two. The latter is the Gaussian/Ashkin-Teller interaction when the corresponding local operator embedding is used. It changes the stiffness or produces a massive flow, within the same complete c=1 continuum theory. Boundary constraints can distinguish the circle and orbifold operator spectra without changing this central-charge count.

For a regular Dirac liquid, bosonize into a single canonical pair phi,theta. Forward scattering changes K and velocities. The leading familiar vertices are

    cos(2 phi), dimension K,
    cos(2 theta), dimension 1/K,
    cos(4 phi), dimension 4K,                              (2)

when momentum, locality, and exact symmetries permit them. They pin combinations of the same fields. Multiple transitions of these fields need not be independent sectors.

[proved conditional on a relativistic continuum flow] If the entire subsequent flow has a unitary Lorentz-invariant ultraviolet continuum fixed point with c<=1, the c-theorem gives c_IR<=1. A common drift can be removed by a frame change. Arbitrary unequal velocities are not covered by that theorem without more input. Weak one-channel bosonization remains the stated physics input in those cases.

[assumed, exact-model continuation] The same argument applies near every regular critical point on the supplied XXZ, twisted XXZ, and XYZ slices after their known continuum theories have been identified. Their regular conformal spectrum has one boson, including the SU(2)_1 endpoint. Gapped interiors have stable gapped neighborhoods. Saturation edges, ferromagnetic multiplets, the factorization line, and flat or cubic free points are separate singular cases. No complete classification of all possible integrable Hamiltonians in P is being assumed.

[proved, limitation] A second independent neutral sector requires an additional threshold becoming soft, a competing vacuum, or failure of the assumed complete continuum matching. It cannot be obtained by perturbing the already complete regular theory while all discarded thresholds remain finite. This local observation is not a global exclusion. A trajectory through bare Hamiltonians is not an RG trajectory.

The bare two-Majorana count does not survive as an interacting theorem. Products of neighboring bonds already have longer support,

    (sigma_i^a sigma_(i+1)^b)(sigma_(i+1)^c sigma_(i+2)^d)
       =delta_bc sigma_i^a sigma_(i+2)^d
         +i epsilon_bce sigma_i^a sigma_(i+1)^e sigma_(i+2)^d.

Effective hopping therefore need not have the Fourier degree of the microscopic quadratic hopping. C15 gives an actual second-order two-site hop in the pure-pairing strong-coupling corner. Its leading minima happen to lie at a ferromagnetic degeneracy, but the example already defeats a general argument based only on bare Fourier degree. Symmetry-breaking spin terms can also become Jordan-Wigner strings, so the generic residual is not a weak local fermion polynomial.

## 3. Symmetries, anomalies, and a sharper conditional reduction

[assumed, standard anomaly input] A single on-site Z2 has no nontrivial projective-cell LSM obstruction of the spin-rotation type. The two pi rotations in D2 can act projectively on a qubit. Physical time reversal has T^2=-1 on each qubit and gives a Kramers-per-cell LSM obstruction with translation. The spin-plane reflection antiunitary has square +1 and has no such Kramers obstruction. Without an internal symmetry, translation alone allows a unique product ground state.

With exact SO(3) and the usual antiferromagnetic translation action, a pure SU(2)_k WZW theory matches the spin-half LSM index when k is odd. This is a parity condition, not k<=1. The actual uniform SU(2)-invariant nearest-neighbor qubit bond is already Heisenberg because two qubits decompose into singlet plus triplet without multiplicities. Its complete microscopic reduction, rather than anomaly matching, closes that slice. An emergent SU(2) in a symmetry-free Hamiltonian cannot be assigned a fundamental microscopic embedding index without an additional derivation.

[proved consequence] An anomaly-free neutral CFT can be tensored with an anomaly-matching sector without changing its anomaly. Therefore none of these anomaly constraints upper-bounds total c. Breaking a protecting symmetry allows additional relevant fields, but isolated fine-tuned critical points are still part of the requested claim.

[proved, conditional operator count] For Ising times one boson, epsilon has dimension one and sigma has dimension 1/8. A uniform local sigma is forbidden if it is odd under a surviving unitary symmetry, odd under an antiunitary, or carries translation momentum pi. Its thermal mass is still allowed. Without such protection, both its magnetic and thermal coefficients must be tuned. These are counts of continuum constraints, not known microscopic codimensions.

If a local boson vertex V has dimension d_V, epsilon V and sigma V have dimensions 1+d_V and 1/8+d_V. Mixed products can be invariant even when their factors separately are not. Allowed sigma times a boson current has dimension 9/8. Epsilon times a current has dimension two and needs a marginal-flow analysis. Bilinear mixing of a neutral Majorana and a fermionic liquid field also needs an embedding-specific locality and momentum test. A mere list of symmetry names does not decide these questions.

If continuum couplings g_alpha(p) vanish at some point p and their Jacobian has rank r, their regular zero set has codimension r. Without a known zero and a matching map p -> g_alpha(p), counting tunings neither proves nor excludes the zero. Four, five, and eight microscopic parameters are not too few for a protected Ising mass crossing inside a liquid.

### A useful charged-sector statement

[assumed, standard unitary CFT input] Suppose the exact lattice U(1) acts nontrivially on the gapless fields and yields a nonzero chiral current J with

    J(z) J(0) ~ k/z^2, k>0.

The Sugawara tensor T_J=:JJ:/(2k) has c_J=1. Subtracting it from the physical stress tensor gives a commuting unitary Virasoro tensor,

    T_neutral=T-T_J, c_neutral=c-1 >=0.                    (3)

For a unitary Virasoro vacuum representation with c_neutral<1, the allowed values are 0 and 1-6/[m(m+1)], m=3,4,... . The smallest positive value is 1/2.

[proved from that input] Consequently a charged compact CFT cannot have 1<c<3/2. If 1<c<=3/2 it must have c=3/2 and an Ising neutral stress tensor. The full local theory need not be a literal tensor product. SU(2)_2, for example, is a U(1)-Ising extension with correlated sectors. This strengthens the test for the staggered-charge family, but is not an exclusion of that endpoint.

The hypothesis matters. A microscopic charge can act trivially on all gapless degrees of freedom while charged excitations remain massive. Then (3) supplies no c=1 sector. Likewise, for several velocities, one must establish the charged and neutral factor decomposition before applying the corresponding sum of central charges. No common-velocity Casimir fit proves it.

[proved, staggered symmetry] For A_st, Q=sum_j (-1)^j Z_j is conserved, and physical translation obeys T Q T^-1=-Q. Alternate-site X rotations map Q to uniform magnetization. Physical translation becomes a one-site shift combined with a global X flip, and translation by two sites is ordinary translation. This semidirect action must be kept in a continuum assignment. The explicit staggered field in the transformed Hamiltonian does not break the original one-site translation. It also does not impose the usual one-site uniform-filling LSM constraint on the transformed two-site cell.

## 4. Universality-class catalog and its limits

[assumed, CFT identifications; heuristic, microscopic compatibility] There is no known exhaustive classification of unitary compact CFTs in 1<c<=3/2 with all possible microscopic embeddings. Nor does the local Hilbert dimension provide a theorem selecting a finite catalog. The following is a useful catalog of candidate families, not a claim that all rows occur in range-one qubit chains. Calling an unknown realization "compatible" can mean only that the listed symmetry arguments do not rule it out.

| Continuum family | Central charge | Range-one qubit status and relevant issue |
| --- | --- | --- |
| Ising | 1/2 | Established in transverse-field Ising and regular XY critical surfaces |
| Unitary Virasoro minimal models M(m,m+1) | 1-6/[m(m+1)], m>=3 | The c<1 central-charge series is complete under the stated CFT hypotheses; its qubit NN realization list is not. Beyond Ising, no realization is established by this work |
| Tricritical Ising, three-state Potts, higher minimal models | 7/10, 4/5, ... | Special members of the preceding series. Extra relevant tunings or emergent symmetries are possible; the standard larger-spin or extended-range realizations do not establish this microscopic class |
| Compact boson, Gaussian/eight-vertex points, SU(2)_1 | 1 | Established on XXZ/XYZ and at controlled or numerically identified floating phases |
| Two Ising stress tensors, Ashkin-Teller/orbifold-related theories | 1 | Two Majorana sectors occur in XY. The detailed local operator algebra and orbifold identification require the specific mapping; no blanket claim that every c=1 theory is realized |
| Ising times a compact boson, or three Majoranas | 3/2 | Principal neutral-sector target. Extended-range qubit examples in the supplied literature are not range-one examples |
| SU(2)_2, equivalently the three-Majorana affine theory | 3/2 | Allowed as an emergent possibility in a residual stratum; exact microscopic SU(2) instead reduces to Heisenberg. The usual spin-half anomaly excludes the pure even-level theory only with the specified SO(3) and translation action |
| Products of minimal models | c_m+c_n<=3/2 | Ising times any finite minimal model, TCI times TCI (7/5), and TCI times Potts (3/2) are concrete targets. Three nontrivial minimal factors fit only as three Ising factors at 3/2. None is constructed here beyond the two-Majorana c=1 case |
| Z_k parafermions | 2(k-1)/(k+2), 2<=k<=10 | k=5,...,10 give 8/7, 5/4, 4/3, 7/5, 16/11, 3/2. Emergent parafermions are not ruled out by having only microscopic Z2 or no symmetry; a clock-model on-site dimension is not a proof of a qubit lower bound |
| Unitary N=1 superconformal minimal series | (3/2)[1-8/(m(m+2))], m>=3 | Includes TCI, c=1, and an infinite sequence in (1,3/2). Supersymmetry can be emergent. The needed local bosonic realization or spin structure must be specified |
| Finite orbifolds, extensions, and other cosets | Inherited or model-dependent c | They can change locality and the relevant-operator list without changing a known central charge. The preceding rows are not exhaustive above one |

[proved, significance] Even a proof excluding all simultaneous Ising and Gaussian sectors would leave parafermionic and other non-product theories. Conversely, none of these abstract CFT possibilities supplies a lattice counterexample. For the charged case of Section 3, the nontrivial-current hypothesis removes every value strictly between one and 3/2, leaving the Ising-coset endpoint. For the symmetry-free residual it does not.

[assumed, literature provenance] The supplied numerical summaries report no known strictly nearest-neighbor, one-site-cell qubit chain with c>1. Their cited Ising-plus-liquid and cluster-Ising realizations use longer-range or three-site terms. That is the inherited literature status, not a new exhaustive search. No outside source was fetched in this run.

## 5. The staggered-field test

### 5.1 Weak field and the strongest singular candidate

[proved, exact lattice map] On even rings, alternate-site X rotations give

    H'=sum_j [X_j X_(j+1)+Y_j Y_(j+1)
                -u Z_j Z_(j+1)+(-1)^j b Z_j].            (4)

This is spinless nearest-neighbor hopping of magnitude 2, density interaction -4u, and alternating chemical potential -2b(-1)^j, plus constants and the half-filling centering term. The doubled cell is a representation of the staggered charge, not two independent low-energy fluids.

At u=0 the two folded one-particle bands are

    E_+/- (p)=+/-sqrt(4b^2+16 cos^2 p),                   (5)

with reduced physical momentum p. They are gapped for b!=0. At b=0 they join the single free band. Folding its Brillouin zone does not double c.

[assumed, exact XXZ parameters and bosonization] For b=0 and -1<u<1, the anisotropy of (4) is Delta=-u and

    K=pi/[2 arccos(u)],
    v=2 pi sqrt(1-u^2)/arccos(-u).                        (6)

These Pauli-normalized formulas give K=1,v=4 at u=0; K=3/2,v=3sqrt(3)/2 at u=1/2; and K=2 at u=1/sqrt(2). The field b couples to staggered density, hence to cos(2 phi) up to a shift of phi. Its dimension is K. Infinitesimal b is relevant for K<2, and irrelevant for K>2, with a BKT separatrix near K=2. The half-filled umklapp cos(4 phi) has dimension 4K. Exact charge conservation forbids theta sources. This establishes a regular weak-field one-boson description, not the full finite-b phase diagram.

[heuristic, target selection] The strongest supplied mechanism for an extra sector in this family is the approach to the upper boundary u=r, especially at larger b. Virtual pairing has two low-momentum minima in the original variables, and the putative liquid interval becomes close to a highly degenerate line. C15 already shows that its leading strong-coupling theory is a ferromagnetic degeneracy, not two independent liquids. The present test probes both a weak-field liquid and finite-field points between the apparently massive side and this boundary. It is not an optimization over all eight parameters.

### 5.2 What was computed

[numerical] `ed_staggered.py` constructs the sparse matrix of (4) at fixed N, with Z=1-2n, in the computational basis. Off-diagonal exchanges of unequal neighboring spins have matrix element 2. The diagonal is

    -u sum_j z_j z_(j+1)+b sum_j (-1)^j z_j.

The main scan uses L=8,10,12,14,16,18 at nine points. Four points are extended to L=20. At each point it obtains fourteen low levels in N=L/2 and the ground energies in N=L/2+1 and N=L/2+2. Every conserved-charge sector was independently checked at L=18; N=L/2 is the global minimizer at all nine points. The combined spin-flip and translation symmetry relates N and L-N. The global-sector check at L=18 is not a theorem at other sizes.

Translation by two sites is diagonalized within numerically degenerate eigenspaces. Momentum labels are used only when the projected translation eigenvalue has modulus greater than 0.999. Truncation through a degenerate multiplet can otherwise spoil a label. Define the candidate oscillator speed from the lowest resolved relative T2 phase +/-4pi/L by

    v_L=L [E_ph(L)-E0(L)]/(2pi).                          (7)

The physical momentum is 2pi/L, so the denominator in (7) is 2pi, not 4pi. Calling this state a boson oscillator is a spectral hypothesis. If a slower additional sector exists, the lowest state at this momentum need not be the charge oscillator.

For entropies, average the two possible block starting parities and fit

    S(l)=s0+(c_chord/3) log[(L/pi) sin(pi l/L)]

over floor(L/4)<=l<=L/2. The size estimator c_size is three times the slope of the averaged half-chain entropy versus log L, using available sizes L>=12. These are distinct fits to the same entanglement observable, so they are supplemented by spectra and energy scaling.

The energy fit uses E0/L=e_inf+A/L^2+B/L^4 over L>=10. Its direct output is the velocity-weighted quantity W=-6A/pi. A separate fit v_L=v_inf+a/L^2+d/L^4 gives c_E=W/v_inf only under the one-velocity CFT hypothesis. It is not a general central-charge measurement for a multi-velocity product. The raw coefficient and both fits are saved in `summary.json`. No c is assigned when the scaling assumptions fail.

### 5.3 Calibration and findings

[numerical] Calibration used exact XX and XXZ one-boson points and a separately constructed transverse-field Ising chain. An additive Ising times XX control has local dimension four and unequal velocities 2 and 4. Its entropies and energies can be added without constructing the enlarged many-body Hilbert space. It is a control, not an in-scope candidate.

| Control | Lmax | c_chord | c_size | Energy diagnostic |
| --- | ---: | ---: | ---: | --- |
| XX, true c=1 | 18 | 1.00610 | 0.98892 | W=3.99944, c_E=0.99986 |
| XXZ Delta=-1/2, true c=1 | 18 | 1.00826 | 0.99126 | v_inf=2.59820, c_E=1.00003 |
| Ising, true c=1/2 | 18 | 0.50063 | 0.49863 | W=0.999991, divided by exact v=2 gives 0.499996 |
| Ising times XX, true c=3/2 | 18 | 1.50672 | 1.48755 | W=4.99943, approaching 2(1/2)+4(1)=5 |

The product control demonstrates why W must not be divided by a single velocity when testing c=3/2. Its value would give 2.5 when divided by the slower speed or 1.25 when divided by the faster speed, neither of which is its central charge.

[numerical] The unresolved-family points give the following finite-size diagnostics. c_E has the conditional meaning just specified. The table does not identify an asymptotic central charge at every row.

| (u,b) | Lmax | c_chord | c_size | c_E assuming one velocity | gap(q=2)/gap(q=1) at Lmax |
| --- | ---: | ---: | ---: | ---: | ---: |
| (0.75,0.5) | 18 | 0.99377 | 0.97198 | 0.92704 | 3.65126 |
| (0.9,0.5) | 20 | 1.02458 | 1.02540 | 1.00391 | 3.95562 |
| (1.3,2) | 18 | 0.67078 | 0.50173 | 0.41016 | 2.63237 |
| (1.35,2) | 20 | 0.99390 | 0.99964 | 0.91957 | 3.26623 |
| (1.4,2) | 20 | 1.11421 | 1.24686 | 0.75143 | 3.90777 |
| (2.2,4) | 18 | 0.50246 | 0.27258 | 0.28324 | 2.35086 |
| (2.22,4) | 20 | 0.95522 | 1.02778 | 0.85397 | 2.87216 |

For a one-boson charge tower,

    E(q)-E0 = (2pi v/L) q^2/(4K)+o(1/L),
    K_L=pi v_L/[2L(E(1)-E0)].                            (8)

Thus the ratio in the last column tends to four in a regular liquid. At (0.9,0.5), K_L=3.069 and the lowest neutral state at the ground state's T2 momentum has energy 1.905 times the oscillator unit 2pi v_L/L. These are consistent with one boson with K>2, whose lowest same-momentum neutral descendant has dimension two. Together with both entropy fits and the energy coefficient, this is the cleanest tested one-boson point.

At (1.35,2), K_L=1.946 and the same-momentum neutral ratio is 1.193. At (2.22,4), these numbers are 2.078 and 1.028. Charge ratios and neutral spectra still have substantial corrections in this part of the apparent BKT region. In particular, a low neutral state with a ratio near one is not, by itself, an Ising energy primary. A whole independently identified tower, locality assignments, and its velocity would be required. No such extra tower was identified here. These two points do not have a c>1 entropy signal.

At (1.4,2), r-u=0.01421356. The entropy profile estimate falls from 1.20945 at L=10 to 1.11421 at L=20. The oscillator-speed estimate falls from 0.86852 to 0.54567, and c_size and c_E do not agree. The charge ratio approaches four, while K_L is large, 8.668 at L=20, and the same-momentum neutral ratio is 1.870. This is consistent with a one-boson crossover near the ferromagnetic boundary; it is not sufficient to extract its asymptotic c from these rings. The exact boundary calculation below explains why a large entropy coefficient is expected there.

At (1.3,2) and (2.2,4), the entropy estimates decrease markedly with size, and charge gaps flatten. This is finite-size evidence for a massive regime, not a proof of a gap. The pilot points (1.2,2) and (2.1,4) lie farther into the same apparent massive side and behave similarly. The weak-field point (0.75,0.5) is close enough to the expected BKT region that these sizes do not distinguish a very long correlation length from a liquid. None of these distinctions affects the negative result of the c>1 test.

[numerical, checks] Independent dense Pauli construction at L=6 agrees with the sparse charge blocks to 1.8e-15, and its full spectrum agrees across the alternate-site rotation to 1.8e-14. The exact upper-wing gap is reproduced to 1.2e-15. The XX energy agrees with -4/sin(pi/L), and the Ising energy with -2/sin(pi/(2L)); the latter check stays within 2.9e-14 through L=18. The largest recorded eigenvector residual in the scan, extensions, and pilot is 8.5e-10. These are numerical checks, not new thermodynamic proofs.

[heuristic, scope of the test] The strongest apparent excess tested here is (1.4,2), and it fails consistency between entropy, velocity, and energy diagnostics. The larger-field near-BKT point (2.22,4) gives approximately one in the entropy fits, with unresolved spectral corrections. This rejects neither a narrowly tuned additional point nor an arbitrarily long crossover. It supplies no candidate that can honestly be labeled c>1.

## 6. An exact c_size=3/2 false positive at the boundary

This section adds an explicit Schmidt-spectrum calculation to the inherited factorization-line analysis. It explains the most dangerous signal in Section 5 without assuming a phase diagram.

### 6.1 State and Schmidt spectrum

[proved] At u=r=sqrt(1+b^2/4), split the original field evenly between bonds. The shifted bond is

    h+r I=4r |v><v|,
    |v>=A|00>+B|11>,
    A^2=(2r+b)/(4r), B^2=(2r-b)/(4r), AB=1/(2r).          (9)

The odd-spin block is zero, and the even-spin block has eigenvalues zero and 4r. After alternate-site X rotations, the forbidden vectors alternate between A|01>+B|10> and A|10>+B|01>.

Let s=A/B=r+b/2, w=s^2, and let n_e count occupied even sites in the rotated computational basis. In each fixed-N sector,

    |psi_N> = Z_L(N)^(-1/2)
                sum_{|n|=N} (-s)^(n_e) |n>,
    Z_L(N)=[t^N](1+w t)^(L/2)(1+t)^(L/2).               (10)

On each bond the two amplitudes paired by its forbidden bra have ratio -A/B or its inverse, and cancel. Equation (9) is positive, so (10) is an exact ground state with energy -rL. Its sign is a product of on-site phases and does not affect the Schmidt eigenvalues. The inherited kernel argument shows that these L+1 states span the entire ground space.

For a block A with l sites, l_e even sites and l_o odd sites, define

    Z_A(n)=[t^n](1+w t)^(l_e)(1+t)^(l_o),
    Z_B(N-n)=[t^(N-n)](1+w t)^(L/2-l_e)(1+t)^(L/2-l_o).

At fixed block occupancy n the amplitudes factor into one vector on A and one on its complement. Different n are orthogonal. The nonzero Schmidt eigenvalues are therefore

    lambda_n=Z_A(n) Z_B(N-n)/Z_L(N).                     (11)

Normalization follows by multiplying the two generating polynomials. No tensor-network approximation or large-Hilbert-space calculation is needed.

### 6.2 Entropy asymptotics

[proved, standard saddle expansion of (11)] Take N=L/2, fixed finite b, and l=xL with x in a compact subinterval of (0,1). A fugacity z0=1/s centers the total occupancy at L/2. Multiplying the coefficients by z0^n changes neither (11) nor its normalization. The two normalized factors then have occupancy parameters

    q=s/(1+s), 1-q=1/(1+s),
    a_b=q(1-q)=1/[2(1+r)].                               (12)

The log of the coefficient-generating function has second derivative L a_b at its saddle. For a block and its complement the corresponding second derivatives are l a_b and (L-l)a_b, up to bounded parity corrections. Multiplying their Gaussian coefficient asymptotics with the total N fixed gives

    variance(lambda)=a_b l(L-l)/L+O(1).

The discrete Gaussian entropy is one half the logarithm of 2pi e times this variance. Thus

    S(l)= (1/2) log[2pi e a_b l(L-l)/L]+o(1)
        = (1/2) log L+(1/2) log[x(1-x)]
            +(1/2) log(2pi e a_b)+o(1).                 (13)

The saddle is nondegenerate for every fixed finite b. This statement is not uniform in b tending to infinity with L.

A periodic-CFT size fit would give

    c_size=3 dS(L/2)/d log L -> 3/2.                     (14)

But the block shape in (13) is log[x(1-x)], not log sin(pi x). Fitting that shape over x in [1/4,1/2] gives c_chord about 1.241, not 3/2. The mismatch persists at arbitrarily large L.

### 6.3 Spectrum and finite checks

[proved, inherited exact defect branch] The same Hamiltonian has a closed one-defect subspace above either original Neel state. Its two-site-cell momentum k gives

    E_-(k)=4r-4sqrt[r^2-sin^2(k/2)]
           =k^2/(2r)+O(k^4).                            (15)

For k=4pi/L this is 8pi^2/(r L^2)+O(L^-4). The would-be speed L E_-(4pi/L)/(2pi) tends to zero as 4pi/(rL). Together with the L+1 ground states, this rules out the ordinary finite-velocity CFT interpretation. The branch is an exact excitation in the appropriate charge sector; it is not an assertion that every low state above the N=L/2 vector has this dispersion.

[numerical] `boundary_entropy.py` evaluates (11) through L=512. These are generating-polynomial calculations, not large ED. At L=8 it also constructs (10), verifies its ground-state residual, and compares the Schmidt formula against a direct reduced density matrix. The largest residual is 6.1e-15 and the largest entropy discrepancy is 4.5e-16.

| b | L | c_size from L/2 to L | c_chord | Exact defect speed estimator |
| --- | ---: | ---: | ---: | ---: |
| 2 | 32 | 1.46751 | 1.24443 | 0.27545 |
| 2 | 128 | 1.49043 | 1.24149 | 0.06939 |
| 2 | 512 | 1.49754 | 1.24092 | 0.01735 |
| 4 | 32 | 1.56330 | 1.24543 | 0.17371 |
| 4 | 128 | 1.50581 | 1.24145 | 0.04388 |
| 4 | 512 | 1.50110 | 1.24092 | 0.01098 |

[proved, interpretation] The half-charge ground state has logarithmic entanglement from a constrained collective degree of freedom in a degenerate ferromagnetic-type manifold. One can also choose product ground states at the same Hamiltonian. Its c_size=3/2 is not an invariant conformal central charge. This supplies an explicit failure mode for exactly the numerical signature that could be mistaken for Ising times Luttinger in this family.

## 7. What would finish the argument

[proved, remaining obligations] The first four settled families in Section 1 and their controlled neighborhoods do not exhaust the unitary-Z2, antiunitary, or symmetry-free branches. Even an exhaustive solution of A_st^o would close only the exact additive-spin-charge route. The inherited numerical junction and two-fluid studies do not give open neighborhoods covering the other branches. C34 explicitly records this limitation.

A sufficient additional theorem would say that every conformal point outside the known gapped and singular sets admits a complete continuum matching with at most one compact boson or at most two Majoranas, with all other thresholds uniformly massive in its matching regime. Given a Lorentz-invariant matching, the c-theorem would then give the desired bound; given a weak multi-velocity matching, its separate controlled RG argument would be needed. This statement is [assumed] if used globally. It contains the missing dynamical content and is not implied by two microscopic Majoranas, a single microscopic charge, or anomaly matching.

[heuristic, focused next tests] Within A_st^o, the charged-sector result (3) makes the next question concrete: can a neutral Ising stress tensor become critical while the charged oscillator retains finite speed? The useful measurements are charge-resolved large-system spectra or symmetry-resolved tensor-network data, a converged bulk entropy slope, and an independently resolved neutral tower. The region around (u,b)=(1.35,2) and (2.22,4) merits such a boundary study, but neither point currently has evidence for c>1. Small-ring ratios near one need not be an Ising primary.

For the full family, one must additionally address irreducible c>1 theories and neutral sectors on the TR-only, reflection-only, and no-symmetry rows. A finite list of Ising and Gaussian junctions does not do that. The supplied points 220 and 1282 are inherited numerical c=1 floating phases, not counterexamples; this run did not rerun their DMRG.

The result is a status-tagged covering with a concrete staggered-family test, not the requested global exclusion. Two-sector coincidence is allowed at c=1. No c>1 candidate survives the present checks, and no claim is made that the unresolved rows are empty.

## 8. Files, reproduction, and provenance

`PROGRESS.md` records the working status. `STATE.md` is the local memory for this work directory; the project-level memory is outside the authorized write scope and was not changed.

New maintained scripts and outputs are

| File | Purpose |
| --- | --- |
| `ed_staggered.py` | Charge-resolved sparse ED, translation labels, entropy, checks |
| `scan.jsonl` | Nine points at six sizes, L=8 through 18 |
| `extended.jsonl` | Four selected points at L=20 |
| `pilot.jsonl` | Initial small-ring probes, including two farther-gapped-side points |
| `extra_checks.py` | All-charge-sector minima at L=18 and Ising/product controls |
| `all_sector_checks.json`, `controls.json` | Results of those checks |
| `analyze.py`, `summary.json` | Entropy-size, energy, velocity, charge, and neutral-gap diagnostics |
| `boundary_entropy.py` | Exact filtered-state Schmidt spectrum and direct state checks |
| `boundary_entropy.json`, `boundary_checks.json` | Large-L coefficient results and small-L convention checks |

Run each calculation from this directory. The ED script appends to its chosen output, so use fresh output names when repeating a scan. `analyze.py` deduplicates by point and size. Limit BLAS threads to one for these small matrices.

```sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 python3 ed_staggered.py --out repeat_scan.jsonl --L 8 10 12 14 16 18 --points 0,0 0.5,0 0.75,0.5 0.9,0.5 1.3,2 1.35,2 1.4,2 2.2,4 2.22,4 --k 14
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 python3 ed_staggered.py --out repeat_extended.jsonl --L 20 --points 0.9,0.5 1.35,2 1.4,2 2.22,4 --k 14
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 python3 extra_checks.py
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 python3 analyze.py
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 python3 boundary_entropy.py
```

`analyze.py` reads the recorded `scan.jsonl` and `extended.jsonl`; to analyze a repeat, change those input names or call its `summarize` function on the repeated rows. The largest ED Hilbert block was binomial(20,10)=184756. No DMRG, network, cluster job, or large local many-body calculation was used.

The principal local sources were `../../STATUS.md` (d=2 entries), `../work_C6_d2/c6.tex`, `../work_C15/RESULT.md`, `../work_C16/RESULT.md`, `../work_C17/RESULT.md`, `../work_C34/RESULT.md`, the README and companion algebra/gap notes in `../../contrib/d2_gap_reduction/`, and the README summaries in `d2_classification`, `d2_multicritical`, `d2_two_fluid`, and `d2_followup_220_1282`. The project memory was read before these sources. Existing constants, spectra, and exclusions were cited rather than claimed as new derivations. C6 supplies the invariant strata and free theorem; C15 the controlled perturbation theory and factorization boundary; C16 the covering and operator-count limitations; C17 the dual two-Ising example and anomaly limitations; C34 the additional gap certificate and residual branches. The new calculations are the staggered-family ED data and the explicit filtered-state entropy diagnostic.

