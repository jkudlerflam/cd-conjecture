# C67. The nearest-neighbor qubit stratum with one unitary Z2

2026-10-28

The stratum can be parametrized completely, and several of its open regions and critical neighborhoods are controlled. The remaining global step is not closed. The new small-ring calculations do not identify a supported c > 1 point. They test four flagged grid candidates and an interacting Ising/chiral-liquid junction. The largest chord-entropy signal, c_chord = 2.259 at L = 20, disagrees with entropy growth across sizes. The junction instead exhibits a collapsing descendant velocity.

[assumed, scope] Central charge is the usual nonchiral value, with c = 1/2 for Ising. A conformal interpretation requires nonzero limiting mode velocities. The exact-model continuum identifications, unitary compact CFT framework, and regular weak-coupling matching are inherited from C63 and C6. This is an offline analysis of the supplied sources, not a new literature survey. No independent referee has checked this calculation.

The labels [proved], [numerical], [heuristic], and [assumed] specify the status of the following claims. Numerical estimates below are finite-size diagnostics, not error-controlled thermodynamic central charges.

## 1. Complete parametrization and the meaning of “only Z2”

[proved, inherited algebra with an explicit strict-symmetry refinement] Remove a constant and telescoping bond fields. A uniform on-site pi rotation can be put along z. Its action on Pauli vectors is R = diag(-1,-1,1). The conditions R J R^T = J and R B = B force

    J = [[a, q+D, 0],
         [q-D, c, 0],
         [0,   0, u]],       B = (0,0,b).

A uniform rotation in the xy plane diagonalizes the real symmetric transverse block. Every Hamiltonian in the row therefore has the representative

    H = sum_j [Jx X_j X_(j+1) + Jy Y_j Y_(j+1)
               + u Z_j Z_(j+1) + D(X_j Y_(j+1)-Y_j X_(j+1))
               + b Z_j].                                            (1)

Its parity is P = product_j Z_j. A symmetric bond includes b(Z_j+Z_(j+1))/2. There are five real parameters before positive scale, hence four after scale. One possible normalization is

    Jx^2 + Jy^2 + u^2 + 2D^2 + b^2 = 1.                              (2)

The discrete identifications do not change this count. Interchanging x and y by a proper rotation changes the sign of Jx-Jy. A uniform pi rotation about x reverses D and b together. Spatial reflection reverses D, and time reversal relates b to -b. The latter two are equivalences between parameter points, not generally internal symmetries of a point. An alternate-site z rotation reverses Jx,Jy,D on even rings. Negative overall scale is not an equivalence of ground-state problems.

[proved] Exactly one nonidentity uniform proper rotation survives if and only if

    Jx != Jy,           (D,b) != (0,0).                              (3)

If Jx=Jy, axial U(1) is restored. If D=b=0, the other principal-axis pi rotations survive. These conclusions do not require u to differ from Jx or Jy.

There is a further distinction if “only symmetry” excludes uniform internal antiunitaries as well. Their Pauli action O obeys

    det O=-1,   O S O^T=S,   O D_vector=-D_vector,   O B=B.

When D and b are both nonzero and parallel to z, the last two equations would require O z=-z and O z=z simultaneously. No such O exists. If D=0 and b!=0, reflections of x or y survive. If b=0 and D!=0, time reversal and reflection of z survive. Thus the strict internal-symmetry stratum is

    Jx != Jy,                   D b != 0.                           (4)

The broader C63 unitary-C2 row also includes those two enhanced faces. They are retained below as useful limits. Spatially combined antiunitaries still exist in (4). For example, spatial reflection times complex conjugation in the Z basis is a symmetry. “Only Z2” here refers to uniform internal operations, not every spatial or nonlocal symmetry.

[proved, inherited charge classification with a direct check here] The strict portion (4) has no nonzero additive spin charge. Commuting an additive Q=sum_j n_j.sigma_j with the nonzero field forces n_j=q_j z. In the fermion representation below, its hopping coefficient is nonzero because D!=0, and charge conservation requires q_(j+1)=q_j. Its pairing coefficient is nonzero because Jx!=Jy, and charge conservation instead requires q_(j+1)=-q_j. Hence every q_j vanishes. With b!=0, the only extra charge class outside uniform U(1) is Jx=-Jy, D=0, where sum_j (-1)^j Z_j is conserved. This lies on the enhanced face, not in (4).

## 2. Fermionic representation and what it does not prove

[proved] Set t=Jx+Jy and delta=Jx-Jy, and use n_j=(1-Z_j)/2. Up to the parity-dependent boundary condition, (1) becomes

    H = sum_j [(t-2iD)c_j^dagger c_(j+1) + h.c.
               + delta(c_j^dagger c_(j+1)^dagger+c_(j+1)c_j)
               + 4u(n_j-1/2)(n_(j+1)-1/2)
               - 2b(n_j-1/2)].                                     (5)

This is one spinless band with pairing and nearest-neighbor density interaction. The transformation is exact. The statement that every interacting infrared limit has only one propagating channel does not follow from it. Eliminating high-energy states generates longer-range effective hopping and interactions. The number of bare fermion operators is not a bound on emergent conformal sectors.

[proved conditional on C63's CFT input] If a nontrivial U(1) current actually exists in the continuum, its Sugawara tensor has c=1. The commuting neutral tensor has c_neutral=c-1. Unitary Virasoro representations with 0<c_neutral<1 begin at c_neutral=1/2. Hence a charged compact theory cannot have 1<c<3/2. This applies to the staggered-charge face when its charge remains gapless. It does not apply to a generic point of (4) merely because (5) uses fermions. Fermion parity supplies no U(1) current. Even with an emergent current, its existence and the commuting-sector hypotheses have to be established first.

## 3. Exactly solvable loci and controlled neighborhoods

### 3.1 The entire free face

[proved, inherited free-mode theorem] At u=0, the quadratic eigenvalues are, up to reversing k,

    lambda_+/-(k) = 4D sin k
        +/- 2 sqrt[(t cos k-b)^2+delta^2 sin^2 k].                    (6)

The range-one Majorana theorem supplied in C6 gives c<=1 for every conformal point on this face. The following direct analysis also locates the singular junctions. Put s^2=4D^2-delta^2.

If s^2<0, the zero condition is a sum of nonnegative terms,

    (t cos k-b)^2 + (delta^2-4D^2) sin^2 k = 0.

For t!=0 it gives a single Ising point at b=+/-t, with speeds

    v_R,L = 2|delta| +/- 4D,

after orienting the momentum. Both must be positive in this labeling. The exceptional t=b=0 point has two Majorana sectors, totaling c=1.

If s^2>0, the zeros satisfy

    (t x-b)^2=s^2(1-x^2),  x=cos k,
    x=[tb +/- s sqrt(t^2+s^2-b^2)]/(t^2+s^2).                       (7)

The regular liquid has one quadratic-fermion channel and c=1. BdG particle-hole partners must not be counted as independent fluids. Its outer boundary |b|=sqrt(t^2+s^2) is a pocket closure with a quadratic soft edge, except at further degenerate limits.

[proved] At |2D|=|delta| and b=t!=0, the apparent Ising/Gaussian meeting is not Ising times a boson. Take delta>0, D=delta/2 and k>0. Expanding (6),

    t(cos k-1)=-t k^2/2+O(k^4),
    sqrt[t^2(cos k-1)^2+delta^2 sin^2 k]
       =delta k+[-delta/6+t^2/(8delta)] k^3+O(k^5).

The slow branch therefore has

    |lambda_slow(k)| = t^2 k^3/(4delta)+O(k^5),                     (8)

while the other branch is linear with speed 4delta. At |b|<|t| on the same tilt boundary, the generic touching at cos k=b/t is quadratic instead. Flat exceptions occur when further coefficients vanish. None of these is a regular finite-velocity c=3/2 theory.

[assumed, controlled continuum matching] Near a regular free Ising point, retain its Majorana fields and require all discarded momentum regions to stay separated by a positive matching scale. A small u shifts the thermal mass and velocities. The local derivative-free quartic of a single Majorana vanishes after normal ordering. Remaining irrelevant operators cannot create an extra low-energy sector while the matching separation and both velocities remain positive. The Z2 symmetry forbids a local odd magnetic source. Retuning the mass gives an Ising critical surface in this controlled regime.

Near a regular liquid, the retained fields form one compact boson. In a convention with dimensions dim cos(2phi)=K and dim cos(2theta)=1/K, forward scattering changes K and the speeds. Pairing, umklapp, and their commensurate versions pin or perturb these same fields. An incommensurate vertex averages away when its momentum mismatch remains finite on the matching scale. At a two-Majorana point the interspecies energy interaction has dimension two and produces the familiar Gaussian interaction. These are local descriptions with a specified domain, not an infrared classification at arbitrary u.

[proved conditional on the continuum assumptions] When the complete matching theory and subsequent flow are unitary and Lorentz invariant, the c-theorem bounds c_IR by the matched c<=1. A common drift can be removed. Unrelated velocities require a separate weak-coupling argument. In particular, the small parameter includes interaction strength relative to the smallest velocity. The above neighborhood shrinks as that velocity tends to zero. It does not cover (8) by continuity.

### 3.2 Enhanced exact-model faces

[proved mappings; assumed exact-model infrared input] At delta=0, a spiral rotation gauges the hopping phase into a boundary twist. The bulk is uniform-field XXZ with transverse exchange sqrt(Jx^2+D^2). Its complete standard ground-state solution has at most one massless Bethe sea and c=1. Boundary twists can shift the finite-size vacuum energy without changing c. Small symmetry-breaking anisotropy near a regular XXZ liquid has the one-boson description just stated, subject to momentum selection and a nonzero matching scale.

At D=b=0 the XYZ gap classification reduces every regular critical locus to a Gaussian c=1 theory, including the antiferromagnetic SU(2)_1 endpoint. Its gapped interiors have gapped neighborhoods under the standard stability input. Saturation edges and the ferromagnetic endpoint are separate nonrelativistic limits.

At Jx=-Jy=a and D=0, alternate-site X rotations on even rings give

    H' = sum_j [a(XX+YY)-u ZZ+(-1)^j b Z].                         (9)

The gap reduction and C63 already establish the gapped wings u<=-|a| and u>r, r=sqrt(a^2+b^2/4), and the L+1-fold ground manifold with a quadratic branch at u=r. For a=1, b=0, -1<u<1, the inherited stiffness is K=pi/[2 arccos u]. A small staggered field couples to cos(2phi) of dimension K, hence is relevant for K<2 and irrelevant for K>2. This controls weak-field regimes, not all finite b. C63's exact boundary entropy has c_size tending to 3/2 despite nonconformal dynamics. Those calculations are inherited, not rerun here.

### 3.3 A quantitative gapped cone inside the strict stratum

[proved using the imported gap theorem] In (1) the transverse exchange matrix is

    M = [[Jx,D],[-D,Jy]].

Its squared singular values follow from the trace and determinant of MM^T. Writing a=(Jx+Jy)/2 and d=(Jx-Jy)/2 gives singular values

    |sqrt(a^2+D^2)+|d||,       |sqrt(a^2+D^2)-|d||.

Consequently the gap criterion from the supplied gap-reduction README becomes

    u <= -sqrt((Jx+Jy)^2/4+D^2)-|Jx-Jy|/2,
    b != 0   =>   gap_L >= 2|b|.                                  (10)

The theorem is Bravyi-Gosset-Liu-Wong as quoted in that README. It was not independently fetched here. Strict inequality in (10) gives an open gapped region within the four-dimensional stratum. For example (Jx,Jy,u,D,b)=(-1,-0.5,-2,0.2,0.3) has only internal Z2 and gap at least 0.6. Large fields also have product-phase neighborhoods under ordinary gapped-phase stability, but no additional numerical radius is claimed for them.

### 3.4 Factorized ordered points and their neighborhood

[proved, elementary lattice construction] Another useful exact family is

    Jx=-A, Jy=-B, u=C,
    A>B>0, C>-B,
    b=-h_f,  h_f=2 sqrt[(A+C)(B+C)],
    cos theta=sqrt[(B+C)/(A+C)].                                  (11)

At D=0 split the field evenly over bonds and subtract e_*=-A-B-C. The odd-spin block of h_bond-e_*I is

    [[A+B,-(A+B)],[-(A+B),A+B]],

with eigenvalues 0 and 2(A+B). The even-spin block is

    [[A+B+2C-h_f,-(A-B)],[-(A-B),A+B+2C+h_f]].

Its determinant vanishes because

    (A+B+2C)^2-h_f^2-(A-B)^2=0,

and its eigenvalues are 0 and 2(A+B+2C), both nonnegative. Its full kernel is spanned by |+ +> and |- ->, where

    |+> = cos(theta/2)|0>+sin(theta/2)|1>,
    |-> = cos(theta/2)|0>-sin(theta/2)|1>.

These vectors are independent since 0<cos theta<1. An invertible single-site map taking the computational basis to |+>,|-> maps the common bond-kernel constraints to “all neighboring labels equal.” Thus the periodic ground space at D=0 is exactly the span of |+>^L and |->^L, with energy L e_*.

The uniform DM sum annihilates each homogeneous product state. For example, choose the orthogonal |perp> to |+> so that

    X|+>=sin theta |+>+cos theta |perp>,
    Y|+>=i|perp>.

Then

    (XY-YX)|+ +>
       =i sin theta (|+ perp>-|perp +>).

Summing over a ring cancels each one-site defect against its neighbor. The same cancellation holds for |->^L. The two states in (11) therefore remain exact eigenstates at energy L e_* for every D. This statement alone does not say they remain ground states for arbitrarily large D.

[assumed, standard gap and stability inputs] At D=0 this is the parent Hamiltonian of two distinct product-state sectors. The standard one-dimensional parent-Hamiltonian gap result supplies a size-independent gap at fixed 0<cos theta<1. The overlap of the two block states is (cos theta)^ell, so the distinction between sectors becomes exponentially sharp under blocking; the singular endpoint cos theta=1 is excluded. Stability of this ordered gapped phase under sufficiently small parity-preserving local perturbations supplies a nonzero neighborhood. Inside it the two exact eigenstates above remain the ground sector. Thus sufficiently small nonzero D in (11) gives controlled factorized ground states in the strict stratum (4). No explicit thermodynamic upper bound on this allowed D is derived here. Nearby generic parity-preserving perturbations remain gapped but need not remain factorized.

[numerical, convention checks] At A=1, B=0.5, C=0.25, h_f=1.9364916731, the L=6 product-state residual is at most 2.1e-15 for D=0,0.1,0.2. Full diagonalization gives two ground states at -10.5 in all three cases, and gaps above them 1.143808, 1.143808, and 0.987155. These finite checks do not set the thermodynamic stability radius.

## 4. Why these neighborhoods do not yet cover the row

[proved mapping; heuristic consequence] On the D=0 face, Kramers-Wannier variables obey

    XX -> tau_j^z,
    YY -> -tau_(j-1)^x tau_j^z tau_(j+1)^x,
    ZZ -> tau_(j-1)^x tau_(j+1)^x,
    Z  -> tau_(j-1)^x tau_j^x.

For Jy=0 this is a quantum ANNNI-type chain, with transverse field Jx, nearest-neighbor coupling b and next-nearest-neighbor coupling u. Its frustration and possible floating regimes are real dynamical questions, despite the original spin bond having range one. The nonzero-D term maps to a local dual chirality term as in C6. This observation supplies a useful boundary model, not an exact solution of the strict stratum.

[heuristic] A continuous transition of a discrete order parameter inside a pre-existing liquid is the most direct c=3/2 mechanism. But one microscopic Z2 need not protect an independent Ising sector. The putative Ising spin may be a vertex of the same boson, or the “new” branch may be a velocity-zero threshold. These alternatives require spectral and operator assignments.

[proved conditional operator statement] In Ising times a boson, epsilon has dimension one and sigma has dimension 1/8. If sigma is odd under the microscopic Z2, its magnetic source is forbidden while its thermal source remains allowed. If a boson vertex V is also odd, sigma V can be invariant and has dimension 1/8+dim V. The product epsilon V has dimension 1+dim V. Symmetry alone does not determine whether these products are local or momentum zero in the lattice embedding. Parameter counting cannot exclude simultaneous zeros of their microscopic coefficients in a four-parameter family.

[proved limitation] Even excluding that product mechanism would not exclude irreducible c>1 theories without a U(1) current. C63 explicitly retains that obligation. Neither the connectedness of (1), the number of lattice Majoranas, nor the absence of a projective-cell obstruction proves a global central-charge bound. A path through bare Hamiltonians is not an RG trajectory. The neighborhoods in Section 3 have no established finite subcover of the unresolved region.

## 5. Small exact diagonalization

### 5.1 Method, controls, and reproducibility

[numerical] `ed_core.py` copies the corrected momentum/parity kernel from C6. `ed_z2.py` computes four low levels in every momentum m=0,...,L-1 and both parities P=(-1)^p. Translation has eigenvalue exp(2pi i m/L). The lowest state across all blocks is expanded in the full spin basis for entropy. The largest L is 20 and the largest block has dimension 26272. There is no large-system local computation or tensor-network run in C67.

For each ground state the chord fit is

    S_L(ell)=s0+(c_chord/3) log[(L/pi)sin(pi ell/L)],
    max(3,floor(L/4)) <= ell <= L/2.                              (12)

The size fit is three times the slope of S_L(L/2) against log L. Separate L mod 4 sequences expose sector changes. Such fits are only diagnostics when the available sizes belong to different finite-size sectors. Eigenvectors within a degenerate ground space can have different entropies; the parity/momentum choice is recorded. No single velocity is used to turn a Casimir coefficient into c.

[numerical, independent checks] `analytic_checks.py` constructs a dense Pauli Hamiltonian at L=6 independently of the orbit kernel. At a generic strict-Z2 point its complete spectrum agrees with the union of orbit blocks to 6.3e-15, and the basis-intertwining residual is below 9.5e-15. The Ising and XX ground energies agree with -2/sin(pi/(2L)) and -4/sin(pi/L), respectively, to 4.7e-14 through L=16.

Across all 70 point-size calculations, the largest retained eigenvector residual is 9.3e-10. Full low spectra are saved, but no exhaustive conformal-tower identification was performed for the four grid points. Their entropy inconsistencies already prevent the proposed high-c interpretation of these sizes.

| Control | c_chord at L=16 | c_size, L=10,12,14,16 | Interpretation |
| :--- | ---: | ---: | :--- |
| Critical Ising | 0.500721 | 0.498124 | Known c=1/2 |
| XX | 1.007044 | 0.984833 | Known c=1 |
| Free tilt junction, (-1,-0.5,0,0.25,-1.5) | 0.507049 | 0.552948 | Cubic slow branch from (8), no ordinary finite-velocity c assignment |

The third control shows why an apparently good Ising entropy fit alone does not certify regular Ising dynamics at the velocity-collapse endpoint.

[numerical, spectral consistency check] In that free junction control, L^2 times the slow descendant estimate (15) is 37.0902,38.9197,40.1642,41.0417 at L=10,12,14,16. Equation (8) predicts the limit 4.5 pi^2=44.4132. The other descendant estimate decreases from 2.24188 to 2.10931, approaching the predicted nonzero speed 2. Thus the spectra distinguish the cubic branch from a second regular conformal sector even while the entropy fit stays near 1/2.

### 5.2 Four candidates from the in-progress grid

[numerical, selection and provenance] The grid's `candidates.json` contained four Z2 entries when read. All four were copied to `grid_candidates.json` and recomputed at L=10,12,14,16,18,20. The tuples below are ordered (Jx,Jy,u,D,b). Full precision is retained in the JSON files.

| Label | Grid name | Couplings |
| :--- | :--- | :--- |
| A | Z2_1_358 | (-0.660398, 0.523541, 0.234266, 0.482843, 0.041960) |
| B | Z2_21_295 | (0.575442, -0.618578, 0.243667, -0.460508, 0.121600) |
| C | Z2_22_54_p01 | (-0.150676, -0.710538, -0.427524, 0.500794, 0.197132) |
| D | Z2_11_128_p05 | (-0.397435, 0.372111, -0.390416, 0.700001, -0.247291) |

All four obey (4) and are interacting. The L=12 and L=16 chord fits at A and B reproduce the inherited grid values within 1.1e-11. Thus the large slopes are not a mismatch between the two ED implementations.

| Point | c_chord(12) | c_chord(16) | c_chord(20) | c_size, L=12,16,20 | c_size, L=10,14,18 |
| :--- | ---: | ---: | ---: | ---: | ---: |
| A | 2.268926 | 1.964341 | 1.745468 | 0.364482 | 2.189231 |
| B | 2.406142 | 2.350631 | 2.259132 | 0.960061 | 0.867243 |
| C | 1.185015 | 1.168927 | 1.161689 | 0.642899 | 2.359829 |
| D | 1.281890 | 1.268930 | 1.215727 | 0.246193 | -0.001871 |

[numerical] B is the largest surviving fixed-size signal, but both separate size-growth fits are near one. Its half-chain entropies are 1.626334, 1.710860, 1.790472 at L=12,16,20, and 1.351125, 1.447890, 1.521094 at L=10,14,18. Its parity is even throughout, while momentum changes. The discrepancy with c_chord above two rules out interpreting these rings as an already resolved c>1 CFT. It does not prove the thermodynamic bound.

At A the two size sequences disagree sharply. The L=14 to 18 half entropy jumps from 1.292146 to 1.616495, while the L=12,16,20 sequence grows slowly. The full L>=12 fit happens to give 1.015701, but averaging incompatible sector behavior is not a controlled extrapolation. A remains a crossover target rather than a measured high-c point.

At C the L=10 ground state has momentum zero and S_half=0.768287. From L=12 onward the ground momentum is m=L-2 and S_half grows from 1.128977 to 1.238233. Restricting to L>=12 gives c_size=0.641414. The apparent 2.359829 in the last column is caused by including the L=10 branch change. The momentum approaches zero as -4pi/L, a useful warning about a finite-particle threshold. This is not a determination of the dynamical exponent.

At D the L=10,14,18 entropies are 0.498695, 0.498088, 0.498369, with even parity and momentum zero. The L=12,16,20 ground states instead have odd parity and m=L/4, and entropies 1.088700, 1.114991, 1.130385. A near-crossing can produce a large entropy offset between sequences without logarithmic growth within either sequence. Its L>=12 combined slope is only 0.014075. No high-c interpretation is supported.

[numerical, inherited grid cross-check] The available one-site-cell finite-entanglement runs for A and B have residuals of order 0.09 to 0.19, enormous transfer lengths, and some energy densities of order 10^13 to 10^15. They fail the grid's own residual filter and cannot be used as evidence. A separate 11-site-cell run for A has chi=8,12,16,24, residuals below 2.4e-6 and correlation lengths 4.39,8.42,11.65,20.94. Its four-point fit 6 dS/dlog xi is 1.020698, with consecutive slopes 0.971224,1.178399,0.970080. This is consistent with a one-boson interpretation, but the small chi range and cell choice prevent a final thermodynamic conclusion. C67 did not rerun that tensor-network calculation. No completed FES output for the two patch points was present in the inspected snapshot.

### 5.3 The interacting Ising/chiral-liquid junction

[numerical] The targeted family is

    (Jx,Jy,u,D,b)=(-1,-0.5,0.25,D,-g).                             (13)

C6 already found Ising behavior near D=0.1, g=1.9891. To follow the massless branch, at each L and D solve the finite-size spectral condition

    E_(even,m=0,1)-E_(even,m=0,0)
       =8 [E_(odd,m=0,0)-E_(even,m=0,0)].                        (14)

The factor eight is the Ising ratio x_epsilon/x_sigma=1/(1/8). This is a phenomenological tuning condition, not an independent test of that same ratio. All momenta and parities are then checked. The actual global ground state remains even and m=0 at the tuned points through L=20.

Define the two odd-sector descendant differences

    v_+(L)=L[E_(odd,1,0)-E_(odd,0,0)]/(2pi),
    v_-(L)=L[E_(odd,L-1,0)-E_(odd,0,0)]/(2pi).                   (15)

They are velocity estimators only while this tower assignment holds. A negative value diagnoses its breakdown.

| L | g_L at D=0.20 | c_chord | v_+, D=0.16 | v_+, D=0.20 | v_+, D=0.24 | v_-, D=0.20 |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 12 | 1.99930056 | 0.482044 | 0.476882 | 0.324108 | 0.171335 | 1.851815 |
| 14 | 1.99547910 | 0.493877 | 0.403611 | 0.248945 | 0.094279 | 1.795575 |
| 16 | 1.99326193 | 0.499494 | 0.352733 | 0.196832 | 0.040930 | 1.755821 |
| 18 | 1.99189745 | 0.502432 | 0.316123 | 0.159370 | 0.002618 | 1.726868 |
| 20 | 1.99101535 | 0.502816 | 0.288984 | 0.131622 | -0.025741 | 1.705225 |

The D-dependence of g_L and c_chord is very small at these sizes, but no exact D-independence is claimed. At D=0.24 and L=20 the m=1 odd state is below the m=0 odd state. The latter is no longer the bottom of an Ising primary family. Both remain above the even ground state at this size.

[heuristic informed by these numerics] Extrapolating v_+ over L=12,...,20 with a constant plus L^-2 gives a zero near D=0.2065. Adding L^-4 moves it to approximately D=0.2016. The difference is model dependence of the fit, not an error bar. The useful target is therefore the neighborhood g near 1.989, D near 0.20 in (13), not a claimed critical coordinate to four digits. Its observed mechanism resembles the free velocity-collapse junction (8). The calculation has not established its interacting exponent, nor identified an additional neutral Ising tower with finite speed.

The Hamiltonian changes with L in (14). A size-entropy slope along that tuned sequence is not used as a central-charge measurement.

[numerical, fixed-coupling follow-up] At the fixed field g=1.9891, D=0.16,0.20,0.24 give almost identical finite-size ground-state entropies through L=20. Their c_chord values at L=12,16,20 are 0.510020,0.509836,0.506764. The size fit over these three rings is 0.563481, and the L=16 to 20 step is 0.541663. These are consistent with the previously found Ising finite-size crossover, but the longer-scale instability at D=0.24 cannot be excluded by its still-Ising-like ground state on these rings.

At D=0.30 the ground state changes from m=0 at L=12 to m=2 at L=16 and 20, with even parity throughout. Its half entropies are 0.458659,0.977097,0.994146. The L=20 chord fit is 1.193648, whereas the L=16 to 20 size step is only 0.229217. Fitting all three sizes would give 3.251540 because it includes the branch change. This is another resolved finite-size ambiguity, not evidence for c>1. These fixed-coupling records are in `fixed_ed.jsonl`.

## 6. Verdict and remaining obligation

[proved, scope of the result] Equations (1) to (4) give the complete algebraic row and distinguish its strict part from the antiunitary-enhanced faces. The free family, the imported exact-model faces, the quantitative cone (10), and the factorized construction (11) provide controlled points and neighborhoods. Their dynamical claims rely only on the explicitly separated exact-model, continuum, or stability inputs.

[numerical] The four selected grid candidates do not supply a consistent c>1 measurement through L=20. B is the strongest fixed-size candidate and fails the size-growth check. A has conflicting sector sequences and an inherited short FES sequence consistent with c=1. The interacting junction has a decreasing and eventually negative candidate descendant difference, not evidence for an extra independent finite-speed sector.

[heuristic] The next useful junction test is to track the slow branch near (13), g near 1.989, D near 0.20, and distinguish its critical dispersion from a new neutral tower while the charged or Gaussian branch retains finite speed. This is a concrete test location, not a supported c>1 point. Larger converged calculations would also resolve the pronounced sector dependence at A and B. Broadening a small-ring entropy maximization would not address the identified ambiguity.

[assumed if used globally] A complete infrared matching with at most one boson or two Majoranas at every remaining conformal point would finish the bound. No such matching theorem or exhaustive phase diagram is established here. Additional soft bound states, strong-coupling multicriticality, and irreducible c>1 theories remain logically unexcluded. The evidence favors the known one-channel mechanisms in the tested regions, but does not justify stating c<=1 throughout this stratum, even with the claim labeled as a theoretical-physics argument.

## 7. Files and reproduction

All scripts run from this directory, with Python, NumPy, SciPy, and Numba. Output files append, so use fresh names for repeated ED runs. `analyze.py` deduplicates by point and size.

    OPENBLAS_NUM_THREADS=1 python3 ed_z2.py grid --sizes 10 12 14 16 18 20 --out repeat_grid.jsonl
    OPENBLAS_NUM_THREADS=1 python3 ed_z2.py junction --sizes 10 12 14 16 --Ds .1 .16 .20 .24 --out repeat_junction.jsonl
    OPENBLAS_NUM_THREADS=1 python3 ed_z2.py junction --sizes 18 20 --Ds .16 .20 .24 --out repeat_junction.jsonl
    OPENBLAS_NUM_THREADS=1 python3 ed_z2.py controls --sizes 10 12 14 16 --out repeat_controls.jsonl
    OPENBLAS_NUM_THREADS=1 python3 ed_z2.py fixed --sizes 12 16 20 --Ds .16 .20 .24 .30 --out repeat_fixed.jsonl
    OPENBLAS_NUM_THREADS=1 python3 analytic_checks.py
    OPENBLAS_NUM_THREADS=1 python3 analyze.py

The last two scripts read the original saved filenames, not the `repeat_` filenames. Their input paths must be changed to analyze a repeat. `summary.json` contains the processed diagnostics. `analytic_checks.json` records independent finite checks. `grid_candidates.json` freezes the four selected parameter vectors. `PROGRESS.md` records the work, and local `STATE.md` carries the handoff because the project-level memory is outside the writable scope.

Primary local sources are `../work_C63/RESULT.md`, `../work_C6_d2/c6.tex` and its corrected ED kernel, `../../contrib/d2_gap_reduction/README.md`, and the stratum, ED, open-chain, FES, and analysis scripts in `../../contrib/d2_grid/`. The grid has no README in the inspected snapshot. Selected grid outputs were read to check the four named candidates, rather than to infer a global phase classification. No files outside C67 were modified. No network commands, cluster jobs, or new tensor-network calculations were run.
