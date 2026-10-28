# C28. A frustration-free charged bus, doped two-mode reference, and a field-rank obstruction

## 1. Outcome and scope

[proved] There is a positive, repeated nearest-neighbor 17-letter initialized core with an exact 2^n-dimensional code, exact total U(1), and C22's nonscalar source resolvent. Replacing the vacuum walk by a Laplacian makes the core frustration-free. This removes an obstacle to applying C8's initialized-core wrapper. The cell is specified completely below, including illegal sectors.

[proved] A uniform density gives an exact common chemical potential, and endpoint densities give two separately adjustable fields at n=2. These operators preserve the entire code, not only its compression. They are retained perturbations of order lambda squared. The positive core continues to annihilate the entire code before those perturbations are added.

[numerical, n=2; heuristic infrared interpretation] Doping the bus repairs the zero-field band gap. The reference below has four simple Fermi points, with two distinct positive Fermi momenta. Its generated inter-flavor pair tunnelling is on the irrelevant side of the one-loop flow. Total U(1) forbids anomalous pairing to every order. A small change of the internal flux removes an additional backscattering resonance. This gives a concrete two-mode chain for central-charge numerics, not an all-n construction or a measured central charge.

[numerical, n=1,2 and two/three cells] Complete legal-sector ED, dimensions 256, 4096, 2304 and 110592, verifies the field and full second-order source-resolvent Hamiltonian. The remainder is O(lambda^4), hence also O(lambda^3), because odd orders vanish. Measured halving exponents are 3.9993 to 4.0012. The calculations include excited vacuum-walker states.

[proved obstruction, restricted scope] Repeated diagonal two-site readouts on this uniform-history code provide at most four independent spatial field profiles, including arbitrary diagonal endpoint terms. They cannot independently specify n fields for n>=5. Scalar Peierls phases on C21's path do not change its return resolvents. Unequal charge-conversion strengths in C22 multiply, rather than repair, its leading hopping matrix. None of these statements excludes a new code, additional noncommuting transport, or a genuine chiral loop construction.

[assumed/open] The requested all-n fixed-alphabet family is not supplied. Nor is there a no-go theorem for every fixed-alphabet implementation. What is established is a positive two-mode candidate and a precise obstruction to its simplest field extension. In particular, the result does not close G1 or kill route A.

## 2. Inputs and conventions

[assumed inputs, attribution] The requested C21, C22, C25, C26 results, C21's cursor table, the route-A background and assembly Section 7 were read first. C8's colored-readout construction and initialized-core wrapper were then read. C22 supplies the 17 letters, arrow grammar, neutral-rail Jordan-Wigner gauge and port maps. C25 supplies the distinction between a positive base and retained payload perturbations, and the ordered fast-wrapper/small-vacancy limit. C26 supplies the requirement to compare against the entire second-order payload and to avoid an extensive isolated-band argument. Claims inherited from those sources are identified below.

[proved, conventions] Write c_i=Z_0...Z_(i-1)|0><1|_i, n_i=c_i^dagger c_i, and P=product_i Z_i. A cell has n data sites followed by one marker site, so k=n+1 physical core sites. Inter-cell fermions use one Jordan-Wigner ordering, first by cell and then by data site. Energies below use t=1 and excited diagonal d=4. The conversion strength g and inter-cell strength lambda are independent small parameters. A coefficient of g^r in the order-lambda^2 bond is denoted C_r, to distinguish it from an order-lambda^r effective Hamiltonian.

## 3. Complete positive core on 17 letters

[proved, definition] The letters are

    L_a, R_a                         a=0,1
    X_(s,a)                         s=v,0,+,-; a=0,1
    M_empty, M_s                    s=v,0,+,-.

They have arrows L=(L,L), R=(R,R), X_s=(L,R), M_empty=(R,L), M_s=(L,L). All M letters are markers. Put a positive penalty K on each arrow mismatch and on every pair whose left letter is a marker. The latter forces the unique marker to the last site once the endpoint conditions are imposed.

[proved, definition] Let Q project onto excited-token letters X_0,X_+,X_-,M_0,M_+,M_-. It vanishes on all other letters. The only conversion entries of C are

    C = g |X_(+,0)><X_(0,1)|
      + g |X_(-,1)><X_(0,0)| + h.c.

Define o=(d-2)Q+C. For |g|<d-2, o is positive, with o >= (d-2-|g|)Q. This follows from its two disjoint 2 by 2 conversion blocks, each with diagonal d-2 and off-diagonal g, and its uncoupled diagonal entries.

[proved, definition] Use the same forward moves for all four species,

    X_(s,a) R_b       -> L_a X_(s,b),
    X_(s,a) M_empty   -> L_a M_s.

For s=v replace each edge by its positive Laplacian (|initial>-|final>)(<initial|-<final|), summed over a,b. Call this sum L_v. For excited species the forward matrix element is -exp(i phi_s), with an additional factor (-1)^a only for s=0. Include its Hermitian conjugate and call the resulting hopping matrix b_exc. The phases used in the final calculation are

    (phi_0,phi_+,phi_-)=(pi/2+0.005,0,0).

[proved, definition] The repeated internal bond and endpoint terms are

    h = (o/2+Q) tensor I + I tensor (o/2+Q)
        + b_exc + L_v + K Pi_bad,
    a = o/2+Q + K Pi_left_bad,
    b = o/2+Q + K Pi_right_bad.                         (1)

Pi_bad is the union of the two forbidden-pair conditions. Pi_left_bad excludes every letter except L and X. Pi_right_bad excludes every nonmarker. Take K=2 in the checks. `ff_bus.literal` specifies every entry of these matrices. Unlisted entries vanish.

[proved, positivity] Each excited hopping edge acts on a disjoint pair of basis states, with one excited token on each side of that edge in configuration space. On that pair,

    Q tensor I + I tensor Q + b_exc
        = [[1,-exp(-i phi)],[-exp(i phi),1]],

with the neutral data sign absorbed in the phase. Its eigenvalues are 0 and 2. On unused basis states the Q term is nonnegative. Thus this operator, o on each site, L_v, and every legality penalty are positive. Equation (1) is a sum of positive operators. The endpoint terms are positive too.

[proved, complete kernel] Zero energy forces every legality penalty to vanish. A legal word starts with L or X, ends with a marker, and contains no earlier marker. Arrow concatenation gives precisely

    L^j X_s R^(n-j-1) M_empty, j=0,...,n-1,
    L^n M_s.

The strictly positive coefficient of Q excludes all excited s. The vacuum Laplacian connects all n+1 remaining positions and acts as the identity on the data. It forces equal amplitudes at all positions for each fixed bit string. Conversely these equal-amplitude vectors annihilate every term. Hence the complete kernel is

    G|psi> = (n+1)^(-1/2) sum_(j=0)^n |j,v> tensor |psi>,
    dim ker H_core = 2^n.                              (2)

This proof includes illegal words and sectors with multiple tokens. No one-token restriction is assumed to prove completeness.

[proved, gap] The complete legal sector has dimension D_n=4(n+1)2^n. Its vacuum block is the path Laplacian, with first positive energy 2[1-cos(pi/(n+1))]. Its excited block has diagonal d, hopping norm at most 2 and conversion norm at most |g|, hence gap at least d-2-|g|. The positive penalties commute with the moves and conversions, since external arrows, marker positions and endpoint legality are preserved. An illegal sector costs at least K. Therefore a complete-core gap bound is

    Delta_n >= min{K, 2[1-cos(pi/(n+1))], d-2-|g|}.     (3)

[proved] Give each letter its data occupation, if any, plus bus charge 0 for v,0 and +/-1 for +,-. The two conversions, all moves and the diagonal terms conserve the sum. On (2) this charge is sum_i n_i. Thus total U(1), and its parity subgroup, are exact microscopic symmetries.

## 4. Fields, including the limitation of diagonal readouts

[proved] The repeated physical density that equals a on every L_a,R_a,X_(s,a), and zero on markers, obeys D G=G sum_i n_i. An endpoint density on the first core site obeys D_L G=G n_0. The repeated two-site diagonal operator that measures the data on its left site only when its right site is a marker obeys D_R G=G n_(n-1). The ordered grammar has only one such bond. These are operator intertwining identities, so there is no leakage out of the code.

[proved] At n=2, D_L and D_R give arbitrary independent coefficients of n_0 and n_1. At n=1 they give the same field. The n=2 hopping below is nonnormal, so these stored-data occupations are not separately conserved band occupations. They provide two independent spatial fields, not an unproved common on-cell rotation into two independent hopping flavors.

[proved] Add lambda^2 alpha_n^4 g^2 sum_i mu_i D_i as a retained perturbation, where alpha_n=(n+1)^(-1/2). The core code remains an invariant subspace of the isolated field-perturbed cell, split only by these fields. A finite-cluster separation follows whenever twice the perturbation norm is below the gap (3). The field terms need not annihilate every code vector. The frustration-free wrapper is applied to (1), before the retained fields and inter-cell exchange are added. Putting a nonscalar field into the common zero-kernel hypothesis would instead contradict that hypothesis.

[proved, field-rank theorem for this code] Allow an arbitrary repeated diagonal two-site operator and arbitrary diagonal one-site operators at the two physical endpoints. Project its field component onto the n linearly independent operators Z_i. In the interior, each neighboring bond sees either LL or RR over a number of cursor positions that is affine in its position, or LX_v or X_vR at one cursor position each. All positions have weight 1/(n+1). A diagonal function of two bits expands uniquely into I,Z_i,Z_(i+1),Z_i Z_(i+1). Therefore the coefficient of Z_i in the sum of the two adjacent bonds is affine in i away from the endpoints. The last bond to M and the left endpoint term can add endpoint corrections. Consequently

    mu_i = A+B i+C delta_(i,0)+D delta_(i,n-1).          (4)

Diagonal nearest-neighbor interactions do not evade this result, since their Z_i Z_(i+1) components are orthogonal to every Z_j under the logical trace. Cancelling those interactions cannot enlarge the projected field span. Thus at most four independent fields are available, and arbitrary n-tuples are impossible for n>=5 in this class. This remains true after a fixed invertible change of coordinates on that n-dimensional field space.

[numerical, n=1,...,8] Direct compression of every diagonal repeated bond entry and endpoint entry gives ranks 1,2,3,4,4,4,4,4. This verifies the analytic upper bound and its saturation at these sizes.

[proved, separate endpoint bound] On a fixed alphabet of size q, a fixed diagonal two-site endpoint readout belongs to a real vector space of dimension q^2. Compression is linear. It therefore cannot yield arbitrarily many independently tunable logical densities from that one endpoint, for any fixed encoding. C8's colored construction uses growing port/color spaces to assign disjoint slots to the signed flavors. Its equation (14) cannot be imported unchanged as an all-n fixed-alphabet readout. This bound concerns independent continuous controls for a fixed code. It does not forbid selecting one useful incommensurate list or changing the encoding with n.

[proved, general linear-control bound] More generally, a repeated Hermitian two-site perturbation and two Hermitian endpoint perturbations have at most q^4+2q^2 real coefficients. For a fixed code their compression is a linear map from that space. Exact arbitrary independent n-field control by these additions is impossible once n exceeds this bound. This is not an obstruction to unbounded c, which needs one suitable member for each n, rather than an open n-dimensional set of exact field controls on the same code.

## 5. The source contraction and two failed phase repairs

[proved, inherited gauge with changed vacuum amplitude] Gauge the neutral rail at position j by Z_0...Z_(j-1). The complete excited block is exactly C22's three-rail H1, with diagonal d, scalar phase-dressed rail hopping and conversions H_(+,0)(j)=g c_j, H_(-,0)(j)=g c_j^dagger. The left source maps are S_L,s=|0,s> tensor I. On the right they are S_R,0=|n,0> tensor P and S_R,+/-=|n,+/-> tensor I. The vacuum source amplitude is now alpha_n=1/sqrt(n+1), rather than C22's sine endpoint amplitude.

[proved, definition] Use C22's local charge-conserving ports

    F0=|0><v|+|v><0|,
    F+=|+><v|+|v><-|, F-=F+^dagger,
    V=F0_R F0_L+i F+_R F-_L-i F-_R F+_L.              (5)

The norm of its active token matrix is sqrt(3). It changes a vacuum token to an excited species, so G^dagger F_s G=0. With

    B=S_R,0 tensor S_L,0+i S_R,+ tensor S_L,--i S_R,- tensor S_L,+,

the exact second-order coefficient is

    K_bond(g)=-alpha_n^4 B^dagger
                (H1 tensor I+I tensor H1)^(-1) B.       (6)

Distinct bonds have zero second-order cross term. An outer cell touched once remains in an excited species until the final projection, which kills it. This is valid on arbitrary open clusters and locally on rings.

[proved, expansion and coefficients] Write the pair excitation Hamiltonian as A+g C, R=A^(-1), and remove the factor alpha_n^4. For |g| ||R C||<1 its Neumann expansion gives

    C_r=-(-1)^r B^dagger R(C R)^r B.
    C_0=-B^dagger R B,
    C_2=-B^dagger R C R C R B,
    C_4=-B^dagger R C R C R C R C R B.                  (7)

Every conversion changes neutral to charged or vice versa. B creates either two neutrals or two charged species, so a nonzero return has an even total number of conversions. Odd coefficients vanish. At order g^2, C22's equal-rail anticommutator cancellation makes the on-cell contribution scalar, and leaves an exactly quadratic hopping. At order g^4 there are at most four logical ladder factors, so normal ordering gives only scalar, quadratic and quartic operators. Equation (7), rather than a finite difference of nearly scalar matrices, computes all the coefficients used here.

[proved] Uniform rail flux still does not give independently adjustable sea centers. C22's formula remains

    T_ij=2i exp(i phi_link) F_ij exp[i chi(n-i+j)]
                       sin[kappa(n-i-j)],              (8)
    kappa=phi_0-(phi_++phi_-)/2, chi=(phi_--phi_+)/2,

where F_ij is C22's positive scalar heat-kernel integral. Our choice phi_link=pi/2 and chi=0 makes T real. The U(1) route, not independent complex phases, is used below.

[proved, failed repair] Replace the two conversions by g_+ and g_-. Each hopping process converts one positive and one negative source across the two cells. Its prefactor is g_+ g_-, for both interfering paths. Same-cell returns still pair into anticommutators separately for each charged rail. Thus the entire order-two hopping (8) is multiplied by g_+ g_-/g^2, with no new field or hopping phase. An imbalance is not a repair of the reversal geometry.

[proved, failed C21 repair] On C21's cursor path, replace scalar edge hopping by -t exp(i phi_j) times its existing unitary decoration. Choose scalar vertex phases recursively to remove every phi_j. These phases commute with every operator-valued on-site potential. An endpoint source only acquires a scalar phase, which cancels in S^dagger f(H) S. Both return measures, and hence C21's J=-2 R G_f, are unchanged. This proof covers arbitrary scalar phases on this path, including a choice called an oriented clock. A directed non-Hermitian shift is outside the problem's Hermitian Hamiltonians. A loop with gauge-invariant flux is outside the theorem.

## 6. The doped n=2 chain

[proved, analytic reference] At kappa=pi/2, n=2, equation (8) has T=[[0,b],[c,0]], with b,c nonzero real numbers and |c|>|b|. Add a common positive density mu. Its single-particle bands are

    epsilon_+/-(p)=mu +/- sqrt(b^2+c^2+2bc cos(2p)).     (9)

Choose 0<k0<pi/2 and mu=sqrt(b^2+c^2+2bc cos(2k0)). The lower band has four simple roots +/-k0 and +/-(pi-k0), since bc sin(2k0) is nonzero. The upper band is positive. The four roots represent two conducting right-left pairs in one partially occupied band, not two separately diagonalized on-site flavors. Simplicity ensures persistence under sufficiently small changes of flux and field, by the implicit-function theorem. Thus the zero-field gap in C22 is not an obstruction at arbitrary chemical potential.

[proved, resonance warning] At that symmetric reference the two positive roots add to pi. This permits an additional inter-pair backscattering process. The reference is not used unchanged for the stability test.

[proved, certificate at the chosen nonzero flux shift] The four crossings also have a finite algebraic certificate at epsilon=0.005. Set A=8I-adjacency(P_3 square P_3), R=A^(-1), s=(2,0), and v=R|s>. The scalar heat-kernel integral in (8) equals

    F_ij=2 v^T(P_i tensor I) R(I tensor P_j)v.

This is the sum of the two conversion orders in the pair resolvent; the two terms coincide by transposition. Exact rational inversion gives

    F_00=1431489/334833680384,
    F_01=108377/167416840192,
    F_10=26801113/167416840192.

Write a=2F_00 sin(2epsilon), b=-2F_01 cos(epsilon), c=-2F_10 cos(epsilon), and use the stated mu^2=b^2+c^2+2bc cos(1.4). The determinant of the Bloch matrix is, with x=cos(p),

    D(x)=-4bc x^2+2mu a x+mu^2-(b-c)^2.

Its value at zero is positive. Its values at both endpoints +/-1 are negative if mu|a|<bc[1-cos(1.4)]. The elementary inequalities cos(1.4)<0.18, |sin(2epsilon)|<=2epsilon and cos^2(epsilon)>=1-epsilon^2 reduce this to the rational inequality

    8F_00(F_01+F_10)/200
        < (82/100)4F_01 F_10(1-1/40000).

The ratio of the right to the left side is greater than 12.36. Alternating Taylor bounds certify 0<cos(1.4)<0.18. Therefore the concave quadratic D has one simple root in (-1,0) and one in (0,1). Each gives two simple momenta. The other band stays positive because the trace 2mu+2a cos(p) is positive, as follows from the same displayed bounds. `certificate.py` implements these exact rational steps. Root coordinates and velocities in the table are numerical evaluations of this proved four-root situation.

[numerical, n=2] Set k0=0.7 and kappa=pi/2+0.005. Evaluate b,c at this flux and keep the same formula for mu. In the leading g^2 units of (7), the resulting hopping and field are

    T = [[ 8.55030557456e-8, -1.29468033402e-6],
         [-3.20168245393e-4, 0]],
    mu = 3.2039083881467e-4.                             (10)

Imaginary entries and same-cell quadratic terms at this order are zero analytically, with numerical remnants below 3e-18. The roots and velocities are

| pair | right momentum | left momentum | positive speed |
| :-- | --: | --: | --: |
| 1 | -2.416591538066 | 2.416591538066 | 2.62593151981e-6 |
| 2 | 0.673661697190 | -0.673661697190 | 2.46936224408e-6 |

[proved, chain convention] Let W=lambda^2 g^2/9, since alpha_2^4=1/9. The leading logical chain is

    H_lead/W = sum_j [c_j^dagger T c_(j+1)+h.c.
                     +mu(n_(j,0)+n_(j,1))].            (11)

Here c_j is the two-component column of globally graded fermions. Add W g^2 sum_j C_4(j,j+1), with the scalar omitted, for the first generated interaction. The quadratic part of C_4 must also be included. The normal-order convention is given in Section 7. These formulas specify the n=2 chain without inferring a central charge from hopping rank.

[numerical, handoff] `logical_chain_n2.npz` contains the actual 16 by 16 nearest-neighbor spin bonds in the cell basis |00>,|01>,|10>,|11>. `leading_bond` is (11) with half the density on each adjacent bond. `g4_bond` is C_4 with its scalar subtracted. `full_g002_bond` is the full-in-g order-lambda^2 contraction at g=0.02, divided by g^2, with the same fields and scalar subtraction. The common omitted energy scale is W. Thus an ordinary dimension-four chain can use these matrices directly. For open boundaries, add the missing half field at each end if the intended field is uniform. No central-charge estimator has been run here.

## 7. Generated pairing, pair transfer and its one-loop sign

[proved] Exact total U(1) forces every anomalous quadratic term c c or c^dagger c^dagger to vanish in a symmetry-covariant canonical effective Hamiltonian. The same holds for every number-changing higher vertex. This statement is to all orders, independent of perturbative power counting. Pair transfer conserves total number and is not forbidden.

[proved, extraction convention] For the four modes on a bond, subtract the vacuum scalar and the one-particle matrix from C_4. In the ordered two-particle basis c_i^dagger c_j^dagger|0>, i<j, the remaining matrix U fixes the quartic operator uniquely as

    sum_(i<j,k<l) U_(ij,kl) c_i^dagger c_j^dagger c_l c_k. (12)

For a normalized Bloch spinor u_r at p_r, let w_r=(u_r,exp(i p_r)u_r) on the two-cell bond, and (w_r wedge w_s)_(ij)=w_(r,i)w_(s,j)-w_(r,j)w_(s,i). The antisymmetrized continuum scattering coefficient is

    A_(rs;tu)=(w_r wedge w_s)^dagger U(w_t wedge w_u).    (13)

This follows by substituting c_(j,a)=sum_r exp(i p_r j)u_(r,a) psi_r into (12). The sum over bond origins enforces the net momentum condition. It includes all density and exchange contractions, including terms supported within just one of the two cells.

[numerical, n=1] C_4 has the density interaction 1.96112889419e-5 n_j n_(j+1), in addition to its scalar and quadratic terms. There is only one conducting pair, so no inter-flavor pair transfer. The n=1 field is 1.11234141966e-4 and the leading roots are +/-1.384724464464.

[numerical, n=2] Order the pairs as in the table above. In the spinor gauge returned by the recorded eigendecomposition, equation (13) gives, before multiplying by lambda^2 alpha_2^4 g^4,

    G_ab = A_(R_a L_b;R_a L_b)
         = [[ 8.67874850129e-7, -1.64536810754e-7],
            [ 2.94385647532e-7,  7.73651048524e-7]],
    p = A_(R_2 L_2;R_1 L_1)
      = -6.45074367714e-7 + 8.04111022848e-8 i.           (14)

All quartic coefficients, and the quadratic C_4 matrix, are saved in `final_algebra.json` and `coefficients_n2.npz`. The finite-dimensional normal-order residual above degree four is 8.4e-20. The order-g^2 quartic residual is 6.0e-18.

[assumed, continuum input] Linearize at the four simple crossings and use the usual weak-coupling fermion-shell RG, with densities rho_R,L and pair operator psi_R2^dagger psi_L2^dagger psi_L1 psi_R1. Its coefficient p can be made real by a constant relative rephasing of the two flavors. Its sign alone therefore cannot decide relevance. At equal speeds v, bosonization gives

    H_G=(1/(2pi))[(partial phi)^T(v I+G_s/(2pi))(partial phi)
                    +(partial theta)^T(v I-G_s/(2pi))(partial theta)]

up to the mixed derivative term from the antisymmetric part of G. To first order K=I-G_s/(2pi v). The pair vertex is exp[2i(theta_1-theta_2)], with dimension

    Delta_pair=2+(G_11+G_22-G_12-G_21)/(2pi v)+O(G^2).  (15)

The antisymmetric forward part contributes to this dimension only beyond first order. The physical forward couplings in this equation are lambda^2 alpha_2^4 g^4 times (14), and the speeds are W times the table entries.

[assumed, one-loop prescription with unequal speeds] The opposite-chirality shell integral with speeds v_a,v_b is d ell/[pi(v_a+v_b)]. Keeping the velocities fixed at this order gives

    dot p = x p,
    x = -G_11/(2pi v_1)-G_22/(2pi v_2)
                         +(G_12+G_21)/(pi(v_1+v_2)),
    dot G_11 = -|p|^2/(2pi v_2),
    dot G_22 = -|p|^2/(2pi v_1),
    dot G_12 = dot G_21 = |p|^2/[pi(v_1+v_2)].          (16)

The negative terms are the pair-channel contractions into the same flavor; the positive terms are the crossed contractions between the two flavors. Differentiating the displayed x, rather than assuming a sign for it, gives

    dot x = |p|^2/pi^2 [1/(2v_1 v_2)+2/(v_1+v_2)^2].

Define y=|p|/pi sqrt[1/(2v_1 v_2)+2/(v_1+v_2)^2]. Then dot x=y^2 and dot y=x y to one loop. The invariant is x^2-y^2. For x<-|y|, y decreases to zero and x tends to a negative value. This includes the equal-velocity normalization y=|p|/(pi v), consistent with (15). Equations (15)-(16) are continuum weak-coupling input and its algebraic consequences, not a finite-dimensional lattice theorem.

[numerical matching; heuristic flow conclusion] Substituting the actual microscopic coefficients gives

    x/g^2 = -0.09435238259,
    y/g^2 =  0.08124029522,
    (-x-y)/g^2 = 0.01311208737 > 0.                     (17)

Thus the generated pair transfer is irrelevant in this one-loop theory. The sign is fixed by the computed relative forward stiffness and its magnitude relative to pair transfer, not by negative semidefiniteness of the Schrieffer-Wolff correction. Higher orders in g shift (17) at relative order g^2. At fixed n=2 a strict leading inequality leaves a sufficiently small nonzero g interval in this perturbative description.

[numerical, momentum enumeration] Enumerating the six distinct two-branch pair sums leaves precisely one equality between different pairs, the Cooper transfer between R_1 L_1 and R_2 L_2. The smallest other mismatch modulo 2pi is 0.1026788367. In particular the extra backscattering resonance at kappa=pi/2 is removed. Coincident-branch local pairs vanish by Fermi statistics; their derivative descendants have larger free dimension. This is the U(1) version of the detuning check, not the assertion that all four points satisfy the parity-only Sidon criterion.

[numerical, finite-g coefficient check] Extracting the quartic part of the full second-order resolvent gives p/g^4=(-6.57272543e-7+8.16097321e-8 i) at g=.2 and (-6.48098588e-7+8.07088290e-8 i) at g=.1, before the common lambda^2 alpha_2^4 factor. They approach (14). The full g^4 quadratic term is complex and can shift the pair centers at the next order. The equality of Cooper-pair momenta used in (16) is a leading-order reference property, not an asserted exact symmetry of the finite-g dispersion. Additional detuning would cut off that pair process.

[heuristic, limitations of the flow] To use this reference in the microscopic chain, first take g sufficiently small that its root shifts, higher interactions and corrections to (17) fit inside these strict margins. Then take lambda small relative to the full cell gap and the very small conducting velocities. The lambda^4 correction cannot be bounded relative to a weak channel merely by lambda/Delta_n. The wrapper introduces further ratios and momentum margins. The ED couplings below test perturbative matrix identities; they are not certified points in a c=2 phase. No all-n one-loop sign theorem is claimed.

## 8. Exact diagonalization and error scaling

[proved, finite-cluster symmetry] On alternating cells multiply every excited-token state by -1, and leave all vacuum-token states unchanged. H_core and the density fields commute with this unitary. Every bond in an open chain changes sign, so H(lambda) is unitarily equivalent to H(-lambda). The unitary is the identity on the code. Canonical polar identification of the low subspace is symmetry covariant. Its effective Hamiltonian is therefore even in lambda, also when the fields are proportional to lambda^2. Thus there is no cubic remainder to measure. The next possible term is lambda^4.

[numerical, method] `verify.py` constructs the full product of legal cells, including all vacuum and excited rail states, with the physical density fields. It diagonalizes every total-charge sector intersecting the code and retains all 2^(n M) code levels, plus a next level in each sector. QR/Rayleigh-Ritz resolves near degeneracies. The polar factor of G^dagger U_low identifies the canonical effective matrix. The comparison is with (6) at finite g=0.2 plus the exact compressed fields, not with its g^2 truncation.

[proved, omitted-charge exclusion] A sector with total charge outside 0,...,n M requires at least one excited token. It has unperturbed energy at least d-2-|g|. The exchange lowers this by at most |lambda|sqrt(3)(M-1). The positive fields used here do not lower it. At the largest tested lambda=0.08 this bound exceeds 1.52. The next levels in the relevant sectors are about 1 at n=2. At n=1 the exact excited-rail lower bound can instead use its one-edge norm, d-1-|g|=2.8, so the omitted sectors lie above 2.52 and do not enter the reported low/next band. The positive illegal-sector penalties provide a separate exact unperturbed exclusion; full-ring marker selection at these couplings is not being simulated.

[numerical] The following are operator-norm errors including the scalar. The last column is log_2(error(.08)/error(.04)).

| n | cells | full legal dimension | error at .08 | error at .04 | exponent |
| --: | --: | --: | --: | --: | --: |
| 1 | 2 | 256 | 4.29647497e-7 | 2.68469910e-8 | 4.00032 |
| 1 | 3 | 4096 | 3.72251130e-7 | 2.32761921e-8 | 3.99935 |
| 2 | 2 | 2304 | 3.61189463e-7 | 2.25565447e-8 | 4.00114 |
| 2 | 3 | 110592 | 4.98486962e-7 | 3.11448772e-8 | 4.00049 |

[numerical] Traceless errors at lambda=.08 are 2.0231e-11, 2.0642e-12, 2.2645e-11 and 2.7399e-11 in the same order. Their halving exponents are 3.9995, 3.9952, 4.0016 and 4.0003. Eigenpair residuals are at most 7.14e-13. The smallest code-overlap singular value is 0.9999168. The smallest next gap is 0.9990118. Charge commutators and field-intertwining residuals are zero in stored arithmetic. The smallest traceless errors, particularly at lambda=.04, are close to numerical residuals and should not be assigned more precision than those residuals support.

[numerical, independent local checks] The literal core of (1), dimensions 289 and 4913 at n=1,2, agrees exactly in stored arithmetic with its legal restriction and has zero coupling out of that subspace. Its complete low spectra give kernel dimensions 2 and 4, followed by gaps 2 and 1. The local bond's smallest numerical eigenvalue is -3.5e-15, consistent with its analytic positivity. The Neumann truncation (7) through g^4 has full bond errors proportional to g^6 at g=.2,.1,.05, with halving exponents 5.996 to 6.004.

[numerical, failed solver and correction] An initial unconstrained ARPACK call on the 289-dimensional literal matrix missed its exact zero eigenspace and returned levels starting at 2. That result failed the known-code check. The n=1 check was replaced by dense Hermitian diagonalization. The larger check uses an explicit seeded start, a larger Krylov space and assertions for every expected zero eigenvalue. The final saved spectrum passes. No conclusion uses the failed output.

## 9. Alphabet accounting and what remains open

[proved] The initialized positive core alphabet is q_core=17 for every n. This is a concrete module count, not an alphabet count for a completed unbounded-c family.

[assumed composition; proved arithmetic] If this core is placed inside the supplied C7 frozen-selector architecture, the direct tensor count is 4 q_core=68. C8's full mobile selector and two-helix construction has count 96 q_core+48=1680. These are architecture counts using the inherited selector/wrapper theorems, not newly verified complete finite-coupling models. The frozen length range and any padding for the smallest test cores must follow the selector hypotheses. The isolated n=1,2 ED does not use a selector or padding. C8's selector-wrapper statement uses k>=3, so n=2 already meets its core-length condition.

[proved, interface qualification] Equation (1) meets the positive initialized-core hypothesis required for tagged vacancies. Thus the unperturbed wrapped code and identity payload transport can inherit C8's theorem, unlike C22's original adjacency-walker ground state. Field readouts remain retained perturbations. Two-site endpoint densities require the corresponding two-occupied-site visibility compensator. Exact unwrapped fields do not imply exact finite-vacancy equalization.

[assumed/open, second-order composition] The new bus mixes syndromes. C25's exact vanishing of its particular u lambda term is not inherited. Its scalar length-dependent self-energy also cannot be discarded as one global constant. The fast-wrapper and small-vacancy estimates must be applied to this source Hamiltonian, with the smallest wrapper gap and the compensator norms included. No numerical wrapped-bus calculation was performed here.

[proved, obstruction scope] Formula (4) excludes all-n independent fields produced by repeated diagonal readouts on this code. The path-phase theorem excludes scalar cursor phases as a C21 repair. The conversion-product theorem excludes imbalance of C22's two conversion strengths as a leading phase repair. None excludes a prescribed useful field pattern with fewer controls, non-diagonal readouts, a different entangled code, an additional internal loop, or a new controller. In particular the n=2 result shows that the original zero-field gap alone cannot be used as a no-go against the charged bus.

[assumed/open, requested all-n ingredient] There is no proof that this family has n conducting pairs for arbitrary n, no independent n-field implementation beyond the restricted rank bound, and no all-n sign control of the resonant pair-transfer system. Therefore neither 68 nor 1680 is a demonstrated d0 for unbounded c. The task is only partially resolved. The useful next steps are (i) test whether a prescribed low-parameter density profile can produce n crossings despite (4), (ii) derive the resulting many-pair RG matrix, or (iii) change the code to obtain more spatial readout profiles. A full target should precede another alphabet optimization.

## 10. Reproduction and files

[proved, file accounting] All new files and computation outputs are inside work_C28. `bus_base.py` is a read-only-source copy of C22's `spatial_bus.py`; the new positive local rule table is in `ff_bus.py`. `explore.py` contains the resolvent coefficient and normal-order/scattering routines, including the unsuccessful conversion-imbalance scan. `verify.py`, `checks.py`, `certificate.py` and `handoff.py` produce the final records. They use Python's standard library, NumPy and SciPy; the rational certificate uses the standard-library Fraction type. The project STATE is outside the permitted folder and was not modified. Local `STATE.md` and `PROGRESS.md` carry the handoff. No network, subagents, cluster job, external messages, commit or push were used.

[numerical, reproduction] Run with NumPy/SciPy and one BLAS thread. Every recorded individual job took at most about two minutes.

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B checks.py
    python3 -B certificate.py
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B verify.py --algebra
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B handoff.py

[numerical, ED reproduction] For n=1,2, cells=2,3 and lam=.08,.04, run

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B verify.py --n 2 --cells 3 --lam .08

[proved, deliverable map] `final_algebra.json` stores the fields, hopping, quadratic/quartic g^4 coefficients, roots and one-loop matching. `checks_results.json` stores literal-core checks, field ranks and g^6 remainders. `ed_n*_M*_lam*.json` and matching effective matrices store all eight ED comparisons. `logical_chain_n2.npz` is the direct numerics handoff. `algebra_results.json` is an earlier failed imbalance scan and is not the final reference. No thermodynamic or central-charge result is contained in those files.
