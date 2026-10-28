# C40 derivations

[proved, scope] All new finite claims below have analytic proofs or exact rational/integer checks. Statements about continuum matching and the infinite infrared flow are tagged assumed or heuristic. A deduction conditional on a named hypothesis is tagged proved conditional. Source results retain their source attribution and have not been independently certified here. No network, subagents, external writes, or large computations were used.

## 1. The fixed-n problem and its parameters

[assumed, finite input F] Fix n, all cell parameters, and a legal selected sector. The complete code has n logical fermions, the unperturbed cell Hamiltonian has this complete kernel and a positive excluded gap Delta_n, and boundary perturbations have finite norms. For the mobile construction, both selected lengths have the same code. The full-chain application also requires exclusion of competing illegal, wrong-length, and marker-free ground sectors. These are separate finite/operator inequalities, not a consequence of a Bloch-band calculation.

[proved, inherited] C32 PROOFS sections 1-3 supply the positive 28-letter core, d0=2736 mobile wrapper, and arbitrarily accurate second-order equalization, conditional on the explicitly inherited C8 wrapper theorem. The smallest wrapper gap includes the charge tunneling splitting proportional to U(w/U)^(2n+4). Replacing it by U is invalid. The frozen 23-letter README explicitly leaves some unrestricted-sector exclusion bounds unfinished.

[proved, a legal-cell gap bound for the original 23-letter bond] At t=t_x=t_sw=1, 0<=t_d<=1, nonnegative D_x,D_M, equal conversions g_+=g_-=g and zero carrier splitting, the excited token graph has absolute row sum at most two. It is the folded path with a link joining its endpoints. The departure Z factors and phases do not increase this bound. Thus its kinetic operator is bounded below by -2I. The conversion is a direct sum of two-state blocks of norm |g|. In every legal cell, of any length,

    H_exc >= (Delta_mass-2-|g|) I,

apart from the common marker chemical potential. The vacuum path ground energy is -2 cos(pi/(l+1))<=0 for l token positions. Hence Delta_mass>=2+|g|+gamma, gamma>0, ensures a legal excited-token gap at least gamma for every l. With nonzero splitting delta, replace |g| by |g|+|delta|. This conservative bound does not prove the illegal-arrow exclusion. A data field can be separated into the conserved total charge and a bounded carrier-charge correction before estimating excitation denominators.

[assumed, finite input B] The selected leading quadratic reference has exactly 2n distinct, simple Fermi roots p_i, n of each chirality, positive speed minimum bar_v_*, and a positive margin for every nontrivial number-conserving quadratic/quartic momentum. Also assume a positive curvature/patch isolation scale. All these statements concern finite matrix functions of momentum and require proofs or certificates for a claimed n. Crossings alone are insufficient.

[proved, definition of the margins] With circular distance d on the reciprocal circle, a sufficient payload test is

    delta_2 = min_(i!=j) d(p_i-p_j,2pi Z)>0,
    delta_4 = min_({i,j}!={k,l}) d(p_i+p_j-p_k-p_l,2pi Z)>0,

where i<j and k<l. Allowing repeated points is a stronger sufficient test, used by C33. Equal unordered pairs are excluded from delta_4 because they give densities. Pairs sharing one index reproduce a quadratic mismatch. Exact total U(1) removes sums appropriate to number-changing operators from the required list.

[proved, normalization] Let g denote conversion strength, called tau in C31. A convenient leading normalization is

    W = lambda^2 g^2                     frozen, in fixed core units,
    W = (lambda^2/s) g^2                 C32 mobile core scale s,
    H_ref = W Hbar_pay + W j_l Hbar_length.

Express attenuation, endpoint factors, Gram eigenvalues and phases remain in Hbar_pay. They are not discarded by this definition of W. Thus physical speeds are W times their barred values. Choose u beta_n=W j_l in the mobile case. The j_l used here differs by a factor g^2 from C32's coefficient before conversion is factored out. Either convention is valid if it is used consistently.

[proved] The bandwidth ratio and all barred speeds are independent of lambda after the core and g are fixed. Decreasing lambda lowers every retained scale together. It does not make a lambda^2 g^4 correction small relative to lambda^2 g^2 hopping. That relative error is controlled by g and the very small propagation factors in Hbar_pay.

## 2. Connected locality and size at every fixed order

[proved, inherited representation] Use C8 section 46 and C36 section 2. In rooted cell coordinates the mobile legal sector has fixed tensor factors

    K = direct_sum_l K_l,   A=sum_j A_j,   A_j>=Delta_n(1-P_j).

P_j includes both selected lengths and their entire logical codes. Perturbations act on neighboring cell factors and are bounded. They preserve total physical length. Quotienting the rooted representation gives the physical cyclic gluing. Frozen cells are the simpler special case with one selected length. The following estimates are independent of total ring size and of the upper excited energy.

[proved, finite-cluster analytic estimate] First write H_I(z)=A_I+z V_I on an interval of l cells, with ||V_I||<=l J, A_I>=Delta Q_I and zero on P_I. On the circle |z|=Delta/(16lJ), use the spectral contour |xi|=Delta/2. Then

    ||(xi-A_I)^(-1)||<=2/Delta,
    ||Pi(z)-P_I|| <= [2|z|lJ/Delta]/[1-2|z|lJ/Delta] <=1/7,
    ||H_I(z) Pi(z)|| <= (Delta/2)/[1-2|z|lJ/Delta] <=4Delta/7.

Here Pi is the Riesz projection. These estimates follow by a norm Neumann series and the contour length 2pi(Delta/2). They remain valid for the direct-sum cell space because V_I is bounded and the contour encloses a finite-rank kernel.

[proved] Set D=Pi-P_I and use the analytic direct rotation

    U=(Pi P_I+(1-Pi)(1-P_I)) (1-D^2)^(-1/2).

For real z it is unitary and maps P_I to Pi. The binomial branch at z=0 defines it for complex z. Both U and U^(-1) have norm at most

    (1+1/7)/sqrt(1-1/49)=2/sqrt(3).

Consequently the canonically identified finite-cluster Hamiltonian

    K_I(z)=P_I U^(-1) H_I(z) Pi(z) U P_I

is analytic and has norm less than Delta on this disk. Cauchy's coefficient estimate gives

    ||[z^r] K_I(z)|| <= Delta [16 l J/Delta]^r.       (1)

This proves a finite-dimensional coefficient bound without assuming a gap above an extensive interacting code band.

[proved, connected subtraction] Give different bonds independent variables. If a set of bonds is disconnected, the low subspace and its polar identification on the code factorize over the disconnected intervals. Thus K is a sum of their identified Hamiltonians. Every Taylor coefficient involving disconnected components vanishes. Equivalently, for a connected interval I sum K over all subsets of its l-1 internal bonds with signs (-1)^(number omitted). This picks out terms using every bond of I. It has at most 2^(l-1) summands. A degree-r connected coefficient therefore has l<=r+1.

[proved] Assign split on-site fields to their incident bonds; adding explicit one-site terms only changes the harmless counting constant. On the infinite line a specified cell lies in l intervals of length l. From (1), a deliberately loose local-norm bound for the order-r coefficient is

    ||Phi_r||_mu <= C_r(mu) J^r/Delta^(r-1),
    C_r(mu)=(r+1)^2 2^(r+1) exp(mu(r+1)) [16(r+1)]^r.  (2)

Here ||Phi||_mu=sup_j sum_(I containing j) exp(mu |I|)||Phi_I||. One can increase C_r to include separately represented on-site terms. Formula (2) is a majorant, not a lower bound on actual coefficients. Its growth does not justify summing r at fixed coupling.

[proved, several small couplings] Apply the same construction with independent variables for port, marker and field terms. A coefficient with r port and t marker insertions has a finite bound C_(r,t) B_V^r B_T^t/Delta^(r+t-1), multiplied by lambda^r u^t. Bounds depend on the full compensator norms, not on projected amplitudes. Conserved charges commute with every resolvent, projection and direct rotation, so the coefficients preserve the microscopic U(1), marker number and translation covariance.

[proved, physical spatial support] A connected r-insertion process starts in at most r+1 selected cells. Marker transfers move one boundary site per insertion. Its support in physical coordinates is at most (r+1)(b+1)+2r sites, a conservative bound including moving ends. Thus locality at fixed order is proved also in physical units. The bound is not uniform in n or b.

[proved, fixed-order temporal bounds] On an irreducible unperturbed excursion an intermediate state has at least one cell defect, so every inverted energy is at least Delta. A reduced resolvent power has

    E^(-k) = integral_0^infinity t^(k-1) exp(-E t) dt/(k-1)!.

A time moment of order m gives m!/E^(m+1) for k=1 and the corresponding factorial ratio for general k. In a word with r-1 fast intervals, r>=2, integrating their total-duration m-th power costs at most

    (r-1)(r)...(r+m-2) Delta^(-m)

relative to its product norm majorant. Frequency derivatives have the same finite bounds. Folded normalization terms are included as products of such resolvents and instantaneous retained operators. An interval wholly in P is retained propagation and is never assigned 1/Delta. These are bounds on finite-order irreducible words, not on the exact bare-source correlator.

[proved, explicit fourth-order check of the construction] For PVP=0, R=QA^(-1)Q and B=QVP, let A2=B^dag R B and T2=B^dag R^2 B. The canonical coefficients are

    K2=-A2,
    K3=B^dag R V R B,
    K4=-B^dag R V R V R B + (T2 A2+A2 T2)/2.          (3)

C26 derives (3) from the invariant graph and its metric. The last term is required. It cancels disconnected pieces and incorporates the energy dependence of lower-order self-energies. Dropping scalar source terms before this normalization can give a wrong higher-order logical interaction.

[proved, limitation] None of (1)-(3) proves an exact local effective action at fixed nonzero coupling. The required perturbative radius on l cells decreases with l. In particular there is no permissible inference from (2) to a geometric all-order bound with a constant independent of r.

## 3. The additional connected-action input, stated precisely

[assumed, M2] After locally accounting for vacuum dressing and subtracting propagation reducible through retained fields, there is an exact low-external-frequency Wilsonian action for the retained cell variables. Its vertices preserve the local fermionic operator algebra, the exact charges and the mobile gluing. Their connected spatial kernels are summable with a fixed exponential weight, or sufficiently many spatial moments for the same matching. Time/frequency moments needed for a derivative expansion are summable as well. The action has a uniform weak-coupling asymptotic expansion through the known second-order coefficient, with constants independent of the number of cells. Its norm must control coefficients and field/filling derivatives, not merely total energy density.

[assumed, a concrete norm for M2] One sufficient convention expands the residual action before its time-derivative expansion in antisymmetric fermionic kernels K_(X,p)(t_1,...,t_(p-1);0), with lattice fields dimensionless, p external legs, and one time fixed by translation invariance. Use the coefficient l1 norm in a fixed local fermion basis, and require, for some a>1 and mu>0,

    sup_j sum_(X containing j) exp(mu |X|)
      sum_(p>=2) a^p/p! integral [product dt_i]
      (1+Delta_n sum_i |t_i|)^2 |K_(X,p)(t;0)|_1 < infinity.

Include instantaneous vertices as delta measures, scalar length interactions, and the first two derivatives with respect to background fillings. The reference kinetic term and known second-order matching terms are subtracted. The length hard-core variable can be represented by its usual fermion field in cell coordinates. A smooth low-frequency restriction is part of this Wilsonian definition; this is not an operator norm on all retained many-body energies. Other equivalent matching norms are possible. Establishing this bound for exact irreducible kernels, rather than merely specifying it, is the missing part of M2.

[proved, finite-order compatibility] A local operator on l retained cells of dimension D per cell has an expansion in D^(2l) normalized unitary matrix-basis elements, each coefficient bounded by its operator norm. Its coefficient l1 norm is therefore at most D^(2l) times its operator norm. A fixed change to a fermion-monomial basis adds another finite factor exponential in l at fixed n. Thus the finite-order bounds in section 2 do control coefficient norms after such factors are included. They do not control the infinite sum over support, order and external legs in M2.

[assumed, usable quantitative version] Fix the core, g, wrapper, and n. Write q=|lambda|B_V+|u|B_T+|h_field|B_F, including all physical norms. A sufficient low-frequency bound is

    ||R_>=3||_(LE,mu) <= A_n q^3/Delta_n^2,
    ||K2(omega)-K2(0)||_(LE,mu)
                       <= A_n q^2 |omega|/Delta_n^2,       (4)

for the irreducible remainder and the leading fast kernel, uniformly for external frequencies below a fixed fraction of Delta_n. The same norm controls a sufficient finite number of spatial moments and derivatives in the reference fillings. Mixed lower-order terms already included in the reference are subtracted first. A corresponding little-o(lambda^2) bound would also suffice at fixed g; analyticity or convergence to all orders is not required. Scalar length interactions are part of the norm.

[proved conditional on M2] For u=O(lambda^2) and a field of retained order, q=O(lambda). Since W=O(lambda^2) at fixed g, the first bound in (4) divided by W tends to zero. At retained frequencies O(W) the second bound divided by W is O(lambda^2/Delta_n^2). Thus the full frequency dependence is a small local matching correction. Wavefunction normalization is included. This is the second-order version of C36's first-order estimate; retaining the entire lambda^2 bond before defining the error is essential.

[heuristic, physical content] The mechanism for M2 is the energy cost of a local cell excitation, together with a small intercell coupling. It should produce local irreducible vertices once long retained propagation is excluded. At each fixed n there is a large separation Delta_n/W. The nontrivial part is performing that subtraction consistently in an interacting, highly degenerate retained background. Calling this the usual integration of massive fields states a physical input. It does not prove (4).

[proved, located obstruction] C36 section 3 proves a geometric bound only for histories in which each insertion touches a currently excited cell. A general connected history can create a new defect pair in an unexcited region and join the older branch later. It can also carry a modified logical state after its last local defect has disappeared. The number of available insertions then is not bounded by the number of current defects. The cancellation of these histories against reducible retained propagation is the missing estimate. Reusing the active-defect count for them is incorrect.

[proved, inherited controls] C8 sections 44 and 51 exhibit respectively a closing global bare-Q threshold despite an exactly removable massive sector, and an algebraic tail of a bare off-code source despite a local massive gap. Subtracting a one-point function does not remove the latter tail; an operator-level retained component must be subtracted. C36 also disproves an unrestricted all-retained-energy interaction-picture norm bound. These examples invalidate shortcuts to M2, not M2 itself or route A. C8's smooth spectral filtering gives source bounds, but does not construct the connected action.

[assumed, sector stability component of M2/F] The local reduction must also retain a positive cost of excluded machinery and length sectors and compare invariant illegal sectors with the dressed legal ground energy. Finite-order small defect occupation is supporting evidence, not a replacement for this comparison. The uniform-vacancy modification of C8 is not silently substituted into C32: its second-order equalization would need a new analysis.

## 4. From the lattice through the separated bandwidths

[assumed, R] Use conventional weak-coupling Wilsonian matching and bosonization for a local number-conserving fermionic action with finitely many simple Fermi points. Matching is continuous in the norm in M2, including first momentum derivatives, filling derivatives and the collective irrelevant-vertex norm. Constants may depend on the fixed quadratic dispersion and therefore arbitrarily badly on its velocity ratios. This is the continuum input, not a finite-dimensional theorem.

[proved, finite patch data] For each simple root choose disjoint intervals |k-p_i|<r_i such that the zero eigenvalue branch is isolated from every other eigenvalue by d_*>0 and its derivative has fixed sign and magnitude at least bar_v_*/2. On the compact complement the absolute band energy has minimum e_off>0. Such data exist by simplicity and continuity. A nearest-neighbor n-band symbol has determinant a Laurent polynomial of degree at most n in either direction, so at most 2n roots counting multiplicity unless its determinant vanishes identically. A certified list of 2n simple roots exhausts the zero set.

[proved, root stability criterion] Let M1=sup ||hbar'(k)|| in the root intervals. For a Hermitian perturbation delta h with norm e<d_*/4, the isolated rank-one spectral projector changes in operator norm by at most 2e/d_*. To see this, project the perturbed eigenvector equation onto the complementary eigenspace of the original matrix. Its inverse has norm at most 1/(d_*-e), so the sine of the eigenvector angle is at most e/(d_*-e). The projector trace norm is twice this sine. Hellmann-Feynman then bounds the derivative change by

    ||delta h'|| + 4 M1 e/d_*.

Require this below bar_v_*/4, e<e_off/2 and e<bar_v_* min r_i/4. The perturbed branch remains monotone, has opposite endpoint signs, and has one root in each interval, with displacement at most 4e/bar_v_*. No extra roots occur outside. Taking this displacement below min(delta_2/4,delta_4/8) keeps both margins positive. This supplies the explicit meaning of a small C1 quadratic error. An interacting self-energy must satisfy the same type of bounds under R; finite matrix perturbation theory alone does not prove that property of it.

[proved, energy definitions] A conservative common linearization energy is

    E_lin = W min(e_off/2, min_i(bar_v_* r_i/8)).

Include the length band's corresponding patch data in the mobile case. Let delta be the least dressed mismatch among all dangerous payload and mixed vertices, reduced if necessary to leave a factor-two reserve. Set

    E_det = min(E_lin, W bar_v_* delta/32),
    E_top = W max(1, max_k||hbar(k)||, 4j_l),
    ell_* = max(0, log(E_top/E_det)).                         (5)

All these are strictly positive or finite at fixed n. The ratio E_top/E_det is independent of lambda. In a linear patch, fields with energy at most E have momentum displacement at most 4E/(W bar_v_*), with the factor allowing dressed slopes. Four such fields cannot compensate a mismatch delta if E<W bar_v_* delta/32. This proves the momentum-support cutoff used in (5), with a conservative strict margin.

[heuristic, stages] Above Delta_n the cell and selector spectra are microscopic. Lowering the cutoff through their gaps integrates the massive degrees of freedom using M2. All n retained fermions are already present after this step. Their full lattice dispersions, including tiny hoppings, are kept exactly. Scalar energies and fields are matched at the same time.

[heuristic, stages continued] Order the individual band linearization energies. At energies below the first one, the fastest bands can be represented by chiral fields. Other bands retain their full lattice propagators until their own linearization energies are reached. They are never removed as if a narrow bandwidth were a mass gap. A slow band still has a Fermi surface and low-energy particle-hole states even when the running cutoff is above its whole bandwidth. Its large susceptibility must be included in matching.

[heuristic, control during the stages] Integrating modes above E_lin has a finite infrared cutoff, so its matching constants are finite at fixed n. They can contain powers of E_top/E_lin as well as logarithms; a claim of uniformly order-one constants would be unjustified. Denote their complete contribution by A_band,n. Use the full free lattice propagators during this matching, rather than expanding in the small hopping of a slow flavor. This avoids treating its entire code as unsplit. Below E_lin every branch has a controlled linear patch. The remaining near-resonant quartics run only for the finite interval down to E_det.

[heuristic, a conservative combined bound] Let x0 be a collective dimensionless norm of all marginal deviations after the finite band matching, including forward Gaussian couplings as well as oscillatory non-Gaussian quartics, and Y0 the corresponding irrelevant norm. Relative Gaussian deformation and filling shifts have their own budgets. R gives a bound of the form

    x0,Y0,eta_initial <= A_match,n epsilon/bar_v_*,          (6)

where epsilon measures the normalized residual action, A_match,n includes A_band,n and all multiplicities, and eta_initial is the relative stiffness error. This bound may be very loose. Using the entire ell_* of (5) for the subsequent logarithmic estimate overcounts the preliminary matching and is safe. It includes log of the bandwidth ratio, which is absent from a condition using only log(1/delta).

[proved within the RG majorant] In a weak neighborhood, suppose the conventional RG equations obey

    x(ell) <= X0 + A int_0^ell x(t)^2 dt,
    X0 = x0 + C (Y0+Y0^2)/sigma,

where sigma>0 is the irrelevance margin below and the second term conservatively bounds integrated irrelevant feedback, including possible linear mixing into Gaussian coefficients. Including forward couplings in x is necessary: they change quartic dimensions and give terms of order eta*x, which would be missed by a majorant for oscillatory couplings alone. The comparison solution is

    x(ell) <= X0/[1-A X0 ell].

Therefore A X0 ell_*<1/2 keeps x below 2X0 until every dangerous oscillatory interaction is cut off. A row-sum or coefficient-sum norm is used here. If g4 bounds one matrix element instead, the initial row-sum can be proportional to n g4. All flavor sums must be included in A_match,n or x0.

[heuristic, already-gapless modes] The same bound controls their self-energies and stiffness shifts throughout the finite matching interval. Quadratic terms cannot annihilate their roots while the C1 margins hold. Non-Gaussian dimension-two terms cannot reach strong coupling before the mismatch cutoff while the majorant holds. Higher vertices give finite drift. Thus the modes that become linear first remain gapless while the remaining bands enter their linear regimes. There is no separate hypothesis that each intermediate effective theory is already a complete CFT.

[proved, why bandwidth-only smallness fails] Take a two-component density stiffness [[1,u],[u,v_s]]. Positivity requires u^2<v_s. With v_s=10^(-12) and u=10^(-5), its determinant is -99*10^(-12), although u is tiny compared with the fast speed one. The correct measure is relative to the full stiffness, here |u|/sqrt(v_s). This exact rational example is recorded in checks.json. It diagnoses a possible failure of a proposed hierarchy, not a failure at parameters satisfying (6).

## 5. Which operators remain after matching

[proved, exact symmetry] Total microscopic U(1) acts on the code as total logical fermion number. Every symmetry-covariant elimination therefore forbids anomalous bilinears psi_i psi_j and their adjoints at every order. This includes third order. Number-conserving interflavor Cooper transfer contains two creations and two annihilations, so it is allowed by U(1). It needs momentum detuning or a separately established irrelevant basin.

[proved, C31 order accounting] C31 uses a different central-switch cell. Its finite-order conversion proof gives

    quadratic: lambda^2 tau^2, lambda^3 tau^2, ...,
    non-Gaussian: lambda^2 tau^4, lambda^3 tau^4, ... .

Two conversions insert at most two data fermion operators, whereas a normal-ordered quartic needs four. The lambda^3 tau^2 term changes bands and roots; it is not a pairing mass or a quartic. C31's numerical favorable centered pair-transfer flow is for n=2,3,4, not an all-n sign theorem. Its rho^2 analogue in C26 uses different reflection symmetries and must not be imported without checking them.

[proved, inherited 23-letter/C32 selection rule] The token excitation parity on alternate ordered cells sends lambda to -lambda and fixes the unperturbed Hamiltonian and marker motion. Thus connected bulk coefficients contain even numbers of port insertions. At each fixed order, sufficiently large rings give the same local coefficients as the line. An odd finite ring can have winding processes at order growing with ring length; it does not generate a local cubic bulk term. A single port insertion with arbitrary marker insertions cannot return to the code, so u lambda vanishes. Pure marker processes are scalar on payload but can interact in the length sector.

[proved, scope of the conversion estimate] The supplied source contractions give a quadratic leading lambda^2 g^2 term and a second-order non-Gaussian remainder O(lambda^2 g^4) for the folded bus. C32 bounds the full higher-conversion remainder by 6g^4/m^5 before restoring physical scales. No all-order g^2 factor for every correction is needed here. A generic higher-port bound divided by the tiny W is allowed to lose 1/g^2 and all propagation factors. Lambda can be made smaller after g is fixed.

[assumed, bosonization dictionary] Let z_i be the signed integer chiral-fermion charges of a local payload primary. U(1) gives sum_i z_i=0. At the free point its dimension and momentum are

    D0(z)=sum_i z_i^2/2,     Q(z)=sum_i z_i p_i mod 2pi.       (7)

For mobile cells add the length harmonic h in Z. Marker conservation forbids a length dual-field vertex. In cell coordinates,

    D0(h,z)=h^2+sum_i z_i^2/2,
    Q(h,z)=sum_i z_i p_i-2pi h nu mod 2pi.                  (8)

The exact identity transport and cyclic gluing of C32/C8 are required for this dictionary. The physical coordinate transformation changes the momenta and field metric together. It does not supply an extra local fractional-charge vertex. Arbitrary gluing would require redoing (8).

[proved conditional on (7)] The nonzero neutral vectors of D0<=2 are e_i-e_j and vectors with two distinct +1 and two distinct -1 entries. The first class is excluded by delta_2, and the second by delta_4. Equal incoming/outgoing multisets have z=0 and give densities and their products. A quartic with a repeated identical field requires point splitting; its nonzero primary charge has D0>=3, or it is a descendant of an already excluded bilinear. The Pauli cost must not be counted twice.

[proved conditional on (8)] The extra possibilities with D0<=2 are h=+/-1,z=0 and h=+/-1,z=e_i-e_j. Choose nu in a fixed interior interval, excluding the finitely many solutions 2pi nu=+/-(p_i-p_j) modulo 2pi. The remaining set contains an open interval with a positive mixed margin. No all-order Diophantine condition is required. Higher pure length harmonics have D0>=4.

[proved, exact count] Write N=2n. Up to conjugation, the nonzero payload primary count through dimension two is

    C(N,2)+3 C(N,4).

The mobile count is

    1+3 C(N,2)+3 C(N,4).

These count charge vectors rather than distinct momentum surfaces or descendants. Derivative descendants of a dimension-one primary share its excluded momentum. Exact enumeration through n=6, at most 3102 oriented payload vectors, agrees with these formulas. This is combinatorics, not evidence for a lattice phase.

[proved conditional on bosonization, forward and exchange] A surviving quartic with equal branch multisets reduces by fermionic anticommutation to a product rho_i rho_j, plus contact terms absorbed in lower-degree coefficients. For distinct branches, an apparent exchange psi_i^dag psi_j psi_j^dag psi_i has this form after normal ordering. Its contribution changes the Gaussian density/current matrices. A different outgoing pair, including R_a L_a to R_b L_b, is not a density interaction; its coefficient oscillates and stops below E_det. Attractive signs are harmless only when the pre-cutoff bound holds. Without detuning an arbitrarily small marginally relevant exchange can generate a gap.

[heuristic] After E_det, all non-Gaussian dimension-at-most-two operators have been removed by symmetry or momentum support. Products of oscillatory vertices can generate resonant operators. Those of degree at most four either are Gaussian or obey the same excluded momentum list; higher ones fall under the dimension bound below. Thus closure of the operator classification is required, not a one-loop inspection of one exchange amplitude.

## 6. Gaussian stability with arbitrary speed ratios

[assumed, continuum normalization] Write the Gaussian energy as

    H_G=(1/(2pi)) integral dx (partial_x X)^T H (partial_x X),
    J=[[0,I],[I,0]].

The free H0 is positive and includes different right/left speeds if needed. In normalized free compact coordinates its scaling metric is I, even with unequal speeds. Cross density-current blocks are allowed.

[proved within the Gaussian theory] The scaling metric is

    M(H)=H^(-1/2) |H^(1/2) J H^(1/2)| H^(-1/2)
        = H^(-1) # (J H J),

where # is the positive matrix geometric mean. A primary of coefficient vector l has dimension l^T M(H) l/4. If

    eta=||H0^(-1/2)(H-H0)H0^(-1/2)||<1,

then matrix order, inverse order reversal and joint monotonicity of # imply

    kappa I <= M(H) <= I/kappa,
    kappa=sqrt((1-eta)/(1+eta)).                            (9)

This is C36's matrix bound, applied to the actual unequal-speed H0. It is invariant under the canonical change between physical and cell coordinates. Positivity alone does not give (9).

[proved conditional] Since sum z_i=0 implies sum z_i^2 is even, (8) is an integer. All nonzero primaries not excluded above have reference dimension at least three. If eta<=1/4, (9) gives

    D >=3 sqrt(3/5)>2,
    sigma=3 sqrt(3/5)-2>0.

The less restrictive threshold is eta<5/13. Derivative descendants add integer dimension and cannot invalidate this bound. Neutral nonlinear density terms of degree at least three have at least dimension three. A collective small irrelevant norm is still required to control the sum of vertices; checking one representative vertex is insufficient.

[heuristic, full bootstrap] Split an eta<=1/4 budget among bare matching, finite pre-detuning running and the integrated irrelevant feedback. Require root shifts inside the reserved margins and A X0 ell_*<1/2, with X0 small enough that the complete Gaussian drift fits the budget. Standard weak-coupling RG then gives Y(ell) decaying at a positive rate and finite further renormalization of H. The remaining density interactions are Gaussian, so their nonzero limiting values do not gap modes. This proves the conditional mechanism for saying generated terms only renormalize couplings. It is not a claim that every small U(1)-invariant perturbation is safe without the listed bounds.

## 7. Closing the microscopic choices without exchanging limits

[proved conditional on F,B,M2,R] Choose a leading normalized band witness first, with fixed express energy, loop, phases and normalized field. Keep these shape parameters fixed while reducing g. The source expansion then gives a finite constant C_g,n such that the normalized higher-conversion error is at most C_g,n g^2. Choose g>0 so that this error fits the root, stiffness and RG budgets, including ell_*. If the only available witness is at finite g with no leading-g margin, this step needs a new certificate or a bound that its existing non-Gaussian terms already fit the budget. An isolated finite-g band certificate alone does not prove that decreasing g improves everything.

[proved conditional] Next equalize the mobile second-order bond, including its scalar. C32's bound is

    ||s K2_(p,q)-K2_bare||<=C_n(w/U)^2+E_fast(s/delta_fast).

Its contribution in W units is this bound divided by g^2. Choose nonzero w/U small enough, then nonzero s/delta_fast small enough. Recompute Delta_n using the resulting small charge gap. Choose an interior length density with the mixed margins in section 5. Set u beta_n=W j_l. The scalar depending only on cell number shifts the marker chemical potential and is retuned there; the combination c_22-c_21-c_12+c_11 is a length-density interaction and stays in the error budget.

[proved conditional] Finally choose lambda>0 with every microscopic perturbation small compared with Delta_n and the M2 remainder small compared with W times the reserved infrared tolerance. Even the crude cubic bound gives

    (A_n lambda^3 B_V^3/Delta_n^2)/W
       = O_n(lambda/g^2)                  frozen,
       = O_n(s lambda/g^2)                mobile.

It tends to zero at fixed earlier choices. The even-port selection rule improves the first corrections to lambda^4, lambda^2 u and u^2, but is not needed for existence of this final interval. In particular the normalized terms are bounded schematically by

    [C_4 lambda^4/Delta_n^3
      +C_21 lambda^2 u/Delta_n^2+C_02 u^2/Delta_n]/W,

and all vanish with lambda. Constants contain compensator and field norms. No coefficient is claimed to be order one uniformly in n.

[proved, scope of this schedule] A crossing-preserving prescription such as t_d<=constant*g^2 is not itself an infrared schedule: its momentum margin may shrink with g. The safe order above freezes a nonresonant leading shape before choosing g. It requires the stronger band input B. C32's all-n crossing theorem does not provide B, and enlarging a fixed alphabet without an actual phase-control construction does not provide it either.

[assumed, equilibrium matching] The dressed fillings must change continuously with the action. At the free point the finite-flavor compressibility is nonsingular when every speed is positive; its smallest inverse eigenvalue may be tiny. R must control the energy's filling derivatives so that the implicit-function argument applies. Put the associated inverse-compressibility constant into A_match,n and require the resulting shifts below the margins. Individual flavor fields need not be tunable if these shifts are already small. The mobile density also stays in its chosen interior interval, ruling out phase separation within this weak Gaussian neighborhood.

[proved, limits] All these choices are made at fixed n and are then held fixed. Take the thermodynamic limit followed by the infrared limit. Since E_det>0 for that fixed model, the limit eventually lies below every bandwidth, curvature and detuning scale. No lambda(L), finite-order gap comparison with 1/L, or limit n->infinity before the infrared limit is used. The final crossover time can be much larger than 1/W or any practical time.

## 8. Identification and limits of the conclusion

[heuristic, under F,B,M2,R] The frozen selected sector flows to n compact bosons, and the mobile sector to n+1. Bosonization and Gaussian quantization assign c=1 to each nonchiral oscillator pair. Thus c=n or n+1 in the multi-velocity convention used by the project. There need not be a common Lorentz velocity or a single common-metric CFT.

[proved within the Gaussian theory] In cell circumference M the oscillator vacuum contribution is

    E_Casimir= -pi/(12M) sum_a(v_(a,R)+v_(a,L)).

The entropy coefficient is N/3 for a long interval, N=n or n+1, under the usual compact-Gaussian ground-state identification. A Casimir fit without measured speeds measures a velocity-weighted sum. The cyclic mobile quotient changes allowed zero modes and twists, not the local oscillator count. It does not prove exact finite-ring uniqueness. Shell crossings can persist at particular sizes or parameter values.

[proved, frozen scope] Frozen markers give translated ground sectors and incompatible-ring defects. The conditional c=n refers to a selected thermodynamic broken-translation sector and compatible rings. This is not automatically a unique-ground-state realization on every ring.

[heuristic, acceptance] A careful theoretical physicist should accept the conditional infrared inference with these strict margins and an M2-type local massive-sector matching input. The bandwidth hierarchy introduces severe smallness requirements but no additional fixed-n obstruction. The supplied work does not establish the unconditional claim for every n. It still lacks an all-n nonresonant band witness, and a derivation of M2 for the actual interacting retained background. Some frozen global exclusion bounds also remain unverified. These are located gaps, not evidence that the infrared phase is impossible.
