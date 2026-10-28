# C36. Infrared bounds for the species-graph candidate

[heuristic, verdict] The six-field and fixed-p parts can be made quantitative. The massive-sector argument can be reduced to a static cell problem, and the claimed bounds can be derived for a specified class of defect histories. They have not been derived for the full low-frequency connected action of the microscopic model. The remaining step is a cancellation and summability estimate for the interacting retained background, not the motion of the marker projectors. An unrestricted operator-norm interpretation of the time bound is false. Below I prove this with a local finite-chain example, explain the restriction needed for an infrared bound, and state a sufficient remaining estimate.

[heuristic, acceptance] Conditional on that low-frequency estimate and the complete exact core/compiler, the candidate is acceptable at the stated physics standard using conventional weak-coupling bosonization and RG. The Gaussian phase has an open stability margin, so no unspecified sign of a second-order coefficient remains to be determined. This report does not establish an unconditional microscopic refutation. In particular, it does not relabel the still-missing linked-background estimate as a proved consequence of the defect gap.

[proved, terminology and scope] “Proved” denotes finite algebra, an explicit estimate, or a deduction under named hypotheses. “Assumed” identifies an input not established here. “Heuristic” denotes a continuum inference or assessment. “Numerical” specifies the finite sizes checked. All new files and computations are inside work_C36. No network was used. The project STATE.md was read but was not changed because the task forbids writes outside this folder. PROGRESS.md carries the local memory.

## 1. Inputs and normalization

[assumed, finite microscopic input] Fix p, N=(p+1)(p²+p+1), the exact initialized positive core, the two selected lengths b,b+1, and the wrapper/equalizer. Assume the complete finite-cell kernel and port contracts in the requested notes, including exact species and marker conservation. Their autonomous controller composition is not verified here. The C8 bounds, conditional on this input, give a positive gap for selected-cell excitations and a positive energy proportional to the length of a sufficiently long wrong cell.

[proved, normalization] Subtract the unperturbed energy eL and write

    H(s) = H0 + s V,       P V P = H_ref.

H_ref is the displayed first-order target in cell units, with length velocity 1 and payload velocities tan(k_a) in (1,sqrt(3)). The factors 1/(alpha_L alpha_R), history readout weights, marker overlap, helix compensators, and the sum of physical terms are all included in the fixed norm of V. No visibility factor is discarded. The retained scale is s, not the microscopic norm s||V||.

[proved, conserved field] In MICROSCOPIC_HAMILTONIAN equations (7)-(10), the summed chemical field is

    Z = Q_upper + (zeta-4u) M.

It commutes with H0 and every hopping term. Within a fixed total-charge and marker-number sector it is a scalar. It therefore cancels from virtual energy differences in that sector. It must be restored when comparing ground energies between sectors. This distinction avoids treating the sum of single-site charges in an arbitrarily long cell as a bounded boundary interaction.

## 2. Fixed cell coordinates remove the moving-projector difficulty

### 2.1 Exact rooted representation

[proved, conditional on the stated wrapper] On a legal ring with M>0 markers, mark one of them as the origin and retain its physical position x. Read cells in cyclic marker order. A cell of length l has the finite Hilbert space H_l of its legal selector string and all internal wrapper/core/helix states, with its next marker excluded. Put

    H_cell = direct_sum_(l>=2) H_l,
    h = direct_sum_l [H_open(l) - e l],
    P_cell = P_b + P_(b+1).

Here P_b and P_(b+1) are the complete selected ground projectors, including all logical payload states. For a finite ring, restrict the tensor product of M copies to sum_j l_j=L. Only finitely many summands then occur.

[proved] Concatenation gives a bijection between basis words with a distinguished marker and rooted cell words. Removing the distinguished marker imposes the cyclic relabeling relation, with x shifted by the length of the first cell. This is the quotient already used in C8/joint_fluid. On the cover,

    H0 = sum_j h_j,       [h_i,h_j]=0,
    P = product_j P_cell,j

on the fixed-total-length subspace. The cyclic quotient commutes with these operators. These are fixed projectors on fixed tensor factors. Marker motion changes the length summands inside two adjacent factors; it does not change the definition of P_cell.

[proved] A microscopic marker hop is a partial isometry which transfers the boundary site from one of those factors to the other. A port bond acts on the two boundary records. Both are nearest-neighbor cell operators on the direct sums. Their norms are bounded independently of the lengths: their coefficients are fixed, and each length pair has bounded incoming and outgoing degree under a one-site transfer. All their internal matrices have fixed finite norms. Let B_c bound the sum of these boundary interaction norms incident on one cell. B_c can be extremely large, but it does not depend on M or L. Crossing the chosen root implements the known global twist and changes no bulk estimate.

[proved, gap estimate] Let delta_* be the minimum of the selected internal gaps and the wrong-length coefficient in IR_AUDIT (8). In these coordinates

    h_j >= delta_* Q_j,       Q_j=1-P_cell,j.

For wrong lengths one has the stronger weighted bound in that equation. All h_j commute with their P_j,Q_j. Decomposing by defect patterns sigma gives

    H0|sigma >= delta_* |sigma|,
    ||exp(-tau H0)|sigma|| <= exp(-delta_* |sigma| tau).

This is an operator inequality on each finite ring, not a count of classical configurations. Degeneracy of the logical kernel does not affect it.

[proved, physical locality at fixed order] Starting from selected lengths, r marker insertions change any individual length by at most r. A connected sequence of r boundary insertions touches at most r+1 cells. The total length of this interval initially is at most (r+1)(b+1); transfers across its outer boundaries add at most r sites. Consequently a bound exponential in the number of boundary insertions would give exponential spatial decay in physical units of order b+1. This conversion is finite at fixed p,b. It is not uniform as b tends to infinity.

### 2.2 What follows without an all-order assumption

[proved, second order] Set Q=1-P. For a boundary term v_X, the vector Q v_X P has defects only on X, which contains at most two cells. Since H0 is a sum of independent cell operators, propagation cannot move those defects. In

    H_eff^(2) = -s² sum_(X,Y) P v_X Q H0^(-1) Q v_Y P,

a disjoint X,Y pair vanishes: defects created in Y cannot be removed in X. Overlapping pairs have norm at most s²||v_X||||v_Y||/delta_*. There are a bounded number of such pairs through any cell. Hence, for a local interaction norm defined by summing supports containing that cell,

    ||H_eff^(2)||_loc <= C0 s² B_c²/delta_* .             (1)

The time-resolved version obeys

    ||K_2(tau)||_loc <= C0 s² B_c² exp(-delta_* tau).     (2)

This uses no logical level spacing. Integrating (2) gives (1); inserting tau^n gives n!/delta_*^(n+1). All symmetry generators commute with the projectors, propagators, and v_X, so the effective terms preserve the exact species and marker charges.

[proved, finite order] Degenerate perturbation theory on any fixed finite connected cell interval gives finite sums of products of v_X and reduced resolvents of H0. The only inverted energies are positive sums of cell excitation energies. Folded terms subtract disconnected and reducible pieces. At order r the connected support has bounded size, and the reachable lengths are bounded. Thus for each fixed r there is a finite constant C_r, independent of total ring length, such that

    ||H_eff^(r)||_loc <= C_r s^r B_c^r/delta_*^(r-1).     (3)

This proves neither C_r<=C^r nor convergence of the series. A finite-order statement cannot provide a bound at fixed s after arbitrarily many infrared scales.

## 3. Derivation of the active-defect bounds, and their limitation

### 3.1 A class of histories for which (9)-(10) can be proved

[proved, restricted statement] Consider histories whose every internal perturbation insertion touches at least one currently excited cell. Stop a history when its last defect disappears. Do not include arbitrary retained evolution in this auxiliary propagator. Sum the norms of its blocks between defect patterns. If the input pattern has n defects, at most 2n nearest-cell bonds touch it. Splitting the output on a two-cell bond into P/Q patterns adds at most four blocks. Their norms are at most the norm of the original bond. Increasing a fixed coordination constant c if necessary gives

    sum_(outgoing blocks) ||T_(sigma',sigma)|| <= c B_c n. (4)

The constant does not depend on cell Hilbert-space dimension or on n. Internal h_j evolution is included exactly.

[proved] Let R0 be the time integral of the unperturbed propagator, restricted to nonempty defect patterns. On pattern sigma its norm is at most 1/(delta_* n). In the block-column norm, (4) therefore implies

    ||s T R0||_col <= q_s,       q_s=c |s| B_c/delta_*.

Neumann expansion of this restricted propagator yields, with two endpoint insertions,

    ||K_r^active||_loc
      <= C0 s² B_c²/delta_* q_s^(r-2),       r>=2.       (5)

The compensation of n is an actual column-sum estimate here. An exit insertion can remove the last defects only when n<=2. Converting a source-rooted estimate to the interaction norm through any cell adds a support-marking factor at most proportional to r+1, absorbed by increasing the fixed geometric constant. No unbounded tree-counting factor is hidden in c.

[proved] The corresponding positive majorant in time has loss delta_* n and total gain at most c|s|B_c n in every column. If W(tau) is its total weight, then

    dW/dtau <= -(delta_*-c|s|B_c) sum_sigma |sigma| W_sigma
             <= -(delta_*-c|s|B_c) W.

Consequently

    ||K^active(tau)||_loc
       <= C0 s² B_c² exp[-(delta_*-c|s|B_c) tau].        (6)

The same enlargement of c absorbs support marking in the time bound. Weighting each spatial step by exp(mu) replaces q_s by exp(mu)q_s. For exp(mu)q_s<1 this also bounds exponential spatial moments. Equations (5)-(6) have the form of IR_AUDIT (9)-(10), but their restricted domain is essential.

[numerical, finite majorants] checks.py verifies the column-sum and spectral bounds for positive majorant matrices with 8,32,128 states and integer defect weights from 1 to 7. This checks the implementation of the inequality, not its extension to the actual microscopic action.

### 3.2 Terms not included in that derivation

[proved, limitation] A connected microscopic history need not satisfy the premise of (4). An insertion can create defects on a currently unexcited bond, with those defects joining an earlier branch only later. A previously visited cell can return to P and continue to carry a changed logical state. Retained hopping can spread that state before a later insertion meets an active defect. The number of such available insertions is not bounded by the number of currently excited cells. Counting them with c n in (4) would be incorrect.

[proved, limitation] Purely retained insertions commute neither with a general endpoint source nor with a defect-dependent boundary interaction. Dividing by the retained propagator does not turn their imaginary-time conjugations into contractions. Cancellation of disconnected vacuum pieces does not by itself cancel connected pieces which join the defect history at a later time. The fixed cell representation resolves the moving-projector bookkeeping, but does not prove this additional cancellation.

### 3.3 An exact obstruction to an unrestricted norm interpretation

[proved, finite local example] Let low spins on an open chain of length L have occupation n_j and the dimensionless range-three Hamiltonian

    A_L = sum_(j=1)^L n_j
          + t sum_(j=1)^L n_(j-1) X_j (1-n_(j+1)),
    n_0=1,       n_(L+1)=0,       t>0.

Boundary constants are suppressed tensor identities. The domain-wall subspace spanned by |r>=|1^r 0^(L-r)>, r=0,...,L, is invariant. In this subspace A_L is the Jacobi matrix with diagonal 0,1,...,L and both off-diagonals t. All local term norms are bounded independently of L.

[proved] Write its eigenvalues lambda_j in increasing order. The Rayleigh quotients of |0> and |L> give lambda_min<=0 and lambda_max>=L. Every eigenvector has a nonzero |0> component. If its first component vanished, the first row of the eigenvalue equation would force its second component to vanish, and the tridiagonal recurrence would force all components to vanish. Thus the local operator O=1-n_1, which restricts to |0><0| in this subspace, has nonzero matrix elements between the lowest and highest eigenvectors.

[proved] Append an immobile high state with energy delta>0 and a local coupling

    H = delta |1><1|_high + s A_L + s (|1><0|+|0><1|)_high O.

At s=0 the high sector has the product defect gap delta. Keeping retained propagation s A_L unexpanded and factoring it out gives the following second-order interaction-picture memory kernel, or its adjoint depending on the endpoint convention:

    K_rel(tau) = s² exp(tau s A_L) O
                 exp[-tau(delta+s A_L)] O.

In the domain-wall subspace its matrix element is

    <max|K_rel(tau)|0>
      = s² <max|0> exp[(s lambda_max-delta)tau]
          sum_j |<0|j>|² exp(-s lambda_j tau).          (7)

The leading coefficient is nonzero. Its asymptotic exponential rate is

    s(lambda_max-lambda_min)-delta >= sL-delta.

For any fixed s>0 this is positive at sufficiently large L. No bound of the form C s² exp[-(delta-cs)tau], with constants independent of L, can hold for this interaction-picture operator norm.

[proved] Integrating (7) first at sufficiently small s gives denominators

    delta+s(lambda_j-lambda_max).

The j=min term has a nonzero pole at s=delta/(lambda_max-lambda_min). These poles approach zero as L grows. Thus the Taylor coefficients of this all-retained-energy kernel cannot have a uniform geometric majorant either. This is a statement about the displayed kernel, not a claim that every possible effective-Hamiltonian convention has the same poles.

[proved, distinction from an infrared instability] In the low ground state the same kernel instead has

    <min|K_rel(tau)|min>
       = s² sum_j |<min|O|j>|²
          exp[-(delta+s(lambda_j-lambda_min))tau],

which decays at least as exp(-delta tau). The example demonstrates the need to restrict retained external energies. It does not demonstrate a failure of low-energy massive-sector elimination, and it is not a counterexample to the species-graph Hamiltonian. Boundary conditions and the three-site interaction are immaterial to this specific logical obstruction, whose premise was only local boundedness and a defect cost.

[numerical, sizes] For t=0.4, delta=1, s=0.1, checks.py diagonalizes Jacobi matrices of dimensions 5,9,17,33, corresponding to L=4,8,16,32. The rates in (7) are approximately -0.57017,-0.17017,0.62983,2.22983. The nonzero prefactors become small; the recurrence above, rather than their floating-point size, proves nonvanishing. These calculations concern the diagnostic example only.

[proved, a corrected low-energy bound for the example] Let Pi_Omega project onto retained energies s lambda_j within Omega of s lambda_min. Since Pi_Omega commutes with A_L, the separate exponential bounds give

    ||Pi_Omega K_rel(tau) Pi_Omega||
       <= s² exp[-(delta-Omega)tau].

Indeed the left exponential is bounded by exp[(s lambda_min+Omega)tau] on the projected subspace, the middle one by exp[-(delta+s lambda_min)tau], and ||O||=1. Choosing Omega<delta/2 restores the claimed fast decay in this example. The projector is global, so this finite-dimensional estimate alone says nothing about spatial locality of the resulting effective vertices. Both the energy restriction and local connectedness are needed in the microscopic application.

[proved, a second caution about the bare Q projector] Even after restricting external states to low energy, the global bare Q=1-P need not have a positive dressed threshold. Consider independent high two-level systems with h=delta|1><1|+sX at each of L sites, optionally tensored with retained spectator degrees of freedom. Their ground energy is E0=L e_minus, with e_minus=(delta-sqrt(delta²+4s²))/2. If |g> is the one-site ground state, set w=|<0|g>|², so 0<w<1 at s!=0. The product ground vector |G> has squared overlap alpha_L=w^L with the bare vacuum P. The trial vector Q|G>/sqrt(1-alpha_L) gives

    0 <= min spec(Q H Q|Q)-E0
       <= L |e_minus| alpha_L/(1-alpha_L) -> 0.

To verify the numerator, use H|G>=E0|G>, <G|P|G>=alpha_L, and P H P=0, giving <G|Q H Q|G>=E0(1-2alpha_L). Each individual dressed two-level gap remains sqrt(delta²+4s²). A product of on-site rotations integrates these systems out without difficulty. Thus the vanishing global Q threshold diagnoses the inappropriate bare projection, not an instability. In the mobile candidate one must similarly resum the vacuum dressing locally before assuming a dressed fast threshold.

### 3.4 The remaining microscopic input

[assumed, E: low-frequency linked elimination] A sufficient input is a ground-state or low-external-frequency connected expansion, after subtracting reducible retained propagation and vacuum pieces, with local kernels K_X(tau). For one fixed p and some A_p,mu_p>0 independent of L,s, require

    sup_j sum_(X containing j) exp(mu_p diam_cell X)
      integral_0^infinity (1+delta_* tau) ||K_X(tau)||_LE dtau
        <= A_p s² B_c²/delta_* .                       (8)

The norm here is on the low-frequency effective vertices with retained external energies restricted well below delta_*, not the operator norm of K_rel on all retained states. The construction must specify that norm and implement the restriction. It must retain the exact symmetries and the local operator lattice. It must also establish that the excluded sectors do not become a new low-energy phase. Equation (8) is a sufficient assumption, not a result obtained by renaming (5).

[proved, consequence of E] The zero-frequency correction divided by s is bounded by

    epsilon(s) <= A_p s B_c²/delta_* .                  (9)

Moreover |exp(i omega tau)-1|<=|omega|tau gives a first frequency correction bounded by |omega|/delta_* times the same scale. At retained frequencies of order s or smaller, this is further suppressed as s/delta_*. Exponential spatial summability in (8), together with the cell-to-physical conversion above, gives a local low-frequency action on the scale of finitely many cells.

[assumed, stronger sufficient version] One way to establish E would be an irreducible fast spectral representation K_X(tau)=integral exp(-E tau) dmu_X(E), with dressed support E>=delta_d>=delta_*/2 and spatially weighted total variation at most A'_p s² B_c². The spectral measure here is for irreducible fast vertices after the slow subspace is separated; it is not the full microscopic local spectral function, which has gapless contributions. Integrating the representation proves (8) and an exponential counterpart of IR_AUDIT (10). Such a representation, including its total-variation and locality bounds for the interacting slow background, has not been supplied for this model.

[heuristic, located obstruction] The physical argument still needed is therefore a resummation of the histories excluded from (4), or another controlled low-frequency construction establishing E. A bare defect-energy inequality, an ordinary all-band resolvent, and a second-order calculation do not establish E. The direct-sum cell representation makes this a static local problem with a finite selected kernel and bounded boundary perturbations. It is a more specific task than an estimate with moving projectors, but it remains a task.

## 4. Which hexagons are non-oscillating?

### 4.1 Classification without enumerating the full fermion space

[assumed, continuum normalization] Use the usual chiral vertex charges u,v, with z=u+v, m=u-v, and free payload dimension

    Delta_pay = (||m||²+||z||²)/4.

The creation-minus-annihilation vector z obeys Wz=0. For centered bands the payload momentum is m dot k. For a physical resonant length harmonic h, C13's locality dictionary gives

    m dot k - 2 pi h nu = 0 mod 2 pi,
    Delta_ref,res = h² + Delta_pay.                    (10)

This dictionary and bosonization are continuum inputs, not new finite-dimensional theorems.

[proved under this normalization] Decompose a nonzero integral circulation Wz=0 into oriented simple cycles. Girth six gives ||z||_1>=6. Since m_e=z_e modulo 2,

    m_e²+z_e² >= 2|z_e|,
    Delta_pay >= ||z||_1/2 >= 3.

For a genuine six-field flavor conversion equality requires a single primitive hexagon, with alternating z_e=+1,-1 around it and |m_e|=1 on its six edges. The chirality sign m_e need not equal z_e. One independently chooses which chirality supplies each creation or annihilation. Neutral contractions and current descendants are handled by normal ordering and do not evade this bound.

[proved, exact h=0 criterion] Label the edges a=1,...,N in the order used to set

    k_a = pi/4 + a pi/[12(N+1)].

For a hexagon, every k_a lies strictly between pi/4 and pi/3. The absolute sum of six such signed momenta is strictly below 2pi. Thus its zero-mod-2pi condition is ordinary equality to zero. Four positive signs and two negative signs give a value greater than 4pi/4-2pi/3=pi/3; more imbalanced signs are also excluded. Three signs of each kind are necessary. They cancel the pi/4 part, leaving

    sum_(three positive-m edges) a
      = sum_(three negative-m edges) a.                (11)

This is the complete criterion for an undressed, non-oscillating six-field conversion. It depends on the assigned edge labels, not just the abstract graph.

[proved] There is at most one unordered equal-sum partition of six distinct labels into two triples. If two noncomplementary triples of equal half-total sum intersect in two elements, their remaining elements would be equal. If they intersect in one element, their common element would equal the unique element outside their union, since twice the half-total sum is the total. Both contradict distinctness. Thus each resonant hexagon has two choices m,-m for a fixed orientation z. Allowing both orientations of z gives four vertices, paired into two Hermitian-conjugate pairs. Counts of hexagons, chiral assignments, and cosine terms should not be confused.

[proved, explicit example] For cycle labels (1,6,5,8,9,3), take z=(1,-1,1,-1,1,-1) and m=(-1,-1,1,1,-1,1). Then 1+6+9=5+8+3 and m dot k=0. The operator is allowed and has free dimension three. A generic filling shift can detune it, but the stability argument must allow it to remain resonant.

### 4.2 Enumeration and dressed harmonics

[numerical, exact integer checks at stated sizes] checks.py constructs the projective-plane representatives in the order (1,a,b), then (0,1,a), then (0,0,1), and numbers incidence edges in point-major, line-minor order. It enumerates every noncollinear point triple, hence every hexagon. Complete edge-label and equal-partition lists for p=2,3 are in hexagons_p2.json and hexagons_p3.json. No Fock-space diagonalization is involved.

| p | Payload flavors N | Hexagons | Resonant hexagons, h=0 | Chiral assignments for one z orientation |
| --: | --: | --: | --: | --: |
| 2 | 21 | 28 | 9 | 18 |
| 3 | 52 | 234 | 110 | 220 |
| 5 | 186 | 3,875 | 2,290 | 4,580 |
| 7 | 456 | 26,068 | 17,556 | 35,112 |

[proved, enumeration size] Three noncollinear points specify their three joining lines and one hexagon uniquely. Subtracting collinear triples gives

    choose(p²+p+1,3) - (p²+p+1) choose(p+1,3)
      = p³(p+1)(p²+p+1)/6.

The integer count in each row agrees with this formula. For p=2,3 the script also tests all 64 choices of m on each cycle against the original momentum congruence, independently of the equal-partition criterion.

[proved under (10)] Length dressing changes the exact resonance test. In integer form it is

    3(N+1) sum_e m_e + sum_e m_e a_e
          - 24(N+1) h nu = 0 mod 24(N+1).              (12)

Whenever it holds the dimension is 3+h². At nu=1/2, even h reduces to the h=0 congruence and odd h requires payload momentum pi modulo 2pi. Equation (12) is the criterion at any other chosen nu. In all cases with h!=0 the dimension is at least four, so such dressing cannot weaken the six-field irrelevance threshold.

[proved, natural ordering at nu=1/2] In the stated point-major edge ordering there are no odd-h six-cycle resonances. Four positive chirality signs and two negative ones give momentum in (pi/3,5pi/6), three of each give absolute momentum less than pi/4, and six positive signs give momentum in (3pi/2,2pi). Only five positives and one negative could give pi; their constant pi/4 terms already sum to pi, so that would require one edge label to equal the sum of the other five. It cannot. A cycle has two labels from each of three point rows i<j<k. Each row r occupies labels r(p+1)+1 through (r+1)(p+1). The largest label and the other label in row k differ by at most p. The two row-j labels, with j>=1, sum to at least 2p+5, and the two row-i labels sum to at least 3. Hence the other five labels exceed the largest by at least p+8. Any smaller candidate negative label also fails. Reversing signs treats negative pi. Thus at half density the listed hexagons and their even-h dressings exhaust the non-oscillating six-cycle vertices in this ordering.

## 5. Keeping the stiffness above kappa=2/3

[proved, Gaussian matrix statement imported with its proof mechanism] Let G=[[A,C],[C^T,B]] be the positive full energy matrix and J=[[0,I],[I,0]]. In the Gaussian normalization the dimension metric is

    M(G) = G^(-1/2) |G^(1/2) J G^(1/2)| G^(-1/2)
         = G^(-1) # (J G J),
    Delta_G(ell) = ell^T M(G) ell/4.

Here # is the positive matrix geometric mean. Diagonalizing G^(1/2)JG^(1/2) gives independent chiral oscillators and the first expression. If

    eta_G = ||G_ref^(-1/2)(G-G_ref)G_ref^(-1/2)|| < 1,

then inverse order reversal, congruence by J, and joint monotonicity of # give

    kappa M_ref <= M(G) <= M_ref/kappa,
    kappa = sqrt((1-eta_G)/(1+eta_G)).                  (13)

Joint monotonicity follows from operator monotonicity of the square root and symmetry of #. The resolvent integral for sqrt gives the former; the identity Z P^(-1) Z=Q gives the latter. This is the certificate in contrib/digit_compilation/audit/GAUSSIAN_CERTIFICATE.md, rather than a new assumption about interacting coefficients.

[proved under the Gaussian description] Equation (13) implies

    eta_G < 5/13  =>  kappa>2/3,
    eta_G <= 1/4 =>  Delta_cycle >= 3 sqrt(3/5)>2.32379.

The second choice leaves an irrelevance margin sigma_0=3 sqrt(3/5)-2>0.32379. Positive energy alone is insufficient. The bound is on the entire matrix, including length-payload and density-current blocks, and must hold for the renormalized infrared matrix, not only the bare matching matrix.

[proved, coordinates] With a=b+nu, f_alpha=k_alpha/pi, and V0 the diagonal reference velocity matrix,

    T = [[a,0],[-f,I]],
    G_ref = diag(T^T V0 T, T^(-1) V0 T^(-T)).

The canonical congruence S=diag(T^(-1),T^T) sends this to diag(V0,V0). The relative generalized eigenvalues in (13) are invariant under this congruence. Thus the small physical marker entry K_marker=1/a² is not itself a small stability margin. Microscopic corrections must be matched in these normalized coordinates. The reference T must use dressed fillings, or the fillings must be retuned.

[proved under (10),(13)] Every exactly resonant non-Gaussian vertex not excluded below has reference dimension at least three, so its dimension is greater than two. This includes the resonant hexagons. The pure h=2 length umklapp at nu=1/2 has reference dimension four and interacting dimension at least 4 kappa. In a decoupled length notation its separate threshold would be K_length>1/2; the full kappa>2/3 certificate is stronger and remains meaningful when the modes mix.

## 6. RG tolerance and filling shifts

### 6.1 Dangerous channels and the scale at which they stop

[proved under the vertex dictionary] When z=0, m=2u and Delta_pay=||u||². The only nonconstant payload primaries of dimension at most two are one backscattering and a pair on distinct flavors. For the reference momenta,

    min_(a!=b) 2|k_a-k_b| = delta_pay=pi/[6(N+1)].

Single backscattering has momentum in (pi/2,2pi/3); same-sign pairs have momentum in (pi,4pi/3). They are separated from reciprocal vectors by order-one margins. Exact total upper-species U(1) forbids pairing. The only additional dimension-two class is one length backscattering times one payload backscattering. For nu in [5/12,7/12], its mismatch is at least pi/6. Pure h=1 length resonance requires an integer nu and is excluded.

[proved, units] The physical mismatch is the cell mismatch divided by a=b+nu. The effective physical cutoff is of order 1/a, and physical velocities are s a times the cell velocities. Thus the detuning energy is of order s delta_pay and the logarithmic running interval is

    ell_delta = log(Lambda_cell/delta_pay) + O(1),

where Lambda_cell is of order one. There is no extra log b from mixing physical and cell units. Finite constants in matching may still depend badly on b.

### 6.2 Why N g log is the relevant condition

[assumed, weak-coupling RG input] Let g bound a dimensionless individual quartic amplitude at the matching scale, in units of the smallest retained velocity. Bosonization matching of the local action is assumed continuous in the norm of (8). Constants in that matching may depend on p,b and the fixed core. Let F_ab denote the potentially marginal flavor-scattering matrix. One-loop contractions have sums of the form sum_c F_ac F_cb; higher orders are controlled while the same collective norm is small.

[proved within these RG equations] If |F_ab|<=g, the row-sum norm x satisfies x(0)<=C N g. Matrix multiplication is submultiplicative in that norm. A beta-function majorant |dx/dell|<=A x² therefore gives

    x(ell) <= x(0)/(1-A x(0) ell)

up to the first pole of the majorant. There are N² entries but at most N terms in one quadratic matrix contraction. The conservative requirement before the smallest mismatch stops these channels is

    N g log(1/delta_pay) << 1.                         (14)

If g is instead defined as a collective row-sum norm, the explicit N is already included. An entrywise g must not be used in a condition missing this N. The bound is sufficient and need not be sharp for a specially structured scattering matrix.

[heuristic, bounded irrelevant feedback] In the neighborhood eta_G<=1/4, let Y bound the collective norm of the remaining irrelevant vertices, with their finite matching multiplicities included. Standard perturbative RG gives a majorant Y(ell)<=Y0 exp(-sigma_0 ell/2) while the marginal couplings remain small. Its integrated feedback on the Gaussian or marginal coefficients is bounded by A Y0²/sigma_0. This includes six-field conversions. It is finite even for conversions of zero momentum. The matching norm must control the complete irrelevant operator expansion, including its multiplicities, not merely one representative hexagon. It does not require their coefficients to vanish.

[proved within those majorants] Define X0=x(0)+A Y0²/sigma_0. The integral inequality x(ell)<=X0+A integral_0^ell x(t)² dt implies

    x(ell)<=X0/(1-A X0 ell).

Thus A X0 ell_delta<1/2 and a sufficiently small X0 keep the running matrix inside the chosen Gaussian neighborhood until the nonresonant quartics stop. Beyond that scale only Gaussian quadratic terms and irrelevant vertices remain. The extra Gaussian drift from their feedback is finite. This bootstrap explains how the full renormalized kappa can stay above 2/3; positivity at the initial scale would not be enough.

### 6.3 Controlling the equilibrium momenta

[proved, conditional on smooth matching] At the free point the payload inverse compressibilities are pi tan(k_a)>pi. The length energy per physical site, in first-order units, is

    E_ref(nu) = [1-sin(pi nu)]/[pi(b+nu)],
    E_ref''(1/2) = pi/(b+1/2)>0.

These follow by differentiating the free dispersions and MICROSCOPIC_HAMILTONIAN (12). Hence the local equilibrium Hessian is nonsingular. If the matched energy density and its required derivatives change continuously by O(epsilon), the implicit-function theorem bounds the filling shifts by D_p epsilon for a finite constant D_p. Smooth matching is an input; it does not follow just from an extensive operator-norm error.

[proved, sufficient margins] Impose

    max_a |delta k_a| < delta_pay/8,
    |delta nu| < min(1/24, delta_pay/(16 pi)).           (15)

The opposite-sign pair mismatch changes by at most delta_pay/2, so remains nonzero. The mixed mismatch changes by at most 2 max|delta k|+2pi|delta nu|, well below pi/6. The length density stays in the stated interior interval. The small bound on delta nu is conservative. Other finite bounds preserving the same strict mismatches are equally sufficient.

[proved, optional retuning] The bare relation cos k_a=1/(2t_a) gives dk_a/dt_a=1/(2t_a² sin k_a), which is nonzero throughout the reference interval. Together with the nonzero length response to its chemical field, this supplies a nonsingular free tuning map. Under smooth matching one may retune the hoppings and marker chemical potential to fix the dressed fillings. This is optional; (15) suffices for an existence argument.

[heuristic, combined tolerance] Under E and standard weak-coupling matching, one can choose finite constants D_p so that the initial relative Gaussian deformation, collective irrelevant norm, and filling shifts are bounded by D_p epsilon(s), and the quartic entrywise coupling obeys g<=D_p epsilon(s). A sufficient smallness requirement is then simultaneously

    eta_G, including RG feedback, <= 1/4,
    D_p epsilon(s) << delta_pay,
    N D_p epsilon(s) log(1/delta_pay) << 1,
    Y0²/sigma_0 << the remaining Gaussian budget.       (16)

Every multiplicity not represented explicitly by N belongs in D_p or Y0; neither is asserted uniform in p. These conditions tend to zero with s if E holds. They show that kappa can be kept above 2/3 conditionally, without computing the signs of all induced couplings. They do not supply a microscopic value of kappa at any particular s.

## 7. Fixed-p order of limits and acceptance conditions

[proved, consistency of scales conditional on E] Fix p first. Fix the exact core and selected length, all wrapper parameters, complete port weights, and compensator norms. Then delta_*, B_c, A_p, D_p, N and delta_pay are fixed finite quantities with delta_*,delta_pay>0. Equations (9),(16) admit a strictly positive interval of sufficiently small s. For example choose a tolerance epsilon_p inside all margins in (16) and then choose one constant

    0 < s << min(delta_*/(c B_c),
                 epsilon_p delta_*/(A_p B_c²)).         (17)

The constants in this formula can be very unfavorable. No lower bound uniform in p or in core length is needed.

[proved, order of limits] After that choice, s is held fixed. Take L to infinity and then the infrared wavelength to infinity within the resulting model, equivalently examine arbitrarily long distances in its thermodynamic ground state. The detuning crossover length is of order a/delta_pay in physical sites; the corresponding energy is of order s delta_pay. Both are fixed, though possibly enormous. There is no condition sL<<delta_* and no sequence s(L) tending to zero. To compare central charge with local dimension, first define this limit separately for each p; only then examine the family as p increases.

[proved, exact-core qualification] The exact construction has no residual on-cell logical splitting to overwhelm the first-order scale as s decreases. With a genuinely approximate core splitting epsilon_core, an additional term epsilon_core/s appears in the normalized error. Then taking s arbitrarily small need not help. The exact-code contract, or a compatible accuracy interval, is essential to (17).

[assumed, remaining inputs] The complete list of unestablished inputs for an accepted microscopic conclusion is:

1. A complete autonomous initialized core with the stated kernel, gap, charge symmetry, readout maps, and a p-independent finite alphabet prefactor. The task's reviewed finite ingredients do not by themselves verify that composition.
2. The low-frequency linked elimination E, or another estimate implying its small local correction and derivative/locality bounds, for the actual boundary operators. This includes stability against wrong-length, machinery, selector-violation and marker-free sectors. Selector violations are conserved, so this part includes ground-energy comparison between invariant sectors, not just virtual excursions in one legal sector.
3. Continuous weak-coupling matching to the stated compact-boson local operator lattice, including filling derivatives, followed by the standard Gaussian RG stability used in Section 6. This is the ordinary continuum input at the requested physics standard. It includes control of the collective irrelevant-vertex norm, rather than a numerical sign assumption for a particular quartic.
4. The convention that a stable collection of compact bosons with unequal velocities has total c equal to the sum of its mode charges, and that size-dependent zero-mode sectors and the ring gluing are allowed.

[heuristic, conditional conclusion] Once inputs 1 and 2 are established by a controlled physical derivation, input 3 is the usual weak-coupling argument with the explicit strict margins above. Each of the N payload oscillator families and the length oscillator survives. Under input 4 the result is c=N+1 at a fixed nonzero coupling. The ring quotient changes allowed zero modes and twists, not the number of oscillator families. The family then violates c<=d-1 for sufficiently large p by the already-stated finite alphabet bound. A Casimir coefficient alone would measure a velocity-weighted sum and would not test this unweighted c.

[heuristic, final assessment] The current calculation closes the hexagon classification and makes the stiffness, flavor multiplicity, and order of limits quantitative. It also proves the restricted version of the excursion estimate and removes moving projectors as the source of difficulty. It leaves the full linked low-frequency estimate open. A physics-standard acceptance is warranted after that estimate and the core contract are established, without demanding a rigorous interacting lattice-to-CFT theorem. Those requirements have not both been met by this report.

## 8. Reproducibility

[proved, files] Run `python3 checks.py` in this folder. It uses only the Python standard library, NumPy and SciPy, writes checks.json and the two hexagon lists here, and finishes in seconds. The proof of the kernel obstruction and the active-history bounds is in Sections 3.1-3.3, independent of floating-point checks. Source hashes are recorded in source_manifest.json.

[numerical, limits] No complete encoded core was simulated. No microscopic second-order coefficient was evaluated from its compiled cell matrices. No dressed stiffness or microscopic central charge was measured. The largest graph check had 456 flavors and 26,068 hexagons, the largest diagnostic Jacobi matrix had dimension 33, and the largest positive-majorant matrix had dimension 128. These checks do not test assumption E.
