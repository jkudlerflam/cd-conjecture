# Thermal operator-valued matching for arbitrary cell spaces

Status: **proved abstract vacuum and thermal matching theorem**, independently checked by a separate Codex agent. The thermal theorem requires a finite auxiliary heat trace and a sufficiently small, exponentially summable coefficient norm. The microscopic application remains conditional on instantiating every invariant-sector, charge and ensemble hypothesis. No Gaussian infrared phase or central charge is claimed.

Provenance: Codex `infrared_theorem`, 2026-10-28. The parent and a separate agent independently suggested the same abstract auxiliary representation and weighted-decoration summation. The sources read below supply the strategy; the Grassmann-valued extension is derived here, not asserted to be a published theorem. See [the independent audit](../alternatives/INFINITE_CELL_AUDIT.md). The parent independently checked the thermal semigroup argument. This is AI mathematical review, not formal verification or human refereeing.

## 1. A dimension-independent exact representation

[proved] Let K be any separable graded Hilbert space with positive self-adjoint A and a finite code G F_n satisfying

    ker A = G F_n,       A >= g (1-GG^dagger),       g>0.

There is no finite-dimensional assumption on the excited space QK. Define the abstract graded auxiliary space

    E = C Omega direct_sum QK,
    a Omega=0,          a|_(QK)=A|_(QK),

and embed

    J G|u> = |u>_log tensor Omega,
    J|s> = |0>_log tensor |s>_aux,       s in QK.           (1)

No excited eigenbasis is needed. The grading on E is inherited from QK and Omega is even. Then J is even, a is gapped by g, and

    (I_log tensor a) J = J A.

The local valid subspace is the code with auxiliary vacuum plus the logical vacuum with an excited auxiliary. Its complement consists of nonzero logical states with excited auxiliaries and has energy at least g. All embedded interactions vanish on invalid local states; the cut-bond comparison of FINITE_CELL_MATCHING Section 7 applies whenever its interaction-norm condition holds. For thermal statements an excited heat trace is additionally required; a gap alone need not give one.

[proved] Every bounded operator O supported on r cells has an expansion in normal logical fermion monomials with auxiliary-operator coefficients O_I such that

    sum_I H^(degree I) ||O_I|| <= (8 H^2)^(nr) ||O||,
    H>=1.                                               (2)

Proof: take the 4^(nr) logical matrix blocks, each of operator norm at most ||O||, and expand each logical matrix unit using its finite vacuum projector, which contributes at most 2^(nr) terms of degree at most 2nr. Unlike the dedicated-fermion embedding, (2) has **no dependence on dim E**. The auxiliary vacuum projector remains a bounded operator and is never expanded in auxiliary creation operators.

[proved] For a locally block-diagonal reference D=A+B, let H_log be the compression PBP acting on the logical factors in every auxiliary sector and set

    K_D = J B J^dagger-H_log.

Each auxiliary-operator coefficient in K_D annihilates the auxiliary vacuum from both sides. This follows by applying K_D to every logical state with auxiliary vacuum and using independence of the normal logical matrix-unit basis. Consequently its complete auxiliary vacuum effective action is zero at zero temperature, before any logical fields are integrated.

## 2. Norm and the analytic statement

Work first in a finite spatial interval and an imaginary-time slab of duration T, with product auxiliary vacuum at both time boundaries. The unperturbed Hamiltonian is sum_j a_j. An interaction is a sum over spatial intervals X of bounded auxiliary-operator-valued even logical Grassmann polynomials V_X(t). The total parity, including auxiliary operators, is even. Its coefficients are instantaneous before a normal-order point split. All logical fields are external.

For fixed H>=2 and L>0 define the strong input norm

    v_L = sup_j sum_(X contains j) exp(L |X|)
            sum_I H^(degree I) ||V_(X,I)||.               (3)

It is a bound on the density at each time. The interactions relevant here are time independent apart from their external fields. Finite-range cutoffs and finite auxiliary spectral cutoffs may initially be imposed, but all estimates below are independent of those cutoffs.

The output is

    W[V] = log <Omega| T exp[-integral (sum_j a_j+V(t))dt]
                              |Omega>,                 (4)

with the trivial vacuum propagation understood. This is a Grassmann-valued vacuum functional. Its infinite-slab connected kernels are translation invariant in time; in a homogeneous spatial system they are also translation invariant in space.

Fix a desired output spatial exponent a_out>=0 and put

    L >= 4(a_out+log 256+2),
    tau = 2L/g,          kappa=g/8,
    v_* = 1/(4tau).                                     (5)

[proved theorem, proof below] If v_L<=v_*, (4) has an absolutely convergent connected polymer expansion. Its nonconstant kernels have a pinned-coordinate coefficient norm, with field weight H and weight

    exp(a_out diam_x + kappa diam_t),

bounded by a finite B=B(g,L), independent of spatial volume, slab duration, auxiliary dimension, and any logical bandwidth. Scalar terms obey the corresponding bound per spacetime volume. One may take an inefficient constant such as

    B = 100 exp(2 kappa tau)/tau.                        (6)

The precise numerical prefactor is unimportant; the derivation below retains much more slack.

The same Cauchy argument as in the finite-cell note gives

    ||W[K_D+R]-W[K_D]|| <= (4B/v_*) r                   (7)

provided k+r<=v_*/4 in (3). Here scalar densities are included in the input norm and output pressure density. Since W[K_D]=0 with vacuum time boundaries, (7) is a dimension-independent matching bound proportional to the leakage remainder. Two time moments follow from the exponential time weight after a fixed loss in field weight, exactly as in the finite-cell note.

## 3. The time-block estimate

[proved] Divide the slab into blocks of duration tau. For one block, expand the semigroup by Duhamel insertions, leaving the factors exp(-s sum a_j) intact. These are contractions for s>=0. Never estimate an unrestricted interaction-picture operator exp(sa)Vexp(-sa).

For a finite set S of interval labels X, let T_S be the sum of all words in which every label X in S occurs at least once and no other interaction label occurs. Remove the free evolution outside U_S=union_(X in S) X. If

    v_X = sum_I H^(degree I)||V_(X,I)||,

then the auxiliary-operator/Grassmann coefficient norm satisfies

    ||T_S|| <= product_(X in S) [exp(tau v_X)-1].         (8)

Proof. A word containing m_X insertions of each label has norm bounded by the product of its coefficient norms. Its time-simplex volume is tau^p/p!, with p=sum m_X. There are p!/product m_X! label orderings. Summing over every m_X>=1 gives the product in (8). This also applies to measurable time-dependent coefficients with the same uniform norm bound. The vacuum and excited projectors have norm one; no state-counting factor appears.

For S empty, use free evolution. Outside U_S insert, at each site,

    exp(-tau a_j) = P_j + exp(-tau a_j) Q_j,
    P_j=|Omega_j><Omega_j|.

If J is the set of inactive sites using the Q term, the block term has norm at most

    product_(X in S) z_X times b^|J|,
    z_X=exp(tau v_X)-1,      b=exp(-g tau).               (9)

There are no additional P/Q choices *inside* active supports in this representation: T_S retains the entire active operator. Its interfaces are contracted when adjacent blocks are multiplied. Inactive P/Q choices are exactly the set J.

## 4. Exact polymer factorization

[proved] A configuration is the sequence of active label sets S_t and passive excited sets J_t. Mark the spacetime nodes (t,j) with j in U_(S_t) union J_t. Take connected components in the nearest-neighbor square grid, including horizontal neighbors; an active interval fills every node between its endpoints.

A component is a decorated polymer. Different components are separated by at least one empty grid node. Their weights factor exactly in the auxiliary vacuum amplitude:

- Disjoint active spatial supports factor in a single block. Duhamel interleavings sum to the tensor product of the separate propagators.
- At an empty node the rank-one P_j resets that site's auxiliary state to Omega between neighboring occupied time components.
- Spatially disconnected supports commute in the graded tensor product because the total interaction parity is even.

An intentionally coarse connectivity rule can merge otherwise independent pieces; this does not invalidate factorization between distinct components. It only increases the number of possible polymers in a majorant.

Each vacuum-completed polymer has even external logical parity: its auxiliary matrix element is zero unless its auxiliary parity is even, and the full interaction parity was even. Therefore the surviving polymer weights commute as elements of the even Grassmann algebra. This is essential when taking the logarithm and using a scalar hard-core polymer majorant. Open block operators are graded and are not assumed to commute.

The resulting identity is

    Z[V] = sum_(mutually nonadjacent polymers gamma_i)
                         product_i w(gamma_i).          (10)

A polymer's weight is the corresponding auxiliary vacuum contraction, including every active-block Duhamel resummation. Its norm is bounded by (9) multiplied over its time blocks. Configurations made solely of passive Q propagation have zero vacuum weight, but counting them in the majorant does no harm.

## 5. Unbounded interval supports: weighted decorations

[proved given factorization] The finite-range estimate 2^(R N) for decorations is not needed. Set epsilon=exp(-L) and

    z_tilde_X = epsilon^(-|X|) z_X.

By v_X<=v_L and (3),

    eta := sup_j sum_(X contains j) z_tilde_X
        <= tau exp(tau v_L) v_L < 1                    (11)

when v_L<=v_*.

For a fixed occupied spacetime support Y of N nodes, the passive sites in each time block are determined by the active intervals: J_t=Y_t minus U_(S_t). Since b=exp(-2L)<=epsilon and sum_(X in S_t)|X|>=|U_(S_t)|, every decoration has norm at most

    epsilon^N product_(t,X in S_t) z_tilde_X.

Summing all possible active sets contained in Y and dropping the requirement that they have exactly that occupied support gives

    sum_(decorations of Y) ||w||
      <= epsilon^N product_(t,X subset Y_t)(1+z_tilde_X)
      <= epsilon^N exp(eta N)
      <= exp[-(L-1)N].                                  (12)

This is the step that prevents a 2^(N^2) interval-decoration loss. It uses the exponential support norm but not a maximum interval length. An active interval's operator and every logical field in it are assigned to its single time block, so |X| occupied nodes suffice; one does not separately count input/output time-boundary nodes. Using twice that support size would also work after strengthening L in the input norm.

## 6. Connected cluster convergence and weights

[proved] A square-grid connected set of N nodes through a given node can be encoded by a depth-first walk of length at most 2N, so at most 16^N such supports occur. For a polymer or connected cluster of occupied nodes, external spatial diameter is at most its node count, and external time diameter is at most tau times its node count plus two block widths. Thus the desired output weight is bounded by

    exp(2 kappa tau) exp[(a_out+kappa tau) N].

Use the standard hard-core connected-graph expansion in the commutative Banach algebra of even Grassmann kernels. The norm is submultiplicative, so its absolute convergence is bounded by the scalar gas with activities equal to the norms. After inserting the output weight and the usual exp(N) cluster-expansion reserve, the rooted activity sum is bounded by

    sum_(N>=1) [16 exp(-L+1+a_out+kappa tau+1)]^N.

The choices (5) put the bracket below 1/16. Polymers are incompatible if their supports overlap or are adjacent. An N-node polymer has at most 5N nodes in its closed nearest-neighbor neighborhood. Therefore its total incompatible weighted activity is at most 5N/15<N, which is the usual rooted tree criterion for absolute convergence of the logarithm. The same bound controls clusters through a specified spacetime node. It proves a bound substantially smaller than (6) after dividing by tau.

The criterion is the Kotecky-Preiss condition, stated for example in [Fernandez and Procacci, Section 2, equation (2.15)](https://webspace.science.uu.nl/~ferna107/papers/ferpro07.pdf), which was opened in this audit. Its application to the coefficient Banach algebra follows by taking the scalar norm majorant. Equivalently it follows directly from the hard-core tree-graph bound and the rooted-tree recursion: if sum_(gamma incompatible gamma0) z_gamma exp(|gamma|)<=|gamma0|, then the sum of absolute connected cluster weights rooted at gamma0 is at most z_gamma0 exp(|gamma0|). No positivity of the original quantum weights is assumed.

## 7. Continuous-time and pinned-kernel issues

[proved] Normal-order the finite logical operator polynomials before introducing external fields. At degree l, fit the point-split logical times into an interval of length epsilon_time, with creation times later than annihilation times. The spacing can be epsilon_time/(l+1), so the added temporal diameter is uniformly bounded even for unbounded spatial support. The active-set estimate uses the integrated coefficient norm and is uniform under this split. One can first mollify the delta kernels and then take the normal-ordered limit. Fixed-order auxiliary coefficients are bounded matrix elements of contraction semigroups and bounded insertions; coincidence sets for distinct integrated vertices have measure zero. The convergent polymer majorant permits the limiting operation.

The absolute cluster estimate initially bounds the *integrated coefficient norm per spacetime volume*. To obtain a pinned-coordinate norm, use translation invariance of the infinite-slab and infinite-volume kernels. For a nonconstant kernel, integrate all its external coordinates, divide by spacetime volume, and change variables to one external coordinate plus all relative coordinates. This is precisely the relative-coordinate integral with the anchor species summed. It bounds the supremum over the anchor species. The same change of variables works with any chosen external leg and has unit Jacobian; the diameter weight is translation invariant. Thus no uncontrolled degree factor is required to pass to the pinned norm.

To justify continuous, rather than merely block-periodic, time translation invariance, allow an arbitrary origin for the block partition and partial first/last blocks of width at most tau. An inactive Q site in an outer partial block gives zero because it meets the vacuum boundary directly; its lack of a full-block decay cost therefore causes no problem. Active-set bounds only improve when a block is shortened. The majorant is uniform in the partition origin. The exact Duhamel functional is independent of that partition. A cluster contained in its interior has the same weight as in any larger slab. A discrepancy must reach a time boundary; its diameter is at least the distance to that boundary. Reserving half the temporal exponential weight gives an exponentially small interior error. Consequently interior limits obtained with any shifted block origin agree, proving continuous time translation invariance. The per-volume bound therefore transfers to their pinned norms. The analogous spatial argument gives the thermodynamic limit. If an application lacks translation invariance, a separately marked local source estimate is needed; that more general statement is not asserted here.

Scalar coefficients have no external coordinate. They retain the per-spacetime-volume bound, rather than being silently included in a nonconstant seminorm. Sources can be introduced as additional local even vertices in the same norm and differentiated inside the uniform analytic domain. No finite-temperature winding trace has been included by this argument.

## 8. Analyticity, vacuum cancellation, and the leakage bound

[proved] Each active-block amplitude is an entire function of its bounded coefficient vertices. The connected majorant is uniform on the ball v_L<=v_*. Therefore W is analytic there as a map into the weighted kernel norm.

Take a vacuum-reducing baseline K_D and a residual R with coefficient norms k,r and k+r<=v_*/4. For 0<=s<=1, the complex disk

    K_D+sR+zR,          |z|<=v_*/(4r),

stays inside v_L<=v_*/2. The output bound B and Cauchy's estimate give derivative norm at most 4Br/v_*. Integrating in s proves (7). Auxiliary vacuum reduction gives W[K_D]=0 exactly, including its scalar part; it is not an estimate that loses the much larger auxiliary interaction scale k.

With C45's fixed-order local SW input, choose its support exponent large enough for (2)-(3), then choose lambda last. The coefficient norm of K_D is O(lambda) and that of the residual is O(lambda^3), with constants independent of auxiliary dimension. The resulting zero-temperature action error would be O(lambda^3) against the payload scale W=c_n lambda^2. This is the route that preserves the static second-order Hamiltonian correctly. Merely expanding the original finite-range action through order lambda^2 would leave a separate retarded-to-static/contact-contraction matching obligation.

## 9. Application restrictions and remaining work

[proved representation; microscopic application conditional] C45's retained code includes the two selected lengths, so use n+1 logical modes in (1)-(3). Its infinite set of excluded cell states is permitted in E. The original paper's boundedness qualification remains essential: in well-formed invariant sectors the port, transport and readout interactions are bounded; positive unbounded readout contributions in illegal invariant sectors stay in A and are not put into (3). Their zero source connection to the code must be maintained.

[scope of the vacuum theorem] The vacuum functional alone does not identify sectors orthogonal to the global auxiliary vacuum. Section 10 resolves this particular issue by proving a full-thermal-trace matching theorem under an explicit heat-trace hypothesis. A bare gap alone would not suffice. Its use for the complete microscopic physical-length ensemble still requires the separately stated sector and quotient arguments.

[unresolved] Separate payload and length charges, their local operator lattice and the cyclic quotient must be preserved. The abstract embedding preserves the grading, but it does not prove that any convenient common Jordan-Wigner convention leaves the reference exactly Gaussian. The audited identity-transport representation remains necessary.

[unresolved after matching] Even a successful proof of the proposed theorem does not establish the multichannel retarded-action RG, the required physical ensemble transfer, or a microscopic CFT/current/stress scaling limit. The conjecture is not solved by this note.

## Source used

[read] Yarotsky, [math-ph/0412040](https://arxiv.org/pdf/math-ph/0412040), allows infinite on-site Hilbert spaces and proves stability and analyticity of product-vacuum phases. Section 2 builds a coarse-time polymer expansion from operator semigroup bounds. This supports the route used here. Its stated theorem does not itself include arbitrary time-dependent external Grassmann coefficients or the pinned effective-action norm; those are the additional proof steps above and must be checked independently.

## 10. Thermal extension with an explicit heat-trace hypothesis

[proved thermal extension] The vacuum-only restriction can be removed if the on-cell auxiliary heat trace is finite at some positive time. Assume

    Tr_E exp(-tau_0 a) < infinity

for one tau_0>0. Since a has one zero eigenvector and a>=g on its complement, dominated convergence lets us choose tau sufficiently large that

    tau >= 2L/g,        Tr_E exp(-tau a) <= e.            (13)

The same tau works for every spectral or length cutoff if its heat trace is bounded by the full trace. This is a substantive hypothesis; the bare auxiliary gap alone is insufficient.

For the thermal version take

    L >= 4(a_out+log 256+4),
    kappa=L/(8tau),     v_*thermal=1/(8tau),
    B_thermal=100 exp(4 kappa tau)/tau.                  (14)

For beta>=2tau, split the time circle into N_t=floor(beta/tau) blocks of widths tau_t in [tau,2tau]. The active-set bound now uses

    z_X = exp(2tau v_X)-1.

For v_L<=v_*thermal,

    eta=sup_j sum_(X contains j) exp(L|X|)z_X
         <=2tau exp(2tau v_L)v_L <1.                    (15)

The thermal functional is the full auxiliary trace, not a vacuum matrix element:

    Z_beta[V] = Tr_aux T exp[-integral_0^beta(A_aux+V(t))dt],
    W_beta[V] = log Z_beta[V].                          (16)

Its free auxiliary scalar may be kept or normalized out. Neither convention changes differences on adding R. Logical odd fields use their usual antiperiodic time identification. The output temporal diameter is measured with the periodic circle metric d_beta(t,s)=min_(m in Z)|t-s+m beta|, and the spatial diameter uses the corresponding line or finite-ring metric. Total interaction vertices are even in the combined logical/auxiliary grading.

To fix signs and the normal-order convention, first use finite spatial and auxiliary spectral cutoffs and a finite Trotter time lattice. Insert fermionic coherent-state resolutions only in the finite logical Fock factors; retain auxiliary operators and the ordinary auxiliary trace. The logical trace imposes the antiperiodic closing condition, and the graded tensor product fixes every sign from odd auxiliary coefficients. The normal-ordered coherent symbols reproduce the point-split vertices used above. The time-block bounds are uniform in the time step: their discrete insertion sums are bounded by the same exponential active-set majorants. Taking the time step to zero gives the semigroup Duhamel coefficients. Their bounds use contraction and heat-trace estimates, not the spectral-cutoff norm of A, so the spectral cutoff may then be removed. This supplies the thermal influence-functional identity; it is not based on naive cyclicity of odd Grassmann-valued operators.

### 10.1 Trace bound for a completed component

[proved] Let S be the spatial support of one periodic-time component and let A_S=sum_(j in S)a_j. Put P_j at every absent node when evaluating its completed component weight. Expand its active block amplitudes into Duhamel words, but do not expand a_j in an eigenbasis. For each free segment write

    exp(-s A_S)=exp(-s A_S/2) exp(-s A_S/2).

Retain one half of **every** segment for a final Schatten-norm estimate. During a block, a passive site j occurs in no active insertion. Thus all its other half-segments commute with every insertion in that block and combine with its Q_j into

    exp(-tau_t a_j/2) Q_j,

whose operator norm is at most exp(-g tau_t/2). This extraction crosses no noncommuting insertion on j. The remaining first-half evolutions on active sites are contractions. All P factors also have norm one.

The retained global second-half segments have total duration beta/2. Generalized Schatten Holder gives, for arbitrary bounded insertions B_i,

    |Tr product_i [exp(-s_i A_S/2) B_i]|
       <= Tr exp(-beta A_S/2) product_i ||B_i||,
    sum_i s_i=beta.                                    (17)

Indeed assign Schatten exponent beta/s_i to exp(-s_i A_S/2), and exponent infinity to every B_i. Their reciprocal exponents sum to one. A zero-duration segment is omitted. The product of semigroup Schatten norms is exactly Tr exp(-beta A_S/2).

Since beta/2>=tau, (13) gives

    Tr exp(-beta A_S/2)
       = product_(j in S) Tr_E exp(-beta a_j/2)
       <= e^|S|.                                      (18)

Consequently a component with N occupied spacetime nodes and B passive-Q nodes obeys the same active-set majorant as before, multiplied by at most

    e^N exp(-g tau B/2).

This bound is independent of auxiliary dimension. It treats noncommuting active insertions by their operator norms, not by a false bound on exp(tA)Vexp(-tA).

### 10.2 Factorization and convergence on the time circle

[proved] Cutting every inactive P_j wire into |Omega_j><Omega_j| separates components exactly as in Section 4. A component crossing the chosen temporal seam is evaluated by a trace on its own spatial support S, with the above P insertions. If several components use the same spatial site at disjoint times, intervening P wires separate their contractions. No trace over the entire spatial volume is used to estimate an individual component.

An odd auxiliary trace component vanishes by parity, so every nonzero completed activity again lies in the even, commuting logical Grassmann algebra. Graded cyclicity together with antiperiodic logical fields makes the exact trace functional continuously translation invariant on the time circle. In a direct coefficient check, moving a vertex of odd logical degree past the other external fields gives a minus sign, canceled by its antiperiodic time identification. Absolute coefficient kernels are translation invariant regardless of this sign convention.

In (12) replace the passive cost b by exp(-g tau/2)<=exp(-L), and include the extra e^N from (18). The sum over decorations of a fixed support is now at most

    exp[-(L-2)N].                                     (19)

The output temporal diameter costs at most exp[2 kappa tau N+4 kappa tau], and 2 kappa tau=L/4. The choices (14) leave enough room for the animal count 16^N, incompatibility neighborhood 5N, the usual exp(N) rooted-tree reserve, **and one additional exp(N) reserve**. The same convergent hard-core cluster expansion therefore proves the full thermal analogue of (6)-(7), uniformly for beta>=2tau, with B_thermal and v_*thermal.

This applies to arbitrary interval supports satisfying (3), and not only to fixed interaction range. The partial trace has retained all logical fields; no logical correlation function or velocity entered the estimate.

### 10.3 Exponentially small thermal baseline

[proved] If K_D annihilates the auxiliary vacuum from both sides on every support, any K_D-only component leaving a complete empty time layer has zero weight: cut there, obtain a vacuum amplitude, and let its rightmost active coefficient act on the vacuum. Thus a nonzero K_D-only periodic component must occupy every one of the N_t time layers. In particular it has at least N_t nodes. This argument does not require every such component to have a simple topological winding loop.

Use the additional exp(N) convergence reserve above. Every surviving baseline cluster contains at least one component of size >=N_t, so its weighted kernel norm, including scalar density, is at most a fixed multiple of

    B_thermal exp(-N_t)
       <= B_thermal exp[-beta/(2tau)].                  (20)

The harmless fixed multiple can be absorbed into the deliberately loose constant 100 in B_thermal. This proves a thermal source-action bound for the baseline, not just its pressure. At zero coupling it includes the free auxiliary heat trace; subtracting that scalar only improves the relevant difference bound.

For k+r<=v_*thermal/4, the Cauchy estimate therefore gives a full-thermal matching statement of the form

    ||W_beta[K_D+R]||
       <= (4 B_thermal/v_*thermal) r
          + B_thermal exp[-beta/(2tau)],                (21)

with scalar coefficients interpreted as density per beta times spatial volume. This is the properly qualified improvement over the vacuum-only statement. Source derivatives are controlled in a smaller analytic ball; low-frequency time moments follow from the exponential temporal weight.

### 10.4 A beta-independent endpoint and periodization

[proved] Let W_infinity[V] denote the infinite-slab vacuum connected kernels. Periodize each relative external time with the appropriate antiperiodic sign for each odd logical field. The periodization operator is bounded from the exponentially weighted line kernel norm to the circular kernel norm: integrating the absolute sum over one period unfolds into the integral over all relative times, and circular diameter never exceeds the diameter of a chosen lift.

Use the same periodic sequence of block widths to construct the infinite-slab expansion. A connected cluster occupying fewer than N_t time layers can be cut on an empty layer and lifted to the line. Its component trace factors into vacuum contractions at that cut, its incompatibility graph is unchanged, and its weight equals the corresponding vacuum-cluster weight with the graded time-periodization signs. This correspondence counts each lift only modulo common translation by beta. The only unmatched clusters are those covering all N_t time layers, or their infinite-line lifts spanning at least one full period. Their total occupied-node count, with cluster multiplicities, is at least N_t.

The extra exp(N) reserve in the convergence estimate therefore proves

    ||W_beta[V]-Per_beta W_infinity[V]||
          <= C B_thermal exp[-beta/(2tau)],              (21a)

for a fixed numerical C absorbed by enlarging B_thermal if needed. The scalar version compares their coefficients per spacetime volume. This establishes a beta-independent initial endpoint for later RG, with exponentially small thermal corrections. It does not compute the later gapless thermal coefficient.

### 10.5 Invalid states and the original cell model

[proved] The abstract representation has infinitely many invalid states when QK is infinite, but their heat trace is explicit:

    z_invalid(beta)=(2^n-1) Tr_Q exp(-beta A),
    z_valid(beta)=2^n+Tr_Q exp(-beta A).

For interaction norm w=sup_j sum_(X contains j)||interaction_X||<g, resolve the enlarged trace into invariant invalid-site sets B_bad. In that sector every interaction touching B_bad is zero. Its trace is z_invalid^|B_bad| times the remaining cut-system trace, whereas the cut all-valid trace has the factor z_valid^|B_bad|. Restoring the deleted interactions changes log Z by at most beta w |B_bad|, since their total norm is at most w |B_bad|. Thus the ratio to the uncut all-valid partition function is at most [exp(beta w)z_invalid/z_valid]^|B_bad|. Sum over invalid subsets and use z_valid>=2^n. This proves

    0<=f_original-f_enlarged
       <=beta^(-1) log[1+(z_invalid(beta)/2^n) exp(beta w)].
                                                               (22)

For beta>=tau_0,

    Tr_Q exp(-beta A)
       <=exp[-g(beta-tau_0)] Tr_Q exp(-tau_0 A).

Hence (22) is exponentially small with rate g-w. The enlarged action is obtained from a **full auxiliary trace**, so it no longer excludes sectors orthogonal to the global auxiliary vacuum. This removes that specific limitation of Section 9 under the explicit heat-trace and interaction-norm conditions. It does not independently prove the remaining physical-length or cyclic-quotient ensemble transfer.

### 10.6 Remaining status

[checked abstract theorem] The thermal semigroup/Holder argument was independently rederived by the parent and checked by the other agent. The independent audit accepts the representation, coefficient norms, connected expansion, pinned-kernel limit, thermal estimates and invalid-sector comparison under the stated hypotheses. The microscopic application is not certified by that abstract verdict.

[unresolved application] C46 supplies a candidate heat-trace bound for the rooted all-length auxiliary space, and C45 supplies the local SW remainder. Their constants must be inserted into (2)-(3), (13)-(14) with the strengthened support exponent chosen before lambda. The illegal-sector readout qualification remains in force. Physical ensemble transfer and the retarded-action multichannel RG/CFT theorem remain separate obligations. No central-charge conclusion is asserted here.


### 10.7 Explicit candidate heat input from C46

[proved conditional on C46's specified rooted-space counting bound] Let gamma be C46's positive linear wrong-length energy coefficient, d0 its physical alphabet, Delta its complete cell gap, and

    tau_0=log(2 d0)/gamma,
    C0=D_fin+b d0^(b-1)+2 d0^(b+2).

At zero intercell norm, C46 Section 2 bounds the excited cell heat trace by

    Tr_Q exp(-beta A) <= C0 exp(-beta Delta),
    beta>=tau_0.

Thus use g=Delta in the abstract theorem and choose

    tau >= max(2L/g, tau_0, g^(-1) log[C0/(e-1)]).

This satisfies (13). All constants are fixed by the chosen cell and wrapper before the final coupling lambda is decreased. With n+1 retained logical modes, a physical support exponent larger than

    L+(n+1) log(8H^2)

pays (2)-(3); choose that exponent before computing the local SW constants. The resulting coefficient norm is O(lambda) for K_D and O(lambda^3) for R, so a fixed positive lambda can meet the analytic smallness conditions and make the matching error arbitrarily small relative to c_n lambda^2. This is a consistent conditional schedule, not a new numerical parameter certificate.

The representation can preserve commuting payload and length charges by assigning each auxiliary state its original charge relative to the code vacuum and keeping the natural code charges on logical states. This does not by itself supply every required local fermionization identity. Also, the physical length operator is unbounded off the code. An arbitrary length-chemical-potential source is **not** included in the bounded-source norm proved here; it needs a relative-bound or ensemble argument of its own.
