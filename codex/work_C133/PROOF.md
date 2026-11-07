# C133 proof and obstruction audit

[proved, scope] This report does not construct a fixed-d family with unbounded A/M and does not prove a general bound on A/M. It proves restrictions on a specified encoder and on exact quadratic product-block compressions. The remaining microscopic realization problem is stated explicitly. Tags on a paragraph apply to its displayed equations and the ensuing derivation until the next tag.

## 1. Definitions and two distinctions

[assumed] The infrared consists of finitely many unitary conformal sectors for each individual bond. Write A=sum_a c_a v_a and V=4eM, using the inherited microscopic Lieb-Robinson identification. A sector's velocity is in physical site units. The number of sectors may increase along a family of different bonds, while d stays fixed. This is not an assertion of infinitely many sectors in any one member.

[proved, inherited definition] M is the minimum bond half-spread over scalar shifts and h -> h+o tensor I-I tensor o. On a vacuum-compatible sequence,

    E_L=e_inf L-pi A/(6L)+o(L^(-1)),
    Q_L=L(E_(2L)-2E_L)/M -> pi A/(4M).

Here both L and 2L must lie in the prescribed vacuum sequence and the Hamiltonian is fixed as L varies. A nonvacuum ground-sector contribution (2pi/L) sum_a v_a x_a changes this coefficient. In the numerical controls the sequences are specified in NUMERICS.md and CHECKS.md. In a general incommensurate theory the required condition is that the weighted zero-mode dimensions on both sizes tend to zero, not merely that no massive defect is present.

[proved, logical qualification] A uniform bound on all Q_L implies an A/M bound. An A/M bound by itself does not control arbitrary preasymptotic Q_L peaks. Thus the full supremum C(d) in _Cd_target.md is stronger than the asymptotic assertion unless another uniform estimate is supplied. This run attacks the asymptotic assertion, which a counterexample would also refute for C(d).

[proved, logical qualification] Unbounded A/M need not mean an unbounded number of order-M velocities. The finite towers

    c_k=1, v_k=M/k, 1<=k<=N

have A/M=H_N, while for every fixed u>0 their velocity-weighted count satisfies

    N_c(u):=sum_(v_a>=u)c_a <= floor(M/u),

uniformly in N. This is conformal-sector data, not a constructed lattice chain. It is a reason that threshold counting must have an integrable bound at u=0. The fast equal-speed tower in the question is a sufficient counterexample mechanism, not the only one.

## 2. Where route A loses speed

[proved, inherited microscopic inputs] The sources are C39 PROOFS sections 1 and 2, C44 PROOFS sections 1 through 5, and C45 RESULT, especially its normalization and ordered schedule. A cell has n+2 occupied core sites, n data qubits, and physical length b=n+3 or b+1=n+4. Its code dimension is 2^n. Its full physical dimension grows with n. It is the dimension per physical site, d0=4464, that stays fixed.

[proved] Any isometric encoding of n independent qubits in b sites of dimension d obeys

    2^n <= d^b, hence n <= b log_2 d.

This is a capacity bound only. Multiplying a cell velocity by b converts it to a physical velocity. Capacity does not imply that a logical coupling must be exponentially small, or that the physical velocity is divided by b.

[proved, inherited] Let s be the core energy scale and lambda the physical port coefficient. The port changes a vacuum token into an excited token, so PVP=0. The first nonzero payload interaction is second order,

    -(lambda^2/s) P V Q (H_exc/s)^(-1) Q V P.

For the two buses the leading Gaussian hopping matrix is

    T=-W[Gamma(w)+i eta Gamma(w v_bus)],
    W=(2 lambda^2 g^2/s) V_n D^(-p),
    p=2n+4, V_n=q^(2n+2)/S_n^2,
    S_n=sum_(j=0)^(n+1)q^(2j).

The symbol is T exp(ik)+T^dag exp(-ik)+F, where the readout F is onsite in cell coordinates. The leading Gram formula takes D large at the g^2 coefficient, then g small relative to that scale. It is not an identity for a finite, strongly coupled microscopic chain.

[proved, mechanism] The factor lambda^2/s follows from the silent ports. The factor g^2 follows from one conversion in each cell. The power D^(-2n-4) comes from the actual propagation graph. The left source traverses its high endpoint twice, contributing two inverse high energies. The right source traverses its endpoint and n express vertices and returns, contributing 2n+2. C39 proves that no shorter high-vertex path produces the logical hopping. The two endpoint amplitudes supply V_n. None of these factors can be removed from the reference by changing units.

### 2.1 Endpoint visibility bound

[proved] For every q>0,

    V_n <= 1/(n+2)^2.                                      (1)

Indeed set m=n+2. Then

    S_n/q^(m-1)=sum_(j=0)^(m-1)q^(2j-(m-1)) >= m.

The exponents occur in opposite pairs, and x+x^(-1)>=2, with a middle 1 when m is odd. Squaring gives (1). Equality occurs at q=1. For the biased choice q<1 one also has V_n<=q^(2n+2). Removing the bias removes that exponential factor but cannot improve the endpoint factor beyond (1) within this path code.

### 2.2 A uniform bound on the bus Gram trace

[proved] Put B(w)=aI-w adjacency(P_n), delta=a-2w>0, K(t)=exp(-t B), and

    f_j(t)=int_0^t K_(0j)(r)K_(j0)(t-r) dr,
    Gamma_ij=int_0^infinity f_i(t)f_j(t) dt.

The entries of K are nonnegative since exp(t w adjacency) has nonnegative entries. The semigroup identity gives

    sum_j f_j(t)=int_0^t K_(00)(t) dr=t K_(00)(t).

Since B>=delta I,

    sum_j f_j(t)^2 <= [sum_j f_j(t)]^2
                     <= t^2 exp(-2 delta t).

Integrating and using int_0^infinity t^2 exp(-2 delta t)dt=1/(4 delta^3),

    tr Gamma(w) <= 1/(4 delta^3).                          (2)

The same bound holds at w v_bus for 0<v_bus<=1. The Gram matrices are positive, so their nuclear norms equal their traces. Hence

    ||T||_1 <= W(1+eta)/(4 delta^3).                       (3)

The key point is that this trace bound has no factor n. Numerically integrating the defining functions through n=12 verifies the normalization, but positivity and the semigroup identity prove the bound for all n.

### 2.3 Physical velocities and M

[proved, using the inherited free-band trace theorem] For the Gaussian reference, the sum of nonchiral complex-fermion velocities obeys

    A_payload,cell <= 2||T||_1.

This is the trace inequality proved in REPORT 9.71, before replacing ||T||_1 by a bond scale. It allows the onsite matrix F and unequal left/right slopes. The convention is A=(1/2) sum over simple Fermi points |epsilon'(k)|. Thus there is no implicit factor n in (3).

[assumed coordinate dictionary; proved implication] In the frozen reference x=b j, and in the mobile separated reference x=bar_b j+fluctuations with bar_b in [n+3,n+4]. Consequently k_cell=bar_b k_phys and v_phys=bar_b v_cell. This is also the physical Casimir conversion E_cell(N)=eN-pi A_cell/(6N), N=L/bar_b. In the mobile model the joint zero-mode constraints must still be imposed; they do not change this oscillator coefficient. Therefore

    A_payload,phys/M
      <= [bar_b lambda^2 g^2(1+eta)/(s M delta^3)] V_n D^(-p)
      <= [(n+4)/(n+2)^2]
         [lambda^2 g^2(1+eta)/(s M delta^3 D^(2n+4))].     (4)

Each payload speed is bounded by the same right-hand side times M. This is a bound on the largest achievable reference speed in the stated parameter range, not on the optimum among all encodings.

[proved, microscopic norm witness] In the literal C39/C44 bond an enabled express-rail hop has coefficient of modulus s and changes both physical site labels. Neither a scalar shift nor o tensor I-I tensor o changes such a matrix element. If B is any centered representative, ||B|| is at least that matrix element. Hence the gauge-optimized M>=s. This only needs one legal internal selector/tag choice; other positive penalty scales can make M much larger.

[proved, explicit controlled range] Impose the conservative port condition

    lambda sqrt(1+eta) <= kappa s.

This can always be added to the C45 schedule by decreasing lambda after the other parameters have been fixed. The complete compensated port norm and the wrapper gap usually impose much stronger conditions. Equation (4) then becomes

    A_payload,phys/M
      <= kappa^2 g^2 (n+4)
         /[delta^3 (n+2)^2 D^(2n+4)].                     (5)

This explicitly relates channel number, physical cell size and speed in units of M. It even decreases as 1/n if one formally drops the high-bus factor while keeping the same leading formula and bounded a,w. That formal observation is not a justification for using the large-D formula at D=1.

[proved for the separated reference] C45 chooses the length hopping J_l=W, so its free hard-core length mode has cell velocity at most 2W. Adding it gives

    A_ref/M <= kappa^2 g^2 (n+4)/[(n+2)^2 D^(2n+4)]
               * [4/(1+eta)+1/delta^3].                  (6)

An alternative marker schedule could make that one mode faster. Its contribution is still at most 4eM by the inherited per-sector bound, and does not supply unbounded multiplicity.

[assumed S; conditional conclusion] At each fixed n the finite-D, finite-g, wrapper and matching errors can be selected small enough relative to every retained velocity to keep the dressed velocities within, for example, a factor two of the reference. Under that part of S, twice (6) bounds the scheduled member. No unqualified all-coupling inequality is being inferred from the leading coefficient. C45's final requirement

    lambda <= epsilon_n (J0/s)/(16 A_dec,n C_SW)

is even more suppressive because J0 contains the visibility and propagation factors. Lowering lambda cannot repair the relative g^4 interaction, which must first be controlled by lowering g.

[proved diagnosis; heuristic scope] Silent ports, the high-mass propagation path and endpoint visibility force slowness within this design. The numerical exponents quoted in STATUS are not universal, and the conservative nonresonance and remainder budgets are not optimality theorems. The worst bounds can be relaxed in small projected models, as REPORT 9.84 illustrates. That calculation rescales the logical energies for visibility and does not establish a microscopic fast tower. A different nonperturbative encoder is not excluded by (1)-(6). Perturbative elimination alone also is not a universal speed bound: physical distance rescaling can compensate a factor of 1/b in a completely different effective coupling.

### 2.4 A separate stiffness obstruction to speeding up the weak band phase

[proved, conditional on the stiffness dictionary and relative Gaussian matching] The route-A payload fermions each carry one unit of the same physical conserved charge. Suppose, in that band basis, the matched current quadratic form satisfies

    V_J >= (1-epsilon) diag(v_1^0,...,v_n^0),
    A_payload <= (1+epsilon) sum_a v_a^0,

with 0<=epsilon<1. The velocities in these inequalities are physical. For q=(1,...,1),

    D_Q=(1/pi) q^T V_J q >= (1-epsilon) sum_a v_a^0/pi.

Let S_Q=tr_site[(Q-mean Q)^2]. The normalized on-site charge is (Q-mean Q)/sqrt(S_Q), so the proved lattice budget gives D_Q/S_Q<=4M. Combining these facts,

    A_payload <= 4pi S_Q [(1+epsilon)/(1-epsilon)] M.     (6a)

This is an aggregate comparison. It does not require every normal mode to have a nonzero diagonal-charge weight after a change of basis. Such a requirement would incorrectly discard neutral relative combinations of equal-speed bosons.

[proved, finite on-site calculation] In the 46-letter core the multiplicities of charges -1,0,1,2 are 8,19,15,4. Each core letter occurs 96 times in the wrapper and there are 48 neutral hole states. Thus

    tr Q=1440, tr Q^2=3744, d0=4464,
    S_Q=3744-1440^2/4464=101664/31.

Adding the one length boson costs at most 4eM. This yields an n-independent, though very large, bound for any family retaining the stated relative weak-band matching. The local charge range alone would give the weaker S_Q<=9d0/4.

[assumed, scope] The relative Gaussian matching is the kind of quantitative control demanded by C45's smallness schedule and S. It is not implied by existence of an arbitrary interacting Luttinger liquid. A proposed acceleration of route A that remains uniformly close to its charged free-band reference therefore cannot produce an unbounded A/M even if it avoids the particular bus losses. To exploit the U(1) loophole it must leave this regime and suppress the aggregate current weight relative to the sum of velocities. This additional obstruction is conditional on a continuum dictionary, unlike the finite Gram and visibility inequalities.

### 2.5 Attempt to remove the bias and resonantly accelerate the bus

[heuristic construction attempt] One natural modification is q->1, followed by lowering the excited-bus masses until the high-D propagation penalty disappears. The large-D estimate (5) cannot be continued into that regime. There is, however, a different restriction when the modified unperturbed cell still contains the uniform vacuum-path excitation.

[proved, restricted perturbative lemma] Suppose a product-cell construction of size b has a complete auxiliary gap Delta, port norm epsilon_port=lambda||V||, no first-order intercell code coupling, and only a nearest-neighbor second-order payload bond. Its exact second-order bond satisfies

    ||K_eff|| <= epsilon_port^2/Delta.

This follows directly from ||Q H_0^(-1) Q||<=1/Delta in the source-resolvent formula. With epsilon_port<=kappa Delta, its bond scale is at most kappa^2 Delta. Any retained onsite terms can be removed in an onsite interaction picture, which preserves the norms and supports of the intercell bonds. Applying the nearest-neighbor Lieb-Robinson bound there and converting physical distance gives

    v_payload,phys <= 4e b kappa^2 Delta.                 (6b)

This bounds each conformal speed of that effective Hamiltonian even when its second-order bond is interacting. Identification with the microscopic infrared again requires matching.

[proved, application with a stated spectral hypothesis] The unbiased C39 core path at q=1 has first excited vacuum energy

    s gamma_path=2s[1-cos(pi/(n+2))] <= s pi^2/(n+2)^2.

For the frozen product-cell model, or any modified wrapper that retains this excitation without increasing its energy, the complete auxiliary gap obeys Delta<=s gamma_path. Using M>=s and b<=n+4 in (6b),

    v_payload,phys/M <= 4e pi^2 kappa^2(n+4)/(n+2)^2.

Thus for n unit-central-charge payloads their total A/M is bounded by 4e pi^2 kappa^2 n(n+4)/(n+2)^2, even without the high-D suppression. This is a restriction on controlled acceleration through that unbiased path, not an all-wrapper theorem: the retained-excitation hypothesis must be checked for a changed wrapper. A smaller actual wrapper gap only strengthens the inequality.

[heuristic, remaining possibility] This attempt trades the biased endpoint suppression for a path gap of order 1/n^2. Resonant elimination does not help while the port must remain smaller than that gap. A different mediator with a gap of order M/b could evade this particular estimate, since (6b) then permits order-M physical speeds. It would still have to escape the weak charged-band stiffness comparison or employ neutral interacting modes. No such complete encoder was constructed here.

## 3. A first-order boundary lemma

[proved, restricted fixed-d locality lemma] Partition a nearest-neighbor chain into blocks of b sites and choose a product of block isometries. Suppose its compressed Hamiltonian is exactly a number-conserving quadratic fermion Hamiltonian, with arbitrary quadratic onsite terms and only adjacent-block hopping. Let T be that hopping matrix and assume the Fermi points are simple, with no band identically at the Fermi energy. Then

    rank T <= floor((d^2-1)/2),
    c_complex <= rank T.                                 (7)

This holds for any b and number of encoded orbitals. If this exact code is an invariant low-energy sector of the microscopic chain with the inherited local velocity identification, then

    A_code <= 4eM floor((d^2-1)/2).                       (8)

Equation (8) is conditional on that identification and exact quadratic description, not a theorem for generic interacting codes.

[proved] Expand the single microscopic bond crossing the block boundary in the identity and d^2-1 traceless operators on each of the two physical sites. The terms that act on both blocks have an operator-Schmidt rank at most d^2-1. Product compression and subsequent linear projections onto chosen logical operator subspaces cannot increase this rank. A number-conserving hopping term has, in the bases c_i^dag P,c_i P on the left and c_j,c_j^dag on the right, coefficient blocks T and its conjugate, up to invertible sign conventions. Here P is the left-block parity needed for ordinary tensor products. Both bases are linearly independent and traceless. The coefficient rank is 2 rank T. Single-block terms project to zero. This proves the first inequality.

[proved] Write the one-particle symbol as F+zT+z^(-1)T^dag at the Fermi energy. If r=rank T, determinant multilinearity permits at most r columns from T in a nonzero term and at most r from T^dag. Its determinant is therefore a Laurent polynomial with exponents between -r and r. Unless identically zero, it has at most 2r roots on the unit circle, with multiplicity. Each nonchiral complex gapless sector requires two simple Fermi points. This proves the second inequality even if onsite mixing makes many bands disperse. Merely bounding the rank of T without this root argument would not bound the number of dispersive bands.

[proved limitation] This lemma does not extend by counting the number of port insertions in perturbation theory. Already at second order a resolvent produces the kernel 1/(E_i+E_j). For arbitrarily many distinct positive E_i, the Cauchy determinant

    det[1/(E_i+E_j)]
      = prod_(i<j)(E_j-E_i)^2 / prod_(i,j)(E_i+E_j)

is nonzero. Thus a small number of microscopic boundary operators can produce a full-rank effective matrix. Route A exploits a related source-resolvent mechanism. Quantitative control of the singular values, not their mere count, would be needed to extend (7). Also, generic interacting compressions are not quadratic and need not inherit any such Fermi-root description.

## 4. Attempt (a), a hierarchy of emergence scales

[heuristic construction attempt] Assign c_k=1, v_k=v0 of order M and xi_k=b^k. A possible low-temperature entropy profile, as abstract thermodynamic data, is

    s_N(T)=(pi/3) sum_(k=1)^N min(T/v0,1/xi_k).           (9)

It is nondecreasing and has positive heat capacity wherever differentiable. For every finite N its T->0 slope is pi N/(3v0), the desired N-boson slope. Yet at every temperature

    s_N(T) <= pi/[3(b-1)].

Taking b>=1+pi/(3 log d) makes it obey the microscopic entropy capacity. One can smooth the corners while preserving the bound and the small-T slope, and add high-energy entropy to reach log d. Integration df/dT=-s produces a thermodynamically consistent free-energy function. This is not a lattice realization, a partition-function construction or a proof of locality.

[proved consequence of (9)] Entropy capacity and a per-sector speed bound alone do not contradict the desired fast hierarchy. The same formula with T/v_k and v_k=M/k passes the capacity bound and produces the nonsummable velocity tower of section 1. Equal velocities are not the only dangerous case.

[proved in concrete restricted mechanisms] An exactly quadratic nearest-neighbor physical chain is covered by REPORT 9.71, A<=2M, even with arbitrarily late crossover. The C106 family tests the distinction. A direct product of N microscopic layers needs dimension multiplying with N. A first-order quadratic product-block encoder obeys (7). The controlled route-A encoder obeys (5). None realizes the proposed hierarchy.

[heuristic obstruction, not no-go theorem] Assigning xi_k=b^k does not specify a single repeated bond. A periodic hierarchy with an externally imposed period b^N violates one-site translation invariance unless the period is encoded dynamically or blocked, and blocking raises the on-site dimension. A dynamical selector returns to the problem of transferring many logical operators across one physical cut. The missing ingredient is a fixed-d, one-site invariant interacting bond whose matched Hamiltonian carries the desired spectrum with physical speeds preserved. No such bond was obtained here. The absence of a bond is a construction failure, not a contradiction of the hierarchy.

## 5. Attempt (b), small Luttinger weights

[proved within the compact-boson dictionary] Use N bosons with integer local vertex labels m_k,n_k, one diagonal conserved charge Q=sum_k m_k, velocities v_k=v0 and K_k=kappa/k^2. The local dimensions and integer spins are

    x=sum_k [m_k^2/(4K_k)+K_k n_k^2],
    spin=sum_k m_k n_k.

The stiffness in the convention of the question is

    D_Q=(1/pi) sum_k v_k K_k
        <= v0 kappa pi/6,

whereas A=Nv0. Taking kappa<=1/16 and v0=M makes sum vK far below 4pi M for every N. The local vertex lattice remains mutually local. This supplies a continuum counterexample to any inference from the Kohn budget alone, not a fixed-d lattice counterexample.

[heuristic, stability attempt] Small K makes winding vertices relevant if they are translation invariant. At incommensurate momenta choose the finite list of relevant winding combinations away from reciprocal lattice vectors. For each finite N that list is finite because the dimension form is positive definite. Pure charge-transfer vertices must have sum m_k=0 and at least two nonzero integer entries. For K_k<=1/16 they have dimension at least 8. Thus a diagonal U(1) plus suitable momenta need not force a relevant locking term in this continuum ansatz. Forward couplings can be chosen small. This is a bosonization consistency check, not a proof that these momenta, dimensions and coupling tolerances are simultaneously realizable by a fixed-d bond.

[proved] Taking instead v_k=M/k and the same K_k keeps the stiffness proportional to sum 1/k^3 while A/M=H_N. This version survives every fixed positive velocity threshold bound of order M/u. It again lacks a microscopic realization.

[proved, normalization obstruction to a shortcut] Rescaling a conserved charge does not change its normalized stiffness. If Q'=Q/b, then D_(Q')=D_Q/b^2 and tr[(Q'-mean Q')^2]=tr[(Q-mean Q)^2]/b^2. Their ratio is unchanged. Similarly, if a marker boson describes composites of b occupied particles, Q_occ=b Q_marker gives K_occ=b^2 K_marker. All available trace-orthonormal conserved on-site charges must be included. Selecting only the weak marker charge cannot hide a larger stiffness in another conserved charge. This does not rule out small weights for genuinely weakly visible sectors.

[proved, inherited restricted obstruction] For the usual gapless spin-1/2 XXZ branch with Delta=cos gamma, 0<gamma<pi, K=pi/[2(pi-gamma)]>=1/2. Its velocity is pi sin gamma/(2gamma) at exchange one. Tuning Delta along this branch therefore fails to generate K->0 or extra sectors. The numerical scan checks this accessible interacting family. More general fixed-d chains are not covered by that formula. The small-K continuum construction remains an open lattice loophole, as in C104 and C108.

## 6. Attempt (c), neutral and coset sectors

[proved from CFT formulas and the inherited velocity bound] An Ising factor contributes v/2. One unitary minimal-model factor has c<1, and a Z_p parafermion factor has c=2(p-1)/(p+2)<2. Thus increasing p alone cannot make one such factor's cv unbounded when v<=4eM. A stack of many factors or a large-c coset would still be dangerous. Literal stacking multiplies microscopic dimensions; virtual encoding is not covered by that elementary observation.

[proved conditional on its stated dictionary] A coset contained in a full simple WZW theory at the same velocity has c_coset<=c_parent, since the removed unitary Sugawara central charge is nonnegative. If that full parent is matched to a conserved on-site semisimple symmetry, C104 already charges its entire c_parent v to the stiffness budget. This includes the parafermion piece in SU(2)_k. Gapping or retaining only the coset can remove the required current matching, so it cannot automatically be charged to that parent afterward.

[proved, restricted encoder obstruction] Replacing the payload by independent neutral Gaussian Majoranas does not remove the endpoint and virtual-propagation losses if the same bus architecture is used. Changing the exact symmetry label does not change those amplitudes. Equation (7) concerns number-conserving Gaussian complex modes, so a separate Majorana-rank version would be needed for its strongest counting statement; no generic interacting parafermion inference is made from it.

[numerical] The neutral self-dual Ashkin-Teller line at d=4, FZ parafermion chains at d=2 through 6, the spin-1 TB point and ULS point, and a d=4 spin-orbital interpolation were tested with full-ring energies. The Ashkin-Teller tuning changes velocities and compactification data, not an unbounded number of independent Ising factors. The FZ order parameter p equals the microscopic dimension and is bounded by 6 in this scan. The spin-orbital model can reach the c=3 SU(4) point at fixed d=4, but has supplied no N->infinity mechanism. Their finite-size results and known endpoint checks are in NUMERICS.md.

[heuristic, unresolved] A neutral interacting code that evades the controlled bus suppression is still an open possibility. The same is true of a large-c coset not accounted for by any conserved on-site parent. No general neutral-sector multiplicity theorem has been proved.

## 7. Attempt (d), literature candidates

[proved, source audit] The explicit finite-dimensional black-hole chain of Ikhlef, Jacobsen and Saleur, arXiv:1109.1119, is non-Hermitian. It also separates the Virasoro central charge c_BH=2+6/(k-2) from the ground-energy coefficient c_eff=c_BH-24 h_min=2, with h_min=1/[4(k-2)]. The identity is not a normalizable state. Taking k->2 therefore does not give the requested Casimir tower, even before the Hermiticity obstruction. No numerical scan of this excluded bond was run.

[heuristic, source audit] Dorey and Ravanini's staircase models, arXiv:hep-th/9206052 and hep-th/9211115, exhibit successive near-fixed-point regimes. They do not add all plateau central charges into the final infrared theory. For fixed G and k the coset expression G_k times G_l divided by G_(k+l) has central charge less than c(G_k)<=dim G. Taking rank large changes the field content. These continuum constructions supply no fixed-d Hermitian nearest-neighbor fast tower.

[proved, source audit] The standard nearest-neighbor examples collected in Latorre and Sierra, arXiv:2403.17242, include the Ashkin-Teller orbifold, parafermions, SU(N) fundamental chains and spin-S WZW chains. The directly implemented high-rank examples increase d. Their microscopic Hamiltonians also provide the small controls used here. This audit does not treat their original c-d conjecture as a theorem, and does not override the repository's conditional counterexamples.

[heuristic, literature-search limit] No in-scope literature construction was identified by this search. This is not an exhaustive nonexistence theorem. Links and the precise exclusions are collected in LITERATURE.md. Models without a normalizable conformal vacuum must not have c substituted for c_eff in a Casimir test.

## 8. What would close a fixed-d theorem

[proved] Define the central-charge-weighted count N_c(u)=sum_(v_a>=u)c_a. Positivity and v_a<=V give the layer-cake identity

    A=int_0^V N_c(u) du.                                  (10)

For example, the following genuinely microscopic statement would suffice, without any crossover restriction:

[assumed, missing lemma] For every fixed d there are B(d)<infinity and 0<p(d)<1 such that every chain in the class satisfies

    N_c(u) <= B(d) (M/u)^p(d),  0<u<=4eM.                 (11)

[proved implication] Substituting (11) into (10) gives

    A/M <= B(d) (4e)^(1-p(d))/(1-p(d)).                   (12)

Allowing any explicit integrable function of u/M in place of the power would suffice. A bounded count independent of u would incorrectly bound total central charge. A borderline 1/u bound does not suffice, as section 1 demonstrates. Counting factors without weighting by c also misses a single large-c factor.

[proved partial progress toward (11)] Equation (7) proves a finite total count for exact first-order quadratic product-block encodings, and (5) bounds the whole velocity sum for the controlled route-A reference. C104 supplies another restricted estimate for full symmetry-matched WZW factors. None proves (11) for arbitrary interacting neutral or weakly charged modes. The resolvent Cauchy example identifies one precise failure in promoting the first-order lemma to all virtual orders.

[assumed, sharper research target] An extension would have to bound the cumulative contribution of all dressed boundary channels after eliminating arbitrarily long-scale degrees of freedom, by an integrable function of physical speed. It must control their weights or singular values, not just their rank, and must not identify nonzero overlap of a microscopic observable with independent access to every infrared field. No such extension is established here. Proving (11) would close the asymptotic A bound; a uniform control of Q_L through crossover would still be needed for the full finite-size C(d) target.

[proved, relation to lessons 1 through 5] Lesson 1 rules out treating the formal entropy or stiffness towers above as microscopic constructions, or excluding them using crossed-channel positivity alone. The first-order lemma uses a physical d-dimensional cut, which is extra microscopic input. Lesson 2 is respected by using ring ground energies rather than local energy response. Lesson 3 is respected by allowing unbounded xi_k, testing C106, and not using raw low-energy state counts. Lesson 4 is respected by not rebranding a seam estimate equivalent to C(d) as a proof. Lesson 5 is respected by making no eventual-convexity assumption; the XXZ scan includes Delta=-0.9. None of the lemmas asserts a crossover length bounded in terms of d and speed.

## 9. Numerical normalization and estimators

[proved] For the reflection-symmetric bonds in scan.py, let S exchange the two sites and let w(B)=(lambda_max(B)-lambda_min(B))/2. The function w is convex. Since ShS=h and S(o tensor I-I tensor o)S=-(o tensor I-I tensor o),

    w(h+g(o))=w(h-g(o)),
    w(h)<= [w(h+g(o))+w(h-g(o))]/2=w(h+g(o)).

Thus the ungauged half-spread is the minimum over all Hermitian o, not only a diagonal numerical search. The mobile toy lacks this simple symmetry; its M is instead bracketed by gauge-invariant equal-marginal witnesses and an explicit gauge representative, as explained in MOBILE.md.

[proved] A fit y_L=E_L/L=e_inf+a/L^2+b/L^4 identifies A=-6a/pi if the assumed asymptotic form and controlled corrections hold. The doubled pair estimator needs no e_inf fit,

    A_eff(L)/M=(4/pi) Q_L.

A term of order L^(-3) in E_L produces order L^(-2) drift in this estimator, while marginal logarithms can decay more slowly. Fitting small rings does not prove that all modes have reached their conformal regime. The observed ULS and Heisenberg overshoots and the gapped XXZ control illustrate this limitation.

[proved, inherited free formulas used for verification] With h=(XX+YY)/4, E_L=-csc(pi/L) on even rings. With h=-ZZ-(X tensor I+I tensor X)/2, E_L=-2 csc(pi/(2L)). The decoupled Ashkin-Teller point is two copies of the latter. checks.py compares these predictions with independently assembled full-ring sparse matrices. Its separate matrix-free tensor contraction is checked against sparse assembly on three interacting controls, including the periodic seam.

[numerical, work record] During implementation, the summary formatter initially attempted to format an unresolved fit as a number, and the forced matrix-free small-size control initially entered the dense-solver branch. Both errors were corrected. The final summarizer and numerical checks pass. No failing numerical check was retained as a success, and no larger-ring result was generated by the SLURM file in this sandbox.
