# P1. The c–d conjecture

Verdict: not proved in the stated generality. The argument below proves the bound for the full SU(d)-invariant fundamental nearest-neighbor class, uses the supplied result for quasi-free chains, and proves two conditional current-algebra theorems. It also identifies failures in all four proposed general routes. No counterexample is constructed, and none of these failures is evidence that the conjecture is false.

This is a closed-book analysis. No pre-existing files, literature, or numerical data were consulted. Facts explicitly supplied in the problem are inputs, not claims newly proved here. All additional hypotheses needed for a theorem or a continuum calculation are stated below.

Throughout, c means the usual nonchiral central charge, equal to either chiral central charge when there is one sector. For several sectors, C = c_T = sum_a c_a, not the sum of left and right central charges. Lattice spacing and hbar are one. A conformal low-energy description is understood to describe the full low-energy spectrum, with finite, nonzero velocities. This excludes a ferromagnetic quadratic branch coexisting below the purported conformal spectrum.

## 1. An unconditional class in which equality is forced

### Theorem 1. Full fundamental SU(d) symmetry

Let d >= 2 and suppose the bond obeys

    [h, U tensor U] = 0 for every U in SU(d),

where the on-site representation is the defining d-dimensional representation. If the chain satisfies the conformal hypothesis of the problem, then C = d - 1.

Proof. Invariance under SU(d) implies invariance under U(d), since its extra scalar phases act by scalars on two sites. We first show

    h = a I + b P,       a,b real.

For a unit vector u, invariance under diagonal phase rotations in a basis containing u forces h(u tensor u) to be proportional to u tensor u. Unitary transitivity on unit vectors makes this proportionality constant independent of u. The vectors u tensor u span the symmetric subspace, so h is scalar there. Hermiticity makes its orthogonal complement, the antisymmetric subspace, invariant. For two orthonormal vectors u,v, diagonal phases then force u wedge v to be an eigenvector in that complement. Unitary transitivity on orthonormal pairs makes its eigenvalue independent of the pair. Such wedges span the antisymmetric subspace. Thus h is a linear combination of their two projectors, which are (I+P)/2 and (I-P)/2.

If b > 0, H is the ULS Hamiltonian multiplied by b and shifted by aL. The supplied ULS result gives C = d - 1. Multiplication changes velocities, not central charges.

If b = 0, the Hamiltonian is scalar on the entire Hilbert space. It has no nontrivial conformal low-energy spectrum with positive velocity.

If b < 0, write J = -b > 0. The ground energy is

    E0 = L(a-J),

because each P has eigenvalues +/-1 and fully symmetric vectors attain P=1 on every bond. The shifted Hamiltonian is J sum_j (I-P_{j,j+1}). Choose two colors and let |j> have the second color at site j and the first color elsewhere. Direct application of the two bonds adjacent to j gives

    (H-E0)|j> = J(2|j> - |j-1> - |j+1>).

The Fourier state with momentum k=2 pi n/L therefore has energy

    E(k)-E0 = 2J(1-cos k).

For n=1 this is 4J sin^2(pi/L), of order L^(-2), and is positive. Thus the full low-energy spectrum is not a collection of CFT spectra with nonzero velocities, whose nonzero fixed-dimension gaps scale as L^(-1). Only b>0 remains. QED.

The bondwise symmetry hypothesis can be weakened to global SU(d) symmetry of H_L for every L. Averaging the bond over U tensor U then leaves H_L unchanged, and produces an invariant bond to which the proof applies. This also removes symmetry-breaking telescoping terms in a chosen representation of h.

This theorem does not cover a general bond on C^d, or an arbitrary smaller on-site symmetry. In particular, SU(2) symmetry on a spin-s site is full SU(d) symmetry only when s=1/2.

### Other unconditional consequences of the supplied facts

For the quasi-free nearest-neighbor class, supplied fact 2 gives C <= log_2 d <= d-1. The second inequality follows from 2^(d-1) >= d for integer d>=1, proved by induction using 2d >= d+1. The fermionic bound itself is not re-derived here.

The proved classes are closed under independent tensor products at each site. Suppose the local space is tensor_i C^{d_i}, the Hamiltonian is a sum of independent nearest-neighbor Hamiltonians on these factors, and each factor has C_i <= d_i-1. Then its IR tensor product has C=sum_i C_i and

    C <= sum_i(d_i-1) <= product_i d_i - 1 = d-1.

The second inequality follows by expanding product_i[1+(d_i-1)] and retaining its nonnegative linear terms. This statement assumes independent factors, not arbitrary couplings between them.

## 2. Route (a). What on-site density counting proves

### Lemma 2. Rank of commuting on-site charges

Let q_1,...,q_r be mutually commuting Hermitian d by d matrices, linearly independent modulo scalar matrices. Then r <= d-1.

Proof. Simultaneous unitary diagonalization embeds their real span modulo the identity in the traceless diagonal matrices. That vector space has dimension d-1. QED.

If independent uniform charges Q_alpha=sum_j q_alpha(j) commute on the full lattice Hilbert space, their on-site generators commute as well. Indeed,

    [Q_alpha,Q_beta] = sum_j [q_alpha,q_beta](j),

and tracing out all sites but one gives d^(L-1)[q_alpha,q_beta]=0, because the trace of a commutator vanishes. Linear independence modulo scalars can be checked by the same partial-trace argument.

These are statements about exact microscopic operators. Commutation only after projection to low energies does not justify this partial-trace argument. Nor does a charge with a several-site or increasingly extended density define an on-site matrix q_alpha.

### Theorem 3. Complete abelian current sectors

Assume all the following.

1. The IR has r independent holomorphic U(1) currents with positive-definite level matrix K, and the antiholomorphic description has the corresponding nonanomalous central charge.
2. These r currents come from r independent, mutually commuting on-site charges, independent modulo scalars. In a product of CFTs the currents here are counted across the product.
3. The entire holomorphic stress tensor is the abelian Sugawara stress tensor of these currents. There is no residual stress tensor or additional gapless sector.

Then C=r <= d-1.

Proof. Normalize the current OPE by

    J_i(z) J_j(0) ~ K_ij/z^2.

Whitening the positive matrix K gives r currents with unit level. Their Sugawara tensor is

    T_J = (1/2) sum_i :J_i J_i:.

Double contractions in T_J(z)T_J(0) give r/(2z^4), hence central charge r. Hypothesis 3 gives C=r. Lemma 2 gives r<=d-1. QED.

Here is why hypothesis 3 is substantial. In a unitary CFT with weight-one primary currents, the OPE of T_J with each current agrees with that of the full T. Consequently T_rest=T-T_J commutes with the current algebra and is a stress tensor with central charge C-r. The norm of its vacuum level-two state is

    ||L^rest_{-2}|0>||^2 = (C-r)/2 >= 0.

Thus currents alone establish C >= r. Equality requires the residual stress tensor to be null. The elementary counting argument naturally gives a lower contribution to C, not an upper bound on its remainder.

### Two concrete failures of the proposed universal density premise

First, the critical transverse-field Ising chain has no nontrivial uniform continuous on-site charge. In Pauli notation its interaction is proportional to sum_j X_j X_{j+1}, and its nonzero transverse field to sum_j Z_j. Write a candidate generator as q=aI+bX+cY+eZ. Commuting sum_j q_j with the interaction produces distinct two-site Pauli strings, which force c=e=0. Commuting with the field then forces b=0. One-site and two-site strings cannot cancel. Nevertheless its CFT has c=1/2. Hence a nonzero conformal sector does not require a continuous on-site charge.

Second, SU(2)_k has one Cartan current, but

    c = 3k/(k+2) = 1 + 2(k-1)/(k+2).

For k>1 the second term is positive. It is the central charge remaining after removing the Cartan U(1) Sugawara tensor. The spin-1 Takhtajan–Babujian example, k=2, already has c=3/2 with one Cartan direction. Nonabelian currents are not independent bosons: even SU(d)_1 has d^2-1 affine current components while c=d-1.

The missing lemma for this route cannot be “every CFT sector has an on-site U(1) density,” since Ising falsifies it. A viable replacement would need an additive microscopic resource that also counts neutral cosets and minimal-model degrees of freedom, and whose total is at most d-1. No such resource is constructed here. Replacing exact densities by “asymptotically conserved” ones requires an additional theorem controlling their support and their microscopic algebra; that theorem is also absent.

### Theorem 4. Complete reductive current sectors with index ceilings

Suppose a compact reductive Lie algebra

    g = u(1)^r plus g_1 plus ... plus g_m

acts on one site, faithfully at the Lie-algebra level after removing scalar generators, so it embeds in su(d). Let I_i be the levels induced on its simple factors by the single simultaneous embedding in su(d)_1. Assume:

1. The full IR stress tensor is the sum of the Sugawara tensors for this algebra, with positive simple-factor levels k_i and r nondegenerate abelian currents.
2. Every k_i <= I_i.

Then C <= d-1.

Proof. Supplied fact 6 and coset unitarity applied to this simultaneous embedding give

    r + sum_i [I_i dim(g_i)/(I_i+h_i^vee)] <= d-1.

The abelian contribution is r regardless of the positive abelian level matrix. For each simple factor, the derivative of k dim(g)/(k+h^vee) with respect to k is dim(g) h^vee/(k+h^vee)^2 > 0. Replacing I_i by k_i cannot increase the sum. Hypothesis 1 identifies that sum with the full C. QED.

It is not sufficient to bound one visible WZW factor while leaving an uncounted coset. It is also not sufficient to apply several unrelated embeddings using the same d-dimensional space. The product algebra must have the stated simultaneous representation.

There is also an index-free subcase. If hypothesis 1 holds and dim(g)<=d-1, then C<=dim(g)<=d-1, since each simple-factor contribution is less than dim(g_i) and the abelian contribution equals r. For a pure SU(2)_k theory on a spin-s site, s>=3/2, this gives c=3k/(k+2)<3<=2s=d-1. These statements do not cover extra gapless sectors. For smaller representations, the embedding indices listed in fact 6 give the corresponding conditional WZW conclusions once the required level ceilings hold. An embedding of a parton theory does not by itself identify the IR of an arbitrary bond with that theory. No general lattice-to-IR level-ceiling theorem is proved here.

## 3. Route (b). Energy spectral functions and a dimensionless quantity

### Exact lattice identities

Set h_j=h_{j,j+1} and define the energy current

    j_j = i[h_{j-1},h_j].

Only overlapping bonds fail to commute, so

    d h_j/dt = i[H,h_j] = j_j-j_{j+1}.

For H(q)=sum_j exp(-iqj)h_j and J(q)=sum_j exp(-iqj)j_j,

    i[H,H(q)] = (1-exp(iq))J(q).

Let |0> be a ground state and omega_n=E_n-E0>0. Define the positive spectral measure

    S_h(q,omega) = (1/L) sum_{n:omega_n>0}
                   |<n|H(q)|0>|^2 delta(omega-omega_n).

For every positive-energy matrix element the continuity identity implies

    omega^2 S_h(q,omega) = 4 sin^2(q/2) S_j(q,omega).

These identities require no conformal approximation. Ground-state zero modes are excluded from inverse moments below.

The first-moment sum rule is

    [M_1(q)+M_1(-q)]/2
      = <[H(-q),[H,H(q)]]>/(2L),

where M_1(q)=integral_0^infinity omega S_h(q,omega) d omega. This follows by inserting energy eigenstates in the two terms of the double commutator. Locality makes its right-hand side a finite trigonometric polynomial in q for fixed h. It vanishes at q=0, and the symmetrized expression is even, hence is O(q^2). Its coefficients contain microscopic correlators and commutators, not just d.

### CFT normalization and order of limits

The following calculation has an explicit extra hypothesis: the long-wavelength energy density being measured is the sum of the CFT stress-energy densities, without leading derivative improvements or other terms that contribute at the order being extracted. The state is the periodic conformal vacuum. For a sector of velocity v and central charge c, on a circle of length L,

    H_n = (2 pi v/L)(L_n + bar L_{-n}),       n != 0,

up to interchanging the two chiralities according to Fourier convention. For n>=2 only one term creates a state from the vacuum, and the Virasoro commutator gives

    ||H_n|0>||^2 = (2 pi v/L)^2 [c n(n^2-1)/12].

The excitation energy of this state is 2 pi v n/L. Taking the thermodynamic limit first, then q=2 pi n/L to zero, gives the stress contribution

    S_h^stress(q,omega)
       = sum_a [c_a v_a^2 |q|^3/(24 pi)]
                 delta(omega-v_a|q|).

Thus its m-th spectral moment is

    M_m^stress(q) = |q|^(m+3)/(24 pi) sum_a c_a v_a^(m+2).

In particular, the first moment is O(q^4), whereas its exact positive sum rule generally has an O(q^2) ultraviolet contribution. An upper bound of the latter form gives only

    sum_a c_a v_a^3 <= constant/q^2,

at a conformal momentum q. Taking q to zero makes this bound vacuous; choosing its largest valid value introduces the unknown conformal window. In special models the q^2 coefficient can vanish. A resulting q^4 bound controls sum_a c_a v_a^3 instead, and still needs an additional estimate on velocities or another microscopic inequality to bound C.

An equal-time local variance has an integrated stress contribution proportional to sum_a c_a v_a^2 times the fourth power of a small momentum cutoff. Higher positive moments are likewise suppressed by positive powers of that cutoff. These observations explain why those particular inequalities do not remove the window. They do not constitute a no-go theorem for every possible spectral argument.

### A velocity-free inverse moment

There is a dimensionless observable with a better normalization. Define

    F_h(q) = integral_0^infinity S_h(q,omega)/omega^2 d omega.

For the pure stress contribution,

    F_h^stress(q) = C |q|/(24 pi).

All velocities cancel, even if they differ between sectors. Under H -> lambda H, lambda>0, the matrix elements and gaps scale by the same factor, so F_h is unchanged.

At fixed integer mode n>=2 on a large circle, the more precise formula is

    F_h^stress(n,L) = C(n-1/n)/(12L).

For n=1 the vacuum stress contribution vanishes. Therefore one must not use the infinite-line coefficient at fixed n without its factor 1-1/n^2. This is a useful finite-size consistency check.

To infer F_h(q)/|q| -> C/(24 pi) for the full spectral measure, an additional remainder estimate is needed. One sufficient condition is that the non-stress contribution and cross terms are o(|q|) in this inverse moment. For a fixed high-energy cutoff Omega, the exact continuity identity gives

    integral_Omega^infinity S_h/omega^2
       = 4 sin^2(q/2) integral_Omega^infinity S_j/omega^4.

If the integral on the right stays bounded as q->0, that high-energy part is O(q^2). This conditional estimate does not control additional low-energy operators or improvements.

### An exact obstruction: the bond density is not unique

For any Hermitian on-site A and real lambda, replace the bond by

    h'_j = h_j + lambda(A_j-A_{j+1}).

On a periodic chain sum_j h'_j=H, identically for every L. The local dimension, eigenstates, spectrum, velocities, and central charge are unchanged. Yet

    H'(q) = H(q) + lambda(1-exp(iq)) A(q).

For example, writing F_A(q)=(1/L) sum_{n:omega_n>0}|<n|A(q)|0>|^2/omega_n^2, the coefficient of lambda^2 in F_{h'}(q) is |1-exp(iq)|^2 F_A(q)>=0. Therefore F_h is not an invariant of the Hamiltonian.

This is consequential already at leading conformal order. For a conserved U(1) density whose uniform scaling limit contains nonzero left and right currents, a spatial derivative of that density has dimension two. Its inverse-moment contribution is proportional to lambda^2 |q|, of the same order as the stress contribution. The cross term with the vacuum stress state vanishes: in current notation <0|L_n J_{-n}|0>=n<0|J_0|0>=0. The added contribution is positive and can be increased without limit.

ULS supplies a realization at fixed d: take A to be a nonzero Cartan generator. The SU(d)_1 description has its nonzero Cartan current levels. Thus a proposed inequality

    limsup_{q->0} F_h(q)/|q| <= (d-1)/(24 pi)

for every representation of the bond is false, even at saturation of the conjecture. Normalizing h' to operator norm one does not repair it: dividing the entire bond by ||h'|| divides the total Hamiltonian, matrix elements, and gaps by that same number, leaving F_{h'} unchanged. This is a counterexample to that proposed spectral inequality, not to the c–d conjecture.

For current improvements only, their fixed-mode contribution has the form B n/L, with B>=0, whereas the stress contribution is C(n-1/n)/(12L). In that restricted situation the combination

    16L [F_h(2,L)/2 - F_h(1,L)] -> C

eliminates B. This follows by direct substitution. The combination is not a positive spectral sum, and additional lower-dimension derivative improvements need not have this n-dependence. It therefore supplies an extraction formula under the stated restriction, not an upper bound.

One could instead project onto vacuum stress descendants, or minimize over suitably defined density improvements. Neither is currently a microscopic bound. The missing lemma is a dimension-dependent bound on such a gauge-independent inverse response, together with a proof that its IR coefficient is C. A finite on-site Hilbert space alone does not bound susceptibilities: for H(t)=epsilon Z+tX, the ground-state fidelity susceptibility at t=0 is 1/(4 epsilon^2), although the space has dimension two. This example only refutes a dimension-only susceptibility argument; it is not a conformal chain example.

### Casimir energy and moment ratios

Under the usual periodic-vacuum finite-size hypothesis,

    E0(L) = e_infinity L - pi [sum_a c_a v_a]/(6L) + o(L^-1).

Together with stress spectral moments, this determines velocity-weighted sums. Writing A_p=sum_a c_a v_a^p, positivity gives

    A_1^2 <= C A_2,

by Cauchy–Schwarz. Thus A_1^2/A_2 is a lower bound on C, equal to C for one common velocity. It is not a general upper bound. Even in the common-velocity case, identifying a dimensionless ratio does not bound that ratio by d-1. That would require a new microscopic inequality.

Conclusion for route (b): a velocity-free, window-free IR coefficient can be written down. What remains unproved is its microscopic bound, and the raw energy-density version is spoiled by an exact telescoping freedom. Positive-moment locality bounds alone retain the window.

## 4. Route (c). Counting operators and chiral generators

The space of traceless one-site observables has dimension d^2-1. Three separate steps would be required to turn this observation into a central-charge bound.

1. Show that a uniformly bounded set of one-site scaling fields generates the entire relevant IR algebra.
2. Show that this generation respects enough finite-dimensional microscopic algebraic relations after short-distance products are renormalized.
3. Prove that those relations bound the Virasoro central term by d-1.

None of these steps follows from finite local dimension, and the following tests locate the problems.

### Even the stress tensor need not occur in a one-site linear scaling image

Suppose the on-site action of G is irreducible and the ground state and a stress descendant |t> are G-singlets. For a traceless one-site O, group averaging gives

    integral_G U O U^dagger dU = (tr O/d) I = 0.

Singlet invariance then implies <t|O_j|0>=0. The identity has zero matrix element to the orthogonal excited state as well. Thus no one-site observable linearly creates that stress descendant.

This applies to the singlet vacuum sector of the fundamental SU(d) ULS chain, for example along lengths divisible by d. Its nonzero central charge is already an equality example. Products of adjacent observables can create singlets, as the bond energy does. Therefore surjectivity of a one-site linear scaling map onto the stress sector is an untenable missing lemma.

### A generator count does not bound central charge

For N independent unitary free bosons, let T=sum_{i=1}^N T_i. The OPE is

    T(z)T(0) ~ N/(2z^4) + 2T(0)/z^2 + partial T(0)/z.

The chiral subalgebra generated by the single field T is a unitary Virasoro algebra with central charge N. Unitarity follows because its states form a subspace of the unitary N-boson Hilbert space. N is arbitrary. Thus no bound on c follows merely from the number of generators of a unitary chiral algebra, or from having one stress field at dimension two.

This example is not a fixed-d lattice construction, and does not assert that the full N-boson operator algebra is generated by T. It isolates the invalid algebraic inference from a small generator count to a small Virasoro central term.

### Increasing support and renormalized products

Operators supported on m sites form a space of dimension d^(2m). The lattice-to-continuum construction permits m to increase as resolution is improved. There is no supplied bound making m a fixed function of a continuum scaling dimension, independent of the conformal crossover scale. Moreover, a single lattice observable generally has an asymptotic expansion in multiple scaling fields. The dimension of the microscopic vector space does not bound the number of subleading fields in those expansions.

Finally, multiplication at the same lattice site is not automatically mapped to a renormalized coincident-point continuum product. Singular OPE terms must be subtracted. Consequently finite matrix identities such as a degree-d characteristic polynomial cannot be carried over as continuum operator identities without a separate theorem.

The missing lemma for a viable version of this route is a quantitative lattice-to-continuum algebra theorem that bounds the central extension, rather than merely the number of generators. It would have to include products, growing supports, and the normalization of the stress tensor. No such theorem is proved here.

## 5. Route (d). A complete twist estimate and its direction

### Lemma 5. An exact O(1/L) variational bound

Assume q is a Hermitian on-site generator with integer eigenvalues and

    [h, q tensor I + I tensor q] = 0.

An overall scalar shift of q is immaterial, so integer eigenvalue differences suffice. Define

    U_+ = exp[(2 pi i/L) sum_{j=1}^L j q_j],   U_- = U_+^dagger,
    D_q = ||[q_2,[q_2,h]]||.

For any normalized ground state |Omega>, both twisted-state excitation expectations obey

    0 <= <Omega|U_+/-^dagger H U_+/-|Omega>-E0
       <= 4 pi^2 D_q/L.

Proof. Set delta=2 pi/L. Charge conservation on a bond removes the common rotation of its two sites. The relative rotation is delta on every bond, including the boundary bond, because exp(2 pi i q)=I. Consequently the twisted local term is e^(-i delta q_2)h e^(i delta q_2), with delta replaced by -delta for the other twist.

Taylor's formula with integral remainder gives

    e^(-i delta q_2)h e^(i delta q_2)
       = h + delta i[h,q_2] + R_delta,
    ||R_delta|| <= (delta^2/2) D_q.

The two twists have opposite linear terms. Adding their energy changes therefore gives an upper bound L delta^2 D_q = 4 pi^2 D_q/L. Each change is nonnegative by the ground-state variational principle, so each is bounded by their sum. QED.

No assumption about vanishing persistent charge current was needed. The constant D_q depends on h and the chosen charge normalization, not just d.

If |Omega> is also an eigenstate of translation T and Q=sum_j q_j, the convention Tq_jT^-1=q_{j+1} gives

    T U_+ T^-1 = U_+ exp(-2 pi i Q/L).

Therefore the twisted state is orthogonal to |Omega> when Q/L is not an integer. The bound then supplies a low-energy state orthogonal to that ground state, or another ground state in a degenerate ground space. This is the familiar lower/existence direction of an LSM argument.

For several commuting on-site U(1)s, Lemma 2 limits their independent number to d-1. However, neither the energy estimate nor that count says that the twisted states exhaust the low-energy Hilbert space. Different twists can share the same momentum or the same low-energy sector. A single compact-boson sector already has many twist/vertex states. Conversely Ising has a conformal spectrum without this continuous twist family.

The missing upper-bound lemma would need to prove that all low-energy conformal degrees of freedom are generated and controlled by these twists, including neutral sectors. That statement is false if it requires only on-site U(1) twists, as the Ising example shows. Allowing other twists introduces new data with no demonstrated d-1 bound.

## 6. Checks against every requested example

These checks use the example central charges supplied in the question. A numerical inequality check is not evidence that a conditional theorem's extra hypotheses hold.

| Example | d and C | C <= d-1 | Implication for the proposed routes |
| :--- | :--- | :--- | :--- |
| ULS | d, d-1 | Equality | Theorem 1 applies. Its d-1 Cartan currents have the entire central charge, so their residual stress tensor has zero norm. Raw spectral bounds must survive arbitrary telescoping changes of its bond, which the proposed F_h bound does not. |
| Ising | 2, 1/2 | Strict | No nontrivial uniform continuous on-site charge is needed. This rules out universal U(1)-density completeness and twist completeness. |
| Three-state Potts | 3, 4/5 | Strict | A current-exhaustion hypothesis cannot be inferred from conformality or the value of d. The proof here does not require or establish such a hypothesis for Potts. |
| Z_Q parafermions | Q, 2(Q-1)/(Q+2), Q>=2 | Strict, since 2/(Q+2)<1 | Neutral or non-bosonic conformal degrees of freedom must be included in a general resource count. For Q=2,3 the values reproduce Ising and Potts. |
| Spin-s Takhtajan–Babujian | 2s+1, 3s/(s+1), s>=1/2 | Gap is s(2s-1)/(s+1)>=0 | Equality only at s=1/2. For s>1/2 its one Cartan direction does not exhaust c. For s>=3/2 the pure SU(2) WZW bound follows from c<3<=2s. |
| t–J | 3, 2 | Equality | Compatible with all the inequalities. The supplied charge alone does not verify the hypotheses of Theorem 3 or 4. |
| Qutrit XX | 3, 2 | Equality | Same conclusion. No unprovided microscopic or current-algebra identification is assumed. |

Other checks are built into the derivations. The energy spectral coefficient agrees with the circle Virasoro norm, gives zero stress response at n=1, and is invariant under positive rescaling of H. Its sensitivity to a telescoping density change is explicit. The ferromagnetic sign of the SU(d) bond gives a quadratic branch, rather than being mistaken for a conformal equality case. The twist proof handles a possible linear current term by using both signs.

## 7. What would close the gap, ranked

1. A gauge-independent microscopic stress-response inequality. The inverse spectral moment supplies a concrete normalization that cancels velocities. The required new theorem would control its stress component, or a rigorously defined improvement-independent replacement, by (d-1)/(24 pi). It must also justify the IR remainder estimate. The raw bond-density version is false. This is the most concrete target found here, but its upper bound is still as strong as the original conjecture once the CFT coefficient is identified.

2. A theorem bounding the full IR current algebra and its neutral remainder from a single on-site representation. In a current-complete class, simultaneous embedding plus level ceilings already proves the result. To become general, this would need both a derivation of the ceilings and a budget for cosets, emergent symmetries, and CFTs without spin-one currents. The exact-on-site-U(1) version is ruled out, so a different microscopic resource is necessary.

3. A quantitative lattice-to-continuum algebra theorem controlling the Virasoro central term. Such a theorem would have to use more than generator number and more than the dimension of a one-site vector space. A support bound uniform in the conformal crossover scale, together with a central-term inequality, would be substantial progress. Neither is supplied by locality alone in the arguments above.

Reverse LSM is ranked below these: its existing variational inequality points in the wrong direction, and the necessary U(1)-twist exhaustiveness fails in an elementary example.

There is no candidate counterexample from this analysis. As an algebraic stress test only, an on-site spin-1 chain flowing to a pure SU(2)_6 theory would have d=3 and c=9/4>2, exceeding the spin-1 embedding ceiling I=4. No such nearest-neighbor Hamiltonian, tuning, or realization is produced or asserted here. Writing down an allowed unitary CFT is not evidence that this lattice class realizes it.

The general inequality remains open in this attack. The proved results are the symmetric and quasi-free classes, their independent tensor-product closure, the conditional current-algebra bounds, and the exact twist estimate. The most useful new diagnostic is the velocity-free stress susceptibility together with its telescoping-density obstruction.
