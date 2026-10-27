# C8 progress

2026-10-27. Third run finished, nonconclusive. The complete report is `RESULT.md`,
with new proofs in Sections 15-19. All writes stayed in this folder.

## Current best answer and restart point

There is still no verified counterexample, no whole-class c bound, and no
fixed-d construction with unbounded c. Do not redo the earlier joint-code,
second-order, exact-encoding, or frozen-ring calculations. Their proofs
remain in RESULT.md and the supporting notes.

New proved work in this run:

1. `equal_hopping_selector.md`: the d=4 mobile selector at equal hopping
   t=u=a is exactly a free hard-core domain-wall chain with even wall
   number and two boundary twists. K>2a+|mu| excludes every illegal
   sector by an L-independent gap. For |mu|<4a, standard free-fermion
   continuum input gives c=1 on all ring lengths. There is no payload.
2. `padding_transport.md`: any two-port code can be padded to two lengths
   using one extra vacancy symbol, with identical endpoint maps and
   projected transport I/4. This is only an abstract code construction.
3. `vacancy_parent.md`: the exact two-placement cat cannot be the sole
   local ground code. A repeated-bond vacancy model instead preserves
   ANY uniform frustration-free core code exactly in fixed occupied-
   number sectors, with coherent scalar projected transport. Its two
   endpoint strengths approach the same value but differ at finite
   hopping. Its small tunneling gap is essential to scale accounting.
4. `finite_ff_witness.md`: Hilbert basis theorem plus positivity gives a
   finite K_FF(q) detecting uniform open-chain frustration freeness.
   Selecting the longest frustration-free bulk core cannot give
   unbounded length at fixed q. No bound on general CFTs follows.

Latest reproducible checks, all with single-thread BLAS:

- `python3 equal_hopping_selector.py`: full spectra L=5..9, every illegal
  pattern L=4,5,6, physical entropies L=8,12,16,20. The L=20 interval
  estimator gives 1.003981 and the 16-to-20 half-chain estimator 0.992700.
  The independent Casimir fit L=40..160 gives 0.999999704. Entropies equal
  the XX control plus log(2), as proved.
- `python3 padding_transport.py`: both code lengths, both ports, I/4
  transport and four length-pair bonds pass below 5e-16.
- `python3 vacancy_parent.py`: a three-qubit symmetric core (D=4),
  physical lengths 4 and 5, three hopping strengths. Code/port/transport
  errors are below 1e-15, ground splittings below 3e-15. The initial
  unconstrained eigensolver failed a port tolerance at a tiny tunneling
  gap; enforcing the proved even-reflection sector corrected it.
- Python compilation passed for all three new scripts. Earlier checks
  were not rerun. No full large-c microscopic CFT computation was done.

Open problems, ranked:

1. A high-channel uniform core with controlled endpoint maps. For the
   new local vacancy extension, it must be frustration free. The prior
   arbitrary Choi parent does not qualify.
2. Occupied-number selection in both cell lengths. A chemical potential
   alone fails for the explicit ferromagnetic core, and Section 18
   restricts the longest-frustration-free-core mechanism. Boundary or
   frustrated selection remains possible.
3. Volume-uniform effective theory and infrared stability including
   length-dependent port strengths and ring twists. The exact equal-
   hopping selector is not a proof for the small-hopping long-cell phase.
4. Separately, prove the thermodynamic and infrared premises for the
   rank-one d=20992 chain in RESULT.md Section 6.
5. For a bound, control neutral sectors or conformal onset without a
   lower-gap or visibility assumption. None of the new restricted
   theorems does that.

Next action: work on item 1 or 2. A finite-history port-idle interface is
available in `../../contrib/programmable_cell/astra_P1/RESULT.md`, but it
explicitly does NOT prove the translation-invariant embedding claimed
as routine in parts of its parent README. Do not silently assume it.

The remainder is the chronological record. Its initial task list has been superseded by the ranking above.

Read the project memory and the requested README, REPORT, C7, unbounded-c draft, C6 and tools index. No existing C8 progress file was present. All writes stay in this folder, as requested; the project STATE.md is read-only for this run.

Important update from the live contributions: `../../contrib/unbounded_c/G1/G1_result.md` already proves that arbitrary, unrestricted blocks can produce an antiferromagnetic SU(N) swap at second order through a fixed two-dimensional port, even a rank-one bond. First order vanishes and second order has no unwanted lower-order splitting. Its arbitrary internal Hamiltonians are not a nearest-neighbor realization. Do not redo that construction. The parent README has not yet incorporated it and its universal hierarchy language is too strong.

Initial lines of attack, ranked:

1. Explicit simultaneous left/right boundary compression, including the fermion-parity dressing needed for a quadratic Majorana chain. Check whether this is an algebraic obstacle or solely an internal-Hamiltonian realization problem.
2. Extend G1 to identical blocks with distinct left and right ports; identify the actual sign obstruction for a single identical port. This is distinct from its alternating unrestricted-block construction.
3. Investigate local history-state blocks with a mobile selector head and spectator gates. Seek a realizable boundary code, or prove why the simplest one-pass head construction cannot improve channel capacity.
4. Conditional rank bound for general critical chains remains open. The known current completeness and finite-volume hypotheses are already recorded by C7 and will not be re-derived.

Next step: derive a simultaneous-channel dilation and its dimension/scale costs, then test it by small explicit matrices. No CFT or counterexample is claimed from an abstract dilation.

Substantive step 1 completed.

Proved an explicit joint completely positive map from two q-dimensional boundary sites onto a logical space with m even Clifford generators, for every m <= q^2-1. Its images are alpha Gamma_a on the left and alpha i Gamma_a P on the right, with alpha = 1/(4m sqrt(q)). A Stinespring isometry realizes both maps on the SAME code, with environment dimension at most q^2 2^(m/2). Thus the full encoded block fits in q^4 2^(m/2) dimensions. This does not supply a local parent Hamiltonian.

The compressed coupling sum_a F_a(right) F_a(left) gives the genuine free Majorana target sum_a i gamma_(j,a) gamma_(j+1,a), not the generally interacting undressed Gamma.Gamma chain. Its c is m/2 and velocity is 4J. For q=9, m=80, the hypothetical C7 tensor-selector implementation would have d=36 and c=40. This is ONLY resource accounting: no internal uniform h implementing the code exists in this work.

Checks written and run in `c8_checks.py --ed`, outputs `checks.json`. Joint Choi positivity and both compression maps checked at (q,m)=(2,2),(3,4),(3,8), errors <=3.7e-15. The logical target was calibrated using independent covariance-entropy and Casimir estimators on L=40,48,64,80,96,128: c/n = 1.000000235 and 0.999999661. Direct spin-chain ED for n=1,2 and L=4,8 agrees with the exact energy and doubled entropies. These are target checks, NOT counterexample numerics.

Substantive step 2 completed analytically.

The simplest one-way history head does not realize the desired code. For gates U_j on bonds (j,j+1) visited once in increasing order, W_t=U_t...U_1 gives W_t^dagger O_1 W_t = U_1^dagger O_1 U_1 for every t>=1. Hence the left boundary compression on the unrestricted input register acts on only two input sites, independently of block length. Any even family of nonzero scaled Clifford generators there has m <= 4 log_2(q), and more sharply m <= 4 nu_2(q). Numerical random-unitary check at q=2 and five sites passed at 4.5e-15. This excludes the most direct C7 walker-plus-one-pass-gates extension, not general clocks or constrained inputs.

Next: complete identical-block/two-port extension of G1, physically check its simultaneous dilation, and examine what remains of the second-order sign obstruction. Then write the complete proofs and the precise remaining construction problem in RESULT.md.

Substantive step 3 completed.

Extended G1 to a uniform chain of IDENTICAL blocks, each carrying distinct left and right ports. Their excited spectral weights can be prescribed jointly using block-diagonal syndrome weights. Cross terms from two distinct bonds vanish at second order because an outer block is left excited. Therefore the chain's second-order Hamiltonian is a sum of the desired swaps, with no alternating-block requirement and no same-order three-block terms. For the explicit N=3, q=2 multichannel example the single block has dimension 48 and gap 1. Code builds its full H0 and both port operator families, saved in `two_port_N3.npz`. Physical second-order error is 8.4e-16; first-order local compressions vanish exactly. This remains an unrestricted-block result, not a microscopic fixed-d construction.

Also derived the narrower identical-port sign obstruction: for a Hermitian product bond F tensor F and the same internal H0, second order has connected coefficient matrix minus a Gram matrix. It cannot be a positive antiferromagnetic swap. Distinct left/right ports remove this obstruction.

Substantive step 4, proof drafted mentally and being written.

An all-orders symmetry obstruction survives. If a diagonal SU(N) symmetry of identical blocks preserves H0 and the boundary bond, its conjugation normalizes the finite algebra generated by each port's Schmidt operators. That algebra embeds in M_q and its derivation Lie algebra has dimension <=q^2-1. For N>q, simplicity of su(N) forces its action on each boundary algebra to be trivial. The full chain then has independent SU(N) symmetries on every block. If each ground code is a single fundamental representation, every effective operator on the product code is scalar. Thus an exact SU(N)-protected swap with N>q is impossible at EVERY order, even allowing arbitrary block Hamiltonians. This does not exclude unprotected or emergently symmetric targets, or neutral multiplicity spaces.

Remaining work: finish proofs in RESULT.md, audit normalization and assumptions, and rank the unresolved local-code/stability/ring-sector steps. No conclusive outcome under the user's definition has been reached.

Substantive step 5 completed.

A simpler simultaneous code exists for q a power of two. Uniformly and coherently choose one of 2m pairs (boundary, Pauli channel), record the logical Pauli eigenvalue, prepare the corresponding physical Pauli eigenspace state at that boundary and a maximally mixed state at the other, and retain all records. This is an isometry with BOTH desired compressions at scale 1/(2m), larger than the general Choi estimate. It has an explicit finite-circuit description. Its total Hilbert dimension is <=4m q^4 2^(m/2). Example q=16, m=128 gives c=64, hypothetical selector dimension d=64, and fits in 23 spectator sites. The Choi-minimal version fits in 20. No repeated-bond parent Hamiltonian is supplied. Checks at (2,2),(4,8) pass to 6.1e-15.

Also examined target stability beyond the symmetric free point. With n>=9 complex flavors, generic incommensurate Fermi momenta eliminate nonzero-charge relevant and marginal vertices by translation, while total-parity twist operators have free dimension n/4>2. Remaining low-dimension couplings are currents and current products. RESULT.md gives explicit irrational Fermi points and proves their nonresonance. The operator classification assumes standard bosonization; a stable finite-coupling microscopic phase is still an RG inference, not a proved lattice result.

RESULT.md now contains the full proofs, constructions, numerics and limitations. Next is a final audit and a concise handoff. The unresolved core is local implementation of the code (or spectral weights), not its abstract compatibility or second-order rank.

Substantive step 6 completed, conditional rank-conjecture route.

For a target bond of Schmidt rank r, the two-port rank-one gadget needs only L=r+2 excited energies, not N^2+1. With N=2^9=512 and the generic nine-flavor quadratic target, r=18, L=20 and physical site dimension d=20992 suffice. Treating the entire unrestricted block as a physical site gives an actual identical nearest-neighbor bond of connected rank one, explicitly specified in RESULT.md (31). Its second-order logical Hamiltonian is the c=9 target plus a scalar, with the required chemical potentials added as on-site terms.

This would refute c<=r(h) if two missing premises hold: a thermodynamic local effective-Hamiltonian reduction with small relative errors, and stability of the incommensurate c=9 target under those errors. Neither premise has been proved in this run. It is NOT a verified counterexample, and even if validated it would not violate c<=d-1 because d=20992. This is worth pursuing separately because it bypasses the small-d internal-block locality problem.

## Resumed run, 2026-10-27

Read the saved state and the current C9 report. C9 independently reached the joint-code, two-port spectral synthesis, exact large-symmetry obstruction, and missing local-simulation issues. Its work does not supply a repeated-bond realization. The d2 contribution has two provisional numerical candidates awaiting further checks; this run does not duplicate that search.

Substantive step 7, new exact-encoding obstruction (proof being written).

For blocks with disjoint endpoint sites, assume the full Hamiltonian exactly preserves the tensor product of block code spaces. Extract a traceless logical coefficient on the neighboring block from the one-block leakage equation. Every active interaction operator at a boundary then has a Hermitian physical lift that preserves the block code. The left and right lifts commute physically, so their logical restrictions commute. Consequently the left and right interaction algebras of an exactly encoded logical chain must commute. This rules out the Section 2 hopping target as an exact product-code restriction at any port dimension, because [R_a,Gamma_a]=-2iP. It does not rule out perturbative encodings.

For a conventional ladder target with separate left and right Majoranas and an internal mass, this argument instead makes the active Clifford algebra at each port a representation of a subalgebra of M_q. For m=2n port generators, q>=2^n, so c=m/2<=log_2 q in this exact-encoding subclass. Scaled compressions do not avoid that bound once the active operators preserve the code.

Next: write and numerically check this statement and its leakage inequality, then investigate whether finite-order elimination or a quantitative locality restriction can close the perturbative route. The whole-class verdict remains nonconclusive.

Substantive step 8 completed.

`exact_encoding.md` now proves the exact product-code theorem, its Clifford dimension corollary, and a quantitative leakage bound ||E_j||>=|J| for the parity-dressed hopping target. It also derives, from an explicitly stated Lieb-Robinson estimate, the endpoint commutator bound C exp(-mu r/2)+4 exp(-Delta r/(4v)). This forces a gap/length/amplitude tradeoff for local parents, but the C7 shrinking gap is compatible with it.

Substantive step 9, flat-parent second order derived and checked.

Choose a joint code using only physical operators odd under a port grading, so all compressed same-port second moments are delta_ab I and all left/right second moments vanish. For the flat parent H0=Delta(1-P), the full second-order chain term is a scalar plus (lambda^2 alpha^4/Delta)[(3/2)sum_j T_j^2 + sum_j {T_j,T_(j+1)}], where T is the logical Majorana hopping bond. The adjacent-bond term is essential. `c8_resume_checks.py` verifies the physical transition amplitudes on three open blocks at (q,m)=(2,2),(3,4). The largest second-order error is 2.8e-14.

For the symmetric free target, the staggered mass-channel coefficient of this expression is 3/2-2=-1/2. Its interpretation as an attractive, marginally relevant SO(m) Gross-Neveu coupling for m>2 uses continuum perturbation theory, not a proved microscopic gap theorem. This makes the flat parent an unsafe completion of the first-order target.

A new way past this particular obstruction is being checked. Two excited energies can keep the same single-port spectral weights while changing the cross-port weights. The preliminary formula gives a positive mass-channel coefficient for an explicit positive spectral measure. This would remove the forced sign at second order in unrestricted blocks, without solving microscopic locality or all-order stability.

Substantive step 10 completed, explicit spectral repair.

`second_order_stability.md` proves the full adjacent-bond calculation and gives a two-energy parent with exactly the same ground code and first-order hopping, but a positive second-order singlet mass-channel coefficient. For E=(1,3), p=(1/2,1/2), s=(-alpha^2/2,-alpha^2/2), and cross-port weights t=(0,-alpha^2), the coefficients are A=(25/24)alpha^4 on sum T_j^2 and B=(1/3)alpha^4 on adjacent anticommutators, hence A-2B=(3/8)alpha^4. Positivity follows from explicit weight eigenvalues. A parity-preserving dilation gives an actual finite matrix H0.

The q=3,m=4 example is saved as `spectral_parent_q3_m4.npz` (dimension 324, gap 1). Physical kernel, parity, spectral-weight, same-bond, and adjacent-bond checks are all below 1.0e-14. The first parity implementation failed because its Choi factor and assumed environment grading used different bases; switching to the specified positive square root fixed the mismatch. `c8_resume_checks.py` now runs both flat and engineered parents and passes. No finite-lambda CFT was measured. Higher-order marginal anisotropies and microscopic repeated-bond locality remain open.

Substantive step 11 completed, all-ring restriction.

`ring_obstruction.md` proves that a translation-invariant chain with conserved binary markers, a positive energy-density separation from constant marker sectors, and a bounded incompatible-ring fitting cost has at least L ground states on every sufficiently large prime ring of length L. This excludes convergence to a fixed finite collection of compact CFT sectors with finite vacuum multiplicity on all rings. It applies to the frozen C7 selector and any spectator extension retaining the stated energy separation. A separate finite-band proposition proves O(L^-2) defect levels after ordinary finite-range defect hopping; coupling to a gapless bath is outside that proposition.

Substantive step 12 completed, small local-parent search.

`local_parent_scan.py` tests all six placements of the three environment qubits between the boundary ports of the q=2,m=2 five-qubit positive-square-root Choi code. No nonconstant uniform two-site bond has this code as a degenerate eigenspace. Smallest singular values of the exact linear constraint matrices are about 1.408, far from numerical tolerance. Even allowing nonuniform nearest-neighbor couplings gives no solution for five orderings. The sixth has two null directions with code energy zero and traceless nonzero Hamiltonians, so none of their nonzero combinations can make the code a ground space. These are numerical exclusions of these six specific dilation embeddings, not a no-go for all equivalent dilations or local codes.

Substantive step 13 completed, a way around the isolated-defect setup.

`mobile_selector.md` gives an explicit additional nearest-neighbor term -u(|#L><R#|+h.c.) in C7's d=4 selector. It preserves every legality projector and total marker number but allows marker motion. At the crossing between preferred lengths b and b+1, its first-order projection exchanges neighboring short and long cells with amplitude -u A_b^2, A_b=sum_i phi_b(i)phi_(b+1)(i). This is an XX chain of cell lengths, with reference-marker and twist bookkeeping on the physical ring. Every sufficiently large ring admits the two-length manifold. The leading thermodynamic energy favors an extensive mixture, not one forced defect.

`mobile_selector.py` checks local legality exactly by enumeration and diagonalizes the physical marker-pattern graphs. For b=5, half-long-cell sectors at L=33,44,55,66 have connected graphs and L*gap=2.3854,2.3406,2.3203,2.3094. Their ground energies match the projected free-XX formula to 3e-14. For one long cell at L=31,41,51,61, L^2*gap instead approaches 4pi^2. This distinguishes the finite-density fluid from an isolated quadratic defect. It is not a finite-u microscopic CFT determination.

New next step: combine a two-length local payload code with a marker hop that transports its logical state coherently. The code must have compatible boundary maps at both lengths, and the coupled length/payload theory must remain conformal. No such payload code is constructed yet. The earlier assumption that frozen blocks are sufficient even after addressing the first-order code is no longer tenable under the all-ring requirement.

Final audit of the resumed run completed. RESULT.md contains the complete new proofs, explicit spectral data, numerical reproduction instructions, and the updated ranked next steps. Compilation passed for all three new scripts. No conclusive whole-class result was reached. The main directions that stalled were exact product-code embedding (excluded for the stated hopping target), the flat-parent completion (a dangerous singlet sign, repaired only at second order), and six direct local-parent fits (numerically excluded for those embeddings). The mobile two-length selector is the next construction direction, with payload transport and infrared control explicitly open.

## Third run, continuation from the saved state

Read the prior progress before project files. The latest programmable-cell contribution has a proved finite-history port-idle interface, but its own supporting report still explicitly lacks a translation-invariant embedding theorem. Its stability and ring caveats remain open. No existing C8 checks were rerun for orientation.

New directions being developed, without repeating the earlier algebraic dilation work:

1. The mobile selector at equal head and marker hopping appears exactly equivalent, in its legal sector, to even-number hard-core domain walls with two boundary twists. This should give a finite-coupling c=1 selector on all ring lengths, with a direct bound excluding illegal sectors. It does not carry the desired payload yet.
2. A vacuum-padding construction appears to give arbitrary existing two-port codes compatible with TWO cell lengths: equal endpoint compressions, and a projected boundary shift proportional to the identity on both logical payloads. This is an abstract code compatibility construction, not a repeated-bond parent.
3. If these work, examine whether vacancy doping of a hypothetical compatible-ring counterexample avoids frozen defects while preserving its violation. Separate any exact algebra from the spin-charge separation assumption.

Next step: prove and check 1 and 2, then examine the effective coupling between length and payload modes. Whole-class status remains nonconclusive.

Substantive step 14 completed, exact finite-hopping selector.

`equal_hopping_selector.md` proves that the d=4 mobile selector at t=u=a is exactly even-number hard-core wall hopping with two boundary twists. The full legal spectrum is free, epsilon(k)=-2a cos(k)+mu/2. Every illegal pattern with k violations has energy at least E0(legal)+k(K-2a-|mu|); K>2a+|mu| therefore gives a uniform gap to all forbidden sectors. For |mu|<4a the continuum is the standard c=1 free-fermion CFT with a finite parity projection/twist completion, on all ring lengths. There is no long-period logical payload in this theorem.

`equal_hopping_selector.py` passed full-spectrum comparisons at L=5..9 and all nonempty violation patterns at L=4,5,6. At mu=0, L divisible by four, E0=-2a csc(pi/L), v=2a, and the first opposite-twist gap is 2a tan(pi/(2L)). The physical-site interval entropy is EXACTLY the wall-chain entropy plus log(2). Direct entropy fits at L=20 give c=1.003981, with half-chain size slope c=0.992700 from L=16 to 20. The independent Casimir fit at L=40..160 gives c=0.999999704. Outputs are in `equal_hopping_selector_checks.json`.

Substantive step 15 completed, two-length payload algebra.

`padding_transport.md` gives W_p=(V tensor vacuum^p + vacuum^p tensor V)/sqrt(2), p=1,2, for ANY k-site two-port code V. Adding one local vacuum symbol costs q -> q+1. Both lengths have the same endpoint compressions, one half those of V. Moving their shared boundary gives (W_2 tensor W_1)^dagger(W_1 tensor W_2)=I/4. Hence the projected marker hop preserves both payloads and the projected payload bond is independent of lengths. This is NOT an invariant encoding, and no repeated-bond local parent is provided.

`padding_transport.py` checked the actual q=2,m=2 Choi code with k=5, lengths 6 and 7. Isometry, both ports, coherent transport, and all four length-pair bond identities pass below 5e-16. Formal graded resource accounting q=17,m=144 gives padded d=72 and payload c=72, still above d-1, but this is only a target. Data are in `padding_transport_checks.json`.

Next: investigate whether the padded code can have a local parent, and whether mobile vacancies offer a distinct all-ring completion. Do not conflate the exact equal-hopping selector with the small-hopping long-cell regime.

Substantive step 16 completed, locality audit and a conditional local repair.

The exact cat padding of step 15 cannot be the sole ground code of a nearest-neighbor parent. The two placements differ in occupation at both endpoints, so every local cross matrix element vanishes. If their symmetric combination minimizes H, positivity forces each placement separately to minimize H, leaving at least 2D ground states. This theorem is in `vacancy_parent.md`, Section 1. It is why step 15 alone does not solve local realizability.

A constructive repair is proved for ANY frustration-free core with uniform positive bond h0. Add one vacancy symbol and use the same bond U(n_i-n_(i+1))^2 - w*(vacancy-label swap) + g*h0 on both lengths. In the sector with k occupied labels, vacancy hopping preserves their order. The unique positive charge ground state chi_p tensored with the original core code is the COMPLETE ground code, for p=1,2 holes. Endpoint maps are r_p times the original ones, where r_p tends to 1/2 as w/U tends to zero. A cut shift still compresses to beta*I on both payloads, with beta tending to 1/4. Unlike the cat construction, r_1 and r_2 need not match at finite w. The small odd/even tunneling gap must be respected when coupling cells.

`vacancy_parent.py` checks a genuine three-qubit symmetric ferromagnetic core (D=4, noncommuting compressed ports), with repeated physical bonds on lengths 4 and 5, at w=0.15,0.3,0.5. All code, port, and two-cell transport errors are below 1e-15, and ground splittings below 3e-15. The smallest tested charge gap is 5.8573e-6. The initial w=0.15 run failed the 1e-10 port tolerance because the unconstrained eigensolver mixed the very close reflection partners at roundoff scale (error 1.244e-10); projecting onto the mathematically known even-reflection sector resolves this and retains an eigenvector residual below 1e-12. Outputs are `vacancy_parent_checks.json`.

Whole-class status remains nonconclusive. Missing: a high-channel frustration-free uniform core, selecting the required occupied-number sectors in a full ring, and controlling length/payload interactions at finite coupling. Vacancy doping of an entire chain would preserve the formal excess c-(d-1) if spin-charge separation held, but its gapless charge denominators prevent a volume-uniform conclusion from elementary perturbation theory. This line is recorded as conditional, not proved.

Substantive step 17 completed, finite witness for a restricted selection route.

`finite_ff_witness.md` proves that for each q there is a finite K_FF(q) such that a uniform positive nearest-neighbor bond having a zero-energy OPEN-chain state at K_FF(q) has one at every length. Proof: the determinant polynomials det H_N generate a finitely generated ideal in q^4 real variables (Hilbert basis theorem); positivity makes existence at length K imply existence at all shorter lengths. No effective value of K_FF is obtained.

Thus selecting an arbitrarily long core as the LONGEST frustration-free open segment is impossible at fixed q when there are no extra boundary terms. It gives c_target<K_FF(q) log_2(q) for the particular first-order Clifford target with that selection rule. It does not apply to boundary computations or C7's negative kinetic-energy selection and does not establish a whole-class c bound. The same note proves why a vacancy chemical potential alone cannot select the partially occupied examples when U>=w and fully occupied cores are frustration free.

The alternative routes have now been followed beyond the first partial results. The exact selector works, the exact cat parent fails, a finite-hopping local replacement works under a frustration-free-core hypothesis, and its number-sector selection and general high-channel core remain unresolved. Next is to integrate the new proofs and check the final report without promoting these restricted results to a conclusion.

Final audit of run 3 completed. RESULT.md Sections 15-19 contain every new proof, construction, check and limitation, with a fresh ranked handoff. The largest complete-spectrum error is 1.96e-14 and the largest entropy-identity error is 4.41e-13. Python compilation passed for all three new scripts. The restricted frustration-free witness argument was checked for its positivity, nesting and boundary-term hypotheses. No conclusive outcome has been reached. Only this folder was modified; no commit, push, network, cluster computation or independent referee check was performed.
