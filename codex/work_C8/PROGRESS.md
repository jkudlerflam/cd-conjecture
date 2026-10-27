# C8 progress

2026-10-27. No conclusive answer. The completed report is `RESULT.md`.

## Current best answer and restart point

The resumed run is complete and nonconclusive. `RESULT.md` now contains all original results plus the full proofs and constructions from `exact_encoding.md`, `second_order_stability.md`, `ring_obstruction.md`, and `mobile_selector.md` in Sections 10-14.

New proved results: exact invariant product codes force commuting active boundary algebras, excluding the old single-set hopping target as an exact separated-endpoint encoding. There is a leakage lower bound and a Lieb-Robinson gap/length/amplitude tradeoff. The full second-order flat-parent Hamiltonian includes adjacent bonds; an explicit two-energy, parity-preserving parent changes its singlet mass-channel coefficient from -1/2 to +3/8. Frozen-marker chains satisfying the energy-density hypotheses have at least L ground states on large prime rings. Ordinary finite-band defect motion has O(L^-2) levels.

A new construction direction is explicit. The local move #L <-> R# added at the crossing of cell lengths b,b+1 projects to XX hopping of the lengths with amplitude u A_b^2. A finite-density mixture has projected gaps consistent with 1/L, while one long-cell defect has 1/L^2 gaps. This does not yet include a logical payload or establish the finite-u microscopic CFT.

Open problems, ranked:

1. Build a repeated-bond local code for BOTH lengths b and b+1, with compatible port maps and a marker hop that transports the logical state coherently. Do not assume one-length frozen blocks handle all rings.
2. Establish thermodynamic effective-Hamiltonian control and infrared stability of the coupled length/payload theory. The spectral repair fixes one second-order sign, not all marginal couplings.
3. Independently, settle the two conditional premises for the rank-one d=20992 chain in RESULT.md Section 6. It does not require small-d blocks or a selector.
4. For a bound, find a restriction stronger than the new gap/length tradeoff that cannot be evaded by arbitrarily small gaps and exchange scales.

Latest reproducible checks:

- `python3 c8_resume_checks.py`, with single-thread BLAS: new finite-dimensional identities pass, maximum flat-parent second-order residual 4.9e-14; the explicit 324-dimensional spectral parent has gap 1 and kernel/parity/weight/effective-coefficient errors below 1e-14. Writes `resume_checks.json`, `spectral_parent_checks.json`, `spectral_parent_q3_m4.npz`.
- `python3 mobile_selector.py`: all local legality checks pass; b=5 graphs at L=33,44,55,66 match the free gap-chain energies to 3e-14 and have L*gap approaching about 2.3; one-defect L^2*gap approaches 4pi^2. Writes `mobile_selector_checks.json`.
- `python3 local_parent_scan.py`: all six tested five-qubit dilation embeddings lack a uniform nonconstant two-site parent; this is numerical and specific to those embeddings. Writes `local_parent_scan.json`.
- Compilation checks passed for the new scripts. Earlier `c8_checks.py --ed` results remain in the run 1 record and were not repeated.

No full microscopic CFT numerical calculation, network search, cluster computation, or independent referee check was done. All project changes stayed in this folder. An incidental temporary stdout capture in /private/tmp was created and removed.

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
