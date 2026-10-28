# C8 progress, ninth run in progress

2026-10-28. Resumed from the completed eighth-run controller. Read the
project memory, repository status/background, requested historical sources,
C19 and C21-C23 handoffs. No external writes, network or subagents.

The earlier main target is complete in the saved run. This continuation
is auditing the full microscopic composition and the remaining infrared
step. The source audit has found no controller error so far. The C21/C22
local source constructions still have the reported target-band obstruction
and are not being duplicated.

Active work: (i) implement the full periodic bond as a sparse column oracle,
including wrapper, selector, helices, normalized ports and density fields;
(ii) use exact total U(1) to replace exponentially crowded geometric Fermi
points by an explicit rational Sidon construction with finite velocity
bounds; (iii) examine which massive-sector estimates actually follow from
the local defect bounds. Neither a conclusive result nor an all-order
elimination theorem is currently claimed.

Substantive step 41. periodic_clean_bond.py now implements the complete
periodic matrix, including all illegal sectors, all wrapper/filter/helix
terms, marker hopping, the two physical ports and colored fields. Its
q=4 dimension is 214272048, but individual columns remain sparse. Checks
pass on 3232 base and perturbation columns separately, with exact charge,
marker number and selector-legality preservation. A separate dense wrapper
implementation agrees on all 225 columns of a two-state payload fixture.
Actual q=4 history readouts reproduce the 16-dimensional two-cell fermion
Hamiltonian with zero matrix discrepancy. No microscopic CFT measurement.
The first assembly left charge-violating roundoff entries from cancellation
of Hermitian components; direct charge-definite assembly removes them.

Substantive step 42. rational_target.py supplies exact integer Sidon
momenta for the total-U(1) core. All pair sums are unique, so the finite
quadratic/quartic detuning margin is at least pi/(8 A), A=O(m^3).
Payload velocities exceed 2 sin(pi/64), independent of m. Mixed low-shell
margins exceed 7pi/8 at length density 1/2. Integer checks pass through
m=256 (512 branches and 131328 unordered pair sums). This replaces the
exponentially slow geometric target, conditional on the same IR inputs.

Active next steps: prove polynomial bounds on all visibility and wrapper
normalizations; then isolate what is missing in the connected massive
excursion argument. The original task remains nonconclusive.

The previous completed handoff follows and remains the source for all
unchanged core proofs and checks.

# C8 progress

2026-10-28 UTC. Eighth run completed. Nonconclusive for the original
c-d conjecture. The autonomous high-channel initialized core is now
explicit, and its alphabet is counted. Start with this summary.

## Current best answer

The first-order core bottleneck from run 7 is closed by
clean_readout_core.md and RESULT Section 40. For q divisible by four,
m=q^2/8 and n=m+6, the 15-letter-clock core has site dimension
446400(q+1), complete zero-space dimension 2^m, a positive explicit
finite-length gap, both Clifford endpoint maps, colored flavor densities,
and exact on-site U(1) restricting to total logical number. All source
records of the exposed bit are removed before READY, at every partial
clock step. The input pins reject every surplus work/head state.

The machine has 58 direction-fixed controls, 799 prescribed collision
blocks and 928 completed blocks. Its counters live on four existing
q-state tape registers. The cursor selects the logical bit. No generic
compiler, external instruction list or program verification is needed.
The full cycle has 11m^2+105m steps, below the clock's 2^(m+6) sweeps.
The first-order logical target is C19's displaced flavor hopping plus
independent fields; C19's second-order spectral filters are NOT realized.

Composing with the already proved mobile wrapper/equalizer gives
physical d=96*446400(q+1)+48 and projected m+1 modes. A counted resource
example is q=685670400, d=29383993632614448, projected c=58767987179520001.
This is conditional resource accounting, not a measured violation.
run8_microscopic_reduction.md and RESULT Section 41 specify the bond
family, exact symmetries, an explicit nonresonant target, and two remaining
inputs: uniform massive-sector elimination in the mobile critical system,
and Gaussian infrared stability/identification. No conclusive original
answer, full-class bound or fixed-d unbounded-c result is claimed.

## Files to read and run

* clean_readout_core.md, the full main theorem, construction, initialization,
  grading, all-time proof, counts and reproduction instructions.
* clean_readout_machine.py and clean_readout_table.json, actual complete
  controller. q-dependent operations are six fixed record-gate families
  and seven onsite families, each explicitly defined on its whole space.
* clean_core_bond.py, sparse column oracle for the complete microscopic
  positive core bulk/boundary matrices, including every illegal sector.
* clean_readout_checks.py, cycle and full-clock checks.
* clean_readout_physical_checks.py, direct physical port matrices, the
  entire local initialization alphabet, and full-space involution tests.
* clean_core_kernel_checks.py, actual Hamiltonian acting on sparse codes.
* clean_readout_counts.py, exact active counts from one polynomial cycle
  and binary valuation sums, without exponential clock simulation.
* run8_microscopic_reduction.md and clean_readout_resources.json, counted
  periodic candidate and the remaining finite-coupling assumptions.
* endpoint_shuttle.py, endpoint_shuttle.md and associated checks/table,
  an earlier independent 21-control generic transport module. Useful for
  the sparse-route compiler, but not needed inside the new direct machine.

## Verification completed

All eight new Python files compile. All new scripts passed.
Full primitive permutation and exact number checks at q=4,8,12.
All 4/256 inputs at q=4/8 pass complete cycles; 24 q=12 samples pass.
At q=8, 51200 paired active environments agree under flipping the exposed
bit. Full-clock checks cover 5124 vertices for every q=4 input and
524298 vertices at q=8. Direct physical compression at q=4 has zero
integer error. The full local Hamiltonian annihilates all four q=4
histories with zero residual, using 21002 matrix entries per column.
The local input pins reject 148790 onsite q=4 columns, retaining exactly
the stated 1 or 2 by role. Full-space G,F tests cover 12000 samples.
Hermiticity and charge checks test 1600 bulk columns at each q=4,8 and
both boundaries. Polynomial count results agree with both explicit-clock
runs and extend to the 872415252-vertex q=12 history.

No microscopic CFT numerics, cluster jobs, network, subagents, independent
referee, commit or push. Only this folder was modified. C21/C22/C23 PROGRESS files were absent initially, appeared during the run,
and were read in the final audit. They report local second-order buses
with target-band obstructions and a conditional assembly ledger, not a
completed C19 compiler. These are source reports, not audited C8 results.
New sparse-route and C19 material was read but not modified or counted
as new C8 work.

## Remaining work, ranked

1. Independent adversarial audit of the new machine and the full-space
   core bond. Focus on reverse collisions, full initialization, exact
   charge, color support and partial-time flags. Do not repeat the old
   clock, parity toy or unspecified compiler argument.
2. Establish the massive-sector integration input for this explicit
   first-order mobile model. The complete finite code, positive gap and
   additive defect-energy/volume bounds are now available. Finite-order
   SW and prethermal control alone do not establish a CFT at fixed coupling.
3. Apply the number-conserving C13 Gaussian analysis with the explicit
   momenta of Section 41, tracking the smallest velocity, complete physical
   norms, mixed resonances and ring zero modes. A microscopic two-estimator
   study needs a much smaller equivalent implementation or cluster work.
4. Only separately revisit C19's fixed-d spectral filters. The first-order
   core does not fix their source resolvents or prove unbounded c at fixed d.

The next run should begin with an independent audit or the concrete
massive-sector reduction, not try to reconstruct the completed controller.

Final audit: the positive core bulk/boundary norms obey 17,7,6 at unit
scales, independent of q and n. Its sparse column oracle also runs at
q=685670400 without allocating a q-sized or block-sized matrix. All new
prose was checked for em dashes and forbidden confidence language. The
complete current handoff is duplicated as RUN8_RESULT.md and
RUN8_PROGRESS.md to survive the shared-repository synchronization issue
seen at the end of run 7. Only local files were changed.

## Earlier chronological record


## Run 8, active continuation

Substantive step 37, new main result, 2026-10-28 UTC. The autonomous
clean-readout machine is now implemented in clean_readout_machine.py.
It avoids a program generator: a moving logical cursor selects the bit,
and two q-state registers next to each port hold row/column counters.
All counters are on the tape. A constant set of explicit two-site modular
permutations updates them. The head carries only a neutral parity bit
and the single graded data bit being moved; that data bit is zeroed in
the head before READY. Both port maps and colored density maps are clean.

The actual table has Q=58 direction-fixed controls, 799 prescribed
collision blocks and 928 completed blocks, each decorated by one of seven
explicit local permutations. The work dimension is 29760(q+1), and the
15-letter-clock core dimension is 446400(q+1). For q divisible by four,
m=q^2/8 logical data sites plus six endpoint/counter sites give n=m+6.
Static role syntax selects this layout, and role-dependent I-sweep pins
initialize every surplus bit, record, flag and head state. There is no
preloaded program and no position-dependent instruction list.

Checks passed: full primitive permutation/grade tests at q=4,8,12;
complete cycles for all 4 inputs at q=4 and all 256 at q=8, plus 24
q=12 inputs. Every paired active environment was identical under flipping
the exposed bit (51200 pairs at q=8). The full machine resets in
T=11m^2+105m sweeps, below 2^(m+6) for m>=2. A first guessed timing
formula failed at m=8; it was replaced by this derived formula.
Every one of 5124 clock vertices was tested for all q=4 inputs, and every
one of 524298 vertices at q=8 for a representative input, comparing with
the machine after each complete sweep. Nonuniform flavor counts occur
and must be compensated. No central charge measurement is claimed.

Next, urgently: write the all-sector initialization/core theorem with
proof, independently check the full initialized kernel in a small
reachable-state fixture and physical compression matrices rather than
only the clean-state invariants, and audit compatibility with the mobile
wrapper. C19 spectral filters are not implied by this first-order core.
The original conjecture is not yet declared settled.

2026-10-28 UTC. The saved run-7 results are retained below. No conclusive
answer has been obtained in this continuation.

Read STATUS, route-A background, C19 and the new sparse-controller files.
C21, C22 and C23 PROGRESS files are absent in this checkout. The new
controller_instantiation contribution supplies a preloaded dispatcher and
a blank length counter, but explicitly leaves their autonomous composition,
format verifier and endpoint readout incomplete. Its status is narrower
than the earlier existential-core composition claim. No external files
were changed and no old clock calculation was repeated.

Substantive step 36. Wrote endpoint_shuttle.py and its complete 378-entry
collision permutation, with 336 prescribed entries, 21 direction-fixed
controls and 18 neutral tape symbols. It transports four marked output
records to the two physical endpoints and their inner neighbors, opens
READY windows, clears them, and reverses all record moves. The cycle uses
no stored travel history and no q-valued controller. It is conditional on
prepared, marked outputs, not a complete initialized core. Its standalone
15-letter-clock factor is 11610 times the quantum-record dimension.

The first check caught an incorrect proposed cycle-length formula, which
omitted the source positions. Corrected it to 10n+8a-4b-12. Full local
unitarity, complete cycles, all partial-clock vertices, both Clifford maps
and colored density checks are now running in endpoint_shuttle_checks.py.
All checks passed. endpoint_shuttle.md contains the complete proof and
count. Checks covered 75 cycles through length 32, 12096 bond cases, every
clock vertex at lengths 8,9,10, and both graded Clifford/color maps through
(q,m)=(8,8), logical dimension 256. Port and density active counts differ
and must be normalized separately. Next: surplus-input checking and the
C19 spectral interface; the full autonomous encoder is still missing.

2026-10-27. Seventh run finished, nonconclusive. RESULT.md Sections 36-39
contain the new proofs, constructions, checks, and ranked continuation.
All writes stayed in this folder. Read this first and do not redo
completed work.

## Current best answer

No verified microscopic counterexample, full-class bound, or fixed-d
unbounded-c construction has been obtained. The initialized clock is now
explicit, but its high-channel clean-readout work program is still missing.
Finite-coupling infrared control is a separate task, treated by C13.

A nine-letter uniform positive NN clock has exactly one legal history
path at every length n, with C_n=n+4(2^n-1) vertices. The path includes
an initialization sweep. A fifteen-letter extension inserts complete
outward/return sweeps and has C_n=(2n+4)2^n+n-4 vertices, with exactly
2^n full sweeps. Every local legality-defect projector is conserved.
There are no extra penalty-free histories, cycles, or internal lengths.

The full-sweep version executes a partitioned local work update. Any
specified finite injective graded collision table can be completed
canonically to a permutation and executed on this clock. Initialization
and later checks determine the COMPLETE quantum ground code, with a
finite-length gap. The work machine may have arbitrary other orbits;
these cannot create additional ground states after initialization pins.
This removes the need to assume a chosen work-machine sector.

It does not supply the desired collision table. The tested packet
machine computes total logical parity at a boundary, but its off-diagonal
accumulator compression is zero because the logical bits remain recorded
elsewhere. A successful readout must erase all other records of the
exposed bit before activating the port. All active partial-sweep times
must have the intended matrix elements and input-independent counts.

## New proved results and files

1. initialized_binary_clock.md gives the nine-letter alphabet, repeated
   positive bond and boundary terms, complete path proof, initialized
   code, gap, hopping spectrum, and a merged frozen selector. Strict
   convexity selects every interior length n>=2. Full illegal sectors
   are controlled. It is not a conformal ring model.
2. sweeping_clock.md gives the fifteen-letter full-sweep extension,
   a local update identity on the whole work space, canonical collision
   completion, and input/check-kernel theorems. For classical permutation
   programs with diagonal checks, every rejected input has a path pin,
   giving gap >=min(K_penalty,J[1-cos(pi/C)],eta/(2C),J/(4C^2)).
3. readout_machine_contract.md specifies the unfinished clean-readout
   program. If its ordinary tape alphabet and packet state set are fixed
   as q grows, the core dimension is K(q+1), with K=30 s0 r^2. Both
   boundary Clifford maps then follow exactly by counting active times.
   The conditional wrapped dimension is d=96K(q+1)+48 and target c=q^2/8.
   K is NOT supplied. This is a conditional reduction, not a counterexample.
4. run7_audit.md is the final audit and ranked continuation, also included
   as RESULT Section 39.

## Verification this run

All four new scripts passed with single-thread BLAS and compiled.

- binary_clock_checks.py exhausts all clock words through n=5 (59049
  words), checks every defect on every edge, and computes complete small
  illegal-component spectra. Both legal spectra through n=8 (1028
  vertices) agree below 4e-15. A 224-dimensional entangling initialized
  core has exactly four zero states, gap .0062877901, and code residual
  below 1.7e-15; the proved gap bound is .0012755102.
- sweeping_clock_checks.py exhausts graphs through n=4 (50625 words,
  16477 components, 3855 illegal components containing cycles). Full
  legal spectra through 1026 vertices agree below 5.3e-15. On a
  64-dimensional work space all eight completed updates agree with
  powers of the intended local update below 8e-15.
- sweep_machine_checks.py builds a literal eighteen-state collision
  permutation, checks grading, and tests every initialized input through
  n=5. The diagonal parity readout succeeds with zero residual and the
  off-diagonal readout vanishes. This is a controller test, not a
  high-channel code or a CFT measurement.
- clock_ring_checks.py computes every component spectrum on complete
  ten-letter rings of dimensions 1000,10000,100000 at L=3,4,5. The
  selected ground multiplicities are 3,4,5, and the analytic energies
  agree within 3e-16. Every illegal component obeys the ring bound.

Outputs have the corresponding .json names. Earlier suites were not
rerun. No microscopic CFT estimators, network, cluster jobs, subagents,
independent referee, commit, or push were used by this run.

## Earlier interfaces, retained

For a q-state initialized positive FF core of length k>=3, the existing
selector/wrapper has dimension 24q+12 and selects exactly b/(b+1) cells,
b=k+1, on every sufficiently large ring, with the complete unperturbed
manifold proved. Marker hopping transports payloads by a scalar. Two
qubit helices equalize leading couplings at d=96q+48. The full projected
length/payload spectrum and the physical K=(b+nu)^(-2) are known. Positive
approximate-core bounds are also available. These are RESULT Sections
30-35 and their separate notes. They have not been re-derived this run.

C13 finished while this run was active. Its final PROGRESS reports a
conditional positive Gaussian-phase analysis and an all-integer fixed-
vacuum qualification. Generic projected numerics remain too small for
accurate central charges. It still requires the actual microscopic core
and corrections. These are source-reported results, not independently
audited C8 claims.

## Remaining problems, ranked

1. Write the actual clean-readout collision table described in
   readout_machine_contract.md. Begin with the fifteen-letter clock.
   Keep control states independent of q, put counters on the tape,
   preserve parity, test surplus input bits, include colored flavor-density
   readouts, and verify both ports at
   ALL active partial-sweep times. Count s0,r and hence K. A generic
   reversible simulator retaining copies of the exposed bit fails.
2. Put that core through the existing wrapper/equalizer and calculate
   actual finite-coupling corrections below its very small but positive
   cell gap. A cell gap does not isolate an extensive many-cell band.
3. Use C13's resonance and sector analysis for the microscopic model,
   with calibrated estimators or a complete reduction to stated IR
   inputs. Do not duplicate its projected-fluid work.
4. Only then consider a universal quantum work program for route A.
   The clock does not supply the required excited-state spectral weights
   and does not prove unbounded c at fixed d.

Do not redo the nine-letter counter, full-sweep clock, defect conservation,
initialization-gap proof, bare-clock selector, or packet parity toy. Do
not redo the earlier joint dilation, spectral synthesis, exact-code
no-go, frozen-ring no-go, graded-swap core, FF witness, boundary-rank
bound, classical readout code, or mobile-cell interface.

Next action: a concrete finite reversible clean-readout program and its
collision table, not another clock or an unspecified universality claim.

The remainder is the chronological record. Its initial task list has
been superseded by the ranking above.



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

## Fourth run, continuation from the third-run handoff

Read the saved progress and requested context. No newer contribution
supplies the missing repeated-bond embedding or finite-coupling IR theorem.
All work remains confined to this folder, with no network or subagents.

Substantive step 18 completed analytically.

`boundary_rank_obstruction.md` proves a new necessary condition for ANY
scaled Clifford code. If even one code vector has Schmidt rank s across
a q-dimensional active port, the number of channels satisfies
m<=2qs-s^2. In particular, one product vector forces even m<=2q-2 and
c_target<=q-1, closing the first-order violation window for such codes.
The proof counts real orthogonal vectors Gamma_a psi, not operator
Schmidt rank. Small nonzero compression amplitudes do not avoid it.

For a positive uniform bond of matrix rank r<q, every chosen endpoint
vector extends to a product ground state at every open length. Exact
involution compressions then preserve the entire ground code and obey
the physical Clifford dimension bound. Joint parity-dressed ports are
impossible. Quantitative Gram and leakage inequalities show that small
relative Clifford errors do not fix this subclass. This does not apply
to a selected subcode or to extra boundary constraints.

Next: check these statements numerically, classify the q=2 finite FF
witness explicitly, and test a higher-rank uniform quantum core to locate
how it evades the stronger low-rank obstruction. No whole-class conclusion.

Substantive step 19 completed.

`qubit_ff_witness.md` proves the sharp K_FF(2)=4. The proof classifies the
product-vector intersection of a two-dimensional two-qubit allowed
space into a secant, tangent, or all-product projective line. The explicit
positive bond h=I-|01><01|-|1+><1+| has a one-dimensional three-site
kernel and no four-site kernel. Its q-state extension proves
K_FF(q)>=q+2 for all q>=2. No upper bound at q>=3 is claimed.

`boundary_ff_checks.py` passes. Product-extension and compressed spectral
extrema checks at q=2,3,4 and N=4 have errors below 3e-15. All real-image
dimensions match 2qs-s^2. Eigh checks cover six qubit line types and common
invertible changes of physical basis at N=3,4,5. Sparse checks verify the
sharp-witness examples at q=2,3,4 up to physical dimension 4096. These
checks validate finite-dimensional formulas, not a CFT claim.

Next: test a natural higher-rank graded-ferromagnetic core and turn to a
symmetry-based rank bound if its channels only saturate the budget. The
previous universal conclusions remain open.

Substantive step 20 completed, constructive higher-rank test.

`graded_core.md` gives a genuine repeated positive NN bond
h0=(I-P_graded)/2 on q=n+1 states. Its complete N>=n ground code is the
2^n-dimensional signed subset space. Both endpoint maps are derived
exactly, including parity dressing. The normalized Clifford error is
1-sqrt(1-(n-1)/N). Thus long cores approximate 2n channels, but only
c_target=n=q-1, before selector cost. Adding even colors cannot bypass
the product-vector obstruction on the complete code. This natural
higher-rank construction therefore stalls as a violation route.

`graded_core_checks.py` checks five (n,N) pairs up to physical dimension
1024. Both ports and the error formula agree below 3e-16, and code
isometry below 1.5e-15. Ground multiplicities are checked by full spectra
through dimension 256. These are code tests, not CFT measurements.

Substantive step 21 completed, a different conditional bound.

`symmetry_rank_bound.md` proves that the part of a microscopic on-site
Lie algebra acting on the connected bond embeds into direct_sum so(n_j),
where n_j are the multiplicities of the bond's Schmidt singular values.
Its kernel consists of independently conserved single-site generators.
If microscopic dynamical currents exhaust the IR stress tensor, this
gives c<=min(d^2-1,sum_j n_j(n_j-1)/2). For Gaussian modes resolved by
commuting on-site charges it improves to c<=sum_j floor(n_j/2)<=floor(r/2).
The proof is independent of perturbative order. It leaves neutral and
emergent current sectors open, so is not a general proof of c<=r(h).

`schmidt_symmetry_checks.py` verifies the induced algebras and kernel
for Ising, XX, XXZ, XYZ, Heisenberg, and two independent XX chains.
The largest invariance residual is below 4e-16.

Also added two short parent exclusions to `boundary_rank_obstruction.md`.
Commuting FF projectors give commuting compressed disjoint-endpoint
algebras. A standard complete injective MPS edge code with a nonzero
2-site parent has D_virtual<=q-1 and its even Clifford target obeys
c_target<=q-1 by dimension and divisibility. High-rank noncommuting,
non-MPS cores with entangled boundaries remain possible.

Next: audit the new finite-dimensional proofs and their hypotheses,
incorporate them into RESULT.md, and rank the remaining construction
and neutral-sector problems. No whole-class conclusion has been reached.

Final audit of run 4 completed. RESULT.md Sections 20-24 contain every
new proof, construction, numerical specification, and ranked next step.
All three new numerical scripts ran successfully and compiled. The
later edits only clarified proofs and synchronized the reports, so no
unchanged numerical suite was rerun. The commuting-parent and MPS
corollaries were checked algebraically; no independent referee was used.
The earlier Choi code was checked analytically against the new Schmidt
obstruction and has full boundary rank for all its code states.
The verdict is still nonconclusive. Only this folder was modified.

## Fifth run, continuation from the fourth-run handoff

Read the saved state and current requested context. The new ir_control
folder contains separate tasks on effective-theory control and operator
classification, so this run will not duplicate them. The programmable-cell
outline still lacks the complete repeated-bond embedding and all-ring
completion. The spin-3/2 tricritical search has returned a negative result
at the sizes tested. No network, cluster, or subagent work is being used.

New attack in progress: an explicit nonnegative, equal-amplitude joint
Clifford code. Instead of quantum state preparation, use a reversible
classical readout that temporarily moves one logical bit to a boundary
pair of physical letters. Two signed physical pairs per logical flavor
appear sufficient, giving n<=floor(q(q-1)/4) complex modes and scale
1/(2n) at BOTH ports. The intended right operators implement the required
parity dressing. If the algebra checks out, this supplies a stoquastic
flat parent and a classical-history interface, but not yet a uniform NN
parent. Next: derive both maps, check small examples, then examine the
locality and clock obligations without assuming a compiler theorem.

Substantive step 22 completed.

`classical_readout.md` constructs a joint two-port code with nonnegative,
equal amplitudes. It uses two distinct signed physical letter pairs per
logical flavor and realizes n<=floor(q(q-1)/4) complex flavors at scale
1/(2n), with left Gamma and right i Gamma P. The flat parent is stoquastic.
With a physical grading, use bipartite pairs and obtain
n<=floor(q_+ q_-/2), still quadratic. This is an unrestricted block code,
not a uniform NN realization. It evades the low-boundary-rank obstruction
because all code vectors must be entangled once n>q-1.

A second exact construction uses an inactive port letter and classical
permutation readout histories. Its complete zero code and boundary maps
are known, and its gap is 1-cos(pi/(2C)) for C clock values. A scheduled
local reversible circuit implements the clean readout interface using
swaps, parity computations on even ancillas, controlled cycles within
port grades, and a final partner-port exchange. The global clock remains
nonlocal. Converting the schedule to a repeated NN controller, excluding
malformed histories, and proving all-ring IR behavior remain open.

`classical_readout_checks.py` passed six static-code tests, including
(q,n)=(5,5),(6,6) with n>q-1, and two complete global-clock spectra at
dimensions 160 and 900. Port errors below 4e-16, full projected-bond error
below 1e-15, full-spectrum error below 1.2e-14. No CFT claim is made.

Next: audit whether clean readouts can be connected by a repeated local
FF Hamiltonian without introducing extra complete-code states. If that
stalls, examine a distinct limitation of stoquastic FF local parents,
without treating nonnegative amplitudes themselves as an obstruction.

Substantive step 23 completed.

The vacancy extension now also works for a frustration-free core
sum h0 + a_left + b_right with positive one-site boundary penalties.
Give the first/last occupied particles L/R tags and the others B tags,
for local dimension 3q+1. The repeated bond contains the oriented terms
a_L tensor I + I tensor b_R; each boundary penalty is then counted
exactly once even when vacancies touch the physical endpoints. Within
the fixed ordered-tag sector, the complete code is chi_p tensor V,
with the same endpoint compressions and scalar cut transport as before.
Other tag and particle-number sectors are not selected by this theorem.

A quantitative gap bound is proved as well. If delta is the charge gap,
gamma the initialized core gap, p_min the least adjacent-particle
probability, and M bounds the positive core interaction, the lifted gap
is at least delta p_min gamma / ((sqrt(delta)+sqrt(M))^2+p_min gamma).
This removes the previous need to assert a positive gap only by finite
dimension, but it gives no bound uniform in core length.

`tagged_vacancy_checks.py` directly assembles the physical bond in four
ordered-tag sectors, dimensions 108 and 270, using a q=3 symmetric core
with two unequal boundary penalties that remove its third color. The
complete code has dimension four. Code/port/transport errors are below
8e-16; all measured gaps exceed the proved bound. This is a local-parent
check, not a high-channel example or a CFT measurement.

Next: write the proof, audit the initialization-sector limitation, and
investigate whether the cheapest controller can avoid spurious cyclic
histories. The fifth-run results remain nonconclusive.

Substantive step 24 completed.

`length_filter.md` proves a new occupied-number selector. A qubit twisted
ferromagnetic core with a pinned left |+> and a right phase (k-1)phi,
phi=2pi/(k+2), is FF exactly at lengths r=k mod(k+2). It therefore
accepts only r=k among 2<=r<=k+2. A proved pinning/gap estimate gives
an explicit rejection energy epsilon_* of order k^-4.

Tensoring this filter onto any initialized FF payload core and using the
tagged vacancy extension selects exactly k occupied sites on BOTH physical
lengths k+1,k+2, within the union of ordered one-L/one-R sectors. Choose
J_filter>=2U/epsilon_* and 0<w<U/(k+1). Every wrong occupied number is
at least U-w(N-1) above the target code. The local dimension is 6q+1,
and the original endpoint maps and scalar cut transport survive. This
resolves the prior occupied-number gap conditionally on the tag sector,
without solving its selection in an unrestricted ring.

`length_filter_checks.py` passed all lengths and number sectors for
k=3,4, N=k+1,k+2, up to dimension 240. Measured selection margins
.99135,.89036,.99010,.86476 exceed the rigorous .7,.6,.6,.5. Maximum
code residual below 4.4e-13. No CFT calculation was performed.

Next: assess tag-sector enforcement and cyclic-history holonomy. Do not
claim the attractive finite-cell vacancy model is itself a homogeneous
conformal charge fluid. High-channel spatial cores and finite-coupling
all-ring IR control remain the main unresolved issues.

Substantive step 25 completed.

`arrow_wrapper.md` removes the remaining sector assumption INSIDE an
open cell. Use three arrow vertices, occupied letters L:0->1, B:1->1,
R:1->2, and three vacancy loops. Legal paths have at most one L and R.
Vacancy hopping changes its type when it passes an endpoint. Rewards
-mu per endpoint select both; K times the matching-violation projector
excludes every illegal sector. The violation projectors commute with
all hopping terms. With mu=2U,K=6U, the length-filter theorem gives the
COMPLETE desired code on both lengths k+1,k+2 at local dimension 6q+3,
without external sector restrictions or boundary terms. An explicit
full gap bound is included.

Moving the cut now requires |H_after><H_before| on the crossing vacancy.
The identity overlap is zero because those vacancy types are orthogonal.
With the conversion, the overlap is again beta times the identity on
both payloads. This conversion is local and can accompany a two-site
marker hop.

`arrow_wrapper_checks.py` passed complete legal-sector spectra of
sizes 315 and 961, containing all tag and occupation sectors at physical
dimension 9 and lengths 4,5. Ground-energy errors below 7e-13, code
residuals below 2.6e-13, converted transport error below 3e-16.
All length-five words obey the endpoint/violation count; overlapping
violation commutators vanish. Illegal sectors are excluded analytically.

A direct periodic use of this wrapper FAILS, proved: a closed legal
arrow path cannot contain L/R, and illegal sectors cost (K-2mu)v>0.
The all-vacancy ring is a zero-energy ground state. A selector must
supply arrow resets at cell boundaries. Thus this is a complete finite-
cell wrapper, not a microscopic ring counterexample.

Next: final audit and integration into RESULT.md. The leading remaining
construction problem is now a high-channel initialized uniform core,
possibly via the new classical readout controller, followed by combining
the wrapper with the mobile selector and controlling its infrared phase.

Substantive step 26 completed.

The graded classical readout code now also has independent logical
Z_f fields from diagonal TWO-SITE physical operators. Properly color
the selected bipartite port edges by color(a,b)=(a+b) mod(q/2), put
that color in the adjacent readout register, and assign +/-1 to the
two edge endpoints. Proper coloring makes the supports for distinct
flavors disjoint. Compression is Z_f/(2n), and the norm of a combined
field is max_f |mu_f|. Together with the arbitrary real boundary
Majorana hopping matrix this reaches independent hopping phases and
chemical potentials at first order. The spatial controller must gate
the fields to completed readouts. Colored-code checks at (4,2),(6,4)
pass below 3e-17, with port and density identities checked separately.

The live `contrib/ir_control/` report finished during this run. Its
verdict is still heuristic: exact grading and well-separated Fermi
points are useful, but general-cell leakage and a full conformal
spectral conclusion remain open. This run did not independently audit
its literature or all of its proofs. No newly supplied theorem closes
the main question.

Fifth-run audit completed. The new finite-dimensional arguments were
checked for port signs, complete-code hypotheses, initialization,
charge-gap dependence, all competing open-cell sectors, and the
necessary vacancy-type conversion. The periodic wrapper failure is
proved and is recorded explicitly. All four numerical scripts passed;
the classical-readout script was rerun only after adding density
readouts and removing a dtype warning. Compilation passed. Earlier
suites were not rerun. No whole-class conclusion, independent referee,
network use, cluster job, commit, or push. All writes stayed here.


## Sixth run, completion record

Re-read the requested context and live contributions. The first attempted
line, turning the scheduled readout circuit into a repeated-bond controller,
still lacked program selection and exclusion of all unwanted components.
No implementation theorem was assumed. Moved to the independent local
composition and ring-spectrum problems rather than stopping there.

Substantive steps 27-31 are complete: selector/wrapper composition,
full projected marker spectrum, helix port equalization, combined
payload/ring spectral reduction, and positive approximate-core lifting.
The five corresponding notes and scripts are listed at the top. They
are reproduced in RESULT.md Sections 30-34, with audit and ranked next
steps in Section 35. The previous steps remain in the chronology above.

Final source audit checked all sector bounds, transfer orientations,
auxiliary phases, ring momentum constraints, and approximation bounds.
All five new numerical scripts and compilation passed. A final read of
the live ir_control report found new toy numerics but no theorem closing
the general-cell leakage or spectral step. The whole-class verdict is
still nonconclusive. Only this folder was modified.

## Seventh run, in progress

Read the saved state, STATUS, requested context, minimal_d0 and C13's
live progress. C13 is treating finite-coupling length/payload resonances;
this run will not duplicate it. No external source or network used.

Substantive step 32, explicit clock construction derived, checks running.

A nine-letter repeated NN clock now has exactly ONE legal path at every
open length n. Its C_n=n+4(2^n-1) vertices include a dedicated n-step
initialization sweep, then a reversible binary carry/return counter.
All locally legal configurations lie on that path. The local rules
preserve EACH legality-defect projector, including boundary defects.
There is no uninitialized clock cycle or independent internal tape length.
The first implementation is `binary_clock_checks.py`.

Unitary decorations of its transitions can act on a fixed work register.
Pins during the idle initialization sweep give an exact positive core
with the complete initialized history code and an explicit finite-length
gap. This initializes a product workspace; it does NOT yet implement
or validate the high-channel readout program. A counter is not a compiler.

The hopping version has energy -2J cos(pi/(C_n+1)). Its shifted energy
appears strictly convex in n, so it can select cell lengths without a
separate walker. The proof, illegal-sector comparison, and full checks
are being completed before recording that as established.

Substantive step 32 completed.

`initialized_binary_clock.md` proves the complete nine-letter graph,
positive initialized core, gap bound, and an optional merged frozen
selector. The shifted hopping energy is strictly convex, so every
interior length n>=2 is selectable with one chemical potential. The
full ring illegal-sector comparison is included. This does not repair
frozen-ring conformality and is not used as such.

`binary_clock_checks.py` passed. Exhaustive clock graphs at n=1,...,5
cover up to 59049 words and all their components. Every local legality
projector is conserved. Complete illegal spectra through n=4 obey the
bounds; both legal spectra through n=8 (1028 vertices) agree below
4e-15. A 224-dimensional entangling-work example has exactly four
initialized zero states and gap .0062877901, above .0012755102.

Substantive step 33 completed, full-sweep extension.

The sparse counter visits the far end too few times for a straightforward
interpreter. A 15-letter extension inserts an outward and return sweep
before EVERY increment. Its complete legal graph is again one path,
now C_n=(2n+4)2^n+n-4 vertices, with exactly 2^n full sweeps. The
initialization and illegal-sector proofs transfer. Each sweep executes
one partitioned local update U=(tensor F_j) product G_(j,j+1) when G
acts on disjoint subregisters. Arbitrary boundary F_j are available
through boundary rules. This is an exact gate interface, not yet the
readout program.

`sweeping_clock_checks.py` passed complete graphs through n=4 (50625
words), full legal hopping spectra through 1026 vertices, and an
entangling three-site update on a 64-dimensional work space repeated
eight times, with error below 8e-15. It explicitly includes thousands
of illegal components with cycles, all separated by conserved penalties.

Next: give the mechanical reversible-machine-to-sweep reduction, then
try to complete the actual initialized readout program. Avoid treating
an unspecified universal machine as the missing explicit construction.

Substantive step 34 completed, mechanical work interface and readout audit.

`sweeping_clock.md` gives a canonical graded permutation completion for
an injective local collision table. Transport lanes plus onsite
collisions implement a specified reversible packet/tape machine on the
full work space. The clock does not require every work-machine orbit to
be valid: input pins select the complete history kernel. Later diagonal
verification checks for a classical permutation program retain the
explicit path-pin gap bound, even at input-dependent checking times.

`sweep_machine_checks.py` builds a literal 18-state work machine and
checks every initialized input through five sites. A left endpoint
observable compresses to alpha times total logical parity. Its
accumulator-X compression is exactly zero, because the original logical
bits remain in the tape. This confirms the clean-readout requirement;
it is not a high-channel code or a new general boundary-complexity result.

`readout_machine_contract.md` isolates the still missing program. It
states a sufficient finite collision-table contract, proves the ensuing
complete high-channel core conditionally, and gives the dimension count
q_core=K(q+1), d=96K(q+1)+48. K is NOT supplied, and the resource inequality
is NOT a counterexample. A generic reversible simulator retaining copies
of the exposed bit is insufficient. The active flag must follow complete
cleanup, and the schedule must give input-independent readout counts at
partial sweep times as well as full updates.

Substantive step 35 completed, full periodic selector check.

`clock_ring_checks.py` exhausts the entire ten-letter ring Hilbert space
at L=3,4,5, dimensions 1000,10000,100000. Selected cells have interiors
2,3,4 and ground multiplicities 3,4,5. Energies agree with the analytic
formula within 3e-16; all illegal components obey the proved bound.
Some illegal components are longer than the legal clock path, but their
conserved defects exclude them. This is a frozen selector, not a CFT.

Latest C13 PROGRESS was read after its final update. It reports a
conditional positive Gaussian-phase analysis, an all-integer fixed-vacuum
qualification, and finite-size generic numerics that remain too small
for accurate central charges. No part of that analysis was duplicated
or promoted to a complete microscopic theorem here.

Next: audit and integrate the results into RESULT Sections 36-39. The
actual clean-readout collision table remains unfinished. The next run
should start there, not reconstruct the now explicit clock.


Seventh-run final audit completed.

RESULT Sections 36-39 and the top of this file now hold the complete
handoff. The source audit checked the clock boundary seed, all legality
constraints, positive core terms, both convexity calculations, ring
fragment bounds, local update order, and the clean-readout condition.
Four new scripts passed and compiled. No additional machine compiler or
microscopic conformal result was claimed. The result is nonconclusive
under the user's stated definition, with the remaining step localized
in the finite readout program. Only this folder was modified.


Handoff persistence note. A concurrent repository update temporarily
restored older tracked versions of RESULT and PROGRESS during final
integration, then restored the new versions. The final seventh-run
contents are also saved as RUN7_RESULT.md and RUN7_PROGRESS.md. These
snapshots preserve the handoff if another synchronization touches the
shared tracked files. No external repository files were changed here.
