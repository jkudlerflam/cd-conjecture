# C8 progress

2026-10-27. No conclusive answer. The completed report is `RESULT.md`.

## Current best answer and restart point

The abstract first-order Majorana window remains open: both endpoints can realize the required scaled Clifford operators on the same code. Two explicit constructions and their proofs are in RESULT.md Section 1. The simplest one-pass history implementation fails (Section 3). Second-order exchange works with identical unrestricted blocks and independent left/right ports (Section 4), but exact large microscopic symmetry protection is excluded at all orders (Section 5).

No microscopic small-d bond implementing the code was found. No bound c<=f(d) for the full class was proved. The free target is conformal; the encoded microscopic model has not been shown to be conformal.

Open problems, ranked for the next run:

1. Local realization of the explicit joint code using a repeated nearest-neighbor bond. Try multiple sweeps or a programmable clock, with physical boundary operators controlled throughout. Do not redo the joint channel dilation or the one-pass obstruction.
2. Prove the local effective-Hamiltonian reduction and the incommensurate target stability assumptions in RESULT.md Section 6. The rank-one d=20992 construction there would then refute c<=r(h), although not c<=d-1. This route avoids the internal-block locality problem by treating the whole block as one physical site.
3. If retaining the small-d goal, realize the joint boundary spectral weights of Section 4 locally and control thermodynamic errors. G1 and this run already remove the abstract second-order rank and alternating-block problems.
4. Resolve forced-defect sectors on microscopic rings not divisible by the selector period.
5. For a general bound, find a locality restriction on these explicit codes or resolvents. Port dimension, complete positivity and perturbative insertion counting alone are insufficient.

Validation completed: `python3 c8_checks.py --ed`, with single-thread BLAS, reproduces `checks.json` and `two_port_N3.npz`; all reported checks passed. `python3 -m py_compile c8_checks.py` passed. No full microscopic CFT numerics, cluster job, literature search, or independent referee check was performed. The calculations are finite-dimensional consistency checks and exact logical-target calibrations.

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
