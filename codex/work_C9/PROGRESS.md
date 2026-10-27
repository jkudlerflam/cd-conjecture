# C9 progress

Working closed book. No pre-existing files are read or listed. All calculations and files stay in the current directory. The supplied ingredients are treated as premises, not independently rederived.

## Initial plan

1. Derive necessary restrictions on fixed-dimensional boundary ports and compare first- and second-order exchange.
2. Examine whether finite programmable cells can be supplied without assuming a translation-invariant simulation theorem whose boundary-operator guarantees are stronger than stated.
3. Test explicit completely positive scaled-Clifford boundary encodings numerically, if constructible, and analyze their exact first-order dynamics and corrections.
4. State a construction only if microscopic locality, exact degeneracy, error control, and infrared stability all close. Otherwise locate the obstruction and separate proved no-go statements from unresolved steps.

No external references, project files, or network sources will be consulted. The user explicitly requests closed-book work, so project-memory instructions are superseded for this task.

## First deductions

The first-order free-Majorana route needs two different port compressions. With on-cell Clifford generators gamma_a and parity P, the right compression must be i gamma_a P while the left compression is gamma_a (up to orientation/sign). Using gamma_a at both ends gives a spin interaction, not the desired Jordan-Wigner free hopping. An odd number of generators with a genuine on-cell parity introduces an unused Majorana and extensive zero-mode degeneracy unless another mechanism removes it. The direct ladder should use an even number of generators.

A joint, two-port completely positive map can realize both sets at sufficiently small scale. This solves an abstract finite-dimensional encoding problem, not the translation-invariant nearest-neighbor cell construction. I will build this encoding and test its first- and second-order spectra.

There is a prospective exact-symmetry obstruction: an exact on-cell SU(N) symmetry that preserves a full q-dimensional port exchange normalizes the port matrix algebra. For N > q its induced projective q-dimensional representation is trivial. The port then commutes with the on-cell symmetry, which makes all induced interactions scalar on fundamental code spaces. This would rule out a common symmetry-protection strategy, not the requested family in full generality.

## Verified abstract cells

Created checks.py and ran it with Python 3, NumPy, and SciPy. The shell has no `python` executable; `python3` works. Numerical runtime is under one second per run with BLAS threads limited to one.

For two qubit ports, a joint CP map gives exact scaled left gamma_a and right i gamma_a P compressions. The r=2 example has code dimension N=2 and full abstract cell dimension D=32. The r=3 example has N=4 and D=64. Both use scale s=0.1 and H_A=I-GG^dagger, with exact abstract gap 1. Isometry and compression errors are below 1.2e-15. Exact low spectra of two cells (dimensions 1024 and 4096) agree with second-order perturbation theory with O(lambda^3) residuals.

For r=2, the second-order term is explicitly a scalar plus (3 s^4/2) T^2, where T is the desired hopping bond. It is not just a renormalization of the free hopping. For r=3, the second-order correction also violates the total logical parity. Four open r=3 ladder cells have the predicted fourfold ground-state degeneracy from unused Majoranas. This is an obstruction to using that odd-channel code as a pure three-channel ladder.

These cells are finite-dimensional matrix constructions only. No local implementation inside a spontaneous cell has been supplied.

## Further deductions

The exact-symmetry obstruction extends from full port swap to any two-cell interaction supported on q-dimensional ports, if its symmetry is an exact product of identical on-cell SU(N) actions with fundamental code spaces. The minimal operator-Schmidt spaces generate subalgebras of M_q. SU(N), for N>q, acts trivially on these algebras. Independent cell SU(N) symmetries then survive, prohibiting every nonconstant logical interaction.

A single port in the supplied G1 gadget is not automatically two geometric boundary ports. This finite-dimensional issue can be repaired by using two independent maximally entangled port-ancilla pairs and applying the stated spectral-weight lemma to the direct sum of their syndrome spaces. It does not solve nearest-neighbor programmability.

Error smallness at fixed N is insufficient for an infinite critical chain. A finite relevant coupling produces a finite crossover length, and a marginally relevant coupling can produce an exponentially large but finite one. Choosing finite energy scales independently of L cannot remove this issue.

## General second-order bond synthesis

Derived a stronger finite-dimensional result. With two qubit ports and a code of dimension N, any Hermitian logical bond K can be obtained, up to a scalar and an overall positive scale, at second order. Choose k distinct excited energies and the positive Cauchy matrix C_ef=1/(e_e+e_f). A subspace orthogonal to both the all-ones vector and C p cancels POVM normalization errors and first-order-in-weight one-cell fields. Choose C-orthonormal columns Z in this subspace. Opposite independently chosen left/right spectral-weight perturbations then give any real coefficient matrix K in a Hermitian operator basis. Positivity follows by making their common amplitude sufficiently small. A Naimark dilation of the combined left/right syndrome space realizes the weights.

This also removes the antiferromagnetic-sign ambiguity at the abstract-cell level. It does not supply microscopic nearest-neighbor cells or infrared protection.

Created spectral_cell_check.py. A D=40, N=2 cell with two qubit ports, excited energies 1, 4, 16 (and spectator energy 32), realizes the directed logical bond X tensor Z. The predicted second-order coefficient is 0.0006363304136605024 and the scalar is -0.42463235294117646, in units Delta=1. Direct spectral-weight and effective-Hamiltonian errors are below 1.8e-15. Exact two-cell diagonalization gives cubic residuals with error/lambda^3 approaching 0.1106. Runtime was about five seconds.

## Infrared alternatives and qualifications

A stable Narain target would require a lattice realization in an open gapless region and a strict gap above dimension 2 for nonzero-charge vertex operators. The supplied phrase 'no primary of dimension <=2' needs the current-algebra meaning: Narain theories always have dimension-one currents and dimension-two Gaussian deformations. These are vacuum descendants for the extended current algebra, not all absent Virasoro primaries.

Exact logical fermion parity can be built into the abstract CP construction by using only physical port operators odd under a physical grading. Their real dimension is 2 q_+ q_-, giving c=q_+ q_- at first order. This removes the odd-generator and parity-breaking defects of the ungraded maximal construction, but does not remove potentially dangerous interactions or solve programmability.

## Final status

RESULT.md is complete. It contains a proof of arbitrary logical bond synthesis through two qubit ports at second order, a conservative D = 4 ceil((6 N^3 + 13 N)/4) whole-cell dimension bound, an exact on-cell SU(N) protection no-go for N larger than the boundary Hilbert dimension, joint CP Majorana encodings, explicit correction calculations, conditional clock/error scales, commensurability, and a requirement-by-requirement gap statement.

No fixed microscopic d0 or translation-invariant nearest-neighbor program rule has been constructed. No general impossibility theorem and no c-d counterexample are claimed. The finite-dimensional constructions and their checks do not establish thermodynamic criticality.

The completed numerical checks use D=32 and D=64 CP cells and a D=40 spectral-weight cell. All effective-matrix checks agree at roundoff, and all tested exact two-cell spectra have the predicted cubic residual after second order. Eigenpair residuals are at most about 1.1e-12. The r=3 parity violation and unused-mode degeneracy are failures of that candidate route, recorded explicitly rather than treated as successful criticality checks.

Python syntax checks pass for both scripts. Total numerical runtime was seconds; the wall-clock work remained well below 90 minutes. No pre-existing project files, data, state, or literature were read or listed. Project STATE.md was not accessed or modified because the closed-book and current-directory restrictions supersede project-memory instructions.
