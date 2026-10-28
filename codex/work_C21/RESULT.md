# C21. Local source resolvents and the remaining spectral compiler

## Verdict

[proved partial result] A fixed 16-state site now supports an explicit
nearest-neighbor isolated-cell Hamiltonian with an exactly degenerate
2^n-dimensional code, exact parity, and a nonscalar source resolvent.
Its second-order exchange is exactly Gaussian and has odd hopping rank
2n. The source spectrum has between 2n+1 and 6n distinct accessible
energies. This demonstrates that local repeated terms can overcome the
restricted scalar-resolvent obstruction and carry growing spectral rank.

[proved limitation] This cell does not realize C19's requested target.
Its hopping matrix has the form J=-2 R Gf, with R a reversal and Gf
positive definite. At zero field its Fermi points can occur only at
0 or pi; the tested one-flavour example is gapped. Independent displaced
flavours and chemical potentials have not been constructed. The final
translation-invariant rule table, marker-sector exclusion, and d0 for
the desired model therefore remain open. Neither this failure nor the
other obstructions below kills route A.

The full local rule table and proofs are in [LOCAL_CURSOR.md](LOCAL_CURSOR.md).
The cell's 16-state module count must not be used as a completed d0.

## Additional finite results

[proved] C19's original Jacobi weights cannot be realized unchanged
on a product-subspace code in a nearest-neighbor cell of length at
least four. Its left and right first source moments fail to commute,
whereas product-code compressions of operators on the two disjoint
boundary neighborhoods commute. This remains an obstruction after
an arbitrary common logical basis rotation. Entangled history codes
are outside the statement.

[proved] Generalized Jacobi filters repair this particular early-moment
problem. The functions

    x(1-x)^s P_j^(2s,1)(2x-1), x=exp(-u), j=0,...,n-1,

after normalization, use k=n+s poles, have all moments through s-1
zero, and retain C19's Cauchy orthonormality. A delayed field filter
gives the same independent hopping, phases and chemical potentials at
the abstract matrix level, with exact parity. The bound k>=n+s is
necessary for simultaneous delay s. This is an algebraic repair,
not a compiled Hamiltonian. Details are in
[MOMENTS_AND_FILTERS.md](MOMENTS_AND_FILTERS.md).

[proved] An eight-state local excitation ladder has n+1 distinct
endpoint-accessible poles, but its source weights are scalar. Arbitrary
unitary decorations on a path or rooted tree do not change that fact
when the on-site energies and edge magnitudes are scalar. This excludes
a tempting implementation even with an arbitrarily long ladder. The
Clifford cursor escapes it by putting logical operators in the excited
on-site energies. See [SCALAR_LADDER.md](SCALAR_LADDER.md).

## Numerical checks

[numerical] All computations used Python, NumPy and SciPy inside this
folder. The new cursor's physical one-excitation matrices have dimensions
12,48,144 at n=1,2,3. Direct spectral decomposition and two-cell source
contraction give the following results.

| n | Accessible energies | Non-Gaussian second-order residual | Smallest hopping singular value |
| --: | --: | --: | --: |
| 1 | 6 | 2.24e-16 | 4.19e-5 |
| 2 | 12 | 9.61e-16 | 1.55e-9 |
| 3 | 18 | 4.97e-15 | 5.75e-14 |

Full rank is proved by distinct first nonzero Taylor orders of the
source filters. The n=3 smallest singular value is too small to use
as a reliable standalone numerical rank test.

[numerical] Full finite-coupling cursor ED retains all filler and
multiple-cursor states. The comparison uses the canonical effective
matrix from the polar overlap with the exact code. With Delta=3,
t=1, g=0.6 and lambda=0.04,0.02,0.01,

| n | Cells | Full physical dimension | Code levels | Error at lambda=0.01 | Halving exponents |
| --: | --: | --: | --: | --: | :-- |
| 1 | 2 | 256 | 4 | 2.07457e-7 | 3.0035, 3.0018 |
| 1 | 3 | 4096 | 8 | 4.14423e-7 | 3.0007, 3.0005 |
| 2 | 2 | 65536 | 16 | 2.09823e-7 | 3.0037, 3.0019 |

These errors compare with the cursor bond J=-2 R Gf, not C19's target.
Eigenpair residuals are at most 1.82e-12. The n=2 calculation uses both
complete parity sectors with eight code levels in each; it is not a
single-excitation truncation. Its separation from the other bands uses
the analytic gap bound, while the n=1 next gaps were measured directly.

The initial scalar ARPACK attempt at n=2 stalled and was interrupted.
A block LOBPCG calculation initialized with the complete code succeeded.
The final script places a 180-second limit on each cluster calculation.
No n=2 three-cell or n=3 finite-coupling cursor ED is claimed. The
two-cell n=3 source resolvent was evaluated, and the absence of
second-order mixed-bond terms on arbitrary chains is proved analytically.

The cubic residual does not establish that all very weak hopping
channels are resolved at these couplings. In particular, errors can
exceed their small logical energy scales. No finite-coupling CFT or
stability conclusion is drawn.

[numerical control] The scalar ladder was tested on two and three cells
at n=1,2 and on two cells at n=3. Its largest complete symmetry sector
has dimension 65,536. Direct resolvent errors are at most 1.2e-15;
finite-coupling exponents are 3.0021 to 3.0038. These results confirm
the scalar failure. They do not fill the missing three-cell target ED.

[exact arithmetic and numerical] Delayed-filter normalization and
moment identities passed twenty Fraction checks, through n=8,s=8.
Direct target contractions at n=1,2,3 and s=1,2,3 have absolute errors
below 5.6e-16, with exact numerical parity commutators. The worst error
relative to rho^2 is 9.4e-7, reflecting filter conditioning.

## Task accounting and continuation

| Requested task | Result |
| :-- | :-- |
| Local operator-valued interface reproducing the target | Local full-rank Gaussian interface proved; required flavour structure remains open |
| Growing accessible poles on a fixed alphabet | Proved for both local modules; C19's weights remain uncompiled |
| Full homogeneous rule table, initialization, illegal sectors, d0 | Complete isolated-cell tables and all internal-sector bounds supplied; ring selector and target d0 absent |
| Two/three-cell target ED at n=1,2 | Cursor checks completed at (1,2),(1,3),(2,2); scalar controls also at (2,3); no target validation |

[proposed next step] Change the relative left/right source geometry of
the Clifford cursor while retaining its exact Clifford-vector resolvent.
The immediate test is the resulting Majorana hopping matrix, including
its Bloch zeros and available chemical potentials. The form R Gf must
be changed before spending effort on a selector. A candidate that
instead compiles C19 directly must evade the boundary-moment obstruction
through entangled encoding or altered source weights.

[scope] C8's clock and C18's richer boundary checks were reused as
context, not rederived. The collaborator's graded repair supplies the
parity and contraction conventions; the sparse first-order work has a
growing alphabet and does not implement these filters. C22 independently
obtained the tree-gauge observation and is studying a different
charge-conserving bus. Its later progress and C23's completed assembly
audit were read. Their bus and infrared calculations were not duplicated
or treated as proved target implementations here.

## Files and reproduction

The derivations are [LOCAL_CURSOR.md](LOCAL_CURSOR.md),
[MOMENTS_AND_FILTERS.md](MOMENTS_AND_FILTERS.md), and
[SCALAR_LADDER.md](SCALAR_LADDER.md). The runnable checks and saved
records are

- [clifford_cursor_checks.py](clifford_cursor_checks.py), [clifford_cursor_results.json](clifford_cursor_results.json).
- [filter_checks.py](filter_checks.py), [filter_results.json](filter_results.json).
- [ladder_checks.py](ladder_checks.py), [ladder_results.json](ladder_results.json).

Run each script with

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 SCRIPT.py

PROGRESS.md records the investigation and the final handoff. No files
outside work_C21 were changed, including the project STATE.md, because
the task explicitly forbids those writes. No network, cluster job,
subagent, external message, commit, push, or independent referee was
used. This is a partial construction with a proved limitation, not a
completed constant-d counterexample.
