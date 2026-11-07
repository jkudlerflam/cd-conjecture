# C134 result

2026-11-07. The full interval -1<Delta<=1 remains unresolved. An explicit theorem on the attractive half of the interval is proved in [PROOF.md](PROOF.md), including odd lengths. The proof is new in this run and has not had an independent referee pass.

## A compact-parameter theorem

[Proved] For every 0<delta<=1, every Delta in [-1+delta,0], and every integer L>=2, with the half-filled minimizing-flux convention in the question,

    F_L >= L e_0 - 1200 J/(delta^3 L),
    |F_L-L e_0| <= 1200 J/(delta^3 L)              if L is even,
    L(F_(2L)-2F_L) <= 3000 J/delta^3.                         (1)

Thus C(delta)=3000/delta^3 is one explicit answer on this range. PROOF.md Section 1 gives substantially smaller parameter-dependent constants and an explicit compact majorant. All constants are independent of L. The estimate is not uniform as Delta approaches -1.

[Proved, mechanism] For Delta<=0 the scattering kernel is negative and has L1 norm below one. This makes the finite counting function globally increasing. Its roots are midpoint nodes, and the cumulative root error is bounded by B/L using a Neumann inverse. Dressing the energy gives an exact quadrature identity. Clipping the dressed energy where the infinite-volume root density reaches 2D/L produces a function whose derivative has bounded variation independent of L. Midpoint error and the two tails are then O(1/L) in total energy. The finite Bethe ground branch, including odd lengths, is identified by strict convexity of the Yang-Yang functional, a coordinate Bethe calculation near the free point, analytic continuation of the characteristic polynomial, and Perron-Frobenius simplicity.

[Proved] Odd lengths use n=(L-1)/2, or its spin flip, but F_(2L) uses n=L. These are not C129's matched-density sectors (2L,2n). The extra counting-coordinate strips at odd L have the favorable sign for the lower bound.

## Exact points

[Proved, inherited free result with corrected odd specialization] At Delta=0, the sharper constant is C=4-2sqrt(2) in units of J. The odd half-filled doubling defect is zero.

[Proved using a published exact-energy identity] At Delta=1/2, C=2pi^2/9 works for all L. The even ring has energy -3JL/8 at the supersymmetric twist 2pi/3. Comparing this with the minimizing twist costs at most pi^2 J/(9L), while flux minimization gives the needed upper bound on the doubled ring. At odd L the minimizing-flux energy itself is -3JL/8 and the defect is nonpositive. PROOF.md Section 11 checks the literature sign convention and cites the exact identity.

[Proved, inherited] At Delta=-1 the defect is zero. C129's band (Delta+1)L^2<=A also bounds the present half-filled defect by A in units of M, including odd L. This size-dependent band does not close the gap to (1) uniformly in Delta.

## Literature and remaining interval

[Proved, literature audit] The finite-size expansion in Kozlowski 1508.05741 does not include the half-filled massless endpoint. Theorem 3.5 has convergence without a rate; Corollary 4.3 has the energy term below half filling. DKKMT 2012.11675 excludes exact half filling in both of its relevant error estimates. Details, checked citations, and access limitations for the older papers are in [LITERATURE.md](LITERATURE.md).

[Assumed, unresolved] No theorem for a positive-anisotropy interval, or for Delta=1, was obtained. A common bound for Delta<1 would extend to Delta=1 by finite-matrix continuity. Parameter-dependent bounds that diverge at 1 do not suffice. Monotonicity of D(L) alone cannot produce the missing factor 1/L. The new real-axis proof uses a kernel sign unavailable at positive Delta. The remaining tasks are the positive interval, the uniform ferromagnetic crossover, and an independent audit of the new compact theorem.

## Numerical results and computational limit

[Numerical] The requested eight anisotropies were checked in 152 length/anisotropy pairs. ED covers every L=2,...,12, both neighboring sectors at odd L. Bethe calculations cover L=2,...,12 and L=16,17,32,33,64,65,128,129, with doubled rings through 258 sites. There were no solver failures. The maximum Bethe/ED energy discrepancy is 8.44e-15. The largest sampled Q=L(F_(2L)-2F_L)/M is 4 at L=2, Delta=1. All applicable proved ceilings pass. These are floating-point checks, not certified enclosures.

| Delta | Q at L=128 | Q at L=129 |
| :--- | ---: | ---: |
| -0.999 | 0.032583757 | 0.031618077 |
| -0.99 | 0.115254793 | 0.104807419 |
| -0.9 | 0.399580519 | 0.284785026 |
| 0 | 1.570934329 | 0.000000000 |
| 0.5 | 2.040785996 | -0.680192574 |
| 0.9 | 2.387225369 | -1.685155117 |
| 0.99 | 2.463203882 | -2.081140331 |
| 1 | 2.471475390 | -2.134273930 |

[Numerical] Separate checks verify the dressing equation, the exact dressed-energy quadrature identity, the scattering-phase identity, cutoff variation bounds, and supersymmetric twisted energies. See [CHECKS.md](CHECKS.md) for commands and errors.

[Assumed, pending computation] The requested several-thousand-site scan was not run. [run_typhon.slurm](run_typhon.slurm) prepares an eight-job array through L=4097 and doubled length 8194, retaining both parities, solver failures, and cancellation diagnostics. Typhon is unavailable from this sandbox, so no job was submitted. All writes remain in work_C134. [STATE.md](STATE.md) is the local handoff; the project-level memory was not modified under the task's outside-write prohibition.
