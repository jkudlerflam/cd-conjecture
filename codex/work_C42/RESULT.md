# C42 result

[proved, finite-dimensional result] The modified 46-letter initialized payload core has an all-n nonresonance construction. For every n, finite nonzero parameters give exactly 2n simple Fermi points in the quadratic projection of the exact second-order bond, distinct sea centers modulo pi, and a strictly positive margin against every nontrivial number-conserving quartic resonance. Repeated endpoint pairs and umklapp are included. Exact U(1), parity, the complete 2^n-dimensional degenerate unperturbed code, and zero first-order port coupling are preserved.

[proved, attribution] This uses C39's two-bus core and biased-vacuum readout. It is not a solution within the original 23-letter alphabet. C42 gives the detailed singular-pencil argument, a terminating rational parameter schedule, explicit finite-express and finite-conversion error bounds, and local reruns of the finite-cell checks. The proof is in [PROOF.md](PROOF.md). No independent human review is claimed.

[assumed, remaining infrared step] This establishes the requested finite payload input. It does not establish an interacting central charge, a finite-coupling vacancy transfer theorem, or a counterexample to c<=d-1.

## Why the all-n argument works

[proved] The two virtual buses give the reference hopping

    T=-Gamma(w)-i eta Gamma(wv),

where Gamma is a positive filter Gram matrix. After a congruence that retains every cross entry, the zero-field centers are -atan(eta r_j), where r_j are generalized eigenvalues of two different positive Gram matrices. Their scales are r_j=kappa_j v^(4j)(1+o(1)), kappa_j>0. Unlike the single-bus bandwidth hierarchy audited in C38, these scales change the Fermi zeros themselves.

[proved] The biased vacuum compresses one repeated local density readout into a geometric field. In the v->0 eigenvectors its width derivatives have scales s_j=d_j u^(4j)(1+o(1)), d_j>0. Choose u first and then v. Both center magnitudes and width derivatives can have successive ratios below 1/32. A difference of endpoint pairs either has a nonzero center combination, whose first nonzero term dominates, or has the same flavor multiset and a nonzero width derivative, whose first term also dominates. A sufficiently small nonzero field splits every latter coincidence. This proves the nonidentity statement for every n.

[proved, constructive schedule] Fix w=1/4 and enumerate dyadic eta,u,v,mu with increasing rational precision. Accept a tuple only after exact Sturm root isolation and rational pair comparisons certify all 2n roots and all pair sums. The preceding argument proves termination. From the accepted brackets compute the positive rational tolerance tau in PROOF (7). With p=2n+4, r=1/16384 and B0=128, it is sufficient to choose

    D >= max(2/r, 8n B0(1+eta)/(tau r^(p+1))),
    0<g<=1/4,
    g^2 <= tau/[192n(1+eta)D^p].

These bounds preserve every bracket for the exact finite-g second-order quadratic projection. PROOF (8) supplies a positive velocity lower bound. PROOF (11) bounds the complete non-Gaussian remainder, which can be made arbitrarily small by reducing g further while retaining the same reference momenta. Lambda is chosen last for finite-cluster control. No coefficient is varied with the physical ring length.

[proved, computational scope] This is an explicit terminating algorithm in n, not a closed elementary formula for its search depth. The proof of termination is analytic. The implemented rational certificates and budgets were executed at n=1,2,3, not at arbitrary n. Their fixed witness values u=1/64,v=1/1024,eta=1/100 are not asserted to work for every n.

## Exact small-n certificates

[proved, computer-assisted, n=1,2,3] certify.py uses rational Gram matrices, Sturm chains and rigorous angle enclosures. It checks all unordered pairs, with replacement. The following conservative decimal lower bounds hold for both the reference and the finite-D, finite-g parameters produced by budgets.py, since the same root brackets survive.

| n | Simple roots | Endpoint margin greater than | Quartic margin greater than |
| :-- | --: | --: | --: |
| 1 | 2 | 3.14 | 5.11e-16 |
| 2 | 4 | 0.01660 | 1.23e-28 |
| 3 | 6 | 1.42e-14 | 2.36e-41 |

[proved] The fields are mu=10^(-8n-10). The exact parameter files budgets_n1.json, budgets_n2.json and budgets_n3.json contain integer D and dyadic g=2^(-1138), 2^(-1738), 2^(-2761), respectively. These deliberately loose certificates use log10 D approximately 102.15, 120.51 and 155.36. The decimals display the exact formulas, not floating-point sign tests. No microscopic ED is claimed at these extremely small physical scales.

## Two- and three-cell checks

[numerical, n=1,2,3] The literal 2116-by-2116 two-site matrix matches the legal-cell matrix to 5.9e-16, with zero computed charge commutator. The complete legal-cell dimensions are 54,160,424. The full-code second-order bond dimensions are 4,16,64. Its g^2 coefficient, Gram limit and relative bus phase are checked independently of the finite-lambda spectra by direct source contraction. All these numerical tests supplement the exact local and source identities proved in PROOF.md.

[numerical, 12 untruncated symmetry sectors] Two and three cells were diagonalized at n=1,2,3 in total charge Q=0 and Q=1, in the even sector of each bus count containing the code. Each calculation uses the entire invariant legal sector, with no energy truncation. The moderate test parameters are w=.8,v=.65,q=.8,D=.2,g=.4,a=3.5,eta=.6 and a lambda^2*.03 conditional readout. These test the implementation, not the asymptotic nonresonance witness. The table reports Q=1 operator-norm errors against the full-in-g second-order bond plus the compressed readout.

| n | Cells | Q=0 dimension | Q=1 dimension | Error/lambda^4 at .1 | Error/lambda^4 at .05 | Maximum Q=1 eigenpair residual |
| :-- | --: | --: | --: | --: | --: | --: |
| 1 | 2 | 265 | 338 | 0.03046069 | 0.03037385 | 6.63e-13 |
| 1 | 3 | 3483 | 5265 | 0.04485196 | 0.04478197 | 1.54e-12 |
| 2 | 2 | 1240 | 2224 | 0.02062424 | 0.02053909 | 4.99e-13 |
| 2 | 3 | 25984 | 60864 | 0.03370639 | 0.03360378 | 9.77e-13 |
| 3 | 2 | 3865 | 9110 | 0.01399300 | 0.01392244 | 6.81e-13 |
| 3 | 3 | 109565 | 356325 | 0.02373706 | 0.02363903 | 5.96e-11 |

[numerical, outcome] Every computed first-order compression is zero. The error divided by lambda^4 is stable as lambda is halved. Q=0 gives the same scaling, with coefficients between 0.0129 and 0.0438. The largest successful calculation took 126 seconds. Its eigensolver returned the nine-dimensional low-energy code subspace, not an additional excited eigenvalue, so no numerical gap is reported for that run. Numerical agreement is absolute; it does not resolve every tiny hopping coefficient relatively. No full 4464-letter wrapper spectrum or central charge was computed.

[numerical, failed run retained in the record] The first n=3, three-cell Q=1 attempt exceeded its 180-second limit and is not evidence. The accepted rerun uses a stated 8e-11 residual target and warms the second calculation with the first low-energy subspace. Its largest actual residual is 5.96e-11. PROGRESS.md records the failure and the changed solver settings.

## Vacancy completion and alphabet cost

[proved, inherited finite wrapper] The initialized core has 46 letters, but this count does not include a spontaneous selector. The existing positive C8/C32 wrapper gives a conservative complete internal alphabet of 96*46+48=4464. Its hypotheses are verified by the new core. The selector can be detuned to one length using the explicit inequality in PROOF Section 8. The wrapper's finite-resolvent error can be put below any already certified payload tolerance.

[proved, exact vacancy compatibility] Add one species-independent vacancy to that complete internal bond, lifting the whole bond, including all half-onsite terms. The specified vacancy squeezing identity applies to any internal alphabet and gives the same positive leading scalar multiplying all internal couplings. Thus a conservative fixed total alphabet is 4465, independent of n. All payload momentum margins are unchanged by this leading scalar. A complete selected chain with only 47 letters is not claimed.

[assumed, limitation] The charge mode is gapless. Control of finite-coupling mixed operators, fractional marker harmonics, dressing and physical momentum gluing remains the mobile infrared problem. The 23-letter loop family also still lacks an all-n nonidentity proof here. The modified core removes the finite payload obstruction without settling those separate questions.

## Reproduction and handoff

[proved, provenance] model.py, literal.py, ed.py, pencil.py, certify.py and checks.py were copied from C39 into this folder. C42 changed ed.py's large-sector solver tolerance, warm start and off-diagonal diagnostics. budgets.py and leading_checks.py are new. leading_checks.json records an additional C42 full-code g^2 test. This is a rerun and audit of C39's implementation, not an independently authored implementation.

[proved, file map] PROOF.md contains the derivations and parameter inequalities. certificate_n*.json contain rational root brackets and margins; budgets_n*.json contain the finite-parameter witnesses. ed_n*_m*_q*.json retain all numerical errors and residuals. literal_checks.json, bond_checks.json and leading_checks.json retain the local and source checks. PROGRESS.md and STATE.md carry the handoff. All writes were confined to work_C42; the project-level STATE was not modified.

[proved, reproduction commands] Use Python with NumPy and SciPy and one BLAS thread. Each ED process has its own finite time limit. The n=3, three-cell run should be run without competing heavy calculations.

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B literal.py
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B checks.py
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B leading_checks.py
    python3 -B certify.py 1
    python3 -B certify.py 2
    python3 -B certify.py 3
    python3 -B budgets.py 1
    python3 -B budgets.py 2
    python3 -B budgets.py 3
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B ed.py 3 3 1 --seconds 180

[proved, reproduction scope] Repeat the last command for n=1,2,3, cell counts 2,3, and Q=0,1 to regenerate the full ED table. The exact reference certifier and finite-parameter budgets are self-contained within this folder. The general all-n search is specified and proved terminating in PROOF.md, rather than implemented as an unbounded background computation.
