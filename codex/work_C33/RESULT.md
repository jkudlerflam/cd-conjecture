# C33 result

[Proved, computer-assisted] The 23-letter chain passes the requested finite band test at n=23. An explicit phase choice and finite parameters give a positive Hermitian hopping matrix, 46 simple Fermi points, distinct sea centers, and

    delta_2 > 0.0140193,
    delta_4 > 3.1558e-7.

The quartic test includes all distinct unordered pair sums, with repeated Fermi points allowed. These bounds hold for the quadratic projection of the exact second-order interaction at finite g, not only for its leading g^2 coefficient. The proof combines rational matrix-entry enclosures, exact integer conditioning, rational endpoint inertia tests, and exact pair-phase comparisons. See [PROOF.md](PROOF.md), especially Section 9.

[Proved] For every n, a separate explicit schedule D_x(n), g(n), epsilon(n) guarantees positive Hermitian hopping and common-field crossings. The proof replaces the informal Gram-limit argument by an exact positive convolution kernel and a finite-D estimate.

[Open] I do not prove distinct centers or Sidon-compatible margins for every n. No failure of the adjustable-parameter family was found. Thus the all-n question remains open, while the finite n near 23 needed for the proposed d=23 application has been reached. No alphabet enlargement is needed for this result.

## Explicit certified n=23 point

[Definition] Let alpha=sqrt(2/25) sin(pi/25), and define

    d_0=1, d_1=1004, d_(j+1)=1004 d_j-d_(j-1), beta=1/d_24.

Use the original nearest-neighbor 23-letter bond with

    t=t_x=t_switch=1,
    D_x=D_M=1000,
    Delta=4-2 cos(pi/25), delta=0,
    g_+=g_-=2^(-365),
    t_d=beta 2^(-136),
    phi_x,0=pi/46, phi_x,+ = phi_x,- = 0,
    phi_d,0=pi, phi_d,+ = phi_d,- = pi/2,
    phi_port=pi/2.

Stub and switch phases are zero. The common physical data field is

    lambda^2 epsilon,
    epsilon=(2 alpha^4 beta^2 g^2) 2^(-214).

[Proved] The new loop phase choice makes the three scalar carrier graphs identical after gauging the tree. It is a parameter change within the supplied alphabet. The unperturbed code, U(1), parity and zero first order are preserved. The certificate concerns the frozen legal-cell sector; the inherited marker selector is not re-audited here.

[Proved, computer-assisted] In the conditioned metric, the Hermitian hopping margin exceeds 0.5028. The certificate permits an additional normalized Bloch perturbation of norm 2^(-33). An exact rational inequality shows that the entire higher-conversion correction to the quadratic Bloch matrix at g=2^(-365) fits this budget. The same root brackets consequently apply at finite g. The one-particle hopping is positive as well.

[Numerical display of proved formulas and bounds] All logarithms below are base ten, with lambda=1 referring to the coefficient of lambda^2, not a finite-lambda approximation claim.

| Quantity | log10 value or bound | Status |
| :-- | --: | :-- |
| beta | -72.0416 | Exact recurrence, decimal display |
| t_d | -112.9817 | Explicit parameter |
| g | -109.8759 | Explicit parameter |
| Largest entry of leading T | -371.6909 | Numerical value from certified entries |
| Smallest Hermitian eigenvalue, finite g | at least -433.2504 | Proved lower bound, logarithmic display |
| Common field epsilon | -433.7560 | Explicit parameter |
| Smallest leading Fermi speed | -431.6922 | Numerical, 46 roots |
| Norm of all g^4-and-higher conversion terms | at most -445.9074 | Proved operator-norm bound |

[Proved scope of the scale statement] The last row bounds the entire correction to the leading Gaussian interaction, including its quadratic corrections. It is not a separately measured non-Gaussian projection norm. All logical energies and speeds acquire an additional factor lambda^2. The absolute scales fall below float64 range, so the saved integer matrices and logarithms, rather than a float64 microscopic Hamiltonian at these parameters, are the reproducible representation.

## What holds for every n

[Proved] Set the loop to zero and let Gamma be the isolated-stub filter Gram matrix. The exact normalized finite-express matrix obeys

    ||F/beta^2-Gamma|| <= n(n+2)/D.

The mechanism is explicit. The right filter is the left filter convolved with a positive normalized kernel whose mean delay is at most 2(n+1)/(D+2). Its Laplace transform is the product of 2(n+1) exponential-kernel transforms.

[Proved] Equations (8)-(14) of PROOF.md give a positive rational lower bound L_n for Gamma and the explicit schedule

    D(n)=ceil[4n(n+2)/L_n],
    s_n=L_n/2, b_n=[D(n)+6]^(-(n+1)),
    g(n)=b_n s_n/[16(n+1)],
    epsilon(n)=alpha_n^4 g(n)^2 b_n^2 s_n/2.

Both the exact one-particle hopping and the exact quadratic projection have positive Hermitian part. Every ordered quadratic band crosses the common field on both semicircles. Multiplicity and resonance are not resolved by this theorem alone.

[Numerical display of a proved schedule] The analytic bound is intentionally loose. At n=23 it gives log10 D about 8403.34 and log10 g about -210083.13. The rational certificates make the much smaller D=1000 usable at the target size. These numbers are existence bounds, not optimized physical scales.

[Proved obstruction on the tree] At t_d=0 the centers coincide modulo pi and interband Cooper transfer has zero mismatch. Shrinking g cannot cure this resonance. The loop is necessary for the certified displaced target above. The missing all-n statement is that a loop choice breaks every relevant pair-sum identity, with quantitative margins, not that the Gram matrix has full rank.

## Size dependence and computation limits

[Proved, computer-assisted] At a=4, D=1000, the following are rational lower and upper bounds on the smallest eigenvalue of Herm(F/beta^2). The last column is the numerical leading hopping entry scale before multiplying by g^2. Every positivity decision in the first two columns uses integers.

| n | log10 lower bound | log10 upper bound | log10 max entry of T/g^2 |
| --: | --: | --: | --: |
| 2 | -5.198 | -4.897 | -21.264 |
| 4 | -10.766 | -10.164 | -34.224 |
| 8 | -21.895 | -20.992 | -59.518 |
| 12 | -33.022 | -31.943 | -84.396 |
| 16 | -44.148 | -42.943 | -109.060 |
| 23 | -63.615 | -62.253 | -151.939 |
| 32 | -88.644 | -87.139 | -206.774 |

[Numerical, n through 64] The positive scalar-resolvent method evaluates leading matrix entries through n=64 in seconds, without a 2^n space. At D=1000, log10 max|T/g^2| is -303.845 at n=48 and -400.639 at n=64. Double-precision eigenvalue signs are unresolved at these sizes and are not claimed. Integer positivity certificates reached n=32. The n=48 certificate stopped at its 175-second limit.

[Numerical, exact-in-g inverse] The polynomial fixed-charge calculation reaches n=12 at D=10, a=4, g=.05, using 641875 intermediate states. It took 32 seconds in the initial run. At n=8, D=1000, g=.05, the exact-in-g hopping has max|T|=7.58286e-63 and a positive numerical Hermitian margin 2.85625e-20 relative to that entry scale. These calculations invert the finite-g operator, with a bounded Neumann-series tail. Their floating-point roundoff is not interval-certified. Attempts at n=16,D=10 and n=12,D=1000 hit their per-run limits and are not counted as completed results.

[Numerical, failure at fixed parameters] For a=4,D=10,g=.05 the Hermitian part is already indefinite at n=3. This is a failure of that fixed express energy, not of the increasing-D construction. At g^2 order, even D=100 is insufficient by n=12. The delay kernel can rotate the narrow filter combinations enough to change their sign.

[Numerical, loop scans] Separate full finite-loop g^2 calculations give positive quartic margins through n=10 at a=2.2,D=100. Those scans use direct imaginary onsite expressions and conditioned root pencils. They are supporting numerical results; the stronger n=23 statement comes from the rational certificate above.

## Checks and remaining scope

[Numerical checks] The exact-in-g reduction agrees with literal ycell.py at n=1,2,3 at a resolved test point, with relative errors below 7.1e-10. The independent source-block calculation and integer finite-loop calculation agree at n=2,3 to relative errors below 1.1e-15. The tree limit of the general leading calculation agrees with the positive scalar formula. At n=2, the literal common-loop quadratic fit has relative hopping errors 0.0037610, 0.00093848, 0.00023451 at g=.1,.05,.025, confirming the expected g^2 relative scaling.

[Rejected checks, retained] The original eigensystem extraction loses relative accuracy when T is tiny. The first 256-bit full-loop attempt also failed its own error enclosure after narrow-band conditioning. That failed output is preserved as full_loop_256bit_rejected.json. The final 512-bit calculation has conditioned entry error below 1.697e-50 and passes every rational root test. No conclusion uses the rejected eigenvalues.

[Open] The all-n Sidon theorem, frozen-cell length-sector issues, mobility and the interacting infrared limit remain separate. No central charge was measured, and this report does not claim a completed counterexample to c<=d-1. The finite n=23 Gaussian band conditions are now proved for an explicit member of the original alphabet.

## Reproduction and files

[Reproduction] Python, NumPy and SciPy are the only required packages. Python's standard-library integer and Fraction arithmetic supplies the certificates. From this folder, use one BLAS thread. Each listed calculation is bounded to minutes.

    OPENBLAS_NUM_THREADS=1 python3 certify_tree.py 23
    OPENBLAS_NUM_THREADS=1 python3 rational_bands.py 23
    OPENBLAS_NUM_THREADS=1 python3 full_loop_integer.py
    OPENBLAS_NUM_THREADS=1 python3 certify_roots.py
    OPENBLAS_NUM_THREADS=1 python3 summarize.py

[File map] PROOF.md gives the derivations. structured.py and fast_loop.py implement the polynomial leading reductions. finite_g.py implements the exact-in-g one-particle inverse. certify_tree.py gives rational tree enclosures. rational_bands.py supplies fraction-free conditioning. full_loop_integer.py encloses the finite-loop matrices. certify_roots.py proves endpoint inertias, simplicity, momentum margins and the finite-g budget. root_certificate_n23.json summarizes the successful certificate; root_certificate_exact_n23.json retains rational brackets and margins. schedule.py evaluates the all-n analytic schedule through logarithms. PROGRESS.md and the local STATE.md are handoff records. No file outside this folder was modified.
