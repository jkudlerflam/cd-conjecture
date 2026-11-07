# C127 checks

2026-11-06. Numerical findings below are floating-point checks, not certified spectral enclosures. Analytic statements and the assumptions on their use are in PROOF.md.

## Reproduction

[Work record] Run from this folder with NumPy and SciPy installed.

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 checks.py
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 audit_checks.py

[Numerical] The first script writes checks.json. The independent audit writes audit.json and asserts the stated inequalities, while recording the expected failures separately. Both finish in about one second each locally. No external project script is imported, no network command is used, and all outputs remain in this folder.

## Exact XX and scalar bounds

[Numerical] The physical XX trace is evaluated by the exact finite-ring Jordan-Wigner parity recursion. Subtracting L times the infinite-line pressure gives log Ztilde. The grid has L=2,4,6,16,64,256 and beta/L=0.01,0.1,1,4,32, totaling 30 samples. All pass the two proved inequalities

    log Ztilde<=log 3+beta[csc(pi/L)-L/pi],
    log Ztilde<=K_XX+B_2(2pi L/(2e beta)),
    K_XX=log 3-2log(1-exp(-1)).

[Numerical] The minimum respective margins are 1.1038483684 and 2.0109920983. The first minimum need not approach zero because log 3 is slack in the bound. At fixed L the beta coefficient, rather than the constant, is sharp. K=0 in the oscillator comparison fails in the high-temperature samples, for example the S=2 product has log-margin -0.0049704811 at L=2,beta=0.2.

[Proved error control used in the numerical check] The independent crossed product keeps N=20000 odd and even frequency levels. Its omitted log-trace tail is bounded by

    4 (beta/(4pi))^L (N-1/2)^(1-L)/(L-1).

This uses asinh u>=log(2u), log(1+z)<=z and an integral tail. For the parity-projected sector, the finite odd-occupation sum is at most its even-occupation sum, so adding the tail increases that sector by no more than the unrestricted fermion tail factor. The code evaluates the bound logarithmically to avoid overflow.

[Numerical] The largest observed crossed-product discrepancy is 0.0051875795 at L=2,beta=64, below the tail bound 0.0051877743. It is a truncation effect. The audit allows an additional 2e-8 absolute tolerance for subtraction of large physical log traces and quadrature error. Every crossed-product comparison lies in its resulting interval. No claim of machine precision is made for the slowly convergent L=2 product.

## Insertion spectrum and decay

[Numerical] The exact XX branch asinh(pi/beta)+asinh((2k-1)pi/beta) is checked at beta=1,10,100,1000 and k=1,2,10,100,1000. The piecewise LR rate passes all samples. The global linear cone fails, as expected and proved in C125.

[Numerical] For the fixed external frequency k=1, 500 distinct positive masses of the reflected density insertion are enumerated at each beta=10,100,1000. Their counts below crossed energy one are 0,8,82. At beta=100, the first eight distinct energies and their positive weights are saved. This demonstrates more than three detected energies in a low-energy window at d=2. The infinite-dimensional cyclic-span conclusion uses the analytic formula, not a rank extrapolation.

[Numerical] Sixteen independent momentum-integral checks of the fermionic resolvent use nu=-3,-0.2,0.2,3 and spatial separations 1,2,3,5. The maximum residual is 5.27e-16. The bubble checks use beta=2,10,100, k=1,2,5 whenever 0<omega<2e, and odd spatial separations 1,3,5,11. They verify the reflected sign and the full LR amplitude with ||Sz||=1/2. The frequency sum uses n from -10000 to 10000. For separation x its omitted positive tail is bounded by

    2 beta^(2x+1)/[2^(2x) pi^(2x+2) (2x+1) 10000^(2x+1)].

[Numerical] The minimum LR margin after adding that tail is 0.0006542961. These checks concern the exact XX model, not a small-Trotter approximation.

## Interacting transfer matrices

[Numerical] The grid uses Delta=0,0.5,1,2, beta=0.25,1,2,8 and m=2,3,4. There are 48 matrices and 240 comparisons at L=2,4,8,16,64. History dimensions are 16,64,256. Energies are defined per spatial site, E=-(1/2)log(u/u0), after diagonalizing the positive two-site transfer U.

[Numerical] The shift by two time entries is resolved before diagonalizing U. At m=4 its q=0,1,2,3 blocks have dimensions 70,60,66,60. The eigenphase is exp(2pi i q/4). At beta=2,m=4 the lowest energies in nonzero sectors are as follows. q=3 agrees with q=1 to roundoff, and q=0 includes the vacuum.

| Model | q=1 minimum | q=2 minimum | log Z_m(4)-log Z_exact(4) |
| :--- | ---: | ---: | ---: |
| XX | 2.4060043730 | 2.9162198218 | 0.0058712717 |
| XXZ Delta=0.5 | 1.8439275876 | 2.2968506401 | 0.0178289457 |
| Heisenberg | 1.5906009528 | 2.0096084123 | 0.0395091329 |
| XXZ Delta=2 | 1.2966068304 | 1.6791028138 | 0.0995250496 |

[Numerical] The bond scale is M=1/2 for Delta=0,0.5,1 and M=3/4 for Delta=2. All scalar samples pass the candidate

    log Ztilde_m(L)<=log 3+B_6(2pi L/(4eM beta)).

| Model | Smallest scalar margin | beta | m | L |
| :--- | ---: | ---: | ---: | ---: |
| XX | 1.0914049080 | 0.25 | 4 | 2 |
| XXZ Delta=0.5 | 0.8134473410 | 8 | 2 | 16 |
| Heisenberg | 0.0350981831 | 8 | 2 | 16 |
| XXZ Delta=2 | 0.0082541890 | 8 | 2 | 64 |

[Numerical] The reference sector comparison also keeps the momentum aliasing. There are three bosonic species of each sign, with every integer frequency reduced modulo m. Positive cyclic convolutions sum the complete geometric occupation series for each oscillator. Their sector sum is independently checked against the scalar product. This avoids Fourier cancellation in exponentially small sectors.

[Numerical failure] The momentum-sector candidate with the same log 3 constant fails. The worst sampled log-margin is -259.5066468 for Heisenberg at beta=0.25,m=4,L=64,q=2. Failure is retained in audit.json. The exact XX lattice argument, rather than an extrapolation of these data, establishes the global linear momentum-sector obstruction.

[Numerical failure] Applying the continuum low-frequency cone directly to the coarse Heisenberg matrix at beta=8,m=2,q=1 gives gap-omega/(4eM)=-0.0724667252. This is not a violation of the continuum theorem. Finite m has a nonzero time step and only aliased momenta. At beta=8,m=4 the four-site log-trace errors are already 0.75537,1.75153,2.85461,4.77674 for Delta=0,0.5,1,2. No low-temperature Trotter extrapolation is claimed.

[Numerical validation] The independently assembled physical four-site Trotter trace agrees with the column contraction to at most 4.03e-12 in log trace. Time-rotation commutators vanish to the test tolerance. The minimum relative Ke eigenvalue is 3.54e-15. The minimum relative U eigenvalue is -1.51e-16, consistent with roundoff. Values of u/u0<=1e-13 are omitted from the listed energies. Their omitted normalized trace is bounded by the omitted count times (1e-13)^(L/2), and is smaller than every passing scalar margin. These are numerical sign checks rather than interval certificates.

## Other analytic comparisons

[Numerical] The right-integral upper bound B_S(a)<=S pi^2/(6a) and the elementary floor-count lower bound are checked at 101 logarithmically spaced a from 0.001 to 100. The charge inequality log theta(u)<=4u is checked at 101 u from 0.001 to 1000. The directions of both Gibbs bounds are checked on a noncommuting two-dimensional example at L=0.1,1,10. All pass.

[Numerical] The nonlinear reference product using the piecewise LR rate is checked with S=6, L=2,4,8,16 and beta=0.1,1,10,100. An explicit power-law tail is included in the upper-bound comparison. All 16 samples pass, with minimum margin 0.0635304259. This tests the conditional reference estimate, not the missing many-insertion comparison to an interacting chain.

[Numerical] The convex midpoint and right-endpoint bounds used for XX are checked at real L=2,3,4,8,16 and beta=0.1,1,10,100. Independent quadrature verifies the integral [pi csc(pi/L)-L]/2. Six-color partition coefficients through level 100 verify the cumulative square-root count in PROOF equation (30). These checks pass.

[Numerical] The diagonal ferromagnetic chain is checked for d=2,3,5 and beta J=4,8,16,24 along even L near exp(beta J/2). At d=5,beta J=24,L=162754, log Ztilde=1.6094133358, approaching log 5=1.6094379124, while the oscillator log trace is below floating-point resolution. This tests the lattice necessity of K>=log d by exact transfer eigenvalues and does not construct a large physical Hilbert space.

## Corrections and unrun work

[Work record] The first product-tail calculation overflowed by separately computing a large power and a small power; it was replaced by a log-domain evaluation. The first reference momentum-sector calculation used a discrete Fourier transform and lost accuracy in tiny sectors; it was replaced by positive cyclic convolution before the reported sector results and audit. Neither issue affected the analytic arguments. The complete local scripts were rerun after these changes.

[Work record] run_typhon.slurm prepares m=5,6 and the same model/temperature grid. checks.py --cluster requires a SLURM job identifier. This larger calculation was not run locally and was not submitted, because Typhon is unavailable in this sandbox. The scripts use CPU linear algebra on the requested apollo partition; no GPU implementation is claimed.

[Assumed or unresolved] No numerical test establishes a universal dimension-dependent species count, a uniform interacting zero-mode comparison, or a scalar majorant beyond the sampled matrices. The general scalar question remains open.
