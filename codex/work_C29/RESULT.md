# C29. A flux obstruction and a two-bus, two-flavor repair

[proved, partial result] A computation-free 26-state repeated bond gives normal, diagonal hopping with two independently adjustable phases and two local chemical potentials at n=2, at leading order lambda^2 g^2. Its quadratic reference has four distinct displaced Fermi points. The phase separation follows from unequal spatial resolvent weights of two virtual bus types, not a flavor label in the physical alphabet.

[proved, limitation] This does not complete the requested arbitrary-n repair. The same rules are nonnormal at n=3, with an exact rational certificate. Normality is proved for the leading conversion coefficient, not for the complete order-lambda-squared bond at fixed nonzero g. The latter has controlled O(g^4) corrections and is numerically non-Gaussian and nonnormal. No completed mobile construction or arbitrary-n d0 is supplied.

[proved, obstruction] Arbitrary position-dependent flux alone cannot give normal exchange with independent phases at n=2 in C22's balanced, equal-rail family. This excludes more than uniform flux, but it is not a universal obstruction to computation-free cells. The two-bus example already escapes its hypotheses.

## 1. What was reused and what changed

[assumed, inherited inputs] C22 supplies the arrow grammar, vacuum walker, frozen selector, exact charge, neutral parity wire and endpoint visibility. C21 and assembly section 7 supply the reason to test the actual exchange matrix rather than count poles. Their target or compiler claims are not assumed solved. All required initial notes and project STATE.md were read.

[proved, provenance] C28 was empty initially. Its subsequently available PROGRESS.md was read. C28 is examining doping and quartic interactions in the original bus. This run instead derives an arbitrary-flux obstruction and adds a second spectrally distinct bus to obtain a normal two-flavor reference. No subagent or independent referee was used.

[proved, microscopic construction] Replace C22's four token species by seven,

    v, A0, A+, A-, B0, B+, B-.

The alphabet has 4 arrow/data letters, 14 token/data letters and 8 marker letters, totaling 26. Each bus has the same local charge conversion as C22. Its neutral hop has flux pi/2 on bulk-to-bulk bonds and zero flux on the bulk-to-marker bond. These are local letter-dependent rules. Set excited diagonal energies relative to the vacuum ground to d_A=3 and d_B=5, hopping magnitude t=1, inter-cell bus amplitudes w_A=1,w_B=8, and inter-cell phases phi_A=pi/2,phi_B=0. The separate n=1 control uses terminal flux pi/3. The full literal bond is implemented in [literal.py](literal.py), with all rules and proofs in [DERIVATION.md](DERIVATION.md).

[proved] A legal cell contains exactly one token and has dimension 7(n+1)2^n. Its unperturbed code is exactly 2^n-dimensional. The first-order active port compressions vanish, and total physical U(1) and parity are exact. The zero-coupling frozen selector is unchanged. The code-gap bound is the minimum of the vacuum-walker gap and min_r(d_r-2 cos(pi/(n+2))-|g|).

## 2. The obstruction and its escape

[proved] For arbitrary link phases in one balanced bus, introduce cumulative kappa and chi phases K_i,X_i as in DERIVATION section 1. The leading hopping is

    T_ij=2i exp(i phi) F_ij(d) exp[i(X_n-X_i+X_j)]
                              sin(K_n-K_i-K_j),

where F_ij is a positive heat-kernel integral. At n=2, both off-diagonal entries contain sin(K_2-K_1), but F_10>F_01. The (0,0) entry of [T,T^dagger] is |T_01|^2-|T_10|^2. Normality forces both entries to vanish. The two surviving diagonal entries then share a phase modulo pi. The proof allows independently prescribed flux on every link and includes a frozen background that only supplies such flux. It does not cover unequal rail dynamics, fluctuating backgrounds, extra bus types or different excitation graphs.

[proved] With terminal flux zero and bulk flux pi/2, each of the two new buses instead contributes a diagonal matrix,

    T^(r)=2i exp(i phi_r) diag(F00(d_r),-F11(d_r)).

Their sum is normal for all choices of the two inter-cell phases. The spatial weight ratios differ,

    F00(3)/F11(3)=265887/845036,
    F00(5)/F11(5)=5885381/19294436.

More explicitly,

    F00(3) F11(5)-F11(3) F00(5)
      =76553167/21709765967228073000 >0.

This proves distinct diagonal hopping phases modulo pi. The two-phase Jacobian is nonsingular at the stated point, so the phases are independently adjustable locally. Independent first-site and last-site density terms provide the two chemical potentials at n=2 without any new letters.

[proved, exact second order] For bus r let H_r be its complete excited-cell matrix and B_r its three pair-source maps, including right-port parity. With alpha_n=sqrt(2/(n+2)) sin(pi/(n+2)), the full order-lambda-squared bond is

    lambda^2 K2(g)
      =-lambda^2 alpha_n^4 sum_r w_r^2 B_r^dagger
                    (H_r tensor I+I tensor H_r)^-1 B_r.

Different bus types and different spatial bonds have zero second-order cross terms. Expanding this formula gives the stated quadratic term at g^2. DERIVATION section 4 proves an explicit Neumann-series O(g^4) bound. It also proves the cubic correction and a finite-cluster analytic remainder bound.

## 3. The effective n=2 chain for numerics

[proved] After subtracting scalar energies, the leading chain is

    H=J sum_(x,a=0,1) [t_a f_(x,a)^dagger f_(x+1,a)+h.c.
                                  +mu_a n_(x,a)],
    J=lambda^2 alpha_2^4 g^2,  alpha_2=1/2,
    mu_0=|t_0|/2,  mu_1=-7|t_1|/10,
    t_0=-2 F00(3)+128 i F00(5),
    t_1= 2 F11(3)-128 i F11(5).

[Numerical evaluation of exact formulas] In these units,

    t_0=-0.0000913087676989 +0.0000986765771459 i,
    t_1= 0.000290195443257  -0.000323497986357 i,
    mu_0= 0.0000672204543571,
    mu_1=-0.000304209677373.

The Fermi points in (-pi,pi] are approximately

    flavor 0:  2.14227574663, -0.493956396677,
    flavor 1: -0.373614468167, 2.05283597813.

The center difference modulo pi is 0.01545108 in absolute value. The velocities divided by J are 0.00026034370025 and 0.00081419521127.

[proved] The dispersions are mu_a+2|t_a| cos(k+arg t_a). Since |mu_a|<2|t_a|, both bands have simple Fermi crossings. Rational alternating-series bounds certify that all four points are distinct, every nontrivial bilinear momentum mismatch exceeds 0.088, and every mismatch between distinct unordered pairs of distinct branches exceeds 0.030. The certificate is [interval_certificate.py](interval_certificate.py). It uses exact fractions, not a grid search.

[proved, software definition] [effective_chain.py](effective_chain.py) exports the leading free chain and the complete order-lambda-squared interacting chain at fixed g. Both use local logical dimension four. The latter includes the higher-body terms of K2(g); it is not substituted by its quadratic fit. Numerical Hermiticity is restored by averaging the mathematically Hermitian bond with its adjoint. Open chains are implemented. Periodic spin-to-fermion conversion additionally requires the parity boundary twist.

[Numerical, four logical cells] The leading 256-dimensional many-body spectrum agrees with all subset sums of its eight one-particle energies within 9.2e-18. This checks the Jordan-Wigner convention and the exported chain, not a central charge.

## 4. Exact diagonalization

[Numerical, literal and source checks] The actual 676-by-676 physical bond is Hermitian and charge-conserving with zero residual. Its restrictions from physical cell dimensions 676 and 17576 agree with the complete legal matrices of dimensions 28 and 84 within 8.0e-15. All port restrictions agree exactly. Independent complete-cell source contractions at n=1,2 agree with the two-bus second-order formula within 1.2e-14. Nonuniform-flux checks at n=1,2,3 agree between the inverse expansion and heat-kernel formula within 1.9e-18. The largest excited-cell matrix in those coefficient checks has dimension 96.

[proved, tested perturbation] The ED Hamiltonian includes the explicitly local term lambda eta E, where E counts excited tokens, with eta=10. It annihilates the code and leaves second order unchanged. It gives the nonzero cubic bond

    K3=2 eta alpha_n^4 sum_r w_r^2 B_r^dagger
                          (H_r tensor I+I tensor H_r)^-2 B_r.

Without that mass probe, the open-chain Hamiltonian is even in lambda and its leading remainder is quartic. This is why a cubic scaling test must specify the probe instead of calling a quartic law a cubic one.

[Numerical, method] [ed.py](ed.py) retains every state in the exact symmetry sectors containing the code, including all vacuum-walker excitations. The number of tokens of each excited bus type modulo two is conserved. The code lies in the even-even sector. No fixed excitation-number truncation is made. Each total-charge block is diagonalized for all its code levels and the next level. QR and Rayleigh-Ritz refinement precede the polar-overlap construction of the canonical effective matrix. Errors below are operator norms including scalar shifts, comparing against the full-in-g second-order formula and the prescribed fields.

| n | Cells | Full legal dimension | Even-even sector | Error at .02 | Error at .01 | Error at .005 | Halving exponents |
| --: | --: | --: | --: | --: | --: | --: | :-- |
| 1 | 2 | 784 | 304 | 7.79048e-5 | 1.02180e-5 | 1.30869e-6 | 2.9306, 2.9649 |
| 1 | 3 | 21952 | 3520 | 1.60308e-4 | 2.07292e-5 | 2.63609e-6 | 2.9511, 2.9752 |
| 2 | 2 | 7056 | 2736 | 1.91849e-5 | 2.53781e-6 | 3.26396e-7 | 2.9183, 2.9589 |
| 2 | 3 | 592704 | 95040 | 3.89153e-5 | 5.11105e-6 | 6.55054e-7 | 2.9286, 2.9639 |

[Numerical] All rows use g=.2. The largest diagonalized charge block has dimension 21924 and the largest code has 64 levels. Subtracting lambda^3 K3 gives halving exponents between 3.9529 and 3.9806. The largest eigenpair residual is 1.09e-12. The longest final calculation took 129 seconds; all cluster calculations finished within minutes. At n=2, three cells and lambda=.005, the next-level gap is 1.41404 and the minimum code-overlap singular value is 0.9999967. Charge, bus-type parity and Hermiticity residuals are zero.

[proved, excluded symmetry sectors] Every omitted odd-type or out-of-code-range charge sector contains at least one excited token. At n=2 its H0 energy is at least 3-sqrt(2)-|g|. A conservative bond norm is sqrt(3)(1+8). At the largest final coupling .02 on three cells, subtracting two such bond norms times lambda still leaves a lower bound greater than .76, apart from the negligible bounded lambda^2 endpoint field. The mass probe is positive. Thus these sectors cannot hide another low code state. The n=1 bound is stronger. Illegal arrow and other marker patterns are not diagonalized here; these are fixed frozen-marker-sector checks, not ED of the entire 26-site-alphabet ring Hilbert space.

[Numerical, resolution limit] The full error is largely scalar. At n=2, three cells and lambda=.005, the traceless residual is 5.34e-12 against eigenpair residual 4.32e-13. It should not be treated as a high-precision interaction measurement. The cubic comparison uses the full error and its analytic coefficient. No finite-size central charge was extracted.

[Numerical, unsuccessful initial scaling protocol] Earlier eta=.37 runs at .04,.02,.01 had cancellations between cubic and quartic scalar terms and did not show a clean cubic halving law. Their JSON records are retained. The final eta=10 protocol resolves that issue and also checks the derived cubic term. An initial interval-certificate assertion requested a bilinear bound of .0892; its conservative enclosure only gave .08919265. The final certified bound is .088. No failed assertion is reported as a passed check.

## 5. What remains open

[Numerical, n=2 finite g] The nonquadratic norm of the exact K2(g), before endpoint visibility, is 2.77748e-7, 1.72190e-8 and 1.07405e-9 at g=.2,.1,.05. The fitted hopping's relative normality defect is 0.00241739, 0.000601167 and 0.000150094. The results are consistent with O(g^4) absolute corrections, and the analytic bound controls them as g goes to zero. They do not make the finite-g bond exactly normal.

[proved, n=3 failure] The same repeated flux rules have T_02=T_20=0 and |T_10|>|T_01|. The positive differences of the necessary spatial weights are

    F_10(3)-F_01(3)=429761/8546653500,
    F_10(5)-F_01(5)=1180141/1487796700500.

Both contributions add in quadrature, so the diagonal commutator entry is nonzero. [exact_weights.py](exact_weights.py) provides a rational 16-by-16 scalar-resolvent certificate. No arbitrary-n normal exchange follows from this construction.

[proved, alphabet accounting] The frozen n=2 demonstration uses d0=26. This is an actual local rule table, not a matrix dimension interpreted as an alphabet. There is no completed mobile d0 and no successful arbitrary-n frozen d0 in this result.

[assumed, conditional arithmetic] The inherited mobile-wrapper count would be 96*26+48=2544 if a compatible 26-letter core were constructed. Its hypotheses have not been established here. In particular, the present adjacency-walker ground state is not supplied as a positive frustration-free local core, adjacent lengths change the data count, and the second-order readout has not been equalized. The number 2544 must not be cited as a realized mobile model.

[heuristic] The detuned n=2 quadratic reference should survive sufficiently weak local corrections under the usual Gaussian-stability input, with exact total U(1) excluding pairing. The rational momentum margins and finite perturbative bounds make that a specified weak-coupling question. This report does not prove the interacting infrared limit or measure c=2.

[proved, final scope] The task's full arbitrary-n mechanism remains open. What is established is (i) a sharper obstruction to every flux-only repair of C22's balanced n=2 bus, (ii) a fixed-alphabet two-bus escape giving independent phases, local fields and normal leading exchange at n=2, and (iii) exact second-order formulas and all four requested finite-cluster ED checks for that partial construction. Neither the restricted obstruction nor the failure at n=3 excludes all computation-free designs.

## 6. Files and reproduction

[proved, file map] [DERIVATION.md](DERIVATION.md) contains the finite proofs. [model.py](model.py) contains the complete legal cell and source resolvents. [literal.py](literal.py) builds the physical bond. [checks.py](checks.py), [exact_weights.py](exact_weights.py) and [interval_certificate.py](interval_certificate.py) provide independent floating-point, rational and interval checks. [ed.py](ed.py) runs cluster ED. [effective_chain.py](effective_chain.py) exports the leading and full second-order n=2 chains. Their JSON records and [ed_summary.json](ed_summary.json) contain the numerical values.

[proved, reproduction commands] Use one BLAS thread. Each ed.py invocation has a 180-second limit.

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 literal.py
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 checks.py
    OPENBLAS_NUM_THREADS=1 python3 exact_weights.py
    OPENBLAS_NUM_THREADS=1 python3 interval_certificate.py
    OPENBLAS_NUM_THREADS=1 python3 effective_chain.py --level leading
    OPENBLAS_NUM_THREADS=1 python3 effective_chain.py --level full
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 ed.py --n 2 --cells 3 --lam .005

[proved, remaining reproduction] For the ED table, run the last command with n=1,2, cells=2,3, and lam=.02,.01,.005. eta=10 and g=.2 are the defaults. No network or cluster computation was used. Only this folder was modified. The project STATE.md was left unchanged under the explicit write restriction; [PROGRESS.md](PROGRESS.md) is the local handoff memory.
