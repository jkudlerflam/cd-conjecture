# C37. Certified distinct centers and quartic nonresonance through n=23

[Proved, computer-assisted finite-dimensional result] The original 23-letter candidate has a parameter choice with exactly 2n simple Fermi points, distinct sea centers modulo pi, and no nontrivial number-conserving quartic resonance, for every n=1,...,23. The test includes pair transfer and umklapp modulo 2pi. No letters or letter rules were added.

[Proved] The result uses the existing loop with neutral phase pi, a nonzero common field, C33's positive hopping certificates, and a new exact finite-field certificate for every unordered pair of Fermi roots. It avoids resolving the small Fermi-point differences numerically. The largest carrier-pair inverse has dimension 2209. Its pair-product polynomial has degree 1081.

[Proved] Explicit margins, velocities and a conversion schedule are given below. They control the full finite-g, second-order quadratic projection. The complete second-order non-Gaussian remainder is bounded as well. Higher-order thermodynamic transfer remains conditional on the already identified local-elimination and dressing input. The result is not a proof of an interacting central charge.

[Open] There is no all-n nonresonance theorem here. No claim is made for n>23 or for C32's different mobile core.

## The concrete parameter choice

[Proved construction] In the notation of explicit_candidate/README.md, use

    t = t_x = t_sw = 1,
    D_x = D_M = 1000,
    Delta_n = 4 - 2 cos(pi/(n+2)),
    delta = 0,       g_+ = g_- = g_n,
    phi^x_0 = pi/(2n),   phi^x_+ = phi^x_- = 0,
    phi^d_0 = pi,        phi^d_+ = phi^d_- = 0,
    phi_port = pi/2,
    t_d = 10^(-2000).

All stub and switch phases vanish. These coefficients use the same 23-letter bond. D=1000 is a common sufficiently large choice on the finite certified range, not a proposed fixed-D all-n schedule. Delta_n>=2 exceeds the chosen conversion amplitude.

[Proved convention] Put

    a4_n = [2/(n+2)]^2 sin(pi/(n+2))^4,
    W_n = lambda^2 g_n^2 a4_n,
    epsilon_physical = W_n * 10^(-1000).

After removing W_n, the leading symbol is h(k)=T exp(ik)+T^dagger exp(-ik)+M+10^(-1000) I. Its coefficients are in Q(i). The physical field is the uniform data field in the candidate, not a separately adjustable field for each flavor. If the conversion is reduced in the conditional schedule below, replace g_n by that chosen g in both W_n and epsilon_physical.

## Why the certificate answers the resonance question

[Proved] Let z_1,...,z_(2n) be the roots of

    f(z)=det[T z^2+(M+10^(-1000) I)z+T^dagger].

The C33 hopping certificates, with an elementary norm bound for the small loop and field, force 2n crossings on the unit circle. The new calculation proves that

    Q_pair(w)=product_(i<=j)(w-z_i z_j)

has no repeated root. Thus all original roots are simple and every two different unordered pairs have different momentum sums modulo 2pi. Every ordered band has exactly two roots. Comparing the two endpoints of each band gives distinct sea centers modulo pi.

[Proved] The calculation reduces Q(i) modulo p=65521, with i mapped to 24297. It verifies gcd(Q_pair,Q_pair')=1. Nonzero reduction proves that the characteristic-zero discriminant is nonzero. This is an exact algebraic certificate, not evidence from a finite momentum grid. Repeated indices are included in the pair list, which is a conservative quartic test. Identical-pair density/current interactions remain allowed. Exact U(1) forbids number-changing resonances.

[Proved, attribution] Positivity was not rederived. The existing C33 certificates at n=2,4,8,12,16,23 were reused, and its certifier was evaluated only at the missing intermediate sizes. An exact integer comparison gives

    lambda_min[-(T_tree+T_tree^T)/2] > 10^(-250)

for every n<=23, in the normalization without a4_n. The field and loop move the symbol by less than half this bound.

[Proved, certificate coverage] In the table, the gcd degree is an exact finite-field result. The last column is only a rounded display of the stronger inherited integer lower bound for the normalized tree hopping.

| n | Pair inverse dimension | Pair-polynomial degree | gcd degree | log10 hopping lower bound |
| :-- | --: | --: | --: | --: |
| 1 | 9 | 3 | 0 | -14.114 |
| 2 | 25 | 10 | 0 | -22.907 |
| 3 | 49 | 21 | 0 | -31.695 |
| 4 | 81 | 36 | 0 | -40.482 |
| 5 | 121 | 55 | 0 | -49.268 |
| 6 | 169 | 78 | 0 | -58.054 |
| 7 | 225 | 105 | 0 | -66.840 |
| 8 | 289 | 136 | 0 | -75.626 |
| 9 | 361 | 171 | 0 | -84.411 |
| 10 | 441 | 210 | 0 | -93.196 |
| 11 | 529 | 253 | 0 | -101.981 |
| 12 | 625 | 300 | 0 | -110.766 |
| 13 | 729 | 351 | 0 | -119.551 |
| 14 | 841 | 406 | 0 | -128.336 |
| 15 | 961 | 465 | 0 | -137.121 |
| 16 | 1089 | 528 | 0 | -145.906 |
| 17 | 1225 | 595 | 0 | -154.690 |
| 18 | 1369 | 666 | 0 | -163.475 |
| 19 | 1521 | 741 | 0 | -172.259 |
| 20 | 1681 | 820 | 0 | -181.044 |
| 21 | 1849 | 903 | 0 | -189.828 |
| 22 | 2025 | 990 | 0 | -198.613 |
| 23 | 2209 | 1081 | 0 | -207.397 |

## Explicit margins and corrections

[Proved] The following intentionally conservative bounds avoid any root-finding assumption. Define

    B_n = 32n,
    Qbar_n = 10^1000 (2014*10^2000)^(9(2n+1)^2),
    H_n = (B_n Qbar_n)^n,
    N_n = n(2n+1),
    delta_n = (2 H_n)^(-N_n(N_n-1)),
    v_n = Qbar_n^(-n) delta_n^(2n-1)/B_n^(n-1),
    r_n = delta_n/64,
    tau_n = Qbar_n^(-n) r_n^(2n)
            /[2^(n+4) n (B_n+1)^(n-1)],
    g_n^2 = tau_n/(256n).

[Proved] Qbar_n bounds a common denominator of T,M and the field. The nonzero Gaussian-integer discriminant after clearing denominators gives

    delta_2 >= delta_n,        delta_4 >= delta_n,
    center separation modulo pi >= delta_n/2,
    normalized minimum speed >= v_n.

The center comparison is empty at n=1. The physical minimum speed is at least W_n v_n. PROOFS.md gives the denominator bound, discriminant argument and velocity derivation step by step.

[Proved] An analytic Hermitian Laurent-symbol perturbation of norm at most tau_n on 1/2<|z|<2 preserves exactly 2n simple unit-circle roots. Each root moves by less than delta_n/32 in angular distance. In particular the corrected quadratic and quartic margins exceed delta_n/2, and the center separation exceeds delta_n/4. A minimum-modulus argument also gives the corrected normalized speed bound

    vtilde_n = Qbar_n^(-n) delta_n^(2n-1)
               /[2^(3n+1)(B_n+1)^(n-1)].

[Proved] The complete second-order bond beyond its g^2 coefficient has norm at most a4_n g^4. The chosen g_n puts its normalized quadratic correction below tau_n/8. Its normalized non-Gaussian norm is at most (2+16n^2)g_n^2, and obeys

    epsilon_NG/(v_n delta_n) < 10^(-6),
    (epsilon_NG/v_n) log(1/delta_n) < 10^(-6).

These ratios use the leading reference speed. With the corrected speed bound vtilde_n and margin delta_n/2, both ratios are less than 10^(-4). These are lattice coefficient bounds. Converting them to dressed roots and continuum couplings still requires the stated infrared input and its matching constants. The conversion can be decreased further without changing the leading certificate.

[Assumed, thermodynamic inputs] Suppose local elimination bounds the additional bare quadratic symbol by C_n lambda^2/(a4_n g^2), and controlled dressing bounds the retained interaction's quadratic self-energy by A_n(2+16n^2)g^2, in the same analytic norm. A_n,C_n are finite and independent of ring length. The latter bound includes matching constants and velocity denominators. The interaction dressing does not disappear by decreasing lambda alone.

[Proved conditional schedule] The sufficient choices are

    g^2 <= min(g_n^2, tau_n/[4 A_n(2+16n^2)]),
    lambda^2 <= a4_n g^2 tau_n/(4 C_n).

They put the finite-g quadratic error, interaction dressing and higher-order bare error below tau_n/8, tau_n/4 and tau_n/4. The margins therefore survive. Decrease lambda further for the inherited selector/exclusion hierarchy. Finite-cluster analyticity proves a bare higher-order bound with a cluster-dependent constant. It does not prove either thermodynamic assumption. The order-lambda^4 effects of placing the physical data field inside virtual denominators are part of the bare error budget.

[Numerical evaluation of proved bounds] These values are logarithms of the explicit lower bounds and chosen conversion, not estimates of the actual margins. The exact formulas above specify the witness.

| n | log10 delta_n | log10 v_n | log10 g_n |
| :-- | --: | --: | --: |
| 1 | -9.796166e+05 | -1.142884e+06 | -1.061254e+06 |
| 2 | -8.131417e+07 | -2.448460e+08 | -1.630801e+08 |
| 3 | -1.114419e+09 | -5.574746e+09 | -3.344582e+09 |
| 4 | -7.365511e+09 | -5.156442e+10 | -2.946497e+10 |
| 8 | -7.654793e+11 | -1.148223e+13 | -6.123855e+12 |
| 12 | -1.213058e+13 | -2.790036e+14 | -1.455671e+14 |
| 16 | -8.741862e+13 | -2.709978e+15 | -1.398698e+15 |
| 23 | -1.069482e+15 | -4.812667e+16 | -2.459808e+16 |

## Band scales

[Proved] Let beta_n be the reciprocal determinant of the n+1 vertex express block with diagonal 1004 and off-diagonal -1. Then

    1006^(-n-1) <= beta_n <= 1002^(-n-1),
    J_n = W_n ||T||
      <= W_n [n beta_n^2/16 + 32n*10^(-2000)].

The crossing excursion is at least (3/2) W_n*10^(-250). The normalized velocity lower bound is v_n. The target-cell gap is at least the minimum of 2-|g_n| and the vacuum-walker gap 2[cos(pi/(n+2))-cos(2pi/(n+2))]. The latter is of order n^(-2), while a4_n is of order n^(-6).

[Numerical display] At n=23 the tree hopping upper bound before multiplying by lambda^2 g_n^2 is about 10^(-149.73). The C33 sharper tree crossing lower bound, in the same units including a4_n, is about 10^(-213.20). The much smaller conversion schedule above is chosen to make the algebraic separation bound survive, not to optimize an observable energy window.

[Proved limitation] These results close the requested finite-dimensional separation and nonresonance gap through n=23. They do not provide a central-charge measurement, prove local elimination, audit the unrestricted selector, or establish an all-n statement. A successful modular certificate at any additional n can be lifted to larger D and smaller field/loop while preserving its residues, so it is compatible with the required increasing-D hierarchy. Proving that such certificates exist for every n remains open.

## Checks and reproduction

[Proved, exact checks] audit.py verifies primality and the image of i, compares sparse and dense modular inverses through dimension 121, compares the symmetric-square polynomial against direct pair multiplication through degree 1081, and compares independent coefficient contractions at n=1,2,3. All exact checks passed. The main certificate checks each modular solve by its residual and checks its determinant polynomial at extra points.

[Numerical, literal-cell check] literal_coefficient.py directly extracts the g^2 coefficient of the literal candidate at n=1,2,3, without subtracting two finite-g Hamiltonians. Pair-space dimensions are 324,3600,28224. Hopping relative discrepancies are at most 1.4e-12, and onsite discrepancies divided by the hopping scale at most 2.3e-10. The earlier finite-g subtraction lost precision in the onsite term, reaching 4.5e-5 at n=3. It is retained in literal_results.json and is not used as a certificate.

[Proved/numerical controls] The finite-field test detects the expected collisions on removing the loop at n=2,3,4 and on removing the field at n=1,2,3,4. Thus neither positivity alone nor the loop without a field is silently substituted for the full test.

[Proved file map] PROOFS.md contains the full derivations. modular.py builds the leading matrices over the residue field. certify_pairs.py creates pairs_n*.json. positivity_C33.py is the copied C33 certifier, complete_positivity.py reuses its prior results and fills missing sizes, and verify_bounds.py performs the exact uniform-bound comparisons. crossing_bounds.json and positivity_summary.json record the positivity evidence. scales.py displays the proved bounds using logarithms. The copied structured.py and loop_scan.py are C33 code used only for diagnostic floating-point checks.

[Proved reproduction] From this folder, with Python 3, NumPy and SciPy only,

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B audit.py
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B certify_pairs.py 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B complete_positivity.py
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B verify_bounds.py
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B literal_coefficient.py
    python3 -B scales.py

[Proved resource record] All computations completed within minutes. The original n=23 modular calculation took about 12 seconds; the final rerun using only integer modular arithmetic took about 33 seconds. The main algorithm uses carrier matrices, not exponentially large logical spaces. Only the explicitly labeled n<=3 literal cross-check uses those small logical spaces. No network, subagents, external-file edits, commit or push was used. Project STATE.md lies outside the authorized write scope; this folder's PROGRESS.md and STATE.md carry the handoff instead.
