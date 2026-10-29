# C83 result

[Proved construction; conditional infrared conclusion] The twelve-letter
C80 grammar needs two vacancy colors for the diagonal nearest-neighbor
repair preserving all legal vacancy placements. The resulting dimension
is d=14. A newly certified n=13 member meets n=d-1 and supports the
conditional count c=13+1=14>13. No microscopic central charge was measured.

[Proved, restricted] A single uncolored vacancy still hides R_0 V L_0.
Both R_0 V and V L_0 must be allowed in legal embeddings, so a one-state
pair penalty cannot detect the squeezed illegal junction. This excludes
the direct d=13 repair, not every possible d=13 construction.

## Construction

[Proved] V_0 and V_1 carry the preceding material letter's outgoing
grammar state. The incoming rules differ from C82. R_0 and E_0 accept
both states, whereas R_1 and E_1 accept only state 1. L_a,D_s,X_s,M_s
accept state 0. L_a and all markers output state 0; R_a,D_s,X_s output
state 1. V_q accepts and outputs q.

Every illegal squeezed material pair forces a distinct physical mismatch,
including across arbitrarily many vacancies. The legal coloring is unique
at nonzero occupancy. Hopping is

    -t sum_a sum_(q in I(a)) (|a,V_o(a)><V_q,a|+h.c.).

The full bond adds JK times every extended mismatch and J times the
unpenalized C80 occupied bond, with its half-site terms retained inside
occupied contacts. [bond14.py](bond14.py) implements the complete 196 by
196 matrix. [PROOF.md](PROOF.md), Sections 1-3, give every transition and
the coloring argument. Occupied number, marker number and internal U(1)
are exactly conserved. Finite K permits suppressed illegal admixture.

## New n=13 certificate

[Definition] Set D=1000, Delta=4-2cos(pi/15), phi_A=pi/26, phi_B=0,
psi_A=pi, psi_B=pi/2, alpha=sqrt(2/15)sin(pi/15), and

    d_0=1, d_1=1004, d_(j+1)=1004d_j-d_(j-1), beta=1/d_14,
    t_d=beta 2^-74, g=2^-256, lambda=2^-576, K>=2^2304,
    physical field=lambda^2 alpha^4 beta^2 g^2 2^-121.

The b=14 midpoint selector has mu_sel=1.8742920463679156... .

[Proved, computer-assisted] C80's integer/rational pipeline was regenerated
at n=13. The result includes all 351 unordered pair sums with replacement.

| Finite claim | Certified result |
| :-- | :-- |
| Simple Fermi points, exhaustive | 26 |
| Root separation | greater than 0.0149087301 |
| Pair-sum separation | greater than 8.70283e-8 |
| Sea centers | all distinct modulo pi |
| Positive hopping margin in the dyadic metric | greater than 0.5608079 |
| Additional normalized Bloch error allowed | 2^-35 |
| Finite-g error | less than 1.87e-33 |
| Finite-K error including zero-port splitting | less than 1.12e-40 |
| Local two-cell lambda remainder | less than 2.60e-50 |
| All-length selector density gap | greater than 3.8800592e-5 |

The error sum is less than one eighth of the allowed budget. The penalty
estimate divides by the physical lambda^2 band scale, retaining C80's
payload-dependent zero-port splitting. Exact fractions are in
[root_certificate_exact_n13.json](root_certificate_exact_n13.json) and
[budget_exact_n13.json](budget_exact_n13.json). Decimal summaries are in
[root_certificate_n13.json](root_certificate_n13.json) and
[budget_n13.json](budget_n13.json). The root routine's intermediate
finite_g_power=-219 is superseded by the final choice -256.

These statements certify the quadratic projection through second order.
They do not establish Gaussianity of the interacting microscopic model.
The local lambda bound is not a volume-uniform many-cell remainder bound.

## Leakage and transfer of F3 and F4

[Numerical] At J=t, two n=3 zero-bit cells plus one vacancy give the
following illegal material weights. The full reachable closures have
dimension 2304 without color and 4608 with two colors.

| K | Uncolored | Two colors |
| --: | --: | --: |
| 100 | .0267832 | 1.07551e-4 |
| 1000 | .0257054 | 1.08213e-6 |
| 10000 | .0256079 | 1.08279e-8 |

The colored weight scales as K^-2. Literal-bond tests at nonzero port
and conversion, with Q=2, use dimensions 268 and 1045 for one and two
vacancies. Their K=10000 weights are 2.71795e-11 and 1.54444e-11, also
with inverse-square scaling. The largest leakage-test eigenpair residual
is below 1.55e-11. All 38400 nonempty four-site words satisfy the mismatch
inequality. Legal colored and uncolored projected matrices agree entry
by entry in all three tests. See [lift_checks.json](lift_checks.json).

[Proved] F3 transfers to frozen full marker sectors. Unlike C82, a marker
can have a bad pair on its right, but its outgoing state is independent
of which marker it is. That boundary penalty belongs to the right cell.
The only genuine intercell term is the port, with norm sqrt(2)|lambda|.
C80's all-cell inequality and C82's boundary-energy argument then give

    E_wrong-E_selected >= k(gamma-5sqrt(2)|lambda|),
    gamma=delta_selector-1/(K-4)>2^-16.

The schedule makes this positive. A test of complete 1024-dimensional
frozen sectors, including every charge and excitation sector, satisfies
the bound at three resolved port strengths. This proves frozen selection,
not the finite-J mobile dressed-defect comparison.

[Proved] F4 transfers to the commuting integer mismatch penalty, with
kappa=K in the frozen chain and kappa=JK in the lift. At each fixed order,

    U_r H U_r^dagger=kappa N_bad+D_r+R_r,
    [D_r,N_bad]=0,
    ||R_r||_local,exponential <= C_r nu(nu/kappa)^r.

The constants are independent of volume. A 60-dimensional literal-bond
check verifies the second-order sign and K^-2 remainder scaling.
[PROOF.md](PROOF.md), Sections 6-7, supply both derivations.

## Remaining scope and reproduction

[Proved sufficient scale requirement] After fixing positive J and a
material tolerance tau, increase K beyond 2^2304 until JK is large
compared with the unpenalized local strength nu and
C_2 nu^2/(JK)+C'_2 nu^3/(JK)^2<J tau. This is a finite, volume-independent
choice. The frozen lower bound on K alone does not certify the mobile
model for arbitrarily small J. No numerical finite-J matching window
is claimed.

[Assumed] Finite-J vacancy matching, a uniform mobile dressed-defect
comparison, massive-sector matching in the needed norms, stable density
tuning, multichannel RG and the final conformal identification remain
necessary. The two colors repair the hidden-wall mechanism and add no
independent legal color mode. Subject to those inputs, the candidate has
d=14 and c=14. It is not an established microscopic counterexample.

Run from this folder with Python, NumPy and SciPy.

    OPENBLAS_NUM_THREADS=1 python3 -B certify_n13.py
    OPENBLAS_NUM_THREADS=1 python3 -B budget_n13.py
    OPENBLAS_NUM_THREADS=1 python3 -B check_lift.py
    OPENBLAS_NUM_THREADS=1 python3 -B check_frozen.py
    OPENBLAS_NUM_THREADS=1 python3 -B check_sw.py

All new files are confined to work_C83. No network commands, cluster
jobs or subagents were used, and no external files were modified.
[STATE.md](STATE.md) contains the local handoff.
