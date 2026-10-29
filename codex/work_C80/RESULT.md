# C80 result

[Proved finite ingredients; conditional infrared conclusion] Twelve letters
suffice for a finite n=12 member with 24 simple Fermi points, distinct sea
centers, positive margins for all 300 unordered pair sums, and exact U(1).
The cell length is 13. This improves C75's fifteen-letter construction.

The frozen-chain prediction is c=12>11. Adding one vacancy gives a
thirteen-letter candidate with predicted c=13>12. Both predictions retain
the parent construction's infrared and dressed-defect assumptions. No
central charge was measured, and no completed microscopic counterexample
is claimed.

## Construction

[Proved] Use bitless tokens and shift the displaced bit to the next site.
The final bit is stored on the empty marker. Conversion flips the bit
immediately to the token's right, using a two-site term. The folded
carrier graph and the constrained second-order operator are preserved.

This first gives fourteen letters. Identifying the vacuum marker with
the zero-bit empty marker gives thirteen. Identifying the vacuum token
with R_0 gives twelve. The resulting alphabet is

    L_0,L_1,R_0,R_1,D_A,D_B,X_A,X_B,E_0,E_1,M_A,M_B.

The full bond, grammar and charge table are in [PROOF.md](PROOF.md),
Section 1, and [bond_shift.py](bond_shift.py). The bond is 144 by 144,
with 205 nonzero entries at a generic diagnostic point. Its Hermiticity
and charge commutator pass the literal matrix checks.

[Proved limitation] At finite penalty the new vacuum code has a small
payload-dependent splitting. It is bounded relative to the physical
lambda^2 band scale, including the zero-port contribution. Treating this
code as exactly degenerate would invalidate the twelve-letter argument.

## Certified member

[Definition] Set n=12, D=1000, Delta=4-2cos(pi/14), and

    phi_A=pi/24, phi_B=0, psi_A=pi, psi_B=pi/2.

Let d_0=1, d_1=1004, d_(j+1)=1004d_j-d_(j-1), beta=1/d_13,
and alpha=sqrt(2/14)sin(pi/14). A sufficient finite schedule is

    t_d=beta 2^-66,
    g=2^-216,
    K=2^2048,
    lambda=2^-480,
    physical field=lambda^2 alpha^4 beta^2 g^2 2^-112.

The marker chemical potential is the exact midpoint selector in PROOF
(15), numerically 1.8561806413628503. The hierarchy remains severe.

[Proved, computer-assisted] The regenerated n=12 certificate gives

| Quantity | Certified bound |
| :-- | :-- |
| Simple Fermi points | 24 |
| Minimum Fermi-point separation | greater than 0.01509618 |
| Minimum distinct unordered pair-sum separation | greater than 5.5894e-7 |
| Pair sums, including repeated points | 300 |
| Positive Hermitian hopping margin after congruence | greater than 0.5652 |
| Additional normalized Bloch error allowed | 2^-35 |
| Finite-g error | less than 4.93e-18 |
| Finite-K error, including zero-port splitting | less than 1.16e-55 |

The full scalar loop uses 512-bit fixed-point integer enclosures.
Endpoint inertias and pair-phase comparisons use rational arithmetic.
`root_certificate_exact_n12.json` retains exact brackets and margins.
`finite_budget_exact.json` retains the rational transport bounds. These
certify the quadratic projection through second order, not an exactly
Gaussian interacting Hamiltonian or independence of all momentum relations.

[Proved] The zero-port selector survives finite K uniformly in cell
length, including illegal sectors and marker-free rings. For k data sites,

    mu_sel-2cos(pi/(k+2))-max(k-1,0)/(K-4)
      <= E_K(k) <= mu_sel-2cos(pi/(k+2)).

The nearest competing energy-density gap exceeds 0.0000547433.

## Checks and rejected alternatives

[Numerical] Legal-space H_0 and V agree entry for entry with C75's
18-letter parent at n=1,2. The full g^2 pair operator has relative
non-Gaussian residual below 3.5e-16. Exact word enumeration through n=6
checks completeness, injectivity, bit count and charge.

[Numerical] Full invariant-space ED, including illegal states, uses
346 and 8224 states at n=1,2. All 4 or 16 code levels agree with the
quasidegenerate second-order formula with error proportional to lambda^4.
At n=2, error/lambda^4 is 0.00174381, 0.00173938, 0.00173829 for
lambda=0.2,0.1,0.05. Eigenpair residuals are below 1.6e-13.
The zero-port pair-band splitting decreases from 0.00495170 at K=100
to 0.000499514 at K=1000. These checks use resolved couplings.

[Numerical] C75's direct equal-charge data and express merges introduce
non-Gaussian terms already at g=0. A different fold using a fast direct
B carrier loses onsite cancellation and fails the common-field crossing
test at n=2,3 for five speeds. The literal no-fold deletion also fails
at n=2. [Proved, restricted] A two-token order encoding cannot implement
the required charge-changing conversion. Putting bits only on a fixed
number of tokens cannot retain the exponentially large payload code.
The successful design instead uses the marker for displaced-bit storage
and as one port. PROOF Section 8 gives the mechanisms and test scope.

[Heuristic] Twelve is the smallest alphabet supported by a complete
finite construction here. The twelve possible charge-preserving
single-letter identifications were audited. They collapse legal states,
introduce order-one vacuum splitting, or destroy the conserved marker
positions. None supplies an eleven-letter continuation of this proof.
There is no global minimality theorem. A redesigned eleven-letter model
is not excluded.

## Reproduction and remaining work

Run from this folder, with Python, NumPy and SciPy.

    OPENBLAS_NUM_THREADS=1 python3 run_certificate.py
    OPENBLAS_NUM_THREADS=1 python3 checks.py
    OPENBLAS_NUM_THREADS=1 python3 selector_checks.py
    python3 grammar_audit.py
    python3 alternatives.py
    python3 merge_c75.py
    python3 fast_carrier.py

[Open] A common-bulk-subtracted defect comparison at nonzero port and
volume-uniform infrared matching remain necessary. The vacancy lift
also needs its existing matching argument. No subagent, external referee,
network command or cluster job was used. No file outside work_C80 was
modified. [STATE.md](STATE.md) carries the local handoff.
