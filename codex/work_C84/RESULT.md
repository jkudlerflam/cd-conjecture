# C84 result

[Proved finite ingredients; assumed infrared transfer] Eleven letters
suffice, with a certified frozen member at n=11 flavors and cell length
13. This improves C80's twelve-letter construction. The conditional
infrared prediction is c=11>10. No central charge was measured, and the
many-cell infrared and nonzero-port defect arguments remain open.

[Proved] Remove C80's E_1 marker and fix the last displaced bit to zero.
A cell with p data sites retains p-1 payload bits. Taking p=12 gives
eleven flavors. The alphabet is

    L_0,L_1,R_0,R_1,D_A,D_B,X_A,X_B,E,M_A,M_B.

The folded carrier and its last data vertex remain, but that last vertex
does not convert. A penalty on L_1 followed by a marker enforces the
fixed boundary bit. The exact legal graph has (5p+3)2^(p-1) states.
This loses one flavor per cell, not one flavor per unit of alphabet.
The full bond and all derivations are in [PROOF.md](PROOF.md).

[Proved, computer-assisted] The regenerated C80 pipeline gives

| Quantity | Certified result |
| :-- | :-- |
| Alphabet and active flavors | a=n=11 |
| Physical cell length | 13 |
| Simple Fermi points | 22 |
| Minimum root separation | greater than 0.0164202883 |
| Minimum distinct unordered pair-sum separation | greater than 9.11916e-7 |
| Unordered pairs, including repeats | 253 |
| Hermitian hopping margin in the dyadic metric | greater than 0.56526 |
| Additional normalized error allowed | 2^-35 |
| Finite-g error | less than 3.380e-21 |
| Finite-K error, including zero-port splitting | less than 7.946e-59 |

All eleven sea centers are distinct modulo pi. The certificate uses
512-bit integer enclosures, rational endpoint inertias and exact pair
phase comparisons. The bounds cover the quadratic projection through
second order in the port coupling, including finite g and finite K.
They do not assert that the effective interacting Hamiltonian is Gaussian.

[Definition] The frozen member uses D=1000, Delta=4-2cos(pi/14),
phi_A=pi/24, phi_B=0, psi_A=pi, psi_B=pi/2. Set

    d_0=1, d_1=1004, d_(j+1)=1004d_j-d_(j-1), beta=1/d_13,
    alpha=sqrt(2/14)sin(pi/14),
    t_d=beta 2^-55, g=2^-216, K=2^2048, lambda=2^-480,
    physical field=lambda^2 alpha^4 beta^2 g^2 2^-101.

The exact midpoint marker selector is in PROOF (18). The severe coupling
hierarchy remains. [MEMBER.json](MEMBER.json) records the member using
integer powers and exact defining expressions.

[Proved] The zero-port all-length selector survives finite K. The new
Hamiltonian is a compression of C80 plus a positive boundary penalty.
Its length-13 density gap exceeds 5.4743300e-5. Illegal sectors and
marker-free rings are included. The selected finite-K payload band has
a small splitting, which is budgeted in physical lambda^2 E_unit units.

[Numerical] The literal bond has 166 nonzero entries and passes Hermiticity,
U(1), and frozen-marker checks. Independent legal-space matrices agree
entry for entry at n=1,2. Full invariant-space ED, including illegal
states, uses dimensions 2221 and 42424. All 4 or 16 code levels exhibit
fourth-order port errors, with eigenpair residuals below 1.6e-13. Grammar
and charge checks pass through seven data sites. The first attempted
loop exponent 66 failed the rational pair-bracket test; exponent 55
passes. The failed parameter is retained in `rejected_parameters.json`.

[Proved, restricted] The faithful shifted-bit, two-rail encoding requires
at least eleven letters: four binary L/R backgrounds, four charged rail
states, and three marker states. These distinctions are justified in
PROOF Section 8, including the general additive-charge equations. Keeping
the on-site charge roles and the retained local port equations gives the weaker
bound a>=7. Neither is a lower bound for all redesigned route-A models.
All ten literal charge-preserving single identifications fail the inherited
construction. [ALPHABET.md](ALPHABET.md) separates necessary functions
from encoding choices and records the rejected alternatives.

[Heuristic] Eleven is the smallest alphabet found and certified here.
A ten-letter design would need a different encoding or mechanism; it is
not excluded. The express fold was retained, not replaced.

Reproduce from this folder with Python, NumPy and SciPy.

    OPENBLAS_NUM_THREADS=1 python3 run_certificate.py
    OPENBLAS_NUM_THREADS=1 python3 checks.py
    OPENBLAS_NUM_THREADS=1 python3 selector_checks.py
    python3 grammar_audit.py
    python3 alphabet_audit.py

[Open] Nonzero-port selection needs a common-bulk-subtracted dressed-defect
comparison. The final c=11 inference also needs volume-uniform infrared
matching and the multichannel stability input. No independent referee,
network command, subagent, or cluster job was used. No file outside
work_C84 was modified. [STATE.md](STATE.md) carries the local handoff.
