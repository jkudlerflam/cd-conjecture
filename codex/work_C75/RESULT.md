# C75 result

[Proved finite ingredients; conditional infrared conclusion] The alphabet
can be reduced from 23 to 15 letters while retaining a finite n=15 member
with 30 simple Fermi crossings, distinct sea centers, exact U(1), and
strict positive pair-sum margins. The cell has length 16. The finite
band certificate includes the corrections from finite conversion strength
and finite grammar penalty.

The corresponding infrared prediction is c=15>d-1=14 for the frozen
chain, or c=16>d-1=15 after adding one vacancy. These remain conditional
predictions. No central charge was measured, and a complete microscopic
counterexample is not claimed. Fifteen is the smallest alphabet found
here, not a proved minimum.

## Construction

[Definition] The alphabet is

    L_a, R_a, D_(A,a), D_(B,a), X_(A,a), X_(B,a),  a=0,1,
    M_0, M_A, M_B.

The bond is specified in full in [PROOF.md](PROOF.md), Section 1, and
implemented in [bond15.py](bond15.py). Its local dimension is 15, so
the literal bond is 225 by 225. At a generic resolved parameter point
it has 374 nonzero entries, is Hermitian, and conserves the on-site
integer charge exactly in the matrix check.

[Proved] The eight-letter saving has two mechanisms. Excited charges
-1,+1, in units where a data bit has charge 2, replace the neutral,+,-
triplet, saving five letters across the data track, express track and
markers. A vacuum domain wall replaces the explicit vacuum token, and
the two vacuum marker letters merge, saving three more. The express
fold itself remains. [ALPHABET.md](ALPHABET.md) records the retained
groups, failed literal deletions, and unresolved alternatives.

[Proved] On legal cells, the new two-species second-order Gaussian
matrices are one half of C33's common-loop matrices. The finite penalty
produces extra processes, including illegal marker excitations. They
are bounded rather than omitted. A separate uniform estimate proves
the zero-port length selector at finite penalty, for all cell lengths
and illegal sectors.

## Explicit certified member

[Definition] Set n=15, D=1000, Delta=4-2cos(pi/17), all main hoppings
and the switch to one, and

    phi_A=pi/30, phi_B=0, psi_A=pi, psi_B=pi/2,
    port amplitudes exp(-i pi/4), exp(+i pi/4).

Define d_0=1, d_1=1004, d_(j+1)=1004d_j-d_(j-1), beta=1/d_16, and
alpha=sqrt(2/17)sin(pi/17). One sufficient finite schedule is

    t_d=beta 2^-84,
    g=2^-256,
    K=2^1152,
    lambda=2^-640,
    physical data field=lambda^2 alpha^4 beta^2 g^2 2^-140.

The marker chemical potential is the exact midpoint selector in PROOF
equation (5), numerically 1.9015866376241617. The hierarchy is severe.
The alphabet saving is obtained partly by allowing a large finite K.
The small ED checks use resolved parameters, not this literal hierarchy.
The diagnostic matrix builders use floating point. The final hierarchy
is represented by the displayed powers and integer recurrences, not by
a floating-point array of its microscopic bond entries.

## Certificates and checks

[Proved, computer-assisted] The final n=15 certificate gives

| Quantity | Certified result |
| :-- | :-- |
| Simple Fermi points | 30 |
| Minimum Fermi-point separation | greater than 0.0145540 |
| Minimum distinct unordered pair-sum separation | greater than 1.1444e-7 |
| Pair sums tested, including repeated points | 465 |
| Positive Hermitian hopping margin after congruence | greater than 0.5510 |
| Additional normalized Bloch error allowed | 2^-32 |
| Finite-g error upper bound in that metric | less than 1.137e-15 |
| Finite-K error upper bound in that metric | less than 1.293e-28 |

The scalar resolvents use 512-bit integer arithmetic, followed by
rational endpoint inertia tests and exact comparisons of pair phases.
The raw roots and centers are in full_loop_bands_n15.json; the exact
brackets are in root_certificate_exact_n15.json. The bounds refer to
the quadratic projection of the exact second-order effective bond at
finite g,K. They are not a claim that the interacting bond is exactly
Gaussian. Finite low-order nonresonance is also not a proof of complete
rational independence of the Fermi momenta.

[Proved] At zero port coupling, the complete ground energy of a cell
with k data sites obeys

    mu_sel-2cos(pi/(k+2))-max(k-1,0)/(K-4)
      <= E_K(k) <= mu_sel-2cos(pi/(k+2)).

The selector's neighboring energy-density gap exceeds 2^-16. This
uniform bound proves that length 16 remains selected at the stated
finite K, including comparison with marker-free rings. The selected
cell retains its exact 2^15-dimensional code and a positive auxiliary
gap.

[Numerical] Literal-bond and source-resolvent calculations agree through
n=2 to relative error below 1.1e-15. The g^2 coefficient is Gaussian
on the entire pair code to residual below 3.2e-16. Full invariant-space
ED, including illegal states, uses dimensions 436 and 15404 at n=1,2.
All four or sixteen code levels agree with the exact second-order
Hamiltonian with error proportional to lambda^4. At n=2 the ratios
error/lambda^4 are 0.00173389, 0.00172948, 0.00172839 for
lambda=0.2,0.1,0.05. Final eigenpair residuals are below 8.0e-14.

[Rejected check] The initial unshifted Lanczos run missed code-band
multiplicities and failed the energy test. That output is retained as
small_checks_unshifted_rejected.json. The final checks use shift-invert
and a larger Krylov space. No conclusion uses the failed spectrum.

## Remaining scope

[Open] The nonzero-port common-bulk-subtracted defect comparison and
volume-uniform infrared matching are not supplied by a finite band
certificate. The one-vacancy version also inherits the parent candidate's
vacancy matching problem. The reduced alphabet has not undergone an
independent referee check. The c=15 or c=16 conclusion remains conditional
on these steps and multichannel Luttinger-liquid stability.

[Open] No smaller complete bond was found. A different fold, a two-token
cell, or merging equal-charge excited letters might improve on 15, but
their second-order structure is not certified here.

## Reproduction

Use Python with NumPy and SciPy, one BLAS thread, from this folder.

    OPENBLAS_NUM_THREADS=1 python3 run_certificate.py
    OPENBLAS_NUM_THREADS=1 python3 checks.py
    OPENBLAS_NUM_THREADS=1 python3 selector_checks.py

The computations are finite local checks, not large-scale simulations.
No cluster job, network command, subagent, or external referee was used.
No file outside work_C75 was modified. Local STATE.md carries the
handoff because the project-level STATE lies outside the allowed scope.
