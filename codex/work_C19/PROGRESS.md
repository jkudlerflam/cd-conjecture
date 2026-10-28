# Current handoff

2026-10-27. C19 finished with a partial construction. Read RESULT.md
for the full proof and the remaining task requirements.

[proved] Uniform Jacobi filters, n+1 accessible energies, exact parity,
and exchange-generated hopping plus chemical potentials. The explicit
reflection cell has D=16(n+2)2^n. A separate-end version has
D=16(n+2)^2 2^n and never puts both ports in the same energy term.
Neither D is a microscopic-site alphabet.

[numerical] The final exchange-only version passed n=1,2 two-cell
and n=1 three-cell ED. Matrix errors at second order are below 1.6e-14;
ED residuals are below 9.6e-13. Exact rational filter checks pass through
n=12. Saved-data audit passes all 21 ED records across the final cell
and its controls. Higher-order non-Gaussian terms are resolved.

[proved, limitation] A constant number k of accessible poles permits
at most 2n<=4(k-1). C18's independent exact-FF static-Bell-check
obstruction excludes that particular implementation of the new cell.
It does not prove general impossibility. No complete NN program,
exact d0, BCO adaptation, or fixed-d unbounded-c model was established.

Next useful work: compile the separate-end target using a richer
syndrome-conditioned local menu or a parity-preserving non-FF program.
Count the actual work alphabet and control errors relative to rho^2.
Do not repeat the finite-matrix construction or claim a Python loop
as a fixed-alphabet Hamiltonian compiler.

Primary reproduction: python3 exchange_field_checks.py.
Fast saved-data audit: python3 audit_results.py.
All files stayed here; no network, cluster, or subagents were used.

# C19 progress

2026-10-27. Initial context read. No earlier local progress existed.

[proved in source, used as input] Read STATUS, minimal_d0 in full, G1,
C9, programmable_cell, ir_control, C13, and C8 progress. C8 is constructing
nine- and fifteen-letter clocks and their work interpreter. C18 now has
an initial progress file and is treating idle-port boundaries and merged
selection. Those controller tasks will not be duplicated here.

[research direction] The n-mode target needs only 2n odd Majoranas,
not the full 4^n operator basis. Use an odd syndrome matrix between the
X and Z Bell syndromes, paired with an odd logical Majorana. This keeps
physical parity exact. Replace the dense Cauchy whitening of C9 by
orthogonal polynomial Laplace filters with consecutive energies.
Candidate functions are x(1-x) times Jacobi polynomials for weight
x(1-x)^2 on [0,1], x=exp(-t). Their coefficients sum to zero and their
Cauchy Gram matrix is diagonal analytically. This is not yet a spatial
nearest-neighbor construction or an alphabet count.

Next: derive the filters and parity-respecting dilation, compute the
free-fermion bond and three-cell cross terms, then test finite matrices.
All writes and computations remain inside this folder. The project
STATE.md is outside the allowed write scope and will not be edited.

## Step 1 completed, analytic graded filters

[proved] For j=0,...,n-1 use f_j(t)=x(1-x)Q_j(x)/sqrt(h_j),
x=exp(-t), Q_j=P_j^(2,1)(2x-1), h_j=(j+1)/[(2j+4)(j+3)].
Expansion in x gives n+1 consecutive energies. Sum-zero coefficients
and Cauchy orthonormality follow from Rodrigues and integration by parts.
Two orthonormal odd syndrome matrices carry the two Majoranas per
flavour. Their transpose on the left is necessary. The resulting bond
is a scalar plus rho^2 times arbitrary flavour-diagonal hopping,
including phases. Chemical potentials are separate even order-lambda^2
on-cell terms. This qualification preserves the exactly degenerate H0.

[proved] The cell grading is Z_pL Z_pR Z_sL Z_sR P_logical on the
unpadded registers. Syndrome X/Z mixing is odd, so multiplying it by an
odd Majorana gives an even cell Hamiltonian. Heisenberg port exchange
preserves the product of cell gradings. The constructed weights have
strict positivity. Their square roots are quadratic polynomials in the
nonconstant weight, because it has only three distinct eigenvalues.
An explicit reflection dilation avoids an unspecified unitary completion.

[numerical, preliminary] uniform_cell.py builds a smaller equivalent
parity-respecting dilation using numerical null spaces. n=1,2,3 cells
have dimensions 32,76,200. Code and parity errors are below 4e-16.
The second-order bond at n=1,2 agrees within 4e-15. Exact rational
filter identities passed at n=1,2,4,8,12. Full two- and three-cell ED
is running in checks.py. No spatial NN realization has been claimed.

Next: finish ED, prove the bounded-pole obstruction, and audit what
fixed-alphabet rules would still have to compute. A formula uniform
in n is not yet a small-alphabet Hamiltonian program.

## Step 2 completed, finite-cell tests and obstruction scope

[proved] A uniform reflection dilation has whole-cell dimension
16(n+2)2^n. Its only n-dependent scalar data are Jacobi coefficients,
normalizations, and independently specified flavour hopping parameters.
The weight square root is quadratic, so no exponential-size matrix
square root or arbitrary unitary completion is needed. The compact
numerical dilation has dimension 4 ceil(max(16,6n+7)2^n/4).

[numerical] Full compact-cell ED finished for (n,M)=(1,2),(1,3),(2,2),
at lambda=.08,.04,.02, dimensions 1024,32768,5776. All low-code levels
were extracted; residuals were below 2e-12. The error to second order
is cubic. The uniform reflection completion separately passed n=1,2
matrix tests and n=1 two-cell full ED at dimension 9216.

[numerical] Third-order operators are not purely quadratic. After
projecting out every scalar and Majorana bilinear on two cells, the
reflection completion has residuals 2.16e-7 (n=1) and 1.35e-9 (n=2),
for the coefficient of lambda^3. Thus “only coupling renormalization”
is established at second order, not at all orders. Standard weak-coupling
IR control is still required. The compact completion has a larger
n=2 cubic non-Gaussian term, illustrating completion dependence.

[proved] A qubit port with grading Z has only four odd syndrome-matrix
slots. With k accessible positive energies, normalization restricts
each odd slot to k-1 independent logical operators. Thus exact n-flavour
hopping at second order requires 2n<=4(k-1). A fixed number of poles
cannot work. This does not obstruct the new k=n+1 uniform family.

[proved in C18, applied here] C18's live progress now states the stronger
static-check obstruction: positive Bell checks and positive port-idle
energies annihilating the code at every time give only scalar logical
source resolvents. This includes initialized ancillas. Do not duplicate
its controller proof. The present nontrivial filters cannot be compiled
by that exact FF interface. Non-FF programs or richer port-touching
checks remain possible, so no whole-class impossibility is claimed.

Next: finish uniform-reflection three-cell ED using its exact invariant
endpoint sector, write the complete derivation and conditional alphabet
accounting, and leave a concrete compiler specification. No valid d0
or BCO alphabet count has been established.

## Step 3 completed, exchange-generated chemical potentials

[proved] Added fields to the unused left Y syndrome only. With
 d=(1,-1,0,...), g=p^T C d=1/[2(k+2)], and sigma=-rho^2/g,
add sigma d_e |Y><Y| tensor F to M_e^L. Odd hopping weights use
X/Z, so they have zero contraction with this diagonal term. There
are no field-field terms. The exact second-order bond is now
 a0 I + rho^2 (K_hop + I tensor F), in units lambda^2/Delta.
Every cell gets its field once on a ring, entirely from exchange.
The code of H0 remains exactly degenerate. On open chains the first
cell has no chemical field unless a separate boundary term is added.

[proved] Strict positivity follows by taking the minimum in RESULT
Eq. (16). The Y-block square root is a scalar square root of a sum
of commuting occupations, computable by the same accumulation loop
for every n. The X/Z square root remains quadratic. Thus this repair
does not restore arbitrary N-dimensional matrix processing.

[numerical] Direct second-order checks of the exchange-only field
version pass for n=1,2 at cell dimensions 96,256, errors below 1.6e-14.
Its exact invariant endpoint sector gives n=1 two- and three-cell ED
at dimensions 576 and 55296, embedded in full dimensions 9216 and
884736. All code levels are retained. n=2 two-cell ED is finishing.
The endpoint restriction is exact, not a low-energy truncation.

[proved, untested alternative] Separate flags at the two ends give
H_cell=h_L+h_R, dimension 16(n+2)^2 2^n, the same single-source
weights, and gap Delta. Each end Hamiltonian commutes with the other
end's ancillary ground projector. This avoids a term containing
both distant ports, but remains a nonlocal workspace Hamiltonian.

RESULT.md now contains the full derivation, pole-count obstruction,
C18 static-check obstruction, and conditional clock counts. Neither
a complete local energy program nor an exact d0 is claimed. Remaining
work in this run is numerical-data audit, reproduction documentation,
and a final source/claim check.


## Final audit completed

[proved and numerical] RESULT.md includes the full filter proof,
positive weights, analytic square roots and reflection, chemical-field
contraction, exact parity, absence of second-order three-cell terms,
pole-count obstruction, conditional skeleton counts, and ranked next
steps. It distinguishes the final exchange-only cell, earlier explicit-
field controls, and the untested additive-end variant.

[numerical] audit_results.py passes every saved record and the exact
filter checks. The final one-flavour two-cell error ratio under halving
lambda is 8.04965. Compilation of all six Python files passed. No
unchanged ED suite was rerun for the final prose edits.

[proved, scope] Task 3 and the spatial part of task 1 remain unfinished.
The fallback premise of general uniform-gadget impossibility was not
established. No unverified small d0 or exact BCO count is reported.
Only this folder was modified. Project STATE.md remains untouched
under the user's explicit write restriction.
