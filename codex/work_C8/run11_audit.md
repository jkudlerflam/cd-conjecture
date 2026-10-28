# Eleventh-run composition and adversarial audit

Status: nonconclusive for the original c-d conjecture. This run closes
three energy/length obligations of the explicit first-order candidate,
with proofs and finite controls. It does not close the low-frequency
massive-sector elimination or the subsequent Gaussian infrared step.

## 1. Which Hamiltonian is now preferred

Use certified_periodic_bond.CertifiedPeriodicBond, with epsilon_max=1,
and multiply its perturbation by a fixed epsilon in the interval below.
The only change from PeriodicBond(flat_vacancies=True) is the larger
diagonal penalty for a forbidden selector bond. Its site dimension is
still d=96*446400(q+1)+48, with m=q^2/8 projected payload modes and
one projected length mode. The clean-readout machine, all its gates,
ports, densities, the flat wrapper, and the legal perturbation are
identical. No claim is made for the adequacy of the old smaller penalty
under the new finite-coupling sector proof.

The new penalty depends on the fixed q and epsilon_max, not on ring
length. Its full-space meaning is specified in selector_repair.md.
At q=4 the old coefficient was 4503.12845888; the new one is
1.67536545854e13 for the deliberately generous epsilon_max=1.
Its lower defect-gap bound is 8.37682729270e12. A smaller fixed
epsilon_max can reduce this coefficient without changing any legal
matrix element. This large penalty is not a measured many-body gap.

## 2. One explicit positive coupling interval

Fix the inherited positive atomic gap Delta, wrong-length slope gamma,
and full unscaled per-cell perturbation bound J<=2B_flat. The projected
length hopping is J0=1, and the marker coefficient is the existing
kappa=2/pi-e_pay. Set

    C=max(2/Delta,1/gamma), eta=1/5, delta=1/10,
    g_delta=(2/pi)[1-cos(pi/10)],
    D0=2C Delta+2J+4,
    epsilon_* = min{1, Delta/(2D0),
                    g_delta Delta/(4J^2),
                    eta/(4CJ), eta gamma/(16J)}.          (1)

[proved] Every fixed 0<epsilon<epsilon_* meets all new hypotheses.
For sufficiently large L, the full ground state has no selector
defects, has positive marker number, and its sector satisfies

    2/5 < (L-bM)/M < 3/5.

The proof combines the unitary defect repair, the negative-base-energy
trial excluding marker-free rings, and the pressure comparison. No
claim that this interval is sufficient for a CFT is made.

The individual coefficients can be very small or large. Finiteness
and positivity follow after fixing the finite q, b and core. This
ordering of limits agrees with the earlier construction. The estimates
do not allow epsilon or q to vary with L.

## 3. A full microscopic low-window cutoff theorem

[proved] Fix a positive Omega<=gamma/2 below the selector-defect gap.
Increase L until the marker-free sector is also above E0+Omega and

    g_delta L/[2(b+1/2+delta)] > 2beta+B_L+Omega/epsilon.   (2)

Here beta=J0+B_pay and B_L is defined in marker_density_bounds.md.
Every sector meeting [E0,E0+Omega] is then positive-marker legal and
has density in the same interior interval. This follows from the
same lower bound as marker selection, with Omega added to the trial
upper bound. Sectors with fixed M=1,2 are excluded by their linear
wrong-length cost. Each retained sector has its own ground energy
at least E0; its part of the global window lies within Omega of that
sector ground. The canonical spectral-window theorem therefore applies
to each of them.

Let C_R retain exactly those legal configurations all of whose cells
have length less than R=b+r, r>=2. It is a well-defined physical
projection, independent of the choice of reference marker. Define

    a_r=(2epsilon J/gamma)^(r-1)/(r-1)!,
    z_r=sqrt(L/2) a_r.

Then the global physical low-energy projection satisfies

    ||(1-C_R)Pi_[E0,E0+Omega]|| <= z_r.                  (3)

Indeed each allowed sector has M<=L/2 cells, and the cell-wise
operator bound plus the union bound gives sqrt(M) a_r. Direct sums
over M do not add their norms. There is no factor counting the
number of eigenvectors in the window.

For z_r<1, the k-th eigenvalue, counting multiplicity from the bottom,
of the restricted physical Hamiltonian H_R=C_R H C_R obeys

    0 <= e_k(H_R)-e_k(H)
      <= [Omega z_r^2+(L/2)epsilon J z_r]/(1-z_r^2)       (4)

whenever e_k(H)<=E0+Omega. The proof is the same min-max comparison
as canonical_length_bounds.md. The atomic part, including the common
base scalar, commutes with C_R. The interaction norm on each relevant
sector is at most (L/2)epsilon J. Unused illegal or marker-free sectors
are invariant, so removing them creates no extra off-diagonal term.

Taking r=A log L/log log L with any fixed A>5/2, rounded up and with
L sufficiently large, makes the right side o(1/L). Larger A produces
any desired fixed inverse power of L. All energy comparisons use the
one fixed microscopic bond.

This establishes a controlled finite-volume length truncation at the
scale where conformal towers would be resolved. It does not establish
that those towers occur. The truncated cell dimension grows with L,
and the clock/internal excitations inside each allowed length remain.
The cutoff is a cell-coordinate restriction with the original physical
quotient; it is not another microscopic nearest-neighbor chain of a
claimed smaller on-site dimension.

## 4. Canonical bare defect occupation away from half filling

[proved] The chemical-tilt method also strengthens the canonical
occupation estimate. For any rho=N/M in [0,1], choose
alpha=-2J0 cos(pi rho), mu=epsilon alpha. The exact reference has
a physical canonical trial of energy at most
epsilon M[kappa+e_pay+e_len(rho)]+2epsilon beta, by the open-chain
construction. Its grand-canonical reference energy is at least the
same bulk value minus mu N-2epsilon beta. Hence the reference-sector
excess needed in defect_energy_bounds.md is at most 4epsilon beta.
Writing g_mu as in marker_density_bounds.md gives

    <sum_j Q_j>/M
       <=4epsilon^2 J^2/g_mu^2+8epsilon beta/(g_mu M).     (5)

This is a bare atomic defect bound for every such canonical ground
state, and in particular for the full ground sectors just selected.
It does not require rho to be the nearest finite-volume value to 1/2.
The total all-code overlap still need not stay positive as L grows.

## 5. Audit of possible shortcuts and hidden assumptions

The raising operator R between the two code lengths is not a physical
nearest-neighbor gate. It is used only to build a variational trial
in an exact atomic coordinate space. The Hamiltonian itself is still
the already specified repeated microscopic bond.

Changing L changes the physical quotient. The proof never equates
the ground energies of two arbitrary seam twists. It first proves
canonical increments for open cells, then uses a pinned-marker state
and explicit boundary-norm costs. Reference relabelings cannot give
Hamiltonian cross terms in that trial because every cell has length
at least two and the reference marker moves by at most one.

The long-tail inverse is a positive compression of H-E_ground after
a canonical energy comparison. It is not the global bare-Q inverse
from the failed run-9 argument. The spectral-window extension uses a
separated Sylvester equation and therefore controls a whole subspace,
not just individual eigenvectors with an uncounted dimension factor.

The pressure proof does not call the unbounded length operator a
bounded perturbation. Its off-code part is absorbed into A_mu, with
the explicit positive lower bound (3) of marker_density_bounds.md.
The reference pressure is minimized over all length-bit charges and
all reference phases before its Legendre transform is used.

The selector repair compares energies by an exact local unitary and
a norm cost per defect. It needs no extensive remainder estimate and
does not assume that the legal chain is gapped. The strengthened
coefficient is necessary for this proof; it is not silently ascribed
to the old oracle. Its numerical class name certifies this selector
statement, not a central charge or an infrared reduction.

The marker-density conclusion allows sector degeneracies and level
crossings. Neither the norm estimates nor the O(sqrt(epsilon)) density
bound proves a unique ground state on every ring or a positive dressed
compressibility. These cannot be inferred from a convex reference
pressure alone.

## 6. Attempt at the remaining local action, and where it stops

Re-read the connected-excursion argument in
../../contrib/research_continuation/IR_AUDIT.md Section 9 and its
independent review. Its geometric majorants are explicitly not proved
there. The present run supplies exact coordinates, canonical controls
and a full-sector comparison that those arguments lacked, but it
still does not construct their irreducible kernels.

A tempting next step is to apply finite-dimensional elimination after
(4). This does not close the argument. To resolve energies of order
1/L, the cutoff grows with L. A theorem with constants depending on
the cutoff dimension would need a new uniformity proof. More seriously,
even at a fixed finite dimension the entire retained many-body band
is extensive and overlaps excited bands. A global Schrieffer-Wolff
inverse is still unavailable. None of the tail bounds separates
retained propagation from operator mixing in the connected action.

The previous exactly dressed XX source control still applies: an
ordinary connected bare ancillary correlator can have a small algebraic
tail. Removing its scalar expectation does not produce the required
massive kernel. Filtering by the full interacting dynamics gives a
useful operator split, but not independent local retained variables
or a constructive effective action. The new length cutoff does not
change this issue.

Consequently the remaining inference to c=m+1 still requires
(i) a local low-frequency elimination with uniform connected support
and time bounds after retained-operator subtractions, and (ii) the
Gaussian infrared argument inside the rational target's momentum,
filling and stiffness margins. Treating (i) as a conclusion of the
new energy bounds would repeat the gap identified in the prior audit.
No new CFT conclusion, general upper bound or fixed-d unboundedness
result follows from this run.

## 7. Verification and next steps

Run, with single-thread BLAS,

    python3 canonical_length_checks.py
    python3 marker_density_checks.py
    python3 selector_repair_checks.py

The first passes 18 complex atomic sectors, two reference-isometry
checks and four physical scalar rings through dimension 7007, including
the pinned trial, tail recursion, low-state operator norm and spectral
cutoff. The largest residual is below 1.1e-14. Its
eight-state spectral tests use the span of the eight computed lowest
eigenvectors; they do not claim an exhaustive numerical diagonalization
of every state in an entire fixed-width thermodynamic window.
The second passes 34 physical sectors through dimension 7007, five
chemical tilts per sector, 10001 Legendre checks and 967417 marker-count
inequalities. The integer exterior-density test includes equality at
both boundary densities; its smallest margin is -8.9e-16 from roundoff.
The largest eigensolver residual is 1.12e-12. Its small ED rings are
below the conservative analytic large-L threshold, which the output
marks explicitly. The third passes 21840 arrow words, 24 quantum
selector sectors through block dimension 256, 8192 differential
physical-bond columns and 3232 active columns for each of the base and
perturbation. The legal columns are unchanged exactly, and relative
Hermiticity errors are below 3.2e-17.
All are control calculations, not microscopic CFT measurements.

The ranked continuation is (i) the uniform irreducible action, using
the now established full-sector and canonical bounds; (ii) dressed
payload fillings and the complete Gaussian quadratic form; (iii) an
independent audit of the controller and the new proofs; (iv) a smaller
implementation permitting two calibrated CFT estimators. The fixed-d
C19 and later spatial-bus program remains separate. Their current
source files do not supply a completed all-n conformal construction.
