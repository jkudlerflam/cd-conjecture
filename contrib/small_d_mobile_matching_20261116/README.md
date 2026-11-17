# First mobile estimate: a charged-rail pin cost

Date: 2026-11-16. Author and AI assistance: OpenAI Codex, GPT-6 family.
No independent agent or human review is claimed.

## Result and scope

**[Proved, awaiting independent review.]** The literal d=22 mobile
candidate has a positive, volume-independent cost for pinning a selected
cell in its **charged-rail** auxiliary sector, with the full gapless
vacancy dynamics retained. The comparison uses a positive length-tilt
ensemble and also applies to its canonical minimizing vacancy-number
blocks. It does not establish the claim at every prescribed density.

For the parameters below, in the exact legal ordered-particle model,

    F_pin(s)-F(s) >= J d_ex,
    d_ex=(1-|g|-q^2)/1024 >= 7/16384.

Here F denotes a ground energy minimized over vacancy number, including
sV. It is not a thermal free energy. The reference F retains all the
same interacting bulk terms. At zero port this includes a single
charged-rail cell versus vacuum cells. At nonzero port it is a
compressed pin comparison, allowing charged virtual pairs elsewhere.

The resulting charged-cell occupation bound is

    <P_ex> <= (2B/d_ex)^2,  B=sqrt(2+2 eta)|lambda|.

This controls **charged-rail admixture**, not a full dressed auxiliary
gap. Excited vacuum-clock eigenstates remain uncontrolled. Wrong-length
cells and their common-bulk-subtracted matching error remain open.

The Hamiltonian still has one repeated Hermitian nearest-neighbour
bond, one-site translation symmetry, d=22, and size-independent
couplings. The legal subspace is exactly invariant at finite K. No
claim that it is the complete ground sector at finite K is made here.

## What is new in the estimate

A contact lower bound alone does not imply a positive energy for a
moving auxiliary head: the head can be correlated with vacancies.
The proof instead bounds the joint observable

    P_head-contact=sum_(i in cell) Pi_i A_i,

where Pi_i is the quantum head-position projector and A_i says that
the following material particles are adjacent. A controlled
vacancy-removal channel gives <P_head-contact>>=1/1024 in the pinned
interacting ground state. Replacing the charged head by a zero-bit
vacuum head at the same position lowers energy by at least
J(1-|g|-q^2)<P_head-contact>.

The channel retains the entire vacancy marginal and all other cells.
The unchanged bulk cancels exactly. No charge-hopping link is cut, no
contact is replaced by its mean, and no vanishing finite-ring charge
gap is used as a perturbative denominator.

The sufficient assumptions are the parent's frozen parameter conditions,
f>=0, D>=q^2, b=n+1, and

    t/2<=s<=3t/4,
    J(n+2)M<=t/32,
    M=Delta+D+|mu_sel|+f+4+B.

The latter is a conservative local material-bond norm bound. Constants
may deteriorate with n and with the large finite D. Choosing lambda
first and J/t afterwards preserves the parent parameter order.

## Exact representation, scalar checkpoint, and physical ring

[PROOF.md](PROOF.md) derives the representation directly from
`mobile_bond`:

    H=H_charge+J sum_j A_j h_j-mu_occ N.

Every material onsite term is contact-weighted by
(A_(j-1)+A_j)/2, including the marker chemical potential, q^2 terminal
penalty, auxiliary masses, and positive readout. The occupied-number
chemical potential is unweighted. The leading average is
alpha(rho)=rho-sin^2(pi rho)/(pi^2 rho), per **material** site.
The corresponding material correction per physical site has an
additional factor rho. The balanced (n-1,n+1) versus (n,n) first-order
comparison keeps both occupied length and vacancy number fixed.

**[Proved scalar checkpoint.]** For a free canonical gap bath, forbidding
one tagged contact is exactly isometric to removing one vacancy. Its
compressed energy cost Delta_p tends to the positive charge pressure,
even though the full bath gap vanishes. A two-block estimate gives

    E(H_charge+u A_0)-E(H_charge)
      >= u Delta_p^2/(Delta_p^2+64t^2),  0<=u<=Delta_p/2.

For N>=128, density in [1/3,2/3], and 0<=u<=t/16, the bound is at least
u/4097. This is a fixed-vacancy result, distinct from C92's tilted
quantum-source contact theorem. It does not justify replacing u by an
arbitrary extensive operator-valued clock problem.

**[Proved ground-energy boundary estimate.]** The physical cyclic
quotient is retained using its exact total-momentum twist. A large gauge
variational bound limits the difference between physical and unrestricted
rooted ground energies to 48t/N. A comparison has one global allowance
of at most 96t/N. This vanishes as N grows and is not charged per defect.
For example, the physical single-charged-cell cost at zero port is at
least J d_ex-96t/N.

At nonzero port, physical ground states selected by fixed s obey the
cell-averaged occupation version

    <P_ex> <= [2JB/(J d_ex-96t/N)]^2

when the denominator is positive. At these canonical minimizing blocks,
the density stays in (1/3,2/3) for N>=32 under the stated schedule.
This does not prove uniqueness of the limiting density or equivalence
at every specified density.

A compatible sufficient allowance is

    B<=min((1-q^2)^2 q^(2n+4)/10, d_ex/16),
    0<J/t<=1/[32(n+2)M].

For N>=192t/(J d_ex), the physical charged-cell occupation bound is at
most 1/16. These estimates hold within the invariant legal space for
any finite K. Selecting that space against illegal competitors remains
a separate task.

## Numerical checks

**[Numerical corroboration; not interval certificates.]**
[checks.py](checks.py) and [checks.json](checks.json) retain:

- An independent mobile-bond rebuild and complete physical/gap spectrum
  comparisons at n=1,2, dimensions 168 and 1104. The largest difference
  is 2.05e-14. Both onsite contact halves and every cyclic momentum
  sector are included.
- The canonical contact compression isometry and scalar impurity
  inequalities through N=V=6.
- The charged-head replacement operator inequality for all 56 binary
  local contact patterns through n=3, including maximal |g|=1/2 and
  the minimal allowed D=q^2.
- Exact controlled-channel completeness, kinetic, and vacancy-number
  identities, retaining a quantum head and a nonzero hopping phase.
- Full-payload charged-pin checks at n=1,2 with both buses, including
  nonzero lambda. At J/t=0.0005 the observed pin costs exceed 0.00046;
  the conservative theorem lower bound for those parameters is
  3.94e-7. A free-bath norm enclosure rules out the unexamined vacancy
  numbers, so these small tilted minima are not an unchecked cutoff.
- Balanced-pair comparisons in exact zero-bit vacuum blocks at n=1,2.
  Their remainders have both negative sign and approximately J^2
  scaling on the sampled small rings. This is finite-size evidence;
  **no uniform o(J gamma) remainder follows from it**.

The largest recorded eigensolver residual is below 9.7e-11. The tiny
nonzero-port occupation is a numerical diagnostic, not a certified
measurement of its asymptotic coefficient.

Run from the repository root with Python, NumPy, SciPy and mpmath
(the last is imported by the independent frozen helper):

```sh
OPENBLAS_NUM_THREADS=1 python3 -B contrib/small_d_mobile_matching_20261116/checks.py
```

The checks take a few seconds (2.6 seconds in the latest saved run) and
write only their adjacent small JSON. The parent model is used for entrywise checking; the mobile
physical matrix, gap matrix, cyclic projection, and contact-source
assembly are constructed separately.

## Remaining hypotheses

The charged replacement uses the large positive Delta,D contact energy.
It cannot be applied to excited states of the vacuum clock, which have
no charged mass to remove. The readout remains positive and does not
provide a gap uniform in the number of coupled payload cells.

The replacement preserves the marker pattern. It therefore cannot
compare wrong lengths, retune their chemical potential, or establish
the requested common-bulk-subtracted o(J gamma) error. It allows
arbitrarily long **other** cells but pins a cell of fixed length b;
the head-conditioned support constant grows for a pinned long cell.
Arbitrary-wrong-length control has not been silently inferred.

The [route audit](../route_audit_20261109/README.md) and C91–C96 retain
the response, ensemble, and clock-matching limitations used here.
The finite-temperature counterexamples in C94 and C96 are not evaded by
a finite-order band calculation. No thermal all-defect theorem is claimed.

The full interacting multichannel CFT, fast-vacancy vacuum-clock matching,
positive mobile wrong-length costs, and illegal-sector comparison remain
open. Proposed c=n+1 is total nonchiral Virasoro central charge, counting
a compact boson as one. No microscopic c is measured or established.
