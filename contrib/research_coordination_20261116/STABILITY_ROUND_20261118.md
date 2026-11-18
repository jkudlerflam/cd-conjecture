# Long-distance stability of qubits and markerless chains

Date: 2026-11-18. Author and AI assistance: OpenAI Codex, GPT-6 family.
Status: research assignments, not new physical results.

The user requested deployment of at most four project chats following the
discussion of decisive tests for the asymmetric qubit candidate. The two
standing priorities remain large c at microscopic d=2 and unbounded c at
fixed small d. This round reuses exactly the four completed GPT-6 Astra Max
chats. Their previous results are in the
[focused-round synthesis](FOCUSED_RESULTS_20261118.md).

Two chats now work on the same asymmetric qubit Hamiltonian, separating the
locking amplitude from the dressed-mode and response calculation. Further
three-valley exploration is deferred. Two chats pursue the remaining
interaction and ground-sector gaps of the markerless construction. The new
all-n frozen clock/readout and the original mobile model remain distinct.

## Common qubit target and matching conventions

Use the existing p153c physical bond, as defined by
[`work_C195/models.py`](../../codex/work_C195/models.py):

```text
w = exp(1.53 i)/80
Jx = (1 + Re w)/2, Jy = (-1 + Re w)/2, D = Im w/2
Jz = 10, hz = -19.97589125
H = sum_j [Jx XX + Jy YY + D(XY-YX) + Jz ZZ]_(j,j+1)
    + hz sum_j Z_j.
```

All coefficients, including the field and endpoint convention on an open
chain, stay fixed as L grows. In the established convention Delta=1,
V=40 and |t2|=1/80. A neighboring parameter interval is a family of fixed
Hamiltonians, not a license to retune the field separately at each L.
The microscopic symmetry is parity, not particle-number U(1).

Both qubit chats must specify the dressed-field convention, physical-site
momentum units, matching scale, energy and velocity units, and uncertainty
of every extracted quantity. They may read each other's saved checkpoint
results but must not wait for each other or write into each other's folder.
Unknown inputs stay explicitly unknown. The symmetry-line two-field vertex
lattice is not automatically the complete off-symmetry lattice.

## A. Match the asymmetric quartet and compare locking with detuning

Chat: `01a1ef70-d7d6-71f1-9d56-05c96f7513d9`.
Owned output: `contrib/qubit_detuned_quartet_20261118/`.

The previous [quartet calculation](../qubit_stable_two_modes_20261118/README.md)
establishes finite interacting amplitudes and lifting of a special leading
cancellation. The next task is a finite-density connected matching
calculation at the asymmetric bond, not another vacuum crossing.

**First calculation:** express the low-energy quartet vertex as a connected,
amputated spectral or correlation quantity in the actual interacting medium,
with a defined projection and counterterms. Compute a new controlled term,
bound, or finite-density estimate of its amplitude. Include the additional
pair sources admitted by asymmetric hopping and keep their perturbative
orders. Use adequate saved data where available and bounded calculations
to validate the matching; do not infer the thermodynamic coefficient from
a few-particle matrix element or substitute free intermediate denominators.

Then compare the running coupling with the dressed momentum cutoff. The
schematic weak-coupling lengths are a0 |y|^(-1/(2-x)) and 1/|delta k|;
they are not a verdict without the coupling convention, running dimensions
and uncertainty. Supply a defensible stable/unstable/undecided conclusion
for the inputs actually determined. If a full value is inaccessible, a
derived bound or exact matching identity with a quantitative missing input
is useful progress; an order-one guess for C4 is not.

Read `contrib/qubit_stable_two_modes_20261118/`,
`contrib/v4_locking_operators_20261118/`, `codex/work_C179/`,
`codex/work_C191/` and the actual p153c records in `codex/work_C195/`.
The old qubit production registration is a historical protocol for sym40;
it must not silently be applied to p153c or amended in place.

## B. Determine the two dressed modes at the same finite-density point

Chat: `01a1ef71-2c80-74b1-9b8e-4a1ed526b266`.
Owned output: `contrib/qubit_finite_density_modes_20261118/`.

This chat moves from the obstructed three-valley route to the most promising
two-valley candidate. It owns the medium response and momentum calculation,
while chat A owns the connected number-changing amplitude.

**First calculation:** derive and evaluate the asymmetric finite-density
equation of state or response matrix sufficiently to test occupation of
both valleys, positive compressibility and relative stiffness, the dressed
momentum mismatch, and the complete Gaussian covariance including mixing
and unequal chiral velocities. Audit existing long-chain correlations for
adequate bulk windows and controls before extracting exponents. The
symmetry-line K values, bare band minima, raw valley peak weights and raw
spin particle count cannot substitute for the dressed quantities.

Use the correct chiral metric and rederive the allowed operator lattice at
the asymmetric point, including new pair sources and mixed vertices. Return
either a physically matched parameter region with explicit uncertainty, a
concrete instability, or a new response/matching calculation that identifies
the exact missing input. A formal stable covariance chosen independently of
the Hamiltonian is not a result. A numerical lack of convergence remains
undecided.

Where existing data cannot decide the issue, produce a calibrated, dated
measurement protocol for two independent bosonic modes and the minimum new
observables required. Do not restart production. Read the C191 frame and
boundary cautions, the corrected locking audit, the latest quartet
diagnostics, and C190/C193/C195 records. The planned interface to chat A is
the covariance and its errors, actual momentum mismatches, velocities and
field/vertex normalization at the same physical point.

## C. Control the interactions generated by the all-n nine-state bond

Chat: `01a1ef71-7c28-7d71-8156-6188bfdc3250`.
Owned output: `contrib/markerless_uniform_matching_20261118/`.

The [all-n band proof](../markerless_all_n_20261118/README.md) now exists.
The next obstacle is an interacting effective Hamiltonian controlled
uniformly in the number of cells. Work on the modified nine-state alphabet,
with its actual biased clock, L1 readout and positive simple-band margins.

**First calculation:** obtain a local, connected expansion of the eliminated
clock/charged sector and calculate the first nonquadratic terms, including
connected multi-cell port processes. Seek a remainder bound in an explicitly
defined local interaction norm, with constants depending on fixed n and its
couplings but not chain length. Retain all generated terms and include the
actual readout's effect on the clock and source visibility. Independently
check the first connected coefficients on bounded literal cell clusters.

Compare the interaction and remainder scales with the smallest physical
mode velocity and the nonresonance margins, not with the microscopic gap
alone. Use this to give a finite nonzero coupling schedule inside a stated
multichannel stability criterion if justified. Distinguish a finite-order
normal form from an exact thermodynamic elimination. A global Kato bound
on an extensively broad code band is not volume-uniform matching.

State separately what is established within the selected frozen sector and
what needs a nonzero-port ground-sector proof. If the earlier
[bulk-subtracted defect argument](../small_d_defect_stability_20261116/README.md)
transfers, verify its hypotheses for the conserved two-site boundaries and
zero-data cells; do not import the d=20 result by dimension substitution.
No extra states or mobile charge-field shortcut are allowed. Report the
additional hierarchy imposed by each bound.

## D. Establish or refute thermodynamic mobile selection at half occupancy

Chat: `01a1ef71-ce96-7703-aa33-1337b9d16956`.
Owned output: `contrib/markerless_thermodynamic_selector_20261118/`.

Start from the original corrected contact-weighted d=10 model and the
[finite mobile result](../markerless_mobile_control_20261118/README.md).
Keep the gapless vacancy fluid explicitly. The finite N=20,V=2 certificate
does not establish selection in an arbitrarily long bath or at rho=1/2.

**First calculation:** establish or refute the static susceptibility bound
(22) in that contribution's `PROOF.md` at V=N, then attack its relative
curvature bound (24) for the balanced replacement (10,10)->(9,11) and the
weak head-number replacement of nine length-10 cells by ten length-9 cells.
Include the opposite head-number direction when extending the result.
Every comparison has identical material/vacancy totals and the same
interacting background outside the patch.

Analyze potential long-wavelength and commensurate singularities rather
than assuming a volume-independent vacancy gap. If the unweighted
susceptibility bound fails, test the actual source-projected combination
in (21), where cancellations may remain. A uniform integrated second-
difference estimate is acceptable when level crossings prevent pointwise
curvature. Determine what survives at nonzero J, not just the coefficient
at J=0. A counterexample to the proposed bound is valuable if it identifies
the mechanism and a concrete repair within d=10.

Give a justified J/t allowance for a declared comparison if possible;
otherwise isolate the reduced missing inequality and compute a new term
or obstruction. Keep zero-port selection, finite ports, charged gaps and
the separate dressed-band budget distinct. The conserved-charge field
repair applies to this original model, not to the all-n L1 readout.

## Scope, compute and publication

All four chats retain GPT-6 Astra, Max reasoning (`gpt-6-astra`, `max`),
as requested by the user for these research chats. No additional chat or
subagent is to be created during this deployment.

The previous battery-saving stop on long DMRG and background production
remains in force. Use analytic work, existing records and bounded one-worker
checks. Preparing a converged production protocol does not authorize its
launch. Do not impose microscopic U(1) where only parity is exact, fit
unresolved data, or present finite rings as thermodynamic proof.

Each chat owns only its new folder and preserves prior contributions. The
coordinator owns shared README, STATUS, REPORT and coordination files. Read
AGENTS.md and CONTRIBUTING.md; serialize complete Git mutation transactions
with the atomic directory lock `/tmp/cd-conjecture-git-sync.lock`, record the
owner and release only your own lock. If another owner is publishing, keep
reading and deriving rather than deleting the lock. Never stash, discard
or stage another chat's edits. If their edits prevent rebase, fetch and
report the integration still needed.

Simplicity remains a scientific objective. Count physical states, independent
real controls and literal terms separately; remove constants, telescoping
terms and common scale before comparing coefficient ratios. Return a concrete
calculation with hypotheses, uncertainty and reproduction instructions.

## Deployment record

The coordinator will record acknowledged redeployment and observed progress
here. This round reuses existing chats rather than creating additional ones.
