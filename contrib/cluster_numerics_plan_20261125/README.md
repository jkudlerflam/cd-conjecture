# Recommended IAS cluster numerics

2026-11-25. **[Proposed numerical campaign, not a new physical result.]**
Start with the simplest qubit candidate at delta=3/16, h=1.783. Resolve
the quartet's scaling and operator identity before spending the main
allocation on longer-chain entropy fits. Then measure bulk modes and
velocities. A smaller parallel campaign can test the original d=10
mobile family's unresolved contact response.

## Basis for the priorities

This plan reviews GitHub through `cc8afc10`, especially REPORT §§9.216–9.217,
the committed quartet data and analysis, and the Typhon runbook. The fixed
qubit Hamiltonian is

    H = sum_j [(19/32) X_j Y_(j+1) + (13/32) Y_j X_(j+1)
               + Z_j Z_(j+1)] - (1783/1000) sum_j Z_j.

It is Hermitian, one-site translation invariant and nearest-neighbour,
with microscopic d=2 and full fields at both open ends. The intended c
is unitary compact Virasoro central charge. An entropy slope is not by
itself its asymptotic identification.

The November 23 missing-data list has been partly overtaken: complex pair
data are committed, a physical pi component has been measured, and an
independent DMRG implementation now reaches L=256. The direct quartet
estimate is **x4 approximately 2.05 +/- 0.2**, which does not resolve
relevance against x=2. The quartet fits use r=32–56; weaker longer-distance
signals move with chi. Increasing L raises the fitted exponent while
increasing chi lowers it. The source's smaller-quadrature interpretation
assumes that the unwanted neutral contribution has rank one; that
assumption needs testing alongside convergence. No new numerical audit
of that extraction is claimed here.

Sources: [latest results](../../REPORT.md),
[quartet table](../d2_xy_yx_chirality_imps_20261118/tables/quartet_k0.txt),
[extraction](../d2_xy_yx_chirality_imps_20261118/scripts/ana_quartet.py),
[earlier audited questions](../qubit_ir_cft_team_20261123/TODO.md).

## 1 Resolve the quartet at fixed lengths

**First inexpensive step.** Reanalyse the six committed quartet datasets
with varied demodulation widths, harmonic content, reference positions
and fitting windows. Test the extraction on synthetic correlators with
known quartet exponents, multiple neutral components and oscillatory
amplitudes comparable to the data. A 48-site local demodulation window
is substantial relative to the fitted r=32–56 range; measure its bias
instead of absorbing it into a regression error bar. This needs minutes
of CPU work, not new ground-state optimization.

**First production target.** At L=256 and 384, run both parity sectors
through chi=1024,1536,2048, saving states and measurements at each stage.
Use recoverable old states if available; their large checkpoint files
are not on GitHub, so access on IAS must be established before assuming
warm starts. The independent code found opposite ground parities at
L=128 and 256. Report both sectors and the lower-energy choice, with
energy uncertainty small compared with the parity splitting.

Retain the pair and all normal/anomalous quartet quadratures. Extend the
measurement to cross-correlations between quartet shapes and candidate
neutral probes such as density and bond operators. Test whether the
putative quartet remains stable when the neutral basis is enlarged;
the smaller eigenvalue alone does not certify emergent charge four.
Run the h=1.772 one-mode control and a known two-boson control to check
the analysis. A control's larger local dimension is not a qubit example.

**Decision.** Aim for combined chi, window, operator-choice and size
uncertainty below about 0.05 in x4. A value below two supports a relevant
lock only after operator matching; its nonzero target coupling is still
needed. A value above two removes that particular weak-coupling threat.
Persistent marginality requires a flow/amplitude analysis. Add L=512
only after fixed-length convergence improves the useful distance range.

The existing [finite-chain driver](../d2_xy_yx_chirality_imps_20261118/scripts/ofd.py)
and [measurement routine](../d2_xy_yx_chirality_imps_20261118/scripts/meas_quartet.py)
provide the base. The enlarged operator matrix requires a small extension
and dense-vector checks. The driver deletes earlier-stage states within
one ladder; preserve checkpoints or measure each stage before advancing.

## 2 Identify the bulk modes and velocities

Use converged states from campaign 1 for the connected dynamic matrix
S_ab(q,omega), initially on L=128 and 256. Start with density, the two
nearest-neighbour hopping quadratures, and uniform distance-two chirality.
Measure each momentum sign separately and at least three resolved small
absolute momenta. Unfold the computational gauge and the two-site cell.

A CPU pilot can use MPS real-time evolution of local perturbations with
bulk windows; halve the time step and raise chi to assess errors, and
stop before boundary reflections. Increase time and length only if the
frequency resolution separates the candidate branches. An efficient
matrix-free excitation method is an alternative requiring implementation.
The current repository has no validated large-chain version of this
measurement, so its cost must be benchmarked.

**Decision.** Look for two independent current channels with linear
dispersion and stable positive velocities, or one channel developing a
resolved threshold. Several peaks or multiparticle harmonics are not a
mode count. Calibrate on c=1 and on a known two-boson model, preferably
with unequal velocities. Operator algebra and zero-mode sectors remain
necessary for a full CFT identification. This data can also constrain the
normalized response and mixing, addressing the old density/pair mismatch.

Do not simply increase the existing Koo–Saleur calculation's allocation.
Its two-mode control and the target-family c=1 control failed; the dense
periodic implementation costs O(L d^2 D^6) per momentum. Any further use
starts with a bounded successful control calculation, not target readings.

## 3 Test the original mobile model's missing response

For the distinct original d=10 family at vacancy number V=N, extend the
contact response gamma_N(q) and its local diagonal beyond N=10. The q=0
channel already has a uniform analytic bound; the unknown nonzero-q and
local channels deserve the computing time.

Start with the existing N=10 calculation and then N=12, using projected
sparse linear solves, reporting residuals, the q of the maximum and the
small-q versus exactly-zero-q distinction. The gap-coordinate dimensions
are binomial(2N-1,N): 92,378 at N=10, 1,352,078 at N=12 and 20,058,300
at N=14. The present Python basis/transition builder is memory intensive;
do not increase N blindly. Benchmark N=12 before committing to a
matrix-free or tensor-network extension.

Period-10 background harmonics need compatible sizes, such as N=20;
N=12 is a convergence/algorithm pilot, not an exact period-10 calculation.
Reaching those sizes needs an improved algorithm. Do not reuse the
same-grid two-particle-hole truncation away from q=0; the prior work
explicitly showed that it fails there.

**Decision.** Determine whether the relevant responses saturate or grow,
and which channels threaten sector selection. Numerical saturation would
guide a uniform bound, not prove it or establish unbounded microscopic c.
This does not transfer the modified frozen d=9 proof to the mobile model.
The smaller d=9 construction was reviewed first: its next stated gap is
all-scale interacting matching, for which larger finite-cluster runs alone
are less discriminating.

Sources: [original mobile response](../markerless_thermodynamic_selector_20261118/README.md),
[response implementation](../markerless_thermodynamic_selector_20261118/checks.py),
[modified d=9 matching](../markerless_uniform_matching_20261118/README.md).

## Cluster staging and resource gates

Read-only queries on November 25 found 63 Typhon nodes configured with
96 CPUs and 384000 MB each. The short QoS has a 24-hour ceiling; the
user's queue was empty at inspection. These are configured capacities,
not reserved resources or runtime predictions. The verified workflow is
the [IAS CPU driver](../../tools/typhon/README.md).

| Stage | Suggested initial request | Expand only after |
|---|---|---|
| Saved-data extraction audit | 1 CPU, 2 GB, 5 minutes | Measured runtime requires it |
| DMRG portability/timing pilot, L=128, chi up to 512 | 4 CPUs, 16 GB, 1 hour | Correct energy, reliable checkpoints, per-sweep timing and peak memory |
| Fixed-length higher-chi ladder | Initially 8 CPUs, 64–128 GB, up to 12 hours | Pilot predicts the next stage fits; choose memory from measurements |
| Dynamic-response pilot | 4–8 CPUs, 32 GB, 1 hour | Controls pass and time/chi convergence improves resolution |
| N=12 mobile response pilot | 1 CPU, 32 GB, up to 4 hours | N=10 baseline and measured construction/memory scaling |

These are proposed allocation caps, not measured completion times.
Benchmark 1 versus 4 versus 8 numerical threads on a representative
sweep before scaling up. Parallel independent lengths, sectors and probes
are a better initial use of nodes than assigning 96 CPUs to one serial
solver. Make the matrix-building mobile code memory efficient before
attempting its N=14 extension.

For each job use a fresh name with `python3 tools/typhon/typhon.py submit
RUN RECIPE`, listing every actual working-tree source/input. Keep large
restart states on IAS outside fetched results, and retain small raw
correlators, convergence logs and provenance. The current finite-chain
driver needs a wrapper to set its wall-time/checkpoint deadline and route
outputs correctly; a scheduler completion alone is not scientific
convergence. Use `status`, `logs`, then `fetch`, and check accounting,
exit status and numerical residuals before reporting completion.

## Status and review

AI assistance: OpenAI Codex, GPT-6 family, repository and code review.
This proposal was checked against the latest committed results, source
interfaces and live read-only scheduler output. No new research jobs
were submitted, no new physical finding is claimed, and there was no
independent audit of this plan. It updates priorities from the earlier
to-do list without changing the underlying results.
