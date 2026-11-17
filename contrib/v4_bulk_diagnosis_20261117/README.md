# V=4 bulk-phase diagnosis

**Stopped at the user's battery-saving request; central charge UNKNOWN.** This bounded round tests the
large density modulation of the fixed nearest-neighbour qubit Hamiltonian
`sum[(9/16)XY+(7/16)YX+ZZ]-(7/4)sum Z`, with full endpoint fields.
No entropy estimate near 2.5--3.1 is treated as a measurement of c.

Author: OpenAI Codex, GPT-6 family, 2026-11-17. Numerical work reuses the
unchanged C190 parity solver (Claude Opus 5.5 agents). New density, bond
and pair contractions are independently checked against literal L=8 ED;
there is no independent human review or second DMRG implementation.

[Registration](REGISTRATION.md) fixes the size/accuracy sequence,
observable definitions, controls, conditional extensions and resource
cap. [Sequential decisions](DECISIONS.md) preserve adaptive choices.
Saved L=96 states are reused; new production is restricted to two workers
with one BLAS thread each. Outputs are JSON; state checkpoints stay local
and ignored. Recommendation: **missing test**. Finish the interrupted accuracy checks
and establish stable scaling beyond these sizes before identifying c.

## Results through completed chi=256 and interrupted chi=320

Both L=48 sectors pass the existing accuracy gates and the new observable
stability checks. The reused L=96 states also retain their prior local
gate passes. Their central-half density metrics are:

| L | Parity | Density RMS | Profile amplitude | q/pi | Connected peak / M |
|---:|---|---:|---:|---:|---:|
| 48 | even | 0.11211 | 0.14892 | 0.88770 | 0.01678 |
| 48 | odd | 0.11982 | 0.15932 | 0.88721 | 0.01786 |
| 96 | even | 0.10851 | 0.15767 | 0.88428 | 0.01111 |
| 96 | odd | 0.10588 | 0.15385 | 0.88428 | 0.01084 |

Connected density RMS at separations M/3 to M/2 falls from about
0.0115--0.0124 to 0.0064--0.0066. This is distinct from the large
one-point modulation. The wavevector's fine grid is not evidence for
commensurability at that precision; the physical resolution is limited
by the finite observation window.

Both L=144,chi=256 runs are now complete. Density RMS is about 0.0923
and profile amplitude about 0.128 in both sectors, with q near 0.880 pi.
The entropy chi-step gate fails (0.01894/0.01840), while the other local
gates pass. The registered chi=320 continuations pass the resource
benchmark but were stopped to save battery: odd parity has four saved
sweeps, even parity has no new completed sweep. Neither is marked done.
The independent L=96 start is complete:
energy differs by 2.33e-11, density/bond profiles by less than 1.6e-8,
and state overlap has magnitude 0.999999999875. All registered restart
checks pass. The current boundary-probe amplitude trigger is not met.

At the user's subsequent request, a **post-hoc** central-charge estimator
has been added in `entropy_analysis.py`. The lower-variational-energy
branch gives c_ent=2.818 for L=48 to 96 and 2.454 for L=96 to 144
using the established W=6 average and 6 delta S/log(size ratio).
The three-size fit gives 2.698; its geometry-corrected version gives
2.631. These are finite-size entropy slopes, not an identified
thermodynamic c or the non-unitary quantity c_eff. The large-size chi
gate failure and size/window/sector dependence remain explicit.

The same estimator on fixed-parameter exact free controls gives
0.995 and 1.954 for known c=1 and c=2 at L=96 to 144. The new valley
channel XX+YY has decreasing uniform normalized weight, about
0.064--0.065, 0.026, and 0.0157 at L=48,96,144. Small periodic-ring
spectra and correlations are saved separately. Periodic entropy now uses
a pure translation eigenstate, with the projection weights and residuals
recorded explicitly. The ground density jumps from about 0.159 at L=12
to 0.288 at L=16; these rings have large size/sector effects and do not
determine c.

`audit_results.py` checks covariance identities and positivity, pair
subtraction/Cauchy bounds and agreement of total energy reconstructed
from density/bond profiles with the independent energy contraction.
The saved stage/baseline observations pass at this checkpoint. Final
counts and all completed/running tags are recorded in `audit_results.json`.

## Numerical central-charge estimate at stopping

The latest **provisional finite-size estimate is c_ent = 2.57** from
L=96 to 144, using the four saved odd-sector chi=320 sweeps.
The last energy change is 1.15e-8, above the 1e-8 stationarity gate;
the strict stopping rule has not passed. This estimate uses the last
saved sweep entropy profile, not a completed stage measurement.

| Estimate | c_ent | Status |
|---|---:|---|
| L=48 to 96, W=6 | 2.8184 | Completed states pass local gates |
| L=96 to 144, chi=256, W=6 | 2.4542 | Largest-size chi gate fails |
| L=96 to 144, interrupted chi=320, W=6 | 2.5669 | Four saved sweeps; convergence incomplete |
| Same interrupted pair, averaged-log-chord | 2.5406 | Geometry sensitivity |
| Three-size log(L) fit, interrupted chi=320 | 2.7356 | Size drift remains |

At the interrupted highest accuracy, W=1 and W=12 give 2.3372 and
2.6471; a fixed central-quarter fraction gives 2.5790. W=12 with
averaged-log-chord geometry gives 2.5475. These are sensitivity checks,
not confidence intervals. The single midpoint is much more oscillatory
(3.497), so retaining all cuts and windows matters. The report does not
identify these finite-size slopes with a thermodynamic Virasoro c.

The bulk diagnostics through completed chi=256 are compatible with a
decaying boundary-induced pattern on the tested sizes: density profile amplitude falls from about 0.154--0.158
at 96 to 0.128 at 144, and normalized connected density weight also
falls. Long-distance connected density RMS falls from 0.0115--0.0124
at 48 to 0.0064--0.0066 at 96 and 0.0057--0.0058 at 144. Full-bond
correlations likewise decrease. The secondary XX+YY peak weight divided
by window length is about 0.143--0.147, 0.072, 0.048, while its uniform
weight also falls. These finite sizes do not exclude a smaller plateau
or a crossover. Connected pair amplitudes are anchor- and size-dependent;
no pair exponent is fitted. Explicitly sourced disconnected pair weight
is stored separately.

The registered optional boundary probe was not triggered (amplitude
ratios 0.813/0.835, below 0.85), so no altered-boundary or field scan
was performed. All twelve completed observations and twelve secondary
channel records pass the internal consistency audit. The independent
start passes. The missing tests are completion of both interrupted
chi=320 checks, a further converged size sequence, and independent
evidence for positive-velocity conformal modes. No robust critical
parameter neighbourhood or microscopic counterexample is established.

All compute jobs from this chat are stopped. Checkpoints are preserved
locally; production is not scheduled to restart automatically.

## Entropy convention and interpretation

All entropies use natural logarithms. For an open critical chain the
conformal fit is

`S(l,L) = (c/6) log[(L/pi) sin(pi*l/L)] + constant + corrections`.

The coefficient is c/3 for a periodic chain. These conventions follow
[Calabrese and Cardy, equations (1)--(2)](https://arxiv.org/abs/hep-th/0405152).
The primary finite-size diagnostic is
`c_ent(L1,L2) = 6 [Sbar6(L2)-Sbar6(L1)] / log(L2/L1)`.
W=6 averages the thirteen cuts from L/2-6 through L/2+6, inclusive.
`entropy_analysis.json` also retains W=0,1,12; a central-quarter
fraction; averaged-log-chord denominators; fixed-parity sequences; and
profile fits with three boundary exclusions and two smoothing widths.
The exact controls and synthetic conformal profiles validate the
normalization and averaging. Their finite-size errors are not assumed
to be a universal correction that can simply be subtracted from V=4.

A local accuracy pass is not a proof of the global ground state or a
CFT. Parity selection uses the lower converged variational energy at
each size; fixed-parity slopes remain available to expose crossings.
No fit residual is represented as a statistical confidence interval.
Three sizes cannot identify both an asymptotic c and unspecified
finite-size corrections. Positive thermodynamic mode velocities and a
robust critical parameter neighbourhood have not been established.

## Reproduction and provenance

Install `requirements.txt`. From the repository root, regenerate the
small exact checks and the summaries of archived observations with:

```bash
python3 -B contrib/v4_bulk_diagnosis_20261117/validate.py
python3 -B contrib/v4_bulk_diagnosis_20261117/controls.py
python3 -B contrib/v4_bulk_diagnosis_20261117/ring_ed.py
python3 -B contrib/v4_bulk_diagnosis_20261117/analyze.py
python3 -B contrib/v4_bulk_diagnosis_20261117/entropy_analysis.py
python3 -B contrib/v4_bulk_diagnosis_20261117/audit_results.py
```

Those commands regenerate their small report files in this contribution.
Rebuild production states in a separate directory to preserve the archive:

```bash
python3 -B contrib/v4_bulk_diagnosis_20261117/reproduce.py \
  --run-root /tmp/v4-reproduction --tag v4_L144_P1_s1_c320
```

The wrapper rebuilds parent stages sequentially and records the final
energy difference from the archive. Repeat `--tag` for selected cases;
omit it to run the completed configurations listed in
`completed_configs.txt`. Each process uses one BLAS thread.
The baseline L=96 states can be regenerated with the previous
contribution's `reproduce.py`, selecting `legacy_V4_a1_L96_P0` and P1,
then measured with this contribution's `measure_saved.py --parity P
--checkpoint PATH --run-root OUTPUT_DIRECTORY`. `measure_channel.py`
measures the secondary XX+YY channel of a specified saved checkpoint.
`auxiliary.py` sequences exact checks, measurements and an optional
benchmark in one worker slot; run it only when that slot is available.

Configuration, observation and source hashes, package versions, runtime,
peak memory and covariance/energy consistency checks are recorded in
`audit_results.json` and `source_audit.json`. State files remain ignored
in `runtime/`. New contractions are checked against full small-chain
matrices; the DMRG solver itself is shared with the prior campaign,
not an independent implementation. [Decisions](DECISIONS.md) records
the post-hoc entropy analysis, resource-accounting deviation and
translation-degeneracy handling for the small rings.
