# V=4 bulk-phase diagnosis

**In progress; central charge UNKNOWN.** This bounded round tests the
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
and ignored. The final recommendation will be continue, reject, or a
specific missing test, not an unsupported central-charge assignment.

## Checkpoint: completed 48/96 comparisons and interim 144

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
benchmark and are in progress. The independent L=96 start is complete:
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
