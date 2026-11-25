# IAS quartet pilot and estimator audit

2026-11-25. **[Numerical; ongoing.]** The historical six-file quartet
table reproduces exactly, but its quoted `x4 approximately 2.05 +/- .2`
does **not** include a substantial local-filter bias. That estimate is
not a calibrated charge-four dimension. The new data do not establish
quartet locking or asymptotic stability of the qubit candidate.

## Claim and hypotheses

The target is the Hermitian, one-site translation invariant nearest-
neighbour physical-qubit (`d=2`) Hamiltonian

    H = sum_j [(19/32) X_j Y_(j+1) + (13/32) Y_j X_(j+1)
               + Z_j Z_(j+1)] - (1783/1000) sum_j Z_j.

Every open-chain site receives one full field. Delta=3/16 and h=1.783
are fixed. The real computational frame is the site-dependent z rotation
theta_j=j*pi/2; it shifts pair momentum by pi and quartet momentum by zero.
Only spin parity is exact. A lattice quartet monomial need not be an
emergent charge-four primary. The desired CFT is unitary and compact and
its central charge is Virasoro c; finite entropy slopes do not identify it.

## Saved-data audit

The original compressed data and scripts are unchanged. The exact
[reproduced table](analysis/historical_table_reproduced.txt) is followed
by 12,672 deterministic variations in [sensitivity.csv](analysis/sensitivity.csv):
window widths 24,32,40,48,56,64, harmonic orders 1,2,3, each reference
position separately or their mean, and five distance windows. Further
[calibrations](analysis/calibration_sensitivity.json) vary the fitted
momentum by up to .002*pi and allow linear/quadratic harmonic envelopes.
Overlapping windows are **not independent statistical samples**. Fits to
log magnitudes are suppressed on sign-changing data; signed residuals
are retained. The historical table itself used log absolute values even
for some sign-changing control data.

The local Hann filter averages a rapidly varying amplitude over 48
sites before assigning it to the centre. Thus fitting the output to
an unfiltered power law is biased, even at zero numerical noise.

| Signal, with the historical 48-site filter and r=32–56 | True p | Naive output p | Forward-filtered power fit p |
|---|---:|---:|---:|
| Pure power plus rank-one neutral part | 3.68 | 4.8459 | 3.6800 |
| Pure power plus rank-one neutral part | 4.00 | 5.3763 | 4.0000 |
| Exact infinite XX x XX quartet | 4.00 asymptotically | 5.3812 | 4.0040 |
| Exact OBC XX x XX, L=256 | 4.00 asymptotically | 5.5913 | 4.0880 |
| Exact OBC XX x XX, L=384 | 4.00 asymptotically | 5.4684 | 4.0400 |

Here `p=2x` only after a scalar operator is matched. The exact two-boson
control consists of two independent XX chains per site, with microscopic
**d=4**, hopping tA=1,tB=.7 and zero fields (unequal velocities). For
`Q_i=(Sm_A_i Sm_A_(i+1))(Sm_B_i Sm_B_(i+1))`, Wick's theorem gives
`G_Q(r)=[C(r)^2-C(r-1)C(r+1)]^2`, where
`C(r)=sin(pi*r/2)/(pi*r)` in the infinite half-filled chain. Finite OBC
values use the exact occupied sine-wave projector. This is a calibration,
not a d=2 construction. The h=1.772 qubit control is retained from both
original chain lengths.

Applying only a pure-power forward correction to the high-chi L=384
target gives p=3.067 for Q(2,5) and 3.481 for Q(3,6). These are conditional
curve parameters, **not replacement x4 measurements**: a quadratic
harmonic-envelope filter moves Q(3,6) to p=4.035, with sign changes and
larger residuals. Shape, envelope, reference, window and chi dependence
are not controlled to .05 in x. Simply increasing chi cannot repair this
identification problem.

Synthetic tests include data-scale oscillations and multiple independent
neutral components. Two neutral r^-4 components can exactly reproduce
the isotropic two-quadrature covariance of a charged r^-4 vertex. This
is an explicit nonidentifiability example, not evidence that the target
actually has those neutral components. The rank-one neutral assumption
must be tested with enlarged operator measurements and further matching.

## Pilot and validation

The authorized IAS run is `quartet-pilot-20261125a`, Slurm **17764244**,
4 CPUs, 16 GB, one hour, L=128, even parity, chi=64,128,256,384,512.
At this checkpoint it is running. No unavailable historical MPS is assumed.
Each stage saves an atomic HDF5 checkpoint in the remote `states/` directory,
outside fetched `results/`; earlier stages are preserved. The new engine
uses the original `ofd.py` model, initializer and solver settings, with
one-sweep statistics and a conservative deadline. Resuming rebuilds the
environment from the saved MPS and resets mixer/sweep counters; it is a
variational warm restart, not a bitwise engine continuation.

- Literal original/rotated MPO agreement and full end fields: below 2e-15.
- Random complex L=12 parity-state correlators: below 3e-16 (shortened
  quartet shapes make a dense-vector test possible).
- Exact production quartet offsets, L=24 complex parity cat: 704 matrix
  entries agree with independent local-matrix contractions within 7e-16.
- HDF5 roundtrip: vector error below 4e-16, energy unchanged.

The enlarged raw matrix keeps Q(2,5), Q(3,6), their adjoints, Z, ZZ and
both nearest-neighbour hopping orientations. Hermitian quadratures can
be formed afterwards. Candidate-neutral probes are labelled as candidates;
their continuum charge is not assumed. Finite-distance correlation matrices
are not necessarily positive semidefinite Gram matrices.

## Reproducing

Requirements: numpy, scipy, h5py, physics-tenpy==1.1.1, threadpoolctl.
The verified IAS environment is documented in the [runbook](../../tools/typhon/README.md).
Saved-data analysis ran in an isolated temporary Mac environment; sustained
optimization runs only on IAS. The initial six-file audit took 6.8 seconds.

From the repository root:

```sh
python contrib/quartet_typhon_20261125/audit_saved.py
python contrib/quartet_typhon_20261125/audit_calibration.py
python contrib/quartet_typhon_20261125/validate.py --output /tmp/quartet-validation.json
python contrib/quartet_typhon_20261125/validate_production_ops.py --output /tmp/quartet-ops.json
python3 tools/typhon/typhon.py submit FRESH-RUN-NAME contrib/quartet_typhon_20261125/pilot.job.json
python3 tools/typhon/typhon.py status FRESH-RUN-NAME
python3 tools/typhon/typhon.py logs FRESH-RUN-NAME
python3 tools/typhon/typhon.py fetch FRESH-RUN-NAME
```

See [PROGRESS.md](PROGRESS.md), [RUNS.md](RUNS.md) and
[VALIDATION.md](VALIDATION.md) for current status, accounting, source
hashes and convergence. Allocation completion is not scientific convergence.

## Authorship and scope

OpenAI Codex, **gpt-6-astra, maximum reasoning**, verified from this chat's
model metadata; chat `01a21370-d364-7c60-a006-ee9c4074ddf3`. The checks
above are self-checks with distinct implementations, not independent human
or agent review. No subagents. This qualifies REPORT §9.217 without rewriting
its historical account. A matched x4<2 still needs a nonzero target coupling
to imply locking; x4>2 addresses only that weak-coupling perturbation.
