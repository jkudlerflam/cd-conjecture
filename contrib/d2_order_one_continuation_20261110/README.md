# Continued search for a simple nearest-neighbour qubit candidate

OpenAI Codex, GPT-6 family, 2026-11-10. This is a continuation of local
commit `b8b6179e`. Previous raw results are immutable inputs. All new
Hamiltonians are one-site translation invariant, nearest neighbour, and
act on physical qubits; longer-range free controls are labeled controls.

**A promising finite-size candidate was found in the exploratory R<=10
band.** The fixed Hamiltonian below has entropy size and profile estimates
near 1.45 through 192 sites. Both parities, an intermediate length,
independent initialization, higher bond dimensions, order diagnostics and
energy scaling support continued investigation. It is not an established
c>1 CFT or a counterexample. No numerical job remains active. Deferred grids
and the interrupted wider-ratio exploration are not completed results.

## Leading fixed Hamiltonian

With X,Y,Z denoting Pauli matrices and every bond joining j to j+1,

`H = sum_j [(XX-YY)/2 + (XY-YX)/10 + 5 ZZ/8] - 7 sum_j Z_j/8`.

Equivalently,

`40 H = sum_j [20(XX-YY) + 4(XY-YX) + 25 ZZ] - 35 sum_j Z_j`.

The model has four grouped terms, six expanded Pauli strings, and four
distinct nonzero coefficient magnitudes. Its raw coefficient spread is 8.75,
within the stated exploratory 3<R<=10 band, but outside the preferred R<=3
band. The exchange-only spread is 6.25. The bond interaction matrix has
singular values (0.625, 0.600, 0.400), spread 1.5625; the field norm divided by
the exchange Frobenius norm is 0.91693. The antisymmetric/symmetric exchange
Frobenius ratio is 0.14985. These invariant scales are reported alongside
the raw Pauli spread rather than used to disguise a small DM component.

The exact relation Jx+Jy=0 removes real hopping in the chosen fermion
representation. Jx-Jy=1 and the imaginary hopping scale is 2D=0.2. This is a
structural equality of the selected symmetric family, not evidence that
all effective low-energy scales are comparable. The field is the exact
rational -7/8, selected from a registered grid. Robustness over a parameter
interval and the ultimate tuning precision are not established.

## Evidence through 192 sites

Both exact parity sectors were optimized in the full microscopic model.
The lower variational-energy branch is even at 24, 48, 64, 96 and 192 sites. All selected
records pass the registered sweep stopping gate. This branch selection is
not a rigorous certification of the ground state.

The fixed estimators are `6 Delta S/log(L2/L1)`, with natural-log entropy
and open boundaries, using (i) the three central cuts and (ii) a
trapezoidal average over cuts [3L/8,5L/8].

| Sizes | Central-three | Fractional window |
|---|---:|---:|
| 24→48 | 1.6826 | 1.7038 |
| 48→64 | 1.4210 | 1.4474 |
| 48→96 | 1.4185 | 1.4363 |
| 64→96 | 1.4167 | 1.4284 |
| 96→192 | 1.4475 | 1.4530 |

The completed 192-site profile fits give 1.45393 and 1.46296 in the two
registered fit windows, with RMS residuals 0.00201 and 0.00181. Exact free XX
controls give a fractional 96→192 slope of 0.99224, and two independent XX chains
give 1.96969. These controls calibrate the estimator, not a universal bound
on finite-size corrections in the interacting model.

The 96-site chi192→256 change is 8.11e-6 in fractional entropy and 8.70e-6
in central-three entropy. At 192 sites the corresponding changes are
0.0002470 and 0.0002606, moving the final entropy coefficients by only
0.00214 and 0.00226. An independent L48 seed23 calculation at chi192
matches the original seed17 result; the chi128→192 entropy change is below
1e-6. Thus the observed entropy excess is much larger than these numerical
changes. The final even-sector energy variance per site is 1.91e-9 at L192.
Global energy variances and their limits are retained in
`variances.json`; variances alone do not certify a ground state or the sign
of a small parity splitting.

Independent saved-state contractions give central-window transverse order
measure m²=0.44742, 0.36226, 0.28904, 0.22307 at 24, 48, 96, 192 sites.
The far X correlation RMS decreases 0.36169→0.27654→0.20322→0.14257.
Density modulation also decreases: central detrended RMS
0.00910→0.00650→0.00160→0.000620. No order plateau
has been established in these data. Odd-parity results behave similarly.

The exact decomposition of entropy into binary left-block parity entropy
and within-parity entropy assigns only 0.00173 of the 48→96 fractional slope
to changing parity weights; the 96→192 contribution is -0.00470. It excludes that narrow explanation of the
excess, but does not exclude every possible ordered cat. The independent
correlation measurements are essential. The variational parity splitting
has L times its magnitude approximately 0.00661, 0.00607, 0.00642, 0.00721
across the four doubling sizes, broadly compatible with 1/L over this range.
It is not a velocity measurement or proof of a CFT, and the energy standard
deviation exceeds this small splitting. Both fixed-parity entropy sequences
retain excess growth, so the result does not rest on certifying its sign.

A post-hoc energy combination `D(L)=E(4L)-3E(2L)+2E(L)` cancels bulk and
boundary constants. Its effective correction power is 1.0361, versus 0.9497
for the exact XX control on the same sizes. A three-size `eL+es+a/L` fit
predicts E192 within 2.18e-5. This supports consistency with inverse-size
energy corrections; it determines neither c, velocities, nor a dynamical
exponent. See `diagnostics/energy_consistency.json`.

## What is unresolved

The initial entropy coefficient decreases from 24→48 to 48→96, then remains
near 1.45 over 96→192. The L64 check supports a smooth trend between 48 and
96. Longer-scale crossover remains possible. All couplings stay fixed across sizes.

Renyi alpha2 and3 fractional slopes over 96→192 are 1.4809 and 1.3702,
compared with the von Neumann estimate 1.4530. The exact controls also show substantial
alpha-dependent finite-size corrections, so this is a warning rather than
an automatic rejection. A common asymptotic coefficient is not established.

Full-Hilbert-space periodic ED at L14,16,20 has profile slopes close to 1.03
and a 16→20 half-chain size slope close to 1.04. These small rings do not
independently confirm a second mode. They also cannot exclude a longer
crossover scale or a dilute slow sector. See `RINGS.md`; raw and
translation-resolved entropies and gap semantics are preserved. A larger
periodic or infinite-chain calculation and bulk low-energy spectra/velocities
are important next discriminants. Boundary effects, a c=1 crossover and
nonrelativistic critical behavior remain alternatives to an identified CFT.

An exact nearby D=0 weighted-Dicke control has an OBC-normalized effective
entropy coefficient approaching 3 without a relativistic CFT, despite clean individual profiles
and saturated binary parity entropy. Its open-parent end fields differ from
the candidate convention. It demonstrates why entropy alone is insufficient;
it does not classify this nonzero-D candidate. The control and its precise
boundary convention are documented in `diagnostics/README.md`.

The preferred attractive points with coefficient spreads 2.656 and 2.8
reproduce 48→96 slopes near 1.42, but about 31% of their growth comes from
changing parity weights. The reference falls to 1.142 at 96→192, with
sub-one profile fits and increasing endpoint X correlations; the rational
192-site comparison is deferred. The R=3 reserve shows a strong longitudinal-order
warning. The simple integer model `(7XX-YY+5DM-7ZZ)/12` instead gives a
48→96 fractional slope 0.982 and is not promoted. The older R=28 bridge
remains separate from the comparable-coupling band. Its 96-site chi256
checks completed; its optional 192-site runs were interrupted under the
registered focused stopping rule after this candidate passed validation.
Their raw records and last checkpoints are preserved, and no incomplete
192-site value enters a size claim.

## Reproduction and audit trail

- `PROTOCOL.md` and each branch's protocol/amendments distinguish initial
  registrations, adaptive follow-ups, completed runs and deferred work.
- `bridge_scan/REPORT.md`, `attractive_scan/RESULTS.md`, and
  `simple_attractive/PROTOCOL.md` document the separate routes.
- `SIZE_COMPARISONS.md` and `size_comparisons.json` are generated by
  `analyze.py`; every selected record has a source path and SHA256 hash.
- `run_inventory.json` lists current completed/running records and preserved
  lower-chi snapshots. `progress.py` gives compact live status.
- `audit.py --final` verifies quiescent records, solver hashes, geometry,
  convergence, entropy/energy consistency and checkpoint metadata.
- `diagnostics/README.md` explains the exact saved-state contractions,
  symmetry caveats, archived lead audit, and factorization-line controls.
- `renyi.py`, `variance.py`, and `analyze_rings.py` regenerate their
  corresponding summaries. Completed checkpoint snapshots are verified by
  hashes before and after expensive contractions.
- `archive_candidate.py` made 12 hash-verified winning-candidate checkpoint
  copies in `checkpoints_local/` (about 613 MB), including the intermediate
  size, independent seed and pre-upgrade L192 state. These binaries are
  Git-ignored; `local_checkpoint_manifest.json` records their paths/hashes.
  Two unfinished R=28 states are separately preserved there and documented
  in `interrupted_exploration.json`.

One preservation limitation is explicit: the pre-upgrade L96 even chi192
MPS was overwritten before a binary snapshot was made. Its raw JSON and
previously measured, hash-tagged spectra/diagnostics remain; its exact old
MPS is unavailable. The final L96 chi256 state and the L192 pre/post-upgrade
states are archived. See `bridge_scan/STATE_PRESERVATION.md`.

The prior `dmrg_screen.py` wrapper and C190 solver are reused unchanged.
The local scientific environment is `/private/tmp/cd-conjecture-venv`;
each numerical process limits BLAS to one thread. Every result records its
exact checkpoint path. Root-run checkpoints use
`/private/tmp/d2-order-one-continuation-checkpoints`; the bridge and attractive
workers use `/private/tmp/d2-order-one-continuation-bridge` and
`/private/tmp/d2-order-one-attractive-checkpoints`, respectively.
The previous checkpoint collection is not modified. New results are local;
no cluster submission or remote publication is claimed.
