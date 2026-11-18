# What the committed finite-entanglement ladders actually support

This is an independent **reanalysis of rounded numbers**, not a reproduction
of the iDMRG states. Input hashes, every parsed row, every fit window and
branch comparison are in `results/ladder_audit.json`; the compact fit table
is `results/ladder_fits.csv`. Entropies use natural logs and the fitted
quantity is c_slope=6 dS/dlog(xi), following the one-cut scaling convention
of [Calabrese--Cardy](https://arxiv.org/abs/hep-th/0405152) and finite
entanglement scaling of [Pollmann et al.](https://arxiv.org/abs/0812.2903).
It is a Virasoro-c estimator under the relevant scaling assumptions, not
the nonunitary c_eff=c-24h_min.

## Declared windows and results

We fit chi>=64,96,128,160 and the last four entries, requiring at least
three points and strictly increasing xi. Branches are never spliced to
manufacture a ladder. The following are raw slopes; the spread is fit-window
sensitivity, not a statistical confidence interval.

| h, branch | chi maximum | c, chi>=96 | c, chi>=128 | main issue |
| --- | ---: | ---: | ---: | --- |
| 1.738, pin | 200 | 1.065 | 1.055 | odd transfer cat at chi>=64 |
| 1.739, pinB | 320 | 2.158 | 2.357 | large nonmonotone local slopes |
| 1.740, cplx | 200 | 1.963 | 1.899 | competing bond/chirality branches |
| 1.7405, pinB | 320 | 1.919 | 1.929 | no equally long competing branch |
| 1.741, pin | 200 | 1.979 | 1.970 | order orientation changes at chi=200 |
| 1.744, pin | 200 | 1.734 | 1.725 | still substantially above 1.5 |
| 1.746, pin | 200 | 1.685 | 1.552 | significant early-point drift |
| 1.7465, pin | 256 | 1.589 | 1.512 | critical candidate |
| 1.747, pin | 256 | 1.525 | 1.445 | latest local slope is 1.402 |
| 1.749, pin | 200 | 1.345 | 1.264 | flow toward one-boson behavior |
| 1.750, pin | 200 | 1.254 | 1.195 | source's h=7/4 point |
| 1.755, pin | 200 | 1.093 | 1.085 | nearly stable chiral expectation |
| 1.758, pin | 200 | 1.109 | 1.102 | highest supplied target field |

The table extends through h=1.758, beyond the abbreviated source README.
At h=1.7465 the chi>=160 fit is about 1.50 and the last adjacent slope
about 1.51. At h=1.747 those values are lower. This supports examining
the lower part of the inherited h_I=1.7466 +/-0.0006 band, without
constituting an independent location or thermodynamic uncertainty bound.

The proposed order exponent is not a second independent proof. The
submission locates h_I by requiring the finite-chi logarithmic derivative
of |m| to approach -1/8, then quotes that same exponent as Ising evidence.
At h=1.7465 the full chi>=96 slope is -0.168 and the final adjacent slope
about -0.137; at h=1.747 they are -0.141 and -0.102. Window/field dependence
matters. We do not impose -1/8 in any new spectral calculation or fit an
independent eta from data that are absent.

## The apparent c about 2 is real as a finite-window observation

At h=1.7405, chi=96--320 spans xi=250.7--736.5. Adjacent entropy slopes
are about 1.91, 1.86, 1.93, 1.98, 1.91; the c=1.919 fit is quite smooth.
The reported B magnitude falls from 0.04458 to 0.02356 with fitted
dlog|B|/dlog(xi)=-0.591. This is **not** a demonstrated bond-order plateau.
It may be a finite-entanglement symmetry-breaking approximation to a
critical state; an asymptotic plateau below the available scale is not
excluded. A decreasing pinned one-point function does not independently
measure the connected bond exponent.

Replacing xi_max=xi1 by xi0 gives c=1.924 here. The ratio xi1/xi0 is almost
constant at about 4.18--4.20, so a simple accidental change of the longest
parity branch cannot remove this c about 2 signal. At h=1.7465 the same
ratio is about 4.20--4.22 and the two fitted c values differ by only 0.004.
Neither length has an operator label identifying a boson or an Ising field;
the ratio is not a velocity ratio. Raw eigenvalue phases and overlaps are
unavailable. A second boson with a different velocity remains an alternative
to crossover behavior, not a conclusion from two lengths.

Branch sensitivity is substantial in the lower window. At h=1.740,
chi=128, branches differing by only 2.80e-7 in energy density have |m|
between 0 and 0.02026 and |B| between 0.00039 and 0.03253. At chi=200
the energy spread is only 2.10e-8 while the order orientations still differ.
This is concrete evidence that energy agreement alone cannot identify a
stable ordered phase there. It does not prove that the c about 2 slope
is metastable: those branches have similar entropies and lengths.

The defensible alternatives remain a narrow extended two-mode phase,
an isolated higher-c point with a wide finite-scale crossover, and a
finite-entanglement/branch artifact. The saved data do not distinguish
them. There is no basis for assigning one precise boundary to that window.

## Controls and audit of convergence claims

Recomputing the stated chi>=96 fits gives:

| Supplied control | c_slope |
| --- | ---: |
| XX with gapped polarized spectator | 0.998 |
| XX x ordered TFIM, g/J=.8,.97,.98,.99 | 1.004, 1.016, 1.026, 1.071 |
| XX x critical TFIM, J/tA=.35,1,2.5 | 1.429, 1.442, 1.495 |
| **Interacting XXZ x critical TFIM, JzA=-.390** | **1.626** |
| XX x XX, equal hopping | 1.963 (three high-chi points) |

The interacting boson control is outside the source README's summarized
critical-control range 1.42--1.52. Its adjacent high-chi slope reaches
1.914 at chi=200. It cannot be used as a calibrated percent-level error
bar for the target. The tB=.35 unequal-velocity two-boson ladder ends at
chi=128, leaving only two points above 96; our declared fit is rejected.
The tA=0 "pure Ising" row has an exactly free A spectator at zero field
and nonmonotone lengths/energies; it is not a useful unique-vacuum control.

The source code's idmrg_conv stops when successive energy, entropy and
chiral expectation meet its thresholds, or when its block/time limit is
reached. It **does not test convergence of the competing B observable**.
The committed tables omit the convergence boolean, block history, actual
retained chi, truncation errors and full transfer spectrum. The saved
NPZ report dictionary also omits the convergence boolean. The driver can
save a final HDF5/NPZ and later skip that stage even when `ok=False`, unless
its separate near-job-end condition intervenes. The existence of a saved
row or DONE marker is therefore not itself a convergence certificate.

This does **not** establish that any particular supplied row failed: the
missing logs are needed to decide that. It limits what this audit can
certify. The tables alone cannot check residual stationarity, rule out
other unit cells, reconstruct connected correlators, measure velocities,
or exclude weak first-order behavior beyond the observed lengths.

Finite-chi Hellmann--Feynman energy secants against the reported Z are
retained as a normalization check. They support use of -h sum Z with a
single full on-site field, but do not add a new method for the critical
exponents or establish the global variational minimum.
