# Audit of the actual asymmetric finite-density records

## Point, provenance and convergence limits

All entries analyzed by `records.py` have the p153 bond and the same
`hz=-19.97589125`. They are full parity-only qubit states, not the
number-conserving projected model. Physical sites have spacing one and
momenta are radians per physical site. No field is refitted with length.

The script records each input path and SHA-256 hash in
[record_results.json](record_results.json). Relevant saved sets are:

| Set | Largest state represented | Correlations available here |
|---|---|---|
| C190 outG | L=256, chi=768, both parities | Central 48-site one-body matrices |
| C190 outX1024 | L=256, chi=1024, both parities | Same 48-site window |
| C190 outX1536 | L=256, chi=1536, odd parity | Same 48-site window |
| C195 out512 | L=512, chi=1536, both parities | Same 48-site window |
| C195 out_x2048 | L=512, chi=2048, even parity | Same 48-site window |
| C193 checkpoint measurements | L=64,128,256; largest chi=768 | Long one-body rows, connected-density inputs, d=2 and d=3 pair rows and their means |

Neither the high-chi L=256 nor the L=512 set supplies the needed long pair
correlators, density-current response tensor, dynamic poles, or matched
valley charges. The long C193 records are measurements of the earlier
chi=768 states, not a continuation to the later chi values. Their agreement
of local contractions with the parent measurements is not a chi test.

At L=512 the final even state has a small within-stage entropy change, but
the chi=1536 to 2048 increase in the central entropy is **0.01177684**, above
the recorded 0.01 diagnostic. The odd chi=1536 state's preceding increase
is **0.0868665**. At L=256 the odd chi=1024 to 1536 change is **0.00387464**.
Small discarded weight and local eigensolver residuals do not remove these
finite-entanglement or finite-length limitations.

C195's latest saved analysis has a ground-branch entropy-size estimate
near 1.954 and profile estimate near 1.989. These are evidence for a large
finite-scale entropy coefficient. They do not give two independently
normalized bosonic responses or establish the asymptotic Virasoro c.
The open-boundary parity gap cannot fill that gap: C191 shows that different
boundary fixed points exchange the relative density and winding terms in
the gap, and no appropriate asymptotic boundary condition is established
at this asymmetric point.

## New finite-window carrier calculation

For each saved central one-body matrix, average its complex r-th diagonal
to obtain `G_b(r)=<b_i^dagger b_(i+r)>` averaged over the available anchors.
For each trial carrier seed -2.9 and -0.24, demodulate by the seed, apply
four `[1/4,1/2,1/4]` convolutions, unwrap the remaining complex phase and
fit its slope. The filter radius is four sites; edge-contaminated distances
are removed. The windows (8,24), (8,32), (12,36) are exploratory diagnostic
choices, not registered asymptotic scaling windows.

Let `Q_hat=kA_hat+kB_hat`, with `Q_hat+pi` reduced near zero. Representative
results for window (8,32) are:

| L | Parity | chi | Raw central density | Q_hat + pi |
|---:|---|---:|---:|---:|
| 128 | even | 512 | 0.221109 | +0.011110 |
| 128 | odd | 512 | 0.212980 | +0.011361 |
| 256 | even | 1024 | 0.219250 | +0.009165 |
| 256 | odd | 1536 | 0.215345 | +0.009248 |
| 512 | even | 2048 | 0.215383 | +0.008303 |
| 512 | odd | 1536 | 0.217293 | +0.008306 |

At L=512 even, changing among all three windows gives
`Q_hat+pi=0.008303 to 0.008329`, and the corresponding diagnostic
`wrap(2Q_hat)=0.016606 to 0.016659`. At fixed L the displayed carrier sum
is much more chi-stable than the entropy: its shift from chi=1536 to 2048
is below 4e-7 at L=512. The shift between L=256 and 512 is about 9.5e-4,
and from 128 to 256 about 0.002. The same 48-site window is being used as
L grows; this is **not** a growing-distance thermodynamic fit.

For comparison, in units |t2|=1 the bare second-order dispersion is

\[
\varepsilon(k)=2\cos(k-1.53)-2\cos(2k),\quad
\mu_{\rm eff}=(40+2h_z-3/40)\,80=-2.1426.
\]

Its minima are -2.8918126002 and -0.2552137280, with
`Q_bare+pi=-0.00543367464` and `wrap(2Q_bare)=-0.01086734929`.
Thus substituting bare minima for the observed finite-density carrier sum
would even give the opposite sign. This comparison is a numerical statement
about the recorded correlators, **not** a proof that the thermodynamic
dressed quartet mismatch has the observed sign or magnitude.

## Independent connected-pair crosscheck

For `P_j=b_j b_(j+d)`, d=2,3, C193 stores the raw
`<P_i^dagger P_(i+r)>` and both means. Subtract
`<P_i^dagger><P_(i+r)>` before filtering around -pi. In this file format
the stored `mm` is already the creation-pair mean; it is **not** conjugated
a second time. Preserve complex phases throughout.

At L=256, chi=768, using four central anchors and windows (12,32) and
(16,40) gives:

| Parity | Range of connected-pair Q_hat + pi | Centered one-body carrier sum, longer windows |
|---|---:|---:|
| even | 0.00924043 to 0.00970692 | 0.00922396 to 0.00935052 |
| odd | 0.00933357 to 0.00977609 | 0.00935771 to 0.00944565 |

The right endpoint clearance is at least about two times r_max for the
(16,40) pair windows. This is useful geometry, but only one L/chi pair
provides these longer observables. Different anchors retain Friedel
oscillations; d=2 and d=3 may have different subleading overlaps. L=128
permits shorter checks but not a comparable long-distance window.

The combined result supports two finite-window carriers and a pair harmonic
near their sum. It establishes neither separate valley occupations nor two
independent massless modes. Peak amplitudes are operator residues, not
valley particle counts.

## Calibration and honest uncertainty

`calibrate.py` uses the exact pair correlator of two decoupled half-filled
XX legs, multiplied by a specified gauge phase, and adds a known
disconnected expectation plus a distant contaminating carrier. With true
`Q=-pi+0.0085`, the recovered phase errors are:

| Window | Carrier error |
|---|---:|
| 12--32 | 4.24e-5 |
| 16--40 | 2.52e-5 |
| 32--96 | 4.96e-6 |
| 48--160 | 1.93e-6 |

Adding a second harmonic only 0.012 radians away, with relative amplitude
0.45, makes the (16,40) bias **0.00370146** while its formal regression
error is only **4.61e-6**. This is an intentionally synthetic contamination
test of the estimator, not a statement about p153c's amplitudes. It proves
that residual-based phase errors alone are inadequate at the required
precision. The finite-size/window/anchor spread is reported as a spread;
it is not promoted to a confidence interval for the true momentum.

Consequently no p153c exponent is fitted here. In particular, neither the
old balanced `Kc,Ks` inversion nor the raw `Cnn/|q|` slope supplies the
asymmetric G. The range of possible matching scales a0 is also unresolved.
If one *assumed* the L=512 carrier sum were the matched quartet momentum,
`1/|2Q_hat|` would be about 60 sites. That assumption has not been validated,
so this number is not a locking exclusion or an input with a defensible
RG uncertainty bar.
