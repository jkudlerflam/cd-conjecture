# D3. Targeted symmetric nearest-neighbor search

No verified counterexample was found. The search produces many finite-size entropy slopes above 2 or 3, and some energy-based estimates far above 3, but none establishes a collection of unitary CFTs with total charge above the conjectured bound. Several apparent violations can be located and explained. The small-angle Sp(4) region and some spin-3/2 and U(1) points remain unresolved.

This is a finite-size catalogue, not an exhaustive thermodynamic phase diagram. In particular, the requirement to identify every critical point and its central charge is not achieved. The available sizes and finite parameter grids do not support that claim. No inference of the conjecture from a negative search is made.

The main deliverables are [d3.py](d3.py), the full [numerical catalogue](catalogue.csv), [calibrations](calibrations_summary.json), the [derivations and assumptions](METHODS.md), and the [survey figure](survey.png). Reproduction and cluster instructions are in [README.md](README.md). All computations used periodic boundary conditions and every momentum up to the exact k versus -k equivalence.

## Scope and numerical checks

The local campaign finished in approximately 38 minutes, below the 75-minute cap. No cluster job was submitted. Exact timestamps and file hashes are recorded in `provenance.json`.

The retained catalogue has 3,526 finite-size records. These are records, not distinct Hamiltonians. They include repeated symmetry realizations and 768 energy-only points. The principal coverage is as follows.

| Family | Parameters | Sizes |
|---|---|---|
| Sp(4) | 96 angles over the full circle, then 768 energy-only angles at L=8; 23 selected angles | Full coarse scan L=6,8,10; selected L=4,8,10,12 |
| SU(2), spin 3/2 | 192 sphere points; 48 perturbations around SU(4), plus ULS, TB, and Heisenberg anchors | L=6,8,10; 28 selected sphere points at L=11 |
| Spin-orbital | 182 signed product points with x<=y; decoupled controls | L=4,6,8,10; 13 selected points at L=12, plus a velocity-ratio-3 control |
| Spin 1 | 96 bilinear-biquadratic angles and 24 deformations of the c=2 region | L=6,9,12; 14 selected angles additionally at L=8,10,14 |
| Spin-1 U(1) grid | 150 points in (Delta,beta,D), with all magnetizations checked | L=8,10; 10 selected points at L=12,14 |
| Saturation follow-up | 16 points, beta=1, Delta=0,1/2,1,3/2 and D=3,3.5,4,4.5 | L=8,10,12 |

The d=4 L=12 calculations use additional symmetry reduction. The generic spin-3/2 sphere reaches L=11. The Sp(4) L=12 ground-state calculation assumes that the zero-Cartan sector continues to contain a global ground state. Its agreement with the sufficient Sz sector is checked at L=8,10, to less than 4e-10 including the small rounding differences in angles. That continuation assumption is not a theorem supplied by this calculation.

Direct L=4 matrix and full-spectrum tests give errors below 8e-15. The independent ULS realizations agree in energy and, where the singlet is unique, in entropy. The largest recorded production eigenvector residual is 1.91e-8. These checks do not certify that an iterative eigensolver found every low eigenvalue, but no failed production solve remains in the retained data.

Two implementation failures were found and corrected. Odd-L spin 1 initially used Sz=1 rather than Sz=0. Those 120 records were quarantined in `invalid_spin1_L9.jsonl`, recomputed, and excluded from analysis. Empty momentum sectors initially caused a concatenation exception in the all-charge audit. The empty-sector and zero-component cases are now handled. There are also 25 physical replacements of constrained-sector states by lower-energy states in other magnetization sectors. These are not numerical failures.

## Calibration

The energy estimate assumes a single velocity and a dimension-one adjacent-momentum current. Its error below is the half-width between a two-size linear and a three-size quadratic Casimir fit, with corresponding velocity extrapolations. It is a fit-choice envelope, not a confidence interval. Entropy entries are raw finite-size slopes with the maximum shift under changing the first fitted cut from 2 to 1 or 3. Those errors are not bounds on thermodynamic c.

| Calibration | Supplied c | Largest L | Energy c, fit envelope | Raw entropy c at largest L | Soft primary x, observed / expected |
|---|---:|---:|---:|---:|---:|
| SU(3) ULS | 2 | 12 | 2.015 +/- 0.036 | 2.239 +/- 0.060 | 0.6421 / 0.6667 |
| SU(4) ULS | 3 | 12 | 3.038 +/- 0.072 | 3.511 +/- 0.098 | 0.7466 / 0.7500 |
| TB spin 1 | 1.5 | 14 | 1.495 +/- 0.020 | 1.580 +/- 0.070 | 0.3326 / 0.3750 |
| TB spin 3/2 | 1.8 | 10 | 1.811 +/- 0.055 | 1.956 +/- 0.117 | 0.2788 / 0.3000 |

The known charges all lie in these energy-fit envelopes. This is a successful calibration of the chosen procedure at these four points, not a general error theorem. The tower identifications use the supplied WZW descriptions and their primary dimensions as assumptions. The full low spectra are retained, but an exhaustive descendant/multiplicity identification is not claimed.

The entropy regression error at SU(4), L=12 is only 0.014, despite the 0.511 error in c. A good straight-line entropy fit is therefore insufficient evidence of a violation. SU(4) lengths 6 and 10 are not vacuum lengths and were excluded from its Casimir fit. For example, the inferred adjacent-gap velocity is only 0.0420 at L=6, compared with 1.5614 at L=12.

## (a) Sp(4)

The local Hamiltonian is h=cos(theta)P_5+sin(theta)P_10, with the singlet coefficient zero. No Sp(4)_2 tower was identified.

| Location or region | Evidence | Charge assessment |
|---|---|---|
| theta=pi/2 | h=(I+swap)/2 | Supplied SU(4)_1 anchor, c=3 |
| theta=2pi/3 and pi | Adjoint and vector current-momentum states, split velocities, Casimir consistency test below | Consistent with Sp(4)_1 times Ising, c_T=3; continuum identification assumed |
| Other sampled points from roughly 5pi/12 toward 7pi/6 | Commensurate L=8,12 entropy slopes drift toward about 3; strong velocity splitting near the upper end | Critical-compatible finite-size window, boundaries and eventual gaps not established |
| Around pi/12 through 3pi/8 | Entropy decreases with size and low 0/pi doublets appear | Dimer-compatible, no reliable nonzero c assigned |
| Small positive theta, particularly pi/48 | Very large entropy slopes, shrinking adjacent-gap velocity, moving 0/pi crossing | Unresolved slow-mode regime, not an established CFT |
| theta in [5pi/4,2pi], and theta=0 | Exact polarized/symmetric ground states or increased degeneracy; quadratic branch | Exact non-CFT exclusion of the apparent entropy excesses |

The interval described as critical-compatible is not a claim of a thermodynamic phase. A weak gap below the SU(4) point, in particular, is not excluded by L<=12.

At theta=pi/48 the raw entropy slopes at L=8,10,12 are 8.150, 7.992, 8.029. The corresponding adjacent-gap velocities are 0.05608, 0.04103, 0.03114. The effective gap exponents are z=2.400 and 2.513 on successive size pairs, while a CFT requires z=1. The L=12 entropy fits range from 6.370 to 9.167 as the cut window changes. Their raw regression error does not account for this discrepancy. Half-chain entropy differences give effective charges 6.32 and 6.92, also inconsistent with the spatial-slope estimate. This is a located failure of the CFT estimators at the accessible sizes, not a proof of the eventual infrared behavior.

The first crossing of the k=0 and k=pi energies moves substantially.

| L | Crossing theta |
|---:|---:|
| 6 | 0.235160201 |
| 8 | 0.171388684 |
| 10 | 0.138280852 |
| 12 | 0.114416470 |

The numerical root tolerance is 1e-9, but that precision has no bearing on the thermodynamic limit. The approximately 1/L drift argues against calling any of these a fixed critical coupling. At theta=pi/16, where L=10 gives a tempting entropy slope 3.841, the L=12 slope falls to 3.453 and the opposite-momentum doublet splitting is only about 0.0085 times the adjacent gap. There is no consistent Sp(4)_2 identification there.

The split-velocity test gives a more concrete explanation of the energy excesses away from the small-angle region. The exact lattice Casimir labels of the relevant states are C_2=3 for the adjoint 10 and C_2=2 for the vector 5. Assume these become the adjoint current and mixed bilinear of five Majoranas with velocity v_5 and one with velocity v_1. Then c_T=5/2+1/2=3, v_mix=(v_5+v_1)/2, and the predicted Casimir coefficient is A=pi(2v_5+v_mix)/6.

| theta | Extrapolated v_5 | Extrapolated v_1 | A from ground energy | A from two velocities | Naive single-velocity c |
|---|---:|---:|---:|---:|---:|
| 2pi/3 | 0.88449 | 1.67248 | 1.58490 | 1.59565 | 3.42223 |
| pi | 0.42243 | 1.88467 | 1.04388 | 1.04636 | 4.71953 |

The agreements are 0.67 and 0.24 percent. This is evidence for the proposed c_T=3 description, not a proof of that continuum theory. The independent entropy extrapolation from L=8,12, with 1/L as the leading correction, gives 3.025 +/- 0.449 and 2.977 +/- 1.140. Here the error is the maximum change under the cut windows and a 1/L^2 extrapolation. The broad errors are intentional.

## (b) SU(2), spin 3/2

The sphere scan covers the normalized coefficient vector perpendicular to (1,1,1,1). The ULS and TB anchors reproduce the calibrations. The region a_3=min_J a_J has an exact polarized ground state and a quadratic one-magnon branch. It cannot furnish the required CFT counterexample.

The strongest remaining sampled excesses are not accompanied by stable towers.

| Sphere point | a_0,a_1,a_2,a_3 | Raw c at L=10 | Observation |
|---|---|---:|---|
| 70 | -0.368683,-0.262743,0.861465,-0.230038 | 7.698 +/- 1.670 | Total spin 0; adjacent gap has z_eff=2.826 from L=8 to 10 |
| 110 | 0.692909,-0.697790,-0.125925,0.130806 | 5.243 +/- 0.813 | Total spin 0; very small velocity 0.02687 at L=10 |
| 160 | 0.399192,-0.388422,-0.592632,0.581861 | 5.117 +/- 0.771 | Total spin 2; changing ground momenta and strong size dependence |

The errors in this table are cut-window changes, not errors on c_T. Point 70 falls to a slope 3.826 at odd L=11, where the vacuum sector differs. That comparison cannot be used as a controlled extrapolation. The remaining sphere points are catalogued individually, but no complete phase boundary map or new CFT identification is established. This is a substantive limitation of the attack.

## (c) Spin-orbital

The scanned product Hamiltonian differs from sign*(S dot S'+x)*(T dot T'+y) only by a constant. Both signs and the x/y exchange symmetry are included.

| Point or region | Numerical result | Interpretation |
|---|---|---|
| x=y=1/4, positive sign | SU(4) calibration | c=3, supplied anchor |
| x=0,y=1/4 | Raw c falls 4.060 at L=8 to 3.723 at L=12 | Entropy extrapolation 3.048 +/- 0.526, compatible with 3 |
| x=1/8,y=1/2 | Raw c falls 3.905 to 3.611 | Entropy extrapolation 3.024 +/- 0.843, compatible with 3 |
| x=y=3/8 | Raw c falls 3.741 to 3.480 | Entropy extrapolation 2.959 +/- 0.313, compatible with 3 |
| x=y=1 | Raw c falls 1.640 at L=8 to 1.289 at L=12; adjacent-gap velocity grows | Gapped-compatible, no nonzero c assigned |
| x=y=3/4 | Exact alternating spin/orbital singlet ground states; ED gives E0=-9L/16 in the stored convention | Dimerized ground states proved; a bulk gap is not proved here |
| Local triplet-triplet energy is minimal, including several degenerate product lines | Polarized product attains the local lower bound | Exact polarized or degenerate exclusions; see METHODS.md |

The extrapolation errors use the same cut and correction-model variation as in the Sp(4) discussion. They are too large to resolve a small excess above 3.

A deliberate control uses h=A+3B, with A and B the two AF spin-1/2 Heisenberg bonds. Under the usual SU(2)_1 identification of each factor, c_T=2. The ground state is a tensor product and is independent of their positive coupling ratio. The entropy slope at L=12 is 2.092, while the single-velocity Casimir procedure gives 4.010 +/- 0.055. This false positive persists in the thermodynamic formula because it measures c_1+3c_2, not c_1+c_2. No larger-c inference based only on that energy estimator is valid.

## (d) Spin 1

The supplied extended c=2 phase has the expected upward finite-size bias. At theta=pi/3 and 5pi/12, L=12 gives raw c=2.254 +/- 0.062 and 2.303 +/- 0.072. These are not violations of the supplied continuum identification. The TB point gives the c=1.5 calibration. The interval [pi/2,5pi/4] has an exact polarized ground state and a quadratic branch, with additional degeneracy at the endpoints.

The all-charge audit matters for the SU(3)-breaking deformations. At theta=5pi/12, extra SzSz coefficient -1 and D=1, the L=12 M=0 state has raw c=5.907. Its energy is 15.4583. The actual ground state is in M=10, with energy 14.6600 and raw c=1.168. The former cannot be used in a ground-state CFT test.

Several U(1) grid excesses remain finite-size candidates.

| (Delta,beta,D) | Largest L | Raw c and cut variation | Qualification |
|---|---:|---:|---|
| (0,1,2) | 14 | 2.508 +/- 0.111 in M=2 | A degenerate M=6 ground state gives c=1.038; no state-independent CFT inference |
| (1,1,2) | 14 | 2.748 +/- 0.213 | Conserved nonzero-color occupation is only K=6; filling changes and short conformal window remain unresolved |
| (3/2,2,1) | 14 | 2.223 +/- 0.059 | Similar scale of excess to the ULS calibration; no independent larger-c tower identified |
| (3/2,2,2) | 14 | 2.160 +/- 0.081 | Same limitation |

At beta=1 the extra color conservation was checked explicitly. For the first two rows K=6 has variance below 4e-27 in the selected L=14 eigenvectors. This explains why fixed total magnetization alone need not select a representative vacuum sector, but it does not prove a continuum charge bound.

For beta=1 and D>=max(4,1+|Delta|), the all-zero product ground state follows from a local positive-operator bound. At D=4 a single-particle branch is quadratic around pi. The saturation scan verifies these energies and does not produce a relativistic multicritical counterexample. At strict inequality the same local bound gives a nonzero many-body gap. The numerical catalogue does not identify every other multicritical surface in the three-parameter U(1) family.

## Verdict and continuation

There is no demonstrated c>2 at d=3 or c>3 at d=4. There is also no numerical proof that such a point is absent. The strongest retained unresolved cases are (i) the small-angle Sp(4) regime and its SU(2)-breaking neighbors, (ii) the nonpolarized high-entropy SU(2) sphere points, and (iii) selected dilute-color U(1) states. The near-SU(4) spin-orbital results are compatible with c=3 within the available extrapolation uncertainty.

The prepared 154-task SLURM manifest targets these cases and larger compatible calibration lengths. It has not been submitted. Larger L, controlled vacuum/filling sequences, and independent velocity/tower identification are needed before promoting any raw excess to a counterexample. A complete thermodynamic phase catalogue remains open.
