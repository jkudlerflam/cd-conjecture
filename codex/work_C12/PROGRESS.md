# C12 progress

The task is to find a simple one-site translation-invariant nearest-neighbor Hermitian chain with total central charge c > d - 1, preferably d <= 64. This run is closed book. No pre-existing files are read or listed, and no external sources are consulted. Work and generated files stay in the current directory. The supplied facts are taken as input.

## Initial strategy

Screen frustration-free cells, squeezed-space models, emergent constraints, and spontaneously formed ladders before committing to numerical work. A successful construction must have an explicit two-site interaction, a controlled relativistic infrared theory, and two calibrated central-charge estimators. Degeneracy and a quadratic ferromagnetic mode are not central charge. Negative results will be reported as such.

The request to remain closed book supersedes the usual project-memory reading and updating workflow. PROGRESS.md and RESULT.md will carry this run's record.

## First screening

1. A frustration-free ferromagnetic SU(q) cell has symmetric ground space Sym^b(C^q), dimension binomial(b+q-1,q-1). Projection of an end generator gives the collective generator divided by b. This produces a large representation, not b independent flavors. An SU(q)_b WZW point would be attractive because c=b(q^2-1)/(b+q), but the projected first-order bond is only bilinear and does not supply the integrable higher-degree interaction. This is a possible target, not a constructed counterexample.
2. A fixed injective MPS parent has a cell ground space bounded by D^2 once the length exceeds its injectivity scale. Its two edge spaces do not grow with cell length. Colored-path constraints can have growing degeneracy, but boundary access is through stack operations, not independent gapless channels.
3. A rigid rung carrying k independent q-state chains costs d=q^k and gives c=k(q-1), which never exceeds d-1. Unrolling the rung into a long cell moves the problem to transmitting k independent channels through its single boundary bond.
4. Large entropy from an SU(q) ferromagnet, extensively degenerate constraints, or a quadratic band edge cannot be counted as a relativistic central charge. These need explicit z tests.

No candidate has yet passed the screen. The main remaining analytic possibility is an inexpensive large-representation cell tuned to a higher-level WZW point. Its missing ingredients are period selection and independently tunable projected interactions.

## A direct dilute-gas route at d=3

Consider local occupation states 0,1,2, with a conserved total occupation N. A general small family has atomic hopping 10 <-> 01, pair hopping 20 <-> 02, and pair conversion 11 <-> 20,02. Three stable incommensurate minima, one atomic and two molecular, would potentially give three Luttinger modes. Ordinary two-atom conversion is oscillatory when molecular minima occur at nonzero momenta; the momentum-neutral four-atom conversion must also be tested by its scaling dimension. Stability against larger bound clusters and phase separation is essential.

Created pair_screen.py, an exact two-particle relative-coordinate diagonalization. In 3,000 random parameter samples for equal-sign, reflection-even conversion, no interior pair-band global minimum below the two-atom threshold was found. The relative coordinate is truncated at 64 sites, and K is sampled at 65 points. This is a negative numerical screen, not a general exclusion theorem. The first launch failed because this environment provides python3 rather than python; the rerun completed.

An odd-parity conversion (opposite signs for 11 -> 20 and 11 -> 02) changes the pair-channel form factor from cos(K/2) to sin(K/2). It is Hermitian and one-site translation invariant. A second two-body scan is running for this case.

## Explicit two-body candidate

Odd conversion produced 282 interior bound-band minima below the atomic continuum threshold in 5,000 random samples. A simpler analytic threshold family exists. With atomic hopping t=1, V11=0, conversion amplitudes h[20,11]=-g and h[02,11]=+g, pair hopping h[20,02]=-tp, and onsite pair energy U, choose s=sin(Q/2),

    tp = g^2 (2+s) / [4(1+s)^2],
    U = -4 + 2 tp - 4 tp s^2 + 2 g^2 s^2/(1+s).

Then the molecular minima at ±Q touch E=-4, twice the atomic minimum -2. At g=2, s=1/2 this gives tp=10/9, U=-14/9, Q=pi/3. Relative-coordinate diagonalization at cutoff 160 gives Q=1.04719755088, E=-4.00000000000, curvature 0.46666679. The atomic mass is 1/2; the molecular mass is 15/7 to the observed precision. These are quadratic vacuum thresholds, not yet a z=1 phase.

Added qutrit_ed.py with the explicit 9 by 9 Hermitian bond. V12=V22=20 penalizes adjacent atoms and pairs and adjacent pairs; t21=0. Few-body periodic ED through N=6 and L=24 (N=6 reached L=18 so far) found no three- through six-particle energy below -2N. Positive finite-size shifts decrease with L. This only screens small-cluster collapse. It does not establish miscibility or central charge. Odd-even energy shifts suggest that the threshold choice favors pairs at finite density, so U must be detuned before counting an independent atomic mode.

## Finite-density screening and numerical safeguards

At L=16, N=6, varying U by delta=0 through 1.8 produces raw periodic-chain entropy slopes above 2. The largest short-chain value is around delta=0.9, with c_vN=4.19 and c_Renyi2=4.49. These are not central-charge measurements yet. The ground state changes momentum, many low levels are close, and at several couplings the two estimators disagree. A real Lanczos vector can also be a standing wave of degenerate ±K states. Added momentum_eigenstate() to select one momentum and check its translation residual before fitting.

Created simple_dmrg.py, a standalone NumPy/SciPy two-site finite DMRG implementation for nearest-neighbor Hamiltonians. A full-Hilbert-space ED check at L=7 matched its grand-canonical ground-state energy within 7e-15. Calibration runs on the XX chain (c=1, exact free-fermion energy) and SU(3) swap chain (c=2) are running alongside an L=48 candidate run. Both von Neumann and second-Renyi profiles will be retained, with their correlated systematic errors stated. No claim of c>2 is justified at this stage.

## Analytic threshold and larger-system work

DERIVATION.md now gives the explicit bond and the full two-body calculation. At the threshold family, the Schur complement at E=-4 factors as

    g^2 (s-s0)^2 [(2+s0)s+s0] / [(1+s)(1+s0)^2].

It is nonnegative for s in [0,1]. This proves the two-body global minima and yields the exact mass 15/7 for g=2, s0=1/2. It proves nothing by itself about a three-fluid many-body phase.

Added charge-conserving DMRG and translation-orbit ED. The fixed-charge DMRG matched L=8 ED to 3e-11. The momentum implementation initially exposed a normalization bug: NumPy sqrt on int16 periods produced float32. Casting periods to float64 fixed it; reconstructed momentum eigenstates then matched full-Hamiltonian eigenvectors with residuals below 2e-13 in the test.

Periodic spectra now reach L=21, N=7 at density 1/3. For delta=0.9, the k=2pi/L gaps at L=15,18,21 are 0.0192319, 0.0161025, 0.0104060. Their L times gap values are not yet stable. No z=1 conclusion follows. Fixed-density OBC DMRG runs at L=24,N=8 and L=48,N=16 are still converging and currently have much smaller entropy than a naive c=3 interpretation would suggest. XX and SU(3) calibration demonstrate substantial oscillatory corrections; SU(3) calibration is being refined to larger bond dimension before final fits.

## Boundary and convergence caveats

The periodic L=21 entropy fits at delta=0.9 are c_vN=3.055 and c_Renyi2=3.032. Nevertheless, the small-momentum gaps can also be read as approximately quadratic over the limited size window: L^2 gap is 4.33, 5.22, 4.59 at L=15,18,21. The asymptotic dynamical exponent remains unresolved.

Open chains have substantial boundary accumulation of charge, so comparing their total density N/L to a periodic chain can change the bulk filling at these sizes. A periodic-Hamiltonian DMRG implementation is therefore being tested. Its first L=8 validation missed the 1e-8 energy tolerance after six sweeps, with error 3.77e-7. The error decreases monotonically; the same test is being continued to twelve sweeps before this implementation is used for production.

Same- and opposite-valley two-pair energies were compared at L=18,24,30 for eleven interaction choices. Splittings are small and show size-dependent sign changes. These data do not establish a miscible valley mixture. Reducing V12 tends to favor the same-valley state in this screen.

## L=24 result and an independent finite-size warning

Translation-resolved L=24,N=8 ED gives E(k=0)=-14.18107312747 and E(k=2pi/L)=-14.17196929015. The resulting gap is 0.00910384, and L*gap=0.21849 agrees with L=21's 0.21853. Only selected additional sectors k=4,8,12 were tested at L=24; they are higher. All sectors were checked through L=21.

The k=0 L=24 entropy fits are c_vN=2.92389 and c_Renyi2=2.75315. The half-chain entropies are S=2.46463 and S2=1.88008. The two shape estimators still exceed 2 but decrease with size. More importantly, half-chain *size* scaling between L=18 and L=24 gives a much smaller estimate (about 1.23 from von Neumann entropy). Shell changes and few molecular particles can invalidate either fit. The inconsistency prevents declaring a counterexample.

Periodic DMRG at chi=96 significantly underestimates the L=24 exact entropy during convergence. Its energy error is also comparable to the small gaps. Larger-chain results at finite bond dimension must not be advertised as converged CFT data. A separate random fixed-charge initialization is running at L=48 as a metastability check.

Periodic calibration has been added using exact free-fermion correlation matrices and translation-resolved SU(3) swap ED. This uses the same entropy fits as the candidate. Open-chain SU(3) calibration refined to chi=192 gives c=2.031 from half-chain von Neumann size scaling between L=24 and 48, while second-Renyi size scaling is 2.326, illustrating its larger corrections.

## Scope of the candidate and numerical status

All fixed-density tests must still be embedded in a full Hamiltonian with a chemical potential. Positive finite-size compressibility and absence of a lower convex-envelope sector are being checked; fixing N alone is not sufficient. The explicit bond includes -mu(n_i+n_j)/2 for this reason.

A more systematic, momentum-resolved detuning scan at density 1/3 shows large shell effects. At L=18, delta=0.6 has raw shape estimates 4.52 and 5.47, while delta=0.9 gives 3.30 and 3.40. Such large finite-size numbers are not evidence for those central charges. Ground-state momentum changes and commensurate valley superpositions can inflate them. Delta=0.75 is being checked at L=21 and 24 as a nearby comparison, not selected solely because it maximizes a short-chain fit.

The two-pair comparison was extended to L=60. Opposite-valley energies are lower at L=36,48,60 by 3.54e-6, 1.02e-6, and 3.95e-7. The splitting decays very rapidly and does not by itself establish the bulk valley phase.

Periodic DMRG jobs were explicitly interrupted and resumed from their saved states after a contraction bottleneck was identified. NumPy's default einsum memory restriction chose a naive high-cost contraction for the boundary operator environments. Explicit matrix contractions give the same short-chain energy (validation error 6e-12) and permit larger bond dimensions. This is a performance correction, not a change to the Hamiltonian.

## Important negative result: the detuned molecular minima disappear

The two-body dispersion was rechecked at the detunings used for the finite-density entropy screen. At delta=0.15,0.30,0.45, local molecular minima remain, with positive binding below the continuum at the same momentum. At delta=0.60,0.75,0.90 they no longer exist. Solving the simultaneous stationary-point and zero-curvature conditions gives the saddle-node

    delta_crit = 0.56981795454,
    Q_crit = 1.00068764889,
    E_crit = -3.69864239469.

The residual of these two equations is 7.8e-13. Thus the large-entropy delta=0.75 and 0.9 points cannot be justified by the original three-minimum dilute-gas mechanism. A separate many-body mechanism would have to be demonstrated. This is a substantive obstruction, not a proof of their asymptotic central charge.

The fixed-charge stability check also finds parity and shell effects. At L=18, delta=0.9, the N=6 sector is locally on the convex envelope over N=4..8 for approximately -1.58245 < mu < -1.56046. N=7 is skipped in favor of a two-particle jump to N=8. The rational value mu=-1.57 lies inside this finite-size interval. This does not establish the thermodynamic density or absence of remote density jumps.
