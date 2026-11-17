# Sequential decisions

1. Saved L=96,chi=256 checkpoints are present (about 79 MB each).
   The prior raw outputs identify exactly V=4,kappa=a=1, full endpoint
   fields, both parity sectors. They will be read without modification.
2. The actual V=4 fixed-work L=144,chi=256 benchmark uses 0.793 GB RSS
   and about 0.55 seconds for a central 100-matvec update. Prior L=96
   chi=256 full-sweep timings, scaled by 144/96, remain below the
   registered 10-minute cap even at the observed maximum. This supports
   L=144,chi=256 production; it does not pre-approve chi=320.
3. All new density/bond/pair contractions, including adjacent overlapping
   bond products, agree with literal L=8 ED in both parities to below
   2e-9. The validation uses the physical bond in the computational gauge;
   density and entropy are invariant under both onsite changes of basis.
4. The saved L=96 states have density/profile peaks near q=0.884 pi;
   density RMS is 0.1085/0.1059 in even/odd sectors, and normalized
   connected peak weights are 0.01111/0.01084. A fine Fourier grid does
   not give commensurability precision finer than the finite window.
   The new L=48 states pass all local gates, and normalized connected
   weights are larger there (0.01678/0.01786), while the visible profile
   amplitudes remain comparable. This motivates the registered L=144 test.
5. As an explicitly exploratory summary of the already registered full
   correlation matrices, report connected RMS at separations M/3 through
   M/2 within the central-half window. It is not a new registered pass
   threshold. This separates long-distance behaviour from the local
   terms contributing to a structure-factor peak.
6. Add one explicitly secondary, unregistered bond channel P=XX+YY,
   including q=0 in its connected structure-factor search. Real-gauge
   finite states can restore its antiunitary symmetry, so a vanishing
   mean or balanced valley weights would not exclude valley order.
   This uses existing checkpoints only; no bulk parameter or production
   size is added. Validate its complex-gauge contractions against ED
   before measurements, and retain the original phase/accuracy gates.
7. The exact free controls use fixed hopping 1 and fixed chemical
   potential -2 cos(pi*0.296) across all sizes (finite occupation is
   floor(0.296*(leg_length+1))). This keeps their bulk parameters fixed
   while calibrating a density close to the V=4 state. Known pinned and
   translation-restored Neel examples distinguish disconnected from
   connected ordered weight. No control result is a candidate phase test.
8. The user subsequently requested a numerical central-charge estimate.
   Add a clearly post-hoc entropy analysis, without changing the fixed
   Hamiltonian or registered production cap. Report ground-branch and
   both fixed-parity size slopes, W=1,6,12 central averages, a common
   central fraction, and profile-fit sensitivity. Compare the preceding
   chi stage, and calibrate identical estimators on the exact one-/two-leg
   controls at the same sizes. Retain raw log(L) and averaged-log-chord
   versions to expose finite-window geometry effects. Periodic ED entropy
   is an additional small-size diagnostic, not a conformal-mode test.
   These finite-size slope estimates are denoted c_ent, not c_eff (which
   has a separate non-unitary-CFT meaning). No asymptotic c is assigned
   merely by rounding a fit or by fitting three sizes with three free
   parameters. The open-boundary coefficient c/6 and periodic coefficient
   c/3 follow Calabrese and Cardy, hep-th/0405152, equations (1)--(2).
9. Both completed L=144,chi=256 sectors fail the entropy chi-step gate
   (0.018942 even, 0.018397 odd) while passing the other local gates.
   The chi=320 benchmark uses 1.517 GB RSS. Its median fixed-work cost
   is 2.0565 times chi=256; multiplying the observed chi=256 sweep times
   gives 473 seconds at the median and 501 seconds at the maximum,
   within the registered 600-second projection cap. Execute the one
   permitted chi=320 continuation in each sector, starting the lower
   energy odd sector first while the independent L=96 run occupies the
   other worker. No chi above 320 is authorized by this registration.
   At chi=256 the optional boundary-probe amplitude ratios are only
   0.813/0.835, so its trigger is currently not met.
