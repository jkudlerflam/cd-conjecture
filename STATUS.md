# Where the c–d conjecture stands

Last updated 2026-11-24. This page lists only the main results and what is open. [REPORT.md](REPORT.md) is the chronological record, [UPDATES.md](UPDATES.md) keeps the dated announcements of each research round, and [STATUS_DETAIL.md](STATUS_DETAIL.md) is an earlier long version. The standard is theoretical physics (standard RG, controlled perturbation theory, calibrated numerics), and nothing has had independent human review. Conditional infrared predictions are not counted as established microscopic counterexamples.

Latorre and Sierra ([arXiv:2403.17242](https://arxiv.org/abs/2403.17242)) conjecture c ≤ d − 1 for translation-invariant nearest-neighbour chains on (ℂ^d)^⊗L whose low-energy spectrum is a unitary CFT.

## In short

1. At d = 2 there is a candidate counterexample, a family of explicit nearest-neighbour qubit chains that reads c ≈ 2 on open and infinite chains, with a mechanism and an effective model that give the same answer. Its long-distance stability is not settled.
2. At fixed d ≥ 9 there are explicit cell constructions with predicted c growing without bound. Their finite band data are certified, and the step to the microscopic CFT is conditional on stated infrared inputs.
3. A replacement bound, Σ_a c_a v_a ≤ π M log₂ d, holds in every solvable family and is proved in several classes, including the whole U(1) class at d = 2.

Current priority (set 2026-11-17): large c at d = 2, and unbounded c at the smallest fixed d, while tracking term counts, coefficient hierarchy and strength of evidence.

## d = 2: a candidate counterexample with c = 2

The chain is

H = Σ_j [ (X_j X_{j+1} − Y_j Y_{j+1})/2 + (X_j Y_{j+1} − Y_j X_{j+1})/160 + 10 Z_j Z_{j+1} ] + h_z Σ_j Z_j, with h_z = −19.976047,

together with a family of neighbours of it. One of them, proposed by Pace, adds a small XX + YY term (hopping phase 1.53). Off the symmetric line a momentum mismatch may cut off charge locking in the estimate of open item 5, making the phase-1.53 bond the stronger candidate in that analysis.

The strong ZZ term, with the field near its blockade value, leaves a dilute gas of blockaded particles. Their direct hop (the XY − YX term) is tuned down to the size of the next-nearest hop that the pairing term XX − YY generates at second order. The particles then move on a zigzag ladder whose band has two valleys, and the proposed infrared theory is a two-component Luttinger liquid with c = 2.

The spread of couplings is one small parameter ε = Δ/V (pair amplitude Δ = J_x − J_y, blockade energy V = 4J_z), not three tunings. The next-nearest hop exists only at second order, |t₂| = Δ²/2V, and the direct hop 2D is set equal to it (D = 1/160 above), giving the ladder εΔ : Δ : Δ/ε. The field sits at the blockade value plus a remainder of order t₂ that sets the density. As ε → 0 only the hopping phase, the pairing strength and the density matter, and the readings agree at V/Δ = 40 and 80.

The same family also contains a much simpler point at V/Δ = 4, on its tuning condition (REPORT §9.213). After a π/4 rotation it is the four-term chain

H(δ, h) = Σ_j [ ((1+δ)/2) X_j Y_{j+1} + ((1−δ)/2) Y_j X_{j+1} + Z_j Z_{j+1} ] − h Σ_j Z_j,

which reads c ≈ 2 in a narrow wedge below a staggered-chirality ordering line, with δ = 3/16, h = 1.783 the cleanest point. A c = 3/2 reading reported here on 2026-11-18 was withdrawn on 2026-11-19 (REPORT §9.209b).

Evidence, all numerical:

- The large-blockade bonds. Parity-resolved DMRG on open chains of L = 64, 128 and 256, with no projection, converged in bond dimension up to χ = 1024, reads c = 2.00 for the bond above and for the phase-1.53 bond from the entropy difference between L = 128 and 256, which passes the pre-registered test at both. From 256 to 512 at χ = 2048 they read 1.95 and 1.98, where an exact c = 2 control reads 1.95. Four further tuned bonds (two densities, V/Δ = 40 and 80) read 1.94 to 2.03 at χ = 768. Independent CPU and GPU codes agree (REPORT §9.178, §9.186, §9.192, §9.203).
- The four-term chain. At δ = 3/16, h = 1.783 the infinite chain reads 1.97 ± 0.01 and open chains from L = 96 to 384 read 1.96 ± 0.08. A nearby one-boson point reads 1.00 with the same estimator, and across the wedge c(h) reads 1.7 to 2.7 between L = 128 and 256. An independent DMRG code reproduces the open-chain states to L = 128 ([contrib](contrib/d2_xy_yx_chirality_imps_20261118/README.md), [review](contrib/qubit_c2_status_20261123/README.md), REPORT §9.212, §9.216).
- The effective model. The second-order model of the blockaded particles reads c = 1.99 to 2.00 on the infinite chain at three densities, out to correlation lengths of 500 to 600 sites. Its open chains read 1.82, 1.88 and 1.93 over the same three size doublings, approaching 2 from below while the qubit chain approaches from above (REPORT §9.163, §9.179, §9.205).
- The mechanism. One particle from each valley binds into a pair in the dilute limit, which alone would give a pair liquid with c = 1. At the candidate density the inter-valley backscattering responsible is irrelevant if the pseudospin Luttinger parameter K_s exceeds 1. K_s reads 1.2 to 1.5 in the effective model, 1.4 to 1.6 on the large-blockade chain at L = 256, and 1.8 to 2.0 on the four-term chain (REPORT §9.160, §9.161, §9.182, §9.214).

Open, before this is called established:

1. Convergence at the largest sizes. The entropy at L = 512 still moves by 0.010 to 0.012 from χ = 1536 to 2048, just above the registered 0.01 (REPORT §9.199, §9.203). On the four-term chain the L = 512 entries fail the convergence filter, and dropping the smallest length moves the fits appreciably.
2. Entanglement scaling on the infinite large-blockade chain. Two attempts did not converge, because the finite-χ state breaks particle number (REPORT §9.183, §9.187).
3. The pseudospin gap. In the effective model the pair gap is seen directly at low density and closes near n ≈ 0.07, against the candidate's 0.19. The registered test puts the candidate on the gapless side in substance, and falls short of the letter only because no single density shows the gap above three errors on both chains at once (REPORT §9.188, §9.204).
4. One unexplained reading. An older-solver run at V/Δ = 8 reads about 1.4 on the infinite chain.
5. Charge locking. The chain conserves only parity, so particle number is emergent. A number-changing term O₄ (two inter-valley pairs, scaling dimension about 1.25) is relevant when its momentum matches, and on the symmetric line it does. Its coupling is about (16/9) C₄ ε² with C₄ not computed, which for C₄ of order one gives a locking length of several thousand sites, beyond which c ≤ 1. At the phase-1.53 bond the coupling oscillates with a period of about 420 sites and the liquid survives in this estimate. This is an order-of-magnitude estimate, and chains of L ≤ 512 cannot test it (REPORT §9.161, §9.162; [collaborator note](contrib/collaborator_notes/qubit_c2/notes_qubit_c2.pdf), Sec. 6). A finite-ring quartet response has been computed but is not yet a continuum coupling ([quartet](contrib/qubit_detuned_quartet_20261118/README.md), [dressed modes](contrib/qubit_finite_density_modes_20261118/README.md)).
6. Operator identification. On the four-term chain the measured exponents fix the locking question only once one knows which operators carry them. Under one natural two-boson dictionary the quartet is the only relevant generic locking harmonic, while a different assignment reproduces the same exponents with no relevant locking. The pair momentum is now measured to be π, so the quartet sits at momentum zero, and if the pair is a vertex operator its dimension is about 1.84, weakly relevant. Measured directly, the quartet reads x4 ≈ 2.05 ± 0.2, drifting up with L, so it is not resolved against marginal and is now only weak evidence against an asymptotic c = 2. Koo–Saleur on rings does not resolve the chain at reachable bond dimension. The bulk velocities and the compactification are the missing measurements ([audit](contrib/qubit_ir_cft_team_20261123/README.md), [list](contrib/qubit_ir_cft_team_20261123/TODO.md), REPORT §9.215, §9.216, §9.217).

Pace's earlier runs of the phase-1.53 bond read a downturn from L = 64 to 128 on states that were not fully optimized ([record](contrib/d2_tenpy_20261109/RESULT.md), [rational bond](contrib/d2_dmrg_priority_20261109/SIMPLE_PARITY_CANDIDATE.md)). The converged rerun lies lower in energy and reads 2.2, so the downturn does not reproduce (REPORT §9.167). Their [transverse-DM deformation](contrib/counterexample_push_20261109/FACILITATED_ROUTE.md) reads c ≈ 1, and their [symmetry audit](contrib/d2_priority_audit_20261109/UNIFORM_Y_SOURCE.md) allows a small uniform Y source.

## Elsewhere at d = 2

Away from this family no chain with c > 1 was found in about 103,000 sampled chains, targeted searches, or the literature (REPORT §9.110 to §9.150, refereed only where noted). A recalibrated screen of the uniform sample bounds any open c > 1 region at about 0.1 percent of bond space.

- c ≤ 1 is derived at the physics standard at generic points of every symmetric class: U(1), free fermion, XYZ, cyclic, parity (outside the corner above), staggered charge, and the singular strata. On the time-reversal strata T and K_I it holds given a one-species assumption, and the sheet K_II is open. Isolated multicritical points on the strata have not been inventoried.
- The XZ point, the only bond with a second decoupled gapless species, gives no c > 1. Its neighbourhood, which reads c ≈ 1.3 on open chains to L = 256, is a gapped Néel bulk seen through a slowly relaxing cat, confirmed at L = 512.
- Close to the ferromagnetic Heisenberg point every conformal bond has c ≤ 1, given standard semiclassics (refereed, and checked on 48 lattice bonds).
- Where a symmetry makes two soliton valleys degenerate, the two-valley gas does not mix, so c = 1.

Open here: generic bonds with no symmetry, and attractive liquids whose U(1)-breaking is as large as the exchange (REPORT §9.127). This region includes single-species c > 1 theories such as Z₅ parafermions and three-state Potts order next to a floating liquid (c = 9/5). A structured search for exact maps into the nearest-neighbour class found no c > 1 model (REPORT §9.134). Searches for a second gapless sector at strong U(1)-breaking, through a staggered Ising transition inside an incommensurate liquid or accidental two-valley soliton onsets, found neither. The c ≈ 2 wedge they turned up instead belongs to the candidate family above.

## Conditional constructions with c > d − 1

These are explicit chains of gapped cells whose logical particles have 2n simple Fermi points, so the predicted infrared theory has c growing with the cell length n at fixed d. Their finite ingredients are certified. The step to the microscopic CFT rests on stated infrared inputs, namely matching to the logical model, selection of the intended ground sector, and Gaussian RG stability. Not established: a lattice-to-CFT theorem for any of them.

1. d = 9 and 10, the [markerless chains](contrib/markerless_d9_d10/README.md). Nine material states with frozen cells give predicted c = n, and one uncoloured vacancy gives c = n + 1 at d = 10. Every head moves by one rule, so the legal sector and the cell pattern are exact invariants, which removes finite-K leakage. The band data are certified with C84's pipeline at n = 9, 10, 11, 12, 14, 16 and 20 (predicted c up to 21 at d = 10), not for all n. An AI referee confirmed the frozen chain and re-checked the corrected mobile bond.
   - A modified d = 9 bond, with a biased clock and an L1 readout, has [2n simple, pair-sum-separated roots for every n](contrib/markerless_all_n_20261118/README.md), plus [local interaction bounds independent of volume and exact ground-sector selection](contrib/markerless_uniform_matching_20261118/README.md) on compatible rings at nonzero port coupling.
   - For the d = 10 mobile chain, [finite certificates](contrib/markerless_mobile_control_20261118/README.md) compare the intended sector with competitors at N = 20, V = 2, J/t = 1/25, and one uniform susceptibility channel at half occupancy is [proved](contrib/markerless_thermodynamic_selector_20261118/README.md). Thermodynamic mobile selection, the dressed-band J/t budget and the severe coefficient hierarchy remain open.
2. d = 11 (frozen, c = 11) and d = 13 (mobile, c = 14), in [contrib/small_d_note/](contrib/small_d_note/note.pdf). Certified finite members with AI referee checks, and a red-team review of d = 13. The [independent audit](contrib/counterexample_push_20261109/README.md) leaves the d = 13 model-specific hypotheses open. A [d = 12 mobile proposal](contrib/mobile_below_13/README.md) has additional heuristic dynamical steps.
3. d₀ = 4464, with predicted c = n + 1 for every n (C39, C42, C44, C45). The finite inputs are proved and refereed. The conclusion is conditional on the stated Wilsonian input.
4. Growing d with c/d → ∞. A sparse species-graph family has predicted c = (p+1)(p²+p+1) + 1 at d = O(p² log p) ([argument](contrib/research_continuation/PROPOSED_COUNTEREXAMPLE.md), [infrared audit](contrib/research_continuation/IR_AUDIT.md)).
5. Numerics agree. The logical vacancy model gives c = 2.0, 2.9 ± 0.1 and 3.98 ± 0.05 at n = 1, 2, 3, and the projected d₀ = 4464 chain gives c = 3.00 ± 0.02 at n = 2.

## A replacement bound

The proposed replacement is A = Σ_a c_a v_a ≤ π M log₂ d, with M the gauge-optimized bond half-spread. It is proved for quasi-free chains (A ≤ 2M, sharp), for free-fermion and classical chains, and for WZW sectors coupled to on-site symmetries. Every solvable family satisfies it, and an adversarial search found no violation. The general case reduces to C(d) = sup_L L(E_{2L} − 2E_L)/M < ∞.

- d = 2, U(1) class (XXZ plus DM along z plus field). C(2) is finite on every bond, every row proved and refereed by AI agents (REPORT §9.94 to §9.109, §9.207). The ferromagnetic strip is closed by a uniform counting bound on the Bethe roots, and the easy-axis side has C = 10781. On −1 < J_z/J ≤ −1/√2 the constant is explicit, L(E_{2L} − 2E_L)/M ≤ 2π sin η + π²/4 ≤ 6.91 with J_z/J = −cos η, which also lowers the middle range to about 10²⁰ (REPORT §9.210, checked by the coordinator, not yet refereed).
- General d, three theorems (REPORT §9.207). For every chain, (L e₀ − E_L)/M → 0 uniformly in the bond at fixed d, a qualitative result with no content at reachable sizes that does not bound A. Two-species exchange chains with any internal register obey A ≤ 2M at zero interaction. Critical ferromagnetic Potts chains at Q = 3, 4 obey a crossover-free bound inherited from XXZ, which covers the neutral c = 4/5 sector.
- Cell constructions. For gapped cells coupled at second order through one bond per interface, with quadratic effective dynamics, the sum of Fermi speeds is at most C d⁴ times the largest bond norm, uniformly in the cell length and the number of modes. The budget extends to all orders for g/Δ ≤ 1/(1728√10 d) as long as the effective dynamics stays quadratic. With interactions the sum of renormalized velocities is not controlled (REPORT §9.210, §9.211, not yet refereed). The Fermi velocities of the d = 9 chain lie within a factor 76 of each other at n = 9, not on the 10^−2.74 ladder quoted earlier (REPORT §9.211).

A general finite C(d) is still open.

## Small d other than 2

At d = 3 to 7, nothing above d − 1 in symmetric families. A Casimir-type bound Σ c v ≤ 2π‖h‖ is false (SU(N) chains with rank-2 symmetric representations, from N = 7).

## Where help is most useful

1. An independent reproduction of the d = 2 candidate. Either Hamiltonian above, any DMRG code, L = 128 and 256 at bond dimension 768 or more.
2. The missing qubit measurements in the [prioritized list](contrib/qubit_ir_cft_team_20261123/TODO.md), above all the pair momentum and the locking amplitude.
3. Constructive RG for the Wilsonian step behind the cell constructions.
4. Independent review of the markerless d = 9 and d = 10 bonds, especially the corrected mobile selector, its quantitative J/t budget and all-n band separation.
5. A general proof of the replacement bound (a finite C(d) at any d ≥ 2 beyond the classes above), or a reasonable constant on −1/√2 ≤ J_z/J ≤ 0 at d = 2.

All of this was produced with AI assistance (Claude Opus 5.5 and Claude Fable 5 coordinating, with gpt-6-astra runs).
