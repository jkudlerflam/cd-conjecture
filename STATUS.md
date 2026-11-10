# Where the c–d conjecture stands

Last updated 2026-11-10. This page lists only the main results. [REPORT.md](REPORT.md) is the chronological record, and [STATUS_DETAIL.md](STATUS_DETAIL.md) is an earlier long version of this page. The standard throughout is theoretical physics, and nothing has had independent human review.

Latorre and Sierra ([arXiv:2403.17242](https://arxiv.org/abs/2403.17242)) conjecture c ≤ d − 1 for translation-invariant nearest-neighbour chains on (ℂ^d)^⊗L whose low-energy spectrum is a unitary CFT.

## In short

1. The conjecture is false at large d. There are explicit counterexamples at d = 11 and d = 13, and c is unbounded at fixed d.
2. At d = 2 there is a candidate counterexample, an explicit nearest-neighbour qubit chain that reads c ≈ 2 on open chains up to L = 256, with a mechanism and an effective model that give the same answer. What separates it from an established counterexample is a short list of convergence checks, given below, most of them running.
3. A replacement bound, Σ_a c_a v_a ≤ π M log₂ d, holds in every solvable family and is proved in several classes.

## The conjecture is false at large d

1. Unbounded c at fixed local dimension d₀ = 4464, with c = n + 1 for every n (C39, C42, C44, C45). All finite inputs are proved and refereed. The one remaining input is the Wilsonian step from the microscopic chain to its Luttinger liquid, of the same kind as large-U Hubbard to Heisenberg.
2. The smallest counterexamples are explicit bonds at d = 11 (frozen cells, c = 11) and d = 13 (mobile cells, c = 14), written out in [contrib/small_d_note/](contrib/small_d_note/note.pdf). Both have certified finite members and are refereed and red-teamed. They rest on the same infrared inputs as the d = 24 candidate of a collaborator (GitHub: salpace), whose [independent audit](contrib/counterexample_push_20261109/README.md) leaves the model-specific hypotheses of the d = 13 bond in place as hypotheses.
3. Numerics agree. The logical vacancy model gives c = 2.0, 2.9 ± 0.1 and 3.98 ± 0.05 at n = 1, 2, 3, and the projected d₀ = 4464 chain gives c = 3.00 ± 0.02 at n = 2.

Not established: a rigorous lattice-to-CFT theorem for any counterexample.

## d = 2: a candidate counterexample with c = 2

The chain is

H = Σ_j [ (X_j X_{j+1} − Y_j Y_{j+1})/2 + (X_j Y_{j+1} − Y_j X_{j+1})/160 + 10 Z_j Z_{j+1} ] + h_z Σ_j Z_j, with h_z = −19.976047,

together with a family of neighbours of it. One of them is the bond proposed by the collaborator, which adds a small XX + YY term (hopping phase 1.53).

The strong ZZ term, with the field near its blockade value, leaves a dilute gas of blockaded particles. Their direct hop (the XY − YX term) is tuned down to the size of the next-nearest hop that the pairing term XX − YY generates at second order. The particles then move on a zigzag ladder whose band has two valleys, and the proposed infrared theory is a two-component Luttinger liquid with c = 2.

The spread of couplings is one small parameter ε = Δ/V (pair amplitude Δ = J_x − J_y, blockade energy V = 4J_z), not three tunings. The next-nearest hop exists only at second order, |t₂| = Δ²/2V, and the direct hop 2D is set equal to it (D = 1/160 above), giving the ladder εΔ : Δ : Δ/ε. The field sits at the blockade value plus a remainder of order t₂ that sets the density. As ε → 0 only the hopping phase, the pairing strength and the density matter, and the readings agree at V/Δ = 40 and 80. How large ε can be is open.

Evidence, all numerical:

- The qubit chain itself. Parity-resolved DMRG on open chains of L = 64, 128 and 256, with no projection, reads c = 1.95 for the bond above from the entropy difference between L = 128 and 256 (the ground branch is the odd sector, per the collaborator's [parity audit](contrib/d2_status_audit_20261110/README.md)), and 1.94 to 2.03 over six tuned bonds (two densities, three bonds, V/Δ = 40 and 80). The c = 1 and c = 2 controls calibrate both size pairs, the registered fail outcome (c ≤ 1.3) is excluded, and independent CPU and GPU codes agree (REPORT §9.167 to §9.181).
- The effective model. The second-order model of the blockaded particles reads c = 1.99 to 2.00 on the infinite chain at three densities, out to correlation lengths of 500 to 600 sites (REPORT §9.163, §9.179).
- The mechanism. One particle from each valley binds into a pair in the dilute limit, which alone would give a pair liquid with c = 1. At the candidate density the inter-valley backscattering responsible is irrelevant if the pseudospin Luttinger parameter K_s exceeds 1. K_s reads 1.2 to 1.5 in the effective model and 1.4 to 1.6 on the qubit chain at L = 256, where the pair exponent and K_c also match the effective model (REPORT §9.160, §9.161, §9.182).

Open, before this is called established:

1. Truncation at L = 256. The entropy still rises with bond dimension, by about 0.03 from χ = 512 to 768 against a registered tolerance of 0.01, so the registered pass is not declared. The bias lowers c. A χ = 1024 stage is running.
2. A third size. L = 512 is running (REPORT §9.184).
3. Entanglement scaling on the infinite qubit chain. Two attempts did not converge at the candidate, because the finite-χ state breaks particle number (REPORT §9.183, §9.187).
4. The pseudospin gap. In the effective model the pair gap is seen at low density and closes near n ≈ 0.05 to 0.10, below the candidate's 0.19, while n = 1/8 is marginal. The registered verdict is formally the gap outcome, through a fit that does not constrain n_c, and runs at larger χ are queued (REPORT §9.188).
5. Two unexplained readings. An older-solver run at V/Δ = 8 reads about 1.4 on the infinite chain, and the c = 1 control of the infinite-chain calculation sits narrowly outside its registered tolerance.

The collaborator's earlier runs of the phase-1.53 bond read a downturn from L = 64 to 128 on states that were not fully optimized ([record](contrib/d2_tenpy_20261109/RESULT.md), [rational bond](contrib/d2_dmrg_priority_20261109/SIMPLE_PARITY_CANDIDATE.md)). The converged rerun lies lower in energy and reads 2.2, so the downturn does not reproduce (REPORT §9.167). Their [transverse-DM deformation](contrib/counterexample_push_20261109/FACILITATED_ROUTE.md) reads c ≈ 1, and their [symmetry audit](contrib/d2_priority_audit_20261109/UNIFORM_Y_SOURCE.md) allows a small uniform Y source.

## Elsewhere at d = 2

Away from this corner no chain with c > 1 was found in about 103,000 sampled chains, targeted searches, or the literature (REPORT §9.110 to §9.150, refereed only where noted). A recalibrated screen of the uniform sample bounds any open c > 1 region at about 0.1 percent of bond space.

- c ≤ 1 is derived at the physics standard at generic points of every symmetric class: U(1), free fermion, XYZ, cyclic, parity (outside the corner above), staggered charge, and the singular strata. On the time-reversal strata T and K_I it holds given a one-species assumption, and the sheet K_II is open. Isolated multicritical points on the strata have not been inventoried.
- The XZ point, the only bond with a second decoupled gapless species, gives no c > 1. Its induced interaction has the wrong sign for the c = 3/2 phase.
- Close to the ferromagnetic Heisenberg point every conformal bond has c ≤ 1, given standard semiclassics (refereed, and checked on 48 lattice bonds).
- Where a symmetry makes two soliton valleys degenerate, the two-valley gas does not mix, so c = 1.

Open here:

1. Generic bonds with no symmetry, and attractive liquids whose U(1)-breaking is as large as the exchange. The symmetric results do not propagate into this region (REPORT §9.127). It includes single-species c > 1 theories such as Z₅ parafermions and three-state Potts order next to a floating liquid (c = 9/5). A structured search for exact maps into the nearest-neighbour class found no c > 1 model (REPORT §9.134).
2. Two mechanisms for a second gapless sector at strong U(1)-breaking, a staggered Ising transition inside an incommensurate liquid and accidental two-valley soliton onsets. Searches found neither.
3. The XZ neighbourhood, where open-chain DMRG reads c ≈ 1.3 at L = 256. This appears to be a Néel-ordered gapped bulk seen through a slowly relaxing cat state. Open chains to L = 512 are pending.

## A replacement bound

The proposed replacement is A = Σ_a c_a v_a ≤ π M log₂ d, with M the gauge-optimized bond half-spread. It is proved for quasi-free chains (A ≤ 2M, sharp), for free-fermion and classical chains, and for WZW sectors coupled to on-site symmetries. Every solvable family satisfies it, and an adversarial search found no violation. The general case reduces to C(d) = sup L(E_{2L} − 2E_L)/M < ∞.

At d = 2 on the U(1) class (XXZ plus DM along z plus field) this is proved with C = 10781 everywhere except a thin strip −1 < J_z/J < −1 + δ near the ferromagnet. There it rests on one inequality on the Bethe-root density, which holds numerically with wide margin (REPORT §9.94 to §9.109). The parts with J_z/J ≤ 1 are refereed; the easy-axis side is not.

## Small d other than 2

At d = 3 to 7, nothing above d − 1 in symmetric families. A Casimir-type bound Σ c v ≤ 2π‖h‖ is false (SU(N) chains with rank-2 symmetric representations, from N = 7).

## Where help is most useful

1. An independent reproduction of the d = 2 candidate. The Hamiltonian above, any DMRG code, L = 128 and 256 at bond dimension 768 or more.
2. Constructive RG for the Wilsonian step behind the large-d counterexamples.
3. Independent review of the d = 11 and d = 13 bonds.
4. A general proof of the replacement bound, or its d = 2 strip.

All of this was produced with AI assistance (Claude Opus 5.5 and Claude Fable 5 coordinating, with gpt-6-astra runs).
