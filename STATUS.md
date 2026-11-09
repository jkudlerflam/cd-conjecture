# Where the c–d conjecture stands

Last updated 2026-11-09. This page lists only the main results. [REPORT.md](REPORT.md) is the chronological record, and [STATUS_DETAIL.md](STATUS_DETAIL.md) is the earlier long version of this page (construction step by step, per-run history, glossary). The standard throughout is theoretical physics, and nothing has had independent human review.

## The conjecture is false

Latorre and Sierra ([arXiv:2403.17242](https://arxiv.org/abs/2403.17242)) conjecture c ≤ d − 1 for translation-invariant nearest-neighbour chains on (ℂ^d)^⊗L whose low-energy spectrum is a unitary CFT.

1. Unbounded c at fixed local dimension d₀ = 4464, with c = n + 1 for every n (C39, C42, C44, C45). All finite inputs are proved and refereed. The one remaining input is the Wilsonian step from the microscopic chain to its Luttinger liquid, of the same kind as large-U Hubbard to Heisenberg. Its pieces are proved, and what is left is constructive RG at the level of mathematical physics.
2. The smallest counterexamples are explicit bonds at d = 11 (frozen cells, c = 11) and d = 13 (mobile cells, c = 14), written out in [contrib/small_d_note/](contrib/small_d_note/note.pdf). Both have certified finite members and are refereed and red-teamed. They rest on the same infrared inputs as the d = 24 candidate of a collaborator (GitHub: salpace).
3. Numerics agree. The logical vacancy model gives c = 2.0, 2.9 ± 0.1 and 3.98 ± 0.05 at n = 1, 2, 3, and the projected d₀ = 4464 chain gives c = 3.00 ± 0.02 at n = 2.

Not established: a rigorous lattice-to-CFT theorem for any counterexample.

Independent continuation (2026-11-09): [the new audit and calculations](contrib/counterexample_push_20261109/README.md)
retain the d=13 model-specific dressed-clock and relative defect
hypotheses, improve its legal projected pin activation bound to J/135,
and prove an exact two-body binding boundary for the qubit cusp. The
qubit reduction also has a nonzero distance-three pairing term at
relative order Delta/V; its off-resonant elimination requires a separate
argument. Neither result identifies a microscopic central charge.
An explicit [transverse-DM deformation of the qubit cusp](contrib/counterexample_push_20261109/FACILITATED_ROUTE.md)
now gives a constructive test family: its full number-mixing effective
model weakens the pair binding after valley retuning. The finite-ring
crossing drifts strongly with size; unbinding and c=2 remain unproved.

## A replacement bound

The proposed replacement is A = Σ_a c_a v_a ≤ π M log₂ d, with M the gauge-optimized bond half-spread. It is proved for quasi-free chains (A ≤ 2M, sharp), for free-fermion and classical chains, and for WZW sectors coupled to on-site symmetries. Every solvable family satisfies it, and an adversarial search found no violation. The general case reduces to C(d) = sup L(E_{2L} − 2E_L)/M < ∞.

At d = 2 on the U(1) class (XXZ plus DM along z plus field) this is proved with C = 10781 everywhere except a thin strip −1 < J_z/J < −1 + δ near the ferromagnet. There it rests on one inequality on the Bethe-root density, which holds numerically with wide margin (REPORT §9.94 to §9.109). The parts with J_z/J ≤ 1 are refereed; the easy-axis side is not.

## Is c ≤ 1 at d = 2?

No nearest-neighbour qubit chain with c > 1 has been established, in about 103,000 sampled chains, targeted searches, or the literature (REPORT §9.110 to §9.139, refereed only where noted). One unconverged candidate is under test (item 1 below). The early screens were weakly sensitive to c = 2. A recalibrated screen of the uniform sample bounds any open c > 1 region at about 0.1 percent of bond space, pending a second-stage check of two points.

- c ≤ 1 is derived at the physics standard at generic points of every symmetric class: U(1), free fermion, XYZ, cyclic, parity, staggered charge, and the singular strata, except one corner of the parity class near the Fibonacci point (item 1 below). On the time-reversal strata T and K_I it holds given a one-species assumption of the same kind as item 2 below, and the sheet K_II is open. Isolated multicritical points on the strata have not been inventoried.
- The XZ point, the only bond with a second decoupled gapless species, gives no c > 1. Its induced interaction has the wrong sign for the c = 3/2 phase.
- Close to the ferromagnetic Heisenberg point every conformal bond has c ≤ 1, given standard semiclassics. The chain is a semiclassical ferromagnet there with at most one soft mode (refereed, and checked on 48 lattice bonds). The neighbourhood of its twisted images and the bound-state route to a second fluid are also covered.
- Where a symmetry makes two soliton valleys degenerate, the two-valley gas does not mix, so c = 1.

Open, in order of importance:
1. A candidate c = 2 corner near the Fibonacci point (strong ZZ with the field near its blockade value), where the transverse hop is tuned down to the size of the next-nearest hop that the pairing generates at second order. The chain then behaves as a zigzag ladder whose band has two valleys, and DMRG on the derived effective model reads c ≈ 2 at L ≤ 96, against 1.88 for a c = 2 control. The runs are not converged and the full qubit chain has not been run. The exact two-body problem then showed that one particle from each valley binds into a pair, so the dilute gas is a pair liquid with c = 1. At the candidate density the pair gap may be too small to see at the sizes reached, which would make the c ≈ 2 reading a crossover, but a gapless window there is not excluded. Earlier one-species arguments for this region assumed the first-order hop dominates and do not cover it. Elsewhere the ladder picture gives one valley in every regime with a separation of scales (REPORT §9.137 to §9.139).
2. Generic bonds with no symmetry, and attractive liquids whose U(1)-breaking is as large as the exchange, far from every controlled region. The symmetric results do not propagate into this region, so it needs a dynamical argument of its own (REPORT §9.127). Besides two-species phases, this includes single-species c > 1 theories such as Z₅ parafermions, which could sit where commensurate locked phases melt, and three-state Potts order next to a floating liquid (c = 9/5). Both need at least two tunings and amplitudes of order the exchange, and a structured search for exact maps into the nearest-neighbour class found no c > 1 model (REPORT §9.134). The adversarial search (about 150,000 evaluations) has no survivor outside item 4.
3. Two specific mechanisms for a second gapless sector at strong U(1)-breaking: a staggered Ising transition inside an incommensurate liquid (exists classically, but at spin 1/2 the two orders are one SU(2)₁ field and it fails to all orders in weak coupling, and a 252-point quantum screen found no coexistence), and accidental two-valley soliton onsets (a search over about 48,000 classical bonds and 6,000 quantum restarts found none).
4. The XZ neighbourhood, where open-chain DMRG reads c ≈ 1.3 at L = 256. This appears to be a Néel-ordered gapped bulk seen through a slowly relaxing cat state, not a critical one. Edge fields along the order turn the open chain into the gapped VUMPS state at the same energy. Open chains to L = 512, which test this directly, are pending.

## Small d other than 2

At d = 3 to 7, nothing above d − 1 in symmetric families. A Casimir-type bound Σ c v ≤ 2π‖h‖ is false (SU(N) chains with rank-2 symmetric representations, from N = 7).

## Where help is most useful

1. Constructive RG for the Wilsonian step behind the counterexamples.
2. Independent review of the d = 11 and d = 13 bonds.
3. A general proof of the replacement bound, or its d = 2 strip.
4. The generic d = 2 residual.

All of this was produced with AI assistance (Claude Opus 5.5 coordinating, with gpt-6-astra runs).
