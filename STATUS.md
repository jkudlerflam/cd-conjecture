# Where the c–d conjecture stands

Last updated 2026-11-08. This page lists only the main results. [REPORT.md](REPORT.md) is the chronological record, and [STATUS_DETAIL.md](STATUS_DETAIL.md) is the earlier long version of this page (construction step by step, per-run history, glossary). The standard throughout is theoretical physics, and nothing has had independent human review.

## The conjecture is false

Latorre and Sierra ([arXiv:2403.17242](https://arxiv.org/abs/2403.17242)) conjecture c ≤ d − 1 for translation-invariant nearest-neighbour chains on (ℂ^d)^⊗L whose low-energy spectrum is a unitary CFT.

1. Unbounded c at fixed local dimension d₀ = 4464, with c = n + 1 for every n (C39, C42, C44, C45). All finite inputs are proved and refereed. The one remaining input is the Wilsonian step from the microscopic chain to its Luttinger liquid, of the same kind as large-U Hubbard to Heisenberg. Its pieces are proved, and what is left is constructive RG at the level of mathematical physics.
2. The smallest counterexamples are explicit bonds at d = 11 (frozen cells, c = 11) and d = 13 (mobile cells, c = 14), written out in [contrib/small_d_note/](contrib/small_d_note/note.pdf). Both have certified finite members and are refereed and red-teamed. They rest on the same infrared inputs as the d = 24 candidate of a collaborator (GitHub: salpace).
3. Numerics agree. The logical vacancy model gives c = 2.0, 2.9 ± 0.1 and 3.98 ± 0.05 at n = 1, 2, 3, and the projected d₀ = 4464 chain gives c = 3.00 ± 0.02 at n = 2.

Not established: a rigorous lattice-to-CFT theorem for any counterexample.

## A replacement bound

The proposed replacement is A = Σ_a c_a v_a ≤ π M log₂ d, with M the gauge-optimized bond half-spread. It is proved for quasi-free chains (A ≤ 2M, sharp), for free-fermion and classical chains, and for WZW sectors coupled to on-site symmetries. Every solvable family satisfies it, and an adversarial search found no violation. The general case reduces to C(d) = sup L(E_{2L} − 2E_L)/M < ∞.

At d = 2 on the U(1) class (XXZ plus DM along z plus field) this is proved with C = 10781 everywhere except a thin strip −1 < J_z/J < −1 + δ near the ferromagnet. There it rests on one inequality on the Bethe-root density, which holds numerically with wide margin (REPORT §9.94 to §9.109). The parts with J_z/J ≤ 1 are refereed; the easy-axis side is not.

## Is c ≤ 1 at d = 2?

No nearest-neighbour qubit chain with c > 1 has been found, in about 103,000 sampled chains, targeted searches, or the literature (REPORT §9.110 to §9.124, refereed only where noted).

- c ≤ 1 is derived at the physics standard on every symmetric class: U(1), free fermion, XYZ, cyclic, parity, staggered charge, and the singular strata. On the time-reversal strata it holds given one assumption.
- The XZ point, the only bond with a second decoupled gapless species, gives no c > 1. Its induced interaction has the wrong sign for the c = 3/2 phase.
- Close to the ferromagnetic Heisenberg point every conformal bond has c ≤ 1, given standard semiclassics. The chain is a semiclassical ferromagnet there with at most one soft mode (refereed, and checked on 48 lattice bonds). The neighbourhood of its twisted images and the bound-state route to a second fluid are also covered.
- Where a symmetry makes two soliton valleys degenerate, the two-valley gas does not mix, so c = 1.

Open, in order of importance:
1. Generic bonds with no symmetry, and attractive liquids whose U(1)-breaking is as large as the exchange, far from every controlled region (adversarial Typhon search running).
2. Two specific mechanisms for a second gapless sector at strong U(1)-breaking: a staggered Ising transition inside an incommensurate liquid (exists classically, but at spin 1/2 the two orders are one SU(2)₁ field and it fails to all orders in weak coupling), and accidental two-valley soliton onsets (local search null). Large searches for both are pending.
3. Two points near XZ whose c readings still rise at L = 128 (up to 1.27), probably the slow XZ crossover; L = 256 runs are pending.

## Small d other than 2

At d = 3 to 7, nothing above d − 1 in symmetric families. A Casimir-type bound Σ c v ≤ 2π‖h‖ is false (SU(N) chains with rank-2 symmetric representations, from N = 7).

## Where help is most useful

1. Constructive RG for the Wilsonian step behind the counterexamples.
2. Independent review of the d = 11 and d = 13 bonds.
3. A general proof of the replacement bound, or its d = 2 strip.
4. The generic d = 2 residual.

All of this was produced with AI assistance (Claude Opus 5.5 coordinating, with gpt-6-astra runs).
