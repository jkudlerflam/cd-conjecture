# C17: Exclude, one by one, the mechanisms that could give c > 1 in a qubit chain

Read `../../specs/_d2_background.md` first; it defines the problem, the parametrization, what is known and the rules. If `PROGRESS.md` exists in your folder, continue from it.

**Your part: (c) of the d = 2 program, the analytic side.** List every physical mechanism by which a nearest-neighbour chain can produce c > 1, and for each decide whether it can operate at d = 2, with an argument at the physics standard.

Mechanisms to cover at least:
1. **Spontaneous enlargement of the unit cell.** Long spontaneous cells (frustration-free selection is closed by K_FF(2) = 4; what about kinetic selection, quantum-walk energies as in C7's selector, or period-2 and period-3 order?), and whether any period-p order can carry more than one gapless mode at d = 2. Note that blocking p sites gives local dimension 2^p but with p-site translation, so the question is genuinely about spontaneous order.
2. **Emergent non-abelian currents.** SU(2)_k with k ≥ 2, SO(3)₁, or products, for bonds without SU(2) symmetry. Use anomaly and level arguments (Lieb–Schultz–Mattis with whatever symmetry survives, the embedding-index level ceiling of `../../REPORT.md` §3.4) and decide what remains when the bond has no symmetry.
3. **Multi-fluid mixtures.** Unpaired particles plus bound clusters (pairs of adjacent up spins, longer strings), in the presence and absence of a U(1). Conversion vertices, their dimensions and momenta, and the range of Luttinger parameters nearest-neighbour qubit bonds can produce.
4. **Fractionalization and partons.** Spinon–holon or Majorana partons with emergent gauge fields, and whether the parton level bounds of `../../REPORT.md` §8 apply at d = 2.
5. **Decoupled sublattices.** Interleaved chains need range R ≥ 2 (c = R at range R). Explain precisely why one-site translation with range 1 cannot decouple sublattices, and whether a spontaneous version can.
6. **Anything else you find** (for example multicritical points glued from two c ≤ 1 sectors, or non-compact or non-unitary limits that might be mistaken for c > 1).

For each item give the argument, its status tag, and a concrete numerical test of the conclusion where one exists. End with the list of mechanisms that remain open at d = 2, ranked by how plausible they are.
