# C153: c ≤ 1 on the staggered-charge class

Your folder is `../work_C153/`. Read `specs/_qubit_c_context.md` first, then `contrib/d2_gap_reduction/staggered_slice.md` and `additive_charge.md` in full, `codex/work_C15/RESULT.md` §§2.3–2.4, and the C143 staggered-U(1) class (C143 PROOF §2.7).

## Setting

The only bonds with an exact additive charge beyond XXZ and commuting Ising are a(XX − YY) + uZZ + bZ. A π rotation about x on odd sites maps this to the XXZ chain (Δ = −u) in a staggered field b. Exact results exclude everything except the box −1 < u < √(1 + b²/4), b > 0, u ≠ 0, with a = 1. Zhang–Song (arXiv:2609.29728) report a gapless phase there from DMRG. C143 found that the Kurmann–Thomas–Müller frustration-free points lie in this class.

## Tasks

1. Derive the phase diagram inside the box: where the staggered field is relevant or irrelevant at the XXZ fixed point (its scaling dimension depends on K), the gapless region, its boundaries, and the universality class everywhere. Use bosonization, the exact factorization line u = r, the quadratic line u = 0, the exact Néel and gap theorems of d2_gap_reduction, and spin-wave and strong-coupling expansions.
2. Derive c ≤ 1 at every gapless point of the box at physics level, showing in particular that the conserved staggered charge carries a single Luttinger liquid and that no second sector becomes critical (the charge is exact, so the c = 1 boson is protected; the question is whether an extra neutral mode can also close its gap).
3. Calibrate with DMRG and ED (copy `codex/work_C144/dmrg.py`), on a grid of the box and at its boundaries, and record c with error bars. Use Typhon for the grid.
4. State the outcome as a theorem-style summary, tagged.
