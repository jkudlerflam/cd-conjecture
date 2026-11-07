# Progress, rigorous toolbox map

2026-11-06.

Done.
- Read `specs/_Cd_target.md`, `specs/_bound_background.md` with updates, `specs/_reconnection_lemma.md`, `specs/C127_crossed_majorization.md`, REPORT.md sections 9.68 to 9.88, and RESULT.md of `codex/work_C125` to `codex/work_C128`.
- Surveyed families 1 to 7 plus complexity results, reading theorem statements where the PDF text was extractable (Borgs–Kotecký Thm 4.1, Sutter–Berta–Tomamichel Cor. 3.3, Hiai–Petz via Carlen–Lieb Eq. (1.10), Bluhm–Capel–Pérez-Hernández Thms I, II and Prop. 3.2, Kato–Brandão Thm 4 and Cor. 5, Kuwahara 2407.05835 main theorem, Wolf–Verstraete–Hastings–Cirac Eq. (3), Tasaki Thm 0 and Cor. 1a, DKKMT Thms 1, 2, 22, Benfatto–Mastropietro Eq. (1) and Thm 1.1, Kozlowski introduction and theorem) and abstracts otherwise. Tags in TOOLBOX.md say which.
- Scope update from Jonah (C(2) alone counts): added TOOLBOX Section 8 (d = 2 normal form, JW, integrable and stoquastic subclasses, rigorous spin-1/2 results, gaplessness), ranked a d = 2 direction first.
- Numerical sanity checks, numpy and scipy only, local: `d2_sanity.py` (output `d2_sanity_output.txt`, about 45 s) and `twist_check.py` (about 1 min).
- Wrote TOOLBOX.md and three draft specs: DRAFT_C129_d2_u1_bethe.md, DRAFT_C130_constructive_rg_casimir.md, DRAFT_C131_d2_cover.md. Nothing placed in specs/.

Not done or not verified.
- Several items are tagged recalled rather than read (Affleck–Lieb, Lieb–Schultz–Mattis 1961, Kotecký–Preiss, Nachtergaele–Sims, Anshu, Tóth, Björnberg–Ueltschi, the frustration-free gap-stability literature, the Hermitian content of the de Leeuw–Pribytok–Ryan classes, the exact Delta range in Yang–Yang).
- Whether Kozlowski's remainders are locally uniform in (Delta, D) was not checked; it is the first task of DRAFT_C129.
- No quantum-ring version of Giuliani–Mastropietro was found; absence in print is my reading of search results, not a certainty.
- No commits, no cluster use, no Codex.
