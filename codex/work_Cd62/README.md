# Cd62: one GPT-6.2 Astra run on a finite C(d)

One run of `gpt-6.2-astra` at effort max, 2026-11-19, asked to prove or disprove that C(d) = sup L(E_{2L} − 2E_L)/M is finite. The brief is `TASK.md`. It was built from REPORT §9.68 to §9.109 and the specs `_Cd_target.md`, `_bound_background.md` and `_reconnection_lemma.md`, so the run did not know about C200 to C204 (REPORT §9.207). The outcome is neither a proof nor a disproof. Summary and context are in REPORT §9.210.

- `RESULT.md`, `PROOF.md`, `PROGRESS.md`: the run's own output, unedited.
- `check_xxz_doubling.py`: the run's check script. Run `python3 -B check_xxz_doubling.py` (standard library only).
- `claude_check/ed_strip.py`, `claude_check/ed_strip_weak.py`: an independent exact-diagonalization check of the strip bound by the coordinator (numpy and scipy), on random U(1) qubit bonds with field and DM angle, rings L = 2 to 7 against 2L, at strong, weak and zero field. The largest value found is 0.094 of the bound.
