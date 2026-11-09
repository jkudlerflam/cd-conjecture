# C180 progress log

This run was done by a Claude Opus 5.5 agent (2026-11-09), not by Codex (Astra out of credits).

- Read specs/C175_C180_fibonacci_swarm.md in full, specs/_qubit_c_context.md, STATUS.md, REPORT §9.126 to §9.137, codex/work_C173/ (RESULT, PROOF, PROGRESS, twobody.py, valley_ed.py and its output), codex/work_C162/RESULT.md, codex/work_C164/RESULT.md, codex/work_C163/PROOF.md §4.2, codex/work_C156/RESULT.md (earlier literature survey, not repeated here).
- Folders work_C175 to work_C179 did not exist at the start of this run.
- Bosonization and operator table derived (PROOF 0 to 3). optable.py written and run (optable.out): band facts at phi = pi/2 checked, vertex-operator enumeration confirms that at the Sigma point only n_a = n_b in 2Z, p_a = -p_b in 2Z are momentum-zero on an open set, and none off Sigma. Dimension formula from the Gaussian covariance checked against the closed form. Drift-term stability check (lam = 1.2 unstable at the sample parameters).
- Two corrections made while checking the resonance table by hand: the 4kF-dressed pair transfer carries K_+ (not K_-), and the phi_- dressed valley transfer resonates at |w| = 4, not |w| = 0. The corrected table puts a pair-transfer resonance at n = 0.161 for |w| = 1, close to the candidate density.
- Literature (WebSearch/WebFetch, PDFs converted with pdftotext in the session scratchpad): KHGV arXiv:1112.5720, Hikihara et al. arXiv:1004.0622, Mishra-Greschner-Santos arXiv:1411.7270, Verresen-Vishwanath-Pollmann arXiv:1903.09179, Chepiga-Mila arXiv:2206.11754, Corcoran-de Leeuw-Pozsgay arXiv:2405.15848, Latorre-Sierra arXiv:2403.17242 read in full or in part; several others at abstract level (listed in RESULT.md).
- Main literature finding: MGS's blockaded frustrated ladder at t = |t'| is (0.2) without V2, V3 and is two-component (c = 2) at low density. The open question reduces to the effect of the generated attractions.
- No Typhon job, no long local runs. RESULT.md and PROOF.md written.
