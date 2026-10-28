# C31: Displaced Fermi seas and the infrared of the number-conserving switch cell (route A)

Read `../../specs/_routeA_background.md`, `../../contrib/routeA_assembly/README.md` §7, `../../contrib/routeA_compile/README.md`, `../work_C26/RESULT.md` (Sidon target, ρ² suppression) and `../../contrib/ir_control/README.md` first. If `PROGRESS.md` exists in your folder, continue from it.

**Your task.** The central-switch cell of `contrib/routeA_compile/` conserves fermion number exactly, so pairing is forbidden, but an exact symmetry between its ports makes the junction hopping Hermitian and centres every Fermi sea at 0 or π, which leaves inter-flavour pair tunnelling ψ†_{aR}ψ†_{aL}ψ_{bL}ψ_{bR} at zero momentum. At n = 2 a small arm asymmetry plus a carrier phase displaces the seas.
1. For general n, find a fixed-alphabet way to give every flavour a distinct displaced Fermi-sea centre (or distinct Fermi momenta with a Sidon-type design), within the switch cell's structure, with margins as functions of n. Verify on the effective chain at n = 2, 3, 4.
2. Alternatively, with all seas centred, compute at one loop the sign and size of pair tunnelling and the interflavour forward interactions the cell actually generates (from the λ³ and τ² non-Gaussian terms), and decide whether pair tunnelling is irrelevant (relative-mode Luttinger parameter below one) for an open set of parameters.
3. Quantify the scales: the reported λ_min/λ_max falls by about 10⁻³ per flavour; check that the hierarchy still closes at each fixed n and state the resulting energy scale of the logical chain.
4. Conclude whether the switch cell gives a stable c = n phase at fixed small coupling for every n, at the physics standard.

Rules: theoretical-physics standard, finite-dimensional statements proved. Tag every claim as proved, numerical (with sizes), heuristic, or assumed. Do not put probabilities, odds or percent-confidence estimates in any file. Python (numpy and scipy) only, inside your folder; no network; never modify files outside your folder; keep each computation to minutes. Write prose plainly, without em-dashes. Keep `PROGRESS.md` updated and write `RESULT.md` at the end.
