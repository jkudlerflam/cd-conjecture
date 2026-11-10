# C188 follow-up: the corner point β = 0.2, t = 0.85 at L = 96 and 192

Coordinator (Claude Fable 5), 2026-11-09. Jonah asked for this check after REPORT §9.165. This file was written before any output of these runs existed.

## Why

C188's open-chain data leave one reading unsettled. At (α, β, γ) = (1, 0.2, −0.17), next to the boundary (t_c = 0.87 at β = 0.2), c_pair(48, 96) = 1.666, but the L = 96 run was not converged (last energy change 5e-6 at χ = 160) and no L = 192 run exists. The neighbouring points fall with size (t = 0.8: 1.322 → 0.939, t = 0.6: 1.144 → 1.057, from (48, 96) to (96, 192)), and the Heisenberg control reads 0.916 and 0.934.

## Runs

`run085.py` is C188's two-site DMRG (dmrg3.py, same model builder, same random start with seed 0), with the sweep loop copied so that it can log every sweep and continue past the old schedule. Schedule: χ = 16, 32, 64, 128, 160, 160 (the old schedule, so the first six sweeps repeat C188's run), then 240, 240, then 320 for 2 to 6 sweeps, stopping when the energy and S̄₄ are stationary. Sizes L = 48, 96, 192. S̄₄ is the mean entropy of the four central cuts and c_pair(L₁, L₂) = 6 [S̄₄(L₂) − S̄₄(L₁)] / ln(L₂/L₁), as in summarize_dmrg.py.

## Registered rule

Convergence, required at L = 96 and 192 before any reading is quoted as a result.
- K1. |ΔE| ≤ 1e-7 between the last two sweeps (total energy).
- K2. |ΔS̄₄| ≤ 1e-3 between the last two sweeps.
- K3. |S̄₄(χ = 320) − S̄₄(χ = 240)| ≤ 0.01.
- K4. Largest discarded weight in the last sweep ≤ 1e-8.
- K5. Largest |⟨X⟩|, |⟨Z⟩| ≤ 0.05. A larger value means the open chain sits on the ordered side and the reading is not a central charge.

If any of K1 to K5 fails, the outcome is "not converged" or "ordered", with the failing item named, and no verdict on c.

Verdict, from c_pair(96, 192) at χ = 320.
- W1. "Crossover to c = 1" if c_pair(96, 192) ≤ 1.15.
- W2. "Stays elevated" if c_pair(96, 192) ≥ 1.40. That makes the point a candidate that needs L = 384 and a ring check. It is not by itself a counterexample, because the velocity vanishes 0.02 away in t and a vanishing velocity inflates entropy slopes in c = 1 controls.
- W3. Anything between is unresolved at these sizes.

Also reported, without a rule attached: c_pair(48, 96) with the converged L = 96 state, the profile fits c_fit at each size, and whether the first six sweeps reproduce C188's L = 96 energy (−62.4641216654) and S̄₄ (1.3654).

## Files

sweeps_L*.jsonl (one line per sweep), stages_L*.jsonl (a full record at the end of each χ stage from 160 up, same keys as dmrg_tau.py), log_L*.txt.
