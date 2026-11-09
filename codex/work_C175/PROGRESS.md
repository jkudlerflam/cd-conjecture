# C175 progress log

This run was done by a Claude Opus 5.5 agent (2026-11-09), not by Codex (Astra out of credits).

- Read the swarm spec, C173 (RESULT, PROOF 1 to 7, twobody.py, valley_ed.py and .out, ring_ed.py, block_dmrg.py), C162 RESULT and PROOF 0 to 2 and 10, C164 RESULT, C163 PROOF 4.2, the shared context and REPORT 9.137.
- Pre-registered P1 to P3 in RESULT.md before computing anything.
- scat.py: exact lattice threshold solver (free double-root modes plus two boundary rows). Hard-rod check a = 1 passed. The unlike-pair matrix A came out non-Hermitian with anti-Hermitian part diag(∓iγ): γ = 0.0630, 0.1291, 0.2887 at |w| = 0.5, 1, 2, identical with V₂ = V₃ = 0, equal to ε‴/(6ε″). Identified as kinematic (cubic relative dispersion), Hermitian part used.
- First physics output: at |w| = 1, a_s = −3.173, A_d = 10.5, |A_o| = 10.211, so Δ_σ = −1.34. The large A_d suggested a shallow bound state.
- Bound-state scan (scat.bound_scan, relative box 500): lowest level minus 2ε_min over K is minimal at K = ±π, −0.0017, −0.0195, −0.1225 at |w| = 0.5, 1, 2; at K = 0 it is +0.125, +0.500, +2.0 and at K = ±π/4 +0.53, +0.14, +0.11 (above threshold). Bound state below the gas threshold at K = π.
- ring2.py (independent brute-force two-particle ring ED): −0.01958, −0.01953 at L = 120, 240 (|w| = 1), −0.12255 (|w| = 2), −0.00175 at L = 240 (|w| = 0.5). Confirms. ring2.out.
- scan_w.py: |w| from 0.25 to 3.9 (scan_w.out). Bound everywhere, E_b ≈ ε″/a_+² at small |w|. The like-pair a_s diverges near |w| ≈ 3.63.
- analytic_A.py (sympy): closed form A_d = 5/2 + 8/u², a_∓ = [5u² + 16 ∓ √(u⁴ + 160u² + 256)]/(2u²). analytic_res.py: det of the α block ∝ D = 2V₂V₃ + 4V₂ + V₃u² + 4V₃ + 8, so V₃ = 0, V₂ = −2 is an exact zero-energy resonance for all u, and V₃ = −2/3 detunes it to D = −2u²/3.
- nbody_box.py (gap coordinates, hard walls): N = 2 reproduces E_b. Trimers bound at |w| = 1 (6e-4) and 2 (0.018), tetramer and pentamer bound at |w| = 2, tetramer at |w| = 1 not resolved by R = 110. Initial K choice for |w| = 2 trimer was wrong (0.2527, the |w| = 1 value), redone at the right K. Outputs nbody3.out, nbody3b.out, nbody45.out, nbody4_w1_R110.out.
- ring_k.py (momentum-resolved ring ED, hard-core bosons, validated against ring2.py to 1e-14 at N = 2, L = 60, and the configuration count against brute force). ring_eos.py: E0(N, L) and spectra per momentum at |w| = 1 for L = 12, 18, 24, 30, 36 (ring_eos_w1.jsonl, ring_eos_w1b.jsonl). Ground states of even N at K = Nπ/2, polarized sectors 0.46 to 0.66 higher, flip energies, odd-even energies (PROOF 6).
- phi_scan.py: pair binding off φ = π/2 (phi_scan.out).
- smatrix.py: energy-dependent eigenphases (smatrix.out), illustrative only.
- Coordinator note (C180's criterion and calibration) received mid-run. calib.py (calib.out): V₂ = V₃ = 0 is two-component (Δ_σ = 0.125, 0.25, 0.50 at |w| = 0.5, 1, 2), reproducing the MGS limit; V₂ alone is the exact resonance; V₃ alone stays two-component; both together bind the unlike pair and put a_s below the window. Folded into PROOF 1, 2.4 and RESULT.
- ring_clusters.py (ring_clusters.out): ring check of the cluster energies, binding confirmed, finite-ring oscillations of a few 1e-3.
- Machine load average rose to about 75 during the run (other agents). block_dmrg.py and dmrg_d.py copied unmodified from C173 for the Typhon job. Local pilots in local/ (Lb = 12) to set the chemical potentials: at |w| = 2, μ = −2.9 gives n = 0.355 with the full V's against n = 0.066 with V₂ = V₃ = 0, so the attraction pulls the density up sharply (consistent with clustering).
- Typhon job 17604038 submitted (14 tasks, typhon/tasks.txt, prediction registered in RESULT.md before submission). Local pilot: |w| = 1, μ = −2.24 gives n = 0.091 at Lb = 12.
