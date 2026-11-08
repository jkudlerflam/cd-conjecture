This run was done by Claude Opus agents (2026-11-07), not by Codex. It was interrupted and resumed: a first agent did the analytic sections and local ED, a second agent resubmitted the Typhon job, added the large-b scaling analysis and the b = 6, 8 checks, and wrote this file.

# C153 result: c <= 1 on the staggered-charge box

Derivations are in PROOF.md (sections cited as 1 to 8), the log in PROGRESS.md. Scripts: ed.py, kscan.py, kcheck.py, dsk.py, spec_ed.py, boundary_spec.py, lswt.py (ED, spin waves), fes.py and vumps.py (infinite-size FES), drun.py and udmrg.py (U(1) DMRG, from contrib/d2_staggered), fesfit.py (estimators following C158's fesfit.py, a copy of which is fesfit_c158.py), dmrg_c144.py (copy of C144's dmrg.py, unused). Typhon inputs in typhon/.

## Model

H = sum_j [X_j X_{j+1} - Y_j Y_{j+1} + u Z_j Z_{j+1} + b Z_j], b > 0, -1 < u < r = sqrt(1 + b^2/4), u != 0. A pi rotation about x on odd sites maps it to the XXZ chain with Delta = -u (J = 4 in spin-1/2 units) in a staggered field 2b, with the staggered charge Q = sum (-1)^j Z_j becoming 2 S^z_tot.

## Phase diagram of the box

| region | where | ground state | c | standard |
| --- | --- | --- | ---: | --- |
| P, polarized | -1 < u < u_BKT(b) | unique, gapped (pinned boson at small b, polarized product state at large b) | 0 | derived; proved on u = 0 (free, gap 2b) |
| BKT line | u = u_BKT(b) | K = 2, marginal cos(2 phi) | 1 | derived at small b, numerical elsewhere |
| G, liquid | u_BKT(b) < u < r | one compact boson for Q, K from 2 to infinity, v -> 0 at u = r | 1 | derived at small b and near u = r; heuristic in the large-b scaling regime; numerical between |
| factorization line | u = r | L + 1 ground states (one per Q), quadratic defect branch | not a CFT (z = 2) | proved |
| outside | u > r; u <= -1 | Neel; unique | 0 | proved |

The BKT line [derived at small b, numerical values]: u_BKT = 1/sqrt(2) + kappa b + O(b^2), from the dimension K of the staggered-field vertex cos(2 phi) at the XXZ point (K = pi/(2 arccos u) = 2 at u = 1/sqrt(2)), with kappa of about 0.25 from the data. At large b the strip is r - u_BKT = g_c b^{-3}, where g = b^3 (r - u) is the single scaling variable of the dilute two-species problem (PROOF 5.1: kinetic energy, conversion and detuning are all of order b^{-3} per particle at fixed g). Estimates of g_c: 0.025, 0.12, 0.40, 0.64 at b = 0.5, 1, 2, 4 (contrib/d2_staggered), below 0.69 at b = 3 (ring crossing), between 0.65 and 1.7 at b = 6 and between 0.77 and 2.0 at b = 8 (local FES). So g_c seems to saturate between about 0.6 and 1. [numerical]

Near the factorization line [proved at first order in r - u, derived for the liquid]: e''(0) = 8 r (r - u) exactly at first order (from the kernel states as projected product measures, checked by ED to 1 percent), the diamagnetic stiffness is 2/(1 + r), the quantum stiffness D_s(b) is below it (1.00, 0.98, 0.85, 0.59, 0.28 at b = 0, 0.5, 1, 2, 4, by ED), and K = pi sqrt(D_s/(8 r (r - u))), v = sqrt(8 r (r - u) D_s). Since K diverges like (r - u)^{-1/2}, the only symmetry-allowed gapping vertex cos(2 phi) becomes ever more irrelevant, so G is a nonempty open strip below the line at every b > 0 (given D_s(b) > 0, numerical).

## Is c <= 1 on the box?

Yes, at every gapless point found or derived. c = 1 throughout G and on the BKT line, carried by one Luttinger liquid of the exact charge Q. The arguments that no second (neutral) sector becomes critical:

1. [derived] Small b: the infrared theory is the XXZ compact boson perturbed by irrelevant operators; there is nothing else to become critical.
2. [derived, with proved ingredients] Near u = r at fixed b: the unperturbed zero modes are one state per Q-sector, a single U(1) family; the excitations at u = r form one quadratic branch (checked by ring ED against the exact one-defect energy, 68.99 against 69.15 at b = 1, L = 16). The perturbation (u - r) D_L turns that one z = 2 branch into one linear boson. Spin waves give one Goldstone branch; the "optical" branch at 4.0 is the folded zone-boundary value of the same band, not an independent mode.
3. [derived] Symmetry: the discrete symmetries (translation times spin flip, reflection, complex conjugation) have no order parameter outside the boson that could order or give an Ising line inside G; all boson vertices allowed are of dimension K or 4K > 2.
4. [heuristic] Large b: the only candidate second mode is the total density of the two sublattice species, and the third-order conversion e + o <-> vacuum locks its phase. Near u = r both species are condensed (each sublattice chain sits next to its ferromagnetic point), so the conversion vertex has dimension (1/4)(1/K_e + 1/K_o), far below 2. Its gap is of the same order b^{-3} as the liquid scale, so there is no parametric c = 2 window.
5. [numerical] contrib/d2_staggered: FES c between 0.95 and 1.05 at all 16 critical points with b = 0.5 to 4 (Ising x XX control 1.50 with the same estimator), open chains to L = 256 falling toward 1, no neutral level outside the boson tower. This run: 120 ring-ED points (b = 0.25 to 3) show only the gapped and liquid behaviours; local FES at b = 6 and 8 (chi <= 32) gives c of 0.9 to 1.1 where critical (least squares 0.96 at (b, g) = (6, 0.65), 1.05 / 0.88 at (8, 0.77)) and saturation where gapped.

Standard of the overall conclusion: theoretical-physics level (derived) in the small-b and near-line corners, heuristic in the large-b scaling regime g = O(1) with b > 4, numerical in between. The large-b statement does not depend on any one b: by the scale counting of PROOF 5.1 the corner is one universal problem in g, so the b = 4 to 12 data test it as a whole (collapse heuristic, rough check: v b^2 = 3.67, 3.80, 4.26 at (b, g) = (3, 0.54), (6, 0.43), (8, 0.51)).

## Candidates with c > 1

None. The only feature that could be read that way is at b = 6, 8 on rings L <= 20 (PROOF 7.3): a zero-momentum neutral level lies next to the boson oscillator doublet (ratio 0.74 to 1.08, rising with L). These rings contain fewer than three particles per species and are inside the z = 2 crossover (ring c(16, 20) = 0.5 to 0.6, below 1). [heuristic] The level is read as the two-defect continuum of that crossover, which at u = r approaches twice the one-defect level with L, not as a second critical mode. The open-chain strips at L = 64, 128 with two neutral excited states (Typhon, below) decide it. Earlier open-chain values above 1 near u = r (up to 3 at L = 64) are the known crossover from the factorization line and fall monotonically with L.

## Input to C154

C154's conditional theorem needs c <= 1 on the staggered-charge class (A4). This run supplies it at the derived level near b = 0 and near u = r, at the heuristic level for large b, and numerically for b <= 4 (plus low-chi checks at b = 6, 8). The B_2 statement inherits the same standard on this class.

## Pending Typhon job

- 17596442 (c153b): one 48-core node, partition typhon, 20 h limit, 45 tasks (typhon/tasks2.txt, typhon/c153b.slurm). PENDING (Priority) at the end of this run, estimated start 2026-11-12. It replaces 17594700, which was cancelled before running (and whose scratch copy lacked vumps.py, now fixed).
- Outputs in /scratch/lustre/jkudlerflam/c153/: fes_data/*.json (FES to chi = 128: controls xx and isxx; b = 3, 6, 8, 12 spanning g = 0.1 to 1.7; read with python3 fesfit.py fes_data/*.json, which prints g and c), strip_data/u*_b*_L{64,128}.json (open chains at b = 3, 6, 8 near u = r; q = 0 with one or two excited states and q = 1), grid_data/*.json (L = 64 grid at b = 0.25 and 6), logs/, slurm_17596442.out.
- What to read off: c from FES at b = 6, 8, 12 at fixed g (collapse, g_c, and c within 0.05 of 1 if the picture holds); in the strips, whether the second neutral level tends to twice the first (one boson) or grows like L (gapped total-density mode) rather than to a constant below 2 (a second critical sector).
